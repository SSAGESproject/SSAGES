/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Hythem Sidky <hsidky@nd.edu>
 *                Jonathan K. Whitmer <jwhitme1@nd.edu>
 *                Ben Sikora <bsikora906@gmail.com>
 *
 * SSAGES is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SSAGES is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SSAGES.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "Meta.h"
#include <math.h>
#include <iostream>
#include "Drivers/DriverException.h"
#include "CVs/CVManager.h"
#include "Validator/ObjectRequirement.h"
#include "Drivers/DriverException.h"
#include "Snapshot.h"
#include "schema.h"

using namespace Json;

namespace SSAGES
{
	//! Evlauate Gaussian. Helper function.
	/*!
	 * \param dx x-value.
	 * \param sigma Width of Gaussian
	 * \return Value at x of Gaussian with center at zero and width sigma.
	 */
	double gaussian(double dx, double sigma)
	{
		double arg = (dx * dx) / (2. * sigma * sigma);
		return exp(-arg);
	}

	//! Evaluate Gaussian derivative. Helper function.
	/*!
	 * \param dx Value of x.
	 * \param sigma Width of Gaussian.
	 * \return Derivative at x of Gaussian with center at zero and width sigma.
	 */
	double gaussianDerv(double dx, double sigma)
	{
		double arg =  (dx * dx) / (2. * sigma * sigma);
		double pre = - dx / (sigma * sigma);
		return pre * exp(-arg);
	}

	// Pre-simulation hook.
	void Meta::PreSimulation(Snapshot* snapshot, const CVManager& cvmanager)
	{
		auto cvs = cvmanager.GetCVs(cvmask_);
		// Write ouput file header.
		if(IsMasterRank(world_))
		{
			hillsout_.open("hills.out");
			hillsout_ << "#Iteration "; 

			for(size_t i = 0; i < cvs.size(); ++i)
				hillsout_ << "center." << i << " ";
	
			for(size_t i = 0; i < cvs.size(); ++i)
				hillsout_ << "sigma." << i << " ";
			
			hillsout_ << "height" << std::endl;
				
			hillsout_.close();
		}

		// Initialize grid to zero. 
		if(grid_ != nullptr)
		{
			Vector vec = Vector::Zero(cvs.size());
			std::fill(grid_->begin(), grid_->end(), vec);
		} 

		auto n = snapshot->GetTargetIterations();
		n = n ? n : 1e5; // Pre-allocate at least something.
	
		hills_.reserve(n+1);
		widths_.reserve(n+1);
		derivatives_.resize(cvs.size());
		tder_.resize(cvs.size());
		dx_.resize(cvs.size());
	}

	// Post-integration hook.
	void Meta::PostIntegration(Snapshot* snapshot, const CVManager& cvmanager)
	{
		auto cvs = cvmanager.GetCVs(cvmask_);
		// Add hills when needed.
		if(snapshot->GetIteration() % hillfreq_ == 0)
			AddHill(cvs, snapshot->GetIteration());

		// Always calculate the current bias.
		CalcBiasForce(cvs);

		// Take each CV and add its biased forces to the atoms
		// using the chain rule.
		auto& forces = snapshot->GetForces();
		auto& virial = snapshot->GetVirial();

		for(size_t i = 0; i < cvs.size(); ++i)
		{
			auto& grad = cvs[i]->GetGradient();
			auto& boxgrad = cvs[i]->GetBoxGradient();
			
			// Update the forces in snapshot by adding in the force bias from each
			// CV to each atom based on the gradient of the CV.
			for (size_t j = 0; j < forces.size(); ++j)
				for(size_t k = 0; k < 3; ++k)
					forces[j][k] -= derivatives_[i]*grad[j][k];
			
			virial += derivatives_[i]*boxgrad;
		}
	}

	// Post-simulation hook.
	void Meta::PostSimulation(Snapshot*, const CVManager&)
	{
	}

	// Load hills from file. 
	void Meta::LoadHills(const std::string& filename)
	{
		std::ifstream file(filename);
		std::string line; 

		auto dim = widths_.size();
		double iteration = 0, height = 0;
		std::vector<double> width(dim, 0.), center(dim, 0);

		// First line is a comment. 
		std::getline(file, line);
		while(std::getline(file, line))
		{
			std::istringstream iss(line);
			iss >> iteration; // Not really using this. 
			
			// Load centers.
			for(size_t i = 0; i < dim; ++i)
				iss >> center[i];
			
			// Load widths.
			for(size_t i = 0; i < dim; ++i)
				iss >> width[i];
			
			// Load height. 
			iss >> height;

			hills_.emplace_back(center, width, height);
		}
	}

	// Drop a new hill.
	void Meta::AddHill(const CVList& cvs, int iteration)
	{
		size_t n = cvs.size();

		// Assume we have the same number of procs per walker.
		int nwalkers = world_.size()/comm_.size();

		// We need to exchange CV values across the walkers 
		// and to each proc on a walker.	
		std::vector<double> cvals(n*nwalkers, 0);

		if(IsMasterRank(comm_))
		{
			for(size_t i = 0, j = GetWalkerID(world_,comm_)*n; i < n; ++i,++j)
				cvals[j] = cvs[i]->GetValue();
		}

		// Reduce across all processors and add hills.
		MPI_Allreduce(MPI_IN_PLACE, cvals.data(), n*nwalkers, MPI_DOUBLE, MPI_SUM, world_);
		
		for(size_t i = 0; i < n*nwalkers; i += n)
		{
			std::vector<double> cval(cvals.begin() + i, cvals.begin() + i + n);
			hills_.emplace_back(cval, widths_, height_);
			
			// Write hill to file.
			if(IsMasterRank(world_))
				PrintHill(hills_.back(), iteration);
		}

		// If grid is defined, add bias onto grid. 
		if(grid_ != nullptr)
		{
			std::vector<double> dx(n, 0.0), df(n, 1.0);
			auto& hill = hills_.back();
			for(auto it = grid_->begin(); it != grid_->end(); ++it)
			{
				auto& val = *it;
				auto coord = it.coordinates();

				// Compute difference between grid point and current val. 
				for(size_t i = 0; i < n; ++i)
				{
					dx[i] = -cvs[i]->GetDifference(coord[i]);
					df[i] = 1.;
				}

				// Compute derivative.
				for(size_t i = 0; i < n; ++i)
				{
					for(size_t j = 0; j < n; ++j)
					{
						if(j != i) 
							df[i] *= gaussian(dx[j], hill.width[j]);
						else
							df[i] *= gaussianDerv(dx[j], hill.width[j]);
					}
				}

				// Add to grid. 
				for(size_t i = 0; i < n; ++i)
					val[i] += height_*df[i];
			}
		}
	}

	// Writes hill to output file. This should only be called by the 
	// world master node. 
	void Meta::PrintHill(const Hill& hill, int iteration)
	{
		hillsout_.open("hills.out", std::fstream::app);
		
		hillsout_ << iteration << " ";
		hillsout_.precision(8);
		
		for(auto& cv : hill.center)
			hillsout_ << cv << " ";
		
		for(auto& w : hill.width)
			hillsout_ << w << " ";

		hillsout_ << height_ << std::endl;
		hillsout_.close();
	}

	void Meta::CalcBiasForce(const CVList& cvs)
	{	
		// Reset bias and derivatives.
		double bias = 0.;
		auto n = cvs.size();

		// Reset vectors.
		std::fill(derivatives_.begin(), derivatives_.end(), 0);
		
		// Look up and apply grid bias. 
		if(grid_ != nullptr)
		{
			bool inbounds = true;
			std::vector<double> val(n, 0.);
			for(size_t i = 0; i < n; ++i)
			{
				val[i] = cvs[i]->GetValue();
				if(val[i] < grid_->GetLower(i) || val[i]  > grid_->GetUpper(i))
					inbounds = false;
			}

			if(inbounds)
			{
				auto frc = (*grid_)[val];
				for(size_t i = 0; i < n; ++i)
					derivatives_[i] = frc[i];
			}
			else
			{
				if(IsMasterRank(comm_))
				{
					std::cerr << "Metadynamics: out of bounds ( ";
					for(auto& v : val)
						std::cerr << v << " "; 
					std::cerr << ")" << std::endl;
				}
			}
		}
		else
		{
			// Loop through hills and calculate the bias force.
			for(auto& hill : hills_)
			{		
				auto tbias = 1.;
				std::fill(tder_.begin(), tder_.end(), 1.0);
				std::fill(dx_.begin(), dx_.end(), 1.0);
				
				for(size_t i = 0; i < n; ++i)
				{
					dx_[i] = cvs[i]->GetDifference(hill.center[i]);
					tbias *= gaussian(dx_[i], hill.width[i]);
				}

				for(size_t i = 0; i < n; ++i)
					for(size_t j = 0; j < n; ++j)
					{
						if(j != i) 
							tder_[i] *= gaussian(dx_[j], hill.width[j]);
						else
							tder_[i] *= gaussianDerv(dx_[j], hill.width[j]);
					}

				bias += height_ * tbias;
				for(size_t i = 0; i < n; ++i)
					derivatives_[i] += height_*tder_[i];
			}
		}

		// Restraints.
		for(size_t i = 0; i < n; ++i)
		{
			auto cval = cvs[i]->GetValue();
			if(cval < lowerb_[i])
				derivatives_[i] += lowerk_[i]*(cval - lowerb_[i]);
			else if(cval > upperb_[i])
				derivatives_[i] += upperk_[i]*(cval - upperb_[i]);
		}
	}

	Meta* Meta::Build(const Json::Value& json, 
		                  const MPI_Comm& world,
		                  const MPI_Comm& comm,
					      const std::string& path)
	{
		ObjectRequirement validator;
		Value schema;
		CharReaderBuilder rbuilder;
		CharReader* reader = rbuilder.newCharReader();

		reader->parse(JsonSchema::MetadynamicsMethod.c_str(),
		              JsonSchema::MetadynamicsMethod.c_str() + JsonSchema::MetadynamicsMethod.size(),
		              &schema, NULL);
		validator.Parse(schema, path);

		// Validate inputs.
		validator.Validate(json, path);
		if(validator.HasErrors())
			throw BuildException(validator.GetErrors());

		std::vector<double> widths;
		for(auto& s : json["widths"])
			widths.push_back(s.asDouble());

		std::vector<double> lowerb, upperb, lowerk, upperk;
		Grid<Vector>* grid = nullptr;
		if(json.isMember("grid"))
			grid = Grid<Vector>::BuildGrid(json.get("grid", Json::Value()));
		else if(!json.isMember("lower_bounds") || !json.isMember("upper_bounds"))
			throw BuildException({
				"#/Method/Metadynamics: Both upper_bounds and lower_bounds "
				"must be defined if grid is not being used."});

		// Assume all vectors are the same size. 
		for(int i = 0; i < static_cast<int>(json["lower_bound_restraints"].size()); ++i)
		{
			lowerk.push_back(json["lower_bound_restraints"][i].asDouble());
			upperk.push_back(json["upper_bound_restraints"][i].asDouble());
			lowerb.push_back(json["lower_bounds"][i].asDouble());
			upperb.push_back(json["upper_bounds"][i].asDouble());
		}
	
		auto height = json.get("height", 1.0).asDouble();
		auto hillfreq = json.get("hill_frequency", 1).asInt();
		auto freq = json.get("frequency", 1).asInt();

		auto* m = new Meta(
			world, comm, height, widths, 
			lowerb, upperb, lowerk,	upperk,
			grid, hillfreq, freq
		);

		if(json.isMember("load_hills"))
			m->LoadHills(json["load_hills"].asString());

		return m;		
	}
}
