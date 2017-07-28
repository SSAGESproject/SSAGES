/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Emre Sevgen <sesevgen@uchicago.edu>
 *                Hythem Sidky <hsidky@nd.edu>
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
#include "ABF.h"
#include "CVs/CVManager.h"
#include "Drivers/DriverException.h"
#include "Validator/ObjectRequirement.h"
#include "Snapshot.h"
#include "schema.h"
#include <algorithm>
#include <fstream>
#include <iostream>

using namespace Json;
// "Adaptive biasing force method for scalar and vector free energy calculations"
// Darve, Rodriguez-Gomez, Pohorille
// J. Chem. Phys. (2008)
namespace SSAGES
{

	// Pre-simulation hook.
	void ABF::PreSimulation(Snapshot* snapshot, const CVManager& cvmanager)
	{
		
		// Open/close outfile to create it fresh. 
		if(world_.rank() == 0)
		{
			worldout_.open(filename_);
			worldout_.close();
		}
		
		// Convenience. Number of CVs.
		auto cvs = cvmanager.GetCVs(cvmask_);
		dim_ = cvs.size();

		// Size and initialize Fold_
		Fold_.setZero(dim_);
		
		// Size and initialize Fworld_ and Nworld_		
		auto nel = 1;

		// Initialize biases.
		biases_.resize(snapshot->GetPositions().size(), Vector3{0, 0, 0});
		
		// Initialize w \dot p's for finite difference. 
		wdotp1_.setZero(dim_);
		wdotp2_.setZero(dim_);	
	}

	//! Post-integration hook.
	/*!
	 * Post-integration hook is where processes carried out every timestep happen.
	 * First, information from current snapshot is retrieved and stored in variables as necessary.
	 * Then, coordinates in CV space are determined.
	 * Then, for each CV, the time derivative of w.p is calculated.
	 * Then, information is printed out if its enabled.
	 * Finally, bias is applied from current estimate of generalized force.
	 */
	void ABF::PostIntegration(Snapshot* snapshot, const CVManager& cvmanager)
	{
		
		++iteration_;

		// Gather information.
		auto cvs = cvmanager.GetCVs(cvmask_);
		auto& vels = snapshot->GetVelocities();
		auto& mass = snapshot->GetMasses();
		auto& forces = snapshot->GetForces();
		auto n = snapshot->GetNumAtoms();
		

		//! Coord holds where we are in CV space in current timestep.
		//int coord = histCoords(cvs);

		//! Eigen::MatrixXd to hold the CV gradient.
		Eigen::MatrixXd J(dim_, 3*n);
		
		//! std::vector<double> to hold CV values to address grid.
		std::vector<double> cvVals(dim_);

		// Fill J and CV. Each column represents grad(CV) with flattened Cartesian elements. 
		for(int i = 0; i < dim_; ++i)
		{
			auto& grad = cvs[i]->GetGradient();
			for(size_t j = 0; j < n; ++j)
				J.block<1, 3>(i,3*j) = grad[j];
			cvVals[i] = cvs[i]->GetValue();
		}

		//* Calculate W using Darve's approach (http://mc.stanford.edu/cgi-bin/images/0/06/Darve_2008.pdf).
		Eigen::MatrixXd Jmass = J.transpose();
		if(massweigh_)
		{
			for(size_t i = 0; i < forces.size(); ++i)
				Jmass.block(3*i, 0, 3, dim_) = Jmass.block(3*i, 0, 3, dim_)/mass[i];
		}

		Eigen::MatrixXd Minv = J*Jmass;
		MPI_Allreduce(MPI_IN_PLACE, Minv.data(), Minv.size(), MPI_DOUBLE, MPI_SUM, comm_);
		Eigen::MatrixXd Wt = Minv.inverse()*Jmass.transpose();	
		// Fill momenta.
		Eigen::VectorXd momenta(3*vels.size());
		for(size_t i = 0; i < vels.size(); ++i)
			momenta.segment<3>(3*i) = mass[i]*vels[i];
		
		// Compute dot(w,p)
		Eigen::VectorXd wdotp = Wt*momenta;

		
		// Reduce dot product across processors.
		MPI_Allreduce(MPI_IN_PLACE, wdotp.data(), wdotp.size(), MPI_DOUBLE, MPI_SUM, comm_);
		
	
		// Compute d(wdotp)/dt second order backwards finite difference. 
		// Adding old force removes bias. 
		Eigen::VectorXd dwdotpdt = unitconv_*(1.5*wdotp - 2.0*wdotp1_ + 0.5*wdotp2_)/timestep_ + Fold_;

		// If we are in bounds, sum force into running total.
		if(boundsCheck(cvVals))
		{
			for(size_t i=0; i<dim_; ++i)
				F_[i]->at(cvVals) += dwdotpdt[i];
			N_->at(cvVals)++;				
		}
		
		// Reduce data across processors.
		for(size_t i=0; i<dim_; ++i)				
			MPI_Allreduce(F_[i]->data(), Fworld_[i]->data(), (F_[i]->size()), MPI_DOUBLE, MPI_SUM, world_);	
		MPI_Allreduce(N_->data(), Nworld_->data(), N_->size(), MPI_INT, MPI_SUM, world_);

		// If we are in bounds, store the old summed force.
		if(boundsCheck(cvVals))
		{
			for(size_t i=0; i < dim_; ++i)
				Fold_[i] = Fworld_[i]->at(cvVals)/std::max(min_, Nworld_->at(cvVals));
		}
	
		// Update finite difference time derivatives.
		wdotp2_ = wdotp1_;
		wdotp1_ = wdotp;		

		// Write out data to file.
		if(iteration_ % FBackupInterv_ == 0)
			WriteData();

		// Calculate the bias from averaged F at current CV coordinates.
		// Or apply harmonic restraints to return CVs back in bounds.
		CalcBiasForce(snapshot, cvs, cvVals);		

		// Update the forces in snapshot by adding in the force bias from each
		// CV to each atom based on the gradient of the CV.
		for (size_t j = 0; j < forces.size(); ++j)
			forces[j] += biases_[j];	
	}

	// Post-simulation hook.
	void ABF::PostSimulation(Snapshot*, const CVManager&)
	{
		WriteData();
	}

	void ABF::CalcBiasForce(const Snapshot* snapshot, const CVList& cvs, const std::vector<double> &cvVals)
	{
		// Reset the bias.
		biases_.resize(snapshot->GetNumAtoms(), Vector3{0,0,0});
		for(auto& b : biases_)
			b.setZero();
		
		// Compute bias if within bounds.
		if(boundsCheck(cvVals))
		{
			for(int i = 0; i < dim_; ++i)
			{
				auto& grad = cvs[i]->GetGradient();
				for(size_t j = 0; j < biases_.size(); ++j)
					biases_[j] -= (Fworld_[i]->at(cvVals))*grad[j]/std::max(min_,Nworld_->at(cvVals));
			}
		}
		// Otherwise, apply harmonic restraint if enabled.
		else
		{
			for(int i = 0; i < dim_; ++i)
			{
				double cvVal = cvs[i]->GetValue();
				auto k = 0.;
				auto x0 = 0.;
				
				if(isperiodic_[i])
				{
					double periodsize = periodicboundaries_[i][1]-periodicboundaries_[i][0];
					double cvRestrMidpoint = (restraint_[i][1]+restraint_[i][0])/2;
					while((cvVal-cvRestrMidpoint) > periodsize/2)
						cvVal -= periodsize;
					while((cvVal-cvRestrMidpoint) < -periodsize/2)
						cvVal += periodsize;
				}

				if(cvVal < restraint_[i][0] && restraint_[i][2] > 0)
				{
					k = restraint_[i][2];
					x0 = restraint_[i][0];
				}
				else if (cvVal > restraint_[i][1] && restraint_[i][2] > 0)
				{
					k = restraint_[i][2];
					x0 = restraint_[i][1];
				}

				auto& grad = cvs[i]->GetGradient();
				for(size_t j = 0; j < biases_.size(); ++j)
					biases_[j] -= (cvVal - x0)*k*grad[j];
			}	
		}
	}

	// Write out the average generalized force.
	// Also write out Fworld and Nworld backups for restarts.
	void ABF::WriteData()
	{
		// Only one processor should be performing I/O.
		if(world_.rank() != 0)
			return;

		// Backup Fworld and Nworld.
		Nworld_->WriteToFile(Nworld_filename_);
		for(size_t i = 0 ; i < dim_; ++i)
		{
			Fworld_[i]->WriteToFile(Fworld_filename_+std::to_string(i));
		}
		

		// Average the generalized force for each bin and print it out.	
		int gridPoints = 1;			
		for(size_t i = 0 ; i < dim_; ++i)
			gridPoints = gridPoints * N_->GetNumPoints(i);

		worldout_.open(filename_,std::ios::out);
		worldout_ << std::endl;
		worldout_ << "Iteration: " << iteration_ << std::endl;			
		worldout_ << "Printing out the current Adaptive Biasing Vector Field." << std::endl;
		worldout_ << "First (Nr of CVs) columns are the coordinates, the next (Nr of CVs) columns are components of the Adaptive Force vector at that point." << std::endl;
		worldout_ << "The columns are " << gridPoints << " long, mapping out a surface of ";
		
		for(size_t i = 0 ; i < dim_-1; ++i)
			worldout_ << (N_->GetNumPoints())[i] << " by " ;
		worldout_ << (N_->GetNumPoints(dim_-1)) << " points in " << dim_ << " dimensions." << std::endl;
		worldout_ << std::endl;	

		std::vector<int> printCoords(dim_);
		int div = 1;
		int index = 0;
		std::vector<double> tempcoord(dim_);
		for(size_t i = 0; i < gridPoints; ++i)
		{
			printCoords[0] = i%(Nworld_->GetNumPoints(0));
			div = 1;
			for(size_t j = 1; j < dim_; ++j)
			{
				div = div * Nworld_->GetNumPoints(j);
				printCoords[j]=i/div;						
			}
			for(size_t j = 0; j < dim_; ++j)
			{
				worldout_ << std::setw(10) << (printCoords[j]+0.5)*((Nworld_->GetUpper(j)-Nworld_->GetLower(j))/Nworld_->GetNumPoints(j)) + Nworld_->GetLower(j) << " ";
				tempcoord[j] = ((printCoords[j]+0.5)*((Nworld_->GetUpper(j)-Nworld_->GetLower(j))/Nworld_->GetNumPoints(j)) + Nworld_->GetLower(j));
			}
			for(size_t j = 0; j < dim_; ++j)
			{
				worldout_ <<  std::setw(10) << (Fworld_[j]->at(tempcoord))/std::max(Nworld_->at(tempcoord),min_)<< " ";
			}
			worldout_ << std::endl;
			
		}

		worldout_ << std::endl;
		worldout_ << std::endl;
		worldout_.close();
	}
	
	// Checks whether walker is within CV bounds.
	bool ABF::boundsCheck(const std::vector<double> &CVs)
	{
		for(size_t i = 0; i < dim_ ; ++i)
		{
			if((CVs[i] < N_->GetLower(i)) || (CVs[i] > N_->GetUpper(i)))
			{
				return false;
			}
		}
		return true;
	}


	ABF* ABF::Build(const Value& json, 
		                const MPI_Comm& world,
		                const MPI_Comm& comm,
			            const std::string& path)
	{
		ObjectRequirement validator;
		Value schema;
		Reader reader;
		
		reader.parse(JsonSchema::ABFMethod, schema);
		validator.Parse(schema, path);

		// Validate inputs.
		validator.Validate(json, path);
		if(validator.HasErrors())
			throw BuildException(validator.GetErrors());

		uint wid = mxx::comm(world).rank()/mxx::comm(comm).size(); 

		//bool multiwalkerinput = false;

		// Obtain lower bound for each CV.
		std::vector<double> minsCV;
		for(auto& mins : json["CV_lower_bounds"])
			if(mins.isArray())
				{
				throw BuildException({"Separate inputs for multi-walker not fully implemented. Please use global entries for CV ranges"});
				//multiwalkerinput = true;
				//for(auto& bound : mins)
				//	minsCV.push_back(bound.asDouble());
				}
			else
				minsCV.push_back(mins.asDouble());
			
		// Obtain upper bound for each CV.
		std::vector<double> maxsCV;
		for(auto& maxs : json["CV_upper_bounds"])
			/*if(maxs.isArray())
				{
				for(auto& bound : maxs)
					maxsCV.push_back(bound.asDouble());
				}*/
			//else
				maxsCV.push_back(maxs.asDouble());

		// Obtain number of bins for each CV dimension.
		std::vector<int> binsCV;
		for(auto& bins : json["CV_bins"])
			binsCV.push_back(bins.asInt());

		// Obtain lower bounds for restraints for each CV.
		std::vector<double> minsrestCV;
		for(auto& mins : json["CV_restraint_minimums"])
			/*if(mins.isArray())
				{
				for(auto& bound : mins)
					minsrestCV.push_back(bound.asDouble());
				}*/
			//else
				minsrestCV.push_back(mins.asDouble());

		// Obtain upper bounds for restraints for each CV.		
		std::vector<double> maxsrestCV;
		for(auto& maxs : json["CV_restraint_maximums"])
			/*if(maxs.isArray())
				{
				for(auto& bound : maxs)
					maxsrestCV.push_back(bound.asDouble());
				}*/
			//else
				maxsrestCV.push_back(maxs.asDouble());
		
		// Obtain harmonic spring constant for restraints for each CV.
		std::vector<double> springkrestCV;
		for(auto& bins : json["CV_restraint_spring_constants"])
			springkrestCV.push_back(bins.asDouble());

		// Obtain periodicity information for restraints for each CV for the purpose of correctly applying restraints through periodic boundaries.
		std::vector<bool> isperiodic;
		for(auto& isperCV : json["CV_isperiodic"])
			isperiodic.push_back(isperCV.asBool());
		
		// Obtain lower periodic boundary for each CV.
		std::vector<double> minperboundaryCV;
		for(auto& minsperCV : json["CV_periodic_boundary_lower_bounds"])
			minperboundaryCV.push_back(minsperCV.asDouble());

		// Obtain upper periodic boundary for each CV.
		std::vector<double> maxperboundaryCV;
		for(auto& maxsperCV : json["CV_periodic_boundary_upper_bounds"])
			maxperboundaryCV.push_back(maxsperCV.asDouble());

		// Verify inputs are all correct lengths.
		if(!((	 minsCV.size() == maxsCV.size() && 
			 maxsCV.size() == minsrestCV.size() &&
			 minsrestCV.size() == maxsrestCV.size()) &&
			(binsCV.size() == springkrestCV.size() &&
			 springkrestCV.size() == isperiodic.size())))		
			throw BuildException({"CV lower bounds, upper bounds, restrain minimums, restrains maximums must match in size. Bins, spring constants and periodicity info must match in size."});

		// Verify that all periodicity information is provided if CVs are periodic.
		bool anyperiodic=false;
		for(size_t i = 0; i<isperiodic.size(); ++i)
			if(isperiodic[i])
			{
				anyperiodic = true;
				if(!(isperiodic.size() == minperboundaryCV.size() &&
				     minperboundaryCV.size() == maxperboundaryCV.size()))
					throw BuildException({"If any CV is defined as periodic, please define the full upper and lower bound vectors. They should both have the same number of entries as CV lower bounds, upper bounds... Entries corresponding to non-periodic CVs will not be used."});
			}

		int dim = binsCV.size();
		
		// Read in JSON info.
		auto FBackupInterv = json.get("output_frequency", 1000).asInt();
		auto unitconv = json.get("unit_conversion", 0).asDouble();
		auto timestep = json.get("timestep", 2).asDouble();
		auto min = json.get("minimum_count", 200).asDouble();
		auto massweigh = json.get("mass_weighing",false).asBool();			

		std::vector<std::vector<double>> histdetails;
		std::vector<std::vector<double>> restraint;
		std::vector<std::vector<double>> periodicboundaries;
		std::vector<double> temp1(3);
		std::vector<double> temp2(3);
		std::vector<double> temp3(2);

		auto freq = json.get("frequency", 1).asInt();
		auto filename = json.get("output_file", "F_out").asString();
		auto Nworld_filename = json.get("Nworld_output_file", "Nworld").asString();
		auto Fworld_filename = json.get("Fworld_output_file", "Fworld_cv").asString();

		// Generate the grids based on JSON.
		Grid<int> *N;
		Grid<int> *Nworld;
		std::vector<Grid<double>*> F(dim);
		std::vector<Grid<double>*> Fworld(dim);

		// This feature is disabled for now.
		/*if(multiwalkerinput)
		{
			for(size_t i=0; i<dim; ++i)
			{
				temp1 = {minsCV[i+wid*dim], maxsCV[i+wid*dim], binsCV[i]};
				temp2 = {minsrestCV[i+wid*dim], maxsrestCV[i+wid*dim], springkrestCV[i]};
				histdetails.push_back(temp1);
				restraint.push_back(temp2);
				if(anyperiodic)
				{
					temp3 = {minperboundaryCV[i], maxperboundaryCV[i]};
					periodicboundaries.push_back(temp3);
				}
			}
			for(size_t i=0; i<dim; ++i)
			{
				minsCVperwalker[i] = histdetails[i][0];
				maxsCVperwalker[i] = histdetails[i][1];
			}

			N= new Grid<int>(binsCV, minsCVperwalker, maxsCVperwalker, isperiodic);
			Nworld= new Grid<int>(binsCV, minsCVperwalker, maxsCVperwalker, isperiodic);
			
			for(auto& grid : F)
			{
				grid= new Grid<double>(binsCV, minsCVperwalker, maxsCVperwalker, isperiodic);
			}
			for(auto& grid : Fworld)
			{
				grid= new Grid<double>(binsCV, minsCVperwalker, maxsCVperwalker, isperiodic);
			}
			
		}
		*/
		//else
		//{
		
		// Appropriately size the grids.
		
			N= new Grid<int>(binsCV, minsCV, maxsCV, isperiodic);
			Nworld= new Grid<int>(binsCV, minsCV, maxsCV, isperiodic);
			
			for(auto& grid : F)
			{
				grid= new Grid<double>(binsCV, minsCV, maxsCV, isperiodic);
			}
			for(auto& grid : Fworld)
			{
				grid= new Grid<double>(binsCV, minsCV, maxsCV, isperiodic);
			}

			for(size_t i=0; i<dim; ++i)
			{
				temp1 = {minsCV[i], maxsCV[i], binsCV[i]};
				temp2 = {minsrestCV[i], maxsrestCV[i], springkrestCV[i]};
				histdetails.push_back(temp1);
				restraint.push_back(temp2);
				if(anyperiodic)
				{
					temp3 = {minperboundaryCV[i], maxperboundaryCV[i]};
					periodicboundaries.push_back(temp3);
				}
			}
		//}

		// Check if previously saved grids exist. If so, check that data match and load grids.
		if(std::ifstream(Nworld_filename) && std::ifstream(Fworld_filename+std::to_string(0)) && wid == 0)
		{
			std::cout << "Attempting to load data from a previous run of ABF." << std::endl;
			N->LoadFromFile(Nworld_filename);
			for(size_t i=0; i<dim; ++i)
				if(std::ifstream(Fworld_filename+std::to_string(i)))
					F[i]->LoadFromFile(Fworld_filename+std::to_string(i));
				else
					throw BuildException({"Some, but not all Fworld outputs were found. Please check that these are appropriate inputs, or clean the working directory of other Fworld and Nworld inputs."});
		}
	 
		
		auto* m = new ABF(world, comm, N, Nworld, F, Fworld, restraint, isperiodic, periodicboundaries, min, massweigh, filename, Nworld_filename, Fworld_filename, histdetails, FBackupInterv, unitconv, timestep, freq);


		if(json.isMember("iteration"))
			m->SetIteration(json.get("iteration",0).asInt());

		return m;
		
	}

}



































