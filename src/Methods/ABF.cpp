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
	// Function to return histogram coordinates, given CV values. Returns -1 if any CV value is outside bounds.
	// Regardless of number of CVs, the histogram is a one dimensional object.
	// Therefore, this function works in two parts:
	// 1) Mapping from CV coordinates to an N dimensional fictitious object where N is the number of CVs.
	// -> This fictitious object is (Number of CV bins) large in each dimension.
	// -> Example: If you have two CVs, X and Y, and you bin X into three bins and Y into five bins, the fictitious object is a 2 dimensional 3x5 matrix. 
	
	// 2) Mapping from N dimensions to one dimension.
	// -> Simply constructed by writing out all members of the fictitious object in a line. This is a vector of length equal to the product of number of bins in each dimension.
	// -> For the above example, this is 3.5 = 15 members. First 5 members are Y = [1 2 3 4 5] bins with X = 1. Next 5 members are Y = [1 2 3 4 5] with X = 2 ....
	/*int ABF::histCoords(const CVList& cvs)
	{
		// Histogram details: This is a 2 Dimensional object set up as the following:
		// Histdetails is a vector of three vectors, each of those three vectors are (Number of CVs) long.
		// First of these vectors hold the lower bound for the CVs in order.
		// Second vector holds the upper bound for the CVs in order.
		// Third vector holds number of histogram bins for the CVs in order.
	
		// The first mapping starts here.
		for(size_t i = 0; i < histdetails_.size(); ++i)
		{
			// Check if CV is in bounds for each CV.
			if ((cvs[i]->GetValue()<histdetails_[i][0]) || (cvs[i]->GetValue()>histdetails_[i][1]))
				return -1;
		}

		// This vector holds the CV coordinates in the fictitious object.
		std::vector<int> coords(cvs.size());
		
		// Loop over each CV dimension.
		for(size_t i = 0; i < histdetails_.size(); ++i)
		{
			// Loop over all bins in each CV dimensions
			for(int j = 0; j < histdetails_[i][2] ; ++j)
			{
				if((histdetails_[i][0] + j*((histdetails_[i][1]-histdetails_[i][0])/histdetails_[i][2]) < cvs[i]->GetValue()) && 
				   (histdetails_[i][0] + (j+1)*((histdetails_[i][1]-histdetails_[i][0])/histdetails_[i][2]) > cvs[i]->GetValue()))
				{
					coords[i] = j;
				}
			}
		} //First mapping ends here. 

		// Second mapping starts here.
		int finalcoord = 0;
		int indexer = 1;
		// Loop over each CV dimension.
		for(size_t i = 0; i < histdetails_.size(); ++i)
		{
			indexer = 1;
			for(size_t j = i+1; j < histdetails_.size(); ++j)
				indexer = indexer*histdetails_[j][2];

			finalcoord += coords[i]*indexer;
		}

		return finalcoord;			
	}
	*/
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
		for(auto& point : (*F_)){
			point.setZero(dim_);
			}
			
		for(auto& point : (*Fworld_)){
			point.setZero(dim_);
			}

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
			F_->at(cvVals) += dwdotpdt;
			N_->at(cvVals)++;				
		}
		
		// Reduce data across processors.					
		MPI_Allreduce(F_->data(), Fworld_->data(), (F_->size())*dim_, MPI_DOUBLE, MPI_SUM, world_);	
		MPI_Allreduce(N_->data(), Nworld_->data(), N_->size(), MPI_INT, MPI_SUM, world_);
		

		// If we are in bounds, store the old summed force.
		if(boundsCheck(cvVals))
		{
			Fold_ = Fworld_->at(cvVals)/std::max(min_, Nworld_->at(cvVals));
		}
		
	
		// Update finite difference time derivatives.
		wdotp2_ = wdotp1_;
		wdotp1_ = wdotp;		

		// Write out data to file.
		if(iteration_ % FBackupInterv_ == 0)
			WriteData();
		
		// Calculate the bias from averaged F at current CV coordinates
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
		
		// Compute bias if within bounds
		if(boundsCheck(cvVals))
		{
			for(int i = 0; i < dim_; ++i)
			{
				auto& grad = cvs[i]->GetGradient();
				for(size_t j = 0; j < biases_.size(); ++j)
					biases_[j] -= (Fworld_->at(cvVals))[i]*grad[j]/std::max(min_,Nworld_->at(cvVals));
			}
		}
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

	void ABF::WriteData()
	{
		if(world_.rank() != 0)
			return;
			
		int gridPoints = 1;
			
		for(size_t i = 0 ; i < dim_; ++i)
			gridPoints = gridPoints * N_->GetNumPoints(i);

		worldout_.open(filename_, std::ofstream::app);
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
				worldout_ <<  std::setw(10) << (Fworld_->at(tempcoord))[j]/std::max(Nworld_->at(tempcoord),min_)<< " ";
			}
			worldout_ << std::endl;
			
		}

		worldout_ << std::endl;
		worldout_ << std::endl;
		worldout_.close();
	}
	
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

	void ABF::printFworld()
	{
		Fworldout_.open(Fworld_filename_);
		for(auto& point : (*F_))
		{
			Fworldout_ << point << std::endl;
		}
		Fworldout_.close();
	}

	ABF* ABF::Construct(const Value& json, 
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

		std::vector<double> minsCV;
		for(auto& mins : json["CV_lower_bounds"])
			minsCV.push_back(mins.asDouble());
		
		std::vector<double> maxsCV;
		for(auto& maxs : json["CV_upper_bounds"])
			maxsCV.push_back(maxs.asDouble());

		std::vector<int> binsCV;
		for(auto& bins : json["CV_bins"])
			binsCV.push_back(bins.asInt());

		std::vector<double> minsrestCV;
		for(auto& mins : json["CV_restraint_minimums"])
			minsrestCV.push_back(mins.asDouble());
		
		std::vector<double> maxsrestCV;
		for(auto& maxs : json["CV_restraint_maximums"])
			maxsrestCV.push_back(maxs.asDouble());

		std::vector<double> springkrestCV;
		for(auto& bins : json["CV_restraint_spring_constants"])
			springkrestCV.push_back(bins.asDouble());

		std::vector<bool> isperiodic;
		for(auto& isperCV : json["CV_isperiodic"])
			isperiodic.push_back(isperCV.asBool());

		std::vector<double> minperboundaryCV;
		for(auto& minsperCV : json["CV_periodic_boundary_lower_bounds"])
			minperboundaryCV.push_back(minsperCV.asDouble());

		std::vector<double> maxperboundaryCV;
		for(auto& maxsperCV : json["CV_periodic_boundary_upper_bounds"])
			maxperboundaryCV.push_back(maxsperCV.asDouble());

		if(!(minsCV.size() == maxsCV.size() && 
			 maxsCV.size() == binsCV.size() &&
			 binsCV.size() == minsrestCV.size() &&
			 minsrestCV.size() == maxsrestCV.size() &&
			 maxsrestCV.size() == springkrestCV.size() &&
			 springkrestCV.size() == isperiodic.size()))			
			throw BuildException({"CV lower bounds, upper bounds, bins, restrain minimums, restrains maximums, spring constants and periodicity info must all have the size == number of CVs defined."});

		bool anyperiodic=false;
		for(size_t i = 0; i<isperiodic.size(); ++i)
			if(isperiodic[i])
			{
				anyperiodic = true;
				if(!(isperiodic.size() == minperboundaryCV.size() &&
				     minperboundaryCV.size() == maxperboundaryCV.size()))
					throw BuildException({"If any CV is defined as periodic, please define the full upper and lower bound vectors. They should both have the same number of entries as CV lower bounds, upper bounds... Entries corresponding to non-periodic CVs will not be used."});
			}
				
		auto FBackupInterv = json.get("backupF_requency", 1000).asInt();
		auto unitconv = json.get("unit_conversion", 0).asDouble();
		auto timestep = json.get("timestep",2).asDouble();
		auto min = json.get("minimum_count",100).asDouble();
		auto massweigh = json.get("mass_weighing",false).asBool();			

		std::vector<std::vector<double>> histdetails;
		std::vector<std::vector<double>> restraint;
		std::vector<std::vector<double>> periodicboundaries;
		std::vector<double> temp1(3);
		std::vector<double> temp2(3);
		std::vector<double> temp3(2);

		for(size_t i=0; i<minsCV.size(); ++i)
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

		auto freq = json.get("frequency", 1).asInt();
		auto filename = json.get("filename", "F_out").asString();
		auto Fworld_filename_ = json.get("F_world_filename", "F_world").asString();
		auto Nworld_filename_ = json.get("N_world_filename", "N_world").asString();
		
		Grid<int> *N;
		Grid<Eigen::VectorXd> *F;

		N= new Grid<int>(binsCV, minsCV, maxsCV, isperiodic);		
		F= new Grid<Eigen::VectorXd>(binsCV, minsCV, maxsCV, isperiodic);

		Grid<int> *Nworld;
		Grid<Eigen::VectorXd> *Fworld;
		
		//Nworld= new Grid<uint>(binsCV, minsCV, maxsCV, isperiodic);
		//Fworld= new Grid<Eigen::VectorXd>(binsCV, minsCV, maxsCV, isperiodic);
	
		//if(std::ifstream(N_world_filename) && std::ifstream(F_world_filename))
		//{
		//	Nworld->LoadFromFile(Nworld_filename);
		//	Fworld->LoadFromFile(Fworld_filename);
		//}
		//else
		//{
			Nworld= new Grid<int>(binsCV, minsCV, maxsCV, isperiodic);
			Fworld= new Grid<Eigen::VectorXd>(binsCV, minsCV, maxsCV, isperiodic);
		//}
	 
		
		auto* m = new ABF(world, comm, N, Nworld, F, Fworld, restraint, isperiodic, periodicboundaries, min, massweigh, filename, histdetails, FBackupInterv, unitconv, timestep, freq);

		/*if(json.isMember("F") && json.isMember("N"))
		{
			Eigen::VectorXd F;
			std::vector<int> N;
			F.resize(json["F"].size());
			for(int i = 0; i < (int)json["F"].size(); ++i)
				F[i] = json["F"][i].asDouble();

			for(auto& n : json["N"])
				N.push_back(n.asInt());

			m->SetHistogram(F,N);
		
		}
		*/

		if(json.isMember("iteration"))
			m->SetIteration(json.get("iteration",0).asInt());

		return m;
		
	}

	void ABF::Serialize(Value& json) const
	{
		/*json["type"] = "ABF";
		for(size_t i = 0; i < histdetails_.size(); ++i)
		{
			json["CV_lower_bounds"].append(histdetails_[i][0]);				
			json["CV_upper_bounds"].append(histdetails_[i][1]);
			json["CV_bins"].append(histdetails_[i][2]);
		}

		for(size_t i = 0; i < restraint_.size(); ++i)
		{
			json["CV_restraint_minimums"].append(restraint_[i][0]);
			json["CV_restraint_maximums"].append(restraint_[i][1]);
			json["CV_restraint_spring_constants"].append(restraint_[i][2]);
		}

		json["timestep"] = timestep_;
		json["minimum_count"] = min_;
		json["backupF_requency"] = FBackupInterv_;			
		json["unit_conversion"] = unitconv_;
		json["iteration"] = iteration_;
		json["filename"] = filename_;		
			
		for(int i = 0; i < F_->size(); ++i)
			for(int j = 0; j < F_[i].size(); ++j)
				json["F"].append(F_[i][j]);

		for(size_t i = 0; i < N_->size(); ++i)
			json["N"].append(N_[i]);

		*/
	}
}



































