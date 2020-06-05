/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2018 Emre Sevgen <sesevgen@uchicago.edu> and Hythem Sidky <hsidky@nd.edu>
 *           2020 Elizabeth M.Y. Lee <emlee@uchicago.edu> and Boyuan Yu <boyuanyu@uchicago.edu>
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
#include "CFF.h"
#include "schema.h"
#include "Snapshot.h"
#include "mxx/bcast.hpp"
#include "CVs/CVManager.h"
#include "Drivers/DriverException.h"
#include "Validator/ObjectRequirement.h"
#include <ctime>

using namespace Json;

// Implementation of the Combined Force Frequency (CFF) method from
// Ref: Sevgen et al. J. Chem. Theory Comput. 2020, 16, 3, 1448-1455.a

namespace SSAGES
{
	//! Pre-simulation hook.
	/*!
	 * Initialize biasing forces and histogram.
	 */
	void CFF::PreSimulation(Snapshot* snapshot, const CVManager& cvmanager)
	{
		auto cvs = cvmanager.GetCVs(cvmask_);
		dim_ = cvs.size();

		// Size and initialize Fold_.
		Fold_.setZero(dim_);

		// Initialize w \dot p's for finite difference.
		wdotp1_.setZero(dim_);
		wdotp2_.setZero(dim_);

		int ndim = hgrid_->GetDimension();
		kbt_ = snapshot->GetKb()*temp_;

		// Zero out forces and histogram.
		Eigen::VectorXd vec = Eigen::VectorXd::Zero(ndim);
		std::fill(hgrid_->begin(), hgrid_->end(), 0);
		std::fill(ugrid_->begin(), ugrid_->end(), 1.0);
		std::fill(fgrid_->begin(), fgrid_->end(), vec);

		// Initialize Nworld.
		Eigen::Map<Eigen::Array<unsigned int, Eigen::Dynamic, 1>> hist(hgrid_->data(), hgrid_->size());
		Nworld_ = hist.cast<int>();

		// Scale initial bias by 1/2*kT and make them positive.
		bias_.array() = abs(bias_.array())*kbt_*0.5;
	 }

	//! Post-integration hook.
	/*!
	 * First, information from current snapshot is retrieved and stored in v-integration hook.
	 * Then, coordinates in CV space are determined.
	 * Then, for each CV, biasing force is calculated.
	 * Then, neural networks are trained if called.
	 * Then, information is printed out if called.
	 * Finally, bias is applied using either genralized force (for the first sweep) or neural networks.
	 */
	void CFF::PostIntegration(Snapshot* snapshot, const CVManager& cvmanager)
	{

		// Gather information.
		// // Bias energy is kb*T*ln(P) where P = unbiased distribution
		auto cvs = cvmanager.GetCVs(cvmask_);
		auto& vels = snapshot->GetVelocities();
		auto& mass = snapshot->GetMasses();
		auto& forces = snapshot->GetForces();
		auto& virial = snapshot->GetVirial();
		auto n = snapshot->GetNumAtoms();

		using NAtoms = decltype(n);

		if(snapshot->GetIteration() && snapshot->GetIteration() % nsweep_ == 0 && snapshot->GetIteration() >= nsweep_*1 )
		{
			// Switch to full blast.
			if(citers_ && snapshot->GetIteration() > citers_)
				pweight_ = 1.0;

			TrainNetwork();
			if(IsMasterRank(world_))
				WriteBias();
		}

		// Determine if we are in bounds.
		Eigen::RowVectorXd vec(dim_);
		std::vector<double> val(dim_);
		bool inbounds = true;
		for(int i = 0; i < dim_; ++i)
		{
			val[i] = cvs[i]->GetValue();
			vec[i] = cvs[i]->GetValue();
			if(val[i] < hgrid_->GetLower(i) || val[i] > hgrid_->GetUpper(i))
				inbounds = false;
		}

		// Initialize matrix to hold the CV gradient.
		Eigen::MatrixXd J(dim_, 3*n);

		// Fill J and CV. Each column represents grad(CV) with flattened Cartesian elements.
		for(int i = 0; i < dim_; ++i)
		{
			auto& grad = cvs[i]->GetGradient();
			for(NAtoms j = 0; j < n; ++j)
				J.block<1, 3>(i,3*j) = grad[j];
		}

		//* Calculate W using Darve's approach (http://mc.stanford.edu/cgi-bin/images/0/06/Darve_2008.pdf).
		// However, we will not use mass weighing.
		Eigen::MatrixXd Jmass = J.transpose();

		Eigen::MatrixXd Minv = J*Jmass;
		MPI_Allreduce(MPI_IN_PLACE, Minv.data(), Minv.size(), MPI_DOUBLE, MPI_SUM, comm_);
		Eigen::MatrixXd Wt = Minv.inverse()*Jmass.transpose();

		// Fill momenta.
		Eigen::VectorXd momenta(3*n);
		for(NAtoms i = 0; i < n; ++i)
			momenta.segment<3>(3*i) = mass[i]*vels[i];

		// Compute dot(w,p)
		Eigen::VectorXd wdotp = Wt*momenta;

		// Reduce dot product across processors.
		MPI_Allreduce(MPI_IN_PLACE, wdotp.data(), wdotp.size(), MPI_DOUBLE, MPI_SUM, comm_);

		// Compute d(wdotp)/dt second order backwards finite difference.
		// Adding old force removes bias.
		Eigen::VectorXd dwdotpdt = unitconv_*(1.5*wdotp - 2.0*wdotp1_ + 0.5*wdotp2_)/timestep_ + Fold_;
		Eigen::VectorXd derivatives = Eigen::VectorXd::Zero(dim_);

		// If we are in bounds, sum force and frequency into running total.
		if (inbounds)
		{
			if(IsMasterRank(comm_))
			{
				for(int i=0; i<dim_; ++i)
					F_[i]->at(val) += dwdotpdt[i];

				hgrid_->at(val) += 1;
			}

			// Initial sweep is the same as doing adaptive basing force
			// i.e., Calculate biasing forces from averaged F at current CV coodinates

			if(snapshot->GetIteration() < nsweep_)
			{
				for(int i = 0; i < dim_; ++i)
					derivatives[i] = (F_[i]->at(val)/std::max((double(hgrid_->at(val))),double(min_)));
			}
			else
			{
				net_.forward_pass(vec);
				net2_.forward_pass(vec);
				derivatives = net_.get_gradient(0)*ratio_ + net2_.get_gradient(0)*(1.0-ratio_);
			}

		}
		// If out of bounds, apply harmonic restraint
		else
		{
			// Output to screen CV value that is out of bounds
			if(IsMasterRank(comm_))
			{
				std::cerr << "CFF (" << snapshot->GetIteration() << "): out of bounds ( ";
				for(auto& v : val)
					std::cerr << v << " ";
				std::cerr << ")" << std::endl;
			}

			// Apply harmonic restraints
			for(int i = 0; i < dim_; ++i)
			{
				auto cval = cvs[i]->GetValue();
				if(cval < lowerb_[i])
					derivatives[i] += lowerk_[i]*cvs[i]->GetDifference(cval - lowerb_[i]);
				else if(cval > upperb_[i])
					derivatives[i] += upperk_[i]*cvs[i]->GetDifference(cval - upperb_[i]);
			}
		}


		for(int i = 0; i < dim_; ++i)
		{
			auto& grad = cvs[i]->GetGradient();
			auto& boxgrad = cvs[i]->GetBoxGradient();

			// Update the forces in snapshot by adding in the force bias from each
			// CV to each atom based on the gradient of the CV.
			for (NAtoms j = 0; j < n; ++j)
				forces[j] -= derivatives[i]*grad[j];

			// Update virial.
			virial += derivatives[i]*boxgrad;
		}
		// Collect all walkers.
		MPI_Barrier(world_);

		// Store the old summed forces.
		Fold_ = derivatives;

		// Update finite difference time derivatives.
		wdotp2_ = wdotp1_;
		wdotp1_ = wdotp;
	}

	//! Post-simulation hook.
	void CFF::PostSimulation(Snapshot*, const CVManager&)
	{
	}

	//! Train neural networks.
	void CFF::TrainNetwork()
	{
		// Start the clock to measure the neural network training time.
		std::clock_t start = std::clock();

		// Increment cycle counter.
		++sweep_;

		// Reduce histogram across procs and sync in case system is periodic.
		mxx::allreduce(hgrid_->data(), hgrid_->size(), std::plus<unsigned int>(), world_);
		hgrid_->syncGrid();

		// Update visit frequencies.
		Eigen::Map<Eigen::Array<unsigned int, Eigen::Dynamic, 1>> hist(hgrid_->data(), hgrid_->size());
		Nworld_ += hist.cast<int>();

		// Synchronize unbiased histogram.
		ugrid_->syncGrid();
		Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, 1>> uhist(ugrid_->data(), ugrid_->size());

		// Update average biased forces across all processors
		std::vector<Eigen::MatrixXd> Ftrain;
		for(int i=0; i<dim_; ++i)
		{
			MPI_Allreduce(F_[i]->data(), Fworld_[i]->data(), (F_[i]->size()), MPI_DOUBLE, MPI_SUM, world_);
			Fworld_[i]->syncGrid();

			// Damp forces used to train net2_ via minimum number of hits
			Eigen::ArrayXd Nmin = Eigen::ArrayXd::Ones(Fworld_[0]->size())*min_;
			Ftrain.push_back(Eigen::Map<Eigen::Array<double, Eigen::Dynamic, 1>> (Fworld_[i]->data(),Fworld_[i]->size()) /
				( (Nworld_.cast<double>()).max(Nmin) ) );
		}

		// Calculate unbiased histrogram from previous unbiased histogram plus estimates from bias energy.
		uhist.array() = pweight_*uhist.array() + hist.cast<double>()*(bias_/kbt_).array().exp()*weight_;

		// Synchronize unbiased histogram and clear global histogram holder.
		ugrid_->syncGrid();
		hist.setZero();

		// Initialize boolean vector to enable training data.
		// 1 = include specified CV bin for training; 0 = remove from training.
		// Useful for training data partially in the future.
		force_to_val_ratio_ = Eigen::MatrixXd::Zero(hist_.rows(),1);

		// Bias energy is kb*T*ln(P) where P = unbiased distribution.
		bias_.array() = kbt_*uhist.array().log();
		bias_.array() -= bias_.minCoeff();

		if (ratio_ > 0.6)
		{
			net2_.init_weights();
		}

		// Train network.
		net_.autoscale(hist_, bias_);
		net2_.autoscale_w_grad(hist_, bias_, Ftrain);
		if(IsMasterRank(world_))
		{
			SetRatio(1.0);
			double gamma1;
			gamma1 = net_.train(hist_, bias_, true);
			SetRatio(0.0);
			double gamma2 = net2_.train_w_grad(hist_, bias_,Ftrain, force_to_val_ratio_, true);
			std::cout << "gamma1 " << gamma1 << " " << gamma2 << std::endl;
			ratio_ = gamma1/(gamma1+gamma2);

			std::cout << std::endl << "Ratio: " << ratio_ << std::endl;

			// Output file that accumulates information on the value of ratio_ for every sweeps
			std::ofstream gammaout;
			gammaout.open(outfile_ + "_gamma", std::ios_base::app);
			gammaout.precision(16);
			gammaout << sweep_ << " " << gamma1 << " " << gamma2 << " " << ratio_ << std::endl;
			gammaout.close();
		}

		// Send optimal neural net params to all procs.
		nnet::vector_t wb = net_.get_wb();
		mxx::bcast(wb.data(), wb.size(), 0, world_);
		net_.set_wb(wb);

		nnet::vector_t wb2 = net2_.get_wb();
		mxx::bcast(wb2.data(), wb2.size(), 0, world_);
		net2_.set_wb(wb2);

		mxx::bcast(&ratio_, 1, 0, world_);

		// Evaluate and subtract off min value for applied bias.
		net_.forward_pass(hist_);
		net2_.forward_pass(hist_);
		bias_.array() = net_.get_activation().col(0).array() * ratio_ + net2_.get_activation().col(0).array() * (1.0-ratio_);
		bias_.array() -= bias_.minCoeff();

		// Average the generalized force for each bin and output the file.
		if(IsMasterRank(world_))
		{
			int gridPoints = 1;
			for(int i = 0 ; i < dim_; ++i)
				gridPoints = gridPoints * hgrid_->GetNumPoints(i);

			std::string filename;
			filename = overwrite_ ? "F_out" : "F_out_"+std::to_string(sweep_);
			std::ofstream file(filename);
			file << std::endl;
			file << "Sweep: " << sweep_ << std::endl;
			file << "Printing out the current Biasing Vector Field from CFF." << std::endl;
			file << "First (Nr of CVs) columns are the coordinates, the next (Nr of CVs) columns are components of the Generalized Force vector at that point." << std::endl;
			file << "The columns are " << gridPoints << " long, mapping out a surface of ";
			for (int i = 0 ; i < dim_-1; ++i)
				file << (hgrid_->GetNumPoints())[i] << " by " ;
			file << (hgrid_->GetNumPoints(dim_-1)) << " points in " << dim_ << " dimensions." << std::endl;
			file << std::endl;

			for(int j=0; j < (Ftrain[0].array()).size();++j)
			{
			    file << std::setw(14) << std::setprecision(8) << std::fixed << hist_.row(j);
			    for(int i = 0 ; i < dim_; ++i){
					nnet::matrix_t x = Ftrain[i].array();
			   	    file << std::setw(14) << std::setprecision(8) << std::fixed << x(j) << " ";
			    }
			    file << std::endl;
			}
			file << std::endl;
			file << std::endl;
			file.close();
		}

		// Output training time information
		double duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
		std::ofstream file1("traintime.out",std::ofstream::app);
		file1 << duration << std::endl;
		file1.close();

	}

	// Write out neural network data
	void CFF::WriteBias()
	{
		// Write neural network topology and parameters
		std::string filename;
		filename = overwrite_ ? "netstate.dat" : "netstate_"+std::to_string(sweep_)+".dat";
		net_.write(filename.c_str());
		filename = overwrite_ ? "netstate2.dat" : "netstate2_"+std::to_string(sweep_)+".dat";
		net2_.write(filename.c_str());


		// Write CFF output
		filename = overwrite_ ? outfile_ : outfile_ + std::to_string(sweep_);
		std::ofstream file(filename);
		file.precision(8);
		net_.forward_pass(hist_);
		net2_.forward_pass(hist_);
		nnet::matrix_t x = net_.get_activation().array();
		nnet::matrix_t y = net2_.get_activation().array();
		nnet::matrix_t q = bias_;
		nnet::matrix_t m = -q;
		m.array() -= m.minCoeff();

		file << std::endl;
		file << "Sweep: " << sweep_ << std::endl;
		file << "Printing out the current Combined Force Frequency data." << std::endl;
		file << "First (Nr of CVs) columns are the coordinates." << std::endl;
		file << "Next columns (left to right) are: bias(freq_NN) bias(force_NN) bias(both_NN) free_energy" << std::endl;
		file << std::endl;
		for(int i = 0; i < y.rows(); ++i)
		{
			for(int j = 0; j < hist_.cols(); ++j)
				file << std::fixed << hist_(i,j) << " ";
			file << std::fixed<< x(i)<< " " << std::fixed << y(i)<< " " << std::fixed << q(i) <<" " << std::fixed << m(i) << "\n";
		}
		file.close();

	}

	CFF* CFF::Build(
		const Json::Value& json,
		const MPI_Comm& world,
		const MPI_Comm& comm,
		const std::string& path)
	{
		ObjectRequirement validator;
		Value schema;
		CharReaderBuilder rbuilder;
		CharReader* reader = rbuilder.newCharReader();
		reader->parse(
			JsonSchema::CFFMethod.c_str(),
			JsonSchema::CFFMethod.c_str() + JsonSchema::CFFMethod.size(),
			&schema,
			nullptr
		);
		validator.Parse(schema, path);

		// Validate inputs.
		validator.Validate(json, path);
		if(validator.HasErrors())
			throw BuildException(validator.GetErrors());

		// Grid.
		auto jsongrid = json.get("grid", Json::Value());
		auto* fgrid = Grid<Eigen::VectorXd>::BuildGrid(jsongrid);
		auto* hgrid = Grid<unsigned int>::BuildGrid(jsongrid);
		auto* ugrid = Grid<double>::BuildGrid(jsongrid);
		std::vector<Grid<double>*> F(json["grid"]["upper"].size());
		for(auto& grid : F)
			grid =	Grid<double>::BuildGrid(jsongrid);
		std::vector<Grid<double>*> Fworld(json["grid"]["upper"].size());
		for(auto& grid : Fworld)
			grid =	Grid<double>::BuildGrid(jsongrid);

		// Topology.
		auto nlayers = json["topology"].size() + 2;
		Eigen::VectorXi topol(nlayers);
		topol[0] = fgrid->GetDimension();
		topol[nlayers-1] = 1;
		for(decltype(nlayers) i = 0; i < json["topology"].size(); ++i)
			topol[i+1] = json["topology"][i].asInt();

		// CFF parameters
		auto weight = json.get("weight", 1.).asDouble();
		auto temp = json["temperature"].asDouble();
		auto nsweep = json["nsweep"].asUInt();
		auto unitconv = json.get("unit_conversion", 1).asDouble();
		auto timestep = json.get("timestep", 0.002).asDouble();
        auto min = json["minimum_count"].asInt();

		// Assume all vectors are the same size.
		std::vector<double> lowerb, upperb, lowerk, upperk;
		auto lbr_size = json["lower_bound_restraints"].size();
		for(decltype(lbr_size) i = 0; i < lbr_size; ++i)
		{
			lowerk.push_back(json["lower_bound_restraints"][i].asDouble());
			upperk.push_back(json["upper_bound_restraints"][i].asDouble());
			lowerb.push_back(json["lower_bounds"][i].asDouble());
			upperb.push_back(json["upper_bounds"][i].asDouble());
		}

		auto* m = new CFF(world, comm, topol, fgrid, hgrid, ugrid, F, Fworld, lowerb, upperb, lowerk, upperk, temp, unitconv, timestep, weight, nsweep, min);

		// Set optional params.
		m->SetPrevWeight(json.get("prev_weight", 1).asDouble());
		m->SetOutput(json.get("output_file", "CFF.out").asString());
		m->SetOutputOverwrite( json.get("overwrite_output", true).asBool());
		m->SetConvergeIters(json.get("converge_iters", 0).asUInt());
		m->SetMaxIters(json.get("max_iters", 1000).asUInt());
		m->SetMinLoss(json.get("min_loss", 0).asDouble());
		m->SetRatio(json.get("ratio", 0.0).asDouble());

		return m;
	}
}
