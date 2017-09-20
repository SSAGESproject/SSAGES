#include "ANN.h"
#include "mxx/bcast.hpp"
#include "Snapshot.h"
#include "CVs/CVManager.h"

using namespace Eigen; 
using namespace nnet;

namespace SSAGES
{
	ANN::ANN(const MPI_Comm& world, 
		     const MPI_Comm& comm, 
		     const VectorXi& topol,
		     Grid<VectorXd>* fgrid,
		     const std::vector<double>& lowerb,
		     const std::vector<double>& upperb,
		     const std::vector<double>& lowerk,
		     const std::vector<double>& upperk,
		     double temperature,
		     double weight,
		     uint maxiter, 
			 uint sweep) : 
	Method(1, world, comm), topol_(topol), maxiter_(maxiter), sweep_(0),
	nsweep_(sweep),  net_(topol), pweight_(1.), weight_(weight), temp_(temperature), 
	kbt_(0), fgrid_(fgrid), hgrid_(nullptr), hist_(), bias_(), uhist_(),
	lowerb_(lowerb), upperb_(upperb), lowerk_(lowerk), upperk_(upperk)
	{
		// Build histogram grid. 
		hgrid_ = new Grid<uint>(
			fgrid_->GetNumPoints(), 
			fgrid_->GetLower(), 
			fgrid_->GetUpper(), 
			fgrid_->GetPeriodic()
		);

		// Create histogram grid matrix.
		auto points = hgrid_->GetNumPoints();
		auto ntot = std::accumulate(points.begin(), points.end(), 1, std::multiplies<int>());
		hist_.resize(ntot, hgrid_->GetDimension());

		// Fill it up. 
		size_t i = 0;
		for(auto it = hgrid_->begin(); it != hgrid_->end(); ++it)
		{
			auto& val = *it; 
			auto coord = it.coordinates(); 
			for(size_t j = 0; j < coord.size(); ++j)
				hist_(i, j) = coord[j]; 
			++i;
		}

		// Initialize FES vector.
		bias_.resize(ntot, 1);
		uhist_.resize(ntot, 1);
		
		for(i = 0; i < ntot; ++i)
		{
			bias_(i) = 0;
			uhist_(i) = 1;
		}
	}

	void ANN::PreSimulation(Snapshot* snapshot, const CVManager&) 
	{
		auto ndim = hgrid_->GetDimension();
		kbt_ = snapshot->GetKb()*temp_;
		
		// Zero out forces and histogram. 
		VectorXd vec = VectorXd::Zero(ndim);
		std::fill(hgrid_->begin(), hgrid_->end(), 0);
		std::fill(fgrid_->begin(), fgrid_->end(), vec);
	}

	void ANN::PostIntegration(Snapshot* snapshot, const CVManager& cvmanager)
	{
		if(snapshot->GetIteration() && snapshot->GetIteration() % nsweep_ == 0)
		{
			TrainNetwork();
			if(world_.rank() == 0)
				WriteBias();
		}

		// Get CV vals.
		auto cvs = cvmanager.GetCVs(cvmask_);
		auto n = cvs.size();

		// Determine if we are in bounds.
		RowVectorXd vec(n);
		std::vector<double> val(n);
		bool inbounds = true;
		for(size_t i = 0; i < n; ++i)
		{
			val[i] = cvs[i]->GetValue();
			vec[i] = cvs[i]->GetValue();
			if(val[i] < hgrid_->GetLower(i) || val[i] > hgrid_->GetUpper(i))
				inbounds = false;
		}

		// If in bounds, bias. 
		VectorXd derivatives = VectorXd::Zero(n);
		if(inbounds)
		{
			// Record histogram hit and get gradient. 
			hgrid_->at(val) += 1;
			//derivatives = (*fgrid_)[val];
			net_.forward_pass(vec);
			derivatives = net_.get_gradient(0);
		}
		else
		{   
			if(comm_.rank() == 0)
			{
				std::cerr << "ANN (" << snapshot->GetIteration() << "): out of bounds ( ";
				for(auto& v : val)
					std::cerr << v << " "; 
				std::cerr << ")" << std::endl;
			}
		}

		// Restraints.
		for(size_t i = 0; i < n; ++i)
		{
			auto cval = cvs[i]->GetValue();
			if(cval < lowerb_[i])
				derivatives[i] += lowerk_[i]*cvs[i]->GetDifference(cval - lowerb_[i]);
			else if(cval > upperb_[i])
				derivatives[i] += upperk_[i]*cvs[i]->GetDifference(cval - upperb_[i]);
		}

		// Apply bias to atoms. 
		auto& forces = snapshot->GetForces(); 
		auto& virial = snapshot->GetVirial();

		for(size_t i = 0; i < cvs.size(); ++i)
		{
			auto& grad = cvs[i]->GetGradient();
			auto& boxgrad = cvs[i]->GetBoxGradient();
			
			// Update the forces in snapshot by adding in the force bias from each
			// CV to each atom based on the gradient of the CV.
			for (size_t j = 0; j < forces.size(); ++j)
				forces[j] -= derivatives[i]*grad[j];
			
			virial += derivatives[i]*boxgrad;
		}
	}

	void ANN::PostSimulation(Snapshot*, const CVManager&) 
	{
	}

	void ANN::TrainNetwork()
	{
		// Increment cycle counter.
		++sweep_;

		// Synchronize grid in case it's periodic. 
		hgrid_->syncGrid();

		// Update FES estimator.
		Map<Array<uint, Dynamic, 1>> hist(hgrid_->data(), hgrid_->size());
		uhist_.array() = pweight_*uhist_.array() + hist.cast<double>()*(1./kbt_*bias_).array().exp()*weight_/sweep_;
		hist = 0;

		bias_.array() = kbt_*uhist_.array().log();
		bias_.array() -= bias_.minCoeff();
		
		// Train network.
		if(world_.rank() == 0)
		{
			net_.autoscale(hist_, bias_);
			net_.train(hist_, bias_, true);
		}

		// Send optimal nnet params to all procs.
		vector_t wb = net_.get_wb();
		mxx::bcast(wb.data(), wb.size(), 0, world_);
		net_.set_wb(wb);

		// Evaluate and subtract off min value for applied bias.
		net_.forward_pass(hist_);
		bias_.array() = net_.get_activation().col(0).array();
		bias_.array() -= bias_.minCoeff();

		// Calc new bias force.
		for(size_t i = 0; i < fgrid_->size(); ++i)
		{
			MatrixXd forces = net_.get_gradient(i); 
			fgrid_->data()[i] = forces.row(i).transpose();
		}
	}
}