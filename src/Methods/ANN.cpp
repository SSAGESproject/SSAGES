#include "ANN.h"
#include "schema.h"
#include "Snapshot.h"
#include "mxx/bcast.hpp"
#include "CVs/CVManager.h"
#include "Drivers/DriverException.h"
#include "Validator/ObjectRequirement.h"

using namespace Eigen; 
using namespace nnet;
using namespace Json;

namespace SSAGES
{
	ANN::ANN(const MPI_Comm& world, 
		     const MPI_Comm& comm, 
		     const VectorXi& topol,
			 Grid<VectorXd>* fgrid,
			 Grid<uint>* hgrid,
			 Grid<double>* ugrid,
		     const std::vector<double>& lowerb,
		     const std::vector<double>& upperb,
		     const std::vector<double>& lowerk,
		     const std::vector<double>& upperk,
		     double temperature,
		     double weight,
			 uint nsweep) : 
	Method(1, world, comm), topol_(topol), sweep_(0), nsweep_(nsweep),  
	citers_(0), net_(topol), pweight_(1.), weight_(weight), temp_(temperature), 
	kbt_(0), fgrid_(fgrid), hgrid_(hgrid), ugrid_(ugrid), hist_(), bias_(),
	lowerb_(lowerb), upperb_(upperb), lowerk_(lowerk), upperk_(upperk),
	outfile_("ann.out"), overwrite_(true)
	{
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
		for(i = 0; i < ntot; ++i)
			bias_(i) = 0;
	}

	void ANN::PreSimulation(Snapshot* snapshot, const CVManager&) 
	{
		auto ndim = hgrid_->GetDimension();
		kbt_ = snapshot->GetKb()*temp_;
		
		// Zero out forces and histogram. 
		VectorXd vec = VectorXd::Zero(ndim);
		std::fill(hgrid_->begin(), hgrid_->end(), 0);
		std::fill(ugrid_->begin(), ugrid_->end(), 1.);
		std::fill(fgrid_->begin(), fgrid_->end(), vec);
	}

	void ANN::PostIntegration(Snapshot* snapshot, const CVManager& cvmanager)
	{
		if(snapshot->GetIteration() && snapshot->GetIteration() % nsweep_ == 0)
		{
			// Switch to full blast.
			if(citers_ && snapshot->GetIteration() > citers_)
				pweight_ = 1.0;
			
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
			// Only record hits on master processes since we will 
			// reduce later. 
			if(comm_.rank() == 0)
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

		// Reduce histogram across procs. 
		mxx::allreduce(hgrid_->data(), hgrid_->size(), std::plus<uint>(), world_);

		// Synchronize grid in case it's periodic. 
		hgrid_->syncGrid();
	
		// Update FES estimator. Synchronize unbiased histogram.
		Map<Array<uint, Dynamic, 1>> hist(hgrid_->data(), hgrid_->size());
		Map<Matrix<double, Dynamic, 1>> uhist(ugrid_->data(), ugrid_->size());
		uhist.array() = pweight_*uhist.array() + hist.cast<double>()*(1./kbt_*bias_).array().exp()*weight_;
		ugrid_->syncGrid();
		hist.setZero();

		bias_.array() = kbt_*uhist.array().log();
		bias_.array() -= bias_.minCoeff();
		
		// Train network.
		net_.autoscale(hist_, bias_);
		if(world_.rank() == 0)
		{
			net_.train(hist_, bias_, false);
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

	void ANN::WriteBias()
	{
		net_.write("netstate.dat");
		
		std::string filename = overwrite_ ? outfile_ : outfile_ + std::to_string(sweep_);
		std::ofstream file(filename);
		file.precision(16);
		net_.forward_pass(hist_);
		matrix_t y = net_.get_activation();
		for(int i = 0; i < y.rows(); ++i)
		{
			for(int j = 0; j < hist_.cols(); ++j)
				file << std::fixed << hist_(i,j) << " ";
			file << std::fixed << ugrid_->data()[i] << " " << std::fixed << y(i) << "\n";
		}

		file.close();
	}

	ANN* ANN::Build(
		const Json::Value& json, 
		const MPI_Comm& world,
		const MPI_Comm& comm,
		const std::string& path)
	{
		ObjectRequirement validator;
		Value schema;
		Reader reader;
		
		reader.parse(JsonSchema::ANNMethod, schema);
		validator.Parse(schema, path);

		// Validate inputs.
		validator.Validate(json, path);
		if(validator.HasErrors())
			throw BuildException(validator.GetErrors());
		
		// Grid. 
		auto* fgrid = Grid<VectorXd>::BuildGrid(json.get("grid", Json::Value()));
		auto* hgrid = Grid<uint>::BuildGrid(json.get("grid", Json::Value()));
		auto* ugrid = Grid<double>::BuildGrid(json.get("grid", Json::Value()));

		// Topology. 
		auto nlayers = json["topology"].size() + 2;
		VectorXi topol(nlayers);
		topol[0] = fgrid->GetDimension();
		topol[nlayers-1] = 1;
		for(int i = 0; i < json["topology"].size(); ++i)
			topol[i+1] = json["topology"][i].asInt();
		
		auto weight = json.get("weight", 1.).asDouble();
		auto temp = json["temperature"].asDouble();
		auto nsweep = json["nsweep"].asUInt();

		// Assume all vectors are the same size. 
		std::vector<double> lowerb, upperb, lowerk, upperk;
		for(int i = 0; i < json["lower_bound_restraints"].size(); ++i)
		{
			lowerk.push_back(json["lower_bound_restraints"][i].asDouble());
			upperk.push_back(json["upper_bound_restraints"][i].asDouble());
			lowerb.push_back(json["lower_bounds"][i].asDouble());
			upperb.push_back(json["upper_bounds"][i].asDouble());
		}

		auto* m = new ANN(world, comm, topol, fgrid, hgrid, ugrid, lowerb, upperb, lowerk, upperk, temp, weight, nsweep);

		// Set optional params.
		m->SetPrevWeight(json.get("prev_weight", 1).asDouble());
		m->SetOutput(json.get("output_file", "ann.out").asString());
		m->SetOutputOverwrite( json.get("overwrite_output", true).asBool());
		m->SetConvergeIters(json.get("converge_iters", 0).asUInt());
		m->SetMaxIters(json.get("max_iters", 1000).asUInt());
		m->SetMinLoss(json.get("min_loss", 0).asDouble());

		return m;
	}
}
