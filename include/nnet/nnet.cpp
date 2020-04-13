#include <iostream>
#include <fstream>
#include <assert.h>
#include <cmath>
#include <random>
#include <Eigen/Dense>
#include "nnet.h"

using namespace Eigen;

namespace nnet
{
	neural_net::neural_net(const VectorXi& topology) 
	{
		assert(topology.size()>1);
		init_layers(topology);
		init_weights();
		autoscale_reset();
	}

	neural_net::neural_net(const char* filename) 
	{
		std::ifstream file(filename, std::ios::in);
		if(file) 
		{ 
			// number of layers
			int num_layers;
			file >> num_layers; 
			VectorXi topology(num_layers);
			
			// topology
			for(int i = 0; i < topology.size(); ++i)
				file >> topology[i];
			
			init_layers(topology);
			autoscale_reset();
			
			// scaling parameters
			for(int i = 0; i < x_scale_.size(); ++i) file >> x_scale_[i];
			for(int i = 0; i < x_shift_.size(); ++i) file >> x_shift_[i];
			for(int i = 0; i < y_scale_.size(); ++i) file >> y_scale_[i];
			for(int i = 0; i < y_shift_.size(); ++i) file >> y_shift_[i];
			
			// weights
			for(int i = 1; i < layers_.size(); ++i)
			{
				auto& layer = layers_[i];
				for(int j = 0; j < layer.b.size(); ++j) file >> layer.b(j);
				for(int j = 0; j < layer.W.size(); ++j) file >> layer.W(j);
			}
		}
		
		file.close();
    }

	void neural_net::init_layers(const VectorXi& topology)
	{
		// init input layer
		nparam_ = 0; 
		nn_layer l;
		l.size = topology(0);
		// added for gradient learning
		for (int i = 0; i < l.size; ++i)
		{		
			matrix_t dai;
			(l.da).push_back(dai);
			(l.delta2).push_back(dai);
			(l.delta3).push_back(dai);	
		}		
		layers_.push_back(l);

		// init hidden and output layer
		for(int i = 1; i < topology.size(); ++i)
		{
			nn_layer l;
			l.size	= topology(i);
			l.W.setZero(l.size, layers_[i-1].size);
			l.b.setZero(l.size);

			// added for gradient learning
			for (int j = 0; j < layers_.back().size; ++j)
			{		
				matrix_t daj;
				(l.da).push_back(daj);	
				(l.delta2).push_back(daj);
				(l.delta3).push_back(daj);	
			}
			layers_.push_back(l);
			nparam_ += l.W.size() + l.b.size();
		}
		tparams_ = {0.005, 1.e10, 10.0, 1.e-7, 0, 1.0, 1000};
    }

	void neural_net::init_weights()
	{
		std::random_device rd{};
		std::mt19937 gen{rd()};

		for(int i = 1; i < layers_.size(); ++i)
		{
			f_type beta = 0.7*std::pow(layers_[i].size, 1.0/layers_[i-1].size);
			f_type range = 2.0;

			std::uniform_real_distribution<> d{-1, 1};
			layers_[i].W = layers_[i].W.unaryExpr([&](f_type x){ return d(gen); });
			layers_[i].W.rowwise().normalize();
			layers_[i].W *= range*beta;

			if(layers_[i].size == 1)
				layers_[i].b.setZero();
			else
				layers_[i].b = range*beta*layers_[i].b.unaryExpr([&](f_type x){ return d(gen); });
		}
	}

	void neural_net::forward_pass(const matrix_t& X) 
	{
		assert(layers_.front().size == static_cast<int>(X.cols()));
		
		// copy and scale data matrix
		layers_[0].a.noalias() = (X.rowwise() - x_shift_.transpose())*x_scale_.asDiagonal();

		// input layer has const. 1st derivative
		for (int i = 0; i < layers_.front().size; ++i)
		{
			(layers_[0]).da[i].noalias() = matrix_t::Zero(X.rows(),X.cols());
			(layers_[0]).da[i].col(i) = vector_t::Ones(X.rows());
		} 
		
		for(int i = 1; i < layers_.size(); ++i)
		{
			// compute input for current layer
			layers_[i].z.noalias() = layers_[i-1].a * layers_[i].W.transpose();
			
			// add bias
			layers_[i].z.rowwise() += layers_[i].b.transpose(); 
			
			// apply activation function
			bool end = (i >= layers_.size() - 1);
			layers_[i].a = end ? layers_[i].z : activation(layers_[i].z);

			// forward propogate derivative -- added for gradient learning
			for (int j = 0; j < layers_.front().size; ++j)
			{
				if(end)
					(layers_[i]).da[j].noalias() = ((layers_[i-1]).da[j] * layers_[i].W.transpose());
				else
				{
					(layers_[i]).da[j].noalias() = (((((layers_[i-1]).da[j] * layers_[i].W.transpose()).array()) *\
						(activation_gradient(layers_[i].a)).array()).matrix());
				}
			}
		}
	}

	f_type neural_net::loss(const matrix_t& X, const matrix_t& Y)
	{
		assert(layers_.front().size == static_cast<int>(X.cols()));
		assert(layers_.back().size == static_cast<int>(Y.cols()));
		assert(X.rows() == Y.rows());
		
		// number of samples and output dim.
		int Q = Y.rows();
		int S = Y.cols();
		
		// Resize Jacobian and define error.
		je_.resize(nparam_);
		j_.resize(S*Q, nparam_);
		vector_t error(S*Q);
		
		// MSE. 
		f_type mse = 0.;
		
		for(int k = 0; k < Q; ++k)
		{
			// forward pass.
			forward_pass(X.row(k));
					
			// compute error.
			error.segment(k*S, S) = (layers_.back().a*y_scale_.asDiagonal().inverse()).transpose() -\
				(Y.row(k).transpose() - y_shift_);

			// compute loss.
			mse += error.segment(k*S, S).transpose().rowwise().squaredNorm().mean()/S;
			
			// Number of layers. 
			int m = layers_.size();

			// Compute sensitivities. 
			int j = nparam_;
			layers_[m-1].delta = y_scale_.asDiagonal().inverse();
			
			// Pack Jacobian.
			j -= layers_[m-1].W.size();
			for(int p = 0; p < S; ++p)
			{
				layers_[m-1].dEdW.noalias() = (layers_[m-1].delta.col(p)*layers_[m-2].a);
				j_.block(S*k+p, j, 1, layers_[m-1].W.size()) =\
				    Map<vector_t>(layers_[m-1].dEdW.data(), layers_[m-1].dEdW.size()).transpose();
			}
		
			j -= layers_[m-1].b.size();
			j_.block(S*k, j, S, layers_[m-1].b.size()) = layers_[m-1].delta;

			for(int i = layers_.size() - 2; i > 0; --i)
			{
				layers_[i].delta.noalias() = activation_gradient(layers_[i].a).asDiagonal()*\
					layers_[i+1].W.transpose()*layers_[i+1].delta;

				// Pack Jacobian.
				j -= layers_[i].W.size();
				for(int p = 0; p < S; ++p)
				{
					layers_[i].dEdW.noalias() = (layers_[i].delta.col(p)*layers_[i-1].a);						 
					j_.block(S*k+p, j, 1, layers_[i].W.size()) = \
   				        Map<vector_t>(layers_[i].dEdW.data(), layers_[i].dEdW.size()).transpose();
				}

				j -= layers_[i].b.size();
				j_.block(S*k, j, S, layers_[i].b.size()) = layers_[i].delta.transpose();
			}
		}

		jj_.noalias() = j_.transpose()*j_;
		jj_ /= (Q*S);
		j_ /= (Q*S);
		je_.noalias() = j_.transpose()*error;
		return mse/Q;
	}

	f_type neural_net::loss_w_grad(const matrix_t& X, const matrix_t& Y, const std::vector<matrix_t> &Z, const matrix_t& B)
	{
		assert(layers_.front().size == X.cols());
		assert(layers_.back().size == Y.cols());
		assert(X.rows() == Y.rows());

		// number of samples and output dim. 
		int Q = Y.rows();
		int S = Y.cols();
		int P = (X.cols()+1)*S;
		f_type Pd = (X.cols()*(1.0-tparams_.ratio)+1*tparams_.ratio)*S;
		
		// Resize jacobian and define error. 
		je_.resize(nparam_);
		j_.resize(P*Q, nparam_);
		vector_t error(P*Q);
		
		// MSE. 
		f_type mse = 0.;

		// Loop over samples.
		for(int k = 0; k < Q; ++k)
		{
			// forward pass
			forward_pass(X.row(k));

			// compute error
			error.segment(k*P, S) = ((layers_.back().a*y_scale_.asDiagonal().inverse()).transpose() -\
	            (Y.row(k).transpose() - y_shift_))*(1.0-B(k,0));

			// compute error on derivative
			for(int l = 0; l < X.cols(); ++l)
			{
				 error.segment(k*P+(l+1)*S, S) = (((layers_.back().da[l]*y_scale_.asDiagonal().inverse())*\
                     x_scale_.row(l) - Z[l].row(k)).transpose())*(1.0-B(k,0));
			}

			// Compute loss. 
			mse += error.segment(k*P, S).transpose().rowwise().squaredNorm().mean()/S;

			for(int l = 0; l < X.cols(); ++l)
				mse += error.segment(k*P+(l+1)*S, S).transpose().rowwise().squaredNorm().mean()/S;

			// Number of layers. 
			int m = layers_.size();

			// Compute sensitivities. 
			int j = nparam_;

			layers_[m-1].delta = y_scale_.asDiagonal().inverse();
			for(int l = 0; l < X.cols(); ++l)
			{
				//Set correct value to this for correct backpropagation from the linear layer to the first nonlinear layer.
				(layers_[m-1].delta2[l]) = ((y_scale_.asDiagonal().inverse()*x_scale_.row(l)*(1.0-B(k,0))*0.0));
				(layers_[m-1].delta3[l]) = ((y_scale_.asDiagonal().inverse()*x_scale_.row(l)*(1.0-B(k,0))*1.0));
			}

			// Pack Jacobian.
			j -= layers_[m-1].W.size();

			for(int p = 0; p < S; ++p)
			{
				layers_[m-1].dEdW = (layers_[m-1].delta.col(p)*layers_[m-2].a);
				j_.block(P*k+p, j, 1, layers_[m-1].W.size()) =\
					 Map<vector_t>(layers_[m-1].dEdW.data(), layers_[m-1].dEdW.size()).transpose();

				for(int l = 0; l < X.cols(); ++l)
				{
					layers_[m-1].dEdW = layers_[m-1].delta3[l].col(p)*layers_[m-2].da[l];
					j_.block(k*P+(l+1)*S+p, j, 1, layers_[m-1].W.size()) =\
                        Map<vector_t>(layers_[m-1].dEdW.data(), layers_[m-1].dEdW.size()).transpose();
				}
			}

			j -= layers_[m-1].b.size();
			j_.block(P*k, j, S, layers_[m-1].b.size()) = layers_[m-1].delta;

			for(int l = 0; l < X.cols(); ++l)
				j_.block(k*P+(l+1)*S, j, S, layers_[m-1].b.size()) = layers_[m-1].delta2[l].transpose()*0.0;

            for(int i = layers_.size() - 2; i > 0; --i)
            {
        		layers_[i].delta = activation_gradient(layers_[i].a).asDiagonal()*layers_[i+1].W.transpose()*layers_[i+1].delta;
	    	    for(int l = 0; l< X.cols(); ++l)
	        	{
		        	layers_[i].delta += \
			    		activation_gradient(layers_[i].a).asDiagonal()*layers_[i+1].W.transpose()*layers_[i+1].delta2[l];

		    	    (layers_[i].delta2)[l] = \
			    		((activation_secondgradient(layers_[i].a).asDiagonal()*layers_[i+1].W.transpose()*\
				    	layers_[i+1].delta3[l]).array() *(layers_[i].W * layers_[i-1].da[l].transpose()).array());

		    	    (layers_[i].delta3)[l] = \
			    		((activation_gradient(layers_[i].a).asDiagonal()*layers_[i+1].W.transpose()*layers_[i+1].delta3[l]));
        		}

	        	// Pack Jacobian.
	        	j -= layers_[i].W.size();
	        	for(int p = 0; p < S; ++p)
	        	{
		        	layers_[i].dEdW = (layers_[i].delta.col(p)*layers_[i-1].a);
		    	    j_.block(P*k+p, j, 1, layers_[i].W.size()) = Map<vector_t>(layers_[i].dEdW.data(), layers_[i].dEdW.size()).transpose();
			        for(int l = 0; l< X.cols(); ++l)
		        	{
				         layers_[i].dEdW = layers_[i].delta3[l].col(p)*layers_[i-1].da[l] + layers_[i].delta2[l].col(p)*layers_[i-1].a;
				         j_.block(k*P+(l+1)*S+p, j, 1, layers_[i].W.size()) = \
                            Map<vector_t>(layers_[i].dEdW.data(), layers_[i].dEdW.size()).transpose();
		        	}
	        	}

		        j -= layers_[i].b.size();
		        j_.block(P*k, j, S, layers_[i].b.size()) = layers_[i].delta.transpose();
	        	for(int l = 0; l< X.cols(); ++l)
			        j_.block(k*P+(l+1)*S, j, S, layers_[i].b.size()) = layers_[i].delta2[l].transpose();

			    }
		}

		jj_.noalias() = j_.transpose()*j_;

		jj_ /= S*Q;
		j_ /= S*Q;
		je_.noalias() = j_.transpose()*error;

		return mse/Q;

	}

	double neural_net::train(const matrix_t& X, const matrix_t& Y, bool verbose)
	{
		// Reset mu.
		tparams_.mu = 0.005;

		int nex = X.rows();
		
		// Forward and back propogate to compute loss and Jacobian.
		f_type mse = loss(X, Y);
		vector_t wb = get_wb(), optwb = get_wb();
		f_type wse = wb.transpose()*wb;

		// Initialize Bayesian regularization parameters.
		f_type gamma = nparam_;
		f_type beta = 0.5*(nex - gamma)/mse;
		beta = beta <= 0 ? 1. : beta;
		f_type alpha = 0.5*gamma/wse;
		f_type tse = beta*mse + alpha*wse;
		f_type grad = 2.*je_.norm();
		matrix_t eye = matrix_t::Identity(nparam_, nparam_);
		
		// Define iteration tol parameters.
		int iter = 0;
		f_type tse2 = 0, mse2 = 0, wse2 = 0;
		do
		{
			if(verbose)
				std::cout << "iter: " << iter << " mse: " << mse << " gamma: " << gamma << " mu: " << tparams_.mu << " grad: " << grad << std::endl;
			
			matrix_t jjb = jj_; 
			vector_t jeb = je_;

			do
			{
				// Compute new weights and performance.
				optwb = wb - (beta*jjb + (tparams_.mu + alpha)*eye).colPivHouseholderQr().solve(beta*jeb + alpha*wb);
				wse2 = optwb.transpose()*optwb;
				set_wb(optwb);
				mse2 = loss(X, Y);
				tse2 = beta*mse2 + alpha*wse2;
				
				// Exit loop or reset values.
				if(tse2 < tse || tparams_.mu > tparams_.mu_max)
					break;
				else
				{
					set_wb(wb);
					//mse2 = loss(X, Y);
					tparams_.mu *= tparams_.mu_scale;
				}
			} while(true);

			wb = optwb;
			mse = mse2; wse = wse2;
			gamma = (f_type)nparam_ - alpha*(beta*jj_ + alpha*eye).inverse().trace();
			beta = mse == 0 ? 1. : 0.5*((f_type)nex - gamma)/mse;
			alpha = wse == 0 ? 1. : 0.5*gamma/wse;
			tse = beta*mse + alpha*wse;
            f_type grad = 2.*je_.norm();
			
			if(tparams_.mu < tparams_.mu_max)
				tparams_.mu /= tparams_.mu_scale;
			if(tparams_.mu < 1.e-20) tparams_.mu = 1.e-20;

			++iter;				 
		} while(
			tparams_.mu < tparams_.mu_max && 
			grad > tparams_.min_grad && 
			iter <= tparams_.max_iter && 
			mse > tparams_.min_loss &&
			!std::isnan(grad) && 
			!std::isnan(gamma));

		if(verbose)
			std::cout << "iter: " << iter << " mse: " << mse << " gamma: " << gamma << " mu: " << tparams_.mu << " grad: " << grad << std::endl;

		return gamma;
	}

	double neural_net::train_w_grad(const matrix_t& X, const matrix_t& Y, const std::vector<matrix_t> &Z,\
		const matrix_t& B, bool verbose)
	{
        // Reset mu.
		tparams_.mu = 0.005;
		int nex = X.rows()*(X.cols()+1);		
        
		// Forward and back propogate to compute loss and Jacobian.
		f_type mse = loss_w_grad(X, Y, Z, B);
		vector_t wb = get_wb(), optwb = get_wb();
		f_type wse = wb.transpose()*wb;

		// Initialize Bayesian regularization parameters.
		f_type gamma = nparam_;
		f_type beta = 0.5*(nex - gamma)/mse;
		beta = beta <= 0 ? 1. : beta;
		f_type alpha = 0.5*gamma/wse;
		f_type tse = beta*mse + alpha*wse;
		f_type grad = 2.*je_.norm();
		matrix_t eye = matrix_t::Identity(nparam_, nparam_);
		
		// Define iteration tol parameters.
		int iter = 0;
		f_type tse2 = 0, mse2 = 0, wse2 = 0;
		do
		{
			if(verbose)
				std::cout << "iter: " << iter << " mse: " << mse << " gamma: " << gamma << " mu: " << tparams_.mu << " grad: " << grad << std::endl;
			
			matrix_t jjb = jj_; 
			vector_t jeb = je_;

			do
			{
				// Compute new weights and performance.
				optwb = wb - (beta*jjb + (tparams_.mu + alpha)*eye).colPivHouseholderQr().solve(beta*jeb + alpha*wb);
				wse2 = optwb.transpose()*optwb;
				set_wb(optwb);
				mse2 = loss_w_grad(X, Y, Z, B);
				tse2 = beta*mse2 + alpha*wse2;
				
				// Exit loop or reset values.
				if(tse2 < tse || tparams_.mu > tparams_.mu_max)
					break;
				else
				{
					set_wb(wb);
					//mse2 = loss(X, Y);
					tparams_.mu *= tparams_.mu_scale;						
				}
			} while(true);

			wb = optwb;
			mse = mse2; wse = wse2;
			gamma = (f_type)nparam_ - alpha*(beta*jj_ + alpha*eye).inverse().trace();
			beta = mse == 0 ? 1. : 0.5*((f_type)nex - gamma)/mse;
			alpha = wse == 0 ? 1. : 0.5*gamma/wse;
			tse = beta*mse + alpha*wse;
			grad = 2.*je_.norm();
			
			if(tparams_.mu < tparams_.mu_max)
				tparams_.mu /= tparams_.mu_scale;
			if(tparams_.mu < 1.e-20) tparams_.mu = 1.e-20;

			++iter;				 
		} while(
			tparams_.mu < tparams_.mu_max && 
			grad > tparams_.min_grad && 
			iter <= tparams_.max_iter && 
			mse > tparams_.min_loss &&
			!std::isnan(grad) && 
			!std::isnan(gamma));

		if(verbose)
			std::cout << "iter: " << iter << " mse: " << mse << " gamma: " << gamma << " mu: " << tparams_.mu << " grad: " << grad << std::endl;

	    return gamma;
	
	}


	train_param neural_net::get_train_params() const
	{
		return tparams_;
	}

	void neural_net::set_train_params(const train_param& params)
	{
		tparams_ = params;
	}

	matrix_t neural_net::get_activation() 
	{
		return (layers_.back().a*y_scale_.asDiagonal().inverse()).rowwise() + y_shift_.transpose();
	}

	matrix_t neural_net::get_gradient(int index)
	{
		layers_.back().delta = matrix_t::Identity(layers_.back().size, layers_.back().size)*y_scale_.asDiagonal().inverse();
		for(int i = layers_.size() - 2; i > 0; --i)
		{
			matrix_t g = activation_gradient(layers_[i].a);
			layers_[i].delta = (layers_[i+1].delta*layers_[i+1].W)*g.row(index).asDiagonal();
		}

		return layers_[1].delta*layers_[1].W*x_scale_.asDiagonal();
	}

	// this is tanh(x)
	matrix_t neural_net::activation(const matrix_t& x)	
	{
		return (2.*((-2.*x).array().exp() + 1.0).inverse() - 1.0).matrix();
	}

	// derivative is 1-tanh^2(x) = sech^2
	matrix_t neural_net::activation_gradient(const matrix_t& x) 
	{
		return (1.0-x.array().square()).matrix();
	}

	matrix_t neural_net::get_gradient_forwardpass(int index) 
	{
		return ((layers_.back().da[index])*y_scale_.asDiagonal().inverse())*x_scale_.row(index);
	}

	// 2nd derivative is -2*tanh*sech^2 = -2*x*(1-x^2) 
	matrix_t neural_net::activation_secondgradient(const matrix_t& x) 
	{
		return (-2.*(1.0-x.array().square())*(x.array())).matrix();
	}

	void neural_net::set_wb(const vector_t& wb)
	{
		int k = 0;
		for(int i = 1; i < layers_.size(); ++i)
		{
			auto& layer = layers_[i];
			for(int j = 0; j < layer.b.size(); ++j)
			{
				layer.b(j) = wb[k];
				++k;
			}

			for(int j = 0; j < layer.W.size(); ++j)
			{
				layer.W(j) = wb[k];
				++k;
			}
		}
	}

	vector_t neural_net::get_wb() const
	{
		int k = 0;
		vector_t wb(nparam_);
		for(int i = 1; i < layers_.size(); ++i)
		{
			auto& layer = layers_[i];
			for(int j = 0; j < layer.b.size(); ++j)
			{
				wb[k] = layer.b(j);
				++k;
			}

			for(int j = 0; j < layer.W.size(); ++j)
			{
				wb[k] = layer.W(j);
				++k;
			}
		}

		return wb;
	}

	void neural_net::autoscale(const matrix_t& X, const matrix_t& Y) 
	{
		assert(layers_.front().size == static_cast<int>(X.cols()));
		assert(layers_.back().size == static_cast<int>(Y.cols()));
		assert(X.rows() == Y.rows());
		
		x_shift_ = 0.5*(X.colwise().minCoeff().array() + X.colwise().maxCoeff().array());
		x_scale_ = 2.0*(X.colwise().maxCoeff() - X.colwise().minCoeff()).array().inverse();

		y_shift_ = 0.5*(Y.colwise().minCoeff().array() + Y.colwise().maxCoeff().array());
		y_scale_ = 2.0*(Y.colwise().maxCoeff() - Y.colwise().minCoeff()).array().inverse();
	}

	void neural_net::autoscale_w_grad(const matrix_t& X, const matrix_t& Y, const std::vector<matrix_t> & Z)
	{
		assert(layers_.front().size == static_cast<int>(X.cols()));
		assert(layers_.back().size == static_cast<int>(Y.cols()));
		assert(X.rows() == Y.rows());
		
		x_shift_ = 0.5*(X.colwise().minCoeff().array() + X.colwise().maxCoeff().array());
		x_scale_ = 2.0*(X.colwise().maxCoeff() - X.colwise().minCoeff()).array().inverse();

		y_shift_ = 0.5*(Y.colwise().minCoeff().array() + Y.colwise().maxCoeff().array());
		y_scale_ = 2.0*(Y.colwise().maxCoeff() - Y.colwise().minCoeff()).array().inverse();

		for (int i = 0; i < Z.size(); ++i)
		{
			vector_t temp = 2.0*(Z[i].colwise().maxCoeff() - Z[i].colwise().minCoeff()).array().inverse();
			y_scale_ =	temp.array().min(y_scale_.array());
		}
	}
	
	void neural_net::autoscale_reset() 
	{
		x_scale_ = vector_t::Ones(layers_.front().size);
		x_shift_ = vector_t::Zero(layers_.front().size);
		y_scale_ = vector_t::Ones(layers_.back().size);
		y_shift_ = vector_t::Zero(layers_.back().size);
	}

	bool neural_net::write(const char* filename) 
	{
		// open file
		std::ofstream file(filename, std::ios::out);
		
		// write everything to disk
		if(file) 
		{
			file.precision(16);
			file << std::scientific;

			// number of layers
			file << static_cast<int>(layers_.size()) << std::endl;

			// topology
			for(int i = 0; i < layers_.size() - 1; ++i)
				file << static_cast<int>(layers_[i].size) << " ";
			file << static_cast<int>(layers_.back().size) << std::endl;

			for(int i = 0; i < x_scale_.size() - 1; ++i) file << x_scale_[i] << " ";
			file << x_scale_[x_scale_.size()-1] << std::endl;

			for(int i = 0; i < x_shift_.size() - 1; ++i) file << x_shift_[i] << " ";
			file << x_shift_[x_shift_.size()-1] << std::endl;

			for(int i = 0; i < y_scale_.size() - 1; ++i) file << y_scale_[i] << " ";
			file << y_scale_[y_scale_.size()-1] << std::endl;

			for(int i = 0; i < y_shift_.size() - 1; ++i) file << y_shift_[i] << " ";
			file << y_shift_[y_shift_.size()-1] << std::endl;

			// weights
			for(int i = 1; i < layers_.size(); ++i)
			{
				auto& layer = layers_[i];
				for(int j = 0; j < layer.b.size(); ++j) file << layer.b(j) << std::endl;
				for(int j = 0; j < layer.W.size(); ++j) file << layer.W(j) << std::endl;
			}
		} 
		else
			return false;
		
		file.close();
		return true;
	}

	neural_net::~neural_net() 
	{
	}
}
