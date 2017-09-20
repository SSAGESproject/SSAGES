#ifndef __NN_H__
#define __NN_H__

#define f_type float

#include <Eigen/Core>
#include <vector>

namespace nnet
{
	typedef Eigen::Matrix<f_type, Eigen::Dynamic, Eigen::Dynamic> matrix_t;
	typedef Eigen::Matrix<f_type, Eigen::Dynamic, 1> vector_t;
	typedef Eigen::Array<f_type, Eigen::Dynamic, Eigen::Dynamic> array_t;

	struct nn_layer 
	{
		size_t size;
		matrix_t a, z, delta;
		matrix_t W, dEdW;
		vector_t b;
	};

	struct train_param
	{
		f_type mu, mu_max, mu_scale, min_grad;
		int max_iter;
	};

	class neural_net 
	{
	private:
		/** Allocate memory and initialize default values. */
		void init_layers(Eigen::VectorXi &topology);
		
		/** Holds the layers of the neural net. */
		std::vector<nn_layer> layers_;

		/** Training params. */
		train_param tparams_;

		/** Holds the error gradient, jacobian, ... */ 
		matrix_t j_, jj_;
		vector_t je_;

		/** Number of adjustable parameters. */
		uint nparam_;

		/** Scaling parameters. */
		vector_t x_shift_;
		vector_t x_scale_;
		vector_t y_shift_;
		vector_t y_scale_;
		
	public:
		/** Init neural net with given topology. */
		neural_net(Eigen::VectorXi& topology);
		
		/** Read neural net from file. */
		neural_net(const char* filename);

		/** Initial weights randomly (zero mean, standard deviation sd) . */
		void init_weights(f_type sd);
		
		/** Propagate data through the net.
		*  Rows of X are instances, columns are features. */
		void forward_pass(const matrix_t& X);

		/** Compute NN loss w.r.t. input and output data.
		* Also backpropogates error. 
		*/
		f_type loss(const matrix_t& X, const matrix_t& Y);

		void train(const matrix_t& X, const matrix_t& Y, bool verbose = false);

		/** Return activation of output layer. */
		matrix_t get_activation();
		
		/** Get gradient of output(s) w.r.t. input i */
		matrix_t get_gradient(int index);
		
		/** Returns the logistic function values f(x) given x. */
		static matrix_t activation(const matrix_t& x);
		
		/** Returns the gradient f'(x) of the logistic function given f(x). */
		static matrix_t activation_gradient(const matrix_t& x);
		
		/** Set weights and biases. */ 
		void set_wb(const vector_t& wb);

		/** Get weights and biases. */
		vector_t get_wb() const; 
		
		/** Compute autoscale parameters. */
		void autoscale(const matrix_t& X, const matrix_t& Y);

		/** Reset autoscale parameters */
		void autoscale_reset();      
		
		/** Write net parameter to file. */
		bool write(const char* filename);
		
		/** Destructor. */
		~neural_net();
	};
}
#endif
