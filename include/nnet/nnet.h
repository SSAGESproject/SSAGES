#ifndef __NN_H__
#define __NN_H__

#define f_type double

#include <Eigen/Core>
#include <vector>

namespace nnet
{
	typedef Eigen::Matrix<f_type, Eigen::Dynamic, Eigen::Dynamic> matrix_t;
	typedef Eigen::Matrix<f_type, Eigen::Dynamic, 1> vector_t;
	typedef Eigen::Array<f_type, Eigen::Dynamic, Eigen::Dynamic> array_t;

	struct nn_layer
	{
		int size;
		matrix_t a, z, delta;
		std::vector<matrix_t> delta2, delta3;
		matrix_t W, dEdW;
		vector_t b;

		// add vector of matrix to hold d(node)/d(input)
		std::vector<matrix_t> da;
	};

	struct train_param
	{
		f_type mu, mu_max, mu_scale, min_grad, min_loss, ratio;
		int max_iter;
	};

	class neural_net
	{
	private:
		/** Allocate memory and initialize default values. */
		void init_layers(const Eigen::VectorXi &topology);

		/** Holds the layers of the neural net. */
		std::vector<nn_layer> layers_;

		/** Training params. */
		train_param tparams_;

		/** Holds the error gradient, jacobian, ... */
		matrix_t j_, jj_;
		vector_t je_;
		std::vector<matrix_t> jder_;

		/** Number of adjustable parameters. */
		unsigned int nparam_;

		/** Scaling parameters. */
		vector_t x_shift_;
		vector_t x_scale_;
		vector_t y_shift_;
		vector_t y_scale_;

	public:
		/** Init neural net with given topology. */
		neural_net(const Eigen::VectorXi& topology);

		/** Read neural net from file. */
		neural_net(const char* filename);

		/** Initial weights randomly (modified Nguyen-Widrow) . */
		void init_weights();

		/** Propagate data through the net.
		*  Rows of X are instances, columns are features. */
		void forward_pass(const matrix_t& X);

		/** Compute NN loss w.r.t. input and output data.
		* Also backpropogates error.
		*/
		// Loss function for frequency-based learning.
		f_type loss(const matrix_t& X, const matrix_t& Y);
		// Loss function for gradient-based learning.
		f_type loss_w_grad(
			const matrix_t& X, const matrix_t& Y,
			const std::vector<matrix_t> &Z, const matrix_t& B
		);

		/** Train neural network. */
		// Training function for frequency-based learning
		double train(const matrix_t& X, const matrix_t& Y, bool verbose = false);
		// Training function for gradient-based learning
		double train_w_grad(
			const matrix_t& X, const matrix_t& Y,
			const std::vector<matrix_t> &Z, const matrix_t& B, bool verbose = false
		);

		/** Get training parameters. */
		train_param get_train_params() const;

		/** Set training parameters. */
		void set_train_params(const train_param& params);

		/** Return activation of output layer. */
		matrix_t get_activation();

		/** Get gradient of output(s) w.r.t. input i, calculated via forward pass */
		matrix_t get_gradient_forwardpass(int index);

		/** Get gradient of output(s) w.r.t. input i */
		vector_t get_gradient(int index);

		/** Returns the logistic function values f(x) given x. */
		static matrix_t activation(const matrix_t& x);

		/** Returns the gradient f'(x) of the logistic function given f(x). */
		static matrix_t activation_gradient(const matrix_t& x);

		/** Returns the second derivative f''(x) of the logistic function given f(x). */
		static matrix_t activation_secondgradient(const matrix_t& x);

		/** Set weights and biases. */
		void set_wb(const vector_t& wb);

		/** Get weights and biases. */
		vector_t get_wb() const;

		/** Compute autoscale parameters. */
		// For frequency-based learning
		void autoscale(const matrix_t& X, const matrix_t& Y);
		// For gradient-based learning
		void autoscale_w_grad(const matrix_t& X, const matrix_t& Y, const std::vector<matrix_t> &Z);

		/** Reset autoscale parameters */
		void autoscale_reset();

		/** Write net parameter to file. */
		bool write(const char* filename);

		/** Destructor. */
		~neural_net();
	};
}
#endif
