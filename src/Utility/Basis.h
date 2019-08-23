#include <cmath>
#include "json/json.h"
#include "../Grids/Grid.h"

#define UNUSED(x) (void)(x)

namespace SSAGES
{
    //! Look-up table for basis functions.
	/*!
     * The structure that holds the Look-up table for the basis function. To
     * prevent repeated calculations, both the derivatives and values of the
     * Legendre polynomials is stored here. More will be added in future
     * versions.
     */
	struct BasisLUT
	{
		//! The values of the basis sets
		std::vector<double> values;

		//! The values of the derivatives of the basis sets
		std::vector<double> derivs;

        //! Constructor.
        /*!
         * \param values The values of the basis sets.
         * \param derivs The values of the derivatives of the basis sets.
         */
		BasisLUT(const std::vector<double>& values,
			const std::vector<double>& derivs) :
			values(values), derivs(derivs)
		{}
	};

    //! Map for histogram and coefficients.
    /*!
     * A clean mapping structure for both the histogram and the coefficients.
     * All vectors are written as 1D with a row major mapping. In order to make
     * iterating easier, the mapping of the 1D vectors are written here.
     */
    struct Map
    {
        //! The coefficient value
        double value;

        //! The mapping in an array of integers
        std::vector<int> map;

        //! Constructor
        /*!
         * \param map The mapping in an array of integers.
         * \param value The coefficient value.
         */
        Map(const std::vector<int>& map,
            double value) :
            value(value), map(map)
        {}
    };

    //! Abstract class for all BasisFunction inheritance.
    class BasisFunction
    {
    protected:
        unsigned int polyOrd_; //!< Order of the polynomial.
        unsigned int nbins_; //!< Number of bins.
        bool isFinite_; //!< Flag for finite-range polynomials.
        bool zeroOrder_; //!< Flag for constant-order polynomials.
        double boundLow_; //!< Lower bound on CV.
        double boundUp_; //!< Upper bound on CV.

    public:
        //! Constructor
        /*!
         * \param polyOrd Order of polynomial.
         * \param nbins Number of bins.
         * \param isFinite Flag for finite-range polynomials.
         * \param zeroOrder Flag for constant-order polynomials.
         * \param boundLow Lower bounds of restraint springs.
         * \param boundUp Upper bounds of restraint springs.
         *
         * Constructs an instance of the Basis function class, which can
         * have multiple different inherited forms.
         */ 
        BasisFunction(unsigned int polyOrd,
                      unsigned int nbins,
                      bool isFinite,
                      bool zeroOrder,
                      double boundLow,
                      double boundUp) :
        polyOrd_(polyOrd), nbins_(nbins), isFinite_(isFinite), zeroOrder_(zeroOrder),
        boundLow_(boundLow), boundUp_(boundUp)
        {
        }

        //! Gets the order of the current polynomial.
        /*!
         * \return Order of current polynomial.
         */
        unsigned int GetOrder() {return polyOrd_;}

        //! Gets the number of bins for the discretization.
        /*!
         * \return Number of bins.
         */
        unsigned int GetBins() {return nbins_;}

        //! Gets the flag for constant-order polynomials.
        /*!
         * \return Flag for constant-order polynomials.
         */
        bool GetZeroOrder() {return zeroOrder_;}

        //! Gets the lower bound on the CV.
        /*!
         * \return Lower bound on CV.
         */
        double GetLower() {return boundLow_;}

        //! Gets the upper bound on the CV.
        /*!
         * \return Upper bound on CV.
         */
        double GetUpper() {return boundUp_;}

        //! Gets the magnitude of the range of bounds.
        /*!
         * \return Magnitude of range of bounds on CV.
         */
        double GetRange() 
        {
            if(isFinite_) 
                return boundUp_ - boundLow_;
            // No infinitely bounded basis functions are included currently so this is going to return nothing for right now
            else
                return 0.0;
        }

        //! Gets the norm of the basis function.
        /*!
         * \param order Order of the polynomial.
         *
         * \return Norm of a specific order of polynomial.
         */
        virtual double GetNorm(int order) {
            UNUSED(order);
            return 0;
        }

        //! Calculates the output of the basis function.
        /*!
         * \param val Input value for function.
         * \param order Order of the polynomial.
         *
         * \return Output value for function.
         */
        virtual double Evaluate(double val, int order) {
            UNUSED(val);
            UNUSED(order);
            return 0;
        }

        //! Calculates the gradient of the basis function.
        /*!
         * \param val Input value for function.
         * \param order Order of the polynomial.
         *
         * \return Gradient for function.
         */
        virtual double EvalGrad(double val, int order) {
            UNUSED(val);
            UNUSED(order);
            return 0;
        }

        //! Calculates the gradient of the basis function.
        /*!
         * \param val Input value for function.
         *
         * \return Gradient of function.
         */
        virtual double Weight(double val) {
            UNUSED(val);
            return 1;
        }

        //! Build BasisFunction from JSON value. 
        /*!
         * \param json JSON value node.
         * \param path Path for JSON path specification.
         * \param nbins Number of bins.
         * 
         * \return Pointer to new BasisFunction.
         */
        static BasisFunction* Build(const Json::Value& json, const std::string& path, unsigned int nbins);

        //! Destructor
        virtual ~BasisFunction() {}
    };

    //! Defines the class of Chebyshev polynomials
    class Chebyshev : public BasisFunction
    {
    private:  
    protected:
    public:
        //! Constructor
        /*!
         * \param polyOrd Order of Chebyshev polynomial.
         * \param boundLow Lower bounds of restraint springs.
         * \param boundUp Upper bounds of restraint springs.
         * \param nbins Number of bins.
         *
         * Constructs an instance of the Chebyshev function class.
         */ 
        Chebyshev(unsigned int polyOrd, double boundLow, double boundUp, unsigned int nbins) :
        BasisFunction(polyOrd, nbins, true, false, boundLow, boundUp)
        {
        }

        //! Transforms variable from absolute value to [-1,1] between boundaries.
        /*!
         * \param x Variable to transform.
         *
         * \return Relative value between lower (-1) and upper (1) bounds on CV.
         */
        double ChangeVariable(double x)
        {
            return (x-boundLow_)*2.0/(boundUp_ - boundLow_)-1.0;
        }

        double Evaluate(double x, int n)
        {
            double xMod = ChangeVariable(x);
            return n == 0 ? 1.0 : 
                    n == 1 ? xMod :
                    2.0*xMod*Evaluate(x,n-1) - Evaluate(x,n-2);
        }

        double EvalGrad(double x, int n)
        {
            double xMod = ChangeVariable(x);
            return n == 0 ? 0.0 :
                    n == 1 ? 1.0 :
                    2.0*(Evaluate(x,n-1) + xMod * EvalGrad(x,n-1)) - EvalGrad(x,n-2);
        }

        double Weight(double x)
        {
            double xMod = ChangeVariable(x);
            return 1.0/sqrt(1.0-xMod*xMod);
        }

        double GetNorm(int /* n */)
        {
            return 2.0/M_PI;
        }
        
        //! Build the Chebyshev polynomial
        //! \copydoc BasisFunction::Build()
        static Chebyshev* Build(const Json::Value& json, const std::string& path, unsigned int nbins);

    };

    //! Defines the class of Legendre polynomials
    class Legendre : public BasisFunction
    {
    private:  
    protected:
    public:
        //! Constructor
        /*!
         * \param polyOrd Order of Legendre polynomial.
         * \param nbins Number of bins.
         *
         * Constructs an instance of the Legendre function class.
         */ 
        Legendre(unsigned int polyOrd, unsigned int nbins) :
        BasisFunction(polyOrd, nbins, true, false, -1.0, 1.0)
        {
        }

        double Evaluate(double x, int n)
        {
            return n == 0 ? 1.0 : 
                    n == 1 ? x :
                    (2.0*n-1.0)/(double)n*x*Evaluate(x,n-1) - (n-1.0)/(double)n*Evaluate(x,n-2);
        }

        double EvalGrad(double x, int n)
        {
            return n == 0 ? 0.0 :
                    n == 1 ? 1.0 :
                    (2.0*n-1.0)/(double)n*(Evaluate(x,n-1) + x * EvalGrad(x,n-1)) - (n-1.0)/(double)n*EvalGrad(x,n-2);
        }

        double GetNorm(int n)
        {
            return n + 0.5;
        }

        //! Build the Legendre polynomial
        //! \copydoc BasisFunction::Build()
        static Legendre* Build(const Json::Value& json, const std::string& path, unsigned int nbins);
    };

    //! Defines the class of Fourier polynomials
    class Fourier : public BasisFunction
    {
    private:
    protected:
    public:
        //! Constructor
        /*!
         * \param polyOrd Order of Fourier polynomial.
         * \param boundLow Lower bounds of restraint springs.
         * \param boundUp Upper bounds of restraint springs.
         * \param nbins Number of bins.
         *
         * Constructs an instance of the Fourier function class.
         */ 
        Fourier(unsigned int polyOrd, double boundLow, double boundUp, unsigned int nbins) :
        BasisFunction(polyOrd, nbins, true, true, boundLow, boundUp)
        {
        }

        //For the Fourier series, I am choosing to split the coefficients in two
        //The first odd coefficients are the sine series and the even are cosine
        //However when evaluated, they will be treated as n, just stored differently
        double Evaluate(double x, int n)
        {
            return n == 0 ? 0.5 : 
                    n % 2 == 1 ? std::sin(2.0*floor(double(n)*0.5)*M_PI*x/(boundUp_-boundLow_)) :
                      std::cos(2.0*floor(double(n)*0.5)*M_PI*x/(boundUp_-boundLow_));
        }

        double EvalGrad(double x, int n)
        {
            return n == 0 ? 0.0 :
                    n % 2 == 0 ? (-2.0*floor(double(n)*0.5)*M_PI/(boundUp_-boundLow_))*(std::sin(2.0*floor(double(n)*0.5)*M_PI*x/(boundUp_-boundLow_))) : 
                     (2.0*floor(double(n)*0.5)*M_PI/(boundUp_-boundLow_))*std::cos(2.0*floor(double(n)*0.5)*M_PI*x/(boundUp_-boundLow_));
        }

        double GetNorm(int /* n */)
        {
            //Assumes the period of the function is the range of the CV
            return 2.0/(boundUp_-boundLow_);
        }

        //! Build the Fourier polynomial
        //! \copydoc BasisFunction::Build()
        static Fourier* Build(const Json::Value& json, const std::string& path, unsigned int nbins);
    };
    
    //! Calculates the inner product of all the basis functions and the histogram
    class BasisEvaluator
    {
    private:
        std::vector<Map> coeff_; //!< Vector of basis function coefficients.
        std::vector<BasisFunction*> functions_; //!< Vector of basis functions.
        std::vector<BasisLUT> lookup_; //!< Lookup table for basis functions.

    public:
        //! Constructor
        /*!
         * \param functions Vector of BasisFunction members.
         *
         * Initialize the evaluator.
         */
        BasisEvaluator(const std::vector<BasisFunction*>& functions) : 
            functions_(functions)
        {
            CoeffInit();
            BasisInit();
        }
    
        //! For now, when the basis evaluator is called it will store all the values
        //! of the basis functions into a lookup table. [Subject to change]
        void CoeffInit(); //!< Initializes coefficient lookup vector
        void BasisInit(); //!< Creates lookup table for basis functions

        //! Outputs the basis projection at a specific coordinate
        /*!
         * \param bias Grid of current bias
         * \param grad Grid of gradient to be applied
         */
        void UpdateBias(Grid<double> *bias, Grid<std::vector<double>> *grad);

        //! Calculates the inner product of the basis set and the biased histogram
        //! This function then returns the coefficients from this calculation
        /*!
         * \param array Basis set.
         * \param hist Biased histrogram.
         *
         * \return Inner product of the basis set and the biased histogram.
         */
        double UpdateCoeff(const std::vector<double> &array, Grid<unsigned int> *hist);

        //! Gets the coefficient array
        /*!
         * \return Vector of current coefficients
         */
        std::vector<double> GetCoeff(void)
        {
            std::vector<double> coeffArray (coeff_.size(),0);
            for (size_t i=0; i<coeff_.size(); i++)
            {
                coeffArray[i] = coeff_[i].value;
            }
            return coeffArray; 
        }

        //! Set the coefficient array in the event of a restart run
        /*!
         * \param coeff Vector of coefficients to set
         */
        void SetCoeff(const std::vector<double>& coeff)
        {
            size_t ii = 0;
            for(auto& val : coeff_) 
            {    
                val.value = coeff[ii];
                ii++;
            }
        }

        ~BasisEvaluator() {}
    };
} 
