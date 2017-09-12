#include <cmath>
#include "json/json.h"
#include "../Grids/Histogram.h"

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

    class BasisFunction
    {
    protected:
        uint polyOrd_; 
        uint nbins_;
        double boundLow_, boundHigh_;
        bool isFinite_;

    public: 
        BasisFunction(uint polyOrd,
                      uint nbins,
                      bool isFinite,
                      double boundLow,
                      double boundHigh) :
        polyOrd_(polyOrd), isFinite_(isFinite), nbins_(nbins),
        boundLow_(boundLow), boundHigh_(boundHigh)
        {
        }

        uint GetOrder() {return polyOrd_;}
        uint GetBins() {return nbins_;}
        double GetLower() {return boundLow_;}
        double GetUpper() {return boundHigh_;}
        double GetRange() 
        {
            if(isFinite_) 
                return boundHigh_ - boundLow_;
            // No infinitely bounded basis functions are included currently so this is going to return nothing for right now
            else
                return 0.0;
        }

        virtual double Evaluate(double val, int order) {return 0;}
        virtual double EvalGrad(double grad, int order) {return 0;}
        //virtual std::vector<double> ConvBounds(std::vector<double> x, double max, double min) = 0;
		static BasisFunction* Build(const Json::Value& json, const std::string& path, uint nbins);

        //! Destructor
        virtual ~BasisFunction()
        {
        }
    };

    class Chebyshev : public BasisFunction
    {
    private:  
    protected:
    public:
        Chebyshev(unsigned int polyOrd, double boundLow, double boundHigh, unsigned int nbins) :
        BasisFunction(polyOrd, nbins, true, boundLow, boundHigh)
        {
        }

        //Quick recursive relation to evaluate the basis sets
        double Evaluate(double x, int n)
        {
            return n == 0 ? 1.0 : 
                    n == 1 ? x :
                    2.0*x*Evaluate(x,n-1) - Evaluate(x,n-2);
        }
        //Same but for the gradients
        double EvalGrad(double x, int n)
        {
            return n == 0 ? 0.0 :
                    n == 1 ? 1.0 :
                    2.0*(Evaluate(x,n-1) + x * EvalGrad(x,n-1)) - EvalGrad(x,n-2);
        }
        
        //Build the Chebyshev polynomial
        static Chebyshev* Build(const Json::Value& json, const std::string& path, uint nbins);

    };

    class Legendre : public BasisFunction
    {
    private:  
    protected:
    public:
        Legendre(unsigned int polyOrd, unsigned int nbins) :
        BasisFunction(polyOrd, nbins, true, -1.0, 1.0)
        {
        }

        //Quick recursive relation to evaluate the basis sets
        virtual double Evaluate(double x, int n)
        {
            return n == 0 ? 1.0 : 
                    n == 1 ? x :
                    (2.0*n-1.0)/(double)n*x*Evaluate(x,n-1) - (n-1.0)/(double)n*Evaluate(x,n-2);
        }
        //Same but for the gradients
        virtual double EvalGrad(double x, int n)
        {
            return n == 0 ? 0.0 :
                    n == 1 ? 1.0 :
                    (2.0*n-1.0)/(double)n*(Evaluate(x,n-1) + x * EvalGrad(x,n-1)) - (n-1.0)/(double)n*EvalGrad(x,n-2);
        }

        //Build the Legendre polynomial
        static Legendre* Build(const Json::Value& json, const std::string& path, uint nbins);
    };

    //Calculates the inner product of all the basis functions and the histogram
    class BasisEvaluator
    {
    private:
        std::vector<Map> coeff_;
        std::vector<BasisFunction*> functions_;
        std::vector<BasisLUT> lookup_;

    public:
        //Initialize the evaluator
        BasisEvaluator(const std::vector<BasisFunction*>& functions) : 
            functions_(functions)
        {
            CoeffInit();
            BasisInit();
        }
    
        //For now when the basis evaluator is called it will store all the values
        //of the basis functions into a lookup table. Subject to change
        void CoeffInit(void);
        void BasisInit(void);

        //Outputs the basis projection at a specific coordinate
		void UpdateBias(Histogram<double> *bias, Histogram<std::vector<double>> *grad);

        //Calculates the inner product of the basis set and the biased histogram
        //This function then returns the coefficients from this calculation
		double UpdateCoeff(const std::vector<double> &array, Histogram<uint> *hist);

        //Gets the coefficient array
        std::vector<double> GetCoeff(void)
        {
            std::vector<double> coeffArray (coeff_.size(),0);
            for (size_t i=0; i<coeff_.size(); i++)
            {
                coeffArray[i] = coeff_[i].value;
            }
            return coeffArray; 
        }

        ~BasisEvaluator()
        {
        }
    };
} 
