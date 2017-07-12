#include <cmath>
#include "../Grids/Grid.h"

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
        unsigned int poly_ord_; 
        unsigned int nbins_;
        double boundLow_, boundHigh_;
        bool isFinite_;

    public: 
        BasisFunction(unsigned int poly_ord,
                      unsigned int nbins,
                      bool isFinite,
                      double boundLow,
                      double boundHigh) :
        poly_ord_(poly_ord), isFinite_(isFinite),
        boundLow_(boundLow), boundHigh_(boundHigh)
        {
        }

        uint GetOrder() {return poly_ord_;}
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

        virtual double Evaluate(double val, int order){ return 0;}
        virtual double EvalGrad(double grad, int order){ return 0;}
        //virtual std::vector<double> ConvBounds(std::vector<double> x, double max, double min) = 0;
        //static BasisFunction* BuildBasisFunction(const Json::Value& json);
        virtual ~BasisFunction()
        {
        }
    };

    class Chebyshev : public BasisFunction
    {
    private: 
        
    protected:

    public:
        Chebyshev(unsigned int poly_ord, double boundLow, double boundHigh, unsigned int nbins) :
        BasisFunction(poly_ord, nbins, true, boundLow, boundHigh)
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
                    (2*n-1)/(double)n*(Evaluate(x,n-1) + x * EvalGrad(x,n-1)) - (n-1)/n*EvalGrad(x,n-2);
        }
    };

    class Legendre : public BasisFunction
    {
    private: 
        
    protected:

    public:
        Legendre(unsigned int poly_ord, unsigned int nbins) :
        BasisFunction(poly_ord, nbins, true, -1.0, 1.0)
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
                    (2*n-1)/(double)n*(Evaluate(x,n-1) + x * EvalGrad(x,n-1)) - (n-1)/n*EvalGrad(x,n-2);
        }
    };

    //Calculates the inner product of all the basis functions and the histogram
    class BasisEvaluator
    {
    private:
        std::vector<Map> coeff_;
        std::vector<BasisFunction> functions_;
        std::vector<BasisLUT> lookup_;

    public:
        //Initialize the evaluator
        BasisEvaluator(const std::vector<BasisFunction>& functions) : 
            functions_(functions)
        {
            CoeffInit();
            BasisInit();
        }
    
        //For now when the basis evaluator is called it will store all the values
        //of the basis functions into a lookup table. Subject to change
        void CoeffInit(void) {
            
            double coeff_size = 0;
            //Get the size of the number of coefficients + 1
            for(size_t i = 0; i < functions_.size(); ++i)
            {
                coeff_size *= functions_[i].GetOrder()+1;
            }

            std::vector<int> jdx(functions_.size(), 0);
            Map temp_map(jdx,0.0);
            //Initialize the mapping for the coeff function
            for(size_t i = 0; i < coeff_size; ++i)
            {
                for(size_t j = 0; j < jdx.size(); ++j)
                {
                    if(jdx[j] > 0 && jdx[j] % (functions_[j].GetOrder()+1) == 0)
                    {
                        if(j != functions_.size() - 1)
                        { 
                            jdx[j+1]++;
                            jdx[j] = 0;
                        }
                    }
                    temp_map.map[j] = jdx[j];
                    temp_map.value  = 0; 
                }
                coeff_.push_back(temp_map);           
                jdx[0]++;
            }
        }

        void BasisInit(void) {
            for (size_t i=0; i<functions_.size(); i++)
            {
                uint poly = functions_[i].GetOrder();
                uint nbin = functions_[i].GetBins();
                std::vector<double> dervs(nbin*poly,0);
                std::vector<double> vals(nbin*poly,0);

                for (size_t j=0; j<poly; j++) {
                    for(size_t k=0; k<nbin; k++) {
                        double x = ((functions_[i].GetRange())*k
                                - functions_[i].GetLower())/nbin + functions_[i].GetLower();
                        vals[i+j*nbin] = functions_[i].Evaluate(x,j);
                        dervs[i+j*nbin] = functions_[i].EvalGrad(x,j);
                    }
                }
                BasisLUT TempLUT(vals,dervs);
                lookup_.push_back(TempLUT);
            }
        }

        //Outputs the basis convolution at a specific coordinate
        void UpdateBias(Grid<double> bias, Grid<std::vector<double>> grad)
        {
            std::vector<double> tmp_grad (functions_.size(),0);
            double tmp_bias = 0;
            double basis;
            double temp;
            size_t j = 0;
            int nbins;

            for(Grid<double>::iterator it = bias.begin(); it != bias.end(); ++it, ++j)
            {
                for (size_t i = 1; i < coeff_.size(); ++i)
                {
                    for (size_t j = 0; j < functions_.size(); ++j)
                    {
                        temp = 1.0;
                        basis = 1.0;
                        for (size_t k = 0; k < functions_.size(); ++k)
                        {
                            nbins = bias.GetNumPoints(k);
                            //I don't know about this 2.0 here it might come from the Legendre basis
                            temp *= j == k ?  lookup_[k].derivs[it.coordinates()[k] + coeff_[i].map[k]*nbins] * 2.0 / (bias.GetUpper(k) - bias.GetLower(k))
                                           :  lookup_[k].values[it.coordinates()[k] + coeff_[i].map[k]*nbins];
                        }
                        tmp_grad[i] -= coeff_[i].value * temp;
                        basis *= lookup_[j].values[it.coordinates()[j] + coeff_[i].map[j]*nbins];
                    }
                    //Update the bias values
                    tmp_bias += coeff_[i].value * basis;
                    //Store the gradient values
                    grad.at(it.coordinates()) = tmp_grad;
                }
                *it = tmp_bias;
            }
        }

        //Calculates the inner product of the basis set and the biased histogram
        //This function then returns the coefficients from this calculation
        double UpdateCoeff(std::vector<double> h_bias, Grid<uint> hist)
        {
            double coeffTemp;
            double sum = 0;

            for(auto& coeff : coeff_)
            {
                // The method uses a standard integration with trap rule weights
                size_t j = 0;
                coeffTemp = coeff.value;
                for(Grid<uint>::iterator it2 = hist.begin(); it2 != hist.end(); ++it2, ++j)
                {
                    double weight = std::pow(2.0,functions_.size());
                    // This adds in a trap-rule type weighting which lowers error significantly at the boundaries
                    for(size_t k = 0; k < functions_.size(); ++k)
                    {
                        if( it2.index(k) == 0 ||
                            it2.index(k) == hist.GetNumPoints(k)-1)
                            weight /= 2.0;
                    }
                
                    /*The numerical integration of the biased histogram across the entirety of CV space
                     *All calculations include the normalization as well
                     */
                    double basis = 1.0;

                    for(size_t l = 0; l < functions_.size(); l++)
                    {
                        int nbins = hist.GetNumPoints(l);
                        basis *= lookup_[l].values[it2.index(l) + coeff.map[l]*nbins] / nbins;
                        basis *= (functions_[l].GetRange()) * coeff.map[l] + functions_[l].GetLower();
                    }
                    coeff.value += basis * log(h_bias[j]) * weight/std::pow(2.0,functions_.size());
                }
                coeffTemp -= coeff.value;
                sum += coeffTemp;
            }
            return sum;
        }

        //Gets the coefficient array
        std::vector<double> GetCoeff(void)
        {
            std::vector<double> coeff_array (coeff_.size(),0);
            for (size_t i=0; i<coeff_.size(); i++)
            {
                coeff_array[i] = coeff_[i].value;
            }
            return coeff_array; 
        }

        ~BasisEvaluator()
        {
        }
    };
}

            
