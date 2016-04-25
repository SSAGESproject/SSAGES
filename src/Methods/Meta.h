#pragma once

#include "Method.h"
#include "../CVs/CollectiveVariable.h"
#include <fstream>
#include <vector>


namespace SSAGES
{
	// Structure representing a multidimensional hill (Gaussian)
	// which is centered at "center" with widths "width" of height
	// "height". A multidimensional Gaussian has one height but 
	// n centers and widths.
	struct Hill 
	{
		// Hill center.
		std::vector<double> center;

		// Hill width.
		std::vector<double> width;

		// Hill height.
		double height;
		
		// Constructs a multidimensional Hill (Gaussian)
		Hill(const std::vector<double>& center, 
			 const std::vector<double>& sigma, 
			 double height) :
		 center(center), width(sigma), height(height)
		{}
	};

	// Implementation of a "vanilla" multi-dimensional Metadynamics
	// method with no bells and whistles.
	class Meta : public Method
	{
	private:	
		// Hills.
		std::vector<Hill> _hills;

		// Hill height.
		double _height;

		// Hill widths.
		std::vector<double> _widths;

		// Derivatives.	
		std::vector<double> _derivatives;

		// Bias magnitude.
		double _bias;

		// Frequency of new hills
		unsigned int _hillfreq;

		// Adds a new hill.
		void AddHill(const CVList& cvs);

		// Computes the bias force.
		void CalcBiasForce(const CVList& cvs);

		// Prints the new hill to file
		void PrintHill(const Hill& hill);
		
		// Output stream for hill data.
		std::ofstream _hillsout;

	public: 
		// Constructs an instance of Metadynamics method.
		// "height" specifies the hieght of the hills to be deposited. 
		// "widths" specifies the widths of the hills to be deposited 
		// along each dimension. "hillfreq" specifies the frequency of 
		// depositing hills. Note that The size of "widths" should be 
		// commensurate with the number of CV's expected.
		Meta(boost::mpi::communicator& world,
			 boost::mpi::communicator& comm,
			 double height, 
			 const std::vector<double>& widths, 
			 unsigned int hillfreq, 
			 unsigned int frequency) : 
		Method(frequency, world, comm), _hills(), _height(height), _widths(widths), 
		_derivatives(0), _bias(0), _hillfreq(hillfreq)
		{
		}

		// Pre-simulation hook.
		void PreSimulation(Snapshot* snapshot, const CVList& cvs) override;

		// Post-integration hook.
		void PostIntegration(Snapshot* snapshot, const CVList& cvs) override;

		// Post-simulation hook.
		void PostSimulation(Snapshot* snapshot, const CVList& cvs) override;

		void Serialize(Json::Value& json) const override
		{

		}

		~Meta() {}
	};
}
			
