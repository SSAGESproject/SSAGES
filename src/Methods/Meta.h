#pragma once

#include "Method.h"
#include "../CVs/CollectiveVariable.h"
#include <fstream>

namespace SSAGES
{
	struct Hill 
	{
		std::vector<double> center;
		std::vector<double> width;
		double height;
	
		Hill(const std::vector<double>& center, 
			 const std::vector<double>& sigma, 
			 double height) :
		 center(center), width(sigma), height(height)
		{}
	};

	class Meta : public Method
	{
	private:		
		// Contains values of the CV's.
		std::vector<double> _cvs;

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
		void CalcBiasForce();

		// Prints the new hill to file
		void printHill ();
		std::ofstream hillsout;

	public: 
		Meta(double height, 
			 const std::vector<double>& widths, 
			 unsigned int hillfreq, 
			 unsigned int frequency) : 
		Method(frequency), _cvs(0), _hills(), _height(height), _widths(widths), 
		_derivatives(0), _bias(0), _hillfreq(hillfreq)
		{
		}

		// Pre-simulation hook.
		void PreSimulation(Snapshot* snapshot, const CVList& cvs) override;

		// Post-integration hook.
		void PostIntegration(Snapshot* snapshot, const CVList& cvs) override;

		// Post-simulation hook.
		void PostSimulation(Snapshot* snapshot, const CVList& cvs) override;

		~Meta() {}
	};
}
			
