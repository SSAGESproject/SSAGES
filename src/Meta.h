#pragma once

#include "Method.h"
#include "CollectiveVariable.h"

namespace SSAGES
{
	struct Hill {
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

		// Hill widths.
		std::vector<double> _widths;

		// Hill height.
		double _height;

		// Derivatives.	
		std::vector<double> _derivatives;

		// Bias magnitude.
		double _bias;

		// Frequency of new hills
		double _hillfreq;

		// Adds a new hill.
		void AddHill(const CVList& cvs);

		// Computes the bias force.
		void CalcBiasForce();

		// Computes chain rule.
		void ChainRule(const CVList& cvs);

	public: 
		Meta(unsigned int frequency) : 
		Method(frequency)
		{
		}

		// Pre-simulation hook.
		void PreSimulation(Snapshot* snapshot, const CVList& cvs) override;

		// Post-integration hook.
		void PostIntegration(Snapshot* snapshot, const CVList& cvs) override;

		~Meta() {}
	};
}
			
