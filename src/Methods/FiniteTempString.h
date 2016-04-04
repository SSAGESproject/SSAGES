#include "Method.h"
#include "../CVs/CollectiveVariable.h"
#include <fstream>

namespace SSAGES
{
	// Implementation of a multi-walker finite string
	// method with hard wall voronoi cells and running block averages.
	class FiniteTempString : public Method
	{
	private:	
		
		// Number of steps to block average the CV's postions over
		unsigned int _blockiterations;

		// The local method iterator
		unsigned int _iterator;

		// Running averages
		std::vector<double> _runavgs;

		// CV starting location values
		std::vector<double> _centers;

		// Vector of CV values at prevoius step
		std::vector<double> _cv_prev;

		// for reparameterization
		std::vector<double> _alpha;

		// The node this belongs to
		unsigned int _mpiid;

		// The world's strings centers for each CV
		std::vector<std::vector<double> > _worldstring;

		// Time step of string change
		double _tau;

		// String modification parameter
		double _kappa;

		// Previous forces for restarting the position
		std::vector<Vector3> _prev_positions;
		
		// Output stream for string data.
		std::ofstream _stringout;

		// Updates the position of the string.
		void StringUpdate();

		// Prints the new hill to file
		void PrintString();

	public: 
		// Constructs an instance of Finite String method.
		// isteps = Number Iterations per block averaging
		// _tau and _kappa default values of 0.1 (JSON reader for this)
		FiniteTempString(unsigned int isteps,
					const std::vector<double>& centers,
					int NumNodes,
					double kappa,
					double tau,
			 		unsigned int frequency) : 
		Method(frequency), _blockiterations(isteps), _centers(centers), _cv_prev(), _alpha(),
		_mpiid(0), _worldstring(), _tau(tau), _kappa(kappa), _prev_positions()
		{
			_worldstring.resize(NumNodes);
		}

		// Pre-simulation hook.
		void PreSimulation(Snapshot* snapshot, const CVList& cvs) override;

		// Post-integration hook.
		void PostIntegration(Snapshot* snapshot, const CVList& cvs) override;

		// Post-simulation hook.
		void PostSimulation(Snapshot* snapshot, const CVList& cvs) override;

		~FiniteTempString() {}
	};
}
			
