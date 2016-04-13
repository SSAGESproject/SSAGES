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
		double _alpha;

		// The node this belongs to
		unsigned int _mpiid;

		// The world's strings centers for each CV. _worldstring[cv#][node#]
		std::vector<std::vector<double> > _worldstring;

		// Time step of string change
		double _tau;

		// String modification parameter
		double _kappa;

		// Previous forces for restarting the position
		std::vector<Vector3> _prev_positions;

		// Number of nodes on a string
		unsigned int _numnodes;

		// Total iterations run so far
		unsigned int _currentiter;
		
		// Output stream for string data.
		std::ofstream _stringout;

		// Flag for when CV is outside newly paramaterized string voronoi cell
		bool _run_SMD;

		// Iterator to increase the spring stiffness on the umbrella method
		int _cv_inside_iterator;

		// Steered MD centers
		std::vector<double> _SMD_centers;

		// length to move the steered MD every time step
		std::vector<double> _SMD_lengths;

		// Updates the position of the string.
		void StringUpdate();

		// Prints the new hill to file
		void PrintString(const CVList& CV);

		// Check if current cv values are in the nodes voronoi cell
		bool InCell(const CVList& CV);

	public: 
		// Constructs an instance of Finite String method.
		// isteps = Number Iterations per block averaging
		// _tau and _kappa default values of 0.1 (JSON reader for this)
		FiniteTempString(boost::mpi::communicator& world,
					boost::mpi::communicator& com,
					unsigned int isteps,
					const std::vector<double>& centers,
					unsigned int NumNodes,
					double kappa,
					double tau,
			 		unsigned int frequency) : 
		Method(frequency, world, com), _blockiterations(isteps), _centers(centers), _cv_prev(), _alpha(),
		_mpiid(0), _worldstring(), _tau(tau), _kappa(kappa), _prev_positions(), _numnodes(NumNodes), _currentiter(0),
		_run_SMD(true), _cv_inside_iterator(0)
		{
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
			
