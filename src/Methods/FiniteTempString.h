#include "Method.h"
#include "../CVs/CollectiveVariable.h"
#include <fstream>

namespace SSAGES
{
	//! Finite Temperature Spring Method
	/*!
	 * \ingroup Methods
	 *
	 * Implementation of a multi-walker finite string
	 * method with hard wall voronoi cells and running block averages.
	 */
	class FiniteTempString : public Method
	{
	private:	
		
		//! Number of steps to block average the CV's postions over
		unsigned int _blockiterations;

		//! The local method iterator
		unsigned int _iterator;

		//! Running averages
		std::vector<double> _runavgs;

		//! CV starting location values
		std::vector<double> _centers;

		//! Vector of CV values at prevoius step
		std::vector<double> _cv_prev;

		//! for reparameterization
		double _alpha;

		//! The node this belongs to
		unsigned int _mpiid;

		//! The world's strings centers for each CV.
		/*!
		 * _worldstring[cv#][node#]
		 */
		std::vector<std::vector<double> > _worldstring;

		//! Time step of string change
		double _tau;

		//! String modification parameter
		double _kappa;

		//! Spring constant for umbrella potential for staying within new voronoi cell after each update
		double _spring;

		//! Tolerance criteria for determining when to stop string (default 0 if no tolerance criteria)
		double _tol;

		//! Previous iteration's CV centers for checking tolerance criteria
		std::vector<std::vector<double> > _tolcheck;

		//! Previous forces for restarting the position
		std::vector<Vector3> _prev_positions;

		//! Number of nodes on a string
		unsigned int _numnodes;
		
		//! Output stream for string data.
		std::ofstream _stringout;

		//! Flag for when CV is outside newly paramaterized string voronoi cell
		bool _run_SMD;

		//! Iterator to increase the spring stiffness on the umbrella method
		int _cv_inside_iterator;

		//! Steered MD centers
		std::vector<double> _SMD_centers;

		//! Updates the position of the string.
		void StringUpdate();

		//! Prints the new hill to file
		void PrintString(const CVList& CV);

		//! Check if current cv values are in the nodes voronoi cell
		bool InCell(const CVList& CV);

		//! Check whether tolerance criteria has been met
		bool TolCheck();

		//! Maximum cap on number of string method iterations performed
		unsigned int _maxiterator;

		//! Options for restarting string
		bool _restart;

		//! Averages of restart configuration.
		std::vector<double> _restartavgs;

	public:
		//! Constructor
		/*!
		 * \param world MPI global communicator.
		 * \param comm MPI local communicator.
		 * \param isteps Number of iterations per block averaging.
		 * \param centers List of centers.
		 * \param NumNodes Number of nodes.
		 * \param kappa Value of kappa (default: 0.1).
		 * \param tau Value of tau (default: 0.1).
		 * \param spring Spring constant.
		 * \param tol Tolerance value.
		 * \param maxiterator Maximum number of iterations.
		 * \param restart If \c True start from restart file.
		 * \param restartiter Iteration of restart system.
		 * \param restartavgs Averages of restart system.
		 * \param frequency Frequency with which this method is invoked.
		 *
		 * Constructs an instance of Finite String method.
		 * isteps = Number Iterations per block averaging
		 * _tau and _kappa default values of 0.1 (JSON reader for this)
		 */
		FiniteTempString(boost::mpi::communicator& world,
					boost::mpi::communicator& comm,
					unsigned int isteps,
					const std::vector<double>& centers,
					unsigned int NumNodes,
					double kappa,
					double tau,
					double spring,
					double tol,
					unsigned int maxiterator,
					bool restart,
					const std::vector<double>& restartavgs,
			 		unsigned int frequency) : 
		Method(frequency, world, comm), _blockiterations(isteps), _centers(centers),
		_cv_prev(), _alpha(), _mpiid(0), _worldstring(), _tau(tau), _kappa(kappa),
		_spring(spring), _tol(tol), _prev_positions(), _numnodes(NumNodes),
        _run_SMD(true), _cv_inside_iterator(0),
		_maxiterator(maxiterator), _restart(restart), _restartavgs(restartavgs)
		{
		}

		//! Pre-simulation hook.
		void PreSimulation(Snapshot* snapshot, const CVList& cvs) override;

		//! Post-integration hook.
		void PostIntegration(Snapshot* snapshot, const CVList& cvs) override;

		//! Post-simulation hook.
		void PostSimulation(Snapshot* snapshot, const CVList& cvs) override;

        void SetIteration(const int iter)
        {
            _iteration = iter;
        }
        
		void Serialize(Json::Value& json) const override
        {
            json["type"] = "FTS";
            for(size_t i = 0; i < _centers.size(); i++)
            {
                json["centers"].append(_centers[i]);
            }

            json["number samples"] = _numnodes;
            json["spring"] = _spring;
            json["block iterations"] = _blockiterations;
            json["kappa"] = _kappa;
            json["time step"] = _tau;
            json["tol"] = _tol;
            json["max iterations"] = _maxiterator;
            json["restart"] = true; 

            for(size_t i = 0; i < _restartavgs.size(); i++)
            {
                json["previous avgs"].append(_restartavgs[i]);
            }

            json["iteration"] = _iteration;
        }

		//! Destructor
		~FiniteTempString() {}
	};
}
			
