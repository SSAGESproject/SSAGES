#pragma once 

#include "Method.h"
#include "../CVs/CollectiveVariable.h"
#include <fstream>

namespace SSAGES
{
	//! Umbrella sampling method.
	/*!
	 * Umbrella sampling method to constrain an arbitrary number of CVs at
	 * specified equilibrium distances.
	 *
	 * \ingroup Methods
	 */
	class Umbrella : public Method
	{
	private:
		//! Vector of spring constants.
		std::vector<double> _kspring;

		//! Vector of equilibrium distances.
		std::vector<double> _centers;

		//! Iterator for this method
		int _currentiter;

		//! Output stream for umbrella data.
		std::ofstream _umbrella;

	public:
		//! Constructor.
		/*!
		 * \param world MPI global communicator.
		 * \param comm MPI local communicator.
		 * \param kspring List of spring constants.
		 * \param centers List of spring equilibrium positions.
		 * \param frequency Frequency with which this method is applied.
		 *
		 * Create instance of umbrella with spring constants "kspring", and
		 * centers "centers". Note the sizes of the vectors should be
		 * commensurate with the number of CVs.
		 */
		Umbrella(boost::mpi::communicator& world,
				 boost::mpi::communicator& comm,
				 const std::vector<double>& kspring,
				 const std::vector<double>& centers,
				 unsigned int frequency) : 
		Method(frequency, world, comm), _kspring(kspring), _centers(centers)
		{}

		//! Pre-simulation hook.
		/*!
		 * \param snapshot Simulation snapshot.
		 * \param cvs List of CVs.
		 */
		void PreSimulation(Snapshot* snapshot, const CVList& cvs) override;

		//! Post-integration hook.
		/*!
		 * \param snapshot Simulation snapshot.
		 * \param cvs List of CVs.
		 */
		void PostIntegration(Snapshot* snapshot, const CVList& cvs) override;

		//! Post-simulation hook.
		/*!
		 * \param snapshot Simulation snapshot.
		 * \param cvs List of CVs.
		 */
		void PostSimulation(Snapshot* snapshot, const CVList& cvs) override;
		
		//! Print umbrella values.
		/*!
		 * \param cvs List of CVs.
		 */
		void PrintUmbrella(const CVList& cvs);

		//! \copydoc Serializable::Serialize()
		/*!
		 * \warning The serialization is not implemented yet.
		 */
		void Serialize(Json::Value& json) const override
		{

		}

	};
}