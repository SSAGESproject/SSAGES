#pragma once 

#include "../Simulations/SimObserver.h"
#include "json/json.h"
#include <iostream>
#include <fstream>
#include <memory>

namespace SSAGES
{
	//! Observer with JSON output
	/*!
	 * JSON Observer class that creates a snapshot of the simulation and dumps
	 * it to a JSON file.
	 *
	 * \ingroup Core
	 */
	class JSONObserver : public SimObserver
	{
	private:
		std::unique_ptr<std::ofstream> _jsonfs; //!< Output file pointer.
		std::string _prefix; //!< Prefix for output file.
		Json::Value _root; //!< JSON root value.
		bool _writetoother; //!< Indicate which json to write to.

	public:
		//! Constructor
		/*!
		 * \param prefix Prefix for output file.
		 * \param world MPI global communicator.
		 * \param comm MPI local communicator.
		 * \param numwalkers Number of walkers.
		 * \param wid ID of walker to observe.
		 * \param frequency Frequency of observations.
		 */
		JSONObserver(const std::string& prefix, 
						boost::mpi::communicator world,
						boost::mpi::communicator comm,
						int numwalkers,
						int wid,
						unsigned int frequency = 1);

		//! Get Observer name.
		/*!
		 * \return String "JSON"
		 */
		virtual std::string GetName() const override{ return "JSON"; }

		//! Get filename prefix.
		/*!
		 * \return Filename prefix.
		 */
		std::string GetPrefix() const {return _prefix;}

		//! Visit driver
		/*!
		 * \param d Driver to visit.
		 */
		virtual void Visit(const Driver& d) override;

		//! \copydoc Serializable::Serialize()
		virtual void Serialize(Json::Value& json) const override;
		
		//! Operations before visiting.
		virtual void PreVisit() override;

		//! Operations after visition.
		virtual void PostVisit() override;

		//! Destructor.
		~JSONObserver()
		{
			if(_jsonfs)
				_jsonfs->close();
		}
	};
}