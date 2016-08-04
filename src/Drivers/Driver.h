#pragma once

#include "../JSON/Serializable.h"
#include "json/json.h"
#include <vector>
#include <boost/mpi.hpp>
#include <random>
#include "../Hook.h"
#include "../Methods/Method.h"
#include "../Snapshot.h"
#include "../JSON/JSONLoader.h"
#include "../Grids/Grid.h"
#include "../Simulations/SimObservable.h"
#include "../Simulations/SimObserver.h"
#include "../Observers/Visitable.h"

namespace mpi = boost::mpi;
using namespace Json;
namespace SSAGES
{
	//! Abstract driver class for creating driver objects.
	/*!
	 * \ingroup Core
	 */
	class Driver : public SimObservable, public Serializable
	{

	protected:
		boost::mpi::communicator _world; //!< MPI global communicator
		boost::mpi::communicator _comm; //!< MPI local communicator

		//! The node id that this driver belongs to
		const int _wid;

		//! The hook that hooks into the MD simulation engine
		Hook* _hook;

		//! The snapshot of your system
		Snapshot* _snapshot;

		//! The Method that will be used
		Method* _method;

		//! The CVs that will be used
		CVList _CVs;

		//! The observers that will be used for this driver
		ObserverList _observers;

		//! Local input file
		std::string _inputfile;

		//! MD engine Restart file name
		std::string _restartname;

		//! Read a restart file
		bool _readrestart;

		//! Random number generators for setting seeds
		std::random_device _rd;

		//! Alternative random number generator
		std::mt19937 _gen;

	public:
		//! Constructor
		/*!
		 * \param world MPI global communicator.
		 * \param comm MPI local communicator.
		 * \param walkerID The walker node for this driver.
		 */
		Driver(boost::mpi::communicator& world, 
			   boost::mpi::communicator& comm,
			   int walkerID) : 
		_world(world), _comm(comm), _wid(walkerID),
		_hook(nullptr), _snapshot(), _method(), _CVs(),
		_observers(), _inputfile(), _restartname(),
		_readrestart(), _rd(), _gen(_rd())
		 {}

		//! Destructor
		virtual ~Driver()
		{
			for(auto& cv : _CVs)
				delete cv;

			for(auto& o : _observers)
				delete o;

			_CVs.clear();
			_observers.clear();

			delete _snapshot;

			delete _method;
		}

		//! Run simulation
		virtual void Run() = 0;

		//! Build the driver, which will create the hook and so forth
		/*!
		 * \param json JSON value containing input information.
		 * \param path Path for JSON path specification.
		 */
		virtual void BuildDriver(const Json::Value& json, const std::string& path) = 0;

		//! Read in driver specific input file (e.g. Lammps.in)
		/*!
		 * \param contents The contents of the input file.
		 *
		 * Execute the contents of an input file. The input file should be
		 * loaded with SSAGES::GetFileContents() first.
		 */
		virtual void ExecuteInputFile(std::string contents) = 0;

		//! Write a driver specific restart file
		/*!
		 * Write out a restart file using driver specific commands.
		 * This will most likely be executed by an observer to force
		 * a write of the restart file.
		 */
		virtual void WriteRestartFile() const = 0;

		//! Get the input file contents
		/*!
		 * \returns The contents of \c _inputfile.
		 *
		 * This function returns the contents of the input file.
		 *
		 * \note This function does not load the input file. If the input file
		 *       has not been loaded before, an empty string will be returned.
		 */
		std::string GetInputFile(){return _inputfile;}

		//! Set the input file
		/*!
		 * \param filename File name of the input file.
		 *
		 * This function sets the name of the input file for the driver.
		 */
		void SetInputFile(const std::string filename){ _inputfile = filename;}

		//! Get the Driver's method current iteration
		/*!
		 * \returns The drivers \c _method current iterator
		 *
		 * This function returns the Driver's methods iterator.
		 */
		int GetIteration(){return _method->GetIteration();}

		//! Build CVs
		/*!
		 * \param json JSON Value containing input information.
		 * \param path Path for JSON path specification.
		 */
		void BuildCVs(const Json::Value& json, const std::string& path)
		{
			CollectiveVariable::BuildCV(json, _CVs, path);
		}

		//! Build method(s).
		/*!
		 * \param json JSON value containing input information.
		 * \param path Path for JSON path specification.
		 */
		void BuildMethod(const Json::Value& json, const std::string& path)
		{
			_method = Method::BuildMethod(json, _world, _comm, path);
		}

		//! Build observer(s).
		/*!
		 * \param json JSON value containing input information.
		 * \param nwalks Number of walkers.
		 */
		void BuildObservers(const Json::Value& json,
			int nwalks)
		{
			SimObserver::BuildObservers(json, _world, _comm, nwalks, _wid, _observers);

			for(auto& o : _observers)
				this->AddObserver(o);
		}

		//! Build the grid.
		/*!
		 * \param json JSON value containing input information.
		 * \param path Path for JSON path specification.
		 */
		void BuildGrid(const Json::Value& json, const std::string& path)
		{
			_method->BuildGrid(json, path);
		}

		//! Finalize the setup
		/*!
		 * Create the snapshot and put all gathered values into the local hook.
		 * Set up listeners as well.
		 */
		void Finalize()
		{
			// Initialize snapshot. 
			_snapshot = new Snapshot(_comm, _wid);

			/* Remove this later? */
			if(_hook != nullptr)
			{
				// Set the hook to snapshot
				_hook->SetSnapshot(_snapshot);

				//Set the driver in the hook
				_hook->SetMDDriver(this);

				_hook->AddListener(_method);
				for(auto&cv : _CVs)
					_hook->AddCV(cv);
			}
		}

		//! \copydoc Serializable::Serialize()
		virtual void Serialize(Json::Value& json) const override
		{
			SerializeObservers(json["observers"]);

			_method->Serialize(json["method"]);

			auto* Grid = _method->GetGrid();

			if(Grid)
				Grid->Serialize(json["grid"]);

			auto& tmp = json["CVs"];
			for(unsigned int i = 0; i < _CVs.size();i++)
				_CVs[i]->Serialize(tmp[i]);

			json["inputfile"] = _inputfile;
			json["number processors"] = _comm.size();
			json["restart file"] = _restartname;
			json["read restart"] = _readrestart;
		}

     	// Accept a visitor.
		virtual void AcceptVisitor(Visitor& v) const override
		{
			v.Visit(*this);
		}
	};
}
