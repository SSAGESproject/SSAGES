#pragma once
#include <vector>
#include <memory>
#include <boost/mpi.hpp>
#include "json/json.h"
#include "JSON/Serializable.h"

namespace SSAGES
{
	using Vector3 = std::array<double, 3>;
	using Label = std::vector<int>;
	
	// Class containing a snapshot of the current simulation in time. 
	// This contains information on the particle positions, velocities, etc.. 
	// and additional information on the state of the system.
	class Snapshot : public Serializable
	{
	private:
		// Local snapshot (walker) communicator.
		boost::mpi::communicator _comm;

		// Walker ID.
		unsigned _wid;

		std::string _ID;

		std::vector<Vector3> _positions;
		std::vector<Vector3> _velocities;
		std::vector<Vector3> _forces;
		std::vector<double> _masses;
		Label _atomids;
		Label _types;

		int _iteration; 
		double _temperature; 
		double _pressure; 
		double _energy;
		double _volume;

		bool _changed;

	public:
		// Initialize a snapshot with MPI communicator and a
		// correpsonding walker ID.
		Snapshot(boost::mpi::communicator& comm, unsigned wid) :
		_comm(comm), _wid(wid), _positions(0), _velocities(0), 
		_forces(0), _atomids(0), _types(0), 
		_iteration(0), _temperature(0), _pressure(0), 
		_energy(0), _volume(0)
		{}

		int GetIteration() const {return _iteration; }
		double GetTemperature() const {return _temperature; }
		double GetPressure() const { return _pressure; }
		double GetEnergy() const { return _energy; }
		double GetVolume() const { return _volume; }
		
		// Get communicator for group (walker).
		// This contains the set of processors that belong to a 
		// single simulation box.
		const boost::mpi::communicator& GetCommunicator() const
		{
			return _comm;
		}

		// Get walker ID.
		unsigned GetWalkerID() const { return _wid; }

		void SetIteration(int iteration) 
		{
			_iteration = iteration; 
			_changed = true; 
		}

		void SetTemperature(double temperature) 
		{ 
			_temperature = temperature;
			_changed = true;
		}

		void SetPressure(double pressure) 
		{ 
			_pressure = pressure;
			_changed = true;
		}
		void SetEnergy(double energy) 
		{
			_energy = energy;
			_changed = true;
		}
		void SetVolume(double volume) 
		{
			_volume = volume;
			_changed = true;
		}

		const std::vector<Vector3>& GetPositions() const { return _positions; }
		std::vector<Vector3>& GetPositions() 
		{ 
			_changed = true;
			return _positions; 
		}

		const std::vector<Vector3>& GetVelocities() const { return _velocities; }
		std::vector<Vector3>& GetVelocities() 
		{
			_changed = true;
			return _velocities; 
		}

		const std::vector<Vector3>& GetForces() const { return _forces; }
		std::vector<Vector3>& GetForces() 
		{
			_changed = true; 
			return _forces; 
		}


		const std::vector<double>& GetMasses() const { return _masses; }
		std::vector<double>& GetMasses() 
		{
			_changed = true; 
			return _masses; 
		}

		const Label& GetAtomIDs() const { return _atomids; }
	 	Label& GetAtomIDs()
	 	{ 
	 		_changed = true;
	 		return _atomids; 
	 	}
		
		const Label& GetAtomTypes() const { return _types; }
		Label& GetAtomTypes() 
		{
			_changed = true; 
			return _types; 
		}

		const std::string& GetSnapshotID() const { return _ID; }
		std::string& GetSnapshotID() 
		{
			_changed = true; 
			return _ID; 
		}


		bool HasChanged() const { return _changed; }
		void Changed(bool state) { _changed = state; }

		void Serialize(Json::Value& json) const override
		{
			json["walker id"] = _wid;
			json["id"] = _ID;

			for(int i = 0; i<_atomids.size(); i++)
			{
				json["Atom"][i]["ID"] = _atomids[i];
				json["Atom"][i]["type"] = _types[i];
				json["Atom"][i]["mass"] = _masses[i];
				json["Atom"][i]["mass"] = _masses[i];
				for(int j = 0; j<3; j++)
				{
					json["Atom"][i]["positions"][j] = _positions[i][j];
					json["Atom"][i]["velocities"][j] = _velocities[i][j];
					json["Atom"][i]["forces"][j] = _forces[i][j];
				}
			}
			
			json["iteration"] = _iteration; 
			json["temperature"] = _temperature; 
			json["pressure"] = _pressure; 
			json["energy"] = _energy;
			json["volume"] = _volume;
		}
		
		~Snapshot(){}		
	};
}