#pragma once
#include <vector>

namespace SSAGES
{
	using Vector = std::vector<double>;
	using Label = std::vector<int>;
	// Class containing a snapshot of the current simulation in time. 
	// This contains information on the particle positions, velocities, etc.. 
	// and additional information on the state of the system.
	class Snapshot
	{
	private:
		Vector _positions;
		Vector _velocities;
		Vector _forces;
		Label _atomids;
		Label _molids;
		Label _types;

		int _iteration; 
		double _time; 
		double _temperature; 
		double _pressure; 
		double _energy;
		double _volume;

		bool _changed;

	public:
		Snapshot() :
		_positions(0), _velocities(0), _forces(0),
		_atomids(0), _molids(0), _types(0), _iteration(0), 
		_time(0), _temperature(0), _pressure(0), _energy(0),
		_volume(0)
		{}

		int GetIteration() const {return _iteration; }
		double GetTime() const { return _time; }
		double GetTemperature() const {return _temperature; }
		double GetPressure() const { return _pressure; }
		double GetEnergy() const { return _energy; }
		double GetVolume() const { return _volume; }

		const Vector& GetPositions() const { return _positions; }
		Vector& GetPositions() 
		{ 
			_changed = true;
			return _positions; 
		}

		const Vector& GetVelocities() const { return _velocities; }
		Vector& GetVelocities() 
		{
			_changed = true;
			return _velocities; 
		}

		const Vector& GetForces() const { return _forces; }
		Vector& GetForces() 
		{
			_changed = true; 
			return _forces; 
		}

		const Label& GetAtomIDs() const { return _atomids; }

		const Label& GetMoleculeIDs() const { return _molids; }
		
		const Label& GetAtomTypes() const { return _types; }

		bool HasChanged() const { return _changed; }
		void Changed(bool state) { _changed = state; }
		
		~Snapshot(){}		
	};
}