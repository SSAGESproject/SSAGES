#pragma once
#include <vector>
#include <memory>

namespace SSAGES
{
	using Vector3 = std::array<double, 3>;
	using Label = std::vector<int>;
	// Class containing a snapshot of the current simulation in time. 
	// This contains information on the particle positions, velocities, etc.. 
	// and additional information on the state of the system.
	class Snapshot
	{
	private:
		Vector3 _positions;
		Vector3 _velocities;
		Vector3 _forces;
		Label _atomids;
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
		_positions({{0,0,0}}), _velocities({{0,0,0}}), _forces({{0,0,0}}),
		_atomids(0), _types(0), _iteration(0), 
		_time(0), _temperature(0), _pressure(0), _energy(0),
		_volume(0)
		{}

		int GetIteration() const {return _iteration; }
		double GetTime() const { return _time; }
		double GetTemperature() const {return _temperature; }
		double GetPressure() const { return _pressure; }
		double GetEnergy() const { return _energy; }
		double GetVolume() const { return _volume; }

		void SetIteration(int iteration) 
		{
			_iteration = iteration; 
			_changed = true; 
		}

		void SetTime(double time) 
		{
			_time = time; 
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

		const Vector3& GetPositions() const { return _positions; }
		Vector3& GetPositions() 
		{ 
			_changed = true;
			return _positions; 
		}

		const Vector3& GetVelocities() const { return _velocities; }
		Vector3& GetVelocities() 
		{
			_changed = true;
			return _velocities; 
		}

		const Vector3& GetForces() const { return _forces; }
		Vector3& GetForces() 
		{
			_changed = true; 
			return _forces; 
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

		bool HasChanged() const { return _changed; }
		void Changed(bool state) { _changed = state; }
		
		~Snapshot(){}		
	};
}