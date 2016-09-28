/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Hythem Sidky <hsidky@nd.edu>
 *                Ben Sikora <bsikora906@gmail.com>
 *
 * SSAGES is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SSAGES is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SSAGES.  If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once
#include <vector>
#include <memory>
#include <boost/mpi.hpp>
#include "json/json.h"
#include "JSON/Serializable.h"
#include "types.h"

namespace SSAGES
{
	inline double roundf(double x)
	{
    	return ( x >= 0 ) ? floor( x + 0.5 ) : ceil( x - 0.5 );
  	}

	//! Class containing a snapshot of the current simulation in time.
	/*!
	 * \ingroup core
	 *
	 * This contains information on the particle positions, velocities, etc...
	 * and additional information on the state of the system.
	 */
	class Snapshot : public Serializable
	{
	private:
		//! Local snapshot (walker) communicator.
		boost::mpi::communicator _comm;

		unsigned _wid; //!< Walker ID.
		unsigned _nlocal; //!< Number of atoms located on this snapshot

		std::string _ID; //!< ID string

		Matrix3 _H; //!< Parrinello-Rahman box H-matrix.
		Matrix3 _Hinv; //!< Parinello-Rahman box inverse.

		Vector3 _origin; //!< Box origin.

		Bool3 _isperiodic; //!< Periodicity of box.

		std::vector<Vector3> _positions; //!< Positions
		std::vector<Integer3> _images; //!< Unwrapped positions
		std::vector<Vector3> _velocities; //!< Velocities
		std::vector<Vector3> _forces; //!< Forces
		std::vector<double> _masses; //!< Masses
		std::vector<double> _charges; //!< Charges
		Label _atomids; //!< List of Atom IDs
		Label _types; //!< List of Atom types
		std::vector<std::vector<double> > _sigma; //!< Sigma

		int _iteration; //!< Iteration of Simulation
		double _temperature; //!< Current temperature
		double _energy; //!< Average per-particle energy
		double _kb; //!< Kb from the MD driver

		double _dielectric; //!< Dielectric
		double _qqrd2e; //!<qqrd2e
		bool _changed; //!< \c TRUE is Simulation state changed

	public:
		//! Constructor
		/*!
		 * \param comm MPI communicator
		 * \param wid Walker ID
		 *
		 * Initialize a snapshot with MPI communicator and a
		 * correpsonding walker ID.
		 */
		Snapshot(boost::mpi::communicator& comm, unsigned wid) :
		_comm(comm), _wid(wid), _H(), _Hinv(), _origin({0,0,0}), 
		_isperiodic({true, true, true}), _positions(0), _images(0), 
		_velocities(0), _forces(0), _masses(0), _atomids(0), _types(0), 
		_iteration(0), _temperature(0), _energy(0), _kb(0)
		{}

		//! Get the current iteration
		/*!
		 * \return Current Iteration
		 */
		int GetIteration() const {return _iteration; }

		//! Get current temperature
		/*!
		 * \return System temperature
		 */
		double GetTemperature() const {return _temperature; }

		//! Get per-particle energy
		/*!
		 * \return Current average per-particle energy
		 */
		double GetEnergy() const { return _energy; }

		//! Get system H-matrix
		/*! 
		 * \return Parrinello-Rahman H-matrix of simulation box
		 */
		const Matrix3& GetHMatrix() const { return _H; }

		//! Get origin of the system.
		/*! 
		 * \return Vector containing coordinates of box origin.
		 */
		const Vector3& GetOrigin() const { return _origin; }

		//! Get periodicity of three dimensions. 
		/*! 
		 * \return Three dimensional boolean containing periodicity of each dimension.
		 */
		const Bool3& IsPeriodic() const {return _isperiodic; }

		//! Get system volume
		/*!
		 * \return Volume of the current simulation box
		 */
		double GetVolume() const { return _H.determinant(); }

		//! Get system Kb
		/*!
		 * \return Kb of the current simulation box
		 */
		double GetKb() const { return _kb; }

		double GetDielectric() const { return _dielectric; }

		double Getqqrd2e() const { return _qqrd2e; }
		
		//! Get communicator for group (walker).
		/*!
		 * \return Communicator used for MetaDynamics.
		 *
		 * Access the communicator containing the set of processors that belong
		 * to a single simulation box.
		 */
		const boost::mpi::communicator& GetCommunicator() const
		{
			return _comm;
		}

		//! Get walker ID.
		/*!
		 * \return ID of the Walker
		 */
		unsigned GetWalkerID() const { return _wid; }


		//! Get number of atoms in this snapshot.
		/*!
		 * \return Number of atoms in this snapshot
		 */
		unsigned GetNumAtoms() const { return _nlocal; }

		//! Set the iteration
		/*!
		 * \param iteration New value for the iteration
		 */
		void SetIteration(int iteration) 
		{
			_iteration = iteration; 
			_changed = true; 
		}

		//! Change the temperature
		/*!
		 * \param temperature New value for the temperature
		 */
		void SetTemperature(double temperature) 
		{ 
			_temperature = temperature;
			_changed = true;
		}

		//! Change the energy
		/*!
		 * \param energy New value for the energy
		 */
		void SetEnergy(double energy) 
		{
			_energy = energy;
			_changed = true;
		}

		//! Change the Box H-matrix. 
		/*!
		 * \param hmat New H-matrix for the system
		*/
		void SetHMatrix(const Matrix3& hmat)
		{
			_H = hmat;
			_Hinv = hmat.inverse();
			_changed = true;
		}

		//! Change the box origin.
		/*!
		 * \param origin New origin for the system
		 */
		void SetOrigin(const Vector3& origin)
		{
			_origin = origin; 
			_changed = true;
		}

		//! Change the periodicity of the system
		/*!
		 * \param isperiodic Periodicity of three dimensions
		 */
		void SetPeriodicity(const Bool3& isperiodic)
		{
			_isperiodic = isperiodic;
			_changed = true;
		}

		//! Change the kb
		/*!
		 * \param kb New value for the kb
		 */
		void SetKb(double kb) 
		{
			_kb = kb;
			_changed = true;
		}
		
		void SetDielectric(double dielectric) 
		{ 
			_dielectric = dielectric;
			_changed = true;
		}
	
		void Setqqrd2e(double qqrd2e) 
		{ 
			_qqrd2e = qqrd2e;
			_changed = true;
		}

		//! Set number of atoms in this snapshot.
		/*!
		 * \param Number of atoms in this snapshot
		 */
		void SetNumAtoms(unsigned int natoms) { _nlocal = natoms; }

		//! Access the particle positions
		/*!
		 * \return List of particle positions
		 */
		const std::vector<Vector3>& GetPositions() const { return _positions; }

		/*! \copydoc Snapshot::GetPositions() const */
		std::vector<Vector3>& GetPositions() 
		{ 
			_changed = true;
			return _positions; 
		}

		//! Access the particles image flags
		/*!
		 * \return List of particle image flags
		 */
		const std::vector<Integer3>& GetImageFlags() const { return _images; }

		//! \copydoc Snapshot::GetImageFlags() const
		std::vector<Integer3>& GetImageFlags() 
		{ 
			_changed = true;
			return _images; 
		}

		//! Access the particle velocities
		/*!
		 * \return List of particle velocities
		 */
		const std::vector<Vector3>& GetVelocities() const { return _velocities; }

		/*! \copydoc Snapshot::GetVelocities() const */
		std::vector<Vector3>& GetVelocities() 
		{
			_changed = true;
			return _velocities; 
		}

		//! Access the per-particle forces
		/*!
		 * \return List of per-particle forces
		 */
		const std::vector<Vector3>& GetForces() const { return _forces; }

		/*! \copydoc Snapshot::GetForces() const */
		std::vector<Vector3>& GetForces() 
		{
			_changed = true; 
			return _forces; 
		}

		//! Const access to the particle masses
		/*!
		 * \return List of Masses
		 *
		 * Note that the Masses can be either stored as per-atom or per-type
		 * depending on the Lammps Atom type used.
		 */
		const std::vector<double>& GetMasses() const { return _masses; }

		/*! \copydoc Snapshot::GetMasses() const */
		std::vector<double>& GetMasses() 
		{
			_changed = true; 
			return _masses; 
		}

		//! Scale a vector into fractional coordinates
		/*! 
		 * \param v Vector of interest
		 * \return Scaled vector in fractional coordinates
		 */
		Vector3 ScaleVector(const Vector3& v) const
		{
			return _Hinv*(v-_origin);
		}

		//! Unwrap a vector's real coordinates according to its image replica count. 
		/*!
		 * \param v Vector of interest
		 * \param image Integer vector representing mirror images in the three dimensions.
		 * 
		 * \return Vector3 Unwrapped vector in real coordinates.
		 * This function takes a set of (wrapped) coordinates and returns the unwrapped
	     * coordinates. 
	     *
	     * \note This function does not require the initial coordinates to be within
	     *       the simulation box.
		 */
		Vector3 UnwrapVector(const Vector3& v, const Integer3& image) const
		{
			return _H*image.cast<double>()+v;
		}

		//! Apply minimum image to a vector.
		/*!
		 * \param v Vector of interest
		 */
		void ApplyMinimumImage(Vector3* v) const
		{
			Vector3 scaled = _Hinv*(*v);

			for(int i = 0; i < 3; ++i)
				scaled[i] -= _isperiodic[i]*roundf(scaled[i]);

			*v = _H*scaled;
		}
		
		//! Apply minimum image to a vector.
		/*!
		 * \param v Vector of interest
		 */
		Vector3 ApplyMinimumImage(const Vector3& v) const
		{
			Vector3 scaled = _Hinv*v;

			for(int i = 0; i < 3; ++i)
				scaled[i] -= _isperiodic[i]*roundf(scaled[i]);

			return _H*scaled;	
		}

		//! Compute the total mass of a group of particles based on index. 
		/*!
		  * \param indices IDs of particles of interest. 
		  * \return double Total mass of particles. 
		  *
		  * \note This function reduces the mass across the communicator associated
		  *       with the snapshot. 
		 */
		double TotalMass(const Label& indices) const
		{
			auto mtot = 0.;
			for(auto& i : indices)
				mtot += _masses[i];
			MPI_Allreduce(MPI_IN_PLACE, &mtot, 1, MPI_DOUBLE, MPI_SUM, _comm);
			return mtot;
		}

		//! Compute center of mass of a group of atoms based on idex with provided 
		//! Total mass.
		/*!
		  * \param indices IDs of particles of interest. 
		  * \return Vector3 Center of mass of particles.
		  */ 		
		Vector3 CenterOfMass(const Label& indices) const
		{
			// Get total mass.
			auto mtot = TotalMass(indices);

			return CenterOfMass(indices, mtot);			
		}

		//! Compute center of mass of a group of atoms based on idex with provided 
		//! Total mass.
		/*!
		  * \param indices IDs of particles of interest. 
		  * \param mtot Total mass of particle group. 
		  * \return Vector3 Center of mass of particles.
		  * \note Each processor passes in the local indices of the atoms of interest
		  *       and this function will collect the data and compute the center of mass.
		  */ 
		Vector3 CenterOfMass(const Label& indices, double mtot) const
		{
			// Store coorinates and masses in vectors to gather. 
			std::vector<double> pos, mass, gpos, gmass;
			std::vector<int> pcounts(_comm.size(), 0), mcounts(_comm.size(), 0); 
			std::vector<int> pdispls(_comm.size()+1, 0), mdispls(_comm.size()+1, 0);

			pcounts[_comm.rank()] = 3*indices.size();
			mcounts[_comm.rank()] = indices.size();

			// Reduce counts.
			MPI_Allreduce(MPI_IN_PLACE, pcounts.data(), pcounts.size(), MPI_INT, MPI_SUM, _comm);
			MPI_Allreduce(MPI_IN_PLACE, mcounts.data(), mcounts.size(), MPI_INT, MPI_SUM, _comm);

			// Compute displacements.
			std::partial_sum(pcounts.begin(), pcounts.end(), pdispls.begin() + 1);
			std::partial_sum(mcounts.begin(), mcounts.end(), mdispls.begin() + 1);
			
			// Fill up mass and position vectors.
			for(auto& idx : indices)
			{
				auto& p = _positions[idx];
				pos.push_back(p[0]);
				pos.push_back(p[1]);
				pos.push_back(p[2]);
				mass.push_back(_masses[idx]);
			}

			// Re-size receiving vectors. 
			gpos.resize(pdispls.back(), 0);
			gmass.resize(mdispls.back(), 0);

			// All-gather data.
			MPI_Allgatherv(pos.data(), pos.size(), MPI_DOUBLE, gpos.data(), pcounts.data(), pdispls.data(), MPI_DOUBLE, _comm);
			MPI_Allgatherv(mass.data(), mass.size(), MPI_DOUBLE, gmass.data(), mcounts.data(), mdispls.data(), MPI_DOUBLE, _comm);

			// Loop through atoms and compute mass weighted sum. 
			// We march linearly through list and find nearest image
			// to each successive particle to properly unwrap object.
			Vector3 ppos = {gpos[0], gpos[1], gpos[2]}; // Previous unwrapped position.
			Vector3 cpos = ppos;
			Vector3 xcm = gmass[0]*cpos;

			for(size_t i = 1, j = 3; i < gmass.size(); ++i, j += 3)
			{
				cpos = {gpos[j], gpos[j+1], gpos[j+2]};
				cpos = ApplyMinimumImage(cpos - ppos) + ppos;
				xcm += gmass[i]*cpos;
				ppos = cpos;
			}

			return xcm/mtot;
		}

		//! Access the atom IDs
		/*!t
		 * \return List of atom IDs
		 */
		const Label& GetAtomIDs() const { return _atomids; }

		/*! \copydoc Snapshot::GetAtomIDs() const */
		Label& GetAtomIDs()
		{
			_changed = true;
			return _atomids;
		}

		//! Gets the local atom index corresponding to an atom ID.
		/*!
		 * \param Atom ID. 
		 * \return Local atom index or -1 if not found.
		 */
		int GetLocalIndex(int id) const
		{

			auto s = std::find(_atomids.begin(), _atomids.end(), id);
			if(s == _atomids.end())
				return -1;
			else
				return s - _atomids.begin();
		}

		//! Gets the local atom indices corresponding to atom IDs in the 
		//! vector.
		/*!
		 * \param ids Vector of atom ID's. 
		 * \param indices Pointer to container for local atom indices. 
		 *
		 * \note If atom does not exist on processor, it will be ignored. 
		 */
		void GetLocalIndices(const Label& ids, Label* indices) const
		{
			for(auto& id : ids)
			{
				auto idx = GetLocalIndex(id);
				if(idx != -1)
					indices->push_back(idx);
			}
		}

		//! Access the atom charges
		/*!
		 * \return List of atom charges
		 */
		const std::vector<double>& GetCharges() const { return _charges; }
		std::vector<double>& GetCharges()
		{
			_changed = true;
			return _charges;
		}
		
		//! Access the atom types
		/*!
		 * \return List of atom types
		 */
		const Label& GetAtomTypes() const { return _types; }

		/*! \copydoc Snapshot::GetAtomTypes() const */
		Label& GetAtomTypes() 
		{
			_changed = true; 
			return _types; 
		}

		//! Access the atom sigmas
		/*!
		 * \return List of atom sigmas
		 */
		const std::vector<std::vector<double>>& GetSigmas() const { return _sigma; }
		std::vector<std::vector<double>>& GetSigmas() 
		{
			_changed = true; 
			return _sigma;
		}

		//! Access the snapshot ID
		/*!
		 * \return Snapshot ID
		 */
		const std::string& GetSnapshotID() const { return _ID; }

		/*! \copydoc Snapshot::GetSnapshotID() const */
		std::string& GetSnapshotID() 
		{
			_changed = true; 
			return _ID; 
		}

		//! Query if Snapshot was modified
		/*!
		 * \return \c True if Snapshot was modified, else return \c False
		 */
		bool HasChanged() const { return _changed; }

		//! Set the "changed" flag of the Snapshot
		/*!
		 * \param state State to which the "changed" flag is set
		 */
		void Changed(bool state) { _changed = state; }

		//! Serialize the Snapshot
		/*!
		 * \param json Json value to write serialized state to
		 */
		void Serialize(Json::Value& json) const override
		{
			json["walker id"] = _wid;
			json["id"] = _ID;

			for(unsigned int i = 0; i < _atomids.size(); i++)
			{
				json["Atom"][i]["ID"] = _atomids[i];
				json["Atom"][i]["type"] = _types[i];
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
			json["energy"] = _energy;
			json["kb"] = _kb;

			for(int i = 0; i < 3; ++i)
				for(int j = 0; j < 3; ++j)
					json["hmat"][i][j] = _H(i,j);
		}
		
		//! Destructor
		~Snapshot(){}
	};
}

