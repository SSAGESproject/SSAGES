/**
 * This file is part of
 * SSAGES - Software Suite for Advanced General Ensemble Simulations
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
#include <numeric>
#include <vector>
#include <memory>
#include <mxx/comm.hpp>
#include "types.h"

namespace SSAGES
{
	//! Quick helper function to round a double.
	/*!
	 * \param x Double to be rounded.
	 * \return Rounded number as a double.
	 */
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
	class Snapshot
	{
	private:
		//! Local snapshot (walker) communicator.
		mxx::comm comm_;

		unsigned int wid_; //!< Walker ID.
		unsigned int nlocal_; //!< Number of atoms located on this snapshot

		std::string ID_; //!< ID string

		Matrix3 H_; //!< Parrinello-Rahman box H-matrix.
		Matrix3 Hinv_; //!< Parinello-Rahman box inverse.

		Matrix3 virial_; //!< Virial tensor.

		Vector3 origin_; //!< Box origin.

		Bool3 isperiodic_; //!< Periodicity of box.

		std::vector<Vector3> positions_; //!< Positions
		std::vector<Vector3> velocities_; //!< Velocities
		std::vector<Vector3> forces_; //!< Forces
		std::vector<double> masses_; //!< Masses
		std::vector<double> charges_; //!< Charges
		Label atomids_; //!< List of Atom IDs
		Label types_; //!< List of Atom types

		size_t iteration_; //!< Iteration of Simulation
		size_t targetiter_; //!< Iteration target of simulation.
		double temperature_; //!< Current temperature
		double energy_; //!< Average per-particle energy
		double kb_; //!< Kb from the MD driver

		double dielectric_; //!< Dielectric
		double qqrd2e_; //!<qqrd2e
		bool changed_; //!< \c TRUE is Simulation state changed

	public:
		//! Constructor
		/*!
		 * \param comm MPI communicator
		 * \param wid Walker ID
		 *
		 * Initialize a snapshot with MPI communicator and a
		 * correpsonding walker ID.
		 */
		Snapshot(const MPI_Comm& comm, unsigned int wid) :
		comm_(comm), wid_(wid), H_(), Hinv_(), virial_(Matrix3::Zero()), 
		origin_({0,0,0}), isperiodic_({true, true, true}), positions_(0), 
		velocities_(0), forces_(0), masses_(0), atomids_(0), types_(0),
		iteration_(0), temperature_(0), energy_(0), kb_(0)
		{}

		//! Get the current iteration
		/*!
		 * \return Current Iteration
		 */
		size_t GetIteration() const {return iteration_; }
		
		//! Get target iterations
		/*!
		 * \return Target Iterations
		 */
		size_t GetTargetIterations() const {return targetiter_; }

		//! Get current temperature
		/*!
		 * \return System temperature
		 */
		double GetTemperature() const {return temperature_; }

		//! Get per-particle energy
		/*!
		 * \return Current average per-particle energy
		 */
		double GetEnergy() const { return energy_; }

		//! Get system H-matrix
		/*! 
		 * \return Parrinello-Rahman H-matrix of simulation box
		 */
		const Matrix3& GetHMatrix() const { return H_; }

		//! Get box virial
		/*!
		 * \return Virial tensor of the simulation box. 
		 */
		const Matrix3& GetVirial() const { return virial_; }

		//! Get box virial
		/*!
		 * \return Virial tensor of the simulation box. 
		 */
		Matrix3& GetVirial() { return virial_; }

		//! Get origin of the system.
		/*! 
		 * \return Vector containing coordinates of box origin.
		 */
		const Vector3& GetOrigin() const { return origin_; }

		//! Get periodicity of three dimensions. 
		/*! 
		 * \return Three dimensional boolean containing periodicity of each dimension.
		 */
		const Bool3& IsPeriodic() const {return isperiodic_; }

		//! Get system volume
		/*!
		 * \return Volume of the current simulation box
		 */
		double GetVolume() const { return H_.determinant(); }

		//! Get system Kb
		/*!
		 * \return Kb of the current simulation box
		 */
		double GetKb() const { return kb_; }

		//! Get dielectric constant.
		/*!
		 * \return dielectric constant.
		 */
		double GetDielectric() const { return dielectric_; }

		//! Get qqrd2e value.
		/*!
		 * \return Value of qqrd2e.
		 */
		double Getqqrd2e() const { return qqrd2e_; }
		
		//! Get communicator for walker.
		/*!
		 * \return Communicator containing processors for simulation instance (walker).
		 *
		 * Access the communicator containing the set of processors that belong
		 * to a single simulation box.
		 */
		const mxx::comm& GetCommunicator() const
		{
			return comm_;
		}

		//! Get walker ID.
		/*!
		 * \return ID of the Walker
		 */
		unsigned GetWalkerID() const { return wid_; }


		//! Get number of atoms in this snapshot.
		/*!
		 * \return Number of atoms in this snapshot
		 */
		unsigned GetNumAtoms() const { return nlocal_; }

		//! Set the iteration
		/*!
		 * \param iteration New value for the iteration
		 */
		void SetIteration(size_t iteration) 
		{
			iteration_ = iteration; 
			changed_ = true; 
		}
		
		//! Set target iterations
		/*!
		 * \param target New value for target iterations
		 */
		void SetTargetIterations(int target) 
		{
			targetiter_ = target; 
			changed_ = true; 
		}

		//! Change the temperature
		/*!
		 * \param temperature New value for the temperature
		 */
		void SetTemperature(double temperature) 
		{ 
			temperature_ = temperature;
			changed_ = true;
		}

		//! Change the energy
		/*!
		 * \param energy New value for the energy
		 */
		void SetEnergy(double energy) 
		{
			energy_ = energy;
			changed_ = true;
		}

		//! Change the Box H-matrix. 
		/*!
		 * \param hmat New H-matrix for the system
		*/
		void SetHMatrix(const Matrix3& hmat)
		{
			H_ = hmat;
			Hinv_ = hmat.inverse();
			changed_ = true;
		}

		//! Change the box virial. 
		/*!
		 * \param virial New virial tensor for the system. 
		 */
		void SetVirial(const Matrix3& virial)
		{
			virial_ = virial; 
			changed_ = true;
		}

		//! Change the box origin.
		/*!
		 * \param origin New origin for the system
		 */
		void SetOrigin(const Vector3& origin)
		{
			origin_ = origin; 
			changed_ = true;
		}

		//! Change the periodicity of the system
		/*!
		 * \param isperiodic Periodicity of three dimensions
		 */
		void SetPeriodicity(const Bool3& isperiodic)
		{
			isperiodic_ = isperiodic;
			changed_ = true;
		}

		//! Change the kb
		/*!
		 * \param kb New value for the kb
		 */
		void SetKb(double kb) 
		{
			kb_ = kb;
			changed_ = true;
		}

		//! Set the dielectric constant.
		/*!
		 * \param dielectric Value for the dielectric constant.
		 */
		void SetDielectric(double dielectric) 
		{ 
			dielectric_ = dielectric;
			changed_ = true;
		}

		//! Set the value for qqrd2e
		/*!
		 * \param qqrd2e Value for qqrd2e.
		 */
		void Setqqrd2e(double qqrd2e) 
		{ 
			qqrd2e_ = qqrd2e;
			changed_ = true;
		}

		//! Set number of atoms in this snapshot.
		/*!
		 * \param natoms Number of atoms in this snapshot
		 */
		void SetNumAtoms(unsigned int natoms) { nlocal_ = natoms; }

		//! Access the particle positions
		/*!
		 * \return List of particle positions
		 */
		const std::vector<Vector3>& GetPositions() const { return positions_; }

		/*! \copydoc Snapshot::GetPositions() const */
		std::vector<Vector3>& GetPositions() 
		{ 
			changed_ = true;
			return positions_; 
		}

		//! Access the particle velocities
		/*!
		 * \return List of particle velocities
		 */
		const std::vector<Vector3>& GetVelocities() const { return velocities_; }

		/*! \copydoc Snapshot::GetVelocities() const */
		std::vector<Vector3>& GetVelocities() 
		{
			changed_ = true;
			return velocities_; 
		}

		//! Access the per-particle forces
		/*!
		 * \return List of per-particle forces
		 */
		const std::vector<Vector3>& GetForces() const { return forces_; }

		/*! \copydoc Snapshot::GetForces() const */
		std::vector<Vector3>& GetForces() 
		{
			changed_ = true; 
			return forces_; 
		}

		//! Const access to the particle masses
		/*!
		 * \return List of Masses
		 *
		 * Note that the Masses can be either stored as per-atom or per-type
		 * depending on the Lammps Atom type used.
		 */
		const std::vector<double>& GetMasses() const { return masses_; }

		/*! \copydoc Snapshot::GetMasses() const */
		std::vector<double>& GetMasses() 
		{
			changed_ = true; 
			return masses_; 
		}

		//! Scale a vector into fractional coordinates
		/*! 
		 * \param v Vector of interest
		 * \return Scaled vector in fractional coordinates
		 */
		Vector3 ScaleVector(const Vector3& v) const
		{
			return Hinv_*(v-origin_);
		}

		//! Apply minimum image to a vector.
		/*!
		 * \param v Vector of interest
		 */
		void ApplyMinimumImage(Vector3* v) const
		{
			Vector3 scaled = Hinv_*(*v);

			for(int i = 0; i < 3; ++i)
				scaled[i] -= isperiodic_[i]*roundf(scaled[i]);

			*v = H_*scaled;
		}
		
		//! Apply minimum image to a vector.
		/*!
		 * \param v Vector of interest
		 * \return Return new vector.
		 */
		Vector3 ApplyMinimumImage(const Vector3& v) const
		{
			Vector3 scaled = Hinv_*v;

			for(int i = 0; i < 3; ++i)
				scaled[i] -= isperiodic_[i]*roundf(scaled[i]);

			return H_*scaled;	
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
				mtot += masses_[i];
			MPI_Allreduce(MPI_IN_PLACE, &mtot, 1, MPI_DOUBLE, MPI_SUM, comm_);
			return mtot;
		}

		//! Compute center of mass of a group of atoms with implicit total mass.
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

		//! Compute center of mass of a group of atoms based on index with
		//! provided total mass.
		/*!
		 * \param indices IDs of particles of interest.
		 * \param mtot Total mass of particle group.
		 * \return Vector3 Center of mass of particles.
		 * \note Each processor passes in the local indices of the atoms of
		 *       interest and this function will collect the data and compute
		 *       the center of mass.
		 * \note If mtot is zero, then the masses are not taken into account and
		 *       the center of geometry is calculated.
		 */
		Vector3 CenterOfMass(const Label& indices, double mtot) const
		{
			// Store coorinates and masses in vectors to gather. 
			std::vector<double> pos, mass, gpos, gmass;
			std::vector<int> pcounts(comm_.size(), 0), mcounts(comm_.size(), 0); 
			std::vector<int> pdispls(comm_.size()+1, 0), mdispls(comm_.size()+1, 0);

			pcounts[comm_.rank()] = 3*indices.size();
			mcounts[comm_.rank()] = indices.size();

			// Reduce counts.
			MPI_Allreduce(MPI_IN_PLACE, pcounts.data(), pcounts.size(), MPI_INT, MPI_SUM, comm_);
			MPI_Allreduce(MPI_IN_PLACE, mcounts.data(), mcounts.size(), MPI_INT, MPI_SUM, comm_);

			// Compute displacements.
			std::partial_sum(pcounts.begin(), pcounts.end(), pdispls.begin() + 1);
			std::partial_sum(mcounts.begin(), mcounts.end(), mdispls.begin() + 1);
			
			// Fill up mass and position vectors.
			for(auto& idx : indices)
			{
				auto& p = positions_[idx];
				pos.push_back(p[0]);
				pos.push_back(p[1]);
				pos.push_back(p[2]);
				mass.push_back(masses_[idx]);
			}

			// Re-size receiving vectors. 
			gpos.resize(pdispls.back(), 0);
			gmass.resize(mdispls.back(), 0);

			// All-gather data.
			MPI_Allgatherv(pos.data(), pos.size(), MPI_DOUBLE, gpos.data(), pcounts.data(), pdispls.data(), MPI_DOUBLE, comm_);
			MPI_Allgatherv(mass.data(), mass.size(), MPI_DOUBLE, gmass.data(), mcounts.data(), mdispls.data(), MPI_DOUBLE, comm_);

			// Fictitious masses for center of geometry.
			if(mtot == 0)
			{
				for(auto& m : gmass)
				{
					m = 1.0;
					mtot++;
				}
			}

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
		/*!
		 * \return List of atom IDs
		 */
		const Label& GetAtomIDs() const { return atomids_; }

		/*! \copydoc Snapshot::GetAtomIDs() const */
		Label& GetAtomIDs()
		{
			changed_ = true;
			return atomids_;
		}

		//! Gets the local atom index corresponding to an atom ID.
		/*!
		 * \param id Atom ID.
		 * \return Local atom index or -1 if not found.
		 */
		int GetLocalIndex(int id) const
		{

			auto s = std::find(atomids_.begin(), atomids_.end(), id);
			if(s == atomids_.end())
				return -1;
			else
				return s - atomids_.begin();
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
		const std::vector<double>& GetCharges() const { return charges_; }

		/*! \copydoc Snapshot::GetCharges() const */
		std::vector<double>& GetCharges()
		{
			changed_ = true;
			return charges_;
		}
		
		//! Access the atom types
		/*!
		 * \return List of atom types
		 */
		const Label& GetAtomTypes() const { return types_; }

		/*! \copydoc Snapshot::GetAtomTypes() const */
		Label& GetAtomTypes() 
		{
			changed_ = true; 
			return types_; 
		}

		//! Access the snapshot ID
		/*!
		 * \return Snapshot ID
		 */
		const std::string& GetSnapshotID() const { return ID_; }

		/*! \copydoc Snapshot::GetSnapshotID() const */
		std::string& GetSnapshotID() 
		{
			changed_ = true; 
			return ID_; 
		}

		//! Query if Snapshot was modified
		/*!
		 * \return \c True if Snapshot was modified, else return \c False
		 */
		bool HasChanged() const { return changed_; }

		//! Set the "changed" flag of the Snapshot
		/*!
		 * \param state State to which the "changed" flag is set
		 */
		void Changed(bool state) { changed_ = state; }

		//! Return the serialized positions across all local cores
		/*!
		 * \return Positions across all local cores.
		 */
		std::vector<double> SerializePositions()
		{

			std::vector<int> pcounts(comm_.size(), 0); 
			std::vector<int> pdispls(comm_.size()+1, 0);

			pcounts[comm_.rank()] = 3*nlocal_;

			// Reduce counts.
			MPI_Allreduce(MPI_IN_PLACE, pcounts.data(), pcounts.size(), MPI_INT, MPI_SUM, comm_);

			// Compute displacements.
			std::partial_sum(pcounts.begin(), pcounts.end(), pdispls.begin() + 1);

			// Re-size receiving vectors.
			std::vector<double> positions;
			positions.resize(pdispls.back(), 0);

			std::vector<double> ptemp;
			
			for(auto& p : positions_)
			{
				ptemp.push_back(p[0]);
				ptemp.push_back(p[1]);
				ptemp.push_back(p[2]);
			}

			// All-gather data.
			MPI_Allgatherv(ptemp.data(), ptemp.size(), MPI_DOUBLE, positions.data(), pcounts.data(), pdispls.data(), MPI_DOUBLE, comm_);
			return positions;
		}

		//! Return the serialized velocities across all local cores
		/*!
		 * \return Velocities across all local cores.
		 */
		std::vector<double> SerializeVelocities()
		{
			std::vector<int> vcounts(comm_.size(), 0); 
			std::vector<int> vdispls(comm_.size()+1, 0);

			vcounts[comm_.rank()] = 3*nlocal_;

			// Reduce counts.
			MPI_Allreduce(MPI_IN_PLACE, vcounts.data(), vcounts.size(), MPI_INT, MPI_SUM, comm_);

			// Compute displacements.
			std::partial_sum(vcounts.begin(), vcounts.end(), vdispls.begin() + 1);

			// Re-size receiving vectors.
			std::vector<double> velocities;
			velocities.resize(vdispls.back(), 0);

			std::vector<double> vtemp;
			
			for(auto& v : velocities_)
			{
				vtemp.push_back(v[0]);
				vtemp.push_back(v[1]);
				vtemp.push_back(v[2]);
			}

			// All-gather data.
			MPI_Allgatherv(vtemp.data(), vtemp.size(), MPI_DOUBLE, velocities.data(), vcounts.data(), vdispls.data(), MPI_DOUBLE, comm_);
			return velocities;
		}

		//! Return the serialized IDs across all local cores
		/*!
		 * \return IDs across all local cores.
		 */
		std::vector<int> SerializeIDs()
		{
			std::vector<int> mcounts(comm_.size(), 0); 
			std::vector<int> mdispls(comm_.size()+1, 0);

			mcounts[comm_.rank()] = nlocal_;

			// Reduce counts.
			MPI_Allreduce(MPI_IN_PLACE, mcounts.data(), mcounts.size(), MPI_INT, MPI_SUM, comm_);

			// Compute displacements.
			std::partial_sum(mcounts.begin(), mcounts.end(), mdispls.begin() + 1);

			// Re-size receiving vectors.
			std::vector<int> IDs;
			IDs.resize(mdispls.back(), 0);

			// All-gather data.
			MPI_Allgatherv(atomids_.data(), atomids_.size(), MPI_INT, IDs.data(), mcounts.data(), mdispls.data(), MPI_INT, comm_);
			return IDs;
		}

		//! Destructor
		~Snapshot(){}
	};
}

