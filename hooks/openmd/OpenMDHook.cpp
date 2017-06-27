#include "OpenMDHook.h"
#include "Snapshot.h"

#include "brains/SnapshotManager.hpp"
#include "brains/Thermo.hpp"
#include "utils/PhysicalConstants.hpp"
#include "primitives/Molecule.hpp"

using namespace OpenMD;

namespace SSAGES
{
	void OpenMDHook::SyncToEngine() 
	{
		const auto& pos = snapshot_->GetPositions();
		const auto& vel = snapshot_->GetVelocities();
		const auto& frc = snapshot_->GetForces();
		const auto& mass = snapshot_->GetMasses();
		const auto& ids = snapshot_->GetAtomIDs();
		
		// Loop through and set atom properties.
		size_t index = 0;
		SimInfo::MoleculeIterator i;
    	Molecule::AtomIterator  j;
		for (auto* mol = siminfo_->beginMolecule(i); mol != NULL; mol = siminfo_->nextMolecule(i)) 
		{
      		for (auto* atom = mol->beginAtom(j); atom != NULL; atom = mol->nextAtom(j)) 
      		{
      			atom->setMass(mass[index]);
      			atom->setPos({pos[index][0], pos[index][1], pos[index][2]});
      			atom->setVel({vel[index][0], vel[index][1], vel[index][2]});
      			atom->setFrc({frc[index][0], frc[index][1], frc[index][2]});
      			++index;
      		}
      	}
	}

	void OpenMDHook::SyncToSnapshot()
	{
		auto natoms = siminfo_->getNGlobalAtoms();

		// Resize vectors.
		auto& pos = snapshot_->GetPositions();
		pos.resize(natoms);
		auto& vel = snapshot_->GetVelocities();
		vel.resize(natoms);
		auto& frc = snapshot_->GetForces();
		frc.resize(natoms);
		auto& mass = snapshot_->GetMasses();
		mass.resize(natoms);
		auto& ids = snapshot_->GetAtomIDs();
		ids.resize(natoms);
		auto& typs = snapshot_->GetAtomTypes();
		typs.resize(natoms);

		// Get OpenMD snapshot. 
		auto* osnap = siminfo_->getSnapshotManager()->getCurrentSnapshot();

		// Get temperature/energy/pressure.
		auto& comm = snapshot_->GetCommunicator();
		auto thermo = Thermo(siminfo_);

		snapshot_->SetIteration(snapshot_->GetIteration()+1);
		snapshot_->SetTemperature(thermo.getTemperature());
		snapshot_->SetEnergy(thermo.getTotalEnergy());
		snapshot_->SetKb(PhysicalConstants::kB);

		// Loop through and get atom properties.
		size_t index = 0;
		SimInfo::MoleculeIterator i;
    	Molecule::AtomIterator  j;
		for (auto* mol = siminfo_->beginMolecule(i); mol != NULL; mol = siminfo_->nextMolecule(i)) 
		{
      		for (auto* atom = mol->beginAtom(j); atom != NULL; atom = mol->nextAtom(j)) 
      		{
        		auto v = atom->getVel();
        		auto p = atom->getPos();
        		auto f = atom->getFrc();

        		//TODO: Get integer type of stunt double.
        		ids[index] = atom->getGlobalIndex() + 1;
        		mass[index] = atom->getMass();

        		pos[index][0] = p[0];
        		pos[index][1] = p[1];
        		pos[index][2] = p[2];

        		vel[index][0] = v[0];
        		vel[index][1] = v[1];
        		vel[index][2] = v[2];

        		frc[index][0] = f[0];
        		frc[index][1] = f[1];
        		frc[index][2] = f[2];

        		++index;
        	}
        }

        auto Homd = osnap->getHmat();
        Matrix3 H;
		H << Homd(0,0), Homd(0,1), Homd(0,2),
		     Homd(1,0), Homd(1,1), Homd(1,2),
		     Homd(2,0), Homd(2,1), Homd(2,2);
		snapshot_->SetHMatrix(H);

        pos.resize(index);
        vel.resize(index);
        frc.resize(index);
        mass.resize(index);
        ids.resize(index);
        typs.resize(index);
        snapshot_->SetNumAtoms(index);

		Hook::PostStepHook();
	}
}