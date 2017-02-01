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
		const auto& pos = _snapshot->GetPositions();
		const auto& vel = _snapshot->GetVelocities();
		const auto& frc = _snapshot->GetForces();
		const auto& mass = _snapshot->GetMasses();
		const auto& ids = _snapshot->GetAtomIDs();
		
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
		auto& pos = _snapshot->GetPositions();
		pos.resize(natoms);
		auto& vel = _snapshot->GetVelocities();
		vel.resize(natoms);
		auto& frc = _snapshot->GetForces();
		frc.resize(natoms);
		auto& mass = _snapshot->GetMasses();
		mass.resize(natoms);
		auto& ids = _snapshot->GetAtomIDs();
		ids.resize(natoms);
		auto& typs = _snapshot->GetAtomTypes();
		typs.resize(natoms);

		// Get OpenMD snapshot. 
		auto* osnap = siminfo_->getSnapshotManager()->getCurrentSnapshot();

		// Get temperature/energy/pressure.
		auto& comm = _snapshot->GetCommunicator();
		auto thermo = Thermo(siminfo_);

		_snapshot->SetIteration(_snapshot->GetIteration()+1);
		_snapshot->SetTemperature(thermo.getTemperature());
		_snapshot->SetEnergy(thermo.getTotalEnergy());
		_snapshot->SetKb(PhysicalConstants::kB);

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
        		ids[index] = atom->getGlobalIndex();
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
		_snapshot->SetHMatrix(H);

        pos.resize(index);
        vel.resize(index);
        frc.resize(index);
        mass.resize(index);
        ids.resize(index);
        typs.resize(index);
        _snapshot->SetNumAtoms(index);

		Hook::PostStepHook();
	}
}