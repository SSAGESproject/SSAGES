#include <iostream>
#include <numeric>
#include "HOOMDHook.h"
#include "Snapshot.h"

using namespace SSAGES;

HOOMDHook::HOOMDHook() :
Hook(), HalfStepHook()
{
  std::cerr << "Installing SSAGES HOOMD-blue hook." << std::endl;
}


void HOOMDHook::update(unsigned int timestep)
{
    timestep_ = timestep;
    SyncToSnapshot();
    Hook::PostIntegrationHook();
}

void HOOMDHook::SyncToSnapshot() //put HOOMD values -> Snapshot
{
    const auto pdata = sysdef_->getParticleData();
    unsigned int n = pdata->getN();
    snapshot_->SetNumAtoms(n);

    auto& snap_positions = snapshot_->GetPositions();
    snap_positions.resize(n);
    auto& snap_velocities = snapshot_->GetVelocities();
    snap_velocities.resize(n);
    auto& snap_forces = snapshot_->GetForces();
    snap_forces.resize(n);
    auto& snap_masses = snapshot_->GetMasses();
    snap_masses.resize(n);
    auto& snap_ids = snapshot_->GetAtomIDs();
    snap_ids.resize(n);
    auto& snap_types = snapshot_->GetAtomTypes();
    snap_types.resize(n);
    auto& snap_charges = snapshot_->GetCharges();
    snap_charges.resize(n);

    {
        ArrayHandle<Scalar4> positions(pdata->getPositions(), access_location::host, access_mode::read);
        ArrayHandle<Scalar4> velocities(pdata->getVelocities(), access_location::host, access_mode::read);
        ArrayHandle<Scalar4> forces(pdata->getNetForce(), access_location::host, access_mode::read);
        ArrayHandle<Scalar> virial(pdata->getNetVirial(), access_location::host, access_mode::read);
        ArrayHandle<Scalar> charges(pdata->getCharges(), access_location::host, access_mode::read);
        ArrayHandle<unsigned int> tags(pdata->getTags(), access_location::host, access_mode::read);

        BoxDim box(pdata->getGlobalBox());

        // Get H-matrix of the box
        Matrix3 H;
        auto L = box.getL();
        auto xy = box.getTiltFactorXY();
        auto xz = box.getTiltFactorXZ();
        auto yz = box.getTiltFactorYZ();
        H << L.x, xy*L.y, xz*L.z,
               0,    L.y, yz*L.z,
               0,      0,    L.z;
        snapshot_->SetHMatrix(H);

        // Get box origin.
        auto boxLo = box.getLo();
        Vector3 snap_origin = {
            boxLo.x,
            boxLo.y,
            boxLo.z
        };
        snapshot_->SetOrigin(snap_origin);

        // Set periodicity.
        auto boxPeriodicity = box.getPeriodic();
        Bool3 snap_periodicity = {
            boxPeriodicity.x != 0,
            boxPeriodicity.y != 0,
            boxPeriodicity.z != 0
        };
        snapshot_->SetPeriodicity(snap_periodicity);

        for (unsigned int i = 0; i < n; ++i)
        {
            // Positions
            snap_positions[i][0] = positions.data[i].x; //x
            snap_positions[i][1] = positions.data[i].y; //y
            snap_positions[i][2] = positions.data[i].z; //z
            snap_types[i] = __scalar_as_int(positions.data[i].w); //type

            // Velocities
            snap_velocities[i][0] = velocities.data[i].x; //x
            snap_velocities[i][1] = velocities.data[i].y; //y
            snap_velocities[i][2] = velocities.data[i].z; //z
            snap_masses[i] = velocities.data[i].w; //mass

            // Forces
            snap_forces[i][0] = forces.data[i].x;
            snap_forces[i][1] = forces.data[i].y;
            snap_forces[i][2] = forces.data[i].z;

            // Tags
            snap_ids[i] = tags.data[i];

            // Charges
            snap_charges[i] = charges.data[i];
        }
    }

    thermo_->compute(timestep_);
    snapshot_->SetTemperature(thermo_->getTemperature());
    auto etot = thermo_->getKineticEnergy() + thermo_->getPotentialEnergy();
    snapshot_->SetEnergy(etot/n);
    snapshot_->SetIteration(timestep_);
    snapshot_->SetKb(1); // TODO: verify that Kb = 1
    //snapshot_->SetDielectric( /* TODO: Do we need this? */ );
    //snapshot_->Setqqrd2e( /* TODO: Do we need this? */ );

    // TODO: LAMMPS fix sets this to zero, not sure if this is desired for HOOMD
    // Zero the virial - we are only interested in accumulation.
    snapshot_->SetVirial(Matrix3::Zero());

}

void HOOMDHook::SyncToEngine() //put Snapshot values -> HOOMD
{
    // Obtain local const reference to snapshot variables.
    // Const will ensure that snapshot_ variables are
    // not being changed. Only engine side variables should
    // change.
    const auto& snap_positions = snapshot_->GetPositions();
    const auto& snap_velocities = snapshot_->GetVelocities();
    const auto& snap_forces = snapshot_->GetForces();
    const auto& snap_masses = snapshot_->GetMasses();
    const auto& snap_ids = snapshot_->GetAtomIDs();
    const auto& snap_types = snapshot_->GetAtomTypes();
    const auto& snap_charges = snapshot_->GetCharges();

    auto pdata = sysdef_->getParticleData();
    unsigned int n = pdata->getN();
    {
        ArrayHandle<Scalar4> positions(pdata->getPositions(), access_location::host, access_mode::overwrite);
        ArrayHandle<Scalar4> velocities(pdata->getVelocities(), access_location::host, access_mode::overwrite);
        ArrayHandle<Scalar4> forces(pdata->getNetForce(), access_location::host, access_mode::overwrite);
        ArrayHandle<Scalar> virial(pdata->getNetVirial(), access_location::host, access_mode::overwrite);
        ArrayHandle<Scalar> charges(pdata->getCharges(), access_location::host, access_mode::overwrite);
        ArrayHandle<unsigned int> tags(pdata->getTags(), access_location::host, access_mode::overwrite);

        for (unsigned int i = 0; i < n; ++i)
        {
            // Positions
            positions.data[i].x = snap_positions[i][0]; //x
            positions.data[i].y = snap_positions[i][1]; //y
            positions.data[i].z = snap_positions[i][2]; //z
            positions.data[i].w = __int_as_scalar(snap_types[i]); //type

            // Velocities
            velocities.data[i].x = snap_velocities[i][0]; //x
            velocities.data[i].y = snap_velocities[i][1]; //y
            velocities.data[i].z = snap_velocities[i][2]; //z
            velocities.data[i].w = snap_masses[i]; //mass

            // Forces
            forces.data[i].x = snap_forces[i][0]; //x
            forces.data[i].y = snap_forces[i][1]; //y
            forces.data[i].z = snap_forces[i][2]; //z

            // Tags
            tags.data[i] = snap_ids[i];

            // Charges
            charges.data[i] = snap_charges[i];
        }

        // Update the virial (negative contribution to box).
        // TODO: This doesn't alter the virial yet.
        /*
        const auto& snap_virial = snapshot_->GetVirial();
        virial.data[0] = -snap_virial(0,0);
        virial.data[1] = -snap_virial(1,1);
        virial.data[2] = -snap_virial(2,2);
        virial.data[3] = -snap_virial(0,1);
        virial.data[4] = -snap_virial(0,2);
        virial.data[5] = -snap_virial(1,2);
        */

        // Update the simulation box
        Matrix3 H = snapshot_->GetHMatrix();
        auto Lx = H(0, 0);
        auto Ly = H(1, 1);
        auto Lz = H(2, 2);
        auto xy = H(0, 1) / Lx;
        auto xz = H(0, 2) / Lx;
        auto yz = H(1, 2) / Ly;
        BoxDim newBox(Lx, Ly, Lz);
        newBox.setTiltFactors(xy, xz, yz);
        pdata->setGlobalBox(newBox);
    }
}

HOOMDHook::~HOOMDHook() {}
