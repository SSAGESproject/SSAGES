#include "QboxHook.h"
#include "Snapshot.h"
#include "pugixml.hpp"
#include <exception>
#include <fstream>

namespace SSAGES
{
	template<typename T>
	void to_vector3(const std::string& s, T& vec)
	{
		std::istringstream ss(s);
		ss >> vec[0] >> vec[1] >> vec[2];
	}

	void QboxHook::InitializeCommands(const std::string& commandfile)
	{
		std::ofstream fout;
		fout.open(commandfile, std::ofstream::trunc);
		
		auto atomset = xml_.child("fpmd:simulation").child("iteration").child("atomset");
		for(auto& atom : atomset.children("atom"))
		{
			auto name = atom.attribute("name").value();
			fout << "extforce define atomic " 
					<< name << " " << name << " 0 0 0 " << std::endl;
		}

		fout.close();		
	}

	void QboxHook::XMLToSSAGES(const std::string& xmlfile)
	{
		// Load the XML file into property tree. 
		auto result = xml_.load_file(xmlfile.c_str()); 
		if(!result)
			throw std::runtime_error(result.description());
		
		SyncToSnapshot();
	}

	void QboxHook::SSAGESToCommands(const std::string& commandfile)
	{
		SyncToEngine();

		std::ofstream fout; 
		fout.open(commandfile, std::ofstream::trunc);
		fout.precision(10);

		const auto& pos = snapshot_->GetPositions();
		const auto& vel = snapshot_->GetVelocities();
		const auto& frc = snapshot_->GetForces();
		const auto& typ = snapshot_->GetAtomTypes();
		
		auto atomset = xml_.child("fpmd:simulation").child("iteration").child("atomset");

		int i = 0;
		for(auto& atom : atomset.children("atom"))
		{
			auto name = atom.attribute("name").value();
			fout << "extforce set " << name << " "
					<< frc[i][0] - prevforces_[i][0] << " "
					<< frc[i][1] - prevforces_[i][1] << " "
					<< frc[i][2] - prevforces_[i][2] << std::endl;

			++i;
		}
	}

	void QboxHook::BuildSpeciesInfo()
	{
		species_.clear();
		speciesmass_.clear();
		species_.reserve(32);
		speciesmass_.reserve(32);
		
		auto species = xml_.child("fpmd:simulation");
		for(auto& s : species.children("species"))
		{
			species_.push_back(s.attribute("name").value());
			speciesmass_.push_back(std::atof(s.child_value("mass")));
		}
	}

	void QboxHook::SyncToSnapshot()
	{
		auto atomset = xml_.child("fpmd:simulation").child("iteration").child("atomset");

		// First entry in <atomset> is <unit_cell>.
		int natoms = std::distance(atomset.children("atom").begin(), atomset.children("atom").end());

		// Build species index on first iteration. 
		if(snapshot_->GetIteration() == 0)
			BuildSpeciesInfo();
		
		// Load em up. 
		snapshot_->SetIteration(snapshot_->GetIteration() + 1);
		snapshot_->SetNumAtoms(natoms);

		//Set temperature.
		//Set potential energy.
		//Set kinetic energy.
		//Set boltzmann. 

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

		// Resize previous forces. 
		prevforces_.resize(natoms, Vector3{0,0,0});

		int i = 0;
		for(auto& atom : atomset.children("atom"))
		{
			to_vector3<Vector3>(atom.child_value("position"), pos[i]);
			to_vector3<Vector3>(atom.child_value("velocity"), vel[i]);
			to_vector3<Vector3>(atom.child_value("force"), frc[i]);
			prevforces_[i] = frc[i]; // Store "previous forces".
			ids[i] = i + 1;
			
			// Get species types by resolving map. 
			int type = GetSpeciesIndex(atom.attribute("species").value());
			if(type == -1)
				throw std::logic_error("Invalid atom type resolved at index " + std::to_string(i));

			typs[i] = type;
			mass[i] = speciesmass_[type];
			++i;
		}
		
		Matrix3 H; 
		Vector3 Hi; 
		auto box = atomset.child("unit_cell");
		to_vector3<Vector3>(box.attribute("a").value(), Hi);
		H.row(0) = Hi;
		to_vector3<Vector3>(box.attribute("b").value(), Hi);
		H.row(1) = Hi;
		to_vector3<Vector3>(box.attribute("c").value(), Hi);
		H.row(2) = Hi;
		snapshot_->SetHMatrix(H);
	}

	void QboxHook::SyncToEngine()
	{
	}
}