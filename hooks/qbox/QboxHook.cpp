#include "QboxHook.h"
#include <boost/property_tree/xml_parser.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <exception>

namespace SSAGES
{
	template<typename T>
	std::vector<T> to_vector(std::string s)
	{
		using boost::algorithm::trim;
		trim(s);

		std::vector<T> result;
		std::stringstream ss(s);
		std::string item;
		while(std::getline(ss, item, ' ')) 
			result.push_back(boost::lexical_cast<T>(item));
		
		return result;
	}

	template<typename T>
	void to_vector3(const std::string& s, T& vec)
	{
		std::istringstream ss(s);
		ss >> vec[0] >> vec[1] >> vec[2];
	}

	void QboxHook::XMLToSSAGES(const std::string& xmlfile)
	{
		// Load the XML file into property tree. 
		read_xml(xmlfile, pt_);
		SyncToEngine();
	}

	void QboxHook::BuildSpeciesInfo()
	{
		species_.clear();
		speciesmass_.clear();
		species_.reserve(32);
		speciesmass_.reserve(32);

		auto& species = pt_.get_child("fpmd:simulation");
		for(auto& s : species)
		{
			if(s.first == "species")
			{
				species_.push_back(s.second.get<std::string>("<xmlattr>.name"));
				speciesmass_.push_back(s.second.get<double>("mass"));
			}
		}
	}

	void QboxHook::SyncToEngine()
	{

		// Build species list and masses.
		BuildSpeciesInfo();

		auto& atomset = pt_.get_child("fpmd:simulation.iteration.atomset");

		// First entry in <atomset> is <unit_cell>.
		int natoms = atomset.size() - 1;
		
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

		int i = 0;
		for(const auto& v : atomset)
		{
			if(v.first == "atom")
			{
				to_vector3<Vector3>(v.second.get<std::string>("position"), pos[i]);
				to_vector3<Vector3>(v.second.get<std::string>("velocity"), vel[i]);
				to_vector3<Vector3>(v.second.get<std::string>("force"), frc[i]);
				ids[i] = i + 1;

				// Get species types by resolving map. 
				int type = GetSpeciesIndex(v.second.get<std::string>("<xmlattr>.species"));
				if(type == -1)
					throw std::logic_error("Invalid atom type resolved at index " + std::to_string(i));

				typs[i] = type;
				mass[i] = speciesmass_[type];

				++i;
			} 
			else if(v.first == "unit_cell")
			{
				Matrix3 H; 
				Vector3 Hi; 
				to_vector3<Vector3>(v.second.get<std::string>("<xmlattr>.a"), Hi);
				H.row(0) = Hi;
				to_vector3<Vector3>(v.second.get<std::string>("<xmlattr>.b"), Hi);
				H.row(1) = Hi;
				to_vector3<Vector3>(v.second.get<std::string>("<xmlattr>.c"), Hi);
				H.row(2) = Hi;
				_snapshot->SetHMatrix(H);
			}
		}
		
	}

	void QboxHook::SyncToSnapshot()
	{}
}