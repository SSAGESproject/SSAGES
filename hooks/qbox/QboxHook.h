#pragma once 

#include "Hook.h"
#include "pugixml.hpp"

namespace SSAGES
{
	class QboxHook : public Hook
	{
	private: 
		//! XML document
		pugi::xml_document xml_;

		//! Vector of unique species. 
		std::vector<std::string> species_; 

		//! Vector of species' masses.
		std::vector<double> speciesmass_;

		//! Previous forces for velocity verlet integration.
		std::vector<Vector3> prevforces_;

		//! Previous positions.
		std::vector<Vector3> prevpositions_;
		
		//! Boltzmann constant in Hartree/K.
		double kb = 1.0 / ( 11605.0 * 2.0 * 13.6058 );

		//! Conversion between a.m.u. and a.u.
		//! https://physics.nist.gov/cgi-bin/cuu/Value?mpsme
		double unitconv = 1836.152673;

		//! Get species index from species vector.
		int GetSpeciesIndex(const std::string& species)
		{
			auto it = std::find(species_.begin(), species_.end(), species);
			if(it == species_.end())
				return -1; 

			return it - species_.begin();
		}

		void BuildSpeciesInfo();

	protected:
		// Implementation of the SyncToEngine interface.
		void SyncToEngine() override;

		// Implementation of the SyncToSnapshot interface. 
		void SyncToSnapshot() override;

	public: 
		QboxHook() : species_(), speciesmass_(), prevforces_(0), prevpositions_()
		{}

		void XMLToSSAGES(const std::string& xmlfile);

		void SSAGESToCommands(const std::string& commandfile);
		
		void InitializeCommands(const std::string& commandfile);

		~QboxHook()
		{}
	};
}
