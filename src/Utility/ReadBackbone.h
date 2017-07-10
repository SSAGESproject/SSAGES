/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2017 Ashley Guo <azguo@uchicago.edu>
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

#include <fstream>

// Not a full PDB parser!!! Only grabs backbone atoms given residue numbers and
// a reference PDB. This will not parse any other PDB record type other than 
// "ATOM", and does not read any info from a PDB file past the column containing
// residue number. Could be extended to parse a PDB file completely, which would
// be better placed in ReadFile.h, but sticking with minimum requirements for 
// secondary structure CVs here.

namespace SSAGES
{
	//! Utility class to read protein backbone atoms from a reference file
	/*!
	 * This utility class allows for the backbone atoms along a protein 
	 * segment to be extracted from a reference file (most likely in a .pdb).
	 *
	 * Supported filetypes:
	 * * PDB
	 */
	class ReadBackbone
	{
	public:
		//! Read protein backbone atoms from reference PDB file
		/*!
		 * \param refpdb Name of reference PDB file.
		 * \param resids Vector of residue numbers for which to find backbone atoms.
		 * \return atomids Vector of atom numbers corresponding to specified residues. 
		 *
		 * Extract backbone atom from a reference PDB file, corresponding to a
		 * sequence of protein residues of interest. Five backbone atoms are 
		 * extracted for each residue in the order: N CA CB C O. For Glycine 
		 * residues, HA1 is located in place of CB. 
		 */
		//static std::vector<int> GetPdbBackbone(std::string refpdb, std::vector<int> resids)
		static std::vector< std::vector< std::string> > GetPdbBackbone(std::string refpdb, std::vector<int> resids)
		{
			//std::vector<int> atomids;
			std::vector< std::vector< std::string > > atomids(2, std::vector<std::string> (0));
			//std::vector<int> atomnums;
			std::vector<std::string> atomnums;
			std::vector<std::string> atomtypes;
			std::vector<std::string> resnames;
			std::vector<std::string> chains;
			std::vector<int> resnums;
			//std::vector<int> tempatoms(5);
			std::vector<std::string> tempatoms(5);
			std::vector<std::string> backboneAtoms = {"N", "CA", "CB", "C", "O", "OT1"};
			std::vector<std::string> glycineAtoms = {"N", "CA", "HA1", "C", "O", "OT1"};
			std::string atomnum, resnum;
			std::string record, atomtype, resname, chain;
			std::ifstream pdbfile;
			pdbfile.open(refpdb, std::ios::in);
			std::string line;

			while(std::getline(pdbfile, line)){
				if(line.length() < 26) line.append(26 - line.length(), ' ');
				record = line.substr(0, 6);
				atomnum = line.substr(6, 5);
				atomtype = line.substr(12, 4);
				resname = line.substr(17, 3);
				chain = line.substr(21, 1);
				resnum = line.substr(22, 4);
				record.erase( std::remove( record.begin(), record.end(), ' '), record.end());
				if((record == "ATOM") && ((std::stoul(resnum) - resids[0]) <= (resids.back() - resids[0]))){
					atomtype.erase( std::remove( atomtype.begin(), atomtype.end(), ' '), atomtype.end());
					resname.erase( std::remove( resname.begin(), resname.end(), ' '), resname.end());
					if(resname == "GLY" && std::find(glycineAtoms.begin(), glycineAtoms.end(), atomtype) != glycineAtoms.end()){
						//atomnums.push_back(std::stoul(atomnum));
						atomnums.push_back(atomnum);
						atomtypes.push_back(atomtype);
						resnames.push_back(resname);
						chains.push_back(chain);
						resnums.push_back(std::stoul(resnum));
					} else if(std::find(backboneAtoms.begin(), backboneAtoms.end(), atomtype) != backboneAtoms.end()){
						//atomnums.push_back(std::stoul(atomnum));
						atomnums.push_back(atomnum);
						atomtypes.push_back(atomtype);
						resnames.push_back(resname);
						chains.push_back(chain);
						resnums.push_back(std::stoul(resnum));
					}
				}
			}

			// atoms in PDB not necessarily in N CA CB C O order, fix that:
			for(size_t i = 0; i < resids.size(); i++){
				for(size_t j = 5 * i; j < (5 * i + 5); j++){
					if(atomtypes[j] == "N"){
						tempatoms[0] = atomnums[j]; 
						std::cout << "N  - Atom " << tempatoms[0] << " - Chain " << chains[j] << std::endl;
					} else if(atomtypes[j] == "CA"){
						tempatoms[1] = atomnums[j]; 
						std::cout << "CA - Atom " << tempatoms[1] << " - Chain " << chains[j] << std::endl;
					} else if(atomtypes[j] == "C"){
						tempatoms[3] = atomnums[j]; 
						std::cout << "C  - Atom " << tempatoms[3] << " - Chain " << chains[j] << std::endl;
					} else if(atomtypes[j] == "O"){
						tempatoms[4] = atomnums[j]; 
						std::cout << "O  - Atom " << tempatoms[4] << " - Chain " << chains[j] << std::endl;
					} else if(atomtypes[j] == "OT1"){
						tempatoms[4] = atomnums[j];
						std::cout << "OT1- Atom " << tempatoms[4] << " - Chain " << chains[j] << std::endl;
					} else{
						tempatoms[2] = atomnums[j]; 
						std::cout << atomtypes[j] << " - Atom " << tempatoms[2] << " - Chain " << chains[j] << std::endl;
					}
					atomids[1].push_back(chains[j]);
				}
				//atomids.insert(atomids.end(), tempatoms.begin(), tempatoms.end());
				atomids[0].insert(atomids[0].end(), tempatoms.begin(), tempatoms.end());
			}

			return atomids;
		}
	};
}
