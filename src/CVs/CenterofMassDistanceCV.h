#pragma once 

#include "CollectiveVariable.h"
#include <array>
#include <math.h>

namespace SSAGES
{
	//! Collective variable on center of mass for group of atoms.
	/*!
	 * This will return the center of mass for a group of atoms       s
	 *
	 * \ingroup CVs
	 */
	class CenterofMassDistanceCV : public CollectiveVariable
	{
	private:

		//! IDs of atoms of interest
		std::vector<int> _atomids1;
		std::vector<int> _atomids2;

		//! Array to store indexes
		std::vector<double> _pertatoms1;
		std::vector<double> _pertatoms2;

		//! Current value of the CV.
		double _val;

		//! Gradient of the CV, dCOM/dxi.
		std::vector<Vector3> _grad;

		//! Bounds on CV.
		std::array<double, 2> _bounds;

	public:
		//! Constructor.
		/*!
		 * \param atomids IDs of the atoms defining the COM.
		 * \param atom_range If \c True Use a range of atoms for the COM based on two atoms supplied in atomids.
		 *
		 * Construct COM.
		 *
		 * \todo Bounds needs to be an input.
		 */
		CenterofMassDistanceCV(std::vector<int> atomids1, std::vector<int> atomids2, bool use_range1 = false, bool use_range2 = false) ://std::array<bool,2> use_range = {false,false}) :
		_atomids1(atomids1), _atomids2(atomids2), _val(0), _grad(0), _bounds{{0,0}}
		{
			if(use_range1) //[0])
			{
				if(atomids1.size() != 2)
					throw std::runtime_error("COMCV input 1: If using range, must define only two atoms!");

				_atomids1.clear();

				if(atomids1[0] >= atomids1[1])
					throw std::runtime_error("COMCV input 1: Please reverse atom range or check that atom range is not equal!");

				for(int i = atomids1[0]; i <= atomids1[1]; i++)
					_atomids1.push_back(i);
			}

			if(use_range2)//[1])
			{
				if(atomids2.size() != 2)
					throw std::runtime_error("COMCV input 2: If using range, must define only two atoms!");

				_atomids2.clear();

				if(atomids2[0] >= atomids2[1])
					throw std::runtime_error("COMCV input 2: Please reverse atom range or check that atom range is not equal!");

				for(int i = atomids2[0]; i <= atomids2[1]; i++)
					_atomids2.push_back(i);
			}
		}

		// Initialize variables
		void Initialize(const Snapshot& snapshot) override
		{
			// Initialize gradient
			auto n = snapshot.GetPositions().size();
			_grad.resize(n);
		}

		// Evaluate the CV
		void Evaluate(const Snapshot& snapshot) override
		{
			const auto& pos = snapshot.GetPositions();
			const auto& ids = snapshot.GetAtomIDs();
			const auto& mass = snapshot.GetMasses();
			const auto& image_flags = snapshot.GetImageFlags();

			Vector3 mass_pos_prod1 = {0,0,0}; //= std::vector<double>(3,0.0);
			double total_mass1 = 0;
			Vector3 mass_pos_prod2 = {0,0,0};; //= std::vector<double>(3,0.0);
			double total_mass2 = 0;

			// Loop through atom ids
			for( size_t i = 0; i < ids.size(); ++i)
			{
				_grad[i].setZero();

				// Pull out atoms of interst
				for(size_t j = 0; j<_atomids1.size(); j++)
				{
					if(_atomids1[j] == ids[i])
					{
						_pertatoms1.push_back(i);
						auto u_coord = snapshot.UnwrapVector(pos[i], image_flags[i]);

						mass_pos_prod1[0] += mass[i]*u_coord[0];
						mass_pos_prod1[1] += mass[i]*u_coord[1];
						mass_pos_prod1[2] += mass[i]*u_coord[2];
						total_mass1 += mass[i];
						break;
					}
				}

				// Pull out atoms of interst
				for(size_t j = 0; j<_atomids2.size(); j++)
				{
					if(_atomids2[j] == ids[i])
					{
						_pertatoms2.push_back(i);
						auto u_coord = snapshot.UnwrapVector(pos[i], image_flags[i]);
						
						mass_pos_prod2[0] += mass[i]*u_coord[0];
						mass_pos_prod2[1] += mass[i]*u_coord[1];
						mass_pos_prod2[2] += mass[i]*u_coord[2];

						total_mass2 += mass[i];
						break;
					}
				}
			}

			Vector3 COM1; //= mass_pos_prod1/total_mass1;
			Vector3 COM2;//= mass_pos_prod2/total_mass2;

			COM1 = mass_pos_prod1/total_mass1;
			COM2 = mass_pos_prod2/total_mass2;

			auto del = snapshot.ApplyMinimumImage(COM1 - COM2);
			auto r = del.norm();
			_val = r ; 
			int i;
			// Calculate gradient
			for(size_t j = 0; j < _pertatoms1.size(); ++j)
			{
				i = _pertatoms1[j];
				for(size_t k = 0; k<del.size(); ++k)
				{
					_grad[i][k] += del[k]/r * mass[i]/total_mass1;

				}

			}

			for(size_t j = 0; j < _pertatoms2.size(); ++j)
			{
				i = _pertatoms2[j];
				for(size_t k = 0; k<del.size(); ++k)
				{
					_grad[i][k] += -del[k]/r * mass[i]/total_mass2;

				}

			}
		}
			// Return the value of the CV.
		double GetValue() const override 
		{ 
			return _val; 
		}

		double GetPeriodicValue(double Location) const override
		{
			return Location;
		}

		// Return the gradient of the CV.
		const std::vector<Vector3>& GetGradient() const override
		{
			return _grad;
		}

		// Return the boundaries of the CV.
		const std::array<double, 2>& GetBoundaries() const override
		{
			return _bounds;
		}

		double GetDifference(const double Location) const override
		{
			return _val - Location;
		}
		//! Serialize this CV for restart purposes.
		/*!
		 * \param json JSON value
		 */
		virtual void Serialize(Json::Value& json) const override
		{
			json["type"] = "CenterofMassDistance";
			for(size_t i=0; i < _atomids1.size(); ++i)
				json["atom ids1"].append(_atomids1[i]);
			json["type"] = "CenterofMassDistance";
			for(size_t i=0; i < _atomids2.size(); ++i)
				json["atom ids1"].append(_atomids2[i]);
			for(size_t i = 0; i < _bounds.size(); ++i)
				json["bounds"].append(_bounds[i]);
		}
	};
}
