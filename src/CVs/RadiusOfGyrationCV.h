#pragma once 

#include "CollectiveVariable.h"

#include <array>
#include <math.h>

namespace SSAGES
{
	//! Collective variable on the radius of gyration.
	/*!
	 * Collective variable on radius of gyration. This will return the
	 * Rg for a group of atoms.
	 *
	 * \ingroup CVs
	 */
	class RadiusOfGyrationCV : public CollectiveVariable
	{
	private:
		
		std::vector<int> _atomids; //!< IDs of the atoms used for Rg calculation

		std::vector<int> _pertatoms; //!< Array to store indicies of atoms of interest

		//! Current value of the CV.
		double _val;
		
		//! Gradient of the CV, dRg/dxi.
		std::vector<Vector3> _grad;

		//! Bounds on CV.
		std::array<double, 2> _bounds;

	public:
		//! Constructor.
		/*!
		 * \param atomids IDs of the atoms defining Rg.
		 * \param userange If \c True Use range of atoms defined by the two atoms in atomids.
		 *
		 * Construct a radius of gyration CV.
		 *
		 * \todo Bounds needs to be an input and periodic boundary conditions
		 */
		RadiusOfGyrationCV(std::vector<int> atomids, bool use_range = false) :
		_atomids(atomids), _val(0), _grad(0), _bounds{{0,0}}
		{
			if(use_range)
			{
				if(atomids.size() != 2)
				{
					std::cout<<"RgCV: If using range, must define only two atoms!"<<std::endl;
					exit(0);
				}

				_atomids.clear();

				if(atomids[0] >= atomids[1])
				{
					std::cout<<"RgCV: Please reverse atom range or check that atom range is not equal!"<<std::endl;
					exit(0);
				}

				for(int i = atomids[0]; i <= atomids[1];i++)
					_atomids.push_back(i);
			}
		}

		//! Initialize necessary variables.
		/*!
		 * \param snapshot Current simulation snapshot.
		 */
		void Initialize(const Snapshot& snapshot) override
		{
			// Initialize gradient
			auto n = snapshot.GetPositions().size();
			_grad.resize(n);
			auto m = _atomids.size();
			_pertatoms.resize(m);
		}

		//! Evaluate the CV.
		/*!
		 * \param snapshot Current simulation snapshot.
		 */
		void Evaluate(const Snapshot& snapshot) override
		{
			const auto& pos = snapshot.GetPositions();
			const auto& ids = snapshot.GetAtomIDs();
			const auto& mass = snapshot.GetMasses();
			const auto& image_flags = snapshot.GetImageFlags();

			// Calculate center of mass
			Vector3 mass_pos_prod = {0,0,0}; //= std::vector<double>(3,0.0);
			double total_mass = 0;

			// Loop through atom positions
			for( size_t i = 0; i < pos.size(); ++i)
			{
				_grad[i].setZero();
				// Loop through pertinent atoms
				for(size_t j = 0; j < _atomids.size(); j++)
				{
					if(ids[i] == _atomids[j])
					{
						_pertatoms[j] = i;
						auto u_coord = pos[i]; //UnwrapCoordinates(snapshot.GetLatticeConstants(), image_flags[i], pos[i]);

						mass_pos_prod[0] += mass[i]*u_coord[0];
						mass_pos_prod[1] += mass[i]*u_coord[1];
						mass_pos_prod[2] += mass[i]*u_coord[2];
						total_mass += mass[i];
						break;
					}
				}
			}

			Vector3 COM;
			for(size_t i = 0; i<COM.size(); i++)
				COM[i] = mass_pos_prod[i]/total_mass;

			// Compute Radius of Gyration
			double num_gyr=0;
			int i;

			for(size_t j = 0; j < _pertatoms.size(); ++j)
			{
				i = _pertatoms[j];
				auto u_coord = pos[i]; //UnwrapCoordinates(snapshot.GetLatticeConstants(), image_flags[i], pos[i]);

				auto diff = u_coord - COM;
				auto diffnorm2 = diff.norm()*diff.norm();

				num_gyr += mass[i]*diffnorm2;

				auto massratio = mass[i]/total_mass;
				_grad[i] = massratio*diffnorm2*diff*(1.0-massratio);
			}

			//Radius of gyration
			_val = sqrt(num_gyr/total_mass);

			// Calculate gradients
			for(size_t j = 0; j < _pertatoms.size(); ++j)
			{
				i = _pertatoms[j];
				_grad[i] *= (1./_val);
			}
		}

		//! Return the value of the CV.
		/*!
		 * \return Current value of the CV.
		 */
		double GetValue() const override 
		{ 
			return _val; 
		}

		//! Return value taking periodic boundary conditions into account.
		/*!
		 * \return Input value
		 *
		 * The Rg CV does not consider periodic boundary
		 * conditions. Thus, this function always returns the input value.
		 */
		double GetPeriodicValue(double Location) const override
		{
			return Location;
		}

		//! Return the gradient of the CV.
		/*!
		 * \return Gradient of the CV.
		 */
		const std::vector<Vector3>& GetGradient() const override
		{
			return _grad;
		}

		//! Return the boundaries of the CV.
		/*!
		 * \return Values of the lower and upper boundary.
		 */
		const std::array<double, 2>& GetBoundaries() const override
		{
			return _bounds;
		}

		//! Return difference considering periodic boundaries.
		/*!
		 * \return Direct difference.
		 *
		 * As the Rg CV does not consider periodic boundary
		 * conditions, the difference between the current value of the CV and
		 * another value is always the direct difference.
		 */
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
			json["type"] = "RadiusOfGyration";
			for(size_t i=0; i < _atomids.size(); ++i)
				json["atom ids"].append(_atomids[i]);
			for(size_t i = 0; i < _bounds.size(); ++i)
				json["bounds"].append(_bounds[i]);
		}
	};
}
