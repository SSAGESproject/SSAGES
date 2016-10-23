/**
 * This file is part of
 * SSAGES - Suite for Advanced Generalized Ensemble Simulations
 *
 * Copyright 2016 Jiyuan Li <jyli@uchicago.edu>
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

#include "Constraint.h"
#include <iostream>
#include <iomanip>
#include <boost/mpi.hpp>
#include <fstream>
#include <vector>
#include <math.h>

namespace SSAGES
{
	//! Image method
	/*!
	 * Image method to include polarization corrections into electrostatic
	 * interactions for systems where dielectric objects are embedded in
	 * dielectric continuum.
	 */
	class COPSSImage : public Constraint
	{
	private:
		//! Dielectric constant of polarizable particles.
		double _einner;

		//! Where non polarizable particles start.
		int _ion_type_start;
		
		//! Dielectric constant of outside continuum.
		double _eouter;

		//! unit conversion constant.
		double _qqrd2e;
		
		//! Lower value for x.
		double xlo = 0.0;
		
		//! Upper value for x.
		double xhi = 1.0;
	
		//! Number of gaussian integrations.
		int ngauss = 5;

		//! Magic numbers for x.
		double _xg0[5] = {-0.9061798459386640,-0.5384693101056831,0.00000000000000000,0.5384693101056831,0.9061798459386640};

		//! Magic numbers for the weight.
   		double _wg0[5] = {0.2369268850561891,0.4786286704993665,0.5688888888888889,0.4786286704993665,0.2369268850561891};	
		
		//! List of radii for all atom types.
		std::vector<double> _atomTypeRadius;
	
	public:
		//! Constructor.
		/*!
		 * \param comm MPI global communicator.
		 * \param frequency Frequency with which this method is invoked.
		 * \param einner Inner bound for electrostatics.
		 * \param ion_type_start The ion type to start with.
		 * \param atomTypeRadius List of atom radii.
		 */
		COPSSImage(boost::mpi::communicator& comm, 
			   unsigned int frequency, 
			   double einner, 
			   int ion_type_start,
			   std::vector<double> atomTypeRadius):
		Constraint(frequency, comm), _einner(einner), _ion_type_start(ion_type_start), _eouter(0), _qqrd2e(1), _atomTypeRadius(atomTypeRadius){}
		
		//! Gauss integration auxiliary parameter e.
		double _e;

		//! Gauss integration auxiliary parameter g inverse.
		double _ginv;

		//! Auxiliary variable R_kj for image kernel functions (x component).
		double Rxkj;

		//! Auxiliary variable R_kj for image kernel functions (y component).
		double Rykj;

		//! Auxiliary variable R_kj for image kernel functions (z component).
		double Rzkj;

		//! Auxiliary variable R_kj squared.
		double Rkj2;

		//! Auxiliary variable r_kj.
		double rkj;

		//! Auxiliary variable R_ij (x component).
		double Rxij;

		//! Auxiliary variable R_ij (y component).
		double Ryij

		//! Auxiliary variable R_ij (z component).
		double Rzij;

		//! Auxiliary variable R_ij squared.
		double Rij2

		//! Auxiliary variable r_ij.
		double rij;

		//! Auxiliary variable u_kj.
		double ukj;

		//! Auxiliary variable v_kj.
		double vkj;

		//! Auxiliary variable w_kj.
		double wkj;

		//! Auxiliary variable.
		double aux1;

		//! Yet another auxiliary variable.
		double aux2;

		//! Auxiliary variable for integration (x component).
		double auxv_x_integ;

		//! Auxiliary variable for integration (y component).
		double auxv_y_integ;

		//! Auxiliary variable for integration (z component).
		double auxv_z_integ;

		//! Auxiliary variable for integration.
		double aux3_integ;

		//! Auxiliary variable, square root of integration result.
		double aux3Sqrt_integ;

		//! Auxiliary variable, delta (x component).
		double auxv_x_delta;

		//! Auxiliary variable, delta (y component).
		double auxv_y_delta;

		//! Auxiliary variable, delta (z component).
		double auxv_z_delta;

		//! Auxiliary variable, delta.
		double aux3_delta;

		//! Auxiliary variable, square root of delta.
		double aux3Sqrt_delta;

		//! Auxiliary variable.
		double aa;

		//! Image method kernel function, 1st order currently
		void force_pol (Snapshot*, size_t, size_t, size_t);
						
		//! Pre-simulation hook
		void PreSimulation(Snapshot*, const CVList&) override;

		//! Post-integration hook
		void PostIntegration(Snapshot* snapshot, const CVList&) override;

		//! Post-simulation hook
		void PostSimulation(Snapshot*, const CVList&) override;

		//! \copydoc Serializable::Serialize()
		void Serialize(Json::Value& json) const override {};
	
		//! Destructor
		~COPSSImage() {}
	};
}
