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
	//Image method to include polarization corrections into electrostatic interactions
	//for systems where dielectric objects are embedded in dielectric continuum
	class COPSSImage : public Constraint
	{
	private:
		//! Dielectric constant of polarizable particles
		double _einner;

		//! Where non polarizable particles start
		int _ion_type_start;
		
		//! Dielectric constant of outside continuum		
		double _eouter;

		//! unit conversion constant
		double _qqrd2e;
		
		double xlo = 0.0;
		
		double xhi = 1.0;
	
		int ngauss = 5;

		double _xg0[5] = {-0.9061798459386640,-0.5384693101056831,0.00000000000000000,0.5384693101056831,0.9061798459386640};

   		double _wg0[5] = {0.2369268850561891,0.4786286704993665,0.5688888888888889,0.4786286704993665,0.2369268850561891};	
		
		std::vector<double> _atomTypeRadius;
	
	public:
		COPSSImage(boost::mpi::communicator& comm, 
			   unsigned int frequency, 
			   double einner, 
			   int ion_type_start,
			   std::vector<double> atomTypeRadius):
		Constraint(frequency, comm), _einner(einner), _ion_type_start(ion_type_start), _eouter(0), _qqrd2e(1), _atomTypeRadius(atomTypeRadius){}
		
		// gauss integration auxiliary params
	                //auxiliary variables for image kernel functions
         	double _e, _ginv;

                double Rxkj, Rykj, Rzkj, Rkj2, rkj;

                double Rxij, Ryij, Rzij, Rij2, rij;

                double ukj, vkj, wkj;

                double aux1, aux2;

                double auxv_x_integ, auxv_y_integ, auxv_z_integ, aux3_integ, aux3Sqrt_integ;

                double auxv_x_delta, auxv_y_delta, auxv_z_delta, aux3_delta, aux3Sqrt_delta;

                double aa;

		//Image method kernel function, 1st order currently
		void force_pol (Snapshot*, size_t, size_t, size_t);
						
		// Pre-simulation hook
		void PreSimulation(Snapshot*, const CVList&) override;

		// Post-integration hook
		void PostIntegration(Snapshot* snapshot, const CVList&) override;

		// Post-simulation hook
		void PostSimulation(Snapshot*, const CVList&) override;

		void Serialize(Json::Value&) const override {};
	
		// Destructor
		~COPSSImage() {}
	};
}
