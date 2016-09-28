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

#include "COPSSImage.h"
#include <math.h>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
//#define _DEBUG_SNAPSHOT
//#define _DEBUG_TMP
namespace SSAGES
{	
	
	void COPSSImage::PreSimulation(Snapshot*, const CVList&)
	{
	//nothing to be done in presimulation
	}

	void COPSSImage::PostIntegration(Snapshot* snapshot, const CVList&)
	{

		// Gather information
  		auto& types = snapshot->GetAtomTypes();
  	        size_t _nlocal = snapshot->GetNumAtoms();

		#ifdef _DEBUG_SNAPSHOT
		fprintf(stdout, "Image method start runnning ... \n");
	        
		auto& forces = snapshot->GetForces();
		auto& positions = snapshot->GetPositions();
		auto& charges = snapshot->GetCharges();

		fprintf(stdout, "check ids: ");
		for (size_t i = 0; i < _nlocal; i++)
		{
			fprintf(stdout, "%d\n; ", ids[i]);
		}

		fprintf(stdout, "check positions: ");
		for (size_t i=0; i< _nlocal; i++)
		{			
			for (size_t j=0; j<3; j++)
			{
				fprintf(stdout,"%f; ", positions[i][j]);
			}
		}

		fprintf(stdout, "\ncheck forces: ");
		for (size_t i=0; i< _nlocal; i++)
		{			
			for (size_t j=0; j<3; j++)
			{
				fprintf(stdout,"%f; ", forces[i][j]);
			}
		}
		
		fprintf(stdout, "\ncheck types: ");
		for (size_t i=0; i< _nlocal; i++)
		{
			fprintf(stdout, "%d; ", types[i]);
		}

		fprintf(stdout, "\ncheck charges: ");
		for (size_t i=0; i< _nlocal; i++)
		{
			fprintf(stdout, "%f; ", charges[i]);
		}
		
		fprintf(stdout, "\ncheck radius: ");
		for (size_t i=0; i< _nlocal; i++)
		{	
				fprintf(stdout,"%f; ", _atomTypeRadius[types[i]-1]);
		}

		fprintf(stdout, "\ncheck nlocal: %d\n", _nlocal);
		fprintf(stdout, "check eouter: %f\n", _eouter);
		fprintf(stdout, "check qqrd2e: %f\n", _qqrd2e);
		fprintf(stdout, "check einner: %f\n", _einner);
		fprintf(stdout, "check ion-type-start: %d\n", _ion_type_start);
		#endif

  		
		//iterate over all permutations {i;j;k}, j has to be polarizable
	
		for (size_t j=0; j < _nlocal; j++) if( types[j] < _ion_type_start )
		{
			for(size_t i=0; i < _nlocal; i++) if(i != j) 
    			{ 
				 for(size_t k=0; k < _nlocal; k++) if (k != j)
      				 {  
					//Interaction: [I] <---- [J] <-------[k]
					//F[i] = f{i;j;k}
					//F[j] = -(f{i;j;k} + f{k;j;i})
					//F[k] = f{k;j;i}   	
					if (i == k)
					{
						force_pol(snapshot, i, j, k);
					}
					else if (i < k)
					{
						force_pol(snapshot, i, j, k);
						force_pol(snapshot, k, j, i);
					}
 				}	
			}
		}


	}

	void COPSSImage::PostSimulation(Snapshot*, const CVList&)
	{
	//nothing to be done in post simulation
	}
		
	// force acting on ith particle in a kernel {ith; jth; kth}
	void COPSSImage::force_pol (Snapshot* snapshot, 
				    size_t ith,
				    size_t jth,
				    size_t kth)
	{
		//fprintf(stdout, "compute: (%d, %d, %d)\n", ith, jth, kth);
                auto& _types  = snapshot->GetAtomTypes();
                auto& _charges = snapshot->GetCharges();
                auto& _positions = snapshot->GetPositions();
                auto& _forces = snapshot->GetForces();
		 _eouter = snapshot->GetDielectric();
                _qqrd2e = snapshot->Getqqrd2e();	
   	 	aa = _atomTypeRadius[_types[jth]-1];
	 	// helper variables
	 	_e = (1.0 - _einner / _eouter) / (1.0 + _einner / _eouter);
         	_ginv = (1.0 + _einner / _eouter);	
	 	//==================kernel function starts========================
	 	// Rxkj, Rykj, Rzkj: vector points from j-th to k-th, in x, y, z direction respectively.
	 	Rxkj = _positions[kth][0] - _positions[jth][0];
    	 	Rykj = _positions[kth][1] - _positions[jth][1];
    	 	Rzkj = _positions[kth][2] - _positions[jth][2];
         	// rkj: module of the vector between j-th and k-th
		Rkj2 = Rxkj*Rxkj + Rykj*Rykj + Rzkj*Rzkj;
         	rkj = sqrt (Rkj2);
  	 	// ukj, vkj, wkj: unit vector points from j-th to k-th, in x, y, z direction repectively
	 	ukj = Rxkj / rkj;
   	 	vkj = Rykj / rkj;
   	 	wkj = Rzkj / rkj;
	 	// Rxij, Ryij, Rzij: vector points from i-th to j-th, in x, y, z direction respectively
	 	Rxij = _positions[ith][0] - _positions[jth][0];
   	 	Ryij = _positions[ith][1] - _positions[jth][1];
    	 	Rzij = _positions[ith][2] - _positions[jth][2];
   	 	Rij2 = Rxij * Rxij + Ryij * Ryij + Rzij * Rzij;
    	        rij = sqrt(Rij2);
	 	// Auxiliary variables
   	 	aux1 = _e * aa / rkj;
	 	aux2 = aa * aa / rkj;
	 	// Tmp variable for return value
    	 	std::vector<double> force_pol(3);
	 	//===============delta_s,1 term of equation.29 in method paper
	 	auxv_x_delta = Rxij - aux2 * ukj;
    	 	auxv_y_delta = Ryij - aux2 * vkj;
         	auxv_z_delta = Rzij - aux2 * wkj;
   	 	aux3Sqrt_delta = sqrt(auxv_x_delta * auxv_x_delta + 
               		                     auxv_y_delta * auxv_y_delta + 
				             auxv_z_delta * auxv_z_delta);
    	 	aux3_delta = aux3Sqrt_delta * aux3Sqrt_delta * aux3Sqrt_delta;
    	 	force_pol[0] = auxv_x_delta / aux3_delta;
         	force_pol[1] = auxv_y_delta / aux3_delta;
         	force_pol[2] = auxv_z_delta / aux3_delta;
		//===============integration term of equation.29 in method paper
		std::vector<double> _xg(5);
		xlo = 0.0;
	 	xhi = 1.0;
   	 	for (int ig = 0; ig < ngauss; ++ig)
   		{
     	   		_xg[ig] = 0.5 * (xhi - xlo) * _xg0[ig] + 0.5 * (xhi + xlo); 
		      	auxv_x_integ = Rxij - ( pow(_xg[ig] , _ginv) * aux2 * ukj );
      			auxv_y_integ = Ryij - ( pow(_xg[ig] , _ginv) * aux2 * vkj );
      			auxv_z_integ = Rzij - ( pow(_xg[ig] , _ginv) * aux2 * wkj );
      			aux3Sqrt_integ = sqrt(auxv_x_integ * auxv_x_integ +
 					      auxv_y_integ * auxv_y_integ +
				      	      auxv_z_integ * auxv_z_integ);
      			aux3_integ = aux3Sqrt_integ * aux3Sqrt_integ * aux3Sqrt_integ;
			force_pol[0] += 0.5 * (xhi-xlo) * (- auxv_x_integ / aux3_integ * _wg0[ig]); 
      			force_pol[1] += 0.5 * (xhi-xlo) * (- auxv_y_integ / aux3_integ * _wg0[ig]);
      			force_pol[2] += 0.5 * (xhi-xlo) * (- auxv_z_integ / aux3_integ * _wg0[ig]);
    	 	}
		// multiple by prefactor in eq.29
	 	force_pol[0] *= ( 0.5 * _qqrd2e * aux1 * _charges[ith] * _charges[kth] );
   	 	force_pol[1] *= ( 0.5 * _qqrd2e * aux1 * _charges[ith] * _charges[kth] );
   	 	force_pol[2] *= ( 0.5 * _qqrd2e * aux1 * _charges[ith] * _charges[kth] );
		
		for (size_t k = 0; k <3; ++k)
		{
			_forces[ith][k] += 2 * force_pol[k];
			_forces[jth][k] -= 2 * force_pol[k];
		}
	} // end force_pol function		
}//end namespace

