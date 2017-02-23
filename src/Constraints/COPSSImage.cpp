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
//#define DEBUG_SNAPSHOT_
//#define DEBUG_TMP_
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
  	        size_t nlocal_ = snapshot->GetNumAtoms();

		#ifdef DEBUG_SNAPSHOT_
		fprintf(stdout, "Image method start runnning ... \n");
	        
		auto& forces = snapshot->GetForces();
		auto& positions = snapshot->GetPositions();
		auto& charges = snapshot->GetCharges();

		fprintf(stdout, "check ids: ");
		for (size_t i = 0; i < nlocal_; i++)
		{
			fprintf(stdout, "%d\n; ", ids[i]);
		}

		fprintf(stdout, "check positions: ");
		for (size_t i=0; i< nlocal_; i++)
		{			
			for (size_t j=0; j<3; j++)
			{
				fprintf(stdout,"%f; ", positions[i][j]);
			}
		}

		fprintf(stdout, "\ncheck forces: ");
		for (size_t i=0; i< nlocal_; i++)
		{			
			for (size_t j=0; j<3; j++)
			{
				fprintf(stdout,"%f; ", forces[i][j]);
			}
		}
		
		fprintf(stdout, "\ncheck types: ");
		for (size_t i=0; i< nlocal_; i++)
		{
			fprintf(stdout, "%d; ", types[i]);
		}

		fprintf(stdout, "\ncheck charges: ");
		for (size_t i=0; i< nlocal_; i++)
		{
			fprintf(stdout, "%f; ", charges[i]);
		}
		
		fprintf(stdout, "\ncheck radius: ");
		for (size_t i=0; i< nlocal_; i++)
		{	
				fprintf(stdout,"%f; ", atomTypeRadius_[types[i]-1]);
		}

		fprintf(stdout, "\ncheck nlocal: %d\n", nlocal_);
		fprintf(stdout, "check eouter: %f\n", eouter_);
		fprintf(stdout, "check qqrd2e: %f\n", qqrd2e_);
		fprintf(stdout, "check einner: %f\n", einner_);
		fprintf(stdout, "check ion-type-start: %d\n", ion_type_start_);
		#endif

  		
		//iterate over all permutations {i;j;k}, j has to be polarizable
	
		for (size_t j=0; j < nlocal_; j++) if( types[j] < ion_type_start_ )
		{
			for(size_t i=0; i < nlocal_; i++) if(i != j) 
    			{ 
				 for(size_t k=0; k < nlocal_; k++) if (k != j)
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
                auto& types_  = snapshot->GetAtomTypes();
                auto& charges_ = snapshot->GetCharges();
                auto& positions_ = snapshot->GetPositions();
                auto& forces_ = snapshot->GetForces();
		 eouter_ = snapshot->GetDielectric();
                qqrd2e_ = snapshot->Getqqrd2e();	
   	 	aa = atomTypeRadius_[types_[jth]-1];
	 	// helper variables
	 	e_ = (1.0 - einner_ / eouter_) / (1.0 + einner_ / eouter_);
         	ginv_ = (1.0 + einner_ / eouter_);	
	 	//==================kernel function starts========================
	 	// Rxkj, Rykj, Rzkj: vector points from j-th to k-th, in x, y, z direction respectively.
	 	Rxkj = positions_[kth][0] - positions_[jth][0];
    	 	Rykj = positions_[kth][1] - positions_[jth][1];
    	 	Rzkj = positions_[kth][2] - positions_[jth][2];
         	// rkj: module of the vector between j-th and k-th
		Rkj2 = Rxkj*Rxkj + Rykj*Rykj + Rzkj*Rzkj;
         	rkj = sqrt (Rkj2);
  	 	// ukj, vkj, wkj: unit vector points from j-th to k-th, in x, y, z direction repectively
	 	ukj = Rxkj / rkj;
   	 	vkj = Rykj / rkj;
   	 	wkj = Rzkj / rkj;
	 	// Rxij, Ryij, Rzij: vector points from i-th to j-th, in x, y, z direction respectively
	 	Rxij = positions_[ith][0] - positions_[jth][0];
   	 	Ryij = positions_[ith][1] - positions_[jth][1];
    	 	Rzij = positions_[ith][2] - positions_[jth][2];
   	 	Rij2 = Rxij * Rxij + Ryij * Ryij + Rzij * Rzij;
    	        rij = sqrt(Rij2);
	 	// Auxiliary variables
   	 	aux1 = e_ * aa / rkj;
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
		std::vector<double> xg_(5);
		xlo = 0.0;
	 	xhi = 1.0;
   	 	for (int ig = 0; ig < ngauss; ++ig)
   		{
     	   		xg_[ig] = 0.5 * (xhi - xlo) * xg0_[ig] + 0.5 * (xhi + xlo); 
		      	auxv_x_integ = Rxij - ( pow(xg_[ig] , ginv_) * aux2 * ukj );
      			auxv_y_integ = Ryij - ( pow(xg_[ig] , ginv_) * aux2 * vkj );
      			auxv_z_integ = Rzij - ( pow(xg_[ig] , ginv_) * aux2 * wkj );
      			aux3Sqrt_integ = sqrt(auxv_x_integ * auxv_x_integ +
 					      auxv_y_integ * auxv_y_integ +
				      	      auxv_z_integ * auxv_z_integ);
      			aux3_integ = aux3Sqrt_integ * aux3Sqrt_integ * aux3Sqrt_integ;
			force_pol[0] += 0.5 * (xhi-xlo) * (- auxv_x_integ / aux3_integ * wg0_[ig]); 
      			force_pol[1] += 0.5 * (xhi-xlo) * (- auxv_y_integ / aux3_integ * wg0_[ig]);
      			force_pol[2] += 0.5 * (xhi-xlo) * (- auxv_z_integ / aux3_integ * wg0_[ig]);
    	 	}
		// multiple by prefactor in eq.29
	 	force_pol[0] *= ( 0.5 * qqrd2e_ * aux1 * charges_[ith] * charges_[kth] );
   	 	force_pol[1] *= ( 0.5 * qqrd2e_ * aux1 * charges_[ith] * charges_[kth] );
   	 	force_pol[2] *= ( 0.5 * qqrd2e_ * aux1 * charges_[ith] * charges_[kth] );
		
		for (size_t k = 0; k <3; ++k)
		{
			forces_[ith][k] += 2 * force_pol[k];
			forces_[jth][k] -= 2 * force_pol[k];
		}
	} // end force_pol function		
}//end namespace

