#pragma once

#include "Meta.h"

using namespace SSAGES

//constructor: must define _width, _height

//destructor

//essential functions

void Meta::PreSimulation(){

  //initialize my cv
  int ii;
  for(ii = 0; ii < dim; ii++)
    mycv_.push_back( mycvlist_[ii].Evaluate() );

}

void Meta::PostIntegration(){
  
  //if you heard something
  if(event) addHill();
  //always
  calcBiasForce();
  chainRule();
  
}


void Meta::addHill(){

  int ii;
  for(ii = 0; ii < dim; ii++)
    mycv_[ii] = mycvlist_[ii].Evaluate();

  Hill newhill = Hill(mycv_,width_,height_);
  hills_.push_back( newhill )

}

void Meta::calcBiasForce(){
  
  //initialize bias and derivative
  int ii;
  bias_ = 0;
  for(ii = 0; ii < dim; ii++)
    myder_[ii] = 0;

  for(ii = 0; ii < hills_.size(); ii++)
    oneBiasForce(hills[ii].center, mycv_, &bias_, &der[0]_, height_);

}

//cleanup with new nomenclature
void Meta::oneBiasForce(const vector<double>& loc, const vector<double>& cv, double *bias, double *der, double height) {
  
  int ii, jj;
  double tbias = 1;
  double *tder;
  double *dx;
  
  tder = new double[loc.size()];
  dx   = new double[loc.size()];

  for(ii = 0; ii < loc.size(); ii++){
    tder[ii] = 1;
    dx[ii] = cv[ii] - loc[ii];
    tbias *= gaussian(dx[ii],sigmas[ii]);
  }

  for(ii = 0; ii < loc.size(); ii++){
    for(jj = 0; jj < loc.size(); jj++){
      if(jj != ii) tder[ii] *= gaussian(dx[jj],sigmas[jj]);
      else         tder[ii] *= gaussianDerv(dx[jj],sigmas[jj]);
    }
  }

  *bias += height*tbias;
  for(ii = 0; ii < loc.size(); ii++)  der[ii] += height*tder[ii];
  
  delete [] tder;
  delete [] dx;
} 

