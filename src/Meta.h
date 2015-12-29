#pragma once

#include "Method.h"
#include "CollectiveVariable.h"

namespace SSAGES
{
  class Meta : public Method
  {
  private:
    struct Hill {
      vector<double> center;
      vector<double> width;
      double height;
      Hill(const vector<double> & center, const vector<double> & sigma, double height):
       center(center), width(width), height(height);
    };
    vector<double> width_;
    vector<Hill>   hills_;

    void addHill(void);
    void calcBiasForce (void);
    void chainRule (void);
    void oneBiasForce (void);
    
    int dim;
    vector<double> mycv_;
    vector<double> myder_;
    double bias_;
    CVList mycvlist_;

  public: 

    Meta(unsigned int frequency) : Method(frequency) {
      //We need to inherit 
    }
    void PreSimulation(Snapshot* snapshot, const CVList& cvs);
    void PostIntegration(Snapshot* snapshot, const CVList& cvs);
    ~Meta(){
      //not sure what to put here currently.
    }
  }


      
