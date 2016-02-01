#include "mpi.h"
#include "src/lammps.h"
#include "src/library.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "colVar.h"

using namespace LAMMPS_NS;
using namespace BIAS_NS;

colVar::colVar() {

  this->lmp     = NULL;
  this->pos     = this->frc = NULL;
  this->list1   = this->list2 = NULL;
  this->der1    = this->der2  = NULL;

  this->id      = 0;
  this->type    = 0;
  this->flag    = 0;

  this->nlist1  = this->nlist2 = 0;


};

void colVar::assign_colVar(LAMMPS *lmp, int id, int type, int flag, double **pos, double **frc, int nlist1, int *list1, int nlist2, int *list2){

  int ii;

  this->lmp     = lmp;
  this->id      = id;
  this->type    = type;
  this->flag    = flag;
  this->pos     = pos;
  this->frc     = frc;

  this->nlist1  = nlist1;
  this->list1   = list1;
  if(nlist1) {
    // this->list1 = new int[nlist1];
    // for(ii = 0; ii < nlist1; ii++)
    //   this->list1[ii] = list1[ii];

    der1 = new double*[nlist1];
    for(ii = 0; ii < nlist1; ii++)
      der1[ii] = new double[3];
  }

  this->nlist2  = nlist2;
  this->list2   = list2;
  if(nlist2) {
    // this->list2 = new int[nlist2];
    // for(ii = 0; ii < nlist2; ii++)
    //   this->list2[ii] = list2[ii];
    
    der2 = new double*[nlist2];
    for(ii = 0; ii < nlist2; ii++)
      der2[ii] = new double[3];
  }

  assignTypeAndPointer(type,flag);

  boxx[0] = *((double *) lammps_extract_global(lmp,"boxxhi"))
    - *((double *) lammps_extract_global(lmp,"boxxlo"));
  boxx[1] = *((double *) lammps_extract_global(lmp,"boxyhi"))
    - *((double *) lammps_extract_global(lmp,"boxylo"));
  boxx[2] = *((double *) lammps_extract_global(lmp,"boxzhi")) 
    - *((double *) lammps_extract_global(lmp,"boxzlo"));
  hbox[0] = boxx[0]/2.;
  hbox[1] = boxx[1]/2.;
  hbox[2] = boxx[2]/2.;
}

colVar::~colVar(){

  int ii;

  if(nlist1) {
    for(ii = 0; ii < nlist1; ii++)
      delete [] der1[ii];

    delete [] der1;
    delete [] list1;
  }

  if(nlist2) {
    for(ii = 0; ii < nlist2; ii++)
      delete [] der2[ii];

    delete [] der2;
    delete [] list2;
  }

  if(type == 3){
    char cmd[1024];
    sprintf(cmd,"group _egroup_%d delete",id);
    lammps_command(lmp,cmd);
    sprintf(cmd,"uncompute _reducepair_%d",id);
    lammps_command(lmp,cmd);
    sprintf(cmd,"uncompute _paircomp_%d",id);
    lammps_command(lmp,cmd);
  }

}

void colVar::assignTypeAndPointer(int type, int flag){
  this->type = type;
  switch(type){
  case 1:
    switch(flag){
    case 1: 
      xflag[0] = 1; xflag[1] = xflag[2] = 0;
      break;
    case 2:
      xflag[1] = 1; xflag[0] = xflag[2] = 0;
      break;
    case 3:
      xflag[2] = 1; xflag[1] = xflag[0] = 0;
      break;
    default:
      fprintf(stderr,"Bad position CV definition");
    }
    typeCalc = &colVar::position;
    break;
  case 2:
    switch(flag){
    case 1: 
      xflag[0] = 1; xflag[1] = xflag[2] = 0;
      break;
    case 2:
      xflag[1] = 1; xflag[0] = xflag[2] = 0;
      break;
    case 3:
      xflag[2] = 1; xflag[1] = xflag[0] = 0;
      break;
    default:
      xflag[0] = xflag[1] = xflag[2] = 1;
      break;
    }
    typeCalc = &colVar::distance;
    break;
  case 3:
    init_energy();
    typeCalc = &colVar::energy;
    break;
  case 4:
    typeCalc = &colVar::gyration;
    break;
  default:
    fprintf(stderr,"Bad colVar definition\n");
    exit(-1);
    break;
  }
}

void colVar::position(){

  int ii,jj;
  double cv = 0;

  for(ii = 0; ii < nlist1; ii++)
    for(jj = 0; jj < 3; jj++)
	cv += xflag[jj] * pos[list1[ii]][jj];

  val = cv/nlist1;

  //derivatives: dr/dx = delta_i{flag}
  for(ii = 0; ii < nlist1; ii++)
    for(jj = 0; jj < 3; jj++)
      der1[ii][jj] = xflag[jj]/nlist1;

  //list2 ignored!
}

void colVar::distance(){

  int ii, jj;
  double cv = 0;
  double cm1[3] = {0}, cm2[3] = {0}, dx[3] = {0};

  //should use unwrapped coordinates in future
  for(ii = 0; ii < nlist1; ii++)
    for(jj = 0; jj < 3; jj++)
      cm1[jj] += pos[list1[ii]][jj]/nlist1;

  for(ii = 0; ii < nlist2; ii++)
    for(jj = 0; jj < 3; jj++)
      cm2[jj] += pos[list2[ii]][jj]/nlist2;

  //hard coded! make flexible [can use lammps data structures or outside help]
  for(ii = 0; ii < 3; ii++) {
    dx[ii] = cm1[ii] - cm2[ii];
    if(dx[ii] >  hbox[ii]) dx[ii]-=boxx[ii];
    if(dx[ii] < -hbox[ii]) dx[ii]+=boxx[ii];
    cv += dx[ii]*dx[ii];
  }

  cv  = sqrt(cv);
  val = cv;

  for(ii = 0; ii < nlist1; ii++)
    for(jj = 0; jj < 3; jj++)
      der1[ii][jj] = (dx[jj])/cv;
  
  for(ii = 0; ii < nlist2; ii++)
    for(jj = 0; jj < 3; jj++)
      der2[ii][jj] = -(dx[jj])/cv;

}
	
void colVar::init_energy() {

  int ii;

  char  *cmd = new char[16*nlist1 + 256];
  char  *atm = new char[16];
  sprintf(cmd, "group _egroup_%d id ",id);
  for(ii = 0; ii < nlist1; ii++){
    sprintf(atm,"%d ",list1[ii]+1);
    strcat(cmd,atm);
  }
  lammps_command(lmp,cmd);
  
  sprintf(cmd,"compute _paircomp_%d _egroup_%d pe/atom",id,id);
  lammps_command(lmp,cmd);

  sprintf(cmd,"compute _reducepair_%d _egroup_%d reduce sum c__paircomp_%d",id,id,id);
  lammps_command(lmp,cmd);

  sprintf(cmd,"thermo_style custom temp ke pe c__reducepair_%d",id);
  lammps_command(lmp,cmd);

  sprintf(cmd,"thermo_modify norm no");
  lammps_command(lmp,cmd);
  
  delete [] cmd;
  delete [] atm;
  
}

void colVar::energy(){

  int ii, jj;
  char  arg[1024];
  sprintf(arg,"_reducepair_%d",id);
  double energy = *((double *) lammps_extract_compute(lmp,arg,0,0));

  val = energy;

  //check sign of this derivative
  for(ii = 0; ii < nlist1; ii++)
    for(jj = 0; jj < 3; jj++)
      der1[ii][jj] = frc[list1[ii]][jj];
}
    
void colVar::gyration(){


  int ii, jj, kk;
  double cv = 0;
  double dx[3] = {0};
  double rg2 = 0, tmp;

  double nsq = nlist1*nlist1;

  //construct a cluster starting from original site using nearest-neighbor images.
  //gyration
  for(ii = 0; ii < nlist1; ii++) {
    for(jj = 0; jj < nlist1; jj++) {
      for(kk = 0; kk < 3; kk++) {
	dx[kk]=pos[list1[ii]][kk]-pos[list1[jj]][kk];
	if(dx[kk] >  hbox[kk]) dx[kk]-=boxx[kk];
	if(dx[kk] < -hbox[kk]) dx[kk]+=boxx[kk];
	rg2 += 0.5 * dx[kk] * dx[kk];
      }
    }
  }
  
  rg2 /= nsq;

  val = rg2;

  //derivative
  for(ii = 0; ii < nlist1; ii++) {
    for(kk = 0; kk < 3; kk++) {
      der1[ii][kk] = 0;
    }
    for(jj = 0; jj < nlist1; jj++) {
      for(kk = 0; kk < 3; kk++) {
	dx[kk]=pos[list1[ii]][kk]-pos[list1[jj]][kk];
	if(dx[kk] >  hbox[kk]) dx[kk]-=boxx[kk];
	if(dx[kk] < -hbox[kk]) dx[kk]+=boxx[kk];
	der1[ii][kk] += dx[kk]/nsq;
	der1[jj][kk] -= dx[kk]/nsq;
      }
    }
  }
}

