#ifndef __COLVAR_H__
#define __COLVAR_H__

namespace BIAS_NS{
class colVar {

public: 

  colVar();
  ~colVar();
 
  void assign_colVar(LAMMPS_NS::LAMMPS *lmp, int id, int type, int flag, double **pos, double **frc, int nlist1, int *list1, int nlist2, int *list2);

  //these values are obtained in the function pointer, and accessible from the wrapping code
  double val; //xi

  //an array of size nlist1/nlist2 
  //containing differentials with respect 
  //to type 1/type 2 coordinates
  double *fctr[3];
  double **der1; 
  double **der2; 

  void calc(){
    (this->*typeCalc)();
  }

  double getval(){ return val; }
  int l1(int jj){ return list1[jj];};
  int l2(int jj){ return list2[jj];};
  int nl1(void){ return nlist1;};
  int nl2(void){ return nlist2;};
  
private:  
  //cv info
  int id, type, flag;
  int xflag[3];
  int nlist1, nlist2, *list1, *list2;

  //interface parameters
  double **pos;
  double **frc;
  //image?
  double  boxx[3];
  double  hbox[3];
  LAMMPS_NS::LAMMPS *lmp;

  //define computation via function pointer
  void assignTypeAndPointer(int type, int flag);	
  //void (colVar::* typeCalc)() = NULL;
  void (colVar::* typeCalc)();
  //nlist2 = 0 and LIST2 should be NULL (they're ignored anyway);
  //  dimension CV is applied to set in assignment above.
  void position();
  
  //dimension of distance set in the assignement above
  void distance();

  //nlist2 = 0 and LIST2 should be NULL (they're ignored anyway);
  //energy! need to define a LAMMPS compute using atoms in LIST1 and extract the compute
  void init_energy();
  void energy();

  //radius of gyration
  void gyration();
};
}
#endif //__COLVAR_H__
