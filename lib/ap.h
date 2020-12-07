#ifndef ___CACTIONPOTENTIAL_H
#define ___CACTIONPOTENTIAL_H

// #define ___USE_VAR_FOR_CONST //use variables for Gto Gks Gkr etc. instead of constants (more memory)

#include <iostream>
using namespace std;
#define _USE_MATH_DEFINES
#include <cmath>

class CActionPotential{
private:
  double jparam;//tauj*jparam
  void pacex(double stim,double jcal, double jnaca, double icabk, double islcap, double cai);
  static const int N=13;
  static const double vc;
  static const double stim;
  static const double stimduration;
  static const double temp;// temperature (K)
  static const double xxr;//
  static const double xf;// Faraday's constant
  static const double frt;

  double comp_ina (void);
  double comp_ikr(void);
  double comp_iks(double ci);
  double comp_ik1(void);
  double comp_ito(void);
  double comp_inak(void);

  // additional currents


  double comp_inabk(void);
  double comp_icl(double cac);  //***
  double comp_iclbk(void);
  //  double comp_icabk(double cac);->compute in subcell
  //  double comp_islcap(double cac);->compute in subcell


  double vold;
  double hode,hpde;

public:
  void pace(double st,double jcal, double jnaca, double icabk, double islcap, double cai);
  void apclamp(double t, double bcl, double Vmin=-80, double Vmax=15, double apd=0, bool pulse=false);
  double setjparam(double newjp){jparam=newjp;return newjp;}
  double setdt(double dtt){hpde=dtt;return hpde;}
  double getdt(void){return hpde;}
  //  int getdim(void){return N;}
  double getvc(void){return vc;}
  double getstim(void){return stim;}
  double getstimduration(void){return stimduration;}
  CActionPotential(void);
  virtual ~CActionPotential();
  double *y;
  double &xm,&xh,&xj,&xr,&xs1,&xs2,&xtos,&ytos,&v, &nai,&xtof,&ytof,&rtos;

  double _inaca,_ica,_iks,_ikr,_itof,_itos,_ik1,_ina,_inak;
  double _inabk,_icl,_iclbk,_icabk,_islcap;

  double gtos;// ito slow conductance 
  double gtof;// ito fast conductance 
  double gks;
  double gkr;
  double gna;// sodium conductance (mS/micro F) 
  double gk1;// Ik1 conductance
  double gnak;

  double xnao;//mM external Na
  double xki;//mM internal K
  double xko;//mM external K
  double cao;//mM external Ca

  double gnabk;
  double gcl;
  double gclbk;

  CActionPotential& operator=(const CActionPotential& ap);

};
#endif /* ___CACTIONPOTENTIAL_H */
