#ifndef ___CELL_H
#define ___CELL_H

#include <iostream>
using namespace std;
#define _USE_MATH_DEFINES
#include <cmath>
#include "ap.h"
#include "subcell.h"

class CCell
{
private:
  static const double vc;
  static const double stim;
  static const double stimduration;
  double jparam;
  double dt;
public:
  CCell(void);
  virtual ~CCell();
  void pace(double stim=0);
  void pace_vclamp(double clampv);
  void prepare(double bcl=300, int iter=0);
  void apclamp(double t, double bcl, double Vmin=-80.0, double Vmax=30.0, double APD=0);//BCL ms
  CActionPotential *ap;
  CSubcell *sc;
  double &v;//,&ci,&cj;

  double ci,cnsr;
  double computeavecnsr(void){cnsr=sc->computeavecnsr();return cnsr;}
  double computeaveci(void){ci=sc->computeaveci();return ci;}

  double setjparam(double newjp){jparam=newjp;ap->setjparam(jparam);return newjp;}
  double setdt(double dtt){dt=dtt;sc->setdt(dt);ap->setdt(dt);return dt;}
  double getdt(void){return dt;}
  double getvc(void){return vc;}
  double getstim(void){return stim;}
  double getstimduration(void){return stimduration;}
  CCell& operator=(const CCell& cell);


  double getvmax(void){return 10;}
  double getvmin(void){return -80;}

};

#endif /* ___CELL_H */
