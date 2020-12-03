#ifndef ___CELL_H
#define ___CELL_H

#include <iostream>
using namespace std;
#define _USE_MATH_DEFINES
#include <cmath>
#include "ap.h"
#include "subcell.h"

class CCells {
private:
  static const double vc;
  static const double stim;
  static const double stimduration;
  double jparam;
  double dt;
public:
  int n;
  CCells(int nn);
  virtual ~CCells();
  void Pace(void);
  void Pace_VClamp(double clampv);
  void ApClamp(double t, double bcl, double Vmin = -80.0, double Vmax = 30.0, double APD = 0);//BCL ms
  CActionPotential *ap;
  CSubcell *sc;
  double *v;
  double *ci;
  double *cnsr;
  double *st;

  void setjparam(double newjp);
  void setdt(double dtt);
  double getdt(void) { return dt; }
  double getvc(void) { return vc; }
  double getstim(void) { return stim; }
  double getstimduration(void) { return stimduration; }
//  CCells& operator=(const CCells& cell);


  double getvmax(void) { return 10; }
  double getvmin(void) { return -80; }
  void CopyParam(const CCells& cell){return;}//do nothing

};

#endif /* ___CELL_H */
