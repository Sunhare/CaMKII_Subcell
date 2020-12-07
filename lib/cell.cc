#include "cell.h"

const double CCell::vc = -80;
const double CCell::stim = 40;
const double CCell::stimduration = 1;

CCell::CCell(void) : ap(new CActionPotential), v(ap->v) {
  sc = new CSubcell;
  //initial conditions
  dt = 0.1;
  jparam = 1;
  ci = 0;
  cnsr = 0;
}
CCell::~CCell() {
  delete sc;
  delete ap;
}
void CCell::prepare(double bcl, int iter) {

  if (iter == 0) {
    double dciold = 0;
    double dciold2 = 0;
    bool first = false;
    int Tn = bcl * 30 / dt, bcln = bcl / dt, durn = stimduration / dt;
    for (int tn = 0; tn < Tn; tn++) {
      double t = tn*dt;
      if (tn%bcln < durn) {
        if (first) {
          computeaveci();
          if (fabs(ci - dciold2) < 0.05 && t > bcl * 5) {
            break;
          }
          dciold2 = dciold;
          dciold = ci;
          first = false;
        }
        pace(stim);
      }
      else {
        first = true;
        pace();
      }
    }
  }
  else {
    int Tn = bcl*iter / dt, bcln = bcl / dt, durn = stimduration / dt;
    for (int tn = 0; tn < Tn; tn++) {
      if (tn%bcln < durn)
        pace(stim);
      else
        pace();
    }
  }
}

CCell& CCell::operator=(const CCell& cell) {
  if (&cell != this) {
    *ap = *cell.ap;
    *sc = *cell.sc;
    jparam = cell.jparam;
    ci = cell.ci;
    cnsr = cell.cnsr;
    ap->setjparam(jparam);
    dt = cell.dt;
    ap->setdt(dt);
    sc->setdt(dt);
  }
  return(*this);
}

void CCell::apclamp(double t, double T, double Vmin, double Vmax, double APD) {
  ci = 0;
  cnsr = 0;
  ap->apclamp(t, T, Vmin, Vmax, APD);
  sc->pace(ap->v, ap->nai);
}
void CCell::pace(double Istim) {
  ci = 0;
  cnsr = 0;
  sc->pace(ap->v, ap->nai);
  ap->pace(Istim, sc->calcjcalave() * 3, sc->calcjnacaave() * 3, sc->calcicabkave(), sc->calcislcapave(), sc->computeavecs());
}
void CCell::pace_vclamp(double clampv) {
  ci = 0;
  cnsr = 0;
  ap->v = clampv;
  sc->pace(ap->v, ap->nai);
}

