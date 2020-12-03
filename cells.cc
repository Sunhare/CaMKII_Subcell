#include "cells.h"

const double CCells::vc = -80;
const double CCells::stim = 40;
const double CCells::stimduration = 1;

CCells::CCells(int nn) {
  n=nn;
  ap = new CActionPotential[n];
  sc = new CSubcell[n];
  v = new double[n];
  ci = new double[n];
  cnsr = new double[n];
  st = new double[n];
  for (int id = 0; id < n; id++) {
    ci[id] = id*0.5;
    cnsr[id] = id;
    v[id] = -86;
    st[id]=0;
  }
  //initial conditions
  dt = 0.1;
  jparam = 1;
}
CCells::~CCells() {
  delete[] sc;
  delete[] ap;
  delete[] v;
  delete[] ci;
  delete[] cnsr;
  delete[] st;
}

void CCells::ApClamp(double t, double T, double Vmin, double Vmax, double APD) {
  for (int id = 0; id < n; id++) {
    ci[id] = 0;
    cnsr[id] = 0;
    ap[id].apclamp(t, T, Vmin, Vmax, APD);
    sc[id].pace(ap[id].v, ap[id].nai);
    v[id] = ap[id].v;
  }
}
void CCells::Pace(void) {
  for (int id = 0; id < n; id++) {
    ci[id] = 0;
    cnsr[id] = 0;
    sc[id].pace(ap[id].v, ap[id].nai);
    ap[id].pace(st[id], sc[id].calcjcalave() * 3, sc[id].calcjnacaave() * 3, sc[id].calcicabkave(), sc[id].calcislcapave(), sc[id].computeavecs());
    v[id] = ap[id].v;
  }
}
void CCells::Pace_VClamp(double clampv) {
  for (int id = 0; id < n; id++) {
    ci[id] = 0;
    cnsr[id] = 0;
    ap[id].v = clampv;
    sc[id].pace(ap[id].v, ap[id].nai);
    v[id] = ap[id].v;
  }
}
void CCells::setjparam(double newjp) {
  jparam = newjp;
  for (int id = 0; id < n; id++) {
    ap[id].setjparam(jparam);
  }
}
void CCells::setdt(double dtt) {
  dt = dtt;
  for (int id = 0; id < n; id++) {
    sc[id].setdt(dt);
    ap[id].setdt(dt);
  }
}
