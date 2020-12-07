
#ifndef ___RECSUBCELL_H
#define ___RECSUBCELL_H

#include <iostream>
using namespace std;
#define _USE_MATH_DEFINES
#include <cmath>
#include <fstream>
#include <cfloat>
#include <string>
#include <cstdlib>

#include "subcell.h"
#include "ap.h"
#include "cell.h"

#define ALL_LAYER -1
#define Z_LAYER_AVERAGE 999999

class RECSUBCELLPARAM{
private:
public:
  string filenametv;//t,v,ciave,csave,cpave
  string filenameci;// spatial ci
  string filenamecimm;// spatial ci
  string filenamecs;// spatial ci
  string filenamecsmm;// spatial ci
  string filenamecp;// spatial ci
  string filenamecpmm;// spatial ci
  string filenamecnsr;// spatial cnsr
  string filenamecnsrmm;// spatial cnsr
  string filenamecjsr;// spatial cnsr
  string filenamecjsrmm;// spatial cnsr
  string filenameryr1;// spatial cnsr
  string filenameryr2;// spatial cnsr
  string filenameryr3;// spatial cnsr
  string filenameryr4;// spatial cnsr
  string filenamecurrent;// ica, ncx, iup, irel, ileak
  //      string filenameapdcip;// apd, cip, apdamp, cipamp
  int step;
  int interval;
  RECSUBCELLPARAM(void){interval=1;step=1;filenametv="";
  filenameci="";filenamecimm="";filenamecs="";filenamecsmm="";filenamecp="";filenamecpmm="";
  filenamecnsr="";filenamecnsrmm="";filenamecjsr="";filenamecjsrmm="";
  filenameryr1="";filenameryr2="";filenameryr3="";filenameryr4="";
  filenamecurrent="";}
};

class CRecSubcell
{
private:

  CSubcell *sc;
  CCell *cell;

  CActionPotential *ap;
  int n,nn,nx,ny,nz,nnx,nny,nnz,nxny,nnxnny;

  char *rec;
  double *ave;
  ofstream *tvout;bool IsRectvout;
  ofstream *ciout; ofstream *cimmout; bool IsRecciout;
  ofstream *csout; ofstream *csmmout; bool IsReccsout;
  ofstream *cpout; ofstream *cpmmout; bool IsReccpout;
  ofstream *cnsrout; ofstream *cnsrmmout; bool IsReccnsrout;
  ofstream *cjsrout; ofstream *cjsrmmout; bool IsReccjsrout;
  ofstream *ryr1out; bool IsRecryr1out;
  ofstream *ryr2out; bool IsRecryr2out;
  ofstream *ryr3out; bool IsRecryr3out;
  ofstream *ryr4out; bool IsRecryr4out;
  ofstream *currentout;bool IsReccurrentout;
  //  ofstream *apdcipout;bool IsRecapdcipout;

  string cifilename;
  string cimmfilename;
  string csfilename;
  string csmmfilename;
  string cpfilename;
  string cpmmfilename;
  string cnsrfilename;
  string cnsrmmfilename;
  string cjsrfilename;
  string cjsrmmfilename;
  string ryr1filename;
  string ryr2filename;
  string ryr3filename;
  string ryr4filename;
  int finestep;
  int step;
public:
  CRecSubcell(class CSubcell *scx, class CActionPotential *apx, class RECSUBCELLPARAM *param);
  CRecSubcell(class CCell *cellx, class RECSUBCELLPARAM *param);
  CRecSubcell(void);
  void setcell(class CCell *cellx, class RECSUBCELLPARAM *param);

  virtual ~CRecSubcell();

  double computeaveci(void);
  double computeavecs(void);
  double computeavecp(void);
  double computeavesr(void);
  double computeavecjsr(void);
  double computeavecnsr(void);
  void recci(int layer=ALL_LAYER);
  void reccs(int layer=ALL_LAYER);
  void reccp(int layer=ALL_LAYER);
  void reccnsr(int layer=ALL_LAYER);
  void reccjsr(int layer=ALL_LAYER);
#ifndef ___DETERMINISTIC
  void recryr1(int layer=ALL_LAYER);
  void recryr2(int layer=ALL_LAYER);
  void recryr3(int layer=ALL_LAYER);
  void recryr4(int layer=ALL_LAYER);
#endif
  void reccurrent(double t);
  void rectv(double t);
};
#endif
