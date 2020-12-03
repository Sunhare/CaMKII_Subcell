//icpc uc.cpp -o ucc -openmp -O3
#define ___REC_CURRENTS
#define ___USE_VAR_FOR_CONST

#ifdef _OPENMP
#include <omp.h>
#endif

#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;
#include "subcell.h"
#include "recsubcell.h"
#include "ap.h"
#include "log.h"

inline unsigned int xorshift2(unsigned int *xx, unsigned int *yy, unsigned int *zz, unsigned int *ww)
{
  unsigned int t=(*xx^(*xx<<11));*xx=*yy;*yy=*zz;*zz=*ww;
  return( *ww=(*ww^(*ww>>19))^(t^(t>>8)) );
}


int main(int argc, char* argv[]) {
  Logging llog;

  //cell size
//    const int nx=65;    const int ny=27;    const int nz=11;
  const int nx = 10;    const int ny = 10;    const int nz = 10;
  const int nn = nx * ny * nz;



  //random number seed
  int r = 100;
  cout << r << endl;
  char tmpchar[255];
  sprintf(tmpchar, "seed=%d", r);
  llog.Note(tmpchar);

  //init CSubcell
  CSubcell sc(nx, ny, nz, 1, 1);
  sc.srand(r);
  sc.init(0.01, 500);


  //parameter setting
  sc.setgleak(0);
  sc.setvup(sc.getvup() * 1);
  sc.setgca(0.0);
  sc.settautr(sc.gettautr() * 10);

  sc.NCXalpha = 0.11;

  double pcl = 1000.0; //1 Hz


  sc.setdt(0.1);
  double dt = sc.getdt();
  sc.setKu(5.0);
  sc.setKb(0.005);
  sc.settauu(1250 * 1);
  sc.settaub(5 * 0.1);
  sc.settauc1(2.0);
  sc.settauc2(0.3);
  sc.setKcp(100);
  sc.setpedk12(0.000001);
  sc.setpedk43(0.000001);
  sc.sethh(10);
  sc.setKK(1400);


  for (int i = 0; i < nn; i++) {
    sc.nryr[i] = 100;
    sc.vp[i] = 0.00126;//uniform vp
  }


  //generate hetero cluster dist
  std::cout << "Generating hetero cluster distribution... " << std::endl;
  for (int i = 0; i < nn; i++) {
    const double mean = 100;
    const double std = 50;
    double res;
    do {
      //Gaussian
      double x1, x2, w;
      do {
        x1 = 2.0 * xorshift2(&sc.xsx[0], &sc.xsy[0], &sc.xsz[0], &sc.xsw[0])/(double)(UINT_MAX) - 1.0;
        x2 = 2.0 * xorshift2(&sc.xsx[0], &sc.xsy[0], &sc.xsz[0], &sc.xsw[0])/(double)(UINT_MAX) - 1.0;
        w = x1 * x1 + x2 * x2;
      } while (w >= 1.0);
      w = sqrt((-2.0 * log(w)) / w);
      double y1 = x1 * w;
      //double y2 = x2 * w;
      res = y1 * std + mean;
    } while (res < 10 || res>300);
    sc.nryr[i] = res;
  }



  //init CRecSubcell
  CActionPotential ap;
  RECSUBCELLPARAM param;

  CRecSubcell rec(&sc, &ap, &param);



  //get steady state
  int Tn = 10 * pcl / dt;//total simulation time
  int ttn = 1;//global time counter
  // std::cout << "Getting steady state.. " << std::endl;
  // for (int tn = 0; tn < Tn; tn++) {
  //   sc.pace(-80, 5.0);
  //   if (tn % 20000 == 0) {
  //     cout << setprecision(10) << tn * dt / pcl << "\t" << rec.computeavecnsr() << endl;
  //   }
  // }



  double prevaveci[nn];
  double prevavecs[nn];
  double prevavecp[nn];
  for (int i = 0; i < nn; i++) {
    prevaveci[i] = 0;
    prevavecs[i] = 0;
    prevavecp[i] = 0;
  }

  int num_sims = 1;
  for (int itr = 0; itr < num_sims; itr++) {
    cout << itr << endl;

    double aveci[nn];
    double avecs[nn];
    double avecp[nn];
    for (int i = 0; i < nn; i++) {
      aveci[i] = 0;
      avecs[i] = 0;
      avecp[i] = 0;
    }

    //main loop

    double v;
    int nbeats = 10;

    Tn = nbeats * pcl / dt;
    for (int tn = 1; tn <= Tn; tn++, ttn++) {
      double t = ttn * dt;

      double t_relative = (tn%(int(pcl/dt)))*dt;

      // std::cout << "Before VC" << std::endl;
      sc.voltage_clamp(t, 20, true);
      // std::cout << "After VC" << std::endl;

    // sc.pace(v, 12.0);
      // if (tn % 20==0){
      //   sc.printsc0d(tn*dt,v,"0d_results.txt");

    // }
      
      if (tn % 20000 == 0) {
        cout << setprecision(10) << ttn * dt / pcl << "\t" << rec.computeavecnsr() << endl;
      }
      for (int i = 0; i < nn; i++) {
        aveci[i] += sc.ci[i];
        avecs[i] += sc.cs[i];
        avecp[i] += sc.cp[i];
      }



      
    }
    //end main loop

    if (!(sc.cp[0] > 0)) {
      cout << "NaN  itr=" << itr << endl;
      break;
    }

    //record average values
    ofstream osave("ave.txt");
    for (int i = 0; i < nn; i++) {
      prevaveci[i] = prevaveci[i] * itr + aveci[i] / Tn;
      prevaveci[i] /= (itr + 1);
      prevavecs[i] = prevavecs[i] * itr + avecs[i] / Tn;
      prevavecs[i] /= (itr + 1);
      prevavecp[i] = prevavecp[i] * itr + avecp[i] / Tn;
      prevavecp[i] /= (itr + 1);
      cout << sc.nryr[i] << "\t" << sc.vp[i] << "\t" << sc.Jmaxx[i] << "\t" << prevaveci[i] << "\t" << prevavecs[i] << "\t" << prevavecp[i] << endl;
      osave << sc.nryr[i] << "\t" << sc.vp[i] << "\t" << sc.Jmaxx[i] << "\t" << prevaveci[i] << "\t" << prevavecs[i] << "\t" << prevavecp[i] << endl;
    }
  }


  return 0;
}
