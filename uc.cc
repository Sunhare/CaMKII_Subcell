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
#include "lib/subcell.h"
#include "lib/recsubcell.h"
#include "lib/ap.h"
#include "lib/log.h"

//starting again

inline unsigned int xorshift2(unsigned int *xx, unsigned int *yy, unsigned int *zz, unsigned int *ww)
{
  unsigned int t=(*xx^(*xx<<11));*xx=*yy;*yy=*zz;*zz=*ww;
  return( *ww=(*ww^(*ww>>19))^(t^(t>>8)) );
}


int main(int argc, char* argv[]) {
  Logging llog;

  //Sanity check
  // ofstream CaM_out; 
  // CaM_out.open("CaMDyad.txt");

  #ifdef ___PTM
    ofstream CaM_out; 
    ofstream CaMKII_out;
    ofstream BAR_out;
    ofstream dydt_CaMKII_out;

    CaM_out.open("CaMDyad.txt");
    CaMKII_out.open("CaMKII_out.txt");
    dydt_CaMKII_out.open("dydt_CaMKII_out.txt");
    BAR_out.open("BAR_out.txt");
    cout << "PTM Simulation" << std::endl;
  #else
    cout << "No PTM Simulation" << std::endl;
  #endif

  //cell size
  // const int nx=65;  const int ny=27; const int nz=11;
  const int nx=10;  const int ny=10; const int nz=10;
  // const int nx= 5;  const int ny= 5; const int nz= 5;
  // const int nx= 3;  const int ny= 3; const int nz= 3;
  const int nn = nx * ny * nz;

//  omp_set_num_threads(16);

  //random number seed
  int r = 100;
  cout << r << endl;
  char tmpchar[255];
  sprintf(tmpchar, "seed=%d", r);
  // llog.Note(tmpchar);

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


  // sc.setdt(0.02);
  sc.setdt(0.001);
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
  param.filenameci = "ci.dat";//record cytosolic Ca2+ (1byte binary 0-255)
  param.filenamecimm = "cimm.txt";//record max & min values of cytosolic Ca2+
  param.filenamecnsr = "cnsr.dat";//record (network) SR Ca2+ (1byte binary 0-255)
  param.filenamecnsrmm = "cnsrmm.txt";//record max & min values of (network) SR Ca2+
  CRecSubcell rec(&sc, &ap, &param);

  // get steady state
  int Tn = 10 * 1000 / dt;//total simulation time
  int ttn = 1;//global time counter
  for (int tn = 0; tn < Tn; tn++) {
    sc.pace(-80, 5.0);
    if (tn % int(2.0*1000.0/dt) == 0) {
      cout << setprecision(10) << tn * dt / 1000.0 << "\t" << rec.computeavecnsr() << endl;
      rec.recci(Z_LAYER_AVERAGE);//record cytosolic Ca2+
      rec.reccnsr(Z_LAYER_AVERAGE);//record (network) SR Ca2+
    }
  }



  double prevaveci[nn];
  double prevavecs[nn];
  double prevavecp[nn];
  for (int i = 0; i < nn; i++) {
    prevaveci[i] = 0;
    prevavecs[i] = 0;
    prevavecp[i] = 0;
  }

  // int num_sims = 100;
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
    int num_beats = 50;
    // int num_beats = 10;
    Tn = num_beats * 1000 / dt;
    // Tn = num_beats * 10 / dt;
    for (int tn = 1; tn <= Tn; tn++, ttn++) {
      double t = ttn * dt;
      sc.pace(-80, 12.0);



      if (tn % int(50.0/dt) == 0) { 
        cout << setprecision(10) << ttn * dt / 1000.0 << "\t" << rec.computeavecnsr() << endl;
        rec.recci(Z_LAYER_AVERAGE);//record cytosolic Ca2+
        rec.reccnsr(Z_LAYER_AVERAGE);//record (network) SR Ca2+
      }
      for (int i = 0; i < nn; i++) {
        aveci[i] += sc.ci[i];
        avecs[i] += sc.cs[i];
        avecp[i] += sc.cp[i];
      }

     #ifdef ___PTM
        // if (tn % int(2.0*1000.0/dt) == 0) {
      if (tn % int(2.0*1000.0/dt) == 0) {
          CaM_out << 
          t                     << "\t" << 
          sc.cp[0]                << "\t" <<
          sc.CaM_dyad[0]          << "\t" << 
          sc.dydt_CaMDyad[0]      << "\t" << 
          sc.CaM_dyad[0]          << "\t" << 
          sc.Ca2CaM_dyad[0]       << "\t" << 
          sc.Ca4CaM_dyad[0]       << "\t" << 
          sc.CaMB_dyad[0]         << "\t" << 
          sc.Ca2CaMB_dyad[0]      << "\t" << 
          sc.Ca4CaMB_dyad[0]      << "\t" << 
          sc.Pb2_dyad[0]          << "\t" << 
          sc.Pb_dyad[0]           << "\t" << 
          sc.Pt_dyad[0]           << "\t" << 
          sc.Pt2_dyad[0]          << "\t" << 
          sc.Pa_dyad[0]           << "\t" << 
          sc.Ca4CaN_dyad[0]       << "\t" << 
          sc.CaMCa4CaN_dyad[0]    << "\t" << 
          sc.Ca2CaMCa4CaN_dyad[0] << "\t" << 
          sc.Ca4CaMCa4CaN_dyad[0] << "\t" << 
          std::endl;

          CaMKII_out << 
          t                     << "\t" << //1 
          sc.LCC_PKAp[0]          << "\t" << //2
          sc.LCC_CKdyadp[0]       << "\t" << //3
          sc.RyR2809p[0]          << "\t" << //4
          sc.RyR2815p[0]          << "\t" << //5
          sc.PLBT17p[0]           << "\t" << //6
          sc.LCC_CKslp[0]         << "\t" << //7
          sc.RYR_multiplier[0]    << "\t" << //8
          std::endl;

          dydt_CaMKII_out <<
          t               << "\t" << //1 
          sc.dydt_CaMKII[0] << "\t" << //2
          sc.dydt_CaMKII[1] << "\t" << //3
          sc.dydt_CaMKII[2] << "\t" << //4
          sc.dydt_CaMKII[3] << "\t" << //5
          sc.dydt_CaMKII[4] << "\t" << //6
          sc.dydt_CaMKII[5] << "\t" << //7
          std::endl;

          BAR_out <<
          t                     << "\t" << //1 
          sc.L[0]                 << "\t" << //2
          sc.B1AR[0]              << "\t" << //3
          sc.Gs[0]                << "\t" << //4
          sc.B1AR_ACT[0]          << "\t" << //5
          sc.B1AR_S464[0]         << "\t" << //6
          sc.B1AR_S301[0]         << "\t" << //7
          sc.GsaGTPtot[0]         << "\t" << //8
          sc.GsaGDP[0]            << "\t" << //9
          sc.GsBy[0]              << "\t" << //10
          sc.GsaGTP[0]            << "\t" << //11
          sc.Fsk[0]               << "\t" << //12
          sc.AC[0]                << "\t" << //13
          sc.PDE[0]               << "\t" << //14
          sc.IBMX[0]              << "\t" << //15
          sc.cAMPtot[0]           << "\t" << //16
          sc.cAMP[0]              << "\t" << //17
          sc.PKAC_I[0]            << "\t" << //18
          sc.PKAC_II[0]           << "\t" << //19
          sc.PLBp[0]              << "\t" << //20
          sc.Inhib1ptot[0]        << "\t" << //21
          sc.Inhib1p[0]           << "\t" << //22
          sc.PP1[0]               << "\t" << //23
          sc.LCCa_PKAp_whole[0]   << "\t" << //24
          sc.LCCb_PKAp_whole[0]   << "\t" << //25
          sc.RyR_PKAp_whole[0]    << "\t" << //26
          sc.TnI_PKAp_whole[0]    << "\t" << //27
          sc.IKs_PKAn[0]          << "\t" << //28
          sc.Yotiao_KCQN1[0]      << "\t" << //29 
          sc.IKs_PKAp_whole[0]    << "\t" << //30
          sc.ICFTR_PKAp_whole[0]  << "\t" << //31
          sc.PLM_PKAp_whole[0]    << "\t" << //32
          sc.Myo_PKAp_whole[0]    << "\t" << //33
          sc.IKr_PKAn[0]          << "\t" << //34
          sc.Yotiao_hERG[0]       << "\t" << //35
          sc.IKr_PKAp_whole[0]    << "\t" << //36
          sc.IClCa_PKAp_whole[0]  << "\t" << //37
          std::endl;
        }
      #endif
    }
    if (!(sc.cp[0] > 0)) {
      cout << "NaN  itr=" << itr << endl;
      break;
    }

    //record average values
    #ifdef ___PTM
      ofstream osave("ave_ptm.txt");
    #else
      ofstream osave("ave.txt");
    #endif 

    for (int i = 0; i < nn; i++) {
      prevaveci[i] = prevaveci[i] * itr + aveci[i] / Tn;
      prevaveci[i] /= (itr + 1);
      prevavecs[i] = prevavecs[i] * itr + avecs[i] / Tn;
      prevavecs[i] /= (itr + 1);
      prevavecp[i] = prevavecp[i] * itr + avecp[i] / Tn;
      prevavecp[i] /= (itr + 1);
      // cout << sc.nryr[i] << "\t" << sc.vp[i] << "\t" << sc.Jmaxx[i] << "\t" << prevaveci[i] << "\t" << prevavecs[i] << "\t" << prevavecp[i] << endl;
      osave << sc.nryr[i] << "\t" << sc.vp[i] << "\t" << sc.Jmaxx[i] << "\t" << prevaveci[i] << "\t" << prevavecs[i] << "\t" << prevavecp[i] << endl;
    }


  }


  return 0;
}
