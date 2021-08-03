#include "subcell.h"
//starting again

#ifdef ___DEBUG
#include <iomanip>
#include <limits>
#include <string>
#endif

//cell parameters
const double CSubcell::vjsr=0.02;
const double CSubcell::vp_ave=0.00126;



#define WT 0
#define OE 1
#define KO 2

inline unsigned int xorshift(unsigned int *xx, unsigned int *yy, unsigned int *zz, unsigned int *ww)
{
  unsigned int t=(*xx^(*xx<<11));*xx=*yy;*yy=*zz;*zz=*ww;
  return( *ww=(*ww^(*ww>>19))^(t^(t>>8)) );
}


CSubcell::CSubcell(int sizex, int sizey, int sizez, int fmesh, double xii)
{
  if (fmesh == 1)
    dt=0.1;
  else if (fmesh == 3)
    dt=0.01;
  else if (fmesh == 5)
    dt=0.005;
  else{
    cerr <<"fine mesh incorrect !\n";
    exit(1);
  }

  nx=sizex;
  ny=sizey;
  nz=sizez;
  finemesh=fmesh;
  xi=xii;
  n=nn=0;

  //parameters
  cao=1.8;
  vup=0.3;
  vnaca=21.0;
#ifdef ___USE_ORG_PARAM
  Jmax=0.0147;
  Ku=3.8*0.0001;
  Kb=5*0.00001;
  tauu=125.00;
  taub=5.0;
  tauc1=1.0;
  tauc2=1.0;
  hh=23;
  KK=850;
  gcabk=0;
  qslcap=0;
#else
  Jmax=0.0147*18;
  Ku=5.0;
  Kb=0.005;
  tauu=1250.0;
  taub=0.5;
  tauc1=2.0;
  tauc2=0.3;
  hh=10;
  KK=1400;
  gcabk=0.0002513;
  qslcap=2.35;
#endif

  MaxSR=15;
  MinSR=1;
  ec50SR=450;
  hkosrca=2.5;



  kup=0.123;
  KNSR=1700;

  gca=1.0;
  NoCaL=4;

  gleak=1.035*0.00001;


  taup=0.022;
  BCSQN=400;

  tautr=5.0;

  Kc=600.0;
  nM=15;
  nD=35;
  rhoinf=5000;

  BCSQN0=400;

  //#ifdef ___SIGMOID
  Kcp=100;
  pedk12=0.000001;
  pedk43=0.000001;
  //#endif

  seed=0;
  initialized=false;
  bc=0;
#ifdef ___NCX
  NCXalpha=0;
#endif
#ifdef ___DEBUG
  bSTOP=false;
#endif

#ifdef ___EGTA
  BEGTA=350.0;
#endif

}
void CSubcell::init(double initci, double initcj)
{
  nxny=nx*ny;
  n=nx*ny*nz;

  finemesh3=finemesh*finemesh*finemesh;
  nnx=finemesh*nx;
  nny=finemesh*ny;
  nnz=finemesh*nz;
  nnxnny=nnx*nny;
  nn=nnx*nny*nnz;

  //cell parameters
  vi=0.5/finemesh3;
  vs=0.025/finemesh3;
  vnsr=0.025/finemesh3;


  tausT=xi*1.42/(finemesh*finemesh);
  tausL=xi*3.4/(finemesh*finemesh);
  taumninv=2/tausL+4/tausT;

  tauiT=xi*2.93/(finemesh*finemesh);
  tauiL=xi*2.32/(finemesh*finemesh);
  taunsrT=xi*7.2/(finemesh*finemesh);
  taunsrL=xi*24.0/(finemesh*finemesh);


  tausi=0.1/(finemesh*finemesh);

  ci=new double [nn];
  cs=new double [nn];
  cp=new double [n];
  cjsr=new double [n];
  cnsr=new double [nn];
  cati=new double [nn];
  //  cats=new double [nn];

#ifdef ___PTM
  allocate_memory_all_PTM_vars(n);
  init_const_parameters_PTM();
#endif 

#ifdef ___DETERMINISTIC
  c1=new double [n];
  c2=new double [n];
  i1ca=new double [n];
  i1ba=new double [n];
  i2ca=new double [n];
  i2ba=new double [n];

  fryr1=new double [n];
  fryr2=new double [n];
  fryr3=new double [n];
#else
  y=new int [n*NoCaL];
  ryr1=new int [n];
  ryr2=new int [n];
  ryr3=new int [n];
  nryr=new int [n];
#endif

  vp=new double [n];
  Jmaxx=new double [n];
  cscp1=new double [nn];
  cscp2=new double [nn];

  Ici=new double [nn];
  Icnsr=new double [nn];

#ifdef ___NO_CS_BUFFER
  csmn=new double [nn];
#else
  Ics=new double [nn];
  Idps=new double [nn];
#endif

#ifdef ___EGTA
  caEGTAi=new double [nn];
  caEGTAs=new double [nn];
#endif


  Itr=new double [nn];

  crupos=new int [n];

  //random number
  xsx=new unsigned int [n];
  xsy=new unsigned int [n];
  xsz=new unsigned int [n];
  xsw=new unsigned int [n];
#pragma ivdep
#pragma vector always
  for (int id=0;id<n;id++)
  {
    xsx[id]=123456789+id+seed;
    xsy[id]=362436069+id*100+seed*10;
    xsz[id]=521288629+id*1000+seed*100;
    xsw[id]=88675123+id*10000+seed*1000;
    for (int i=0;i<1000;i++)
      xorshift(&xsx[id], &xsy[id], &xsz[id], &xsw[id]);
  }

  //initial conditions

#pragma ivdep
#pragma vector always
  for (int id=0;id<nn;id++)
  {
    ci[id]=initci;
    cs[id]=initci;
    cnsr[id]=initcj;
    Icnsr[id]=0;

  #ifdef ___PTM
      init_state_variables_PTMs(id);
  #endif

  #ifdef ___NO_CS_BUFFER
      csmn[id]=0;
  #else
      Ics[id]=0;
      Idps[id]=0;
  #endif

    Ici[id]=0;

    cscp1[id]=0;
    cscp2[id]=0;
    Itr[id]=0;
  }
  resetBuffer();
#pragma ivdep
#pragma vector always
  for (int id=0;id<n;id++)
  {
    Jmaxx[id]=Jmax;
    cp[id]=initci;
    cjsr[id]=initcj;
#ifdef ___DETERMINISTIC
    c1[id]=1;
    c2[id]=0;
    i1ca[id]=0;
    i1ba[id]=0;
    i2ca[id]=0;
    i2ba[id]=0;

    fryr1[id]=0.03;
    fryr2[id]=0;
    fryr3[id]=0;

#else
    ryr1[id]=0+int(5.0*xorshift(&xsx[id], &xsy[id], &xsz[id], &xsw[id])/(double)(UINT_MAX));
    ryr2[id]=0;
    ryr3[id]=0;
    nryr[id]=100;
    for (int j=0;j<NoCaL;j++)
      y[id*NoCaL+j]=2;
#endif

#ifdef ___UNIFORM
    double r=0;
#else
    double r=calcvp(0,0.3,-0.8,0.8,id);//Gaussian distribution (0,0.3) range(-0.8~0.8)
#endif
    vp[id]=vp_ave*(1+r);

    int kk=id/nxny;
    kk=kk*finemesh+finemesh/2;
    int modi=id%nxny;
    int jj=modi/nx;
    jj=jj*finemesh+finemesh/2;
    int ii=modi%nx;
    ii=ii*finemesh+finemesh/2;
    crupos[id]=ii+jj*nnx+kk*nnxnny;

    cscp2[crupos[id]]=1/taup*vp[id]/vs;

  }

  iupave=icaave=incxave=irave=ileakave=icabkave=islcapave=0;
  initialized=true;
}
void CSubcell::srand(int sed)
{
  seed=sed;
#pragma ivdep
#pragma vector always
  for (int id=0;id<n;id++)
  {
    xsx[id]=123456789+id+seed;
    xsy[id]=362436069+id*100+seed*10;
    xsz[id]=521288629+id*1000+seed*100;
    xsw[id]=88675123+id*10000+seed*1000;
  }
  for (int id=0;id<n;id++)
    for (int i=0;i<1000;i++)
      xorshift(&xsx[id], &xsy[id], &xsz[id], &xsw[id]);
}

CSubcell::~CSubcell()
{
  if (initialized)
  {
    delarray();
  }
}
void CSubcell::delarray(void)
{
  delete [] ci;
  delete [] cs;
  delete [] cp;
  delete [] cjsr;
  delete [] cnsr;

#ifdef ___DETERMINISTIC
  delete [] c1;
  delete [] c2;
  delete [] i1ca;
  delete [] i1ba;
  delete [] i2ca;
  delete [] i2ba;

  delete [] fryr1;
  delete [] fryr2;
  delete [] fryr3;
#else
  delete [] ryr1;
  delete [] ryr2;
  delete [] ryr3;
  delete [] nryr;
  delete [] y;
#endif

  delete [] cati;
  //      delete [] cats;

  delete [] Jmaxx;
  delete [] vp;
  delete [] cscp1;
  delete [] cscp2;
  delete [] Itr;

  delete [] Ici;
  delete [] Icnsr;

#ifdef ___NO_CS_BUFFER
  delete [] csmn;
#else
  delete [] Ics;
  delete [] Idps;
#endif
  delete [] crupos;

  delete [] xsx;
  delete [] xsy;
  delete [] xsz;
  delete [] xsw;


  #ifdef ___PTM

    delete [] dydt_CaMDyad;
    delete [] dydt_CaMSL;
    delete [] dydt_CaMCyt;

    delete [] JCaCyt;
    delete [] JCaSL;
    delete [] JCaDyad;

    delete [] CaM_dyad;
    delete [] Ca2CaM_dyad;
    delete [] Ca4CaM_dyad;
    delete [] CaMB_dyad;
    delete [] Ca2CaMB_dyad;
    delete [] Ca4CaMB_dyad;
    delete [] Pb2_dyad;
    delete [] Pb_dyad;
    delete [] Pt_dyad;
    delete [] Pt2_dyad;
    delete [] Pa_dyad;
    delete [] Ca4CaN_dyad;
    delete [] CaMCa4CaN_dyad;
    delete [] Ca2CaMCa4CaN_dyad;
    delete [] Ca4CaMCa4CaN_dyad;

    delete [] CaM_sl;
    delete [] Ca2CaM_sl;
    delete [] Ca4CaM_sl;
    delete [] CaMB_sl;
    delete [] Ca2CaMB_sl;
    delete [] Ca4CaMB_sl;
    delete [] Pb2_sl;
    delete [] Pb_sl;
    delete [] Pt_sl;
    delete [] Pt2_sl;
    delete [] Pa_sl;
    delete [] Ca4CaN_sl;
    delete [] CaMCa4CaN_sl;
    delete [] Ca2CaMCa4CaN_sl;
    delete [] Ca4CaMCa4CaN_sl;

    delete [] CaM_cyt;
    delete [] Ca2CaM_cyt;
    delete [] Ca4CaM_cyt;
    delete [] CaMB_cyt;
    delete [] Ca2CaMB_cyt;
    delete [] Ca4CaMB_cyt;
    delete [] Pb2_cyt;
    delete [] Pb_cyt;
    delete [] Pt_cyt;
    delete [] Pt2_cyt;
    delete [] Pa_cyt;
    delete [] Ca4CaN_cyt;
    delete [] CaMCa4CaN_cyt;
    delete [] Ca2CaMCa4CaN_cyt;
    delete [] Ca4CaMCa4CaN_cyt;

    delete [] dydt_CaMKII;
    delete [] LCC_CKdyadp;
    delete [] LCC_CKslp;
    delete [] RyR2815p;
    delete [] PLBT17p;
    delete [] LCC_PKAp;
    delete [] RyR2809p;
    delete [] CaMKIIactDyad;
    delete [] CaMKIIactSL;
    delete [] PP1_PLB_avail;
    delete [] RyR_CKp;
    delete [] PLB_CKp;

    delete [] dydt_BAR;
    delete [] L;
    delete [] B1AR;
    delete [] Gs;
    delete [] B1AR_ACT;
    delete [] B1AR_S464;
    delete [] B1AR_S301;
    delete [] GsaGTPtot;
    delete [] GsaGDP;
    delete [] GsBy;
    delete [] GsaGTP;
    delete [] Fsk;
    delete [] AC;
    delete [] PDE;
    delete [] IBMX;
    delete [] cAMPtot;
    delete [] cAMP;
    delete [] PKAC_I;
    delete [] PKAC_II;
    delete [] PLBp;
    delete [] Inhib1ptot;
    delete [] Inhib1p;
    delete [] PP1;

    delete [] LCCa_PKAp_whole;
    delete [] LCCb_PKAp_whole;
    delete [] RyR_PKAp_whole;
    delete [] TnI_PKAp_whole;
    delete [] IKs_PKAn;
    delete [] Yotiao_KCQN1;
    delete [] IKs_PKAp_whole;
    delete [] ICFTR_PKAp_whole;
    delete [] PLM_PKAp_whole;
    delete [] Myo_PKAp_whole;
    delete [] IKr_PKAn;
    delete [] Yotiao_hERG;
    delete [] IKr_PKAp_whole;
    delete [] IClCa_PKAp_whole;

    delete [] LCCa_PKAp;
    delete [] LCCb_PKAp;
    delete [] PLB_PKAn;
    delete [] RyR_PKAp;
    delete [] TnI_PKAp;
    delete [] IKs_PKAp;
    delete [] ICFTR_PKAp;
    delete [] PLM_PKAp;
    delete [] Myo_PKAp;
    delete [] IKr_PKAp;
    delete [] IClCa_PKAp;
  #endif
}

CSubcell& CSubcell::operator=(const CSubcell& sc)
{
  if (&sc==this)return(*this);
  if (initialized)
  {
    delarray();
  }
  //constructor
  dt=sc.dt;
  nx=sc.nx;
  ny=sc.ny;
  nz=sc.nz;
  finemesh=sc.finemesh;
  xi=sc.xi;

  cao=sc.cao;
  vup=sc.vup;
  kup=sc.kup;
  KNSR=sc.KNSR;
  vnaca=sc.vnaca;
  Jmax=sc.Jmax;
  gca=sc.gca;

  gcabk=sc.gcabk;
  qslcap=sc.qslcap;
  gleak=sc.gleak;
  BCSQN=sc.BCSQN;

  Kc=sc.Kc;
  nM=sc.nM;
  nD=sc.nD;
  hh=sc.hh;
  KK=sc.KK;
  rhoinf=sc.rhoinf;

  Ku=sc.Ku;
  Kb=sc.Kb;
  tauu=sc.tauu;
  taub=sc.taub;
  tauc1=sc.tauc1;
  tauc2=sc.tauc2;
  BCSQN0=sc.BCSQN0;

  //#ifdef ___SIGMOID
  Kcp=sc.Kcp;
  pedk12=sc.pedk12;
  pedk43=sc.pedk43;
  //#endif
  NoCaL=sc.NoCaL;
#ifdef ___NCX
  NCXalpha=sc.NCXalpha;
#endif

  MaxSR=sc.MaxSR;
  MinSR=sc.MinSR;
  ec50SR=sc.ec50SR;
  hkosrca=sc.hkosrca;


  init();

  //initial conditions
#pragma ivdep
#pragma vector always
  for (int id=0;id<nn;id++)
  {
    ci[id]=sc.ci[id];
    cs[id]=sc.cs[id];
    cati[id]=sc.cati[id];
    //      cats[id]=sc.cats[id];
    cnsr[id]=sc.cnsr[id];
    Icnsr[id]=sc.Icnsr[id];
#ifdef ___NO_CS_BUFFER
    csmn[id]=sc.csmn[id];
#else
    Ics[id]=sc.Ics[id];
    Idps[id]=sc.Idps[id];
#endif

    cscp1[id]=sc.cscp1[id];
    cscp2[id]=sc.cscp2[id];
    Itr[id]=sc.Itr[id];
  }
#pragma ivdep
#pragma vector always
  for (int id=0;id<n;id++)
  {
    Jmaxx[id]=sc.Jmaxx[id];
    cp[id]=sc.cp[id];
    cjsr[id]=sc.cjsr[id];
#ifdef ___DETERMINISTIC
    c1[id]=sc.c1[id];
    c2[id]=sc.c2[id];
    i1ca[id]=sc.i1ca[id];
    i1ba[id]=sc.i1ba[id];
    i2ca[id]=sc.i2ca[id];
    i2ba[id]=sc.i2ba[id];

    fryr1[id]=sc.fryr1[id];
    fryr2[id]=sc.fryr2[id];
    fryr3[id]=sc.fryr3[id];
#else
    ryr1[id]=sc.ryr1[id];
    ryr2[id]=sc.ryr2[id];
    ryr3[id]=sc.ryr3[id];
    nryr[id]=sc.nryr[id];
    for (int j=0;j<NoCaL;j++)
      y[id*NoCaL+j]=sc.y[id*NoCaL+j];
#endif
    vp[id]=sc.vp[id];
    cscp2[crupos[id]]=sc.cscp2[crupos[id]];
    Ici[id]=sc.Ici[id];
  }

#pragma ivdep
#pragma vector always
  for (int id=0;id<n;id++)
  {
    xsx[id]=sc.xsx[id];
    xsy[id]=sc.xsy[id];
    xsz[id]=sc.xsz[id];
    xsw[id]=sc.xsw[id];
  }
  return(*this);
}

void CSubcell::pace(double v, double nai)
{
  iupave=icaave=incxave=irave=ileakave=icabkave=islcapave=0;

#ifndef ___NO_DIFFUSION
  //set diffusion terms
  computeIci();//diffusion ci
  computeIcnsr();//diffusion cnsr
#ifdef ___NO_CS_BUFFER
  computecsmn();//diffusion cs
#else
  computeIcs();//diffusion cs
#endif
#endif

  const double F=96.5;
  const double R=8.314;
  const double T=308;
  const double rtf=R*T/F;//~26.5
  const double rtf2=R*T/(2*F);
  double z=v*F/(R*T);
  double za=z*2.0;

  //L-double Ca current (ca independent part)
  const double pca=11.9;
  double factor1=4.0*pca*F*F/(R*T);
  double factor=v*factor1;
  const double vth=0.0;

  //  const double s6=8.0;
  const double s6=4.0;//Juan code param

  double poinf=1.0/(1.0+exp(-(v-vth)/s6));
  const double taupo=1.0;
  double alpha=poinf/taupo;
  double beta=(1.0-poinf)/taupo;
  const double vx=-40;
  const double sx=3.0;
  double poi=1.0/(1.0+exp(-(v-vx)/sx));
  const double tau3=3.0;

  double xk3=(1.0-poi)/tau3;
  double xk3t=xk3;

  const double vy=-40.0;
  const double sy=4.0;
  double prv=1.0-1.0/(1.0+exp(-(v-vy)/sy));

  double recov=10.0+4954.0*exp(v/15.6);
  double tauba=(recov-450.0)*prv+450.0;

  const double vyr=-40.0;
  const double syr=11.32;
  double poix=1.0/(1.0+exp(-(v-vyr)/syr));


  //NCX (ca independent part)
  const double Kmnai=12.3;
  const double nao=136;
  const double Kmcai=3.59;
  const double ksat=0.27;
  const double eta=0.35;
  double t1=(Kmcai*0.001)*nao*nao*nao*(1+(nai*nai*nai/(Kmnai*Kmnai*Kmnai)));
  double x1a=exp(eta*z)*nai*nai*nai*cao;
  double x1b=exp((eta-1)*z)*nao*nao*nao;
  double x2=(1+ksat*exp((eta-1)*z));
  double x3=vnaca;


  double sumica=0;
  double sumir=0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+: sumica, sumir)
#endif
#pragma ivdep
#pragma vector always
  for (int id=0;id<n;id++)
  {
    //L-double Ca current
    double ica=(factor/(exp(za)-1.0))*(cp[id]*0.001*exp(za)-cao);
    if(fabs(za)<0.1)
    {
      ica=(factor1/(2.0*F/(R*T)))*(cp[id]*0.001*exp(za)-cao);
    }
    if (ica>0)ica=0;
    const double r1=0.30;
    const double r2=6.0;
    const double cat=0.5;
    double fca=1.0/(1.0+pow(double(cat/cp[id]),3));

    //      double s1=0.0182688*fca;
    double s1=0.02*fca;//Juan code param
    const double s1t=0.00195;
    //      double xk1=0.024168*fca;
    double xk1=0.03*fca;//Juan code param
    const double xk2=1.03615e-4;

    const double xk1t=0.00413;
    const double xk2t=0.00224;
    double s2=s1*(r1/r2)*(xk2/xk1);
    const double s2t=s1t*(r1/r2)*(xk2t/xk1t);


    //      const double tca=78.0329;
    const double tca=114;//Juan code param
    const double cpt=1.5;
    double tau_ca=tca/(1.0+pow((cp[id]/cpt),4))+1;

    double tauca=(recov-tau_ca)*prv+tau_ca;


    double xk6=fca*poix/tauca;
    double xk5=(1.0-poix)/tauca;

    double xk6t=poix/tauba;
    double xk5t=(1.0-poix)/tauba;

    double xk4=xk3*(alpha/beta)*(xk1/xk2)*(xk5/xk6);
    double xk4t=xk3t*(alpha/beta)*(xk1t/xk2t)*(xk5t/xk6t);


#ifdef ___DETERMINISTIC
    double capo=1.0-i1ca[id]-i2ca[id]-i1ba[id]-i2ba[id]-c1[id]-c2[id];
    double dc2= beta*c1[id]+xk5*i2ca[id]+xk5t*i2ba[id]-(xk6+xk6t+alpha)*c2[id];
    double dc1=alpha*c2[id]+xk2*i1ca[id]+xk2t*i1ba[id]+r2*capo-(beta+r1+xk1t+xk1)*c1[id];
    double di1ca=xk1*c1[id]+xk4*i2ca[id]+s1*capo-(xk3+xk2+s2)*i1ca[id];
    double di2ca=xk3*i1ca[id]+xk6*c2[id]-(xk5+xk4)*i2ca[id];
    double di1ba=xk1t*c1[id]+xk4t*i2ba[id]+s1t*capo-(xk3t+xk2t+s2t)*i1ba[id];
    double di2ba=xk3t*i1ba[id]+xk6t*c2[id]-(xk5t+xk4t)*i2ba[id];

    c1[id]+=dc1*dt;
    c2[id]+=dc2*dt;
    i1ca[id]+=di1ca*dt;
    i1ba[id]+=di1ba*dt;
    i2ca[id]+=di2ca*dt;
    i2ba[id]+=di2ba*dt;
    double Ica=gca*ica*capo*NoCaL;

#else
    int NL=0;
#pragma ivdep
#pragma vector always
    for (int j=0;j<NoCaL;j++)
    {
      //          I - I
      //          |   |  |
      //          C - C - O
      //          |   |  |
      //          I - I

      //          6 - 5
      //          |   |  |
      //          2 - 1 - 0
      //          |   |  |
      //          4 - 3
      double rr=xorshift(&xsx[id], &xsy[id], &xsz[id], &xsw[id])/(double)(UINT_MAX);

      rr/=dt;
      switch(y[id*NoCaL+j])
      {
      case 0:
        NL++;
        if (rr<r2)
          y[id*NoCaL+j]=1;
        else if (rr<(s1+r2))
          y[id*NoCaL+j]=5;
        else if (rr<(s1+r2+s1t))
          y[id*NoCaL+j]=3;
        break;
      case 1:
        if (rr<beta)
          y[id*NoCaL+j]=2;
        else if (rr<(r1+beta))
          y[id*NoCaL+j]=0;
        else if (rr<(r1+beta+xk1t))
          y[id*NoCaL+j]=3;
        else if (rr<(r1+beta+xk1t+xk1))
          y[id*NoCaL+j]=5;
        break;
      case 2:
        if (rr<xk6)
          y[id*NoCaL+j]=6;
        else if (rr<(xk6+xk6t))
          y[id*NoCaL+j]=4;
        else if (rr<(alpha+xk6t+xk6))
          y[id*NoCaL+j]=1;
        break;
      case 3:
        if (rr<xk3t)
          y[id*NoCaL+j]=4;
        else if (rr<(xk3t+xk2t))
          y[id*NoCaL+j]=1;
        else if (rr<(s2t+xk2t+xk3t))
          y[id*NoCaL+j]=0;
        break;
      case 4:
        if (rr<xk5t)
          y[id*NoCaL+j]=2;
        else if (rr<(xk5t+xk4t))
          y[id*NoCaL+j]=3;
        break;
      case 5:
        if (rr<xk3)
          y[id*NoCaL+j]=6;
        else if (rr<(xk3+xk2))
          y[id*NoCaL+j]=1;
        else if (rr<(s2+xk2+xk3))
          y[id*NoCaL+j]=0;
        break;
      case 6:
        if (rr<xk5)
          y[id*NoCaL+j]=2;
        else if (rr<(xk5+xk4))
          y[id*NoCaL+j]=5;
        break;
      }
    }
    double Ica=gca*ica*NL;
#endif

    sumica+=Ica;

    //release current Ir
#ifdef ___DETERMINISTIC
    double Po=fryr2[id]+fryr3[id];
#else
    double Po=(ryr2[id]+ryr3[id])/100.0;
#endif
    double Ir=Jmaxx[id]*Po*(cjsr[id]-cp[id])/vp[id];
    sumir+=Ir;
    //Diffusion from proximal space to submembrane space Idps
    double kr=Jmaxx[id]*Po/vp[id];

#ifdef ___NCX
    //NCX cleft space
    const double Kmcao=1.3;
    const double Kmnao=87.5;
    const double Kda=0.11;

    double Ka=1/(1+(Kda*Kda*Kda/(cp[id]*cp[id]*cp[id])));
    double t2=Kmnao*Kmnao*Kmnao*(cp[id]*0.001)*(1+cp[id]/Kmcai);
    double t3=Kmcao*nai*nai*nai+nai*nai*nai*cao+nao*nao*nao*(cp[id]*0.001);
    double jnaca=NCXalpha*Ka*x3*(x1a-x1b*(cp[id]*0.001))/((t1+t2+t3)*x2)*vs/vp[id];
    double newcp=(cs[crupos[id]]+taup*(kr*cjsr[id]-Ica+jnaca))/(1+taup*Jmaxx[id]*Po/vp[id]);
#else
    double newcp=(cs[crupos[id]]+taup*(kr*cjsr[id]-Ica))/(1+taup*Jmaxx[id]*Po/vp[id]);
#endif
    if (newcp<=0)newcp=0.00001;

    //Diffusion from NSR to JSR
    Itr[crupos[id]]=(cnsr[crupos[id]]-cjsr[id])/tautr;

    //Nearest-nighbor diffusive current Ici, Ics, IcNSR
    //RyR
    double rho=rhoinf*pow(double(cjsr[id]/1000.0),hh)/(pow(KK/1000,hh)+pow(double(cjsr[id]/1000),hh));
    if (rho<0.0000001)rho=0.0000001;/// to avoid div small (for cjsr<200)
    double MM=(sqrt(1+8*rho*BCSQN)-1)/(4*rho*BCSQN);
    double ncjsr=MM*nM+(1-MM)*nD;

    //      double rhopri=(hh*rhoinf*pow(double(cjsr[id]),double(hh-1))*(pow(double(KK),double(hh))+pow(double(cjsr[id]),double(hh)))-rhoinf*pow(double(cjsr[id]),double(hh))*hh*pow(double(cjsr[id]),double(hh-1)))/pow((pow(double(KK),double(hh))+pow(double(cjsr[id]),double(hh))),2);

    //      double rhopri=(hh*rhoinf*pow(double(cjsr[id]/1000.0),double(hh-1))*(pow(double(KK/1000.0),double(hh))+pow(double(cjsr[id]/1000.0),double(hh)))-rhoinf*pow(double(cjsr[id]/1000.0),double(hh))*hh*pow(double(cjsr[id]/1000.0),double(hh-1)))/pow((pow(double(KK/1000.0),double(hh))+pow(double(cjsr[id]/1000.0),(hh))),2)/1000.0;

    //gpu code
    double rhopri=(hh*rhoinf*pow(double(cjsr[id]/1000.0),double(hh-1))*(pow(double(KK/1000.0),double(hh))+pow(double(cjsr[id]/1000.0),double(hh)))-rhoinf*pow(double(cjsr[id]/1000.0),double(hh))*hh*pow(double(cjsr[id]/1000.0),double(hh-1)))/pow((pow(double(KK/1000.0),double(hh))+pow(double(cjsr[id]/1000.0),double(hh))),2);
    rhopri*=0.001;


    double dMdc=(((1.0/2.0)*pow(double(1+8*rho*BCSQN),double(-1.0/2.0))*8*rhopri*BCSQN*4*rho*BCSQN)-(sqrt(1+8*rho*BCSQN)-1)*4*rhopri*BCSQN)/pow(4*rho*BCSQN,2);
    double dndc=dMdc*(nM-nD);

    double Betajsr=1/(1+(Kc*BCSQN*ncjsr+dndc*(cjsr[id]*Kc+cjsr[id]*cjsr[id]))/((Kc+cjsr[id])*(Kc+cjsr[id])));
    //double m=MM*(nM*BCSQN-cB);


    //      double rhocp=rhoinf*pow(cp[id],hh)/(pow(KK,hh)+pow(cp[id],hh));
    //      double MMcp=(sqrt(1+8*rhocp*BCSQN)-1)/(4*rhocp*BCSQN);

    double cp2=cp[id]*cp[id];
#ifdef ___USE_ORG_PARAM
    double k12=Ku*cp2;
    double k43=Kb*cp2;
#else

#ifdef ___KOSRCA
    double kCaSR = MaxSR - (MaxSR-MinSR)/(1+pow((ec50SR/cjsr[id]),hkosrca));
    double koSRCa = 1/kCaSR;
#else
    const double koSRCa=1;
#endif
    double sgmd=cp2/(Kcp*Kcp+cp2);
    double k12=koSRCa*Ku*sgmd+pedk12;
    double k43=koSRCa*Kb*sgmd+pedk43;
#endif

#ifdef ___SPTM
    if (cp[id]>1.0) {
      k12=k12*2;
      k43=k43*2;
    }

#endif

#ifdef ___PTM 

    //Soltis Saucerman/Negroni PTM 
    double fCKII_RyR = (20.0 * RyR_CKp[id] / 3.0 - 1.0 / 3.0); //≈ 1.0005, Max = 19/3, 6
    double fPKA_RyR = (RyR_PKAp[id] * 1.025) + 0.9750; // ≈
    
    // if(id%3==0){
      // std::cout << id << ": dPa_dyad: " << dydt_CaMDyad[(id*15)+ 10]*dt << std::endl;;
      // std::cout << id << ": Pa_dyad: " << Pa_dyad[id] << std::endl;
      // std::cout << id << ": CaMKIIactDyad: " << CaMKIIactDyad[id] << std::endl;
      // std::cout << id << ": CaMKIIactSL: " << CaMKIIactSL[id] << std::endl;
      // std::cout << id << ": RyR2815p: " << RyR2815p[id] << std::endl;
      // std::cout << id << ": RyR_CKp: " << RyR_CKp[id] << std::endl;
      // std::cout << id << ": RyR_PKAp: " << RyR_CKp[id] << std::endl;
      // std::cout << id << ": fCKII_RYR: " << fCKII_RyR << std::endl;
      // std::cout << id << ": fPKA_RyR: " << fPKA_RyR << std::endl;
      // std::cout << std::endl;
      
    // }
    

    double FINE_TUNING = 1.0; //Theoretical parameter, choose any number ()
    double BASELINE_PARAMETER = 0.0;

    k12 *= (fCKII_RyR + fPKA_RyR - 1)* FINE_TUNING + BASELINE_PARAMETER; //Multiplier is always bigger than 1 
    k43 *= (fCKII_RyR + fPKA_RyR - 1)* FINE_TUNING + BASELINE_PARAMETER; // 

    //Currently, it's 0.66 to 1.05 (during pacing)
    //Need to adjust so that it's from 1.0 to something larger (based on fine tuning)

    RYR_multiplier[id] = (fCKII_RyR + fPKA_RyR - 1); //Sanity check, check if this is NAN

#endif


    double k14=MM/taub*BCSQN/BCSQN0;
    double k21=1/tauc1;
    double k23=MM/taub*BCSQN/BCSQN0;
    double k41=1/tauu;
    double k34=1/tauc2;
    double k32=k41*k12*k23*k34/(k43*k21*k14);

#ifdef ___DETERMINISTIC

    double fryr4=1.0-fryr1[id]-fryr2[id]-fryr3[id];
    double dfryr1=k21*fryr2[id]+k41*fryr4-(k12+k14)*fryr1[id];
    double dfryr2=k12*fryr1[id]+k32*fryr3[id]-(k21+k23)*fryr2[id];
    double dfryr3=k23*fryr2[id]+k43*fryr4-(k34+k32)*fryr3[id];

    fryr1[id]+=dfryr1*dt;
    fryr2[id]+=dfryr2*dt;
    fryr3[id]+=dfryr3*dt;

#else
    int ryr4=nryr[id]-ryr1[id]-ryr2[id]-ryr3[id];
    int ryr12=bino(ryr1[id],k12*dt,id);
    int ryr14=bino(ryr1[id],k14*dt,id);
    int ryr21=bino(ryr2[id],k21*dt,id);
    int ryr23=bino(ryr2[id],k23*dt,id);
    int ryr43=bino(ryr4,k43*dt,id);
    int ryr41=bino(ryr4,k41*dt,id);
    int ryr34=bino(ryr3[id],k34*dt,id);
    int ryr32=bino(ryr3[id],k32*dt,id);
    ryr1[id]=ryr1[id]-(ryr12+ryr14)+(ryr21+ryr41);
    ryr2[id]=ryr2[id]-(ryr21+ryr23)+(ryr12+ryr32);
    ryr3[id]=ryr3[id]-(ryr34+ryr32)+(ryr43+ryr23);
    if (ryr1[id]<0 ||ryr2[id]<0||ryr3[id]<0||ryr1[id]+ryr2[id]+ryr3[id]>nryr[id])
    {
      if (ryr1[id]<0)
      {
        if (xorshift(&xsx[id], &xsy[id], &xsz[id], &xsw[id])%2)
          ryr2[id]+=ryr1[id];
        ryr1[id]=0;
      }
      if (ryr2[id]<0)
      {
        if (xorshift(&xsx[id], &xsy[id], &xsz[id], &xsw[id])%2)
        {
          ryr1[id]+=ryr2[id];
          if (ryr1[id]<0)ryr1[id]=0;
        }
        else
          ryr3[id]+=ryr2[id];
        ryr2[id]=0;
      }
      if (ryr3[id]<0)
      {
        if (xorshift(&xsx[id], &xsy[id], &xsz[id], &xsw[id])%2)
        {
          ryr2[id]+=ryr3[id];
          if (ryr2[id]<0)ryr2[id]=0;
        }
        ryr3[id]=0;
      }
      if (ryr1[id]+ryr2[id]+ryr3[id]>nryr[id])
      {
        ryr4=nryr[id]-(ryr1[id]+ryr2[id]+ryr3[id]);
        if (xorshift(&xsx[id], &xsy[id], &xsz[id], &xsw[id])%2)
        {
          ryr3[id]+=ryr4;
          if (ryr3[id]<0)
          {
            ryr3[id]-=ryr4;
            ryr1[id]+=ryr4;
            if (ryr1[id]<0)
            {
              ryr1[id]-=ryr4;
              ryr2[id]+=ryr4;
            }
          }
        }
        else
        {
          ryr1[id]+=ryr4;
          if (ryr1[id]<0)
          {
            ryr1[id]-=ryr4;
            ryr3[id]+=ryr4;
            if (ryr3[id]<0)
            {
              ryr3[id]-=ryr4;
              ryr2[id]+=ryr4;
            }
          }
        }
      }
    }
#endif

    //update
    double dcjsr=Betajsr*(Itr[crupos[id]]-Ir*(vp[id]/vjsr));

#ifdef ___NO_CS_BUFFER
    cscp2[crupos[id]]=vp[id]/(vs*taup);
    cscp1[crupos[id]]=cp[id]*cscp2[crupos[id]];
    //      cscp1[crupos[id]]=cp[id]*vp[id]/vs/taup;
    //      cscp2[crupos[id]]=1/taup*vp[id]/vs;
#else
    Idps[crupos[id]]=(cp[id]-cs[crupos[id]])/taup;
#endif


#ifdef ___CONTRACTION //TODO implement Negroni 
#endif 

  //Main PTM ODEs
#ifdef ___PTM

    // TODO
    // for(int dtx = 0; dtx < 100; dtx++){
    //   // solve for signaling part with smaller time step
    //   // check if this works
    // }

  //Finer Timestep

  for(int dtx = 0; dtx < 100; dtx++){

  //CaMKII Equations
    CaMKIIactDyad[id] = CaMKIItotDyad*(Pb_dyad[id]+Pt_dyad[id]+Pt2_dyad[id]+Pa_dyad[id]); // Multiply total by fraction of activated CaMKII states (Pb, Pt, Pt2, Pa)
    CaMKIIactSL[id] = CaMKIItotSL*(Pb_sl[id]+Pt_sl[id]+Pt2_sl[id]+Pa_sl[id]);

    // PP1_PLB_avail[id] = PP1[id]/PP1_PLBtot + .0091;  // Active PP1 near PLB / total PP1 conc + basal value

    calc_dydt_CaMKII_ODEs(id);
    LCC_CKdyadp[id] = LCC_CKdyadp[id]/LCCtotDyad;   //136 fractional CaMKII-dependent LCC dyad phosphorylation
    RyR_CKp[id] = RyR2815p[id]/RyRtot;           //138 fractional CaMKII-dependent RyR phosphorylation
    PLB_CKp[id] = PLBT17p[id]/PLBtot;           //139 fractional CaMKII-dependent PLB phosphorylation
  
  //Beta Adrenergic/PKA equations
    calc_dydt_BAR_ODEs(id);
    LCCa_PKAp[id] = LCCa_PKAp_whole[id]/LCCtotBA;
    LCCb_PKAp[id] = LCCb_PKAp_whole[id]/LCCtotBA;
    PLB_PKAn[id] = (PLBtotBA - PLBp[id])/PLBtotBA; // non-phosphorylated PLB targets
    RyR_PKAp[id] = RyR_PKAp_whole[id]/RyRtotBA;
    TnI_PKAp[id] = TnI_PKAp_whole[id]/TnItotBA;
    IKs_PKAp[id] = IKs_PKAp_whole[id]/IKstotBA;
    ICFTR_PKAp[id] = ICFTR_PKAp_whole[id]/ICFTRtotBA;
    PLM_PKAp[id] = PLM_PKAp_whole[id]/PLMtotBA;
    Myo_PKAp[id] = Myo_PKAp_whole[id]/MyototBA;
    IKr_PKAp[id] = IKr_PKAp_whole[id]/IKrtotBA;
    IClCa_PKAp[id] = IClCa_PKAp_whole[id]/IClCatotBA;

    solve_ODE_CaM(id);
    solve_ODE_CaMKII(id);
    solve_ODE_BAR(id);
  }
#endif 

#ifdef ___CPDIFF
    //Instantaneous buffering functions
    const double KCAM=7.0;
    const double BCAM=24.0;
    const double KSR=0.6;
    const double BSR=47.0;
    const double KMCa=0.033;
    const double BMCa=140.0;
    const double KMMg=3.64;
    const double BMMg=140.0;
    const double KSLH=0.3;
    const double BSLH=13.4;

    double CAM=BCAM*KCAM/((cp[id]+KCAM)*(cp[id]+KCAM));
    double SR=BSR*KSR/((cp[id]+KSR)*(cp[id]+KSR));
    double MCa=BMCa*KMCa/((cp[id]+KMCa)*(cp[id]+KMCa));
    double MMg=BMMg*KMMg/((cp[id]+KMMg)*(cp[id]+KMMg));
    double SLH=BSLH*KSLH/((cp[id]+KSLH)*(cp[id]+KSLH)); // only for cs
    double Betap=1/(1+CAM+SR+MCa+MMg+SLH);


#ifdef ___DEBUG
    if(isnan(Idps[crupos[id]])) //(Idsi != Idsi)
    {
      cout<<setprecision(10)<<id<<"\t"<<Idps[crupos[id]]<<"\t cp="<<cp[id]<<"\t cs="<<cs[crupos[id]]<<endl;
      bSTOP = true;
    }
#endif

#ifdef ___NCX
    double dcp=Betap*(Ir-Ica+jnaca-Idps[crupos[id]]);
#else
    double dcp=Betap*(Ir-Ica-Idps[crupos[id]]);
#endif
    cp[id]+=dcp*dt;
#else
    cp[id]=newcp;
#endif
    cjsr[id]+=dcjsr*dt;

  }



  double sumjup=0;
  double sumjleak=0;
  double sumjnaca=0;
  double sumjcabk=0;
  double sumjslcap=0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+: sumjup, sumjleak,sumjnaca,sumjcabk,sumjslcap)
#endif
#pragma ivdep
#pragma vector always
  for (int id=0;id<nn;id++)
  {
    //SERCA Uptake current Iup
    const double H=1.787;
    double Iup=vup*(pow(ci[id]/kup,H)-pow(cnsr[id]/KNSR,H))/(1+pow(ci[id]/kup,H)+pow(cnsr[id]/KNSR,H));

    //Leak current Ileak
    const double KJSR=500;
    double cjsr2=cnsr[id]*cnsr[id];
    double Ileak=gleak*(cjsr2*(cnsr[id]-ci[id]))/(cjsr2+KJSR*KJSR);

    //NCX
    const double Kmcao=1.3;
    const double Kmnao=87.5;
    const double Kda=0.11;

    double Ka=1/(1+(Kda*Kda*Kda/(cs[id]*cs[id]*cs[id])));
    double t2=Kmnao*Kmnao*Kmnao*(cs[id]*0.001)*(1+cs[id]/Kmcai);
    double t3=Kmcao*nai*nai*nai+nai*nai*nai*cao+nao*nao*nao*(cs[id]*0.001);
#ifdef ___NCX
    double jnaca=(1-NCXalpha)*Ka*x3*(x1a-x1b*(cs[id]*0.001))/((t1+t2+t3)*x2);
#else
    double jnaca=Ka*x3*(x1a-x1b*(cs[id]*0.001))/((t1+t2+t3)*x2);
#endif

    // -------Icabk (background SL Ca flux) following Shannon --------------
    double eca=rtf2*log(cao*1000/cs[id]);
    double icabk=gcabk*(v-eca);// current [A/F]  negative current
    //Icabk[A/F]*(Cm[pF]/(65*27*11*fine3))/(z*F[C/mol]*(vs*10-9/fine3)[ul])=Icabk*3.3286
    double jcabk=icabk*3.3286;//coming in
    // -------Islcap (SL Ca pump) following Shannon --------------
    const double vmax=2.2*0.01;
    const double km=0.5;
    const double h=1.6;
    double islcap=qslcap*vmax/(1+pow(km/cs[id],h));// current [A/F] positive current
    double jslcap=islcap*3.3286;//going out



    //Instantaneous buffering functions
    const double KCAM=7.0;
    const double BCAM=24.0;
    const double KSR=0.6;
    const double BSR=47.0;
    const double KMCa=0.033;
    const double BMCa=140.0;
    const double KMMg=3.64;
    const double BMMg=140.0;
    const double KSLH=0.3;
    const double BSLH=13.4;

    //Troponin C dynamic buffering current ITCi and ITCs
    const double BT=70.0;
    const double kon=0.0327;
    const double koff=0.0196;

    double CAM=BCAM*KCAM/((ci[id]+KCAM)*(ci[id]+KCAM));
    double SR=BSR*KSR/((ci[id]+KSR)*(ci[id]+KSR));
    double MCa=BMCa*KMCa/((ci[id]+KMCa)*(ci[id]+KMCa));
    double MMg=BMMg*KMMg/((ci[id]+KMMg)*(ci[id]+KMMg));
    double Betai=1/(1+CAM+SR+MCa+MMg);
    double ITCi=kon*ci[id]*(BT-cati[id])-koff*cati[id];

#ifdef ___EGTA
    const double konEGTA=4E-3;
    const double koffEGTA=2E-3;
    double IEGTAi=konEGTA*ci[id]*(BEGTA-caEGTAi[id])-koffEGTA*caEGTAi[id];
    double IEGTAs=konEGTA*cs[id]*(BEGTA-caEGTAs[id])-koffEGTA*caEGTAs[id];
#endif



    //      double ITCs=kon*cs[id]*(BT-cats[id])-koff*cats[id];
    double ITCs=0;

    //Diffusion from submembrane to myoplasm Idsi
    double Idsi=(cs[id]-ci[id])/tausi;
#ifdef ___DEBUG
    if(isnan(Idsi)) //(Idsi != Idsi)
    {
      cout<<setprecision(10)<<id<<"\t"<<Idsi<<"\t cs="<<cs[id]<<"\t ci="<<ci[id]<<endl;
      bSTOP = true;
    }
#endif


#ifdef ___NO_CS_BUFFER
#ifdef ___NO_DIFFUSION
    double newcs=(cscp1[id]+jnaca+ci[id]/tausi-ITCs+csmn[id]-jcabk-jslcap)/(cscp2[id]+1/tausi);
#else
    double newcs=(cscp1[id]+jnaca+ci[id]/tausi-ITCs+csmn[id]-jcabk-jslcap)/(cscp2[id]+1/tausi+taumninv);
#endif
    if (newcs<=0)newcs=0.00001; // to avoid div 0
#else
    CAM=BCAM*KCAM/((cs[id]+KCAM)*(cs[id]+KCAM));
    SR=BSR*KSR/((cs[id]+KSR)*(cs[id]+KSR));
    MCa=BMCa*KMCa/((cs[id]+KMCa)*(cs[id]+KMCa));
    MMg=BMMg*KMMg/((cs[id]+KMMg)*(cs[id]+KMMg));
    double SLH=BSLH*KSLH/((cs[id]+KSLH)*(cs[id]+KSLH)); // only for cs
    double Betas=1/(1+CAM+SR+MCa+MMg+SLH);
#ifdef ___EGTA
    double dcs=Betas*(Idps[id]*vp[id]/vs+jnaca-Idsi-ITCs-IEGTAs+Ics[id]-jcabk-jslcap);
#else
    double dcs=Betas*(Idps[id]*vp[id]/vs+jnaca-Idsi-ITCs+Ics[id]-jcabk-jslcap);
#endif
#endif


#ifdef ___EGTA
    double dci=Betai*(Idsi*(vs/vi)-Iup+Ileak-ITCi-IEGTAi+Ici[id]);
#else
    double dci=Betai*(Idsi*(vs/vi)-Iup+Ileak-ITCi+Ici[id]);
#endif
    double dcnsr=((Iup-Ileak)*(vi/vnsr)-Itr[id]*(vjsr/vnsr)+Icnsr[id]);



    ci[id]+=dci*dt;
    cati[id]+=ITCi*dt;
#ifdef ___NO_CS_BUFFER
    cs[id]=newcs;
#else
    cs[id]+=dcs*dt;
#endif
    //      cats[id]+=ITCs*dt;
    cnsr[id]+=dcnsr*dt;

#ifdef ___EGTA
    caEGTAi[id]+=IEGTAi*dt;
    caEGTAs[id]+=IEGTAs*dt;

#endif


    sumjup+=Iup;
    sumjleak+=Ileak;
    sumjnaca+=jnaca;
    sumjcabk+=icabk;
    sumjslcap+=islcap;
  }

  icaave=sumica/n;
  irave=sumir/n;
  incxave=sumjnaca/nn;
  iupave=sumjup/nn;
  ileakave=sumjleak/nn;
  icabkave=sumjcabk/nn;
  islcapave=sumjslcap/nn;
}
int CSubcell::bino(double num, double p, int ii)
{
  int res;
  double lambda=num*p;
  if (lambda>12)
  {
    //Gaussian
    double x1,x2,w;
    do
    {
      x1=2.0*xorshift(&xsx[ii],&xsy[ii],&xsz[ii],&xsw[ii])/(double)(UINT_MAX)-1.0;
      x2=2.0*xorshift(&xsx[ii],&xsy[ii],&xsz[ii],&xsw[ii])/(double)(UINT_MAX)-1.0;
      w=x1*x1+x2*x2;
    }while(w >= 1.0);
    w = sqrt((-2.0*log(w))/w);
    double y1=x1*w;
    //double y2=x2*w;
    res=y1*sqrt(num*p*(1-p))+num*p;// *** ave=num*p , rho^2=num*p*(1-p)
    res=int(res+0.5);//round
  }
  else    if (100*p<6.6+52*pow(num,double(-0.5)))
  {
    //Poisson
    double L=exp(-lambda);
    double k=0;
    double pp=1;
    do
    {
      k++;
      double u=xorshift(&xsx[ii],&xsy[ii],&xsz[ii],&xsw[ii])/(double)(UINT_MAX);
      pp*=u;
    }while(pp>=L);
    res=k-1;
  }
  else
  {
    //Gaussian
    double x1,x2,w;
    do
    {
      x1=2.0*xorshift(&xsx[ii],&xsy[ii],&xsz[ii],&xsw[ii])/(double)(UINT_MAX)-1.0;
      x2=2.0*xorshift(&xsx[ii],&xsy[ii],&xsz[ii],&xsw[ii])/(double)(UINT_MAX)-1.0;
      w=x1*x1+x2*x2;
    }while(w >= 1.0);
    w = sqrt((-2.0*log(w))/w);
    double y1=x1*w;
    //double y2=x2*w;
    res=y1*sqrt(num*p*(1-p))+num*p;// *** ave=num*p , rho^2=num*p*(1-p)
    res=int(res+0.5);//round
  }
  if (res<0)res=0;

  return res;
}

double CSubcell::calcvp(double mean,double std,double lim1, double lim2, int ii)
{
  double res;
  do
  {
    //Gaussian
    double x1,x2,w;
    do
    {
      x1=2.0*xorshift(&xsx[ii],&xsy[ii],&xsz[ii],&xsw[ii])/(double)(UINT_MAX)-1.0;
      x2=2.0*xorshift(&xsx[ii],&xsy[ii],&xsz[ii],&xsw[ii])/(double)(UINT_MAX)-1.0;
      w=x1*x1+x2*x2;
    }while(w >= 1.0);
    w = sqrt((-2.0*log(w))/w );
    double y1 = x1 * w;
    //double y2 = x2 * w;
    res=y1*std+mean;
  }while (res<lim1 || res>lim2);
  return res;
}

double CSubcell::computeaveci(void)
{
  double sum=0;
  for (int id=0;id<nn;id++)
    sum+=ci[id];
  return (sum/nn);
}
double CSubcell::computeavecs(void)
{
  double sum=0;
  for (int id=0;id<nn;id++)
    sum+=cs[id];
  return (sum/nn);
}
double CSubcell::computeavecnsr(void)
{
  double sum=0;
  for (int id=0;id<nn;id++)
    sum+=cnsr[id];
  return (sum/nn);
}

void CSubcell::setboundary(int bcc)
{
  if (bc>0)//corner
  {
    Jmaxx[0+0*nx+0*nxny]=0;
    Jmaxx[0+(ny-1)*nx+(nz-1)*nxny]=0;
    Jmaxx[0+(ny-1)*nx+0*nxny]=0;
    Jmaxx[0+0*nx+(nz-1)*nxny]=0;
    Jmaxx[(nx-1)+(ny-1)*nx+0*nxny]=0;
    Jmaxx[(nx-1)+0*nx+(nz-1)*nxny]=0;
    Jmaxx[(nx-1)+0*nx+0*nxny]=0;
    Jmaxx[(nx-1)+(ny-1)*nx+(nz-1)*nxny]=0;
  }
  if (bc>1)//edge
  {
    //x fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int j=1;j<ny-1;j++)
    {
      Jmaxx[0+j*nx+0*nxny]=0;
      Jmaxx[(nx-1)+j*nx+0*nxny]=0;
      Jmaxx[0+j*nx+(nz-1)*nxny]=0;
      Jmaxx[(nx-1)+j*nx+(nz-1)*nxny]=0;
    }
    //y fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i=1;i<(nx-1);i++)
    {
      Jmaxx[i+0*nx+0*nxny]=0;
      Jmaxx[i+(ny-1)*nx+0*nxny]=0;
      Jmaxx[i+0*nx+(nz-1)*nxny]=0;
      Jmaxx[i+(ny-1)*nx+(nz-1)*nxny]=0;
    }
#pragma ivdep
#pragma vector always
    for (int k=1;k<nz-1;k++)
    {
      Jmaxx[0+0*nx+k*nxny]=0;
      Jmaxx[0+(ny-1)*nx+k*nxny]=0;
      Jmaxx[(nx-1)+0*nx+k*nxny]=0;
      Jmaxx[(nx-1)+(ny-1)*nx+k*nxny]=0;
    }
  }
  if (bc>2)//surface
  {
    //x fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int j=1;j<ny-1;j++)
    {
#pragma ivdep
#pragma vector always
      for (int k=1;k<nz-1;k++)
      {
        Jmaxx[0+j*nx+k*nxny]=0;
        Jmaxx[(nx-1)+j*nx+k*nxny]=0;
      }
    }
    //y fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i=1;i<(nx-1);i++)
    {
#pragma ivdep
#pragma vector always
      for (int k=1;k<nz-1;k++)
      {
        Jmaxx[i+0*nx+k*nxny]=0;
        Jmaxx[i+(ny-1)*nx+k*nxny]=0;
      }
      //z fixed
#pragma ivdep
#pragma vector always
      for (int j=1;j<ny-1;j++)
      {
        Jmaxx[i+j*nx+0*nxny]=0;
        Jmaxx[i+j*nx+(nz-1)*nxny]=0;
      }
    }
  }
  if (bc>3)//corner more
  {
    Jmaxx[1+1*nx+1*nxny]=0;
    Jmaxx[1+(ny-2)*nx+(nz-2)*nxny]=0;
    Jmaxx[1+(ny-2)*nx+1*nxny]=0;
    Jmaxx[1+1*nx+(nz-2)*nxny]=0;
    Jmaxx[(nx-2)+(ny-2)*nx+1*nxny]=0;
    Jmaxx[(nx-2)+1*nx+(nz-2)*nxny]=0;
    Jmaxx[(nx-2)+1*nx+1*nxny]=0;
    Jmaxx[(nx-2)+(ny-2)*nx+(nz-2)*nxny]=0;
  }


}
void CSubcell::setJmax(double newJmax){
  Jmax=newJmax;
  for(int id=0;id<n;id++)Jmaxx[id]=Jmax;
  setboundary(bc);
}
void CSubcell::resetBuffer(void){
  const double BT=70.0;
  const double kon=0.0327;
  const double koff=0.0196;

  const double konEGTA=4E-3;
  const double koffEGTA=2E-3;

#pragma ivdep
#pragma vector always
  for (int id=0;id<nn;id++)
  {
    cati[id]=kon*ci[id]*BT/(kon*ci[id]+koff);
    //      cats[id]=kon*cs[id]*BT/(kon*cs[id]+koff);
#ifdef ___EGTA
    caEGTAi[id]=konEGTA*ci[id]*BEGTA/(konEGTA*ci[id]+koffEGTA);
    caEGTAs[id]=konEGTA*cs[id]*BEGTA/(konEGTA*cs[id]+koffEGTA);
#endif
  }
}

void CSubcell::computeIci(void)
{

#ifdef ___PERIODIC
  Ici[0+0*nnx+0*nnxnny]=(ci[0+0*nnx+0*nnxnny+1]+ci[(nnx-1)+0*nnx+0*nnxnny]-2*ci[0+0*nnx+0*nnxnny])/tauiL+
    (ci[0+(0+1)*nnx+0*nnxnny]+ci[0+(nny-1)*nnx+0*nnxnny]-2*ci[0+0*nnx+0*nnxnny])/tauiT+
    (ci[0+0*nnx+(0+1)*nnxnny]+ci[0+0*nnx+(nnz-1)*nnxnny]-2*ci[0+0*nnx+0*nnxnny])/tauiT;
  Ici[0+(nny-1)*nnx+(nnz-1)*nnxnny]=(ci[0+(nny-1)*nnx+(nnz-1)*nnxnny+1]+ci[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny]-2*ci[0+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiL+
    (ci[0+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+ci[0+(0)*nnx+(nnz-1)*nnxnny]-2*ci[0+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiT+
    (ci[0+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+ci[0+(nny-1)*nnx+(0)*nnxnny]-2*ci[0+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiT;
  Ici[0+(nny-1)*nnx+0*nnxnny]=(ci[0+(nny-1)*nnx+0*nnxnny+1]+ci[(nnx-1)+(nny-1)*nnx+0*nnxnny]-2*ci[0+(nny-1)*nnx+0*nnxnny])/tauiL+
    (ci[0+((nny-1)-1)*nnx+0*nnxnny]+ci[0+(0)*nnx+0*nnxnny]-2*ci[0+(nny-1)*nnx+0*nnxnny])/tauiT+
    (ci[0+(nny-1)*nnx+(0+1)*nnxnny]+ci[0+(nny-1)*nnx+(nnz-1)*nnxnny]-2*ci[0+(nny-1)*nnx+0*nnxnny])/tauiT;
  Ici[0+0*nnx+(nnz-1)*nnxnny]=(ci[0+0*nnx+(nnz-1)*nnxnny+1]+ci[(nnx-1)+0*nnx+(nnz-1)*nnxnny]-2*ci[0+0*nnx+(nnz-1)*nnxnny])/tauiL+
    (ci[0+(0+1)*nnx+(nnz-1)*nnxnny]+ci[0+(nny-1)*nnx+(nnz-1)*nnxnny]-2*ci[0+0*nnx+(nnz-1)*nnxnny])/tauiT+
    (ci[0+0*nnx+((nnz-1)-1)*nnxnny]+ci[0+0*nnx+(0)*nnxnny]-2*ci[0+0*nnx+(nnz-1)*nnxnny])/tauiT;
  Ici[(nnx-1)+(nny-1)*nnx+0*nnxnny]=(ci[(nnx-1)+(nny-1)*nnx+0*nnxnny-1]+ci[(0)+(nny-1)*nnx+0*nnxnny]-2*ci[(nnx-1)+(nny-1)*nnx+0*nnxnny])/tauiL+
    (ci[(nnx-1)+((nny-1)-1)*nnx+0*nnxnny]+ci[(nnx-1)+(0)*nnx+0*nnxnny]-2*ci[(nnx-1)+(nny-1)*nnx+0*nnxnny])/tauiT+
    (ci[(nnx-1)+(nny-1)*nnx+(0+1)*nnxnny]+ci[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny]-2*ci[(nnx-1)+(nny-1)*nnx+0*nnxnny])/tauiT;
  Ici[(nnx-1)+0*nnx+(nnz-1)*nnxnny]=(ci[(nnx-1)+0*nnx+(nnz-1)*nnxnny-1]+ci[(0)+0*nnx+(nnz-1)*nnxnny]-2*ci[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/tauiL+
    (ci[(nnx-1)+(0+1)*nnx+(nnz-1)*nnxnny]+ci[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny]-2*ci[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/tauiT+
    (ci[(nnx-1)+0*nnx+((nnz-1)-1)*nnxnny]+ci[(nnx-1)+0*nnx+(0)*nnxnny]-2*ci[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/tauiT;
  Ici[(nnx-1)+0*nnx+0*nnxnny]=(ci[(nnx-1)+0*nnx+0*nnxnny-1]+ci[(0)+0*nnx+0*nnxnny]-2*ci[(nnx-1)+0*nnx+0*nnxnny])/tauiL+
    (ci[(nnx-1)+(0+1)*nnx+0*nnxnny]+ci[(nnx-1)+(nny-1)*nnx+0*nnxnny]-2*ci[(nnx-1)+0*nnx+0*nnxnny])/tauiT+
    (ci[(nnx-1)+0*nnx+(0+1)*nnxnny]+ci[(nnx-1)+0*nnx+(nnz-1)*nnxnny]-2*ci[(nnx-1)+0*nnx+0*nnxnny])/tauiT;
  Ici[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny]=(ci[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny-1]+ci[(0)+(nny-1)*nnx+(nnz-1)*nnxnny]-2*ci[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiL+
    (ci[(nnx-1)+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+ci[(nnx-1)+(0)*nnx+(nnz-1)*nnxnny]-2*ci[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiT+
    (ci[(nnx-1)+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+ci[(nnx-1)+(nny-1)*nnx+(0)*nnxnny]-2*ci[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiT;
  //x fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int j=1;j<nny-1;j++)
  {
    Ici[0+j*nnx+0*nnxnny]=(ci[0+j*nnx+0*nnxnny+1]+ci[(nnx-1)+j*nnx+0*nnxnny]-2*ci[0+j*nnx+0*nnxnny])/tauiL+
      (ci[0+(j+1)*nnx+0*nnxnny]+ci[0+(j-1)*nnx+0*nnxnny]-2*ci[0+j*nnx+0*nnxnny])/tauiT+
      (ci[0+j*nnx+(0+1)*nnxnny]+ci[0+j*nnx+(nnz-1)*nnxnny]-2*ci[0+j*nnx+0*nnxnny])/tauiT;
    Ici[(nnx-1)+j*nnx+0*nnxnny]=(ci[(nnx-1)+j*nnx+0*nnxnny-1]+ci[(0)+j*nnx+0*nnxnny]-2*ci[(nnx-1)+j*nnx+0*nnxnny])/tauiL+
      (ci[(nnx-1)+(j+1)*nnx+0*nnxnny]+ci[(nnx-1)+(j-1)*nnx+0*nnxnny]-2*ci[(nnx-1)+j*nnx+0*nnxnny])/tauiT+
      (ci[(nnx-1)+j*nnx+(0+1)*nnxnny]+ci[(nnx-1)+j*nnx+(nnz-1)*nnxnny]-2*ci[(nnx-1)+j*nnx+0*nnxnny])/tauiT;
    Ici[0+j*nnx+(nnz-1)*nnxnny]=(ci[0+j*nnx+(nnz-1)*nnxnny+1]+ci[(nnx-1)+j*nnx+(nnz-1)*nnxnny]-2*ci[0+j*nnx+(nnz-1)*nnxnny])/tauiL+
      (ci[0+(j+1)*nnx+(nnz-1)*nnxnny]+ci[0+(j-1)*nnx+(nnz-1)*nnxnny]-2*ci[0+j*nnx+(nnz-1)*nnxnny])/tauiT+
      (ci[0+j*nnx+((nnz-1)-1)*nnxnny]+ci[0+j*nnx+(0)*nnxnny]-2*ci[0+j*nnx+(nnz-1)*nnxnny])/tauiT;
    Ici[(nnx-1)+j*nnx+(nnz-1)*nnxnny]=(ci[(nnx-1)+j*nnx+(nnz-1)*nnxnny-1]+ci[(0)+j*nnx+(nnz-1)*nnxnny]-2*ci[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/tauiL+
      (ci[(nnx-1)+(j+1)*nnx+(nnz-1)*nnxnny]+ci[(nnx-1)+(j-1)*nnx+(nnz-1)*nnxnny]-2*ci[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/tauiT+
      (ci[(nnx-1)+j*nnx+((nnz-1)-1)*nnxnny]+ci[(nnx-1)+j*nnx+(0)*nnxnny]-2*ci[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/tauiT;
#pragma ivdep
#pragma vector always
    for (int k=1;k<nnz-1;k++)
    {
      Ici[0+j*nnx+k*nnxnny]=(ci[0+j*nnx+k*nnxnny+1]+ci[(nnx-1)+j*nnx+k*nnxnny]-2*ci[0+j*nnx+k*nnxnny])/tauiL+
        (ci[0+(j+1)*nnx+k*nnxnny]+ci[0+(j-1)*nnx+k*nnxnny]-2*ci[0+j*nnx+k*nnxnny])/tauiT+
        (ci[0+j*nnx+(k+1)*nnxnny]+ci[0+j*nnx+(k-1)*nnxnny]-2*ci[0+j*nnx+k*nnxnny])/tauiT;
      Ici[(nnx-1)+j*nnx+k*nnxnny]=(ci[(nnx-1)+j*nnx+k*nnxnny-1]+ci[(0)+j*nnx+k*nnxnny]-2*ci[(nnx-1)+j*nnx+k*nnxnny])/tauiL+
        (ci[(nnx-1)+(j+1)*nnx+k*nnxnny]+ci[(nnx-1)+(j-1)*nnx+k*nnxnny]-2*ci[(nnx-1)+j*nnx+k*nnxnny])/tauiT+
        (ci[(nnx-1)+j*nnx+(k+1)*nnxnny]+ci[(nnx-1)+j*nnx+(k-1)*nnxnny]-2*ci[(nnx-1)+j*nnx+k*nnxnny])/tauiT;
    }
  }
  //y fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=1;i<(nnx-1);i++)
  {
    Ici[i+0*nnx+0*nnxnny]=(ci[i+0*nnx+0*nnxnny+1]+ci[i+0*nnx+0*nnxnny-1]-2*ci[i+0*nnx+0*nnxnny])/tauiL+
      (ci[i+(0+1)*nnx+0*nnxnny]+ci[i+(nny-1)*nnx+0*nnxnny]-2*ci[i+0*nnx+0*nnxnny])/tauiT+
      (ci[i+0*nnx+(0+1)*nnxnny]+ci[i+0*nnx+(nnz-1)*nnxnny]-2*ci[i+0*nnx+0*nnxnny])/tauiT;
    Ici[i+(nny-1)*nnx+0*nnxnny]=(ci[i+(nny-1)*nnx+0*nnxnny+1]+ci[i+(nny-1)*nnx+0*nnxnny-1]-2*ci[i+(nny-1)*nnx+0*nnxnny])/tauiL+
      (ci[i+((nny-1)-1)*nnx+0*nnxnny]+ci[i+(0)*nnx+0*nnxnny]-2*ci[i+(nny-1)*nnx+0*nnxnny])/tauiT+
      (ci[i+(nny-1)*nnx+(0+1)*nnxnny]+ci[i+(nny-1)*nnx+(nnz-1)*nnxnny]-2*ci[i+(nny-1)*nnx+0*nnxnny])/tauiT;
    Ici[i+0*nnx+(nnz-1)*nnxnny]=(ci[i+0*nnx+(nnz-1)*nnxnny+1]+ci[i+0*nnx+(nnz-1)*nnxnny-1]-2*ci[i+0*nnx+(nnz-1)*nnxnny])/tauiL+
      (ci[i+(0+1)*nnx+(nnz-1)*nnxnny]+ci[i+(nny-1)*nnx+(nnz-1)*nnxnny]-2*ci[i+0*nnx+(nnz-1)*nnxnny])/tauiT+
      (ci[i+0*nnx+((nnz-1)-1)*nnxnny]+ci[i+0*nnx+(0)*nnxnny]-2*ci[i+0*nnx+(nnz-1)*nnxnny])/tauiT;
    Ici[i+(nny-1)*nnx+(nnz-1)*nnxnny]=(ci[i+(nny-1)*nnx+(nnz-1)*nnxnny+1]+ci[i+(nny-1)*nnx+(nnz-1)*nnxnny-1]-2*ci[i+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiL+
      (ci[i+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+ci[i+(0)*nnx+(nnz-1)*nnxnny]-2*ci[i+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiT+
      (ci[i+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+ci[i+(nny-1)*nnx+(0)*nnxnny]-2*ci[i+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiT;
#pragma ivdep
#pragma vector always
    for (int k=1;k<nnz-1;k++)
    {
      Ici[i+0*nnx+k*nnxnny]=(ci[i+0*nnx+k*nnxnny+1]+ci[i+0*nnx+k*nnxnny-1]-2*ci[i+0*nnx+k*nnxnny])/tauiL+
        (ci[i+(0+1)*nnx+k*nnxnny]+ci[i+(nny-1)*nnx+k*nnxnny]-2*ci[i+0*nnx+k*nnxnny])/tauiT+
        (ci[i+0*nnx+(k+1)*nnxnny]+ci[i+0*nnx+(k-1)*nnxnny]-2*ci[i+0*nnx+k*nnxnny])/tauiT;
      Ici[i+(nny-1)*nnx+k*nnxnny]=(ci[i+(nny-1)*nnx+k*nnxnny+1]+ci[i+(nny-1)*nnx+k*nnxnny-1]-2*ci[i+(nny-1)*nnx+k*nnxnny])/tauiL+
        (ci[i+((nny-1)-1)*nnx+k*nnxnny]+ci[i+(0)*nnx+k*nnxnny]-2*ci[i+(nny-1)*nnx+k*nnxnny])/tauiT+
        (ci[i+(nny-1)*nnx+(k+1)*nnxnny]+ci[i+(nny-1)*nnx+(k-1)*nnxnny]-2*ci[i+(nny-1)*nnx+k*nnxnny])/tauiT;
    }
    //z fixed
#pragma ivdep
#pragma vector always
    for (int j=1;j<nny-1;j++)
    {
      Ici[i+j*nnx+0*nnxnny]=(ci[i+j*nnx+0*nnxnny+1]+ci[i+j*nnx+0*nnxnny-1]-2*ci[i+j*nnx+0*nnxnny])/tauiL+
        (ci[i+(j+1)*nnx+0*nnxnny]+ci[i+(j-1)*nnx+0*nnxnny]-2*ci[i+j*nnx+0*nnxnny])/tauiT+
        (ci[i+j*nnx+(0+1)*nnxnny]+ci[i+j*nnx+(nnz-1)*nnxnny]-2*ci[i+j*nnx+0*nnxnny])/tauiT;
      Ici[i+j*nnx+(nnz-1)*nnxnny]=(ci[i+j*nnx+(nnz-1)*nnxnny+1]+ci[i+j*nnx+(nnz-1)*nnxnny-1]-2*ci[i+j*nnx+(nnz-1)*nnxnny])/tauiL+
        (ci[i+(j+1)*nnx+(nnz-1)*nnxnny]+ci[i+(j-1)*nnx+(nnz-1)*nnxnny]-2*ci[i+j*nnx+(nnz-1)*nnxnny])/tauiT+
        (ci[i+j*nnx+((nnz-1)-1)*nnxnny]+ci[i+j*nnx+(0)*nnxnny]-2*ci[i+j*nnx+(nnz-1)*nnxnny])/tauiT;
    }
  }
#pragma ivdep
#pragma vector always
  for (int k=1;k<nnz-1;k++)
  {
    Ici[0+0*nnx+k*nnxnny]=(ci[0+0*nnx+k*nnxnny+1]+ci[(nnx-1)+0*nnx+k*nnxnny]-2*ci[0+0*nnx+k*nnxnny])/tauiL+
      (ci[0+(0+1)*nnx+k*nnxnny]+ci[0+(nny-1)*nnx+k*nnxnny]-2*ci[0+0*nnx+k*nnxnny])/tauiT+
      (ci[0+0*nnx+(k+1)*nnxnny]+ci[0+0*nnx+(k-1)*nnxnny]-2*ci[0+0*nnx+k*nnxnny])/tauiT;
    Ici[0+(nny-1)*nnx+k*nnxnny]=(ci[0+(nny-1)*nnx+k*nnxnny+1]+ci[(nnx-1)+(nny-1)*nnx+k*nnxnny]-2*ci[0+(nny-1)*nnx+k*nnxnny])/tauiL+
      (ci[0+((nny-1)-1)*nnx+k*nnxnny]+ci[0+(0)*nnx+k*nnxnny]-2*ci[0+(nny-1)*nnx+k*nnxnny])/tauiT+
      (ci[0+(nny-1)*nnx+(k+1)*nnxnny]+ci[0+(nny-1)*nnx+(k-1)*nnxnny]-2*ci[0+(nny-1)*nnx+k*nnxnny])/tauiT;
    Ici[(nnx-1)+0*nnx+k*nnxnny]=(ci[(nnx-1)+0*nnx+k*nnxnny-1]+ci[(0)+0*nnx+k*nnxnny]-2*ci[(nnx-1)+0*nnx+k*nnxnny])/tauiL+
      (ci[(nnx-1)+(0+1)*nnx+k*nnxnny]+ci[(nnx-1)+(nny-1)*nnx+k*nnxnny]-2*ci[(nnx-1)+0*nnx+k*nnxnny])/tauiT+
      (ci[(nnx-1)+0*nnx+(k+1)*nnxnny]+ci[(nnx-1)+0*nnx+(k-1)*nnxnny]-2*ci[(nnx-1)+0*nnx+k*nnxnny])/tauiT;
    Ici[(nnx-1)+(nny-1)*nnx+k*nnxnny]=(ci[(nnx-1)+(nny-1)*nnx+k*nnxnny-1]+ci[(0)+(nny-1)*nnx+k*nnxnny]-2*ci[(nnx-1)+(nny-1)*nnx+k*nnxnny])/tauiL+
      (ci[(nnx-1)+((nny-1)-1)*nnx+k*nnxnny]+ci[(nnx-1)+(0)*nnx+k*nnxnny]-2*ci[(nnx-1)+(nny-1)*nnx+k*nnxnny])/tauiT+
      (ci[(nnx-1)+(nny-1)*nnx+(k+1)*nnxnny]+ci[(nnx-1)+(nny-1)*nnx+(k-1)*nnxnny]-2*ci[(nnx-1)+(nny-1)*nnx+k*nnxnny])/tauiT;
  }
#else
  Ici[0+0*nnx+0*nnxnny]=(ci[0+0*nnx+0*nnxnny+1]+ci[0+0*nnx+0*nnxnny+1]-2*ci[0+0*nnx+0*nnxnny])/tauiL+
    (ci[0+(0+1)*nnx+0*nnxnny]+ci[0+(0+1)*nnx+0*nnxnny]-2*ci[0+0*nnx+0*nnxnny])/tauiT+
    (ci[0+0*nnx+(0+1)*nnxnny]+ci[0+0*nnx+(0+1)*nnxnny]-2*ci[0+0*nnx+0*nnxnny])/tauiT;
  Ici[0+(nny-1)*nnx+(nnz-1)*nnxnny]=(ci[0+(nny-1)*nnx+(nnz-1)*nnxnny+1]+ci[0+(nny-1)*nnx+(nnz-1)*nnxnny+1]-2*ci[0+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiL+
    (ci[0+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+ci[0+((nny-1)-1)*nnx+(nnz-1)*nnxnny]-2*ci[0+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiT+
    (ci[0+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+ci[0+(nny-1)*nnx+((nnz-1)-1)*nnxnny]-2*ci[0+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiT;
  Ici[0+(nny-1)*nnx+0*nnxnny]=(ci[0+(nny-1)*nnx+0*nnxnny+1]+ci[0+(nny-1)*nnx+0*nnxnny+1]-2*ci[0+(nny-1)*nnx+0*nnxnny])/tauiL+
    (ci[0+((nny-1)-1)*nnx+0*nnxnny]+ci[0+((nny-1)-1)*nnx+0*nnxnny]-2*ci[0+(nny-1)*nnx+0*nnxnny])/tauiT+
    (ci[0+(nny-1)*nnx+(0+1)*nnxnny]+ci[0+(nny-1)*nnx+(0+1)*nnxnny]-2*ci[0+(nny-1)*nnx+0*nnxnny])/tauiT;
  Ici[0+0*nnx+(nnz-1)*nnxnny]=(ci[0+0*nnx+(nnz-1)*nnxnny+1]+ci[0+0*nnx+(nnz-1)*nnxnny+1]-2*ci[0+0*nnx+(nnz-1)*nnxnny])/tauiL+
    (ci[0+(0+1)*nnx+(nnz-1)*nnxnny]+ci[0+(0+1)*nnx+(nnz-1)*nnxnny]-2*ci[0+0*nnx+(nnz-1)*nnxnny])/tauiT+
    (ci[0+0*nnx+((nnz-1)-1)*nnxnny]+ci[0+0*nnx+((nnz-1)-1)*nnxnny]-2*ci[0+0*nnx+(nnz-1)*nnxnny])/tauiT;
  Ici[(nnx-1)+(nny-1)*nnx+0*nnxnny]=(ci[(nnx-1)+(nny-1)*nnx+0*nnxnny-1]+ci[(nnx-1)+(nny-1)*nnx+0*nnxnny-1]-2*ci[(nnx-1)+(nny-1)*nnx+0*nnxnny])/tauiL+
    (ci[(nnx-1)+((nny-1)-1)*nnx+0*nnxnny]+ci[(nnx-1)+((nny-1)-1)*nnx+0*nnxnny]-2*ci[(nnx-1)+(nny-1)*nnx+0*nnxnny])/tauiT+
    (ci[(nnx-1)+(nny-1)*nnx+(0+1)*nnxnny]+ci[(nnx-1)+(nny-1)*nnx+(0+1)*nnxnny]-2*ci[(nnx-1)+(nny-1)*nnx+0*nnxnny])/tauiT;
  Ici[(nnx-1)+0*nnx+(nnz-1)*nnxnny]=(ci[(nnx-1)+0*nnx+(nnz-1)*nnxnny-1]+ci[(nnx-1)+0*nnx+(nnz-1)*nnxnny-1]-2*ci[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/tauiL+
    (ci[(nnx-1)+(0+1)*nnx+(nnz-1)*nnxnny]+ci[(nnx-1)+(0+1)*nnx+(nnz-1)*nnxnny]-2*ci[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/tauiT+
    (ci[(nnx-1)+0*nnx+((nnz-1)-1)*nnxnny]+ci[(nnx-1)+0*nnx+((nnz-1)-1)*nnxnny]-2*ci[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/tauiT;
  Ici[(nnx-1)+0*nnx+0*nnxnny]=(ci[(nnx-1)+0*nnx+0*nnxnny-1]+ci[(nnx-1)+0*nnx+0*nnxnny-1]-2*ci[(nnx-1)+0*nnx+0*nnxnny])/tauiL+
    (ci[(nnx-1)+(0+1)*nnx+0*nnxnny]+ci[(nnx-1)+(0+1)*nnx+0*nnxnny]-2*ci[(nnx-1)+0*nnx+0*nnxnny])/tauiT+
    (ci[(nnx-1)+0*nnx+(0+1)*nnxnny]+ci[(nnx-1)+0*nnx+(0+1)*nnxnny]-2*ci[(nnx-1)+0*nnx+0*nnxnny])/tauiT;
  Ici[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny]=(ci[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny-1]+ci[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny-1]-2*ci[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiL+
    (ci[(nnx-1)+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+ci[(nnx-1)+((nny-1)-1)*nnx+(nnz-1)*nnxnny]-2*ci[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiT+
    (ci[(nnx-1)+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+ci[(nnx-1)+(nny-1)*nnx+((nnz-1)-1)*nnxnny]-2*ci[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiT;
  //x fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int j=1;j<nny-1;j++)
  {
    Ici[0+j*nnx+0*nnxnny]=(ci[0+j*nnx+0*nnxnny+1]+ci[0+j*nnx+0*nnxnny+1]-2*ci[0+j*nnx+0*nnxnny])/tauiL+
      (ci[0+(j+1)*nnx+0*nnxnny]+ci[0+(j-1)*nnx+0*nnxnny]-2*ci[0+j*nnx+0*nnxnny])/tauiT+
      (ci[0+j*nnx+(0+1)*nnxnny]+ci[0+j*nnx+(0+1)*nnxnny]-2*ci[0+j*nnx+0*nnxnny])/tauiT;
    Ici[(nnx-1)+j*nnx+0*nnxnny]=(ci[(nnx-1)+j*nnx+0*nnxnny-1]+ci[(nnx-1)+j*nnx+0*nnxnny-1]-2*ci[(nnx-1)+j*nnx+0*nnxnny])/tauiL+
      (ci[(nnx-1)+(j+1)*nnx+0*nnxnny]+ci[(nnx-1)+(j-1)*nnx+0*nnxnny]-2*ci[(nnx-1)+j*nnx+0*nnxnny])/tauiT+
      (ci[(nnx-1)+j*nnx+(0+1)*nnxnny]+ci[(nnx-1)+j*nnx+(0+1)*nnxnny]-2*ci[(nnx-1)+j*nnx+0*nnxnny])/tauiT;
    Ici[0+j*nnx+(nnz-1)*nnxnny]=(ci[0+j*nnx+(nnz-1)*nnxnny+1]+ci[0+j*nnx+(nnz-1)*nnxnny+1]-2*ci[0+j*nnx+(nnz-1)*nnxnny])/tauiL+
      (ci[0+(j+1)*nnx+(nnz-1)*nnxnny]+ci[0+(j-1)*nnx+(nnz-1)*nnxnny]-2*ci[0+j*nnx+(nnz-1)*nnxnny])/tauiT+
      (ci[0+j*nnx+((nnz-1)-1)*nnxnny]+ci[0+j*nnx+((nnz-1)-1)*nnxnny]-2*ci[0+j*nnx+(nnz-1)*nnxnny])/tauiT;
    Ici[(nnx-1)+j*nnx+(nnz-1)*nnxnny]=(ci[(nnx-1)+j*nnx+(nnz-1)*nnxnny-1]+ci[(nnx-1)+j*nnx+(nnz-1)*nnxnny-1]-2*ci[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/tauiL+
      (ci[(nnx-1)+(j+1)*nnx+(nnz-1)*nnxnny]+ci[(nnx-1)+(j-1)*nnx+(nnz-1)*nnxnny]-2*ci[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/tauiT+
      (ci[(nnx-1)+j*nnx+((nnz-1)-1)*nnxnny]+ci[(nnx-1)+j*nnx+((nnz-1)-1)*nnxnny]-2*ci[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/tauiT;
#pragma ivdep
#pragma vector always
    for (int k=1;k<nnz-1;k++)
    {
      Ici[0+j*nnx+k*nnxnny]=(ci[0+j*nnx+k*nnxnny+1]+ci[0+j*nnx+k*nnxnny+1]-2*ci[0+j*nnx+k*nnxnny])/tauiL+
        (ci[0+(j+1)*nnx+k*nnxnny]+ci[0+(j-1)*nnx+k*nnxnny]-2*ci[0+j*nnx+k*nnxnny])/tauiT+
        (ci[0+j*nnx+(k+1)*nnxnny]+ci[0+j*nnx+(k-1)*nnxnny]-2*ci[0+j*nnx+k*nnxnny])/tauiT;
      Ici[(nnx-1)+j*nnx+k*nnxnny]=(ci[(nnx-1)+j*nnx+k*nnxnny-1]+ci[(nnx-1)+j*nnx+k*nnxnny-1]-2*ci[(nnx-1)+j*nnx+k*nnxnny])/tauiL+
        (ci[(nnx-1)+(j+1)*nnx+k*nnxnny]+ci[(nnx-1)+(j-1)*nnx+k*nnxnny]-2*ci[(nnx-1)+j*nnx+k*nnxnny])/tauiT+
        (ci[(nnx-1)+j*nnx+(k+1)*nnxnny]+ci[(nnx-1)+j*nnx+(k-1)*nnxnny]-2*ci[(nnx-1)+j*nnx+k*nnxnny])/tauiT;
    }
  }
  //y fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=1;i<(nnx-1);i++)
  {
    Ici[i+0*nnx+0*nnxnny]=(ci[i+0*nnx+0*nnxnny+1]+ci[i+0*nnx+0*nnxnny-1]-2*ci[i+0*nnx+0*nnxnny])/tauiL+
      (ci[i+(0+1)*nnx+0*nnxnny]+ci[i+(0+1)*nnx+0*nnxnny]-2*ci[i+0*nnx+0*nnxnny])/tauiT+
      (ci[i+0*nnx+(0+1)*nnxnny]+ci[i+0*nnx+(0+1)*nnxnny]-2*ci[i+0*nnx+0*nnxnny])/tauiT;
    Ici[i+(nny-1)*nnx+0*nnxnny]=(ci[i+(nny-1)*nnx+0*nnxnny+1]+ci[i+(nny-1)*nnx+0*nnxnny-1]-2*ci[i+(nny-1)*nnx+0*nnxnny])/tauiL+
      (ci[i+((nny-1)-1)*nnx+0*nnxnny]+ci[i+((nny-1)-1)*nnx+0*nnxnny]-2*ci[i+(nny-1)*nnx+0*nnxnny])/tauiT+
      (ci[i+(nny-1)*nnx+(0+1)*nnxnny]+ci[i+(nny-1)*nnx+(0+1)*nnxnny]-2*ci[i+(nny-1)*nnx+0*nnxnny])/tauiT;
    Ici[i+0*nnx+(nnz-1)*nnxnny]=(ci[i+0*nnx+(nnz-1)*nnxnny+1]+ci[i+0*nnx+(nnz-1)*nnxnny-1]-2*ci[i+0*nnx+(nnz-1)*nnxnny])/tauiL+
      (ci[i+(0+1)*nnx+(nnz-1)*nnxnny]+ci[i+(0+1)*nnx+(nnz-1)*nnxnny]-2*ci[i+0*nnx+(nnz-1)*nnxnny])/tauiT+
      (ci[i+0*nnx+((nnz-1)-1)*nnxnny]+ci[i+0*nnx+((nnz-1)-1)*nnxnny]-2*ci[i+0*nnx+(nnz-1)*nnxnny])/tauiT;
    Ici[i+(nny-1)*nnx+(nnz-1)*nnxnny]=(ci[i+(nny-1)*nnx+(nnz-1)*nnxnny+1]+ci[i+(nny-1)*nnx+(nnz-1)*nnxnny-1]-2*ci[i+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiL+
      (ci[i+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+ci[i+((nny-1)-1)*nnx+(nnz-1)*nnxnny]-2*ci[i+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiT+
      (ci[i+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+ci[i+(nny-1)*nnx+((nnz-1)-1)*nnxnny]-2*ci[i+(nny-1)*nnx+(nnz-1)*nnxnny])/tauiT;
#pragma ivdep
#pragma vector always
    for (int k=1;k<nnz-1;k++)
    {
      Ici[i+0*nnx+k*nnxnny]=(ci[i+0*nnx+k*nnxnny+1]+ci[i+0*nnx+k*nnxnny-1]-2*ci[i+0*nnx+k*nnxnny])/tauiL+
        (ci[i+(0+1)*nnx+k*nnxnny]+ci[i+(0+1)*nnx+k*nnxnny]-2*ci[i+0*nnx+k*nnxnny])/tauiT+
        (ci[i+0*nnx+(k+1)*nnxnny]+ci[i+0*nnx+(k-1)*nnxnny]-2*ci[i+0*nnx+k*nnxnny])/tauiT;
      Ici[i+(nny-1)*nnx+k*nnxnny]=(ci[i+(nny-1)*nnx+k*nnxnny+1]+ci[i+(nny-1)*nnx+k*nnxnny-1]-2*ci[i+(nny-1)*nnx+k*nnxnny])/tauiL+
        (ci[i+((nny-1)-1)*nnx+k*nnxnny]+ci[i+((nny-1)-1)*nnx+k*nnxnny]-2*ci[i+(nny-1)*nnx+k*nnxnny])/tauiT+
        (ci[i+(nny-1)*nnx+(k+1)*nnxnny]+ci[i+(nny-1)*nnx+(k-1)*nnxnny]-2*ci[i+(nny-1)*nnx+k*nnxnny])/tauiT;
    }
    //z fixed
#pragma ivdep
#pragma vector always
    for (int j=1;j<nny-1;j++)
    {
      Ici[i+j*nnx+0*nnxnny]=(ci[i+j*nnx+0*nnxnny+1]+ci[i+j*nnx+0*nnxnny-1]-2*ci[i+j*nnx+0*nnxnny])/tauiL+
        (ci[i+(j+1)*nnx+0*nnxnny]+ci[i+(j-1)*nnx+0*nnxnny]-2*ci[i+j*nnx+0*nnxnny])/tauiT+
        (ci[i+j*nnx+(0+1)*nnxnny]+ci[i+j*nnx+(0+1)*nnxnny]-2*ci[i+j*nnx+0*nnxnny])/tauiT;
      Ici[i+j*nnx+(nnz-1)*nnxnny]=(ci[i+j*nnx+(nnz-1)*nnxnny+1]+ci[i+j*nnx+(nnz-1)*nnxnny-1]-2*ci[i+j*nnx+(nnz-1)*nnxnny])/tauiL+
        (ci[i+(j+1)*nnx+(nnz-1)*nnxnny]+ci[i+(j-1)*nnx+(nnz-1)*nnxnny]-2*ci[i+j*nnx+(nnz-1)*nnxnny])/tauiT+
        (ci[i+j*nnx+((nnz-1)-1)*nnxnny]+ci[i+j*nnx+((nnz-1)-1)*nnxnny]-2*ci[i+j*nnx+(nnz-1)*nnxnny])/tauiT;
    }
  }
#pragma ivdep
#pragma vector always
  for (int k=1;k<nnz-1;k++)
  {
    Ici[0+0*nnx+k*nnxnny]=(ci[0+0*nnx+k*nnxnny+1]+ci[0+0*nnx+k*nnxnny+1]-2*ci[0+0*nnx+k*nnxnny])/tauiL+
      (ci[0+(0+1)*nnx+k*nnxnny]+ci[0+(0+1)*nnx+k*nnxnny]-2*ci[0+0*nnx+k*nnxnny])/tauiT+
      (ci[0+0*nnx+(k+1)*nnxnny]+ci[0+0*nnx+(k-1)*nnxnny]-2*ci[0+0*nnx+k*nnxnny])/tauiT;
    Ici[0+(nny-1)*nnx+k*nnxnny]=(ci[0+(nny-1)*nnx+k*nnxnny+1]+ci[0+(nny-1)*nnx+k*nnxnny+1]-2*ci[0+(nny-1)*nnx+k*nnxnny])/tauiL+
      (ci[0+((nny-1)-1)*nnx+k*nnxnny]+ci[0+((nny-1)-1)*nnx+k*nnxnny]-2*ci[0+(nny-1)*nnx+k*nnxnny])/tauiT+
      (ci[0+(nny-1)*nnx+(k+1)*nnxnny]+ci[0+(nny-1)*nnx+(k-1)*nnxnny]-2*ci[0+(nny-1)*nnx+k*nnxnny])/tauiT;
    Ici[(nnx-1)+0*nnx+k*nnxnny]=(ci[(nnx-1)+0*nnx+k*nnxnny-1]+ci[(nnx-1)+0*nnx+k*nnxnny-1]-2*ci[(nnx-1)+0*nnx+k*nnxnny])/tauiL+
      (ci[(nnx-1)+(0+1)*nnx+k*nnxnny]+ci[(nnx-1)+(0+1)*nnx+k*nnxnny]-2*ci[(nnx-1)+0*nnx+k*nnxnny])/tauiT+
      (ci[(nnx-1)+0*nnx+(k+1)*nnxnny]+ci[(nnx-1)+0*nnx+(k-1)*nnxnny]-2*ci[(nnx-1)+0*nnx+k*nnxnny])/tauiT;
    Ici[(nnx-1)+(nny-1)*nnx+k*nnxnny]=(ci[(nnx-1)+(nny-1)*nnx+k*nnxnny-1]+ci[(nnx-1)+(nny-1)*nnx+k*nnxnny-1]-2*ci[(nnx-1)+(nny-1)*nnx+k*nnxnny])/tauiL+
      (ci[(nnx-1)+((nny-1)-1)*nnx+k*nnxnny]+ci[(nnx-1)+((nny-1)-1)*nnx+k*nnxnny]-2*ci[(nnx-1)+(nny-1)*nnx+k*nnxnny])/tauiT+
      (ci[(nnx-1)+(nny-1)*nnx+(k+1)*nnxnny]+ci[(nnx-1)+(nny-1)*nnx+(k-1)*nnxnny]-2*ci[(nnx-1)+(nny-1)*nnx+k*nnxnny])/tauiT;
  }

#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int k=1;k<nnz-1;k++)
  {
    for (int j=1;j<nny-1;j++)
    {
#pragma ivdep
#pragma vector always
      for (int i=1;i<nnx-1;i++)
      {
        Ici[i+j*nnx+k*nnxnny]=(ci[i+j*nnx+k*nnxnny+1]+ci[i+j*nnx+k*nnxnny-1]-2*ci[i+j*nnx+k*nnxnny])/tauiL+
          (ci[i+(j+1)*nnx+k*nnxnny]+ci[i+(j-1)*nnx+k*nnxnny]-2*ci[i+j*nnx+k*nnxnny])/tauiT+
          (ci[i+j*nnx+(k+1)*nnxnny]+ci[i+j*nnx+(k-1)*nnxnny]-2*ci[i+j*nnx+k*nnxnny])/tauiT;
      }
    }
  }
}
void CSubcell::computeIcnsr(void)
{
#ifdef ___PERIODIC
  //boundary
  Icnsr[0+0*nnx+0*nnxnny]=(cnsr[0+0*nnx+0*nnxnny+1]+cnsr[(nnx-1)+0*nnx+0*nnxnny]-2*cnsr[0+0*nnx+0*nnxnny])/taunsrL+
    (cnsr[0+(0+1)*nnx+0*nnxnny]+cnsr[0+(nny-1)*nnx+0*nnxnny]-2*cnsr[0+0*nnx+0*nnxnny])/taunsrT+
    (cnsr[0+0*nnx+(0+1)*nnxnny]+cnsr[0+0*nnx+(nnz-1)*nnxnny]-2*cnsr[0+0*nnx+0*nnxnny])/taunsrT;
  Icnsr[0+(nny-1)*nnx+(nnz-1)*nnxnny]=(cnsr[0+(nny-1)*nnx+(nnz-1)*nnxnny+1]+cnsr[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny]-2*cnsr[0+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrL+
    (cnsr[0+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cnsr[0+(0)*nnx+(nnz-1)*nnxnny]-2*cnsr[0+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrT+
    (cnsr[0+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cnsr[0+(nny-1)*nnx+(0)*nnxnny]-2*cnsr[0+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrT;
  Icnsr[0+(nny-1)*nnx+0*nnxnny]=(cnsr[0+(nny-1)*nnx+0*nnxnny+1]+cnsr[(nnx-1)+(nny-1)*nnx+0*nnxnny]-2*cnsr[0+(nny-1)*nnx+0*nnxnny])/taunsrL+
    (cnsr[0+((nny-1)-1)*nnx+0*nnxnny]+cnsr[0+(0)*nnx+0*nnxnny]-2*cnsr[0+(nny-1)*nnx+0*nnxnny])/taunsrT+
    (cnsr[0+(nny-1)*nnx+(0+1)*nnxnny]+cnsr[0+(nny-1)*nnx+(nnz-1)*nnxnny]-2*cnsr[0+(nny-1)*nnx+0*nnxnny])/taunsrT;
  Icnsr[0+0*nnx+(nnz-1)*nnxnny]=(cnsr[0+0*nnx+(nnz-1)*nnxnny+1]+cnsr[(nnx-1)+0*nnx+(nnz-1)*nnxnny]-2*cnsr[0+0*nnx+(nnz-1)*nnxnny])/taunsrL+
    (cnsr[0+(0+1)*nnx+(nnz-1)*nnxnny]+cnsr[0+(nny-1)*nnx+(nnz-1)*nnxnny]-2*cnsr[0+0*nnx+(nnz-1)*nnxnny])/taunsrT+
    (cnsr[0+0*nnx+((nnz-1)-1)*nnxnny]+cnsr[0+0*nnx+(0)*nnxnny]-2*cnsr[0+0*nnx+(nnz-1)*nnxnny])/taunsrT;
  Icnsr[(nnx-1)+(nny-1)*nnx+0*nnxnny]=(cnsr[(nnx-1)+(nny-1)*nnx+0*nnxnny-1]+cnsr[(0)+(nny-1)*nnx+0*nnxnny]-2*cnsr[(nnx-1)+(nny-1)*nnx+0*nnxnny])/taunsrL+
    (cnsr[(nnx-1)+((nny-1)-1)*nnx+0*nnxnny]+cnsr[(nnx-1)+(0)*nnx+0*nnxnny]-2*cnsr[(nnx-1)+(nny-1)*nnx+0*nnxnny])/taunsrT+
    (cnsr[(nnx-1)+(nny-1)*nnx+(0+1)*nnxnny]+cnsr[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny]-2*cnsr[(nnx-1)+(nny-1)*nnx+0*nnxnny])/taunsrT;
  Icnsr[(nnx-1)+0*nnx+(nnz-1)*nnxnny]=(cnsr[(nnx-1)+0*nnx+(nnz-1)*nnxnny-1]+cnsr[(0)+0*nnx+(nnz-1)*nnxnny]-2*cnsr[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/taunsrL+
    (cnsr[(nnx-1)+(0+1)*nnx+(nnz-1)*nnxnny]+cnsr[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny]-2*cnsr[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/taunsrT+
    (cnsr[(nnx-1)+0*nnx+((nnz-1)-1)*nnxnny]+cnsr[(nnx-1)+0*nnx+(0)*nnxnny]-2*cnsr[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/taunsrT;
  Icnsr[(nnx-1)+0*nnx+0*nnxnny]=(cnsr[(nnx-1)+0*nnx+0*nnxnny-1]+cnsr[(0)+0*nnx+0*nnxnny]-2*cnsr[(nnx-1)+0*nnx+0*nnxnny])/taunsrL+
    (cnsr[(nnx-1)+(0+1)*nnx+0*nnxnny]+cnsr[(nnx-1)+(nny-1)*nnx+0*nnxnny]-2*cnsr[(nnx-1)+0*nnx+0*nnxnny])/taunsrT+
    (cnsr[(nnx-1)+0*nnx+(0+1)*nnxnny]+cnsr[(nnx-1)+0*nnx+(nnz-1)*nnxnny]-2*cnsr[(nnx-1)+0*nnx+0*nnxnny])/taunsrT;
  Icnsr[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny]=(cnsr[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny-1]+cnsr[(0)+(nny-1)*nnx+(nnz-1)*nnxnny]-2*cnsr[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrL+
    (cnsr[(nnx-1)+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cnsr[(nnx-1)+(0)*nnx+(nnz-1)*nnxnny]-2*cnsr[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrT+
    (cnsr[(nnx-1)+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cnsr[(nnx-1)+(nny-1)*nnx+(0)*nnxnny]-2*cnsr[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrT;
  //x fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int j=1;j<nny-1;j++)
  {
    Icnsr[0+j*nnx+0*nnxnny]=(cnsr[0+j*nnx+0*nnxnny+1]+cnsr[(nnx-1)+j*nnx+0*nnxnny]-2*cnsr[0+j*nnx+0*nnxnny])/taunsrL+
      (cnsr[0+(j+1)*nnx+0*nnxnny]+cnsr[0+(j-1)*nnx+0*nnxnny]-2*cnsr[0+j*nnx+0*nnxnny])/taunsrT+
      (cnsr[0+j*nnx+(0+1)*nnxnny]+cnsr[0+j*nnx+(nnz-1)*nnxnny]-2*cnsr[0+j*nnx+0*nnxnny])/taunsrT;
    Icnsr[(nnx-1)+j*nnx+0*nnxnny]=(cnsr[(nnx-1)+j*nnx+0*nnxnny-1]+cnsr[(0)+j*nnx+0*nnxnny]-2*cnsr[(nnx-1)+j*nnx+0*nnxnny])/taunsrL+
      (cnsr[(nnx-1)+(j+1)*nnx+0*nnxnny]+cnsr[(nnx-1)+(j-1)*nnx+0*nnxnny]-2*cnsr[(nnx-1)+j*nnx+0*nnxnny])/taunsrT+
      (cnsr[(nnx-1)+j*nnx+(0+1)*nnxnny]+cnsr[(nnx-1)+j*nnx+(nnz-1)*nnxnny]-2*cnsr[(nnx-1)+j*nnx+0*nnxnny])/taunsrT;
    Icnsr[0+j*nnx+(nnz-1)*nnxnny]=(cnsr[0+j*nnx+(nnz-1)*nnxnny+1]+cnsr[(nnx-1)+j*nnx+(nnz-1)*nnxnny]-2*cnsr[0+j*nnx+(nnz-1)*nnxnny])/taunsrL+
      (cnsr[0+(j+1)*nnx+(nnz-1)*nnxnny]+cnsr[0+(j-1)*nnx+(nnz-1)*nnxnny]-2*cnsr[0+j*nnx+(nnz-1)*nnxnny])/taunsrT+
      (cnsr[0+j*nnx+((nnz-1)-1)*nnxnny]+cnsr[0+j*nnx+(0)*nnxnny]-2*cnsr[0+j*nnx+(nnz-1)*nnxnny])/taunsrT;
    Icnsr[(nnx-1)+j*nnx+(nnz-1)*nnxnny]=(cnsr[(nnx-1)+j*nnx+(nnz-1)*nnxnny-1]+cnsr[(0)+j*nnx+(nnz-1)*nnxnny]-2*cnsr[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/taunsrL+
      (cnsr[(nnx-1)+(j+1)*nnx+(nnz-1)*nnxnny]+cnsr[(nnx-1)+(j-1)*nnx+(nnz-1)*nnxnny]-2*cnsr[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/taunsrT+
      (cnsr[(nnx-1)+j*nnx+((nnz-1)-1)*nnxnny]+cnsr[(nnx-1)+j*nnx+(0)*nnxnny]-2*cnsr[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/taunsrT;
#pragma ivdep
#pragma vector always
    for (int k=1;k<nnz-1;k++)
    {
      Icnsr[0+j*nnx+k*nnxnny]=(cnsr[0+j*nnx+k*nnxnny+1]+cnsr[(nnx-1)+j*nnx+k*nnxnny]-2*cnsr[0+j*nnx+k*nnxnny])/taunsrL+
        (cnsr[0+(j+1)*nnx+k*nnxnny]+cnsr[0+(j-1)*nnx+k*nnxnny]-2*cnsr[0+j*nnx+k*nnxnny])/taunsrT+
        (cnsr[0+j*nnx+(k+1)*nnxnny]+cnsr[0+j*nnx+(k-1)*nnxnny]-2*cnsr[0+j*nnx+k*nnxnny])/taunsrT;
      Icnsr[(nnx-1)+j*nnx+k*nnxnny]=(cnsr[(nnx-1)+j*nnx+k*nnxnny-1]+cnsr[(0)+j*nnx+k*nnxnny]-2*cnsr[(nnx-1)+j*nnx+k*nnxnny])/taunsrL+
        (cnsr[(nnx-1)+(j+1)*nnx+k*nnxnny]+cnsr[(nnx-1)+(j-1)*nnx+k*nnxnny]-2*cnsr[(nnx-1)+j*nnx+k*nnxnny])/taunsrT+
        (cnsr[(nnx-1)+j*nnx+(k+1)*nnxnny]+cnsr[(nnx-1)+j*nnx+(k-1)*nnxnny]-2*cnsr[(nnx-1)+j*nnx+k*nnxnny])/taunsrT;
    }
  }
  //y fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=1;i<(nnx-1);i++)
  {
    Icnsr[i+0*nnx+0*nnxnny]=(cnsr[i+0*nnx+0*nnxnny+1]+cnsr[i+0*nnx+0*nnxnny-1]-2*cnsr[i+0*nnx+0*nnxnny])/taunsrL+
      (cnsr[i+(0+1)*nnx+0*nnxnny]+cnsr[i+(nny-1)*nnx+0*nnxnny]-2*cnsr[i+0*nnx+0*nnxnny])/taunsrT+
      (cnsr[i+0*nnx+(0+1)*nnxnny]+cnsr[i+0*nnx+(nnz-1)*nnxnny]-2*cnsr[i+0*nnx+0*nnxnny])/taunsrT;
    Icnsr[i+(nny-1)*nnx+0*nnxnny]=(cnsr[i+(nny-1)*nnx+0*nnxnny+1]+cnsr[i+(nny-1)*nnx+0*nnxnny-1]-2*cnsr[i+(nny-1)*nnx+0*nnxnny])/taunsrL+
      (cnsr[i+((nny-1)-1)*nnx+0*nnxnny]+cnsr[i+(0)*nnx+0*nnxnny]-2*cnsr[i+(nny-1)*nnx+0*nnxnny])/taunsrT+
      (cnsr[i+(nny-1)*nnx+(0+1)*nnxnny]+cnsr[i+(nny-1)*nnx+(nnz-1)*nnxnny]-2*cnsr[i+(nny-1)*nnx+0*nnxnny])/taunsrT;
    Icnsr[i+0*nnx+(nnz-1)*nnxnny]=(cnsr[i+0*nnx+(nnz-1)*nnxnny+1]+cnsr[i+0*nnx+(nnz-1)*nnxnny-1]-2*cnsr[i+0*nnx+(nnz-1)*nnxnny])/taunsrL+
      (cnsr[i+(0+1)*nnx+(nnz-1)*nnxnny]+cnsr[i+(nny-1)*nnx+(nnz-1)*nnxnny]-2*cnsr[i+0*nnx+(nnz-1)*nnxnny])/taunsrT+
      (cnsr[i+0*nnx+((nnz-1)-1)*nnxnny]+cnsr[i+0*nnx+(0)*nnxnny]-2*cnsr[i+0*nnx+(nnz-1)*nnxnny])/taunsrT;
    Icnsr[i+(nny-1)*nnx+(nnz-1)*nnxnny]=(cnsr[i+(nny-1)*nnx+(nnz-1)*nnxnny+1]+cnsr[i+(nny-1)*nnx+(nnz-1)*nnxnny-1]-2*cnsr[i+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrL+
      (cnsr[i+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cnsr[i+(0)*nnx+(nnz-1)*nnxnny]-2*cnsr[i+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrT+
      (cnsr[i+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cnsr[i+(nny-1)*nnx+(0)*nnxnny]-2*cnsr[i+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrT;
#pragma ivdep
#pragma vector always
    for (int k=1;k<nnz-1;k++)
    {
      Icnsr[i+0*nnx+k*nnxnny]=(cnsr[i+0*nnx+k*nnxnny+1]+cnsr[i+0*nnx+k*nnxnny-1]-2*cnsr[i+0*nnx+k*nnxnny])/taunsrL+
        (cnsr[i+(0+1)*nnx+k*nnxnny]+cnsr[i+(nny-1)*nnx+k*nnxnny]-2*cnsr[i+0*nnx+k*nnxnny])/taunsrT+
        (cnsr[i+0*nnx+(k+1)*nnxnny]+cnsr[i+0*nnx+(k-1)*nnxnny]-2*cnsr[i+0*nnx+k*nnxnny])/taunsrT;
      Icnsr[i+(nny-1)*nnx+k*nnxnny]=(cnsr[i+(nny-1)*nnx+k*nnxnny+1]+cnsr[i+(nny-1)*nnx+k*nnxnny-1]-2*cnsr[i+(nny-1)*nnx+k*nnxnny])/taunsrL+
        (cnsr[i+((nny-1)-1)*nnx+k*nnxnny]+cnsr[i+(0)*nnx+k*nnxnny]-2*cnsr[i+(nny-1)*nnx+k*nnxnny])/taunsrT+
        (cnsr[i+(nny-1)*nnx+(k+1)*nnxnny]+cnsr[i+(nny-1)*nnx+(k-1)*nnxnny]-2*cnsr[i+(nny-1)*nnx+k*nnxnny])/taunsrT;
    }
    //z fixed
#pragma ivdep
#pragma vector always
    for (int j=1;j<nny-1;j++)
    {
      Icnsr[i+j*nnx+0*nnxnny]=(cnsr[i+j*nnx+0*nnxnny+1]+cnsr[i+j*nnx+0*nnxnny-1]-2*cnsr[i+j*nnx+0*nnxnny])/taunsrL+
        (cnsr[i+(j+1)*nnx+0*nnxnny]+cnsr[i+(j-1)*nnx+0*nnxnny]-2*cnsr[i+j*nnx+0*nnxnny])/taunsrT+
        (cnsr[i+j*nnx+(0+1)*nnxnny]+cnsr[i+j*nnx+(nnz-1)*nnxnny]-2*cnsr[i+j*nnx+0*nnxnny])/taunsrT;
      Icnsr[i+j*nnx+(nnz-1)*nnxnny]=(cnsr[i+j*nnx+(nnz-1)*nnxnny+1]+cnsr[i+j*nnx+(nnz-1)*nnxnny-1]-2*cnsr[i+j*nnx+(nnz-1)*nnxnny])/taunsrL+
        (cnsr[i+(j+1)*nnx+(nnz-1)*nnxnny]+cnsr[i+(j-1)*nnx+(nnz-1)*nnxnny]-2*cnsr[i+j*nnx+(nnz-1)*nnxnny])/taunsrT+
        (cnsr[i+j*nnx+((nnz-1)-1)*nnxnny]+cnsr[i+j*nnx+(0)*nnxnny]-2*cnsr[i+j*nnx+(nnz-1)*nnxnny])/taunsrT;
    }
  }
#pragma ivdep
#pragma vector always
  for (int k=1;k<nnz-1;k++)
  {
    Icnsr[0+0*nnx+k*nnxnny]=(cnsr[0+0*nnx+k*nnxnny+1]+cnsr[(nnx-1)+0*nnx+k*nnxnny]-2*cnsr[0+0*nnx+k*nnxnny])/taunsrL+
      (cnsr[0+(0+1)*nnx+k*nnxnny]+cnsr[0+(nny-1)*nnx+k*nnxnny]-2*cnsr[0+0*nnx+k*nnxnny])/taunsrT+
      (cnsr[0+0*nnx+(k+1)*nnxnny]+cnsr[0+0*nnx+(k-1)*nnxnny]-2*cnsr[0+0*nnx+k*nnxnny])/taunsrT;
    Icnsr[0+(nny-1)*nnx+k*nnxnny]=(cnsr[0+(nny-1)*nnx+k*nnxnny+1]+cnsr[(nnx-1)+(nny-1)*nnx+k*nnxnny]-2*cnsr[0+(nny-1)*nnx+k*nnxnny])/taunsrL+
      (cnsr[0+((nny-1)-1)*nnx+k*nnxnny]+cnsr[0+(0)*nnx+k*nnxnny]-2*cnsr[0+(nny-1)*nnx+k*nnxnny])/taunsrT+
      (cnsr[0+(nny-1)*nnx+(k+1)*nnxnny]+cnsr[0+(nny-1)*nnx+(k-1)*nnxnny]-2*cnsr[0+(nny-1)*nnx+k*nnxnny])/taunsrT;
    Icnsr[(nnx-1)+0*nnx+k*nnxnny]=(cnsr[(nnx-1)+0*nnx+k*nnxnny-1]+cnsr[(0)+0*nnx+k*nnxnny]-2*cnsr[(nnx-1)+0*nnx+k*nnxnny])/taunsrL+
      (cnsr[(nnx-1)+(0+1)*nnx+k*nnxnny]+cnsr[(nnx-1)+(nny-1)*nnx+k*nnxnny]-2*cnsr[(nnx-1)+0*nnx+k*nnxnny])/taunsrT+
      (cnsr[(nnx-1)+0*nnx+(k+1)*nnxnny]+cnsr[(nnx-1)+0*nnx+(k-1)*nnxnny]-2*cnsr[(nnx-1)+0*nnx+k*nnxnny])/taunsrT;
    Icnsr[(nnx-1)+(nny-1)*nnx+k*nnxnny]=(cnsr[(nnx-1)+(nny-1)*nnx+k*nnxnny-1]+cnsr[(0)+(nny-1)*nnx+k*nnxnny]-2*cnsr[(nnx-1)+(nny-1)*nnx+k*nnxnny])/taunsrL+
      (cnsr[(nnx-1)+((nny-1)-1)*nnx+k*nnxnny]+cnsr[(nnx-1)+(0)*nnx+k*nnxnny]-2*cnsr[(nnx-1)+(nny-1)*nnx+k*nnxnny])/taunsrT+
      (cnsr[(nnx-1)+(nny-1)*nnx+(k+1)*nnxnny]+cnsr[(nnx-1)+(nny-1)*nnx+(k-1)*nnxnny]-2*cnsr[(nnx-1)+(nny-1)*nnx+k*nnxnny])/taunsrT;
  }
#else
  //boundary

  Icnsr[0+0*nnx+0*nnxnny]=(cnsr[0+0*nnx+0*nnxnny+1]+cnsr[0+0*nnx+0*nnxnny+1]-2*cnsr[0+0*nnx+0*nnxnny])/taunsrL+
    (cnsr[0+(0+1)*nnx+0*nnxnny]+cnsr[0+(0+1)*nnx+0*nnxnny]-2*cnsr[0+0*nnx+0*nnxnny])/taunsrT+
    (cnsr[0+0*nnx+(0+1)*nnxnny]+cnsr[0+0*nnx+(0+1)*nnxnny]-2*cnsr[0+0*nnx+0*nnxnny])/taunsrT;
  Icnsr[0+(nny-1)*nnx+(nnz-1)*nnxnny]=(cnsr[0+(nny-1)*nnx+(nnz-1)*nnxnny+1]+cnsr[0+(nny-1)*nnx+(nnz-1)*nnxnny+1]-2*cnsr[0+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrL+
    (cnsr[0+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cnsr[0+((nny-1)-1)*nnx+(nnz-1)*nnxnny]-2*cnsr[0+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrT+
    (cnsr[0+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cnsr[0+(nny-1)*nnx+((nnz-1)-1)*nnxnny]-2*cnsr[0+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrT;
  Icnsr[0+(nny-1)*nnx+0*nnxnny]=(cnsr[0+(nny-1)*nnx+0*nnxnny+1]+cnsr[0+(nny-1)*nnx+0*nnxnny+1]-2*cnsr[0+(nny-1)*nnx+0*nnxnny])/taunsrL+
    (cnsr[0+((nny-1)-1)*nnx+0*nnxnny]+cnsr[0+((nny-1)-1)*nnx+0*nnxnny]-2*cnsr[0+(nny-1)*nnx+0*nnxnny])/taunsrT+
    (cnsr[0+(nny-1)*nnx+(0+1)*nnxnny]+cnsr[0+(nny-1)*nnx+(0+1)*nnxnny]-2*cnsr[0+(nny-1)*nnx+0*nnxnny])/taunsrT;
  Icnsr[0+0*nnx+(nnz-1)*nnxnny]=(cnsr[0+0*nnx+(nnz-1)*nnxnny+1]+cnsr[0+0*nnx+(nnz-1)*nnxnny+1]-2*cnsr[0+0*nnx+(nnz-1)*nnxnny])/taunsrL+
    (cnsr[0+(0+1)*nnx+(nnz-1)*nnxnny]+cnsr[0+(0+1)*nnx+(nnz-1)*nnxnny]-2*cnsr[0+0*nnx+(nnz-1)*nnxnny])/taunsrT+
    (cnsr[0+0*nnx+((nnz-1)-1)*nnxnny]+cnsr[0+0*nnx+((nnz-1)-1)*nnxnny]-2*cnsr[0+0*nnx+(nnz-1)*nnxnny])/taunsrT;
  Icnsr[(nnx-1)+(nny-1)*nnx+0*nnxnny]=(cnsr[(nnx-1)+(nny-1)*nnx+0*nnxnny-1]+cnsr[(nnx-1)+(nny-1)*nnx+0*nnxnny-1]-2*cnsr[(nnx-1)+(nny-1)*nnx+0*nnxnny])/taunsrL+
    (cnsr[(nnx-1)+((nny-1)-1)*nnx+0*nnxnny]+cnsr[(nnx-1)+((nny-1)-1)*nnx+0*nnxnny]-2*cnsr[(nnx-1)+(nny-1)*nnx+0*nnxnny])/taunsrT+
    (cnsr[(nnx-1)+(nny-1)*nnx+(0+1)*nnxnny]+cnsr[(nnx-1)+(nny-1)*nnx+(0+1)*nnxnny]-2*cnsr[(nnx-1)+(nny-1)*nnx+0*nnxnny])/taunsrT;
  Icnsr[(nnx-1)+0*nnx+(nnz-1)*nnxnny]=(cnsr[(nnx-1)+0*nnx+(nnz-1)*nnxnny-1]+cnsr[(nnx-1)+0*nnx+(nnz-1)*nnxnny-1]-2*cnsr[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/taunsrL+
    (cnsr[(nnx-1)+(0+1)*nnx+(nnz-1)*nnxnny]+cnsr[(nnx-1)+(0+1)*nnx+(nnz-1)*nnxnny]-2*cnsr[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/taunsrT+
    (cnsr[(nnx-1)+0*nnx+((nnz-1)-1)*nnxnny]+cnsr[(nnx-1)+0*nnx+((nnz-1)-1)*nnxnny]-2*cnsr[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/taunsrT;
  Icnsr[(nnx-1)+0*nnx+0*nnxnny]=(cnsr[(nnx-1)+0*nnx+0*nnxnny-1]+cnsr[(nnx-1)+0*nnx+0*nnxnny-1]-2*cnsr[(nnx-1)+0*nnx+0*nnxnny])/taunsrL+
    (cnsr[(nnx-1)+(0+1)*nnx+0*nnxnny]+cnsr[(nnx-1)+(0+1)*nnx+0*nnxnny]-2*cnsr[(nnx-1)+0*nnx+0*nnxnny])/taunsrT+
    (cnsr[(nnx-1)+0*nnx+(0+1)*nnxnny]+cnsr[(nnx-1)+0*nnx+(0+1)*nnxnny]-2*cnsr[(nnx-1)+0*nnx+0*nnxnny])/taunsrT;
  Icnsr[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny]=(cnsr[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny-1]+cnsr[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny-1]-2*cnsr[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrL+
    (cnsr[(nnx-1)+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cnsr[(nnx-1)+((nny-1)-1)*nnx+(nnz-1)*nnxnny]-2*cnsr[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrT+
    (cnsr[(nnx-1)+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cnsr[(nnx-1)+(nny-1)*nnx+((nnz-1)-1)*nnxnny]-2*cnsr[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrT;
  //x fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int j=1;j<nny-1;j++)
  {
    Icnsr[0+j*nnx+0*nnxnny]=(cnsr[0+j*nnx+0*nnxnny+1]+cnsr[0+j*nnx+0*nnxnny+1]-2*cnsr[0+j*nnx+0*nnxnny])/taunsrL+
      (cnsr[0+(j+1)*nnx+0*nnxnny]+cnsr[0+(j-1)*nnx+0*nnxnny]-2*cnsr[0+j*nnx+0*nnxnny])/taunsrT+
      (cnsr[0+j*nnx+(0+1)*nnxnny]+cnsr[0+j*nnx+(0+1)*nnxnny]-2*cnsr[0+j*nnx+0*nnxnny])/taunsrT;
    Icnsr[(nnx-1)+j*nnx+0*nnxnny]=(cnsr[(nnx-1)+j*nnx+0*nnxnny-1]+cnsr[(nnx-1)+j*nnx+0*nnxnny-1]-2*cnsr[(nnx-1)+j*nnx+0*nnxnny])/taunsrL+
      (cnsr[(nnx-1)+(j+1)*nnx+0*nnxnny]+cnsr[(nnx-1)+(j-1)*nnx+0*nnxnny]-2*cnsr[(nnx-1)+j*nnx+0*nnxnny])/taunsrT+
      (cnsr[(nnx-1)+j*nnx+(0+1)*nnxnny]+cnsr[(nnx-1)+j*nnx+(0+1)*nnxnny]-2*cnsr[(nnx-1)+j*nnx+0*nnxnny])/taunsrT;
    Icnsr[0+j*nnx+(nnz-1)*nnxnny]=(cnsr[0+j*nnx+(nnz-1)*nnxnny+1]+cnsr[0+j*nnx+(nnz-1)*nnxnny+1]-2*cnsr[0+j*nnx+(nnz-1)*nnxnny])/taunsrL+
      (cnsr[0+(j+1)*nnx+(nnz-1)*nnxnny]+cnsr[0+(j-1)*nnx+(nnz-1)*nnxnny]-2*cnsr[0+j*nnx+(nnz-1)*nnxnny])/taunsrT+
      (cnsr[0+j*nnx+((nnz-1)-1)*nnxnny]+cnsr[0+j*nnx+((nnz-1)-1)*nnxnny]-2*cnsr[0+j*nnx+(nnz-1)*nnxnny])/taunsrT;
    Icnsr[(nnx-1)+j*nnx+(nnz-1)*nnxnny]=(cnsr[(nnx-1)+j*nnx+(nnz-1)*nnxnny-1]+cnsr[(nnx-1)+j*nnx+(nnz-1)*nnxnny-1]-2*cnsr[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/taunsrL+
      (cnsr[(nnx-1)+(j+1)*nnx+(nnz-1)*nnxnny]+cnsr[(nnx-1)+(j-1)*nnx+(nnz-1)*nnxnny]-2*cnsr[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/taunsrT+
      (cnsr[(nnx-1)+j*nnx+((nnz-1)-1)*nnxnny]+cnsr[(nnx-1)+j*nnx+((nnz-1)-1)*nnxnny]-2*cnsr[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/taunsrT;
#pragma ivdep
#pragma vector always
    for (int k=1;k<nnz-1;k++)
    {
      Icnsr[0+j*nnx+k*nnxnny]=(cnsr[0+j*nnx+k*nnxnny+1]+cnsr[0+j*nnx+k*nnxnny+1]-2*cnsr[0+j*nnx+k*nnxnny])/taunsrL+
        (cnsr[0+(j+1)*nnx+k*nnxnny]+cnsr[0+(j-1)*nnx+k*nnxnny]-2*cnsr[0+j*nnx+k*nnxnny])/taunsrT+
        (cnsr[0+j*nnx+(k+1)*nnxnny]+cnsr[0+j*nnx+(k-1)*nnxnny]-2*cnsr[0+j*nnx+k*nnxnny])/taunsrT;
      Icnsr[(nnx-1)+j*nnx+k*nnxnny]=(cnsr[(nnx-1)+j*nnx+k*nnxnny-1]+cnsr[(nnx-1)+j*nnx+k*nnxnny-1]-2*cnsr[(nnx-1)+j*nnx+k*nnxnny])/taunsrL+
        (cnsr[(nnx-1)+(j+1)*nnx+k*nnxnny]+cnsr[(nnx-1)+(j-1)*nnx+k*nnxnny]-2*cnsr[(nnx-1)+j*nnx+k*nnxnny])/taunsrT+
        (cnsr[(nnx-1)+j*nnx+(k+1)*nnxnny]+cnsr[(nnx-1)+j*nnx+(k-1)*nnxnny]-2*cnsr[(nnx-1)+j*nnx+k*nnxnny])/taunsrT;
    }
  }
  //y fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=1;i<(nnx-1);i++)
  {
    Icnsr[i+0*nnx+0*nnxnny]=(cnsr[i+0*nnx+0*nnxnny+1]+cnsr[i+0*nnx+0*nnxnny-1]-2*cnsr[i+0*nnx+0*nnxnny])/taunsrL+
      (cnsr[i+(0+1)*nnx+0*nnxnny]+cnsr[i+(0+1)*nnx+0*nnxnny]-2*cnsr[i+0*nnx+0*nnxnny])/taunsrT+
      (cnsr[i+0*nnx+(0+1)*nnxnny]+cnsr[i+0*nnx+(0+1)*nnxnny]-2*cnsr[i+0*nnx+0*nnxnny])/taunsrT;
    Icnsr[i+(nny-1)*nnx+0*nnxnny]=(cnsr[i+(nny-1)*nnx+0*nnxnny+1]+cnsr[i+(nny-1)*nnx+0*nnxnny-1]-2*cnsr[i+(nny-1)*nnx+0*nnxnny])/taunsrL+
      (cnsr[i+((nny-1)-1)*nnx+0*nnxnny]+cnsr[i+((nny-1)-1)*nnx+0*nnxnny]-2*cnsr[i+(nny-1)*nnx+0*nnxnny])/taunsrT+
      (cnsr[i+(nny-1)*nnx+(0+1)*nnxnny]+cnsr[i+(nny-1)*nnx+(0+1)*nnxnny]-2*cnsr[i+(nny-1)*nnx+0*nnxnny])/taunsrT;
    Icnsr[i+0*nnx+(nnz-1)*nnxnny]=(cnsr[i+0*nnx+(nnz-1)*nnxnny+1]+cnsr[i+0*nnx+(nnz-1)*nnxnny-1]-2*cnsr[i+0*nnx+(nnz-1)*nnxnny])/taunsrL+
      (cnsr[i+(0+1)*nnx+(nnz-1)*nnxnny]+cnsr[i+(0+1)*nnx+(nnz-1)*nnxnny]-2*cnsr[i+0*nnx+(nnz-1)*nnxnny])/taunsrT+
      (cnsr[i+0*nnx+((nnz-1)-1)*nnxnny]+cnsr[i+0*nnx+((nnz-1)-1)*nnxnny]-2*cnsr[i+0*nnx+(nnz-1)*nnxnny])/taunsrT;
    Icnsr[i+(nny-1)*nnx+(nnz-1)*nnxnny]=(cnsr[i+(nny-1)*nnx+(nnz-1)*nnxnny+1]+cnsr[i+(nny-1)*nnx+(nnz-1)*nnxnny-1]-2*cnsr[i+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrL+
      (cnsr[i+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cnsr[i+((nny-1)-1)*nnx+(nnz-1)*nnxnny]-2*cnsr[i+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrT+
      (cnsr[i+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cnsr[i+(nny-1)*nnx+((nnz-1)-1)*nnxnny]-2*cnsr[i+(nny-1)*nnx+(nnz-1)*nnxnny])/taunsrT;
#pragma ivdep
#pragma vector always
    for (int k=1;k<nnz-1;k++)
    {
      Icnsr[i+0*nnx+k*nnxnny]=(cnsr[i+0*nnx+k*nnxnny+1]+cnsr[i+0*nnx+k*nnxnny-1]-2*cnsr[i+0*nnx+k*nnxnny])/taunsrL+
        (cnsr[i+(0+1)*nnx+k*nnxnny]+cnsr[i+(0+1)*nnx+k*nnxnny]-2*cnsr[i+0*nnx+k*nnxnny])/taunsrT+
        (cnsr[i+0*nnx+(k+1)*nnxnny]+cnsr[i+0*nnx+(k-1)*nnxnny]-2*cnsr[i+0*nnx+k*nnxnny])/taunsrT;
      Icnsr[i+(nny-1)*nnx+k*nnxnny]=(cnsr[i+(nny-1)*nnx+k*nnxnny+1]+cnsr[i+(nny-1)*nnx+k*nnxnny-1]-2*cnsr[i+(nny-1)*nnx+k*nnxnny])/taunsrL+
        (cnsr[i+((nny-1)-1)*nnx+k*nnxnny]+cnsr[i+((nny-1)-1)*nnx+k*nnxnny]-2*cnsr[i+(nny-1)*nnx+k*nnxnny])/taunsrT+
        (cnsr[i+(nny-1)*nnx+(k+1)*nnxnny]+cnsr[i+(nny-1)*nnx+(k-1)*nnxnny]-2*cnsr[i+(nny-1)*nnx+k*nnxnny])/taunsrT;
    }
    //z fixed
#pragma ivdep
#pragma vector always
    for (int j=1;j<nny-1;j++)
    {
      Icnsr[i+j*nnx+0*nnxnny]=(cnsr[i+j*nnx+0*nnxnny+1]+cnsr[i+j*nnx+0*nnxnny-1]-2*cnsr[i+j*nnx+0*nnxnny])/taunsrL+
        (cnsr[i+(j+1)*nnx+0*nnxnny]+cnsr[i+(j-1)*nnx+0*nnxnny]-2*cnsr[i+j*nnx+0*nnxnny])/taunsrT+
        (cnsr[i+j*nnx+(0+1)*nnxnny]+cnsr[i+j*nnx+(0+1)*nnxnny]-2*cnsr[i+j*nnx+0*nnxnny])/taunsrT;
      Icnsr[i+j*nnx+(nnz-1)*nnxnny]=(cnsr[i+j*nnx+(nnz-1)*nnxnny+1]+cnsr[i+j*nnx+(nnz-1)*nnxnny-1]-2*cnsr[i+j*nnx+(nnz-1)*nnxnny])/taunsrL+
        (cnsr[i+(j+1)*nnx+(nnz-1)*nnxnny]+cnsr[i+(j-1)*nnx+(nnz-1)*nnxnny]-2*cnsr[i+j*nnx+(nnz-1)*nnxnny])/taunsrT+
        (cnsr[i+j*nnx+((nnz-1)-1)*nnxnny]+cnsr[i+j*nnx+((nnz-1)-1)*nnxnny]-2*cnsr[i+j*nnx+(nnz-1)*nnxnny])/taunsrT;
    }
  }
#pragma ivdep
#pragma vector always
  for (int k=1;k<nnz-1;k++)
  {
    Icnsr[0+0*nnx+k*nnxnny]=(cnsr[0+0*nnx+k*nnxnny+1]+cnsr[0+0*nnx+k*nnxnny+1]-2*cnsr[0+0*nnx+k*nnxnny])/taunsrL+
      (cnsr[0+(0+1)*nnx+k*nnxnny]+cnsr[0+(0+1)*nnx+k*nnxnny]-2*cnsr[0+0*nnx+k*nnxnny])/taunsrT+
      (cnsr[0+0*nnx+(k+1)*nnxnny]+cnsr[0+0*nnx+(k-1)*nnxnny]-2*cnsr[0+0*nnx+k*nnxnny])/taunsrT;
    Icnsr[0+(nny-1)*nnx+k*nnxnny]=(cnsr[0+(nny-1)*nnx+k*nnxnny+1]+cnsr[0+(nny-1)*nnx+k*nnxnny+1]-2*cnsr[0+(nny-1)*nnx+k*nnxnny])/taunsrL+
      (cnsr[0+((nny-1)-1)*nnx+k*nnxnny]+cnsr[0+((nny-1)-1)*nnx+k*nnxnny]-2*cnsr[0+(nny-1)*nnx+k*nnxnny])/taunsrT+
      (cnsr[0+(nny-1)*nnx+(k+1)*nnxnny]+cnsr[0+(nny-1)*nnx+(k-1)*nnxnny]-2*cnsr[0+(nny-1)*nnx+k*nnxnny])/taunsrT;
    Icnsr[(nnx-1)+0*nnx+k*nnxnny]=(cnsr[(nnx-1)+0*nnx+k*nnxnny-1]+cnsr[(nnx-1)+0*nnx+k*nnxnny-1]-2*cnsr[(nnx-1)+0*nnx+k*nnxnny])/taunsrL+
      (cnsr[(nnx-1)+(0+1)*nnx+k*nnxnny]+cnsr[(nnx-1)+(0+1)*nnx+k*nnxnny]-2*cnsr[(nnx-1)+0*nnx+k*nnxnny])/taunsrT+
      (cnsr[(nnx-1)+0*nnx+(k+1)*nnxnny]+cnsr[(nnx-1)+0*nnx+(k-1)*nnxnny]-2*cnsr[(nnx-1)+0*nnx+k*nnxnny])/taunsrT;
    Icnsr[(nnx-1)+(nny-1)*nnx+k*nnxnny]=(cnsr[(nnx-1)+(nny-1)*nnx+k*nnxnny-1]+cnsr[(nnx-1)+(nny-1)*nnx+k*nnxnny-1]-2*cnsr[(nnx-1)+(nny-1)*nnx+k*nnxnny])/taunsrL+
      (cnsr[(nnx-1)+((nny-1)-1)*nnx+k*nnxnny]+cnsr[(nnx-1)+((nny-1)-1)*nnx+k*nnxnny]-2*cnsr[(nnx-1)+(nny-1)*nnx+k*nnxnny])/taunsrT+
      (cnsr[(nnx-1)+(nny-1)*nnx+(k+1)*nnxnny]+cnsr[(nnx-1)+(nny-1)*nnx+(k-1)*nnxnny]-2*cnsr[(nnx-1)+(nny-1)*nnx+k*nnxnny])/taunsrT;
  }
#endif
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int k=1;k<nnz-1;k++)
  {
    for (int j=1;j<nny-1;j++)
    {
#pragma ivdep
#pragma vector always
      for (int i=1;i<nnx-1;i++)
      {
        Icnsr[i+j*nnx+k*nnxnny]=(cnsr[i+j*nnx+k*nnxnny+1]+cnsr[i+j*nnx+k*nnxnny-1]-2*cnsr[i+j*nnx+k*nnxnny])/taunsrL+
          (cnsr[i+(j+1)*nnx+k*nnxnny]+cnsr[i+(j-1)*nnx+k*nnxnny]-2*cnsr[i+j*nnx+k*nnxnny])/taunsrT+
          (cnsr[i+j*nnx+(k+1)*nnxnny]+cnsr[i+j*nnx+(k-1)*nnxnny]-2*cnsr[i+j*nnx+k*nnxnny])/taunsrT;
      }
    }
  }
}


#ifdef ___NO_CS_BUFFER
void CSubcell::computecsmn(void)
{
#ifdef ___PERIODIC
  csmn[0+0*nnx+0*nnxnny]=(cs[0+0*nnx+0*nnxnny+1]+cs[(nnx-1)+0*nnx+0*nnxnny])/tausL+
    (cs[0+(0+1)*nnx+0*nnxnny]+cs[0+(nny-1)*nnx+0*nnxnny])/tausT+
    (cs[0+0*nnx+(0+1)*nnxnny]+cs[0+0*nnx+(nnz-1)*nnxnny])/tausT;
  csmn[0+(nny-1)*nnx+(nnz-1)*nnxnny]=(cs[0+(nny-1)*nnx+(nnz-1)*nnxnny+1]+cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/tausL+
    (cs[0+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cs[0+(0)*nnx+(nnz-1)*nnxnny])/tausT+
    (cs[0+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cs[0+(nny-1)*nnx+(0)*nnxnny])/tausT;
  csmn[0+(nny-1)*nnx+0*nnxnny]=(cs[0+(nny-1)*nnx+0*nnxnny+1]+cs[(nnx-1)+(nny-1)*nnx+0*nnxnny])/tausL+
    (cs[0+((nny-1)-1)*nnx+0*nnxnny]+cs[0+(0)*nnx+0*nnxnny])/tausT+
    (cs[0+(nny-1)*nnx+(0+1)*nnxnny]+cs[0+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT;
  csmn[0+0*nnx+(nnz-1)*nnxnny]=(cs[0+0*nnx+(nnz-1)*nnxnny+1]+cs[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/tausL+
    (cs[0+(0+1)*nnx+(nnz-1)*nnxnny]+cs[0+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT+
    (cs[0+0*nnx+((nnz-1)-1)*nnxnny]+cs[0+0*nnx+(0)*nnxnny])/tausT;
  csmn[(nnx-1)+(nny-1)*nnx+0*nnxnny]=(cs[(nnx-1)+(nny-1)*nnx+0*nnxnny-1]+cs[(0)+(nny-1)*nnx+0*nnxnny])/tausL+
    (cs[(nnx-1)+((nny-1)-1)*nnx+0*nnxnny]+cs[(nnx-1)+(0)*nnx+0*nnxnny])/tausT+
    (cs[(nnx-1)+(nny-1)*nnx+(0+1)*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT;
  csmn[(nnx-1)+0*nnx+(nnz-1)*nnxnny]=(cs[(nnx-1)+0*nnx+(nnz-1)*nnxnny-1]+cs[(0)+0*nnx+(nnz-1)*nnxnny])/tausL+
    (cs[(nnx-1)+(0+1)*nnx+(nnz-1)*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT+
    (cs[(nnx-1)+0*nnx+((nnz-1)-1)*nnxnny]+cs[(nnx-1)+0*nnx+(0)*nnxnny])/tausT;
  csmn[(nnx-1)+0*nnx+0*nnxnny]=(cs[(nnx-1)+0*nnx+0*nnxnny-1]+cs[(0)+0*nnx+0*nnxnny])/tausL+
    (cs[(nnx-1)+(0+1)*nnx+0*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+0*nnxnny])/tausT+
    (cs[(nnx-1)+0*nnx+(0+1)*nnxnny]+cs[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/tausT;
  csmn[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny]=(cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny-1]+cs[(0)+(nny-1)*nnx+(nnz-1)*nnxnny])/tausL+
    (cs[(nnx-1)+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cs[(nnx-1)+(0)*nnx+(nnz-1)*nnxnny])/tausT+
    (cs[(nnx-1)+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+(0)*nnxnny])/tausT;
  //x fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int j=1;j<nny-1;j++)
  {
    csmn[0+j*nnx+0*nnxnny]=(cs[0+j*nnx+0*nnxnny+1]+cs[(nnx-1)+j*nnx+0*nnxnny])/tausL+
      (cs[0+(j+1)*nnx+0*nnxnny]+cs[0+(j-1)*nnx+0*nnxnny])/tausT+
      (cs[0+j*nnx+(0+1)*nnxnny]+cs[0+j*nnx+(nnz-1)*nnxnny])/tausT;
    csmn[(nnx-1)+j*nnx+0*nnxnny]=(cs[(nnx-1)+j*nnx+0*nnxnny-1]+cs[(0)+j*nnx+0*nnxnny])/tausL+
      (cs[(nnx-1)+(j+1)*nnx+0*nnxnny]+cs[(nnx-1)+(j-1)*nnx+0*nnxnny])/tausT+
      (cs[(nnx-1)+j*nnx+(0+1)*nnxnny]+cs[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/tausT;
    csmn[0+j*nnx+(nnz-1)*nnxnny]=(cs[0+j*nnx+(nnz-1)*nnxnny+1]+cs[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/tausL+
      (cs[0+(j+1)*nnx+(nnz-1)*nnxnny]+cs[0+(j-1)*nnx+(nnz-1)*nnxnny])/tausT+
      (cs[0+j*nnx+((nnz-1)-1)*nnxnny]+cs[0+j*nnx+(0)*nnxnny])/tausT;
    csmn[(nnx-1)+j*nnx+(nnz-1)*nnxnny]=(cs[(nnx-1)+j*nnx+(nnz-1)*nnxnny-1]+cs[(0)+j*nnx+(nnz-1)*nnxnny])/tausL+
      (cs[(nnx-1)+(j+1)*nnx+(nnz-1)*nnxnny]+cs[(nnx-1)+(j-1)*nnx+(nnz-1)*nnxnny])/tausT+
      (cs[(nnx-1)+j*nnx+((nnz-1)-1)*nnxnny]+cs[(nnx-1)+j*nnx+(0)*nnxnny])/tausT;
#pragma ivdep
#pragma vector always
    for (int k=1;k<nnz-1;k++)
    {
      csmn[0+j*nnx+k*nnxnny]=(cs[0+j*nnx+k*nnxnny+1]+cs[(nnx-1)+j*nnx+k*nnxnny])/tausL+
        (cs[0+(j+1)*nnx+k*nnxnny]+cs[0+(j-1)*nnx+k*nnxnny])/tausT+
        (cs[0+j*nnx+(k+1)*nnxnny]+cs[0+j*nnx+(k-1)*nnxnny])/tausT;
      csmn[(nnx-1)+j*nnx+k*nnxnny]=(cs[(nnx-1)+j*nnx+k*nnxnny-1]+cs[(0)+j*nnx+k*nnxnny])/tausL+
        (cs[(nnx-1)+(j+1)*nnx+k*nnxnny]+cs[(nnx-1)+(j-1)*nnx+k*nnxnny])/tausT+
        (cs[(nnx-1)+j*nnx+(k+1)*nnxnny]+cs[(nnx-1)+j*nnx+(k-1)*nnxnny])/tausT;
    }
  }
  //y fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=1;i<(nnx-1);i++)
  {
    csmn[i+0*nnx+0*nnxnny]=(cs[i+0*nnx+0*nnxnny+1]+cs[i+0*nnx+0*nnxnny-1])/tausL+
      (cs[i+(0+1)*nnx+0*nnxnny]+cs[i+(nny-1)*nnx+0*nnxnny])/tausT+
      (cs[i+0*nnx+(0+1)*nnxnny]+cs[i+0*nnx+(nnz-1)*nnxnny])/tausT;
    csmn[i+(nny-1)*nnx+0*nnxnny]=(cs[i+(nny-1)*nnx+0*nnxnny+1]+cs[i+(nny-1)*nnx+0*nnxnny-1])/tausL+
      (cs[i+((nny-1)-1)*nnx+0*nnxnny]+cs[i+(0)*nnx+0*nnxnny])/tausT+
      (cs[i+(nny-1)*nnx+(0+1)*nnxnny]+cs[i+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT;
    csmn[i+0*nnx+(nnz-1)*nnxnny]=(cs[i+0*nnx+(nnz-1)*nnxnny+1]+cs[i+0*nnx+(nnz-1)*nnxnny-1])/tausL+
      (cs[i+(0+1)*nnx+(nnz-1)*nnxnny]+cs[i+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT+
      (cs[i+0*nnx+((nnz-1)-1)*nnxnny]+cs[i+0*nnx+(0)*nnxnny])/tausT;
    csmn[i+(nny-1)*nnx+(nnz-1)*nnxnny]=(cs[i+(nny-1)*nnx+(nnz-1)*nnxnny+1]+cs[i+(nny-1)*nnx+(nnz-1)*nnxnny-1])/tausL+
      (cs[i+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cs[i+(0)*nnx+(nnz-1)*nnxnny])/tausT+
      (cs[i+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cs[i+(nny-1)*nnx+(0)*nnxnny])/tausT;
#pragma ivdep
#pragma vector always
    for (int k=1;k<nnz-1;k++)
    {
      csmn[i+0*nnx+k*nnxnny]=(cs[i+0*nnx+k*nnxnny+1]+cs[i+0*nnx+k*nnxnny-1])/tausL+
        (cs[i+(0+1)*nnx+k*nnxnny]+cs[i+(nny-1)*nnx+k*nnxnny])/tausT+
        (cs[i+0*nnx+(k+1)*nnxnny]+cs[i+0*nnx+(k-1)*nnxnny])/tausT;
      csmn[i+(nny-1)*nnx+k*nnxnny]=(cs[i+(nny-1)*nnx+k*nnxnny+1]+cs[i+(nny-1)*nnx+k*nnxnny-1])/tausL+
        (cs[i+((nny-1)-1)*nnx+k*nnxnny]+cs[i+(0)*nnx+k*nnxnny])/tausT+
        (cs[i+(nny-1)*nnx+(k+1)*nnxnny]+cs[i+(nny-1)*nnx+(k-1)*nnxnny])/tausT;
    }
    //z fixed
#pragma ivdep
#pragma vector always
    for (int j=1;j<nny-1;j++)
    {
      csmn[i+j*nnx+0*nnxnny]=(cs[i+j*nnx+0*nnxnny+1]+cs[i+j*nnx+0*nnxnny-1])/tausL+
        (cs[i+(j+1)*nnx+0*nnxnny]+cs[i+(j-1)*nnx+0*nnxnny])/tausT+
        (cs[i+j*nnx+(0+1)*nnxnny]+cs[i+j*nnx+(nnz-1)*nnxnny])/tausT;
      csmn[i+j*nnx+(nnz-1)*nnxnny]=(cs[i+j*nnx+(nnz-1)*nnxnny+1]+cs[i+j*nnx+(nnz-1)*nnxnny-1])/tausL+
        (cs[i+(j+1)*nnx+(nnz-1)*nnxnny]+cs[i+(j-1)*nnx+(nnz-1)*nnxnny])/tausT+
        (cs[i+j*nnx+((nnz-1)-1)*nnxnny]+cs[i+j*nnx+(0)*nnxnny])/tausT;
    }
  }
#pragma ivdep
#pragma vector always
  for (int k=1;k<nnz-1;k++)
  {
    csmn[0+0*nnx+k*nnxnny]=(cs[0+0*nnx+k*nnxnny+1]+cs[(nnx-1)+0*nnx+k*nnxnny])/tausL+
      (cs[0+(0+1)*nnx+k*nnxnny]+cs[0+(nny-1)*nnx+k*nnxnny])/tausT+
      (cs[0+0*nnx+(k+1)*nnxnny]+cs[0+0*nnx+(k-1)*nnxnny])/tausT;
    csmn[0+(nny-1)*nnx+k*nnxnny]=(cs[0+(nny-1)*nnx+k*nnxnny+1]+cs[(nnx-1)+(nny-1)*nnx+k*nnxnny])/tausL+
      (cs[0+((nny-1)-1)*nnx+k*nnxnny]+cs[0+(0)*nnx+k*nnxnny])/tausT+
      (cs[0+(nny-1)*nnx+(k+1)*nnxnny]+cs[0+(nny-1)*nnx+(k-1)*nnxnny])/tausT;
    csmn[(nnx-1)+0*nnx+k*nnxnny]=(cs[(nnx-1)+0*nnx+k*nnxnny-1]+cs[(0)+0*nnx+k*nnxnny])/tausL+
      (cs[(nnx-1)+(0+1)*nnx+k*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+k*nnxnny])/tausT+
      (cs[(nnx-1)+0*nnx+(k+1)*nnxnny]+cs[(nnx-1)+0*nnx+(k-1)*nnxnny])/tausT;
    csmn[(nnx-1)+(nny-1)*nnx+k*nnxnny]=(cs[(nnx-1)+(nny-1)*nnx+k*nnxnny-1]+cs[0+(nny-1)*nnx+k*nnxnny])/tausL+
      (cs[(nnx-1)+((nny-1)-1)*nnx+k*nnxnny]+cs[(nnx-1)+(0)*nnx+k*nnxnny])/tausT+
      (cs[(nnx-1)+(nny-1)*nnx+(k+1)*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+(k-1)*nnxnny])/tausT;
  }
#else
  csmn[0+0*nnx+0*nnxnny]=(cs[0+0*nnx+0*nnxnny+1]+cs[0+0*nnx+0*nnxnny+1])/tausL+
    (cs[0+(0+1)*nnx+0*nnxnny]+cs[0+(0+1)*nnx+0*nnxnny])/tausT+
    (cs[0+0*nnx+(0+1)*nnxnny]+cs[0+0*nnx+(0+1)*nnxnny])/tausT;
  csmn[0+(nny-1)*nnx+(nnz-1)*nnxnny]=(cs[0+(nny-1)*nnx+(nnz-1)*nnxnny+1]+cs[0+(nny-1)*nnx+(nnz-1)*nnxnny+1])/tausL+
    (cs[0+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cs[0+((nny-1)-1)*nnx+(nnz-1)*nnxnny])/tausT+
    (cs[0+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cs[0+(nny-1)*nnx+((nnz-1)-1)*nnxnny])/tausT;
  csmn[0+(nny-1)*nnx+0*nnxnny]=(cs[0+(nny-1)*nnx+0*nnxnny+1]+cs[0+(nny-1)*nnx+0*nnxnny+1])/tausL+
    (cs[0+((nny-1)-1)*nnx+0*nnxnny]+cs[0+((nny-1)-1)*nnx+0*nnxnny])/tausT+
    (cs[0+(nny-1)*nnx+(0+1)*nnxnny]+cs[0+(nny-1)*nnx+(0+1)*nnxnny])/tausT;
  csmn[0+0*nnx+(nnz-1)*nnxnny]=(cs[0+0*nnx+(nnz-1)*nnxnny+1]+cs[0+0*nnx+(nnz-1)*nnxnny+1])/tausL+
    (cs[0+(0+1)*nnx+(nnz-1)*nnxnny]+cs[0+(0+1)*nnx+(nnz-1)*nnxnny])/tausT+
    (cs[0+0*nnx+((nnz-1)-1)*nnxnny]+cs[0+0*nnx+((nnz-1)-1)*nnxnny])/tausT;
  csmn[(nnx-1)+(nny-1)*nnx+0*nnxnny]=(cs[(nnx-1)+(nny-1)*nnx+0*nnxnny-1]+cs[(nnx-1)+(nny-1)*nnx+0*nnxnny-1])/tausL+
    (cs[(nnx-1)+((nny-1)-1)*nnx+0*nnxnny]+cs[(nnx-1)+((nny-1)-1)*nnx+0*nnxnny])/tausT+
    (cs[(nnx-1)+(nny-1)*nnx+(0+1)*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+(0+1)*nnxnny])/tausT;
  csmn[(nnx-1)+0*nnx+(nnz-1)*nnxnny]=(cs[(nnx-1)+0*nnx+(nnz-1)*nnxnny-1]+cs[(nnx-1)+0*nnx+(nnz-1)*nnxnny-1])/tausL+
    (cs[(nnx-1)+(0+1)*nnx+(nnz-1)*nnxnny]+cs[(nnx-1)+(0+1)*nnx+(nnz-1)*nnxnny])/tausT+
    (cs[(nnx-1)+0*nnx+((nnz-1)-1)*nnxnny]+cs[(nnx-1)+0*nnx+((nnz-1)-1)*nnxnny])/tausT;
  csmn[(nnx-1)+0*nnx+0*nnxnny]=(cs[(nnx-1)+0*nnx+0*nnxnny-1]+cs[(nnx-1)+0*nnx+0*nnxnny-1])/tausL+
    (cs[(nnx-1)+(0+1)*nnx+0*nnxnny]+cs[(nnx-1)+(0+1)*nnx+0*nnxnny])/tausT+
    (cs[(nnx-1)+0*nnx+(0+1)*nnxnny]+cs[(nnx-1)+0*nnx+(0+1)*nnxnny])/tausT;
  csmn[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny]=(cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny-1]+cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny-1])/tausL+
    (cs[(nnx-1)+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cs[(nnx-1)+((nny-1)-1)*nnx+(nnz-1)*nnxnny])/tausT+
    (cs[(nnx-1)+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+((nnz-1)-1)*nnxnny])/tausT;
  //x fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int j=1;j<nny-1;j++)
  {
    csmn[0+j*nnx+0*nnxnny]=(cs[0+j*nnx+0*nnxnny+1]+cs[0+j*nnx+0*nnxnny+1])/tausL+
      (cs[0+(j+1)*nnx+0*nnxnny]+cs[0+(j-1)*nnx+0*nnxnny])/tausT+
      (cs[0+j*nnx+(0+1)*nnxnny]+cs[0+j*nnx+(0+1)*nnxnny])/tausT;
    csmn[(nnx-1)+j*nnx+0*nnxnny]=(cs[(nnx-1)+j*nnx+0*nnxnny-1]+cs[(nnx-1)+j*nnx+0*nnxnny-1])/tausL+
      (cs[(nnx-1)+(j+1)*nnx+0*nnxnny]+cs[(nnx-1)+(j-1)*nnx+0*nnxnny])/tausT+
      (cs[(nnx-1)+j*nnx+(0+1)*nnxnny]+cs[(nnx-1)+j*nnx+(0+1)*nnxnny])/tausT;
    csmn[0+j*nnx+(nnz-1)*nnxnny]=(cs[0+j*nnx+(nnz-1)*nnxnny+1]+cs[0+j*nnx+(nnz-1)*nnxnny+1])/tausL+
      (cs[0+(j+1)*nnx+(nnz-1)*nnxnny]+cs[0+(j-1)*nnx+(nnz-1)*nnxnny])/tausT+
      (cs[0+j*nnx+((nnz-1)-1)*nnxnny]+cs[0+j*nnx+((nnz-1)-1)*nnxnny])/tausT;
    csmn[(nnx-1)+j*nnx+(nnz-1)*nnxnny]=(cs[(nnx-1)+j*nnx+(nnz-1)*nnxnny-1]+cs[(nnx-1)+j*nnx+(nnz-1)*nnxnny-1])/tausL+
      (cs[(nnx-1)+(j+1)*nnx+(nnz-1)*nnxnny]+cs[(nnx-1)+(j-1)*nnx+(nnz-1)*nnxnny])/tausT+
      (cs[(nnx-1)+j*nnx+((nnz-1)-1)*nnxnny]+cs[(nnx-1)+j*nnx+((nnz-1)-1)*nnxnny])/tausT;
#pragma ivdep
#pragma vector always
    for (int k=1;k<nnz-1;k++)
    {
      csmn[0+j*nnx+k*nnxnny]=(cs[0+j*nnx+k*nnxnny+1]+cs[0+j*nnx+k*nnxnny+1])/tausL+
        (cs[0+(j+1)*nnx+k*nnxnny]+cs[0+(j-1)*nnx+k*nnxnny])/tausT+
        (cs[0+j*nnx+(k+1)*nnxnny]+cs[0+j*nnx+(k-1)*nnxnny])/tausT;
      csmn[(nnx-1)+j*nnx+k*nnxnny]=(cs[(nnx-1)+j*nnx+k*nnxnny-1]+cs[(nnx-1)+j*nnx+k*nnxnny-1])/tausL+
        (cs[(nnx-1)+(j+1)*nnx+k*nnxnny]+cs[(nnx-1)+(j-1)*nnx+k*nnxnny])/tausT+
        (cs[(nnx-1)+j*nnx+(k+1)*nnxnny]+cs[(nnx-1)+j*nnx+(k-1)*nnxnny])/tausT;
    }
  }
  //y fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=1;i<(nnx-1);i++)
  {
    csmn[i+0*nnx+0*nnxnny]=(cs[i+0*nnx+0*nnxnny+1]+cs[i+0*nnx+0*nnxnny-1])/tausL+
      (cs[i+(0+1)*nnx+0*nnxnny]+cs[i+(0+1)*nnx+0*nnxnny])/tausT+
      (cs[i+0*nnx+(0+1)*nnxnny]+cs[i+0*nnx+(0+1)*nnxnny])/tausT;
    csmn[i+(nny-1)*nnx+0*nnxnny]=(cs[i+(nny-1)*nnx+0*nnxnny+1]+cs[i+(nny-1)*nnx+0*nnxnny-1])/tausL+
      (cs[i+((nny-1)-1)*nnx+0*nnxnny]+cs[i+((nny-1)-1)*nnx+0*nnxnny])/tausT+
      (cs[i+(nny-1)*nnx+(0+1)*nnxnny]+cs[i+(nny-1)*nnx+(0+1)*nnxnny])/tausT;
    csmn[i+0*nnx+(nnz-1)*nnxnny]=(cs[i+0*nnx+(nnz-1)*nnxnny+1]+cs[i+0*nnx+(nnz-1)*nnxnny-1])/tausL+
      (cs[i+(0+1)*nnx+(nnz-1)*nnxnny]+cs[i+(0+1)*nnx+(nnz-1)*nnxnny])/tausT+
      (cs[i+0*nnx+((nnz-1)-1)*nnxnny]+cs[i+0*nnx+((nnz-1)-1)*nnxnny])/tausT;
    csmn[i+(nny-1)*nnx+(nnz-1)*nnxnny]=(cs[i+(nny-1)*nnx+(nnz-1)*nnxnny+1]+cs[i+(nny-1)*nnx+(nnz-1)*nnxnny-1])/tausL+
      (cs[i+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cs[i+((nny-1)-1)*nnx+(nnz-1)*nnxnny])/tausT+
      (cs[i+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cs[i+(nny-1)*nnx+((nnz-1)-1)*nnxnny])/tausT;
#pragma ivdep
#pragma vector always
    for (int k=1;k<nnz-1;k++)
    {
      csmn[i+0*nnx+k*nnxnny]=(cs[i+0*nnx+k*nnxnny+1]+cs[i+0*nnx+k*nnxnny-1])/tausL+
        (cs[i+(0+1)*nnx+k*nnxnny]+cs[i+(0+1)*nnx+k*nnxnny])/tausT+
        (cs[i+0*nnx+(k+1)*nnxnny]+cs[i+0*nnx+(k-1)*nnxnny])/tausT;
      csmn[i+(nny-1)*nnx+k*nnxnny]=(cs[i+(nny-1)*nnx+k*nnxnny+1]+cs[i+(nny-1)*nnx+k*nnxnny-1])/tausL+
        (cs[i+((nny-1)-1)*nnx+k*nnxnny]+cs[i+((nny-1)-1)*nnx+k*nnxnny])/tausT+
        (cs[i+(nny-1)*nnx+(k+1)*nnxnny]+cs[i+(nny-1)*nnx+(k-1)*nnxnny])/tausT;
    }
    //z fixed
#pragma ivdep
#pragma vector always
    for (int j=1;j<nny-1;j++)
    {
      csmn[i+j*nnx+0*nnxnny]=(cs[i+j*nnx+0*nnxnny+1]+cs[i+j*nnx+0*nnxnny-1])/tausL+
        (cs[i+(j+1)*nnx+0*nnxnny]+cs[i+(j-1)*nnx+0*nnxnny])/tausT+
        (cs[i+j*nnx+(0+1)*nnxnny]+cs[i+j*nnx+(0+1)*nnxnny])/tausT;
      csmn[i+j*nnx+(nnz-1)*nnxnny]=(cs[i+j*nnx+(nnz-1)*nnxnny+1]+cs[i+j*nnx+(nnz-1)*nnxnny-1])/tausL+
        (cs[i+(j+1)*nnx+(nnz-1)*nnxnny]+cs[i+(j-1)*nnx+(nnz-1)*nnxnny])/tausT+
        (cs[i+j*nnx+((nnz-1)-1)*nnxnny]+cs[i+j*nnx+((nnz-1)-1)*nnxnny])/tausT;
    }
  }
#pragma ivdep
#pragma vector always
  for (int k=1;k<nnz-1;k++)
  {
    csmn[0+0*nnx+k*nnxnny]=(cs[0+0*nnx+k*nnxnny+1]+cs[0+0*nnx+k*nnxnny+1])/tausL+
      (cs[0+(0+1)*nnx+k*nnxnny]+cs[0+(0+1)*nnx+k*nnxnny])/tausT+
      (cs[0+0*nnx+(k+1)*nnxnny]+cs[0+0*nnx+(k-1)*nnxnny])/tausT;
    csmn[0+(nny-1)*nnx+k*nnxnny]=(cs[0+(nny-1)*nnx+k*nnxnny+1]+cs[0+(nny-1)*nnx+k*nnxnny+1])/tausL+
      (cs[0+((nny-1)-1)*nnx+k*nnxnny]+cs[0+((nny-1)-1)*nnx+k*nnxnny])/tausT+
      (cs[0+(nny-1)*nnx+(k+1)*nnxnny]+cs[0+(nny-1)*nnx+(k-1)*nnxnny])/tausT;
    csmn[(nnx-1)+0*nnx+k*nnxnny]=(cs[(nnx-1)+0*nnx+k*nnxnny-1]+cs[(nnx-1)+0*nnx+k*nnxnny-1])/tausL+
      (cs[(nnx-1)+(0+1)*nnx+k*nnxnny]+cs[(nnx-1)+(0+1)*nnx+k*nnxnny])/tausT+
      (cs[(nnx-1)+0*nnx+(k+1)*nnxnny]+cs[(nnx-1)+0*nnx+(k-1)*nnxnny])/tausT;
    csmn[(nnx-1)+(nny-1)*nnx+k*nnxnny]=(cs[(nnx-1)+(nny-1)*nnx+k*nnxnny-1]+cs[(nnx-1)+(nny-1)*nnx+k*nnxnny-1])/tausL+
      (cs[(nnx-1)+((nny-1)-1)*nnx+k*nnxnny]+cs[(nnx-1)+((nny-1)-1)*nnx+k*nnxnny])/tausT+
      (cs[(nnx-1)+(nny-1)*nnx+(k+1)*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+(k-1)*nnxnny])/tausT;
  }
#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int k=1;k<nnz-1;k++)
  {
    for (int j=1;j<nny-1;j++)
    {
#pragma ivdep
#pragma vector always
      for (int i=1;i<nnx-1;i++)
      {
        csmn[i+j*nnx+k*nnxnny]=(cs[i+j*nnx+k*nnxnny+1]+cs[i+j*nnx+k*nnxnny-1])/tausL+
          (cs[i+(j+1)*nnx+k*nnxnny]+cs[i+(j-1)*nnx+k*nnxnny])/tausT+
          (cs[i+j*nnx+(k+1)*nnxnny]+cs[i+j*nnx+(k-1)*nnxnny])/tausT;
      }
    }
  }
}
#else
void CSubcell::computeIcs(void)
{

#ifdef ___PERIODIC
  Ics[0+0*nnx+0*nnxnny]=(cs[0+0*nnx+0*nnxnny+1]+cs[(nnx-1)+0*nnx+0*nnxnny]-2*cs[0+0*nnx+0*nnxnny])/tausL+
    (cs[0+(0+1)*nnx+0*nnxnny]+cs[0+(nny-1)*nnx+0*nnxnny]-2*cs[0+0*nnx+0*nnxnny])/tausT+
    (cs[0+0*nnx+(0+1)*nnxnny]+cs[0+0*nnx+(nnz-1)*nnxnny]-2*cs[0+0*nnx+0*nnxnny])/tausT;
  Ics[0+(nny-1)*nnx+(nnz-1)*nnxnny]=(cs[0+(nny-1)*nnx+(nnz-1)*nnxnny+1]+cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny]-2*cs[0+(nny-1)*nnx+(nnz-1)*nnxnny])/tausL+
    (cs[0+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cs[0+(0)*nnx+(nnz-1)*nnxnny]-2*cs[0+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT+
    (cs[0+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cs[0+(nny-1)*nnx+(0)*nnxnny]-2*cs[0+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT;
  Ics[0+(nny-1)*nnx+0*nnxnny]=(cs[0+(nny-1)*nnx+0*nnxnny+1]+cs[(nnx-1)+(nny-1)*nnx+0*nnxnny]-2*cs[0+(nny-1)*nnx+0*nnxnny])/tausL+
    (cs[0+((nny-1)-1)*nnx+0*nnxnny]+cs[0+(0)*nnx+0*nnxnny]-2*cs[0+(nny-1)*nnx+0*nnxnny])/tausT+
    (cs[0+(nny-1)*nnx+(0+1)*nnxnny]+cs[0+(nny-1)*nnx+(nnz-1)*nnxnny]-2*cs[0+(nny-1)*nnx+0*nnxnny])/tausT;
  Ics[0+0*nnx+(nnz-1)*nnxnny]=(cs[0+0*nnx+(nnz-1)*nnxnny+1]+cs[(nnx-1)+0*nnx+(nnz-1)*nnxnny]-2*cs[0+0*nnx+(nnz-1)*nnxnny])/tausL+
    (cs[0+(0+1)*nnx+(nnz-1)*nnxnny]+cs[0+(nny-1)*nnx+(nnz-1)*nnxnny]-2*cs[0+0*nnx+(nnz-1)*nnxnny])/tausT+
    (cs[0+0*nnx+((nnz-1)-1)*nnxnny]+cs[0+0*nnx+(0)*nnxnny]-2*cs[0+0*nnx+(nnz-1)*nnxnny])/tausT;
  Ics[(nnx-1)+(nny-1)*nnx+0*nnxnny]=(cs[(nnx-1)+(nny-1)*nnx+0*nnxnny-1]+cs[(0)+(nny-1)*nnx+0*nnxnny]-2*cs[(nnx-1)+(nny-1)*nnx+0*nnxnny])/tausL+
    (cs[(nnx-1)+((nny-1)-1)*nnx+0*nnxnny]+cs[(nnx-1)+(0)*nnx+0*nnxnny]-2*cs[(nnx-1)+(nny-1)*nnx+0*nnxnny])/tausT+
    (cs[(nnx-1)+(nny-1)*nnx+(0+1)*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny]-2*cs[(nnx-1)+(nny-1)*nnx+0*nnxnny])/tausT;
  Ics[(nnx-1)+0*nnx+(nnz-1)*nnxnny]=(cs[(nnx-1)+0*nnx+(nnz-1)*nnxnny-1]+cs[(0)+0*nnx+(nnz-1)*nnxnny]-2*cs[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/tausL+
    (cs[(nnx-1)+(0+1)*nnx+(nnz-1)*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny]-2*cs[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/tausT+
    (cs[(nnx-1)+0*nnx+((nnz-1)-1)*nnxnny]+cs[(nnx-1)+0*nnx+(0)*nnxnny]-2*cs[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/tausT;
  Ics[(nnx-1)+0*nnx+0*nnxnny]=(cs[(nnx-1)+0*nnx+0*nnxnny-1]+cs[(0)+0*nnx+0*nnxnny]-2*cs[(nnx-1)+0*nnx+0*nnxnny])/tausL+
    (cs[(nnx-1)+(0+1)*nnx+0*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+0*nnxnny]-2*cs[(nnx-1)+0*nnx+0*nnxnny])/tausT+
    (cs[(nnx-1)+0*nnx+(0+1)*nnxnny]+cs[(nnx-1)+0*nnx+(nnz-1)*nnxnny]-2*cs[(nnx-1)+0*nnx+0*nnxnny])/tausT;
  Ics[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny]=(cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny-1]+cs[(0)+(nny-1)*nnx+(nnz-1)*nnxnny]-2*cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/tausL+
    (cs[(nnx-1)+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cs[(nnx-1)+(0)*nnx+(nnz-1)*nnxnny]-2*cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT+
    (cs[(nnx-1)+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+(0)*nnxnny]-2*cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT;
  //x fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int j=1;j<nny-1;j++)
  {
    Ics[0+j*nnx+0*nnxnny]=(cs[0+j*nnx+0*nnxnny+1]+cs[(nnx-1)+j*nnx+0*nnxnny]-2*cs[0+j*nnx+0*nnxnny])/tausL+
      (cs[0+(j+1)*nnx+0*nnxnny]+cs[0+(j-1)*nnx+0*nnxnny]-2*cs[0+j*nnx+0*nnxnny])/tausT+
      (cs[0+j*nnx+(0+1)*nnxnny]+cs[0+j*nnx+(nnz-1)*nnxnny]-2*cs[0+j*nnx+0*nnxnny])/tausT;
    Ics[(nnx-1)+j*nnx+0*nnxnny]=(cs[(nnx-1)+j*nnx+0*nnxnny-1]+cs[(0)+j*nnx+0*nnxnny]-2*cs[(nnx-1)+j*nnx+0*nnxnny])/tausL+
      (cs[(nnx-1)+(j+1)*nnx+0*nnxnny]+cs[(nnx-1)+(j-1)*nnx+0*nnxnny]-2*cs[(nnx-1)+j*nnx+0*nnxnny])/tausT+
      (cs[(nnx-1)+j*nnx+(0+1)*nnxnny]+cs[(nnx-1)+j*nnx+(nnz-1)*nnxnny]-2*cs[(nnx-1)+j*nnx+0*nnxnny])/tausT;
    Ics[0+j*nnx+(nnz-1)*nnxnny]=(cs[0+j*nnx+(nnz-1)*nnxnny+1]+cs[(nnx-1)+j*nnx+(nnz-1)*nnxnny]-2*cs[0+j*nnx+(nnz-1)*nnxnny])/tausL+
      (cs[0+(j+1)*nnx+(nnz-1)*nnxnny]+cs[0+(j-1)*nnx+(nnz-1)*nnxnny]-2*cs[0+j*nnx+(nnz-1)*nnxnny])/tausT+
      (cs[0+j*nnx+((nnz-1)-1)*nnxnny]+cs[0+j*nnx+(0)*nnxnny]-2*cs[0+j*nnx+(nnz-1)*nnxnny])/tausT;
    Ics[(nnx-1)+j*nnx+(nnz-1)*nnxnny]=(cs[(nnx-1)+j*nnx+(nnz-1)*nnxnny-1]+cs[(0)+j*nnx+(nnz-1)*nnxnny]-2*cs[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/tausL+
      (cs[(nnx-1)+(j+1)*nnx+(nnz-1)*nnxnny]+cs[(nnx-1)+(j-1)*nnx+(nnz-1)*nnxnny]-2*cs[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/tausT+
      (cs[(nnx-1)+j*nnx+((nnz-1)-1)*nnxnny]+cs[(nnx-1)+j*nnx+(0)*nnxnny]-2*cs[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/tausT;
#pragma ivdep
#pragma vector always
    for (int k=1;k<nnz-1;k++)
    {
      Ics[0+j*nnx+k*nnxnny]=(cs[0+j*nnx+k*nnxnny+1]+cs[(nnx-1)+j*nnx+k*nnxnny]-2*cs[0+j*nnx+k*nnxnny])/tausL+
        (cs[0+(j+1)*nnx+k*nnxnny]+cs[0+(j-1)*nnx+k*nnxnny]-2*cs[0+j*nnx+k*nnxnny])/tausT+
        (cs[0+j*nnx+(k+1)*nnxnny]+cs[0+j*nnx+(k-1)*nnxnny]-2*cs[0+j*nnx+k*nnxnny])/tausT;
      Ics[(nnx-1)+j*nnx+k*nnxnny]=(cs[(nnx-1)+j*nnx+k*nnxnny-1]+cs[(0)+j*nnx+k*nnxnny]-2*cs[(nnx-1)+j*nnx+k*nnxnny])/tausL+
        (cs[(nnx-1)+(j+1)*nnx+k*nnxnny]+cs[(nnx-1)+(j-1)*nnx+k*nnxnny]-2*cs[(nnx-1)+j*nnx+k*nnxnny])/tausT+
        (cs[(nnx-1)+j*nnx+(k+1)*nnxnny]+cs[(nnx-1)+j*nnx+(k-1)*nnxnny]-2*cs[(nnx-1)+j*nnx+k*nnxnny])/tausT;
    }
  }
  //y fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=1;i<(nnx-1);i++)
  {
    Ics[i+0*nnx+0*nnxnny]=(cs[i+0*nnx+0*nnxnny+1]+cs[i+0*nnx+0*nnxnny-1]-2*cs[i+0*nnx+0*nnxnny])/tausL+
      (cs[i+(0+1)*nnx+0*nnxnny]+cs[i+(nny-1)*nnx+0*nnxnny]-2*cs[i+0*nnx+0*nnxnny])/tausT+
      (cs[i+0*nnx+(0+1)*nnxnny]+cs[i+0*nnx+(nnz-1)*nnxnny]-2*cs[i+0*nnx+0*nnxnny])/tausT;
    Ics[i+(nny-1)*nnx+0*nnxnny]=(cs[i+(nny-1)*nnx+0*nnxnny+1]+cs[i+(nny-1)*nnx+0*nnxnny-1]-2*cs[i+(nny-1)*nnx+0*nnxnny])/tausL+
      (cs[i+((nny-1)-1)*nnx+0*nnxnny]+cs[i+(0)*nnx+0*nnxnny]-2*cs[i+(nny-1)*nnx+0*nnxnny])/tausT+
      (cs[i+(nny-1)*nnx+(0+1)*nnxnny]+cs[i+(nny-1)*nnx+(nnz-1)*nnxnny]-2*cs[i+(nny-1)*nnx+0*nnxnny])/tausT;
    Ics[i+0*nnx+(nnz-1)*nnxnny]=(cs[i+0*nnx+(nnz-1)*nnxnny+1]+cs[i+0*nnx+(nnz-1)*nnxnny-1]-2*cs[i+0*nnx+(nnz-1)*nnxnny])/tausL+
      (cs[i+(0+1)*nnx+(nnz-1)*nnxnny]+cs[i+(nny-1)*nnx+(nnz-1)*nnxnny]-2*cs[i+0*nnx+(nnz-1)*nnxnny])/tausT+
      (cs[i+0*nnx+((nnz-1)-1)*nnxnny]+cs[i+0*nnx+(0)*nnxnny]-2*cs[i+0*nnx+(nnz-1)*nnxnny])/tausT;
    Ics[i+(nny-1)*nnx+(nnz-1)*nnxnny]=(cs[i+(nny-1)*nnx+(nnz-1)*nnxnny+1]+cs[i+(nny-1)*nnx+(nnz-1)*nnxnny-1]-2*cs[i+(nny-1)*nnx+(nnz-1)*nnxnny])/tausL+
      (cs[i+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cs[i+(0)*nnx+(nnz-1)*nnxnny]-2*cs[i+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT+
      (cs[i+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cs[i+(nny-1)*nnx+(0)*nnxnny]-2*cs[i+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT;
#pragma ivdep
#pragma vector always
    for (int k=1;k<nnz-1;k++)
    {
      Ics[i+0*nnx+k*nnxnny]=(cs[i+0*nnx+k*nnxnny+1]+cs[i+0*nnx+k*nnxnny-1]-2*cs[i+0*nnx+k*nnxnny])/tausL+
        (cs[i+(0+1)*nnx+k*nnxnny]+cs[i+(nny-1)*nnx+k*nnxnny]-2*cs[i+0*nnx+k*nnxnny])/tausT+
        (cs[i+0*nnx+(k+1)*nnxnny]+cs[i+0*nnx+(k-1)*nnxnny]-2*cs[i+0*nnx+k*nnxnny])/tausT;
      Ics[i+(nny-1)*nnx+k*nnxnny]=(cs[i+(nny-1)*nnx+k*nnxnny+1]+cs[i+(nny-1)*nnx+k*nnxnny-1]-2*cs[i+(nny-1)*nnx+k*nnxnny])/tausL+
        (cs[i+((nny-1)-1)*nnx+k*nnxnny]+cs[i+(0)*nnx+k*nnxnny]-2*cs[i+(nny-1)*nnx+k*nnxnny])/tausT+
        (cs[i+(nny-1)*nnx+(k+1)*nnxnny]+cs[i+(nny-1)*nnx+(k-1)*nnxnny]-2*cs[i+(nny-1)*nnx+k*nnxnny])/tausT;
    }
    //z fixed
#pragma ivdep
#pragma vector always
    for (int j=1;j<nny-1;j++)
    {
      Ics[i+j*nnx+0*nnxnny]=(cs[i+j*nnx+0*nnxnny+1]+cs[i+j*nnx+0*nnxnny-1]-2*cs[i+j*nnx+0*nnxnny])/tausL+
        (cs[i+(j+1)*nnx+0*nnxnny]+cs[i+(j-1)*nnx+0*nnxnny]-2*cs[i+j*nnx+0*nnxnny])/tausT+
        (cs[i+j*nnx+(0+1)*nnxnny]+cs[i+j*nnx+(nnz-1)*nnxnny]-2*cs[i+j*nnx+0*nnxnny])/tausT;
      Ics[i+j*nnx+(nnz-1)*nnxnny]=(cs[i+j*nnx+(nnz-1)*nnxnny+1]+cs[i+j*nnx+(nnz-1)*nnxnny-1]-2*cs[i+j*nnx+(nnz-1)*nnxnny])/tausL+
        (cs[i+(j+1)*nnx+(nnz-1)*nnxnny]+cs[i+(j-1)*nnx+(nnz-1)*nnxnny]-2*cs[i+j*nnx+(nnz-1)*nnxnny])/tausT+
        (cs[i+j*nnx+((nnz-1)-1)*nnxnny]+cs[i+j*nnx+(0)*nnxnny]-2*cs[i+j*nnx+(nnz-1)*nnxnny])/tausT;
    }
  }
#pragma ivdep
#pragma vector always
  for (int k=1;k<nnz-1;k++)
  {
    Ics[0+0*nnx+k*nnxnny]=(cs[0+0*nnx+k*nnxnny+1]+cs[(nnx-1)+0*nnx+k*nnxnny]-2*cs[0+0*nnx+k*nnxnny])/tausL+
      (cs[0+(0+1)*nnx+k*nnxnny]+cs[0+(nny-1)*nnx+k*nnxnny]-2*cs[0+0*nnx+k*nnxnny])/tausT+
      (cs[0+0*nnx+(k+1)*nnxnny]+cs[0+0*nnx+(k-1)*nnxnny]-2*cs[0+0*nnx+k*nnxnny])/tausT;
    Ics[0+(nny-1)*nnx+k*nnxnny]=(cs[0+(nny-1)*nnx+k*nnxnny+1]+cs[(nnx-1)+(nny-1)*nnx+k*nnxnny]-2*cs[0+(nny-1)*nnx+k*nnxnny])/tausL+
      (cs[0+((nny-1)-1)*nnx+k*nnxnny]+cs[0+(0)*nnx+k*nnxnny]-2*cs[0+(nny-1)*nnx+k*nnxnny])/tausT+
      (cs[0+(nny-1)*nnx+(k+1)*nnxnny]+cs[0+(nny-1)*nnx+(k-1)*nnxnny]-2*cs[0+(nny-1)*nnx+k*nnxnny])/tausT;
    Ics[(nnx-1)+0*nnx+k*nnxnny]=(cs[(nnx-1)+0*nnx+k*nnxnny-1]+cs[(0)+0*nnx+k*nnxnny]-2*cs[(nnx-1)+0*nnx+k*nnxnny])/tausL+
      (cs[(nnx-1)+(0+1)*nnx+k*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+k*nnxnny]-2*cs[(nnx-1)+0*nnx+k*nnxnny])/tausT+
      (cs[(nnx-1)+0*nnx+(k+1)*nnxnny]+cs[(nnx-1)+0*nnx+(k-1)*nnxnny]-2*cs[(nnx-1)+0*nnx+k*nnxnny])/tausT;
    Ics[(nnx-1)+(nny-1)*nnx+k*nnxnny]=(cs[(nnx-1)+(nny-1)*nnx+k*nnxnny-1]+cs[(0)+(nny-1)*nnx+k*nnxnny]-2*cs[(nnx-1)+(nny-1)*nnx+k*nnxnny])/tausL+
      (cs[(nnx-1)+((nny-1)-1)*nnx+k*nnxnny]+cs[(nnx-1)+(0)*nnx+k*nnxnny]-2*cs[(nnx-1)+(nny-1)*nnx+k*nnxnny])/tausT+
      (cs[(nnx-1)+(nny-1)*nnx+(k+1)*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+(k-1)*nnxnny]-2*cs[(nnx-1)+(nny-1)*nnx+k*nnxnny])/tausT;
  }
#else
  Ics[0+0*nnx+0*nnxnny]=(cs[0+0*nnx+0*nnxnny+1]+cs[0+0*nnx+0*nnxnny+1]-2*cs[0+0*nnx+0*nnxnny])/tausL+
    (cs[0+(0+1)*nnx+0*nnxnny]+cs[0+(0+1)*nnx+0*nnxnny]-2*cs[0+0*nnx+0*nnxnny])/tausT+
    (cs[0+0*nnx+(0+1)*nnxnny]+cs[0+0*nnx+(0+1)*nnxnny]-2*cs[0+0*nnx+0*nnxnny])/tausT;
  Ics[0+(nny-1)*nnx+(nnz-1)*nnxnny]=(cs[0+(nny-1)*nnx+(nnz-1)*nnxnny+1]+cs[0+(nny-1)*nnx+(nnz-1)*nnxnny+1]-2*cs[0+(nny-1)*nnx+(nnz-1)*nnxnny])/tausL+
    (cs[0+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cs[0+((nny-1)-1)*nnx+(nnz-1)*nnxnny]-2*cs[0+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT+
    (cs[0+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cs[0+(nny-1)*nnx+((nnz-1)-1)*nnxnny]-2*cs[0+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT;
  Ics[0+(nny-1)*nnx+0*nnxnny]=(cs[0+(nny-1)*nnx+0*nnxnny+1]+cs[0+(nny-1)*nnx+0*nnxnny+1]-2*cs[0+(nny-1)*nnx+0*nnxnny])/tausL+
    (cs[0+((nny-1)-1)*nnx+0*nnxnny]+cs[0+((nny-1)-1)*nnx+0*nnxnny]-2*cs[0+(nny-1)*nnx+0*nnxnny])/tausT+
    (cs[0+(nny-1)*nnx+(0+1)*nnxnny]+cs[0+(nny-1)*nnx+(0+1)*nnxnny]-2*cs[0+(nny-1)*nnx+0*nnxnny])/tausT;
  Ics[0+0*nnx+(nnz-1)*nnxnny]=(cs[0+0*nnx+(nnz-1)*nnxnny+1]+cs[0+0*nnx+(nnz-1)*nnxnny+1]-2*cs[0+0*nnx+(nnz-1)*nnxnny])/tausL+
    (cs[0+(0+1)*nnx+(nnz-1)*nnxnny]+cs[0+(0+1)*nnx+(nnz-1)*nnxnny]-2*cs[0+0*nnx+(nnz-1)*nnxnny])/tausT+
    (cs[0+0*nnx+((nnz-1)-1)*nnxnny]+cs[0+0*nnx+((nnz-1)-1)*nnxnny]-2*cs[0+0*nnx+(nnz-1)*nnxnny])/tausT;
  Ics[(nnx-1)+(nny-1)*nnx+0*nnxnny]=(cs[(nnx-1)+(nny-1)*nnx+0*nnxnny-1]+cs[(nnx-1)+(nny-1)*nnx+0*nnxnny-1]-2*cs[(nnx-1)+(nny-1)*nnx+0*nnxnny])/tausL+
    (cs[(nnx-1)+((nny-1)-1)*nnx+0*nnxnny]+cs[(nnx-1)+((nny-1)-1)*nnx+0*nnxnny]-2*cs[(nnx-1)+(nny-1)*nnx+0*nnxnny])/tausT+
    (cs[(nnx-1)+(nny-1)*nnx+(0+1)*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+(0+1)*nnxnny]-2*cs[(nnx-1)+(nny-1)*nnx+0*nnxnny])/tausT;
  Ics[(nnx-1)+0*nnx+(nnz-1)*nnxnny]=(cs[(nnx-1)+0*nnx+(nnz-1)*nnxnny-1]+cs[(nnx-1)+0*nnx+(nnz-1)*nnxnny-1]-2*cs[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/tausL+
    (cs[(nnx-1)+(0+1)*nnx+(nnz-1)*nnxnny]+cs[(nnx-1)+(0+1)*nnx+(nnz-1)*nnxnny]-2*cs[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/tausT+
    (cs[(nnx-1)+0*nnx+((nnz-1)-1)*nnxnny]+cs[(nnx-1)+0*nnx+((nnz-1)-1)*nnxnny]-2*cs[(nnx-1)+0*nnx+(nnz-1)*nnxnny])/tausT;
  Ics[(nnx-1)+0*nnx+0*nnxnny]=(cs[(nnx-1)+0*nnx+0*nnxnny-1]+cs[(nnx-1)+0*nnx+0*nnxnny-1]-2*cs[(nnx-1)+0*nnx+0*nnxnny])/tausL+
    (cs[(nnx-1)+(0+1)*nnx+0*nnxnny]+cs[(nnx-1)+(0+1)*nnx+0*nnxnny]-2*cs[(nnx-1)+0*nnx+0*nnxnny])/tausT+
    (cs[(nnx-1)+0*nnx+(0+1)*nnxnny]+cs[(nnx-1)+0*nnx+(0+1)*nnxnny]-2*cs[(nnx-1)+0*nnx+0*nnxnny])/tausT;
  Ics[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny]=(cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny-1]+cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny-1]-2*cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/tausL+
    (cs[(nnx-1)+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cs[(nnx-1)+((nny-1)-1)*nnx+(nnz-1)*nnxnny]-2*cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT+
    (cs[(nnx-1)+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+((nnz-1)-1)*nnxnny]-2*cs[(nnx-1)+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT;
  //x fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int j=1;j<nny-1;j++)
  {
    Ics[0+j*nnx+0*nnxnny]=(cs[0+j*nnx+0*nnxnny+1]+cs[0+j*nnx+0*nnxnny+1]-2*cs[0+j*nnx+0*nnxnny])/tausL+
      (cs[0+(j+1)*nnx+0*nnxnny]+cs[0+(j-1)*nnx+0*nnxnny]-2*cs[0+j*nnx+0*nnxnny])/tausT+
      (cs[0+j*nnx+(0+1)*nnxnny]+cs[0+j*nnx+(0+1)*nnxnny]-2*cs[0+j*nnx+0*nnxnny])/tausT;
    Ics[(nnx-1)+j*nnx+0*nnxnny]=(cs[(nnx-1)+j*nnx+0*nnxnny-1]+cs[(nnx-1)+j*nnx+0*nnxnny-1]-2*cs[(nnx-1)+j*nnx+0*nnxnny])/tausL+
      (cs[(nnx-1)+(j+1)*nnx+0*nnxnny]+cs[(nnx-1)+(j-1)*nnx+0*nnxnny]-2*cs[(nnx-1)+j*nnx+0*nnxnny])/tausT+
      (cs[(nnx-1)+j*nnx+(0+1)*nnxnny]+cs[(nnx-1)+j*nnx+(0+1)*nnxnny]-2*cs[(nnx-1)+j*nnx+0*nnxnny])/tausT;
    Ics[0+j*nnx+(nnz-1)*nnxnny]=(cs[0+j*nnx+(nnz-1)*nnxnny+1]+cs[0+j*nnx+(nnz-1)*nnxnny+1]-2*cs[0+j*nnx+(nnz-1)*nnxnny])/tausL+
      (cs[0+(j+1)*nnx+(nnz-1)*nnxnny]+cs[0+(j-1)*nnx+(nnz-1)*nnxnny]-2*cs[0+j*nnx+(nnz-1)*nnxnny])/tausT+
      (cs[0+j*nnx+((nnz-1)-1)*nnxnny]+cs[0+j*nnx+((nnz-1)-1)*nnxnny]-2*cs[0+j*nnx+(nnz-1)*nnxnny])/tausT;
    Ics[(nnx-1)+j*nnx+(nnz-1)*nnxnny]=(cs[(nnx-1)+j*nnx+(nnz-1)*nnxnny-1]+cs[(nnx-1)+j*nnx+(nnz-1)*nnxnny-1]-2*cs[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/tausL+
      (cs[(nnx-1)+(j+1)*nnx+(nnz-1)*nnxnny]+cs[(nnx-1)+(j-1)*nnx+(nnz-1)*nnxnny]-2*cs[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/tausT+
      (cs[(nnx-1)+j*nnx+((nnz-1)-1)*nnxnny]+cs[(nnx-1)+j*nnx+((nnz-1)-1)*nnxnny]-2*cs[(nnx-1)+j*nnx+(nnz-1)*nnxnny])/tausT;
#pragma ivdep
#pragma vector always
    for (int k=1;k<nnz-1;k++)
    {
      Ics[0+j*nnx+k*nnxnny]=(cs[0+j*nnx+k*nnxnny+1]+cs[0+j*nnx+k*nnxnny+1]-2*cs[0+j*nnx+k*nnxnny])/tausL+
        (cs[0+(j+1)*nnx+k*nnxnny]+cs[0+(j-1)*nnx+k*nnxnny]-2*cs[0+j*nnx+k*nnxnny])/tausT+
        (cs[0+j*nnx+(k+1)*nnxnny]+cs[0+j*nnx+(k-1)*nnxnny]-2*cs[0+j*nnx+k*nnxnny])/tausT;
      Ics[(nnx-1)+j*nnx+k*nnxnny]=(cs[(nnx-1)+j*nnx+k*nnxnny-1]+cs[(nnx-1)+j*nnx+k*nnxnny-1]-2*cs[(nnx-1)+j*nnx+k*nnxnny])/tausL+
        (cs[(nnx-1)+(j+1)*nnx+k*nnxnny]+cs[(nnx-1)+(j-1)*nnx+k*nnxnny]-2*cs[(nnx-1)+j*nnx+k*nnxnny])/tausT+
        (cs[(nnx-1)+j*nnx+(k+1)*nnxnny]+cs[(nnx-1)+j*nnx+(k-1)*nnxnny]-2*cs[(nnx-1)+j*nnx+k*nnxnny])/tausT;
    }
  }
  //y fixed
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int i=1;i<(nnx-1);i++)
  {
    Ics[i+0*nnx+0*nnxnny]=(cs[i+0*nnx+0*nnxnny+1]+cs[i+0*nnx+0*nnxnny-1]-2*cs[i+0*nnx+0*nnxnny])/tausL+
      (cs[i+(0+1)*nnx+0*nnxnny]+cs[i+(0+1)*nnx+0*nnxnny]-2*cs[i+0*nnx+0*nnxnny])/tausT+
      (cs[i+0*nnx+(0+1)*nnxnny]+cs[i+0*nnx+(0+1)*nnxnny]-2*cs[i+0*nnx+0*nnxnny])/tausT;
    Ics[i+(nny-1)*nnx+0*nnxnny]=(cs[i+(nny-1)*nnx+0*nnxnny+1]+cs[i+(nny-1)*nnx+0*nnxnny-1]-2*cs[i+(nny-1)*nnx+0*nnxnny])/tausL+
      (cs[i+((nny-1)-1)*nnx+0*nnxnny]+cs[i+((nny-1)-1)*nnx+0*nnxnny]-2*cs[i+(nny-1)*nnx+0*nnxnny])/tausT+
      (cs[i+(nny-1)*nnx+(0+1)*nnxnny]+cs[i+(nny-1)*nnx+(0+1)*nnxnny]-2*cs[i+(nny-1)*nnx+0*nnxnny])/tausT;
    Ics[i+0*nnx+(nnz-1)*nnxnny]=(cs[i+0*nnx+(nnz-1)*nnxnny+1]+cs[i+0*nnx+(nnz-1)*nnxnny-1]-2*cs[i+0*nnx+(nnz-1)*nnxnny])/tausL+
      (cs[i+(0+1)*nnx+(nnz-1)*nnxnny]+cs[i+(0+1)*nnx+(nnz-1)*nnxnny]-2*cs[i+0*nnx+(nnz-1)*nnxnny])/tausT+
      (cs[i+0*nnx+((nnz-1)-1)*nnxnny]+cs[i+0*nnx+((nnz-1)-1)*nnxnny]-2*cs[i+0*nnx+(nnz-1)*nnxnny])/tausT;
    Ics[i+(nny-1)*nnx+(nnz-1)*nnxnny]=(cs[i+(nny-1)*nnx+(nnz-1)*nnxnny+1]+cs[i+(nny-1)*nnx+(nnz-1)*nnxnny-1]-2*cs[i+(nny-1)*nnx+(nnz-1)*nnxnny])/tausL+
      (cs[i+((nny-1)-1)*nnx+(nnz-1)*nnxnny]+cs[i+((nny-1)-1)*nnx+(nnz-1)*nnxnny]-2*cs[i+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT+
      (cs[i+(nny-1)*nnx+((nnz-1)-1)*nnxnny]+cs[i+(nny-1)*nnx+((nnz-1)-1)*nnxnny]-2*cs[i+(nny-1)*nnx+(nnz-1)*nnxnny])/tausT;
#pragma ivdep
#pragma vector always
    for (int k=1;k<nnz-1;k++)
    {
      Ics[i+0*nnx+k*nnxnny]=(cs[i+0*nnx+k*nnxnny+1]+cs[i+0*nnx+k*nnxnny-1]-2*cs[i+0*nnx+k*nnxnny])/tausL+
        (cs[i+(0+1)*nnx+k*nnxnny]+cs[i+(0+1)*nnx+k*nnxnny]-2*cs[i+0*nnx+k*nnxnny])/tausT+
        (cs[i+0*nnx+(k+1)*nnxnny]+cs[i+0*nnx+(k-1)*nnxnny]-2*cs[i+0*nnx+k*nnxnny])/tausT;
      Ics[i+(nny-1)*nnx+k*nnxnny]=(cs[i+(nny-1)*nnx+k*nnxnny+1]+cs[i+(nny-1)*nnx+k*nnxnny-1]-2*cs[i+(nny-1)*nnx+k*nnxnny])/tausL+
        (cs[i+((nny-1)-1)*nnx+k*nnxnny]+cs[i+((nny-1)-1)*nnx+k*nnxnny]-2*cs[i+(nny-1)*nnx+k*nnxnny])/tausT+
        (cs[i+(nny-1)*nnx+(k+1)*nnxnny]+cs[i+(nny-1)*nnx+(k-1)*nnxnny]-2*cs[i+(nny-1)*nnx+k*nnxnny])/tausT;
    }
    //z fixed
#pragma ivdep
#pragma vector always
    for (int j=1;j<nny-1;j++)
    {
      Ics[i+j*nnx+0*nnxnny]=(cs[i+j*nnx+0*nnxnny+1]+cs[i+j*nnx+0*nnxnny-1]-2*cs[i+j*nnx+0*nnxnny])/tausL+
        (cs[i+(j+1)*nnx+0*nnxnny]+cs[i+(j-1)*nnx+0*nnxnny]-2*cs[i+j*nnx+0*nnxnny])/tausT+
        (cs[i+j*nnx+(0+1)*nnxnny]+cs[i+j*nnx+(0+1)*nnxnny]-2*cs[i+j*nnx+0*nnxnny])/tausT;
      Ics[i+j*nnx+(nnz-1)*nnxnny]=(cs[i+j*nnx+(nnz-1)*nnxnny+1]+cs[i+j*nnx+(nnz-1)*nnxnny-1]-2*cs[i+j*nnx+(nnz-1)*nnxnny])/tausL+
        (cs[i+(j+1)*nnx+(nnz-1)*nnxnny]+cs[i+(j-1)*nnx+(nnz-1)*nnxnny]-2*cs[i+j*nnx+(nnz-1)*nnxnny])/tausT+
        (cs[i+j*nnx+((nnz-1)-1)*nnxnny]+cs[i+j*nnx+((nnz-1)-1)*nnxnny]-2*cs[i+j*nnx+(nnz-1)*nnxnny])/tausT;
    }
  }
#pragma ivdep
#pragma vector always
  for (int k=1;k<nnz-1;k++)
  {
    Ics[0+0*nnx+k*nnxnny]=(cs[0+0*nnx+k*nnxnny+1]+cs[0+0*nnx+k*nnxnny+1]-2*cs[0+0*nnx+k*nnxnny])/tausL+
      (cs[0+(0+1)*nnx+k*nnxnny]+cs[0+(0+1)*nnx+k*nnxnny]-2*cs[0+0*nnx+k*nnxnny])/tausT+
      (cs[0+0*nnx+(k+1)*nnxnny]+cs[0+0*nnx+(k-1)*nnxnny]-2*cs[0+0*nnx+k*nnxnny])/tausT;
    Ics[0+(nny-1)*nnx+k*nnxnny]=(cs[0+(nny-1)*nnx+k*nnxnny+1]+cs[0+(nny-1)*nnx+k*nnxnny+1]-2*cs[0+(nny-1)*nnx+k*nnxnny])/tausL+
      (cs[0+((nny-1)-1)*nnx+k*nnxnny]+cs[0+((nny-1)-1)*nnx+k*nnxnny]-2*cs[0+(nny-1)*nnx+k*nnxnny])/tausT+
      (cs[0+(nny-1)*nnx+(k+1)*nnxnny]+cs[0+(nny-1)*nnx+(k-1)*nnxnny]-2*cs[0+(nny-1)*nnx+k*nnxnny])/tausT;
    Ics[(nnx-1)+0*nnx+k*nnxnny]=(cs[(nnx-1)+0*nnx+k*nnxnny-1]+cs[(nnx-1)+0*nnx+k*nnxnny-1]-2*cs[(nnx-1)+0*nnx+k*nnxnny])/tausL+
      (cs[(nnx-1)+(0+1)*nnx+k*nnxnny]+cs[(nnx-1)+(0+1)*nnx+k*nnxnny]-2*cs[(nnx-1)+0*nnx+k*nnxnny])/tausT+
      (cs[(nnx-1)+0*nnx+(k+1)*nnxnny]+cs[(nnx-1)+0*nnx+(k-1)*nnxnny]-2*cs[(nnx-1)+0*nnx+k*nnxnny])/tausT;
    Ics[(nnx-1)+(nny-1)*nnx+k*nnxnny]=(cs[(nnx-1)+(nny-1)*nnx+k*nnxnny-1]+cs[(nnx-1)+(nny-1)*nnx+k*nnxnny-1]-2*cs[(nnx-1)+(nny-1)*nnx+k*nnxnny])/tausL+
      (cs[(nnx-1)+((nny-1)-1)*nnx+k*nnxnny]+cs[(nnx-1)+((nny-1)-1)*nnx+k*nnxnny]-2*cs[(nnx-1)+(nny-1)*nnx+k*nnxnny])/tausT+
      (cs[(nnx-1)+(nny-1)*nnx+(k+1)*nnxnny]+cs[(nnx-1)+(nny-1)*nnx+(k-1)*nnxnny]-2*cs[(nnx-1)+(nny-1)*nnx+k*nnxnny])/tausT;
  }

#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int k=1;k<nnz-1;k++)
  {
    for (int j=1;j<nny-1;j++)
    {
#pragma ivdep
#pragma vector always
      for (int i=1;i<nnx-1;i++)
      {
        Ics[i+j*nnx+k*nnxnny]=(cs[i+j*nnx+k*nnxnny+1]+cs[i+j*nnx+k*nnxnny-1]-2*cs[i+j*nnx+k*nnxnny])/tausL+
          (cs[i+(j+1)*nnx+k*nnxnny]+cs[i+(j-1)*nnx+k*nnxnny]-2*cs[i+j*nnx+k*nnxnny])/tausT+
          (cs[i+j*nnx+(k+1)*nnxnny]+cs[i+j*nnx+(k-1)*nnxnny]-2*cs[i+j*nnx+k*nnxnny])/tausT;
      }
    }
  }
}

#endif


#ifdef ___PTM
void CSubcell::allocate_memory_all_PTM_vars(int n){

  //CaM Differentials
  dydt_CaMDyad = new double[15*n];
  dydt_CaMSL = new double[15*n];
  dydt_CaMCyt = new double[15*n];
  dydt_CaMKII = new double[6*n];
  dydt_BAR = new double[36*n];

  //CaM State Variables
  //Ca Fluxes from CaM_dyad
  JCaCyt = new double[n];
  JCaSL = new double[n];
  JCaDyad = new double[n];

  //CaM_Dyad
  CaM_dyad = new double[n];
  Ca2CaM_dyad = new double[n];
  Ca4CaM_dyad = new double[n];
  CaMB_dyad = new double[n];
  Ca2CaMB_dyad = new double[n];
  Ca4CaMB_dyad = new double[n];
  Pb2_dyad = new double[n];
  Pb_dyad = new double[n];
  Pt_dyad = new double[n];
  Pt2_dyad = new double[n];
  Pa_dyad = new double[n];
  Ca4CaN_dyad = new double[n];
  CaMCa4CaN_dyad = new double[n];
  Ca2CaMCa4CaN_dyad = new double[n];
  Ca4CaMCa4CaN_dyad = new double[n];

  //CaM_SL
  CaM_sl = new double[n];
  Ca2CaM_sl = new double[n];
  Ca4CaM_sl = new double[n];
  CaMB_sl = new double[n];
  Ca2CaMB_sl = new double[n];
  Ca4CaMB_sl = new double[n];
  Pb2_sl = new double[n];
  Pb_sl = new double[n];
  Pt_sl = new double[n];
  Pt2_sl = new double[n];
  Pa_sl = new double[n];
  Ca4CaN_sl = new double[n];
  CaMCa4CaN_sl = new double[n];
  Ca2CaMCa4CaN_sl = new double[n];
  Ca4CaMCa4CaN_sl = new double[n];

  //CaM_cytosol
  CaM_cyt = new double[n];
  Ca2CaM_cyt = new double[n];
  Ca4CaM_cyt = new double[n];
  CaMB_cyt = new double[n];
  Ca2CaMB_cyt = new double[n];
  Ca4CaMB_cyt = new double[n];
  Pb2_cyt = new double[n];
  Pb_cyt = new double[n];
  Pt_cyt = new double[n];
  Pt2_cyt = new double[n];
  Pa_cyt = new double[n];
  Ca4CaN_cyt = new double[n];
  CaMCa4CaN_cyt = new double[n];
  Ca2CaMCa4CaN_cyt = new double[n];
  Ca4CaMCa4CaN_cyt = new double[n];

  // //CaMKII// //
  
  CaMKIIactDyad = new double[n]; 
  CaMKIIactSL = new double[n];
  PP1_PLB_avail = new double[n];

  LCC_PKAp = new double[n];
  LCC_CKdyadp = new double[n]; //TODO this variable name is used twice, fix please
  LCC_CKslp = new double[n];
  RyR2809p = new double[n];
  RyR2815p = new double[n];
  PLBT17p = new double[n];

  RyR_CKp = new double[n];
  PLB_CKp = new double[n];

  // //BAR// //
  L = new double[n];
  B1AR = new double[n];
  Gs = new double[n];
  B1AR_ACT = new double[n];
  B1AR_S464 = new double[n];

  B1AR_S301 = new double[n];
  GsaGTPtot = new double[n];
  GsaGDP = new double[n];
  GsBy = new double[n];

  GsaGTP = new double[n];
  Fsk = new double[n];
  AC = new double[n];
  PDE = new double[n];
  IBMX = new double[n];

  cAMPtot = new double[n];
  cAMP = new double[n];
  PKAC_I = new double[n];
  PKAC_II = new double[n];

  PLBp = new double[n];
  Inhib1ptot = new double[n];
  Inhib1p = new double[n];
  PP1 = new double[n];

  LCCa_PKAp_whole = new double[n];
  LCCb_PKAp_whole = new double[n];
  RyR_PKAp_whole = new double[n];
  TnI_PKAp_whole = new double[n];
  IKs_PKAn = new double[n]; 
  Yotiao_KCQN1 = new double[n]; 
  IKs_PKAp_whole = new double[n];
  ICFTR_PKAp_whole = new double[n];
  PLM_PKAp_whole = new double[n];
  Myo_PKAp_whole = new double[n];
  IKr_PKAn = new double[n]; 
  Yotiao_hERG = new double[n]; 
  IKr_PKAp_whole = new double[n];
  IClCa_PKAp_whole = new double[n];

  LCCa_PKAp = new double[n];
  LCCb_PKAp = new double[n];
  PLB_PKAn = new double[n];
  RyR_PKAp = new double[n];
  TnI_PKAp = new double[n];
  IKs_PKAp = new double[n];
  ICFTR_PKAp = new double[n];
  PLM_PKAp = new double[n];
  Myo_PKAp = new double[n];
  IKr_PKAp = new double[n];
  IClCa_PKAp = new double[n];


  RYR_multiplier = new double[n];
}
void CSubcell::init_const_parameters_PTM(){
  // freq = 1.0;                 // [Hz] CHANGE DEPENDING ON FREQUENCY
  // cycleLength = 1e3/freq;     // [ms]
  CaMtotDyad = 418;             //[uM]
  BtotDyad = 1.54/8.293e-4;        //[uM]
  CaMKIItotDyad = 120;        //[uM]
  CaNtotDyad = 3e-3/8.293e-4;    //[uM]
  PP1totDyad = 96.5;             //[uM]
  CaMtotSL = 5.65;                 //[uM]
  BtotSL = 24.2;                     //[uM]
  CaMKIItotSL = 120*8.293e-4;   //[uM]
  CaNtotSL = 3e-3;                 //[uM]
  PP1totSL = 0.57;                 //[uM]
  CaMtotCyt = 5.65;               //[uM]
  BtotCyt = 24.2;                   //[uM]
  CaMKIItotCyt = 120*8.293e-4; //[uM]
  CaNtotCyt = 3e-3;               //[uM]
  PP1totCyt = 0.57;               //[uM]

  // ADJUST CAMKII ACTIVITY LEVELS (expression = 'WT', 'OE', or 'KO')
  expression = WT;
  CKIIOE = 0; // Should be zero during "WT" and "KO" runs
  if(expression == OE){
    CKIIOE = 1; // Flag for CKII OE in ECC file (0=WT, 1=OE) - for Ito and INa
    CaMKIItotDyad = 120*6;          // [uM] 
    CaMKIItotSL = 120*8.293e-4*6;   // [uM]
    CaMKIItotCyt = 120*8.293e-4*6;  // [uM]
  }
  else if(expression == KO){
    CaMKIItotDyad = 0;              // [uM] 
    CaMKIItotSL = 0;                // [uM]
    CaMKIItotCyt = 0;               // [uM]
  }

  // Parameters for CaMKII module
  LCCtotDyad = 31.4*.9;       // [uM] - Total Dyadic [LCC] - (umol/l dyad)
  LCCtotSL = 0.0846;          // [uM] - Total Subsarcolemmal [LCC] (umol/l sl)
  RyRtot = 382.6;             // [uM] - Total RyR (in Dyad)
  PP1_dyad = 95.7;            // [uM] - Total dyadic [PP1]
  PP1_SL = 0.57;              // [uM] - Total Subsarcolemmal [PP1]
  PP2A_dyad = 95.76;          // [uM] - Total dyadic PP2A
  OA = 0;                     // [uM] - PP1/PP2A inhibitor Okadaic Acid
  PLBtot = 38;                // [uM] - Total [PLB] in cytosolic units

  // Parameters for BAR module
  Ligtot = 0;// 0.1 or 0.02    // [uM] - SET LIGAND CONCENTRATION HERE
  LCCtotBA = 0.025;           // [uM] - [umol/L cytosol]
  RyRtotBA = 0.135;           // [uM] - [umol/L cytosol]
  PLBtotBA = PLBtot;          // [uM] - [umol/L cytosol]
  TnItotBA = 70;              // [uM] - [umol/L cytosol]
  IKstotBA = 0.025;           // [uM] - [umol/L cytosol]
  ICFTRtotBA = 0.025;         // [uM] - [umol/L cytosol]
  PP1_PLBtot = 0.89;          // [uM] - [umol/L cytosol]
  PLMtotBA = 48;              // [uM] - [umol/L cytosol]
  MyototBA = 70;              // [uM] - [umol/L cytosol]
  IKrtotBA = 0.025;           // [uM] - [umol/L cytosol]
  IClCatotBA = 0.025;         // [uM] - [umol/L cytosol]

  // For Recovery from inactivation of LCC
  recoveryTime = 10;  // initialize to smallest value
}
void CSubcell::init_state_variables_PTMs(int id){
    // double *dydt_CaMDyad, *dydt_CaMSL, *dydt_CaMCyt;
    // double *JCaCyt, *JCaSL, *JCaSLDyad;

    //CaM_Dyad
    CaM_dyad[id]= 356.109460603695;
    Ca2CaM_dyad[id]= 5.97654735341164;
    Ca4CaM_dyad[id]= 0.00223486677258696;
    CaMB_dyad[id]= 0;
    Ca2CaMB_dyad[id]= 0;
    Ca4CaMB_dyad[id]= 0;
    Pb2_dyad[id]= 0.527046816424795;
    Pb_dyad[id]= 0.0328341389969648;
    Pt_dyad[id]= 7.55606853909437e-06;
    Pt2_dyad[id]= 2.99351358946522e-09;
    Pa_dyad[id]= 2.49374553021126e-09;
    Ca4CaN_dyad[id]= 0.000136542384253690;
    CaMCa4CaN_dyad[id]= 0.00414331349458221;
    Ca2CaMCa4CaN_dyad[id]= 0.0125769071945199;
    Ca4CaMCa4CaN_dyad[id]= 3.60060936227458;

    //CaM_SL
    
    CaM_sl[id]= 0.0369987905424850;
    Ca2CaM_sl[id]= 4.96044557788166e-05;
    Ca4CaM_sl[id]= 1.25071429861889e-08;
    CaMB_sl[id]= 1.99877548871130;
    Ca2CaMB_sl[id]= 13.7576843921400;
    Ca4CaMB_sl[id]= 0.000294277468518282;
    Pb2_sl[id]= 1.00823339154849e-05;
    Pb_sl[id]= 1.49023312731869e-05;
    Pt_sl[id]= 3.28512350734024e-07;
    Pt2_sl[id]= 3.56376700614912e-12;
    Pa_sl[id]= 2.13842749748637e-08;
    Ca4CaN_sl[id]= 0.000363360135943152;
    CaMCa4CaN_sl[id]= 1.15201811885199e-06;
    Ca2CaMCa4CaN_sl[id]= 5.16530421197686e-06;
    Ca4CaMCa4CaN_sl[id]= 0.00204745211775456;

    //CaM_cytosol
    CaM_cyt[id]= 0.0367160202068570;
    Ca2CaM_cyt[id]= 3.48648648558715e-05;
    Ca4CaM_cyt[id]= 5.72196307711396e-10;
    CaMB_cyt[id]= 3.50223317045483;
    Ca2CaMB_cyt[id]= 1.80060664830752;
    Ca4CaMB_cyt[id]= 2.60549585613380e-05;
    Pb2_cyt[id]= 6.97710649793748e-06;
    Pb_cyt[id]= 1.00475277347346e-07;
    Pt_cyt[id]= 3.78363623175676e-12;
    Pt2_cyt[id]= 4.03604680784057e-17;
    Pa_cyt[id]= 2.52209309010385e-13;
    Ca4CaN_cyt[id]= 0.000201018368785850;
    CaMCa4CaN_cyt[id]= 6.32018257082229e-07;
    Ca2CaMCa4CaN_cyt[id]= 3.95696769885873e-07;
    Ca4CaMCa4CaN_cyt[id]= 7.89429746296487e-06;

    //* * * * * CaM Variables * * * * * //

    //* * * * * CaMKII Variables * * * * *  //

    //CaMKII State Variables
    
    LCC_PKAp[id] = 16.4544;
    LCC_CKdyadp[id] = 15.6555;
    RyR2809p[id] = 297.3572;
    RyR2815p[id] = 74.9548;
    PLBT17p[id] = 0.4628;
    LCC_CKslp[id] = 1.90123528218208e-05;


    CaMKIIactDyad[id] = CaMKIItotDyad*(Pb_dyad[id]+Pt_dyad[id]+Pt2_dyad[id]+Pa_dyad[id]); // Multiply total by fraction of activated CaMKII states (Pb, Pt, Pt2, Pa)

    

    CaMKIIactSL[id] = CaMKIItotSL*(Pb_sl[id]+Pt_sl[id]+Pt2_sl[id]+Pa_sl[id]);
    PP1_PLB_avail[id] = 0.8819/PP1_PLBtot + .0091;  // Active PP1 near PLB / total PP1 conc + basal value

    LCC_CKdyadp[id] = LCC_CKdyadp[id]/LCCtotDyad;   //136 fractional CaMKII-dependent LCC dyad phosphorylation
    RyR_CKp[id] = RyR2815p[id]/RyRtot;           //138 fractional CaMKII-dependent RyR phosphorylation
    PLB_CKp[id] = PLBT17p[id]/PLBtot; 
    //* * * * * CaMKII Variables * * * * *  

    //* * * * * B-AR Variables * * * * *  //
    L[id] = 0;
    B1AR[id] = 0.0247551234189601;
    Gs[id] = 3.82712905902435;
    B1AR_ACT[id] = 0.0276260643946097;
    B1AR_S464[id] = 1.47940139415962e-21;
    B1AR_S301[id] = 0.000373935605630215;
    GsaGTPtot[id] = 0.0574188195042173;
    GsaGDP[id] = 0.000654316122780689;
    GsBy[id] = 0.0580191356270164;
    GsaGTP[id] = 0.0520107506131518;
    Fsk[id] = 0;
    AC[id] = 0.0415919311089345;
    PDE[id] = 0.0360000000000000;
    IBMX[id] = 0;
    cAMPtot[id] = 0.486632223132263;
    cAMP[id] = 0.0716671431739411;
    PKAC_I[id] = 0.00827172716628298;
    PKAC_II[id] = 0.00317804750275709;
    PLBp[id] = 0.278595379571156;
    Inhib1ptot[id] = 0.00807506817947479;
    Inhib1p[id] = 9.14572037129325e-06;
    PP1[id] = 0.881934077540897;
    LCCa_PKAp_whole[id] = 0.000695895330926661;
    LCCb_PKAp_whole[id] = 0.000821087078294546;
    RyR_PKAp_whole[id] = 0.00329852374692194;
    TnI_PKAp_whole[id] = 0.219146073024094;
    IKs_PKAn[id] = 0.00153192920195564;
    Yotiao_KCQN1[id] = 0.00153192920195564;
    IKs_PKAp_whole[id] = 0.00274417781311398;
    ICFTR_PKAp_whole[id] = 0.00405926811380601;
    PLM_PKAp_whole[id] = 0.302091135992497;
    Myo_PKAp_whole[id] = 0.219146073024094;
    IKr_PKAn[id] = 0.00153192920195564;
    Yotiao_hERG[id] = 0.00153192920195564;
    IKr_PKAp_whole[id] = 0.00274417781311045;
    IClCa_PKAp_whole[id] = 0.00405926811380601;

    LCCa_PKAp[id] = LCCa_PKAp_whole[id]/LCCtotBA;
    LCCb_PKAp[id] = LCCb_PKAp_whole[id]/LCCtotBA;
    PLB_PKAn[id] = (PLBtotBA - PLBp[id])/PLBtotBA; // non-phosphorylated PLB targets
    RyR_PKAp[id] = RyR_PKAp_whole[id]/RyRtotBA;
    TnI_PKAp[id] = TnI_PKAp_whole[id]/TnItotBA;
    IKs_PKAp[id] = IKs_PKAp_whole[id]/IKstotBA;
    ICFTR_PKAp[id] = ICFTR_PKAp_whole[id]/ICFTRtotBA;
    PLM_PKAp[id] = PLM_PKAp_whole[id]/PLMtotBA;
    Myo_PKAp[id] = Myo_PKAp_whole[id]/MyototBA;
    IKr_PKAp[id] = IKr_PKAp_whole[id]/IKrtotBA;
    IClCa_PKAp[id] = IClCa_PKAp_whole[id]/IClCatotBA;

    //* * * * * B-AR Variables * * * * *  //
}

void CSubcell::calc_dydt_CaM_Dyad_ODEs(int id){
  double K = 135.0; // mM
  double Mg = 1.0; //mM


  // Negroni et al model - CaM_dyad module

  //// Descriptions for state variables
  // Ca2CaM_dyad[id] = y(1);  // 2 Ca bound to C terminal sites
  // Ca4CaM_dyad[id] = y(2);  // 4 Ca bound
  // CaMB_dyad[id] = y(3);   
  // Ca2CaMB_dyad[id] = y(4);
  // Ca4CaMB_dyad[id] = y(5);
  // Pb2_dyad[id] = y(6);     // probability of a Ca2CaM_dyad[id] bound CaMKII subunit
  // Pb_dyad[id] = y(7);      // probability of a Ca4CaM_dyad[id] bound CaMKII subunit
  // Pt_dyad[id] = y(8);      // probability of a Ca4CaM_dyad[id] bound autophosphorylated CaMKII subunit
  // Pt2_dyad[id] = y(9);     // probability of a Ca2CaMxx bound autophosphorylated CaMKII subunit
  // Pa_dyad[id] = y(10);     // probability of an autonomous autophosphorylated CaMKII subunit
  // CaMCaN = y(11);
  // Ca2CaMCaN = y(12);
  // Ca4CaMCaN = y(13);  // active calcineurin
  // description of intermediate variables
  // CaM_dyad- Ca free CaM_dyad
  ////

  // Parameters
  // Ca/CaM_dyad parameters
  double Kd02, Kd24;
  if(Mg <= 1){
      Kd02 = 0.0025*(1.0+K/0.94-Mg/0.012)*(1.0+K/8.1+Mg/0.022);  // [uM^2]
      Kd24 = 0.128*(1.0+K/0.64+Mg/0.0014)*(1.0+K/13.0-Mg/0.153); // [uM^2]
  }
  else{
      Kd02 = 0.0025*(1.0+K/0.94-1.0/0.012+(Mg-1.0)/0.060)*(1.0+K/8.1+1.0/0.022+(Mg-1)/0.068);   // [uM^2]
      Kd24 = 0.128*(1.0+K/0.64+1.0/0.0014+(Mg-1.0)/0.005)*(1.0+K/13.0-1.0/0.153+(Mg-1)/0.150);  // [uM^2]
  }

  double k20 = 10.0;               // [s^-1]      
  double k02 = k20/Kd02;         // [uM^-2 s^-1]
  double k42 = 500.0;              // [s^-1]      
  double k24 = k42/Kd24;         // [uM^-2 s^-1]

  // CaM_dyad buffering (B) parameters
  double k0Boff = 0.0014;        // [s^-1] 
  double k0Bon = k0Boff/0.2;   // [uM^-1 s^-1] kon = koff/Kd
  double k2Boff = k0Boff/100.0;    // [s^-1] 
  double k2Bon = k0Bon;          // [uM^-1 s^-1]
  double k4Boff = k2Boff;        // [s^-1]
  double k4Bon = k0Bon;          // [uM^-1 s^-1]
  // using thermodynamic constraints
  double k20B = k20/100.0; // [s^-1] thermo constraint on loop 1
  double k02B = k02;     // [uM^-2 s^-1] 
  double k42B = k42;     // [s^-1] thermo constraint on loop 2
  double k24B = k24;     // [uM^-2 s^-1]

  // CaMKII parameters
  // Wi Wa Wt Wp
  double kbi = 2.2;      // [s^-1] (Ca4CaM_dyad[id] dissocation from Wb)
  double kib = kbi/33.5e-3; // [uM^-1 s^-1]
  double kib2 = kib;
  double kb2i = kib2*5.0;
  double kb24 = k24;
  double kb42 = k42*33.5e-3/5;
  double kpp1 = 1.72;    // [s^-1] (PP1-dep dephosphorylation rates)
  double Kmpp1 = 11.5;   // [uM]
  double kta = kbi/1000; // [s^-1] (Ca4CaM_dyad[id] dissociation from Wt)
  double kat = kib;      // [uM^-1 s^-1] (Ca4CaM_dyad[id] reassociation with Wa)
  double kt42 = k42*33.5e-6/5;
  double kt24 = k24;
  double kat2 = kib;
  double kt2a = kib*5;

  // CaN parameters
  double kcanCaoff = 1.0;              // [s^-1] 
  double kcanCaon = kcanCaoff/0.5;   // [uM^-1 s^-1] 
  double kcanCaM4on = 46.0;            // [uM^-1 s^-1]
  double kcanCaM4off = 1.3e-3;       // [s^-1]
  double kcanCaM2on = kcanCaM4on;
  double kcanCaM2off = 2508.0*kcanCaM4off;
  double kcanCaM0on = kcanCaM4on;
  double kcanCaM0off = 165.0*kcanCaM2off;
  double k02can = k02;
  double k20can = k20/165.0;
  double k24can = k24;
  double k42can = k20/2508.0;

  // CaM_dyad Reaction fluxes
  double rcn02 = k02*pow(cp[id],2)*CaM_dyad[id] - k20*Ca2CaM_dyad[id];
  double rcn24 = k24*pow(cp[id],2)*Ca2CaM_dyad[id] - k42*Ca4CaM_dyad[id];
  
  // CaM_dyad buffer fluxes
  // double B = BtotDyad - CaMB_dyad[id] - Ca2CaMB_dyad[id] - Ca4CaMB_dyad[id];
  //Note, BtotDyad is set to 0
  double B = 0 - CaMB_dyad[id] - Ca2CaMB_dyad[id] - Ca4CaMB_dyad[id];
  double rcn02B = k02B*pow(cp[id],2)*CaMB_dyad[id] - k20B*Ca2CaMB_dyad[id];
  double rcn24B = k24B*pow(cp[id],2)*Ca2CaMB_dyad[id] - k42B*Ca4CaMB_dyad[id];
  double rcn0B = k0Bon*CaM_dyad[id]*B - k0Boff*CaMB_dyad[id];
  double rcn2B = k2Bon*Ca2CaM_dyad[id]*B - k2Boff*Ca2CaMB_dyad[id];
  double rcn4B = k4Bon*Ca4CaM_dyad[id]*B - k4Boff*Ca4CaMB_dyad[id];
  // CaN reaction fluxes 
  double Ca2CaN = CaNtotDyad - Ca4CaN_dyad[id] - CaMCa4CaN_dyad[id] - Ca2CaMCa4CaN_dyad[id] - Ca4CaMCa4CaN_dyad[id];
  double rcnCa4CaN = kcanCaon*pow(cp[id],2)*Ca2CaN - kcanCaoff*Ca4CaN_dyad[id];
  double rcn02CaN = k02can*pow(cp[id],2)*CaMCa4CaN_dyad[id] - k20can*Ca2CaMCa4CaN_dyad[id]; 
  double rcn24CaN = k24can*pow(cp[id],2)*Ca2CaMCa4CaN_dyad[id] - k42can*Ca4CaMCa4CaN_dyad[id];
  double rcn0CaN = kcanCaM0on*CaM_dyad[id]*Ca4CaN_dyad[id] - kcanCaM0off*CaMCa4CaN_dyad[id];
  double rcn2CaN = kcanCaM2on*Ca2CaM_dyad[id]*Ca4CaN_dyad[id] - kcanCaM2off*Ca2CaMCa4CaN_dyad[id];
  double rcn4CaN = kcanCaM4on*Ca4CaM_dyad[id]*Ca4CaN_dyad[id] - kcanCaM4off*Ca4CaMCa4CaN_dyad[id];
  // CaMKII reaction fluxes
  double Pi = 1.0 - Pb2_dyad[id] - Pb_dyad[id] - Pt_dyad[id] - Pt2_dyad[id] - Pa_dyad[id];
  double rcnCKib2 = kib2*Ca2CaM_dyad[id]*Pi - kb2i*Pb2_dyad[id];
  double rcnCKb2b = kb24*pow(cp[id],2)*Pb2_dyad[id] - kb42*Pb_dyad[id];
  double rcnCKib = kib*Ca4CaM_dyad[id]*Pi - kbi*Pb_dyad[id];
  double T = Pb_dyad[id] + Pt_dyad[id] + Pt2_dyad[id] + Pa_dyad[id];
  double kbt = 0.055*T + .0074*pow(T,2) + 0.015*pow(T,3);
  double rcnCKbt = kbt*Pb_dyad[id] - kpp1*PP1totDyad*Pt_dyad[id]/(Kmpp1+CaMKIItotDyad*Pt_dyad[id]);
  double rcnCKtt2 = kt42*Pt_dyad[id] - kt24*pow(cp[id],2)*Pt2_dyad[id];
  double rcnCKta = kta*Pt_dyad[id] - kat*Ca4CaM_dyad[id]*Pa_dyad[id];
  double rcnCKt2a = kt2a*Pt2_dyad[id] - kat2*Ca2CaM_dyad[id]*Pa_dyad[id];
  double rcnCKt2b2 = kpp1*PP1totDyad*Pt2_dyad[id]/(Kmpp1+CaMKIItotDyad*Pt2_dyad[id]);
  double rcnCKai = kpp1*PP1totDyad*Pa_dyad[id]/(Kmpp1+CaMKIItotDyad*Pa_dyad[id]);

  // CaM_dyad equations
  double dCaM = 1e-3*(-rcn02 - rcn0B - rcn0CaN);
  // if(id == 0){
  //   printf("B: %f, BtotDyad: %f, CaMB_dyad[id]: %f, Ca2CaMB_dyad[id]: %f, Ca4CaMB_dyad[id]: %f\n", B, BtotDyad, CaMB_dyad[id], Ca2CaMB_dyad[id], Ca4CaMB_dyad[id]);
  //   printf("rcn02: %f, rcn0B: %f, rcn0CaN: %f\n", rcn02, rcn0B, rcn0CaN);
  //   printf("cp[id]: %f\n\n", cp[id]);
  // }
  double dCa2CaM = 1e-3*(rcn02 - rcn24 - rcn2B - rcn2CaN + CaMKIItotDyad*(-rcnCKib2 + rcnCKt2a) );
  double dCa4CaM = 1e-3*(rcn24 - rcn4B - rcn4CaN + CaMKIItotDyad*(-rcnCKib+rcnCKta) );
  double dCaMB = 1e-3*(rcn0B-rcn02B);
  double dCa2CaMB = 1e-3*(rcn02B + rcn2B - rcn24B);
  double dCa4CaMB = 1e-3*(rcn24B + rcn4B);

  // CaMKII equations
  double dPb2 = 1e-3*(rcnCKib2 - rcnCKb2b + rcnCKt2b2); // Pb2_dyad[id]
  double dPb = 1e-3*(rcnCKib + rcnCKb2b - rcnCKbt);    // Pb_dyad[id]
  double dPt = 1e-3*(rcnCKbt-rcnCKta-rcnCKtt2);        // Pt_dyad[id]
  double dPt2 = 1e-3*(rcnCKtt2-rcnCKt2a-rcnCKt2b2);     // Pt2_dyad[id]
  double dPa = 1e-3*(rcnCKta+rcnCKt2a-rcnCKai);       // Pa_dyad[id]

  // CaN equations
  double dCa4CaN = 1e-3*(rcnCa4CaN - rcn0CaN - rcn2CaN - rcn4CaN);                       // Ca4CaN_dyad[id]
  double dCaMCa4CaN = 1e-3*(rcn0CaN - rcn02CaN);           // CaMCa4CaN_dyad[id]
  double dCa2CaMCa4CaN = 1e-3*(rcn2CaN+rcn02CaN-rcn24CaN);    // Ca2CaMCa4CaN_dyad[id]
  double dCa4CaMCa4CaN = 1e-3*(rcn4CaN+rcn24CaN);             // Ca4CaMCa4CaN_dyad[id]

  // double dydt[]={dCaM,dCa2CaM,dCa4CaM,dCaMB,dCa2CaMB,dCa4CaMB,dPb2,dPb,dPt,dPt2,dPa,dCa4CaN,dCaMCa4CaN,dCa2CaMCa4CaN,dCa4CaMCa4CaN};
  
  //Assign dydt
  // dydt_CaMDyad = dydt; //dydt for CaMDyad
  // for(int i=0;i<15;i++){
  //   dydt_CaMDyad[(id*15)+i]=dydt[i];
  // }

  // write to global variables for adjusting Ca dyad buffering in EC coupling model
  JCaDyad[id] = 1e-3*(2.0*CaMKIItotDyad*(rcnCKtt2-rcnCKb2b) - 2.0*(rcn02+rcn24+rcn02B+rcn24B+rcnCa4CaN+rcn02CaN+rcn24CaN)); // [uM/msec]

  dydt_CaMDyad[(id*15)+  0] = dCaM;
  dydt_CaMDyad[(id*15)+  1] = dCa2CaM;
  dydt_CaMDyad[(id*15)+  2] = dCa4CaM;
  dydt_CaMDyad[(id*15)+  3] = dCaMB;
  dydt_CaMDyad[(id*15)+  4] = dCa2CaMB;
  dydt_CaMDyad[(id*15)+  5] = dCa4CaMB;
  dydt_CaMDyad[(id*15)+  6] = dPb2;
  dydt_CaMDyad[(id*15)+  7] = dPb;
  dydt_CaMDyad[(id*15)+  8] = dPt;
  dydt_CaMDyad[(id*15)+  9] = dPt2;
  dydt_CaMDyad[(id*15)+ 10] = dPa;
  dydt_CaMDyad[(id*15)+ 11] = dCa4CaN;
  dydt_CaMDyad[(id*15)+ 12] = dCaMCa4CaN;
  dydt_CaMDyad[(id*15)+ 13] = dCa2CaMCa4CaN;
  dydt_CaMDyad[(id*15)+ 14] = dCa4CaMCa4CaN; 
}
void CSubcell::calc_dydt_CaM_SL_ODEs(int id){
  double K = 135.0; // mM
  double Mg = 1.0; //mM


  // Negroni et al model - CaM_sl module

  //// Descriptions for state variables
  // Ca2CaM_sl[id] = y(1);  // 2 Ca bound to C terminal sites
  // Ca4CaM_sl[id] = y(2);  // 4 Ca bound
  // CaMB_sl[id] = y(3);   
  // Ca2CaMB_sl[id] = y(4);
  // Ca4CaMB_sl[id] = y(5);
  // Pb2_sl[id] = y(6);     // probability of a Ca2CaM_sl[id] bound CaMKII subunit
  // Pb_sl[id] = y(7);      // probability of a Ca4CaM_sl[id] bound CaMKII subunit
  // Pt_sl[id] = y(8);      // probability of a Ca4CaM_sl[id] bound autophosphorylated CaMKII subunit
  // Pt2_sl[id] = y(9);     // probability of a Ca2CaMxx bound autophosphorylated CaMKII subunit
  // Pa_sl[id] = y(10);     // probability of an autonomous autophosphorylated CaMKII subunit
  // CaMCaN = y(11);
  // Ca2CaMCaN = y(12);
  // Ca4CaMCaN = y(13);  // active calcineurin
  // description of intermediate variables
  // CaM_sl- Ca free CaM_sl
  ////

  // Parameters
  // Ca/CaM_sl parameters
  double Kd02, Kd24;
  if(Mg <= 1.0){
      Kd02 = 0.0025*(1.0+K/0.94-Mg/0.012)*(1.0+K/8.1+Mg/0.022);  // [uM^2]
      Kd24 = 0.128*(1.0+K/0.64+Mg/0.0014)*(1.0+K/13.0-Mg/0.153); // [uM^2]
  }
  else{
      Kd02 = 0.0025*(1.0+K/0.94-1.0/0.012+(Mg-1.0)/0.060)*(1.0+K/8.1+1.0/0.022+(Mg-1.0)/0.068);   // [uM^2]
      Kd24 = 0.128*(1.0+K/0.64+1.0/0.0014+(Mg-1.0)/0.005)*(1.0+K/13.0-1.0/0.153+(Mg-1.0)/0.150);  // [uM^2]
  }

  double k20 = 10.0;               // [s^-1]      
  double k02 = k20/Kd02;         // [uM^-2 s^-1]
  double k42 = 500.0;              // [s^-1]      
  double k24 = k42/Kd24;         // [uM^-2 s^-1]

  // CaM_sl buffering (B) parameters
  double k0Boff = 0.0014;        // [s^-1] 
  double k0Bon = k0Boff/0.2;   // [uM^-1 s^-1] kon = koff/Kd
  double k2Boff = k0Boff/100.0;    // [s^-1] 
  double k2Bon = k0Bon;          // [uM^-1 s^-1]
  double k4Boff = k2Boff;        // [s^-1]
  double k4Bon = k0Bon;          // [uM^-1 s^-1]
  // using thermodynamic constraints
  double k20B = k20/100.0; // [s^-1] thermo constraint on loop 1
  double k02B = k02;     // [uM^-2 s^-1] 
  double k42B = k42;     // [s^-1] thermo constraint on loop 2
  double k24B = k24;     // [uM^-2 s^-1]

  // CaMKII parameters
  // Wi Wa Wt Wp
  double kbi = 2.2;      // [s^-1] (Ca4CaM_sl[id] dissocation from Wb)
  double kib = kbi/33.5e-3; // [uM^-1 s^-1]
  double kib2 = kib;
  double kb2i = kib2*5.0;
  double kb24 = k24;
  double kb42 = k42*33.5e-3/5.0;
  double kpp1 = 1.72;    // [s^-1] (PP1-dep dephosphorylation rates)
  double Kmpp1 = 11.5;   // [uM]
  double kta = kbi/1000.0; // [s^-1] (Ca4CaM_sl[id] dissociation from Wt)
  double kat = kib;      // [uM^-1 s^-1] (Ca4CaM_sl[id] reassociation with Wa)
  double kt42 = k42*33.5e-6/5.0;
  double kt24 = k24;
  double kat2 = kib;
  double kt2a = kib*5.0;

  // CaN parameters
  double kcanCaoff = 1.0;              // [s^-1] 
  double kcanCaon = kcanCaoff/0.5;   // [uM^-1 s^-1] 
  double kcanCaM4on = 46.0;            // [uM^-1 s^-1]
  double kcanCaM4off = 1.3e-3;       // [s^-1]
  double kcanCaM2on = kcanCaM4on;
  double kcanCaM2off = 2508.0*kcanCaM4off;
  double kcanCaM0on = kcanCaM4on;
  double kcanCaM0off = 165.0*kcanCaM2off;
  double k02can = k02;
  double k20can = k20/165.0;
  double k24can = k24;
  double k42can = k20/2508;

  // CaM_sl Reaction fluxes
  double rcn02 = k02*pow(cs[id],2)*CaM_sl[id] - k20*Ca2CaM_sl[id];
  double rcn24 = k24*pow(cs[id],2)*Ca2CaM_sl[id] - k42*Ca4CaM_sl[id];
  // CaM_sl buffer fluxes
  double B = BtotSL - CaMB_sl[id] - Ca2CaMB_sl[id] - Ca4CaMB_sl[id];
  double rcn02B = k02B*pow(cs[id],2)*CaMB_sl[id] - k20B*Ca2CaMB_sl[id];
  double rcn24B = k24B*pow(cs[id],2)*Ca2CaMB_sl[id] - k42B*Ca4CaMB_sl[id];
  double rcn0B = k0Bon*CaM_sl[id]*B - k0Boff*CaMB_sl[id];
  double rcn2B = k2Bon*Ca2CaM_sl[id]*B - k2Boff*Ca2CaMB_sl[id];
  double rcn4B = k4Bon*Ca4CaM_sl[id]*B - k4Boff*Ca4CaMB_sl[id];
  // CaN reaction fluxes 
  double Ca2CaN = CaNtotSL - Ca4CaN_sl[id] - CaMCa4CaN_sl[id] - Ca2CaMCa4CaN_sl[id] - Ca4CaMCa4CaN_sl[id];
  double rcnCa4CaN = kcanCaon*pow(cs[id],2)*Ca2CaN - kcanCaoff*Ca4CaN_sl[id];
  double rcn02CaN = k02can*pow(cs[id],2)*CaMCa4CaN_sl[id] - k20can*Ca2CaMCa4CaN_sl[id]; 
  double rcn24CaN = k24can*pow(cs[id],2)*Ca2CaMCa4CaN_sl[id] - k42can*Ca4CaMCa4CaN_sl[id];
  double rcn0CaN = kcanCaM0on*CaM_sl[id]*Ca4CaN_sl[id] - kcanCaM0off*CaMCa4CaN_sl[id];
  double rcn2CaN = kcanCaM2on*Ca2CaM_sl[id]*Ca4CaN_sl[id] - kcanCaM2off*Ca2CaMCa4CaN_sl[id];
  double rcn4CaN = kcanCaM4on*Ca4CaM_sl[id]*Ca4CaN_sl[id] - kcanCaM4off*Ca4CaMCa4CaN_sl[id];
  // CaMKII reaction fluxes
  double Pi = 1.0 - Pb2_sl[id] - Pb_sl[id] - Pt_sl[id] - Pt2_sl[id] - Pa_sl[id];
  double rcnCKib2 = kib2*Ca2CaM_sl[id]*Pi - kb2i*Pb2_sl[id];
  double rcnCKb2b = kb24*pow(cs[id],2)*Pb2_sl[id] - kb42*Pb_sl[id];
  double rcnCKib = kib*Ca4CaM_sl[id]*Pi - kbi*Pb_sl[id];
  double T = Pb_sl[id] + Pt_sl[id] + Pt2_sl[id] + Pa_sl[id];
  double kbt = 0.055*T + .0074*pow(T,2) + 0.015*pow(T,3);
  double rcnCKbt = kbt*Pb_sl[id] - kpp1*PP1totSL*Pt_sl[id]/(Kmpp1+CaMKIItotSL*Pt_sl[id]);
  double rcnCKtt2 = kt42*Pt_sl[id] - kt24*pow(cs[id],2)*Pt2_sl[id];
  double rcnCKta = kta*Pt_sl[id] - kat*Ca4CaM_sl[id]*Pa_sl[id];
  double rcnCKt2a = kt2a*Pt2_sl[id] - kat2*Ca2CaM_sl[id]*Pa_sl[id];
  double rcnCKt2b2 = kpp1*PP1totSL*Pt2_sl[id]/(Kmpp1+CaMKIItotSL*Pt2_sl[id]);
  double rcnCKai = kpp1*PP1totSL*Pa_sl[id]/(Kmpp1+CaMKIItotSL*Pa_sl[id]);

  // CaM_sl equations
  double dCaM = 1e-3*(-rcn02 - rcn0B - rcn0CaN);
  double dCa2CaM = 1e-3*(rcn02 - rcn24 - rcn2B - rcn2CaN + CaMKIItotSL*(-rcnCKib2 + rcnCKt2a) );
  double dCa4CaM = 1e-3*(rcn24 - rcn4B - rcn4CaN + CaMKIItotSL*(-rcnCKib+rcnCKta) );
  double dCaMB = 1e-3*(rcn0B-rcn02B);
  double dCa2CaMB = 1e-3*(rcn02B + rcn2B - rcn24B);
  double dCa4CaMB = 1e-3*(rcn24B + rcn4B);

  // CaMKII equations
  double dPb2 = 1e-3*(rcnCKib2 - rcnCKb2b + rcnCKt2b2); // Pb2_sl[id]
  double dPb = 1e-3*(rcnCKib + rcnCKb2b - rcnCKbt);    // Pb_sl[id]
  double dPt = 1e-3*(rcnCKbt-rcnCKta-rcnCKtt2);        // Pt_sl[id]
  double dPt2 = 1e-3*(rcnCKtt2-rcnCKt2a-rcnCKt2b2);     // Pt2_sl[id]
  double dPa = 1e-3*(rcnCKta+rcnCKt2a-rcnCKai);       // Pa_sl[id]

  // CaN equations
  double dCa4CaN = 1e-3*(rcnCa4CaN - rcn0CaN - rcn2CaN - rcn4CaN);                       // Ca4CaN_sl[id]
  double dCaMCa4CaN = 1e-3*(rcn0CaN - rcn02CaN);           // CaMCa4CaN_sl[id]
  double dCa2CaMCa4CaN = 1e-3*(rcn2CaN+rcn02CaN-rcn24CaN);    // Ca2CaMCa4CaN_sl[id]
  double dCa4CaMCa4CaN = 1e-3*(rcn4CaN+rcn24CaN);             // Ca4CaMCa4CaN_sl[id]

  // double dydt[15]={dCaM,dCa2CaM,dCa4CaM,dCaMB,dCa2CaMB,dCa4CaMB,dPb2,dPb,dPt,dPt2,dPa,dCa4CaN,dCaMCa4CaN,dCa2CaMCa4CaN,dCa4CaMCa4CaN};
  
  
  // write to global variables for adjusting Ca dyad buffering in EC coupling model
  JCaSL[id] = 1e-3*(2.0*CaMKIItotSL*(rcnCKtt2-rcnCKb2b) - 2.0*(rcn02+rcn24+rcn02B+rcn24B+rcnCa4CaN+rcn02CaN+rcn24CaN)); // [uM/msec]
  
  //Assign dydt
  dydt_CaMSL[(id*15)+  0] = dCaM;
  dydt_CaMSL[(id*15)+  1] = dCa2CaM;
  dydt_CaMSL[(id*15)+  2] = dCa4CaM;
  dydt_CaMSL[(id*15)+  3] = dCaMB;
  dydt_CaMSL[(id*15)+  4] = dCa2CaMB;
  dydt_CaMSL[(id*15)+  5] = dCa4CaMB;
  dydt_CaMSL[(id*15)+  6] = dPb2;
  dydt_CaMSL[(id*15)+  7] = dPb;
  dydt_CaMSL[(id*15)+  8] = dPt;
  dydt_CaMSL[(id*15)+  9] = dPt2;
  dydt_CaMSL[(id*15)+ 10] = dPa;
  dydt_CaMSL[(id*15)+ 11] = dCa4CaN;
  dydt_CaMSL[(id*15)+ 12] = dCaMCa4CaN;
  dydt_CaMSL[(id*15)+ 13] = dCa2CaMCa4CaN;
  dydt_CaMSL[(id*15)+ 14] = dCa4CaMCa4CaN; 
}
void CSubcell::calc_dydt_CaM_Cyt_ODEs(int id){
  double K = 135.0; // mM
  double Mg = 1.0; //mM


  // Negroni et al model - CaM_cyt module

  //// Descriptions for state variables
  // Ca2CaM_cyt[id] = y(1);  // 2 Ca bound to C terminal sites
  // Ca4CaM_cyt[id] = y(2);  // 4 Ca bound
  // CaMB_cyt[id] = y(3);   
  // Ca2CaMB_cyt[id] = y(4);
  // Ca4CaMB_cyt[id] = y(5);
  // Pb2_cyt[id] = y(6);     // probability of a Ca2CaM_cyt[id] bound CaMKII subunit
  // Pb_cyt[id] = y(7);      // probability of a Ca4CaM_cyt[id] bound CaMKII subunit
  // Pt_cyt[id] = y(8);      // probability of a Ca4CaM_cyt[id] bound autophosphorylated CaMKII subunit
  // Pt2_cyt[id] = y(9);     // probability of a Ca2CaMxx bound autophosphorylated CaMKII subunit
  // Pa_cyt[id] = y(10);     // probability of an autonomous autophosphorylated CaMKII subunit
  // CaMCaN = y(11);
  // Ca2CaMCaN = y(12);
  // Ca4CaMCaN = y(13);  // active calcineurin
  // description of intermediate variables
  // CaM_cyt- Ca free CaM_cyt
  ////

  // Parameters
  // Ca/CaM_cyt parameters
  double Kd02, Kd24;
  if(Mg <= 1.0){
      Kd02 = 0.0025*(1.0+K/0.94-Mg/0.012)*(1.0+K/8.1+Mg/0.022);  // [uM^2]
      Kd24 = 0.128*(1.0+K/0.64+Mg/0.0014)*(1.0+K/13.0-Mg/0.153); // [uM^2]
  }
  else{
      Kd02 = 0.0025*(1.0+K/0.94-1.0/0.012+(Mg-1.0)/0.060)*(1.0+K/8.1+1.0/0.022+(Mg-1.0)/0.068);   // [uM^2]
      Kd24 = 0.128*(1.0+K/0.64+1.0/0.0014+(Mg-1.0)/0.005)*(1.0+K/13.0-1.0/0.153+(Mg-1.0)/0.150);  // [uM^2]
  }

  double k20 = 10.0;               // [s^-1]      
  double k02 = k20/Kd02;         // [uM^-2 s^-1]
  double k42 = 500.0;              // [s^-1]      
  double k24 = k42/Kd24;         // [uM^-2 s^-1]

  // CaM_cyt buffering (B) parameters
  double k0Boff = 0.0014;        // [s^-1] 
  double k0Bon = k0Boff/0.2;   // [uM^-1 s^-1] kon = koff/Kd
  double k2Boff = k0Boff/100.0;    // [s^-1] 
  double k2Bon = k0Bon;          // [uM^-1 s^-1]
  double k4Boff = k2Boff;        // [s^-1]
  double k4Bon = k0Bon;          // [uM^-1 s^-1]
  // using thermodynamic constraints
  double k20B = k20/100.0; // [s^-1] thermo constraint on loop 1
  double k02B = k02;     // [uM^-2 s^-1] 
  double k42B = k42;     // [s^-1] thermo constraint on loop 2
  double k24B = k24;     // [uM^-2 s^-1]

  // CaMKII parameters
  // Wi Wa Wt Wp
  double kbi = 2.2;      // [s^-1] (Ca4CaM_cyt[id] dissocation from Wb)
  double kib = kbi/33.5e-3; // [uM^-1 s^-1]
  double kib2 = kib;
  double kb2i = kib2*5.0;
  double kb24 = k24;
  double kb42 = k42*33.5e-3/5;
  double kpp1 = 1.72;    // [s^-1] (PP1-dep dephosphorylation rates)
  double Kmpp1 = 11.5;   // [uM]
  double kta = kbi/1000; // [s^-1] (Ca4CaM_cyt[id] dissociation from Wt)
  double kat = kib;      // [uM^-1 s^-1] (Ca4CaM_cyt[id] reassociation with Wa)
  double kt42 = k42*33.5e-6/5.0;
  double kt24 = k24;
  double kat2 = kib;
  double kt2a = kib*5.0;

  // CaN parameters
  double kcanCaoff = 1.0;              // [s^-1] 
  double kcanCaon = kcanCaoff/0.5;   // [uM^-1 s^-1] 
  double kcanCaM4on = 46.0;            // [uM^-1 s^-1]
  double kcanCaM4off = 1.3e-3;       // [s^-1]
  double kcanCaM2on = kcanCaM4on;
  double kcanCaM2off = 2508.0*kcanCaM4off;
  double kcanCaM0on = kcanCaM4on;
  double kcanCaM0off = 165.0*kcanCaM2off;
  double k02can = k02;
  double k20can = k20/165.0;
  double k24can = k24;
  double k42can = k20/2508.0;

  // CaM_cyt Reaction fluxes
  double rcn02 = k02*pow(ci[id],2)*CaM_cyt[id] - k20*Ca2CaM_cyt[id];
  double rcn24 = k24*pow(ci[id],2)*Ca2CaM_cyt[id] - k42*Ca4CaM_cyt[id];
  // CaM_cyt buffer fluxes
  double B = BtotCyt - CaMB_cyt[id] - Ca2CaMB_cyt[id] - Ca4CaMB_cyt[id];
  double rcn02B = k02B*pow(ci[id],2)*CaMB_cyt[id] - k20B*Ca2CaMB_cyt[id];
  double rcn24B = k24B*pow(ci[id],2)*Ca2CaMB_cyt[id] - k42B*Ca4CaMB_cyt[id];
  double rcn0B = k0Bon*CaM_cyt[id]*B - k0Boff*CaMB_cyt[id];
  double rcn2B = k2Bon*Ca2CaM_cyt[id]*B - k2Boff*Ca2CaMB_cyt[id];
  double rcn4B = k4Bon*Ca4CaM_cyt[id]*B - k4Boff*Ca4CaMB_cyt[id];
  // CaN reaction fluxes 
  double Ca2CaN = CaNtotCyt - Ca4CaN_cyt[id] - CaMCa4CaN_cyt[id] - Ca2CaMCa4CaN_cyt[id] - Ca4CaMCa4CaN_cyt[id];
  double rcnCa4CaN = kcanCaon*pow(ci[id],2)*Ca2CaN - kcanCaoff*Ca4CaN_cyt[id];
  double rcn02CaN = k02can*pow(ci[id],2)*CaMCa4CaN_cyt[id] - k20can*Ca2CaMCa4CaN_cyt[id]; 
  double rcn24CaN = k24can*pow(ci[id],2)*Ca2CaMCa4CaN_cyt[id] - k42can*Ca4CaMCa4CaN_cyt[id];
  double rcn0CaN = kcanCaM0on*CaM_cyt[id]*Ca4CaN_cyt[id] - kcanCaM0off*CaMCa4CaN_cyt[id];
  double rcn2CaN = kcanCaM2on*Ca2CaM_cyt[id]*Ca4CaN_cyt[id] - kcanCaM2off*Ca2CaMCa4CaN_cyt[id];
  double rcn4CaN = kcanCaM4on*Ca4CaM_cyt[id]*Ca4CaN_cyt[id] - kcanCaM4off*Ca4CaMCa4CaN_cyt[id];
  // CaMKII reaction fluxes
  double Pi = 1.0 - Pb2_cyt[id] - Pb_cyt[id] - Pt_cyt[id] - Pt2_cyt[id] - Pa_cyt[id];
  double rcnCKib2 = kib2*Ca2CaM_cyt[id]*Pi - kb2i*Pb2_cyt[id];
  double rcnCKb2b = kb24*pow(ci[id],2)*Pb2_cyt[id] - kb42*Pb_cyt[id];
  double rcnCKib = kib*Ca4CaM_cyt[id]*Pi - kbi*Pb_cyt[id];
  double T = Pb_cyt[id] + Pt_cyt[id] + Pt2_cyt[id] + Pa_cyt[id];
  double kbt = 0.055*T + .0074*pow(T,2) + 0.015*pow(T,3);
  double rcnCKbt = kbt*Pb_cyt[id] - kpp1*PP1totCyt*Pt_cyt[id]/(Kmpp1+CaMKIItotCyt*Pt_cyt[id]);
  double rcnCKtt2 = kt42*Pt_cyt[id] - kt24*pow(ci[id],2)*Pt2_cyt[id];
  double rcnCKta = kta*Pt_cyt[id] - kat*Ca4CaM_cyt[id]*Pa_cyt[id];
  double rcnCKt2a = kt2a*Pt2_cyt[id] - kat2*Ca2CaM_cyt[id]*Pa_cyt[id];
  double rcnCKt2b2 = kpp1*PP1totCyt*Pt2_cyt[id]/(Kmpp1+CaMKIItotCyt*Pt2_cyt[id]);
  double rcnCKai = kpp1*PP1totCyt*Pa_cyt[id]/(Kmpp1+CaMKIItotCyt*Pa_cyt[id]);

  // CaM_cyt equations
  double dCaM = 1e-3*(-rcn02 - rcn0B - rcn0CaN);
  double dCa2CaM = 1e-3*(rcn02 - rcn24 - rcn2B - rcn2CaN + CaMKIItotCyt*(-rcnCKib2 + rcnCKt2a) );
  double dCa4CaM = 1e-3*(rcn24 - rcn4B - rcn4CaN + CaMKIItotCyt*(-rcnCKib+rcnCKta) );
  double dCaMB = 1e-3*(rcn0B-rcn02B);
  double dCa2CaMB = 1e-3*(rcn02B + rcn2B - rcn24B);
  double dCa4CaMB = 1e-3*(rcn24B + rcn4B);

  // CaMKII equations
  double dPb2 = 1e-3*(rcnCKib2 - rcnCKb2b + rcnCKt2b2); // Pb2_cyt[id]
  double dPb = 1e-3*(rcnCKib + rcnCKb2b - rcnCKbt);    // Pb_cyt[id]
  double dPt = 1e-3*(rcnCKbt-rcnCKta-rcnCKtt2);        // Pt_cyt[id]
  double dPt2 = 1e-3*(rcnCKtt2-rcnCKt2a-rcnCKt2b2);     // Pt2_cyt[id]
  double dPa = 1e-3*(rcnCKta+rcnCKt2a-rcnCKai);       // Pa_cyt[id]



  // CaN equations
  double dCa4CaN = 1e-3*(rcnCa4CaN - rcn0CaN - rcn2CaN - rcn4CaN);                       // Ca4CaN_cyt[id]
  double dCaMCa4CaN = 1e-3*(rcn0CaN - rcn02CaN);           // CaMCa4CaN_cyt[id]
  double dCa2CaMCa4CaN = 1e-3*(rcn2CaN+rcn02CaN-rcn24CaN);    // Ca2CaMCa4CaN_cyt[id]
  double dCa4CaMCa4CaN = 1e-3*(rcn4CaN+rcn24CaN);             // Ca4CaMCa4CaN_cyt[id]

  // double dydt[15]={dCaM,dCa2CaM,dCa4CaM,dCaMB,dCa2CaMB,dCa4CaMB,dPb2,dPb,dPt,dPt2,dPa,dCa4CaN,dCaMCa4CaN,dCa2CaMCa4CaN,dCa4CaMCa4CaN};
  
  //Assign dydt

  // for(int i=0;i<15;i++){
  //   dydt_CaMCyt[(id*15)+i]=dydt[i];
  // }

  // write to global variables for adjusting Ca dyad buffering in EC coupling model
  JCaCyt[id] = 1e-3*(2.0*CaMKIItotCyt*(rcnCKtt2-rcnCKb2b) - 2.0*(rcn02+rcn24+rcn02B+rcn24B+rcnCa4CaN+rcn02CaN+rcn24CaN)); // [uM/msec]

  dydt_CaMCyt[(id*15)+  0] = dCaM;
  dydt_CaMCyt[(id*15)+  1] = dCa2CaM;
  dydt_CaMCyt[(id*15)+  2] = dCa4CaM;
  dydt_CaMCyt[(id*15)+  3] = dCaMB;
  dydt_CaMCyt[(id*15)+  4] = dCa2CaMB;
  dydt_CaMCyt[(id*15)+  5] = dCa4CaMB;
  dydt_CaMCyt[(id*15)+  6] = dPb2;
  dydt_CaMCyt[(id*15)+  7] = dPb;
  dydt_CaMCyt[(id*15)+  8] = dPt;
  dydt_CaMCyt[(id*15)+  9] = dPt2;
  dydt_CaMCyt[(id*15)+ 10] = dPa;
  dydt_CaMCyt[(id*15)+ 11] = dCa4CaN;
  dydt_CaMCyt[(id*15)+ 12] = dCaMCa4CaN;
  dydt_CaMCyt[(id*15)+ 13] = dCa2CaMCa4CaN;
  dydt_CaMCyt[(id*15)+ 14] = dCa4CaMCa4CaN;
}
void CSubcell::solve_ODE_CaM(int id, double dtx){
  // double dydt[]={dCaM,dCa2CaM,dCa4CaM,dCaMB,dCa2CaMB,dCa4CaMB,dPb2,dPb,dPt,dPt2,dPa,dCa4CaN,dCaMCa4CaN,dCa2CaMCa4CaN,dCa4CaMCa4CaN};
  
  //CaM Cytosol
    CaM_cyt[id]           += dydt_CaMCyt[(id*15) +  0]*dtx;
    Ca2CaM_cyt[id]        += dydt_CaMCyt[(id*15) +  1]*dtx;
    Ca4CaM_cyt[id]        += dydt_CaMCyt[(id*15) +  2]*dtx;
    CaMB_cyt[id]          += dydt_CaMCyt[(id*15) +  3]*dtx;
    Ca2CaMB_cyt[id]       += dydt_CaMCyt[(id*15) +  4]*dtx;
    Ca4CaMB_cyt[id]       += dydt_CaMCyt[(id*15) +  5]*dtx;
    Pb2_cyt[id]           += dydt_CaMCyt[(id*15) +  6]*dtx;
    Pb_cyt[id]            += dydt_CaMCyt[(id*15) +  7]*dtx;
    Pt_cyt[id]            += dydt_CaMCyt[(id*15) +  8]*dtx;
    Pt2_cyt[id]           += dydt_CaMCyt[(id*15) +  9]*dtx;
    Pa_cyt[id]            += dydt_CaMCyt[(id*15) + 10]*dtx;
    Ca4CaN_cyt[id]        += dydt_CaMCyt[(id*15) + 11]*dtx;
    CaMCa4CaN_cyt[id]     += dydt_CaMCyt[(id*15) + 12]*dtx;
    Ca2CaMCa4CaN_cyt[id]  += dydt_CaMCyt[(id*15) + 13]*dtx;
    Ca4CaMCa4CaN_cyt[id]  += dydt_CaMCyt[(id*15) + 14]*dtx;
  //CaM Sarcolemmal

    CaM_sl[id]            += dydt_CaMSL[(id*15) +  0]*dtx;
    Ca2CaM_sl[id]         += dydt_CaMSL[(id*15) +  1]*dtx;
    Ca4CaM_sl[id]         += dydt_CaMSL[(id*15) +  2]*dtx;
    CaMB_sl[id]           += dydt_CaMSL[(id*15) +  3]*dtx;
    Ca2CaMB_sl[id]        += dydt_CaMSL[(id*15) +  4]*dtx;
    Ca4CaMB_sl[id]        += dydt_CaMSL[(id*15) +  5]*dtx;
    Pb2_sl[id]            += dydt_CaMSL[(id*15) +  6]*dtx;
    Pb_sl[id]             += dydt_CaMSL[(id*15) +  7]*dtx;
    Pt_sl[id]             += dydt_CaMSL[(id*15) +  8]*dtx;
    Pt2_sl[id]            += dydt_CaMSL[(id*15) +  9]*dtx;
    Pa_sl[id]             += dydt_CaMSL[(id*15) + 10]*dtx;
    Ca4CaN_sl[id]         += dydt_CaMSL[(id*15) + 11]*dtx;
    CaMCa4CaN_sl[id]      += dydt_CaMSL[(id*15) + 12]*dtx;
    Ca2CaMCa4CaN_sl[id]   += dydt_CaMSL[(id*15) + 13]*dtx;
    Ca4CaMCa4CaN_sl[id]   += dydt_CaMSL[(id*15) + 14]*dtx;
  
  //CaM Dyad/Cleft space
    CaM_dyad[id]          += dydt_CaMDyad[(id*15)+  0]*dtx;
    Ca2CaM_dyad[id]       += dydt_CaMDyad[(id*15)+  1]*dtx;
    Ca4CaM_dyad[id]       += dydt_CaMDyad[(id*15)+  2]*dtx;
    CaMB_dyad[id]         += dydt_CaMDyad[(id*15)+  3]*dtx;
    Ca2CaMB_dyad[id]      += dydt_CaMDyad[(id*15)+  4]*dtx;
    Ca4CaMB_dyad[id]      += dydt_CaMDyad[(id*15)+  5]*dtx;
    Pb2_dyad[id]          += dydt_CaMDyad[(id*15)+  6]*dtx;
    Pb_dyad[id]           += dydt_CaMDyad[(id*15)+  7]*dtx;
    Pt_dyad[id]           += dydt_CaMDyad[(id*15)+  8]*dtx;
    Pt2_dyad[id]          += dydt_CaMDyad[(id*15)+  9]*dtx;
    Pa_dyad[id]           += dydt_CaMDyad[(id*15)+ 10]*dtx;
    Ca4CaN_dyad[id]       += dydt_CaMDyad[(id*15)+ 11]*dtx;
    CaMCa4CaN_dyad[id]    += dydt_CaMDyad[(id*15)+ 12]*dtx;
    Ca2CaMCa4CaN_dyad[id] += dydt_CaMDyad[(id*15)+ 13]*dtx;
    Ca4CaMCa4CaN_dyad[id] += dydt_CaMDyad[(id*15)+ 14]*dtx;
}


void CSubcell::calc_dydt_CaMKII_ODEs(int id)
  {
  //// Description of state variables
  // LCCp_PKA = y(1);        // [LCCp] by PKA (currently unused anywhere else)
  // LCCp-CaMKIIdyad = y(2); // Dyadic [LCCp] by dyadic CaMKII
  // RyR-Ser2809p = y(3);    // [RyR-Ser2809p] by PKA (currently unused anywhere else)
  // RyR-Ser2815p = y(4);    // [RyR-Ser2815p] by CaMKII 
  // PLB-Thr17p = y(5);      // [PLB-Thr17p] by CaMKII
  // LCCp-CaMKIIsl = y(6);   // Subsarcolemmal [LCCp] by subsarcolemmal CaMKII
  //// RATE CONSTANTS and KM VALUES
  // L-Type Ca Channel (LCC) parameters
  double k_ckLCC = 0.4;                  // [s^-1]
  double k_pp1LCC = 0.1103;              // [s^-1] 
  double k_pkaLCC = 13.5;                // [s^-1] 
  double k_pp2aLCC = 10.1;               // [s^-1] 

  double KmCK_LCC = 12;                  // [uM] 
  double KmPKA_LCC = 21;                 // [uM] 
  double KmPP2A_LCC = 47;                // [uM] 
  double KmPP1_LCC = 9;                  // [uM] 

  // Ryanodine Receptor (RyR) parameters
  double k_ckRyR = 0.4;                  // [s^-1] 
  double k_pkaRyR = 1.35;                // [s^-1] 
  double k_pp1RyR = 1.07;                // [s^-1] 
  double k_pp2aRyR = 0.481;              // [s^-1] 

  // Basal RyR phosphorylation (numbers based on param estimation)
  double kb_2809 = 0.51;                 // [uM/s] - PKA site
  double kb_2815 = 0.35;                 // [uM/s] - CaMKII site

  double KmCK_RyR = 12;                  // [uM] 
  double KmPKA_RyR = 21;                 // [uM] 
  double KmPP1_RyR = 9;                  // [uM] 
  double KmPP2A_RyR = 47;                // [uM] 

  // Phospholamban (PLB) parameters
  double k_ckPLB = 8e-3;                 // [s^-1]
  double k_pp1PLB = .0428;               // [s^-1]

  double KmCK_PLB = 12;
  double KmPP1_PLB = 9;

  // Okadaic Acid inhibition params (based on Huke/Bers [2008])
  // Want to treat OA as non-competitive inhibitor of PP1 and PP2A
  double Ki_OA_PP1 = 0.78;        // [uM] - Values from fit
  double Ki_OA_PP2A = 0.037;      // [uM] - Values from fit

  // Default PKA level
  double PKAc = 95.6*.54;

  //// OA inhibition term (non-competitive) for PP1 and PP2A
  double OA_PP1 = 1/(1 + pow(OA/Ki_OA_PP1,3));
  double OA_PP2A = 1/(1 + pow(OA/Ki_OA_PP2A,3));

  //// ODE EQUATIONS
  //// LCC states (note: PP2A is acting on PKA site and PP1 on CKII site)
  // CaMKII phosphorylation of Dyadic LCCs
  double LCC_CKdyadn = LCCtotDyad - LCC_CKdyadp[id];
  double LCCDyad_PHOS = (k_ckLCC*CaMKIIactDyad[id]*LCC_CKdyadn)/(KmCK_LCC+LCC_CKdyadn);
  double LCCDyad_DEPHOS = (k_pp1LCC*PP1_dyad*LCC_CKdyadp[id])/(KmPP1_LCC+LCC_CKdyadp[id])*OA_PP1;
  double dLCC_CKdyadp = LCCDyad_PHOS - LCCDyad_DEPHOS;

  // CaMKII phosphorylation of Sub-sarcolemmal LCCs
  double LCC_CKsln = LCCtotSL - LCC_CKslp[id];
  double LCCSL_PHOS = (k_ckLCC*CaMKIIactSL[id]*LCC_CKsln)/(KmCK_LCC+LCC_CKsln); 
  double LCCSL_DEPHOS = (k_pp1LCC*PP1_SL*LCC_CKslp[id])/(KmPP1_LCC+LCC_CKslp[id])*OA_PP1;
  double dLCC_CKslp = LCCSL_PHOS - LCCSL_DEPHOS; 

  // PKA phosphorylation (currently unused elsewhere, obsolete with B-AR module)
  double LCC_PKAn = LCCtotDyad - LCC_PKAp[id];
  double dLCC_PKAp = (k_pkaLCC*PKAc*LCC_PKAn)/(KmPKA_LCC+LCC_PKAn) - 
              (k_pp2aLCC*PP2A_dyad*LCC_PKAp[id])/(KmPP2A_LCC+LCC_PKAp[id])*OA_PP2A;
  //// RyR states
  double RyR2815n = RyRtot - RyR2815p[id];
  double RyR_BASAL = kb_2815*RyR2815n;
  double RyR_PHOS = (k_ckRyR*CaMKIIactDyad[id]*RyR2815n)/(KmCK_RyR+RyR2815n);
  double RyR_PP1_DEPHOS = (k_pp1RyR*PP1_dyad*RyR2815p[id])/(KmPP1_RyR+RyR2815p[id])*OA_PP1;
  double RyR_PP2A_DEPHOS = (k_pp2aRyR*PP2A_dyad*RyR2815p[id])/(KmPP2A_RyR+RyR2815p[id])*OA_PP2A;
  double dRyR2815p = RyR_BASAL + RyR_PHOS - RyR_PP1_DEPHOS - RyR_PP2A_DEPHOS;

  // PKA phosphorylation of Ser 2809 on RyR (currently unused elsewhere)
  double RyR2809n = RyRtot - RyR2809p[id];
  double dRyR2809p = kb_2809*RyR2809n + (k_pkaRyR*PKAc*RyR2809n)/(KmPKA_RyR+RyR2809n) - 
              (k_pp1RyR*PP1_dyad*RyR2809p[id])/(KmPP1_RyR+RyR2809p[id])*OA_PP1;        
  //// PLB states
  double PP1_PLB = PP1_dyad*PP1_PLB_avail[id];    // Inhibitor-1 regulation of PP1_dyad included here
  double PLBT17n = PLBtot - PLBT17p[id];
  double PLB_PHOS = (k_ckPLB*PLBT17n*CaMKIIactDyad[id])/(KmCK_PLB+PLBT17n);
  double PLB_DEPHOS = (k_pp1PLB*PP1_PLB*PLBT17p[id])/(KmPP1_PLB+PLBT17p[id])*OA_PP1;
  double dPLBT17p = PLB_PHOS - PLB_DEPHOS; 

  if(id == 0){
   printf(
      "PLB_PHOS: %f\tk_ckPLB: %f\tPLBT17n: %f\t CaMKIIactDyad[id]: %f\tKmCK_PLB: %f\tk_ckPLB[id]: %f\tPLB_DEPHOS: %f\tk_pp1PLB: %f\tPP1_PLB: %f\tPLBT17p[id]: %f\tKmPP1_PLB: %f\tOA_PP1: %f\t\n", 
      PLB_PHOS,  
      k_ckPLB, 
      PLBT17n,
      CaMKIIactDyad[id],
      KmCK_PLB,
      k_ckPLB,
      PLB_DEPHOS,
      k_pp1PLB,
      PP1_PLB,
      PLBT17p[id],
      KmPP1_PLB,
      OA_PP1);

  }
  
  //// Collect ODEs and convert to uM/ms
  dydt_CaMKII[(id*6) + 0] = dLCC_PKAp*1e-3;
  dydt_CaMKII[(id*6) + 1] = dLCC_CKdyadp*1e-3;
  dydt_CaMKII[(id*6) + 2] = dRyR2809p*1e-3;
  dydt_CaMKII[(id*6) + 3] = dRyR2815p*1e-3;
  dydt_CaMKII[(id*6) + 4] = dPLBT17p*1e-3;
  dydt_CaMKII[(id*6) + 5] = dLCC_CKslp*1e-3;

}
void CSubcell::solve_ODE_CaMKII(int id, double dtx){
  LCC_PKAp[id]    += dydt_CaMKII[(id*6)+0]*dtx; //unused
  LCC_CKdyadp[id] += dydt_CaMKII[(id*6)+1]*dtx;
  RyR2809p[id]    += dydt_CaMKII[(id*6)+2]*dtx; //unused
  RyR2815p[id]    += dydt_CaMKII[(id*6)+3]*dtx;
  PLBT17p[id]     += dydt_CaMKII[(id*6)+4]*dtx;
  LCC_CKslp[id]   += dydt_CaMKII[(id*6)+5]*dtx;
}

void CSubcell::calc_dydt_BAR_ODEs(int id){
  // Negroni et al model - PKA module

  //Matlab State Variables 
  // y(1) = L;
  // y(2) = B1AR;
  // y(3) = Gs;
  // y(4) = B1AR_ACT;
  // y(5) = B1AR_S464;
  // y(6) = B1AR_S301;
  // y(7) = GsaGTPtot;
  // y(8) = GsaGDP;
  // y(9) = GsBy;
  // y(10) = GsaGTP;
  // y(11) = Fsk;
  // y(12) = AC;
  // y(13) = PDE;
  // y(14) = IBMX;
  // y(15) = cAMPtot;
  // y(16) = cAMP;
  // y(17) = PKAC_I;
  // y(18) = PKAC_II;
  // y(19) = PLBp;
  // y(20) = Inhib1ptot;
  // y(21) = Inhib1p;
  // y(22) = PP1;
  // y(23) = LCCa_PKAp_whole;
  // y(24) = LCCb_PKAp_whole;
  // y(25) = RyR_PKAp_whole;
  // y(26) = TnI_PKAp_whole;
  // y(27) = IKs_PKAn;
  // y(28) = Yotiao_KCQN1;
  // y(29) = IKs_PKAp_whole;
  // y(30) = ICFTR_PKAp_whole;
  // y(31) = PLM_PKAp_whole;
  // y(32) = Myo_PKAp_whole;
  // y(33) = IKr_PKAn;
  // y(34) = Yotiao_hERG;
  // y(35) = IKr_PKAp_whole;
  // y(36) = IClCa_PKAp_whole;


  //// Parameters
  //// ----- Signaling model parameters -------
  // b-AR/Gs module
  double Ltotmax = Ligtot;    // Ltotmax   [uM] ** apply agonist concentration here **
  double sumb1AR = 0.028;   // sumb1AR   [uM]
  double Gstot = 3.83;    // Gstot     [uM]
  double Kl = 0.285;   // Kl        [uM]
  double Kr = 0.062;   // Kr        [uM]
  double Kc = 33.0;    // Kc        [uM]
  double k_barkp = 1.1e-3;  // k_barkp   [1/sec]
  double k_barkm = 2.2e-3;  // k_barkm   [1/sec]
  double k_pkap = 3.6e-3;  // k_pkap    [1/sec/uM]
  double k_pkam = 2.2e-3; // k_pkam    [1/sec]
  double k_gact = 16.0;   // k_gact    [1/sec]
  double k_hyd = 0.8;    // k_hyd     [1/sec]
  double k_reassoc = 1.21e3; // k_reassoc [1/sec/uM]
  // cAMP module
  double AC_tot = 0.047;  // AC_tot    [uM]
  double ATP = 5.0e3;  // ATP       [uM]
  double PDE3tot = 0.036;  // PDE3tot   [uM] // Changed from .06 to .036
  double PDE4tot = 0.036;  // PDE4tot   [uM]
  double IBMXtot = 0.0;    // IBMXtot   [uM]
  double Fsktot = 0.0;    // Fsktot    [uM] (10 uM when used)
  double k_ac_basal = 0.2;    // k_ac_basal[1/sec]
  double k_ac_gsa = 8.5;    // k_ac_gsa  [1/sec]
  double k_ac_fsk = 7.3;    // k_ac_fsk  [1/sec]
  double Km_basal = 1.03e3; // Km_basal  [uM]
  double Km_gsa = 315.0;  // Km_gsa    [uM]
  double Km_fsk = 860.0;  // Km_fsk    [uM]
  double Kgsa = 0.4;    // Kgsa      [uM]
  double Kfsk = 44.0;   // Kfsk      [uM]
  double k_pde3 = 3.5;    // k_pde3    [1/sec]
  double Km_pde3 = 0.15;   // Km_pde3   [uM]
  double k_pde4 = 5.0;    // k_pde4    [1/sec]
  double Km_pde4 = 1.3;    // Km_pde4   [uM]
  double Ki_ibmx = 30.0;   // Ki_ibmx   [uM]
  // PKA module
  double PKAItot = 0.46;   // PKAItot   [uM]
  double PKAIItot = 0.084;  // PKAIItot  [uM]
  double PKItot = 0.18;   // PKItot    [uM]
  double Ka = 9.14;   // Ka        [uM]
  double Kb = 1.64;   // Kb        [uM]
  double Kd = 4.375;  // Kd        [uM]
  double Ki_pki = 0.2e-3; // Ki_pki    [uM]
  // PLB & PP1 module
  double epsilon = 10;     // epsilon       [none]
  double PLBtot = PLBtotBA; // PLBtot        [uM]
  double PP1tot = PP1_PLBtot; // PP1tot    [uM]
  double Inhib1tot = 0.3;    // Inhib1tot     [uM]
  double k_pka_plb = 54;     // k_pka_plb     [1/sec]
  double Km_pka_plb = 21;     // Km_pka_plb    [uM]
  double k_pp1_plb = 8.5;    // k_pp1_plb     [1/sec]
  double Km_pp1_plb = 7.0;    // Km_pp1_plb    [uM]
  double k_pka_i1 = 60;     // k_pka_i1      [1/sec]
  double Km_pka_i1 = 1.0;    // Km_pka_i1     [uM]
  double Vmax_pp2a_i1 = 14.0;   // Vmax_pp2a_i1  [uM/sec]
  double Km_pp2a_i1 = 1.0;    // Km_pp2a_i1    [uM]
  double Ki_inhib1 = 1.0e-3; // Ki_inhib1     [uM]
  // LCC module
  double LCCtot = LCCtotBA; // LCCtot        [uM]
  double PKAIIlcctot = 0.025;  // PKAIIlcctot   [uM]
  double PP1lcctot = 0.025;  // PP1lcctot     [uM]
  double PP2Alcctot = 0.025;  // PP2Alcctot    [uM]
  double k_pka_lcc = 54;     // k_pka_lcc     [1/sec]
  double Km_pka_lcc = 21;     // Km_pka_lcc    [uM]
  double k_pp1_lcc = 8.52;   // k_pp1_lcc     [1/sec]
  double Km_pp1_lcc = 3;      // Km_pp1_lcc    [uM]
  double k_pp2a_lcc = 10.1;   // k_pp2a_lcc    [1/sec]
  double Km_pp2a_lcc = 3;      // Km_pp2a_lcc   [uM]
  // RyR module
  double RyRtot = RyRtotBA; // RyRtot        [uM]
  double PKAIIryrtot = 0.034;  // PKAIIryrtot   [uM]
  double PP1ryr = 0.034;  // PP1ryr        [uM]
  double PP2Aryr = 0.034;  // PP2Aryr       [uM]
  double kcat_pka_ryr = 54;     // kcat_pka_ryr  [1/sec]
  double Km_pka_ryr = 21;     // Km_pka_ryr    [uM]
  double kcat_pp1_ryr = 8.52;   // kcat_pp1_ryr  [1/sec]
  double Km_pp1_ryr = 7;      // Km_pp1_ryr    [uM]
  double kcat_pp2a_ryr = 10.1;   // kcat_pp2a_ryr [1/sec]
  double Km_pp2a_ryr = 4.1;    // Km_pp2a_ryr   [uM]
  // TnI module
  double TnItot = TnItotBA; // TnItot        [uM]
  double PP2Atni = 0.67;   // PP2Atni       [uM]
  double kcat_pka_tni = 54;     // kcat_pka_tni  [1/sec]
  double Km_pka_tni = 21;     // Km_pka_tni    [uM]
  double kcat_pp2a_tni = 10.1;   // kcat_pp2a_tni [1/sec]
  double Km_pp2a_tni = 4.1;    // Km_pp2a_tni   [uM]
  // Iks module
  double Iks_tot = IKstotBA; // Iks_tot       [uM]
  double Yotiao_tot = 0.025;  // Yotiao_tot    [uM]
  double K_yotiao = 0.1e-3; // K_yotiao      [uM] ** apply G589D mutation here **
  double PKAII_ikstot = 0.025;  // PKAII_ikstot  [uM]
  double PP1_ikstot = 0.025;  // PP1_ikstot    [uM]
  double k_pka_iks = 1.87;//54; // k_pka_iks   [1/sec] // adjusted as in Xie et al 2013
  double Km_pka_iks = 21;     // Km_pka_iks    [uM]
  double k_pp1_iks   = 0.19;//8.52; // k_pp1_iks [1/sec] // adjusted as in Xie et al 2013
  double Km_pp1_iks = 7;      // Km_pp1_iks    [uM]
  // Icftr Module - Added 04/30/10 by Anthony Soltis
  double CFTR_tot = ICFTRtotBA; // CFTR_tot    [uM]
  double PKAII_CFTRtot = 0.025;  // PKAII_CFTRtot [uM]
  double PP1_CFTRtot = 0.025;  // PP1_CFTRtot   [uM]
  double k_pka_CFTR = 54;     // k_pka_CFTR    [1/sec]
  double Km_pka_CFTR = 8.5;    // Km_pka_CFTR   [uM]
  double k_pp1_CFTR = 8.52;   // k_pp1_CFTR    [1/sec]
  double Km_pp1_CFTR = 7;      // Km_pp1_CFTR   [uM]
  // PLM module (from PLB)
  double PLMtot = PLMtotBA; // PLM tot       [uM]
  double kcat_pka_plm = 54;     // kcat_pka_plm  [1/sec]
  double Km_pka_plm = 21;     // Km_pka_plm    [uM]
  double kcat_pp2a_plm = 8.5;    // kcat_pp2a_plm [1/sec]
  double Km_pp2a_plm = 7.0;    // Km_pp2a_plm   [uM]
  // Myofilament module (from TnI)
  double Myo_tot = MyototBA; // Myo tot      [uM]
  double PP2Amyo = 0.67;  // PP2Amyo       [uM]
  double kcat_pka_myo = 54;    // kcat_pka_myo  [1/sec]
  double Km_pka_myo = 21;    // Km_pka_myo    [uM]
  double kcat_pp2a_myo = 10.1;  // kcat_pp2a_myo [1/sec]
  double Km_pp2a_myo = 4.1;   // Km_pp2a_myo   [uM]
  // Ikr module (from Iks)
  double Ikr_tot = IKrtotBA; // Iks_tot       [uM]
  // double Yotiao_tot = 0.025;  // Yotiao_tot    [uM]
  // double K_yotiao = 0.1e-3; // K_yotiao      [uM] ** apply G589D mutation here **
  double PKAII_ikrtot = 0.025;  // PKAII_ikstot  [uM]
  double PP1_ikrtot = 0.025;  // PP1_ikstot    [uM]
  double k_pka_ikr = 1.87;//54; // k_pka_iks   [1/sec] // adjusted as in Xie et al 2013
  double Km_pka_ikr = 21;     // Km_pka_iks    [uM]
  double k_pp1_ikr = 0.19;//8.52; // k_pp1_iks [1/sec] // adjusted as in Xie et al 2013
  double Km_pp1_ikr = 7;      // Km_pp1_iks    [uM]
  // Iclca module (from Icftr)
  double Iclca_tot = IClCatotBA; // Iclca_tot   [uM]
  double PKAII_ClCatot = 0.025;  // PKAII_ClCatot [uM]
  double PP1_ClCatot = 0.025;  // PP1_ClCatot   [uM]
  double k_pka_ClCa = 54;     // k_pka_ClCa    [1/sec]
  double Km_pka_ClCa = 8.5;    // Km_pka_ClCa   [uM]
  double k_pp1_ClCa = 8.52;   // k_pp1_ClCa    [1/sec]
  double Km_pp1_ClCa = 7;      // Km_pp1_ClCa   [uM]
  //// -------- SIGNALING MODEL -----------

  //// b-AR module
  double LR = L[id]*B1AR[id]/Kl;
  double LRG = LR*Gs[id]/Kr;
  double RG = B1AR[id]*Gs[id]/Kc;
  double BARKDESENS = k_barkp*(LR+LRG);
  double BARKRESENS = k_barkm*B1AR_S464[id];
  double PKADESENS = k_pkap*PKAC_I[id]*B1AR_ACT[id];  
  double PKARESENS = k_pkam*B1AR_S301[id];
  double GACT = k_gact*(RG+LRG);
  double HYD = k_hyd*GsaGTPtot[id];
  double REASSOC = k_reassoc*GsaGDP[id]*GsBy[id];
  double dL = Ltotmax-LR-LRG-L[id];
  double dB1AR = B1AR_ACT[id]-LR-LRG-RG-B1AR[id];
  double dGs = Gstot-LRG-RG-Gs[id];
  double dB1AR_ACT = (BARKRESENS-BARKDESENS)+(PKARESENS-PKADESENS);
  double dB1AR_S464 = BARKDESENS-BARKRESENS;
  double dB1AR_S301 = PKADESENS-PKARESENS;
  double dGsaGTPtot = GACT-HYD;
  double dGsaGDP = HYD-REASSOC;
  double dGsBy = GACT-REASSOC;
  // end b-AR module

  //// cAMP module
  double Gsa_gtp_AC = GsaGTP[id]*AC[id]/Kgsa;
  double Fsk_AC = Fsk[id]*AC[id]/Kfsk;
  double AC_ACT_BASAL = k_ac_basal*AC[id]*ATP/(Km_basal+ATP);     
  double AC_ACT_GSA = k_ac_gsa*Gsa_gtp_AC*ATP/(Km_gsa+ATP); 
  double AC_ACT_FSK = k_ac_fsk*Fsk_AC*ATP/(Km_fsk+ATP);    
  // PDE3_ACT = k_pde3*PDE[id]*cAMP[id]/(Km_pde3+cAMP[id]); 
  // PDE4_ACT = k_pde4*PDE[id]*cAMP[id]/(Km_pde4+cAMP[id]); 
  double PDE3_ACT = k_pde3*PDE3tot*cAMP[id]/(Km_pde3*(1+IBMXtot/Ki_ibmx)+cAMP[id]); // new PDE3 term w IBMX
  double PDE4_ACT = k_pde4*PDE4tot*cAMP[id]/(Km_pde4*(1+IBMXtot/Ki_ibmx)+cAMP[id]); // new PDE4 term w IBMX
  double PDE_IBMX = PDE[id]*IBMX[id]/Ki_ibmx;
  double dGsaGTP = GsaGTPtot[id]-Gsa_gtp_AC-GsaGTP[id];
  double dFsk = Fsktot-Fsk_AC-Fsk[id];
  double dAC = AC_tot-Gsa_gtp_AC-AC[id];  // note: assumes Fsk = 0.  Change Gsa_gtp_AC to Fsk_AC for Forskolin.
  // dPDE = PDE4tot-PDE_IBMX-PDE[id];
  // dIBMX = IBMXtot-PDE_IBMX-IBMX[id];
  double dPDE = 0;
  double dIBMX = 0;
  double dcAMPtot = AC_ACT_BASAL+AC_ACT_GSA+AC_ACT_FSK-PDE3_ACT-PDE4_ACT;
  // end cAMP module

  //// PKA module
  double PKI = PKItot*Ki_pki/(Ki_pki+PKAC_I[id]+PKAC_II[id]);
  double A2RC_I = (PKAC_I[id]/Kd)*PKAC_I[id]*(1+PKI/Ki_pki);
  double A2R_I = PKAC_I[id]*(1+PKI/Ki_pki);
  double A2RC_II = (PKAC_II[id]/Kd)*PKAC_II[id]*(1+PKI/Ki_pki);
  double A2R_II = PKAC_II[id]*(1+PKI/Ki_pki);
  double ARC_I = (Ka/cAMP[id])*A2RC_I;
  double ARC_II = (Ka/cAMP[id])*A2RC_II;
  double dcAMP = cAMPtot[id]-(ARC_I+2*A2RC_I+2*A2R_I)-(ARC_II+2*A2RC_II+2*A2R_II)-cAMP[id];
  double PKAtemp = Ka*Kb/Kd+Ka*cAMP[id]/Kd+(cAMP[id]*cAMP[id])/Kd;
  double dPKAC_I = 2*PKAItot*(cAMP[id]*cAMP[id])-PKAC_I[id]*(1+PKI/Ki_pki)*(PKAtemp*PKAC_I[id]+(cAMP[id]*cAMP[id]));
  double dPKAC_II = 2*PKAIItot*(cAMP[id]*cAMP[id])-PKAC_II[id]*(1+PKI/Ki_pki)*(PKAtemp*PKAC_II[id]+(cAMP[id]*cAMP[id]));
  // end PKA module

  //// PLB & PP1 module
  double PLB = PLBtot-PLBp[id];
  double PLB_PHOSPH = k_pka_plb*PKAC_I[id]*PLB/(Km_pka_plb+PLB);
  double PLB_DEPHOSPH = k_pp1_plb*PP1[id]*PLBp[id]/(Km_pp1_plb+PLBp[id]);
  double dPLBp = PLB_PHOSPH-PLB_DEPHOSPH;
   
  double Inhib1 = Inhib1tot-Inhib1ptot[id];
  double Inhib1p_PP1 = Inhib1p[id]*PP1[id]/Ki_inhib1;
  double Inhib1_PHOSPH = k_pka_i1*PKAC_I[id]*Inhib1/(Km_pka_i1+Inhib1); 
  double Inhib1_DEPHOSPH = Vmax_pp2a_i1*Inhib1ptot[id]/(Km_pp2a_i1+Inhib1ptot[id]);
  double dInhib1ptot = Inhib1_PHOSPH-Inhib1_DEPHOSPH;
  double dInhib1p = Inhib1ptot[id]-Inhib1p_PP1-Inhib1p[id];
  double dPP1 = PP1tot-Inhib1p_PP1-PP1[id];
  // end PLB & PP1 module

  //// LCC module
  double PKAClcc = (PKAIIlcctot/PKAIItot)*PKAC_II[id];
  double LCCa = LCCtot-LCCa_PKAp_whole[id];
  double LCCa_PHOSPH = epsilon*k_pka_lcc*PKAClcc*LCCa/(Km_pka_lcc + epsilon*LCCa);
  double LCCa_DEPHOSPH = epsilon*k_pp2a_lcc*PP2Alcctot*LCCa_PKAp_whole[id]/(Km_pp2a_lcc+epsilon*LCCa_PKAp_whole[id]);
  double dLCCa_PKAp_whole = LCCa_PHOSPH - LCCa_DEPHOSPH;
   
  double LCCb = LCCtot-LCCb_PKAp_whole[id];
  double LCCb_PHOSPH = epsilon*k_pka_lcc*PKAClcc*LCCb/(Km_pka_lcc+epsilon*LCCb);   
  double LCCb_DEPHOSPH = epsilon*k_pp1_lcc*PP1lcctot*LCCb_PKAp_whole[id]/(Km_pp1_lcc+epsilon*LCCb_PKAp_whole[id]);
  double dLCCb_PKAp_whole = LCCb_PHOSPH-LCCb_DEPHOSPH;
  // end LCC module

  //// RyR module
  double PKACryr = (PKAIIryrtot/PKAIItot)*PKAC_II[id];
  double RyR = RyRtot-RyR_PKAp_whole[id];
  double RyRPHOSPH = epsilon*kcat_pka_ryr*PKACryr*RyR/(Km_pka_ryr+epsilon*RyR);
  double RyRDEPHOSPH1 = epsilon*kcat_pp1_ryr*PP1ryr*RyR_PKAp_whole[id]/(Km_pp1_ryr+epsilon*RyR_PKAp_whole[id]);
  double RyRDEPHOSPH2A = epsilon*kcat_pp2a_ryr*PP2Aryr*RyR_PKAp_whole[id]/(Km_pp2a_ryr+epsilon*RyR_PKAp_whole[id]);
  double dRyR_PKAp_whole = RyRPHOSPH-RyRDEPHOSPH1-RyRDEPHOSPH2A;
  // end RyR module

  //// TnI module
  double TnI = TnItot-TnI_PKAp_whole[id];
  double TnIPHOSPH = kcat_pka_tni*PKAC_I[id]*TnI/(Km_pka_tni+TnI);
  double TnIDEPHOSPH = kcat_pp2a_tni*PP2Atni*TnI_PKAp_whole[id]/(Km_pp2a_tni+TnI_PKAp_whole[id]);
  double dTnI_PKAp_whole = TnIPHOSPH-TnIDEPHOSPH;
  // end TnI module

  //// Iks module
  double IksYot = IKs_PKAn[id]*Yotiao_KCQN1[id]/K_yotiao;           // [uM]
  double dIKs_PKAn = Iks_tot - IksYot - IKs_PKAn[id];    // [uM]
  double dYotiao_KCQN1 = Yotiao_tot - IksYot - Yotiao_KCQN1[id];    // [uM]
  double PKACiks = (IksYot/Iks_tot)*(PKAII_ikstot/PKAIItot)*PKAC_II[id];
  double PP1iks = (IksYot/Iks_tot)*PP1_ikstot;
  double Iks = Iks_tot-IKs_PKAp_whole[id];
  double IKS_PHOSPH = epsilon*k_pka_iks*PKACiks*Iks/(Km_pka_iks+epsilon*Iks);
  double IKS_DEPHOSPH = epsilon*k_pp1_iks *PP1iks*IKs_PKAp_whole[id]/(Km_pp1_iks+epsilon*IKs_PKAp_whole[id]);
  double dIKs_PKAp_whole = IKS_PHOSPH-IKS_DEPHOSPH;
  // end Iks module

  //// CFTR module (included 04/30/10)
  double CFTRn = CFTR_tot - ICFTR_PKAp_whole[id];  // Non-phos = tot - phos
  double PKAC_CFTR = (PKAII_CFTRtot/PKAIItot)*PKAC_II[id];    // (PKACFTRtot/PKAIItot)*PKAIIact
  double CFTRphos = epsilon*CFTRn*PKAC_CFTR*k_pka_CFTR/(Km_pka_CFTR+epsilon*CFTRn);
  double CFTRdephos = PP1_CFTRtot*k_pp1_CFTR*epsilon*ICFTR_PKAp_whole[id]/(Km_pp1_CFTR+epsilon*ICFTR_PKAp_whole[id]);
  double dICFTR_PKAp_whole = CFTRphos - CFTRdephos;
  // end CFTR module

  //// PLM module (from PLB)
  double PLM = PLMtot-PLM_PKAp_whole[id];
  double PLM_PHOSPH = kcat_pka_plm*PKAC_I[id]*PLM/(Km_pka_plm+PLM);
  double PLM_DEPHOSPH = kcat_pp2a_plm*PP1[id]*PLM_PKAp_whole[id]/(Km_pp2a_plm+PLM_PKAp_whole[id]);
  double dPLM_PKAp_whole = PLM_PHOSPH-PLM_DEPHOSPH;
  // end PLM module

  //// Myofilament module (from TnI)
  double Myo = Myo_tot-Myo_PKAp_whole[id];
  double MyoPHOSPH = kcat_pka_myo*PKAC_I[id]*Myo/(Km_pka_myo+Myo);
  double MyoDEPHOSPH = kcat_pp2a_myo*PP2Amyo*Myo_PKAp_whole[id]/(Km_pp2a_myo+Myo_PKAp_whole[id]);
  double dMyo_PKAp_whole = MyoPHOSPH-MyoDEPHOSPH;
  // end myofilament module

  //// Ikr module (from Iks)
  double IkrYot = IKr_PKAn[id]*Yotiao_hERG[id]/K_yotiao;           // [uM]
  double dIKr_PKAn = Ikr_tot - IkrYot - IKr_PKAn[id];    // [uM]
  double dYotiao_hERG = Yotiao_tot - IkrYot - Yotiao_hERG[id];    // [uM]
  double PKACikr = (IkrYot/Ikr_tot)*(PKAII_ikrtot/PKAIItot)*PKAC_II[id];
  double PP1ikr = (IkrYot/Ikr_tot)*PP1_ikrtot;
  double Ikr = Ikr_tot-IKr_PKAp_whole[id];
  double IKR_PHOSPH = epsilon*k_pka_ikr*PKACikr*Ikr/(Km_pka_ikr+epsilon*Ikr);
  double IKR_DEPHOSPH = epsilon*k_pp1_ikr*PP1ikr*IKr_PKAp_whole[id]/(Km_pp1_ikr+epsilon*IKr_PKAp_whole[id]);
  double dIKr_PKAp_whole = IKR_PHOSPH-IKR_DEPHOSPH;
  // end Ikr module

  //// ICl(Ca) module
  // 88->115 89->116 90-117 91-118 92-119 93-12o 94-121
  double ClCan = CFTR_tot - IClCa_PKAp_whole[id];  // Non-phos = tot - phos
  double PKAC_ClCa = (PKAII_ClCatot/PKAIItot)*PKAC_II[id];    // (PKACFTRtot/PKAIItot)*PKAIIact
  double ClCaphos = epsilon*ClCan*PKAC_ClCa*k_pka_ClCa/(Km_pka_ClCa+epsilon*ClCan);
  double ClCadephos = PP1_ClCatot*k_pp1_ClCa*epsilon*IClCa_PKAp_whole[id]/(Km_pp1_ClCa+epsilon*IClCa_PKAp_whole[id]);
  double dIClCa_PKAp_whole = ClCaphos - ClCadephos;
  // end ICl(Ca) module

  //// ISO-target (set ydot(x) = 0 to prevent ISO effect on a specific target)
  //dPLBp = 0; // PLB
  //dLCCa_PKAp_whole = 0; dLCCb_PKAp_whole = 0; // LCCa and LCCb
  //dRyR_PKAp_whole = 0; // RyR
  // //dTnI_PKAp_whole = 0; // TnI // not used
  //dIKs_PKAp_whole = 0; // IKs
  //dICFTR_PKAp_whole = 0; // ICFTR
  //dPLM_PKAp_whole = 0; // PLM
  //dIKr_PKAp_whole = 0; // IKr
  //dIClCa_PKAp_whole = 0; // IClCa

  //dMyo_PKAp_whole = 0; // Myofilament

  //// Gather odes
  // Need to convert all ydot terms that are ODEs (not DAEs) to miliseconds
  // odes = [4,5,6,7,8,9,13,14,15,19,20,23,24,25,26,29,30,31,32,35,36];
  // ydot(odes) = ydot(odes).*1e-3;
  dB1AR_ACT         *= 1e-3;
  dB1AR_S464        *= 1e-3;
  dB1AR_S301        *= 1e-3;
  dGsaGTPtot        *= 1e-3;
  dGsaGDP           *= 1e-3;
  dGsBy             *= 1e-3;
  dPDE              *= 1e-3;
  dIBMX             *= 1e-3;
  dcAMPtot          *= 1e-3;
  dPLBp             *= 1e-3;
  dInhib1ptot       *= 1e-3;
  dLCCa_PKAp_whole  *= 1e-3;
  dLCCb_PKAp_whole  *= 1e-3;
  dRyR_PKAp_whole   *= 1e-3;
  dTnI_PKAp_whole   *= 1e-3;
  dIKs_PKAp_whole   *= 1e-3;
  dICFTR_PKAp_whole *= 1e-3;
  dPLM_PKAp_whole   *= 1e-3;
  dMyo_PKAp_whole   *= 1e-3;
  dIKr_PKAp_whole   *= 1e-3;
  dIClCa_PKAp_whole *= 1e-3;

  // double dydt[36] = {dL,dB1AR,dGs,dB1AR_ACT,dB1AR_S464,dB1AR_S301,dGsaGTPtot,dGsaGDP,dGsBy,dGsaGTP,dFsk,dAC,dPDE,dIBMX,dcAMPtot,dcAMP,dPKAC_I,dPKAC_II,dPLBp,dInhib1ptot,dInhib1p,dPP1,dLCCa_PKAp_whole,dLCCb_PKAp_whole,dRyR_PKAp_whole,dTnI_PKAp_whole,dIKs_PKAn,dYotiao_KCQN1,dIKs_PKAp_whole,dICFTR_PKAp_whole,dPLM_PKAp_whole,dMyo_PKAp_whole,dIKr_PKAn,dYotiao_hERG,dIKr_PKAp_whole,dIClCa_PKAp_whole};

  // for(int i=0;i<36;i++){
  //   dydt_BAR[(id*36)+i] = dydt[i];
  // }
  dydt_BAR[(id*36)+  0] = dL;
  dydt_BAR[(id*36)+  1] = dB1AR;
  dydt_BAR[(id*36)+  2] = dGs;
  dydt_BAR[(id*36)+  3] = dB1AR_ACT;
  dydt_BAR[(id*36)+  4] = dB1AR_S464;
  dydt_BAR[(id*36)+  5] = dB1AR_S301;
  dydt_BAR[(id*36)+  6] = dGsaGTPtot;
  dydt_BAR[(id*36)+  7] = dGsaGDP;
  dydt_BAR[(id*36)+  8] = dGsBy;
  dydt_BAR[(id*36)+  9] = dGsaGTP;
  dydt_BAR[(id*36)+ 10] = dFsk;
  dydt_BAR[(id*36)+ 11] = dAC;
  dydt_BAR[(id*36)+ 12] = dPDE;
  dydt_BAR[(id*36)+ 13] = dIBMX;
  dydt_BAR[(id*36)+ 14] = dcAMPtot;
  dydt_BAR[(id*36)+ 15] = dcAMP;
  dydt_BAR[(id*36)+ 16] = dPKAC_I;
  dydt_BAR[(id*36)+ 17] = dPKAC_II;
  dydt_BAR[(id*36)+ 18] = dPLBp;
  dydt_BAR[(id*36)+ 19] = dInhib1ptot;
  dydt_BAR[(id*36)+ 20] = dInhib1p;
  dydt_BAR[(id*36)+ 21] = dPP1;
  dydt_BAR[(id*36)+ 22] = dLCCa_PKAp_whole;
  dydt_BAR[(id*36)+ 23] = dLCCb_PKAp_whole;
  dydt_BAR[(id*36)+ 24] = dRyR_PKAp_whole;
  dydt_BAR[(id*36)+ 25] = dTnI_PKAp_whole;
  dydt_BAR[(id*36)+ 26] = dIKs_PKAn;
  dydt_BAR[(id*36)+ 27] = dYotiao_KCQN1;
  dydt_BAR[(id*36)+ 28] = dIKs_PKAp_whole;
  dydt_BAR[(id*36)+ 29] = dICFTR_PKAp_whole;
  dydt_BAR[(id*36)+ 30] = dPLM_PKAp_whole;
  dydt_BAR[(id*36)+ 31] = dMyo_PKAp_whole;
  dydt_BAR[(id*36)+ 32] = dIKr_PKAn;
  dydt_BAR[(id*36)+ 33] = dYotiao_hERG;
  dydt_BAR[(id*36)+ 34] = dIKr_PKAp_whole;
  dydt_BAR[(id*36)+ 35] = dIClCa_PKAp_whole;
}
void CSubcell::solve_ODE_BAR(int id, double dtx){
  L[id]                   += dydt_BAR[ (id*36) +  0]*dtx;
  B1AR[id]                += dydt_BAR[ (id*36) +  1]*dtx;   
  Gs[id]                  += dydt_BAR[ (id*36) +  2]*dtx; 
  B1AR_ACT[id]            += dydt_BAR[ (id*36) +  3]*dtx;       
  B1AR_S464[id]           += dydt_BAR[ (id*36) +  4]*dtx;        
  B1AR_S301[id]           += dydt_BAR[ (id*36) +  5]*dtx;        
  GsaGTPtot[id]           += dydt_BAR[ (id*36) +  6]*dtx;        
  GsaGDP[id]              += dydt_BAR[ (id*36) +  7]*dtx;     
  GsBy[id]                += dydt_BAR[ (id*36) +  8]*dtx;   
  GsaGTP[id]              += dydt_BAR[ (id*36) +  9]*dtx;     
  Fsk[id]                 += dydt_BAR[ (id*36) + 10]*dtx; 
  AC[id]                  += dydt_BAR[ (id*36) + 11]*dtx;
  PDE[id]                 += dydt_BAR[ (id*36) + 12]*dtx; 
  IBMX[id]                += dydt_BAR[ (id*36) + 13]*dtx;  
  cAMPtot[id]             += dydt_BAR[ (id*36) + 14]*dtx;     
  cAMP[id]                += dydt_BAR[ (id*36) + 15]*dtx;  
  PKAC_I[id]              += dydt_BAR[ (id*36) + 16]*dtx;    
  PKAC_II[id]             += dydt_BAR[ (id*36) + 17]*dtx;     
  PLBp[id]                += dydt_BAR[ (id*36) + 18]*dtx;  
  Inhib1ptot[id]          += dydt_BAR[ (id*36) + 19]*dtx;        
  Inhib1p[id]             += dydt_BAR[ (id*36) + 20]*dtx;     
  PP1[id]                 += dydt_BAR[ (id*36) + 21]*dtx; 
  LCCa_PKAp_whole[id]     += dydt_BAR[ (id*36) + 22]*dtx;             
  LCCb_PKAp_whole[id]     += dydt_BAR[ (id*36) + 23]*dtx;             
  RyR_PKAp_whole[id]      += dydt_BAR[ (id*36) + 24]*dtx;            
  TnI_PKAp_whole[id]      += dydt_BAR[ (id*36) + 25]*dtx;            
  IKs_PKAn[id]            += dydt_BAR[ (id*36) + 26]*dtx;      
  Yotiao_KCQN1[id]        += dydt_BAR[ (id*36) + 27]*dtx;          
  IKs_PKAp_whole[id]      += dydt_BAR[ (id*36) + 28]*dtx;            
  ICFTR_PKAp_whole[id]    += dydt_BAR[ (id*36) + 29]*dtx;              
  PLM_PKAp_whole[id]      += dydt_BAR[ (id*36) + 30]*dtx;            
  Myo_PKAp_whole[id]      += dydt_BAR[ (id*36) + 31]*dtx;            
  IKr_PKAn[id]            += dydt_BAR[ (id*36) + 32]*dtx;      
  Yotiao_hERG[id]         += dydt_BAR[ (id*36) + 33]*dtx;         
  IKr_PKAp_whole[id]      += dydt_BAR[ (id*36) + 34]*dtx;            
  IClCa_PKAp_whole[id]    += dydt_BAR[ (id*36) + 35]*dtx;              
}

double CSubcell::computeave_RyRp_Multiplier(void)
{

  double sum=0;

  for (int id=0;id<n;id++){
    double fCKII_RyR = (20 * RyR_CKp[id] / 3.0 - 1 / 3.0); //
    double fPKA_RyR = (RyR_PKAp[id] * 1.025) + 0.9750; //
  
    sum+=(fCKII_RyR + fPKA_RyR - 1);
  }
  return (sum/nn);
}
#endif