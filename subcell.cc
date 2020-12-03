#include "subcell.h"

#ifdef ___DEBUG
#include <iomanip>
#include <limits>
#include <string>
#endif

//cell parameters
const double CSubcell::vjsr=0.02;
const double CSubcell::vp_ave=0.00126;

#include "diffusion.cc"


inline unsigned int xorshift(unsigned int *xx, unsigned int *yy, unsigned int *zz, unsigned int *ww)
{
  unsigned int t=(*xx^(*xx<<11));*xx=*yy;*yy=*zz;*zz=*ww;
  return( *ww=(*ww^(*ww>>19))^(t^(t>>8)) );
}


CSubcell::CSubcell(int sizex, int sizey, int sizez, int fmesh, double xii)
{
  if (fmesh == 1)
    dt=0.1;
  else if (fmesh == 5)
    dt=0.005;
  else if (fmesh == 3)
    dt=0.01;
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

#ifdef ___PTM
    if (cp[id]>1.0) {
      k12=k12*2;
      k43=k43*2;
    }

#endif


    double k14=MM/taub*BCSQN/BCSQN0;
    double k21=1/tauc1;
    double k23=MM/taub*BCSQN/BCSQN0;
    double k41=1/tauu;
    double k34=1/tauc2;
    double k32=k41*k12/k43;

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
      //          cout<<"RyR is negative "<<ryr1[id]<<"\t"<<ryr2[id]<<"\t"<<ryr3[id]<<endl;
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
  std::cout << "icaave" << std::end;
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

//Bardia 2020
void CSubcell::printsc0d(double t, double v,std::string filename){
  ofstream results0d;
  results0d.open(filename.c_str(), std::ofstream::app);

  results0d << t << "\t" 
            << v << "\t" 
            << this->computeaveci() << "\t" 
            << this->calcjcalave() << "\t" 
            << this->icaave << "\t"
            << this->calcjnacaave() << "\t"

            << std::endl;

  results0d.close();
}


// void CSubcell::voltage_clamp(double t, double peak_v, bool output, double dt=0.1, double min_v = -80.0, double APD = 250.0, double pcl=1000.0, double na_conc=12.0){
void CSubcell::voltage_clamp(double t, double peak_v, bool output, double dt, double min_v, double APD, double pcl, double na_conc){
 
  double relative_t = (int(t/dt)%int(pcl/dt))*dt;

  double v;

  if(relative_t < APD){
    v=peak_v;
  }
  else if(relative_t >= APD && relative_t <= pcl){
    v=min_v;
  }
  else{
    std::cout << "Error: unexpected time value: "<< relative_t << " (Voltage Clamp Error)";
    v=-80;
  }

  pace(v, na_conc);

  // std::cout << "Before output"  << std::endl;
  if(output){
    // char* fname_voltage_clamp;
    // sprintf(fname_voltage_clamp, "0d_results_vc_%02.0f.txt", peak_v);
    std::string fname_voltage_clamp = "0d_results_vc_20.txt";
    printsc0d(t, v, fname_voltage_clamp);
  }
  // std::cout << "AFter output" <<  std::endl;

}
