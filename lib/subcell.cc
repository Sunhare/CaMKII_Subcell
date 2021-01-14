#include "subcell.h"

#ifdef ___DEBUG
#include <iomanip>
#include <limits>
#include <string>
#endif

//cell parameters
const double CSubcell::vjsr=0.02;
const double CSubcell::vp_ave=0.00126;

// #include "diffusion.cc"


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

  //TODO Create arrays for Soltis-Saucerman/Negroni model

  #ifdef ___PTM
    allocate_memory_PTM_CaM_vars(n);
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

  #ifdef ___NO_CS_BUFFER
      csmn[id]=0;
  #else
      Ics[id]=0;
      Idps[id]=0;
  #endif

  #ifdef ___PTM
    init_parameters_PTM_CaM(id);
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
    for (int j=0;j<NoCaL;j++) //Mahajan model
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

    sumica+=Ica; //TODO Change ICa to include Soltis Saucerman PTMs

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
    double koSRCa = 10/kCaSR;
#else
    const double koSRCa=1;
#endif
    double sgmd=cp2/(Kcp*Kcp+cp2);
    double k12=koSRCa*Ku*sgmd+pedk12;
    double k43=koSRCa*Kb*sgmd+pedk43;
#endif

#ifdef ___PTM // TODO: Change PTMs

    // //D. Sato Quick Demo
    // if (cp[id]>1.0) {
    //   k12=k12*2;
    //   k43=k43*2;
    // }

    //Soltis Saucerman PTM 
    //0 < RyR_CKP < 1
    //0 < RYR_PKAp < 1
    double fCKII_RyR = (20 * RyR_CKp / 3.0 - 1 / 3.0); //≈ 1.0005, Max = 19/3, 6
    double fPKA_RyR = RyR_PKAp * 1.025 + 0.9750; //1.025*0.025 + 0.9750= 1.00065, Max = 2
    k12 *= (fCKII_RyR + fPKA_RyR - 1); //*= 1.0005, Max 7
    k43 *= (fCKII_RyR + fPKA_RyR - 1); // Max 7

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


#ifdef ___CONTRACTION //TODO implement Negroni 
#endif 

  //Main PTM ODEs
#ifdef ___PTM

    //CaM Dyad
    double* dydt_CaMDyad;
    dydt_CaMDyad = calc_dydt_CaM_ODEs(id, dydt_CaMDyad,
      //CaM Dyad/Cleft State Variables
      Ca2CaM_dyad[id],Ca4CaM_dyad[id],CaMB_dyad[id],Ca2CaMB_dyad[id],Ca4CaMB_dyad[id],
      Pb2_dyad[id],Pb_dyad[id],Pt_dyad[id],Pt2_dyad[id],Pa_dyad[id],
      CaMCaN_dyad[id],Ca2CaMCaN_dyad[id],Ca4CaMCaN_dyad[id], CaM_dyad[id],
      Ca4CaN_dyad[id], CaMCa4CaN_dyad[id],Ca2CaMCa4CaN_dyad[id],Ca4CaMCa4CaN_dyad[id],


      //CaM Dyad/Cleft Parameters
      CaMtotDyad[id],BtotDyad[id],CaMKIItotDyad[id],CaNtotDyad[id],PP1totDyad[id],cp[id], 2);

    //CaM SL
    double* dydt_CaMSL;
    dydt_CaMSL = calc_dydt_CaM_ODEs(id, dydt_CaMSL,
      //CaM Sarcolemmal State Variables
      Ca2CaM_sl[id],Ca4CaM_sl[id],CaMB_sl[id],Ca2CaMB_sl[id],Ca4CaMB_sl[id],
      Pb2_sl[id],Pb_sl[id],Pt_sl[id],Pt2_sl[id],Pa_sl[id],
      CaMCaN_sl[id],Ca2CaMCaN_sl[id],Ca4CaMCaN_sl[id], CaM_sl[id],
      Ca4CaN_sl[id], CaMCa4CaN_sl[id],Ca2CaMCa4CaN_sl[id],Ca4CaMCa4CaN_sl[id],


      //CaM Sarcolemmal Parameters
      CaMtotSL[id],BtotSL[id],CaMKIItotSL[id],CaNtotSL[id],PP1totSL[id],cs[id], 1);

    //CaM Cytosol
    double* dydt_CaMCyt;
    dydt_CaMCyt = calc_dydt_CaM_ODEs(id, dydt_CaMCyt,
      //CaM Cytosol State Variables
      Ca2CaM_cyt[id],Ca4CaM_cyt[id],CaMB_cyt[id],Ca2CaMB_cyt[id],Ca4CaMB_cyt[id],
      Pb2_cyt[id],Pb_cyt[id],Pt_cyt[id],Pt2_cyt[id],Pa_cyt[id],
      CaMCaN_cyt[id],Ca2CaMCaN_cyt[id],Ca4CaMCaN_cyt[id], CaM_cyt[id],
      Ca4CaN_cyt[id], CaMCa4CaN_cyt[id],Ca2CaMCa4CaN_cyt[id],Ca4CaMCa4CaN_cyt[id],


      //CaM Cytosol Parameters
      CaMtotCyt[id],BtotCyt[id],CaMKIItotCyt[id],CaNtotCyt[id],PP1totCyt[id],ci[id], 0);

    solve_ODE_CaM(id, dt, dydt_CaMDyad, 2);
    solve_ODE_CaM(id, dt, dydt_CaMSL, 1);
    solve_ODE_CaM(id, dt, dydt_CaMCyt, 0);


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

//Signaling network (ODEs go here, before dci) TODO
    // dRYR_CKp = ci[id]*
    // RyR_CKp[id] += dRYR_CKp*dt;

    ci[id]+=dci*dt; // Update the  
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

void CSubcell::setboundary(int bcc){
  if (bc>0){//corner

    Jmaxx[0+0*nx+0*nxny]=0;
    Jmaxx[0+(ny-1)*nx+(nz-1)*nxny]=0;
    Jmaxx[0+(ny-1)*nx+0*nxny]=0;
    Jmaxx[0+0*nx+(nz-1)*nxny]=0;
    Jmaxx[(nx-1)+(ny-1)*nx+0*nxny]=0;
    Jmaxx[(nx-1)+0*nx+(nz-1)*nxny]=0;
    Jmaxx[(nx-1)+0*nx+0*nxny]=0;
    Jmaxx[(nx-1)+(ny-1)*nx+(nz-1)*nxny]=0;
  }
  if (bc>1){//edge
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

//PTMs 2020-2021
void CSubcell::printsc0d(double t, double v,std::string filename){
  ofstream results0d;
  results0d.open(filename.c_str(), std::ofstream::app);

  results0d << t << "\t" 
            << v << "\t" 
            << this->computeaveci() << "\t" 
            << this->calcjcalave() << "\t" 
            << this->icaave << "\t"
            << this->calcjnacaave() << "\t"

            //SR Ca concentration
            //SR Release
            //

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
    std::string fname_voltage_clamp = "results/0d_results_vc_20.txt";
    printsc0d(t, v, fname_voltage_clamp);
  }

}


void CSubcell::allocate_memory_PTM_CaM_vars(int n){

  RyR_CKp = new double(n);
  RyR_PKAp = new double(n);
  //All Parameters
  CaMtotDyad = new double(n);
  BtotDyad = new double(n);
  CaMKIItotDyad = new double(n);
  CaNtotDyad = new double(n);
  PP1totDyad = new double(n);
  CaMtotSL = new double(n);
  BtotSL = new double(n);
  CaMKIItotSL = new double(n);
  CaNtotSL = new double(n);
  PP1totSL = new double(n);
  CaMtotCyt = new double(n);
  BtotCyt = new double(n);
  CaMKIItotCyt = new double(n);
  CaNtotCyt = new double(n);
  PP1totCyt = new double(n);
  LCCtotDyad = new double(n);
  RyRtot = new double(n);
  PP1_dyad = new double(n);
  PP2A_dyad = new double(n);
  OA = new double(n);
  LCCtotSL = new double(n);
  PP1_SL = new double(n);
  PLBtot = new double(n);
  Ligtot = new double(n);
  LCCtotBA = new double(n);
  RyRtotBA = new double(n);
  PLBtotBA = new double(n);
  TnItotBA = new double(n);
  IKstotBA = new double(n);
  ICFTRtotBA = new double(n);
  PP1_PLBtot = new double(n);
  PLMtotBA = new double(n);
  MyototBA = new double(n);
  IKrtotBA = new double(n);
  IClCatotBA = new double(n);
  CKIIOE = new double(n);
  recoveryTime = new double(n);

  ////CaM State Variables

  //Ca Fluxes from CaM
  JCaCyt = new double(n);
  JCaSL = new double(n);
  JCaDyad = new double(n);

  //CaM_Dyad
  CaM_dyad = new double(n);
  Ca2CaM_dyad = new double(n);
  Ca4CaM_dyad = new double(n);
  CaMB_dyad = new double(n);
  Ca2CaMB_dyad = new double(n);
  Ca4CaMB_dyad = new double(n);
  Pb2_dyad = new double(n);
  Pb_dyad = new double(n);
  Pt_dyad = new double(n);
  Pt2_dyad = new double(n);
  Pa_dyad = new double(n);
  Ca4CaN_dyad = new double(n);
  CaMCa4CaN_dyad = new double(n);
  Ca2CaMCa4CaN_dyad = new double(n);
  Ca4CaMCa4CaN_dyad = new double(n);

  //CaM_SL
  CaM_sl = new double(n);
  Ca2CaM_sl = new double(n);
  Ca4CaM_sl = new double(n);
  CaMB_sl = new double(n);
  Ca2CaMB_sl = new double(n);
  Ca4CaMB_sl = new double(n);
  Pb2_sl = new double(n);
  Pb_sl = new double(n);
  Pt_sl = new double(n);
  Pt2_sl = new double(n);
  Pa_sl = new double(n);
  Ca4CaN_sl = new double(n);
  CaMCa4CaN_sl = new double(n);
  Ca2CaMCa4CaN_sl = new double(n);
  Ca4CaMCa4CaN_sl = new double(n);

  //CaM_cytosol
  CaM_cyt = new double(n);
  Ca2CaM_cyt = new double(n);
  Ca4CaM_cyt = new double(n);
  CaMB_cyt = new double(n);
  Ca2CaMB_cyt = new double(n);
  Ca4CaMB_cyt = new double(n);
  Pb2_cyt = new double(n);
  Pb_cyt = new double(n);
  Pt_cyt = new double(n);
  Pt2_cyt = new double(n);
  Pa_cyt = new double(n);
  Ca4CaN_cyt = new double(n);
  CaMCa4CaN_cyt = new double(n);
  Ca2CaMCa4CaN_cyt = new double(n);
  Ca4CaMCa4CaN_cyt = new double(n);
}


void CSubcell::init_parameters_PTM_CaM(int n){
  // freq = 1.0;                 // [Hz] CHANGE DEPENDING ON FREQUENCY
  // cycleLength = 1e3/freq;     // [ms]
  CaMtotDyad[n] = 418;             //[uM]
  BtotDyad[n] = 1.54/8.293e-4;        //[uM]
  CaMKIItotDyad[n] = 120;        //[uM]
  CaNtotDyad[n] = 3e-3/8.293e-4;    //[uM]
  PP1totDyad[n] = 96.5;             //[uM]
  CaMtotSL[n] = 5.65;                 //[uM]
  BtotSL[n] = 24.2;                     //[uM]
  CaMKIItotSL[n] = 120*8.293e-4;   //[uM]
  CaNtotSL[n] = 3e-3;                 //[uM]
  PP1totSL[n] = 0.57;                 //[uM]
  CaMtotCyt[n] = 5.65;               //[uM]
  BtotCyt[n] = 24.2;                   //[uM]
  CaMKIItotCyt[n] = 120*8.293e-4; //[uM]
  CaNtotCyt[n] = 3e-3;               //[uM]
  PP1totCyt[n] = 0.57;               //[uM]

  // ADJUST CAMKII ACTIVITY LEVELS (expression = 'WT', 'OE', or 'KO')
  expression[n] = (char* )"WT";
  CKIIOE[n] = 0; // Should be zero during "WT" and "KO" runs

  if(expression[n] =="OE"){
    CKIIOE[n] = 1; // Flag for CKII OE in ECC file (0=WT, 1=OE) - for Ito and INa
    CaMKIItotDyad[n] = 120*6;          // [uM] 
    CaMKIItotSL[n] = 120*8.293e-4*6;   // [uM]
    CaMKIItotCyt[n] = 120*8.293e-4*6;  // [uM]
  }
  else if(expression[n] == "KO"){
    CaMKIItotDyad[n] = 0;              // [uM] 
    CaMKIItotSL[n] = 0;                // [uM]
    CaMKIItotCyt[n] = 0;               // [uM]
  }
  else{
    expression[n] = (char* )"WT";
  }

  // Parameters for CaMKII module
  LCCtotDyad[n] = 31.4*.9;       // [uM] - Total Dyadic [LCC] - (umol/l dyad)
  LCCtotSL[n] = 0.0846;          // [uM] - Total Subsarcolemmal [LCC] (umol/l sl)
  RyRtot[n] = 382.6;             // [uM] - Total RyR (in Dyad)
  PP1_dyad[n] = 95.7;            // [uM] - Total dyadic [PP1]
  PP1_SL[n] = 0.57;              // [uM] - Total Subsarcolemmal [PP1]
  PP2A_dyad[n] = 95.76;          // [uM] - Total dyadic PP2A
  OA[n] = 0;                     // [uM] - PP1/PP2A inhibitor Okadaic Acid
  PLBtot[n] = 38;                // [uM] - Total [PLB] in cytosolic units

  // Parameters for BAR module
  Ligtot[n] = 0;// 0.1 or 0.02    // [uM] - SET LIGAND CONCENTRATION HERE
  LCCtotBA[n] = 0.025;           // [uM] - [umol/L cytosol]
  RyRtotBA[n] = 0.135;           // [uM] - [umol/L cytosol]
  PLBtotBA[n] = PLBtot[n];          // [uM] - [umol/L cytosol]
  TnItotBA[n] = 70;              // [uM] - [umol/L cytosol]
  IKstotBA[n] = 0.025;           // [uM] - [umol/L cytosol]
  ICFTRtotBA[n] = 0.025;         // [uM] - [umol/L cytosol]
  PP1_PLBtot[n] = 0.89;          // [uM] - [umol/L cytosol]
  PLMtotBA[n] = 48;              // [uM] - [umol/L cytosol]
  MyototBA[n] = 70;              // [uM] - [umol/L cytosol]
  IKrtotBA[n] = 0.025;           // [uM] - [umol/L cytosol]
  IClCatotBA[n] = 0.025;         // [uM] - [umol/L cytosol]

  // For Recovery from inactivation of LCC
  recoveryTime[n] = 10;  // initialize to smallest value

}


  // pCaMDyad = [K, Mg, CaMtotDyad, 0, CaMKIItotDyad, CaNtotDyad, PP1totDyad, CaDyad, cycleLength, compart_dyad];
void CSubcell::calc_dydt_CaM_ODEs(
  int id, double* dydt_CaM,
  //CaM State Variables
  double& Ca2CaM, double& Ca4CaM, double& CaMB, double& Ca2CaMB, double& Ca4CaMB,
  double& Pb2, double& Pb, double& Pt, double& Pt2, double& Pa,
  double& CaMCaN, double& Ca2CaMCaN, double& Ca4CaMCaN, double& CaM,
  double& Ca4CaN, double& CaMCa4CaN, double& Ca2CaMCa4CaN, double& Ca4CaMCa4CaN,

  //CaM Parameters 
  double CaMtot, double Btot, double CaMKIItot, double CaNtot, double PP1tot, double Ca, 
  int CaM_compartment)

{
  double K = 135; // mM
  double Mg = 1; //mM


  // Negroni et al model - CaM module

  //// Descriptions for state variables
  // Ca2CaM = y(1);  // 2 Ca bound to C terminal sites
  // Ca4CaM = y(2);  // 4 Ca bound
  // CaMB = y(3);   
  // Ca2CaMB = y(4);
  // Ca4CaMB = y(5);
  // Pb2 = y(6);     // probability of a Ca2CaM bound CaMKII subunit
  // Pb = y(7);      // probability of a Ca4CaM bound CaMKII subunit
  // Pt = y(8);      // probability of a Ca4CaM bound autophosphorylated CaMKII subunit
  // Pt2 = y(9);     // probability of a Ca2CaMxx bound autophosphorylated CaMKII subunit
  // Pa = y(10);     // probability of an autonomous autophosphorylated CaMKII subunit
  // CaMCaN = y(11);
  // Ca2CaMCaN = y(12);
  // Ca4CaMCaN = y(13);  // active calcineurin
  // description of intermediate variables
  // CaM- Ca free CaM
  ////

  // Parameters
  // Ca/CaM parameters
  double Kd02, Kd24;
  if(Mg <= 1){
      Kd02 = 0.0025*(1+K/0.94-Mg/0.012)*(1+K/8.1+Mg/0.022);  // [uM^2]
      Kd24 = 0.128*(1+K/0.64+Mg/0.0014)*(1+K/13.0-Mg/0.153); // [uM^2]
  }
  else{
      Kd02 = 0.0025*(1+K/0.94-1/0.012+(Mg-1)/0.060)*(1+K/8.1+1/0.022+(Mg-1)/0.068);   // [uM^2]
      Kd24 = 0.128*(1+K/0.64+1/0.0014+(Mg-1)/0.005)*(1+K/13.0-1/0.153+(Mg-1)/0.150);  // [uM^2]
  }

  double k20 = 10;               // [s^-1]      
  double k02 = k20/Kd02;         // [uM^-2 s^-1]
  double k42 = 500;              // [s^-1]      
  double k24 = k42/Kd24;         // [uM^-2 s^-1]

  // CaM buffering (B) parameters
  double k0Boff = 0.0014;        // [s^-1] 
  double k0Bon = k0Boff/0.2;   // [uM^-1 s^-1] kon = koff/Kd
  double k2Boff = k0Boff/100;    // [s^-1] 
  double k2Bon = k0Bon;          // [uM^-1 s^-1]
  double k4Boff = k2Boff;        // [s^-1]
  double k4Bon = k0Bon;          // [uM^-1 s^-1]
  // using thermodynamic constraints
  double k20B = k20/100; // [s^-1] thermo constraint on loop 1
  double k02B = k02;     // [uM^-2 s^-1] 
  double k42B = k42;     // [s^-1] thermo constraint on loop 2
  double k24B = k24;     // [uM^-2 s^-1]

  // CaMKII parameters
  // Wi Wa Wt Wp
  double kbi = 2.2;      // [s^-1] (Ca4CaM dissocation from Wb)
  double kib = kbi/33.5e-3; // [uM^-1 s^-1]
  double kib2 = kib;
  double kb2i = kib2*5;
  double kb24 = k24;
  double kb42 = k42*33.5e-3/5;
  double kpp1 = 1.72;    // [s^-1] (PP1-dep dephosphorylation rates)
  double Kmpp1 = 11.5;   // [uM]
  double kta = kbi/1000; // [s^-1] (Ca4CaM dissociation from Wt)
  double kat = kib;      // [uM^-1 s^-1] (Ca4CaM reassociation with Wa)
  double kt42 = k42*33.5e-6/5;
  double kt24 = k24;
  double kat2 = kib;
  double kt2a = kib*5;

  // CaN parameters
  double kcanCaoff = 1;              // [s^-1] 
  double kcanCaon = kcanCaoff/0.5;   // [uM^-1 s^-1] 
  double kcanCaM4on = 46;            // [uM^-1 s^-1]
  double kcanCaM4off = 1.3e-3;       // [s^-1]
  double kcanCaM2on = kcanCaM4on;
  double kcanCaM2off = 2508*kcanCaM4off;
  double kcanCaM0on = kcanCaM4on;
  double kcanCaM0off = 165*kcanCaM2off;
  double k02can = k02;
  double k20can = k20/165;
  double k24can = k24;
  double k42can = k20/2508;

  // CaM Reaction fluxes
  double rcn02 = k02*pow(Ca,2)*CaM - k20*Ca2CaM;
  double rcn24 = k24*pow(Ca,2)*Ca2CaM - k42*Ca4CaM;
  // CaM buffer fluxes
  double B = Btot - CaMB - Ca2CaMB - Ca4CaMB;
  double rcn02B = k02B*pow(Ca,2)*CaMB - k20B*Ca2CaMB;
  double rcn24B = k24B*pow(Ca,2)*Ca2CaMB - k42B*Ca4CaMB;
  double rcn0B = k0Bon*CaM*B - k0Boff*CaMB;
  double rcn2B = k2Bon*Ca2CaM*B - k2Boff*Ca2CaMB;
  double rcn4B = k4Bon*Ca4CaM*B - k4Boff*Ca4CaMB;
  // CaN reaction fluxes 
  double Ca2CaN = CaNtot - Ca4CaN - CaMCa4CaN - Ca2CaMCa4CaN - Ca4CaMCa4CaN;
  double rcnCa4CaN = kcanCaon*pow(Ca,2)*Ca2CaN - kcanCaoff*Ca4CaN;
  double rcn02CaN = k02can*pow(Ca,2)*CaMCa4CaN - k20can*Ca2CaMCa4CaN; 
  double rcn24CaN = k24can*pow(Ca,2)*Ca2CaMCa4CaN - k42can*Ca4CaMCa4CaN;
  double rcn0CaN = kcanCaM0on*CaM*Ca4CaN - kcanCaM0off*CaMCa4CaN;
  double rcn2CaN = kcanCaM2on*Ca2CaM*Ca4CaN - kcanCaM2off*Ca2CaMCa4CaN;
  double rcn4CaN = kcanCaM4on*Ca4CaM*Ca4CaN - kcanCaM4off*Ca4CaMCa4CaN;
  // CaMKII reaction fluxes
  double Pi = 1 - Pb2 - Pb - Pt - Pt2 - Pa;
  double rcnCKib2 = kib2*Ca2CaM*Pi - kb2i*Pb2;
  double rcnCKb2b = kb24*pow(Ca,2)*Pb2 - kb42*Pb;
  double rcnCKib = kib*Ca4CaM*Pi - kbi*Pb;
  double T = Pb + Pt + Pt2 + Pa;
  double kbt = 0.055*T + .0074*pow(T,2) + 0.015*pow(T,3);
  double rcnCKbt = kbt*Pb - kpp1*PP1tot*Pt/(Kmpp1+CaMKIItot*Pt);
  double rcnCKtt2 = kt42*Pt - kt24*pow(Ca,2)*Pt2;
  double rcnCKta = kta*Pt - kat*Ca4CaM*Pa;
  double rcnCKt2a = kt2a*Pt2 - kat2*Ca2CaM*Pa;
  double rcnCKt2b2 = kpp1*PP1tot*Pt2/(Kmpp1+CaMKIItot*Pt2);
  double rcnCKai = kpp1*PP1tot*Pa/(Kmpp1+CaMKIItot*Pa);

  // CaM equations
  double dCaM = 1e-3*(-rcn02 - rcn0B - rcn0CaN);
  double dCa2CaM = 1e-3*(rcn02 - rcn24 - rcn2B - rcn2CaN + CaMKIItot*(-rcnCKib2 + rcnCKt2a) );
  double dCa4CaM = 1e-3*(rcn24 - rcn4B - rcn4CaN + CaMKIItot*(-rcnCKib+rcnCKta) );
  double dCaMB = 1e-3*(rcn0B-rcn02B);
  double dCa2CaMB = 1e-3*(rcn02B + rcn2B - rcn24B);
  double dCa4CaMB = 1e-3*(rcn24B + rcn4B);

  // CaMKII equations
  double dPb2 = 1e-3*(rcnCKib2 - rcnCKb2b + rcnCKt2b2); // Pb2
  double dPb = 1e-3*(rcnCKib + rcnCKb2b - rcnCKbt);    // Pb
  double dPt = 1e-3*(rcnCKbt-rcnCKta-rcnCKtt2);        // Pt
  double dPt2 = 1e-3*(rcnCKtt2-rcnCKt2a-rcnCKt2b2);     // Pt2
  double dPa = 1e-3*(rcnCKta+rcnCKt2a-rcnCKai);       // Pa

  // CaN equations
  double dCa4CaN = 1e-3*(rcnCa4CaN - rcn0CaN - rcn2CaN - rcn4CaN);                       // Ca4CaN
  double dCaMCa4CaN = 1e-3*(rcn0CaN - rcn02CaN);           // CaMCa4CaN
  double dCa2CaMCa4CaN = 1e-3*(rcn2CaN+rcn02CaN-rcn24CaN);    // Ca2CaMCa4CaN
  double dCa4CaMCa4CaN = 1e-3*(rcn4CaN+rcn24CaN);             // Ca4CaMCa4CaN

  double dydt[]={dCaM,dCa2CaM,dCa4CaM,dCaMB,dCa2CaMB,dCa4CaMB,dPb2,dPb,dPt,dPt2,dPa,dCa4CaN,dCaMCa4CaN,dCa2CaMCa4CaN,dCa4CaMCa4CaN};
  
  //Assign dydt
  dydt_CaM = dydt;

  // write to global variables for adjusting Ca buffering in EC coupling model
  double JCa = 1e-3*(2*CaMKIItot*(rcnCKtt2-rcnCKb2b) - 2*(rcn02+rcn24+rcn02B+rcn24B+rcnCa4CaN+rcn02CaN+rcn24CaN)); // [uM/msec]

  if(CaM_compartment == 0){
    JCaCyt[id] = JCa;
  }
  else if(CaM_compartment == 1){
    JCaSL[id] = JCa;
  }
  else if(CaM_compartment == 2){
    JCaDyad[id] = JCa;
  }
  else{
    std::cout <<"Error: Unrecognized CaM_compartment in subcell.cc: " << CaM_compartment << std::endl;
    std::exit(1);
  }

}


void CSubcell::solve_ODE_CaM(int id, double dt, double* dydt_CaM, int CaM_compartment){
  int n_CaM_params = 15;
  // double dydt[]={dCaM,dCa2CaM,dCa4CaM,dCaMB,dCa2CaMB,dCa4CaMB,dPb2,dPb,dPt,dPt2,dPa,dCa4CaN,dCaMCa4CaN,dCa2CaMCa4CaN,dCa4CaMCa4CaN};
  //CaM Cytosol
  if(CaM_compartment == 0){
    CaM_cyt[id] += dydt_CaM[0]*dt;
    Ca2CaM_cyt[id] += dydt_CaM[1]*dt;
    Ca4CaM_cyt[id] += dydt_CaM[2]*dt;
    CaMB_cyt[id] += dydt_CaM[3]*dt;
    Ca2CaMB_cyt[id] += dydt_CaM[4]*dt;
    Ca4CaMB_cyt[id] += dydt_CaM[5]*dt;
    Pb2_cyt[id] += dydt_CaM[6]*dt;
    Pb_cyt[id] += dydt_CaM[7]*dt;
    Pt_cyt[id] += dydt_CaM[8]*dt;
    Pt2_cyt[id] += dydt_CaM[9]*dt;
    Pa_cyt[id] += dydt_CaM[10]*dt;
    Ca4CaN_cyt[id] += dydt_CaM[11]*dt;
    CaMCa4CaN_cyt[id] += dydt_CaM[12]*dt;
    Ca2CaMCa4CaN_cyt[id] += dydt_CaM[13]*dt;
    Ca4CaMCa4CaN_cyt[id] += dydt_CaM[14]*dt;
  } 
  //CaM Sarcolemmal
  else if(CaM_compartment == 1){
    CaM_sl[id] += dydt_CaM[0]*dt;
    Ca2CaM_sl[id] += dydt_CaM[1]*dt;
    Ca4CaM_sl[id] += dydt_CaM[2]*dt;
    CaMB_sl[id] += dydt_CaM[3]*dt;
    Ca2CaMB_sl[id] += dydt_CaM[4]*dt;
    Ca4CaMB_sl[id] += dydt_CaM[5]*dt;
    Pb2_sl[id] += dydt_CaM[6]*dt;
    Pb_sl[id] += dydt_CaM[7]*dt;
    Pt_sl[id] += dydt_CaM[8]*dt;
    Pt2_sl[id] += dydt_CaM[9]*dt;
    Pa_sl[id] += dydt_CaM[10]*dt;
    Ca4CaN_sl[id] += dydt_CaM[11]*dt;
    CaMCa4CaN_sl[id] += dydt_CaM[12]*dt;
    Ca2CaMCa4CaN_sl[id] += dydt_CaM[13]*dt;
    Ca4CaMCa4CaN_sl[id] += dydt_CaM[14]*dt;
  }
  //CaM Dyad/Cleft space
  else if(CaM_compartment == 2){
    CaM_dyad[id] += dydt_CaM[0]*dt;
    Ca2CaM_dyad[id] += dydt_CaM[1]*dt;
    Ca4CaM_dyad[id] += dydt_CaM[2]*dt;
    CaMB_dyad[id] += dydt_CaM[3]*dt;
    Ca2CaMB_dyad[id] += dydt_CaM[4]*dt;
    Ca4CaMB_dyad[id] += dydt_CaM[5]*dt;
    Pb2_dyad[id] += dydt_CaM[6]*dt;
    Pb_dyad[id] += dydt_CaM[7]*dt;
    Pt_dyad[id] += dydt_CaM[8]*dt;
    Pt2_dyad[id] += dydt_CaM[9]*dt;
    Pa_dyad[id] += dydt_CaM[10]*dt;
    Ca4CaN_dyad[id] += dydt_CaM[11]*dt;
    CaMCa4CaN_dyad[id] += dydt_CaM[12]*dt;
    Ca2CaMCa4CaN_dyad[id] += dydt_CaM[13]*dt;
    Ca4CaMCa4CaN_dyad[id] += dydt_CaM[14]*dt;
  }
  else{
    std::cout << "Error: Unreognized CaM_compartment in solve_ODE_CaM: " << CaM_compartment << std::endl;
    exit(0);
  }
}























