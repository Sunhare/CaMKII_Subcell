#include "ap.h"


const double CActionPotential::vc=-80;
const double CActionPotential::stim=40;
const double CActionPotential::stimduration=1;


// ---------------constant parameters ------------------------------
const double CActionPotential::temp=308.0;// temperature (K)
const double CActionPotential::xxr=8.314;//
const double CActionPotential::xf=96.485;// Faraday's constant
const double CActionPotential::frt=xf/(xxr*temp);


CActionPotential::CActionPotential(void) : y(new double[N]),
  xm(y[0]), xh(y[1]), xj(y[2]), xr(y[3]), 
  xs1(y[4]), xs2(y[5]), xtos(y[6]), ytos(y[7]),
  v(y[8]), nai(y[9]), xtof(y[10]), ytof(y[11]), 
  rtos(y[12])
{
  // initial conditions
  xm=0.001;// sodium m-gate
  xh=0.99;// sodium h-gate
  xj=0.99;// soiumj-gate

  xr=0.01;// ikr gate variable 
  xs1=0.1;// iks gate variable
  xs2=0.14;// iks gate varaible 


  v=-86.8; // voltage

  // Markov gate variables 

  nai=14.0;// internal Na conc.

  xtof=0.003;// ito fast activation
  ytof=0.98;// ito fast inactivation
  xtos=0.003;// ito slow activation
  ytos=0.155;// ito slow inactivation
  rtos=0.98;// ito slow inactivation

  hpde=0.1;
  vold = v;
  jparam=1;

  _inaca=_ica=_iks=_ikr=_ik1=_ina=_inak=_itof=_itos=0;
  _inabk=_icl=_iclbk=_icabk=_islcap=0;

  xnao=136.0;//mM   external Na
  xki=140.0;// mM  internal K
  xko=5.40;//mM  external K
  cao=1.8;// mM  external Ca

  gtos=0.04;// ito slow conductance 
  gtof=0.11;// ito fast conductance 
  gkr=0.0125;// Ikr conductance 
  gks=0.32;
  gk1=0.3;// Ik1 conductance
  gnak=1.5;
  gna=12.0;// sodium conductance (mS/micro F) 

  gnabk=0.297*0.001;
  gcl=0.109625*0.1;
  gclbk=0.009*0.1;

}
CActionPotential::~CActionPotential()
{
  delete[] y;
}

CActionPotential& CActionPotential::operator=(const CActionPotential& ap)
{
  if (&ap==this)return(*this);

  //constructor
  xm=ap.xm;// sodium m-gate
  xh=ap.xh;// sodium h-gate
  xj=ap.xj;// soiumj-gate

  xr=ap.xr;// ikr gate variable 
  xs1=ap.xs1;// iks gate variable
  xs2=ap.xs2;// iks gate varaible 


  v=ap.v; // voltage

  // Markov gate variables 

  nai=ap.nai;// internal Na conc.

  xtof=ap.xtof;// ito fast activation
  ytof=ap.ytof;// ito fast inactivation
  xtos=ap.xtos;// ito slow activation
  ytos=ap.ytos;// ito slow inactivation
  rtos=ap.rtos;// ito slow inactivation

  hpde=ap.hpde;
  vold = ap.vold;
  jparam=ap.jparam;

  _inaca=_ica=_iks=_ikr=_ik1=_ina=_inak=_itof=_itos=0;
  _inabk=_icl=_iclbk=_icabk=_islcap=0;

  xnao=ap.xnao;//mM   external Na
  xki=ap.xki;// mM   internal K
  xko=ap.xko;//mM  external K
  cao=ap.cao;// mM   external Ca

  gtos=ap.gtos;// ito slow conductance 
  gtof=ap.gtof;// ito fast conductance 
  gkr=ap.gkr;// Ikr conductance 
  gks=ap.gks;
  gk1=ap.gk1;// Ik1 conductance
  gnak=ap.gnak;
  gna=ap.gna;// sodium conductance (mS/micro F) 


  gnabk=ap.gnabk;
  gcl=ap.gcl;
  gclbk=ap.gclbk;

  return(*this);
}


void CActionPotential::apclamp(double t, double pcl, double vmin, double vmax, double apd, bool pulse)
{
  if (apd==0)
  {
    const double a=2.0/3.0*1000;
    double x=a/(a+pcl);
    int m=(int)(t/pcl);
    if (m*pcl+x*pcl>t)
    {
      if (pulse)
        v=vmax;
      else
        v=vmin+(vmax-vmin)*sqrt(1-((t-m*pcl)/x/pcl)*((t-m*pcl)/x/pcl));
    }
    else
    {
      v=vmin;
    }
    nai=78.0/(1+10.0*sqrt(pcl/1000));
  }
  else
  {
    double x=apd/pcl;
    int m=(int)(t/pcl);
    if (m*pcl+x*pcl>t)
    {
      if (pulse)
        v=vmax;
      else
        v=vmin+(vmax-vmin)*sqrt(1-((t-m*pcl)/x/pcl)*((t-m*pcl)/x/pcl));
    }
    else
    {
      v=vmin;
    }
  }
}


void CActionPotential::pace(double st,double jcal, double jnaca, double icabk, double islcap, double cai)
{
  double dv=(vold-v)/hpde;
  vold=v;
  double Itotal;
  if(fabs(dv)>25.0)// then finer time step when dv/dt large
  {
    hode=hpde/10;
    for (int iii=0;iii<10;iii++)
    {
      pacex(st, jcal, jnaca, icabk, islcap, cai);
    }
  }
  else
  {
    hode=hpde;
    pacex(st, jcal, jnaca, icabk, islcap, cai);
  }
}

void CActionPotential::pacex(double st,double jcal, double jnaca, double icabk, double islcap, double cai)
{
  double ik1=comp_ik1();
  double ito=comp_ito();//itos and itof
  double inak=comp_inak();

  double ina=comp_ina();
  double ikr=comp_ikr();
  double iks=comp_iks(cai);

  //additional currents
  double inabk=comp_inabk();
  double icl=comp_icl(cai);
  double iclbk=comp_iclbk();
  //  double icabk=comp_icabk(cai);
  //  double islcap=comp_islcap(cai);


  //-------convert ion flow to current---------
  const double wca=8.0;//conversion factor between micro molar/ms to micro amps/ micro farads
  //  double xinaca=wca*jnaca;
  //  double xica=-2.0*wca*jcal;
  //--------sodium dynamics -------------------------
  const double xrr=(1.0/wca)/1000.0;// note: sodium is in m molar so need to divide by 1000


  // conversion Jca (sum(Jca*vp)/vp_ave/n)  Jncx (sum(Jncx)/n
  //  const double wca1=vs*1000*10^-9*1*96485*(65*27*11)/310; =0.1502 see below
  const double wca1=0.1502;
  //  const double wca2=vp_ave*1000*10^-9*2*96485*(65*27*11)/310; =0.0151 see below
  const double wca2=0.0151;
  double inaca=wca1*jnaca;
  double ica=wca2*jcal;


  nai+=(-xrr*(ina+3.0*inak+3.0*inaca+inabk))*hode;
  //  nai+=(-xrr*(ina+3.0*inak+3.0*inaca))*hode;
  // -------- dV/dt ------------------------------------
  double Itotal=(-(ina+ik1+ikr+iks+ito+inaca+ica+inak+(inabk+icl+iclbk+icabk+islcap))+ st);
  //  double Itotal=(-(ina+xik1+ikr+iks+xito+inaca+ica+xinak)+ st);
  v+=Itotal*hode;

  _inaca=inaca;_ica=ica;_iks=iks;_ikr=ikr;_ik1=ik1;_ina=ina;_inak=inak;
  _inabk=inabk;_icl=icl;_iclbk=iclbk;_icabk=icabk;_islcap=islcap;

}
//----------- sodium current following Hund-Rudy -------------------
double CActionPotential::comp_ina(void)
{
  double ena = (1.0/frt)*log(xnao/nai);
  double am=3.2;
  if (fabs(v+47.13)>0.001/0.1)
    am = 0.32*(v+47.13)/(1.0-exp(-0.1*(v+47.13)));
  double bm = 0.08*exp(-v/11.0);
  double ah,bh,aj,bj;
  if(v<(-40.0))
  {
    ah=0.135*exp((80.0+v)/(-6.8));
    bh=3.56*exp(0.079*v)+310000.0*exp(0.35*v);
    aj=((-127140.0*exp(0.2444*v)-0.00003474*exp(-0.04391*v))*(v+37.78))/(1.0+exp(0.311*(v+79.23)));
    bj=(0.1212*exp(-0.01052*v))/(1.0+exp(-0.1378*(v+40.14)));
  }
  else
  {
    ah=0.0;
    bh=1.0/(0.13*(1.0+exp((v+10.66)/(-11.1))));
    aj=0.0;
    bj=(0.3*exp(-0.0000002535*v))/(1.0+exp(-0.1*(v+32.0)));
  }

  double tauh=1.0/(ah+bh);
  double tauj=1.0/(aj+bj)*jparam;
  double taum=1.0/(am+bm);
  double xina= gna*xh*xj*xm*xm*xm*(v-ena);

  xh = ah/(ah+bh)-((ah/(ah+bh))-xh)*exp(-hode/tauh);
  xj = aj/(aj+bj)-((aj/(aj+bj))-xj)*exp(-hode/tauj);
  xm = am/(am+bm)-((am/(am+bm))-xm)*exp(-hode/taum);
  return xina;
}
//-------------- Ikr following Shannon------------------ 
double CActionPotential::comp_ikr(void)
{
  double ek = (1.0/frt)*log(xko/xki);// K reversal potential
  const double gss=sqrt(xko/5.4);
  double xkrv1;
  if (fabs(v+7.0)<0.001/0.123)
    xkrv1=0.00138/0.123;
  else
    xkrv1=0.00138*(v+7.0)/( 1.-exp(-0.123*(v+7.0)));
  double xkrv2;
  if (fabs(v+10.0)<0.001/0.145)
    xkrv2=0.00061/0.145;
  else
    xkrv2=0.00061*(v+10.0)/(exp( 0.145*(v+10.0))-1.0);
  double taukr=1.0/(xkrv1+xkrv2);
  double xkrinf=1.0/(1.0+exp(-(v+50.0)/7.5));
  double rg=1.0/(1.0+exp((v+33.0)/22.4));
  double xikr=gkr*gss*xr*rg*(v-ek);
  xr=xkrinf-(xkrinf-xr)*exp(-hode/taukr);
  return xikr;
}
// ----- Iks modified from Shannon, with new Ca dependence------------
double CActionPotential::comp_iks(double ci)
{
  const double prnak=0.018330;
  double eks=(1.0/frt)*log((xko+prnak*xnao)/(xki+prnak*nai));
  double xs1ss=1.0/(1.0+exp(-(v-1.50)/16.70));
  double xs2ss=xs1ss;
  double tauxs1;
  if (fabs(v+30.0)<0.001/0.0687)
    tauxs1=1/(0.0000719/0.148+0.000131/0.0687);
  else
    tauxs1=1.0/(0.0000719*(v+30.0)/(1.0-exp(-0.148*(v+30.0)))+0.000131*(v+30.0)/(exp(0.0687*(v+30.0))-1.0));
  double tauxs2=4*tauxs1;
  double gksx=0.433*(1+0.8/(1+pow((0.5/ci),3)));
  double xiks=gks*gksx*xs1*xs2*(v-eks);
  xs1=xs1ss-(xs1ss-xs1)*exp(double(-hode/tauxs1));
  xs2=xs2ss-(xs2ss-xs2)*exp(double(-hode/tauxs2));
  return xiks;
}
//------Ik1 following Luo-Rudy formulation (from Shannon model) ------
double CActionPotential::comp_ik1(void)
{
  double ek = (1.0/frt)*log(xko/xki);// K reversal potential
  const double gki=(sqrt(xko/5.4));
  double aki=1.02/(1.0+exp(0.2385*(v-ek-59.215)));
  double bki=(0.49124*exp(0.08032*(v-ek+5.476))+exp(0.061750*(v-ek-594.31)))/(1.0+exp(-0.5143*(v-ek+4.753)));
  double xkin=aki/(aki+bki);
  double xik1=gk1*gki*xkin*(v-ek);
  return xik1;
}
//------- Ito slow following Shannon et. al. 2005 -----------
//------- Ito fast following Shannon et. al. 2005 -----------
double CActionPotential::comp_ito(void)
{
  double ek = (1.0/frt)*log(xko/xki);// K reversal potential
  double rt1=-(v+3.0)/15.0;
  double rt2=(v+33.5)/10.0;
  double rt3=(v+60.0)/10.0;
  double xtos_inf=1.0/(1.0+exp(rt1));
  double ytos_inf=1.0/(1.0+exp(rt2));
  double rs_inf=1.0/(1.0+exp(rt2));
  double trs=((2800-500)/(1+exp((v+60)/10)))+220+500;
  double txs=9.0/(1.0+exp(-rt1)) + 0.5;
  double tys=3000.0/(1.0+exp(rt3)) + 30.0;
  double xitos=gtos*xtos*(0.5*ytos+0.5*rtos)*(v-ek);// ito slow
  xtos = xtos_inf-(xtos_inf-xtos)*exp(-hode/txs);
  ytos = ytos_inf-(ytos_inf-ytos)*exp(-hode/tys);
  rtos = rs_inf-(rs_inf-rtos)*exp(-hode/(trs));

  double xtof_inf=xtos_inf;
  double ytof_inf=ytos_inf;
  double rt4=-(v/30.0)*(v/30.0);
  double rt5=(v+33.5)/10.0;
  double txf=3.5*exp(rt4)+1.5;
  double tyf=20.0/(1.0+exp(rt5))+20.0;
  double xitof=gtof*xtof*ytof*(v-ek);// ito fast
  xtof = xtof_inf-(xtof_inf-xtof)*exp(-hode/txf);
  ytof = ytof_inf-(ytof_inf-ytof)*exp(-hode/tyf);

  _itof=xitof;_itos=xitos;

  return xitos+xitof;
}
// -------Inak (sodium-potassium exchanger) following Shannon --------------
double CActionPotential::comp_inak(void)
{
  const double xkmko=1.5; //these are Inak constants adjusted to fit
  //the experimentally measured dynamic restitution curve
  const double xkmnai=12.0;
  const double sigma = (exp(xnao/67.3)-1.0)/7.0;
  double fnak = 1.0/(1+0.1245*exp(-0.1*v*frt)+0.0365*sigma*exp(-v*frt));
  double xinak = gnak*fnak*(1./(1.+(xkmnai/nai)))*xko/(xko+xkmko);
  return xinak;
}


//
// additional currents from Shannon model
//

// -------Inabk (Na leak current) following Shannon --------------
double CActionPotential::comp_inabk(void)
{
  double ena = (1.0/frt)*log(xnao/nai);
  return gnabk*(v-ena);
}
// -------ICl(Ca) (Ca-dependent Cl current) following Shannon --------------

// Cas dependent
double CActionPotential::comp_icl(double cac)
{
  const double cli=15.0;
  const double clo=150.0;
  const double ecl=(1.0/frt)*log(cli/clo);
  const double kdclca=100;
  return gcl*(v-ecl)/(1+kdclca/cac);
}
// -------Iclbk (background SL Cl flux) following Shannon --------------
double CActionPotential::comp_iclbk(void)
{
  const double cli=15.0;
  const double clo=150.0;
  const double ecl=(1.0/frt)*log(cli/clo);
  return gclbk*(v-ecl);
}

