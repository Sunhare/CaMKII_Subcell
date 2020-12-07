#include "recsubcell.h"

CRecSubcell::CRecSubcell(class CSubcell *scx, class CActionPotential *apx, class RECSUBCELLPARAM *param)
{
  sc=scx;
  ap=apx;
  n=sc->n;
  nn=sc->nn;
  nx=sc->nx;
  ny=sc->ny;
  nz=sc->nz;
  nxny=nx*ny;
  nnx=nx*sc->finemesh;
  nny=ny*sc->finemesh;
  nnz=nz*sc->finemesh;
  nnxnny=nnx*nny;
  rec=new char [nn];
  ave=new double [nn];

  step=param->step;

  IsRecciout=IsReccsout=IsReccpout=IsReccnsrout=IsReccjsrout=IsRecryr1out=IsRecryr2out=IsRecryr3out=IsRecryr4out=IsRectvout=IsReccurrentout=false;
  if (param->filenameci!="" && param->filenamecimm!="")
  {
    ciout=new ofstream(param->filenameci.c_str(),ios::out|ios::binary);
    cimmout=new ofstream(param->filenamecimm.c_str());
    IsRecciout=true;
    cifilename=param->filenameci;
    cimmfilename=param->filenamecimm;
  }
  if (param->filenamecs!="" && param->filenamecsmm!="")
  {
    csout=new ofstream(param->filenamecs.c_str(),ios::out|ios::binary);
    csmmout=new ofstream(param->filenamecsmm.c_str());
    IsReccsout=true;
    csfilename=param->filenamecs;
    csmmfilename=param->filenamecsmm;
  }
  if (param->filenamecp!="" && param->filenamecpmm!="")
  {
    cpout=new ofstream(param->filenamecp.c_str(),ios::out|ios::binary);
    cpmmout=new ofstream(param->filenamecpmm.c_str());
    IsReccpout=true;
    cpfilename=param->filenamecp;
    cpmmfilename=param->filenamecpmm;
  }
  if (param->filenamecnsr!="" && param->filenamecnsrmm!="")
  {
    cnsrout=new ofstream(param->filenamecnsr.c_str(),ios::out|ios::binary);
    cnsrmmout=new ofstream(param->filenamecnsrmm.c_str());
    IsReccnsrout=true;
    cnsrfilename=param->filenamecnsr;
    cnsrmmfilename=param->filenamecnsrmm;
  }
  if (param->filenamecjsr!="" && param->filenamecjsrmm!="")
  {
    cjsrout=new ofstream(param->filenamecjsr.c_str(),ios::out|ios::binary);
    cjsrmmout=new ofstream(param->filenamecjsrmm.c_str());
    IsReccjsrout=true;
    cjsrfilename=param->filenamecjsr;
    cjsrmmfilename=param->filenamecjsrmm;
  }
  if (param->filenameryr1!="")
  {
    ryr1out=new ofstream(param->filenameryr1.c_str(),ios::out|ios::binary);
    IsRecryr1out=true;
    ryr1filename=param->filenameryr1;
  }
  if (param->filenameryr2!="")
  {
    ryr2out=new ofstream(param->filenameryr2.c_str(),ios::out|ios::binary);
    IsRecryr2out=true;
    ryr2filename=param->filenameryr2;
  }
  if (param->filenameryr3!="")
  {
    ryr3out=new ofstream(param->filenameryr3.c_str(),ios::out|ios::binary);
    IsRecryr3out=true;
    ryr3filename=param->filenameryr3;
  }
  if (param->filenameryr4!="")
  {
    ryr4out=new ofstream(param->filenameryr4.c_str(),ios::out|ios::binary);
    IsRecryr4out=true;
    ryr4filename=param->filenameryr4;
  }
  if (param->filenametv!="")
  {
    tvout=new ofstream(param->filenametv.c_str());
    IsRectvout=true;
  }
  if (param->filenamecurrent!="")
  {
    currentout=new ofstream(param->filenamecurrent.c_str());
    IsReccurrentout=true;
  }
}

CRecSubcell::CRecSubcell(class CCell *cellx, class RECSUBCELLPARAM *param)
{
  cell=cellx;
  sc=cellx->sc;
  ap=cellx->ap;
  n=sc->n;
  nn=sc->nn;
  rec=new char [nn];

  IsRecciout=IsReccnsrout=IsRectvout=IsReccurrentout=false;
  if (param->filenameci!="" && param->filenamecimm!="")
  {
    ciout=new ofstream(param->filenameci.c_str(),ios::out|ios::binary);
    cimmout=new ofstream(param->filenamecimm.c_str());
    IsRecciout=true;
    cifilename=param->filenameci;
    cimmfilename=param->filenamecimm;
  }
  if (param->filenamecnsr!="" && param->filenamecnsrmm!="")
  {
    cnsrout=new ofstream(param->filenamecnsr.c_str(),ios::out|ios::binary);
    cnsrmmout=new ofstream(param->filenamecnsrmm.c_str());
    IsReccnsrout=true;
  }
  if (param->filenametv!="")
  {
    tvout=new ofstream(param->filenametv.c_str());
    IsRectvout=true;
  }
  if (param->filenamecurrent!="")
  {
    currentout=new ofstream(param->filenamecurrent.c_str());
    IsReccurrentout=true;
  }
}
CRecSubcell::CRecSubcell(void)
{
}

void CRecSubcell::setcell(class CCell *cellx, class RECSUBCELLPARAM *param)
{
  cell=cellx;
  sc=cellx->sc;
  ap=cellx->ap;
  n=sc->n;
  nn=sc->nn;
  rec=new char [nn];

  IsRecciout=IsReccnsrout=IsRectvout=IsReccurrentout=false;
  if (param->filenameci!="" && param->filenamecimm!="")
  {
    ciout=new ofstream(param->filenameci.c_str(),ios::out|ios::binary);
    cimmout=new ofstream(param->filenamecimm.c_str());
    IsRecciout=true;
    cifilename=param->filenameci;
    cimmfilename=param->filenamecimm;
  }
  if (param->filenamecnsr!="" && param->filenamecnsrmm!="")
  {
    cnsrout=new ofstream(param->filenamecnsr.c_str(),ios::out|ios::binary);
    cnsrmmout=new ofstream(param->filenamecnsrmm.c_str());
    IsReccnsrout=true;
  }
  if (param->filenametv!="")
  {
    tvout=new ofstream(param->filenametv.c_str());
    IsRectvout=true;
  }
  if (param->filenamecurrent!="")
  {
    currentout=new ofstream(param->filenamecurrent.c_str());
    IsReccurrentout=true;
  }
}



CRecSubcell::~CRecSubcell()
{
  delete [] rec;
  delete [] ave;
  char temp[255];
  if (IsRectvout){tvout->close();delete tvout;}
  if (IsRecciout){
    ciout->close();cimmout->close();delete ciout;delete cimmout;
#ifdef _WIN32
    sprintf_s(temp,"xz -9f %s",cifilename.c_str());
#else
    sprintf(temp,"xz -9f %s",cifilename.c_str());
#endif
    int res=system(temp);if(res)cout<<res<<endl;
  }
  if (IsReccsout){
    csout->close();csmmout->close();delete csout;delete csmmout;
#ifdef _WIN32
    sprintf_s(temp,"xz -9f %s",csfilename.c_str());int res=system(temp);if(res)cout<<res<<endl;
#else
    sprintf(temp,"xz -9f %s",csfilename.c_str());int res=system(temp);if(res)cout<<res<<endl;
#endif
  }
  if (IsReccpout){
    cpout->close();cpmmout->close();delete cpout;delete cpmmout;
#ifdef _WIN32
    sprintf_s(temp,"xz -9f %s",cpfilename.c_str());int res=system(temp);if(res)cout<<res<<endl;
#else
    sprintf(temp,"xz -9f %s",cpfilename.c_str());int res=system(temp);if(res)cout<<res<<endl;
#endif
  }
  if (IsReccnsrout){
    cnsrout->close();cnsrmmout->close();delete cnsrout;delete cnsrmmout;
#ifdef _WIN32
    sprintf_s(temp,"xz -9f %s",cnsrfilename.c_str());int res=system(temp);if(res)cout<<res<<endl;
#else
    sprintf(temp,"xz -9f %s",cnsrfilename.c_str());int res=system(temp);if(res)cout<<res<<endl;
#endif
  }
  if (IsReccjsrout){
    cjsrout->close();cjsrmmout->close();delete cjsrout;delete cjsrmmout;
#ifdef _WIN32
    sprintf_s(temp,"xz -9f %s",cjsrfilename.c_str());int res=system(temp);if(res)cout<<res<<endl;
#else
    sprintf(temp,"xz -9f %s",cjsrfilename.c_str());int res=system(temp);if(res)cout<<res<<endl;
#endif
  }
  if (IsRecryr1out){
    ryr1out->close();delete ryr1out;
#ifdef _WIN32
    sprintf_s(temp,"xz -9f %s",ryr1filename.c_str());int res=system(temp);if(res)cout<<res<<endl;
#else
    sprintf(temp,"xz -9f %s",ryr1filename.c_str());int res=system(temp);if(res)cout<<res<<endl;
#endif
  }
  if (IsRecryr2out){
    ryr2out->close();delete ryr2out;
#ifdef _WIN32
    sprintf_s(temp,"xz -9f %s",ryr2filename.c_str());int res=system(temp);if(res)cout<<res<<endl;
#else
    sprintf(temp,"xz -9f %s",ryr2filename.c_str());int res=system(temp);if(res)cout<<res<<endl;
#endif
  }
  if (IsRecryr3out){
    ryr3out->close();delete ryr3out;
#ifdef _WIN32
    sprintf_s(temp,"xz -9f %s",ryr3filename.c_str());int res=system(temp);if(res)cout<<res<<endl;
#else
    sprintf(temp,"xz -9f %s",ryr3filename.c_str());int res=system(temp);if(res)cout<<res<<endl;
#endif
  }
  if (IsRecryr4out){
    ryr4out->close();delete ryr4out;
#ifdef _WIN32
    sprintf_s(temp,"xz -9f %s",ryr4filename.c_str());int res=system(temp);if(res)cout<<res<<endl;
#else
    sprintf(temp,"xz -9f %s",ryr4filename.c_str());int res=system(temp);if(res)cout<<res<<endl;
#endif
  }
  if (IsReccurrentout){currentout->close();delete currentout;}
}
void CRecSubcell::recci(int layer)
{
#ifndef CPU
  sc->writetohost(sc->ci,sc->d_ci,nn);
#endif
  if (IsRecciout)
  {
    double cimin=FLT_MAX;
    double cimax=-FLT_MAX;
    if (layer==ALL_LAYER)
    {
      for (int i=0;i<nn;i++)
      {
        if (sc->ci[i]<cimin)cimin=sc->ci[i];
        if (sc->ci[i]>cimax)cimax=sc->ci[i];
      }
      double cidiff=cimax-cimin;
      if (cidiff<0.01)
      {
        cimax=cimin+0.01;
        cidiff=0.01;
      }
      *cimmout<<cimax<<"\t"<<cimin<<endl;
      for (int i=0;i<nn;i++)
      {
        rec[i]=(unsigned char)((sc->ci[i]-cimin)*(255/cidiff));
      }
      ciout->write( (char * ) rec, sizeof(unsigned char)*(nn) );
    }
    else if (layer==Z_LAYER_AVERAGE){
      int l=nnxnny/step/step;
      for (int i=0;i<l;i++)ave[i]=0;
      l=0;
      for (int j=0;j<nny;j+=step){
        for (int i=0;i<nnx;i+=step){
#pragma ivdep
#pragma vector always
          for (int k=0;k<nnz;k++){
            ave[l]+=sc->ci[i+j*nnx+k*nnxnny];
          }
          l++;
        }
      }
      for (int i=0;i<l;i++)ave[i]/=nnz;

#pragma ivdep
#pragma vector always
      for (int i=0;i<l;i++){
        if (ave[i]<cimin)cimin=ave[i];
        if (ave[i]>cimax)cimax=ave[i];
      }


      double cidiff=cimax-cimin;
      if (cidiff<0.01){
        cimax=cimin+0.01;
        cidiff=0.01;
      }
      *cimmout<<cimax<<"\t"<<cimin<<endl;
      for (int i=0;i<l;i++){
        rec[i]=(unsigned char)((ave[i]-cimin)*(255/cidiff));
      }
      ciout->write( (char * ) rec, sizeof(unsigned char)*(l) );
    }
    else{
      int k=layer;
      int l=0;
      for (int j=0;j<nny;j+=step){
        for (int i=0;i<nnx;i+=step){
          ave[l++]=sc->ci[i+j*nnx+k*nnxnny];
        }
      }
#pragma ivdep
#pragma vector always
      for (int i=0;i<l;i++){
        if (ave[i]<cimin)cimin=ave[i];
        if (ave[i]>cimax)cimax=ave[i];
      }
      double cidiff=cimax-cimin;
      if (cidiff<0.01){
        cimax=cimin+0.01;
        cidiff=0.01;
      }
      *cimmout<<cimax<<"\t"<<cimin<<endl;
      for (int i=0;i<l;i++){
        rec[i]=(unsigned char)((ave[i]-cimin)*(255/cidiff));
      }
      ciout->write( (char * ) rec, sizeof(unsigned char)*(l) );
    }
    ciout->flush();
  }
}
void CRecSubcell::reccs(int layer)
{
#ifndef CPU
  sc->writetohost(sc->cs,sc->d_cs,nn);
#endif
  if (IsReccsout)
  {
    double csmin=FLT_MAX;
    double csmax=-FLT_MAX;
    if (layer==ALL_LAYER)
    {
      for (int i=0;i<nn;i++)
      {
        if (sc->cs[i]<csmin)csmin=sc->cs[i];
        if (sc->cs[i]>csmax)csmax=sc->cs[i];
      }
      double csdiff=csmax-csmin;
      if (csdiff<0.01)
      {
        csmax=csmin+0.01;
        csdiff=0.01;
      }
      *csmmout<<csmax<<"\t"<<csmin<<endl;
      for (int i=0;i<nn;i++)
      {
        rec[i]=(unsigned char)((sc->cs[i]-csmin)*(255/csdiff));
      }
      csout->write( (char * ) rec, sizeof(unsigned char)*(nn) );
    }
    else if (layer==Z_LAYER_AVERAGE){
      int l=nnxnny/step/step;
      for (int i=0;i<l;i++)ave[i]=0;
      l=0;
      for (int j=0;j<nny;j+=step){
        for (int i=0;i<nnx;i+=step){
#pragma ivdep
#pragma vector always
          for (int k=0;k<nnz;k++){
            ave[l]+=sc->cs[i+j*nnx+k*nnxnny];
          }
          l++;
        }
      }
      for (int i=0;i<l;i++)ave[i]/=nnz;

#pragma ivdep
#pragma vector always
      for (int i=0;i<l;i++){
        if (ave[i]<csmin)csmin=ave[i];
        if (ave[i]>csmax)csmax=ave[i];
      }


      double csdiff=csmax-csmin;
      if (csdiff<0.01){
        csmax=csmin+0.01;
        csdiff=0.01;
      }
      *csmmout<<csmax<<"\t"<<csmin<<endl;
      for (int i=0;i<l;i++){
        rec[i]=(unsigned char)((ave[i]-csmin)*(255/csdiff));
      }
      csout->write( (char * ) rec, sizeof(unsigned char)*(l) );
    }
    else{
      int k=layer;
      int l=0;
      for (int j=0;j<nny;j+=step){
        for (int i=0;i<nnx;i+=step){
          ave[l++]=sc->cs[i+j*nnx+k*nnxnny];
        }
      }
#pragma ivdep
#pragma vector always
      for (int i=0;i<l;i++){
        if (ave[i]<csmin)csmin=ave[i];
        if (ave[i]>csmax)csmax=ave[i];
      }
      double csdiff=csmax-csmin;
      if (csdiff<0.01){
        csmax=csmin+0.01;
        csdiff=0.01;
      }
      *csmmout<<csmax<<"\t"<<csmin<<endl;
      for (int i=0;i<l;i++){
        rec[i]=(unsigned char)((ave[i]-csmin)*(255/csdiff));
      }
      csout->write( (char * ) rec, sizeof(unsigned char)*(l) );
    }
    csout->flush();
  }
}
void CRecSubcell::reccp(int layer)
{
#ifndef CPU
  sc->writetohost(sc->cp,sc->d_cp,nn);
#endif
  if (IsReccpout)
  {
    double cpmin=FLT_MAX;
    double cpmax=-FLT_MAX;
    if (layer==ALL_LAYER)
    {
      for (int i=0;i<n;i++)
      {
        if (sc->cp[i]<cpmin)cpmin=sc->cp[i];
        if (sc->cp[i]>cpmax)cpmax=sc->cp[i];
      }
      double cpdiff=cpmax-cpmin;
      if (cpdiff<0.01)
      {
        cpmax=cpmin+0.01;
        cpdiff=0.01;
      }
      *cpmmout<<cpmax<<"\t"<<cpmin<<endl;
      for (int i=0;i<n;i++)
      {
        rec[i]=(unsigned char)((sc->cp[i]-cpmin)*(255/cpdiff));
      }
      cpout->write( (char * ) rec, sizeof(unsigned char)*(n) );
    }
    else if (layer==Z_LAYER_AVERAGE){
      int l=nxny/step/step;
      for (int i=0;i<l;i++)ave[i]=0;
      l=0;
      for (int j=0;j<ny;j+=step){
        for (int i=0;i<nx;i+=step){
#pragma ivdep
#pragma vector always
          for (int k=0;k<nz;k++){
            ave[l]+=sc->cp[i+j*nx+k*nxny];
          }
          l++;
        }
      }
      for (int i=0;i<l;i++)ave[i]/=nz;

#pragma ivdep
#pragma vector always
      for (int i=0;i<l;i++){
        if (ave[i]<cpmin)cpmin=ave[i];
        if (ave[i]>cpmax)cpmax=ave[i];
      }

      double cpdiff=cpmax-cpmin;
      if (cpdiff<0.01){
        cpmax=cpmin+0.01;
        cpdiff=0.01;
      }
      *cpmmout<<cpmax<<"\t"<<cpmin<<endl;
      for (int i=0;i<l;i++){
        rec[i]=(unsigned char)((ave[i]-cpmin)*(255/cpdiff));
      }
      cpout->write( (char * ) rec, sizeof(unsigned char)*(l) );
    }
    else{
      int k=layer;
      int l=0;
      for (int j=0;j<ny;j+=step){
        for (int i=0;i<nx;i+=step){
          ave[l++]=sc->cp[i+j*nx+k*nxny];
        }
      }
#pragma ivdep
#pragma vector always
      for (int i=0;i<l;i++){
        if (ave[i]<cpmin)cpmin=ave[i];
        if (ave[i]>cpmax)cpmax=ave[i];
      }
      double cpdiff=cpmax-cpmin;
      if (cpdiff<0.01){
        cpmax=cpmin+0.01;
        cpdiff=0.01;
      }
      *cpmmout<<cpmax<<"\t"<<cpmin<<endl;
      for (int i=0;i<l;i++){
        rec[i]=(unsigned char)((ave[i]-cpmin)*(255/cpdiff));
      }
      cpout->write( (char * ) rec, sizeof(unsigned char)*(l) );
    }
    cpout->flush();
  }
}
void CRecSubcell::reccnsr(int layer)
{
#ifndef CPU
  sc->writetohost(sc->cnsr,sc->d_cnsr,nn);
#endif
  if (IsReccnsrout)
  {
    double cnsrmin=FLT_MAX;
    double cnsrmax=-FLT_MAX;
    if (layer==ALL_LAYER)
    {
      for (int i=0;i<nn;i++)
      {
        if (sc->cnsr[i]<cnsrmin)cnsrmin=sc->cnsr[i];
        if (sc->cnsr[i]>cnsrmax)cnsrmax=sc->cnsr[i];
      }
      double cnsrdiff=cnsrmax-cnsrmin;
      if (cnsrdiff<0.01)
      {
        cnsrmax=cnsrmin+0.01;
        cnsrdiff=0.01;
      }
      *cnsrmmout<<cnsrmax<<"\t"<<cnsrmin<<endl;
      for (int i=0;i<nn;i++)
      {
        rec[i]=(unsigned char)((sc->cnsr[i]-cnsrmin)*(255/cnsrdiff));
      }
      cnsrout->write( (char * ) rec, sizeof(unsigned char)*(nn) );
    }
    else if (layer==Z_LAYER_AVERAGE){
      int l=nnxnny/step/step;
      for (int i=0;i<l;i++)ave[i]=0;
      l=0;
      for (int j=0;j<nny;j+=step){
        for (int i=0;i<nnx;i+=step){
#pragma ivdep
#pragma vector always
          for (int k=0;k<nnz;k++){
            ave[l]+=sc->cnsr[i+j*nnx+k*nnxnny];
          }
          l++;
        }
      }
      for (int i=0;i<l;i++)ave[i]/=nnz;

#pragma ivdep
#pragma vector always
      for (int i=0;i<l;i++){
        if (ave[i]<cnsrmin)cnsrmin=ave[i];
        if (ave[i]>cnsrmax)cnsrmax=ave[i];
      }


      double cnsrdiff=cnsrmax-cnsrmin;
      if (cnsrdiff<0.01){
        cnsrmax=cnsrmin+0.01;
        cnsrdiff=0.01;
      }
      *cnsrmmout<<cnsrmax<<"\t"<<cnsrmin<<endl;
      for (int i=0;i<l;i++){
        rec[i]=(unsigned char)((ave[i]-cnsrmin)*(255/cnsrdiff));
      }
      cnsrout->write( (char * ) rec, sizeof(unsigned char)*(l) );
    }
    else{
      int k=layer;
      int l=0;
      for (int j=0;j<nny;j+=step){
        for (int i=0;i<nnx;i+=step){
          ave[l++]=sc->cnsr[i+j*nnx+k*nnxnny];
        }
      }
#pragma ivdep
#pragma vector always
      for (int i=0;i<l;i++){
        if (ave[i]<cnsrmin)cnsrmin=ave[i];
        if (ave[i]>cnsrmax)cnsrmax=ave[i];
      }
      double cnsrdiff=cnsrmax-cnsrmin;
      if (cnsrdiff<0.01){
        cnsrmax=cnsrmin+0.01;
        cnsrdiff=0.01;
      }
      *cnsrmmout<<cnsrmax<<"\t"<<cnsrmin<<endl;
      for (int i=0;i<l;i++){
        rec[i]=(unsigned char)((ave[i]-cnsrmin)*(255/cnsrdiff));
      }
      cnsrout->write( (char * ) rec, sizeof(unsigned char)*(l) );
    }
    cnsrout->flush();
  }
}
void CRecSubcell::reccjsr(int layer)
{
#ifndef CPU
  sc->writetohost(sc->cjsr,sc->d_cjsr,nn);
#endif
  if (IsReccjsrout)
  {
    double cjsrmin=FLT_MAX;
    double cjsrmax=-FLT_MAX;
    if (layer==ALL_LAYER)
    {
      for (int i=0;i<n;i++)
      {
        if (sc->cjsr[i]<cjsrmin)cjsrmin=sc->cjsr[i];
        if (sc->cjsr[i]>cjsrmax)cjsrmax=sc->cjsr[i];
      }
      double cjsrdiff=cjsrmax-cjsrmin;
      if (cjsrdiff<0.01)
      {
        cjsrmax=cjsrmin+0.01;
        cjsrdiff=0.01;
      }
      *cjsrmmout<<cjsrmax<<"\t"<<cjsrmin<<endl;
      for (int i=0;i<n;i++)
      {
        rec[i]=(unsigned char)((sc->cjsr[i]-cjsrmin)*(255/cjsrdiff));
      }
      cjsrout->write( (char * ) rec, sizeof(unsigned char)*(n) );
    }
    else if (layer==Z_LAYER_AVERAGE){
      int l=nxny/step/step;
      for (int i=0;i<l;i++)ave[i]=0;
      l=0;
      for (int j=0;j<ny;j+=step){
        for (int i=0;i<nx;i+=step){
#pragma ivdep
#pragma vector always
          for (int k=0;k<nz;k++){
            ave[l]+=sc->cjsr[i+j*nx+k*nxny];
          }
          l++;
        }
      }
      for (int i=0;i<l;i++)ave[i]/=nz;

#pragma ivdep
#pragma vector always
      for (int i=0;i<l;i++){
        if (ave[i]<cjsrmin)cjsrmin=ave[i];
        if (ave[i]>cjsrmax)cjsrmax=ave[i];
      }

      double cjsrdiff=cjsrmax-cjsrmin;
      if (cjsrdiff<0.01){
        cjsrmax=cjsrmin+0.01;
        cjsrdiff=0.01;
      }
      *cjsrmmout<<cjsrmax<<"\t"<<cjsrmin<<endl;
      for (int i=0;i<l;i++){
        rec[i]=(unsigned char)((ave[i]-cjsrmin)*(255/cjsrdiff));
      }
      cjsrout->write( (char * ) rec, sizeof(unsigned char)*(l) );
    }
    else{
      int k=layer;
      int l=0;
      for (int j=0;j<ny;j+=step){
        for (int i=0;i<nx;i+=step){
          ave[l++]=sc->cjsr[i+j*nx+k*nxny];
        }
      }
#pragma ivdep
#pragma vector always
      for (int i=0;i<l;i++){
        if (ave[i]<cjsrmin)cjsrmin=ave[i];
        if (ave[i]>cjsrmax)cjsrmax=ave[i];
      }
      double cjsrdiff=cjsrmax-cjsrmin;
      if (cjsrdiff<0.01){
        cjsrmax=cjsrmin+0.01;
        cjsrdiff=0.01;
      }
      *cjsrmmout<<cjsrmax<<"\t"<<cjsrmin<<endl;
      for (int i=0;i<l;i++){
        rec[i]=(unsigned char)((ave[i]-cjsrmin)*(255/cjsrdiff));
      }
      cjsrout->write( (char * ) rec, sizeof(unsigned char)*(l) );
    }
    cjsrout->flush();
  }
}
#ifndef ___DETERMINISTIC
//nryr must be less than 255
void CRecSubcell::recryr1(int layer)
{
#ifndef CPU
  sc->writetohost(sc->ryr1,sc->d_ryr1,n);
#endif
  if (IsRecryr1out)
  {
    if (layer==ALL_LAYER)
    {
      for (int i=0;i<n;i++){
        rec[i]=(unsigned char)(sc->ryr1[i]);
      }
      ryr1out->write( (char * ) rec, sizeof(unsigned char)*(n) );
    }
    else if (layer==Z_LAYER_AVERAGE){
      int l=nxny/step/step;
      for (int i=0;i<l;i++)ave[i]=0;
      l=0;
      for (int j=0;j<ny;j+=step){
        for (int i=0;i<nx;i+=step){
#pragma ivdep
#pragma vector always
          for (int k=0;k<nz;k++){
            ave[l]+=sc->ryr1[i+j*nx+k*nxny];
          }
          l++;
        }
      }

      for (int i=0;i<l;i++){
        rec[i]=(unsigned char)(ave[i]/nz);
      }
      ryr1out->write( (char * ) rec, sizeof(unsigned char)*(l) );
    }
    else{
      int l=0;
      int k=layer;
      for (int j=0;j<ny;j+=step){
#pragma ivdep
#pragma vector always
        for (int i=0;i<nx;i+=step)
        {
          rec[l++]=(unsigned char)(sc->ryr1[i+j*nx+k*nxny]);
        }
      }
      ryr1out->write( (char * ) rec, sizeof(unsigned char)*(nxny) );
    }
    ryr1out->flush();
  }
}
void CRecSubcell::recryr2(int layer)
{
#ifndef CPU
  sc->writetohost(sc->ryr12,sc->d_ryr2,n);
#endif
  if (IsRecryr2out)
  {
    if (layer==ALL_LAYER)
    {
      for (int i=0;i<n;i++){
        rec[i]=(unsigned char)(sc->ryr2[i]);
      }
      ryr2out->write( (char * ) rec, sizeof(unsigned char)*(n) );
    }
    else if (layer==Z_LAYER_AVERAGE){
      int l=nxny/step/step;
      for (int i=0;i<l;i++)ave[i]=0;
      l=0;
      for (int j=0;j<ny;j+=step){
        for (int i=0;i<nx;i+=step){
#pragma ivdep
#pragma vector always
          for (int k=0;k<nz;k++){
            ave[l]+=sc->ryr2[i+j*nx+k*nxny];
          }
          l++;
        }
      }

      for (int i=0;i<l;i++){
        rec[i]=(unsigned char)(ave[i]/nz);
      }
      ryr2out->write( (char * ) rec, sizeof(unsigned char)*(l) );
    }
    else{
      int l=0;
      int k=layer;
      for (int j=0;j<ny;j+=step){
#pragma ivdep
#pragma vector always
        for (int i=0;i<nx;i+=step)
        {
          rec[l++]=(unsigned char)(sc->ryr2[i+j*nx+k*nxny]);
        }
      }
      ryr2out->write( (char * ) rec, sizeof(unsigned char)*(nxny) );
    }
    ryr2out->flush();
  }
}
void CRecSubcell::recryr3(int layer)
{
#ifndef CPU
  sc->writetohost(sc->ryr3,sc->d_ryr3,n);
#endif
  if (IsRecryr3out)
  {
    if (layer==ALL_LAYER)
    {
      for (int i=0;i<n;i++){
        rec[i]=(unsigned char)(sc->ryr3[i]);
      }
      ryr3out->write( (char * ) rec, sizeof(unsigned char)*(n) );
    }
    else if (layer==Z_LAYER_AVERAGE){
      int l=nxny/step/step;
      for (int i=0;i<l;i++)ave[i]=0;
      l=0;
      for (int j=0;j<ny;j+=step){
        for (int i=0;i<nx;i+=step){
#pragma ivdep
#pragma vector always
          for (int k=0;k<nz;k++){
            ave[l]+=sc->ryr3[i+j*nx+k*nxny];
          }
          l++;
        }
      }

      for (int i=0;i<l;i++){
        rec[i]=(unsigned char)(ave[i]/nz);
      }
      ryr3out->write( (char * ) rec, sizeof(unsigned char)*(l) );
    }
    else{
      int l=0;
      int k=layer;
      for (int j=0;j<ny;j+=step){
#pragma ivdep
#pragma vector always
        for (int i=0;i<nx;i+=step)
        {
          rec[l++]=(unsigned char)(sc->ryr3[i+j*nx+k*nxny]);
        }
      }
      ryr3out->write( (char * ) rec, sizeof(unsigned char)*(nxny) );
    }
    ryr3out->flush();
  }
}
void CRecSubcell::recryr4(int layer)
{
#ifndef CPU
  sc->writetohost(sc->ryr1,sc->d_ryr1,n);
  sc->writetohost(sc->ryr2,sc->d_ryr2,n);
  sc->writetohost(sc->ryr3,sc->d_ryr3,n);
#endif
  if (IsRecryr4out)
  {
    if (layer==ALL_LAYER)
    {
      for (int i=0;i<n;i++){
        rec[i]=(unsigned char)(sc->nryr[i]-sc->ryr1[i]-sc->ryr2[i]-sc->ryr3[i]);
      }
      ryr4out->write( (char * ) rec, sizeof(unsigned char)*(n) );
    }
    else if (layer==Z_LAYER_AVERAGE){
      int l=nxny/step/step;
      for (int i=0;i<l;i++)ave[i]=0;
      l=0;
      for (int j=0;j<ny;j+=step){
        for (int i=0;i<nx;i+=step){
#pragma ivdep
#pragma vector always
          for (int k=0;k<nz;k++){
            ave[l]+=sc->nryr[i+j*nx+k*nxny]-sc->ryr1[i+j*nx+k*nxny]-sc->ryr2[i+j*nx+k*nxny]-sc->ryr3[i+j*nx+k*nxny];
          }
          l++;
        }
      }

      for (int i=0;i<l;i++){
        rec[i]=(unsigned char)(ave[i]/nz);
      }
      ryr4out->write( (char * ) rec, sizeof(unsigned char)*(l) );
    }
    else{
      int l=0;
      int k=layer;
      for (int j=0;j<ny;j+=step){
#pragma ivdep
#pragma vector always
        for (int i=0;i<nx;i+=step)
        {
          rec[l++]=(unsigned char)(sc->nryr[i+j*nx+k*nxny]-sc->ryr1[i+j*nx+k*nxny]-sc->ryr2[i+j*nx+k*nxny]-sc->ryr3[i+j*nx+k*nxny]);
        }
      }
      ryr4out->write( (char * ) rec, sizeof(unsigned char)*(nxny) );
    }
    ryr4out->flush();
  }
}
#endif
void CRecSubcell::reccurrent(double t)
{
  if (IsReccurrentout)
    *currentout<<t<<"\t"<<sc->icaave<<"\t"<<sc->incxave<<"\t"<<sc->iupave<<"\t"<<sc->irave<<"\t"<<sc->ileakave<<endl;
}
void CRecSubcell::rectv(double t)
{
  if (IsRectvout)
    *tvout<<t<<"\t"<<ap->v<<"\t"<<computeaveci()<<"\t"<<computeavecs()<<"\t"<<computeavecp()<<"\t"<<computeavecnsr()<<"\t"<<computeavecjsr()<<"\t"<<ap->nai<<endl;
}

double CRecSubcell::computeaveci(void)
{
#ifndef CPU
  sc->writetohost(sc->ci,sc->d_ci,nn);
#endif
  double sum=0;
  for (int i=0;i<nn;i++)
    sum+=sc->ci[i];
  return (sum/nn);
}
double CRecSubcell::computeavecs(void)
{
#ifndef CPU
  sc->writetohost(sc->cs,sc->d_cs,nn);
#endif
  double sum=0;
  for (int i=0;i<nn;i++)
    sum+=sc->cs[i];
  return (sum/nn);
}
double CRecSubcell::computeavecp(void)
{
#ifndef CPU
  sc->writetohost(sc->cp,sc->d_cp,n);
#endif
  double sum=0;
  for (int i=0;i<n;i++)
    sum+=sc->cp[i];
  return (sum/n);
}
double CRecSubcell::computeavesr(void) // incorrect     ****
{
  double sum=0;
  for (int i=0;i<n;i++)
    sum+=(sc->vjsr*sc->cjsr[i]+sc->vnsr*sc->cnsr[i]);
  return (sum/(sc->vjsr+sc->vnsr)/n);
}
double CRecSubcell::computeavecjsr(void)
{
#ifndef CPU
  sc->writetohost(sc->cjsr,sc->d_cjsr,n);
#endif
  double sum=0;
  for (int i=0;i<n;i++)
    sum+=sc->cjsr[i];
  return (sum/n);
}
double CRecSubcell::computeavecnsr(void)
{
#ifndef CPU
  sc->writetohost(sc->cnsr,sc->d_cnsr,nn);
#endif
  double sum=0;
  for (int i=0;i<nn;i++)
    sum+=sc->cnsr[i];
  return (sum/nn);
}
