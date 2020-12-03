#ifndef ___LOGGINGCLASS_H
#define ___LOGGINGCLASS_H
#include <iostream>
#include <fstream>
#include <time.h>
#include <string>
#include <cstdlib>
#include <cmath>
#include <iomanip>
using namespace std;
#ifdef ___USE_MPI
#include <mpi.h>
#endif

#ifndef _WIN32

//link advapi32.lib Pdh.lib User32.lib winmm.lib

#include <unistd.h>
#include <pwd.h>
#include <sys/utsname.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/resource.h>
#else
#include <pdh.h>
#include <windows.h>
#include <mmsystem.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

//for VC warning
//#define _CRT_SECURE_NO_DEPRECATE
#pragma warning ( disable : 4996 )

#define LOG_DATE 0
#define LOG_UNIQUE 1

#define ___HIGHEST_PRIORITY -4
#define ___HIGH_PRIORITY -3
#define ___LOWEST_PRIORITY -2
#define ___LOW_PRIORITY -1


class Logging {
public:
  //      void Alive(void);
  //      void Alive(string str);
  Logging(void);
  #ifndef ___USE_MPI
  Logging(string filename);
  Logging(int i);
  #endif
  ~Logging(void);
  void Note(string note);
  friend ostream &operator<<(ostream &stream, Logging *llog);
  void Wait(int percent = 50);
  int Kill(void);
  int LoadAdjust(int maxload = 0, int minload = 0);
  int Renice(int val = 19);
  void Idleprocess(void);
  void MemoryUsage(void);

private:
  long long int GetTime(char *tim);
  char filename[256 + 1];
  ofstream *logfile;
  //      ofstream *alivefile;
  void Start(void);
  void End(void);
  long long int t1, t2;//cpu time
  long long int a1, a2;//physical time
  long long int t1s, a1s;//cpu and physical time after sleeping
  long long int n1, n2;//physical time for note
  long long int Cputime();
  bool firstnote;

  int node;
  int numofnode;

  int num_procs;
  time_t previous_loadave_time, last_change_time;
  int nextchecktime;
  time_t nextinterval;

  void WriteHostInfo(void);
  void WriteProcessInfo(void);
  void WriteCompilerInfo(void);

  unsigned long cpu_prev[10];
};


#endif /* ___LOGGINGCLASS_H */
