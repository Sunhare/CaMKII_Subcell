#include "log.h"

#ifdef _WIN32
#include <windows.h>
#endif

#ifdef __NVCC__
#include <cuda_runtime.h>
inline int _ConvertSMVer2Cores(int major, int minor) {
  // Defines for GPU Architecture types (using the SM version to determine the # of cores per SM
  typedef struct
  {
    int SM; // 0xMm (hexidecimal notation), M = SM Major version, and m = SM minor version
    int Cores;
  } sSMtoCores;

  sSMtoCores nGpuArchCoresPerSM[] =
  {
    { 0x20, 32 }, // Fermi Generation (SM 2.0) GF100 class
    { 0x21, 48 }, // Fermi Generation (SM 2.1) GF10x class
    { 0x30, 192 }, // Kepler Generation (SM 3.0) GK10x class
    { 0x32, 192 }, // Kepler Generation (SM 3.2) GK10x class
    { 0x35, 192 }, // Kepler Generation (SM 3.5) GK11x class
    { 0x37, 192 }, // Kepler Generation (SM 3.7) GK21x class
    { 0x50, 128 }, // Maxwell Generation (SM 5.0) GM10x class
    { 0x52, 128 }, // Maxwell Generation (SM 5.2) GM20x class
    { -1, -1 }
  };

  int index = 0;

  while (nGpuArchCoresPerSM[index].SM != -1)
  {
    if (nGpuArchCoresPerSM[index].SM == ((major << 4) + minor))
    {
      return nGpuArchCoresPerSM[index].Cores;
    }

    index++;
  }

  // If we don't find the values, we default use the previous one to run properly
  printf("MapSMtoCores for SM %d.%d is undefined.  Default to use %d Cores/SM\n", major, minor, nGpuArchCoresPerSM[index - 1].Cores);
  return nGpuArchCoresPerSM[index - 1].Cores;
}
#endif

static inline char *
skip_token(const char *p) {
  while (isspace(*p)) p++;
  while (*p && !isspace(*p)) p++;
  return (char *)p;
}

ostream &operator<<(ostream &stream, Logging *llog) {
#ifndef _WIN32
  struct passwd *pwdptr;
  pwdptr = getpwuid(getuid());
  stream << "User : " << pwdptr->pw_name << " (UID=" << pwdptr->pw_uid << ", GID=" << pwdptr->pw_gid << ", HOME=" << pwdptr->pw_dir << ")" << endl;
  struct utsname uts;
  uname(&uts);
  stream << "Host name : " << uts.nodename << endl;
  stream << "OS name and version : " << uts.sysname << "   " << uts.release << "  " << uts.version << endl;
  pid_t pid = getpid();
  stream << "PID : " << pid << endl;
  char dir[255];
  if (getcwd(dir, 255))
    stream << "Execute Directory : " << dir << endl;
#endif

#ifdef _OPENMP
  stream << "OpenMP : Enabled" << endl;
#endif

  char ttt[256];
  llog->a2 = llog->GetTime(ttt);
  stream << "Current Time : " << ttt << endl;

#ifndef _WIN32
  llog->t2 = llog->Cputime();
  if ((llog->a2 - llog->a1) / 1000000.0 < 60 * 60) {
    stream << "CPU time      : " << (llog->t2 - llog->t1) / 1000000.0 / 60 << " min (" << (llog->t2 - llog->t1) / 1000000.0 << " sec)" << endl;
    stream << "Physical time : " << (llog->a2 - llog->a1) / 1000000.0 / 60 << " min (" << (llog->a2, llog->a1) / 1000000.0 << " sec)" << endl;
  }
  else if ((llog->a2 - llog->a1) / 1000000.0 < 60 * 60 * 24 * 2) {
    stream << "CPU time      : " << (llog->t2 - llog->t1) / 1000000.0 / (60 * 60) << " hour (" << (llog->t2 - llog->t1) / 1000000.0 << " sec)" << endl;
    stream << "Physical time : " << (llog->a2 - llog->a1) / 1000000.0 / (60 * 60) << " hour (" << (llog->a2 - llog->a1) / 1000000.0 << " sec)" << endl;
  }
  else {
    stream << "CPU time      : " << (llog->t2 - llog->t1) / 1000000.0 / (60 * 60 * 24) << " days (" << (llog->t2 - llog->t1) / 1000000.0 << " sec)" << endl;
    stream << "Physical time : " << (llog->a2 - llog->a1) / 1000000.0 / (60 * 60 * 24) << " days (" << (llog->a2 - llog->a1) / 1000000.0 << " sec)" << endl;
  }
  stream << "CPU usgae     : " << (llog->t2 - llog->t1)*100.0 / (llog->a2, llog->a1) << " %  (priority nice = " << getpriority(PRIO_PROCESS, 0) << ")" << endl;
#endif
  return stream;
}
Logging::Logging(void) {
  node = 0;
  numofnode = 1;
#ifdef ___USE_MPI
  char dummy[255];
  MPI_Comm_rank(MPI_COMM_WORLD, &node);
  MPI_Comm_size(MPI_COMM_WORLD, &numofnode);
#ifdef ___MPI_OLDSTYLE
  sprintf(dummy, "logfile%d.log", node);
  logfile = new ofstream(dummy);
#else
  if (node == 0) {
    logfile = new ofstream("logfile.log");
  }
#endif
#else
  logfile = new ofstream("logfile.log");
#endif
  Start();
}

#ifndef ___USE_MPI
Logging::Logging(string filename) {
  logfile = new ofstream(filename.c_str());
  Start();
}
Logging::Logging(int i) {
  if (i == LOG_DATE) {
    struct tm *tptr;
    time_t now = time(NULL);
    tptr = localtime(&now);
    char dummy[255];
    sprintf(dummy, "logfile%4d%02d%02d%02d%02d%02d.log",
      (tptr->tm_year) + 1900, tptr->tm_mon + 1, tptr->tm_mday,
      tptr->tm_hour, tptr->tm_min, tptr->tm_sec);
    logfile = new ofstream(dummy);
    Start();
  }
  if (i == LOG_UNIQUE) {
    char dummy[255];

    struct tm *tptr;
    time_t now = time(NULL);
    tptr = localtime(&now);
#ifndef _WIN32
    struct utsname uts;
    uname(&uts);
    pid_t pid = getpid();
    sprintf(dummy, "logfile_%4d%02d%02d%02d%02d%02d_%s_%05d.log",
      (tptr->tm_year) + 1900, tptr->tm_mon + 1, tptr->tm_mday,
      tptr->tm_hour, tptr->tm_min, tptr->tm_sec, uts.nodename, pid);
#else
    sprintf(dummy, "logfile%4d%02d%02d%02d%02d%02d.log",
      (tptr->tm_year) + 1900, tptr->tm_mon + 1, tptr->tm_mday,
      tptr->tm_hour, tptr->tm_min, tptr->tm_sec);
#endif
    logfile = new ofstream(dummy);
    Start();
  }
  else {
    logfile = new ofstream("logfile.log");
    Start();
  }
}
#endif

Logging::~Logging(void) {
  End();
#ifdef ___MPI_OLDSTYLE
  delete logfile;
#else
  if (node == 0) {
    delete logfile;
  }
#endif
}
void Logging::Start(void) {
  nextchecktime = 5 + rand() % 55;
  nextinterval = time(NULL);

  FILE *fid = fopen("/proc/stat", "r");
  if (fid) {
    char *p;
    const int len = 256;
    char buffer[len];
    if (fgets(buffer, len, fid)) {
      p = skip_token(buffer);
      for (int j = 0; j < 10; j++) {
        cpu_prev[j] = strtoul(p, &p, 0);
      }
    }
    fclose(fid);
  }


  firstnote = true;
  char ttt[256];
  a1s = n1 = a1 = GetTime(ttt);
  t1s = t1 = Cputime();
#ifdef _OPENMP
  num_procs = omp_get_num_procs();
#endif
  previous_loadave_time = last_change_time = time(NULL);
#ifdef ___USE_MPI

  pid_t pid = getpid();
  struct utsname uts;
  uname(&uts);
  if (node == 0) {
    WriteHostInfo();
    WriteProcessInfo();
    WriteCompilerInfo();
    *logfile << "MPI : Enabled (# of processes = " << numofnode << ")" << endl;
    *logfile << "Headnode : " << uts.nodename << endl;
    *logfile << "Node\thost\tPID\tStart time" << endl;
    *logfile << node << "\t\t" << uts.nodename << "\t\t" << pid << "\t" << ttt << endl;
    for (int i = 1; i < numofnode; i++) {
      MPI_Status status;
      MPI_Recv(ttt, 255, MPI_CHAR, i, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(&pid, sizeof(pid_t), MPI_CHAR, i, 0, MPI_COMM_WORLD, &status);
      MPI_Recv(uts.nodename, sizeof(uts.nodename), MPI_CHAR, i, 0, MPI_COMM_WORLD, &status);
      *logfile << i << "\t\t" << uts.nodename << "\t\t" << pid << "\t" << ttt << endl;
    }
  }
  else {
    MPI_Send(ttt, 255, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&pid, sizeof(pid_t), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    MPI_Send(uts.nodename, sizeof(uts.nodename), MPI_CHAR, 0, 0, MPI_COMM_WORLD);
  }

#else
  *logfile << "Start Time : " << ttt << endl;
  WriteHostInfo();
  WriteProcessInfo();
  WriteCompilerInfo();
#endif
}
/*
void Logging::Alive(void) {
char ttt[256];
GetTime(ttt);
*alivefile<<"alive : "<<ttt<<endl;
}
void Logging::Alive(string str) {
char ttt[256];
a2=GetTime(ttt);
*alivefile<<"alive : "<<str<<" "<<ttt<<"("<< difftime(a2,a1)<<" sec)"<<endl;
}
*/

void Logging::Wait(int percent) {
  int sleeptime = 15;
#ifndef _WIN32
  double numofprocs = 1;
#ifdef _OPENMP
  numofprocs = omp_get_num_procs();
#endif

  FILE *fid = fopen("/proc/stat", "r");
  if (fid)
  {
    double threshold = percent / 100.0;
    //    cp_diff[0]/total; user 0~1
    //    cp_diff[1]/total; niced 0~1
    //    cp_diff[2]/total; system 0~1
    //    cp_diff[3]/total; idle 0~1

    unsigned long cp_time[10];
    unsigned long cp_prev[10];
    unsigned long cp_diff[10];
    char *p;

    const int len = 256;
    char buffer[len];
    if (fgets(buffer, len, fid)) {
      p = skip_token(buffer);
      for (int j = 0; j < 10; j++) {
        cpu_prev[j] = strtoul(p, &p, 0);
      }
    }
    fclose(fid);
    sleep(1);
    fid = fopen("/proc/stat", "r");
    if (fgets(buffer, len, fid))
    {
      p = skip_token(buffer);
      for (int j = 0; j < 10; j++) {
        cp_time[j] = strtoul(p, &p, 0);
      }
    }
    fclose(fid);
    double total = 0;
    for (int j = 0; j < 10; j++) {
      cp_diff[j] = cp_time[j] - cp_prev[j];
      total += cp_diff[j];
    }
    if (cp_diff[3] / total > threshold)return;

    *logfile << "\nwaiting for previous jobs ... (CPU idle: " << cp_diff[3] * 100.0 / total << " %)" << endl;
    ofstream osslp("sleepbeacon.txt");
    do {
      char tttt[256];
      GetTime(tttt);
      osslp << tttt << endl;
      sleep(sleeptime);
      for (int j = 0; j < 10; j++) {
        cp_prev[j] = cp_time[j];
      }
      fid = fopen("/proc/stat", "r");
      if (fgets(buffer, len, fid)) {
        p = skip_token(buffer);
        for (int j = 0; j < 10; j++) {
          cp_time[j] = strtoul(p, &p, 0);
        }
      }
      fclose(fid);
      total = 0;
      for (int j = 0; j < 10; j++) {
        cp_diff[j] = cp_time[j] - cp_prev[j];
        total += cp_diff[j];
      }
    } while (cp_diff[3] / total < threshold);
    remove("sleepbeacon.txt");
    char ttt[256];
    a1s = GetTime(ttt);
    t1s = Cputime();
    *logfile << "CPU is idle more than " << percent << " percent (" << cp_diff[3] * 100.0 / total << " %)\nprogram start : " << ttt << endl;
  }
  else {
    double threshold = numofprocs*percent / 100.0;
    double loadavg[3];
    if (getloadavg(loadavg, 3) == -1) {
      cout << "Can't get load average" << endl;
      return;
    }
    if (loadavg[0]<threshold)return;
    *logfile << "\nwaiting for previous jobs ... (" << loadavg[0] << ")" << endl;
    ofstream osslp("sleepbeacon.txt");
    do {
      char tttt[256];
      GetTime(tttt);
      osslp << tttt << endl;
      sleep(sleeptime);
      getloadavg(loadavg, 3);
    } while (loadavg[0]>threshold);
    remove("sleepbeacon.txt");
    char ttt[256];
    a1s = GetTime(ttt);
    t1s = Cputime();
    *logfile << "load average is less than " << percent << " percent(" << loadavg[0] << ")\nprogram start : " << ttt << endl;
  }
  *logfile << endl;
#else
  HQUERY hQuery;//handle to the Query
  PdhOpenQuery(NULL, 0, &hQuery);
  HCOUNTER hCounter;//Handle to the counter
  PdhAddCounter(hQuery, "\\Processor(_Total)\\% Processor Time", 0, &hCounter);
  PdhCollectQueryData(hQuery);
  Sleep(1000);
  PdhCollectQueryData(hQuery);
  PDH_FMT_COUNTERVALUE fntValue;
  PdhGetFormattedCounterValue(hCounter, PDH_FMT_DOUBLE, NULL, &fntValue);
  if (fntValue.doubleValue > percent) {
    *logfile << "\nwaiting for previous jobs ... (" << fntValue.doubleValue << ")" << endl;
    ofstream osslp("sleepbeacon.txt");
    do {
      char tttt[256];
      GetTime(tttt);
      osslp << tttt << endl;
      Sleep(sleeptime*1000.0);
      PdhCollectQueryData(hQuery);
      PdhGetFormattedCounterValue(hCounter, PDH_FMT_DOUBLE, NULL, &fntValue);
    } while (fntValue.doubleValue > percent);
    remove("sleepbeacon.txt");
    char ttt[256];
    a1s = GetTime(ttt);
    t1s = Cputime();
    *logfile << "load average is less than " << percent << " percent(" << fntValue.doubleValue << ")\nprogram start : " << ttt << endl;
  }

  PdhCloseQuery(hQuery);
#endif
  return;
}
int Logging::Kill(void) {
  //  bool IsAbort=false;
  FILE *fpp = fopen("kill.txt", "r");
  if (fpp != NULL) {
    //      IsAbort=true;
    fclose(fpp);
    char ttt[256];
    GetTime(ttt);
    *logfile << "kill.txt was found : " << ttt << endl;
    *logfile << endl;
    return 1;
  }
#ifndef _WIN32
  FILE *fpp2 = fopen("sleep.txt", "r");
  if (fpp2 != NULL) {
    //      IsAbort=true;
    fclose(fpp2);
    char ttt[256];
    GetTime(ttt);
    *logfile << "sleep.txt was found : " << ttt << endl;
    *logfile << endl;
    remove("sleep.txt");
    ofstream osslp("sleepbeacon.txt");
    double loadavg[3];
    do {
      char tttt[256];
      GetTime(tttt);
      osslp << tttt << endl;
      sleep(10);
      getloadavg(loadavg, 3);
    } while (loadavg[0] > 0.8);
    remove("sleepbeacon.txt");
    GetTime(ttt);
    *logfile << "load average is less than 0.8\nprogram continue : " << ttt << endl;
    *logfile << endl;
  }
#endif
  return 0;
}
void Logging::Idleprocess(void)
{
#ifndef _WIN32
  ofstream osslp("sleepbeacon.txt");
  char tttt[256];
  GetTime(tttt);
  osslp << tttt << endl;
  sleep(30);
  remove("sleepbeacon.txt");
#ifdef _OPENMP
  int currentnum = omp_get_max_threads();
  currentnum /= 2;
  if (currentnum >= 1)
    omp_set_num_threads(currentnum);
  nextchecktime = 120;
  previous_loadave_time = time(NULL);
#endif
#endif
}
int Logging::LoadAdjust(int maxload, int minload) {
  time_t current_time = time(NULL);
  if (difftime(current_time, previous_loadave_time) < nextchecktime)return -1;
  previous_loadave_time = current_time;
  nextchecktime = 5 + rand() % 55;

  bool idle = false;
  if (minload == ___LOWEST_PRIORITY || minload == ___LOW_PRIORITY) {
    if (minload == ___LOWEST_PRIORITY && difftime(current_time, nextinterval) > 3600) {
#ifdef _WIN32
      Sleep(60 * 1000);
#else
      sleep(60);
#endif
      nextinterval = current_time;
    }
    idle = true;
    minload = 0;
  }


  int numofthread = 1;
#ifndef _WIN32
  double loadavg[3];
  getloadavg(loadavg, 3);
#ifdef _OPENMP
  if (maxload <= 0) {
    maxload = num_procs;
  }
  int currentnum = omp_get_max_threads();

  if (minload == ___HIGH_PRIORITY || minload == ___HIGHEST_PRIORITY) {
    if (difftime(current_time, nextinterval) > 3600) {
      if (difftime(current_time, nextinterval) > 3600 + 180) {
        nextinterval = current_time;
      }
      else {
        omp_set_num_threads(maxload);
        return maxload;
      }
    }
    if (minload == ___HIGHEST_PRIORITY)
      minload = 1;
    else
      minload = 0;
  }

  FILE *fid = fopen("/proc/stat", "r");
  if (fid) {
    //    cpu_diff[0]/total; user 0~1
    //    cpu_diff[1]/total; niced 0~1
    //    cpu_diff[2]/total; system 0~1
    //    cpu_diff[3]/total; idle 0~1

    unsigned long cpu_time[10];
    unsigned long cpu_diff[10];

    const int len = 256;
    char buffer[len];
    if (fgets(buffer, len, fid)) {
      char *p = skip_token(buffer);
      for (int j = 0; j < 10; j++) {
        cpu_time[j] = strtoul(p, &p, 0);
      }
    }
    fclose(fid);
    double total = 0;
    for (int j = 0; j < 10; j++) {
      cpu_diff[j] = cpu_time[j] - cpu_prev[j];
      total += cpu_diff[j];
      cpu_prev[j] = cpu_time[j];
    }
    if (cpu_diff[2] / total > 0.02){
      //        cout<<cpu_diff[2]/total<<" take 1"<<endl;
      numofthread = currentnum - 1;

      int overload = int(loadavg[0] + 0.5) - num_procs;
      if (overload > 0) {
        numofthread = numofthread - overload;
      }
      //        cout<<"overload="<<overload<<endl;

    }
    else {
      int overload = int(loadavg[1] + 0.5) - num_procs;
      if (overload > 0) {
        numofthread = numofthread - overload;
      }
      else {
        numofthread = currentnum + ceil(num_procs*(cpu_diff[3] / total - 0.8 / num_procs));
        //        cout<<"current="<<currentnum<<"\t"<<num_procs*(cpu_diff[3]/total-0.8/num_procs)<<"\t"<<cpu_diff[3]/total<<"\t"<<cpu_diff[2]<<"\t"<<cpu_diff[3]<<"\t"<<total<<endl;
      }
    }

    if (numofthread > maxload)numofthread = maxload;
    if (numofthread < minload)numofthread = minload;
    if (numofthread < 1) {
      ofstream osslp("sleepbeacon.txt");
      do {
        nextchecktime = 5 + rand() % 25;
        //                cout<<"nextchecksleeptime="<<nextchecktime<<endl;
        char tttt[256];
        GetTime(tttt);
        osslp << tttt << endl;
        sleep(nextchecktime);
        for (int j = 0; j < 10; j++) {
          cpu_prev[j] = cpu_time[j];
        }
        fid = fopen("/proc/stat", "r");
        if (fgets(buffer, len, fid)) {
          char *p = skip_token(buffer);
          for (int j = 0; j < 10; j++) {
            cpu_time[j] = strtoul(p, &p, 0);
          }
        }
        fclose(fid);
        total = 0;
        for (int j = 0; j < 10; j++) {
          cpu_diff[j] = cpu_time[j] - cpu_prev[j];
          total += cpu_diff[j];
        }
        getloadavg(loadavg, 3);
      } while (cpu_diff[3] / total<1.0 / num_procs || loadavg[0]>num_procs);
      remove("sleepbeacon.txt");
      numofthread = 1;
      previous_loadave_time = time(NULL);
      nextchecktime = 5 + rand() % 55;
    }
  }
  else {
    do{
      double loadavg[3];
      getloadavg(loadavg, 3);
      int overload = int(loadavg[0] + 0.5) - num_procs;
      numofthread = currentnum - overload;
      //      cout<<currentnum<<" "<<loadavg[0]<<" "<<overload<<" "<<numofthread<<endl;
      if (numofthread > maxload)numofthread = maxload;
      if (numofthread < minload)numofthread = minload;
      if (numofthread < 1) {
        //          cout<<"sleep 10 sec"<<endl;
        sleep(10);
        currentnum = 0;
      }
    } while (numofthread < 1);
  }
  if (currentnum != numofthread) {
    //        cout<<"change "<<currentnum<<" to "<<numofthread<<endl;
    if (idle&&numofthread < maxload)
    {
      ofstream osslp("sleepbeacon.txt");
      unsigned long cpu_diff[10];
      double total = 0;
      do
      {
        nextchecktime = 20 + rand() % 60;
        char tttt[256];
        GetTime(tttt);
        osslp << tttt << endl;
        sleep(nextchecktime);

        fid = fopen("/proc/stat", "r");
        const int len = 256;
        char buffer[len];
        unsigned long cpu_time[10];
        if (fgets(buffer, len, fid)) {
          char *p = skip_token(buffer);
          for (int j = 0; j < 10; j++) {
            cpu_time[j] = strtoul(p, &p, 0);
          }
        }
        fclose(fid);
        total = 0;
        for (int j = 0; j < 10; j++) {
          cpu_diff[j] = cpu_time[j] - cpu_prev[j];
          total += cpu_diff[j];
        }
        for (int j = 0; j < 10; j++) {
          cpu_prev[j] = cpu_time[j];
        }
        getloadavg(loadavg, 3);
      } while (cpu_diff[3] / total<(maxload - 0.5) / num_procs || loadavg[0]>num_procs);
      remove("sleepbeacon.txt");
      numofthread = maxload;
      previous_loadave_time = time(NULL);
      nextchecktime = 5 + rand() % 55;
    }
    omp_set_num_threads(numofthread);
  }
#endif
#else
  HANDLE hCurrentProcess = GetCurrentProcess();
  SetPriorityClass(hCurrentProcess, IDLE_PRIORITY_CLASS);
#endif
  return numofthread;
}
int Logging::Renice(int val)
{
#ifdef _WIN32
  HANDLE hCurrentProcess = GetCurrentProcess();
  if (val == 0)
    SetPriorityClass(hCurrentProcess, NORMAL_PRIORITY_CLASS);
  else if (val < 0)
    SetPriorityClass(hCurrentProcess, ABOVE_NORMAL_PRIORITY_CLASS);
  else
    SetPriorityClass(hCurrentProcess, IDLE_PRIORITY_CLASS);
#else
  val = setpriority(PRIO_PROCESS, 0, val);
#endif
  return val;
}

void Logging::End(void)
{
#ifdef ___USE_MPI

  char ttt[256];
  a2 = GetTime(ttt);
  char tmpstring[1025];

  t2 = Cputime();
  if ((a2 - a1) / 1000000.0 < 60 * 60) {
    sprintf(tmpstring, "%s\t%f min (%f sec)\t%f min (%f sec)\t%f (%d)", ttt, (t2 - t1) / 1000000.0 / 60, (t2 - t1) / 1000000.0, (a2 - a1) / 1000000.0 / 60, (a2 - a1) / 1000000.0, (t2 - t1)*100.0 / (a2 - a1), getpriority(PRIO_PROCESS, 0));
  }
  else if ((a2 - a1) / 1000000.0 < 60 * 60 * 24 * 2) {
    sprintf(tmpstring, "%s\t%d hour %f min (%f sec)\t%d hour %f min (%f sec)\t%f (%d)", ttt, int((t2 - t1) / 1000000.0 / (60 * 60)), 60.0*((t2 - t1) / 1000000.0 / (60 * 60) - int((t2 - t1) / 1000000.0 / (60 * 60))), (t2 - t1) / 1000000.0, int((a2 - a1) / 1000000.0 / (60 * 60)), 60.0*((a2 - a1) / 1000000.0 / (60 * 60) - int((a2 - a1) / 1000000.0 / (60 * 60))), (a2 - a1) / 1000000.0, (t2 - t1)*100.0 / (a2 - a1), getpriority(PRIO_PROCESS, 0));
  }
  else {
    sprintf(tmpstring, "%s\t%d days %f hour (%f sec)\t%d days %f hour (%f sec)\t%f (%d)", ttt, int((t2 - t1) / 1000000.0 / (60 * 60 * 24)), 24.0*((t2 - t1) / 1000000.0 / (60 * 60 * 24) - int((t2 - t1) / 1000000.0 / (60 * 60 * 24))), (t2 - t1) / 1000000.0, int((a2 - a1) / 1000000.0 / (60 * 60 * 24)), 24.0*((a2 - a1) / 1000000.0 / (60 * 60 * 24) - int((a2 - a1) / 1000000.0 / (60 * 60 * 24))), (a2 - a1) / 1000000.0, (t2 - t1)*100.0 / (a2 - a1), getpriority(PRIO_PROCESS, 0));
  }

  if (node == 0) {
    *logfile << "\tEnd Time\t\t\t\t\t\t\t\tCPU Time\t\t\t\t\tPhysical time\t\t\t\tCPU usgae (nice)" << endl;
    *logfile << node << "\t" << tmpstring << endl;
    for (int i = 1; i < numofnode; i++) {
      MPI_Status status;
      MPI_Recv(tmpstring, 1024, MPI_CHAR, i, 0, MPI_COMM_WORLD, &status);
      *logfile << i << "\t" << tmpstring << endl;
    }

    FILE *fp = fopen("/proc/cpuinfo", "r");
    if (fp != NULL) {
      *logfile << endl;
      char tmp[255];
      while (fgets(tmp, 255, fp) != NULL) {
        *logfile << tmp;
      }
      fclose(fp);
      *logfile << endl;
    }

#define MAXLINE 1024
    char s[MAXLINE];
    FILE *fpin;
    fpin = popen("sensors", "r");
    while (fgets(s, MAXLINE, fpin) != NULL){
      *logfile << s;
    }
    pclose(fpin);
    *logfile << endl;

    fpin = popen("ls -la --full-time", "r");
    while (fgets(s, MAXLINE, fpin) != NULL){
      *logfile << s;
    }
    pclose(fpin);
    *logfile << endl;
  }
  else {
    MPI_Send(tmpstring, 1024, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
  }
#else
  char ttt[256];
  a2 = GetTime(ttt);
  *logfile << "End Time : " << ttt << endl;
#ifndef _WIN32
  t2 = Cputime();
  if ((a2 - a1) / 1000000.0 < 60) {
    *logfile << "CPU time      : " << (t2 - t1) / 1000000.0 << " sec" << endl;
    *logfile << "Physical time : " << (a2 - a1) / 1000000.0 << " sec" << endl;
  }
  else if ((a2 - a1) / 1000000.0 < 60 * 60) {
    *logfile << "CPU time      : " << (t2 - t1) / 1000000.0 / 60 << " min (" << (t2 - t1) / 1000000.0 << " sec)" << endl;
    *logfile << "Physical time : " << (a2 - a1) / 1000000.0 / 60 << " min (" << (a2 - a1) / 1000000.0 << " sec)" << endl;
  }
  else if ((a2 - a1) / 1000000.0 < 60 * 60 * 24 * 2) {
    *logfile << "CPU time      : " << int((t2 - t1) / 1000000.0 / (60 * 60)) << " hour " << 60.0*((t2 - t1) / 1000000.0 / (60 * 60) - int((t2 - t1) / 1000000.0 / (60 * 60))) << " min (" << (t2 - t1) / 1000000.0 << " sec)" << endl;
    *logfile << "Physical time : " << int((a2 - a1) / 1000000.0 / (60 * 60)) << " hour " << 60.0*((a2 - a1) / 1000000.0 / (60 * 60) - int((a2 - a1) / 1000000.0 / (60 * 60))) << " min (" << (a2 - a1) / 1000000.0 << " sec)" << endl;
  }
  else
  {
    *logfile << "CPU time      : " << int((t2 - t1) / 1000000.0 / (60 * 60 * 24)) << " days " << 24.0*((t2 - t1) / 1000000.0 / (60 * 60 * 24) - int((t2 - t1) / 1000000.0 / (60 * 60 * 24))) << " hour (" << setprecision(8) << (t2 - t1) / 1000000.0 << " sec)" << endl;
    *logfile << "Physical time : " << int((a2 - a1) / 1000000.0 / (60 * 60 * 24)) << " days " << 24.0*((a2 - a1) / 1000000.0 / (60 * 60 * 24) - int((a2 - a1) / 1000000.0 / (60 * 60 * 24))) << " hour (" << setprecision(8) << (a2 - a1) / 1000000.0 << " sec)" << endl;
  }
  *logfile << "CPU usgae     : " << (t2 - t1)*100.0 / (a2 - a1) << " %  (priority nice = " << getpriority(PRIO_PROCESS, 0) << ")" << endl;
  if (a1 != a1s) {
    *logfile << "Time after sleeping" << endl;

    if ((a2 - a1) / 1000000.0 < 60) {
      *logfile << "CPU time      : " << (t2 - t1s) / 1000000.0 << " sec" << endl;
      *logfile << "Physical time : " << (a2 - a1s) / 1000000.0 << " sec" << endl;
    }
    else if ((a2 - a1) / 1000000.0 < 60 * 60) {
      *logfile << "CPU time      : " << (t2 - t1s) / 1000000.0 / 60 << " min (" << (t2 - t1s) / 1000000.0 << " sec)" << endl;
      *logfile << "Physical time : " << (a2 - a1s) / 1000000.0 / 60 << " min (" << (a2 - a1s) / 1000000.0 << " sec)" << endl;
    }
    else if ((a2 - a1) / 1000000.0 < 60 * 60 * 24 * 2) {
      *logfile << "CPU time      : " << int((t2 - t1s) / 1000000.0 / (60 * 60)) << " hour " << 60.0*((t2 - t1s) / 1000000.0 / (60 * 60) - int((t2 - t1s) / 1000000.0 / (60 * 60))) << " min (" << (t2 - t1s) / 1000000.0 << " sec)" << endl;
      *logfile << "Physical time : " << int((a2 - a1s) / 1000000.0 / (60 * 60)) << " hour " << 60.0*((a2 - a1s) / 1000000.0 / (60 * 60) - int((a2 - a1s) / 1000000.0 / (60 * 60))) << " min (" << (a2 - a1s) / 1000000.0 << " sec)" << endl;
    }
    else {
      *logfile << "CPU time      : " << int((t2 - t1s) / 1000000.0 / (60 * 60 * 24)) << " days " << 24.0*((t2 - t1s) / 1000000.0 / (60 * 60 * 24) - int((t2 - t1s) / 1000000.0 / (60 * 60 * 24))) << " hour (" << setprecision(8) << (t2 - t1s) / 1000000.0 << " sec)" << endl;
      *logfile << "Physical time : " << int((a2 - a1s) / 1000000.0 / (60 * 60 * 24)) << " days " << 24.0*((a2 - a1s) / 1000000.0 / (60 * 60 * 24) - int((a2 - a1s) / 1000000.0 / (60 * 60 * 24))) << " hour (" << setprecision(8) << (a2 - a1s) / 1000000.0 << " sec)" << endl;
    }
    *logfile << "CPU usgae     : " << (t2 - t1s)*100.0 / (a2 - a1s) << " %  (priority nice = " << getpriority(PRIO_PROCESS, 0) << ")" << endl;
  }

#ifdef _OPENMP
  *logfile << "CPU efficiency: " << 100.0*(t2 - t1) / (a2 - a1) / omp_get_max_threads() << " %" << endl;
#endif

  FILE *fp = fopen("/proc/cpuinfo", "r");
  if (fp != NULL) {
    *logfile << endl;
    char tmp[255];
    while (fgets(tmp, 255, fp) != NULL) {
      *logfile << tmp;
    }
    fclose(fp);
    *logfile << endl;
  }


#define MAXLINE 1024
  char s[MAXLINE];
  FILE *fpin;
  fpin = popen("sensors", "r");
  if (fpin != NULL) {
    while (fgets(s, MAXLINE, fpin) != NULL){
      *logfile << s;
    }
    pclose(fpin);
    *logfile << endl;
  }

  fpin = popen("lspci", "r");
  if (fpin != NULL) {
    while (fgets(s, MAXLINE, fpin) != NULL){
      *logfile << s;
    }
    pclose(fpin);
    *logfile << endl;
  }

  fp = fopen("/proc/meminfo", "r");
  if (fp != NULL) {
    char tmp[255];
    while (fgets(tmp, 255, fp) != NULL) {
      *logfile << tmp;
    }
    fclose(fp);
    *logfile << endl;
  }


  fpin = popen("df -h", "r");
  if (fpin != NULL) {
    while (fgets(s, MAXLINE, fpin) != NULL){
      *logfile << s;
    }
    pclose(fpin);
    *logfile << endl;
  }

  fpin = popen("ls -la --full-time", "r");
  if (fpin != NULL) {
    while (fgets(s, MAXLINE, fpin) != NULL){
      *logfile << s;
    }
    pclose(fpin);
    *logfile << endl;
  }
#else
  if ((a2 - a1) / 1000000.0 < 60) {
    *logfile << "Physical time : " << (a2 - a1) / 1000000.0 << " sec" << endl;
  }
  else if ((a2 - a1) / 1000000.0 < 60 * 60) {
    *logfile << "Physical time : " << (a2 - a1) / 1000000.0 / 60 << " min (" << (a2 - a1) / 1000000.0 << " sec)" << endl;
  }
  else if ((a2 - a1) / 1000000.0 < 60 * 60 * 24 * 2) {
    *logfile << "Physical time : " << (a2 - a1) / 1000000.0 / (60 * 60) << " hour (" << (a2 - a1) / 1000000.0 << " sec)" << endl;
  }
  else {
    *logfile << "Physical time : " << (a2 - a1) / 1000000.0 / (60 * 60 * 24) << " days (" << setprecision(8) << (a2 - a1) / 1000000.0 << " sec)" << endl;
  }
#endif

#ifdef __NVCC__
  //GPU info

  int deviceCount = 0;
  cudaError_t error_id = cudaGetDeviceCount(&deviceCount);
  if (error_id != cudaSuccess)return;
  *logfile << endl;
  *logfile << "GPU info" << endl;
  *logfile << "CUDA Capable device(s) : " << deviceCount << endl;

  int dev, driverVersion = 0, runtimeVersion = 0;
  for (dev = 0; dev < deviceCount; ++dev)
  {
    cudaSetDevice(dev);
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, dev);
    *logfile << "Device " << dev << " : " << deviceProp.name << endl;

    cudaDriverGetVersion(&driverVersion);
    cudaRuntimeGetVersion(&runtimeVersion);
    *logfile << "  CUDA Driver Version : " << driverVersion / 1000 << "." << (driverVersion % 100) / 10 << endl;
    *logfile << "  CUDA Runtime Version : " << runtimeVersion / 1000 << "." << (runtimeVersion % 100) / 10 << endl;
    *logfile << "  CUDA Capability : " << deviceProp.major << "." << deviceProp.minor << endl;

    *logfile << "  Total amount of global memory : " << (float)deviceProp.totalGlobalMem / 1048576.0f << " MB" << endl;

    *logfile << "  " << deviceProp.multiProcessorCount << " Multiprocessors, " << _ConvertSMVer2Cores(deviceProp.major, deviceProp.minor) << " CUDA Cores/MP: " << _ConvertSMVer2Cores(deviceProp.major, deviceProp.minor) * deviceProp.multiProcessorCount << "CUDA Cores" << endl;

    *logfile << "  GPU Max Clock rate : " << deviceProp.clockRate * 1e-3f << " MHz" << endl;

    *logfile << "  Memory Clock rate : " << deviceProp.memoryClockRate * 1e-3f << " MHz" << endl;
    *logfile << "  Memory Bus Width : " << deviceProp.memoryBusWidth << "-bit" << endl;
    if (deviceProp.l2CacheSize)
      *logfile << "  L2 Cache Size : " << deviceProp.l2CacheSize << " bytes" << endl;
    *logfile << "  Maximum Texture Dimension Size (x,y,z)" << endl;
    *logfile << "    1D=" << deviceProp.maxTexture1D << endl;
    *logfile << "    2D=" << deviceProp.maxTexture2D[0] << " , " << deviceProp.maxTexture2D[1] << endl;
    *logfile << "    3D=" << deviceProp.maxTexture3D[0] << " , " << deviceProp.maxTexture3D[1] << " , " << deviceProp.maxTexture3D[2] << endl;
    *logfile << "    Maximum Layered 1D Texture Size, (num) layers  1D = (" << deviceProp.maxTexture1DLayered[0] << "), " << deviceProp.maxTexture1DLayered[1] << " layers" << endl;
    *logfile << "    Maximum Layered 2D Texture Size, (num) layers  2D=(" << deviceProp.maxTexture2DLayered[0] << "," << deviceProp.maxTexture2DLayered[1] << "), " << deviceProp.maxTexture2DLayered[2] << " layers" << endl;
    *logfile << "  Total amount of constant memory : " << deviceProp.totalConstMem << " bytes" << endl;
    *logfile << "  Total amount of shared memory per block : " << deviceProp.sharedMemPerBlock << " bytes" << endl;
    *logfile << "  Total number of registers available per block : " << deviceProp.regsPerBlock << endl;
    *logfile << "  Warp size : " << deviceProp.warpSize << endl;
    *logfile << "  Maximum number of threads per multiprocessor : " << deviceProp.maxThreadsPerMultiProcessor << endl;
    *logfile << "  Maximum number of threads per block : " << deviceProp.maxThreadsPerBlock << endl;
    *logfile << "  Max dimension size of a thread block (x,y,z) : (" << deviceProp.maxThreadsDim[0] << ", " << deviceProp.maxThreadsDim[1] << ", " << deviceProp.maxThreadsDim[2] << ")" << endl;
    *logfile << "  Max dimension size of a grid size    (x,y,z): (" << deviceProp.maxGridSize[0] << ", " << deviceProp.maxGridSize[1] << ", " << deviceProp.maxGridSize[2] << ")" << endl;
    *logfile << "  Maximum memory pitch : " << deviceProp.memPitch << " bytes" << endl;
    *logfile << "  Texture alignment : " << deviceProp.textureAlignment << " bytes" << endl;
    *logfile << "  Concurrent copy and kernel execution : " << (deviceProp.deviceOverlap ? "Yes" : "No") << " with " << deviceProp.asyncEngineCount << " copy engine(s)" << endl;
    *logfile << "  Run time limit on kernels : " << (deviceProp.kernelExecTimeoutEnabled ? "Yes" : "No") << endl;
    *logfile << "  Integrated GPU sharing Host Memory : " << (deviceProp.integrated ? "Yes" : "No") << endl;
    *logfile << "  Support host page-locked memory mapping : " << (deviceProp.canMapHostMemory ? "Yes" : "No") << endl;
    *logfile << "  Alignment requirement for Surfaces : " << (deviceProp.surfaceAlignment ? "Yes" : "No") << endl;
    *logfile << "  Device has ECC support : " << (deviceProp.ECCEnabled ? "Enabled" : "Disabled") << endl;
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
    *logfile << "  CUDA Device Driver Mode (TCC or WDDM) : " << (deviceProp.tccDriver ? "TCC (Tesla Compute Cluster Driver)" : "WDDM (Windows Display Driver Model)") << endl;
#endif
    *logfile << "  Device supports Unified Addressing (UVA) : " << (deviceProp.unifiedAddressing ? "Yes" : "No") << endl;
    *logfile << "  Device PCI Domain ID / Bus ID / location ID : " << deviceProp.pciDomainID << " / " << deviceProp.pciBusID << " / " << deviceProp.pciDeviceID << endl;

    const char *sComputeMode[] =
    {
      "Default (multiple host threads can use ::cudaSetDevice() with device simultaneously)",
      "Exclusive (only one host thread in one process is able to use ::cudaSetDevice() with this device)",
      "Prohibited (no host thread can use ::cudaSetDevice() with this device)",
      "Exclusive Process (many threads in one process is able to use ::cudaSetDevice() with this device)",
      "Unknown",
      NULL
    };
    *logfile << "  Compute Mode : " << sComputeMode[deviceProp.computeMode] << endl;
  }
  // If there are 2 or more GPUs, query to determine whether RDMA is supported
  if (deviceCount >= 2)
  {
    cudaDeviceProp prop[64];
    int gpuid[64]; // we want to find the first two GPUs that can support P2P
    int gpu_p2p_count = 0;

    for (int i = 0; i < deviceCount; i++)
    {
      cudaGetDeviceProperties(&prop[i], i);

      // Only boards based on Fermi or later can support P2P
      if ((prop[i].major >= 2)
#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
        // on Windows (64-bit), the Tesla Compute Cluster driver for windows must be enabled to support this
        && prop[i].tccDriver
#endif
        )
      {
        // This is an array of P2P capable GPUs
        gpuid[gpu_p2p_count++] = i;
      }
    }
    // Show all the combinations of support P2P GPUs
    int can_access_peer;
    if (gpu_p2p_count >= 2)
    {
      for (int i = 0; i < gpu_p2p_count; i++)
      {
        for (int j = 0; j < gpu_p2p_count; j++)
        {
          if (gpuid[i] == gpuid[j])
          {
            continue;
          }
          cudaDeviceCanAccessPeer(&can_access_peer, gpuid[i], gpuid[j]);
          *logfile << "Peer access from " << prop[gpuid[i]].name << " (GPU" << gpuid[i] << ") -> " << prop[gpuid[j]].name << " (GPU" << gpuid[j] << ") : " << (can_access_peer ? "Yes" : "No") << endl;
        }
      }
    }
  }
  cudaDeviceReset();
#endif
#endif
}
long long int Logging::GetTime(char *tim) {
  time_t now;
  struct tm *tptr;
  char *weekday[] = { (char *)"Sun", (char *)"Mon", (char *)"Tue", (char *)"Wed", (char *)"Thr", (char *)"Fri", (char *)"Sat" };
  now = time(NULL);
  tptr = localtime(&now);
#ifndef _WIN32
  struct timeval tv;
  struct timezone tz;
  gettimeofday(&tv, &tz);

  time_t timer;
  struct tm *timerlocal;
  timer = time(NULL);
  timerlocal = localtime(&timer);
  strftime(tim, 255, "%Y/%m/%d (%a)  %H:%M:%S %z (%Z)", timerlocal);
  return (long long int)((long long int)tv.tv_sec * 1000000 + tv.tv_usec);
#else
  TIME_ZONE_INFORMATION tzi;
  DWORD res = GetTimeZoneInformation(&tzi);
  char tzname[255];
  if (res == TIME_ZONE_ID_STANDARD) {
    wcstombs(tzname, tzi.StandardName, 255);
  }
  else if (res == TIME_ZONE_ID_DAYLIGHT) {
    wcstombs(tzname, tzi.DaylightName, 255);
  }
  else {
    tzname[0] = 0;
  }
  sprintf(tim, "%4d/%02d/%02d (%s)  %02d:%02d:%02d %2.2d00 (%s)",
    (tptr->tm_year) + 1900, tptr->tm_mon + 1, tptr->tm_mday,
    weekday[tptr->tm_wday], tptr->tm_hour, tptr->tm_min, tptr->tm_sec, (-tzi.Bias - tzi.DaylightBias) / 60, tzname);
  return (long long int)timeGetTime() * 1000;
#endif
}
void Logging::Note(string note) {
  if (node == 0) {
    if (firstnote) {
      *logfile << "\n";
    }
    char ttt[256];
    n2 = GetTime(ttt);
    *logfile << "Note : " << ttt;

    if ((n2 - n1) / 1000000.0 < 60) {
      *logfile << "  [Physical time interval : " << (n2 - n1) / 1000000.0 << " sec]" << endl;
    }
    else if ((n2 - n1) / 1000000.0 < 60 * 60) {
      *logfile << "  [Physical time interval : " << (n2 - n1) / 1000000.0 / 60 << " min (" << (n2 - n1) / 1000000.0 << " sec)]" << endl;
    }
    else if ((n2 - n1) / 1000000.0 < 60 * 60 * 24 * 2) {
      *logfile << "  [Physical time interval : " << int((n2 - n1) / 1000000.0 / (60 * 60)) << " hour " << 60.0*((n2 - n1) / 1000000.0 / (60 * 60) - int((n2 - n1) / 1000000.0 / (60 * 60))) << " min (" << (n2 - n1) / 1000000.0 << " sec)]" << endl;
    }
    else {
      *logfile << "  [Physical time interval : " << int((n2 - n1) / 1000000.0 / (60 * 60 * 24)) << " days " << int(24.0*((n2 - n1) / 1000000.0 / (60 * 60 * 24) - int((n2 - n1) / 1000000.0 / (60 * 60 * 24)))) << " hour " << 60.0*((n2 - n1) / 1000000.0 / (60 * 60) - int((n2 - n1) / 1000000.0 / (60 * 60))) << " min (" << setprecision(8) << (n2 - n1) / 1000000.0 << " sec)]" << endl;
    }
    n1 = n2;
    *logfile << note << "\n" << endl;
  }
}
long long int Logging::Cputime() {
#ifndef _WIN32
  struct rusage RU;
  getrusage(RUSAGE_SELF, &RU);
  return (long long int)((long long int)(RU.ru_utime.tv_sec) * 1000000.0 + (double)RU.ru_utime.tv_usec);
#else
  return 0;
#endif
}


void Logging::WriteHostInfo(void) {
  //
  // write hostname, OS version etc
  //
#ifndef _WIN32
  struct utsname uts;
  uname(&uts);
  *logfile << "Host name : " << uts.nodename << endl;
  *logfile << "OS name and version : " << uts.sysname << "   " << uts.release << "  " << uts.version << endl;
  FILE *fp = fopen("/etc/lsb-release", "r");
  if (fp != NULL) {
    char tmp[255];
    while (fgets(tmp, 255, fp) != NULL) {
      *logfile << tmp;
    }
    fclose(fp);
  }
  FILE *fp2 = fopen("/etc/redhat-release", "r");// for CentOS
  if (fp2 != NULL) {
    char tmp[255];
    while (fgets(tmp, 255, fp2) != NULL) {
      *logfile << tmp;
    }
    fclose(fp2);
  }

#endif

#ifdef _WIN32
  char acComputerName[32767];
  DWORD nComputerName = sizeof(acComputerName);
  if (GetComputerName(acComputerName, &nComputerName)) {
    *logfile << "Computer Name : " << acComputerName << endl;
  }
  else {
    cerr << "Failed to lookup computer name, error code " <<
      GetLastError() << "." << endl;
  }
  OSVERSIONINFO osvi;
  ZeroMemory(&osvi, sizeof(OSVERSIONINFO));
  osvi.dwOSVersionInfoSize = sizeof(OSVERSIONINFO);
  GetVersionEx(&osvi);
  *logfile << "OS Version number : " << osvi.dwMajorVersion << "." << osvi.dwMinorVersion << " build " << osvi.dwBuildNumber << " " << osvi.szCSDVersion << endl;


  {
    typedef void (WINAPI *PGNSI)(LPSYSTEM_INFO);
    typedef BOOL(WINAPI *PGPI)(DWORD, DWORD, DWORD, DWORD, PDWORD);

    OSVERSIONINFOEX osvi;
    SYSTEM_INFO si;
    PGNSI pGNSI;
    PGPI pGPI;
    BOOL bOsVersionInfoEx;
    DWORD dwType;

    ZeroMemory(&si, sizeof(SYSTEM_INFO));
    ZeroMemory(&osvi, sizeof(OSVERSIONINFOEX));

    osvi.dwOSVersionInfoSize = sizeof(OSVERSIONINFOEX);
    bOsVersionInfoEx = GetVersionEx((OSVERSIONINFO*)&osvi);

    if (bOsVersionInfoEx == 0) return;

    // Call GetNativeSystemInfo if supported or GetSystemInfo otherwise.

    pGNSI = (PGNSI)GetProcAddress(
      GetModuleHandle(TEXT("kernel32.dll")),
      "GetNativeSystemInfo");
    if (NULL != pGNSI)
      pGNSI(&si);
    else GetSystemInfo(&si);

    if (VER_PLATFORM_WIN32_NT == osvi.dwPlatformId &&
      osvi.dwMajorVersion > 4) {
      *logfile << "Microsoft ";

      // Test for the specific product.

      if (osvi.dwMajorVersion == 6) {
        if (osvi.dwMinorVersion == 0) {
          if (osvi.wProductType == VER_NT_WORKSTATION)
            *logfile << "Windows Vista ";
          else *logfile << "Windows Server 2008 ";
        }

        if (osvi.dwMinorVersion == 1) {
          if (osvi.wProductType == VER_NT_WORKSTATION)
            *logfile << "Windows 7 ";
          else *logfile << "Windows Server 2008 R2 ";
        }

        if (osvi.dwMinorVersion == 2) {
          if (osvi.wProductType == VER_NT_WORKSTATION)
            *logfile << "Windows 8 ";
          else *logfile << "Windows Server 2012 ";
        }

        if (osvi.dwMinorVersion == 3) {
          if (osvi.wProductType == VER_NT_WORKSTATION)
            *logfile << "Windows 8.1 ";
          else *logfile << "Windows Server 2012 R2 ";
        }

        pGPI = (PGPI)GetProcAddress(
          GetModuleHandle(TEXT("kernel32.dll")),
          "GetProductInfo");

        pGPI(osvi.dwMajorVersion, osvi.dwMinorVersion, 0, 0, &dwType);


        //
#ifndef PRODUCT_PROFESSIONAL
#define PRODUCT_PROFESSIONAL 0x00000030
#endif
        //
        switch (dwType)
        {
        case PRODUCT_ULTIMATE:
          *logfile << "Ultimate Edition";
          break;
        case PRODUCT_PROFESSIONAL:
          *logfile << "Professional";
          break;
        case PRODUCT_HOME_PREMIUM:
          *logfile << "Home Premium Edition";
          break;
        case PRODUCT_HOME_BASIC:
          *logfile << "Home Basic Edition";
          break;
        case PRODUCT_ENTERPRISE:
          *logfile << "Enterprise Edition";
          break;
        case PRODUCT_BUSINESS:
          *logfile << "Business Edition";
          break;
        case PRODUCT_STARTER:
          *logfile << "Starter Edition";
          break;
        case PRODUCT_CLUSTER_SERVER:
          *logfile << "Cluster Server Edition";
          break;
        case PRODUCT_DATACENTER_SERVER:
          *logfile << "Datacenter Edition";
          break;
        case PRODUCT_DATACENTER_SERVER_CORE:
          *logfile << "Datacenter Edition (core installation)";
          break;
        case PRODUCT_ENTERPRISE_SERVER:
          *logfile << "Enterprise Edition";
          break;
        case PRODUCT_ENTERPRISE_SERVER_CORE:
          *logfile << "Enterprise Edition (core installation)";
          break;
        case PRODUCT_ENTERPRISE_SERVER_IA64:
          *logfile << "Enterprise Edition for Itanium-based Systems";
          break;
        case PRODUCT_SMALLBUSINESS_SERVER:
          *logfile << "Small Business Server";
          break;
        case PRODUCT_SMALLBUSINESS_SERVER_PREMIUM:
          *logfile << "Small Business Server Premium Edition";
          break;
        case PRODUCT_STANDARD_SERVER:
          *logfile << "Standard Edition";
          break;
        case PRODUCT_STANDARD_SERVER_CORE:
          *logfile << "Standard Edition (core installation)";
          break;
        case PRODUCT_WEB_SERVER:
          *logfile << "Web Server Edition";
          break;
        }
      }

      if (osvi.dwMajorVersion == 5 && osvi.dwMinorVersion == 2) {
        if (GetSystemMetrics(SM_SERVERR2))
          *logfile << "Windows Server 2003 R2, ";
        else if (osvi.wSuiteMask & VER_SUITE_STORAGE_SERVER)
          *logfile << "Windows Storage Server 2003";
        else if (osvi.wSuiteMask & VER_SUITE_WH_SERVER)
          *logfile << "Windows Home Server";
        else if (osvi.wProductType == VER_NT_WORKSTATION &&
          si.wProcessorArchitecture == PROCESSOR_ARCHITECTURE_AMD64) {
          *logfile << "Windows XP Professional x64 Edition";
        }
        else *logfile << "Windows Server 2003, ";

        // Test for the server type.
        if (osvi.wProductType != VER_NT_WORKSTATION) {
          if (si.wProcessorArchitecture == PROCESSOR_ARCHITECTURE_IA64) {
            if (osvi.wSuiteMask & VER_SUITE_DATACENTER)
              *logfile << "Datacenter Edition for Itanium-based Systems";
            else if (osvi.wSuiteMask & VER_SUITE_ENTERPRISE)
              *logfile << "Enterprise Edition for Itanium-based Systems";
          }
          else if (si.wProcessorArchitecture == PROCESSOR_ARCHITECTURE_AMD64) {
            if (osvi.wSuiteMask & VER_SUITE_DATACENTER)
              *logfile << "Datacenter x64 Edition";
            else if (osvi.wSuiteMask & VER_SUITE_ENTERPRISE)
              *logfile << "Enterprise x64 Edition";
            else *logfile << "Standard x64 Edition";
          }
          else {
            if (osvi.wSuiteMask & VER_SUITE_COMPUTE_SERVER)
              *logfile << "Compute Cluster Edition";
            else if (osvi.wSuiteMask & VER_SUITE_DATACENTER)
              *logfile << "Datacenter Edition";
            else if (osvi.wSuiteMask & VER_SUITE_ENTERPRISE)
              *logfile << "Enterprise Edition";
            else if (osvi.wSuiteMask & VER_SUITE_BLADE)
              *logfile << "Web Edition";
            else *logfile << "Standard Edition";
          }
        }
      }

      if (osvi.dwMajorVersion == 5 && osvi.dwMinorVersion == 1) {
        *logfile << "Windows XP ";
        if (osvi.wSuiteMask & VER_SUITE_PERSONAL)
          *logfile << "Home Edition";
        else *logfile << "Professional";
      }

      if (osvi.dwMajorVersion == 5 && osvi.dwMinorVersion == 0) {
        *logfile << "Windows 2000 ";

        if (osvi.wProductType == VER_NT_WORKSTATION) {
          *logfile << "Professional";
        }
        else {
          if (osvi.wSuiteMask & VER_SUITE_DATACENTER)
            *logfile << "Datacenter Server";
          else if (osvi.wSuiteMask & VER_SUITE_ENTERPRISE)
            *logfile << "Advanced Server";
          else *logfile << "Server";
        }
      }

      // Include service pack (if any) and build number.

      if (osvi.szCSDVersion > 0) {
        *logfile << " ";
        *logfile << osvi.szCSDVersion;
      }

      *logfile << " (build " << osvi.dwBuildNumber << ")";

      if (osvi.dwMajorVersion >= 6) {
        if (si.wProcessorArchitecture == PROCESSOR_ARCHITECTURE_AMD64)
          *logfile << ", 64-bit";
        else if (si.wProcessorArchitecture == PROCESSOR_ARCHITECTURE_INTEL)
          *logfile << ", 32-bit";
      }
      *logfile << endl;
    }
  }
#endif
}

void Logging::WriteProcessInfo(void)
{
  //
  // write PID, username, dir etc
  //

#ifndef _WIN32
  struct passwd *pwdptr;
  pwdptr = getpwuid(getuid());
  *logfile << "User : " << pwdptr->pw_name << " (UID=" << pwdptr->pw_uid << ", GID=" << pwdptr->pw_gid << ", HOME=" << pwdptr->pw_dir << ")" << endl;
  pid_t pid = getpid();
  *logfile << "PID : " << pid << endl;
  char dir[255];
  if (getcwd(dir, 255))
    *logfile << "Execute Directory : " << dir << endl;
#endif

#ifdef _WIN32
  char acUserName[32767];
  DWORD nUserName = sizeof(acUserName);
  if (GetUserName(acUserName, &nUserName)) { //need advapi32.lib
    *logfile << "User : " << acUserName << endl;
  }
  else {
    cerr << "Failed to lookup user name, error code " <<
      GetLastError() << "." << endl;
  }
  if (DWORD PID = GetCurrentProcessId()) {
    *logfile << "Process ID : " << PID << endl;
  }
  char dir[255];
  GetCurrentDirectory(255, dir);
  *logfile << "Execute Directory : " << dir << endl;
#endif
}

void Logging::WriteCompilerInfo(void)
{
#ifdef __INTEL_COMPILER
  *logfile << "Intel Compiler " << __INTEL_COMPILER / 100.0 << " build " << __INTEL_COMPILER_BUILD_DATE << endl;
#endif
#ifdef _MSC_VER
  *logfile << "Microsoft Visual C++ " << _MSC_VER;
#ifdef _MSC_BUILD
  *logfile << " build " << _MSC_BUILD;
#endif
  *logfile << endl;
#endif
#ifdef __GNUC__
  *logfile << "GNU C/C++ " << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__ << endl;
#endif
#ifdef __BORLANDC__
  *logfile << "Borland C++ " << __BORLANDC__ << endl;
#endif
#ifdef __NVCC__
  *logfile << "nvcc : yes" << endl;
#endif

#ifdef _OPENMP
  *logfile << "OpenMP (" << _OPENMP << ") : Enabled (Max # of threads = " << omp_get_max_threads() << ")" << endl;
#endif
}

void Logging::MemoryUsage(void) {
  struct utsname uts;
  uname(&uts);
  pid_t pid = getpid();
  char tmp[255];
  sprintf(tmp, "cp -p /proc/%d/status ./memusg_%s_%05d.log", pid, uts.nodename, pid);
  int ret = system(tmp);
  if (ret)cout << "system function error" << endl;
}
