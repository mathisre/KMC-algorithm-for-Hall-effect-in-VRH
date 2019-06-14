#include "params.h"
#include "paramfile.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <iostream>

#pragma hdrstop


#include "systemclass.h"
#include "CGanalysis.h"
#include "fileutils.h"
#include "stringutils.h"
#include "paramfile.h"
#include "MKcurrent.h"
#include "params.h"
#include "run.h"

void params::readparams(int n, char* cmdline[])
{
  int a;
  ctype = 5; // 1 normal, 2 heat map, 3 from state file, 4 from state spes file, 5 rejection algo, 6 rejection maps, 7 rejection random lattice, 8 dynamic random lattice
  cutoffexp = 7;
  d2el = 5;
  d = 1;
  dim = 2;
  disU = 1.0; // Magnitude of disorder potential
  doTriJumps = true;
  eonsteps = 256; //# of time step in which the ext field is on
  extf = 0;
  Ex = 0.040;
  Ey = 0.0;
  Ez = 0.0;
  fill = 0.5;
  fillPassive = 0.5;
  htype = 3;    //1: NN, 2: random,
  Hx = 0.0;
  Hy = 0.0;
  Hz = 0.2;
  initruns = 1;
  inittemp = 0.02;
  istate = 1; //random state
  lattdis = 0.0; //lattice distortion
  len = 100;
  limitOfNoReturn = 100.0;
  loc = 3.0; // loc = 2/a
  maxExpMax = 20;
  maxExpMin = 3;
  maxJL2el = 5;
  maxJL = 5;
  maxProbability = 20;
  omega  =  0;
  ratefun = 0;
  rhopps = 100;
  rsconfig = 30;    //the disorder config
  rsstate = 2;     //the initial state config
  running = false;
  runs = 1;
  screen = 0.5; //screening length: < 0 sets: len/2 (max: len-1)
  seed2 = 1;
  skipsteps = 0; //all timesteps are recorded
  spesfilename = "";
  spesfilenamePassive = "";
  startrun = 0;
  starttime = 0;
  statefilename = "";
  statefilenamePassive = "";
  stopfilename = "stop";
  stopfilename = "stop";
  temp = 0.40;
  tposn = 0;
  trace = false;
  tracenum = 0; // = all particles
  xpbc = true; //periodic boundary condition in x-direction can be lifted
  outputprefix = "/home/mathisre/Master/data/";
  randShift = 0.0;


  char Hz_s [10], Ex_s[10], temp_s[10];
  sprintf(Hz_s, "%.5f", Hz);
  sprintf(Ex_s, "%.5f", Ex);
  sprintf(temp_s, "%.5f", temp);

  outputendfix = "Hz_" + string(Hz_s)+  "_Ex_" + string(Ex_s) + "_T_" + string(temp_s);
  writelines = 10000;
  timesteps = pow(10,4);
  if(n >1)
  {

    paramfilereader *pf = new paramfilereader();
    a = pf->cmdlineopenread(n,cmdline);
    if(a < 1) a = pf->openread("esystem.ini");
    if(a > 0)
    {  //traceing settings
       if(outputprefix == "") outputprefix = "out";
       if(pf->getstring("trace") == "yes") trace = true;
       if(pf->getstring("xpbc") == "no") xpbc = false;
       if(stopfilename == "") stopfilename = "stop";
       doTriJumps = bool(pf->getint("triJumps", doTriJumps));
       Ex = pf->getdouble("Ex",Ex);
       Ey = pf->getdouble("Ey",Ey);
       Ez = pf->getdouble("Ez",Ez);
       Hx = pf->getdouble("Hx",Hx);
       Hy = pf->getdouble("Hy",Hy);
       Hz = pf->getdouble("Hz",Hz);
       ctype = pf->getint("ctype",ctype);
       cutoffexp = pf->getdouble("cutoffexp",cutoffexp);
       d = pf->getdouble("d",d);
       d2el = pf->getint("d2el",d2el);
       dim = pf->getint("dim",dim);
       disU = pf->getdouble("U",disU);
       eonsteps = pf->getint("NEon",eonsteps);
       extf = pf->getdouble("E",extf);
       fill = pf->getdouble("fill",fill);
       fillPassive = fill;
       fillPassive = pf->getdouble("fillPassive",fillPassive);
       htype = pf->getint("htype",htype);
       initruns = pf->getint("initruns",initruns);
       inittemp = pf->getdouble("initT",inittemp);
       istate = pf->getint("istate",istate);
       lattdis = pf->getdouble("lattdis",lattdis);
       len = pf->getint("L",len);
       limitOfNoReturn = pf->getdouble("limitOfNoReturn",limitOfNoReturn);
       loc = pf->getdouble("loc",loc);
       maxExpMax = pf->getint("maxExpMax", maxExpMax);
       maxExpMin = pf->getint("maxExpMin", maxExpMin);
       maxJL = pf->getint("maxJL",maxJL);
       maxJL2el = pf->getint("maxJL2el",maxJL2el);
       maxProbability = pf->getdouble("maxProbability",maxProbability);
       omega = pf->getdouble("omega",omega);
       outputprefix = pf->getstring("outpre");
       ratefun = pf->getint("ratefun",ratefun);
       rhopps = pf->getint("rhopps",rhopps);
       rsconfig = pf->getint("rsconfig",rsconfig);
       rsstate = pf->getint("rsstate",rsstate);
       runs = pf->getint("runs",runs);
       screen = pf->getdouble("screen",screen);
       seed2 = pf->getint("seed2",seed2);
       skipsteps = pf->getint("skipsteps",skipsteps);
       spesfilename = pf->getstring("spesfilename");
       spesfilenamePassive = pf->getstring("spesfilenamePassive");
       startrun = pf->getint("startrun",startrun);
       starttime = pf->getdouble("starttime",starttime);
       statefilename = pf->getstring("statefilename");
       statefilenamePassive = pf->getstring("statefilenamePassive");
       stopfilename = pf->getstring("stopfn");
       temp = pf->getdouble("T",temp);
       timesteps = pf->getint("Nt",timesteps);
       tposn = pf->getintarray("tpos",tpos,maxtpos);
       tracenum = pf->getint("tracenum",tracenum);
       writelines = pf->getint("writelines",writelines);
       randShift = pf->getdouble("randShift",randShift);
   }

    char Hz_s [10], Ex_s[10], temp_s[10];
    sprintf(Hz_s, "%.5f", Hz);
    sprintf(Ex_s, "%.5f", Ex);
    sprintf(temp_s, "%.5f", temp);


    if (doTriJumps==true) outputendfix = "Hz_" + string(Hz_s)+  "_Ex_" + string(Ex_s) + "_T_" + string(temp_s);
    else outputendfix = "noTri_Hz_" + string(Hz_s)+  "_Ex_" + string(Ex_s) + "_T_" + string(temp_s);
   delete pf;
  }
}
