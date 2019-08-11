//---------------------------------------------------------------------------

#include "fileutils.h"
#include "stringutils.h"
#include "CGanalysis.h"
#include "systemclass.h"
#include "treeutils.h"
#include "MKcurrentRandPos.h"



// Remove eventually
#include <iostream>
#include <fstream>

MKcurrentRandPos::MKcurrentRandPos(int steps)
{
  from = new int[steps];
  to = new int[steps];
  energy = new double[steps];
  MCsteps = new int[steps];
  dxs = new double[steps];
  dys = new double[steps];
  ts = new double[steps];
  des = new double[steps];
}

MKcurrentRandPos::~MKcurrentRandPos()
{
    if (from!=NULL) delete[] from;
    if (to!=NULL) delete[] to;
    if (energy!=NULL) delete[] energy;
    if (MCsteps!=NULL) delete[] MCsteps;
    if (dxs!=NULL) delete[] dxs;
    if (dys!=NULL) delete[] dys;
    if (ts!=NULL) delete[] ts;
    if (des!=NULL) delete[] des;

    for (int i=0; i<N; i++){
        if (GammaT[i]!=NULL) delete[] GammaT[i];
        if (sitePositions[i]!=NULL) delete[] sitePositions[i];
        if (distanceMatrix[i]!=NULL) delete[] distanceMatrix[i];
    }
    if (GammaT!=NULL) delete[] GammaT;
    if (sitePositions!=NULL) delete[] sitePositions;
    if (distanceMatrix!=NULL) delete[] distanceMatrix;

    if (Nmem!=NULL) delete[] Nmem;
    if (GammaTtot!=NULL) delete[] GammaTtot;
    if (GammaTtotForN!=NULL) delete[] GammaTtotForN;
}


void MKcurrentRandPos::setE(double Ex1, double Ey1, double Ez1)
{
  Ex=Ex1; Ey=Ey1; Ez=Ez1;
}

void MKcurrentRandPos::setH(double Hx1, double Hy1, double Hz1)
{
  Hx=Hx1;
  Hy=Hy1;
  Hz=Hz1;
}

void MKcurrentRandPos::setE0(double Ex1, double Ey1, double Ez1)
{
  E0x=Ex1; E0y=Ey1; E0z=Ez1;
}

void MKcurrentRandPos::setOmega(double omega1)
{
  omega=omega1;
}

void MKcurrentRandPos::setWritelines(int writelines1)
{
  writelines=writelines1;
}

void MKcurrentRandPos::setRatefun(int ratefun1)
{
  ratefun = ratefun1;

}

void MKcurrentRandPos::setES(ESystem *e)
{
  es = e;
}


void MKcurrentRandPos::init(int D1, int L1, int N1, int maxJL1, double a1, double beta1, double maxProbability1, bool doTriJumps, double shift)
{ //  (D, L, N, p->maxJL, p->loc, 1.0/p->temp,p->maxProbability)

  int  n,  i, j;
  double dx, dy;
  double gamma;
  double overlap_ij;
  double t0, tau1;

  RanGen = new CRandomMersenne(1);


  D=D1; L=L1; N=N1; A=a1; beta = beta1;
  maxJL = maxJL1; JL2 = 2*maxJL+1;

  maxProbability=maxProbability1/beta;
  randShift = shift;
  t0   = 1;
  tau1 = 1;


//  Nmem = JL2*JL2; // Number of sites we can jump to
//  jl            = new double[Nmem];

  sitePositions   = new double*[N];
  GammaT          = new double*[N];
  GammaTtot       = new double [N];
  GammaTtotForN       = new double [N];

  finalSite    = new int*[N];
  jumpArea = new double*[N];
  Nmem = new int[N];


  dxToSite = new double*[N];
  dyToSite = new double*[N];


  meanTotalGamma = 0;


  // Initializing site positions
  for (i=0; i<N; i++){
    sitePositions[i] = new double[2];
    //sitePositions[i][0] = i%L + (es->ran2(0)-0.5)*randShift;
    //sitePositions[i][1] = i/L + (es->ran2(0)-0.5)*randShift;

    sitePositions[i][0] = es->ran2(0)*L;
    sitePositions[i][1] = es->ran2(0)*L;
  }


  // Setting up arrays to be used
  distanceMatrix = new double*[N];
  for (i=0; i<N; i++){
    distanceMatrix[i] = new double[N];
    Nmem[i] = 0;
    for (j=0; j<N; j++){
        dx = sitePositions[j][0] - sitePositions[i][0];
        dy = sitePositions[j][1] - sitePositions[i][1];

        if      (dx > L/2) dx -= L;
        else if (dx < -L/2) dx += L;

        if      (dy > L/2) dy -= L;
        else if (dy < -L/2) dy += L;


        distanceMatrix[i][j] = sqrt(dx*dx + dy*dy);
        if (distanceMatrix[i][j] < maxJL) Nmem[i]++;
    }
    GammaT[i]  = new double[Nmem[i]];
    finalSite[i] = new int[Nmem[i]];
    dxToSite[i] = new double[Nmem[i]];
    dyToSite[i] = new double[Nmem[i]];


  }


  totGammaTtot  = 0;


  meandx = 0;
  meandy = 0;

  // Calculate all transition rates
  for (i=0; i<N; i++){

      gamma  = 0.0;
      n  = 0;

      for (j=0; j<N; j++){
          overlap_ij = distanceMatrix[i][j];
          if (overlap_ij < maxJL){
              if (i != j){
                  dxToSite[i][n] = sitePositions[j][0] - sitePositions[i][0];
                  dyToSite[i][n] = sitePositions[j][1] - sitePositions[i][1];

                  if      (dxToSite[i][n] >  L/2) dxToSite[i][n] -= L;
                  else if (dxToSite[i][n] < -L/2) dxToSite[i][n] += L;

                  if      (dyToSite[i][n] >  L/2) dyToSite[i][n] -= L;
                  else if (dyToSite[i][n] < -L/2) dyToSite[i][n] += L;

                  gamma += exp(-A*overlap_ij/tau1);
                  GammaT[i][n] = gamma;
                  finalSite[i][n] = j;

                  n++;


              }
              else {
                  GammaT[i][n] = gamma;
                  n++;
              }
          }

      }
      GammaTtot[i] = gamma;
      totGammaTtot += gamma;
      GammaTtotForN[i] = totGammaTtot ;
  }
  meanTotalGamma = totGammaTtot/ N;

  printf("Mean site rate: %.5f\n",totGammaTtot /N);




}

void MKcurrentRandPos::setMTseed(int seed)
{
  RanGen->RandomInit(seed);
}


void inline MKcurrentRandPos::getjump(int &i, int &j, double &dx, double &dy)
{
  double r2;
  int h,l,step, occu, tries;


  tries = 0;
  occu = 0;

  while (occu == 0){
      r2 = RanGen->Random()*totGammaTtot;
      l = -1;
      h = N-1;
      step = N/2;
      while (step > 0)  {
        if (GammaTtotForN[l+step] >= r2) {
            h = l+step;
        }
        else l = l+step;
        step = (h-l)/2;
      }
      tries++;
      occu = es->getocci(h);
  }
  i = h;


  r2 = RanGen->Random()*GammaTtot[i];
  l = -1;
  h = Nmem[i]-1;
  step = Nmem[i]/2;

  while (step > 0)
  {

    if (GammaT[i][l+step] >= r2) {
        h = l+step;
    }
    else {
        l = l+step;
    }
    step = (h-l)/2;
  }
  //at end h=correct jump;

  j = finalSite[i][h];

  dx = dxToSite[i][h];
  dy = dyToSite[i][h];

//  printf("gj: %d %d %d %d h: %d r2: %le\n ",i,j,dx,dy,h,r2);
}


bool inline MKcurrentRandPos::testjump(int i, int j, double dx, double dy)
{
  double dE1, dE2, rate, r2, exponent;

   if (es->getocci(j)!=0){
      return false;
    }

//  printf("dE: %f\n ", dE1);
   dE1 = es->hoppEdiffij(i,j);
  dE2 = dx*Ex + dy*Ey + dE1;
//  printf("dx: %f\n", dx);

  exponent = beta*dE2;

  if(ratefun == 0){   // Approximate rate
    if (dE2 < 0) rate = 1;
    else rate = exp(-exponent);
    r2 = es->ran2(0);
//    printf("r2: %f, rate %f", r2, rate);
    return (r2 < rate);
  }

  if(ratefun == 1){   // Approximate rate
    if (dE2 < 0) rate = 1;
    else rate = exp(-exponent);
    r2 = es->ran2(0);
    printf("rate: %.4f\n",rate);
    return (r2 < rate);
  }
  else if(ratefun == 2){  //Rate according to ES:
   if (dE2 < 0) rate = -dE2*(1+1/(exp(-exponent)-1));
   else if (dE2 == 0) rate = 1/beta;
   else rate = dE2/(exp(exponent)-1);
   if (rate > maxProbability) printf("Rate larger than max rate: %f\n",rate); 

   r2 = es->ran2(0);
   return (maxProbability*r2<rate);

  }
  else if(ratefun == 3){ // Rate according to ES with quadratic prefactor:
    if (dE2<0) rate = dE2*dE2*(1+1/(exp(-exponent)-1)); 
    else if (dE2==0) rate = 0;
    else rate= dE2*dE2/(exp(exponent)-1);
    if (rate > maxProbability) printf("Rate larger than max rate: %f\n",rate); 

    r2 = es->ran2(0);
    return (maxProbability*r2<rate);
  }
}




void MKcurrentRandPos::runCurrent(int steps,double &E, double &t)
{

  int s, MCs, i=0, j=0;
  double dE,dt, dx=0, dy=0, dxI=0, dyI=0;
  bool jump;

  tMC = 1/(               N * es->nu() * meanTotalGamma);

  meanArea = 0;

  meanJumpLength = 0;
  meandx = 0;
  meandy = 0;
  double actualmeanMCs = 0;

  for (s = 0; s < steps; s++)
  {
      MCs = 0;
      jump = false;
      while (!jump)
      {
          MCs++;
          getjump(i,j,dx,dy);
          jump=testjump(i,j,dx,dy);
      }

      dE = es->hopp(i,-1,0,j,0,0);

      dt = MCs*tMC;
      E += dE;
      t += dt;

      actualmeanMCs += MCs;

      meanJumpLength += distanceMatrix[i][j];
      meandx += dx;
      meandy += dy;
      meandxI += dxI;
      meandyI += dyI;



      from[s] = i;
      to[s]   = j;

      ts[s] = t;
      energy[s] = E;
      des[s] = dE;
      dxs[s] = dx;
      dys[s] = dy;

      MCsteps[s]=MCs;

      if (s%100000==0) {
          printf("%6d E=%le MCs=%.3f\n",s,E,actualmeanMCs/100000);
          actualmeanMCs =0;
      }


    }
  printf("\n-------Finishing timesteps--------\n\n");


  meanJumpLength /= steps;
  meandx /= steps;
  meandy /= steps;


  printf("\nMean jump length:              %f\n", meanJumpLength);
  printf("Mean dx:      %8.5f,  Mean dy:      %8.5f\n", meandx, meandy);

}

void MKcurrentRandPos::closePositions(FILE* f)
{
    FileClose(f);
}



void MKcurrentRandPos::jumpsToFileSmall(string filename, int steps)
{
  int s;
  double cumdx, cumdy;
  FILE* f;
  string st;

  cumdx=0;
  cumdy=0;
  f=FileCreate(filename);

  st="t  \t \t E \t \t from \t to \t dx \t dy \t MCs \t dE \t tMC="+DoubleToStr(tMC)+ "\t <jl>=" + DoubleToStr(meanJumpLength) +"\t <dx>=" + DoubleToStr(meandx) +"\t <dy>=" + DoubleToStr(meandy) + "\t <jl_I>=" +DoubleToStr(meanInterJumpLength) +"\t <dx_I>=" + DoubleToStr(meandxI) +"\t <dy_I>=" + DoubleToStr(meandyI) + "\t Acc 2=" +  "\t Acc 3=" +  "\n";
  FileWrite(f,st.c_str(),st.length());

  for(s = 0; s < steps; s++){
      cumdx += dxs[s];
      cumdy += dys[s];
      if (s%writelines == 0){
        st = FloatToStr(ts[s])+"\t"+DoubleToStr(energy[s])+"\t"+IntToStr(from[s])+"\t"+IntToStr(to[s])+"\t"+DoubleToStr(cumdx)+"\t"+DoubleToStr(cumdy)+"\t"+IntToStr(MCsteps[s])+"\t"+FloatToStr(des[s])+"\n";
        FileWrite(f,st.c_str(),st.length());
      }
    }

  FileClose(f);
}

void MKcurrentRandPos::jumpsToFileAC(string filename, int steps)
{
  int s, x1, x2, y1, y2, dx, dy, cumdx;
  FILE* f;
  double cumc, r, dt, t, e, de;
  string st;

  cumdx=0;
  f=FileCreate(filename);

  st="t  \t E \t from \t to \t dx \t MCs \t dE  tMC="+DoubleToStr(tMC)+"\n";
  FileWrite(f,st.c_str(),st.length());

  for(s=0;s<steps;s++)
    {
      cumdx+=dxs[s];
      if (s%writelines==0){
    st=FloatToStr(ts[s])+"\t"+DoubleToStr(energy[s])+"\t"+IntToStr(from[s])+"\t"+IntToStr(to[s])+"\t"+IntToStr(dxs[s])+"\t"+IntToStr(MCsteps[s])+"\t"+FloatToStr(des[s])+"\n";
    FileWrite(f,st.c_str(),st.length());
      }
    }

  FileClose(f);
}


void MKcurrentRandPos::heatMapToFile(string filename, int steps)
{
  int s;
  double def, deph, *relEn;
  string st;
  FILE *f;

  relEn = new double[L*L];

  for(s=0;s<steps;s++)
    {
      def = dxs[s]*Ex + dys[s]*Ey;
      deph = des[s] + def;
      relEn[from[s]] += deph/2;
      relEn[to[s]] += deph/2;
    }

  f=FileCreate(filename);

  for(s=0;s<L*L;s++)
    {
      if ((s+1)%L == 0)
    st = FloatToStr(relEn[s])+"\n";
      else
    st = FloatToStr(relEn[s])+"\t";
      FileWrite(f,st.c_str(),st.length());
    }

  FileClose(f);
}
