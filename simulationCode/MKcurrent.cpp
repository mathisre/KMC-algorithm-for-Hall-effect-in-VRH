//---------------------------------------------------------------------------

#include "fileutils.h"
#include "stringutils.h"
#include "CGanalysis.h"
#include "systemclass.h"
#include "treeutils.h"
#include "MKcurrent.h"

#include <math.h>

// Remove eventually
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>

MKcurrent::MKcurrent(int steps)
{
  from = new int[steps];
  to = new int[steps];
  energy = new double[steps];
  MCsteps = new int[steps];
  dxs = new int[steps];
  dys = new int[steps];
  dxIs = new int[steps];
  dyIs = new int[steps];
  dy3s = new int[steps];
  dx3s = new int[steps];
  ts = new double[steps];
  des = new double[steps];

  meanWeightedAreaPer = new double[steps];
  dyPer = new double[steps];


  positions  = new double*[L];
  movement   = new double*[L];
  xmovement  = new double*[L];
  ymovement  = new double*[L];
  twoSite    = new double*[L];
  x2Site     = new double*[L];
  y2Site     = new double*[L];
  triangles  = new double*[L];
  xTriangles = new double*[L];
  yTriangles = new double*[L];

  for (int x = 0; x<L; x++){
      positions[x] = new double[L];
      movement[x] = new double[L];
      xmovement[x] = new double[L];
      ymovement[x] = new double[L];
      twoSite[x] = new double[L];
      x2Site[x] = new double[L];
      y2Site[x] = new double[L];
      triangles[x] = new double[L];
      xTriangles[x] = new double[L];
      yTriangles[x] = new double[L];
  }



}

MKcurrent::~MKcurrent()
{
    if (from!=NULL) delete[] from;
      if (to!=NULL) delete[] to;
      if (energy!=NULL) delete[] energy;
      if (MCsteps!=NULL) delete[] MCsteps;
      if (dxs!=NULL) delete[] dxs;
      if (dys!=NULL) delete[] dys;
      if (dxIs!=NULL) delete[] dxIs;
      if (dyIs!=NULL) delete[] dyIs;
      if (ts!=NULL) delete[] ts;
      if (des!=NULL) delete[] des;
      if (jl!=NULL) delete[] jl;

      if (RanGen!=NULL) delete RanGen;

      if (dxInter!=NULL) delete[] dxInter;
      if (dyInter!=NULL) delete[] dyInter;
      if (dxFinal!=NULL) delete[] dxFinal;
      if (dyFinal!=NULL) delete[] dyFinal;

      if (GammaT!=NULL) delete[] GammaT;
      if (GammaT3Sites!=NULL) delete[] GammaT3Sites;
      if (jumpArea!=NULL) delete[] jumpArea;
      if (areaHistogram!=NULL) delete[] areaHistogram;


      if (testdxs!=NULL) delete[]  testdxs;
      if (testdys!=NULL) delete[]  testdys;
      if (testdxIs!=NULL) delete[] testdxIs;
      if (testdyIs!=NULL) delete[] testdyIs;


      if (discdxs!=NULL) delete[]  discdxs;
      if (discdys!=NULL) delete[]  discdys;
      if (discdxIs!=NULL) delete[] discdxIs;
      if (discdyIs!=NULL) delete[] discdyIs;

      if (transitionTracer!=NULL)  free(transitionTracer);

      for (int i=0; i < JL2; ++i){
          if (rateMatrix[i]!=NULL) delete[] rateMatrix[i];
      }
      if (rateMatrix!=NULL) delete[] rateMatrix;


      for (int x = 0; x<L; x++){
          if (positions[x] !=NULL) delete[] positions[x];
          if (movement[x]  !=NULL) delete[] movement[x];
          if (xmovement[x] !=NULL) delete[] xmovement[x];
          if (ymovement[x] !=NULL) delete[] ymovement[x];
          if (twoSite[x]   !=NULL) delete[] twoSite[x];
          if (x2Site[x]    !=NULL) delete[] x2Site[x];
          if (y2Site[x]    !=NULL) delete[] y2Site[x];
          if (triangles[x] !=NULL) delete[] triangles[x];
          if (xTriangles[x]!=NULL) delete[] xTriangles[x];
          if (yTriangles[x]!=NULL) delete[] yTriangles[x];
      }


      positions  = new double*[L];
      movement   = new double*[L];
      xmovement  = new double*[L];
      ymovement  = new double*[L];
      twoSite    = new double*[L];
      x2Site     = new double*[L];
      y2Site     = new double*[L];
      triangles  = new double*[L];
      xTriangles = new double*[L];
      yTriangles = new double*[L];


}


void MKcurrent::setE(double Ex1, double Ey1, double Ez1)
{
  Ex=Ex1; Ey=Ey1; Ez=Ez1;
}

void MKcurrent::setH(double Hx1, double Hy1, double Hz1)
{
  Hx=Hx1;
  Hy=Hy1;
  Hz=Hz1;
}

void MKcurrent::setE0(double Ex1, double Ey1, double Ez1)
{
  E0x=Ex1; E0y=Ey1; E0z=Ez1;
}

void MKcurrent::setOmega(double omega1)
{
  omega=omega1;
}

void MKcurrent::setWritelines(int writelines1)
{
  writelines=writelines1;
}

void MKcurrent::setRatefun(int ratefun1)
{
  ratefun = ratefun1;

}

void MKcurrent::setES(ESystem *e)
{
  es = e;
}


void MKcurrent::init(int D1, int L1, int N1, int maxJL1, double a1, double beta1, double maxProbability1, bool doTriJumps)
{ //  (D, L, N, p->maxJL, p->loc, 1.0/p->temp,p->maxProbability)
  int dx, dy, n, dxF, dyF, dxI, dyI;
  double gamma;
  double overlap_3, overlap_ik, overlap_jk, overlap_ij;
  double area, jlInterFinalSite, t0, tau1;
  double onePlusHA, oneMinusHA;
  RanGen = new CRandomMersenne(1);

  D=D1; L=L1; N=N1; A=a1; beta = beta1;
  maxJL = maxJL1; JL2 = 2*maxJL+1;

  maxProbability=maxProbability1/beta;
//  writelines = 10000;

//  tau1
  t0   = 1;
  tau1 = 1;

  //printf("L: %d N: %d maxJL: %d\n", L, N, maxJL);
  if (D==2)
    {
      Nmem = JL2*JL2; // Number of sites we can jump to
      jl     = new double[Nmem];
      GammaT = new double[Nmem];
      for (dx = -maxJL; dx < maxJL + 1; dx++) // This double loop works, dw
      for (dy = -maxJL; dy < maxJL + 1; dy++)
      {
        jl[(dx+maxJL) + JL2*(dy+maxJL)] = sqrt(dx*dx + dy*dy); // jump length, distance between sites
      }

        // if you read it out, then it reads like a matrix that's been put in an array


      jl[maxJL + JL2*maxJL] = 0.0;

      gamma = 0.0;
      for (n = 0; n < Nmem/2; n++)
      {
        // HERE IS WHERE TO PUT 2SITE RATE CALC
        gamma += exp(-A*jl[n])/tau1;
        GammaT[n] = gamma;
//        std::cout << GammaT[n] - GammaT[n-1] << std::endl;
        //printf("%le ",gamma); if ((n+1)%JL2==0) printf("\n");
      }
      GammaT[Nmem/2] = gamma;
      for (n = Nmem/2 + 1; n < Nmem; n++)
      {
        gamma += exp(-A*jl[n])/tau1;
        GammaT[n] = gamma;
//        std::cout << GammaT[n] - GammaT[n-1] << std::endl;
      }
      GammaTtot = gamma;

      // 3 sites rates
      Nmem3Sites = (Nmem-1)*(Nmem-1) - (Nmem-1);
      Nmem3Sites *= 2; // ta med occupied elecs


      GammaT3Sites = new double[Nmem3Sites];

      jumpArea = new double[Nmem3Sites];
      dxInter =  new int[Nmem3Sites];
      dyInter =  new int[Nmem3Sites];
      dxFinal =  new int[Nmem3Sites];
      dyFinal =  new int[Nmem3Sites];
      occupied = new bool[Nmem3Sites];

      gamma = 0.0;
      n = 0;
      double testsum = 0;
      double testsum2 = 0;
      rateMatrix = new double*[JL2];


      transitionTracer= (double*)calloc( JL2*JL2*JL2*JL2,sizeof(double));


      // dxF, dyF for final site, dxI, dyI for intermediate site
      if (doTriJumps == true)
      for (dyF= -maxJL; dyF< maxJL + 1; dyF++){
          rateMatrix[dyF + maxJL] = new double[JL2];
          for (dxF= -maxJL; dxF< maxJL + 1; dxF++){
              if( (dxF!= 0 || dyF!=0)) {
                  for (dxI= -maxJL; dxI< maxJL + 1; dxI++){
                      for (dyI= -maxJL; dyI< maxJL + 1; dyI++){
                          if ( (dxI!= 0 || dyI!=0) ){
                              if (dxF!= dxI || dyF!= dyI){



                                  area = 0.5*(dxI*dyF- dxF*dyI); // anticlockwise triangles are positive

                                  jlInterFinalSite = sqrt((dxF-dxI)*(dxF-dxI) + (dyF-dyI)*(dyF-dyI));
                                  overlap_ij = jl[(dxF+maxJL) + JL2*(dyF+maxJL)];
                                  overlap_ik = jl[(dxI+maxJL) + JL2*(dyI+maxJL)];

                                  overlap_3 = exp(-0.5*A*(overlap_ij + overlap_ik + jlInterFinalSite ));

                                  onePlusHA  = 1 + Hz*area;
                                  oneMinusHA = 1 - Hz*area;

                                  if (onePlusHA < 0){
                                      onePlusHA = 0;
                                  }
                                  if (oneMinusHA < 0){
                                      oneMinusHA = 0;
                                  }

                                  gamma += overlap_3*onePlusHA/(tau1*t0);

                                  occupied[n] = false;

                                  GammaT3Sites[n] = gamma;

                                  jumpArea[n] = area;
                                  dxFinal[n] = dxF;
                                  dyFinal[n] = dyF;
                                  dxInter[n] = dxI;
                                  dyInter[n] = dyI;
                                  n++;


                                  gamma += overlap_3*oneMinusHA/(tau1*t0);


                                  occupied[n] = true;

                                  GammaT3Sites[n] = gamma;

                                  jumpArea[n] = area;
                                  dxFinal[n] = dxF;
                                  dyFinal[n] = dyF;
                                  dxInter[n] = dxI;
                                  dyInter[n] = dyI;
                                  n++;


                                  transitionTracer[ (dyF+maxJL)*(2*maxJL+1)*(2*maxJL+1)*(2*maxJL+1) + (dxF+maxJL)*(2*maxJL+1)*(2*maxJL+1)
                                                   +(dyI+maxJL)*(2*maxJL+1)+(dxI+maxJL)] = overlap_3*( 1+ Hz*area)/(tau1*t0);


                                  rateMatrix[dyF+maxJL][dxF+maxJL] +=  overlap_3*( 1+ Hz*area)/(tau1*t0);

                              }
                          }
                      }
                  }
              }
          }
      }
      nGamma3Sites = n;
      GammaTtot3Sites = gamma;
      totGammaTtot = GammaTtot + GammaTtot3Sites;
      printf("Total 2 site rate: %.5f\n",GammaTtot);
      printf("Total 3 site rate: %.5f\n",GammaTtot3Sites);
      printf("Total rate: %.5f\n",totGammaTtot);
    }


}


void MKcurrent::setMTseed(int seed)
{
  RanGen->RandomInit(seed);
}


void inline MKcurrent::getjump(int &i, int &j, int &dx, int &dy)
{
  double r2;
  int h,l,step, x, y, x2, y2;

  i = RanGen->IRandom(0,N-1);
  while (es->getocci(i) == 0){
    i = RanGen->IRandom(0,N-1);
  }

  // r2=es->ran2(0)*GammaTtot;
  r2 = RanGen->Random()*GammaTtot;
  l = -1;
  h = Nmem-1;
  step = Nmem/2;

  while (step > 0)
  {

    if (GammaT[l+step] >= r2) {
        h = l+step;
    }
    else {
        l = l+step;
    }
    step = (h-l)/2;
  }
  //at end h=correct jump;

  x = i%L;
  y = i/L;

  dx = h%JL2 - maxJL;
  dy = h/JL2 - maxJL;

  x2 = x+dx; if (x2 >= L) x2 -= L; else if (x2 < 0) x2 += L;
  y2 = y+dy; if (y2 >= L) y2 -= L; else if (y2 < 0) y2 += L;

  j = x2 + y2 * L;
}

void inline MKcurrent::getjump3sites(int &i, int &j, int &k, int &dx, int &dy, int &dxI, int &dyI, int &jumpNumber)
{
  double r2;
  int h,l,step, x, y, x2, y2, xI, yI;


  i = RanGen->IRandom(0,N-1);
  while (es->getocci(i) == 0){
    i = RanGen->IRandom(0,N-1);
  }

  r2 = RanGen->Random()*GammaTtot3Sites;
  l = -1;
  h = Nmem3Sites-1;
  step = Nmem3Sites/2;

  while (step > 0)  {
    if (GammaT3Sites[l+step] >= r2) {
        h = l+step;
    }
    else {
        l = l+step;
    }
    step = (h-l)/2;
  }
  //at end h=correct jump;
  jumpNumber = h;

  x = i%L; // rando site
  y = i/L;

  dx = dxFinal[h];
  dy = dyFinal[h];

  x2 = x+dx; if (x2 >= L) x2 -= L; else if (x2 < 0) x2 += L;
  y2 = y+dy; if (y2 >= L) y2 -= L; else if (y2 < 0) y2 += L;

  j = x2 + y2 * L; // end position

  dxI = dxInter[h];
  dyI = dyInter[h];

  xI = x+dxI; if (xI >= L) xI -= L; else if (xI < 0) xI += L;
  yI = y+dyI; if (yI >= L) yI -= L; else if (yI < 0) yI += L;

  k = xI + yI * L; // intermediate position

//  printf("gj: %d %d %d %d h: %d r2: %le\n",i,j,dx,dy,h,r2);
}



bool inline MKcurrent::testjump(int i, int j, int dx, int dy)
{
  double dE1, dE2, rate, r2, exponent;

   if (es->getocci(j)!=0)
    {
      return false;
    }
  dE1 = es->hoppEdiffij(i,j);
  dE2 = dx*Ex + dy*Ey + dE1;
  exponent = beta*dE2;

  if(ratefun == 0){   // Approximate rate
    if (dE2 < 0) rate = 1;
    else rate = exp(-exponent);
    r2 = es->ran2(0);
    return (r2 < rate);
  }

  else if(ratefun == 1){  //Rate according to ES:
   if (dE2 < 0) rate = -dE2*(1+1/(exp(-exponent)-1));
   else if (dE2 == 0) rate = 1/beta;
   else rate = dE2/(exp(exponent)-1);
//   if (rate > maxProbability) printf("Rate: %.4f larger than max rate: %.4f\n",rate, maxProbability);

   r2 = es->ran2(0);
//   return (maxProbability*r2<rate);
   return (r2<rate);

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

bool inline MKcurrent::testjump3sites(int i, int j, int k, int dx, int dy, int dxI, int dyI, int jumpNumber)
{
  double rate, r2;
  double dEij, dEjk, dEik, dEkj, dEki;
  double rinv_ij, rinv_ik, rinv_jk;
  bool performJump;

   bool occuIntermediate = es->getocci(k) == 1;
   bool shouldBeOccupied = occupied[jumpNumber];

   if (es->getocci(j)!=0 || shouldBeOccupied != occuIntermediate)
    {
      return false;
    }
    // hoppEdiffij(i,j) gir spe(i) - spe(j) tilbake
    // min dEik betyr Shumilins dEkj

   rinv_ij = 1 / jl[(  dx+maxJL) + JL2*( dy+maxJL)];
   rinv_ik = 1 / jl[( dxI+maxJL) + JL2*(dyI+maxJL)];
   rinv_jk = 1 / jl[((dxI-dx)+maxJL) + JL2*((dyI-dy)+maxJL)];

   if (rinv_ij < es->rmaxi) rinv_ij = 0;
   if (rinv_ik < es->rmaxi) rinv_ik = 0;
   if (rinv_jk < es->rmaxi) rinv_jk = 0;

   if (!occuIntermediate)
  {
      dEij = es->hoppEdiffij(i,j) + dx*Ex + dy*Ey;
      dEik = es->hoppEdiffij(i,k) + dxI*Ex + dyI*Ey;
      dEjk = es->hoppEdiffij(j,k) + (dxI-dx)*Ex + (dyI-dy)*Ey + rinv_ij + rinv_jk - rinv_ik;
      dEkj = es->hoppEdiffij(k,j) + (dx-dxI)*Ex + (dy-dyI)*Ey + rinv_ik + rinv_jk - rinv_ij;

      (dEij < 0) ? dEij = 1 : dEij = exp(-beta*dEij);
      (dEik < 0) ? dEik = 1 : dEik = exp(-beta*dEik);
      (dEjk < 0) ? dEjk = 1 : dEjk = exp(-beta*dEjk);
      (dEkj < 0) ? dEkj = 1 : dEkj = exp(-beta*dEkj);

      rate = dEik*dEkj + dEij*dEik + dEij*dEjk;

      r2 = 3*es->ran2(0);
      performJump = r2 < rate;

      notOccupiedJumpCounter += performJump;
      dxNotOccupied+= dx*performJump;
      dyNotOccupied+= dy*performJump;
      return performJump;
  }
  else
  {
      dEij = es->hoppEdiffij(i,j) + dx*Ex + dy*Ey;
      dEik = es->hoppEdiffij(i,k) + dxI*Ex + dyI*Ey + rinv_jk + rinv_ik - rinv_ij;
      dEki = es->hoppEdiffij(k,i) - dxI*Ex - dyI*Ey - rinv_jk + rinv_ik + rinv_ij;
      dEkj = es->hoppEdiffij(k,j) + (dx-dxI)*Ex + (dy-dyI)*Ey;

      (dEij < 0) ? dEij = 1 : dEij = exp(-beta*dEij);
      (dEik < 0) ? dEik = 1 : dEik = exp(-beta*dEik);
      (dEki < 0) ? dEki = 1 : dEki = exp(-beta*dEki);
      (dEkj < 0) ? dEkj = 1 : dEkj = exp(-beta*dEkj);

      rate = dEik*dEkj + dEij*dEki + dEij*dEkj;
      r2 = 3*es->ran2(0);
      performJump = r2 < rate;




      occupiedJumpCounter += performJump;
      dxOccupied += dx*performJump;
      dyOccupied += dy*performJump;
      return performJump;
  }

}




void MKcurrent::runCurrent(int steps,double &E, double &t)
{

  int s, MCs, i=0, j=0, dx=0, dy=0, jumpNumber = 0;
  int k, dxI=0, dyI=0;
  double dE,dt, jumpLength;
  bool jump;
  double prob2Site = GammaTtot/totGammaTtot;
  printf("prob2site: %.3f\n", prob2Site);


    tMC = 1/(               N * es->nu() * totGammaTtot);


  meanArea = 0;  meanWeightedArea = 0;  meanJumpLength = 0;
  meanInterJumpLength = 0;  meanWeightedJL_x = 0;  meanWeightedJL_y = 0;

  occupiedJumpCounter =  0;  notOccupiedJumpCounter = 0;
  dyOccupied = 0; dyNotOccupied = 0;  dxOccupied = 0; dxNotOccupied = 0;

  meandx = 0;  meandy = 0;  meandxI = 0;  meandyI = 0;
  testedNumberOf2Site = 0;  testedNumberOf3Site = 0;
  numberOf2Site = 0;  numberOf3Site = 0;

  zeroAreaCounter = 0;  meanSomething = 0;

  meanInterFinalJumpLength = 0;
  meanDiscdx =  0; meanDiscdy  = 0;
  meanAbsdx =   0; meanAbsdy   = 0;
  meanDiscdxI = 0; meanDiscdyI = 0;
  meanAbsdxI =  0; meanAbsdyI  = 0;

  double meanMCs = 0;
  double actualmeanMCs = 0;


  dy3 = 0;  dx3 = 0;  dy2 = 0;  dx2 = 0;
  meanJumpLength3 = 0;  meanJumpLength2 = 0;

  test3JumpCounter = 0;
  disc3JumpCounter = 0;

  bool threeSiteJump = false;

  for (s = 0; s < steps; s++)
  {
      MCs = 0;
      jump = false;
      while (!jump)
      {
              MCs++;
              if ((es->ran2(0))<prob2Site){
                testedNumberOf2Site++;
                getjump(i,j,dx,dy);
                jump=testjump(i,j,dx,dy);
                sample2SiteJump(dx,dy,jump);
                dxI = 0; dyI = 0;
                threeSiteJump = false;
              }
              else{
                getjump3sites(i,j,k,dx,dy,dxI,dyI, jumpNumber);
                jump=testjump3sites(i,j,k,dx,dy,dxI,dyI,jumpNumber);
                sample3SiteJump(dx, dy, dxI, dyI, jumpNumber, jump);
                threeSiteJump = true;
              }
      }

      dE = es->hopp(i,-1,0,j,0,0);

      dt = MCs*tMC;
      E += dE;
      t += dt;
      meanMCs += MCs;

      jumpLength = jl[(dx+maxJL) + JL2*(dy+maxJL)];

      meanJumpLength      += jumpLength;
      meanWeightedJL_x    += jumpLength*dx;

      meandx  += dx;
      meandy  += dy;
      meandxI += dxI;
      meandyI += dyI;

      actualmeanMCs += MCs;

      from[s] = i;
      to[s]   = j;

      ts[s] = t;
      energy[s] = E;
      des[s] = dE;
      dxs[s] = dx;
      dys[s] = dy;
      dx3s[s] = dx*threeSiteJump;
      dy3s[s] = dy*threeSiteJump;

      dxIs[s] = dxI;
      dyIs[s] = dyI;

      MCsteps[s]=MCs;

      if (s%100000==0) {
          printf("%6d E=%le MCs=%.3f\n",s,E,actualmeanMCs/100000);
          actualmeanMCs =0;
      }

     // printf("Step: %d: %d %d dE: %le I:%le MCs: %d\n", s, i, j,energy[s],dx,MCsteps[s]);
    }

  printf("\n-------Finishing timesteps--------\n\n");
  finishTimeStepsAndSample(steps);


}


void MKcurrent::finishTimeStepsAndSample(int steps)
{


    meanJumpLength /= steps;
    meanJumpLength3 /= numberOf3Site;
    meanJumpLength2 /= numberOf2Site;
    meanInterJumpLength /= numberOf3Site;
    meanInterFinalJumpLength /= numberOf3Site;

    meanWeightedArea /= dy3;


    zeroAreaCounter /= numberOf3Site;
    meanArea /= numberOf3Site-zeroAreaCounter;

    meanWeightedJL_x    /= meandx;
    meanWeightedJL_y    /= dy3;

    meandx /= steps;  meandy /= steps;
    meandxI /= steps;  meandyI /= steps;

    meanAbsdx /= testedNumberOf2Site + testedNumberOf3Site;
    meanAbsdy /= testedNumberOf2Site + testedNumberOf3Site;

    meanDiscdx /= (testedNumberOf2Site - numberOf2Site) + (testedNumberOf3Site - numberOf3Site);
    meanDiscdy /= (testedNumberOf2Site - numberOf2Site) + (testedNumberOf3Site - numberOf3Site);

    meanAbsdxI /=testedNumberOf3Site;
    meanAbsdyI /=testedNumberOf3Site;

    meanDiscdxI /= (testedNumberOf3Site - numberOf3Site);
    meanDiscdyI /= (testedNumberOf3Site - numberOf3Site);

    double meandxOccupiedJump = dxOccupied/occupiedJumpCounter;
    double meandyOccupiedJump = dyOccupied/occupiedJumpCounter;

    double meandxNotOccupiedJump = dxNotOccupied/notOccupiedJumpCounter;
    double meandyNotOccupiedJump = dyNotOccupied/notOccupiedJumpCounter;

    occupiedJumpCounter /= numberOf3Site;
    notOccupiedJumpCounter /= numberOf3Site;

    dxOccupied /= numberOf3Site;
    dyOccupied /= numberOf3Site;
    dxNotOccupied /= numberOf3Site;
    dyNotOccupied /= numberOf3Site;

    driftAngle = atan2(meandy, meandx)*180/M_PI;

    printf("Number of 3 site jumps tested   : %d\n", testedNumberOf3Site);
    printf("Number of 3 site jumps performed: %d\n", numberOf3Site);

    printf("Mean 2-site dx: %f, dy: %f\n", dx2/numberOf2Site, dy2/numberOf2Site);
    printf("Mean 3-site dx: %f, dy: %f\n", dx3/numberOf3Site, dy3/numberOf3Site);


    printf("Mean dE: %f \n", meanSomething/testedNumberOf2Site);

    printf("\nProposed  dx:  %8.5f, dy:  %8.5f\n", meanAbsdx , meanAbsdy );
    printf("Discarded dx:  %8.5f, dy:  %8.5f\n", meanDiscdx, meanDiscdy);
    printf("Mean dx:       %8.5f, dy:  %8.5f\n", meandx, meandy);


    printf("\nProposed  dxI: %8.5f, dyI:  %8.5f\n", meanAbsdxI , meanAbsdyI );
    printf("Discarded dxI: %8.5f, dyI:  %8.5f\n",  meanDiscdxI, meanDiscdyI);
    printf("Mean dxI:      %8.5f, dyI:  %8.5f\n", meandxI, meandyI);

    printf("\nMean jump length:                    %8.5f\n", meanJumpLength);
    printf("Mean intermediate jump length:       %8.5f\n", meanInterJumpLength);
    printf("Mean intermediate-final jump length: %8.5f\n", meanInterFinalJumpLength);
    printf("Mean jump length, only 2 site:       %8.5f\n", meanJumpLength2);
    printf("Mean jump length, only 3 site:       %8.5f\n", meanJumpLength3);

    printf("\nPercentage zero area:      %8.5f\n", zeroAreaCounter);
    printf("Mean area:                 %8.5f\n", meanArea);
    printf("Mean current area:         %8.5f\n", meanWeightedArea);

    printf("\nIntermediate occupied: %5.4f, unoccupied: %5.4f\n", occupiedJumpCounter, notOccupiedJumpCounter);
    printf("Used occupied     dx: %5.4f, dy: %5.4f\n", dxOccupied, dyOccupied);
    printf("Used not occupied dx: %5.4f, dy: %5.4f\n", dxNotOccupied, dyNotOccupied);
    printf("Mean occupied     dx: %5.4f, dy: %5.4f\n", meandxOccupiedJump, meandyOccupiedJump);
    printf("Mean not occupied dx: %5.4f, dy: %5.4f\n", meandxNotOccupiedJump, meandyNotOccupiedJump);

    printf("\n Drift angle: %4f\n", driftAngle);

    printf("Acceptance ratio for 2 site: %.3f, 3 site: %.3f\n", double(numberOf2Site)/testedNumberOf2Site, double(numberOf3Site)/testedNumberOf3Site);
    printf("Percentage number of 3 site jumps: %f\n", double(numberOf3Site)/(numberOf2Site+numberOf3Site));


}

void MKcurrent::initPositions(string filename)
{
    int i,j=0;
    // i is position number, j is electron number
    for (i = 0; i<N; i++){
        if (es->getocci(i)!=0){
            electronPositions[j][0] = j;
            electronPositions[j][1] = i%L;
            electronPositions[j][2] = i/L;
            electronPositions[j][3] = i;
            j++;
        }
    }
    string st;
    FILE* f;
    f=FileCreate(filename);
    FileClose(f);
}

void MKcurrent::updatePositions(int i, int j, int dx, int dy)
{
    int k=0;
    while (electronPositions[k][3] != i) k++;

    electronPositions[k][1] += dx;
    electronPositions[k][2] += dy;
    electronPositions[k][3] = j;

    if      (electronPositions[k][1] >  L) electronPositions[k][1] -= L;
    else if (electronPositions[k][1] <  0) electronPositions[k][1] += L;

    if      (electronPositions[k][2] >  L) electronPositions[k][2] -= L;
    else if (electronPositions[k][2] <  0) electronPositions[k][2] += L;
}

void MKcurrent::writePositions(string filename, int s)
{

    string st;
    FILE* g;
    g = FileContinue(filename);
    int i;
    st="ITEM: TIMESTEP\n" + IntToStr(s/100) + "\nITEM: NUMBER OF ATOMS\n5000\nITEM: BOX BOUNDS pp pp pp\n0 100\n0 100\n0 100\nITEM: ATOMS id x y z\n";
    FileWrite(g,st.c_str(),st.length());

    for(i = 0; i < N/2; i++){
        st= IntToStr(electronPositions[i][0]) + " " + IntToStr(electronPositions[i][1])
                + " " + IntToStr(electronPositions[i][2]) + " " + DoubleToStr(50*(es->spe[electronPositions[i][3]]+1)) +"\n";
        FileWrite(g,st.c_str(),st.length());
    }
    FileClose(g);

}

void MKcurrent::closePositions(FILE* f)
{
    FileClose(f);
}

void MKcurrent::sample2SiteJump(int dx, int dy, bool jump)
{
    meanAbsdx +=dx; meanAbsdy +=dy;
    if (jump == true) {
        numberOf2Site++;

        meanJumpLength2 += jl[(dx+maxJL) + JL2*(dy+maxJL)];
        dx2 += dx;
        dy2 += dy;
    }
    else meanDiscdx += dx; meanDiscdy += dy;
}

void MKcurrent::sample3SiteJump(int dx, int dy, int dxI, int dyI, int jumpNumber, bool jump)
{
    int jumpLength;

    meanAbsdx += dx; meanAbsdy += dy;
    meanAbsdxI += dxI; meanAbsdyI += dyI;

//    testdxs[ testedNumberOf3Site] = dx;
//    testdys[ testedNumberOf3Site] = dy;
//    testdyIs[testedNumberOf3Site] = dyI;
//    testdxIs[testedNumberOf3Site] = dxI;
    testedNumberOf3Site++;

    if (jump == true) {

        meanInterFinalJumpLength += sqrt((dx-dxI)*(dx-dxI) + (dy-dyI)*(dy-dyI));


        if (jumpArea[jumpNumber] == 0) zeroAreaCounter++;
        else {
            meanArea += fabs(jumpArea[jumpNumber]);
            meanWeightedArea += fabs(jumpArea[jumpNumber])*dy;
        }
        jumpLength = jl[(dx+maxJL) + JL2*(dy+maxJL)];
        meanJumpLength3  += jumpLength;
        meanWeightedJL_y += jumpLength*dy;

        meanInterJumpLength += jl[(dxI+maxJL) + JL2*(dyI+maxJL)];
        dx3 += dx;
        dy3 += dy;

        numberOf3Site++;
         }
    else {
        meanDiscdx += dx; meanDiscdy += dy;
        meanDiscdxI += dxI; meanDiscdyI += dyI;

//        discdxs[ disc3JumpCounter] = dx;
//        discdys[ disc3JumpCounter] = dy;
//        discdxIs[disc3JumpCounter] = dxI;
//        discdyIs[disc3JumpCounter] = dyI;
        disc3JumpCounter++;
    }
}


void MKcurrent::normalizeMap(double** &vector, int steps)
{
    int x, y, i;
    for (i = 0; i < N; i++){
        x = i%L;
        y = i/L;
        vector[y][x] /= steps;
    }
}

void MKcurrent::updateMovement(double** &vector, int i, int j, int dWhat)
{
    int x0, y0, x1, y1;
    x0 = i%L; y0 = i/L; x1 = j%L; y1 = j/L;
    vector[y0][x0] += dWhat;
    vector[y1][x1] += dWhat;
}



void MKcurrent::updateTriangles(double** &vector, int i, int k, int j, int dWhat)
{
    int x0, y0, x1, y1, x2, y2;
    x0 = i%L; y0 = i/L; x1 = j%L; y1 = j/L; x2 = k%L; y2 = k/L;
    vector[y0][x0] += dWhat;
    vector[y1][x1] += dWhat;
    vector[y2][x2] += dWhat;

}


void MKcurrent::update2Site(double** &vector, int i, int j, int dWhat)
{
    int x0, y0, x1, y1;
    x0 = i%L; y0 = i/L; x1 = j%L; y1 = j/L;

    vector[y0][x0] += dWhat;
    vector[y1][x1] += dWhat;

}


void MKcurrent::writeMapToFile(double **vector, int size, string filename)
{

    unsigned int x, y, L1;
    L1 = unsigned(size);
    ofstream outfileMap;
    outfileMap.open(filename.c_str());
    char string[50];

    for (x = 0; x < L1; x++){
        for (y = 0; y < L1; y++){
            sprintf(string, "%.4e",vector[x][y]);
            outfileMap << string << " ";
        }
        outfileMap << "\n";
    }
}

void MKcurrent::saveMaps(string outputprefix, string outputendfix, int steps)
{
    normalizeMap(triangles,  numberOf3Site);
    normalizeMap(xTriangles, numberOf3Site);
    normalizeMap(yTriangles, numberOf3Site);

    normalizeMap(twoSite,    numberOf2Site);
    normalizeMap(x2Site,     numberOf2Site);
    normalizeMap(y2Site,     numberOf2Site);

    normalizeMap(movement,   steps);
    normalizeMap(xmovement,   steps);
    normalizeMap(ymovement,   steps);
    normalizeMap(positions,  steps);

    writeMapToFile(triangles,  L, outputprefix   +  "3SiteMap_" + outputendfix    + ".dat");
    writeMapToFile(xTriangles, L, outputprefix   + "x3SiteMap_" + outputendfix    + ".dat");
    writeMapToFile(yTriangles, L, outputprefix   + "y3SiteMap_" + outputendfix    + ".dat");

    writeMapToFile(x2Site, L,  outputprefix      + "x2SiteMap_" + outputendfix    + ".dat");
    writeMapToFile(y2Site, L,  outputprefix      + "y2SiteMap_" + outputendfix    + ".dat");
    writeMapToFile(twoSite,  L,  outputprefix    + "2SiteMap_" +  outputendfix    + ".dat");

    writeMapToFile(movement,  L,  outputprefix   +   "allJumps_" +  outputendfix  + ".dat");
    writeMapToFile(xmovement,  L,  outputprefix  +  "xallJumps_" +  outputendfix  + ".dat");
    writeMapToFile(ymovement,  L,  outputprefix  +  "yallJumps_" +  outputendfix  + ".dat");
    writeMapToFile(positions, L,  outputprefix   + "electronMap_" + outputendfix  + ".dat");


}

void MKcurrent::jumpsToFileSmall(string filename, int steps)
{
  int s, cumdx, cumdy;
  int cumdx3, cumdy3;

  FILE* f;
  string st;

  cumdx  = 0;  cumdy =  0;
  cumdx3 = 0;  cumdy3 = 0;
  f=FileCreate(filename);



  st="t  \t \t E \t \t from \t to \t dx \t dy \t dx3 \t dy3 \t MCs \t dE \t tMC="+DoubleToStr(tMC)+ "\t <jl>="
          + DoubleToStr(meanJumpLength) +"\t <dx>=" + DoubleToStr(meandx) +"\t <dy>="
          + DoubleToStr(meandy) + "\t <jl_I>=" +DoubleToStr(meanInterJumpLength)
          +"\t <dx_I>=" + DoubleToStr(meandxI) +"\t <dy_I>=" + DoubleToStr(meandyI)
          + "\t Acc 2=" + DoubleToStr(double(numberOf2Site)/testedNumberOf2Site)
          + "\t Acc 3=" + DoubleToStr(double(numberOf3Site)/testedNumberOf3Site)
          + "\t meanArea=" + DoubleToStr(meanArea)+ "\t meanCurrentArea=" + DoubleToStr(meanWeightedArea)
          + "\t mean_x_CurrWeightedJL=" + DoubleToStr(meanWeightedJL_x)
          + "\t mean_y_CurrWeightedJL=" + DoubleToStr(meanWeightedJL_y) + "\n";
  FileWrite(f,st.c_str(),st.length());



  for(s = 0; s < steps; s++){
      cumdx += dxs[s];
      cumdy += dys[s];
      cumdx3 += dx3s[s];
      cumdy3 += dy3s[s];

      if (s%writelines == 0){
        st = FloatToStr(ts[s])+"\t"+DoubleToStr(energy[s])+"\t"+IntToStr(from[s])+"\t"+IntToStr(to[s])+"\t"+IntToStr(cumdx)+"\t"+IntToStr(cumdy)+"\t"+IntToStr(cumdx3)+"\t"+IntToStr(cumdy3)+"\t"+IntToStr(MCsteps[s])+"\t"+FloatToStr(des[s])+"\n";
        FileWrite(f,st.c_str(),st.length());
      }
    }

  printf("Drift Angle: %f\n", driftAngle);

  printf("Rough estimate for conductivity_x: %8.5f\n", double(cumdx)/steps);
  printf("Rough estimate for conductivity_y: %8.5f\n", double(cumdy)/steps);

  FileClose(f);
}

void MKcurrent::jumpsToFileAC(string filename, int steps)
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


void MKcurrent::heatMapToFile(string filename, int steps)
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


void MKcurrent::sample3SiteJumpRates(int steps, string proposedFilename, string discardedFilename, string performedFilename)
{
    int s, i,j,iI,jI;
    double **jumpsPerformed,**counterPerformed, **tracerSumPerformed;
    double **jumpsTested, **counterTested, **tracerSumTested;
    double **jumpsDiscarded, **counterDiscarded, **tracerSumDiscarded;

    int L=3;
    // Tracer tracks what intermediate jump went to what final site
        jumpsPerformed = new double*[2*L+1];
      counterPerformed = new double*[2*L+1];
    tracerSumPerformed = new double*[2*L+1];

        jumpsDiscarded = new double*[2*L+1];
      counterDiscarded = new double*[2*L+1];
    tracerSumDiscarded = new double*[2*L+1];

        jumpsTested = new double*[2*L+1];
      counterTested = new double*[2*L+1];
    tracerSumTested = new double*[2*L+1];



    for (i = 0; i < 2*L+1; i++) {
            jumpsPerformed[i] = new double[2*L+1];
          counterPerformed[i] = new double[2*L+1];
        tracerSumPerformed[i] = new double[2*L+1];

            jumpsDiscarded[i] = new double[2*L+1];
          counterDiscarded[i] = new double[2*L+1];
        tracerSumDiscarded[i] = new double[2*L+1];

            jumpsTested[i] = new double[2*L+1];
          counterTested[i] = new double[2*L+1];
        tracerSumTested[i] = new double[2*L+1];

        for (j =0; j<2*L+1; ++j){

            jumpsPerformed[i][j] = 0;
          counterPerformed[i][j] = 0;
        tracerSumPerformed[i][j] = 0;

            jumpsDiscarded[i][j] = 0;
          counterDiscarded[i][j] = 0;
        tracerSumDiscarded[i][j] = 0;

            jumpsTested[i][j] = 0;
          counterTested[i][j] = 0;
        tracerSumTested[i][j] = 0;
        }
    }

    double *tracerPerformed = (double*)calloc( (2*L+1) * (2*L+1)* (2*L+1)* (2*L+1),sizeof(double));
    double *tracerDiscarded = (double*)calloc( (2*L+1) * (2*L+1)* (2*L+1)* (2*L+1),sizeof(double));
    double *tracerTested    = (double*)calloc( (2*L+1) * (2*L+1)* (2*L+1)* (2*L+1),sizeof(double));

    // Trace performed jumps


    // Jumps is a 2d matrice with the normalized number of jumps going to individual final sites. Jumps only considers final sites within L, but also considers intermediates sites outside.
    // Counter is the same but for the intermediate sites. How many intermediate jumps within L go to final site
    // Tracer is the 4d matrice for each final site containing the matrice of what intermediate sites were used with what rate.
    for (s=0; s < steps; s++){
        if (dxIs[s] != 0 || dyIs[s] != 0){ // If we did 3 jump in step s
            if (fabs(dxs[s]) <= L && fabs(dys[s]) <= L ){
                jumpsPerformed[dys[s]+L][dxs[s]+L]++;
                if (fabs(dxIs[s]) <= L && fabs(dyIs[s]) <= L ){
                    counterPerformed[dys[s]+L][dxs[s]+L]++;
                    tracerPerformed[  (dys[s]+L)*(2*L+1)*(2*L+1)*(2*L+1) +
                             (dxs[s]+L)*(2*L+1)*(2*L+1) +
                            (dyIs[s]+L)*(2*L+1)+
                            (dxIs[s]+L)]++;
                }
            }
        }
    }
    // Trace proposed (tested) jumps
    for (s=0; s < testedNumberOf3Site; s++){
            if (testdxIs[s] != 0 || testdyIs[s] != 0){ // If we did 3 jump in step s
                if (fabs(testdxs[s]) <= L && fabs(testdys[s]) <= L ){
                    jumpsTested[testdys[s]+L][testdxs[s]+L]++;
                    if (fabs(testdxIs[s]) <= L && fabs(testdyIs[s]) <= L ){
                        counterTested[testdys[s]+L][testdxs[s]+L]++;
                        tracerTested[  (testdys[s]+L)*(2*L+1)*(2*L+1)*(2*L+1) +
                                 (testdxs[s]+L)*(2*L+1)*(2*L+1) +
                                (testdyIs[s]+L)*(2*L+1)+
                                (testdxIs[s]+L)]++;
                    }
                }
            }
        }

    // Trace discarded jumps
    for (s=0; s < disc3JumpCounter; s++){
            if (discdxIs[s] != 0 || discdyIs[s] != 0){ // If we did 3 jump in step s
                if (fabs(discdxs[s]) <= L && fabs(discdys[s]) <= L ){
                    jumpsDiscarded[discdys[s]+L][discdxs[s]+L]++;
                    if (fabs(discdxIs[s]) <= L && fabs(discdyIs[s]) <= L ){
                        counterDiscarded[discdys[s]+L][discdxs[s]+L]++;
                        tracerDiscarded[  (discdys[s]+L)*(2*L+1)*(2*L+1)*(2*L+1) +
                                 (discdxs[s]+L)*(2*L+1)*(2*L+1) +
                                (discdyIs[s]+L)*(2*L+1)+
                                (discdxIs[s]+L)]++;
                    }
                }
            }
        }



    fflush(stdout);
    printf("\n");
    char *st = "FINAL";
    char *stc = "START";


    printf("\nUsed 3-site jump probability matrices\n");

    // Tested, discarded used
    for (i=-L; i<L+1; i++){
            for (j = -L; j < L+1; j++) {
                printf(" Proposed \t \t\t \t \t \t\t Discarded \t  \t \t \t \t  \t \t Used\n");
                for (jI=L; jI>-L-1; jI--) {
                    for (iI=-L; iI<L+1; iI++){
                        if (i+L!=L || j+L!=L){
                            tracerTested[(j+L)*(2*L+1)*(2*L+1)*(2*L+1) +
                                                              (i+L)*(2*L+1)*(2*L+1) +
                                                             (jI+L)*(2*L+1)+
                                                             (iI+L)] /= jumpsTested[j+L][i+L];
                        }
                        if (i==iI && j==jI){
                            tracerSumTested[j+L][i+L] += jumpsTested[j+L][i+L]/testedNumberOf3Site;
                            printf("%8s",st);
                        }
                        else if (iI==0 && jI==0) printf("%8s", stc);
                        else printf("%8.5f",tracerTested[(j+L)*(2*L+1)*(2*L+1)*(2*L+1) +
                                                            (i+L)*(2*L+1)*(2*L+1) +
                                                           (jI+L)*(2*L+1)+
                                                           (iI+L)]);
                    }

                    printf("\t");
                    for (iI=-L; iI<L+1; iI++){
                        if (i+L!=L || j+L!=L){
                            tracerDiscarded[(j+L)*(2*L+1)*(2*L+1)*(2*L+1) +
                                                              (i+L)*(2*L+1)*(2*L+1) +
                                                             (jI+L)*(2*L+1)+
                                                             (iI+L)] /= jumpsDiscarded[j+L][i+L];
                        }
                        if (i==iI && j==jI){
                            tracerSumDiscarded[j+L][i+L] += jumpsDiscarded[j+L][i+L]/disc3JumpCounter;
                            printf("%8s",st);
                        }
                        else if (iI==0 && jI==0) printf("%8s", stc);
                        else printf("%8.5f",tracerDiscarded[(j+L)*(2*L+1)*(2*L+1)*(2*L+1) +
                                                            (i+L)*(2*L+1)*(2*L+1) +
                                                           (jI+L)*(2*L+1)+
                                                           (iI+L)]);
                    }
                    printf("\t");
                    for (iI=-L; iI<L+1; iI++){
                        if (i+L!=L || j+L!=L){
                            tracerPerformed[(j+L)*(2*L+1)*(2*L+1)*(2*L+1) +
                                                              (i+L)*(2*L+1)*(2*L+1) +
                                                             (jI+L)*(2*L+1)+
                                                             (iI+L)] /= jumpsPerformed[j+L][i+L];
                        }
                        if (i==iI && j==jI){
                            tracerSumPerformed[j+L][i+L] += jumpsPerformed[j+L][i+L]/numberOf3Site;
                            printf("%8s",st);
                        }
                        else if (iI==0 && jI==0) printf("%8s", stc);
                        else printf("%8.5f",tracerPerformed[(j+L)*(2*L+1)*(2*L+1)*(2*L+1) +
                                                            (i+L)*(2*L+1)*(2*L+1) +
                                                           (jI+L)*(2*L+1)+
                                                           (iI+L)]);
                    }
                    printf("\n");
                }
                printf("Jumps included:   %8.5f",  1-(jumpsTested[j+L][i+L]-counterTested[j+L][i+L])/jumpsTested[j+L][i+L]  );
                printf("\t \t   \t  \t \t");
                printf("Jumps included:   %8.5f",  1-(jumpsDiscarded[j+L][i+L]-counterDiscarded[j+L][i+L])/jumpsDiscarded[j+L][i+L]  );
                printf("\t \t   \t  \t \t");
                printf("Jumps included:   %8.5f",  1-(jumpsPerformed[j+L][i+L]-counterPerformed[j+L][i+L])/jumpsPerformed[j+L][i+L]  );
                printf("\n");

                printf("Jumps going here: %8.5f",tracerSumTested[j+L][i+L] );
                printf("\t \t \t   \t \t");
                printf("Jumps going here: %8.5f",tracerSumDiscarded[j+L][i+L] );
                printf("\t \t \t  \t \t");
                printf("Jumps going here: %8.5f",tracerSumPerformed[j+L][i+L] );
                printf("\n\n");
            }
        }

    printf("\n\n---------------------------------------------------------------------------------------------------------------------------------------------------------\n");
    // Up down left right
    double sumPerformed=0, sum_upPerformed=0, sum_downPerformed=0, sum_leftPerformed=0, sum_rightPerformed=0;
    double sumTested=0, sum_upTested=0, sum_downTested=0, sum_leftTested=0, sum_rightTested=0;
    double sumDiscarded=0, sum_upDiscarded=0, sum_downDiscarded=0, sum_leftDiscarded=0, sum_rightDiscarded=0;
    // Corners
    double sum_upleftTested=0, sum_uprightTested=0, sum_downleftTested=0, sum_downrightTested=0;
    double sum_upleftPerformed=0, sum_uprightPerformed=0, sum_downleftPerformed=0, sum_downrightPerformed=0;
    double sum_upleftDiscarded=0, sum_uprightDiscarded=0, sum_downleftDiscarded=0, sum_downrightDiscarded=0;
    // Straight up down left right
    double sum_onlydownTested=0, sum_onlyupTested=0, sum_onlyleftTested=0, sum_onlyrightTested=0;
    double sum_onlydownDiscarded=0, sum_onlyupDiscarded=0, sum_onlyleftDiscarded=0, sum_onlyrightDiscarded=0;
    double sum_onlydownPerformed=0, sum_onlyupPerformed=0, sum_onlyleftPerformed=0, sum_onlyrightPerformed=0;

    FILE* proposedFile;
    FILE* discardedFile;
    FILE* performedFile;

    proposedFile =  FileCreate(proposedFilename);
    discardedFile = FileCreate(discardedFilename);
    performedFile = FileCreate(performedFilename);
    string proposedText, discardedText, performedText;



    printf("\n3-site jump probability matrix: \n");
    printf(" Proposed \t \t\t \t \t \t\t Discarded  \t \t \t \t \t \t \t Used\n");
    for (j=L; j>-L-1; j--) {
        proposedText=""; discardedText=""; performedText="";
        for (i=-L; i<L+1; i++){

             jumpsTested[j+L][i+L] /= testedNumberOf3Site;
             if (i==0 && j==0){
                 printf("%8s", stc);
                 proposedText += "0.000 ";

             }
             else{
                 printf("%8.5f", jumpsTested[j+L][i+L]);
                 proposedText += DoubleToStr(jumpsTested[j+L][i+L]) + " ";
             }

             sumTested += jumpsTested[j+L][i+L];
             if (j>0) sum_upTested         += jumpsTested[j+L][i+L];
             else if (j<0) sum_downTested  += jumpsTested[j+L][i+L];
             if (i<0) sum_leftTested       += jumpsTested[j+L][i+L];
             else if (i>0) sum_rightTested += jumpsTested[j+L][i+L];

             if (j>0 && i < 0) sum_upleftTested    += jumpsTested[j+L][i+L];
             if (j>0 && i > 0) sum_uprightTested   += jumpsTested[j+L][i+L];
             if (j<0 && i > 0) sum_downrightTested += jumpsTested[j+L][i+L];
             if (j<0 && i < 0) sum_downleftTested  += jumpsTested[j+L][i+L];

             if (j == 0 && i < 0) sum_onlyleftTested  += jumpsTested[j+L][i+L];
             if (j == 0 && i > 0) sum_onlyrightTested += jumpsTested[j+L][i+L];
             if (i == 0 && j > 0) sum_onlyupTested    += jumpsTested[j+L][i+L];
             if (i == 0 && j < 0) sum_onlydownTested  += jumpsTested[j+L][i+L];

         }


        printf("\t");
        for (i=-L; i<L+1; i++){

             jumpsDiscarded[j+L][i+L] /= disc3JumpCounter;
             sumDiscarded += jumpsDiscarded[j+L][i+L];
             if (i==0 && j==0){
                 printf("%8s", stc);
                 discardedText += "0.000 ";

             }
             else{
                 printf("%8.5f", jumpsDiscarded[j+L][i+L]);
                 discardedText += DoubleToStr(jumpsDiscarded[j+L][i+L]) + " ";
             }
             if (j>0) sum_upDiscarded         += jumpsDiscarded[j+L][i+L];
             else if (j<0) sum_downDiscarded  += jumpsDiscarded[j+L][i+L];
             if (i<0) sum_leftDiscarded       += jumpsDiscarded[j+L][i+L];
             else if (i>0) sum_rightDiscarded += jumpsDiscarded[j+L][i+L];

             if (j>0 && i < 0) sum_upleftDiscarded    += jumpsDiscarded[j+L][i+L];
             if (j>0 && i > 0) sum_uprightDiscarded   += jumpsDiscarded[j+L][i+L];
             if (j<0 && i > 0) sum_downrightDiscarded += jumpsDiscarded[j+L][i+L];
             if (j<0 && i < 0) sum_downleftDiscarded   += jumpsDiscarded[j+L][i+L];


             if (j == 0 && i < 0) sum_onlyleftDiscarded  += jumpsDiscarded[j+L][i+L];
             if (j == 0 && i > 0) sum_onlyrightDiscarded += jumpsDiscarded[j+L][i+L];
             if (i == 0 && j > 0) sum_onlyupDiscarded    += jumpsDiscarded[j+L][i+L];
             if (i == 0 && j < 0) sum_onlydownDiscarded  += jumpsDiscarded[j+L][i+L];
         }

        printf("\t");
         for (i=-L; i<L+1; i++){

             jumpsPerformed[j+L][i+L] /= numberOf3Site;
             sumPerformed += jumpsPerformed[j+L][i+L];
             if (i==0 && j==0){
                 printf("%8s", stc);
                 performedText += "0.000 ";

             }
             else{
                 printf("%8.5f", jumpsPerformed[j+L][i+L]);
                 performedText += DoubleToStr(jumpsPerformed[j+L][i+L]) + " ";
             }
             if (j>0) sum_upPerformed         += jumpsPerformed[j+L][i+L];
             else if (j<0) sum_downPerformed  += jumpsPerformed[j+L][i+L];
             if (i<0) sum_leftPerformed       += jumpsPerformed[j+L][i+L];
             else if (i>0) sum_rightPerformed += jumpsPerformed[j+L][i+L];

             if (j>0 && i < 0) sum_upleftPerformed     += jumpsPerformed[j+L][i+L];
             if (j>0 && i > 0) sum_uprightPerformed    += jumpsPerformed[j+L][i+L];
             if (j<0 && i > 0) sum_downrightPerformed  += jumpsPerformed[j+L][i+L];
             if (j<0 && i < 0) sum_downleftPerformed    += jumpsPerformed[j+L][i+L];


             if (j == 0 && i < 0) sum_onlyleftPerformed  += jumpsPerformed[j+L][i+L];
             if (j == 0 && i > 0) sum_onlyrightPerformed += jumpsPerformed[j+L][i+L];
             if (i == 0 && j > 0) sum_onlyupPerformed    += jumpsPerformed[j+L][i+L];
             if (i == 0 && j < 0) sum_onlydownPerformed  += jumpsPerformed[j+L][i+L];

         }

         proposedText += "\n";
         discardedText += "\n";
         performedText += "\n";


         FileWrite(proposedFile, proposedText.c_str(), proposedText.length());
         FileWrite(discardedFile,discardedText.c_str(),discardedText.length());
         FileWrite(performedFile,performedText.c_str(),performedText.length());
         printf("\n");


        }


    FileClose(proposedFile );
    FileClose(discardedFile);
    FileClose(performedFile);

    printf("\n");

    printf("Prob 3 sum: %8.5f",sumTested);
    printf("\t \t   \t \t  \t \t ");
    printf("Prob 3 sum: %8.5f",sumDiscarded);
    printf("\t \t   \t \t \t \t ");
    printf("Prob 3 sum: %8.5f",sumPerformed);
    printf("\n");

    printf("Prob left:  %8.5f Prob up:   %8.5f", sum_leftTested, sum_upTested);
    printf("\t \t     \t ");
    printf("Prob left:  %8.5f Prob up:   %8.5f", sum_leftDiscarded, sum_upDiscarded);
    printf("\t \t     \t ");
    printf("Prob left:  %8.5f Prob up:   %8.5f", sum_leftPerformed, sum_upPerformed);
    printf("\n");

    printf("Prob right: %8.5f Prob down: %8.5f", sum_rightTested, sum_downTested);
    printf("\t  \t \t ");
    printf("Prob right: %8.5f Prob down: %8.5f", sum_rightDiscarded, sum_downDiscarded);
    printf("\t  \t \t ");
    printf("Prob right: %8.5f Prob down: %8.5f", sum_rightPerformed, sum_downPerformed);
    printf("\n");


    printf("Prob up left:    %8.5f Prob up right:     %8.5f", sum_upleftTested, sum_uprightTested);
    printf("\t         ");
    printf("Prob up left:    %8.5f Prob up right:     %8.5f", sum_upleftDiscarded, sum_uprightDiscarded);
    printf("\t         ");
    printf("Prob up left:    %8.5f Prob up right:     %8.5f", sum_upleftPerformed, sum_uprightPerformed);
    printf("\n");

    printf("Prob down left:  %8.5f Prob down right:   %8.5f", sum_downleftTested, sum_downrightTested);
    printf("\t         ");
    printf("Prob down left:  %8.5f Prob down right:   %8.5f", sum_downleftDiscarded, sum_downrightDiscarded);
    printf("\t         ");
    printf("Prob down left:  %8.5f Prob down right:   %8.5f", sum_downleftPerformed, sum_downrightPerformed);
    printf("\n");



    printf("Prob straight left:  %8.5f Prob straight up:   %8.5f", sum_onlyleftTested, sum_onlyupTested);
    printf("       ");
    printf("Prob straight left:  %8.5f Prob straight up:   %8.5f", sum_onlyleftDiscarded, sum_onlyupDiscarded);
    printf("      ");
    printf("Prob straight left:  %8.5f Prob straight up:   %8.5f", sum_onlyleftPerformed, sum_onlyupPerformed);
    printf("\n");

    printf("Prob straight right: %8.5f Prob straight down: %8.5f", sum_onlyrightTested, sum_onlydownTested);
    printf("       ");
    printf("Prob straight right: %8.5f Prob straight down: %8.5f", sum_onlyrightDiscarded, sum_onlydownDiscarded);
    printf("      ");
    printf("Prob straight right: %8.5f Prob straight down: %8.5f", sum_onlyrightPerformed, sum_onlydownPerformed);
    printf("\n");

    printf("\n\n---------------------------------------------------------------------------------------------------------------------------------------------------------\n");
    printf("\n\n");


    for (i = 0; i < 2*L+1; i++) {
        if (jumpsPerformed[i] != NULL) delete[] jumpsPerformed[i];
        if (tracerSumPerformed[i] != NULL) delete[] tracerSumPerformed[i];
        if (counterPerformed[i] != NULL) delete[] counterPerformed[i];
        if (jumpsTested[i] != NULL) delete[] jumpsTested[i];
        if (tracerSumTested[i] != NULL) delete[] tracerSumTested[i];
        if (counterTested[i] != NULL) delete[] counterTested[i];
        if (jumpsDiscarded[i] != NULL) delete[] jumpsDiscarded[i];
        if (tracerSumDiscarded[i] != NULL) delete[] tracerSumDiscarded[i];
        if (counterDiscarded[i] != NULL) delete[] counterDiscarded[i];
    }
    delete[] jumpsPerformed;
    delete[] jumpsDiscarded;
    delete[] jumpsTested;

    delete[] tracerSumPerformed;
    delete[] tracerSumDiscarded;
    delete[] tracerSumTested;

    delete[] counterPerformed;
    delete[] counterDiscarded;
    delete[] counterTested;

    delete tracerDiscarded;
    delete tracerPerformed;
    delete tracerTested;

}


void MKcurrent::sample2SiteJumpRates(int steps)
{
    int s,i,j;
    double **jumps;
    int L=3;
    jumps = new double*[2*L+1];
    for (i = 0; i < 2*L+1; i++) {
        jumps[i] = new double[2*L+1];
        for (j = 0; j < 2*L+1; ++j){
            jumps[i][j] = 0;
        }
    }
    for (s=0; s < steps; s++){
        if (dxIs[s] == 0 && dyIs[s] == 0){ // If we did 2 jump
            if (fabs(dxs[s]) <= L && fabs(dys[s]) <= L ){
                jumps[dys[s]+L][dxs[s]+L]++;
            }
        }
    }
    double sum=0, sum_up=0, sum_down=0, sum_left=0, sum_right=0;

    printf("2-site jump probability matrix:\n");

    for (j=L; j>-L-1; j--) {
        for (i=-L; i<L+1; i++){
            jumps[j+L][i+L] /= numberOf2Site;
            sum+= jumps[j+L][i+L];
            if (i<0) sum_left += jumps[j+L][i+L];
            else if (i>0) sum_right += jumps[j+L][i+L];
            if (j>0) sum_up += jumps[j+L][i+L];
            else if (j<0) sum_down += jumps[j+L][i+L];
            printf("%8.5f ", jumps[j+L][i+L]);
        }
        printf("\n");
    }

    printf("Prob 2 sum: %8.5f\n",sum);
    printf("Prob left:  %8.5f Prob up:   %8.5f\n", sum_left, sum_up);
    printf("Prob right: %8.5f Prob down: %8.5f\n", sum_right, sum_down);
    printf("\n");

    for (i = 0; i < 2*L+1; i++) {
        if (jumps[i] != NULL) delete[] jumps[i];
    }
    delete[] jumps;

}

void MKcurrent::comparePerformedWithInitialized(int steps)
{
    int s,i,j,iI,jI;
    int L = maxJL;

    double **jumpsPerformed,**counterPerformed, **tracerSumPerformed;

    double **jumpsTested, **counterTested, **tracerSumTested;
    double *tracerPerformed = (double*)calloc( (2*L+1) * (2*L+1)* (2*L+1)* (2*L+1),sizeof(double));
    jumpsPerformed = new double*[2*L+1];
  counterPerformed = new double*[2*L+1];
tracerSumPerformed = new double*[2*L+1];


jumpsTested = new double*[2*L+1];
counterTested = new double*[2*L+1];
tracerSumTested = new double*[2*L+1];
    for (i = 0; i < 2*L+1; i++) {
            jumpsPerformed[i] = new double[2*L+1];
          counterPerformed[i] = new double[2*L+1];
        tracerSumPerformed[i] = new double[2*L+1];

        jumpsTested[i] = new double[2*L+1];
      counterTested[i] = new double[2*L+1];
    tracerSumTested[i] = new double[2*L+1];
    }

    double *tracerTested = (double*)calloc( (2*L+1) * (2*L+1)* (2*L+1)* (2*L+1),sizeof(double));

    for (s=0; s < steps; s++){
        if (dxIs[s] != 0 || dyIs[s] != 0){ // If we did 3 jump in step s
            if (fabs(dxs[s]) <= L && fabs(dys[s]) <= L ){
                jumpsPerformed[dys[s]+L][dxs[s]+L]++;
                if (fabs(dxIs[s]) <= L && fabs(dyIs[s]) <= L ){
                    counterPerformed[dys[s]+L][dxs[s]+L]++;
                    tracerPerformed[  (dys[s]+L)*(2*L+1)*(2*L+1)*(2*L+1) +
                             (dxs[s]+L)*(2*L+1)*(2*L+1) +
                            (dyIs[s]+L)*(2*L+1)+
                            (dxIs[s]+L)]++;
                }
            }
        }
    }
    // Trace proposed (tested) jumps
    for (s=0; s < testedNumberOf3Site; s++){
            if (testdxIs[s] != 0 || testdyIs[s] != 0){ // If we did 3 jump in step s
                if (fabs(testdxs[s]) <= L && fabs(testdys[s]) <= L ){
                    jumpsTested[testdys[s]+L][testdxs[s]+L]++;
                    if (fabs(testdxIs[s]) <= L && fabs(testdyIs[s]) <= L ){
                        counterTested[testdys[s]+L][testdxs[s]+L]++;
                        tracerTested[  (testdys[s]+L)*(2*L+1)*(2*L+1)*(2*L+1) +
                                 (testdxs[s]+L)*(2*L+1)*(2*L+1) +
                                (testdyIs[s]+L)*(2*L+1)+
                                (testdxIs[s]+L)]++;
                    }
                }
            }
        }


    printf("\n");
    char *st = "FINAL";
    char *stc = "START";

    printf("\nUsed 3-site jump probability matrices\n");

    // Initialized numeric rate, tested monte carlo rate, performed jump rate
    for (i=-L; i<L+1; i++){
            for (j = -L; j < L+1; j++) {
                printf(" Initialized \t \t  \t \t  \t \t \t  \t \t  \t \t \t  \t Performed \n");
                for (jI=L; jI>-L-1; jI--) {
                    for (iI=-L; iI<L+1; iI++){
//                        if (fabs(iI) <= L && fabs(jI) <= L && fabs(i) <= L && fabs(j) <= L)
                        transitionTracer[ (j+L)*(2*L+1)*(2*L+1)*(2*L+1) +
                                          (i+L)*(2*L+1)*(2*L+1) +
                                         (jI+L)*(2*L+1)+
                                         (iI+L)] /= rateMatrix[j+L][i+L]; // Normalize

                        if (i==iI && j==jI){
                            printf("%9s",st);
                        }
                        else if (iI==0 && jI==0) printf("%9s", stc);
                        else printf("%9.5f",transitionTracer[(j+L)*(2*L+1)*(2*L+1)*(2*L+1) +
                                                            (i+L)*(2*L+1)*(2*L+1) +
                                                            (jI+L)*(2*L+1)+
                                                           (iI+L)]);
                    }
                    printf("\t");


                    for (iI=-L; iI<L+1; iI++){
                        if (i+L!=L || j+L!=L){
                            tracerPerformed[(j+L)*(2*L+1)*(2*L+1)*(2*L+1) +
                                                              (i+L)*(2*L+1)*(2*L+1) +
                                                             (jI+L)*(2*L+1)+
                                                             (iI+L)] /= jumpsPerformed[j+L][i+L];
                        }
                        if (i==iI && j==jI){
                            tracerSumPerformed[j+L][i+L] += jumpsPerformed[j+L][i+L]/numberOf3Site;
                            printf("%9s",st);
                        }
                        else if (iI==0 && jI==0) printf("%9s", stc);
                        else printf("%9.5f",tracerPerformed[(j+L)*(2*L+1)*(2*L+1)*(2*L+1) +
                                                            (i+L)*(2*L+1)*(2*L+1) +
                                                           (jI+L)*(2*L+1)+
                                                           (iI+L)]);
                    }
                    printf("\n");
                }

                printf("Jumps going here: %8.6f",rateMatrix[j+L][i+L]/GammaTtot3Sites );
                printf("\t \t \t  \t \t \t \t \t \t   \t \t ");
                printf("Jumps going here: %8.6f",tracerSumPerformed[j+L][i+L] );

                printf("\n\n");
                printf("Tested jumps\n");


                for (jI=L; jI>-L-1; jI--) {
                for (iI=-L; iI<L+1; iI++){
                    if (i+L!=L || j+L!=L){ // To avoid dividing by zero
                        tracerTested[(j+L)*(2*L+1)*(2*L+1)*(2*L+1) +
                                                          (i+L)*(2*L+1)*(2*L+1) +
                                                         (jI+L)*(2*L+1)+
                                                         (iI+L)] /= jumpsTested[j+L][i+L];
                    }
                    if (i==iI && j==jI){
                        tracerSumTested[j+L][i+L] += jumpsTested[j+L][i+L]/testedNumberOf3Site;
                        printf("%9s",st);
                    }
                    else if (iI==0 && jI==0) printf("%9s", stc);
                    else printf("%9.5f",tracerTested[(j+L)*(2*L+1)*(2*L+1)*(2*L+1) +
                                                        (i+L)*(2*L+1)*(2*L+1) +
                                                       (jI+L)*(2*L+1)+
                                                       (iI+L)]);
                }
                printf("\n");
                }
                printf("Jumps going here: %8.6f",tracerSumTested[j+L][i+L] );
                printf("\n-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
                printf("\n\n");

            }
        }

    int dxF, dyF;
    printf("\n----------------------------------------------------------------------------------------------------------------------------------------");
    printf("\n----------------------------------------------------------------------------------------------------------------------------------------");
    printf("\n----------------------------------------------------------------------------------------------------------------------------------------");
    printf("\n\n Initialized jump rates: \t \t \t \t \t \t \t \t \t \t Performed jump rates: \n");

    for (dyF= maxJL; dyF> -maxJL - 1; dyF--){
            for (dxF= -maxJL; dxF< maxJL + 1; dxF++){
                printf("%9.5f", rateMatrix[dyF+maxJL][dxF+maxJL]/GammaTtot3Sites);
            }
            printf("\t");
            for (dxF= -maxJL; dxF< maxJL + 1; dxF++){
                printf("%9.5f", jumpsPerformed[dyF+maxJL][dxF+maxJL]/numberOf3Site);
            }
            printf("\n");
    }
    printf("\nTested jump rates: \n");
    for (dyF= maxJL; dyF> -maxJL - 1; dyF--){
            for (dxF= -maxJL; dxF< maxJL + 1; dxF++){
                printf("%9.5f", jumpsTested[dyF+maxJL][dxF+maxJL]/testedNumberOf3Site);
            }
            printf("\n");
    }
    printf("\n----------------------------------------------------------------------------------------------------------------------------------------");
    printf("\n----------------------------------------------------------------------------------------------------------------------------------------");
    printf("\n----------------------------------------------------------------------------------------------------------------------------------------");

    for (i = 0; i < 2*L+1; i++) {
        if (jumpsPerformed[i] != NULL) delete[] jumpsPerformed[i];
        if (tracerSumPerformed[i] != NULL) delete[] tracerSumPerformed[i];
        if (counterPerformed[i] != NULL) delete[] counterPerformed[i];
        if (jumpsTested[i] != NULL) delete[] jumpsTested[i];
        if (tracerSumTested[i] != NULL) delete[] tracerSumTested[i];
        if (counterTested[i] != NULL) delete[] counterTested[i];
    }
    delete[] jumpsPerformed;
    delete[] jumpsTested;

    delete[] tracerSumPerformed;
    delete[] tracerSumTested;

    delete[] counterPerformed;
    delete[] counterTested;

    delete tracerPerformed;
    delete tracerTested;
}

void MKcurrent::sampleProposed3SiteJumps(int steps)
{
    int s, i,j,iI,jI;
    double **jumps,**counter, **tracerSum;

    int L=3;
    // Tracer tracks what intermediate jump went to what final site
    jumps = new double*[2*L+1];
    counter = new double*[2*L+1];
    tracerSum = new double*[2*L+1];

    for (i = 0; i < 2*L+1; i++) {
        jumps[i] = new double[2*L+1];
        counter[i] = new double[2*L+1];
        tracerSum[i] = new double[2*L+1];
    }

    double *tracer = (double*)calloc( (2*L+1) * (2*L+1)* (2*L+1)* (2*L+1),sizeof(double));

    for (s=0; s < test3JumpCounter; s++){
        if (testdxIs[s] != 0 || testdyIs[s] != 0){ // If we did 3 jump in step s
            for (i=-L; i<L+1; i++){ // Potential final site in jump s
                for (j = -L; j < L+1; j++) {
                    if (testdxs[s] == i && testdys[s] == j){
                        jumps[j+L][i+L]++;
                        for (iI=-L; iI<L+1; iI++){ // Potential intermediate site
                            for (jI=-L; jI<L+1; jI++) {
                                if (testdxIs[s] == iI && testdyIs[s] == jI){
                                    counter[j+L][i+L]++;
                                    tracer[  (j+L)*(2*L+1)*(2*L+1)*(2*L+1) +
                                             (i+L)*(2*L+1)*(2*L+1) +
                                            (jI+L)*(2*L+1)+
                                            (iI+L)]++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    fflush(stdout);
    printf("\n");
    char *st = "FINAL";
    char *stc = "START";

    printf("\nProposed 3-site jump probability matrices\n");

    for (i=-L; i<L+1; i++){
        for (j = -L; j < L+1; j++) {
            for (jI=L; jI>-L-1; jI--) {
                for (iI=-L; iI<L+1; iI++){
                    if (i+L!=L || j+L!=L) tracer[(j+L)*(2*L+1)*(2*L+1)*(2*L+1) +
                                                 (i+L)*(2*L+1)*(2*L+1) +
                                                (jI+L)*(2*L+1)+
                                                (iI+L)] /= jumps[j+L][i+L];
                    if (i==iI && j==jI){
                        tracerSum[j+L][i+L] += jumps[j+L][i+L]/test3JumpCounter;
                        printf("%8s",st);
                    }
                    else if (iI==0 && jI==0) printf("%8s", stc);
                    else printf("%8.5f",tracer[(j+L)*(2*L+1)*(2*L+1)*(2*L+1) +
                                               (i+L)*(2*L+1)*(2*L+1) +
                                              (jI+L)*(2*L+1)+
                                              (iI+L)]);
                }
                printf("\n");
            }
            printf("Jumps tested included:   %8.5f\n",  1-(jumps[j+L][i+L]-counter[j+L][i+L])/jumps[j+L][i+L]  );
            printf("Jumps tested going here: %8.5f\n\n",tracerSum[j+L][i+L] );

        }
    }

    printf("\n\n");
    double sum=0, sum_up=0, sum_down=0, sum_left=0, sum_right=0;
    printf("3-site jump test probability matrix:\n");

    for (j=L; j>-L-1; j--) {
        for (i=-L; i<L+1; i++){
            jumps[j+L][i+L] /= test3JumpCounter;
            sum += jumps[j+L][i+L];
            if (i==0 && j==0) printf("%8s", stc);
            else printf("%8.5f ", jumps[j+L][i+L]);
            if (j>0) sum_up += jumps[j+L][i+L];
            else if (j<0) sum_down += jumps[j+L][i+L];
            if (i<0) sum_left += jumps[j+L][i+L];
            else if (i>0) sum_right += jumps[j+L][i+L];
        }
        printf("\n");
    }

    printf("Prob 3 sum: %8.5f\n",sum);
    printf("Prob left:  %8.5f Prob up:   %8.5f\n", sum_left, sum_up);
    printf("Prob right: %8.5f Prob down: %8.5f\n", sum_right, sum_down);
    printf("\n");
    delete tracer;
}

void MKcurrent::writeAreaHistogramToFile(string filename)
{

    int s, cumdx, cumdy;
    FILE* f;
    string st;

    f=FileCreate(filename);
    double area;

    st = "Area \t probability \t bin #\n";
    FileWrite(f,st.c_str(),st.length());
    double sum = 0;
    for(s = 0; s < bins; s++){
        areaHistogram[s] /= numberOf3Site;
        sum += areaHistogram[s];

//        if (s > bins/2)area = (-bins/2 + s)*bucketSize - bins/4;
        if (s > bins/2)area = -(-bins/2 + s)*bucketSize;
        else area = s*bucketSize;
        st = FloatToStr(area) + "\t" + DoubleToStr(areaHistogram[s])  + "\t" + IntToStr(s)+ "\n";
        FileWrite(f,st.c_str(),st.length());
    }

    FileClose(f);
}

void MKcurrent::printWeightedAreaTimeData(){

    int s;
    FILE* f;
    string st;

    f = FileCreate("weightedAreaData.dat");
    double cumulativeWeightedArea = 0;
    int cumdy = 0;
    for (s = 0; s<numberOf3Site; s++){
        cumulativeWeightedArea += meanWeightedAreaPer[s];
        cumdy += dyPer[s];

        st = IntToStr(s) + "\t" +  FloatToStr(cumulativeWeightedArea)  + "\t"
                + DoubleToStr(cumdy)+ "\t" + IntToStr(dy3) + "\n";
        FileWrite(f,st.c_str(),st.length());
    }

    FileClose(f);


}
