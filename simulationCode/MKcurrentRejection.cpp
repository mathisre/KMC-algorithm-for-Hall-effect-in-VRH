//---------------------------------------------------------------------------

#include "fileutils.h"
#include "stringutils.h"
#include "CGanalysis.h"
#include "systemclass.h"
#include "treeutils.h"
#include "MKcurrentRejection.h"




#include <fstream>
MKcurrentRejection::MKcurrentRejection(int steps)
{
  from = new int[steps];
  to = new int[steps];
  energy = new double[steps];
  MCsteps = new int[steps];
  dxs = new int[steps];
  dys = new int[steps];
  dxsTri = new double[steps];
  dysTri = new double[steps];
  ts = new double[steps];
  des = new double[steps];
}

void MKcurrentRejection::initMaps(){

    movement   = new double*[L];
    xmovement  = new double*[L];
    ymovement  = new double*[L];

    for (int x = 0; x<L; x++){
        movement[x] = new double[L];
        xmovement[x] = new double[L];
        ymovement[x] = new double[L];
    }
}



MKcurrentRejection::~MKcurrentRejection()
{

  if (from!=NULL) delete[] from;
  if (to!=NULL) delete[] to;
  if (jl!=NULL) delete[] jl;
  if (GammaT!=NULL) delete[] GammaT;
  if (energy!=NULL) delete[] energy;
  if (MCsteps!=NULL) delete[] MCsteps;
  if (dxs!=NULL) delete[] dxs;
  if (dys!=NULL) delete[] dys;
  if (dxsTri!=NULL) delete[] dxsTri;
  if (dysTri!=NULL) delete[] dysTri;
  if (ts!=NULL) delete[] ts;
  if (des!=NULL) delete[] des;
  if (dxInter!=NULL) delete[] dxInter;
  if (dyInter!=NULL) delete[] dyInter;
  if (dxFinal!=NULL) delete[] dxFinal;
  if (dyFinal!=NULL) delete[] dyFinal;
  if (GammaT3Sites!=NULL) delete[] GammaT3Sites;
}


void MKcurrentRejection::setE(double Ex1, double Ey1, double Ez1)
{
  Ex=Ex1; Ey=Ey1; Ez=Ez1;
}

void MKcurrentRejection::setH(double Hx1, double Hy1, double Hz1)
{
  Hx=Hx1;
  Hy=Hy1;
  Hz=Hz1;
}

void MKcurrentRejection::setE0(double Ex1, double Ey1, double Ez1)
{
  E0x=Ex1; E0y=Ey1; E0z=Ez1;
}

void MKcurrentRejection::setOmega(double omega1)
{
  omega=omega1;
}

void MKcurrentRejection::setWritelines(int writelines1)
{
  writelines=writelines1;
}

void MKcurrentRejection::setRatefun(int ratefun1)
{
  ratefun = ratefun1;

}

void MKcurrentRejection::setES(ESystem *e)
{
  es = e;
}


void MKcurrentRejection::init(int D1, int L1, int N1, int maxJL1, double a1, double beta1, double maxProbability1)
{ //  (D, L, N, p->maxJL, p->loc, 1.0/p->temp,p->maxProbability)
  int dx, dy, n;
  RanGen = new CRandomMersenne(1); //seed??

  D=D1; L=L1; N=N1; A=a1; beta = beta1;
  maxJL = maxJL1; JL2 = 2*maxJL+1;

  maxProbability=maxProbability1/beta;
  maxAccept = 2*exp(-A);
//  writelines = 10000;


  //printf("L: %d N: %d maxJL: %d\n", L, N, maxJL);
  if (D==2)
    {
      Nmem = JL2*JL2; // Number of sites we can jump to
      int Nmem2 = (JL2+2)*(JL2+2);
      jl     = new double[Nmem2];
      GammaT = new double[Nmem];
      dxIndex = new int[Nmem];
      dyIndex = new int[Nmem];
      n = 0;
      for (dx = -maxJL-1; dx < maxJL + 2; dx++)
      for (dy = -maxJL-1; dy < maxJL + 2; dy++)
      {
        jl[(dx+maxJL+1) + JL2*(dy+maxJL+1)] = sqrt(dx*dx + dy*dy);
      }
      jl[maxJL + 1 + JL2*(maxJL+1)] = 0.0;

      for (dx = -maxJL; dx < maxJL + 1; dx++)
      for (dy = -maxJL; dy < maxJL + 1; dy++)
      {
        dxIndex[(dx+maxJL) + JL2*(dy+maxJL)] = dx;
        dyIndex[(dx+maxJL) + JL2*(dy+maxJL)] = dy;
      }
    }

}




void MKcurrentRejection::setMTseed(int seed)
{
  RanGen->RandomInit(seed);
}


void inline MKcurrentRejection::getjump(int &i, int &j, int &dx, int &dy)
{
  int x, y, x2, y2;

  i = RanGen->IRandom(0,N-1);
  while (es->getocci(i) == 0){ // i is occupied
    i = RanGen->IRandom(0,N-1);
  }

  x = i%L;
  y = i/L;

  j = RanGen->IRandom(0,Nmem);

  dx = dxIndex[j];
  dy = dyIndex[j];

  x2 = x+dx; if (x2 >= L) x2 -= L; else if (x2 < 0) x2 += L;
  y2 = y+dy; if (y2 >= L) y2 -= L; else if (y2 < 0) y2 += L;

  j = x2 + y2 * L; // end position

  while (es->getocci(j) == 1){ // j is free
    j = RanGen->IRandom(0, Nmem);
    dx = dxIndex[j];
    dy = dyIndex[j];

    x2 = x+dx; if (x2 >= L) x2 -= L; else if (x2 < 0) x2 += L;
    y2 = y+dy; if (y2 >= L) y2 -= L; else if (y2 < 0) y2 += L;

    j = x2 + y2 * L; // end position
  }

}


bool inline MKcurrentRejection::testjump(int i, int j, int dx, int dy)
{
  double r2, gamma;
  double r_ij, dEij;

  r_ij =  jl[(dx+maxJL+1) + JL2*(dy+maxJL+1)];
  if (r_ij > maxJL) return false;

  dEij = es->hoppEdiffij(i,j) + dx*Ex + dy*Ey;
  (dEij < 0 ) ? dEij = 1 : dEij = exp(-beta*dEij);
  gamma = exp(-A*r_ij)*dEij;

  if (Hz != 0) gamma += threeSiteCorrectionRates(i, j, dx, dy, dEij, r_ij);

//  double testsum = gamma - (exp(-A*r_ij)*dEij);
//  if (dx == -4 && dy == 0)printf("dx: %2d, dy: %2d, full rate: %7.4e, 2_rate %7.4e, correction %7.4f\n", dx, dy, gamma, exp(-A*r_ij)*dEij, testsum / (exp(-A*r_ij)*dEij));
//  if (dx == -1 && dy == 0 && dEij == 1 ) printf("g/g_max: %6.4e\n", gamma/maxAccept);

  if (gamma >= maxAccept) printf("Rate larger than max accept!!! gamma: %.4f, maxAccept: %.4f\n", gamma, maxAccept);
  r2 = maxAccept*es->ran2(0);

  return (r2 < gamma);

}

double MKcurrentRejection::threeSiteCorrectionRates(int i, int j, int dx, int dy, double dEij, double r_ij){
    double dEik, dEjk, dEkj, dEki, area, overlap3;
    double acceptance3, r_ik, r_jk, rinv_ij, rinv_ik, rinv_jk;
    int x, y, dxI, dyI, k, xI, yI;
    int startDx, startDy, endDx, endDy;
    double gammaThree = 0;

    x = i%L; y = i/L;

    rinv_ij = 1 / r_ij;;
    if (rinv_ij < es->rmaxi) rinv_ij = 0;

    // Optimization to only consider sites  k along the path from i to j
    startDx = 0;  startDy = 0;
    endDx   = dx; endDy   = dy;

    if (dx < 0){
        startDx = dx;
        endDx = 0;
    }
    if (dy < 0){
        startDy = dy;
        endDy = 0;
    }
    if (dx == 0){
        startDx = -1;
        endDx = 1;
    }
    if (dy == 0){
        startDy = -1;
        endDy = 1;
    }


    for (dxI = startDx; dxI <= endDx; dxI+=1){
        for (dyI = startDy; dyI <= endDy; dyI+=1){


            area = 0.5*(dxI*dy - dx*dyI);
            if (area != 0){ // makes it faster

                r_ik =  jl[(dxI+maxJL+1) + JL2*(dyI+maxJL+1)];
                r_jk = sqrt((dx-dxI)*(dx-dxI) + (dy-dyI)*(dy-dyI));
                overlap3 = exp(-0.5*A*(r_ij + r_ik + r_jk));

                rinv_ik = 1 / r_ik;
                rinv_jk = 1 / r_jk;
                if (rinv_ik < es->rmaxi) rinv_ik = 0;
                if (rinv_jk < es->rmaxi) rinv_jk = 0;

                xI = x+dxI; if (xI >= L) xI -= L; else if (xI < 0) xI += L;
                yI = y+dyI; if (yI >= L) yI -= L; else if (yI < 0) yI += L;

                k = xI + yI * L; // intermediate site

                if (es->getocci(k) == 0){

                    dEik = es->hoppEdiffij(i,k) + dxI*Ex + dyI*Ey;
                    dEjk = es->hoppEdiffij(j,k) + (dxI-dx)*Ex + (dyI-dy)*Ey + rinv_ij + rinv_jk - rinv_ik;
                    dEkj = es->hoppEdiffij(k,j) + (dx-dxI)*Ex + (dy-dyI)*Ey + rinv_ik + rinv_jk - rinv_ij;

                    (dEik < 0) ? dEik = 1 : dEik = exp(-beta*dEik);
                    (dEjk < 0) ? dEjk = 1 : dEjk = exp(-beta*dEjk);
                    (dEkj < 0) ? dEkj = 1 : dEkj = exp(-beta*dEkj);

                    acceptance3 = dEik*dEkj + dEij*dEik + dEij*dEjk;
                    gammaThree += overlap3*(Hz*area)*acceptance3;

                }
                else {

                    dEik = es->hoppEdiffij(i,k) + dxI*Ex + dyI*Ey + rinv_jk + rinv_ik - rinv_ij;
                    dEki = es->hoppEdiffij(k,i) - dxI*Ex - dyI*Ey - rinv_jk + rinv_ik + rinv_ij;
                    dEkj = es->hoppEdiffij(k,j) + (dx-dxI)*Ex + (dy-dyI)*Ey;

                    (dEik < 0) ? dEik = 1 : dEik = exp(-beta*dEik);
                    (dEki < 0) ? dEki = 1 : dEki = exp(-beta*dEki);
                    (dEkj < 0) ? dEkj = 1 : dEkj = exp(-beta*dEkj);

                    acceptance3 = dEik*dEkj + dEij*dEki + dEij*dEkj;
                    gammaThree -= overlap3*(Hz*area)*acceptance3;
                }
            }
        }
    }
    return gammaThree;
}


void MKcurrentRejection::runCurrent(int steps,double &E, double &t)
{

  int s, MCs, i=0, j=0, dx=0, dy=0, jumpNumber = 0;
  int k, dxI=0, dyI=0;
  double dE,dt;
  bool jump;

  tMC = 1/(L*L*(Nmem-1)*es->occnum*(1-es->occnum)*maxAccept);

  meanArea = 0;

  meanJumpLength = 0;
  meanInterJumpLength = 0;
  meandx = 0;
  meandy = 0;
  meandxI = 0;
  meandyI = 0;
  testedNumberOf2Site = 0;
  testedNumberOf3Site = 0;
  numberOf2Site = 0;
  numberOf3Site = 0;

  meanSomething = 0;

  for (s = 0; s < steps; s++)
  {
      MCs = 0;
      jump = false;
      while (!jump)
      {
              MCs++;

                jumpNumber = 11; dxI = 0; dyI = 0;
                testedNumberOf2Site++;
                getjump(i,j,dx,dy);
                jump=testjump(i,j,dx,dy);

      }
      dE = es->hopp(i,-1,0,j,0,0);

      dt = MCs*tMC;
      E += dE;
      t += dt;

      from[s] = i;
      to[s]   = j;

      ts[s] = t;
      energy[s] = E;
      des[s] = dE;
      dxs[s] = dx;
      dys[s] = dy;
      MCsteps[s]=MCs;

      if (s%100000==0) {
          printf("%6d E=%le MCs=%d\n",s,E,MCs);
      }

    }
  printf("\n-------Finishing timesteps--------\n\n");

//  printf("2 site %d, 3 site %d\n",numberOf2Site,numberOf3Site);
//  printf("Tested 2 site %d, tested 3 site %d\n",testedNumberOf2Site,testedNumberOf3Site);

}

void MKcurrentRejection::runCurrentMeasure(int steps,double &E, double &t)
{

  int s, MCs, i=0, j=0, dx=0, dy=0;
  double dE,dt;
  bool jump;
  tMC = 1/(L*L*(Nmem-1)*es->occnum*(1-es->occnum)*maxAccept);
  printf("t: %.7f\n", tMC);


  meanJumpLength = 0;
  meandx = 0;
  meandy = 0;
  testedNumberOf2Site = 0;
  numberOf2Site = 0;



  double meanDiscdx = 0, meanDiscdy = 0;
  double meanAbsdx = 0, meanAbsdy = 0;

  double meanMCs = 0;
  double actualmeanMCs = 0;


  for (s = 0; s < steps; s++)
  {
      MCs = 0;
      jump = false;
      while (!jump)
      {

          MCs++;

          testedNumberOf2Site++;
          getjump(i,j,dx,dy);
          jump=testjump(i,j,dx,dy);

          meanAbsdx +=dx; meanAbsdy +=dy;
          if (jump == true) {
              numberOf2Site++;
          }
          else meanDiscdx += dx; meanDiscdy += dy;

      }

      dE = es->hopp(i,-1,0,j,0,0);

      dt = MCs*tMC;
      E += dE;
      t += dt;
      meanMCs += MCs;


      meanJumpLength += jl[(dx+maxJL+1) + JL2*(dy+maxJL+1)];
      meandx  += dx;
      meandy  += dy;

      actualmeanMCs += MCs;

      from[s] = i;
      to[s]   = j;

      ts[s] = t;
      energy[s] = E;
      des[s] = dE;
      dxs[s] = dx;
      dys[s] = dy;
      MCsteps[s]=MCs;

      if (s%10000==0) {
          printf("%6d E=%le t=%.4e MCs=%.3f\n",s,E,t,actualmeanMCs/10000);
          actualmeanMCs =0;
      }

     // printf("Step: %d: %d %d dE: %le I:%le MCs: %d\n", s, i, j,energy[s],dx,MCsteps[s]);
    }
  printf("\n-------Finishing timesteps--------\n\n");


  meanJumpLength /= steps;
  meandx /= steps;
  meandy /= steps;

  printf("\nProposed  dx: %8.5f, dy:  %8.5f\n", meanAbsdx /testedNumberOf2Site, meanAbsdy /testedNumberOf2Site);
  printf("Discarded dx: %8.5f, dy:  %8.5f\n", meanDiscdx/testedNumberOf2Site, meanDiscdy/testedNumberOf2Site);
  printf("Mean dx:      %8.5f, dy:  %8.5f\n", meandx, meandy);
  printf("\nMean jump length:                    %8.5f\n", meanJumpLength);
  printf("Acceptance ratio: %.5f\n", double(numberOf2Site)/testedNumberOf2Site);

}

void MKcurrentRejection::jumpsToFileSmall(string filename, int steps)
{
  int s, cumdx, cumdy;
  FILE* f;
  string st;

  cumdx=0;
  cumdy=0;
  f=FileCreate(filename);

  st="t  \t \t E \t \t from \t to \t dx \t dy \t MCs \t dE \t tMC="+DoubleToStr(tMC)+ "\t <jl>=" + DoubleToStr(meanJumpLength) +"\t <dx>=" + DoubleToStr(meandx) +"\t <dy>=" + DoubleToStr(meandy) + "\t <jl_I>=" +DoubleToStr(meanInterJumpLength) +"\t <dx_I>=" + DoubleToStr(meandxI) +"\t <dy_I>=" + DoubleToStr(meandyI) + "\t Acc 2=" + DoubleToStr(double(numberOf2Site)/testedNumberOf2Site)+ "\t Acc 3=" + DoubleToStr(double(numberOf3Site)/testedNumberOf3Site) + "\n";
  FileWrite(f,st.c_str(),st.length());

  for(s = 0; s < steps; s++){
      cumdx += dxs[s];
      cumdy += dys[s];
      if (s%writelines == 0){
        st = FloatToStr(ts[s])+"\t"+DoubleToStr(energy[s])+"\t"+IntToStr(from[s])+"\t"+IntToStr(to[s])+"\t"+IntToStr(cumdx)+"\t"+IntToStr(cumdy)+"\t"+IntToStr(MCsteps[s])+"\t"+FloatToStr(des[s])+"\n";
        FileWrite(f,st.c_str(),st.length());
      }
    }
  printf("Rough estimate for conductivity_x: %8.5f\n", cumdx/ts[steps-1]);
  printf("Rough estimate for conductivity_y: %8.5f\n", cumdy/ts[steps-1]);

  FileClose(f);
}




void MKcurrentRejection::runCurrentRandomLattice(int steps,double &E, double &t)
{

  int s, MCs, i=0, j=0,  neighborNumber=0;
  double dx=0, dy=0;
  double dE,dt;
  bool jump;

  tMC = 1/(L*L*(averageNmem-1)*es->occnum*(1-es->occnum)*maxAccept);
  printf("t: %.7f\n", tMC);



  dxsReal = new double[steps];
  dysReal = new double[steps];

  double meanMCs = 0;
  double actualmeanMCs = 0;


  for (s = 0; s < steps; s++)
  {
      MCs = 0;
      jump = false;
      while (!jump)
      {

          MCs++;;
          getjumpRandomLattice(i,j,neighborNumber);
          jump=testjumpRandomLattice(i,j,neighborNumber);
      }
      dE = es->hoppRandomPositions(i,-1,0,j,0,0);

      dt = MCs*tMC;
      E += dE;
      t += dt;
      meanMCs += MCs;

      meanJumpLength += es->getDistanceMatrix(i,j);

      dx = positiondx[i][neighborNumber];
      dy = positiondy[i][neighborNumber];

      actualmeanMCs += MCs;

      from[s] = i;
      to[s]   = j;

      ts[s] = t;
      energy[s] = E;
      des[s] = dE;
      dxsReal[s] =  dx;
      dysReal[s] =  dy;
      MCsteps[s]=MCs;

      if (s%100==0) {
          printf("%6d E=%le t=%.4e MCs=%.3f\n",s,E,t,actualmeanMCs/100);
          actualmeanMCs =0;
      }

     // printf("Step: %d: %d %d dE: %le I:%le MCs: %d\n", s, i, j,energy[s],dx,MCsteps[s]);
    }
  printf("\n-------Finishing timesteps--------\n\n");


}






void MKcurrentRejection::initRandomLattice(int D1, int L1, int N1, int maxJL1, double a1, double beta1, double maxProbability1, bool CintOn)
{ //  (D, L, N, p->maxJL, p->loc, 1.0/p->temp,p->maxProbability)


    int  n, i, j;
    double dx, dy;
    RanGen = new CRandomMersenne(1);
    Cint = CintOn;




    D=D1; L=L1; N=N1; A=a1; beta = beta1;
    maxJL = maxJL1; JL2 = 2*maxJL+1;

    maxProbability=maxProbability1/beta;
    maxAccept = 1.8*2*6*exp(-A);


  //  Nmem = JL2*JL2; // Number of sites we can jump to
  //  jl            = new double[Nmem];

    sitePositions   = new double*[N];
    randShift = 0.5;

    finalSite    = new int*[N];
    NmemArray = new int[N];

    indexdx = new int*[N];
    indexdy = new int*[N];

    positiondx = new double*[N];
    positiondy = new double*[N];

    neighborList = new int*[N];


    // Initializing site positions
    for (i=0; i<N; i++){
      sitePositions[i] = new double[2];
//      sitePositions[i][0] = i%L + (es->ran2(0)-0.5)*randShift;
//      sitePositions[i][1] = i/L + (es->ran2(0)-0.5)*randShift;

      sitePositions[i][0] = es->ran2(0)*L;
      sitePositions[i][1] = es->ran2(0)*L;
    }


   averageNmem = 0;

double distance = 0;


    // Setting up arrays to be used
    distanceMatrix = new double*[N];
    for (i=0; i<N; i++){
      distanceMatrix[i] = new double[i];
      NmemArray[i] = 0;
      for (j=0; j<N; j++){
          dx = sitePositions[j][0] - sitePositions[i][0];
          dy = sitePositions[j][1] - sitePositions[i][1];

          if      (dx > L/2) dx -= L;
          else if (dx < -L/2) dx += L;

          if      (dy > L/2) dy -= L;
          else if (dy < -L/2) dy += L;


          distance = sqrt(dx*dx + dy*dy);
          if (distance < maxJL) NmemArray[i]++;
          if (j < i) distanceMatrix[i][j] = distance;

      }
      averageNmem += NmemArray[i];
      neighborList[i] = new int[NmemArray[i]];
      indexdx[i] = new int[NmemArray[i]];
      indexdy[i] = new int[NmemArray[i]];
      positiondx[i] = new double[NmemArray[i]];
      positiondy[i] = new double[NmemArray[i]];
    }

    es->setDistanceMatrix(distanceMatrix, N);

    averageNmem /= N;

    double distanceFromMatrix;
    for (i=0; i<N; i++){
        n = 0;
        for (j=0; j<N; j++){

             distance = es->getDistanceMatrix(i,j);
             (j < i) ? distanceFromMatrix = distanceMatrix[i][j] : distanceFromMatrix = distanceMatrix[j][i];
//             printf("distance es: %.3f, distanceMatrix: %.3f\n", es->getDistanceMatrix(i,j), distanceFromMatrix);

             if (distance < maxJL) {
                 neighborList[i][n] = j;
//                 indexdx[i][n] = dxSiteIndex;
//                 indexdy[i][n] = dySiteIndex;


                 dx = sitePositions[j][0] - sitePositions[i][0];
                 dy = sitePositions[j][1] - sitePositions[i][1];


                 dx = sitePositions[j][0] - sitePositions[i][0];
                 dy = sitePositions[j][1] - sitePositions[i][1];

                 if      (dx > L/2) dx -= L;
                 else if (dx < -L/2) dx += L;

                 if      (dy > L/2) dy -= L;
                 else if (dy < -L/2) dy += L;


                 positiondx[i][n] = dx;
                 positiondy[i][n] = dy;

                 n++;
             }


        }
    }

}


void inline MKcurrentRejection::getjumpRandomLattice(int &i, int &j, int &neighborNumber)
{
  int x, y;

  i = RanGen->IRandom(0,N-1);
  while (es->getocci(i) == 0){ // i is occupied
    i = RanGen->IRandom(0,N-1);
  }

  x = i%L;
  y = i/L;

  neighborNumber = RanGen->IRandom(0,NmemArray[i]);
  j = neighborList[i][neighborNumber];

  while (es->getocci(j) == 1){ // j is free
    neighborNumber = RanGen->IRandom(0,NmemArray[i]);
    j = neighborList[i][neighborNumber];

  }

}

#include <stdio.h>

bool inline MKcurrentRejection::testjumpRandomLattice(int i, int j, int neighborNumber)
{
  double r2, gamma, acceptance3;
  double r_ij, r_ik, r_jk, dEij, dEik, dEjk, dEkj, dEki, area, overlap3;
  double rinv_ij, rinv_ik, rinv_jk;
  double physicaldx, physicaldy;
  int indexdxToFinal, indexdyToFinal;
  double box_x_max, box_x_min, box_y_max, box_y_min;
  int neighbor;
  double dxNeighbor, dyNeighbor;

  physicaldx = positiondx[i][neighborNumber];
  physicaldy = positiondy[i][neighborNumber];

  indexdxToFinal = indexdx[i][neighborNumber];
  indexdyToFinal = indexdy[i][neighborNumber];

  r_ij =  es->getDistanceMatrix(i,j);

  dEij = es->hoppEdiffRandomPositions(i,j) + physicaldx*Ex + physicaldy*Ey;
  if (dEij < 0) dEij = 1;
  else dEij = exp(-beta*dEij);
  gamma = exp(-A*r_ij)*dEij;

  rinv_ij = 1 / r_ij;;
  if (!Cint) rinv_ij = 0;

  // BBbox optimization
  makeBBbox(physicaldx, physicaldy, box_x_max, box_x_min, box_y_max, box_y_min);

//  printf("dx: %.3f, dy: %.3f, BBbox: top: %.3f, bot: %.3f, left: %.3f, right: %.3f\n", physicaldx, physicaldy, box_y_max, box_y_min, box_x_min, box_x_max);


  int neighborSite;
  if (Hz != 0) {
      for (neighbor = 0; neighbor<NmemArray[i]; neighbor++){
          dxNeighbor = positiondx[i][neighbor];
          dyNeighbor = positiondy[i][neighbor];
          area = 0.5*(dxNeighbor*physicaldy - physicaldx*dyNeighbor);

          neighborSite = neighborList[i][neighbor];
          if ( neighborSite != j &&  isInsideBox(dxNeighbor, dyNeighbor,box_x_max,box_x_min,box_y_max,box_y_min) ){
              r_ik =  es->getDistanceMatrix(i,neighborSite);
              r_jk =  es->getDistanceMatrix(j,neighborSite);


              overlap3 = exp(-0.5*A*(r_ij + r_ik + r_jk));

              rinv_ik = 1 / r_ik;
              rinv_jk = 1 / r_jk;
              if ( !Cint ) rinv_ik = 0;
              if ( !Cint ) rinv_jk = 0;

              if (es->getocci(neighbor) == 0){

                  dEik = es->hoppEdiffRandomPositions(i,neighborSite) + dxNeighbor*Ex + dyNeighbor*Ey;
                  dEjk = es->hoppEdiffRandomPositions(j,neighborSite) + (physicaldx-dxNeighbor)*Ex + (physicaldy-dyNeighbor)*Ey + rinv_ij + rinv_jk - rinv_ik;
                  dEkj = es->hoppEdiffRandomPositions(neighborSite,j) + (dxNeighbor-physicaldx)*Ex + (dxNeighbor-physicaldy)*Ey + rinv_ik + rinv_jk - rinv_ij;

                  (dEik < 0) ? dEik = 1 : dEik = exp(-beta*dEik);
                  (dEjk < 0) ? dEjk = 1 : dEjk = exp(-beta*dEjk);
                  (dEkj < 0) ? dEkj = 1 : dEkj = exp(-beta*dEkj);

                  acceptance3 = dEik*dEkj + dEij*dEik + dEij*dEjk;
                  gamma += overlap3*(Hz*area)*acceptance3;
              }
              else {

                  dEik = es->hoppEdiffRandomPositions(i,neighborSite) + dxNeighbor*Ex + dyNeighbor*Ey + rinv_jk + rinv_ik - rinv_ij;
                  dEki = es->hoppEdiffRandomPositions(neighborSite,i) - dxNeighbor*Ex - dyNeighbor*Ey - rinv_jk + rinv_ik + rinv_ij;
                  dEkj = es->hoppEdiffRandomPositions(neighborSite,j) + (physicaldy-dxNeighbor)*Ex + (physicaldy-dyNeighbor)*Ey;

                  (dEik < 0) ? dEik = 1 : dEik = exp(-beta*dEik);
                  (dEki < 0) ? dEki = 1 : dEki = exp(-beta*dEki);
                  (dEkj < 0) ? dEkj = 1 : dEkj = exp(-beta*dEkj);

                  acceptance3 = dEik*dEkj + dEij*dEki + dEij*dEkj;
                  gamma -= overlap3*(Hz*area)*acceptance3;
              }
          }
      }
  }



//  double testsum = gamma - exp(-A*r_ij)*dEij;
//  printf("total rate: %.4e, direct rate: %.4e\n", gamma, exp(-A*r_ij)*dEij);
//  printf("dx: %.3f, dy: %.3f, r_ij: %.4f, 2_rate %7.4e, correction %7.4e\n", physicaldx, physicaldy, r_ij, exp(-A*r_ij)*dEij, testsum / (exp(-A*r_ij)*dEij));
//  if (dx == -1 && dy == 0 && dEij == 1 ) printf("g/g_max: %6.4e\n", gamma/maxAccept);

//  printf("\n----------------------------------------------------------------\n");

//  printf("Ratio of used neighbors: %.3f\n", double(insideBox)/NmemArray[i]);


  if (gamma >= maxAccept) printf("Rate larger than max accept!!! gamma: %.4f, maxAccept: %.4f\n", gamma, maxAccept);
  r2 = maxAccept*es->ran2(0);
  return (r2 < gamma);


}

void MKcurrentRejection::makeBBbox(double dx, double dy, double &box_x_max, double &box_x_min, double &box_y_max, double &box_y_min){

    if (dx < 0){
        box_x_min = dx;
        box_x_max = 0;
    }
    if (dx > 0){
        box_x_max = dx;
        box_x_min = 0;
    }

    if (dy < 0){
        box_y_min = dy;
        box_y_max = 0;
    }
    if (dy > 0){
        box_y_max = dy;
        box_y_min = 0;
    }

    if (fabs(dx) < 1){
        box_x_max += 1;
        box_x_min -= 1;
    }

    if (fabs(dy) < 1){
        box_y_max += 1;
        box_y_min -= 1;
    }
}

bool MKcurrentRejection::isInsideBox(double dx, double dy, double box_x_max, double box_x_min, double box_y_max, double box_y_min){
    if (dx > box_x_max) return false;
    if (dx < box_x_min) return false;
    if (dy > box_y_max) return false;
    if (dy < box_y_min) return false;
    return true;
}


void MKcurrentRejection::jumpsToFileSmallRandomLattice(string filename, int steps)
{
  int s;
  double cumdx, cumdy;
  FILE* f;
  string st;

  cumdx=0;
  cumdy=0;
  f=FileCreate(filename);

  st="t  \t \t E \t \t from \t to \t dx \t dy \t MCs \t dE \t tMC="+DoubleToStr(tMC)+ "\t <jl>=" + DoubleToStr(meanJumpLength) +"\t <dx>=" + DoubleToStr(meandx) +"\t <dy>=" + DoubleToStr(meandy) + "\t <jl_I>=" +DoubleToStr(meanInterJumpLength) +"\t <dx_I>=" + DoubleToStr(meandxI) +"\t <dy_I>=" + DoubleToStr(meandyI) + "\t Acc 2=" + DoubleToStr(double(numberOf2Site)/testedNumberOf2Site)+ "\t Acc 3=" + DoubleToStr(double(numberOf3Site)/testedNumberOf3Site) + "\n";
  FileWrite(f,st.c_str(),st.length());

  for(s = 0; s < steps; s++){
      cumdx += dxsReal[s];
      cumdy += dysReal[s];
      if (s%writelines == 0){
        st = FloatToStr(ts[s])+"\t"+DoubleToStr(energy[s])+"\t"+IntToStr(from[s])+"\t"+IntToStr(to[s])+"\t"+DoubleToStr(cumdx)+"\t"+DoubleToStr(cumdy)+"\t"+IntToStr(MCsteps[s])+"\t"+FloatToStr(des[s])+"\n";
        FileWrite(f,st.c_str(),st.length());
      }
    }
  printf("Rough estimate for conductivity_x: %8.5f\n", cumdx/ts[steps-1]);
  printf("Rough estimate for conductivity_y: %8.5f\n", cumdy/ts[steps-1]);

  FileClose(f);
}

void MKcurrentRejection::sample2SiteJumps(int steps)
{

    int s,i,j;
    double **jumps;
    int L=3;
    jumps = new double*[2*L+1];
    for (i = 0; i < 2*L+1; i++) {
        jumps[i] = new double[2*L+1];
    }

    for (s=0; s < steps; s++){
            for (i=-L; i<L+1; i++){
                for (j = -L; j < L+1; j++) {
                    if (dxs[s] == i && dys[s] == j){
                        jumps[j+L][i+L]++;
                    }
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


}


void MKcurrentRejection::runCurrentMaps(int steps,double &E, double &t)
{


  int s, MCs, i=0, j=0, dx=0, dy=0;
  double dE,dt;
  bool jump;
  tMC = 1/(L*L*(Nmem-1)*es->occnum*(1-es->occnum)*maxAccept);
  printf("t: %.7f\n", tMC);


  meanJumpLength = 0;
  meandx = 0;
  meandy = 0;
  testedNumberOf2Site = 0;
  numberOf2Site = 0;



  double meanDiscdx = 0, meanDiscdy = 0;
  double meanAbsdx = 0, meanAbsdy = 0;

  double meanMCs = 0;
  double actualmeanMCs = 0;



  for (s = 0; s < steps; s++)
  {
      MCs = 0;
      jump = false;
      while (!jump)
      {

          MCs++;

          testedNumberOf2Site++;
          getjump(i,j,dx,dy);
          jump=testjump(i,j,dx,dy);

          meanAbsdx +=dx; meanAbsdy +=dy;
          if (jump == true) {
              numberOf2Site++;
          }
          else meanDiscdx += dx; meanDiscdy += dy;

      }
      if (s > 1000000){
          updateMovement(movement, i,j, 1);
          updateMovement(xmovement,i,j, dx);
          updateMovement(ymovement,i,j,-dy);
      }

      dE = es->hopp(i,-1,0,j,0,0);

      dt = MCs*tMC;
      E += dE;
      t += dt;
      meanMCs += MCs;


      meanJumpLength += jl[(dx+maxJL+1) + JL2*(dy+maxJL+1)];
      meandx  += dx;
      meandy  += dy;

      actualmeanMCs += MCs;

      from[s] = i;
      to[s]   = j;

      ts[s] = t;
      energy[s] = E;
      des[s] = dE;
      dxs[s] = dx;
      dys[s] = dy;
      MCsteps[s]=MCs;

      if (s%10000==0) {
          printf("%6d E=%le t=%.4e MCs=%.3f\n",s,E,t,actualmeanMCs/10000);
          actualmeanMCs =0;
      }

     // printf("Step: %d: %d %d dE: %le I:%le MCs: %d\n", s, i, j,energy[s],dx,MCsteps[s]);
    }
  printf("\n-------Finishing timesteps--------\n\n");


  meanJumpLength /= steps;
  meandx /= steps;
  meandy /= steps;

  printf("\nProposed  dx: %8.5f, dy:  %8.5f\n", meanAbsdx /testedNumberOf2Site, meanAbsdy /testedNumberOf2Site);
  printf("Discarded dx: %8.5f, dy:  %8.5f\n", meanDiscdx/testedNumberOf2Site, meanDiscdy/testedNumberOf2Site);
  printf("Mean dx:      %8.5f, dy:  %8.5f\n", meandx, meandy);
  printf("\nMean jump length:                    %8.5f\n", meanJumpLength);
  printf("Acceptance ratio: %.5f\n", double(numberOf2Site)/testedNumberOf2Site);

}

void MKcurrentRejection::normalizeMap(double** &vector, int steps)
{
    int x, y, i;
    for (i = 0; i < N; i++){
        x = i%L;
        y = i/L;
        vector[y][x] /= steps;
    }
}

void MKcurrentRejection::updateMovement(double** &vector, int i, int j, int dWhat)
{
    int x0, y0, x1, y1;
    x0 = i%L; y0 = i/L; x1 = j%L; y1 = j/L;
    vector[y0][x0] += dWhat;
    vector[y1][x1] += dWhat;
}



void MKcurrentRejection::writeMapToFile(double **vector, int size, string filename)
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

void MKcurrentRejection::saveMaps(string outputprefix, string outputendfix, int steps)
{


    normalizeMap(movement,    steps-1000000);
    normalizeMap(xmovement,   steps-1000000);
    normalizeMap(ymovement,   steps-1000000);

    writeMapToFile(movement,  L,  outputprefix   +   "allJumps_" +  outputendfix  + ".dat");
    writeMapToFile(xmovement,  L,  outputprefix  +  "xallJumps_" +  outputendfix  + ".dat");
    writeMapToFile(ymovement,  L,  outputprefix  +  "yallJumps_" +  outputendfix  + ".dat");



}


void MKcurrentRejection::heatMapToFile(string filename, int steps)
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

