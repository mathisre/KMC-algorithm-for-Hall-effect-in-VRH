//---------------------------------------------------------------------------
#ifndef MKcurrentRandPosH
#define MKcurrentRandPosH
//---------------------------------------------------------------------------
#include "systemclass.h"
#include "treeutils.h"
#include "randomc.h"
#include "mersenne.h"
#include "vector"

using namespace std;

#define MAXFINDMIN 1000

//class for current within a systemclass


class MKcurrentRandPos {

public:
  MKcurrentRandPos(int steps);
  ~MKcurrentRandPos();

	void setES(ESystem *e);

    void init(int D1, int L1, int N1, int maxJL1, double a1, double beta1, double maxProbability1, bool doTriJumps, double shift);

    void setE(double Ex1, double Ey1, double Ez1);
    void setH(double Hx1, double Hy1, double Hz1);
	void setE0(double Ex1, double Ey1, double Ez1);
	void setOmega(double omega1);
	void setWritelines(int writelines1);
	void setRatefun(int ratefun1);

    void runCurrent(int steps, double &E, double &t);
    void runCurrentTrace(int steps, double &E, double &t, string filename);
	void runCurrentAC(int steps, double &E, double &t, int &dx);

	void jumpsToFileSmall(string filename, int steps);
	void jumpsToFileAC(string filename, int steps);
	void heatMapToFile(string filename, int steps);
	void setMTseed(int seed);
    void updateMap();
    void updateMovement(int i, int j);
    void writeGammaToFile(double Hz);
    void writeGamma2SiteToFile(double Hz);
    void initPositions(string filename);
    void updatePositions(int i, int j,int dx, int dy);
    void writePositions(string filename, int s);
    void closePositions(FILE* f);






private:
    CRandomMersenne *RanGen;
    void inline getjump(int &i, int &j, double &dx, double &dy);
    bool inline testjump(int i, int j, double dx, double dy);

	//	double inline calcRate2d(int n1, int n2, int dx, int dy);

	int maxJL, JL2;
    double *jl, maxProbability; //jumpLength
    double  totGammaTtot, abstotGammaTtot;
    int *Nmem;
    double **dxToSite, **dyToSite;
	double A, beta; //2/a, for the localization...
    double meanArea;

	double deltat;

	int *from;
	int *to;
	double *energy,*ts,*des;
	int *MCsteps;
    double *dxs;
    double *dys;

    double randShift;


    double Ex, Ey, Ez, E0x, E0y, E0z, Hx, Hy, Hz;
	double tMC,omega;
	int writelines,ratefun;
    double meanJumpLength, meanInterJumpLength;
    double meandx, meandy,meandxI, meandyI;
    int L,N,D;
    int electronPositions[5000][4];
	//	char *n;

	ESystem *es;
    double **sitePositions, **GammaT, *GammaTtot, **jumpArea, **distanceMatrix;
    int  **finalSite;
    double meanTotalGamma;
    double *GammaTtotForN;


    double meanSomething;

};

#endif
