//---------------------------------------------------------------------------
#ifndef MKcurrentH
#define MKcurrentH
//---------------------------------------------------------------------------
#include "systemclass.h"
#include "treeutils.h"
#include "randomc.h"
#include "mersenne.h"
#include "vector"

using namespace std;

#define MAXFINDMIN 1000

//class for current within a systemclass


class MKcurrent {

public:
  MKcurrent(int steps);
  ~MKcurrent();

	void setES(ESystem *e);

    void init(int D1, int L1, int N1, int maxJL1, double a1, double beta1, double maxProbability1, bool doTriJumps);

    void setE(double Ex1, double Ey1, double Ez1);
    void setH(double Hx1, double Hy1, double Hz1);
	void setE0(double Ex1, double Ey1, double Ez1);
	void setOmega(double omega1);
	void setWritelines(int writelines1);
	void setRatefun(int ratefun1);

    void runCurrent(int steps, double &E, double &t);
    void runCurrentMeasure(int steps, double &E, double &t);
    void runCurrentTrace(int steps, double &E, double &t, string filename);
//    void runCurrentMaps(int steps, double &E, double &t);
	void runCurrentAC(int steps, double &E, double &t, int &dx);

	//	double getDeltat(){return deltat;}

	//	void toFile(string filename, int steps);
	void jumpsToFileSmall(string filename, int steps);
	void jumpsToFileAC(string filename, int steps);
	void heatMapToFile(string filename, int steps);
    void setMTseed(int seed);
    void initPositions(string filename);
    void updatePositions(int i, int j,int dx, int dy);
    void writePositions(string filename, int s);
    void closePositions(FILE* f);
    void sample3SiteJumpRates(int steps, string proposedFilename, string discardedFilename, string performedFilename);
    void sample2SiteJumpRates(int steps);
    void writeAreaHistogramToFile(string filename);
    void sampleProposed3SiteJumps(int steps);
    void comparePerformedWithInitialized(int steps);

    void sample3SiteJump(int dx, int dy, int dxI, int dyI, int jumpNumber, bool jump);
    void sample2SiteJump(int dx, int dy, bool jump);

    void printWeightedAreaTimeData();
    void finishTimeStepsAndSample(int steps);


    void updateMovement(double** &vector, int i, int j, int dWhat);
    void update2Site(double** &vector, int i, int j, int dWhat);
    void updateTriangles(double** &vector, int i, int k, int j, int dWhat);

    void normalizeMap(double** &vector, int steps);
    void writeMapToFile(double** vector, int size, string filename);
    void saveMaps(string outputprefix, string inputprefix, int steps);




private:
    CRandomMersenne *RanGen;
    //	double inline getJL(int dx, int dy, int dz);
    //      savedstate inline *checkstateexistance();
    void inline getjump(int &i, int &j, int &dx, int &dy);
    void inline getjump3sites(int &i, int &j, int &k, int &dx, int &dy, int &dxI, int &dyI, int &jumpNumber);
    bool inline testjump(int i, int j, int dx, int dy);
    bool inline testjump3sites(int i, int j, int k, int dx, int dy, int dxI, int dyI, int jumpNumber);


	//	double inline calcRate2d(int n1, int n2, int dx, int dy);

	int maxJL, JL2;
	double *jl; //jumpLength
    double *GammaT,*GammaT3Sites; //Tunneling Gamma (see Tsiganskov, Pazy, Laikhmann, Efros
    double GammaTtot, maxProbability, GammaTtot3Sites, totGammaTtot;
    int Nmem,Nmem3Sites, nGamma3Sites;
    int *dxInter,*dyInter,*dxFinal,*dyFinal;
    bool *occupied;
	double A, beta; //2/a, for the localization...
    double *jumpArea;
    double occupiedJumpCounter, notOccupiedJumpCounter;
    double dyOccupied, dyNotOccupied;
    double dxOccupied, dxNotOccupied;

    double deltat;

	int *from;
	int *to;
	double *energy,*ts,*des;
	int *MCsteps;
    int *dxs, *dxIs, *dys, *dyIs;
    int *dx3s, *dy3s;


    double Ex, Ey, Ez, E0x, E0y, E0z, Hx, Hy, Hz;
	double tMC,omega;
	int writelines,ratefun;
    double meanJumpLength, meanInterJumpLength, meanWeightedJL_x, meanWeightedJL_y;
    double meanArea, meanWeightedArea;
    double meandx, meandy,meandxI, meandyI;
    int L,N,D;
    int electronPositions[5000][4];
    double meanSomething;

    double *areaHistogram, bucketSize;
    int bins;

    int *testdxs, *testdys, *testdxIs, *testdyIs;
    int test3JumpCounter;
    int *discdxs, *discdys, *discdxIs, *discdyIs;
    int disc3JumpCounter;
    double *transitionTracer, **transitionJumps, **transitionCounter;
    double **rateMatrix;

    double driftAngle;

    double *meanWeightedAreaPer;
    double *dyPer;
    double dy3, dx3, dy2, dx2;

    double meanInterFinalJumpLength ;
    double meanDiscdx, meanDiscdy;
    double meanAbsdx, meanAbsdy;

    double meanDiscdxI , meanDiscdyI ;
    double meanAbsdxI, meanAbsdyI;

    double zeroAreaCounter;
    double meanJumpLength3,meanJumpLength2;


    int testedNumberOf2Site, testedNumberOf3Site, numberOf2Site, numberOf3Site;

    double **positions ;
    double **movement  ;
    double **xmovement ;
    double **ymovement ;
    double **twoSite   ;
    double **x2Site    ;
    double **y2Site    ;
    double **triangles ;
    double **xTriangles;
    double **yTriangles;


	ESystem *es;

};

#endif
