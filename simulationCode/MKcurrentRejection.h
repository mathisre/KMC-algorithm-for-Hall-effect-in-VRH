//---------------------------------------------------------------------------
#ifndef MKcurrentRejectionH
#define MKcurrentRejectionH
//---------------------------------------------------------------------------
#include "systemclass.h"
#include "treeutils.h"
#include "randomc.h"
#include "mersenne.h"
#include "vector"

using namespace std;

#define MAXFINDMIN 1000

//class for current within a systemclass


class MKcurrentRejection {

public:
  MKcurrentRejection(int steps);
  ~MKcurrentRejection();

    void setES(ESystem *e);

    void init(int D1, int L1, int N1, int maxJL1, double a1, double beta1, double maxProbability1);
    void initRandomLattice(int D1, int L1, int N1, int maxJL1, double a1, double beta1, double maxProbability1, double Bz, bool doTriJumps);

    void setE(double Ex1, double Ey1, double Ez1);
    void setH(double Hx1, double Hy1, double Hz1);
    void setE0(double Ex1, double Ey1, double Ez1);
    void setOmega(double omega1);
    void setWritelines(int writelines1);
    void setRatefun(int ratefun1);

    void runCurrent(int steps, double &E, double &t);
    void runCurrentMeasure(int steps, double &E, double &t);
    void runCurrentRandomLattice(int steps, double &E, double &t);
    void runCurrentAC(int steps, double &E, double &t, int &dx);
    void runCurrentMaps(int steps, double &E, double &t);

    //	double getDeltat(){return deltat;}

    //	void toFile(string filename, int steps);
    void jumpsToFileSmall(string filename, int steps);
    void jumpsToFileSmallTriLattice(string filename, int steps);
    void jumpsToFileSmallRandomLattice(string filename, int steps);


    void jumpsToFileAC(string filename, int steps);
    void heatMapToFile(string filename, int steps);
    void setMTseed(int seed);
    void writeGammaToFile(double Hz);
    void writeGamma2SiteToFile(double Hz);
    void initPositions(string filename);
    void updatePositions(int i, int j,int dx, int dy);
    void writePositions(string filename, int s);
    void closePositions(FILE* f);
    void sample3SiteJumps(int steps);
    void sample2SiteJumps(int steps);
    void writeAreaHistogramToFile(string filename);



    void initMaps();
    void updateMovement(double** &vector, int i, int j, int dWhat);

    void normalizeMap(double** &vector, int steps);
    void writeMapToFile(double** vector, int size, string filename);
    void saveMaps(string outputprefix, string inputprefix, int steps);






private:
    CRandomMersenne *RanGen;
    //	double inline getJL(int dx, int dy, int dz);
    //      savedstate inline *checkstateexistance();
    void inline getjump(int &i, int &j, int &dx, int &dy);
    bool inline testjump(int i, int j, int dx, int dy);
    double threeSiteCorrectionRates(int i, int j, int dx, int dy, double dEij, double r_ij);

    void inline getjumpTriLattice(int &i, int &j, int &dx, int &dy);
    bool inline testjumpTriLattice(int i, int j, int dx, int dy);

    void inline getjumpRandomLattice(int &i, int &j, int &jumpNumber);
    bool inline testjumpRandomLattice(int i, int j, int jumpNumber);

    //	double inline calcRate2d(int n1, int n2, int dx, int dy);

    int maxJL, JL2;
    double *jl; //jumpLength
    double *GammaT,*GammaT3Sites; //Tunneling Gamma (see Tsiganskov, Pazy, Laikhmann, Efros
    double GammaTtot, maxProbability, GammaTtot3Sites, totGammaTtot;
    int Nmem,Nmem3Sites, nGamma3Sites;
    int testedNumberOf2Site, testedNumberOf3Site, numberOf2Site, numberOf3Site;
    double A, beta; //2/a, for the localization...
    double *jumpArea;

    double deltat;



    int *from;
    int *to;
    double *energy,*ts,*des;
    int *MCsteps;
    int *dxs,*dys;
    double *dxsTri,*dysTri;

    double *dxArray, *dyArray;
    int *dxIndex, *dyIndex;


    double Ex, Ey, Ez, E0x, E0y, E0z, Hx, Hy, Hz;
    double tMC,omega;
    int writelines,ratefun;
    double meanJumpLength, meanInterJumpLength;
    double meanArea;
    double meandx, meandy,meandxI, meandyI;
    int L,N,D;
    int electronPositions[5000][4];
    //	char *n;
    double meanSomething;

    double *areaHistogram, bucketSize;
    int bins;
    double maxAccept;
    double randShift;


    ESystem *es;

    double **sitePositions,  **distanceMatrix;
    double **dxInter,**dyInter,**dxFinal,**dyFinal;
    double **dxR, **dyR;
    int **dxRIndex, **dyRIndex;
    int **inter3Site, **final3Site,  **finalSite;
    double *dxsR, *dysR;
    int *NmemArray;
    double averageNmem;


    double **movement  ;
    double **xmovement ;
    double **ymovement ;


};

#endif
