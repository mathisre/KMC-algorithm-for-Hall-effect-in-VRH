//--------------
#ifndef RUN_H
#define RUN_H
//--------------

#include <stdio.h>
#include <math.h>
#include <string.h>

#pragma hdrstop


#include "systemclass.h"
#include "CGanalysis.h"
#include "fileutils.h"
#include "stringutils.h"
#include "paramfile.h"
#include "MKcurrent.h"
#include "paramfile.h"


class run
{
public:
    run();
    void FMcb(int step,int Xn,double E,int type,unsigned int secs,void *data);
    void printHistogram(int *histogram, int nrBoxes);
    static void printES(ESystem *es);

    static void runCurrent(class params *p);
    static void runCurrentMeasure(class params *p);
    static void runCurrentRandPos(class params *p);
    static void runCurrentTrace(class params *p);
    static void runCurrentHeatMap(params *p);
    static void runCurrentStateFile(params *p);
    static void runCurrentStateSpesFile(params *p);
    static void runCurrentRejection(params *p);
    static void runCurrentRejectionTriLattice(params *p);
    static void runCurrentRejectionRandomLattice(params *p);
    static void runCurrentRejectionMaps(params *p);

//    static void runCurrentMaps(params *p);
};

#endif // RUN_H
