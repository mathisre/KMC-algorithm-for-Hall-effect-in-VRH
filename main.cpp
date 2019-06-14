//---------------------------------------------------------------------------
#include <stdio.h>
#include <math.h>
#include <string.h>

#pragma hdrstop


#include "systemclass.h"
#include "CGanalysis.h"
#include "fileutils.h"
#include "stringutils.h"
#include "treeutils.h"
#include "paramfile.h"
#include "MKcurrent.h"
#include "params.h"
#include "run.h"
#include "randomc.h"



int main(int argc, char* argv[])
{
  params *p=new params();

  p->readparams(argc,argv);

  switch (p->ctype) {
         case 1: run::runCurrent(p);                       break;
         case 2: run::runCurrentHeatMap(p);                break;
         case 3: run::runCurrentStateFile(p);              break;
         case 4: run::runCurrentStateSpesFile(p);          break;
         case 5: run::runCurrentRejection(p);              break;
         case 6: run::runCurrentRejectionMaps(p);          break;
         case 7: run::runCurrentRejectionRandomLattice(p); break;
         case 8: run::runCurrentRandPos(p);                break;
  }
  delete p;
  return 0;
}
