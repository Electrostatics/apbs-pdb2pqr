/* ///////////////////////////////////////////////////////////////////////////
/// APBS -- Adaptive Poisson-Boltzmann Solver
///
///  Nathan A. Baker (nbaker@wasabi.ucsd.edu)
///  Dept. of Chemistry and Biochemistry
///  Dept. of Mathematics, Scientific Computing Group
///  University of California, San Diego 
///
///  Additional contributing authors listed in the code documentation.
///
/// Copyright © 1999. The Regents of the University of California (Regents).
/// All Rights Reserved. 
/// 
/// Permission to use, copy, modify, and distribute this software and its
/// documentation for educational, research, and not-for-profit purposes,
/// without fee and without a signed licensing agreement, is hereby granted,
/// provided that the above copyright notice, this paragraph and the
/// following two paragraphs appear in all copies, modifications, and
/// distributions.
/// 
/// IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
/// SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
/// ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
/// REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  
/// 
/// REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
/// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
/// PARTICULAR PURPOSE.  THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF
/// ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS".  REGENTS HAS NO OBLIGATION
/// TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
/// MODIFICATIONS. 
///
/// rcsid="$Id$"
//////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// File:     routines.h
//
// Purpose:  APBS ``front end" using formatted input files auxiliary routines
//           header file
//
// Author:   Nathan Baker
//
// rcsid="$Id$"
/////////////////////////////////////////////////////////////////////////// */

#ifndef _APBSROUTINES_H_
#define _APBSROUTINES_H_

#include "apbscfg.h"
#include "apbs/apbs.h"  
#include "apbs/nosh.h"  
#include "apbs/mgparm.h"  
#include "apbs/pbeparm.h"  
#include "apbs/femparm.h"  

#define APBSRC 13

/* ///////////////////////////////////////////////////////////////////////////
// Class:    AtomForce
//
// Purpose:  Container class for atomic forces
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
typedef struct AtomForce {
   double ibForce[3];
   double qfForce[3];
   double dbForce[3];
   double npForce[3];
} AtomForce;


/* ///////////////////////////////////////////////////////////////////////////
// PUBLIC ROUTINES
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VEXTERNC int loadMolecules(Vcom *com, NOsh *nosh, Valist *alist[NOSH_MAXMOL]);
VEXTERNC void printPBEPARM(Vcom *com, PBEparm *pbeparm);
VEXTERNC void printMGPARM(Vcom *com, MGparm *mgparm, double realCenter[3]);
VEXTERNC int initMG(Vcom *com, int i, NOsh *nosh, MGparm *mgparm,
  PBEparm *pbeparm, double realCenter[3], Vpbe *pbe[NOSH_MAXCALC],
  Valist *alist[NOSH_MAXMOL], Vpmgp *pmgp[NOSH_MAXCALC], 
  Vpmg *pmg[NOSH_MAXCALC]);
VEXTERNC int solveMG(Vcom *com, Vpmg *pmg);
VEXTERNC int setPartMG(Vcom *com, MGparm *mgparm, Vpmg *pmg);
VEXTERNC int energyMG(Vcom *com, NOsh* nosh, int icalc, Vpmg *pmg,
  int *nenergy, double *totEnergy, double *qfEnergy, double *qmEnergy,
  double *dielEnergy);
VEXTERNC int forceMG(Vcom *com, Vmem *mem, NOsh *nosh, PBEparm *pbeparm, 
  Vpmg *pmg, int *nforce, AtomForce **atomForce, Valist *alist[NOSH_MAXMOL]);
VEXTERNC int writepotMG(Vcom *com, NOsh *nosh, PBEparm *pbeparm, Vpmg *pmg);
VEXTERNC int writeaccMG(Vcom *com, NOsh *nosh, PBEparm *pbeparm, Vpmg *pmg);
VEXTERNC int printEnergy(Vcom *com, NOsh *nosh, double totEnergy[NOSH_MAXCALC],
  int i);

#endif
