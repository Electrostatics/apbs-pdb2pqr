/**
 *  @file    routines.h
 *  @author  Nathan Baker
 *  @brief   Header file for front end auxiliary routines
 *  @version $Id$
 *  @attention
 *  @verbatim
 *
 * APBS -- Adaptive Poisson-Boltzmann Solver
 *
 * Nathan A. Baker (baker@biochem.wustl.edu)
 * Dept. of Biochemistry and Molecular Biophysics
 * Center for Computational Biology
 * Washington University in St. Louis
 *
 * Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 2003.  Washington University in St. Louis.
 * All Rights Reserved.
 * Portions Copyright (c) 1999-2003.  The Regents of the University of
 * California.  
 * Portions Copyright (c) 1995.  Michael Holst.
 *
 * This file is part of APBS.
 *
 * APBS is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * APBS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with APBS; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 *
 * @endverbatim
 */

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
VEXTERNC int loadMolecules(NOsh *nosh, Valist *alist[NOSH_MAXMOL]);
VEXTERNC void killMolecules(NOsh *nosh, Valist *alist[NOSH_MAXMOL]);
VEXTERNC int loadDielMaps(NOsh *nosh, Vgrid *dielXMap[NOSH_MAXMOL],
Vgrid *dielYMap[NOSH_MAXMOL], Vgrid *dielZMap[NOSH_MAXMOL]);
VEXTERNC void killDielMaps(NOsh *nosh, Vgrid *dielXMap[NOSH_MAXMOL],
Vgrid *dielYMap[NOSH_MAXMOL], Vgrid *dielZMap[NOSH_MAXMOL]);
VEXTERNC int loadKappaMaps(NOsh *nosh, Vgrid *kappa[NOSH_MAXMOL]);
VEXTERNC void killKappaMaps(NOsh *nosh, Vgrid *kappa[NOSH_MAXMOL]);
VEXTERNC int loadChargeMaps(NOsh *nosh, Vgrid *charge[NOSH_MAXMOL]);
VEXTERNC void killChargeMaps(NOsh *nosh, Vgrid *charge[NOSH_MAXMOL]);
VEXTERNC void printPBEPARM(PBEparm *pbeparm);
VEXTERNC void printMGPARM(MGparm *mgparm, double realCenter[3]);
VEXTERNC int initMG(int i, NOsh *nosh, MGparm *mgparm,
  PBEparm *pbeparm, double realCenter[3], Vpbe *pbe[NOSH_MAXCALC],
  Valist *alist[NOSH_MAXMOL], Vgrid *dielXMap[NOSH_MAXMOL], 
  Vgrid *dielYMap[NOSH_MAXMOL], Vgrid *dielZMap[NOSH_MAXMOL], 
  Vgrid *kappaMap[NOSH_MAXMOL], Vgrid *chargeMap[NOSH_MAXMOL], 
  Vpmgp *pmgp[NOSH_MAXCALC], Vpmg *pmg[NOSH_MAXCALC]);
VEXTERNC void killMG(NOsh *nosh, Vpbe *pbe[NOSH_MAXCALC],
  Vpmgp *pmgp[NOSH_MAXCALC], Vpmg *pmg[NOSH_MAXCALC]);
VEXTERNC int solveMG(NOsh *nosh, Vpmg *pmg, int type);
VEXTERNC int setPartMG(NOsh *nosh, MGparm *mgparm, Vpmg *pmg);
VEXTERNC int energyMG(NOsh* nosh, int icalc, Vpmg *pmg,
  int *nenergy, double *totEnergy, double *qfEnergy, double *qmEnergy,
  double *dielEnergy);
VEXTERNC int npenergyMG(NOsh* nosh, int icalc, Vpmg *pmg, int *nenergy, double *npEnergy);
VEXTERNC void killEnergy();
VEXTERNC int forceMG(Vmem *mem, NOsh *nosh, PBEparm *pbeparm, 
  Vpmg *pmg, int *nforce, AtomForce **atomForce, Valist *alist[NOSH_MAXMOL]);
VEXTERNC void killForce(Vmem *mem, NOsh *nosh, int nforce[NOSH_MAXCALC],
  AtomForce *atomForce[NOSH_MAXCALC]);
VEXTERNC int writedataMG(int rank, NOsh *nosh, PBEparm *pbeparm, Vpmg *pmg);
VEXTERNC int writematMG(int rank, NOsh *nosh, PBEparm *pbeparm, Vpmg *pmg);
VEXTERNC int printEnergy(Vcom *com, NOsh *nosh, double totEnergy[NOSH_MAXCALC],
  int i);
VEXTERNC int printForce(Vcom *com, NOsh *nosh, int nforce[NOSH_MAXCALC],
  AtomForce *atomForce[NOSH_MAXCALC], int i);
VEXTERNC void startVio();

#endif
