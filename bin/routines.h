/**
 * @defgroup  Frontend  High-level front-end routines
 */

/**
 *  @file    routines.h
 *  @author  Nathan Baker
 *  @brief   Header file for front end auxiliary routines
 *  @ingroup  Frontend
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
 * Copyright (c) 2002-2006.  Washington University in St. Louis.
 * All Rights Reserved.
 * Portions Copyright (c) 1999-2002.  The Regents of the University of
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
 * Linking APBS statically or dynamically with other modules is making a
 * combined work based on APBS. Thus, the terms and conditions of the GNU
 * General Public License cover the whole combination.
 * 
 * SPECIAL GPL EXCEPTION
 * In addition, as a special exception, the copyright holders of APBS
 * give you permission to combine the APBS program with free software
 * programs and libraries that are released under the GNU LGPL or with
 * code included in releases of ISIM, Ion Simulator Interface, PMV, PyMOL
 * SMOL, VMD, and Vision. Such combined software may be linked with APBS and 
 * redistributed together in original or modified form as mere aggregation
 * without requirement that the entire work be under the scope of the GNU 
 * General Public License. This special exception permission is also extended
 * to any software listed in the SPECIAL GPL EXCEPTION clauses by the PMG,
 * FEtk, MC, or MALOC libraries.
 * 
 * Note that people who make modified versions of APBS are not obligated
 * to grant this special exception for their modified versions; it is
 * their choice whether to do so. The GNU General Public License gives
 * permission to release a modified version without this exception; this
 * exception also makes it possible to release a modified version which
 * carries forward this exception.
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

/**
 * @brief  Return code for APBS during failure
 * @ingroup  Frontend */
#define APBSRC 13

/** 
 * @brief  Set this macro to 1 for hierarchical basis, 0 for normal solver
 * @ingroup  Frontend */
#define USEHB 1

/**
 * @brief  Structure to hold atomic forces
 * @ingroup  Frontend
 * @author  Nathan Baker */
struct AtomForce {
   double ibForce[3];  /**< Ion-boundary force */
   double qfForce[3];  /**< Charge-field force */
   double dbForce[3];  /**< Dielectric boundary force */
   double npForce[3];  /**< Apolar force */
};

/**
 * @brief  Define AtomForce type
 * @ingroup  Frontend */
typedef struct AtomForce AtomForce;

/**
 * @brief  Load the molecules given in NOsh into atom lists
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  nosh  NOsh object with input file information
 * @param  alist  List of atom list objects
 * @returns  1 if successful, 0 otherwise */
VEXTERNC int loadMolecules(NOsh *nosh, Valist *alist[NOSH_MAXMOL]);

/**
 * @brief  Destroy the loaded molecules
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  nosh  NOsh object with input file information
 * @param  alist  List of atom list objects */
VEXTERNC void killMolecules(NOsh *nosh, Valist *alist[NOSH_MAXMOL]);

/**
 * @brief  Load the dielectric maps given in NOsh into grid objects
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  nosh  NOsh object with input file information
 * @param  dielXMap  List of x-shifted dielectric maps
 * @param  dielYMap  List of y-shifted dielectric maps
 * @param  dielZMap  List of x-shifted dielectric maps
 * @returns  1 if successful, 0 otherwise */
VEXTERNC int loadDielMaps(NOsh *nosh, Vgrid *dielXMap[NOSH_MAXMOL],
  Vgrid *dielYMap[NOSH_MAXMOL], Vgrid *dielZMap[NOSH_MAXMOL]);

/**
 * @brief  Destroy the loaded dielectric
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  nosh  NOsh object with input file information
 * @param  dielXMap  List of x-shifted dielectric maps
 * @param  dielYMap  List of y-shifted dielectric maps
 * @param  dielZMap  List of x-shifted dielectric maps */
VEXTERNC void killDielMaps(NOsh *nosh, Vgrid *dielXMap[NOSH_MAXMOL],
  Vgrid *dielYMap[NOSH_MAXMOL], Vgrid *dielZMap[NOSH_MAXMOL]);

/**
 * @brief  Load the kappa maps given in NOsh into grid objects
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  nosh  NOsh object with input file information
 * @param  kappa  List of kappa maps
 * @returns  1 if successful, 0 otherwise */
VEXTERNC int loadKappaMaps(NOsh *nosh, Vgrid *kappa[NOSH_MAXMOL]);

/**
 * @brief  Destroy the loaded kappa maps 
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  nosh  NOsh object with input file information
 * @param  kappa  List of kappa maps */
VEXTERNC void killKappaMaps(NOsh *nosh, Vgrid *kappa[NOSH_MAXMOL]);

/**
 * @brief  Load the charge maps given in NOsh into grid objects
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  nosh  NOsh object with input file information
 * @param  charge  List of kappa maps
 * @returns  1 if successful, 0 otherwise */
VEXTERNC int loadChargeMaps(NOsh *nosh, Vgrid *charge[NOSH_MAXMOL]);

/**
 * @brief  Destroy the loaded charge maps 
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  nosh  NOsh object with input file information
 * @param  charge  List of charge maps */
VEXTERNC void killChargeMaps(NOsh *nosh, Vgrid *charge[NOSH_MAXMOL]);

/**
 * @brief  Print out generic PBE params loaded from input
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  pbeparm  PBEparm object */
VEXTERNC void printPBEPARM(PBEparm *pbeparm);

/**
 * @brief  Print out MG-specific params loaded from input
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  realCenter  Center of mesh for actual calculation
 * @param  mgparm  MGparm object */
VEXTERNC void printMGPARM(MGparm *mgparm, double realCenter[3]);

/**
 * @brief  Initialize an MG calculation
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param icalc  Index of calculation in pmg/pmpg arrays
 * @param nosh  Object with parsed input file parameters
 * @param mgparm  Object with MG-specific parameters
 * @param pbeparm  Object with generic PBE parameters 
 * @param realCenter  The actual center of the current mesh
 * @param pbe  Array of Vpbe objects (one for each calc)
 * @param alist  Array of atom lists
 * @param dielXMap  Array of x-shifted dielectric maps 
 * @param dielYMap  Array of y-shifted dielectric maps 
 * @param dielZMap  Array of z-shifted dielectric maps 
 * @param kappaMap  Array of kappa maps 
 * @param chargeMap  Array of charge maps 
 * @param pmgp  Array of MG parameter objects (one for each calc)
 * @param pmg  Array of MG objects (one for each calc)
 * @return  1 if succesful, 0 otherwise */
VEXTERNC int initMG(int icalc, NOsh *nosh, MGparm *mgparm,
  PBEparm *pbeparm, double realCenter[3], Vpbe *pbe[NOSH_MAXCALC],
  Valist *alist[NOSH_MAXMOL], Vgrid *dielXMap[NOSH_MAXMOL], 
  Vgrid *dielYMap[NOSH_MAXMOL], Vgrid *dielZMap[NOSH_MAXMOL], 
  Vgrid *kappaMap[NOSH_MAXMOL], Vgrid *chargeMap[NOSH_MAXMOL], 
  Vpmgp *pmgp[NOSH_MAXCALC], Vpmg *pmg[NOSH_MAXCALC]);

/**
 * @brief  Kill structures initialized during an MG calculation
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param nosh  Object with parsed input file parameters
 * @param pbe  Array of Vpbe objects (one for each calc)
 * @param pmgp  Array of MG parameter objects (one for each calc)
 * @param pmg  Array of MG objects (one for each calc) */
VEXTERNC void killMG(NOsh *nosh, Vpbe *pbe[NOSH_MAXCALC],
  Vpmgp *pmgp[NOSH_MAXCALC], Vpmg *pmg[NOSH_MAXCALC]);

/**
 * @brief  Solve the PBE with MG
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param nosh  Object with parsed input file parameters
 * @param pmg  MG objects for this calculation
 * @param type  Type of MG calculation
 * @return  1 if successful, 0 otherwise */
VEXTERNC int solveMG(NOsh *nosh, Vpmg *pmg, MGparm_CalcType type);

/**
 * @brief  Set MG partitions for calculating observables and performing I/O
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param nosh  Object with parsed input file parameters
 * @param mgparm  MG parameters from input file
 * @param pmg  MG object 
 * @return  1 if successful, 0 otherwise */
VEXTERNC int setPartMG(NOsh *nosh, MGparm *mgparm, Vpmg *pmg);

/**
 * @brief  Calculate electrostatic energies from MG solution
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param nosh  Object with parsed input file parameters
 * @param icalc  Index of calculation 
 * @param pmg  MG object 
 * @param nenergy  Set to number of entries in energy arrays
 * @param totEnergy  Set to total energy (in kT)
 * @param qfEnergy  Set to charge-potential energy (in kT)
 * @param qmEnergy  Set to mobile ion energy (in kT)
 * @param dielEnergy  Set to polarization energy (in kT)
 * @return  1 if successful, 0 otherwise */
VEXTERNC int energyMG(NOsh* nosh, int icalc, Vpmg *pmg,
  int *nenergy, double *totEnergy, double *qfEnergy, double *qmEnergy,
  double *dielEnergy);

/**
 * @brief  Calculate apolar energies from MG solution
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param nosh  Object with parsed input file parameters
 * @param icalc  Index of calculation 
 * @param pmg  MG object 
 * @param nenergy  Set to number of entries in energy arrays
 * @param npEnergy  Set to apolar energy (in kT)
 * @return  1 if successful, 0 otherwise */
VEXTERNC int npenergyMG(NOsh* nosh, int icalc, Vpmg *pmg, int *nenergy, 
  double *npEnergy);

/** 
 * @brief  Kill arrays allocated for energies
 * @ingroup  Frontend
 * @author  Nathan Baker */
VEXTERNC void killEnergy();

/**
 * @brief  Calculate forces from MG solution
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  mem  Memory management object
 * @param  nosh  Parameters from input file
 * @param  pbeparm  Generic PBE parameters
 * @param  mgparm  MG-specific parmaeters
 * @param pmg  MG object
 * @param nforce  Set to number of forces in arrays
 * @param  atomForce  List of atom forces
 * @param  alist  List of atom lists
 * @return  1 if successful, 0 otherwise */
VEXTERNC int forceMG(Vmem *mem, NOsh *nosh, PBEparm *pbeparm,  MGparm *mgparm,
  Vpmg *pmg, int *nforce, AtomForce **atomForce, Valist *alist[NOSH_MAXMOL]);

/**
 * @brief  Free memory from MG force calculation
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  mem  Memory management object
 * @param  nosh  Parameters from input file
 * @param  nforce  Number of forces in arrays
 * @param  atomForce  List of atom forces */
VEXTERNC void killForce(Vmem *mem, NOsh *nosh, int nforce[NOSH_MAXCALC],
  AtomForce *atomForce[NOSH_MAXCALC]);

/**
 * @brief Store energy in arrays for future use
 * @ingroup Frontend
 * @author Todd Dolinsky */
VEXTERNC void storeAtomEnergy(
		Vpmg *pmg, /**< MG object */
		int icalc, /**< Calculation number */
		double **atomEnergy, /**< Pointer to storage array of doubles */
		int *nenergy /**< Stores number of atoms per calc */
		);

/**
 * @brief Write out information to a flat file
 * @ingroup Frontend
 * @author  Todd Dolinsky
 * @param nosh  Parameters from input file
 * @param com   The communications object
 * @param fname The target XML file name
 * @param totEnergy An array with per-calc total energies (in kT)
 * @param qfEnergy  An array with per-calc charge-potential energies (in kT)
 * @param qmEnergy  An array with per-calc mobile energies (in kT)
 * @param dielEnergy  An array with per-calc polarization energies (in kT)
 * @param nenergy  An array containing the number of atoms per-calc
 * @param atomEnergy An array containing per-atom energies (in KT) per calc
 * @param nforce  An array containing the number of forces calculated per-calc
 * @param atomForce An array containing per-atom forces per calc
 * @return 1 if successful, 0 otherwise */
VEXTERNC int writedataFlat(NOsh *nosh, Vcom *com, const char *fname, 
  double totEnergy[NOSH_MAXCALC], double qfEnergy[NOSH_MAXCALC], 
  double qmEnergy[NOSH_MAXCALC], double dielEnergy[NOSH_MAXCALC],
  int nenergy[NOSH_MAXCALC], double *atomEnergy[NOSH_MAXCALC],
  int nforce[NOSH_MAXCALC], AtomForce *atomForce[NOSH_MAXCALC]);

/**
 * @brief Write out information to an XML file
 * @ingroup Frontend
 * @author  Todd Dolinsky
 * @param nosh  Parameters from input file
 * @param com   The communications object 
 * @param fname The target XML file name
 * @param totEnergy An array with per-calc total energies (in kT)
 * @param qfEnergy  An array with per-calc charge-potential energies (in kT)
 * @param qmEnergy  An array with per-calc mobile energies (in kT)
 * @param dielEnergy  An array with per-calc polarization energies (in kT)
 * @param nenergy  An array containing the number of atoms per-calc
 * @param atomEnergy An array containing per-atom energies (in KT) per calc
 * @param nforce  An array containing the number of forces calculated per-calc
 * @param atomForce An array containing per-atom forces per calc
 * @return 1 if successful, 0 otherwise */
VEXTERNC int writedataXML(NOsh *nosh, Vcom *com, const char *fname, 
  double totEnergy[NOSH_MAXCALC], double qfEnergy[NOSH_MAXCALC], 
  double qmEnergy[NOSH_MAXCALC], double dielEnergy[NOSH_MAXCALC],
  int nenergy[NOSH_MAXCALC], double *atomEnergy[NOSH_MAXCALC],
  int nforce[NOSH_MAXCALC], AtomForce *atomForce[NOSH_MAXCALC]);

/**
 * @brief  Write out observables from MG calculation to file
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  rank  Processor rank (if parallel calculation)
 * @param  nosh  Parameters from input file
 * @param  pbeparm  Generic PBE parameters
 * @param pmg  MG object
 * @return  1 if successful, 0 otherwise */
VEXTERNC int writedataMG(int rank, NOsh *nosh, PBEparm *pbeparm, Vpmg *pmg);

/**
 * @brief  Write out operator matrix from MG calculation to file
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  rank  Processor rank (if parallel calculation)
 * @param  nosh  Parameters from input file
 * @param  pbeparm  Generic PBE parameters
 * @param pmg  MG object
 * @return  1 if successful, 0 otherwise */
VEXTERNC int writematMG(int rank, NOsh *nosh, PBEparm *pbeparm, Vpmg *pmg);

/** 
 * @brief  Combine and pretty-print energy data
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  com  Communications object
 * @param  nosh  Parameters from input file
 * @param  totEnergy  Array of energies from different calculations
 * @param  iprint  Index of energy statement to print
 * @return  1 if successful, 0 otherwise */
VEXTERNC int printEnergy(
						 Vcom *com,
						 NOsh *nosh,
						 double totEnergy[NOSH_MAXCALC], 
						 int iprint
						 );

/** 
 * @brief  Combine and pretty-print force data
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  com  Communications object
 * @param  nosh  Parameters from input file
 * @param  nforce  Number of forces calculated
 * @param  atomForce  Array of force structures
 * @param  i  Index of force statement to print
 * @return  1 if successful, 0 otherwise */
VEXTERNC int printForce(Vcom *com, NOsh *nosh, int nforce[NOSH_MAXCALC],
  AtomForce *atomForce[NOSH_MAXCALC], int i);

/**
 * @brief  Wrapper to start MALOC Vio layer
 * @ingroup  Frontend
 * @author  Nathan Baker and Robert Konecny */
VEXTERNC void startVio();

#ifdef HAVE_MC_H
/**
 * @brief  Print out FE-specific params loaded from input
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  icalc  Calculation index
 * @param  nosh  Master parameter object
 * @param feparm  FE-specific parameters 
 * @param fetk  Array of FE solver objects  */
VEXTERNC void printFEPARM(int icalc, NOsh *nosh, FEMparm *feparm,
  Vfetk *fetk[NOSH_MAXCALC]);

/**
 * @brief  Calculate electrostatic energies from FE solution
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param nosh  Object with parsed input file parameters
 * @param icalc  Index of calculation 
 * @param fetk  FE object  array
 * @param nenergy  Set to number of entries in energy arrays
 * @param totEnergy  Set to total energy (in kT)
 * @param qfEnergy  Set to charge-potential energy (in kT)
 * @param qmEnergy  Set to mobile ion energy (in kT)
 * @param dielEnergy  Set to polarization energy (in kT)
 * @bug  "calcenergy 2" does not work
 * @return  1 if successful, 0 otherwise */
VEXTERNC int energyFE(NOsh* nosh, int icalc, Vfetk *fetk[NOSH_MAXCALC],
  int *nenergy, double *totEnergy, double *qfEnergy, double *qmEnergy,
  double *dielEnergy);

/**
 * @brief  Initialize FE solver objects
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  icalc  Index in pbe, fetk to initialize -- calculation index
 * @param nosh  Master parameter object
 * @param feparm  FE-specific parameters 
 * @param pbeparm  Generic PBE parameters
 * @param pbe  Array of PBE objects 
 * @param alist Array of atom lists 
 * @param fetk  Array of FE solver objects 
 * @bug  THIS FUNCTION IS HARD-CODED TO SOLVE LRPBE
 * @todo  THIS FUNCTION IS HARD-CODED TO SOLVE LRPBE
 * @return  1 if successful, 0 otherwise */
VEXTERNC int initFE(int icalc, NOsh *nosh, FEMparm *feparm, PBEparm *pbeparm,
  Vpbe *pbe[NOSH_MAXCALC], Valist *alist[NOSH_MAXMOL], 
  Vfetk *fetk[NOSH_MAXCALC]);

/**
 * @brief  Pre-refine mesh before solve
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  i  Calculation index
 * @param  nosh  Master parameter object
 * @param feparm  FE-specific parameters 
 * @param fetk  Array of FE solver objects 
 * @return  1 if successful, 0 otherwise */
VEXTERNC int preRefineFE(int i, NOsh *nosh, FEMparm *feparm,
  Vfetk *fetk[NOSH_MAXCALC]);

/**
 * @brief  Partition mesh (if applicable)
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  i  Calculation index
 * @param  nosh  Master parameter object
 * @param feparm  FE-specific parameters 
 * @param fetk  Array of FE solver objects 
 * @return  1 if successful, 0 otherwise */
VEXTERNC int partFE(int i, NOsh *nosh, FEMparm *feparm,
  Vfetk *fetk[NOSH_MAXCALC]);

/**
 * @brief  Solve-estimate-refine
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  i  Calculation index
 * @param  nosh  Master parameter object
 * @param feparm  FE-specific parameters 
 * @param pbeparm  Generic PBE parameters
 * @param fetk  Array of FE solver objects 
 * @return  1 if successful, 0 otherwise */
VEXTERNC int solveFE(int i, NOsh *nosh, PBEparm *pbeparm, FEMparm *feparm,
  Vfetk *fetk[NOSH_MAXCALC]);

/**
 * @brief  Estimate error, mark mesh, and refine mesh after solve
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  icalc  Calculation index
 * @param  nosh  Master parameter object
 * @param feparm  FE-specific parameters 
 * @param fetk  Array of FE solver objects 
 * @return  1 if successful, 0 otherwise -- note that a 0 will likely imply
 * that either the max number of vertices have been met or no vertices were
 * marked for refinement.  In either case, this should not be treated as a
 * fatal error.  */
VEXTERNC int postRefineFE(int icalc, NOsh *nosh, FEMparm *feparm,
  Vfetk *fetk[NOSH_MAXCALC]);

/**
 * @brief  Write FEM data to files
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  rank  Rank of processor (for parallel runs)
 * @param  nosh  NOsh object
 * @param  pbeparm  PBEparm object
 * @param  fetk  FEtk object (with solution)
 * @return  1 if successful, 0 otherwise */
VEXTERNC int writedataFE(int rank, NOsh *nosh, PBEparm *pbeparm, Vfetk *fetk);
#endif

#endif
