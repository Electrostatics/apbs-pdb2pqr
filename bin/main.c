/**
 *  @file    main.c
 *  @author  Nathan Baker
 *  @brief   APBS "front end" program using formatted input files.
 * 
 *           This driver program represents a mish-mash of
 *           instructions for
 *           calculating electrostatic potentials, as well as free energies of
 *           binding and solvation.  It is invoked as:
 *
 *               apbs apbs.in
 *
 *           where apbs.in is a formatted input file (see documentation and
 *           examples).
 * 
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
 * Copyright (c) 2002-2007.  Washington University in St. Louis.
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

#include "apbscfg.h"
#include "apbs/apbs.h"  
#include "apbs/nosh.h"  
#include "apbs/mgparm.h"  
#include "apbs/pbeparm.h"  
#include "apbs/femparm.h"  

#include "routines.h"

VEMBED(rcsid="$Id$")

int main(
		 int argc, 
		 char **argv
		 ) 
{
	
	NOsh *nosh = VNULL;
	
	MGparm *mgparm = VNULL;
	FEMparm *feparm = VNULL;
	PBEparm *pbeparm = VNULL;
	APOLparm *apolparm = VNULL;
	Vparam *param = VNULL;
	
	Vmem *mem = VNULL;
	Vcom *com = VNULL;
	Vio *sock = VNULL;
#ifdef HAVE_MC_H
	Vfetk *fetk[NOSH_MAXCALC];
#else
	void *fetk[NOSH_MAXCALC];
#endif
	Vpmg *pmg[NOSH_MAXCALC];
	Vpmgp *pmgp[NOSH_MAXCALC];
	Vpbe *pbe[NOSH_MAXCALC];
	Valist *alist[NOSH_MAXMOL];
	Vgrid *dielXMap[NOSH_MAXMOL],*dielYMap[NOSH_MAXMOL],*dielZMap[NOSH_MAXMOL];
	Vgrid *kappaMap[NOSH_MAXMOL];
	Vgrid *chargeMap[NOSH_MAXMOL];
	char *input_path = VNULL;
	char *output_path = VNULL;
	int i, rank, size, isolve, k;
	size_t bytesTotal, highWater;
	Voutput_Format outputformat;
	
	int rc = 0;
	
	/* These variables require some explaining... The energy double arrays
		* store energies from the various calculations.  The energy int array
		* stores either a flag (0,1) displaying whether energies were calculated
		* or if PCE_COMPS is used, the number of atom energies stored
		* for the given calculation.  Likewise, the
		* force double arrays store forces from the various calcualtions.  The
		* force int array stores an integer which either says no calculation was
		* performed (0) or gives the number of entries in the force array for each
		* calculation */
	double qfEnergy[NOSH_MAXCALC], qmEnergy[NOSH_MAXCALC];
	double dielEnergy[NOSH_MAXCALC], totEnergy[NOSH_MAXCALC];
	AtomForce *atomForce[NOSH_MAXCALC];
	double *atomEnergy[NOSH_MAXCALC];
	int nenergy[NOSH_MAXCALC], nforce[NOSH_MAXCALC];
	/* THe real partition centers */
	double realCenter[3];
	
	/* Instructions: */
	char header[] = {"\n\n\
----------------------------------------------------------------------\n\
	APBS -- Adaptive Poisson-Boltzmann Solver\n\
	Version 0.5.1\n\
	\n\
	Nathan A. Baker (baker@biochem.wustl.edu)\n\
	Dept. Biochemistry and Molecular Biophysics\n\
	Center for Computational Biology\n\
	Washington University in St. Louis\n\
	\n\
	Additional contributing authors listed in the code documentation.\n\
	\n\
	Copyright (c) 2002-2007.  Washington University in St. Louis.\n\
	All Rights Reserved.\n\
	Portions Copyright (c) 1999-2002.  The Regents of the University of \n\
	California.\n\
	Portions Copyright (c) 1995.  Michael Holst.\n\
	\n\
	This program is free software; you can redistribute it and/or modify\n\
	it under the terms of the GNU General Public License as published by\n\
	the Free Software Foundation; either version 2 of the License, or\n\
	(at your option) any later version.\n\
	\n\
	This program is distributed in the hope that it will be useful,\n\
	but WITHOUT ANY WARRANTY; without even the implied warranty of\n\
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n\
	GNU General Public License for more details.\n\
	\n\
	You should have received a copy of the GNU General Public License\n\
	along with this program; if not, write to the Free Software\n\
	Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA\n\
	\n\
	Linking APBS statically or dynamically with other modules is making a\n\
	combined work based on APBS. Thus, the terms and conditions of the GNU\n\
	General Public License cover the whole combination.\n\
	\n\
	SPECIAL GPL EXCEPTION\n\
	In addition, as a special exception, the copyright holders of APBS\n\
	give you permission to combine the APBS program with free software\n\
	programs and libraries that are released under the GNU LGPL or with\n\
	code included in releases of ISIM, Ion Simulator Interface, PMV,\n\
	SMOL, VMD, and Vision. Such combined software may be linked with APBS and\n\
	redistributed together in original or modified form as mere aggregation\n\
	without requirement that the entire work be under the scope of the GNU\n\
	General Public License. This special exception permission is also extended\n\
	to any software listed in the SPECIAL GPL EXCEPTION clauses by the PMG,\n\
	FEtk, MC, or MALOC libraries.\n\
	\n\
	Note that people who make modified versions of APBS are not obligated\n\
	to grant this special exception for their modified versions; it is\n\
	their choice whether to do so. The GNU General Public License gives\n\
	permission to release a modified version without this exception; this\n\
	exception also makes it possible to release a modified version which\n\
	carries forward this exception.\n\
----------------------------------------------------------------------\n\
	APBS uses FETK (the Finite Element ToolKit) to solve the\n\
	Poisson-Boltzmann equation numerically.  FETK is a portable collection\n\
	of finite element modeling class libraries written in an object-oriented\n\
	version of C.  It is designed to solve general coupled systems of nonlinear\n\
	partial differential equations using adaptive finite element methods,\n\
	inexact Newton methods, and algebraic multilevel methods.  More information\n\
	about FEtk may be found at <http://www.FEtk.ORG>.\n\
----------------------------------------------------------------------\n\
Please cite your use of APBS as:\n\n\
	Baker NA, Sept D, Joseph S, Holst MJ, McCammon JA. Electrostatics of\n\
	nanosystems: application to microtubules and the ribosome. Proc.\n\
	Natl. Acad. Sci. USA 98, 10037-10041 2001.\n\
	\n\n"};
	char *usage = 
{"\n\n\
----------------------------------------------------------------------\n\
	This driver program calculates electrostatic potentials, energies,\n\
	and forces using both multigrid and finite element methods.\n\
        It is invoked as:\n\n\
	apbs [options] apbs.in\n\n\
	where apbs.in is a formatted input file and [options] are:\n\n\
--output-file=<name>     Enables output logging to the path\n\
	listed in <name>.  Uses flat-file\n\
	format is --output-format is not used.\n\
--output-format=<type>   Specifies format for logging.  Options\n\
	for type are either \"xml\" or \"flat\".\n\
--help                   Display this help information.\n\
--version                Display the current APBS version.\n\
----------------------------------------------------------------------\n\n"};

	/* ************** CHECK PARALLEL STATUS *************** */
	VASSERT(Vcom_init(&argc, &argv));
	com = Vcom_ctor(1);
	rank = Vcom_rank(com);
	size = Vcom_size(com);
	startVio(); 
	Vnm_setIoTag(rank, size);
	Vnm_tprint( 0, "Hello world from PE %d\n", rank);

	/* A bit of array/pointer initialization */
	mem = Vmem_ctor("MAIN");
	for (i=0; i<NOSH_MAXCALC; i++) {
		pmg[i] = VNULL;
		pmgp[i] = VNULL;
		fetk[i] = VNULL;
		pbe[i] = VNULL;
		qfEnergy[i] = 0;
		qmEnergy[i] = 0;
		dielEnergy[i] = 0;
		totEnergy[i] = 0;
		atomForce[i] = VNULL;
		nenergy[i] = 0;
		nforce[i] = 0;
	}
	for (i=0; i<NOSH_MAXMOL; i++) {
		alist[i] = VNULL;
		dielXMap[i] = VNULL;
		dielYMap[i] = VNULL;
		dielZMap[i] = VNULL;
		kappaMap[i] = VNULL;
		chargeMap[i] = VNULL;
	}

	/* ********* CHECK INVOCATION AND OPTIONS ************* */
	Vnm_tstart(APBS_TIMER_WALL_CLOCK, "APBS WALL CLOCK");
	Vnm_tprint( 1, "%s", header);

#ifdef APBS_FAST
	printf("WARNING: APBS was compiled with the --enable-fast option.\n"
		   "WARNING: This mode is experimental and subject to change in future releases.\n"
		   "WARNING: The fast mode enables: Gauess-Seidel Smoothing and \n"
		   "WARNING:   Conjugate Gradient Multigrid methods.\n\n");
#endif

	Vnm_tprint( 1, "This executable compiled on %s at %s\n\n", __DATE__, __TIME__);

	i=0;
	outputformat = OUTPUT_NULL;
	while (i<argc){
		if (strncmp(argv[i], "--", 2) == 0) {
			
			/* Long Options */
			if (Vstring_strcasecmp("--version", argv[i]) == 0){
				Vnm_tprint(2, "%s\n", PACKAGE_STRING);
				VJMPERR1(0);
			} else if (Vstring_strcasecmp("--help", argv[i]) == 0){
				Vnm_tprint(2, "%s\n", usage);
				VJMPERR1(0); 
			} else if (strncmp(argv[i], "--output-format", 15) == 0) {
				if (strstr(argv[i], "xml") != NULL) outputformat = OUTPUT_XML;
				else if (strstr(argv[i], "flat") != NULL) {
					outputformat = OUTPUT_FLAT;
				} else {
					Vnm_tprint(2, "Invalid output-format type!\n");
					VJMPERR1(0);
				}
			} else if (strncmp(argv[i], "--output-file=", 14) == 0){
				output_path = strstr(argv[i], "=");
				++output_path;
				if (outputformat == OUTPUT_NULL) outputformat = OUTPUT_FLAT;
			} else {
				Vnm_tprint(2, "UNRECOGNIZED COMMAND LINE OPTION %s!\n", argv[i]);
				Vnm_tprint(2, "%s\n", usage);
				VJMPERR1(0); 
			}  
		} else {
			
			/* Set the path to the input file */
			if ((input_path == VNULL) && (i != 0)) input_path = argv[i];
			else if (i != 0) {
				Vnm_tprint(2, "ERROR -- CALLED WITH TOO MANY ARGUMENTS!\n", \
						   argc);
				Vnm_tprint(2, "%s\n", usage);
				VJMPERR1(0);
			}
		}
		i++; 
	}

	if ((outputformat != 0) && (output_path == NULL)) {
		Vnm_tprint(2, "The --output-path variable must be set when using --output-format!\n");
		VJMPERR1(0);
	}

	if (input_path == NULL) {
		Vnm_tprint(2, "ERROR -- APBS input file not specified!\n", argc);
		Vnm_tprint(2, "%s\n", usage);
		VJMPERR1(0);
	} 

	/* Append rank info if a parallel run */
	if ((size > 1) && (output_path != NULL))
		printf(output_path, "%s_%d", output_path, rank);

	/* *************** PARSE INPUT FILE ******************* */
	nosh = NOsh_ctor(rank, size);
	Vnm_tprint( 1, "Parsing input file %s...\n", input_path);
	sock = Vio_ctor("FILE", "ASC", VNULL, input_path, "r");
	if (sock == VNULL) {
		Vnm_tprint(2, "Error while opening input file %s!\n", input_path);
		VJMPERR1(0);
	}
	if (!NOsh_parseInput(nosh, sock)) {
		Vnm_tprint( 2, "Error while parsing input file.\n");
		VJMPERR1(0);
	} else Vnm_tprint( 1, "Parsed input file.\n");
	Vio_dtor(&sock);
	
	/* *************** LOAD PARAMETERS AND MOLECULES ******************* */	
	param = loadParameter(nosh);
	if (loadMolecules(nosh, param, alist) != 1) {
		Vnm_tprint(2, "Error reading molecules!\n");
		VJMPERR1(0);
	}

	/* *************** SETUP CALCULATIONS *************** */
	if (NOsh_setupElecCalc(nosh, alist) != 1) {
		Vnm_tprint(2, "Error setting up ELEC calculations\n");
		VJMPERR1(0);
	}
	
	if ((rc = NOsh_setupApolCalc(nosh, alist)) == ACD_ERROR) {
		Vnm_tprint(2, "Error setting up APOL calculations\n");
		VJMPERR1(0);
	}
	
	/* ******************* CHECK APOL********************** */
	//if((nosh->gotparm == 0) && (rc == ACD_YES)){
	//	Vnm_print(1,"\nError you must provide a parameter file if you\n" \
	//				"     are performing an APOLAR calculation\n");
	//	VJMPERR1(0);
	//}

	/* *************** LOAD MAPS ******************* */
	if (loadDielMaps(nosh, dielXMap, dielYMap, dielZMap) != 1) {
		Vnm_tprint(2, "Error reading dielectric maps!\n");
		VJMPERR1(0);
	}
	if (loadKappaMaps(nosh, kappaMap) != 1) {
		Vnm_tprint(2, "Error reading kappa maps!\n");
		VJMPERR1(0);
	}
	if (loadChargeMaps(nosh, chargeMap) != 1) {
		Vnm_tprint(2, "Error reading charge maps!\n");
		VJMPERR1(0);
	}

	/* *************** DO THE CALCULATIONS ******************* */
	Vnm_tprint( 1, "Preparing to run %d PBE calculations.\n",
				nosh->ncalc);
	for (i=0; i<nosh->ncalc; i++) {
		Vnm_tprint( 1, "----------------------------------------\n");
		
		switch (nosh->calc[i]->calctype) {  
			case NCT_MG:
				/* What is this?  This seems like a very awkward way to find 
				the right ELEC statement... */
				for (k=0; k<nosh->nelec; k++) {
					if (nosh->elec2calc[k] >= i) {
						break;
					}
				}
				if (Vstring_strcasecmp(nosh->elecname[k], "") == 0) {
					Vnm_tprint( 1, "CALCULATION #%d: MULTIGRID\n", i+1);
				} else {
					Vnm_tprint( 1, "CALCULATION #%d (%s): MULTIGRID\n", 
								i+1, nosh->elecname[k]);
				}
				/* Useful local variables */
				mgparm = nosh->calc[i]->mgparm;
				pbeparm = nosh->calc[i]->pbeparm;
				
				/* Set up problem */
				Vnm_tprint( 1, "  Setting up problem...\n");
				if (!initMG(i, nosh, mgparm, pbeparm, realCenter, pbe, 
							alist, dielXMap, dielYMap, dielZMap, kappaMap, chargeMap, 
							pmgp, pmg)) {
					Vnm_tprint( 2, "Error setting up MG calculation!\n");
					VJMPERR1(0);
				}
					
				/* Print problem parameters */
				printMGPARM(mgparm, realCenter);
				printPBEPARM(pbeparm);
				
				/* Solve PDE */
				if (solveMG(nosh, pmg[i], mgparm->type) != 1) {
					Vnm_tprint(2, "Error solving PDE!\n");
					VJMPERR1(0);
				}
					
				/* Set partition information for observables and I/O */
				if (setPartMG(nosh, mgparm, pmg[i]) != 1) {
					Vnm_tprint(2, "Error setting partition info!\n");
					VJMPERR1(0);
				}
					
				/* Write out energies */
				energyMG(nosh, i, pmg[i], 
						&(nenergy[i]), &(totEnergy[i]), &(qfEnergy[i]), 
						&(qmEnergy[i]), &(dielEnergy[i]));
				
				/* Write out forces */
				forceMG(mem, nosh, pbeparm, mgparm, pmg[i], &(nforce[i]), 
						&(atomForce[i]), alist);
				
				/* Write out data folks might want */
				writedataMG(rank, nosh, pbeparm, pmg[i]);
				
				/* Write matrix */
				writematMG(rank, nosh, pbeparm, pmg[i]);
				
				/* If needed, cache atom energies */				
				nenergy[i] = 0;
				if ((pbeparm->calcenergy == PCE_COMPS) && (outputformat != OUTPUT_NULL)){
					storeAtomEnergy(pmg[i], i, &(atomEnergy[i]), &(nenergy[i]));
				}
					
				fflush(stdout);
				fflush(stderr);
				
				break;
				
				/* ***** Do FEM calculation ***** */
			case NCT_FEM:
#ifdef HAVE_MC_H
				for (k=0; k<nosh->nelec; k++) {
					if (nosh->elec2calc[k] >= i) break;
				}
				if (Vstring_strcasecmp(nosh->elecname[i+1], "") == 0) {
					Vnm_tprint( 1, "CALCULATION #%d: FINITE ELEMENT\n", i+1);
				} else {
					Vnm_tprint( 1, "CALCULATION #%d (%s): FINITE ELEMENT\n", i+1, nosh->elecname[k+1]);
				}
				
				/* Useful local variables */
				feparm = nosh->calc[i]->femparm;
				pbeparm = nosh->calc[i]->pbeparm;
				
				/* Warn the user about some things */
				Vnm_tprint(2, "#################### WARNING ###################\n");
				Vnm_tprint(2, "## FE support is currently very experimental! ##\n");
				Vnm_tprint(2, "#################### WARNING ###################\n");
				
				/* Set up problem */
				Vnm_tprint( 1, "  Setting up problem...\n");
				if (!initFE(i, nosh, feparm, pbeparm, pbe, alist, fetk)) {
					Vnm_tprint( 2, "Error setting up FE calculation!\n");
					VJMPERR1(0);
				}
					
					/* Print problem parameters */
					printFEPARM(i, nosh, feparm, fetk);
				printPBEPARM(pbeparm);
				
				/* Refine mesh */
				if (!preRefineFE(i, nosh, feparm, fetk)) {
					Vnm_tprint( 2, "Error pre-refining mesh!\n");
					VJMPERR1(0);
				}
					
				/* Solve-estimate-refine */
				Vnm_tprint(2, "\n\nWARNING!  DO NOT EXPECT PERFORMANCE OUT OF THE APBS/FEtk\n");
				Vnm_tprint(2, "INTERFACE AT THIS TIME.  THE FINITE ELEMENT SOLVER IS\n");
				Vnm_tprint(2, "CURRENTLY NOT OPTIMIZED FOR THE PB EQUATION.  IF YOU WANT\n");
				Vnm_tprint(2, "PERFORMANCE, PLEASE USE THE MULTIGRID-BASED METHODS, E.G.\n");
				Vnm_tprint(2, "MG-AUTO, MG-PARA, and MG-MANUAL (SEE DOCS.)\n\n");
				Vnm_tprint(1, "  Beginning solve-estimate-refine cycle:\n");
				for (isolve=0; isolve<feparm->maxsolve; isolve++) {
					Vnm_tprint(1, "    Solve #%d...\n", isolve);
					if (!solveFE(i, nosh, pbeparm, feparm, fetk)) {
						Vnm_tprint(2, "ERROR SOLVING EQUATION!\n");
						VJMPERR1(0);
					}
					if (!energyFE(nosh, i, fetk, &(nenergy[i]), 
								  &(totEnergy[i]), &(qfEnergy[i]), 
								  &(qmEnergy[i]), &(dielEnergy[i]))) {
						Vnm_tprint(2, "ERROR SOLVING EQUATION!\n");
						VJMPERR1(0);
					}
					/* We're not going to refine if we've hit the max number
						* of solves */
					if (isolve < (feparm->maxsolve)-1) {
						if (!postRefineFE(i, nosh, feparm, fetk)) break;
					}
					bytesTotal = Vmem_bytesTotal();
					highWater = Vmem_highWaterTotal();
					Vnm_tprint(1, "      Currently memory use:  %g MB\n", 
							   ((double)bytesTotal/(1024.)/(1024.)));
					Vnm_tprint(1, "      High-water memory use:  %g MB\n", 
							   ((double)highWater/(1024.)/(1024.)));
				}
					
					Vnm_tprint(1, "  Writing FEM data to files.\n");
				if (!writedataFE(rank, nosh, pbeparm, fetk[i])) {
					Vnm_tprint(2, "  Error while writing FEM data!\n");
				}
#else /* ifdef HAVE_MC_H */
					Vnm_print(2, "Error!  APBS not compiled with FEtk!\n");
				exit(2);
#endif /* ifdef HAVE_MC_H */
				break;
				
			/* Do an apolar calculation */
			case NCT_APOL:
				/* Copied from NCT_MG. See the note above (top of loop) for
					information about this loop.
				*/
				for (k=0; k<nosh->napol; k++) {
					if (nosh->apol2calc[k] >= i) {
						break;
					}
				}
				
				if (Vstring_strcasecmp(nosh->apolname[k], "") == 0) {
					Vnm_tprint( 1, "CALCULATION #%d: APOLAR\n", i+1);
				} else {
					Vnm_tprint( 1, "CALCULATION #%d (%s): APOLAR\n", 
								i+1, nosh->apolname[k]);
				}

				apolparm = nosh->calc[i]->apolparm;
				rc = initAPOL(nosh, mem, param, apolparm, &(nforce[i]), &(atomForce[i]), 
						 alist[(apolparm->molid)-1]);
				if(rc == 0) {
					Vnm_tprint(2, "Error calculating apolar solvation quantities!\n");
					VJMPERR1(0);
				}
				break;
			default:
				Vnm_tprint(2, "  Unknown calculation type (%d)!\n", 
						   nosh->calc[i]->calctype);
				exit(2);
		}
	}
	
	//Clear out the parameter file memory
	if(param != VNULL) Vparam_dtor(&param);
	
	/* *************** HANDLE PRINT STATEMENTS ******************* */
	if (nosh->nprint > 0) {
		Vnm_tprint( 1, "----------------------------------------\n");
		Vnm_tprint( 1, "PRINT STATEMENTS\n");
	}
	for (i=0; i<nosh->nprint; i++) {
		/* Print energy */
		if (nosh->printwhat[i] == NPT_ENERGY) {
			printEnergy(com, nosh, totEnergy, i);
			/* Print force */
		} else if (nosh->printwhat[i] == NPT_FORCE) {
			printForce(com, nosh, nforce, atomForce, i);
		} else if (nosh->printwhat[i] == NPT_ELECENERGY) {
			printElecEnergy(com, nosh, totEnergy, i);
		} else if (nosh->printwhat[i] == NPT_ELECFORCE) {
			printElecForce(com, nosh, nforce, atomForce, i);
		} else if (nosh->printwhat[i] == NPT_APOLENERGY) {
			printApolEnergy(nosh, i);
		} else if (nosh->printwhat[i] == NPT_APOLFORCE) {
			printApolForce(com, nosh, nforce, atomForce, i);
		} else {
			Vnm_tprint( 2, "Undefined PRINT keyword!\n");
			break;
		}
	} 
	Vnm_tprint( 1, "----------------------------------------\n");
	
	/* *************** HANDLE LOGGING *********************** */

	if (outputformat == OUTPUT_XML) {
		Vnm_tprint(2, "  Writing data to XML file %s...\n\n", output_path);
		writedataXML(nosh, com, output_path, totEnergy, qfEnergy, qmEnergy,
					 dielEnergy, nenergy, atomEnergy, nforce, atomForce);
		
	} else if (outputformat == OUTPUT_FLAT) {
		Vnm_tprint(2," Writing data to flat file %s...\n\n", output_path);
		writedataFlat(nosh, com, output_path, totEnergy, qfEnergy, qmEnergy,
					  dielEnergy, nenergy, atomEnergy, nforce, atomForce);
	}

	/* Destroy energy arrays if they still exist */

	for (i=0; i<nosh->ncalc; i++) {
		if (nenergy[i] > 0) Vmem_free(mem, nenergy[i], sizeof(double),
									  (void **)&(atomEnergy[i]));    
	}

	/* *************** GARBAGE COLLECTION ******************* */

	Vnm_tprint( 1, "CLEANING UP AND SHUTTING DOWN...\n");
	/* Clean up APBS structures */
	killForce(mem, nosh, nforce, atomForce);
	killEnergy();
	killMG(nosh, pbe, pmgp, pmg);
	killChargeMaps(nosh, chargeMap);
	killKappaMaps(nosh, kappaMap);
	killDielMaps(nosh, dielXMap, dielYMap, dielZMap);
	killMolecules(nosh, alist);
	NOsh_dtor(&nosh);

	/* Memory statistics */
	bytesTotal = Vmem_bytesTotal();
	highWater = Vmem_highWaterTotal();
	Vnm_tprint( 1, "Final memory usage:  %4.3f MB total, %4.3f MB high water\n", 
				(double)(bytesTotal)/(1024.*1024.), 
				(double)(highWater)/(1024.*1024.));

	/* Clean up MALOC structures */
	Vcom_dtor(&com);
	Vmem_dtor(&mem);

	/* And now it's time to so "so long"... */
	Vnm_tprint(1, "\n\n");
	Vnm_tprint( 1, "Thanks for using APBS!\n\n");

	/* This should be last */
	Vnm_tstop(APBS_TIMER_WALL_CLOCK, "APBS WALL CLOCK");
	Vnm_flush(1);
	Vnm_flush(2);
	Vcom_finalize();

	fflush(NULL);
	
	return 0;

	VERROR1:
	Vcom_finalize();
	Vcom_dtor(&com);
	Vmem_dtor(&mem);
	return APBSRC;
}
