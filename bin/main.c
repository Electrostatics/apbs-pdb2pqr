/**
 *  @file    main.c
 *  @ingroup Frontend
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
 * Copyright (c) 2002-2010, Washington University in St. Louis.
 * Portions Copyright (c) 2002-2010.  Nathan A. Baker
 * Portions Copyright (c) 1999-2002.  The Regents of the University of California.
 * Portions Copyright (c) 1995.  Michael Holst
 *
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met: 
 *
 * -  Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.  
 * 
 * - Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 * 
 * - Neither the name of Washington University in St. Louis nor the names of its
 * contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
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

/** 
 * @brief The main APBS function
 * @ingroup  Frontend
 * @author  Nathan Baker, Dave Gohara, Todd Dolinsky
 * @returns Status code (0 for success)
 */
int main(
		 int argc,  /**< Number of arguments */
		 char **argv  /**< Argument strings */
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
	Gem *gm[NOSH_MAXMOL];
#else
	void *fetk[NOSH_MAXCALC];
	void *gm[NOSH_MAXMOL];
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
	Version 1.2.1\n\
	\n\
	Nathan A. Baker (baker@biochem.wustl.edu)\n\
	Dept. Biochemistry and Molecular Biophysics\n\
	Center for Computational Biology\n\
	Washington University in St. Louis\n\
	\n\
	Additional contributing authors listed in the code documentation.\n\
	\n\
	Copyright (c) 2002-2010, Washington University in St. Louis.\n\
	Portions Copyright (c) 2002-2010.  Nathan A. Baker\n\
	Portions Copyright (c) 1999-2002.  The Regents of the University of California.\n\
	Portions Copyright (c) 1995.  Michael Holst\n\
	\n\
	All rights reserved.\n\
	\n\
	Redistribution and use in source and binary forms, with or without\n\
	modification, are permitted provided that the following conditions are met: \n\
	\n\
	* Redistributions of source code must retain the above copyright notice, this\n\
	list of conditions and the following disclaimer.  \n\
	\n\
	* Redistributions in binary form must reproduce the above copyright notice,\n\
	this list of conditions and the following disclaimer in the documentation\n\
	and/or other materials provided with the distribution.\n\
	\n\
	* Neither the name of Washington University in St. Louis nor the names of its\n\
	contributors may be used to endorse or promote products derived from this\n\
	software without specific prior written permission.\n\
	\n\
	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS\n\
	\"AS IS\" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT\n\
	LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR\n\
	A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR\n\
	CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,\n\
	EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,\n\
	PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR\n\
	PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF\n\
	LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING\n\
	NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS\n\
	SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n\
----------------------------------------------------------------------\n\
	APBS uses FETK (the Finite Element ToolKit) to solve the\n\
	Poisson-Boltzmann equation numerically.  FETK is a portable collection\n\
	of finite element modeling class libraries developed by the Michael Holst\n\
	research group and written in an object-oriented form of C.  FEtk is\n\
	designed to solve general coupled systems of nonlinear partial differential\n\
	equations using adaptive finite element methods, inexact Newton methods,\n\
	and algebraic multilevel methods.  More information about FEtk may be found\n\
	at <http://www.FEtk.ORG>.\n\
----------------------------------------------------------------------\n\
	APBS also uses Aqua to solve the Poisson-Boltzmann equation numerically.  \n\
	Aqua is a modified form of the Holst group PMG library <http://www.FEtk.ORG>\n\
	which has been modified by Patrice Koehl\n\
	<http://koehllab.genomecenter.ucdavis.edu/> for improved efficiency and\n\
	memory usage when solving the Poisson-Boltzmann equation.\n\
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
				if (strstr(argv[i], "xml") != NULL) {
					Vnm_tprint(2, "XML output format is now deprecated, please use --output-format=flat instead!\n\n");
					VJMPERR1(0); 
				}
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
	/* if((nosh->gotparm == 0) && (rc == ACD_YES)){
	 	Vnm_print(1,"\nError you must provide a parameter file if you\n" \
	 				"     are performing an APOLAR calculation\n");
	 	VJMPERR1(0);
	} */
	
#if defined(DEBUG_MAC_OSX_OCL)
#include "mach_chud.h"
#include <stdint.h>
	uint64_t mbeg;
	machm_(&mbeg);
	
	if(clFinish != NULL)
	{
		int ret = initOpenCL();
        printf("OpenCL runtime present - initialized = %i\n",ret);
    }
    else
	{
		setkOpenCLAvailable_(0);
        printf("OpenCL is not present!\n");
	}
#endif
	
#if defined(DEBUG_MAC_OSX_STANDARD)
#include "mach_chud.h"
#include <stdint.h>
	uint64_t mbeg;
	machm_(&mbeg);
#endif
	
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
				
#if defined(INCLUDE_MULTI)
				/* HACK: Write multivalue values for each atom */
				writeMultivalue(pbeparm,pmg[i],alist[i]);
				/* HACK: End hack */
#endif					
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
				if (initFE(i, nosh, feparm, pbeparm, pbe, alist, fetk, gm) != VRC_SUCCESS) {
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

	if (outputformat == OUTPUT_FLAT) {
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
#ifdef HAVE_MC_H
	killFE(nosh, pbe, fetk, gm);
#endif
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

#if defined(DEBUG_MAC_OSX_OCL)
	mets_(&mbeg, "Main Program CL");
#endif
#if defined(DEBUG_MAC_OSX_STANDARD)
	mets_(&mbeg, "Main Program Standard");
#endif
	
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
