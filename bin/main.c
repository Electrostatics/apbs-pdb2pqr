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
// File:     main.c
//
// Purpose:  APBS ``front end" using formatted input files
//
// Arguments:    
//           This driver program represents a mish-mash of instructions for
//           calculating electrostatic potentials, as well as free energies of 
//           binding and solvation.  It is invoked as:
//
//               apbs apbs.in
//
//           where apbs.in is a formatted input file (see example).
//
// Notes:    Doesn't have parallel implementation yet
//
// Returns:  APBSRC (defined below) on failure, 0 otherwise
//
// Author:   Nathan Baker
//
// rcsid="$Id$"
/////////////////////////////////////////////////////////////////////////// */

#include "apbscfg.h"
#include "apbs/apbs.h"  
#include "apbs/nosh.h"  
#include "apbs/mgparm.h"  
#include "apbs/pbeparm.h"  
#include "apbs/femparm.h"  

#include "routines.h"

VEMBED(rcsid="$Id$")

int main(int argc, char **argv) {

    NOsh *nosh = VNULL;
    MGparm *mgparm = VNULL;
    PBEparm *pbeparm = VNULL;
    FEMparm *femparm = VNULL;
    Vmem *mem = VNULL;
    Vcom *com = VNULL;
    Vio *sock = VNULL;
    Vpmg *pmg[NOSH_MAXCALC];
    Vpmgp *pmgp[NOSH_MAXCALC];
    Vpbe *pbe[NOSH_MAXCALC];
    Valist *alist[NOSH_MAXMOL];
    char *input_path = VNULL;
    char outpath[VMAX_ARGLEN];
    double iparm, sparm;
    int i, j, rank, size;

    /* These variables require some explaining... The energy double arrays
     * store energies from the various calculations.  The energy int array
     * stores a flag (0, 1, 2) that indicates whether no, total, or all
     * energies were evaluated for a particular calculation.  Likewise, the
     * force double arrays store forces from the various calcualtions.  The
     * force int array stores an integer which either says no calculation was
     * performed (0) or gives the number of entries in the force array for each
     * calculation */
    double qfEnergy[NOSH_MAXCALC], qmEnergy[NOSH_MAXCALC];
    double dielEnergy[NOSH_MAXCALC], totEnergy[NOSH_MAXCALC];
    double partMin[3], partMax[3];
    AtomForce *atomForce[NOSH_MAXCALC];
    double ibForce[3], qfForce[3], dbForce[3], npForce[3], tenergy;
    int nenergy[NOSH_MAXCALC], nforce[NOSH_MAXCALC], bytesTotal, highWater;
    int calcid;
    /* THe real partition centers */
    double realCenter[3];

    /* Instructions: */
    char *header = "\n\n\
    ----------------------------------------------------------------------\n\
    Adaptive Poisson-Boltzmann Solver (APBS)\n\
    Version 0.1.6 (November 17, 2001)\n\
    \n\
    Nathan A. Baker (nbaker@wasabi.ucsd.edu)\n\
    Dept. of Chemistry and Biochemistry\n\
    Dept. of Mathematics, Scientific Computing Group\n\
    University of California, San Diego \n\n\
    Additional contributing authors listed in the code documentation.\n\n\
    Copyright (c) 1999-2001.\n\
    The Regents of the University of California (Regents).\n\
    All Rights Reserved.\n\n\
    Permission to use, copy, modify, and distribute this software and its\n\
    documentation for educational, research, and not-for-profit purposes,\n\
    without fee and without a signed licensing agreement, is hereby granted,\n\
    provided that the above copyright notice, this paragraph and the\n\
    following two paragraphs appear in all copies, modifications, and\n\
    distributions.\n\n\
    IN NO EVENT SHALL REGENTS OR THE AUTHORS BE LIABLE TO ANY PARTY FOR\n\
    DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES,\n\
    INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS\n\
    DOCUMENTATION, EVEN IF REGENTS OR THE AUTHORS HAVE BEEN ADVISED OF THE\n\
    POSSIBILITY OF SUCH DAMAGE.\n\n\
    REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT\n\
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A\n\
    PARTICULAR PURPOSE.  THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF\n\
    ANY, PROVIDED HEREUNDER IS PROVIDED \"AS IS\".  REGENTS HAS NO OBLIGATION\n\
    TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR\n\
    MODIFICATIONS.\n\
    ----------------------------------------------------------------------\n\
    \n\n";
    char *usage = "\n\n\
    ----------------------------------------------------------------------\n\
    This driver program calculates electrostatic potentials, energies,\n\
    and forces using both multigrid and finite element methods.\n\
    It is invoked as:\n\n\
      apbs apbs.in\n\n\
    where apbs.in is a formatted input file.\n\
    ----------------------------------------------------------------------\n\n";


    /* ************** CHECK PARALLEL STATUS *************** */
    Vcom_init(&argc, &argv);
    com = Vcom_ctor(1);
    rank = Vcom_rank(com);
    size = Vcom_size(com);
    Vio_start();
    Vnm_setIoTag(rank, size);
    Vnm_tprint( 0, "main:  Hello world from PE %d\n", rank);

    /* A bit of array/pointer initialization */
    mem = Vmem_ctor("MAIN");
    for (i=0; i<NOSH_MAXCALC; i++) {
        atomForce[i] = VNULL;
        nenergy[i] = 0;
        nforce[i] = 0;
    }

    /* *************** CHECK INVOCATION ******************* */
    Vnm_tstart(26, "APBS WALL CLOCK");
    Vnm_tprint( 1, "%s", header);
    Vnm_tprint( 1, "This executable compiled on %s at %s\n\n", __DATE__, 
      __TIME__);
    if (argc != 2) {
        Vnm_tprint( 2,"%s\n", usage);
        return APBSRC;
    } 
    input_path = argv[1];


    /* *************** PARSE INPUT FILE ******************* */
    nosh = NOsh_ctor(com);
    sock = Vio_ctor("FILE", "ASC", VNULL, input_path, "r");
    Vnm_tprint( 1, "main:  Parsing input file %s...\n", input_path);
    if (!NOsh_parse(nosh, sock)) {
        Vnm_tprint( 2, "main:  Error while parsing input file.\n");
        return APBSRC;
    } else Vnm_tprint( 1, "main:  Parsed input file.\n");
    Vio_dtor(&sock);

    /* *************** LOAD MOLECULES ******************* */
    if (loadMolecules(com, nosh, alist) != 1) {
        Vnm_tprint(2, "main:  Error reading molecules!\n");
        return APBSRC;
    }

    /* *************** DO THE CALCULATIONS ******************* */
    Vnm_tprint( 1, "main:  Preparing to run %d PBE calculations.\n",
      nosh->ncalc);
    for (i=0; i<nosh->ncalc; i++) {
        Vnm_tprint( 1, "main:  ----------------------------------------\n");

        /* ***** Do MG calculation ***** */
        if (nosh->calc[i].calctype == 0) {

            Vnm_tprint( 1, "main:  CALCULATION #%d: MULTIGRID\n", i+1);

            /* Useful local variables */
            mgparm = nosh->calc[i].mgparm;
            pbeparm = nosh->calc[i].pbeparm;

            /* Set up problem */
            Vnm_tprint( 1, "main:    Setting up problem...\n");
            if (!initMG(com, i, nosh, mgparm, pbeparm, realCenter, pbe, 
              alist, pmgp, pmg)) {
                Vnm_tprint( 2, "main:  Error setting up MG calculation!\n");
                return APBSRC;
            }

            /* Print problem parameters */
            printMGPARM(com, mgparm, realCenter);
            printPBEPARM(com, pbeparm);

            /* Solve PDE */
            if (solveMG(com, pmg[i]) != 1) {
                Vnm_tprint(2, "main:  Error solving PDE!\n");
                return APBSRC;
            }

            /* Set partition information for observables and I/O */
            if (setPartMG(com, mgparm, pmg[i]) != 1) {
                Vnm_tprint(2, "main:  Error setting partition info!\n");
                return APBSRC;
            }

            /* Write out energies */
            energyMG(com, nosh, i, pmg[i], &(nenergy[i]), 
              &(totEnergy[i]), &(qfEnergy[i]), &(qmEnergy[i]), 
              &(dielEnergy[i]));

            /* Write out forces */
            forceMG(com, mem, nosh, pbeparm, pmg[i], &(nforce[i]), 
              &(atomForce[i]), alist);

            /* Write out potential */
            writepotMG(com, pbeparm, pmg[i]);
            
            /* Write accessibility */
            writeaccMG(com, pbeparm, pmg[i]);

            fflush(stdout);
            fflush(stderr);

        /* ***** Do FEM calculation ***** */
        } else {
            Vnm_tprint( 2, "main: FEM shell support not implemented yet\n");
            return APBSRC;
        }
    } 

    
    /* *************** HANDLE PRINT STATEMENTS ******************* */
    if (nosh->nprint > 0) {
        Vnm_tprint( 1, "main:  ----------------------------------------\n");
        Vnm_tprint( 1, "main:  PRINT STATEMENTS\n");
    }
    for (i=0; i<nosh->nprint; i++) {
        /* Print energy */
        if (nosh->printwhat[i] == 0) {
            printEnergy(com, nosh, totEnergy, i);
        } else {
            Vnm_tprint( 2, "main:  Undefined PRINT keyword!\n");
            break;
        }
    }
 

    /* *************** GARBAGE COLLECTION ******************* */
    Vnm_tprint( 1, "main:  ----------------------------------------\n");
    Vnm_tprint( 1, "main:  CLEANING UP AND SHUTTING DOWN...\n");
    for (i=0; i<nosh->ncalc; i++) {
        if (nforce[i] > 0) 
          Vmem_free(mem, nforce[i], sizeof(AtomForce), 
            (void **)&(atomForce[i])); 
    }
    for (i=0; i<nosh->nmol; i++) Valist_dtor(&(alist[i]));
    NOsh_dtor(&nosh);
    Vmem_dtor(&mem);

    /* Stop wall clock timer */
    Vnm_tstop(26, "APBS WALL CLOCK");

    Vnm_print(1, "\n\n");
    Vnm_tprint( 1, "Thanks for using APBS!\n\n");

    Vcom_finalize();
    Vcom_dtor(&com);

    return 0;

}
