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
//////////////////////////////////////////////////////////////////////////// 
/// rcsid="$Id$"
//////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// File:     nosh.h    
//
// Purpose:  No shell class (i.e., fixed format input files)
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */

#ifndef _NOSH_H_
#define _NOSH_H_

#define NOSH_MAXMOL 20
#define NOSH_MAXPATH 256
#define NOSH_MAXMGPARM 20
#define NOSH_MAXFEMPARM 20
#define NOSH_MAXCALC 20
#define NOSH_MAXPRINT 20
#define NOSH_MAXPOP 20

#include "apbs/apbs.h"
#include "maloc/maloc.h"

/* ///////////////////////////////////////////////////////////////////////////
// Class NOsh_femparm: Definition
/////////////////////////////////////////////////////////////////////////// */
typedef struct NOsh_femparm {

   int dummy;

} NOsh_femparm;

/* ///////////////////////////////////////////////////////////////////////////
// Class NOsh_mgparm: Definition
/////////////////////////////////////////////////////////////////////////// */
typedef struct NOsh_mgparm {
   
    int dime[3];               /* Grid dimensions */
    int setdime;
    int nlev;                  /* Levels in multigrid hierarchy */
    int setnlev;
    double grid[3];            /* Grid spacings */
    int setgrid;
    double glen[3];            /* Grid side lengths */
    int setglen;
    int cmeth;                 /* Centering method: 0 => center on point, 
                                * 1 => center on molecule */
    double center[3];          /* Grid center */
    int setgcent;  
    int centmol;               /* Particular molecule on which we want to
                                * center the grid */
    int molid;                 /* Molecule ID to perform calculation on */
    int setmolid;
    int nonlin;                /* 0 => LPBE, 1 => NPBE */
    int setnonlin;
    int bcfl;                  /* Boundary condition: 0 => zero, 1 => single
                                * Debye-Huckel sphere, 2 => multiple Debye-
                                * Huckel spheres, 4 => focusing */
    int setbcfl;
    int nion;                  /* Number of counterion species */
    int setnion;
    double ionq[MAXION];  /* Counterion charges (in e) */
    double ionc[MAXION];  /* Counterion concentrations (in M) */
    double ionr[MAXION];  /* Counterion radii (in A) */
    int setion[MAXION];
    double pdie;               /* Solute dielectric */
    int setpdie;
    double sdie;               /* Solvent dielectric */
    int setsdie;
    int srfm;                  /* Surface calculation method
				* 0 => Mol surface for epsilon; inflated VdW
				* for kappa; no smoothing 
                                * 1 => As 0 with harmoic average
                                * smoothing
                                * 2 => Cubic spline */
    int setsrfm;
    double srad;               /* Solvent radius */
    int setsrad;
    double swin;               /* Cubic spline window */
    int setswin; 
    double temp;               /* Temperarture (in K) */
    int settemp;
    double gamma;              /* Surface tension for apolar energies/forces
                                * (in kJ/mol/A^2) */
    int setgamma;
    int calcenergy;            /* Energy calculation
				* 0 => don't calculate out energy
                                * 1 => calculate total energy 
                                * 2 => calculate total energy and all energy
                                *      components*/
    int setcalcenergy;      
    int calcforce;             /* Atomic forces I/O 
                                * 0 => don't calculate forces
                                * 1 => calculate net forces on molecule
                                * 2 => calculate atom-level forces */
    int setcalcforce;       
    int writepot;              /* 0 => no, 1 => yes */
    int setwritepot; 
    char writepotstem[NOSH_MAXPATH];    /* File stem to write pot */
    int writepotfmt;           /* Potential file formats: 0 => dx, 1 => avs, 
                                * 2 => UHBD */
    int writeacc;              /* 0 => no, 1 => yes */
    int setwriteacc; 
    char writeaccstem[256];    /* File stem to write pot */
    int writeaccfmt;           /* Potential file formats: 0 => dx, 1 => avs, 
                                * 2 => UHBD */
} NOsh_mgparm;

/* ///////////////////////////////////////////////////////////////////////////
// Class NOsh: Definition
/////////////////////////////////////////////////////////////////////////// */

typedef struct NOsh {

    NOsh_mgparm mgparm[NOSH_MAXMGPARM];       /* Parameter objects for MG 
                                               * types of calculations */
    NOsh_femparm femparm[NOSH_MAXFEMPARM];    /* Parameter objects for FEM 
                                               * types of calculations */
    int ncalc;                                /* The number of calculations to 
                                               * be done */
    int imgcalc, ifemcalc, nmgcalc, nfemcalc; /* Counters for the various MG 
                                               * and FEM calculations */
    int calctype[NOSH_MAXCALC];               /* The list of calculations: 0 =>
                                               * multigrid and 1 => FEM */
    int nmol;                                 /* Number of molecules */
    char molpath[NOSH_MAXMOL][NOSH_MAXPATH];  /* Paths to mol files */
    int nprint;                               /* How many print sections? */
    int printwhat[NOSH_MAXPRINT];
                                              /* What do we print (0=>energy) */
    int printnarg[NOSH_MAXPRINT];             /* How many arguments in energy 
                                               * list */
    int printcalc[NOSH_MAXPRINT][NOSH_MAXPOP];/* Calculation id */
    int printop[NOSH_MAXPRINT][NOSH_MAXPOP];  /* Operation id (0 = add, 1 = 
                                               * subtract) */
  int parsed;                                 /* Have we parsed an input file
                                               * yet? */

} NOsh;

/* ///////////////////////////////////////////////////////////////////////////
// Class NOsh: Non-inlineable methods (mcsh.c)
/////////////////////////////////////////////////////////////////////////// */

VEXTERNC NOsh* NOsh_ctor();
VEXTERNC int   NOsh_ctor2(NOsh *thee);
VEXTERNC void  NOsh_dtor(NOsh **thee);
VEXTERNC void  NOsh_dtor2(NOsh *thee);
VEXTERNC int   NOsh_parse(NOsh *thee, Vio *sock);

#endif 

