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
#define NOSH_MAXCALC 20
#define NOSH_MAXPRINT 20
#define NOSH_MAXPOP 20

#include "apbs/apbs.h"
#include "maloc/maloc.h"
#include "apbs/femparm.h"
#include "apbs/mgparm.h"
#include "apbs/pbeparm.h"

/* ///////////////////////////////////////////////////////////////////////////
// Class NOsh_calc: A small class which allows NOsh to keep track of various
//                  calculations in a relatively straightforward way.
/////////////////////////////////////////////////////////////////////////// */
typedef struct NOsh_calc {

    MGparm *mgparm;                      /* Multigrid parameters */
    FEMparm *femparm;                    /* Finite element parameters */
    PBEparm *pbeparm;                    /* Generic PBE parameters */
    int calctype;                        /* 0 => multigrid, 1 => FEM */

} NOsh_calc;


/* ///////////////////////////////////////////////////////////////////////////
// Class NOsh: Definition
/////////////////////////////////////////////////////////////////////////// */

typedef struct NOsh {

    NOsh_calc calc[NOSH_MAXCALC];        /* The array of calculation objects */

    int ncalc;                           /* The number of calculations in the
                                          * calc array */
    int nelec;
    int ispara;                          /* (1 => is a parallel calculation, 
                                          *  0 => is not) */
    Vcom *com;                           /* Communications object for parallel
                                          * focusing calculations */
    int bogus;                           /* A flag which tells routines using
                                          * NOsh that this particular NOsh is
                                          * broken -- useful for parallel
                                          * focusing calculations where the
                                          * user gave us too many processors 
                                          * (1 => ignore this NOsh; 
                                          *  0 => this NOsh is OK) */
    int elec2calc[NOSH_MAXCALC];         /* A mapping between ELEC statements
					  * which appear in the input file and
					  * calc objects stored above.  Since
					  * we allow both normal and focused
					  * multigrid, there isn't a 1-to-1
					  * correspondence between ELEC
					  * statements and actual calcualtions.
					  * This can really confuse operations
					  * which work on specific calculations
					  * further down the road (like PRINT).
					  * Therefore this array is the initial
					  * point of entry for any
					  * calculation-specific operation.  It
					  * points to a specific entry in the
					  * calc array. */
    int nmol;                            /* Number of molecules */
    char molpath[NOSH_MAXMOL][VMAX_ARGLEN];   
                                         /* Paths to mol files */
    int nprint;                          /* How many print sections? */
    int printwhat[NOSH_MAXPRINT];        /* What do we print (0=>energy) */
    int printnarg[NOSH_MAXPRINT];        /* How many arguments in energy 
                                          * list */
    int printcalc[NOSH_MAXPRINT][NOSH_MAXPOP];
                                         /* ELEC id (see elec2calc) */
    int printop[NOSH_MAXPRINT][NOSH_MAXPOP];  
                                         /* Operation id (0 = add, 1 = 
                                          * subtract) */
  int parsed;                            /* Have we parsed an input file
                                          * yet? */

} NOsh;

/* ///////////////////////////////////////////////////////////////////////////
// Class NOsh: Non-inlineable methods (mcsh.c)
/////////////////////////////////////////////////////////////////////////// */

VEXTERNC NOsh* NOsh_ctor(Vcom *com);
VEXTERNC int   NOsh_ctor2(NOsh *thee, Vcom *com);
VEXTERNC void  NOsh_dtor(NOsh **thee);
VEXTERNC void  NOsh_dtor2(NOsh *thee);
VEXTERNC int   NOsh_parse(NOsh *thee, Vio *sock);

#endif 

