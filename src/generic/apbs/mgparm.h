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
// File:     mgparm.h    
//
// Purpose:  A set of useful parameters for a generic multigrid calculation
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */

#ifndef _MGPARM_H_
#define _MGPARM_H_

#include "apbs/apbs.h"
#include "maloc/maloc.h"

/* ///////////////////////////////////////////////////////////////////////////
// Class MGparm: Definition
//
// If you add/change something here, it must be added/changed in MGparm_copy
/////////////////////////////////////////////////////////////////////////// */
typedef struct MGparm {

    int type;                            /* What type of MG calculation?
                                          *    0 => sequential manual
                                          *    1 => sequential auto-focus
                                          *    2 => parallel auto-focus */
    int parsed;

    /* *** GENERIC PARAMETERS *** */
    int dime[3];               /* Grid dimensions */
    int setdime;

    /* *** TYPE 0 PARAMETERS (SEQUENTIAL MANUAL) *** */
    int nlev;                  /* DEPRECATED!!!! (Ignored now)
                                * Levels in multigrid hierarchy */
    int setnlev;
    double grid[3];            /* Grid spacings */
    int setgrid;
    double glen[3];            /* Grid side lengths. */
    int setglen;
    int cmeth;                 /* Centering method: 0 => center on point, 
                                * 1 => center on molecule */
    double center[3];          /* Grid center. If ispart = 0, then this is only
                                * meaningful if cmeth = 0.  However, if ispart
                                * = 1 and cmeth = 0, then this is the center of
                                * the non-disjoint (overlapping) partition.  If
                                * ispart = 1 and cmeth = 1, then this is the
                                * vector that must be added to the center of
                                * the molecule to give the center of the
                                * non-disjoint partition.  */
    int centmol;               /* Particular molecule on which we want to
                                * center the grid */
    int setgcent;  

    /* ******** TYPE 1 & 2 PARAMETERS (SEQUENTIAL & PARALLEL AUTO-FOCUS) *** */
    double cglen[3];           /* Coarse grid side lengths */
    int setcglen;
    double fglen[3];           /* Fine grid side lengths */
    int setfglen;
    int ccmeth;                /* Coarse grid centering method:  0 => center on
                                * point, 1 => center on molecule */
    double ccenter[3];         /* Coarse grid center.  */
    int ccentmol;              /* Particular molecule on which we want to
                                * center the coarse grid */
    int setcgcent;
    int fcmeth;                /* Fine grid centering method:  0 => center on
                                * point, 1 => center on molecule */
    double fcenter[3];         /* Fine grid center.  */
    int fcentmol;              /* Particular molecule on which we want to
                                * center the fine grid */
    int setfgcent;


    /* ********* TYPE 2 PARAMETERS (PARALLEL AUTO-FOCUS) ******** */
    double partDisjCenterShift[3];       /* When added to the actual (local)
                                          * mesh center, this gives the center
                                          * of the disjoint partitions */
    double partDisjLength[3];            /* This gives the lengths of the
                                          * disjoint partitions */
    int partDisjOwnSide[6];              /* Tells whether the boundary points
                                          * are ours (1) or not (0) */
    double partOlapCenterShift[3];       /* When added to the actual (local)
                                          * mesh center, this gives the center
                                          * of the overlapping partitions */
    double partOlapLength[3];            /* This gives the lengths of the
                                          * overlapping partitions */

    int pdime[3];                        /* Grid of processors to be used in
                                          * calculation */
    int setpdime;
    Vcom *com;                           /* Communications object for this 
                                          * processor */
    int setcom;
    double ofrac;                        /* Overlap fraction between procs */
    int setofrac; 

} MGparm;

/* ///////////////////////////////////////////////////////////////////////////
// Class NOsh: Non-inlineable methods (mcsh.c)
/////////////////////////////////////////////////////////////////////////// */

VEXTERNC MGparm*  MGparm_ctor(int type);
VEXTERNC int      MGparm_ctor2(MGparm *thee, int type);
VEXTERNC void     MGparm_dtor(MGparm **thee);
VEXTERNC void     MGparm_dtor2(MGparm *thee);
VEXTERNC int      MGparm_check(MGparm *thee);
VEXTERNC void     MGparm_copy(MGparm *thee, MGparm *parm);
VEXTERNC int      MGparm_parseToken(MGparm *thee, char tok[VMAX_BUFSIZE], 
                    Vio *sock);


#endif 

