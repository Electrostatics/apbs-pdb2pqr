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
// File:     vhal.h   
//
// Purpose:  
//      Provides macros for general configuration of APBS.
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */

#ifndef _VAPBSHAL_H_
#define _VAPBSHAL_H_

/* ///////////////////////////////////////////////////////////////////////////
// Macros for space allocation.
/////////////////////////////////////////////////////////////////////////// */
#define MAXMOL 5             /* The maximum number of molecules that can be
                                involved in a single PBE calculation */
#define MAXION 10            /* The maximum number of ion species */
#define MAXFOCUS 5           /* The maximum number of times an MG calculation 
                                can be focused */
#define VMGNLEV 4            /* Desired number of levels in multigrid
                              * hierarchies */
#define VREDFRAC 0.25        /* Max reduction of grid spacing in focusing
                              * calculations */
#define VAPBS_RIGHT 0        /* Directions for looking at parallel focusing */
#define VAPBS_FRONT 1        /* partitions; consistent with PMG if RIGHT = */
#define VAPBS_UP    2        /* EAST, BACK = SOUTH */
#define VAPBS_LEFT  3
#define VAPBS_BACK  4
#define VAPBS_DOWN  5

/* ///////////////////////////////////////////////////////////////////////////
// Inlining via macros for speed.
// 
// If you want to debug, do not define these
/////////////////////////////////////////////////////////////////////////// */
#if !defined(VDEBUG)
#   define VINLINE_VACC
#   define VINLINE_VATOM
#   define VINLINE_VCSM
#   define VINLINE_VPBE
#   define VINLINE_VPEE
#   define VINLINE_VGREEN
#   define VINLINE_VFETK
#   define VINLINE_VPMG
#endif


#endif
