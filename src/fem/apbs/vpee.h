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
// File:     vpee.h    < vpee.c >
//
// Purpose: Parallel error estimator class
//
//    Class Vpee:  This class provides some functionality for error esimation
//      in parallel.  The purpose is to modulate the error returned by some
//      external error estimator according to the partitioning of the mesh.
//      For example, the Bank/Holst parallel refinement routine essentially
//      reduces the error outside the ``local" partition to zero.  However, 
//      this leads to the need for a few final overlapping Schwarz solves to
//      smooth out the errors near partition boundaries.  Supposedly, if the
//      region in which we allow error-based refinement includes the ``local"
//      partition and an external buffer zone approximately equal in size to
//      the local region, then the solution will asymptotically approach the
//      solution obtained via more typical methods.
//
//      This is essentially a more flexible parallel implementation of MC's
//      AM_markRefine.
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */

#ifndef _VPEE_H
#define _VPEE_H

#include "mc/vgm.h"
#include "mc/am.h"
#include "mc/alg.h"
#include "mc/ves.h"
#include "mc/vmem.h"

#include "apbs/vhal.h"

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpee: Parameters and datatypes
/////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpee: Definition
/////////////////////////////////////////////////////////////////////////// */

typedef struct Vpee {

  Vgm *gm;                     /* Grid manager */
  int localPartID;             /* The local partition ID: i.e. the partition 
                                * whose boundary simplices we're keeping
                                * track of */
  double localPartCenter[3];   /* The coordinates of the center of the local
                                * partition */
  double localPartRadius;      /* The radius of the circle/sphere which
                                * circumscribes the local partition */
  int killFlag;                /* A flag indicating the method we're using to
                                * artificially decrease the error esimate
                                * outside the local partition */
  double killParam;            /* A parameter for the error estimate
                                * attenuation method */
  Vmem *mem;                   /* Memory manager */

} Vpee;

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpee Inlineable methods 
/////////////////////////////////////////////////////////////////////////// */

#if !defined(VINLINE_VPEE)
#else /* if defined(VINLINE_VPEE) */
#endif /* if !defined(VINLINE_VPEE) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpee: Non-Inlineable methods (vpee.c)
/////////////////////////////////////////////////////////////////////////// */

VEXTERNC Vpee* Vpee_ctor(Vgm *gm, int localPartID, int killFlag, 
                 double killParam);
VEXTERNC int   Vpee_ctor2(Vpee *thee, Vgm *gm, int localPartID, int killFlag,
                 double killParam);
VEXTERNC void  Vpee_dtor(Vpee **thee);
VEXTERNC void  Vpee_dtor2(Vpee *thee);

VEXTERNC int   Vpee_numSS(Vpee *thee);
VEXTERNC void  Vpee_estimate(Vpee *thee, AM *am, int level, int akey); 
VEXTERNC int   Vpee_markRefine(Vpee *thee, AM *am, int level, int akey, 
                 int rcol, double etol);
VEXTERNC int   Vpee_memChk(Vpee *thee);

#endif    /* ifndef _VPEE_H_ */
