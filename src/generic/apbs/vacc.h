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
// File:     vacc.h    < vaacc.c >
//
// Purpose:
//    Class Vacc:
//      Determines (protein) accessibility of point with respect probe.
//      Uses list of atoms.
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */

#ifndef _VACC_H_
#define _VACC_H_

#include "mc/vhal.h"
#include "mc/vec3.h"
#include "mc/vset.h"
#include "mc/vmem.h"
#include "mc/vgm.h"

#include "apbs/valist.h"
#include "apbs/vunit.h"

/* ///////////////////////////////////////////////////////////////////////////
// Class Vacc: Parameters and datatypes
/////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vacc: Definition
/////////////////////////////////////////////////////////////////////////// */

typedef struct Vacc {

  Vmem *vmem;                    /* Memory management object for this class */
  Vatom ***atoms;                /* An array of arrays of pointers to atoms */
  int *natoms;                   /* An array telling how many pointers are 
                                  * stored in atoms[i] */
  double **sphere;               /* An array of points on the surface of a 
                                  * sphere */
  int nsphere;
  Vset acc;                      /* An integer array (to be treated as 
                                  * bitfields) of Vset type with length equal 
                                  * to the number of vertices in the mesh */
  Vec3 grid_lower_corner;        /* Hash table grid corner */
  double hx, hy, hzed;           /* Hash table grid spacings */
  int nx, ny, nz, n;             /* Hash table grid dimensions, 
                                  * n = nx*nz*ny */
  double max_radius;             /* Maximum probe radius */

} Vacc;

/* ///////////////////////////////////////////////////////////////////////////
// Class Vacc: Inlineable methods (vacc.c)
/////////////////////////////////////////////////////////////////////////// */

#if !defined(VINLINE_VACC)
#else /* if defined(VINLINE_VACC) */
#endif /* if !defined(VINLINE_VACC) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vacc: Non-Inlineable methods (vacc.c)
/////////////////////////////////////////////////////////////////////////// */

VEXTERNC Vacc* Vacc_ctor(Valist *alist, double max_radius, int nx, 
               int ny, int nz, int nsphere);
VEXTERNC int Vacc_ctor2(Vacc *thee, Valist *alist, double max_radius, 
             int nx, int ny, int nz, int nsphere);
VEXTERNC void Vacc_dtor(Vacc **thee);
VEXTERNC void Vacc_dtor2(Vacc *thee);

VEXTERNC int Vacc_memChk(Vacc *thee);

VEXTERNC int Vacc_vdwAcc(Vacc *thee, Vec3 center);
VEXTERNC int Vacc_ivdwAcc(Vacc *thee, Vec3 center, double radius);
VEXTERNC int Vacc_molAcc(Vacc *thee, Vec3 center, double radius);
VEXTERNC double** Vacc_sphere(Vacc *thee, int *npts);
VEXTERNC void Vacc_writeGMV(Vacc *thee, double radius, int meth, Vgm *gm,
  char *iodev, char *iofmt, char *iohost, char *iofile);


#endif    /* ifndef _VACC_H_ */
