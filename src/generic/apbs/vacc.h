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

#include "maloc/maloc.h"
#include "apbs/vhal.h"
#include "apbs/valist.h"
#include "apbs/vunit.h"

/* ///////////////////////////////////////////////////////////////////////////
// Class Vacc: Parameters and datatypes
/////////////////////////////////////////////////////////////////////////// */
#define VACCMAXNBOR 20

/* ///////////////////////////////////////////////////////////////////////////
// Class Vacc: Definition
/////////////////////////////////////////////////////////////////////////// */

typedef struct Vacc {

  Vmem *vmem;                    /* Memory management object for this class */
  Valist *alist;                 /* Valist structure for list of atoms */
  int **atomIDs;                 /* An array of arrays of atom IDs in alist */
  int *natoms;                   /* An array telling how many pointers are 
                                  * stored in atoms[i] */
  int *atomFlags;                /* Flags for keeping track of atoms */
  double **sphere;               /* An array of points on the surface of a 
                                  * sphere */
  int nsphere;                   /* The number of points in thee->sphere */
  Vset acc;                      /* An integer array (to be treated as 
                                  * bitfields) of Vset type with length equal 
                                  * to the number of vertices in the mesh */
  double grid_lower_corner[3];   /* Hash table grid corner */
  double hx, hy, hzed;           /* Hash table grid spacings */
  int nx, ny, nz, n;             /* Hash table grid dimensions, 
                                  * n = nx*nz*ny */
  double max_radius;             /* Maximum probe radius */
  double *area;                  /* The contribution to the solvent-accessible
                                  * surface area from each atom.  This array is
                                  * filled by Vacc_totalSASA */


} Vacc;

/* ///////////////////////////////////////////////////////////////////////////
// Class Vacc: Inlineable methods (vacc.c)
/////////////////////////////////////////////////////////////////////////// */

#if !defined(VINLINE_VACC)
    VEXTERNC int Vacc_memChk(Vacc *thee);
#else /* if defined(VINLINE_VACC) */
#   define Vacc_memChk(thee) (Vmem_bytes((thee)->vmem))
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

VEXTERNC double Vacc_vdwAcc(Vacc *thee, double center[3]);
VEXTERNC double Vacc_ivdwAcc(Vacc *thee, double center[3], double radius);
VEXTERNC double Vacc_molAcc(Vacc *thee, double center[3], double radius);
VEXTERNC double Vacc_fastMolAcc(Vacc *thee, double center[3], double radius);
VEXTERNC double Vacc_splineAcc(Vacc *thee, double center[3], double window,
  double infrad);
VEXTERNC void Vacc_splineAccGrad(Vacc *thee, double center[3], double win,
  double infrad, int atomID, double *force);
VEXTERNC double** Vacc_sphere(Vacc *thee, int *npts);
VEXTERNC double Vacc_totalSASA(Vacc *thee, double radius);
VEXTERNC double Vacc_atomSASA(Vacc *thee, double radius, int iatom);
VEXTERNC void Vacc_getNborAtoms(Vacc *thee, double coord[3], int *atomIDs[6]);
/* ** FORCE EVALUATION STUFF ** */
VEXTERNC double** Vacc_qrule(Vacc *thee, int *npts);
VEXTERNC int Vacc_contact(Vacc *thee, double center[3], double radius, 
  int *ncontact, int *contacts);
VEXTERNC double** Vacc_srf(Vacc *thee, double probe_radius, int *npts);
VEXTERNC double** Vacc_SRsrf(Vacc *thee, double probe_radius, int *npts);

#endif    /* ifndef _VACC_H_ */
