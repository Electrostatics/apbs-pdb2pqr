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
/// PARTICULAR PURPOSE.  THE SOFTWARE AND PMGOMPANYING DOCUMENTATION, IF
/// ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS".  REGENTS HAS NO OBLIGATION
/// TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
/// MODIFICATIONS. 
////////////////////////////////////////////////////////////////////////////
/// rcsid="$Id$"
//////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// File:     vpmg.h    < vpmg.c >
//
// Purpose:
//    Class Vpmg:
//      A wrapper for Mike Holst's PMG multigrid code.  Many of the routines
//      and macros are borrowed from the main.c driver (written by Mike Holst)
//      provided with the PMG code.
//     
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */

#ifndef _VPMG_H_
#define _VPMG_H_

#include "mc/mc.h"
#include "apbs/vhal.h"
#include "apbs/vpmgp.h"
#include "apbs/vunit.h"
#include "apbs/vpbe.h"
#include "apbs/valist.h"
#include "apbs/vacc.h"

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpmg: Parameters and datatypes
/////////////////////////////////////////////////////////////////////////// */
#define VPMGMAXPART 2000  /* The max number of partitions we can divide the
                           * mesh into */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpmg: Definition
/////////////////////////////////////////////////////////////////////////// */

typedef struct Vpmg {

  Vmem *vmem;                    /* Memory management object for this class */
  Vpmgp *pmgp;                   /* Parameters */
  Vpbe *pbe;                     /* Information about the PBE system */

  int *iparm;                    /* Passing int parameters to FORTRAN */
  double *rparm;                 /* Passing real parameters to FORTRAN */
  int *iwork;                    /* Work array */
  double *rwork;                 /* Work array */
  double *a1cf, *a2cf, *a3cf;    /* Operator coefficient values */
  double *ccf;                   /* Helmholtz term */
  double *fcf;                   /* Right-hand side */
  double *tcf;                   /* True solution */
  double *u;                     /* Solution */
  double *xf, *yf, *zf;          /* Mesh point coordinates */
  double *gxcf, *gycf, *gzcf;    /* Boundary conditions */
  int partFlag;                  /* Flag which restricts calculation of
                                  * energies to the specified partition */
  double partLower[3];           /* Lower corner of partition */
  double partUpper[3];           /* Upper corner of partition */
  double extEnergy;              /* Storing the contribution to the energy from
                                  * regions outside the immediate problem
                                  * domain */

} Vpmg;

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpmg: Inlineable methods (vpmg.c)
/////////////////////////////////////////////////////////////////////////// */

#if !defined(VINLINE_VPMG)
    VEXTERNC int Vpmg_memChk(Vpmg *thee);
#else /* if defined(VINLINE_VPMG) */
#   define Vpmg_memChk(thee) (Vmem_bytes((thee)->vmem))
#endif /* if !defined(VINLINE_VPMG) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpmg: Non-Inlineable methods (vpmg.c)
/////////////////////////////////////////////////////////////////////////// */

VEXTERNC Vpmg* Vpmg_ctor(Vpmgp *parms, Vpbe *pbe);
VEXTERNC int Vpmg_ctor2(Vpmg *thee, Vpmgp *parms, Vpbe *pbe);
VEXTERNC Vpmg* Vpmg_ctorFocus(Vpmgp *parms, Vpbe *pbe, Vpmg *pmgOLD);
VEXTERNC int Vpmg_ctor2Focus(Vpmg *thee, Vpmgp *parms, Vpbe *pbe, Vpmg *pmgOLD);
VEXTERNC void Vpmg_dtor(Vpmg **thee);
VEXTERNC void Vpmg_dtor2(Vpmg *thee);

VEXTERNC void Vpmg_fillco(Vpmg *thee, int epsmeth, double epsparm);
VEXTERNC void Vpmg_solve(Vpmg *thee);
VEXTERNC double Vpmg_getLinearEnergy1(Vpmg *thee, int extFlag);
VEXTERNC double Vpmg_getLinearEnergy2(Vpmg *thee, int extFlag);
VEXTERNC double Vpmg_getAnorm(Vpmg *thee, double *u);
VEXTERNC void Vpmg_getAdet(Vpmg *thee, double det[2]);
VEXTERNC void Vpmg_writeUHBD(Vpmg *thee, const char *iodev, const char *iofmt,
  const char *thost, const char *fname, char *title, double *data);
VEXTERNC void Vpmg_writeDX(Vpmg *thee, const char *iodev, const char *iofmt,
  const char *thost, const char *fname, char *title, double *data);
VPUBLIC int Vpmg_readArrayDX(const char *iodev, const char *iofmt,
  const char *thost, const char *fname, double *xlow, double *ylow, 
  double *zlow, double *hx, double *hy, double *hzed, int *nx, int *ny, 
  int *nz, double **data);
VEXTERNC void Vpmg_setPart(Vpmg *thee, double xmin, double ymin, double zmin,
           double xmax, double ymax, double zmax);
VEXTERNC void Vpmg_fillAcc(Vpmg *thee, double *vec, int meth, double parm);

#endif    /* ifndef _VPMG_H_ */
