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
/// PARTICULAR PURPOSE.  THE SOFTWARE AND PMGPOMPANYING DOCUMENTATION, IF
/// ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS".  REGENTS HAS NO OBLIGATION
/// TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
/// MODIFICATIONS. 
////////////////////////////////////////////////////////////////////////////
/// rcsid="$Id$"
//////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// File:     vpmgpp.h    < vpmgpp.c >
//
// Purpose:
//    Class Vpmgpp:
//      Parameter structure for Mike Holst's PMGP code
//     
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */

#ifndef _VPMGP_H_
#define _VPMGP_H_

#include "maloc/maloc.h"
#include "apbs/vhal.h"
#include "apbs/valist.h"
#include "apbs/vunit.h"

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpmgp: Parameters and datatypes
/////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpmgp: Definition
/////////////////////////////////////////////////////////////////////////// */

typedef struct Vpmgp {

    /* ********** USER-SPECIFIED PARAMETERS ********** */
    int nx, ny, nz;              /* Grid dimensions [no default]  */
    int nlev;                    /* Number of mesh levels [no default] */
    double hx, hy, hzed;         /* Grid spacings [no default]  */
    int nonlin;                  /* Problem type [no default]
                                  *   0 => linear
                                  *   1 => nonlinear
                                  *   2 => linear then nonlinear */

    /* ********** DERIVED PARAMETERS ********** */
    int nrwk;                    /* Real work storage */
    int niwk;                    /* Integer work storage */
    int narr;                    /* Array work storage */
    int ipkey;                   /* Toggles nonlinearity (set by nonlin)
                                  *  -1 => Linearized PBE
                                  *   0 => Nonlinear PBE with capped sinh 
                                  *        term [default]
                                  *  >1 => Polynomial approximation to sinh, 
                                  *        note that ipkey must be odd       */

    /* ********** PARAMETERS WITH DEFAULT VALUES ********** */
    double xcent, ycent, zcent;  /* Grid center [0, 0, 0]  */
    double errtol;               /* Desired error tolerance [default = 1e-9] */
    int itmax;                   /* Maximum number of iters [default = 100] */
    int istop;                   /* Stopping criterion [default = 1]
                                  *   0 => residual
                                  *   1 => relative residual
                                  *   2 => diff
                                  *   3 => errc
                                  *   4 => errd
                                  *   5 => aerrd */
    int iinfo;                   /* Runtime status messages [default = 1]
                                  *   0 => none
                                  *   1 => some
                                  *   2 => lots
                                  *   3 => more */
    int bcfl;                    /* Boundary condition method [default = 1]
                                  *   0 => zero boundary conditions
                                  *   1 => boundary condition approximated by
                                  *        single Debye-Huckel sphere for
                                  *        entire molecule
                                  *   2 => boundary condition approximated by 
                                  *        Debye-Huckel spheres for each atom 
                                  *   4 => boundary condition determined from 
                                  *        previously calculated solution 
                                  *        (i.e., for focusing calculations)
                                  */
    int key;                     /* Print solution to file [default = 0] 
                                  *   0 => no
                                  *   1 => yes */
    int iperf;                   /* Analysis of the operator [default = 0]
                                  *   0 => no
                                  *   1 => condition number
                                  *   2 => spectral radius
                                  *   3 => cond. number & spectral radius */
    int meth;                    /* Solution method [default = 2]
                                  *   0 => conjugate gradient multigrid
                                  *   1 => newton
                                  *   2 => multigrid
                                  *   3 => conjugate gradient
                                  *   4 => sucessive overrelaxation
                                  *   5 => red-black gauss-seidel
                                  *   6 => weighted jacobi
                                  *   7 => richardson */
    int mgkey;                   /* Multigrid method [default = 0]
                                  *   0 => variable v-cycle
                                  *   1 => nested iteration */
    int nu1;                     /* Number of pre-smoothings [default = 2] */
    int nu2;                     /* Number of post-smoothings [default = 2] */
    int mgsmoo;                  /* Smoothing method [default = 1]
                                  *   0 => weighted jacobi
                                  *   1 => gauss-seidel
                                  *   2 => SOR
                                  *   3 => richardson
                                  *   4 => cghs */
    int mgprol;                  /* Prolongation method [default = 0]
                                  *   0 => trilinear
                                  *   1 => operator-based
                                  *   2 => mod. operator-based */
    int mgcoar;                  /* Coarsening method [default = 2]
                                  *   0 => standard
                                  *   1 => harmonic
                                  *   2 => galerkin */
    int mgsolv;                  /* Coarse equation solve method [default = 1]
                                  *   0 => cghs
                                  *   1 => banded linpack */
    int mgdisc;                  /* Discretization method [default = 0]
                                  *   0 => finite volume
                                  *   1 => finite element */
    double omegal;               /* Linear relax parameter [default = 8e-1] */
    double omegan;               /* Nonlin relax parameter [default = 9e-1] */
    int irite;                   /* FORTRAN output unit [default = 8] */
    int ipcon;                   /* Preconditioning method [default = 3]
                                  *   0 => diagonal
                                  *   1 => ICCG 
                                  *   2 => ICCGDW
                                  *   3 => MICCGDW
                                  *   4 => none */
    double xlen, ylen, zlen;     /* Domain dimensions */
    double xmin, ymin, zmin;     /* Domain dimensions */
    double xmax, ymax, zmax;     /* Domain dimensions */
} Vpmgp;

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpmgp: Inlineable methods (vpmgp.c)
/////////////////////////////////////////////////////////////////////////// */

#if !defined(VINLINE_VPMGP)
#else /* if defined(VINLINE_VPMGP) */
#endif /* if !defined(VINLINE_VPMGP) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpmgp: Non-Inlineable methods (vpmgp.c)
/////////////////////////////////////////////////////////////////////////// */

VEXTERNC Vpmgp* Vpmgp_ctor(int nx, int ny, int nz, int nlev, 
  double hx, double hy, double hzed, int nonlin);
VEXTERNC int Vpmgp_ctor2(Vpmgp *thee, int nx, int ny, int nz, int nlev, 
  double hx, double hy, double hzed, int nonlin);
VEXTERNC void Vpmgp_dtor(Vpmgp **thee);
VEXTERNC void Vpmgp_dtor2(Vpmgp *thee);

#endif    /* ifndef _VPMGP_H_ */
