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
// File:     vpmgp.c
//
// Purpose:  Class Vpmgp: methods.  This class wraps options which get passed
//                                  to the FORTRAN PMG code by Mike Holst.
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */

#include "apbscfg.h"
#include "apbs/vpmgp.h"

VEMBED(rcsid="$Id$")

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpmgp: Inlineable methods
/////////////////////////////////////////////////////////////////////////// */
#if !defined(VINLINE_VACC)
#endif /* if !defined(VINLINE_VACC) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpmgp: Non-inlineable methods
/////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmgp_ctor
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Vpmgp* Vpmgp_ctor(int nx, int ny, int nz, int nlev, double hx, 
  double hy, double hzed, int nonlin) {

    Vpmgp *thee = VNULL;

    /* Set up the structure */
    thee = Vmem_malloc(VNULL, 1, sizeof(Vpmgp) );
    VASSERT( thee != VNULL);
    VASSERT(Vpmgp_ctor2(thee, nx, ny, nz, nlev, hx, hy, hzed, nonlin));

    return thee;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmgp_ctor2
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vpmgp_ctor2(Vpmgp *thee, int nx, int ny, int nz, int nlev,
  double hx, double hy, double hzed, int nonlin) {

    /* Specified parameters */
    thee->nx = nx;
    thee->ny = ny;
    thee->nz = nz;
    thee->hx = hx;
    thee->hy = hy;
    thee->hzed = hzed;
    thee->xlen = ((double)(nx-1))*hx;
    thee->ylen = ((double)(ny-1))*hy;
    thee->zlen = ((double)(nz-1))*hzed;
    thee->nlev = nlev; 
    thee->nonlin = nonlin;
    if (nonlin == 0) thee->ipkey = -1;
    else thee->ipkey = 0;


    /* Default parameters */
    thee->errtol = 1.0e-6;   /* Here are a few comments.  Mike had this set to
                              * 1e-9; convential wisdom sets this at 1e-6 for
                              * the PBE; Ray Luo sets this at 1e-3 for his
                              * accelerated PBE (for dynamics, etc.) */
    thee->itmax = 200;
    thee->istop = 1;
    thee->iinfo = 1;         /* I'd recommend either 1 (for debugging LPBE) or 
                              * 2 (for debugging NPBE), higher values give 
                              * too much output */
    thee->bcfl = 1;
    thee->key = 0;
    thee->iperf = 0;
    if (thee->nonlin == 1) { 
        Vnm_print(0, "Vpmp_ctor2:  Using meth = 1, mgcoar = 2, mgsolv = 0\n");
        thee->meth = 1;
        thee->mgcoar = 2;
        thee->mgsolv = 0;
    } else {                 
        Vnm_print(0, "Vpmp_ctor2:  Using meth = 0, mgcoar = 2, mgsolv = 0\n");
#if 0                               /* Fastest convergence, but lots of mem */
        thee->meth = 0;
        thee->mgcoar = 2;
        thee->mgsolv = 0;
#else                               /* Most rigorous (good for testing) */
        thee->meth = 2;
        thee->mgcoar = 2;
        thee->mgsolv = 1;
#endif
    }
    thee->mgkey = 0;
    thee->nu1 = 2;
    thee->nu2 = 2;
    thee->mgsmoo = 1;
    thee->mgprol = 0;
    thee->mgdisc = 0;
    thee->omegal = 8.0e-1;
    thee->omegan = 9.0e-1;
    thee->ipcon = 3;
    thee->irite = 8;
    thee->xcent = 0.0;
    thee->ycent = 0.0;
    thee->zcent = 0.0;

    return 1;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmgp_dtor
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpmgp_dtor(Vpmgp **thee) {
    
    if ((*thee) != VNULL) {
        Vpmgp_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(Vpmgp), (void **)thee);
        (*thee) = VNULL;
    }

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmgp_dtor2
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpmgp_dtor2(Vpmgp *thee) { ; }
