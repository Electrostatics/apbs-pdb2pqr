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
// File:     vgreen.c
//
// Purpose:  Class Vgreen: methods. 
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */

#include "apbscfg.h"
#include "apbs/vgreen.h"
#include "cxxfmm.h"

#define DTORCXXFMM      Vgreen_dtorCXXFMM__FP6Vgreen
#define FIELDCXXFMM     Vgreen_fieldCXXFMM__FP6VgreenPdT1
#define INITCXXFMM      Vgreen_initCXXFMM__FP6Vgreendiiiddd
#define POTCXXFMM       Vgreen_potentialCXXFMM__FP6VgreenPd
#define UPDATECXXFMM    Vgreen_updateCXXFMM__FP6Vgreen

/* ///////////////////////////////////////////////////////////////////////////
// Class Vgreen: Inlineable methods
/////////////////////////////////////////////////////////////////////////// */
#if !defined(VINLINE_VGREEN)

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vgreen_getValist
//
// Purpose:  Get a pointer to the Valist (atom list) object
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Valist* Vgreen_getValist(Vgreen *thee) { 

   VASSERT(thee != VNULL);
   return thee->alist;

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vgreen_memChk
//
// Purpose:  Return number of bytes used by this object
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vgreen_memChk(Vgreen *thee) {
    if (thee == VNULL) return 0;
    return Vmem_bytes(thee->vmem);
}

#endif /* if !defined(VINLINE_VCSM) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vgreen: Non-inlineable methods
/////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vgreen_ctor
//
// Purpose:  Construct the Green's function oracle
//
// Notes:    This doesn't have any multipole method implementation yet, so it's
//           most likeley slow
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Vgreen* Vgreen_ctor(Valist *alist) {

    /* Set up the structure */
    Vgreen *thee = VNULL;
    thee = (Vgreen *)Vmem_malloc(VNULL, 1, sizeof(Vgreen) );
    VASSERT( thee != VNULL);
    VASSERT( Vgreen_ctor2(thee, alist));

    return thee;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vgreen_ctor2
//
// Purpose:  Construct the Vgreen object
//
// Notes:    Constructor broken into two parts for FORTRAN users.
//
// Returns:  1 if sucessful, 0 otherwise
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vgreen_ctor2(Vgreen *thee, Valist *alist) { 
 
    VASSERT( thee != VNULL );

    /* Memory management object */
    thee->vmem = Vmem_ctor("APBS:VGREEN");

    /* Set up the atom list and grid manager */
    if (alist == VNULL) {
        Vnm_print(2,"Vgreen_ctor2: got null pointer to Valist object!\n");
    }

    thee->alist = alist;
    thee->initFlagCXXFMM = 0;
   
    return 1;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vgreen_initFMM
//
// Purpose:  Construct the Vgreen FMM object
//
// Args:     spacing   Spacing for the FMM cells
//           n[xyz]    Number of FMM cells in the specified direction
//           [xyz]low  Lower corner of the cells
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vgreen_initFMM(Vgreen *thee, double spacing, int nx, int ny,
  int nz, double xlow, double ylow, double zlow) {

    VASSERT( thee != VNULL );

#if defined(USE_CXX_FMM)
    INITCXXFMM(thee, spacing, nx, ny, nz, xlow, ylow, zlow);
    thee->initFlagCXXFMM = 1;
#else
    Vnm_print(2, "Vgreen_initFMM: Not compiled with FMM support!\n");
    VASSERT(0);
#endif 

}


/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vgreen_dtor
//
// Purpose:  Destroy the charge-simplex map.
// 
// Notes:    Since the grid manager and atom list were allocated outside of
//           the Vgreen routines, they are not destroyed.
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vgreen_dtor(Vgreen **thee) {
    if ((*thee) != VNULL) {
        Vgreen_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(Vgreen), (void **)thee);
        (*thee) = VNULL;
    }
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vgreen_dtor2
//
// Purpose:  Destroy the atom object
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vgreen_dtor2(Vgreen *thee) { 

    Vmem_dtor(&(thee->vmem));

#if defined(USE_CXX_FMM)
    if (thee->initFlagCXXFMM) DTORCXXFMM(thee);
#endif

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vgreen_helmholtz
//
// Purpose:  Return the Green's function for Helmholtz's equation
//           integrated over the atomic point charges
//
//           G(r) = \sum_i \frac{q_i e^{-\kappa r_i}}{r_i}
//
//           where \kappa is the inverse screening length (in Angstroms),
//           q_i is the atomic charge (in e), and r_i is the distance from atom
//           i to the observation point.  The potential returned is in units of
//           J/C.
//
// Note:     This quantity needs to be multiplied by (1/\epsilon)
//           give the actual potential.
//
// Args:     position = vector containing the position of the observation pt
//           dim      = number of elements in position
//           kappa    = Helmholtz coefficient (units of inverse length)
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vgreen_helmholtz(Vgreen *thee, double *position, double dim, 
  double kappa) { 

    VASSERT(0);
    return 0.;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vgreen_helmholtzD
//
// Purpose:  Return the gradient of the Green's function for Helmholtz's
//           equation integrated over the atomic point charges
//
//           G(r) = \nabla (\sum_i \frac{q_i e^{-\kappa r_i}}{r_i})
//
//           where \kappa is the inverse screening length (in Angstroms),
//           q_i is the atomic charge (in e), and r_i is the distance from atom
//           i to the observation point.  The field is returned in units of
//           J/C/A.
//
// Note:     This quantity needs to be multiplied by (1/\epsilon) to give the
//           actual field.
//
// Args:     position = vector containing the position of the observation pt
//           dim      = number of elements in position
//           kappa    = Helmholtz coefficient (units of inverse length)
//           grad     = storage for the gradient of the Green's function (dim
//                      elements)
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vgreen_helmholtzD(Vgreen *thee, double *position, 
  double dim, double kappa, double *grad) {
    VASSERT(0);
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vgreen_coulomb
//
// Purpose:  Return the Coulomb's Law Green's function for Poisson's equation
//           integrated over the atomic point charges
//
//           G(r) = \sum_i \frac{q_i}{r_i}
//
//           where q_i is the atomic charge (in e), and r_i is the distance
//           from atom i to the observation point.  The potential is returned
//           in units of J/C.
//
// Note:     This quantity needs to be multiplied by (1/\epsilon) to give the
//           potential.
//
// Args:     position = vector containing the position of the observation pt
//           dim      = number of elements in position
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vgreen_coulomb(Vgreen *thee, double *position, double dim) {

    Vatom *atom;
    double *x, pot, charge, dist;
    int iatom, j;

    VASSERT(dim < 4);
    pot = 0;
  
#if defined(USE_CXX_FMM)
    VASSERT(thee->initFlagCXXFMM);
    Vnm_print(2, "Vgreen_coulomb: pos = (%g, %g, %g)\n", position[0],
position[1], position[2]);
    pot = POTCXXFMM(thee, position);
#else
    for (iatom=0; iatom<Valist_getNumberAtoms(thee->alist); iatom++) {
        atom = Valist_getAtom(thee->alist, iatom);
        x = Vatom_getPosition(atom);
        charge = Vatom_getCharge(atom);
        dist = 0;
        for (j=0; j<dim; j++) 
          dist += ((x[j] - position[j])*(x[j] - position[j]));
        dist = 1.0e-10*VSQRT(dist);
        pot += (charge/dist);
    }
#endif

    return pot*Vunit_ec/(4*Vunit_pi*Vunit_eps0);
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vgreen_coulombD
//
// Purpose:  Return the gradient for the Coulomb's Law Green's function for
//           Poisson's equation integrated over the atomic point charges
//
//           G(r) = \sum_i \frac{q_i}{r_i}
//
//           where q_i is the atomic charge (in e), and r_i is the distance
//           from atom i to the observation point.  The field is returned in
//           units of J/C/A.
//
// Note:     This quantity needs to be scaled by (1/\epsilon) to give the
//           actual field.
//
// Args:     position = vector containing the position of the observation pt
//           dim      = number of elements in position
//           grad     = storage for the gradient (at least dim elems)
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vgreen_coulombD(Vgreen *thee, double *position, double dim,
  double *grad) {

    Vatom *atom;
    double *x, tpot[MAXV], charge, dist;
    int iatom, j;

    VASSERT(dim < 4);
    
    for (j=0; j<dim; j++) grad[j] = 0;
  
#if defined(USE_CXX_FMM)
    VASSERT(thee->initFlagCXXFMM);
    FIELDCXXFMM(thee, position, grad);
#else
    for (iatom=0; iatom<Valist_getNumberAtoms(thee->alist); iatom++) {
        atom = Valist_getAtom(thee->alist, iatom);
        x = Vatom_getPosition(atom);
        charge = Vatom_getCharge(atom);
        dist = 0;
        for (j=0; j<dim; j++) 
          dist += ((x[j] - position[j])*(x[j] - position[j]));
        for (j=0; j<dim; j++) 
          grad[j] -= (charge*(position[j] - x[j])/dist/VSQRT(dist));
    }
#endif

    for (j=0; j<dim; j++) 
      grad[j] = grad[j]*Vunit_ec/(4*Vunit_pi*Vunit_eps0*(1.0e-10));
}

