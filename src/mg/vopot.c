/**
 *  @file    vopot.c
 *  @author  Nathan Baker
 *  @brief   Class Vopot methods
 *  @ingroup Vopot
 *  @version $Id$
 *  @attention
 *  @verbatim
 *
 * APBS -- Adaptive Poisson-Boltzmann Solver
 *
 * Nathan A. Baker (nbaker@wasabi.ucsd.edu)
 * Dept. of Chemistry and Biochemistry
 * University of California, San Diego 
 *
 * Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 1999-2002.  The Regents of the University of California.
 * Portions Copyright (c) 1995.  Michael Holst.
 *
 * Permission to use, copy, modify, and distribute this software and its
 * documentation for educational, research, and not-for-profit purposes,
 * without fee and without a signed licensing agreement, is hereby granted,
 * provided that the above copyright notice, this paragraph and the
 * following two paragraphs appear in all copies, modifications, and
 * distributions.
 *
 * IN NO EVENT SHALL THE AUTHORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE
 * AUTHORS HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * THE AUTHORS SPECIFICALLY DISCLAIM ANY WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE.  THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED
 * HEREUNDER IS PROVIDED "AS IS".  THE AUTHORS HAVE NO OBLIGATION TO PROVIDE
 * MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
 * @endverbatim
 */

#include "apbscfg.h"
#include "apbs/vopot.h"

VEMBED(rcsid="$Id$")

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vopot_ctor
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Vopot* Vopot_ctor(Vgrid *grid, Vpbe *pbe, int bcfl) {

    Vopot *thee = VNULL;

    thee = Vmem_malloc(VNULL, 1, sizeof(Vopot));
    VASSERT(thee != VNULL);
    VASSERT(Vopot_ctor2(thee, grid, pbe, bcfl));

    return thee;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vopot_ctor2
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vopot_ctor2(Vopot *thee, Vgrid *grid, Vpbe *pbe, int bcfl) {

    if (thee == VNULL) return 0;
    if ((bcfl < 0) || (bcfl > 2)) {
        Vnm_print(2, "Vopot_ctor:  Bogus bcfl flag (%d)!\n", bcfl);
        return 0;
    }
    thee->bcfl = bcfl;
    thee->grid = grid;
    thee->pbe = pbe;

    return 1;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vopot_dtor
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vopot_dtor(Vopot **thee) {

    if ((*thee) != VNULL) {
        Vopot_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(Vopot), (void **)thee);
        (*thee) = VNULL;
    }
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vopot_dtor2
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vopot_dtor2(Vopot *thee) { return; }

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vopot_pot
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
#define IJK(i,j,k)  (((k)*(nx)*(ny))+((j)*(nx))+(i))
VPUBLIC int Vopot_pot(Vopot *thee, double pt[3], double *value) {

    Vatom *atom;
    int i, iatom;
    double u, T, charge, eps_w, xkappa, dist, size, val, *position;
    Valist *alist;

    VASSERT(thee != VNULL);

    eps_w = Vpbe_getSolventDiel(thee->pbe);
    xkappa = (1.0e10)*Vpbe_getXkappa(thee->pbe);
    T = Vpbe_getTemperature(thee->pbe);
    alist = Vpbe_getValist(thee->pbe);

    u = 0;

    /* See if we're on the mesh */
    if (Vgrid_value(thee->grid, pt, &u)) {

        *value = u;

    } else {

        switch (thee->bcfl) {

            case 0:
                u = 0;
                break;

            case 1:
                size = (1.0e-10)*Vpbe_getSoluteRadius(thee->pbe);
                position = Vpbe_getSoluteCenter(thee->pbe);
                charge = Vunit_ec*Vpbe_getSoluteCharge(thee->pbe);
                dist = 0;
                for (i=0; i<3; i++)
                  dist += VSQR(position[i] - pt[i]);
                dist = (1.0e-10)*VSQRT(dist);
                val = (charge)/(4*VPI*Vunit_eps0*eps_w*dist);
                if (xkappa != 0.0) 
                  val = val*(exp(-xkappa*(dist-size))/(1+xkappa*size));
                val = val*Vunit_ec/(Vunit_kb*T);
                u = val;
                break;

            case 2:
                u = 0;
                for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {
                    atom = Valist_getAtom(alist, iatom);
                    position = Vatom_getPosition(atom);
                    charge = Vunit_ec*Vatom_getCharge(atom);
                    size = (1e-10)*Vatom_getRadius(atom);
                    dist = 0;
                    for (i=0; i<3; i++)
                      dist += VSQR(position[i] - pt[i]);
                    dist = (1.0e-10)*VSQRT(dist);
                    val = (charge)/(4*VPI*Vunit_eps0*eps_w*dist);
                    if (xkappa != 0.0)
                      val = val*(exp(-xkappa*(dist-size))/(1+xkappa*size));
                    val = val*Vunit_ec/(Vunit_kb*T);
                    u = u + val;
                }
                break;

            default:
                Vnm_print(2, "Vopot_pot:  Bogus thee->bcfl flag (%d)!\n", 
                  thee->bcfl);
                return 0;        
                break;
        }

        *value = u;

    }

    return 1;

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vopot_curvature
//
//   Notes:  cflag=0 ==> Reduced Maximal Curvature
//           cflag=1 ==> Mean Curvature (Laplace)
//           cflag=2 ==> Gauss Curvature
//           cflag=3 ==> True Maximal Curvature
//
// Authors:  Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vopot_curvature(Vopot *thee, double pt[3], int cflag, 
  double *value) {

    if (!Vgrid_curvature(thee->grid, pt, cflag, value)) {
        Vnm_print(2, "Vopot_curvature:  Off mesh!\n");
        return 0;
    } 

    return 1;

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vopot_gradient
//
// Authors:  Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vopot_gradient(Vopot *thee, double pt[3], double grad[3]) {

    if (!Vgrid_gradient(thee->grid, pt, grad)) {
        Vnm_print(2, "Vopot_curvature:  Off mesh!\n");
        return 0;
    } 

    return 1;


}

