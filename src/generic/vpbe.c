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
// File:     vpbe.c
//
// Purpose:  Class Vpbe: methods. 
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */

#include "apbscfg.h"
#include "apbs/vpbe.h"

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpbe: Private method declaration
/////////////////////////////////////////////////////////////////////////// */
#define MAX_SPLINE_WINDOW 0.5
#define MAX_HASH_DIM 75
#define VACC_SPHERE 200

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpbe: Inlineable methods
/////////////////////////////////////////////////////////////////////////// */
#if !defined(VINLINE_VPBE)

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpbe_getValist
//
// Purpose:  Get a pointer to the Valist (atom list) object
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Valist* Vpbe_getValist(Vpbe *thee) { 

   VASSERT(thee != VNULL);
   return thee->alist;

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpbe_getVgreen
//
// Purpose:  Get a pointer to the Vgreen (grid manager) object
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Gem* Vpbe_getVgreen(Vpbe *thee) {

   VASSERT(thee != VNULL);
   return thee->green;

}


/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpbe_getVacc
//
// Purpose:  Get a pointer to the Vacc accessibility object 
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Vacc* Vpbe_getVacc(Vpbe *thee) { 

   VASSERT(thee != VNULL);
   VASSERT(thee->paramFlag);
   return thee->acc; 

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpbe_getIonConc
//
// Purpose:  Get the ionic strength in M
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vpbe_getIonConc(Vpbe *thee) { 

   VASSERT(thee != VNULL);
   VASSERT(thee->paramFlag);
   return thee->ionConc; 

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpbe_getTemperature
//
// Purpose:  Get the temperature in K
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vpbe_getTemperature(Vpbe *thee) { 

   VASSERT(thee != VNULL);
   VASSERT(thee->paramFlag);
   return thee->T; 

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpbe_getSoluteDiel
//
// Purpose:  Get the solute dielectric
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vpbe_getSoluteDiel(Vpbe *thee) { 

   VASSERT(thee != VNULL);
   VASSERT(thee->paramFlag);
   return thee->soluteDiel; 

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpbe_getSoluteCenter
//
// Purpose:  Get the center of the solute molecule
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double* Vpbe_getSoluteCenter(Vpbe *thee) { 

   VASSERT(thee != VNULL);
   return thee->soluteCenter; 
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpbe_getSolventDiel
//
// Purpose:  Get the solvent dielectric
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vpbe_getSolventDiel(Vpbe *thee) { 

   VASSERT(thee != VNULL);
   VASSERT(thee->paramFlag);
   return thee->solventDiel; 
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpbe_getSolventRadius
//
// Purpose:  Get the solvent radius in angstroms
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vpbe_getSolventRadius(Vpbe *thee) { 

   VASSERT(thee != VNULL);
   VASSERT(thee->paramFlag);
   return thee->solventRadius; 
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpbe_getIonRadius
//
// Purpose:  Get the ion probe radius in angstroms 
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vpbe_getIonRadius(Vpbe *thee) { 

   VASSERT(thee != VNULL);
   VASSERT(thee->paramFlag);
   return thee->ionRadius; 
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpbe_getXkappa
//
// Purpose:  Get the Debye-Huckel parameter in reciprocal angstroms
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vpbe_getXkappa(Vpbe *thee) { 

   VASSERT(thee != VNULL);
   VASSERT(thee->paramFlag);
   return thee->xkappa; 
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpbe_getDeblen
//
// Purpose:  Get the Debye length in angstroms
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vpbe_getDeblen(Vpbe *thee) { 

   VASSERT(thee != VNULL);
   VASSERT(thee->paramFlag);
   return thee->deblen; 
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpbe_getZkappa2
//
// Purpose:  Get the squared modified Debye-Huckel parameter
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vpbe_getZkappa2(Vpbe *thee) { 

   VASSERT(thee != VNULL);
   VASSERT(thee->paramFlag);
   return thee->zkappa2; 
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpbe_getZmagic
//
// Purpose:  Get the delta function scaling factor
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vpbe_getZmagic(Vpbe *thee) { 

   VASSERT(thee != VNULL);
   VASSERT(thee->paramFlag);
   return thee->zmagic; 
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpbe_getSoluteRadius
//
// Purpose:  Get the radius of the solute molecule
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vpbe_getSoluteRadius(Vpbe *thee) { 

   VASSERT(thee != VNULL);
   return thee->soluteRadius; 
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpbe_getSoluteXlen
//
// Purpose:  Solute length in x-direction
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vpbe_getSoluteXlen(Vpbe *thee) { 

   VASSERT(thee != VNULL);
   return thee->soluteXlen; 
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpbe_getSoluteYlen
//
// Purpose:  Solute length in Y-direction
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vpbe_getSoluteYlen(Vpbe *thee) { 

   VASSERT(thee != VNULL);
   return thee->soluteYlen; 
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpbe_getSoluteZlen
//
// Purpose:  Solute length in Z-direction
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vpbe_getSoluteZlen(Vpbe *thee) { 

   VASSERT(thee != VNULL);
   return thee->soluteZlen; 
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpbe_getSoluteCharge
//
// Purpose:  Get the charge of the solute molecule
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vpbe_getSoluteCharge(Vpbe *thee) { 

   VASSERT(thee != VNULL);
   return thee->soluteCharge; 
}
#endif /* if !defined(VINLINE_VPBE) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpbe: Non-inlineable methods
/////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpbe_ctor
//
// Purpose:  Set up parameters, Vacc objects, and charge-simplex map
//
// Arguments: ionConc       = ionic strength in M
//            ionRadius     = ionic probe radius in A
//            T             = temperature in K
//            soluteDiel    = solute dielectric (unitless)
//            solventDiel   = solvent dielectric (unitless)
//            solventRadius = solvent radius in Angstroms
//
// Notes:  Here's the original function coments from Mike
//
//    kappa is defined as follows:
//
//       kappa^2 = (8 pi N_A e_c^2) I_s / (1000 eps_w k_B T)
//
//    note that the units here are:  esu^2 / erg-mole.
//    to obtain angstroms^{-2}, we multiply by 10^{-16}.
//
//    thus, in angstroms^{-2}, where k_B and e_c are in gaussian
//    rather than mks units, the proper value for kappa is:
//
//       kappa^2 = 8 pi N_A e_c^2 I_s / (1000 eps_w k_b T)  * 1.0e-16
//
//    where:
//
//       e_c = 4.803242384e-10 esu  (rather than 1.6021892e-19 coulombs)
//       k_B = 1.380662000e-16 ers  (rather than 1.380662000e-23 joules)
//
//    and the factor of 1.0e-16 results from converting cm^2 to anstroms^2,
//    noting that the 1000 in the denominator has converted m^3 to cm^3,
//    since the ionic strength I_s is assumed to have been provided in
//    moles
//    per liter, which is moles per 1000 cm^3.
//
// Author: Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Vpbe* Vpbe_ctor(Valist *alist, double ionConc, double ionRadius,
  double T, double soluteDiel, double solventDiel, double solventRadius) {

    /* Set up the structure */
    Vpbe *thee = VNULL;
    thee = Vmem_malloc(VNULL, 1, sizeof(Vpbe) );
    VASSERT( thee != VNULL);
    VASSERT( Vpbe_ctor2(thee, alist, ionConc, ionRadius, T, soluteDiel, 
      solventDiel, solventRadius) );

    return thee;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpbe_ctor2
//
// Purpose:  Construct the Vpbe object
//
// Notes:    Constructor broken into two parts for FORTRAN users.
//
// Returns:  1 if sucessful, 0 otherwise
//
// Author:   Nathan Baker and Mike Holst
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vpbe_ctor2(Vpbe *thee, Valist *alist, double ionConc, 
  double ionRadius, double T, double soluteDiel, double solventDiel, 
  double solventRadius) {

    int iatom;
    double atomRadius;
    Vatom *atom;
    double center[3] = {0.0, 0.0, 0.0};
    double disp[3], dist, radius, charge, xmin, xmax, ymin, ymax, zmin, zmax;
    double x, y, z;
    double nhash;
    const double N_A = 6.022045000e+23;
    const double e_c = 4.803242384e-10;
    const double k_B = 1.380662000e-16;
    const double pi  = 4. * VATAN(1.);

    /* Set up memory management object */
    thee->vmem = Vmem_ctor("APBS::VPBE");

    VASSERT(thee != VNULL);
    if (alist == VNULL) {
        Vnm_print(2, "Vpbe_ctor2: Got null pointer to Valist object!\n");
        return 0;
    }

    /* **** STUFF THAT GETS DONE FOR EVERYONE **** */
    /* Set pointers */
    thee->alist = alist;
    thee->paramFlag = 0;

    /* Set up Green's function oracle */
    thee->green = Vgreen_ctor(alist);

    /* Determine solute center */
    center[0] = thee->alist->center[0];
    center[1] = thee->alist->center[1];
    center[2] = thee->alist->center[2];
    thee->soluteCenter[0] = center[0];
    thee->soluteCenter[1] = center[1];
    thee->soluteCenter[2] = center[2];

    /* Determine solute length and charge*/
    radius = 0;
    atom = Valist_getAtom(thee->alist, 0);
    xmin = Vatom_getPosition(atom)[0];
    xmax = Vatom_getPosition(atom)[0];
    ymin = Vatom_getPosition(atom)[1];
    ymax = Vatom_getPosition(atom)[1];
    zmin = Vatom_getPosition(atom)[2];
    zmax = Vatom_getPosition(atom)[2];
    charge = 0;
    for (iatom=0; iatom<Valist_getNumberAtoms(thee->alist); iatom++) {
        atom = Valist_getAtom(thee->alist, iatom);
        atomRadius = Vatom_getRadius(atom);
        x = Vatom_getPosition(atom)[0];
        y = Vatom_getPosition(atom)[1];
        z = Vatom_getPosition(atom)[2];
        if ((x+atomRadius) > xmax) xmax = x + atomRadius;
        if ((x-atomRadius) < xmin) xmin = x - atomRadius;
        if ((y+atomRadius) > ymax) ymax = y + atomRadius;
        if ((y-atomRadius) < ymin) ymin = y - atomRadius;
        if ((z+atomRadius) > zmax) zmax = z + atomRadius;
        if ((z-atomRadius) < zmin) zmin = z - atomRadius;
        disp[0] = (x - center[0]);
        disp[1] = (y - center[1]);
        disp[2] = (z - center[2]);
        dist = (disp[0]*disp[0]) + (disp[1]*disp[1]) + (disp[2]*disp[2]); 
        dist = VSQRT(dist) + atomRadius;
        if (dist > radius) radius = dist;
        charge += Vatom_getCharge(Valist_getAtom(thee->alist, iatom));
    }
    thee->soluteRadius = radius;
    thee->soluteXlen = xmax - xmin;
    thee->soluteYlen = ymax - ymin;
    thee->soluteZlen = zmax - zmin;
    thee->soluteCharge = charge;

    /* Set parameters */
    thee->ionConc = ionConc;
    thee->ionRadius = ionRadius;
    thee->T = T;
    thee->soluteDiel = soluteDiel;
    thee->solventDiel = solventDiel;
    thee->solventRadius = solventRadius;

    /* Compute parameters: 
     *
     * kappa^2 = (8 pi N_A e_c^2) I_s / (1000 eps_w k_B T)
     * kappa   = 0.325567 * I_s^{1/2}   angstroms^{-1}
     * deblen  = 1 / kappa
     *         = 3.071564378 * I_s^{1/2}   angstroms
     * \bar{kappa}^2 = eps_w * kappa^2 
     * zmagic  = (4 * pi * e_c^2) / (k_B T)   (we scale the diagonal later)
     *         = 7046.528838
     */
    if ((thee->T == 0.) || (thee->ionConc == 0.)) {
        thee->xkappa  = 0.;
        thee->deblen  = 0.;
        thee->zkappa2 = 0.;
    } else {
        thee->xkappa  = VSQRT( thee->ionConc * 1.0e-16 *
            ((8.0 * pi * N_A * e_c*e_c) / 
            (1000.0 * thee->solventDiel * k_B * T))
        );
        thee->deblen  = 1. / thee->xkappa;
        thee->zkappa2 = thee->solventDiel * VSQR(thee->xkappa);
    }
    thee->zmagic  = ((4.0 * pi * e_c*e_c) / (k_B * thee->T)) * 1.0e+8;

    /* Compute accessibility objects:
     *   - Allow for extra room in the case of spline windowing 
     *   - Place some limits on the size of the hash table in the case of very
     *     large molecules
     */
    if (thee->ionRadius > thee->solventRadius) 
      radius = thee->ionRadius + MAX_SPLINE_WINDOW;
    else radius = thee->solventRadius + MAX_SPLINE_WINDOW;
    nhash = VPOW(8.0*(double)Valist_getNumberAtoms(thee->alist), 1.0/3.0);
    if (((int)nhash) < 3) nhash = 3;
    if (((int)nhash) > MAX_HASH_DIM) nhash = MAX_HASH_DIM;
    Vnm_print(0, "Vpbe_ctor2: Started constructing Vacc object with %d^3 hash table\n",
      (int)nhash); 
    thee->acc = Vacc_ctor(thee->alist, radius, (int)(nhash), (int)(nhash),
      (int)(nhash), VACC_SPHERE);
    Vnm_print(0, "Vpbe_ctor2: Done constructing Vacc object...\n"); 
    VASSERT(thee->acc != VNULL);

    return 1; 
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpbe_dtor
//
// Purpose:  Destroy the charge-simplex map.
// 
// Notes:    Since the grid manager and atom list were allocated outside of
//           the Vpbe routines, they are not destroyed.
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpbe_dtor(Vpbe **thee) {
    if ((*thee) != VNULL) {
        Vpbe_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(Vpbe), (void **)thee);
        (*thee) = VNULL;
    }
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpbe_dtor2
//
// Purpose:  Destroy the atom object
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpbe_dtor2(Vpbe *thee) { 
    Vacc_dtor(&(thee->acc));
    Vgreen_dtor(&(thee->green));
    Vmem_dtor(&(thee->vmem));
}


/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpbe_getCoulombEnergy1
//
// Purpose:  Perform an inefficient double sum to calculate the Coulombic
//           energy of a set of charges in a homogeneous dielectric (with
//           permittivity equal to the protein interior) and zero ionic
//           strength.  Result is returned in units of k_B T.  The sum can be
//           restriction to charges present in simplices of specified color
//           (pcolor); if (color == -1) no restrictions are used.
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vpbe_getCoulombEnergy1(Vpbe *thee) {

    int i, j, k, natoms;

    double dist, *ipos, *jpos, icharge, jcharge;
    double energy = 0.0;
    double eps, T;
    Vatom *iatom, *jatom;
    Valist *alist;
 
    VASSERT(thee != VNULL);
    alist = Vpbe_getValist(thee);
    VASSERT(alist != VNULL);
    natoms = Valist_getNumberAtoms(alist);
  
    /* Do the sum */ 
    for (i=0; i<natoms; i++) {
        iatom = Valist_getAtom(alist,i);
        icharge = Vatom_getCharge(iatom);
        ipos = Vatom_getPosition(iatom);
        for (j=i+1; j<natoms; j++) {
            jatom = Valist_getAtom(alist,j);
            jcharge = Vatom_getCharge(jatom);
            jpos = Vatom_getPosition(jatom);
            dist = 0;
            for (k=0; k<3; k++) dist += ((ipos[k]-jpos[k])*(ipos[k]-jpos[k]));
            dist = VSQRT(dist);
            energy = energy + icharge*jcharge/dist;
        }
    }

    /* Convert the result to J */
    T = Vpbe_getTemperature(thee);
    eps = Vpbe_getSoluteDiel(thee);
    energy = energy*Vunit_ec*Vunit_ec/(4*Vunit_pi*Vunit_eps0*eps*(1.0e-10));
   
    /* Scale by Boltzmann energy */
    energy = energy/(Vunit_kb*T);

    return energy;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpbe_memChk
//
// Purpose:  Returns the bytes used by the specified object
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vpbe_memChk(Vpbe *thee) {
   
    int memUse = 0;

    if (thee == VNULL) return 0;

    memUse = memUse + sizeof(Vpbe);
    memUse = memUse + Vacc_memChk(thee->acc);

    return memUse;
}
