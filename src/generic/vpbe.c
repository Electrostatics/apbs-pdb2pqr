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

#include "apbs/vpbe.h"

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpbe: Inlineable methods
/////////////////////////////////////////////////////////////////////////// */
#if !defined(VINLINE_VPBE)
#endif /* if !defined(VINLINE_VPBE) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpbe: Non-inlineable methods
/////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpbe_ctor
//
// Purpose:  Construct the charge-vertex map, assign atoms to vertices,
//           and assign vertices to atoms
//
// Notes:    The initial mesh must be sufficiently coarse for the
//           assignment procedures to be efficient.  
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Vpbe* Vpbe_ctor(Valist *alist, Vgm *gm) {

    /* Set up the structure */
    Vpbe *thee = VNULL;
    thee = Vram_ctor( 1, sizeof(Vpbe) );
    VASSERT( thee != VNULL);
    VASSERT( Vpbe_ctor2(thee, alist, gm));

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
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vpbe_ctor2(Vpbe *thee, Valist *alist, Vgm *gm) { 

    int iatom;
    double atomRadius;
    Vatom *atom;
    double center[3] = {0.0, 0.0, 0.0};
    double disp[3], dist, radius, charge;

    VASSERT(thee != VNULL);
    if (alist == VNULL) {
        Vnm_print(1,"Vpbe_ctor2: Got null pointer to Valist object!\n");
        return 0;
    }
    if (gm == VNULL) {
        Vnm_print(1,"Vpbe_ctor2: Got null pointer to Vgm object!\n");
        return 0;
    }

    /* Set pointers */
    thee->alist = alist;
    thee->gm = gm;
    thee->paramFlag = 0;

    /* Set up charge-simplex map */
    VASSERT((thee->csm = Vcsm_ctor(thee->alist, thee->gm)) != VNULL);

    /* Determine solute center */
    for (iatom=0; iatom<Valist_getNumberAtoms(thee->alist); iatom++) {
        center[0] += Vatom_getPosition(Valist_getAtom(thee->alist, iatom))[0];
        center[1] += Vatom_getPosition(Valist_getAtom(thee->alist, iatom))[1];
        center[2] += Vatom_getPosition(Valist_getAtom(thee->alist, iatom))[2];
    }
    center[0] = center[0]/((double)(Valist_getNumberAtoms(thee->alist)));
    center[1] = center[1]/((double)(Valist_getNumberAtoms(thee->alist)));
    center[2] = center[2]/((double)(Valist_getNumberAtoms(thee->alist)));
    thee->soluteCenter[0] = center[0];
    thee->soluteCenter[1] = center[1];
    thee->soluteCenter[2] = center[2];

    /* Determine solute radius and charge*/
    radius = 0;
    charge = 0;
    for (iatom=0; iatom<Valist_getNumberAtoms(thee->alist); iatom++) {
        atom = Valist_getAtom(thee->alist, iatom);
        disp[0] = (Vatom_getPosition(atom)[0] - center[0]);
        disp[1] = (Vatom_getPosition(atom)[1] - center[1]);
        disp[2] = (Vatom_getPosition(atom)[2] - center[2]);
        atomRadius = Vatom_getRadius(atom);
        dist = (disp[0]*disp[0]) + (disp[1]*disp[1]) + (disp[2]*disp[2]); 
        dist = VSQRT(dist) + atomRadius;
        if (dist > radius) radius = dist;
        charge += Vatom_getCharge(Valist_getAtom(thee->alist, iatom));
    }
    thee->soluteRadius = radius;
    thee->soluteCharge = charge;

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
        Vcsm_dtor(&((*thee)->csm));
        if ((*thee)->paramFlag) {
            Vhash_dtor(&((*thee)->solvHash));
            Vhash_dtor(&((*thee)->ionHash));
        }
        Vpbe_dtor2(*thee);
        Vram_dtor((Vram *)thee, 1, sizeof(Vpbe) );
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
VPUBLIC void Vpbe_dtor2(Vpbe *thee) { ; }

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpbe_initialize
//
// Purpose:  Set up parameters, hash tables, and charge-simplex map
//
// Arguments: ionConc       = ionic strength in M
//            ionRadius     = ionic probe radius in A
//            T             = temperature in K
//            soluteDiel    = solute dielectric (unitless)
//            solventDiel   = solvent dielectric (unitless)
//            solventRadius = solvent radius in Angstroms
//
// Notes:  Here's the original function coments from Mike
// notes:
// ------
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
//    some reference numbers in delphi:
//
//       salt = 1.0, deblen * I_s^{1/2} = 3.047000 angstroms
//       salt = 0.1, deblen * I_s^{1/2} = 9.635460 angstroms
//
//    delphi:  (kappa=0.328191663 * I_s^{1/2}, must have computed with
//    T!=298)
//    -------
//
//       deblen        = 1 / kappa
//                     = 1 / (0.328191663 * I_s^{1/2})
//                     = 3.047 / I_s^{1/2}    angstroms
//
//       kappa         = 1 / deblen
//                     = 0.328191663 * I_s^{1/2}    angstroms^{-1}
//
//       \bar{kappa}^2 = eps_w * kappa^2    angstroms^{-2}
//
//       debfact       = \bar{kappa}^2 * h^2
//                     = eps_w / (deblen * scale)**2
//
//    mike:  (kappa=0.325567 * I_s^{1/2}, with T=298)
//    -----
//
//       deblen        = 1 / kappa
//                     = 1 / (0.325567 * I_s^{1/2})
//                     = 3.071564378 * I_s^{1/2}   angstroms
//
//       kappa         = 1 / deblen
//                     = 0.325567 * I_s^{1/2}   angstroms^{-1}
//
//       \bar{kappa}^2 = eps_w * kappa^2   angstroms^{-2}
//
//       zkappa2       = \bar{kappa}^2  angstroms^{-2}
//
//    notes on scaling for the charges:
//    ---------------------------------
//
//       delphi:  (the 7049.484 seems to correspond to T=297.875)
//       -------
//
//          zmagic  = (4 * pi * e_c^2) * scale / (k_B T)
//                     / 6  (for diag scale of laplacean)
//                  = (4 * pi * e_c^2) / (6 h k_B T)
//                  = 7049.484 / (6 h)
//                  = 1174.914 / h
//
//       mike:    (the 7046.528838 corresponds to T=298)
//       -----
//
//          zmagic  = (4 * pi * e_c^2) / (k_B T)   (we scale the diagonal
//          later)
//                  = 7046.528838
//
//       since the units are esu^2 / erg, when converting to
//       angstroms^{-2},
//       we multiply by 1.0e8, yielding the 7046.528838 above.
//
//       important note:  in delphi, rhs_sca is multiplied against
//       the charges to yield directly the integrated function values
//       at grid points.  this is correct since the charges are represented
//       as *delta functions*, the "h^3" term does NOT appear after the
//       volume integrals are done as a result.  thus, since our generic
//       discretization module integrates this as a normal function,
//       yielding the h^3 term, we must account for it in the function
//       (i.e., divide by h^3 here).  thus, our scaling is:
//
//          zmagic  = 7046.528838 / h^3.
//
//    notes on scale and meshsize in delphi:
//    --------------------------------------
//
//       scale = 1 / h.
//       to properly discretize the original pde on the domain, we
//       reconstruct the domain from the mesh size and the number of
//       points.
//
//          meshsize = h = 1 / scale
//          hx = hy = hz = h
//          xmin = ymin = zmin = 0.0
//          xmax = ymax = zmax = 65.0 * h
//
//
// Author:   Nathan Baker and Michael Holst
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpbe_initialize(Vpbe *thee, double ionConc, double ionRadius,
                    double T, double soluteDiel, double solventDiel, 
                    double solventRadius) {

    const double N_A = 6.022045000e+23;
    const double e_c = 4.803242384e-10;
    const double k_B = 1.380662000e-16;
    const double pi  = 4. * VATAN(1.);
 
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
        thee->zkappa2 = thee->solventDiel * VPOW(thee->xkappa,2.);
    }
    thee->zmagic  = ((4.0 * pi * e_c*e_c) / (k_B * thee->T)) * 1.0e+8;

    /* Compute atomic hash tables */
    VASSERT( (thee->solvHash = Vhash_ctor(thee->alist, thee->solventRadius, 
                                      110, 110, 110)) != VNULL);
    VASSERT( (thee->ionHash = Vhash_ctor(thee->alist, thee->ionRadius, 
                                      110, 110, 110)) != VNULL);

    /* Compute charge-simplex map */
    Vcsm_init(thee->csm);

    thee->paramFlag = 1;
}

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
// Routine:  Vpbe_getVgm
//
// Purpose:  Get a pointer to the Vgm (grid manager) object
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Vgm* Vpbe_getVgm(Vpbe *thee) { 

   VASSERT(thee != VNULL);
   return thee->gm; 

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpbe_getSolventHash
//
// Purpose:  Get a pointer to the Vhash (atomic hash table) object for solvent
//           accessiblity
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Vhash* Vpbe_getSolventHash(Vpbe *thee) { 

   VASSERT(thee != VNULL);
   VASSERT(thee->paramFlag);
   return thee->solvHash; 

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpbe_getIonHash
//
// Purpose:  Get a pointer to the Vhash (atomic hash table) object for ionic
//           accessibility
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Vhash* Vpbe_getIonHash(Vpbe *thee) { 

   VASSERT(thee != VNULL);
   VASSERT(thee->paramFlag);
   return thee->ionHash; 

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpbe_getVcsm
//
// Purpose:  Get a pointer to the Vcsm (charge-simplex map) object
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Vcsm* Vpbe_getVcsm(Vpbe *thee) { 

   VASSERT(thee != VNULL);
   return thee->csm; 

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

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpbe_getSolution
//
// Purpose:  Get the electrostatic potential (in units of kT/e) at the
//           finest level of the passed AM object as a (newly allocated) array
//           of doubles and store the length in *length.  You'd better destroy
//           the returned array later!
//
// Author:   Nathan Baker and Michael Holst
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double* Vpbe_getSolution(Vpbe *thee, AM *am, int *length) { 

   int level, i;
   double *solution;
   double *theAnswer;

   VASSERT(thee != VNULL);
   VASSERT(am != VNULL);
   VASSERT(thee->gm != VNULL);

   /* Get the max level from AM */
   level = AM_maxLevel(am);

   /* Copy the solution into the w0 vector */
   Bvec_copy(AM_alg(am, level)->W[W_w0], 
     AM_alg(am, level)->W[W_u]);
   /* Add the Dirichlet conditions */
   Bvec_axpy(AM_alg(am, level)->W[W_w0], 
     AM_alg(am, level)->W[W_ud], 1.);
   /* Get the data from the Bvec */
   solution = Bvec_data(AM_alg(am, level)->W[W_w0], 0);
   /* Get the length of the data from the Bvec */
   *length = Bvec_numRT(AM_alg(am, level)->W[W_w0]);
   /* Make sure that we got scalar data (only one block) for the solution
    * to the PBE */
   VASSERT(1 == Bvec_numB(AM_alg(am, level)->W[W_w0]));
   /* Allocate space for the returned vector and copy the solution into it */
   theAnswer = VNULL;
   theAnswer = Vram_ctor(*length, sizeof(double));
   VASSERT(theAnswer != VNULL);
   for (i=0; i<*length; i++) theAnswer[i] = solution[i];
   
   return theAnswer;
}


/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpbe_getLinearEnergy1
//
// Purpose:  using the solution at the finest mesh level, get the 
//           electrostatic energy using the free energy functional for the 
//           linearized Poisson-Boltzmann equation without removing any 
//           self-interaction terms (i.e., removing the reference state of
//           isolated charges present in an infinite dielectric continuum with 
//           the same relative permittivity as the interior of the protein).
//           In other words, we calculate
//             \[ G = \frac{1}{2} \sum_i q_i u(r_i) \]
//           and return the result in units of $k_B T$.  The argument color
//           allows the user to control the partition on which this energy
//           is calculated; if (color == -1) no restrictions are used.
//           The solution is obtained from the finest level of the passed AM
//           object, but atomic data from the Vpbe object is used to
//           calculate the energy
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vpbe_getLinearEnergy1(Vpbe *thee, AM *am, int color) { 

   double *sol; int nsol;
   double charge;
   double phi[4], phix[4][3], *position;
   int iatom, natoms;
   int isimp, nsimps;
   int ivert;
   SS *simp;
   double energy = 0.0;
   double uval;

   VASSERT(thee != VNULL);
   VASSERT(am != VNULL);
   VASSERT(thee->gm != VNULL);
   VASSERT(thee->alist != VNULL);
   VASSERT(thee->csm != VNULL);

   /* Get the finest level solution */
   sol= VNULL;
   sol = Vpbe_getSolution(thee, am, &nsol);
   VASSERT(sol != VNULL);

   /* Make sure the number of entries in the solution array matches the
    * number of vertices currently in the mesh */
   if (nsol != Vgm_numVV(thee->gm)) {
      Vnm_print(2, "Vpbe_getLinearEnergy1: Number of unknowns in solution does not match\n");
      Vnm_print(2, "Vpbe_getLinearEnergy1: number of vertices in mesh!!!  Bailing out!\n");
      VASSERT(0);
   }

   /* Now we do the sum over atoms... */
   natoms = Valist_getNumberAtoms(thee->alist);
   for (iatom=0; iatom<natoms; iatom++) {
       charge = Vatom_getCharge(Valist_getAtom(thee->alist, iatom));
       position = Vatom_getPosition(Valist_getAtom(thee->alist, iatom));
       /* Loop over the simps associated with this atom */
       nsimps =  Vcsm_getNumberSimplices(thee->csm, iatom);
       /* Don't overcount q*phi contributions if found in more than one
        * simplex.  Since the solution is piecewise linear (i.e.,
        * continuous), FOR THE ENERGY ONLY, we should be able to integrate
        * the delta function in only one of the simplices.  Note that this
        * is a completely different story for the forces...
        */
       for (isimp=0; isimp<nsimps; isimp++) {
           simp = Vcsm_getSimplex(thee->csm, isimp, iatom);
           /* Loop through the list of simplices to see if any belong to
            * our partition */
           if ((SS_chart(simp)==color)||(color<0)) {
               /* Get the value of each basis function evaluated at this
                * point */
               Vgm_pointInSimplexVal(thee->gm, simp, position, phi, phix);
               for (ivert=0; ivert<SS_dimVV(simp); ivert++) {
                   uval = sol[VV_id(SS_vertex(simp,ivert))];
                   energy += (charge*phi[ivert]*uval);
               } /* end for ivert */
               break;
           } /* endif (color) */
       } /* end for isimp */
   } /* end for iatom */

   /* Destroy the finest level solution */
   Vram_dtor((Vram **)&sol, nsol, sizeof(double));

   /* Return the energy */
   return 0.5*energy;
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
