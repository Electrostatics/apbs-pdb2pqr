/**
 *  @file    vfetk.c
 *  @ingroup Vfetk
 *  @author  Nathan Baker
 *  @brief   Class Vfetk methods
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
 * Copyright (c) 1999-2002.  Nathan A. Baker.  All Rights Reserved.
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
 *
 * @endverbatim
 */


#include "apbscfg.h"

#ifdef HAVE_MC_H

#include "apbs/vfetk.h"


/* ///////////////////////////////////////////////////////////////////////////
// Class Vfetk: Private method declaration
/////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vfetk: Inlineable methods
/////////////////////////////////////////////////////////////////////////// */
#if !defined(VINLINE_VFETK)

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vfetk_getGem
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Gem* Vfetk_getGem(Vfetk *thee) { 

   VASSERT(thee != VNULL);
   return thee->gm; 

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vfetk_getAM
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC AM* Vfetk_getAM(Vfetk *thee) { 

   VASSERT(thee != VNULL);
   return thee->am; 
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vfetk_getVpbe
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Vpbe* Vfetk_getVpbe(Vfetk *thee) { 

   VASSERT(thee != VNULL);
   return thee->pbe; 

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vfetk_getVcsm
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Vcsm* Vfetk_getVcsm(Vfetk *thee) { 

   VASSERT(thee != VNULL);
   return thee->csm; 

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vfetk_getAtomColor
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vfetk_getAtomColor(Vfetk *thee, int iatom) {

    int natoms;

    VASSERT(thee != VNULL);

    natoms = Valist_getNumberAtoms(Vpbe_getValist(thee->pbe));
    VASSERT(iatom < natoms);

    return Vatom_getPartID(Valist_getAtom(Vpbe_getValist(thee->pbe), iatom));
}
#endif /* if !defined(VINLINE_VFETK) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vfetk: Non-inlineable methods
/////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vfetk_ctor
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Vfetk* Vfetk_ctor(Vpbe *pbe, Gem *gm, AM *am) {

    /* Set up the structure */
    Vfetk *thee = VNULL;
    thee = Vmem_malloc(VNULL, 1, sizeof(Vfetk) );
    VASSERT(thee != VNULL);
    VASSERT(Vfetk_ctor2(thee, pbe, gm, am));

    return thee;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vfetk_ctor2
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vfetk_ctor2(Vfetk *thee, Vpbe *pbe, Gem *gm, AM *am) { 

    /* Make sure things have been properly initialized & store them */
    VASSERT(pbe != VNULL);
    VASSERT(pbe->alist != VNULL);
    VASSERT(pbe->acc != VNULL);
    VASSERT(gm != VNULL);
    VASSERT(am != VNULL);
    thee->pbe = pbe;
    thee->gm = gm;
    thee->am = am;

    /* Set up memory management object */
    thee->vmem = Vmem_ctor("APBS::VFETK");

    /* Set up charge-simplex map */
    thee->csm = Vcsm_ctor(Vpbe_getValist(thee->pbe), thee->gm);
    VASSERT(thee->csm != VNULL);
    Vcsm_init(thee->csm);

    return 1; 
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vfetk_dtor
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vfetk_dtor(Vfetk **thee) {
    if ((*thee) != VNULL) {
        Vfetk_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(Vfetk), (void **)thee);
        (*thee) = VNULL;
    }
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vfetk_dtor2
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vfetk_dtor2(Vfetk *thee) { 
    Vmem_dtor(&(thee->vmem));
    Vcsm_dtor(&(thee->csm));
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vfetk_getSolution
//
// Author:   Nathan Baker and Michael Holst
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double* Vfetk_getSolution(Vfetk *thee, int *length) { 

   int level, i;
   double *solution;
   double *theAnswer;
   Alg *alg;

   VASSERT(thee != VNULL);

   /* Get the max level from AM */
   level = AM_maxLevel(thee->am);
   /* Get the alg object at that level */
   alg = AM_alg(thee->am, level);

   /* Copy the solution into the w0 vector */
   Bvec_copy(alg->W[W_w0], alg->W[W_u]);
   /* Add the Dirichlet conditions */
   Bvec_axpy(alg->W[W_w0], alg->W[W_ud], 1.);
   /* Get the data from the Bvec */
   solution = Bvec_addr(alg->W[W_w0], 0);
   /* Get the length of the data from the Bvec */
   *length = Bvec_numRT(alg->W[W_w0]);
   /* Make sure that we got scalar data (only one block) for the solution
    * to the FETK */
   VASSERT(1 == Bvec_numB(alg->W[W_w0]));
   /* Allocate space for the returned vector and copy the solution into it */
   theAnswer = VNULL;
   theAnswer = Vmem_malloc(VNULL, *length, sizeof(double));
   VASSERT(theAnswer != VNULL);
   for (i=0; i<(*length); i++) theAnswer[i] = solution[i];
   
   return theAnswer;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vfetk_energy
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vfetk_energy(Vfetk *thee, int color, int nonlin) {

    double totEnergy = 0.0;
    double qfEnergy = 0.0;
    double dqmEnergy = 0.0;

    VASSERT(thee != VNULL);

    if (nonlin && (Vpbe_getBulkIonicStrength(thee->pbe) > 0.)) {
        Vnm_print(0, "Vfetk_energy:  calculating full PBE energy\n");
        Vnm_print(0, "Vfetk_energy:  bulk ionic strength = %g M\n",
          Vpbe_getBulkIonicStrength(thee->pbe));
        dqmEnergy = Vfetk_dqmEnergy(thee, color);
        Vnm_print(0, "Vfetk_energy:  dqmEnergy = %g kT\n", dqmEnergy);
        qfEnergy = Vfetk_qfEnergy(thee, color);
        Vnm_print(0, "Vfetk_energy:  qfEnergy = %g kT\n", qfEnergy);

        totEnergy = qfEnergy - dqmEnergy;
    } else {
        Vnm_print(0, "Vfetk_energy:  calculating only q-phi energy\n");
        dqmEnergy = Vfetk_dqmEnergy(thee, color);
        Vnm_print(0, "Vfetk_energy:  dqmEnergy = %g kT (NOT USED)\n", dqmEnergy);
        qfEnergy = Vfetk_qfEnergy(thee, color);
        Vnm_print(0, "Vfetk_energy:  qfEnergy = %g kT\n", qfEnergy);
        totEnergy = 0.5*qfEnergy;
    }

    return totEnergy;

}


/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vfetk_qfEnergy
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vfetk_qfEnergy(Vfetk *thee, int color) { 

   double *sol; int nsol;
   double charge;
   double phi[4], phix[4][3], *position;
   int iatom, natoms;
   int isimp, nsimps;
   int icolor;
   int ivert;
   SS *simp;
   double energy = 0.0;
   double uval;

   AM *am;

   VASSERT(thee != VNULL);
   am = thee->am;
    
   /* Get the finest level solution */
   sol= VNULL;
   sol = Vfetk_getSolution(thee, &nsol);
   VASSERT(sol != VNULL);

   /* Make sure the number of entries in the solution array matches the
    * number of vertices currently in the mesh */
   if (nsol != Gem_numVV(thee->gm)) {
      Vnm_print(2, "Vfetk_getLinearEnergy1: Number of unknowns in solution does not match\n");
      Vnm_print(2, "Vfetk_getLinearEnergy1: number of vertices in mesh!!!  Bailing out!\n");
      VASSERT(0);
   }
    
   /* Now we do the sum over atoms... */
   natoms = Valist_getNumberAtoms(thee->pbe->alist);
   for (iatom=0; iatom<natoms; iatom++) {
       /* Get atom information */
       icolor = Vfetk_getAtomColor(thee, iatom);
       charge = Vatom_getCharge(Valist_getAtom(thee->pbe->alist, iatom));
       position = Vatom_getPosition(Valist_getAtom(thee->pbe->alist, iatom));
       /* Check if this atom belongs to the specified partition */
       if ((color>=0) && (icolor<0)) {
           Vnm_print(2, "Vfetk_getLinearEnergy1: Atom colors not set!\n");
           VASSERT(0);
       }
       if ((icolor==color) || (color<0)) { 
           /* Loop over the simps associated with this atom */
           nsimps =  Vcsm_getNumberSimplices(thee->csm, iatom);
           /* Get the first simp of the correct color; we can use just one
            * simplex for energy evaluations, but not for force 
            * evaluations */
           for (isimp=0; isimp<nsimps; isimp++) {
               simp = Vcsm_getSimplex(thee->csm, isimp, iatom);
               /* If we've asked for a particular partition AND if the atom 
                * is our partition, then compute the energy */
               if ((SS_chart(simp)==color)||(color<0)) {
                   /* Get the value of each basis function evaluated at this
                    * point */
                   Gem_pointInSimplexVal(thee->gm, simp, position, phi, phix);
                   for (ivert=0; ivert<SS_dimVV(simp); ivert++) {
                       uval = sol[VV_id(SS_vertex(simp,ivert))];
                       energy += (charge*phi[ivert]*uval);
                   } /* end for ivert */
                   /* We only use one simplex of the appropriate color for
                    * energy calculations, so break here */
                   break;
               } /* endif (color) */
           } /* end for isimp */
       } 
   } /* end for iatom */

   /* Destroy the finest level solution */
   Vmem_free(VNULL, nsol, sizeof(double), (void **)&sol);
    
   /* Return the energy */
   return energy;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vfetk_dqmEnergy
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vfetk_dqmEnergy(Vfetk *thee, int color) { 

    return AM_evalJ(thee->am, AM_maxLevel(thee->am));

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vfetk_lnDet
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vfetk_lnDet(Vfetk *thee, int color, int flag) {

    Bmat *A;
    Zslu *slu;

    Alg *alg; 
    AM *am;
    int level, ip[10];
    int evalKey, tangKey, energyKey, residKey, massKey;
    double lndet, rp[10];

    VASSERT(thee != VNULL);
    am = thee->am;
    VASSERT(am != VNULL);

    if (color>=0) Vnm_print(2,"Vfetk_lnDet: color argument ignored!\n");

    /* Figure out key settings */
    evalKey = flag;
    energyKey = 1;
    residKey = 0;
    tangKey = 1;
    massKey = 0;
    if (flag == 0) tangKey = 0;
    else if (flag == 1) tangKey = 1;

    /* Get the max level from AM */
    level = AM_maxLevel(thee->am);

    /* Assemble the requested operator */
    Vnm_print(1, "Vfetk_lnDet: assembling operator...\n");
    AM_zeroA(thee->am, level);
    AM_init(am, level, W_f, 0.);
    AM_assem(am, level, evalKey, energyKey, residKey, tangKey, massKey, 
      W_u, W_ud, W_f, ip, rp);

    /* Au = A u */
    alg = AM_alg(thee->am, level);
    A = alg->A;
    Vnm_print(1, "Vfetk_lnDet: factoring matrix...\n");
    fflush(stdout);
    if (Bmat_sluFactor(A) == 0) {
        Vnm_print(2, "Vfetk_lnDet:  Error factoring matrix!\n");
        Vnm_print(2, "Vfetk_lnDet:  Last state = %d\n", A->AG->state);
        return 0.0;
    }
    slu = A->AG->slu;

    lndet = Zslu_lnDet(slu);

    return lndet;
}
    
    
/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vfetk_setAtomColors
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vfetk_setAtomColors(Vfetk *thee) {

#define VMAXLOCALCOLORSDONTREUSETHISVARIABLE 1024
    SS *simp;
    Vatom *atom;
    int i, natoms;

    VASSERT(thee != VNULL);

    natoms = Valist_getNumberAtoms(thee->pbe->alist);
    for (i=0; i<natoms; i++) {
        atom = Valist_getAtom(thee->pbe->alist, i);
        simp = Vcsm_getSimplex(thee->csm, 0, i);
        Vatom_setPartID(atom, SS_chart(simp));
    }

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vfetk_memChk
//
// Purpose:  Returns the bytes used by the specified object
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vfetk_memChk(Vfetk *thee) {
   
    int memUse = 0;

    if (thee == VNULL) return 0;

    memUse = memUse + sizeof(Vfetk);
    memUse = memUse + Vcsm_memChk(thee->csm);

    return memUse;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vfetk_genIcosGem
//  
// Purpose:   Given a non-NULL (but with 0 vertices) Gem object, create an
//            icosahedral domain with the outer boundary set to Dirichlet.
//
// Arguments: radius -- distance of outer vertices from center
//            center -- center of mesh
//  
// Returns:  0 if successful
//
// Author:   Tongye Shen and Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vfetk_genIcosGem(Gem *gm, double radius, double center[3]) {

    VV *vx;
    SS *sm;
    int i, j, theDim, theDimII;
    int topdata[80] = {1,6,12,0, 1,12,4,0, 1,4,8,0, 1,8,10,0, 1,10,6,0,
                       2,9,11,0, 2,5,9,0,  2,3,5,0, 2,7,3,0,  2,11,7,0,
                       7,6,10,0, 3,10,8,0, 5,8,4,0, 9,4,12,0, 11,12,6,0,
                       8,5,3,0,  4,9,5,0,  12,11,9,0, 6,7,11,0, 10,3,7,0 };
    double xyzdata[39] = { 0.000000e+00,  0.000000e+00,  0.000000e+00,
                           5.257311e-01,  0.000000e+00,  8.506508e-01,
                          -5.257311e-01,  0.000000e+00, -8.506508e-01,
                           5.257311e-01,  0.000000e+00, -8.506508e-01,
                           0.000000e+00,  8.506508e-01,  5.257311e-01,
                           0.000000e+00,  8.506508e-01, -5.257311e-01,
                           0.000000e+00, -8.506508e-01,  5.257311e-01,
                           0.000000e+00, -8.506508e-01, -5.257311e-01,
                           8.506508e-01,  5.257311e-01,  0.000000e+00,
                          -8.506508e-01,  5.257311e-01,  0.000000e+00,
                           8.506508e-01, -5.257311e-01,  0.000000e+00,
                          -8.506508e-01, -5.257311e-01,  0.000000e+00,
                          -5.257311e-01,  0.000000e+00,  8.506508e-01
                         };
    int numVV, chartV, numSS;
    int vnum, vtp, vtpI; 
    int fnum[3], ftp[3], ftpB[3];

    /* THIS ROUTINE IS BROKEN!!! */
    VASSERT(0);

    /* Make sure this is an empty Gem object */
    if (Gem_numVV(gm) != 0) {
        Vnm_print(2, "Vbnd_genIcosGem:  Error! Gem object has non-zero \
number of vertices!\n");
        return -1;
    }


    /* Set up parameters for 3D domain */
    theDim = 3;
    theDimII = 3;
    numVV = 13;
    numSS = 20;
    chartV = 0;

    /* Dump parameters for debugging info */
    Vnm_print(0, "Vbnd_genIcosGem:  radius = %g, center = (%g, %g, %g)\n",
      radius, center[0], center[1], center[2]);
    Vnm_print(0, "Vbnd_genIcosGem:  theDim=%d, theDimII=%d, numVV=%d, \
numSS=%d\n", theDim, theDimII, numVV, numSS);
    Vnm_print(0, "Vbnd_genIcosGem: Reseting manifold structures.\n");
    Gem_reset(gm, theDim, theDimII);

    
    /* Create the vertices */
    for (i=0; i<numVV; i++) {
        vx = Gem_createAndInitVV(gm);
        VV_setReality(vx, 0);
        VV_setDim(vx, theDim);
        VV_setClass(vx, 0);
        VV_setType(vx, 0);
        VV_setId(vx, i);
        VV_setChart(vx, chartV);
    
        /* set the vertex coordinates */
        VV_setCoord(vx, 0, radius*xyzdata[3*i]+center[0]);
        VV_setCoord(vx, 1, radius*xyzdata[3*i+1]+center[1]);
        VV_setCoord(vx, 2, radius*xyzdata[3*i+2]+center[2]);
    }
    
    /* Create the simplices */
    for (i=0; i<numSS; i++) {
    
        /* create the new simplex */
        sm = Gem_createAndInitSS(gm);
        SS_setReality(sm, 0);
        SS_setDim(sm, theDim);
        SS_setClass(sm, theDim);
        SS_setType(sm, 0);
        SS_setId(sm, i);
        SS_setChart(sm, 0);

        /* set the simplex face types and vertex labels */
        SS_setFaceType( sm, 0, 0 );
        SS_setFaceType( sm, 1, 0 );
        SS_setFaceType( sm, 2, 0 );
        SS_setFaceType( sm, 3, 1 );
        (gm->numBF)++;
        for (j=0; j<theDim+1; j++) {
            vx = Gem_VV(gm, topdata[4*i+j] );
            SS_setVertex( sm, j, vx );
        }

        /* calculate (our contribution to) vertex types from our face types */
        for (j=0; j<theDim+1; j++) {
            /* get the vertex in question */
            vnum = j;
            vx = SS_vertex( sm, vnum );
            /* get face numbers of two/three faces which touch vertex vnum */
            fnum[0] = vmapOV3[vnum][0];
            fnum[1] = vmapOV3[vnum][1];  
            fnum[2] = vmapOV3[vnum][2];  /* 2D: third face always interior */
            /* some shorthand notation... */
            vtp     = VV_type(vx);
            vtpI    = VINTERIOR( VV_type(vx) );
            ftp[0]  = SS_faceType(sm,fnum[0]);
            ftp[1]  = SS_faceType(sm,fnum[1]);
            ftp[2]  = SS_faceType(sm,fnum[2]);
            ftpB[0] = VBOUNDARY( SS_faceType(sm,fnum[0]) );
            ftpB[1] = VBOUNDARY( SS_faceType(sm,fnum[1]) );
            ftpB[2] = VBOUNDARY( SS_faceType(sm,fnum[2]) );
            /* if any of the faces are Boundary, then mark vertex Boundary */
            if ( ftpB[0] || ftpB[1] || ftpB[2] ) {

                /* deal with existing vertex type */
                if (vtpI) (gm->numBV)++;

                /* okay, determine max boundary flag (including vtp) */
                if (ftpB[0]) vtp = VMAX2(vtp,ftp[0]);
                if (ftpB[1]) vtp = VMAX2(vtp,ftp[1]);
                if (ftpB[2]) vtp = VMAX2(vtp,ftp[2]);

                /* set the type */
                VV_setType(vx, vtp);
            }

            /* build the ringed vertex datastructure */
        }
        SS_buildRing(sm);
    }

    /* create initial edge markings in the simplices */
    Gem_markEdges(gm);
    /* count v/e/f/s and check the mesh */
    Gem_countChk(gm);

    return 0;
}

#endif
