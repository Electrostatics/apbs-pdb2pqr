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
 * Nathan A. Baker (baker@biochem.wustl.edu)
 * Dept. of Biochemistry and Molecular Biophysics
 * Washington University in St. Louis
 *
 * Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 2002.  Washington University in St. Louis.
 * All Rights Reserved.
 *
 * Portions Copyright (c) 1999-2002.  The Regents of the University of
 * California.  
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
#if 0
    Vcsm_init(thee->csm); /* Catch 22 */
#endif

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

   int i;
   double *solution;
   double *theAnswer;
   AM *am;

   VASSERT(thee != VNULL);

   /* Get the AM object */
   am = thee->am;
   /* Copy the solution into the w0 vector */
   Bvec_copy(am->w0, am->u);
   /* Add the Dirichlet conditions */
   Bvec_axpy(am->w0, am->ud, 1.);
   /* Get the data from the Bvec */
   solution = Bvec_addr(am->w0, 0);
   /* Get the length of the data from the Bvec */
   *length = Bvec_numRT(am->w0);
   /* Make sure that we got scalar data (only one block) for the solution
    * to the FETK */
   VASSERT(1 == Bvec_numB(am->w0));
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

    return AM_evalJ(thee->am);

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vfetk_lnDet
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vfetk_lnDet(Vfetk *thee, int color, int oflag, int mflag) {

    Bmat *A;
    Zslu *slu;

    AM *am;
    int ip[10];
    int evalKey, tangKey, energyKey, residKey, massKey;
    double lndet = 0, rp[10];

    VASSERT(thee != VNULL);
    am = thee->am;
    VASSERT(am != VNULL);

    if (color>=0) Vnm_print(2,"Vfetk_lnDet: color argument ignored!\n");

    /* Figure out key settings */
    evalKey = oflag;
    energyKey = 1;
    residKey = 0;
    tangKey = 1;
    massKey = 0;
    if (oflag == 0) tangKey = 0;
    else if (oflag == 1) tangKey = 1;

    /* Assemble the requested operator */
    Vnm_print(1, "Vfetk_lnDet: assembling operator...\n");
    Bmat_zero(am->A);
    Bvec_init(am->f, 0.);
    AM_assem(am, evalKey, energyKey, residKey, tangKey, massKey,
        am->u, am->ud, am->f, ip, rp);

    /* Au = A u */
    A = am->A;
    Vnm_print(1, "Vfetk_lnDet: factoring matrix...\n");
    fflush(stdout);
    switch(mflag) {
    case 0: /* Calculate using ZSLU */
        if (Bmat_sluFactor(A) == 0) {
            Vnm_print(2, "Vfetk_lnDet:  Error factoring matrix!\n");
            if (A->AG == VNULL) {
                Vnm_print(2, "Vfetk_lnDet:  NULL A->AG: ");
                Vnm_print(2, "Mike -- what the hell are you doing???\n");
            }
            Vnm_print(2, "Vfetk_lnDet:  Last state = %d\n", Mat_state(A->AG));
            return 0.0;
        }
        slu = A->AG->slu;

        lndet = Zslu_lnDet(slu);
        break;
    case 1:
        if (Bmat_choleskyFactor(A,1) == 0) {
            Vnm_print(2, "Vfetk_lnDet:  Error Cholesky factoring matrix!\n");
            Vnm_print(2, "              Go complain to Steve?!\n");
        }

        lndet = 2*Mat_lnDetDiag(A->AG);
        Mat_dtor( &(A->AG) );
        break;
    default:
        Vnm_print(2, "Vfetk_lnDet:  Unsupported flag, doing nothing!\n");
    }

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

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  mContribRecycle
//
// Purpose:  Set or add a value to a linked matrix entry list.
//
// Notes:    key==0 ==> Set the value
//           key==1 ==> Add the value
//
//           rclnk==VNULL  ==> Allocate a new link if necessary
//           rclnk!=VNULL  ==> If necessary use rclnk, returning rclnk->next
//
// Author:   Stephen Bond and Michael Holst
/////////////////////////////////////////////////////////////////////////// */
VPRIVATE void mContribRecycle(Vset *mtpool,
    LinkA **rclnk, int key, int *count, int i, int j, double val)
{
    int done;
    LinkA *curr, *mnew;

    mnew=VNULL;
    curr=(LinkA*)Vset_access(mtpool,i);
    VASSERT( curr != VNULL );

    /* we have to look through the row(col) */
    done=0;
    while (!done) {

        /* if first guy in row(col) is "blank", fill him with the info */
        if (curr->idx == -1) {
            (*count)++;
            curr->idx  = j;
            curr->val  = val;
            curr->next = VNULL;
            done = 1;

        /* we found the position; just add in the contribution */
        } else if (curr->idx == j) {
            if (key == 0)
                curr->val = val;
            else
                curr->val += val;
            done = 1;

        /* it is not in the list; insert new struct AFTER it */
        /* KEY PLAY: this always inserts into the head of list */
        } else if (curr->idx > j) {
            (*count)++;
            if ( (*rclnk) == VNULL) {
                mnew = (LinkA*)Vset_create(mtpool);
            } else {
                mnew = (*rclnk);
                (*rclnk) = mnew->next;
            }
            mnew->idx  = curr->idx;
            mnew->val  = curr->val;
            mnew->next = curr->next;
            curr->idx  = j;
            curr->val  = val;
            curr->next = mnew;
            done = 1;

        /* not found; no one left */
        } else if (curr->next == VNULL) {
            (*count)++;
            if ( (*rclnk) == VNULL) {
                mnew = (LinkA*)Vset_create(mtpool);
            } else {
                mnew = (*rclnk);
                (*rclnk) = mnew->next;
            }
            mnew->idx  = j;
            mnew->val  = val;
            mnew->next = VNULL;
            curr->next = mnew;
            done = 1;

        /* not found yet; still hope */
        } else {
            curr=curr->next;
        }
    }
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Bmat_choleskyCreate
//
// Purpose:  Create the global matrix from the blocks in CLN storage format.
//           This is required for preparing a single global matrix for input
//           to the Cholesky Factor routine.
//
//    Note:  It is assumed implicitly that the Bmat is symmetric, and hence
//           we only copy and store the lower triangular part!
//
// Author:   Stephen Bond
/////////////////////////////////////////////////////////////////////////// */
VPRIVATE void Bmat_choleskyCreate(Bmat *thee)
{
    int i, j, k, p, q, count, istart, jstart;
    Vset  *lnk;
    LinkA *mt;
    Mat *blk, *gmat;
    MATmirror mirror;
    MATformat format;
    /* Bound on the size of L */
    int maxnZ = Bmat_numRT(thee)*1000;
    double *diag, *offU, *offL;

    /* initialize the global matrix datastructure */
    gmat = Mat_ctor(thee->vmem, "AG", Bmat_numRT(thee), Bmat_numCT(thee));
    Mat_initStructure(gmat, CLN_FORMAT, IS_SYM, 0, VNULL);

    /* initialize/clear the dynamic array */
    gmat->lnkL = Vset_ctor( thee->vmem, "lnk", sizeof( LinkA ), maxnZ, 0 );
    Vset_reset( gmat->lnkL );
    lnk = gmat->lnkL;

    /* create an empty entry to start each row of global matrix */
    for (i=0; i<Bmat_numRT(thee); i++) {
        mt=(LinkA*)Vset_create(lnk);
        mt->idx  = -1;
        mt->val  = 0.;
        mt->next = VNULL;
    }

    /* now get the COL structure of the matrix */
    count = 0;

    istart = 0;
    for (p=0; p<thee->numB; p++) {
        jstart = 0;
        for (q=0; q<=p; q++) {
            mirror = thee->mirror[p][q];
            blk    = thee->AD[p][q];
            format = Mat_format(blk);
            diag   = blk->diag;
            offU   = blk->offU;
            offL   = blk->offL;
            if (mirror) {
                if (format == ROW_FORMAT) {
                    format = COL_FORMAT;
                } else if (format == COL_FORMAT) {
                    format = ROW_FORMAT;
                } else {
                    VASSERT(0);
                }
            }

            for (j=0; j<Mat_numR(blk); j++) {
                if (format == DRC_FORMAT) {
                    i = j;
                    mContrib(lnk,0,&count,jstart+j,istart+i,diag[i]);
                    for (k=blk->IA[j]; k<blk->IA[j+1]; k++) {
                        i = blk->JA[k];
                        mContrib(lnk,0,&count,jstart+j,istart+i,offL[k]);
                    }
                } else if (format == ROW_FORMAT) {
                    for (k=blk->IA[j]; k<blk->IA[j+1]; k++) {
                        i = blk->JA[k];
                        mContrib(lnk,0,&count,jstart+i,istart+j,blk->A[k]);
                    }
                } else if (format == COL_FORMAT) {
                    for (k=blk->IA[j]; k<blk->IA[j+1]; k++) {
                        i = blk->JA[k];
                        mContrib(lnk,0,&count,jstart+j,istart+i,blk->A[k]);
                    }
                } else { VASSERT(0); }

            }
            /* increment index for the column block */
            jstart += Bmat_numC(thee,q,q);
        }
        /* increment index for the row block */
        istart += Bmat_numR(thee,p,p);
    }
    gmat->numO = count;

    thee->AG = gmat;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Mat_choleskyFactor
//
// Purpose:  Create the sparse L L^T factors for the system matrix.
//
//    Note:  The matrix is overwriten with the result!
//
//  Method:  In OCTAVE pseudocode
//
//           for k = 1:n
//             A(k,k) = sqrt(A(k,k));
//             for i = (k+1):n;
//               A(i,k) = A(i,k)/A(k,k);
//             end
//             for j = (k+1):n
//               for i = j:n
//                 A(i,j) = A(i,j) - A(i,k)*A(j,k);
//               end
//             end
//           end
//
//
// Author:   Stephen Bond
/////////////////////////////////////////////////////////////////////////// */
VPRIVATE int Mat_choleskyFactor(Mat *thee, int flag)
{
    int count, k, n;
    int recyc_idx;
    double Akk;
    LinkA *mt, *mti, *mtj, *mt_recyc, *mt_tmp;

    Vset *lnk = thee->lnkL;

    VASSERT( Mat_format(thee) == CLN_FORMAT );
    VASSERT( Mat_sym(thee) == IS_SYM );
    VASSERT( Mat_numR(thee) == Mat_numC(thee) );

    /*  matrix has already been sparse factored */
    if ( Mat_state(thee) == FACTORED_STATE) {
        return 1;
    }

    /*  groovy, it hasn't been factored yet!  */
    n = Mat_numR(thee);
    count = Mat_numO(thee);

    switch( flag ) {
    case 0:  /* CASE 0:  Keep full factoring and No pivoting */
        for (k=0; k<n; k++) {
            mt = (LinkA*)Vset_access(lnk,k); /* A(k,k) */
            if ( mt->val > 0 && mt->idx >= 0 ) {
                mt->val = sqrt(mt->val);
                Akk = mt->val;
                for ( mti=mt->next; mti!=VNULL; mti=mti->next ) {
                    mti->val /= Akk; /* A(i,k)/A(k,k) */
                }

                for ( mtj=mt->next; mtj!=VNULL; mtj=mtj->next ) {
                    for ( mti=mtj; mti!=VNULL; mti=mti->next ) {
                        mContrib( lnk, 1, &count, mtj->idx, mti->idx,
                                  -1*mti->val*mtj->val );
                        /* A(i,j) -= A(i,k)*A(j,k) */
                    }
                }
            } else {
                return 0;
            }
        }
        break;
    case 1:  /* CASE 1:  Only get the diagonal correct */
        recyc_idx = 0;
        mt_recyc = VNULL;
        for (k=0; k<n; k++) {
            mt = (LinkA*)Vset_access(lnk,k); /* A(k,k) */
            if ( mt->val > 0 && mt->idx >= 0 ) {
                mt->val = sqrt(mt->val);
                Akk = mt->val;
                for ( mti=mt->next; mti!=VNULL; mti=mti->next ) {
                    mti->val /= Akk; /* A(i,k)/A(k,k) */
                }

                for ( mtj=mt->next; mtj!=VNULL; mtj=mtj->next ) {
                    for ( mti=mtj; mti!=VNULL; mti=mti->next ) {
                        while( (mt_recyc == VNULL) && (recyc_idx < k) ) {
                            mt_tmp = (LinkA*)Vset_access(lnk,recyc_idx);
                            mt_recyc = mt_tmp->next;
                            mt_tmp = VNULL;
                            recyc_idx++;
                        }
                        mContribRecycle( lnk, &(mt_recyc), 1, &count,
                            mtj->idx, mti->idx, -1*mti->val*mtj->val );
                        /* A(i,j) -= A(i,k)*A(j,k) */
                    }
                }
            } else {
                return 0;
            }
        }
        break;
    default:
        VASSERT( 0 );
    }

    thee->state = FACTORED_STATE;

    return 1;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Bmat_choleskyFactor
//
// Purpose:  Create the sparse L L^T factors for global matrix.
//
//    Note:  Will only work if A is symmetric positive definite!
//
// Author:   Stephen Bond
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Bmat_choleskyFactor(Bmat *thee, int flag)
{
    if ( thee->AG == VNULL ) {
        Bmat_choleskyCreate(thee);
    }

    if ( Mat_state(thee->AG) == FACTORED_STATE ) {
        return 1;
    } else {
        return Mat_choleskyFactor(thee->AG, flag);
    }
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Mat_lnDetDiag
//
// Purpose:  Calculate the log(det(A)), where A is diagonal or triangular.
//
// Notes:    Use another algorithm first to reduce to the triangular.
//
// Author:   Stephen Bond
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Mat_lnDetDiag(Mat *thee)
{
    int k,n;
    LinkA *mt;
    Vset *lnk;
    double lndet;

    VASSERT( Mat_numR(thee) == Mat_numC(thee) );
    n = Mat_numR(thee);

    lndet = 0.0;
    switch( Mat_format(thee) ) {
    case RLN_FORMAT:
    case CLN_FORMAT:
        lnk = thee->lnkL;
        for ( k=0; k<n; k++ ) {
            mt=(LinkA*)Vset_access(lnk,k);
            lndet += log(VABS(mt->val));
        }
        break;
    case DRC_FORMAT:
        for ( k=0; k<n; k++ ) {
            lndet += log(VABS(thee->diag[k]));
        }
        break;
    default:
        VASSERT(0); /* Currently Unsupported */
    }

    return lndet;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Bmat_printHB
//
// Purpose:  Prints a Bmat in sparse Harwell-Boeing format
//
//  Author:  Stephen Bond (following HB matrix user's guide)
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Bmat_printHB( Bmat *thee, char *fname )
{
    Mat *Ablock;
    MATsym pqsym;
    int i, j, jj;
    int *IA, *JA;
    double *D, *L, *U;
    FILE *fp;

    char mmtitle[72];
    char mmkey[] = {"8charkey"};
    int totc = 0, ptrc = 0, indc = 0, valc = 0;
    char mxtyp[] = {"RUA"}; /* Real Unsymmetric Assembled */
    int nrow = 0, ncol = 0, numZ = 0;
    int numZdigits = 0, nrowdigits = 0;
    int nptrline = 8, nindline = 8, nvalline = 5;
    char ptrfmt[] = {"(8I10)          "}, ptrfmtstr[] = {"%10d"};
    char indfmt[] = {"(8I10)          "}, indfmtstr[] = {"%10d"};
    char valfmt[] = {"(5E16.8)            "}, valfmtstr[] = {"%16.8E"};

    VASSERT( thee->numB == 1 );             /* HARDWIRE FOR NOW */
    Ablock = thee->AD[0][0];

    VASSERT( Mat_format( Ablock ) == DRC_FORMAT );  /* HARDWIRE FOR NOW */

    pqsym = Mat_sym( Ablock );

    if ( pqsym == IS_SYM ) {
        mxtyp[1] = 'S';
    } else if ( pqsym == ISNOT_SYM ) {
        mxtyp[1] = 'U';
    } else {
        VASSERT( 0 ); /* NOT VALID */
    }

    nrow = Bmat_numRT( thee ); /* Number of rows */
    ncol = Bmat_numCT( thee ); /* Number of cols */
    numZ = Bmat_numZT( thee ); /* Number of entries */

    nrowdigits = (int) (log( nrow )/log( 10 )) + 1;
    numZdigits = (int) (log( numZ )/log( 10 )) + 1;

    nptrline = (int) ( 80 / (numZdigits + 1) );
    nindline = (int) ( 80 / (nrowdigits + 1) );

    sprintf(ptrfmt,"(%dI%d)",nptrline,numZdigits+1);
    sprintf(ptrfmtstr,"%%%dd",numZdigits+1);
    sprintf(indfmt,"(%dI%d)",nindline,nrowdigits+1);
    sprintf(indfmtstr,"%%%dd",nrowdigits+1);

    ptrc = (int) ( ( (ncol + 1) - 1 ) / nptrline ) + 1;
    indc = (int) ( (numZ - 1) / nindline ) + 1;
    valc = (int) ( (numZ - 1) / nvalline ) + 1;

    totc = ptrc + indc + valc;

    sprintf( mmtitle, "Sparse '%s' Matrix - Harwell-Boeing Format - '%s'",
             thee->name, fname );

   /* Step 0:  Open the file for writing */

    fp = fopen( fname, "w" );
    if (fp == VNULL) {
        Vnm_print(2,"Bmat_printHB:  Ouch couldn't open file <%s>\n",fname);
        return;
    }

    /* Step 1:  Print the header information */

    fprintf( fp, "%-72s%-8s\n", mmtitle, mmkey );
    fprintf( fp, "%14d%14d%14d%14d%14d\n", totc, ptrc, indc, valc, 0 );
    fprintf( fp, "%3s%11s%14d%14d%14d\n", mxtyp, " ", nrow, ncol, numZ );
    fprintf( fp, "%-16s%-16s%-20s%-20s\n", ptrfmt, indfmt, valfmt, "6E13.5" );

    IA = Ablock->IA;
    JA = Ablock->JA;
    D = Ablock->diag;
    L = Ablock->offL;
    U = Ablock->offU;

    if ( pqsym == IS_SYM ) {

        /* Step 2:  Print the pointer information */

        for (i=0; i<(ncol+1); i++) {
            fprintf( fp, ptrfmtstr, Ablock->IA[i] + (i+1) );
            if ( ( (i+1) % nptrline ) == 0 ) {
                fprintf( fp, "\n" );
            }
        }

        if ( ( (ncol+1) % nptrline ) != 0 ) {
            fprintf( fp, "\n" );
        }

        /* Step 3:  Print the index information */

        j = 0;
        for (i=0; i<ncol; i++) {
            fprintf( fp, indfmtstr, i+1); /* diagonal */
            if ( ( (j+1) % nindline ) == 0 ) {
                fprintf( fp, "\n" );
            }
            j++;
            for (jj=IA[i]; jj<IA[i+1]; jj++) {
                fprintf( fp, indfmtstr, JA[jj] + 1 ); /* lower triangle */
                if ( ( (j+1) % nindline ) == 0 ) {
                    fprintf( fp, "\n" );
                }
                j++;
            }
        }

        if ( ( j % nindline ) != 0 ) {
            fprintf( fp, "\n" );
        }

        /* Step 4:  Print the value information */

        j = 0;
        for (i=0; i<ncol; i++) {
            fprintf( fp, valfmtstr, D[i] );
            if ( ( (j+1) % nvalline ) == 0 ) {
                fprintf( fp, "\n" );
            }
            j++;
            for (jj=IA[i]; jj<IA[i+1]; jj++) {
                fprintf( fp, valfmtstr, L[jj] );
                if ( ( (j+1) % nvalline ) == 0 ) {
                    fprintf( fp, "\n" );
                }
                j++;
            }
        }

        if ( ( j % nvalline ) != 0 ) {
            fprintf( fp, "\n" );
        }

    } else { /* ISNOT_SYM */

        VASSERT( 0 ); /* NOT CODED YET */
    }

    /* Step 5:  Close the file */
    fclose( fp );

}

#endif
