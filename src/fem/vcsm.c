/* ///////////////////////////////////////////////////////////////////////////
/// MC = < Manifold Code >
/// Copyright (C) 1993  Michael Holst
/// Copyright (C) 1999  Contributing authors listed in the code documentation
///                     retain the copyrights to their code contributions.
///
/// The program MC is copyrighted software, and there are restrictions on its
/// use, modification, and distribution; please refer to the document COPYING
/// which accompanies this software.
///
/// This program is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
/// See the document COPYING for more details.
///
/// You should have received a copy of the COPYING file with this program;
/// if not, write to the author Michael Holst at:  mholst@math.ucsd.edu
///
/// rcsid="$Id$"
/// /////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// File:     vcsm.c
//
// Purpose:  Class Vcsm: methods. 
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */

#include "mc/vcsm.h"

/* ///////////////////////////////////////////////////////////////////////////
// Class Vcsm: Inlineable methods
/////////////////////////////////////////////////////////////////////////// */
#if !defined(VINLINE_VCSM)
#endif /* if !defined(VINLINE_VCSM) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vcsm: Non-inlineable methods
//
/////////////////////////////////////////////////////////////////////////// */
VPRIVATE void Vcsm_freeArrays(Vcsm *thee);

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vcsm_ctor
//
// Purpose:  Construct the charge-vertex map, assign atoms to vertices,
//           and assign vertices to atoms
//
// Notes:    The initial mesh must be sufficiently coarse for the
//           assignment procedures to be efficient
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Vcsm* Vcsm_ctor(Valist *alist, Vgm *gm) {

    /* Set up the structure */
    Vcsm *thee = VNULL;
    thee = (Vcsm*)calloc( 1, sizeof(Vcsm) );
    VASSERT( thee != VNULL);
    VASSERT( Vcsm_ctor2(thee, alist, gm));

    return thee;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vcsm_ctor2
//
// Purpose:  Construct the Vcsm object
//
// Notes:    Constructor broken into two parts for FORTRAN users.
//
// Returns:  1 if sucessful, 0 otherwise
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vcsm_ctor2(Vcsm *thee, Valist *alist, Vgm *gm) { 
 
    VASSERT( thee != VNULL );

    /* Set up the atom list and grid manager */
    if( alist == VNULL) {
        Vnm_print(2,"Vcsm_ctor2: got null pointer to Valist object!\n");
        return 0;
    }
    thee->alist = alist;
    if( gm == VNULL) {
        Vnm_print(2,"Vcsm_ctor2: got a null pointer to the Vgm object!\n");
        return 0;
    }
    thee->gm = gm;
   
    thee->initFlag = 0;
    return 1;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vcsm_init
//
// Purpose:  Initialize the charge-simplex map with mesh and charge data.
// 
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vcsm_init(Vcsm *thee) {
 
    /* Counters */
    int iatom, isimp, jsimp;
    int nsimps, gotSimp;
    /* Atomic information */
    Vatom *atom;
    double *position, charge, nsimpsf;
    /* Simplex/Vertex information */
    SS *simplex;

    VASSERT(thee != NULL);
    thee->natom = Valist_getNumberAtoms(thee->alist);
    thee->nsimp = thee->gm->numSS;
    VASSERT(thee->nsimp > 0);

    /* Allocate and initialize space for the first dimensions of the 
     * simplex-charge map, the simplex array, and the counters */
    VASSERT( (thee->sqm = Vram_ctor(thee->nsimp, sizeof(int *))) != VNULL);
    VASSERT( (thee->nsqm = Vram_ctor(thee->nsimp, sizeof(int))) != VNULL);
    for (isimp=0; isimp<thee->nsimp; isimp++) (thee->nsqm)[isimp] = 0;

    /* Count the number of charges per simplex.
     * IF AN ATOM IS IN MORE THAN ONE SIMPLEX, COUNT THE TOTAL NUMBER OF
     * SIMPLICES IT RESIDES IN AND DIVIDE THE ATOMIC CHARGE BY THAT NUMBER.
     * WE ASSUME THAT SIMPLICES ARE NEVER UNREFINED, SO ONCE AN ATOM'S
     * CHARGE IS DIVIDED IT WILL NEVER BE REINTEGRATED */
    for (iatom=0; iatom<thee->natom; iatom++) {
        atom = Valist_getAtom(thee->alist, iatom);
        position = Vatom_getPosition(atom);
        gotSimp = 0;
        for (isimp=0; isimp<thee->nsimp; isimp++) {
            simplex = Vgm_SS(thee->gm, isimp);
            if (Vgm_pointInSimplex(thee->gm, simplex, position)) {
                (thee->nsqm)[isimp]++;
                gotSimp = 1;
             }
        }
        if (!gotSimp) {
            Vnm_print(2, "Vcsm_init: Atom #%d (%4.3f, %4.3f, %4.3f) was not located in a simplex!\n", 
                iatom, position[0], position[1], position[2]);
            Vnm_print(2, "Vcsm_init: Confirming...\n");
            for (isimp=0; isimp<thee->nsimp; isimp++) {
                simplex = Vgm_SS(thee->gm, isimp);
                if (Vgm_pointInSimplex(thee->gm, simplex, position)) {
                    Vnm_print(2, "Vcsm_init:     Atom %d IN simplex %d!!!\n", iatom, isimp);
                 } else Vnm_print(2, "Vcsm_init:     Atom %d not in simplex %d\n", iatom, isimp);
            }
            VASSERT(0);
        }
    }

    /* Allocate the space for the simplex-charge map */
    for (isimp=0; isimp<thee->nsimp; isimp++) {
        if ((thee->nsqm)[isimp] > 0) {
            VASSERT(((thee->sqm)[isimp] = Vram_ctor((thee->nsqm)[isimp],
                                            sizeof(int)) ) != VNULL);
        }
    }

    /* Finally, set up the map */
    for (isimp=0; isimp<thee->nsimp; isimp++) {
        jsimp = 0;
        simplex = Vgm_SS(thee->gm, isimp);
        for (iatom=0; iatom<thee->natom; iatom++) {
            atom = Valist_getAtom(thee->alist, iatom);
            position = Vatom_getPosition(atom);
            /* Check to see if the atom's in this simplex */
            if (Vgm_pointInSimplex(thee->gm, simplex, position)) {
                /* Assign the entries in the next vacant spot */
                (thee->sqm)[isimp][jsimp] = iatom;
                jsimp++;
            }
        }
    }

    thee->initFlag = 1;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vcsm_dtor
//
// Purpose:  Destroy the charge-simplex map.
// 
// Notes:    Since the grid manager and atom list were allocated outside of
//           the Vcsm routines, they are not destroyed.
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vcsm_dtor(Vcsm **thee) {
    Vcsm_freeArrays(*thee);
    if ((*thee) != VNULL) {
        Vcsm_dtor2(*thee);
        free( *thee );
        (*thee) = VNULL;
    }
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vcsm_dtor2
//
// Purpose:  Destroy the atom object
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vcsm_dtor2(Vcsm *thee) { ; }

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vcsm_freeArrays
//
// Purpose:  Frees the memory allocated to the map arrays
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPRIVATE void Vcsm_freeArrays(Vcsm *thee) {
    int i;

    if ((thee != VNULL) && thee->initFlag) {

        Vnm_print(0,"Vcsm_freeArrays: freeing sqm entries.\n"); 
        for (i=0; i<thee->nsimp; i++) {
            if (thee->nsqm[i] > 0) free(thee->sqm[i]);
        }
        Vnm_print(0,"Vcsm_freeArrays: freeing sqm.\n"); 
        free(thee->sqm);
        Vnm_print(0,"Vcsm_freeArrays: freeing nsqm.\n"); 
        free(thee->nsqm);

    }
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vcsm_getValist
//
// Purpose:  Get a pointer to the Valist (atom list) object
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Valist* Vcsm_getValist(Vcsm *thee) { 

   VASSERT(thee != VNULL);
   return thee->alist;

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vcsm_getVgm
//
// Purpose:  Get a pointer to the Vgm (grid manager) object
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Vgm* Vcsm_getVgm(Vcsm *thee) { 

   VASSERT(thee != VNULL);
   return thee->gm; 

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vcsm_getNumberAtoms
//
// Purpose:  Get the number of atoms associated with simplex isimp.  
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vcsm_getNumberAtoms(Vcsm *thee, int isimp) {

   VASSERT(thee != VNULL);
   VASSERT(thee->initFlag);
   return thee->nsqm[isimp];

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vcsm_getAtom
//
// Purpose:  Get the iatom-th atom associated with isimp-th simplex.  
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Vatom* Vcsm_getAtom(Vcsm *thee, int iatom, int isimp) {


   VASSERT(thee != VNULL);
   VASSERT(thee->initFlag);

   VASSERT(iatom < (thee->nsqm)[isimp]);
   return Valist_getAtom(thee->alist, (thee->sqm)[isimp][iatom]);

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vcsm_getAtomIndex
//
// Purpose:  Get the iatom-th atom associated with isimp-th simplex.
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vcsm_getAtomIndex(Vcsm *thee, int iatom, int isimp) {


   VASSERT(thee != VNULL);
   VASSERT(thee->initFlag);

   VASSERT(iatom < (thee->nsqm)[isimp]);
   return (thee->sqm)[isimp][iatom];

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vcsm_update
//
// Purpose:  Update the charge-vertex and vertex-charge maps after
//           refinement
//
// Args:     Takes list of simplices.  The first simplex is derived from
//           the parent simplex and therefore has the same ID.  The
//           remaining simplices are the children and should represent new
//           entries in the charge-simplex maps.
//
// Returns:  1 if successful, 0 otherwise
//
// Author:   Nathan Baker
///////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vcsm_update(Vcsm *thee, SS **simps, int num) {

    /* Counters */
    int isimp, jsimp, iatom, atomID, simpID;
    int nsimps, gotMem;
    /* Object info */
    Vatom *atom;
    SS *simplex;
    double *position, charge, nsimpsf;
    /* Lists */
    int *qParent; int nqParent;
    int **sqmNew; int *nsqmNew;

    VASSERT(thee != VNULL);
    VASSERT(thee->initFlag);

    /* If we don't have enough memory to accommodate the new entries, 
     * add more by doubling the existing amount */
    isimp = thee->nsimp + num - 1;
    gotMem = 0;
    while (!gotMem) {
        if (isimp > thee->msimp) {
            isimp = 2 * isimp;
            printf("Vcsm_update: going to allocate %d entries\n", isimp);
            VASSERT( (thee->nsqm = 
              realloc(thee->nsqm, isimp * sizeof(int))) != VNULL); 
            VASSERT( (thee->sqm = 
              realloc(thee->sqm, isimp * sizeof(int *))) != VNULL); 
            thee->msimp = isimp;
            printf("Vcsm_update: reallocated %d entries\n",thee->msimp);
        } else gotMem = 1;
    }
    /* Initialize the nsqm entires we just allocated */
    for (isimp = thee->nsimp; isimp<thee->nsimp+num-1 ; isimp++) 
      thee->nsqm[isimp] = 0;
    thee->nsimp = thee->nsimp + num - 1;

    /* There's a simple case to deal with:  if simps[0] didn't have a
     * charge in the first place */
    isimp = SS_id(simps[0]);
    if (thee->nsqm[isimp] == 0) {
        for (isimp=1; isimp<num; isimp++) {
            thee->nsqm[SS_id(simps[isimp])] = 0;
        }
        return 1;
    }

    /*
     * printf("Vcsm_update: Updating Simp %d and children\n", SS_id(simps[0]));
     */

    /* The more complicated case has occured; the parent simplex had one or
     * more charges.  First, generate the list of affected charges. */
    isimp = SS_id(simps[0]);
    nqParent = thee->nsqm[isimp];
    qParent = thee->sqm[isimp];

    VASSERT( (sqmNew = Vram_ctor(num, sizeof(int *))) != VNULL);
    VASSERT( (nsqmNew = Vram_ctor(num, sizeof(int))) != VNULL);
    for (isimp=0; isimp<num; isimp++) nsqmNew[isimp] = 0;

    /* Loop throught the affected atoms to determine how many atoms each
     * simplex will get.  
     * IF AN ATOM WILL BE ASSIGNED TO MORE THAN ONE SIMPLEX, IT'S CHARGE IS
     * DIVIDED BY THE TOTAL NUMBER OF SIMPLICES TO WHICH IT WILL BE
     * ASSIGNED.  WE MAKE NO PROVISIONS FOR UNREFINEMENT OF SIMPLICES, SO
     * THIS DIVISION IS PERMANENT. */
    for (iatom=0; iatom<nqParent; iatom++) {
        atomID = qParent[iatom];
        atom = Valist_getAtom(thee->alist, atomID);
        position = Vatom_getPosition(atom);
        nsimps = 0;

        for (isimp=0; isimp<num; isimp++) {
            simplex = simps[isimp];
            if (Vgm_pointInSimplex(thee->gm, simplex, position)) {
                nsqmNew[isimp]++;
            }
        }
    }

    /* Allocate the storage */
    for (isimp=0; isimp<num; isimp++) {
        if (nsqmNew[isimp] > 0) {
            sqmNew[isimp] = Vram_ctor(nsqmNew[isimp], sizeof(int));
            VASSERT(sqmNew[isimp] != VNULL);
        }
    }

    /* Assign charges to simplices */
    for (isimp=0; isimp<num; isimp++) {

        jsimp = 0;
        simplex = simps[isimp];

        /* Loop over the atoms associated with the parent simplex */
        for (iatom=0; iatom<nqParent; iatom++) {

            atomID = qParent[iatom];
            atom = Valist_getAtom(thee->alist, atomID);
            position = Vatom_getPosition(atom);
            if (Vgm_pointInSimplex(thee->gm, simplex, position)) {
                sqmNew[isimp][jsimp] = atomID;
                jsimp++;
            }
        }
    }

    /* Sanity check that we didn't lose any atoms... */
    iatom = 0;
    for (isimp=0; isimp<num; isimp++) {
        if (nsqmNew[isimp] > 0) {
            /*
             * printf("Vcsm_update: Simp %d has %d atoms.\n",
             *   SS_id(simps[isimp]), nsqmNew[isimp]);
             */
        } 
        iatom += nsqmNew[isimp];
    }
    if (iatom < nqParent) {
        Vnm_print(2,"Vcsm_update: Lost %d (of %d) atoms!\n", 
            nqParent - iatom, nqParent);
        VASSERT(0);
    }

    /* Replace the existing entries in the table */
    for (isimp=0; isimp<num; isimp++) {
        simpID = SS_id(simps[isimp]);
        if (thee->nsqm[simpID] > 0) free(thee->sqm[simpID]);
        thee->sqm[simpID] = sqmNew[isimp];
        thee->nsqm[simpID] = nsqmNew[isimp];
    }


    return 1;


}

