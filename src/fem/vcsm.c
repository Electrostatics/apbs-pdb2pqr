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
VPRIVATE Vram* Vcsm_realloc(Vram **thee, int num, int size, int newNum);

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
    thee = Vram_ctor( 1, sizeof(Vcsm) );
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

    /* Set up the pool of Vlink objects */
    thee->pool = Vpool_ctor(Valist_getNumberAtoms(alist));
   
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
    int iatom, isimp;
    int gotSimp;
    /* Atomic information */
    Vatom *atom;
    double *position;
    /* Simplex/Vertex information */
    SS *simplex;

    VASSERT(thee != NULL);
    thee->natom = Valist_getNumberAtoms(thee->alist);
    thee->nsimp = Vgm_numSS(thee->gm);
    VASSERT(thee->nsimp > 0);

    /* Allocate and initialize space for the first dimensions of the 
     * simplex-charge map, the simplex array, and the counters */
    VASSERT( (thee->sqm = Vram_ctor(thee->nsimp, sizeof(Vlink *))) != VNULL);
    for (isimp=0; isimp<thee->nsimp; isimp++) (thee->sqm)[isimp] = VNULL;

    /* Assign charges to simplices */
    for (iatom=0; iatom<thee->natom; iatom++) {

        atom = Valist_getAtom(thee->alist, iatom);
        position = Vatom_getPosition(atom);

        gotSimp = 0;
        for (isimp=0; isimp<thee->nsimp; isimp++) {
            simplex = Vgm_SS(thee->gm, isimp);
            if (Vgm_pointInSimplex(thee->gm, simplex, position)) {

                /* Start the linked list if it hasn't been initialized */
                if (thee->sqm[isimp] == VNULL) 
                  thee->sqm[isimp] = (Vlink *)Vpool_create(thee->pool, atom);

                /* Add link for this atom to list */
                Vlink_setNext(Vlink_getLast(thee->sqm[isimp]), 
                   (Vlink *)Vpool_create(thee->pool, atom));
                gotSimp = 1;
             }
        }
        if (!gotSimp) {
            Vnm_print(2, "Vcsm_init: Atom #%d (%4.3f,%4.3f,%4.3f) not in simplex!\n", 
                iatom, position[0], position[1], position[2]);
            VASSERT(0);
        }
    }

    thee->msimp = thee->nsimp;
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
VPUBLIC void Vcsm_dtor(Vcsm **thee) { Vcsm_dtor2(*thee); }

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vcsm_dtor2
//
// Purpose:  Destroy the atom object
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vcsm_dtor2(Vcsm *thee) { 
    VASSERT ( thee != NULL);
    Vpool_dtor(&(thee->pool));
    Vram_dtor((Vram **)thee, 1, sizeof(Vcsm));
    thee = VNULL;
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

   Vlink *link1;
   int natoms = 0;

   VASSERT(thee != VNULL);
   VASSERT(thee->initFlag);

   if (thee->sqm[isimp] != VNULL) {
       link1 = thee->sqm[isimp];
       while (link1 != VNULL) {
           natoms++;
           link1 = Vlink_getNext(link1);
       }
   }

   return natoms;

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vcsm_getAtom
//
// Purpose:  Get the iatom-th atom associated with isimp-th simplex.  
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Vatom* Vcsm_getAtom(Vcsm *thee, int iatom, int isimp) {

   Vlink *link;
   int i;

   VASSERT(thee != VNULL);
   VASSERT(thee->initFlag);
   VASSERT(thee->sqm[isimp] != VNULL);
   
   link = thee->sqm[isimp];
   for (i=0; i<iatom; i++) VASSERT((link = Vlink_getNext(link)) != VNULL);
   return (Vatom *)Vlink_getData(link);

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
    int isimp, simpID, gotSimp, gotMem;
    /* Object info */
    Vatom *atom;
    SS *simplex;
    double *position;
    /* Linked list variables */
    Vlink **sqmNew, *qParent, *link;

    VASSERT(thee != VNULL);
    VASSERT(thee->initFlag);

    /* If we don't have enough memory to accommodate the new entries, 
     * add more by doubling the existing amount */
    isimp = thee->nsimp + num - 1;
    gotMem = 0;
    while (!gotMem) {
        if (isimp > thee->msimp) {
            isimp = 2 * isimp;
            thee->sqm = Vcsm_realloc((Vram **)&(thee->sqm), thee->msimp, 
              sizeof(int *), isimp);
            VASSERT(thee->sqm != VNULL);
            thee->msimp = isimp;
        } else gotMem = 1;
    }

    /* Update the number of simplices in the map */
    thee->nsimp = thee->nsimp + num - 1;

    /* There's a simple case to deal with:  if simps[0] didn't have a
     * charge in the first place */
    isimp = SS_id(simps[0]);
    if (thee->sqm[isimp] == VNULL) {
        for (isimp=1; isimp<num; isimp++) 
          thee->sqm[SS_id(simps[isimp])] = VNULL;
        return 1;
    }
    /* The more complicated case has occured; the parent simplex had one or
     * more charges.  First, generate the list of affected charges. */
    isimp = SS_id(simps[0]);
    qParent = thee->sqm[isimp];

    /* Here the new array of linked lists */
    sqmNew = (Vlink **)Vram_ctor(num, sizeof(Vlink *));
    VASSERT(sqmNew != VNULL);
    for (isimp=0; isimp<num; isimp++) sqmNew[isimp] = VNULL;

    /* Assign atoms to simplices */
    link = qParent;
    while (link != VNULL) {
        atom = (Vatom *)Vlink_getData(link);
        position = Vatom_getPosition(atom);

        gotSimp = 0;
        for (isimp=0; isimp<num; isimp++) {
            simplex = simps[isimp];

            if (Vgm_pointInSimplex(thee->gm, simplex, position)) {
                if (sqmNew[isimp] == VNULL) {
                    sqmNew[isimp] = (Vlink *)Vpool_create(thee->pool, atom);
                } else {
                    Vlink_setNext(Vlink_getLast(sqmNew[isimp]), 
                       (Vlink *)Vpool_create(thee->pool, atom));
                }
                gotSimp = 1;
             }
        }
        if (!gotSimp) {
            Vnm_print(2, "Vcsm_init: Atom (%4.3f,%4.3f,%4.3f) not in simplex!\n", 
                position[0], position[1], position[2]);
            VASSERT(0);
        }

        link = Vlink_getNext(link);
    }


    /* Replace the existing entries in the table */
    for (isimp=0; isimp<num; isimp++) {
        simpID = SS_id(simps[isimp]);
        if (thee->sqm[simpID] != VNULL) Vpool_destroyList(thee->pool, 
          thee->sqm[simpID]);
        thee->sqm[simpID] = sqmNew[isimp];
    }
    Vram_dtor((Vram **)&sqmNew, num, sizeof(Vlink *));

    return 1;


}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vcsm_realloc
//
// Purpose:  A logged version of realloc (using this is usually a bad idea)
//
// Author:   Michael Holst
/////////////////////////////////////////////////////////////////////////// */
VPRIVATE Vram* Vcsm_realloc(Vram **thee, int num, int size, int newNum)
{
    Vram *tee = Vram_ctor(newNum, size);
    memcpy(tee, (*thee), size*VMIN2(num,newNum));
    Vram_dtor((Vram **)thee, num, size);                  
    return tee;                
}


