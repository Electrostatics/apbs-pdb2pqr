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
    thee = (Vcsm*)Vram_ctor( 1, sizeof(Vcsm) );
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
    int iatom, isimp;
    int gotSimp;
    char setid[80];
    /* Atomic information */
    Vatom *atom;
    double *position;
    /* Simplex/Vertex information */
    SS *simplex;

    VASSERT(thee != NULL);
    Vnm_tstart(65, "charge-simplex map init");
    thee->natom = Valist_getNumberAtoms(thee->alist);

    /* Allocate a Vset (dynamic array) which contains pointers to
     * additional Vset objects (the rows of this array) */
    Vnm_print(0, "Vcsm_init: Constructing initial sqm object for %d simps\n", 
      Vgm_numSS(thee->gm));
    thee->sqm = Vset_ctor("CSM.SQM", sizeof(Vset *), VMAX_OBJECTS);
    for (isimp=0; isimp<Vgm_numSS(thee->gm); isimp++) {
        sprintf(setid, "CSM.SQM[%d]", isimp);
        /* Obfuscated C at its best... Allocate a new (Vset *) entry for
         * integer objects in 
         * the dynamic array.  This gets returned as a (Vset **).
         * Dereference this and use it as the pointer to a new Vset object
         */
        *(Vset **)Vset_create(thee->sqm) = \
           Vset_ctor(setid, sizeof(int), VMAX_OBJECTS);
    }

    /* Place charges in simplices...  */
    Vnm_print(0, "Vcsm_init: Placing charges in simplices.\n");
    for (iatom=0; iatom<thee->natom; iatom++) {
        /* Get the atomic position */
        atom = Valist_getAtom(thee->alist, iatom);
        position = Vatom_getPosition(atom);

        gotSimp = 0;
        for (isimp=0; isimp<Vset_num(thee->sqm); isimp++) {
            simplex = Vgm_SS(thee->gm, isimp);
            if (Vgm_pointInSimplex(thee->gm, simplex, position)) {
                *(int *)Vset_create(*(Vset **)Vset_access(thee->sqm,
                    isimp)) = iatom;
                gotSimp = 1;
             }
        }

        if (!gotSimp) {
            Vnm_print(2, "Vcsm_init: Atom #%d", iatom);
            Vnm_print(2, " (%4.3f, %4.3f, %4.3f)", position[0], position[1], 
              position[2]);
            Vnm_print(2, " was not located in a simplex!\n"); 
            VASSERT(0);
        }
    }

    Vnm_tstop(65, "charge-simplex map init");
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
        for (i=0; i<Vset_num(thee->sqm); i++) 
          Vset_dtor((Vset **)Vset_access(thee->sqm, i));
        Vnm_print(0,"Vcsm_freeArrays: freeing sqm.\n"); 
        Vset_dtor(&(thee->sqm));

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
  
   Vset *set;

   VASSERT(thee != VNULL);
   VASSERT(thee->initFlag);
   set = *(Vset **)Vset_access(thee->sqm, isimp);
   return Vset_num(set);

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vcsm_getAtom
//
// Purpose:  Get the iatom-th atom associated with isimp-th simplex.  
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Vatom* Vcsm_getAtom(Vcsm *thee, int iatom, int isimp) {

   Vset *set;

   VASSERT(thee != VNULL);
   VASSERT(thee->initFlag);

   set = *(Vset **)Vset_access(thee->sqm, isimp);
   VASSERT(iatom < Vset_num(set));
 
   return Valist_getAtom(thee->alist, *Vset_access(set, iatom));

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vcsm_getAtomIndex
//
// Purpose:  Get the iatom-th atom associated with isimp-th simplex.
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vcsm_getAtomIndex(Vcsm *thee, int iatom, int isimp) {

   Vset *set;

   VASSERT(thee != VNULL);
   VASSERT(thee->initFlag);

   set = *(Vset **)Vset_access(thee->sqm, isimp);
   VASSERT(iatom < Vset_num(set));
 
   return *Vset_access(set, iatom);

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
    Vatom *atom;
    SS *simplex;
    double *position;
    char setid[80];
    Vset *qParent;
    Vset **sqmNew;

    VASSERT(thee != VNULL);
    VASSERT(thee->initFlag);

    /* Get the number new simplices and initialize */
    isimp = Vset_num(thee->sqm) + num - 1;

    Vnm_print(0, "Vcsm_init: Adding %d entries to CSM (%d entries)\n",
       num, Vset_num(thee->sqm));
    for (isimp = 0; isimp<num; isimp++) {
        sprintf(setid, "CSM.SQM[%d]", isimp);
        /* Obfuscated C at its best... Allocate a new (Vset *) entry for
         * integer objects in 
         * the dynamic array.  This gets returned as a (Vset **).
         * Dereference this and use it as the pointer to a new Vset object
         */
        *(Vset **)Vset_create(thee->sqm) = \
           Vset_ctor(setid, sizeof(int), VMAX_OBJECTS);
    }

    /* There's a simple case to deal with:  if simps[0] didn't have a
     * charge in the first place */
    isimp = SS_id(simps[0]);
    if (Vset_num(*(Vset **)Vset_access(thee->sqm, isimp)) == 0) return 1;

    /* The more complicated case has occured; the parent simplex had one or
     * more charges.  First, generate the list of affected charges. */
    isimp = SS_id(simps[0]);
    qParent = *(Vset **)Vset_access(thee->sqm, isimp);

    VASSERT( (sqmNew = Vram_ctor(num, sizeof(Vset *))) != VNULL);
    for (isimp=0; isimp<num; isimp++) {
        sprintf(setid, "CSM.TMP[%d]", isimp);
        sqmNew[isimp] = Vset_ctor(setid, sizeof(int), VMAX_OBJECTS);
    }

    /* Assign charges to simplices */
    for (isimp=0; isimp<num; isimp++) {

        jsimp = 0;
        simplex = simps[isimp];

        /* Loop over the atoms associated with the parent simplex */
        for (iatom=0; iatom<Vset_num(qParent); iatom++) {

            atomID = *(int *)Vset_access(qParent, iatom);
            atom = Valist_getAtom(thee->alist, atomID);
            position = Vatom_getPosition(atom);
            if (Vgm_pointInSimplex(thee->gm, simplex, position)) {
                *(int *)Vset_create(sqmNew[isimp]) = atomID;
            }
        }
    }

    /* Sanity check that we didn't lose any atoms... */
    iatom = 0;
    for (isimp=0; isimp<num; isimp++) {
        iatom += Vset_num(sqmNew[isimp]);
    }
    if (iatom < Vset_num(qParent)) {
        Vnm_print(2,"Vcsm_update: Lost %d (of %d) atoms!\n", 
            Vset_num(qParent) - iatom, Vset_num(qParent));
        VASSERT(0);
    }

    /* Replace the existing entries in the table */
    for (isimp=0; isimp<num; isimp++) {
        simpID = SS_id(simps[isimp]);
        Vset_dtor((Vset **)Vset_access(thee->sqm, simpID));
        *(Vset **)Vset_access(thee->sqm, simpID) = sqmNew[isimp];
    }

    Vram_dtor((Vram **)&sqmNew, num, sizeof(Vset *));

    Vnm_tstop(65, "charge-simplex map update");
    return 1;
}

