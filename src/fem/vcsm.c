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
// File:     vcsm.c
//
// Purpose:  Class Vcsm: methods. 
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */

#include "apbs/vcsm.h"

/* ///////////////////////////////////////////////////////////////////////////
// Class Vcsm: Inlineable methods
/////////////////////////////////////////////////////////////////////////// */
#if !defined(VINLINE_VCSM)
#endif /* if !defined(VINLINE_VCSM) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vcsm: Non-inlineable methods
/////////////////////////////////////////////////////////////////////////// */

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
    int iatom, jatom, isimp, jsimp;
    int gotSimp;
    /* Atomic information */
    Vatom *atom;
    double *position;
    /* Simplex/Vertex information */
    SS *simplex;

    VASSERT(thee != NULL);
    thee->natom = Valist_getNumberAtoms(thee->alist);
    thee->nsimp = thee->gm->numSS;
    VASSERT(thee->nsimp > 0);

    /* Set up the array of colors and initialize all to -1 */
    thee->colors = Vram_ctor(thee->natom, sizeof(int));
    VASSERT(thee->colors != VNULL);
    for (iatom=0; iatom<thee->natom; iatom++) thee->colors[iatom] = -1;

    /* Allocate and initialize space for the first dimensions of the 
     * simplex-charge map, the simplex array, and the counters */
    VASSERT( (thee->sqm = Vram_ctor(thee->nsimp, sizeof(int *))) != VNULL);
    VASSERT( (thee->nsqm = Vram_ctor(thee->nsimp, sizeof(int))) != VNULL);
    for (isimp=0; isimp<thee->nsimp; isimp++) (thee->nsqm)[isimp] = 0;

    /* Count the number of charges per simplex. */
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
    }

    /* Allocate the space for the simplex-charge map */
    for (isimp=0; isimp<thee->nsimp; isimp++) {
        if ((thee->nsqm)[isimp] > 0) {
            thee->sqm[isimp] = Vram_ctor((thee->nsqm)[isimp], sizeof(int));
            VASSERT(thee->sqm[isimp] != VNULL);
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

    thee->msimp = thee->nsimp;

    /* Allocate space for the charge-simplex map */
    VASSERT( (thee->qsm = Vram_ctor(thee->natom, sizeof(int *))) != VNULL);
    VASSERT( (thee->nqsm = Vram_ctor(thee->natom, sizeof(int))) != VNULL);
    for (iatom=0; iatom<thee->natom; iatom++) (thee->nqsm)[iatom] = 0;
    /* Loop through the list of simplices and count the number of times
     * each atom appears */
    for (isimp=0; isimp<thee->nsimp; isimp++) {
        for (iatom=0; iatom<thee->nsqm[isimp]; iatom++) {
            jatom = thee->sqm[isimp][iatom];
            thee->nqsm[jatom]++;
        }
    }
    /* Do a TIME-CONSUMING SANITY CHECK to make sure that each atom was
     * placed in at simplex */
    for (iatom=0; iatom<thee->natom; iatom++) {
        if (thee->nqsm[iatom] == 0) {
            Vnm_print(2,"Vcsm_init: atom %d not placed in simplex!\n");
            VASSERT(0);
        }
    }
    /* Allocate the appropriate amount of space for each entry in the
     * charge-simplex map and clear the counter for re-use in assignment */
    for (iatom=0; iatom<thee->natom; iatom++) {
        VASSERT(((thee->qsm)[iatom] = Vram_ctor((thee->nqsm)[iatom],
                                            sizeof(int)) ) != VNULL);
        thee->nqsm[iatom] = 0;
    }
    /* Assign the simplices to atoms */
    for (isimp=0; isimp<thee->nsimp; isimp++) {
        for (iatom=0; iatom<thee->nsqm[isimp]; iatom++) {
            jatom = thee->sqm[isimp][iatom];
            thee->qsm[jatom][thee->nqsm[jatom]] = isimp;
            thee->nqsm[jatom]++;
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
    if ((*thee) != VNULL) {
        Vcsm_dtor2(*thee);
        Vram_dtor((Vram **)thee, 1, sizeof(Vcsm));
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
VPUBLIC void Vcsm_dtor2(Vcsm *thee) { 
    int i;

    if ((thee != VNULL) && thee->initFlag) {

        for (i=0; i<thee->msimp; i++) {
            if (thee->nsqm[i] > 0) Vram_dtor((Vram **)&(thee->sqm[i]),
               thee->nsqm[i], sizeof(int));
        }
        for (i=0; i<thee->natom; i++) {
            if (thee->nqsm[i] > 0) Vram_dtor((Vram **)&(thee->qsm[i]),
               thee->nqsm[i], sizeof(int));
        }
        Vram_dtor((Vram **)&(thee->sqm), thee->msimp, sizeof(int *));
        Vram_dtor((Vram **)&(thee->nsqm), thee->msimp, sizeof(int));
        Vram_dtor((Vram **)&(thee->qsm), thee->natom, sizeof(int *));
        Vram_dtor((Vram **)&(thee->nqsm), thee->natom, sizeof(int));
        Vram_dtor((Vram **)&(thee->colors), thee->natom, sizeof(int));

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
// Routine:  Vcsm_getNumberSimplices
//
// Purpose:  Get the number of simplices associated with atom iatom
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vcsm_getNumberSimplices(Vcsm *thee, int iatom) {


   VASSERT(thee != VNULL);
   VASSERT(thee->initFlag);

   return (thee->nqsm)[iatom];

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vcsm_getSimplex
//
// Purpose:  Get simplex isimp associated with atom iatom
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC SS* Vcsm_getSimplex(Vcsm *thee, int isimp, int iatom) {


   VASSERT(thee != VNULL);
   VASSERT(thee->initFlag);

   return Vgm_SS(thee->gm, (thee->qsm)[iatom][isimp]);

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vcsm_getSimplexIndex
//
// Purpose:  Get simplex isimp associated with atom iatom
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vcsm_getSimplexIndex(Vcsm *thee, int isimp, int iatom) {


   VASSERT(thee != VNULL);
   VASSERT(thee->initFlag);

   return (thee->qsm)[iatom][isimp];

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
    int isimp, jsimp, iatom, jatom, atomID, simpID;
    int nsimps, gotMem;
    /* Object info */
    Vatom *atom;
    SS *simplex;
    double *position;
    /* Lists */
    int *qParent, nqParent;
    int **sqmNew, *nsqmNew;
    int *affAtoms, nAffAtoms;
    int *dnqsm, *nqsmNew, **qsmNew;

    VASSERT(thee != VNULL);
    VASSERT(thee->initFlag);

    /* If we don't have enough memory to accommodate the new entries, 
     * add more by doubling the existing amount */
    isimp = thee->nsimp + num - 1;
    gotMem = 0;
    while (!gotMem) {
        if (isimp > thee->msimp) {
            isimp = 2 * isimp;
            VASSERT( (thee->nsqm = 
              Vram_realloc((Vram **)&(thee->nsqm), thee->msimp, sizeof(int), 
              isimp)) != VNULL); 
            VASSERT( (thee->sqm = 
              Vram_realloc((Vram **)&(thee->sqm), thee->msimp, sizeof(int *), 
              isimp)) != VNULL); 
            thee->msimp = isimp;
        } else gotMem = 1;
    }
    /* Initialize the nsqm entires we just allocated */
    for (isimp = thee->nsimp; isimp<thee->nsimp+num-1 ; isimp++) {
       thee->nsqm[isimp] = 0;
    }
    
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
    for (isimp=0; isimp<num; isimp++) iatom += nsqmNew[isimp];
    if (iatom < nqParent) {
        Vnm_print(2,"Vcsm_update: Lost %d (of %d) atoms!\n", 
            nqParent - iatom, nqParent);
        VASSERT(0);
    }

    /* Update the QSM map using the old and new SQM lists */
    /* The affected atoms are those contained in the parent simplex; i.e.
     * thee->sqm[SS_id(simps[0])] */
    affAtoms = thee->sqm[SS_id(simps[0])];
    nAffAtoms = thee->nsqm[SS_id(simps[0])];
    /* Each of these atoms will go somewhere else; i.e., the entries in
     * thee->qsm are never destroyed and thee->nqsm never decreases.
     * However, it is possible that a subdivision could cause an atom to be
     * shared by two child simplices.  Here we record the change, if any,
     * in the number of simplices associated with each atom. */
    dnqsm = Vram_ctor(nAffAtoms, sizeof(int));
    VASSERT(dnqsm != VNULL);
    nqsmNew = Vram_ctor(nAffAtoms, sizeof(int));
    VASSERT(nqsmNew != VNULL);
    qsmNew = Vram_ctor(nAffAtoms, sizeof(int*));
    VASSERT(qsmNew != VNULL);
    for (iatom=0; iatom<nAffAtoms; iatom++) {
        dnqsm[iatom] = -1;
        atomID = affAtoms[iatom];
        for (isimp=0; isimp<num; isimp++) {
            for (jatom=0; jatom<nsqmNew[isimp]; jatom++) {
                if (sqmNew[isimp][jatom] == atomID) dnqsm[iatom]++;
            }
        }
        VASSERT(dnqsm[iatom] > -1);
    }
    /* Setup the new entries in the array */
    for (iatom=0;iatom<nAffAtoms; iatom++) {
        atomID = affAtoms[iatom];
        qsmNew[iatom] = Vram_ctor(dnqsm[iatom] + thee->nqsm[atomID], 
                                  sizeof(int));
        nqsmNew[iatom] = 0;
        VASSERT(qsmNew[iatom] != VNULL);
    }
    /* Fill the new entries in the array */
    /* First, do the modified entries */
    for (isimp=0; isimp<num; isimp++) {
        simpID = SS_id(simps[isimp]);
        for (iatom=0; iatom<nsqmNew[isimp]; iatom++) {
            atomID = sqmNew[isimp][iatom];
            for (jatom=0; jatom<nAffAtoms; jatom++) {
                if (atomID == affAtoms[jatom]) break;
            }
            if (jatom < nAffAtoms) {
                qsmNew[jatom][nqsmNew[jatom]] = simpID;
                nqsmNew[jatom]++;
            } 
        }
    }
    /* Now do the unmodified entries */
    for (iatom=0; iatom<nAffAtoms; iatom++) {
        atomID = affAtoms[iatom];
        for (isimp=0; isimp<thee->nqsm[atomID]; isimp++) {
            for (jsimp=0; jsimp<num; jsimp++) {
                simpID = SS_id(simps[jsimp]);
                if (thee->qsm[atomID][isimp] == simpID) break;
            }
            if (jsimp == num) {
                qsmNew[iatom][nqsmNew[iatom]] = thee->qsm[atomID][isimp];
                nqsmNew[iatom]++;
            }
        }
    }

    /* Replace the existing entries in the table.  Do the QSM entires
     * first, since they require affAtoms = thee->sqm[simps[0]] */
    for (iatom=0; iatom<nAffAtoms; iatom++) {
        atomID = affAtoms[iatom]; 
        Vram_dtor((Vram **)&(thee->qsm[atomID]), thee->nqsm[atomID], 
          sizeof(int));
        thee->qsm[atomID] = qsmNew[iatom];
        thee->nqsm[atomID] = nqsmNew[iatom];
    }
    for (isimp=0; isimp<num; isimp++) {
        simpID = SS_id(simps[isimp]);
        if (thee->nsqm[simpID] > 0) Vram_dtor((Vram **)&(thee->sqm[simpID]),
          thee->nsqm[simpID], sizeof(int));
        thee->sqm[simpID] = sqmNew[isimp];
        thee->nsqm[simpID] = nsqmNew[isimp];
    }

    Vram_dtor((Vram **)&sqmNew, num, sizeof(int *));
    Vram_dtor((Vram **)&nsqmNew, num, sizeof(int));
    Vram_dtor((Vram **)&qsmNew, nAffAtoms, sizeof(int *));
    Vram_dtor((Vram **)&nqsmNew, nAffAtoms, sizeof(int));
    Vram_dtor((Vram **)&dnqsm, nAffAtoms, sizeof(int));


    return 1;


}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vcsm_memChk
//
// Purpose:  Return number of bytes used by this object
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vcsm_memChk(Vcsm *thee)
{
    int i;
    int memUse = 0;

    VASSERT(thee != VNULL);

    memUse = memUse + sizeof(Vcsm);
    /* thee->sqm first dimension */
    memUse = memUse + (thee->msimp)*sizeof(int *);
    /* thee->sqm second dimension */
    for (i=0; i<thee->nsimp; i++) 
      memUse = memUse + (thee->nsqm[i])*sizeof(int);
    /* thee->nsqm */
    memUse = memUse + (thee->msimp)*sizeof(int);
    /* thee->qsm first dimension */
    memUse = memUse + (thee->natom)*sizeof(int *);
    /* thee->qsm second dimension */
    for (i=0; i<thee->natom; i++)
      memUse = memUse + (thee->nqsm[i])*sizeof(int);
    /* thee->nqsm */
    memUse = memUse + (thee->natom)*sizeof(int);
    /* thee->colors */
    memUse = memUse + (thee->natom)*sizeof(int);

    return memUse;
}

