/**
 *  @file    vcsm.c
 *  @ingroup Vcsm
 *  @author  Nathan Baker
 *  @brief   Class Vcsm methods
 *  @version $Id$
 *  @attention
 *  @verbatim
 *
 * APBS -- Adaptive Poisson-Boltzmann Solver
 *
 * Nathan A. Baker (baker@biochem.wustl.edu)
 * Dept. of Biochemistry and Molecular Biophysics
 * Center for Computational Biology
 * Washington University in St. Louis
 *
 * Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 2002-2006.  Washington University in St. Louis.
 * All Rights Reserved.
 * Portions Copyright (c) 1999-2002.  The Regents of the University of
 * California.  
 * Portions Copyright (c) 1995.  Michael Holst.
 *
 * This file is part of APBS.
 *
 * APBS is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * APBS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with APBS; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 *
 * Linking APBS statically or dynamically with other modules is making a
 * combined work based on APBS. Thus, the terms and conditions of the GNU
 * General Public License cover the whole combination.
 * 
 * SPECIAL GPL EXCEPTION
 * In addition, as a special exception, the copyright holders of APBS
 * give you permission to combine the APBS program with free software
 * programs and libraries that are released under the GNU LGPL or with
 * code included in releases of ISIM, PMV, PyMOL, SMOL, VMD, and Vision.
 * Such combined software may be linked with APBS and redistributed together 
 * in original or modified form as mere aggregation without requirement that 
 * the entire work be under the scope of the GNU General Public License.
 * This special exception permission is also extended to any software listed
 * in the SPECIAL GPL EXCEPTION clauses by the PMG, FEtk, MC, or MALOC
 * libraries.
 * 
 * Note that people who make modified versions of APBS are not obligated
 * to grant this special exception for their modified versions; it is
 * their choice whether to do so. The GNU General Public License gives
 * permission to release a modified version without this exception; this
 * exception also makes it possible to release a modified version which
 * carries forward this exception.
 *
 * @endverbatim
 */


#include "apbscfg.h"

#if defined(HAVE_MC_H)
#include "apbs/vcsm.h"

/* Inlineable methods */
#if !defined(VINLINE_VCSM)

VPUBLIC Valist* Vcsm_getValist(Vcsm *thee) { 

   VASSERT(thee != VNULL);
   return thee->alist;

}

VPUBLIC int Vcsm_getNumberAtoms(Vcsm *thee, int isimp) {

   VASSERT(thee != VNULL);
   VASSERT(thee->initFlag);
   return thee->nsqm[isimp];

}

VPUBLIC Vatom* Vcsm_getAtom(Vcsm *thee, int iatom, int isimp) {


   VASSERT(thee != VNULL);
   VASSERT(thee->initFlag);

   VASSERT(iatom < (thee->nsqm)[isimp]);
   return Valist_getAtom(thee->alist, (thee->sqm)[isimp][iatom]);

}

VPUBLIC int Vcsm_getAtomIndex(Vcsm *thee, int iatom, int isimp) {


   VASSERT(thee != VNULL);
   VASSERT(thee->initFlag);

   VASSERT(iatom < (thee->nsqm)[isimp]);
   return (thee->sqm)[isimp][iatom];

}

VPUBLIC int Vcsm_getNumberSimplices(Vcsm *thee, int iatom) {


   VASSERT(thee != VNULL);
   VASSERT(thee->initFlag);

   return (thee->nqsm)[iatom];

}

VPUBLIC SS* Vcsm_getSimplex(Vcsm *thee, int isimp, int iatom) {


   VASSERT(thee != VNULL);
   VASSERT(thee->initFlag);

   return Gem_SS(thee->gm, (thee->qsm)[iatom][isimp]);

}

VPUBLIC int Vcsm_getSimplexIndex(Vcsm *thee, int isimp, int iatom) {


   VASSERT(thee != VNULL);
   VASSERT(thee->initFlag);

   return (thee->qsm)[iatom][isimp];

}

VPUBLIC unsigned long int Vcsm_memChk(Vcsm *thee) {
    if (thee == VNULL) return 0;
    return Vmem_bytes(thee->vmem);
}

#endif /* if !defined(VINLINE_VCSM) */

VPUBLIC Vcsm* Vcsm_ctor(Valist *alist, Gem *gm) {

    /* Set up the structure */
    Vcsm *thee = VNULL;
    thee = Vmem_malloc(VNULL, 1, sizeof(Vcsm) );
    VASSERT( thee != VNULL);
    VASSERT( Vcsm_ctor2(thee, alist, gm));

    return thee;
}

VPUBLIC int Vcsm_ctor2(Vcsm *thee, Valist *alist, Gem *gm) { 
 
    VASSERT( thee != VNULL );

    /* Memory management object */
    thee->vmem = Vmem_ctor("APBS:VCSM");

    /* Set up the atom list and grid manager */
    if( alist == VNULL) {
        Vnm_print(2,"Vcsm_ctor2: got null pointer to Valist object!\n");
        return 0;
    }
    thee->alist = alist;
    if( gm == VNULL) {
        Vnm_print(2,"Vcsm_ctor2: got a null pointer to the Gem object!\n");
        return 0;
    }
    thee->gm = gm;
   
    thee->initFlag = 0;
    return 1;
}

VPUBLIC void Vcsm_init(Vcsm *thee) {
 
    /* Counters */
    int iatom, jatom, isimp, jsimp, gotSimp;
    /* Atomic information */
    Vatom *atom;
    double *position;
    /* Simplex/Vertex information */
    SS *simplex;
    /* Basis function values */

    if (thee == VNULL) {
        Vnm_print(2, "Vcsm_init:  Error!  Got NULL thee!\n");
        VASSERT(0);
    }
    if (thee->gm == VNULL) {
    VASSERT(thee->gm != VNULL);
        Vnm_print(2, "Vcsm_init:  Error!  Got NULL thee->gm!\n");
        VASSERT(0);
    }
    thee->nsimp = Gem_numSS(thee->gm);
    if (thee->nsimp <= 0) {
        Vnm_print(2, "Vcsm_init:  Error!  Got %d simplices!\n", thee->nsimp);
        VASSERT(0);
    }
    thee->natom = Valist_getNumberAtoms(thee->alist);

    /* Allocate and initialize space for the first dimensions of the 
     * simplex-charge map, the simplex array, and the counters */
    thee->sqm = Vmem_malloc(thee->vmem, thee->nsimp, sizeof(int *));
    VASSERT(thee->sqm != VNULL);
    thee->nsqm = Vmem_malloc(thee->vmem, thee->nsimp, sizeof(int));
    VASSERT(thee->nsqm != VNULL);
    for (isimp=0; isimp<thee->nsimp; isimp++) (thee->nsqm)[isimp] = 0;

    /* Count the number of charges per simplex. */
    for (iatom=0; iatom<thee->natom; iatom++) {
        atom = Valist_getAtom(thee->alist, iatom);
        position = Vatom_getPosition(atom);
        gotSimp = 0;
        for (isimp=0; isimp<thee->nsimp; isimp++) {
            simplex = Gem_SS(thee->gm, isimp);
            if (Gem_pointInSimplex(thee->gm, simplex, position)) {
                (thee->nsqm)[isimp]++;
                gotSimp = 1;
             }
        }
    }

    /* Allocate the space for the simplex-charge map */
    for (isimp=0; isimp<thee->nsimp; isimp++) {
        if ((thee->nsqm)[isimp] > 0) {
            thee->sqm[isimp] = Vmem_malloc(thee->vmem, (thee->nsqm)[isimp], 
              sizeof(int));
            VASSERT(thee->sqm[isimp] != VNULL);
        }
    }

    /* Finally, set up the map */
    for (isimp=0; isimp<thee->nsimp; isimp++) {
        jsimp = 0;
        simplex = Gem_SS(thee->gm, isimp);
        for (iatom=0; iatom<thee->natom; iatom++) {
            atom = Valist_getAtom(thee->alist, iatom);
            position = Vatom_getPosition(atom);
            /* Check to see if the atom's in this simplex */
            if (Gem_pointInSimplex(thee->gm, simplex, position)) {
                /* Assign the entries in the next vacant spot */
                (thee->sqm)[isimp][jsimp] = iatom;
                jsimp++;
            }
        }
    }

    thee->msimp = thee->nsimp;

    /* Allocate space for the charge-simplex map */
    thee->qsm = Vmem_malloc(thee->vmem, thee->natom, sizeof(int *));
    VASSERT(thee->qsm != VNULL);
    thee->nqsm = Vmem_malloc(thee->vmem, thee->natom, sizeof(int));
    VASSERT(thee->nqsm != VNULL);
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
            Vnm_print(2, "Vcsm_init: Atom %d not placed in simplex!\n", iatom);
            VASSERT(0);
        }
    }
    /* Allocate the appropriate amount of space for each entry in the
     * charge-simplex map and clear the counter for re-use in assignment */
    for (iatom=0; iatom<thee->natom; iatom++) {
        thee->qsm[iatom] = Vmem_malloc(thee->vmem, (thee->nqsm)[iatom],
          sizeof(int));
        VASSERT(thee->qsm[iatom] != VNULL);
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

VPUBLIC void Vcsm_dtor(Vcsm **thee) {
    if ((*thee) != VNULL) {
        Vcsm_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(Vcsm), (void **)thee);
        (*thee) = VNULL;
    }
}

VPUBLIC void Vcsm_dtor2(Vcsm *thee) { 
    int i;

    if ((thee != VNULL) && thee->initFlag) {

        for (i=0; i<thee->msimp; i++) {
            if (thee->nsqm[i] > 0) Vmem_free(thee->vmem, thee->nsqm[i], 
              sizeof(int), (void **)&(thee->sqm[i]));
        }
        for (i=0; i<thee->natom; i++) {
            if (thee->nqsm[i] > 0) Vmem_free(thee->vmem, thee->nqsm[i], 
              sizeof(int), (void **)&(thee->qsm[i]));
        }
        Vmem_free(thee->vmem, thee->msimp, sizeof(int *), 
          (void **)&(thee->sqm));
        Vmem_free(thee->vmem, thee->msimp, sizeof(int),
          (void **)&(thee->nsqm));
        Vmem_free(thee->vmem, thee->natom, sizeof(int *), 
          (void **)&(thee->qsm));
        Vmem_free(thee->vmem, thee->natom, sizeof(int),
          (void **)&(thee->nqsm));

    }
    Vmem_dtor(&(thee->vmem));
}

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
            thee->nsqm = Vmem_realloc(thee->vmem, thee->msimp, sizeof(int), 
              (void **)&(thee->nsqm), isimp);
            VASSERT(thee->nsqm != VNULL);
            thee->sqm = Vmem_realloc(thee->vmem, thee->msimp, sizeof(int *), 
              (void **)&(thee->sqm), isimp);
            VASSERT(thee->sqm != VNULL);
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

    sqmNew = Vmem_malloc(thee->vmem, num, sizeof(int *));
    VASSERT(sqmNew != VNULL);
    nsqmNew = Vmem_malloc(thee->vmem, num, sizeof(int));
    VASSERT(nsqmNew != VNULL);
    for (isimp=0; isimp<num; isimp++) nsqmNew[isimp] = 0;

    /* Loop throught the affected atoms to determine how many atoms each
     * simplex will get. */
    for (iatom=0; iatom<nqParent; iatom++) {

        atomID = qParent[iatom];
        atom = Valist_getAtom(thee->alist, atomID);
        position = Vatom_getPosition(atom);
        nsimps = 0;

        jsimp = 0;

        for (isimp=0; isimp<num; isimp++) {
            simplex = simps[isimp];
            if (Gem_pointInSimplex(thee->gm, simplex, position)) {
                nsqmNew[isimp]++;
                jsimp = 1;
            }
        }
 
        VASSERT(jsimp != 0);
    }

    /* Sanity check that we didn't lose any atoms... */
    iatom = 0;
    for (isimp=0; isimp<num; isimp++) iatom += nsqmNew[isimp];
    if (iatom < nqParent) {
        Vnm_print(2,"Vcsm_update: Lost %d (of %d) atoms!\n", 
            nqParent - iatom, nqParent);
        VASSERT(0);
    }

    /* Allocate the storage */
    for (isimp=0; isimp<num; isimp++) {
        if (nsqmNew[isimp] > 0) {
            sqmNew[isimp] = Vmem_malloc(thee->vmem, nsqmNew[isimp], 
              sizeof(int));
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
            if (Gem_pointInSimplex(thee->gm, simplex, position)) {
                sqmNew[isimp][jsimp] = atomID;
                jsimp++;
            }
        }
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
    dnqsm = Vmem_malloc(thee->vmem, nAffAtoms, sizeof(int));
    VASSERT(dnqsm != VNULL);
    nqsmNew = Vmem_malloc(thee->vmem, nAffAtoms, sizeof(int));
    VASSERT(nqsmNew != VNULL);
    qsmNew = Vmem_malloc(thee->vmem, nAffAtoms, sizeof(int*));
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
        qsmNew[iatom] = Vmem_malloc(thee->vmem, 
                (dnqsm[iatom] + thee->nqsm[atomID]), 
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
        Vmem_free(thee->vmem, thee->nqsm[atomID], sizeof(int), 
          (void **)&(thee->qsm[atomID]));
        thee->qsm[atomID] = qsmNew[iatom];
        thee->nqsm[atomID] = nqsmNew[iatom];
    }
    for (isimp=0; isimp<num; isimp++) {
        simpID = SS_id(simps[isimp]);
        if (thee->nsqm[simpID] > 0) Vmem_free(thee->vmem, thee->nsqm[simpID], 
          sizeof(int), (void **)&(thee->sqm[simpID]));
        thee->sqm[simpID] = sqmNew[isimp];
        thee->nsqm[simpID] = nsqmNew[isimp];
    }

    Vmem_free(thee->vmem, num, sizeof(int *), (void **)&sqmNew);
    Vmem_free(thee->vmem, num, sizeof(int), (void **)&nsqmNew);
    Vmem_free(thee->vmem, nAffAtoms, sizeof(int *), (void **)&qsmNew);
    Vmem_free(thee->vmem, nAffAtoms, sizeof(int), (void **)&nqsmNew);
    Vmem_free(thee->vmem, nAffAtoms, sizeof(int), (void **)&dnqsm);


    return 1;


}

#endif
