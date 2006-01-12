/**
 *  @file    valist.c
 *  @author  Nathan Baker
 *  @brief   Class Valist methods
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
 * Copyright (c) 2002-2005.  Washington University in St. Louis.
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
#include "apbs/valist.h"

VEMBED(rcsid="$Id$")

VPRIVATE char *Valist_whiteChars = " \t\n";
VPRIVATE char *Valist_commChars  = "#%";

/** Get statistics on a newly read-in atom list  */
VPRIVATE int Valist_getStatistics(Valist *thee);

#if !defined(VINLINE_VATOM)

VPUBLIC double Valist_getCenterX(Valist *thee) {
 
    VASSERT(thee != NULL);
    return thee->center[0];

}

VPUBLIC double Valist_getCenterY(Valist *thee) {

    VASSERT(thee != NULL);
    return thee->center[1];

}
VPUBLIC double Valist_getCenterZ(Valist *thee) {

    VASSERT(thee != NULL);
    return thee->center[2];

}

VPUBLIC Vatom* Valist_getAtomList(Valist *thee) {

    VASSERT(thee != NULL);
    return thee->atoms;

}

VPUBLIC int Valist_getNumberAtoms(Valist *thee) {

    VASSERT(thee != NULL);
    return thee->number;

}

VPUBLIC Vatom* Valist_getAtom(Valist *thee, int i) {

    VASSERT(thee != NULL);
    VASSERT(i < thee->number);
    return &(thee->atoms[i]);

}

VPUBLIC unsigned long int Valist_memChk(Valist *thee) {

    if (thee == NULL) return 0;
    return Vmem_bytes(thee->vmem);

}

#endif /* if !defined(VINLINE_VATOM) */

VPUBLIC Valist* Valist_ctor() {

    /* Set up the structure */
    Valist *thee = VNULL;
    thee = Vmem_malloc(VNULL, 1, sizeof(Valist));
    VASSERT( thee != VNULL);
    VASSERT( Valist_ctor2(thee));
 
    return thee;
}

VPUBLIC int Valist_ctor2(Valist *thee) {
  
    thee->atoms = VNULL;
    thee->number = 0;

    /* Initialize the memory management object */
    thee->vmem = Vmem_ctor("APBS:VALIST");

    return 1;    

}

VPUBLIC void Valist_dtor(Valist **thee)
{
    if ((*thee) != VNULL) {
        Valist_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(Valist), (void **)thee);
        (*thee) = VNULL;
    }
}

VPUBLIC void Valist_dtor2(Valist *thee) {

    Vmem_free(thee->vmem, thee->number, sizeof(Vatom), (void **)&(thee->atoms));
    thee->atoms = VNULL;
    thee->number = 0;

    Vmem_dtor(&(thee->vmem));
} 

/** Read serial number from PDB ATOM/HETATM field */
VPRIVATE int Valist_readPDBSerial(Valist *thee, Vio *sock, int *serial) {

    char tok[VMAX_BUFSIZE];
    int ti = 0;

    if (Vio_scanf(sock, "%s", tok) != 1) {
        Vnm_print(2, "Valist_readPDB:  Ran out of tokens while parsing serial!\n");
        return 0;
    } 
    if (sscanf(tok, "%d", &ti) != 1) {
        Vnm_print(2, "Valist_readPDB:  Unable to parse serial token (%s) as int!\n",
                tok);
        return 0;
    } 
    *serial = ti;

    return 1;
}

/** Read atom name from PDB ATOM/HETATM field */
VPRIVATE int Valist_readPDBAtomName(Valist *thee, Vio *sock, 
        char atomName[VMAX_ARGLEN]) {

    char tok[VMAX_BUFSIZE];

    if (Vio_scanf(sock, "%s", tok) != 1) {
        Vnm_print(2, "Valist_readPDB:  Ran out of tokens while parsing atom name!\n");
        return 0;
    }
    if (strlen(tok) < VMAX_ARGLEN) strcpy(atomName, tok);
    else {
        Vnm_print(2, "Valist_readPDB:  Atom name (%s) too long!\n", tok);
        return 0;
    }
    return 1;
}

/** Read residue name from PDB ATOM/HETATM field */
VPRIVATE int Valist_readPDBResidueName(Valist *thee, Vio *sock, 
        char resName[VMAX_ARGLEN]) {

    char tok[VMAX_BUFSIZE];

    if (Vio_scanf(sock, "%s", tok) != 1) {
        Vnm_print(2, "Valist_readPDB:  Ran out of tokens while parsing residue name!\n");
        return 0;
    }
    if (strlen(tok) < VMAX_ARGLEN) strcpy(resName, tok);
    else {
        Vnm_print(2, "Valist_readPDB:  Residue name (%s) too long!\n", tok);
        return 0;
    }
    return 1;
}

/** Read residue number from PDB ATOM/HETATM field */
VPRIVATE int Valist_readPDBResidueNumber(
        Valist *thee, Vio *sock, int *resSeq) {

    char tok[VMAX_BUFSIZE];
    int ti = 0;

    if (Vio_scanf(sock, "%s", tok) != 1) {
        Vnm_print(2, "Valist_readPDB:  Ran out of tokens while parsing resSeq!\n");
        return 0;
    } 
    if (sscanf(tok, "%d", &ti) != 1) {
        Vnm_print(2, "Valist_readPDB:  Unable to parse resSeq token (%s) as int!\n",
                tok);
        return 0;
    } 
    *resSeq = ti;

    return 1;
}

/** Read atom coordinate from PDB ATOM/HETATM field */
VPRIVATE int Valist_readPDBAtomCoord(Valist *thee, Vio *sock, double *coord) {

    char tok[VMAX_BUFSIZE];
    double tf = 0;

    if (Vio_scanf(sock, "%s", tok) != 1) {
        Vnm_print(2, "Valist_readPDB:  Ran out of tokens while parsing atom coordinate!\n");
        return 0;
    } 
    if (sscanf(tok, "%lf", &tf) != 1) {
        return 0;
    } 
    *coord = tf;

    return 1;
}

/** Read charge and radius from PQR ATOM/HETATM field */
VPRIVATE int Valist_readPDBChargeRadius(Valist *thee, Vio *sock, 
        double *charge, double *radius) {

    char tok[VMAX_BUFSIZE];
    double tf = 0;

    if (Vio_scanf(sock, "%s", tok) != 1) {
        Vnm_print(2, "Valist_readPQR:  Ran out of tokens while parsing charge!\n");
        return 0;
    } 
    if (sscanf(tok, "%lf", &tf) != 1) {
        return 0;
    } 
    *charge = tf;

    if (Vio_scanf(sock, "%s", tok) != 1) {
        Vnm_print(2, "Valist_readPQR:  Ran out of tokens while parsing radius!\n");
        return 0;
    } 
    if (sscanf(tok, "%lf", &tf) != 1) {
        return 0;
    } 
    *radius = tf;

    return 1;
}

/** Read ATOM/HETATM field of PDB through the X/Y/Z fields */
VPRIVATE int Valist_readPDB_throughXYZ(
        Valist *thee, 
        Vio *sock, /** Socket ready for reading */
        int *serial, /** Set to atom number */
        char atomName[VMAX_ARGLEN], /** Set to atom name */
        char resName[VMAX_ARGLEN], /** Set to residue name */
        int *resSeq, /** Set to residue number */
        double *x, /** Set to x-coordinate */
        double *y, /** Set to y-coordinate */
        double *z  /** Set to z-coordinate */
        ) {


    int i, njunk, gotit;

    /* Grab serial */
    if (!Valist_readPDBSerial(thee, sock, serial)) {
        Vnm_print(2, "Valist_readPDB:  Error while parsing serial!\n");
    }

    /* Grab atom name */
    if (!Valist_readPDBAtomName(thee, sock, atomName)) {
        Vnm_print(2, "Valist_readPDB:  Error while parsing atom name!\n");
        return 0;
    }

    /* Grab residue name */
    if (!Valist_readPDBResidueName(thee, sock, resName)) {
        Vnm_print(2, "Valist_readPDB:  Error while parsing residue name!\n");
        return 0;
    }


    /* Grab residue number */
    if (!Valist_readPDBResidueNumber(thee, sock, resSeq)) {
        Vnm_print(2, "Valist_readPDB:  Error while parsing residue name!\n");
        return 0;
    }


    /* Read tokens until we find one that can be parsed as an atom
     * x-coordinate.  We will allow njunk=1 intervening field that
     * cannot be parsed as a coordinate */
    njunk = 1;
    gotit = 0;
    for (i=0; i<(njunk+1); i++) {
        if (Valist_readPDBAtomCoord(thee, sock, x)) {
            gotit = 1;
            break;
        }
    }
    if (!gotit) {
        Vnm_print(2, "Valist_readPDB:  Can't find x!\n");
        return 0;
    }
    /* Read y-coordinate */
    if (!Valist_readPDBAtomCoord(thee, sock, y)) {
        Vnm_print(2, "Valist_readPDB:  Can't find y!\n");
        return 0;
    }
    /* Read z-coordinate */
    if (!Valist_readPDBAtomCoord(thee, sock, z)) {
        Vnm_print(2, "Valist_readPDB:  Can't find z!\n");
        return 0;
    }

#if 0 /* Set to 1 if you want to debug */
    Vnm_print(1, "Valist_readPDB:  serial = %d\n", serial);
    Vnm_print(1, "Valist_readPDB:  atomName = %s\n", atomName);
    Vnm_print(1, "Valist_readPDB:  resName = %s\n", resName);
    Vnm_print(1, "Valist_readPDB:  resSeq = %d\n", resSeq);
    Vnm_print(1, "Valist_readPDB:  pos = (%g, %g, %g)\n", 
            pos[0], pos[1], pos[2]);
#endif

    return 1;
}

/** Get a the next available atom storage location, increasing the storage
 * space if necessary.  Return VNULL if something goes wrong. */
VPRIVATE Vatom* Valist_getAtomStorage(
        Valist *thee,
        Vatom **plist, /** Pointer to existing list of atoms */
        int *pnlist, /** Size of existing list, may be changed */
        int *pnatoms /** Existing number of atoms in list; incremented 
                       before exit */
        ) {

    Vatom *oldList, *newList, *theList;
    Vatom *oldAtom, *newAtom;
    int iatom, inext, oldLength, newLength, natoms;

    newList = VNULL;

    /* See if we need more space */
    if (*pnatoms >= *pnlist) {

        /* Double the storage space */
        oldLength = *pnlist;
        newLength = 2*oldLength;
        newList = Vmem_malloc(thee->vmem, newLength, sizeof(Vatom));
        oldList = *plist;

        /* Check the allocation */
        if (newList == VNULL) {
            Vnm_print(2, "Valist_readPDB:  failed to allocate space for %d (Vatom)s!\n", newLength);
            return VNULL;
        }

        /* Copy the atoms over */
        natoms = *pnatoms;
        for (iatom=0; iatom<natoms; iatom++) { 
            oldAtom = &(oldList[iatom]);
            newAtom = &(newList[iatom]);
            Vatom_copyTo(oldAtom, newAtom);
            Vatom_dtor2(oldAtom);
        }

        /* Free the old list */
        Vmem_free(thee->vmem, oldLength, sizeof(Vatom), (void **)plist);

        /* Copy new list to plist */
        *plist = newList;
        *pnlist = newLength;
    }

    theList = *plist;
    inext = *pnatoms;

    /* Get the next available spot and increment counters */
    newAtom = &(theList[inext]);
    *pnatoms = inext + 1;

    return newAtom;
}

VPRIVATE int Valist_setAtomArray(Valist *thee, 
        Vatom **plist, /** Pointer to list of atoms to store */
        int nlist, /** Length of list */
        int natoms /** Number of real atom entries in list */
        ) {

    Vatom *list, *newAtom, *oldAtom;
    int i;

    list = *plist;

    /* Allocate necessary space */
    thee->number = 0;
    thee->atoms = Vmem_malloc(thee->vmem, natoms, sizeof(Vatom));
    if (thee->atoms == VNULL) {
        Vnm_print(2, "Valist_readPDB:  Unable to allocate space for %d (Vatom)s!\n", 
                natoms);
        return 0;
    }
    thee->number = natoms;

    /* Copy over data */
    for (i=0; i<thee->number; i++) {
        newAtom = &(thee->atoms[i]);
        oldAtom = &(list[i]);
        Vatom_copyTo(oldAtom, newAtom);
        Vatom_dtor2(oldAtom);
    }

    /* Free old array */
    Vmem_free(thee->vmem, nlist, sizeof(Vatom), (void **)plist);

    return 1;
}

VPUBLIC int Valist_readPDB(Valist *thee, Vparam *param, Vio *sock) {

    /* WE DO NOT DIRECTLY CONFORM TO PDB STANDARDS -- TO ALLOW LARGER FILES, WE
     * REQUIRE ALL FIELDS TO BE WHITESPACE DELIMITED */


    Vatom *atoms = VNULL;
    Vatom *nextAtom = VNULL;
    Vparam_AtomData *atomData = VNULL;
    char tok[VMAX_BUFSIZE];
    char atomName[VMAX_ARGLEN], resName[VMAX_ARGLEN]; 
    int nlist, natoms, serial, resSeq;
    double x, y, z, charge, radius;
    double pos[3];
 
    VASSERT(thee != VNULL);
    thee->number = 0;

    Vio_setWhiteChars(sock, Valist_whiteChars);
    Vio_setCommChars(sock, Valist_commChars);

    /* Allocate some initial space for the atoms */
    nlist = 200;
    atoms = Vmem_malloc(thee->vmem, nlist, sizeof(Vatom));

    natoms = 0;
    /* Read until we run out of lines */
    while (Vio_scanf(sock, "%s", tok) == 1) {

        /* Parse only ATOM/HETATOM fields */
        if ((Vstring_strcasecmp(tok, "ATOM") == 0) || 
            (Vstring_strcasecmp(tok, "HETATM") == 0)) {

            /* Read ATOM/HETATM field of PDB through the X/Y/Z fields */
            if (!Valist_readPDB_throughXYZ(thee, sock, &serial, atomName, 
                        resName, &resSeq, &x, &y, &z)) {
                Vnm_print(2, "Valist_readPDB:  Error parsing ATOM field!\n");
                return 0;
            }

            /* Try to find the parameters. */
            atomData = Vparam_getAtomData(param, resName, atomName);
            if (atomData == VNULL) {
                Vnm_print(2, "Valist_readPDB:  Couldn't find parameters for \
atom = %s, residue = %s\n", atomName, resName);
                return 0;
            }
            charge = atomData->charge;
            radius = atomData->radius;

            /* Get pointer to next available atom position */
            nextAtom = Valist_getAtomStorage(thee, &atoms, &nlist, &natoms);
            if (nextAtom == VNULL) {
                Vnm_print(2, "Valist_readPDB:  Error in allocating spacing for atoms!\n");
                return 0;
            }

            /* Store the information */
            pos[0] = x; pos[1] = y; pos[2] = z; 
            Vatom_setPosition(nextAtom, pos);
            Vatom_setCharge(nextAtom, charge);
            Vatom_setRadius(nextAtom, radius);
            Vatom_setAtomID(nextAtom, natoms);

        } /* if ATOM or HETATM */
    } /* while we haven't run out of tokens */

    Vnm_print(0, "Valist_readPDB: Counted %d atoms\n", natoms);
    fflush(stdout);

    /* Store atoms internally */
    if (!Valist_setAtomArray(thee, &atoms, nlist, natoms)) {
        Vnm_print(2, "Valist_readPDB:  unable to store atoms!\n");
        return 0;
    }

    return Valist_getStatistics(thee);


}

VPUBLIC int Valist_readPQR(Valist *thee, Vio *sock) {

    /* WE DO NOT DIRECTLY CONFORM TO PDB STANDARDS -- TO ALLOW LARGER FILES, WE
     * REQUIRE ALL FIELDS TO BE WHITESPACE DELIMITED */


    Vatom *atoms = VNULL;
    Vatom *nextAtom = VNULL;
    char tok[VMAX_BUFSIZE];
    char atomName[VMAX_ARGLEN], resName[VMAX_ARGLEN]; 
    int nlist, natoms, serial, resSeq;
    double x, y, z, charge, radius;
    double pos[3];
 
    VASSERT(thee != VNULL);
    thee->number = 0;

    Vio_setWhiteChars(sock, Valist_whiteChars);
    Vio_setCommChars(sock, Valist_commChars);

    /* Allocate some initial space for the atoms */
    nlist = 200;
    atoms = Vmem_malloc(thee->vmem, nlist, sizeof(Vatom));

    natoms = 0;
    /* Read until we run out of lines */
    while (Vio_scanf(sock, "%s", tok) == 1) {

        /* Parse only ATOM/HETATOM fields */
        if ((Vstring_strcasecmp(tok, "ATOM") == 0) || 
            (Vstring_strcasecmp(tok, "HETATM") == 0)) {

            /* Read ATOM/HETATM field of PDB through the X/Y/Z fields */
            if (!Valist_readPDB_throughXYZ(thee, sock, &serial, atomName, 
                        resName, &resSeq, &x, &y, &z)) {
                Vnm_print(2, "Valist_readPQR:  Error parsing ATOM field!\n");
                return 0;
            }

            /* Read Q/R fields */
            if (!Valist_readPDBChargeRadius(thee, sock, &charge, &radius)) {
                Vnm_print(2, "Valist_readPQR:  Error parsing ATOM field!\n");
                return 0;
            }


            /* Get pointer to next available atom position */
            nextAtom = Valist_getAtomStorage(thee, &atoms, &nlist, &natoms);
            if (nextAtom == VNULL) {
                Vnm_print(2, "Valist_readPQR:  Error in allocating spacing for atoms!\n");
                return 0;
            }

            /* Store the information */
            pos[0] = x; pos[1] = y; pos[2] = z; 
            Vatom_setPosition(nextAtom, pos);
            Vatom_setCharge(nextAtom, charge);
            Vatom_setRadius(nextAtom, radius);
            Vatom_setAtomID(nextAtom, natoms-1);

        } /* if ATOM or HETATM */
    } /* while we haven't run out of tokens */

    Vnm_print(0, "Valist_readPQR: Counted %d atoms\n", natoms);
    fflush(stdout);

    /* Store atoms internally */
    if (!Valist_setAtomArray(thee, &atoms, nlist, natoms)) {
        Vnm_print(2, "Valist_readPDB:  unable to store atoms!\n");
        return 0;
    }

    return Valist_getStatistics(thee);


}

/** Load up Valist with various statistics */
VPRIVATE int Valist_getStatistics(Valist *thee) {

    Vatom *atom;
    int i, j;

    VASSERT(thee != VNULL);

    thee->center[0] = 0.;
    thee->center[1] = 0.;
    thee->center[2] = 0.;
    thee->maxrad = 0.;
    thee->charge = 0.;

    if (thee->number == 0) return 0;

    /* Reset stat variables */
    atom = &(thee->atoms[0]);
    for (i=0; i<3; i++) {
        thee->maxcrd[i] = thee->mincrd[i] = atom->position[i];
    }
    thee->maxrad = atom->radius;
    thee->charge = 0.0;

    for (i=0; i<thee->number; i++) {

        atom = &(thee->atoms[i]);
        for (j=0; j<3; j++) {
            if (atom->position[j] < thee->mincrd[j]) 
              thee->mincrd[j] = atom->position[j];
            if (atom->position[j] > thee->maxcrd[j]) 
              thee->maxcrd[j] = atom->position[j];
        }
        if (atom->radius > thee->maxrad) thee->maxrad = atom->radius;
        thee->charge = thee->charge + atom->charge;
    } 

    thee->center[0] = 0.5*(thee->maxcrd[0] + thee->mincrd[0]);
    thee->center[1] = 0.5*(thee->maxcrd[1] + thee->mincrd[1]);
    thee->center[2] = 0.5*(thee->maxcrd[2] + thee->mincrd[2]);

    return 1;
}
