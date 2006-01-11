/**
 *  @file    vparam.c
 *  @ingroup Vparam
 *  @author  Nathan Baker
 *  @brief   Class Vparam methods
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
 * programs or libraries that are released under the GNU LGPL and with
 * code included in releases of ISIM, PMV, PyMOL, SMOL, VMD.  This
 * special exception permission is also extended to any software listed
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
#include "apbs/vparam.h"

#if defined(HAVE_MC_H)
#include "mc/mc.h"
#endif

VEMBED(rcsid="$Id$")

/**
 * @brief  Whitespace characters for socket reads 
 * @ingroup  Vparam
 */
VPRIVATE char *MCwhiteChars = " =,;\t\n\r";

/**
 * @brief  Comment characters for socket reads 
 * @ingroup  Vparam
 */
VPRIVATE char *MCcommChars  = "#%";

/**
 * @brief  Read a single line of the flat file database
 * @author  Nathan Baker
 * @ingroup  Vparam
 * @param  sock  Socket ready for reading
 * @param  atom  Atom to hold parsed data
 * @returns 1 if successful, 0 otherwise
 */
VPRIVATE int readFlatFileLine(Vio *sock, Vparam_AtomData *atom);

#if !defined(VINLINE_VPARAM)

VPUBLIC unsigned long int Vparam_memChk(Vparam *thee) {
    if (thee == VNULL) return 0;
    return Vmem_bytes(thee->vmem);
}

#endif /* if !defined(VINLINE_VPARAM) */

VPUBLIC Vparam_AtomData* Vparam_AtomData_ctor() {

    Vparam_AtomData *thee = VNULL;

    /* Set up the structure */
    thee = Vmem_malloc(VNULL, 1, sizeof(Vparam_AtomData) );
    VASSERT(thee != VNULL);
    VASSERT(Vparam_AtomData_ctor2(thee));

    return thee;
}

VPUBLIC int Vparam_AtomData_ctor2(Vparam_AtomData *thee) { return 1; }

VPUBLIC void Vparam_AtomData_dtor(Vparam_AtomData **thee) {
    
    if ((*thee) != VNULL) {
        Vparam_AtomData_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(Vparam_AtomData), (void **)thee);
        (*thee) = VNULL;
    }

}

VPUBLIC void Vparam_AtomData_dtor2(Vparam_AtomData *thee) { ; }

VPUBLIC Vparam_ResData* Vparam_ResData_ctor(Vmem *mem) {

    Vparam_ResData *thee = VNULL;

    /* Set up the structure */
    thee = Vmem_malloc(mem, 1, sizeof(Vparam_ResData) );
    VASSERT(thee != VNULL);
    VASSERT(Vparam_ResData_ctor2(thee, mem));

    return thee;
}

VPUBLIC int Vparam_ResData_ctor2(Vparam_ResData *thee, Vmem *mem) { 
    
    if (thee == VNULL) {
        Vnm_print(2, "Vparam_ResData_ctor2:  Got VNULL thee!\n");
        return 0;
    }
    thee->vmem = mem;
    thee->nAtomData = 0;
    thee->atomData = VNULL;

    return 1;
}

VPUBLIC void Vparam_ResData_dtor(Vparam_ResData **thee) {
    
    if ((*thee) != VNULL) {
        Vparam_ResData_dtor2(*thee);
        Vmem_free((*thee)->vmem, 1, sizeof(Vparam_ResData), (void **)thee);
        (*thee) = VNULL;
    }

}

VPUBLIC void Vparam_ResData_dtor2(Vparam_ResData *thee) { 
    
    if (thee == VNULL) return; 
    if (thee->nAtomData > 0) {
        Vmem_free(thee->vmem, thee->nAtomData, sizeof(Vparam_AtomData), 
          (void **)&(thee->atomData));
    }
    thee->nAtomData = 0;
    thee->atomData = VNULL;
}

VPUBLIC Vparam* Vparam_ctor() {

    Vparam *thee = VNULL;

    /* Set up the structure */
    thee = Vmem_malloc(VNULL, 1, sizeof(Vparam) );
    VASSERT(thee != VNULL);
    VASSERT(Vparam_ctor2(thee));

    return thee;
}

VPUBLIC int Vparam_ctor2(Vparam *thee) {

    if (thee == VNULL) {
        Vnm_print(2, "Vparam_ctor2: got VNULL thee!\n");
        return 0;
    }

    thee->vmem = VNULL;
    thee->vmem = Vmem_ctor("APBS:VPARAM");
    if (thee->vmem == VNULL) {
        Vnm_print(2, "Vparam_ctor2:  failed to init Vmem!\n");
        return 0;
    }

    thee->nResData = 0;
    thee->resData = VNULL;

    return 1;
}

VPUBLIC void Vparam_dtor(Vparam **thee) {
    
    if ((*thee) != VNULL) {
        Vparam_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(Vparam), (void **)thee);
        (*thee) = VNULL;
    }

}

VPUBLIC void Vparam_dtor2(Vparam *thee) {

    int i;

    if (thee == VNULL) return;

    /* Destroy the residue data */
    for (i=0; i<thee->nResData; i++) Vparam_ResData_dtor2(&(thee->resData[i]));
    if (thee->nResData > 0) Vmem_free(thee->vmem, thee->nResData, 
      sizeof(Vparam_ResData), (void **)&(thee->resData));
    thee->nResData = 0;
    thee->resData = VNULL;

    if (thee->vmem != VNULL) Vmem_dtor(&(thee->vmem));
    thee->vmem = VNULL;

}

VPUBLIC Vparam_ResData* Vparam_getResData(Vparam *thee, 
  char resName[VMAX_ARGLEN]) {

    int i;
    Vparam_ResData *res = VNULL;

    VASSERT(thee != VNULL);

    if ((thee->nResData == 0) || (thee->resData == VNULL)) {
        res = VNULL;
        return res;
    }

    /* Look for the matching residue */
    for (i=0; i<thee->nResData; i++) {
        res = &(thee->resData[i]);
        if (Vstring_strcasecmp(resName, res->name) == 0) return res;

    }

    /* Didn't find a matching residue */
    res = VNULL;
    Vnm_print(2, "Vparam_getResData:  unable to find res=%s\n", resName);
    return res;
}

VPUBLIC Vparam_AtomData* Vparam_getAtomData(Vparam *thee, 
  char resName[VMAX_ARGLEN], char atomName[VMAX_ARGLEN]) {

    int i;
    Vparam_ResData *res = VNULL;
    Vparam_AtomData *atom = VNULL;

    VASSERT(thee != VNULL);

    if ((thee->nResData == 0) || (thee->resData == VNULL)) {
        atom = VNULL;
        return atom;
    }

    /* Look for the matching residue */
    res = Vparam_getResData(thee, resName);
    if (res == VNULL) {
        atom = VNULL;
        return atom;
    }
    for (i=0; i<res->nAtomData; i++) {
        atom = &(res->atomData[i]);
        if (Vstring_strcasecmp(atomName, atom->atomName) == 0) {
            return atom;
        }
    }

    /* Didn't find a matching atom/residue */
    atom = VNULL;
    Vnm_print(2, "Vparam_getAtomData:  unable to find atom=%s, res=%s\n",
      atomName, resName);
    return atom;
}

VPUBLIC int Vparam_readFlatFile(Vparam *thee, const char *iodev,
  const char *iofmt, const char *thost, const char *fname) {

    int i, iatom, jatom, ires, natoms, nalloc;
    Vparam_AtomData *atoms = VNULL;
    Vparam_AtomData *tatoms = VNULL;
    Vparam_AtomData *atom = VNULL;
    Vparam_ResData *res = VNULL;
    Vio *sock = VNULL;
    char currResName[VMAX_ARGLEN];

    VASSERT(thee != VNULL);

    /* Setup communication */
    sock = Vio_ctor(iodev,iofmt,thost,fname,"r");
    if (sock == VNULL) {
        Vnm_print(2, "Vparam_readFlatFile: Problem opening virtual socket %s\n",
          fname);
        return 0;
    }
    if (Vio_accept(sock, 0) < 0) {
        Vnm_print(2, "Vparam_readFlatFile: Problem accepting virtual socket %s\n",
          fname);
        return 0;
    }
    Vio_setWhiteChars(sock, MCwhiteChars);
    Vio_setCommChars(sock, MCcommChars);

    /* Clear existing parameters */
    if (thee->nResData > 0) {
        Vnm_print(2, "WARNING -- CLEARING PARAMETER DATABASE!\n");
        for (i=0; i<thee->nResData; i++) {
            Vparam_ResData_dtor2(&(thee->resData[i]));
        }
        Vmem_free(thee->vmem, thee->nResData, 
          sizeof(Vparam_ResData), (void **)&(thee->resData));
    }

    /* Initial space for atoms */
    nalloc = 200;
    natoms = 0;
    atoms = Vmem_malloc(thee->vmem, nalloc, sizeof(Vparam_AtomData));

    /* Read until we run out of entries, allocating space as needed */
    while (1) {
        if (natoms >= nalloc) {
            tatoms = Vmem_malloc(thee->vmem, 2*nalloc, sizeof(Vparam_AtomData));
            VASSERT(tatoms != VNULL);
            for (i=0; i<natoms; i++) {
                Vparam_AtomData_copyTo(&(atoms[i]), &(tatoms[i]));
            }
            Vmem_free(thee->vmem, nalloc, sizeof(Vparam_AtomData), 
              (void **)&(atoms));
            atoms = tatoms; 
            tatoms = VNULL;
            nalloc = 2*nalloc;
        }
        atom = &(atoms[natoms]);
        if (!readFlatFileLine(sock, atom)) break;
        natoms++;
    }
    if (natoms == 0) return 0;

    /* Count the number of residues */
    thee->nResData = 1;
    strcpy(currResName, atoms[0].resName);
    for (i=1; i<natoms; i++) {
        if (Vstring_strcasecmp(atoms[i].resName, currResName) != 0) {
            strcpy(currResName, atoms[i].resName);
            (thee->nResData)++;
        }
    }

    /* Create the residues */
    thee->resData = Vmem_malloc(thee->vmem, thee->nResData, 
      sizeof(Vparam_ResData));
    VASSERT(thee->resData != VNULL);
    for (i=0; i<(thee->nResData); i++) {
        res = &(thee->resData[i]);
        Vparam_ResData_ctor2(res, thee->vmem);
    }

    /* Count the number of atoms per residue */
    ires = 0;
    res = &(thee->resData[ires]);
    res->nAtomData = 1;
    strcpy(res->name, atoms[0].resName);
    for (i=1; i<natoms; i++) {
        if (Vstring_strcasecmp(atoms[i].resName, res->name) != 0) {
            (ires)++;
            res = &(thee->resData[ires]);
            res->nAtomData = 1;
            strcpy(res->name, atoms[i].resName);
        } else (res->nAtomData)++;
    }

    /* Allocate per-residue space for atoms */
    for (ires=0; ires<thee->nResData; ires++) {
        res = &(thee->resData[ires]);
        res->atomData = Vmem_malloc(thee->vmem, res->nAtomData, 
          sizeof(Vparam_AtomData));
    }

    /* Copy atoms into residues */
    iatom = 0;
    Vparam_AtomData_copyTo(&(atoms[0]), &(res->atomData[iatom]));
    for (ires=0; ires<thee->nResData; ires++) {
        res = &(thee->resData[ires]);
        for (jatom=0; jatom<res->nAtomData; jatom++) {
            Vparam_AtomData_copyTo(&(atoms[iatom]), &(res->atomData[jatom]));
            iatom++;
        }
    }

    /* Shut down communication */
    Vio_acceptFree(sock);
    Vio_dtor(&sock);

    /* Destroy temporary atom space */
    Vmem_free(thee->vmem, nalloc, sizeof(Vparam_AtomData), (void **)&(atoms));

    return 1;

}

VEXTERNC void Vparam_AtomData_copyTo(Vparam_AtomData *thee,
  Vparam_AtomData *dest) {

    VASSERT(thee != VNULL);
    VASSERT(dest != VNULL);

    strcpy(dest->atomName, thee->atomName);
    strcpy(dest->resName, thee->resName);
    dest->charge = thee->charge;
    dest->radius = thee->radius;
    dest->epsilon = thee->epsilon;

}

VEXTERNC void Vparam_AtomData_copyFrom(Vparam_AtomData *thee,
  Vparam_AtomData *src) {  Vparam_AtomData_copyTo(src, thee); }

VPRIVATE int readFlatFileLine(Vio *sock, Vparam_AtomData *atom) {

    double dtmp;
    char tok[VMAX_BUFSIZE];

    VASSERT(atom != VNULL);

    if (Vio_scanf(sock, "%s", tok) != 1) return 0;
    if (strlen(tok) > VMAX_ARGLEN) {
        Vnm_print(2, "Vparam_readFlatFile:  string (%s) too long (%d)!\n", 
          tok, strlen(tok));
        return 0;
    }
    strcpy(atom->resName, tok);
    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (strlen(tok) > VMAX_ARGLEN) {
        Vnm_print(2, "Vparam_readFlatFile:  string (%s) too long (%d)!\n", 
          tok, strlen(tok));
        return 0;
    }
    strcpy(atom->atomName, tok);
    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%lf", &dtmp) != 1) {
        Vnm_print(2, "Vparam_readFlatFile:  Unexpected token (%s) while \
parsing charge!\n", tok);
        return 0;
    }
    atom->charge = dtmp;
    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%lf", &dtmp) != 1) {
        Vnm_print(2, "Vparam_readFlatFile:  Unexpected token (%s) while \
parsing radius!\n", tok);
        return 0;
    }
    atom->radius = dtmp;
    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%lf", &dtmp) != 1) {
        Vnm_print(2, "Vparam_readFlatFile:  Unexpected token (%s) while \
parsing radius!\n", tok);
        return 0;
    }
    atom->epsilon = dtmp;

    return 1;

VERROR1:
    Vnm_print(2, "Vparam_readFlatFile: Got unexpected EOF reading parameter file!\n");
    return 0;
}
