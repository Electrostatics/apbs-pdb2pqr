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
 * Copyright (c) 2003.  Washington University in St. Louis.
 * All Rights Reserved.
 * Portions Copyright (c) 1999-2003.  The Regents of the University of
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
 * @endverbatim
 */

#include "apbscfg.h"
#include "apbs/valist.h"

VEMBED(rcsid="$Id$")

VPRIVATE char *MCwhiteChars = " =,;\t\n";
VPRIVATE char *MCcommChars  = "#%";

/**
 * @brief  Get statistics on a newly read-in atom list 
 * @author  Nathan Baker
 * @ingroup  Valist
 * @param  thee  Valist object
 * @return  1 if successful, 0 otherwise
 */
VPRIVATE int getStatistics(Valist *thee);

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

VPUBLIC int Valist_memChk(Valist *thee) {

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

VPUBLIC int Valist_readPDB(Valist *thee, Vparam *param, const char *iodev, 
  const char *iofmt, const char *thost, const char *fname) {

    /* WE DO NOT DIRECTLY CONFORM TO PDB STANDARDS -- TO ALLOW LARGER FILES, WE
     * REQUIRE ALL FIELDS TO BE WHITESPACE DELIMITED */

    Vio *sock = VNULL;
    Vatom *atoms = VNULL;
    Vatom *tatoms = VNULL;
    Vparam_AtomData *atomData = VNULL;
    char tok[VMAX_BUFSIZE], tokArray[4][VMAX_BUFSIZE], stmp[VMAX_BUFSIZE];
    char atomName[VMAX_ARGLEN], resName[VMAX_ARGLEN]; 
    int nalloc, itmp, i, gotit, ntok;
    double dtmp, x, y, z, charge, radius;
    double pos[3];
 
    VASSERT(thee != VNULL);
    thee->number = 0;

    /* Open socket for reading */
    sock = Vio_ctor(iodev,iofmt,thost,fname,"r");
    if (sock == VNULL) {
        Vnm_print(2, "Valist_readPDB: Problem opening virtual socket %s\n",
          fname);
        return 0;
    }
    if (Vio_accept(sock, 0) < 0) {
        Vnm_print(2, "Valist_readPDB: Problem accepting virtual socket %s\n",
          fname);
        return 0;
    }
    Vio_setWhiteChars(sock, MCwhiteChars);
    Vio_setCommChars(sock, MCcommChars);

    /* Allocate some initial space for the atoms */
    nalloc = 200;
    atoms = Vmem_malloc(thee->vmem, nalloc, sizeof(Vatom));

    /* Read until we run out of lines */
    thee->number = 0;
    while (1) {

        if (Vio_scanf(sock, "%s", tok) != 1) break;
        if ((Vstring_strcasecmp(tok, "ATOM") == 0) || 
            (Vstring_strcasecmp(tok, "HETATM") == 0)) {

            /* Grab serial */
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%d", &itmp) != 1) {
                Vnm_print(2, "Valist_readPDB:  Error while parsing serial!\n");
                return 0;
            }

            /* Grab name */
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (strlen(tok) < VMAX_ARGLEN) strcpy(atomName, tok);
            else {
                Vnm_print(2, "Valist_readPDB:  Atom name (%s) too long!\n",
                  tok);
                return 0;
            }

            /* We don't care about any of the next 1-3 fields; the next thing
             * we're looking for is resSeq (integer) */
            gotit = 0;
            ntok = 0;
            for (i=0; i<4; i++) {
                VJMPERR1(Vio_scanf(sock, "%s", tokArray[i]) == 1);
                ntok++;
                if ((sscanf(tok, "%d", &itmp) == 1) && 
                        (sscanf(tok, "%s%d%s", stmp, &itmp, stmp) == 1) &&
                        (sscanf(tok, "%d%s", &itmp, stmp) == 1) &&
                        (sscanf(tok, "%s%d", &itmp) == 1)) {
                    gotit = 1;
                    break;
                }
            }
            if (!gotit) {
                Vnm_print(2, "Valist_readPDB:  Can't find resSeq!\n");
                return 0;
            }

            /* We don't care about any of the next 1 fields; the next thing
             * we're looking for is x (float) */
            gotit = 0;
            for (i=0; i<2; i++) {
                VASSERT(Vio_scanf(sock, "%s", tok) == 1);
                if (sscanf(tok, "%lf", &dtmp) == 1) {
                    gotit = 1;
                    break;
                }
            }
            if (!gotit) {
                Vnm_print(2, "Valist_readPDB:  Can't find x!\n");
                return 0;
            }
            pos[0] = dtmp;
            VASSERT(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &dtmp) != 1) {
                Vnm_print(2, "Valist_readPDB:  Can't find y!\n");
                return 0;
            }
            pos[1] = dtmp;
            VASSERT(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &dtmp) != 1) {
                Vnm_print(2, "Valist_readPDB:  Can't find z!\n");
                return 0;
            }
            pos[2] = dtmp;

            /* Try to find the parameters.  We have to loop through all the
             * entries in tokArray because we can't be sure which is the
             * actual resName */
            gotit = 0;
            for (i=0; i<ntok; i++) {
                if (strlen(tokArray[i]) < VMAX_ARGLEN) {
                    strcpy(resName, tokArray[i]);
                    atomData = Vparam_getAtomData(param, resName, atomName);
                    if (atomData != VNULL) {
                        gotit = 1;
                        charge = atomData->charge;
                        radius = atomData->radius;
                        break;
                    }
                }
            }
            if (!gotit) {
                Vnm_print(2, "Valist_readPDB:  Couldn't find parameters for \
atom=%s using following residue names as guesses: ");
                for (i=0; i<ntok; i++) {
                    if (strlen(tokArray[i]) < VMAX_ARGLEN) {
                        Vnm_print(2, " %s,", tokArray[i]);
                    }
                }
                Vnm_print(2, "\n");
                return 0;
            }

            /* Allocate more space for the new atom (if necessary) */
            if (thee->number == (nalloc-1)) {
                tatoms = Vmem_malloc(thee->vmem, 2*nalloc, sizeof(Vatom));
                VASSERT(tatoms != VNULL);
                for (i=0; i<thee->number; i++) {
                    Vatom_copyTo(&(atoms[i]), &(tatoms[i]));
                    Vatom_dtor2(&(atoms[i]));
                }
                Vmem_free(thee->vmem, nalloc, sizeof(Vatom), (void **)&atoms);
                nalloc = 2*nalloc;
                atoms = tatoms;
                tatoms = VNULL;
            }
            Vatom_setPosition(&(atoms[thee->number]), pos);
            Vatom_setCharge(&(atoms[thee->number]), charge);
            Vatom_setRadius(&(atoms[thee->number]), radius);
            (thee->number)++;

        } /* if ATOM or HETATM */
    } /* while we haven't run out of tokens */


    Vnm_print(0, "Valist_readPDB: Counted %d atoms\n", thee->number);
    fflush(stdout);


    /* Allocate the necessary space for the actual atom array */
    thee->atoms = Vmem_malloc(thee->vmem, thee->number, (sizeof(Vatom)));
    VASSERT(thee->atoms != VNULL);
    for (i=0; i<thee->number; i++) {
        Vatom_copyTo(&(atoms[i]), &(thee->atoms[i]));
        Vatom_dtor2(&(atoms[i]));
    }
    Vmem_free(thee->vmem, nalloc, sizeof(Vatom), (void **)&atoms);

    /* Close socket */
    Vio_acceptFree(sock);
    Vio_dtor(&sock);
      
    return getStatistics(thee);

VERROR1:
    Vnm_print(2, "Valist_readPDB:  Ran out of tokens!\n");
    return 0;

}

VPUBLIC int Valist_readPQR(Valist *thee, const char *iodev, const char *iofmt,
  const char *thost, const char *fname) {

    /* WE DO NOT DIRECTLY CONFORM TO PDB STANDARDS -- TO ALLOW LARGER FILES, WE
     * REQUIRE ALL FIELDS TO BE WHITESPACE DELIMITED */

    Vio *sock = VNULL;
    Vatom *atoms = VNULL;
    Vatom *tatoms = VNULL;
    char tok[VMAX_BUFSIZE], stmp[VMAX_BUFSIZE]; 
    int nalloc, itmp, i, gotit;
    double dtmp, x, y, z, charge, radius;
    double pos[3];
 
    VASSERT(thee != VNULL);
    thee->number = 0;

    /* Open socket for reading */
    sock = Vio_ctor(iodev,iofmt,thost,fname,"r");
    if (sock == VNULL) {
        Vnm_print(2, "Valist_readPQR: Problem opening virtual socket %s\n",
          fname);
        return 0;
    }
    if (Vio_accept(sock, 0) < 0) {
        Vnm_print(2, "Valist_readPQR: Problem accepting virtual socket %s\n",
          fname);
        return 0;
    }
    Vio_setWhiteChars(sock, MCwhiteChars);
    Vio_setCommChars(sock, MCcommChars);

    /* Allocate some initial space for the atoms */
    nalloc = 200;
    atoms = Vmem_malloc(thee->vmem, nalloc, sizeof(Vatom));

    /* Read until we run out of lines */
    thee->number = 0;
    while (1) {

        if (Vio_scanf(sock, "%s", tok) != 1) break;
        if ((Vstring_strcasecmp(tok, "ATOM") == 0) || 
            (Vstring_strcasecmp(tok, "HETATM") == 0)) {

            /* Grab serial */
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%d", &itmp) != 1) {
                Vnm_print(2, "Valist_readPQR:  Error while parsing serial!\n");
                return 0;
            }

            /* Grab name */
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);

            /* We don't care about any of the next 1-3 fields; the next thing
             * we're looking for is resSeq (integer) */
            gotit = 0;
            for (i=0; i<4; i++) {
                VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
                if ((sscanf(tok, "%d", &itmp) == 1) && 
                        (sscanf(tok, "%s%d%s", stmp, &itmp, stmp) == 1) &&
                        (sscanf(tok, "%d%s", &itmp, stmp) == 1) &&
                        (sscanf(tok, "%s%d", &itmp) == 1)) {
                    /* Vnm_print(1, "DEBUG:  parsed %s as integer.\n", tok); */
                    gotit = 1;
                    break;
                }
            }
            if (!gotit) {
                Vnm_print(2, "Valist_readPQR:  Can't find resSeq!\n");
                return 0;
            }

            /* We don't care about any of the next 1 fields; the next thing
             * we're looking for is x (float) */
            gotit = 0;
            for (i=0; i<2; i++) {
                VASSERT(Vio_scanf(sock, "%s", tok) == 1);
                if (sscanf(tok, "%lf", &dtmp) == 1) {
                    gotit = 1;
                    break;
                }
            }
            if (!gotit) {
                Vnm_print(2, "Valist_readPQR:  Can't find x!\n");
                return 0;
            }
            pos[0] = dtmp;
            VASSERT(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &dtmp) != 1) {
                Vnm_print(2, "Valist_readPQR:  Can't find y!\n");
                return 0;
            }
            pos[1] = dtmp;
            VASSERT(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &dtmp) != 1) {
                Vnm_print(2, "Valist_readPQR:  Can't find z!\n");
                return 0;
            }
            pos[2] = dtmp;
            VASSERT(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &dtmp) != 1) {
                Vnm_print(2, "Valist_readPQR:  Can't find charge!\n");
                return 0;
            }
            charge = dtmp;
            VASSERT(Vio_scanf(sock, "%s", tok) == 1);
            if (sscanf(tok, "%lf", &dtmp) != 1) {
                Vnm_print(2, "Valist_readPQR:  Can't find radius!\n");
                return 0;
            }
            radius = dtmp;
            if (radius < 0.0) {
                Vnm_print(2, "Valist_readPQR:  radii can't be negative (%g)!\n",
                        radius);
                return 0;
            }
            /* Vnm_print(1, "DEBUG:  x = %g, y = %g, z = %g, charge = %g, radius = %g\n", pos[0], pos[1], pos[2], charge, radius); */

            /* Allocate more space for the new atom (if necessary) */
            if (thee->number == (nalloc-1)) {
                tatoms = Vmem_malloc(thee->vmem, 2*nalloc, sizeof(Vatom));
                VASSERT(tatoms != VNULL);
                for (i=0; i<thee->number; i++) {
                    Vatom_copyTo(&(atoms[i]), &(tatoms[i]));
                    Vatom_dtor2(&(atoms[i]));
                }
                Vmem_free(thee->vmem, nalloc, sizeof(Vatom), (void **)&atoms);
                nalloc = 2*nalloc;
                atoms = tatoms;
                tatoms = VNULL;
            }
            Vatom_setCharge(&(atoms[thee->number]), charge);
            Vatom_setRadius(&(atoms[thee->number]), radius);
            Vatom_setPosition(&(atoms[thee->number]), pos);
            (thee->number)++;
        } /* if ATOM or HETATM */
    } /* while we haven't run out of tokens */


#if defined(VDEBUG)
    Vnm_print(1, "Valist_readPQR: Counted %d atoms\n",thee->number);
    fflush(stdout);
#endif

    /* Allocate the necessary space for the actual atom array */
    thee->atoms = Vmem_malloc(thee->vmem, thee->number,(sizeof(Vatom)));
    VASSERT(thee->atoms != VNULL);
    for (i=0; i<thee->number; i++) {
        Vatom_copyTo(&(atoms[i]), &(thee->atoms[i]));
        Vatom_dtor2(&(atoms[i]));
    }
    Vmem_free(thee->vmem, nalloc, sizeof(Vatom), (void **)&atoms);

    /* Close socket */
    Vio_acceptFree(sock);
    Vio_dtor(&sock);
      
    return getStatistics(thee);

VERROR1:
    Vnm_print(2, "Valist_readPQR:  Ran out of tokens!\n");
    return 0;

}

VPRIVATE int getStatistics(Valist *thee) {

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

VPUBLIC void Valist_buildMesh(Valist *thee, double size, const char *iodev,
  const char *iofmt, const char *thost, const char *fname) {

    Vatom *atom;
    Vio *sock;
    double centC[3], minC[3], maxC[3], pos, rad;
    double x0, x1, y0, y1, z0, z1;
    int i, j;

    VASSERT(thee != VNULL);

    /* Get the protein center and dimensions */
    atom = &((thee->atoms)[0]);
    for (j=0; j<3; j++) {
       rad = Vatom_getRadius(atom);
       pos = (Vatom_getPosition(atom))[j];
       minC[j] = pos - rad;
       maxC[j] = pos + rad;
       centC[j] = pos;
    }
    for (i=1; i<thee->number; i++) {
        for (j=0; j<3; j++) {
            pos = (Vatom_getPosition(atom))[j];
            centC[j] += pos;
            if ((pos+rad) > maxC[j]) maxC[j] = (pos + rad);
            if ((pos-rad) < minC[j]) minC[j] = (pos - rad);
        }
    }
    for (j=0; j<3; j++) centC[j] = centC[j]/((double)(thee->number));

    /* Determine the box corner positions */
    x0 = centC[0] - size*(maxC[0]-minC[0])/2;
    x1 = centC[0] + size*(maxC[0]-minC[0])/2;
    y0 = centC[1] - size*(maxC[1]-minC[1])/2;
    y1 = centC[1] + size*(maxC[1]-minC[1])/2;
    z0 = centC[2] - size*(maxC[2]-minC[2])/2;
    z1 = centC[2] + size*(maxC[2]-minC[2])/2;
   
    /* Open up a socket for writing out the mesh */ 
    sock = Vio_ctor(iodev, iofmt, thost, fname, "w");
    VASSERT(sock != VNULL);

    /* Write out the MCSF header */
    Vio_printf(sock, "mcsf_begin=1;\n");
    Vio_printf(sock, "dim=3;\n");
    Vio_printf(sock, "dimii=3;\n");
    Vio_printf(sock, "vertices=8;\n");
    Vio_printf(sock, "simplices=6;\n");

    /* Write out the vertices */
    Vio_printf(sock, "vert=[\n");
    Vio_printf(sock, "0 0  %11.10e %11.10e %11.10e\n", x0, y0, z0);
    Vio_printf(sock, "1 0  %11.10e %11.10e %11.10e\n", x1, y0, z0);
    Vio_printf(sock, "2 0  %11.10e %11.10e %11.10e\n", x0, y1, z0);
    Vio_printf(sock, "3 0  %11.10e %11.10e %11.10e\n", x1, y1, z0);
    Vio_printf(sock, "4 0  %11.10e %11.10e %11.10e\n", x0, y0, z1);
    Vio_printf(sock, "5 0  %11.10e %11.10e %11.10e\n", x1, y0, z1);
    Vio_printf(sock, "6 0  %11.10e %11.10e %11.10e\n", x0, y1, z1);
    Vio_printf(sock, "7 0  %11.10e %11.10e %11.10e\n", x1, y1, z1);
    Vio_printf(sock, "];\n");

    /* Write out the simplices */
    Vio_printf(sock, "simp=[\n");
    Vio_printf(sock, "0 0 0    0 1 0 1   0 5 1 2\n");
    Vio_printf(sock, "1 0 0    0 1 1 0   0 5 2 4\n");
    Vio_printf(sock, "2 0 0    0 1 0 1   1 5 3 2\n");
    Vio_printf(sock, "3 0 0    0 1 0 1   3 5 7 2\n");
    Vio_printf(sock, "4 0 0    1 1 0 0   2 5 7 6\n");
    Vio_printf(sock, "5 0 0    1 1 0 0   2 5 6 4\n");
    Vio_printf(sock, "];\n");

    /* Write out the MCSF footer */
    Vio_printf(sock, "mcsf_end=1;\n");

    /* Close the socket */
    Vio_dtor(&sock);
}

