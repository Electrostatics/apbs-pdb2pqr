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
 * Washington University in St. Louis
 *
 * Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 2002-2003.  Washington University in St. Louis.
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
#include "apbs/valist.h"

VEMBED(rcsid="$Id$")

/* ///////////////////////////////////////////////////////////////////////////
// Class Valist: Inlineable methods
/////////////////////////////////////////////////////////////////////////// */
#if !defined(VINLINE_VATOM)

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Valist_getCenterCoord
// 
// Author:   Nathan Baker 
/////////////////////////////////////////////////////////////////////////// */
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



/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Valist_getAtomList
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Vatom* Valist_getAtomList(Valist *thee) {

  VASSERT(thee != NULL);
  return thee->atoms;

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Valist_getNumberAtoms
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Valist_getNumberAtoms(Valist *thee) {

  VASSERT(thee != NULL);
  return thee->number;

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Valist_getAtom
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Vatom* Valist_getAtom(Valist *thee, int i) {

  VASSERT(thee != NULL);
  VASSERT(i < thee->number);
  return &(thee->atoms[i]);

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Valist_memChk
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Valist_memChk(Valist *thee) {

  if (thee == NULL) return 0;
  return Vmem_bytes(thee->vmem);

}

#endif /* if !defined(VINLINE_VATOM) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Valist: Non-inlineable methods
/////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Valist_ctor
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Valist* Valist_ctor() {

    /* Set up the structure */
    Valist *thee = VNULL;
    thee = Vmem_malloc(VNULL, 1, sizeof(Valist));
    VASSERT( thee != VNULL);
    VASSERT( Valist_ctor2(thee));
 
    return thee;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Valist_ctor2
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Valist_ctor2(Valist *thee) {
  
    thee->atoms = VNULL;
    thee->number = 0;

    /* Initialize the memory management object */
    thee->vmem = Vmem_ctor("APBS:VALIST");

    return 1;    

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Valist_dtor
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Valist_dtor(Valist **thee)
{
    if ((*thee) != VNULL) {
        Valist_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(Valist), (void **)thee);
        (*thee) = VNULL;
    }
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Valist_dtor2
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Valist_dtor2(Valist *thee) {

    Vmem_free(thee->vmem, thee->number, sizeof(Vatom), (void **)&(thee->atoms));
    thee->atoms = VNULL;
    thee->number = 0;

    Vmem_dtor(&(thee->vmem));
} 

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Valist_readPQR
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Valist_readPQR(Valist *thee, const char *iodev, const char *iofmt,
  const char *thost, const char *fname) {

    FILE *pqrf;                 /* PQR file pointer */
    char line[101];             /* To hold lines from the PQR file */
    int maxl = 100;
    /* Counters */
    int i;
    /* Information to be gleaned from PQR file */
    double x, y, z, charge, radius;
    double pos[3];
 
    VASSERT(thee != VNULL);
    thee->number = 0;

    /* Make sure we're reading in an ASCII file */
    if ((!strcmp(iodev,"FILE")) && (!strcmp(iodev,"file"))) {
        Vnm_print(2, "Valist_readPQR:  This routine cannot read from type %s.\n", 
          iodev);
        return 0;
    }
    if ((!strcmp(iofmt,"ASC")) && (!strcmp(iofmt,"asc"))) {
        Vnm_print(2, "Valist_readPQR:  This routine cannot read type %s.\n",
          iofmt);
        return 0;
    }

    /* Open data files */
    pqrf = fopen(fname,"r");
    if (pqrf == NULL) {
        Vnm_print(2, "Valist_readPQR: Error opening %s\n",fname);
        return 0;
    }

    thee->center[0] = 0.;
    thee->center[1] = 0.;
    thee->center[2] = 0.;
    thee->maxcrd[0] = -VLARGE;
    thee->maxcrd[1] = -VLARGE;
    thee->maxcrd[2] = -VLARGE;
    thee->mincrd[0] = VLARGE;
    thee->mincrd[1] = VLARGE;
    thee->mincrd[2] = VLARGE;
    thee->maxrad = 0.;
    thee->charge = 0.;

    /* Now we read some lines and count the atoms. */
    while (1) {

        if (fgets(line,maxl,pqrf) == NULL) break;

        /* Check to see if we got an ATOM line */
        if ( !strncmp(line,"ATOM",4) ) (thee->number)++;
    
    }

#if defined(VDEBUG)
    Vnm_print(1, "Valist_readPQR: Counted %d atoms\n",thee->number);
    fflush(stdout);
#endif

    /* Allocate the necessary space for the atom array */
    thee->atoms = Vmem_malloc(thee->vmem, thee->number,(sizeof(Vatom)));
    VASSERT(thee->atoms != VNULL);
      

    /* Rewind the file pointer (use rewind to clear the error
       indicator, too) */
    rewind(pqrf);

    i = 0;
    while (1) {

        if (fgets(line,maxl,pqrf) == NULL) {
            Vnm_print(2, "Valist_readPQR: Read EOF instead of atom\n");
            fflush(stderr);
            return 0;
        }

        /* Check to see if we got an ATOM line */
        if ( !strncmp(line,"ATOM",4)) {
#if 0
            /* Try to parse a line with chain IDs and with an integer
             * residue ID */
	    if ((sscanf(line,"ATOM%*7d  %*4s%*4s%*5d %*2s %lf%lf%lf%lf%lf",
              &x,&y,&z,&charge,&radius) == 5)) {;}
            /* Try to parse a line with chain IDs and with an string
             * residue ID */
	    else if ((sscanf(line,"ATOM%*7d  %*4s%*4s%*5s %*2s %lf%lf%lf%lf%lf",
              &x,&y,&z,&charge,&radius) == 5)) {;}
            /* Try to parse a line without chain IDs and with an integer
             * residue ID */
#endif
	    if ((sscanf(line,"ATOM%*7d  %*4s%*4s%*5d    %lf%lf%lf%lf%lf",
              &x,&y,&z,&charge,&radius) == 5)) {;}
            /* Try to parse a line without chain IDs and with a string
             * residue ID */
            else if ((sscanf(line,"ATOM%*7d  %*4s%*4s%*5s    %lf%lf%lf%lf%lf",
              &x,&y,&z,&charge,&radius) == 5)) {;}
            else {
                Vnm_print(2, "Valist_readPQR:  FATAL sscanf (formatting) \
error reading: \n    %s\n", line);
                return 0;
            }

            if (x < thee->mincrd[0]) thee->mincrd[0] = x;
            if (y < thee->mincrd[1]) thee->mincrd[1] = y;
            if (z < thee->mincrd[2]) thee->mincrd[2] = z;
            if (x > thee->maxcrd[0]) thee->maxcrd[0] = x;
            if (y > thee->maxcrd[1]) thee->maxcrd[1] = y;
            if (z > thee->maxcrd[2]) thee->maxcrd[2] = z;
            if (radius > thee->maxrad) thee->maxrad = radius;
            thee->charge = thee->charge + charge;

            /* Put it in the atom list */
            pos[0] = x;
            pos[1] = y;
            pos[2] = z;
            Vatom_setPosition(&(thee->atoms)[i],pos);
            Vatom_setCharge(&(thee->atoms)[i],charge);
            Vatom_setRadius(&(thee->atoms)[i],radius);

            /* Update the number of atoms we've found */
            i++;

            /* Have we gotten all the entries? */
            if (i == thee->number)  break;

        } else {
            Vnm_print(2, "Valist_readPQR:  IGNORED line: \n    %s\n", line);
        }
    } /* while(1) */

    thee->center[0] = 0.5*(thee->maxcrd[0] + thee->mincrd[0]);
    thee->center[1] = 0.5*(thee->maxcrd[1] + thee->mincrd[1]);
    thee->center[2] = 0.5*(thee->maxcrd[2] + thee->mincrd[2]);

    return 1;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Valist_buildMesh
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
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

