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
// File:     valist.c
//
// Purpose:  Class Valist: methods.
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */

#include "apbscfg.h"
#include "apbs/valist.h"

VEMBED(rcsid="$Id$")

/* ///////////////////////////////////////////////////////////////////////////
// Class Valist: Inlineable methods
/////////////////////////////////////////////////////////////////////////// */
#if !defined(VINLINE_VATOM)

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Valist_getAtomList
//
// Purpose:  Get atom list
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
// Purpose:  Get number of atoms in atom list
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
// Purpose:  Get pointer to atom i
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
// Purpose:  Get total memory (in bytes) allocated for this object
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
// Purpose:  Construct the atom list object
//
// Notes:    This routine sets up data members via file I/O.
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
// Purpose:  Construct the atom list object
//
// Notes:    This routine sets up data members of the class.
//           Broken into two parts for FORTRAN users.
//
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
// Purpose:  Destroy the atom list object
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
// Purpose:  Destroy the atom list object
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
// Purpose:  Fill atom list with information from a PQR file
//           A PQR file has PDB structure with charge and radius in the
//           last two columns (instead of weight and occupancy)
//
// Returns:  1 if successful
//
// Notes:    The PDB reader routine was borrowed from Phil Hunenberger
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
        fprintf(stderr,"Valist_readPQR: Error opening %s\n",fname);
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
    printf("Valist_readPQR: Counted %d atoms\n",thee->number);
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
            fprintf(stderr,"Valist_readPQR: Read EOF instead of atom\n");
            fflush(stderr);
            return 0;
        }

        /* Check to see if we got an ATOM line */
        if ( !strncmp(line,"ATOM",4)) {
            if ( (sscanf(line,"ATOM%*7d  %*4s%*4s%*5d    %lf%lf%lf%lf%lf", 
                   &x,&y,&z,&charge,&radius) == 5) != 1) {
                if ( (sscanf(line,"ATOM%*7d  %*4s%*4s%*5s    %lf%lf%lf%lf%lf",
                   &x,&y,&z,&charge,&radius) == 5) != 1) {
                    fprintf(stderr,"Valist_readPQR:  FATAL sscanf (formatting) error reading: \n    %s\n", line);
                    return 0;
                }
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

        }  /* !strncmp(line,"ATOM",4) */
    } /* while(1) */

    thee->center[0] = 0.5*(thee->maxcrd[0] + thee->mincrd[0]);
    thee->center[1] = 0.5*(thee->maxcrd[1] + thee->mincrd[1]);
    thee->center[2] = 0.5*(thee->maxcrd[2] + thee->mincrd[2]);

    return 1;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Valist_buildMesh
//
// Purpose:   Build a cuboid 3D mesh to surround the molecule contained in this
//            Valist.  The mesh will have rectagular sides and consist of 6
//            simplices/8 vertices with all boundaries Dirichlet.  The mesh
//            will be written, in MCSF format, to the output "path" specified
//            by the arguments below.
//
// Arguments: size  The factor by which the mesh is larger than the
//                  biomolecule.  In other words, if the smallest box
//                  containing the protein is dx x dy x dz, then the mesh will
//                  be (size*dx) x (size*dy) x (size*dz).  Clearly, size > 1.
//            iodev Where we write the data: 
//                    FILE -- to some file
//                    INET -- to an INET socket
//                    BUFF -- to a buffer in memory
//                    UNIX -- to a UNIX domain socket
//            iofmt Data format that we write:
//                    ASC -- ASCII
//                    XDR -- Network byte order format
//            thost The hostname to which we may connect the socket
//            fname The file name
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

