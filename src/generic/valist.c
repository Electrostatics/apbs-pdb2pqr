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
// Notes:    The PDB reader routine was borrowed from Phil Hunenberger
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Valist_readPQR(Valist *thee, char *path) {

    FILE *pqrf;                 /* PQR file pointer */
    char line[101];             /* To hold lines from the PQR file */
    int maxl = 100;
    /* Counters */
    int i;
    /* Information to be gleaned from PQR file */
    double x, y, z, charge, radius;
    Vec3 pos;

    VASSERT(thee != VNULL);
    thee->number = 0;

    /* Open data files */
    pqrf = fopen(path,"r");
    if (pqrf == NULL) {
        fprintf(stderr,"Valist_readPQR: Error opening %s\n",path);
        return 0;
    }

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
                fprintf(stderr,"Valist_readPQR: sscanf (formatting) error.\n");
                return 0;
            }
            

            pos[0] = x;
            pos[1] = y;
            pos[2] = z;

            /* Put it in the atom list */
            Vatom_setPosition(&(thee->atoms)[i],pos);
            Vatom_setCharge(&(thee->atoms)[i],charge);
            Vatom_setRadius(&(thee->atoms)[i],radius);

            /* Update the number of atoms we've found */
            i++;

            /* Have we gotten all the entries? */
            if (i == thee->number)  break;

        }  /* !strncmp(line,"ATOM",4) */
    } /* while(1) */

    return 1;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Valist_readPDB
//
// Purpose:  Fill atom list with information from a PDB file
//
// Notes:    The PDB reader routine was borrowed from Phil Hunenberger
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Valist_readPDB(Valist *thee, char *path, char *parameter_path) {

  VASSERT(thee != VNULL);

  fprintf(stderr,"Valist_readPDB: I haven't gotten around to writing the PDB reader yet.\n");
  fprintf(stderr,"Valist_readPDB: Until then, use awk and Valist_readPQR.\n");
  return 0;

}

