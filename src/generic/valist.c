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
    thee = Vram_ctor( 1, sizeof(Valist) );
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
        Vram_dtor((Vram **)thee, 1, sizeof(Valist));
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

    Vram_dtor((Vram **)&(thee->atoms), thee->number, sizeof(Vatom));
    thee->atoms = VNULL;
    thee->number = 0;

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
    if ((thee->atoms = Vram_ctor(thee->number,(sizeof(Vatom)))) == VNULL) {
        fprintf(stderr,"Valist_readPQR: Failed to allocate atom array for %d atoms!\n",thee->number);
        fflush(stderr);
        return 0;
    }
      

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

