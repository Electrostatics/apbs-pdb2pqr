/**
 *  @file    vatom.c
 *  @ingroup Vatom
 *  @author  Nathan Baker
 *  @brief   Class Vatom methods
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
 * Copyright (c) 2002-2007.  Washington University in St. Louis.
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
 * code included in releases of ISIM, Ion Simulator Interface, PMV, PyMOL
 * SMOL, VMD, and Vision. Such combined software may be linked with APBS and 
 * redistributed together in original or modified form as mere aggregation
 * without requirement that the entire work be under the scope of the GNU 
 * General Public License. This special exception permission is also extended
 * to any software listed in the SPECIAL GPL EXCEPTION clauses by the PMG,
 * FEtk, MC, or MALOC libraries.
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
#include "apbs/vatom.h"

VEMBED(rcsid="$Id$")

#if !defined(VINLINE_VATOM)

VPUBLIC double *Vatom_getPosition(Vatom *thee) { 

   VASSERT(thee != VNULL);
   return thee->position; 

}

VPUBLIC double Vatom_getPartID(Vatom *thee) {

   VASSERT(thee != VNULL);
   return thee->partID;

}

VPUBLIC void Vatom_setPartID(Vatom *thee, int partID) {

   VASSERT(thee != VNULL);
   thee->partID = (double)partID;

}

VPUBLIC double Vatom_getAtomID(Vatom *thee) {

   VASSERT(thee != VNULL);
   return thee->atomID;

}

VPUBLIC void Vatom_setAtomID(Vatom *thee, int atomID) {

   VASSERT(thee != VNULL);
   thee->atomID = atomID;

}

VPUBLIC void Vatom_setRadius(Vatom *thee, double radius) { 

   VASSERT(thee != VNULL);
   thee->radius = radius; 

}

VPUBLIC double Vatom_getRadius(Vatom *thee) { 

   VASSERT(thee != VNULL);
   return thee->radius; 

}

VPUBLIC void Vatom_setCharge(Vatom *thee, double charge) { 

   VASSERT(thee != VNULL);
   thee->charge = charge; 

}

VPUBLIC double Vatom_getCharge(Vatom *thee) { 

   VASSERT(thee != VNULL);
   return thee->charge; 

}

VPUBLIC unsigned long int Vatom_memChk(Vatom *thee) { return sizeof(Vatom); }

#endif /* if !defined(VINLINE_VATOM) */

VPUBLIC Vatom* Vatom_ctor() {

    /* Set up the structure */
    Vatom *thee = VNULL;
    thee = (Vatom *)Vmem_malloc( VNULL, 1, sizeof(Vatom) );
    VASSERT( thee != VNULL);
    VASSERT( Vatom_ctor2(thee));

    return thee;
}

VPUBLIC int Vatom_ctor2(Vatom *thee) { 
    thee->partID = -1;
    return 1; 
}

VPUBLIC void Vatom_dtor(Vatom **thee) {
    if ((*thee) != VNULL) {
        Vatom_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(Vatom), (void **)thee);
        (*thee) = VNULL;
    }
}

VPUBLIC void Vatom_dtor2(Vatom *thee) { ; }

VPUBLIC void Vatom_setPosition(Vatom *thee, double position[3]) { 

   VASSERT(thee != VNULL);
   (thee->position)[0] = position[0]; 
   (thee->position)[1] = position[1]; 
   (thee->position)[2] = position[2]; 

}

VPUBLIC void Vatom_copyTo(Vatom *thee, Vatom *dest) {

    int i;

    VASSERT(thee != VNULL);
    VASSERT(dest != VNULL);

    memcpy(dest, thee, sizeof(Vatom));

}

VPUBLIC void Vatom_copyFrom(Vatom *thee, Vatom *src) {

    Vatom_copyTo(src, thee);

}

VPUBLIC void Vatom_setResName(Vatom *thee, char resName[VMAX_RECLEN]) { 
	
	VASSERT(thee != VNULL);
	strcpy(thee->resName, resName);
	
}

VPUBLIC void Vatom_getResName(Vatom *thee, char resName[VMAX_RECLEN]) { 
	
	
	VASSERT(thee != VNULL);
	strcpy(resName,thee->resName);
	
}

VPUBLIC void Vatom_setAtomName(Vatom *thee, char atomName[VMAX_RECLEN]) { 
	
	VASSERT(thee != VNULL);
	strcpy(thee->atomName, atomName);
	
}

VPUBLIC void Vatom_getAtomName(Vatom *thee, char atomName[VMAX_RECLEN]) { 
	
	VASSERT(thee != VNULL);
	strcpy(atomName,thee->atomName);
	
}

#if defined(WITH_TINKER)

VPUBLIC void Vatom_setDipole(Vatom *thee, double dipole[3]) { 

   VASSERT(thee != VNULL);
   (thee->dipole)[0] = dipole[0]; 
   (thee->dipole)[1] = dipole[1]; 
   (thee->dipole)[2] = dipole[2]; 

}

VPUBLIC void Vatom_setQuadrupole(Vatom *thee, double quadrupole[9]) { 

   int i;
   VASSERT(thee != VNULL);
   for (i=0; i<9; i++)  (thee->quadrupole)[i] = quadrupole[i];
}

VPUBLIC void Vatom_setInducedDipole(Vatom *thee, double dipole[3]) { 

   VASSERT(thee != VNULL);
   (thee->inducedDipole)[0] = dipole[0]; 
   (thee->inducedDipole)[1] = dipole[1]; 
   (thee->inducedDipole)[2] = dipole[2]; 
}

VPUBLIC void Vatom_setNLInducedDipole(Vatom *thee, double dipole[3]) { 

   VASSERT(thee != VNULL);
   (thee->nlInducedDipole)[0] = dipole[0]; 
   (thee->nlInducedDipole)[1] = dipole[1]; 
   (thee->nlInducedDipole)[2] = dipole[2]; 

}

VPUBLIC double *Vatom_getDipole(Vatom *thee) { 

   VASSERT(thee != VNULL);
   return thee->dipole;

}

VPUBLIC double *Vatom_getQuadrupole(Vatom *thee) { 

   VASSERT(thee != VNULL);
   return thee->quadrupole;

}

VPUBLIC double *Vatom_getInducedDipole(Vatom *thee) { 

   VASSERT(thee != VNULL);
   return thee->inducedDipole;

}
 
VPUBLIC double *Vatom_getNLInducedDipole(Vatom *thee) { 

   VASSERT(thee != VNULL);
   return thee->nlInducedDipole;

}

VPUBLIC void Vatom_setResName(Vatom *thee, char resName[VMAX_ARGLEN]) { 
	
	VASSERT(thee != VNULL);
	strcpy(thee->resName, resName);
	
}

VPUBLIC void Vatom_getResName(Vatom *thee, char resName[VMAX_ARGLEN]) { 
	
	
	VASSERT(thee != VNULL);
	strcpy(resName,thee->resName);
	return 1; 
	
}

VPUBLIC void Vatom_setAtomName(Vatom *thee, char atomName[VMAX_ARGLEN]) { 
	
	VASSERT(thee != VNULL);
	strcpy(thee->atomName, atomName);
	
}

VPUBLIC void Vatom_getAtomName(Vatom *thee, char atomName[VMAX_ARGLEN]) { 
	
	VASSERT(thee != VNULL);
	strcpy(atomName,thee->atomName);
	return 1; 
	
}

#endif /* if defined(WITH_TINKER) */
