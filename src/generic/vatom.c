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
 * Copyright (c) 2002-2004.  Washington University in St. Louis.
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

VPUBLIC int Vatom_getPartID(Vatom *thee) {

   VASSERT(thee != VNULL);
   return thee->partID;

}

VPUBLIC void Vatom_setPartID(Vatom *thee, int partID) {

   VASSERT(thee != VNULL);
   thee->partID = partID;

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

VPUBLIC int Vatom_memChk(Vatom *thee) { return sizeof(Vatom); }

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

    for (i=0; i<3; i++) dest->position[i] = thee->position[i];
    dest->charge = thee->charge;
    dest->radius = thee->radius;
}

VPUBLIC void Vatom_copyFrom(Vatom *thee, Vatom *src) {

    Vatom_copyTo(src, thee);

}
