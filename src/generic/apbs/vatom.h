/** @defgroup Vatom Vatom class
 *  @brief  Atom class for interfacing APBS with PDB files
 */

/**
 *  @file     vatom.h
 *  @ingroup  Vatom
 *  @brief    Contains declarations for class Vatom
 *  @version  $Id$
 *  @author   Nathan A. Baker
 *
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

#ifndef _VATOM_H_
#define _VATOM_H_

#include "maloc/maloc.h"
#include "apbs/vhal.h"

/** 
 *  @ingroup Vatom
 *  @author  Nathan Baker
 *  @brief   Contains public data members for Vatom class/module
 */
struct sVatom {

    double position[3];  /**< Atomic position */
    double radius;  /**< Atomic radius   */
    double charge;  /**< Atomic charge   */
    double partID;  /**< Partition value for assigning atoms to particular
                     * processors and/or partitions   */
    int id;  /**< Atomic ID; this should be a unique non-negative integer
              * assigned based on the index of the atom in a Valist atom 
              * array */

};

/** 
 *  @ingroup Vatom
 *  @brief   Declaration of the Vatom class as the Vatom structure
 */
typedef struct sVatom Vatom;

#if !defined(VINLINE_VATOM)

    /** @brief   Get atomic position
     *  @ingroup Vatom
     *  @author  Nathan Baker
     *  @param   thee  Vatom object
     *  @returns Pointer to 3*double array of atomic coordinates (in &Aring;)
     */
    VEXTERNC double* Vatom_getPosition(Vatom *thee);

    /** @brief   Set atomic radius
     *  @ingroup Vatom
     *  @author  Nathan Baker
     *  @param   thee    Vatom object
     *  @param   radius  Atomic radius (in &Aring;)
     */
    VEXTERNC void    Vatom_setRadius(Vatom *thee, double radius);

    /** @brief   Get atomic position
     *  @ingroup Vatom
     *  @author  Nathan Baker
     *  @param   thee  Vatom object
     *  @returns Atomic radius (in &Aring;)
     */
    VEXTERNC double  Vatom_getRadius(Vatom *thee);

    /** @brief   Set partition ID
     *  @ingroup Vatom
     *  @author  Nathan Baker
     *  @param   thee    Vatom object
     *  @param   partID  Partition ID; a negative value means this atom is not
     *                   assigned to any partition
     */
    VEXTERNC void    Vatom_setPartID(Vatom *thee, int partID);

    /** @brief   Get partition ID
     *  @ingroup Vatom
     *  @author  Nathan Baker
     *  @param   thee    Vatom object
     *  @return  Partition ID; a negative value means this atom is not
     *           assigned to any partition
     */
    VEXTERNC double     Vatom_getPartID(Vatom *thee);

    /** @brief   Set atom ID
     *  @ingroup Vatom
     *  @author  Nathan Baker
     *  @param   thee    Vatom object
     *  @param   id  Unique non-negative number
     */
    VEXTERNC void Vatom_setAtomID(Vatom *thee, int id);

    /** @brief   Get atom ID
     *  @ingroup Vatom
     *  @author  Nathan Baker
     *  @param   thee    Vatom object
     *  @return  Unique non-negative number
     */
    VEXTERNC double Vatom_getPartID(Vatom *thee);

    /** @brief   Set atomic charge
     *  @ingroup Vatom
     *  @author  Nathan Baker
     *  @param   thee    Vatom object
     *  @param   charge  Atom partial charge (in e)
     */
    VEXTERNC void    Vatom_setCharge(Vatom *thee, double charge);

    /** @brief   Get atomic charge
     *  @ingroup Vatom
     *  @author  Nathan Baker
     *  @param   thee    Vatom object
     *  @return  Atom partial charge (in e)
     */
    VEXTERNC double  Vatom_getCharge(Vatom *thee);

    /** @brief   Return the memory used by this structure (and its contents)
     *           in bytes
     *  @ingroup Vatom
     *  @author  Nathan Baker
     *  @param   thee  Vpmg object
     *  @return  The memory used by this structure and its contents in bytes
     */
    VEXTERNC unsigned long int Vatom_memChk(Vatom *thee);

#else /* if defined(VINLINE_VATOM) */
#   define Vatom_getPosition(thee) ((thee)->position)
#   define Vatom_setRadius(thee, tRadius) ((thee)->radius = (tRadius))
#   define Vatom_getRadius(thee) ((thee)->radius)
#   define Vatom_setPartID(thee, tpartID) ((thee)->partID = (double)(tpartID))
#   define Vatom_getPartID(thee) ((thee)->partID)
#   define Vatom_setAtomID(thee, tatomID) ((thee)->id = (tatomID))
#   define Vatom_getAtomID(thee) ((thee)->id)
#   define Vatom_setCharge(thee, tCharge) ((thee)->charge = (tCharge))
#   define Vatom_getCharge(thee) ((thee)->charge)
#   define Vatom_memChk(thee) (sizeof(Vatom))
#endif /* if !defined(VINLINE_VATOM) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vatom: Non-Inlineable methods (vatom.c)
/////////////////////////////////////////////////////////////////////////// */

/** @brief   Constructor for the Vatom class 
 *  @author  Nathan Baker
 *  @ingroup Vatom
 *  @returns Pointer to newly allocated Vatom object
 */
VEXTERNC Vatom* Vatom_ctor();

/** @brief   FORTRAN stub constructor for the Vatom class 
 *  @author  Nathan Baker
 *  @ingroup Vatom
 *  @param   thee Pointer to Vatom allocated memory location
 *  @returns 1 if succesful, 0 otherwise
 */
VEXTERNC int     Vatom_ctor2(Vatom *thee);

/** @brief   Object destructor
 *  @ingroup Vatom
 *  @author  Nathan Baker
 *  @param   thee   Pointer to memory location of object to be destroyed
 */
VEXTERNC void    Vatom_dtor(Vatom **thee);

/** @brief   FORTRAN stub object destructor
 *  @ingroup Vatom
 *  @author  Nathan Baker
 *  @param   thee   Pointer to object to be destroyed
 */
VEXTERNC void    Vatom_dtor2(Vatom *thee);

/** @brief   Set the atomic position
 *  @ingroup Vatom
 *  @author  Nathan Baker
 *  @param   thee   Vatom object to be modified
 *  @param   position  Coordinates (in &Aring;)
 */
VEXTERNC void   Vatom_setPosition(Vatom *thee, double position[3]);

/**
 * @brief  Copy information to another atom
 * @ingroup  Vatom
 * @author  Nathan Baker
 * @param  thee Source for atom information
 * @param  dest Destination for atom information
 */
VEXTERNC void Vatom_copyTo(Vatom *thee, Vatom *dest);

/**
 * @brief  Copy information to another atom
 * @ingroup  Vatom
 * @author  Nathan Baker
 * @param  thee Destination for atom information
 * @param  src Source for atom information
 */
VEXTERNC void Vatom_copyFrom(Vatom *thee, Vatom *src);

#endif /* ifndef _VATOM_H_ */
