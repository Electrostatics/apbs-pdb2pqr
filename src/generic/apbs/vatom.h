/** @defgroup Vatom Vatom class
 *  @brief  Atom class for interfacing APBS with PDB files
 */

/**
 *  @file     vatom.h
 *  @ingroup  Vatom
 *  @brief    Contains declarations for class Vatom
 *  @version  $Id$
 *  @author   Nathan A. Baker
 *  @attention
 *  @verbatim
 *
 * APBS -- Adaptive Poisson-Boltzmann Solver
 *
 * Nathan A. Baker (nbaker@wasabi.ucsd.edu)
 * Dept. of Chemistry and Biochemistry
 * University of California, San Diego
 *
 * Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 1999-2002. The Regents of the University of California 
 *                          (Regents).  All Rights Reserved.
 *
 * Permission to use, copy, modify, and distribute this software and its
 * documentation for educational, research, and not-for-profit purposes,
 * without fee and without a signed licensing agreement, is hereby granted,
 * provided that the above copyright notice, this paragraph and the
 * following two paragraphs appear in all copies, modifications, and
 * distributions.
 *
 * IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
 * REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE.  THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF
 * ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS".  REGENTS HAS NO OBLIGATION
 * TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 * MODIFICATIONS.
 *
 * @endverbatim
 */

#ifndef _VATOM_H_
#define _VATOM_H_

#include "maloc/maloc.h"
#include "apbs/vhal.h"

/** 
 *  @struct  Vatom
 *  @ingroup Vatom
 *  @author  Nathan Baker
 *  @brief   Contains public data members for Vatom class/module
 */
struct Vatom {

    double position[3];     /**< Atomic position */
    double radius;          /**< Atomic radius   */
    double charge;          /**< Atomic charge   */
    int partID;             /**< Partition ID for assigning atoms to particular
                             *   processors and/or partitions   */

};

/** @typedef Vatom
 *  @ingroup Vatom
 *  @brief   Declaration of the Vatom class as the Vatom structure
 */
typedef struct Vatom Vatom;

/* ///////////////////////////////////////////////////////////////////////////
// Class Vatom: Inlineable methods (vatom.c)
/////////////////////////////////////////////////////////////////////////// */

#if !defined(VINLINE_VATOM)

    /** @brief   Get atomic position
     *  @ingroup Vatom
     *  @author  Nathan Baker
     *  @param   thee  Vatom object
     *  @returns Pointer to 3*double array of atomic coordinates (in \f$\AA\f$)
     */
    VEXTERNC double* Vatom_getPosition(Vatom *thee);

    /** @brief   Set atomic radius
     *  @ingroup Vatom
     *  @author  Nathan Baker
     *  @param   thee    Vatom object
     *  @param   radius  Atomic radius (in \f$\AA\f$)
     */
    VEXTERNC void    Vatom_setRadius(Vatom *thee, double radius);

    /** @brief   Get atomic position
     *  @ingroup Vatom
     *  @author  Nathan Baker
     *  @param   thee  Vatom object
     *  @returns Atomic radius (in \f$\AA\f$)
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
    VEXTERNC int     Vatom_getPartID(Vatom *thee);

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
    VEXTERNC int     Vatom_memChk(Vatom *thee);

#else /* if defined(VINLINE_VATOM) */
#   define Vatom_getPosition(thee) ((thee)->position)
#   define Vatom_setRadius(thee, tRadius) ((thee)->radius = (tRadius))
#   define Vatom_getRadius(thee) ((thee)->radius)
#   define Vatom_setPartID(thee, tpartID) ((thee)->partID = (tpartID))
#   define Vatom_getPartID(thee) ((thee)->partID)
#   define Vatom_setCharge(thee, charge) ((thee)->charge = charge)
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
VEXTERNC Vatom * Vatom_ctor();

/** @brief   FORTRAN stub constructor for the Vatom class 
 *  @author  Nathan Baker
 *  @ingroup Vatom
 *  @param   Pointer to Vatom allocated memory location
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
 *  @param   position  Coordinates (in \f$\AA\f$)
 */
VEXTERNC void   Vatom_setPosition(Vatom *thee, double position[3]);

#endif /* ifndef _VATOM_H_ */
