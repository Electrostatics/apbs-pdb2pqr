/** @defgroup FEMparm FEMparm class
 *  @brief    Parameter structure for FEM-specific variables from input files
 */

/**
 *  @file     femparm.h
 *  @ingroup  FEMparm
 *  @brief    Contains declarations for class FEMparm
 *  @version  $Id$
 *  @author   Nathan A. Baker
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
 * Copyright (c) 2002.  Washington University in St. Louis.
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


#ifndef _FEMPARM_H_
#define _FEMPARM_H_

#include "apbs/apbs.h"
#include "maloc/maloc.h"

/**
 *  @struct  FEMparm
 *  @ingroup FEMparm
 *  @author  Nathan Baker
 *  @brief   Parameter structure for FEM-specific variables from input files
 */
struct FEMparm {
   
    int parsed;         /**< Flag:  Has this structure been filled with
                         * anything other than * the default values? (0 = no,
                         * 1 = yes) */
};

/** @typedef FEMparm
 *  @ingroup FEMparm
 *  @brief   Declaration of the FEMparm class as the FEMparm structure
 */
typedef struct FEMparm FEMparm;

/* ///////////////////////////////////////////////////////////////////////////
// Class NOsh: Non-inlineable methods (nosh.c)
/////////////////////////////////////////////////////////////////////////// */

/** @brief   Construct FEMparm
 *  @ingroup FEMparm
 *  @author  Nathan Baker
 *  @returns Newly allocated and initialized Vpmgp object
 */
VEXTERNC FEMparm* FEMparm_ctor();

/** @brief   FORTRAN stub to construct FEMparm
 *  @ingroup FEMparm
 *  @author  Nathan Baker
 *  @param   thee Pointer to allocated FEMparm object
 *  @returns 1 if successful, 0 otherwise
 */
VEXTERNC int       FEMparm_ctor2(FEMparm *thee);

/** @brief   Object destructor
 *  @ingroup FEMparm
 *  @author  Nathan Baker
 *  @param   thee  Pointer to memory location of FEMparm object
 */
VEXTERNC void      FEMparm_dtor(FEMparm **thee);

/** @brief   FORTRAN stub for object destructor
 *  @ingroup FEMparm
 *  @author  Nathan Baker
 *  @param   thee  Pointer to FEMparm object
 */
VEXTERNC void      FEMparm_dtor2(FEMparm *thee);

/** @brief   Consistency check for parameter values stored in object
 *  @ingroup FEMparm
 *  @author  Nathan Baker
 *  @param   thee   FEMparm object
 *  @returns 1 if OK, 0 otherwise
 */
VEXTERNC int       FEMparm_check(FEMparm *thee);

#endif 

