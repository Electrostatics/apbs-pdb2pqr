/** @defgroup FEMparm FEMparm class
 *  @brief    Parameter structure for FEM-specific variables from input files
 */

/**
 *  @file     femparm.h
 *  @ingroup  FEMparm
 *  @brief    Contains declarations for class FEMparm
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


#ifndef _FEMPARM_H_
#define _FEMPARM_H_

/* Generic header files */
#include "maloc/maloc.h"
#include "apbs/vhal.h"

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

