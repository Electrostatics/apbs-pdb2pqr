/** @defgroup Vopot Vopot class
 *  @brief  Potential oracle for Cartesian mesh data
 */

/**
 *  @file    vopot.h
 *  @ingroup Vopot
 *  @author  Nathan Baker
 *  @brief   Potential oracle for Cartesian mesh data
 *  @version $Id$
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

#ifndef _VOPOT_H_
#define _VOPOT_H_

/* Generic headers */
#include "maloc/maloc.h"
#include "apbs/vhal.h"

/* Specific headers */
#include "apbs/vatom.h"
#include "apbs/valist.h"
#include "apbs/vmgrid.h"
#include "apbs/vunit.h"
#include "apbs/vpbe.h"
#include "apbs/pbeparm.h"

/**
 *  @ingroup Vopot
 *  @author  Nathan Baker
 *  @brief   Electrostatic potential oracle for Cartesian mesh data
 */
struct sVopot {

    Vmgrid *mgrid;  /**< Multiple grid object containing potential data (in
                     * units kT/e) */
    Vpbe   *pbe;  /**< Pointer to PBE object */
    Vbcfl bcfl;  /**< Boundary condition flag for returning potential
                  * values at points off the grid. */
};

/** 
 *  @ingroup Vopot
 *  @brief   Declaration of the Vopot class as the Vopot structure
 */
typedef struct sVopot Vopot;

/** @brief   Construct Vopot object with values obtained from Vpmg_readDX (for
 *           example)
 *  @ingroup Vopot
 *  @author  Nathan Baker
 *  @param   mgrid  Multiple grid object containing potential data (in units
 *                  kT/e)
 *  @param   pbe   Pointer to Vpbe object for parameters
 *  @param   bcfl  Boundary condition to use for potential values off the grid
 *  @returns Newly allocated and initialized Vopot object
 */
VEXTERNC Vopot*  Vopot_ctor(Vmgrid *mgrid, Vpbe *pbe, Vbcfl bcfl);

/** @brief   Initialize Vopot object with values obtained from Vpmg_readDX (for
 *           example)
 *  @ingroup Vopot
 *  @author  Nathan Baker
 *  @param   thee  Pointer to newly allocated Vopot object
 *  @param   mgrid  Multiple grid object containing potential data (in units
 *                 kT/e)
 *  @param   pbe   Pointer to Vpbe object for parameters
 *  @param   bcfl  Boundary condition to use for potential values off the grid
 *  @returns 1 if successful, 0 otherwise
 */
VEXTERNC int Vopot_ctor2(Vopot *thee, Vmgrid *mgrid, Vpbe *pbe, Vbcfl bcfl);

/** @brief   Get potential value (from mesh or approximation) at a point
 *  @ingroup Vopot
 *  @author  Nathan Baker
 *  @param   thee  Vopot obejct
 *  @param   x     Point at which to evaluate potential
 *  @param   pot   Set to dimensionless potential (units kT/e) at point x
 *  @returns        1 if successful, 0 otherwise
 */
VEXTERNC int Vopot_pot(Vopot *thee, double x[3], double *pot);

/** @brief   Object destructor
 *  @ingroup Vopot
 *  @author  Nathan Baker
 *  @param   thee   Pointer to memory location of object to be destroyed
 */
VEXTERNC void Vopot_dtor(Vopot **thee);

/** @brief   FORTRAN stub object destructor
 *  @ingroup Vopot
 *  @author  Nathan Baker
 *  @param   thee   Pointer to object to be destroyed
 */
VEXTERNC void Vopot_dtor2(Vopot *thee);

/** @brief   Get second derivative values at a point
 *  @ingroup Vopot
 *  @author  Nathan Baker
 *  @param   thee   Pointer to Vopot object
 *  @param   pt     Location to evaluate second derivative
 *  @param   cflag  
 *             \li  0:  Reduced Maximal Curvature
 *             \li  1:  Mean Curvature (Laplace)
 *             \li  2:  Gauss Curvature
 *             \li  3:  True Maximal Curvature
 *  @param   curv   Set to specified curvature value
 *  @returns        1 if successful, 0 otherwise
 */
VEXTERNC int Vopot_curvature(Vopot *thee, double pt[3], int cflag, double
  *curv);

/** @brief   Get first derivative values at a point
 *  @ingroup Vopot
 *  @author  Nathan Baker
 *  @param   thee   Pointer to Vopot object
 *  @param   pt     Location to evaluate gradient
 *  @param   grad   Gradient
 *  @returns        1 if successful, 0 otherwise
 */
VEXTERNC int Vopot_gradient(Vopot *thee, double pt[3], double grad[3] );


#endif
