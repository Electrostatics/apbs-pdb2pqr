/** @defgroup Vopot Vopot class
 *  @brief  Potential oracle for Cartesian mesh data
 */

/**
 *  @file    vopot.h
 *  @ingroup Vopot
 *  @author  Nathan Baker
 *  @brief   Potential oracle for Cartesian mesh data
 *  @version $Id$
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

#ifndef _VOPOT_H_
#define _VOPOT_H_

#include "maloc/maloc.h"
#include "apbs/vunit.h"
#include "apbs/vgrid.h"
#include "apbs/vpbe.h"

/**
 *  @struct  Vopot
 *  @ingroup Vopot
 *  @author  Nathan Baker
 *  @brief   Electrostatic potential oracle for Cartesian mesh data
 */
struct Vopot {

    Vgrid *grid;  /**< Grid object containing potential data (in units kT/e) */
    Vpbe   *pbe;  /**< Pointer to PBE object */
    int bcfl;     /**< Boundary condition flag for returning potential values
                   * at points off the grid.  0 is zero potential, 1 is single
                   * sphere Debye-Huckel approximation, and 2 is multiple
                   * sphere Debye-Huckel approximation */
};

/** @typedef Vopot
 *  @ingroup Vopot
 *  @brief   Declaration of the Vopot class as the Vopot structure
 */
typedef struct Vopot Vopot;

/** @brief   Construct Vopot object with values obtained from Vpmg_readDX (for
 *           example)
 *  @ingroup Vopot
 *  @author  Nathan Baker
 *  @param   grid  Grid object containing potential data (in units kT/e)
 *  @param   pbe   Pointer to Vpbe object for parameters
 *  @param   bcfl  Boundary condition to use for potential values off the grid
 *                 \li 0:  Zero potential
 *                 \li 1:  Single sphere Debye-Huckel approximation
 *                 \li 2:  Multiple sphere Debye-Huckel approximation
 *  @returns Newly allocated and initialized Vopot object
 */
VEXTERNC Vopot*  Vopot_ctor(Vgrid *grid, Vpbe *pbe, int bcfl);

/** @brief   Initialize Vopot object with values obtained from Vpmg_readDX (for
 *           example)
 *  @ingroup Vopot
 *  @author  Nathan Baker
 *  @param   thee  Pointer to newly allocated Vopot object
 *  @param   grid  Grid object containing potential data (in units kT/e)
 *  @param   pbe   Pointer to Vpbe object for parameters
 *  @param   bcfl  Boundary condition to use for potential values off the grid
 *                 \li 0:  Zero potential
 *                 \li 1:  Single sphere Debye-Huckel approximation
 *                 \li 2:  Multiple sphere Debye-Huckel approximation
 *  @returns 1 if successful, 0 otherwise
 */
VEXTERNC int Vopot_ctor2(Vopot *thee, Vgrid *grid, Vpbe *pbe, int bcfl);

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
