/** @defgroup Vmgrid Vmgrid class
 *  @brief    Oracle for Cartesian mesh data
 */

/**
 *  @file    vmgrid.h
 *  @ingroup Vmgrid
 *  @author  Nathan Baker
 *  @brief   Multiresolution oracle for Cartesian mesh data
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

#ifndef _VMGRID_H_
#define _VMGRID_H_

#include "maloc/maloc.h"
#include "apbs/vgrid.h"
#include "apbs/apbs.h"

/** @def VMGRIDMAX   The maximum number of levels in the grid hiearchy
 *  @ingroup Vmgrid
 */
#define VMGRIDMAX 20


/**
 *  @struct  Vmgrid
 *  @ingroup Vmgrid
 *  @author  Nathan Baker
 *  @brief   Multiresoltion oracle for Cartesian mesh data
 */
struct Vmgrid {

    int ngrids;                /**< Number of grids in hiearchy */
    Vgrid *grids[VMGRIDMAX];   /**< Grids in hiearchy.  Our convention will be
                                *   to have the finest grid first, however,
                                *   this will not be enforced as it may be
                                *   useful to search multiple grids for
                                *   parallel datasets, etc. */
};

/** @typedef Vmgrid
 *  @ingroup Vmgrid
 *  @brief   Declaration of the Vmgrid class as the Vgmrid structure
 */
typedef struct Vmgrid Vmgrid;

/** @brief   Construct Vmgrid object 
 *  @ingroup Vmgrid
 *  @author  Nathan Baker
 *  @returns Newly allocated and initialized Vmgrid object
 */
VEXTERNC Vmgrid*  Vmgrid_ctor();

/** @brief   Initialize Vmgrid object 
 *  @ingroup Vmgrid
 *  @author  Nathan Baker
 *  @param   thee Newly allocated Vmgrid object
 *  @returns Newly allocated and initialized Vmgrid object
 */
VEXTERNC int Vmgrid_ctor2(Vmgrid *thee);

/** @brief   Get potential value (from mesh or approximation) at a point
 *  @ingroup Vmgrid
 *  @author  Nathan Baker
 *  @param   thee  Vmgrid obejct
 *  @param   x     Point at which to evaluate potential
 *  @param   value Value of data at point x
 *  @return  1 if successful, 0 if off grid
 */
VEXTERNC int Vmgrid_value(Vmgrid *thee, double x[3], double *value);

/** @brief   Object destructor
 *  @ingroup Vmgrid
 *  @author  Nathan Baker
 *  @param   thee   Pointer to memory location of object to be destroyed
 */
VEXTERNC void Vmgrid_dtor(Vmgrid **thee);

/** @brief   FORTRAN stub object destructor
 *  @ingroup Vmgrid
 *  @author  Nathan Baker
 *  @param   thee   Pointer to object to be destroyed
 */
VEXTERNC void Vmgrid_dtor2(Vmgrid *thee);

/** @brief   Add a grid to the hierarchy
 *  @ingroup Vmgrid
 *  @author  Nathan Baker
 *  @param   thee   Pointer to object to be destroyed
 *  @param   grid   Grid to be added.  As mentioned above, we would prefer to
 *           have the finest grid added first, next-finest second, ...,
 *           coarsest last -- this is how the grid will be searched when
 *           looking up values for points.  However, this is not enforced to
 *           provide flexibility for cases where the dataset is decomposed into
 *           disjoint partitions, etc.
 *  @returns 1 if successful, 0 otherwise
 */
VEXTERNC int Vmgrid_addGrid(Vmgrid *thee, Vgrid *grid);


/** @brief   Get second derivative values at a point
 *  @ingroup Vmgrid
 *  @author  Nathan Baker (wrapper for Vgrid routine by Steve Bond)
 *  @param   thee   Pointer to Vmgrid object
 *  @param   pt     Location to evaluate second derivative
 *  @param   cflag  
 *             \li  0:  Reduced Maximal Curvature
 *             \li  1:  Mean Curvature (Laplace)
 *             \li  2:  Gauss Curvature
 *             \li  3:  True Maximal Curvature
 *  @param   curv Specified curvature value
 *  @return  1 if successful, 0 if off grid
 */
VEXTERNC int Vmgrid_curvature(Vmgrid *thee, double pt[3], int cflag, 
  double *curv);

/** @brief   Get first derivative values at a point
 *  @ingroup Vmgrid
 *  @author  Nathan Baker and Steve Bond
 *  @param   thee   Pointer to Vmgrid object
 *  @param   pt     Location to evaluate gradient
 *  @param   grad   Gradient
 *  @return  1 if successful, 0 if off grid
 */
VEXTERNC int Vmgrid_gradient(Vmgrid *thee, double pt[3], double grad[3] );

/** @brief   Get specific grid in hiearchy
 *  @ingroup Vmgrid
 *  @author  Nathan Baker 
 *  @param   thee   Pointer to Vmgrid object
 *  @param   num    Number of grid in hiearchy 
 *  @return  Pointer to specified grid
 */
VEXTERNC Vgrid* Vmgrid_getGridByNum(Vmgrid *thee, int num);

/** @brief   Get grid in hiearchy which contains specified point or VNULL
 *  @ingroup Vmgrid
 *  @author  Nathan Baker 
 *  @param   thee   Pointer to Vmgrid object
 *  @param   pt     Point to check
 *  @return  Pointer to specified grid
 */
VEXTERNC Vgrid* Vmgrid_getGridByPoint(Vmgrid *thee, double pt[3]);

#endif

