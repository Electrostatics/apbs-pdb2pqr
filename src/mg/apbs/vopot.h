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
 * Nathan A. Baker (nbaker@wasabi.ucsd.edu)
 * Dept. of Chemistry and Biochemistry
 * University of California, San Diego 
 *
 * Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 1999-2002.  Nathan A. Baker.  All Rights Reserved.
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
 *
 * @endverbatim
 */

#ifndef _VOPOT_H_
#define _VOPOT_H_

#include "maloc/maloc.h"
#include "apbs/vunit.h"
#include "apbs/vpbe.h"

VEMBED(rcsid="$Id$")

/**
 *  @struct  Vopot
 *  @ingroup Vopot
 *  @author  Nathan Baker
 *  @brief   Electrostatic potential oracle for Cartesian mesh data
 */
struct Vopot {

    int nx;       /**< Number grid points in x direction */
    int ny;       /**< Number grid points in y direction */
    int nz;       /**< Number grid points in z direction */
    double hx;    /**< Grid spacing in x direction */
    double hy;    /**< Grid spacing in y direction */
    double hzed;  /**< Grid spacing in z direction */
    double xmin;  /**< x coordinate of lower grid corner */
    double ymin;  /**< y coordinate of lower grid corner */
    double zmin;  /**< z coordinate of lower grid corner */
    double *data; /**< nx*ny*nz array of potential data (in units kT/e) */
    Vpbe   *pbe;  /**< Pointer to PBE object */
};

/** @typedef Vopot
 *  @ingroup Vopot
 *  @brief   Declaration of the Vopot class as the Vopot structure
 */
typedef struct Vopot Vopot;

/** @def     VOPOT_BCFL  Controls behavior of potential outside mesh, 0 =>
 *                       zero potential, 1 => single Debye-Huckel sphere
 *                       approximation, 2 => multiple Debye-Huckel sphere
 *                       approximation (slow)
 *  @ingroup Vopot
 */
#define VOPOT_BCFL 1

/** @brief   Construct Vopot object
 *  @ingroup Vopot
 *  @author  Nathan Baker
 *  @param   nx    Number grid points in x direction
 *  @param   ny    Number grid points in y direction
 *  @param   nz    Number grid points in z direction
 *  @param   hx    Grid spacing in x direction
 *  @param   hy    Grid spacing in y direction
 *  @param   hzed  Grid spacing in z direction
 *  @param   xmin  x coordinate of lower grid corner
 *  @param   ymin  y coordinate of lower grid corner
 *  @param   zmin  z coordinate of lower grid corner
 *  @param   data  nx*ny*nz array of potential data (in units kT/e)
 *  @param   pbe   Pointer to Vpbe object for parameters
 *  @returns Newly allocated and initialized Vopot object
 */
VEXTERNC Vopot*  Vopot_ctor(int nx, int ny, int nz, 
                  double hx, double hy, double hzed,
                  double xmin, double ymin, double zmin,
                  double *data, Vpbe *pbe);

/** @brief   Initialize Vopot object
 *  @ingroup Vopot
 *  @author  Nathan Baker
 *  @param   thee  Pointer to newly allocated Vopot object
 *  @param   nx    Number grid points in x direction
 *  @param   ny    Number grid points in y direction
 *  @param   nz    Number grid points in z direction
 *  @param   hx    Grid spacing in x direction
 *  @param   hy    Grid spacing in y direction
 *  @param   hzed  Grid spacing in z direction
 *  @param   xmin  x coordinate of lower grid corner
 *  @param   ymin  y coordinate of lower grid corner
 *  @param   zmin  z coordinate of lower grid corner
 *  @param   data  nx*ny*nz array of potential data (in units kT/e)
 *  @param   pbe   Pointer to Vpbe object for parameters
 *  @returns Newly allocated and initialized Vopot object
 */
VEXTERNC int Vopot_ctor2(Vopot *thee, int nx, int ny, int nz, 
                  double hx, double hy, double hzed,
                  double xmin, double ymin, double zmin,
                  double *data, Vpbe *pbe);

/** @brief   Get potential value (from mesh or approximation) at a point
 *  @ingroup Vopot
 *  @author  Nathan Baker
 *  @param   thee  Vopot obejct
 *  @param   x     Point at which to evaluate potential
 *  @returns pot   Dimensionless potential (units kT/e) at point x
 */
VEXTERNC double Vopot_pot(Vopot *thee, double x[3]);

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

#endif
