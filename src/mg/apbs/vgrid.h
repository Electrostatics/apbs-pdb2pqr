/** @defgroup Vgrid Vgrid class
 *  @brief    Oracle for Cartesian mesh data
 */

/**
 *  @file    vgrid.h
 *  @ingroup Vgrid
 *  @author  Nathan Baker and Steve Bond
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

#ifndef _VGRID_H_
#define _VGRID_H_

#include "maloc/maloc.h"

/**
 *  @struct  Vgrid
 *  @ingroup Vgrid
 *  @author  Nathan Baker
 *  @brief   Electrostatic potential oracle for Cartesian mesh data
 */
struct Vgrid {

    int nx;       /**< Number grid points in x direction */
    int ny;       /**< Number grid points in y direction */
    int nz;       /**< Number grid points in z direction */
    double hx;    /**< Grid spacing in x direction */
    double hy;    /**< Grid spacing in y direction */
    double hzed;  /**< Grid spacing in z direction */
    double xmin;  /**< x coordinate of lower grid corner */
    double ymin;  /**< y coordinate of lower grid corner */
    double zmin;  /**< z coordinate of lower grid corner */
    double xmax;  /**< x coordinate of upper grid corner */
    double ymax;  /**< y coordinate of upper grid corner */
    double zmax;  /**< z coordinate of upper grid corner */
    double *data; /**< nx*ny*nz array of data */
    int readdata; /**< flag indicating whether data was read from file */
    int ctordata; /**< flag indicating whether data was included at
                   *   construction */
    Vmem *mem;    /**< Memory manager object */
};

/** @typedef Vgrid
 *  @ingroup Vgrid
 *  @brief   Declaration of the Vgrid class as the Vgrid structure
 */
typedef struct Vgrid Vgrid;

#if !defined(VINLINE_VGRID)

    /** @brief   Return the memory used by this structure (and its contents)
     *           in bytes
     *  @ingroup Vgrid
     *  @author  Nathan Baker
     *  @param   thee  Vgrid object
     *  @return  The memory used by this structure and its contents in bytes
     */
    VEXTERNC int Vgrid_memChk(Vgrid *thee);

#else /* if defined(VINLINE_VGRID) */

    /** @brief   Return the memory used by this structure (and its contents)
     *           in bytes
     *  @ingroup Vgrid
     *  @author  Nathan Baker
     *  @param   thee  Vgrid object
     *  @return  The memory used by this structure and its contents in bytes
     */
#   define Vgrid_memChk(thee) (Vmem_bytes((thee)->vmem))

#endif /* if !defined(VINLINE_VPMG) */

/** @brief   Construct Vgrid object with values obtained from Vpmg_readDX (for
 *           example)
 *  @ingroup Vgrid
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
 *  @param   data  nx*ny*nz array of data.  This can be VNULL if you are
 *                 planning to read in data later with one of the read routines
 *  @returns Newly allocated and initialized Vgrid object
 */
VEXTERNC Vgrid*  Vgrid_ctor(int nx, int ny, int nz, 
                  double hx, double hy, double hzed,
                  double xmin, double ymin, double zmin,
                  double *data);

/** @brief   Initialize Vgrid object with values obtained from Vpmg_readDX (for
 *           example)
 *  @ingroup Vgrid
 *  @author  Nathan Baker
 *  @param   thee  Pointer to newly allocated Vgrid object
 *  @param   nx    Number grid points in x direction
 *  @param   ny    Number grid points in y direction
 *  @param   nz    Number grid points in z direction
 *  @param   hx    Grid spacing in x direction
 *  @param   hy    Grid spacing in y direction
 *  @param   hzed  Grid spacing in z direction
 *  @param   xmin  x coordinate of lower grid corner
 *  @param   ymin  y coordinate of lower grid corner
 *  @param   zmin  z coordinate of lower grid corner
 *  @param   data  nx*ny*nz array of data.  This can be VNULL if you are
 *                 planning to read in data later with one of the read routines
 *  @returns Newly allocated and initialized Vgrid object
 */
VEXTERNC int Vgrid_ctor2(Vgrid *thee, int nx, int ny, int nz, 
                  double hx, double hy, double hzed,
                  double xmin, double ymin, double zmin,
                  double *data);

/** @brief   Get potential value (from mesh or approximation) at a point
 *  @ingroup Vgrid
 *  @author  Nathan Baker
 *  @param   thee  Vgrid obejct
 *  @param   x     Point at which to evaluate potential
 *  @param   value Value of data at point x
 *  @return  1 if successful, 0 if off grid
 */
VEXTERNC int Vgrid_value(Vgrid *thee, double x[3], double *value);

/** @brief   Object destructor
 *  @ingroup Vgrid
 *  @author  Nathan Baker
 *  @param   thee   Pointer to memory location of object to be destroyed
 */
VEXTERNC void Vgrid_dtor(Vgrid **thee);

/** @brief   FORTRAN stub object destructor
 *  @ingroup Vgrid
 *  @author  Nathan Baker
 *  @param   thee   Pointer to object to be destroyed
 */
VEXTERNC void Vgrid_dtor2(Vgrid *thee);

/** @brief   Get second derivative values at a point
 *  @ingroup Vgrid
 *  @author  Steve Bond and Nathan Baker
 *  @param   thee   Pointer to Vgrid object
 *  @param   pt     Location to evaluate second derivative
 *  @param   cflag  
 *             \li  0:  Reduced Maximal Curvature
 *             \li  1:  Mean Curvature (Laplace)
 *             \li  2:  Gauss Curvature
 *             \li  3:  True Maximal Curvature
 *  @param   curv Specified curvature value
 *  @return  1 if successful, 0 if off grid
 */
VEXTERNC int Vgrid_curvature(Vgrid *thee, double pt[3], int cflag, 
  double *curv);

/** @brief   Get first derivative values at a point
 *  @ingroup Vgrid
 *  @author  Nathan Baker and Steve Bond
 *  @param   thee   Pointer to Vgrid object
 *  @param   pt     Location to evaluate gradient
 *  @param   grad   Gradient
 *  @return  1 if successful, 0 if off grid
 */
VEXTERNC int Vgrid_gradient(Vgrid *thee, double pt[3], double grad[3] );

/** @brief Write out the data in UHBD grid format 
 *  @note   \li The mesh spacing should be uniform
 *          \li Format changed from %12.6E to %12.5E
 * @ingroup Vgrid
 * @author  Nathan Baker
 * @param   thee   Grid object
 * @param   iodev  Output device type (FILE/BUFF/UNIX/INET)
 * @param   iofmt  Output device format (ASCII/XDR)
 * @param   thost  Output hostname (for sockets)
 * @param   fname  Output FILE/BUFF/UNIX/INET name
 * @param   title  Title to be inserted in grid file
 * @param   pvec   Partition information (1=>point in current partition, 
 *                 0=>point not in current partition)
 * @bug     This routine does not respect partition information
 */
VEXTERNC void Vgrid_writeUHBD(Vgrid *thee, const char *iodev, 
  const char *iofmt, const char *thost, const char *fname, char *title, 
  int *pvec);

/** @brief  Write out the data in OpenDX grid format 
 * @ingroup Vgrid
 * @author  Nathan Baker
 * @param   thee   Grid object
 * @param   iodev  Output device type (FILE/BUFF/UNIX/INET)
 * @param   iofmt  Output device format (ASCII/XDR)
 * @param   thost  Output hostname (for sockets)
 * @param   fname  Output FILE/BUFF/UNIX/INET name
 * @param   title  Title to be inserted in grid file
 * @param   pvec   Partition information (1=>point in current partition,
 *                 0=>point not in current partition)
 * @bug     This routine does not respect partition information
 */
VEXTERNC void Vgrid_writeDX(Vgrid *thee, const char *iodev, 
  const char *iofmt,  const char *thost, const char *fname, char *title, 
  int *pvec);

/** @brief   Read in data in OpenDX grid format
 *  @note    All dimension information is given in order: z, y, x
 *  @ingroup Vgrid
 *  @author  Nathan Baker
 *  @param   thee   Vgrid object
 *  @param   iodev  Input device type (FILE/BUFF/UNIX/INET)
 *  @param   iofmt  Input device format (ASCII/XDR)
 *  @param   thost  Input hostname (for sockets)
 *  @param   fname  Input FILE/BUFF/UNIX/INET name
 *  @returns 1 if sucessful, 0 otherwise
 */
VEXTERNC int Vgrid_readDX(Vgrid *thee, const char *iodev, const char *iofmt,
  const char *thost, const char *fname);

#endif
