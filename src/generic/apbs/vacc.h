/** @defgroup Vacc Vacc class
 *  @brief    Solvent- and ion-accessibility oracle
 */

/**
 *  @file     vacc.h
 *  @ingroup  Vacc
 *  @brief    Contains declarations for class Vacc
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
 * Copyright (c) 1999-2002.  The Regents of the University of California.
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
 *
 * @endverbatim
 */

#ifndef _VACC_H_
#define _VACC_H_

#include "maloc/maloc.h"
#include "apbs/vhal.h"
#include "apbs/valist.h"
#include "apbs/vunit.h"

/** @brief Maximum number of neighbors used in an accessibility calculation 
 *  @ingroup Vacc
 */
#define VACCMAXNBOR 20

/**
 *  @struct  Vacc
 *  @ingroup Vacc
 *  @author  Nathan Baker
 *  @brief   Oracle for solvent- and ion-accessibility around a biomolecule
 */
struct Vacc {

  Vmem *vmem;               /**< Memory management object for this class */
  Valist *alist;            /**< Valist structure for list of atoms */
  int **atomIDs;            /**< An array of arrays of atom IDs in alist */
  int *natoms;              /**< An array telling how many pointers are 
                             * stored in atoms[i] */
  int *atomFlags;           /**< Flags for keeping track of atoms */
  double **sphere;          /**< An array of points on the surface of a 
                             * sphere */
  int nsphere;              /**< The number of points in thee->sphere */
  Vset acc;                 /**< An integer array (to be treated as 
                             * bitfields) of Vset type with length equal 
                             * to the number of vertices in the mesh */
  double grid_lower_corner[3]; /**< Hash table grid corner */
  double hx;                /**< Hash table grid spacings */
  double hy;                /**< Hash table grid spacings */
  double hzed;              /**< Hash table grid spacings */
  int nx;                   /**< Hash table grid dimensions */
  int ny;                   /**< Hash table grid dimensions */
  int nz;                   /**< Hash table grid dimensions */
  int n;                    /**< n = nx*nz*ny */
  double max_radius;        /**< Maximum probe radius */
  double *area;             /**< The contribution to the solvent-accessible
                             * surface area from each atom.  This array is
                             * filled by Vacc_totalSASA */
};

/** @typedef Vacc
 *  @ingroup Vacc
 *  @brief   Declaration of the Vacc class as the Vacc structure
 */
typedef struct Vacc Vacc;

/* ///////////////////////////////////////////////////////////////////////////
// Class Vacc: Inlineable methods (vacc.c)
/////////////////////////////////////////////////////////////////////////// */

#if !defined(VINLINE_VACC)

    /** @brief   Get number of bytes in this object and its members
     *  @ingroup Vacc
     *  @author  Nathan Baker
     *  @param   thee  Vacc object
     *  @returns Number of bytes allocated for object
     */
    VEXTERNC int Vacc_memChk(Vacc *thee);

#else /* if defined(VINLINE_VACC) */

#   define Vacc_memChk(thee) (Vmem_bytes((thee)->vmem))

#endif /* if !defined(VINLINE_VACC) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vacc: Non-Inlineable methods (vacc.c)
/////////////////////////////////////////////////////////////////////////// */

/** @brief   Construct the accessibility object
 *  @ingroup Vacc
 *  @author  Nathan Baker
 *  @param   alist      List of atom objects
 *  @param   max_radius Maximum probe radius to be queried
 *  @param   nx         Number of grid points in hash table (x)
 *  @param   ny         Number of grid points in hash table (y)
 *  @param   nz         Number of grid points in hash table (z)
 *  @param   nsphere    Number of points on surface of reference sphere
 *  @returns Newly allocated Vacc object */
VEXTERNC Vacc* Vacc_ctor(Valist *alist, double max_radius, int nx, 
               int ny, int nz, int nsphere);

/** @brief   FORTRAN stub to construct the accessibility object
 *  @ingroup Vacc
 *  @author  Nathan Baker
 *  @param   thee       Allocated Vacc memory
 *  @param   alist      List of atom objects
 *  @param   max_radius Maximum probe radius to be queried
 *  @param   nx         Number of grid points in hash table (x)
 *  @param   ny         Number of grid points in hash table (y)
 *  @param   nz         Number of grid points in hash table (z)
 *  @param   nsphere    Number of points on surface of reference sphere
 *  @returns 1 if successful, 0 otherwise */
VEXTERNC int Vacc_ctor2(Vacc *thee, Valist *alist, double max_radius, 
             int nx, int ny, int nz, int nsphere);

/** @brief   Destroy object
 *  @ingroup Vacc
 *  @author  Nathan Baker
 *  @param   thee  Pointer to memory location of object */
VEXTERNC void Vacc_dtor(Vacc **thee);

/** @brief   FORTRAN stub to destroy object
 *  @ingroup Vacc
 *  @author  Nathan Baker
 *  @param   thee  Pointer to object */
VEXTERNC void Vacc_dtor2(Vacc *thee);

/** @brief   Report van der Waals accessibility
 *  
 *  Determines if a point is within the union of the atomic spheres (with
 *  radii equal to their van der Waals radii).  
 * 
 *  @ingroup Vacc
 *  @author  Nathan Baker
 *  @param   thee   Vacc object
 *  @param   center Probe center coordinates (point to test)
 *  @returns Characteristic function value between 1.0 (accessible) and 0.0
 *          (inaccessible)
 */
VEXTERNC double Vacc_vdwAcc(Vacc *thee, double center[3]);

/** @brief   Report inflated van der Waals accessibility
 *
 *  Determines if a point is within the union of the spheres centered at the
 *  atomic centers with radii equal to the sum of the atomic van der Waals
 *  radius and the probe radius. 
 *
 *  @ingroup Vacc
 *  @author  Nathan Baker
 *  @param   thee   Vacc object
 *  @param   center Probe center coordinates (point to test)
 *  @param   radius Probe radius
 *  @returns Characteristic function value between 1.0 (accessible) and 0.0
 *          (inaccessible)
 */
VEXTERNC double Vacc_ivdwAcc(Vacc *thee, double center[3], double radius);

/** @brief   Report molecular accessibility
 *
 * Determine accessibility of a probe (of radius radius) at a given point,
 * given a collection of atomic spheres.  Uses molecular (Connolly) surface
 * definition.
 *
 *  @ingroup Vacc
 *  @author  Nathan Baker
 *  @param   thee   Vacc object
 *  @param   center Probe center coordinates (point to test)
 *  @param   radius Probe radius
 *  @returns Characteristic function value between 1.0 (accessible) and 0.0
 *          (inaccessible)
 */
VEXTERNC double Vacc_molAcc(Vacc *thee, double center[3], double radius);

/** @brief   Report molecular accessibility quickly
 *
 *  Given a point which is INSIDE the collection of inflated van der Waals
 *  spheres, but OUTSIDE the collection of non-inflated van der Waals spheres,
 *  determine accessibility of a probe (of radius radius) at a given point,
 *  given a collection of atomic spheres.  Uses molecular (Connolly) surface
 *  definition.  
 *
 *  @note    THIS ASSUMES YOU HAVE TESTED THAT THIS POINT IS DEFINITELY INSIDE
 *           THE INFLATED AND NON-INFLATED VAN DER WAALS SURFACES!
 *  @ingroup Vacc
 *  @author  Nathan Baker
 *  @param   thee   Vacc object
 *  @param   center Probe center coordinates (point to test)
 *  @param   radius Probe radius
 *  @returns Characteristic function value between 1.0 (accessible) and 0.0
 *          (inaccessible)
 */
VEXTERNC double Vacc_fastMolAcc(Vacc *thee, double center[3], double radius);

/** @brief   Report spline-based accessibility
 *
 *  Determine accessibility at a given point, given a collection of atomic
 *  spheres.  Uses Benoit Roux (Im et al, Comp Phys Comm, 111, 59--75, 1998)
 *  definition suitable for force evalation; basically a cubic spline.  
 *
 *  @ingroup Vacc
 *  @author  Nathan Baker
 *  @param   thee   Vacc object
 *  @param   center Probe center coordinates (point to test)
 *  @param   win    Spline window
 *  @param   infrad Inflation radius (for ion-accessibility) 
 *  @returns Characteristic function value between 1.0 (accessible) and 0.0
 *          (inaccessible)
 */
VEXTERNC double Vacc_splineAcc(Vacc *thee, double center[3], double window,
  double infrad);

/** @brief   Report spline-based accessibility for a given atom
 *
 *  Determine accessibility at a given point for a given atomic
 *  spheres.  Uses Benoit Roux (Im et al, Comp Phys Comm, 111, 59--75, 1998)
 *  definition suitable for force evalation; basically a cubic spline.
 *
 *  @ingroup Vacc
 *  @author  Nathan Baker
 *  @param   thee   Vacc object
 *  @param   center Probe center coordinates (point to test)
 *  @param   win    Spline window
 *  @param   infrad Inflation radius (for ion-accessibility)
 *  @returns Characteristic function value between 1.0 (accessible) and 0.0
 *          (inaccessible)
 */
VEXTERNC double Vacc_splineAccAtom(Vacc *thee, double center[3], double window,
  double infrad, int atomID);

/** @brief   Report gradient of spline-based accessibility with respect to a
 *           particular atom normalized by the accessibility value due to that
 *           atom at that point (see Vpmg_splineAccAtom)
 *
 *  Determine accessibility at a given point, given a collection of atomic
 *  spheres.  Uses Benoit Roux (Im et al, Comp Phys Comm, 111, 59--75, 1998)
 *  definition suitable for force evalation; basically a cubic spline. 
 *
 *  @ingroup Vacc
 *  @author  Nathan Baker
 *  @param   thee   Vacc object
 *  @param   center Probe center coordinates (point to test)
 *  @param   win    Spline window
 *  @param   infrad Inflation radius (for ion-accessibility)
 *  @param   atomID Valist ID of atom used in calculation
 *  @param   force  Stores gradient of accesibility (3-vector)
 */
VEXTERNC void Vacc_splineAccGradAtom(Vacc *thee, double center[3], double win,
  double infrad, int atomID, double *force);

/** @brief Set up an array of points for the reference sphere
 *  
 *  Generates approximately npts # of points (actual number stored in npts)
 *  somewhat uniformly distributed across a sphere of unit radius centered at
 *  the origin.  Returns an (npts x 3) double array, which the user is
 *  responsible for destroying.
 * 
 *  @note  This routine was shamelessly ripped off from sphere.f from UHBD as
 *         developed by Michael K. Gilson.
 * 
 *  @ingroup Vacc
 *  @author  Nathan Baker (original FORTRAN code by Mike Gilson)
 *  @param   thee  Vacc object
 *  @param   npts  Requested number of points on sphere, later set to actual
 *                 number of points placed on sphere
 *  @return  Pointer to array of reference sphere points' coordinates
 */ 
VEXTERNC double** Vacc_sphere(Vacc *thee, int *npts);

/** @brief   Calculates the solvent-accessible area of the entire molecules
 *  @ingroup Vacc
 *  @note    Shamelessly ripped off from UHBD FORTRAN routine by Brock Luty
 *  @author  Nathan Baker (original FORTRAN routine by Brock Luty)
 *  @param   thee   Vacc object
 *  @param   radius Radius of probe molecule used to calculate solvent
 *           accessible area
 *  @return  Solvent accessible area
 */
VEXTERNC double Vacc_totalSASA(Vacc *thee, double radius);

/** @brief   Calculates the solvent-accessible area due to the specified atom
 *  @ingroup Vacc
 *  @note    Shamelessly ripped off from UHBD FORTRAN routine by Brock Luty
 *  @author  Nathan Baker (original FORTRAN routine by Brock Luty)
 *  @param   thee   Vacc object
 *  @param   radius Radius of probe molecule used to calculate solvent
 *           accessible area
 *  @param   iatom  Atom ID in Valist list
 *  @return  Solvent accessible area
 */
VEXTERNC double Vacc_atomSASA(Vacc *thee, double radius, int iatom);

#endif    /* ifndef _VACC_H_ */
