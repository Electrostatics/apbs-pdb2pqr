/** @defgroup Vacc Vacc class
 *  @brief    Solvent- and ion-accessibility oracle
 */

/**
 *  @file     vacc.h
 *  @ingroup  Vacc
 *  @brief    Contains declarations for class Vacc
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

#ifndef _VACC_H_
#define _VACC_H_

/* Generic headers */
#include "maloc/maloc.h"
#include "apbs/vhal.h"

/* Headers specific to this file */
#include "apbs/valist.h"
#include "apbs/vatom.h"
#include "apbs/vunit.h"

/** @brief Maximum number of neighbors used in an accessibility calculation 
 *  @ingroup Vacc
 */
#define VACCMAXNBOR 20

/**
 *  @ingroup Vacc
 *  @author  Nathan Baker
 *  @brief   Oracle for solvent- and ion-accessibility around a biomolecule
 */
struct sVacc {

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

/** 
 *  @ingroup Vacc
 *  @brief   Declaration of the Vacc class as the Vacc structure
 */
typedef struct sVacc Vacc;

#if !defined(VINLINE_VACC)

    /** @brief   Get number of bytes in this object and its members
     *  @ingroup Vacc
     *  @author  Nathan Baker
     *  @returns Number of bytes allocated for object
     */
    VEXTERNC unsigned long int Vacc_memChk(
            Vacc *thee /** Object for memory check */
            );

#else /* if defined(VINLINE_VACC) */

#   define Vacc_memChk(thee) (Vmem_bytes((thee)->vmem))

#endif /* if !defined(VINLINE_VACC) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vacc: Non-Inlineable methods (vacc.c)
/////////////////////////////////////////////////////////////////////////// */

/** @brief   Construct the accessibility object
 *  @ingroup Vacc
 *  @author  Nathan Baker
 *  @returns Newly allocated Vacc object */
VEXTERNC Vacc* Vacc_ctor(
        Valist *alist, /** Molecule for accessibility queries */
        double max_radius, /** Max probe radius (&Aring;) to be queried */ 
        int nx, /** Number of x-points in hash table */ 
        int ny, /** Number of y-points in hash table */
        int nz, /** Number of y-points in hash table */
        int nsphere /** Number of points on surface of reference sphere */
        );

/** @brief   FORTRAN stub to construct the accessibility object
 *  @ingroup Vacc
 *  @author  Nathan Baker
 *  @returns 1 if successful, 0 otherwise */
VEXTERNC int Vacc_ctor2(
        Vacc *thee, /** Memory for Vacc objet */
        Valist *alist, /** Molecule for accessibility queries */
        double max_radius, /** Max probe radius (&Aring;) to be queried */ 
        int nx, /** Number of x-points in hash table */ 
        int ny, /** Number of y-points in hash table */
        int nz, /** Number of y-points in hash table */
        int nsphere /** Number of points on surface of reference sphere */
        );

/** @brief Constructor for an accessibility object for a subdomain of a
 * molecule
 * @ingroup Vacc
 * @author  Nathan Baker, Todd Dolinsky
 * @returns Newly allocated Vacc object */
VEXTERNC Vacc* Vacc_ctorFocus(
        Valist *alist, /** Molecule for accessibility queries */
        double max_radius, /** Max probe radius (&Aring;) to be queried */ 
        int nx, /** Number of x-points in hash table */ 
        int ny, /** Number of y-points in hash table */
        int nz, /** Number of y-points in hash table */
        int nsphere, /** Number of points on surface of reference sphere */
        double x_min, /** Min x-coord of subdomain */
        double y_min, /** Min y-coord of subdomain */ 
        double z_min, /** Min z-coord of subdomain */
        double x_max, /** Max x-coord of subdomain */
        double y_max, /** Max y-coord of subdomain */
        double z_max /** Max z-coord of subdomain */
        );

/** @brief   FORTRAN stub for constructor for an accessibility object for a
 * subdomain of a molecule
 *  @ingroup Vacc
 *  @author  Nathan Baker, Todd Dolinsky
 *  @returns 1 if successful, 0 otherwise */
VEXTERNC int Vacc_ctor2Focus(
        Vacc *thee, /** Previously-allocated Vacc memory */
        Valist *alist, /** Molecule for accessibility queries */
        double max_radius, /** Max probe radius (&Aring;) to be queried */ 
        int nx, /** Number of x-points in hash table */ 
        int ny, /** Number of y-points in hash table */
        int nz, /** Number of y-points in hash table */
        int nsphere, /** Number of points on surface of reference sphere */
        double x_min, /** Min x-coord of subdomain */
        double y_min, /** Min y-coord of subdomain */ 
        double z_min, /** Min z-coord of subdomain */
        double x_max, /** Max x-coord of subdomain */
        double y_max, /** Max y-coord of subdomain */
        double z_max /** Max z-coord of subdomain */
        );

/** @brief   Destroy object
 *  @ingroup Vacc
 *  @author  Nathan Baker
 */
VEXTERNC void Vacc_dtor(
        Vacc **thee /** Pointer to memory location of object */
        );

/** @brief   FORTRAN stub to destroy object
 *  @ingroup Vacc
 *  @author  Nathan Baker
 */
VEXTERNC void Vacc_dtor2(
        Vacc *thee /** Pointer to object */
        );

/** @brief   Report van der Waals accessibility
 *  
 *  Determines if a point is within the union of the atomic spheres (with
 *  radii equal to their van der Waals radii).  
 * 
 *  @ingroup Vacc
 *  @author  Nathan Baker
 *  @returns Characteristic function value between 1.0 (accessible) and 0.0
 *          (inaccessible)
 */
VEXTERNC double Vacc_vdwAcc(
        Vacc *thee,  /** Accessibility object */
        double center[3] /** Probe center coordinates */
        );

/** @brief   Report inflated van der Waals accessibility
 *
 *  Determines if a point is within the union of the spheres centered at the
 *  atomic centers with radii equal to the sum of the atomic van der Waals
 *  radius and the probe radius. 
 *
 *  @ingroup Vacc
 *  @author  Nathan Baker
 *  @returns Characteristic function value between 1.0 (accessible) and 0.0
 *          (inaccessible)
 */
VEXTERNC double Vacc_ivdwAcc(
        Vacc *thee, /** Accessibility object */
        double center[3], /** Probe center coordinates */
        double radius /** Probe radius (&Aring;) */
        );

/** @brief   Report molecular accessibility
 *
 * Determine accessibility of a probe (of radius radius) at a given point,
 * given a collection of atomic spheres.  Uses molecular (Connolly) surface
 * definition.
 *
 *  @ingroup Vacc
 *  @author  Nathan Baker
 *  @returns Characteristic function value between 1.0 (accessible) and 0.0
 *          (inaccessible)
 */
VEXTERNC double Vacc_molAcc(
        Vacc *thee, /** Accessibility object */
        double center[3], /** Probe center coordinates */
        double radius /** Probe radius (in &Aring;) */
        );

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
 *  @returns Characteristic function value between 1.0 (accessible) and 0.0
 *          (inaccessible)
 */
VEXTERNC double Vacc_fastMolAcc(
        Vacc *thee,  /** Accessibility object */
        double center[3],  /** Probe center coordinates */
        double radius /** Probe radius (in &Aring;) */
        );

/** @brief   Report spline-based accessibility
 *
 *  Determine accessibility at a given point, given a collection of atomic
 *  spheres.  Uses Benoit Roux (Im et al, Comp Phys Comm, 111, 59--75, 1998)
 *  definition suitable for force evalation; basically a cubic spline.  
 *
 *  @ingroup Vacc
 *  @author  Nathan Baker
 *  @returns Characteristic function value between 1.0 (accessible) and 0.0
 *          (inaccessible)
 */
VEXTERNC double Vacc_splineAcc(
        Vacc *thee, /** Accessibility object */
        double center[3], /** Probe center coordinates */
        double win, /** Spline window (&Aring;) */ 
        double infrad  /** Inflation radius (&Aring;) for ion access. */
        );

/** @brief   Report spline-based accessibility for a given atom
 *
 *  Determine accessibility at a given point for a given atomic
 *  spheres.  Uses Benoit Roux (Im et al, Comp Phys Comm, 111, 59--75, 1998)
 *  definition suitable for force evalation; basically a cubic spline.
 *
 *  @ingroup Vacc
 *  @author  Nathan Baker
 *  @returns Characteristic function value between 1.0 (accessible) and 0.0
 *          (inaccessible)
 */
VEXTERNC double Vacc_splineAccAtom(
        Vacc *thee, /** Accessibility object */
        double center[3], /** Probe center coordinates */
        double win, /** Spline window (&Aring;) */ 
        double infrad, /** Inflation radius (&Aring;) for ion access. */
        int atomID /** Index of atom in Valist object */
        );

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
 */
VEXTERNC void Vacc_splineAccGradAtom(
        Vacc *thee, /** Accessibility object */
        double center[3], /** Probe center coordinates */
        double win, /** Spline window (&Aring;) */ 
        double infrad, /** Inflation radius (&Aring;) for ion access. */
        int atomID, /** Index of atom in Valist object */ 
        double *force /** 3-vector set to gradient of accessibility */
        );

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
 *  @return  Pointer to array of reference sphere points' coordinates
 */ 
VEXTERNC double** Vacc_sphere(
        Vacc *thee,  /** Accessibility object */
        int *npts /** Input:  requested number of points on sphere; Output:
                    set to actual number of points placed on sphere */
        );

/** @brief   Calculates the solvent-accessible area of the entire molecules
 *  @ingroup Vacc
 *  @note    Shamelessly ripped off from UHBD FORTRAN routine by Brock Luty
 *  @author  Nathan Baker (original FORTRAN routine by Brock Luty)
 *  @return  Solvent accessible area
 */
VEXTERNC double Vacc_totalSASA(
        Vacc *thee,  /** Accessibility object */
        double radius  /** Probe molecule radius (&Aring;) */
        );

/** @brief   Calculates the solvent-accessible area due to the specified atom
 *  @ingroup Vacc
 *  @note    Shamelessly ripped off from UHBD FORTRAN routine by Brock Luty
 *  @author  Nathan Baker (original FORTRAN routine by Brock Luty)
 *  @return  Solvent accessible area
 */
VEXTERNC double Vacc_atomSASA(
        Vacc *thee, /** Acessibility object */
        double radius, /** Probe radius (&Aring;) */
        int iatom /** Atom index in Valist object */
        );

#endif    /* ifndef _VACC_H_ */
