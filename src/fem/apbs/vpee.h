/** @defgroup Vpee Vpee class
 *  @brief  This class provides some functionality for error esimation
 *          in parallel. 
 * 
 *    This class provides some functionality for error esimation in parallel.
 *    The purpose is to modulate the error returned by some external error
 *    estimator according to the partitioning of the mesh.  For example, the
 *    Bank/Holst parallel refinement routine essentially reduces the error
 *    outside the ``local" partition to zero.  However,  this leads to the need
 *    for a few final overlapping Schwarz solves to smooth out the errors near
 *    partition boundaries.  Supposedly, if the region in which we allow
 *    error-based refinement includes the ``local" partition and an external
 *    buffer zone approximately equal in size to the local region, then the
 *    solution will asymptotically approach the solution obtained via more
 *    typical methods.  This is essentially a more flexible parallel
 *    implementation of MC's AM_markRefine.
 */

/**
 *  @file     vpee.h
 *  @ingroup  Vpee
 *  @brief    Contains declarations for class Vpee
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

#ifndef _VPEE_H
#define _VPEE_H

/* Generic headers */
#include "maloc/maloc.h"
#include "mc/mc.h"

/**
 *  @ingroup Vpee
 *  @author  Nathan Baker
 *  @brief   Contains public data members for Vpee class/module
 */
struct sVpee {

  Gem *gm;                     /**< Grid manager */
  int localPartID;             /**< The local partition ID: i.e. the partition 
                                * whose boundary simplices we're keeping
                                * track of */
  double localPartCenter[3];   /**< The coordinates of the center of the local
                                * partition */
  double localPartRadius;      /**< The radius of the circle/sphere which
                                * circumscribes the local partition */
  int killFlag;                /**< A flag indicating the method we're using to
                                * artificially decrease the error esimate
                                * outside the local partition */
  double killParam;            /**< A parameter for the error estimate
                                * attenuation method */
  Vmem *mem;                   /**< Memory manager */

};

/** 
 *  @ingroup Vpee
 *  @brief   Declaration of the Vpee class as the Vpee structure
 */
typedef struct sVpee Vpee;

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpee Inlineable methods 
/////////////////////////////////////////////////////////////////////////// */

#if !defined(VINLINE_VPEE)
#else /* if defined(VINLINE_VPEE) */
#endif /* if !defined(VINLINE_VPEE) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpee: Non-Inlineable methods (vpee.c)
/////////////////////////////////////////////////////////////////////////// */

/** @brief   Construct the Vpee object
 *  @ingroup Vpee
 *  @author  Nathan Baker
 *  @param   gm  Gem (geometry manager) object
 *  @param   localPartID The ID of the local partition
 *  @param   killFlag  A flag to indicate how error estimates are to be
 *                     attenuated outside the local partition:
 *                     \li 0:  no attenuation
 *                     \li 1:  all error outside the local partition set to
 *                           zero
 *                     \li 2:  all error is set to zero outside a sphere of
 *                           radius (killParam*partRadius), where
 *                           partRadius is the radius of the sphere
 *                           circumscribing the local partition
 *                     \li 3:  all error is set to zero except for the local
 *                           partition and its immediate neighbors
 * @param    killParam @see killFlag for usage
 * @return   Newly constructed Vpee object
 */
VEXTERNC Vpee* Vpee_ctor(Gem *gm, int localPartID, int killFlag, 
                 double killParam);

/** @brief   FORTRAN stub to construct the Vpee object
 *  @ingroup Vpee
 *  @author  Nathan Baker
 *  @param   thee  Memory location for new object
 *  @param   gm  Gem (geometry manager) object
 *  @param   localPartID The ID of the local partition
 *  @param   killFlag  A flag to indicate how error estimates are to be
 *                     attenuated outside the local partition:
 *                     \li 0:  no attenuation
 *                     \li 1:  all error outside the local partition set to
 *                           zero
 *                     \li 2:  all error is set to zero outside a sphere of
 *                           radius (killParam*partRadius), where
 *                           partRadius is the radius of the sphere
 *                           circumscribing the local partition
 *                     \li 3:  all error is set to zero except for the local
 *                           partition and its immediate neighbors
 * @param    killParam @see killFlag for usage
 * @return   1 if successful, 0 otherwise
 */
VEXTERNC int   Vpee_ctor2(Vpee *thee, Gem *gm, int localPartID, int killFlag,
                 double killParam);

/** @brief   Object destructor
 *  @ingroup Vpee
 *  @author  Nathan Baker
 *  @param   thee   Pointer to memory location of object to be destroyed
 */
VEXTERNC void  Vpee_dtor(Vpee **thee);

/** @brief   FORTRAN stub object destructor
 *  @ingroup Vpee
 *  @author  Nathan Baker
 *  @param   thee   Pointer to object to be destroyed
 */
VEXTERNC void  Vpee_dtor2(Vpee *thee);

/** @brief   Mark simplices for refinement based on attenuated error estimates.
 *  
 *  A wrapper/reimplementation of AM_markRefine that allows for more flexible
 *  attenuation of error-based markings outside the local partition.  The error
 *  in each simplex is modified by the method (see killFlag) specified in the
 *  Vpee constructor.  This allows the user to confine refinement to an
 *  arbitrary area around the local partition.
 * 
 *  @ingroup Vpee
 *  @author  Nathan Baker and Mike Holst
 *  @note  This routine borrows very heavily from FEtk routines by Mike Holst.
 *  @param thee The Vpee object
 *  @param am   The AM (algebra manager) object currently in use for solving
 *              the PBE
 *  @param level Current level of multigrid hierarchy
 *  @param akey  The marking method:
 *               \li -1:  Reset markings  --> killFlag has no effect.
 *               \li 0:  Uniform.
 *               \li 1:  User defined (geometry-based).
 *               \li >1:  A numerical estimate for the error has already been
 *                        set in am and should be attenuated according to
 *                        killFlag and used, in conjunction with etol, to mark
 *                        simplices for refinement.
 *  @param rcol The ID of the main partition on which to mark (or -1 if all
 *              partitions should be marked).  Note that we should have (rcol
 *              == thee->localPartID) for (thee->killFlag == 2 or 3)
 *  @param etol The error tolerance criterion for marking
 *  @param bkey How the error tolerance is interpreted:
 *              \li 0:  Simplex marked if error > \f$\eta\f$.
 *              \li 1:  Simplex marked if error >.
 *                \f$\left(\eta^2/L\right)^{1/2}\f$
 *         where \f$\eta\f$ is the error tolerance and \f$L\f$ is the number of
 *         simplices
 *  @return The number of simplices marked for refinement.
 */
VEXTERNC int  Vpee_markRefine(Vpee *thee, AM *am, int level, int akey, 
                 int rcol, double etol, int bkey);

/** @brief   Returns the number of simplices in the local partition
 *  @ingroup Vpee
 *  @author  Nathan Baker
 *  @param   thee Vpee object
 *  @return  Number of simplices in the local partition
 */
VEXTERNC int   Vpee_numSS(Vpee *thee);

#endif    /* ifndef _VPEE_H_ */
