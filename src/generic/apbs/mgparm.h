/** @defgroup MGparm MGparm class
 *  @brief    Parameter which holds useful parameters for generic multigrid
 *            calculations
 */

/**
 *  @file     mgparm.h
 *  @ingroup  MGparm
 *  @brief    Contains declarations for class MGparm
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


#ifndef _MGPARM_H_
#define _MGPARM_H_

#include "apbs/apbs.h"
#include "maloc/maloc.h"

/**
 *  @struct  MGparm
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @brief   Parameter structure for MG-specific variables from input files
 *  @note    If you add/delete/change something in this class, the member
 *           functions -- especially MGparm_copy -- must be modified
 *           accordingly
 */
struct MGparm {

    int type;                   /**< What type of MG calculation?
                                 *   \li 0: sequential manual
                                 *   \li 1: sequential auto-focus
                                 *   \li 2: parallel auto-focus 
                                 *   \li 3: dummy calculation for coefficient
                                 *          I/O */
    int parsed;                 /**< Has this structure been filled? (0 = no,
                                 * 1 = yes) */

    /* *** GENERIC PARAMETERS *** */
    int dime[3];               /**< Grid dimensions */
    int setdime;               /**< Flag, @see dime */

    /* *** TYPE 0 PARAMETERS (SEQUENTIAL MANUAL) *** */
    int nlev;                  /**< Levels in multigrid hierarchy 
                                *   @deprecated Just ignored now */
    int setnlev;               /**< Flag, @see nlev */
    double grid[3];            /**< Grid spacings */
    int setgrid;               /**< Flag, @see grid */
    double glen[3];            /**< Grid side lengths. */
    int setglen;               /**< Flag, @see glen */
    int cmeth;                 /**< Centering method:  
                                *   \li 0: center on point, 
                                *   \li 1: center on molecule */
    double center[3];          /**< Grid center. If ispart = 0, then this is
				* only meaningful if cmeth = 0.  However, if
				* ispart = 1 and cmeth = 0, then this is the
				* center of the non-disjoint (overlapping)
				* partition.  If ispart = 1 and cmeth = 1, then
				* this is the vector that must be added to the
				* center of the molecule to give the center of
				* the non-disjoint partition.  */
    int centmol;               /**< Particular molecule on which we want to
                                * center the grid */
    int setgcent;              /**< Flag, @see cmeth */

    /* ******** TYPE 1 & 2 PARAMETERS (SEQUENTIAL & PARALLEL AUTO-FOCUS) *** */
    double cglen[3];           /**< Coarse grid side lengths */
    int setcglen;              /**< Flag, @see cglen */
    double fglen[3];           /**< Fine grid side lengths */
    int setfglen;              /**< Flag, @see fglen */
    int ccmeth;                /**< Coarse grid centering method:  0 => center
				* on point, 1 => center on molecule */
    double ccenter[3];         /**< Coarse grid center.  */
    int ccentmol;              /**< Particular molecule on which we want to
                                * center the coarse grid */
    int setcgcent;             /**< Flag, @see ccmeth */
    int fcmeth;                /**< Fine grid centering method:  0 => center on
                                * point, 1 => center on molecule */
    double fcenter[3];         /**< Fine grid center.  */
    int fcentmol;              /**< Particular molecule on which we want to
                                * center the fine grid */
    int setfgcent;             /**< Flag, @see fcmeth */


    /* ********* TYPE 2 PARAMETERS (PARALLEL AUTO-FOCUS) ******** */
    double partDisjCenterShift[3];       /**< When added to the actual (local)
                                          * mesh center, this gives the center
                                          * of the disjoint partitions */
    double partDisjLength[3];            /**< This gives the lengths of the
                                          * disjoint partitions */
    int partDisjOwnSide[6];              /**< Tells whether the boundary points
                                          * are ours (1) or not (0) */
    double partOlapCenterShift[3];       /**< When added to the actual (local)
                                          * mesh center, this gives the center
                                          * of the overlapping partitions */
    double partOlapLength[3];            /**< This gives the lengths of the
                                          * overlapping partitions */

    int pdime[3];                        /**< Grid of processors to be used in
                                          * calculation */
    int setpdime;                        /**< Flag, @see pdime */
    int proc_rank;                       /**< Rank of this processor */
    int setrank;                         /**< Flag, @see proc_rank */
    int proc_size;                       /**< Total number of processors */
    int setsize;                         /**< Flag, @see proc_size */
    double ofrac;                        /**< Overlap fraction between procs */
    int setofrac;                        /**< Flag, @see ofrac */

};

/** @typedef MGparm
 *  @ingroup MGparm
 *  @brief   Declaration of the MGparm class as the MGparm structure
 */
typedef struct MGparm MGparm;

/** @brief   Get number of grid points in x direction
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @returns Number of grid points in the x direction
 */
VEXTERNC int MGparm_getNx(MGparm *thee);

/** @brief   Get number of grid points in y direction
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @returns Number of grid points in the y direction
 */
VEXTERNC int MGparm_getNy(MGparm *thee);

/** @brief   Get number of grid points in z direction
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @returns Number of grid points in the z direction
 */
VEXTERNC int MGparm_getNz(MGparm *thee);

/** @brief   Get grid spacing in x direction (\f$\AA\f$)
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @returns Grid spacing in the x direction
 */
VEXTERNC double MGparm_getHx(MGparm *thee);

/** @brief   Get grid spacing in y direction (\f$\AA\f$)
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @returns Grid spacing in the y direction
 */
VEXTERNC double MGparm_getHy(MGparm *thee);

/** @brief   Get grid spacing in z direction (\f$\AA\f$)
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @returns Grid spacing in the z direction
 */
VEXTERNC double MGparm_getHz(MGparm *thee);

/** @brief   Set center x-coordinate
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @param   thee  x-coordinate
 */
VEXTERNC void MGparm_setCenterX(MGparm *thee, double x);

/** @brief   Set center y-coordinate
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @param   thee  y-coordinate
 */  
VEXTERNC void MGparm_setCenterY(MGparm *thee, double y);

/** @brief   Set center z-coordinate
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @param   thee  z-coordinate
 */
VEXTERNC void MGparm_setCenterZ(MGparm *thee, double z);

/** @brief   Get center x-coordinate
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @returns  x-coordinate
 */
VEXTERNC double MGparm_getCenterX(MGparm *thee);

/** @brief   Get center y-coordinate
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @returns  y-coordinate
 */
VEXTERNC double MGparm_getCenterY(MGparm *thee);

/** @brief   Get center z-coordinate
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @returns  z-coordinate
 */
VEXTERNC double MGparm_getCenterZ(MGparm *thee);

/** @brief   Get x-coordinate shift of partition center in parallel calculation
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @returns  x-coordinate shift of partition center in parallel calculation
 */
VEXTERNC double MGparm_getPartOlapCenterShiftX(MGparm *thee);

/** @brief   Get y-coordinate shift of partition center in parallel calculation
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @returns  y-coordinate
 */
VEXTERNC double MGparm_getPartOlapCenterShiftY(MGparm *thee);

/** @brief   Get z-coordinate shift of partition center in parallel calculation
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @returns  z-coordinate shift of partition center in parallel calculation
 */
VEXTERNC double MGparm_getPartOlapCenterShiftZ(MGparm *thee);

/** @brief   Construct MGparm object
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   type Type of MG calculation
 *                \li 0: sequential manual
 *                \li 1: sequential auto-focus
 *                \li 2: parallel auto-focus
 *  @returns Newly allocated and initialized MGparm object
 */
VEXTERNC MGparm*  MGparm_ctor(int type);

/** @brief   FORTRAN stub to construct MGparm object
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee Space for MGparm object
 *  @param   type Type of MG calculation
 *                \li 0: sequential manual
 *                \li 1: sequential auto-focus
 *                \li 2: parallel auto-focus
 *  @returns 1 if succesful, 0 otherwise
 */
VEXTERNC int      MGparm_ctor2(MGparm *thee, int type);

/** @brief   Object destructor
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  Pointer to memory location of MGparm object
 */
VEXTERNC void     MGparm_dtor(MGparm **thee);

/** @brief   FORTRAN stub for object destructor
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  Pointer to MGparm object
 */
VEXTERNC void     MGparm_dtor2(MGparm *thee);

/** @brief   Consistency check for parameter values stored in object
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee   MGparm object
 *  @returns 1 if OK, 0 otherwise
 */
VEXTERNC int      MGparm_check(MGparm *thee);

/** @brief   Copy MGparm object into thee
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee   MGparm object (target for copy)
 *  @param   parm   MGparm object (source for copy)
 */
VEXTERNC void     MGparm_copy(MGparm *thee, MGparm *parm);

/** @brief   Parse an MG keyword from an input file
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee   MGparm object 
 *  @param   tok    Token to parse
 *  @param   sock   Stream for more tokens
 *  @return   1 if matched and assigned; -1 if matched, but there's some sort
 *            of error (i.e., too few args); 0 if not matched
 */
VEXTERNC int      MGparm_parseToken(MGparm *thee, char tok[VMAX_BUFSIZE], 
                    Vio *sock);

#endif 

