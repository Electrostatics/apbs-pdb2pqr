/** @defgroup Vhal Vhal class
 *  @brief    A "class" which consists solely of macro definitions which are
 *            used by several other classes
 */

/**
 *  @file       vhal.h
 *  @ingroup    Vhal
 *  @brief      Contains generic macro definitions for APBS
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


#ifndef _VAPBSHAL_H_
#define _VAPBSHAL_H_

/**
 *  @ingroup Vhal
 *  @author  Nathan Baker
 *  @brief   Types of (scalar) data that can be written out of APBS
 */
enum Vdata_Type {
    VDT_CHARGE, /**< Charge distribution (e) */
    VDT_POT,    /**< Potential (kT/e) */
    VDT_SMOL,   /**< Solvent accessibility defined by molecular/Connolly
                 * surface definition (1 = accessible, 0 = inaccessible) */
    VDT_SSPL,   /**< Spline-based solvent accessibility (1 = accessible, 0 =
                 * inaccessible) */
    VDT_VDW,    /**< van der Waals-based accessibility (1 = accessible, 0 =
                 * inaccessible) */
    VDT_IVDW,   /**< Ion accessibility/inflated van der Waals (1 =
                 * accessible, 0 = inaccessible) */
    VDT_LAP,    /**< Laplacian of potential (kT/e/A^2) */
    VDT_EDENS,  /**< Energy density \f$\epsilon (\nabla u)^2\f$, where \f$u\f$
                 * is potential (kT/e/A)^2 */
    VDT_NDENS,  /**< Ion number density \f$\sum c_i \exp (-q_i u)^2\f$, 
	         * where $u$ is potential (output in M) */
    VDT_QDENS,  /**< Ion charge density $\sum q_i c_i \exp (-q_i u)^2$,
	         * where $u$ is potential (output in \f$e_c M\f$) */
    VDT_DIELX,  /**< Dielectric x-shifted map as calculated with the currently
                 * specified scheme (dimensionless) */
    VDT_DIELY,  /**< Dielectric y-shifted map as calculated with the currently
                 * specified scheme (dimensionless) */
    VDT_DIELZ,  /**< Dielectric y-shifted map as calculated with the currently
                 * specified scheme (dimensionless) */
    VDT_KAPPA   /**< Kappa map as calculated with the currently
                 * specified scheme (\f$A^{-3}\f$) */
};

/** @typedef Vdata_Type
 *  @ingroup Vhal
 *  @brief   Declaration of the Vdata_Type type as the Vdata_Type enum
 */
typedef enum Vdata_Type Vdata_Type;

/**
 *  @ingroup Vhal
 *  @author  Nathan Baker
 *  @brief   Format of data that can be written out of APBS
 */
enum Vdata_Format {
    VDF_DX,   /**< OpenDX (Data Explorer) format */
    VDF_UHBD, /**< UHBD format */
    VDF_AVS   /**< AVS UCD format */
};

/** @typedef Vdata_Format
 *  @ingroup Vhal
 *  @brief   Declaration of the Vdata_Format type as the Vdata_Format enum
 */
typedef enum Vdata_Format Vdata_Format;


/** @brief The maximum number of molecules that can be involved in a single 
 *         PBE calculation
 *  @ingroup Vhal 
 */
#define MAXMOL 5

/** @brief The maximum number of ion species that can be involved in a single 
 *         PBE calculation
 *  @ingroup Vhal
 */
#define MAXION 10

/** @brief The maximum number of times an MG calculation can be focused
 *  @ingroup Vhal
 */
#define MAXFOCUS 5

/** @brief   Minimum number of levels in a multigrid calculations
 *  @ingroup Vhal
 */
#define VMGNLEV 4

/** @brief   Maximum reduction of grid spacing during a focusing calculation 
 *  @ingroup Vhal
 */ 
#define VREDFRAC 0.25

/** @brief   Face definition for a volume
 *  @note    Consistent with PMG if RIGHT = EAST, BACK = SOUTH 
 *  @ingroup Vhal
 */
#define VAPBS_RIGHT 0

/** @brief   Face definition for a volume
 *  @note    Consistent with PMG if RIGHT = EAST, BACK = SOUTH 
 *  @ingroup Vhal
 */
#define VAPBS_FRONT 1

/** @brief   Face definition for a volume
 *  @note    Consistent with PMG if RIGHT = EAST, BACK = SOUTH 
 *  @ingroup Vhal
 */
#define VAPBS_UP    2

/** @brief   Face definition for a volume
 *  @note    Consistent with PMG if RIGHT = EAST, BACK = SOUTH 
 *  @ingroup Vhal
 */
#define VAPBS_LEFT  3

/** @brief   Face definition for a volume
 *  @note    Consistent with PMG if RIGHT = EAST, BACK = SOUTH 
 *  @ingroup Vhal
 */
#define VAPBS_BACK  4

/** @brief   Face definition for a volume
 *  @note    Consistent with PMG if RIGHT = EAST, BACK = SOUTH 
 *  @ingroup Vhal
 */
#define VAPBS_DOWN  5

#if !defined(VDEBUG)

/** @brief   Turns on inlining macros in Vacc class if defined
 *  @ingroup Vhal
 */
#   define VINLINE_VACC

/** @brief   Turns on inlining macros in Vatom class if defined
 *  @ingroup Vhal
 */
#   define VINLINE_VATOM

/** @brief   Turns on inlining macros in Vcsm class if defined
 *  @ingroup Vhal
 */
#   define VINLINE_VCSM

/** @brief   Turns on inlining macros in Vpbe class if defined
 *  @ingroup Vhal
 */
#   define VINLINE_VPBE

/** @brief   Turns on inlining macros in Vpee class if defined
 *  @ingroup Vhal
 */
#   define VINLINE_VPEE

/** @brief   Turns on inlining macros in Vgreen class if defined
 *  @ingroup Vhal
 */
#   define VINLINE_VGREEN

/** @brief   Turns on inlining macros in Vfetk class if defined
 *  @ingroup Vhal
 */
#   define VINLINE_VFETK

/** @brief   Turns on inlining macros in Vpmg class if defined
 *  @ingroup Vhal
 */
#   define VINLINE_VPMG

#endif


#endif
