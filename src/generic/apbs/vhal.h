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
 * Copyright (c) 1999-2002. The Regents of the University of California 
 *                          (Regents).  All Rights Reserved.
 *
 * Permission to use, copy, modify, and distribute this software and its
 * documentation for educational, research, and not-for-profit purposes,
 * without fee and without a signed licensing agreement, is hereby granted,
 * provided that the above copyright notice, this paragraph and the
 * following two paragraphs appear in all copies, modifications, and
 * distributions.
 *
 * IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
 * REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE.  THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF
 * ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS".  REGENTS HAS NO OBLIGATION
 * TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 * MODIFICATIONS.
 *
 * @endverbatim
 */


#ifndef _VAPBSHAL_H_
#define _VAPBSHAL_H_

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
