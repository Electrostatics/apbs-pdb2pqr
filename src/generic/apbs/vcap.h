/** @defgroup Vcap Vcap class
 *  @brief  Collection of routines which cap certain exponential and hyperbolic
 *          functions
 *  @note   These routines are based on FORTRAN code by Mike Holst
 */

/**
 *  @file     vcap.h
 *  @ingroup  Vcap
 *  @brief    Contains declarations for class Vcap
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

#ifndef _VCAP_H_
#define _VCAP_H_

/** @brief   Maximum argument for exp(), sinh(), or cosh()
 *  @ingroup Vcap
 */
#define EXPMAX  85.00

/** @brief   Minimum argument for exp(), sinh(), or cosh()
 *  @ingroup Vcap
 */
#define EXPMIN -85.00

#include "maloc/maloc.h"

/* ///////////////////////////////////////////////////////////////////////////
// Class Vcap: Non-Inlineable methods (vcap.c)
/////////////////////////////////////////////////////////////////////////// */

/** @brief   Provide a capped exp() function
 *
 *           If the argument x of Vcap_exp() exceeds EXPMAX or EXPMIN, then we
 *           return exp(EXPMAX) or exp(EXPMIN) rather than exp(x).
 *
 *  @note    Original FORTRAN routine from PMG library by Mike Holst
 *           Original notes:
 *           to control overflow in the hyperbolic and exp functions, note
 *           that the following are the argument limits of the various 
 *           functions on various machines after which overflow occurs:
 *           Convex C240, Sun 3/60, Sun SPARC, IBM RS/6000:
 *           sinh, cosh, exp:     maximal argument (abs value) =  88.0d0
 *           dsinh, dcosh, dexp:  maximal argument (abs value) = 709.0d0
 *
 *  @ingroup Vcap
 *  @author  Nathan Baker (based on FORTRAN code by Mike Holst)
 *  @param   x     Argument to exp()
 *  @param   ichop set to 1 if exp function capped
 *  @return  exp(x) or capped equivalent
 */
VEXTERNC double Vcap_exp(double x, int *ichop);


/** @brief   Provide a capped sinh() function
 *
 *           If the argument x of Vcap_sinh() exceeds EXPMAX or EXPMIN, then we
 *           return sinh(EXPMAX) or sinh(EXPMIN) rather than sinh(x).
 *
 *  @note    Original FORTRAN routine from PMG library by Mike Holst
 *           Original notes:
 *           to control overflow in the hyperbolic and exp functions, note
 *           that the following are the argument limits of the various 
 *           functions on various machines after which overflow occurs:
 *           Convex C240, Sun 3/60, Sun SPARC, IBM RS/6000:
 *           sinh, cosh, exp:     maximal argument (abs value) =  88.0d0
 *           dsinh, dcosh, dexp:  maximal argument (abs value) = 709.0d0
 *
 *  @ingroup Vcap
 *  @author  Nathan Baker (based on FORTRAN code by Mike Holst)
 *  @param   x     Argument to sinh()
 *  @param   ichop set to 1 if sinh function capped
 *  @return  sinh(x) or capped equivalent
 */
VEXTERNC double Vcap_sinh(double x, int *ichop);

/** @brief   Provide a capped cosh() function
 *
 *           If the argument x of Vcap_cosh() exceeds EXPMAX or EXPMIN, then we
 *           return cosh(EXPMAX) or cosh(EXPMIN) rather than cosh(x).
 *
 *  @note    Original FORTRAN routine from PMG library by Mike Holst
 *           Original notes:
 *           to control overflow in the hyperbolic and exp functions, note
 *           that the following are the argument limits of the various
 *           functions on various machines after which overflow occurs:
 *           Convex C240, Sun 3/60, Sun SPARC, IBM RS/6000:
 *           sinh, cosh, exp:     maximal argument (abs value) =  88.0d0
 *           dsinh, dcosh, dexp:  maximal argument (abs value) = 709.0d0
 *
 *  @ingroup Vcap
 *  @author  Nathan Baker (based on FORTRAN code by Mike Holst)
 *  @param   x     Argument to cosh()
 *  @param   ichop set to 1 if cosh function capped
 *  @return  cosh(x) or capped equivalent
 */
VEXTERNC double Vcap_cosh(double x, int *ichop);

#endif    /* ifndef _VCAP_H_ */
