/**
 *  @file    vunit.h
 *  @ingroup Vunit
 *  @author  Nathan Baker
 *  @brief   Contains a collection of useful constants and conversion factors
 *  @author  Nathan A. Baker
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

/** @defgroup Vunit Vunit class
 *  @brief    Collection of constants and conversion factors
 */

#ifndef _VUNIT_H_
#define _VUNIT_H_

/** @brief   Multiply by this to convert J to cal 
 *  @ingroup Vunit */
#define Vunit_J_to_cal		4.1840000e+00

/** @brief   Multiply by this to convert cal to J
 *  @ingroup Vunit */
#define Vunit_cal_to_J		2.3900574e-01

/** @brief   Multiply by this to convert amu to kg
 *  @ingroup Vunit */
#define Vunit_amu_to_kg 	1.6605402e-27

/** @brief   Multiply by this to convert kg to amu
 *  @ingroup Vunit */
#define Vunit_kg_to_amu 	6.0221367e+26

/** @brief   Multiply by this to convert ec to C
 *  @ingroup Vunit */
#define Vunit_ec_to_C		1.6021773e-19

/** @brief   Multiply by this to convert C to ec
 *  @ingroup Vunit */
#define Vunit_C_to_ec		6.2415065e+18

/** @brief   Charge of an electron in C
 *  @ingroup Vunit */
#define Vunit_ec		1.6021773e-19

/** @brief   Boltzmann constant
 *  @ingroup Vunit */
#define Vunit_kb		1.3806581e-23

/** @brief   Avogadro's number
 *  @ingroup Vunit */
#define Vunit_Na		6.0221367e+23

/** @brief   Pi
 *  @ingroup Vunit */
#define Vunit_pi		VPI

/** @brief   Vacuum permittivity
 *  @ingroup Vunit */
#define Vunit_eps0		8.8541878e-12

/** @brief \f${e_c}^2/\AA\f$ in ESU units => kcal/mol 
 *  @ingroup Vunit */
#define Vunit_esu_ec2A		3.3206364e+02

/** @brief \f$k_b\f$ in ESU units => kcal/mol 
 *  @ingroup Vunit */
#define Vunit_esu_kb            1.9871913e-03

#endif /* ifndef _VUNIT_H_ */
