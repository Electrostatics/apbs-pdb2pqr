/**
 *  @file    vcap.c
 *  @ingroup Vcap
 *  @author  Nathan Baker
 *  @brief   Class Vcap methods
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


#include "apbscfg.h"
#include "apbs/vcap.h"

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vcap_exp
//
// Author:  Nathan Baker (based on PMG by Mike Holst)
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vcap_exp(double x, int *ichop) {

    /* The two chopped arguments */
    if (x > EXPMAX) {
       (*ichop) = 1;
       return VEXP(EXPMAX);
    } else if (x < EXPMIN) {
       (*ichop) = 1;
       return VEXP(EXPMIN);
    } 

    /* The normal EXP */
    (*ichop) = 0;
    return VEXP(x);
}
    
/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vcap_sinh
//
// Author:  Nathan Baker (based on PMG by Mike Holst)
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vcap_sinh(double x, int *ichop) {

    /* The two chopped arguments */
    if (x > EXPMAX) {
       (*ichop) = 1;
       return VSINH(EXPMAX);
    } else if (x < EXPMIN) {
       (*ichop) = 1;
       return VSINH(EXPMIN);
    } 

    /* The normal SINH */
    (*ichop) = 0;
    return VSINH(x);
}
    
/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vcap_cosh
//
// Author:  Nathan Baker (based on PMG by Mike Holst)
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vcap_cosh(double x, int *ichop) {

    /* The two chopped arguments */
    if (x > EXPMAX) {
       (*ichop) = 1;
       return VCOSH(EXPMAX);
    } else if (x < EXPMIN) {
       (*ichop) = 1;
       return VCOSH(EXPMIN);
    } 

    /* The normal COSH */
    (*ichop) = 0;
    return VCOSH(x);
}
    
