/* ///////////////////////////////////////////////////////////////////////////
/// APBS -- Adaptive Poisson-Boltzmann Solver
///
///  Nathan A. Baker (nbaker@wasabi.ucsd.edu)
///  Dept. of Chemistry and Biochemistry
///  Dept. of Mathematics, Scientific Computing Group
///  University of California, San Diego 
///
///  Additional contributing authors listed in the code documentation.
///
/// Copyright © 1999. The Regents of the University of California (Regents).
/// All Rights Reserved. 
/// 
/// Permission to use, copy, modify, and distribute this software and its
/// documentation for educational, research, and not-for-profit purposes,
/// without fee and without a signed licensing agreement, is hereby granted,
/// provided that the above copyright notice, this paragraph and the
/// following two paragraphs appear in all copies, modifications, and
/// distributions.
/// 
/// IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
/// SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
/// ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
/// REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  
/// 
/// REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
/// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
/// PARTICULAR PURPOSE.  THE SOFTWARE AND PMGOMPANYING DOCUMENTATION, IF
/// ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS".  REGENTS HAS NO OBLIGATION
/// TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
/// MODIFICATIONS. 
////////////////////////////////////////////////////////////////////////////
/// rcsid="$Id$"
//////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// File:     vcap.c 
//
// Purpose:  Capped exponential and hyperbolic functions
//     
// Author:   Nathan Baker (based on FORTRAN code by Mike Holst)
/////////////////////////////////////////////////////////////////////////// */

#include "apbscfg.h"
#include "apbs/vcap.h"

#define EXPMAX  85.00
#define EXPMIN -85.00

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vcap_exp
//
// Purpose:  Provide a capped exp() function
//
// Notes: From Mike's PMG FORTRAN CODE:
//
//   c* to control overflow in the hyperbolic and exp functions, note
//   c* that the following are the argument limits of the various 
//   c* functions on various machines after which overflow occurs:
//   c* 
//   c*  Convex C240, Sun 3/60, Sun SPARC, IBM RS/6000:
//   c* 
//   c*    sinh, cosh, exp:     maximal argument (abs value) =  88.0d0
//   c*    dsinh, dcosh, dexp:  maximal argument (abs value) = 709.0d0
//
//   If the argument of Vcap_exp() exceeds EXPMAX or EXPMIN, then we replace it
//   with exp(EXPMAX or EXPMIN).
//
// Arguments:
//   x = argument
//   ichop = set to 1 if EXP is chopped, zero otherwise
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
// Purpose:  Provide a capped exp() function
//
// Notes: From Mike's PMG FORTRAN CODE:
//
//   c* to control overflow in the hyperbolic and exp functions, note
//   c* that the following are the argument limits of the various 
//   c* functions on various machines after which overflow occurs:
//   c* 
//   c*  Convex C240, Sun 3/60, Sun SPARC, IBM RS/6000:
//   c* 
//   c*    sinh, cosh, exp:     maximal argument (abs value) =  88.0d0
//   c*    dsinh, dcosh, dexp:  maximal argument (abs value) = 709.0d0
//
//   If the argument of Vcap_sinh() exceeds EXPMAX or EXPMIN, then we replace it
//   with sinh(EXPMAX or EXPMIN).
//
// Arguments:
//   x = argument
//   ichop = set to 1 if SINH is chopped, zero otherwise
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
// Purpose:  Provide a capped cosh() function
//
// Notes: From Mike's PMG FORTRAN CODE:
//
//   c* to control overflow in the hyperbolic and exp functions, note
//   c* that the following are the argument limits of the various 
//   c* functions on various machines after which overflow occurs:
//   c* 
//   c*  Convex C240, Sun 3/60, Sun SPARC, IBM RS/6000:
//   c* 
//   c*    sinh, cosh, exp:     maximal argument (abs value) =  88.0d0
//   c*    dsinh, dcosh, dexp:  maximal argument (abs value) = 709.0d0
//
//   If the argument of Vcap_cosh() exceeds EXPMAX or EXPMIN, then we replace it
//   with cosh(EXPMAX or EXPMIN).
//
// Arguments:
//   x = argument
//   ichop = set to 1 if COSH is chopped, zero otherwise
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
    
