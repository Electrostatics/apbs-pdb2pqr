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
 * Nathan A. Baker (baker@biochem.wustl.edu)
 * Dept. of Biochemistry and Molecular Biophysics
 * Center for Computational Biology
 * Washington University in St. Louis
 *
 * Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 2003.  Washington University in St. Louis.
 * All Rights Reserved.
 * Portions Copyright (c) 1999-2003.  The Regents of the University of
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
    
