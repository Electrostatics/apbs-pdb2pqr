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
 * Copyright (c) 2002-2005.  Washington University in St. Louis.
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
 * Linking APBS statically or dynamically with other modules is making a
 * combined work based on APBS. Thus, the terms and conditions of the GNU
 * General Public License cover the whole combination.
 * 
 * SPECIAL GPL EXCEPTION
 * In addition, as a special exception, the copyright holders of APBS
 * give you permission to combine the APBS program with free software
 * programs and libraries that are released under the GNU LGPL or with
 * code included in releases of ISIM, PMV, PyMOL, SMOL, VMD, and Vision.
 * Such combined software may be linked with APBS and redistributed together 
 * in original or modified form as mere aggregation without requirement that 
 * the entire work be under the scope of the GNU General Public License.
 * This special exception permission is also extended to any software listed
 * in the SPECIAL GPL EXCEPTION clauses by the PMG, FEtk, MC, or MALOC
 * libraries.
 * 
 * Note that people who make modified versions of APBS are not obligated
 * to grant this special exception for their modified versions; it is
 * their choice whether to do so. The GNU General Public License gives
 * permission to release a modified version without this exception; this
 * exception also makes it possible to release a modified version which
 * carries forward this exception.
 *
 * @endverbatim
 */

#include "apbscfg.h"
#include "apbs/vcap.h"

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
    
