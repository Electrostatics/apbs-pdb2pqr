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
/// PARTICULAR PURPOSE.  THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF
/// ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS".  REGENTS HAS NO OBLIGATION
/// TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
/// MODIFICATIONS. 
//////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// File:     vunit.h    < vunit.c >
//
// Purpose:  Class Vunit:  
// A unit conversion class.
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */

#ifndef _VUNIT_H_
#define _VUNIT_H_

/* ///////////////////////////////////////////////////////////////////////////
// Class Vunit: Conversion macros
/////////////////////////////////////////////////////////////////////////// */
#define Vunit_J_to_cal		4.1840000e+00
#define Vunit_cal_to_J		2.3900574e-01
#define Vunit_amu_to_kg 	1.6605402e-27
#define Vunit_kg_to_amu 	6.0221367e+26
#define Vunit_ec_to_C		1.6021773e-19
#define Vunit_C_to_ec		6.2415065e+18

/* ///////////////////////////////////////////////////////////////////////////
// Class Vunit: Constant macros
/////////////////////////////////////////////////////////////////////////// */
#define Vunit_ec		1.6021773e-19
#define Vunit_kb		1.3806581e-23
#define Vunit_Na		6.0221367e+23
#define Vunit_eps0		8.8541878e-12

#endif /* ifndef _VUNIT_H_ */
