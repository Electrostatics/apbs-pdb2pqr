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
////////////////////////////////////////////////////////////////////////////
/// rcsid="$Id$"
//////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// File:     vgreen.h    < vgreen.c >
//
// Purpose:  Class Vgreen:  
//   Provides capabilities for pointwise evaluation of free space Green's
//   function for point charges in a uniform dielectric.
//
// Notes:  Right now, this is a very slow method without any fast multipole 
//   method acceleration
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */

#ifndef _VGREEN_H_
#define _VGREEN_H_

#include "mc/mc.h"
#include "apbs/vatom.h"
#include "apbs/valist.h"
#include "apbs/vunit.h"

/* ///////////////////////////////////////////////////////////////////////////
// Class Vgreen: Parameters and datatypes
/////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vgreen: Definition
/////////////////////////////////////////////////////////////////////////// */

typedef struct Vgreen { 

  Valist *alist;      /* Atom (charge) list */
  Vmem *vmem;         /* Memory management object */
  int initFlagCXXFMM; /* Flag to indicate whether the C++ FMM object has been
                       * initialized */

} Vgreen;

/* ///////////////////////////////////////////////////////////////////////////
// Class Vgreen: Inlineable methods (vgreen.c)
/////////////////////////////////////////////////////////////////////////// */

#if !defined(VINLINE_VGREEN)
    VEXTERNC Valist* Vgreen_getValist(Vgreen *thee);
    VEXTERNC int     Vgreen_memChk(Vgreen *thee);
#else /* if defined(VINLINE_VGREEN) */
#   define Vgreen_getValist(thee) ((thee)->alist)
#   define Vgreen_memChk(thee) (Vmem_bytes((thee)->vmem))
#endif /* if !defined(VINLINE_VGREEN) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vgreen: Non-Inlineable methods (vgreen.c)
/////////////////////////////////////////////////////////////////////////// */

VEXTERNC Vgreen* Vgreen_ctor(Valist *alist);
VEXTERNC int     Vgreen_ctor2(Vgreen *thee, Valist *alist);
VEXTERNC void    Vgreen_dtor(Vgreen **thee);
VEXTERNC void    Vgreen_dtor2(Vgreen *thee);
VEXTERNC double  Vgreen_helmholtz(Vgreen *thee, double *position, double dim,
                   double kappa);
VEXTERNC void    Vgreen_helmholtzD(Vgreen *thee, double *position, 
                   double dim, double kappa, double *grad);
VEXTERNC double  Vgreen_coulomb(Vgreen *thee, double *position, double dim);
VEXTERNC void    Vgreen_coulombD(Vgreen *thee, double *position, double dim,
                   double *grad);

#endif /* ifndef _VGREEN_H_ */
