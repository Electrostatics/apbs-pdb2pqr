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
// File:     vfetk.h    < vfetk.c >
//
// Purpose:  
//    Class Vfetk:  
//      Provides the interface to FEtk for solution of the PBE
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */

#ifndef _VFETK_H_
#define _VFETK_H_

#include "maloc/maloc.h"
#include "mc/mc.h"
#include "apbs/vhal.h"
#include "apbs/vatom.h"
#include "apbs/valist.h"
#include "apbs/vcsm.h"
#include "apbs/vpbe.h"
#include "apbs/vunit.h"
#include "apbs/vgreen.h"

/* ///////////////////////////////////////////////////////////////////////////
// Class Vfetk: Parameters and datatypes
/////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vfetk: Definition
/////////////////////////////////////////////////////////////////////////// */

typedef struct Vfetk { 

  Vmem *vmem;         /* Memory management object */
  Gem *gm;            /* Grid manager (container class for master vertex
                       * and simplex lists as well as prolongation
                       * operator for updating after refinement ) */
  AM *am;             /* Multilevel algebra manager */

  Vpbe *pbe;          /* Poisson-Boltzmann object */
  Vcsm *csm;          /* Charge-simplex map */
  
} Vfetk;

/* ///////////////////////////////////////////////////////////////////////////
// Class Vfetk: Inlineable methods (vfetk.c)
/////////////////////////////////////////////////////////////////////////// */

#if !defined(VINLINE_VFETK)
    VEXTERNC Gem*    Vfetk_getGem(Vfetk *thee);
    VEXTERNC AM*     Vfetk_getAM(Vfetk *thee);
    VEXTERNC Vpbe*   Vfetk_getVpbe(Vfetk *thee);
    VEXTERNC Vcsm*   Vfetk_getVcsm(Vfetk *thee);
    VEXTERNC int     Vfetk_getAtomColor(Vfetk *thee, int iatom);
#else /* if defined(VINLINE_VFETK) */
#   define Vfetk_getGem(thee) ((thee)->gm)
#   define Vfetk_getAM(thee) ((thee)->am)
#   define Vfetk_getVpbe(thee) ((thee)->pbe)
#   define Vfetk_getVcsm(thee) ((thee)->csm)
#   define Vfetk_getAtomColor(thee, iatom) (Vatom_getPartID(Valist_getAtom(Vpbe_getValist(thee->pbe), iatom)))
#endif /* if !defined(VINLINE_VFETK) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vfetk: Non-Inlineable methods (vfetk.c)
/////////////////////////////////////////////////////////////////////////// */

VEXTERNC Vfetk*  Vfetk_ctor(Vpbe *pbe, Gem *gm, AM *am);
VEXTERNC int     Vfetk_ctor2(Vfetk *thee, Vpbe *apbe, Gem *gm, AM *am);
VEXTERNC void    Vfetk_dtor(Vfetk **thee);
VEXTERNC void    Vfetk_dtor2(Vfetk *thee);

VEXTERNC double* Vfetk_getSolution(Vfetk *thee, int *length);
VEXTERNC double  Vfetk_getLinearEnergy1(Vfetk *thee, int color);
VEXTERNC double  Vfetk_getLinearEnergy2(Vfetk *thee, int color);
VEXTERNC double  Vfetk_getEnergyNorm2(Vfetk *thee, int color);
VEXTERNC int     Vfetk_memChk(Vfetk *thee);
VEXTERNC void    Vfetk_setAtomColors(Vfetk *thee);

#endif /* ifndef _VFETK_H_ */
