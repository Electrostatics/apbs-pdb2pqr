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
// File:     vpbe.h    < vpbe.c >
//
// Purpose:  
//    Class Vpbe:  
//      The Poisson-Boltzmann equation master class
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */

#ifndef _VPBE_H_
#define _VPBE_H_

#include "mc/vhal.h"
#include "mc/vgm.h"
#include "mc/vram.h"
#include "mc/ves.h"
#include "mc/am.h"
#include "mc/alg.h"
#include "mc/bvec.h"

#include "apbs/vhal.h"
#include "apbs/vatom.h"
#include "apbs/valist.h"
#include "apbs/vcsm.h"
#include "apbs/vacc.h"
#include "apbs/vunit.h"

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpbe: Parameters and datatypes
/////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpbe: Definition
/////////////////////////////////////////////////////////////////////////// */

typedef struct Vpbe { 

  Valist *alist;      /* Atom (charge) list */
  Vgm *gm;            /* Grid manager (container class for master vertex
                       * and simplex lists as well as prolongation
                       * operator for updating after refinement ) */


  Vacc *acc;          /* Accessibility object */
  Vcsm *csm;          /* Charge-simplex map */
  Vmem *vmem;         /* Memory management object */

  int methFlag;       /* Method of solution
                       *  0 ==> MC (adaptive multilevel FEM) 
                       *  1 ==> PMGC (multigrid) */

  double ionConc;     /* Ionic strength (M) */
  double T;           /* Temperature (K) */
  double soluteDiel;  /* Solute dielectric constant (unitless) */
  double solventDiel; /* Solvent dielectric constant (unitless) */
  double solventRadius;
                      /* Solvent probe radius (angstroms) for accessibility;
                       * determining defining volumes for the dielectric
                       * coefficient */
  double ionRadius;   /* Ion probe radius (angstroms) for accessibility;
                       * determining defining volumes for the ionic strength
                       * coefficient */

  double xkappa;      /* Debye-Huckel parameter */
  double deblen;      /* Debye length */
  double zkappa2;     /* Square of modified Debye-Huckel parameter */
  double zmagic;      /* Delta function scaling parameter */

  double soluteCenter[3];
                      /* Center of solute molecule (A) */
  double soluteRadius;/* Radius of solute molecule (A) */
  double soluteMaxX;  /* Max dist from solute center in x-direction */
  double soluteMaxY;  /* Max dist from solute center in y-direction */
  double soluteMaxZ;  /* Max dist from solute center in z-direction */
  double soluteCharge;
                      /* Charge of solute molecule (e) */

  int paramFlag;      /* Check to see if the parameters have been set */
  
} Vpbe;

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpbe: Inlineable methods (vpbe.c)
/////////////////////////////////////////////////////////////////////////// */

#if !defined(VINLINE_VPBE)
    VEXTERNC Valist* Vpbe_getValist(Vpbe *thee);
    VEXTERNC Vgm*    Vpbe_getVgm(Vpbe *thee);
    VEXTERNC Vacc*   Vpbe_getVacc(Vpbe *thee);
    VEXTERNC Vcsm*   Vpbe_getVcsm(Vpbe *thee);
    VEXTERNC double  Vpbe_getIonConc(Vpbe *thee);
    VEXTERNC double  Vpbe_getTemperature(Vpbe *thee);           
    VEXTERNC double  Vpbe_getSoluteDiel(Vpbe *thee); 
    VEXTERNC double  Vpbe_getSoluteRadius(Vpbe *thee);
    VEXTERNC double  Vpbe_getSoluteMaxX(Vpbe *thee);
    VEXTERNC double  Vpbe_getSoluteMaxY(Vpbe *thee);
    VEXTERNC double  Vpbe_getSoluteMaxZ(Vpbe *thee);
    VEXTERNC double* Vpbe_getSoluteCenter(Vpbe *thee);
    VEXTERNC double  Vpbe_getSoluteCharge(Vpbe *thee);
    VEXTERNC double  Vpbe_getSolventDiel(Vpbe *thee);
    VEXTERNC double  Vpbe_getSolventRadius(Vpbe *thee);
    VEXTERNC double  Vpbe_getIonRadius(Vpbe *thee);
    VEXTERNC double  Vpbe_getXkappa(Vpbe *thee);
    VEXTERNC double  Vpbe_getDeblen(Vpbe *thee);
    VEXTERNC double  Vpbe_getZkappa2(Vpbe *thee);
    VEXTERNC double  Vpbe_getZmagic(Vpbe *thee);
    VEXTERNC int     Vpbe_getAtomColor(Vpbe *thee, int iatom);
#else /* if defined(VINLINE_VPBE) */
#   define Vpbe_getValist(thee) ((thee)->alist)
#   define Vpbe_getVgm(thee) ((thee)->gm)
#   define Vpbe_getVacc(thee) ((thee)->acc)
#   define Vpbe_getVcsm(thee) ((thee)->csm)
#   define Vpbe_getIonConc(thee) ((thee)->ionConc)
#   define Vpbe_getTemperature(thee) ((thee)->T)           
#   define Vpbe_getSoluteDiel(thee) ((thee)->soluteDiel) 
#   define Vpbe_getSoluteCenter(thee) ((thee)->soluteCenter)
#   define Vpbe_getSoluteRadius(thee) ((thee)->soluteRadius)
#   define Vpbe_getSoluteMaxX(thee) ((thee)->soluteMaxX)
#   define Vpbe_getSoluteMaxY(thee) ((thee)->soluteMaxY)
#   define Vpbe_getSoluteMaxZ(thee) ((thee)->soluteMaxZ)
#   define Vpbe_getSoluteCharge(thee) ((thee)->soluteCharge)
#   define Vpbe_getSolventDiel(thee) ((thee)->solventDiel)
#   define Vpbe_getSolventRadius(thee) ((thee)->solventRadius)
#   define Vpbe_getIonRadius(thee) ((thee)->ionRadius)
#   define Vpbe_getXkappa(thee) ((thee)->xkappa)
#   define Vpbe_getDeblen(thee) ((thee)->deblen)
#   define Vpbe_getZkappa2(thee) ((thee)->zkappa2)
#   define Vpbe_getZmagic(thee) ((thee)->zmagic)
#   define Vpbe_getAtomColor(thee, iatom) ((thee)->csm->colors[iatom])
#endif /* if !defined(VINLINE_VPBE) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpbe: Non-Inlineable methods (vpbe.c)
/////////////////////////////////////////////////////////////////////////// */

VEXTERNC Vpbe*   Vpbe_ctor(Valist *alist, Vgm *gm, int methFlag);
VEXTERNC int     Vpbe_ctor2(Vpbe *thee, Valist *alist, Vgm *gm, int methFlag);
VEXTERNC void    Vpbe_dtor(Vpbe **thee);
VEXTERNC void    Vpbe_dtor2(Vpbe *thee);

VEXTERNC void    Vpbe_initialize(Vpbe *thee, double ionConc, double ionRadius, 
                    double T, double soluteDiel, double solventDiel, 
                    double solventRadius); 
VEXTERNC double* Vpbe_getSolution(Vpbe *thee, AM *am, int *length);
VEXTERNC double  Vpbe_getLinearEnergy1(Vpbe *thee, AM *am, int color);
VEXTERNC double  Vpbe_getLinearEnergy2(Vpbe *thee, AM *am, int color);
VEXTERNC double  Vpbe_getEnergyNorm2(Vpbe *thee, Alg *am, int color);
VEXTERNC double  Vpbe_getCoulombEnergy1(Vpbe *thee);
VEXTERNC int     Vpbe_memChk(Vpbe *thee);
VEXTERNC void    Vpbe_setAtomColors(Vpbe *thee);


#endif /* ifndef _VALIST_H_ */
