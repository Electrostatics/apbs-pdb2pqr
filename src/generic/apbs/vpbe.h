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
#include "mc/bvec.h"

#include "apbs/vatom.h"
#include "apbs/valist.h"
#include "apbs/vcsm.h"
#include "apbs/vhash.h"
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
  AM *am;             /* Algebraic multilevel structure; wraps solver
                       * portion of MC.  am is the AM object for the problem to
                       * be solved.  amRef is an optional AM object for a
                       * reference problem to be solved on the same mesh. */


  Vhash *solvHash;    /* Atomic hash table for solvent accessibility */
  Vhash *ionHash;     /* Atomic hash table for ionic accessibility */
  Vcsm *csm;          /* Charge-simplex map */

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
  double soluteRadius;
                      /* Radius of solute molecule (A) */
  double soluteCharge;
                      /* Charge of solute molecule (e) */

  int paramFlag;      /* Check to see if the parameters have been set */
  
} Vpbe;

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpbe: Inlineable methods (vpbe.c)
/////////////////////////////////////////////////////////////////////////// */

#if !defined(VINLINE_VPBE)
#else /* if defined(VINLINE_VPBE) */
#endif /* if !defined(VINLINE_VPBE) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpbe: Non-Inlineable methods (vpbe.c)
/////////////////////////////////////////////////////////////////////////// */

VEXTERNC Vpbe*   Vpbe_ctor(Valist *alist, Vgm *gm, AM *am);
VEXTERNC int     Vpbe_ctor2(Vpbe *thee, Valist *alist, Vgm *gm, AM *am);
VEXTERNC void    Vpbe_dtor(Vpbe **thee);
VEXTERNC void    Vpbe_dtor2(Vpbe *thee);

VEXTERNC void    Vpbe_initialize(Vpbe *thee, double ionConc, double ionRadius, 
                    double T, double soluteDiel, double solventDiel, 
                    double solventRadius); 

VEXTERNC Valist* Vpbe_getValist(Vpbe *thee);
VEXTERNC Vgm*    Vpbe_getVgm(Vpbe *thee);
VEXTERNC AM*     Vpbe_getAM(Vpbe *thee);
VEXTERNC Vhash*  Vpbe_getSolventHash(Vpbe *thee);
VEXTERNC Vhash*  Vpbe_getIonHash(Vpbe *thee);
VEXTERNC Vcsm*   Vpbe_getVcsm(Vpbe *thee);

VEXTERNC double Vpbe_getIonConc(Vpbe *thee);
VEXTERNC double Vpbe_getIonRadius(Vpbe *thee);
VEXTERNC double Vpbe_getTemperature(Vpbe *thee);           
VEXTERNC double Vpbe_getSoluteDiel(Vpbe *thee); 
VEXTERNC double* Vpbe_getSoluteCenter(Vpbe *thee);
VEXTERNC double Vpbe_getSoluteRadius(Vpbe *thee);
VEXTERNC double Vpbe_getSoluteCharge(Vpbe *thee);
VEXTERNC double Vpbe_getSolventDiel(Vpbe *thee);
VEXTERNC double Vpbe_getSolventRadius(Vpbe *thee);
VEXTERNC double Vpbe_getXkappa(Vpbe *thee);
VEXTERNC double Vpbe_getDeblen(Vpbe *thee);
VEXTERNC double Vpbe_getZkappa2(Vpbe *thee);
VEXTERNC double Vpbe_getZmagic(Vpbe *thee);
VEXTERNC double* Vpbe_getSolution(Vpbe *thee, int *length);
VEXTERNC double Vpbe_getLinearEnergy1(Vpbe *thee, int color);
VEXTERNC double Vpbe_getCoulombEnergy1(Vpbe *thee);

#endif /* ifndef _VALIST_H_ */
