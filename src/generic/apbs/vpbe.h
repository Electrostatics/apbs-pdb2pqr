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

#include "mc/mc.h"
#include "apbs/vhal.h"
#include "apbs/vatom.h"
#include "apbs/valist.h"
#include "apbs/vcsm.h"
#include "apbs/vacc.h"
#include "apbs/vunit.h"
#include "apbs/vgreen.h"

#if defined(HAVE_PMGC_H)
#   include "pmgc/mgmlsys.h"
#   include "pmgc/mgarray.h"
#endif

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpbe: Parameters and datatypes
/////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpbe: Definition
/////////////////////////////////////////////////////////////////////////// */

typedef struct Vpbe { 

  Vmem *vmem;         /* Memory management object */

  Valist *alist;      /* Atom (charge) list */
  Vacc *acc;          /* Accessibility object */
  Vgreen *green;      /* Green's function oracle */

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
  double soluteXlen;  /* Solute length in x-direction */
  double soluteYlen;  /* Solute length in y-direction */
  double soluteZlen;  /* Solute length in z-direction */
  double soluteCharge;
                      /* Charge of solute molecule (e) */

  int paramFlag;      /* Check to see if the parameters have been set */
  
} Vpbe;

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpbe: Inlineable methods (vpbe.c)
/////////////////////////////////////////////////////////////////////////// */

#if !defined(VINLINE_VPBE)
    VEXTERNC Valist* Vpbe_getValist(Vpbe *thee);
    VEXTERNC Vacc*   Vpbe_getVacc(Vpbe *thee);
    VEXTERNC Vgreen* Vpbe_getVgreen(Vpbe *thee);
    VEXTERNC double  Vpbe_getIonConc(Vpbe *thee);
    VEXTERNC double  Vpbe_getTemperature(Vpbe *thee);           
    VEXTERNC double  Vpbe_getSoluteDiel(Vpbe *thee); 
    VEXTERNC double  Vpbe_getSoluteRadius(Vpbe *thee);
    VEXTERNC double  Vpbe_getSoluteXlen(Vpbe *thee);
    VEXTERNC double  Vpbe_getSoluteYlen(Vpbe *thee);
    VEXTERNC double  Vpbe_getSoluteZlen(Vpbe *thee);
    VEXTERNC double* Vpbe_getSoluteCenter(Vpbe *thee);
    VEXTERNC double  Vpbe_getSoluteCharge(Vpbe *thee);
    VEXTERNC double  Vpbe_getSolventDiel(Vpbe *thee);
    VEXTERNC double  Vpbe_getSolventRadius(Vpbe *thee);
    VEXTERNC double  Vpbe_getIonRadius(Vpbe *thee);
    VEXTERNC double  Vpbe_getXkappa(Vpbe *thee);
    VEXTERNC double  Vpbe_getDeblen(Vpbe *thee);
    VEXTERNC double  Vpbe_getZkappa2(Vpbe *thee);
    VEXTERNC double  Vpbe_getZmagic(Vpbe *thee);
#else /* if defined(VINLINE_VPBE) */
#   define Vpbe_getValist(thee) ((thee)->alist)
#   define Vpbe_getVacc(thee) ((thee)->acc)
#   define Vpbe_getVgreen(thee) ((thee)->green)
#   define Vpbe_getIonConc(thee) ((thee)->ionConc)
#   define Vpbe_getTemperature(thee) ((thee)->T)           
#   define Vpbe_getSoluteDiel(thee) ((thee)->soluteDiel) 
#   define Vpbe_getSoluteCenter(thee) ((thee)->soluteCenter)
#   define Vpbe_getSoluteRadius(thee) ((thee)->soluteRadius)
#   define Vpbe_getSoluteXlen(thee) ((thee)->soluteXlen)
#   define Vpbe_getSoluteYlen(thee) ((thee)->soluteYlen)
#   define Vpbe_getSoluteZlen(thee) ((thee)->soluteZlen)
#   define Vpbe_getSoluteCharge(thee) ((thee)->soluteCharge)
#   define Vpbe_getSolventDiel(thee) ((thee)->solventDiel)
#   define Vpbe_getSolventRadius(thee) ((thee)->solventRadius)
#   define Vpbe_getIonRadius(thee) ((thee)->ionRadius)
#   define Vpbe_getXkappa(thee) ((thee)->xkappa)
#   define Vpbe_getDeblen(thee) ((thee)->deblen)
#   define Vpbe_getZkappa2(thee) ((thee)->zkappa2)
#   define Vpbe_getZmagic(thee) ((thee)->zmagic)
#endif /* if !defined(VINLINE_VPBE) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpbe: Non-Inlineable methods (vpbe.c)
/////////////////////////////////////////////////////////////////////////// */

VEXTERNC Vpbe*   Vpbe_ctor(Valist *alist, double ionConc, double ionRadius,
                    double T, double soluteDiel, double solventDiel,
                    double solventRadius);
VEXTERNC int     Vpbe_ctor2(Vpbe *thee, Valist *alist, double ionConc, 
		    double ionRadius, double T, double soluteDiel, double
                    solventDiel, double solventRadius);

VEXTERNC void    Vpbe_dtor(Vpbe **thee);
VEXTERNC void    Vpbe_dtor2(Vpbe *thee);

VEXTERNC double  Vpbe_getCoulombEnergy1(Vpbe *thee);
VEXTERNC int     Vpbe_memChk(Vpbe *thee);


#endif /* ifndef _VALIST_H_ */
