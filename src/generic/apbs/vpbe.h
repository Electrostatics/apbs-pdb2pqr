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

#include "maloc/maloc.h"
#include "apbs/vhal.h"
#include "apbs/vatom.h"
#include "apbs/valist.h"
#include "apbs/vacc.h"
#include "apbs/vunit.h"
#include "apbs/vgreen.h"

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

  double T;           /* Temperature (K) */
  double soluteDiel;  /* Solute dielectric constant (unitless) */
  double solventDiel; /* Solvent dielectric constant (unitless) */
  double solventRadius;
                      /* Solvent probe radius (angstroms) for accessibility;
                       * determining defining volumes for the dielectric
                       * coefficient */
  double bulkIonicStrength; /* Bulk ionic strength (M) */
  double maxIonRadius;      /* Max ion radius (A; used for calculating
                             * accessiblity and defining volumes for ionic
                             * strength coeffcients) */
  int    numIon;            /* Total number of ion species */
  double ionConc[MAXION];   /* Concentration (M) of each species */
  double ionRadii[MAXION];  /* Ionic radius (A) of each species */
  double ionQ[MAXION];      /* Charge (e) of each species */

  double xkappa;      /* Debye-Huckel parameter (bulk) */
  double deblen;      /* Debye length (bulk) */
  double zkappa2;     /* Square of modified Debye-Huckel parameter (bulk) */
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
    VEXTERNC double  Vpbe_getBulkIonicStrength(Vpbe *thee);
    VEXTERNC double  Vpbe_getMaxIonRadius(Vpbe *thee);
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
    VEXTERNC double  Vpbe_getSplineWin(Vpbe *thee);
    VEXTERNC double  Vpbe_getXkappa(Vpbe *thee);
    VEXTERNC double  Vpbe_getDeblen(Vpbe *thee);
    VEXTERNC double  Vpbe_getZkappa2(Vpbe *thee);
    VEXTERNC double  Vpbe_getZmagic(Vpbe *thee);
#else /* if defined(VINLINE_VPBE) */
#   define Vpbe_getValist(thee) ((thee)->alist)
#   define Vpbe_getVacc(thee) ((thee)->acc)
#   define Vpbe_getVgreen(thee) ((thee)->green)
#   define Vpbe_getBulkIonicStrength(thee) ((thee)->bulkIonicStrength)
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
#   define Vpbe_getMaxIonRadius(thee) ((thee)->maxIonRadius)
#   define Vpbe_getXkappa(thee) ((thee)->xkappa)
#   define Vpbe_getDeblen(thee) ((thee)->deblen)
#   define Vpbe_getZkappa2(thee) ((thee)->zkappa2)
#   define Vpbe_getZmagic(thee) ((thee)->zmagic)
#endif /* if !defined(VINLINE_VPBE) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpbe: Non-Inlineable methods (vpbe.c)
/////////////////////////////////////////////////////////////////////////// */

VEXTERNC Vpbe*   Vpbe_ctor(Valist *alist, int ionNum, double *ionConc, 
		    double *ionRadii, double *ionQ, double T,
                    double soluteDiel, double solventDiel,  
                    double solventRadius);
VEXTERNC int    Vpbe_ctor2(Vpbe *thee, Valist *alist, int ionNum, 
		    double *ionConc, double *ionRadii, double *ionQ, 
                    double T, double soluteDiel, 
                    double solventDiel, double solventRadius);
VEXTERNC int     Vpbe_getIons(Vpbe *thee, int *nion, double ionConc[MAXION],
                    double ionRadii[MAXION], double ionQ[MAXION]);

VEXTERNC void    Vpbe_dtor(Vpbe **thee);
VEXTERNC void    Vpbe_dtor2(Vpbe *thee);

VEXTERNC double  Vpbe_getCoulombEnergy1(Vpbe *thee);
VEXTERNC int     Vpbe_memChk(Vpbe *thee);


#endif /* ifndef _VALIST_H_ */
