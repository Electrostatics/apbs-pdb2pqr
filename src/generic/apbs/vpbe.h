/* ///////////////////////////////////////////////////////////////////////////
/// MC = < Manifold Code >
///
/// Copyright (C) 1993  Michael Holst
/// Copyright (C) 1999  Contributing authors listed in the code documentation
///                     retain the copyrights to their code contributions.
///
/// The program MC is copyrighted software, and there are restrictions on its
/// use, modification, and distribution; please refer to the document COPYING
/// which accompanies this software.
///
/// This program is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
/// See the document COPYING for more details.
///
/// You should have received a copy of the COPYING file with this program;
/// if not, write to the author Michael Holst at:  mholst@math.ucsd.edu
///
/// rcsid="$Id$"
/// /////////////////////////////////////////////////////////////////////// */

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

#include "apbs/vatom.h"
#include "apbs/valist.h"
#include "apbs/vcsm.h"
#include "apbs/vhash.h"

/* ///////////////////////////////////////////////////////////////////////////
// Class Vcsm: Parameters and datatypes
/////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vcsm: Definition
/////////////////////////////////////////////////////////////////////////// */

typedef struct Vpbe { 

  Valist *alist;      /* Atom (charge) list */
  Vgm *gm;            /* Grid manager (container class for master vertex
                       * and simplex lists as well as prolongation
                       * operator for updating after refinement ) */

  Vhash *hash;        /* Atomic hash table */
  Vcsm *csm;          /* Charge-simplex map */

  double ionConc;     /* Ionic strength (M) */
  double T;           /* Temperature (K) */
  double soluteDiel;  /* Solute dielectric constant (unitless) */
  double solventDiel; /* Solvent dielectric constant (unitless) */
  double solventRadius;
                      /* Solvent radius (angstroms) */

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
  int csmFlag;        /* Check to see if the charge-simplex map has been
                       * constructued */
  
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

VEXTERNC Vpbe*   Vpbe_ctor(Valist *alist, Vgm *gm);
VEXTERNC int     Vpbe_ctor2(Vpbe *thee, Valist *alist, Vgm *gm);
VEXTERNC void    Vpbe_dtor(Vpbe **thee);
VEXTERNC void    Vpbe_dtor2(Vpbe *thee);

VEXTERNC void    Vpbe_initialize(Vpbe *thee, 
                    double ionConc, double T, double soluteDiel, 
                    double solventDiel, double solventRadius); 

VEXTERNC Valist* Vpbe_getValist(Vpbe *thee);
VEXTERNC Vgm*    Vpbe_getVgm(Vpbe *thee);
VEXTERNC Vhash*  Vpbe_getVhash(Vpbe *thee);
VEXTERNC Vcsm*   Vpbe_getVcsm(Vpbe *thee);

VEXTERNC double Vpbe_getIonConc(Vpbe *thee);
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

#endif /* ifndef _VALIST_H_ */
