/* ///////////////////////////////////////////////////////////////////////////
/// MC = < Manifold Code >
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
// File:     valist.h    < valist.c >
//
// Purpose:  
//    Class Valist:  
//      List of atoms (Vatom) = a molecule.
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */

#ifndef _VALIST_H_
#define _VALIST_H_

#include "apbs/vatom.h"
#include "mc/vhal.h"
#include "mc/vram.h"

/* ///////////////////////////////////////////////////////////////////////////
// Class Valist: Parameters and datatypes
/////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// Class Valist: Definition
/////////////////////////////////////////////////////////////////////////// */

typedef struct Valist { 

  int number;         /* Number of atoms in list */
  Vatom *atoms;       /* Atom list */

} Valist;

/* ///////////////////////////////////////////////////////////////////////////
// Class Valist: Inlineable methods (valist.c)
/////////////////////////////////////////////////////////////////////////// */

#if !defined(VINLINE_VATOM)
#else /* if defined(VINLINE_VATOM) */
#endif /* if !defined(VINLINE_VATOM) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Valist: Non-Inlineable methods (valist.c)
/////////////////////////////////////////////////////////////////////////// */

VEXTERNC Valist* Valist_ctor();
VEXTERNC int    Valist_ctor2(Valist *thee);
VEXTERNC void   Valist_dtor(Valist **thee);
VEXTERNC void   Valist_dtor2(Valist *thee);

VEXTERNC int    Valist_readPQR(Valist *thee, char *path);
VEXTERNC int    Valist_readPDB(Valist *thee, char *path, 
                  char *parameter_path);
VEXTERNC Vatom* Valist_getAtomList(Valist *thee);
VEXTERNC Vatom* Valist_getAtom(Valist *thee, int i);
VEXTERNC int    Valist_getNumberAtoms(Valist *thee);

#endif /* ifndef _VALIST_H_ */
