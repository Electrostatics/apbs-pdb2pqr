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

#include "maloc/maloc.h"
#include "apbs/vhal.h"
#include "apbs/vatom.h"

/* ///////////////////////////////////////////////////////////////////////////
// Class Valist: Parameters and datatypes
/////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// Class Valist: Definition
/////////////////////////////////////////////////////////////////////////// */

typedef struct Valist { 

  int number;         /* Number of atoms in list */
  Vatom *atoms;       /* Atom list */
  Vmem *vmem;         /* Memory management object */

} Valist;

/* ///////////////////////////////////////////////////////////////////////////
// Class Valist: Inlineable methods (valist.c)
/////////////////////////////////////////////////////////////////////////// */

#if !defined(VINLINE_VATOM)
    VEXTERNC Vatom* Valist_getAtomList(Valist *thee);
    VEXTERNC int    Valist_getNumberAtoms(Valist *thee);
    VEXTERNC Vatom* Valist_getAtom(Valist *thee, int i);
    VEXTERNC int    Valist_memChk(Valist *thee);
#else /* if defined(VINLINE_VATOM) */
#   define Valist_getAtomList(thee) ((thee)->atoms)
#   define Valist_getNumberAtoms(thee) ((thee)->number)
#   define Valist_getAtom(thee, i) (&((thee)->atoms[i]))
#   define Valist_memChk(thee) (Vmem_bytes((thee)->vmem))
#endif /* if !defined(VINLINE_VATOM) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Valist: Non-Inlineable methods (valist.c)
/////////////////////////////////////////////////////////////////////////// */

VEXTERNC Valist* Valist_ctor();
VEXTERNC int     Valist_ctor2(Valist *thee);
VEXTERNC void    Valist_dtor(Valist **thee);
VEXTERNC void    Valist_dtor2(Valist *thee);

VEXTERNC int     Valist_readPQR(Valist *thee, char *path);
VEXTERNC int     Valist_readPDB(Valist *thee, char *path, char *parameter_path);
VEXTERNC void    Valist_buildMesh(Valist *thee, double size, const char *iodev,
                   const char *iofmt, const char *thost, const char *fname);

#endif /* ifndef _VALIST_H_ */
