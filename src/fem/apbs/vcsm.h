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
// File:     vcsm.h    < vcsm.c >
//
// Purpose:  
//    Class Vcsm:  
//      A charge-simplex map for evaluating integrals of delta functions
//      with a Lagrange finite element basis.
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */

#ifndef _VCSM_H_
#define _VCSM_H_

#include "mc/mc.h"
#include "apbs/vhal.h"
#include "apbs/vatom.h"
#include "apbs/valist.h"

/* ///////////////////////////////////////////////////////////////////////////
// Class Vcsm: Parameters and datatypes
/////////////////////////////////////////////////////////////////////////// */

VEXTERNC void Gem_setExternalUpdateFunction(Gem *thee,
    void (*externalUpdate)(SS **simps, int num));

/* ///////////////////////////////////////////////////////////////////////////
// Class Vcsm: Definition
/////////////////////////////////////////////////////////////////////////// */

typedef struct Vcsm { 

  Valist *alist;      /* Atom (charge) list */
  int natom;          /* Size of thee->alist; redundant, but useful for
                       * convenience */
  Gem *gm;            /* Grid manager (container class for master vertex
                       * and simplex lists as well as prolongation
                       * operator for updating after refinement ) */
  int **sqm;          /* The map which gives the list charges associated with
                       * each simplex in gm->simplices.  The indices of
                       * the first dimension are associated with the
                       * simplex ID's in Vgm.  Each charge list (second 
                       * dimension) contains entries corresponding to
                       * indicies in thee->alist with lengths given in 
                       * thee->nsqm */
  int *nsqm;          /* The length of the charge lists in thee->sqm */
  int nsimp;          /* The _currently used) length of sqm, nsqm -- may not 
                       * always be up-to-date with Gem */
  int msimp;          /* The maximum number of entries that can be 
                       * accomodated by sqm or nsqm  -- saves on realloc's */
  int **qsm;          /* The inverse of sqm; the list of simplices
                       * associated with a given charge */
  int *nqsm;          /* The length of the simplex lists in thee->qsm */
  int *colors;        /* A list of atom colors for association with specific
                       * mesh partitions. */
  int initFlag;       /* Indicates whether the maps have been initialized
                       * yet */
  Vmem *vmem;         /* Memory management object */

} Vcsm;

/* ///////////////////////////////////////////////////////////////////////////
// Class Vcsm: Inlineable methods (vcsm.c)
/////////////////////////////////////////////////////////////////////////// */

#if !defined(VINLINE_VCSM)
    VEXTERNC Valist* Vcsm_getValist(Vcsm *thee);
    VEXTERNC Vgm*    Vcsm_getVgm(Vcsm *thee);
    VEXTERNC int     Vcsm_getNumberAtoms(Vcsm *thee, int isimp);
    VEXTERNC Vatom*  Vcsm_getAtom(Vcsm *thee, int iatom, int isimp);
    VEXTERNC int     Vcsm_getAtomIndex(Vcsm *thee, int iatom, int isimp);
    VEXTERNC int     Vcsm_getNumberSimplices(Vcsm *thee, int iatom);
    VEXTERNC SS*     Vcsm_getSimplex(Vcsm *thee, int isimp, int iatom);
    VEXTERNC int     Vcsm_getSimplexIndex(Vcsm *thee, int isimp, int iatom);
    VEXTERNC int     Vcsm_memChk(Vcsm *thee);
#else /* if defined(VINLINE_VCSM) */
#   define Vcsm_getValist(thee) ((thee)->alist)
#   define Vcsm_getVgm(thee) ((thee)->gm)
#   define Vcsm_getNumberAtoms(thee, isimp) ((thee)->nsqm[isimp])
#   define Vcsm_getAtom(thee, iatom, isimp) (Valist_getAtom((thee)->alist, ((thee)->sqm)[isimp][iatom]))
#   define Vcsm_getAtomIndex(thee, iatom, isimp) (((thee)->sqm)[isimp][iatom])
#   define Vcsm_getNumberSimplices(thee, iatom) (((thee)->nqsm)[iatom])
#   define Vcsm_getSimplex(thee, isimp, iatom) (Gem_SS((thee)->gm, ((thee)->qsm)[iatom][isimp]))
#   define Vcsm_getSimplexIndex(thee, isimp, iatom) (((thee)->qsm)[iatom][isimp])
#   define Vcsm_memChk(thee) (Vmem_bytes((thee)->vmem))
#endif /* if !defined(VINLINE_VCSM) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vcsm: Non-Inlineable methods (vcsm.c)
/////////////////////////////////////////////////////////////////////////// */

VEXTERNC Vcsm*   Vcsm_ctor(Valist *alist, Gem *gm);
VEXTERNC int     Vcsm_ctor2(Vcsm *thee, Valist *alist, Gem *gm);
VEXTERNC void    Vcsm_dtor(Vcsm **thee);
VEXTERNC void    Vcsm_dtor2(Vcsm *thee);


VEXTERNC void    Vcsm_init(Vcsm *thee);
VEXTERNC int     Vcsm_update(Vcsm *thee, SS **simps, int num);

#endif /* ifndef _VALIST_H_ */
