/** @defgroup Vcsm Vcsm class
 *  @brief  A charge-simplex map for evaluating integrals of delta functions
 *          in a finite element setting
 */

/**
 *  @file      vcsm.h
 *  @brief     Contains declarations for the Vcsm class
 *  @ingroup   Vcsm
 *  @version   $Id$
 *  @author    Nathan A. Baker
 *  @attention
 *  @verbatim
 *
 * APBS -- Adaptive Poisson-Boltzmann Solver
 *
 * Nathan A. Baker (nbaker@wasabi.ucsd.edu)
 * Dept. of Chemistry and Biochemistry
 * University of California, San Diego
 *
 * Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 1999-2002. The Regents of the University of California 
 *                          (Regents).  All Rights Reserved.
 *
 * Permission to use, copy, modify, and distribute this software and its
 * documentation for educational, research, and not-for-profit purposes,
 * without fee and without a signed licensing agreement, is hereby granted,
 * provided that the above copyright notice, this paragraph and the
 * following two paragraphs appear in all copies, modifications, and
 * distributions.
 *
 * IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
 * REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE.  THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF
 * ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS".  REGENTS HAS NO OBLIGATION
 * TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
 * MODIFICATIONS.
 *
 * @endverbatim
 */

#ifndef _VCSM_H_
#define _VCSM_H_

#include "maloc/maloc.h"
#include "mc/mc.h"
#include "apbs/vhal.h"
#include "apbs/vatom.h"
#include "apbs/valist.h"

/** @brief   External function for FEtk Gem class to use during mesh refinement
 *  @ingroup Vcsm
 *  @author  Nathan Baker
 *  @param   thee   Gem (FEtk geometry manager) object
 *  @param   externalUpdate  Function pointer for function to call during 
 *           mesh refinement
 */
VEXTERNC void Gem_setExternalUpdateFunction(Gem *thee,
    void (*externalUpdate)(SS **simps, int num));

/** @brief   Charge-simplex map class
 *  @ingroup Vcsm
 *  @author  Nathan Baker
 */
struct Vcsm { 

  Valist *alist;      /**< Atom (charge) list */
  int natom;          /**< Size of thee->alist; redundant, but useful for
                       * convenience */
  Gem *gm;            /**< Grid manager (container class for master vertex
                       * and simplex lists as well as prolongation
                       * operator for updating after refinement ) */
  int **sqm;          /**< The map which gives the list charges associated with
                       * each simplex in gm->simplices.  The indices of
                       * the first dimension are associated with the
                       * simplex ID's in Vgm.  Each charge list (second 
                       * dimension) contains entries corresponding to
                       * indicies in thee->alist with lengths given in 
                       * thee->nsqm */
  int *nsqm;          /**< The length of the charge lists in thee->sqm */
  int nsimp;          /**< The _currently used) length of sqm, nsqm -- may not 
                       * always be up-to-date with Gem */
  int msimp;          /**< The maximum number of entries that can be 
                       * accomodated by sqm or nsqm  -- saves on realloc's */
  int **qsm;          /**< The inverse of sqm; the list of simplices
                       * associated with a given charge */
  int *nqsm;          /**< The length of the simplex lists in thee->qsm */
  int initFlag;       /**< Indicates whether the maps have been initialized
                       * yet */
  Vmem *vmem;         /**< Memory management object */

};

/** @typedef Vcsm
 *  @ingroup Vcsm
 *  @brief   Declaration of the Vcsm class as the Vcsm structure
 */
typedef struct Vcsm Vcsm;

/* ///////////////////////////////////////////////////////////////////////////
// Class Vcsm: Inlineable methods (vcsm.c)
/////////////////////////////////////////////////////////////////////////// */

#if !defined(VINLINE_VCSM)

    /** @brief   Get atom list
     *  @ingroup Vcsm
     *  @author  Nathan Baker
     *  @param   thee  Vcsm object
     *  @return  Pointer to Valist atom list
     */
    VEXTERNC Valist* Vcsm_getValist(Vcsm *thee);

    /** @brief   Get number of atoms associated with a simplex
     *  @ingroup Vcsm
     *  @author  Nathan Baker
     *  @param   thee  Vcsm object
     *  @param   isimp Simplex ID
     *  @return  Number of atoms associated with a simplex
     */
    VEXTERNC int     Vcsm_getNumberAtoms(Vcsm *thee, int isimp);

    /** @brief   Get list of atoms associated with a simplex
     *  @ingroup Vcsm
     *  @author  Nathan Baker
     *  @param   thee  Vcsm object
     *  @param   isimp Simplex ID
     *  @return  Array of atoms associated with a simplex
     */
    VEXTERNC Vatom*  Vcsm_getAtom(Vcsm *thee, int iatom, int isimp);

    /** @brief   Get ID of particular atom in a simplex
     *  @ingroup Vcsm
     *  @author  Nathan Baker
     *  @param   thee  Vcsm object
     *  @param   isimp Simplex ID
     *  @oaram   iatom Index of atom in Vcsm list
     *  @return  Index of atom in Valist object
     */
    VEXTERNC int     Vcsm_getAtomIndex(Vcsm *thee, int iatom, int isimp);

    /** @brief   Get number of simplices associated with an atom
     *  @ingroup Vcsm
     *  @author  Nathan Baker
     *  @param   thee  Vcsm object
     *  @param   iatom  Valist atom index
     *  @return  Number of simplices associated with an atom
     */
    VEXTERNC int     Vcsm_getNumberSimplices(Vcsm *thee, int iatom);

    /** @brief   Get particular simplex associated with an atom
     *  @ingroup Vcsm
     *  @author  Nathan Baker
     *  @param   thee  Vcsm object
     *  @param   iatom  Valist atom index
     *  @param   isimp  Index of simplex in Vcsm list
     *  @return  Pointer to simplex object
     */ 
    VEXTERNC SS*     Vcsm_getSimplex(Vcsm *thee, int isimp, int iatom);

    /** @brief   Get index particular simplex associated with an atom
     *  @ingroup Vcsm
     *  @author  Nathan Baker
     *  @param   thee  Vcsm object
     *  @param   iatom  Valist atom index
     *  @param   isimp  Index of simplex in Vcsm list
     *  @return  Gem index of specified simplex
     */ 
    VEXTERNC int     Vcsm_getSimplexIndex(Vcsm *thee, int isimp, int iatom);

    /** @brief   Return the memory used by this structure (and its contents)
     *           in bytes
     *  @ingroup Vcsm
     *  @author  Nathan Baker
     *  @param   thee  Vcsm object
     *  @return  The memory used by this structure and its contents in bytes
     */
    VEXTERNC int     Vcsm_memChk(Vcsm *thee);

#else /* if defined(VINLINE_VCSM) */
#   define Vcsm_getValist(thee) ((thee)->alist)
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

/** @brief   Construct Vcsm object
 *  @ingroup Vcsm 
 *  @author  Nathan Baker
 *  @note    \li The initial mesh must be sufficiently coarse for the assignment
 *             procedures to be efficient
 *           \li The map is not built until Vcsm_init is called
 *  @param   alist  List of atoms
 *  @param   gm     FEtk geometry manager (defines mesh)
 *  @return  Pointer to newly allocated Vcsm object 
 */
VEXTERNC Vcsm*   Vcsm_ctor(Valist *alist, Gem *gm);

/** @brief   FORTRAN stub to construct Vcsm object
 *  @ingroup Vcsm 
 *  @author  Nathan Baker
 *  @note    \li The initial mesh must be sufficiently coarse for the assignment
 *             procedures to be efficient
 *           \li The map is not built until Vcsm_init is called
 *  @param   thee   Pointer to space for Vcsm object
 *  @param   alist  List of atoms
 *  @param   gm     FEtk geometry manager (defines mesh)
 *  @return  1 if successful, 0 otherwise 
 */
VEXTERNC int     Vcsm_ctor2(Vcsm *thee, Valist *alist, Gem *gm);

/** @brief   Destroy Vcsm object
 *  @ingroup Vcsm
 *  @author  Nathan Baker
 *  @param   thee Pointer to memory location for Vcsm object
 */
VEXTERNC void    Vcsm_dtor(Vcsm **thee);

/** @brief   FORTRAN stub to destroy Vcsm object
 *  @ingroup Vcsm
 *  @author  Nathan Baker
 *  @param   thee Pointer to Vcsm object
 */
VEXTERNC void    Vcsm_dtor2(Vcsm *thee);

/** @brief   Initialize charge-simplex map with mesh and atom data
 *  @ingroup Vcsm
 *  @author  Nathan Baker
 *  @param   thee Vcsm object to be initialized
 *  @note    The initial mesh must be sufficiently coarse for the assignment
 *            procedures to be efficient
 */
VEXTERNC void    Vcsm_init(Vcsm *thee);

/** @brief   Update the charge-simplex and simplex-charge maps after
 *           refinement
 *  @ingroup Vcsm
 *  @author  Nathan Baker
 *  @param   thee  Vcsm object to be updated
 *  @param   simps List of pointer to newly created (by refinement) simplex
 *           objects.  The first simplex is expected to be derived from the
 *           parent simplex and therefore have the same ID.  The remaining
 *           simplices are the children and should represent new entries in the
 *           charge-simplex map.
 *  @return  1 if successful, 0 otherwise 
 */
VEXTERNC int     Vcsm_update(Vcsm *thee, SS **simps, int num);

#endif /* ifndef _VCSM_H_ */
