/** @defgroup Valist Valist class
 *  @brief    Container class for list of atom objects
 */

/** 
 *  @file     valist.h
 *  @ingroup  Valist
 *  @brief    Contains declarations for class Valist
 *  @version  $Id$
 *  @author   Nathan A. Baker
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
 * Copyright (c) 1999-2002.  Nathan A. Baker.  All Rights Reserved.
 *
 * Permission to use, copy, modify, and distribute this software and its
 * documentation for educational, research, and not-for-profit purposes,
 * without fee and without a signed licensing agreement, is hereby granted,
 * provided that the above copyright notice, this paragraph and the
 * following two paragraphs appear in all copies, modifications, and
 * distributions.
 *
 * IN NO EVENT SHALL THE AUTHORS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
 * SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE
 * AUTHORS HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * THE AUTHORS SPECIFICALLY DISCLAIM ANY WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE.  THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED
 * HEREUNDER IS PROVIDED "AS IS".  THE AUTHORS HAVE NO OBLIGATION TO PROVIDE
 * MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
 *
 * @endverbatim
 */

#ifndef _VALIST_H_
#define _VALIST_H_

#include "maloc/maloc.h"
#include "apbs/vhal.h"
#include "apbs/vatom.h"

/** 
 *  @struct  Valist
 *  @ingroup Valist
 *  @author  Nathan Baker
 *  @brief   Container class for list of atom objects
 */
struct Valist { 

  int number;         /**< Number of atoms in list */
  double center[3];   /**< Molecule center (xmin - xmax)/2, etc.*/
  double mincrd[3];   /**< Minimum coordinates */
  double maxcrd[3];   /**< Maximum coordinates */
  double maxrad;      /**< Maximum radius */
  double charge;      /**< Net charge */
  Vatom *atoms;       /**< Atom list */
  Vmem *vmem;         /**< Memory management object */

};

/** @typedef Valist
 *  @ingroup Valist
 *  @brief Declaration of the Valist class as the Valist structure
 */
typedef struct Valist Valist;

/* ///////////////////////////////////////////////////////////////////////////
// Class Valist: Inlineable methods (valist.c)
/////////////////////////////////////////////////////////////////////////// */

#if !defined(VINLINE_VATOM)

    /** @brief   Get actual array of atom objects from the list
     *  @ingroup Valist
     *  @author  Nathan Baker
     *  @param   thee Valist object
     *  @return  Array of atom objects 
     */
    VEXTERNC Vatom* Valist_getAtomList(Valist *thee);

    /** @brief   Get number of atoms in the list
     *  @ingroup Valist
     *  @author  Nathan Baker
     *  @param   thee Valist object
     *  @return  Number of atoms in list 
     */
    VEXTERNC int    Valist_getNumberAtoms(Valist *thee);

    /** @brief   Get pointer to particular atom in list
     *  @ingroup Valist
     *  @author  Nathan Baker
     *  @param   thee Valist object
     *  @param   i    Index of atom of interest
     *  @return  Pointer to atom object i
     */
    VEXTERNC Vatom* Valist_getAtom(Valist *thee, int i);

    /** @brief   Get total memory allocated for this object and its members
     *  @ingroup Valist
     *  @author  Nathan Baker
     *  @param   thee Valist object
     *  @return  Total memory in bytes
     */
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

/** @brief   Construct the atom list object
 *  @ingroup Valist
 *  @author  Nathan Baker
 *  @returns Pointer to newly allocated (empty) atom list 
 */
VEXTERNC Valist* Valist_ctor();

/** @brief   FORTRAN stub to construct the atom list object
 *  @ingroup Valist
 *  @author  Nathan Baker
 *  @param   thee Storage for new atom list
 *  @returns 1 if successful, 0 otherwise
 */
VEXTERNC int     Valist_ctor2(Valist *thee);

/** @brief   Destroys atom list object
 *  @ingroup Valist
 *  @author  Nathan Baker
 *  @param   thee Pointer to memory location for atom list object
 */
VEXTERNC void    Valist_dtor(Valist **thee);

/** @brief   FORTRAN stub to destroy atom list object
 *  @ingroup Valist
 *  @author  Nathan Baker
 *  @param   thee Pointer to atom list object
 */
VEXTERNC void    Valist_dtor2(Valist *thee);

/** @brief   Fill atom list with information from a PQR file
 *  @note     \liA PQR file has PDB structure with charge and radius in the last
 *             two columns instead of weight and occupancy
 *            \liThe original PDB reading algorithm was based on routines by Phil
 *             Hunenberger; it has since been modified by Nathan Baker
 *  @ingroup Valist
 *  @author  Nathan Baker
 *  @param   thee  Valist object to be filled
 *  @param   iodev  Input device type (FILE/BUFF/UNIX/INET)
 *  @param   iofmt  Input device format (ASCII/XDR)
 *  @param   thost  Input hostname (for sockets)
 *  @param   fname  Input FILE/BUFF/UNIX/INET name
 *  @return  1 if successful, 0 otherwise
 */
VEXTERNC int     Valist_readPQR(Valist *thee, const char *iodev, 
                 const char *iofmt, const char *thost, const char *fname);

/** @brief   Build rectangular prismatic finite element mesh which surrounds
 *           molecule contained in Valist object.
 * 
 *           The mesh will have rectangular sides and consist of 6 simplices
 *           and 8 vertices will all boundaries Dirichlet.  The mesh will be
 *           written to the sock specified by the parameters below.
 * 
 *  @bug     This routine has not been tested.
 *  @ingroup Valist
 *  @author  Nathan Baker (mesh topology based on MCSF file by Mike Holst)
 *  @param   thee   Valist object to use as reference for the mesh
 *  @param   size   The factor by which the mesh is larger than the
 *                  biomolecule.  In other words, if the smallest box
 *                  containing the biomolecule is dx x dy x dz, then the mesh
 *                  built by this routine will have dimensions (size*dx) x
 *                  (size*dy) x (size*dz).
 *  @param   iodev  Output device type (FILE/BUFF/UNIX/INET)
 *  @param   iofmt  Output device format (ASCII/XDR)
 *  @param   thost  Output hostname (for sockets)
 *  @param   fname  Output FILE/BUFF/UNIX/INET name
 */
VEXTERNC void    Valist_buildMesh(Valist *thee, double size, const char *iodev,
                   const char *iofmt, const char *thost, const char *fname);

#endif /* ifndef _VALIST_H_ */
