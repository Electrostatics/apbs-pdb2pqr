/** @defgroup Valist Valist class
 *  @brief    Container class for list of atom objects
 */

/** 
 *  @file     valist.h
 *  @ingroup  Valist
 *  @brief    Contains declarations for class Valist
 *  @version  $Id$
 *  @author   Nathan A. Baker
 *
 *  @attention
 *  @verbatim
 *
 * APBS -- Adaptive Poisson-Boltzmann Solver
 *
 * Nathan A. Baker (baker@biochem.wustl.edu)
 * Dept. of Biochemistry and Molecular Biophysics
 * Center for Computational Biology
 * Washington University in St. Louis
 *
 * Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 2003.  Washington University in St. Louis.
 * All Rights Reserved.
 * Portions Copyright (c) 1999-2003.  The Regents of the University of
 * California.  
 * Portions Copyright (c) 1995.  Michael Holst.
 *
 * This file is part of APBS.
 *
 * APBS is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * APBS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with APBS; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
 *
 * @endverbatim
 */

#ifndef _VALIST_H_
#define _VALIST_H_

/* Generic headers */
#include "maloc/maloc.h"
#include "apbs/vhal.h"

/* Headers specific to this file */
#include "apbs/vatom.h"
#include "apbs/vparam.h"

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

    /** @brief   Get x-coordinate of molecule center
     *  @ingroup Valist
     *  @author  Nathan Baker
     *  @param   thee Valist object
     *  @return  X-coordinate of molecule center
     */
    VEXTERNC double Valist_getCenterX(Valist *thee);

    /** @brief   Get y-coordinate of molecule center
     *  @ingroup Valist
     *  @author  Nathan Baker
     *  @param   thee Valist object
     *  @return  Y-coordinate of molecule center
     */ 
    VEXTERNC double Valist_getCenterY(Valist *thee);

    /** @brief   Get z-coordinate of molecule center
     *  @ingroup Valist
     *  @author  Nathan Baker
     *  @param   thee Valist object
     *  @return  Z-coordinate of molecule center
     */ 
    VEXTERNC double Valist_getCenterZ(Valist *thee);


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
#   define Valist_getCenterX(thee) ((thee)->center[0])
#   define Valist_getCenterY(thee) ((thee)->center[1])
#   define Valist_getCenterZ(thee) ((thee)->center[2])
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

/** 
 * @brief  Fill atom list with information from a PQR file
 * @ingroup Valist
 * @author  Nathan Baker
 * @param   thee  Valist object to be filled
 * @param   iodev  Input device type (FILE/BUFF/UNIX/INET)
 * @param   iofmt  Input device format (ASCII/XDR)
 * @param   thost  Input hostname (for sockets)
 * @param   fname  Input FILE/BUFF/UNIX/INET name
 * @return  1 if successful, 0 otherwise
 * @note  \li A PQR file has PDB structure with charge and radius in the last
 *            two columns instead of weight and occupancy
 *        \li We don't actually respect PDB format; instead recognize
 *            whitespace- or tab-delimited fields which allows us to deal with
 *            structures with coordinates > 999 or < -999.
 */
VEXTERNC int Valist_readPQR(Valist *thee, const char *iodev, 
  const char *iofmt, const char *thost, const char *fname);

/** 
 * @brief  Fill atom list with information from a PDB file
 * @ingroup Valist
 * @author  Nathan Baker
 * @param   thee  Valist object to be filled
 * @param   param  A pre-initialized parameter object
 * @param   iodev  Input device type (FILE/BUFF/UNIX/INET)
 * @param   iofmt  Input device format (ASCII/XDR)
 * @param   thost  Input hostname (for sockets)
 * @param   fname  Input FILE/BUFF/UNIX/INET name
 * @return  1 if successful, 0 otherwise
 * @note  \li A PQR file has PDB structure with charge and radius in the last
 *            two columns instead of weight and occupancy
 *        \li We don't actually respect PDB format; instead recognize
 *            whitespace- or tab-delimited fields which allows us to deal with
 *            structures with coordinates > 999 or < -999.
 */
VEXTERNC int Valist_readPDB(Valist *thee, Vparam *param, const char *iodev,
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
