/** @defgroup Vfetk Vfetk class
 *  @brief    FEtk master class (interface between FEtk and APBS)
 */

/**
 *  @file     vfetk.h
 *  @ingroup  Vfetk
 *  @brief    Contains declarations for class Vfetk
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

#ifndef _VFETK_H_
#define _VFETK_H_

#include "maloc/maloc.h"
#include "mc/mc.h"
#include "apbs/vhal.h"
#include "apbs/vatom.h"
/* #include "apbs/valist.h" */
#include "apbs/vcsm.h"
#include "apbs/vpbe.h"
#include "apbs/vunit.h"
#include "apbs/vgreen.h"
#include "apbs/vcap.h"
#include "apbs/pbeparm.h"
#include "apbs/femparm.h"

/**
 * @brief  Linear solver type 
 * @ingroup Vfetk
 * @note  Do not change these values; they correspond to settings in FEtk
 */
enum eVfetk_LsolvType {
    VLT_SLU=0,  /**< SuperLU direct solve */
    VLT_MG=1,  /**< Multigrid */
    VLT_CG=2,  /**< Conjugate gradient */
    VLT_BCG=3  /**< BiCGStab */
};

/**
 * @brief  Declare FEMparm_LsolvType type
 * @ingroup  Vfetk
 */
typedef enum eVfetk_LsolvType Vfetk_LsolvType;

/**
 * @brief  Non-linear solver type 
 * @ingroup Vfetk
 * @note  Do not change these values; they correspond to settings in FEtk
 */
enum eVfetk_NsolvType {
    VNT_NEW=0,  /**< Newton solver */
    VNT_INC=1,  /**< Incremental */
    VNT_ARC=2  /**< Psuedo-arclength */
};

/**
 * @brief  Declare FEMparm_NsolvType type
 * @ingroup  Vfetk
 */
typedef enum eVfetk_NsolvType Vfetk_NsolvType;

/**
 * @brief  Initial guess type
 * @ingroup Vfetk
 * @note  Do not change these values; they correspond to settings in FEtk
 */
enum eVfetk_GuessType {
    VGT_ZERO=0,  /**< Zero initial guess */
    VGT_DIRI=1,  /**< Dirichlet boundary condition initial guess */
    VGT_PREV=2  /**< Previous level initial guess */
};

/**
 * @brief  Declare FEMparm_GuessType type
 * @ingroup  Vfetk
 */
typedef enum eVfetk_GuessType Vfetk_GuessType;

/**
 * @brief  Preconditioner type
 * @ingroup Vfetk
 * @note  Do not change these values; they correspond to settings in FEtk
 */
enum eVfetk_PrecType {
    VPT_IDEN=0,  /**< Identity matrix */
    VPT_DIAG=1,  /**< Diagonal scaling */
    VPT_MG=2  /**< Multigrid */
};

/**
 * @brief  Declare FEMparm_GuessType type
 * @ingroup  Vfetk
 */
typedef enum eVfetk_PrecType Vfetk_PrecType;

/**
 *  @ingroup Vfetk
 *  @author  Nathan Baker
 *  @brief   Contains public data members for Vfetk class/module
 *
 *  Many of the routines and macros are borrowed from the main.c driver
 *  (written by Mike Holst) provided with the PMG code.
 *
 */
struct sVfetk { 

  Vmem *vmem;  /**< Memory management object */
  Gem *gm;  /**< Grid manager (container class for master vertex
             * and simplex lists as well as prolongation operator for updating
             * after refinement).  */
  AM *am;  /**< Multilevel algebra manager. */
  Aprx *aprx;  /**< Approximation manager. */
  PDE *pde;  /**< FEtk PDE object */
  Vpbe *pbe;  /**< Poisson-Boltzmann object */
  Vcsm *csm;  /**< Charge-simplex map */
  Vfetk_LsolvType lkey;  /**< Linear solver method */
  int lmax;  /**< Maximum number of linear solver iterations */
  double ltol;  /**< Residual tolerance for linear solver */
  Vfetk_NsolvType nkey;  /**< Nonlinear solver method */
  int nmax;  /**< Maximum number of nonlinear solver iterations */
  double ntol;  /**< Residual tolerance for nonlinear solver */
  Vfetk_GuessType gues;  /**< Initial guess method */
  Vfetk_PrecType lprec;  /**< Linear preconditioner */
  int pjac;  /**< Flag to print the jacobians (usually set this to -1, 
              * please) */
  PBEparm *pbeparm;  /**<  Generic PB parameters */
  FEMparm *feparm;  /**<  FEM-specific parameters */
  Vhal_PBEType type;  /**< Version of PBE to solve */
  int level;  /**< Refinement level (starts at 0) */

};

/** @typedef Vfetk
 *  @ingroup Vfetk
 *  @brief   Declaration of the Vfetk class as the Vfetk structure */
typedef struct sVfetk Vfetk;

/**
 * @brief  Vfetk LocalVar subclass
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @brief  Contains variables used when solving the PDE with FEtk
 */
struct sVfetk_LocalVar {
    double nvec[3];  /**< Normal vector for a simplex face */
    double vx[4][3];  /**< Vertex coordinates */
    double xq[3];  /**< Quadrature pt */
    double U[MAXV];  /**< Solution value */
    double dU[MAXV][3];  /**< Solution gradient */
    double W;  /**< Coulomb regularization term scalar value */
    double dW[3];  /**< Coulomb regularization term gradient */
    double d2W;  /**< Coulomb regularization term Laplacia */
    int sType;  /**< Simplex type */
    int fType;  /**< Face type */
    double diel;  /**< Dielectric value */
    double kappa2;  /**< Kappa^2 value */
    double A;  /**< Second-order differential term */
    double F;  /**< RHS characteristic function value */
    double B;  /**< Entire ionic strength term */
    double DB;  /**< Entire ionic strength term derivative */
    double jumpDiel;  /**< Dielectric value on one side of a simplex face */
    Vfetk *fetk;  /**< Pointer to the VFETK object */
    Vgreen *green;  /**< Pointer to a Green's function object */
    int initGreen;  /**< Boolean to designate whether Green's function 
                     * has been initialized */
    SS *simp;  /**< Pointer to the latest simplex object; set in initElement() 
                *  and delta() */
    VV *verts[4];  /**< Pointer to the latest vertices; set in initElement */
    int nverts;  /**< number of vertices in the simplex */
    double ionConc[MAXION];  /**< Counterion species' concentrations */
    double ionQ[MAXION];  /**< Counterion species' valencies */
    double ionRadii[MAXION];  /**< Counterion species' radii */
    double zkappa2; /**< Ionic strength parameters */
    double zks2; /**< Ionic strength parameters */
    double ionstr; /**< Ionic strength parameters (M) */
    int nion;  /**<  Number of ion species */
    double Fu_v;  /**< Store Fu_v value */
    double DFu_wv;  /**< Store DFu_wv value */
    double delta;  /**< Store delta value */
    double u_D;  /**< Store Dirichlet value */
    double u_T;  /**< Store true value */
};

/** 
 *  @ingroup Vfetk
 *  @brief   Declaration of the Vfetk_LocalVar subclass as the Vfetk_LocalVar 
 *           structure */
typedef struct sVfetk_LocalVar Vfetk_LocalVar;

#if !defined(VINLINE_VFETK)

    /** @brief   Get a pointer to the Gem (grid manager) object
     *  @ingroup Vfetk
     *  @author  Nathan Baker
     *  @param   thee  Vfetk object
     *  @return  Pointer to the Gem (grid manager) object
     */
    VEXTERNC Gem*    Vfetk_getGem(Vfetk *thee);

    /** @brief   Get a pointer to the AM (algebra manager) object
     *  @ingroup Vfetk
     *  @author  Nathan Baker
     *  @param   thee  Vfetk object
     *  @return  Pointer to the AM (algebra manager) object
     */
    VEXTERNC AM*     Vfetk_getAM(Vfetk *thee);

    /** @brief   Get a pointer to the Vpbe (PBE manager) object
     *  @ingroup Vfetk
     *  @author  Nathan Baker
     *  @param   thee  Vfetk object
     *  @return  Pointer to the Vpbe (PBE manager) object
     */
    VEXTERNC Vpbe*   Vfetk_getVpbe(Vfetk *thee);

    /** @brief   Get a pointer to the Vcsm (charge-simplex map) object
     *  @ingroup Vfetk
     *  @author  Nathan Baker
     *  @param   thee  Vfetk object
     *  @return  Pointer to the Vcsm (charge-simplex map) object
     */
    VEXTERNC Vcsm*   Vfetk_getVcsm(Vfetk *thee);

    /** @brief   Get the partition information for a particular atom
     *  @ingroup Vfetk
     *  @author  Nathan Baker
     *  @note    Friend function of Vatom
     *  @param   thee  Vfetk object
     *  @param   iatom Valist atom ID
     *  @returns Partition ID 
     */
    VEXTERNC int     Vfetk_getAtomColor(Vfetk *thee, int iatom);

#else /* if defined(VINLINE_VFETK) */
#   define Vfetk_getGem(thee) ((thee)->gm)
#   define Vfetk_getAM(thee) ((thee)->am)
#   define Vfetk_getVpbe(thee) ((thee)->pbe)
#   define Vfetk_getVcsm(thee) ((thee)->csm)
#   define Vfetk_getAtomColor(thee, iatom) (Vatom_getPartID(Valist_getAtom(Vpbe_getValist(thee->pbe), iatom)))
#endif /* if !defined(VINLINE_VFETK) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vfetk: Non-Inlineable methods (vfetk.c)
/////////////////////////////////////////////////////////////////////////// */

/** 
 * @brief  Constructor for Vfetk object
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @param  pbe  Vpbe (PBE manager) object
 * @param  type  Version of PBE to solve
 * @return  Pointer to newly allocated Vfetk object 
 * @note  This sets up the Gem, AM, and Aprx FEtk objects but does not create
 *         a mesh.  The easiest way to create a mesh is to then call
 *         Vfetk_genCube */
VEXTERNC Vfetk*  Vfetk_ctor(Vpbe *pbe, Vhal_PBEType type);

/** 
 * @brief  FORTRAN stub constructor for Vfetk object
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @param  thee Vfetk obeject memory address
 * @param  pbe  Vpbe (PBE manager) object
 * @param  type  Version of PBE to solve
 * @return  1 if successful, 0 otherwise
 * @note  This sets up the Gem, AM, and Aprx FEtk objects but does not create
 *         a mesh.  The easiest way to create a mesh is to then call
 *         Vfetk_genCube */
VEXTERNC int     Vfetk_ctor2(Vfetk *thee, Vpbe *pbe, Vhal_PBEType type);

/** @brief   Object destructor
 *  @ingroup Vfetk
 *  @author  Nathan Baker
 *  @param   thee   Pointer to memory location of object to be destroyed
 */
VEXTERNC void    Vfetk_dtor(Vfetk **thee);

/** @brief   FORTRAN stub object destructor
 *  @ingroup Vfetk
 *  @author  Nathan Baker
 *  @param   thee   Pointer to object to be destroyed
 */
VEXTERNC void    Vfetk_dtor2(Vfetk *thee);

/** @brief   Create an array containing the solution (electrostatic potential
 *           in units of \f$k_B T/e\f$) at the finest mesh level. 
 *  @ingroup Vfetk
 *  @author  Nathan Baker and Michael Holst
 *  @note    The user is responsible for destroying the newly created array
 *  @param   thee   Vfetk object
 *  @param   length Set to the length of the newly created array
 *  @return  Newly created array of length "length" (see above); the user is
 *           responsible for destruction
 */
VEXTERNC double* Vfetk_getSolution(Vfetk *thee, int *length);

/**
 * @brief  Set the parameter objects
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @param  thee  Vfetk object
 * @param  pbeparm  PBE generic parameters
 * @param  feparm  FEM-specific parameters
 */
VEXTERNC void Vfetk_setParameters(Vfetk *thee, PBEparm *pbeparm,
  FEMparm *feparm);

/** @brief   Return the total electrostatic energy
 * 
 *   Using the solution at the finest mesh level, get the electrostatic energy
 *   using the free energy functional for the Poisson-Boltzmann equation
 *   without removing any self-interaction terms (i.e., removing the reference
 *   state of isolated charges present in an infinite dielectric continuum with
 *   the same relative permittivity as the interior of the protein) and return
 *   the result in units of \f$k_B T\f$.  The argument color allows the user to
 *   control the partition on which this energy is calculated; if (color == -1)
 *   no restrictions are used.  The solution is obtained from the finest level
 *   of the passed AM object, but atomic data from the Vfetk object is used to
 *   calculate the energy.
 * 
 *  @ingroup Vfetk
 *  @author  Nathan Baker
 *  @param thee Vfetk object
 *  @param color Partition restriction; if non-negative, energy calculation is
 *               restricted to the specified partition.
 *  @param nonlin If 1, the NPBE energy functional is used, if 0 then the LPBE
 *                energy functional is used.
 *  @return Total electrostatic energy in units of \f$k_B T\f$.
 */
VEXTERNC double  Vfetk_energy(Vfetk *thee, int color, int nonlin);

/** @brief   Get the "mobile charge" and "polarization" contributions to the
 *           electrostatic energy.
 * 
 *           Using the solution at the finest mesh level, get the
 *           electrostatic energy due to the interaction of the mobile charges
 *           with the potential and polarization of the dielectric medium:
 *              \f[ G = \frac{1}{4 I_s} \sum_i c_i q_i^2 \int
 *              \overline{\kappa}^2(x) e^{-q_i u(x)} dx + \frac{1}{2} \int 
 *              \epsilon ( \nabla u )^2 dx \f]
 *           for the NPBE and
 *              \f[ G = \frac{1}{2} \int \overline{\kappa}^2(x) u^2(x) dx + 
 *              \frac{1}{2} \int \epsilon ( \nabla u )^2 dx \f]
 *           for the LPBE.  Here \f$i\f$ denotes the counterion species,
 *           \f$I_s\f$ is the bulk ionic strength, \f$\overline{\kappa}^2(x)\f$
 *           is the modified Debye-Huckel parameter, \f$c_i\f$ is the  
 *           concentration of species \f$i\f$, \f$q_i\f$ is the charge of
 *           species \f$i\f$, \f$\epsilon\f$ is the dielectric function, and
 *           \f$u(x)\f$ is the dimensionless electrostatic potential.  The
 *           energy is scaled to units of \f$k_b T\f$.
 *
 *  @ingroup Vfetk
 *  @author  Nathan Baker
 *  @param thee  Vfetk object
 *  @param color Partition restriction for energy evaluation, only used if
 *               non-negative
 *  @return The "mobile charge" and "polarization" contributions to the
 *           electrostatic energy in units of \f$k_B T\f$.
 */
VEXTERNC double  Vfetk_dqmEnergy(Vfetk *thee, int color);

/** @brief   Get the "fixed charge" contribution to the electrostatic energy
 *
 *           Using the solution at the finest mesh level, get the
 *           electrostatic energy due to the interaction of the fixed charges
 *           with the potential: \f[ G = \sum_i q_i u(r_i) \f]
 *           and return the result in units of \f$k_B T\f$.  Clearly, no
 *           self-interaction terms are removed.  A factor a 1/2 has to be
 *           included to convert this to a real energy.
 *
 *  @ingroup Vfetk
 *  @author  Nathan Baker
 *  @param   thee   Vfetk object
 *  @param   color Partition restriction for energy evaluation, only used if
 *               non-negative
 *  @returns The fixed charge electrostatic energy in units of \f$k_B T\f$.
 */
VEXTERNC double  Vfetk_qfEnergy(Vfetk *thee, int color);

/** @brief   Calculate the log determinant of the specified operator.
 *  @ingroup Vfetk
 *  @author  Nathan Baker and Stephen Bond
 *  @note    \li Only works with symmetric positive definite matrices
 *           \li Uses LU or Recycled Cholesky factorization and can be 
 *           very memory- and time-consuming for large matrices.
 *  @param   thee  Vfetk object
 *  @param   color Partition to evaluate over (ignored if <0)
 *  @param   oflag  Operator to evaluate:
 *           \li 0:  Helmholtz operator (NPBE tangent operator evaluated at zero
 *                 solution)
 *           \li 1:  Response function (NPBE tangent operator evaluated at NPBE
 *                 solution)
 *  @param   mflag  Method to use:
 *           \li 0:  Full nonsymmetric SuperLU factor with ROW/COL reordering
 *           \li 1:  Recycled symmetric Cholesky factor with no reordering
 *  @return  The log determinant of the specified operator
 *  @bug     color argument ignored
 */
VEXTERNC double Vfetk_lnDet(Vfetk *thee, int color, int oflag, int mflag);

/** @brief   Return the memory used by this structure (and its contents)
 *           in bytes
 *  @ingroup Vfetk
 *  @author  Nathan Baker
 *  @param   thee  Vfetk object
 *  @return  The memory used by this structure and its contents in bytes
 */
VEXTERNC int     Vfetk_memChk(Vfetk *thee);

/** @brief   Transfer color (partition ID) information frmo a partitioned mesh
 *           to the atoms.
 * 
 *           Transfer color information from partitioned mesh to the atoms.
 *           In the case that a charge is shared between two partitions, the
 *           partition color of the first simplex is selected.  Due to the
 *           arbitrary nature of this selection, THIS METHOD SHOULD ONLY BE
 *           USED IMMEDIATELY AFTER PARTITIONING!!!
 *  @warning This function should only be used immediately after mesh
 *           partitioning
 *  @ingroup Vfetk
 *  @author  Nathan Baker
 *  @note    This is a friend function of Vcsm
 *  @param   thee Vfetk object
 */
VEXTERNC void    Vfetk_setAtomColors(Vfetk *thee);

/** @brief   Writes a Bmat to disk in Harwell-Boeing sparse matrix format.
 * 
 *  @ingroup Vfetk
 *  @author  Stephen Bond
 *  @note    This is a friend function of Bmat
 *  @param   thee Bmat object
 *  @param   fname char Output filename
 *  @bug     Hardwired to only handle the single block symmetric case.
 */
VEXTERNC void    Bmat_printHB(Bmat *thee, char *fname);

/** @brief   Assembles the Cholesky factorization of a Bmat.
 * 
 *  @ingroup Vfetk
 *  @author  Stephen Bond
 *  @note    This is a friend function of Bmat
 *  @param   thee Bmat object
 *  @param   flag  Type of factor to be assembled:
 *           \li 0:  Full Cholesky factor
 *           \li 1:  Diagonal of the Cholesky factor
 */
VEXTERNC int     Bmat_choleskyFactor(Bmat *thee, int flag);

/** 
 * @brief   Returns the log(abs(det(D))) of a diagonal Mat, D.
 * @ingroup Vfetk
 * @author  Stephen Bond
 * @note    This is a friend function of Mat
 * @param   thee Mat object
 * @bug     Only works for Mat's of type RLN, CLN, or DRC.
 */
VEXTERNC double  Mat_lnDetDiag(Mat *thee);

/** 
 * @brief   Generates a simple cubic tetrahedral mesh
 * @ingroup Vfetk
 * @author  Nathan Baker (based on Mike Holst's Gem_makeCube code)
 * @param   thee Vfetk object
 * @param   center  desired mesh center (in &Aring;)
 * @param   length  desired mesh length (in &Aring;)
 * @return  1 if successful, 0 otherwise
 */
VEXTERNC int Vfetk_genCube(Vfetk *thee, double center[3], double length[3]);

/**
 * @brief  Constructs the FEtk PDE object 
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @param  fetk  Vfetk (parent) object
 * @return  Newly-allocated PDE object 
 * @bug  Not thread-safe */
VEXTERNC PDE* Vfetk_PDE_ctor(Vfetk *fetk);

/**
 * @brief  Intializes the FEtk PDE object 
 * @ingroup  Vfetk
 * @author  Nathan Baker (with code by Mike Holst)
 * @param  thee newly-allocated PDE object
 * @param  fetk  Vfetk (parent) object
 * @return  1 if successful, 0 otherwise
 * @bug  Not thread-safe */
VEXTERNC int Vfetk_PDE_ctor2(PDE *thee, Vfetk *fetk);

/**
 * @brief  Destroys FEtk PDE object
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @param  thee Pointer to PDE object memory
 * @note  Thread-safe
 */
VEXTERNC void Vfetk_PDE_dtor(PDE **thee);

/**
 * @brief  FORTRAN stub:  destroys FEtk PDE object
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @param  thee Pointer to PDE object
 * @note  Thread-safe
 */
VEXTERNC void Vfetk_PDE_dtor2(PDE *thee);

/** 
 * @brief  Do once-per-assembly initialization
 * @ingroup  Vfetk
 * @author  Nathan Baker and Mike Holst
 * @param  thee PDE object
 * @param ip  Integer parameter array (not used)
 * @param rp  Double parameter array (not used)
 * @note  Thread-safe */
VEXTERNC void Vfetk_PDE_initAssemble(PDE *thee, int ip[], double rp[]);

/** 
 * @brief  Do once-per-element initialization
 * @ingroup  Vfetk
 * @author  Nathan Baker and Mike Holst
 * @param  thee PDE object
 * @param  elementType  Material type of element
 * @param  chart  Chart in which vertex coordinates are provided (used here as
 *                bitfield for molecular accessibiliity)
 * @param tvx  Vertex coordinates
 * @param data  Simplex pointer
 * @todo  Jump term is not implemented
 * @bug This function is not thread-safe */
VEXTERNC void Vfetk_PDE_initElement(PDE *thee, int elementType, int chart,
  double tvx[][3], void *data);

/** 
 * @brief  Do once-per-face initialization
 * @ingroup  Vfetk
 * @author  Nathan Baker and Mike Holst
 * @param  thee PDE object
 * @param  faceType  (interior or various boundary types)
 * @param  chart  Chart in which vertex coordinates are provided (used here as
 *                bitfield for molecular accessibiliity)
 * @param tnvec  Coordinates of outward normal vector for face
 * @bug This function is not thread-safe */
VEXTERNC void Vfetk_PDE_initFace(PDE *thee, int faceType, int chart,
  double tnvec[]);

/** 
 * @brief  Do once-per-point initialization
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @param  thee PDE object
 * @param  pointType  (interior or various boundary types)
 * @param  chart  Chart in which vertex coordinates are provided (used here as
 *                bitfield for molecular accessibiliity)
 * @param txq  Coordinates of point
 * @param tU  Solution value at point
 * @param tdU  Solution gradient at point
 * @bug This function is not thread-safe 
 * @bug This function uses pre-defined boudnary definitions for the molecular
 *      surface. */
VEXTERNC void Vfetk_PDE_initPoint(PDE *thee, int pointType, int chart, 
  double txq[], double tU[], double tdU[][3]);

/** 
 * @brief  Evaluate strong form of PBE.  For interior points, this is:
 *  \f[ -\nabla \cdot \epsilon \nabla u + b(u) - f \f]
 *  where \f$b(u)\f$ is the (possibly nonlinear) mobile ion term and \f$f\f$ is
 *  the source charge distribution term (for PBE) or the induced surface charge
 *  distribution (for RPBE).  For an interior-boundary (simplex face) point,
 *  this is: 
 *  \f[ [\epsilon(x) \nabla u(x) \cdot n(x)]_{x=0^+} - [\epsilon(x) \nabla u(x)
 *  \cdot n(x)]_{x=0^-} \f]
 *  where \f$n(x)\f$ is the normal to the simplex face and the term represents
 *  the jump in dielectric displacement across the face.  There is no
 *  outer-boundary contribution for this problem.
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @param  thee PDE object
 * @param  key type of point (0 = interior, 1 = boundary, 2 = interior
 * boundary)
 * @param  F value of residual
 * @bug This function is not thread-safe 
 * @bug This function is not implemented (sets error to zero)
 */
VEXTERNC void Vfetk_PDE_Fu(PDE *thee, int key, double F[]);

/**
 * @brief  This is the weak form of the PBE; i.e. the strong form integrated
 * with a test function to give:
 * \f[ \int_\Omega \left[ \epsilon \nabla u \cdot \nabla v + b(u) v - f v
 * \right] dx \f]
 * where \f$b(u)\f$ denotes the mobile ion term.
 * @ingroup  Vfetk
 * @author  Nathan Baker and Mike Holst
 * @param  thee  PDE object
 * @param  key  Integrand to evaluate
 *   \li 0 => Interior weak form (see above)
 *   \li 1 => Boundary weak form (not required)
 * @param  V  test function at current point
 * @param  dV  test function gradient at current point
 * @return  Integrand value
 * @bug  This function is not thread-safe */
VEXTERNC double Vfetk_PDE_Fu_v(PDE *thee, int key, double V[], double dV[][3]);

/**
 * @brief  This is the linearization of the weak form of the PBE; e.g., for use
 * in a Newton iteration. This is the functional linearization of the strong
 * form integrated with a test function to give:
 * \f[ \int_\Omega \left[ \epsilon \nabla w \cdot \nabla v + b'(u) w v - f v
 * \right] dx \f]
 * where \f$b'(u)\f$ denotes the functional derivation of the mobile ion term.
 * @ingroup  Vfetk
 * @author  Nathan Baker and Mike Holst
 * @param  thee PDE object
 * @param  key  Integrand to evaluate
 *   \li  0 => Interior form (see above)
 *   \li  1 => Boundary form (not required)
 * @param  W trial function at current point
 * @param  dW trial function gradient at current point
 * @param  V test function at current point
 * @param  dV test function gradient at current point
 * @return  Integrand value
 * @bug  This function is not thread-safe */
VEXTERNC double Vfetk_PDE_DFu_wv(PDE *thee, int key, double W[], double dW[][3],
  double V[], double dV[][3]);

/**
 * @brief  Evaluate a (discretized) delta function source term at the given
 * point
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @param  thee  PDE object
 * @param  chart  Chart for point coordinates
 * @param  type  Vertex type
 * @param  txq  Point coordinates
 * @param  user  User-defined data (cast to vertex of interest)
 * @param  F  set to delta function value
 * @bug  This function is not thread-safe */
VEXTERNC void Vfetk_PDE_delta(PDE *thee, int type, int chart, double txq[],
  void *user, double F[]);

/**
 * @brief  Evaluate the Dirichlet boundary condition at the given point
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @param  thee  PDE object
 * @param  type  Vertex type (i.e., different types of boundary points)
 * @param  chart  Chart in which vertex coordinates are provided
 * @param  txq  Vertex coordinates
 * @param  F  Set to boundary condition value
 * @bug  This function is hard-coded to call only multiple-sphere
 * Debye-H&uuml; functions. 
 * @bug  This function is not thread-safe. */
VEXTERNC void Vfetk_PDE_u_D(PDE *thee, int type, int chart, double txq[],
  double F[]);

/**
 * @brief  Evaluate the "true solution" at the given point for comparison with
 * the numerical solution
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @param  thee  PDE object
 * @param  type  Vertex type (i.e., different types of boundary points)
 * @param  chart  Chart in which vertex coordinates are provided
 * @param  txq  Vertex coordinates
 * @param  F  Set to true solution value
 * @note  This function only returns zero.
 * @bug  This function is not thread-safe. */
VEXTERNC void Vfetk_PDE_u_T(PDE *thee, int type, int chart, double txq[],
  double F[]);

/**
 * @brief  Define the way manifold edges are bisected
 * @ingroup  Vfetk
 * @author  Nathan Baker and Mike Holst
 * @param  dim  intrinsic dimension of the manifold
 * @param  dimII  imbedding dimension of the manifold
 * @param  edgeType  edge type being refined
 * @param  chart  Accessibility bitfields
 *   \li chart[0]  Accessibility bitfield for first vertex of edge
 *   \li chart[1]  Accessibility bitfield for second vertex of edge
 *   \li chart[2]  Set to accessibility bitfield for new vertex
 * @param  vx  Vertex coordinates
 *   \li  vx[0]  First vertex coordinates 
 *   \li  vx[1]  Second vertex coordinates
 *   \li  vx[2]  Set to new vertex coordinates 
 * @note  This function is thread-safe. */
VEXTERNC void Vfetk_PDE_bisectEdge(int dim, int dimII, int edgeType, 
  int chart[], double vx[][3]);

/**
 * @brief  Map a boundary point to some pre-defined shape
 * @ingroup  Vfetk
 * @author  Nathan Baker and Mike Holst
 * @param  dim  intrinsic dimension of the manifold
 * @param  dimII  imbedding dimension of the manifold
 * @param  vertexType  type of boundary vertex
 * @param  chart  Vertex chart
 * @param  vx  Vertex coordinates
 * @note  This function is thread-safe and is a no-op */
VEXTERNC void Vfetk_PDE_mapBoundary(int dim, int dimII, int vertexType, 
  int chart, double vx[3]);

/**
 * @brief  User-defined error estimator -- in our case, a geometry-based
 * refinement method; forcing simplex refinement at the dielectric boundary and
 * (for non-regularized PBE) the charges.
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @param  dim  intrinsic dimension of the manifold
 * @param  dimII  imbedding dimension of the manifold
 * @param  simplexType  type of simplex 
 * @param  faceType  types of simplex faces
 * @param  vertexType  types of simplex vertices
 * @param  chart  Vertices' charts
 * @param  vx  Vertex coordinates
 * @param  simplex  Simplex pointer
 * @return 1 if mark simplex for refinement, 0 otherwise
 * @bug  This function is not thread-safe */
VEXTERNC int Vfetk_PDE_markSimplex(int dim, int dimII, int simplexType, 
  int faceType[4], int vertexType[4], int chart[], double vx[][3],
  void *simplex);

/**
 * @brief  Unify the chart for different coordinate systems -- a no-op for us.
 * @ingroup  Vfetk
 * @author Nathan Baker
 * @param  dim  Dimension of manifold 
 * @param  dimII  Imbedding dimension 
 * @param  objType  Don't know 
 * @param  chart  Charts of vertices
 * @param  vx  Vertex coordinates 
 * @param  dimV  Number of vertices 
 * @note  Thread-safe; a no-op */
VEXTERNC void Vfetk_PDE_oneChart(int dim, int dimII, int objType, int chart[],
  double vx[][3], int dimV);

/**
 * @brief  Energy functional.  This returns the energy (less delta function
 * terms) in the form:
 *   \f[ c^{-1}/2 \int (\epsilon (\nabla u)^2 + \kappa^2 (cosh u - 1)) dx \f]
 * for a 1:1 electrolyte where \f$c\f$ is the output from Vpbe_getZmagic.
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @param  thee  The PDE object
 * @param  key  interior (0) or boundary (1)
 * @returns  Energy value (in kT)
 * @bug  This function is not thread-safe. */
VEXTERNC double Vfetk_PDE_Ju(PDE *thee, int key);

/**
 * @brief  External hook to simplex subdivision routines in Gem.  Called each
 * time a simplex is subdivided (we use it to update the charge-simplex map)
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @param  simps  List of parent and children simplices
 * @param  num  Number of simplices in list
 * @bug  This function is not thread-safe.
 */
VEXTERNC void Vfetk_externalUpdateFunction(SS **simps, int num);


/**
 * @brief Initialize the bases for the trial or the test space, for a
 * particular component of the system, at all quadrature points on the master
 * simplex element.
 * @ingroup Vfetk
 * @author Mike Holst
 * @param key basis type to evaluate
 *        \li 0 => trial
 *        \li 1 => test
 *        \li 2 => trialB
 *        \li 3 => testB
 * @param dim  spatial dimension
 * @param comp which component of elliptic system to produce basis for
 * @param ndof set to number of degrees of freedom
 * @param dof set to degree of freedom per v/e/f/s
 * @note
 *   @verbatim
 *   The basis ordering is important.  For a fixed quadrature
 *   point iq, you must follow the following ordering in p[iq][],
 *   based on how you specify the degrees of freedom in dof[]:
 *   
 *   <v_0 vDF_0>,      <v_1 vDF_0>,      ..., <v_{nv} vDF_0>
 *   <v_0 vDF_1>,      <v_1 vDF_1>,      ..., <v_{nv} vDF_1>
 *                           ...
 *   <v_0 vDF_{nvDF}>, <v_0 vDF_{nvDF}>, ..., <v_{nv} vDF_{nvDF}>
 *   
 *   <e_0 eDF_0>,      <e_1 eDF_0>,      ..., <e_{ne} eDF_0>
 *   <e_0 eDF_1>,      <e_1 eDF_1>,      ..., <e_{ne} eDF_1>
 *                           ...
 *   <e_0 eDF_{neDF}>, <e_1 eDF_{neDF}>, ..., <e_{ne} eDF_{neDF}>
 *
 *   <f_0 fDF_0>,      <f_1 fDF_0>,      ..., <f_{nf} fDF_0>
 *   <f_0 fDF_1>,      <f_1 fDF_1>,      ..., <f_{nf} fDF_1>
 *                           ...
 *   <f_0 fDF_{nfDF}>, <f_1 fDF_{nfDF}>, ..., <f_{nf} fDF_{nfDF}>
 *
 *   <s_0 sDF_0>,      <s_1 sDF_0>,      ..., <s_{ns} sDF_0>
 *   <s_0 sDF_1>,      <s_1 sDF_1>,      ..., <s_{ns} sDF_1>
 *                           ...
 *   <s_0 sDF_{nsDF}>, <s_1 sDF_{nsDF}>, ..., <s_{ns} sDF_{nsDF}>
 *
 *   For example, linear elements in R^3, with one degree of freedom at each *
 *   vertex, would use the following ordering: 
 *
 *     <v_0 vDF_0>, <v_1 vDF_0>, <v_2 vDF_0>, <v_3 vDF_0> 
 *
 *   Quadratic elements in R^2, with one degree of freedom at each vertex and
 *   edge, would use the following ordering: 
 * 
 *     <v_0 vDF_0>, <v_1 vDF_0>, <v_2 vDF_0> 
 *     <e_0 eDF_0>, <e_1 eDF_0>, <e_2 eDF_0> 
 *
 *   You can use different trial and test spaces for each component of the
 *   elliptic system, thereby allowing for the use of Petrov-Galerkin methods.
 *   You MUST then tag the bilinear form symmetry entries as nonsymmetric in
 *   your PDE constructor to reflect that DF(u)(w,v) will be different from
 *   DF(u)(v,w), even if your form acts symmetrically when the same basis is
 *   used for w and v.
 * 
 *   You can also use different trial spaces for each component of the elliptic
 *   system, and different test spaces for each component of the elliptic
 *   system.  This allows you to e.g.  use a basis which is vertex-based for 
 *   one component, and a basis which is edge-based for another.  This is
 *   useful in fluid mechanics, eletromagnetics, or simply to play around with
 *   different elements.  
 *   
 *   This function is called by MC to build new master elements whenever it
 *   reads in a new mesh.  Therefore, this function does not have to be all
 *   that fast, and e.g.  could involve symbolic computation.
 *   @endverbatim
 */
VEXTERNC int Vfetk_PDE_simplexBasisInit(int key, int dim, int comp, int *ndof, 
  int dof[]);

/**
 * @brief Evaluate the bases for the trial or test space, for a particular
 * component of the system, at all quadrature points on the master simplex
 * element.
 * @ingroup Vfetk
 * @author Mike Holst
 * @param key basis type to evaluate
 *        \li 0 => trial
 *        \li 1 => test
 *        \li 2 => trialB
 *        \li 3 => testB
 * @param dim  spatial dimension
 * @param comp which component of elliptic system to produce basis for
 * @param pdkey basis partial differential equation evaluation key:
 *        \li 0 = evaluate basis(x,y,z)
 *        \li 1 = evaluate basis_x(x,y,z)
 *        \li 2 = evaluate basis_y(x,y,z)
 *        \li 3 = evaluate basis_z(x,y,z)
 *        \li 4 = evaluate basis_xx(x,y,z)
 *        \li 5 = evaluate basis_yy(x,y,z)
 *        \li 6 = evaluate basis_zz(x,y,z)
 *        \li 7 = etc...
 * @param xq  set to quadrature point coordinate
 * @param basis set to all basis functions evaluated at all quadrature points
 */
VEXTERNC void Vfetk_PDE_simplexBasisForm(int key, int dim, int comp,
    int pdkey, double xq[], double basis[]);

/**
 * @brief  Read in mesh and initialize associated internal structures
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @param  thee  Vfetk object
 * @param  sock  Socket ready for reading
 * @param  skey  Format key (0 => simplex format)
 * @note  @see Vfetk_genCube */
VEXTERNC void Vfetk_readMesh(Vfetk *thee, int skey, Vio *sock);

/**
 * @brief  Debugging routine to print out local variables used by PDE object
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @bug  This function is not thread-safe */
VEXTERNC void Vfetk_dumpLocalVar();

/**
 * @brief  Fill an array with the specified data
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @param  thee  Vfetk object
 * @param  vec  FEtk Bvec vector to use
 * @param  type  Type of data to write
 * @note  This function is thread-safe 
 * @bug  Several values of type are not implemented
 * @return  1 if successful, 0 otherwise */  
VEXTERNC int Vfetk_fillArray(Vfetk *thee, Bvec *vec, Vdata_Type type);

/**
 * @brief  Write out data
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @param  thee  Vfetk object
 * @param  vec  FEtk Bvec vector to use
 * @param  format  Format for data
 * @param iodev  Output device type (FILE/BUFF/UNIX/INET)
 * @param iofmt  Output device format (ASCII/XDR)
 * @param thost  Output hostname (for sockets)
 * @param fname  Output FILE/BUFF/UNIX/INET name
 * @note  This function is thread-safe 
 * @bug  Some values of format are not implemented
 * @return  1 if successful, 0 otherwise */  
VEXTERNC int Vfetk_write(Vfetk *thee, const char *iodev, const char *iofmt,
  const char *thost, const char *fname, Bvec *vec, Vdata_Format format);


#endif /* ifndef _VFETK_H_ */
