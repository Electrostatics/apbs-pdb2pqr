/** @defgroup Vmultigrid Multi-grid class
 *  @brief    Multi-grid class used to encapsulate values and methods for
 *            multi-grid computations
 */

/**
 *  @file    vmultigrid.h
 *  @ingroup Vmultigrid
 *  @author  Tucker Beck
 *  @brief
 *  @version
 *
 *  @attention
 *  @verbatim
 *
 * APBS -- Adaptive Poisson-Boltzmann Solver
 *
 * Nathan A. Baker (nathan.baker@pnl.gov)
 * Pacific Northwest National Laboratory
 *
 * Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 2010-2011 Battelle Memorial Institute. Developed at the Pacific Northwest National Laboratory, operated by Battelle Memorial Institute, Pacific Northwest Division for the U.S. Department Energy.  Portions Copyright (c) 2002-2010, Washington University in St. Louis.  Portions Copyright (c) 2002-2010, Nathan A. Baker.  Portions Copyright (c) 1999-2002, The Regents of the University of California. Portions Copyright (c) 1995, Michael Holst.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * -  Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * - Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * - Neither the name of Washington University in St. Louis nor the names of its
 * contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF   THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @endverbatim
 */

#ifndef _VMULTIGRID_H_
#define _VMULTIGRID_H_

#include "maloc/maloc.h"

#include "apbs/vhal.h"
#include "apbs/vsize3d.h"
#include "apbs/vmatrix3dfast.h"



/// Defines the method to be used for the computation
enum eVmultigrid_method {
    VMULTIGRID_CGMG,      ///< Conjugate Gradient Multi-grid
    VMULTIGRID_NEW,       ///< Newton's
    VMULTIGRID_MG,        ///< Standard Multigrid
    VMULTIGRID_CG,        ///< Conjugate Gradient
    VMULTIGRID_SOR,       ///< Successive Overrelaxation
    VMULTIGRID_RBGS,      ///< Red-bklack Gauss-Seidel
    VMULTIGRID_WJ,        ///< Weighted Jacobian
    VMULTIGRID_RICH,      ///< Richardson
    VMULTIGRID_CGMG_AQUA, ///< Conjugate Gradient Multigrid Aqua
    VMULTIGRID_NEW_AQUA   ///< Newton's Aqua
};

typedef enum eVmultigrid_method Vmultigrid_method;



/// Defines the boundary condition method
enum eVmultigrid_boundary_method {
    VMULTIGRID_BOUNDARY_ZERO,   ///< Zero Dirichlet
    VMULTIGRID_BOUNDARY_SDH,    ///< Single-sphere Debye-Huckel Dirichlet
    VMULTIGRID_BOUNDARY_MDH,    ///< Multiple-sphere Debye-Huckel Dirichlet
    VMULTIGRID_BOUNDARY_FOCUS,  ///< Focusing Dirichlet
    VMULTIGRID_BOUNDARY_MEM,    ///< Focusing membrane
    VMULTIGRID_BOUNDARY_MAP,    ///< Skip first level focusing. Use external map
    VMULTIGRID_BOUNDARY_UNUSED  ///< Unused mehtod (placeholder) @todo Remove
};

typedef enum eVmultigrid_boundary_method Vmultigrid_boundary_method;



/// Controls the Prolongation method
enum eVmultigrid_prolong_method {
    VMULTIGRID_PROLONG_TRILINEAR, ///< Trilinear prolongation
    VMULTIGRID_PROLONG_OPERATOR,  ///< Operator based prolongation
    VMULTIGRID_PROLONG_MODULO,    ///< Modulo operator based prolongation
};

typedef enum eVmultigrid_prolong_method Vmultigrid_prolong_method;



/// Controls the coarsening method
enum eVmultigrid_coarsen_method {
    VMULTIGRID_COARSEN_STANDARD, ///< Standard Discretization
    VMULTIGRID_COARSEN_AVERAGED, ///< Averaged Coef and Std Discretization
    VMULTIGRID_COARSEN_GALERKIN  ///< Algebraic Galerkin Coarsening
};

typedef enum eVmultigrid_coarsen_method Vmultigrid_coarsen_method;



/// Controls the discretization method
enum eVmultigrid_discrete_method {
    VMULTIGRID_DISCRETE_VOLUME, ///< Finite Volume
    VMULTIGRID_DISCRETE_ELEMENT ///< Finite Element Method
};

typedef enum eVmultigrid_discrete_method Vmultigrid_discrete_method;



/// Controls the smoothing method
enum eVmultigrid_smooth_method {
    VMULTIGRID_SMOOTH_JACOBIAN,   ///< Weighted Jacobian method
    VMULTIGRID_SMOOTH_GAUSS,      ///< Gauss-seidel method
    VMULTIGRID_SMOOTH_SOR,        ///< Successive overrelaxtion
    VMULTIGRID_SMOOTH_RICHARDSON, ///< Richardson's method
    VMULTIGRID_SMOOTH_CGHS        ///< @todo Find name: Conjugate Gradient?
};

typedef enum eVmultigrid_smooth_method Vmultigrid_smooth_method;



/// Controls the coarse grid solver
enum eVmultigrid_coarse_solver {
    VMULTIGRID_SOLVER_CG,       ///< Conjugate Gradients solver
    VMULTIGRID_SOLVER_SYMMETRIC ///< Smmetric banded linpack solver
};

typedef enum eVmultigrid_coarse_solver Vmultigrid_coarse_solver;



/// Controls the multigrid iteration method
enum eVmultigrid_iterate_method {
    VMULTIGRID_ITERATE_VARIABLE, ///< Variable v-cycle
    VMULTIGRID_ITERATE_NESTED,   ///< Nested iteration
};

typedef enum eVmultigrid_iterate_method Vmultigrid_iterate_method;



/// Controls run-time status message display
enum eVmultigrid_verbosity {
    VMULTIGRID_VERBOSE_NONE, ///< Print no status messages
    VMULTIGRID_VERBOSE_SOME, ///< Print some status messages
    VMULTIGRID_VERBOSE_MOST, ///< Print most status messages
    VMULTIGRID_VERBOSE_ALL,  ///< Print all status messages
};

typedef enum eVmultigrid_verbosity Vmultigrid_verbosity;




/// Controls the stopping criterion
enum eVmultigrid_stop_criterion {
    VMULTIGRID_STOP_RESIDUAL, ///< Residual
    VMULTIGRID_STOP_RELATIVE, ///< Relative Residual
    VMULTIGRID_STOP_DIFF,     ///< @todo Decipher and un-abbreviate
    VMULTIGRID_STOP_ERRC,     ///< @todo Decipher and un-abbreviate
    VMULTIGRID_STOP_ERRD,     ///< @todo Decipher and un-abbreviate
    VMULTIGRID_STOP_AERRD     ///< @todo Decipher and un-abbreviate
};

typedef enum eVmultigrid_stop_criterion Vmultigrid_stop_criterion;



/// Controls the use of non-linearity
enum eVmultigrid_linear {
    VMULTIGRID_LINEAR_SIZE,      ///< Size-modified PBE
    VMULTIGRID_LINEAR_NORMAL,    ///< @todo Document Me
    VMULTIGRID_LINEAR_NON,       ///< Nonlinear PBE with capped sinh term
    VMULTIGRID_LINEAR_POLYNOMIAL ///< Polynomial approximation to sinh
};

typedef enum eVmultigrid_linear Vmultigrid_linear;



/// Controls method for analysis of the operator
enum eVmultigrid_operator_analysis {
    VMULTIGRID_ANALYSIS_NONE,        ///< No analysis
    VMULTIGRID_ANALYSIS_CONDITIONAL, ///< Conditional number analysis
    VMULTIGRID_ANALYSIS_SPECTRAL,    ///< Spectral radius analysis
    VMULTIGRID_ANALYSIS_BOTH         ///< Use both previous analysis methods
};

typedef enum eVmultigrid_operator_analysis Vmultigrid_operator_analysis;



/**
*  @ingroup Vmultigrid
*  @author  Tucker Beck
*  @brief   Declaration of Vmultigrid class
*/
struct sVmultigrid {

    // Algorithmic control flags

    /// The computation method to be used
    Vmultigrid_method method;

    /// The method to be used for handling boundary conditions
    Vmultigrid_boundary_method boundary_method;

    /// The prolongation method to be used
    /// @note Replaces mgprol from mgdrvd.c
    Vmultigrid_prolong_method prolong_method;

    /// The coarsening method to be used
    /// @note Replaces mgcoar from mgdrvd.c
    Vmultigrid_coarsen_method coarsen_method;

    /// The discretization method to be used
    /// @note Replaces mgdisc from mgdrvd.c
    Vmultigrid_discrete_method discrete_method;

    /// The smoothing method to be used
    /// @note Replaces mgsmoo from mgdrvd.c
    Vmultigrid_smooth_method smooth_method;

    /// The coarse solver to be used
    /// @note Replaces mgsolv from mgdrvd.c
    Vmultigrid_coarse_solver coarse_solver;

    /// The multigrid iteration method to be used
    Vmultigrid_iterate_method iterate_method;

    /// The verbosity of status messages to display
    /// @note Replaces iinfo from mgdrvd.c
    Vmultigrid_verbosity verbosity;

    /// The stopping criterion to cease comutation iterations
    /// @note Replaces istop from mgdrvd.c
    Vmultigrid_stop_criterion stop_criterion;

    /// The linear or non-linear approximation to use
    Vmultigrid_linear linear;

    /// The operator analysis method to use
    Vmultigrid_operator_analysis operator_analysis;



    // Workspace boundary totals

    /// The upper bound on the size of the real work space
    /// @note Replaces nrwk from mgdrvd.c
    int real_workspace_limit;

    /// The upper bound on the size of the integer work space
    /// @note Replaces irwk from mgdrvd.c
    int integer_workspace_limit;

    /// The total size of the real work space
    /// @note Replaces iretot from mgdrvd.c
    int real_workspace_size;

    /// The total size of the integer work space
    /// @note Replaces iintot from mgdrvd.c
    int integer_workspace_size;

    
    
    // Grid dimensions
    
    /// The size of the grid in the x dimension
    /// @note Replaces nx from mgdrvd.h
    int x_size;
    
    /// The size of the grid in the y dimension
    /// @note Replaces ny from mgdrvd.h
    int y_size;
    
    /// The size of the grid in the z dimension
    /// @note Replaces nz from mgdrvd.h
    int z_size;
    
    /// The number of levels of the multi-grid
    /// @note Replaces nlev from mgdrvd.h
    int level_count;
    
    /// The number of unknowns on the finest mesh
    /// @note Replaces nf from mgdrvd.h
    int fine_count;
    
    /// The number of unknowns on the coarsest mesh
    /// @note Replaces nc from mgdrvd.h
    int coarse_count;
    

    

    // Iteration control

    /// The maximum number of iterations to perform
    int max_iterations;

    
    
    // Computation Matrices
    Vm


    // Undocumented
    // @todo  Finish investigation, refactor, and documentation

    /// @todo Document Me
    int ipcon;

    /// @todo Document Me
    int irite;

    /// @todo Document Me
    int nonlin;

};

/**
 *  @ingroup Vmultigrid
 *  @brief   Declaration of the Vmultigrid class as the Vmultigrid structure
 */
typedef struct sVmultigrid Vmultigrid;

#endif // _VMULTIGRID_H_

