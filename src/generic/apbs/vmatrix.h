/** @defgroup Vmatrix  Matrix wrapper class
 *  @brief    A header for including data wrapping matrices
 */

/**
 *  @file     vmatrix.h
 *  @ingroup  Vmatrix
 *  @brief    Contains inclusions for matrix data wrappers
 *  @version  @todo  figure out how the version identification works
 *  @author   Tucker A. Beck
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
 * Copyright (c) 2010-2011 Battelle Memorial Institute. Developed at the
 * Pacific Northwest National Laboratory,
 * operated by Battelle Memorial Institute,
 * Pacific Northwest Division for the U.S. Department Energy.
 * Portions Copyright (c) 2002-2010, Washington University in St. Louis.
 * Portions Copyright (c) 2002-2010, Nathan A. Baker.
 * Portions Copyright (c) 1999-2002, The Regents of the University of California
 * Portions Copyright (c) 1995, Michael Holst.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * Neither the name of Washington University in St. Louis nor the names of its
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
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @endverbatim
 */

#ifndef _VMATRIX_H_
#define _VMATRIX_H_

/* Generic headers */
#include "maloc/maloc.h"

#define MAT2(mat, rows, cols) \
    int rows_##mat = rows;    \
    int cols_##mat = rows

#define RAT2(mat, i, j) \
    (mat + ((j - 1) * rows_##mat + (i - 1)))

#define VAT2(mat, i, j) \
    mat[(j - 1) * rows_##mat + (i - 1)]

#define MAT3(mat, rows, cols, levs) \
    int rows_##mat = rows;          \
    int cols_##mat = cols;          \
    int levs_##mat = levs

#define RAT3(mat, i, j, k) \
    (mat + ((k - 1) * rows_##mat * cols_##mat + (j - 1) * rows_##mat + (i - 1)))

#define VAT3(mat, i, j, k) \
    mat[(k - 1) * rows_##mat * cols_##mat + (j - 1) * rows_##mat + (i - 1)]

#endif /* _VMATRIX_H_ */
