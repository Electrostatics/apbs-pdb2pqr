/** @defgroup Vstring Vstring class
 *  @brief    Provides a collection of useful non-ANSI string functions
 */

/**
 *  @file     vstring.h
 *  @ingroup  Vstring
 *  @brief    Contains declarations for class Vstring
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

#ifndef _VSTRING_H_
#define _VSTRING_H_

#include "maloc/maloc.h"

/** @brief   Case-insensitive string comparison (BSD standard)
 *  @ingroup Vstring
 *  @author  Copyright (c) 1988-1993 The Regents of the University of
 *           California.  Copyright (c) 1995-1996 Sun Microsystems, Inc.
 *  @note    Copyright (c) 1988-1993 The Regents of the University of
 *           California.  Copyright (c) 1995-1996 Sun Microsystems, Inc.
 *  @param   s1   First string for comparison
 *  @param   s2   Second string for comparison
 *  @return  An integer less than, equal to, or greater than zero if s1 is
 *           found,  respectively,  to  be  less  than, to match, or be greater
 *           than s2. (Source:  Linux man pages)
 */
VEXTERNC int Vstring_strcasecmp(const char *s1, const char *s2);

#endif    /* ifndef _VSTRING_H_ */
