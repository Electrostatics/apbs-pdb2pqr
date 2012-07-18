/*
 * ***************************************************************************
 * MALOC = < Minimal Abstraction Layer for Object-oriented C >
 * Copyright (C) 1994--2008 Michael Holst
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * rcsid="$Id: vpup.h,v 1.13 2008/03/12 05:13:58 fetk Exp $"
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * File:     vpup.h    < vpup.c >
 *
 * Purpose:  Class Vpup: provides shell redirection, pipes, and execv/execvp.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */

#ifndef _VPUP_H_
#define _VPUP_H_

#include <maloc/maloc_base.h>
#include "maloccf.h"

/*
 * ***************************************************************************
 * Class Vpup: Parameters and datatypes
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * Class Vpup: Definition
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * Class Vpup: Inlineable methods (vpup.c)
 * ***************************************************************************
 */

#if !defined(VINLINE_MALOC)
#else /* if defined(VINLINE_MALOC) */
#endif /* if !defined(VINLINE_MALOC) */

/*
 * ***************************************************************************
 * Class Vpup: Non-Inlineable methods (vpup.c)
 * ***************************************************************************
 */

VEXTERNC void Vpup_execCmd(const char *PR, int argc, char **argv, char *inbuf);

#endif /* _VPUP_H_ */

