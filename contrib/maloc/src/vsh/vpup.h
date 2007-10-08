/*
 * ***************************************************************************
 * MALOC = < Minimal Abstraction Layer for Object-oriented C >
 * Copyright (C) 1994--2000  Michael Holst
 * 
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 675 Mass Ave, Cambridge, MA 02139, USA.
 * 
 * rcsid="$Id$"
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

