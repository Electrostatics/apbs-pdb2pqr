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
 * File:     maloc_base.c
 *
 * Purpose:  MALOC linkage assistance.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */

#include "maloc_base_p.h"

VEMBED(rcsid="$Id$")

/*
 * ***************************************************************************
 * Routine:  maloc_link, maloc_needs_XXX
 *
 * Purpose:  Autoconf linkage assistance for packages built on top of MALOC.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void maloc_link(void)
{
}

#if (defined(HAVE_READLINE_READLINE_H) || defined(HAVE_READLINE_HISTORY_H))
    VPUBLIC void maloc_needs_rl(void)
    {
    }
#endif

#if defined(HAVE_MPI_H)
    VPUBLIC void maloc_needs_mpi(void)
    {
    }
#endif

