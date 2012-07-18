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
 * rcsid="$Id: maloc_base.c,v 1.16 2008/03/12 05:13:57 fetk Exp $"
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

VEMBED(rcsid="$Id: maloc_base.c,v 1.16 2008/03/12 05:13:57 fetk Exp $")

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

