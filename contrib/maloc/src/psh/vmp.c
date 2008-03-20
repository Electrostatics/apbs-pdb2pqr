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
 * rcsid="$Id: vmp.c,v 1.10 2008/03/12 05:13:58 fetk Exp $"
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * File:     vmp.c
 *
 * Purpose:  Class Vmp: methods.
 *
 * Notes:    Thin MPI abstraction layer on top of VCOM and VMPI.
 *           This layer is going to dissappear completely when
 *           VCOM and VMPI are merged.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */

#include "vmp_p.h"

VEMBED(rcsid="$Id: vmp.c,v 1.10 2008/03/12 05:13:58 fetk Exp $")

#define USE_VCOM_NOT 1

/*
 * ***************************************************************************
 * Class Vmp: Inlineable methods
 * ***************************************************************************
 */
#if !defined(VINLINE_MALOC)

#endif /* if !defined(VINLINE_MALOC) */
/*
 * ***************************************************************************
 * Class Vmp: Non-inlineable methods
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * Routine:  Vmp_init
 *
 * Purpose:  The Vmp initializer.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vmp_init(int *argc, char ***argv)
{
#   if defined(USE_VCOM)
        return Vcom_init(argc,argv);
#   else
        return Vmpi_init(argc,argv);
#   endif
}

/*
 * ***************************************************************************
 * Routine:  Vmp_finalize
 *
 * Purpose:  The Vmp finalizer.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vmp_finalize(void)
{
#   if defined(USE_VCOM)
        return Vcom_finalize();
#   else
        return Vmpi_finalize();
#   endif
}

/*
 * ***************************************************************************
 * Routine:  Vmp_ctor
 *
 * Purpose:  The Vmp constructor.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC Vmp* Vmp_ctor(void)
{
#   if defined(USE_VCOM)
        return (Vmp*)Vcom_ctor(1);
#   else
        return (Vmp*)Vmpi_ctor();
#   endif
}

/*
 * ***************************************************************************
 * Routine:  Vmp_dtor
 *
 * Purpose:  The Vmp destructor.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void Vmp_dtor(Vmp **thee)
{
#   if defined(USE_VCOM)
        Vcom_dtor( (Vcom**)thee );
#   else
        Vmpi_dtor( (Vmpi**)thee );
#   endif
}

/*
 * ***************************************************************************
 * Routine:  Vmp_rank
 *
 * Purpose:  Return my processor ID.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vmp_rank(Vmp *thee)
{
#   if defined(USE_VCOM)
        return Vcom_rank( (Vcom*)thee );
#   else
        return Vmpi_rank( (Vmpi*)thee );
#   endif
}

/*
 * ***************************************************************************
 * Routine:  Vmp_size
 *
 * Purpose:  Return the number of processors involved.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vmp_size(Vmp *thee)
{
#   if defined(USE_VCOM)
        return Vcom_size( (Vcom*)thee );
#   else
        return Vmpi_size( (Vmpi*)thee );
#   endif
}

/*
 * ***************************************************************************
 * Routine:  Vmp_barr
 *
 * Purpose:  An MPI barrier.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vmp_barr(Vmp *thee)
{
#   if defined(USE_VCOM)
        return Vcom_barr( (Vcom*)thee );
#   else
        return Vmpi_barr( (Vmpi*)thee );
#   endif
}

/*
 * ***************************************************************************
 * Routine:  Vmp_send
 *
 * Purpose:  An MPI blocking send.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vmp_send(Vmp *thee, int des, char *buf, int bufsize)
{
#   if defined(USE_VCOM)
        return Vcom_send( (Vcom*)thee, des, buf, bufsize, 0, 1 );
#   else
        return Vmpi_send( (Vmpi*)thee, des, buf, bufsize );
#   endif
}

/*
 * ***************************************************************************
 * Routine:  Vmp_recv
 *
 * Purpose:  An MPI blocking receive.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vmp_recv(Vmp *thee, int src, char *buf, int bufsize)
{
#   if defined(USE_VCOM)
        return Vcom_recv( (Vcom*)thee, src, buf, bufsize, 0, 1 );
#   else
        return Vmpi_recv( (Vmpi*)thee, src, buf, bufsize );
#   endif
}

