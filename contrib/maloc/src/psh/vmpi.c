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
 * rcsid="$Id: vmpi.c,v 1.21 2008/03/12 05:13:58 fetk Exp $"
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * File:     vmpi.c
 *
 * Purpose:  Class Vmpi: methods.
 *
 * Notes:    Class Vmpi is a thin object-oriented Clean C layer on top of the
 *           MPI communication library.  Vmpi provides access to the minimal
 *           set of ten MPI primitives required to implement the Bank-Holst
 *           parallel adaptive algorithm, using either the Bank-Holst Oracle
 *           library, or directly.
 *
 *           Vmpi is also just about the smallest possible communication
 *           library that still contains just enough functionality to solve
 *           an elliptic equation in parallel using a typical parallel
 *           algorithm such as domain decomposition or parallel multigrid.
 *
 *           This minimal functionality is provided by Vmpi in ten core
 *           methods: startup, shutdown, size, rank, send, receive,
 *           broadcast, reduce, barrier, and non-blocking send (sometimes
 *           very useful).  That is it, and it all fits in about 90 lines
 *           of object-oriented Clean C code (not counting comments...).
 *           To fit into such a small amount of code, Vmpi assumes that the
 *           user does all of his own datapacking into pure character arrays.
 *           The best way to make use of Vmpi is to pack every message into
 *           a character buffer in XDR format (RPC's answer to
 *           machine-independent binary format).  You can accomplish this
 *           e.g. using the Vio library that sits at the bottom of MC.
 *
 *           The 6 primary MPI routines (according to Ian Foster):
 *           -----------------------------------------------------
 *
 *           MPI_Init       : Initiate an MPI computation
 *           MPI_Finalize   : Terminate a computation
 *           MPI_Comm_size  : Determine number of processes
 *           MPI_Comm_rank  : Determine my process identifier
 *           MPI_Send       : Send a message (blocking)
 *           MPI_Recv       : Receive a message (blocking)
 *
 *           The 4 additional useful MPI routines (according to me):
 *           -------------------------------------------------------
 *
 *           MPI_Barrier    : Synchronizes all processes
 *           MPI_Bcast      : Sends data from one process to all processes
 *           MPI_Reduce     : Sums, etc., distributed data (result in root)
 *           MPI_Isend      : Send a message (non-blocking)
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */

#include "vmpi_p.h"

#if defined(HAVE_MPI_H)
#   include <mpi.h>
#endif

VEMBED(rcsid="$Id: vmpi.c,v 1.21 2008/03/12 05:13:58 fetk Exp $")

/*
 * ***************************************************************************
 * Class Vmpi: Inlineable methods
 * ***************************************************************************
 */
#if !defined(VINLINE_MALOC)

#endif /* if !defined(VINLINE_MALOC) */
/*
 * ***************************************************************************
 * Class Vmpi: Non-inlineable methods
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * Routine:  Vmpi_init
 *
 * Purpose:  The Vmp initializer.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vmpi_init(int *argc, char ***argv)
{
#if defined(HAVE_MPI_H)
    return (MPI_SUCCESS == MPI_Init(argc,argv));
#else
    return 1;
#endif
}

/*
 * ***************************************************************************
 * Routine:  Vmpi_finalize
 *
 * Purpose:  The Vmp finalizerr.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vmpi_finalize(void)
{
#if defined(HAVE_MPI_H)
    return (MPI_SUCCESS == MPI_Finalize());
#else
    return 1;
#endif
}

/*
 * ***************************************************************************
 * Routine:  Vmpi_ctor
 *
 * Purpose:  The Vmpi constructor.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC Vmpi* Vmpi_ctor(void)
{
#if defined(HAVE_MPI_H)
    int dummy;
#endif
    Vmpi *thee = VNULL;
    VDEBUGIO("Vmpi_ctor: CREATING object..");
    thee = Vmem_malloc( VNULL, 1, sizeof(Vmpi) );
    thee->mpi_rank = 0;
    thee->mpi_size = 0;
#if defined(HAVE_MPI_H)
    VASSERT( MPI_SUCCESS == MPI_Initialized(&dummy) );
    VASSERT( MPI_SUCCESS == MPI_Comm_rank(MPI_COMM_WORLD, &(thee->mpi_rank)) );
    VASSERT( MPI_SUCCESS == MPI_Comm_size(MPI_COMM_WORLD, &(thee->mpi_size)) );
    Vnm_setIoTag(thee->mpi_rank, thee->mpi_size);
    Vnm_print(2,"Vmpi_ctor: process %d of %d is ALIVE!\n",
        thee->mpi_rank, thee->mpi_size);
#endif
    VDEBUGIO("..done.\n");
    return thee;
}

/*
 * ***************************************************************************
 * Routine:  Vmpi_dtor
 *
 * Purpose:  The Vmpi destructor.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void Vmpi_dtor(Vmpi **thee)
{
    VDEBUGIO("Vmpi_dtor: DESTROYING object..");
    Vmem_free( VNULL, 1, sizeof(Vmpi), (void**)thee );
#if defined(HAVE_MPI_H)
#if 0
    VASSERT( MPI_SUCCESS == MPI_Finalize() );
#endif
#endif
    VDEBUGIO("..done.\n");
}

/*
 * ***************************************************************************
 * Routine:  Vmpi_rank
 *
 * Purpose:  Return my processor ID.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vmpi_rank(Vmpi *thee)
{
    return thee->mpi_rank;
}

/*
 * ***************************************************************************
 * Routine:  Vmpi_size
 *
 * Purpose:  Return the number of processors involved.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vmpi_size(Vmpi *thee)
{
    return thee->mpi_size;
}

/*
 * ***************************************************************************
 * Routine:  Vmpi_barr
 *
 * Purpose:  An MPI barrier.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vmpi_barr(Vmpi *thee)
{
    int rc=0;
#if defined(HAVE_MPI_H)
    rc = (MPI_SUCCESS == MPI_Barrier(MPI_COMM_WORLD));
#endif
    return rc;
}

/*
 * ***************************************************************************
 * Routine:  Vmpi_send
 *
 * Purpose:  An MPI blocking send.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vmpi_send(Vmpi *thee, int des, char *buf, int bufsize)
{
    int rc=0;
#if defined(HAVE_MPI_H)
    int tag=0;
    rc = (MPI_SUCCESS == MPI_Send(buf,bufsize,MPI_CHAR,des,tag,MPI_COMM_WORLD));
#endif
    return rc;
}

/*
 * ***************************************************************************
 * Routine:  Vmpi_recv
 *
 * Purpose:  An MPI blocking receive.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vmpi_recv(Vmpi *thee, int src, char *buf, int bufsize)
{
    int rc=0;
#if defined(HAVE_MPI_H)
    int rsrc=0;
    int tag=0;
    MPI_Status stat;
    if (src < 0) {
        rsrc = MPI_ANY_SOURCE;
    } else {
        rsrc = src;
    }
    rc = (MPI_SUCCESS == MPI_Recv(buf,bufsize,MPI_CHAR,rsrc,tag,
        MPI_COMM_WORLD,&stat));
#endif
    return rc;
}

/*
 * ***************************************************************************
 * Routine:  Vmpi_bcast
 *
 * Purpose:  An MPI broadcast.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vmpi_bcast(Vmpi *thee, char *buf, int bufsize)
{
    int rc=0;
#if defined(HAVE_MPI_H)
    rc = (MPI_SUCCESS == MPI_Bcast(buf,bufsize,MPI_CHAR,0,MPI_COMM_WORLD));
#endif
    return rc;
}

/*
 * ***************************************************************************
 * Routine:  Vmpi_reduce
 *
 * Purpose:  An MPI reduce.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vmpi_reduce(Vmpi *thee, char *sbuf, char *rbuf, int bufsize)
{
    int rc=0;
#if defined(HAVE_MPI_H)
    rc = (MPI_SUCCESS ==
        MPI_Reduce(sbuf,rbuf,bufsize,MPI_CHAR,MPI_SUM,0,MPI_COMM_WORLD));
#endif
    return rc;
}

/*
 * ***************************************************************************
 * Routine:  Vmpi_isend
 *
 * Purpose:  An MPI non-blocking send.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vmpi_isend(Vmpi *thee, int des, char *buf, int bufsize)
{
    int rc=0;
#if defined(HAVE_MPI_H)
    int tag=0;
    MPI_Request requ;
    rc = (MPI_SUCCESS == MPI_Isend(buf,bufsize,MPI_CHAR,des,tag,MPI_COMM_WORLD,
        &requ));
#endif
    return rc;
}

