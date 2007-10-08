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
 * File:     vcom_p.h
 *
 * Purpose:  PRIVATE header.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */

#ifndef _VCOM_P_H_
#define _VCOM_P_H_

#include <maloc/vcom.h>
#include "maloccf.h"

#if defined(HAVE_MPI_H)
#   include <mpi.h>
#endif

typedef struct Vcom_core {
#if defined(HAVE_MPI_H)
    MPI_Status  mpi_status;       /* MPI_Status object   (4-int struct)     */
    MPI_Request mpi_request;      /* MPI_Request object  (union of structs) */
    MPI_Comm    mpi_comm;         /* MPI_Comm object     (int)              */
#else
    int mpi_status;
    int mpi_request;
    int mpi_comm;
#endif
} Vcom_core;

#endif /* _VCOM_P_H_ */

