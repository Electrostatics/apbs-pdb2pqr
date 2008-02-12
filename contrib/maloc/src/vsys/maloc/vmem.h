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
 * File:     vmem.h    < vmem.c >
 *
 * Purpose:  Class Vmem: A safer, object-oriented, malloc/free object.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */

#ifndef _VMEM_H_
#define _VMEM_H_

#include <maloc/maloc_base.h>

/*
 * ***************************************************************************
 * Class Vmem: Parameters and datatypes
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * Class Vmem: Definition
 * ***************************************************************************
 */

typedef struct Vmem {

    char name[80];      /* name of class we are managing malloc areas for   */

    size_t mallocBytes;    /* total size of all current malloc areas of class  */
    size_t freeBytes;      /* total size of all freed malloc areas of class    */
    size_t highWater;      /* high-water malloc bytemark for this class        */
    size_t mallocAreas;    /* total number of individual malloc areas          */

} Vmem;

/*
 * ***************************************************************************
 * Class Vmem: Inlineable methods (vmem.c)
 * ***************************************************************************
 */

#if !defined(VINLINE_MALOC)
#else /* if defined(VINLINE_MALOC) */
#endif /* if !defined(VINLINE_MALOC) */

/*
 * ***************************************************************************
 * Class Vmem: Non-Inlineable methods (vmem.c)
 * ***************************************************************************
 */

VEXTERNC size_t Vmem_bytesTotal(void);
VEXTERNC size_t Vmem_mallocBytesTotal(void);
VEXTERNC size_t Vmem_freeBytesTotal(void);
VEXTERNC size_t  Vmem_highWaterTotal(void);
VEXTERNC size_t Vmem_mallocAreasTotal(void);
VEXTERNC void Vmem_prsize_tTotal(void);

VEXTERNC Vmem* Vmem_ctor(char *name);
VEXTERNC void Vmem_dtor(Vmem **thee);

VEXTERNC void *Vmem_malloc(Vmem *thee, size_t num, size_t size);
VEXTERNC void Vmem_free(Vmem *thee, size_t num, size_t size, void **ram);
VEXTERNC void *Vmem_realloc(Vmem *thee, size_t num, size_t size, void **ram,
    size_t newNum);

VEXTERNC size_t Vmem_bytes(Vmem *thee);
VEXTERNC size_t Vmem_mallocBytes(Vmem *thee);
VEXTERNC size_t Vmem_freeBytes(Vmem *thee);
VEXTERNC size_t Vmem_highWater(Vmem *thee);
VEXTERNC size_t Vmem_mallocAreas(Vmem *thee);
VEXTERNC void Vmem_prsize_t(Vmem *thee);

#endif /* _VMEM_H_ */

