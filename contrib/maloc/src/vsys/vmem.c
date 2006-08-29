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
 * rcsid="$Id: vmem.c,v 1.12 2006/07/28 02:11:50 mholst Exp $"
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * File:     vmem.c
 *
 * Purpose:  Class Vmem: methods.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */

#include "vmem_p.h"

VEMBED(rcsid="$Id: vmem.c,v 1.12 2006/07/28 02:11:50 mholst Exp $")

/* total and misc (default) malloc/free tracking */
VPRIVATE Vmem vmemTotal, vmemMisc;
VPRIVATE int vmemInit=0;

/*
 * ***************************************************************************
 * Class Vmem: Inlineable methods
 * ***************************************************************************
 */
#if !defined(VINLINE_MALOC)

#endif /* if !defined(VINLINE_MALOC) */
/*
 * ***************************************************************************
 * Class Vmem: Non-inlineable methods
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * Routine:  Vmem_init
 *
 * Purpose:  Initialize the total log and the catch-all misc log.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPRIVATE void Vmem_init(void)
{
    if (!vmemInit) {

        strncpy(vmemTotal.name, "TOTAL", 80);
        vmemTotal.mallocBytes = 0;
        vmemTotal.freeBytes = 0;
        vmemTotal.highWater = 0;
        vmemTotal.mallocAreas = 0;

        strncpy(vmemMisc.name, "MISC", 80);
        vmemMisc.mallocBytes = 0;
        vmemMisc.freeBytes = 0;
        vmemMisc.highWater = 0;
        vmemMisc.mallocAreas = 0;

        vmemInit = 1;
    }
}

/*
 * ***************************************************************************
 * Routine:  Vmem_bytesTotal
 *
 * Purpose:  Return total size of ALL currently ACTIVE malloc areas that
 *           went through Vmem_malloc.  This is the current memory footprint.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vmem_bytesTotal(void)
{
    Vmem_init();
    return (vmemTotal.mallocBytes - vmemTotal.freeBytes);
}

/*
 * ***************************************************************************
 * Routine:  Vmem_mallocBytesTotal
 *
 * Purpose:  Return total size of ALL malloc areas that went through
 *           Vmem_malloc (even if they have been subsequently freed).
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vmem_mallocBytesTotal(void)
{
    Vmem_init();
    return vmemTotal.mallocBytes;
}

/*
 * ***************************************************************************
 * Routine:  Vmem_freeBytesTotal
 *
 * Purpose:  Return total size of ALL freed malloc areas that
 *           went through Vmem_free.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vmem_freeBytesTotal(void)
{
    Vmem_init();
    return vmemTotal.freeBytes;
}

/*
 * ***************************************************************************
 * Routine:  Vmem_highWaterTotal
 *
 * Purpose:  Return the high-water malloc bytemark hit by ALL active
 *           malloc areas; this is the largest active malloc total hit at
 *           one time by areas obtained through Vmem_malloc.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vmem_highWaterTotal(void)
{
    Vmem_init();
    return vmemTotal.highWater;
}

/*
 * ***************************************************************************
 * Routine:  Vmem_mallocAreasTotal
 *
 * Purpose:  Return the total number of ALL individual active malloc areas.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vmem_mallocAreasTotal(void)
{
    Vmem_init();
    return vmemTotal.mallocAreas;
}

/*
 * ***************************************************************************
 * Routine:  Vmem_printTotal
 *
 * Purpose:  Print current memory statistics for all malloc areas that have
 *           gone through Vmem_malloc and Vmem_free.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void Vmem_printTotal(void)
{
    Vmem_init();
    fprintf(stderr,"%12d %12d %12d %12d %12d %% %s\n",
        (vmemTotal.mallocBytes-vmemTotal.freeBytes),
        vmemTotal.mallocAreas, vmemTotal.mallocBytes,
        vmemTotal.freeBytes, vmemTotal.highWater,
        vmemTotal.name);
}

/*
 * ***************************************************************************
 * Routine:  Vmem_ctor
 *
 * Purpose:  Construct the dynamic memory allocation logging object.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC Vmem *Vmem_ctor(char *name)
{
    Vmem *thee;

    thee = Vmem_malloc( VNULL, 1, sizeof(Vmem) );
    VASSERT( thee != VNULL );

    strncpy( thee->name, name, 80 );
    thee->mallocBytes = 0;
    thee->freeBytes = 0;
    thee->highWater = 0;
    thee->mallocAreas = 0;

    return thee;
}

/*
 * ***************************************************************************
 * Routine:  Vmem_dtor
 *
 * Purpose:  Destruct the dynamic memory allocation logging object.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void Vmem_dtor(Vmem **thee)
{
    Vmem_free( VNULL, 1, sizeof(Vmem), (void**)thee );
}

/*
 * ***************************************************************************
 * Routine:  Vmem_malloc
 *
 * Purpose:  A logged version of malloc.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void *Vmem_malloc(Vmem *thee, int num, int size)
{
    int btmp;
    void *ram = VNULL;

    Vmem_init();

    /* VWARN( (num > 0) && (size > 0) ); */
    VASSERT( (num > 0) && (size > 0) );
    if ( (num > 0) && (size > 0) ) {

        ram = (void*)calloc((unsigned int)num, (unsigned int)size);
        VASSERT( ram != VNULL );

        vmemTotal.mallocBytes += (num * size);
        btmp = (vmemTotal.mallocBytes - vmemTotal.freeBytes);
        if ( vmemTotal.highWater < btmp ) vmemTotal.highWater = btmp;
        vmemTotal.mallocAreas++;

        if (thee != VNULL) {
            thee->mallocBytes += (num * size);
            btmp = (thee->mallocBytes - thee->freeBytes);
            if ( thee->highWater < btmp ) thee->highWater = btmp;
            thee->mallocAreas++;
        } else {
            vmemMisc.mallocBytes += (num * size);
            btmp = (vmemMisc.mallocBytes - vmemMisc.freeBytes);
            if ( vmemMisc.highWater < btmp ) vmemMisc.highWater = btmp;
            vmemMisc.mallocAreas++;
        }

    }

    return ram;
}

/*
 * ***************************************************************************
 * Routine:  Vmem_free
 *
 * Purpose:  A logged version of free.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void Vmem_free(Vmem *thee, int num, int size, void **ram)
{
    Vmem_init();

    /* VWARN( (*ram) != VNULL ); */
    VASSERT( (*ram) != VNULL );
    if ((*ram) != VNULL) {

        free(*ram);
        (*ram) = VNULL;

        vmemTotal.freeBytes += (num * size);
        vmemTotal.mallocAreas--;

        if (thee != VNULL) {
            thee->freeBytes += (num * size);
            thee->mallocAreas--;
        } else {
            vmemMisc.freeBytes += (num * size);
            vmemMisc.mallocAreas--;
        }
    }
}

/*
 * ***************************************************************************
 * Routine:  Vmem_realloc
 *
 * Purpose:  A logged version of realloc (using this is usually a bad idea).
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void *Vmem_realloc(Vmem *thee, int num, int size, void **ram,
    int newNum)
{
    void *tee = Vmem_malloc(thee, newNum, size);
    memcpy(tee, (*ram), size*VMIN2(num,newNum));
    Vmem_free(thee, num, size, ram);
    return tee;
}

/*
 * ***************************************************************************
 * Routine:  Vmem_bytes
 *
 * Purpose:  Return total of ACTIVE malloc areas used by the Vmem object.
 *           (If Vmem is VNULL, return the misc catch-all malloc total).
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vmem_bytes(Vmem *thee)
{
    Vmem_init();

    if (thee != VNULL) {
        return (thee->mallocBytes - thee->freeBytes);
    } else {
        return (vmemMisc.mallocBytes - vmemMisc.freeBytes);
    }
}

/*
 * ***************************************************************************
 * Routine:  Vmem_mallocBytes
 *
 * Purpose:  Return total of all mallocs performed by the Vmem object.
 *           (If Vmem is VNULL, return the misc catch-all malloc total).
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vmem_mallocBytes(Vmem *thee)
{
    Vmem_init();

    if (thee != VNULL) {
        return thee->mallocBytes;
    } else {
        return vmemMisc.mallocBytes;
    }
}

/*
 * ***************************************************************************
 * Routine:  Vmem_freeBytes
 *
 * Purpose:  Return total of the frees performed by the Vmem object.
 *           (If Vmem is VNULL, return misc catch-all free total).
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vmem_freeBytes(Vmem *thee)
{
    Vmem_init();

    if (thee != VNULL) {
        return thee->freeBytes;
    } else {
        return vmemMisc.freeBytes;
    }
}

/*
 * ***************************************************************************
 * Routine:  Vmem_highWater
 *
 * Purpose:  Return the high-water malloc bytemark hit by the Vmem object.
 *           (If Vmem is VNULL, return misc catch-all malloc highwater).
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vmem_highWater(Vmem *thee)
{
    Vmem_init();

    if (thee != VNULL) {
        return thee->highWater;
    } else {
        return vmemMisc.highWater;
    }
}

/*
 * ***************************************************************************
 * Routine:  Vmem_mallocAreas
 *
 * Purpose:  Return the total number of individual active malloc areas.
 *           (If Vmem is VNULL, return misc catch-all malloc areas).
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vmem_mallocAreas(Vmem *thee)
{
    Vmem_init();

    if (thee != VNULL) {
        return thee->mallocAreas;
    } else {
        return vmemMisc.mallocAreas;
    }
}

/*
 * ***************************************************************************
 * Routine:  Vmem_print
 *
 * Purpose:  Print current memory statistics associated with this Vmem object.
 *           (If Vmem is VNULL, print info for the catch-all ovject).
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void Vmem_print(Vmem *thee)
{
    Vmem_init();

    if (thee != VNULL) {
        fprintf(stderr,"%12d %12d %12d %12d %12d %% %s\n",
            (thee->mallocBytes-thee->freeBytes),
            thee->mallocAreas, thee->mallocBytes,
            thee->freeBytes, thee->highWater,
            thee->name);
    } else {
        fprintf(stderr,"%12d %12d %12d %12d %12d %% %s\n",
            (vmemMisc.mallocBytes-vmemMisc.freeBytes),
            vmemMisc.mallocAreas, vmemMisc.mallocBytes,
            vmemMisc.freeBytes, vmemMisc.highWater,
            vmemMisc.name);
    }
}

