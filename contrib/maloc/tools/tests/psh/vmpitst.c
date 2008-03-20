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
 * rcsid="$Id: vmpitst.c,v 1.13 2008/03/12 05:14:00 fetk Exp $"
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * File:     main.c
 *
 * Purpose:  Test main driver for the VMPI layer.
 *
 * Author:   Michael Holst 060197
 * ***************************************************************************
 */

#include <maloc/maloc.h>

#define VEMBED(rctag) VPRIVATE const char* rctag; \
    static void* use_rcsid=(0 ? &use_rcsid : (void*)&rcsid);
VEMBED(rcsid="$Id: vmpitst.c,v 1.13 2008/03/12 05:14:00 fetk Exp $")

int main(int argc, char *argv[])
{
    /*
     * *************************************************************
     * variables
     * **************************************************************
     */

    /* mpi variables */
    char buffer1, buffer2;
    int key;

    /* the mpi objects */
    Vmpi *vmpi = VNULL;

    /*
     * *************************************************************
     * mpi setup
     * **************************************************************
     */

    /* construct the vmpi object */
    VASSERT( Vmpi_init(&argc, &argv) );
    vmpi = Vmpi_ctor();

    /* setup -- root guy */
    if (Vmpi_rank(vmpi) == 0) {
        fprintf(stderr,"<Process #%d (of %d) -- STARTUP(ROOT)>\n",
            Vmpi_rank(vmpi), Vmpi_size(vmpi));
        buffer1 = 2;

    /* setup -- all non-root guys */
    } else {
        fprintf(stderr,"<Process #%d (of %d) -- STARTUP(DRONE)>\n",
            Vmpi_rank(vmpi), Vmpi_size(vmpi));
    }

    /*
     * *************************************************************
     * computations...
     * **************************************************************
     */

    /* send root guy's value to everyone; check to see that everyone got it */
    key = 1;
    Vmpi_bcast(vmpi,&buffer1,key);
    if (buffer1 != 2) fprintf(stderr,"Problem!\n");

    /* now add the broadcasted value up on every process */
    key = 1;
    buffer2 = buffer1;
    Vmpi_reduce(vmpi,&buffer2,&buffer1,key);

    /*
     * *************************************************************
     * mpi shutdown
     * **************************************************************
     */

    /* cleanup -- root guy */
    if (Vmpi_rank(vmpi) == 0) {
        /* default parameters */

        /* finish up */
        fprintf(stderr,"<Process #%d (of %d) -- SHUTDOWN(ROOT)>\n",
            Vmpi_rank(vmpi), Vmpi_size(vmpi));
        fprintf(stderr,"<RESULT = %d>\n", buffer1);

    /* cleanup -- all non-root guys */
    } else {
        fprintf(stderr,"<Process #%d (of %d) -- SHUTDOWN(DRONE)>\n",
            Vmpi_rank(vmpi), Vmpi_size(vmpi));
    }

    /* keep everyone waiting here until root is ready */
    Vmpi_barr(vmpi);

    /* destroy the vmpi object */
    Vmpi_dtor(&vmpi);
    VASSERT( Vmpi_finalize() );

    /* return */
    return 0;
}

