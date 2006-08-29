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
 * rcsid="$Id: mainc.c,v 1.10 2006/05/08 04:39:44 mholst Exp $"
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * File:     main.c
 *
 * Purpose:  Test the C interface to the VIO library.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */

#include <maloc/maloc.h>

#define VEMBED(rctag) VPRIVATE const char* rctag; \
    static void* use_rcsid=(0 ? &use_rcsid : (void*)&rcsid);
VEMBED(rcsid="$Id: mainc.c,v 1.10 2006/05/08 04:39:44 mholst Exp $")

extern void outputMyStuff(Vio *sock);

int main(int argc, char *argv[])
{
    Vio *sock;

    char *oskey = "INET";      /* FILE, PIPE, UNIX, INET, MPI1 */
    char *osfmt = "ASC";       /* ASC, XDR */
    char *osnam = "1";         /* socket number or filename */
    char *ohval = "localhost"; /* hostname to send things to */

    /* open and connect to the socket */
    Vio_start();
    VJMPERR1( VNULL != (sock=Vio_ctor(oskey,osfmt,ohval,osnam,"w")) );
    VJMPERR1( 0 <= Vio_connect(sock,0) );

    /* write some data to the socket */
    outputMyStuff(sock);

    /* disconnect and close the socket */
    Vio_connectFree(sock);
    Vio_dtor(&sock);
    Vio_stop();

    return 1;

  VERROR1:
    fprintf(stderr,"Problem occurred.\n");
    return 0;
}

void outputMyStuff(Vio *sock)
{
    Vio_printf(sock,"%s\n","bhsingle");

    Vio_printf(sock,"%s %d\n","putl",4);
    Vio_printf(sock,"%e %e %e %e\n",0.0,1.0,0.0,0.0);
    Vio_printf(sock,"%e %e %e %e\n",0.0,0.0,1.0,0.0);
    Vio_printf(sock,"%e %e %e %e\n",0.0,0.0,0.0,1.0);

    Vio_printf(sock,"%s %d\n","list",7);
    Vio_printf(sock,"%s %d %d\n","fill",3,3);
    Vio_printf(sock,"%e %e %e\n",0.1,0.9,0.1);
    Vio_printf(sock,"%e %e %e\n",0.1,0.1,0.9);
    Vio_printf(sock,"%e %e %e\n",0.1,0.1,0.1);
    Vio_printf(sock,"%s %d %d\n","line",1,4);
    Vio_printf(sock,"%e %e %e\n",0.1,0.9,0.1);
    Vio_printf(sock,"%e %e %e\n",0.1,0.1,0.9);
    Vio_printf(sock,"%e %e %e\n",0.1,0.1,0.1);
    Vio_printf(sock,"%e %e %e\n",0.1,0.9,0.1);
    Vio_printf(sock,"%s %d\n","list",-7);

    Vio_printf(sock,"%s %d\n","list",2);
    Vio_printf(sock,"%s %d %d\n","fill",2,3);
    Vio_printf(sock,"%e %e %e\n",0.9,0.9,0.1);
    Vio_printf(sock,"%e %e %e\n",0.1,0.9,0.9);
    Vio_printf(sock,"%e %e %e\n",0.1,0.1,0.1);
    Vio_printf(sock,"%s %d %d\n","line",1,4);
    Vio_printf(sock,"%e %e %e\n",0.9,0.9,0.1);
    Vio_printf(sock,"%e %e %e\n",0.1,0.9,0.9);
    Vio_printf(sock,"%e %e %e\n",0.1,0.1,0.1);
    Vio_printf(sock,"%e %e %e\n",0.9,0.9,0.1);
    Vio_printf(sock,"%s %d\n","list",-2);

    Vio_printf(sock,"%s %d\n","list",3);
    Vio_printf(sock,"%s %d %d\n","fill",4,3);
    Vio_printf(sock,"%e %e %e\n",0.9,0.9,0.1);
    Vio_printf(sock,"%e %e %e\n",0.1,0.9,0.9);
    Vio_printf(sock,"%e %e %e\n",0.1,0.1,0.1);
    Vio_printf(sock,"%s %d %d\n","line",1,4);
    Vio_printf(sock,"%e %e %e\n",0.9,0.9,0.1);
    Vio_printf(sock,"%e %e %e\n",0.1,0.9,0.9);
    Vio_printf(sock,"%e %e %e\n",0.1,0.1,0.1);
    Vio_printf(sock,"%e %e %e\n",0.9,0.9,0.1);
    Vio_printf(sock,"%s %d\n","list",-3);

    Vio_printf(sock,"%s %d\n","putl",-1);
}

