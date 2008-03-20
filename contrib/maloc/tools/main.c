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
 * rcsid="$Id: main.c,v 1.12 2008/02/22 23:06:00 fetk Exp $"
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * File:     bridge.c
 *
 * Purpose:  Socket-to-socket bridge.
 *
 * Notes:    The source and destination sockets can be any combination of
 *           UNIX domain (AF_UNIX, AF_LOCAL) and internet (AF_INET) sockets.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */

#include "maloc/maloc.h"

#define VEMBED(rctag) VPRIVATE const char* rctag; \
    static void* use_rcsid=(0 ? &use_rcsid : (void*)&rcsid);
VEMBED(rcsid="$Id: main.c,v 1.12 2008/02/22 23:06:00 fetk Exp $")

/* int MAIN__(void) {}; */

int main(int argc, char *argv[])
{
    char inType[80], outType[80], hostname[64], buf[VMAX_BUFSIZE];
    int  ok, inputWasFile, done, firstTime, bufsize;
    Vio   *inSock, *outSock;

    inputWasFile = 0;
    ok = 1;
    if (argc == 4) {
        strncpy(hostname,"localhost",sizeof(hostname));
        fprintf(stderr,"%s: Starting up on host <%s>\n", argv[0], hostname);
        if (!strcmp(argv[1],"-u2i")) {
            strcpy(inType,"UNIX");
            strcpy(outType,"INET");
        } else if (!strcmp(argv[1],"-i2u")) {
            strcpy(inType,"INET");
            strcpy(outType,"UNIX");
        } else if (!strcmp(argv[1],"-u2u")) {
            strcpy(inType,"UNIX");
            strcpy(outType,"UNIX");
        } else if (!strcmp(argv[1],"-i2i")) {
            strcpy(inType,"INET");
            strcpy(outType,"INET");
        } else if (!strcmp(argv[1],"-f2u")) {
            strcpy(inType,"FILE");
            strcpy(outType,"UNIX");
            inputWasFile = 1;
        } else if (!strcmp(argv[1],"-f2i")) {
            strcpy(inType,"FILE");
            strcpy(outType,"INET");
            inputWasFile = 1;
        } else ok = 0;
    } else ok = 0;

    VJMPERR1( ok );

    /* start the Vio communication layer */
    Vio_start();

    /* open the input socket */
    if ( VNULL==(inSock=Vio_ctor(inType,"ASC",hostname,argv[2],"r")) ) {
        fprintf(stderr,"%s: Problem open input socket <%s>\n",
            argv[0],argv[2]);
        VJMPERR1(0);
    } else {
        firstTime = 1;
        done = 0;
        while (!done) {

            /* do some extra stuff the first time */
            if (firstTime) {
                firstTime = 0;

                /* if input was a file, we just do one shot */
                if (inputWasFile) done = 1;

            /* must be non-file input; go to sleep for a bit */
            } else {
                /* (don't need this since now do BLOCKing accept/connect) */
                /* Vsig_sleep(100); */
            }

            /* call (blocking) accept on input socket */
            if ( 0 > Vio_accept(inSock,0) ) {
                fprintf(stderr,"%s: Problem accept input socket <%s>\n",
                    argv[0], argv[3]);
                VJMPERR1(0);
            } else {

                /* some i/o */
                fprintf(stderr,"%s: Hit on %s port <%s> from <%s>\n",
                    argv[0],inType,inSock->file,inSock->rhost);

                /* now that we have some data, open the output socket */
                if ( VNULL==(outSock=Vio_ctor(outType,"ASC",
                  hostname,argv[3],"w")) ) {
                    fprintf(stderr,"%s: Problem open output socket <%s>\n",
                        argv[0], argv[3]);
                    VJMPERR1(0);
                } else {

                    /* call (blocking) connect on output socket */
                    if ( 0 > Vio_connect(outSock,0) ) {
                        fprintf(stderr,"%s: Problem conn output socket <%s>\n",
                            argv[0], argv[3]);
                        VJMPERR1(0);
                    } else {

                        /* some i/o */
                        fprintf(stderr,"%s: Pass to %s port <%s> on <%s>\n",
                            argv[0],outType,outSock->file,outSock->lhost);

                        /* grab input and write to output until done */
                        memset(buf, '\0', sizeof(buf));
                        while (0<(bufsize=Vio_read(inSock,buf,sizeof(buf)))) {
                            VJMPERR2(bufsize==Vio_write(outSock,buf,bufsize));
                            memset(buf, '\0', sizeof(buf));
                        }

                        /* release the input subsocket and the output socket */
                      VERROR2:
                        Vio_connectFree(outSock);
                    }
                    /* close output socket */
                    Vio_dtor(&outSock);
                }
                /* free the input subsocket */
                Vio_acceptFree(inSock);
            }
        }
    }
    /* close input socket */
    Vio_dtor(&inSock);

    /* stop the Vio communication layer */
    Vio_stop();

    /* return with no errors */
    fprintf(stderr,"%s: Shutting down on host <%s>\n", argv[0], hostname);
    return 1;

  VERROR1:
    fprintf(stderr,"Usage: %s -u2i unixSocket inetSocket\n", argv[0]);
    fprintf(stderr,"   or: %s -i2u inetSocket unixSocket\n", argv[0]);
    fprintf(stderr,"   or: %s -u2u unixSocket unixSocket\n", argv[0]);
    fprintf(stderr,"   or: %s -i2i inetSocket inetSocket\n", argv[0]);
    fprintf(stderr,"   or: %s -f2u fileName   unixSocket\n", argv[0]);
    fprintf(stderr,"   or: %s -f2i fileName   inetSocket\n", argv[0]);
    return 0;
}

