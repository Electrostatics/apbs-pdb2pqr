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
 * rcsid="$Id: vcomtst.c,v 1.12 2008/03/12 05:14:00 fetk Exp $"
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * File:     vcomtst.c
 *
 * Purpose:  Low-level test of MPI routines in vcom.c
 *
 * Author:   Nathan Baker and Michael Holst
 * ***************************************************************************
 */

#include <maloc/maloc.h>

#define VEMBED(rctag) VPRIVATE const char* rctag; \
    static void* use_rcsid=(0 ? &use_rcsid : (void*)&rcsid);
VEMBED(rcsid="$Id: vcomtst.c,v 1.12 2008/03/12 05:14:00 fetk Exp $")

int main(int argc, char *argv[]) {

    int i, testi, testj;
    int myrank;
    int nproc;
    char string[100];
    Vcom *vcom;

    VASSERT( Vcom_init(&argc,&argv) );
    vcom = Vcom_ctor(1);
    myrank = Vcom_rank(vcom);
    nproc = Vcom_size(vcom);

    printf("PE %d: Starting test program with %d total PEs.\n",myrank,nproc);

    /* Character send/recv test */
    if (myrank == 0) {
        printf("\n\nPE %d: STARTING SEND/RECV TEST.\n",myrank);
        printf("PE %d: Sending non-blocked messages.\n",myrank);
        for (i=1; i<nproc; i++) {
            sprintf(&(string[0]), "non-blocked message (PE %d --> %d)", 
              myrank, i);
            printf( "PE %d: Sent non-blocked message to PE %d\n", myrank, i);
            assert(Vcom_send(vcom,i,string,100,3,0));
        }
        printf("PE %d: Sending blocked messages.\n",myrank);
        for (i=1; i<nproc; i++) {
            sprintf(&(string[0]),"blocked message (PE %d --> %d)",i,nproc);
            printf( "PE %d: Sent blocked message to PE %d\n", myrank, i);
            assert(Vcom_send(vcom,i,string,100,3,1));
        }
    } else {
        printf("PE %d: Recving blocked message from PE 0.\n",myrank);
        /* Get blocked sent messages */
        assert(Vcom_recv(vcom,0,string,100,3,1));
        printf( "PE %d: Blocked message is: '%s'\n",myrank,string);
        /* Get blocked sent messages */
        printf("PE %d: Recving non-blocked message from PE 0.\n",myrank);
        assert(Vcom_recv(vcom,0,string,100,3,1));
        printf( "PE %d: Non-blocked message is: '%s'\n",myrank,string);
    }
    fflush(stdout);
    Vcom_barr(vcom);

    /* Barrier test */
    if (myrank == 0) printf("\n\nPE %d: STARTING BARRIER TEST.\n",myrank);
    printf("PE %d: Hit barrier\n", myrank);
    Vcom_barr(vcom);
    printf("PE %d: Passed barrier\n", myrank);
    fflush(stdout);
    Vcom_barr(vcom);

    /* getCount test */
    if (myrank == 0) {
        printf("\n\nPE %d: STARTING PROBE TEST.\n",myrank);
        printf("PE %d: Sending string of length 100\n", myrank);
        for (i=1; i<nproc; i++) assert(Vcom_send(vcom,i,string,100,3,0));
    } else {
        Vcom_getCount(vcom, 0, &testi, 3);
        printf("PE %d: Probed for string of length %d\n", myrank, testi);
    }
    fflush(stdout);
    Vcom_barr(vcom);
    
    /* Reduction test */
    testi = 4;
    if (myrank == 0) {
        printf("\n\nPE %d: STARTING REDUCTION TEST.\n",myrank);
        printf("PE %d: All %d PEs have the number %d.\n", myrank, nproc, testi);
        printf("PE %d: Performing global sum reduction.\n", myrank);
    }
    Vcom_reduce(vcom, &testi, &testj, 1, 1, 0);
    printf("PE %d: Global sum = %d\n", myrank, testj);
    fflush(stdout);
    Vcom_barr(vcom);
    testi = 4;
    if (myrank == 0) {
        testj = 8;
        printf("PE %d: PE %d has the number %d; all others have %d.\n", myrank,
          myrank, testj, testi);
        testi = testj;
        printf("PE %d: Performing global max reduction.\n", myrank);
    }
    Vcom_reduce(vcom, &testi, &testj, 1, 1, 3);
    printf("PE %d: Global max = %d\n", myrank, testj);
    fflush(stdout);
    Vcom_barr(vcom);

    /* Resize the communicator */
    printf("PE %d: Resizing communicator from %d to %d.\n", myrank, nproc,
      (int)(nproc/2));
    Vcom_resize(vcom, (int)(nproc/2));
    myrank = Vcom_rank(vcom);
    if (Vcom_rank(vcom) != -1) printf("PE %d: Hello world from resized communicator\n", myrank);
    Vcom_barr(vcom);

    /* Finish up */
    printf("PE %d: Exiting test program.\n",myrank);
    fflush(stdout);
    Vcom_dtor(&vcom);
    VASSERT( Vcom_finalize() );

    return 0;
}

