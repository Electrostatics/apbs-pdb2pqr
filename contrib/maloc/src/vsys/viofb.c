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
 * rcsid="$Id: viofb.c,v 1.18 2008/03/12 05:13:59 fetk Exp $"
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * File:     viofb.c
 *
 * Purpose:  FORTRAN bindings for the Vio class methods.
 *
 * Notes:    We provide FORTRAN stubs for the following manglings:
 *
 *               vrnd   --> no underscore,     lowercase (default)
 *               VRND   --> no underscore,     uppercase
 *               vrnd_  --> single underscore, lowercase
 *               VRND_  --> single underscore, uppercase
 *               vrnd__ --> double underscore, lowercase
 *               VRND__ --> double underscore, uppercase
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */

#include "vio_p.h"

VEMBED(rcsid="$Id: viofb.c,v 1.18 2008/03/12 05:13:59 fetk Exp $")

#define MAXVIO 10
VPRIVATE Vio theVio[MAXVIO];
VPRIVATE int stack[MAXVIO];
VPRIVATE int stackPtr = 0;

/*
 * ***************************************************************************
 * Class Vio FORTRAN binding prototypes (default: no underscore, lowercase)
 * ***************************************************************************
 */

VEXTERNC void viosta(void);
VEXTERNC void viostp(void);

VEXTERNC int vioctr(char type[4], char frmt[3], 
    char *host, int *lenh, char *file, int *lenf,
    char mode[1]);
VEXTERNC int viodtr(int *socknum);
VEXTERNC int vioutl(int *socknum, char mode[1]);

VEXTERNC void vioint(int *socknum, int *ival, int *len);
VEXTERNC void vioflt(int *socknum, float *fval, int *len);
VEXTERNC void viodbl(int *socknum, double *dval, int *len);
VEXTERNC void viostr(int *socknum, char *sval, int *len);

/*
 * ***************************************************************************
 * Class Vio FORTRAN bindings
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * Routine:  viosta
 *
 * Purpose:  Start the Vio communication layer.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void viosta(void)
{
    int i;

    for (i=0; i<MAXVIO; i++) {
        stack[i] = i+1;
    }
    stack[MAXVIO-1] = -1;
    stackPtr = 0;

    Vio_start();
}

/*
 * ***************************************************************************
 * Routine:  viostp
 *
 * Purpose:  Stop the Vio communication layer.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void viostp(void)
{
    Vio_stop();
}

/*
 * ***************************************************************************
 * Routine:  vioctr
 *
 * Purpose:  Construct the Vio object.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int vioctr(char type[4], char frmt[3],
    char *host, int *lenh, char *file, int *lenf,
    char mode[1])
{
    int i, socknum;
    char phost[VMAX_ARGLEN], pfile[VMAX_ARGLEN], ptype[VMAX_ARGLEN];
    char pfrmt[VMAX_ARGLEN], pmode[VMAX_ARGLEN];
    Vio *sock;

#if 0
    Vio sockSize;
    fprintf(stderr,"vioctr: Vio structure size is exactly <%d> bytes.\n",
        sizeof(sockSize) );
#endif

    for (i=0; i<4; i++) ptype[i] = type[i];
    ptype[4] = '\0';
    for (i=0; i<3; i++) pfrmt[i] = frmt[i];
    pfrmt[3] = '\0';
    for (i=0; i<*lenh; i++) phost[i] = host[i];
    phost[*lenh] = '\0';
    for (i=0; i<*lenf; i++) pfile[i] = file[i];
    pfile[*lenf] = '\0';
    pmode[0] = mode[0];
    pmode[1] = '\0';

    VASSERT( (0 <= stackPtr) && (stackPtr < MAXVIO) );

    socknum = stackPtr;
    stackPtr = stack[socknum];
    sock = &theVio[socknum];
    VJMPERR1(0 != Vio_ctor2(sock, ptype, pfrmt, phost, pfile, pmode));

    return socknum;

  VERROR1:
    return -1;
}

/*
 * ***************************************************************************
 * Routine:  viodtr
 *
 * Purpose:  Destruct the Vio object.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int viodtr(int *socknum)
{
    Vio *sock = &theVio[*socknum];

    VASSERT( (0 <= *socknum) && (*socknum < MAXVIO) );

    Vio_dtor2(sock);

    stack[*socknum] = stackPtr;
    stackPtr = *socknum;

    return 0;
}

/*
 * ***************************************************************************
 * Routine:  vioutl
 *
 * Purpose:  Vio state utility.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int vioutl(int *socknum, char mode[1])
{
    Vio *sock = &theVio[*socknum];
    char pmode[VMAX_ARGLEN];

    VASSERT( (0 <= *socknum) && (*socknum < MAXVIO) );

    pmode[0] = mode[0];
    pmode[1] = '\0';

    if ( !strcmp(pmode,"o") ) {

        if ( sock->rwkey == VIO_R ) {
            /* BLOCKING READ (blocking accept) */
            VJMPERR1( 0 <= Vio_accept(sock,0) );
        } else if ( sock->rwkey == VIO_W ) {
            /* BLOCKING WRITE (blocking connect) */
            VJMPERR1( 0 <= Vio_connect(sock,0) );
        } else { VJMPERR1(0); }

        return 0;

    } else if (!strcmp(pmode,"c")) {

        if ( sock->rwkey == VIO_R ) {
            Vio_acceptFree(sock);
        } else if ( sock->rwkey == VIO_W ) {
            Vio_connectFree(sock);
        } else { VJMPERR1(0); }

        return 0;

    } else { VJMPERR1(0); }

  VERROR1:
    return 1;
}

/*
 * ***************************************************************************
 * Routine:  vioint
 *
 * Purpose:  Integer READ/WRITE.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void vioint(int *socknum, int *ival, int *len)
{
    Vio *sock = &theVio[*socknum];
    int i;

    VASSERT( (0 <= *socknum) && (*socknum < MAXVIO) );

    if ( sock->rwkey == VIO_R ) {
        for (i=0; i<*len; i++)
            Vio_scanf(sock,"%d",&(ival[i]));
    } else if ( sock->rwkey == VIO_W ) {
        for (i=0; i<*len; i++)
            Vio_printf(sock,"%d ",ival[i]);
        Vio_printf(sock,"\n");
    }
}

/*
 * ***************************************************************************
 * Routine:  vioflt
 *
 * Purpose:  Float READ/WRITE.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void vioflt(int *socknum, float *fval, int *len)
{
    Vio *sock = &theVio[*socknum];
    int i;

    VASSERT( (0 <= *socknum) && (*socknum < MAXVIO) );

    if ( sock->rwkey == VIO_R ) {
        for (i=0; i<*len; i++)
            Vio_scanf(sock,"%e",&(fval[i]));
    } else if ( sock->rwkey == VIO_W ) {
        for (i=0; i<*len; i++)
            Vio_printf(sock,"%e ",fval[i]);
        Vio_printf(sock,"\n");
    }
}

/*
 * ***************************************************************************
 * Routine:  viodbl
 *
 * Purpose:  Double READ/WRITE.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void viodbl(int *socknum, double *dval, int *len)
{
    Vio *sock = &theVio[*socknum];
    int i;

    VASSERT( (0 <= *socknum) && (*socknum < MAXVIO) );

    if ( sock->rwkey == VIO_R ) {
        for (i=0; i<*len; i++)
            Vio_scanf(sock,"%le",&(dval[i]));
    } else if ( sock->rwkey == VIO_W ) {
        for (i=0; i<*len; i++)
            Vio_printf(sock,"%le ",dval[i]);
        Vio_printf(sock,"\n");
    }
}

/*
 * ***************************************************************************
 * Routine:  viostr
 *
 * Purpose:  String READ/WRITE.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void viostr(int *socknum, char *sval, int *len)
{
    Vio *sock = &theVio[*socknum];
    int i;
    char buf[VMAX_BUFSIZE];

    VASSERT( (0 <= *socknum) && (*socknum < MAXVIO) );

    if ( sock->rwkey == VIO_R ) {
        Vio_scanf(sock,"%s",buf);
        VASSERT( (int)strlen(buf) == *len );
        for (i=0; i<*len; i++) sval[i] = buf[i];
    } else if ( sock->rwkey == VIO_W ) {
        for (i=0; i<*len; i++) buf[i] = sval[i];
        buf[*len] = '\0';
        Vio_printf(sock,"%s\n",buf);
    }
}

/*
 * ***************************************************************************
 * Class Vio FORTRAN binding STUBS (no underscore, uppercase)
 * ***************************************************************************
 */

VPUBLIC void VIOSTA(void)
{
    viosta();
}

VPUBLIC void VIOSTP(void)
{
    viostp();
}

VPUBLIC int VIOCTR(char type[4], char frmt[3], 
    char *host, int *lenh, char *file, int *lenf,
    char mode[1])
{
    return vioctr(type, frmt, host, lenh, file, lenf, mode);
}

VPUBLIC int VIODTR(int *socknum)
{
    return viodtr(socknum);
}

VPUBLIC int VIOUTL(int *socknum, char mode[1])
{
    return vioutl(socknum, mode);
}

VPUBLIC void VIOINT(int *socknum, int *ival, int *len)
{
    vioint(socknum, ival, len);
}

VPUBLIC void VIOFLT(int *socknum, float *fval, int *len)
{
    vioflt(socknum, fval, len);
}

VPUBLIC void VIODBL(int *socknum, double *dval, int *len)
{
    viodbl(socknum, dval, len);
}

VPUBLIC void VIOSTR(int *socknum, char *sval, int *len)
{
    viostr(socknum, sval, len);
}

/*
 * ***************************************************************************
 * Class Vio FORTRAN binding STUBS (single underscore, lowercase)
 * ***************************************************************************
 */

VPUBLIC void viosta_(void)
{
    viosta();
}

VPUBLIC void viostp_(void)
{
    viostp();
}

VPUBLIC int vioctr_(char type[4], char frmt[3], 
    char *host, int *lenh, char *file, int *lenf,
    char mode[1])
{
    return vioctr(type, frmt, host, lenh, file, lenf, mode);
}

VPUBLIC int viodtr_(int *socknum)
{
    return viodtr(socknum);
}

VPUBLIC int vioutl_(int *socknum, char mode[1])
{
    return vioutl(socknum, mode);
}

VPUBLIC void vioint_(int *socknum, int *ival, int *len)
{
    vioint(socknum, ival, len);
}

VPUBLIC void vioflt_(int *socknum, float *fval, int *len)
{
    vioflt(socknum, fval, len);
}

VPUBLIC void viodbl_(int *socknum, double *dval, int *len)
{
    viodbl(socknum, dval, len);
}

VPUBLIC void viostr_(int *socknum, char *sval, int *len)
{
    viostr(socknum, sval, len);
}

/*
 * ***************************************************************************
 * Class Vio FORTRAN binding STUBS (double underscore, lowercase)
 * ***************************************************************************
 */

VPUBLIC void viosta__(void)
{
    viosta();
}

VPUBLIC void viostp__(void)
{
    viostp();
}

VPUBLIC int vioctr__(char type[4], char frmt[3], 
    char *host, int *lenh, char *file, int *lenf,
    char mode[1])
{
    return vioctr(type, frmt, host, lenh, file, lenf, mode);
}

VPUBLIC int viodtr__(int *socknum)
{
    return viodtr(socknum);
}

VPUBLIC int vioutl__(int *socknum, char mode[1])
{
    return vioutl(socknum, mode);
}

VPUBLIC void vioint__(int *socknum, int *ival, int *len)
{
    vioint(socknum, ival, len);
}

VPUBLIC void vioflt__(int *socknum, float *fval, int *len)
{
    vioflt(socknum, fval, len);
}

VPUBLIC void viodbl__(int *socknum, double *dval, int *len)
{
    viodbl(socknum, dval, len);
}

VPUBLIC void viostr__(int *socknum, char *sval, int *len)
{
    viostr(socknum, sval, len);
}

/*
 * ***************************************************************************
 * Class Vio FORTRAN binding STUBS (single underscore, uppercase)
 * ***************************************************************************
 */

VPUBLIC void VIOSTA_(void)
{
    viosta();
}

VPUBLIC void VIOSTP_(void)
{
    viostp();
}

VPUBLIC int VIOCTR_(char type[4], char frmt[3], 
    char *host, int *lenh, char *file, int *lenf,
    char mode[1])
{
    return vioctr(type, frmt, host, lenh, file, lenf, mode);
}

VPUBLIC int VIODTR_(int *socknum)
{
    return viodtr(socknum);
}

VPUBLIC int VIOUTL_(int *socknum, char mode[1])
{
    return vioutl(socknum, mode);
}

VPUBLIC void VIOINT_(int *socknum, int *ival, int *len)
{
    vioint(socknum, ival, len);
}

VPUBLIC void VIOFLT_(int *socknum, float *fval, int *len)
{
    vioflt(socknum, fval, len);
}

VPUBLIC void VIODBL_(int *socknum, double *dval, int *len)
{
    viodbl(socknum, dval, len);
}

VPUBLIC void VIOSTR_(int *socknum, char *sval, int *len)
{
    viostr(socknum, sval, len);
}

/*
 * ***************************************************************************
 * Class Vio FORTRAN binding STUBS (double underscore, uppercase)
 * ***************************************************************************
 */

VPUBLIC void VIOSTA__(void)
{
    viosta();
}

VPUBLIC void VIOSTP__(void)
{
    viostp();
}

VPUBLIC int VIOCTR__(char type[4], char frmt[3], 
    char *host, int *lenh, char *file, int *lenf,
    char mode[1])
{
    return vioctr(type, frmt, host, lenh, file, lenf, mode);
}

VPUBLIC int VIODTR__(int *socknum)
{
    return viodtr(socknum);
}

VPUBLIC int VIOUTL__(int *socknum, char mode[1])
{
    return vioutl(socknum, mode);
}

VPUBLIC void VIOINT__(int *socknum, int *ival, int *len)
{
    vioint(socknum, ival, len);
}

VPUBLIC void VIOFLT__(int *socknum, float *fval, int *len)
{
    vioflt(socknum, fval, len);
}

VPUBLIC void VIODBL__(int *socknum, double *dval, int *len)
{
    viodbl(socknum, dval, len);
}

VPUBLIC void VIOSTR__(int *socknum, char *sval, int *len)
{
    viostr(socknum, sval, len);
}

