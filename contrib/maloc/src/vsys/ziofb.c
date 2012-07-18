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
 * rcsid="$Id: ziofb.c,v 1.16 2008/03/12 05:13:59 fetk Exp $"
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * File:     ziofb.c
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

VEMBED(rcsid="$Id: ziofb.c,v 1.16 2008/03/12 05:13:59 fetk Exp $")

/*
 * ***************************************************************************
 * Class Vio FORTRAN binding prototypes (default: no underscore, lowercase)
 * ***************************************************************************
 */

VEXTERNC void ziosta(void);
VEXTERNC void ziostp(void);
VEXTERNC void zioctr(Vio *sock,
    char type[4], char frmt[3], 
    char *host, int *lenh, char *file, int *lenf,
    char mode[1], int *iflag);
VEXTERNC void ziodtr(Vio *sock, int *iflag);
VEXTERNC void zioutl(Vio *sock, char mode[1], int *iflag);
VEXTERNC void zioint(Vio *sock, int *ival, int *len);
VEXTERNC void zioflt(Vio *sock, float *fval, int *len);
VEXTERNC void ziodbl(Vio *sock, double *dval, int *len);
VEXTERNC void ziostr(Vio *sock, char *sval, int *len);

/*
 * ***************************************************************************
 * Class Vio FORTRAN bindings
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * Routine:  ziosta
 *
 * Purpose:  Start the Vio communication layer.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void ziosta(void)
{
    Vio_start();
}

/*
 * ***************************************************************************
 * Routine:  ziostp
 *
 * Purpose:  Stop the Vio communication layer.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void ziostp(void)
{
    Vio_stop();
}

/*
 * ***************************************************************************
 * Routine:  zioctr
 *
 * Purpose:  Construct the Vio object.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void zioctr(Vio *sock,
    char type[4], char frmt[3],
    char *host, int *lenh, char *file, int *lenf,
    char mode[1], int *iflag)
{
    int i;
    char phost[VMAX_ARGLEN], pfile[VMAX_ARGLEN], ptype[VMAX_ARGLEN];
    char pfrmt[VMAX_ARGLEN], pmode[VMAX_ARGLEN];

#if 0
    Vio sockSize;
    fprintf(stderr,"zioctr: Vio structure size is exactly <%d> bytes.\n",
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

    VJMPERR1(0 != Vio_ctor2(sock, ptype, pfrmt, phost, pfile, pmode));

    *iflag = 0;
    return;

  VERROR1:
    *iflag = 1;
    return;
}

/*
 * ***************************************************************************
 * Routine:  ziodtr
 *
 * Purpose:  Destruct the Vio object.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void ziodtr(Vio *sock, int *iflag)
{
    Vio_dtor2(sock);
    *iflag = 0;
    return;
}

/*
 * ***************************************************************************
 * Routine:  zioutl
 *
 * Purpose:  Vio state utility.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void zioutl(Vio *sock, char mode[1], int *iflag)
{
    char pmode[VMAX_ARGLEN];

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

        *iflag = 0;
        return;

    } else if (!strcmp(pmode,"c")) {

        if ( sock->rwkey == VIO_R ) {
            Vio_acceptFree(sock);
        } else if ( sock->rwkey == VIO_W ) {
            Vio_connectFree(sock);
        } else { VJMPERR1(0); }

        *iflag = 0;
        return;

    } else { VJMPERR1(0); }

  VERROR1:
    *iflag = 1;
    return;
}

/*
 * ***************************************************************************
 * Routine:  zioint
 *
 * Purpose:  Integer READ/WRITE.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void zioint(Vio *sock, int *ival, int *len)
{
    int i;

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
 * Routine:  zioflt
 *
 * Purpose:  Float READ/WRITE.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void zioflt(Vio *sock, float *fval, int *len)
{
    int i;

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
 * Routine:  ziodbl
 *
 * Purpose:  Double READ/WRITE.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void ziodbl(Vio *sock, double *dval, int *len)
{
    int i;

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
 * Routine:  ziostr
 *
 * Purpose:  String READ/WRITE.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void ziostr(Vio *sock, char *sval, int *len)
{
    int i;
    char buf[VMAX_BUFSIZE];

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

VPUBLIC void ZIOSTA(void)
{
    ziosta();
}

VPUBLIC void ZIOSTP(void)
{
    ziostp();
}

VPUBLIC void ZIOCTR(Vio *sock,
    char type[4], char frmt[3], 
    char *host, int *lenh, char *file, int *lenf,
    char mode[1], int *iflag)
{
    zioctr(sock, type, frmt, host, lenh, file, lenf, mode, iflag);
}

VPUBLIC void ZIODTR(Vio *sock, int *iflag)
{
    ziodtr(sock, iflag);
}

VPUBLIC void ZIOUTL(Vio *sock, char mode[1], int *iflag)
{
    zioutl(sock, mode, iflag);
}

VPUBLIC void ZIOINT(Vio *sock, int *ival, int *len)
{
    zioint(sock, ival, len);
}

VPUBLIC void ZIOFLT(Vio *sock, float *fval, int *len)
{
    zioflt(sock, fval, len);
}

VPUBLIC void ZIODBL(Vio *sock, double *dval, int *len)
{
    ziodbl(sock, dval, len);
}

VPUBLIC void ZIOSTR(Vio *sock, char *sval, int *len)
{
    ziostr(sock, sval, len);
}

/*
 * ***************************************************************************
 * Class Vio FORTRAN binding STUBS (single underscore, lowercase)
 * ***************************************************************************
 */

VPUBLIC void ziosta_(void)
{
    ziosta();
}

VPUBLIC void ziostp_(void)
{
    ziostp();
}

VPUBLIC void zioctr_(Vio *sock,
    char type[4], char frmt[3], 
    char *host, int *lenh, char *file, int *lenf,
    char mode[1], int *iflag)
{
    zioctr(sock, type, frmt, host, lenh, file, lenf, mode, iflag);
}

VPUBLIC void ziodtr_(Vio *sock, int *iflag)
{
    ziodtr(sock, iflag);
}

VPUBLIC void zioutl_(Vio *sock, char mode[1], int *iflag)
{
    zioutl(sock, mode, iflag);
}

VPUBLIC void zioint_(Vio *sock, int *ival, int *len)
{
    zioint(sock, ival, len);
}

VPUBLIC void zioflt_(Vio *sock, float *fval, int *len)
{
    zioflt(sock, fval, len);
}

VPUBLIC void ziodbl_(Vio *sock, double *dval, int *len)
{
    ziodbl(sock, dval, len);
}

VPUBLIC void ziostr_(Vio *sock, char *sval, int *len)
{
    ziostr(sock, sval, len);
}

/*
 * ***************************************************************************
 * Class Vio FORTRAN binding STUBS (single underscore, uppercase)
 * ***************************************************************************
 */

VPUBLIC void ZIOSTA_(void)
{
    ziosta();
}

VPUBLIC void ZIOSTP_(void)
{
    ziostp();
}

VPUBLIC void ZIOCTR_(Vio *sock,
    char type[4], char frmt[3], 
    char *host, int *lenh, char *file, int *lenf,
    char mode[1], int *iflag)
{
    zioctr(sock, type, frmt, host, lenh, file, lenf, mode, iflag);
}

VPUBLIC void ZIODTR_(Vio *sock, int *iflag)
{
    ziodtr(sock, iflag);
}

VPUBLIC void ZIOUTL_(Vio *sock, char mode[1], int *iflag)
{
    zioutl(sock, mode, iflag);
}

VPUBLIC void ZIOINT_(Vio *sock, int *ival, int *len)
{
    zioint(sock, ival, len);
}

VPUBLIC void ZIOFLT_(Vio *sock, float *fval, int *len)
{
    zioflt(sock, fval, len);
}

VPUBLIC void ZIODBL_(Vio *sock, double *dval, int *len)
{
    ziodbl(sock, dval, len);
}

VPUBLIC void ZIOSTR_(Vio *sock, char *sval, int *len)
{
    ziostr(sock, sval, len);
}

/*
 * ***************************************************************************
 * Class Vio FORTRAN binding STUBS (double underscore, lowercase)
 * ***************************************************************************
 */

VPUBLIC void ziosta__(void)
{
    ziosta();
}

VPUBLIC void ziostp__(void)
{
    ziostp();
}

VPUBLIC void zioctr__(Vio *sock,
    char type[4], char frmt[3], 
    char *host, int *lenh, char *file, int *lenf,
    char mode[1], int *iflag)
{
    zioctr(sock, type, frmt, host, lenh, file, lenf, mode, iflag);
}

VPUBLIC void ziodtr__(Vio *sock, int *iflag)
{
    ziodtr(sock, iflag);
}

VPUBLIC void zioutl__(Vio *sock, char mode[1], int *iflag)
{
    zioutl(sock, mode, iflag);
}

VPUBLIC void zioint__(Vio *sock, int *ival, int *len)
{
    zioint(sock, ival, len);
}

VPUBLIC void zioflt__(Vio *sock, float *fval, int *len)
{
    zioflt(sock, fval, len);
}

VPUBLIC void ziodbl__(Vio *sock, double *dval, int *len)
{
    ziodbl(sock, dval, len);
}

VPUBLIC void ziostr__(Vio *sock, char *sval, int *len)
{
    ziostr(sock, sval, len);
}

/*
 * ***************************************************************************
 * Class Vio FORTRAN binding STUBS (double underscore, uppercase)
 * ***************************************************************************
 */

VPUBLIC void ZIOSTA__(void)
{
    ziosta();
}

VPUBLIC void ZIOSTP__(void)
{
    ziostp();
}

VPUBLIC void ZIOCTR__(Vio *sock,
    char type[4], char frmt[3], 
    char *host, int *lenh, char *file, int *lenf,
    char mode[1], int *iflag)
{
    zioctr(sock, type, frmt, host, lenh, file, lenf, mode, iflag);
}

VPUBLIC void ZIODTR__(Vio *sock, int *iflag)
{
    ziodtr(sock, iflag);
}

VPUBLIC void ZIOUTL__(Vio *sock, char mode[1], int *iflag)
{
    zioutl(sock, mode, iflag);
}

VPUBLIC void ZIOINT__(Vio *sock, int *ival, int *len)
{
    zioint(sock, ival, len);
}

VPUBLIC void ZIOFLT__(Vio *sock, float *fval, int *len)
{
    zioflt(sock, fval, len);
}

VPUBLIC void ZIODBL__(Vio *sock, double *dval, int *len)
{
    ziodbl(sock, dval, len);
}

VPUBLIC void ZIOSTR__(Vio *sock, char *sval, int *len)
{
    ziostr(sock, sval, len);
}

