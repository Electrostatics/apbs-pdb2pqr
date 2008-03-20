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
 * rcsid="$Id: vnmfb.c,v 1.8 2008/03/12 05:13:59 fetk Exp $"
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * File:     vnmfb.c
 *
 * Purpose:  FORTRAN bindings for the Vnm class methods.
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

#include "vnm_p.h"

VEMBED(rcsid="$Id: vnmfb.c,v 1.8 2008/03/12 05:13:59 fetk Exp $")

/*
 * ***************************************************************************
 * Class Vnm FORTRAN binding prototypes (default: no underscore, lowercase)
 * ***************************************************************************
 */

VEXTERNC int vrnd(void);
VEXTERNC int vrndmx(void);
VEXTERNC double vepsma(void);
VEXTERNC void vtstrt(int *timer, char *name, int *len);
VEXTERNC void vtstop(int *timer, char *name, int *len);
VEXTERNC int vsystm(char *cmd, int *len);
VEXTERNC void vnmprt(int *unit, char *strng, int *len);
VEXTERNC void vnmpri(int *unit, char *strng, int *len, int *val);
VEXTERNC void vnmprd(int *unit, char *strng, int *len, double *val);

/*
 * ***************************************************************************
 * Class Vnm FORTRAN bindings
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * Routine:  vrnd
 *
 * Purpose:  Returns a pseudo-random integer between 0 and VRANDMAX.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int vrnd(void)
{
    return VRAND;
}

/*
 * ***************************************************************************
 * Routine:  vrndmx
 *
 * Purpose:  Returns VRANDMAX.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int vrndmx(void)
{
    return VRANDMAX; 
}

/*
 * ***************************************************************************
 * Routine:  vepsma
 *
 * Purpose:  Computes the unit roundoff of the machine in single precision.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC double vepsma(void)
{
    return Vnm_epsmac();
}

/*
 * ***************************************************************************
 * Routine:  vtstrt
 *
 * Purpose:  Starts the timer on the particular machine.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void vtstrt(int *timer, char *name, int *len)
{
    int i;
    char buf[VMAX_ARGLEN];

    VASSERT( VMAX_ARGLEN > *len );
    for (i=0; i<*len; i++) {
        buf[i] = name[i];
    }
    buf[*len] = '\0';

    Vnm_tstart(*timer, buf);
}

/*
 * ***************************************************************************
 * Routine:  vtstop
 *
 * Purpose:  Stops the timer on the particular machine.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void vtstop(int *timer, char *name, int *len)
{
    int i;
    char buf[VMAX_ARGLEN];

    VASSERT( VMAX_ARGLEN > *len );
    for (i=0; i<*len; i++) {
        buf[i] = name[i];
    }
    buf[*len] = '\0';

    Vnm_tstop(*timer, buf);
}

/*
 * ***************************************************************************
 * Routine:  vsystm
 *
 * Purpose:  An improved ANSI-C "system" call.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int vsystm(char *cmd, int *len)
{
    int i;
    char buf[VMAX_ARGLEN];

    VASSERT( VMAX_ARGLEN > *len );
    for (i=0; i<*len; i++) {
        buf[i] = cmd[i];
    }
    buf[*len] = '\0';

    return Vnm_system(buf);
}   

/*
 * ***************************************************************************
 * Routine:  vnmprt
 *
 * Purpose:  External interface to the console i/o routine.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void vnmprt(int *unit, char *strng, int *len)
{
    int i;
    char buf[VMAX_ARGLEN];

    VASSERT( VMAX_ARGLEN > *len );
    for (i=0; i<*len; i++) {
        buf[i] = strng[i];
    }
    buf[*len] = '\0';

    Vnm_print(*unit, "%s\n", buf);
}

/*
 * ***************************************************************************
 * Routine:  vnmpri
 *
 * Purpose:  External interface to the console i/o routine.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void vnmpri(int *unit, char *strng, int *len, int *val)
{
    int i;
    char buf[VMAX_ARGLEN];

    VASSERT( VMAX_ARGLEN > *len );
    for (i=0; i<*len; i++) {
        buf[i] = strng[i];
    }
    buf[*len] = '\0';

    Vnm_print(*unit, "%s %d\n", buf, *val);
}

/*
 * ***************************************************************************
 * Routine:  vnmprd
 *
 * Purpose:  External interface to the console i/o routine.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void vnmprd(int *unit, char *strng, int *len, double *val)
{
    int i;
    char buf[VMAX_ARGLEN];

    VASSERT( VMAX_ARGLEN > *len );
    for (i=0; i<*len; i++) {
        buf[i] = strng[i];
    }
    buf[*len] = '\0';

    Vnm_print(*unit, "%s %e\n", buf, *val);
}

/*
 * ***************************************************************************
 * Class Vnm FORTRAN binding STUBS (no underscore, uppercase)
 * ***************************************************************************
 */

VPUBLIC int VRND(void)
{
    return vrnd();
}

VPUBLIC int VRNDMX(void)
{
    return vrndmx();
}

VPUBLIC double VEPSMA(void)
{
    return vepsma();
}

VPUBLIC void VTSTRT(int *timer, char *name, int *len)
{
    vtstrt(timer, name, len);
}

VPUBLIC void VTSTOP(int *timer, char *name, int *len)
{
    vtstop(timer, name, len);
}

VPUBLIC int VSYSTM(char *cmd, int *len)
{
    return vsystm(cmd, len);
}

VPUBLIC void VNMPRT(int *unit, char *strng, int *len)
{
    vnmprt(unit, strng, len);
}

VPUBLIC void VNMPRI(int *unit, char *strng, int *len, int *val)
{
    vnmpri(unit, strng, len, val);
}

VPUBLIC void VNMPRD(int *unit, char *strng, int *len, double *val)
{
    vnmprd(unit, strng, len, val);
}

/*
 * ***************************************************************************
 * Class Vnm FORTRAN binding STUBS (one underscore, lowercase)
 * ***************************************************************************
 */

VPUBLIC int vrnd_(void)
{
    return vrnd();
}

VPUBLIC int vrndmx_(void)
{
    return vrndmx();
}

VPUBLIC double vepsma_(void)
{
    return vepsma();
}

VPUBLIC void vtstrt_(int *timer, char *name, int *len)
{
    vtstrt(timer, name, len);
}

VPUBLIC void vtstop_(int *timer, char *name, int *len)
{
    vtstop(timer, name, len);
}

VPUBLIC int vsystm_(char *cmd, int *len)
{
    return vsystm(cmd, len);
}

VPUBLIC void vnmprt_(int *unit, char *strng, int *len)
{
    vnmprt(unit, strng, len);
}

VPUBLIC void vnmpri_(int *unit, char *strng, int *len, int *val)
{
    vnmpri(unit, strng, len, val);
}

VPUBLIC void vnmprd_(int *unit, char *strng, int *len, double *val)
{
    vnmprd(unit, strng, len, val);
}

/*
 * ***************************************************************************
 * Class Vnm FORTRAN binding STUBS (one underscore, uppercase)
 * ***************************************************************************
 */

VPUBLIC int VRND_(void)
{
    return vrnd();
}

VPUBLIC int VRNDMX_(void)
{
    return vrndmx();
}

VPUBLIC double VEPSMA_(void)
{
    return vepsma();
}

VPUBLIC void VTSTRT_(int *timer, char *name, int *len)
{
    vtstrt(timer, name, len);
}

VPUBLIC void VTSTOP_(int *timer, char *name, int *len)
{
    vtstop(timer, name, len);
}

VPUBLIC int VSYSTM_(char *cmd, int *len)
{
    return vsystm(cmd, len);
}

VPUBLIC void VNMPRT_(int *unit, char *strng, int *len)
{
    vnmprt(unit, strng, len);
}

VPUBLIC void VNMPRI_(int *unit, char *strng, int *len, int *val)
{
    vnmpri(unit, strng, len, val);
}

VPUBLIC void VNMPRD_(int *unit, char *strng, int *len, double *val)
{
    vnmprd(unit, strng, len, val);
}

/*
 * ***************************************************************************
 * Class Vnm FORTRAN binding STUBS (double underscore, lowercase)
 * ***************************************************************************
 */

VPUBLIC int vrnd__(void)
{
    return vrnd();
}

VPUBLIC int vrndmx__(void)
{
    return vrndmx();
}

VPUBLIC double vepsma__(void)
{
    return vepsma();
}

VPUBLIC void vtstrt__(int *timer, char *name, int *len)
{
    vtstrt(timer, name, len);
}

VPUBLIC void vtstop__(int *timer, char *name, int *len)
{
    vtstop(timer, name, len);
}

VPUBLIC int vsystm__(char *cmd, int *len)
{
    return vsystm(cmd, len);
}

VPUBLIC void vnmprt__(int *unit, char *strng, int *len)
{
    vnmprt(unit, strng, len);
}

VPUBLIC void vnmpri__(int *unit, char *strng, int *len, int *val)
{
    vnmpri(unit, strng, len, val);
}

VPUBLIC void vnmprd__(int *unit, char *strng, int *len, double *val)
{
    vnmprd(unit, strng, len, val);
}

/*
 * ***************************************************************************
 * Class Vnm FORTRAN binding STUBS (double underscore, uppercase)
 * ***************************************************************************
 */

VPUBLIC int VRND__(void)
{
    return vrnd();
}

VPUBLIC int VRNDMX__(void)
{
    return vrndmx();
}

VPUBLIC double VEPSMA__(void)
{
    return vepsma();
}

VPUBLIC void VTSTRT__(int *timer, char *name, int *len)
{
    vtstrt(timer, name, len);
}

VPUBLIC void VTSTOP__(int *timer, char *name, int *len)
{
    vtstop(timer, name, len);
}

VPUBLIC int VSYSTM__(char *cmd, int *len)
{
    return vsystm(cmd, len);
}

VPUBLIC void VNMPRT__(int *unit, char *strng, int *len)
{
    vnmprt(unit, strng, len);
}

VPUBLIC void VNMPRI__(int *unit, char *strng, int *len, int *val)
{
    vnmpri(unit, strng, len, val);
}

VPUBLIC void VNMPRD__(int *unit, char *strng, int *len, double *val)
{
    vnmprd(unit, strng, len, val);
}

