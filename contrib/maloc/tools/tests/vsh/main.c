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
 * rcsid="$Id: main.c,v 1.12 2008/03/12 05:14:00 fetk Exp $"
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * File:     main.c
 *
 * Purpose:  Main driver for testing the Vsh_shell.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */

#include <maloc/maloc.h>

#define VEMBED(rctag) VPRIVATE const char* rctag; \
    static void* use_rcsid=(0 ? &use_rcsid : (void*)&rcsid);
VEMBED(rcsid="$Id: main.c,v 1.12 2008/03/12 05:14:00 fetk Exp $")

/* Some help to organize the application-specific shell commands */
typedef enum APPcommand {
    app_none,
    app_help,
    app_stat,
    app_hello,
    app_bye
} APPcommand;

/*
 * ***************************************************************************
 * Routine:  APPgetCmd
 *
 * Purpose: The application-specific shell parser
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPRIVATE APPcommand APPgetCmd(void *thee, int argc, char *argv[])
{
    APPcommand theCmd = app_none;
    if (!strcmp(argv[0],"")) {
        theCmd = app_none;
    } else if (!strcmp(argv[0],"help")) {
        theCmd = app_help;
    } else if (!strcmp(argv[0],"stat")) {
        theCmd = app_stat;
    } else if (!strcmp(argv[0],"hello")) {
        theCmd = app_hello;
    } else if (!strcmp(argv[0],"bye")) {
        theCmd = app_bye;
    } else {
        theCmd = app_none;
    }
    return theCmd;
}

/*
 * ***************************************************************************
 * Routine:  APPsh
 *
 * Purpose: The application-specific shell for enriching the vsh_shell
 *          or for overriding the builtin vsh_shell commands.
 *
 * Return codes (required):
 *
 *     rc=0   --> APPsh does not know about this command
 *     rc=1   --> APPsh handled this command sucessfully
 *     rc=2   --> APPsh handled this command sucessfully and wants to exit!
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPRIVATE int APPsh(void *pthee, int argc, char *argv[])
{
    Vsh *thee = (Vsh*)pthee;

    int rc;
    APPcommand theCmd;
    static int init=0;

    /* one-time intialization */
    if (!init) {
       init=1;
    }

    /* default return code (success) */
    rc = 1;

    /* get the command */
    theCmd = APPgetCmd(thee, argc, argv);

    /* decode and execute the command */
    switch (theCmd) {

      case app_help:
        if (argc==1) {
            Vnm_print(1,"%s: Application-layer Help Menu:\n",
                Vsh_getenv(thee,"SHELL"));
            Vnm_print(1,"    help      --> Print this menu.\n");
            Vnm_print(1,"    stat      --> Print some environ variables.\n");
            Vnm_print(1,"    hello     --> Print 'Hello, World!'.\n");
            Vnm_print(1,"    bye       --> Print 'Bye, World!' and exit.\n");
            rc = 0;  /* pretend we didn't see it so subshell can help too */
        } else {
            rc = 0;  /* pretend we didn't see it so subshell can help too */
        }
        break;

      case app_stat:
        Vnm_print(1,"%s(APPsh): IFNAM=<%s>  LMAX=<%d>  LTOL=<%e>\n",
            Vsh_getenv(thee,"SHELL"),
            Vsh_getenv(thee,"IFNAM"),
            Vsh_getenvInt(thee,"LMAX"),
            Vsh_getenvReal(thee,"LTOL"));
        break;

      case app_hello:
        Vnm_print(1,"%s(APPsh): Hello, World!\n",
            Vsh_getenv(thee,"SHELL"));
        break;

      case app_bye:
        Vnm_print(1,"%s(APPsh): Bye, World!\n",
            Vsh_getenv(thee,"SHELL"));
        rc = 2;
        break;

      default:
        rc = 0;
        break;
    }
    return rc;
}

/*
 * ***************************************************************************
 * Routine:  main
 *
 * Purpose: The main driver for initiating the Vsh_shell.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
int main(int argc, char **argv)
{
    int rc;
    Vsh *vsh;

    /* ********************************************************************* */
    /* NOTE: avoiding Vnm_print before Vsh_ctor() allows output file tagging */
    /* ********************************************************************* */

    vsh = Vsh_ctor(VNULL, argc, argv);
    rc = Vsh_shell(vsh, VNULL, vsh, &APPsh);
    Vsh_dtor(&vsh);

    /* some final i/o */
    Vnm_print(1,"\n");
    Vnm_print(1,"maloc_leaks = [\n");
    Vnm_print(1,"%% --------------------------------------"
                "--------------------------------------\n");
    Vnm_print(1,"%%  Footprint        Areas       Malloc         Free"
                "    Highwater   Class\n"),
    Vnm_print(1,"%% --------------------------------------"
                "--------------------------------------\n");
    Vmem_print(VNULL);
    Vmem_printTotal();
    Vnm_print(1,"%% --------------------------------------"
                "--------------------------------------\n");
    Vnm_print(1,"];\n");

    /* normal return */
    return rc;
}

