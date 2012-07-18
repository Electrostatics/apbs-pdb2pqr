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
 * rcsid="$Id: psh.c,v 1.40 2008/03/12 05:13:58 fetk Exp $"
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * File:     psh.c
 *
 * Purpose:  Vsh_pshell w/ infrastructure (parallel extension of Vsh_shell).
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */

#include "psh_p.h"

VEMBED(rcsid="$Id: psh.c,v 1.40 2008/03/12 05:13:58 fetk Exp $")

/* parallel shell commands */
typedef enum PSH_command {
    pshcom_none,
    pshcom_ignore,
    pshcom_set,
    pshcom_help,
    pshcom_vmp_snd,
    pshcom_vmp_rcv,
    pshcom_vmp_bar
} PSH_command;

/* local variables we need (THIS CODE IS NOT REENTRANT!) */
VPRIVATE Vmp *theeVMP = VNULL;
VPRIVATE Vsh *theePSH = VNULL;
VPRIVATE int (*theeFunc)(void *thee, int argc, char **argv) = VNULL;

/*
 * ***************************************************************************
 * Routine:  Vsh_publishVars
 *
 * Purpose:  Publish environment variables.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPRIVATE void PSH_publishVars(Vsh *thee)
{
    int i, numVars = 4;
    typedef struct vshVars {
        char envi[VMAX_ARGLEN];
        char valu[VMAX_ARGLEN];
        char info[VMAX_ARGLEN];
    } vshVars;
    vshVars envVars[] = {
        /* --------   -----       ----------- */
        /* VARIABLE   VALUE       EXPLANATION */
        /* --------   -----       ----------- */
        /* ===[ VMP=4 ]=== */
        { "VMP_I",    "0",
            "VMP id (my VMP process number)" },
        { "VMP_N",    "1",
            "VMP nproc (number of VMP processes in this execution)" },
        { "VMP_P",    "0",
            "VMP send/recv partner (current VMP send/recv partner)" },
        { "VMP_F",    "-1",
            "VMP cmd effect (-1=all,0=proc0,1=proc1,...,N=procN)" }
    };  

    /* publish remaining variables */
    for (i=0; i<numVars; i++) {
        VASSERT( Vsh_putenv(     thee, envVars[i].envi, envVars[i].valu )
              && Vsh_putenvInfo( thee, envVars[i].envi, envVars[i].info ) );
    }
}

/*
 * ***************************************************************************
 * Routine:  PSH_getCmd
 *
 * Purpose:  Decode the input string into a legal command.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPRIVATE PSH_command PSH_getCmd(int argc, char **argv)
{
    PSH_command theCmd = pshcom_none;
    if (!strcmp(argv[0],"")) {
        theCmd = pshcom_none;
    } else if (!strcmp(argv[0],"set")) {
        theCmd = pshcom_set;
    } else if (!strcmp(argv[0],"help")) {
        theCmd = pshcom_help;
    } else if (!strcmp(argv[0],"vmp_snd")) {
        theCmd = pshcom_vmp_snd;
    } else if (!strcmp(argv[0],"vmp_rcv")) {
        theCmd = pshcom_vmp_rcv;
    } else if (!strcmp(argv[0],"vmp_bar")) {
        theCmd = pshcom_vmp_bar;
    } else {
        theCmd = pshcom_none;
    }
    return theCmd;
}

/*
 * ***************************************************************************
 * Routine:  PSH_builtin
 *
 * Purpose: Parallel extensions to the vsh.
 *
 * Return codes (required):
 *
 *     rc=0   --> Psh does not know about this command
 *     rc=1   --> Psh handled this command sucessfully
 *     rc=2   --> Psh handled this command sucessfully and wants to exit!
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPRIVATE int PSH_builtin(void *pthee, int argc, char **argv)
{
    int rc, me, src, des, bufLen;
    unsigned int bufLenMessage;
    char *bufPtr;
    PSH_command theCmd;

    static int init=0;
    static char vmp[VMAX_BUFSIZE], vmp_min[VMAX_BUFSIZE];
    const char *stmp;

    /* one-time intialization */
    if (!init) {
        init=1;

        /* make the minimal vmp message (%s slots = 2) */
        stmp = "%s: pVsh-layer Help Menu: \n"
            "    help vmp  --> Help on %s communication commands\n";
        sprintf(vmp_min,stmp,theePSH->PR,theePSH->PR);

        /* make the full vmp message (%s slots = 1) */
        stmp = "%s: Parallel shell extensions: \n"
            "    vmp_snd   --> VMP send local buffer to selected proc\n"
            "    vmp_rcv   --> VMP recv into local buffer\n"
            "    vmp_bar   --> VMP synchronization barrier\n";
        sprintf(vmp,stmp,theePSH->PR);

        /* publish VMP variables we need */
        PSH_publishVars(theePSH);

        /* initialize VMP variables */
        if (theeVMP != VNULL) {
            Vsh_putenvInt(theePSH,"VMP_I",Vmp_rank(theeVMP));
            Vsh_putenvInt(theePSH,"VMP_N",Vmp_size(theeVMP));
        } else {
            Vsh_putenvInt(theePSH,"VMP_I",0);
            Vsh_putenvInt(theePSH,"VMP_N",1);
        }
    }

    /* get the command */
    theCmd = PSH_getCmd(argc, argv);

    /*
     * see if we are supposed to execute the command
     * parallel mode: we execute the command in only three situations:
     *     (1) everyone has the focus (VMP_F==-1)
     *     (2) we alone have the focus (VMP_F==VMP_I)
     *     (3) we are setting the focus (set VMP_F X)
     */
    if ( Vsh_getenvInt(theePSH,"VMP_F") == -1 ) {
        /*
         * everyone is supposed to execute it
         */
        rc = 0;
    } else if ( Vsh_getenvInt(theePSH,"VMP_F") 
             == Vsh_getenvInt(theePSH,"VMP_I") ) {
        /*
         * only we are supposed to execute it
         */
        rc = 0;
    } else if ( (theCmd == pshcom_set)
             && (!strcmp(argv[1],"VMP_F"))
             && (argc == 3) ) {
        /*
         * OK, it is a legal "set VMP_F X" command.
         * SO, we ARE supposed to execute it.
         * NOTE: we won't recognize "set" as one of our
         *       legal commands, so it will fall through
         *       to the base VSH builtin which will deal with it.
         */
        rc = 0;
    } else {
        /*
         * we are supposed to ignore it
         */
        theCmd = pshcom_ignore;
        rc = 1;
    }

    /* have a look at the command if we are not ignoring it */
    if (theCmd != pshcom_ignore) {

        /* the normal vsh shell gets first shot at the command */
        if (theeFunc != VNULL) {
            rc = (*(theeFunc))(pthee,argc,argv);
            if (rc != 0) return rc;
        }

        /* decode and execute the command */
        switch (theCmd) {

          case pshcom_help:
            if (argc==1) {
                Vnm_print(1,"%s",vmp_min);
                rc = 0;  /* pretend we didn't see it so subshell can help too */
            } else if ((argc==2) && (!strcmp(argv[1],"vmp"))) {
                Vnm_print(1,"%s",vmp);
                rc = 1;
            } else {
                rc = 0;  /* pretend we didn't see it so subshell can help too */
            }
            break;

          case pshcom_vmp_snd:
            me     = Vsh_getenvInt(theePSH,"VMP_I");
            des    = Vsh_getenvInt(theePSH,"VMP_P");
            bufLen = theePSH->bufsize;
            bufPtr = theePSH->buf;

            Vnm_print(2,"Vsh_builtIn: [%d --> %d] sending mesg size=<%d>\n",
                me, des, bufLen);
            bufLenMessage = (unsigned int)bufLen;
            Vmp_send(theeVMP, des, (char*)&bufLenMessage, 4);

            Vnm_print(2,"Vsh_builtIn: [%d --> %d] sending the real mesg.\n",
                me, des);
            Vmp_send(theeVMP, des, bufPtr, bufLen);
            rc = 1;
            break;

          case pshcom_vmp_rcv:
            me  = Vsh_getenvInt(theePSH,"VMP_I");
            src = Vsh_getenvInt(theePSH,"VMP_P");

            Vmp_recv(theeVMP, src, (char*)&bufLenMessage, 4);
            bufLen = (int)bufLenMessage;
            Vnm_print(2,"Vsh_builtIn: [%d <-- %d] received mesg size=<%d>\n",
                me, src, bufLen);

            /* bufPtr = Vmem_malloc( theePSH->vmem, bufLen, sizeof(char) ); */
            bufPtr = calloc( bufLen, sizeof(char) );

            Vmp_recv(theeVMP, src, bufPtr, bufLen);
            Vnm_print(2,"Vsh_builtIn: [%d <-- %d] received the real mesg.\n",
                me, src);
            theePSH->bufsize = bufLen;
            theePSH->buf = bufPtr;
            rc = 1;
            break;

          case pshcom_vmp_bar:
            Vmp_barr(theeVMP);
            rc = 1;
            break;

          case pshcom_ignore:
            /* truly ignore; keep subsequent layers from trying to execute! */
            rc = 1;
            break;

          default:
            rc = 0;
            break;
        }
    }

    return rc;
}

/*
 * ***************************************************************************
 * Routine:  Vsh_pshell
 *
 * Purpose:  Drop-in replacement for <Vsh_shell> giving parallel extensions.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vsh_pshell(Vsh *thee, char *pPR, void *pthee,
    int (*builtin)(void *thee, int argc, char **argv))
{
    int rc;

    /* don't let shell process args (e.g. MPI wants to use argc/argv itself) */
    thee->processArgs = 0;
    
    /* save the function pointer to be called by our parallel layer */
    theeFunc = builtin;
    theePSH = thee;

    /* create communication object (do first so file i/o gets tagged) */
    theeVMP = Vmp_ctor();

    /* start the vsh, slipping our extra parallel layer in */
    rc = Vsh_shell(thee, pPR, pthee, &PSH_builtin);

    /* destroy communication object */
    Vmp_dtor( &(theeVMP) );

    /* return */
    return rc;
}

