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
 * rcsid="$Id: vpars.c,v 1.7 2008/03/12 05:13:58 fetk Exp $"
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * File:     vpars.c
 *
 * Purpose:  Class Vsh: methods.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */

#include "vsh_p.h"

VEMBED(rcsid="$Id: vpars.c,v 1.7 2008/03/12 05:13:58 fetk Exp $")

VPRIVATE char inbuf[VMAX_BUFSIZE];
VPRIVATE int numRead;

/*
 * ***************************************************************************
 * Routine:  Vsh_parse
 *
 * Purpose:  Parser.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
void Vsh_parse(void)
{
    /* get an input line */
    /* numRead = Vsh_input(inbuf,VMAX_BUFSIZE); */
    VSH_INPUT(inbuf,numRead,VMAX_BUFSIZE);
}

/*
 * ***************************************************************************
 * Routine:  Vsh_parseHandoff
 *
 * Purpose:  Fake parser.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
void Vsh_parseHandoff(char *buf)
{
    strcpy(inbuf,buf);
    numRead = strlen(inbuf);
}

/*
 * ***************************************************************************
 * Routine:  Vsh_execute
 *
 * Purpose:  Executor (execv of a command).
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
void Vsh_execute(void)
{
    char inbuf_argv[VMAX_BUFSIZE], *argv[VMAX_ARGNUM];
    int i, argc;

    /*
     * Init: "cmdKey=0" AFTER possible setJmp/longJmp
     *    cmdKey==0  ==> This WAS NOT a builtin command
     *    cmdKey==1  ==> This WAS a builtin command (non-exit)
     *    cmdKey==2  ==> This WAS the builtin EXIT command
     */
    cmdKey = 0;
    if (numRead == 0) {
        if (Vsh_thee->cinUnit == stdin) {
            Vnm_print(1,"%s",Vsh_thee->PR_EXIT);
            Vnm_print(1,"%s",VNEWLINE_STRING);
        }
        cmdKey = 2;
    } else {
        if (*(inbuf+strlen(inbuf)-1) == VNEWLINE_SYMBOL)
            *(inbuf+strlen(inbuf)-1) = VNULL_SYMBOL;
        strcpy(inbuf_argv,inbuf);
        argc = Vnm_gentokens(inbuf_argv,argv,VMAX_ARGNUM," ","#");

        /*
         * Vnm_print(1,"   TOKENS: (%d)",argc);
         * for (i=0; i<=argc; i++) Vnm_print(1,"   <%s>",argv[i]);
         * Vnm_print(1,"\n");
         */

        if (argc > 0) {

            /* write out the input to our log file */
            /* (We don't want to log "." commands; otherwise running the */
            /* history won't recreate same situation it documented... */
            /* ...so we will put these in the log file as comments... */
            for (i=0; i<argc; i++) {
                if (i==0) {
                    if (!strcmp(argv[0],".")) {
                        Vnm_print(3,"# ");
                    }
                }
                Vnm_print(3,"%s", argv[i]);
                if (i==(argc-1)) {
                    Vnm_print(3,"\n");
                } else {
                    Vnm_print(3," ");
                }
            }
            Vsh_addhist(inbuf,strlen(inbuf));
            cmdKey = Vsh_builtIn(Vsh_thee,argc,argv);
            if (!cmdKey) Vsh_execCmd(Vsh_thee->PR,argc,argv,inbuf);
        }
    }
}

/*
 * ***************************************************************************
 * Routine:  Vsh_yyexecute
 *
 * Purpose:  Fork a child to exec a command, wait for child to finish.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void Vsh_yyexecute(COMMAND *cmd)
{
}

