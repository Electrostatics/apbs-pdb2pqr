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
 * rcsid="$Id: vpup.c,v 1.10 2008/03/12 05:13:58 fetk Exp $"
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * File:     vpup.c
 *
 * Purpose:  Handle the following useful features of a shell like bash:
 *
 *               (1) input/output redirection from/into files
 *               (2) pipes of arbitrary length and complexity
 *               (3) execv/execvp of the actual user commands
 *
 * Notes:    These five routines were twisted somewhat from the examples
 *           in the book "Practical Unix Programming" by Robbins and Robbins.
 *           The five routines are:
 *
 *               VPUBLIC  Vpup_execCmd
 *               VPRIVATE Vpup_makeargv
 *               VPRIVATE Vpup_parseDelim
 *               VPRIVATE Vpup_reDir
 *               VPRIVATE Vpup_connPipe
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */

#include "vsh_p.h"
#include "vpup.h"

VEMBED(rcsid="$Id: vpup.c,v 1.10 2008/03/12 05:13:58 fetk Exp $")

/* define VPUP_SHELL 1 */

#if !defined(HAVE_WINSOCK_H) && defined(VPUP_SHELL)
    VPRIVATE int Vpup_makeargv(char *s, const char *delimiters, char ***argvp);
    VPRIVATE int Vpup_parseDelim(char *s, char delimiter, char **v);
    VPRIVATE int Vpup_reDir(char *infilename, char *outfilename);
    VPRIVATE int Vpup_connPipe(char *cmd, int *frontfd, int *backfd);
#endif

/*
 * ***************************************************************************
 * Class Vsig: Inlineable methods
 * ***************************************************************************
 */
#if !defined(VINLINE_MALOC)

#endif /* if !defined(VINLINE_MALOC) */
/*
 * ***************************************************************************
 * Class Vsig: Non-inlineable methods
 * ***************************************************************************
 */

#if !defined(HAVE_WINSOCK_H) && defined(VPUP_SHELL)
/*
 * ***************************************************************************
 * Routine:  Vpup_execCmd
 *
 * Purpose:  Execute a shell command.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void Vpup_execCmd(const char *PR, int argc, char **argv, char *inbuf)
{
    char **chargv;
    pid_t child_pid;
    char *cmd;
    char *nextcmd;
    int frontfd[2];
    int backfd[2];
  
    frontfd[0] = -1;
    frontfd[1] = -1;
    backfd[0]  = -1;
    backfd[1]  = -1;
    child_pid  = 0;

    if ((nextcmd = inbuf) == VNULL) exit(1);

    while (1) {

        cmd = nextcmd;
        if (cmd == VNULL) break;

        /* if last in pipeline, do not fork another */
        if ((nextcmd = strchr(nextcmd, VPIPE_SYMBOL)) == VNULL) {
            backfd[1] = -1;
            child_pid = 0;

        /* else fork a child to execute next pipeline command */
        } else {
            *nextcmd = VNULL_SYMBOL;
            nextcmd++;
            if (pipe(backfd)== -1) {
                perror(PR);
                exit(1);
            } else if ((child_pid = fork()) == -1) {
                perror(PR);
                exit(1);
            }
        }

        /* the child execs the command */
        if (child_pid == 0) {
            if (Vpup_connPipe(cmd, frontfd, backfd) == -1) {
                perror(PR);
                exit(1);
            } else if (Vpup_makeargv(cmd, VBLANK_STRING, &chargv) > 0) {
                if (execvp(chargv[0], chargv) == -1) {
                    perror(PR);
                    /* Vnm_system("play sorry.au"); */
                }
            }
            exit(1);
        }

        /* the parent closes front pipe and makes back pipe, front */
        close(frontfd[0]);
        close(frontfd[1]);
        frontfd[0] = backfd[0];
        frontfd[1] = backfd[1];
    }
    close(backfd[0]);
    close(backfd[1]);
    exit(1);
}

/*
 * ***************************************************************************
 * Routine:  Vpup_makeargv
 *
 * Purpose:  Create a set of tokens from the input buffer.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPRIVATE int Vpup_makeargv(char *s, const char *delimiters, char ***argvp)
{
    char         *t;
    char         *snew;
    int          numtokens;
    int          i;
    unsigned int size;

    /* snew is real start of string after skipping leading delimiters */
    snew = s + strspn(s, delimiters);

    /* try to create space for a copy of snew in t */
    if ((t = calloc(strlen(snew) + 1, sizeof(char))) == VNULL) {
        *argvp = VNULL;
        numtokens = -1;

    /* if space creation successful, parse the string */
    } else {

        /* count the number of tokens in snew */
        strcpy(t, snew);
        if (strtok(t, delimiters) == VNULL)
            numtokens = 0;
        else
            for (numtokens=1; strtok(VNULL,delimiters) != VNULL; numtokens++);  

        /* create an argument array to contain ptrs to tokens */
        size = (unsigned int)(numtokens + 1);
        if ((*argvp = calloc(size, sizeof(char *))) == VNULL) {
            free(t);
            numtokens = -1;

        /* if successful, insert pointers to tokens into the array */
        } else {
            if (numtokens > 0) {
                strcpy(t, snew);
                **argvp = strtok(t, delimiters);
                for (i=1; i<numtokens+1; i++)
                   *((*argvp)+i) = strtok(VNULL, delimiters);
            } else {
                **argvp = VNULL;
                free(t);
            }
        }
    }   
    return numtokens;
}

/*
 * ***************************************************************************
 * Routine:  Vpup_parseDelim
 *
 * Purpose:  buffer.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPRIVATE int Vpup_parseDelim(char *s, char delimiter, char **v)
{
    char *p;
    char *q;
    int offset;
    int error = 0;

    /* Find position of the delimiting character */
    *v = NULL;
    if ((p = strchr(s, delimiter)) != NULL)  {

        /* Split off the token following delimiter */
        if ((q = (char *)malloc(strlen(p + 1) + 1)) == NULL)
            error = -1;

        else {
            strcpy(q, p + 1);
            if ((*v = strtok(q, VDELIM_SET)) == NULL) error = -1;
            offset = strlen(q);
            strcpy(p, p + offset + 1);
        } 
    }
    return error;
}

/*
 * ***************************************************************************
 * Routine:  Vpup_reDir
 *
 * Purpose:  Do a redirection.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPRIVATE int Vpup_reDir(char *infilename, char *outfilename)
{
    int indes;
    int outdes;

    /* redirect standard in to infilename */
    if (infilename != VNULL) {
        if ((indes = open(infilename, O_RDONLY, VSTDMODE)) == -1)
            return -1;
        if (dup2(indes, STDIN_FILENO) == -1) {
            close(indes);
            return -1;
        }  
        close(indes);
    }

    /* redirect standard out to outfilename */
    if (outfilename != VNULL) {
        if ((outdes=open(outfilename,O_WRONLY|O_CREAT|O_TRUNC,VSTDMODE))==-1)
            return -1;
        if (dup2(outdes, STDOUT_FILENO) == -1) {
            close(outdes);
            return -1;
        }
        close(outdes);
    }   
    return 0;
}

/*
 * ***************************************************************************
 * Routine:  Vpup_connPipe
 *
 * Purpose:  Connect up some pipes.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPRIVATE int Vpup_connPipe(char *cmd, int *frontfd, int *backfd)
{
    int error = 0;
    char *infilename, *outfilename;

    /* look for stdin redirect */
    if (Vpup_parseDelim(cmd, VRDIN_SYMBOL, &infilename) == -1)
        error = -1;

    /* no redirection allowed at front of pipeline */
    else if (infilename != VNULL && frontfd[0] != -1)
        error = -1;

    /* look for stdout redirect */
    else if (Vpup_parseDelim(cmd, VRDOUT_SYMBOL, &outfilename) == -1)
        error = -1;

    /* no redirection allowed at back of pipeline */
    else if (outfilename != VNULL && backfd[1] != -1)
        error = -1;

    /* do any required redirections */
    else if (Vpup_reDir(infilename, outfilename) == -1)
        error = -1;

    /* now connect up appropriate pipes */
    else { 
        if (frontfd[0] != -1) {
            if (dup2(frontfd[0], STDIN_FILENO) == -1)
            error = -1;
        } 
        if (backfd[1] != -1) {
            if (dup2(backfd[1], STDOUT_FILENO) == -1)
                error = -1;
        } 
    }

    /* close unneeded file descriptors */
    close (frontfd[0]);
    close (frontfd[1]);
    close (backfd[0]);
    close (backfd[1]);
    return error;
}
#else

/*
 * ***************************************************************************
 * Routine:  Vpup_execCmd
 *
 * Purpose:  A simple shell command exec.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void Vpup_execCmd(const char *PR, int argc, char **argv, char *inbuf)
{
    Vnm_exec(argc,argv);
}
#endif /* if !defined(HAVE_WINSOCK_H) */

