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
 * rcsid="$Id: vsh.c,v 1.26 2008/03/12 05:13:58 fetk Exp $"
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * File:     vsh.c
 *
 * Purpose:  Class Vsh: methods.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */

#include "vsh_p.h"
#include "vpup.h"

VEMBED(rcsid="$Id: vsh.c,v 1.26 2008/03/12 05:13:58 fetk Exp $")

/* use lex/yacc or not */
#define VSH_LEX_YACC_NOT 1

/* lex/yacc support */
VPUBLIC int cmdKey = 0;
VPUBLIC Vsh *Vsh_thee = VNULL;
VPUBLIC COMMAND *global_command = VNULL;

/* command structure */
typedef enum VSH_command {
    vshcom_none,
    vshcom_clear,
    vshcom_help,
    vshcom_pause,
    vshcom_delay,
    vshcom_set,
    vshcom_penv,
    vshcom_pinfo,
    vshcom_cd,
    vshcom_cdw,
    vshcom_io,
    vshcom_noio,
    vshcom_exit,
    vshcom_dot,
    vshcom_sockg,
    vshcom_sockm,
    vshcom_sockk,
    vshcom_sorry
} VSH_command;

VPRIVATE void Vsh_publishVars(Vsh *thee, int argc, char **argv);
VPRIVATE void Vsh_readlineReset(void);

/*
 * ***************************************************************************
 * Class Vsh: Inlineable methods
 * ***************************************************************************
 */
#if !defined(VINLINE_MALOC)
#endif /* if !defined(VINLINE_MALOC) */

/*
 * ***************************************************************************
 * Class Vsh: Non-inlineable methods
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * Routine:  Vsh_ctor
 *
 * Purpose:  Create the shell.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC Vsh* Vsh_ctor(Vmem *vmem, int argc, char **argv)
{
    Vsh *thee = VNULL;

    VDEBUGIO("Vsh_ctor: CREATING object..");

    thee = Vmem_malloc( VNULL, 1, sizeof(Vsh) );
    if (vmem == VNULL) {
        thee->vmem = Vmem_ctor( "Vsh" );
        thee->iMadeVmem = 1;
    } else {
        thee->vmem = vmem;
        thee->iMadeVmem = 0;
    }

    VDEBUGIO("..done.\n");

    /* start i/o layer */
    Vio_start();

    /* initialization */
    thee->processArgs = 1;
    thee->inUnit      = VNULL;
    thee->scUnit      = VNULL;
    thee->clUnit      = VNULL;
    thee->cinUnit     = VNULL;
    thee->cinName[0]  = '\0';
    thee->PR[0]       = '\0';
    thee->PR_PATH[0]  = '\0';
    strcpy(thee->PR_EXIT,"exit");
    thee->envValuLen = 0;
    thee->envInfoLen = 0;
    thee->envValu = Vmem_malloc(thee->vmem,1,sizeof(char*));
    thee->envInfo = Vmem_malloc(thee->vmem,1,sizeof(char*));
    thee->envValu[0] = VNULL;
    thee->envInfo[0] = VNULL;
    thee->buf        = VNULL;
    thee->bufsize    = 0;

    /* set the builtin/thee pointers */
    thee->Ext_thee    = VNULL;
    thee->Ext_builtin = VNULL;
    Vsh_thee = thee;

    /* publish other variables */
    Vsh_publishVars(thee, argc, argv);

    /* return the object */
    return thee;
}

/*
 * ***************************************************************************
 * Routine:  Vsh_dtor
 *
 * Purpose:  Destroy the shell.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void Vsh_dtor(Vsh **thee)
{
    VASSERT( (*thee) != VNULL );
    if ((*thee) != VNULL) {

        /* wipe the environment */
        Vsh_wipe( (*thee) );

        /* stop i/o layer */
        Vio_stop();

        VDEBUGIO("Vsh_dtor: DESTROYING object..");
        if ((*thee)->iMadeVmem) Vmem_dtor( &((*thee)->vmem) );
        Vmem_free( VNULL, 1, sizeof(Vsh), (void**)thee );
        VDEBUGIO("..done.\n");

        (*thee) = VNULL;
    }
}

/*
 * ***************************************************************************
 * Routine:  Vsh_publishVars
 *
 * Purpose:  Publish environment variables.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPRIVATE void Vsh_publishVars(Vsh *thee, int argc, char **argv)
{
    char homeDirectory[VMAX_ARGLEN];
    char workDirectory[VMAX_ARGLEN];
    char userName[VMAX_ARGLEN];
    char hostName[VMAX_ARGLEN];
    char osName[VMAX_ARGLEN];
    char configFile[VMAX_ARGLEN];
    char buf1[VMAX_ARGLEN], buf2[VMAX_ARGLEN];
    char *term, *cterm;
    int i, numVars = 11;
    typedef struct vshVars {
        char envi[VMAX_ARGLEN];
        char valu[VMAX_ARGLEN];
        char info[VMAX_ARGLEN];
    } vshVars;
    vshVars envVars[] = {
        /* --------   -----       ----------- */
        /* VARIABLE   VALUE       EXPLANATION */
        /* --------   -----       ----------- */
        /* ===[ SOCKET=1 ]=== */
        { "GVLO",     "2x1",
            "socket layout (1x1,...,4x1,1x1s,...,4x1s)" },

        /* ===[ INPUTDEV=5 ]=== */
        { "ISKEY",    "file",
            "VIO INPUT DEV type (sdio,file,buff,unix,inet)" },
        { "ISFMT",    "asc",
            "VIO INPUT DEV format (asc,xdr)" },
        { "IFNAM",    "mcin.m",
            "VIO INPUT DEV file (filename for file I/O)" },
        { "ISNAM",    "0",
            "VIO INPUT DEV name ([ buff | unix | inet ] number)" },
        { "IHVAL",    "localhost",
            "VIO INPUT DEV host (INET hostname or IP address)" },

        /* ===[ OUTPUTDEV=5 ]=== */
        { "OSKEY",    "inet",
            "VIO OUTPUT DEV type (sdio,file,buff,unix,inet)" },
        { "OSFMT",    "asc",
            "VIO OUTPUT DEV format (asc,xdr)" },
        { "OFNAM",    "mcout.m",
            "VIO OUTPUT DEV file (filename for file I/O)" },
        { "OSNAM",    "1",
            "VIO OUTPUT DEV name ([ buff | unix | inet ] number)" },
        { "OHVAL",    "localhost",
            "VIO OUTPUT DEV host (INET hostname or IP address)" },
    };

    /* get username, hostname, and osname */
    VASSERT( Vnm_getuser(userName,sizeof(userName)) );
    VASSERT( Vnm_gethost(hostName,sizeof(hostName)) );
    VASSERT( Vnm_getos(osName,sizeof(osName)) );

    /* get the home directory */
    /*
     * NOTES: the first call to Vnm_gethome fixes the home directory
     *        for all time, regardless of who makes the call.
     */
    VASSERT( Vnm_gethome(homeDirectory,sizeof(homeDirectory)) );

    /* get the working directory */
    /*
     * NOTES: each call to Vnm_getcwd may return a different value;
     *        it returns the current working directory, which may be
     *        modified by calls to Vnm_chdir.
     */
    VASSERT( Vnm_getcwd(workDirectory,sizeof(workDirectory)) );

    /* get some other stuff (may be null) */
    term  = getenv("TERM");
    cterm = getenv("COLORTERM");

    /* config file */
    sprintf(configFile,"%s/%s",homeDirectory,"rc.mcsh");

    /* export to the variables */
    VASSERT( Vsh_putenv(     thee, "USER",      userName               )
          && Vsh_putenvInfo( thee, "USER",      "user name"            ) );
    VASSERT( Vsh_putenv(     thee, "HOSTNAME",  hostName               )
          && Vsh_putenvInfo( thee, "HOSTNAME",  "host name"            ) );
    VASSERT( Vsh_putenv(     thee, "OSTYPE",    osName                 )
          && Vsh_putenvInfo( thee, "OSTYPE",    "operating system"     ) );
    VASSERT( Vsh_putenv(     thee, "HOME",      homeDirectory          )
          && Vsh_putenvInfo( thee, "HOME",      "home directory"       ) );
    VASSERT( Vsh_putenv(     thee, "CWD",       workDirectory          )
          && Vsh_putenvInfo( thee, "CWD",       "working directory"    ) );
    VASSERT( Vsh_putenv(     thee, "TERM",      term                   )
          && Vsh_putenvInfo( thee, "TERM",      "terminal type"        ) );
    VASSERT( Vsh_putenv(     thee, "COLORTERM", cterm                  )
          && Vsh_putenvInfo( thee, "COLORTERM", "color terminal type"  ) );

    /* init file name, shell name, prompt */
    VASSERT( Vsh_putenv(     thee, "ENV",       configFile             )
          && Vsh_putenvInfo( thee, "ENV",       "environ file"         ) );
    VASSERT( Vsh_putenv(     thee, "SHELL",     thee->PR               )
          && Vsh_putenvInfo( thee, "SHELL",     "command shell"        ) );
    VASSERT( Vsh_putenv(     thee, "PROMPT",    thee->PR_PATH          )
          && Vsh_putenvInfo( thee, "PROMPT",    "command shell prompt" ) );

    /* publish remaining variables */
    for (i=0; i<numVars; i++) {
        VASSERT( Vsh_putenv(     thee, envVars[i].envi, envVars[i].valu )
              && Vsh_putenvInfo( thee, envVars[i].envi, envVars[i].info ) );
    }

    /* publish argc/argv variables */
    VASSERT( Vsh_putenvInt(thee,"ARGC",argc)
          && Vsh_putenvInfo(thee,"ARGC","Number of command line parameters") );
    for (i=0; i<argc; i++) {
        sprintf(buf1,"ARG%d",i);
        sprintf(buf2,"Command line parameter <%d>",i);
        VASSERT( Vsh_putenv(     thee, buf1, argv[i] )
              && Vsh_putenvInfo( thee, buf1, buf2    ) );
    }
}

/*
 * ***************************************************************************
 * Routine:  Vsh_trace
 *
 * Purpose:  Trace of token stream.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void Vsh_trace(char *from, char *arg) {
#if defined(VSH_TRACE)
#   if defined(VSH_LEX_YACC)
        fprintf(stderr, "%s token=<%s>, yytext=<%s>\n", from, arg, yytext);
#   else
        fprintf(stderr, "%s token=<%s>\n", from, arg);
#   endif
#endif
}

/*
 * ***************************************************************************
 * Routine:  Vsh_isInteractive
 *
 * Purpose:  Is this an interactive shell.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vsh_isInteractive(Vsh *thee) {
    return ((thee->cinUnit == stdin) && isatty(fileno(stdin)));
}

/*
 * ***************************************************************************
 * Routine:  Vsh_findVar
 *
 * Purpose:  Find a variable in the environment.
 *
 * Notes:    The parameters are:
 *
 *               env    --> the environment variable array
 *               envLen --> the environment variable array length
 *               var    --> the variable we are looking for
 *               term   --> the character terminator (usually "=" or ":")
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vsh_findVar(char **env, int envLen,
    const char *var, const char term)
{
    int i, j, len, ifnd, foundEq;
    char varBuf[VMAX_BUFSIZE];

    /* look for variable in the environment */
    ifnd = -1;
    i = 0;
    while ((ifnd < 0) && (i<envLen)) {

        /* grab the complete string */
        strcpy(varBuf,env[i]);

        /* strip out the variable and the value */
        len = strlen(varBuf);
        foundEq = 0;
        for (j=0; j<len; j++) {
            if (!foundEq) {
                if (varBuf[j] == term) {
                    varBuf[j] = '\0';
                    foundEq = 1;
                }
            } else {
               varBuf[j] = '\0';
            }
        }

        /* now check for match */
        if (!strcmp(varBuf,var)) {
            ifnd = i;
        }
        
        /* next iteration */
        i++;
    }

    return ifnd;
}

/*
 * ***************************************************************************
 * Routine:  Vsh_putenv
 *
 * Purpose:  Place a variable with a value in the environment.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vsh_putenv(Vsh *thee, const char *envi, const char *valu)
{
    int i, len, ifnd;
    char *newValu, **newValuList, valuLoc[VMAX_BUFSIZE];

    VASSERT( envi != VNULL );
    if (valu == VNULL ) {
        valuLoc[0] = '\0';
    } else {
        strcpy(valuLoc,valu);
    }

    /* make the variable=value string */
    len = strlen(envi) + 1 + strlen(valuLoc) + 1;
    newValu = Vmem_malloc(thee->vmem,len,sizeof(char));
    sprintf(newValu,"%s=%s",envi,valuLoc);

    /* look for variable in the environment */
    ifnd = Vsh_findVar(thee->envValu,thee->envValuLen,envi,'=');

    /* if variable exists, just change it */
    if (ifnd >= 0) {

        if (valuLoc[0] != '\0') {

            /* free old VALU */
            len = strlen(thee->envValu[ifnd]) + 1;
            Vmem_free(thee->vmem,len,sizeof(char),
                (void**)&(thee->envValu[ifnd]) );

            /* point to new VALU */
            thee->envValu[ifnd] = newValu;
        }

    /* variable doesn't exist; must create it */
    } else {

        /* expand the environment */
        thee->envValuLen++;

        /* make a new VALU list with one more slot */
        len = thee->envValuLen + 1;
        newValuList = Vmem_malloc(thee->vmem,len,sizeof(char*));

        /* copy the old VALU list into the new VALU list */
        for (i=0; i<thee->envValuLen-1; i++) {
            newValuList[i] = thee->envValu[i];
        }
        newValuList[thee->envValuLen-1] = newValu;
        newValuList[thee->envValuLen] = VNULL;

        /* free old VALU list */
        len = thee->envValuLen;
        Vmem_free(thee->vmem,len,sizeof(char*), (void**)&(thee->envValu) );

        /* setup the new VALU list */
        thee->envValu = newValuList;
    }

    return 1;
}

/*
 * ***************************************************************************
 * Routine:  Vsh_putenvInfo
 *
 * Purpose:  Place a variable with an info string in the environment.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vsh_putenvInfo(Vsh *thee, const char *envi, const char *info)
{
    int i, len, ifnd;
    char *newInfo, **newInfoList, infoLoc[VMAX_BUFSIZE];

    VASSERT( envi != VNULL );
    if (info == VNULL ) {
        infoLoc[0] = '\0';
    } else {
        strcpy(infoLoc,info);
    }

    /* make the variable=info string */
    len = strlen(envi) + 2 + strlen(infoLoc) + 1;
    newInfo = Vmem_malloc(thee->vmem,len,sizeof(char));
    sprintf(newInfo,"%s: %s",envi,infoLoc);

    /* look for variable in the environment */
    ifnd = Vsh_findVar(thee->envInfo,thee->envInfoLen,envi,':');

    /* if variable exists, just change it */
    if (ifnd >= 0) {

        if (infoLoc[0] != '\0') {

            /* free old INFO */
            len = strlen(thee->envInfo[ifnd]) + 1;
            Vmem_free(thee->vmem,len,sizeof(char),
                (void**)&(thee->envInfo[ifnd]) );

            /* point to new INFO */
            thee->envInfo[ifnd] = newInfo;
        }

    /* variable doesn't exist; must create it */
    } else {

        /* expand the environment */
        thee->envInfoLen++;

        /* make a new INFO list with one more slot */
        len = thee->envInfoLen + 1;
        newInfoList = Vmem_malloc(thee->vmem,len,sizeof(char*));

        /* copy the old INFO list into the new INFO list */
        for (i=0; i<thee->envInfoLen-1; i++) {
            newInfoList[i] = thee->envInfo[i];
        }
        newInfoList[thee->envInfoLen-1] = newInfo;
        newInfoList[thee->envInfoLen] = VNULL;

        /* free old INFO list */
        len = thee->envInfoLen;
        Vmem_free(thee->vmem,len,sizeof(char*), (void**)&(thee->envInfo) );

        /* setup the new INFO list */
        thee->envInfo = newInfoList;
    }

    return 1;
}

/*
 * ***************************************************************************
 * Routine:  Vsh_putenvInt
 *
 * Purpose:  Place a variable with a value (integer) in the environment.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vsh_putenvInt(Vsh *thee, const char *envi, const int valu)
{
    char buf[VMAX_BUFSIZE];

    sprintf(buf,"%d",valu);
    Vsh_putenv(thee,envi,buf);

    return 1;
}

/*
 * ***************************************************************************
 * Routine:  Vsh_putenvReal
 *
 * Purpose:  Place a variable with a value (real) in the environment.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vsh_putenvReal(Vsh *thee, const char *envi, const double valu)
{
    char buf[VMAX_BUFSIZE];

    sprintf(buf,"%e",valu);
    Vsh_putenv(thee,envi,buf);

    return 1;
}


/*
 * ***************************************************************************
 * Routine:  Vsh_getenv
 *
 * Purpose:  Get a value of variable in the environment.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC char* Vsh_getenv(Vsh *thee, const char *envi)
{
    int ifnd;

    ifnd = Vsh_findVar(thee->envValu,thee->envValuLen,envi,'=');
    if (ifnd >= 0) {
        return (thee->envValu[ifnd]+strlen(envi)+1);
    } else {
        return VNULL;
    }
}

/*
 * ***************************************************************************
 * Routine:  Vsh_getenvInfo
 *
 * Purpose:  Get info associated with a variable in the environment.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC char* Vsh_getenvInfo(Vsh *thee, const char *envi)
{
    int ifnd;

    ifnd = Vsh_findVar(thee->envInfo,thee->envInfoLen,envi,':');
    if (ifnd >= 0) {
        return (thee->envInfo[ifnd]+strlen(envi)+2);
    } else {
        return VNULL;
    }
}

/*
 * ***************************************************************************
 * Routine:  Vsh_getenvInt
 *
 * Purpose:  Get a value of variable in the environment as an integer.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vsh_getenvInt(Vsh *thee, const char *envi)
{
    int ifnd;

    ifnd = Vsh_findVar(thee->envValu,thee->envValuLen,envi,'=');
    if (ifnd >= 0) {
        return atoi(thee->envValu[ifnd]+strlen(envi)+1);
    } else {
        return 0;
    }
}

/*
 * ***************************************************************************
 * Routine:  Vsh_getenvReal
 *
 * Purpose:  Get a value of variable in the environment as a real.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC double Vsh_getenvReal(Vsh *thee, const char *envi)
{
    int ifnd;

    ifnd = Vsh_findVar(thee->envValu,thee->envValuLen,envi,'=');
    if (ifnd >= 0) {
        return atof(thee->envValu[ifnd]+strlen(envi)+1);
    } else {
        return 0.0;
    }
}

/*
 * ***************************************************************************
 * Routine:  Vsh_remove
 *
 * Purpose:  Remove a variable from the environment.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void Vsh_remove(Vsh *thee, const char *envi)
{
    /* unsetenv(envi); */
}

/*
 * ***************************************************************************
 * Routine:  Vsh_wipe
 *
 * Purpose:  Wipe the environment.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void Vsh_wipe(Vsh *thee)
{
    int i, len;

    VASSERT( thee->envValu != VNULL );

    /* wipe the entire environment */
    for (i=0; i<thee->envValuLen; i++) {
        len = strlen(thee->envValu[i]) + 1;
        Vmem_free(thee->vmem,len,sizeof(char), (void**)&(thee->envValu[i]) );
    }
    len = thee->envValuLen + 1;
    Vmem_free(thee->vmem,len,sizeof(char*), (void**)&(thee->envValu) );
    for (i=0; i<thee->envInfoLen; i++) {
        len = strlen(thee->envInfo[i]) + 1;
        Vmem_free(thee->vmem,len,sizeof(char), (void**)&(thee->envInfo[i]) );
    }
    len = thee->envInfoLen + 1;
    Vmem_free(thee->vmem,len,sizeof(char*), (void**)&(thee->envInfo) );
}

/*
 * ***************************************************************************
 * Routine:  Vsh_printenv
 *
 * Purpose:  Print the environment.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void Vsh_printenv(Vsh *thee)
{
    int i;

    for (i=0; i<thee->envValuLen; i++) {
        Vnm_print(1,"%s\n",thee->envValu[i]);
    }
}

/*
 * ***************************************************************************
 * Routine:  Vsh_printenvInfo
 *
 * Purpose:  Print the environment info.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void Vsh_printenvInfo(Vsh *thee)
{
    int i;

    for (i=0; i<thee->envInfoLen; i++) {
        Vnm_print(1,"%s\n",thee->envInfo[i]);
    }
}

/*
 * ***************************************************************************
 * Routine:  Vsh_completion
 *
 * Purpose:  Command completion action.
 *
 * Notes:    The return type and argument list is mandated by readline:
 *
 *               int func(int count, int key)
 *
 *           Useful tips: the readline library uses the following conventions
 *           for simplifying the notation for function pointers:
 *          
 *               typedef int IFunction ();
 *               typedef void VFunction ();
 *               typedef char *CPFunction ();
 *               typedef char **CPPFunction ();
 *
 *           This allows one to replace something like:
 *
 *               int (*)()func;
 *
 *           with simply:
 *
 *               IFunction *func;
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
#if defined(HAVE_READLINE_READLINE_H)
VPRIVATE int Vsh_completion(int count, int key)
{
    int cargc = 1; 
    char *cargv[2] = { "help", VNULL }; 

    Vnm_print(1,"\n");
    Vsh_builtIn(Vsh_thee, cargc, cargv);
    Vnm_print(1,"%s",Vsh_thee->PR_PATH);

    Vsh_readlineReset();

    return 0;
}
#endif

/*
 * ***************************************************************************
 * Routine:  Vsh_readlineInit
 *
 * Purpose:  Initialize the readline library.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPRIVATE void Vsh_readlineInit(void)
{
    static int init=0;

    if (!init) {
        init = 1;

#if defined(HAVE_READLINE_READLINE_H)
#if 0
        rl_catch_signals = 0;
        rl_catch_sigwinch = 0;
#endif
        rl_bind_key(9, &Vsh_completion);
#endif
    }
}

/*
 * ***************************************************************************
 * Routine:  Vsh_readlineReset
 *
 * Purpose:  Reset the readline library (e.g. free partial input string).
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPRIVATE void Vsh_readlineReset(void)
{
    Vsh_readlineInit();

#if defined(HAVE_READLINE_READLINE_H)
#if 0
    rl_free_line_state();
    rl_resize_terminal();
#endif
#endif
}

/*
 * ***************************************************************************
 * Routine:  Vsh_addhist
 *
 * Purpose:  Put an input line into the history list.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void Vsh_addhist(char *buf, int buflen)
{
    Vsh_readlineInit();

#if defined(HAVE_READLINE_HISTORY_H)
    add_history(buf);
#endif
}

/*
 * ***************************************************************************
 * Routine:  Vsh_readline
 *
 * Purpose:  Get an input line.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC char *Vsh_readline(char *prompt, char *buf, int buflen, FILE *stream)
{
    char *key;

    if (stream != stdin) {
        key = fgets(buf, buflen, stream);
    } else {

#if defined(HAVE_READLINE_READLINE_H)
        Vsh_readlineInit();
        key = readline(prompt);
        if (key == VNULL) {
            buf[0] = '\n';
            buf[1] = '\0';
        } else if (key[0] == '\0') {
            buf[0] = '\n';
            buf[1] = '\0';
            free(key);
        } else {
            strncpy(buf,key,buflen);
            free(key);
        }
#else
        Vnm_print(1,"%s",prompt);
        key = fgets(buf, buflen, stream);
#endif

    }

    return key;
}

/*
 * ***************************************************************************
 * Routine:  Vsh_shell
 *
 * Purpose:  A bash-like shell with user-definable extensions.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vsh_shell(Vsh *thee, char *pPR, void *pthee,
    int (*builtin)(void *thee, int argc, char **argv))
{
    int i, argc, offset;
    char **argvPtr, buf[VMAX_ARGLEN];
    char *argvNULL = "\0";
    struct stat fInfo;
 
    /* we will need argc and argv[] */
    argc = Vsh_getenvInt(thee, "ARGC");
    argvPtr = Vmem_malloc(thee->vmem, VMAX_ARGLEN, sizeof(char*));
    for (i=0; i<argc; i++) {
        sprintf(buf,"ARG%d",i);
        argvPtr[i] = Vsh_getenv(thee, buf);
    }
    argvPtr[argc] = argvNULL;

    /* paranoia: check type sizes on the machine */
    Vnm_typeChk();

    /* the externally supplied builtin object pointer and function */
    thee->Ext_thee    = pthee;
    thee->Ext_builtin = builtin;

    /* construct a reasonable shell command prompt if not given as argument */
    if (pPR != VNULL) {
        if (pPR[0] != '\0') {
            sprintf(thee->PR,"%s",pPR);
        }
    }
    if (thee->PR[0] == '\0') {
        VASSERT( argc > 0 );
        strncpy(buf,argvPtr[0],VMAX_ARGLEN);
        offset = 0;
        if (strlen(buf) >= 2) {
            /* remove the "./" if it is there */
            if ((buf[0] == '.') && (buf[1] == '/')) {
                offset = 2;
            /* remove "-" if there; happens when vsh is a login shell */
            } else if (buf[0] == '-') {
                offset = 1;
            }
        }
        sprintf(thee->PR,"%s",buf+offset);
    }
    Vsh_putenv(thee,"SHELL",thee->PR);

    /*
     * if filename given on command line, try to take input from there.
     */
    thee->inUnit = stdin;
    if (thee->processArgs) {
        for (i=1; i<argc; i++) {

            /* is the argument the "-h" option */
            if (!strcmp(argvPtr[i],"-h")) {
                VJMPERR1(1);

            /* is the argument the "-io" option */
            } else if (!strcmp(argvPtr[i],"-io")) {
                Vnm_redirect(0);

            /* is the argument the "-noio" option */
            } else if (!strcmp(argvPtr[i],"-noio")) {
                Vnm_redirect(1);

            /* try to open the argument as a script file */
            } else {
                thee->clUnit = fopen(argvPtr[i], "r");
                if (thee->clUnit == VNULL) {
                    Vnm_print(2,"%s: Problem opening file <%s>\n",
                        thee->PR, argvPtr[i]);
                    thee->inUnit = stdin;
                    VJMPERR1(1);
                } else thee->inUnit = thee->clUnit;
            }
        }
    }

    /* the current input unit starts out as stdin or a script */
    thee->cinUnit = thee->inUnit;

    /* we first execute the user's configuration file */
    if ( !stat(Vsh_getenv(thee,"ENV"), &fInfo) ) {
        if (VS_ISREG(fInfo.st_mode)) {
            thee->scUnit = fopen(Vsh_getenv(thee,"ENV"), "r");
            if (thee->scUnit != VNULL) {
                thee->cinUnit = thee->scUnit;
                strncpy(thee->cinName,Vsh_getenv(thee,"ENV"),VMAX_ARGLEN);
                Vnm_print(0,"Starting <%s> script named <%s>\n",
                    thee->PR,thee->cinName);
            }
        }
    }

    /* start the command shell parsing loop */
    cmdKey = 0;
    while (cmdKey != 2) {

        /*
         * Parse a single input unit.
         * An input unit is a complete shell statement
         * (e.g. one-line command, if-then-else, while, case, etc)
         * which may span multiple lines of input.
         */
#if defined(VSH_LEX_YACC)
        yyparse();  /* for complex Bourne-shell compatible input units */
        Vsh_yyexecute(global_command);
#else
        Vsh_parse();  /* for simple one-line commands only */
        Vsh_execute();
#endif

    }

    /* close the input unit for a command-line file if we had one */
    if (thee->clUnit != VNULL) VASSERT( !fclose(thee->clUnit) );

    /* free the argv storage */
    Vmem_free( thee->vmem, VMAX_ARGLEN, sizeof(char*), (void**)&argvPtr );

    /* no error */
    return 1;

  VERROR1:
    /* close the input unit for a command-line file if we had one */
    if (thee->clUnit != VNULL) VASSERT( !fclose(thee->clUnit) );

    /* free the argv storage */
    Vmem_free( thee->vmem, VMAX_ARGLEN, sizeof(char*), (void**)argvPtr );

    /* print a usage menu */
    Vnm_print(2,"usage: ./%s [ -h | -io | -noio ]"
                " [ shellScriptFile ]\n",thee->PR);

    /* error, but not fatal */
    return 1;
}

/*
 * ***************************************************************************
 * Routine:  Vsh_input
 *
 * Purpose:  Read a single newline-terminated line of input.
 *
 * Notes:    We output a prompt if appropriate to do so.
 *           We catch any SIGINTs generated during input just like
 *           any normal command shell.  If one is caught, we jump
 *           back to the input prompt via setjmp/longjmp.
 *           We also handle killing of scripts via SIGINT, and
 *           correctly handle prompt output on redirection.
 *
 *           This routine is used with the following macro:
 *
 *               define VSH_INPUT(buf,result,max_size) { \  
 *                   result = Vsh_input(buf,max_size); \  
 *               }
 *
 *           which has exactly the same functionality as the YY_INPUT
 *           macro used by LEX.  Forcing LEX to use VSH_INPUT in place
 *           of YY_INPUT allows all of the above to work for a
 *           LEX/YACC-generated command shell parser.
 *
 * Input:    buf[0..buflen-1] = character buffer to read in to
 *           buflen           = length of buf
 *
 * Output:   returned = number of characters actually read in.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vsh_input(char *buf, int buflen)
{
    int numRead;
    char *key, currDirectory[VMAX_ARGLEN];
    jmp_buf *jbuf;

    /* setup for the jump */
    jbuf = Vnm_signalInit();
    if (setjmp(*jbuf)) {

        /* FIRST: kill any script that may have been executing */
        Vsh_thee->cinUnit = Vsh_thee->inUnit;

        /* SECOND: reset the readline input buffer */
        Vsh_readlineReset();

        /* THIRD: restart lex on the (possibly) new input unit */
#if defined(VSH_LEX_YACC)
            yyrestart(Vsh_thee->cinUnit);
#endif

        /* FOURTH: print a newline if we are in interactive mode */
        if (Vsh_isInteractive(Vsh_thee)) {
            Vnm_print(1,"%s",VNEWLINE_STRING);
        }
    }

    /* close script unit if open -- must have killed it */
    if ( ((Vsh_thee->cinUnit != Vsh_thee->scUnit) 
         || feof(Vsh_thee->scUnit)) && (Vsh_thee->scUnit != VNULL) ) {
        VASSERT( !fclose(Vsh_thee->scUnit) );
        Vsh_thee->scUnit = VNULL;
        Vnm_print(0,"Stopping <%s> script named <%s>\n",
            Vsh_thee->PR,Vsh_thee->cinName);
        strncpy(Vsh_thee->cinName," ",VMAX_ARGLEN);
    }

    /* OKAY-TO-JUMP back to the input prompt now */
    Vnm_jmpOkSet();

    /* deal with different I/O units */
    if (Vsh_isInteractive(Vsh_thee)) {
        VASSERT( Vnm_getcwd(currDirectory,sizeof(currDirectory)) );
        sprintf(Vsh_thee->PR_PATH,"%s@%s<%s+%s>%% ",
            Vsh_getenv(Vsh_thee,"USER"),
            Vsh_getenv(Vsh_thee,"HOSTNAME"),
            Vsh_getenv(Vsh_thee,"OSTYPE"),
            Vsh_thee->PR);
        Vsh_putenv(Vsh_thee,"PROMPT",Vsh_thee->PR_PATH);
    }

    /* get the input line */
    memset(buf, VNULL_SYMBOL, buflen);
    key = Vsh_readline(Vsh_thee->PR_PATH, buf, buflen, Vsh_thee->cinUnit);
    if ((key == VNULL) && (Vsh_thee->cinUnit==Vsh_thee->scUnit)) {

        /* shut down the script handling */
        VASSERT( Vsh_thee->scUnit != VNULL );
        VASSERT( feof(Vsh_thee->scUnit) );
        VASSERT( !fclose(Vsh_thee->scUnit) );
        Vsh_thee->scUnit = VNULL;
        Vnm_print(0,"Stopping <%s> script named <%s>\n",
            Vsh_thee->PR,Vsh_thee->cinName);
        strncpy(Vsh_thee->cinName," ",VMAX_ARGLEN);

        /* deal with different I/O units */
        Vsh_thee->cinUnit = Vsh_thee->inUnit;
        if (Vsh_isInteractive(Vsh_thee)) {
            VASSERT( Vnm_getcwd(currDirectory,sizeof(currDirectory)) );
            sprintf(Vsh_thee->PR_PATH,"%s@%s<%s+%s>%% ",
                Vsh_getenv(Vsh_thee,"USER"),
                Vsh_getenv(Vsh_thee,"HOSTNAME"),
                Vsh_getenv(Vsh_thee,"OSTYPE"),
                Vsh_thee->PR);
            Vsh_putenv(Vsh_thee,"PROMPT",Vsh_thee->PR_PATH);
        }

        /* get the input line */
        memset(buf, VNULL_SYMBOL, buflen);
        key = Vsh_readline(Vsh_thee->PR_PATH,
            buf, buflen, Vsh_thee->cinUnit);
    }

    /* calculate the number of characters actually read */
    if (key == VNULL) {
        numRead = 0;
    } else {
        numRead = strlen(buf);
    }
    VASSERT( numRead <= buflen );

    /* NOT-OKAY-TO-JUMP back to the input prompt now */
    Vnm_jmpOkClear();

    /* clear the numerical loop signal before executing the command */
    Vnm_sigIntClear();

    /* return num chars read */
    return numRead;
}

/*
 * ***************************************************************************
 * Routine:  Vsh_getCmd
 *
 * Purpose:  Decode the input string into a legal command.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPRIVATE VSH_command Vsh_getCmd(int argc, char **argv)
{
    VSH_command theCmd = vshcom_none;
    if (!strcmp(argv[0],"")) {
        theCmd = vshcom_none;
    } else if ( (!strcmp(argv[0],"c")) || (!strcmp(argv[0],"clear")) ) { 
        theCmd = vshcom_clear;
    } else if (!strcmp(argv[0],"help")) { 
        theCmd = vshcom_help;
    } else if (!strcmp(argv[0],"pause")) {
        theCmd = vshcom_pause;
    } else if (!strcmp(argv[0],"delay")) {
        theCmd = vshcom_delay;
    } else if (!strcmp(argv[0],"set")) {
        theCmd = vshcom_set;
    } else if (!strcmp(argv[0],"penv")) {
        theCmd = vshcom_penv;
    } else if (!strcmp(argv[0],"pinfo")) {
        theCmd = vshcom_pinfo;
    } else if (!strcmp(argv[0],"cd")) {
        theCmd = vshcom_cd;
    } else if (!strcmp(argv[0],"cdw")) {
        theCmd = vshcom_cdw;
    } else if (!strcmp(argv[0],"io")) {
        theCmd = vshcom_io;
    } else if (!strcmp(argv[0],"noio")) {
        theCmd = vshcom_noio;
    } else if (!strcmp(argv[0],"exit")) {
        theCmd = vshcom_exit;
    } else if (!strcmp(argv[0],".")) {
        theCmd = vshcom_dot;
    } else if (!strcmp(argv[0],"sockg")) {
        theCmd = vshcom_sockg;
    } else if (!strcmp(argv[0],"sockm")) {
        theCmd = vshcom_sockm;
    } else if (!strcmp(argv[0],"sockk")) {
        theCmd = vshcom_sockk;
    } else {
        theCmd = vshcom_none;
    }
    return theCmd;
}

/*
 * ***************************************************************************
 * Routine:  Vsh_execCmd
 *
 * Purpose:  Fork a child to exec a command, wait for child to finish.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void Vsh_execCmd(const char *PR, int argc, char **argv, char *inbuf)
{
    /* fork a child to do the work (real shell behavior) */
#if !defined(HAVE_WINSOCK_H)

    static pid_t child_pid;
    char PR_TMP[VMAX_ARGLEN];

    VASSERT( argc > 0 );
    sprintf(PR_TMP,"%s: %s",PR,argv[0]);

    if ((child_pid=fork()) == 0) {
        /* NOTE: child should NOT return except on error */
        Vpup_execCmd(PR_TMP,argc,argv,inbuf);
        perror(PR_TMP);
        exit(1);
    } else if (child_pid > 0) {
        wait(VNULL);
    } else {
        perror(PR_TMP);
    }

    /* fork() does not exist on Win32 */
    /* (also fork() is BROKEN in Linux 2.1.85 believe it or not...) */
    /* SO, we fake it by passing command to underlying shell -- bummer! */
#else

    Vnm_system(inbuf);

#endif
}

/*
 * ***************************************************************************
 * Routine:  Vsh_memChk
 *
 * Purpose:  Print the exact current malloc usage.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void Vsh_memChk(Vsh *thee)
{
    if (thee->iMadeVmem) Vmem_print(thee->vmem);
}

/*
 * ***************************************************************************
 * Routine:  Vsh_killSockets
 *
 * Purpose:  Kill any socket graphics.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPRIVATE void Vsh_killSockets(Vsh *thee)
{
    Vnm_systemKill("gvx");
    Vnm_systemKill("mcsg");
    Vnm_systemKill("mcbridge");
    Vnm_system("sleep 1");
}

/*
 * ***************************************************************************
 * Routine:  Vsh_builtIn
 *
 * Purpose:  Vsh_shell built-in commands.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC int Vsh_builtIn(Vsh *thee, int argc, char **argv)
{
    int i, rc;
    VSH_command theCmd;
    struct stat fInfo;
    jmp_buf *jbuf;

    char sysc[VMAX_BUFSIZE];
    char sofmt[VMAX_ARGLEN], sokey[VMAX_ARGLEN];

    static int  init=0;
    static char buf[VMAX_BUFSIZE], workDirectory[VMAX_ARGLEN];
    static char PR_TMP[VMAX_ARGLEN];
    static char env[VMAX_BUFSIZE], com[VMAX_BUFSIZE];
    static char sock[VMAX_BUFSIZE];
    const char *stmp;

    /* one-time intialization */
    if (!init) {
        init=1;

        /* make the env message (%s slots = 12) */
        stmp =
          "%s: Execution environment, directories, and files:\n"
          "    Shell                --> <%s>\n"
          "    User name            --> <%s>\n"
          "    Host name            --> <%s>\n"
          "    Operating system     --> <%s>\n"
          "    Home directory       --> <%s>\n"
          "    Work directory       --> <%s>\n"
          "    Startup script       --> <%s/%s>\n"
          "    Command history file --> <%s/%s>\n"
          "    I/O capture file     --> <%s/%s>\n";
        sprintf(env,stmp,thee->PR,
            thee->PR,
            Vsh_getenv(thee,"USER"),
            Vsh_getenv(thee,"HOSTNAME"),
            Vsh_getenv(thee,"OSTYPE"),
            Vsh_getenv(thee,"HOME"),
            Vsh_getenv(thee,"CWD"),
            Vsh_getenv(thee,"HOME"),"rc.mcsh",
            Vsh_getenv(thee,"HOME"),"hist.mcsh",
            Vsh_getenv(thee,"HOME"),"io.mc");

        /* make the com message (%s slots = 2) */
        stmp =
          "%s: Shell interaction commands: \n"
          "    help [ env|com|sock ] --> print help messages\n"
          "    c | clear                 --> clear the screen\n"
          "    pause                     --> pause until carriage return\n"
          "    delay                     --> delay for three seconds\n"
          "    cd | cdw                  --> cd to home or work directory\n";
        sprintf(com,stmp,thee->PR);
        stmp =
          "    io | noio                 --> display extra io or not\n"
          "    set [var [val]]           --> set or print variables\n"
          "    penv [var]                --> print variable\n"
          "    pinfo [var]               --> print variable information\n"
          "    . scriptfile              --> execute an vsh scriptfile\n"
          "    exit | CTRL-D             --> exit the <%s> shell\n";
        sprintf(buf,stmp,thee->PR);
        strcat(com,buf);

        /* make the socket message (%s slots = 1) */
        stmp = "%s: Socket Graphics manipulation commands: \n"
            "    sockk --> kill all socket graphics processes\n"
            "    sockg --> start geomview socket displays\n"
            "    sockm --> start mcsg socket displays\n";
        sprintf(sock,stmp,thee->PR);
    }

    /* the user-defined shell gets first shot at the command */
    if (thee->Ext_builtin != VNULL) {
        rc = (*(thee->Ext_builtin))(thee->Ext_thee,argc,argv);
        if (rc != 0) return rc;
    }

    /* now it is our turn; set default return code (success) */
    rc = 1;

    /* get the command */
    theCmd = Vsh_getCmd(argc, argv);

    /* decode and execute the command */
    switch (theCmd) {

      case vshcom_pause:
        sprintf(PR_TMP,"%s: Press <return> to continue..", thee->PR);
        memset(buf, VNULL_SYMBOL, sizeof(buf));

        /* setup for the jump */
        jbuf = Vnm_signalInit();
        if (setjmp(*jbuf)) {

            /* FIRST: kill any script that may have been executing */
            thee->cinUnit = thee->inUnit;

            /* SECOND: reset the readline input buffer */
            Vsh_readlineReset();

            /* THIRD: restart lex on the (possibly) new input unit */
#if defined(VSH_LEX_YACC)
                yyrestart(thee->cinUnit);
#endif

            /* FOURTH: print a newline if we are in interactive mode */
            if (Vsh_isInteractive(thee)) {
                Vnm_print(1,"%s",VNEWLINE_STRING);
            }

            Vnm_jmpOkClear();
            Vnm_sigIntClear();
        } else {
            Vnm_jmpOkSet();
            Vsh_readline(PR_TMP, buf, sizeof(buf), stdin);
            Vnm_jmpOkClear();
            Vnm_sigIntClear();
        }
        rc = 1;
        break;

      case vshcom_delay:
        Vnm_sleep(3000000);
        break;

      case vshcom_clear:
        memset(buf, VNULL_SYMBOL, sizeof(buf));
#if !defined(HAVE_WINSOCK_H)
        strcpy(buf,"clear");
#else
        strcpy(buf,"cls");
#endif
        Vnm_system(buf);
        break;

      case vshcom_help:
        if (argc==1) {
            Vnm_print(1,"%s: Vsh-layer Help Menu:\n",thee->PR);
            Vnm_print(1,"    help env  --> Help on %s environment\n",
                thee->PR);
            Vnm_print(1,"    help com  --> Help on %s shell commands\n",
                thee->PR);
            Vnm_print(1,"    help sock --> Help on %s socket graphics\n",
                thee->PR);
        } else {
            if (!strcmp(argv[1],"env")) { 
                Vnm_print(1, "%s", env);
            } else if (!strcmp(argv[1],"com")) { 
                Vnm_print(1, "%s", com);
            } else if (!strcmp(argv[1],"sock")) { 
                Vnm_print(1, "%s", sock);
            } else {
                /* we are the last shell layer; nothing to defer to */
            }
        }
        break;

      case vshcom_set:
        if (argc <= 1) {
            Vsh_printenv(thee);
        } else if (argc == 2) {
            Vnm_print(2,"%s: %s=%s\n",
                thee->PR,argv[1],Vsh_getenv(thee,argv[1]));
        } else if (argc == 3) {
            Vsh_putenv(thee,argv[1],argv[2]);
        } else if (argc > 3) {
            Vsh_putenv(thee,argv[1],argv[2]);
            buf[0] = '\0';
            for (i=3; i<argc; i++) {
                if (i>3) strcat(buf," ");
                strcat(buf,argv[i]);
            }
            Vsh_putenvInfo(thee,argv[1],buf);
        } else VASSERT(0);
        break;

      case vshcom_penv:
        if (argc <= 1) {
            Vsh_printenv(thee);
        } else if (argc == 2) {
            Vnm_print(2,"%s: %s=%s\n",
                thee->PR,argv[1],Vsh_getenv(thee,argv[1]));
        } else Vnm_print(2,"%s: %s: Too many arguments\n",thee->PR,argv[0]);
        break;

      case vshcom_pinfo:
        if (argc <= 1) {
            Vsh_printenvInfo(thee);
        } else if (argc == 2) {
            Vnm_print(2,"%s: %s: %s\n",
                thee->PR,argv[1],Vsh_getenvInfo(thee,argv[1]));
        } else Vnm_print(2,"%s: %s: Too many arguments\n",thee->PR,argv[0]);
        break;

      case vshcom_cd:
        if (argc==1) {
            if (Vnm_chdir(Vsh_getenv(thee,"HOME")) == -1) {
                Vnm_print(2,"%s: %s: %s: No such directory\n",
                    thee->PR, argv[0], Vsh_getenv(thee,"HOME"));
            }
            VASSERT( Vnm_getcwd(workDirectory,sizeof(workDirectory)) );
            Vsh_putenv(thee,"CWD", workDirectory);
        } else {
            if (Vnm_chdir(argv[1]) == -1) {
                Vnm_print(2,"%s: %s: %s: No such directory\n",
                    thee->PR, argv[0], argv[1]);
            }
            VASSERT( Vnm_getcwd(workDirectory,sizeof(workDirectory)) );
            Vsh_putenv(thee,"CWD", workDirectory);
        }
        break;

      case vshcom_cdw:
        if (argc==1) {
            if (Vnm_chdir(Vsh_getenv(thee,"MCSH_HOME")) == -1) {
                Vnm_print(2,"%s: %s: %s: No such directory\n",
                    thee->PR, argv[0], Vsh_getenv(thee,"MCSH_HOME"));
            }
            VASSERT( Vnm_getcwd(workDirectory,sizeof(workDirectory)) );
            Vsh_putenv(thee,"CWD", workDirectory);
        } else {
            if (Vnm_chdir(argv[1]) == -1) {
                Vnm_print(2,"%s: %s: %s: No such directory\n",
                    thee->PR, argv[0], argv[1]);
            }
            VASSERT( Vnm_getcwd(workDirectory,sizeof(workDirectory)) );
            Vsh_putenv(thee,"CWD", workDirectory);
        }
        break;

      case vshcom_io:
        Vnm_redirect(0);
        break;

      case vshcom_noio:
        Vnm_redirect(1);
        break;

      case vshcom_exit:
        if (Vnm_chdir(Vsh_getenv(thee,"HOME")) == -1)
            Vnm_print(2,"%s: %s: %s: No such directory\n",
                thee->PR, "cd", Vsh_getenv(thee,"HOME"));
        rc = 2;                                                           
        break;

      case vshcom_dot:
        if (argc <= 1) {
            Vnm_print(2,"%s: Filename argument required\n", thee->PR);
        } else if ( !stat(argv[1], &fInfo) ) {
            if (VS_ISREG(fInfo.st_mode)) {
                thee->scUnit = fopen(argv[1], "r");
                if (thee->scUnit != VNULL) {
                    thee->cinUnit = thee->scUnit;
                    strncpy(thee->cinName,argv[1],VMAX_ARGLEN);
                    Vnm_print(0,"Starting <%s> script named <%s>\n",
                        thee->PR,thee->cinName);
                } else {
                    Vnm_print(2,"%s: Problem opening <%s>\n",
                        thee->PR,argv[1]);
                }
            } else {
                Vnm_print(2,"%s: File <%s> not normal\n",
                    thee->PR,argv[1]);
            }
        } else {
            Vnm_print(2,"%s: File <%s> not found\n",thee->PR,argv[1]);
        }
        break;

      case vshcom_sockg:
        Vsh_killSockets(thee);
        strncpy(sofmt,Vsh_getenv(thee,"OSFMT"),VMAX_ARGLEN);
        if (!strcmp("unix",Vsh_getenv(thee,"OSKEY"))) {
            strncpy(sokey,"Mcs",VMAX_ARGLEN);
        } else if (!strcmp("inet",Vsh_getenv(thee,"OSKEY"))) {
            strncpy(sokey,"Mcs",VMAX_ARGLEN);
        } else {
            strncpy(sokey,"Mcs",VMAX_ARGLEN);
        }
        if (!strcmp("1x1",Vsh_getenv(thee,"GVLO"))) {
            sprintf(sysc,"geomview -nopanels -wpos 436,445@55,520 -%s 0",sokey);
            Vnm_systemBack(sysc);
            if (!strcmp("inet",Vsh_getenv(thee,"OSKEY"))) {
                Vnm_systemBack("mcbridge -i2u 0 0");
            }
        } else if (!strcmp("2x1",Vsh_getenv(thee,"GVLO"))) {
            sprintf(sysc,"geomview -nopanels -wpos 436,445@55,520 -%s 0",sokey);
            Vnm_systemBack(sysc);
            sprintf(sysc,"geomview -nopanels -wpos 436,445@55,50  -%s 1",sokey);
            Vnm_systemBack(sysc);
            if (!strcmp("inet",Vsh_getenv(thee,"OSKEY"))) {
                Vnm_systemBack("mcbridge -i2u 0 0");
                Vnm_systemBack("mcbridge -i2u 1 1");
            }
        } else if (!strcmp("3x1",Vsh_getenv(thee,"GVLO"))) {
            sprintf(sysc,"geomview -nopanels -wpos 436,287@55,684 -%s 0",sokey);
            Vnm_systemBack(sysc);
            sprintf(sysc,"geomview -nopanels -wpos 436,287@55,367 -%s 1",sokey);
            Vnm_systemBack(sysc);
            sprintf(sysc,"geomview -nopanels -wpos 436,287@55,50  -%s 2",sokey);
            Vnm_systemBack(sysc);
            if (!strcmp("inet",Vsh_getenv(thee,"OSKEY"))) {
                Vnm_systemBack("mcbridge -i2u 0 0");
                Vnm_systemBack("mcbridge -i2u 1 1");
                Vnm_systemBack("mcbridge -i2u 2 2");
            }
        } else if (!strcmp("4x1",Vsh_getenv(thee,"GVLO"))) {
            sprintf(sysc,"geomview -nopanels -wpos 255,255 -%s 0",sokey);
            Vnm_systemBack(sysc);
            sprintf(sysc,"geomview -nopanels -wpos 255,255 -%s 1",sokey);
            Vnm_systemBack(sysc);
            sprintf(sysc,"geomview -nopanels -wpos 255,255 -%s 2",sokey);
            Vnm_systemBack(sysc);
            sprintf(sysc,"geomview -nopanels -wpos 255,255 -%s 3",sokey);
            Vnm_systemBack(sysc);
            if (!strcmp("inet",Vsh_getenv(thee,"OSKEY"))) {
                Vnm_systemBack("mcbridge -i2u 0 0");
                Vnm_systemBack("mcbridge -i2u 1 1");
                Vnm_systemBack("mcbridge -i2u 2 2");
                Vnm_systemBack("mcbridge -i2u 3 3");
            }
        } else if (!strcmp("1x1s",Vsh_getenv(thee,"GVLO"))) {
            sprintf(sysc,"geomview -nopanels -wpos 282,245 -%s 0",sokey);
            Vnm_systemBack(sysc);
            if (!strcmp("inet",Vsh_getenv(thee,"OSKEY"))) {
                Vnm_systemBack("mcbridge -i2u 0 0");
            }
        } else if (!strcmp("2x1s",Vsh_getenv(thee,"GVLO"))) {
            sprintf(sysc,"geomview -nopanels -wpos 282,245 -%s 0",sokey);
            Vnm_systemBack(sysc);
            sprintf(sysc,"geomview -nopanels -wpos 282,245 -%s 1",sokey);
            Vnm_systemBack(sysc);
            if (!strcmp("inet",Vsh_getenv(thee,"OSKEY"))) {
                Vnm_systemBack("mcbridge -i2u 0 0");
                Vnm_systemBack("mcbridge -i2u 1 1");
            }
        } else if (!strcmp("3x1s",Vsh_getenv(thee,"GVLO"))) {
            sprintf(sysc,"geomview -nopanels -wpos 282,152 -%s 0",sokey);
            Vnm_systemBack(sysc);
            sprintf(sysc,"geomview -nopanels -wpos 282,152 -%s 1",sokey);
            Vnm_systemBack(sysc);
            sprintf(sysc,"geomview -nopanels -wpos 282,152 -%s 2",sokey);
            Vnm_systemBack(sysc);
            if (!strcmp("inet",Vsh_getenv(thee,"OSKEY"))) {
                Vnm_systemBack("mcbridge -i2u 0 0");
                Vnm_systemBack("mcbridge -i2u 1 1");
                Vnm_systemBack("mcbridge -i2u 2 2");
            }
        } else if (!strcmp("4x1s",Vsh_getenv(thee,"GVLO"))) {
            sprintf(sysc,"geomview -nopanels -wpos 255,255 -%s 0",sokey);
            Vnm_systemBack(sysc);
            sprintf(sysc,"geomview -nopanels -wpos 255,255 -%s 1",sokey);
            Vnm_systemBack(sysc);
            sprintf(sysc,"geomview -nopanels -wpos 255,255 -%s 2",sokey);
            Vnm_systemBack(sysc);
            sprintf(sysc,"geomview -nopanels -wpos 255,255 -%s 3",sokey);
            Vnm_systemBack(sysc);
            if (!strcmp("inet",Vsh_getenv(thee,"OSKEY"))) {
                Vnm_systemBack("mcbridge -i2u 0 0");
                Vnm_systemBack("mcbridge -i2u 1 1");
                Vnm_systemBack("mcbridge -i2u 2 2");
                Vnm_systemBack("mcbridge -i2u 3 3");
            }
        } else {
            Vnm_print(2,"%s: %s: Incorrect argument <%s>\n",
                thee->PR,"GVLO",Vsh_getenv(thee,"GVLO"));
        }
        Vnm_system("sleep 3");
        break;

      case vshcom_sockm:
        Vsh_killSockets(thee);
        strncpy(sofmt,Vsh_getenv(thee,"OSFMT"),VMAX_ARGLEN);
        if (!strcmp("unix",Vsh_getenv(thee,"OSKEY"))) {
            strncpy(sokey,"Mcs",VMAX_ARGLEN);
        } else if (!strcmp("inet",Vsh_getenv(thee,"OSKEY"))) {
            strncpy(sokey,"Mci",VMAX_ARGLEN);
        } else {
            strncpy(sokey,"Mci",VMAX_ARGLEN);
        }
        if (!strcmp("1x1",Vsh_getenv(thee,"GVLO"))) {
            sprintf(sysc,"mcsg -%s -wpos 436,445@55,520 -%s 0",sofmt,sokey);
            Vnm_systemBack(sysc);
        } else if (!strcmp("2x1",Vsh_getenv(thee,"GVLO"))) {
            sprintf(sysc,"mcsg -%s -wpos 436,445@55,520 -%s 0",sofmt,sokey);
            Vnm_systemBack(sysc);
            sprintf(sysc,"mcsg -%s -wpos 436,445@55,50  -%s 1",sofmt,sokey);
            Vnm_systemBack(sysc);
        } else if (!strcmp("3x1",Vsh_getenv(thee,"GVLO"))) {
            sprintf(sysc,"mcsg -%s -wpos 436,287@55,684 -%s 0",sofmt,sokey);
            Vnm_systemBack(sysc);
            sprintf(sysc,"mcsg -%s -wpos 436,287@55,367 -%s 1",sofmt,sokey);
            Vnm_systemBack(sysc);
            sprintf(sysc,"mcsg -%s -wpos 436,287@55,50  -%s 2",sofmt,sokey);
            Vnm_systemBack(sysc);
        } else if (!strcmp("4x1",Vsh_getenv(thee,"GVLO"))) {
            sprintf(sysc,"mcsg -%s -wpos 255,255 -%s 0",sofmt,sokey);
            Vnm_systemBack(sysc);
            sprintf(sysc,"mcsg -%s -wpos 255,255 -%s 1",sofmt,sokey);
            Vnm_systemBack(sysc);
            sprintf(sysc,"mcsg -%s -wpos 255,255 -%s 2",sofmt,sokey);
            Vnm_systemBack(sysc);
            sprintf(sysc,"mcsg -%s -wpos 255,255 -%s 3",sofmt,sokey);
            Vnm_systemBack(sysc);
        } else if (!strcmp("1x1s",Vsh_getenv(thee,"GVLO"))) {
            sprintf(sysc,"mcsg -%s -wpos 282,245 -%s 0",sofmt,sokey);
            Vnm_systemBack(sysc);
        } else if (!strcmp("2x1s",Vsh_getenv(thee,"GVLO"))) {
            sprintf(sysc,"mcsg -%s -wpos 282,245 -%s 0",sofmt,sokey);
            Vnm_systemBack(sysc);
            sprintf(sysc,"mcsg -%s -wpos 282,245 -%s 1",sofmt,sokey);
            Vnm_systemBack(sysc);
        } else if (!strcmp("3x1s",Vsh_getenv(thee,"GVLO"))) {
            sprintf(sysc,"mcsg -%s -wpos 282,152 -%s 0",sofmt,sokey);
            Vnm_systemBack(sysc);
            sprintf(sysc,"mcsg -%s -wpos 282,152 -%s 1",sofmt,sokey);
            Vnm_systemBack(sysc);
            sprintf(sysc,"mcsg -%s -wpos 282,152 -%s 2",sofmt,sokey);
            Vnm_systemBack(sysc);
        } else if (!strcmp("4x1s",Vsh_getenv(thee,"GVLO"))) {
            sprintf(sysc,"mcsg -%s -wpos 500,300 -%s 0",sofmt,sokey);
            Vnm_systemBack(sysc);
            sprintf(sysc,"mcsg -%s -wpos 500,300 -%s 1",sofmt,sokey);
            Vnm_systemBack(sysc);
            sprintf(sysc,"mcsg -%s -wpos 500,300 -%s 2",sofmt,sokey);
            Vnm_systemBack(sysc);
            sprintf(sysc,"mcsg -%s -wpos 500,300 -%s 3",sofmt,sokey);
            Vnm_systemBack(sysc);
        } else if (!strcmp("2x1r",Vsh_getenv(thee,"GVLO"))) {
            sprintf(sysc,"mcsg -%s -wpos 425,345@55,0 -%s 0",sofmt,sokey);
            Vnm_systemBack(sysc);
            sprintf(sysc,"mcsg -%s -wpos 425,345@55,369 -%s 1",sofmt,sokey);
            Vnm_systemBack(sysc);
        } else {
            Vnm_print(2,"%s: %s: Incorrect argument <%s>\n",
                thee->PR,"GVLO",Vsh_getenv(thee,"GVLO"));
        }
        Vnm_system("sleep 3");
        break;

      case vshcom_sockk:
        Vsh_killSockets(thee);
        break;

      case vshcom_sorry:
        Vnm_system("play sorry.au");
        break;

      default:
        rc = 0;
        break;
    }
    return rc;
}

/*
 * ***************************************************************************
 * Routine:  Vsh_ioSetup
 *
 * Purpose:  Setup for an I/O command.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC Vio *Vsh_ioSetup(Vsh *thee, char *key)
{
    char iodev[VMAX_BUFSIZE], iofmt[VMAX_BUFSIZE];
    char iohost[VMAX_BUFSIZE], iofile[VMAX_BUFSIZE];
    Vio *sock;

    /* setup for a read */
    if (!strcmp("r",key)) {

        strncpy(iohost,Vsh_getenv(thee,"IHVAL"),VMAX_BUFSIZE);

        if (!strcmp("sdio",Vsh_getenv(thee,"ISKEY"))) {
            strncpy(iodev,"SDIO",VMAX_BUFSIZE);
            strncpy(iofile,"console",VMAX_BUFSIZE);
        } else if (!strcmp("file",Vsh_getenv(thee,"ISKEY"))) {
            strncpy(iodev,"FILE",VMAX_BUFSIZE);
            strncpy(iofile,Vsh_getenv(thee,"IFNAM"),VMAX_BUFSIZE);
        } else if (!strcmp("buff",Vsh_getenv(thee,"ISKEY"))) {
            strncpy(iodev,"BUFF",VMAX_BUFSIZE);
            strncpy(iofile,Vsh_getenv(thee,"ISNAM"),VMAX_BUFSIZE);
        } else if (!strcmp("unix",Vsh_getenv(thee,"ISKEY"))) {
            strncpy(iodev,"UNIX",VMAX_BUFSIZE);
            strncpy(iofile,Vsh_getenv(thee,"ISNAM"),VMAX_BUFSIZE);
        } else if (!strcmp("inet",Vsh_getenv(thee,"ISKEY"))) {
            strncpy(iodev,"INET",VMAX_BUFSIZE);
            strncpy(iofile,Vsh_getenv(thee,"ISNAM"),VMAX_BUFSIZE);
        } else {
            Vnm_print(2,"Vsh_ioSetup: Internal logic error.\n");
            VJMPERR1( 0 );
        }

        if (!strcmp("asc",Vsh_getenv(thee,"ISFMT"))) {
            strncpy(iofmt,"ASC", VMAX_BUFSIZE);
        } else if (!strcmp("xdr",Vsh_getenv(thee,"ISFMT"))) {
            strncpy(iofmt,"XDR", VMAX_BUFSIZE);
        } else {
            Vnm_print(2,"Vsh_ioSetup: Internal logic error.\n");
            VJMPERR1( 0 );
        }

    /* setup for a write */
    } else if (!strcmp("w",key)) {

        strncpy(iohost,Vsh_getenv(thee,"OHVAL"),VMAX_BUFSIZE);

        if (!strcmp("sdio",Vsh_getenv(thee,"OSKEY"))) {
            strncpy(iodev,"SDIO",VMAX_BUFSIZE);
            strncpy(iofile,"console",VMAX_BUFSIZE);
        } else if (!strcmp("file",Vsh_getenv(thee,"OSKEY"))) {
            strncpy(iodev,"FILE",VMAX_BUFSIZE);
            strncpy(iofile,Vsh_getenv(thee,"OFNAM"),VMAX_BUFSIZE);
        } else if (!strcmp("buff",Vsh_getenv(thee,"OSKEY"))) {
            strncpy(iodev,"BUFF",VMAX_BUFSIZE);
            strncpy(iofile,Vsh_getenv(thee,"OSNAM"),VMAX_BUFSIZE);
        } else if (!strcmp("unix",Vsh_getenv(thee,"OSKEY"))) {
            strncpy(iodev,"UNIX",VMAX_BUFSIZE);
            strncpy(iofile,Vsh_getenv(thee,"OSNAM"),VMAX_BUFSIZE);
        } else if (!strcmp("inet",Vsh_getenv(thee,"OSKEY"))) {
            strncpy(iodev,"INET",VMAX_BUFSIZE);
            strncpy(iofile,Vsh_getenv(thee,"OSNAM"),VMAX_BUFSIZE);
        } else {
            Vnm_print(2,"Vsh_ioSetup: Internal logic error.\n");
            VJMPERR1( 0 );
        }

        if (!strcmp("asc",Vsh_getenv(thee,"OSFMT"))) {
            strncpy(iofmt,"ASC", VMAX_BUFSIZE);
        } else if (!strcmp("xdr",Vsh_getenv(thee,"OSFMT"))) {
            strncpy(iofmt,"XDR", VMAX_BUFSIZE);
        } else {
            Vnm_print(2,"Vsh_ioSetup: Internal logic error.\n");
            VJMPERR1( 0 );
        }

    } else {
        Vnm_print(2,"Vsh_ioSetup: Internal logic error.\n");
        VJMPERR1( 0 );
    }

    /* create socket and associate the buffer */
    VJMPERR1( VNULL != (sock=Vio_socketOpen(key,iodev,iofmt,iohost,iofile)) );
    Vio_bufTake(sock, thee->buf, thee->bufsize);
    thee->bufsize = 0;
    thee->buf     = VNULL;

    /* return without error */
    return sock;

  VERROR1: 
    Vnm_print(2,"Vsh_ioSetup: bailing out.\n");
    return VNULL;
}

/*
 * ***************************************************************************
 * Routine:  Vsh_ioCleanup
 *
 * Purpose:  Cleanup an I/O command.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */
VPUBLIC void Vsh_ioCleanup(Vsh *thee, Vio **sock)
{
    VJMPERR1( VNULL != thee );
    VJMPERR1( VNULL != *sock );

    /* snag the buffer before destroying the socket */
    thee->bufsize = Vio_bufSize(*sock);
    thee->buf     = Vio_bufGive(*sock);

    /* return without error */
    Vio_socketClose( sock );
    return;

  VERROR1: 
    Vnm_print(2,"Vsh_ioCleanup: bailing out.\n");
    return; 
}

