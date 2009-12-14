dnl ##########################################################################
dnl MALOC = < Minimal Abstraction Layer for Object-oriented C >
dnl Copyright (C) 1994--2008 Michael Holst
dnl
dnl This library is free software; you can redistribute it and/or
dnl modify it under the terms of the GNU Lesser General Public
dnl License as published by the Free Software Foundation; either
dnl version 2.1 of the License, or (at your option) any later version.
dnl
dnl This library is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
dnl Lesser General Public License for more details.
dnl
dnl You should have received a copy of the GNU Lesser General Public
dnl License along with this library; if not, write to the Free Software
dnl Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
dnl ##########################################################################

dnl ##########################################################################
dnl File:    configure.ac
dnl
dnl Purpose: configure.ac script for MALOC (Minimal Abstraction Layer
dnl          for Object-oriented C).
dnl
dnl Notes:   You process this file with autoconf to produce a "configure"
dnl          shell script by using several tools in the autoconf suite in
dnl          succession:
dnl
dnl             rm -rf config.cache autom4te.cache
dnl             aclocal \
dnl             && automake --gnu --add-missing --copy \
dnl             && autoconf \
dnl             && autoheader \
dnl             && libtoolize --automake --copy --force
dnl             rm -rf config.cache autom4te.cache
dnl
dnl          The above generates the script "configure", the package
dnl          pre-header "config.h", and all required "Makefile.in"
dnl          pre-makefiles.  To build the package you type the following:
dnl
dnl             rm -rf config.cache autom4te.cache
dnl             ./configure
dnl             make
dnl             make install
dnl
dnl          Running the "configure" script generates the package header
dnl          "config.h" from the pre-header "config.h.in", and generates
dnl          all "Makefiles" from all pre-makefiles "Makefile.in".
dnl
dnl Details: Below are some diagrams twisted from the GNU autoconf info
dnl          pages, with some automake paths thrown in, which show the
dnl          relationships between the various autoconf utilities.  The
dnl          names followed by "*" indicate the execution of a program
dnl          or a script (possibly produced at an earlier stage).  The
dnl          names in brackets "[ ]" represent optional steps, such as
dnl          the use of automake and/or autoheader.
dnl
dnl
dnl          Preparing a software package for distribution using autoconf:
dnl
dnl
dnl             source files-->[autoscan*]-->[configure.scan]-->configure.ac
dnl
dnl             [Makefile.am] -+---------> [automake*] --------> Makefile.in
dnl                            |
dnl             configure.ac --.   .--------> autoconf* ---------> configure
dnl                            +---+
dnl             [aclocal.m4] --+   `---.
dnl             [acsite.m4] ---'       |
dnl                                    +--> [autoheader*] ---> [config.h.in]
dnl                                    |
dnl             - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
dnl                                    |
dnl             [acconfig.h] ----.     |  BELOW THIS LINE IS NOW DEPRECATED;
dnl                              +-----'  AUTOCONF NOW SUPPORTS A COMPLETELY
dnl             [config.h.top] --+        SELF-CONTAINED configure.ac FILE.
dnl             [config.h.bot] --'
dnl
dnl
dnl
dnl          Building the resulting package using "configure" and "make":
dnl
dnl
dnl                                   .------------> config.cache
dnl             configure* -----------+------------> config.log
dnl                                   |
dnl             [config.h.in] -.      v           .-> [config.h] -.
dnl                            +-->config.status*-+               +--> make*
dnl             Makefile.in ---'                  `-> Makefile ---'
dnl
dnl Author:  Michael Holst
dnl ##########################################################################

dnl ##########################################################################
dnl 0. autoconf initialization
dnl ##########################################################################

dnl # Autoconf init, autoheader output location, config and prefix directories
AC_INIT(maloc, 1.0, mholst@math.ucsd.edu)
AC_CONFIG_SRCDIR(src/aaa_inc/maloc/maloc.h)
AC_CONFIG_HEADERS(src/aaa_inc/maloccf.h)
AC_CONFIG_AUX_DIR(config)
AC_PREFIX_DEFAULT(${HOME})

dnl # Get information we need later for the CPU-VENDOR-OS triplet
AC_CANONICAL_BUILD
AC_CANONICAL_HOST
AC_CANONICAL_TARGET

dnl ##########################################################################
dnl 1. automake initialization
dnl ##########################################################################

dnl -------------------------------------------
dnl AM_INIT_AUTOMAKE does:    PACKAGE="fetk";
dnl                           VERSION="1.0";
dnl                           AC_SUBST(PACKAGE)
dnl                           AC_SUBST(VERSION)
dnl                           AC_ARG_PROGRAM
dnl                           AC_PROG_MAKE_SET
dnl                           AM_SANITY_CHECK
dnl                           AC_PROG_INSTALL
dnl                           AM_MISSING_PROG
dnl So, we do this all in one shot
dnl AM_INIT_AUTOMAKE(maloc, 1.0)
dnl -------------------------------------------
dnl The above is deprecated; now we only do:
AM_INIT_AUTOMAKE([foreign]) 

dnl We need a separate FETK_VERSION macro for use by libtool
FETK_VERSION="1:0:0"
AC_SUBST(FETK_VERSION)

AM_MAINTAINER_MODE

dnl ##########################################################################
dnl 2. checks for programs
dnl ##########################################################################

dnl # Basic compiler and related tools
AC_LANG_C
AC_PROG_CC
AC_PROG_CPP
AM_PROG_LEX
AC_PROG_YACC

dnl # Libtool stuff
AC_DISABLE_SHARED
AM_PROG_LIBTOOL

dnl ##########################################################################
dnl 3. checks for libraries
dnl ##########################################################################

dnl # Check for the math library.
AC_CHECK_LIB(m,pow,[math_lib="-lm"],[math_lib=""])
AC_SUBST(math_lib)

dnl # Check for the Liberty library (useful for debugging on Linux).
AC_CHECK_LIB(iberty,strerrno,[liberty_lib="-liberty"],[liberty_lib=""])
AC_SUBST(liberty_lib)

dnl # Check for "gethostbyname" and nsl library
nsl_lib="";
AC_CHECK_FUNC(gethostbyname,[try_nsl=""],[try_nsl="yes"])
if test -n "${try_nsl}"; then
    AC_CHECK_LIB(nsl,gethostbyname,[nsl_lib="-lnsl"],[nsl_lib=""])
fi
AC_SUBST(nsl_lib)

dnl # Check for "connect" and socket library
socket_lib="";
AC_CHECK_FUNC(connect,[try_socket=""],[try_socket="yes"])
if test -n "${try_socket}"; then
    AC_CHECK_LIB(socket,connect,[socket_lib="-lsocket"],[socket_lib=""])
fi
AC_SUBST(socket_lib)

dnl # Check for "thr_create" and thread library
thread_lib="";
AC_CHECK_FUNC(thr_create,[try_thread=""],[try_thread="yes"])
if test -n "${try_thread}"; then
    AC_CHECK_LIB(thread,thr_create,[thread_lib="-lthread"],[thread_lib=""])
fi
AC_SUBST(thread_lib)

dnl ##########################################################################
dnl 4. checks for header files
dnl ##########################################################################

dnl # Check for the 15 Standard C (ANSI-C/ISO-C) headers.
AC_HEADER_STDC
dnl AC_CHECK_HEADER(assert.h,, AC_MSG_ERROR([STDC <assert.h> not found.]))
dnl AC_CHECK_HEADER(ctype.h,,  AC_MSG_ERROR([STDC <ctype.h> not found.]))
dnl AC_CHECK_HEADER(errno.h,,  AC_MSG_ERROR([STDC <errno.h> not found.]))
dnl AC_CHECK_HEADER(float.h,,  AC_MSG_ERROR([STDC <float.h> not found.]))
dnl AC_CHECK_HEADER(limits.h,, AC_MSG_ERROR([STDC <limits.h> not found.]))
dnl AC_CHECK_HEADER(locale.h,, AC_MSG_ERROR([STDC <locale.h> not found.]))
dnl AC_CHECK_HEADER(math.h,,   AC_MSG_ERROR([STDC <math.h> not found.]))
dnl AC_CHECK_HEADER(setjmp.h,, AC_MSG_ERROR([STDC <setjmp.h> not found.]))
dnl AC_CHECK_HEADER(signal.h,, AC_MSG_ERROR([STDC <signal.h> not found.]))
dnl AC_CHECK_HEADER(stdarg.h,, AC_MSG_ERROR([STDC <stdarg.h> not found.]))
dnl AC_CHECK_HEADER(stddef.h,, AC_MSG_ERROR([STDC <stddef.h> not found.]))
dnl AC_CHECK_HEADER(stdio.h,,  AC_MSG_ERROR([STDC <stdio.h> not found.]))
dnl AC_CHECK_HEADER(stdlib.h,, AC_MSG_ERROR([STDC <stdlib.h> not found.]))
dnl AC_CHECK_HEADER(string.h,, AC_MSG_ERROR([STDC <string.h> not found.]))
dnl AC_CHECK_HEADER(time.h,,   AC_MSG_ERROR([STDC <time.h> not found.]))

dnl # Check for some additional UNIX-like headers that we need.
AC_HEADER_SYS_WAIT

dnl ##########################################################################
dnl 5. checks for typedefs
dnl ##########################################################################

AC_TYPE_MODE_T
AC_TYPE_PID_T
AC_TYPE_SIZE_T

dnl ##########################################################################
dnl 6. checks for structures
dnl ##########################################################################

dnl ##########################################################################
dnl 7. checks for compiler characteristics
dnl ##########################################################################

dnl ##########################################################################
dnl 8. checks for library functions
dnl ##########################################################################

AC_HEADER_STAT

AC_MSG_CHECKING([for xdrmem_create in <rpc/rpc.h>])
AC_HEADER_EGREP(xdrmem_create, rpc/rpc.h,
    [xdrmem_create_ok=yes;AC_DEFINE(HAVE_XDR)],[xdrmem_create_ok=no;])
AC_MSG_RESULT(${xdrmem_create_ok})

AC_MSG_CHECKING([for getcwd in <unistd.h>])
AC_HEADER_EGREP(getcwd, unistd.h,
    [getcwd_ok=yes;AC_DEFINE(HAVE_GETCWD)],[getcwd_ok=no;])
AC_MSG_RESULT(${getcwd_ok})

AC_MSG_CHECKING([for O_NONBLOCK in <fcntl.h>])
AC_TRY_COMPILE(
    [#include <fcntl.h>],
    [unsigned short a; a && O_NONBLOCK;],
    [o_nonblock_ok=yes;AC_DEFINE(HAVE_O_NONBLOCK)],[o_nonblock_ok=no])
AC_MSG_RESULT(${o_nonblock_ok})

AC_MSG_CHECKING([whether accept() in <sys/socket.h> wants ulong, uint, or int])
dnl AC_DEFINE(ACCEPT_USES_ULONG)
dnl AC_DEFINE(ACCEPT_USES_UINT)
AC_MSG_RESULT([unclear])

dnl ##########################################################################
dnl 9. checks for system services
dnl ##########################################################################

dnl ##########################################################################
dnl 10. checks for anything else
dnl ##########################################################################

dnl # ------------------------------------------------
dnl # CHECK WHETHER TO DISABLE CPU-VENDOR-HOST TRIPLET
dnl # ------------------------------------------------

fetk_cpu_vendor_os="";
AC_SUBST(fetk_cpu_vendor_os)

dnl fetk_cpu_vendor_os="";
dnl AC_MSG_CHECKING([whether you want to disable the CPU-VENDOR-OS install triplet])
dnl AC_ARG_ENABLE(triplet,
dnl     [  --disable-triplet       disable CPU-VENDOR-OS install triplet [[default=no]]],
dnl     [ triplet_kill="yes"],
dnl     [ triplet_kill="" ]
dnl )
dnl if test "${triplet_kill}" == yes; then
dnl     AC_MSG_RESULT([yes])
dnl     fetk_cpu_vendor_os="";
dnl else
dnl     AC_MSG_RESULT([no])
dnl     fetk_cpu_vendor_os="${host_cpu}-${host_vendor}-${host_os}";
dnl fi
dnl AC_SUBST(fetk_cpu_vendor_os)

dnl # ---------------------------------------------------
dnl # HANDLE VARIOUS ARCHITECTURE-SPECIFIC CONFIGURATIONS
dnl # ---------------------------------------------------
if [[ "${host_os}" = "cygwin" ] && [ "${CC}" != "icl" ]]; then
    AC_DEFINE(MSDOS, [], [MSDOS Macro for Cygwin and gcc])
    AC_DEFINE(WIN32, [], [WIN32 Macro for Cygwin and gcc])
elif [[ "${host_os}" = "cygwin" ] && [ "${CC}" = "icl" ]]; then
   AC_DEFINE(MSDOS, [], [MSDOS Macro for Cygwin and icl - obsolete])
   AC_DEFINE(WIN32, [], [WIN32 Macro for Cygwin and icl])
#   AC_CHECK_HEADERS(\
#    unistd.h \
#    sys/types.h \
#    sys/time.h \
#    sys/times.h \
#    sys/stat.h \
#    fcntl.h \
#    rpc/rpc.h \
#    sys/socket.h \
#    sys/un.h \
#    netinet/in.h \
#    arpa/inet.h \
#    netdb.h \
#   )
else
   AC_CHECK_HEADERS(\
    unistd.h \
    sys/types.h \
    sys/time.h \
    sys/times.h \
    sys/stat.h \
    fcntl.h \
    rpc/rpc.h \
    sys/socket.h \
    sys/un.h \
    netinet/in.h \
    arpa/inet.h \
    netdb.h \ 
   )

dnl # ---------------------------
dnl # HANDLE THE READLINE LIBRARY
dnl # ---------------------------
rl_lib="";
rl_inc="";

dnl # READLINE header location ENV specification
AC_MSG_CHECKING([whether your environment defines FETK_RL_INCLUDE])
if test -n "${FETK_RL_INCLUDE}"; then
    AC_MSG_RESULT([yes])
    rl_inc_path="-I${FETK_RL_INCLUDE} ";
else
    AC_MSG_RESULT([no])
    rl_inc_path="";
fi

dnl # READLINE library location ENV specification
AC_MSG_CHECKING([whether your environment defines FETK_RL_LIBRARY])
if test -n "${FETK_RL_LIBRARY}"; then
    AC_MSG_RESULT([yes])
    rl_lib_path="-L${FETK_RL_LIBRARY} ";
else
    AC_MSG_RESULT([no])
    rl_lib_path="";
fi

dnl # Setup temp LDFLAGS and CPPFLAGS and look for READLINE library/header
LDFLAGS_SAVE=${LDFLAGS};
CPPFLAGS_SAVE=${CPPFLAGS};
LDFLAGS=${rl_lib_path};
CPPFLAGS=${rl_inc_path};
AC_CHECK_HEADER(readline/readline.h,
    [rl_use="yes";rl_inc="${rl_inc_path}"],
    [rl_use="";rl_inc=""])
if test -n "${rl_use}"; then
    AC_CHECK_HEADER(readline/history.h,
        [rl_use="yes";rl_inc="${rl_inc_path}"],
        [rl_use="";rl_inc=""])
    if test -n "${rl_use}"; then
        AC_CHECK_LIB(ncurses,tgetnum,
            [ncurses_lib="-lncurses"],
            [ncurses_lib=""],
            [${math_lib}])
        AC_CHECK_LIB(readline,readline,
            [rl_use="yes";rl_lib="${rl_lib_path}-lreadline ${ncurses_lib}"],
            [rl_use="";rl_lib=""],
            [${rl_lib_path} ${ncurses_lib} ${math_lib}])
        if test -n "${rl_use}"; then
            AC_CHECK_LIB(readline,add_history,
                [rl_use="yes";rl_lib="${rl_lib_path}-lreadline ${ncurses_lib}"],
                [rl_use="";rl_lib=""],
                [${rl_lib_path} ${ncurses_lib} ${math_lib}])
        fi
        if test -n "${rl_use}"; then
            AC_DEFINE(HAVE_READLINE_READLINE_H)
            AC_DEFINE(HAVE_READLINE_HISTORY_H)
        fi
    fi
fi

dnl # Return LDFLAGS and CPPFLAGS to normal
LDFLAGS=${LDFLAGS_SAVE};
CPPFLAGS=${CPPFLAGS_SAVE};

AC_SUBST(rl_lib)
AC_SUBST(rl_inc)

fi

dnl # ----------------------
dnl # HANDLE THE MPI LIBRARY
dnl # ----------------------
AC_ARG_WITH(mpi,
    [AC_HELP_STRING([--with-mpi], [Use MPI pre-configured from APBS])],
    [mpi_from_apbs="yes"],
    [mpi_from_apbs=""]
)

if test -n "${mpi_from_apbs}"; then

	echo "Using MPI enabled settings from APBS"
	AC_DEFINE([HAVE_MPI_H], [], [])
	mpi_inc=${APBS_MPI_INC}
	mpi_lib=${APBS_MPI_LIB}
	
else

	AC_MSG_CHECKING([whether you want MPI])
	AC_ARG_ENABLE(mpi,
		[  --enable-mpi            enable MPI layer [[default=no]]],
		[ mpi_use="yes" ],
		[ mpi_use=""    ]
	)
	mpi_lib="";
	mpi_inc="";
	if test -z "${mpi_use}"; then
		AC_MSG_RESULT([no])
	else
		AC_MSG_RESULT([yes ... poking around])

		dnl # MPI header location ENV specification
		AC_MSG_CHECKING([whether your environment defines FETK_MPI_INCLUDE])
		if test -n "${FETK_MPI_INCLUDE}"; then
			AC_MSG_RESULT([yes])
			mpi_inc_path="-I${FETK_MPI_INCLUDE} ";
		else
			AC_MSG_RESULT([no])
			mpi_inc_path="";
		fi

		dnl # MPI library location ENV specification
		AC_MSG_CHECKING([whether your environment defines FETK_MPI_LIBRARY])
		if test -n "${FETK_MPI_LIBRARY}"; then
			AC_MSG_RESULT([yes])
			mpi_lib_path="-L${FETK_MPI_LIBRARY} ";
		else
			AC_MSG_RESULT([no])
			mpi_lib_path="";
		fi

		dnl # Setup temp LDFLAGS and CPPFLAGS and look for MPI library/header
		LDFLAGS_SAVE=${LDFLAGS};
		CPPFLAGS_SAVE=${CPPFLAGS};
		LDFLAGS=${mpi_lib_path};
		CPPFLAGS=${mpi_inc_path};
		mpi_use="";
		AC_CHECK_HEADER(mpi.h,
			[mpi_use="yes";mpi_inc="${mpi_inc_path}"],
			[mpi_use="";mpi_inc=""])
		if test -n "${mpi_use}"; then
			AC_CHECK_LIB(mpi,main,
            [mpi_use="yes";mpi_lib="${mpi_lib_path}-lmpi"],
            [mpi_use="";mpi_lib=""],${mpi_lib_path})
			if test -z "${mpi_use}"; then
				AC_CHECK_LIB(mpich,main,
                [mpi_use="yes";mpi_lib="${mpi_lib_path}-lmpich"],
                [mpi_use="";mpi_lib=""],${mpi_lib_path})
			fi
			if test -n "${mpi_use}"; then
				AC_DEFINE([HAVE_MPI_H], [], [We have and want to use MPI])
			fi
		fi

		dnl # Return LDFLAGS and CPPFLAGS to normal
		LDFLAGS=${LDFLAGS_SAVE};
		CPPFLAGS=${CPPFLAGS_SAVE};
	fi
fi

AC_SUBST(mpi_lib)
AC_SUBST(mpi_inc)

dnl # -----------------------
dnl # HANDLE SYMBOL EMBEDDING
dnl # -----------------------
AC_MSG_CHECKING([whether EMBED macro can be used without compilation errors])
AC_TRY_COMPILE(,
    [#define EMBED(rctag) \
     static const char* rctag; \
     static void* use_rcsid=(0 ? &use_rcsid : (void*)&rcsid);
     EMBED(rcsid)],
    [embed_ok=yes;AC_DEFINE(HAVE_EMBED)],[embed_ok=no])
AC_MSG_RESULT(${embed_ok})

AC_MSG_RESULT([----------- begin processing MALOC configure options ---------])

dnl # ------------------------
dnl # HANDLE DEBUG COMPILATION
dnl # ------------------------
AC_MSG_CHECKING([whether you want debug compilation])
AC_ARG_ENABLE(vdebug,
    [  --enable-vdebug         enable debug compilation [[default=no]]],
    [ vdebug_use="yes" ],
    [ vdebug_use=""    ]
)
if test -z "${vdebug_use}"; then
    AC_MSG_RESULT([no])
else
    AC_MSG_RESULT([yes])
    AC_DEFINE(HAVE_DEBUG)
fi

dnl # --------------------------------
dnl # HANDLE ANSI/PEDANTIC COMPILATION
dnl # --------------------------------
AC_MSG_CHECKING([whether you want pedantic ANSI compilation])
AC_ARG_ENABLE(pansi,
    [  --enable-pansi          enable pedantic ANSI compilation [[default=no]]],
    [ pansi_use="yes" ],
    [ pansi_use=""    ]
)
if test -z "${pansi_use}"; then
    AC_MSG_RESULT([no])
    wall=""
    wall_pedantic=""
    wall_pedantic_ansi=""
    pedantic_ansi=""
else
    AC_MSG_RESULT([yes])
    wall="-Wall"
    wall_pedantic="-Wall -pedantic"
    wall_pedantic_ansi="-Wall -pedantic -ansi"
    pedantic_ansi="-pedantic -ansi"
fi
AC_SUBST(wall)
AC_SUBST(wall_pedantic)
AC_SUBST(wall_pedantic_ansi)
AC_SUBST(pedantic_ansi)

dnl # ----------------
dnl # HANDLE PROFILING
dnl # ----------------
AC_MSG_CHECKING([whether you want profiling])
AC_ARG_ENABLE(gprof,
    [  --enable-gprof          enable gprof profiling [[default=no]]],
    [ gprof_use="yes" ],
    [ gprof_use=""    ]
)
if test -z "${gprof_use}"; then
    AC_MSG_RESULT([no])
    profile=""
    profile_lib=""
else
    AC_MSG_RESULT([yes])
    profile="-pg"
    profile_lib="-pg"
fi
AC_SUBST(profile)
AC_SUBST(profile_lib)

dnl # ---------------------------------
dnl # HANDLE THE ELECTRIC FENCE LIBRARY
dnl # ---------------------------------
AC_MSG_CHECKING([whether you want to build ElectricFence])
AC_ARG_ENABLE(efence,
    [  --enable-efence         build ElectricFence ourself [[default=no]]],
    [ efence_build="yes" ],
    [ efence_build="" ]
)
if test -n "${efence_build}"; then
    AC_MSG_RESULT([yes])
    efence="efence"
else
    AC_MSG_RESULT([no])
    efence=""
fi
AC_SUBST(efence)
AM_CONDITIONAL(BUILD_EFENCE, test x$efence = xefence)

AC_MSG_RESULT([------------ end processing MALOC configure options ----------])

dnl ##########################################################################
dnl 11. templates
dnl ##########################################################################

AH_TEMPLATE([HAVE_XDR],
            [Do I have the XDR datastructure in the RPC package?])

AH_TEMPLATE([HAVE_O_NONBLOCK],
            [Do I have the O_NONBLOCK macro?])

AH_TEMPLATE([ACCEPT_USES_ULONG],
            [Does accept() use unsigned long?])

AH_TEMPLATE([ACCEPT_USES_UINT],
            [Does accept() use unsigned int?])

AH_TEMPLATE([HAVE_GETCWD],
            [Do I have the getcwd routine?])

AH_TEMPLATE([HAVE_CYGWIN],
            [Am I running in a Cygwin/Win32 environment?])

AH_TEMPLATE([HAVE_EMBED],
            [Does EMBED macro for embedding rcsid symbols into binaries work?])

AH_TEMPLATE([HAVE_DEBUG],
            [Do I compile as a debug version?])

AH_TEMPLATE([HAVE_MPI_H],
            [Do I have MPI header?])

AH_TEMPLATE([HAVE_READLINE_READLINE_H],
            [Do I have readline.h header?])

AH_TEMPLATE([HAVE_READLINE_HISTORY_H],
            [Do I have history.h header?])

AH_TOP([
/*
 * ***************************************************************************
 * MALOC = < Minimal Abstraction Layer for Object-oriented C >
 * Copyright (C) 1994--2000  Michael Holst
 * 
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 675 Mass Ave, Cambridge, MA 02139, USA.
 * 
 * rcsid="$Id: configure.ac,v 1.94 2006/07/14 06:23:15 mholst Exp $"
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * File:     acconfig.h
 *
 * Purpose:  Generates the main configuration header "maloccf.h" for MALOC.
 *
 * Notes:    See the comments at the top of the file "configure.ac" for
 *           an outline of the sequence of steps that turns acconfig.h
 *           into <src/aaa_inc/maloccf.h.in> and then eventually into
 *           <src/aaa_inc/maloccf.h> when you are using GNU autoconf.
 *
 *           This file can also form the basis for a manually-produced
 *           maloccf.h file.  In fact, a correct Win32 maloccf.h file can be
 *           generated simply by removing the two lines containing the
 *           GNU autoconf tags "TOP" and "BOTTOM".
 *
 *           The final autoconf (or manually) generated "maloccf.h" attempts
 *           to produce a correct header file layout for various UNIX-like
 *           and Win32 machines, giving access to several things beyond ISO
 *           C/C++, including BSD Signals, UNIX Domain sockets, INET TCP/IP
 *           sockets, and the WINSOCK implementation of INET TCP/IP sockets.
 *
 *           The MALOC library then provides a portable abstract interface
 *           to UNIX domain sockets, INET sockets, pipes, signals, and other
 *           system-dependent things that one usually wants to get to in a
 *           fairly standard C or C++ scientific software package.  Once
 *           MALOC is built, "maloccf.h" is no longer needed (it is not
 *           included in the set of API headers that are copied into the
 *           specified header install directory.  In other words, none of
 *           the MALOC headers forming the API include the config file
 *           "maloccf.h"; it is only included by the source files.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */

#ifndef _MALOCCF_H_
#define _MALOCCF_H_
])

AH_BOTTOM([
/*
 * ***************************************************************************
 * Handle some additional things manually (Win32, NeXT, etc)
 * ***************************************************************************
 */

/* Win32 configuration (non-CygWin) */
#if !defined(HAVE_CYGWIN)
#   if defined(WIN32) || defined(__WATCOMC__)

        /* Set the main key for specifying WIN32 code */
#       define HAVE_WIN32

        /* Deal with some basic problems with UNIX/WIN32 compatibility */
#       define HAVE_O_NONBLOCK 1
#       define HAVE_MODE_T 1
#       define HAVE_GETCWD 1

        /* WATCOM does STAT macros right; Microsoft does not */
#       if !defined(__WATCOMC__)
#           define STAT_MACROS_BROKEN 1
#       endif

        /* Required headers that exist in both UNIX and WIN32 */
#       define HAVE_SYS_TYPES_H 1
#       define HAVE_SYS_STAT_H 1
#       define HAVE_FCNTL_H 1
#       define HAVE_RPC_H 1

        /* Required headers that exist only in WIN32 (non-CygWin) */
#       define HAVE_DIRECT_H 1
#       define HAVE_PROCESS_H 1
#       define HAVE_WINSOCK_H 1
#       define HAVE_IO_H 1

#   endif
#endif

#if defined(NeXT) || defined(__NeXT__)
#   define HAVE_NEXT
#endif

/*
 * ***************************************************************************
 * Deal with macros we need that are sometimes missing
 * ***************************************************************************
 */

/* Deal with broken stat macros on some platforms */
#if !defined(STAT_MACROS_BROKEN)
#   define VS_ISREG(a) (((a) & S_IFMT) == S_IFREG)
#else
#   define VS_ISREG(a) (0)
#endif

/* Deal a missing macro on some unix platforms (NeXT, etc) */
#if !defined(HAVE_O_NONBLOCK)
#   define VO_NONBLOCK 00004
#else
#   define VO_NONBLOCK O_NONBLOCK
#endif

/*
 * ***************************************************************************
 * Define some RCS tag embedding and debug I/O macros
 * ***************************************************************************
 */

/* Embedded RCS tags ("ident filename" prints module versions in filename) */
#if defined(HAVE_EMBED)
#    define VEMBED(rctag) \
         VPRIVATE const char* rctag; \
         static void* use_rcsid=(0 ? &use_rcsid : (void*)&rcsid);
#else
#    define VEMBED(rctag)
#endif

/* Produce additional debugging I/O */
#if defined(HAVE_DEBUG)
#    define VDEBUGIO(str) fprintf(stderr,str)
#else
#    define VDEBUGIO(str)
#endif

#endif /* _MALOCCF_H_ */
])

dnl ##########################################################################
dnl 12. generate the makefiles
dnl ##########################################################################

AC_OUTPUT([
    Makefile
    src/Makefile

    src/base/Makefile

    src/efence/Makefile

    src/vsys/Makefile
    src/vsh/Makefile
    src/psh/Makefile

    src/aaa_inc/Makefile
    src/aaa_inc/maloc/Makefile

    src/aaa_lib/Makefile
])

