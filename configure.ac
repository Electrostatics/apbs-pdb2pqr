dnl ##########################################################################
dnl APBS -- Adaptive Poisson-Boltzmann Solver
dnl 
dnl Nathan A. Baker (baker@biochem.wustl.edu)
dnl Dept. of Biochemistry and Molecular Biophysics
dnl Center for Computational Biology
dnl Washington University in St. Louis
dnl 
dnl Additional contributing authors listed in the code documentation.
dnl 
dnl Copyright (c) 2006.  Washington University in St. Louis.
dnl All Rights Reserved.
dnl Portions Copyright (c) 1999-2002.  The Regents of the University of
dnl California.  
dnl Portions Copyright (c) 1995.  Michael Holst.
dnl 
dnl This file is part of APBS.
dnl 
dnl APBS is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 2 of the License, or
dnl (at your option) any later version.
dnl 
dnl APBS is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl 
dnl You should have received a copy of the GNU General Public License
dnl along with APBS; if not, write to the Free Software
dnl Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA
dnl
dnl Linking APBS statically or dynamically with other modules is making a
dnl combined work based on APBS. Thus, the terms and conditions of the GNU
dnl General Public License cover the whole combination.
dnl 
dnl SPECIAL GPL EXCEPTION
dnl In addition, as a special exception, the copyright holders of APBS
dnl give you permission to combine the APBS program with free software
dnl programs and libraries that are released under the GNU LGPL or with
dnl code included in releases of ISIM, Ion Simulator Interface, PMV, PyMOL
dnl SMOL, VMD, and Vision. Such combined software may be linked with APBS and 
dnl redistributed together in original or modified form as mere aggregation
dnl without requirement that the entire work be under the scope of the GNU 
dnl General Public License. This special exception permission is also extended
dnl to any software listed in the SPECIAL GPL EXCEPTION clauses by the PMG,
dnl FEtk, MC, or MALOC libraries.
dnl 
dnl Note that people who make modified versions of APBS are not obligated
dnl to grant this special exception for their modified versions; it is
dnl their choice whether to do so. The GNU General Public License gives
dnl permission to release a modified version without this exception; this
dnl exception also makes it possible to release a modified version which
dnl carries forward this exception.
dnl
dnl ##########################################################################

dnl ##########################################################################
dnl File: configure.in
dnl Purpose: configure.in script for APBS
dnl Note:  For use with autoconf-2.56 and above
dnl Author:  Michael Holst and Nathan Baker
dnl ##########################################################################

dnl ##########################################################################
dnl 0. autoconf initialization
dnl ##########################################################################
dnl # AC_INIT(src/aaa_inc/apbs/apbs.h)
AC_INIT([APBS], [0.4.0], [baker@biochem.wustl.edu])
AC_PREREQ(2.56)
AM_CONFIG_HEADER(src/aaa_inc/apbscfg.h)
AC_CONFIG_AUX_DIR(config)
AC_PREFIX_DEFAULT(${HOME})
AC_CANONICAL_SYSTEM
fetk_cpu_vendor_os="${host_cpu}-${host_vendor}-${host_os}";
AC_SUBST(fetk_cpu_vendor_os)

dnl ##########################################################################
dnl 1. automake initialization (uncomment only if rebuilding Makefile.in's)
dnl ##########################################################################
dnl # THIS IS THE VERSION INFORMATION FOR THE APBS LIBRARY.  ACCORDING TO 
dnl # THE LIBTOOL MANUAL:
dnl #   This flag accepts an argument of the form `current[:revision[:age]]'.
dnl #   So, passing `-version-info 3:12:1' sets current to 3, revision to 12,
dnl #   and age to 1.  If either revision or age are omitted, they default to 
dnl #   0. Also note that age must be less than or equal to the current 
dnl #   interface number.  Here are a set of rules to help you update your
dnl #   library version information:
dnl #      1. Start with version information of `0:0:0' for each libtool 
dnl #         library.
dnl #      2. Update the version information only immediately before a public
dnl #         release of your software. More frequent updates are unnecessary,
dnl #         and only guarantee that the current interface number gets larger
dnl #         faster.  
dnl #      3. If the library source code has changed at all since the last
dnl #         update, then increment revision (`c:r:a' becomes `c:r+1:a').
dnl #      4. If any interfaces have been added, removed, or changed since the
dnl #         last update, increment current, and set revision to 0.
dnl #      5. If any interfaces have been added since the last public release,
dnl #         then increment age.
dnl #      6. If any interfaces have been removed since the last public
dnl #         release, then set age to 0.
dnl #
dnl #   Never try to set the interface numbers so that they correspond to the
dnl #   release number of your package. This is an abuse that only fosters
dnl #   misunderstanding of the purpose of library versions.
AM_INIT_AUTOMAKE(apbs, 0.4.0)
APBS_VERSION="1:0:0"
AC_SUBST(APBS_VERSION)
ACLOCAL=":";
AUTOCONF=":";
AUTOMAKE=":";
AUTOHEADER=":";
AC_SUBST(ACLOCAL)
AC_SUBST(AUTOCONF)
AC_SUBST(AUTOMAKE)
AC_SUBST(AUTOHEADER)

dnl ##########################################################################
dnl 2. checks for programs
dnl ##########################################################################
AC_LANG_C
AC_PROG_CC
AM_PROG_CC_STDC
AC_PROG_CXX
AC_PROG_F77
AC_LANG_FORTRAN77
dnl AC_F77_LIBRARY_LDFLAGS 
dnl # -----------------------
dnl # HANDLE OSF WEIRDNESS:
dnl #   This includes the inability to link FORTRAN and C code compiled with
dnl #   debugger flags and a g77 that breaks on some FORTRAN code with the
dnl #   -O2 flag set.
dnl # -----------------------
short_os=`echo ${host_os} | sed 's/osf.*/osf/g'`;
if test -z "`echo ${short_os} | sed 's/osf//g'`"; then
  AC_MSG_RESULT([Disabling debugger compiler flags to avert OSF weirdness])
  if test -z ${FFLAGS-''}; then
      FFLAGS="";
  else
     FFLAGS=`echo ${FFLAGS} | sed 's/-g//g'`;
     if test -z "`echo ${F77-f77} | sed 's/g77//g'`"; then
         FFLAGS=`echo ${FFLAGS} | sed 's/-O2/-O/g'`;
     fi
  fi
  if test -z ${CFLAGS-''}; then
      CFLAGS="";
  else
     CFLAGS=`echo ${CFLAGS} | sed 's/-g//g'`;
  fi
  if test -z ${CPPFLAGS-''}; then
      CPPFLAGS="";
  else
     CPPFLAGS=`echo ${CPPFLAGS} | sed 's/-g//g'`;
  fi
  if test -z ${CXXFLAGS-''}; then
      CXXFLAGS="";
  else
     CXXFLAGS=`echo ${CXXFLAGS} | sed 's/-g//g'`;
  fi
fi
AC_F77_NOFORMAIN
FFLAGS="${FLIBS} ${FFLAGS}"
AC_F77_FCN_MANGLE()
AC_LANG_C
AC_PROG_CPP
AC_DISABLE_SHARED
AM_PROG_LIBTOOL

dnl ##########################################################################
dnl 3. checks for libraries
dnl ##########################################################################
dnl # Check for the math library.
AC_LANG_C
AC_CHECK_LIB(m,pow,[math_lib="-lm"],[math_lib=""])
AC_SUBST(math_lib)
dnl # Check for the liberty library
AC_LANG_C
AC_CHECK_LIB(iberty,strerrno,[liberty_lib="-liberty"],[liberty_lib=""])
AC_SUBST(liberty_lib)
dnl # Check for "gethostbyname" and nsl library
nsl_lib="";
AC_LANG_C
AC_CHECK_FUNC(gethostbyname,[try_nsl=""],[try_nsl="yes"])
if test -n "${try_nsl}"; then
    AC_CHECK_LIB(nsl,gethostbyname,[nsl_lib="-lnsl"],[nsl_lib=""])
fi
AC_SUBST(nsl_lib)
dnl # Check for "connect" and socket library
socket_lib="";
AC_LANG_C
AC_CHECK_FUNC(connect,[try_socket=""],[try_socket="yes"])
if test -n "${try_socket}"; then
    AC_CHECK_LIB(socket,connect,[socket_lib="-lsocket"],[socket_lib=""])
fi
AC_SUBST(socket_lib)
dnl # Check for "thr_create()" and thread library
thread_lib="";
AC_LANG_C
AC_CHECK_FUNC(thr_create,[try_thread=""],[try_thread="yes"])
if test -n "${try_thread}"; then
    AC_CHECK_LIB(thread,thr_create,[thread_lib="-lthread"],[thread_lib=""])
fi
AC_SUBST(thread_lib)
dnl # CHECK FOR TIME
try_time=""
AC_LANG_C
AC_CHECK_FUNC(time,[AC_DEFINE([HAVE_TIME_FUNC], [], [have the time() function])], [])
dnl # CHECK FOR RAND
try_rand=""
AC_LANG_C
AC_CHECK_FUNC(rand,[AC_DEFINE([HAVE_RAND_FUNC], [], [have the rand() function])], [])
dnl # CHECK FOR SRAND
try_srand=""
AC_LANG_C
AC_CHECK_FUNC(srand,[AC_DEFINE([HAVE_SRAND_FUNC], [], [have the srand() function])], [])

dnl ##########################################################################
dnl 4. checks for header files
dnl ##########################################################################

dnl ##########################################################################
dnl 5. checks for typedefs
dnl ##########################################################################

dnl ##########################################################################
dnl 6. checks for structures
dnl ##########################################################################

dnl ##########################################################################
dnl 7. checks for compiler characteristics
dnl ##########################################################################

dnl # ------------------------------------
dnl # HANDLE COMPILER FLOATING POINT ERROR
dnl # ------------------------------------

AC_FPERROR

dnl ##########################################################################
dnl 8. checks for library functions
dnl ##########################################################################

dnl ##########################################################################
dnl 9. checks for system services
dnl ##########################################################################

dnl ##########################################################################
dnl 10. checks for anything else
dnl ##########################################################################

AC_MSG_RESULT([---------- begin configuring for MALOC linkage ------------])

dnl # ----------------------
dnl # HANDLE THE MALOC LIBRARY
dnl # ----------------------

AC_MSG_CHECKING([for FETK_INCLUDE])
if test -n "${FETK_INCLUDE}"; then
    AC_MSG_RESULT([yes])
    maloc_inc_path="-I${FETK_INCLUDE} ";
else
    if test "${prefix}" != "NONE"; then
        AC_MSG_RESULT([no; will use ${prefix}/include])
        maloc_inc_path="-I${prefix}/include ";
    else
        AC_MSG_RESULT([no; will use ${HOME}/include])
        maloc_inc_path="-I${HOME}/include ";
    fi
fi

dnl # MALOC library location ENV specification
AC_MSG_CHECKING([for FETK_LIBRARY])
if test -n "${FETK_LIBRARY}"; then
    AC_MSG_RESULT([yes])
    maloc_lib_path="-L${FETK_LIBRARY}/${fetk_cpu_vendor_os} ";
else
    if test "${prefix}" != "NONE"; then
        AC_MSG_RESULT([no; will use ${prefix}/lib])
        maloc_lib_path="-L${prefix}/lib/${fetk_cpu_vendor_os} ";
    else
        AC_MSG_RESULT([no; will use ${HOME}/lib])
        maloc_lib_path="-L${HOME}/lib/${fetk_cpu_vendor_os} ";
    fi
fi

LDFLAGS_SAVE=${LDFLAGS};
CPPFLAGS_SAVE=${CPPFLAGS};
LDFLAGS=${maloc_lib_path};
CPPFLAGS=${maloc_inc_path};
maloc_found="";
AC_CHECK_HEADER(maloc/maloc.h,
    [maloc_found="yes";maloc_inc="${maloc_inc_path}"],
    [maloc_found="";maloc_inc=""])
if test -n "${maloc_found}"; then
    AC_CHECK_LIB(maloc,maloc_link,
        [maloc_found="yes";maloc_lib="${maloc_lib_path} -lmaloc"],
        [maloc_found="";maloc_lib=""],
        [${maloc_lib_path}-lmaloc -lvf2c ${mpi_lib} ${rl_lib} ${liberty_lib} ${socket_lib} ${nsl_lib} ${math_lib}])
fi
if test -z "${maloc_found}"; then
    AC_MSG_ERROR([MALOC libraries or headers missing:  exiting on error!])
fi
dnl # Return LDFLAGS and CPPFLAGS to normal
LDFLAGS=${LDFLAGS_SAVE};
CPPFLAGS=${CPPFLAGS_SAVE};
AC_SUBST(maloc_lib)
AC_SUBST(maloc_inc)

dnl # ----------------------
dnl # HANDLE THE MPI LIBRARY
dnl # ----------------------
LDFLAGS_SAVE=${LDFLAGS};
CPPFLAGS_SAVE=${CPPFLAGS};
LDFLAGS=${maloc_lib_path};
CPPFLAGS=${maloc_inc_path};
mpi_lib="";
AC_CHECK_LIB(maloc,maloc_needs_mpi,
    [mpi_use="yes"],
    [mpi_use=""],
    [${maloc_lib_path}-lmaloc -lvf2c ${liberty_lib} ${socket_lib} ${nsl_lib} ${math_lib}])
LDFLAGS=${LDFLAGS_SAVE};
CPPFLAGS=${CPPFLAGS_SAVE};
mpi_inc="";
if test -z "${mpi_use}"; then
    AC_MSG_RESULT([MALOC not compiled with MPI:  MPI disabled])
else

    dnl # MPI header location ENV specification
    AC_MSG_CHECKING([for FETK_MPI_INCLUDE environmental variable])
    if test -n "${FETK_MPI_INCLUDE}"; then
        mpi_inc_path_env="${FETK_MPI_INCLUDE} ";
        AC_MSG_RESULT([yes])
    else
        mpi_inc_path_env="";
        AC_MSG_RESULT([no])
    fi

    dnl # MPI library location ENV specification
    AC_MSG_CHECKING([for FETK_MPI_LIBRARY environmental variable])
    if test -n "${FETK_MPI_LIBRARY}"; then
        mpi_lib_path_env="${FETK_MPI_LIBRARY}";
        AC_MSG_RESULT([yes])
    else
        mpi_lib_path_env="";
        AC_MSG_RESULT([no])
    fi

    dnl # MPI header location ARG specification (override environment)
    AC_MSG_CHECKING([for --with-mpiinc option])
    AC_ARG_WITH(mpiinc,
        [AC_HELP_STRING([--with-mpiinc=PATH], [directory containing <mpi.h>])],
        [mpi_inc_path="-I${with_mpiinc} "; AC_MSG_RESULT([yes])],
        [mpi_inc_path="-I${mpi_inc_path_env}"; AC_MSG_RESULT([no])]
    )

    dnl # MPI library location ARG specification (override environment)
    AC_MSG_CHECKING([for --with-mpilib option])
    AC_ARG_WITH(mpilib,
        [AC_HELP_STRING([--with-mpilib=PATH], 
                        [directory containing <libmpi.a>])],
        [mpi_lib_path="${with_mpilib}"; AC_MSG_RESULT([yes])],
        [mpi_lib_path="${mpi_lib_path_env}"; AC_MSG_RESULT([no])]
    )
    dnl # Setup temp LDFLAGS and look for MPI library/header
    LDFLAGS_SAVE=${LDFLAGS}
    CPPFLAGS_SAVE=${CPPFLAGS}
    LDFLAGS="-L${mpi_lib_path} ${nsl_lib} ${socket_lib} -L${mpi_lib_path}";
    CPPFLAGS="${mpi_inc_path}";
    mpi_use="";
    AC_MSG_RESULT([Looking for MPI header with ${mpi_inc_path}...])
    AC_CHECK_HEADER(mpi.h,
        [mpi_use="yes";mpi_inc="${mpi_inc_path}"],
        [mpi_use="";mpi_inc=""])
    if test -n "${mpi_use}"; then
        AC_MSG_RESULT([Looking for MPI library with ${mpi_lib_path}...])
        
        AC_CHECK_FILE(${mpi_lib_path}/libmpich.a, 
           [HAS_MPICH="mpich"], [])
        AC_CHECK_FILE(${mpi_lib_path}/libmpich.so, 
           [HAS_MPICH="mpich"], [])   

        if test -n "${HAS_MPICH}"; then
          AC_CHECK_LIB(mpich,main,
          [mpi_use="yes";mpi_lib="-L${mpi_lib_path} -lmpich"],
          [mpi_use="";mpi_lib=""])
         if test -z "${mpi_use}"; then   
          AC_CHECK_LIB(mpi,main,
          [mpi_use="yes";mpi_lib="-L${mpi_lib_path} -lmpi";mpi_bin="pnoshell"],
          [mpi_use="";mpi_lib="";mpi_bin=""])
         fi

        else
          AC_CHECK_LIB(mpi,main,
          [mpi_use="yes";mpi_lib="-L${mpi_lib_path} -lmpi";mpi_bin="pnoshell"],
          [mpi_use="";mpi_lib="";mpi_bin=""])
         if test -z "${mpi_use}"; then
          AC_CHECK_LIB(mpich,main,
          [mpi_use="yes";mpi_lib="-L${mpi_lib_path} -lmpich"],
          [mpi_use="";mpi_lib=""])
         fi
        fi
    fi
    if test -z "${mpi_use}"; then
        AC_MSG_RESULT([MALOC was compiled with MPI enabled])
        AC_MSG_RESULT([MALOC needs to link against the MPI library])
        AC_MSG_RESULT([Try defining FETK_MPI_INCLUDE and FETK_MPI_LIBRARY])
        AC_ERROR([MPI library or header is missing!!])

    else
        AC_DEFINE([HAVE_MPI_H], [], [have the MPI library])
    fi

    dnl # Return LDFLAGS and CPPFLAGS to normal
fi
LDFLAGS=${LDFLAGS_SAVE};
CPPFLAGS=${CPPFLAGS_SAVE};
AC_SUBST(mpi_lib)
AC_SUBST(mpi_inc)
AC_SUBST(mpi_bin)

dnl # ---------------------------
dnl # HANDLE THE READLINE LIBRARY
dnl # ---------------------------

LDFLAGS_SAVE=${LDFLAGS};
CPPFLAGS_SAVE=${CPPFLAGS};
LDFLAGS=${maloc_lib_path};
CPPFLAGS=${maloc_inc_path};
rl_lib="";
AC_CHECK_LIB(maloc,maloc_needs_rl,
    [rl_use="yes"],
    [rl_use=""],
    [${maloc_lib_path}-lmaloc -lvf2c ${liberty_lib} ${socket_lib} ${nsl_lib} ${math_lib}])
if test -n "${rl_use}"; then

    dnl # READLINE library location ENV specification
    AC_MSG_CHECKING([whether your environment defines FETK_RL_LIBRARY])
    if test -n "${FETK_RL_LIBRARY}"; then
        AC_MSG_RESULT([yes])
        rl_lib_path="-L${FETK_RL_LIBRARY} ";
    else
        AC_MSG_RESULT([no])
        rl_lib_path="";
    fi

    dnl # Setup temp LDFLAGS and look for READLINE library
  
    LDFLAGS=${rl_lib_path};
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

    dnl # Return LDFLAGS to normal
    LDFLAGS=${LDFLAGS2_SAVE};

    dnl # Final decision
    if test -z "${rl_use}"; then
        AC_MSG_RESULT([MALOC was compiled with READLINE enabled])
        AC_MSG_RESULT([MALOC needs to link against the READLINE library])
        AC_MSG_ERROR([READLINE cannot be found!])
    fi
fi
LDFLAGS=${LDFLAGS_SAVE};
CPPFLAGS=${CPPFLAGS_SAVE};
AC_SUBST(rl_lib)

AC_MSG_RESULT([---------- end configuring for MALOC linkage --------------])

AC_MSG_RESULT([---------- begin configuring for Python wrappers --------------])

dnl # Check for Python (for wrappers)
AC_PATH_PROG(py_path, python)
if test -n "${py_path}"; then
   changequote(<<, >>)dnl
   PY_VERSION=`$py_path -c 'import sys; print sys.version[:3]'`   
   changequote([, ])dnl
   PY_PREFIX=`$py_path -c 'import sys; print sys.prefix'`
   PY_EXEC_PREFIX=`$py_path -c 'import sys; print sys.exec_prefix'`
   if test -f $PY_PREFIX/include/python$PY_VERSION/Python.h; then
      PY_CFLAGS="-I$PY_PREFIX/include/python$PY_VERSION -I$PY_EXEC_PREFIX/include/python$PY_VERSION"
   else
	  PY_CFLAGS=""
   fi
else
   PY_CFLAGS=""
fi

PY_LINKFLAGS=${CFLAGS}
PY_CC=${CC}

dnl # Platform specific linking flags
PY_LDFLAGS=""
PY_SHARED=""
case $host_os in
   linux*) 
      PY_SHARED="-shared";;
   irix*) 
      PY_SHARED="-shared"
      PY_LDFLAGS="-lftn -lm";;
   sun|solaris**) 
      PY_SHARED="-G"
      PY_LDFLAGS="-lF77 -lm -lM77 -lsunmath"
      PY_CFLAGS="";;
   cygwin|mingw32|windows**)
      PY_CFLAGS="";;
   darwin**)
      PY_CC=${F77}
      PY_LINKFLAGS="-Wl,-framework -Wl,Python -bundle";;
   *) 
      PY_CFLAGS="";;
esac

dnl # The Python wrappers use different flags
dnl # On different architectures

IS_STATIC=`echo ${LDFLAGS} | grep static`
case $host_cpu in
   x86_64*)
      if test "${F77}" = "pgf77"; then
	  PY_CC=${F77}
      else
	  PY_CFLAGS=""
      fi
      if test -n "${IS_STATIC}"; then
	  PY_CFLAGS=""
      fi;;
   powerpc64*)
      PY_CFLAGS="";;
   i*86)
      if test "${CC}" = "icc"; then
          AC_CHECK_LIB(ifcore_pic,d_int_val,[ifcore_lib="-lifcore_pic"],[ifcore_lib=""])
	      if test -n "${ifcore_lib}"; then
	          PY_LDFLAGS=${ifcore_lib}
	      else
	          PY_CFLAGS=""	
	      fi
	      if test -n "${IS_STATIC}"; then
	          PY_LINKFLAGS="$CFLAGS -static-libcxa -Wl,-static"
          fi
      else
          if test "$CC" != "gcc"; then
              PY_CFLAGS=""
          fi
      fi;; 
   ia64*)
      if test -n "${mpi_use}"; then
          PY_CFLAGS=""
      fi
      if test "${CC}" = "icc"; then
          AC_CHECK_LIB(ifcore_pic,d_int_val,[ifcore_lib="-lifcore_pic"],[ifcore_lib=""])
	      if test -n "${ifcore_lib}"; then
	          PY_LDFLAGS=${ifcore_lib}
	      else
	          PY_CFLAGS=""	
	      fi
	      if test -n "${IS_STATIC}"; then
	          if (test -f /usr/lib/libc.so) && (test -f /usr/lib/libm.so); then
	              PY_LINKFLAGS="$CFLAGS -static-libcxa -Wl,-static,/usr/lib/libc.so,/usr/lib/libm.so"
	          else
                      PY_CFLAGS=""
	          fi 
          fi
      else
          PY_CFLAGS=""
      fi 
esac

if test -n "${PY_CFLAGS}"; then
   AC_MSG_RESULT(Python wrappers enabled)
else
   AC_MSG_RESULT(Python wrappers disabled)
fi

AC_SUBST(PY_SHARED,[$PY_SHARED])
AC_SUBST(PY_LDFLAGS,[$PY_LDFLAGS])
AC_SUBST(PY_CFLAGS,[$PY_CFLAGS])
AC_SUBST(PY_LINKFLAGS,[$PY_LINKFLAGS])
AC_SUBST(PY_CC,[$PY_CC])
AM_CONDITIONAL(PY_PRESENT, [test -n "${PY_CFLAGS}"])
AC_MSG_RESULT([---------- end configuring for Python wrappers --------------])


dnl # BEGIN FETK LIBRARY/HEADER TEST BLOCK
dnl # -----------------------
dnl # HANDLE THE FEtk LIBRARIES
dnl #
dnl # We require MC to be present.  If this library is not found, then 
dnl # we simply disable compilation of the portions of APBS
dnl # which depend on FEtk.  Compilation can also be disabled by the
dnl # '--disable-fetk' flag.
dnl # -----------------------
dnl # MC header location ENV specification
AC_MSG_RESULT([---------- begin configuring for FEtk linkage --------------])
dnl # BEGIN FETK LIBRARY/HEADER TEST BLOCK
use_fetk=""
mc_lib=""
mc_inc=""
mc_do_src=""
mc_do_lib=""
mc_do_test=""
AC_MSG_CHECKING([for --disable-fetk])
AC_ARG_ENABLE(fetk,
    [AC_HELP_STRING([--disable-fetk],[disable FEtk features])],
    [ with_fetk="no" ],
    [ with_fetk=""    ]
)
if test -z "${with_fetk}"; then
    AC_MSG_RESULT([no])
    AC_MSG_CHECKING([for FETK_INCLUDE])
    if test -n "${FETK_INCLUDE}"; then
        mc_inc_path="-I${FETK_INCLUDE} ";
        AC_MSG_RESULT([yes])
    else
        AC_MSG_RESULT([no])
        mc_inc_path="-I${HOME}/include ";
    fi
    AC_MSG_CHECKING([for FETK_LIBRARY])
    if test -n "${FETK_LIBRARY}"; then
        AC_MSG_RESULT([yes])
        mc_lib_path="-L${FETK_LIBRARY}/${fetk_cpu_vendor_os} ";
    else
        AC_MSG_RESULT([no])
        mc_lib_path="-L${HOME}/lib/${fetk_cpu_vendor_os} ";
    fi
    dnl # Setup temp LDFLAGS and look for MC library/header
    LDFLAGS_SAVE=${LDFLAGS};
    CPPFLAGS_SAVE=${CPPFLAGS};
    LDFLAGS=${mc_lib_path};
    CPPFLAGS=${mc_inc_path};
    mc_found="";
    AC_CHECK_HEADER(mc/mc.h,
        [mc_found="yes";mc_inc="${mc_inc_path}"],
        [mc_found="";mc_inc=""])
    AC_CHECK_LIB(mc,mc_link,
        [mc_found="yes";mc_lib="${mc_lib_path}-lmc"],
        [mc_found="";mc_lib=""],
        [${mc_lib_path}-lmc ${maloc_lib} ${mpi_lib} ${rl_lib} ${liberty_lib} ${socket_lib} ${nsl_lib} ${math_lib}])
    dnl # Setup temp LDFLAGS and look for MCX library/header
    LDFLAGS=${mc_lib_path};
    CPPFLAGS=${mc_inc_path};
    mcx_found="";
    AC_CHECK_HEADER(mcx/mcx.h,
        [mcx_found="yes";mcx_inc="${mc_inc_path}"],
        [mcx_found="";mc_inc=""])
    AC_CHECK_LIB(mcx,mcx_link,
        [mcx_found="yes";mcx_lib="${mc_lib_path}-lmcx"],
        [mcx_found="";mcx_lib=""],
        [${mc_lib_path}-lmc ${maloc_lib} ${mpi_lib} ${rl_lib} ${liberty_lib} ${socket_lib} ${nsl_lib} ${math_lib}])
    

    dnl # Return LDFLAGS and CPPFLAGS to normal
    LDFLAGS=${LDFLAGS_SAVE}
    CPPFLAGS=${CPPFLAGS_SAVE}

    if test -n "${mc_found}"; then
        AC_MSG_RESULT([FEtk component MC found:  enabled.])
        AC_DEFINE([HAVE_MC_H], [], [have FEtk component MC])
    else
        AC_MSG_RESULT([FEtk component MC not found:  disabled.])
    fi
    if test -n "${mcx_found}"; then
        AC_MSG_RESULT([FEtk component MCX found:  enabled.])
        AC_DEFINE([HAVE_MCX_H], [], [have FEtk component MCX])
    else
        AC_MSG_RESULT([FEtk component MCX not found:  disabled.])
    fi
        

else
    AC_MSG_RESULT([yes])
    AC_MSG_RESULT([FEtk disabled at user request.])
fi
AC_SUBST(mcx_lib)
AC_SUBST(mcx_inc)
AC_SUBST(mc_lib)
AC_SUBST(mc_inc)

dnl # Check for various optional MCX functions 
AC_CHECK_LIB(mcx_vor,Vor_delaunay,
        [vor_found="yes"],
        [vor_found=""],
        [${mcx_lib_path}-lmcx ${mc_lib_path}-lmc ${maloc_lib} ${maloc_lib_path} -lblas -lvf2c ${mpi_lib} ${rl_lib} ${liberty_lib} ${socket_lib} ${nsl_lib} ${math_lib}])
if test -n "${vor_found}"; then
    vor_lib="${mcx_lib_path} -lmcx_vor"
    AC_DEFINE([HAVE_VOR], [], [have the Vor library])
fi
AC_CHECK_LIB(mcx_vmetis,Vmetis_ctor,
        [vmetis_found="yes"],
        [vmetis_found=""],
        [${mcx_lib_path}-lmcx ${mc_lib_path}-lmc ${maloc_lib} ${maloc_lib_path} -lblas -lvf2c ${mpi_lib} ${rl_lib} ${liberty_lib} ${socket_lib} ${nsl_lib} ${math_lib}])
if test -n "${vmetis_found}"; then
    vmetis_lib="${mcx_lib_path} -lmcx_vmetis"
    AC_DEFINE([HAVE_VMETIS], [], [have the Vmetis library])
fi

AC_MSG_RESULT([---------- end configuring for FEtk linkage --------------])

AC_MSG_RESULT([------------ begin processing APBS configure options ----------])

dnl # -----------------------
dnl # HANDLE SYMBOL EMBEDDING
dnl # -----------------------
AC_MSG_CHECKING([for EMBED macro can be used])
AC_TRY_COMPILE(,
    [#define EMBED(rctag) \
     static const char* rctag; \
     static void* use_rcsid=(0 ? &use_rcsid : (void**)&rcsid);
     EMBED(rcsid)],
    [embed_ok=yes; AC_DEFINE([HAVE_EMBED], [], [EMBED macro can be used])],
    [embed_ok=no])
AC_MSG_RESULT(${embed_ok})


dnl # --------------------------------
dnl # HANDLE ANSI/PEDANTIC COMPILATION
dnl # --------------------------------
AC_MSG_CHECKING([--enable-pansi])
AC_ARG_ENABLE(pansi,
    [AC_HELP_STRING([--enable-pansi], [enable pedantic ANSI compilation])],
    [pansi_use="yes"],
    [pansi_use=""   ]
)
if test -z "${pansi_use}"; then
    AC_MSG_RESULT([no])
    wall=""
    wall_pedantic=""
    wall_pedantic_ansi=""
    pedantic_ansi=""
    dnl # --------------------------
    dnl # We only look for the non-ANSI strcasecmp if configure wasn't invoked
    dnl # with --enable-pansi
    dnl # --------------------------
    try_strcasecmp="";
    AC_CHECK_FUNC(strcasecmp, [try_strcasecmp=""], [try_strcasecmp="yes"])
    if test -n "${try_strcasecmp}"; then
        AC_DEFINE([HAVE_STRCASECMP], [], [have machine-supplied strcasecmp])
    fi
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

dnl # -----------------------
dnl # HANDLE THE TOOLS DIRECTORY
dnl # Disable compilation of executables in the tools directory with the
dnl # '--disable-tools' flag (for Alpha).
dnl # -----------------------
AC_MSG_CHECKING([for --disable-tools])
AC_ARG_ENABLE(tools,
    [AC_HELP_STRING([--disable-tools],[disable tool compilation])],
    [ with_tools="no" ],
    [ with_tools=""    ]
)
apbs_do_tools=""
if test -z "${with_tools}"; then
    AC_MSG_RESULT([yes])
    apbs_do_tools="tools"
else
    AC_MSG_RESULT([no])
fi
AC_SUBST(apbs_do_tools)


dnl # -----------------------
dnl # HANDLE THE BLAS LIBRARIES
dnl # -----------------------
blas_lib=""
AC_ARG_WITH(blas,
    [AC_HELP_STRING([--with-blas=LDFLAGS],
                  [link to user-specified BLAS (specify flags to pass to ld)])],
    [ blas_lib="${withval}"],
    [ blas_lib="" ]
)
AC_MSG_CHECKING([for --with-blas])
if test -n "${blas_lib}"; then
    AC_MSG_RESULT([yes])

    AC_ARG_WITH(blas_name,
        [AC_HELP_STRING([--with-blas-name=NAME],
                  [name of BLAS library, defaults to blas])],
        [ blas_name="${withval}"],
        [ blas_name=blas]
    )
    AC_MSG_RESULT([Looking for BLAS using ${blas_lib}])
    AC_LANG_FORTRAN77
    LDFLAGS_SAVE=${LDFLAGS};
    CPPFLAGS_SAVE=${CPPFLAGS};
    LDFLAGS="-L${maloc_lib_path} ${blas_lib}";
    AC_CHECK_LIB(${blas_name},dscal,
        [blas_found="yes"],
        [blas_found=""],
        [${blas_lib} ${socket_lib} ${nsl_lib} ${math_lib}])
    LDFLAGS=${LDFLAGS_SAVE};
    CPPFLAGS=${CPPFLAGS_SAVE};
    if test -z "${blas_found}"; then
        AC_MSG_ERROR([BLAS libraries or headers missing:  exiting on error!])
    fi
    AC_LANG_C

else
    AC_MSG_RESULT([no])
    blas_lib_path="${maloc_lib_path}"
    blas_lib="${blas_lib_path} -lblas ${blas_lib_path} -lvf2c"
fi
AC_SUBST(blas_lib)

dnl # -----------------------
dnl # If using Portland compilers, do not compile
dnl # tools/mesh/uhbd_asc2bin.f
dnl # -----------------------
case $F77 in
   pgf77*)
      UHBD_PORTLAND=""
      FFLAGS="${FFLAGS} -Mnomain";;
   *)
      UHBD_PORTLAND="Y";;
esac
AM_CONDITIONAL(UHBD_COMPILE, [test -n "${UHBD_PORTLAND}"])


dnl # -----------------------
dnl # CHECK WHETHER OR NOT TO LINK IN ARPACK FUNCTIONS
dnl # -----------------------
AC_MSG_CHECKING([for --with-arpack])
AC_ARG_WITH(arpack,
    [AC_HELP_STRING([--with-arpack=PATH],
                    [enable ARPACK functions (specify path to lib)])],
    [ arpack_lib="-larpack -L${withval}"],
    [ arpack_lib="" ]
)
if test -n "${arpack_lib}"; then
    arpack_do_tools="arpack"
    AC_MSG_RESULT([yes])
else
    arpack_do_tools=""
    AC_MSG_RESULT([no])
fi
AC_SUBST(arpack_lib)
AC_SUBST(arpack_do_tools)

dnl # -----------------------
dnl # CHECK FOR FMM LIBRARY
dnl # -----------------------
AC_LANG_FORTRAN77
LDFLAGS_SAVE=${LDFLAGS};
AC_CHECK_LIB(tree,treepeforce,
        [tree_found="yes"],
        [tree_found=""],
        [${mcx_lib_path} ${liberty_lib} ${socket_lib} ${nsl_lib} ${math_lib}])
if test -n "${tree_found}"; then
    tree_lib="${mc_lib_path} -ltree"
    AC_DEFINE([HAVE_TREE], [], [have the FMM library])
else
    tree_lib=""
fi
AC_LANG_C
AC_SUBST(tree_lib)

dnl # -----------------------
dnl # CHECK WHETHER YOU WANT INLINE FUNCTIONS
dnl # -----------------------
AC_MSG_CHECKING([for --disable-inline])
AC_ARG_ENABLE(inline,
    [AC_HELP_STRING([--disable-inline],[disable inline functions (macros)])],
    [ use_inline="${enableval}"],
    [ use_inline="" ]
)
if test -z "${use_inline}"; then
    AC_MSG_RESULT([yes])
else
    if test "${use_inline}" == no; then
        AC_MSG_RESULT([no])
        AC_DEFINE([APBS_NOINLINE], [1], [do not inline functions])
    else
        AC_MSG_RESULT([yes])
    fi
fi

dnl # -----------------------
dnl # CHECK WHETHER YOU WANT TO INCLUDE TINKER
dnl # -----------------------
AC_MSG_CHECKING([for --with-tinker])
AC_ARG_WITH(tinker,
    [AC_HELP_STRING([--with-tinker],[include TINKER support])],
    [ use_tinker="yes"],
    [ use_tinker="" ]
)
if test -n "${use_tinker}"; then
    AC_MSG_RESULT([yes])
    AC_DEFINE([WITH_TINKER], [], [include TINKER support])
else
    AC_MSG_RESULT([no]) 
fi

dnl # -----------------------
dnl # CHECK FOR WINDOWS SPECIFIC LIBRARIES
dnl # -----------------------
case $host_os in
   mingw32*)
      AC_CHECK_LIB(wsock32,main,
        [wsock_found="yes"],
        [wsock_found=""])
      if test -n "${wsock_found}"; then
           socket_lib="-lwsock32"
      else
           AC_MSG_ERROR([Unable to find wsock32 libraries!])
      fi
      ;;
   cygwin*)
      AC_CHECK_LIB(rpclib,main,
        [rpc_found="yes"],
        [rpc_found=""])
      if test -n "${rpc_found}"; then
           socket_lib="-lrpclib"
      fi
      ;;
esac
AC_MSG_RESULT([------------- end processing APBS configure options -----------])

dnl # -----------------------
dnl # ASSEMBLE A LIST OF LIBRARIES AND INCLUDE FILES THAT WILL BE USED IN ALL
dnl # COMPILATIONS
dnl # -----------------------
top_srcdir="`pwd`"
lib_list="${FLIBS} -L${top_srcdir}/lib ${tree_lib} ${vmetis_lib} ${vor_lib} ${mcx_lib} ${mc_lib} ${maloc_lib} ${blas_lib} ${thread_lib} ${mpi_lib} ${math_lib} ${rl_lib} ${socket_lib} ${nsl_lib}"
inc_list="-I. ${maloc_inc} ${mc_inc} ${mcx_inc} ${mpi_inc}"
AC_SUBST(lib_list)
AC_SUBST(inc_list)

dnl ##########################################################################
dnl 11. generate the makefiles
dnl ##########################################################################

AC_OUTPUT([
    Makefile

    tools/Makefile
    tools/mesh/Makefile
    tools/manip/Makefile
    tools/python/Makefile
    tools/arpack/Makefile
    tools/python/vgrid/Makefile

    src/Makefile
    src/generic/Makefile
    src/mg/Makefile
    src/fem/Makefile
    src/blas/Makefile
    src/aaa_inc/Makefile
    src/aaa_inc/apbs/Makefile
    src/aaa_lib/Makefile

    examples/Makefile
    examples/actin-dimer/Makefile
    examples/born/Makefile
    examples/FKBP/Makefile
    examples/hca-bind/Makefile
    examples/ion-protein/Makefile
    examples/ionize/Makefile
    examples/pka-lig/Makefile
    examples/point-pmf/Makefile
    examples/solv/Makefile
    
    bin/Makefile


])

AC_MSG_RESULT([])
AC_MSG_RESULT([APBS to be linked with following libraries:  ${lib_list}])
AC_MSG_RESULT([APBS to be compiled with following includes:  ${inc_list}])
AC_MSG_RESULT([])
AC_MSG_RESULT([Please type "make; make install"])

