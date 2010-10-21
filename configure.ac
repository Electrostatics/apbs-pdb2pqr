1.3dnl ##########################################################################
dnl APBS -- Adaptive Poisson-Boltzmann Solver
dnl 
dnl Nathan A. Baker (nathan.baker@pnl.gov)
dnl Pacific Northwest National Laboratory
dnl 
dnl Additional contributing authors listed in the code documentation.
dnl 
dnl Copyright (c) 2010.  Pacific Northwest National Laboratory.
dnl Portions Copyright (c) 2002-2010.  Washington University in St. Louis.
dnl Portions Copyright (c) 1999-2002.  The Regents of the University of
dnl California.  
dnl Portions Copyright (c) 1995.  Michael Holst.
dnl All Rights Reserved.
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
dnl code included in releases of ISIM, Ion Simulator Interface, PMV, 
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
AC_INIT([APBS], [1.3], [baker@biochem.wustl.edu])
AC_PREREQ(2.56)
AM_CONFIG_HEADER(src/aaa_inc/apbscfg.h)
AC_CONFIG_AUX_DIR(config)
AC_PREFIX_DEFAULT(${HOME})
AC_CANONICAL_SYSTEM

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
AM_INIT_AUTOMAKE(apbs, 1.3)
APBS_VERSION="1:3"
AC_SUBST(APBS_VERSION)
ACLOCAL=":";
AUTOCONF=":";
AUTOMAKE=":";
AUTOHEADER=":";
AC_SUBST(ACLOCAL)
AC_SUBST(AUTOCONF)
AC_SUBST(AUTOMAKE)
AC_SUBST(AUTOHEADER)

dnl # ---------------------------
dnl # SET CONVENIENCE PATHS
dnl # ---------------------------
top_apbsdir=`pwd`
contrib_path=`pwd`/contrib
AC_SUBST(top_apbsdir)
AC_SUBST(contrib_path)

dnl # ---------------------------
dnl # HANDLE THE FAST OPTION
dnl # ---------------------------

AC_MSG_CHECKING([for --enable-fast option])
AC_ARG_ENABLE(
	[fast],
	AS_HELP_STRING([--enable-fast], [enables fast computation]),,
	[enable_fast=no]
)

if test x"$enable_fast" != "xno"; then
	AC_DEFINE([APBS_FAST], [], [using the fast APBS parameters])
	AC_MSG_RESULT([yes])
else
	AC_MSG_RESULT([no])
fi

dnl # -------------------------------------------------------------------------
dnl # Set the prefix value explicitly if it isn't set
dnl # -------------------------------------------------------------------------
if test [ "${prefix}" = "NONE" ] ; then
	prefix=$HOME
fi

dnl ##########################################################################
dnl 2. checks for programs
dnl ##########################################################################
AC_MSG_CHECKING([for --enable-python])
AC_ARG_ENABLE(
	[python],
	AS_HELP_STRING([--enable-python],[enables Python wrapper compilation]),,
	[enable_python=no]
)
if test x"$enable_python" != "xno"; then
	AC_MSG_RESULT([yes])
	AC_MSG_NOTICE([Python linkage requires GNU compilers])
	cclist="gcc"
	f77list="gfortran g95 g77"
	apbs_do_python="python"
else
	AC_MSG_RESULT([no])
	cclist="icc pgcc gcc"
	f77list="ifort pgf77 gfortran g95 g77"
	apbs_do_python=""
fi
AC_SUBST(apbs_do_python)

AC_LANG_C
AC_PROG_CC(${cclist})
AM_PROG_CC_STDC
AC_PROG_CXX
AC_PROG_F77(${f77list})

dnl # -------------------------------------------------------------------------
dnl # Test to see if the compiler combination is valid for gcc/g77 or
dnl # gcc/gfortran/g95. Since code compiled with gcc 4.0 is only link
dnl # compatible with gfortran/g95. If they are using g77 then they need to
dnl # use gcc < 4.0 - D. Gohara 8/30/2006
dnl # -------------------------------------------------------------------------
if test ["$CC" = "gcc"]; then
	gcc_version=`gcc --version`
	gcc_complete_version=`echo ${gcc_version} | awk '{printf("%s",$3)}'`
	gcc_major_version=`echo ${gcc_complete_version} | awk '{printf("%i",$1)}'`

	if test ["$F77" = "g77"]; then
		if test ["${gcc_major_version}" -ge 4]; then
			AC_MSG_ERROR([Error, your GCC version is: ${gcc_complete_version} which is not compatible])
			AC_MSG_ERROR([with g77. You need to either use gcc < 4.0 as your C/CC/CXX compiler])
			AC_MSG_ERROR([or gfortran/g95 as your Fortran compiler ... exiting.])
		fi
	fi
fi

dnl AC_F77_LIBRARY_LDFLAGS 
dnl # -----------------------
dnl # HANDLE OSF WEIRDNESS:
dnl #   This includes the inability to link FORTRAN and C code compiled with
dnl #   debugger flags and a g77 that breaks on some FORTRAN code with the
dnl #   -O2 flag set.
dnl # -----------------------
short_os=`echo ${host_os} | sed 's/osf.*/osf/g'`;
if test -z "`echo ${short_os} | sed 's/osf//g'`"; then
	AC_MSG_NOTICE([Disabling debugger compiler flags to avert OSF weirdness])
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
if test "${CC}" != "icl"; then
  AC_F77_FCN_MANGLE()
fi
AC_LANG_C
AC_PROG_CPP
AC_DISABLE_SHARED
AM_PROG_LIBTOOL
FFLAGS="${FFLAGS}"

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

fetk_cpu_vendor_os="";
AC_SUBST(fetk_cpu_vendor_os)

AC_MSG_NOTICE([-----------------------------------------------------------])
AC_MSG_NOTICE([---------- begin configuring for MALOC linkage ------------])
AC_MSG_NOTICE([-----------------------------------------------------------])

dnl # ----------------------
dnl # HANDLE THE MALOC AND FETK LIBRARIES
dnl # ----------------------

contrib_path=`pwd`/contrib

dnl First check to see if they are using FETK_INCLUDE or FETK_LIBRARY
dnl if they are abort. We no longer support this...
if test -n "${FETK_INCLUDE}"; then
	AC_MSG_WARN([Error... FETK_INCLUDE is no longer supported.])
	AC_MSG_WARN([Use --with-fetk-include=PATH instead.])
	AC_MSG_WARN([Ignoring these variables...])
fi

if test -n "${FETK_LIBRARY}"; then
	AC_MSG_WARN([Error... FETK_LIBRARY is no longer supported.])
	AC_MSG_WARN([Use --with-fetk-library=PATH instead.])
	AC_MSG_WARN([Ignoring these variables...])
fi

dnl Next see if the user specified fetk-include or fetk-lib if they did
dnl process those and use those values instead of our own MALOC
AC_MSG_CHECKING([for --with-fetk-include option])
AC_ARG_WITH(
	[fetk-include],
	AS_HELP_STRING([--with-fetk-include=PATH], [full path to the FETk INCLUDE directory]),
	[fetk_include_path="${with_fetk_include}"]
)
if test x"$fetk_include_path" != "xno" -a -n "$fetk_include_path"; then
	AC_MSG_RESULT([${fetk_include_path}])
	USE_FETK_INCLUDE="yes"
else
	AC_MSG_RESULT([no])
fi

AC_MSG_CHECKING([for --with-fetk-library option])
AC_ARG_WITH(fetk-library,
	AS_HELP_STRING([--with-fetk-library=PATH], [full path to the FETk LIBRARY directory]),
	[fetk_library_path="${with_fetk_library}"]
)
if test x"$fetk_library_path" != "xno" -a -n "$fetk_library_path"; then
	AC_MSG_RESULT([${fetk_library_path}])
	USE_FETK_LIBRARY="yes"
else
	AC_MSG_RESULT([no])
fi

USE_FETK=
if test -n "${USE_FETK_INCLUDE}"; then
	maloc_inc_path="${fetk_include_path}"
	AC_CHECK_FILE(${maloc_inc_path}/maloc/maloc.h,[maloc_inc="yes"],[maloc_inc=""])	
	if test -z "${maloc_inc}"; then
		AC_MSG_ERROR([[Error... Couldn't find FETK include files in the location you specified]])
	fi
else
	maloc_inc_path="${contrib_path}/include";
fi
AC_MSG_NOTICE([Using MALOC headers in ${maloc_inc_path}])

if test -n "${USE_FETK_LIBRARY}"; then
	maloc_lib_path="${fetk_library_path}"
	AC_CHECK_FILE(${maloc_lib_path}/libmaloc.a,[maloc_lib="yes"],[maloc_lib=""])
	if test -z "${maloc_lib}"; then
		AC_MSG_ERROR([Error... Couldn't find FETK library files in the location you specified])
	fi
	USE_FETK="yes"
else
	maloc_lib_path="${contrib_path}/lib";
fi
AC_MSG_NOTICE([Using MALOC libraries in ${maloc_lib_path}])

maloc_inc="-I${maloc_inc_path}"
maloc_lib="-L${maloc_lib_path} -lmaloc"

dnl # Return LDFLAGS and CPPFLAGS to normal
AC_SUBST(maloc_lib)
AC_SUBST(maloc_inc)

dnl # ---------------------------
dnl # HANDLE THE MPI LIBRARY
dnl # ---------------------------

dnl # MPI header location ARG specification (override environment)
AC_MSG_CHECKING([for --with-mpich option])
AC_ARG_WITH(
	[mpich],
	AS_HELP_STRING([--with-mpich=PATH], [toplevel MPICH directory]),
	[mpi_path="${with_mpich}"]
)

if test x"$mpi_path" != "xno" -a -n "$mpi_path"; then
	AC_MSG_RESULT([${mpi_path}])
	USE_MPICH="yes"
else
	AC_MSG_RESULT([no])
	mpi_path=
fi

if test -z "${mpi_path}"; then
	AC_MSG_CHECKING([for --with-mpich2 option])
	AC_ARG_WITH(
		[mpich2],
		AS_HELP_STRING([--with-mpich2=PATH], [toplevel MPICH2 directory]),
		[mpi_path="${with_mpich2}"]
	)
	if test x"$mpi_path" != "xno" -a -n "$mpi_path"; then
		AC_MSG_RESULT([${mpi_path}])
		USE_MPICH2="yes"
	else
		AC_MSG_RESULT([no])
		mpi_path=
	fi
fi

if test -z "${mpi_path}"; then
	AC_MSG_CHECKING([for --with-lam option])
	AC_ARG_WITH(
		[lam],
		AS_HELP_STRING([--with-lam=PATH], [toplevel LAM-MPI directory]),
		[mpi_path="${with_lam}"]
	)
	if test x"$mpi_path" != "xno" -a -n "$mpi_path"; then
		AC_MSG_RESULT([${mpi_path}])
		USE_LAM="yes"
	else
		AC_MSG_RESULT([no])
		mpi_path=
	fi
fi

if test -z "${mpi_path}"; then
	AC_MSG_CHECKING([for --with-openmpi option])
	AC_ARG_WITH(
		[openmpi],
		AS_HELP_STRING([--with-openmpi=PATH], [enable OpenMPI compilation]),
		[mpi_path="${with_openmpi}"],
	)
	if test x"$mpi_path" != "xno" -a -n "$mpi_path"; then
		AC_MSG_RESULT([${mpi_path}])
		USE_OPMI="yes"
	else
		AC_MSG_RESULT([no])
		mpi_path=
	fi
fi

dnl # MPI library location ARG specification (override environment)
USE_MPI=
apbs_do_mpi=""
if test -n "${mpi_path}"; then
	AC_CHECK_FILE(${mpi_path}/include/mpi.h,[HAS_MPI_INC="yes"], [])
	if test -n "${HAS_MPI_INC}"; then
		mpi_inc="-I${mpi_path}/include"
	else
		AC_MSG_ERROR([Error... mpi.h missing in ${mpi_path}/include], 0)
	fi
	
	if test -n "${USE_MPICH}"; then
		AC_CHECK_FILE(${mpi_path}/lib/libmpich.a, [HAS_MPI_LIB="yes"], [])
		AC_CHECK_FILE(${mpi_path}/lib/libmpich.so, [HAS_MPI_LIB="yes"], [])
		AC_CHECK_FILE(${mpi_path}/lib/libmpich.dylib, [HAS_MPI_LIB="yes"], [])
		
		if test -n "${HAS_MPI_LIB}"; then
			mpi_lib="-L${mpi_path}/lib -lmpich -lpmpich"
			mpi_bin="pnoshell"
		else
			AC_MSG_ERROR([Error...Could not find MPICH library in ${mpi_path}/lib], 0)
		fi
	fi
	
	if test -n "${USE_MPICH2}"; then
		AC_CHECK_FILE(${mpi_path}/lib/libmpich.a, [HAS_MPI_LIB="yes"], [])
		AC_CHECK_FILE(${mpi_path}/lib/libmpich.so, [HAS_MPI_LIB="yes"], [])
		AC_CHECK_FILE(${mpi_path}/lib/libmpich.dylib, [HAS_MPI_LIB="yes"], [])
		if test -n "${HAS_MPI_LIB}"; then
			mpi_lib="-L${mpi_path}/lib -lmpich"
		else
			AC_MSG_ERROR([Error...Could not find MPICH library in ${mpi_path}/lib], 0)
		fi
	fi
	
	if test -n "${USE_LAM}"; then
		AC_CHECK_FILE(${mpi_path}/lib/liblam.a, [HAS_MPI_LIB="yes"], [])
		AC_CHECK_FILE(${mpi_path}/lib/liblam.so, [HAS_MPI_LIB="yes"], [])
		AC_CHECK_FILE(${mpi_path}/lib/liblam.dylib, [HAS_MPI_LIB="yes"], [])
		if test -n "${HAS_MPI_LIB}"; then
			mpi_lib="-L${mpi_path}/lib -lmpi -llam"
		else
			AC_MSG_ERROR([Error...Could not find LAM-MPI libraries in ${mpi_path}/lib], 0)
		fi
	fi
	
	if test -n "${USE_OPMI}"; then
		AC_CHECK_FILE(${mpi_path}/lib/libmpi.a, [HAS_MPI_LIB="yes"], [])
		AC_CHECK_FILE(${mpi_path}/lib/libmpi.so, [HAS_MPI_LIB="yes"], [])
		AC_CHECK_FILE(${mpi_path}/lib/libmpi.dylib, [HAS_MPI_LIB="yes"], [])
		if test -n "${HAS_MPI_LIB}"; then
			mpi_lib="-L${mpi_path}/lib -lmpi"
		else
			AC_MSG_ERROR([Error...Could not find OpenMPI libraries in ${mpi_path}/lib], 0)
		fi
	fi

	USE_MPI="yes"
fi

if test -n "${USE_MPI}"; then
	AC_DEFINE([HAVE_MPI_H], [], [have the MPI library])
	if test -n "${apbs_do_python}"; then
		AC_MSG_NOTICE([Sorry -- Python wrappers do not work with MPI.])
		AC_MSG_NOTICE([Disabling Python support])
		apbs_do_python=""
	fi
fi

AC_SUBST(mpi_lib)
AC_SUBST(mpi_inc)
dnl AC_SUBST(mpi_bin)

AC_MSG_NOTICE([---------------------------------------------------------------])
AC_MSG_NOTICE([---------- begin configuring for Python wrappers --------------])
AC_MSG_NOTICE([---------------------------------------------------------------])

PY_ENABLE=""
AC_MSG_CHECKING([for --enable-python])
if test x"$enable_python" != "xno"; then
	AC_MSG_RESULT([yes])
	AC_PATH_PROG(py_path, python)
	PY_ENABLE="yes"    

	if test -n "${py_path}"; then
		changequote(<<, >>)dnl
		PY_VERSION=`$py_path -c 'import sys; print sys.version[:3]'`   
		changequote([, ])dnl
		PY_PREFIX=`$py_path -c 'import sys; print sys.prefix'`
		PY_EXEC_PREFIX=`$py_path -c 'import sys; print sys.exec_prefix'`
		if test -f $PY_PREFIX/include/python$PY_VERSION/Python.h; then
			PY_CFLAGS="-I$PY_PREFIX/include/python$PY_VERSION -I$PY_EXEC_PREFIX/include/python$PY_VERSION"
		fi
	fi

	PY_LINKFLAGS=${CFLAGS}
	PY_CC=${CC}

	dnl # Platform specific linking flags

	PY_LDFLAGS=""
	PY_SHARED=""
	case $host_os in
		linux*) 
			PY_SHARED="-shared"
			CC="${CC} -fPIC"
			F77="${F77} -fPIC"
			;;
		irix*) 
			PY_SHARED="-shared"
			PY_LDFLAGS="-lftn -lm"
			;;
		sun|solaris**) 
			PY_SHARED="-G"
			PY_LDFLAGS="-lF77 -lm -lM77 -lsunmath"
			;;
		darwin**)
	 	PY_CC=${F77}
                       PY_SHARED="-shared"
                       if test "${F77}" == "g95"; then
                         PY_SHARED="-dynamiclib"
                       fi
                       PY_LINKFLAGS="-Wl,-framework -Wl,Python"
                       PY_LDFLAGS="-Wl,-framework -Wl,Accelerate"
                       ;;
		mingw*)
			PY_PREFIX=`$py_path -c 'import sys, string; a=sys.prefix; b=string.replace(a,"c:\\\","/c/"); print b;'`
			PY_EXEC_PREFIX=`$py_path -c 'import sys, string; a=sys.exec_prefix; b=string.replace(a,"c:\\\","/c/"); print b;'`
			changequote(<<, >>)dnl
			PY_NOSPACE_VER=`$py_path -c 'import sys; a=sys.version[0]; b=sys.version[2]; print a+b;'`
			changequote([, ])dnl
			PY_CC=${F77}
			PY_LDFLAGS="-L$PY_PREFIX/libs -lpython$PY_NOSPACE_VER"
			PY_SHARED="-shared"
			if test -f $PY_PREFIX/include/Python.h; then
				PY_CFLAGS="-I$PY_PREFIX/include -I$PY_EXEC_PREFIX/include"
			fi
			;;
	esac

	dnl # The Python wrappers use different flags on different arches

	IS_STATIC=`echo ${LDFLAGS} | grep static`
	case $host_cpu in
		x86_64*)
			PY_CC="${F77}";;
		i*86)
			PY_CC="${F77}"
			if test "${CC}" = "icc"; then
				AC_CHECK_LIB(ifcore_pic,d_int_val,[ifcore_lib="-lifcore_pic"],[ifcore_lib=""])
				if test -n "${ifcore_lib}"; then
					PY_LDFLAGS=${ifcore_lib}
				else
					AC_MSG_WARN([Unable to find Intel's ifcore_pic library - Python wrappers may fail looking for d_int_val.])
				fi
				if test -n "${IS_STATIC}"; then
					PY_LINKFLAGS="$CFLAGS -static-libcxa -Wl,-static"
				fi
			fi;; 
		ia64*)
			if test "${CC}" = "icc"; then
				AC_CHECK_LIB(ifcore_pic,d_int_val,[ifcore_lib="-lifcore_pic"],[ifcore_lib=""])
				if test -n "${ifcore_lib}"; then
					PY_LDFLAGS=${ifcore_lib}
				else
					AC_MSG_WARN([Unable to find Intel's ifcore_pic library - Python wrappers may fail looking for d_int_val.])
				fi
				if test -n "${IS_STATIC}"; then
					if (test -f /usr/lib/libc.so) && (test -f /usr/lib/libm.so); then
						PY_LINKFLAGS="$CFLAGS -static-libcxa -Wl,-static,/usr/lib/libc.so,/usr/lib/libm.so"
					else
						AC_MSG_WARN([Could not find libc.so or libm.so in /usr/lib])
					fi 
				fi
			fi;; 
	esac

else
	AC_MSG_RESULT([no])
fi


if test -n "${PY_ENABLE}"; then
	AC_MSG_NOTICE(Python wrappers enabled)
else
	AC_MSG_NOTICE(Python wrappers disabled)
fi

AC_SUBST(PY_SHARED,[$PY_SHARED])
AC_SUBST(PY_LDFLAGS,[$PY_LDFLAGS])
AC_SUBST(PY_CFLAGS,[$PY_CFLAGS])
AC_SUBST(PY_LINKFLAGS,[$PY_LINKFLAGS])
AC_SUBST(PY_CC,[$PY_CC])
AM_CONDITIONAL(PY_PRESENT, [test -n "${PY_ENABLE}"])

dnl # ---------------------------
dnl # BUILD PYTHON OPAL FILES
dnl # ---------------------------

AC_ARG_WITH(python,
  [  --with-python=PYTHONPATH           Specify python location.],
  [PYTHONPATH=$withval],
  [PYTHONPATH=""])

if test "$PYTHONPATH" = "" ; then
  AC_PATH_PROG(py_path, python)
  if test -z "${py_path}"; then
	 AC_MSG_RESULT([no])
     AC_MSG_ERROR([Unable to find python!  Python location can be specified with --with-python])
  fi
else
  py_path=$PYTHONPATH
fi

changequote(<<, >>)dnl
PY_VERSION=`$py_path -c 'import sys; print sys.version[:3]'`
changequote([, ])dnl

SERVICEURL='http://kryptonite.nbcr.net/opal2/services/apbs_1.3'
PARALLELSERVICEURL='http://oolite.calit2.optiputer.net/opal2/services/apbs-parallel-1.3'
SITEPKGS=${prefix}/lib/python${PY_VERSION}/site-packages

AC_ARG_WITH(maxmem,
  [  --with-maxmem=MAXMEM                    Maximum memory (in MB) to be used for a job.],
  [MAXMEM=$withval],
  [MAXMEM="-1"])

sed -e 's%@MAXMEM@%'${MAXMEM}'%' \
    -e 's%@SITEPKGS@%'${SITEPKGS}'%' \
    -e 's%@PACKAGE_VERSION@%'${PACKAGE_VERSION}'%' \
    -e 's%@SERVICEURL@%'${SERVICEURL}'%' \
    -e 's%@PARALLELSERVICEURL@%'${PARALLELSERVICEURL}'%' ${top_apbsdir}/bin/ApbsClient.py.in > ${top_apbsdir}/bin/ApbsClient.py

AC_SUBST(SITEPKGS,[$SITEPKGS])

dnl # ---------------------------
dnl # BUILD MALOC IF NECESSARY
dnl # ---------------------------
if test -z "${USE_FETK}";  then
	AC_MSG_CHECKING([for --disable-maloc-rebuild])
	AC_ARG_ENABLE(
		[maloc-rebuild],
		AS_HELP_STRING([--disable-maloc-rebuild], [disable rebuilding maloc]),,
		[enable_maloc_rebuild=yes]
	)
	
	if test x"$enable_maloc_rebuild" = "xno"; then
		AC_MSG_RESULT([yes])
		AC_MSG_NOTICE([Will not rebuild MALOC...skipping])
		apbs_do_maloc=""
	else
		AC_MSG_RESULT([no])
		AC_MSG_NOTICE([-----------------------------------------------------------])
		AC_MSG_NOTICE([-------------------- building MALOC -----------------------])
		AC_MSG_NOTICE([-----------------------------------------------------------])
		apbs_do_maloc="maloc"
		
		export CC=${CC}
		export CFLAGS=${CFLAGS}
		export LDFLAGS=${LDFLAGS}
		cd ${contrib_path}/maloc
	
		make clean
	
		if test -n "${USE_MPI}";then
			export APBS_MPI_INC=${mpi_inc}
			export APBS_MPI_LIB=${mpi_lib}
			
			./configure --prefix=${contrib_path} --with-mpi --enable-vdebug
		else
			./configure --prefix=${contrib_path} --enable-vdebug
		fi
		
		make; make install;
		cd ../../

		AC_MSG_NOTICE([-----------------------------------------------------------])
		AC_MSG_NOTICE([----- MALOC installation finished continuing with APBS ----])
		AC_MSG_NOTICE([-----------------------------------------------------------])
	fi
	
else
	AC_MSG_NOTICE([Skipping MALOC installation; will use MALOC provided with FEtk])
fi
AC_SUBST(apbs_do_maloc)

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
	[${maloc_lib_path} -lmaloc -lvf2c ${liberty_lib} ${socket_lib} ${nsl_lib} ${math_lib}])
if test -n "${rl_use}"; then

	dnl # READLINE library location ENV specification
	AC_MSG_CHECKING([whether your environment defines FETK_RL_LIBRARY])
	if test -n "${FETK_RL_LIBRARY}"; then
		AC_MSG_RESULT([${FETK_RL_LIBRARY}])
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
		[rl_use="yes";rl_lib="${rl_lib_path} -lreadline ${ncurses_lib}"],
		[rl_use="";rl_lib=""],
		[${rl_lib_path} ${ncurses_lib} ${math_lib}])
	if test -n "${rl_use}"; then
		AC_CHECK_LIB(readline,add_history,
			[rl_use="yes";rl_lib="${rl_lib_path} -lreadline ${ncurses_lib}"],
			[rl_use="";rl_lib=""],
			[${rl_lib_path} ${ncurses_lib} ${math_lib}])
	fi

	dnl # Return LDFLAGS to normal
	LDFLAGS=${LDFLAGS2_SAVE};

	dnl # Final decision
	if test -z "${rl_use}"; then
		AC_MSG_NOTICE([MALOC was compiled with READLINE enabled])
		AC_MSG_NOTICE([MALOC needs to link against the READLINE library])
		AC_MSG_ERROR([READLINE cannot be found!])
	fi
fi
LDFLAGS=${LDFLAGS_SAVE};
CPPFLAGS=${CPPFLAGS_SAVE};
AC_SUBST(rl_lib)

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

AC_MSG_NOTICE([------------------------------------------------------------])
AC_MSG_NOTICE([---------- begin configuring for FEtk linkage --------------])
AC_MSG_NOTICE([------------------------------------------------------------])

dnl # BEGIN FETK LIBRARY/HEADER TEST BLOCK
mc_lib=""
mc_inc=""
mc_do_src=""
mc_do_lib=""
mc_do_test=""

AC_MSG_CHECKING([for --disable-fetk])
AC_ARG_ENABLE(
	[fetk],
	AS_HELP_STRING([--disable-fetk],[disable FEtk features]),,
	[enable_fetk=yes]
)

if test x"$enable_fetk" != "xno"; then
	AC_MSG_RESULT([no])
	AC_MSG_CHECKING([for user-specified FEtk include path])
	if test -n "${fetk_include_path}"; then
		AC_MSG_RESULT([${fetk_include_path}])
		mc_inc_path="-I${fetk_include_path} ";
	else
		AC_MSG_RESULT([no])
	fi
	
	AC_MSG_CHECKING([for user-specified FEtk library path])
	if test -n "${fetk_library_path}"; then
		AC_MSG_RESULT([${fetk_library_path}])
		mc_lib_path="-L${fetk_library_path}";
	else
		AC_MSG_RESULT([no])
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
		[mc_found="yes";mc_lib="${mc_lib_path} -lmc"],
		[mc_found="";mc_lib=""],
		[${mc_lib_path} -lmc ${maloc_lib} ${mpi_lib} ${rl_lib} ${liberty_lib} ${socket_lib} ${nsl_lib} ${math_lib}])
	dnl # Setup temp LDFLAGS and look for MCX library/header
	LDFLAGS=${mc_lib_path};
	CPPFLAGS=${mc_inc_path};
	mcx_found="";
	AC_CHECK_HEADER(mcx/mcx.h,
		[mcx_found="yes";mcx_inc="${mc_inc_path}"],
		[mcx_found="";mc_inc=""])
	AC_CHECK_LIB(mcx,mcx_link,
		[mcx_found="yes";mcx_lib="${mc_lib_path} -lmcx"],
		[mcx_found="";mcx_lib=""],
		[${mc_lib_path} -lmc ${maloc_lib} ${mpi_lib} ${rl_lib} ${liberty_lib} ${socket_lib} ${nsl_lib} ${math_lib}])
    

	dnl # Return LDFLAGS and CPPFLAGS to normal
	LDFLAGS=${LDFLAGS_SAVE}
	CPPFLAGS=${CPPFLAGS_SAVE}
	
	if test -n "${mc_found}"; then
		AC_MSG_NOTICE([FEtk component MC found:  enabled.])
		AC_DEFINE([HAVE_MC_H], [], [have FEtk component MC])
	else
		AC_MSG_WARN([FEtk component MC not found:  finite element solver disabled.])
	fi
	if test -n "${mcx_found}"; then
		AC_MSG_NOTICE([FEtk component MCX found:  enabled.])
		AC_DEFINE([HAVE_MCX_H], [], [have FEtk component MCX])\
		AC_DEFINE([USE_HB], [], [use hierarchical basis method])
	else
		AC_MSG_WARN([FEtk component MCX not found:  hierarchical basis methods disabled.])
	fi
else
	AC_MSG_RESULT([yes])
	AC_MSG_NOTICE([FEtk disabled at user request.])
fi
AC_SUBST(mcx_lib)
AC_SUBST(mcx_inc)
AC_SUBST(mc_lib)
AC_SUBST(mc_inc)

AC_MSG_NOTICE([---------------------------------------------------------------])
AC_MSG_NOTICE([------------ begin processing APBS configure options ----------])
AC_MSG_NOTICE([---------------------------------------------------------------])

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
AC_ARG_ENABLE(
	[pansi],
	AS_HELP_STRING([--enable-pansi], [enable pedantic ANSI compilation]),,
	[enable_pansi=no]
)
if test x"$enable_pansi" != "xno"; then
	AC_MSG_RESULT([yes])
	wall="-Wall"
	wall_pedantic="-Wall -pedantic"
	wall_pedantic_ansi="-Wall -pedantic -ansi"
	pedantic_ansi="-pedantic -ansi"
else
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
		if test "${CC}" != "icl"; then
		  AC_DEFINE([HAVE_STRCASECMP], [], [have machine-supplied strcasecmp])
		fi
	fi
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
AC_ARG_ENABLE(
	[tools],
	AS_HELP_STRING([--disable-tools],[disable tool compilation]),,
	[enable_tools=yes]
)
apbs_do_tools=""
if test x"$enable_tools" != "xno"; then
	AC_MSG_RESULT([no])
	AC_MSG_NOTICE([tools will be compiled])
	apbs_do_tools="tools"
else
	AC_MSG_RESULT([yes])
	AC_MSG_NOTICE([tools will not be compiled])
fi
AC_SUBST(apbs_do_tools)


dnl # -----------------------
dnl # HANDLE THE BLAS LIBRARIES
dnl # -----------------------
blas_lib=""
AC_ARG_WITH(
	    [blas], 
	    AS_HELP_STRING([--with-blas=LDFLAGS], [link to user-specified BLAS (specify flags to pass to ld)])
	    )
AC_MSG_CHECKING([for --with-blas])
if test x"$with_blas" != "xno" -a -n "$with_blas"; then
	blas_lib="${with_blas}"
	AC_MSG_RESULT([${blas_lib}])
 	AC_ARG_WITH(
		    [blas_name], 
		    AS_HELP_STRING([--with-blas-name=NAME], [name of BLAS library, defaults to blas]),, 
		    [with_blas_name="blas"]
		    )
	if test x"$with_blas_name" != "xno" -a -n "$with_blas_name"; then
		blas_name="${with_blas_name}"
	fi
 	AC_MSG_NOTICE([checking ${blas_name} using ${blas_lib}])
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
	apbs_do_blas=""
	AC_LANG_C
else
	AC_MSG_RESULT([no])
	apbs_do_blas="blas"
	blas_lib=${contrib_path}/blas/.libs/libapbsblas.a
fi
AC_SUBST(apbs_do_blas)

dnl # --------------------------------
dnl # HANDLE ZLIB COMPILATION
dnl # --------------------------------
disable_zlib=yes
z_lib=
AC_MSG_CHECKING([--disable-zlib])
AC_ARG_ENABLE(
	[zlib],
	AS_HELP_STRING([--disable-zlib], [disable zlib compilation]),,
	[disable_zlib=no]
)
if test x"$disable_zlib" == "xno"; then
	AC_MSG_RESULT([yes])
	zlib="zlib"
	z_lib=${contrib_path}/zlib/.libs/libz.a
	AC_DEFINE([HAVE_ZLIB], [], [have zlib compression enabled])
fi
AC_SUBST(zlib)
AC_SUBST(z_lib)

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
dnl # If using Intel compilers add in -nofor_main
dnl # -----------------------
case $F77 in
	ifort)
			FFLAGS="${FFLAGS} -nofor_main" ;;
	*)
		;;
esac

dnl # -----------------------
dnl # CHECK WHETHER OR NOT TO LINK IN ARPACK FUNCTIONS
dnl # -----------------------
AC_MSG_CHECKING([for --with-arpack])
AC_ARG_WITH(
	[arpack],
	AS_HELP_STRING([--with-arpack=PATH], [enable ARPACK functions (specify path to lib)])
)

if test x"$with_arpack" != "xno" -a -n "$with_arpack"; then
	arpack_lib="-larpack -L${with_arpack}"
	arpack_do_tools="arpack"
	AC_MSG_RESULT([${arpack_lib}])
else
	arpack_do_tools=""
	AC_MSG_RESULT([no])
fi
AC_SUBST(arpack_lib)
AC_SUBST(arpack_do_tools)

dnl # -----------------------
dnl # CHECK WHETHER YOU WANT INLINE FUNCTIONS
dnl # -----------------------
AC_MSG_CHECKING([for --disable-inline])
AC_ARG_ENABLE(
	[inline],
	AS_HELP_STRING([--disable-inline],[disable inline functions (macros)]),,
	[enable_inline=yes]
)
if test x"$enable_inline" != "xno"; then
	AC_MSG_RESULT([no])
	AC_MSG_NOTICE([inline functions will be used])
else
	AC_MSG_RESULT([yes])
	AC_MSG_NOTICE([inline functions will not be used])
	AC_DEFINE([APBS_NOINLINE], [1], [do not inline functions])
fi

dnl # -----------------------
dnl # CHECK WHETHER YOU WANT TO INCLUDE TINKER
dnl # -----------------------
AC_MSG_CHECKING([for --enable-tinker])
AC_ARG_ENABLE(
	[tinker],
	AS_HELP_STRING([--enable-tinker],[include TINKER support]),,
	[enable_tinker=no]
)
if test x"$enable_tinker" != "xno"; then
    AC_MSG_RESULT([yes])
    AC_DEFINE([WITH_TINKER], [], [include TINKER support])
	AC_DEFINE([VAPBSQUIET], [], [APBS diagnostic messages will not display])
else
    AC_MSG_RESULT([no]) 
fi

dnl # -----------------------
dnl # CHECK WHETHER YOU WANT TO DISABLE APBS MESSAGES
dnl # -----------------------
AC_MSG_CHECKING([for --enable-abps-quiet])
AC_ARG_ENABLE(
	[abps-quiet],
	AS_HELP_STRING([--enable-abps-quiet],[prevent APBS diagnostic messages from displaying to stdout]),,
	[enable_apbs_quiet=no]
)
if test x"$enable_apbs_quiet" != "xno"; then
    AC_MSG_RESULT([yes])
    AC_DEFINE([VAPBSQUIET], [], [APBS diagnostic messages will not display])
else
    AC_MSG_RESULT([no]) 
fi

dnl # -----------------------
dnl # CHECK WHETHER YOU WANT TO DISABLE OPENMP SUPPORT
dnl # -----------------------
AC_MSG_CHECKING([for --disable-openmp])
AC_ARG_ENABLE(
	[openmp],
	AS_HELP_STRING([--disable-openmp],[disables OpenMP compiler support]),,
	[disable_openmp=yes]
)
if test x"$disable_openmp" != "xyes"; then
    AC_MSG_RESULT([yes])
else
    AC_MSG_RESULT([no])
	if [[ "${CC}" = "icc" ] && [ "${F77}" = "ifort" ]]; then
		CFLAGS="${CFLAGS} -openmp"
		FFLAGS="${FFLAGS} -openmp -assume cc_omp"
		LDFLAGS="${LDFLAGS} -static-intel"
	fi
	
	if [[ "${CC}" = "gcc" ] && [ "${F77}" = "gfortran" ]]; then
		gcc_version=`gcc --version`
		gcc_complete_version=`echo ${gcc_version} | awk '{printf("%s",$3)}' | awk -F'.' '{print $1$2}'`
		
		if test ["${gcc_complete_version}" -ge 42]; then
			CFLAGS="${CFLAGS} -fopenmp"
			FFLAGS="${FFLAGS} -fopenmp"
		else
		    AC_MSG_RESULT([WARNING])
			AC_MSG_RESULT([WARNING Version of gcc/gfortran does not support OpenMP. Disabling.])
			AC_MSG_RESULT([WARNING])
		fi
	fi
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

dnl # -----------------------
dnl # ASSEMBLE A LIST OF LIBRARIES AND INCLUDE FILES THAT WILL BE USED IN ALL
dnl # COMPILATIONS
dnl # -----------------------
top_srcdir="`pwd`"
base_list="${FLIBS} ${tree_lib} ${vmetis_lib} ${vor_lib} ${mcx_lib} ${mc_lib} \ 
	   ${thread_lib} ${mpi_lib} ${math_lib} ${rl_lib} ${socket_lib} ${nsl_lib} ${blas_lib} ${z_lib}"

case $host_os in

	darwin8*)
		extra_list1="-lSystemStubs"
		;;
	darwin*)
		LDFLAGS="${LDFLAGS} -mmacosx-version-min=10.4"
		;;
esac

dnl # -------------------------
dnl # Need to add in a case statement to put in lcc_dynamic for g77
dnl # This is for OS X only
dnl # -------------------------
extra_list2=""
case $host_os in

	darwin*)
		if test "$F77" == "g77"; then
			extra_list2="-lcc_dynamic"
		else
			extra_list2=""
		fi
		;;
		
esac

dnl # -----------------------
dnl # ADD ON THE LIBRARIES FOR VERSION 1.40 OF FETK
dnl # -----------------------
if test -n "${mc_found}"; then
	case $host_os in

	*)
		extra_list3="-lsuperlu -lumfpack -lamd -lpunc -lblas"
		;;
	
	esac
else

	extra_list3=""

fi

inc_list="-I. ${maloc_inc} ${mc_inc} ${mcx_inc} ${mpi_inc}"
lib_list="${base_list} ${extra_list1} ${extra_list2} ${extra_list3}"
AC_SUBST(lib_list)
AC_SUBST(inc_list)

apbs_build_libraries=""
if [[ "${host_os}" = "cygwin" ] && [ "${CC}" = "icl" ]]; then
  maloc_lib=""
  AC_SUBST(maloc_lib)
  apbs_build_libraries="../bin/libapbsmainroutines.la \
                        ../src/aaa_lib/libapbs.la \
                        ../contrib/blas/libapbsblas.la \
			../contrib/maloc/src/aaa_lib/libmaloc.la \
			../../Program\ Files/Microsoft\ Visual\ Studio\ 8/VC/PlatformSDK/Lib/WSock32.Lib"
else
  apbs_build_libraries="${top_apbsdir}/bin/libapbsmainroutines.la \
			${top_apbsdir}/src/aaa_lib/libapbs.la \
			${maloc_lib} ${lib_list}"
fi
AC_SUBST(apbs_build_libraries)

dnl ##########################################################################
dnl 11. generate the makefiles
dnl ##########################################################################

AC_CONFIG_FILES([
    Makefile

    contrib/Makefile
    contrib/blas/Makefile
    contrib/pmgZ/Makefile
    contrib/aqua/Makefile
    contrib/opal/Makefile
    contrib/opal/ZSI/Makefile
    contrib/zlib/Makefile	
	
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
    src/aaa_inc/Makefile
    src/aaa_inc/apbs/Makefile
    src/aaa_lib/Makefile

    examples/Makefile
    examples/actin-dimer/Makefile
    examples/alkanes/Makefile
    examples/born/Makefile
    examples/FKBP/Makefile
    examples/hca-bind/Makefile
    examples/ion-pmf/Makefile
    examples/ion-protein/Makefile
    examples/ionize/Makefile
    examples/pka-lig/Makefile
    examples/point-pmf/Makefile
    examples/protein-rna/Makefile
    examples/smpbe/Makefile
    examples/solv/Makefile
    examples/opal/Makefile
    
    bin/Makefile

])
AC_OUTPUT

AC_MSG_NOTICE([

APBS to be compiled with following settings:  ${CC} ${CFLAGS}
APBS to be compiled with following settings:  ${F77} ${FFLAGS}
APBS to be linked with following libraries :  ${lib_list}
APBS to be compiled with following includes:  ${inc_list}

Please type "make; make install"
])


