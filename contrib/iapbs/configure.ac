AC_INIT([iAPBS],[2.0.0],[rok@ucsd.edu])
AC_PREREQ(2.65)

AC_CONFIG_AUX_DIR(config)
AC_PREFIX_DEFAULT(${HOME})
AC_CANONICAL_TARGET

dnl AM_INIT_AUTOMAKE(iapbs, 2.0.0)
AM_INIT_AUTOMAKE
dnl AM_CONFIG_HEADER(src/iapbscfg.h)

AC_SUBST(IAPBS_VERSION, 2.0.0)
ACLOCAL=":";
AUTOCONF=":";
AUTOMAKE=":";
AUTOHEADER=":";
AC_SUBST(ACLOCAL)
AC_SUBST(AUTOCONF)
AC_SUBST(AUTOMAKE)
AC_SUBST(AUTOHEADER)
AC_CONFIG_MACRO_DIR([m4])

cclist="icc gcc mpicc"
fclist="ifort gfortran g95 g77 mpif90"

AC_PROG_CC(${cclist})
AC_PROG_FC(${fclist})
AC_PROG_F77(${fclist})
AC_PROG_INSTALL
AC_PROG_CPP
#AC_PROG_RANLIB
# libtool needs to be installed
##AM_PROG_LIBTOOL
#AC_LIBTOOL_DLOPEN
#AC_PROG_LIBTOOL
LT_INIT

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

AC_MSG_CHECKING([if F77 is Intel compiler])
is_intel=`$F77 --version| grep ifort`

# options for F77 preprocessor
if test  -n "${is_intel}"; then
 AC_SUBST([FPP_FLAGS], ["-cpp"])
 AC_MSG_RESULT([yes])
else
 AC_MSG_RESULT([no])
fi

AC_MSG_CHECKING([if F77 is GNU compiler])
if test  $ac_cv_f77_compiler_gnu = yes; then
 AC_SUBST([FPP_FLAGS], ["-x f77-cpp-input"])
 AC_MSG_RESULT([yes])
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




#if test $ac_cv_fc_compiler_gnu = yes; then
#   AC_SUBST([MY_CPPFLAGS], ["-x f95-cpp-input"])
#else
#   AC_SUBST([CPPFLAGS], ["-x c"])
#AC_SUBST([AM_FFLAGS], ["-cpp -DHAVE_MPI_H"])
#fi


# check for --with-mpi option
AC_MSG_CHECKING([for --with-openmpi option])
AC_ARG_WITH(
   [openmpi],
   AS_HELP_STRING([--with-openmpi], [enable OpenMPI compilation]),
      [mpi_path="${with_openmpi}"],
   )
   if test x"$mpi_path" != "xno" -a -n "$mpi_path"; then
      AC_MSG_RESULT([${mpi_path}])
      USE_OPMI="yes"
   else
      AC_MSG_RESULT([no])
      mpi_path=
   fi

if test -n "${mpi_path}"; then
 AC_DEFINE([HAVE_MPI_H], [1], [have the MPI library])
# AC_SUBST([AM_FFLAGS], ["-cpp -DHAVE_MPI_H"])
fi


# check for MALOC/APBS libs and includes
AC_MSG_CHECKING([for APBS_PREFIX])
if test -n "${APBS_PREFIX}"; then
    AC_MSG_RESULT([yes])
    apbs_inc_path="-I${APBS_PREFIX}/include ";
    apbs_lib_path="-L${APBS_PREFIX}/lib"
    apbs_libs=" -lapbsmainroutines -lapbs -lmaloc"
else
    AC_MSG_ERROR([MALOC/APBS libraries or headers missing: exiting on error!])
fi

# check for BLAS
AC_MSG_CHECKING([for APBS_BLAS])
if test -n "${APBS_BLAS}"; then
    AC_MSG_RESULT([yes])
    apbs_blas="${APBS_BLAS}"
else
    AC_MSG_RESULT([no, using MALOC BLAS (libapbsblas.a)])
#    apbs_blas="-L${APBS_PREFIX}/lib -lapbsblas"
    apbs_blas="-lapbsblas"
fi

# here I should actually check if libs/ includes and BLAS are working

AC_SUBST(apbs_inc_path)
AC_SUBST(apbs_lib_path)
AC_SUBST(apbs_libs)
AC_SUBST(apbs_blas)

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
#	z_lib=${contrib_path}/zlib/.libs/libz.a
#	z_lib="-L${APBS_PREFIX}/lib -lz"
	z_lib="-lz"
	AC_DEFINE([HAVE_ZLIB], [], [have zlib compression enabled])
fi

AC_SUBST(zlib)
AC_SUBST(z_lib)


#AC_CONFIG_FILES([Makefile src/Makefile test/Makefile])
AC_OUTPUT
AC_CONFIG_FILES([
	Makefile
	src/Makefile
	test/Makefile
])

AC_OUTPUT

AC_MSG_NOTICE([

iAPBS to be compiled with following settings:  ${CC} ${CFLAGS}
iAPBS to be compiled with following settings:  ${F77} ${FFLAGS}
iAPBS to be linked with following libraries :  ${apbs_lib_path} ${apbs_libs} ${apbs_blas} ${z_lib}
iAPBS to be compiled with following includes:  ${apbs_inc_path}

Please type "make; make install"
])
