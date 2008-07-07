dnl acinclude.m4 generated manually by Nathan Baker

dnl Test for the name mangling scheme used by the Fortran 77 compiler.
dnl Two variables are set by this macro:
dnl
dnl        f77_case: Set to either "upper" or "lower", depending on the
dnl                  case of the name mangling.
dnl
dnl  f77_underscore: Set to either "no", "single" or "double", depending
dnl                  on how underscores (i.e. "_") are appended to
dnl                  identifiers, if at all.
dnl
dnl                  If no underscores are appended, then the value is
dnl                  "no".
dnl
dnl                  If a single underscore is appended, even with
dnl                  identifiers which already contain an underscore
dnl                  somewhere in their name, then the value is
dnl                  "single".
dnl
dnl                  If a single underscore is appended *and* two
dnl                  underscores are appended to identifiers which
dnl                  already contain an underscore somewhere in their
dnl                  name, then the value is "double".
dnl
dnl   This was adapted by Nathan Baker from an e-mail by Steven G. Johnson on
dnl   the Autoconf mailing list.  
dnl
AC_DEFUN([AC_F77_FCN_MANGLE],
[
AC_REQUIRE([AC_PROG_CC])
AC_REQUIRE([AC_PROG_F77])
AC_MSG_CHECKING(how f77 mangles function names)
cat > mangle-func.f <<EOF
      subroutine foobar()
      return
      end
      subroutine foo_bar()
      return
      end
EOF
ac_try='$F77 -c $FFLAGS mangle-func.f 1>&AC_FD_CC'
if AC_TRY_EVAL(ac_try); then
  ac_try=""
else
  echo "configure: failed program was:" >&AC_FD_CC
  cat mangle-func.f >&AC_FD_CC
  rm -f mangle-func*
  AC_MSG_ERROR(failed to compile fortran test program)
fi

ac_f77_mangle_type=unknown
AC_LANG_SAVE
AC_LANG_C
ac_save_LIBS="$LIBS"
LIBS="mangle-func.o $LIBS"
AC_TRY_LINK(,foobar();,
     ac_f77_mangle_type=lowercase,
     AC_TRY_LINK(,foobar_();,
          ac_f77_mangle_type=lowercase-underscore,
          AC_TRY_LINK(,FOOBAR();,
               ac_f77_mangle_type=uppercase,
               AC_TRY_LINK(,FOOBAR_();,
                    ac_f77_mangle_type=uppercase-underscore))))
LIBS="$ac_save_LIBS"
AC_LANG_RESTORE
AC_MSG_RESULT($ac_f77_mangle_type)

case $ac_f77_mangle_type in
        unknown)
                AC_MSG_ERROR(unknown fortran name-mangling scheme)
                ;;
        lowercase)
                AC_DEFINE([VF77_NOUNDERSCORE], [], [no F77 underscore])
                AC_DEFINE([VF77_LOWERCASE], [], [F77 lowercase])
                mangle_try=foo_bar_
                ;;
        lowercase-underscore)
                AC_DEFINE([VF77_ONEUNDERSCORE], [], [one F77 underscore])
                AC_DEFINE([VF77_LOWERCASE], [], [F77 lowercase])
                mangle_try=foo_bar__
                ;;
        uppercase)
                AC_DEFINE([VF77_NOUNDERSCORE], [], [no F77 underscore])
                AC_DEFINE([VF77_UPPERCASE], [], [F77 uppercase])
                mangle_try=FOO_BAR_
                ;;
        uppercase-underscore)
                AC_DEFINE([VF77_ONEUNDERSCORE], [], [one F77 underscore])
                AC_DEFINE([VF77_UPPERCASE], [], [F77 uppercase])
                mangle_try=FOO_BAR__
                ;;
esac

])

rm -f mangle-func.f mangle-func.o

dnl Test for floating point error found in Intel compilers without
dnl   -mp flag (and maybe other compilers as well?)
dnl
dnl   by Todd Dolinsky
dnl   modified by David Gohara
dnl   NOTE: The update just calculates the machine epsilon and uses
dnl		for all systems when calling VFLOOR (only done in nosh)
dnl

AC_DEFUN([AC_FPERROR], [
AC_REQUIRE([AC_PROG_CC])
AC_MSG_CHECKING(the machine epsilon for VFLOOR)

cat > fperror.c <<EOF
#include <stdio.h>

int main( int argc, char **argv )
{
        float machEps = 1.0f;
        
        while(1)
        {
                machEps /= 2.0;
		/*
                 If next epsilon yields 1, then break, because current
                 epsilon is the machine epsilon.
		*/
                if( (float)(1.0 + (machEps/2.0)) == 1.0 )
                        break;  
        }
        
        printf( "%G\n", machEps );
        return 0;
}
EOF 
actry='$CC $CFLAGS -c fperror.c 1>&AC_FD_CC'
if AC_TRY_EVAL(actry); then
    if test "${CC}" = "icc"; then
        $CC -mp fperror.c -o fperror.test -lm
    else
        $CC fperror.c -o fperror.test -lm
    fi
    FPRESULTS=`./fperror.test`  
    if test -n "$FPRESULTS"; then
        AC_MSG_RESULT([yes])
        echo "The machine epsilon for this system is: $FPRESULTS"
         if test "${CC}" = "icc"; then	 
            echo "*** For icc add -mp flag if desired (see icc manual) ***"	 
         fi
        AC_DEFINE([FLOAT_ERRORS], 1, [have floating point errors])
        AC_DEFINE_UNQUOTED([MACHINE_EPS], $FPRESULTS, [machine error])
    else
        AC_MSG_RESULT([no])
        AC_DEFINE([FLOAT_ERRORS], 0, [have floating point errors])
    fi
    rm -f fperror*
else
    AC_MSG_RESULT([unknown])
    echo "configure: failed program was:" >&AC_FD_CC
    cat fperror.c >&AC_FD_CC
    rm -f fperror*
    AC_MSG_ERROR(failed to compile test C program)
fi
])

dnl Test for the use of the -nofor_main option used by Alpha FORTRAN
dnl
dnl   by Nathan Baker 
dnl
AC_DEFUN([AC_F77_NOFORMAIN], [
AC_REQUIRE([AC_PROG_CC])
AC_REQUIRE([AC_PROG_F77])
dnl AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
AC_REQUIRE([AC_CANONICAL_HOST])
dnl First, compile the C code with the main routine
cat > noformain.c <<EOF
    int main(int argc, char **argv) {
        return 0;
    }
EOF
AC_LANG_SAVE
AC_LANG_C
ac_try='$CC -c $CFLAGS noformain.c 1>&AC_FD_CC'
if AC_TRY_EVAL(ac_try); then
  ac_try=""
else
  echo "configure: failed program was:" >&AC_FD_CC
  cat noformain.c >&AC_FD_CC
  rm -f noformain*
  AC_MSG_ERROR(failed to compile test C program)
fi
AC_MSG_CHECKING([if Fortran 77 compiler accepts -nofor_main])
AC_LANG_FORTRAN77
myflibs="$FLIBS"
ac_try='$F77 -nofor_main $FFLAGS $FLIBS noformain.o -o noformain 1>&AC_FD_CC'
if AC_TRY_EVAL(ac_try); then
  ac_try=""
  AC_MSG_RESULT([yes])
  AC_MSG_CHECKING([if Fortran 77 compilers requires -nofor_main])
  AC_LANG_FORTRAN77
  myflibs="$FLIBS"
  ac_try='$F77 $FFLAGS $FLIBS noformain.o -o noformain
  1>&AC_FD_CC'
  if AC_TRY_EVAL(ac_try); then
    ac_try=""     
    AC_MSG_RESULT([no])
  else
    AC_MSG_RESULT([yes])
    FFLAGS="$FFLAGS -nofor_main"
    AC_SUBST(FFLAGS)
  fi
else
  AC_MSG_RESULT([no])
fi
rm -f noformain*
AC_LANG_RESTORE

])

