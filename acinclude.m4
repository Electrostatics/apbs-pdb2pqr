dnl acinclude.m4 generated manually by Nathan Baker

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

