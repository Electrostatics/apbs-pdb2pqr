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
AC_DEFUN(AC_F77_FCN_MANGLE,
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
                AC_DEFINE(VF77_NOUNDERSCORE)
                AC_DEFINE(VF77_LOWERCASE)
                mangle_try=foo_bar_
                ;;
        lowercase-underscore)
                AC_DEFINE(VF77_ONEUNDERSCORE)
                AC_DEFINE(VF77_LOWERCASE)
                mangle_try=foo_bar__
                ;;
        uppercase)
                AC_DEFINE(VF77_NOUNDERSCORE)
                AC_DEFINE(VF77_UPPERCASE)
                mangle_try=FOO_BAR_
                ;;
        uppercase-underscore)
                AC_DEFINE(VF77_ONEUNDERSCORE)
                AC_DEFINE(VF77_UPPERCASE)
                mangle_try=FOO_BAR__
                ;;
esac

])




