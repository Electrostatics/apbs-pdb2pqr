## ###########################################################################
## File:    Makefile.am
##
## Purpose: Automatically generate Makefile.in for autoconf using automake.
##
## Author:  Nathan Baker, Todd Dolinsky
## ###########################################################################

AUTOMAKE_OPTIONS = no-dependencies
SUFFIXES         = .c .C .o .so
INCLUDES         = @inc_list@ $(PY_CFLAGS) -I${top_srcdir}/bin 
CCLD             = $(CC)
CLEANFILES       = apbslibc.so

if PY_PRESENT
   noinst_PROGRAMS       = apbslibc.so
   apbslibc_so_SOURCES   = apbslib.c
   apbslibc_so_LINK      = $(CC) $(PY_LINKFLAGS)
   apbslibc_so_LDFLAGS   = -o $@ $(PY_SHARED)
   apbslibc_so_LDADD     = ${top_builddir}/bin/.libs/libapbsmainroutines.a ${top_builddir}/src/aaa_lib/.libs/libapbs.a $(PY_LDFLAGS) @lib_list@
endif
