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
CCLD             = $(PY_CC)
CLEANFILES       = vgridc.so

if PY_PRESENT
   noinst_PROGRAMS       = vgridc.so
   vgridc_so_SOURCES   = vgridlib.c
   vgridc_so_LINK      = $(PY_CC) $(PY_LINKFLAGS)
   vgridc_so_LDFLAGS   = -o $@ $(PY_SHARED)
   vgridc_so_LDADD     = ${top_builddir}/bin/.libs/libapbsmainroutines.a ${top_builddir}/src/aaa_lib/.libs/libapbs.a $(PY_LDFLAGS) @lib_list@
endif
