## ###########################################################################
## File:    Makefile.am
##
## Purpose: Automatically generate Makefile.in for autoconf using automake.
##
## Author:  Nathan Baker, Todd Dolinsky
## ###########################################################################

AUTOMAKE_OPTIONS = no-dependencies
SUFFIXES         = .c .C .o .so
INCLUDES         = @inc_list@ $(PY_CFLAGS) -I../../bin 
CCLD             = $(CC)
CLEANFILES       = apbslibc.so

if PY_PRESENT
   noinst_PROGRAMS       = apbslibc.so
   apbslibc_so_SOURCES   = apbslib.c
   apbslibc_so_LDADD     = ../../bin/libapbsmainroutines.la ../../src/aaa_lib/libapbs.la @lib_list@   
   apbslibc_so_LDFLAGS   = -shared 
endif
