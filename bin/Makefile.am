## ###########################################################################
## File:    Makefile.am
##
## Purpose: Automatically generate Makefile.in for autoconf using automake.
##
## Author:  Nathan Baker
## ###########################################################################

AUTOMAKE_OPTIONS = no-dependencies
SUFFIXES = .H .h .C .cc .c .f .o
CCLD = $(F77)

bindir = ${prefix}/bin/${fetk_cpu_vendor_os}
bin_PROGRAMS     = apbs
apbs_SOURCES     = main.c routines.c
apbs_LDADD       = ../src/aaa_lib/libapbs.la @lib_list@

libdir = ${prefix}/lib/${fetk_cpu_vendor_os}
lib_LTLIBRARIES = libapbsmainroutines.la
libapbsmainroutines_la_SOURCES = routines.c
libapbsmainroutines_la_LDFLAGS = -version-info ${APBS_VERSION}

INCLUDES         = @inc_list@
AM_CFLAGS        = @wall_pedantic@
AM_LDFLAGS          = @LDFLAGS@ @FLIBS@
