## ###########################################################################
## File:    Makefile.am
##
## Purpose: Automatically generate Makefile.in for autoconf using automake.
##
## Author:  Nathan Baker
## ###########################################################################

AUTOMAKE_OPTIONS = no-dependencies
SUFFIXES = .y .l .H .h .m .C .cc .c .f .java .o
always_built_SUBDIRS =
SUBDIRS = $(always_built_SUBDIRS) .

GEN_LIBS = ../generic/libapbsgen.la
FEM_LIBS  = ../fem/libapbsfem.la 
PMG_LIBS  = ../mg/libapbsmg.la
BLAS_LIBS = @blas_lib@

libdir = ${prefix}/lib/${fetk_cpu_vendor_os}
lib_LTLIBRARIES = libapbs.la

libapbs_la_SOURCES = apbs_link.c
libapbs_la_LIBADD  = ${GEN_LIBS} ${FEM_LIBS} ${PMG_LIBS} ${BLAS_LIBS}
libapbs_la_LDFLAGS = -version-info ${APBS_VERSION}

INCLUDES         = -I${top_srcdir}/src/generic -I${top_srcdir}/src/mg -I${top_srcdir}/src/fem @inc_list@
AM_CFLAGS        = @wall_pedantic@
AM_LDFLAGS          = @LDFLAGS@ @FLIBS@

