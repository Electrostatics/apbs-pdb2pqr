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

libdir = ${prefix}/lib/${fetk_cpu_vendor_os}
lib_LTLIBRARIES = libapbs.la

libapbs_la_SOURCES = apbs_link.c
libapbs_la_LIBADD  = ${GEN_LIBS} ${FEM_LIBS} ${PMG_LIBS}
libapbs_la_LDFLAGS = -version-info ${APBS_VERSION}

INCLUDES         = -I../generic -I../mg -I../fem @inc_list@
AM_CFLAGS        = @wall_pedantic@
LDFLAGS          = @LDFLAGS@ @FLIBS@

