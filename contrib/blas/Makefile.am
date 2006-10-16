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
THISLIB = apbsblas

noinst_LTLIBRARIES = libapbsblas.la
libapbsblas_la_SOURCES  = mblasd.f

all-local:
	mkdir -p $(DESTDIR)$(libdir)
	cp .libs/libapbsblas.a $(DESTDIR)$(libdir)/libapbsblas.a; \
	ranlib $(DESTDIR)$(libdir)/libapbsblas.a;

