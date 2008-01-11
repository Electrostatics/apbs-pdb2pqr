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

libdir = ${contrib_path}/lib/${fetk_cpu_vendor_os}
lib_LTLIBRARIES = libapbsblas.la
libapbsblas_la_SOURCES = mblasd.f
libapbsblas_la_LDFLAGS = -version-info ${APBS_VERSION}

#all-local:
#	mkdir -p $(DESTDIR)$(libdir)
#	cp .libs/libapbsblas.a $(DESTDIR)$(libdir)/libapbsblas.a; \
#	ranlib $(DESTDIR)$(libdir)/libapbsblas.a;

