
SUFFIXES = .H .h .C .cc .c .F .f .o

SUBDIRS = .

F77COMPILE = $(F77) $(DEFS) $(FPP_FLAGS) $(AM_FFLAGS) $(FFLAGS)

libdir = ${prefix}/lib


lib_LIBRARIES = libiapbs.a
libiapbs_a_SOURCES = apbs_driver.c apbs_driver.h

includedir = $(prefix)/include/apbs
include_HEADERS = apbs_driver.h

AM_CPPFLAGS = $(apbs_inc_path)

EXTRA_PROGRAMS = wrapper
wrapper_SOURCES = wrapper.f apbs_driver.c apbs_driver.h

wrapper_LDADD = libiapbs.a
wrapper_LDFLAGS = $(apbs_lib_path) $(apbs_libs) $(apbs_blas) $(z_lib)


