## ###########################################################################
## File:    Makefile.am
##
## Purpose: Automatically generate Makefile.in for autoconf using automake.
##
## Author:  David Gohara
## ###########################################################################

AUTOMAKE_OPTIONS = no-dependencies
SUFFIXES = .y .l .H .h .m .C .cc .c .f .java .o
always_built_SUBDIRS =
SUBDIRS = $(always_built_SUBDIRS) .
THISLIB = aqua

noinst_LTLIBRARIES = libapbsaqua.la
libapbsaqua_la_SOURCES = Blaswrapper.f Cartesian_grid.f \
  CGMG_Solver.f Galerkin.f LinMG.f MatVect.f Newton_Solver_pbe.f \
  Smoothing.f PrintInfo.f
INCLUDES           = -I${top_srcdir}/src/generic -I${top_srcdir}/src/fem @inc_list@
AM_CFLAGS          = @wall_pedantic_ansi@ @CFLAGS@


all-local:
#	-cp -p ${top_srcdir}/src/${THISLIB}/apbs/*.h ${top_srcdir}/src/aaa_inc/apbs/.

