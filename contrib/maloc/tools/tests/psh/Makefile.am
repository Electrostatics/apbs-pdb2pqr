### ##########################################################################
### LAF = < Linear Algebra Framework >
### Copyright (C) 1994--2006  Michael Holst
###
### This program is free software; you can redistribute it and/or modify it
### under the terms of the GNU General Public License as published by the
### Free Software Foundation; either version 2 of the License, or (at your
### option) any later version.
###
### This program is distributed in the hope that it will be useful,
### but WITHOUT ANY WARRANTY; without even the implied warranty of
### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
### See the GNU General Public License for more details.
###
### You should have received a copy of the GNU General Public License along
### with this program; if not, write to the Free Software Foundation, Inc.,
### 675 Mass Ave, Cambridge, MA 02139, USA.
###
### rcsid="$Id: Makefile.am,v 1.25 2006/06/03 07:22:30 mholst Exp $"
### ##########################################################################

## ###########################################################################
## File:    Makefile.am
##
## Purpose: Automatically generate Makefile.in for autoconf using automake.
##
## Notes:   The final makefile produced by autoconf will use the following:
##
##     COMPILE = $(CC) $(DEFS) $(INCLUDES) $(CPPFLAGS) $(CFLAGS)
##     LINK = $(CC) $(CFLAGS) $(LDFLAGS) -o $@
##
## Author:  Michael Holst
## ###########################################################################

AUTOMAKE_OPTIONS = no-dependencies
SUFFIXES = .y .l .H .h .m .C .cc .c .f .java .o
always_built_SUBDIRS =
SUBDIRS = $(always_built_SUBDIRS) .

noinst_PROGRAMS  = govmpi govcom

govmpi_SOURCES = vmpitst.c
govmpi_LDADD   = @maloc_lib@ @mpi_lib@ @rl_lib@ @liberty_lib@ @socket_lib@ @nsl_lib@ @thread_lib@ @math_lib@ @zfence_lib@ @profile_lib@

govcom_SOURCES = vcomtst.c
govcom_LDADD   = @maloc_lib@ @mpi_lib@ @rl_lib@ @liberty_lib@ @socket_lib@ @nsl_lib@ @thread_lib@ @math_lib@ @zfence_lib@ @profile_lib@

INCLUDES       = @maloc_inc@
AM_CFLAGS      = @profile@ @wall_pedantic_ansi@

