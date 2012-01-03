/* src/aaa_inc/mccf.h.in.  Generated from configure.ac by autoheader.  */


/*
 * ***************************************************************************
 * MC = < Manifold Code >
 * Copyright (C) 1994-- Michael Holst
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 * 
 * rcsid="INTENTIONALLY LEFT BLANK"
 * ***************************************************************************
 */

/*
 * ***************************************************************************
 * File:     acconfig.h
 *
 * Purpose:  Generates the main configuration header "mccf.h" for MC.
 *
 * Author:   Michael Holst
 * ***************************************************************************
 */

#ifndef _MCCF_H_
#define _MCCF_H_


/* Does the AMD library exist? */
#undef HAVE_AMD

/* Does the ARPACK library exist? */
#undef HAVE_ARPACK

/* Does the BLAS library exist? */
#undef HAVE_BLAS

/* Does the CGCODE library exist? */
#undef HAVE_CGCODE

/* Am I running in a Cygwin/Win32 environment? */
#undef HAVE_CYGWIN

/* Do I compile as a debug version? */
#undef HAVE_DEBUG

/* Define to 1 if you have the <dlfcn.h> header file. */
#undef HAVE_DLFCN_H

/* Does EMBED macro for embedding rcsid symbols into binaries work? */
#define HAVE_EMBED

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H

/* Does the LAPACK library exist? */
#define HAVE_LAPACK

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H

/* Does the PMG library exist? */
#undef HAVE_PMG

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H

/* Does the SUPERLU library exist? */
#undef HAVE_SUPERLU

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H

/* Does the UMFPACK library exist? */
#undef HAVE_UMFPACK

/* Define to 1 if you have the <unistd.h> header file. */
#undef HAVE_UNISTD_H

/* Does the VF2C library exist? */
#undef HAVE_VF2C

/* Define to the sub-directory in which libtool stores uninstalled libraries.
   */
#undef LT_OBJDIR

/* Name of package */
#undef PACKAGE

/* Define to the address where bug reports for this package should be sent. */
#undef PACKAGE_BUGREPORT

/* Define to the full name of this package. */
#undef PACKAGE_NAME

/* Define to the full name and version of this package. */
#undef PACKAGE_STRING

/* Define to the one symbol short name of this package. */
#undef PACKAGE_TARNAME

/* Define to the version of this package. */
#undef PACKAGE_VERSION

/* Define to 1 if you have the ANSI C header files. */
#undef STDC_HEADERS

/* Version number of package */
#undef VERSION


/*
 * ***************************************************************************
 * Undo the damage done by the autoconf tools
 * ***************************************************************************
 */
#define HAVE_DLFCN_H
#define HAVE_INTTYPES_H
#define HAVE_MEMORY_H
#define HAVE_STDINT_H
#define HAVE_STDLIB_H
#define HAVE_STRINGS_H
#define HAVE_STRING_H
#define HAVE_SYS_STAT_H
#define HAVE_SYS_TYPES_H
#undef HAVE_UNISTD_H

/*  
 * ***************************************************************************
 * Define some RCS tag embedding and debug I/O macros
 * ***************************************************************************
 */

/* Embedded RCS tags ("ident filename" prints module versions in filename) */
#if defined(HAVE_EMBED)
#    define VEMBED(rctag) \
         VPRIVATE const char* rctag; \
         static void* use_rcsid=(0 ? &use_rcsid : (void*)&rcsid);
#else
#    define VEMBED(rctag)
#endif

/* Produce additional debugging I/O */
#if defined(HAVE_DEBUG)
#    define VDEBUGIO(str) fprintf(stderr,str)
#else
#    define VDEBUGIO(str)
#endif

#endif /* _MCCF_H_ */

