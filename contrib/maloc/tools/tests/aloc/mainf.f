c* ***************************************************************************
c* MALOC = < Minimal Abstraction Layer for Object-oriented C >
c* Copyright (C) 1994--2008 Michael Holst
c*
c* This library is free software; you can redistribute it and/or
c* modify it under the terms of the GNU Lesser General Public
c* License as published by the Free Software Foundation; either
c* version 2.1 of the License, or (at your option) any later version.
c*
c* This library is distributed in the hope that it will be useful,
c* but WITHOUT ANY WARRANTY; without even the implied warranty of
c* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
c* Lesser General Public License for more details.
c*
c* You should have received a copy of the GNU Lesser General Public
c* License along with this library; if not, write to the Free Software
c* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
c*
c* rcsid="$Id: mainf.f,v 1.9 2008/03/12 05:14:00 fetk Exp $"
c* ***************************************************************************

c*****************************************************************************
c* FORTRAN interface to vnm library:
c*
c*     int function vrand();
c*     int function vrndmx();
c*     double function vepsma();
c*     subroutine vtstrt(timer, name, len);
c*     subroutine vtstop(timer, name, len);
c*     int function vsystm(cmd, len);
c*     subroutine vnmprt(int *unit, char *strng, int *len);
c*
c*         timer  --- integer timer number
c*         name   --- char timer name of length eq len (with no '\0')
c*         cmd    --- char command of length eq len (with no '\0')
c*         val    --- integer tag number
c*         unit   --- integer I/O unit number
c*         strng  --- char message of length eq len (with no '\0')
c*
c* FORTRAN interface to vio library:
c*
c*     int function vioctr(type, frmt, host, lenh, file, lenf, cmode)
c*     int function viodtr(sockn)
c*     int function vioutl(sockn, mode)
c*     subroutine vioint(sockn, ival, len)
c*     subroutine vioflt(sockn, rval, len)
c*     subroutine viodbl(sockn, dval, len)
c*     subroutine viostr(sockn, sval, len)
c*
c*         sockn  --- integer socket number
c*
c*         type   --- char*4 vio device type (FILE,PIPE,UNIX,INET,MPI1)
c*         frmt   --- char*3 datatype (XDR,ASC)
c*
c*         host   --- char hostname of length eq lenh (with no '\0')
c*         file   --- char filename of length eq lenf (with no '\0')
c*
c*         cmode  --- char*1 'r' create for read, or 'w' create for write
c*         mode   --- char*1 'o' open, or 'c' close (both for read/write)
c*
c*         ival   --- int array of length len
c*         rval   --- float array of length len
c*         dval   --- double array of length len
c*         sval   --- char array of length len (no terminating '\0')
c*
c*****************************************************************************

      program main
c* ***************************************************************************
c* File:     mainf.f
c*
c* Purpose:  Test the FORTRAN interface to the VNM and VIO libraries.
c*
c* Author:   Michael Holst
c* ***************************************************************************
      implicit none
      integer vio1, vio2, rc
      integer vioctr, viodtr, vioutl
      integer vrand, vrndmx
      double precision vepsma

      call vtstrt(10, "Test Program", 12);

c*    *** some VNM testing ***
      call vnmprt(0, '0: This message is being sent to VNM unit 0',43)
      call vnmprt(1, '1: This message is being sent to VNM unit 1',43)
      call vnmprt(2, '2: This message is being sent to VNM unit 2',43)
      call vnmpri(2, '2: A call to vrand()  produces:   ',34,vrand())
      call vnmpri(2, '2: A call to vrndmx() produces:   ',34,vrndmx())
      call vnmprd(2, '2: A call to vepsma() produces:   ',34,vepsma())
      call vnmprt(2, '2: Sending "df" command to system ',34)
      call vsystm('df', 2)
      call vnmprt(2, '2: Returned from the "df" command ',34)

c*    *** start (initialize) the VIO library ***
      call viosta()

c*    *** you can call the VIO initializer multiple times w/o error ***
      call viosta()
      call viosta()
      call viosta()

c*    *** construct VIO1 (socket 1) and VIO2 (socket 2) ***
      vio1 = vioctr('INET', 'ASC', 'localhost',9, '1',1, 'w')
      if (vio1 .lt. 0) call vnmprt(2,'2: vioctr vio1 failed',21)
      vio2 = vioctr('INET', 'ASC', 'localhost',9, '2',1, 'w')
      if (vio2 .lt. 0) call vnmprt(2,'2: vioctr vio2 failed',21)

c*    *** open VIO1 and VIO2 connections for writing ***
      rc = vioutl(vio1,'o')
      if (rc .ne. 0) call vnmprt(2,'2: vioutl vio1 failed',21)
      rc = vioutl(vio2,'o')
      if (rc .ne. 0) call vnmprt(2,'2: vioutl vio2 failed',21)

c*    *** write things to VIO1 and VIO2 ***
      call myoutput(vio1)
      call myoutput(vio2)

c*    *** close VIO1 and VIO2 connections ***
      rc = vioutl(vio1,'c')
      if (rc .ne. 0) call vnmprt(2,'2: vioutl vio1 failed',21)
      rc = vioutl(vio2,'c')
      if (rc .ne. 0) call vnmprt(2,'2: vioutl vio2 failed',21)

c*    *** destroy VIO1 and VIO2 ***
      rc = viodtr(vio1)
      if (rc .ne. 0) call vnmprt(2,'2: viodtr vio1 failed',21)
      rc = viodtr(vio2)
      if (rc .ne. 0) call vnmprt(2,'2: viodtr vio2 failed',21)

c*    *** stop (de-initialize) the VIO library ***
      call viostp()

c*    *** you can call the VIO de-initializer multiple times w/o error ***
      call viostp()
      call viostp()
      call viostp()

c*    ************************************************************************
c*    *** NOTE: once you call "viostp()" one time, you MUST call "viosta()"
c*    ***       BEFORE you can use any VIO routines again...otherwise an
c*    ***       error message is generated and you don't get any I/O.
c*    ***       I.e., any calls to the VIO routine must be preceeded by
c*    ***       at least one call to "viosta()", with no intervening calls
c*    ***       to "viostp()".
c*    ************************************************************************

      call vtstop(10, "Test Program", 12);
      call vnmprt(2,'2: Stopping the program.                   ',43)
      stop
      end
c* ***************************************************************************
      subroutine myoutput(vio)
      integer vio

c*    *** write out a single red triangle in BH format ***
      call viostr(vio,'bhsingle',8)

      call viostr(vio,'putl',4)
      call vioint(vio,1,1)
      call vioflt(vio,1.0,1)
      call vioflt(vio,0.0,1)
      call vioflt(vio,0.0,1)

      call viostr(vio,'list',4)
      call vioint(vio,7,1)
      call viostr(vio,'fill',4)
      call vioint(vio,1,1)
      call vioint(vio,3,1)
      call vioflt(vio,0.1,1)
      call vioflt(vio,0.9,1)
      call vioflt(vio,0.1,1)
      call vioflt(vio,0.1,1)
      call vioflt(vio,0.1,1)
      call vioflt(vio,0.9,1)
      call vioflt(vio,0.1,1)
      call vioflt(vio,0.1,1)
      call vioflt(vio,0.1,1)
      call viostr(vio,'list',4)
      call vioint(vio,-7,1)

      call viostr(vio,'putl',4)
      call vioint(vio,-1,1)

      return
      end

