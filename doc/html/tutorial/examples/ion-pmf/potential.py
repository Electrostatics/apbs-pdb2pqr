""" Vgrid: an OpenDX data manipulation class.

    A module to read and manipulate OpenDX data with strong similarities to 
    the APBS class with the same name. 
   
    Nathan A. Baker (baker@biochem.wustl.edu)
    Dept. of Biochemistry and Molecular Biophysics
    Center for Computational Biology
    Washington University in St. Louis
    Additional contributing authors listed in the code documentation.
    
    Copyright (c) 2003. Washington University in St. Louis
    All Rights Reserved.

    Portions copyright (c) 1999-2002.  University of California.
    Portions copyright (c) 1995.  Michael Holst.

    Permission to use, copy, modify, and distribute this software and its
    documentation for educational, research, and not-for-profit purposes,
    without fee and without a signed licensing agreement, is hereby granted,
    provided that the above copyright notice, this paragraph and the
    following two paragraphs appear in all copies, modifications, and
    distributions.

    IN NO EVENT SHALL THE AUTHORS BE LIABLE TO ANY PARTY FOR DIRECT,
    INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST
    PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION,
    EVEN IF THE AUTHORS HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH
    DAMAGE.

    THE AUTHORS SPECIFICALLY DISCLAIM ANY WARRANTIES, INCLUDING, BUT NOT
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE.  THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF
    ANY, PROVIDED HEREUNDER IS PROVIDED \"AS IS\".  THE AUTHORS HAVE NO
    OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
    MODIFICATIONS.
"""

import math, string 

class Vgrid:
    def __init__(self, dims, spac, origin, data, colMajor=1):
        """ Initialize the structure with:
            dims -- the number of grid points in each dimension 
                    (a 3-tuple of integers)
            spac -- the grid spacing (in A) in each dimension 
                    (a 3-tuple of floats)
            origin -- the coordinates (in A) of the grid lower corner
                      (a 3-tuple of floats)
            data -- array of data with x increasing most quickly
                    (an array of floats)
            colMajor -- flag indicating whether the data is in column major (z
                        increases most quickly) or row major (x increases most
                        quickly) format.  If colMajor == 1, then we assume
                        column-major; if colMajor = 0, then we assume row-major
        """
        self.dims = dims
        self.spac = spac
        self.origin = origin
        self.data = data
        self.colMajor = colMajor

    def ijk2xyz(self, index):
        """ Transform a data array index tuple into coordinate tuple
            index -- the index of interest (3-tuple of integers) 
        """
        coord = []
        for i in range(0, 3):
            ii = index[i]
            if (ii >= self.dims[i]) or (i < 0):
                errstr = "Index element %d (%d) out of range!\n" % (i, ii)
                raise IndexError, errstr
            x = float(ii)*self.spac[i] + self.origin[i]
            coord.append(x)
        return (coord[0], coord[1], coord[2])        
        
    def ijk2u(self, index):
        """ Transform a data array index tuple into the natural ordering index
            of the data array
            index -- the index of interest (3-tuple of integers) 
        """
        for i in range(0, 3):
            ii = index[i]
            if (ii >= self.dims[i]) or (i < 0):
                errstr = "Index element %d (%d) out of range!\n" % (i, ii)
                raise IndexError, errstr
        if self.colMajor:
            u = (index[2]*self.dims[0]*self.dims[1]) + (index[1]*self.dims[0]) + index[0]
        else:
            u = (index[0]*self.dims[2]*self.dims[1]) + (index[1]*self.dims[2]) + index[2]
        return u

    def value(self, pt):
        """ Get the data value at a specified point. 
            pt -- the coordinates of the point of interest (3-tuple of floats)
        """
        nx = self.dims[0]
        ny = self.dims[1]
        nz = self.dims[2]
        hx = self.spac[0]
        hy = self.spac[1]
        hz = self.spac[2]
        xmin = self.origin[0]
        ymin = self.origin[1]
        zmin = self.origin[2]
        xmax = xmin+ float(nx-1)*hx
        ymax = ymin+ float(ny-1)*hy
        zmax = zmin+ float(nz-1)*hz

        ifloat = (pt[0] - xmin)/hx
        jfloat = (pt[1] - ymin)/hy
        kfloat = (pt[2] - zmin)/hz
        # print
        # print ifloat, jfloat, kfloat
        # print 
        ihi = int(math.ceil(ifloat))
        jhi = int(math.ceil(jfloat))
        khi = int(math.ceil(kfloat))
        ilo = int(math.floor(ifloat))
        jlo = int(math.floor(jfloat))
        klo = int(math.floor(kfloat))
        if math.fabs(pt[0] - xmin) < 0.0:  ilo = 0
        if math.fabs(pt[1] - ymin) < 0.0:  jlo = 0
        if math.fabs(pt[2] - zmin) < 0.0:  klo = 0
        if math.fabs(pt[0] - xmax) < 0.0:  ihi = nx - 1
        if math.fabs(pt[1] - ymax) < 0.0:  jhi = ny - 1
        if math.fabs(pt[2] - zmax) < 0.0:  khi = nz - 1

        dx = ifloat - float(ilo);
        dy = jfloat - float(jlo);
        dz = kfloat - float(klo);

        u =      dx *     dy *     dx *self.data[self.ijk2u((ihi,jhi,khi))] \
          +      dx *(1.0-dy)*     dx *self.data[self.ijk2u((ihi,jlo,khi))] \
          +      dx *     dy *(1.0-dx)*self.data[self.ijk2u((ihi,jhi,klo))] \
          +      dx *(1.0-dy)*(1.0-dx)*self.data[self.ijk2u((ihi,jlo,klo))] \
          + (1.0-dx)*     dy *     dx *self.data[self.ijk2u((ilo,jhi,khi))] \
          + (1.0-dx)*(1.0-dy)*     dx *self.data[self.ijk2u((ilo,jlo,khi))] \
          + (1.0-dx)*     dy *(1.0-dx)*self.data[self.ijk2u((ilo,jhi,klo))] \
          + (1.0-dx)*(1.0-dy)*(1.0-dx)*self.data[self.ijk2u((ilo,jlo,klo))]

        return u

    def integrate(self, vmin, vmax, nquad):
        """ Integrate over a rectangular volume with the specified spacing:
            vmin -- coordinates of volume lower corner (3-tuple of floats)
            vmax -- coordinates of volume upper corner (3-tuple of floats)
            nquad -- number of qudrature points (3-tuple of ints)
        """
        xlen = vmax[0] - vmin[0]
        ylen = vmax[1] - vmin[1]
        zlen = vmax[2] - vmin[2]
        nx = nquad[0]
        ny = nquad[1]
        nz = nquad[2]
        hx = xlen/float(nx-1)
        hy = ylen/float(ny-1)
        hz = zlen/float(nz-1)
        xmin = self.origin[0]
        ymin = self.origin[1]
        zmin = self.origin[2]

        u = 0.0
        for k in range(0, nz):
            if (k == 0) or (k == (nz-1)):  wz = 0.5
            else:  wz = 1.0
            z = float(k)*hz + zmin
            for j in range(0, ny):
                if (j == 0) or (j == (ny-1)):  wy = 0.5
                else:  wy = 1.0
                y = float(j)*hy + ymin
                for i in range(0, nx):
                    if (i == 0) or (i == (nx-1)):  wx = 0.5
                    else:  wx = 1.0
                    x = float(i)*hx + xmin
                    u = u + wx*wy*wz*self.value((x, y, z))
        u = hx*hy*hz*u
        return u

    def readOpenDX(self, file):
        """ Read OpenDX-format data from the open file object, replacing all
            existing internal data """
        data = []
        nx = None
        ny = None
        nz = None
        hx = None
        hy = None
        hz = None
        ox = None
        oy = None
        oz = None
        data = []
        while (1):
            line = file.readline()
            if line == "":  break
            words = string.split(string.strip(line))
            if len(words[0]) > 0:
                if words[0] != "#":
                    if words[0] == "object":
                        if words[3] == "gridpositions":
                            nx = int(words[5])
                            ny = int(words[6])
                            nz = int(words[7])
                    elif words[0] == "delta":  
                        dx = float(words[1])
                        dy = float(words[2])
                        dz = float(words[3])
                        if dx != 0.0:  
                            hx = dx
                            if hz == None:  colMajor = 0
                        elif dy != 0.0:  hy = dy
                        elif dz != 0.0:  
                            hz = dz
                            if hx == None:  colMajor = 1
                    elif words[0] == "origin":  
                        ox = float(words[1])
                        oy = float(words[2])
                        oz = float(words[3])
                    else: 
                        try:
                            for word in words: data.append(float(word))
                        except ValueError, details:  pass
        if len(data) != (nx*ny*nz):
            errstr = "Read %d values, expected %d!" % (len(data), (nx*ny*nz))
        self.dims = (nx, ny, nz)
        self.spac = (hx, hy, hz)
        self.origin = (ox, oy, oz)
        self.data = data
        self.colMajor = colMajor

def main():
    """ A main driver routine, mainly for testing """

    import sys, random
 
    # Read or fabricate data
    if len(sys.argv) == 2:
        path = sys.argv[1]
        file = open(path, "r")
        dims = (None, None, None)
        spac = (None, None, None)
        origin = (None, None, None)
        n = None
        data = []
        vgrid = Vgrid(dims, spac, origin, data)
        vgrid.readOpenDX(file)
        file.close()
    else:
        dims = (4, 4, 4)
        spac = (0.5, 0.5, 0.5)
        origin = (0.5, 0.5, 0.5)
        n = dims[0]*dims[1]*dims[2]
        data = []
        for i in range(0, n):  data.append(1.0)
        vgrid = Vgrid(dims, spac, origin, data)
    # Derived
    nx = vgrid.dims[0]
    ny = vgrid.dims[1]
    nz = vgrid.dims[2]
    hx = vgrid.spac[0]
    hy = vgrid.spac[1]
    hz = vgrid.spac[2]
    xmin = vgrid.origin[0]
    ymin = vgrid.origin[1]
    zmin = vgrid.origin[2]
    xlen = float(nx-1)*hx
    ylen = float(ny-1)*hy
    zlen = float(nz-1)*hz
    xmax = xmin + xlen
    ymax = ymin + ylen
    zmax = zmin + zlen
    xcent = xmin + 0.5*xlen
    ycent = ymin + 0.5*ylen
    zcent = zmin + 0.5*zlen
    for i in range(1, (nx-1)):
        x = xmin + float(i)*hx
        y = ycent
        z = zcent 
        print x, vgrid.value((x, y, z))

if __name__ == "__main__":  main()
    
