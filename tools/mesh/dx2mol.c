
/* ///////////////////////////////////////////////////////////////////////////
// File:     dx2mol.c
//
// Purpose:  Convert OpenDx format potential to MolMol format
//
// Author:   Jung-Hsin Lin (bits added/modified by Nathan Baker)
//           Additional changes by Fred Damberger
//
// rcsid="$Id$"
/////////////////////////////////////////////////////////////////////////// */

#include "apbscfg.h"
#include "maloc/maloc.h"
#include "apbs/apbs.h"  

VEMBED(rcsid="$Id$")

int main(int argc, char **argv) {

    /* *************** VARIABLES ******************* */
    int u, i, j, k, nx, ny, nz;
    int n;
    double avg;
    double hy, hx, hzed, xmin, ymin, zmin;
    Vio *sock;
    Vgrid *grid;
    char *inpath = VNULL;
    char *outpath = VNULL;
    char *iodev = "FILE";
    char *iofmt = "ASC";
    char *thost = VNULL;
    char *usage = "\n\n\
    -----------------------------------------------------------------------\n\
    dx2mol (Contributed by Jung-Hsin Lin)\n\
    \n\
    For converting the OpenDX format of the electrostatic potential to the\n\
    MOLMOL format. MOLMOL is a popular free molecular display program\n\
    (http://www.mol.biol.ethz.ch/wuthrich/software/molmol/).\n\
    \n\
    Usage:  dx2mol file1.dx file2.pot\n\
            where file1.dx is a file in OpenDX format and file2.pot is the\n\
            file to be written in MOLMOL format.\n\
    -----------------------------------------------------------------------\n\
    \n";
 
 
    /* *************** CHECK INVOCATION ******************* */
    Vio_start();
    if (argc != 3) {
        Vnm_print(2, "\n*** Syntax error: got %d arguments, expected 3.\n\n",argc);
        Vnm_print(2,"%s\n", usage);
        return -1;
    } else {
        inpath = argv[1];
        outpath = argv[2];
    }

    /* Read DX format file */ 
    grid = Vgrid_ctor(0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, VNULL);
    Vgrid_readDX(grid, "FILE", "ASC", VNULL, inpath);
    nx = grid->nx;
    ny = grid->ny;
    nz = grid->nz;
    hx = grid->hx;
    hy = grid->hy;
    hzed = grid->hzed;
    xmin = grid->xmin;
    ymin = grid->ymin;
    zmin = grid->zmin;

    /* Intialize socket for writing */
    sock = Vio_ctor(iodev,iofmt,thost, outpath, "w");
    if (sock == VNULL) {
        Vnm_print(2, "Problem opening virtual socket %s\n",
          outpath);
        return -1;
    }
    if (Vio_connect(sock, 0) < 0) {
        Vnm_print(2, "Problem connecting virtual socket %s\n",
          outpath);
        return -1;
    }

    /* Write out data */
    Vio_printf(sock, "%8.3f %4d %5.3f\n",xmin,nx,hx); 
    Vio_printf(sock, "%8.3f %4d %5.3f\n",ymin,ny,hy); 
    Vio_printf(sock, "%8.3f %4d %5.3f\n",zmin,nz,hzed); 

    n = 0;
    for (k=0; k<nz; k++) {
        for (j=0; j<ny; j++) {
            for (i=0; i<nz; i++) {
                n++;
                u = nx*ny*k + nx*j + i;
                Vio_printf(sock, "%10.3e ", grid->data[u]);
                if (n == 10) {
                    Vio_printf(sock, "\n");
                    n = 0;
                }
            }
        }
    }

    /* Close off the socket */
    Vio_connectFree(sock);
    Vio_dtor(&sock);

    return 0;
}
