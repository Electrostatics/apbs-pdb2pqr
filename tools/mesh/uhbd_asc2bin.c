
#include "apbscfg.h"
#include "maloc/maloc.h"
#include "apbs/apbs.h"
#include "apbs/vhal.h"
#include "apbs/vmatrix.h"

#include <stdio.h>

#define MXGRD 129

int main(int argc, char **argv) {

    int im, jm, km, grdflg, one, kk;
    double h, ox, oy, oz;
    double scale;
    double dum2, dum3, dum4, dum5, dum6, dum7, dum8;
    int idum2, idum3, idum4;
    int i, j, k;

    char newfile[256];
    char flnm[256];
    char title[72];

    FILE *inFile = VNULL;
    FILE *outFile = VNULL;

    double *grid = (double*)malloc(MXGRD * MXGRD * MXGRD * sizeof(double));
    MAT3(grid, MXGRD, MXGRD, MXGRD);

    printf("Convert UHBD ascii grid to binary\n");
    printf("---------------------------------\n");

    printf("Enter grid map file name: ");
    scanf("%s", flnm);

    inFile = fopen(flnm, "r");
    VASSERT_MSG0(inFile != VNULL, "Couldn't open output file");



    printf("Enter output file name: ");
    scanf("%s", newfile);

    outFile = fopen(newfile, "wb");
    VASSERT_MSG0(outFile != VNULL, "Couldn't open input file");


    /* Read in the DX regular positions */
    /* Get "object" */
    fscanf(inFile, "%s", title);
    fwrite(title, 1, 72, outFile);

    fscanf(inFile, "%f", &scale);
    fwrite(&scale, 1, sizeof(double), outFile);

    fscanf(inFile, "%f", &dum2);
    fwrite(&dum2, 1, sizeof(double), outFile);

    fscanf(inFile, "%d", &grdflg);
    fwrite(&grdflg, 1, sizeof(int), outFile);

    fscanf(inFile, "%d", &idum2);
    fwrite(&idum2, 1, sizeof(int), outFile);

    fscanf(inFile, "%d", &km);
    fwrite(&km, 1, sizeof(int), outFile);

    fscanf(inFile, "%d", &one);
    fwrite(&one, 1, sizeof(int), outFile);

    fscanf(inFile, "%d", &km);
    fwrite(&km, 1, sizeof(int), outFile);

    fscanf(inFile, "%d", &im);
    fwrite(&im, 1, sizeof(int), outFile);

    fscanf(inFile, "%d", &jm);
    fwrite(&jm, 1, sizeof(int), outFile);

    fscanf(inFile, "%d", &km);
    fwrite(&km, 1, sizeof(int), outFile);

    fscanf(inFile, "%f", &h);
    fwrite(&h, 1, sizeof(double), outFile);

    fscanf(inFile, "%f", &ox);
    fwrite(&ox, 1, sizeof(double), outFile);

    fscanf(inFile, "%f", &oy);
    fwrite(&oy, 1, sizeof(double), outFile);

    fscanf(inFile, "%f", &oz);
    fwrite(&oz, 1, sizeof(double), outFile);

    fscanf(inFile, "%f", &dum3);
    fwrite(&dum3, 1, sizeof(double), outFile);

    fscanf(inFile, "%f", &dum4);
    fwrite(&dum4, 1, sizeof(double), outFile);

    fscanf(inFile, "%f", &dum5);
    fwrite(&dum5, 1, sizeof(double), outFile);

    fscanf(inFile, "%f", &dum6);
    fwrite(&dum6, 1, sizeof(double), outFile);

    fscanf(inFile, "%f", &dum7);
    fwrite(&dum7, 1, sizeof(double), outFile);

    fscanf(inFile, "%f", &dum8);
    fwrite(&dum8, 1, sizeof(double), outFile);

    fscanf(inFile, "%d", &idum3);
    fwrite(&idum3, 1, sizeof(int), outFile);

    fscanf(inFile, "%d", &idum4);
    fwrite(&idum4, 1, sizeof(int), outFile);

    for (k=0; k<km; k++) {

        fscanf(inFile, "%d", &kk);
        fwrite(&kk, 1, sizeof(int), outFile);

        fscanf(inFile, "%d", &im);
        fwrite(&im, 1, sizeof(int), outFile);

        fscanf(inFile, "%d", &jm);
        fwrite(&jm, 1, sizeof(int), outFile);

        for (j=0; j<jm; j++) {

            for (i=0; i<im; i++ ) {

                fscanf(inFile, "%f", RAT3(grid, i, j, k));
                fwrite(RAT3(grid, i, j, k), 1, sizeof(double), outFile);

            }
        }
    }

    return 0;
}