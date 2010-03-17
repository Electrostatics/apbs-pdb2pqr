/*
 *  del2dx.c
 *  apbs
 *
 *  Created by David Gohara on 3/17/10.
 *  Copyright 2010 Washington University. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

char *usage = "\n\n\
-----------------------------------------------------------------------\n\
del2dx\n\
\n\
For converting the DelPhi format of the electrostatic potential to the\n\
OpenDX format.\n\
\n\
Usage:  del2dx delphi_file opendx_file\n\
where delphi_file is a file in DelPhi format and opendx_file is the\n\
file to be written in OpenDX format.\n\
-----------------------------------------------------------------------\n\
\n";

int main(int argc, char **argv) {
	
	int i,j,k,index;
	char *inpath = NULL;
    char *outpath = NULL;
	
	char buffer[1024];
	
	int igrid, tot_grid;
	float val, scale, oldmid[3];
	
	FILE * pfile = NULL;
	
    if (argc != 3) {
        printf("\n*** Syntax error: got %d arguments, expected 3.\n\n",argc);
        printf("%s\n", usage);
        return -1;
    } else {
        inpath = argv[1];
        outpath = argv[2];
    }
	
	pfile = fopen(inpath, "r+b");
	
	//Get the grid length (all dimensions are the same)
	fseek(pfile, sizeof(int) * -1, SEEK_END);
	fread(&igrid, 1, sizeof(int), pfile);
	tot_grid = igrid*igrid*igrid;
	
	float * data = (float *)calloc(tot_grid, sizeof(float));
	
	rewind(pfile);
	
	//Skip the first 20 characters, they simply say
	// "now starting phimap "
	fseek(pfile, sizeof(char) * 20, SEEK_CUR);
	
	//get the map file time: potential or concentrat
	memset(buffer, 0, 1024);
	fread(buffer, 1, sizeof(char) * 10, pfile);
	printf("DelPhi %s map grid dimensions are: %i %i %i\n",
		   buffer,igrid,igrid,igrid);
	
	//Now skip ahead 60 characters to the data set.
	fseek(pfile, sizeof(char) * 60, SEEK_CUR);
	
	//Read in the actual data (float)
	fread(data, tot_grid, sizeof(float), pfile);
	
	//Skip over the text " end of phimap "
	fseek(pfile, sizeof(char) * 16, SEEK_CUR);
	
	//Get the scale
	fread(&scale, 1, sizeof(float), pfile);
	
	//Get the old midpoints?
	fread(oldmid, 3, sizeof(float), pfile);
	
	printf("The DelPhi %s map scale is: %f\n" \
		   "Old midpoint values: %f %f %f\n",
		   buffer,scale,oldmid[0],oldmid[1],oldmid[2]);
	
	fclose(pfile);
	
	//Now calculate some numbers needed for writing the DX header
	float temp = 0.;
	float xmax = 0.;
	for(i=0;i<3;i++){
		temp = fabsf(oldmid[i]);
		xmax = (xmax < temp) ? temp : xmax;
	}
	
	float range = ((float)igrid - 1.) / (2. * scale);
	float extent = range + xmax;
	float origin[3], delta[3];
	
	for(i=0;i<3;i++) origin[i] = oldmid[i] - range;
	for(i=0;i<3;i++) delta[i] = (range * 2.) / (float)igrid;
	
	//Now write the OpenDX formatted file
	pfile = fopen(outpath, "w");
	
	fprintf(pfile,
		   "# Data from del2dx (APBS 1.2.1)\n"	\
		   "# \n"							\
		   "# POTENTIAL (kT/e)\n"			\
		   "# \n"							\
		   "object 1 class gridpositions counts %i %i %i\n"	\
		   "origin %1.6e %1.6e %1.6e\n"	\
		   "delta %1.6e 0.000000e+00 0.000000e+00\n"		\
		   "delta 0.000000e+00 %1.6e 0.000000e+00\n"		\
		   "delta 0.000000e+00 0.000000e+00 %1.6e\n"		\
		   "object 2 class gridconnections counts %i %i %i\n"\
		   "object 3 class array type double rank 0 items %i data follows\n",
			igrid,igrid,igrid,
			origin[0],origin[1],origin[2],
			delta[0],delta[1],delta[2],
			igrid,igrid,igrid,tot_grid);
	
	//For the moment I'm assuming the data for DelPhi is row major
	//Write out the data
	int icol = 0;
	for (i=0; i<igrid; i++) {
		for (j=0; j<igrid; j++) { 
			for (k=0; k<igrid; k++) {
				index = k*(igrid)*(igrid)+j*(igrid)+i;
				fprintf(pfile, "%12.6e ", data[index]);
				icol++;
				if (icol == 3) {
					icol = 0;
					fprintf(pfile, "\n");
				}
			}
		}
	}
	if (icol != 0) fprintf(pfile, "\n");
	fprintf(pfile, "attribute \"dep\" string \"positions\"\n" \
			"object \"regular positions regular connections\" class field\n" \
			"component \"positions\" value 1\n" \
			"component \"connections\" value 2\n" \
			"component \"data\" value 3\n");
	
	free(data);
	
	return 0;
}