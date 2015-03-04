/*
 * =====================================================================================
 *
 *       Filename:  hmm.c
 *
 *    Description:  Hidden Markov Model
 *
 *        Version:  1.0
 *        Created:  03/02/2015 16:40:50
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME Hyunkyu Lee 
 *   Organization:  Han Lab 
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int main(int argc, char **argv) {
	FILE *fin;   
   int i,j,h,k; //global variables that we may want to eliminate.
   char *gmapf; // genetic map data file name
   char *snpdataf; //reference snp data filename
	int nsnp; //total number of snps
	int nsamp; //total number of individuals in snp data
	int nhap; //total number of chromosomes (nsamp*2
	int *snpdata; //reference snp data.
	double *gmap; // genetic map data.
	int *target; //target genome
	int *ref; // reference data
	int nsampref; // total number of individuals in reference
	int N; // total number of chromosomes in reference (nsmapref * 2)
	char *tagf; // tag SNPs of a chipset
	
	int ntag; // number of tags
	int *is_tag; // is tag ?
	
	enum Mode {T} mode;
 
	if (!strcmp(argv[1], "-T")) { mode = T; }//support the first argument as mode selection.
	else { printf("only mode -T is available now");
		break;
	}
	i = 2; //since we already used argv[1]
	snpdataf = argv[i++]; //just file name we should call it by FILE function in c++ lib (stdio.h)
	gmapf = argv[i++];
	nsnp = atoi(argv[i++]);
	nsamp = atoi(argv[i++]);
   tagf = argv[i++];
	ntag = atoi(argv[i++]);
	nhap = 2*nsamp;

	snpdata = (int*)malloc(nsnp*nhap*sizeof(int));
	gmap = (double*)malloc(nsnp*sizeof(double));
	is_tag = (int*)malloc(nsnp*sizeof(int));
	
	// file read
	fin = fopen(snpdataf, "r");
	for (i = 0; i < nsnp; i++) { 
		for (h = 0; h < nhap; h++) {
		fscanf(fin, "%d", &snpdata[i*nhap+h]); //one at a time to address at the snpdata size
		}	
	}	
	fclose(fin);

	fin = fopen(gmapf, "r");
	for (i = 0; i < nsnp; i++) {
		fscanf(fin, "%lf", &gmap[i]);
		gmap[i] *= 0.01; // unit conversion factor (cM)
	}
	fclose(fin);

	if (mode == T) {
		fin = fopen(tagf, "r");
		for (i = 0; i < nsnp; i++)
			is_tag[i] = 0;
		for (u = 0; i < ntag; i++)
				fscanf(fin, "%d", &j);
				j--;
				is_tag[j] = 1;	
		}
	}
	nsampref = nsamp-1;
	N = nsampref*2;
	target = (int*)malloc(nsnp*sizeof(int));
	ref = (int*)malloc(nsnp*N*sizeof(int));
	for (i = 0; i < nsnp; i++) {
		target[i] = snpdata[i*nhap] + snpdata[i*nhap+1];
	}
	for (i = 0; i < nsnp; i++) {
		for (h = 0; h < N; h++) {
			ref[i*N*h] = snpdata[i*nhap+(h+2)];
		}
	}

	// HMM


		
}
