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

//#include <stdio.h>
#include <iostream>
#include <string.h>
#include <math.h>
#include <string.h>

int main(int argc, char **argv) {
   
   int i,j,h,k; //global variables that we may want to eliminate.
   char *gmapf; // genetic map data file name
   char *snpdataf; //reference snp data filename
	//
	int nsnp;
	int nsamp;

// exclude seqf 
	if (!strcmp(argv[1], "-T")) { mode = T; }//support the first argument as mode selection.
	else {cout << "only mode -T is available now"
	break;
	}
	i=2;
	snpdataf = argv[i++]; //is it saving filename? 
	gmapf = argv[i++];
	//
	nsnp = atoi(argv[i++]);
	nsamp = atoi(argv[i++]);
   tagf = argv[i++];
	ntag = atoi(argv[i++]);
	nhap = 2*nsamp;

	snpdata = (int*)malloc(nsnp*nhap*sizeof(int));

}








