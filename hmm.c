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
	const double Ne = 11418;
	const double e = 0.01; //sequencing error rate
	double theta;
	double lambda;
	double **emission;
	double **alpha;
	double **beta;
	double trans_allsame;
	double trans_alldiff;
	double trans_onesame;
	int *arnrm;
	int *brnrm;
	double rho;
	double e_rhoN;
	double a, b, c, sum;
	double asum, bsum;
	double *ch;
	double *ck;
	double *beta_times_b;
	double *pstate;
	double *postp;
	double maxpstate;
	int maxh, maxk;
	int incorrect_allele_cnt;
	double accuracy;
	double allele_dosage_diff;
	const double BIG = pow(2,66);
	const double BIGI = pow(2,-66);
	
	theta = 0.;
	for (i = 1; i < N; i++) theta += 1./i;
	theta = 1./theta;
	lambda = theta / (2*(theta+N));
	emission = (double**)malloc(3*sizeof(double*));
	for (i = 0; i < 3; i++)
		emission[i] = (double*)malloc(3*sizeof(double));
	emission[0][0] = emission[2][2] = (1-lambda)*(1-lambda);
	emission[0][2] = emission[2][0] = lambda*lambda;
	emission[0][1] = emission[2][1] = 2*lambda*(1-lambda);
	emission[1][0] = emission[1][2] = lambda*(1-lambda);
	emission[1][1] = lambda*lambda+(1-lambda)*(1-lambda);

	alpha = (double**)malloc(nsnp*sizeof(double*));
	beta = (double**)malloc(nsnp*sizeof(double*));
	for (i = 0; i < nsnp; i++)	{
		alpha[i] = (double*)malloc(N*N*sizeof(double));
		beta[i] = (double*)malloc(N*N*sizeof(double));
	}
	pstate = (double*)malloc(N*N*sizeof(double));
	arnrm = (int*)malloc(nsnp*sizeof(int));
	brnrm = (int*)malloc(nsnp*sizeof(int));
	ch = (double*)malloc(N*sizeof(double));
	ck = (double*)malloc(N*sizeof(double));
	postp = (double*)malloc(3*sizeof(double));
	beta_times_b = (double*)malloc(N*N*sizeof(double));
	
	//forward probability	
	for (h = 0; h < N; h++) {
		for (k = 0; k < N; k++) { // this indicates for N*N states
			if (mode == T) {
				alpha[0][h*N+k] = (is_tag[0] == 1)? emission[ref[0*N+h]+ref[0*N+k]][target[0]] : 1.; // isn't observation tag? 
			}
		}
	}
	arnrm[0] = 0;
	for (i = 1; i < nsnp; i++) { // for all snps
		if (i % 1000 == 0) fprintf(stderr, "forward &d\n", i);
		rho = 4*Ne*(gmap[i]-gmap[i-1]);
		e_rhoN = exp(-rho/N);
		trans_allsame = pow(e_rhoN+(1-e_rhoN)/N, 2);
		trans_alldiff = pow((1-e_rhoN)/N, 2);
		trans_onesame = (e_rhoN+(1-e_rhoN)/N) * ((1-e_rhoN)/N);

// for efficiency, introduce ch = sum_k(alpha[i-1]), and similarly ck
// also, c = sum_{hk}(alpha[i-1])

		c = 0.;
		for (h = 0; h < N; h++) {
			ch[h] = ck[h] = 0.;
		}
		for (h = 0; h < N; h++) {
			for (k = 0; k < N; k++) {
				ch[h] += alpha[i-1][h*N+k];
				ck[k] += alpha[i-1][h*N+k];
				c += alpha[i-1][h*N+k];
			}	
		}
	
		// calculate alpha for each state
		asum = 0.;
		for (h = 0; h < N; h++) {
			for (k = 0; k < N; k++) { // for all state N*N
				if (mode == T) {
					b = (is_tag[i] == 1)? emission[ref[i*N+h]+ref[i*N+k]][target[i]] : 1.;  // isn't observation tag?
				}
				a = alpha[i-1][h*N+k];
				alpha[i][h*N+k] = ( (c-ch[h]-ck[k]+a)*trans_alldiff + (ch[h]+ck[k]-a)*trans_onesame + a*trans_allsame ) * b;
				asum += alpha[i][h*N+k];
			}	
		}
		arnrm[i] = arnrm[i-1];
		if (asum < BIGI) { // renormalize if necessary
			++arnrm[i];
			for (h = 0; h < N; h++) {
				for (k = 0; k < N; k++) { //for all N*N states
				alpha[i][h*N+k] *= BIG;
				}
			}
		}
		if (asum > BIG) { 
			--arnrm[i];
			for (h = 0; h < N; h++) {
				for (k = 0; k < N; k++) {
					alpha[i][h*N+k] *= BIGI;
				}
			}
		}
	}

	// backward probability
	for (h = 0; h < N; h++) {
		for (k = 0; k < N; k++) { //for all N*N
			beta[nsnp-1][h*N+k] = 1.;
		}
	}
	brnrm[nsnp-1] = 0;
	for (i = nsnp-2; i >= 0; i--) {
		if (i%1000 == 0) fprintf(stderr, "backward %d\n", i);
		rho = 4*Ne*(gmap[i+1]-gmap[i]);
		e_rhoN = exp(-rho/N);
		trans_allsame = pow(e_rhoN+(1-e_rhoN)/N, 2);
		trans_alldiff = pow((1-e_rhoN)/N, 2);
		trans_onesame = (e_rhoN+(1-e_rhoN)/N) * ((1-e_rhoN)/N);
		// introduce ch, ck, c
		// store beta*b first
		for (h = 0; h < N; h++) {
			for (k = 0; k < N; k++) {
				if (mode == T) {
					b = (is_tag[i+1] == 1)? emission[ref[(i+1)*N+h]+ref[(i+1)*N+k]][target[i+1]] : 1.; //isn't observation tag?
				}
				beta_times_b[h*N+k] = b * beta[i+1][h*N+k];
			}
		}
		c = 0.;
		for (h = 0; h < N; h++) {
			ch[h] = ck[h] = 0.;
		}
		for (h = 0; h < N; h++) {
			for (k = 0; k < N; k++) {
				ch[h] += beta_times_b[h*N+k];
				ck[k] += beta_times_b[h*N+k];
				c += beta_times_b[h*N+k];
			}
		}
		// calculate beta for each state
		bsum = 0.; 
		for (h = 0; h < N; h++) {
			for (k = 0; k < N; k++) { // for all previous N*N
				a = beta_times_b[h*N+k];
				beta[i][h*N+k] = (c-ch[h]-ck[k]+a)*trans_alldiff + (ch[h]+ck[k]-a)*trans)_onesame + a*trans_allsame;
				bsum += beta[i][h*N+k];
			}	
		}
		brnrm[i] = brnrm[i+1];
		if (bsum < BIGI) { // renormalize if necessary
			++brnrm[i];
			for (h = 0; h < N; h++) {
				for (k = 0; k < N; k++) { //for all N*N states
					beta[i][h*N+k] *= BIG;
				}
			}
		}
		if (bsum > BIG) {
			--brnrm[i];
			for (h = 0; h < N; h++) {
				for (k = 0; k < N; k++) { // for all N*N states 
					beta[i][h*N+k] *= BIGI;
				}
			}
		}
	}

	// pstate calculation
	incorrect_allele_cnt = 0;
	allele_dosage_diff = 0.;
	for (i = 0; i < nsnp; i++) {
		sum = 0.;
		for (h = 0; h < N; h++) {
			for (k = 0; k < N; k++) {
				sum += (pstate[h*N+k] = alpha[i][h*N+k]*beta[i][h*N+k]);
			}
		}
		for (h = 0; h < N; h++) {
			for (k = 0; k < N; k++) {
				pstate[h*N+k] /= sum;
			}
		}
		postp[0] = post[1] = post[2] = 0.;
		maxpstate = 0.;
		for (h = 0; h < N; h++) {
			for (k = 0; k < N; k++) {
				postp[ref[i*N+h]+ref[i*N+k]] += pstate[h*N+k];
				if (pstate[h*N+k] > maxpstate) {
				maxh = h;
				maxk = k;
				}
			}
		}
		if (postp[0] >= postp[1] && postp[0] >= postp[2]) {
			incorrect_allele_cnt += abs(target[i] - 0);
		} else if (postp[1] >= postp[2]) {
			incorrect_allele_cnt += abs(target[i] - 1);
		} else {
			incorrect_allele_cnt += abs(target[i] - 2);
		}
		allele_dosage_diff += fabs(postp[1]+postp[2]*2 - target[i]);
		
	}
	accuracy = (2.0*nsnp - incorrect_allele_cnt)/(2.0*nsnp);
	printf("%.5lf ", accuracy);
	accuracy = 1. - allele_dosage_diff/(2.0*nsnp);
	printf("%.5lf ", accuracy);
	printf("\n");
	return 0;
}


