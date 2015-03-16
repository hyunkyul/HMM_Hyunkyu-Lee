/*
* =====================================================================================
*
*			Filename:   example.cpp
*		Description: 	Beagle HMM model (node)
*		
*			 Version:	1.0
*			 Created:	03/ 16/ 2015	08:20:00
*			Revision:	none
*			Compiler:	gcc
*
*			  Author:	Hyunkyu Lee
*    Organization:	Han Lab
*
* ======================================================================================
*/			

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>

int main(int argc, char **argv)
{

	FILE *fin;
	char *nodecountf;
	char *edgepathf;

	int i,j,h,k; //universal variables
	int **nodecount; //node count
	int **edgepath; //edge path to next edge
	int D; //number of levels

	i=1;
	nodecountf= argv[i++];
	edgepathf= argv[i++];
	D=4;

	nodecount = (int**)malloc((D+1)*sizeof(int*));
	for (i=0;i<D+1;i++)
			nodecount[i] = (int*)malloc(pow(2, D)*sizeof(int));

	edgepath = (int**)malloc((D+1)*sizeof(int*));
	for (i=0;i<D+1;i++)
			edgepath[i] = (int*)malloc(pow(2, D)*sizeof(int));

// file read

	fin = fopen(nodecountf,"r");
	for (i=0;i<D+1;i++)
	{
		for (j=0;j<16; j++)
		{
		fscanf(fin, "%d", &nodecount[i][j]); //be at a tune to address at the nodecount size
		printf("%d ",nodecount[i][j]);
		}
		printf("\n");

   }
 
	fin = fopen(edgepathf,"r");
 	for (i=0;i<D+1;i++)
	{
		for (j=0;j<16; j++)
		{
		fscanf(fin, "%d", &edgepath[i][j]); //be at a time to address at the nodecount size
		printf("%d ",edgepath[i][j]);
		}
	printf("\n");
	}

//Beagle HMM!
	
	int *missing;
	missing = (int*)malloc(D*sizeof(int));

	printf("\nmissing i : ");
	for (i=0;i<D;i++) 
	{
		if(i==0) missing[i]=1;
		else missing[i]=0;
		printf("%d ",missing[i]);
	}
	printf("\n");


//forward algorithm  

	double **alpha;
	int **boolean;

   boolean = (int**)malloc(D*sizeof(int*));
	for (i=0;i<D;i++)
		boolean[i] = (int*)malloc(pow(2,D)*sizeof(int));

	alpha = (double**)malloc(D*sizeof(double*));
	for (i=0;i<D;i++)
		alpha[i] = (double*)malloc(pow(D,2)*sizeof(double));

	printf("\nboolean : \n"); //build boolean to calculate P(ei)*P(ei)
	for (i=0;i<D;i++)
	{
		for (j=0;j<pow(2,D);j++)
		{
			if(nodecount[i+1][j]!=0) boolean[i][j]=1;
			else boolean[i][j]=0;
			printf("%d",boolean[i][j]);
		}
		printf("\n");
	}	

	// calculate alpha
	
	for (i=0;i<D;i++)// for D*D states
	{	
		for(j=0;j<D;j++) alpha[0][D*i+j]= boolean[1][i]*boolean[1][j]*(nodecount[1][i]*nodecount[1][j]/pow(nodecount[0][0],2));
	}

	for (i=0;i<D;i++)
	{
		for (j=0;j<pow(D,2); j++)
		{
		printf("%.2f ",alpha[i][j]);
		}
	printf("\n");
	}



	for (i=1;i<D;i++)// for D*D states
	{	
		for (j=0;j<D;j++)
		{
			for(k=0;k<d;k++)
			{
				if (boolean[i][j]*boolean[i][k]==1) //only if both booleans are 1 
				{
					for (h=0;h<pow(D,2);h++) 
					{
						if(edgepath[j][k]==edgepath[i][h/4]&&edgepath[j][k]==edgepath[i][h%4]) 
							alpha[i][D*j+k]+=alpha[i-1][h]*nodecount[i][h%4]/nodecount[i-1][h/4 
					}
				}
			}
		}
	}

}






