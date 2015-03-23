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
	char *samplef;

	int i,j,h,k; //universal variables
	int **sample; //sample data
	int **nodecount; //node count
	int **edgepath; //edge path to next edge
	int D; //number of levels

	i=1;
	samplef= argv[i++];
	nodecountf= argv[i++];
	edgepathf= argv[i++];
	D=4;

// allocate memories on the matrices

	sample = (int**)malloc(D*sizeof(int*));
	for (i=0;i<D;i++)
			sample[i] = (int*)malloc(2*sizeof(int));

	nodecount = (int**)malloc((D+1)*sizeof(int*));
	for (i=0;i<D+1;i++)
			nodecount[i] = (int*)malloc(pow(2, D)*sizeof(int));

	edgepath = (int**)malloc((D+1)*sizeof(int*));
	for (i=0;i<D+1;i++)
			edgepath[i] = (int*)malloc(pow(2, D)*sizeof(int));

// file read
	
	printf("nodecount matrix\n\n");
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
 	
	printf("\n\nedgepath matrix\n\n");
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
	
	printf("\n\nsample status\n\n");
	fin = fopen(samplef,"r");
	for (i=0;i<2;i++)
	{
		for (j=0;j<D;j++)
		{
		fscanf(fin, "%d", &sample[i][j]); // be at a time to address at the nodecount size
		printf("%d ",sample[i][j]);
		}
	printf("\n");
	}

//Beagle HMM!

	int **opennode;
	double **alpha;
	int *edgenum;


	printf("\n\n edgenumbers[i]\n\n");	
	edgenum = (int*)malloc(D*sizeof(int*));
	for (i=1;i<D+1;i++) //nodecount[][] has fake edge 0
	{	
		for (j=pow(2,D);j>0;j--)
		{
			if (nodecount[i][j-1]!=0)	
			{		
			edgenum[i-1]=j; 
			j=0;
			}
		}
		printf("%d",edgenum[i-1]);
	}
	printf("\n");

	opennode = (int**)malloc(D*sizeof(int*));  			//matrix size of [D][edgenum[i]^2]
	for (i=0;i<D;i++)
			opennode[i] = (int*)malloc(pow(edgenum[i],2)*sizeof(int));
	
	alpha = (double**)malloc(D*sizeof(double**));		//matrix size of [D][edgenum[i]^2]
	for (i=0;i<D;i++)
			alpha[i] = (double*)malloc(pow(edgenum[i],2)*sizeof(double));

	// open node calculation : calculating paths which we are interested	and vanishing nodes which are not useful




	for (i=0;i<D;i++)
	{
		for (j=0;j<edgenum[i];j++)
		{
			for (k=0;k<edgenum[i];k++)
			{
				if (sample[i][0] == 0 || sample [i][1] == 0) 
				{
					if(sample[i][0] == 0 && sample[i][1] != 0)		opennode[i][edgenum[i]*j+k] = ((k%2)+1 == sample[i][1] && nodecount[i+1][k] != 0 )? 1: 0;
					else if(sample[i][0] != 0 && sample[i][1] == 0)	opennode[i][edgenum[i]*j+k] = ((j%2)+1 == sample[i][0] && nodecount[i+1][j] != 0 )? 1: 0;
					else 															opennode[i][edgenum[i]*j+k] = 1;
				} // if sample has a missing node for given i phase ::: and remained opennodes are about {1,2}


				else if (sample[i][0] == 1 || sample[i][1] == 1)
				{
					if(sample[i][0] == 1 && sample[i][1] != 1) 		opennode[i][edgenum[i]*j+k] = ((j%2)+1 == sample[i][0] && nodecount[i+1][j] != 0)? 1: 0;
					else if(sample[i][0] != 1 && sample[i][1] == 1) opennode[i][edgenum[i]*j+k] = ((k%2)+1 == sample[i][1] && nodecount[i+1][k] != 0)? 1: 0;
					else 	opennode[i][edgenum[i]*j+k] = ((j%2)+1 == sample[i][0] && (k%2)+1 == sample[i][1] && nodecount[i+1][k] != 0 && nodecount[i+1][j] != 0)? 1: 0;
				}
				else if (sample[i][0] == 2 || sample[i][1] == 2)
				{
					if(sample[i][0] == 2 && sample[i][1] != 2) 		opennode[i][edgenum[i]*j+k] = ((j%2)+1 == sample[i][0] && nodecount[i+1][j] != 0)? 1: 0;
					else if(sample[i][0] != 2 && sample[i][1] == 2) opennode[i][edgenum[i]*j+k] = ((k%2)+1 == sample[i][1] && nodecount[i+1][k] != 0)? 1: 0;
					else 	opennode[i][edgenum[i]*j+k] = ((j%2)+1 == sample[i][0] && (k%2)+1 == sample[i][1] && nodecount[i+1][k] != 0 && nodecount[i+1][j] != 0)? 1: 0;
				}
			printf("%d	",opennode[i][edgenum[i]*j+k]);
	 		}
		}

	printf("\n");
	}
	
		


//forward algorithm  








/*
	double **alpha;
	int **boolean;
   int edgenum[4]={2,3,3,5};


   boolean = (int**)malloc(D*sizeof(int*));
	for (i=0;i<D;i++)
		boolean[i] = (int*)malloc(pow(2,D)*sizeof(int));

	alpha = (double**)malloc(D*sizeof(double*));
	for (i=0;i<D;i++)
		alpha[i] = (double*)malloc(pow(D,2)*sizeof(double));

	printf("\nboolean : \n"); //build boolean to calculate P(ei)*P(ei)
	for (i=0;i<D;i++)
	{
		for (j=0;j<edgenum[i];j++)
		{
			if(nodecount[i+1][j]!=0) boolean[i][j]=1;
			else boolean[i][j]=0;
			printf("%d",boolean[i][j]);
		}
		printf("\n");
	}	

	// calculate alpha
	
	for (i=0;i<edgenum[0];i++)// for D*D states
	{	
		for(j=0;j<edgenum[i];j++)
		{
			alpha[0][D*i+j]= boolean[1][i]*boolean[1][j]*(nodecount[1][(edgepath[0][i]-1)*2]*nodecount[1][(edgepath[0][j]-1)*2]/pow(nodecount[0][0],2));
			printf("\n %d * %d \n ------------- \n %d^2",nodecount[1][(edgepath[0][i]-1)*2],nodecount[1][(edgepath[0][j]-1)*2],nodecount[0][0]);
			printf("\n %d   %d   \n",(edgepath[0][i]-1)*2,(edgepath[0][j]-1)*2);
		}
	}




	for (i=1;i<D;i++)// for D*D states
	{	
		for (j=0;j<D;j++)
		{
			for(k=0;k<D;k++)
			{
				if (boolean[i][j]*boolean[i][k]==1) //only if both booleans are 1 
				{
					for (h=0;h<pow(D,2);h++) 
					{
						if(edgepath[j][k]==edgepath[i][h/4] && edgepath[j][k]==edgepath[i][h%4]){ 
							alpha[i][D*j+k]+=alpha[i-1][h]*nodecount[i][h%4]/nodecount[i-1][((h%4)-1)/2]*nodecount[i][h/4]/nodecount[i-1][((h/4)-1)/2];
						}
					}
				}
			
			printf("%.2f ",alpha[i][D*j+k]);
			}
			printf("\n");
		}
	}
*/


}






