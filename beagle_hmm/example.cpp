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
	char *f1;
	char *f2;

	int i,j,h,k; //universal variables
	int **nodecount; //node count
	int **edgepath; //edge path to next edge
	int D; //number of levels

	i=1;
	f1= argv[i++];
	f2= argv[i++];
	D=4;

	nodecount = (int**)malloc((D+1)*sizeof(int*));
	for (i=0;i<D+1;i++)
			nodecount[i] = (int*)malloc(pow(2, D)*sizeof(int));

	edgepath = (int**)malloc((D+1)*sizeof(int*));
	for (i=0;i<D+1;i++)
			edgepath[i] = (int*)malloc(pow(2, D)*sizeof(int));

// file read

	fin = fopen(f1,"r");
	for (i=0;i<D+1;i++)
	{
		for (j=0;j<16; j++)
		{
		fscanf(fin, "%d", &nodecount[i][j]); //ibe at a tune to address at the nodecount size
		printf("%d ",nodecount[i][j]);
		}
		printf("\n");

   }
 
	fin = fopen(f2,"r");
 	for (i=0;i<D+1;i++)
	{
		for (j=0;j<16; j++)
		{
		fscanf(fin, "%d", &edgepath[i][j]); //ibe at a tune to address at the nodecount size
		printf("%d ",edgepath[i][j]);
		}
	printf("\n");
	}

//file read end  


//file input

//print out file 
	for (i=0;i<5;i++)
	{
		for (j=0;j<16; j++)
		{
		printf("%d ",nodecount[i][j]);
		}
	printf("\n");
	}

}







