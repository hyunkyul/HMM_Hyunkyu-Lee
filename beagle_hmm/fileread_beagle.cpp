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

int main()
{
//	FILE *fin;
//	char *f1;
//	char *f2;


	int i,j,h,k; //universal variables
	int nodecount[5][16] = {{600, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
						  {311, 289, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
						  {195, 116, 289, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
						  {237, 247, 116, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
						  {46, 191, 247, 0, 116, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}; // node matrix which has counts at the node
	int edgepath[5][16] = {{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
						  {1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
						  {1, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
						  {1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
						  {1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
//	int D;
/* 
	nodecount = (int**)malloc((D+1)*sizeof(int*));
	for (i=0;i<D+1;i++)
			nodecount[i] = (int*)malloc(pow(2, D)*sizeof(int));

	edgepath = (int**)malloc((D+1)*sizeof(int*));
	for (i=0;i<D+1;i++)
			edgepath[i] = (int*)malloc(pow(2, D)*sizeof(int));
*/
/*

// file read

	fin = fopen(f1,"r");
	for (i=0;i<D+1;i++)
	{
		for (j=0;j<16; j++)
		{
		fscanf(fin, "%d", &nodecount[i][j]); //ibe at a tune to address at the nodecount size
		printf("%d",nodecount[i][j]);
		}
		return 0;

	}
 
	fin = fopen(f2,"r");
 	for (i=0;i<D+1;i++)
	{
		for (j=0;j<16; j++)
		{
		fscanf(fin, "%d", &edgepath[i][j]); //ibe at a tune to address at the nodecount size
		}
	}

//file read end  
return 0;
*/


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







