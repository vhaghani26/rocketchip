#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
#include "limits.h"
#include "time.h"

#include "StringLib.h"
#include "MatrixLib.h"
#include "RandomLib.h"
#include "MathLib.h"
#include "MotifLib.h"
#include "SequenceLib.h"
#include "GenomeLib.h"
#include "MicroarrayLib.h"
#include "AffyLib.h"
#include "WorkLib.h"

int menu_motifmap_getcluster(int argv, char **argc);

int main(int argv, char **argc)
{
	int nLen;
	int nseed;

	/* init rand */
	srand( (unsigned)time( NULL ) );
	if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
	{
		nseed = (int)(rand()*1000/RAND_MAX);
	}
	else
	{
		nseed = rand()%1000;
	}
	rand_u_init(nseed);


	/* ---- */
	/* menu */
	/* ---- */
	menu_motifmap_getcluster(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_motifmap_getcluster(int argv, char **argc)
{
	/* define */
	int iOK,wOK,oOK,rOK;
	char strInputPath[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	int nW = 500;
	int nInputType = 0;
	int ni;

	int nResult;
	
	/* ------------------------------- */
	/*    menu_motifmap_getcluster     */
	/* -i input coordinates file       */
	/* -w window size                  */ 
	/* -o output file                  */
	/* -r input type, 0:cod, 1:bed,    */
	/*                2:codp,3:bedp    */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    menu_motifmap_getcluster  \n");
		printf(" -i input coordinates file   \n");
		printf(" -w window size [default = 500] \n");
		printf(" -o output file \n");
		printf(" -r input type (0:cod, 1:bed, 2:codp,3:bedp) [default = 0] \n");
		printf(" example: \n");
		printf("    menu_motifmap_getcluster -i Gli_map.txt -w 200 -o Gli_clusteredsites.txt -r 0\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	
	iOK = 0;
	wOK = 0;
	oOK = 0;
	rOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-w") == 0)
		{
			ni++;
			nW = atoi(argc[ni]);
			wOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-r") == 0)
		{
			ni++;
			nInputType = atoi(argc[ni]);
			rOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((iOK == 0) || (oOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = MotifMap_getcluster_Main(strInputPath, strOutputPath, nW, nInputType);
	}

	/* return */
	return nResult;
}