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

int menu_motifmap_getsitearound(int argv, char **argc);

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
	menu_motifmap_getsitearound(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_motifmap_getsitearound(int argv, char **argc)
{
	/* define */
	int iOK,wOK,dOK,sOK,oOK,cnOK,aOK;
	char strInputPath[MED_LINE_LENGTH];
	char strGenomePath[MED_LINE_LENGTH];
	char strSpecies[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	int nW,nCN,nA;
	int ni;

	int nResult;
	
	/* ------------------------------- */
	/*     motifmap_getsitearound      */
	/* -i input coordinates file       */
	/* -w half flanking window size    */ 
	/* -d genome sequence path         */
	/* -s species                      */
	/* -cn numerical chromosome id     */
	/* -o output file                  */
	/* -a alias type                   */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    motifmap_getsitearound  \n");
		printf(" -i input coordinates file   \n");
		printf(" -w half flanking window size\n");
		printf(" -d genome sequence path     \n");
		printf(" -s species                  \n");
		printf(" -cn numerical chromosome id (0: string chr id; 1: integer chr id\n");
		printf(" -o output file \n");
		printf(" -a alias type (0: use original, 1: reorder) \n");
		printf(" example: \n");
		printf("    motifmap_getsitearound -i Gli_map.txt -w 200 -d /data/genomes/mouse/mm6/ -s mouse -cn 1 -a 1 -o Gli_sitearound.txt\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	
	iOK = 0;
	wOK = 0;
	dOK = 0;
	sOK = 0;
	cnOK = 0;
	oOK = 0;
	aOK = 0;
	
	nW = 0;
	nCN = 0;
	nA = 0;
	
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
		
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-a") == 0)
		{
			ni++;
			nA = atoi(argc[ni]);
			aOK = 1;
		}
		else if(strcmp(argc[ni], "-cn") == 0)
		{
			ni++;
			nCN = atoi(argc[ni]);
			cnOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((iOK == 0) || (wOK == 0) || (dOK == 0) || (sOK == 0) || (oOK == 0) || (aOK == 0) || (cnOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = MotifMap_GetSiteAround_Main(strInputPath, strGenomePath, 
			strOutputPath, strSpecies, nW, nCN, nA);
	}

	/* return */
	return nResult;
}