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

int menu_motifmap_matrixscan_genomebackground(int argv, char **argc);

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
	menu_motifmap_matrixscan_genomebackground(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_motifmap_matrixscan_genomebackground(int argv, char **argc)
{
	/* define */
	int dOK,oOK,bOK,sOK,wOK;
	char strGenomePath[MED_LINE_LENGTH];
	char strOutPath[MED_LINE_LENGTH];
	int nBGOrder;
	int nS;
	int nW;
	
	int ni;
	int nResult;
	
	/* ------------------------------- */
	/*  motifmap_matrixscan_genomebg   */
	/* -d genome database path         */
	/* -o exporting directory          */
	/* -b background markov order      */
	/* -s step size                    */
	/* -w window size                  */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("   motifmap_matrixscan_genomebg    \n");
		printf(" -d genome database path \n");
		printf(" -o exporting directory \n");
		printf(" -b background markov order  \n");
		printf(" -s step size (should <= 1000000)\n");
		printf(" -w window size\n");
		printf(" [note: w|s must be 0, i.e. w=a*s, where a is an integer]");
		printf(" example: \n");
		printf("    motifmap_matrixscan_genomebg -d /data/mm6 -o /data/mm6/mcbg/3/ -b 3 -s 100000 -w 1000000\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

		
	dOK = 0;
	oOK = 0;
	bOK = 0;
	sOK = 0;
	wOK = 0;
	
	nBGOrder = 0;
	nS = 100000;
	nW = 1000000;
	
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-b") == 0)
		{
			ni++;
			nBGOrder = atoi(argc[ni]);
			bOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			nS = atoi(argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-w") == 0)
		{
			ni++;
			nW = atoi(argc[ni]);
			wOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((dOK == 0) || (oOK == 0) || (bOK == 0) || (sOK == 0) || (wOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = MotifMap_ScanMatrix_GenomeBackground_Main(strGenomePath, 
					strOutPath, nBGOrder, nS, nW);
	}

	/* return */
	return nResult;
}