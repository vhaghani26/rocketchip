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
#include "HTSequencingLib.h"

int menu_hts_onesample_enrich(int argv, char **argc);

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
	menu_hts_onesample_enrich(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_hts_onesample_enrich(int argv, char **argc)
{
	char strInputPath[MED_LINE_LENGTH];
	char strOutputFolder[MED_LINE_LENGTH];
	char strOutputTitle[MED_LINE_LENGTH];
	int nW = 100;
	int nS = 25;
	int nCut = 10;
	int nMinLen = 1;
	int nMaxGap = 0;

	int nResult;
	int iOK;
	int dOK;
	int oOK;
	int wOK;
	int sOK;
	int cOK;
	int gOK;
	int lOK;
	int ni;

	/* ------------------------------- */
	/*    hts_onesample_enrich         */
	/* -i input                        */
	/* -d output folder                */
	/* -o output title                 */
	/* -w window size (default = 100)  */
	/* -s step size (default = 25)     */
	/* -c cutoff (default = 10)        */
	/* -g max gap (default = 0)        */
	/* -l min region length (default = 1) */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("       hts_peakdetector            \n");
		printf(" -i input      \n");
		printf(" -d output folder             \n");
		printf(" -o output title              \n");
		printf(" -w window size (default = 100)  \n");
		printf(" -s step size (default = 25)     \n");
		printf(" -c cutoff (default = 10)        \n");
		printf(" -g max gap (default = 0)        \n");
		printf(" -l min region length (default = 1) \n");
		printf(" example: \n");
		printf("    hts_peakdetector -i input.bar -d . -o output -w 200 -s 25 -c 10 -g 50 -l 500\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	dOK = 0;
	oOK = 0;
	wOK = 0;
	sOK = 0;
	cOK = 0;
	gOK = 0;
	lOK = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strOutputFolder, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputTitle, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-w") == 0)
		{
			ni++;
			nW = atoi(argc[ni]);
			wOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			nS = atoi(argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			nCut = atoi(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-g") == 0)
		{
			ni++;
			nMaxGap = atoi(argc[ni]);
			gOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			nMinLen = atoi(argc[ni]);
			lOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (dOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		if(nW <= 0)
		{
			nW = 100;
			printf("Warning: Window size<=0! Proceed with default window size = 100.\n");
		}
		if(nS <= 0)
		{
			nS = 25;
			printf("Warning: Step size<=0! Proceed with default step size = 25.\n");
		}
		if(nCut <= 0)
		{
			nCut = 10;
			printf("Warning: Cutoff<=0! Proceed with default cutoff = 10.\n");
		}

		nResult =  HTS_Enrich_OneSample_Main(strInputPath, nW, nS, nCut,
					  nMinLen, nMaxGap, strOutputFolder, strOutputTitle);
	}

	return nResult;
}