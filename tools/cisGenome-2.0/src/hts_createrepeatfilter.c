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

int menu_hts_createrepeatfilter(int argv, char **argc);

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
	menu_hts_createrepeatfilter(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_hts_createrepeatfilter(int argv, char **argc)
{
	char strChrListFile[MED_LINE_LENGTH];
	char strChrLenFile[MED_LINE_LENGTH];
	char strGenomeFolder[MED_LINE_LENGTH];
	char strOutputFolder[MED_LINE_LENGTH];
	char strOutputTitle[MED_LINE_LENGTH];
	int nW = 100;
	double dR = 0.5;

	int nResult;
	int iOK;
	int lOK;
	int gdOK;
	int dOK;
	int oOK;
	int wOK;
	int rOK;
	int ni;

	/* ------------------------------- */
	/*   hts_createrepeatfilter        */
	/* -i input chrlist file           */
	/* -l input chrlen file			   */
	/* -gd genome sequence folder      */
	/* -d output folder                */
	/* -o output title                 */
	/* -w window size (default = 100)  */
	/* -r repeat cutoff, if repeat percentage>=r, then discard. */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("     hts_createrepeatfilter        \n");
		printf(" -i input chrlist file     \n");
		printf(" -l input chrlen file           \n");
		printf(" -gd genome sequence folder           \n");
		printf(" -d output folder  \n");
		printf(" -o output title    \n");
		printf(" -w window size (default = 100)        \n");
		printf(" -r repeat cutoff, if repeat percentage>=r, then discard. \n");
		printf(" example: \n");
		printf("    hts_createrepeatfilter -i chrlist.txt -l chrlen.txt -gd /home/hg17 -d . -o output -w 100 -r 0.5\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	lOK = 0;
	gdOK = 0;
	dOK = 0;
	oOK = 0;
	wOK = 0;
	rOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strChrListFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			strcpy(strChrLenFile, argc[ni]);
			lOK = 1;
		}
		else if(strcmp(argc[ni], "-gd") == 0)
		{
			ni++;
			strcpy(strGenomeFolder, argc[ni]);
			gdOK = 1;
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
		else if(strcmp(argc[ni], "-r") == 0)
		{
			ni++;
			dR = atof(argc[ni]);
			rOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (lOK == 0) || (gdOK == 0) || (dOK == 0) || (oOK == 0) )
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
		
		nResult = HTS_CreateRepeatFilter(strGenomeFolder, strChrListFile, strChrLenFile, 
					 nW, dR, strOutputFolder, strOutputTitle);
	}

	return nResult;
}