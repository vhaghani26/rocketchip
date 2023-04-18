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

int menu_hts_aln2winbar(int argv, char **argc);

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
	menu_hts_aln2winbar(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_hts_aln2winbar(int argv, char **argc)
{
	char strInputFile[MED_LINE_LENGTH];
	char strOutputFolder[MED_LINE_LENGTH];
	char strBARHeader[LINE_LENGTH];
	int nExtLen = 0;
	int nBinSize = 25;
	char strSpecies[LINE_LENGTH];
	char strChrLenFile[MED_LINE_LENGTH];
	int nStrand = 0;
	int nResult;
	
	int iOK;
	int dOK;
	int oOK;
	int sOK;
	int lOK;

	int ni;

	/* ------------------------------- */
	/*        hts_aln2winbar           */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("          hts_aln2winbar          \n");
		printf(" -i input aln file \n");
		printf(" -d output folder  \n");
		printf(" -o output bar file header \n");
		printf(" -e extension length \n");
		printf(" -b bin size \n");
		printf(" -strand 1: produce strand specific curves; 0 (default): do not produce strand profiles. \n");
		printf(" -s species \n");
		printf(" -l chromosome length list \n");
		printf(" example: \n");
		printf("    hts_aln2winbar -i input.aln -d . -o output -e 100 -b 25 -s human -l chrlen.txt\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	dOK = 0;
	oOK = 0;
	sOK = 0;
	lOK = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputFile, argc[ni]);
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
			strcpy(strBARHeader, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-e") == 0)
		{
			ni++;
			nExtLen = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-b") == 0)
		{
			ni++;
			nBinSize = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-strand") == 0)
		{
			ni++;
			nStrand = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			strcpy(strChrLenFile, argc[ni]);
			lOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (dOK == 0) || (oOK == 0) || (sOK == 0) || (lOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = HTS_Aln2WinBAR(strInputFile, strOutputFolder, strBARHeader, 
				   nExtLen, nBinSize, nStrand,
				   strSpecies, strChrLenFile);
	}

	return nResult;
}

