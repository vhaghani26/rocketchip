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

int menu_hts_collectprofile(int argv, char **argc);

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
	menu_hts_collectprofile(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_hts_collectprofile(int argv, char **argc)
{
	char strInputFile[MED_LINE_LENGTH];
	char strOutputFile[MED_LINE_LENGTH];
	char strBARFile[MED_LINE_LENGTH];
	int nW = 200;
	int nS = 20;

	int nResult;
	int iOK;
	int dOK;
	int oOK;
	int ni;

	/* ------------------------------- */
	/*      hts_collectprofile         */
	/* -i input region file            */
	/* -d input read bar file.         */
	/* -o output title                 */
	/* -w half window size             */
	/* -s step size                    */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("     hts_collectprofile           \n");
		printf(" -i input region file     \n");
		printf(" -d input alignment bar file        \n");
		printf(" -o output file           \n");
		printf(" -w half window size (default = 200 bp)     \n");
		printf(" -s step size (default = 20 bp)            \n");
		printf(" example: \n");
		printf("    hts_collectprofile -i input.txt -o output.txt -d chip.aln.bar -w 200 -s 50\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	dOK = 0;
	oOK = 0;
	
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
			strcpy(strBARFile, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputFile, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-w") == 0)
		{
			ni++;
			nW = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			nS = atoi(argc[ni]);
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
		if(nW < nS)
		{
			printf("Warning: nW cannot be smaller than nS. nW is set to nS!\n");
			nW = nS;
		}
		nResult = HTS_CollectProfile_Main(strInputFile, strOutputFile, strBARFile, nW, nS);
	}

	return nResult;
}