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

int menu_hts_filterrepeatreads(int argv, char **argc);

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
	menu_hts_filterrepeatreads(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_hts_filterrepeatreads(int argv, char **argc)
{
	char strInputFile[MED_LINE_LENGTH];
	char strOutputFile[MED_LINE_LENGTH];
	char strGenomeFolder[MED_LINE_LENGTH];

	int nResult;
	int iOK;
	int gdOK;
	int oOK;
	int ni;

	/* ------------------------------- */
	/*   hts_filterrepeatreads         */
	/* -i input file                   */
	/* -o output title                 */
	/* -gd genome sequence folder      */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("     hts_filterrepeatreads         \n");
		printf(" -i input file     \n");
		printf(" -o output file    \n");
		printf(" -gd genome sequence folder           \n");
		printf(" example: \n");
		printf("    hts_filterrepeatreads -i input.txt -o output.txt -gd /home/hg17\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	gdOK = 0;
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
		else if(strcmp(argc[ni], "-gd") == 0)
		{
			ni++;
			strcpy(strGenomeFolder, argc[ni]);
			gdOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputFile, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (gdOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = HTS_FilterRepeatReads(strInputFile, strOutputFile, strGenomeFolder);
	}

	return nResult;
}