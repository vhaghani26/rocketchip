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

int menu_hts_windowsummaryv2(int argv, char **argc);

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
	menu_hts_windowsummaryv2(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_hts_windowsummaryv2(int argv, char **argc)
{
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	int nW = 100;
	int nZ = 0;
	char strChrList[MED_LINE_LENGTH];
	char strChrLen[MED_LINE_LENGTH];
	char strRepeatFile[MED_LINE_LENGTH];
	int nResult;
	int iOK;
	int oOK;
	int wOK;
	int gOK;
	int lOK;
	int zOK;
	int mOK;
	int ni;

	/* ------------------------------- */
	/*     hts_windowsummaryv2         */
	/* -i input                        */
	/* -o output                       */
	/* -w window size (default = 100)  */
	/* -g chromosome list              */
	/* -l chromosome length            */
	/* -z use combined data after shifting */
	/* -m mask repeats                 */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("        hts_windowsummaryv2        \n");
		printf(" -i input      \n");
		printf(" -o output     \n");
		printf(" -w window size (default = 100)  \n");
		printf(" -g chromosome list              \n");
		printf(" -l chromosome length            \n");
		printf(" -z use combined data after shifting (default = 0) \n");
		printf(" -m mask windows specified in a BAR file which usually represents repeat regions (can be skipped) \n");
		printf(" example: \n");
		printf("    hts_windowsummaryv2 -i input.bar -g chrlist.txt -l chrlen.txt -w 200 -o output.txt -z 1 -m repeat.bar\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	wOK = 0;
	gOK = 0;
	lOK = 0;
	mOK = 0;
	strcpy(strRepeatFile, "");

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-g") == 0)
		{
			ni++;
			strcpy(strChrList, argc[ni]);
			gOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			strcpy(strChrLen, argc[ni]);
			lOK = 1;
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
		else if(strcmp(argc[ni], "-z") == 0)
		{
			ni++;
			nZ = atoi(argc[ni]);
			zOK = 1;
		}
		else if(strcmp(argc[ni], "-m") == 0)
		{
			ni++;
			strcpy(strRepeatFile, argc[ni]);
			mOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (oOK == 0) || (gOK == 0) || (lOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = HTS_WindowSummaryv2(strInputPath, strChrList, strChrLen, nW, strOutputPath, nZ, strRepeatFile);
	}

	return nResult;
}