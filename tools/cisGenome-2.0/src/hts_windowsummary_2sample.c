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

int menu_hts_twosample_windowsummary(int argv, char **argc);

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
	menu_hts_twosample_windowsummary(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_hts_twosample_windowsummary(int argv, char **argc)
{
	char strInputPath[LINE_LENGTH];
	char strNegInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	int nW = 100;
	char strChrList[MED_LINE_LENGTH];
	char strChrLen[MED_LINE_LENGTH];
	int nResult;
	int iOK;
	int nOK;
	int oOK;
	int wOK;
	int gOK;
	int lOK;
	int ni;

	/* ------------------------------- */
	/*  hts_windowsummary_twosample    */
	/* -i input1 (positive sample)     */
	/* -n input2 (negative sample)     */
	/* -o output                       */
	/* -w window size (default = 100)  */
	/* -g chromosome list              */
	/* -l chromosome length            */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    hts_windowsummary_twosample    \n");
		printf(" -i input (positive)     \n");
		printf(" -n input (negative)     \n");
		printf(" -o output     \n");
		printf(" -w window size (default = 100)  \n");
		printf(" -g chromosome list              \n");
		printf(" -l chromosome length            \n");
		printf(" example: \n");
		printf("    hts_windowsummary_twosample -i posinput.bar -n neginput.bar -g chrlist.txt -l chrlen.txt -w 200 -o output.txt \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	nOK = 0;
	oOK = 0;
	wOK = 0;
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
		else if(strcmp(argc[ni], "-n") == 0)
		{
			ni++;
			strcpy(strNegInputPath, argc[ni]);
			nOK = 1;
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
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (nOK == 0) || (oOK == 0) || (gOK == 0) || (lOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = HTS_TwoSample_WindowSummary(strInputPath, strNegInputPath, strChrList, strChrLen, nW, strOutputPath);
	}

	return nResult;
}