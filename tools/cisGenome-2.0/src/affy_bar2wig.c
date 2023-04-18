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

int menu_affy_bar2wig(int argv, char **argc);

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
	menu_affy_bar2wig(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_affy_bar2wig(int argv, char **argc)
{
	/* ------------------------------- */
	/*        affy_bar2wig             */
	/* convert bar file to wig file    */
	/* ------------------------------- */
	char strInputFile[LINE_LENGTH];
	char strOutputFile[LINE_LENGTH];
	int ni;
	int iOK,oOK;
	int nResult;

	/* ------------------------------- */
	/*        affy_bar2wig             */
	/* -i input file                   */
	/* -o output file                  */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("        affy_bar2wig               \n");
		printf(" -i input file                     \n");
		printf(" -o output file                    \n");
		printf(" example: \n");
		printf("    affy_bar2wig -i input.bar -o output.wig\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;

	/* set default */
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputFile, argc[ni]);
			iOK = 1;
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

	/* int nCount; */
	if( (iOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = Affy_BAR2WIG_Main(strInputFile, strOutputFile); 
	}

	return nResult;
}