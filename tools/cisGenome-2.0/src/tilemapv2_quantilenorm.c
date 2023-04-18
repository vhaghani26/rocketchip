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
#include "TilingArrayLib.h"
#include "AffyLib.h"
#include "WorkLib.h"

int menu_tilemapv2_quantilenorm(int argv, char **argc);

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
	menu_tilemapv2_quantilenorm(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_tilemapv2_quantilenorm(int argv, char **argc)
{
	/* ------------------------------- */
	/*   tilemapv2_quantilenorm        */
	/* quantile normalizatio of a txt  */
	/* ------------------------------- */
	char strInFile[LINE_LENGTH];
	char strOutFile[LINE_LENGTH];
	int ni;
	int iOK,oOK,cOK,lOK,uOK,tOK;
	int nResult;
	int nSkipCol = 2;
	double dL = -1000000.0;
	double dU = 1000000.0;
	int nTransform = 0;

	/* ------------------------------- */
	/*  tilemapv2_quantilenorm         */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    tilemapv2_quantilenorm         \n");
		printf(" -i input file (full path)         \n");
		printf(" -o output file name               \n");
		printf(" -c number of columns to skip (default = 2) \n");
		printf(" -l truncation lower bound, values smaller than the lower bound will be truncated to the lower bound (default = -1000000) \n");
		printf(" -u truncation upper bound, values bigger than the upper bound will be truncated to the upper bound (default = 1000000) \n");
		printf(" -t transformation (0:identity(default), 1:log2). Transformation will be performed after truncation \n");
		printf(" example: \n");
		printf("   tilemapv2_quantilenorm -i input.txt -o output.txt -c 1 -l 1 -t l\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	cOK = 0;
	lOK = 0;
	uOK = 0;
	tOK = 0;
	
	/* set default */
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutFile, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			nSkipCol = atoi(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			dL = atof(argc[ni]);
			lOK = 1;
		}
		else if(strcmp(argc[ni], "-u") == 0)
		{
			ni++;
			dU = atof(argc[ni]);
			uOK = 1;
		}
		else if(strcmp(argc[ni], "-t") == 0)
		{
			ni++;
			nTransform = atoi(argc[ni]);
			tOK = 1;
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
		nResult = TileMapv2_TXT_QuantileNormalization(strInFile, strOutFile,
				nSkipCol, nTransform, dL, dU);

	}

	return nResult;
}
