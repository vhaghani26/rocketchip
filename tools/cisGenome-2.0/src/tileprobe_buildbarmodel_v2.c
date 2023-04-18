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

int menu_tileprobe_buildbarmodel_v2(int argv, char **argc);

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
	menu_tileprobe_buildbarmodel_v2(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_tileprobe_buildbarmodel_v2(int argv, char **argc)
{
	/* ------------------------------- */
	/*   tileprobe_buildbarmodel          */
	/* build probe effect model        */
	/* ------------------------------- */
	char strInFile[MED_LINE_LENGTH];
	char strOutFile[MED_LINE_LENGTH];
	double dL = -100000000.0;
	int nTransform = 0;
	int nShrink = 1;
	int nTest = 0;

	int ni;
	int iOK,oOK,lOK,tOK,sOK,eOK;
	int nResult;

	/* ------------------------------- */
	/*  tileprobe_buildbarmodel           */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    tileprobe_buildbarmodel_v2     \n");
		printf(" -i input file that contains the array list for building the background probe model. \n");
		printf(" -o output file title               \n");
		printf(" -l truncation lower bound, values smaller than the lower bound will be truncated to the lower bound (default = -100000000.0) \n");
		printf(" -t transformation (0=default:identity, 1:log2). Transformation will be performed after truncation \n");
		printf(" -s 1: shrink variance (default); 0: not shrink \n");
		printf(" -e 1: used by developers only, for test purpose; 0: default \n");
		printf(" example: \n");
		printf("   tileprobe_buildbarmodel_v2 -i input.txt -o MouseProm1.0R\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	lOK = 0;
	tOK = 0;
	sOK = 0;
	eOK = 0;
	
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
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			dL = atof(argc[ni]);
			lOK = 1;
		}
		else if(strcmp(argc[ni], "-t") == 0)
		{
			ni++;
			nTransform = atoi(argc[ni]);
			tOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			nShrink = atoi(argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-e") == 0)
		{
			ni++;
			nTest = atoi(argc[ni]);
			eOK = 1;
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
		nResult = TileProbe_BARBuildv2_Main(strInFile, strOutFile, nTransform, dL, nShrink, nTest);
	}

	return nResult;
}