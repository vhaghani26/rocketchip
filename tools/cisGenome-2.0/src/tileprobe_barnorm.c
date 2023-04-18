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

int menu_tileprobe_barnorm(int argv, char **argc);

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
	menu_tileprobe_barnorm(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_tileprobe_barnorm(int argv, char **argc)
{
	/* ------------------------------- */
	/*   tileprobe_barnorm             */
	/* normalization based on probe effect model        */
	/* ------------------------------- */
	char strInFile[MED_LINE_LENGTH];
	char strOutFile[MED_LINE_LENGTH];
	char strModelFile[MED_LINE_LENGTH];
	int nInputType = 0;
	double dL = -100000000.0;
	int nTransform = 0;
	double dB = 0.0;

	int ni;
	int iOK,oOK,mOK,cOK,lOK,tOK,bOK;
	int nResult;

	/* ------------------------------- */
	/*  tileprobe_norm                 */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    tileprobe_barnorm           \n");
		printf(" -i input bar file or file list \n");
		printf(" -o output folder               \n");
		printf(" -m path and title of the file that contains probe effect models \n");
		printf(" -c 0 (default): input is a single bar file; 1: input is a file that contains a list of array bar files. \n");
		printf(" -l truncation lower bound, values smaller than the lower bound will be truncated to the lower bound (default = -100000000.0) \n");
		printf(" -t transformation (0:identity (default), 1:log2). Transformation will be performed after truncation \n");
		printf(" -b shrinkage factor (if <0, no shrinking; 0~1, shrinking; >1 (for future support, automatic shrinking) \n");
		printf(" example: \n");
		printf("   tileprobe_barnorm -i input.bar -o /home/normalized/ -m MouseProm\n");
		printf("   tileprobe_barnorm -i input_list.txt -c 1 -o /home/normalized/ -m MouseProm -l -100000 -t 0  -b 0.1\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	mOK = 0;
	cOK = 0;
	lOK = 0;
	tOK = 0;
	bOK = 0;
	
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
		else if(strcmp(argc[ni], "-m") == 0)
		{
			ni++;
			strcpy(strModelFile, argc[ni]);
			mOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			nInputType = atoi(argc[ni]);
			cOK = 1;
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
		else if(strcmp(argc[ni], "-b") == 0)
		{
			ni++;
			dB = atof(argc[ni]);
			bOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	/* int nCount; */
	if( (iOK == 0) || (oOK == 0) || (mOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = TileProbe_BARNorm_Main(strInFile, nInputType, strOutFile, strModelFile, nTransform, dL, dB);
	}

	return nResult;
}
