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

int menu_tilemapv2_txt2bar(int argv, char **argc);

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
	menu_tilemapv2_txt2bar(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_tilemapv2_txt2bar(int argv, char **argc)
{
	/* ------------------------------- */
	/*     tilemapv2_txt2bar           */
	/* convert bar file to txt file    */
	/* ------------------------------- */
	char strInputFile[LINE_LENGTH];
	char strOutputFolder[LINE_LENGTH];
	char strOutputFile[LINE_LENGTH];
	int ni;
	int iOK,oOK,dOK;
	int nResult;

	/* ------------------------------- */
	/*   tilemapv2_txt2bar             */
	/* -i input file                   */
	/* -d output folder                */
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
		printf("    tilemapv2_txt2bar               \n");
		printf(" -i input file (full path)         \n");
		printf(" -d output folder (default = current folder) \n");
		printf(" -o output file name               \n");
		printf(" example: \n");
		printf("   tilemapv2_txt2bar -i input.txt -d ./data -o output.cgw\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	dOK = 0;
	oOK = 0;
	
	strcpy(strOutputFolder, ".");

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
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strOutputFolder, argc[ni]);
			dOK = 1;
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
		nResult = TileMapv2_TXT2BAR(strInputFile, strOutputFolder, strOutputFile); 
	}

	return nResult;
}
