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

int menu_tilemapv2_regioninfo(int argv, char **argc);

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
	menu_tilemapv2_regioninfo(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_tilemapv2_regioninfo(int argv, char **argc)
{
	/* define */
	char strParamPath[MED_LINE_LENGTH];	
	char strRegionPath[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	int nResult;
	int ni;
	int iOK,oOK,dOK;

	/* ------------------------------- */
	/*    tilemapv2                    */
	/* -i input file                   */
	/* -o output file                  */
	/* -d track info file              */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    tilemapv2_regioninfo           \n");
		printf(" -i input file                     \n");
		printf(" -o output file                    \n");
		printf(" -d track info file                \n");
		printf(" example: \n");
		printf("    tilemapv2_regioninfo -i region.txt -o output.txt -d tilemap_info_arg.txt \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	dOK = 0;

	/* set default */
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strRegionPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strParamPath, argc[ni]);
			dOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	/* int nCount; */
	if( (iOK == 0) || (oOK == 0) || (dOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = TileMapv2_RegionInfo_Main(strRegionPath, strParamPath, strOutputPath);
	}
	
	/* return */
	return nResult;
}