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
#include "TilingArrayLib.h"
#include "HTSequencingLib.h"

int menu_seqclust_count(int argv, char **argc);

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
	menu_seqclust_count(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_seqclust_count(int argv, char **argc)
{
	char strInputPath[LINE_LENGTH];
	char strDataPath[LINE_LENGTH];
	char strOutputFile[LINE_LENGTH];
	int nExtendLen = 0;
	
	int iOK;
	int dOK;
	int oOK;
	int ni;
	int nResult;

	/* ------------------------------- */
	/*     menu_seqclust_count         */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    seqclust_count                \n");
		printf(" -i input coordinate file \n");
		printf(" -d sample file list  \n");
		printf(" -o output file  \n");
		printf(" -e extension length (default = 0 bp), typically equal to DNA fragment length.\n");
		printf(" Example: \n");
		printf("    seqclust_count -i input_file.cod -d sample_filelist.txt -o mycountdata.txt -e 100\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	dOK = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strDataPath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputFile, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-e") == 0)
		{
			ni++;
			nExtendLen = atoi(argc[ni]);
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (oOK == 0) || (dOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = SeqClust_Count_Main(strInputPath, strDataPath, strOutputFile, nExtendLen);
	}

	return nResult;
}