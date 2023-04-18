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

int menu_rnaseq_countgeneread(int argv, char **argc);

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
	menu_rnaseq_countgeneread(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_rnaseq_countgeneread(int argv, char **argc)
{
	char strInputPath[LINE_LENGTH];
	int nInputType = 0;
	char strOutputPath[LINE_LENGTH];
	char strDatabasePath[LINE_LENGTH];
	int nDatabaseType = 1;
	char strSpecies[LINE_LENGTH];
	int nStandardize = 1;
	int nResult;

	int iOK;
	int tOK;
	int oOK;
	int dOK;
	int dtOK;
	int sOK;
	int zOK;
	int ni;

	/* ------------------------------- */
	/*     rnaseq_countgeneread        */
	/* -i input                        */
	/* -t 0(default): input is an alignment file  */
	/*    1: input is a file that contains a list of alignment files */
	/* -o output                       */
	/* -d gene annotatation database   */
	/* -dt annotation type. 0=refGene; 1=refFlat; 2=refLocus */
	/* -s species                      */
	/* -z 1(default): standardize by total read count and length, report # reads/1M reads/1kb; 0: not standardize, report raw counts */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    rnaseq_countgeneread               \n");
		printf(" -i input \n");
		printf(" -t 0(default): input is an alignment file; \n");
		printf("    1: input is a file that contains a list of alignment files. \n");
		printf(" -o output                       \n");
		printf(" -d gene annotatation database             \n");
		printf(" -dt annotation type. 0=refGene; 1=refFlat (default); 2=refLocus \n");
		printf(" -s species \n");
		printf(" -z 1(default): standardize by total read count and length, report # reads/1M reads/1kb; 0: not standardize, report raw counts. \n");
		printf(" example: \n");
		printf("    rnaseq_countgeneread -i input.aln -o output.txt -d /data/hg18/annotation/refFlat_sorted.txt -dt 1 -s human\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	tOK = 0;
	oOK = 0;
	dOK = 0;
	dtOK = 0;
	sOK = 0;
	zOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-t") == 0)
		{
			ni++;
			nInputType = atoi(argc[ni]);
			tOK = 1;
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
			strcpy(strDatabasePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-dt") == 0)
		{
			ni++;
			nDatabaseType = atoi(argc[ni]);
			dtOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-z") == 0)
		{
			ni++;
			nStandardize = atoi(argc[ni]);
			zOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (oOK == 0) || (dOK == 0) || (sOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = RNASEQ_CountReadPerTranscript_Main(strInputPath, nInputType,
						strOutputPath, strDatabasePath, 
						nDatabaseType, strSpecies, nStandardize);
	}

	return nResult;
}
