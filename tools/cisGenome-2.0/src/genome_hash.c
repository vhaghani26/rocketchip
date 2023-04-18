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

int menu_genome_hash(int argv, char **argc);

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
	menu_genome_hash(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_genome_hash(int argv, char **argc)
{
	/* define */
	char strGenomePath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	int nK = 13;

	int ni;
	int dOK,oOK,kOK;
	int nResult;

	/* ------------------------------- */
	/* genome_hash                     */
	/* -d database                     */
	/* -o output path                  */
	/* -k key length                 */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("     genome_hash                 \n");
		printf(" -d path of genome database      \n");
		printf(" -o output path \n");
		printf(" -k key length \n");
		printf(" example: \n");
		printf("    genome_hash -d /data/genomes/human/b35_hg17/ -o ./ -k 14\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}


	dOK = 0;
	oOK = 0;
	kOK = 0;

	/* set default */
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-k") == 0)
		{
			ni++;
			nK = atoi(argc[ni]);
			kOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if( (dOK == 0) || (oOK == 0) || (kOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}

	nResult = Genome_CreateHash_Main(strGenomePath, strOutputPath, 
						   nK);

	return nResult;
}