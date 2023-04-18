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

int menu_codinggenome(int argv, char **argc);

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
	menu_codinggenome(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_codinggenome(int argv, char **argc)
{
	/* define */
	char strGenomePath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	int nResult;
	int dOK;
	int oOK;
	int ni;

	/* ------------------------------- */
	/*        codinggenome             */
	/* -d database                     */
	/* -o output                       */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("          genome_encode         \n");
		printf(" -d path of genome database      \n");
		printf(" -o path for saving coded genome (*.sq files) \n");
		printf(" example: \n");
		printf("    genome_encode -d /data/genomes/human/b35_hg17/ -o /data/genomes/human/b35_hg17/\n\n");
		printf(" [note]                           \n");
		printf("    Before coding, there must be a chrlist.txt file to list all chromosomes that need to be coded. \n");
		printf("    After coding, a new file chrlen.txt will be created to record the length for each chromosome. \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	dOK = 0;
	oOK = 0;

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
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((dOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = Genome_Fasta_To_Code_4bit_Main(strGenomePath, strOutputPath);
	}

	return nResult;
}