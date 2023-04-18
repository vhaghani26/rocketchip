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

int menu_getseqfromgenome(int argv, char **argc);

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
	menu_getseqfromgenome(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_getseqfromgenome(int argv, char **argc)
{
	/* ------------------------------- */
	/*        getseqfromgenome         */
	/* get sequence from genome module */
	/* ------------------------------- */
	char strGenomePath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	char strTargetFile[LINE_LENGTH];
	char strOutputFile[LINE_LENGTH];
	char strStrandType[LINE_LENGTH];
	int nStrandType;
	int ni;
	int dOK,sOK,iOK,oOK,rOK;
	int nResult;

	/* ------------------------------- */
	/*        getseqfromgenome         */
	/* -d database                     */
	/* -s species                      */
	/* -i target region                */
	/* -o output                       */
	/* -r strand                       */
	/* get sequence from genome module */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("        genome_getseq            \n");
		printf(" -d path of genome database      \n");
		printf(" -s species name                 \n");
		printf("		currently supporting: human, mouse, dog, D_melanogaster\n");
		printf(" -i tab-delimited file for specifying target region \n");
		printf("		col1: seqid; \n");
		printf("		col2: chromosome in numbers, e.g. use 23 for human chrX, 24 for human chrY \n");
		printf("		col3: zero-based start position in the assembly \n");
		printf("		col4: zero-based end position in the assembly \n");
		printf("		col5: strand in the assembly, '+' or '-' \n");
		printf(" -o file for saving the sequences\n");
		printf(" -r method to deal with +/- strand        \n");
		printf("		genebase: always save DNA from 5' to 3' in relative to a gene; \n");
		printf("		assemblybase [default]: always save the + strand from genome assembly. \n");
		printf(" example: \n");
		printf("    genome_getseq -d /data/genomes/human/b35_hg17/ -s human -i /data/genomes/human/b35_hg17/testid.txt -o testseq.fa -r genebase \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}


	dOK = 0;
	sOK = 0;
	iOK = 0;
	oOK = 0;
	rOK = 0;

	/* set default */
	nStrandType = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strTargetFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputFile, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-r") == 0)
		{
			ni++;
			strcpy(strStrandType, argc[ni]);
			if(strcmp(strStrandType, "genebase") == 0)
			{
				nStrandType = 1;
			}
			else if(strcmp(strStrandType, "assemblybase") == 0)
			{
				nStrandType = 0;
			}
			rOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if( (dOK == 0) || (sOK == 0) || (iOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = Genome_Code_4bit_GetSeq_Main(strGenomePath, strSpecies, strTargetFile, strOutputFile, nStrandType); 
	}

	return nResult;
}