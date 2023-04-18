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

int menu_getseqcsfromgenome_c(int argv, char **argc);

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
	menu_getseqcsfromgenome_c(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_getseqcsfromgenome_c(int argv, char **argc)
{
	/* define */
	char strGenomePath[LINE_LENGTH];
	char strConservePath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	char strTargetFile[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strSeqFile[LINE_LENGTH];
	char strStrandType[LINE_LENGTH];
	char strCSFormat[LINE_LENGTH];

	int nStrandType;
	int ni;
	int dOK,cOK,sOK,iOK,oOK,aOK,rOK,fOK;
	int nResult;

	/* ------------------------------- */
	/*   genome_getseqcs_c             */
	/* -d database                     */
	/* -c conservation                 */
	/* -s species                      */
	/* -i target region                */
	/* -o output path                  */
	/* -a output sequence file name    */
	/* -r strand                       */
	/* -f conservation format          */
	/* get sequence and conservation   */
	/* from genome                     */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    genome_getseqcs_c            \n");
		printf(" -d path of genome database      \n");
		printf(" -c path of conservation database \n");
		printf(" -s species name                 \n");
		printf("		currently supporting: human, mouse, dog, cow, chicken, zebrafish, D_melanogaster, C_elegans, yeast, arabidopsis\n");
		printf(" -i tab-delimited file for specifying target region \n");
		printf("		col1: seqid; \n");
		printf("		col2: chromosome in string, e.g. use chrX, chr1, etc. \n");
		printf("		col3: zero-based start position in the assembly \n");
		printf("		col4: zero-based end position in the assembly \n");
		printf("		col5: strand in the assembly, '+' or '-' \n");
		printf(" -o output path \n");
		printf(" -a file for saving the sequences \n");
		printf(" -r method to deal with +/- strand        \n");
		printf("		genebase: always save DNA from 5' to 3' in relative to a gene; \n");
		printf("		assemblybase [default]: always save the + strand from genome assembly. \n");
		printf(" -f conservation score format \n");
		printf("        three possible formats: cs, txt, bed\n");
		printf(" example: \n");
		printf("    genome_getseqcs_c -d /data/genomes/human/b35_hg17/ -c /data/genomes/human/b35_hg17/conservation/phastcons/ -s human -i /data/genomes/human/b35_hg17/testid.txt -o ./ -a testseq -r genebase -f txt\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}


	dOK = 0;
	cOK = 0;
	sOK = 0;
	iOK = 0;
	oOK = 0;
	aOK = 0;
	rOK = 0;
	fOK = 0;

	/* set default */
	nStrandType = 0;
	strcpy(strCSFormat, "cs");

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			strcpy(strConservePath, argc[ni]);
			cOK = 1;
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
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-a") == 0)
		{
			ni++;
			strcpy(strSeqFile, argc[ni]);
			aOK = 1;
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
		else if(strcmp(argc[ni], "-f") == 0)
		{
			ni++;
			strcpy(strCSFormat, argc[ni]);
			fOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if( (dOK == 0) || (cOK == 0) ||  (sOK == 0) || (iOK == 0) || (oOK == 0) || (aOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = Genome_Code_4bit_GetSeqCS_C_Main(strGenomePath, strConservePath, strSpecies, 
			strTargetFile, strOutputPath, strSeqFile, nStrandType, strCSFormat); 
	}

	return nResult;
}
