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

int menu_fasta_soft2hardmask(int argv, char **argc);

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
	menu_fasta_soft2hardmask(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_fasta_soft2hardmask(int argv, char **argc)
{
	/* define */
	char strSeqFile[LINE_LENGTH];
	char strOutFile[LINE_LENGTH];
	int ni;
	int iOK,oOK;
	int nResult;

	/* ------------------------------- */
	/* fasta_soft2hardmask             */
	/* -i sequence                     */
	/* -o output path                  */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    fasta_soft2hardmask          \n");
		printf(" -i FASTA sequences                \n");
		printf(" -o output path \n");
		printf(" example: \n");
		printf("    fasta_soft2hardmask -i test.fa -o testm.fa\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}


	iOK = 0;
	oOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strSeqFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutFile, argc[ni]);
			oOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if( (iOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = FastaSequenceSoft2HardMask_Main(strSeqFile, strOutFile); 
	}

	return nResult;
}
