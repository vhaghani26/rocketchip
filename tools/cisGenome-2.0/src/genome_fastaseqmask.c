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

int menu_seqmask(int argv, char **argc);

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
	menu_seqmask(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_seqmask(int argv, char **argc)
{
	/* define */
	char strSeqFile[LINE_LENGTH];
	char strMaskFile[LINE_LENGTH];
	char strOutFile[LINE_LENGTH];
	int ni;
	int nMaskType = 0;
	int iOK,mOK,mtOK,oOK;
	int nResult;

	/* ------------------------------- */
	/* menu_fastaseqmask               */
	/* -i sequence                     */
	/* -m masks                        */
	/* -mt mask type:                  */
	/*     0-soft mask, 1-hard mask    */
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
		printf("     genome_fastaseqmask           \n");
		printf(" -i FASTA sequences                \n");
		printf(" -m masks                          \n");
		printf(" -mt mask type (default 0)         \n");
		printf("     0: soft mask, masking with small letters a, c, g, t.\n");
		printf("     1: hard mask, masking with N \n");
		printf(" -o output path \n");
		printf(" example: \n");
		printf("    genome_fastaseqmask -i test.fa -m test.mask -mt 1 -o testm.fa\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}


	iOK = 0;
	mOK = 0;
	mtOK = 0;
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
		else if(strcmp(argc[ni], "-m") == 0)
		{
			ni++;
			strcpy(strMaskFile, argc[ni]);
			mOK = 1;
		}
		else if(strcmp(argc[ni], "-mt") == 0)
		{
			ni++;
			nMaskType = atoi(argc[ni]);
			mtOK = 1;
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

	if( (iOK == 0) ||  (mOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = FastaSequenceMask_Main(strSeqFile, strMaskFile,
			nMaskType, strOutFile); 
	}

	return nResult;
}
