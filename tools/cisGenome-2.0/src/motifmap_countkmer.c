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

int menu_motifmap_countkmer(int argv, char **argc);

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
	menu_motifmap_countkmer(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_motifmap_countkmer(int argv, char **argc)
{
	/* define */
	int dOK,iOK,oOK,kOK,cOK,chOK;
	char strInPath[MED_LINE_LENGTH];
	char strSeqFile[MED_LINE_LENGTH];
	char strOutPath[MED_LINE_LENGTH];
	char strCSAlias[MED_LINE_LENGTH];
	double dC;
	int nK = 4;
	int ni;
	int nResult;

	/* ------------------------------- */
	/*    menu_motifmap_countkmer      */
	/* -d sequence & conservation score data path */
	/* -i sequence file                */
	/* -o output file                  */
	/* -k kmer length                  */
	/* -c conservation cutoff\n")      */
	/* -ch head tag for conservation score *.cs file */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    menu_motifmap_countkmer        \n");
		printf(" -d sequence & conservation score data path\n");
		printf(" -i sequence file                  \n");
		printf(" -o output file (full path)        \n");
		printf(" -k kmer length                    \n");
		printf(" -c conservation cutoff\n");
		printf(" -ch head tag for conservation score *.cs file \n");
		printf(" example: \n");
		printf("    motifmap_countkmer -d . -i inputseq.fa -o kmer.txt -k 5 -c 200 -ch es_cluster5\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	
	dOK = 0;
	iOK = 0;
	oOK = 0;
	kOK = 0;
	cOK = 0;
	chOK = 0;
	nK = 4;
	dC = 0.0;
	strcpy(strCSAlias, "");
	

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strInPath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strSeqFile, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-k") == 0)
		{
			ni++;
			nK = atoi(argc[ni]);
			kOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			dC = atof(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-ch") == 0)
		{
			ni++;
			strcpy(strCSAlias, argc[ni]);
			chOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((dOK == 0) || (iOK == 0) || (oOK == 0) || (kOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		if( cOK == 1 )
		{
			if(chOK == 1)
			{
				nResult = MotifMap_CountKmer_Sequential_Main(strInPath, strSeqFile,  
						strOutPath, nK, cOK, dC, strCSAlias);
			}
			else
			{
				nResult = MotifMap_CountKmer_Sequential_Main(strInPath, strSeqFile,  
						strOutPath, nK, cOK, dC, "");
			}
		}
		else
		{
			nResult = MotifMap_CountKmer_Sequential_Main(strInPath, strSeqFile,
					strOutPath, nK, cOK, 0.0, "");
		}
	}

	/* return */
	return nResult;
}