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

int menu_motifmap_consensusscan(int argv, char **argc);

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
	menu_motifmap_consensusscan(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_motifmap_consensusscan(int argv, char **argc)
{
	/* define */
	int dOK,mOK,iOK,oOK,mcOK,mdOK,cOK,chOK;
	char strMotifPath[MED_LINE_LENGTH];
	char strInputPath[MED_LINE_LENGTH];
	char strSeqFile[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	char strCSAlias[MED_LINE_LENGTH];
	int nMC,nMD;
	double dC;
	int ni;

	int nResult;
	
	/* ------------------------------- */
	/*    motifmap_consensusscan       */
	/* -m motif consensus              */
	/* -d sequence & conservation score data path */
	/* -i sequence file                */
	/* -o output file                  */
	/* -mc number of stringent         */
	/*    consensus mismatches allowed */
	/* -md number of relaxed mismatches*/
	/*    allowed                      */
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
		printf("    motifmap_consensusscan         \n");
		printf(" -m motif consensus (full path)    \n");
		printf(" -d sequence & conservation score data path\n");
		printf(" -i sequence file    \n");
		printf(" -o output file (full path) \n");
		printf(" -mc number of stringent consensus mismatches allowed\n");
		printf(" -md number of relaxed consensus mismatches allowed\n");
		printf(" -c conservation cutoff\n");
		printf(" -ch head tag for conservation score *.cs file \n");
		printf(" example: \n");
		printf("    motifmap_consensusscan -m E2F.txt -d . -i inputseq.fa -o outputseq.txt -mc 2 -md 0 -c 200 -ch es_cluster5\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	
	mOK = 0;
	dOK = 0;
	iOK = 0;
	oOK = 0;
	mcOK = 0;
	mdOK = 0;
	cOK = 0;
	chOK = 0;
	nMC = 0;
	nMD = 0;
	dC = 0.0;
	strcpy(strCSAlias, "");
	

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-m") == 0)
		{
			ni++;
			strcpy(strMotifPath, argc[ni]);
			mOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
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
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-mc") == 0)
		{
			ni++;
			nMC = atoi(argc[ni]);
			mcOK = 1;
		}
		else if(strcmp(argc[ni], "-md") == 0)
		{
			ni++;
			nMD = atoi(argc[ni]);
			mdOK = 1;
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

	

	if((mOK == 0) || (dOK == 0) || (iOK == 0) || (oOK == 0) || (mcOK == 0) || (mdOK == 0))
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
				nResult = MotifMap_ScanConsensus_Sequential_Main(strMotifPath, strInputPath, strSeqFile,  
						strOutputPath, nMC, nMD, cOK, dC, strCSAlias);
			}
			else
			{
				nResult = MotifMap_ScanConsensus_Sequential_Main(strMotifPath, strInputPath, strSeqFile,  
						strOutputPath, nMC, nMD, cOK, dC, "");
			}
		}
		else
		{
			nResult = MotifMap_ScanConsensus_Sequential_Main(strMotifPath, strInputPath, strSeqFile,
					strOutputPath, nMC, nMD, cOK, 0.0, "");
		}
	}

	/* return */
	return nResult;
}