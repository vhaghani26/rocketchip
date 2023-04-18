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

int menu_motifmap_matrixscan(int argv, char **argc);

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
	menu_motifmap_matrixscan(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_motifmap_matrixscan(int argv, char **argc)
{
	/* define */
	int dOK,mOK,iOK,oOK,rOK,bOK,cOK,chOK;
	char strMotifPath[MED_LINE_LENGTH];
	char strInputPath[MED_LINE_LENGTH];
	char strSeqFile[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	char strCSAlias[MED_LINE_LENGTH];
	double dC;
	double dR;
	int nBGOrder;
	int ni;
	int nUseCS;

	int nResult;
	
	/* ------------------------------- */
	/*    motifmap_matrixscan          */
	/* -m motif matrix                 */
	/* -d sequence & conservation path */
	/* -i sequence file                */
	/* -o output file                  */
	/* -r likelihood ratio             */
	/* -b background markov order      */
	/* -c conservation cutoff          */
	/* -ch head tag for *.cs file      */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    motifmap_matrixscan            \n");
		printf(" -m motif matrix (full path)    \n");
		printf(" -d sequence & conservation score data path\n");
		printf(" -i sequence file    \n");
		printf(" -o output file (full path) \n");
		printf(" -r likelihood ratio \n");
		printf(" -b background markov order  \n");
		printf(" -c conservation cutoff \n");
		printf(" -ch head tag for conservation score *.cs file \n");
		printf(" example: \n");
		printf("    motifmap_matrixscan -m E2F.txt -d . -i inputseq.fa -o outputseq.txt -r 100 -b 3 -c 100 -ch es_cluster5\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	
	mOK = 0;
	dOK = 0;
	iOK = 0;
	oOK = 0;
	rOK = 0;
	bOK = 0;
	cOK = 0;
	chOK = 0;
	
	dC = 0.0;
	dR = 100.0;
	nBGOrder = 0;
	

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
		else if(strcmp(argc[ni], "-r") == 0)
		{
			ni++;
			dR = atof(argc[ni]);
			rOK = 1;
		}
		else if(strcmp(argc[ni], "-b") == 0)
		{
			ni++;
			nBGOrder = atoi(argc[ni]);
			bOK = 1;
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

	

	if((mOK == 0) || (dOK == 0) || (iOK == 0) || (oOK == 0) || (rOK == 0) || (bOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		if( cOK == 1 )
		{
			nUseCS = 1;
			if(chOK == 1)
			{
				nResult = MotifMap_ScanMatrix_Sequential_Main(strMotifPath, strInputPath, strSeqFile,  
					strOutputPath, dR, nBGOrder, nUseCS, dC, strCSAlias);
			}
			else
			{
				nResult = MotifMap_ScanMatrix_Sequential_Main(strMotifPath, strInputPath, strSeqFile,  
					strOutputPath, dR, nBGOrder, nUseCS, dC, "");
			}
		}
		else
		{
			nUseCS = 0;
			nResult = MotifMap_ScanMatrix_Sequential_Main(strMotifPath, strInputPath, strSeqFile,
				strOutputPath, dR, nBGOrder, nUseCS, 0.0, "");
		}
	}

	/* return */
	return nResult;
}