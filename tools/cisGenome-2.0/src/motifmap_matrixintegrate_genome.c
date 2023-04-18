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

int menu_motifmap_matrixintegrate_genome(int argv, char **argc);

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
	menu_motifmap_matrixintegrate_genome(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_motifmap_matrixintegrate_genome(int argv, char **argc)
{
	/* define */
	int mOK,gdOK,iOK,sampOK,oOK,sOK,wOK,rOK,bOK,btOK,bdOK,bsOK,cOK,cdOK,repOK;
	char strMotifPath[MED_LINE_LENGTH];
	char strGenomePath[MED_LINE_LENGTH];
	char strCodPath[MED_LINE_LENGTH];
	char strSampPath[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	char strSpecies[MED_LINE_LENGTH];
	int nW;
	double dR;
	int nBGOrder;
	char strBGType[MED_LINE_LENGTH];
	char strBGPath[MED_LINE_LENGTH];
	int nBGStepSize;
	int nUseCS;
	double dC;
	int nRepeat;
	char strCSPath[MED_LINE_LENGTH];
	
	int ni;
	int nResult = 0;
	
	/* ------------------------------- */
	/* motifmap_matrixintegrate_genome */
	/* -m motif matrix                 */
	/* -gd genome sequence path        */
	/* -i coordinates file             */
	/* -samp sampling file             */
	/* -o output file                  */
	/* -s species                      */
	/* -w window size                  */
	/* -repeat include repeat          */
	/* -r likelihood ratio             */
	/* -b background markov order      */
	/* -bt type of background          */
	/* -bd background mc path          */
	/* -c conservation cutoff          */
	/* -cd conservation score data path*/
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("   motifmap_matrixintegrate_genome   \n");
		printf(" -m motif matrix (full path)    \n");
		printf(" -gd genome sequence path\n");
		printf(" -i coordinates file    \n");
		printf(" -samp sampling file    \n");
		printf(" -o output file (full path) \n");
		printf(" -s species             \n");
		printf(" -w window size         \n");
		printf(" -repeat include repeat \n");
		printf(" -r likelihood ratio \n");
		printf(" -b background markov order  \n");
		printf(" -bt type of background, [genome or region] \n");
		printf(" -bd if use genome-wide background, the path for precomputed MC matrices \n");
		printf(" -bs if use genome-wide background, the stepsize used for computing MC matrices \n");
		printf(" -c conservation cutoff \n");
		printf(" -cd conservation score path\n");
		printf(" example: \n");
		printf("    motifmap_matrixintegrate_genome -m Gli_matrix.txt -gd /data/mm6 -i inputseq.cod -samp affytile.txt -o Gli_map.txt -s mouse -w 200 -repeat 1 -r -1.0 -b 3 -bt genome -bd /data/mm6/MC/3 -bs 100000 -c -1 -cd /data/mm6/conservation/\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

		
	mOK = 0;
	gdOK = 0;
	iOK = 0;
	oOK = 0;
	rOK = 0;
	bOK = 0;
	btOK = 0;
	bdOK = 0;
	bsOK = 0;
	cOK = 0;
	cdOK = 0;
	sampOK = 0;
	sOK = 0;
	wOK = 0;
	repOK = 0;
	
	dR = 100.0;
	dC = 0.0;
	nW = 100;
	nBGOrder = 0;
	nBGStepSize = 100000;
	nUseCS = 0;
	nRepeat = 0;
	strcpy(strBGType, "REGION");
	strcpy(strBGPath, "");
	strcpy(strSpecies, "");

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-m") == 0)
		{
			ni++;
			strcpy(strMotifPath, argc[ni]);
			mOK = 1;
		}
		else if(strcmp(argc[ni], "-gd") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			gdOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strCodPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-samp") == 0)
		{
			ni++;
			strcpy(strSampPath, argc[ni]);
			sampOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-w") == 0)
		{
			ni++;
			nW = atoi(argc[ni]);
			wOK = 1;
		}
		else if(strcmp(argc[ni], "-repeat") == 0)
		{
			ni++;
			nRepeat = atoi(argc[ni]);
			repOK = 1;
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
		else if(strcmp(argc[ni], "-bt") == 0)
		{
			ni++;
			strcpy(strBGType, argc[ni]);
			btOK = 1;
		}
		else if(strcmp(argc[ni], "-bd") == 0)
		{
			ni++;
			strcpy(strBGPath, argc[ni]);
			bdOK = 1;
		}
		else if(strcmp(argc[ni], "-bs") == 0)
		{
			ni++;
			nBGStepSize = atoi(argc[ni]);
			bsOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			dC = atof(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-cd") == 0)
		{
			ni++;
			strcpy(strCSPath, argc[ni]);
			cdOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	StrMakeUpper(strBGType);
	if(strcmp(strBGType, "GENOME") == 0)
	{
		if(bdOK == 0)
		{
			printf("Error: the background MC path is not specified!\n");
			exit(EXIT_FAILURE);
		}
		if(bsOK == 0)
		{
			printf("Error: the background MC step size is not specified!\n");
			exit(EXIT_FAILURE);
		}
	}
	else if(strcmp(strBGType, "REGION") == 0)
	{
	}
	else
	{
		printf("Error: unknown background type parameter!\n");
		exit(EXIT_FAILURE);
	}
	
	if((mOK == 0) || (gdOK == 0) || (iOK == 0) || (oOK == 0) || (rOK == 0) || (bOK == 0)
		|| (sampOK == 0) || (sOK == 0) || (wOK == 0) || (nW < 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		if( cOK == 1 )
		{
			if(cdOK == 1)
			{
				nResult = MotifMap_ScanMatrix_Integrate_Genome_Main(strMotifPath, 
					strGenomePath, strCodPath, 
					strSampPath, strOutputPath, 
					strSpecies, nW, nRepeat, dR, 
					nBGOrder, strBGType, 
					strBGPath, nBGStepSize,
					cOK, dC, strCSPath);
			}
			else
			{
				printf("Error: the conservation score path is not specified!\n");
				exit(EXIT_FAILURE);
			}
		}
		else
		{
			nResult = MotifMap_ScanMatrix_Integrate_Genome_Main(strMotifPath, 
					strGenomePath, strCodPath, 
					strSampPath, strOutputPath, 
					strSpecies, nW, nRepeat, dR, 
					nBGOrder, strBGType, 
					strBGPath, nBGStepSize,
					cOK, dC, strCSPath);
		}
	}

	/* return */
	return nResult;
}