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

int menu_motifmap_multiplematrixscan_genome(int argv, char **argc);

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
	menu_motifmap_multiplematrixscan_genome(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_motifmap_multiplematrixscan_genome(int argv, char **argc)
{
	/* define */
	int sOK,gdOK,dOK,iOK,oOK,bOK,btOK,bdOK,bsOK,cOK,cdOK,wOK;
	char strSpecies[MED_LINE_LENGTH];
	char strGenomePath[MED_LINE_LENGTH];
	char strWorkPath[MED_LINE_LENGTH];
	char strInputPath[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	int nBGOrder;
	char strBGType[MED_LINE_LENGTH];
	char strBGPath[MED_LINE_LENGTH];
	int nBGStepSize;
	int nUseCS;
	char strCSPath[MED_LINE_LENGTH];
	int nW;
	int nIncludeRepeat = 0;
	
	int ni;
	int nResult;
	
	/* --------------------------------------- */
	/*   motifmap_multiplematrixscan_genome    */
	/* -s species                              */
	/* -gd genome sequence path                */
	/* -d working directory                    */
	/* -i data files                           */
	/* -o output file                          */
	/* -b background markov order              */
	/* -bt type of background                  */
	/* -bd background mc path                  */
	/* -bs the stepsize used for computing MC  */
	/* -c use conservation score?              */
	/* -cd conservation score data path        */
	/* -w window half width for testing        */
	/*    clustering                           */
	/* -u include repeat [0: no (default); 1: using unmasked genome] */
	/* --------------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ---------------------------------- */\n");
		printf("   motifmap_multiplematrixscan_genome    \n");
		printf(" -s species\n");
		printf(" -gd genome sequence path\n");
		printf(" -d working directory \n");
		printf(" -i data (coordinates&motif) file    \n");
		printf(" -o output file (full path) \n");
		printf(" -b background markov order  \n");
		printf(" -bt type of background, [genome or region] \n");
		printf(" -bd if use genome-wide background, the path for precomputed MC matrices \n");
		printf(" -bs if use genome-wide background, the stepsize used for computing MC matrices \n");
		printf(" -c use conservation? (0: no; 1: yes) \n");
		printf(" -cd conservation score path \n");
		printf(" -w window half width for testing clustering \n");
		printf(" -u include repeat [0: no (default); 1: using unmasked genome]\n");
		printf(" example: \n");
		printf("    motifmap_multiplematrixscan_genome -s human -gd /data/genomes/human/b35_hg17/ -d . -i ChIPmotif.txt -o ChIPmap -b 3 -bt genome -bd /data/mm6/markovbg/S100000_W1000000 -bs 100000 -c 1 -cd /data/mm6/conservation/phastcons -w 100\n");
		printf("/* ---------------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	sOK = 0;	
	gdOK = 0;
	dOK = 0;
	iOK = 0;
	oOK = 0;
	bOK = 0;
	btOK = 0;
	bdOK = 0;
	bsOK = 0;
	cOK = 0;
	cdOK = 0;
	wOK = 0;

	nBGOrder = 0;
	nBGStepSize = 100000;
	nUseCS = 0;
	strcpy(strBGType, "REGION");
	strcpy(strBGPath, "");
	strcpy(strCSPath, "");
	nW = 100;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			strcpy(strSpecies, argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-gd") == 0)
		{
			ni++;
			strcpy(strGenomePath, argc[ni]);
			gdOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strWorkPath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
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
			nUseCS = atoi(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-cd") == 0)
		{
			ni++;
			strcpy(strCSPath, argc[ni]);
			cdOK = 1;
		}
		else if(strcmp(argc[ni], "-u") == 0)
		{
			ni++;
			nIncludeRepeat = atoi(argc[ni]);
			if(nIncludeRepeat != 1)
				nIncludeRepeat = 0;
		}
		else if(strcmp(argc[ni], "-w") == 0)
		{
			ni++;
			nW = atoi(argc[ni]);
			wOK = 1;
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

	if((sOK == 0) || (gdOK == 0) || (dOK == 0) || (iOK == 0) || (oOK == 0) || (bOK == 0) || (cOK == 0) || (wOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		if(nUseCS == 1)
		{
			if(cdOK == 0)
			{
				printf("Error: the conservation score path is not specified!\n");
				exit(EXIT_FAILURE);
			}
		}

		nResult = MotifMap_ScanMultipleMatrix_Genome_Main(strSpecies, strGenomePath,  
				strWorkPath, strInputPath, strOutputPath,  
				nBGOrder, strBGType, strBGPath, nBGStepSize,
				nUseCS, strCSPath, nIncludeRepeat, nW);
	}

	/* return */
	return nResult;
}