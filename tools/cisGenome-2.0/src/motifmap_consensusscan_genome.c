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

int menu_motifmap_consensusscan_genome(int argv, char **argc);

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
	menu_motifmap_consensusscan_genome(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_motifmap_consensusscan_genome(int argv, char **argc)
{
	/* define */
	int gdOK,cdOK,mOK,iOK,oOK,mcOK,mdOK,cOK;
	char strMotifPath[MED_LINE_LENGTH];
	char strGenomePath[MED_LINE_LENGTH];
	char strCSPath[MED_LINE_LENGTH];
	char strCodPath[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	int nMC,nMD;
	double dC;
	int nIncludeRepeat = 0;
	int ni;

	int nResult;
	
	/* ------------------------------- */
	/*   motifmap_consensusscan_genome */
	/* -m motif consensus              */
	/* -gd genome sequence path        */
	/* -i coordinates file             */
	/* -o output file                  */
	/* -mc number of stringent         */
	/*    consensus mismatches allowed */
	/* -md number of relaxed mismatches*/
	/*    allowed                      */
	/* -c conservation cutoff\n")      */
	/* -cd conservation score data path*/
	/* -u include repeat [0: no (default); 1: using unmasked genome] */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("  motifmap_consensusscan_genome \n");
		printf(" -m motif consensus (full path)    \n");
		printf(" -gd genome sequence path\n");
		printf(" -i coordinates file (full path)   \n");
		printf(" -o output file (full path) \n");
		printf(" -mc number of stringent consensus mismatches allowed\n");
		printf(" -md number of relaxed consensus mismatches allowed\n");
		printf(" -c conservation cutoff\n");
		printf(" -cd conservation score path\n");
		printf(" -u include repeat [0: no (default); 1: using unmasked genome]\n");
		printf(" example: \n");
		printf("    motifmap_consensusscan_genome -m E2F.txt -gd mm6 -i inputseq.cod -o outputseq.txt -mc 2 -md 0 -c 200 -cd /mm6/conservation/phastcons\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	
	mOK = 0;
	gdOK = 0;
	iOK = 0;
	oOK = 0;
	mcOK = 0;
	mdOK = 0;
	cOK = 0;
	cdOK = 0;
	nMC = 0;
	nMD = 0;
	dC = 0.0;
	

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
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((mOK == 0) || (gdOK == 0) || (iOK == 0) || (oOK == 0) || (mcOK == 0) || (mdOK == 0))
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
				nResult = MotifMap_ScanConsensus_Genome_Main(strMotifPath, strGenomePath, strCodPath,  
						strOutputPath, nMC, nMD, cOK, dC, strCSPath, nIncludeRepeat);
			}
			else
			{
				printf("Error: the path where conservation score is stored need to be specified!\n");
				exit(EXIT_FAILURE);
			}
		}
		else
		{
			nResult = MotifMap_ScanConsensus_Genome_Main(strMotifPath, strGenomePath, strCodPath,
					strOutputPath, nMC, nMD, cOK, 0.0, "", nIncludeRepeat);
		}
	}

	/* return */
	return nResult;
}
