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

int menu_genome_regionextend(int argv, char **argc);

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
	menu_genome_regionextend(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_genome_regionextend(int argv, char **argc)
{
	/* define */
	int iOK,lOK,rOK,dOK,sOK,oOK,cnOK,aOK,tOK;
	char strInputPath[MED_LINE_LENGTH];
	char strGenomePath[MED_LINE_LENGTH];
	char strSpecies[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	int nL,nR,nCN,nA,nUseStrand;
	int ni;

	int nResult;
	
	/* ------------------------------- */
	/*    genome_regionextend          */
	/* -i input coordinates file       */
	/* -l left extension length        */
	/* -r right extension length       */
	/* -t left/right based on strand info */
	/*    0: ignore strand, assemblybase (default)*/
	/*    1: +/- base                   */
	/* -d genome sequence path         */
	/* -s species                      */
	/* -cn numerical chromosome id     */
	/* -o output file                  */
	/* -a alias type                   */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    genome_regionextend  \n");
		printf(" -i input coordinates file   \n");
		printf(" -l left extension length (default = 0) \n");
		printf(" -r right extension length (default = 0) \n");
		printf(" -t left/right based on strand info   \n");
		printf("    0: ignore strand, assemblybase (default)  \n");
		printf("    1: 5' left, 3' right; 5'/3' determined by +/- \n");
		printf(" -d genome sequence path     \n");
		printf(" -s species  (current support: human, mouse, dog, cow, chicken, zebrafish, D_melanogaster, C_elegans, yeast, arabidopsis) \n");
		printf(" -cn numerical chromosome id [0=string chr name (default); 1=integer chr id]\n");
		printf(" -o output file \n");
		printf(" -a alias type [0=use original (default); 1=reorder] \n");
		printf(" example: \n");
		printf("    genome_regionextend -i Gli_map.txt -l 100 -r 200 -t 1 -d /data/genomes/mouse/mm6/ -s mouse -cn 0 -a 1 -o Gli_sitearound.txt\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	
	iOK = 0;
	lOK = 0;
	rOK = 0;
	dOK = 0;
	sOK = 0;
	cnOK = 0;
	oOK = 0;
	aOK = 0;
	tOK = 0;
	
	nL = 0;
	nR = 0;
	nCN = 0;
	nA = 0;
	nUseStrand = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			nL = atoi(argc[ni]);
			lOK = 1;
		}
		else if(strcmp(argc[ni], "-r") == 0)
		{
			ni++;
			nR = atoi(argc[ni]);
			rOK = 1;
		}
		else if(strcmp(argc[ni], "-t") == 0)
		{
			ni++;
			nUseStrand = atoi(argc[ni]);
			tOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
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
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-a") == 0)
		{
			ni++;
			nA = atoi(argc[ni]);
			aOK = 1;
		}
		else if(strcmp(argc[ni], "-cn") == 0)
		{
			ni++;
			nCN = atoi(argc[ni]);
			cnOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((iOK == 0) || (dOK == 0) || (sOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = Genome_RegionExtend_Main(strInputPath, strGenomePath, 
			strOutputPath, strSpecies, nL, nR, nCN, nA, nUseStrand);
	}

	/* return */
	return nResult;
}
