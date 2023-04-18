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

int menu_motifmap_getsitearoundcs_genome(int argv, char **argc);

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
	menu_motifmap_getsitearoundcs_genome(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_motifmap_getsitearoundcs_genome(int argv, char **argc)
{
	/* define */
	int iOK,oOK,lOK,wOK,sOK,gdOK,cdOK;
	char strInputPath[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	char strGenomePath[MED_LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	int nW,nMotifLen;
	char strCSPath[MED_LINE_LENGTH];
	
	int ni;

	int nResult;
	
	/* ------------------------------- */
	/* motifmap_getsitearoundcs        */
	/* -i input site file              */
	/* -o output file                  */
	/* -l motif length                 */
	/* -w flanking window size         */
	/* -s species                      */
	/* -gd genome sequence data path   */
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
		printf("  motifmap_getsitearoundcs        \n");
		printf(" -i input site file (full path)   \n");
		printf(" -o output file (full path) \n");
		printf(" -l motif length \n");
		printf(" -w flanking window size\n");
		printf(" -s species (current support: human, mouse, dog, cow, chicken, zebrafish, D_melanogaster, C_elegans, yeast, arabidopsis) \n");
		printf(" -gd genome sequence data path\n");
		printf(" -cd conservation score data path\n");
		printf(" example: \n");
		printf("    motifmap_getsitearoundcs -i Gli.map -o Gli_cscurve.txt -l 12 -w 50 -s mouse -gd mm6/ -cd mm6/conservation\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	
	iOK = 0;
	oOK = 0;
	lOK = 0;
	wOK = 0;
	sOK = 0;
	gdOK = 0;
	cdOK = 0;
	nW = 0;
	nMotifLen = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
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
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			nMotifLen = atoi(argc[ni]);
			lOK = 1;
		}
		else if(strcmp(argc[ni], "-w") == 0)
		{
			ni++;
			nW = atoi(argc[ni]);
			wOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
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

	if((iOK == 0) || (oOK == 0) || (lOK == 0) || (wOK == 0) || (sOK == 0) || (gdOK == 0) || (cdOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	

	nResult = MotifMap_GetSiteAroundCS_Genome_Main(strInputPath, strOutputPath, 
		nMotifLen, nW, strSpecies, strGenomePath, strCSPath);
	
	/* return */
	return nResult;
}