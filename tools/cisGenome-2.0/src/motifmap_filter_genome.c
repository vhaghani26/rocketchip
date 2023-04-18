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

int menu_motifmap_filter_genome(int argv, char **argc);

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
	menu_motifmap_filter_genome(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_motifmap_filter_genome(int argv, char **argc)
{
	/* define */
	int iOK,oOK,cOK,cmOK,cdOK,cdsOK,cdsdOK;
	char strInputPath[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	char strMaskPath[MED_LINE_LENGTH];
	double dC;
	char strCSPath[MED_LINE_LENGTH];
	double dCds;
	char strCdsPath[MED_LINE_LENGTH];
	
	int ni;

	int nResult;
	
	/* ------------------------------- */
	/*   motifmap_filter_genome        */
	/* -i input site file              */
	/* -o output file                  */
	/* -c conservation cutoff          */
	/* -cm conservation mask           */
	/* -cd conservation score data path*/
	/* -cds coding sequence cutoff     */
	/* -cdsd coding sequence data path */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("  motifmap_filter_genome \n");
		printf(" -i input site file (full path)   \n");
		printf(" -o output file (full path) \n");
		printf(" -c conservation cutoff\n");
		printf(" -cm conservation mask\n");
		printf(" -cd conservation score data path\n");
		printf(" -cds coding sequence cutoff     \n");
		printf(" -cdsd coding sequence data path \n");
		printf(" example: \n");
		printf("    motifmap_filter_genome -i Gli.map -o Gli_filtered.map -c 1 -cm Gli_mask.txt -cd mm6/conservation -cds 0.9 -cdsd mm6/cds\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	
	iOK = 0;
	oOK = 0;
	cOK = 0;
	cmOK = 0;
	cdOK = 0;
	cdsOK = 0;
	cdsdOK = 0;
	dC = 0.0;
	dCds = 0.0;

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
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			dC = atof(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-cm") == 0)
		{
			ni++;
			strcpy(strMaskPath, argc[ni]);
			cmOK = 1;
		}
		else if(strcmp(argc[ni], "-cd") == 0)
		{
			ni++;
			strcpy(strCSPath, argc[ni]);
			cdOK = 1;
		}
		else if(strcmp(argc[ni], "-cds") == 0)
		{
			ni++;
			dCds = atof(argc[ni]);
			cdsOK = 1;
		}
		else if(strcmp(argc[ni], "-cdsd") == 0)
		{
			ni++;
			strcpy(strCdsPath, argc[ni]);
			cdsdOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((iOK == 0) || (oOK == 0) || ((cOK == 0) && (cdsOK == 0)) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	

	if(cOK == 1)
	{
		if( (cdOK == 0) || (cmOK == 0) )
		{
			printf("Error: Input Parameter not correct!\n");
			exit(EXIT_FAILURE);
		}
	}

	if(cdsOK == 1)
	{
		if(cdsdOK == 0)
		{
			printf("Error: Input Parameter not correct!\n");
			exit(EXIT_FAILURE);
		}
	}
		
	nResult = MotifMap_Filter_Genome_Main(strInputPath, strOutputPath, 
		cOK, dC, strMaskPath, strCSPath, cdsOK, dCds, strCdsPath);
	
	/* return */
	return nResult;
}
