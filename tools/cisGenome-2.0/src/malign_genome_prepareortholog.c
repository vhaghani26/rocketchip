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

int menu_malign_genome_prepareortholog(int argv, char **argc);

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
	menu_malign_genome_prepareortholog(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_malign_genome_prepareortholog(int argv, char **argc)
{
	/* define */
	char strInPath[LINE_LENGTH];
	char strOutPath[LINE_LENGTH];
	int nSkip = 0;
	int nRefType = 0;
	int nUp = 5000;
	int nDown = 5000;
	int nSpeciesNum = 0;

	int iOK,oOK,nOK,sfOK,rOK,upOK,downOK;
	int ni;
	int nResult;
	
	/* ------------------------------- */
	/*  malign_genome_prepareortholog  */
	/* ------------------------------- */
	if((argv < 1))
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("   malign_genome_prepareortholog   \n");
		printf("    -i input nstmap or nsomap file \n");
		printf("    -o output file                 \n");
		printf("    -n species number              \n");
		printf("    -sf skip the first species (0: No, 1: Yes) \n");
		printf("    -r reference type              \n");
		printf("       0: TSS-up, TES-down         \n");
		printf("       1: TSS-up, TSS-down         \n");
		printf("       2: TES-up, TES-down         \n");
		printf("       3: CDSS-up, CDSE-down       \n");
		printf("       4: CDSS-up, CDSS-down       \n");
		printf("       5: CDSE-up, CDSE-down       \n");
		printf("    -up up distance                \n");
		printf("    -down down distance down       \n");
		printf(" example: \n");
		printf("   malign_genome_prepareortholog -i input.txt -o output.txt -n 4 -sf 1 -r 0 -up 10000 -down 10000 \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	nOK = 0;
	sfOK = 0;
	rOK = 0;
	upOK = 0;
	downOK = 0;


	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-n") == 0)
		{
			ni++;
			nSpeciesNum = atoi(argc[ni]);
			nOK = 1;
		}
		else if(strcmp(argc[ni], "-sf") == 0)
		{
			ni++;
			nSkip = atoi(argc[ni]);
			sfOK = 1;
		}
		else if(strcmp(argc[ni], "-r") == 0)
		{
			ni++;
			nRefType = atoi(argc[ni]);
			rOK = 1;
		}
		else if(strcmp(argc[ni], "-up") == 0)
		{
			ni++;
			nUp = atoi(argc[ni]);
			if(nUp < 0)
			{
				printf("Error: -up must >=0! \n");
				exit(EXIT_FAILURE);
			}

			upOK = 1;
		}
		else if(strcmp(argc[ni], "-down") == 0)
		{
			ni++;
			nDown = atoi(argc[ni]);
			if(nDown < 0)
			{
				printf("Error: -down must >=0! \n");
				exit(EXIT_FAILURE);
			}
			downOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((iOK == 0) || (oOK == 0) || (nOK == 0) || (upOK == 0) || (downOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = MAlign_Genome_PrepareOrtholog (strInPath, strOutPath,
			nSpeciesNum, nSkip, nRefType, nUp, nDown);
	}
	
	/* return */
	return nResult;
}
