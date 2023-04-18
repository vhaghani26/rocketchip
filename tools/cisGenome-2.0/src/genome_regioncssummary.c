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

int menu_regioncs_summary(int argv, char **argc);

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
	menu_regioncs_summary(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_regioncs_summary(int argv, char **argc)
{
	/* define */
	int gdOK,iOK,oOK,cOK,cdOK;
	char strGenomePath[MED_LINE_LENGTH];
	char strCodPath[MED_LINE_LENGTH];
	char strOutputPath[MED_LINE_LENGTH];
	int nUseCS;
	double dC;
	char strCSPath[MED_LINE_LENGTH];
	
	int ni;
	int nResult;
	
	/* ------------------------------- */
	/*   genome_regioncssummary        */
	/* -gd genome sequence path        */
	/* -i coordinates file             */
	/* -o output file                  */
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
		printf("   genome_regioncssummary   \n");
		printf(" -gd genome sequence path \n");
		printf(" -i coordinates file    \n");
		printf(" -o output file (full path) \n");
		printf(" -c conservation cutoff \n");
		printf(" -cd conservation score path\n");
		printf(" example: \n");
		printf("    genome_regioncssummary -gd /data/mm6 -i inputseq.cod -o inputcs.txt -c 100 -cd /data/mm6/conservation/\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

		
	gdOK = 0;
	iOK = 0;
	oOK = 0;
	cOK = 0;
	cdOK = 0;
	
	dC = -0.01;
	nUseCS = 1;
	strcpy(strCSPath, "");
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-gd") == 0)
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

	if((gdOK == 0) || (iOK == 0) || (oOK == 0) || (cdOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		cOK = 1;
		nResult = Genome_GetRegionCS_Summary_Main(strGenomePath, 
				strCodPath,	strOutputPath, 
				cOK, dC, strCSPath);
	}

	/* return */
	return nResult;
}