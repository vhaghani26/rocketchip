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

int menu_refgene_getmultiorthologtssaround(int argv, char **argc);

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
	menu_refgene_getmultiorthologtssaround(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_refgene_getmultiorthologtssaround(int argv, char **argc)
{
	/* define */
	char strTargetPath[LINE_LENGTH];
	char strOutPath[LINE_LENGTH];
	char strParamPath[LINE_LENGTH];	
	char strSeqFile[LINE_LENGTH];
	int nUp;
	int nDown;

	int nResult;
	int iOK;
	int dOK;
	int oOK;
	int upOK;
	int downOK;
	int aOK;
	int ni;

	/* ------------------------------- */
	/*    refgene_getmultiorthologtss  */
	/* -i target list                  */
	/* -d database info                */
	/* -o output                       */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    refgene_getmultiorthologtssaround      \n");
		printf(" -i path of target list           \n");
		printf(" -d database infomation           \n");
		printf(" -up TSS up                       \n");
		printf(" -down TSS down                   \n");
		printf(" -o path for saving retrieved orthologs \n");
		printf(" -a file for saving the sequences \n");
		printf(" example: \n");
		printf("    refgene_getmultiorthologtssaround -i testortho.msomap -up -5000 -down 1000 -d genomesetting.txt -o ./ -a testortho\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	dOK = 0;
	oOK = 0;
	upOK = 0;
	downOK = 0;
	aOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strTargetPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-up") == 0)
		{
			ni++;
			nUp = atoi(argc[ni]);
			upOK = 1;
		}
		else if(strcmp(argc[ni], "-down") == 0)
		{
			ni++;
			nDown = atoi(argc[ni]);
			downOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strParamPath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-a") == 0)
		{
			ni++;
			strcpy(strSeqFile, argc[ni]);
			aOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	

	if((iOK == 0) || (dOK == 0) || (oOK == 0) || (upOK == 0) || (downOK == 0) || (aOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = RefGene_MultiOrthologGetTSSAround(strTargetPath, nUp, nDown, strParamPath, strOutPath, strSeqFile);
	}

	/* return */
	return nResult;
}