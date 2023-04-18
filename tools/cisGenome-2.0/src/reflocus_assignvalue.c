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

int menu_reflocus_assignvalue(int argv, char **argc);

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
	menu_reflocus_assignvalue(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_reflocus_assignvalue(int argv, char **argc)
{
	/* define */
	char strRefLocusDatabasePath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	char strOutPath[LINE_LENGTH];
	char strScorePath[LINE_LENGTH];	
	double dTruncateLowerBound; 
	double dTruncateUpperBound;
	char strTransform[LINE_LENGTH];
	int nNormalize = 0;
	double dPostLowerBound; 
	double dPostUpperBound;
	char strPostTransform[LINE_LENGTH]; 
	int nTakeAbsoluteValue = 0;
	char strRefToProbeRowIDPath[LINE_LENGTH];
	char strRefToProbeNamePath[LINE_LENGTH];
	char strNetworkPath[LINE_LENGTH];
	char strNetworkAnnotationPath[LINE_LENGTH];
	int nNetDepth = 1;
		
	int nResult;
	int dOK = 0;
	int sOK = 0;
	int oOK = 0;
	int vOK = 0;
	int vlOK = 0;
	int vuOK = 0;
	int vtOK = 0;
	int nOK = 0;
	int plOK = 0;
	int puOK = 0;
	int ptOK = 0;
	int absOK = 0;
	int d2rOK = 0;
	int d2aOK = 0;
	int netOK = 0;
	int netAOK = 0;
	int netDOK = 0;

	int ni;

	/* ------------------------------- */
	/* menu_reflocus_assignvalue       */
	/* -d reflocus database path       */
	/* -s species                      */
	/* -o output path                  */
	/* -v score path                   */
	/* -vl score truncation lower bound*/
	/* -vu score truncation upper bound*/
	/* -vt score transformation        */
	/* -nOK normalization indicator    */
	/* -pl post normalization truncation lower bound*/
	/* -pu post normalization truncation upper bound*/
	/* -pt post normalization transformation        */
	/* -absOK take absolute value      */
	/* -d2r reflocus to score rowid map*/
	/* -d2a reflocus to probe name map */
	/* -net use physical network to enhance inference*/
	/* -netAOK network annotation      */
	/* -netDOK network search depth    */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("   reflocus_assignvalue            \n");
		printf(" -d reflocus database path       \n");
		printf(" -s species                      \n");
		printf(" -o output path                  \n");
		printf(" -v score path                   \n");
		printf(" -vl score truncation lower bound\n");
		printf(" -vu score truncation upper bound\n");
		printf(" -vt score transformation        \n");
		printf(" -nOK normalization indicator    \n");
		printf(" -pl post normalization truncation lower bound\n");
		printf(" -pu post normalization truncation upper bound\n");
		printf(" -pt post normalization transformation        \n");
		printf(" -abs take absolute value        \n");
		printf(" -d2r reflocus to score rowid map\n");
		printf(" -d2a reflocus to probe name map \n");
		printf(" -net use physical network to enhance inference\n");
		printf(" -netAOK network annotation      \n");
		printf(" -netDOK network search depth    \n");
		printf(" example: \n");
		printf("    reflocus_assignvalue -d refLocus_sorted.txt -s mouse -o refLocus_POSPOST.txt -v Positive_Post.ori -vl 0.0005 -vu 0.9995 -vt logit -n 0 -pl -1e20 -pu 1e20 -pt -1 -abs 0 -d2r refFlat_exonarray_rowidmap.txt -d2a refFlat_exonarray_transcriptidmap.txt -net HPRD_v6_mouse_enhance.txt -netA geneid2genename_nr.txt -netD 2\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	dTruncateLowerBound = -1e20; 
	dTruncateUpperBound = 1e20;
	strcpy(strTransform, "Identity");
	nNormalize = 0;
	dPostLowerBound = -1e20; 
	dPostUpperBound = 1e20;
	strcpy(strPostTransform, "Identity"); 
	nTakeAbsoluteValue = 0;
	strcpy(strNetworkPath, "NULL");
	strcpy(strNetworkAnnotationPath, "NULL");
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strRefLocusDatabasePath, argc[ni]);
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
			strcpy(strOutPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-v") == 0)
		{
			ni++;
			strcpy(strScorePath, argc[ni]);
			vOK = 1;
		}
		else if(strcmp(argc[ni], "-vl") == 0)
		{
			ni++;
			dTruncateLowerBound = atof(argc[ni]);
			vlOK = 1;
		}
		else if(strcmp(argc[ni], "-vu") == 0)
		{
			ni++;
			dTruncateUpperBound = atof(argc[ni]);
			vuOK = 1;
		}
		else if(strcmp(argc[ni], "-vt") == 0)
		{
			ni++;
			strcpy(strTransform, argc[ni]);
			vtOK = 1;
		}
		else if(strcmp(argc[ni], "-n") == 0)
		{
			ni++;
			nNormalize = atoi(argc[ni]);
			nOK = 1;
		}
		else if(strcmp(argc[ni], "-pl") == 0)
		{
			ni++;
			dPostLowerBound = atof(argc[ni]);
			plOK = 1;
		}
		else if(strcmp(argc[ni], "-pu") == 0)
		{
			ni++;
			dPostUpperBound = atof(argc[ni]);
			puOK = 1;
		}
		else if(strcmp(argc[ni], "-pt") == 0)
		{
			ni++;
			strcpy(strPostTransform, argc[ni]);
			ptOK = 1;
		}
		else if(strcmp(argc[ni], "-abs") == 0)
		{
			ni++;
			nTakeAbsoluteValue = atoi(argc[ni]);
			absOK = 1;
		}
		else if(strcmp(argc[ni], "-d2r") == 0)
		{
			ni++;
			strcpy(strRefToProbeRowIDPath, argc[ni]);
			d2rOK = 1;
		}
		else if(strcmp(argc[ni], "-d2a") == 0)
		{
			ni++;
			strcpy(strRefToProbeNamePath, argc[ni]);
			d2aOK = 1;
		}
		else if(strcmp(argc[ni], "-net") == 0)
		{
			ni++;
			strcpy(strNetworkPath, argc[ni]);
			netOK = 1;
		}
		else if(strcmp(argc[ni], "-netA") == 0)
		{
			ni++;
			strcpy(strNetworkAnnotationPath, argc[ni]);
			netAOK = 1;
		}
		else if(strcmp(argc[ni], "-netD") == 0)
		{
			ni++;
			nNetDepth = atoi(argc[ni]);
			netDOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((dOK == 0) || (sOK == 0) || (oOK == 0) || (vOK == 0) || (vlOK == 0)
		|| (vuOK == 0) || (vtOK == 0) || (nOK == 0) || (plOK == 0) || (puOK == 0)
		|| (ptOK == 0) || (absOK == 0) || (d2rOK == 0) || (d2aOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = RefGene_Database_AssignValues_Main(strRefLocusDatabasePath, 
			strSpecies, strOutPath, strScorePath, nNormalize,
			dTruncateLowerBound, dTruncateUpperBound, strTransform, 
			dPostLowerBound, dPostUpperBound, strPostTransform, nTakeAbsoluteValue,
			strRefToProbeRowIDPath, strRefToProbeNamePath,
			strNetworkPath, strNetworkAnnotationPath, nNetDepth);
	}

	/* return */
	return nResult;
}