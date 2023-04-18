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
#include "HTSequencingLib.h"

int menu_hts_twosample_enrich(int argv, char **argc);

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
	menu_hts_twosample_enrich(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_hts_twosample_enrich(int argv, char **argc)
{
	char strPosInputPath[MED_LINE_LENGTH];
	char strNegInputPath[MED_LINE_LENGTH];
	char strOutputFolder[MED_LINE_LENGTH];
	char strOutputTitle[MED_LINE_LENGTH];
	char strFDRPath[MED_LINE_LENGTH];
	int nOneSide = 1;
	int nW = 100;
	int nS = 25;
	double dFDRCut = 0.1;
	int nMinLen = 1;
	int nMaxGap = 0;
	int nTCut = 10;
	double dP0 = 0.5;
	
	int nResult;
	int iOK;
	int nOK;
	int dOK;
	int oOK;
	int tOK;
	int wOK;
	int sOK;
	int fOK;
	int cOK;
	int mOK;
	int pOK;
	int gOK;
	int lOK;
	int ni;

	  
	/* ------------------------------- */
	/*    hts_peakdetector_2sample     */
	/* -i input (positive)             */
	/* -n input (negative)             */
	/* -d output folder                */
	/* -o output title                 */
	/* -t one-sided or two-sided test  */
	/* -w window size (default = 100)  */
	/* -s step size (default = 25)     */
	/* -f FDR file                     */
	/* -c fdr cutoff (default = 0.1)   */
	/* -m min window read count        */
	/* -p p0 = pos/neg background ratio*/
	/* -g max gap (default = 0)        */
	/* -l min region length (default = 1) */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("     hts_peakdetector_2sample      \n");
		printf(" -i input (positive)          \n");
		printf(" -n input (negative)          \n");
		printf(" -d output folder             \n");
		printf(" -o output title              \n");
		printf(" -t 1: one-sided or 0: two-sided test (default=1) \n");
		printf(" -w window size (default = 100)  \n");
		printf(" -s step size (default = 25)     \n");
		printf(" -f FDR file                     \n");
		printf(" -c fdr cutoff (default = 0.1)   \n");
		printf(" -m min window read_num (default = 10) \n");
		printf(" -p p0 = pos/neg background ratio \n");
		printf(" -g max gap (default = 0)        \n");
		printf(" -l min region length (default = 1) \n");
		printf(" example: \n");
		printf("    hts_peakdetector_2sample -i pos.bar -n neg.bar -d . -o output -f FDR.txt -c 0.2 -p 0.5\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	nOK = 0;
	dOK = 0;
	oOK = 0;
	tOK = 0;
	wOK = 0;
	sOK = 0;
	fOK = 0;
	cOK = 0;
	mOK = 0;
	pOK = 0;
	gOK = 0;
	lOK = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strPosInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-n") == 0)
		{
			ni++;
			strcpy(strNegInputPath, argc[ni]);
			nOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strOutputFolder, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputTitle, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-t") == 0)
		{
			ni++;
			nOneSide = atoi(argc[ni]);
			tOK = 1;

			if(nOneSide != 0)
				nOneSide = 1;
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
			nS = atoi(argc[ni]);
			sOK = 1;
		}
		else if(strcmp(argc[ni], "-f") == 0)
		{
			ni++;
			strcpy(strFDRPath, argc[ni]);
			fOK = 1;
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			dFDRCut = atof(argc[ni]);
			cOK = 1;
		}
		else if(strcmp(argc[ni], "-m") == 0)
		{
			ni++;
			nTCut = atoi(argc[ni]);
			mOK = 1;
		}
		else if(strcmp(argc[ni], "-p") == 0)
		{
			ni++;
			dP0 = atof(argc[ni]);
			pOK = 1;

			if( (dP0 <= 0.0) || (dP0 >= 1.0) )
			{
				printf("Error: -p option must be set to a value between 0 and 1 (i.e. 0<p<1)!\n");
				exit(EXIT_FAILURE);
			}
		}
		else if(strcmp(argc[ni], "-g") == 0)
		{
			ni++;
			nMaxGap = atoi(argc[ni]);
			gOK = 1;
		}
		else if(strcmp(argc[ni], "-l") == 0)
		{
			ni++;
			nMinLen = atoi(argc[ni]);
			lOK = 1;
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (nOK == 0) || (dOK == 0) || (oOK == 0) || (fOK == 0) || (pOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		if(nW <= 0)
		{
			nW = 100;
			printf("Warning: Window size<=0! Proceed with default window size = 100\n");
		}
		if(nS <= 0)
		{
			nS = 25;
			printf("Warning: Step size<=0! Proceed with default step size = 25\n");
		}
		if(nTCut <= 0)
		{
			nTCut = 10;
			printf("Warning: Window read_num cutoff <= 0! Proceed with default window read_num cutoff = 10\n");
		}
		if(dFDRCut <= 0.0)
		{
			dFDRCut = 0.1;
			printf("Warning: FDR cutoff <= 0! Proceed with default fdr cutoff = 0.1\n");
		}
	

		nResult =  HTS_Enrich_TwoSample_Main(strPosInputPath, strNegInputPath, 
					  nOneSide, nW, nS, nTCut, strFDRPath, dFDRCut,
					  nMinLen, nMaxGap, dP0,
					  strOutputFolder, strOutputTitle);
	}

	return nResult;
}