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

int menu_seqpeak(int argv, char **argc);

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
	menu_seqpeak(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_seqpeak(int argv, char **argc)
{
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strOutputFile[LINE_LENGTH];
	int nPaired = 0;
	int nBinSize = 50;
	int nExtLen = 150;
	int nSegType = 0;
	int nWinSize = 1;
	int nCutType = 0;
	double dCutoff = 3.0;
	int nMaxGap = 50;
	int nMinLen = 100;
	int nExportBAR = 0;
	int nKeepTempFiles = 0;
	int nBoundaryRefine = 1;
	int nBRWin = 5;
	int nCollectRawData = 0;
	int nPoisFilter = 1;
	int nPoisWin = 10000;
	double dPoisCut = 1.0e-5;
	int nTStandardize = 1;
	int nLFCAdj = 0;
	int nLFCWin = 10000;

	int iOK;
	int oOK;
	int dOK;

	int bOK;
	int wOK;
	int ni;
	int nResult;

	/* ------------------------------- */
	/*     menu_seqpeak                */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    seqpeak                \n");
		printf(" -i input \n");
		printf(" -d output folder  \n");
		printf(" -o output file  \n");
		printf(" -b bin size, the resolution of peak detection (default = 50 bp) \n");
		printf(" -w half window size, 2*w+1 bins will be used to determine enrichment of the center bin (default = 1) \n");
		printf(" -e read extension length  (default = 150 bp) \n");
		/* printf(" -st approach for genome segmentation         \n");
		printf("     0 (default): automatic segmentation      \n");
		printf(" -ct cuttype.(default = 0) \n");
		*/
		printf(" -ts 1 (default): standardize window statistics. 0: do not standardize. \n");
		printf(" -c cutoff for defining peak boundaries. (default = 3.0) \n");
		printf(" -maxgap (default = 50 bp): two peaks will be merged if the gap (in bp) between them is smaller than this value \n");
		printf(" -minlen minimal region size (default = 100 bp). Regions shorter than this size will not be reported\n");
		printf(" -br 1 (default): refine peak boundaries using +/- strand peak offsets; 0: no boundary refinement \n");
		printf(" -bw boundary refinement resolution (in bp). (default = 5). The smaller the number (i.e. the higher the resolution), the more memory is required. \n");
		printf(" -bar 1: export bin counts in bar file (may take significantly more time, but bar file can be used for visualization); 0 (default): do not export bar file. \n");
		printf(" -dat 1: report normalized read counts for each peak and each sample; 0 (default): not report the normalized raw data. \n");
		printf(" -keeptemp 1: keep intermediate files (may take big disk space); 0 (default): remove intermediate files. \n");
		printf(" -lpois 1 (default): apply local poisson rate filter (need to be careful when only IP sample is available and the signal is long); 0: no local poisson rate filter. \n");
		printf(" -lpwin window size (in bp) used to estimate local poisson rate (default = 10000). \n");
		printf(" -lpcut p-value cutoff for local poisson rate correction (default = 1e-5). \n");
		printf(" -lfc 1: local fold change adjustment (usually do not use it as it may reduce sensitivity); 0 (default): no adjustment. \n");
		printf(" -lfcwin window size (in bp) used to estimate average local coverage for lfc adjustment (default = 10000). \n");
		printf(" Example: \n");
		printf("    seqpeak -i input_filelist.txt -d . -o mypeak -b 25 -w 5 -e 100\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	dOK = 0;
	bOK = 0;
	wOK = 0;
	
	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-d") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			dOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputFile, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-b") == 0)
		{
			ni++;
			nBinSize = atoi(argc[ni]);
			bOK = 1;
		}
		else if(strcmp(argc[ni], "-w") == 0)
		{
			ni++;
			nWinSize = atoi(argc[ni]);
			wOK = 1;
		}
		else if(strcmp(argc[ni], "-e") == 0)
		{
			ni++;
			nExtLen = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-st") == 0)
		{
			ni++;
			nSegType = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-ct") == 0)
		{
			ni++;
			nCutType = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-ts") == 0)
		{
			ni++;
			nTStandardize = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-c") == 0)
		{
			ni++;
			dCutoff = atof(argc[ni]);
		}
		else if(strcmp(argc[ni], "-br") == 0)
		{
			ni++;
			nBoundaryRefine = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-bw") == 0)
		{
			ni++;
			nBRWin = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-bar") == 0)
		{
			ni++;
			nExportBAR = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-keeptemp") == 0)
		{
			ni++;
			nKeepTempFiles = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-maxgap") == 0)
		{
			ni++;
			nMaxGap = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-minlen") == 0)
		{
			ni++;
			nMinLen = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-dat") == 0)
		{
			ni++;
			nCollectRawData = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-lpois") == 0)
		{
			ni++;
			nPoisFilter = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-lpwin") == 0)
		{
			ni++;
			nPoisWin = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-lpcut") == 0)
		{
			ni++;
			dPoisCut = atof(argc[ni]);
		}
		else if(strcmp(argc[ni], "-lfc") == 0)
		{
			ni++;
			nLFCAdj = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-lfcwin") == 0)
		{
			ni++;
			nLFCWin = atoi(argc[ni]);
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if(nBinSize < 0)
		nBinSize = 50;

	if(nWinSize < 0)
		nWinSize = 0;
	
	if((iOK == 0) || (oOK == 0) || (dOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = SeqPeak_Main(strInputPath, strOutputPath, strOutputFile, 
			nPaired, nBinSize, nExtLen, 
			nSegType, nWinSize, nCutType, 
			nTStandardize, dCutoff, nMaxGap, nMinLen, 
			nBoundaryRefine, nBRWin,
			nExportBAR, nKeepTempFiles, nCollectRawData, 
			nPoisFilter, nPoisWin, dPoisCut, 
			nLFCAdj, nLFCWin);
	}

	return nResult;
}

