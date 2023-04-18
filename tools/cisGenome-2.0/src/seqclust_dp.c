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
#include "TilingArrayLib.h"
#include "HTSequencingLib.h"

int menu_seqclust_dp(int argv, char **argc);

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
	menu_seqclust_dp(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_seqclust_dp(int argv, char **argc)
{
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strOutputFile[LINE_LENGTH];
	int nBinSize = 25;
	int nKernelType = 1;
	int nKernelStep = 1;
	int nKernelLen = 300;
	int nKernelBand = 100;
	int nSegType = 0;
	char strSegFile[MED_LINE_LENGTH];
	int nDistType = 0;
	int nUp = 0;
	int nDown = 0;
	int nDatabaseType = 1;
	int nMemBlock = 100000;
	int nCorrBlock = 10;
	int nGridNum = 10000;
	int nCutType = 0;
	double dCutL = 0.99;
	double dCutH = 0.999;
	int nBlockLenCut = 100;
	int nExportBAR = 0;


	int iOK;
	int oOK;
	int dOK;
	int bOK;
	int ksOK;
	int siOK;
	int rOK;
	int upOK;
	int downOK;
	int dtOK;
	int ctOK;
	int clOK;
	int chOK;
	int ni;
	int nResult;

	/* ------------------------------- */
	/*     menu_seqclust_dp            */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    seqclust_dp                \n");
		printf(" -i input \n");
		printf(" -d output folder  \n");
		printf(" -o output file  \n");
		printf(" -b bin size (default = 25 bp) \n");
		printf(" -k smoothing kernel \n");
		printf("    0: no kernel, use raw read count \n");
		printf("    1 (default): one-sided exponential \n");
		printf("    2: two-sided exponential \n"); 
		printf("    3: one-sided gaussian \n");
		printf("    4: two-sided gaussian \n");
		printf(" -ks step size for kernel (a positive integer, default = bin_size) \n");
		printf(" -kw window size for kernel smoothing (a positive integer, default = 300, \n");
		printf("     which means smoothing is only carried out for positions within 300 base \n");
		printf("     pairs if one uses one-sided kernel, or within 300 bp each side if two-sided \n");
		printf("     kernel is used (i.e. 600 bp in total).\n");
		printf(" -kb kernel bandwidth (a positive integer, default = 100 (bp)) \n");
		printf(" -st approach for genome segmentation         \n");
		printf("     0 (default): automatic segmentation (no need to specify -si)      \n");
		printf("     1: user supplied coordinates (*.COD) file for genomic intervals (-si needs to be specified) \n");
		printf("     2: intervals defined based on gene structures (-si is needed; -r, -up, -down, -dt are optional) \n");
		printf(" -si if st=1, si specifies a COD file for genomic intervals \n");
		printf("     if st=2, si specifies a gene annotation file (refFlat_sorted.txt from CisGenome website) \n");
		printf(" -r if st=2, r specifies the reference points to define genomic intervals \n");
		printf("     0: TSS-up, TES-down (default)    \n");
		printf("     1: TSS-up, TSS-down     \n");
		printf("     2: TES-up, TES-down     \n");
		printf("     3: CDSS-up, CDSE-down   \n");
		printf("     4: CDSS-up, CDSS-down   \n");
		printf("     5: CDSE-up, CDSE-down   \n");
		printf(" -up if st=2, up specifies the distance 5' to the first reference point (default = 0) \n");
		printf(" -down if st=2, down specifies the distance 3' to the second reference point (default = 0) \n");
		printf("     Example: -st 2 -si refFlat_sorted.txt -r 0 -up 1000 -down 500 means get intervals \n");
		printf("	   that start from 1000 bp 5' upstream of TSS to 500 bp 3' downstream of TES of    \n");
		printf("       genes specified in the file refFlat_sorted.txt \n");
		printf(" -dt optional if st=2, annotation type: 0=refGene; 1=refFlat (default); 2=refLocus \n");
		printf(" -mb memory block size (default = 100000, which means each processing cycle will handle 100000 bins simultaneously \n");
		printf(" -cb correlation block size (default = 10, which means correlation coef of a bin is determined using its corr with 20 flanking bins (10 from left and 10 from right) \n");
		printf(" -ct cuttype. 0 (default): cut by percentiles of bin statistics; 1: cut by FDR; 2: cut by user supplied cutoffs. \n");
		printf(" -cl lower cutoff, default: 0.99 if cuttype = 0; 0.25 if cuttype = 1. This is the cutoff for defining peak boundaries. \n");
		printf(" -ch higher cutoff, default: 0.999 if cuttype = 0; 0.05 if cuttype = 1. This is the cutoff for initiating peaks. \n");
		printf(" -grid Number of grid points for computing FDR or percentile cutoffs (default = 10000). \n");
		printf(" -bar 1: export bin counts in bar file (may take significantly more time, but bar file can be used for visualization); 0: do not export bar file (default = 0). \n");
		printf(" -minlen minimal region size in base pairs, regions shorter than this size will not be reported\n");
		printf(" Example: \n");
		printf("    seqclust_dp -i input_filelist.txt -d /workingfolder -o mycluster\n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;
	dOK = 0;
	bOK = 0;
	ksOK = 0;
	siOK = 0;
	rOK = 0;
	upOK = 0;
	downOK = 0;
	dtOK = 0;
	ctOK = 0;
	chOK = 0;
	clOK = 0;
	strcpy(strSegFile, "");
	
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
		else if(strcmp(argc[ni], "-k") == 0)
		{
			ni++;
			nKernelType = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-ks") == 0)
		{
			ni++;
			nKernelStep = atoi(argc[ni]);
			ksOK = 1;
		}
		else if(strcmp(argc[ni], "-kw") == 0)
		{
			ni++;
			nKernelLen = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-kb") == 0)
		{
			ni++;
			nKernelBand = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-st") == 0)
		{
			ni++;
			nSegType = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-si") == 0)
		{
			ni++;
			strcpy(strSegFile, argc[ni]);
			siOK = 1;
		}
		else if(strcmp(argc[ni], "-r") == 0)
		{
			ni++;
			nDistType = atoi(argc[ni]);
			rOK = 1;
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
		else if(strcmp(argc[ni], "-dt") == 0)
		{
			ni++;
			nDatabaseType = atoi(argc[ni]);
			dtOK = 1;
		}
		else if(strcmp(argc[ni], "-mb") == 0)
		{
			ni++;
			nMemBlock = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-cb") == 0)
		{
			ni++;
			nCorrBlock = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-ct") == 0)
		{
			ni++;
			nCutType = atoi(argc[ni]);
			ctOK = 1;
		}
		else if(strcmp(argc[ni], "-cl") == 0)
		{
			ni++;
			dCutL = atof(argc[ni]);
			clOK = 1;
		}
		else if(strcmp(argc[ni], "-ch") == 0)
		{
			ni++;
			dCutH = atof(argc[ni]);
			chOK = 1;
		}
		else if(strcmp(argc[ni], "-grid") == 0)
		{
			ni++;
			nGridNum = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-bar") == 0)
		{
			ni++;
			nExportBAR = atoi(argc[ni]);
		}
		else if(strcmp(argc[ni], "-minlen") == 0)
		{
			ni++;
			nBlockLenCut = atoi(argc[ni]);
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if(ksOK == 0)
	{
		nKernelStep = nBinSize;
		if(nKernelStep == 0)
			nKernelStep = 1;
	}
	if(nSegType != 0)
	{
		if(siOK == 0)
		{
			printf("Error: please use -si to specify a genomic interval (if -st 1) or gene annotation file (if -st 2)!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	if(nCutType == 0)
	{
		if(clOK == 0)
			dCutL = 0.99;
		if(chOK == 0)
			dCutH = 0.999;
	}
	else if(nCutType == 1)
	{
		if(clOK == 0)
			dCutL = 0.25;
		if(chOK == 0)
			dCutH = 0.05;
	}
	else
	{
		if((clOK == 0) || (chOK == 0))
		{
			printf("Error: please specify segmentation cutoffs using -cl and/or -ch options!\n");
			exit(EXIT_FAILURE);
		}
	}

	if((iOK == 0) || (oOK == 0) || (dOK == 0))
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = SeqClust_SegDP_Main(strInputPath, strOutputPath, strOutputFile, nBinSize,
			nKernelType, nKernelStep, nKernelLen, nKernelBand, 
			nSegType, strSegFile, nDistType, nUp, nDown, nDatabaseType, 
			nMemBlock, nCorrBlock, nGridNum, nCutType, dCutL, dCutH, nBlockLenCut, nExportBAR);
	}

	return nResult;
}
