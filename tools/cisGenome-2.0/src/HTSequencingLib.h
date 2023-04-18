/* ----------------------------------------------------------------------- */
/*                                 Macro                                   */
/* ----------------------------------------------------------------------- */


/* ----------------------------------------------------------------------- */ 
/*                                 Struct                                  */
/* ----------------------------------------------------------------------- */ 

struct tagProbeGenomeInfo;

/* HTS aln2window paramters */
struct tagHTSAln2WinParam
{
	/* project name */
	char strProjectName[MED_LINE_LENGTH];

	/* chromosome list */
	char strChrList[MED_LINE_LENGTH];

	/* chromosome length */
	char strChrLen[MED_LINE_LENGTH];

	/* data directory */
	char strDataFolder[MED_LINE_LENGTH];

	/* No. of Libraries, Samples, Arrays and Groups */
	int nSampleNum;
	int nGroupNum;
	
	/* Sample name */
	struct INTMATRIX *vGroupLabel;
	struct tagString **vSampleAlias;
	struct tagString **vSampleFile;

	/* scaling parameter */
	double dScaling;

	/* window size */
	int nW;

	/* step size */
	int nS;

	/* export folder */
	char strExportFolder[MED_LINE_LENGTH];
};

/* node for hierarchical clustering */
struct tagHCNode
{
	/* node id */
	int nNodeId;

	/* child number */
	int nChildNum;

	/* child vector */
	struct tagHCNode **vChildNode;

	/* leaf number */
	int nLeafNum;

	/* leaf vector */
	int *vLeafId;

	/* height of the node */
	int nHeight;
};

/* ----------------------------------------------------------------------- */
/*                              Declaration                                */
/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Unique()                                                       */
/*  Take a sorted aln file as input, remove redundant reads.               */
/*  At each position, keep at most nMax reads.                             */
/* ----------------------------------------------------------------------- */ 
int HTS_Aln2Unique(char strInputPath[], char strOutputPath[], int nMax);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Window_Main()                                                  */
/*  Convert high throughput sequencing alignment to windowed bar tiling    */
/*  data. Each sample will be scaled to have the same number of aligned    */
/*  reads if specified.                                                    */
/* ----------------------------------------------------------------------- */ 
int HTS_Aln2Window_Main(char strParamPath[]);

/* ----------------------------------------------------------------------- */ 
/*  HTSAln2WinParam_Create()                                               */
/*  Create hts_aln2window paramter.                                        */
/* ----------------------------------------------------------------------- */ 
struct tagHTSAln2WinParam *HTSAln2WinParam_Create();

/* ----------------------------------------------------------------------- */ 
/*  HTSAln2WinParam_Destroy()                                              */
/*  Destroy hts_aln2window paramter.                                       */
/* ----------------------------------------------------------------------- */ 
void HTSAln2WinParam_Destroy(struct tagHTSAln2WinParam **pParam);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Window_LoadParamter()                                          */
/*  Load hts_aln2window paramter.                                          */
/* ----------------------------------------------------------------------- */ 
struct tagHTSAln2WinParam *HTS_Aln2Window_LoadParamter(char strParamPath[]);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Window_PrepareBARData()                                        */
/*  Prepare the bar data object.                                           */
/* ----------------------------------------------------------------------- */ 
struct tagBARData * HTS_Aln2Window_PrepareBARData(struct tagHTSAln2WinParam *pParam, 
				struct INTMATRIX *pChrLen);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Window_PrepareBARData()                                        */
/*  Count hits in each window.                                             */
/* ----------------------------------------------------------------------- */ 
double HTS_Aln2Window_CountHits(struct tagBARData *pBARData, char strInFile[], int nW);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Window_Scaling()                                               */
/*  Scaling HTS data.                                                      */
/* ----------------------------------------------------------------------- */ 
int HTS_Aln2Window_Scaling(struct tagBARData *pBARData, double dScale);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Window_WriteCGW()                                              */
/*  Write to CGW file.                                                     */
/* ----------------------------------------------------------------------- */ 
int HTS_Aln2Window_WriteCGW(struct tagHTSAln2WinParam *pParam);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Window_WriteCGB()                                              */
/*  Write to CGB file.                                                     */
/* ----------------------------------------------------------------------- */
int HTS_Aln2Window_WriteCGB(struct tagHTSAln2WinParam *pParam);


/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Diff_Main()                                                    */
/*  Detect differentially aligned regions.                                 */
/* ----------------------------------------------------------------------- */ 
int HTS_Aln2Diff_Main(char strParamPath[]);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Diff_PrepareBARData()                                          */
/*  Prepare the bar data object.                                           */
/* ----------------------------------------------------------------------- */ 
struct tagBARData * HTS_Aln2Diff_PrepareBARData(struct tagHTSAln2WinParam *pParam, 
				struct INTMATRIX *pChrLen);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Diff_PrepareRecountingData()                                   */
/*  Prepare the bar data object.                                           */
/* ----------------------------------------------------------------------- */ 
struct tagBARData * HTS_Aln2Diff_PrepareRecountingData(struct tagHTSAln2WinParam *pParam, 
				struct INTMATRIX *pChrLen);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Diff_CountHits()                                               */
/*  Count hits in each window.                                             */
/* ----------------------------------------------------------------------- */ 
double HTS_Aln2Diff_CountHits(struct tagBARData *pBARData, char strInFile[], int nW, int nSampleID);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Diff_EstimateP()                                               */
/*  Estimate hyperparameter.                                               */
/* ----------------------------------------------------------------------- */ 
int HTS_Aln2Diff_EstimateP(struct tagBARData *pBARData, double *vP, int nSampleNum);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Diff_EstimatePrior2Col()                                       */
/*  Estimate prior probability of difference                               */
/* ----------------------------------------------------------------------- */ 
int HTS_Aln2Diff_EstimatePrior2Col(struct tagBARData *pBARData, int nSampleNum, 
								   double *vP, double *dPriorP);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Diff_EstimatePBinom()                                          */
/*  Estimate hyperparamter for mixture binominal                           */
/* ----------------------------------------------------------------------- */ 
int HTS_Aln2Diff_EstimatePBinom(struct tagBARData *pBARData, int nSampleNum, 
								   double *vP);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Diff_EstimateP3()                                              */
/*  Estimate hyperparamter for mixture binominal                           */
/* ----------------------------------------------------------------------- */ 
int HTS_Aln2Diff_EstimateP3(struct tagBARData *pBARData, int nSampleNum, 
								   double *vP);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Diff_CallRegion()                                              */
/*  Call differentially expressed regions.                                 */
/* ----------------------------------------------------------------------- */ 
int HTS_Aln2Diff_CallRegion(struct tagHTSAln2WinParam *pParam, double *vP);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Enrich_Main()                                                  */
/*  Detect enriched alignment regions.                                     */
/* ----------------------------------------------------------------------- */ 
int HTS_Aln2Enrich_Main(char strParamPath[]);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Enrich_EstimateP()                                             */
/*  Estimate hyperparamter for mixture binominal                           */
/* ----------------------------------------------------------------------- */ 
int HTS_Aln2Enrich_EstimateP(struct tagBARData *pBARData, int nSampleNum, 
								   double *vP);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2BAR()                                                          */
/*  Convert high throughput sequencing alignment to bar file.              */
/* ----------------------------------------------------------------------- */ 
int HTS_Aln2BAR(char strTXTFile[], char strBARFile[]);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2BARv2()                                                        */
/*  Convert high throughput sequencing alignment to bar file.              */
/* ----------------------------------------------------------------------- */ 
int HTS_Aln2BARv2(char strTXTFile[], char strBARFile[]);

/* ----------------------------------------------------------------------- */ 
/*  HTS_WindowSummary()                                                    */
/*  Summarize window counts.                                               */
/* ----------------------------------------------------------------------- */ 
int HTS_WindowSummary(char strBARFile[], char strChrList[], char strChrLen[], 
					  int nW, char strOutFile[], struct tagBARData *pRepeatData);

/* ----------------------------------------------------------------------- */ 
/*  HTS_WindowSummaryv2()                                                  */
/*  Summarize window counts.                                               */
/* ----------------------------------------------------------------------- */ 
int HTS_WindowSummaryv2(char strBARFile[], char strChrList[], char strChrLen[], 
					  int nW, char strOutFile[], int nCombineShift, char strRepeatFile[]);

/* ----------------------------------------------------------------------- */ 
/*  HTS_WindowSummaryPaper()                                               */
/*  Summarize window counts.                                               */
/* ----------------------------------------------------------------------- */ 
int HTS_WindowSummaryPaper(char strBARFile[], char strChrList[], char strChrLen[], 
					  int nW, char strOutFile[], int nCombineShift);


/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_OneSample_Main()                                            */
/*  Find enriched region from a sequencing data set.                       */
/* ----------------------------------------------------------------------- */ 
int HTS_Enrich_OneSample_Main(char strBARFile[], int nW, int nS, int nCutoff,
					  int nMinLen, int nMaxGap,
					  char strExportFolder[], char strOutFileTitle[]);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_OneSamplev2_Main()                                          */
/*  Find enriched region from a sequencing data set.                       */
/* ----------------------------------------------------------------------- */ 
int HTS_Enrich_OneSamplev2_Main_Old(char strBARFile[], int nW, int nS, int nCutoff,
					  int nCutoffF, int nCutoffR, int nMinLen, int nMaxGap,
					  char strExportFolder[], char strOutFileTitle[]);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_OneSamplev2_Main()                                          */
/*  Find enriched region from a sequencing data set.                       */
/* ----------------------------------------------------------------------- */ 
int HTS_Enrich_OneSamplev2_Main(char strBARFile[], int nW, int nS, int nCutoff, 
					  int nCutoffF, int nCutoffR, int nMinLen, int nMaxGap, 
					  char strExportFolder[], char strOutFileTitle[],
					  int nBR, int nBRL, int nSSF, int nCombineShift);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_OneSample_CallRegion_Initial()                              */
/*  Search for enriched windows.                                           */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *HTS_Enrich_OneSample_CallRegion_Initial(struct tagBARData *pBARData,
			int nW, int nS, int nCutoff, char strExportFolder[], char strOutFileTitle[]);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_OneSample_MergeRegion()                                     */
/*  Merge and filter regions.                                              */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *HTS_Enrich_OneSample_MergeRegion(struct DOUBLEMATRIX *pRegion0, 
	int nMinLen, int nMaxGap);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_OneSample_MergeRegion()                                     */
/*  Merge and filter regions.                                              */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *HTS_Enrich_OneSamplev2_MatchFRRegions(struct DOUBLEMATRIX *pRegion, 
		struct DOUBLEMATRIX *pRegionF, struct DOUBLEMATRIX *pRegionR);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_OneSamplev2_FindOverlap()                                   */
/*  Find overlap region.                                                   */
/* ----------------------------------------------------------------------- */ 
int HTS_Enrich_OneSamplev2_FindOverlap(struct DOUBLEMATRIX *pRegion, 
					int nChr, int nStart, int nEnd);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_OneSample_ExportResults()                                   */
/*  Export results.                                                        */
/* ----------------------------------------------------------------------- */ 
int	HTS_Enrich_OneSample_ExportResults(struct DOUBLEMATRIX *pRegion, 
							struct tagBARData *pBARData, char strFileName[]);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_OneSamplev2_ExportResults()                                 */
/*  Export results.                                                        */
/* ----------------------------------------------------------------------- */ 
int	HTS_Enrich_OneSamplev2_ExportResults_Old(struct DOUBLEMATRIX *pRegion, 
				struct tagBARData *pBARData, char strFileName[]);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_OneSamplev2_ExportResults()                                 */
/*  Export results.                                                        */
/* ----------------------------------------------------------------------- */ 
int	HTS_Enrich_OneSamplev2_ExportResults(struct DOUBLEMATRIX *pRegion, 
				struct tagBARData *pBARData, char strFileName[], int nBR, int nBRL,
				int nSSF, int nSSFF, int nSSFR);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_OneSample_RegionCollectInfo()                               */
/*  Collect reads number.                                                  */
/* ----------------------------------------------------------------------- */ 
int HTS_Enrich_OneSample_RegionCollectInfo(struct DOUBLEMATRIX *pRegion, struct tagBARData *pBARData);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_OneSamplev2_RegionCollectInfo()                             */
/*  Collect forward/reverse reads info                                     */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *HTS_Enrich_OneSamplev2_RegionCollectInfo(struct DOUBLEMATRIX *pRegion, 
				struct tagBARData *pBARData, int nW);

/* ----------------------------------------------------------------------- */ 
/*  HTS_TwoSample_WindowSummary()                                          */
/*  Summarize window counts.                                               */
/* ----------------------------------------------------------------------- */ 
int HTS_TwoSample_WindowSummary(char strPosBARFile[], char strNegBARFile[],
					char strChrList[], char strChrLen[], 
					int nW, char strOutFile[]);

/* ----------------------------------------------------------------------- */ 
/*  HTS_TwoSample_WindowSummaryv2()                                        */
/*  Summarize window counts.                                               */
/* ----------------------------------------------------------------------- */ 
int HTS_TwoSample_WindowSummaryv2(char strPosBARFile[], char strNegBARFile[],
					char strChrList[], char strChrLen[], 
					int nW, char strOutFile[], int nCombineShift);

/* ----------------------------------------------------------------------- */ 
/*  HTS_TwoSample_FDR()                                                    */
/*  Compute FDR for two sample comparison.                                 */
/* ----------------------------------------------------------------------- */ 
int HTS_TwoSample_FDR(struct DOUBLEMATRIX *pCount, struct DOUBLEMATRIX *pCount2D, 
					  int nMaxC, int nMaxCPos, struct DOUBLEMATRIX **pFDR, 
					  double *dP0);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_TwoSample_Main()                                            */
/*  Find enriched region from a sequencing data set.                       */
/* ----------------------------------------------------------------------- */ 
int HTS_Enrich_TwoSample_Main(char strPosBARFile[], char strNegBARFile[], 
					  int nOneSide, int nW, int nS, int nTCut, char strFDRFile[], double dFDRCut,
					  int nMinLen, int nMaxGap, double dP0,
					  char strExportFolder[], char strOutFileTitle[]);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_TwoSamplev2_Main()                                          */
/*  Find enriched region from a sequencing data set.                       */
/* ----------------------------------------------------------------------- */ 
int HTS_Enrich_TwoSamplev2_Main(char strPosBARFile[], char strNegBARFile[], 
					  int nOneSide, int nW, int nS, int nTCut, char strFDRFile[], double dFDRCut,
					  int nMinLen, int nMaxGap, double dP0,
					  char strExportFolder[], char strOutFileTitle[],
					  int nBR, int nBRL, int nSSF, int nSSFF, int nSSFR, int nCombineShift,
					  double dFC, double dTFC);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_TwoSample_CallRegion_Initial()                              */
/*  Search for differentially expressed windows.                           */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *HTS_Enrich_TwoSample_CallRegion_Initial(struct tagBARData *pBARDataPos,
			struct tagBARData *pBARDataNeg, int nW, int nS, int nTCut,
			struct DOUBLEMATRIX *pFDR, double dFDRCut, double dP0, int nOneSide,
			char strExportFolder[], char strOutFileTitle[]);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_TwoSample_MergeRegion()                                     */
/*  Merge and filter regions.                                              */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *HTS_Enrich_TwoSample_MergeRegion(struct DOUBLEMATRIX *pRegion0, 
	int nMinLen, int nMaxGap);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_TwoSample_ExportResults()                                   */
/*  Export results.                                                        */
/* ----------------------------------------------------------------------- */ 
int	HTS_Enrich_TwoSample_ExportResults(struct DOUBLEMATRIX *pRegion, 
			struct tagBARData *pBARDataPos, struct tagBARData *pBARDataNeg,
			char strFileName[]);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_TwoSamplev2_ExportResults()                                 */
/*  Export results.                                                        */
/* ----------------------------------------------------------------------- */ 
int	HTS_Enrich_TwoSamplev2_ExportResults(struct DOUBLEMATRIX *pRegion, 
			struct tagBARData *pBARDataPos, struct tagBARData *pBARDataNeg,
			char strFileName[], int nBR, int nBRL, int nSSF, int nSSFF, int nSSFR,
			double dFC, double dTFC);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_TwoSample_RegionCollectInfo()                               */
/*  Collect reads number.                                                  */
/* ----------------------------------------------------------------------- */ 
int HTS_Enrich_TwoSample_RegionCollectInfo(struct DOUBLEMATRIX *pRegion, 
		struct tagBARData *pBARDataPos, struct tagBARData *pBARDataNeg);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_TwoSamplev2_RegionCollectInfo()                             */
/*  Collect reads number.                                                  */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *HTS_Enrich_TwoSamplev2_RegionCollectInfo(struct DOUBLEMATRIX *pRegion, 
		struct tagBARData *pBARDataPos, struct tagBARData *pBARDataNeg, int nW,
		struct DOUBLEMATRIX *pFDR, double dP0, int nOneSide);

/* ----------------------------------------------------------------------- */ 
/*  HTS_AlnShift2BAR()                                                     */
/*  shift base pairs.                                                      */
/* ----------------------------------------------------------------------- */ 
int HTS_AlnShift2BAR(char strInputPath[], char strOutputPath[], int nS);

/* ----------------------------------------------------------------------- */ 
/*  HTS_CreateRepeatFilter()                                               */
/*  Create repeat filter.                                                  */
/* ----------------------------------------------------------------------- */ 
int HTS_CreateRepeatFilter(char strGenomePath[], char strChrList[], char strChrLen[], 
					 int nW, double dR, char strExportFolder[], char strOutFileTitle[]);

/* ----------------------------------------------------------------------- */ 
/*  HTS_FilterRepeatReads()                                                */
/*  Filter repeat reads.                                                   */
/* ----------------------------------------------------------------------- */ 
int HTS_FilterRepeatReads(char strInputFile[], char strOutputFile[], char strGenomePath[]);

/* ----------------------------------------------------------------------- */ 
/*  HTS_CollectReads()                                                     */
/*  Collecte reads in specific regions.                                    */
/* ----------------------------------------------------------------------- */ 
int HTS_CollectReads(char strInputFile[], char strOutputFile[], char strBARFile[]);

/* ----------------------------------------------------------------------- */ 
/*  HTS_CollectProfile_Main()                                              */
/*  Collecte read profile in specific regions.                             */
/* ----------------------------------------------------------------------- */ 
int HTS_CollectProfile_Main(char strInputFile[], char strOutputFile[], char strBARFile[], int nW, int nS);

/* ----------------------------------------------------------------------- */ 
/*  HTS_CollectProfile()                                                   */
/*  Collecte read profile in specific regions.                             */
/* ----------------------------------------------------------------------- */ 
struct INTMATRIX* HTS_CollectProfile(char strChr[], int nPos, int nW, int nS, struct tagBARData *pBARData);

/* ----------------------------------------------------------------------- */ 
/*  CNV_Aln2Window_Main()                                                  */
/*  Convert copy number sequencing alignment to windowed bar tiling        */
/*  data.                                                                  */
/* ----------------------------------------------------------------------- */ 
int CNV_Aln2Window_Main(char strInFile[], char strOutFile[], 
						char strChrListFile[], char strChrLenFile[], 
						int nW, int nL, int nN);

/* ----------------------------------------------------------------------- */ 
/*  CNV_Aln2WindowC_Main()                                                 */
/*  Convert copy number sequencing alignment to windowed bar tiling        */
/*  data.                                                                  */
/* ----------------------------------------------------------------------- */ 
int CNV_Aln2WindowC_Main(char strInFile[], char strOutFile[], 
						char strChrListFile[], char strChrLenFile[], 
						int nW, int nL, int nN);

/* ----------------------------------------------------------------------- */ 
/*  CNV_Aln2Window_PrepareBARData()                                        */
/*  Prepare the bar data object.                                           */
/* ----------------------------------------------------------------------- */ 
struct tagBARData * CNV_Aln2Window_PrepareBARData(int nW, struct INTMATRIX *pChrLen, 
				char strChrListFile[]);

/* ----------------------------------------------------------------------- */ 
/*  CNV_Repeat2Window_Main()                                               */
/*  Convert repeat percentage to windowed bar tiling data.                 */
/* ----------------------------------------------------------------------- */ 
int CNV_Repeat2Window_Main(char strGenomePath[], char strOutFile[], 
						char strChrListFile[], char strChrLenFile[], 
						int nW, int nL, int nN);

/* ----------------------------------------------------------------------- */ 
/*  RNASEQ_CountReadPerTranscript_Main()                                   */
/*  Count reads for all transcripts.                                       */
/*  nInputType = 0: input is a single file specified by strInFile.         */
/*  nInputType = 1: input is a list of files specified by strInFile.       */
/* ----------------------------------------------------------------------- */ 
int RNASEQ_CountReadPerTranscript_Main(char strInFile[], int nInputType,
						char strOutFile[], char strDatabaseFile[], 
						int nDatabaseType, char strSpecies[], int nStandardize);

/* ----------------------------------------------------------------------- */ 
/*  RNASEQ_HitTest()                                                       */
/*  Test if a read is aligned to an exon of a gene.                        */
/* ----------------------------------------------------------------------- */ 
int RNASEQ_HitTest(int nChr, int nPos, struct tagRefGene *pRefGene);

/* ----------------------------------------------------------------------- */ 
/*  HTS_CountReads4RefGene_Main()                                          */
/*  Count reads for all transcripts.                                       */
/*  nInputType = 0: input is a single file specified by strInFile.         */
/*  nInputType = 1: input is a list of files specified by strInFile.       */
/* ----------------------------------------------------------------------- */ 
int HTS_CountReads4RefGene_Main(char strDatabasePath[], int nDatabaseType,
			char strSpecies[], char strInputPath[], int nInputType, char strOutputPath[],
			int nRefType, int nUP, int nDOWN, int nStandardizebyT, int nStandardizebyL);

/* ----------------------------------------------------------------------- */ 
/*  HTS_CountReads4Region()                                                */
/*  Count reads for a region.                                              */
/* ----------------------------------------------------------------------- */ 
double HTS_CountReads4Region(char strChr[], int nStart, int nEnd, struct tagBARData *pBARData);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_Seg_Main()                                                    */
/*  Clustering analysis of ChIP-seq.                                       */
/* ----------------------------------------------------------------------- */ 
int SeqClust_Seg_Main(char strInputPath[], char strOutputPath[], char strOutputFile[],
				  int nBinSize, int nKernelType, int nKernelStep, 
				  int nKernelLen, int nKernelBand,
				  int nSegType, char strSegFile[], int nDistType, 
				  int nUp, int nDown, int nDatabaseType, int nMemBlock,
				  int nCorrBlock, int nGridNum,
				  int nCutType, double dCutL, double dCutH,
				  int nExportBAR);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_InitializeKernel()                                            */
/*  Clustering analysis of ChIP-seq: Initialize smoothing kernel vector    */
/* ----------------------------------------------------------------------- */ 
int SeqClust_InitializeKernel(struct DOUBLEMATRIX **ppSmoothK, int nKernelType, 
			int nKernelStep, int nKernelLen, int nKernelBand, int *pKStart, int *pKNum);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_Initialize()                                                  */
/*  Clustering analysis of ChIP-seq: Initialize                            */
/* ----------------------------------------------------------------------- */ 
int SeqClust_Initialize(char strInputPath[], int *pFileNum, struct tagString ***vFileName,
						struct DOUBLEMATRIX **ppFileReadCount, 
						struct DOUBLEMATRIX **ppFileNormFactor, int *pBaseLineFileId, 
						int *pChrNum, 
						struct tagString ***vvChrName, struct INTMATRIX **ppChrLen);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_ReadPos_From_Aln()                                            */
/*  Clustering analysis of ChIP-seq: Parse coordinates from aln file       */
/* ----------------------------------------------------------------------- */ 
int SeqClust_ReadPos_From_Aln(char strLine[], char strChr[], int *pPos, int *pStrand);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_ReadPos_From_Cod()                                            */
/*  Clustering analysis of ChIP-seq: Parse coordinates from cod file       */
/* ----------------------------------------------------------------------- */ 
int SeqClust_ReadPos_From_Cod(char strLine[], char strAlias[], char strChr[], 
							  int *pPos1, int *pPos2, int *pStrand);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_FindChr()                                                     */
/*  Clustering analysis of ChIP-seq: Find matching chromosome              */
/*  If no matching chromosome is found, the chromosome will be added to    */
/*  the chromosome list. If there are more than nMaxChrNum, return -1.     */
/*  otherwise return the index of the matching chromosome.                 */
/* ----------------------------------------------------------------------- */ 
int SeqClust_FindChrInfo(struct tagProbeGenomeInfo **vChrList, int nMaxChrNum, 
						 int *pChrNum, char strChr[]);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_CountBinReads()                                               */
/*  Count reads for genomic bins                                           */
/* ----------------------------------------------------------------------- */ 
int SeqClust_CountBinReads(char strInputPath[], char strOutputPath[],
		char strOutputFile[], int nFileNum, int nChrNum, 
		struct tagString **vChrName, struct INTMATRIX *pChrLen, int nBinSize,
		int nKernelType, struct DOUBLEMATRIX *pSmoothK, int nSmoothStart, 
		int nKernelStep, int nSmoothNum, struct BYTEMATRIX *pIncludeChr,
		int nExportBAR);


/* ----------------------------------------------------------------------- */ 
/*  SeqClust_CountBinReads_SingleFile()                                    */
/*  Count reads for genomic bins for a single sample file                  */
/* ----------------------------------------------------------------------- */ 
int SeqClust_CountBinReads_SingleFile(char strInputFile[], char strOutputPath[],
				char strOutputFile[], int nChrNum, struct tagString **vChrName, 
				struct INTMATRIX *pChrLen, int nBinSize,
				int nKernelType, struct DOUBLEMATRIX *pSmoothK, int nSmoothStart, 
				int nKernelStep, int nSmoothNum, struct BYTEMATRIX *pIncludeChr,
				int nExportBAR);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_FindChr()                                                     */
/*  Find matching chromosome                                               */
/*  If no matching chromosome is found, return -1.                         */
/*  otherwise return the index of the matching chromosome.                 */
/* ----------------------------------------------------------------------- */ 
int SeqClust_FindChr(struct tagString **vChrName, int nChrNum, char strChr[]);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_ComputeBinStats()                                             */
/*  Compute summary statistics for genomic bins.                           */
/* ----------------------------------------------------------------------- */ 
int SeqClust_ComputeBinStats(char strOutputPath[], char strOutputFile[], 
			int nFileNum, struct tagString **vFileName, struct DOUBLEMATRIX *pFileNormFactor, 
			int nChrNum, struct tagString **vChrName, struct INTMATRIX *pChrLen, 
			int nBinSize, int nMemBlock, int nCorrBlock, double *pBinMax, double *pBinMin);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_ComputeBinStats_Chr()                                         */
/*  Compute summary statistics for genomic bins for a single chromosome    */
/* ----------------------------------------------------------------------- */ 
int SeqClust_ComputeBinStats_Chr(char strOutputPath[], char strOutputFile[], 
			int nFileNum, struct tagString **vFileName, struct DOUBLEMATRIX *pFileNormFactor, 
			char strChr[], int nBinNum, int nBinSize, int nMemBlock, int nCorrBlock,
			double *pBinMax, double *pBinMin);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_QuantileCutoff()                                              */
/*  Find quantile cutoffs, save to dCL and dCH.                            */
/* ----------------------------------------------------------------------- */ 
int SeqClust_QuantileCutoff(char strOutputPath[], char strOutputFile[], 
		int nChrNum, struct tagString **vChrName, struct INTMATRIX *pChrLen, 
		int nBinSize, double dBinMax, double dBinMin, 
		int nGridNum, double dCutL, double dCutH, double *pCL, double *pCH,
		struct BYTEMATRIX *pIncludeChr);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_FDRCutoff()                                                   */
/*  Find FDR cutoffs, save to dCL and dCH.                                 */
/* ----------------------------------------------------------------------- */ 
int SeqClust_FDRCutoff(char strOutputPath[], char strOutputFile[], 
		int nChrNum, struct tagString **vChrName, struct INTMATRIX *pChrLen, 
		int nBinSize, double dBinMax, double dBinMin, 
		int nGridNum, double dCutL, double dCutH, double *pCL, double *pCH,
		struct BYTEMATRIX *pIncludeChr);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_Segmentation()                                                */
/*  Genome segmentation. Divide genome into correlation blocks.            */
/* ----------------------------------------------------------------------- */ 
int SeqClust_Segmentation(char strOutputPath[], char strOutputFile[], 
				int nChrNum, struct tagString **vChrName, struct INTMATRIX *pChrLen, 
				int nBinSize, double dCL, double dCH);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_CollectData_ForAutoSeg()                                      */
/*  collect data from automatically determined blocks for clustering.      */
/* ----------------------------------------------------------------------- */ 
int SeqClust_CollectData_ForAutoSeg(char strOutputPath[], char strOutputFile[], 
				int nFileNum, struct tagString **vFileName, 
				struct DOUBLEMATRIX *pFileNormFactor, double dSmoothFactor,
				int nChrNum, struct tagString **vChrName, 
				struct INTMATRIX *pChrLen, int nBinSize,
				struct BYTEMATRIX *pIncludeChr);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_CollectData_ForAutoSeg_SingleSample()                         */
/*  collect data from a single sample (file) for a region list.            */
/* ----------------------------------------------------------------------- */ 
int SeqClust_CollectData_ForAutoSeg_SingleSample(struct DOUBLEMATRIX *pReg, 
					double *pData, char strFileName[], int nBinNum);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_CollectData_RefineFileNormFactor()                            */
/*  Recompute file normalizing factor after excluding variable blocks.     */
/* ----------------------------------------------------------------------- */ 
int SeqClust_CollectData_RefineFileNormFactor(char strOutputPath[], char strOutputFile[], 
				int nFileNum, struct tagString **vFileName, int nBaseFileId,
				struct DOUBLEMATRIX *pFileReadCount,
				struct DOUBLEMATRIX *pFileNormFactor, double dSmoothFactor,
				int nChrNum, struct tagString **vChrName, 
				struct INTMATRIX *pChrLen, int nBinSize);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_Main()                                                        */
/*  Model Based Clustering of ChIP-seq.                                    */
/* ----------------------------------------------------------------------- */ 
int SeqClust_Main(char strInputPath[], char strOutputPath[], char strOutputFile[],
			int nSkipCol, int nKmin, int nKmax, int nKr, int nMethod, 
			int nBmin, int nBmax, int nSeed,
			int nTransform, double dTL, int nRowStandardize, double dCut,
			int nMaxIter, double dTol);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_FixedBlock()                                                  */
/*  Fixed covariance block clustering.                                     */
/* ----------------------------------------------------------------------- */ 
int SeqClust_FixedBlock(int nGeneNum, int nSampleNum, double **vData, 
			int nKmin, int nKmax, int nKr, int nMethod, 
			int nMaxIter, double dTol, 
			char strOutputPath[], char strOutputFile[],	
			int nGroupNum, int *vGroupSize, int **vGroupId,
			int *pBestK, int *pBestId);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_AdaptiveBlock()                                               */
/*  Adaptive covariance block clustering.                                  */
/* ----------------------------------------------------------------------- */ 
int SeqClust_AdaptiveBlock(int nGeneNum, int nSampleNum, double **vData, 
			int nKmin, int nKmax, int nKr, int nMethod, int nBmin, int nBmax,
			int nMaxIter, double dTol, 
			char strOutputPath[], char strOutputFile[],	
			int *pBestB, int *pBestK, int *pBestId);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_AdaptiveFactor()                                              */
/*  Adaptive mixture factor model based clustering.                        */
/* ----------------------------------------------------------------------- */ 
int SeqClust_AdaptiveFactor(int nGeneNum, int nSampleNum, double **vData, 
			int nKmin, int nKmax, int nKr, int nMethod, int nBmin, int nBmax,
			int nMaxIter, double dTol, 
			char strOutputPath[], char strOutputFile[],	
			int *pBestB, int *pBestK, int *pBestId);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_ReadCovBlock()                                                */
/*  Load covariance blocks for clustering analysis.                        */
/* ----------------------------------------------------------------------- */ 
int SeqClust_ReadCovBlock(char strFileName[], int *pGroupNum, 
						  int **vGroupSize, int ***vGroupId);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_LoadData_Kmeans()                                             */
/*  Load data for K-means clustering.                                      */
/* ----------------------------------------------------------------------- */ 
int SeqClust_LoadData_Kmeans(char strInputPath[], int nSkipCol, 
		int *pGeneNum, int *pSampleNum, 
		struct tagString ***vInfo, double ***vData, 
		int nTransform, double dTL, int nRowStandardize,
		int *pGroupNum, int **vGroupSize, int ***vGroupId);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_Knorm()                                                       */
/*  Model based clustering (K normal distributions with diagonal COV       */
/* ----------------------------------------------------------------------- */ 
double SeqClust_Knorm(int nGeneNum, int nSampleNum, double **vData, 
				int nK, int nMaxIter, double dTol, 
				char strOutputPath[], char strOutputFile[]);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_KMVnorm()                                                     */
/*  Model based clustering (K normal distributions with block covariance   */
/*  matrix.                                                                */
/* ----------------------------------------------------------------------- */ 
double SeqClust_KMVnorm(int nGeneNum, int nSampleNum, double **vData, 
				int nK, int nMaxIter, double dTol, 
				char strOutputPath[], char strOutputFile[],
				int nGroupNum, int *vGroupSize, int **vGroupId);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_KMFnorm()                                                     */
/*  Mixture factor model based clustering (K clusters with covariance      */
/*  matrix modeled by nF factors).                                         */
/* ----------------------------------------------------------------------- */ 
double SeqClust_KMFnorm(int nGeneNum, int nSampleNum, double **vData, 
				int nF, int nK, int nMaxIter, double dTol, 
				char strOutputPath[], char strOutputFile[]);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_KMVnorm_Quadrature()                                          */
/*  compute normal quadrature -(x-mu)'invV(x-mu)/2                         */
/* ----------------------------------------------------------------------- */ 
double SeqClust_KMVnorm_Quadrature(double *vX, double *vMu, 
					struct DOUBLEMATRIX *pinvV, int *vGroupId, int nGroupSize);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_KMFnorm_Quadrature()                                          */
/*  compute normal quadrature -(x-mu)'invV(x-mu)/2                         */
/* ----------------------------------------------------------------------- */ 
double SeqClust_KMFnorm_Quadrature(double *vX, double *vMu, struct DOUBLEMATRIX *pinvV);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_Kpolya()                                                      */
/*  Model based clustering (K polya distributions)                         */
/* ----------------------------------------------------------------------- */ 
double SeqClust_Kpolya(int nGeneNum, int nSampleNum, double **vData, 
				int nK, int nMaxIter, double dTol, 
				char strOutputPath[], char strOutputFile[]);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_Kpolya_MM()                                                   */
/*  MM algorithm for fitting polya distribution.                           */
/* ----------------------------------------------------------------------- */ 
int SeqClust_Kpolya_MM(int nSampleNum, int nMMax, int nXMax, 
					   double *vMHist, double **vXHist, 
					   double *vAold, double *vAnew);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_Kpolya_MM_LogLike()                                           */
/*  Compute Polya loglikelihood for MM algorithm                           */
/* ----------------------------------------------------------------------- */ 
double SeqClust_Kpolya_MM_LogLike(int nSampleNum, int nMMax, int nXMax, 
					   double *vMHist, double **vXHist, double *vA);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_Count_Main()                                                  */
/*  Get read count for genomic regions.                                    */
/* ----------------------------------------------------------------------- */ 
int SeqClust_Count_Main(char strInputPath[], char strDataPath[], char strOutputFile[], 
						int nExtendLen);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_Count_LoadRegion()                                            */
/*  Load input coordinates.                                                */
/* ----------------------------------------------------------------------- */ 
int SeqClust_Count_LoadRegion(char strInputPath[], struct DOUBLEMATRIX **pRegion,
			    struct tagString ***vvRegionName,
				struct DOUBLEMATRIX **pSortRegion, struct LONGMATRIX **pSortId,
				struct tagString ***vvChrName, int *pChrNum);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_Count_SingleSample()                                          */
/*  Count read for a single file.                                          */
/* ----------------------------------------------------------------------- */ 
int SeqClust_Count_SingleSample(char strFileName[], struct DOUBLEMATRIX *pRegion, 
							struct tagString **vChrName, int nChrNum, 
							double **vData, int nDataCol, int nExtendLen);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_Count_ExportCod()                                             */
/*  Save count data to output files.                                       */
/* ----------------------------------------------------------------------- */ 
int SeqClust_Count_ExportCod(char strOutputFile[], struct DOUBLEMATRIX *pRegion, 
							 struct tagString **vRegionName, 
							 struct tagString **vChrName, int nChrNum, 
							 double **vDataS, int nFileNum,
							 struct tagString **vFileName);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_Count_RegCompare()                                            */
/*  Compare two regions.                                                   */
/* ----------------------------------------------------------------------- */ 
int SeqClust_Count_RegCompare(int nChrId, int nStart, int nEnd, 
							  struct DOUBLEMATRIX *pRegion, int nRegionId);

/* ----------------------------------------------------------------------- */ 
/*  HierarchicalCluster_Column()                                           */
/*  Hierarchical clustering. Cluster columns of a data matrix.             */
/*  nDistanceType: 0=correlation, 1=absolute correlation.                  */
/*  nMergeType: 0=average linkage.                                         */
/*  nExportCluster: 0=save cluster to files                                */
/*  nMinClustNum, nMaxClustNum: cut so that there are min<=n<=max clusters */
/* ----------------------------------------------------------------------- */ 
int HierarchicalCluster_Column(int nRowNum, int nColNum, double **vData, 
		int nDistanceType, int nMergeType, 
		int nExportCluster, char strOutputPath[], char strOutputFile[],
		int nMinClustNum, int nMaxClustNum);

/* ----------------------------------------------------------------------- */ 
/*  HierarchicalCluster_Column_Dist_Corr()                                 */
/*  Compute correlation distance.                                          */
/* ----------------------------------------------------------------------- */ 
int HierarchicalCluster_Column_Dist_Corr(int nRowNum, int nColNum, 
					double **vData, double **vDist);

/* ----------------------------------------------------------------------- */ 
/*  HierarchicalCluster_Column_Dist_AbsCorr()                              */
/*  Compute absolute correlation distance.                                 */
/* ----------------------------------------------------------------------- */ 
int HierarchicalCluster_Column_Dist_AbsCorr(int nRowNum, int nColNum, 
					double **vData, double **vDist);

/* ----------------------------------------------------------------------- */ 
/*  HCNode_Create()                                                        */
/*  Create node for hierarchical clustering.                               */
/* ----------------------------------------------------------------------- */ 
struct tagHCNode *HCNode_Create(int nNodeId, int nHeight, int nChildNum, int nLeafNum);

/* ----------------------------------------------------------------------- */ 
/*  HCNode_Delete()                                                        */
/*  delete a node for hierarchical clustering.                             */
/* ----------------------------------------------------------------------- */ 
int HCNode_Delete(struct tagHCNode **pNode);

/* ----------------------------------------------------------------------- */ 
/*  HCNode_ClearTree()                                                     */
/*  delete all nodes in a tree.                                            */
/* ----------------------------------------------------------------------- */ 
int HCNode_ClearTree(struct tagHCNode **pNode);

/* ----------------------------------------------------------------------- */ 
/*  HierarchicalCluster_UpdateDist_AverageLinkage()                        */
/*  compute average linkage distance.        .                             */
/* ----------------------------------------------------------------------- */ 
double HierarchicalCluster_UpdateDist_AverageLinkage(struct tagHCNode *pNode1,
				struct tagHCNode *pNode2, double **vDist);

/* ----------------------------------------------------------------------- */ 
/*  HierarchicalCluster_UpdateDist_AverageLinkage()                        */
/*  compute average linkage distance.        .                             */
/* ----------------------------------------------------------------------- */ 
int HierarchicalCluster_ExportCluster(struct tagHCNode **vHCNode, 
					int nNodeNum, char strOutFile[]);

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_SegDP_Main()                                                  */
/*  Segmentation by Dynamic Programming.                                   */
/* ----------------------------------------------------------------------- */ 
int SeqClust_SegDP_Main(char strInputPath[], char strOutputPath[], char strOutputFile[],
				  int nBinSize, int nKernelType, int nKernelStep, 
				  int nKernelLen, int nKernelBand,
				  int nSegType, char strSegFile[], int nDistType, 
				  int nUp, int nDown, int nDatabaseType, 
				  int nMemBlock, int nCorrBlock, int nGridNum,
				  int nCutType, double dCutL, double dCutH,
				  int nBlockLenCut, int nExportBAR);

/* ----------------------------------------------------------------------- */ 
/*  SeqClustDP_ComputeBinStats()                                           */
/*  Compute summary statistics for genomic bins.                           */
/* ----------------------------------------------------------------------- */ 
int SeqClustDP_ComputeBinStats(char strOutputPath[], char strOutputFile[], 
			int nFileNum, struct tagString **vFileName, struct DOUBLEMATRIX *pFileNormFactor, 
			int nChrNum, struct tagString **vChrName, struct INTMATRIX *pChrLen,
			int nBinSize, int nMemBlock, int nCorrBlock, double *pBinMax, double *pBinMin,
			int nExportBAR);

/* ----------------------------------------------------------------------- */ 
/*  SeqClustDP_ComputeBinStats_Chr()                                       */
/*  Compute summary statistics for genomic bins for a single chromosome    */
/* ----------------------------------------------------------------------- */ 
int SeqClustDP_ComputeBinStats_Chr(char strOutputPath[], char strOutputFile[], 
			int nFileNum, struct tagString **vFileName, struct DOUBLEMATRIX *pFileNormFactor, 
			char strChr[], int nBinNum, int nBinSize, int nMemBlock, int nCorrBlock,
			double *pBinMax, double *pBinMin, int nExportBAR);

/* ----------------------------------------------------------------------- */ 
/*  SeqClustDP_QuantileCutoff()                                            */
/*  Find quantile cutoffs, save to dCL and dCH.                            */
/* ----------------------------------------------------------------------- */ 
int SeqClustDP_QuantileCutoff(char strOutputPath[], char strOutputFile[], 
		int nChrNum, struct tagString **vChrName, struct INTMATRIX *pChrLen, 
		int nBinSize, double dBinMax, double dBinMin, 
		int nGridNum, double dCutL, double dCutH, double *pCL, double *pCH,
		struct BYTEMATRIX *pIncludeChr);

/* ----------------------------------------------------------------------- */ 
/*  SeqClustDP_Segmentation()                                              */
/*  Genome segmentation. Divide genome into correlation blocks.            */
/* ----------------------------------------------------------------------- */ 
int SeqClustDP_Segmentation(char strOutputPath[], char strOutputFile[], 
				int nChrNum, struct tagString **vChrName, struct INTMATRIX *pChrLen, 
				int nBinSize, double dCL, double dCH);

/* ----------------------------------------------------------------------- */ 
/*  SeqClustDP_FindBlock()                                                 */
/*  Find correlation blocks.                                               */
/* ----------------------------------------------------------------------- */ 
int SeqClustDP_FindBlock(char strOutputPath[], char strOutputFile[], 
				int nFileNum, struct tagString **vFileName, 
				struct DOUBLEMATRIX *pFileNormFactor, double dSmoothFactor,
				int nChrNum, struct tagString **vChrName, 
				struct INTMATRIX *pChrLen, int nBinSize, 
				struct BYTEMATRIX *pIncludeChr, int nBlockLenCut);

/* ----------------------------------------------------------------------- */ 
/*  SeqClustDP_FindBlock_SingleRegion()                                    */
/*  Using dynamic programming to find blocks in a single region.           */
/* ----------------------------------------------------------------------- */ 
int SeqClustDP_FindBlock_SingleRegion(char strChr[], int nStart, int nEnd, int nFileNum, 
				char strOutputPath[], char strOutputFile[], 
				struct tagString **vFileName, struct DOUBLEMATRIX *pFileNormFactor, 
				double dSmoothFactor, FILE *fpOut, int nBlockLenCut);

/* ----------------------------------------------------------------------- */ 
/*  SeqClustDP_FindBlock_SingleRegion_BinStandardize()                     */
/*  Using dynamic programming to find blocks in a single region.           */
/*  Standardize bin first.                                                 */
/* ----------------------------------------------------------------------- */ 
int SeqClustDP_FindBlock_SingleRegion_BinStandardize(char strChr[], int nStart, int nEnd, int nFileNum, 
				char strOutputPath[], char strOutputFile[], 
				struct tagString **vFileName, struct DOUBLEMATRIX *pFileNormFactor, 
				double dSmoothFactor, FILE *fpOut, int nBlockLenCut);

/* ----------------------------------------------------------------------- */ 
/*  SeqClustDP_BlockBIC_Normal()                                           */
/*  Compute BIC for a block using independent normal.                      */
/* ----------------------------------------------------------------------- */ 
double SeqClustDP_BlockBIC_Normal(float **vData, int nFileNum, int nStart, int nEnd, 
				int nRegionNum, double *vMu, double *vVar);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2WinBAR()                                                       */
/*  Convert alignment to bar file after read extension.                    */
/* ----------------------------------------------------------------------- */ 
int HTS_Aln2WinBAR(char strTXTFile[], char strOutputFolder[], char strBARHeader[], 
				   int nExtLen, int nBinSize, int nExportStrandProfile, 
				   char strSpecies[], char strChrLenFile[]);

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2WinBAR_Export()                                                */
/*  Convert alignment to bar file after read extension.                    */
/*  If nExportType == 0: only export + reads, reads extended both sides.   */
/*                 == 1: only export - reads, reads extended both sides.   */
/*                 == 2: export both + and - reads, reads extended         */
/*                       on single side toward DNA fragment center.        */
/* ----------------------------------------------------------------------- */ 
int HTS_Aln2WinBAR_Export(char strTXTFile[], char strOutputFolder[], 
				char strBARHeader[], int nExtLen, int nBinSize, 
				 char strSpecies[], int nChrNum, struct INTMATRIX *pChrLen,
				 int nExportType);

/* ----------------------------------------------------------------------- */ 
/*  SeqPeak_Main()                                                         */
/*  CisGenome Peak detection v2.                                           */
/* ----------------------------------------------------------------------- */ 
int SeqPeak_Main(char strInputPath[], char strOutputPath[], char strOutputFile[],
				  int nPaired, int nBinSize, int nExtLen,
				  int nSegType, int nWinSize,
				  int nCutType, int nTStandardize, double dCutoff, 
				  int nMaxGap, int nMinLen,
				  int nBoundaryRefine, int nBRWin,
				  int nExportBAR, int nKeepTempFiles, int nCollectRawData,
				  int nPoisFilter, int nPoisWin, double dPoisCut,
				  int nLFCAdj, int nLFCWin);

/* ----------------------------------------------------------------------- */ 
/*  SeqPeak_Initialize()                                                   */
/*  SeqPeak peak detection: Initialize                                     */
/* ----------------------------------------------------------------------- */ 
int SeqPeak_Initialize(char strInputPath[], int *pIPNum, int *pCTNum, int *pFileNum, 
						struct tagString ***vFileName, int **vGroupId,
						struct DOUBLEMATRIX **ppFileReadCount, 
						struct DOUBLEMATRIX **ppFileNormFactor, int *pBaseLineFileId, 
						int *pChrNum, 
						struct tagString ***vvChrName, struct INTMATRIX **ppChrLen);

/* ----------------------------------------------------------------------- */ 
/*  SeqPeak_CountBinReads()                                                */
/*  Count reads for genomic bins                                           */
/* ----------------------------------------------------------------------- */ 
int SeqPeak_CountBinReads(char strInputPath[], char strOutputPath[],
		char strOutputFile[], int nFileNum, int nChrNum, 
		struct tagString **vChrName, struct INTMATRIX *pChrLen, int nBinSize,
		int nExtLen, struct BYTEMATRIX *pIncludeChr,
		int nExportBAR);

/* ----------------------------------------------------------------------- */ 
/*  SeqPeak_CountBinReads_SingleFile()                                     */
/*  Count reads for genomic bins for a single sample file                  */
/* ----------------------------------------------------------------------- */ 
int SeqPeak_CountBinReads_SingleFile(char strInputFile[], char strOutputPath[],
				char strOutputFile[], int nChrNum, struct tagString **vChrName, 
				struct INTMATRIX *pChrLen, int nBinSize, 		
				int nExtLen, struct BYTEMATRIX *pIncludeChr,
				int nExportBAR);

/* ----------------------------------------------------------------------- */ 
/*  SeqPeak_ComputeBinStats()                                              */
/*  Compute summary statistics for genomic bins.                           */
/* ----------------------------------------------------------------------- */ 
int SeqPeak_ComputeBinStats(char strOutputPath[], char strOutputFile[], 
			int nFileNum, int nIPNum, int nCTNum, 
			struct tagString **vFileName, int *vGroupId, 
			struct DOUBLEMATRIX *pFileNormFactor, 
			int nChrNum, struct tagString **vChrName, struct INTMATRIX *pChrLen,
			int nBinSize, int nWinSize, double dSmoothFactor,
			struct DOUBLEMATRIX *pFileReadCount, int nTStandardize, 
			int nPoisFilter, int nPoisWin, 
			int nLFCAdj, int nLFCWin);

/* ----------------------------------------------------------------------- */ 
/*  SeqPeak_ComputeBinStats_Chr_Initial()                                  */
/*  Compute initial summary statistics for genomic bins for a single       */
/*  chromosome.                                                            */
/* ----------------------------------------------------------------------- */ 
int SeqPeak_ComputeBinStats_Chr_Initial(char strOutputPath[], char strOutputFile[], 
			int nFileNum, int nIPNum, int nCTNum,
			struct tagString **vFileName, int *vGroupId,
			struct DOUBLEMATRIX *pFileNormFactor, 
			char strChr[], int nBinNum, int nBinSize, int nWinSize, 
			double dSmoothFactor, struct DOUBLEMATRIX *pFileReadCount, 
			int nPoisFilter, int nPoisWin, 
			int nLFCAdj, int nLFCWin,
			double *pS2M, double *pS2S, int *pEffectBinNum);

/* ----------------------------------------------------------------------- */ 
/*  SeqPeak_ComputeBinStats_Chr_Tstat()                                    */
/*  Compute t-statistics.                                                  */
/* ----------------------------------------------------------------------- */ 
int SeqPeak_ComputeBinStats_Chr_Tstat(char strOutputPath[], char strOutputFile[], 
			char strChr[], int nBinNum, double dB, double dVM, int nIPNum, int nCTNum, 
			double *pTM, double *pTS);

/* ----------------------------------------------------------------------- */ 
/*  SeqPeak_ComputeBinStats_Chr_Tstandard()                                */
/*  Standardize t-statistics.                                              */
/* ----------------------------------------------------------------------- */
int SeqPeak_ComputeBinStats_Chr_Tstandard(char strOutputPath[], char strOutputFile[],
			char strChr[], int nBinNum, double dTM, double dTS);

/* ----------------------------------------------------------------------- */ 
/*  SeqPeak_CleanFiles()                                                   */
/*  Export results to BAR files and clean intermediate results.            */
/* ----------------------------------------------------------------------- */
int SeqPeak_CleanFiles(char strOutputPath[], char strOutputFile[],
		int nFileNum, struct tagString **vFileName,
		int nChrNum, struct tagString **vChrName, struct INTMATRIX * pChrLen, int nBinSize,
		int nExportBAR, int nKeepTempFiles, int nPoisFilter);

/* ----------------------------------------------------------------------- */ 
/*  SeqPeak_Segmentation()                                                 */
/*  Genome segmentation. Find peaks.                                       */
/* ----------------------------------------------------------------------- */ 
int SeqPeak_Segmentation(char strOutputPath[], char strOutputFile[], 
				int nChrNum, struct tagString **vChrName, struct INTMATRIX *pChrLen, 
				int nBinSize, double dCutoff, int nMaxGap, int nMinLen, 
				int nFlipGroupLabel, int nPoisFilter, 
				int nPoisFilterLabel, double dPoisCut);

/* ----------------------------------------------------------------------- */ 
/*  SeqPeak_FDR_Flip()                                                     */
/*  Compute FDR.                                                           */
/* ----------------------------------------------------------------------- */ 
int SeqPeak_FDR_Flip(char strOutputPath[], char strOutputFile[], 
			struct DOUBLEMATRIX **ppRegion,
			struct DOUBLEMATRIX **ppRegionSort, struct LONGMATRIX **ppRegionSid,
			struct DOUBLEMATRIX **ppRegionFDR);

/* ----------------------------------------------------------------------- */ 
/*  SeqPeak_ExportPeaks()                                                  */
/*  Export peaks.                                                          */
/* ----------------------------------------------------------------------- */ 
int SeqPeak_ExportPeaks(char strOutputPath[], char strOutputFile[], 
		struct DOUBLEMATRIX *pRegionSort, struct LONGMATRIX *pRegionSid, 
		struct DOUBLEMATRIX *pRegionFDR, int nBinSize, 
		int nChrNum, struct tagString **vChrName,
		int nBoundaryRefine, struct INTMATRIX *pPeakBoundary,
		int nFileNum, struct tagString **vFileName, 
		int nCollectRawData, struct DOUBLEMATRIX *pPeakRawData,
		int nPoisFilter, double dPoisCut);

/* ----------------------------------------------------------------------- */ 
/*  SeqPeak_RefineBoundary()                                               */
/*  Refine peak boundaries.                                                */
/* ----------------------------------------------------------------------- */ 
int SeqPeak_RefineBoundary(struct DOUBLEMATRIX *pRegion, int nBinSize, int nBRWin,
			char strInputPath[], char strOutputPath[], char strOutputFile[], 
			int nIPNum, int nCTNum, int nFileNum, struct tagString **vFileName, 
			int *vGroupId, struct DOUBLEMATRIX *pFileNormFactor, 
			int nChrNum, struct tagString **vChrName, 
			struct INTMATRIX *pChrLen, int nExtLen,
			struct INTMATRIX **ppPeakBoundary);

/* ----------------------------------------------------------------------- */ 
/*  SeqPeak_RefineBoundary_RegionProfile()                                 */
/*  Read region profile for boundary refinement.                           */
/* ----------------------------------------------------------------------- */ 
int SeqPeak_RefineBoundary_RegionProfile(int nRegionNum, int *vChrId, 
		int *vStart, int *vEnd, int *vBRWinNum, 
		float **vPosCount, float **vNegCount, float **vTotCount,
	    int nBRWin,	char strInputPath[], int nIPNum, int nCTNum, 
		int nFileNum, int *vGroupId, struct DOUBLEMATRIX *pFileNormFactor, 
		int nChrNum, struct tagString **vChrName, 
		struct INTMATRIX *pChrLen, int nExtLen);

/* ----------------------------------------------------------------------- */ 
/*  SeqPeak_RefineBoundary_FindOverlapRegions()                            */
/*  Find regions overlapping with an input region. Return region index.    */
/* ----------------------------------------------------------------------- */ 
int SeqPeak_RefineBoundary_FindOverlapRegions(int nChrId, int nP1, int nP2, 
		int nRegionNum, int *vChrId, int *vStart, int *vEnd, int *px, int *py);

/* ----------------------------------------------------------------------- */ 
/*  SeqPeak_RegionCompare()                                                */
/*  Compare two regions.                                                   */
/*  Return  -1 if region1 < region2;                                       */
/*           0 if region1 = region2;                                       */
/*           1 if region1 > region2;                                       */
/* ----------------------------------------------------------------------- */ 
int SeqPeak_RegionCompare(int nChr1, int nStart1, int nEnd1,
						  int nChr2, int nStart2, int nEnd2);

/* ----------------------------------------------------------------------- */ 
/*  SeqPeak_RefineBoundary_DetectBoundary()                                */
/*  Find boundaries.                                                       */
/* ----------------------------------------------------------------------- */ 
int SeqPeak_RefineBoundary_DetectBoundary(int nRegionNum, 
		int *vChrId, int *vStart, int *vEnd, int *vBRWinNum, 
		float **vPosCount, float **vNegCount, float **vTotCount,
		int nBRWin, int nExtLen, struct INTMATRIX *pPeakBoundary);

/* ----------------------------------------------------------------------- */ 
/*  SeqPeak_CollectRawData()                                               */
/*  Collect data for reporting.                                            */
/* ----------------------------------------------------------------------- */ 
int SeqPeak_CollectRawData(struct DOUBLEMATRIX *pRegion, 
		char strOutputPath[], char strOutputFile[], 
		int nFileNum, struct tagString **vFileName, 
		struct DOUBLEMATRIX *pFileNormFactor, double dSmoothFactor,
		int nChrNum, struct tagString **vChrName, 
		struct INTMATRIX *pChrLen, int nBinSize, 
		struct DOUBLEMATRIX **ppPeakRawData);


/* ----------------------------------------------------------------------- */ 
/*  SeqPeak_CollectRawData_SingleFile()                                    */
/*  Collect data for one file for reporting.                               */
/* ----------------------------------------------------------------------- */ 
int SeqPeak_CollectRawData_SingleFile(int nRegionNum, 
		int *vChrId, int *vStart, int *vEnd,
		char strOutputPath[], char strOutputFile[], 
		char strFileName[], double dFileNormFactor,
		double dSmoothFactor, int nChrNum, 
		struct tagString **vChrName, struct INTMATRIX *pChrLen, 
		int nBinSize, struct DOUBLEMATRIX *pPeakRawData, 
		int nCol);











