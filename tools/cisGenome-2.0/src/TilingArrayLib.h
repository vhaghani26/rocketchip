/* ----------------------------------------------------------------------- */
/*                                 Macro                                   */
/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------- */ 
/*                                 Struct                                  */
/* ----------------------------------------------------------------------- */ 

/* probe genome mapping */
struct tagProbeGenomeInfo
{
	/* array coordinates */
	int nX;
	int nY;

	/* probe sequence */
	struct tagString *pProbe;

	/* genomic coordinates */
	int nChr;
	int nPos;

	/* next */
	struct tagProbeGenomeInfo *pNext;
};


/* genome pos */
struct tagChrPos
{
	/* position */
	int nPos;

	/* next */
	struct tagChrPos *pNext;
};

/* genome pos */
struct tagChrPosList
{
	/* number of nodes */
	int nNodeNum;

	/* next */
	struct tagChrPos *pPosList;
};

/* tilemap parameters */
struct tagTileMapv2Param
{
	/* comparison type */
	int nComparisonType;

	/* working directory */
	char strWorkPath[MED_LINE_LENGTH];

	/* project title */
	char strProjectTitle[LINE_LENGTH];

	/* No. of Libraries, Samples, Arrays and Groups */
	int nLibNum;
	int nSampleNum;
	int nGroupNum;
	int nArrayNum;

	/* Sample name */
	struct INTMATRIX *pGroupLabel;
	struct INTMATRIX *pGroupSize;
	struct tagString **vSampleAlias;
	struct tagString **vArrayFile;

	/* Pattern of Interest */
	int nPatternType;
	char strPattern[MED_LINE_LENGTH];

	/* Preprocessing steps */
	int nNoiseMask;
	double dLower;
	double dUpper;
	int nTransformType;

	/* Probe summary */
	int nMonteCarloNum;
	int nVargroupNum;
	struct INTMATRIX *pVargroupSize;
	struct INTMATRIX **vVargroupMap;

	/* Region Summary */
	int nRegionSummaryType;

	/* MA parameters */
	int nW;
	int nWindowBoundary;
	int nMAStandardize;
	double dBaseStd;
	int nFDRType;

	/* HMM parameters */
	int nExpLen;
	double dPostCut;
	int nHMMParamUserSpecified;
	char strTransitionPath[MED_LINE_LENGTH];
	char strEmissionPath[MED_LINE_LENGTH];
	
	/* UMS parameters */
	double dTp;
	double dTq;
	int nOffset;
	int nGridSize;
	
	/* permutation parameters */
	int nPermutationNum;
	int nPermgroupNum;
	struct INTMATRIX *pPermgroupSize;
	struct INTMATRIX **vPermgroupMap;

	/* postfiltering parameters */
	int nGap;
	int nGapW;
	int nMinRegLen;
	int nMinRegProbeNum;
};

struct tagTileMapv2InfoTrack
{
	/* track name */
	char strTrackName[MED_LINE_LENGTH];

	/* data file */
	char strDataFile[MED_LINE_LENGTH];

	/* filter number */
	int nFilterNum;

	/* filter file */
	struct tagString **vFilterFile;

	/* filter type */
	struct BYTEMATRIX *pFilterType;

	/* filter value */
	struct DOUBLEMATRIX *pFilterValue;

	/* next */
	struct tagTileMapv2InfoTrack *pNext;
};

/* ----------------------------------------------------------------------- */
/*                              Declaration                                */
/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ProbeReMapping_Main()                                        */
/*  Remap probes to a new genome version.                                  */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_ProbeReMapping_Main(char strBpmapFile[], char strGenomePath[],
						char strOutputFile[], char strSpecies[],
						int nKeyLen);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ProbeReMapping_ProbeAlign()                                  */
/*  Align probes to genome.                                                */
/* ----------------------------------------------------------------------- */ 
struct tagAffyBpMapUnit *TileMapv2_ProbeReMapping_ProbeAlign(char strBpmapFile[], 
				char strSpecies[], char strSeqFile[], int nChrLen, 
				int nKeyLen, int nBaseTypeNum,
				int *vGenomeHashIndex, int nIndexNum,
				int *vGenomeHashTable, int nTableNum,
				int *pProbeNum);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ProbeReMapping_HashAlign()                                   */
/*  Align a probe to a chromosome.                                         */
/* ----------------------------------------------------------------------- */ 
struct tagAffyBpMapUnit *TileMapv2_ProbeReMapping_HashAlign(struct tagAffyBpMapUnit *pUnit,
			char strSeqFile[], int nChrLen, int nKeyLen, int nBaseTypeNum,
			int *vGenomeHashIndex, int nIndexNum,
			int *vGenomeHashTable, int nTableNum,
			int nForward);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ProbeReMapping_SortProbeList()                               */
/*  Resorting probes according to their position.                          */
/* ----------------------------------------------------------------------- */ 
struct tagAffyBpMapUnit **TileMapv2_ProbeReMapping_SortProbeList(int nProbeNum, 
		struct tagAffyBpMapUnit **ppProbeList);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ProbeReMapping_WriteBpmapFile()                              */
/*  Write remapped probe to a new bpmap file.                              */
/* ----------------------------------------------------------------------- */
int TileMapv2_ProbeReMapping_WriteBpmapFile(char strOutputFile[], 
			char strSpecies[], int nChrNum, struct tagString **vChrName, 
			struct tagAffyBpMapUnit ***vvProbeMap, struct INTMATRIX *pProbeNum);

/* ----------------------------------------------------------------------- */ 
/*  ProbeGenomeInfoCreate()                                                */
/*  Create Probe Genome Info object.                                       */
/* ----------------------------------------------------------------------- */ 
struct tagProbeGenomeInfo *ProbeGenomeInfoCreate();

/* ----------------------------------------------------------------------- */ 
/*  ProbeGenomeInfoDestroy()                                               */
/*  Destroy Probe Genome Info object.                                      */
/* ----------------------------------------------------------------------- */ 
void ProbeGenomeInfoDestroy(struct tagProbeGenomeInfo **pInfo);

/* ----------------------------------------------------------------------- */ 
/*  ChrPosCreate()                                                         */
/*  Create ChrPos object.                                                  */
/* ----------------------------------------------------------------------- */ 
struct tagChrPos *ChrPosCreate();

/* ----------------------------------------------------------------------- */ 
/*  ChrPosDestroy()                                                        */
/*  Destroy ChrPos object.                                                 */
/* ----------------------------------------------------------------------- */ 
void ChrPosDestroy(struct tagChrPos **pChrPos);

/* ----------------------------------------------------------------------- */ 
/*  ChrPosListCreate()                                                     */
/*  Create ChrPos List object.                                             */
/* ----------------------------------------------------------------------- */ 
struct tagChrPosList *ChrPosListCreate();

/* ----------------------------------------------------------------------- */ 
/*  ChrPosListDestroy()                                                    */
/*  Destroy ChrPosList object.                                             */
/* ----------------------------------------------------------------------- */ 
void ChrPosListDestroy(struct tagChrPosList **pChrPosList);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ImportAffy_Main()                                            */
/*  Create a tiling array analysis project by loading affymetrix CEL data  */
/*  The data will be normalized, mapped to chromosomes, and saved to *.bar */
/*  files. If specified by users, intensities will also be exported to a   */
/*  combined .txt file.                                                    */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_ImportAffy_Main(char strParamPath[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_LoadCEL()                                                    */
/*  Load CEL files.                                                        */
/* ----------------------------------------------------------------------- */ 
struct tagCELData *TileMapv2_LoadCEL(char strFileName[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_GetMask()                                                    */
/*  Get mask indicators.                                                   */
/* ----------------------------------------------------------------------- */ 
struct BYTEMATRIX *TileMapv2_GetMask(struct tagCELData *pCELData, 
	int nIncludeMasked, int nIncludeOutlier,  int nIncludeModified);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_SaveToBinaryFile()                                           */
/*  Write data to binary files.                                            */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_SaveToBinaryFile(const void *buffer, size_t size, size_t count, char strFileName[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_LoadFromBinaryFile()                                           */
/*  Load data from binary files.                                            */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_LoadFromBinaryFile(void *buffer, size_t size, size_t count, char strFileName[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_QuantileNorm_AddQuantile()                                   */
/*  Add quantiles for quantile normalization.                              */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *TileMapv2_QuantileNorm_AddQuantile(struct DOUBLEMATRIX *pSortMean, 
			struct DOUBLEMATRIX *pArray, struct LONGMATRIX *pSortId, 
			struct BYTEMATRIX *pMask, int nMaskedCells);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_QuantileNorm_Rescale()                                       */
/*  Rescale arrays according to the quantile normalized values.            */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_QuantileNorm_Rescale(char strPrcPath[], struct DOUBLEMATRIX *pSortMean, 
								   char strOutPath[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_BpmapToBAR()                                                 */
/*  Create BAR template from bpmap file. The function will create local    */
/*  repeat mask at the same time.                                          */
/* ----------------------------------------------------------------------- */
struct tagBARData *TileMapv2_BpmapToBAR(char strBpmapFile[], char strMaskPath[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ExportAffyIntensity()                                        */
/*  Export intensity data to files.                                        */
/*  nIntensityType = 0: PMonly; 1:PM-MM.                                   */
/*  nApplyMask = 0: no mask; 1: using masks.                               */
/*  nMode = 0: separate bar files; 1: a combined bar file; 2: a text file  */
/* ----------------------------------------------------------------------- */
int TileMapv2_ExportAffyIntensity(char strWorkPath[], char strProjectTitle[],
			struct tagString **vCELPath,
			int nCELNumberCells, int nCELTotalX, int nCELTotalY,
			int nSampleNum, struct tagString **vSampleAlias,  
			int nLibNum, int nLibId, struct tagBARData *pBARPos, 
			int nIntensityType, double dIntLowerBound, int nIntLogTransform,
			int nApplyMask, int nExportMode);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ExportAffyIntensity_MultiBAR()                               */
/*  Export intensity data to separate bar files.                           */
/* ----------------------------------------------------------------------- */
int TileMapv2_ExportAffyIntensity_MultiBAR(struct tagString **vExportPath, 
			int nSampleNum, struct tagString **vSampleAlias, FILE **vfpIn, 
			int nCELNumberCells, int nCELTotalX, int nCELTotalY,
			int nApplyMask, FILE **vfpMask,	struct tagBARData *pBARPos, 
			int nIntensityType,	double dIntLowerBound, int nIntLogTransform);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ExportAffyIntensity_SingleBAR()                              */
/*  Export intensity data to a single bar files.                           */
/* ----------------------------------------------------------------------- */
int TileMapv2_ExportAffyIntensity_SingleBAR(char strExportPath[], 
			int nSampleNum, struct tagString **vSampleAlias, FILE **vfpIn,
			int nCELNumberCells, int nCELTotalX, int nCELTotalY,
			int nApplyMask, FILE **vfpMask,	struct tagBARData *pBARPos, 
			int nIntensityType,	double dIntLowerBound, int nIntLogTransform);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ExportAffyIntensity_SingleTXT()                              */
/*  Export intensity data to a single txt files.                           */
/* ----------------------------------------------------------------------- */
int TileMapv2_ExportAffyIntensity_SingleTXT(char strExportPath[], 
			int nSampleNum, struct tagString **vSampleAlias, FILE **vfpIn,
			int nCELNumberCells, int nCELTotalX, int nCELTotalY,
			int nApplyMask, FILE **vfpMask,	struct tagBARData *pBARPos, 
			int nIntensityType,	double dIntLowerBound, int nIntLogTransform);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2Param_Create()                                                */
/*  Create TileMapv2Param object.                                          */
/* ----------------------------------------------------------------------- */ 
struct tagTileMapv2Param *TileMapv2Param_Create();

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2Param_Destroy()                                               */
/*  Create TileMapv2Param object.                                          */
/* ----------------------------------------------------------------------- */ 
int TileMapv2Param_Destroy(struct tagTileMapv2Param **pParam);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_Main()                                                       */
/*  TileMapv2 pipeline                                                     */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_Main(char strParamPath[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_LoadParameters()                                             */
/*  Load tilemap parameters.                                               */
/* ----------------------------------------------------------------------- */ 
struct tagTileMapv2Param *TileMapv2_LoadParameters(char strParamFile[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2()                                                            */
/*  TileMap for two sample comparisons.                                    */
/* ----------------------------------------------------------------------- */ 
int TileMapv2(struct tagTileMapv2Param *pParam);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ProbeSelection_Main()                                        */
/*  TileMapv2 probe level summary.                                         */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_ProbeSelection_Main(struct tagTileMapv2Param *pParam, 
						int nLibId);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ProbeSelection_DataTransformation()                          */
/*  TileMapv2 transform raw data.                                          */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_ProbeSelection_DataTransformation(struct tagBARData *pData,
						struct tagTileMapv2Param *pParam);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ProbeSelection_OneSample_WithMask()                          */
/*  TileMapv2 probe level summary: two sample comparisons.                 */
/*  this function can handle masked cells.                                 */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_ProbeSelection_OneSample_WithMask(struct tagBARData *pData, 
						struct tagBARData *pMask, struct tagTileMapv2Param *pParam, 
						char strOutFile[], char strFCFile[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ProbeSelection_TwoSample()                                   */
/*  TileMapv2 probe level summary: two sample comparisons.                 */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_ProbeSelection_TwoSample(struct tagBARData *pData,
						struct tagTileMapv2Param *pParam, char strOutFile[],
						char strFCFile[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ProbeSelection_TwoSample_WithMask()                          */
/*  TileMapv2 probe level summary: two sample comparisons.                 */
/*  this function can handle masked cells.                                 */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_ProbeSelection_TwoSample_WithMask(struct tagBARData *pData, 
						struct tagBARData *pMask, struct tagTileMapv2Param *pParam, 
						char strOutFile[], char strFCFile[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ProbeSelection_MultiSample()                                 */
/*  TileMapv2 probe level summary: multiple sample comparisons.            */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_ProbeSelection_MultiSample(struct tagBARData *pData,
						struct tagTileMapv2Param *pParam, char strOutFile[],
						char strFCFile[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ProbeSelection_MultiSample_Fast()                            */
/*  TileMapv2 probe level summary: multiple sample comparisons.            */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_ProbeSelection_MultiSample_Fast(struct tagBARData *pData,
						struct tagTileMapv2Param *pParam, char strOutFile[],
						char strFCFile[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ProbeSelection_MultiSample_WithMask_Fast()                   */
/*  TileMapv2 probe level summary: multiple sample comparisons.            */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_ProbeSelection_MultiSample_WithMask_Fast(struct tagBARData *pData,
						struct tagBARData *pMask, struct tagTileMapv2Param *pParam, 
						char strOutFile[], char strFCFile[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ProbeSelection_RemoveMean()                                  */
/*  Remove mean from each group, this is a preparation for permutations.   */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_ProbeSelection_RemoveMean(struct tagBARData *pData, 
				struct tagTileMapv2Param *pParam, int nNoiseMask, 
				struct tagBARData *pMask);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_GetShrinkingVar_UnequalDF()                                  */
/*  Variance shrinking: unequal d.f. case.                                 */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_GetShrinkingVar_UnequalDF(struct DOUBLEMATRIX *pV, 
					struct DOUBLEMATRIX *pDf, int nMaxDf, int nProbeNum);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_GetShrinkingVar_FromSSR()                                    */
/*  Variance shrinking: unequal d.f. case.                                 */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_GetShrinkingVar_FromSSR(struct tagBARData *pData, int nMaxDf,
						int nVarId, int nDfId);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_RegionDetection_MA_Main()                                    */
/*  TileMapv2 region detection, MA.                                        */
/*  return number of regions detected.                                     */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_RegionDetection_MA_Main(struct tagTileMapv2Param *pParam, int nLibId);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_RegionDetection_MA_PrepareData()                             */
/*  TileMapv2 MA preparation: loading probe level statistics.              */
/* ----------------------------------------------------------------------- */ 
struct tagBARData *TileMapv2_RegionDetection_MA_PrepareData(char strProbeFile[],
							char strMaskFile[], char strFCFile[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_MA()                                                         */
/*  TileMapv2 MA.                                                          */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_MA(struct tagBARData *pData, struct tagTileMapv2Param *pParam);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_MA_CallRegion()                                              */
/*  TileMapv2 call binding regions based on MA statistics.                 */
/* ----------------------------------------------------------------------- */
struct DOUBLEMATRIX *TileMapv2_MA_CallRegion(struct tagBARData *pData, 
					struct tagTileMapv2Param *pParam, int nLibId);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_RegionDetection_MA_ResetData()                               */
/*  TileMapv2 MA reset data.                                               */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_RegionDetection_MA_ResetData(struct tagBARData *pData, 
					struct INTMATRIX *pCol);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_RegionFDR_Count()                                            */
/*  Count random observations for computing FDR based on column nCol.      */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_RegionFDR_Count(struct DOUBLEMATRIX *pRegionSort, 
							  struct DOUBLEMATRIX *pRegionCTSort, 
							  struct DOUBLEMATRIX *pRegionFDR, int nCol);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_RegionFDR_Compute()                                          */
/*  Compute FDR from the count.                                            */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_RegionFDR_Compute(struct DOUBLEMATRIX *pFDR);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_RegionFDR_Reorganize()                                       */
/*  Reorganize regions by adding FDR.                                      */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *TileMapv2_RegionFDR_Reorganize(struct DOUBLEMATRIX *pRegion, 
				struct DOUBLEMATRIX *pRegionMaxFDR, struct LONGMATRIX *pRegionMaxSid, 
				struct DOUBLEMATRIX *pRegionSumFDR, struct LONGMATRIX *pRegionSumSid);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_RegionFDR_LocalFDR()                                         */
/*  Compute local fdr.                                                     */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX* TileMapv2_RegionFDR_LocalFDR(struct DOUBLEMATRIX* pFDR);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_MA_ExportRegion()                                            */
/*  Export regions.                                                        */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_MA_ExportRegion(struct DOUBLEMATRIX *pRegion, struct tagBARData *pData, 
							  struct tagTileMapv2Param *pParam, int nLibId);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_MA_RegionFDR_Permutation()                                   */
/*  TileMap Region FDR.                                                    */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_MA_RegionFDR_Permutation(struct INTMATRIX *pType,
				struct DOUBLEMATRIX *pRegionMaxSort, struct DOUBLEMATRIX *pRegionMaxFDR,
				struct DOUBLEMATRIX *pRegionSumSort, struct DOUBLEMATRIX *pRegionSumFDR, 
				struct tagTileMapv2Param *pParam, int nLibId);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_MA_RegionFDR_LeftTail()                                      */
/*  TileMap Region FDR.                                                    */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_MA_RegionFDR_LeftTail(struct INTMATRIX *pType,
				struct DOUBLEMATRIX *pRegionMaxSort, struct DOUBLEMATRIX *pRegionMaxFDR,
				struct DOUBLEMATRIX *pRegionSumSort, struct DOUBLEMATRIX *pRegionSumFDR, 
				struct tagBARData *pData, struct tagTileMapv2Param *pParam, int nLibId);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_MA_RegionFDR_UMS()                                           */
/*  TileMap Region FDR: UMS                                                */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *TileMapv2_MA_RegionFDR_UMS(struct DOUBLEMATRIX *pRegion,
				struct tagBARData *pData, struct tagTileMapv2Param *pParam, 
				int nLibId);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_SummarizeResults_Main()                                      */
/*  Summarize region detection results.                                    */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_SummarizeResults_Main(struct tagTileMapv2Param *pParam, int nRegionNum);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_RegionDetection_HMM_Main()                                   */
/*  TileMapv2 region detection, HMM.                                       */
/*  return number of regions detected.                                     */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_RegionDetection_HMM_Main(struct tagTileMapv2Param *pParam, int nLibId);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_RegionDetection_HMM_PrepareData()                            */
/*  TileMapv2 HMM preparation: loading probe level statistics.             */
/* ----------------------------------------------------------------------- */ 
struct tagBARData *TileMapv2_RegionDetection_HMM_PrepareData(char strProbeFile[],
							char strMaskFile[], char strFCFile[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_HMM()                                                        */
/*  TileMapv2 HMM.                                                         */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_HMM(struct tagBARData *pData, struct tagTileMapv2Param *pParam, 
				  int nLibId);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_UMS_HMM_Main()                                               */
/*  TileMapv2 UMS HMM.                                                     */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_UMS_HMM_Main(struct tagBARData *pData, struct tagTileMapv2Param *pParam,
				 int nLibId);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_HMM_Decoding()                                               */
/*  TileMapv2 HMM decoding.                                                */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_HMM_Decoding(struct tagBARData *pData, struct tagTileMapv2Param *pParam, 
						   int nLibId);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_HMM_Decoding_Chr()                                           */
/*  TileMapv2 HMM decoding.                                                */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_HMM_Decoding_Chr(int nProbeNum, int nStateNum, 
			struct DOUBLEMATRIX *pStationary, struct DOUBLEMATRIX *pTransition, 
			struct DOUBLEMATRIX *pEmission, double dGapDist, 
			struct DOUBLEMATRIX *pPosition, struct DOUBLEMATRIX *pScore,
			struct DOUBLEMATRIX *pPosterior, struct DOUBLEMATRIX *pMask);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_HMM_CallRegion()                                             */
/*  TileMapv2 call binding regions based on HMM statistics.                */
/* ----------------------------------------------------------------------- */
struct DOUBLEMATRIX *TileMapv2_HMM_CallRegion(struct tagBARData *pData, 
					struct tagTileMapv2Param *pParam, int nLibId);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_HMM_RegionFDR_Permutation()                                  */
/*  TileMap Region FDR.                                                    */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_HMM_RegionFDR_Permutation(struct INTMATRIX *pType,
				struct DOUBLEMATRIX *pRegionMaxSort, struct DOUBLEMATRIX *pRegionMaxFDR,
				struct DOUBLEMATRIX *pRegionSumSort, struct DOUBLEMATRIX *pRegionSumFDR, 
				struct tagTileMapv2Param *pParam, int nLibId);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_HMM_RegionFDR_LeftTail()                                     */
/*  TileMap Region FDR.                                                    */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_HMM_RegionFDR_LeftTail(struct INTMATRIX *pType,
				struct DOUBLEMATRIX *pRegionMaxSort, struct DOUBLEMATRIX *pRegionMaxFDR,
				struct DOUBLEMATRIX *pRegionSumSort, struct DOUBLEMATRIX *pRegionSumFDR, 
				struct tagBARData *pData, struct tagTileMapv2Param *pParam, int nLibId);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_DataTransform()                                              */
/*  Transform column nCol according to strTransform.                       */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_DataTransform(struct tagBARData *pData, int nCol, 
							char strTransform[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_RegionInfo_Main()                                            */
/*  TileMapv2 get region information                                       */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_RegionInfo_Main(char strRegionPath[], char strInfoPath[], char strOutPath[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_InfoTrack_Create()                                           */
/*  TileMapv2 create information track.                                    */
/* ----------------------------------------------------------------------- */ 
struct tagTileMapv2InfoTrack *TileMapv2_InfoTrack_Create(int nFilterNum);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_InfoTrack_Destroy()                                           */
/*  TileMapv2 destroy information track.                                    */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_InfoTrack_Destroy(struct tagTileMapv2InfoTrack **pTrack);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_RegionInfo_LoadTrackInfo()                                   */
/*  TileMapv2 load info track information                                  */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_RegionInfo_LoadTrackInfo(char strInfoPath[], 
				struct tagTileMapv2InfoTrack **pTrackList);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_RegionInfo_CollectTrackInfo()                                */
/*  TileMapv2 collect track information                                    */
/* ----------------------------------------------------------------------- */
int TileMapv2_RegionInfo_CollectTrackInfo(char strRegionPath[], 
			struct tagTileMapv2InfoTrack *pTrack, 
			struct DOUBLEMATRIX *pInfo, int nCol);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_RegionInfo_Integral_Main()                                   */
/*  TileMapv2 get region information                                       */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_RegionInfo_Integral_Main(char strRegionPath[], char strInfoPath[], char strOutPath[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_RegionInfo_Integral_CollectTrackInfo()                       */
/*  TileMapv2 collect track information                                    */
/* ----------------------------------------------------------------------- */
int TileMapv2_RegionInfo_Integral_CollectTrackInfo(char strRegionPath[], 
			struct tagTileMapv2InfoTrack *pTrack, 
			struct DOUBLEMATRIX *pInfo, int nCol);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_CollectProbes_Main()                                         */
/*  TileMapv2 get probe information                                        */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_CollectProbes_Main(char strRegionPath[], char strInfoPath[], char strOutPath[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_CollectProbes_ProcessRegion()                                */
/*  TileMapv2 get probe information                                        */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_CollectProbes_ProcessRegion(char strChr[], int nStart, int nEnd, 
							struct tagBARData **vData, int nTrackNum, FILE *fpOut);


/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_TXT2BAR()                                                    */
/*  Convert a text file to BAR tiling array project.                       */
/* ----------------------------------------------------------------------- */
int TileMapv2_TXT2BAR(char strTXTFile[], char strExportFolder[], char strProjectTitle[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_TXT_QuantileNormalization()                                  */
/*  Quantile normalization for a tab-delimited txt file.                   */
/*  Return PROC_SUCCESS if success.                                        */
/*  nTransform: 0=identity (default); 1=log2.                              */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_TXT_QuantileNormalization(char strDataFile[], char strOutFile[],
			int nSkipColNum, int nTransform, double dTruncLow, double dTruncHigh);

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_Build_Main()                                                 */
/*  Build a probe background model for microarray CEL files.               */
/* ----------------------------------------------------------------------- */ 
int TileProbe_Build_Main(char strInputPath[], char strOutputPath[], 
						 int nTakeLog, double dNormLowerBound, int nShrink, int nTest);

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_Norm_Main()                                                  */
/*  Adjust probe intensities based on the probe background model to remove */
/*  probe effects.                                                         */
/*  nInputType: 0=single input cel file; 1=a file of array list            */
/*  strOutputPath: folder of output                                        */
/* ----------------------------------------------------------------------- */ 
int TileProbe_Norm_Main(char strInputPath[], int nInputType,
						char strOutputPath[], char strModelPath[],
						int nTakeLog, double dNormLowerBound,
						double dB, int nLogAfterNorm);

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_Peak_Main()                                                  */
/*  Find peaks.                                                            */
/* ----------------------------------------------------------------------- */ 
int TileProbe_Peak_Main(char strInputPath[]);

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_Peak()                                                       */
/*  Find peaks.                                                            */
/* ----------------------------------------------------------------------- */ 
int TileProbe_Peak(int nIPNum, int nCTNum, struct tagBARData *vIP, struct tagBARData *vCT,
				   int nBandWidth, int nMaxGap, int nMinProbe, int nUseVar, double dCut,
				   char strOutputPath[], char strGrpName[]);

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_Peak_Call()                                                  */
/*  Call peaks.                                                            */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *TileProbe_Peak_Call(int nTCol, struct tagBARData *vIP, 
	double dCut, int nMaxGap, int nPosVal, char strOutputPath[]);

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_MAT_Main()                                                   */
/*  MAT background correction for tileprobe.                               */
/* ----------------------------------------------------------------------- */ 
int TileProbe_MAT_Main(char strParamPath[]);


/* ----------------------------------------------------------------------- */ 
/*  TileProbe_BpmapToBAR()                                                 */
/*  Create BAR template from bpmap file. The function will create local    */
/*  repeat mask at the same time.                                          */
/* ----------------------------------------------------------------------- */
struct tagBARData *TileProbe_BpmapToBAR(char strBpmapFile[], char strMaskPath[]);

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_MAT_XX()                                                     */
/*  Obtain X'X.                                                            */
/* ----------------------------------------------------------------------- */
struct DOUBLEMATRIX *TileProbe_MAT_XX(struct tagBARData *pBARPos);

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_MAT_XY()                                                     */
/*  Obtain X'Y.                                                            */
/* ----------------------------------------------------------------------- */
struct DOUBLEMATRIX *TileProbe_MAT_XY(struct tagBARData *pBARPos, struct tagCELData *pCELData);

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_MAT_Correction()                                             */
/*  Background correction.                                                 */
/* ----------------------------------------------------------------------- */
int TileProbe_MAT_Correction(struct tagBARData *pBARPos, struct DOUBLEMATRIX *pBeta,
							 char strOutFile[]);

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_BARBuild_Main()                                              */
/*  Build a probe background model for microarray BAR files.               */
/* ----------------------------------------------------------------------- */ 
int TileProbe_BARBuild_Main(char strInputPath[], char strOutputPath[], 
						 int nTakeLog, double dNormLowerBound, int nShrink, int nTest);

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_BARNorm_Main()                                               */
/*  Adjust probe intensities based on the probe background model to remove */
/*  probe effects.                                                         */
/*  nInputType: 0=single input bar file; 1=a file of array list            */
/*  strOutputPath: folder of output                                        */
/*  dB shrinkage factor for Q50-Q25.                                       */
/* ----------------------------------------------------------------------- */ 
int TileProbe_BARNorm_Main(char strInputPath[], int nInputType,
						char strOutputPath[], char strModelPath[],
						int nTakeLog, double dNormLowerBound,
						double dB);

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_MATv2_Main()                                                 */
/*  MAT background correction for tileprobe.                               */
/* ----------------------------------------------------------------------- */ 
int TileProbe_MATv2_Main(char strParamPath[]);

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_BpmapToBARv2()                                               */
/*  Create BAR template from bpmap file. The function will create local    */
/*  repeat mask at the same time.                                          */
/* ----------------------------------------------------------------------- */
struct tagBARData *TileProbe_BpmapToBARv2(char strBpmapFile[], char strMaskPath[],
				char strGenomeGrp[], int nIncludeNonGrp);

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_MATv2_XX()                                                   */
/*  Obtain X'X.                                                            */
/* ----------------------------------------------------------------------- */
struct DOUBLEMATRIX *TileProbe_MATv2_XX(struct tagBARData *pBARPos, char strGenomeGrp[],
				int nIncludeNonGrp);

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_MATv2_XY()                                                   */
/*  Obtain X'Y.                                                            */
/* ----------------------------------------------------------------------- */
struct DOUBLEMATRIX *TileProbe_MATv2_XY(struct tagBARData *pBARPos, struct tagCELData *pCELData,
				char strGenomeGrp[], int nIncludeNonGrp);

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_MATv2_Correction()                                           */
/*  Background correction.                                                 */
/* ----------------------------------------------------------------------- */
int TileProbe_MATv2_Correction(struct tagBARData *pBARPos, struct DOUBLEMATRIX *pBeta,
				char strOutFile[], char strGenomeGrp[], int nIncludeNonGrp);

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_BARBuildv2_Main()                                            */
/*  Build a probe background model for microarray BAR files.               */
/*  The variance of each probe is the mean within group variance           */
/* ----------------------------------------------------------------------- */ 
int TileProbe_BARBuildv2_Main(char strInputPath[], char strOutputPath[], 
						 int nTakeLog, double dNormLowerBound, int nShrink, int nTest);

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_Buildv2_Main()                                               */
/*  Build a probe background model for microarray CEL files.               */
/* ----------------------------------------------------------------------- */ 
int TileProbe_Buildv2_Main(char strInputPath[], char strOutputPath[], 
						 int nTakeLog, double dNormLowerBound, int nShrink, 
						 int nLogAfterNorm, int nTest);

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_BuildHMMT_Main()                                             */
/*  Build a probe background model for microarray CEL files using HMM      */
/*  tiling algorithm.                                                      */
/* ----------------------------------------------------------------------- */ 
int TileProbe_BuildHMMT_Main(char strInputPath[], char strOutputPath[], 
						 int nTakeLog, double dNormLowerBound, 
						 double dFilterTopPrc, int nTest);

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_BuildHMMB_Main()                                             */
/*  Build a probe background model for microarray CEL files using HMM      */
/*  tiling algorithm with variance shrinking.                              */
/* ----------------------------------------------------------------------- */ 
int TileProbe_BuildHMMB_Main(char strInputPath[], char strOutputPath[], 
						 int nTakeLog, double dNormLowerBound, 
						 double dFilterTopPrc, int nTest);

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_BuildHMMM_Main()                                             */
/*  Build a probe background model for microarray CEL files using HMM      */
/*  tiling algorithm with variance shrinking.                              */
/* ----------------------------------------------------------------------- */ 
int TileProbe_BuildHMMM_Main(char strInputPath[], char strOutputPath[], 
						 int nTakeLog, double dNormLowerBound, 
						 double dFilterTopPrc, int nTest);

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_BuildHMMW_Main()                                             */
/*  Build a probe background model for microarray CEL files using HMM      */
/*  tiling algorithm with variance shrinking.                              */
/* ----------------------------------------------------------------------- */ 
int TileProbe_BuildHMMW_Main(char strInputPath[], char strOutputPath[], 
						 int nTakeLog, double dNormLowerBound, 
						 double dFilterTopPrc, int nShrink, int nTest);

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_HuberData_Main()                                             */
/*  Get data for testing Huber et al.                                      */
/* ----------------------------------------------------------------------- */ 
int TileProbe_HuberData_Main(char strParamPath[]);

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_BpmapToBARv2()                                               */
/*  Create BAR template from bpmap file. The function will create local    */
/*  repeat mask at the same time.                                          */
/* ----------------------------------------------------------------------- */
struct tagBARData *TileProbe_BpmapToBARHuber(char strBpmapFile[], char strMaskPath[],
			char strGenomeGrp[], int nIncludeNonGrp);

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_Resample_Main()                                              */
/*  Test robustness by resampling                                          */
/* ----------------------------------------------------------------------- */
int TileProbe_Resample_Main(char strParamFile[]);

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_Resample_Model()                                             */
/*  Build probe model                                                      */
/* ----------------------------------------------------------------------- */
int TileProbe_Resample_Model(int nArrayNum, int nProbeNum, int nShrink,
							 int *vGroupID, struct tagBARData **vBARData, 
							 struct DOUBLEMATRIX *pM, struct DOUBLEMATRIX *pV);

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_Resample_MATPeak()                                           */
/*  Find peaks.                                                            */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *TileProbe_Resample_MATPeak(int nIPNum, int nCTNum, int nProbeNum,
				   struct tagBARData *vIP, struct tagBARData *vCT,
				   int nBandWidth, int nMaxGap, int nMinProbe, int nUseVar,
				   struct DOUBLEMATRIX *pM, struct DOUBLEMATRIX *pV, int nDivideV);

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_Resample_Residual()                                          */
/*  Build probe model                                                      */
/* ----------------------------------------------------------------------- */
int TileProbe_Resample_Residual(int nArrayNum, int nProbeNum,
							 int *vGroupID, struct tagBARData **vBARData);

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_Resample_VarModel()                                          */
/*  Build probe model                                                      */
/* ----------------------------------------------------------------------- */
int TileProbe_Resample_VarModel(int nArrayNum, int nProbeNum, int nShrink,
							 int *vGroupID, struct tagBARData **vBARData, 
							 struct DOUBLEMATRIX *pV);