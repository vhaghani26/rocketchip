/* ----------------------------------------------------------------------- */
/*  MicroarrayLib.h : interface of the microarray library                  */
/*  Author : Ji HongKai ; Time: 2004.08                                    */
/* ----------------------------------------------------------------------- */

#define POWEXPRESSLOC_MINWIN 10

#define TRANSLOC_PREGENENUM 5
#define TRANSLOC_REFGENEMATCH_TH 0.001
#define TRANSLOC_REPEATMASK_WIN 5000000
#define TRANSLOC_GENECORR_WIN 10
#define TRANSLOC_DIST_UPPERBOUND 5000000

/* ----------------------------------------------------------------------- */ 
/*                                Functions                                */
/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------- */ 
/*  Expression_GetSpecificProbe_Main()                                     */
/*  Get specific probe from the expression data.                           */
/*  The first row of the raw data will be ignored.                         */
/*  Return PROC_SUCCESS if success.                                        */
/* ----------------------------------------------------------------------- */ 
int Expression_GetSpecificProbe_Main(char strDatabasePath[], 
			char strInputPath[], int nColumn, char strOutputPath[]);

/* ----------------------------------------------------------------------- */ 
/*  Expression_GetNonRedundantProbe_Main()                                 */
/*  Get non-redundant probe from the expression data.                      */
/*  The first row of the raw data will be ignored.                         */
/*  Return PROC_SUCCESS if success.                                        */
/* ----------------------------------------------------------------------- */ 
int Expression_GetNonRedundantProbe_Main(char strInputPath[], char strOutputPath[]);

/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneSelection_Main()                                        */
/*  Gene selection based on criteria specified in the criteria file        */
/*  Return PROC_SUCCESS if success.                                        */
/* ----------------------------------------------------------------------- */ 
int Expression_GeneSelection_Main(char strDataFile[], char strGeneInfoFile[],
								  char strCriteriaFile[], char strOutFile[]);

/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneSelection_LoadGroupId()                                 */
/*  Load group ids                                                         */
/* ----------------------------------------------------------------------- */ 
struct INTMATRIX *Expression_GeneSelection_LoadGroupId(char strInLine[]);

/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneSelection_MonteCarlo()                                  */
/*  Get a score for every gene.                                            */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *Expression_GeneSelection_MonteCarlo(int nProbeNum, 
					int nArrayNum, struct DOUBLEMATRIX **vArray,
					int nClassNum, struct INTMATRIX *pClassSize, 
					struct INTMATRIX *pClassID,
					int nVargroupNum, struct INTMATRIX *pVargroupSize, 
					struct INTMATRIX **vVargroupMap, 
					int nCycPermNum, char strComparisons[]);

/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneSelection_t-Test()                                      */
/*  Get a score for every gene.                                            */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *Expression_GeneSelection_tTest(int nProbeNum, 
					int nArrayNum, struct DOUBLEMATRIX **vArray,
					int nClassNum, struct INTMATRIX *pClassSize, 
					struct INTMATRIX *pClassID,
					int nVargroupNum, struct INTMATRIX *pVargroupSize, 
					struct INTMATRIX **vVargroupMap, 
					char strComparisons[]);

/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneSelection_SubtractMean()                                */
/*  Get a score for every gene.                                            */
/* ----------------------------------------------------------------------- */ 
int Expression_GeneSelection_SubtractMean(int nProbeNum, 
					int nArrayNum, struct DOUBLEMATRIX **vArray,
					int nClassNum, struct INTMATRIX *pClassSize, 
					struct INTMATRIX *pClassID);

/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneSelection_Evaluate()                                    */
/*  Evaluate a expression.                                                 */
/*  return 1 if the expression is True; 0 if the expression is false.      */
/*  if nSimple = 0, evaluate () first; if nSimple = 1, no () evaluation.   */ 
/* ----------------------------------------------------------------------- */ 
int Expression_GeneSelection_Evaluate(struct DOUBLEMATRIX *pDraws, 
									  char vLogic[], double vGid[], 
									  int nLogicLen, int nSimple);


/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneSelection_EvaluateVec()                                 */
/*  Evaluate a expression.                                                 */
/*  return 0 if any error happened.                                        */
/*  return 1 if the result is a value vector pointer which is not new and  */
/*  shouldn't be deleted.                                                  */
/*  return 2 if the result is a value vector pointer which is new and      */
/*  need to be deleted later.                                              */
/*  return 3 if the result is a logic vector pointer which is new and      */
/*  need to be deleted later.                                              */
/* ----------------------------------------------------------------------- */ 
int Expression_GeneSelection_EvaluateVec(int nProbeNum,
									  char vLogic[], 
									  struct DOUBLEMATRIX *vVid[], 
									  struct BYTEMATRIX *vLid[],
									  int nLogicLen, int nSimple,
									  struct BYTEMATRIX **pResultVec, 
									  struct DOUBLEMATRIX **pResultVal);

/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneSelection_ClassPerm()                                   */
/*  permute class labels.                                                  */
/*  return permuted class labels.                                          */
/* ----------------------------------------------------------------------- */ 
struct INTMATRIX *Expression_GeneSelection_ClassPerm(int nClassNum,
				struct INTMATRIX *pClassID, struct INTMATRIX *pClassSize, 
				int nPermgroupNum, struct INTMATRIX *pPermgroupSize, 
				struct INTMATRIX *vPermgroupMap[]);

/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneSelection_Output()                                      */
/*  Output results.                                                        */
/* ----------------------------------------------------------------------- */ 
int Expression_GeneSelection_Output(int nProbeNum, int nOutputNum, 
				struct DOUBLEMATRIX *pScore, struct LONGMATRIX *pSortID, 
				struct DOUBLEMATRIX *pSortScore, struct DOUBLEMATRIX *pFDR,
				char strProbeFile[], char strGeneInfoFile[], char strOutFile[]);


/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneSelection_Output_1()                                    */
/*  Output results.                                                        */
/* ----------------------------------------------------------------------- */ 
int Expression_GeneSelection_Output_1(int nProbeNum, int nOutputNum, 
				struct DOUBLEMATRIX *pScore, struct LONGMATRIX *pSortID, 
				struct DOUBLEMATRIX *pSortScore, struct DOUBLEMATRIX *pFDR,
				char strProbeFile[], char strGeneInfoFile[], char strOutFile[]);

/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneSelection_Output_WithRandomControl()                    */
/*  Output results, a group of random control will be selected.            */
/* ----------------------------------------------------------------------- */ 
int Expression_GeneSelection_Output_WithRandomControl(int nProbeNum, int nOutputNum, 
				struct DOUBLEMATRIX *pScore, struct LONGMATRIX *pSortID, 
				struct DOUBLEMATRIX *pSortScore, struct DOUBLEMATRIX *pFDR,
				char strProbeFile[], char strGeneInfoFile[], char strOutFile[]);

/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneSelection_Output_WithRandomControl_Fast()               */
/*  Fast version of Expression_GeneSelection_Output_WithRandomControl.     */
/*  Requires annotation ordering to be the same as probeset odering.       */
/*  Output results, a group of random control will be selected.            */
/* ----------------------------------------------------------------------- */ 
int Expression_GeneSelection_Output_WithRandomControl_Fast(int nProbeNum, int nOutputNum, 
				struct DOUBLEMATRIX *pScore, struct LONGMATRIX *pSortID, 
				struct DOUBLEMATRIX *pSortScore, struct DOUBLEMATRIX *pFDR,
				char strProbeFile[], char strGeneInfoFile[], char strOutFile[]);

/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneRankByLocusLink_Main()                                  */
/*  Gene selection based on combining multiple probesets by Locuslink ID   */
/*  Return PROC_SUCCESS if success.                                        */
/* ----------------------------------------------------------------------- */ 
int Expression_GeneRankByLocusLink_Main(char strParamFile[]);

/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneRankByLocusLink()                                       */
/*  Gene selection based on combining multiple probesets by Locuslink ID   */
/*  Return PROC_SUCCESS if success.                                        */
/* ----------------------------------------------------------------------- */ 
int Expression_GeneRankByLocusLink(int nProbeNum, char strScorePath[],
								   char strTransform[], double dResolution,
								   char strMapPath[], int nMaxIter, 
								   double dPrecision, char strOutPath[],
								   int nOutputNum);

/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneRankByLocusLink()                                       */
/*  Gene selection based on combining multiple probesets by Locuslink ID   */
/*  Return PROC_SUCCESS if success.                                        */
/* ----------------------------------------------------------------------- */ 
int Expression_GeneRankByLocusLink_TestPaper(int nProbeNum, char strScorePath[],
								   char strTransform[], double dResolution,
								   char strMapPath[], int nMaxIter, 
								   double dPrecision, char strOutPath[],
								   int nOutputNum);

/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneRankByLocusLink_ScoreTransform()                        */
/*  transform original scores to working scores so that normality is       */
/*  appropriate. Acceptable transformations are logit, rank, identity      */
/*  return the transformed matrix.                                         */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *Expression_GeneRankByLocusLink_ScoreTransform(struct DOUBLEMATRIX *pScore, 
								  double dResolution, char strTransform[]);

/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneRankByLocusLink_ScoreInverseTransform()                 */
/* ----------------------------------------------------------------------- */ 
int Expression_GeneRankByLocusLink_ScoreInverseTransform(struct DOUBLEMATRIX *pScore, 
														 char strTransform[]);

/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneRankByLocusLink_EstimateVar()                           */
/*  estimate the variance                                                  */
/* ----------------------------------------------------------------------- */ 
int Expression_GeneRankByLocusLink_EstimateVar(struct DOUBLEMATRIX *pS, 
							struct INTMATRIX *pN);

/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneRankByLocusLink_EstimateVar0()                           */
/*  estimate the variance                                                  */
/* ----------------------------------------------------------------------- */ 
int Expression_GeneRankByLocusLink_EstimateVar0(struct DOUBLEMATRIX *pS, 
							struct INTMATRIX *pN);

/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneRankByLocusLink_EstimateMean()                          */
/*  estimate the variance                                                  */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *Expression_GeneRankByLocusLink_EstimateMean(struct DOUBLEMATRIX *pX, 
							struct DOUBLEMATRIX *pS, struct INTMATRIX *pN,
							int nMaxIter, double dPrecision);

/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneRankByLocusLink_Output()                                */
/*  output the selected gene                                               */
/* ----------------------------------------------------------------------- */ 
int Expression_GeneRankByLocusLink_Output(int nLocusNum, 
		struct DOUBLEMATRIX *pSortScore, struct LONGMATRIX *pSortIndex, 
		struct INTMATRIX *pLocus, struct INTMATRIX **vLocusMap, 
		int nProbeNum, struct tagString **vProbeName, 
		struct DOUBLEMATRIX *pScore,
		char strMapPath[], char strOutPath[], int nOutputNum);


/* ----------------------------------------------------------------------- */ 
/*  Expression_PowExpressLocusLink_Main()                                  */
/*  This is the algorithm used in powerexpress.                            */
/*  Gene selection based on combining multiple probesets by Locuslink ID   */
/*  Return PROC_SUCCESS if success.                                        */
/* ----------------------------------------------------------------------- */ 
int Expression_PowExpressLocusLink_Main(char strParamFile[]);

/* ----------------------------------------------------------------------- */ 
/*  Expression_PowExpressLocusLink()                                       */
/*  This is the algorithm used by PowerExpress.                            */
/*  Gene selection based on combining multiple probesets by Locuslink ID   */
/*  Return PROC_SUCCESS if success.                                        */
/* ----------------------------------------------------------------------- */ 
int Expression_PowExpressLocusLink(int nProbeNum, char strScorePath[],
								   char strTransform[], double dResolution,
								   char strMapPath[], char strOutPath[],
								   int nOutputNum);

/* ----------------------------------------------------------------------- */ 
/*  Expression_Normalization_Quantile_Main()                               */
/*  Quantile normalization.                                                */
/*  Return PROC_SUCCESS if success.                                        */
/* ----------------------------------------------------------------------- */ 
int Expression_Normalization_Quantile_Main(int nArrayNum, int nProbeNum,
								char strDataFile[], char strOutFile[],
								int nTakeLog, double dTruncLow);

/* ----------------------------------------------------------------------- */ 
/*  Expression_Normalization_Quantile()                                    */
/*  Quantile normalization.                                                */
/*  Return PROC_SUCCESS if success.                                        */
/* ----------------------------------------------------------------------- */ 
int Expression_Normalization_Quantile(int nArrayNum, int nProbeNum, 
								struct DOUBLEMATRIX **vArray);

/* ----------------------------------------------------------------------- */ 
/*  Tiling_ProbeSelection_Main()                                           */
/*  Probe selection based on criteria specified in the criteria file       */
/*  Return PROC_SUCCESS if success.                                        */
/* ----------------------------------------------------------------------- */ 
int Tiling_ProbeSelection_Main(char strDataFile[], char strGeneInfoFile[],
							   char strCriteriaFile[], char strOutFile[]);


/* ----------------------------------------------------------------------- */ 
/*  Tiling_ProbeSelection_Output_WithRandomControl()                       */
/*  Output results, a group of random control will be selected.            */
/* ----------------------------------------------------------------------- */ 
int Tiling_ProbeSelection_Output_WithRandomControl(int nProbeNum, int nOutputNum, 
				struct DOUBLEMATRIX *pScore, struct LONGMATRIX *pSortID, 
				struct DOUBLEMATRIX *pSortScore, struct DOUBLEMATRIX *pFDR,
				char strProbeFile[], char strGeneInfoFile[], char strOutFile[]);

/* ----------------------------------------------------------------------- */ 
/*  Tiling_BindingRegionSelection_HMM_Main()                               */
/*  HMM Forward-Backward algorithm for calculating the posterior prob. of  */
/*  binding for each position.                                             */
/* ----------------------------------------------------------------------- */ 
int Tiling_BindingRegionSelection_HMM_Main(int nProbeNum, char strDataPath[],
				char strTransformType[], double dScoreResolution,
				char strTransitionPath[], char strEmissionPath[], 
				double dGapDist, double dPosteriorCutoff,
				char strOutPath[]);

/* ----------------------------------------------------------------------- */ 
/*  Tiling_BindingRegionSelection_HMM()                                    */
/*  HMM Forward-Backward algorithm for calculating the posterior prob. of  */
/*  binding for each position.                                             */
/* ----------------------------------------------------------------------- */ 
int Tiling_BindingRegionSelection_HMM(int nProbeNum, int nStateNum,
				struct DOUBLEMATRIX *pStationary, struct DOUBLEMATRIX *pTransition, 
				struct DOUBLEMATRIX *pEmission, double dGapDist, 
				struct DOUBLEMATRIX *pScore, struct DOUBLEMATRIX *pPosition, 
				struct DOUBLEMATRIX **vPosterior);

/* ----------------------------------------------------------------------- */ 
/*  Tiling_BindingRegionSelection_ScoreTransform()                         */
/*  transform original scores to working scores                            */
/*  appropriate. Acceptable transformations are logit, rank, identity      */
/*  return the transformed matrix.                                         */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *Tiling_BindingRegionSelection_ScoreTransform(struct DOUBLEMATRIX *pScore, 
										double dResolution, char strTransform[]);

/* ----------------------------------------------------------------------- */ 
/*  Tiling_BindingRegionSelection_HMM_GetTransition()                      */
/*  Get transition probability for HMM                                     */
/* ----------------------------------------------------------------------- */ 
double Tiling_BindingRegionSelection_HMM_GetTransition(struct DOUBLEMATRIX *pStationary,
					struct DOUBLEMATRIX *pTransition, 
					int nFromS, int nToS, double dDist, double dGapDist);

/* ----------------------------------------------------------------------- */ 
/*  Tiling_BindingRegionSelection_HMM_GetEmission()                        */
/*  Get emission probability for HMM                                       */
/* ----------------------------------------------------------------------- */ 
double Tiling_BindingRegionSelection_HMM_GetEmission(struct DOUBLEMATRIX *pEmission, 
					int nState, double dScore, int nEqualLenInterval, 
					double dIntS, double dIntE, double dIntStep);

/* ----------------------------------------------------------------------- */ 
/*  Tiling_ProbeSelection_t-Test()                                         */
/*  Get a score for every probe.                                           */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *Tiling_ProbeSelection_tTest(int nProbeNum, 
					int nArrayNum, struct DOUBLEMATRIX **vArray,
					int nClassNum, struct INTMATRIX *pClassSize, 
					struct INTMATRIX *pClassID,
					int nVargroupNum, struct INTMATRIX *pVargroupSize, 
					struct INTMATRIX **vVargroupMap, 
					char strComparisons[]);

/* ----------------------------------------------------------------------- */ 
/*  Tiling_ProbeSelection_MonteCarlo()                                     */
/*  Get a score for every probe.                                           */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *Tiling_ProbeSelection_MonteCarlo(int nProbeNum, 
					int nArrayNum, struct DOUBLEMATRIX **vArray,
					int nClassNum, struct INTMATRIX *pClassSize, 
					struct INTMATRIX *pClassID,
					int nVargroupNum, struct INTMATRIX *pVargroupSize, 
					struct INTMATRIX **vVargroupMap, 
					int nCycPermNum, char strComparisons[]);

/* ----------------------------------------------------------------------- */ 
/*  Tiling_ProbeSelection_Output_WithRandomControl_ForTest()               */
/*  Output results, a group of random control will be selected.            */
/* ----------------------------------------------------------------------- */ 
int Tiling_ProbeSelection_Output_WithRandomControl_ForTest(int nProbeNum, int nOutputNum, 
				struct DOUBLEMATRIX *pScore, struct LONGMATRIX *pSortID, 
				struct DOUBLEMATRIX *pSortScore, struct DOUBLEMATRIX *pFDR,
				char strProbeFile[], char strGeneInfoFile[], char strOutFile[]);

/* ----------------------------------------------------------------------- */ 
/*  Tiling_BindingRegionSelection_OutputToBed()                            */
/*  Output binding region to a bed file.                                   */
/* ----------------------------------------------------------------------- */ 
int Tiling_BindingRegionSelection_OutputToBed(int nProbeNum, 
					struct DOUBLEMATRIX *pPosition, struct DOUBLEMATRIX **vPosterior, 
					int nStateId, double dPosteriorCutoff, double dGapDist, 
					char strOutPath[]);

/* ----------------------------------------------------------------------- */ 
/*  Tiling_BindingRegionSelection_HMM_BaumWelch_Main()                     */
/*  HMM Forward-Backward algorithm for calculating the posterior prob. of  */
/*  binding for each position. Baum-Welch algorithm is used to estimate    */
/*  unknown parameters.                                                    */
/* ----------------------------------------------------------------------- */ 
int Tiling_BindingRegionSelection_HMM_BaumWelch_Main(int nProbeNum, char strDataPath[],
				char strTransformType[], double dScoreResolution,
				char strTransitionPath[], char strEmissionPath[], 
				double dGapDist, double dPosteriorCutoff,
				double dStopCut, int nMaxIter,
				char strOutPath[]);

/* ----------------------------------------------------------------------- */ 
/*  Tiling_BindingRegionSelection_HMM_BaumWelch()                          */
/*  HMM Baum-Welch algorithm for estimating unknown parameters.            */
/* ----------------------------------------------------------------------- */ 
int Tiling_BindingRegionSelection_HMM_BaumWelch(int nProbeNum, int nStateNum,
				struct DOUBLEMATRIX *pStationary, struct DOUBLEMATRIX *pTransition, 
				struct DOUBLEMATRIX *pEmission, double dGapDist, 
				struct DOUBLEMATRIX *pScore, struct DOUBLEMATRIX *pPosition, 
				struct DOUBLEMATRIX *pStationaryNew, struct DOUBLEMATRIX *pTransitionNew, 
				struct DOUBLEMATRIX *pEmissionNew, struct DOUBLEMATRIX *pStationaryTrue,
				double dStopCut, int nMaxIter);

/* ----------------------------------------------------------------------- */ 
/*  Tiling_BindingRegionSelection_HMM_GetEmission()                        */
/*  Get emission probability for HMM                                       */
/* ----------------------------------------------------------------------- */ 
int Tiling_BindingRegionSelection_HMM_GetEmissionId(struct DOUBLEMATRIX *pEmission, 
					double dScore, int nEqualLenInterval, 
					double dIntS, double dIntE, double dIntStep);

/* ----------------------------------------------------------------------- */ 
/*  Tiling_BindingRegionSelection_HMM_ExpLen_Main()                        */
/*  HMM Forward-Backward algorithm for calculating the posterior prob. of  */
/*  binding for each position. Transition probability is gap length        */
/*  dependent and exponetially decaying.                                   */
/* ----------------------------------------------------------------------- */ 
int Tiling_BindingRegionSelection_HMM_ExpLen_Main(int nProbeNum, char strDataPath[],
				char strTransformType[], double dScoreResolution,
				char strTransitionPath[], char strEmissionPath[], 
				double dGapDist, double dPosteriorCutoff,
				char strOutPath[]);

/* ----------------------------------------------------------------------- */ 
/*  Tiling_BindingRegionSelection_HMM_ExpLen()                             */
/*  HMM Forward-Backward algorithm for calculating the posterior prob. of  */
/*  binding for each position. Transition probability is gap length        */
/*  dependent and exponentially decaying.                                  */
/* ----------------------------------------------------------------------- */ 
int Tiling_BindingRegionSelection_HMM_ExpLen(int nProbeNum, int nStateNum,
				struct DOUBLEMATRIX *pStationary, struct DOUBLEMATRIX *pTransition, 
				struct DOUBLEMATRIX *pEmission, double dGapDist, 
				struct DOUBLEMATRIX *pScore, struct DOUBLEMATRIX *pPosition, 
				struct DOUBLEMATRIX **vPosterior);

/* ----------------------------------------------------------------------- */ 
/*  Tiling_BindingRegionSelection_HMM_GetTransition_ExpLen()               */
/*  Get transition probability for HMM. Transition probability is gap      */
/*  length dependent and exponentially decaying.                           */
/* ----------------------------------------------------------------------- */ 
double Tiling_BindingRegionSelection_HMM_GetTransition_ExpLen(struct DOUBLEMATRIX *pStationary,
					struct DOUBLEMATRIX *pTransition, 
					int nFromS, int nToS, double dDist, double dGapDist);

/* ----------------------------------------------------------------------- */ 
/*  Tiling_BindingRegionSelection_HMM_ConstLen_Main()                      */
/*  HMM Forward-Backward algorithm for calculating the posterior prob. of  */
/*  binding for each position. The binding region has a fixed length.      */
/* ----------------------------------------------------------------------- */ 
int Tiling_BindingRegionSelection_HMM_ConstLen_Main(int nProbeNum, char strDataPath[],
				char strTransformType[], double dScoreResolution,
				char strTransitionPath[], char strEmissionPath[], 
				double dGapDist, double dPosteriorCutoff,
				char strOutPath[]);

/* ----------------------------------------------------------------------- */ 
/*  Tiling_BindingRegionSelection_HMM_ConstLen()                           */
/*  HMM Forward-Backward algorithm for calculating the posterior prob. of  */
/*  binding for each position. The length of the binding region is fixed.  */
/* ----------------------------------------------------------------------- */ 
int Tiling_BindingRegionSelection_HMM_ConstLen(int nProbeNum, int nStateNum,
				struct DOUBLEMATRIX *pStationary, struct DOUBLEMATRIX *pTransition, 
				struct DOUBLEMATRIX *pEmission, double dGapDist, 
				struct DOUBLEMATRIX *pScore, struct DOUBLEMATRIX *pPosition, 
				struct DOUBLEMATRIX **vPosterior);

/* ----------------------------------------------------------------------- */ 
/*  Tiling_UMS_FDR_Main()                                                  */
/*  UMS for FDR estimation.                                                */
/* ----------------------------------------------------------------------- */ 
int Tiling_UMS_FDR_Main(char strSelectPath[], char strScorePath[], 
						double dPcut, double dQcut,
						int nStepSize, int nIntervalNum,
						char strOutPath[]);

/* ----------------------------------------------------------------------- */ 
/*  Tiling_UMS_FDR()                                                       */
/*  UMS for FDR estimation.                                                */
/* ----------------------------------------------------------------------- */ 
int Tiling_UMS_FDR(struct DOUBLEMATRIX *pSelect, double dPcut, double dQcut,
				   int nStepSize, struct DOUBLEMATRIX *pScore, 
				   int nIntervalNum, struct DOUBLEMATRIX *pFDR);

/* ----------------------------------------------------------------------- */ 
/*  TileMap_ImportAffy_Main()                                              */
/*  TileMap loading data from affymetrix's *.CEL files                     */
/* ----------------------------------------------------------------------- */ 
int TileMap_ImportAffy_Main(char strParamPath[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMap_ImportCELIntensity()                                           */
/*  TileMap loading data from a single affymetrix's *.CEL file             */
/* ----------------------------------------------------------------------- */ 
int TileMap_ImportCELIntensity(char strWorkPath[], char strBpmapPath[], 
			char strCelFile[], char strAlias[],
			int nIntensityType, double dLowerBound, int nLogTransform);

/* ----------------------------------------------------------------------- */ 
/*  TileMap_CombineIntensity()                                             */
/*  TileMap combine intensities into one file.                             */
/* ----------------------------------------------------------------------- */ 
int TileMap_CombineIntensity(int nArrayNum, int nTotalProbeNum, char strWorkPath[], 
							 char strBpmapPath[], struct tagString **vAlias, 
							 char strOutFile[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMap_ImportAffy_Normalization_Main()                                */
/*  TileMap loading data from affymetrix's *.CEL files, do normalizations, */
/*  and compute intensities.                                               */
/* ----------------------------------------------------------------------- */ 
int TileMap_ImportAffy_Normalization_Main(char strParamPath[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMap_LoadCEL()                                                      */
/*  TileMap loading raw data from a single affymetrix's *.CEL file         */
/* ----------------------------------------------------------------------- */ 
int TileMap_LoadCEL(char strWorkPath[], char strCelFile[], 
					int *pTotalX, int *pTotalY,
					struct DOUBLEMATRIX **pMean, struct DOUBLEMATRIX **pSD);

/* ----------------------------------------------------------------------- */ 
/*  TileMap_ExportAffyIntensity()                                          */
/*  TileMap export *.CEL intensity to a single file according to *.bpmap   */
/* ----------------------------------------------------------------------- */ 
int TileMap_ExportAffyIntensity(int nArrayNum, int nCellNum, int *pnTotalProbeNum, 
								struct DOUBLEMATRIX **vArray, 
								struct tagString **vAlias,
								int nTotalX, int nTotalY, 
								char strBpmapPath[], char strMaskPath[],
								char strOutPath[], int nIntensityType, 
								double dLowerBound, int nLogTransform);

/* ----------------------------------------------------------------------- */ 
/*  TileMap_Normalization_Main()                                           */
/*  TileMap normalization module.                                          */
/* ----------------------------------------------------------------------- */ 
int TileMap_Normalization_Main(char strParamPath[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMap_Normalization_Quantile()                                       */
/*  Quantile normalization for tilemap.                                    */
/*  Return PROC_SUCCESS if success.                                        */
/* ----------------------------------------------------------------------- */ 
int TileMap_Normalization_Quantile(int nArrayNum, int nProbeNum, 
								char strDataFile[], char strOutFile[],
								int nTakeLog, double dTruncLow);

/* ----------------------------------------------------------------------- */ 
/*  TileMap_Main()                                                         */
/*  TileMap pipeline                                                       */
/* ----------------------------------------------------------------------- */ 
int TileMap_Main(char strParamPath[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMap_ProbeSelection_Main()                                          */
/*  Probe selection based on criteria specified in the criteria file       */
/*  Return PROC_SUCCESS if success.                                        */
/* ----------------------------------------------------------------------- */ 
int TileMap_ProbeSelection_Main(char strDataFile[], char strGeneInfoFile[],
							   char strCriteriaFile[], char strOutFile[],
							   int nProbeNum, int nApplyPerm, int *nPermNum,
							   int *nProbeSummaryRange,	double dZeroCut,
							   int nApplyFilter);

/* ----------------------------------------------------------------------- */ 
/*  TileMap_ProbeSelection_t-Test()                                        */
/*  Get a score for every probe.                                           */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *TileMap_ProbeSelection_tTest(int nProbeNum, 
					int nArrayNum, struct DOUBLEMATRIX **vArray,
					int nClassNum, struct INTMATRIX *pClassSize, 
					struct INTMATRIX *pClassID,
					int nVargroupNum, struct INTMATRIX *pVargroupSize, 
					struct INTMATRIX **vVargroupMap, 
					char strComparisons[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMap_ProbeSelection_MonteCarlo()                                    */
/*  Get a score for every probe.                                           */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *TileMap_ProbeSelection_MonteCarlo(int nProbeNum, 
					int nArrayNum, struct DOUBLEMATRIX **vArray,
					int nClassNum, struct INTMATRIX *pClassSize, 
					struct INTMATRIX *pClassID,
					int nVargroupNum, struct INTMATRIX *pVargroupSize, 
					struct INTMATRIX **vVargroupMap, 
					int nCycPermNum, char strComparisons[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMap_ProbeSelection_Output()                                        */
/*  Output results, a group of random control will be selected.            */
/* ----------------------------------------------------------------------- */ 
int TileMap_ProbeSelection_Output(int nProbeNum, struct DOUBLEMATRIX *pScore, 
				char strProbeFile[], char strOutFile[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMap_BPMAPFilter_Pbsum()                                            */
/*  Filter out repeat and bad probes with bpmap informaiton.               */
/* ----------------------------------------------------------------------- */ 
int TileMap_BPMAPFilter_Pbsum(char strInFile[], char strRefFile[], char strOutFile[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMap_UMS_HMM_Main()                                                 */
/*  UMS for estimating HMM parameters.                                     */
/* ----------------------------------------------------------------------- */ 
int TileMap_UMS_HMM_Main(char strSelectPath[], char strScorePath[], 
						int nProbeNum,  int nScoreRange,
						double dPcut, double dQcut,
						int nStepSize, int nIntervalNum,
						int nFragLen, char strOutPath[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMap_UMS()                                                          */
/*  Unbalanced Mixture Subtraction.                                        */
/* ----------------------------------------------------------------------- */ 
int TileMap_UMS(struct DOUBLEMATRIX *pSelect, double dPcut, double dQcut,
				   int nStepSize, struct DOUBLEMATRIX *pScore, 
				   int nIntervalNum, struct DOUBLEMATRIX *pIntx,
				   struct DOUBLEMATRIX *pFhat0, struct DOUBLEMATRIX *pFhat1,
				   double *dTheta);

/* ----------------------------------------------------------------------- */ 
/*  TileMap_BindingRegionSelection_HMM_Main()                              */
/*  HMM Forward-Backward algorithm for calculating the posterior prob. of  */
/*  binding for each position.                                             */
/* ----------------------------------------------------------------------- */ 
int TileMap_BindingRegionSelection_HMM_Main(int nProbeNum, char strDataPath[],
				char strTransformType[], double dScoreResolution,
				char strTransitionPath[], char strEmissionPath[], 
				double dGapDist, double dPosteriorCutoff,
				char strOutPath[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMap_BindingRegionSelection_HMM()                                   */
/*  HMM Forward-Backward algorithm for calculating the posterior prob. of  */
/*  binding for each position.                                             */
/* ----------------------------------------------------------------------- */ 
int TileMap_BindingRegionSelection_HMM(int nProbeNum, int nStateNum,
				struct DOUBLEMATRIX *pStationary, struct DOUBLEMATRIX *pTransition, 
				struct DOUBLEMATRIX *pEmission, double dGapDist, 
				struct DOUBLEMATRIX *pScore, struct INTMATRIX *pChr,
				struct DOUBLEMATRIX *pPosition, struct DOUBLEMATRIX **vPosterior);

/* ----------------------------------------------------------------------- */ 
/*  TileMap_BindingRegionSelection_MA_UMS_Main()                           */
/*  Moving Average algorithm for calling binding regions in TileMap.       */
/*  FDR will be estimated by UMS.                                          */
/* ----------------------------------------------------------------------- */ 
int TileMap_BindingRegionSelection_MA_UMS_Main(int nProbeNum, 
				char strScorePath[], char strWorkPath[], char strProjectTitle[],
				int nMAW, int nScoreRange,	double dZeroCut, char strUMSSelectPath[], 
				double dTp, double dTq, int nOffset, int nIntervalNum, 
				int nMaxGap, double dFDRCut);

/* ----------------------------------------------------------------------- */ 
/*  TileMap_MA()                                                           */
/*  Moving Average of TileMap.                                             */
/* ----------------------------------------------------------------------- */ 
int TileMap_MA(struct DOUBLEMATRIX *pScore, struct INTMATRIX * pChr, 
			   struct DOUBLEMATRIX *pPosition, int nW, 
			   struct DOUBLEMATRIX *pMA);

/* ----------------------------------------------------------------------- */ 
/*  TileMap_BindingRegionSelection_MA_PERM_Main()                          */
/*  Moving Average algorithm for calling binding regions in TileMap.       */
/*  FDR will be estimated by permutation test.                             */
/* ----------------------------------------------------------------------- */ 
int TileMap_BindingRegionSelection_MA_PERM_Main(int nProbeNum, 
				char strScorePath[], char strWorkPath[], char strProjectTitle[],
				int nW, int nScoreRange, double dZeroCut, int nIntervalNum,
				int nPermNum, int nMaxGap, double dFDRCut);

/* ----------------------------------------------------------------------- */ 
/*  TileMap_CountScore()                                                   */
/*  Count empirical distributions of scores.                               */
/* ----------------------------------------------------------------------- */ 
int TileMap_CountScore(struct DOUBLEMATRIX *pPermMA, int nIntervalNum, 
					   struct DOUBLEMATRIX *pF);

/* ----------------------------------------------------------------------- */ 
/*  TileMap_CallBindingRegion_HMM()                                        */
/*  Call binding regions based on HMM posterior probability.               */
/* ----------------------------------------------------------------------- */ 
int TileMap_CallBindingRegion_HMM(int nProbeNum, struct DOUBLEMATRIX *pScore,
				struct INTMATRIX *pChr, struct DOUBLEMATRIX *pPosition, 
				struct tagString **vChrName, double dCutoff, double dGapDist, 
				char strOutPath[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMap_CallBindingRegion_MA()                                         */
/*  Call binding regions based on MA FDR.                                  */
/* ----------------------------------------------------------------------- */ 
int TileMap_CallBindingRegion_MA(int nProbeNum, struct DOUBLEMATRIX *pScore,
				struct INTMATRIX *pChr, struct DOUBLEMATRIX *pPosition, 
				struct tagString **vChrName, double dCutoff, double dGapDist, 
				char strOutPath[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMap_Extract_Main()                                                 */
/*  TileMap extract data for plot.                                         */
/* ----------------------------------------------------------------------- */ 
int TileMap_Extract_Main(char strParamPath[]);

/* ----------------------------------------------------------------------- */ 
/*  TileMap_Extract()                                                      */
/*  TileMap extract data for plot.                                         */
/*  position, probe-sum, hmm, ma-stat, ma-fdr, raw data.                   */
/* ----------------------------------------------------------------------- */ 
int TileMap_Extract(char strWorkPath[], char strProjectTitle[], 
					char strProbeFile[], char strRawData[],
					char strChr[], int nStart, int nEnd);


/* ----------------------------------------------------------------------- */ 
/*  TransLoc_Main()                                                        */
/*  Transloc Main function.                                                */
/*  get clusters of expression units.                                      */
/* ----------------------------------------------------------------------- */ 
int TransLoc_Main(char strParamPath[]);

/* ----------------------------------------------------------------------- */ 
/*  TransLoc_GetProbeCoordinates()                                         */
/*  get probe coordinates.                                                 */
/* ----------------------------------------------------------------------- */ 
int TransLoc_GetProbeCoordinates(char strArrayAnnotPath[], char strWorkPath[], 
								 char strOutputFile[], char strSpecies[]);

/* ----------------------------------------------------------------------- */ 
/*  TransLoc_FilterProbes()                                                */
/*  filter probes and sort them along chromosomes.                         */
/* ----------------------------------------------------------------------- */ 
int TransLoc_FilterProbes(char strWorkPath[], char strOutputFile[], 
						  char strArrayDataPath[], int nArrayNum,
						  double dTruncLow, int nTakeLog, double dCVCutoff,
						  int *pGeneNum, int *pFinalProbeNum);

/* ----------------------------------------------------------------------- */ 
/*  TransLoc_ComputeCV()                                                   */
/*  Compute coefficient of variation.                                      */
/* ----------------------------------------------------------------------- */ 
double TransLoc_ComputeCV(struct DOUBLEMATRIX *pData);

/* ----------------------------------------------------------------------- */ 
/*  TransLoc_ImportScore()                                                 */
/*  import scores and sort them along chromosomes.                         */
/* ----------------------------------------------------------------------- */ 
int TransLoc_ImportScore(char strWorkPath[], char strOutputFile[], char strArrayDataPath[]);

/* ----------------------------------------------------------------------- */ 
/*  TransLoc_GroupRefGene()                                                */
/*  group refgene annotations.                                             */
/* ----------------------------------------------------------------------- */ 
int TransLoc_GroupRefGene(char strGenomeAnnotPath[], char strWorkPath[], 
						  char strOutputFile[], char strSpecies[]);

/* ----------------------------------------------------------------------- */ 
/*  TransLoc_NeighborDistance()                                            */
/*  compute distance of neighboring genes.                                 */
/* ----------------------------------------------------------------------- */ 
int TransLoc_NeighborDistance(char strWorkPath[], char strOutputFile[], 
		int nArrayNum, int nGeneNum, int nProbeNum);

/* ----------------------------------------------------------------------- */ 
/*  TransLoc_ComputeCorrelation()                                          */
/*  compute correlation of two expression vectors.                         */
/* ----------------------------------------------------------------------- */ 
double TransLoc_ComputeCorrelation(struct DOUBLEMATRIX *pVec1, struct DOUBLEMATRIX *pVec2);
