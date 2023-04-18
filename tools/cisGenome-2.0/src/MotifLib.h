/* ----------------------------------------------------------------------- */
/*  MotifLib.h : interface of the motif library                            */
/*  Author : Ji HongKai ; Time: 2004.07                                    */
/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------- */
/*                              Declaration                                */
/* ----------------------------------------------------------------------- */
#define FASTA_LINE_LEN 60
#define BGMC_MAX_ORDER 10

/* MAX_HMMSTATUS_NUM = 2*(MAX_MOTIF_NUM*2+1) */
#define MAX_MOTIF_NUM 10
#define MAX_HMMSTATUS_NUM 42
#define MAX_TYPE_NUM 12

/* motif probability in the end of sequences */
#define MARGIN_MOTIF_PROB 1e-10
#define MODULE_MOTIF_DENSITY_LOWERBOUND 0.02
#define BG_MOTIF_DENSITY_UPPERBOUND 0.001

#define MOTIF_CALL_CUTOFF 0.5
#define MODULE_SAMPLE_NUM 100

#define QUALITY_HALFWINDOW 50
#define QUALITY_OVERLAP_PENALTY 25

/* known motif mapping */
#define MOTIFMAP_MAXMTFNUM 200 
#define MOTIFMAP_CONSERVE_CUTOFF 150
#define MAX_MOTIF_LEN 10000
#define MOTIFMAP_GENOMECONTIG_SIZE 1000000

/* flex module */
#define FLEXMODULE_SCOREBIN 256
#define FLEXMODULE_MOTIFLEN_LOWERBOUND 6
#define FLEXMODULE_CUMPROB_WID 4

/* ----------------------------------------------------------------------- */ 
/*                                 Struct                                  */
/* ----------------------------------------------------------------------- */ 
struct tagSequence;

/* ----------------------------------------------------------------------- */ 
/*              MOTIFMATRIX: for recording PWM of motifs                   */
/* ----------------------------------------------------------------------- */ 
struct MOTIFMATRIX
{
	/* motif identity */
	int nMotifType;

	/* prior probability for the whole motif */
	double dPt;
	/* prior probability for the + strand motif */
	double dPp;
	/* total count of occurence */
	double dCt;
	/* count of + strand occurence */
	double dCp;
	
	/* nx4 matrices: ROWS are positions, COLUMNs are A,C,G,T's */
	/* position specific weight matrix for + strand. */
	struct DOUBLEMATRIX *pPWMp;
	/* position specific weight matrix for - strand. */
	struct DOUBLEMATRIX *pPWMn;
	/* sample count matrix for + strand. */
	struct DOUBLEMATRIX *pCOUNTp;
	/* sample count matrix for - strand. */
	struct DOUBLEMATRIX *pCOUNTn;

	/* for constructing linear list */
	struct MOTIFMATRIX *pNext;
};

/* ----------------------------------------------------------------------- */ 
/*                MOTIFSITE: for recording motif sites.                    */
/* ----------------------------------------------------------------------- */ 
struct MOTIFSITE
{
	/* motif identity */
	int nMotifType;
	/* sequence id */
	int nSeqId;
	/* start position */
	int nStartPos;
	/* length */
	int nMotifLen;
	/* strand: 0: +; 1: -. */
	int nStrand;
	/* score for the motif */
	double dScore;

	/* for constructing linear list */
	struct MOTIFSITE *pNext;
};

/* ----------------------------------------------------------------------- */ 
/*           FLEXMOTIFMATRIX: for recording PWM of motifs                  */
/* ----------------------------------------------------------------------- */ 
struct FLEXMOTIFMATRIX
{
	/* motif identity */
	int nMotifType;

	/* prior probability for the + strand motif */
	double dPp;
	/* prior probability for the - strand motif */
	double dPm;

	/* count of + strand occurence */
	double dCp;
	/* count of - strand occurence */
	double dCm;
	
	/* nx4 matrices: ROWS are positions, COLUMNs are A,C,G,T's */
	/* prior count matrix for + strand. */
	struct DOUBLEMATRIX *pPriorCount;
	/* sample count matrix for + strand. */
	struct DOUBLEMATRIX *pSampleCount;
	/* position specific weight matrix for + strand. */
	struct DOUBLEMATRIX *pPWM;

	/* for storing extra prior counts */
	/* extra prior for the beginning of the motif */
	int nHeadExtraNum;
	struct DOUBLEMATRIX *pPriorHeadExtra;
	/* extra prior for the end of the motif */
	int nTailExtraNum;
	struct DOUBLEMATRIX *pPriorTailExtra;
};

/* ----------------------------------------------------------------------- */ 
/*                FLEXMOTIFSITE: for recording motif sites.                */
/* ----------------------------------------------------------------------- */ 
struct FLEXMOTIFSITE
{
	/* motif identity */
	int nMotifType;
	/* sequence id */
	int nSeqId;
	/* start position */
	int nStartPos;
	/* strand: 0: +; 1: -. */
	int nStrand;
	/* score for the motif */
	double dScore;
	/* posterior probability of the site */
	double dProb;

	/* context matrix */
	struct INTMATRIX *pContextPos;
	struct INTMATRIX *pContextMotif;
	struct INTMATRIX *pContextStrand;


	/* for constructing bi-linear list */
	struct FLEXMOTIFSITE *pPrev;
	struct FLEXMOTIFSITE *pNext;
};

/* ----------------------------------------------------------------------- */ 
/*         FLEXSEQMOTIF: for organizing sequences for flexmodule.          */
/* ----------------------------------------------------------------------- */ 
struct FLEXSEQMOTIF
{
	/* seqmotif id */
	int nId;
	/* seqmotif name */
	char strAlias[255];
	/* reference sequence length */
	int nSeqLen;

	/* sequences: human,mouse,dog... crossspecies sequences. */
	/* sequence number */
	int nSeqNum;
	/* sequence ids */
	int *pSeqId;
	/* coded sequences */
	struct BYTEMATRIX **vSeq;

	/* all information below are for reference sequence, have length = nSeqLen */

	/* scores: conservation score, CHIP-Chip score... information scores. */
	/* additional information score sequence number */
	int nScoreNum;
	/* information score type */
	int *pScoreId;
	/* information score */
	struct BYTEMATRIX **vScore;

	/* sampler status: A, B status. */
	/* status type number */
	int nStatusNum;
	/* status type */
	int *pStatusId;
	/* status */
	struct BYTEMATRIX **vStatus;

	/* MC counts */
	/* monitor number */
	int nMonitorNum;
	/* monitor id */
	int *pMonitorId;
	/* monitor count */
	struct DOUBLEMATRIX **vMonitor;

	/* motifs */
	struct FLEXMOTIFSITE **vMotif;
	/* site number */
	int nSiteNum;
	/* site info */
	struct FLEXMOTIFSITE *pMotifList;
	/* trial site info */
	struct FLEXMOTIFSITE *pTrialMotifList;
};

/* ----------------------------------------------------------------------- */ 
/*           SEQMOTIF: for organizing sequences for motif discovery.       */
/* ----------------------------------------------------------------------- */ 
struct SEQMOTIF
{
	/* seqmotif id */
	int nId;

	/* sequences: human,mouse,dog... crossspecies sequences. */
	/* sequence number */
	int nSeqNum;
	/* sequence ids */
	int *pSeqId;
	/* coded sequences */
	struct BYTEMATRIX **ppSeq;

	/* scores: conservation score, CHIP-Chip score... information scores. */
	/* additional information score sequence number */
	int nScoreNum;
	/* information score type */
	int *pScoreId;
	/* information score */
	struct DOUBLEMATRIX **ppScore;

	/* sample paths */
	/* HMM Path*/
	struct DOUBLEMATRIX *pHMM;
	/* recorded path number */
	int nPathNum;
	/* path id */
	int *pPathId;
	/* path sample */
	struct BYTEMATRIX **ppSamplePath;

	/* motifs */
	/* site number */
	int nSiteNum;
	/* site info */
	struct MOTIFSITE *pMotifSite;
};

/* ----------------------------------------------------------------------- */ 
/*     MARKOVCHAINMOVE: for organizing markov chain calculation.           */
/* ----------------------------------------------------------------------- */ 
struct MARKOVCHAINMOVE
{
	/* current base or simulated base */
	int nCurrentBase;
	/* oldest base or simulated base */
	int nFirstPreBase;
	/* previous word */
	int nPreWordId;
	/* is current base simulated or not? 1-simulated; 0-good */
	int nIsBad;
	/* transition probability */
	double dPTrans;
	/* next */
	struct MARKOVCHAINMOVE *pNext;
};

/* ----------------------------------------------------------------------- */ 
/*     SITEPOS: for organizing neighboring sites.                          */
/* ----------------------------------------------------------------------- */ 
struct SITEPOS
{
	/* position */
	int nPos;
	/* next */
	struct SITEPOS *pNext;
};

/* ----------------------------------------------------------------------- */ 
/*  MOTIFSITEGROUP: for recording all motif sites in a region.             */
/* ----------------------------------------------------------------------- */ 
struct MOTIFSITEGROUP
{
	/* Group Id */
	int nGroupId;
	/* number of different motifs */
	int nMotifNum;
	/* vector of motif site lists */
	struct MOTIFSITE **pSites;
};

/* ----------------------------------------------------------------------- */ 
/*                                Functions                                */
/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_Main: flexmodule main function                              */
/*  sample cis-regulatory module with conservation scores.                 */
/* ----------------------------------------------------------------------- */ 
int FlexModule_Main(char strParamFile[]);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_Sampler: flexmodule sampler                                 */
/*  sample cis-regulatory module with conservation scores.                 */
/* ----------------------------------------------------------------------- */ 
int FlexModule_Sampler(char strWorkPath[], char strSeqFile[], char strOutFile[],
		/* motif number K, motif length poisson parameter, motif len distribution bin number */
		int nMotifNum, int nMeanMotifLen, int nMotifLenUpperBound,
		/* init motif lengths, init motif pseudocounts (PWMs) */
		struct INTMATRIX *pMotifLen, struct tagString **vPriorPWMFile, 
		/* init module score parameter d, module length */
		double dInitModuleD, int nModuleLen, int nSampleModuleLen,
		/* order of background markov chain */
		int nBGOrder, int nUseFittedBG, char strFittedBGPrefix[],
		/* MC draw number, using conservation indicator, conservation file prefix, conservation likelihood file */
		int nMCNum, int nUseCS, char strCSPrefix[], char strCSLike[],
		/* prior abundance pseudocounts */
		char strPriorAbundance[]);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_InitScoreLike: prepare likelihood function for scores       */
/*  sample cis-regulatory module with conservation scores.                 */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX **FlexModule_InitScoreLike(struct DOUBLEMATRIX *pScoreLike);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_InitSite: Initialize motif sites at random.                 */
/* ----------------------------------------------------------------------- */ 
int FlexModule_InitSite(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,
			struct DOUBLEMATRIX *pFreqPrior0, struct DOUBLEMATRIX *pFreqCount,
			int nBGOrder, struct DOUBLEMATRIX *pLogBG, struct DOUBLEMATRIX *pLogBG0,  
			int nMotifNum, struct FLEXMOTIFMATRIX **vMotif, 
			struct INTMATRIX *pMotifLen, int nMaxMotifLen,
			int nUseCS, int nCSLikeNum, struct DOUBLEMATRIX **vCSLike);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_InitContext: create initial context.                        */
/* ----------------------------------------------------------------------- */ 
int FlexModule_InitContext(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, double dT,
			struct DOUBLEMATRIX *pFreqCount, double dInitModuleD, 
			double dModuleL, double dModuleA, int nModuleLen,
			int nMotifNum, struct FLEXMOTIFMATRIX **vMotif, 
			struct INTMATRIX *pMotifLen, int nMaxMotifLen);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_SampleBA: sample B and A status in flex module sampler.     */
/* ----------------------------------------------------------------------- */ 
int FlexModule_SampleBA(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, double dT, 
				         struct DOUBLEMATRIX *pFreqCount, 
						 double dModuleD, double dModuleL, double dModuleA, int nModuleLen, 
						 int nBGOrder, struct DOUBLEMATRIX *pLogBG, struct DOUBLEMATRIX *pLogBG0,
						 int nMotifNum, struct FLEXMOTIFMATRIX **vMotif, 
						 struct INTMATRIX *pMotifLen, int nMaxMotifLen, 
						 int nUseCS, int nCSLikeNum, struct DOUBLEMATRIX **vCSLike,
						 int nIsRecording, struct DOUBLEMATRIX *pMotifNSample);


/* ----------------------------------------------------------------------- */ 
/*  FlexModule_SampleL: sample module length.                              */
/* ----------------------------------------------------------------------- */ 
int FlexModule_SampleL(struct FLEXSEQMOTIF **vSeqMtf, int nSeqCount, 
				int nIndex, double dT, 
				double dModuleD, double dModuleL, double dModuleA, 
				int *pModuleLen, int nMinModuleLen, int nMaxModuleLen, 
				int nMotifNum, struct FLEXMOTIFMATRIX **vMotif, 
				struct INTMATRIX *pMotifLen, int nMaxMotifLen, 
				int nIsRecording);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_LocalShift: shift motif sites by 1 or -1 or 0.              */
/* ----------------------------------------------------------------------- */ 
int FlexModule_LocalShift(struct FLEXSEQMOTIF **vSeqMtf, int nSeqCount, 
				int nIndex, double dT, 
				struct DOUBLEMATRIX *pFreqCount, double dModuleD, 
				double dModuleL, double dModuleA, int nModuleLen,
				int nBGOrder, struct DOUBLEMATRIX *pLogBG, struct DOUBLEMATRIX *pLogBG0,
				int nMotifId, struct FLEXMOTIFMATRIX **vMotif, 
				struct INTMATRIX *pMotifLen, int nMaxMotifLen,
				int nUseCS, int nCSLikeNum, struct DOUBLEMATRIX **vCSLike,
				int nIsRecording);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_ClearOldStatus: clear old status of a given position.       */
/* ----------------------------------------------------------------------- */ 
double FlexModule_ShiftContextLogLike(struct INTMATRIX *pPos, struct INTMATRIX *pMotif, 
				struct INTMATRIX *pStrand, struct BYTEMATRIX *pA,
				int nLen, struct INTMATRIX *pMotifLen, 
				int nModuleLen, double dModuleD, 
				double dModuleL, double dModuleA, 
				double dT, int *pMask);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_ShiftMotifLogLike: compute shifted motif loglikelihood      */
/* ----------------------------------------------------------------------- */ 
double FlexModule_ShiftMotifLogLike(struct FLEXMOTIFMATRIX *pMotif, 
							 struct DOUBLEMATRIX *pMCount, int nStep);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_ShiftUpdateMotif: update shifted motif matrix               */
/* ----------------------------------------------------------------------- */ 
int FlexModule_ShiftUpdateMotif(struct FLEXMOTIFMATRIX *pMotif, 
							 struct DOUBLEMATRIX *pMCount, int nStep);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_ShiftUpdate: update all indicators and context after shift  */
/* ----------------------------------------------------------------------- */ 
int FlexModule_ShiftUpdate(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, 
				int nMotifId, int nStep,
				int nModuleLen, struct INTMATRIX *pMotifLen, int nMaxMotifLen);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_SampleW: sample motif width W.                              */
/* ----------------------------------------------------------------------- */ 
int FlexModule_SampleW(struct FLEXSEQMOTIF **vSeqMtf, int nSeqCount, 
				int nIndex, double dT, 
				struct DOUBLEMATRIX *pFreqCount, double dModuleD, 
				double dModuleL, double dModuleA, int nModuleLen,
				int nBGOrder, struct DOUBLEMATRIX *pLogBG, struct DOUBLEMATRIX *pLogBG0,
				int nMotifId, int nMotifNum, struct FLEXMOTIFMATRIX **vMotif, 
				struct INTMATRIX *pMotifLen, int *pMaxMotifLen,
				int nUseCS, int nCSLikeNum, struct DOUBLEMATRIX **vCSLike,
				int nIsRecording, int nMeanMotifLen, 
				int nMotifLenUpperBound, int nMotifLenLowerBound,
				double dDefaultPriorCount);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_IncreaseW: try to increase motif width W by 1.              */
/* ----------------------------------------------------------------------- */ 
int FlexModule_IncreaseW(struct FLEXSEQMOTIF **vSeqMtf, int nSeqCount, 
				int nIndex, double dT, 
				struct DOUBLEMATRIX *pFreqCount, double dModuleD, 
				double dModuleL, double dModuleA, int nModuleLen,
				int nBGOrder, struct DOUBLEMATRIX *pLogBG, struct DOUBLEMATRIX *pLogBG0,
				int nMotifId, int nMotifNum, struct FLEXMOTIFMATRIX **vMotif, 
				struct INTMATRIX *pMotifLen, int *pMaxMotifLen,
				int nUseCS, int nCSLikeNum, struct DOUBLEMATRIX **vCSLike,
				int nIsRecording, int nDirec, int nMeanMotifLen, 
				int nMotifLenUpperBound, int nMotifLenLowerBound,
				double dDefaultPriorCount);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_IncreaseWMotifLogLike: compute length-increased motif's     */
/*  loglikelihood.                                                         */
/* ----------------------------------------------------------------------- */ 
double FlexModule_IncreaseWMotifLogLike(struct FLEXMOTIFMATRIX *pMotif, 
							 struct DOUBLEMATRIX *pMCount, int nDirec,
							 double dDefaultPriorCount);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_IncreaseWUpdate: update all indicators and context after    */
/*  increasing motif length by 1.                                          */
/* ----------------------------------------------------------------------- */ 
int FlexModule_IncreaseWUpdate(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, 
				int nMotifId, int nDirec,
				int nModuleLen, struct INTMATRIX *pMotifLen, int *pMaxMotifLen);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_IncreaseWUpdateMotif: update motif after increasing its     */
/*  length by one.                                                         */
/* ----------------------------------------------------------------------- */ 
int FlexModule_IncreaseWUpdateMotif(struct FLEXMOTIFMATRIX *pMotif, 
							 struct DOUBLEMATRIX *pMCount, int nDirec,
							 double dDefaultPriorCount);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_DecreaseW: try to decrease motif width W by 1.              */
/* ----------------------------------------------------------------------- */ 
int FlexModule_DecreaseW(struct FLEXSEQMOTIF **vSeqMtf, int nSeqCount, 
				int nIndex, double dT, 
				struct DOUBLEMATRIX *pFreqCount, double dModuleD, 
				double dModuleL, double dModuleA, int nModuleLen,
				int nBGOrder, struct DOUBLEMATRIX *pLogBG, struct DOUBLEMATRIX *pLogBG0,
				int nMotifId, int nMotifNum, struct FLEXMOTIFMATRIX **vMotif, 
				struct INTMATRIX *pMotifLen, int *pMaxMotifLen,
				int nUseCS, int nCSLikeNum, struct DOUBLEMATRIX **vCSLike,
				int nIsRecording, int nDirec, int nMeanMotifLen, 
				int nMotifLenUpperBound, int nMotifLenLowerBound,
				double dDefaultPriorCount);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_DecreaseWMotifLogLike: compute length-decreased motif's     */
/*  loglikelihood.                                                         */
/* ----------------------------------------------------------------------- */ 
double FlexModule_DecreaseWMotifLogLike(struct FLEXMOTIFMATRIX *pMotif, 
							 struct DOUBLEMATRIX *pMCount, int nDirec,
							 double dDefaultPriorCount);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_DecreaseWUpdate: update all indicators and context after    */
/*  decreasing motif length by 1.                                          */
/* ----------------------------------------------------------------------- */ 
int FlexModule_DecreaseWUpdate(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, 
				int nMotifId, int nDirec,
				int nModuleLen, struct INTMATRIX *pMotifLen, int *pMaxMotifLen);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_DecreaseWUpdateMotif: update motif after decreasing its     */
/*  length by one.                                                         */
/* ----------------------------------------------------------------------- */ 
int FlexModule_DecreaseWUpdateMotif(struct FLEXMOTIFMATRIX *pMotif, 
							 struct DOUBLEMATRIX *pMCount, int nDirec,
							 double dDefaultPriorCount);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_ClearOldStatus: clear old status of a given position.       */
/* ----------------------------------------------------------------------- */ 
int FlexModule_ClearOldStatus(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, int nPos, 
			struct FLEXMOTIFSITE *pContext, int *pContextNum, int nBaseTypeNum,
			struct FLEXMOTIFSITE **pPrev, struct FLEXMOTIFSITE **pNext, 
			struct DOUBLEMATRIX *pFreqCount, int nHalfWin, 
			int nMotifNum, struct FLEXMOTIFMATRIX **vMotif, 
			struct INTMATRIX *pMotifLen, int nMaxMotifLen);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_UpdateNewStatus: update new context status.                 */
/* ----------------------------------------------------------------------- */ 
int FlexModule_UpdateNewStatus(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, int nPos,
			int nMotifId, int nStrand,
			struct FLEXMOTIFSITE *pContext, int *pContextNum, int nBaseTypeNum,
			struct FLEXMOTIFSITE **pPrev, struct FLEXMOTIFSITE **pNext, 
			int nHalfWin, 
			int nMotifNum, struct FLEXMOTIFMATRIX **vMotif, 
			struct INTMATRIX *pMotifLen, int nMaxMotifLen);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_InitBGLogLike: Initialize background log likelihood         */
/* ----------------------------------------------------------------------- */ 
int FlexModule_InitBGLogLike(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,
			int nBGOrder, struct DOUBLEMATRIX *pLogBG, struct DOUBLEMATRIX *pLogBG0,
			int nUseCS, int nCSLikeNum, struct DOUBLEMATRIX **vCSLike);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_FromFittedBGFile: Initialize background log likelihood from */
/*  fitted background files.                                               */
/* ----------------------------------------------------------------------- */ 
int FlexModule_InitBGLogLike_FromFittedBGFile(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,
			char strFilePath[],	int nUseCS, int nCSLikeNum, struct DOUBLEMATRIX **vCSLike);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_ComputeContextLogLike: compute context log likelihood.      */
/* ----------------------------------------------------------------------- */ 
double FlexModule_ComputeContextLogLike(struct FLEXSEQMOTIF *pSeqMtf, struct FLEXMOTIFSITE *pContext,
					int nPos, int nContextNum, int nAddPos, int nAddMotifId, int nAddStrand,
					int nMotifNum, double dModuleD, double dModuleL, double dModuleA, double dT);

/* ----------------------------------------------------------------------- */ 
/* FlexModule_AddMotifMatrixCount: Add count to motif matrix.              */
/* ----------------------------------------------------------------------- */ 
int FlexModule_AddMotifMatrixCount(struct DOUBLEMATRIX *pCountMat, 
						  unsigned char pSite[], int nMotifLen, char chStrand);

/* ----------------------------------------------------------------------- */ 
/* FlexModule_SubtractMotifMatrixCount: Subtract count to motif matrix.    */
/* ----------------------------------------------------------------------- */ 
int FlexModule_SubtractMotifMatrixCount(struct DOUBLEMATRIX *pCountMat, 
						  unsigned char pSite[], int nMotifLen, char chStrand);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_NormalizeFreq: normalize the frequency matrix.              */
/* ----------------------------------------------------------------------- */ 
int FlexModule_NormalizeFreq(struct DOUBLEMATRIX *pFreq);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_LogFreq: take log-transformation of the frequency matrix.   */
/* ----------------------------------------------------------------------- */ 
int FlexModule_LogFreq(struct DOUBLEMATRIX *pFreq);

/* ----------------------------------------------------------------------- */ 
/*  FlexBaseLogLikelihood_Motif: get loglikelihood for motif.              */
/* ----------------------------------------------------------------------- */ 
double FlexBaseLogLikelihood_Motif(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, int nPos, 
					struct DOUBLEMATRIX *pLogPWM, char chStrand, int nMaxMotifLen,
					int nBGOrder, struct DOUBLEMATRIX *pLogBG, struct DOUBLEMATRIX *pLogBG0, 
					int nUseCS, int nCSLikeNum, struct DOUBLEMATRIX **vCSLike,
					int *pMask);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_NormalizePost: normalize log frequency.                     */
/* ----------------------------------------------------------------------- */ 
int FlexModule_NormalizePost(struct DOUBLEMATRIX *pFreq, struct BYTEMATRIX *pValid);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_SamplePost: sample from a frequency matrix.                 */
/* ----------------------------------------------------------------------- */ 
int FlexModule_SamplePost(struct DOUBLEMATRIX *pFreq, struct BYTEMATRIX *pValid,
						  int *pI, int *pJ);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_CallMotifSitesInFinalSample: call motif site according to   */
/*  the last sample.                                                       */
/* ----------------------------------------------------------------------- */ 
int FlexModule_CallMotifSitesInFinalSample(struct FLEXSEQMOTIF *pSeqMtf, 
					struct INTMATRIX *pMotifLen);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_WriteMotifToFile: write motif sites to file.                */
/* ----------------------------------------------------------------------- */ 
int FlexModule_WriteMotifToFile(char strFilePath[], int nMotifNum, 
					struct FLEXMOTIFMATRIX **vMotif, struct DOUBLEMATRIX *pBG0, 
					struct INTMATRIX *pMotifLen,
					int nSeqCount, struct FLEXSEQMOTIF **vSeqMtf, int nIndex);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_ComputeSiteScore: compute motif score for a site.           */
/* ----------------------------------------------------------------------- */ 
double FlexModule_ComputeSiteScore(unsigned char *pBase, 
					struct DOUBLEMATRIX *pPWM, struct DOUBLEMATRIX *pBG, 
					char chStrand);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_WritePostMotifToFile: write motif sites to file.            */
/* ----------------------------------------------------------------------- */ 
int FlexModule_WritePostMotifToFile(char strFilePath[], int nMotifNum, 
					struct FLEXMOTIFMATRIX **vMotif, 
					struct DOUBLEMATRIX *pBG0, struct DOUBLEMATRIX *pLogBG0, 
					struct INTMATRIX *pMotifLen, struct INTMATRIX *pMeanMotifLen, 
					struct DOUBLEMATRIX *pMotifNSample, 
					int nSeqCount, struct FLEXSEQMOTIF **vSeqMtf, 
					int nIndex, int nRecordNum, double dDefaultPrior);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_BGFit_Main: flexmodule fit background likelihood function   */
/* ----------------------------------------------------------------------- */ 
int FlexModule_BGFit_Main(char strParamFile[]);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_BGFit: flexmodule fit background likelihood                 */
/*  fit a mixture background model, and export log mean likelihood.        */
/* ----------------------------------------------------------------------- */ 
int FlexModule_BGFit(char strWorkPath[], char strSeqFile[], char strOutFile[], int nExportBGFit,
		/* order of background markov chain, motif number K, motif PWM (count matrix) */
		int nBGOrder, int nMotifNum, struct tagString **vPWMFile,
		/* prior count for background components, control motif concentration */
		struct DOUBLEMATRIX *pFreqPrior, int nUsePostEnrich, struct DOUBLEMATRIX *pControlRate,	
		/* MC draw number, using conservation indicator, conservation file prefix, conservation likelihood file */
		int nMCNum, int nUseCS, char strCSPrefix[], char strCSLike[]);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_BGFit_InitMotifLogLike: Initialize motif log likelihood     */
/* ----------------------------------------------------------------------- */ 
int FlexModule_BGFit_InitMotifLogLike(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,
			int nBGOrder, struct DOUBLEMATRIX *pLogBG, struct DOUBLEMATRIX *pLogBG0, 
			int nMotifNum, struct FLEXMOTIFMATRIX **vMotif, struct INTMATRIX *pMotifLen, 
			int nMaxMotifLen,
			int nUseCS, int nCSLikeNum, struct DOUBLEMATRIX **vCSLike);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_BGFit_InitSite: Initialize motif sites at random for        */
/*  background fitting.                                                    */
/* ----------------------------------------------------------------------- */ 
int FlexModule_BGFit_InitSite(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,
			struct DOUBLEMATRIX *pFreqPrior0, struct DOUBLEMATRIX *pFreqCount,
			int nBGOrder, struct DOUBLEMATRIX *pLogBG, struct DOUBLEMATRIX *pLogBG0,  
			int nMotifNum, struct FLEXMOTIFMATRIX **vMotif, 
			struct INTMATRIX *pMotifLen, int nMaxMotifLen,
			int nUseCS, int nCSLikeNum, struct DOUBLEMATRIX **vCSLike);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_BGFit_SampleBA: sample motifs in flexmodule background fit. */
/* ----------------------------------------------------------------------- */ 
int FlexModule_BGFit_SampleBA(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, 
						 int nUsePostEnrich, struct DOUBLEMATRIX *pFreqPrior, 
				         struct DOUBLEMATRIX *pFreqCount,
						 int nBGOrder, struct DOUBLEMATRIX *pLogBG, struct DOUBLEMATRIX *pLogBG0,
						 int nMotifNum, struct FLEXMOTIFMATRIX **vMotif, 
						 struct INTMATRIX *pMotifLen, int nMaxMotifLen, 
						 int nUseCS, int nCSLikeNum, struct DOUBLEMATRIX **vCSLike,
						 int nIsRecording);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_BGFit_ClearOldStatus: clear old status of a given position. */
/* ----------------------------------------------------------------------- */ 
int FlexModule_BGFit_ClearOldStatus(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, int nPos, 
			int nBaseTypeNum,
			struct FLEXMOTIFSITE **pPrev, struct FLEXMOTIFSITE **pNext, 
			struct DOUBLEMATRIX *pFreqCount, 
			int nMotifNum, struct FLEXMOTIFMATRIX **vMotif, 
			struct INTMATRIX *pMotifLen, int nMaxMotifLen);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_BGFit_UpdateNewStatus: update new context status.           */
/* ----------------------------------------------------------------------- */ 
int FlexModule_BGFit_UpdateNewStatus(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, int nPos,
			int nMotifId, int nStrand,
			int nBaseTypeNum,
			struct FLEXMOTIFSITE **pPrev, struct FLEXMOTIFSITE **pNext, 
			int nMotifNum, struct FLEXMOTIFMATRIX **vMotif, 
			struct INTMATRIX *pMotifLen, int nMaxMotifLen);


/* ----------------------------------------------------------------------- */ 
/*  FlexModule_BGFit_UpdateLikelihood: update likelihood.                  */
/* ----------------------------------------------------------------------- */
int FlexModule_BGFit_UpdateLikelihood(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,  
				         struct DOUBLEMATRIX *pTDR, int nBGOrder, 
						 struct DOUBLEMATRIX *pLogBG, struct DOUBLEMATRIX *pLogBG0,
						 int nMotifNum, struct FLEXMOTIFMATRIX **vMotif, 
						 struct INTMATRIX *pMotifLen, int nMaxMotifLen, 
						 int nUseCS, int nCSLikeNum, struct DOUBLEMATRIX **vCSLike);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_BGFit_UpdateLikelihood: update likelihood.                  */
/* ----------------------------------------------------------------------- */
int FlexModule_BGFit_ExportLogLikelihood(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, 
			int nRecordNum, char strFilePath[]);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_BGFit_WritePostMotifToFile: write motif sites to file.      */
/* ----------------------------------------------------------------------- */ 
int FlexModule_BGFit_WritePostMotifToFile(char strFilePath[], int nMotifNum, 
					struct FLEXMOTIFMATRIX **vMotif, 
					struct DOUBLEMATRIX *pBG0, struct DOUBLEMATRIX *pLogBG0, 
					struct INTMATRIX *pMotifLen, struct DOUBLEMATRIX *pFreqCount, 
					int nUsePostEnrich, struct DOUBLEMATRIX *pFreqPrior,
					struct DOUBLEMATRIX *pMotifSiteNum,
					int nSeqCount, struct FLEXSEQMOTIF **vSeqMtf, 
					int nIndex, int nRecordNum);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_WriteCisGenomeIniFiles: prepare cisgenome ini file.         */
/* ----------------------------------------------------------------------- */ 
int FlexModule_WriteCisGenomeIniFiles(char strParamFile[], char strWorkPath[], 
					char strOutFile[], int nMotifNum);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_ExtractMotifsFromResult: prepare cisgenome ini file.        */
/* ----------------------------------------------------------------------- */ 
int FlexModule_ExtractMotifsFromResult(char strInFile[], char strOutFile[]);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanConsensus_Genome_Main: map motif consensus pattern to     */
/*  Genome sequences.                                                      */
/*  nCM = # of mismatches to the consensus allowed.                        */
/*  nDM = # of mismatches to the degenerated consensus allowed.            */
/*  dC  = the cutoff for average conservation scores.                      */
/*  Example: if we are searching for TGGGT[A]GG, and encouter TGGAAGG,     */
/*   then the actual CM is 2, and the actual DM is 1.                      */  
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanConsensus_Genome_Main(char strMotifPath[], 
					char strGenomePath[], char strCodPath[], char strOutputPath[],  
					int nMC, int nMD, int nUseCS, double dC, char strCSPath[],
					int nIncludeRepeat);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanConsensus_Export_Genome: save motif sites to a file.      */
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanConsensus_Export_Genome(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,
				int nMotifLen, char strSeqAlias[], char strChr[], int nOffset,
				FILE *fpOut);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanConsensus_Sequential_Main: map motif consensus pattern to */
/*  FASTA sequences.                                                       */
/*  nCM = # of mismatches to the consensus allowed.                        */
/*  nDM = # of mismatches to the degenerated consensus allowed.            */
/*  dC  = the cutoff for average conservation scores.                      */
/*  Example: if we are searching for TGGGT[A]GG, and encouter TGGAAGG,     */
/*   then the actual CM is 2, and the actual DM is 1.                      */  
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanConsensus_Sequential_Main(char strMotifPath[], 
					char strInputPath[], char strSeqFile[], char strOutputPath[], 
					int nMC, int nMD, int nUseCS, double dC, char strCSPrefix[]);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanConsensus_In_FlexSeqMtf: call motif sites.                */
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanConsensus_In_FlexSeqMtf(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,
			int nMotifLen, struct BYTEMATRIX *pCMotif, struct BYTEMATRIX *pDMotif, 
			int nMC, int nMD, int nUseCS, double dC,
			int *pEffecLen, int *pConsEffecLen, int *pSiteNum);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanConsensus_In_FlexSeqMtf: count mismatches to a motif      */
/*  consensus.                                                             */
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanConsensus_MotifCountMisMatch(struct FLEXSEQMOTIF *pSeqMtf, 
					int nIndex, int nPos, char chStrand, int nMotifLen, 
					struct BYTEMATRIX *pCMotif, struct BYTEMATRIX *pDMotif, 
					int *pMisC, int *pMisD, int *pMask);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanConsensus_Export_Single: save motif sites to a file.      */
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanConsensus_Export_Single(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,
						int nMotifLen, char strSeqAlias[], FILE *fpOut);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_GenomeBackground_Main: compute local markov        */
/*  background matrices. These matrices will be used by genome-wide        */
/*  matrix scan routines.                                                  */ 
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanMatrix_GenomeBackground_Main(char strGenomePath[], 
					char strOutPath[], int nBGOrder, int nS, int nW);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_GenomeBackground_Chr: compute local markov         */
/*  background matrices for a single chromosome. These matrices will be    */
/*  used by genome-wide matrix scan routines.                              */ 
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanMatrix_GenomeBackground_Chr(char strGenomePath[], char strChr[],
			char strExportPath[], int nChrLen, int nS, int nW, int nBGOrder);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_Genome_Main: map motif probability matrix to       */
/*  genomic regions, using nBGOrder markov chain as background and dR as   */
/*  likelihood ratio cutoff. dC is the conservation score cutoff.          */
/*  strBGType specifies whether region-derived or genomic location         */
/*  specific background should be used.                                    */ 
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanMatrix_Genome_Main(char strMotifPath[], char strGenomePath[], 
					char strCodPath[], char strOutputPath[], double dR, 
					int nBGOrder, char strBGType[], 
					char strBGPath[], int nBGStepSize,
					int nUseCS, double dC, char strCSPath[],
					int nIncludeRepeat);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_Genome_ComputeBackground_Region: construct         */
/*  background markov matrices for specified regions.                      */ 
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanMatrix_Genome_ComputeBackground_Region(char strGenomePath[],
				char strCodPath[], int nBGOrder, 
				struct DOUBLEMATRIX *pBGF, struct DOUBLEMATRIX *pBGB, 
				struct DOUBLEMATRIX *pBG0);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_Genome_CountBackground: count empirical            */
/*  observations for constructing background markov matrices.              */ 
/* ----------------------------------------------------------------------- */
int MotifMap_ScanMatrix_Genome_CountBackground(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,
					int nFrom, int nTo, int nBGOrder, 
					struct DOUBLEMATRIX *pBGF, struct DOUBLEMATRIX *pBGB);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_Export_Genome: save motif sites to a file.         */
/* ----------------------------------------------------------------------- */
int MotifMap_ScanMatrix_Export_Genome(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,
						int nMotifLen, char strSeqAlias[], 
						char strChr[], int nOffset, FILE *fpOut);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_Sequential_Main: map motif probability matrix to   */
/*  FASTA sequences, using nBGOrder markov chain as background and dR as   */
/*  likelihood ratio cutoff. dC is the conservation score cutoff.          */ 
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanMatrix_Sequential_Main(char  strMotifPath[], char strInputPath[],
					char strSeqFile[], char strOutputPath[], 
					double dR, int nBGOrder,
					int nUseCS, double dC, char strCSAlias[]);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_Sequential_ComputeBackground: construct background */
/*  markov matrices.                                                       */ 
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanMatrix_Sequential_ComputeBackground(char strInputPath[],
				char strSeqFile[], int nBGOrder, 
				struct DOUBLEMATRIX *pBGF, struct DOUBLEMATRIX *pBGB, 
				struct DOUBLEMATRIX *pBG0);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_CountBackground: count empirical observations for  */
/*  constructing background markov matrices.                               */ 
/* ----------------------------------------------------------------------- */
int MotifMap_ScanMatrix_CountBackground(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,
					int nBGOrder, struct DOUBLEMATRIX *pBGF, struct DOUBLEMATRIX *pBGB);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_Group_Main: map motif probability matrix to        */
/*  FASTA sequences, using nBGOrder markov chain as background and dR as   */
/*  likelihood ratio cutoff. dC is the conservation score cutoff.          */
/*  All the sequences will be loaded into the memory before motif scan,    */
/*  therefore this function is faster but requires more memory than        */
/*  MotifMap_ScanMatrix_Sequential_Main.                                   */  
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanMatrix_Group_Main(char  strMotifPath[], char strInputPath[],
					char strSeqFile[], char strOutputPath[], 
					double dR, int nBGOrder,
					int nUseCS, double dC, char strCSAlias[]);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_InitBGLogLike: Initialize background loglikelihood */
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanMatrix_InitBGLogLike(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,
			int nBGOrder, struct DOUBLEMATRIX *pLogBGF, struct DOUBLEMATRIX *pLogBGB,
			struct DOUBLEMATRIX *pLogBG0);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_InitBGLogLike_Genome: Initialize background        */
/*  loglikelihood using genomic local MC matrices.                         */
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanMatrix_InitBGLogLike_Genome(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,
					char strChr[], int nChrStart, int nChrEnd, 
					int nBGOrder, char strBGPath[], int nBGStepSize);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_In_FlexSeqMtf: call motif sites.                   */
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanMatrix_In_FlexSeqMtf(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,
			int nFrom, int nTo, 
			struct DOUBLEMATRIX *pLogPWM, double dR, 
			int nUseCS, double dC, 
			int *pEffecLen, int *pConsEffecLen, int *pSiteNum);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_Export_Single: save motif sites to a file.         */
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanMatrix_Export_Single(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,
						int nMotifLen, char strSeqAlias[], FILE *fpOut);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_Export: export matrix mapping result.              */
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanMatrix_Export(int nSeqCount, struct FLEXSEQMOTIF **vSeqMtf, int nIndex,
			int nMotifLen, double dTotLen, double dTotSite, char strOutputPath[]);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_MotifLikeRatio: get loglikelihood ratio for motif. */
/* ----------------------------------------------------------------------- */ 
double MotifMap_ScanMatrix_MotifLikeRatio(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, int nPos, 
					struct DOUBLEMATRIX *pLogPWM, char chStrand,
					int *pMask);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_Genome_GetCutoff_Main: compute cutoff to calibrate */
/*  motif occurence rate.                                                  */
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanMatrix_Genome_GetCutoff_Main(char strMotifListPath[], char strGenomePath[], 
					char strCodPath[], char strOutputPath[],
					int nBGOrder, char strBGType[], char strBGPath[], int nBGStepSize,
					int nUseCS, double dC, char strCSPath[],
					int nIncludeRepeat);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_Genome_Summary_Main: compute motif enrichment      */
/*  for targeted genomic regions. All motifs in the motif list will be     */
/*  examined.                                                              */
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanMatrix_Genome_Summary_Main(char strMotifListPath[], char strGenomePath[], 
					char strCodPath[], char strNegCodPath[], char strOutputPath[],
					int nBGOrder, char strBGType[], char strBGPath[], int nBGStepSize,
					int nUseCS, double dC, char strCSPath[],
					int nIncludeRepeat);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanConsensus_Genome_Summary_Main: compute motif enrichment   */
/*  for targeted genomic regions. All motifs in the motif list will be     */
/*  examined.                                                              */
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanConsensus_Genome_Summary_Main(char strMotifListPath[], char strGenomePath[], 
					char strCodPath[], char strNegCodPath[], char strOutputPath[],
					int nUseCS, double dC, char strCSPath[],
					int nIncludeRepeat);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_Genome_Enrich_Main: compute motif enrichment       */
/*  for targeted and ordered genomic regions. Target regions will be       */
/*  grouped into several tiers according to nTierSize. The enrichment      */
/*  level of each tier will be compared to negative control regions.       */
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanMatrix_Genome_Enrich_Main(char strMotifPath[], char strGenomePath[], 
					char strCodPath[], char strNegCodPath[], int nTierSize, 
					char strOutputPath[], double dR, 
					int nBGOrder, char strBGType[], char strBGPath[], int nBGStepSize,
					int nUseCS, double dC, char strCSPath[],
					int nIncludeRepeat);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanConsensus_Genome_Enrich_Main: compute motif enrichment    */
/*  for targeted and ordered genomic regions. Target regions will be       */
/*  grouped into several tiers according to nTierSize. The enrichment      */
/*  level of each tier will be compared to negative control regions.       */
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanConsensus_Genome_Enrich_Main(char strMotifPath[], char strGenomePath[], 
					char strCodPath[], char strNegCodPath[], int nTierSize, 
					char strOutputPath[], int nMC, int nMD, 
					int nUseCS, double dC, char strCSPath[],
					int nIncludeRepeat);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_GetSiteAround_Main: get regions surrounding sites derived     */
/*  from motif mapping                                                     */
/* ----------------------------------------------------------------------- */ 
int MotifMap_GetSiteAround_Main(char strInputPath[], char strGenomePath[], 
			char strOutputPath[], char strSpecies[], int nW, int nCN, int nA);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_getcluster_Main: get clustered sites.                         */
/*  If distance between 2 sites <= nW, they will be report as clustered.   */
/*  nInputType is the input file type, 0=cod, 1=bed, 2=codp, 3=bedp.       */
/* ----------------------------------------------------------------------- */ 
int MotifMap_getcluster_Main(char strInputPath[], char strOutputPath[],
							 int nW, int nInputType);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_GetSiteAroundCS_Genome_Main: get mean conservation score      */
/*  TFBS sites.                                                            */
/* ----------------------------------------------------------------------- */ 
int MotifMap_GetSiteAroundCS_Genome_Main(char strInputPath[], char strOutputPath[], 
		int nMotifLen, int nW, char strSpecies[], char strGenomePath[], char strCSPath[]);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_Filter_Genome_Main: fiter TFBS by footprint conservation      */
/*  and mask TFBS in coding sequences.                                     */
/* ----------------------------------------------------------------------- */ 
int MotifMap_Filter_Genome_Main(char strInputPath[], char strOutputPath[], 
		int nUseCS, double dC, char strMaskPath[], char strCSPath[], 
		int nUseCds, double dCds, char strCdsPath[]);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_SimuScoreDistribution_typeI: get motif mapping score          */
/*  distribution through simulation. This function will tell you if a      */
/*  sequence looks more like a motif than the background or not.           */ 
/* ----------------------------------------------------------------------- */ 
int MotifMap_SimuScoreDistribution_typeI(char strMotifName[], struct DOUBLEMATRIX *pMotifPWM, 
								   int nBGOrder, struct DOUBLEMATRIX *pBG0, 
								   struct DOUBLEMATRIX *pBGF, struct DOUBLEMATRIX *pBGR, 
								   int nSimuNum, struct DOUBLEMATRIX *pQ, 
								   char strCutoffPath[]);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_SimuScore_typeI: simulate a background sequence and calculate */
/*  its motif mapping score.                                               */
/* ----------------------------------------------------------------------- */ 
double MotifMap_SimuScore_typeI(struct BYTEMATRIX *pSeq, int nBurninRound,
						  struct DOUBLEMATRIX *pLogPWM, int nBGOrder, 
						  struct DOUBLEMATRIX *pBG0, 
						  struct DOUBLEMATRIX *pBGF, struct DOUBLEMATRIX *pBGR,
						  struct DOUBLEMATRIX *pLogBGF, struct DOUBLEMATRIX *pLogBGR);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_SimuScoreDistribution_typeII: get motif score distribution    */
/*  through simulation. Motifs will be simulated from PWM and quantiles    */
/*  for log(Scores) will be estimated. This function will give you a sense */
/*  of to what extent a sequence looks like a motif instead of telling     */
/*  if that sequence looks more like a motif than like the background.     */  
/* ----------------------------------------------------------------------- */ 
int MotifMap_SimuScoreDistribution_typeII(char strMotifName[], 
				struct DOUBLEMATRIX *pMotifPWM, int nSimuNum, 
				struct DOUBLEMATRIX *pQ, char strCutoffPath[]);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_SimuScore_typeII: simulate a background sequence and          */
/*  calculate its motif mapping score.                                     */
/* ----------------------------------------------------------------------- */ 
double MotifMap_SimuScore_typeII(struct BYTEMATRIX *pSeq, struct DOUBLEMATRIX *pPWM, 
						  struct DOUBLEMATRIX *pLogPWM);


/* ----------------------------------------------------------------------- */ 
/*  MotifMap_SimuScoreDistribution_typeIII: get motif mapping score        */
/*  distribution through simulation. This function will tell you the distn */
/*  of likelihood ratio under the motif model.                             */ 
/* ----------------------------------------------------------------------- */ 
int MotifMap_SimuScoreDistribution_typeIII(char strMotifName[], struct DOUBLEMATRIX *pMotifPWM, 
								   int nBGOrder, struct DOUBLEMATRIX *pBG0, 
								   struct DOUBLEMATRIX *pBGF, struct DOUBLEMATRIX *pBGR, 
								   int nSimuNum, struct DOUBLEMATRIX *pQ, 
								   char strCutoffPath[]);


/* ----------------------------------------------------------------------- */ 
/*  MotifMap_SimuScore_typeIII: simulate a motif sequence and calculate    */
/*  its motif mapping score.                                               */
/* ----------------------------------------------------------------------- */ 
double MotifMap_SimuScore_typeIII(struct BYTEMATRIX *pSeq, int nBurninRound, 
						  struct DOUBLEMATRIX *pLogPWM, int nBGOrder, 
						  struct DOUBLEMATRIX *pPWM, struct DOUBLEMATRIX *pBG0,
						  struct DOUBLEMATRIX *pBGF, struct DOUBLEMATRIX *pBGR,
						  struct DOUBLEMATRIX *pLogBGF, struct DOUBLEMATRIX *pLogBGR);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScoreRefGene: Motif mapping based on already mapped sites     */
/*  This function will score refgene based on known genomic motif mapping. */ 
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScoreRefGene(char strSpecies[], char strVersion[], 
		char strGenomePath[], char strCScorePath[], 
		char strTargetMotifList[], char strMotifMapPath[], double dCSCutoff, int nRepeatMask,
		char strRefGeneFile[], int nTSSUP, int nTSSDOWN, int nTESUP, int nTESDOWN, 
		int nIncludeIntron, int nIncludeExon,
		char strOutPath[]);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_GetContigIndex: get contig index for a position.              */ 
/* ----------------------------------------------------------------------- */ 
int MotifMap_GetContigIndex(int nChrId, int nStart, 
							struct INTMATRIX *pChrContigNum, int nContigSize);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_GetStatistics: read summary statistics from motifmap.         */ 
/* ----------------------------------------------------------------------- */ 
int MotifMap_GetStatistics(char strFileName[], double *pEffecLen, double *pSiteNum);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMultipleMatrix_Genome_Main: map multiple motif matrix to  */
/*  genomic regions, using nBGOrder markov chain as background.            */
/*  strBGType specifies whether region-derived or genomic location         */
/*  specific background should be used. After mapping, each pair of motifs */
/*  will be tested to see if they are clustered together.                  */ 
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanMultipleMatrix_Genome_Main(char strSpecies[], char strGenomePath[], 
				char strWorkPath[], char strInputPath[], char strOutputPath[],  
				int nBGOrder, char strBGType[], char strBGPath[], int nBGStepSize,
				int nUseCS, char strCSPath[], int nIncludeRepeat, int nW);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMultipleMatrix_ClusterTest: test clustering of motif site */ 
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanMultipleMatrix_ClusterTest(int nMotifNum, struct tagString **vMotifName, 
			char strSortedSitePath[], int nTotSiteNum, int nW, struct DOUBLEMATRIX *pNeighborNum);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMultipleMatrix_Export: save summary statistics.           */ 
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanMultipleMatrix_Export(char strFileName[], int nMotifNum, 
		struct tagString **vMotifName, int nPosOnly,
		struct DOUBLEMATRIX *pEffecLenPos, struct DOUBLEMATRIX *pSiteNumPos, 
		struct DOUBLEMATRIX *pNeighborNumPos, struct DOUBLEMATRIX *pNeighborLenPos, 
		struct DOUBLEMATRIX *pNeighborSiteNumPos, struct DOUBLEMATRIX *pEffecLenNeg, 
		struct DOUBLEMATRIX *pSiteNumNeg, struct DOUBLEMATRIX *pNeighborNumNeg);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_CountKmer_Sequential_Main: count k-mers in a set of           */
/*  FASTA sequences.                                                       */
/*  dC  = the cutoff for average conservation scores.                      */
/* ----------------------------------------------------------------------- */ 
int MotifMap_CountKmer_Sequential_Main(char strInputPath[], char strSeqFile[],  
					char strOutputPath[], int nK, int nUseCS, 
					double dC, char strCSPrefix[]);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_CountKmer_In_FlexSeqMtf: count kmer                           */
/* ----------------------------------------------------------------------- */ 
int MotifMap_CountKmer_In_FlexSeqMtf(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,
				int nUseCS, double dC, int nK, struct DOUBLEMATRIX *pKmerCount);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_Rescore_Main: rescore the TFBS using the background model     */
/*  learnt from current sites.                                             */
/* ----------------------------------------------------------------------- */ 
int MotifMap_Rescore_Main(char strMotifPath[], char strInputPath[], char strOutputPath[]);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_FilterOverlappingSite:                                        */
/*  Filter overlapping sites. Return the number of non-overlapping sites.  */
/* ----------------------------------------------------------------------- */ 
int MotifMap_FilterOverlappingSite(char strInFile[], char strOutFile[]);

/* ----------------------------------------------------------------------- */ 
/*  MotifSampler_Quality: Motif sampler based on two stage A/B quality     */
/*  control.                                                               */
/* ----------------------------------------------------------------------- */ 
int MotifSampler_Quality(char strInFile[], char strOutFile[], int nIteration, 
						   int nBGOrder, struct DOUBLEMATRIX *pAFreqPrior,
						   int nMotifNum, 
						   struct DOUBLEMATRIX *pMatrixPrior[],
						   int nSeedNum,
						   struct DOUBLEMATRIX *pSeedMotif[]);

/* ----------------------------------------------------------------------- */ 
/*  Quality_Sampler: Quality sampler for one sequence.                     */
/* ----------------------------------------------------------------------- */ 
int Quality_Sampler(struct SEQMOTIF *pSeqMtf, double vFreqCount[MAX_TYPE_NUM][2],
					int nBGOrder, struct DOUBLEMATRIX *pLogBG, struct DOUBLEMATRIX *pLogBG0,
					int nMotifNum, int vMotifLen[], int nMaxMotifLen,
					struct MOTIFMATRIX *vMotif[],
					int niter, int nIteration);

/* ----------------------------------------------------------------------- */ 
/*  Quality_EstimateTruePosProb: Get the true positive probability.        */
/* ----------------------------------------------------------------------- */ 
double Quality_EstimateTruePosProb(int nNeighborCount, int niter, int nIteration);

/* ----------------------------------------------------------------------- */ 
/*  Quality_Neighbor_Likelihood: Get the likelihood for neighbors.         */
/* ----------------------------------------------------------------------- */ 
int Quality_Neighbor_Likelihood(struct SITEPOS *pNeighbors, unsigned char *pCount,
								unsigned char *pB, double *pLn0, double *pLn1,
								int niter, int nIteration);

/* ----------------------------------------------------------------------- */ 
/*  Quality_BaseLogLikelihood_Background: get loglikelihood for background */
/* ----------------------------------------------------------------------- */ 
double Quality_BaseLogLikelihood_Background(unsigned char *pSite, int nMaxMotifLen, 
											struct DOUBLEMATRIX *pLogBG, int nBGOrder, 
											struct DOUBLEMATRIX *pLogBG0, int nPreNum, 
											char chStrand, int *pMask);

/* ----------------------------------------------------------------------- */ 
/*  Quality_BaseLogLikelihood_Motif: get loglikelihood for motif.          */
/* ----------------------------------------------------------------------- */ 
double Quality_BaseLogLikelihood_Motif(unsigned char *pSite, int nMaxMotifLen, 
									   int nMotifLen, struct DOUBLEMATRIX *pLogPWMp, 
									   struct DOUBLEMATRIX *pLogBG, int nBGOrder, 
									   struct DOUBLEMATRIX *pLogBG0, 
									   char chStrand, int *pMask);

/* ----------------------------------------------------------------------- */ 
/*  Quality_CallMotifSites: Call motif sites.                              */
/* ----------------------------------------------------------------------- */ 
int Quality_CallMotifSites(struct SEQMOTIF *pSeqMtf, int vMotifLen[]);

/* ----------------------------------------------------------------------- */ 
/*  QUALITY_RandInit: Random init.                                         */
/* ----------------------------------------------------------------------- */ 
int Quality_RandInit(struct SEQMOTIF *pSeqMtf, int nBGOrder, 
					 struct DOUBLEMATRIX *pBG, 
					 int nMotifNum, struct MOTIFMATRIX *vMotif[], 
					 int vMotifLen[], int nMaxMotifLen, int nSeedNum, 
					 struct DOUBLEMATRIX *pSeedMotif[], 
					 double vFreqCount[MAX_TYPE_NUM][2]);

/* ----------------------------------------------------------------------- */ 
/*  QUALITY_UpdateNeighborCount: Update neighboring count                  */
/* ----------------------------------------------------------------------- */ 
int Quality_UpdateNeighborCount(struct SEQMOTIF *pSeqMtf);


/* ----------------------------------------------------------------------- */ 
/*  MotifSampler_Collapsed: Motif sampler based on collapsed gibbs sampler */
/*  Find one motif per call                                                */
/* ----------------------------------------------------------------------- */ 
int MotifSampler_Collapsed(char strInFile[], char strOutFile[], int nIteration, 
						   int nBGOrder, struct DOUBLEMATRIX *pFreqPrior,
						   struct DOUBLEMATRIX *pMatrixPrior);

/* ----------------------------------------------------------------------- */ 
/*  AddNucleicMatrixCount: Add count to motif matrix, for nucleic acid.    */
/* ----------------------------------------------------------------------- */ 
int AddNucleicMatrixCount(struct DOUBLEMATRIX *pCountMat, char chStrand, 
						  unsigned char pSite[], int nMotifLen);

/* ----------------------------------------------------------------------- */ 
/*  BaseLikelihood_Motif: get collapsed likelihood for motif.              */
/* ----------------------------------------------------------------------- */ 
double BaseLikelihood_Motif(unsigned char pSite[], int nMotifLen, struct DOUBLEMATRIX *pCountMat, char chStrand);

/* ----------------------------------------------------------------------- */ 
/*  BaseLikelihood_Background: get likelihood for background.              */
/* ----------------------------------------------------------------------- */ 
double BaseLikelihood_Background(unsigned char pSite[], int nMotifLen, struct DOUBLEMATRIX *pBG, int nBGOrder, int nPreNum, char chStrand);

/* ----------------------------------------------------------------------- */ 
/*  SubtractNucleicMatrixCount: Add count to motif matrix, nucleic acid.   */
/* ----------------------------------------------------------------------- */ 
int SubtractNucleicMatrixCount(struct DOUBLEMATRIX *pCountMat, char chStrand, 
						  unsigned char pSite[], int nMotifLen);

/* ----------------------------------------------------------------------- */ 
/*  MotifSampler_HMM: Motif sampler based on HMM.                          */
/*  This sampler incorporates module, tiling, conservation information     */
/* ----------------------------------------------------------------------- */ 
int MotifSampler_HMM(char strInFile[], char strOutFile[], int nIteration, 
						   int nBGOrder, int nBGTypeNum,
						   struct DOUBLEMATRIX *pModuleTransP, 
						   struct DOUBLEMATRIX *pFreqPrior[],
						   int nMotifNum, 
						   struct DOUBLEMATRIX *pMatrixPrior[],
						   int nSeedNum,
						   struct DOUBLEMATRIX *pSeedMotif[]);

/* ----------------------------------------------------------------------- */ 
/*  HMM_RandInit: Randomly initialize the motifs.                          */
/* ----------------------------------------------------------------------- */ 
int HMM_RandInit(struct SEQMOTIF *pSeqMtf, 
				 int nBGOrder, struct DOUBLEMATRIX *pBG, 
				 int nMotifNum, struct MOTIFMATRIX *vMotif[], int vMotifLen[],
				 int nSeedNum, struct DOUBLEMATRIX *pSeedMotif[], 
				 double vFreqPrior[]);

/* ----------------------------------------------------------------------- */ 
/*  HMM_CallMotifSites: Call motif sites.                                  */
/* ----------------------------------------------------------------------- */ 
int HMM_CallMotifSites(struct SEQMOTIF *pSeqMtf,
					   int nBGOrder, struct DOUBLEMATRIX *pBG,
					   int nMotifNum, struct MOTIFMATRIX *vMotif[], int vMotifLen[],
					   double vFreqPrior[2][MAX_HMMSTATUS_NUM]);

/* ----------------------------------------------------------------------- */ 
/*  HMM_ForwardSummation: Forward integration procedure of HMM.            */
/* ----------------------------------------------------------------------- */ 
int HMM_ForwardSummation(struct SEQMOTIF *pSeqMtf, int nHMMStatusNum,
						 struct DOUBLEMATRIX *pModuleTransP,
						 double vFreqPrior[2][MAX_HMMSTATUS_NUM],
						 int nBGOrder, struct DOUBLEMATRIX *pBG,
						 int nMotifNum, int vMotifLen[], 
						 struct MOTIFMATRIX *vMotif[]);

/* ----------------------------------------------------------------------- */ 
/*  HMM_BackwardSampling: Backward sampling procedure of HMM.              */
/* ----------------------------------------------------------------------- */ 
int HMM_BackwardSampling(struct SEQMOTIF *pSeqMtf, int nHMMStatusNum,
						 struct DOUBLEMATRIX *pModuleTransP, 
						 double vFreqCount[2][MAX_HMMSTATUS_NUM],
						 int nMotifNum, int vMotifLen[], 
						 struct MOTIFMATRIX *vMotif[],
						 int nRecordModuleOn);

/* ----------------------------------------------------------------------- */ 
/*  WriteMotifToFile: Report motif discovery results.                      */
/*  Find one motif per call                                                */
/* ----------------------------------------------------------------------- */ 
int WriteMotifToFile(char strOutFile[], 
					 int nMotifNum, struct MOTIFMATRIX *vMotif[],
					 struct DOUBLEMATRIX *pBG0,
					 int nCount, struct SEQMOTIF *vSeqMtf[]);

/* ----------------------------------------------------------------------- */ 
/*  BaseLogLikelihood_MotifPWM: get loglikelihood for motif.               */
/* ----------------------------------------------------------------------- */ 
double BaseLogLikelihood_MotifPWM(unsigned char pSite[], int nMotifLen, struct DOUBLEMATRIX *pLogPWMMat, char chStrand);

/* ----------------------------------------------------------------------- */ 
/*  BaseLikelihood_MotifPWM: get likelihood for motif.                     */
/* ----------------------------------------------------------------------- */ 
double BaseLikelihood_MotifPWM(unsigned char pSite[], int nMotifLen, struct DOUBLEMATRIX *pPWMMat, char chStrand);

/* ----------------------------------------------------------------------- */ 
/*                     MTFMCREATE: create PWM of motifs                    */
/* ----------------------------------------------------------------------- */ 
struct MOTIFMATRIX *MTFMCREATE();

/* ----------------------------------------------------------------------- */ 
/*                    MTFMDESTROY: destroy PWM of motifs                   */
/* ----------------------------------------------------------------------- */ 
int MTFMDESTROY(struct MOTIFMATRIX *pM);

/* ----------------------------------------------------------------------- */ 
/*                  MTFMDESTROYLIST: destroy PWM list of motifs            */
/* ----------------------------------------------------------------------- */ 
int MTFMDESTROYLIST(struct MOTIFMATRIX **pMList);

/* ----------------------------------------------------------------------- */ 
/*                  MTFMSAMPLEPWM: sample PWM from counts                  */
/* ----------------------------------------------------------------------- */ 
int MTFMSAMPLEPWM(struct MOTIFMATRIX *pM);


/* ----------------------------------------------------------------------- */ 
/*                  MTFMSCORE: get motif scores                            */
/* ----------------------------------------------------------------------- */ 
double MTFMSCORE(struct MOTIFMATRIX *pM, struct DOUBLEMATRIX *pBG);

/* ----------------------------------------------------------------------- */ 
/*           MTFMREFRESH_QUALITY: set PWMp = log(COUNTp/sum(COUNTp)        */
/* ----------------------------------------------------------------------- */ 
int MTFMREFRESH_QUALITY(struct MOTIFMATRIX *pM);

/* ----------------------------------------------------------------------- */ 
/*                     MTFSCREATE: create site of motifs                   */
/* ----------------------------------------------------------------------- */ 
struct MOTIFSITE *MTFSCREATE();

/* ----------------------------------------------------------------------- */ 
/*                    MTFSDESTROY: destroy motif site                      */
/* ----------------------------------------------------------------------- */ 
int MTFSDESTROY(struct MOTIFSITE *pS);

/* ----------------------------------------------------------------------- */ 
/*                  MTFSDESTROYLIST: destroy list of motif sites           */
/* ----------------------------------------------------------------------- */ 
int MTFSDESTROYLIST(struct MOTIFSITE **pSList);

/* ----------------------------------------------------------------------- */ 
/*                  SEQMTFCREATE: create seq/motif complex                 */
/* ----------------------------------------------------------------------- */ 
struct SEQMOTIF *SEQMTFCREATE(int nId, int nSeqNum, int nScoreNum, int nPathNum);

/* ----------------------------------------------------------------------- */ 
/*                  SEQMTFDESTROY: destroy seq/motif complex               */
/* ----------------------------------------------------------------------- */ 
int SEQMTFDESTROY(struct SEQMOTIF *pSeqMtf);

/* ----------------------------------------------------------------------- */ 
/*              SEQMTFWRITESEQTOFASTA: write sequence to fasta file        */
/* ----------------------------------------------------------------------- */ 
int SEQMTFWRITESEQTOFASTA(FILE *fpOut, struct SEQMOTIF *pSeqMtf);

/* ----------------------------------------------------------------------- */ 
/*          SEQMTFWRITESAMPPATHTOFASTA: write sample path to fasta file    */
/* ----------------------------------------------------------------------- */ 
int SEQMTFWRITESAMPPATHTOFASTA(FILE *fpOut, struct SEQMOTIF *pSeqMtf);

/* ----------------------------------------------------------------------- */ 
/*          SEQMTFCREATESEQ: allocate memory required for a sequence       */
/* ----------------------------------------------------------------------- */ 
int SEQMTFCREATESEQ(struct SEQMOTIF *pSeqMtf, int nIndex, int nLen);

/* ----------------------------------------------------------------------- */ 
/*          SEQMTFDESTROYSEQ: release memory required for a sequence       */
/* ----------------------------------------------------------------------- */ 
int SEQMTFDESTROYSEQ(struct SEQMOTIF *pSeqMtf, int nIndex);

/* ----------------------------------------------------------------------- */ 
/*          SEQMTFCREATESCORE: allocate memory required for a score        */
/* ----------------------------------------------------------------------- */ 
int SEQMTFCREATESCORE(struct SEQMOTIF *pSeqMtf, int nIndex, int nLen);

/* ----------------------------------------------------------------------- */ 
/*          SEQMTFDESTROYSCORE: release memory required for a score        */
/* ----------------------------------------------------------------------- */ 
int SEQMTFDESTROYSCORE(struct SEQMOTIF *pSeqMtf, int nIndex);

/* ----------------------------------------------------------------------- */ 
/*          SEQMTFCREATEPATH: allocate memory required for a sample path   */
/* ----------------------------------------------------------------------- */ 
int SEQMTFCREATEPATH(struct SEQMOTIF *pSeqMtf, int nIndex, int nLen);

/* ----------------------------------------------------------------------- */ 
/*          SEQMTFDESTROYPATH: release memory required for a sample path   */
/* ----------------------------------------------------------------------- */ 
int SEQMTFDESTROYPATH(struct SEQMOTIF *pSeqMtf, int nIndex);

/* ----------------------------------------------------------------------- */ 
/*          SEQMTFCREATEHMM: allocate memory required for HMM.             */
/* ----------------------------------------------------------------------- */ 
int SEQMTFCREATEHMM(struct SEQMOTIF *pSeqMtf, int nStatusNum, int nLen);

/* ----------------------------------------------------------------------- */ 
/*          SEQMTFDESTROYHMM: release memory required for HMM.             */
/* ----------------------------------------------------------------------- */ 
int SEQMTFDESTROYHMM(struct SEQMOTIF *pSeqMtf);

/* ----------------------------------------------------------------------- */ 
/*          SEQMTFCODENUCLEICSEQ: code nucleic acid sequences.             */
/* ----------------------------------------------------------------------- */ 
int SEQMTFCODENUCLEICSEQ(struct SEQMOTIF *pSeqMtf, int nIndex, int nLen, char strSeq[]);

/* ----------------------------------------------------------------------- */ 
/*       SEQMTFLOADPRIMARYSEQFROMFASTA: load sequence from fasta file.     */
/* Return a vector of pointers to newly created motif/seq complexes, and   */
/* *pCount, the number of loaded complexes.                                */
/* ----------------------------------------------------------------------- */
struct SEQMOTIF **SEQMTFLOADPRIMARYSEQFROMFASTA(int *pCount, char strFilePath[]);

/* ----------------------------------------------------------------------- */ 
/* SEQMTFLOADPRIMARYSEQFROMFASTA_FORQUALITY: load sequence from fasta file */
/* Return a vector of pointers to newly created motif/seq complexes, and   */
/* nCount, the number of loaded complexes.                                 */
/* ----------------------------------------------------------------------- */
struct SEQMOTIF **SEQMTFLOADPRIMARYSEQFROMFASTA_FORQUALITY(int *pCount, char strFilePath[]);

/* ----------------------------------------------------------------------- */ 
/* SEQMTFLOADPRIMARYSEQFROMFASTA_FORREGRESSOR: load sequence from fasta    */
/* file and corresponding conservation scores.                             */
/* Return a vector of pointers to newly created motif/seq complexes, and   */
/* nCount, the number of loaded complexes.                                 */
/* ----------------------------------------------------------------------- */
struct SEQMOTIF **SEQMTFLOADPRIMARYSEQFROMFASTA_FORREGRESSOR(int *pCount, char strFilePath[], char strFastaName[]);

/* ----------------------------------------------------------------------- */ 
/*       SEQMTFESTIMATENUCLEICBGMC: estimate markov background.            */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *SEQMTFESTIMATENUCLEICBGMC(struct SEQMOTIF **pSeqMtf, int nCount, int nIndex, int nOrder);

/* ----------------------------------------------------------------------- */ 
/*       SEQMTFESTIMATENUCLEICBGMC: estimate markov background.            */
/* ----------------------------------------------------------------------- */ 
double SEQMTF_MOTIFMAP_SCORE(struct SEQMOTIF *pSeqMtf, int nIndex, 
							 int nBGOrder, struct DOUBLEMATRIX *pBG,
							 struct DOUBLEMATRIX *pMotif,
							 double dRatioCutoff);

/* ----------------------------------------------------------------------- */ 
/*      SEQMTFADDCOUNTTONUCLEICBGMC: add count to markov background.       */
/* ----------------------------------------------------------------------- */ 
int SEQMTFADDCOUNTTONUCLEICBGMC(struct DOUBLEMATRIX *pBG, struct SEQMOTIF *pSeqMtf, int nIndex, int nOrder);

/* ----------------------------------------------------------------------- */ 
/*      SEQMTF_LOADSEQFROMGENOME: load sequences from genome.              */
/* This function will load sequences from fpSeq (from nStart base to nEnd  */
/* base), and write it to the nSeqId th sequence of pSeqMtf, starting from */
/* nOffset th bp.                                                          */
/* ----------------------------------------------------------------------- */ 
int SEQMTF_LOADSEQFROMGENOME(struct SEQMOTIF *pSeqMtf, int nSeqId, int nOffset,
							 FILE *fpSeq, int nStart, int nEnd, char chStrand);

/* ----------------------------------------------------------------------- */ 
/* SEQMTF_ADDCOUNTTONUCLEICBGMC_BOTHSTRAND: add count to markov background */
/* both strands will be considered and two backgroud matrix will be        */
/* generated, one for forward background, one for backward backgroud       */
/* Forward Background: P( b[j+3] | b[j], b[j+1], b[j+2] )                  */
/* Backward Background: P( b[j-1] | b[j], b[j+1], b[j+2] )                 */
/* the sequences nIndex: nFromPos-nToPos will be considered.               */
/* ----------------------------------------------------------------------- */ 
int SEQMTF_ADDCOUNTTONUCLEICBGMC_BOTHSTRAND(struct DOUBLEMATRIX *pBGF, struct DOUBLEMATRIX *pBGR, 
											struct SEQMOTIF *pSeqMtf, int nIndex, int nEffectiveLen,
											int nFromPos, int nToPos, int nOrder);

/* ----------------------------------------------------------------------- */ 
/* SEQMTF_SETBGMCSCORE_BOTHSTRAND: get background log-likelihood for       */
/* every base in the sequences.                                            */
/* ----------------------------------------------------------------------- */ 
int SEQMTF_SETBGMCSCORE_BOTHSTRAND(struct SEQMOTIF *pSeqMtf, int nIndex, int nEffectiveLen, 
								   int nFromPos, int nToPos, 
								   struct DOUBLEMATRIX *pBGF, struct DOUBLEMATRIX *pBGR, 
								   int nOrder);

/* ----------------------------------------------------------------------- */ 
/* SEQMTF_MAPMOTIF_BOTHSTRANDASBG: map a motif to a genomic sequence.      */
/* the map locations will be recorded in motif site list.                  */
/* ----------------------------------------------------------------------- */ 
int SEQMTF_MAPMOTIF_BOTHSTRANDASBG(struct DOUBLEMATRIX *pMotifPWM, int nMotifId, double dCutoff, 
								   struct SEQMOTIF *pSeqMtf, int nIndex, int nEffectiveLen, 
								   int nFromPos, int nToPos);

/* ----------------------------------------------------------------------- */ 
/* SEQMTF_WRITEMOTIFSITETOBED_THENCLEARSITES: write motifs to bed files.   */
/* ----------------------------------------------------------------------- */ 
int SEQMTF_WRITEMOTIFSITETOBED_THENCLEARSITES(struct SEQMOTIF *pSeqMtf, int nMaxMotifNum, 
											  struct tagString **vMotifName, 
											  char strChrName[], int nOffset,
											  FILE **vfpMapOut);

/* ----------------------------------------------------------------------- */ 
/* SEQMTF_LOADSEQCSFROMGENOME_FORMOTIFFILTER: load sequence and            */
/* conservation score.                                                     */
/* ----------------------------------------------------------------------- */ 
int SEQMTF_LOADSEQCSFROMGENOME_FORMOTIFFILTER(struct SEQMOTIF *pSeqMtf, int nSeqId, 
											  char strChrName[], int nStart, int nEnd, 
											  char strGenomePath[], char strCScorePath[]);

/* ----------------------------------------------------------------------- */ 
/*             MCMOVECREATE: create markov chain move elements             */
/* ----------------------------------------------------------------------- */ 
struct MARKOVCHAINMOVE *MCMOVECREATE();

/* ----------------------------------------------------------------------- */ 
/*             MCMOVEDESTROY: destroy markov chain move elements           */
/* ----------------------------------------------------------------------- */ 
int MCMOVEDESTROY(struct MARKOVCHAINMOVE *pM);

/* ----------------------------------------------------------------------- */ 
/*             MCMOVEDESTROYCIRCLE: destroy markov chain circle            */
/* ----------------------------------------------------------------------- */ 
int MCMOVEDESTROYCIRCLE(struct MARKOVCHAINMOVE **pMList);

/* ----------------------------------------------------------------------- */ 
/*             SITEPOSCREATE: create site pos elements                     */
/* ----------------------------------------------------------------------- */ 
struct SITEPOS *SITEPOSCREATE();

/* ----------------------------------------------------------------------- */ 
/*             SITEPOSDESTROY: destroy site pos elements                   */
/* ----------------------------------------------------------------------- */ 
int SITEPOSDESTROY(struct SITEPOS *pS);

/* ----------------------------------------------------------------------- */ 
/*             SITEPOSDESTROYLIST: destroy position list                   */
/* ----------------------------------------------------------------------- */ 
int SITEPOSDESTROYLIST(struct SITEPOS **pSList);

/* ----------------------------------------------------------------------- */ 
/*             SITEPOSLISTGETCOUNT: get # of site pos elements             */
/* ----------------------------------------------------------------------- */ 
int SITEPOSLISTGETCOUNT(struct SITEPOS *pSList);

/* ----------------------------------------------------------------------- */ 
/*              MTFSGRPCREATE: create motif site group                     */
/* ----------------------------------------------------------------------- */ 
struct MOTIFSITEGROUP *MTFSGRPCREATE(int nMotifNum);

/* ----------------------------------------------------------------------- */ 
/*             MTFSGRPDESTROY: destroy motif site group                    */
/* ----------------------------------------------------------------------- */ 
int MTFSGRPDESTROY(struct MOTIFSITEGROUP *pG);


/* ----------------------------------------------------------------------- */ 
/*             FLEXMTFSCREATE: create site of flexmotifs                   */
/* ----------------------------------------------------------------------- */ 
struct FLEXMOTIFSITE *FLEXMTFSCREATE();

/* ----------------------------------------------------------------------- */ 
/*             FLEXMTFSDESTROY: destroy flexmotif site                     */
/* ----------------------------------------------------------------------- */ 
int FLEXMTFSDESTROY(struct FLEXMOTIFSITE *pS);

/* ----------------------------------------------------------------------- */ 
/*             FLEXMTFSDESTROYLIST: destroy list of flex motif sites       */
/* ----------------------------------------------------------------------- */ 
int FLEXMTFSDESTROYLIST(struct FLEXMOTIFSITE **pSList);

/* ----------------------------------------------------------------------- */ 
/* FLEXMTFSADDTOLIST_SORTBYPROB: add flex motif sites to list, and sort    */
/* them according to the posterior probability.                            */
/* ----------------------------------------------------------------------- */ 
int FLEXMTFSADDTOLIST_SORTBYPROB(struct FLEXMOTIFSITE **pSList, struct FLEXMOTIFSITE *pSite);

/* ----------------------------------------------------------------------- */ 
/* FLEXMTFSADDTOLIST_SORTBYPOS: add flex motif sites to list, and sort     */
/* them according to their position.                                       */
/* ----------------------------------------------------------------------- */ 
int FLEXMTFSADDTOLIST_SORTBYPOS(struct FLEXMOTIFSITE **pSList, struct FLEXMOTIFSITE *pSite);

/* ----------------------------------------------------------------------- */ 
/*             FLEXMTFSREMOVESITE: remove one site from the matrix         */
/* ----------------------------------------------------------------------- */ 
int FLEXMTFSREMOVESITE(struct FLEXMOTIFSITE *pS, int nPos);

/* ----------------------------------------------------------------------- */ 
/*            FLEXMTFSADDSITE: add one site to the context matrix.         */
/* ----------------------------------------------------------------------- */ 
int FLEXMTFSADDSITE(struct FLEXMOTIFSITE *pS, int nPos, int nMotifId, int nStrand);

/* ----------------------------------------------------------------------- */ 
/*             FLEXSEQMTFCREATE: create flex seq/motif complex             */
/* ----------------------------------------------------------------------- */ 
struct FLEXSEQMOTIF *FLEXSEQMTFCREATE(int nId, int nSeqNum, int nScoreNum, 
									  int nStatusNum, int nMonitorNum);

/* ----------------------------------------------------------------------- */ 
/*             FLEXSEQMTFDESTROY: destroy flex seq/motif complex           */
/* ----------------------------------------------------------------------- */ 
int FLEXSEQMTFDESTROY(struct FLEXSEQMOTIF *pSeqMtf);

/* ----------------------------------------------------------------------- */ 
/*      FLEXSEQMTFCREATESEQ: allocate memory required for a sequence       */
/* ----------------------------------------------------------------------- */ 
int FLEXSEQMTFCREATESEQ(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, int nLen);

/* ----------------------------------------------------------------------- */ 
/*      FLEXSEQMTFDESTROYSEQ: release memory required for a sequence       */
/* ----------------------------------------------------------------------- */ 
int FLEXSEQMTFDESTROYSEQ(struct FLEXSEQMOTIF *pSeqMtf, int nIndex);

/* ----------------------------------------------------------------------- */ 
/*      FLEXSEQMTFCREATESCORE: allocate memory required for a score        */
/* ----------------------------------------------------------------------- */ 
int FLEXSEQMTFCREATESCORE(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, int nLen);

/* ----------------------------------------------------------------------- */ 
/*      FLEXSEQMTFDESTROYSCORE: release memory required for a score        */
/* ----------------------------------------------------------------------- */ 
int FLEXSEQMTFDESTROYSCORE(struct FLEXSEQMOTIF *pSeqMtf, int nIndex);

/* ----------------------------------------------------------------------- */ 
/*      FLEXSEQMTFCREATESTATUS: allocate memory required for a status      */
/* ----------------------------------------------------------------------- */ 
int FLEXSEQMTFCREATESTATUS(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, int nLen);

/* ----------------------------------------------------------------------- */ 
/*      FLEXSEQMTFDESTROYSTATUS: release memory required for a status      */
/* ----------------------------------------------------------------------- */ 
int FLEXSEQMTFDESTROYSTATUS(struct FLEXSEQMOTIF *pSeqMtf, int nIndex);

/* ----------------------------------------------------------------------- */ 
/*      FLEXSEQMTFCREATEMONITOR: allocate memory required for a monitor    */
/* ----------------------------------------------------------------------- */ 
int FLEXSEQMTFCREATEMONITOR(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, int nLen);

/* ----------------------------------------------------------------------- */ 
/*      FLEXSEQMTFDESTROYMONITOR: release memory required for a monitor    */
/* ----------------------------------------------------------------------- */ 
int FLEXSEQMTFDESTROYMONITOR(struct FLEXSEQMOTIF *pSeqMtf, int nIndex);

/* ----------------------------------------------------------------------- */ 
/*      FLEXSEQMTFCREATEMOTIF: allocate memory required for motifs         */
/* ----------------------------------------------------------------------- */ 
int FLEXSEQMTFCREATEMOTIF(struct FLEXSEQMOTIF *pSeqMtf, int nLen);

/* ----------------------------------------------------------------------- */ 
/*      FLEXSEQMTFDESTROYMOTIF: release memory required for motifs         */
/* ----------------------------------------------------------------------- */ 
int FLEXSEQMTFDESTROYMOTIF(struct FLEXSEQMOTIF *pSeqMtf);

/* ----------------------------------------------------------------------- */ 
/* FLEXSEQMTFLOADSEQ: load sequence and scores from input files            */
/* Return a vector of pointers to newly created flexmotif/seq complexes,   */
/* and nSeqCount, the number of loaded complexes.                          */
/* ----------------------------------------------------------------------- */
struct FLEXSEQMOTIF **FLEXSEQMTFLOADSEQ(int *nSeqCount, char strWorkPath[], 
										char strSeqFile[], 
										int nUseCS, char strCSPrefix[],
										int nMotifNum);

/* ----------------------------------------------------------------------- */ 
/* FLEXSEQMTFLOADSEQ_FOR_BGFIT: load sequence and scores from input files  */
/* This loading is for fitting background likelihood.                      */
/* Return a vector of pointers to newly created flexmotif/seq complexes,   */
/* and nSeqCount, the number of loaded complexes.                          */
/* ----------------------------------------------------------------------- */
struct FLEXSEQMOTIF **FLEXSEQMTFLOADSEQ_FOR_BGFIT(int *nSeqCount, char strWorkPath[], 
										char strSeqFile[],
										int nUseCS, char strCSPrefix[],
										int nMotifNum);

/* ----------------------------------------------------------------------- */ 
/* FLEXSEQMTFLOADSEQ_FOR_SCANCONSENSUS_GENOME: load single sequence and    */
/* scores from genome. Return a pointer to newly created                   */
/* flexmotif/seq complexes.                                                */
/* ----------------------------------------------------------------------- */
struct FLEXSEQMOTIF *FLEXSEQMTFLOADSEQ_FOR_SCANCONSENSUS_GENOME(int nSeqId,
			char strGenomePath[], char strChr[], int nStart, int nEnd, 
			int nUseCS, char strCSPath[]);

/* ----------------------------------------------------------------------- */ 
/* FLEXSEQMTFLOADSEQ_FOR_SCANCONSENSUS_SINGLE: load single sequence and    */
/* scores from input sequences. Return a pointer to newly created          */
/* flexmotif/seq complexes and nSeqCount, the number of loaded complexes.  */
/* ----------------------------------------------------------------------- */
struct FLEXSEQMOTIF *FLEXSEQMTFLOADSEQ_FOR_SCANCONSENSUS_SINGLE(struct tagSequence *pSeq,
					char strWorkPath[], int nUseCS, char strCSPrefix[]);

/* ----------------------------------------------------------------------- */ 
/* FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_SINGLE: load single sequence and       */
/* scores from input sequences. Return a pointer to newly created          */
/* flexmotif/seq complexes.                                                */
/* ----------------------------------------------------------------------- */
struct FLEXSEQMOTIF *FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_SINGLE(struct tagSequence *pSeq,
					char strWorkPath[], int nUseCS, char strCSPrefix[]);

/* ----------------------------------------------------------------------- */ 
/* FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_GENOME: load single sequence and       */
/* scores from input sequences. Return a pointer to newly created          */
/* flexmotif/seq complexes.                                                */
/* ----------------------------------------------------------------------- */
struct FLEXSEQMOTIF *FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_GENOME(int nSeqId,
				char strGenomePath[], char strChr[], int nStart, int nEnd,
				int nUseCS, char strCSPath[]);

/* ----------------------------------------------------------------------- */ 
/* FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_SINGLE: load single sequence and       */
/* scores from input sequences. Return a pointer to newly created          */
/* flexmotif/seq complexes.                                                */
/* ----------------------------------------------------------------------- */
struct FLEXSEQMOTIF *FLEXSEQMTFLOADSEQ_FOR_COUNTMCBACKGROUND_GENOME(int nSeqId,
				char strGenomePath[], char strChr[], int nStart, int nEnd,
				int nUseCS, char strCSPath[]);

/* ----------------------------------------------------------------------- */ 
/* FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_GROUP: load sequence&scores from input */
/* files. Return a vector of pointers to newly created flexmotif/seq       */
/* complexes and nSeqCount, the number of loaded complexes.                */
/* ----------------------------------------------------------------------- */
struct FLEXSEQMOTIF **FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_GROUP(int *nSeqCount, 
				char strWorkPath[], char strSeqFile[], int nUseCS, char strCSPrefix[]);

/* ----------------------------------------------------------------------- */ 
/* FLEXSEQMTFLOADSEQ_FOR_HASHGENOME: load single sequence and              */
/* scores from input sequences. Return a pointer to newly created          */
/* flexmotif/seq complexes.                                                */
/* ----------------------------------------------------------------------- */
struct FLEXSEQMOTIF *FLEXSEQMTFLOADSEQ_FOR_HASHGENOME(int nSeqId,
				char strGenomePath[], char strChr[], int nStart, int nEnd,
				int nUseCS, char strCSPath[]);

/* ----------------------------------------------------------------------- */ 
/*      FLEXSEQMTFCODENUCLEICSEQ: code nucleic acid sequences.             */
/* ----------------------------------------------------------------------- */ 
int FLEXSEQMTFCODENUCLEICSEQ(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, int nLen, char strSeq[]);

/* ----------------------------------------------------------------------- */ 
/*      FLEXSEQMTFLOADCS: load conservation score.                         */
/* ----------------------------------------------------------------------- */ 
int FLEXSEQMTFLOADCS(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, int nSeqLen, char strConsFile[]);

/* ----------------------------------------------------------------------- */ 
/*    FLEXSEQMTFESTIMATENUCLEICBGMC: estimate markov background.           */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *FLEXSEQMTFESTIMATENUCLEICBGMC(struct FLEXSEQMOTIF **vSeqMtf, int nCount, 
								int nIndex, int nOrder, int nBaseTypeNum);

/* ----------------------------------------------------------------------- */ 
/*      FLEXMTFMCREATE: create PWM of flex motifs                          */
/* ----------------------------------------------------------------------- */ 
struct FLEXMOTIFMATRIX *FLEXMTFMCREATE();

/* ----------------------------------------------------------------------- */ 
/*      FLEXMTFMDESTROY: destroy PWM of flex motifs                        */
/* ----------------------------------------------------------------------- */ 
int FLEXMTFMDESTROY(struct FLEXMOTIFMATRIX *pM);

/* ----------------------------------------------------------------------- */ 
/*  FLEXMTFMREFRESH: update PWM based on current count matrix.             */
/* ----------------------------------------------------------------------- */ 
int FLEXMTFMREFRESH(struct FLEXMOTIFMATRIX *pM);

/* ----------------------------------------------------------------------- */ 
/*  FLEXMTFMSCORE: get motif score.                                        */
/* ----------------------------------------------------------------------- */ 
double FLEXMTFMSCORE(struct FLEXMOTIFMATRIX *pM, struct DOUBLEMATRIX *pBG);

/* ----------------------------------------------------------------------- */ 
/*  FLEXMTFMSHIFTPRIORCOUNT: shift prior count.                            */
/* ----------------------------------------------------------------------- */
int FLEXMTFMSHIFTPRIORCOUNT(struct FLEXMOTIFMATRIX *pM, int nMotifLen, 
							int nOffset, double dDefaultPrior);


/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_Integrate_Genome_Main: map motif probability       */
/*  matrix to genomic regions, and then integrate signal strength and      */
/*  export the strength at specific positions specified by strSampPath.    */
/*  Using nBGOrder markov chain as background and dR as                    */
/*  likelihood ratio cutoff. dC is the conservation score cutoff.          */
/*  strBGType specifies whether region-derived or genomic location         */
/*  specific background should be used.                                    */ 
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanMatrix_Integrate_Genome_Main(char strMotifPath[], 
					char strGenomePath[], char strCodPath[], 
					char strSampPath[], char strOutputPath[], 
					char strSpecies[], int nW, int nRepeat, double dR, 
					int nBGOrder, char strBGType[], 
					char strBGPath[], int nBGStepSize,
					int nUseCS, double dC, char strCSPath[]);

/* ----------------------------------------------------------------------- */ 
/* FLEXSEQMTFLOADSEQ_FOR_INTEGRATEMATRIX_GENOME: load single sequence and  */
/* scores from input sequences. Return a pointer to newly created          */
/* flexmotif/seq complexes.                                                */
/* ----------------------------------------------------------------------- */
struct FLEXSEQMOTIF *FLEXSEQMTFLOADSEQ_FOR_INTEGRATEMATRIX_GENOME(int nSeqId,
				char strGenomePath[], char strChr[], int nStart, int nEnd,
				int nUseCS, char strCSPath[]);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_IntegrateMatrix_In_FlexSeqMtf: integrate motif signals.       */
/* ----------------------------------------------------------------------- */ 
int MotifMap_IntegrateMatrix_In_FlexSeqMtf(struct FLEXSEQMOTIF *pSeqMtf, 
			int nIndex, int nW, 
			struct DOUBLEMATRIX *pLogPWM, double dR, 
			int nUseCS, double dC, 
			int *pEffecLen, int *pConsEffecLen, int *pSiteNum);

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_IntegrateMatrix_Export_Genome: export integrated singal.      */
/* ----------------------------------------------------------------------- */ 
int MotifMap_IntegrateMatrix_Export_Genome(struct FLEXSEQMOTIF *pSeqMtf, 
				int nIndex, int nChr, int nActualStart, int nFrom, int nTo, 
				int *nCurrentChr, int *nCurrentPos, FILE *fpSamp, FILE *fpOut,
				int *pnNewSampLine);

/* ----------------------------------------------------------------------- */ 
/*  FLEXSEQMTF_REPEATUNMASK: unmask repeat sequences.                      */
/* ----------------------------------------------------------------------- */ 
int FLEXSEQMTF_REPEATUNMASK(struct FLEXSEQMOTIF *pSeqMtf, int nIndex);

/* ----------------------------------------------------------------------- */
/* tred2 Function                                                          */
/* Householder reduction of a real, symmetric matrix a[1:n][1:n].          */
/* On output, a is replaced by the orthogonal matrix Q effecting the       */
/* transformation. d[1:n] returns the diagonal elements of the             */
/* resulting tridiagonal matrix, and e[1:n] the off-diagonal elements,     */
/* with e[1] = 0.                                                          */
/* ----------------------------------------------------------------------- */
void tred2(double **a, int n, double d[], double e[]);


