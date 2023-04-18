/* ----------------------------------------------------------------------- */
/*  FlexModuleInterface.h : interface of the flexmodule library            */
/*  Author : Ji HongKai ; Time: 2005.07                                    */
/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------- */
/*                              Declaration                                */
/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------- */ 
/*                                Functions                                */
/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_ComputeModuleScore: compute context score.                  */
/* ----------------------------------------------------------------------- */ 
double FlexModule_ComputeModuleScore(int nContextNum, int *pPos,
			int *pMotifId, int *pStrand,
			int nAddone, int nAddPos, int nAddMotifId, int nAddStrand,
			int nMotifNum);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_ComputeTrueProb: compute probability for a motif to be      */
/*  true.                                                                  */
/* ----------------------------------------------------------------------- */ 
double FlexModule_ComputeTrueProb(double dContextScore, double dModuleD, 
								  double dModuleL, double dModuleA, double dT);

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_SampleD: sample module parameter D.                         */
/* ----------------------------------------------------------------------- */ 
int FlexModule_SampleD(struct FLEXSEQMOTIF **vSeqMtf, int nSeqCount, 
				int nIndex, double dT, double dInitModuleD,
				double *pModuleD, double *pModuleL, double *pModuleA, int nModuleLen, 
				int nMotifNum, struct FLEXMOTIFMATRIX **vMotif, 
				struct INTMATRIX *pMotifLen, int nMaxMotifLen, 
				int nIsRecording);