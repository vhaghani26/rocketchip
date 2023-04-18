/* ----------------------------------------------------------------------- */
/*  FlexModuleInterface.c : user interface of flexmodule functions         */
/*  Author : Ji HongKai ; Time: 2005.07                                    */
/* ----------------------------------------------------------------------- */

#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
#include "limits.h"
#include "time.h"

#include "MathLib.h"
#include "MatrixLib.h"
#include "StringLib.h"
#include "RandomLib.h"
#include "SequenceLib.h"
#include "MotifLib.h"
#include "GenomeLib.h"
#include "AffyLib.h"
#include "FlexModuleInterface.h"

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_ComputeModuleScore: compute context score.                  */
/* ----------------------------------------------------------------------- */ 
double FlexModule_ComputeModuleScore(int nContextNum, int *pPos,
			int *pMotifId, int *pStrand,
			int nAddone, int nAddPos, int nAddMotifId, int nAddStrand,
			int nMotifNum)
{
	/* define */
	double dScore = 0.0;
	
	/* compute score */
	dScore = nContextNum+nAddone-1.0;
	if(dScore < 0.0)
		dScore = 0.0;
	
	/* return */
	return dScore;
}

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_ComputeTrueProb: compute probability for a motif to be      */
/*  true.                                                                  */
/* ----------------------------------------------------------------------- */ 
double FlexModule_ComputeTrueProb(double dContextScore, double dModuleD, 
								  double dModuleL, double dModuleA, double dT)
{
	/* define */
	double dR = 1.0;
	double dA;
	double dX;

	/* compute probability */
	dA = 1.0-(1.0-dModuleL)*dT;
	dX = dModuleA*(dContextScore*2.0/dModuleD-1.0);
	dX = exp(dX);
	dR = dA+(1.0-dA)*dX/(1.0+dX);
	if(dR >= 0.999)
		dR = 0.999;

	/* return */
	return dR;
}

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_SampleD: sample module parameter D.                         */
/* ----------------------------------------------------------------------- */ 
int FlexModule_SampleD(struct FLEXSEQMOTIF **vSeqMtf, int nSeqCount, 
				int nIndex, double dT, double dInitModuleD,
				double *pModuleD, double *pModuleL, double *pModuleA, int nModuleLen, 
				int nMotifNum, struct FLEXMOTIFMATRIX **vMotif, 
				struct INTMATRIX *pMotifLen, int nMaxMotifLen, 
				int nIsRecording)
{
	/* define */
	double dNewD,dNewL,dNewA;
	double dLike0,dLike1,dR,dRand;
	int ni,nj;
	unsigned char *pAi;
	struct FLEXMOTIFSITE *pSite;
	double dScore,dPAi0,dPAi1;
	double dMu0 = 1000000.0;

	/* get trial */
	dNewD = *pModuleD + normrnd(0.0, 0.2);
	if( (dNewD < dInitModuleD) || (dNewD >= (2.0*dInitModuleD)) )
		dNewD = *pModuleD;
	
	/* dNewL = rand_u(); */
	dNewL = *pModuleL + 0.05*(rand_u()-0.5);
	if( (dNewL<0.001) || (dNewL >= 0.99) )
		dNewL = *pModuleL;
	dNewA = *pModuleA + 2.0*(rand_u()-0.5);
	if( dNewA<1e-6 )
		dNewA = *pModuleA;
	/* dNewA = *pModuleA; */

	/* get likelihood */
	dLike0 = 0.0;
	dLike1 = 0.0;
	
	for(ni=0; ni<nSeqCount; ni++)
	{
		pAi = vSeqMtf[ni]->vStatus[0]->pMatElement;
		pSite = vSeqMtf[ni]->pMotifList;
		while(pSite != NULL)
		{
			dScore = FlexModule_ComputeModuleScore(pSite->pContextPos->nWidth, pSite->pContextPos->pMatElement,
				pSite->pContextMotif->pMatElement, pSite->pContextStrand->pMatElement, 0, 0, 0, 0,
				nMotifNum);
			
			dPAi0 = FlexModule_ComputeTrueProb(dScore, *pModuleD, *pModuleL, *pModuleA, dT);
			dPAi1 = FlexModule_ComputeTrueProb(dScore, dNewD, dNewL, dNewA, dT);

			nj = pSite->nStartPos;
			if(pAi[nj] == 1)
			{
				dLike0 += log(dPAi0);
				dLike1 += log(dPAi1);
			}
			else
			{
				dLike0 += log(1.0-dPAi0);
				dLike1 += log(1.0-dPAi1);
			}

			/* get next */
			pSite = pSite->pNext;
		}
	}

	dLike0 -= ((*pModuleD-dInitModuleD)/dMu0+(*pModuleA)/dMu0);
	dLike1 -= ((dNewD-dInitModuleD)/dMu0+dNewA/dMu0);

	dR = exp(dLike1-dLike0);

	dRand = rand_u();
	if(dRand < dR)
	{
		*pModuleD = dNewD;
		*pModuleL = dNewL;
		*pModuleA = dNewA;
	}

	/* return */
	return PROC_SUCCESS;
}