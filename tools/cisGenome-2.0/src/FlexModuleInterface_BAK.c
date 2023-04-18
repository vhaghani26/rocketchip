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
	int ni,nj,nk;
	struct INTMATRIX *pAPos;
	struct INTMATRIX *pOK;
	struct DOUBLEMATRIX *pMScore;
	int nAnchorNum;
	double dDist,dMinDist;
	double dModLen = 200.0;

	/* int nZeroType,nSum; */

	/* compute score */
	/* dScore = nContextNum+nAddone-1.0;
	if(dScore < 0.0)
		dScore = 0.0; */
	
	/* compute score */
	/* ------------------------------------ */
	/* anchor distance                      */
	/* ------------------------------------ */
	pOK = NULL;
	pOK = CreateIntMatrix(1, nMotifNum);
	if(pOK == NULL)
	{
		printf("Error: FlexModule_ComputeModuleScore, cannot allocate memory for computing module score!\n");
		exit(EXIT_FAILURE);
	}

	nAnchorNum = 0;
	for(ni=0; ni<nContextNum; ni++)
	{
		nk = pMotifId[ni]-1;
		pOK->pMatElement[nk] += 1;
	}

	if(nAddone == 1)
	{
		nk = nAddMotifId-1;
		pOK->pMatElement[nk] += 1;
	}

	nAnchorNum = pOK->pMatElement[0];

	if(nAnchorNum == 0)
	{
		dScore = 0.0;
		DestroyIntMatrix(pOK);
		return dScore;
	}
	
	dScore = 0.0;
	for(ni=0; ni<nMotifNum; ni++)
	{
		if(pOK->pMatElement[ni] > 0)
		{
			dScore = dScore+1.0+0.5*(pOK->pMatElement[ni]-1.0);
		}
	}
	
	DestroyIntMatrix(pOK);
	
	/* pMScore = NULL;
	pMScore = CreateDoubleMatrix(1, nMotifNum);
	if(pMScore == NULL)
	{
		printf("Error: FlexModule_ComputeModuleScore, cannot allocate memory for computing module score!\n");
		exit(EXIT_FAILURE);
	}

	pAPos = NULL;
	pAPos = CreateIntMatrix(1, nAnchorNum);
	if(pAPos == NULL)
	{
		printf("Error: FlexModule_ComputeModuleScore, cannot create anchor pos!\n");
		exit(EXIT_FAILURE);
	}

	nj = 0;
	for(ni=0; ni<nContextNum; ni++)
	{
		if(pMotifId[ni] == 1)
		{
			pAPos->pMatElement[nj] = pPos[ni];
			nj++;
		}
	}

	if(nAddone == 1)
	{
		if(nAddMotifId == 1)
		{
			pAPos->pMatElement[nj] = nAddPos;
			nj++;
		}
	}

	for(ni=0; ni<nContextNum; ni++)
	{
		if(pMotifId[ni] != 1)
		{
			nk = pMotifId[ni]-1;

			dMinDist = dModLen;
			for(nj=0; nj<nAnchorNum; nj++)
			{
				dDist = fabs(pPos[ni]-pAPos->pMatElement[nj]);
				if(dMinDist > dDist)
				{
					dMinDist = dDist;
				}
			}
			
			dDist = 1.0-dMinDist/dModLen;
			if(dDist < 0.0)
				dDist = 0.0;

			pMScore->pMatElement[nk] += dDist;
		}
	}

	if(nAddone == 1)
	{
		if(nAddMotifId != 1)
		{
			nk = nAddMotifId-1;

			dMinDist = dModLen;
			for(nj=0; nj<nAnchorNum; nj++)
			{
				dDist = fabs(nAddPos-pAPos->pMatElement[nj]);
				if(dMinDist > dDist)
				{
					dMinDist = dDist;
				}
			}
			
			dDist = 1.0-dMinDist/dModLen;
			if(dDist < 0.0)
				dDist = 0.0;

			pMScore->pMatElement[nk] += dDist;
		}
	}

	dScore = 0.0;
	for(ni=0; ni<nMotifNum; ni++)
	{
		dScore += pMScore->pMatElement[ni];
	}
	dScore += (double)nAnchorNum;

	DestroyIntMatrix(pAPos);
	DestroyIntMatrix(pOK);
	DestroyDoubleMatrix(pMScore);
	*/
	
	
	/* ------------------------------------ */
	/* motif mode                           */
	/* ------------------------------------ */
	/* dScore = 1000.0; */
	
	/* ------------------------------------ */
	/* each motif has at least one site     */
	/* >1 motif sites is equivalent to 1    */
	/* site                                 */
	/* ------------------------------------ */
	/* pOK = NULL;
	pOK = CreateIntMatrix(1, nMotifNum);
	if(pOK == NULL)
	{
		printf("Error: FlexModule_ComputeModuleScore, cannot allocate memory for computing module score!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nContextNum; ni++)
	{
		nk = pMotifId[ni]-1;
		pOK->pMatElement[nk] += 1;
	}

	if(nAddone == 1)
	{
		pOK->pMatElement[nAddMotifId-1] += 1;
	}
	
	dScore = 0.0;
	for(ni=0; ni<nMotifNum; ni++)
	{
		if(pOK->pMatElement[ni] > 0)
			dScore += 1.0;
	}

	if(dScore < ((double)nMotifNum-0.5))
	{
		dScore = 0.0;
	}
	else
	{
		dScore = 0.0;
		for(ni=0; ni<nMotifNum; ni++)
		{
			dScore += pOK->pMatElement[ni];
		}
		if(dScore < 0.0)
			dScore = 0.0;
	}

	DestroyIntMatrix(pOK);
	*/

	
	
	/* nZeroType = 0;
	nOK1 = 0;
	nOK2 = 0;
	nOK3 = 0;
	for(ni=0; ni<nContextNum; ni++)
	{
		if(pPos[ni] == 0)
		{
			nZeroType = pMotifId[ni];
		}

		if(pMotifId[ni] == 1)
		{
			nOK1 += 1;
		}
		else if(pMotifId[ni] == 2)
		{
			nOK2 += 1;
		}
		else if(pMotifId[ni] == 3)
		{
			nOK3 += 1;
		}
	}

	if(nAddone == 1)
	{
		if(nAddPos == 0)
		{
			nZeroType = nAddMotifId;
		}

		if(nAddMotifId == 1)
		{
			nOK1 += 1;
		}
		else if(nAddMotifId == 2)
		{
			nOK2 += 1;
		}
		else if(nAddMotifId == 3)
		{
			nOK3 += 1;
		}
	}

	if(nZeroType == 1)
		dScore = nOK1-1.0;
	else if(nZeroType == 2)
		dScore = nOK2-1.0;
	else if(nZeroType == 3)
		dScore = nOK3-1.0;
	else
		dScore = 0.0;

	if(dScore < 0.0)
		dScore = 0.0;
	*/

	/* dScore = 0.0;
	for(ni=0; ni<nContextNum; ni++)
	{
		if( (pPos[ni] == 0) && (pStrand[ni] == 1) )
		{
			return 0.0;
		}

		if(pStrand[ni] == 0)
		{
			dScore += 1.0;
		}
	}

	if(nAddone == 1)
	{
		if( (nAddPos == 0) && (nAddStrand == 1) )
		{
			return 0.0;
		}

		if(nAddStrand == 0)
		{
			dScore += 1.0;
		}
	}

	dScore -= 1.0;
	if(dScore < 0.0)
		dScore = 0.0;

	*/

	/* nOK1 = 0;
	nOK2 = 0;
	nOK3 = 0;
	for(ni=0; ni<nContextNum; ni++)
	{
		if(pMotifId[ni] == 1)
		{
			nOK1 += 1;
		}
		else if(pMotifId[ni] == 2)
		{
			nOK2 += 1;
		}
		else if(pMotifId[ni] == 3)
		{
			nOK3 += 1;
		}
	}

	if(nAddone == 1)
	{
		if(nAddMotifId == 1)
		{
			nOK1 += 1;
		}
		else if(nAddMotifId == 2)
		{
			nOK2 += 1;
		}
		else if(nAddMotifId == 3)
		{
			nOK3 += 1;
		}
	}

	dScore = nOK1+nOK2+nOK3-1.0;

	if(nOK1 > 0)
		nOK1 = 1;
	if(nOK2 > 0)
		nOK2 = 1;
	if(nOK3 > 0)
		nOK3 = 1;
	if( (nOK1+nOK2+nOK3) < 2)
	{
		dScore = 0.0;
	}

	if(dScore < 0.0)
		dScore = 0.0;
	*/


	/* ------------------------------------ */
	/* each motif has one and only one site */
	/* ------------------------------------ */
	/* nOK1 = 0;
	nOK2 = 0;
	nOK3 = 0;
	nZeroType = 0;
	for(ni=0; ni<nContextNum; ni++)
	{
		if(pPos[ni] == 0)
		{
			if(pMotifId[ni] == 1)
			{
				nZeroType = 1;
			}
			else if(pMotifId[ni] == 2)
			{
				nZeroType = 2;
			}
			else if(pMotifId[ni] == 3)
			{
				nZeroType = 3;
			}

		}
		else
		{
			if(pMotifId[ni] == 1)
			{
				nOK1 += 1;
			}
			else if(pMotifId[ni] == 2)
			{
				nOK2 += 1;
			}
			else if(pMotifId[ni] == 3)
			{
				nOK3 += 1;
			}
		}
	}

	if(nAddone == 1)
	{
		if(nAddPos == 0)
		{
			if(nAddMotifId == 1)
			{
				nZeroType = 1;
			}
			else if(nAddMotifId == 2)
			{
				nZeroType = 2;
			}
			else if(nAddMotifId == 3)
			{
				nZeroType = 3;
			}
		}
		else
		{
			if(nAddMotifId == 1)
			{
				nOK1 += 1;
			}
			else if(nAddMotifId == 2)
			{
				nOK2 += 1;
			}
			else if(nAddMotifId == 3)
			{
				nOK3 += 1;
			}
		}
	}

	if(nZeroType == 0)
	{
		if(nOK1 > 0)
			nOK1 = 1;
		if(nOK2 > 0)
			nOK2 = 1;
		if(nOK3 > 0)
			nOK3 = 1;
		nSum = nOK1+nOK2+nOK3;
		if(nSum == 0)
			dScore = 0.0;
		else if(nSum == 1)
			dScore = 1.0;
		else if(nSum == 2)
			dScore = 2.0;
		else if(nSum == 3)
			dScore = 0.0;
	}
	else if(nZeroType == 1)
	{
		if(nOK2 > 0)
			nOK2 = 1;
		if(nOK3 > 0)
			nOK3 = 1;
		nSum = nOK2+nOK3;

		if(nOK1 > 0)
			dScore = 0.0;
		else if(nSum == 0)
			dScore = 0.0;
		else if(nSum == 1)
			dScore = 1.0;
		else if(nSum == 2)
			dScore = 2.0;
	}
	else if(nZeroType == 2)
	{
		if(nOK1 > 0)
			nOK1 = 1;
		if(nOK3 > 0)
			nOK3 = 1;
		nSum = nOK1+nOK3;

		if(nOK2 > 0)
			dScore = 0.0;
		else if(nSum == 0)
			dScore = 0.0;
		else if(nSum == 1)
			dScore = 1.0;
		else if(nSum == 2)
			dScore = 2.0;
	}
	else if(nZeroType == 3)
	{
		if(nOK1 > 0)
			nOK1 = 1;
		if(nOK2 > 0)
			nOK2 = 1;
		nSum = nOK1+nOK2;

		if(nOK3 > 0)
			dScore = 0.0;
		else if(nSum == 0)
			dScore = 0.0;
		else if(nSum == 1)
			dScore = 1.0;
		else if(nSum == 2)
			dScore = 2.0;
	}


	if(dScore < 0.0)
		dScore = 0.0;
	*/

	/* if(nOK1 <= 0)
		dScore /= 2.0; */

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
	/* dA = 1.0-dT;
	dR = dA+(1.0-dA)*dContextScore/dModuleD;
	if(dR >= 0.95)
		dR = 0.95;
	if(dR < 0.04)
		dR = 0.04; */
	
	dA = 1.0-(1.0-dModuleL)*dT;
	dX = dModuleA*(dContextScore*2.0/dModuleD-1.0);
	dX = exp(dX);
	dR = dA+(1.0-dA)*dX/(1.0+dX);
	if(dR >= 0.999)
		dR = 0.999;
	
	/* double d1,d2,d3;
	d1 = 0.04;
	d2 = 0.4;
	d3 = 0.95;
	
	/* get pdf */
	/*if(dContextScore <= 1.1)
	{
		dR = 0.999-(1.0-d1)*dT;
		if(dR < d1)
			dR = d1;
	}
	else if(dContextScore <= 2.1)
	{
		dR = 0.999-(1.0-d2)*dT;
		if(dR < d2)
			dR = d2;
	}
	else if(dContextScore <= 3.1)
	{
		dR = 0.999-(1.0-d3)*dT;
		if(dR < d3)
			dR = d3;
	}
	else if(dContextScore >= 3.9)
		dR = 0.999; */

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
	/* dNewD = *pModuleD + normrnd(0.0, 0.2);
	if( (dNewD < dInitModuleD) || (dNewD >= (2.0*dInitModuleD)) )
		dNewD = *pModuleD;
	*/
	
	/* dNewL = rand_u(); */
	/* dNewL = *pModuleL + 0.05*(rand_u()-0.5);
	if( (dNewL<0.001) || (dNewL >= 0.5) )
		dNewL = *pModuleL;
	dNewA = *pModuleA + 2.0*(rand_u()-0.5);
	if( dNewA<1e-6 )
		dNewA = *pModuleA; */

	dNewD = *pModuleD + normrnd(0.0, 0.2);
	if( (dNewD < dInitModuleD) || (dNewD >= (2.0*dInitModuleD)) )
		dNewD = *pModuleD;
	
	dNewL = *pModuleL + 0.0005*(rand_u()-0.5);
	if( (dNewL<0.001) || (dNewL >= 0.01) )
		dNewL = *pModuleL;
	dNewA = *pModuleA + 2.0*(rand_u()-0.5);
	if( dNewA<1e-6 )
		dNewA = *pModuleA;


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