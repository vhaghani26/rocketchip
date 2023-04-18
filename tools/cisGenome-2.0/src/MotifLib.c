/* ----------------------------------------------------------------------- */
/*  MotifLib.c : implementation of the motif library                       */
/*  Author : Ji HongKai ; Time: 2004.07                                    */
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
#include "GenomeLib.h"
#include "MotifLib.h"
#include "AffyLib.h"
#include "FlexModuleInterface.h"


/* ----------------------------------------------------------------------- */ 
/*  FlexModule_Main: flexmodule main function                              */
/*  sample cis-regulatory module with conservation scores.                 */
/* ----------------------------------------------------------------------- */ 
int FlexModule_Main(char strParamFile[])
{
	/* define */

	/* for loading parameters */
	FILE *fpIn;
	char strLine[MED_LINE_LENGTH];
	char *chSep,*chSep2;
	int nlen;
	int nError = 0;
	int ni;

	/* parameters */
	char strParamPath[MED_LINE_LENGTH];
	char strWorkPath[MED_LINE_LENGTH];
	char strSeqFile[LINE_LENGTH];
	char strOutFile[LINE_LENGTH];
	int nMotifNum = 1;
	int nMeanMotifLen = 10;
	int nMotifLenUpperBound = 255;
	struct INTMATRIX *pMotifLen = NULL;
	struct tagString **vPriorPWMFile = NULL;
	double dInitModuleD = 1.0;
	int nModuleLen = 200;
	int nSampleModuleLen = 0;
	int nMCNum = 500;
	int nUseCS = 0;
	int nBGOrder = 3;
	int nUseFittedBG = 0;
	char strFittedBGPrefix[LINE_LENGTH];
	char strCSPrefix[LINE_LENGTH];
	char strCSLike[LINE_LENGTH];
	char strPriorAbundance[LINE_LENGTH];

	/* init */
	strcpy(strWorkPath, "");
	
	/* --------------- */
	/* load parameters */
	/* --------------- */
	fpIn = NULL;
	fpIn = fopen(strParamFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: FlexModule_Main, cannot open the parameter file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strstr(strLine, "[Working Directory]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n%s\n", strLine, chSep);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strWorkPath, chSep);
			
			AdjustDirectoryPath(strWorkPath);
		}


		else if(strstr(strLine, "[FASTA Sequence]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n%s\n", strLine, chSep);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep == '\0')
			{
				printf("Error: FlexModule_Main, there are no sequences!\n");
				nError = 1;
				break;
			}
			else
			{
				strcpy(strSeqFile, chSep);
			}
		}

		else if(strstr(strLine, "[Output File]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n%s\n", strLine, chSep);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep == '\0')
			{
				printf("Error: FlexModule_Main, you have to specify a file for exporting the results!\n");
				nError = 1;
				break;
			}
			else
			{
				strcpy(strOutFile, chSep);
			}
		}
		
		else if(strstr(strLine, "[Motif Number K]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n%s\n", strLine, chSep);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nMotifNum = atoi(chSep);
			}
		}

		else if(strstr(strLine, "[Mean Motif Length Lamda]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n%s\n", strLine, chSep);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nMeanMotifLen = atoi(chSep);
			}
		}

		else if(strstr(strLine, "[Maximal Motif Length Allowed]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n%s\n", strLine, chSep);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nMotifLenUpperBound = atoi(chSep);
			}
		}

		else if(strstr(strLine, "[Init Motif Length L]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n%s\n", strLine, chSep);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep == '\0')
			{
				printf("Error: FlexModule_Main, please specify the motif lengths!\n");
				nError = 1;
				break;
			}
			else
			{
				StrTrimLeft(chSep);

				pMotifLen = NULL;
				pMotifLen = CreateIntMatrix(1,nMotifNum);
				if(pMotifLen == NULL)
				{
					printf("Error: FlexModule_Main, cannot load initial motif lengths!\n");
					exit(EXIT_FAILURE);
				}

				ni = 0;
				chSep2 = strchr(chSep, ',');
				while(chSep2 != NULL)
				{
					*chSep2 = '\0';
					pMotifLen->pMatElement[ni] = atoi(chSep);
					chSep = chSep2+1;
					StrTrimLeft(chSep);
					chSep2 = strchr(chSep, ',');
					ni++;
				}

				pMotifLen->pMatElement[ni] = atoi(chSep);
				ni++;

				if(ni != nMotifNum)
				{
					printf("Error: FlexModule_Main, motif number not match!\n");
					exit(EXIT_FAILURE);
				}
			}
		}


		else if(strstr(strLine, "[Init Motif Matrix]") == strLine)
		{
			vPriorPWMFile = NULL;
			vPriorPWMFile = (struct tagString **)calloc(nMotifNum, sizeof(struct tagString *));
			if(vPriorPWMFile == NULL)
			{
				printf("Error: FlexModule_Main, cannot allocate memory for storing prior PWM path!\n");
				exit(EXIT_FAILURE);
			}

			for(ni=0; ni<nMotifNum; ni++)
			{
				fgets(strLine, MED_LINE_LENGTH, fpIn);
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				StringAddTail(vPriorPWMFile+ni, strLine);
			}

		}

		else if(strstr(strLine, "[Init Motif Abundance]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n%s\n", strLine, chSep);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep == '\0')
			{
				strcpy(strPriorAbundance, "NULL");
			}
			else
			{
				strcpy(strPriorAbundance, chSep);
			}
		}
		
		else if(strstr(strLine, "[Init Module Size D]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n%s\n", strLine, chSep);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				dInitModuleD = atof(chSep);
			}
		}

		else if(strstr(strLine, "[Module Length]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n%s\n", strLine, chSep);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nModuleLen = atoi(chSep);
			}
		}

		else if(strstr(strLine, "[Sample Module Length?]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n%s\n", strLine, chSep);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nSampleModuleLen = atoi(chSep);
			}
		}

		

		else if(strstr(strLine, "[Order of Background Markov Chain]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n%s\n", strLine, chSep);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nBGOrder = atoi(chSep);
			}
		}

		else if(strstr(strLine, "[Use Fitted Background?]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n%s\n", strLine, chSep);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nUseFittedBG = atoi(chSep);
			}
		}

		else if(strstr(strLine, "[Fitted Background Prefix]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n%s\n", strLine, chSep);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strFittedBGPrefix, chSep);
		}

		else if(strstr(strLine, "[MCMC Iteration]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n%s\n", strLine, chSep);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nMCNum = atoi(chSep);
			}
		}

		else if(strstr(strLine, "[Use *.cs?]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n%s\n", strLine, chSep);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nUseCS = atoi(chSep);
			}
		}

		else if(strstr(strLine, "[CS Prefix]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n%s\n", strLine, chSep);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strCSPrefix, chSep);
		}
		
		else if(strstr(strLine, "[CS Likelihood f]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n%s\n", strLine, chSep);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep == '\0')
			{
				strcpy(strCSLike, "NULL");
			}
			else
			{
				strcpy(strCSLike, chSep);
			}
		}

		else if(strLine[0] == '#')
		{
		}

		else
		{
			printf("Error: FlexModule_Main, unknown parameters!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	fclose(fpIn);

	if(nError > 0)
	{
		printf("Error: FlexModule_Main, cannot load the parameter file correctly!\n");
		exit(EXIT_FAILURE);
	}


	/* flexmodule */
	FlexModule_Sampler(strWorkPath, strSeqFile, strOutFile,
		/* motif number K, motif length poisson parameter, motif len distribution bin number */
		nMotifNum, nMeanMotifLen, nMotifLenUpperBound,
		/* init motif lengths, init motif pseudocounts (PWMs) */
		pMotifLen, vPriorPWMFile, 
		/* init module score parameter d, module length */
		dInitModuleD, nModuleLen, nSampleModuleLen, 
		/* order of background markov chain */
		nBGOrder, nUseFittedBG, strFittedBGPrefix,
		/* MC draw number, using conservation indicator, conservation file prefix, conservation likelihood file */
		nMCNum, nUseCS, strCSPrefix, strCSLike, 
		/* prior abundance pseudocounts */
		strPriorAbundance);

	/* release memory */
	DestroyIntMatrix(pMotifLen);
	for(ni=0; ni<nMotifNum; ni++)
	{
		DeleteString(vPriorPWMFile[ni]);
		vPriorPWMFile[ni] = NULL;
	}
	free(vPriorPWMFile);

	/* #################################### */
	/* Prepare CisGenome Ini files          */
	/* #################################### */
	strcpy(strParamPath, strParamFile);
	FlexModule_WriteCisGenomeIniFiles(strParamPath, strWorkPath, strOutFile, nMotifNum);

	/* return */
	return PROC_SUCCESS;
}


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
		char strPriorAbundance[])
{
	/* define */
	char strFilePath[MED_LINE_LENGTH];
	
	/* sequences and conservation */
	/* all sequences */
	struct FLEXSEQMOTIF **vSeqMtf = NULL;
	/* sequence number */
	int nSeqCount = 0;
	
	/* background */
	struct DOUBLEMATRIX *pBG = NULL;
	/* log background */
	struct DOUBLEMATRIX *pLogBG = NULL;
	/* background0 */
	struct DOUBLEMATRIX *pBG0 = NULL;
	/* log background0 */
	struct DOUBLEMATRIX *pLogBG0 = NULL;

	/* conservation likelihood */
	struct DOUBLEMATRIX *pCSLike = NULL;
	struct DOUBLEMATRIX **vCSLike = NULL;
	int nCSLikeNum = 0;

	/* motif matrix */
	struct FLEXMOTIFMATRIX **vMotif = NULL;
	/* maximum motif length */
	int nMaxMotifLen = 0;
	int nMotifLenLowerBound = FLEXMODULE_MOTIFLEN_LOWERBOUND;

	/* motif length distribution */
	struct INTMATRIX **vMotifLenSample = NULL;
	struct INTMATRIX *pMeanMotifLen = NULL;
	/* module size D,L,A distribution */
	struct DOUBLEMATRIX *pModuleDSample = NULL;
	struct DOUBLEMATRIX *pModuleLSample = NULL;
	struct DOUBLEMATRIX *pModuleASample = NULL;
	struct INTMATRIX *pModuleLenSample = NULL;
	/* motif number sample */
	struct DOUBLEMATRIX *pMotifNSample = NULL;

	/* number of base types */
	int nBaseTypeNum = 4;
	/* background base transition */
	int nScale = 0;

	/* module parameters */
	double dModuleL = 0.005;
	double dModuleA = 1.0;
	double dModuleD = dInitModuleD;
	int nMinModuleLen = 100;
	int nMaxModuleLen = 1000;

	/* frequency count */
	struct DOUBLEMATRIX *pFreqCount = NULL;
	struct DOUBLEMATRIX *pFreqPrior = NULL;
	double dDefaultPriorCount = 0.5;

	/* others */
	int ni,nj,nk,niter;
	double dT;
	double *pElement;
	/*double *pScoreEle;*/
	/*FILE *fpSeg;*/
	int nRecordNum,nWRelaxNum,nIsRecording;
	int nTempLen,nTotSamp;
	double dMeanLen;
	char strMotifFile[255];

	/* #################################### */
	/* initialize                           */
	/* #################################### */
	nBaseTypeNum = 4;
	nScale = (int)(pow((double)nBaseTypeNum, (double)(nBGOrder-1)));
	nRecordNum = (int)(nMCNum/4);
	nWRelaxNum = (int)(nMCNum/2);
	
	/* #################################### */
	/* load initial parameters              */
	/* #################################### */

	/* prior frequency */
	sprintf(strFilePath, "%s%s", strWorkPath, strPriorAbundance);
	pFreqPrior = NULL;
	pFreqPrior = DMLOAD(strFilePath);
	if(pFreqPrior == NULL)
	{
		printf("Error: FlexModule_Sampler, cannot load prior pseudocount of observing motifs!\n");
		exit(EXIT_FAILURE);
	}
	if( (pFreqPrior->nHeight != (nMotifNum+1)) || (pFreqPrior->nWidth != 2) )
	{
		printf("Error: FlexModule_Sampler, prior count matrix dimension does not match motif number!\n");
		exit(EXIT_FAILURE);
	}

	pFreqCount = NULL;
	pFreqCount = DMCLONE(pFreqPrior);
	if(pFreqCount == NULL)
	{
		printf("Error: FlexModule_Sampler, cannot initialize motif site counting matrix!\n");
		exit(EXIT_FAILURE);
	}

	/* prior conservation likelihood */
	pCSLike = NULL;
	vCSLike = NULL;
	if(nUseCS == 1)
	{
		sprintf(strFilePath, "%s%s", strWorkPath, strCSLike);
		pCSLike = DMLOAD(strFilePath);
		if(pCSLike == NULL)
		{
			printf("Error: FlexModule_Sampler, cannot load conservation likelihood matrix!\n");
			exit(EXIT_FAILURE);
		}

		vCSLike = FlexModule_InitScoreLike(pCSLike);
		if(vCSLike == NULL)
		{
			printf("Error: FlexModule_Sampler, cannot create conservation score likelihood function!\n");
			exit(EXIT_FAILURE);
		}

		nCSLikeNum = pCSLike->nHeight-1;
	}

	/* space for counting sampled module length */
	vMotifLenSample = NULL;
	vMotifLenSample = (struct INTMATRIX **)calloc(nMotifNum, sizeof(struct INTMATRIX *));
	if(vMotifLenSample == NULL)
	{
		printf("Error: FlexModule_Sampler, cannot create space for counting sampled motif length!\n");
		exit(EXIT_FAILURE);
	}

	pMeanMotifLen = NULL;
	pMeanMotifLen = CreateIntMatrix(1, nMotifNum);
	if(pMeanMotifLen == NULL)
	{
		printf("Error: FlexModule_Sampler, cannot create space for counting sampled motif length!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nMotifNum; ni++)
	{
		vMotifLenSample[ni] = NULL;
		vMotifLenSample[ni] = CreateIntMatrix(1, nMotifLenUpperBound);
		if(vMotifLenSample[ni] == NULL)
		{
			printf("Error: FlexModule_Sampler, cannot create space for counting sampled motif length!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* space for storing samples of module size parameter D */
	pModuleDSample = NULL;
	pModuleDSample = CreateDoubleMatrix(1, nRecordNum);
	if(pModuleDSample == NULL)
	{
		printf("Error: FlexModule_Sampler, cannot create space for storing sampled module parameter D!\n");
		exit(EXIT_FAILURE);
	}

	pModuleLSample = NULL;
	pModuleLSample = CreateDoubleMatrix(1, nRecordNum);
	if(pModuleLSample == NULL)
	{
		printf("Error: FlexModule_Sampler, cannot create space for storing sampled module parameter L!\n");
		exit(EXIT_FAILURE);
	}

	pModuleASample = NULL;
	pModuleASample = CreateDoubleMatrix(1, nRecordNum);
	if(pModuleASample == NULL)
	{
		printf("Error: FlexModule_Sampler, cannot create space for storing sampled module parameter A!\n");
		exit(EXIT_FAILURE);
	}

	pModuleLenSample = NULL;
	pModuleLenSample = CreateIntMatrix(1, nRecordNum);
	if(pModuleLenSample == NULL)
	{
		printf("Error: FlexModule_Sampler, cannot create space for storing sampled module length!\n");
		exit(EXIT_FAILURE);
	}

	/* space for counting true motif site number */
	pMotifNSample = NULL;
	pMotifNSample = CreateDoubleMatrix(1, nMotifNum);
	if(pMotifNSample == NULL)
	{
		printf("Error: FlexModule_Sampler, cannot create space for storing sampled motif site number!\n");
		exit(EXIT_FAILURE);
	}

	/* prior motif PWM */
	vMotif = NULL;
	vMotif = (struct FLEXMOTIFMATRIX **)calloc(nMotifNum, sizeof(struct FLEXMOTIFMATRIX *));
	if(vMotif == NULL)
	{
		printf("Error: FlexModule_Sampler, cannot create motif PWM matrix!\n");
		exit(EXIT_FAILURE);
	}

	nMaxMotifLen = 0;
	for(ni=0; ni<nMotifNum; ni++)
	{
		/* create motif */
		vMotif[ni] = NULL;
		vMotif[ni] = FLEXMTFMCREATE();
		if( vMotif[ni] == NULL)
		{
			printf("Error: FlexModule_Sampler, cannot load init motif matrices!\n");
			exit(EXIT_FAILURE);
		}

		/* motif id */
		vMotif[ni]->nMotifType = (ni+1);

		/* if no user-specified prior */
		if(	strcmp(vPriorPWMFile[ni]->m_pString, "NULL") == 0 )
		{
			vMotif[ni]->pPriorCount = NULL;
			vMotif[ni]->pPriorCount = CreateDoubleMatrix(pMotifLen->pMatElement[ni], nBaseTypeNum);
			if(vMotif[ni]->pPriorCount == NULL)
			{
				printf("Error: FlexModule_Sampler, cannot create prior count matrix!\n");
				exit(EXIT_FAILURE);
			}

			pElement = vMotif[ni]->pPriorCount->pMatElement;
			for(nj=0; nj<vMotif[ni]->pPriorCount->nHeight; nj++)
			{
				for(nk=0; nk<vMotif[ni]->pPriorCount->nWidth; nk++)
				{
					*pElement = dDefaultPriorCount;
					pElement++;
				}
			}
		}
		/* if user-specified prior */
		else
		{			
			/* init prior matrix */
			sprintf(strFilePath, "%s%s", strWorkPath, vPriorPWMFile[ni]->m_pString);
			vMotif[ni]->pPriorCount = NULL;
			vMotif[ni]->pPriorCount = DMLOAD(strFilePath);
			if( vMotif[ni]->pPriorCount == NULL)
			{
				printf("Error: FlexModule_Sampler, cannot load init motif pseudocount!\n");
				exit(EXIT_FAILURE);
			}
			if(	vMotif[ni]->pPriorCount->nHeight != pMotifLen->pMatElement[ni] )
			{
				printf("Error: FlexModule_Sampler, init motif length does not match the length of init motif matrix!\n");
				exit(EXIT_FAILURE);
			}
			if(	vMotif[ni]->pPriorCount->nWidth != nBaseTypeNum )
			{
				printf("Error: FlexModule_Sampler, the width of init motif matrix shoule be %d!\n", nBaseTypeNum);
				exit(EXIT_FAILURE);
			}
		}

		/* init sample count matrix */
		vMotif[ni]->pSampleCount = NULL;
		vMotif[ni]->pSampleCount = CreateDoubleMatrix(pMotifLen->pMatElement[ni], nBaseTypeNum);
		if(vMotif[ni]->pSampleCount == NULL)
		{
			printf("Error: FlexModule_Sampler, cannot create sample count matrix!\n");
			exit(EXIT_FAILURE);
		}

		/* init PWM matrix */
		vMotif[ni]->pPWM = NULL;
		vMotif[ni]->pPWM = CreateDoubleMatrix(pMotifLen->pMatElement[ni], nBaseTypeNum);
		if(vMotif[ni]->pPWM == NULL)
		{
			printf("Error: FlexModule_Sampler, cannot create PWM matrix!\n");
			exit(EXIT_FAILURE);
		}
		FLEXMTFMREFRESH(*(vMotif+ni));

		/* update the maximal motif length */
		if(pMotifLen->pMatElement[ni] > nMaxMotifLen)
			nMaxMotifLen = pMotifLen->pMatElement[ni];
	}


	/* #################################### */
	/* load sequences                       */
	/* #################################### */
	vSeqMtf = NULL;
	vSeqMtf = FLEXSEQMTFLOADSEQ(&nSeqCount, strWorkPath, strSeqFile, nUseCS, strCSPrefix, nMotifNum);

	/* #################################### */
	/* compute background markov chain      */
	/* #################################### */
	pBG = NULL;
	pBG = FLEXSEQMTFESTIMATENUCLEICBGMC(vSeqMtf, nSeqCount, 0, nBGOrder, nBaseTypeNum);
	if(pBG == NULL)
	{
		printf("Error: FlexModule_Sampler, failure to estimate background markov matrix!\n");
		exit(EXIT_FAILURE);
	}
	pLogBG = NULL;
	pLogBG = DMCLONE(pBG);
	pElement = pLogBG->pMatElement;
	for(ni=0; ni<pLogBG->nHeight; ni++)
	{
		for(nj=0; nj<pLogBG->nWidth; nj++)
		{
			*pElement = log(*pElement);
			pElement++;
		}
	}
	/* DMSAVE(pBG, "bgestimate.txt"); */
	
	pBG0 = NULL;
	pBG0 = FLEXSEQMTFESTIMATENUCLEICBGMC(vSeqMtf, nSeqCount, 0, 0, nBaseTypeNum);
	if(pBG0 == NULL)
	{
		printf("Error: FlexModule_Sampler, failure to estimate background markov matrix!\n");
		exit(EXIT_FAILURE);
	}
	pLogBG0 = NULL;
	pLogBG0 = DMCLONE(pBG0);
	pElement = pLogBG0->pMatElement;
	for(ni=0; ni<pLogBG0->nHeight; ni++)
	{
		for(nj=0; nj<pLogBG0->nWidth; nj++)
		{
			*pElement = log(*pElement);
			pElement++;
		}
	}
	/* DMSAVE(pBG0, "bgestimate0.txt"); */
	
	/* #################################### */
	/* initialize the sampler               */
	/* #################################### */
	/* init background loglikelihood */
	if(nUseFittedBG == 1)
	{
		 for(ni=0; ni<nSeqCount; ni++)
		 {
			sprintf(strFilePath, "%s%s_%d_%s.bgf", strWorkPath, 
				strFittedBGPrefix, ni, vSeqMtf[ni]->strAlias); 
			FlexModule_InitBGLogLike_FromFittedBGFile(vSeqMtf[ni], 0,
				strFilePath, nUseCS, nCSLikeNum, vCSLike); 
		 }
	}
	else
	{
		for(ni=0; ni<nSeqCount; ni++)
		{
			FlexModule_InitBGLogLike(vSeqMtf[ni], 0,
				nBGOrder, pLogBG, pLogBG0, 
				nUseCS, nCSLikeNum, vCSLike);
		}
	}

	/* init sites */
	for(ni=0; ni<nSeqCount; ni++)
	{
 		FlexModule_InitSite(vSeqMtf[ni], 0,
			pFreqPrior, pFreqCount,
			nBGOrder, pLogBG, pLogBG0,  
			nMotifNum, vMotif, pMotifLen, nMaxMotifLen,
			nUseCS, nCSLikeNum, vCSLike);

		FlexModule_InitContext(vSeqMtf[ni], 0, 0.5,
			pFreqCount, dModuleD, dModuleL, dModuleA, nModuleLen,
			nMotifNum, vMotif, pMotifLen, nMaxMotifLen);
	}

	/* DMSAVE(pFreqCount, "temp.txt");
	DMSAVE(vMotif[0]->pSampleCount, "motif0.txt");
	DMSAVE(vMotif[1]->pSampleCount, "motif1.txt");
	DMSAVE(vMotif[2]->pSampleCount, "motif2.txt"); */

	/* update motif matrix */
	for(ni=0; ni<nMotifNum; ni++)
	{
		FLEXMTFMREFRESH(vMotif[ni]);
	}

	/* #################################### */
	/* flexmodule sampler                   */
	/* #################################### */
	for(niter=0; niter<nMCNum; niter++)
	{
		if(niter%50 == 0)
		{
			printf("iter=%d...\n", niter);
		}

		dT = 2.0*(double)niter/(double)nMCNum;
		if(dT > 1.0)
			dT = 1.0;

		if((niter+nWRelaxNum) < nMCNum)
			nIsRecording = 0;
		else if((niter+nRecordNum) < nMCNum)
			nIsRecording = 1;
		else
			nIsRecording = 2;

		/* sequence-wise site sampling */
		for(ni=0; ni<nSeqCount; ni++)
		{
			FlexModule_SampleBA(vSeqMtf[ni], 0, dT, 
				         pFreqCount, dModuleD, dModuleL, dModuleA, nModuleLen, 
						 nBGOrder, pLogBG, pLogBG0,
						 nMotifNum, vMotif, pMotifLen, nMaxMotifLen, 
						 nUseCS, nCSLikeNum, vCSLike,
						 nIsRecording, pMotifNSample);
		}

		/* recording current samples */
		if(nIsRecording == 2)
		{
			for(ni=0; ni<nMotifNum; ni++)
			{
				nTempLen = pMotifLen->pMatElement[ni]-1;
				vMotifLenSample[ni]->pMatElement[nTempLen] += 1;
			}
		}

		/* adjust sampled motifs */
		if(nIsRecording > 0)
		{
			/* sample module size factor D */
			FlexModule_SampleD(vSeqMtf, nSeqCount, 0, dT, dInitModuleD,
				&dModuleD, &dModuleL, &dModuleA, 
				nModuleLen, nMotifNum, vMotif, pMotifLen,
				nMaxMotifLen, nIsRecording);
			

			/* sample module length */
			if(nSampleModuleLen == 1)
			{
				FlexModule_SampleL(vSeqMtf, nSeqCount, 0, dT, 
				dModuleD, dModuleL, dModuleA, 
				&nModuleLen, nMinModuleLen, nMaxModuleLen, 
				nMotifNum, vMotif, pMotifLen,
				nMaxMotifLen, nIsRecording);
			}
		
			if(nIsRecording == 2)
			{
				pModuleDSample->pMatElement[niter+nRecordNum-nMCNum] = dModuleD;
				pModuleLSample->pMatElement[niter+nRecordNum-nMCNum] = dModuleL;
				pModuleASample->pMatElement[niter+nRecordNum-nMCNum] = dModuleA;
				pModuleLenSample->pMatElement[niter+nRecordNum-nMCNum] = nModuleLen;
			}

			/* sample motif length */
			if(niter%2 == 0)
			{
				for(ni=0; ni<nMotifNum; ni++)
				{
					/* change motif length */
					FlexModule_SampleW(vSeqMtf, nSeqCount, 0, dT, 
						pFreqCount, dModuleD, dModuleL, dModuleA, 
						nModuleLen,
						nBGOrder, pLogBG, pLogBG0,
						(ni+1), nMotifNum, vMotif, pMotifLen, &nMaxMotifLen,
						nUseCS, nCSLikeNum, vCSLike,
						nIsRecording, 
						nMeanMotifLen, nMotifLenUpperBound, nMotifLenLowerBound,
						dDefaultPriorCount);
				}
			}

			/* local shift */
			if(niter%2 == 1)
			{
				for(ni=0; ni<nMotifNum; ni++)
				{
					/* local mode shift */
					FlexModule_LocalShift(vSeqMtf, nSeqCount, 0, dT, 
						pFreqCount, dModuleD, dModuleL, dModuleA, 
						nModuleLen,
						nBGOrder, pLogBG, pLogBG0,
						(ni+1), vMotif, pMotifLen, nMaxMotifLen,
						nUseCS, nCSLikeNum, vCSLike,
						nIsRecording);
				}
			}
		}
	}

	/* #################################### */
	/* get posterior mean lengths of motifs */
	/* #################################### */
	for(ni=0; ni<nMotifNum; ni++)
	{
		pMotifNSample->pMatElement[ni] /= (double)nRecordNum;
		nTempLen = (int)pMotifNSample->pMatElement[ni];
		if(pMotifNSample->pMatElement[ni] - nTempLen >= 0.5)
		{
			pMotifNSample->pMatElement[ni] = nTempLen+1.0;
		}
		else
		{
			pMotifNSample->pMatElement[ni] = nTempLen;
		}

		nTempLen = 0;
		nTotSamp = 0;
		for(nj=0; nj<nMotifLenUpperBound; nj++)
		{
			nTempLen += (vMotifLenSample[ni]->pMatElement[nj])*(nj+1);
			nTotSamp += vMotifLenSample[ni]->pMatElement[nj];
		}
		
		dMeanLen = (double)nTempLen/(double)nTotSamp;
		pMeanMotifLen->pMatElement[ni] = (int)dMeanLen;
		if(dMeanLen - pMeanMotifLen->pMatElement[ni] >= 0.5)
		{
			pMeanMotifLen->pMatElement[ni] += 1;
		}
	}

	/* #################################### */
	/* export sampled parameters            */
	/* #################################### */
	sprintf(strFilePath, "%s%s_modD.txt", strWorkPath, strOutFile);
	DMSAVE(pModuleDSample, strFilePath);
	sprintf(strFilePath, "%s%s_modL.txt", strWorkPath, strOutFile);
	DMSAVE(pModuleLSample, strFilePath);
	sprintf(strFilePath, "%s%s_modA.txt", strWorkPath, strOutFile);
	DMSAVE(pModuleASample, strFilePath);
	sprintf(strFilePath, "%s%s_modLen.txt", strWorkPath, strOutFile);
	IMSAVE(pModuleLenSample, strFilePath);

	/* #################################### */
	/* export last sample                   */
	/* #################################### */
	for(ni=0; ni<nSeqCount; ni++)
	{
		FlexModule_CallMotifSitesInFinalSample(vSeqMtf[ni], pMotifLen);
	}
	sprintf(strFilePath, "%s%s_l.txt", strWorkPath, strOutFile);
	FlexModule_WriteMotifToFile(strFilePath, nMotifNum, vMotif, pBG0, pMotifLen,
		nSeqCount, vSeqMtf, 0);

	/* #################################### */
	/* export posterior sample              */
	/* #################################### */
	sprintf(strFilePath, "%s%s_p.txt", strWorkPath, strOutFile);
	FlexModule_WritePostMotifToFile(strFilePath, nMotifNum, vMotif, pBG0, pLogBG0, 
		pMotifLen, pMeanMotifLen, pMotifNSample, 
		nSeqCount, vSeqMtf, 0, nRecordNum, dDefaultPriorCount);

	/* #################################### */
	/* release memory                       */
	/* #################################### */
	DestroyDoubleMatrix(pFreqPrior);
	DestroyDoubleMatrix(pFreqCount);
	
	if(nUseCS == 1)
	{
		for(ni=0; ni<(pCSLike->nHeight-1); ni++)
		{
			DestroyDoubleMatrix(vCSLike[ni]);
			vCSLike[ni] = NULL;
		}
		free(vCSLike);
		DestroyDoubleMatrix(pCSLike);
	}
	
	for(ni=0; ni<nMotifNum; ni++)
	{
		FLEXMTFMDESTROY(vMotif[ni]);
		vMotif[ni] = NULL;
	}
	free(vMotif);

	for(ni=0; ni<nMotifNum; ni++)
	{
		DestroyIntMatrix(vMotifLenSample[ni]);
		vMotifLenSample[ni] = NULL;
	}
	free(vMotifLenSample);

	DestroyIntMatrix(pMeanMotifLen);

	DestroyDoubleMatrix(pModuleDSample);
	DestroyDoubleMatrix(pModuleLSample);
	DestroyDoubleMatrix(pModuleASample);
	DestroyIntMatrix(pModuleLenSample);

	DestroyDoubleMatrix(pMotifNSample);

	for(ni=0; ni<nSeqCount; ni++)
	{
		FLEXSEQMTFDESTROY(vSeqMtf[ni]);
		vSeqMtf[ni] = NULL;
	}
	free(vSeqMtf);

	DestroyDoubleMatrix(pBG);
	DestroyDoubleMatrix(pLogBG);
	DestroyDoubleMatrix(pBG0);
	DestroyDoubleMatrix(pLogBG0);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_InitScoreLike: prepare likelihood function for scores       */
/*  sample cis-regulatory module with conservation scores.                 */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX **FlexModule_InitScoreLike(struct DOUBLEMATRIX *pScoreLike)
{
	/* define */
	struct DOUBLEMATRIX **vL;
	int nScoreNum;
	int ni,nj,nid;
	double dIntS,dIntE,dInitStep;

	/* init */
	if(pScoreLike->nWidth <= 2)
	{
		printf("Error: FlexModule_InitScoreLike, score likelihood function needs to be set for more than two intervals!\n");
		exit(EXIT_FAILURE);
	}
	dIntS = DMGETAT(pScoreLike, 0, 0);
	dIntE = DMGETAT(pScoreLike, 0, (pScoreLike->nWidth-2));
	
	dInitStep = (dIntE-dIntS)/(double)(pScoreLike->nWidth-2);
	if(dInitStep < 0.0)
	{
		printf("Error: FlexModule_InitScoreLike, negative interval step length!\n");
		exit(EXIT_FAILURE);
	}

	vL = NULL;
	vL = (struct DOUBLEMATRIX **)calloc(2, sizeof(struct DOUBLEMATRIX *));
	if(vL == NULL)
	{
		printf("Error: FlexModule_InitScoreLike, cannot create the matrix for likelihood function!\n");
		exit(EXIT_FAILURE);
	}

	/* create matrix */
	nScoreNum = pScoreLike->nHeight-1;
	for(ni=0; ni<nScoreNum; ni++)
	{
		vL[ni] = NULL;
		vL[ni] = CreateDoubleMatrix(1, FLEXMODULE_SCOREBIN);
		if(vL[ni] == NULL)
		{
			printf("Error: FlexModule_InitScoreLike, cannot create the matrix for likelihood function!\n");
			exit(EXIT_FAILURE);
		}

		for(nj=0; nj<FLEXMODULE_SCOREBIN; nj++)
		{
			nid = (int)(nj/dInitStep);
			if(nid >= pScoreLike->nWidth)
				nid = pScoreLike->nWidth-1;

			vL[ni]->pMatElement[nj] = DMGETAT(pScoreLike, (ni+1), nid);
			if(vL[ni]->pMatElement[nj] > 0.0)
			{
				vL[ni]->pMatElement[nj] = log(vL[ni]->pMatElement[nj]);
			}
			else
			{
				printf("Error: FlexModule_InitScoreLike, log(0), likelihood cannot be zero!\n");
				exit(EXIT_FAILURE);
			}
		}

	}

	/* return */
	return vL;
}

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_InitBGLogLike: Initialize background log likelihood         */
/* ----------------------------------------------------------------------- */ 
int FlexModule_InitBGLogLike(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,
			int nBGOrder, struct DOUBLEMATRIX *pLogBG, struct DOUBLEMATRIX *pLogBG0,
			int nUseCS, int nCSLikeNum, struct DOUBLEMATRIX **vCSLike)
{
	/* define */
	int nBaseTypeNum,nScale;
	int nWordId,nWordBadLen;
	int nW,nCSId;
	int nLen;
	unsigned char *pBase;
	/* unsigned char *pCS; */
	double *pBLike;
	int ni;

	/* init */
	nLen = pSeqMtf->vSeq[nIndex]->nWidth;
	pBase = pSeqMtf->vSeq[nIndex]->pMatElement;
	/* pCS = pSeqMtf->vScore[0]->pMatElement; */
	pBLike = pSeqMtf->vMonitor[0]->pMatElement;

	nBaseTypeNum = pLogBG->nWidth;
	if( nBGOrder == 0)
	{
		nScale = 0;
	}
	else
	{
		nScale = (int)(pow((double)nBaseTypeNum, (double)(nBGOrder-1)));
	}
	
	/* if zero order Markov chain as background */
	if(nBGOrder == 0)
	{
		for(ni=0; ni<nLen; ni++)
		{
			nW = (int)pBase[ni];
			if(nW < nBaseTypeNum)
			{
				pBLike[ni] = DMGETAT(pLogBG, 0, nW);
				/* if(nUseCS == 1)
				{
					nCSId = pCS[ni];
					pBLike[ni] += vCSLike[0]->pMatElement[nCSId];
				} */
			}
			else
			{
				pBLike[ni] = 0.0;
			}
		}
	}

	/* if higher order Markov chain as background */
	else
	{
		nWordId = 0;
		nWordBadLen = 0;
		
		/* initial words */
		for(ni=0; ni<nBGOrder; ni++)
		{
			if(ni >= nLen)
				break;

			nW = (int)pBase[ni];
			if(nW < nBaseTypeNum)
			{
				pBLike[ni] = DMGETAT(pLogBG0, 0, nW);
				/* if(nUseCS == 1)
				{
					nCSId = pCS[ni];
					pBLike[ni] += vCSLike[0]->pMatElement[nCSId];
				} */
				nWordId = nWordId*nBaseTypeNum+nW;
			}
			else
			{
				pBLike[ni] = 0.0;
				nWordId = nWordId*nBaseTypeNum;
				nWordBadLen++;
			}
		}

		/* continuing words */
		for(; ni<nLen; ni++)
		{
			/* get likelihood */
			nW = (int)pBase[ni];
			if(nW < nBaseTypeNum)
			{
				if(nWordBadLen == 0)
					pBLike[ni] = DMGETAT(pLogBG, nWordId, nW);
				else
					pBLike[ni] = DMGETAT(pLogBG0, 0, nW);
				
				/* if(nUseCS == 1)
				{
					nCSId = pCS[ni];
					pBLike[ni] += vCSLike[0]->pMatElement[nCSId];
				} */
			}
			else
			{
				pBLike[ni] = 0.0;
			}
			
			/* update word id */
			/* minus old letter */
			nW = (int)(pBase[ni-nBGOrder]);
			if(nW < nBaseTypeNum)
			{
				nWordId -= nW*nScale;
			}
			else
			{
				nWordBadLen--;
			} 
			
			/* add new letter */
			nW = (int)pBase[ni];
			if(nW < nBaseTypeNum)
			{
				nWordId = nWordId*nBaseTypeNum+nW;
			}
			else
			{
				nWordId = nWordId*nBaseTypeNum;
				nWordBadLen++;
			}
		}
	}


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_FromFittedBGFile: Initialize background log likelihood from */
/*  fitted background files.                                               */
/* ----------------------------------------------------------------------- */ 
int FlexModule_InitBGLogLike_FromFittedBGFile(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,
			char strFilePath[],	int nUseCS, int nCSLikeNum, struct DOUBLEMATRIX **vCSLike)
{
	/* define */
	int nLen;
	double *pBLike;
	int ni,nj;
	FILE *fpIn;
	int nTemp;
	char strITemp[INT_SIZE];
	char strITemp2[INT_SIZE];
	char strDTemp[DOUBLE_SIZE];
	char strDTemp2[DOUBLE_SIZE];
	
	/* FILE *fpOut;
	char strTempFile[LINE_LENGTH];
	sprintf(strTempFile, "%s.txt", strFilePath);
	fpOut = NULL;
	fpOut = fopen(strTempFile, "w");
	*/

	/* init */
	nLen = pSeqMtf->vSeq[nIndex]->nWidth;
	pBLike = pSeqMtf->vMonitor[0]->pMatElement;

	/* open file */
	fpIn = NULL;
	fpIn = fopen(strFilePath, "rb");
	if(fpIn == NULL)
	{
		printf("Error: FlexModule_InitBGLogLike_FromFittedBGFile, cannot open input file\n");
		exit(EXIT_FAILURE);
	}

	/* check length */
	fread(&nTemp, sizeof(int), 1, fpIn);
	if(nTemp == nLen)
	{
		fread(pBLike, sizeof(double), nLen, fpIn);

		/* for(nj=0; nj<nLen; nj++)
		{
			fprintf(fpOut, "%f\n", pBLike[nj]);
		} */
	}
	else
	{
		if(fseek( fpIn, 0, SEEK_SET ) != 0)
		{
			printf("Error: FlexModule_InitBGLogLike_FromFittedBGFile, cannot locate the required sequence!\n");
			exit(EXIT_FAILURE);
		}

		fread(strITemp, sizeof(char), INT_SIZE, fpIn);
		for(ni=0; ni<INT_SIZE; ni++)
			strITemp2[ni] = strITemp[INT_SIZE-1-ni];
		memcpy(&nTemp, strITemp2, INT_SIZE);

		if(nTemp != nLen)
		{
			printf("Error: FlexModule_InitBGLogLike_FromFittedBGFile, cannot load fitted background!\n");
			exit(EXIT_FAILURE);
		}

		for(nj=0; nj<nLen; nj++)
		{
			fread(strDTemp, sizeof(char), DOUBLE_SIZE, fpIn);
			for(ni=0; ni<DOUBLE_SIZE; ni++)
				strDTemp2[ni] = strDTemp[DOUBLE_SIZE-1-ni];
			memcpy(pBLike+nj, strDTemp2, DOUBLE_SIZE);

			/* fprintf(fpOut, "%f\n", pBLike[nj]); */
		}
	}
	
	/* close files */
	fclose(fpIn);
	/* fclose(fpOut); */

	/* return */
	return PROC_SUCCESS;
}



/* ----------------------------------------------------------------------- */ 
/*  FlexModule_InitSite: Initialize motif sites at random.                 */
/* ----------------------------------------------------------------------- */ 
int FlexModule_InitSite(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,
			struct DOUBLEMATRIX *pFreqPrior0, struct DOUBLEMATRIX *pFreqCount,
			int nBGOrder, struct DOUBLEMATRIX *pLogBG, struct DOUBLEMATRIX *pLogBG0,  
			int nMotifNum, struct FLEXMOTIFMATRIX **vMotif, 
			struct INTMATRIX *pMotifLen, int nMaxMotifLen,
			int nUseCS, int nCSLikeNum, struct DOUBLEMATRIX **vCSLike)
{
	/* define */
	/* sampler probability */
	struct DOUBLEMATRIX *pFreqPrior = NULL;
	struct DOUBLEMATRIX *pFreqPost = NULL;
	struct BYTEMATRIX *pFreqValid = NULL;
	
	/* sampler status */
	int ni,nj,nx,nLen;
	unsigned char *pSi,*pBase;
	struct FLEXMOTIFSITE *pNewSite,*pCurrentSite;
	struct FLEXMOTIFSITE **vSite;

	/* sampler likelihood */
	int nBaseTypeNum;
	int nStateId, nStrand;
	int nMotifId;
	double dTemp;
	double dBaseLike;
	double dLike;
	int nMask;

	/* init */
	pFreqPrior = NULL;
	pFreqPrior = DMCLONE(pFreqPrior0); 
	if(pFreqPrior == NULL)
	{
		printf("Error: FlexModule_InitSite, cannot create prior probability matrix!\n");
		exit(EXIT_FAILURE);
	}

	pFreqPost = NULL;
	pFreqPost = CreateDoubleMatrix(pFreqCount->nHeight, pFreqCount->nWidth); 
	if(pFreqPost == NULL)
	{
		printf("Error: FlexModule_InitSite, cannot create posterior probability matrix!\n");
		exit(EXIT_FAILURE);
	}

	pFreqValid = NULL;
	pFreqValid = CreateByteMatrix(pFreqCount->nHeight, pFreqCount->nWidth); 
	if(pFreqValid == NULL)
	{
		printf("Error: FlexModule_InitSite, cannot create sequence status indicator matrix!\n");
		exit(EXIT_FAILURE);
	}

	/* get prior probability */
	FlexModule_NormalizeFreq(pFreqPrior);
	FlexModule_LogFreq(pFreqPrior);
	
	/* init */
	nLen = pSeqMtf->vSeq[nIndex]->nWidth;
	pBase = pSeqMtf->vSeq[nIndex]->pMatElement;
	pSi = pSeqMtf->vStatus[1]->pMatElement;
	vSite = pSeqMtf->vMotif;
	pCurrentSite = pSeqMtf->pMotifList;
	
	nBaseTypeNum = pLogBG->nWidth;
	
	/* ######################### */
	/* scan the whole sequence   */
	/* ######################### */
	
	/* set background likelihood */
	dBaseLike = DMGETAT(pFreqPrior, 0, 0);
	BMSETAT(pFreqValid, 0, 0, 1);
		
	for(ni=0; ni<nLen; ni++)
	{
		/* if repeat or 'N', skip it */
		if(pBase[ni] >= nBaseTypeNum)
			continue;

		/* set background likelihood */
		DMSETAT(pFreqPost, 0, 0, dBaseLike);
	
		/* get motif likelihood */
		for(nj=0; nj<nMotifNum; nj++)
		{
			/* '+' strand */
			dLike = FlexBaseLogLikelihood_Motif(pSeqMtf, nIndex, ni, 
				vMotif[nj]->pPWM, '+', nMaxMotifLen,
				nBGOrder, pLogBG, pLogBG0, 
				nUseCS, nCSLikeNum, vCSLike,
				&nMask);

			if(nMask == 0)
				dLike += DMGETAT(pFreqPrior, (nj+1), 0);
			BMSETAT(pFreqValid, nj+1, 0, (unsigned char)(1-nMask));
			DMSETAT(pFreqPost, nj+1, 0, dLike);

			/* '-' strand */
			dLike = FlexBaseLogLikelihood_Motif(pSeqMtf, nIndex, ni, 
				vMotif[nj]->pPWM, '-', nMaxMotifLen,
				nBGOrder, pLogBG, pLogBG0, 
				nUseCS, nCSLikeNum, vCSLike,
				&nMask);

			if(nMask == 0)
				dLike += DMGETAT(pFreqPrior, (nj+1), 1);
			BMSETAT(pFreqValid, nj+1, 1, (unsigned char)(1-nMask));
			DMSETAT(pFreqPost, nj+1, 1, dLike);
		}

		/* normalize posterior probability */
		FlexModule_NormalizePost(pFreqPost, pFreqValid);
		
		/* sample motifs */
		FlexModule_SamplePost(pFreqPost, pFreqValid, &nStateId, &nStrand);

		if(nStateId > 0)
		{
			/* generate pNewSite */
			pNewSite = NULL;
			pNewSite = FLEXMTFSCREATE();
			if(pNewSite == NULL)
			{
				printf("Error: FlexModule_InitSite, cannot create new binding site!\n");
				exit(EXIT_FAILURE);
			}
			pNewSite->nMotifType = nStateId;
			pNewSite->nSeqId = pSeqMtf->nId;
			pNewSite->nStartPos = ni;
			pNewSite->nStrand = nStrand;
			
			vSite[ni] = pNewSite;
			pSeqMtf->nSiteNum += 1;
			if(pSeqMtf->pMotifList == NULL)
			{
				pSeqMtf->pMotifList = pNewSite;
				pCurrentSite = pNewSite;
			}
			else
			{
				pCurrentSite->pNext = pNewSite;
				pNewSite->pPrev = pCurrentSite;
				pCurrentSite = pNewSite;
			}

			
			/* update pSi */
			nMotifId = nStateId-1;
			for(nx=0; nx <pMotifLen->pMatElement[nMotifId]; nx++)
			{
				pSi[ni+nx] = nStateId;
			}
		}

		/* update count */
		dTemp = DMGETAT(pFreqCount, nStateId, nStrand)+1.0;
		DMSETAT(pFreqCount, nStateId, nStrand, dTemp);
	}
	
	/* release memory */
	DestroyDoubleMatrix(pFreqPrior);
	DestroyDoubleMatrix(pFreqPost);
	DestroyByteMatrix(pFreqValid);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_InitContext: create initial context.                        */
/* ----------------------------------------------------------------------- */ 
int FlexModule_InitContext(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, double dT,
			struct DOUBLEMATRIX *pFreqCount, double dInitModuleD, 
			double dModuleL, double dModuleA, int nModuleLen,
			int nMotifNum, struct FLEXMOTIFMATRIX **vMotif, 
			struct INTMATRIX *pMotifLen, int nMaxMotifLen)
{
	/* define */
	struct FLEXMOTIFSITE *pSite;
	struct FLEXMOTIFSITE *pHead,*pTail,*pTest;
	unsigned char *vAi;
	int nContextNum;
	double dScore,dRand;
	int ni,nMotifId;
	int nHalfWin;
	double dPAi;
	
	/* init */
	nHalfWin = (int)(nModuleLen/2);
	vAi = pSeqMtf->vStatus[0]->pMatElement;

	/* create context one by one */
	pSite = pSeqMtf->pMotifList;
	while(pSite != NULL)
	{
		/* establish context limit */
		nContextNum = 1;
		pHead = pSite;
		pTail = pSite;

		pTest = pHead->pPrev;
		while(pTest != NULL)
		{
			if( (pSite->nStartPos-pTest->nStartPos) > nHalfWin )
				break;

			pHead = pTest;
			pTest = pHead->pPrev;
			nContextNum++;
		}

		pTest = pTail->pNext;
		while(pTest != NULL)
		{
			if( (pTest->nStartPos-pSite->nStartPos) > nHalfWin )
				break;

			pTail = pTest;
			pTest = pTail->pNext;
			nContextNum++;
		}

		/* create context matrix */
		pSite->pContextPos = NULL;
		pSite->pContextPos = CreateIntMatrix(1, nContextNum);
		if(pSite->pContextPos == NULL)
		{
			printf("Error: FlexModule_InitContext, cannot create context matrix correctly!\n");
			exit(EXIT_FAILURE);
		}
		pSite->pContextMotif = NULL;
		pSite->pContextMotif = CreateIntMatrix(1, nContextNum);
		if(pSite->pContextMotif == NULL)
		{
			printf("Error: FlexModule_InitContext, cannot create context matrix correctly!\n");
			exit(EXIT_FAILURE);
		}
		pSite->pContextStrand = NULL;
		pSite->pContextStrand = CreateIntMatrix(1, nContextNum);
		if(pSite->pContextStrand == NULL)
		{
			printf("Error: FlexModule_InitContext, cannot create context matrix correctly!\n");
			exit(EXIT_FAILURE);
		}

		/* fill out the context */
		pTest = pHead;
		ni = 0;
		while(pTest != NULL)
		{
			if(pTest->nStartPos - pTail->nStartPos > 0)
				break;

			pSite->pContextPos->pMatElement[ni] = pTest->nStartPos-pSite->nStartPos;
			pSite->pContextMotif->pMatElement[ni] = pTest->nMotifType;
			pSite->pContextStrand->pMatElement[ni] = pTest->nStrand;
			ni++;
			pTest = pTest->pNext;
		}

		if(ni != nContextNum)
		{
			printf("Error: FlexModule_InitContext, number of motif sites in the context do not match!\n");
			exit(EXIT_FAILURE);
		}

		/* get Di */
		dScore = FlexModule_ComputeModuleScore(nContextNum, pSite->pContextPos->pMatElement,
			pSite->pContextMotif->pMatElement, pSite->pContextStrand->pMatElement,
			0, 0, 0, 0, nMotifNum);

		/* sample Ai */
		dPAi = FlexModule_ComputeTrueProb(dScore, dInitModuleD, dModuleL, dModuleA, dT);
		dRand = rand_u();

		/* update motif PWM count */
		if(dRand <= dPAi)
		{
			vAi[pSite->nStartPos] = 1;
			nMotifId = pSite->nMotifType-1;
			if(pSite->nStrand == 1)
			{
				FlexModule_AddMotifMatrixCount(vMotif[nMotifId]->pSampleCount, 
						  pSeqMtf->vSeq[nIndex]->pMatElement+pSite->nStartPos, 
						  pMotifLen->pMatElement[nMotifId], '-');
			}
			else
			{
				FlexModule_AddMotifMatrixCount(vMotif[nMotifId]->pSampleCount, 
						  pSeqMtf->vSeq[nIndex]->pMatElement+pSite->nStartPos, 
						  pMotifLen->pMatElement[nMotifId], '+');
			}
		}
		else
		{
			vAi[pSite->nStartPos] = 0;
		}

		/* get next */
		pSite = pSite->pNext;
	}

	/* return */
	return PROC_SUCCESS;
}

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
						 int nIsRecording, struct DOUBLEMATRIX *pMotifNSample)
{
	/* define */
	/* sampler probability */
	struct DOUBLEMATRIX *pFreqPost = NULL;
	struct BYTEMATRIX *pFreqValid = NULL;

	/* context sites */
	struct FLEXMOTIFSITE *pPrev,*pNext,*pHead,*pTail,*pTest,*pCurrent;

	/* sampler status */
	int ni,nj,nk,nLen;
	unsigned char *pAi,*pSi,*pBase;
	struct FLEXMOTIFSITE **vSite;
	struct FLEXMOTIFSITE *pContext;

	/* sampler likelihood */
	int nBaseTypeNum;
	double dTemp;
	int nMotifId;
	int nHalfWin;
	int nStateId, nStrand;
	int nContextNum;
	double dRand;
	
	double dLike,dScore,dPAi,dContextLike;
	int nMask;

	/* init */
	nHalfWin = (int)(nModuleLen/2);
	pFreqPost = NULL;
	pFreqPost = CreateDoubleMatrix(pFreqCount->nHeight, pFreqCount->nWidth); 
	if(pFreqPost == NULL)
	{
		printf("Error: FlexModule_SampleBA, cannot create posterior probability matrix!\n");
		exit(EXIT_FAILURE);
	}

	pFreqValid = NULL;
	pFreqValid = CreateByteMatrix(pFreqCount->nHeight, pFreqCount->nWidth); 
	if(pFreqValid == NULL)
	{
		printf("Error: FlexModule_SampleBA, cannot create sequence status indicator matrix!\n");
		exit(EXIT_FAILURE);
	}

	/* init */
	nLen = pSeqMtf->vSeq[nIndex]->nWidth;
	pBase = pSeqMtf->vSeq[nIndex]->pMatElement;
	pAi = pSeqMtf->vStatus[0]->pMatElement;
	pSi = pSeqMtf->vStatus[1]->pMatElement;
	vSite = pSeqMtf->vMotif;
	pCurrent = pSeqMtf->pMotifList;
	
	nBaseTypeNum = pLogBG->nWidth;
	
	/* ######################### */
	/* scan the whole sequence   */
	/* ######################### */
	
	/* initialize pPrev, pNext */
	pPrev = NULL;
	pNext = pSeqMtf->pMotifList;
	
	/* initialize context */
	nContextNum = 0;
	pContext = NULL;
	pContext = FLEXMTFSCREATE();
	pHead = NULL;
	pTail = NULL;
	pTest = pSeqMtf->pMotifList;
	while(pTest != NULL)
	{
		if(pTest->nStartPos >= nHalfWin)
		{
			break;
		}

		if(pHead == NULL)
			pHead = pTest;
		pTail = pTest;
		nContextNum++;
		pTest = pTest->pNext;
	}

	pContext->pContextPos = NULL;
	pContext->pContextMotif = NULL;
	pContext->pContextStrand = NULL;
	if(nContextNum > 0)
	{
		pContext->pContextPos = CreateIntMatrix(1, nContextNum);
		if(pContext->pContextPos == NULL)
		{
			printf("Error: FlexModule_SampleBA, cannot create the context matrix!\n");
			exit(EXIT_FAILURE);
		}
		pContext->pContextMotif = CreateIntMatrix(1, nContextNum);
		if(pContext->pContextMotif == NULL)
		{
			printf("Error: FlexModule_SampleBA, cannot create the context matrix!\n");
			exit(EXIT_FAILURE);
		}
		pContext->pContextStrand = CreateIntMatrix(1, nContextNum);
		if(pContext->pContextStrand == NULL)
		{
			printf("Error: FlexModule_SampleBA, cannot create the context matrix!\n");
			exit(EXIT_FAILURE);
		}
		
		/* get initial context */
		ni = 0;
		pTest = pHead;
		while(pTest != NULL)
		{
			if(pTest->nStartPos >= nHalfWin)
				break;

			pContext->pContextPos->pMatElement[ni] = pTest->nStartPos+1;
			pContext->pContextMotif->pMatElement[ni] = pTest->nMotifType;
			pContext->pContextStrand->pMatElement[ni] = pTest->nStrand;

			ni++;
			pTest = pTest->pNext;
		}

		if(ni != nContextNum)
		{
			printf("Error: FlexModule_SampleBA, context not synchronized!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* sample base by base */
	BMSETAT(pFreqValid, 0, 0, 1);
	for(ni=0; ni<nLen; ni++)
	{
		/* clear old motif status, and get updated context */
		FlexModule_ClearOldStatus(pSeqMtf, nIndex, ni, 
			pContext, &nContextNum, nBaseTypeNum,
			&pPrev, &pNext, 
			pFreqCount, nHalfWin, 
			nMotifNum, vMotif, 
			pMotifLen, nMaxMotifLen);


		/* if repeat or 'N', skip it */
		if(pBase[ni] >= nBaseTypeNum)
		{
			continue;
		}

		/* get motif likelihood and set posterior probability */
		dLike = log(DMGETAT(pFreqCount, 0, 0));
		DMSETAT(pFreqPost, 0, 0, dLike);
		
		for(nj=0; nj<nMotifNum; nj++)
		{
			nStateId = nj+1;
	
			/* '+' '-' strand */
			for(nk=0; nk<2; nk++)
			{
				/* compute context score & probabiity */
				if(nk == 0)
				{
					dLike =	FlexBaseLogLikelihood_Motif(pSeqMtf, nIndex, ni, 
						vMotif[nj]->pPWM, '+', nMaxMotifLen,
						nBGOrder, pLogBG, pLogBG0, 
						nUseCS, nCSLikeNum, vCSLike,
						&nMask);
				}
				else
				{
					dLike =	FlexBaseLogLikelihood_Motif(pSeqMtf, nIndex, ni, 
						vMotif[nj]->pPWM, '-', nMaxMotifLen,
						nBGOrder, pLogBG, pLogBG0, 
						nUseCS, nCSLikeNum, vCSLike,
						&nMask);
				}

				if(nMask == 1)
				{
					DMSETAT(pFreqPost, nStateId, nk, 0.0);
					BMSETAT(pFreqValid, nStateId, nk, 0);
				}
				else
				{
					if(nContextNum == 0)
					{
						dScore = FlexModule_ComputeModuleScore(nContextNum, NULL,
							NULL, NULL,
							1, 0, nStateId, nk, 
							nMotifNum);
					}
					else
					{
						dScore = FlexModule_ComputeModuleScore(nContextNum, pContext->pContextPos->pMatElement,
							pContext->pContextMotif->pMatElement, pContext->pContextStrand->pMatElement,
							1, 0, nStateId, nk,
							nMotifNum);
					}

					dPAi = FlexModule_ComputeTrueProb(dScore, dModuleD, dModuleL, dModuleA, dT);
					dLike = log((1.0-dPAi)+dPAi*exp(dLike));

					dContextLike = FlexModule_ComputeContextLogLike(pSeqMtf, pContext, ni, 
						nContextNum, 0, nStateId, nk, nMotifNum, dModuleD, dModuleL, dModuleA, dT);
					
					dLike += log(DMGETAT(pFreqCount, nStateId, nk))+dContextLike;

					DMSETAT(pFreqPost, nStateId, nk, dLike);
					BMSETAT(pFreqValid, nStateId, nk, 1);
				}
			}
		}

		 
		/* normalize posterior probability */
		FlexModule_NormalizePost(pFreqPost, pFreqValid);
		
		/* sample motifs */
		FlexModule_SamplePost(pFreqPost, pFreqValid, &nStateId, &nStrand);

		/* update motif frequency count */
		dTemp = DMGETAT(pFreqCount, nStateId, nStrand)+1.0;
		DMSETAT(pFreqCount, nStateId, nStrand, dTemp);

		/* update the new context */
		if(nStateId > 0)
		{
			nMotifId = nStateId-1;

			/* sample Ai */
			if(nStrand == 0)
			{
				dLike =	FlexBaseLogLikelihood_Motif(pSeqMtf, nIndex, ni, 
					vMotif[nMotifId]->pPWM, '+', nMaxMotifLen,
					nBGOrder, pLogBG, pLogBG0, 
					nUseCS, nCSLikeNum, vCSLike,
					&nMask);
			}
			else
			{
				dLike =	FlexBaseLogLikelihood_Motif(pSeqMtf, nIndex, ni, 
					vMotif[nMotifId]->pPWM, '-', nMaxMotifLen,
					nBGOrder, pLogBG, pLogBG0, 
					nUseCS, nCSLikeNum, vCSLike,
					&nMask);
			}

			if(nContextNum == 0)
			{
				dScore = FlexModule_ComputeModuleScore(nContextNum, NULL,
					NULL, NULL,
					1, 0, nStateId, nStrand,
					nMotifNum);
			}
			else
			{
				dScore = FlexModule_ComputeModuleScore(nContextNum, pContext->pContextPos->pMatElement,
					pContext->pContextMotif->pMatElement, pContext->pContextStrand->pMatElement,
					1, 0, nStateId, nStrand,
					nMotifNum);
			}
			
			dPAi = FlexModule_ComputeTrueProb(dScore, dModuleD, dModuleL, dModuleA, dT);
			dLike += (log(dPAi)-log(1.0-dPAi));
			dLike = exp(dLike);
			dLike = dLike/(1.0+dLike);

			dRand = rand_u();
			if(dRand <= dLike)
			{
				pAi[ni] = 1;
			}
			else
			{
				pAi[ni] = 0;
			}

			/* update PWM count */
			if(pAi[ni] == 1)
			{
				if(nStrand == 1)
				{
					FlexModule_AddMotifMatrixCount(vMotif[nMotifId]->pSampleCount, 
							  pBase+ni, pMotifLen->pMatElement[nMotifId], '-');
				}
				else
				{
					FlexModule_AddMotifMatrixCount(vMotif[nMotifId]->pSampleCount, 
							  pBase+ni, pMotifLen->pMatElement[nMotifId], '+');
				}

				FLEXMTFMREFRESH(vMotif[nMotifId]);

				if(nIsRecording == 2)
				{
					pSeqMtf->vMonitor[nStateId]->pMatElement[ni] += 1.0;
					pMotifNSample->pMatElement[nMotifId] += 1.0;
				}
			}

			/* add new Bi site */
			pSeqMtf->nSiteNum += 1;
			FlexModule_UpdateNewStatus(pSeqMtf, nIndex, ni,
				nStateId, nStrand,
				pContext, &nContextNum, nBaseTypeNum,
				&pPrev, &pNext, 
				nHalfWin, 
				nMotifNum, vMotif, 
				pMotifLen, nMaxMotifLen);

		}
	}

	/* release memory */
	DestroyDoubleMatrix(pFreqPost);
	DestroyByteMatrix(pFreqValid);

	if(nContextNum == 0)
	{
		if(pContext->pContextPos != NULL)
		{
			printf("Error: FlexModule_SampleBA, context not synchronized!\n");
			exit(EXIT_FAILURE);
		}
	}
	else
	{
		if(pContext->pContextPos->nWidth != nContextNum)
		{
			printf("Error: FlexModule_SampleBA, context not synchronized!\n");
			exit(EXIT_FAILURE);
		}
	}

	FLEXMTFSDESTROY(pContext);

	/* return */
	return PROC_SUCCESS;
}



/* ----------------------------------------------------------------------- */ 
/*  FlexModule_SampleL: sample module length.                              */
/* ----------------------------------------------------------------------- */ 
int FlexModule_SampleL(struct FLEXSEQMOTIF **vSeqMtf, int nSeqCount, 
				int nIndex, double dT, 
				double dModuleD, double dModuleL, double dModuleA, 
				int *pModuleLen, int nMinModuleLen, int nMaxModuleLen, 
				int nMotifNum, struct FLEXMOTIFMATRIX **vMotif, 
				struct INTMATRIX *pMotifLen, int nMaxMotifLen, 
				int nIsRecording)
{
	/* define */
	int nNewModuleLen;
	double dLike0,dLike1,dR,dRand;
	int ni,nj,nk;
	struct FLEXMOTIFSITE *pSite;
	unsigned char *pAi;
	struct FLEXMOTIFSITE *pNewSite,*pPrevSite;
	struct FLEXMOTIFSITE *pHead,*pTail,*pTest;
	int nContextNum,nHalfWin;
	double dScore,dPAi0,dPAi1;
	struct INTMATRIX *pTempPos,*pTempMotif,*pTempStrand;
	
	/* init trial */
	dRand = rand_u();
	if(dRand < 0.5)
		nNewModuleLen = *pModuleLen-20;
	else if(dRand > 0.5)
		nNewModuleLen = *pModuleLen+20;

	if( (nNewModuleLen < nMinModuleLen) || (nNewModuleLen > nMaxModuleLen) )
	{
		/* return */
		return PROC_SUCCESS;
	}

	/* get likelihood */
	nHalfWin = (int)(nNewModuleLen/2);
	dLike0 = 0.0;
	dLike1 = 0.0;

	/* compute likelihood */
	for(ni=0; ni<nSeqCount; ni++)
	{
		pAi = vSeqMtf[ni]->vStatus[0]->pMatElement;
		pSite = vSeqMtf[ni]->pMotifList;
		while(pSite != NULL)
		{
			/* create new site */
			pNewSite = NULL;
			pNewSite = FLEXMTFSCREATE();
			if(pNewSite == NULL)
			{
				printf("Error: FlexModule_SampleL, cannot create trial site!\n");
				exit(EXIT_FAILURE);
			}
			pNewSite->pNext = NULL;
			pNewSite->pPrev = NULL;

			/* get context */
			/* establish context limit */
			nContextNum = 1;
			pHead = pSite;
			pTail = pSite;

			pTest = pHead->pPrev;
			while(pTest != NULL)
			{
				if( (pSite->nStartPos-pTest->nStartPos) > nHalfWin )
					break;

				pHead = pTest;
				pTest = pHead->pPrev;
				nContextNum++;
			}

			pTest = pTail->pNext;
			while(pTest != NULL)
			{
				if( (pTest->nStartPos-pSite->nStartPos) > nHalfWin )
					break;

				pTail = pTest;
				pTest = pTail->pNext;
				nContextNum++;
			}

			/* create context matrix */
			pNewSite->pContextPos = NULL;
			pNewSite->pContextPos = CreateIntMatrix(1, nContextNum);
			if(pNewSite->pContextPos == NULL)
			{
				printf("Error: FlexModule_SampleL, cannot create context matrix correctly!\n");
				exit(EXIT_FAILURE);
			}
			pNewSite->pContextMotif = NULL;
			pNewSite->pContextMotif = CreateIntMatrix(1, nContextNum);
			if(pNewSite->pContextMotif == NULL)
			{
				printf("Error: FlexModule_SampleL, cannot create context matrix correctly!\n");
				exit(EXIT_FAILURE);
			}
			pNewSite->pContextStrand = NULL;
			pNewSite->pContextStrand = CreateIntMatrix(1, nContextNum);
			if(pNewSite->pContextStrand == NULL)
			{
				printf("Error: FlexModule_SampleL, cannot create context matrix correctly!\n");
				exit(EXIT_FAILURE);
			}

			/* fill out the context */
			pTest = pHead;
			nk = 0;
			while(pTest != NULL)
			{
				if(pTest->nStartPos - pTail->nStartPos > 0)
					break;

				pNewSite->pContextPos->pMatElement[nk] = pTest->nStartPos-pSite->nStartPos;
				pNewSite->pContextMotif->pMatElement[nk] = pTest->nMotifType;
				pNewSite->pContextStrand->pMatElement[nk] = pTest->nStrand;
				nk++;
				pTest = pTest->pNext;
			}

			if(nk != nContextNum)
			{
				printf("Error: FlexModule_SampleL, number of motif sites in the context do not match!\n");
				exit(EXIT_FAILURE);
			}

			/* add to list */
			if(vSeqMtf[ni]->pTrialMotifList == NULL)
			{
				vSeqMtf[ni]->pTrialMotifList = pNewSite;
				pPrevSite = pNewSite;
			}
			else
			{
				pPrevSite->pNext = pNewSite;
				pNewSite->pPrev = pPrevSite;
				pPrevSite = pNewSite;
			}

			/* compute likelihood */
			dScore = FlexModule_ComputeModuleScore(pSite->pContextPos->nWidth, pSite->pContextPos->pMatElement,
				pSite->pContextMotif->pMatElement, pSite->pContextStrand->pMatElement, 
				0, 0, 0, 0, nMotifNum);
			dPAi0 = FlexModule_ComputeTrueProb(dScore, dModuleD, dModuleL, dModuleA, dT);
			
			dScore = FlexModule_ComputeModuleScore(nContextNum, pNewSite->pContextPos->pMatElement,
				pNewSite->pContextMotif->pMatElement, pNewSite->pContextStrand->pMatElement, 
				0, 0, 0, 0, nMotifNum);
			dPAi1 = FlexModule_ComputeTrueProb(dScore, dModuleD, dModuleL, dModuleA, dT);

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

	/* sample */
	dR = exp(dLike1-dLike0);

	dRand = rand_u();
	/* update new context and length */
	if(dRand < dR)
	{
		/* modify module length */
		*pModuleLen = nNewModuleLen;
		
		/* change context */
		for(ni=0; ni<nSeqCount; ni++)
		{
			pSite = vSeqMtf[ni]->pMotifList;
			pNewSite = vSeqMtf[ni]->pTrialMotifList;
			while(pSite != NULL)
			{
				/* update context */
				pTempPos = pSite->pContextPos;
				pTempMotif = pSite->pContextMotif;
				pTempStrand = pSite->pContextStrand;
				pSite->pContextPos = pNewSite->pContextPos;
				pSite->pContextMotif = pNewSite->pContextMotif;
				pSite->pContextStrand = pNewSite->pContextStrand;
				pNewSite->pContextPos = pTempPos;
				pNewSite->pContextMotif = pTempMotif;
				pNewSite->pContextStrand = pTempStrand;

				/* get next */
				pSite = pSite->pNext;
				vSeqMtf[ni]->pTrialMotifList = pNewSite->pNext;
				FLEXMTFSDESTROY(pNewSite);
				pNewSite = vSeqMtf[ni]->pTrialMotifList;
			}

			vSeqMtf[ni]->pTrialMotifList = NULL;
		}
	}
	else
	{
		for(ni=0; ni<nSeqCount; ni++)
		{
			FLEXMTFSDESTROYLIST(&(vSeqMtf[ni]->pTrialMotifList));
			vSeqMtf[ni]->pTrialMotifList = NULL;
		}
	}

	/* return */
	return PROC_SUCCESS;
}

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
				int nIsRecording)
{
	/* define */
	/* for data access */
	unsigned char *pBase,*pAi,*pSi,*pCS;
	double *pBL;
	struct FLEXMOTIFSITE **vSite;
	struct FLEXMOTIFSITE *pSite;
	int ni,nj,nPos,nNewPos,nBaseId;
	
	/* new context */
	struct INTMATRIX *pTestPos;
	struct INTMATRIX *pTestMotif;
	struct INTMATRIX *pTestStrand;
	struct BYTEMATRIX *pTestA;
	double dScore, dRatio, dPAi, dTemp;
	
	/* motif parameters */
	int nMId,nLen;
	int nMotifLen;
	
	/* for local move */
	double dProposalProb[3];
	double dRand;
	int nStep;

	/* for likelihood calculation */
	struct DOUBLEMATRIX *pMCount;
	int nMask;
	double dLike[2];
	int nCSId;
	int nCSShift;
	int nBaseTypeNum;

	/* choose local move */
	for(ni=0; ni<3; ni++)
	{
		dProposalProb[ni] = (double)(ni+1)/3.0;
	}
	dRand = rand_u();
	for(ni=0; ni<3; ni++)
	{
		if(dRand <= dProposalProb[ni])
		{
			break;
		}
	}
	if(ni == 1)
	{
		nStep = 1;
	}
	else if(ni==2)
	{
		nStep = -1;
	}
	else
	{
		/* return */
		return PROC_SUCCESS;
	}

	/* if need to check move, do initialization first */
	nBaseTypeNum = pLogBG->nWidth;
	nMId = nMotifId-1;
	nMotifLen = pMotifLen->pMatElement[nMId];
	nCSShift = (int)((nMotifLen-1)/2);

	/* vectors below are arranged in the order of 0, 1/-1 */
	for(ni=0; ni<2; ni++)
	{
		dLike[ni] = 0.0;
	}

	pMCount = NULL;
	pMCount = CreateDoubleMatrix(1, nBaseTypeNum);
	if(pMCount == NULL)
	{
		printf("Error: FlexModule_LocalShift, cannot prepare space for likelihood computation!\n");
		exit(EXIT_FAILURE);
	}


	/* scan sequence by sequence */
	nMask = 0;
	for(ni=0; ni<nSeqCount; ni++)
	{
		/* if there are no sites, jump to the next sequence */
		if(vSeqMtf[ni]->nSiteNum == 0)
		{
			if(vSeqMtf[ni]->pMotifList != NULL)
			{
				printf("Error: FlexModule_LocalShift, motif site counting error!\n");
				exit(EXIT_FAILURE);
			}
			continue;
		}

		/* if already masked, break */
		if(nMask == 1)
		{
			break;
		}

		/* sequence wise init */
		nLen = vSeqMtf[ni]->vSeq[nIndex]->nWidth;
		pBase = vSeqMtf[ni]->vSeq[nIndex]->pMatElement;
		if(nUseCS == 1)
		{
			pCS = vSeqMtf[ni]->vScore[0]->pMatElement;
		}
		pAi = vSeqMtf[ni]->vStatus[0]->pMatElement;
		pSi = vSeqMtf[ni]->vStatus[1]->pMatElement;
		vSite = vSeqMtf[ni]->vMotif;
		pBL = vSeqMtf[ni]->vMonitor[0]->pMatElement;

		/* init whole sequence context */		
		pTestPos = NULL;
		pTestPos = CreateIntMatrix(1, vSeqMtf[ni]->nSiteNum);
		pTestMotif = NULL;
		pTestMotif = CreateIntMatrix(1, vSeqMtf[ni]->nSiteNum);
		pTestStrand = NULL;
		pTestStrand = CreateIntMatrix(1, vSeqMtf[ni]->nSiteNum);
		pTestA = NULL;
		pTestA = CreateByteMatrix(1, vSeqMtf[ni]->nSiteNum);

		/* scan candidate site one by one */
		pSite = vSeqMtf[ni]->pMotifList;
		nj = 0;
		while(pSite != NULL)
		{
			/* store sites */
			nPos = pSite->nStartPos;
			pTestPos->pMatElement[nj] = nPos;
			pTestMotif->pMatElement[nj] = pSite->nMotifType;
			pTestStrand->pMatElement[nj] = pSite->nStrand;
			pTestA->pMatElement[nj] = pAi[nPos];

			/* get context score */
			dScore = FlexModule_ComputeModuleScore(pSite->pContextPos->nWidth,
				pSite->pContextPos->pMatElement, pSite->pContextMotif->pMatElement,
				pSite->pContextStrand->pMatElement, 0, 0, 0, 0,
				pMotifLen->nWidth);
			dPAi = FlexModule_ComputeTrueProb(dScore, dModuleD, dModuleL, dModuleA, dT);
			if(pAi[nPos] == 1)
			{
				dLike[0] += log(dPAi);
			}
			else
			{
				dLike[0] += log(1.0-dPAi);
			}

			/* if the motif site in question */
			if( pSite->nMotifType == nMotifId ) 
			{
				/* count bases */
				/* '-' strand */
				if(pSite->nStrand == 1)
				{
					/* +1 */
					if(nStep == 1)
					{
						/* move +nk */
						pTestPos->pMatElement[nj] = nPos-nStep;
						nNewPos = nPos-nStep;
						if( nNewPos < 0 )
						{
							nMask = 1;
							break;
						}

						if(pBase[nNewPos] >= nBaseTypeNum)
						{
							nMask = 1;
							break;
						}

						if(pAi[nPos] == 1)
						{
							dTemp = pBL[nNewPos+nMotifLen]-pBL[nNewPos];
							dLike[1] += dTemp;
							
							if(nUseCS == 1)
							{
								nCSId = pCS[nNewPos+nCSShift];
								dLike[1] += (vCSLike[1]->pMatElement[nCSId]-vCSLike[0]->pMatElement[nCSId]);
							}

							/* count bases */
							nBaseId = nBaseTypeNum-1-(int)pBase[nNewPos];
							pMCount->pMatElement[nBaseId] += 1.0;
						}
					}
					/* -1 */
					else if(nStep == -1)
					{
						/* move -nk */
						pTestPos->pMatElement[nj] = nPos-nStep;
						nNewPos = nPos-nStep+nMotifLen-1;
						if( nNewPos >= nLen )
						{
							nMask = 1;
							break;
						}

						if(pBase[nNewPos] >= nBaseTypeNum)
						{
							nMask = 1;
							break;
						}
						
						if(pAi[nPos] == 1)
						{
							dTemp = pBL[nNewPos-nMotifLen]-pBL[nNewPos];
							dLike[1] += dTemp;
							
							if(nUseCS == 1)
							{
								nCSId = pCS[nPos-nStep+nCSShift];
								dLike[1] += (vCSLike[1]->pMatElement[nCSId]-vCSLike[0]->pMatElement[nCSId]);
							}

							/* count bases */
							nBaseId = nBaseTypeNum-1-(int)pBase[nNewPos];
							pMCount->pMatElement[nBaseId] += 1.0;
						}
					}
					/* otherwise */
					else
					{
						printf("Error: FlexModule_LocalShift, sample strand wrong!\n");
						exit(EXIT_FAILURE);
					}
				}
				
				/* '+' strand */
				else
				{
					/* +1 */
					if(nStep == 1)
					{
						/* move +nk */
						pTestPos->pMatElement[nj] = nPos+nStep;
						nNewPos = nPos+nStep+nMotifLen-1;
						if( nNewPos >= nLen )
						{
							nMask = 1;
							break;
						}

						if(pBase[nNewPos] >= nBaseTypeNum)
						{
							nMask = 1;
							break;
						}
						
						if(pAi[nPos] == 1)
						{
							dTemp = pBL[nNewPos-nMotifLen]-pBL[nNewPos];
							dLike[1] += dTemp;
							
							if(nUseCS == 1)
							{
								nCSId = pCS[nPos+nStep+nCSShift];
								dLike[1] += (vCSLike[1]->pMatElement[nCSId]-vCSLike[0]->pMatElement[nCSId]);
							}

							/* count bases */
							nBaseId = (int)pBase[nNewPos];
							pMCount->pMatElement[nBaseId] += 1.0;
						}
					}
					/* -1 */
					else if(nStep == -1)
					{
						/* move -nk */
						pTestPos->pMatElement[nj] = nPos+nStep;
						nNewPos = nPos+nStep;
						if( nNewPos < 0 )
						{
							nMask = 1;
							break;
						}

						if(pBase[nNewPos] >= nBaseTypeNum)
						{
							nMask = 1;
							break;
						}

						if(pAi[nPos] == 1)
						{
							dTemp = pBL[nNewPos+nMotifLen]-pBL[nNewPos];
							dLike[1] += dTemp;
							
							if(nUseCS == 1)
							{
								nCSId = pCS[nNewPos+nCSShift];
								dLike[1] += (vCSLike[1]->pMatElement[nCSId]-vCSLike[0]->pMatElement[nCSId]);
							}

							/* count bases */
							nBaseId = (int)pBase[nNewPos];
							pMCount->pMatElement[nBaseId] += 1.0;
						}
					}
					/* otherwise */
					else
					{
						printf("Error: FlexModule_LocalShift, sample strand wrong!\n");
						exit(EXIT_FAILURE);
					}
				}

				/* change dLike0 */
				if( (nUseCS == 1) && (pAi[nPos] == 1) )
				{
					nCSId = pCS[nPos+nCSShift];
					dLike[0] += (vCSLike[1]->pMatElement[nCSId]-vCSLike[0]->pMatElement[nCSId]);
				}
			}

			/* get next */
			nj++;
			pSite = pSite->pNext;
		}

		/* if already masked, stop scan */
		if(nMask == 1)
		{
			/* destroy whole sequence context */
			DestroyIntMatrix(pTestPos);
			DestroyIntMatrix(pTestMotif);
			DestroyIntMatrix(pTestStrand);
			DestroyByteMatrix(pTestA);
			break;
		}

		if(nj != vSeqMtf[ni]->nSiteNum)
		{
			printf("Error: FlexModule_LocalShift, motif site counting error!\n");
			exit(EXIT_FAILURE);
		}


		/* get new context likelihood */
		dScore = FlexModule_ShiftContextLogLike(pTestPos, pTestMotif, pTestStrand, pTestA,
			nLen, pMotifLen, nModuleLen, dModuleD, dModuleL, dModuleA, dT, &nMask);

		/* if already masked, stop scan */
		if(nMask == 1)
		{
			/* destroy whole sequence context */
			DestroyIntMatrix(pTestPos);
			DestroyIntMatrix(pTestMotif);
			DestroyIntMatrix(pTestStrand);
			DestroyByteMatrix(pTestA);
			break;
		}

		dLike[1] += dScore;
			
		/* destroy whole sequence context */
		DestroyIntMatrix(pTestPos);
		DestroyIntMatrix(pTestMotif);
		DestroyIntMatrix(pTestStrand);
		DestroyByteMatrix(pTestA);
	}

	/* if masked, return */
	if(nMask == 1)
	{
		DestroyDoubleMatrix(pMCount);
		return PROC_SUCCESS;
	}

	/* update the whole count likelihood */
	dScore = FlexModule_ShiftMotifLogLike(vMotif[nMId], pMCount, nStep);
	dLike[1] += dScore;
	dScore = FlexModule_ShiftMotifLogLike(vMotif[nMId], pMCount, 0);
	dLike[0] += dScore;

	/* sample shift */
	dRatio = exp(dLike[1]-dLike[0]);
	dRand = rand_u();

	/* update status */
	if(dRand <= dRatio)
	{
		/* update sequence by sequence */
		for(ni=0; ni<nSeqCount; ni++)
		{
			FlexModule_ShiftUpdate(vSeqMtf[ni], nIndex, nMotifId, nStep,
				nModuleLen, pMotifLen, nMaxMotifLen);
		}

		/* update PWM */
		FlexModule_ShiftUpdateMotif(vMotif[nMId], pMCount, nStep);
	}

	/* release memory */
	DestroyDoubleMatrix(pMCount);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_ShiftUpdate: update all indicators and context after shift  */
/* ----------------------------------------------------------------------- */ 
int FlexModule_ShiftUpdate(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, 
				int nMotifId, int nStep,
				int nModuleLen, struct INTMATRIX *pMotifLen, int nMaxMotifLen)
{
	/* define */
	unsigned char *pBase,*pAi,*pSi;
	struct FLEXMOTIFSITE **vSite;
	struct FLEXMOTIFSITE *pSite;
	int nLen,nMId,nMotifLen;
	int ni,nSiteNum,nOldPos,nNewPos,nContextNum;
	int nHalfWin;
	struct FLEXMOTIFSITE *pHead,*pTail,*pTest;
	

	/* init */
	nHalfWin = (int)(nModuleLen/2);
	nMId = nMotifId-1;
	nMotifLen = pMotifLen->pMatElement[nMId];
	nLen = pSeqMtf->vSeq[nIndex]->nWidth;
	pBase = pSeqMtf->vSeq[nIndex]->pMatElement;
	pAi = pSeqMtf->vStatus[0]->pMatElement;
	pSi = pSeqMtf->vStatus[1]->pMatElement;
	vSite = pSeqMtf->vMotif;

	/* scan all sites */
	pSite = pSeqMtf->pMotifList;
	nSiteNum = 0;
	while(pSite != NULL)
	{
		/* if motif in question, update */
		if(pSite->nMotifType == nMotifId)
		{
			/* if - strand */
			if(pSite->nStrand == 1)
			{
				nOldPos = pSite->nStartPos;
				nNewPos = nOldPos-nStep;
				
				/* update Ai */
				pAi[nNewPos] = pAi[nOldPos];
				pAi[nOldPos] = 0;

				/* update Si */
				if(nStep > 0)
				{
					for(ni=0; ni<nStep; ni++)
					{
						pSi[nNewPos-ni] = nMotifId;
						pSi[nNewPos+nMotifLen-ni] = 0; 
					}
				}
				else
				{
					for(ni=0; ni>nStep; ni--)
					{
						pSi[nOldPos+nMotifLen-ni] = nMotifId;
						pSi[nOldPos-ni] = 0; 
					}
				}

				/* update vSite */
				vSite[nNewPos] = vSite[nOldPos];
				vSite[nOldPos] = NULL;

				/* update site itself */
				pSite->nStartPos = nNewPos;
			}
			/* if + strand */
			else
			{
				nOldPos = pSite->nStartPos;
				nNewPos = nOldPos+nStep;
				
				/* update Ai */
				pAi[nNewPos] = pAi[nOldPos];
				pAi[nOldPos] = 0;

				/* update Si */
				if(nStep > 0)
				{
					for(ni=0; ni<nStep; ni++)
					{
						pSi[nOldPos+nMotifLen+ni] = nMotifId;
						pSi[nOldPos+ni] = 0; 
					}
				}
				else
				{
					for(ni=0; ni>nStep; ni--)
					{
						pSi[nNewPos+ni] = nMotifId;
						pSi[nNewPos+nMotifLen+ni] = 0; 
					}
				}

				/* update vSite */
				vSite[nNewPos] = vSite[nOldPos];
				vSite[nOldPos] = NULL;

				/* update site itself */
				pSite->nStartPos = nNewPos;
			}
		}

		/* get next */
		nSiteNum++;
		pSite = pSite->pNext;
	}

	if(nSiteNum != pSeqMtf->nSiteNum)
	{
		printf("Error: FlexModule_ShiftUpdate, site number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* update context */
	pSite = pSeqMtf->pMotifList;
	while(pSite != NULL)
	{
		/* establish context limit */
		nContextNum = 1;
		pHead = pSite;
		pTail = pSite;

		pTest = pHead->pPrev;
		while(pTest != NULL)
		{
			if( (pSite->nStartPos-pTest->nStartPos) > nHalfWin )
				break;

			pHead = pTest;
			pTest = pHead->pPrev;
			nContextNum++;
		}

		pTest = pTail->pNext;
		while(pTest != NULL)
		{
			if( (pTest->nStartPos-pSite->nStartPos) > nHalfWin )
				break;

			pTail = pTest;
			pTest = pTail->pNext;
			nContextNum++;
		}

		/* create context matrix */
		DestroyIntMatrix(pSite->pContextPos);
		pSite->pContextPos = NULL;
		pSite->pContextPos = CreateIntMatrix(1, nContextNum);
		if(pSite->pContextPos == NULL)
		{
			printf("Error: FlexModule_ShiftUpdate, cannot create context matrix correctly!\n");
			exit(EXIT_FAILURE);
		}
		DestroyIntMatrix(pSite->pContextMotif);
		pSite->pContextMotif = NULL;
		pSite->pContextMotif = CreateIntMatrix(1, nContextNum);
		if(pSite->pContextMotif == NULL)
		{
			printf("Error: FlexModule_ShiftUpdate, cannot create context matrix correctly!\n");
			exit(EXIT_FAILURE);
		}
		DestroyIntMatrix(pSite->pContextStrand);
		pSite->pContextStrand = NULL;
		pSite->pContextStrand = CreateIntMatrix(1, nContextNum);
		if(pSite->pContextStrand == NULL)
		{
			printf("Error: FlexModule_ShiftUpdate, cannot create context matrix correctly!\n");
			exit(EXIT_FAILURE);
		}

		/* fill out the context */
		pTest = pHead;
		ni = 0;
		while(pTest != NULL)
		{
			if(pTest->nStartPos - pTail->nStartPos > 0)
				break;

			pSite->pContextPos->pMatElement[ni] = pTest->nStartPos-pSite->nStartPos;
			pSite->pContextMotif->pMatElement[ni] = pTest->nMotifType;
			pSite->pContextStrand->pMatElement[ni] = pTest->nStrand;
			ni++;
			pTest = pTest->pNext;
		}

		if(ni != nContextNum)
		{
			printf("Error: FlexModule_ShiftUpdate, number of motif sites in the context do not match!\n");
			exit(EXIT_FAILURE);
		}

		/* get next */
		pSite = pSite->pNext;
	}
	

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_ClearOldStatus: clear old status of a given position.       */
/* ----------------------------------------------------------------------- */ 
double FlexModule_ShiftContextLogLike(struct INTMATRIX *pPos, struct INTMATRIX *pMotif, 
				struct INTMATRIX *pStrand, struct BYTEMATRIX *pA,
				int nLen, struct INTMATRIX *pMotifLen, 
				int nModuleLen, double dModuleD, double dModuleL, double dModuleA,
				double dT, int *pMask)
{
	/* define */
	double dL = 0.0;
	struct INTMATRIX *pLocalPos;
	int nHead,nTail;
	int ni,nj,nk;
	int nSiteNum;
	int nContextNum;
	int nHalfWin;
	double dScore,dPAi;

	/* init */
	*pMask = 0;
	nHalfWin = (int)(nModuleLen/2);
	nSiteNum = pPos->nWidth;

	if(nSiteNum <= 0)
	{
		printf("Error: FlexModule_ShiftContextLogLike, no sites available!\n");
		exit(EXIT_FAILURE);
	}

	/* check boundary */
	if(pPos->pMatElement[0] < 0)
	{
		*pMask = 1;
		dL = 0.0;
		return dL;
	}

	ni = nSiteNum-1;
	nj = pMotif->pMatElement[ni]-1;
	if( (pPos->pMatElement[ni]+pMotifLen->pMatElement[nj]) > nLen )
	{
		*pMask = 1;
		dL = 0.0;
		return dL;
	}

	pLocalPos = NULL;
	pLocalPos = IMCLONE(pPos);
	if(pLocalPos == NULL)
	{
		printf("Error: FlexModule_ShiftContextLogLike, cannot compute the context log-likelihood!\n");
		exit(EXIT_FAILURE);
	}

	/* start */
	nHead = 0;
	nTail = 0;
	nContextNum = 1;
	for(ni=0; ni<nSiteNum; ni++)
	{
		/* check motif overlap */
		if(ni > 0)
		{
			nj = pMotif->pMatElement[ni-1]-1;
			if( (pPos->pMatElement[ni]-pPos->pMatElement[ni-1]) < pMotifLen->pMatElement[nj] )
			{
				*pMask = 1;
				dL = 0.0;
				DestroyIntMatrix(pLocalPos);
				return dL;
			}
		}
		
		/* get context */
		while(nHead < nSiteNum)
		{
			if( (pPos->pMatElement[nHead] - pPos->pMatElement[ni]) >= -nHalfWin )
			{
				break;
			}
			nHead++;
		}
		while(nTail < nSiteNum)
		{
			if( (pPos->pMatElement[nTail] - pPos->pMatElement[ni]) > nHalfWin )
			{
				break;
			}
			nTail++;
		}

		nContextNum = nTail-nHead;
		if(nContextNum <= 0)
		{
			printf("Error: FlexModule_ShiftContextLogLike, context not correctly synchronized!\n");
			exit(EXIT_FAILURE);
		}

		for(nk=nHead; nk<nTail; nk++)
		{
			pLocalPos->pMatElement[nk] = pPos->pMatElement[nk]-pPos->pMatElement[ni];
		}

		/* get context score */
		dScore = FlexModule_ComputeModuleScore(nContextNum, pLocalPos->pMatElement+nHead, 
			pMotif->pMatElement+nHead, pStrand->pMatElement+nHead, 
			0, 0, 0, 0, pMotifLen->nWidth);
		dPAi = FlexModule_ComputeTrueProb(dScore, dModuleD, dModuleL, dModuleA, dT);
		
		if(pA->pMatElement[ni] == 1)
		{
			dL += log(dPAi);
		}
		else
		{
			dL += log(1.0-dPAi);
		}
	}


	/* destroy */
	DestroyIntMatrix(pLocalPos);

	/* return */
	return dL;
}

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_ShiftMotifLogLike: compute shifted motif loglikelihood      */
/* ----------------------------------------------------------------------- */ 
double FlexModule_ShiftMotifLogLike(struct FLEXMOTIFMATRIX *pMotif, 
							 struct DOUBLEMATRIX *pMCount, int nStep)
{
	/* define */
	double dL = 0.0;
	double dTotal,dTemp;
	double *pEle1,*pEle2;
	int ni,nj,nMotifLen,nBaseTypeNum;

	/* get loglikelihood */
	nMotifLen = pMotif->pPWM->nHeight;
	nBaseTypeNum = pMotif->pPWM->nWidth;
	if(nStep == 1)
	{
		
		pEle1 = pMotif->pSampleCount->pMatElement+nBaseTypeNum;
		pEle2 = pMotif->pPriorCount->pMatElement;
		for(ni=1; ni<nMotifLen; ni++)
		{
			dTotal = 0.0;
			for(nj=0; nj<nBaseTypeNum; nj++)
			{
				dTemp = (*pEle1)+(*pEle2);
				dTotal += dTemp;
				dL += gammaln(dTemp);
				pEle1++;
				pEle2++;
			}
			dL -= gammaln(dTotal);
		}
		
		pEle1 = pMCount->pMatElement;
		dTotal = 0.0;
		for(nj=0; nj<nBaseTypeNum; nj++)
		{
			dTemp = (*pEle1)+(*pEle2);
			dTotal += dTemp;
			dL += gammaln(dTemp);
			pEle1++;
			pEle2++;
		}
		dL -= gammaln(dTotal);
	}
	else if(nStep == -1)
	{
		pEle1 = pMCount->pMatElement;
		pEle2 = pMotif->pPriorCount->pMatElement;
		dTotal = 0.0;
		for(nj=0; nj<nBaseTypeNum; nj++)
		{
			dTemp = (*pEle1)+(*pEle2);
			dTotal += dTemp;
			dL += gammaln(dTemp);
			pEle1++;
			pEle2++;
		}
		dL -= gammaln(dTotal);

		pEle1 = pMotif->pSampleCount->pMatElement;
		for(ni=1; ni<nMotifLen; ni++)
		{
			dTotal = 0.0;
			for(nj=0; nj<nBaseTypeNum; nj++)
			{
				dTemp = (*pEle1)+(*pEle2);
				dTotal += dTemp;
				dL += gammaln(dTemp);
				pEle1++;
				pEle2++;
			}
			dL -= gammaln(dTotal);
		}
	}
	else if(nStep == 0)
	{
		pEle1 = pMotif->pSampleCount->pMatElement;
		pEle2 = pMotif->pPriorCount->pMatElement;
		for(ni=0; ni<nMotifLen; ni++)
		{
			dTotal = 0.0;
			for(nj=0; nj<nBaseTypeNum; nj++)
			{
				dTemp = (*pEle1)+(*pEle2);
				dTotal += dTemp;
				dL += gammaln(dTemp);
				pEle1++;
				pEle2++;
			}
			dL -= gammaln(dTotal);
		}
	}
	else
	{
		printf("Error: FlexModule_ShiftMotifLogLike, cannot get the shifted loglikelihood!\n");
		exit(EXIT_FAILURE);
	}


	/* return */
	return dL;
}

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_ShiftUpdateMotif: update shifted motif matrix               */
/* ----------------------------------------------------------------------- */ 
int FlexModule_ShiftUpdateMotif(struct FLEXMOTIFMATRIX *pMotif, 
							 struct DOUBLEMATRIX *pMCount, int nStep)
{
	/* define */
	double *pEle1,*pEle2;
	int ni,nj,nMotifLen,nBaseTypeNum;

	/* get loglikelihood */
	nMotifLen = pMotif->pPWM->nHeight;
	nBaseTypeNum = pMotif->pPWM->nWidth;
	if(nStep == 1)
	{
		pEle1 = pMotif->pSampleCount->pMatElement;
		pEle2 = pEle1+nBaseTypeNum;
		for(ni=1; ni<nMotifLen; ni++)
		{
			for(nj=0; nj<nBaseTypeNum; nj++)
			{
				*pEle1 = *pEle2;
				pEle1++;
				pEle2++;
			}
		}
		
		pEle2 = pMCount->pMatElement;
		for(nj=0; nj<nBaseTypeNum; nj++)
		{
			*pEle1 = *pEle2;
			pEle1++;
			pEle2++;
		}

		FLEXMTFMREFRESH(pMotif);
	}
	else if(nStep == -1)
	{
		pEle1 = pMotif->pSampleCount->pMatElement+nBaseTypeNum*nMotifLen-1;
		pEle2 = pEle1-nBaseTypeNum;
		for(ni=1; ni<nMotifLen; ni++)
		{
			for(nj=0; nj<nBaseTypeNum; nj++)
			{
				*pEle1 = *pEle2;
				pEle1--;
				pEle2--;
			}
		}
		
		pEle2 = pMCount->pMatElement+nBaseTypeNum-1;
		for(nj=0; nj<nBaseTypeNum; nj++)
		{
			*pEle1 = *pEle2;
			pEle1--;
			pEle2--;
		}

		FLEXMTFMREFRESH(pMotif);

	}
	else
	{
		printf("Error: FlexModule_ShiftUpdateMotif, cannot get the shifted loglikelihood!\n");
		exit(EXIT_FAILURE);
	}

	/* return */
	return PROC_SUCCESS;
}

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
				double dDefaultPriorCount)
{
	/* define */
	/* for local move */
	double dProposalProb[5];
	double dRand;
	int nDirec;
	int ni;

	
	/* propose local move */
	for(ni=0; ni<5; ni++)
	{
		dProposalProb[ni] = (double)(ni+1)/5.0;
	}
	dRand = rand_u();
	for(ni=0; ni<5; ni++)
	{
		if(dRand <= dProposalProb[ni])
		{
			break;
		}
	}

	if(ni == 1)
	{
		nDirec = 1;
		FlexModule_IncreaseW(vSeqMtf, nSeqCount, 
				nIndex, dT, 
				pFreqCount, dModuleD, 
				dModuleL, dModuleA, nModuleLen,
				nBGOrder, pLogBG, pLogBG0,
				nMotifId, nMotifNum, vMotif, 
				pMotifLen, pMaxMotifLen,
				nUseCS, nCSLikeNum, vCSLike,
				nIsRecording, nDirec, nMeanMotifLen, 
				nMotifLenUpperBound, nMotifLenLowerBound,
				dDefaultPriorCount);
	}
	else if(ni==2)
	{
		nDirec = -1;
		FlexModule_IncreaseW(vSeqMtf, nSeqCount, 
				nIndex, dT, 
				pFreqCount, dModuleD, 
				dModuleL, dModuleA, nModuleLen,
				nBGOrder, pLogBG, pLogBG0,
				nMotifId, nMotifNum, vMotif, 
				pMotifLen, pMaxMotifLen,
				nUseCS, nCSLikeNum, vCSLike,
				nIsRecording, nDirec, nMeanMotifLen, 
				nMotifLenUpperBound, nMotifLenLowerBound,
				dDefaultPriorCount);
	}
	else if(ni == 3)
	{
		nDirec = 1;
		FlexModule_DecreaseW(vSeqMtf, nSeqCount, 
				nIndex, dT, 
				pFreqCount, dModuleD, 
				dModuleL, dModuleA, nModuleLen,
				nBGOrder, pLogBG, pLogBG0,
				nMotifId, nMotifNum, vMotif, 
				pMotifLen, pMaxMotifLen,
				nUseCS, nCSLikeNum, vCSLike,
				nIsRecording, nDirec, nMeanMotifLen, 
				nMotifLenUpperBound, nMotifLenLowerBound,
				dDefaultPriorCount);
	}
	else if(ni == 4)
	{
		nDirec = -1;
		FlexModule_DecreaseW(vSeqMtf, nSeqCount, 
				nIndex, dT, 
				pFreqCount, dModuleD, 
				dModuleL, dModuleA, nModuleLen,
				nBGOrder, pLogBG, pLogBG0,
				nMotifId, nMotifNum, vMotif, 
				pMotifLen, pMaxMotifLen,
				nUseCS, nCSLikeNum, vCSLike,
				nIsRecording, nDirec, nMeanMotifLen, 
				nMotifLenUpperBound, nMotifLenLowerBound,
				dDefaultPriorCount);
	}
	else
	{
		/* return */
		return PROC_SUCCESS;
	}

	/* return */
	return PROC_SUCCESS;
}

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
				double dDefaultPriorCount)
{
	/* define */
	/* for data access */
	unsigned char *pBase,*pAi,*pSi,*pCS;
	double *pBL;
	struct FLEXMOTIFSITE **vSite;
	struct FLEXMOTIFSITE *pSite;
	int ni,nj,nPos,nNewPos,nBaseId;
	
	/* new context */
	struct INTMATRIX *pTestPos;
	struct INTMATRIX *pTestMotif;
	struct INTMATRIX *pTestStrand;
	struct BYTEMATRIX *pTestA;
	double dScore, dRatio, dPAi;
	
	/* motif parameters */
	int nMId,nLen;
	int nOldMotifLen,nNewMotifLen;
	int nOldCSShift,nNewCSShift;
	
	/* for local move */
	double dRand;
	
	/* for likelihood calculation */
	double dLike[2];
	struct DOUBLEMATRIX *pMCount;
	int nMask;
	int nCSId;
	int nBaseTypeNum;

	/* do initialization first */
	nBaseTypeNum = pLogBG->nWidth;
	nMId = nMotifId-1;
	nOldMotifLen = pMotifLen->pMatElement[nMId];
	nOldCSShift = (int)((nOldMotifLen-1)/2);
	nNewMotifLen = nOldMotifLen+1;
	nNewCSShift = (int)((nNewMotifLen-1)/2);

	if(nNewMotifLen > nMotifLenUpperBound)
	{
		return PROC_SUCCESS;
	}
	if(nNewMotifLen < nMotifLenLowerBound)
	{
		printf("Error: FlexModule_IncreaseW, motif length out of range!\n");
		exit(EXIT_FAILURE);
	}

	/* vectors below are arranged in the order of 0, 1/-1 */
	for(ni=0; ni<2; ni++)
	{
		dLike[ni] = 0.0;
	}

	pMCount = NULL;
	pMCount = CreateDoubleMatrix(1, nBaseTypeNum);
	if(pMCount == NULL)
	{
		printf("Error: FlexModule_IncreaseW, cannot prepare space for likelihood computation!\n");
		exit(EXIT_FAILURE);
	}


	/* scan sequence by sequence */
	nMask = 0;
	for(ni=0; ni<nSeqCount; ni++)
	{
		/* if there are no sites, jump to the next sequence */
		if(vSeqMtf[ni]->nSiteNum == 0)
		{
			if(vSeqMtf[ni]->pMotifList != NULL)
			{
				printf("Error: FlexModule_IncreaseW, motif site counting error!\n");
				exit(EXIT_FAILURE);
			}
			continue;
		}

		/* if already masked, break */
		if(nMask == 1)
		{
			break;
		}

		/* sequence wise init */
		nLen = vSeqMtf[ni]->vSeq[nIndex]->nWidth;
		pBase = vSeqMtf[ni]->vSeq[nIndex]->pMatElement;
		if(nUseCS == 1)
		{
			pCS = vSeqMtf[ni]->vScore[0]->pMatElement;
		}
		pAi = vSeqMtf[ni]->vStatus[0]->pMatElement;
		pSi = vSeqMtf[ni]->vStatus[1]->pMatElement;
		vSite = vSeqMtf[ni]->vMotif;
		pBL = vSeqMtf[ni]->vMonitor[0]->pMatElement;

		/* init whole sequence context */		
		pTestPos = NULL;
		pTestPos = CreateIntMatrix(1, vSeqMtf[ni]->nSiteNum);
		pTestMotif = NULL;
		pTestMotif = CreateIntMatrix(1, vSeqMtf[ni]->nSiteNum);
		pTestStrand = NULL;
		pTestStrand = CreateIntMatrix(1, vSeqMtf[ni]->nSiteNum);
		pTestA = NULL;
		pTestA = CreateByteMatrix(1, vSeqMtf[ni]->nSiteNum);

		/* scan candidate site one by one */
		pSite = vSeqMtf[ni]->pMotifList;
		nj = 0;
		while(pSite != NULL)
		{
			/* store sites */
			nPos = pSite->nStartPos;
			pTestPos->pMatElement[nj] = nPos;
			pTestMotif->pMatElement[nj] = pSite->nMotifType;
			pTestStrand->pMatElement[nj] = pSite->nStrand;
			pTestA->pMatElement[nj] = pAi[nPos];

			/* get context score */
			dScore = FlexModule_ComputeModuleScore(pSite->pContextPos->nWidth,
				pSite->pContextPos->pMatElement, pSite->pContextMotif->pMatElement,
				pSite->pContextStrand->pMatElement, 0, 0, 0, 0, nMotifNum);
			dPAi = FlexModule_ComputeTrueProb(dScore, dModuleD, dModuleL, dModuleA, dT);
			if(pAi[nPos] == 1)
			{
				dLike[0] += log(dPAi);
			}
			else
			{
				dLike[0] += log(1.0-dPAi);
			}

			/* if the motif site in question */
			if( pSite->nMotifType == nMotifId ) 
			{
				/* count bases */
				/* '-' strand */
				if(pSite->nStrand == 1)
				{
					/* +1 */
					if(nDirec == 1)
					{
						/* move +nk */
						pTestPos->pMatElement[nj] = nPos-1;
						nNewPos = nPos-1;
						if( nNewPos < 0 )
						{
							nMask = 1;
							break;
						}

						if(pBase[nNewPos] >= nBaseTypeNum)
						{
							nMask = 1;
							break;
						}

						if(pAi[nPos] == 1)
						{
							dLike[1] -= pBL[nNewPos];
							
							if(nUseCS == 1)
							{
								nCSId = pCS[nNewPos+nNewCSShift];
								dLike[1] += (vCSLike[1]->pMatElement[nCSId]-vCSLike[0]->pMatElement[nCSId]);
							}

							/* count bases */
							nBaseId = nBaseTypeNum-1-(int)pBase[nNewPos];
							pMCount->pMatElement[nBaseId] += 1.0;
						}
					}
					/* -1 */
					else if(nDirec == -1)
					{
						/* move -nk */
						pTestPos->pMatElement[nj] = nPos;
						nNewPos = nPos+nNewMotifLen-1;
						if( nNewPos >= nLen )
						{
							nMask = 1;
							break;
						}

						if(pBase[nNewPos] >= nBaseTypeNum)
						{
							nMask = 1;
							break;
						}
						
						if(pAi[nPos] == 1)
						{
							dLike[1] -= pBL[nNewPos];
							
							if(nUseCS == 1)
							{
								nCSId = pCS[nPos+nNewCSShift];
								dLike[1] += (vCSLike[1]->pMatElement[nCSId]-vCSLike[0]->pMatElement[nCSId]);
							}

							/* count bases */
							nBaseId = nBaseTypeNum-1-(int)pBase[nNewPos];
							pMCount->pMatElement[nBaseId] += 1.0;
						}
					}
					/* otherwise */
					else
					{
						printf("Error: FlexModule_IncreaseW, sample strand wrong!\n");
						exit(EXIT_FAILURE);
					}
				}
				
				/* '+' strand */
				else
				{
					/* +1 */
					if(nDirec == 1)
					{
						/* move +nk */
						pTestPos->pMatElement[nj] = nPos;
						nNewPos = nPos+nNewMotifLen-1;
						if( nNewPos >= nLen )
						{
							nMask = 1;
							break;
						}

						if(pBase[nNewPos] >= nBaseTypeNum)
						{
							nMask = 1;
							break;
						}
						
						if(pAi[nPos] == 1)
						{
							dLike[1] -= pBL[nNewPos];
							
							if(nUseCS == 1)
							{
								nCSId = pCS[nPos+nNewCSShift];
								dLike[1] += (vCSLike[1]->pMatElement[nCSId]-vCSLike[0]->pMatElement[nCSId]);
							}

							/* count bases */
							nBaseId = (int)pBase[nNewPos];
							pMCount->pMatElement[nBaseId] += 1.0;
						}
					}
					/* -1 */
					else if(nDirec == -1)
					{
						/* move -nk */
						pTestPos->pMatElement[nj] = nPos-1;
						nNewPos = nPos-1;
						if( nNewPos < 0 )
						{
							nMask = 1;
							break;
						}

						if(pBase[nNewPos] >= nBaseTypeNum)
						{
							nMask = 1;
							break;
						}

						if(pAi[nPos] == 1)
						{
							dLike[1] -= pBL[nNewPos];
							
							if(nUseCS == 1)
							{
								nCSId = pCS[nNewPos+nNewCSShift];
								dLike[1] += (vCSLike[1]->pMatElement[nCSId]-vCSLike[0]->pMatElement[nCSId]);
							}

							/* count bases */
							nBaseId = (int)pBase[nNewPos];
							pMCount->pMatElement[nBaseId] += 1.0;
						}
					}
					/* otherwise */
					else
					{
						printf("Error: FlexModule_IncreaseW, sample strand wrong!\n");
						exit(EXIT_FAILURE);
					}
				}

				/* change dLike0 */
				if( (nUseCS == 1) && (pAi[nPos] == 1) )
				{
					nCSId = pCS[nPos+nOldCSShift];
					dLike[0] += (vCSLike[1]->pMatElement[nCSId]-vCSLike[0]->pMatElement[nCSId]);
				}
			}

			/* get next */
			nj++;
			pSite = pSite->pNext;
		}

		/* if already masked, stop scan */
		if(nMask == 1)
		{
			/* destroy whole sequence context */
			DestroyIntMatrix(pTestPos);
			DestroyIntMatrix(pTestMotif);
			DestroyIntMatrix(pTestStrand);
			DestroyByteMatrix(pTestA);
			break;
		}

		if(nj != vSeqMtf[ni]->nSiteNum)
		{
			printf("Error: FlexModule_IncreaseW, motif site counting error!\n");
			exit(EXIT_FAILURE);
		}


		/* get new context likelihood */
		pMotifLen->pMatElement[nMId] = nNewMotifLen;
		dScore = FlexModule_ShiftContextLogLike(pTestPos, pTestMotif, pTestStrand, pTestA,
			nLen, pMotifLen, nModuleLen, dModuleD, dModuleL, dModuleA, dT, &nMask);
		pMotifLen->pMatElement[nMId] = nOldMotifLen;

		/* if already masked, stop scan */
		if(nMask == 1)
		{
			/* destroy whole sequence context */
			DestroyIntMatrix(pTestPos);
			DestroyIntMatrix(pTestMotif);
			DestroyIntMatrix(pTestStrand);
			DestroyByteMatrix(pTestA);
			break;
		}

		dLike[1] += dScore;
			
		/* destroy whole sequence context */
		DestroyIntMatrix(pTestPos);
		DestroyIntMatrix(pTestMotif);
		DestroyIntMatrix(pTestStrand);
		DestroyByteMatrix(pTestA);
	}

	/* if masked, return */
	if(nMask == 1)
	{
		DestroyDoubleMatrix(pMCount);
		return PROC_SUCCESS;
	}

	/* update the whole count likelihood */
	dScore = FlexModule_IncreaseWMotifLogLike(vMotif[nMId], pMCount, 
		nDirec, dDefaultPriorCount);
	dLike[1] += dScore;
	dLike[1] += log((double)nMeanMotifLen)-log((double)nNewMotifLen);
	
	/* sample shift */
	dRatio = exp(dLike[1]-dLike[0]);
	dRand = rand_u();

	/* update status */
	if(dRand <= dRatio)
	{
		/* update sequence by sequence */
		for(ni=0; ni<nSeqCount; ni++)
		{
			FlexModule_IncreaseWUpdate(vSeqMtf[ni], nIndex, nMotifId, nDirec,
				nModuleLen, pMotifLen, pMaxMotifLen);
		}

		/* update motif length */
		pMotifLen->pMatElement[nMId] += 1;
		if(pMotifLen->pMatElement[nMId] > (*pMaxMotifLen))
			*pMaxMotifLen = pMotifLen->pMatElement[nMId];

		/* update PWM */
		FlexModule_IncreaseWUpdateMotif(vMotif[nMId], pMCount, 
			nDirec, dDefaultPriorCount);
	}

	/* release memory */
	DestroyDoubleMatrix(pMCount);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_IncreaseWMotifLogLike: compute length-increased motif's     */
/*  loglikelihood.                                                         */
/* ----------------------------------------------------------------------- */ 
double FlexModule_IncreaseWMotifLogLike(struct FLEXMOTIFMATRIX *pMotif, 
							 struct DOUBLEMATRIX *pMCount, int nDirec,
							 double dDefaultPriorCount)
{
	/* define */
	double dL = 0.0;
	double dTotal,dTemp,dTotal2;
	double *pEle1,*pEle2;
	int nj,nBaseTypeNum;

	/* get loglikelihood */
	nBaseTypeNum = pMotif->pPWM->nWidth;
	/* if extending the tail */
	if(nDirec == 1)
	{
		if(pMotif->nTailExtraNum > 0)
		{
			if(pMotif->nTailExtraNum != pMotif->pPriorTailExtra->nHeight)
			{
				printf("Error: FlexModule_IncreaseWMotifLogLike, extra prior vector not match!\n");
				exit(EXIT_FAILURE);
			}

			pEle2 = pMotif->pPriorTailExtra->pMatElement+nBaseTypeNum*(pMotif->nTailExtraNum-1);

			pEle1 = pMCount->pMatElement;
			dTotal = 0.0;
			dTotal2 = 0.0;
			for(nj=0; nj<nBaseTypeNum; nj++)
			{
				dTemp = (*pEle1)+(*pEle2);
				dTotal += dTemp;
				dTotal2 += (*pEle2);
				dL += gammaln(dTemp);
				dL -= gammaln(*pEle2);
				pEle1++;
				pEle2++;
			}
			dL -= gammaln(dTotal);
			dL += gammaln(dTotal2);
		}
		else
		{
			pEle1 = pMCount->pMatElement;
			dTotal = 0.0;
			dTotal2 = 0.0;
			for(nj=0; nj<nBaseTypeNum; nj++)
			{
				dTemp = (*pEle1)+dDefaultPriorCount;
				dTotal += dTemp;
				dTotal2 += dDefaultPriorCount;
				dL += gammaln(dTemp);
				dL -= gammaln(dDefaultPriorCount);
				pEle1++;
			}
			dL -= gammaln(dTotal);
			dL += gammaln(dTotal2);
		}
	}
	/* if extending the head */
	else if(nDirec == -1)
	{
		if(pMotif->nHeadExtraNum > 0)
		{
			if(pMotif->nHeadExtraNum != pMotif->pPriorHeadExtra->nHeight)
			{
				printf("Error: FlexModule_IncreaseWMotifLogLike, extra prior vector not match!\n");
				exit(EXIT_FAILURE);
			}

			pEle2 = pMotif->pPriorHeadExtra->pMatElement+nBaseTypeNum*(pMotif->nHeadExtraNum-1);

			pEle1 = pMCount->pMatElement;
			dTotal = 0.0;
			dTotal2 = 0.0;
			for(nj=0; nj<nBaseTypeNum; nj++)
			{
				dTemp = (*pEle1)+(*pEle2);
				dTotal += dTemp;
				dTotal2 += (*pEle2);
				dL += gammaln(dTemp);
				dL -= gammaln(*pEle2);
				pEle1++;
				pEle2++;
			}
			dL -= gammaln(dTotal);
			dL += gammaln(dTotal2);
		}
		else
		{
			pEle1 = pMCount->pMatElement;
			dTotal = 0.0;
			dTotal2 = 0.0;
			for(nj=0; nj<nBaseTypeNum; nj++)
			{
				dTemp = (*pEle1)+dDefaultPriorCount;
				dTotal += dTemp;
				dTotal2 += dDefaultPriorCount;
				dL += gammaln(dTemp);
				dL -= gammaln(dDefaultPriorCount);
				pEle1++;
			}
			dL -= gammaln(dTotal);
			dL += gammaln(dTotal2);
		}
	}
	/* otherwise error */
	else
	{
		printf("Error: FlexModule_IncreaseWMotifLogLike, the direction of extending motif length is not clear!\n");
		exit(EXIT_FAILURE);
	}

	/* return */
	return dL;
}


/* ----------------------------------------------------------------------- */ 
/*  FlexModule_IncreaseWUpdate: update all indicators and context after    */
/*  increasing motif length by 1.                                          */
/* ----------------------------------------------------------------------- */ 
int FlexModule_IncreaseWUpdate(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, 
				int nMotifId, int nDirec,
				int nModuleLen, struct INTMATRIX *pMotifLen, int *pMaxMotifLen)
{
	/* define */
	unsigned char *pBase,*pAi,*pSi;
	struct FLEXMOTIFSITE **vSite;
	struct FLEXMOTIFSITE *pSite;
	int nLen,nMId,nMotifLen;
	int ni,nSiteNum,nOldPos,nNewPos,nContextNum;
	int nHalfWin;
	struct FLEXMOTIFSITE *pHead,*pTail,*pTest;
	

	/* init */
	nHalfWin = (int)(nModuleLen/2);
	nMId = nMotifId-1;
	nMotifLen = pMotifLen->pMatElement[nMId];
	nLen = pSeqMtf->vSeq[nIndex]->nWidth;
	pBase = pSeqMtf->vSeq[nIndex]->pMatElement;
	pAi = pSeqMtf->vStatus[0]->pMatElement;
	pSi = pSeqMtf->vStatus[1]->pMatElement;
	vSite = pSeqMtf->vMotif;

	/* scan all sites */
	pSite = pSeqMtf->pMotifList;
	nSiteNum = 0;
	while(pSite != NULL)
	{
		/* if motif in question, update */
		if(pSite->nMotifType == nMotifId)
		{
			/* if - strand */
			if(pSite->nStrand == 1)
			{
				nOldPos = pSite->nStartPos;

				if(nDirec == 1)
				{
					nNewPos = nOldPos-1;
				
					/* update Ai */
					pAi[nNewPos] = pAi[nOldPos];
					pAi[nOldPos] = 0;

					/* update Si */
					pSi[nNewPos] = nMotifId;
					
					/* update vSite */
					vSite[nNewPos] = vSite[nOldPos];
					vSite[nOldPos] = NULL;

					/* update site itself */
					pSite->nStartPos = nNewPos;
				}
				else if(nDirec == -1)
				{
					/* update Ai */
					
					/* update Si */
					pSi[nOldPos+nMotifLen] = nMotifId;
					
					/* update vSite */
					
					/* update site itself */
				}
				else
				{
					printf("Error: FlexModule_IncreaseWUpdate, direction of extending motif length is not clear!\n");
					exit(EXIT_FAILURE);
				}
			}
			/* if + strand */
			else
			{
				nOldPos = pSite->nStartPos;

				if(nDirec == 1)
				{
					/* update Ai */
					
					/* update Si */
					pSi[nOldPos+nMotifLen] = nMotifId;
					
					/* update vSite */
					
					/* update site itself */
				}
				else if(nDirec == -1)
				{
					nNewPos = nOldPos-1;
				
					/* update Ai */
					pAi[nNewPos] = pAi[nOldPos];
					pAi[nOldPos] = 0;

					/* update Si */
					pSi[nNewPos] = nMotifId;
					
					/* update vSite */
					vSite[nNewPos] = vSite[nOldPos];
					vSite[nOldPos] = NULL;

					/* update site itself */
					pSite->nStartPos = nNewPos;
				}
				else
				{
					printf("Error: FlexModule_IncreaseWUpdate, direction of extending motif length is not clear!\n");
					exit(EXIT_FAILURE);
				}
			}
		}

		/* get next */
		nSiteNum++;
		pSite = pSite->pNext;
	}

	if(nSiteNum != pSeqMtf->nSiteNum)
	{
		printf("Error: FlexModule_IncreaseWUpdate, site number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* update context */
	pSite = pSeqMtf->pMotifList;
	while(pSite != NULL)
	{
		/* establish context limit */
		nContextNum = 1;
		pHead = pSite;
		pTail = pSite;

		pTest = pHead->pPrev;
		while(pTest != NULL)
		{
			if( (pSite->nStartPos-pTest->nStartPos) > nHalfWin )
				break;

			pHead = pTest;
			pTest = pHead->pPrev;
			nContextNum++;
		}

		pTest = pTail->pNext;
		while(pTest != NULL)
		{
			if( (pTest->nStartPos-pSite->nStartPos) > nHalfWin )
				break;

			pTail = pTest;
			pTest = pTail->pNext;
			nContextNum++;
		}

		/* create context matrix */
		DestroyIntMatrix(pSite->pContextPos);
		pSite->pContextPos = NULL;
		pSite->pContextPos = CreateIntMatrix(1, nContextNum);
		if(pSite->pContextPos == NULL)
		{
			printf("Error: FlexModule_IncreaseWUpdate, cannot create context matrix correctly!\n");
			exit(EXIT_FAILURE);
		}
		DestroyIntMatrix(pSite->pContextMotif);
		pSite->pContextMotif = NULL;
		pSite->pContextMotif = CreateIntMatrix(1, nContextNum);
		if(pSite->pContextMotif == NULL)
		{
			printf("Error: FlexModule_IncreaseWUpdate, cannot create context matrix correctly!\n");
			exit(EXIT_FAILURE);
		}
		DestroyIntMatrix(pSite->pContextStrand);
		pSite->pContextStrand = NULL;
		pSite->pContextStrand = CreateIntMatrix(1, nContextNum);
		if(pSite->pContextStrand == NULL)
		{
			printf("Error: FlexModule_IncreaseWUpdate, cannot create context matrix correctly!\n");
			exit(EXIT_FAILURE);
		}

		/* fill out the context */
		pTest = pHead;
		ni = 0;
		while(pTest != NULL)
		{
			if(pTest->nStartPos - pTail->nStartPos > 0)
				break;

			pSite->pContextPos->pMatElement[ni] = pTest->nStartPos-pSite->nStartPos;
			pSite->pContextMotif->pMatElement[ni] = pTest->nMotifType;
			pSite->pContextStrand->pMatElement[ni] = pTest->nStrand;
			ni++;
			pTest = pTest->pNext;
		}

		if(ni != nContextNum)
		{
			printf("Error: FlexModule_IncreaseWUpdate, number of motif sites in the context do not match!\n");
			exit(EXIT_FAILURE);
		}

		/* get next */
		pSite = pSite->pNext;
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_IncreaseWUpdateMotif: update motif after increasing its     */
/*  length by one.                                                         */
/* ----------------------------------------------------------------------- */ 
int FlexModule_IncreaseWUpdateMotif(struct FLEXMOTIFMATRIX *pMotif, 
							 struct DOUBLEMATRIX *pMCount, int nDirec,
							 double dDefaultPriorCount)
{
	/* define */
	double dL = 0.0;
	double *pEle1,*pEle2;
	int nj,nMotifLen,nBaseTypeNum;
	struct DOUBLEMATRIX *pNewPriorCount;
	struct DOUBLEMATRIX *pNewSampleCount;
	struct DOUBLEMATRIX *pNewExtraPrior;

	/* init */
	nMotifLen = pMotif->pPWM->nHeight;
	nBaseTypeNum = pMotif->pPWM->nWidth;

	/* get loglikelihood */
	/* if extending the tail */
	if(nDirec == 1)
	{
		/* copy the old counts */
		pNewPriorCount = NULL;
		pNewPriorCount = CreateDoubleMatrix((nMotifLen+1), nBaseTypeNum);
		if(pNewPriorCount == NULL)
		{
			printf("Error: FlexModule_IncreaseWUpdateMotif, cannot create new prior count!\n");
			exit(EXIT_FAILURE);
		}
		pNewSampleCount = NULL;
		pNewSampleCount = CreateDoubleMatrix((nMotifLen+1), nBaseTypeNum);
		if(pNewSampleCount == NULL)
		{
			printf("Error: FlexModule_IncreaseWUpdateMotif, cannot create new sample count!\n");
			exit(EXIT_FAILURE);
		}

		pEle1 = pNewPriorCount->pMatElement;
		pEle2 = pMotif->pPriorCount->pMatElement;
		memcpy(pEle1, pEle2, nMotifLen*nBaseTypeNum*sizeof(double));

		pEle1 = pNewSampleCount->pMatElement;
		pEle2 = pMotif->pSampleCount->pMatElement;
		memcpy(pEle1, pEle2, nMotifLen*nBaseTypeNum*sizeof(double));

		pEle1 = pNewSampleCount->pMatElement+nMotifLen*nBaseTypeNum;
		pEle2 = pMCount->pMatElement;
		memcpy(pEle1, pEle2, nBaseTypeNum*sizeof(double));

		/* add new counts */
		if(pMotif->nTailExtraNum > 0)
		{
			if(pMotif->nTailExtraNum != pMotif->pPriorTailExtra->nHeight)
			{
				printf("Error: FlexModule_IncreaseWUpdateMotif, extra prior vector not match!\n");
				exit(EXIT_FAILURE);
			}

			pEle1 = pNewPriorCount->pMatElement+nMotifLen*nBaseTypeNum;
			pEle2 = pMotif->pPriorTailExtra->pMatElement+nBaseTypeNum*(pMotif->nTailExtraNum-1);
			memcpy(pEle1, pEle2, nBaseTypeNum*sizeof(double));
		}
		else
		{
			pEle1 = pNewPriorCount->pMatElement+nMotifLen*nBaseTypeNum;
			for(nj=0; nj<nBaseTypeNum; nj++)
			{
				*pEle1 = dDefaultPriorCount;
				pEle1++;
			}
		}

		/* update PWM */
		DestroyDoubleMatrix(pMotif->pPriorCount);
		pMotif->pPriorCount = pNewPriorCount;
		DestroyDoubleMatrix(pMotif->pSampleCount);
		pMotif->pSampleCount = pNewSampleCount;
		DestroyDoubleMatrix(pMotif->pPWM);
		pMotif->pPWM = NULL;
		pMotif->pPWM = CreateDoubleMatrix((nMotifLen+1), nBaseTypeNum);
		if(pMotif->pPWM == NULL)
		{
			printf("Error: FlexModule_IncreaseWUpdateMotif, cannot create new PWM count!\n");
			exit(EXIT_FAILURE);
		}
		FLEXMTFMREFRESH(pMotif);

		/* update extra prior pool */
		pMotif->nTailExtraNum -= 1;
		if(pMotif->nTailExtraNum > 0)
		{
			pNewExtraPrior = NULL;
			pNewExtraPrior = CreateDoubleMatrix(pMotif->nTailExtraNum, nBaseTypeNum);
			if(pNewExtraPrior == NULL)
			{
				printf("Error: FlexModule_IncreaseWUpdateMotif, cannot create new extra prior count!\n");
				exit(EXIT_FAILURE);
			}

			pEle1 = pNewExtraPrior->pMatElement;
			pEle2 = pMotif->pPriorTailExtra->pMatElement;
			memcpy(pEle1, pEle2, (pMotif->nTailExtraNum)*nBaseTypeNum*sizeof(double));
			DestroyDoubleMatrix(pMotif->pPriorTailExtra);
			pMotif->pPriorTailExtra = pNewExtraPrior;

		}
		else if(pMotif->nTailExtraNum == 0)
		{
			DestroyDoubleMatrix(pMotif->pPriorTailExtra);
			pMotif->pPriorTailExtra = NULL;
		}
		else
		{
		}
	}
	/* if extending the head */
	else if(nDirec == -1)
	{
		/* copy the old counts */
		pNewPriorCount = NULL;
		pNewPriorCount = CreateDoubleMatrix((nMotifLen+1), nBaseTypeNum);
		if(pNewPriorCount == NULL)
		{
			printf("Error: FlexModule_IncreaseWUpdateMotif, cannot create new prior count!\n");
			exit(EXIT_FAILURE);
		}
		pNewSampleCount = NULL;
		pNewSampleCount = CreateDoubleMatrix((nMotifLen+1), nBaseTypeNum);
		if(pNewSampleCount == NULL)
		{
			printf("Error: FlexModule_IncreaseWUpdateMotif, cannot create new sample count!\n");
			exit(EXIT_FAILURE);
		}

		pEle1 = pNewPriorCount->pMatElement+nBaseTypeNum;
		pEle2 = pMotif->pPriorCount->pMatElement;
		memcpy(pEle1, pEle2, nMotifLen*nBaseTypeNum*sizeof(double));

		pEle1 = pNewSampleCount->pMatElement;
		pEle2 = pMCount->pMatElement;
		memcpy(pEle1, pEle2, nBaseTypeNum*sizeof(double));

		pEle1 = pNewSampleCount->pMatElement+nBaseTypeNum;
		pEle2 = pMotif->pSampleCount->pMatElement;
		memcpy(pEle1, pEle2, nMotifLen*nBaseTypeNum*sizeof(double));

		/* add new counts */
		if(pMotif->nHeadExtraNum > 0)
		{
			if(pMotif->nHeadExtraNum != pMotif->pPriorHeadExtra->nHeight)
			{
				printf("Error: FlexModule_IncreaseWUpdateMotif, extra prior vector not match!\n");
				exit(EXIT_FAILURE);
			}

			pEle1 = pNewPriorCount->pMatElement;
			pEle2 = pMotif->pPriorHeadExtra->pMatElement+nBaseTypeNum*(pMotif->nHeadExtraNum-1);
			memcpy(pEle1, pEle2, nBaseTypeNum*sizeof(double));
		}
		else
		{
			pEle1 = pNewPriorCount->pMatElement;
			for(nj=0; nj<nBaseTypeNum; nj++)
			{
				*pEle1 = dDefaultPriorCount;
				pEle1++;
			}
		}

		/* update PWM */
		DestroyDoubleMatrix(pMotif->pPriorCount);
		pMotif->pPriorCount = pNewPriorCount;
		DestroyDoubleMatrix(pMotif->pSampleCount);
		pMotif->pSampleCount = pNewSampleCount;
		DestroyDoubleMatrix(pMotif->pPWM);
		pMotif->pPWM = NULL;
		pMotif->pPWM = CreateDoubleMatrix((nMotifLen+1), nBaseTypeNum);
		if(pMotif->pPWM == NULL)
		{
			printf("Error: FlexModule_IncreaseWUpdateMotif, cannot create new PWM count!\n");
			exit(EXIT_FAILURE);
		}
		FLEXMTFMREFRESH(pMotif);

		/* update extra prior pool */
		pMotif->nHeadExtraNum -= 1;
		if(pMotif->nHeadExtraNum > 0)
		{
			pNewExtraPrior = NULL;
			pNewExtraPrior = CreateDoubleMatrix(pMotif->nHeadExtraNum, nBaseTypeNum);
			if(pNewExtraPrior == NULL)
			{
				printf("Error: FlexModule_IncreaseWUpdateMotif, cannot create new extra prior count!\n");
				exit(EXIT_FAILURE);
			}

			pEle1 = pNewExtraPrior->pMatElement;
			pEle2 = pMotif->pPriorHeadExtra->pMatElement;
			memcpy(pEle1, pEle2, (pMotif->nHeadExtraNum)*nBaseTypeNum*sizeof(double));
			DestroyDoubleMatrix(pMotif->pPriorHeadExtra);
			pMotif->pPriorHeadExtra = pNewExtraPrior;
		}
		else if(pMotif->nHeadExtraNum == 0)
		{
			DestroyDoubleMatrix(pMotif->pPriorHeadExtra);
			pMotif->pPriorHeadExtra = NULL;
		}
		else
		{
		}
	}
	/* otherwise error */
	else
	{
		printf("Error: FlexModule_IncreaseWUpdateMotif, the direction of extending motif length is not clear!\n");
		exit(EXIT_FAILURE);
	}

	/* return */
	return PROC_SUCCESS;
}


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
				double dDefaultPriorCount)
{
	/* for data access */
	unsigned char *pBase,*pAi,*pSi,*pCS;
	double *pBL;
	struct FLEXMOTIFSITE **vSite;
	struct FLEXMOTIFSITE *pSite;
	int ni,nj,nPos,nNewPos,nBaseId;
	
	/* new context */
	struct INTMATRIX *pTestPos;
	struct INTMATRIX *pTestMotif;
	struct INTMATRIX *pTestStrand;
	struct BYTEMATRIX *pTestA;
	double dScore, dRatio, dPAi;
	
	/* motif parameters */
	int nMId,nLen;
	int nOldMotifLen,nNewMotifLen;
	int nOldCSShift,nNewCSShift;
	
	/* for local move */
	double dRand;
	
	/* for likelihood calculation */
	double dLike[2];
	struct DOUBLEMATRIX *pMCount;
	int nMask;
	int nCSId;
	int nBaseTypeNum;

	/* do initialization first */
	nBaseTypeNum = pLogBG->nWidth;
	nMId = nMotifId-1;
	nOldMotifLen = pMotifLen->pMatElement[nMId];
	nOldCSShift = (int)((nOldMotifLen-1)/2);
	nNewMotifLen = nOldMotifLen-1;
	nNewCSShift = (int)((nNewMotifLen-1)/2);

	if(nNewMotifLen < nMotifLenLowerBound)
	{
		return PROC_SUCCESS;
	}
	if(nNewMotifLen > nMotifLenUpperBound)
	{
		printf("Error: FlexModule_DecreaseW, motif length out of range!\n");
		exit(EXIT_FAILURE);
	}

	/* vectors below are arranged in the order of 0, 1/-1 */
	for(ni=0; ni<2; ni++)
	{
		dLike[ni] = 0.0;
	}

	pMCount = NULL;
	pMCount = CreateDoubleMatrix(1, nBaseTypeNum);
	if(pMCount == NULL)
	{
		printf("Error: FlexModule_DecreaseW, cannot prepare space for likelihood computation!\n");
		exit(EXIT_FAILURE);
	}


	/* scan sequence by sequence */
	nMask = 0;
	for(ni=0; ni<nSeqCount; ni++)
	{
		/* if there are no sites, jump to the next sequence */
		if(vSeqMtf[ni]->nSiteNum == 0)
		{
			if(vSeqMtf[ni]->pMotifList != NULL)
			{
				printf("Error: FlexModule_DecreaseW, motif site counting error!\n");
				exit(EXIT_FAILURE);
			}
			continue;
		}

		/* if already masked, break */
		if(nMask == 1)
		{
			break;
		}

		/* sequence wise init */
		nLen = vSeqMtf[ni]->vSeq[nIndex]->nWidth;
		pBase = vSeqMtf[ni]->vSeq[nIndex]->pMatElement;
		if(nUseCS == 1)
		{
			pCS = vSeqMtf[ni]->vScore[0]->pMatElement;
		}
		pAi = vSeqMtf[ni]->vStatus[0]->pMatElement;
		pSi = vSeqMtf[ni]->vStatus[1]->pMatElement;
		vSite = vSeqMtf[ni]->vMotif;
		pBL = vSeqMtf[ni]->vMonitor[0]->pMatElement;

		/* init whole sequence context */		
		pTestPos = NULL;
		pTestPos = CreateIntMatrix(1, vSeqMtf[ni]->nSiteNum);
		pTestMotif = NULL;
		pTestMotif = CreateIntMatrix(1, vSeqMtf[ni]->nSiteNum);
		pTestStrand = NULL;
		pTestStrand = CreateIntMatrix(1, vSeqMtf[ni]->nSiteNum);
		pTestA = NULL;
		pTestA = CreateByteMatrix(1, vSeqMtf[ni]->nSiteNum);

		/* scan candidate site one by one */
		pSite = vSeqMtf[ni]->pMotifList;
		nj = 0;
		while(pSite != NULL)
		{
			/* store sites */
			nPos = pSite->nStartPos;
			pTestPos->pMatElement[nj] = nPos;
			pTestMotif->pMatElement[nj] = pSite->nMotifType;
			pTestStrand->pMatElement[nj] = pSite->nStrand;
			pTestA->pMatElement[nj] = pAi[nPos];

			/* get context score */
			dScore = FlexModule_ComputeModuleScore(pSite->pContextPos->nWidth,
				pSite->pContextPos->pMatElement, pSite->pContextMotif->pMatElement,
				pSite->pContextStrand->pMatElement, 0, 0, 0, 0, nMotifNum);
			dPAi = FlexModule_ComputeTrueProb(dScore, dModuleD, dModuleL, dModuleA, dT);
			if(pAi[nPos] == 1)
			{
				dLike[0] += log(dPAi);
			}
			else
			{
				dLike[0] += log(1.0-dPAi);
			}

			/* if the motif site in question */
			if( pSite->nMotifType == nMotifId ) 
			{
				/* count bases */
				/* '-' strand */
				if(pSite->nStrand == 1)
				{
					/* +1 */
					if(nDirec == 1)
					{
						/* move +nk */
						pTestPos->pMatElement[nj] = nPos+1;
						nNewPos = nPos+1;
						
						if( nNewPos >= nLen )
						{
							nMask = 1;
							break;
						}

						if(pBase[nNewPos] >= nBaseTypeNum)
						{
							nMask = 1;
							break;
						}

						if(pAi[nPos] == 1)
						{
							dLike[1] += pBL[nPos];
							
							if(nUseCS == 1)
							{
								nCSId = pCS[nNewPos+nNewCSShift];
								dLike[1] += (vCSLike[1]->pMatElement[nCSId]-vCSLike[0]->pMatElement[nCSId]);
							}

							/* count bases */
							nBaseId = nBaseTypeNum-1-(int)pBase[nPos];
							pMCount->pMatElement[nBaseId] += 1.0;
						}
					}
					/* -1 */
					else if(nDirec == -1)
					{
						/* move -nk */
						pTestPos->pMatElement[nj] = nPos;
						nNewPos = nPos+nOldMotifLen-1;
						if( nNewPos >= nLen )
						{
							nMask = 1;
							break;
						}

						if(pBase[nNewPos] >= nBaseTypeNum)
						{
							nMask = 1;
							break;
						}
						
						if(pAi[nPos] == 1)
						{
							dLike[1] += pBL[nNewPos];
							
							if(nUseCS == 1)
							{
								nCSId = pCS[nPos+nNewCSShift];
								dLike[1] += (vCSLike[1]->pMatElement[nCSId]-vCSLike[0]->pMatElement[nCSId]);
							}

							/* count bases */
							nBaseId = nBaseTypeNum-1-(int)pBase[nNewPos];
							pMCount->pMatElement[nBaseId] += 1.0;
						}
					}
					/* otherwise */
					else
					{
						printf("Error: FlexModule_DecreaseW, sample strand wrong!\n");
						exit(EXIT_FAILURE);
					}
				}
				
				/* '+' strand */
				else
				{
					/* +1 */
					if(nDirec == 1)
					{
						/* move +nk */
						pTestPos->pMatElement[nj] = nPos;
						nNewPos = nPos+nOldMotifLen-1;
						if( nNewPos >= nLen )
						{
							nMask = 1;
							break;
						}

						if(pBase[nNewPos] >= nBaseTypeNum)
						{
							nMask = 1;
							break;
						}
						
						if(pAi[nPos] == 1)
						{
							dLike[1] += pBL[nNewPos];
							
							if(nUseCS == 1)
							{
								nCSId = pCS[nPos+nNewCSShift];
								dLike[1] += (vCSLike[1]->pMatElement[nCSId]-vCSLike[0]->pMatElement[nCSId]);
							}

							/* count bases */
							nBaseId = (int)pBase[nNewPos];
							pMCount->pMatElement[nBaseId] += 1.0;
						}
					}
					/* -1 */
					else if(nDirec == -1)
					{
						/* move -nk */
						pTestPos->pMatElement[nj] = nPos+1;
						nNewPos = nPos+1;
						if( nNewPos >= nLen )
						{
							nMask = 1;
							break;
						}

						if(pBase[nPos] >= nBaseTypeNum)
						{
							nMask = 1;
							break;
						}

						if(pAi[nPos] == 1)
						{
							dLike[1] += pBL[nPos];
							
							if(nUseCS == 1)
							{
								nCSId = pCS[nNewPos+nNewCSShift];
								dLike[1] += (vCSLike[1]->pMatElement[nCSId]-vCSLike[0]->pMatElement[nCSId]);
							}

							/* count bases */
							nBaseId = (int)pBase[nPos];
							pMCount->pMatElement[nBaseId] += 1.0;
						}
					}
					/* otherwise */
					else
					{
						printf("Error: FlexModule_DecreaseW, sample strand wrong!\n");
						exit(EXIT_FAILURE);
					}
				}

				/* change dLike0 */
				if( (nUseCS == 1) && (pAi[nPos] == 1) )
				{
					nCSId = pCS[nPos+nOldCSShift];
					dLike[0] += (vCSLike[1]->pMatElement[nCSId]-vCSLike[0]->pMatElement[nCSId]);
				}
			}

			/* get next */
			nj++;
			pSite = pSite->pNext;
		}

		/* if already masked, stop scan */
		if(nMask == 1)
		{
			/* destroy whole sequence context */
			DestroyIntMatrix(pTestPos);
			DestroyIntMatrix(pTestMotif);
			DestroyIntMatrix(pTestStrand);
			DestroyByteMatrix(pTestA);
			break;
		}

		if(nj != vSeqMtf[ni]->nSiteNum)
		{
			printf("Error: FlexModule_IncreaseW, motif site counting error!\n");
			exit(EXIT_FAILURE);
		}


		/* get new context likelihood */
		pMotifLen->pMatElement[nMId] = nNewMotifLen;
		dScore = FlexModule_ShiftContextLogLike(pTestPos, pTestMotif, pTestStrand, pTestA,
			nLen, pMotifLen, nModuleLen, dModuleD, dModuleL, dModuleA, dT, &nMask);
		pMotifLen->pMatElement[nMId] = nOldMotifLen;

		/* if already masked, stop scan */
		if(nMask == 1)
		{
			/* destroy whole sequence context */
			DestroyIntMatrix(pTestPos);
			DestroyIntMatrix(pTestMotif);
			DestroyIntMatrix(pTestStrand);
			DestroyByteMatrix(pTestA);
			break;
		}

		dLike[1] += dScore;
			
		/* destroy whole sequence context */
		DestroyIntMatrix(pTestPos);
		DestroyIntMatrix(pTestMotif);
		DestroyIntMatrix(pTestStrand);
		DestroyByteMatrix(pTestA);
	}

	/* if masked, return */
	if(nMask == 1)
	{
		DestroyDoubleMatrix(pMCount);
		return PROC_SUCCESS;
	}

	/* update the whole count likelihood */
	dScore = FlexModule_DecreaseWMotifLogLike(vMotif[nMId], pMCount, 
		nDirec, dDefaultPriorCount);
	dLike[1] -= dScore;
	dLike[1] += log((double)nOldMotifLen)-log((double)nMeanMotifLen);
	
	/* sample shift */
	dRatio = exp(dLike[1]-dLike[0]);
	dRand = rand_u();

	/* update status */
	if(dRand <= dRatio)
	{
		/* update sequence by sequence */
		for(ni=0; ni<nSeqCount; ni++)
		{
			FlexModule_DecreaseWUpdate(vSeqMtf[ni], nIndex, nMotifId, nDirec,
				nModuleLen, pMotifLen, pMaxMotifLen);
		}

		/* update motif length */
		pMotifLen->pMatElement[nMId] -= 1;
		*pMaxMotifLen = 0;
		for(ni=0; ni<pMotifLen->nWidth; ni++)
		{
			if(pMotifLen->pMatElement[ni] > (*pMaxMotifLen))
				*pMaxMotifLen = pMotifLen->pMatElement[ni];
		}
	
		/* update PWM */
		FlexModule_DecreaseWUpdateMotif(vMotif[nMId], pMCount, 
			nDirec, dDefaultPriorCount);
	}

	/* release memory */
	DestroyDoubleMatrix(pMCount);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_DecreaseWMotifLogLike: compute length-decreased motif's     */
/*  loglikelihood.                                                         */
/* ----------------------------------------------------------------------- */ 
double FlexModule_DecreaseWMotifLogLike(struct FLEXMOTIFMATRIX *pMotif, 
							 struct DOUBLEMATRIX *pMCount, int nDirec,
							 double dDefaultPriorCount)
{
	/* define */
	double dL = 0.0;
	double dTotal,dTemp,dTotal2;
	double *pEle1,*pEle2,*pEle3;
	int nj,nBaseTypeNum;
	int nOffset;

	/* get loglikelihood */
	nBaseTypeNum = pMotif->pPWM->nWidth;
	/* if shrinking the tail */
	if(nDirec == 1)
	{
		nOffset = nBaseTypeNum*(pMotif->pPWM->nHeight-1);
		pEle3 = pMCount->pMatElement;
		pEle1 = pMotif->pSampleCount->pMatElement+nOffset;
		pEle2 = pMotif->pPriorCount->pMatElement+nOffset;

		dTotal = 0.0;
		dTotal2 = 0.0;
		for(nj=0; nj<nBaseTypeNum; nj++)
		{
			if(*pEle1 != *pEle3)
			{
				printf("Error: FlexModule_DecreaseWMotifLogLike, motif sample count not match!\n");
				exit(EXIT_FAILURE);
			}
			dTemp = (*pEle1)+(*pEle2);
			dTotal += dTemp;
			dTotal2 += (*pEle2);
			dL += gammaln(dTemp);
			dL -= gammaln(*pEle2);
			pEle1++;
			pEle2++;
			pEle3++;
		}
		dL -= gammaln(dTotal);
		dL += gammaln(dTotal2);
	}
	/* if shrinking the head */
	else if(nDirec == -1)
	{
		pEle3 = pMCount->pMatElement;
		pEle1 = pMotif->pSampleCount->pMatElement;
		pEle2 = pMotif->pPriorCount->pMatElement;

		dTotal = 0.0;
		dTotal2 = 0.0;
		for(nj=0; nj<nBaseTypeNum; nj++)
		{
			if(*pEle1 != *pEle3)
			{
				printf("Error: FlexModule_DecreaseWMotifLogLike, motif sample count not match!\n");
				exit(EXIT_FAILURE);
			}
			dTemp = (*pEle1)+(*pEle2);
			dTotal += dTemp;
			dTotal2 += (*pEle2);
			dL += gammaln(dTemp);
			dL -= gammaln(*pEle2);
			pEle1++;
			pEle2++;
			pEle3++;
		}
		dL -= gammaln(dTotal);
		dL += gammaln(dTotal2);
	}
	/* otherwise error */
	else
	{
		printf("Error: FlexModule_DecreaseWMotifLogLike, the direction of shrinking motif length is not clear!\n");
		exit(EXIT_FAILURE);
	}

	/* return */
	return dL;
}

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_DecreaseWUpdate: update all indicators and context after    */
/*  decreasing motif length by 1.                                          */
/* ----------------------------------------------------------------------- */ 
int FlexModule_DecreaseWUpdate(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, 
				int nMotifId, int nDirec,
				int nModuleLen, struct INTMATRIX *pMotifLen, int *pMaxMotifLen)
{
	/* define */
	unsigned char *pBase,*pAi,*pSi;
	struct FLEXMOTIFSITE **vSite;
	struct FLEXMOTIFSITE *pSite;
	int nLen,nMId,nMotifLen;
	int ni,nSiteNum,nOldPos,nNewPos,nContextNum;
	int nHalfWin;
	struct FLEXMOTIFSITE *pHead,*pTail,*pTest;
	

	/* init */
	nHalfWin = (int)(nModuleLen/2);
	nMId = nMotifId-1;
	nMotifLen = pMotifLen->pMatElement[nMId];
	nLen = pSeqMtf->vSeq[nIndex]->nWidth;
	pBase = pSeqMtf->vSeq[nIndex]->pMatElement;
	pAi = pSeqMtf->vStatus[0]->pMatElement;
	pSi = pSeqMtf->vStatus[1]->pMatElement;
	vSite = pSeqMtf->vMotif;

	/* scan all sites */
	pSite = pSeqMtf->pMotifList;
	nSiteNum = 0;
	while(pSite != NULL)
	{
		/* if motif in question, update */
		if(pSite->nMotifType == nMotifId)
		{
			/* if - strand */
			if(pSite->nStrand == 1)
			{
				nOldPos = pSite->nStartPos;

				if(nDirec == 1)
				{
					nNewPos = nOldPos+1;
				
					/* update Ai */
					pAi[nNewPos] = pAi[nOldPos];
					pAi[nOldPos] = 0;

					/* update Si */
					pSi[nOldPos] = 0;
					
					/* update vSite */
					vSite[nNewPos] = vSite[nOldPos];
					vSite[nOldPos] = NULL;

					/* update site itself */
					pSite->nStartPos = nNewPos;
				}
				else if(nDirec == -1)
				{
					/* update Ai */
					
					/* update Si */
					pSi[nOldPos+nMotifLen-1] = 0;
					
					/* update vSite */
					
					/* update site itself */
				}
				else
				{
					printf("Error: FlexModule_DecreaseWUpdate, direction of extending motif length is not clear!\n");
					exit(EXIT_FAILURE);
				}
			}
			/* if + strand */
			else
			{
				nOldPos = pSite->nStartPos;

				if(nDirec == 1)
				{
					/* update Ai */
					
					/* update Si */
					pSi[nOldPos+nMotifLen-1] = 0;
					
					/* update vSite */
					
					/* update site itself */
				}
				else if(nDirec == -1)
				{
					nNewPos = nOldPos+1;
				
					/* update Ai */
					pAi[nNewPos] = pAi[nOldPos];
					pAi[nOldPos] = 0;

					/* update Si */
					pSi[nOldPos] = 0;
					
					/* update vSite */
					vSite[nNewPos] = vSite[nOldPos];
					vSite[nOldPos] = NULL;

					/* update site itself */
					pSite->nStartPos = nNewPos;
				}
				else
				{
					printf("Error: FlexModule_DecreaseWUpdate, direction of extending motif length is not clear!\n");
					exit(EXIT_FAILURE);
				}
			}
		}

		/* get next */
		nSiteNum++;
		pSite = pSite->pNext;
	}

	if(nSiteNum != pSeqMtf->nSiteNum)
	{
		printf("Error: FlexModule_DecreaseWUpdate, site number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* update context */
	pSite = pSeqMtf->pMotifList;
	while(pSite != NULL)
	{
		/* establish context limit */
		nContextNum = 1;
		pHead = pSite;
		pTail = pSite;

		pTest = pHead->pPrev;
		while(pTest != NULL)
		{
			if( (pSite->nStartPos-pTest->nStartPos) > nHalfWin )
				break;

			pHead = pTest;
			pTest = pHead->pPrev;
			nContextNum++;
		}

		pTest = pTail->pNext;
		while(pTest != NULL)
		{
			if( (pTest->nStartPos-pSite->nStartPos) > nHalfWin )
				break;

			pTail = pTest;
			pTest = pTail->pNext;
			nContextNum++;
		}

		/* create context matrix */
		DestroyIntMatrix(pSite->pContextPos);
		pSite->pContextPos = NULL;
		pSite->pContextPos = CreateIntMatrix(1, nContextNum);
		if(pSite->pContextPos == NULL)
		{
			printf("Error: FlexModule_DecreaseWUpdate, cannot create context matrix correctly!\n");
			exit(EXIT_FAILURE);
		}
		DestroyIntMatrix(pSite->pContextMotif);
		pSite->pContextMotif = NULL;
		pSite->pContextMotif = CreateIntMatrix(1, nContextNum);
		if(pSite->pContextMotif == NULL)
		{
			printf("Error: FlexModule_DecreaseWUpdate, cannot create context matrix correctly!\n");
			exit(EXIT_FAILURE);
		}
		DestroyIntMatrix(pSite->pContextStrand);
		pSite->pContextStrand = NULL;
		pSite->pContextStrand = CreateIntMatrix(1, nContextNum);
		if(pSite->pContextStrand == NULL)
		{
			printf("Error: FlexModule_DecreaseWUpdate, cannot create context matrix correctly!\n");
			exit(EXIT_FAILURE);
		}

		/* fill out the context */
		pTest = pHead;
		ni = 0;
		while(pTest != NULL)
		{
			if(pTest->nStartPos - pTail->nStartPos > 0)
				break;

			pSite->pContextPos->pMatElement[ni] = pTest->nStartPos-pSite->nStartPos;
			pSite->pContextMotif->pMatElement[ni] = pTest->nMotifType;
			pSite->pContextStrand->pMatElement[ni] = pTest->nStrand;
			ni++;
			pTest = pTest->pNext;
		}

		if(ni != nContextNum)
		{
			printf("Error: FlexModule_DecreaseWUpdate, number of motif sites in the context do not match!\n");
			exit(EXIT_FAILURE);
		}

		/* get next */
		pSite = pSite->pNext;
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_DecreaseWUpdateMotif: update motif after decreasing its     */
/*  length by one.                                                         */
/* ----------------------------------------------------------------------- */ 
int FlexModule_DecreaseWUpdateMotif(struct FLEXMOTIFMATRIX *pMotif, 
							 struct DOUBLEMATRIX *pMCount, int nDirec,
							 double dDefaultPriorCount)
{
	/* define */
	double dL = 0.0;
	double *pEle1,*pEle2;
	int nMotifLen,nBaseTypeNum;
	struct DOUBLEMATRIX *pNewPriorCount;
	struct DOUBLEMATRIX *pNewSampleCount;
	struct DOUBLEMATRIX *pNewExtraPrior;

	/* init */
	nMotifLen = pMotif->pPWM->nHeight-1;
	nBaseTypeNum = pMotif->pPWM->nWidth;

	/* get loglikelihood */
	/* if shrinking the tail */
	if(nDirec == 1)
	{
		/* update extra prior pool */
		pMotif->nTailExtraNum += 1;
		if(pMotif->nTailExtraNum > 1)
		{
			pNewExtraPrior = NULL;
			pNewExtraPrior = CreateDoubleMatrix(pMotif->nTailExtraNum, nBaseTypeNum);
			if(pNewExtraPrior == NULL)
			{
				printf("Error: FlexModule_DecreaseWUpdateMotif, cannot create new extra prior count!\n");
				exit(EXIT_FAILURE);
			}

			pEle1 = pNewExtraPrior->pMatElement;
			pEle2 = pMotif->pPriorTailExtra->pMatElement;
			memcpy(pEle1, pEle2, (pMotif->nTailExtraNum-1)*nBaseTypeNum*sizeof(double));

			pEle1 = pNewExtraPrior->pMatElement+(pMotif->nTailExtraNum-1)*nBaseTypeNum;
			pEle2 = pMotif->pPriorCount->pMatElement+nMotifLen*nBaseTypeNum;
			memcpy(pEle1, pEle2, nBaseTypeNum*sizeof(double));

			DestroyDoubleMatrix(pMotif->pPriorTailExtra);
			pMotif->pPriorTailExtra = pNewExtraPrior;
		}
		else if(pMotif->nTailExtraNum == 1)
		{
			if(pMotif->pPriorTailExtra != NULL)
			{
				printf("Error: FlexModule_DecreaseWUpdateMotif, wrong extra prior count!\n");
				exit(EXIT_FAILURE);
			}

			pNewExtraPrior = NULL;
			pNewExtraPrior = CreateDoubleMatrix(1, nBaseTypeNum);
			if(pNewExtraPrior == NULL)
			{
				printf("Error: FlexModule_DecreaseWUpdateMotif, cannot create new extra prior count!\n");
				exit(EXIT_FAILURE);
			}

			pEle1 = pNewExtraPrior->pMatElement;
			pEle2 = pMotif->pPriorCount->pMatElement+nMotifLen*nBaseTypeNum;
			memcpy(pEle1, pEle2, nBaseTypeNum*sizeof(double));

			pMotif->pPriorTailExtra = pNewExtraPrior;
		}
		else
		{
		}

		/* copy the old counts */
		pNewPriorCount = NULL;
		pNewPriorCount = CreateDoubleMatrix(nMotifLen, nBaseTypeNum);
		if(pNewPriorCount == NULL)
		{
			printf("Error: FlexModule_DecreaseWUpdateMotif, cannot create new prior count!\n");
			exit(EXIT_FAILURE);
		}
		pNewSampleCount = NULL;
		pNewSampleCount = CreateDoubleMatrix(nMotifLen, nBaseTypeNum);
		if(pNewSampleCount == NULL)
		{
			printf("Error: FlexModule_DecreaseWUpdateMotif, cannot create new sample count!\n");
			exit(EXIT_FAILURE);
		}

		pEle1 = pNewPriorCount->pMatElement;
		pEle2 = pMotif->pPriorCount->pMatElement;
		memcpy(pEle1, pEle2, nMotifLen*nBaseTypeNum*sizeof(double));

		pEle1 = pNewSampleCount->pMatElement;
		pEle2 = pMotif->pSampleCount->pMatElement;
		memcpy(pEle1, pEle2, nMotifLen*nBaseTypeNum*sizeof(double));

		
		/* update PWM */
		DestroyDoubleMatrix(pMotif->pPriorCount);
		pMotif->pPriorCount = pNewPriorCount;
		DestroyDoubleMatrix(pMotif->pSampleCount);
		pMotif->pSampleCount = pNewSampleCount;
		DestroyDoubleMatrix(pMotif->pPWM);
		pMotif->pPWM = NULL;
		pMotif->pPWM = CreateDoubleMatrix(nMotifLen, nBaseTypeNum);
		if(pMotif->pPWM == NULL)
		{
			printf("Error: FlexModule_DecreaseWUpdateMotif, cannot create new PWM count!\n");
			exit(EXIT_FAILURE);
		}
		FLEXMTFMREFRESH(pMotif);
	}
	/* if extending the head */
	else if(nDirec == -1)
	{
		/* update extra prior pool */
		pMotif->nHeadExtraNum += 1;
		if(pMotif->nHeadExtraNum > 1)
		{
			pNewExtraPrior = NULL;
			pNewExtraPrior = CreateDoubleMatrix(pMotif->nHeadExtraNum, nBaseTypeNum);
			if(pNewExtraPrior == NULL)
			{
				printf("Error: FlexModule_DecreaseWUpdateMotif, cannot create new extra prior count!\n");
				exit(EXIT_FAILURE);
			}

			pEle1 = pNewExtraPrior->pMatElement;
			pEle2 = pMotif->pPriorHeadExtra->pMatElement;
			memcpy(pEle1, pEle2, (pMotif->nHeadExtraNum-1)*nBaseTypeNum*sizeof(double));

			pEle1 = pNewExtraPrior->pMatElement+(pMotif->nHeadExtraNum-1)*nBaseTypeNum;
			pEle2 = pMotif->pPriorCount->pMatElement;
			memcpy(pEle1, pEle2, nBaseTypeNum*sizeof(double));

			DestroyDoubleMatrix(pMotif->pPriorHeadExtra);
			pMotif->pPriorHeadExtra = pNewExtraPrior;
		}
		else if(pMotif->nHeadExtraNum == 1)
		{
			if(pMotif->pPriorHeadExtra != NULL)
			{
				printf("Error: FlexModule_DecreaseWUpdateMotif, wrong extra prior count!\n");
				exit(EXIT_FAILURE);
			}

			pNewExtraPrior = NULL;
			pNewExtraPrior = CreateDoubleMatrix(1, nBaseTypeNum);
			if(pNewExtraPrior == NULL)
			{
				printf("Error: FlexModule_DecreaseWUpdateMotif, cannot create new extra prior count!\n");
				exit(EXIT_FAILURE);
			}

			pEle1 = pNewExtraPrior->pMatElement;
			pEle2 = pMotif->pPriorCount->pMatElement;
			memcpy(pEle1, pEle2, nBaseTypeNum*sizeof(double));

			pMotif->pPriorHeadExtra = pNewExtraPrior;
		}
		else
		{
		}

		/* copy the old counts */
		pNewPriorCount = NULL;
		pNewPriorCount = CreateDoubleMatrix(nMotifLen, nBaseTypeNum);
		if(pNewPriorCount == NULL)
		{
			printf("Error: FlexModule_DecreaseWUpdateMotif, cannot create new prior count!\n");
			exit(EXIT_FAILURE);
		}
		pNewSampleCount = NULL;
		pNewSampleCount = CreateDoubleMatrix(nMotifLen, nBaseTypeNum);
		if(pNewSampleCount == NULL)
		{
			printf("Error: FlexModule_DecreaseWUpdateMotif, cannot create new sample count!\n");
			exit(EXIT_FAILURE);
		}

		pEle1 = pNewPriorCount->pMatElement;
		pEle2 = pMotif->pPriorCount->pMatElement+nBaseTypeNum;
		memcpy(pEle1, pEle2, nMotifLen*nBaseTypeNum*sizeof(double));

		pEle1 = pNewSampleCount->pMatElement;
		pEle2 = pMotif->pSampleCount->pMatElement+nBaseTypeNum;
		memcpy(pEle1, pEle2, nMotifLen*nBaseTypeNum*sizeof(double));

		
		/* update PWM */
		DestroyDoubleMatrix(pMotif->pPriorCount);
		pMotif->pPriorCount = pNewPriorCount;
		DestroyDoubleMatrix(pMotif->pSampleCount);
		pMotif->pSampleCount = pNewSampleCount;
		DestroyDoubleMatrix(pMotif->pPWM);
		pMotif->pPWM = NULL;
		pMotif->pPWM = CreateDoubleMatrix(nMotifLen, nBaseTypeNum);
		if(pMotif->pPWM == NULL)
		{
			printf("Error: FlexModule_DecreaseWUpdateMotif, cannot create new PWM count!\n");
			exit(EXIT_FAILURE);
		}
		FLEXMTFMREFRESH(pMotif);
	}
	/* otherwise error */
	else
	{
		printf("Error: FlexModule_DecreaseWUpdateMotif, the direction of extending motif length is not clear!\n");
		exit(EXIT_FAILURE);
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_ClearOldStatus: clear old status of a given position.       */
/* ----------------------------------------------------------------------- */ 
int FlexModule_ClearOldStatus(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, int nPos, 
			struct FLEXMOTIFSITE *pContext, int *pContextNum, int nBaseTypeNum,
			struct FLEXMOTIFSITE **pPrev, struct FLEXMOTIFSITE **pNext, 
			struct DOUBLEMATRIX *pFreqCount, int nHalfWin, 
			int nMotifNum, struct FLEXMOTIFMATRIX **vMotif, 
			struct INTMATRIX *pMotifLen, int nMaxMotifLen)
{
	/* define */
	int ni,nj,nLen,nRelativePos,nx;
	unsigned char *pBase,*pAi,*pSi;
	struct FLEXMOTIFSITE **vSite;
	struct FLEXMOTIFSITE *pSite;
	int nMotifId,nMotifLen,nMotifStrand;
	double dTemp;

	int nRemoveOld,nAddNew;
	struct INTMATRIX *pContextPos;
	struct INTMATRIX *pContextMotif;
	struct INTMATRIX *pContextStrand;

	/* init */
	nLen = pSeqMtf->vSeq[nIndex]->nWidth;
	pBase = pSeqMtf->vSeq[nIndex]->pMatElement;
	pAi = pSeqMtf->vStatus[0]->pMatElement;
	pSi = pSeqMtf->vStatus[1]->pMatElement;
	vSite = pSeqMtf->vMotif;

	/* adjust context position */
	nRemoveOld = 0;
	if(*pContextNum > 0)
	{
		if(pContext->pContextPos->nWidth != *pContextNum)
		{
			printf("Error: FlexModule_ClearOldStatus, context not synchronized!\n");
			exit(EXIT_FAILURE);
		}

		for(ni=0; ni<(*pContextNum); ni++)
		{
			pContext->pContextPos->pMatElement[ni] -= 1;
		}

		/* remove the old one */
		if( pContext->pContextPos->pMatElement[0] < (-nHalfWin) )
		{
			*pContextNum -= 1;
			nRemoveOld = 1;
		}
	}
	/* add the new one */
	nAddNew = 0;
	if(nPos+nHalfWin < nLen)
	{
		if( vSite[nPos+nHalfWin] != NULL )
		{
			*pContextNum += 1;
			nAddNew = 1;
		}
	}

	if((nRemoveOld == 1) || (nAddNew == 1))
	{
		if(*pContextNum == 0)
		{
			DestroyIntMatrix(pContext->pContextPos);
			pContext->pContextPos = NULL;
			DestroyIntMatrix(pContext->pContextMotif);
			pContext->pContextMotif = NULL;
			DestroyIntMatrix(pContext->pContextStrand);
			pContext->pContextStrand = NULL;
		}
		else if((*pContextNum) > 0)
		{
			pContextPos = NULL;
			pContextPos = CreateIntMatrix(1, (*pContextNum));
			if(pContextPos == NULL)
			{
				printf("Error: FlexModule_ClearOldStatus, cannot create context matrix correctly!\n");
				exit(EXIT_FAILURE);
			}
			pContextMotif = NULL;
			pContextMotif = CreateIntMatrix(1, (*pContextNum));
			if(pContextMotif == NULL)
			{
				printf("Error: FlexModule_ClearOldStatus, cannot create context matrix correctly!\n");
				exit(EXIT_FAILURE);
			}
			pContextStrand = NULL;
			pContextStrand = CreateIntMatrix(1, (*pContextNum));
			if(pContextStrand == NULL)
			{
				printf("Error: FlexModule_ClearOldStatus, cannot create context matrix correctly!\n");
				exit(EXIT_FAILURE);
			}

			/* update new context */
			if( (nRemoveOld == 1) && (nAddNew == 1) )
			{
				nx = (*pContextNum)-1;
				if(nx > 0)
				{
					memcpy(pContextPos->pMatElement, pContext->pContextPos->pMatElement+1, nx*sizeof(int));
					memcpy(pContextMotif->pMatElement, pContext->pContextMotif->pMatElement+1, nx*sizeof(int));
					memcpy(pContextStrand->pMatElement, pContext->pContextStrand->pMatElement+1, nx*sizeof(int));
				}
				pContextPos->pMatElement[nx] = nHalfWin;
				pContextMotif->pMatElement[nx] = vSite[nPos+nHalfWin]->nMotifType;
				pContextStrand->pMatElement[nx] = vSite[nPos+nHalfWin]->nStrand;
			}
			else if(nRemoveOld == 1)
			{
				memcpy(pContextPos->pMatElement, pContext->pContextPos->pMatElement+1, (*pContextNum)*sizeof(int));
				memcpy(pContextMotif->pMatElement, pContext->pContextMotif->pMatElement+1, (*pContextNum)*sizeof(int));
				memcpy(pContextStrand->pMatElement, pContext->pContextStrand->pMatElement+1, (*pContextNum)*sizeof(int));
			}
			else if(nAddNew == 1)
			{
				nx = (*pContextNum)-1;
				if(nx > 0)
				{
					memcpy(pContextPos->pMatElement, pContext->pContextPos->pMatElement, nx*sizeof(int));
					memcpy(pContextMotif->pMatElement, pContext->pContextMotif->pMatElement, nx*sizeof(int));
					memcpy(pContextStrand->pMatElement, pContext->pContextStrand->pMatElement, nx*sizeof(int));
				}
				pContextPos->pMatElement[nx] = nHalfWin;
				pContextMotif->pMatElement[nx] = vSite[nPos+nHalfWin]->nMotifType;
				pContextStrand->pMatElement[nx] = vSite[nPos+nHalfWin]->nStrand;
			}
			else
			{
				printf("Error: FlexModule_ClearOldStatus, logic error!\n");
				exit(EXIT_FAILURE);
			}

			/* destroy old context */
			DestroyIntMatrix(pContext->pContextPos);
			pContext->pContextPos = pContextPos;
			DestroyIntMatrix(pContext->pContextMotif);
			pContext->pContextMotif = pContextMotif;
			DestroyIntMatrix(pContext->pContextStrand);
			pContext->pContextStrand = pContextStrand;
		}
		else
		{
			printf("Error: FlexModule_ClearOldStatus, negative context count!\n");
			exit(EXIT_FAILURE);
		}
	}


	/* check the position in question ! */
	/* if repeat or 'N', skip it */
	if(pBase[nPos] >= nBaseTypeNum)
	{
		/* return */
		return PROC_SUCCESS;
	}

	/* if effective base */

	/* if Bi>0, subtract motif count, motif matrix, change context */
	if(vSite[nPos] != NULL)
	{
		/* get pointer */
		pSite = vSite[nPos];
		nMotifId = pSite->nMotifType;
		nMotifLen = pMotifLen->pMatElement[nMotifId-1];
		nMotifStrand = pSite->nStrand;

		/* subtract frequency */
		dTemp = DMGETAT(pFreqCount, nMotifId, nMotifStrand) - 1.0;
		DMSETAT(pFreqCount, nMotifId, nMotifStrand, dTemp);

		/* if Ai == 1, subtract from PWM count */
		if(pAi[nPos] != 0)
		{
			pAi[nPos] = 0;
			if(nMotifStrand == 1)
			{
				FlexModule_SubtractMotifMatrixCount(vMotif[nMotifId-1]->pSampleCount, 
						  pBase+nPos, nMotifLen, '-');
			}
			else
			{
				FlexModule_SubtractMotifMatrixCount(vMotif[nMotifId-1]->pSampleCount, 
						  pBase+nPos, nMotifLen, '+');
			}

			FLEXMTFMREFRESH(vMotif[nMotifId-1]);
		}

		/* clear site vSite */
		vSite[nPos] = NULL;

		/* clear coverage Si */
		for(nj=0; nj<nMotifLen; nj++)
		{
			if(nPos+nj >= nLen)
			{
				printf("Error: FlexModule_ClearOldStatus, base index out of range!\n");
				exit(EXIT_FAILURE);
			}
			pSi[nPos+nj] = 0;
		}

		/* reconnect the motif site list */
		if( (*pNext != pSite) || (*pPrev != pSite->pPrev) )
		{
			printf("Error: FlexModule_ClearOldStatus, motif list and context not synchronized!\n");
			exit(EXIT_FAILURE);
		}
		*pNext = pSite->pNext;

		if((*pPrev == NULL) && (*pNext == NULL))
		{
			pSeqMtf->pMotifList = NULL;
		}
		else if(*pPrev == NULL)
		{
			(*pNext)->pPrev = *pPrev;
			pSeqMtf->pMotifList = *pNext;
		}
		else if(*pNext == NULL)
		{
			(*pPrev)->pNext = *pNext;
		}
		else
		{
			(*pPrev)->pNext = *pNext;
			(*pNext)->pPrev = *pPrev;
		}

		/* update pContext */
		nx = (*pContextNum)-1;
		if(nx == 0)
		{
			DestroyIntMatrix(pContext->pContextPos);
			pContext->pContextPos = NULL;
			DestroyIntMatrix(pContext->pContextMotif);
			pContext->pContextMotif = NULL;
			DestroyIntMatrix(pContext->pContextStrand);
			pContext->pContextStrand = NULL;
		}
		else
		{
			pContextPos = NULL;
			pContextPos = CreateIntMatrix(1, nx);
			if(pContextPos == NULL)
			{
				printf("Error: FlexModule_ClearOldStatus, cannot create context matrix correctly!\n");
				exit(EXIT_FAILURE);
			}
			pContextMotif = NULL;
			pContextMotif = CreateIntMatrix(1, nx);
			if(pContextMotif == NULL)
			{
				printf("Error: FlexModule_ClearOldStatus, cannot create context matrix correctly!\n");
				exit(EXIT_FAILURE);
			}
			pContextStrand = NULL;
			pContextStrand = CreateIntMatrix(1, nx);
			if(pContextStrand == NULL)
			{
				printf("Error: FlexModule_ClearOldStatus, cannot create context matrix correctly!\n");
				exit(EXIT_FAILURE);
			}

			nj = 0;
			for(ni=0; ni<(*pContextNum); ni++)
			{
				if(pContext->pContextPos->pMatElement[ni] == 0)
					continue;

				pContextPos->pMatElement[nj] = pContext->pContextPos->pMatElement[ni];
				pContextMotif->pMatElement[nj] = pContext->pContextMotif->pMatElement[ni];
				pContextStrand->pMatElement[nj] = pContext->pContextStrand->pMatElement[ni];
				nj++;
			}

			DestroyIntMatrix(pContext->pContextPos);
			pContext->pContextPos = pContextPos;
			DestroyIntMatrix(pContext->pContextMotif);
			pContext->pContextMotif = pContextMotif;
			DestroyIntMatrix(pContext->pContextStrand);
			pContext->pContextStrand = pContextStrand;
		}
		
		/* update other sites in the context */
		for(ni=0; ni<nx; ni++)
		{
			nj = pContext->pContextPos->pMatElement[ni];
			nRelativePos = -nj;
			if(vSite[nPos+nj] == NULL)
			{
				printf("Error: FlexModule_ClearOldStatus, context not synchronized!\n");
				exit(EXIT_FAILURE);
			}

			FLEXMTFSREMOVESITE(vSite[nPos+nj], nRelativePos);
		}
		
		/* destroy the old site */
		FLEXMTFSDESTROY(pSite);
		pSeqMtf->nSiteNum -= 1;
		*pContextNum = nx;
	}
	/* if Bi=0, subtract background count */
	else
	{
		if(pAi[nPos] != 0)
		{
			printf("Error: FlexModule_ClearOldStatus, motif indicator wrong !\n");
			exit(EXIT_FAILURE);
		}

		/* subtract frequency */
		dTemp = DMGETAT(pFreqCount, 0, 0) - 1.0;
		DMSETAT(pFreqCount, 0, 0, dTemp);

		/* there is no need to update pPrev, pNext, pContext */		
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_UpdateNewStatus: update new context status.                 */
/* ----------------------------------------------------------------------- */ 
int FlexModule_UpdateNewStatus(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, int nPos,
			int nMotifId, int nStrand,
			struct FLEXMOTIFSITE *pContext, int *pContextNum, int nBaseTypeNum,
			struct FLEXMOTIFSITE **pPrev, struct FLEXMOTIFSITE **pNext, 
			int nHalfWin, 
			int nMotifNum, struct FLEXMOTIFMATRIX **vMotif, 
			struct INTMATRIX *pMotifLen, int nMaxMotifLen)
{
	/* define */
	int ni,nj,nRelativePos,nNum,nLen;
	int nAdd;
	int nMId;
	int nMotifLen;
	unsigned char *vSi;
	struct FLEXMOTIFSITE *pSite;
	struct FLEXMOTIFSITE **vSite;
	struct INTMATRIX *pPos, *pMotif, *pStrand;

	/* init */
	nMId = nMotifId-1;
	nMotifLen = pMotifLen->pMatElement[nMId];
	vSi = pSeqMtf->vStatus[1]->pMatElement;
	vSite = pSeqMtf->vMotif;
	nLen = pSeqMtf->vSeq[nIndex]->nWidth;

	/* update context sites */
	for(ni=0; ni<(*pContextNum); ni++)
	{
		nRelativePos = pContext->pContextPos->pMatElement[ni];
		pSite = vSite[nPos+nRelativePos];
		if(pSite == NULL)
		{
			printf("Error: FlexModule_UpdateNewStatus, context not synchronized!\n");
			exit(EXIT_FAILURE);
		}

		FLEXMTFSADDSITE(pSite, -nRelativePos, nMotifId, nStrand);
	}

	/* update context number */
	nNum = *pContextNum;
	*pContextNum = (*pContextNum)+1;

	/* add new site */
	pPos = NULL;
	pPos = CreateIntMatrix(1, *pContextNum);
	if(pPos == NULL)
	{
		printf("Error: FlexModule_UpdateNewStatus, cannot add context sites!\n");
		exit(EXIT_FAILURE);
	}

	pMotif = NULL;
	pMotif = CreateIntMatrix(1, *pContextNum);
	if(pMotif == NULL)
	{
		printf("Error: FlexModule_UpdateNewStatus, cannot add context sites!\n");
		exit(EXIT_FAILURE);
	}

	pStrand = NULL;
	pStrand = CreateIntMatrix(1, *pContextNum);
	if(pStrand == NULL)
	{
		printf("Error: FlexModule_UpdateNewStatus, cannot add context sites!\n");
		exit(EXIT_FAILURE);
	}

	nAdd = 0;
	nj = 0;
	for(ni=0; ni<nNum; ni++)
	{
		if((pContext->pContextPos->pMatElement[ni] > 0) && (nAdd == 0))
		{
			pPos->pMatElement[nj] = 0;
			pMotif->pMatElement[nj] = nMotifId;
			pStrand->pMatElement[nj] = nStrand;
			nj++;
			nAdd = 1;
		}
		pPos->pMatElement[nj] = pContext->pContextPos->pMatElement[ni];
		pMotif->pMatElement[nj] = pContext->pContextMotif->pMatElement[ni];
		pStrand->pMatElement[nj] = pContext->pContextStrand->pMatElement[ni];
		nj++;
	}
	if(nAdd == 0)
	{
		pPos->pMatElement[nj] = 0;
		pMotif->pMatElement[nj] = nMotifId;
		pStrand->pMatElement[nj] = nStrand;
	}

	DestroyIntMatrix(pContext->pContextPos);
	pContext->pContextPos = pPos;
	DestroyIntMatrix(pContext->pContextMotif);
	pContext->pContextMotif = pMotif;
	DestroyIntMatrix(pContext->pContextStrand);
	pContext->pContextStrand = pStrand;

	/* create new site */
	pSite = NULL;
	pSite = FLEXMTFSCREATE();
	if(pSite == NULL)
	{
		printf("Error: FlexModule_UpdateNewStatus, cannot create new site!\n");
		exit(EXIT_FAILURE);
	}
	pSite->nMotifType = nMotifId;
	pSite->nSeqId = pSeqMtf->nId;
	pSite->nStartPos = nPos;
	pSite->nStrand = nStrand;
	pSite->pContextPos = NULL;
	pSite->pContextPos = IMCLONE(pContext->pContextPos);
	if(pSite->pContextPos == NULL)
	{
		printf("Error: FlexModule_UpdateNewStatus, cannot create new site!\n");
		exit(EXIT_FAILURE);
	}
	pSite->pContextMotif = NULL;
	pSite->pContextMotif = IMCLONE(pContext->pContextMotif);
	if(pSite->pContextMotif == NULL)
	{
		printf("Error: FlexModule_UpdateNewStatus, cannot create new site!\n");
		exit(EXIT_FAILURE);
	}
	pSite->pContextStrand = NULL;
	pSite->pContextStrand = IMCLONE(pContext->pContextStrand);
	if(pSite->pContextStrand == NULL)
	{
		printf("Error: FlexModule_UpdateNewStatus, cannot create new site!\n");
		exit(EXIT_FAILURE);
	}
	pSite->pNext = NULL;
	pSite->pPrev = NULL;

	/* add the new site to the current list */
	vSite[nPos] = pSite;

	if( (*pPrev == NULL) && (*pNext == NULL) )
	{
		pSeqMtf->pMotifList = pSite;
		*pPrev = pSite;
	}
	else if(*pPrev == NULL)
	{
		pSite->pNext = *pNext;
		(*pNext)->pPrev = pSite;
		pSeqMtf->pMotifList = pSite;
		*pPrev = pSite;
	}
	else if(*pNext == NULL)
	{
		pSite->pPrev = *pPrev;
		(*pPrev)->pNext = pSite;
		*pPrev = pSite;
	}
	else
	{
		(*pPrev)->pNext = pSite;
		pSite->pPrev = (*pPrev);
		(*pNext)->pPrev = pSite;
		pSite->pNext = *pNext;
		*pPrev = pSite;
	}

	/* update vSi */
	for(ni=0; ni<nMotifLen; ni++)
	{
		if( (nPos+ni) >= nLen)
		{
			printf("Error: FlexModule_UpdateNewStatus, base index out of range!\n");
			exit(EXIT_FAILURE);
		}
		vSi[nPos+ni] = nMotifId;
	}

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  FlexModule_ComputeContextLogLike: compute context log likelihood.      */
/* ----------------------------------------------------------------------- */ 
double FlexModule_ComputeContextLogLike(struct FLEXSEQMOTIF *pSeqMtf, struct FLEXMOTIFSITE *pContext,
					int nPos, int nContextNum, int nAddPos, int nAddMotifId, int nAddStrand,
					int nMotifNum, double dModuleD, double dModuleL, double dModuleA, double dT)
{
	/* define */
	double dLogLike = 0.0;
	double dScore,dPAi;
	int ni,nj;
	struct FLEXMOTIFSITE *pSite;
	struct FLEXMOTIFSITE **vSite;
	unsigned char *vAi;

	/* init */
	vAi = pSeqMtf->vStatus[0]->pMatElement;

	/* compute one by one */
	vSite = pSeqMtf->vMotif;
	for(ni=0; ni<nContextNum; ni++)
	{
		nj = nPos+pContext->pContextPos->pMatElement[ni];
		pSite = vSite[nj];
		if(pSite == NULL)
		{
			printf("Error: FlexModule_ComputeContextLogLike, context not synchronized!\n");
			exit(EXIT_FAILURE);
		}

		/* add new site */
		dScore = FlexModule_ComputeModuleScore(pSite->pContextPos->nWidth,
			pSite->pContextPos->pMatElement, pSite->pContextMotif->pMatElement,
			pSite->pContextStrand->pMatElement, 1, -(pContext->pContextPos->pMatElement[ni]), 
			nAddMotifId, nAddStrand, nMotifNum);
		dPAi = FlexModule_ComputeTrueProb(dScore, dModuleD, dModuleL, dModuleA, dT);
		if(vAi[nj] == 1)
			dLogLike += log(dPAi);
		else
			dLogLike += log(1.0-dPAi);

		/* do not add new site */
		dScore = FlexModule_ComputeModuleScore(pSite->pContextPos->nWidth,
			pSite->pContextPos->pMatElement, pSite->pContextMotif->pMatElement,
			pSite->pContextStrand->pMatElement, 0, 0, 0, 0, nMotifNum);
		dPAi = FlexModule_ComputeTrueProb(dScore, dModuleD, dModuleL, dModuleA, dT);
		if(vAi[nj] == 1)
			dLogLike -= log(dPAi);
		else
			dLogLike -= log(1.0-dPAi);
	}

	/* return */
	return dLogLike;
}


/* ----------------------------------------------------------------------- */ 
/* FlexModule_AddMotifMatrixCount: Add count to motif matrix.              */
/* ----------------------------------------------------------------------- */ 
int FlexModule_AddMotifMatrixCount(struct DOUBLEMATRIX *pCountMat, 
						  unsigned char pSite[], int nMotifLen, char chStrand)
{
	/* define */
	int ni,nj,nk;
	int nBaseTypeNum;
	double dTemp;

	/* check */
	if(pCountMat == NULL)
	{
		printf("Error: FlexModule_AddNucleicMatrixCount, no count matrix exist for counting!\n");
		exit(EXIT_FAILURE);
	}
	if(pCountMat->nHeight != nMotifLen)
	{
		printf("Error: FlexModule_AddNucleicMatrixCount, motif lengths not match!\n");
		exit(EXIT_FAILURE);
	}
	nBaseTypeNum = pCountMat->nWidth;

	/* add count from a '+' strand site */
	if(chStrand == '+')
	{
		for(ni=0; ni<nMotifLen; ni++)
		{
			if(pSite[ni] >= nBaseTypeNum)
				continue;
			nj = (int)pSite[ni];
			dTemp = DMGETAT(pCountMat, ni, nj) + 1.0;
			DMSETAT(pCountMat, ni, nj, dTemp);
		}
	}
	/* add count from a '-' strand site */
	else if(chStrand == '-')
	{
		for(ni=0; ni<nMotifLen; ni++)
		{
			if(pSite[ni] >= nBaseTypeNum)
				continue;
			nj = nBaseTypeNum-1-(int)pSite[ni];
			nk = nMotifLen-1-ni;
			dTemp = DMGETAT(pCountMat, nk, nj) + 1.0;
			DMSETAT(pCountMat, nk, nj, dTemp);
		}
	}
	/* do nothing */
	else
	{
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/* FlexModule_SubtractMotifMatrixCount: Subtract count to motif matrix.    */
/* ----------------------------------------------------------------------- */ 
int FlexModule_SubtractMotifMatrixCount(struct DOUBLEMATRIX *pCountMat, 
						  unsigned char pSite[], int nMotifLen, char chStrand)
{
	/* define */
	int ni,nj,nk;
	int nBaseTypeNum;
	double dTemp;

	/* check */
	if(pCountMat == NULL)
	{
		printf("Error: FlexModule_SubtractMotifMatrixCount, no count matrix exist for counting!\n");
		exit(EXIT_FAILURE);
	}
	if(pCountMat->nHeight != nMotifLen)
	{
		printf("Error: FlexModule_SubtractMotifMatrixCount, motif lengths not match!\n");
		exit(EXIT_FAILURE);
	}
	nBaseTypeNum = pCountMat->nWidth;

	/* add count from a '+' strand site */
	if(chStrand == '+')
	{
		for(ni=0; ni<nMotifLen; ni++)
		{
			if(pSite[ni] >= nBaseTypeNum)
				continue;
			nj = (int)pSite[ni];
			dTemp = DMGETAT(pCountMat, ni, nj) - 1.0;
			DMSETAT(pCountMat, ni, nj, dTemp);
		}
	}
	/* add count from a '-' strand site */
	else if(chStrand == '-')
	{
		for(ni=0; ni<nMotifLen; ni++)
		{
			if(pSite[ni] >= nBaseTypeNum)
				continue;
			nj = nBaseTypeNum-1-(int)pSite[ni];
			nk = nMotifLen-1-ni;
			dTemp = DMGETAT(pCountMat, nk, nj) - 1.0;
			DMSETAT(pCountMat, nk, nj, dTemp);
		}
	}
	/* do nothing */
	else
	{
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_NormalizeFreq: normalize the frequency matrix.              */
/* ----------------------------------------------------------------------- */ 
int FlexModule_NormalizeFreq(struct DOUBLEMATRIX *pFreq)
{
	/* define */
	double *pElement;
	int ni,nj;
	double dTotal;

	/* normalize */
	pElement = pFreq->pMatElement;
	dTotal = *pElement;
	pElement += pFreq->nWidth;
	for(ni=1; ni<pFreq->nHeight; ni++)
	{
		for(nj=0; nj<pFreq->nWidth; nj++)
		{
			dTotal += *pElement;
			pElement++;
		}
	}

	pElement = pFreq->pMatElement;
	*pElement = *pElement/dTotal;
	pElement += pFreq->nWidth;
	for(ni=1; ni<pFreq->nHeight; ni++)
	{
		for(nj=0; nj<pFreq->nWidth; nj++)
		{
			*pElement = (*pElement)/dTotal;
			pElement++;
		}
	}

	/* return */
	return PROC_SUCCESS;
}



/* ----------------------------------------------------------------------- */ 
/*  FlexModule_LogFreq: take log-transformation of the frequency matrix.   */
/* ----------------------------------------------------------------------- */ 
int FlexModule_LogFreq(struct DOUBLEMATRIX *pFreq)
{
	/* define */
	double *pElement;
	int ni,nj;
	
	/* log */
	pElement = pFreq->pMatElement;
	if(*pElement <= 0.0)
	{
		printf("Error: FlexModule_LogFreq, log(0)!\n");
		exit(EXIT_FAILURE);
	}
	*pElement = log(*pElement);
	pElement += pFreq->nWidth;
	for(ni=1; ni<pFreq->nHeight; ni++)
	{
		for(nj=0; nj<pFreq->nWidth; nj++)
		{
			if(*pElement <= 0.0)
			{
				printf("Error: FlexModule_LogFreq, log(0)!\n");
				exit(EXIT_FAILURE);
			}
			*pElement = log(*pElement);
			pElement++;
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  FlexBaseLogLikelihood_Motif: get loglikelihood for motif.              */
/* ----------------------------------------------------------------------- */ 
double FlexBaseLogLikelihood_Motif(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, int nPos, 
					struct DOUBLEMATRIX *pLogPWM, char chStrand, int nMaxMotifLen,
					int nBGOrder, struct DOUBLEMATRIX *pLogBG, struct DOUBLEMATRIX *pLogBG0, 
					int nUseCS, int nCSLikeNum, struct DOUBLEMATRIX **vCSLike,
					int *pMask)
{
	/* define */
	double dL;
	int ni,nMotifLen,nW;
	int nBaseTypeNum;
	int nSeqLen;
	unsigned char *pBase,*pCS,*pSi;
	double *pBL;
	int nCSId;
	int nMidPos;

	/* init */
	*pMask = 0;
	dL = 0.0;
	nBaseTypeNum = pLogBG->nWidth;

	/* check */
	if(pLogPWM == NULL)
	{
		printf("Error: FlexBaseLogLikelihood_Motif, no PWM exist for likelihood calculation!\n");
		exit(EXIT_FAILURE);
	}
	
	/* get likelihood */
	nMotifLen = pLogPWM->nHeight;
	nMidPos = (int)((nMotifLen-1)/2);
	nSeqLen = pSeqMtf->vSeq[nIndex]->nWidth;
	if((nPos+nMotifLen) > nSeqLen)
	{
		*pMask = 1;
		return dL;
	}

	pBase = pSeqMtf->vSeq[nIndex]->pMatElement+nPos;
	if(nUseCS == 1)
	{
		pCS = pSeqMtf->vScore[0]->pMatElement+nPos;
	}
	pSi = pSeqMtf->vStatus[1]->pMatElement+nPos;
	pBL = pSeqMtf->vMonitor[0]->pMatElement+nPos;

	if(chStrand == '+')
	{
		for(ni=0; ni<nMotifLen; ni++)
		{
			nW = (int)(pBase[ni]);
			
			if( (pSi[ni] > 0) || (nW >= nBaseTypeNum) )
			{
				*pMask = 1;
				dL = 0.0;
				return dL;
			}

			dL += DMGETAT(pLogPWM, ni, nW);
			dL -= pBL[ni];

			/* if(dL < -1e10)
			{
				dL = dL;
			} */

		}

		if(nUseCS == 1)
		{	
			nCSId = (int)(pCS[nMidPos]);
			dL += (vCSLike[1]->pMatElement[nCSId]-vCSLike[0]->pMatElement[nCSId]);

			/* if(dL < -1e10)
			{
				dL = dL;
			} */
		}
	}
	else if(chStrand == '-')
	{
		for(ni=0; ni<nMotifLen; ni++)
		{
			nW = (int)(pBase[ni]);
			
			if( (pSi[ni] > 0) || (nW >= nBaseTypeNum) )
			{
				*pMask = 1;
				dL = 0.0;
				return dL;
			}

			nW = nBaseTypeNum-1-nW;
			dL += DMGETAT(pLogPWM, (nMotifLen-1-ni), nW);
			dL -= pBL[ni];

			/* if(dL < -1e10)
			{
				dL = dL;
			} */
		}

		if(nUseCS == 1)
		{	
			nCSId = (int)(pCS[nMidPos]);
			dL += (vCSLike[1]->pMatElement[nCSId]-vCSLike[0]->pMatElement[nCSId]);

			/* if(dL < -1e10)
			{
				dL = dL;
			} */
		}
	}
	else
	{
		printf("Error: FlexBaseLogLikelihood_Motif, strand information wrong!\n");
		exit(EXIT_FAILURE);
	}

	/* return */
	return dL;
}

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_NormalizePost: normalize log frequency.                     */
/* ----------------------------------------------------------------------- */ 
int FlexModule_NormalizePost(struct DOUBLEMATRIX *pFreq, struct BYTEMATRIX *pValid)
{
	/* define */
	int ni,nj;
	double dMax,dSum;
	double *pElement;
	unsigned char *pVi;

	/* get max */
	pElement = pFreq->pMatElement;
	pVi = pValid->pMatElement;
	dMax = *pElement;
	pElement += pFreq->nWidth;
	pVi += pValid->nWidth;
	for(ni=1; ni<pFreq->nHeight; ni++)
	{
		for(nj=0; nj<pFreq->nWidth; nj++)
		{
			if((*pVi == 1) && (*pElement > dMax))
				dMax = *pElement;
			pElement++;
			pVi++;
		}
	}

	/* get exp(x-dmax) */
	dSum = 0.0;
	pElement = pFreq->pMatElement;
	pVi = pValid->pMatElement;
	*pElement = exp(*pElement-dMax);
	dSum += *pElement;
	pElement += pFreq->nWidth;
	pVi += pValid->nWidth;
	for(ni=1; ni<pFreq->nHeight; ni++)
	{
		for(nj=0; nj<pFreq->nWidth; nj++)
		{
			if(*pVi == 1)
			{
				*pElement = exp(*pElement-dMax);
				dSum += *pElement;
			}
			pElement++;
			pVi++;
		}
	}

	/* normalize */
	if(dSum <= 0.0)
	{
		printf("Error: FlexModule_NormalizePost, divide by zero!\n");
		exit(EXIT_FAILURE);
	}
	pElement = pFreq->pMatElement;
	pVi = pValid->pMatElement;
	*pElement /= dSum;
	pElement += pFreq->nWidth;
	pVi += pValid->nWidth;
	for(ni=1; ni<pFreq->nHeight; ni++)
	{
		for(nj=0; nj<pFreq->nWidth; nj++)
		{
			if(*pVi == 1)
			{
				*pElement /= dSum;
			}
			
			pElement++;
			pVi++;
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_SamplePost: sample from a frequency matrix.                 */
/* ----------------------------------------------------------------------- */ 
int FlexModule_SamplePost(struct DOUBLEMATRIX *pFreq, struct BYTEMATRIX *pValid,
						  int *pI, int *pJ)
{
	/* define */
	double dRand;
	double dSum;
	int ni,nj;
	double *pElement;
	unsigned char *pVi;


	/* sample */
	dRand = rand_u();
	
	/* get max */
	pElement = pFreq->pMatElement;
	pVi = pValid->pMatElement;
	dSum = *pElement;
	if(dRand <= dSum)
	{
		*pI = 0;
		*pJ = 0;
		return PROC_SUCCESS;
	}

	pElement += pFreq->nWidth;
	pVi += pValid->nWidth;
	for(ni=1; ni<pFreq->nHeight; ni++)
	{
		for(nj=0; nj<pFreq->nWidth; nj++)
		{
			if(*pVi == 1)
				dSum += *pElement;
			if(dRand <= dSum)
			{
				*pI = ni;
				*pJ = nj;
				return PROC_SUCCESS;
			}
			pElement++;
			pVi++;
		}
	}

	*pI = 0;
	*pJ = 0;
	
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_CallMotifSitesInFinalSample: call motif site according to   */
/*  the last sample.                                                       */
/* ----------------------------------------------------------------------- */ 
int FlexModule_CallMotifSitesInFinalSample(struct FLEXSEQMOTIF *pSeqMtf, 
					struct INTMATRIX *pMotifLen)
{
	/* define */
	int ni,nLen;
	unsigned char *pAi;
	struct FLEXMOTIFSITE *pMtfSite,*pPrevSite;
	struct FLEXMOTIFSITE **vSite;

	/* init */
	if( pSeqMtf->pMotifList != NULL)
	{
		pSeqMtf->pMotifList = NULL;
	}
	pMtfSite = NULL;
	pPrevSite = NULL;

	nLen = pSeqMtf->nSeqLen;
	pAi = pSeqMtf->vStatus[0]->pMatElement;
	vSite = pSeqMtf->vMotif;

	for(ni=0; ni<nLen; ni++)
	{
		if(vSite[ni] != NULL)
		{
			vSite[ni]->pNext = NULL;
			vSite[ni]->pPrev = NULL;

			if(pAi[ni] == 1)
			{			
				if(pPrevSite == NULL)
				{
					pSeqMtf->pMotifList = vSite[ni];
					pPrevSite = vSite[ni];
				}
				else
				{
					pPrevSite->pNext = vSite[ni];
					pPrevSite = vSite[ni];
				}
			}
			else
			{
				FLEXMTFSDESTROY(vSite[ni]);
			}
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_WriteMotifToFile: write motif sites to file.                */
/* ----------------------------------------------------------------------- */ 
int FlexModule_WriteMotifToFile(char strFilePath[], int nMotifNum, 
					struct FLEXMOTIFMATRIX **vMotif, struct DOUBLEMATRIX *pBG0, 
					struct INTMATRIX *pMotifLen,
					int nSeqCount, struct FLEXSEQMOTIF **vSeqMtf, int nIndex)
{
	/* define */
	FILE *fpOut;
	int ni,nj,nk,nSeqId,nMId;
	double *pElement1,*pElement2;
	char strConsensus[LINE_LENGTH];
	struct FLEXMOTIFSITE *pSite;
	unsigned char *pBi;
	double dScore;
	int nMax;
	double dMax;
	int nBaseTypeNum;
	double dTemp;

	/* report results */
	nBaseTypeNum = pBG0->nWidth;

	fpOut = NULL;
	fpOut = fopen(strFilePath, "w");
	if(fpOut == NULL)
	{
		printf("Error: cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	/* write motif one by one */
	for(ni=0; ni<nMotifNum; ni++)
	{
		/* get motif score */
		dScore = FLEXMTFMSCORE(vMotif[ni], pBG0);

		/* write head */
		fprintf(fpOut, "****** Motif%d ******\n", ni);
		fprintf(fpOut, "Motif Score: %f\n", dScore);
		fprintf(fpOut, "Motif Matrix: \n");
		pElement1 = vMotif[ni]->pPriorCount->pMatElement;
		pElement2 = vMotif[ni]->pSampleCount->pMatElement;
		for(nj=0; nj<vMotif[ni]->pPriorCount->nHeight; nj++)
		{
			dMax = 0.0;
			nMax = nBaseTypeNum;
			for(nk=0; nk<vMotif[ni]->pPriorCount->nWidth; nk++)
			{
				dTemp = (*pElement1)+(*pElement2);
				fprintf(fpOut, "% 9.7e ", dTemp);
				if(dTemp > dMax)
				{
					dMax = dTemp;
					nMax = nk;
				}
				pElement1++;
				pElement2++;
			}

			fprintf(fpOut, "\n");
			
			switch(nMax)
			{
				case 0: strConsensus[nj] = 'A';
					break;
				case 1: strConsensus[nj] = 'C';
					break;
				case 2: strConsensus[nj] = 'G';
					break;
				case 3: strConsensus[nj] = 'T';
					break;
				default: strConsensus[nj] = 'N';
			}
		}
		strConsensus[nj] = '\0';

		fprintf(fpOut, "\n");
		fprintf(fpOut, "Consensus:\n");
		fprintf(fpOut, "%s\n", strConsensus);

		fprintf(fpOut, "\nMotif Sites:\n");
		for(nk=0; nk<nSeqCount; nk++)
		{
			pSite = vSeqMtf[nk]->pMotifList;
			while(pSite != NULL)
			{
				/* if match, then write */
				if((pSite->nMotifType-1) != ni)
				{
					pSite = pSite->pNext;
					continue;
				}

				/* head information */
				nSeqId = pSite->nSeqId;
				nMId = pSite->nMotifType-1;
				/* fprintf(fpOut, ">%d:%d-%d", vSeqMtf[nk]->nId, 
						pSite->nStartPos, (pSite->nStartPos+pMotifLen->pMatElement[nMId]-1));
				if(pSite->nStrand == 0)
				{
					fprintf(fpOut, "(+)\n");
				}
				else if(pSite->nStrand == 1)
				{
					fprintf(fpOut, "(-)\n");
				}
				else
				{
					fprintf(fpOut, "(?)\n");
				}*/

				fprintf(fpOut, "%d\t%d\t%d\t", vSeqMtf[nk]->nId, 
					pSite->nStartPos, (pSite->nStartPos+pMotifLen->pMatElement[nMId]-1));
				if(pSite->nStrand == 0)
				{
					fprintf(fpOut, "+\t");
				}
				else if(pSite->nStrand == 1)
				{
					fprintf(fpOut, "-\t");
				}
				else
				{
					fprintf(fpOut, "?\t");
				}

				/* write seq */
				if(pSite->nStrand == 0)
				{
					pBi = vSeqMtf[nk]->vSeq[nIndex]->pMatElement + pSite->nStartPos;
					for(nj=0; nj<pMotifLen->pMatElement[nMId]; nj++)
					{
						switch((int)(*pBi))
						{
							case 0: fprintf(fpOut, "A");
								break;
							case 1: fprintf(fpOut, "C");
								break;
							case 2: fprintf(fpOut, "G");
								break;
							case 3: fprintf(fpOut, "T");
								break;
							case 10: fprintf(fpOut, "a");
								break;
							case 11: fprintf(fpOut, "c");
								break;
							case 12: fprintf(fpOut, "g");
								break;
							case 13: fprintf(fpOut, "t");
								break;
							default: fprintf(fpOut, "N");
						}
						pBi++;
					}
				}
				else if(pSite->nStrand == 1)
				{
					pBi = vSeqMtf[nk]->vSeq[nIndex]->pMatElement + pSite->nStartPos + pMotifLen->pMatElement[nMId] - 1;
					for(nj=0; nj<pMotifLen->pMatElement[nMId]; nj++)
					{
						switch((int)(*pBi))
						{
							case 0: fprintf(fpOut, "T");
								break;
							case 1: fprintf(fpOut, "G");
								break;
							case 2: fprintf(fpOut, "C");
								break;
							case 3: fprintf(fpOut, "A");
								break;
							case 10: fprintf(fpOut, "t");
								break;
							case 11: fprintf(fpOut, "g");
								break;
							case 12: fprintf(fpOut, "c");
								break;
							case 13: fprintf(fpOut, "a");
								break;
							default: fprintf(fpOut, "N");
						}
						pBi--;
					}
				}
				fprintf(fpOut, "\n");

				/* get next */
				pSite = pSite->pNext;
			}
		}

		fprintf(fpOut, "\n\n");
	}

	/* close file */
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;

}

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_ComputeSiteScore: compute motif score for a site.           */
/* ----------------------------------------------------------------------- */ 
double FlexModule_ComputeSiteScore(unsigned char *pBase, 
					struct DOUBLEMATRIX *pPWM, struct DOUBLEMATRIX *pBG, 
					char chStrand)
{
	/* define */
	double dScore = 0.0;
	int nBaseTypeNum;
	int ni,nj,nLen;

	/* compute */
	nBaseTypeNum = pBG->nWidth;
	nLen = pPWM->nHeight;

	if(chStrand == '+')
	{
		for(ni=0; ni<nLen; ni++)
		{
			nj = (int)pBase[ni];
			if(nj < nBaseTypeNum)
			{
				dScore += DMGETAT(pPWM, ni, nj);
				dScore -= DMGETAT(pBG, 0, nj);
			}
		}
	}
	else if(chStrand == '-')
	{
		for(ni=0; ni<nLen; ni++)
		{
			nj = (int)pBase[ni];
			if(nj < nBaseTypeNum)
			{
				nj = nBaseTypeNum-1-nj;
				dScore += DMGETAT(pPWM, (nLen-1-ni), nj);
				dScore -= DMGETAT(pBG, 0, nj);
			}
		}
	}
	else
	{
		printf("Error: FlexModule_ComputeSiteScore, strand not clear!\n");
		exit(EXIT_FAILURE);
	}

	/* return */
	return dScore;
}

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_WritePostMotifToFile: write motif sites to file.            */
/* ----------------------------------------------------------------------- */ 
int FlexModule_WritePostMotifToFile(char strFilePath[], int nMotifNum, 
					struct FLEXMOTIFMATRIX **vMotif, 
					struct DOUBLEMATRIX *pBG0, struct DOUBLEMATRIX *pLogBG0, 
					struct INTMATRIX *pMotifLen, struct INTMATRIX *pMeanMotifLen, 
					struct DOUBLEMATRIX *pMotifNSample, 
					int nSeqCount, struct FLEXSEQMOTIF **vSeqMtf, 
					int nIndex, int nRecordNum, double dDefaultPrior)
{
	/* define */
	FILE *fpOut;
	int nBaseTypeNum;
	int nLen;
	int ni,nj,nk,nx;
	double *vPostSamp;
	char strLine[MED_LINE_LENGTH];
	char strConsensus[LINE_LENGTH];
	struct DOUBLEMATRIX *pCumProb;
	struct DOUBLEMATRIX *pShiftProb;
	double dLowCut = 0.05;
	struct FLEXMOTIFSITE *pMtfSite,*pPrevSite,*pTempSite,*pSite;
	struct FLEXMOTIFSITE **vMotifList;
	double dLike1,dLike2;
	unsigned char *pBase,*pBi;
	double *pMoniter;
	int nMId;
	double dScore;
	int nMax;
	double dMax;
	double dTemp;
	double *pElement1,*pElement2;
	int nSeqId,nPos;


	/* init */
	nBaseTypeNum = pBG0->nWidth;

	if(nRecordNum <= 0)
	{
		printf("Error: FlexModule_WritePostMotifToFile, no posterior samples available!\n");
		exit(EXIT_FAILURE);
	}

	pCumProb = NULL;
	pCumProb = CreateDoubleMatrix(1, nMotifNum);
	if(pCumProb == NULL)
	{
		printf("Error: FlexModule_WritePostMotifToFile, cannot create space for cumulant probability computation!\n");
		exit(EXIT_FAILURE);
	}

	vMotifList = NULL;
	vMotifList = (struct FLEXMOTIFSITE **)calloc(nMotifNum, sizeof(struct FLEXMOTIFSITE *));
	if(vMotifList == NULL)
	{
		printf("Error: FlexModule_WritePostMotifToFile, cannot create motif site list!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nMotifNum; ni++)
	{
		vMotifList[ni] = NULL;
	}

	/* initial call of motif accoriding to posterior samples, */
	/* motifs will be aligned according to current matrix */
	for(ni=0; ni<nSeqCount; ni++)
	{
		/* init */
		pBase = vSeqMtf[ni]->vSeq[nIndex]->pMatElement;

		/* destroy old lists */
		FLEXMTFSDESTROYLIST(&(vSeqMtf[ni]->pMotifList));
		vSeqMtf[ni]->pMotifList = NULL;
		pPrevSite = NULL;

		/* get posterior probability */
		nLen = vSeqMtf[ni]->nSeqLen;
		for(nj=1; nj<=nMotifNum; nj++)
		{
			vPostSamp = vSeqMtf[ni]->vMonitor[nj]->pMatElement;
			for(nk=0; nk<nLen; nk++)
			{
				vPostSamp[nk] /= (double)nRecordNum;
			}
		}

		/* create new sites according to posterior probability */		
		for(nk=0; nk<FLEXMODULE_CUMPROB_WID; nk++)
		{
			for(nj=1; nj<=nMotifNum; nj++)
			{
				pCumProb->pMatElement[nj-1] += vSeqMtf[ni]->vMonitor[nj]->pMatElement[nk];
			}
		}

		for(nk=0; nk<nLen; nk++)
		{
			/* compute probability */
			nx = nk-FLEXMODULE_CUMPROB_WID-1;
			if(nx >= 0)
			{
				for(nj=1; nj<=nMotifNum; nj++)
				{
					pCumProb->pMatElement[nj-1] -= vSeqMtf[ni]->vMonitor[nj]->pMatElement[nx];
				}
			}

			nx = nk+FLEXMODULE_CUMPROB_WID;
			if(nx < nLen)
			{
				for(nj=1; nj<=nMotifNum; nj++)
				{
					pCumProb->pMatElement[nj-1] += vSeqMtf[ni]->vMonitor[nj]->pMatElement[nx];
				}
			}

			/* create new site */
			for(nj=0; nj<nMotifNum; nj++)
			{
				if((pCumProb->pMatElement[nj] > dLowCut) && (nk+pMeanMotifLen->pMatElement[nj] <= nLen)
					&& (nk+pMotifLen->pMatElement[nj] <= nLen) && 
					(nk+pMotifLen->pMatElement[nj]-pMeanMotifLen->pMatElement[nj] >= 0) )
				{
					pMtfSite = NULL;
					pMtfSite = FLEXMTFSCREATE();
					if(pMtfSite == NULL)
					{
						printf("Error: FlexModule_WritePostMotifToFile, cannot create space for new site!\n");
						exit(EXIT_FAILURE);
					}

					pMtfSite->dProb = pCumProb->pMatElement[nj];
					pMtfSite->nMotifType = nj+1;
					pMtfSite->nSeqId = vSeqMtf[ni]->nId;
					pMtfSite->nStartPos = nk;
					pMtfSite->pNext = NULL;
					pMtfSite->pPrev = NULL;

					dLike1 = FlexModule_ComputeSiteScore(pBase+nk, vMotif[nj]->pPWM, pLogBG0, '+');
					dLike2 = FlexModule_ComputeSiteScore(pBase+nk, vMotif[nj]->pPWM, pLogBG0, '-');
					if(dLike1 > dLike2)
					{
						pMtfSite->nStrand = 0;
						pMtfSite->dScore = dLike1;
						pMtfSite->nStartPos = nk;
					}
					else
					{
						pMtfSite->nStrand = 1;
						pMtfSite->dScore = dLike2;
						pMtfSite->nStartPos = nk+vMotif[nj]->pPWM->nHeight-pMeanMotifLen->pMatElement[nj];
					}

					/* add to site list */
					if(vSeqMtf[ni]->pMotifList == NULL)
					{
						vSeqMtf[ni]->pMotifList = pMtfSite;
						pPrevSite = pMtfSite;
					}
					else
					{
						/* check if overlap. */ 
						nMId = pPrevSite->nMotifType-1;
						
						/* if not overlap, add to motif site list */
						if(pPrevSite->nStartPos+pMeanMotifLen->pMatElement[nMId] <= pMtfSite->nStartPos)
						{
							pPrevSite->pNext = pMtfSite;
							pMtfSite->pPrev = pPrevSite;
							pPrevSite = pMtfSite;
						}

						/* If overlap, choose the best one */
						else
						{
							/* if different motif type */
							if(pPrevSite->nMotifType != pMtfSite->nMotifType)
							{
								if(pMtfSite->dProb <= pPrevSite->dProb)
								{
									FLEXMTFSDESTROY(pMtfSite);
								}
								else
								{
									pTempSite = pPrevSite;
									pPrevSite = pPrevSite->pPrev;
									FLEXMTFSDESTROY(pTempSite);
									if(pPrevSite == NULL)
									{
										vSeqMtf[ni]->pMotifList = pMtfSite;
										pPrevSite = pMtfSite;
									}
									else
									{
										pPrevSite->pNext = pMtfSite;
										pMtfSite->pPrev = pPrevSite;
										pPrevSite = pMtfSite;
									}
								}
							}
							/* if the same motif type */
							else
							{
								if(pMtfSite->dScore <= pPrevSite->dScore)
								{
									FLEXMTFSDESTROY(pMtfSite);
								}
								else
								{
									pTempSite = pPrevSite;
									pPrevSite = pPrevSite->pPrev;
									FLEXMTFSDESTROY(pTempSite);
									if(pPrevSite == NULL)
									{
										vSeqMtf[ni]->pMotifList = pMtfSite;
										pPrevSite = pMtfSite;
									}
									else
									{
										pPrevSite->pNext = pMtfSite;
										pMtfSite->pPrev = pPrevSite;
										pPrevSite = pMtfSite;
									}
								}
							}
						}
					}
				}
			}
		}

		/* FOR TEST: */
		/* sprintf(strLine, "%s_post%d.txt", strFilePath, ni);
		fpOut = NULL;
		fpOut = fopen(strLine, "w");
		for(nk=0; nk<nLen; nk++)
		{
			for(nj=1; nj<=nMotifNum; nj++)
			{
				fprintf(fpOut, "%f\t", vSeqMtf[ni]->vMonitor[nj]->pMatElement[nk]);
			}
			fprintf(fpOut, "\n");
		}
		fclose(fpOut); */
	}
	

	/* sort motifs according to motif type and posterior probability */
	for(ni=0; ni<nSeqCount; ni++)
	{
		/* init */
		while(vSeqMtf[ni]->pMotifList != NULL)
		{
			pMtfSite = vSeqMtf[ni]->pMotifList;
			pTempSite = pMtfSite->pNext;
			if(pTempSite != NULL)
			{
				pTempSite->pPrev = NULL;
			}
			vSeqMtf[ni]->pMotifList = pTempSite;
			pMtfSite->pNext = NULL;
			nMId = pMtfSite->nMotifType-1;

			FLEXMTFSADDTOLIST_SORTBYPROB(vMotifList+nMId, pMtfSite);
		}
	}

	/* only keep the top N motifs */
	for(ni=0; ni<nMotifNum; ni++)
	{
		/* FOR TEST: */
		/* sprintf(strLine, "%s_topN%d.txt", strFilePath, ni);
		fpOut = NULL;
		fpOut = fopen(strLine, "w"); */
		/* END */

		nj = 0;
		pMtfSite = vMotifList[ni];
		while(pMtfSite != NULL)
		{
			/* FOR TEST */
			/* nMId = pMtfSite->nMotifType-1;
			fprintf(fpOut, "%d\t%d\t%d\t", pMtfSite->nSeqId, 
				pMtfSite->nStartPos, (pMtfSite->nStartPos+pMeanMotifLen->pMatElement[nMId]-1));
			if(pMtfSite->nStrand == 0)
			{
				fprintf(fpOut, "+\t");
			}
			else if(pMtfSite->nStrand == 1)
			{
				fprintf(fpOut, "-\t");
			}
			else
			{
				fprintf(fpOut, "?\t");
			}

			nSeqId = pMtfSite->nSeqId;
			if(pMtfSite->nStrand == 0)
			{
				pBi = vSeqMtf[nSeqId]->vSeq[nIndex]->pMatElement + pMtfSite->nStartPos;
				for(nj=0; nj<pMeanMotifLen->pMatElement[nMId]; nj++)
				{
					switch((int)(*pBi))
					{
						case 0: fprintf(fpOut, "A");
							break;
						case 1: fprintf(fpOut, "C");
							break;
						case 2: fprintf(fpOut, "G");
							break;
						case 3: fprintf(fpOut, "T");
							break;
						case 10: fprintf(fpOut, "a");
							break;
						case 11: fprintf(fpOut, "c");
							break;
						case 12: fprintf(fpOut, "g");
							break;
						case 13: fprintf(fpOut, "t");
							break;
						default: fprintf(fpOut, "N");
					}
					pBi++;
				}
			}
			else if(pMtfSite->nStrand == 1)
			{
				pBi = vSeqMtf[nSeqId]->vSeq[nIndex]->pMatElement + pMtfSite->nStartPos + pMeanMotifLen->pMatElement[nMId] - 1;
				for(nj=0; nj<pMeanMotifLen->pMatElement[nMId]; nj++)
				{
					switch((int)(*pBi))
					{
						case 0: fprintf(fpOut, "T");
							break;
						case 1: fprintf(fpOut, "G");
							break;
						case 2: fprintf(fpOut, "C");
							break;
						case 3: fprintf(fpOut, "A");
							break;
						case 10: fprintf(fpOut, "t");
							break;
						case 11: fprintf(fpOut, "g");
							break;
						case 12: fprintf(fpOut, "c");
							break;
						case 13: fprintf(fpOut, "a");
							break;
						default: fprintf(fpOut, "N");
					}
					pBi--;
				}
			}
			fprintf(fpOut, "\t%f\n", pMtfSite->dProb); */
			/* END TEST */


			nj++;
			pMtfSite = pMtfSite->pNext;
			/* FOR TEST */
			if(nj == (int)(pMotifNSample->pMatElement[ni])+5)
			{
				break;
			}
			/* END TEST */
			/* if(nj == (int)(pMotifNSample->pMatElement[ni]))
			{
				break;
			} */
		}

		if(pMtfSite != NULL)
		{
			pPrevSite = pMtfSite->pPrev;
			if(pPrevSite != NULL)
			{
				pPrevSite->pNext = NULL;
			}
			else
			{
				vMotifList[ni] = NULL;
			}
			pMtfSite->pPrev = NULL;
			FLEXMTFSDESTROYLIST(&pMtfSite);
		}


		/* FOR TEST */
		/* fclose(fpOut); */
		/* END */
	}

	/* shift top N motif sites to get the highest sampling rate */
	for(ni=0; ni<nMotifNum; ni++)
	{
		/* create probability space */
		pShiftProb = NULL;
		pShiftProb = CreateDoubleMatrix(1, (2*FLEXMODULE_CUMPROB_WID+1));
		if(pShiftProb == NULL)
		{
			printf("Error: FlexModule_WritePostMotifToFile, cannot create shift probability!\n");
			exit(EXIT_FAILURE);
		}

		/* compute probability */
		pMtfSite = vMotifList[ni];
		while(pMtfSite != NULL)
		{
			if(pMtfSite->nStrand == 0)
			{
				nSeqId = pMtfSite->nSeqId;
				nPos = pMtfSite->nStartPos;
				pMoniter = vSeqMtf[nSeqId]->vMonitor[ni+1]->pMatElement;
				nLen = vSeqMtf[nSeqId]->nSeqLen;
				for(nj=-FLEXMODULE_CUMPROB_WID; nj<=FLEXMODULE_CUMPROB_WID; nj++)
				{
					nk = nj+nPos;
					if( (nk >= 0) && (nk < nLen) )
					{
						pShiftProb->pMatElement[FLEXMODULE_CUMPROB_WID+nj] += pMoniter[nk];
					}
				}
			}
			else if(pMtfSite->nStrand == 1)
			{
				nSeqId = pMtfSite->nSeqId;
				nPos = pMtfSite->nStartPos;
				pMoniter = vSeqMtf[nSeqId]->vMonitor[ni+1]->pMatElement;
				nLen = vSeqMtf[nSeqId]->nSeqLen;
				for(nj=-FLEXMODULE_CUMPROB_WID; nj<=FLEXMODULE_CUMPROB_WID; nj++)
				{
					nk = nj+nPos;
					if( (nk >= 0) && (nk < nLen) )
					{
						pShiftProb->pMatElement[FLEXMODULE_CUMPROB_WID-nj] += pMoniter[nk];
					}
				}
			}
			else
			{
				printf("Error: FlexModule_WritePostMotifToFile, motif site strand wrong!\n");
				exit(EXIT_FAILURE);
			}

			pMtfSite = pMtfSite->pNext;
		}


		/* get shift phase */
		dMax = 0.0;
		nMax = 0;
		for(nj=-FLEXMODULE_CUMPROB_WID; nj<=FLEXMODULE_CUMPROB_WID; nj++)
		{
			if(pShiftProb->pMatElement[FLEXMODULE_CUMPROB_WID+nj] > dMax)
			{
				dMax = pShiftProb->pMatElement[FLEXMODULE_CUMPROB_WID+nj];
				nMax = nj;
			}
		}

		/* shift & update motif PWM */
		DestroyDoubleMatrix(vMotif[ni]->pSampleCount);
		vMotif[ni]->pSampleCount = CreateDoubleMatrix(pMeanMotifLen->pMatElement[ni], nBaseTypeNum);
		if(vMotif[ni]->pSampleCount == NULL)
		{
			printf("Error: FlexModule_WritePostMotifToFile, cannot create sample count!\n");
			exit(EXIT_FAILURE);
		}

		pMtfSite = vMotifList[ni];
		while(pMtfSite != NULL)
		{
			nLen = vSeqMtf[nSeqId]->nSeqLen;
			nSeqId = pMtfSite->nSeqId;

			/* move */

			if(pMtfSite->nStrand == 0)
			{
				pMtfSite->nStartPos += nMax;	
			}
			else if(pMtfSite->nStrand == 1)
			{
				pMtfSite->nStartPos -= nMax;			
			}
			else
			{
				printf("Error: FlexModule_WritePostMotifToFile, motif site strand wrong!\n");
				exit(EXIT_FAILURE);
			}

			/* if out of range, then remove the site */
			if( (pMtfSite->nStartPos < 0) || (pMtfSite->nStartPos+pMeanMotifLen->pMatElement[ni] > nLen) )
			{
				pPrevSite = pMtfSite->pPrev;
				pTempSite = pMtfSite->pNext;

				if( (pPrevSite != NULL) && (pTempSite != NULL) )
				{
					pPrevSite->pNext = pTempSite;
					pTempSite->pPrev = pPrevSite;
				}
				else if(pPrevSite != NULL)
				{
					pPrevSite->pNext = pTempSite;
				}
				else if(pTempSite != NULL)
				{
					pTempSite->pPrev = pPrevSite;
					vMotifList[ni] = pTempSite;
				}

				pMtfSite->pNext = NULL;
				pMtfSite->pPrev = NULL;
				FLEXMTFSDESTROY(pMtfSite);
				pMtfSite = pTempSite;
				continue;
			}

			pBi = vSeqMtf[nSeqId]->vSeq[nIndex]->pMatElement+pMtfSite->nStartPos;
			if(pMtfSite->nStrand == 0)
			{
				FlexModule_AddMotifMatrixCount(vMotif[ni]->pSampleCount, pBi, pMeanMotifLen->pMatElement[ni], '+');
			}
			else
			{
				FlexModule_AddMotifMatrixCount(vMotif[ni]->pSampleCount, pBi, pMeanMotifLen->pMatElement[ni], '-');
			}

			pMtfSite = pMtfSite->pNext;
		}

		/* shift prior count */
		FLEXMTFMSHIFTPRIORCOUNT(vMotif[ni], pMeanMotifLen->pMatElement[ni], nMax, dDefaultPrior);

		/* refresh PWM */
		DestroyDoubleMatrix(vMotif[ni]->pPWM);
		vMotif[ni]->pPWM = CreateDoubleMatrix(pMeanMotifLen->pMatElement[ni], nBaseTypeNum);
		if(vMotif[ni]->pPWM == NULL)
		{
			printf("Error: FlexModule_WritePostMotifToFile, cannot create PWM!\n");
			exit(EXIT_FAILURE);
		}
		FLEXMTFMREFRESH(vMotif[ni]);

		/* release memory */
		DestroyDoubleMatrix(pShiftProb);
	}

	/* resort sites by sequence id and position */
	for(ni=0; ni<nMotifNum; ni++)
	{
		pTempSite = NULL;
		while(vMotifList[ni] != NULL)
		{
			pMtfSite = vMotifList[ni];
			vMotifList[ni] = pMtfSite->pNext;
			pMtfSite->pNext = NULL;
			pMtfSite->pPrev = NULL;
			FLEXMTFSADDTOLIST_SORTBYPOS(&pTempSite, pMtfSite);
		}

		vMotifList[ni] = pTempSite;
	}

	/* export */
	fpOut = NULL;
	fpOut = fopen(strFilePath, "w");
	if(fpOut == NULL)
	{
		printf("Error: cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	/* write motif one by one */
	for(ni=0; ni<nMotifNum; ni++)
	{
		/* get motif score */
		dScore = FLEXMTFMSCORE(vMotif[ni], pBG0);

		/* write head */
		fprintf(fpOut, "****** Motif%d ******\n", ni);
		fprintf(fpOut, "Motif Score: %f\n", dScore);
		fprintf(fpOut, "Motif Matrix: \n");
		pElement1 = vMotif[ni]->pPriorCount->pMatElement;
		pElement2 = vMotif[ni]->pSampleCount->pMatElement;
		for(nj=0; nj<vMotif[ni]->pPriorCount->nHeight; nj++)
		{
			dMax = 0.0;
			nMax = nBaseTypeNum;
			for(nk=0; nk<vMotif[ni]->pPriorCount->nWidth; nk++)
			{
				dTemp = (*pElement1)+(*pElement2);
				fprintf(fpOut, "% 9.7e ", dTemp);
				if(dTemp > dMax)
				{
					dMax = dTemp;
					nMax = nk;
				}
				pElement1++;
				pElement2++;
			}

			fprintf(fpOut, "\n");
			
			switch(nMax)
			{
				case 0: strConsensus[nj] = 'A';
					break;
				case 1: strConsensus[nj] = 'C';
					break;
				case 2: strConsensus[nj] = 'G';
					break;
				case 3: strConsensus[nj] = 'T';
					break;
				default: strConsensus[nj] = 'N';
			}
		}
		strConsensus[nj] = '\0';

		fprintf(fpOut, "\n");
		fprintf(fpOut, "Consensus:\n");
		fprintf(fpOut, "%s\n", strConsensus);

		fprintf(fpOut, "\nMotif Sites:\n");
		pSite = vMotifList[ni];
		while(pSite != NULL)
		{
			/* head information */
			nSeqId = pSite->nSeqId;
			nMId = pSite->nMotifType-1;
			/* fprintf(fpOut, ">%d:%d-%d", pSite->nSeqId, 
				pSite->nStartPos, (pSite->nStartPos+pMeanMotifLen->pMatElement[nMId]-1));
			if(pSite->nStrand == 0)
			{
				fprintf(fpOut, "(+)\n");
			}
			else if(pSite->nStrand == 1)
			{
				fprintf(fpOut, "(-)\n");
			}
			else
			{
				fprintf(fpOut, "(?)\n");
			}
			*/

			fprintf(fpOut, "%d\t%d\t%d\t", pSite->nSeqId, 
				pSite->nStartPos, (pSite->nStartPos+pMeanMotifLen->pMatElement[nMId]-1));
			if(pSite->nStrand == 0)
			{
				fprintf(fpOut, "+\t");
			}
			else if(pSite->nStrand == 1)
			{
				fprintf(fpOut, "-\t");
			}
			else
			{
				fprintf(fpOut, "?\t");
			}

			/* write seq */
			if(pSite->nStrand == 0)
			{
				pBi = vSeqMtf[nSeqId]->vSeq[nIndex]->pMatElement + pSite->nStartPos;
				for(nj=0; nj<pMeanMotifLen->pMatElement[nMId]; nj++)
				{
					switch((int)(*pBi))
					{
						case 0: fprintf(fpOut, "A");
							break;
						case 1: fprintf(fpOut, "C");
							break;
						case 2: fprintf(fpOut, "G");
							break;
						case 3: fprintf(fpOut, "T");
							break;
						case 10: fprintf(fpOut, "a");
							break;
						case 11: fprintf(fpOut, "c");
							break;
						case 12: fprintf(fpOut, "g");
							break;
						case 13: fprintf(fpOut, "t");
							break;
						default: fprintf(fpOut, "N");
					}
					pBi++;
				}
			}
			else if(pSite->nStrand == 1)
			{
				pBi = vSeqMtf[nSeqId]->vSeq[nIndex]->pMatElement + pSite->nStartPos + pMeanMotifLen->pMatElement[nMId] - 1;
				for(nj=0; nj<pMeanMotifLen->pMatElement[nMId]; nj++)
				{
					switch((int)(*pBi))
					{
						case 0: fprintf(fpOut, "T");
							break;
						case 1: fprintf(fpOut, "G");
							break;
						case 2: fprintf(fpOut, "C");
							break;
						case 3: fprintf(fpOut, "A");
							break;
						case 10: fprintf(fpOut, "t");
							break;
						case 11: fprintf(fpOut, "g");
							break;
						case 12: fprintf(fpOut, "c");
							break;
						case 13: fprintf(fpOut, "a");
							break;
						default: fprintf(fpOut, "N");
					}
					pBi--;
				}
			}
			fprintf(fpOut, "\t%f\n", pSite->dProb);

			/* get next */
			pSite = pSite->pNext;
		}

		fprintf(fpOut, "\n\n");
	}

	fclose(fpOut);

	/* release memory */
	DestroyDoubleMatrix(pCumProb);
	for(ni=0; ni<nMotifNum; ni++)
	{
		FLEXMTFSDESTROYLIST(vMotifList+ni);
	}
	free(vMotifList);


	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  FlexModule_BGFit_Main: flexmodule fit background likelihood function   */
/* ----------------------------------------------------------------------- */ 
int FlexModule_BGFit_Main(char strParamFile[])
{
	/* define */

	/* for loading parameters */
	FILE *fpIn;
	char strLine[MED_LINE_LENGTH];
	char *chSep;
	int nlen;
	int nError;
	int ni;

	/* parameters */
	char strWorkPath[MED_LINE_LENGTH];
	char strSeqFile[LINE_LENGTH];
	char strOutFile[LINE_LENGTH];
	int nMotifNum = 0;
	struct tagString **vPWMFile = NULL;
	struct DOUBLEMATRIX *pFreqPrior = NULL;
	struct DOUBLEMATRIX *pControlRate = NULL;

	int nUsePostEnrich = 1;
	int nExportBGFit = 1;
	int nMCNum = 500;
	int nUseCS = 0;
	int nBGOrder = 3;
	char strCSPrefix[LINE_LENGTH];
	char strCSLike[LINE_LENGTH];
	
	/* other variables */
	char strTemp[LINE_LENGTH];
	double dC1,dC2,dC3;

	

	/* --------------- */
	/* load parameters */
	/* --------------- */
	fpIn = NULL;
	fpIn = fopen(strParamFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: FlexModule_BGFit_Main, cannot open the parameter file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strstr(strLine, "[Working Directory]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n%s\n", strLine, chSep);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strWorkPath, chSep);
			
			AdjustDirectoryPath(strWorkPath);
		}


		else if(strstr(strLine, "[FASTA Sequence]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n%s\n", strLine, chSep);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep == '\0')
			{
				printf("Error: FlexModule_BGFit_Main, there are no sequences!\n");
				nError = 1;
				break;
			}
			else
			{
				strcpy(strSeqFile, chSep);
			}
		}

		else if(strstr(strLine, "[Output File]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n%s\n", strLine, chSep);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep == '\0')
			{
				printf("Error: FlexModule_BGFit_Main, you have to specify a file for exporting the results!\n");
				nError = 1;
				break;
			}
			else
			{
				strcpy(strOutFile, chSep);
			}
		}

		else if(strstr(strLine, "[Export Fitted Background]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n%s\n", strLine, chSep);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			nExportBGFit = atoi(chSep);
		}
		
		else if(strstr(strLine, "[Order of Markov Background]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n%s\n", strLine, chSep);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nBGOrder = atoi(chSep);
			}
		}

		else if(strstr(strLine, "[Number of Other Background Components]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n%s\n", strLine, chSep);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nMotifNum = atoi(chSep);
			}
		}

		
		else if(strstr(strLine, "[Background Components]") == strLine)
		{
			/* prepare background count matrix */
			pFreqPrior = NULL;
			pFreqPrior = CreateDoubleMatrix((nMotifNum+1), 2);
			if(pFreqPrior == NULL)
			{
				printf("Error: FlexModule_BGFit_Main, cannot allocate memory for storing prior counts!\n");
				exit(EXIT_FAILURE);
			}

			pControlRate = NULL;
			pControlRate = CreateDoubleMatrix(nMotifNum, 1);
			if(pControlRate == NULL)
			{
				printf("Error: FlexModule_BGFit_Main, cannot allocate memory for storing motif target concentration!\n");
				exit(EXIT_FAILURE);
			}

			vPWMFile = NULL;
			vPWMFile = (struct tagString **)calloc(nMotifNum, sizeof(struct tagString *));
			if(vPWMFile == NULL)
			{
				printf("Error: FlexModule_BGFit_Main, cannot allocate memory for storing PWM path!\n");
				exit(EXIT_FAILURE);
			}

			/* load markov background */
			while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
			{
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] != '\0')
					break;
			}
			sscanf(strLine, "%s %lf %lf", strTemp, &dC1, &dC2);
			DMSETAT(pFreqPrior, 0, 0, dC1);
			DMSETAT(pFreqPrior, 0, 1, dC2);

			/* load motif background */
			for(ni=0; ni<nMotifNum; ni++)
			{
				while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
				{
					StrTrimLeft(strLine);
					StrTrimRight(strLine);
					if(strLine[0] != '\0')
						break;
				}

				sscanf(strLine, "%s %lf %lf %lf", strTemp, &dC1, &dC2, &dC3);
				StringAddTail(vPWMFile+ni, strTemp);
				DMSETAT(pFreqPrior, (ni+1), 0, dC1);
				DMSETAT(pFreqPrior, (ni+1), 1, dC2);
				DMSETAT(pControlRate, ni, 0, dC3);
			}

		}

		else if(strstr(strLine, "[Enrichment Self-determined]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n%s\n", strLine, chSep);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nUsePostEnrich = atoi(chSep);
			}
		}
		
		else if(strstr(strLine, "[MCMC Iteration]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n%s\n", strLine, chSep);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nMCNum = atoi(chSep);
			}
		}

		else if(strstr(strLine, "[Use *.cs?]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n%s\n", strLine, chSep);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nUseCS = atoi(chSep);
			}
		}

		else if(strstr(strLine, "[CS Prefix]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n%s\n", strLine, chSep);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strCSPrefix, chSep);
		}
		
		else if(strstr(strLine, "[CS Likelihood f]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n%s\n", strLine, chSep);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep == '\0')
			{
				strcpy(strCSLike, "NULL");
			}
			else
			{
				strcpy(strCSLike, chSep);
			}
		}

		else if(strLine[0] == '#')
		{
		}

		else
		{
			printf("Error: FlexModule_BGFit_Main, unknown parameters!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	fclose(fpIn);

	if(nError > 0)
	{
		printf("Error: FlexModule_BGFit_Main, cannot load the parameter file correctly!\n");
		exit(EXIT_FAILURE);
	}


	/* flexmodule_bgfit */
	FlexModule_BGFit(strWorkPath, strSeqFile, strOutFile, nExportBGFit,
		/* order of background markov chain, motif number K, motif PWM (count matrix) */
		nBGOrder, nMotifNum, vPWMFile,
		/* prior count for background components, control motif concentration */
		pFreqPrior, nUsePostEnrich, pControlRate,
		/* MC draw number, using conservation indicator, conservation file prefix, conservation likelihood file */
		nMCNum, nUseCS, strCSPrefix, strCSLike);

	/* release memory */
	for(ni=0; ni<nMotifNum; ni++)
	{
		DeleteString(vPWMFile[ni]);
		vPWMFile[ni] = NULL;
	}
	free(vPWMFile);
	DestroyDoubleMatrix(pFreqPrior);
	DestroyDoubleMatrix(pControlRate);

	/* return */
	return PROC_SUCCESS;
}


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
		int nMCNum, int nUseCS, char strCSPrefix[], char strCSLike[])
{
	/* define */
	/* define */
	char strFilePath[MED_LINE_LENGTH];
	
	/* sequences and conservation */
	/* all sequences */
	struct FLEXSEQMOTIF **vSeqMtf = NULL;
	/* sequence number */
	int nSeqCount = 0;
	
	/* background */
	struct DOUBLEMATRIX *pBG = NULL;
	/* log background */
	struct DOUBLEMATRIX *pLogBG = NULL;
	/* background0 */
	struct DOUBLEMATRIX *pBG0 = NULL;
	/* log background0 */
	struct DOUBLEMATRIX *pLogBG0 = NULL;

	/* conservation likelihood */
	struct DOUBLEMATRIX *pCSLike = NULL;
	struct DOUBLEMATRIX **vCSLike = NULL;
	int nCSLikeNum = 0;

	/* motif matrix */
	struct FLEXMOTIFMATRIX **vMotif = NULL;
	struct INTMATRIX *pMotifLen = NULL;
	int nMaxMotifLen = 0;

	/* motif number sample */
	int nEffectLen = 0;
	struct DOUBLEMATRIX *pMotifSiteNum = NULL;
	struct DOUBLEMATRIX *pMotifSiteCount = NULL;
	double dCount0,dCount1;
	double dCountT,dCountP;
	struct DOUBLEMATRIX *pTDR = NULL;
	FILE *fpOut;

	/* number of base types */
	int nBaseTypeNum = 4;
	/* background base transition */
	int nScale = 0;

	/* frequency count */
	struct DOUBLEMATRIX *pFreqCount = NULL;
	double dDefaultPriorCount = 0.5;

	/* others */
	int ni,nj,niter;
	double *pElement;
	/*double *pScoreEle;*/
	/*FILE *fpSeg;*/
	int nRecordNum,nIsRecording;

	/* #################################### */
	/* initialize                           */
	/* #################################### */
	nBaseTypeNum = 4;
	nScale = (int)(pow((double)nBaseTypeNum, (double)(nBGOrder-1)));
	nRecordNum = (int)(nMCNum/2);
	
	/* #################################### */
	/* load initial parameters              */
	/* #################################### */

	/* prior frequency */
	if( (pFreqPrior->nHeight != (nMotifNum+1)) || (pFreqPrior->nWidth != 2) )
	{
		printf("Error: FlexModule_BGFit, prior count matrix dimension does not match motif number!\n");
		exit(EXIT_FAILURE);
	}

	pFreqCount = NULL;
	pFreqCount = DMCLONE(pFreqPrior);
	if(pFreqCount == NULL)
	{
		printf("Error: FlexModule_BGFit, cannot initialize motif site counting matrix!\n");
		exit(EXIT_FAILURE);
	}

	/* prior conservation likelihood */
	pCSLike = NULL;
	vCSLike = NULL;
	if(nUseCS == 1)
	{
		sprintf(strFilePath, "%s%s", strWorkPath, strCSLike);
		pCSLike = DMLOAD(strFilePath);
		if(pCSLike == NULL)
		{
			printf("Error: FlexModule_BGFit, cannot load conservation likelihood matrix!\n");
			exit(EXIT_FAILURE);
		}

		vCSLike = FlexModule_InitScoreLike(pCSLike);
		if(vCSLike == NULL)
		{
			printf("Error: FlexModule_BGFit, cannot create conservation score likelihood function!\n");
			exit(EXIT_FAILURE);
		}

		nCSLikeNum = pCSLike->nHeight-1;
	}

	/* space for counting true motif site number */
	pMotifSiteCount = NULL;
	pMotifSiteCount = CreateDoubleMatrix((nMotifNum+1), 2);
	if(pMotifSiteCount == NULL)
	{
		printf("Error: FlexModule_BGFit, cannot create space for storing sampled motif site count!\n");
		exit(EXIT_FAILURE);
	}

	pMotifSiteNum = NULL;
	pMotifSiteNum = CreateDoubleMatrix(nMotifNum, 1);
	if(pMotifSiteNum == NULL)
	{
		printf("Error: FlexModule_BGFit, cannot create space for storing sampled motif site number!\n");
		exit(EXIT_FAILURE);
	}

	pTDR = NULL;
	pTDR = CreateDoubleMatrix(nMotifNum, 1);
	if(pTDR == NULL)
	{
		printf("Error: FlexModule_BGFit, cannot create space for storing true positive rate!\n");
		exit(EXIT_FAILURE);
	}

	/* prior motif PWM */
	pMotifLen = NULL;
	pMotifLen = CreateIntMatrix(nMotifNum, 1);
	if(pMotifLen == NULL)
	{
		printf("Error: FlexModule_BGFit, cannot create motif length matrix!\n");
		exit(EXIT_FAILURE);
	}

	vMotif = NULL;
	vMotif = (struct FLEXMOTIFMATRIX **)calloc(nMotifNum, sizeof(struct FLEXMOTIFMATRIX *));
	if(vMotif == NULL)
	{
		printf("Error: FlexModule_BGFit, cannot create motif PWM matrix!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nMotifNum; ni++)
	{
		/* create motif */
		vMotif[ni] = NULL;
		vMotif[ni] = FLEXMTFMCREATE();
		if( vMotif[ni] == NULL)
		{
			printf("Error: FlexModule_BGFit, cannot load init motif matrices!\n");
			exit(EXIT_FAILURE);
		}

		/* motif id */
		vMotif[ni]->nMotifType = (ni+1);

		/* init prior matrix */
		sprintf(strFilePath, "%s%s", strWorkPath, vPWMFile[ni]->m_pString);
		vMotif[ni]->pPriorCount = NULL;
		vMotif[ni]->pPriorCount = DMLOAD(strFilePath);
		if( vMotif[ni]->pPriorCount == NULL)
		{
			printf("Error: FlexModule_BGFit, cannot load init motif pseudocount!\n");
			exit(EXIT_FAILURE);
		}
		if(	vMotif[ni]->pPriorCount->nWidth != nBaseTypeNum )
		{
			printf("Error: FlexModule_BGFit, the width of init motif matrix shoule be %d!\n", nBaseTypeNum);
			exit(EXIT_FAILURE);
		}
		pMotifLen->pMatElement[ni] = vMotif[ni]->pPriorCount->nHeight;

		/* init sample count matrix */
		vMotif[ni]->pSampleCount = NULL;
		vMotif[ni]->pSampleCount = CreateDoubleMatrix(pMotifLen->pMatElement[ni], nBaseTypeNum);
		if(vMotif[ni]->pSampleCount == NULL)
		{
			printf("Error: FlexModule_BGFit, cannot create sample count matrix!\n");
			exit(EXIT_FAILURE);
		}

		/* init PWM matrix */
		vMotif[ni]->pPWM = NULL;
		vMotif[ni]->pPWM = CreateDoubleMatrix(pMotifLen->pMatElement[ni], nBaseTypeNum);
		if(vMotif[ni]->pPWM == NULL)
		{
			printf("Error: FlexModule_BGFit, cannot create PWM matrix!\n");
			exit(EXIT_FAILURE);
		}
		FLEXMTFMREFRESH(*(vMotif+ni));

		/* update the maximal motif length */
		if(pMotifLen->pMatElement[ni] > nMaxMotifLen)
			nMaxMotifLen = pMotifLen->pMatElement[ni];
	}


	/* #################################### */
	/* load sequences                       */
	/* #################################### */
	vSeqMtf = NULL;
	vSeqMtf = FLEXSEQMTFLOADSEQ_FOR_BGFIT(&nSeqCount, strWorkPath, strSeqFile, nUseCS, strCSPrefix, nMotifNum);

	/* #################################### */
	/* compute background markov chain      */
	/* #################################### */
	pBG = NULL;
	pBG = FLEXSEQMTFESTIMATENUCLEICBGMC(vSeqMtf, nSeqCount, 0, nBGOrder, nBaseTypeNum);
	if(pBG == NULL)
	{
		printf("Error: FlexModule_BGFit, failure to estimate background markov matrix!\n");
		exit(EXIT_FAILURE);
	}
	pLogBG = NULL;
	pLogBG = DMCLONE(pBG);
	pElement = pLogBG->pMatElement;
	for(ni=0; ni<pLogBG->nHeight; ni++)
	{
		for(nj=0; nj<pLogBG->nWidth; nj++)
		{
			*pElement = log(*pElement);
			pElement++;
		}
	}
	/* DMSAVE(pBG, "bgestimate.txt"); */
	
	pBG0 = NULL;
	pBG0 = FLEXSEQMTFESTIMATENUCLEICBGMC(vSeqMtf, nSeqCount, 0, 0, nBaseTypeNum);
	if(pBG0 == NULL)
	{
		printf("Error: FlexModule_BGFit, failure to estimate background markov matrix!\n");
		exit(EXIT_FAILURE);
	}
	pLogBG0 = NULL;
	pLogBG0 = DMCLONE(pBG0);
	pElement = pLogBG0->pMatElement;
	for(ni=0; ni<pLogBG0->nHeight; ni++)
	{
		for(nj=0; nj<pLogBG0->nWidth; nj++)
		{
			*pElement = log(*pElement);
			pElement++;
		}
	}
	/* DMSAVE(pBG0, "bgestimate0.txt"); */
	
	/* #################################### */
	/* initialize the sampler               */
	/* #################################### */
	/* init background loglikelihood */
	for(ni=0; ni<nSeqCount; ni++)
	{
		FlexModule_InitBGLogLike(vSeqMtf[ni], 0,
			nBGOrder, pLogBG, pLogBG0, 
			nUseCS, nCSLikeNum, vCSLike);

		nEffectLen += FlexModule_BGFit_InitMotifLogLike(vSeqMtf[ni], 0,
			nBGOrder, pLogBG, pLogBG0,
			nMotifNum, vMotif, pMotifLen, nMaxMotifLen, 
			nUseCS, nCSLikeNum, vCSLike);
	}

	/* init sites */
	for(ni=0; ni<nSeqCount; ni++)
	{
 		FlexModule_BGFit_InitSite(vSeqMtf[ni], 0,
			pFreqPrior, pFreqCount,
			nBGOrder, pLogBG, pLogBG0,  
			nMotifNum, vMotif, pMotifLen, nMaxMotifLen,
			nUseCS, nCSLikeNum, vCSLike);
	}

	/* DMSAVE(pFreqCount, "temp.txt");
	DMSAVE(vMotif[0]->pSampleCount, "motif0.txt");
	DMSAVE(vMotif[1]->pSampleCount, "motif1.txt");
	DMSAVE(vMotif[2]->pSampleCount, "motif2.txt"); */

	/* #################################### */
	/* flexmodule sampler                   */
	/* #################################### */
	for(niter=0; niter<nMCNum; niter++)
	{
		if(niter%50 == 0)
		{
			printf("iter=%d...\n", niter);
		}

		if((niter+nRecordNum) < nMCNum)
			nIsRecording = 0;
		else
			nIsRecording = 1;

		/* sequence-wise site sampling */
		for(ni=0; ni<nSeqCount; ni++)
		{
			/* Sample potential site */
			FlexModule_BGFit_SampleBA(vSeqMtf[ni], 0,  
				         nUsePostEnrich, pFreqPrior, pFreqCount,
						 nBGOrder, pLogBG, pLogBG0,
						 nMotifNum, vMotif, pMotifLen, nMaxMotifLen, 
						 nUseCS, nCSLikeNum, vCSLike, 
						 nIsRecording);
		}

		/* get true positive rate and recording */
		if(nIsRecording > 0)
		{
			dCount0 = DMGETAT(pFreqCount, 0, 0)-DMGETAT(pFreqPrior, 0, 0);
			dCount1 = DMGETAT(pFreqCount, 0, 1)-DMGETAT(pFreqPrior, 0, 1);
			
			dCount0 += DMGETAT(pMotifSiteCount, 0, 0);
			DMSETAT(pMotifSiteCount, 0, 0, dCount0);
			dCount1 += DMGETAT(pMotifSiteCount, 0, 1);
			DMSETAT(pMotifSiteCount, 0, 1, dCount1);
		
			for(nj=0; nj<nMotifNum; nj++)
			{
				dCount0 = DMGETAT(pFreqCount, (nj+1), 0)-DMGETAT(pFreqPrior, (nj+1), 0);
				dCount1 = DMGETAT(pFreqCount, (nj+1), 1)-DMGETAT(pFreqPrior, (nj+1), 1);
				
				pMotifSiteNum->pMatElement[nj] = dCount0+dCount1;
				dCount0 += DMGETAT(pMotifSiteCount, (nj+1), 0);
				DMSETAT(pMotifSiteCount, (nj+1), 0, dCount0);
				dCount1 += DMGETAT(pMotifSiteCount, (nj+1), 1);
				DMSETAT(pMotifSiteCount, (nj+1), 1, dCount1);

				if(pControlRate->pMatElement[nj] < 0.0)
				{
					pTDR->pMatElement[nj] = 1.0;
				}
				else
				{
					if( pMotifSiteNum->pMatElement[nj] == 0)
					{
						pTDR->pMatElement[nj] = 1.0;
					}
					else
					{
						pTDR->pMatElement[nj] = (double)nEffectLen*pControlRate->pMatElement[nj]/pMotifSiteNum->pMatElement[nj];
						if(pTDR->pMatElement[nj] > 1.0)
							pTDR->pMatElement[nj] = 1.0;
					}
				}
			}

		
			/* set likelihood */
			for(ni=0; ni<nSeqCount; ni++)
			{
				FlexModule_BGFit_UpdateLikelihood(vSeqMtf[ni], 0,  
				         pTDR, nBGOrder, pLogBG, pLogBG0,
						 nMotifNum, vMotif, pMotifLen, nMaxMotifLen, 
						 nUseCS, nCSLikeNum, vCSLike);
			}
		}
	}

	/* #################################### */
	/* export posterior fitted background   */
	/* #################################### */
	if(nExportBGFit == 1)
	{
		for(ni=0; ni<nSeqCount; ni++)
		{
			sprintf(strFilePath, "%s%s_%d_%s.bgf", strWorkPath, strOutFile,
				ni, vSeqMtf[ni]->strAlias);
			FlexModule_BGFit_ExportLogLikelihood(vSeqMtf[ni], 0, 
				nRecordNum, strFilePath);
		}
	}

	DMPDIVTS(pMotifSiteCount, (double)nRecordNum);

	sprintf(strFilePath, "%s%s.freq", strWorkPath, strOutFile);
	fpOut = NULL;
	fpOut = fopen(strFilePath, "w");
	if(fpOut == NULL)
	{
		printf("Error: FlexModule_BGFit, cannot open file for saving posterior counts!\n");
		exit(EXIT_FAILURE);
	}

	dCount0 = DMGETAT(pMotifSiteCount, 0, 0);
	dCount1 = DMGETAT(pMotifSiteCount, 0, 1);
	dCountT = dCount0+dCount1;
	
	fprintf(fpOut, "%f\t%f\t", dCount0, dCount1);

	dCount0 += DMGETAT(pFreqPrior, 0, 0);
	dCount1 += DMGETAT(pFreqPrior, 0, 1);
	dCountP = dCount0+dCount1;
	DMSETAT(pMotifSiteCount, 0, 0, dCount0);
	DMSETAT(pMotifSiteCount, 0, 1, dCount1);

	fprintf(fpOut, "%f\n", dCountP);

	for(nj=0; nj<nMotifNum; nj++)
	{
		dCount0 = DMGETAT(pMotifSiteCount, (nj+1), 0);
		dCount1 = DMGETAT(pMotifSiteCount, (nj+1), 1);

		fprintf(fpOut, "%f\t%f\t", dCount0, dCount1);
		dCount0 += DMGETAT(pFreqPrior, (nj+1), 0);
		dCount1 += DMGETAT(pFreqPrior, (nj+1), 1);
		DMSETAT(pMotifSiteCount, (nj+1), 0, dCount0);
		DMSETAT(pMotifSiteCount, (nj+1), 1, dCount1);
		pMotifSiteNum->pMatElement[nj] = dCount0+dCount1;

		if(dCountP > 0.0)
		{
			fprintf(fpOut, "%f\t%9.7e\n",  
				pMotifSiteNum->pMatElement[nj], pMotifSiteNum->pMatElement[nj]/dCountP);
		}
		else
		{
			fprintf(fpOut, "%f\t%9.7e\n",  
				pMotifSiteNum->pMatElement[nj], 0.0);
		}
	}

	fclose(fpOut);

	/* #################################### */
	/* export posterior sample              */
	/* #################################### */
	for(nj=0; nj<nMotifNum; nj++)
	{
		pMotifSiteNum->pMatElement[nj] -= (DMGETAT(pFreqPrior, (nj+1), 0)+DMGETAT(pFreqPrior, (nj+1), 1));
	}
	sprintf(strFilePath, "%s%s.map", strWorkPath, strOutFile);
	FlexModule_BGFit_WritePostMotifToFile(strFilePath, nMotifNum, 
					vMotif, pBG0, pLogBG0, 
					pMotifLen, pMotifSiteCount, 
					nUsePostEnrich, pFreqPrior, 
					pMotifSiteNum, nSeqCount, 
					vSeqMtf, 
					0, nRecordNum);

	/* #################################### */
	/* release memory                       */
	/* #################################### */
	DestroyDoubleMatrix(pFreqCount);
	DestroyIntMatrix(pMotifLen);
	
	if(nUseCS == 1)
	{
		for(ni=0; ni<(pCSLike->nHeight-1); ni++)
		{
			DestroyDoubleMatrix(vCSLike[ni]);
			vCSLike[ni] = NULL;
		}
		free(vCSLike);
		DestroyDoubleMatrix(pCSLike);
	}
	
	for(ni=0; ni<nMotifNum; ni++)
	{
		FLEXMTFMDESTROY(vMotif[ni]);
		vMotif[ni] = NULL;
	}
	free(vMotif);

	DestroyDoubleMatrix(pMotifSiteCount);
	DestroyDoubleMatrix(pMotifSiteNum);
	DestroyDoubleMatrix(pTDR);

	for(ni=0; ni<nSeqCount; ni++)
	{
		FLEXSEQMTFDESTROY(vSeqMtf[ni]);
		vSeqMtf[ni] = NULL;
	}
	free(vSeqMtf);

	DestroyDoubleMatrix(pBG);
	DestroyDoubleMatrix(pLogBG);
	DestroyDoubleMatrix(pBG0);
	DestroyDoubleMatrix(pLogBG0);
	
	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  FlexModule_BGFit_InitMotifLogLike: Initialize motif log likelihood     */
/*  return effective length of base pairs.                                 */
/* ----------------------------------------------------------------------- */ 
int FlexModule_BGFit_InitMotifLogLike(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,
			int nBGOrder, struct DOUBLEMATRIX *pLogBG, struct DOUBLEMATRIX *pLogBG0, 
			int nMotifNum, struct FLEXMOTIFMATRIX **vMotif, struct INTMATRIX *pMotifLen, 
			int nMaxMotifLen,
			int nUseCS, int nCSLikeNum, struct DOUBLEMATRIX **vCSLike)
{
	/* define */
	int ni,nj,nk;
	int nLen,nEffectLen;
	unsigned char *pBase;
	int nBaseTypeNum;
	int nId;
	double dLike;
	int nMask;

	nLen = pSeqMtf->vSeq[nIndex]->nWidth;
	pBase = pSeqMtf->vSeq[nIndex]->pMatElement;
	nBaseTypeNum = pLogBG->nWidth;
	
	/* ######################### */
	/* scan the whole sequence   */
	/* ######################### */
	/* sample base by base */
	nEffectLen = 0;
	for(ni=0; ni<nLen; ni++)
	{
		/* if repeat or 'N', skip it */
		if(pBase[ni] >= nBaseTypeNum)
		{
			continue;
		}

		nEffectLen++;

		/* get motif likelihood */
		for(nj=0; nj<nMotifNum; nj++)
		{
			/* '+' '-' strand */
			for(nk=0; nk<2; nk++)
			{
				/* compute context score & probabiity */
				if(nk == 0)
				{
					dLike =	FlexBaseLogLikelihood_Motif(pSeqMtf, nIndex, ni,
						vMotif[nj]->pPWM, '+', nMaxMotifLen,
						nBGOrder, pLogBG, pLogBG0, 
						nUseCS, nCSLikeNum, vCSLike,
						&nMask);
				}
				else
				{
					dLike =	FlexBaseLogLikelihood_Motif(pSeqMtf, nIndex, ni, 
						vMotif[nj]->pPWM, '-', nMaxMotifLen,
						nBGOrder, pLogBG, pLogBG0, 
						nUseCS, nCSLikeNum, vCSLike,
						&nMask);
				}

				if(nMask == 1)
				{
					nId = 2+2*nj+nk;
					pSeqMtf->vStatus[nId]->pMatElement[ni] = 1;
				}
				else
				{
					nId = 2+2*nj+nk;
					pSeqMtf->vStatus[nId]->pMatElement[ni] = 0;
					nId = 5+2*nj+nk;
					pSeqMtf->vMonitor[nId]->pMatElement[ni] = dLike;
				}
			}
		}
	} 
		
	/* return */
	return nEffectLen;
}

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_BGFit_InitSite: Initialize motif sites at random for        */
/*  background fitting.                                                    */
/* ----------------------------------------------------------------------- */ 
int FlexModule_BGFit_InitSite(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,
			struct DOUBLEMATRIX *pFreqPrior0, struct DOUBLEMATRIX *pFreqCount,
			int nBGOrder, struct DOUBLEMATRIX *pLogBG, struct DOUBLEMATRIX *pLogBG0,  
			int nMotifNum, struct FLEXMOTIFMATRIX **vMotif, 
			struct INTMATRIX *pMotifLen, int nMaxMotifLen,
			int nUseCS, int nCSLikeNum, struct DOUBLEMATRIX **vCSLike)
{
	/* define */
	/* sampler probability */
	struct DOUBLEMATRIX *pFreqPrior = NULL;
	struct DOUBLEMATRIX *pFreqPost = NULL;
	struct BYTEMATRIX *pFreqValid = NULL;
	
	/* sampler status */
	int ni,nj,nx,nLen;
	unsigned char *pAi,*pSi,*pBase;
	struct FLEXMOTIFSITE *pNewSite,*pCurrentSite;
	struct FLEXMOTIFSITE **vSite;

	/* sampler likelihood */
	int nBaseTypeNum;
	int nStateId, nStrand;
	int nMotifId;
	double dTemp;
	double dBaseLike;
	double dLike;
	int nMask;
	int nId1,nId2;

	/* init */
	pFreqPrior = NULL;
	pFreqPrior = DMCLONE(pFreqPrior0); 
	if(pFreqPrior == NULL)
	{
		printf("Error: FlexModule_BGFit_InitSite, cannot create prior probability matrix!\n");
		exit(EXIT_FAILURE);
	}

	pFreqPost = NULL;
	pFreqPost = CreateDoubleMatrix(pFreqCount->nHeight, pFreqCount->nWidth); 
	if(pFreqPost == NULL)
	{
		printf("Error: FlexModule_BGFit_InitSite, cannot create posterior probability matrix!\n");
		exit(EXIT_FAILURE);
	}

	pFreqValid = NULL;
	pFreqValid = CreateByteMatrix(pFreqCount->nHeight, pFreqCount->nWidth); 
	if(pFreqValid == NULL)
	{
		printf("Error: FlexModule_BGFit_InitSite, cannot create sequence status indicator matrix!\n");
		exit(EXIT_FAILURE);
	}

	/* get prior probability */
	FlexModule_NormalizeFreq(pFreqPrior);
	FlexModule_LogFreq(pFreqPrior);
	
	/* init */
	nLen = pSeqMtf->vSeq[nIndex]->nWidth;
	pBase = pSeqMtf->vSeq[nIndex]->pMatElement;
	pAi = pSeqMtf->vStatus[0]->pMatElement;
	pSi = pSeqMtf->vStatus[1]->pMatElement;
	vSite = pSeqMtf->vMotif;
	pCurrentSite = pSeqMtf->pMotifList;
	
	nBaseTypeNum = pLogBG->nWidth;
	
	/* ######################### */
	/* scan the whole sequence   */
	/* ######################### */
	
	/* set background likelihood */
	dBaseLike = DMGETAT(pFreqPrior, 0, 0);
	BMSETAT(pFreqValid, 0, 0, 1);
		
	for(ni=0; ni<nLen; ni++)
	{
		/* if repeat or 'N', skip it */
		if(pBase[ni] >= nBaseTypeNum)
			continue;

		/* set background likelihood */
		DMSETAT(pFreqPost, 0, 0, dBaseLike);
	
		/* get motif likelihood */
		for(nj=0; nj<nMotifNum; nj++)
		{
			nMask = 0;
			for(nx=0; nx<pMotifLen->pMatElement[nj]; nx++)
			{
				if(pSi[ni+nx] > 0)
				{
					nMask = 1;
					break;
				}
			}

			/* '+' strand */
			nId1 = 2+2*nj;
			nId2 = 5+2*nj;
			if( (nMask == 0) && (pSeqMtf->vStatus[nId1]->pMatElement[ni] == 0) )
			{
				dLike = DMGETAT(pFreqPrior,(nj+1),0)+pSeqMtf->vMonitor[nId2]->pMatElement[ni];
				DMSETAT(pFreqPost, nj+1, 0, dLike);
				BMSETAT(pFreqValid, nj+1, 0, 1);
			}
			else
			{
				BMSETAT(pFreqValid, nj+1, 0, 0);
			}

			
			/* '-' strand */
			nId1 = 2+2*nj+1;
			nId2 = 5+2*nj+1;
			if( (nMask == 0) && (pSeqMtf->vStatus[nId1]->pMatElement[ni] == 0) )
			{
				dLike = DMGETAT(pFreqPrior,(nj+1),1)+pSeqMtf->vMonitor[nId2]->pMatElement[ni];
				DMSETAT(pFreqPost, nj+1, 1, dLike);
				BMSETAT(pFreqValid, nj+1, 1, 1);
			}
			else
			{
				BMSETAT(pFreqValid, nj+1, 1, 0);
			}
		}

		/* normalize posterior probability */
		FlexModule_NormalizePost(pFreqPost, pFreqValid);
		
		/* sample motifs */
		FlexModule_SamplePost(pFreqPost, pFreqValid, &nStateId, &nStrand);

		if(nStateId > 0)
		{
			/* generate pNewSite */
			pNewSite = NULL;
			pNewSite = FLEXMTFSCREATE();
			if(pNewSite == NULL)
			{
				printf("Error: FlexModule_InitSite, cannot create new binding site!\n");
				exit(EXIT_FAILURE);
			}
			pNewSite->nMotifType = nStateId;
			pNewSite->nSeqId = pSeqMtf->nId;
			pNewSite->nStartPos = ni;
			pNewSite->nStrand = nStrand;
			
			vSite[ni] = pNewSite;
			pSeqMtf->nSiteNum += 1;
			if(pSeqMtf->pMotifList == NULL)
			{
				pSeqMtf->pMotifList = pNewSite;
				pCurrentSite = pNewSite;
			}
			else
			{
				pCurrentSite->pNext = pNewSite;
				pNewSite->pPrev = pCurrentSite;
				pCurrentSite = pNewSite;
			}

			
			/* update pAi and pSi */
			pAi[ni] = 1;
			nMotifId = nStateId-1;
			for(nx=0; nx < pMotifLen->pMatElement[nMotifId]; nx++)
			{
				pSi[ni+nx] = nStateId;
			}
		}

		/* update count */
		dTemp = DMGETAT(pFreqCount, nStateId, nStrand)+1.0;
		DMSETAT(pFreqCount, nStateId, nStrand, dTemp);
	}
	
	/* release memory */
	DestroyDoubleMatrix(pFreqPrior);
	DestroyDoubleMatrix(pFreqPost);
	DestroyByteMatrix(pFreqValid);

	/* return */
	return PROC_SUCCESS;
}

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
						 int nIsRecording)
{
	/* define */
	/* sampler probability */
	struct DOUBLEMATRIX *pFreqPost = NULL;
	struct BYTEMATRIX *pFreqValid = NULL;

	/* context sites */
	struct FLEXMOTIFSITE *pPrev,*pNext,*pCurrent;

	/* sampler status */
	int ni,nj,nk,nLen,nId1,nId2,nx;
	unsigned char *pAi,*pSi,*pBase;
	struct FLEXMOTIFSITE **vSite;
	double *pPost;

	/* sampler likelihood */
	int nBaseTypeNum;
	double dTemp;
	int nMotifId;
	int nStateId, nStrand;
	
	double dLike;
	int nMask;

	/* init */
	pFreqPost = NULL;
	pFreqPost = CreateDoubleMatrix(pFreqCount->nHeight, pFreqCount->nWidth); 
	if(pFreqPost == NULL)
	{
		printf("Error: FlexModule_BGFit_SampleBA, cannot create posterior probability matrix!\n");
		exit(EXIT_FAILURE);
	}

	pFreqValid = NULL;
	pFreqValid = CreateByteMatrix(pFreqCount->nHeight, pFreqCount->nWidth); 
	if(pFreqValid == NULL)
	{
		printf("Error: FlexModule_BGFit_SampleBA, cannot create sequence status indicator matrix!\n");
		exit(EXIT_FAILURE);
	}

	/* init */
	nLen = pSeqMtf->vSeq[nIndex]->nWidth;
	pBase = pSeqMtf->vSeq[nIndex]->pMatElement;
	pAi = pSeqMtf->vStatus[0]->pMatElement;
	pSi = pSeqMtf->vStatus[1]->pMatElement;
	pPost = pSeqMtf->vMonitor[1]->pMatElement;
	vSite = pSeqMtf->vMotif;
	pCurrent = pSeqMtf->pMotifList;
	
	nBaseTypeNum = pLogBG->nWidth;
	
	/* ######################### */
	/* scan the whole sequence   */
	/* ######################### */
	
	/* initialize pPrev, pNext */
	pPrev = NULL;
	pNext = pSeqMtf->pMotifList;
	
	/* sample base by base */
	BMSETAT(pFreqValid, 0, 0, 1);
	for(ni=0; ni<nLen; ni++)
	{
		/* clear old motif status, and get updated context */
		FlexModule_BGFit_ClearOldStatus(pSeqMtf, nIndex, ni, 
			nBaseTypeNum,
			&pPrev, &pNext, 
			pFreqCount, 
			nMotifNum, vMotif, 
			pMotifLen, nMaxMotifLen);


		/* if repeat or 'N', skip it */
		if(pBase[ni] >= nBaseTypeNum)
		{
			continue;
		}

		/* get motif likelihood and set posterior probability */
		if( nUsePostEnrich == 1)
		{
			dLike = log(DMGETAT(pFreqCount, 0, 0));
		}
		else
		{
			dLike = log(DMGETAT(pFreqPrior, 0, 0));
		}

		DMSETAT(pFreqPost, 0, 0, dLike);
		
		for(nj=0; nj<nMotifNum; nj++)
		{
			nStateId = nj+1;

			nMask = 0;
			for(nx=0; nx<pMotifLen->pMatElement[nj]; nx++)
			{
				if(pSi[ni+nx] > 0)
				{
					nMask = 1;
					break;
				}
			}

			/* '+' '-' strand */
			for(nk=0; nk<2; nk++)
			{
				/* compute context score & probabiity */
				nId1 = 2+2*nj+nk;
				nId2 = 5+2*nj+nk;
				if( nUsePostEnrich == 1)
				{
					dLike = log(DMGETAT(pFreqCount, nStateId, nk))+pSeqMtf->vMonitor[nId2]->pMatElement[ni];
				}
				else
				{
					dLike = log(DMGETAT(pFreqPrior, nStateId, nk))+pSeqMtf->vMonitor[nId2]->pMatElement[ni];
				}

				if( (nMask == 0) && (pSeqMtf->vStatus[nId1]->pMatElement[ni] == 0) )
				{
					DMSETAT(pFreqPost, nStateId, nk, dLike);
					BMSETAT(pFreqValid, nStateId, nk, 1);
				}
				else
				{
					DMSETAT(pFreqPost, nStateId, nk, 0.0);
					BMSETAT(pFreqValid, nStateId, nk, 0);
				}
			}
		}
 
		/* normalize posterior probability */
		FlexModule_NormalizePost(pFreqPost, pFreqValid);
		
		/* sample motifs */
		FlexModule_SamplePost(pFreqPost, pFreqValid, &nStateId, &nStrand);

		/* update motif frequency count */
		dTemp = DMGETAT(pFreqCount, nStateId, nStrand)+1.0;
		DMSETAT(pFreqCount, nStateId, nStrand, dTemp);

		/* update the new context */
		if(nStateId > 0)
		{
			nMotifId = nStateId-1;

			/* sample Ai */
			pAi[ni] = 1;
				
			/* add new Bi site */
			pSeqMtf->nSiteNum += 1;
			FlexModule_BGFit_UpdateNewStatus(pSeqMtf, nIndex, ni,
				nStateId, nStrand, nBaseTypeNum,
				&pPrev, &pNext, 
				nMotifNum, vMotif, 
				pMotifLen, nMaxMotifLen);
		}
		else
		{
		
		}
	}

	/* release memory */
	DestroyDoubleMatrix(pFreqPost);
	DestroyByteMatrix(pFreqValid);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_BGFit_ClearOldStatus: clear old status of a given position. */
/* ----------------------------------------------------------------------- */ 
int FlexModule_BGFit_ClearOldStatus(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, int nPos, 
			int nBaseTypeNum,
			struct FLEXMOTIFSITE **pPrev, struct FLEXMOTIFSITE **pNext, 
			struct DOUBLEMATRIX *pFreqCount, 
			int nMotifNum, struct FLEXMOTIFMATRIX **vMotif, 
			struct INTMATRIX *pMotifLen, int nMaxMotifLen)
{
	/* define */
	int nj,nLen;
	unsigned char *pBase,*pAi,*pSi;
	struct FLEXMOTIFSITE **vSite;
	struct FLEXMOTIFSITE *pSite;
	int nMotifId,nMotifLen,nMotifStrand;
	double dTemp;
	
	/* init */
	nLen = pSeqMtf->vSeq[nIndex]->nWidth;
	pBase = pSeqMtf->vSeq[nIndex]->pMatElement;
	pAi = pSeqMtf->vStatus[0]->pMatElement;
	pSi = pSeqMtf->vStatus[1]->pMatElement;
	vSite = pSeqMtf->vMotif;

	/* check the position in question ! */
	/* if repeat or 'N', skip it */
	if(pBase[nPos] >= nBaseTypeNum)
	{
		/* return */
		return PROC_SUCCESS;
	}

	/* if effective base */

	/* if Bi>0, subtract motif count, motif matrix, change context */
	if(vSite[nPos] != NULL)
	{
		/* get pointer */
		pSite = vSite[nPos];
		nMotifId = pSite->nMotifType;
		nMotifLen = pMotifLen->pMatElement[nMotifId-1];
		nMotifStrand = pSite->nStrand;

		/* subtract frequency */
		dTemp = DMGETAT(pFreqCount, nMotifId, nMotifStrand) - 1.0;
		DMSETAT(pFreqCount, nMotifId, nMotifStrand, dTemp);

		/* check if Ai == 1, subtract from PWM count */
		if(pAi[nPos] == 0)
		{
			printf("Error: FlexModule_BGFit_ClearOldStatus, motif status out of track!\n");
			exit(EXIT_FAILURE);
		}
		
		pAi[nPos] = 0;
		vSite[nPos] = NULL;

		/* clear coverage Si */
		for(nj=0; nj<nMotifLen; nj++)
		{
			if(nPos+nj >= nLen)
			{
				printf("Error: FlexModule_BGFit_ClearOldStatus, base index out of range!\n");
				exit(EXIT_FAILURE);
			}
			pSi[nPos+nj] = 0;
		}

		/* reconnect the motif site list */
		if( (*pNext != pSite) || (*pPrev != pSite->pPrev) )
		{
			printf("Error: FlexModule_BGFit_ClearOldStatus, motif list and context not synchronized!\n");
			exit(EXIT_FAILURE);
		}
		*pNext = pSite->pNext;

		if((*pPrev == NULL) && (*pNext == NULL))
		{
			pSeqMtf->pMotifList = NULL;
		}
		else if(*pPrev == NULL)
		{
			(*pNext)->pPrev = *pPrev;
			pSeqMtf->pMotifList = *pNext;
		}
		else if(*pNext == NULL)
		{
			(*pPrev)->pNext = *pNext;
		}
		else
		{
			(*pPrev)->pNext = *pNext;
			(*pNext)->pPrev = *pPrev;
		}

		/* destroy the old site */
		FLEXMTFSDESTROY(pSite);
		pSeqMtf->nSiteNum -= 1;
	}
	/* if Bi=0, subtract background count */
	else
	{
		if(pAi[nPos] != 0)
		{
			printf("Error: FlexModule_BGFit_ClearOldStatus, motif indicator wrong !\n");
			exit(EXIT_FAILURE);
		}

		/* subtract frequency */
		dTemp = DMGETAT(pFreqCount, 0, 0) - 1.0;
		DMSETAT(pFreqCount, 0, 0, dTemp);

		/* there is no need to update pPrev, pNext, pContext */		
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_BGFit_UpdateNewStatus: update new context status.           */
/* ----------------------------------------------------------------------- */
int FlexModule_BGFit_UpdateNewStatus(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, int nPos,
			int nMotifId, int nStrand,
			int nBaseTypeNum,
			struct FLEXMOTIFSITE **pPrev, struct FLEXMOTIFSITE **pNext, 
			int nMotifNum, struct FLEXMOTIFMATRIX **vMotif, 
			struct INTMATRIX *pMotifLen, int nMaxMotifLen)
{
	/* define */
	int ni,nLen;
	int nMId;
	int nMotifLen;
	unsigned char *vSi;
	struct FLEXMOTIFSITE *pSite;
	struct FLEXMOTIFSITE **vSite;

	/* init */
	nMId = nMotifId-1;
	nMotifLen = pMotifLen->pMatElement[nMId];
	vSi = pSeqMtf->vStatus[1]->pMatElement;
	vSite = pSeqMtf->vMotif;
	nLen = pSeqMtf->vSeq[nIndex]->nWidth;

	/* add new site */
	/* create new site */
	pSite = NULL;
	pSite = FLEXMTFSCREATE();
	if(pSite == NULL)
	{
		printf("Error: FlexModule_BGFit_UpdateNewStatus, cannot create new site!\n");
		exit(EXIT_FAILURE);
	}
	pSite->nMotifType = nMotifId;
	pSite->nSeqId = pSeqMtf->nId;
	pSite->nStartPos = nPos;
	pSite->nStrand = nStrand;
	pSite->pContextPos = NULL;
	pSite->pContextMotif = NULL;
	pSite->pContextStrand = NULL;
	pSite->pNext = NULL;
	pSite->pPrev = NULL;

	/* add the new site to the current list */
	vSite[nPos] = pSite;

	if( (*pPrev == NULL) && (*pNext == NULL) )
	{
		pSeqMtf->pMotifList = pSite;
		*pPrev = pSite;
	}
	else if(*pPrev == NULL)
	{
		pSite->pNext = *pNext;
		(*pNext)->pPrev = pSite;
		pSeqMtf->pMotifList = pSite;
		*pPrev = pSite;
	}
	else if(*pNext == NULL)
	{
		pSite->pPrev = *pPrev;
		(*pPrev)->pNext = pSite;
		*pPrev = pSite;
	}
	else
	{
		(*pPrev)->pNext = pSite;
		pSite->pPrev = (*pPrev);
		(*pNext)->pPrev = pSite;
		pSite->pNext = *pNext;
		*pPrev = pSite;
	}

	/* update vSi */
	for(ni=0; ni<nMotifLen; ni++)
	{
		if( (nPos+ni) >= nLen)
		{
			printf("Error: FlexModule_BGFit_UpdateNewStatus, base index out of range!\n");
			exit(EXIT_FAILURE);
		}
		vSi[nPos+ni] = nMotifId;
	}

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  FlexModule_BGFit_UpdateLikelihood: update likelihood.                  */
/* ----------------------------------------------------------------------- */
int FlexModule_BGFit_UpdateLikelihood(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,  
				         struct DOUBLEMATRIX *pTDR, int nBGOrder, 
						 struct DOUBLEMATRIX *pLogBG, struct DOUBLEMATRIX *pLogBG0,
						 int nMotifNum, struct FLEXMOTIFMATRIX **vMotif, 
						 struct INTMATRIX *pMotifLen, int nMaxMotifLen, 
						 int nUseCS, int nCSLikeNum, struct DOUBLEMATRIX **vCSLike)
{
	/* define */
	int ni,nk,nLen,nPos;
	int nMId;
	int nBaseTypeNum;
	int nMotifLen;
	struct FLEXMOTIFSITE *pSite;
	unsigned char *pBase;
	double *pMLike,*pMCount,*pOccur,*pBL,*pPost;
	double dTDR,dTemp;
	
	/* init */
	nLen = pSeqMtf->vSeq[nIndex]->nWidth;
	pBase = pSeqMtf->vSeq[nIndex]->pMatElement;
	pBL = pSeqMtf->vMonitor[0]->pMatElement;
	pPost = pSeqMtf->vMonitor[1]->pMatElement;
	pMLike = pSeqMtf->vMonitor[2]->pMatElement;
	pMCount = pSeqMtf->vMonitor[3]->pMatElement;
	pOccur = pSeqMtf->vMonitor[4]->pMatElement;
	nBaseTypeNum = pLogBG->nWidth;

	/* update likelihood one by one */
	pSite = pSeqMtf->pMotifList;
	while(pSite != NULL)
	{
		nMId = pSite->nMotifType-1;
		nMotifLen = pMotifLen->pMatElement[nMId];
		nPos = pSite->nStartPos;
		dTDR = pTDR->pMatElement[nMId];
		pOccur[nPos] += dTDR;

		if(pSite->nStrand == 0)
		{
			for(ni=0; ni<nMotifLen; ni++)
			{
				pMCount[nPos+ni] += dTDR;
				nk = (int)pBase[nPos+ni];
				dTemp = DMGETAT(vMotif[nMId]->pPWM, ni, nk);
				/* if(pBL[nPos+ni] > dTemp)
					dTemp = pBL[nPos+ni]; */
				pMLike[nPos+ni] += dTDR*exp(dTemp);
				pPost[nPos+ni] -= 1.0;
			}
		}
		else
		{
			for(ni=0; ni<nMotifLen; ni++)
			{
				pMCount[nPos+ni] += dTDR;
				nk = nBaseTypeNum-1-(int)pBase[nPos+ni];
				dTemp = DMGETAT(vMotif[nMId]->pPWM, nMotifLen-1-ni, nk);
				/* if(pBL[nPos+ni] > dTemp)
					dTemp = pBL[nPos+ni]; */
				pMLike[nPos+ni] += dTDR*exp(dTemp);
				pPost[nPos+ni] -= 1.0;
			}
		}

		/* get next */
		pSite = pSite->pNext;
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_BGFit_UpdateLikelihood: update likelihood.                  */
/* ----------------------------------------------------------------------- */
int FlexModule_BGFit_ExportLogLikelihood(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, 
			int nRecordNum, char strFilePath[])
{
	/* define */
	FILE *fpOut;
	int ni,nLen;
	double *pPost,*pBL,*pML,*pMN;
	int numwritten;
	double dTemp,dTemp0;

	/* init */
	pBL = pSeqMtf->vMonitor[0]->pMatElement;
	pPost = pSeqMtf->vMonitor[1]->pMatElement;
	pML = pSeqMtf->vMonitor[2]->pMatElement;
	pMN = pSeqMtf->vMonitor[3]->pMatElement;
	nLen = pSeqMtf->nSeqLen;

	/* compute */
	/* fpOut = NULL;
	fpOut = fopen(strFilePath, "w");
	if(fpOut == NULL)
	{
		printf("Error: FlexModule_BGFit_ExportLogLikelihood, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	
	for(ni=0; ni<nLen; ni++)
	{
		pPost[ni] += nRecordNum;
		dTemp0 = (pML[ni]+exp(pBL[ni])*pPost[ni])/(pPost[ni]+pMN[ni]);
		dTemp0 = log(dTemp0);
		dTemp = (pML[ni]+exp(pBL[ni])*((double)nRecordNum-pMN[ni]))/(double)nRecordNum;
		dTemp = log(dTemp);
		fprintf(fpOut, "%f\t%f\t%f\t%f\t%f\t%f\t%f\n", pBL[ni], dTemp0, dTemp, pPost[ni], pMN[ni], pMN[ni]/(pPost[ni]+pMN[ni]), pMN[ni]/(double)nRecordNum);
		pPost[ni] = dTemp;
	}
	
	fclose(fpOut);
	*/

	for(ni=0; ni<nLen; ni++)
	{
		/* pPost[ni] += nRecordnum;
		pPost[ni] = (pML[ni]+exp(pBL[ni])*pPost[ni])/(pPost[ni]+pMN[ni]); */
		pPost[ni] = (pML[ni]+exp(pBL[ni])*((double)nRecordNum-pMN[ni]))/(double)nRecordNum;
		pPost[ni] = log(pPost[ni]);
	}

	fpOut = NULL;
	fpOut = fopen(strFilePath, "wb");
	if(fpOut == NULL)
	{
		printf("Error: FlexModule_BGFit_ExportLogLikelihood, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	
	numwritten = fwrite(&nLen, sizeof(int), 1, fpOut);
	numwritten = fwrite(pPost, sizeof(double), nLen, fpOut);
	if(numwritten != nLen)
	{
		printf("Error: coordinate wrong!\n");
		exit(EXIT_FAILURE);
	}

	fclose(fpOut);

	
	/* return */
	return PROC_SUCCESS;
}

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
					int nIndex, int nRecordNum)
{
	/* define */
	FILE *fpOut;
	int nBaseTypeNum;
	int nLen;
	int ni,nj,nk,nx;
	char strConsensus[LINE_LENGTH];
	struct DOUBLEMATRIX *pFreqPost;
	double dLowCut = 0.01;
	struct FLEXMOTIFSITE *pMtfSite,*pPrevSite,*pTempSite,*pSite;
	struct FLEXMOTIFSITE **vMotifList;
	unsigned char *pBase,*pBi;
	double *pOccur;
	int nMId;
	double dScore;
	int nMax;
	double dMax;
	double *pElement1,*pElement2;
	int nSeqId;
	int nId1,nId2;
	int nStateId,nStrand;
	double dTemp,dMaxPost;


	/* init */
	nBaseTypeNum = pBG0->nWidth;

	if(nRecordNum <= 0)
	{
		printf("Error: FlexModule_BGFit_WritePostMotifToFile, no posterior samples available!\n");
		exit(EXIT_FAILURE);
	}

	pFreqPost = NULL;
	pFreqPost = CreateDoubleMatrix(pFreqCount->nHeight, pFreqCount->nWidth); 
	if(pFreqPost == NULL)
	{
		printf("Error: FlexModule_BGFit_WritePostMotifToFile, cannot create posterior probability matrix!\n");
		exit(EXIT_FAILURE);
	}


	vMotifList = NULL;
	vMotifList = (struct FLEXMOTIFSITE **)calloc(nMotifNum, sizeof(struct FLEXMOTIFSITE *));
	if(vMotifList == NULL)
	{
		printf("Error: FlexModule_WritePostMotifToFile, cannot create motif site list!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nMotifNum; ni++)
	{
		vMotifList[ni] = NULL;
	}

	/* initial call of motif accoriding to posterior samples, */
	/* motifs will be aligned according to current matrix */
	for(ni=0; ni<nSeqCount; ni++)
	{
		/* init */
		pBase = vSeqMtf[ni]->vSeq[nIndex]->pMatElement;
		pOccur = vSeqMtf[ni]->vMonitor[4]->pMatElement;

		/* destroy old lists */
		FLEXMTFSDESTROYLIST(&(vSeqMtf[ni]->pMotifList));
		vSeqMtf[ni]->pMotifList = NULL;
		pPrevSite = NULL;

		/* get posterior probability */
		nLen = vSeqMtf[ni]->nSeqLen;
		for(nk=0; nk<nLen; nk++)
		{
			pOccur[nk] /= (double)nRecordNum;
		
			/* create new sites according to posterior probability */		
			if(pOccur[nk] > dLowCut)
			{
				/* get motif */
				nStateId = 0;
				nStrand = 0;
				dMaxPost = -1e100;
				for(nj=0; nj<nMotifNum; nj++)
				{
					for(nx=0; nx<2; nx++)
					{
						nId1 = 2+2*nj+nx;
						nId2 = 5+2*nj+nx;
						if(vSeqMtf[ni]->vStatus[nId1]->pMatElement[nk] == 0)
						{
							if(nUsePostEnrich == 1)
							{
								dTemp = log(DMGETAT(pFreqCount, (nj+1), nx))+vSeqMtf[ni]->vMonitor[nId2]->pMatElement[nk];
							}
							else
							{
								dTemp = log(DMGETAT(pFreqPrior, (nj+1), nx))+vSeqMtf[ni]->vMonitor[nId2]->pMatElement[nk];
							}
							
							if(dTemp > dMaxPost)
							{
								dMaxPost = dTemp;
								nStateId = nj+1;
								nStrand = nx;
							}
						}
					}
				}

				if(nStateId == 0)
					continue;

				/* create new site */
				pMtfSite = NULL;
				pMtfSite = FLEXMTFSCREATE();
				if(pMtfSite == NULL)
				{
					printf("Error: FlexModule_WritePostMotifToFile, cannot create space for new site!\n");
					exit(EXIT_FAILURE);
				}

				pMtfSite->dProb = pOccur[nk];
				pMtfSite->nMotifType = nStateId;
				pMtfSite->nSeqId = vSeqMtf[ni]->nId;
				pMtfSite->nStartPos = nk;
				pMtfSite->pNext = NULL;
				pMtfSite->pPrev = NULL;

				if(nStrand = 0)
				{
					pMtfSite->nStrand = 0;
					nId2 = 5+2*(nStateId-1);
					pMtfSite->dScore = vSeqMtf[ni]->vMonitor[nId2]->pMatElement[nk];
				}
				else
				{
					pMtfSite->nStrand = 1;
					nId2 = 5+2*(nStateId-1)+1;
					pMtfSite->dScore = vSeqMtf[ni]->vMonitor[nId2]->pMatElement[nk];
				}

				/* add to site list */
				if(vSeqMtf[ni]->pMotifList == NULL)
				{
					vSeqMtf[ni]->pMotifList = pMtfSite;
					pPrevSite = pMtfSite;
				}
				else
				{
					/* check if overlap. */ 
					nMId = pPrevSite->nMotifType-1;
							
					/* if not overlap, add to motif site list */
					if(pPrevSite->nStartPos+pMotifLen->pMatElement[nMId] <= pMtfSite->nStartPos)
					{
						pPrevSite->pNext = pMtfSite;
						pMtfSite->pPrev = pPrevSite;
						pPrevSite = pMtfSite;
					}

					/* If overlap, choose the best one */
					else
					{
						/* if different motif type */
						if(pPrevSite->nMotifType != pMtfSite->nMotifType)
						{
							if(pMtfSite->dProb <= pPrevSite->dProb)
							{
								FLEXMTFSDESTROY(pMtfSite);
							}
							else
							{
								pTempSite = pPrevSite;
								pPrevSite = pPrevSite->pPrev;
								FLEXMTFSDESTROY(pTempSite);
								if(pPrevSite == NULL)
								{
									vSeqMtf[ni]->pMotifList = pMtfSite;
									pPrevSite = pMtfSite;
								}
								else
								{
									pPrevSite->pNext = pMtfSite;
									pMtfSite->pPrev = pPrevSite;
									pPrevSite = pMtfSite;
								}
							}
						}
						/* if the same motif type */
						else
						{
							if(pMtfSite->dScore <= pPrevSite->dScore)
							{
								FLEXMTFSDESTROY(pMtfSite);
							}
							else
							{
								pTempSite = pPrevSite;
								pPrevSite = pPrevSite->pPrev;
								FLEXMTFSDESTROY(pTempSite);
								if(pPrevSite == NULL)
								{
									vSeqMtf[ni]->pMotifList = pMtfSite;
									pPrevSite = pMtfSite;
								}
								else
								{
									pPrevSite->pNext = pMtfSite;
									pMtfSite->pPrev = pPrevSite;
									pPrevSite = pMtfSite;
								}
							}
						}
					}
				}
			}
		}

		/* FOR TEST: */
		/* sprintf(strLine, "%s_post%d.txt", strFilePath, ni);
		fpOut = NULL;
		fpOut = fopen(strLine, "w");
		for(nk=0; nk<nLen; nk++)
		{
			for(nj=1; nj<=nMotifNum; nj++)
			{
				fprintf(fpOut, "%f\t", vSeqMtf[ni]->vMonitor[nj]->pMatElement[nk]);
			}
			fprintf(fpOut, "\n");
		}
		fclose(fpOut); */
	}
	

	/* sort motifs according to motif type and posterior probability */
	for(ni=0; ni<nSeqCount; ni++)
	{
		/* init */
		while(vSeqMtf[ni]->pMotifList != NULL)
		{
			pMtfSite = vSeqMtf[ni]->pMotifList;
			pTempSite = pMtfSite->pNext;
			if(pTempSite != NULL)
			{
				pTempSite->pPrev = NULL;
			}
			vSeqMtf[ni]->pMotifList = pTempSite;
			pMtfSite->pNext = NULL;
			nMId = pMtfSite->nMotifType-1;

			FLEXMTFSADDTOLIST_SORTBYPROB(vMotifList+nMId, pMtfSite);
		}
	}

	/* only keep the top N motifs */
	for(ni=0; ni<nMotifNum; ni++)
	{
		nj = 0;
		pMtfSite = vMotifList[ni];
		while(pMtfSite != NULL)
		{
			nj++;
			pMtfSite = pMtfSite->pNext;
			if(nj == (int)(pMotifSiteNum->pMatElement[ni]))
			{
				break;
			}
		}

		if(pMtfSite != NULL)
		{
			pPrevSite = pMtfSite->pPrev;
			if(pPrevSite != NULL)
			{
				pPrevSite->pNext = NULL;
			}
			else
			{
				vMotifList[ni] = NULL;
			}
			pMtfSite->pPrev = NULL;
			FLEXMTFSDESTROYLIST(&pMtfSite);
		}
	}

	/* resort sites by sequence id and position */
	for(ni=0; ni<nMotifNum; ni++)
	{
		pTempSite = NULL;
		while(vMotifList[ni] != NULL)
		{
			pMtfSite = vMotifList[ni];
			vMotifList[ni] = pMtfSite->pNext;
			pMtfSite->pNext = NULL;
			pMtfSite->pPrev = NULL;
			FLEXMTFSADDTOLIST_SORTBYPOS(&pTempSite, pMtfSite);
		}

		vMotifList[ni] = pTempSite;
	}

	/* export */
	fpOut = NULL;
	fpOut = fopen(strFilePath, "w");
	if(fpOut == NULL)
	{
		printf("Error: cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	/* write motif one by one */
	for(ni=0; ni<nMotifNum; ni++)
	{
		/* get motif score */
		dScore = FLEXMTFMSCORE(vMotif[ni], pBG0);

		/* write head */
		fprintf(fpOut, "****** Motif%d ******\n", ni);
		fprintf(fpOut, "Motif Score: %f\n", dScore);
		fprintf(fpOut, "Motif Matrix: \n");
		pElement1 = vMotif[ni]->pPriorCount->pMatElement;
		pElement2 = vMotif[ni]->pSampleCount->pMatElement;
		for(nj=0; nj<vMotif[ni]->pPriorCount->nHeight; nj++)
		{
			dMax = 0.0;
			nMax = nBaseTypeNum;
			for(nk=0; nk<vMotif[ni]->pPriorCount->nWidth; nk++)
			{
				dTemp = (*pElement1)+(*pElement2);
				fprintf(fpOut, "% 9.7e ", dTemp);
				if(dTemp > dMax)
				{
					dMax = dTemp;
					nMax = nk;
				}
				pElement1++;
				pElement2++;
			}

			fprintf(fpOut, "\n");
			
			switch(nMax)
			{
				case 0: strConsensus[nj] = 'A';
					break;
				case 1: strConsensus[nj] = 'C';
					break;
				case 2: strConsensus[nj] = 'G';
					break;
				case 3: strConsensus[nj] = 'T';
					break;
				default: strConsensus[nj] = 'N';
			}
		}
		strConsensus[nj] = '\0';

		fprintf(fpOut, "\n");
		fprintf(fpOut, "Consensus:\n");
		fprintf(fpOut, "%s\n", strConsensus);

		fprintf(fpOut, "\nMotif Sites:\n");
		pSite = vMotifList[ni];
		while(pSite != NULL)
		{
			/* head information */
			nSeqId = pSite->nSeqId;
			nMId = pSite->nMotifType-1;
			fprintf(fpOut, "%d\t%d\t%d\t", pSite->nSeqId, 
				pSite->nStartPos, (pSite->nStartPos+pMotifLen->pMatElement[nMId]-1));
			if(pSite->nStrand == 0)
			{
				fprintf(fpOut, "+\t");
			}
			else if(pSite->nStrand == 1)
			{
				fprintf(fpOut, "-\t");
			}
			else
			{
				fprintf(fpOut, "?\t");
			}

			/* write seq */
			if(pSite->nStrand == 0)
			{
				pBi = vSeqMtf[nSeqId]->vSeq[nIndex]->pMatElement + pSite->nStartPos;
				for(nj=0; nj<pMotifLen->pMatElement[nMId]; nj++)
				{
					switch((int)(*pBi))
					{
						case 0: fprintf(fpOut, "A");
							break;
						case 1: fprintf(fpOut, "C");
							break;
						case 2: fprintf(fpOut, "G");
							break;
						case 3: fprintf(fpOut, "T");
							break;
						case 10: fprintf(fpOut, "a");
							break;
						case 11: fprintf(fpOut, "c");
							break;
						case 12: fprintf(fpOut, "g");
							break;
						case 13: fprintf(fpOut, "t");
							break;
						default: fprintf(fpOut, "N");
					}
					pBi++;
				}
			}
			else if(pSite->nStrand == 1)
			{
				pBi = vSeqMtf[nSeqId]->vSeq[nIndex]->pMatElement + pSite->nStartPos + pMotifLen->pMatElement[nMId] - 1;
				for(nj=0; nj<pMotifLen->pMatElement[nMId]; nj++)
				{
					switch((int)(*pBi))
					{
						case 0: fprintf(fpOut, "T");
							break;
						case 1: fprintf(fpOut, "G");
							break;
						case 2: fprintf(fpOut, "C");
							break;
						case 3: fprintf(fpOut, "A");
							break;
						case 10: fprintf(fpOut, "t");
							break;
						case 11: fprintf(fpOut, "g");
							break;
						case 12: fprintf(fpOut, "c");
							break;
						case 13: fprintf(fpOut, "a");
							break;
						default: fprintf(fpOut, "N");
					}
					pBi--;
				}
			}
			fprintf(fpOut, "\n");

			/* get next */
			pSite = pSite->pNext;
		}

		fprintf(fpOut, "\n\n");
	}

	fclose(fpOut);

	/* release memory */
	DestroyDoubleMatrix(pFreqPost);
	for(ni=0; ni<nMotifNum; ni++)
	{
		FLEXMTFSDESTROYLIST(vMotifList+ni);
	}
	free(vMotifList);


	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  FlexModule_WriteCisGenomeIniFiles: prepare cisgenome ini file.         */
/* ----------------------------------------------------------------------- */ 
int FlexModule_WriteCisGenomeIniFiles(char strParamFile[], char strWorkPath[],
							char strOutFile[], int nMotifNum)
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	int ni;
	char strFileName[MED_LINE_LENGTH];
	char strLOutName[MED_LINE_LENGTH];
	char strTemp[MED_LINE_LENGTH];
	char strLine[MED_LINE_LENGTH];

	/* write individual motif file */
	sprintf(strFileName, "%s%s_p.txt", strWorkPath, strOutFile);
	fpIn = NULL;
	fpIn = fopen(strFileName, "r");
	if(fpIn == NULL)
	{
		printf("Error: FlexModule_WriteCisGenomeIniFiles, cannot open result file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strstr(strLine, "Motif Matrix:") == strLine)
		{
			sprintf(strFileName, "%s%s_%d.mat", strWorkPath, strOutFile, ni);
			fpOut = NULL;
			fpOut = fopen(strFileName, "w");
			if(fpOut == NULL)
			{
				printf("Error: FlexModule_WriteCisGenomeIniFiles, cannot open result file!\n");
				exit(EXIT_FAILURE);
			}

			while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
			{
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] == '\0')
					break;

				fprintf(fpOut, "%s\n", strLine);
			}

			fclose(fpOut);
			ni++;
		}
	}

	if(ni != nMotifNum)
	{
		printf("Error: FlexModule_WriteCisGenomeIniFiles, motif number not consistent!\n");
		exit(EXIT_FAILURE);
	}

	fclose(fpIn);

	/* write motif list file */
	sprintf(strFileName, "%s%s.matl", strWorkPath, strOutFile);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: FlexModule_WriteCisGenomeIniFiles, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nMotifNum; ni++)
	{
		fprintf(fpOut, "%s%s_%d.mat\n", strWorkPath, strOutFile, ni);
	}
	
	fclose(fpOut);

	/* write _l.txt mat */
	sprintf(strFileName, "%s%s_l.txt", strWorkPath, strOutFile);
	sprintf(strLOutName, "%s%s_l", strWorkPath, strOutFile);
	FlexModule_ExtractMotifsFromResult(strFileName, strLOutName);

	/* write ini file */
	sprintf(strFileName, "%s%s.cgw", strWorkPath, strOutFile);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: FlexModule_WriteCisGenomeIniFiles, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	/* GetFileName(strParamFile, strTemp); */
	fprintf(fpOut, "[item1]\n");
	fprintf(fpOut, "type=txt\n");
	fprintf(fpOut, "file_path=%s\n", strParamFile);
	
	sprintf(strTemp, "%s%s_p.txt", strWorkPath, strOutFile);
	fprintf(fpOut, "[item2]\n");
	fprintf(fpOut, "type=txt\n");
	fprintf(fpOut, "file_path=%s\n", strTemp);

	sprintf(strTemp, "%s%s.matl", strWorkPath, strOutFile);
	fprintf(fpOut, "[item3]\n");
	fprintf(fpOut, "type=matl\n");
	fprintf(fpOut, "file_path=%s\n", strTemp);

	for(ni=0; ni<nMotifNum; ni++)
	{
		sprintf(strTemp, "%s%s_%d.mat", strWorkPath, strOutFile, ni);
		fprintf(fpOut, "[item%d]\n", 4+ni);
		fprintf(fpOut, "type=mat\n");
		fprintf(fpOut, "file_path=%s\n", strTemp);
	}
	
	sprintf(strTemp, "%s%s_l.txt", strWorkPath, strOutFile);
	fprintf(fpOut, "[item%d]\n", 4+nMotifNum);
	fprintf(fpOut, "type=txt\n");
	fprintf(fpOut, "file_path=%s\n", strTemp);

	sprintf(strTemp, "%s%s_l.matl", strWorkPath, strOutFile);
	fprintf(fpOut, "[item%d]\n", 5+nMotifNum);
	fprintf(fpOut, "type=matl\n");
	fprintf(fpOut, "file_path=%s\n", strTemp);

	for(ni=0; ni<nMotifNum; ni++)
	{
		sprintf(strTemp, "%s%s_l_%d.mat", strWorkPath, strOutFile, ni);
		fprintf(fpOut, "[item%d]\n", 6+nMotifNum+ni);
		fprintf(fpOut, "type=mat\n");
		fprintf(fpOut, "file_path=%s\n", strTemp);
	}
	
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  FlexModule_ExtractMotifsFromResult: prepare cisgenome ini file.        */
/* ----------------------------------------------------------------------- */ 
int FlexModule_ExtractMotifsFromResult(char strInFile[], char strOutFile[])
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	int ni;
	char strFileName[MED_LINE_LENGTH];
	char strTemp[MED_LINE_LENGTH];
	char strLine[MED_LINE_LENGTH];
	int nMotifNum;

	/* write individual motif file */
	fpIn = NULL;
	fpIn = fopen(strInFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: FlexModule_ExtractMotifsFromResult, cannot open result file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strstr(strLine, "Motif Matrix:") == strLine)
		{
			sprintf(strFileName, "%s_%d.mat", strOutFile, ni);
			fpOut = NULL;
			fpOut = fopen(strFileName, "w");
			if(fpOut == NULL)
			{
				printf("Error: FlexModule_ExtractMotifsFromResult, cannot open result file!\n");
				exit(EXIT_FAILURE);
			}

			while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
			{
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] == '\0')
					break;

				fprintf(fpOut, "%s\n", strLine);
			}

			fclose(fpOut);
			ni++;
		}
	}


	fclose(fpIn);
	nMotifNum = ni;

	/* write motif list file */
	sprintf(strFileName, "%s.matl", strOutFile);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: FlexModule_ExtractMotifsFromResult, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nMotifNum; ni++)
	{
		fprintf(fpOut, "%s_%d.mat\n", strOutFile, ni);
	}
	
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}

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
					int nIncludeRepeat)
{
	/* sequences and conservation */
	/* sequences */
	struct FLEXSEQMOTIF *pSeqMtf = NULL;
	/* sequence number */
	int nSeqCount = 0;
	/* site number */
	int nEffecLen = 0;
	int nConsEffecLen = 0;
	int nSiteNum = 0;
	/* total sequence length and total site number */
	double dTotLen = 0.0;
	double dTotConsLen = 0.0;
	double dTotSite = 0.0;

	/* segment position */
	int nM1,nM2;
	int nActualStart,nActualEnd;
	
	/* motif consensus */
	FILE *fpIn;
	FILE *fpOut;
	FILE *fpOut2;
	char strLine[LONG_LINE_LENGTH];
	struct BYTEMATRIX *pCMotif = NULL;
	struct BYTEMATRIX *pDMotif = NULL;
	struct BYTEMATRIX *pTempCMotif = NULL;
	struct BYTEMATRIX *pTempDMotif = NULL;
	int nMotifLen;
	unsigned char bStatus;
	
	/* number of base types */
	int nBaseTypeNum = 4;
	
	/* load sequence */
	char strSeqAlias[LINE_LENGTH];
	char strChr[LINE_LENGTH];
	int nStart,nEnd;
		
	/* other variables */
	int nlen;
	int ni,nj;
	int nIsDegenerate;
	
	/* #################################### */
	/* initialize                           */
	/* #################################### */
	
	/* #################################### */
	/* load initial parameters              */
	/* #################################### */

	if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
	{
		nlen = strlen(strGenomePath);
		if(nlen == 0)
		{
			sprintf(strGenomePath, ".\\"); 
		}
		else if(strGenomePath[nlen-1] != '\\')
		{
			strGenomePath[nlen] = '\\';
			strGenomePath[nlen+1] = '\0';
		}

		if(nUseCS == 1)
		{
			nlen = strlen(strCSPath);
			if(nlen == 0)
			{
				sprintf(strCSPath, ".\\"); 
			}
			else if(strCSPath[nlen-1] != '\\')
			{
				strCSPath[nlen] = '\\';
				strCSPath[nlen+1] = '\0';
			}
		}
	}
	else
	{
		nlen = strlen(strGenomePath);
		if(nlen == 0)
		{
			sprintf(strGenomePath, "./"); 
		}
		else if(strGenomePath[nlen-1] != '/')
		{
			strGenomePath[nlen] = '/';
			strGenomePath[nlen+1] = '\0';
		}

		if(nUseCS == 1)
		{
			nlen = strlen(strCSPath);
			if(nlen == 0)
			{
				sprintf(strCSPath, "./"); 
			}
			else if(strCSPath[nlen-1] != '/')
			{
				strCSPath[nlen] = '/';
				strCSPath[nlen+1] = '\0';
			}
		}
	}

	/* #################################### */
	/* load motif                           */
	/* #################################### */
	fpIn = NULL;
	fpIn = fopen(strMotifPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: MotifMap_ScanConsensus_Genome_Main, cannot load motif consensus pattern!\n");
		exit(EXIT_FAILURE);
	}
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] != '\0')
			break;
	}
	fclose(fpIn);

	StrMakeUpper(strLine);
	nlen = strlen(strLine);
	pTempCMotif = CreateByteMatrix(nlen, 1);
	pTempDMotif = CreateByteMatrix(nlen, nBaseTypeNum);
	nMotifLen = 0;
	nIsDegenerate = 0;
	for(ni=0; ni<nlen; ni++)
	{
		if((nIsDegenerate > 1) || (nIsDegenerate < 0))
		{
			printf("Error: MotifMap_ScanConsensus_Genome_Main, number of  '[' and ']' does not match!\n");
			exit(EXIT_FAILURE);
		}

		if(strLine[ni] == '[')
		{
			nIsDegenerate += 1;
			nMotifLen--;
		}
		else if(strLine[ni] == ']')
		{
			nIsDegenerate -= 1;
			nMotifLen++;
		}
		else if(strLine[ni] == 'A')
		{
			BMSETAT(pTempDMotif, nMotifLen, 0, 1);
			if(nIsDegenerate == 0)
			{
				pTempCMotif->pMatElement[nMotifLen] = 0;
				nMotifLen++;
			}
		}
		else if(strLine[ni] == 'C')
		{
			BMSETAT(pTempDMotif, nMotifLen, 1, 1);
			if(nIsDegenerate == 0)
			{
				pTempCMotif->pMatElement[nMotifLen] = 1;
				nMotifLen++;
			}
		}
		else if(strLine[ni] == 'G')
		{
			BMSETAT(pTempDMotif, nMotifLen, 2, 1);
			if(nIsDegenerate == 0)
			{
				pTempCMotif->pMatElement[nMotifLen] = 2;
				nMotifLen++;
			}
		}
		else if(strLine[ni] == 'T')
		{
			BMSETAT(pTempDMotif, nMotifLen, 3, 1);
			if(nIsDegenerate == 0)
			{
				pTempCMotif->pMatElement[nMotifLen] = 3;
				nMotifLen++;
			}
		}
		else if(strLine[ni] == 'N')
		{
			for(nj=0; nj<nBaseTypeNum; nj++)
			{
				BMSETAT(pTempDMotif, nMotifLen, nj, 1);
			}
			
			if(nIsDegenerate == 0)
			{
				pTempCMotif->pMatElement[nMotifLen] = nBaseTypeNum;
				nMotifLen++;
			}
		}
		else if(strLine[ni] == ',')
		{
		}
		else if(strLine[ni] == ' ')
		{
		}
		else if(strLine[ni] == '\t')
		{
		}
		else
		{
			printf("Error: MotifMap_ScanConsensus_Genome_Main, %c is not accepted to represent consensus pattern!\n", strLine[ni]);
			exit(EXIT_FAILURE);
		}
	}

	if(nIsDegenerate != 0)
	{
		printf("Error: MotifMap_ScanConsensus_Genome_Main, number of  '[' and ']' does not match!\n");
		exit(EXIT_FAILURE);
	}

	pCMotif = CreateByteMatrix(nMotifLen, 1);
	pDMotif = CreateByteMatrix(nMotifLen, nBaseTypeNum);
	for(ni=0; ni<nMotifLen; ni++)
	{
		pCMotif->pMatElement[ni] = pTempCMotif->pMatElement[ni];
		for(nj=0; nj<nBaseTypeNum; nj++)
		{
			bStatus = BMGETAT(pTempDMotif, ni, nj);
			BMSETAT(pDMotif, ni, nj, bStatus);
		}
	}

	DestroyByteMatrix(pTempCMotif);
	DestroyByteMatrix(pTempDMotif);
	/* 
	BMSAVE(pCMotif, "testMC.txt");
	BMSAVE(pDMotif, "testMD.txt");
	*/


	/* #################################### */
	/* load sequences and conservation      */
	/* #################################### */
	fpIn = NULL;
	fpIn = fopen(strCodPath, "r");
	if(fpIn == NULL)
	{
        printf("Error: MotifMap_ScanConsensus_Genome_Main, cannot open coordinates file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
        printf("Error: MotifMap_ScanConsensus_Genome_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	fprintf(fpOut, "#seq_id\tchromosome\tstart\tend\tstrand\tscore\tsite\n");
	
	sprintf(strLine, "%s.stat", strOutputPath);
	fpOut2 = NULL;
	fpOut2 = fopen(strLine, "w");
	if(fpOut2 == NULL)
	{
        printf("Error: MotifMap_ScanConsensus_Genome_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	/* read */
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		sscanf(strLine, "%s %s %d %d", strSeqAlias, strChr,
			&nStart, &nEnd);
		
		nM1 = nStart;
		while(nM1 <= nEnd)
		{
			nM2 = nM1+GENOME_CONTIG_LEN-1;
			if(nM2 > nEnd)
				nM2 = nEnd;

			nActualStart = nM1;
			nActualEnd = nM2+nMotifLen-1;
			if(nActualEnd > nEnd)
				nActualEnd = nEnd;

			/* #################################### */
			/* load sequence and score              */
			/* #################################### */
			pSeqMtf = NULL;
			pSeqMtf = FLEXSEQMTFLOADSEQ_FOR_SCANCONSENSUS_GENOME(nSeqCount,
				strGenomePath, strChr, nActualStart, nActualEnd, 
				nUseCS, strCSPath);

			if(pSeqMtf == NULL)
			{
				nM1 = nM2+1;
				continue;
			}

			if(nIncludeRepeat == 1)
			{
				FLEXSEQMTF_REPEATUNMASK(pSeqMtf, 0);
			}

			/* #################################### */
			/* map consensus                        */
			/* #################################### */
			MotifMap_ScanConsensus_In_FlexSeqMtf(pSeqMtf, 0,
					nMotifLen, pCMotif, pDMotif, nMC, nMD, 
					nUseCS, dC, &nEffecLen, &nConsEffecLen, &nSiteNum);
			dTotLen += nEffecLen;
			dTotConsLen += nConsEffecLen;
			dTotSite += nSiteNum;

			/* #################################### */
			/* export results                       */
			/* #################################### */
			MotifMap_ScanConsensus_Export_Genome(pSeqMtf, 0,
					nMotifLen, strSeqAlias, strChr, nActualStart, fpOut);

			/* #################################### */
			/* release memory                       */
			/* #################################### */
			FLEXSEQMTFDESTROY(pSeqMtf);
			nM1 = nM2+1;
		}
		nSeqCount++;
	}

	/* save statistics */
	fprintf(fpOut2, "EffecLen= %d\n", (int)dTotLen);
	fprintf(fpOut2, "TotalSite= %d\n", (int)dTotSite);
	fprintf(fpOut2, "ConsLen= %d\n", (int)dTotConsLen);
	
	/* close files */
	fclose(fpIn);
	fclose(fpOut);
	fclose(fpOut2);

	/* #################################### */
	/* release memory                       */
	/* #################################### */
	DestroyByteMatrix(pCMotif);
	DestroyByteMatrix(pDMotif);

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanConsensus_Export_Genome: save motif sites to a file.      */
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanConsensus_Export_Genome(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,
				int nMotifLen, char strSeqAlias[], char strChr[], int nOffset,
				FILE *fpOut)
{
	/* define */
	int nj;
	struct FLEXMOTIFSITE *pSite;
	char chStrand;
	char chBase;
	unsigned char *pBase;

	/* export */
	pSite = pSeqMtf->pMotifList;
	while(pSite != NULL)
	{
		if(pSite->nStrand == 0)
			chStrand = '+';
		else
			chStrand = '-';

		/* fprintf(fpOut, ">%d:%d-%d(%c)\n", pSeqMtf->nId, pSite->nStartPos, pSite->nStartPos+nMotifLen-1,
			chStrand); */

		fprintf(fpOut, "%s\t%s\t%d\t%d\t%c\t%f\t", strSeqAlias, strChr, 
			nOffset+pSite->nStartPos, nOffset+pSite->nStartPos+nMotifLen-1,
			chStrand, -(pSite->dScore));

		if(pSite->nStrand == 0)
		{
			pBase = pSeqMtf->vSeq[nIndex]->pMatElement+pSite->nStartPos;
			for(nj=0; nj<nMotifLen; nj++)
			{
				switch(*pBase)
				{
					case 0: chBase = 'A';
						break;
					case 1: chBase = 'C';
						break;
					case 2: chBase = 'G';
						break;
					case 3: chBase = 'T';
						break;
					case 10: chBase = 'a';
						break;
					case 11: chBase = 'c';
						break;
					case 12: chBase = 'g';
						break;
					case 13: chBase = 't';
						break;
					default: chBase = 'N';
				}

				fprintf(fpOut, "%c", chBase);
				pBase++;
			}
			fprintf(fpOut, "\n");
		}
		else
		{
			pBase = pSeqMtf->vSeq[nIndex]->pMatElement+pSite->nStartPos+nMotifLen-1;
			for(nj=0; nj<nMotifLen; nj++)
			{
				switch(*pBase)
				{
					case 0: chBase = 'T';
						break;
					case 1: chBase = 'G';
						break;
					case 2: chBase = 'C';
						break;
					case 3: chBase = 'A';
						break;
					case 10: chBase = 't';
						break;
					case 11: chBase = 'g';
						break;
					case 12: chBase = 'c';
						break;
					case 13: chBase = 'a';
						break;
					default: chBase = 'N';
				}

				fprintf(fpOut, "%c", chBase);
				pBase--;
			}
			fprintf(fpOut, "\n");
		}

		pSite = pSite->pNext;
	}

	/* return */
	return PROC_SUCCESS;
}

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
					int nMC, int nMD, int nUseCS, double dC, char strCSPrefix[])
{
	/* sequences and conservation */
	/* sequences */
	struct FLEXSEQMOTIF *pSeqMtf = NULL;
	/* sequence number */
	int nSeqCount = 0;
	/* site number */
	int nEffecLen = 0;
	int nConsEffecLen = 0;
	int nSiteNum = 0;
	/* total sequence length and total site number */
	double dTotLen = 0.0;
	double dTotConsLen = 0.0;
	double dTotSite = 0.0;
	
	/* motif consensus */
	FILE *fpIn;
	FILE *fpOut;
	FILE *fpOut2;
	char strLine[MED_LINE_LENGTH];
	struct BYTEMATRIX *pCMotif = NULL;
	struct BYTEMATRIX *pDMotif = NULL;
	struct BYTEMATRIX *pTempCMotif = NULL;
	struct BYTEMATRIX *pTempDMotif = NULL;
	int nMotifLen;
	unsigned char bStatus;
	
	/* number of base types */
	int nBaseTypeNum = 4;
	
	/* load sequence */
	char strSeqLine[LONG_LINE_LENGTH];
	struct tagSequence *pNewSeq;
	char strSeqAlias[LINE_LENGTH];
	
	/* other variables */
	int nlen;
	int ni,nj;
	int nIsDegenerate;
	
	/* #################################### */
	/* initialize                           */
	/* #################################### */
	
	/* #################################### */
	/* load initial parameters              */
	/* #################################### */

	if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
	{
		nlen = strlen(strInputPath);
		if(nlen == 0)
		{
			sprintf(strInputPath, ".\\"); 
		}
		else if(strInputPath[nlen-1] != '\\')
		{
			strInputPath[nlen] = '\\';
			strInputPath[nlen+1] = '\0';
		}
	}
	else
	{
		nlen = strlen(strInputPath);
		if(nlen == 0)
		{
			sprintf(strInputPath, "./"); 
		}
		else if(strInputPath[nlen-1] != '/')
		{
			strInputPath[nlen] = '/';
			strInputPath[nlen+1] = '\0';
		}
	}

	/* #################################### */
	/* load motif                           */
	/* #################################### */
	fpIn = NULL;
	fpIn = fopen(strMotifPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: MotifMap_ScanConsensus_Main, cannot load motif consensus pattern!\n");
		exit(EXIT_FAILURE);
	}
	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] != '\0')
			break;
	}
	fclose(fpIn);

	StrMakeUpper(strLine);
	nlen = strlen(strLine);
	pTempCMotif = CreateByteMatrix(nlen, 1);
	pTempDMotif = CreateByteMatrix(nlen, nBaseTypeNum);
	nMotifLen = 0;
	nIsDegenerate = 0;
	for(ni=0; ni<nlen; ni++)
	{
		if((nIsDegenerate > 1) || (nIsDegenerate < 0))
		{
			printf("Error: MotifMap_ScanConsensus_Main, number of  '[' and ']' does not match!\n");
			exit(EXIT_FAILURE);
		}

		if(strLine[ni] == '[')
		{
			nIsDegenerate += 1;
			nMotifLen--;
		}
		else if(strLine[ni] == ']')
		{
			nIsDegenerate -= 1;
			nMotifLen++;
		}
		else if(strLine[ni] == 'A')
		{
			BMSETAT(pTempDMotif, nMotifLen, 0, 1);
			if(nIsDegenerate == 0)
			{
				pTempCMotif->pMatElement[nMotifLen] = 0;
				nMotifLen++;
			}
		}
		else if(strLine[ni] == 'C')
		{
			BMSETAT(pTempDMotif, nMotifLen, 1, 1);
			if(nIsDegenerate == 0)
			{
				pTempCMotif->pMatElement[nMotifLen] = 1;
				nMotifLen++;
			}
		}
		else if(strLine[ni] == 'G')
		{
			BMSETAT(pTempDMotif, nMotifLen, 2, 1);
			if(nIsDegenerate == 0)
			{
				pTempCMotif->pMatElement[nMotifLen] = 2;
				nMotifLen++;
			}
		}
		else if(strLine[ni] == 'T')
		{
			BMSETAT(pTempDMotif, nMotifLen, 3, 1);
			if(nIsDegenerate == 0)
			{
				pTempCMotif->pMatElement[nMotifLen] = 3;
				nMotifLen++;
			}
		}
		else if(strLine[ni] == 'N')
		{
			for(nj=0; nj<nBaseTypeNum; nj++)
			{
				BMSETAT(pTempDMotif, nMotifLen, nj, 1);
			}
			
			if(nIsDegenerate == 0)
			{
				pTempCMotif->pMatElement[nMotifLen] = nBaseTypeNum;
				nMotifLen++;
			}
		}
		else if(strLine[ni] == ',')
		{
		}
		else if(strLine[ni] == ' ')
		{
		}
		else if(strLine[ni] == '\t')
		{
		}
		else
		{
			printf("Error: MotifMap_ScanConsensus_Main, %c is not accepted to represent consensus pattern!\n", strLine[ni]);
			exit(EXIT_FAILURE);
		}
	}

	if(nIsDegenerate != 0)
	{
		printf("Error: MotifMap_ScanConsensus_Main, number of  '[' and ']' does not match!\n");
		exit(EXIT_FAILURE);
	}

	pCMotif = CreateByteMatrix(nMotifLen, 1);
	pDMotif = CreateByteMatrix(nMotifLen, nBaseTypeNum);
	for(ni=0; ni<nMotifLen; ni++)
	{
		pCMotif->pMatElement[ni] = pTempCMotif->pMatElement[ni];
		for(nj=0; nj<nBaseTypeNum; nj++)
		{
			bStatus = BMGETAT(pTempDMotif, ni, nj);
			BMSETAT(pDMotif, ni, nj, bStatus);
		}
	}

	DestroyByteMatrix(pTempCMotif);
	DestroyByteMatrix(pTempDMotif);
	/* 
	BMSAVE(pCMotif, "testMC.txt");
	BMSAVE(pDMotif, "testMD.txt");
	*/


	/* #################################### */
	/* load sequences and conservation      */
	/* #################################### */
	sprintf(strLine, "%s%s", strInputPath, strSeqFile);
	fpIn = NULL;
	fpIn = fopen(strLine, "r");
	if(fpIn == NULL)
	{
        printf("Error: MotifMap_ScanConsensus_Main, cannot open sequence file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
        printf("Error: MotifMap_ScanConsensus_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	fprintf(fpOut, "#seq_name\tseq_id\tstart\tend\tstrand\tscore\tsite\n");
	
	sprintf(strLine, "%s.stat", strOutputPath);
	fpOut2 = NULL;
	fpOut2 = fopen(strLine, "w");
	if(fpOut2 == NULL)
	{
        printf("Error: MotifMap_ScanConsensus_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	/* read */
	pNewSeq = NULL;
	while(fgets(strSeqLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strSeqLine);
		StrTrimRight(strSeqLine);
		if(strSeqLine[0] == '\0')
			continue;

		if(strSeqLine[0] == '>')
		{
			if(pNewSeq != NULL)
			{
				/* #################################### */
				/* load sequence and score              */
				/* #################################### */
				pSeqMtf = NULL;
				pSeqMtf = FLEXSEQMTFLOADSEQ_FOR_SCANCONSENSUS_SINGLE(pNewSeq,
					strInputPath, nUseCS, strCSPrefix);
				SequenceDelete(pNewSeq);

				/* #################################### */
				/* map consensus                        */
				/* #################################### */
				MotifMap_ScanConsensus_In_FlexSeqMtf(pSeqMtf, 0,
						nMotifLen, pCMotif, pDMotif, nMC, nMD, 
						nUseCS, dC, &nEffecLen, &nConsEffecLen, &nSiteNum);
				dTotLen += nEffecLen;
				dTotConsLen += nConsEffecLen;
				dTotSite += nSiteNum;

				/* #################################### */
				/* export results                       */
				/* #################################### */
				MotifMap_ScanConsensus_Export_Single(pSeqMtf, 0,
						nMotifLen, strSeqAlias, fpOut);

				/* #################################### */
				/* release memory                       */
				/* #################################### */
				FLEXSEQMTFDESTROY(pSeqMtf);
				nSeqCount++;
			}

			
			pNewSeq = NULL;
			pNewSeq = SequenceCreate();
			if(pNewSeq == NULL)
			{
				fclose(fpIn);
				printf("Error: MotifMap_ScanConsensus_Main, cannot allocate memory for loading sequences!\n");
				exit(EXIT_FAILURE);
			}
			pNewSeq->m_nIndex = nSeqCount;
			strcpy(pNewSeq->m_strAlias, (strSeqLine+1));
			strcpy(strSeqAlias, (strSeqLine+1));
		}
		else
		{
			if(SequenceAddTail(pNewSeq, strSeqLine) == PROC_FAILURE)
			{
				fclose(fpIn);
				printf("Error: MotifMap_ScanConsensus_Main, cannot allocate memory for loading sequences!\n");
				exit(EXIT_FAILURE);
			}			
		}
	}

	/* the last sequence */
	if(pNewSeq != NULL)
	{
		/* #################################### */
		/* load sequence and score              */
		/* #################################### */
		pSeqMtf = NULL;
		pSeqMtf = FLEXSEQMTFLOADSEQ_FOR_SCANCONSENSUS_SINGLE(pNewSeq,
			strInputPath, nUseCS, strCSPrefix);
		SequenceDelete(pNewSeq);

		/* #################################### */
		/* map consensus                        */
		/* #################################### */
		MotifMap_ScanConsensus_In_FlexSeqMtf(pSeqMtf, 0,
				nMotifLen, pCMotif, pDMotif, nMC, nMD, 
				nUseCS, dC, &nEffecLen, &nConsEffecLen, &nSiteNum);
		dTotLen += nEffecLen;
		dTotConsLen += nConsEffecLen;
		dTotSite += nSiteNum;

		/* #################################### */
		/* export results                       */
		/* #################################### */
		MotifMap_ScanConsensus_Export_Single(pSeqMtf, 0,
				nMotifLen, strSeqAlias, fpOut);

		/* #################################### */
		/* release memory                       */
		/* #################################### */
		FLEXSEQMTFDESTROY(pSeqMtf);
		nSeqCount++;
	}

	/* save statistics */
	fprintf(fpOut2, "EffecLen= %d\n", (int)dTotLen);
	fprintf(fpOut2, "TotalSite= %d\n", (int)dTotSite);
	fprintf(fpOut2, "ConsLen= %d\n", (int)dTotConsLen);
	
	/* close files */
	fclose(fpIn);
	fclose(fpOut);
	fclose(fpOut2);

	/* #################################### */
	/* release memory                       */
	/* #################################### */
	DestroyByteMatrix(pCMotif);
	DestroyByteMatrix(pDMotif);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanConsensus_In_FlexSeqMtf: call motif sites.                */
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanConsensus_In_FlexSeqMtf(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,
			int nMotifLen, struct BYTEMATRIX *pCMotif, struct BYTEMATRIX *pDMotif, 
			int nMC, int nMD, int nUseCS, double dC,
			int *pEffecLen, int *pConsEffecLen, int *pSiteNum)
{
	/* define */
	int nBaseTypeNum;
	int ni,nj,nLen;
	unsigned char *pBase;
	unsigned char *pCS;
	struct FLEXMOTIFSITE *pPrev,*pSite;
	
	int nPass[2];
	int nMisC[2];
	int nMisD[2];
	int nMask[2];
	double dTotConserve;
	double dAveConserve;
	int nStrand;
	double dScore;

	/* init */
	nBaseTypeNum = pDMotif->nWidth;
	*pEffecLen = 0;
	*pConsEffecLen = 0;
	*pSiteNum = 0;
	pPrev = NULL;

	nLen = pSeqMtf->vSeq[nIndex]->nWidth-nMotifLen+1;
	if(nLen <= 0)
	{
		return PROC_SUCCESS;
	}
	pBase = pSeqMtf->vSeq[nIndex]->pMatElement;
	dTotConserve = 0.0;
	if(nUseCS == 1)
	{
		pCS = pSeqMtf->vScore[0]->pMatElement;
		for(nj=0; nj<(nMotifLen-1); nj++)
		{
			dTotConserve += pCS[nj];
		}
	}

	for(ni=0; ni<nLen; ni++)
	{
		if(nUseCS == 1)
		{
			dTotConserve += pCS[ni+nMotifLen-1];
		}

		/* if repeat or 'N', skip it */
		if(pBase[ni] >= nBaseTypeNum)
		{
			if(nUseCS == 1)
			{
				dTotConserve -= pCS[ni];
			}
			continue;
		}
		
		/* add effect length */
		*pEffecLen += 1;
		
		/* compute conservation */
		if(nUseCS == 1)
		{
			dAveConserve = dTotConserve/(double)nMotifLen;
			if(dAveConserve < dC)
			{
				dTotConserve -= pCS[ni];
				continue;
			}
			else
			{
				*pConsEffecLen += 1;
			}
		}
		else
		{
			*pConsEffecLen += 1;
		}

		/* get motif matching score */
		MotifMap_ScanConsensus_MotifCountMisMatch(pSeqMtf, nIndex, ni, 
					'+', nMotifLen, pCMotif, pDMotif,
					(nMisC+0), (nMisD+0), (nMask+0));
		MotifMap_ScanConsensus_MotifCountMisMatch(pSeqMtf, nIndex, ni, 
					'-', nMotifLen, pCMotif, pDMotif,
					(nMisC+1), (nMisD+1), (nMask+1));
		
		if( (nMask[0]==0) && (nMisC[0]<=nMC) && (nMisD[0]<=nMD))
		{
			nPass[0] = 1;
		}
		else
		{
			nPass[0] = 0;
		}

		if( (nMask[1]==0) && (nMisC[1]<=nMC) && (nMisD[1]<=nMD))
		{
			nPass[1] = 1;
		}
		else
		{
			nPass[1] = 0;
		}

		if( (nPass[0]==1) && (nPass[1]==1) )
		{
			if(nMisD[0] < nMisD[1])
			{
				nStrand = 0;
				dScore = nMisD[0]+0.1*nMisC[0];
			}
			else if(nMisD[1] < nMisD[0])
			{
				nStrand = 1;
				dScore = nMisD[1]+0.1*nMisC[1]; 
			}
			else
			{
				if(nMisC[0] <= nMisC[1])
				{
					nStrand = 0;
					dScore = nMisD[0]+0.1*nMisC[0];
				}
				else
				{
					nStrand = 1;
					dScore = nMisD[1]+0.1*nMisC[1];
				}
			}
		}
		else if(nPass[0]==1)
		{
			nStrand = 0;
			dScore = nMisD[0]+0.1*nMisC[0];
		}
		else if(nPass[1]==1)
		{
			nStrand = 1;
			dScore = nMisD[1]+0.1*nMisC[1];
		}
		else
		{
			nStrand = -1;
		}

		/* create new site if above threshold */
		if( (nStrand == 0) || (nStrand == 1) )
		{
			pSite = NULL;
			pSite = FLEXMTFSCREATE();
			if(pSite == NULL)
			{
				printf("Error: MotifMap_ScanConsensus_In_FlexSeqMtf, cannot create new site!\n");
				exit(EXIT_FAILURE);
			}

			pSite->nSeqId = pSeqMtf->nId;
			pSite->nStartPos = ni;
			pSite->nStrand = nStrand;
			pSite->pNext = NULL;
			pSite->pPrev = NULL;
			pSite->dScore = dScore;

			if(pPrev == NULL)
			{
				pSeqMtf->pMotifList = pSite;
				pPrev = pSite;
			}
			else
			{
				pPrev->pNext = pSite;
				pSite->pPrev = pPrev;
				pPrev = pSite;
			}

			pSeqMtf->nSiteNum += 1;
			*pSiteNum += 1;
		}



		if(nUseCS == 1)
		{
			dTotConserve -= pCS[ni];
		}
	}


	/* return */
	return PROC_SUCCESS;
}

		
/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanConsensus_MotifCountMisMatch: count mismatches to a motif */
/*  consensus.                                                             */
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanConsensus_MotifCountMisMatch(struct FLEXSEQMOTIF *pSeqMtf, 
					int nIndex, int nPos, char chStrand, int nMotifLen, 
					struct BYTEMATRIX *pCMotif, struct BYTEMATRIX *pDMotif, 
					int *pMisC, int *pMisD, int *pMask)
{

	/* define */
	int ni,nj;
	int nW;
	int nBaseTypeNum;
	int nSeqLen;
	unsigned char *pBase;
	
	/* init */
	*pMisC = 0;
	*pMisD = 0;
	*pMask = 0;
	nBaseTypeNum = pDMotif->nWidth;

	/* check */
	if( (pCMotif == NULL) || (pDMotif == NULL) )
	{
		printf("Error: MotifMap_ScanConsensus_MotifCountMisMatch, no motif pattern available!\n");
		exit(EXIT_FAILURE);
	}
	
	/* get likelihood */
	nMotifLen = pDMotif->nHeight;
	nSeqLen = pSeqMtf->vSeq[nIndex]->nWidth;
	if((nPos+nMotifLen) > nSeqLen)
	{
		*pMask = 1;
		return PROC_SUCCESS;
	}

	pBase = pSeqMtf->vSeq[nIndex]->pMatElement+nPos;
	
	if(chStrand == '+')
	{
		for(ni=0; ni<nMotifLen; ni++)
		{
			nW = (int)(pBase[ni]);
			
			if(nW >= nBaseTypeNum)
			{
				*pMask = 1;
				return PROC_SUCCESS;
			}

			/* if 'N' in consensus, tolerate */
			if((int)(pCMotif->pMatElement[ni]) == nBaseTypeNum)
			{
			}
			else
			{
				if(nW != (int)(pCMotif->pMatElement[ni]) )
				{
					*pMisC += 1;
				}
			}

			if(BMGETAT(pDMotif, ni, nW) != 1)
			{
				*pMisD += 1;
			}
		}
	}
	else if(chStrand == '-')
	{
		for(ni=0; ni<nMotifLen; ni++)
		{
			nW = (int)(pBase[ni]);
			
			if(nW >= nBaseTypeNum)
			{
				*pMask = 1;
				return PROC_SUCCESS;
			}

			nW = nBaseTypeNum-1-nW;
			nj = nMotifLen-1-ni;
			/* if 'N' in consensus, tolerate */
			if((int)(pCMotif->pMatElement[nj]) == nBaseTypeNum)
			{
			}
			else
			{
				if( nW != (int)(pCMotif->pMatElement[nj]) )
				{
					*pMisC += 1;
				}
			}

			if(BMGETAT(pDMotif, nj, nW) != 1)
			{
				*pMisD += 1;
			}
		}
	}
	else
	{
		printf("Error: MotifMap_ScanConsensus_MotifCountMisMatch, strand information wrong!\n");
		exit(EXIT_FAILURE);
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanConsensus_Export_Single: save motif sites to a file.      */
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanConsensus_Export_Single(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,
						int nMotifLen, char strSeqAlias[], FILE *fpOut)
{
	/* define */
	int nj;
	struct FLEXMOTIFSITE *pSite;
	char chStrand;
	char chBase;
	unsigned char *pBase;

	/* export */
	pSite = pSeqMtf->pMotifList;
	while(pSite != NULL)
	{
		if(pSite->nStrand == 0)
			chStrand = '+';
		else
			chStrand = '-';

		/* fprintf(fpOut, ">%d:%d-%d(%c)\n", pSeqMtf->nId, pSite->nStartPos, pSite->nStartPos+nMotifLen-1,
			chStrand); */

		fprintf(fpOut, "%s\t%d\t%d\t%d\t%c\t%f\t", strSeqAlias, pSeqMtf->nId, pSite->nStartPos, pSite->nStartPos+nMotifLen-1,
			chStrand, -(pSite->dScore));

		if(pSite->nStrand == 0)
		{
			pBase = pSeqMtf->vSeq[nIndex]->pMatElement+pSite->nStartPos;
			for(nj=0; nj<nMotifLen; nj++)
			{
				switch(*pBase)
				{
					case 0: chBase = 'A';
						break;
					case 1: chBase = 'C';
						break;
					case 2: chBase = 'G';
						break;
					case 3: chBase = 'T';
						break;
					case 10: chBase = 'a';
						break;
					case 11: chBase = 'c';
						break;
					case 12: chBase = 'g';
						break;
					case 13: chBase = 't';
						break;
					default: chBase = 'N';
				}

				fprintf(fpOut, "%c", chBase);
				pBase++;
			}
			fprintf(fpOut, "\n");
		}
		else
		{
			pBase = pSeqMtf->vSeq[nIndex]->pMatElement+pSite->nStartPos+nMotifLen-1;
			for(nj=0; nj<nMotifLen; nj++)
			{
				switch(*pBase)
				{
					case 0: chBase = 'T';
						break;
					case 1: chBase = 'G';
						break;
					case 2: chBase = 'C';
						break;
					case 3: chBase = 'A';
						break;
					case 10: chBase = 't';
						break;
					case 11: chBase = 'g';
						break;
					case 12: chBase = 'c';
						break;
					case 13: chBase = 'a';
						break;
					default: chBase = 'N';
				}

				fprintf(fpOut, "%c", chBase);
				pBase--;
			}
			fprintf(fpOut, "\n");
		}

		pSite = pSite->pNext;
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_GenomeBackground_Main: compute local markov        */
/*  background matrices. These matrices will be used by genome-wide        */
/*  matrix scan routines.                                                  */ 
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanMatrix_GenomeBackground_Main(char strGenomePath[], 
					char strOutPath[], int nBGOrder, int nS, int nW)
{
	/* define */
	char strFileName[MED_LINE_LENGTH];
	char strExportPath[MED_LINE_LENGTH];
	struct INTMATRIX *pChrSize;
	char strCommand[MED_LINE_LENGTH];
	char strLine[LINE_LENGTH];
	int ni;
	FILE *fpIn;
	int nChrLen;
	int nlen;

	/* ####################### */
	/* load initial parameters */
	/* ####################### */
	if(nW%nS != 0)
	{
		printf("Error: MotifMap_ScanMatrix_GenomeBackground_Main, W (window size) must be a multiple of S (step size)!\n");
		exit(EXIT_FAILURE);
	}

	if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
	{
		nlen = strlen(strGenomePath);
		if(nlen == 0)
		{
			sprintf(strGenomePath, ".\\"); 
		}
		else if(strGenomePath[nlen-1] != '\\')
		{
			strGenomePath[nlen] = '\\';
			strGenomePath[nlen+1] = '\0';
		}

		nlen = strlen(strOutPath);
		if(nlen == 0)
		{
			sprintf(strOutPath, ".\\"); 
		}
		else if(strOutPath[nlen-1] != '\\')
		{
			strOutPath[nlen] = '\\';
			strOutPath[nlen+1] = '\0';
		}
	}
	else
	{
		nlen = strlen(strGenomePath);
		if(nlen == 0)
		{
			sprintf(strGenomePath, "./"); 
		}
		else if(strGenomePath[nlen-1] != '/')
		{
			strGenomePath[nlen] = '/';
			strGenomePath[nlen+1] = '\0';
		}

		nlen = strlen(strOutPath);
		if(nlen == 0)
		{
			sprintf(strOutPath, "./"); 
		}
		else if(strOutPath[nlen-1] != '/')
		{
			strOutPath[nlen] = '/';
			strOutPath[nlen+1] = '\0';
		}
	}

	/* ###################### */
	/* Load chromosome size   */
	/* ###################### */
	sprintf(strFileName, "%schrlen.txt", strGenomePath);
	pChrSize = NULL;
	pChrSize = IMLOAD(strFileName);
	if(pChrSize == NULL)
	{
		printf("Warning: MotifMap_ScanMatrix_GenomeBackground_Main, no chromosomes are processed since chromosome length were not loaded!\n");
		return PROC_SUCCESS;
	}

	/* ###################### */
	/* process chromosomes    */
	/* ###################### */
	sprintf(strFileName, "%schrlist.txt", strGenomePath);
	fpIn = NULL;
	fpIn = fopen(strFileName, "r");
	if(fpIn == NULL)
	{
		printf("Error: MotifMap_ScanMatrix_GenomeBackground_Main, cannot open chrlist.txt file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
		{
			sprintf(strCommand, "md %s%s", strOutPath, strLine);
			system(strCommand);
			sprintf(strExportPath, "%s%s\\", strOutPath, strLine);
		}
		else
		{
			sprintf(strCommand, "mkdir %s%s", strOutPath, strLine);
			system(strCommand);
			sprintf(strExportPath, "%s%s/", strOutPath, strLine);
		}

		nChrLen = pChrSize->pMatElement[ni];
				
		/* compute background */
		MotifMap_ScanMatrix_GenomeBackground_Chr(strGenomePath, strLine, 
			strExportPath, nChrLen, nS, nW, nBGOrder); 

		ni++;
	}

	fclose(fpIn);
	if(ni != pChrSize->nHeight)
	{
		printf("Error: MotifMap_ScanMatrix_GenomeBackground_Main, chromosome number correct!\n");
		exit(EXIT_FAILURE);
	}

	/* ###################### */
	/* Release memory         */
	/* ###################### */
	DestroyIntMatrix(pChrSize);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_GenomeBackground_Chr: compute local markov         */
/*  background matrices for a single chromosome. These matrices will be    */
/*  used by genome-wide matrix scan routines.                              */ 
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanMatrix_GenomeBackground_Chr(char strGenomePath[], char strChr[],
			char strExportPath[], int nChrLen, int nS, int nW, int nBGOrder)
{
	/* define */
	struct DOUBLEMATRIX **vBGF;
	struct DOUBLEMATRIX **vBGB;
	struct DOUBLEMATRIX *pBGWF;
	struct DOUBLEMATRIX *pBGWB;
	struct DOUBLEMATRIX *pBGNF;
	struct DOUBLEMATRIX *pBGNB;

	/* all sequences */
	struct FLEXSEQMOTIF *pSeqMtf = NULL;
	int nSeqCount = 0;

	/* background order */
	int nBaseTypeNum = 4;
	int nBGHeight;

	/* others */
	int nSN,nWN,nHN;
	int nM1,nM2;
	int ni,nj,nk;
	char strFileName[MED_LINE_LENGTH];

	/* initialize */
	nBGHeight = (int)(pow((double)nBaseTypeNum, (double)nBGOrder));
	nSN = nChrLen/nS;
	if(nChrLen%nS != 0)
	{
		nSN += 1;
	}
	nWN = nW/nS;
	if(nWN%2 == 0)
	{
		nHN = nWN/2;
		nWN += 1;
	}
	else
	{
		nHN = (nWN-1)/2;
	}

	vBGF = NULL;
	vBGF = (struct DOUBLEMATRIX **)calloc(nSN, sizeof(struct DOUBLEMATRIX *));
	if(vBGF == NULL)
	{
		printf("Error: MotifMap_ScanMatrix_GenomeBackground_Chr, cannot allocate memory for local counting matrix!\n");
		exit(EXIT_FAILURE);
	}
	vBGB = NULL;
	vBGB = (struct DOUBLEMATRIX **)calloc(nSN, sizeof(struct DOUBLEMATRIX *));
	if(vBGB == NULL)
	{
		printf("Error: MotifMap_ScanMatrix_GenomeBackground_Chr, cannot allocate memory for local counting matrix!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nSN; ni++)
	{
		vBGF[ni] = CreateDoubleMatrix(nBGHeight, nBaseTypeNum);
		if(vBGF[ni] == NULL)
		{
			printf("Error: MotifMap_ScanMatrix_GenomeBackground_Chr, cannot allocate memory for local counting matrix!\n");
			exit(EXIT_FAILURE);
		}
		vBGB[ni] = CreateDoubleMatrix(nBGHeight, nBaseTypeNum);
		if(vBGB[ni] == NULL)
		{
			printf("Error: MotifMap_ScanMatrix_GenomeBackground_Chr, cannot allocate memory for local counting matrix!\n");
			exit(EXIT_FAILURE);
		}

		for(nj=0; nj<nBGHeight; nj++)
		{
			for(nk=0; nk<nBaseTypeNum; nk++)
			{
				DMSETAT(vBGF[ni], nj, nk, 1.0);
				DMSETAT(vBGB[ni], nj, nk, 1.0);
			}
		}
	}

	/* compute matrix for each step */
	nM1 = 0;
	ni = 0;
	while(nM1<nChrLen)
	{
		/* get coordinates */
		nM2 = nM1+nS-1;
		if(nM2 >= nChrLen)
			nM2 = nChrLen-1;

		/* load sequences */
		pSeqMtf = NULL;
		pSeqMtf = FLEXSEQMTFLOADSEQ_FOR_COUNTMCBACKGROUND_GENOME(ni, strGenomePath, strChr,
			nM1, nM2, 0, "");

		/* count background */
		MotifMap_ScanMatrix_Genome_CountBackground(pSeqMtf, 0, 0, pSeqMtf->nSeqLen-1, nBGOrder,
			vBGF[ni], vBGB[ni]);

		/* get next */
		FLEXSEQMTFDESTROY(pSeqMtf);
		nM1 = nM2+1;
		ni++;
	}
	
	if(ni != nSN)
	{
		printf("Error: MotifMap_ScanMatrix_GenomeBackground_Chr, chromosome length not match!\n");
		exit(EXIT_FAILURE);
	}

	/* compute matrix for each window */
	pBGWF = NULL;
	pBGWF = CreateDoubleMatrix(nBGHeight, nBaseTypeNum);
	if(pBGWF == NULL)
	{
		printf("Error: MotifMap_ScanMatrix_GenomeBackground_Chr, cannot create window count matrix!\n");
		exit(EXIT_FAILURE);
	}
	pBGWB = NULL;
	pBGWB = CreateDoubleMatrix(nBGHeight, nBaseTypeNum);
	if(pBGWB == NULL)
	{
		printf("Error: MotifMap_ScanMatrix_GenomeBackground_Chr, cannot create window count matrix!\n");
		exit(EXIT_FAILURE);
	}

	pBGNF = NULL;
	pBGNF = CreateDoubleMatrix(nBGHeight, nBaseTypeNum);
	if(pBGNF == NULL)
	{
		printf("Error: MotifMap_ScanMatrix_GenomeBackground_Chr, cannot create window count matrix!\n");
		exit(EXIT_FAILURE);
	}
	pBGNB = NULL;
	pBGNB = CreateDoubleMatrix(nBGHeight, nBaseTypeNum);
	if(pBGNB == NULL)
	{
		printf("Error: MotifMap_ScanMatrix_GenomeBackground_Chr, cannot create window count matrix!\n");
		exit(EXIT_FAILURE);
	}

	if(nBGOrder > 0)
	{
		for(ni=0; ni<nWN; ni++)
		{
			if(ni == nSN)
				break;

			DMADDTF(pBGWF, vBGF[ni]);
			DMADDTF(pBGWB, vBGB[ni]);
		}
		DMCOPY(pBGNF, pBGWF);
		DMCOPY(pBGNB, pBGWB);
		DMROWNORM(pBGNF);
		DMROWNORM(pBGNB);

		for(nj=0; nj<=nHN; nj++)
		{
			if(nj == nSN)
				break;

			sprintf(strFileName, "%s%d_f.txt", strExportPath, nj);
			DMSAVE(pBGNF, strFileName);
			sprintf(strFileName, "%s%d_b.txt", strExportPath, nj);
			DMSAVE(pBGNB, strFileName);
		}

		for(; ni<nSN; ni++)
		{
			DMSUBTF(pBGWF, vBGF[ni-nWN]);
			DMSUBTF(pBGWB, vBGB[ni-nWN]);
			DMADDTF(pBGWF, vBGF[ni]);
			DMADDTF(pBGWB, vBGB[ni]);

			DMCOPY(pBGNF, pBGWF);
			DMCOPY(pBGNB, pBGWB);
			DMROWNORM(pBGNF);
			DMROWNORM(pBGNB);

			sprintf(strFileName, "%s%d_f.txt", strExportPath, nj);
			DMSAVE(pBGNF, strFileName);
			sprintf(strFileName, "%s%d_b.txt", strExportPath, nj);
			DMSAVE(pBGNB, strFileName);
			nj++;
		}

		for(; nj<nSN; nj++)
		{
			sprintf(strFileName, "%s%d_f.txt", strExportPath, nj);
			DMSAVE(pBGNF, strFileName);
			sprintf(strFileName, "%s%d_b.txt", strExportPath, nj);
			DMSAVE(pBGNB, strFileName);
		}
	}
	else
	{
		for(ni=0; ni<nWN; ni++)
		{
			if(ni == nSN)
				break;

			DMADDTF(pBGWF, vBGF[ni]);
		}
		DMCOPY(pBGNF, pBGWF);
		DMROWNORM(pBGNF);
		
		for(nj=0; nj<=nHN; nj++)
		{
			if(nj == nSN)
				break;

			sprintf(strFileName, "%s%d.txt", strExportPath, nj);
			DMSAVE(pBGNF, strFileName);
		}

		for(; ni<nSN; ni++)
		{
			DMSUBTF(pBGWF, vBGF[ni-nWN]);
			DMADDTF(pBGWF, vBGF[ni]);
		
			DMCOPY(pBGNF, pBGWF);
			DMROWNORM(pBGNF);
		
			sprintf(strFileName, "%s%d.txt", strExportPath, nj);
			DMSAVE(pBGNF, strFileName);
			nj++;
		}

		for(; nj<nSN; nj++)
		{
			sprintf(strFileName, "%s%d.txt", strExportPath, nj);
			DMSAVE(pBGNF, strFileName);
		}
	}


	/* release memory */
	for(ni=0; ni<nSN; ni++)
	{
		DestroyDoubleMatrix(vBGF[ni]);
		vBGF[ni] = NULL;
		DestroyDoubleMatrix(vBGB[ni]);
		vBGB[ni] = NULL;
	}
	free(vBGF);
	free(vBGB);

	DestroyDoubleMatrix(pBGWF);
	DestroyDoubleMatrix(pBGWB);
	DestroyDoubleMatrix(pBGNF);
	DestroyDoubleMatrix(pBGNB);

	/* return */
	return PROC_SUCCESS;
}

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
					int nIncludeRepeat)
{
	/* background type: 0=region; 1=genome */
	int nBGType = 0;

	/* sequences and conservation */
	/* all sequences */
	struct FLEXSEQMOTIF *pSeqMtf = NULL;
	/* sequence number */
	int nSeqCount = 0;
	/* site number */
	int nEffecLen = 0;
	int nConsEffecLen = 0;
	int nSiteNum = 0;
	/* total sequence length and total site number */
	double dTotLen = 0.0;
	double dTotConsLen = 0.0;
	double dTotSite = 0.0;
	
	/* background */
	struct DOUBLEMATRIX *pBGF = NULL;
	struct DOUBLEMATRIX *pBGB = NULL;
	/* log background */
	struct DOUBLEMATRIX *pLogBGF = NULL;
	struct DOUBLEMATRIX *pLogBGB = NULL;
	/* background0 */
	struct DOUBLEMATRIX *pBG0 = NULL;
	/* log background0 */
	struct DOUBLEMATRIX *pLogBG0 = NULL;
	
	/* motif matrix */
	struct DOUBLEMATRIX *pMotif = NULL;
	
	/* number of base types */
	int nBaseTypeNum = 4;
	int nScale = 0;
	int nBGHeight;

	/* motif consensus */
	FILE *fpIn;
	FILE *fpOut;
	FILE *fpOut2;
	char strLine[LONG_LINE_LENGTH];
	char strSeqAlias[LINE_LENGTH];

	/* coordinates */
	int nMotifLen;
	char strChr[LINE_LENGTH];
	int nStart,nEnd;
	int nM1,nM2;
	int nActualStart,nActualEnd;
	int nFrom,nTo;
	
	/* other variables */
	int ni,nj;
	double *pEle1,*pEle2,*pElement;
	double dSum;
	int nlen;
	
	/* #################################### */
	/* initialize                           */
	/* #################################### */
	nBaseTypeNum = 4;
	nScale = (int)(pow((double)nBaseTypeNum, (double)(nBGOrder-1)));
	nBGHeight = (int)(pow((double)nBaseTypeNum, (double)nBGOrder));
	
	if(strcmp(strBGType, "GENOME") == 0)
	{
		nBGType = 1;
	}
	else
	{
		nBGType = 0;
	}
	
	/* #################################### */
	/* load initial parameters              */
	/* #################################### */

	if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
	{
		nlen = strlen(strGenomePath);
		if(nlen == 0)
		{
			sprintf(strGenomePath, ".\\"); 
		}
		else if(strGenomePath[nlen-1] != '\\')
		{
			strGenomePath[nlen] = '\\';
			strGenomePath[nlen+1] = '\0';
		}

		if(nUseCS == 1)
		{
			nlen = strlen(strCSPath);
			if(nlen == 0)
			{
				sprintf(strCSPath, ".\\"); 
			}
			else if(strCSPath[nlen-1] != '\\')
			{
				strCSPath[nlen] = '\\';
				strCSPath[nlen+1] = '\0';
			}
		}
		if(nBGType == 1)
		{
			nlen = strlen(strBGPath);
			if(nlen == 0)
			{
				sprintf(strBGPath, ".\\"); 
			}
			else if(strBGPath[nlen-1] != '\\')
			{
				strBGPath[nlen] = '\\';
				strBGPath[nlen+1] = '\0';
			}
		}
	}
	else
	{
		nlen = strlen(strGenomePath);
		if(nlen == 0)
		{
			sprintf(strGenomePath, "./"); 
		}
		else if(strGenomePath[nlen-1] != '/')
		{
			strGenomePath[nlen] = '/';
			strGenomePath[nlen+1] = '\0';
		}

		if(nUseCS == 1)
		{
			nlen = strlen(strCSPath);
			if(nlen == 0)
			{
				sprintf(strCSPath, "./"); 
			}
			else if(strCSPath[nlen-1] != '/')
			{
				strCSPath[nlen] = '/';
				strCSPath[nlen+1] = '\0';
			}
		}

		if(nBGType == 1)
		{
			nlen = strlen(strBGPath);
			if(nlen == 0)
			{
				sprintf(strBGPath, "./"); 
			}
			else if(strBGPath[nlen-1] != '/')
			{
				strBGPath[nlen] = '/';
				strBGPath[nlen+1] = '\0';
			}
		}
	}

	/* #################################### */
	/* load motif                           */
	/* #################################### */
	pMotif = NULL;
	pMotif = DMLOAD(strMotifPath);
	if( pMotif == NULL)
	{
		printf("Error: MotifMap_ScanMatrix_Genome_Main, cannot load init motif pseudocount!\n");
		exit(EXIT_FAILURE);
	}

	nMotifLen = pMotif->nHeight;

	pEle1 = pMotif->pMatElement;
	for(ni=0; ni<pMotif->nHeight; ni++)
	{
		pEle2 = pEle1;
		dSum = 0.0;
		for(nj=0; nj<pMotif->nWidth; nj++)
		{
			*pEle2 += 1e-3;
			dSum += (*pEle2);
			pEle2++;
		}

		for(nj=0; nj<pMotif->nWidth; nj++)
		{
			*pEle1 = log(*pEle1/dSum);
			pEle1++;
		}
	}


	/* #################################### */
	/* compute background if needed         */
	/* #################################### */
	
	/* genomic local background */
	if(nBGType == 1)
	{
	}
	/* region based background */
	else
	{
		pBGF = NULL;
		pBGF = CreateDoubleMatrix(nBGHeight, nBaseTypeNum);
		pBGB = NULL;
		pBGB = CreateDoubleMatrix(nBGHeight, nBaseTypeNum);
		pBG0 = NULL;
		pBG0 = CreateDoubleMatrix(1, nBaseTypeNum);
		for(ni=0; ni<nBGHeight; ni++)
		{
			for(nj=0; nj<nBaseTypeNum; nj++)
			{
				DMSETAT(pBGF, ni, nj, 1.0);
				DMSETAT(pBGB, ni, nj, 1.0);
			}
		}
		for(nj=0; nj<nBaseTypeNum; nj++)
		{
			DMSETAT(pBG0, 0, nj, 1.0);
		}
		
		if( (pBGF == NULL) || (pBGB == NULL) || (pBG0 == NULL) )
		{
			printf("Error: MotifMap_ScanMatrix_Genome_Main, failure to create background markov matrix!\n");
			exit(EXIT_FAILURE);
		}

		/* compute background */
		MotifMap_ScanMatrix_Genome_ComputeBackground_Region(strGenomePath,
					strCodPath, nBGOrder, pBGF, pBGB, pBG0);
		
		pLogBGF = NULL;
		pLogBGF = DMCLONE(pBGF);
		pElement = pLogBGF->pMatElement;
		for(ni=0; ni<pLogBGF->nHeight; ni++)
		{
			for(nj=0; nj<pLogBGF->nWidth; nj++)
			{
				*pElement = log(*pElement);
				pElement++;
			}
		}

		pLogBGB = NULL;
		pLogBGB = DMCLONE(pBGB);
		pElement = pLogBGB->pMatElement;
		for(ni=0; ni<pLogBGB->nHeight; ni++)
		{
			for(nj=0; nj<pLogBGB->nWidth; nj++)
			{
				*pElement = log(*pElement);
				pElement++;
			}
		}

		pLogBG0 = NULL;
		pLogBG0 = DMCLONE(pBG0);
		pElement = pLogBG0->pMatElement;
		for(ni=0; ni<pLogBG0->nHeight; ni++)
		{
			for(nj=0; nj<pLogBG0->nWidth; nj++)
			{
				*pElement = log(*pElement);
				pElement++;
			}
		}

		sprintf(strLine, "%s.bg0", strOutputPath);
		DMSAVE(pBG0, strLine);
		sprintf(strLine, "%s.bgf", strOutputPath);
		DMSAVE(pBGF, strLine);
		sprintf(strLine, "%s.bgb", strOutputPath);
		DMSAVE(pBGB, strLine);
	}

	/* #################################### */
	/* load sequences and conservation      */
	/* #################################### */
	fpIn = NULL;
	fpIn = fopen(strCodPath, "r");
	if(fpIn == NULL)
	{
        printf("Error: MotifMap_ScanMatrix_Genome_Main, cannot open coordinates file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
        printf("Error: MotifMap_ScanMatrix_Genome_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	fprintf(fpOut, "#seq_id\tchromosome\tstart\tend\tstrand\tlog10(Likelihood_Ratio)\tsite\n");
	
	sprintf(strLine, "%s.stat", strOutputPath);
	fpOut2 = NULL;
	fpOut2 = fopen(strLine, "w");
	if(fpOut2 == NULL)
	{
        printf("Error: MotifMap_ScanMatrix_Genome_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	/* read */
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		sscanf(strLine, "%s %s %d %d", strSeqAlias, strChr,
			&nStart, &nEnd);
		
		/* process segment by segment */
		nM1 = nStart;
		while(nM1 <= nEnd)
		{
			nM2 = nM1+GENOME_CONTIG_LEN-1;
			if(nM2 > nEnd)
				nM2 = nEnd;

			nActualStart = nM1-nBGOrder;
			if(nActualStart >= nStart)
			{
				nFrom = nBGOrder;
			}
			else
			{
				nActualStart = nStart;
				nFrom = 0;
			}

			nActualEnd = nM2+nMotifLen-1;
			if(nActualEnd > nEnd)
			{
				nActualEnd = nEnd;
				nTo = nActualEnd-nActualStart+1-nMotifLen;
			}
			else
			{
				nTo = nM2-nM1+nFrom;
				nActualEnd += nBGOrder;
				if(nActualEnd > nEnd)
				{
					nActualEnd = nEnd;
				}
			}

			if(nTo < nFrom)
			{
				break;
			}

			/* #################################### */ 
			/* load sequence and score              */
			/* #################################### */
			pSeqMtf = NULL;
			pSeqMtf = FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_GENOME(nSeqCount,
				strGenomePath, strChr, nActualStart, nActualEnd, 
				nUseCS, strCSPath);
			if(pSeqMtf == NULL)
			{
				break;
			}
			if(nIncludeRepeat == 1)
			{
				FLEXSEQMTF_REPEATUNMASK(pSeqMtf, 0);
			}

			/* #################################### */
			/* set background                       */
			/* #################################### */
			/* genomic local background */
			if(nBGType == 1)
			{
				MotifMap_ScanMatrix_InitBGLogLike_Genome(pSeqMtf, 0,
					strChr, nActualStart, nActualEnd, 
					nBGOrder, strBGPath, nBGStepSize);
			}
			/* region based background */
			else
			{
				MotifMap_ScanMatrix_InitBGLogLike(pSeqMtf, 0,
					nBGOrder, pLogBGF, pLogBGB, pLogBG0);
			}
			
			/* #################################### */
			/* scan motif                           */
			/* #################################### */
			MotifMap_ScanMatrix_In_FlexSeqMtf(pSeqMtf, 0,
				nFrom, nTo, pMotif, dR, 
				nUseCS, dC, &nEffecLen, &nConsEffecLen, &nSiteNum);
			dTotLen += nEffecLen;
			dTotConsLen += nConsEffecLen;
			dTotSite += nSiteNum;

			/* #################################### */
			/* export results                       */
			/* #################################### */
			MotifMap_ScanMatrix_Export_Genome(pSeqMtf, 0,
				pMotif->nHeight, strSeqAlias, strChr, 
				nActualStart, fpOut);

			/* #################################### */
			/* release memory                       */
			/* #################################### */
			FLEXSEQMTFDESTROY(pSeqMtf);
			nM1 = nM2+1;
		}
		nSeqCount++;
	}


	/* save statistics */
	fprintf(fpOut2, "EffecLen= %d\n", (int)dTotLen);
	fprintf(fpOut2, "TotalSite= %d\n", (int)dTotSite);
	fprintf(fpOut2, "ConsLen= %d\n", (int)dTotConsLen);
	
	/* close files */
	fclose(fpIn);
	fclose(fpOut);
	fclose(fpOut2);

	/* #################################### */
	/* release memory                       */
	/* #################################### */
	DestroyDoubleMatrix(pMotif);
	if(nBGType == 0)
	{
		DestroyDoubleMatrix(pBGF);
		DestroyDoubleMatrix(pBGB);
		DestroyDoubleMatrix(pLogBGF);
		DestroyDoubleMatrix(pLogBGB);
		DestroyDoubleMatrix(pBG0);
		DestroyDoubleMatrix(pLogBG0);
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_Genome_ComputeBackground_Region: construct         */
/*  background markov matrices for specified regions.                      */ 
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanMatrix_Genome_ComputeBackground_Region(char strGenomePath[],
				char strCodPath[], int nBGOrder, 
				struct DOUBLEMATRIX *pBGF, struct DOUBLEMATRIX *pBGB, 
				struct DOUBLEMATRIX *pBG0)
{
	/* all sequences */
	struct FLEXSEQMOTIF *pSeqMtf = NULL;
	/* sequence number */
	int nSeqCount = 0;
	
	/* number of base types */
	int nBaseTypeNum = 4;
	double *pEle1,*pEle2;
	double dSum;

	/* motif consensus */
	FILE *fpIn;
	char strLine[LONG_LINE_LENGTH];
	char strSeqAlias[LINE_LENGTH];
	int ni,nj;
	
	/* other variables */
	char strChr[LINE_LENGTH];
	int nStart,nEnd;
	int nM1,nM2;
	int nActualStart,nActualEnd;
	int nFrom,nTo;


	/* #################################### */
	/* load sequences and conservation      */
	/* #################################### */
	fpIn = NULL;
	fpIn = fopen(strCodPath, "r");
	if(fpIn == NULL)
	{
        printf("Error: MotifMap_ScanMatrix_Genome_ComputeBackground_Region, cannot open sequence file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		sscanf(strLine, "%s %s %d %d", strSeqAlias, strChr,
			&nStart, &nEnd);
		
		/* process segment by segment */
		nM1 = nStart;
		while(nM1 <= nEnd)
		{
			nM2 = nM1+GENOME_CONTIG_LEN-1;
			if(nM2 > nEnd)
				nM2 = nEnd;

			nActualStart = nM1;
			nFrom = 0;
			nTo = nM2-nM1;

			nActualEnd = nM2+nBGOrder;
			if(nActualEnd > nEnd)
			{
				nActualEnd = nEnd;
			}
			
			if(nTo < nFrom)
			{
				break;
			}

			/* #################################### */ 
			/* load sequence and score              */
			/* #################################### */
			pSeqMtf = NULL;
			pSeqMtf = FLEXSEQMTFLOADSEQ_FOR_COUNTMCBACKGROUND_GENOME(nSeqCount,
				strGenomePath, strChr, nActualStart, nActualEnd,
				0, "");

			if(pSeqMtf == NULL)
			{
				break;
			}

			/* #################################### */
			/* compute background                   */
			/* #################################### */
			MotifMap_ScanMatrix_Genome_CountBackground(pSeqMtf, 0,
				nFrom, nTo, 0, pBG0, NULL);
			if(nBGOrder > 0)
			{
				MotifMap_ScanMatrix_Genome_CountBackground(pSeqMtf, 0,
					0, pSeqMtf->nSeqLen-1, nBGOrder, pBGF, pBGB);
			}

			/* #################################### */
			/* release memory                       */
			/* #################################### */
			FLEXSEQMTFDESTROY(pSeqMtf);
			nM1 = nM2+1;
		}
		nSeqCount++;
	}


	/* close files */
	fclose(fpIn);
	
	/* #################################### */
	/* normalization                        */
	/* #################################### */
	pEle1 = pBG0->pMatElement;
	for(ni=0; ni<pBG0->nHeight; ni++)
	{
		pEle2 = pEle1;
		dSum = 0.0;
		for(nj=0; nj<pBG0->nWidth; nj++)
		{
			dSum += (*pEle2);
			pEle2++;
		}

		for(nj=0; nj<pBG0->nWidth; nj++)
		{
			*pEle1 = *pEle1/dSum;
			pEle1++;
		}
	}

	
	if(nBGOrder > 0)
	{
		pEle1 = pBGF->pMatElement;
		for(ni=0; ni<pBGF->nHeight; ni++)
		{
			pEle2 = pEle1;
			dSum = 0.0;
			for(nj=0; nj<pBGF->nWidth; nj++)
			{
				dSum += (*pEle2);
				pEle2++;
			}

			for(nj=0; nj<pBGF->nWidth; nj++)
			{
				*pEle1 = *pEle1/dSum;
				pEle1++;
			}
		}

		pEle1 = pBGB->pMatElement;
		for(ni=0; ni<pBGB->nHeight; ni++)
		{
			pEle2 = pEle1;
			dSum = 0.0;
			for(nj=0; nj<pBGB->nWidth; nj++)
			{
				dSum += (*pEle2);
				pEle2++;
			}

			for(nj=0; nj<pBGB->nWidth; nj++)
			{
				*pEle1 = *pEle1/dSum;
				pEle1++;
			}
		}
	}
	else
	{
		DMCOPY(pBGF, pBG0);
		DMCOPY(pBGB, pBG0);
	}

	/* return */
	return nSeqCount;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_Genome_CountBackground: count empirical            */
/*  observations for constructing background markov matrices.              */ 
/* ----------------------------------------------------------------------- */
int MotifMap_ScanMatrix_Genome_CountBackground(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,
					int nFrom, int nTo, int nBGOrder, 
					struct DOUBLEMATRIX *pBGF, struct DOUBLEMATRIX *pBGB)
{
	/* define */
	int nBaseTypeNum,nScale;
	int nWordId,nWordBadLen;
	int nW,nWB;
	int nLen;
	unsigned char *pBase;
	/* unsigned char *pCS; */
	int ni,nj;
	double dTemp;

	/* init */
	nLen = pSeqMtf->vSeq[nIndex]->nWidth;
	pBase = pSeqMtf->vSeq[nIndex]->pMatElement;
	nBaseTypeNum = pBGF->nWidth;
	if( nBGOrder == 0)
	{
		nScale = 0;
	}
	else
	{
		nScale = (int)(pow((double)nBaseTypeNum, (double)(nBGOrder-1)));
	}
	
	/* if zero order Markov chain as background */
	if(nBGOrder == 0)
	{
		for(ni=nFrom; ni<=nTo; ni++)
		{
			nW = (int)pBase[ni];
			if(nW < nBaseTypeNum)
			{
				dTemp = DMGETAT(pBGF, 0, nW)+1.0;
				DMSETAT(pBGF, 0, nW, dTemp);
			}
		}
	}

	/* if higher order Markov chain as background */
	else if(nLen > (nBGOrder+nFrom))
	{
		nWordId = 0;
		nWordBadLen = 0;
		
		/* initial words */
		for(ni=nFrom; ni<(nFrom+nBGOrder); ni++)
		{
			nW = (int)pBase[ni];
			if(nW < nBaseTypeNum)
			{
				nWordId = nWordId*nBaseTypeNum+nW;
			}
			else
			{
				nWordId = nWordId*nBaseTypeNum;
				nWordBadLen++;
			}
		}

		/* continuing words */
		for(; ni<=nTo; ni++)
		{
			/* get forward bg */
			nW = (int)pBase[ni];
			if(nW < nBaseTypeNum)
			{
				if(nWordBadLen == 0)
				{
					dTemp = DMGETAT(pBGF, nWordId, nW)+1.0;
					DMSETAT(pBGF, nWordId, nW, dTemp);
				}
				else
				{		
				}
			}
			else
			{
			}
			
			/* update word id */
			/* minus old letter */
			nj = ni-nBGOrder;
			nWB = (int)(pBase[nj]);
			if(nWB < nBaseTypeNum)
			{
				nWordId -= nWB*nScale;
			}
			else
			{
				nWordBadLen--;
			}
						
			/* add new letter */
			nW = (int)pBase[ni];
			if(nW < nBaseTypeNum)
			{
				nWordId = nWordId*nBaseTypeNum+nW;
			}
			else
			{
				nWordId = nWordId*nBaseTypeNum;
				nWordBadLen++;
			}

			/* get backward likelihood */
			if(nWB < nBaseTypeNum)
			{
				if(nWordBadLen == 0)
				{	
					dTemp = DMGETAT(pBGB, nWordId, nWB)+1.0;
					DMSETAT(pBGB, nWordId, nWB, dTemp);
				}
				else
				{
				}
			}
			else
			{
			}
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_Export_Genome: save motif sites to a file.         */
/* ----------------------------------------------------------------------- */
int MotifMap_ScanMatrix_Export_Genome(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,
						int nMotifLen, char strSeqAlias[], 
						char strChr[], int nOffset, FILE *fpOut)
{
	/* define */
	int nj;
	struct FLEXMOTIFSITE *pSite;
	char chStrand;
	char chBase;
	unsigned char *pBase;

	/* export */
	pSite = pSeqMtf->pMotifList;
	while(pSite != NULL)
	{
		if(pSite->nStrand == 0)
			chStrand = '+';
		else
			chStrand = '-';

		/* fprintf(fpOut, ">%d:%d-%d(%c)\n", pSeqMtf->nId, pSite->nStartPos, pSite->nStartPos+nMotifLen-1,
			chStrand); */

		fprintf(fpOut, "%s\t%s\t%d\t%d\t%c\t%f\t", strSeqAlias, strChr, 
			nOffset+pSite->nStartPos, nOffset+pSite->nStartPos+nMotifLen-1,
			chStrand, pSite->dScore);

		if(pSite->nStrand == 0)
		{
			pBase = pSeqMtf->vSeq[nIndex]->pMatElement+pSite->nStartPos;
			for(nj=0; nj<nMotifLen; nj++)
			{
				switch(*pBase)
				{
					case 0: chBase = 'A';
						break;
					case 1: chBase = 'C';
						break;
					case 2: chBase = 'G';
						break;
					case 3: chBase = 'T';
						break;
					case 10: chBase = 'a';
						break;
					case 11: chBase = 'c';
						break;
					case 12: chBase = 'g';
						break;
					case 13: chBase = 't';
						break;
					default: chBase = 'N';
				}

				fprintf(fpOut, "%c", chBase);
				pBase++;
			}
			fprintf(fpOut, "\n");
		}
		else
		{
			pBase = pSeqMtf->vSeq[nIndex]->pMatElement+pSite->nStartPos+nMotifLen-1;
			for(nj=0; nj<nMotifLen; nj++)
			{
				switch(*pBase)
				{
					case 0: chBase = 'T';
						break;
					case 1: chBase = 'G';
						break;
					case 2: chBase = 'C';
						break;
					case 3: chBase = 'A';
						break;
					case 10: chBase = 't';
						break;
					case 11: chBase = 'g';
						break;
					case 12: chBase = 'c';
						break;
					case 13: chBase = 'a';
						break;
					default: chBase = 'N';
				}

				fprintf(fpOut, "%c", chBase);
				pBase--;
			}
			fprintf(fpOut, "\n");
		}

		pSite = pSite->pNext;
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_Sequential_Main: map motif probability matrix to   */
/*  FASTA sequences, using nBGOrder markov chain as background and dR as   */
/*  likelihood ratio cutoff. dC is the conservation score cutoff.          */ 
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanMatrix_Sequential_Main(char  strMotifPath[], char strInputPath[],
					char strSeqFile[], char strOutputPath[], 
					double dR, int nBGOrder,
					int nUseCS, double dC, char strCSAlias[])
{
	/* sequences and conservation */
	/* all sequences */
	struct FLEXSEQMOTIF *pSeqMtf = NULL;
	/* sequence number */
	int nSeqCount = 0;
	/* site number */
	int nEffecLen = 0;
	int nConsEffecLen = 0;
	int nSiteNum = 0;
	/* total sequence length and total site number */
	double dTotLen = 0.0;
	double dTotConsLen = 0.0;
	double dTotSite = 0.0;
	
	/* background */
	struct DOUBLEMATRIX *pBGF = NULL;
	struct DOUBLEMATRIX *pBGB = NULL;
	/* log background */
	struct DOUBLEMATRIX *pLogBGF = NULL;
	struct DOUBLEMATRIX *pLogBGB = NULL;
	/* background0 */
	struct DOUBLEMATRIX *pBG0 = NULL;
	/* log background0 */
	struct DOUBLEMATRIX *pLogBG0 = NULL;
	
	/* motif matrix */
	struct DOUBLEMATRIX *pMotif = NULL;
	
	/* number of base types */
	int nBaseTypeNum = 4;
	int nScale = 0;
	int nBGHeight;

	/* motif consensus */
	FILE *fpIn;
	FILE *fpOut;
	FILE *fpOut2;
	char strLine[MED_LINE_LENGTH];
	char strSeqLine[LONG_LINE_LENGTH];
	struct tagSequence *pNewSeq;
	char strSeqAlias[LINE_LENGTH];
	
	/* other variables */
	int ni,nj;
	double *pEle1,*pEle2,*pElement;
	double dSum;
	int nlen;
	
	/* #################################### */
	/* initialize                           */
	/* #################################### */
	nBaseTypeNum = 4;
	nScale = (int)(pow((double)nBaseTypeNum, (double)(nBGOrder-1)));
	nBGHeight = (int)(pow((double)nBaseTypeNum, (double)nBGOrder));
	
	/* #################################### */
	/* load initial parameters              */
	/* #################################### */
	if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
	{
		nlen = strlen(strInputPath);
		if(nlen == 0)
		{
			sprintf(strInputPath, ".\\"); 
		}
		else if(strInputPath[nlen-1] != '\\')
		{
			strInputPath[nlen] = '\\';
			strInputPath[nlen+1] = '\0';
		}
	}
	else
	{
		nlen = strlen(strInputPath);
		if(nlen == 0)
		{
			sprintf(strInputPath, "./"); 
		}
		else if(strInputPath[nlen-1] != '/')
		{
			strInputPath[nlen] = '/';
			strInputPath[nlen+1] = '\0';
		}
	}

	/* motif PWM */
	pMotif = NULL;
	pMotif = DMLOAD(strMotifPath);
	if( pMotif == NULL)
	{
		printf("Error: MotifMap_ScanMatrix_Sequential_Main, cannot load init motif pseudocount!\n");
		exit(EXIT_FAILURE);
	}
	pEle1 = pMotif->pMatElement;
	for(ni=0; ni<pMotif->nHeight; ni++)
	{
		pEle2 = pEle1;
		dSum = 0.0;
		for(nj=0; nj<pMotif->nWidth; nj++)
		{
			*pEle2 += 1e-3;
			dSum += (*pEle2);
			pEle2++;
		}

		for(nj=0; nj<pMotif->nWidth; nj++)
		{
			*pEle1 = log(*pEle1/dSum);
			pEle1++;
		}
	}

	/* #################################### */
	/* compute background markov chain      */
	/* #################################### */
	pBGF = NULL;
	pBGF = CreateDoubleMatrix(nBGHeight, nBaseTypeNum);
	pBGB = NULL;
	pBGB = CreateDoubleMatrix(nBGHeight, nBaseTypeNum);
	pBG0 = NULL;
	pBG0 = CreateDoubleMatrix(1, nBaseTypeNum);
	for(ni=0; ni<nBGHeight; ni++)
	{
		for(nj=0; nj<nBaseTypeNum; nj++)
		{
			DMSETAT(pBGF, ni, nj, 1.0);
			DMSETAT(pBGB, ni, nj, 1.0);
		}
	}
	for(nj=0; nj<nBaseTypeNum; nj++)
	{
		DMSETAT(pBG0, 0, nj, 1.0);
	}
	
	if( (pBGF == NULL) || (pBGB == NULL) || (pBG0 == NULL) )
	{
		printf("Error: MotifMap_ScanMatrix_Sequential_Main, failure to create background markov matrix!\n");
		exit(EXIT_FAILURE);
	}

	/* compute background */
	MotifMap_ScanMatrix_Sequential_ComputeBackground(strInputPath,
				strSeqFile, nBGOrder, pBGF, pBGB, pBG0);
	
	pLogBGF = NULL;
	pLogBGF = DMCLONE(pBGF);
	pElement = pLogBGF->pMatElement;
	for(ni=0; ni<pLogBGF->nHeight; ni++)
	{
		for(nj=0; nj<pLogBGF->nWidth; nj++)
		{
			*pElement = log(*pElement);
			pElement++;
		}
	}

	pLogBGB = NULL;
	pLogBGB = DMCLONE(pBGB);
	pElement = pLogBGB->pMatElement;
	for(ni=0; ni<pLogBGB->nHeight; ni++)
	{
		for(nj=0; nj<pLogBGB->nWidth; nj++)
		{
			*pElement = log(*pElement);
			pElement++;
		}
	}

	pLogBG0 = NULL;
	pLogBG0 = DMCLONE(pBG0);
	pElement = pLogBG0->pMatElement;
	for(ni=0; ni<pLogBG0->nHeight; ni++)
	{
		for(nj=0; nj<pLogBG0->nWidth; nj++)
		{
			*pElement = log(*pElement);
			pElement++;
		}
	}

	sprintf(strLine, "%s.bg0", strOutputPath);
	DMSAVE(pBG0, strLine);
	sprintf(strLine, "%s.bgf", strOutputPath);
	DMSAVE(pBGF, strLine);
	sprintf(strLine, "%s.bgb", strOutputPath);
	DMSAVE(pBGB, strLine);

	/* #################################### */
	/* load sequences and conservation      */
	/* #################################### */
	sprintf(strLine, "%s%s", strInputPath, strSeqFile);
	fpIn = NULL;
	fpIn = fopen(strLine, "r");
	if(fpIn == NULL)
	{
        printf("Error: MotifMap_ScanMatrix_Sequential_Main, cannot open sequence file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
        printf("Error: MotifMap_ScanMatrix_Sequential_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	fprintf(fpOut, "#seq_name\tseq_id\tstart\tend\tstrand\tlog10(Likelihood_Ratio)\tsite\n");
	
	sprintf(strLine, "%s.stat", strOutputPath);
	fpOut2 = NULL;
	fpOut2 = fopen(strLine, "w");
	if(fpOut2 == NULL)
	{
        printf("Error: MotifMap_ScanMatrix_Sequential_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	/* read */
	pNewSeq = NULL;
	while(fgets(strSeqLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strSeqLine);
		StrTrimRight(strSeqLine);
		if(strSeqLine[0] == '\0')
			continue;

		if(strSeqLine[0] == '>')
		{
			if(pNewSeq != NULL)
			{
				/* #################################### */
				/* load sequence and score              */
				/* #################################### */
				pSeqMtf = NULL;
				pSeqMtf = FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_SINGLE(pNewSeq,
					strInputPath, nUseCS, strCSAlias);
				SequenceDelete(pNewSeq);

				/* #################################### */
				/* init background loglikelihood        */
				/* #################################### */
				MotifMap_ScanMatrix_InitBGLogLike(pSeqMtf, 0,
					nBGOrder, pLogBGF, pLogBGB, pLogBG0); 
				
				/* #################################### */
				/* scan motif                           */
				/* #################################### */
				MotifMap_ScanMatrix_In_FlexSeqMtf(pSeqMtf, 0,
					0, pSeqMtf->nSeqLen-1, pMotif, dR, 
					nUseCS, dC, &nEffecLen, &nConsEffecLen, &nSiteNum);
				dTotLen += nEffecLen;
				dTotConsLen += nConsEffecLen;
				dTotSite += nSiteNum;

				/* #################################### */
				/* export results                       */
				/* #################################### */
				MotifMap_ScanMatrix_Export_Single(pSeqMtf, 0,
					pMotif->nHeight, strSeqAlias, fpOut);

				/* #################################### */
				/* release memory                       */
				/* #################################### */
				FLEXSEQMTFDESTROY(pSeqMtf);
				nSeqCount++;
			}

			pNewSeq = NULL;
			pNewSeq = SequenceCreate();
			if(pNewSeq == NULL)
			{
				fclose(fpIn);
				printf("Error: MotifMap_ScanMatrix_Sequential_Main, cannot allocate memory for loading sequences!\n");
				exit(EXIT_FAILURE);
			}
			pNewSeq->m_nIndex = nSeqCount;
			strcpy(pNewSeq->m_strAlias, (strSeqLine+1));
			strcpy(strSeqAlias, (strSeqLine+1));
		}
		else
		{
			if(SequenceAddTail(pNewSeq, strSeqLine) == PROC_FAILURE)
			{
				fclose(fpIn);
				printf("Error: MotifMap_ScanMatrix_Sequential_Main, cannot allocate memory for loading sequences!\n");
				exit(EXIT_FAILURE);
			}			
		}
	}

	/* the last sequence */
	if(pNewSeq != NULL)
	{
		/* #################################### */
		/* load sequence and score              */
		/* #################################### */
		pSeqMtf = NULL;
		pSeqMtf = FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_SINGLE(pNewSeq,
			strInputPath, nUseCS, strCSAlias);
		SequenceDelete(pNewSeq);

		/* #################################### */
		/* init background loglikelihood        */
		/* #################################### */
		MotifMap_ScanMatrix_InitBGLogLike(pSeqMtf, 0,
			nBGOrder, pLogBGF, pLogBGB, pLogBG0); 
		
		/* #################################### */
		/* scan motif                           */
		/* #################################### */
		MotifMap_ScanMatrix_In_FlexSeqMtf(pSeqMtf, 0,
			0, pSeqMtf->nSeqLen-1, pMotif, dR, 
			nUseCS, dC, &nEffecLen, &nConsEffecLen, &nSiteNum);
		dTotLen += nEffecLen;
		dTotConsLen += nConsEffecLen;
		dTotSite += nSiteNum;

		/* #################################### */
		/* export results                       */
		/* #################################### */
		MotifMap_ScanMatrix_Export_Single(pSeqMtf, 0,
			pMotif->nHeight, strSeqAlias, fpOut);

		/* #################################### */
		/* release memory                       */
		/* #################################### */
		FLEXSEQMTFDESTROY(pSeqMtf);
		nSeqCount++;
	}

	/* save statistics */
	fprintf(fpOut2, "EffecLen= %d\n", (int)dTotLen);
	fprintf(fpOut2, "TotalSite= %d\n", (int)dTotSite);
	fprintf(fpOut2, "ConsLen= %d\n", (int)dTotConsLen);
		
	/* close files */
	fclose(fpIn);
	fclose(fpOut);
	fclose(fpOut2);

	/* #################################### */
	/* release memory                       */
	/* #################################### */
	
	DestroyDoubleMatrix(pMotif);
	DestroyDoubleMatrix(pBGF);
	DestroyDoubleMatrix(pBGB);
	DestroyDoubleMatrix(pLogBGF);
	DestroyDoubleMatrix(pLogBGB);
	DestroyDoubleMatrix(pBG0);
	DestroyDoubleMatrix(pLogBG0);
	

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_Sequential_ComputeBackground: construct background */
/*  markov matrices.                                                       */ 
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanMatrix_Sequential_ComputeBackground(char strInputPath[],
				char strSeqFile[], int nBGOrder, 
				struct DOUBLEMATRIX *pBGF, struct DOUBLEMATRIX *pBGB, 
				struct DOUBLEMATRIX *pBG0)
{
	/* all sequences */
	struct FLEXSEQMOTIF *pSeqMtf = NULL;
	/* sequence number */
	int nSeqCount = 0;
	/* site number */
	int nEffecLen = 0;
	int nSiteNum = 0;
	/* total sequence length and total site number */
	double dTotLen = 0.0;
	double dTotSite = 0.0;
	
	/* number of base types */
	int nBaseTypeNum = 4;
	int nScale = 0;
	double *pEle1,*pEle2;
	double dSum;

	/* motif consensus */
	FILE *fpIn;
	char strLine[MED_LINE_LENGTH];
	char strSeqLine[LONG_LINE_LENGTH];
	struct tagSequence *pNewSeq;
	char strSeqAlias[LINE_LENGTH];
	int ni,nj;
	
	/* other variables */

	/* #################################### */
	/* load sequences and conservation      */
	/* #################################### */
	sprintf(strLine, "%s%s", strInputPath, strSeqFile);
	fpIn = NULL;
	fpIn = fopen(strLine, "r");
	if(fpIn == NULL)
	{
        printf("Error: MotifMap_ScanMatrix_Sequential_ComputeBackground, cannot open sequence file!\n");
		exit(EXIT_FAILURE);
	}

	/* read */
	pNewSeq = NULL;
	while(fgets(strSeqLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strSeqLine);
		StrTrimRight(strSeqLine);
		if(strSeqLine[0] == '\0')
			continue;

		if(strSeqLine[0] == '>')
		{
			if(pNewSeq != NULL)
			{
				/* #################################### */
				/* load sequence and score              */
				/* #################################### */
				pSeqMtf = NULL;
				pSeqMtf = FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_SINGLE(pNewSeq,
					strInputPath, 0, "");
				SequenceDelete(pNewSeq);

				/* #################################### */
				/* compute background                   */
				/* #################################### */
				MotifMap_ScanMatrix_CountBackground(pSeqMtf, 0,
					0, pBG0, NULL);
				if(nBGOrder > 0)
				{
					MotifMap_ScanMatrix_CountBackground(pSeqMtf, 0,
						nBGOrder, pBGF, pBGB);
				}

				/* #################################### */
				/* destroy sequences.                   */
				/* #################################### */
				FLEXSEQMTFDESTROY(pSeqMtf);
				nSeqCount++;
			}

			pNewSeq = NULL;
			pNewSeq = SequenceCreate();
			if(pNewSeq == NULL)
			{
				fclose(fpIn);
				printf("Error: MotifMap_ScanMatrix_Sequential_Main, cannot allocate memory for loading sequences!\n");
				exit(EXIT_FAILURE);
			}
			pNewSeq->m_nIndex = nSeqCount;
			strcpy(pNewSeq->m_strAlias, (strSeqLine+1));
			strcpy(strSeqAlias, (strSeqLine+1));
		}
		else
		{
			if(SequenceAddTail(pNewSeq, strSeqLine) == PROC_FAILURE)
			{
				fclose(fpIn);
				printf("Error: MotifMap_ScanMatrix_Sequential_Main, cannot allocate memory for loading sequences!\n");
				exit(EXIT_FAILURE);
			}			
		}
	}

	/* the last sequence */
	if(pNewSeq != NULL)
	{
		/* #################################### */
		/* load sequence and score              */
		/* #################################### */
		pSeqMtf = NULL;
		pSeqMtf = FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_SINGLE(pNewSeq,
			strInputPath, 0, "");
		SequenceDelete(pNewSeq);

		/* #################################### */
		/* compute background                   */
		/* #################################### */
		MotifMap_ScanMatrix_CountBackground(pSeqMtf, 0,
			0, pBG0, NULL);
		if(nBGOrder > 0)
		{
			MotifMap_ScanMatrix_CountBackground(pSeqMtf, 0,
				nBGOrder, pBGF, pBGB);
		}

		/* #################################### */
		/* destroy sequences.                   */
		/* #################################### */
		FLEXSEQMTFDESTROY(pSeqMtf);
		nSeqCount++;
	}

	/* close files */
	fclose(fpIn);
	
	/* #################################### */
	/* normalization                        */
	/* #################################### */
	pEle1 = pBG0->pMatElement;
	for(ni=0; ni<pBG0->nHeight; ni++)
	{
		pEle2 = pEle1;
		dSum = 0.0;
		for(nj=0; nj<pBG0->nWidth; nj++)
		{
			dSum += (*pEle2);
			pEle2++;
		}

		for(nj=0; nj<pBG0->nWidth; nj++)
		{
			*pEle1 = *pEle1/dSum;
			pEle1++;
		}
	}

	
	if(nBGOrder > 0)
	{
		pEle1 = pBGF->pMatElement;
		for(ni=0; ni<pBGF->nHeight; ni++)
		{
			pEle2 = pEle1;
			dSum = 0.0;
			for(nj=0; nj<pBGF->nWidth; nj++)
			{
				dSum += (*pEle2);
				pEle2++;
			}

			for(nj=0; nj<pBGF->nWidth; nj++)
			{
				*pEle1 = *pEle1/dSum;
				pEle1++;
			}
		}

		pEle1 = pBGB->pMatElement;
		for(ni=0; ni<pBGB->nHeight; ni++)
		{
			pEle2 = pEle1;
			dSum = 0.0;
			for(nj=0; nj<pBGB->nWidth; nj++)
			{
				dSum += (*pEle2);
				pEle2++;
			}

			for(nj=0; nj<pBGB->nWidth; nj++)
			{
				*pEle1 = *pEle1/dSum;
				pEle1++;
			}
		}
	}
	else
	{
		DMCOPY(pBGF, pBG0);
		DMCOPY(pBGB, pBG0);
	}

	/* return */
	return nSeqCount;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_CountBackground: count empirical observations for  */
/*  constructing background markov matrices.                               */ 
/* ----------------------------------------------------------------------- */
int MotifMap_ScanMatrix_CountBackground(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,
					int nBGOrder, struct DOUBLEMATRIX *pBGF, struct DOUBLEMATRIX *pBGB)
{
	/* define */
	int nBaseTypeNum,nScale;
	int nWordId,nWordBadLen;
	int nW,nWB;
	int nLen;
	unsigned char *pBase;
	/* unsigned char *pCS; */
	int ni,nj;
	double dTemp;

	/* init */
	nLen = pSeqMtf->vSeq[nIndex]->nWidth;
	pBase = pSeqMtf->vSeq[nIndex]->pMatElement;
	nBaseTypeNum = pBGF->nWidth;
	if( nBGOrder == 0)
	{
		nScale = 0;
	}
	else
	{
		nScale = (int)(pow((double)nBaseTypeNum, (double)(nBGOrder-1)));
	}
	
	/* if zero order Markov chain as background */
	if(nBGOrder == 0)
	{
		for(ni=0; ni<nLen; ni++)
		{
			nW = (int)pBase[ni];
			if(nW < nBaseTypeNum)
			{
				dTemp = DMGETAT(pBGF, 0, nW)+1.0;
				DMSETAT(pBGF, 0, nW, dTemp);
			}
		}
	}

	/* if higher order Markov chain as background */
	else if(nLen > nBGOrder)
	{
		nWordId = 0;
		nWordBadLen = 0;
		
		/* initial words */
		for(ni=0; ni<nBGOrder; ni++)
		{
			nW = (int)pBase[ni];
			if(nW < nBaseTypeNum)
			{
				nWordId = nWordId*nBaseTypeNum+nW;
			}
			else
			{
				nWordId = nWordId*nBaseTypeNum;
				nWordBadLen++;
			}
		}

		/* continuing words */
		for(; ni<nLen; ni++)
		{
			/* get forward bg */
			nW = (int)pBase[ni];
			if(nW < nBaseTypeNum)
			{
				if(nWordBadLen == 0)
				{
					dTemp = DMGETAT(pBGF, nWordId, nW)+1.0;
					DMSETAT(pBGF, nWordId, nW, dTemp);
				}
				else
				{		
				}
			}
			else
			{
			}
			
			/* update word id */
			/* minus old letter */
			nj = ni-nBGOrder;
			nWB = (int)(pBase[nj]);
			if(nWB < nBaseTypeNum)
			{
				nWordId -= nWB*nScale;
			}
			else
			{
				nWordBadLen--;
			}
						
			/* add new letter */
			nW = (int)pBase[ni];
			if(nW < nBaseTypeNum)
			{
				nWordId = nWordId*nBaseTypeNum+nW;
			}
			else
			{
				nWordId = nWordId*nBaseTypeNum;
				nWordBadLen++;
			}

			/* get backward likelihood */
			if(nWB < nBaseTypeNum)
			{
				if(nWordBadLen == 0)
				{	
					dTemp = DMGETAT(pBGB, nWordId, nWB)+1.0;
					DMSETAT(pBGB, nWordId, nWB, dTemp);
				}
				else
				{
				}
			}
			else
			{
			}
		}
	}

	/* return */
	return PROC_SUCCESS;
}

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
					int nUseCS, double dC, char strCSAlias[])
{
	/* sequences and conservation */
	/* all sequences */
	struct FLEXSEQMOTIF **vSeqMtf = NULL;
	/* sequence number */
	int nSeqCount = 0;
	/* site number */
	int nEffecLen = 0;
	int nConsEffecLen = 0;
	int nSiteNum = 0;
	/* total sequence length and total site number */
	double dTotLen = 0.0;
	double dTotConsLen = 0.0;
	double dTotSite = 0.0;
	
	/* background */
	struct DOUBLEMATRIX *pBGF = NULL;
	struct DOUBLEMATRIX *pBGB = NULL;
	/* log background */
	struct DOUBLEMATRIX *pLogBGF = NULL;
	struct DOUBLEMATRIX *pLogBGB = NULL;
	/* background0 */
	struct DOUBLEMATRIX *pBG0 = NULL;
	/* log background0 */
	struct DOUBLEMATRIX *pLogBG0 = NULL;
	
	/* motif matrix */
	struct DOUBLEMATRIX *pMotif = NULL;
	
	/* number of base types */
	int nBaseTypeNum = 4;
	int nScale = 0;
	int nBGHeight;

	/* motif consensus */
	FILE *fpOut;
	FILE *fpOut2;
	char strLine[MED_LINE_LENGTH];
	
	/* other variables */
	int ni,nj;
	double *pEle1,*pEle2;
	double dSum;
	int nlen;
	

	/* #################################### */
	/* initialize                           */
	/* #################################### */
	nBaseTypeNum = 4;
	nScale = (int)(pow((double)nBaseTypeNum, (double)(nBGOrder-1)));
	nBGHeight = (int)(pow((double)nBaseTypeNum, (double)nBGOrder));
	
	/* #################################### */
	/* load initial parameters              */
	/* #################################### */
	if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
	{
		nlen = strlen(strInputPath);
		if(nlen == 0)
		{
			sprintf(strInputPath, ".\\"); 
		}
		else if(strInputPath[nlen-1] != '\\')
		{
			strInputPath[nlen] = '\\';
			strInputPath[nlen+1] = '\0';
		}
	}
	else
	{
		nlen = strlen(strInputPath);
		if(nlen == 0)
		{
			sprintf(strInputPath, "./"); 
		}
		else if(strInputPath[nlen-1] != '/')
		{
			strInputPath[nlen] = '/';
			strInputPath[nlen+1] = '\0';
		}
	}

	/* motif PWM */
	pMotif = NULL;
	pMotif = DMLOAD(strMotifPath);
	if( pMotif == NULL)
	{
		printf("Error: MotifMap_ScanMatrix_Sequential_Main, cannot load init motif pseudocount!\n");
		exit(EXIT_FAILURE);
	}
	pEle1 = pMotif->pMatElement;
	for(ni=0; ni<pMotif->nHeight; ni++)
	{
		pEle2 = pEle1;
		dSum = 0.0;
		for(nj=0; nj<pMotif->nWidth; nj++)
		{
			*pEle2 += 1e-3;
			dSum += (*pEle2);
			pEle2++;
		}

		for(nj=0; nj<pMotif->nWidth; nj++)
		{
			*pEle1 = log(*pEle1/dSum);
			pEle1++;
		}
	}

	/* #################################### */
	/* load sequences                       */
	/* #################################### */
	vSeqMtf = NULL;
	vSeqMtf = FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_GROUP(&nSeqCount, strInputPath, strSeqFile, nUseCS, strCSAlias);

	/* #################################### */
	/* compute background markov chain      */
	/* #################################### */
	pBGF = NULL;
	pBGF = CreateDoubleMatrix(nBGHeight, nBaseTypeNum);
	pBGB = NULL;
	pBGB = CreateDoubleMatrix(nBGHeight, nBaseTypeNum);
	pBG0 = NULL;
	pBG0 = CreateDoubleMatrix(1, nBaseTypeNum);
	for(ni=0; ni<nBGHeight; ni++)
	{
		for(nj=0; nj<nBaseTypeNum; nj++)
		{
			DMSETAT(pBGF, ni, nj, 1.0);
			DMSETAT(pBGB, ni, nj, 1.0);
		}
	}
	for(nj=0; nj<nBaseTypeNum; nj++)
	{
		DMSETAT(pBG0, 0, nj, 1.0);
	}
	
	if( (pBGF == NULL) || (pBGB == NULL) || (pBG0 == NULL) )
	{
		printf("Error: MotifMap_ScanMatrix_Sequential_Main, failure to create background markov matrix!\n");
		exit(EXIT_FAILURE);
	}

	/* compute background */
	for(ni=0; ni<nSeqCount; ni++)
	{
		MotifMap_ScanMatrix_CountBackground(vSeqMtf[ni], 0,
				0, pBG0, NULL);
		if(nBGOrder > 0)
		{
			MotifMap_ScanMatrix_CountBackground(vSeqMtf[ni], 0,
				nBGOrder, pBGF, pBGB);
		}
	}
	
	DMROWNORM(pBG0);

	if(nBGOrder > 0)
	{
		DMROWNORM(pBGF);
		DMROWNORM(pBGB);
	}
	else
	{
		DMCOPY(pBGF, pBG0);
		DMCOPY(pBGB, pBG0);
	}

	pLogBGF = NULL;
	pLogBGF = DMCLONE(pBGF);
	DMLOGTS(pLogBGF);

	pLogBGB = NULL;
	pLogBGB = DMCLONE(pBGB);
	DMLOGTS(pLogBGB);

	pLogBG0 = NULL;
	pLogBG0 = DMCLONE(pBG0);
	DMLOGTS(pLogBG0);

	sprintf(strLine, "%s.bg0", strOutputPath);
	DMSAVE(pBG0, strLine);
	sprintf(strLine, "%s.bgf", strOutputPath);
	DMSAVE(pBGF, strLine);
	sprintf(strLine, "%s.bgb", strOutputPath);
	DMSAVE(pBGB, strLine);

	fpOut = NULL;
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
        printf("Error: MotifMap_ScanMatrix_Sequential_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	
	sprintf(strLine, "%s.stat", strOutputPath);
	fpOut2 = NULL;
	fpOut2 = fopen(strLine, "w");
	if(fpOut2 == NULL)
	{
        printf("Error: MotifMap_ScanMatrix_Sequential_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	/* #################################### */
	/* scan for motif matrices              */
	/* #################################### */
	for(ni=0; ni<nSeqCount; ni++)
	{
		/* #################################### */
		/* init background loglikelihood        */
		/* #################################### */
		MotifMap_ScanMatrix_InitBGLogLike(vSeqMtf[ni], 0,
			nBGOrder, pLogBGF, pLogBGB, pLogBG0); 
		
		/* #################################### */
		/* scan motif                           */
		/* #################################### */
		MotifMap_ScanMatrix_In_FlexSeqMtf(vSeqMtf[ni], 0,
			0, vSeqMtf[ni]->nSeqLen-1, pMotif, dR, 
			nUseCS, dC, &nEffecLen, &nConsEffecLen, &nSiteNum);
		dTotLen += nEffecLen;
		dTotConsLen += nConsEffecLen;
		dTotSite += nSiteNum;

		/* #################################### */
		/* export results                       */
		/* #################################### */
		MotifMap_ScanMatrix_Export_Single(vSeqMtf[ni], 0,
			pMotif->nHeight, vSeqMtf[ni]->strAlias, fpOut);
	}
	
	/* save statistics */
	fprintf(fpOut2, "EffecLen= %d\n", (int)dTotLen);
	fprintf(fpOut2, "TotalSite= %d\n", (int)dTotSite);
	fprintf(fpOut2, "ConsLen= %d\n", (int)dTotConsLen);
	
	/* close files */
	fclose(fpOut);
	fclose(fpOut2);

	/* #################################### */
	/* release memory                       */
	/* #################################### */
	for(ni=0; ni<nSeqCount; ni++)
	{
		FLEXSEQMTFDESTROY(vSeqMtf[ni]);
		vSeqMtf[ni] = NULL;
	}
	free(vSeqMtf);

	DestroyDoubleMatrix(pMotif);
	DestroyDoubleMatrix(pBGF);
	DestroyDoubleMatrix(pBGB);
	DestroyDoubleMatrix(pLogBGF);
	DestroyDoubleMatrix(pLogBGB);
	DestroyDoubleMatrix(pBG0);
	DestroyDoubleMatrix(pLogBG0);
	

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_InitBGLogLike_Genome: Initialize background        */
/*  loglikelihood using genomic local MC matrices.                         */
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanMatrix_InitBGLogLike_Genome(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,
					char strChr[], int nChrStart, int nChrEnd, 
					int nBGOrder, char strBGPath[], int nBGStepSize)
{
	/* define */
	int nBaseTypeNum,nScale;
	int nWordId,nWordBadLen;
	int nW,nWB;
	int nLen;
	unsigned char *pBase;
	/* unsigned char *pCS; */
	double *pBLikeF,*pBLikeB;
	int ni,nj;
	/* for accessing background matrices */
	int nBinidx,nBinidxF,nBinidxB,nMatidx,nMatidxF,nMatidxB;
	char strFileName[MED_LINE_LENGTH];
	char strFileName1[MED_LINE_LENGTH];
	char strFileName2[MED_LINE_LENGTH];
	char strFileName3[MED_LINE_LENGTH];
	struct DOUBLEMATRIX *pLogBGF,*pLogBGB,*pLogBG0,*pLogBG0B;

	/* init */
	nLen = pSeqMtf->vSeq[nIndex]->nWidth;
	pBase = pSeqMtf->vSeq[nIndex]->pMatElement;
	pBLikeF = pSeqMtf->vMonitor[0]->pMatElement;
	pBLikeB = pSeqMtf->vMonitor[1]->pMatElement;

	nBaseTypeNum = 4;
	if( nBGOrder == 0)
	{
		nScale = 0;
	}
	else
	{
		nScale = (int)(pow((double)nBaseTypeNum, (double)(nBGOrder-1)));
	}
	
	/* if zero order Markov chain as background */
	if((nBGOrder == 0) || (nLen <= nBGOrder))
	{
		nMatidx = nChrStart/nBGStepSize;
		if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
		{
			sprintf(strFileName, "%s0\\%s\\%d.txt", strBGPath, strChr, nMatidx);
		}
		else
		{
			sprintf(strFileName, "%s0/%s/%d.txt", strBGPath, strChr, nMatidx);
		}
		pLogBG0 = NULL;
		pLogBG0 = DMLOAD(strFileName);
		DMLOGTS(pLogBG0);

		for(ni=0; ni<nLen; ni++)
		{
			nBinidx = (nChrStart+ni)/nBGStepSize;
			if(nBinidx != nMatidx)
			{
				DestroyDoubleMatrix(pLogBG0);
				nMatidx = nBinidx;
				if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
				{
					sprintf(strFileName, "%s0\\%s\\%d.txt", strBGPath, strChr, nMatidx);
				}
				else
				{
					sprintf(strFileName, "%s0/%s/%d.txt", strBGPath, strChr, nMatidx);
				}
				pLogBG0 = NULL;
				pLogBG0 = DMLOAD(strFileName);
				DMLOGTS(pLogBG0);
			}

			nW = (int)pBase[ni];
			if(nW < nBaseTypeNum)
			{
				pBLikeF[ni] = DMGETAT(pLogBG0, 0, nW);
				pBLikeB[ni] = pBLikeF[ni];
			}
			else
			{
				pBLikeF[ni] = 0.0;
				pBLikeB[ni] = 0.0;
			}
		}

		DestroyDoubleMatrix(pLogBG0);

		if( (nChrStart+ni-1) != nChrEnd )
		{
			printf("Error: MotifMap_ScanMatrix_InitBGLogLike_Genome, coordinates do not match!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* if higher order Markov chain as background */
	else
	{
		nMatidxF = nChrStart/nBGStepSize;
		nMatidxB = nMatidxF;
		if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
		{
			sprintf(strFileName, "%s0\\%s\\%d.txt", strBGPath, strChr, nMatidxF);
			sprintf(strFileName3, "%s0\\%s\\%d.txt", strBGPath, strChr, nMatidxB);
			sprintf(strFileName1, "%s%d\\%s\\%d_f.txt", strBGPath, nBGOrder, strChr, nMatidxF);
			sprintf(strFileName2, "%s%d\\%s\\%d_b.txt", strBGPath, nBGOrder, strChr, nMatidxB);
		}
		else
		{
			sprintf(strFileName, "%s0/%s/%d.txt", strBGPath, strChr, nMatidxF);
			sprintf(strFileName3, "%s0/%s/%d.txt", strBGPath, strChr, nMatidxB);
			sprintf(strFileName1, "%s%d/%s/%d_f.txt", strBGPath, nBGOrder, strChr, nMatidxF);
			sprintf(strFileName2, "%s%d/%s/%d_b.txt", strBGPath, nBGOrder, strChr, nMatidxB);
		}
		pLogBG0 = NULL;
		pLogBG0 = DMLOAD(strFileName);
		DMLOGTS(pLogBG0);
		pLogBG0B = NULL;
		pLogBG0B = DMLOAD(strFileName3);
		DMLOGTS(pLogBG0B);
		pLogBGF = NULL;
		pLogBGF = DMLOAD(strFileName1);
		DMLOGTS(pLogBGF);
		pLogBGB = NULL;
		pLogBGB = DMLOAD(strFileName2);
		DMLOGTS(pLogBGB);


		nWordId = 0;
		nWordBadLen = 0;
		
		/* initial words */
		for(ni=0; ni<nBGOrder; ni++)
		{
			nBinidxF = (nChrStart+ni)/nBGStepSize;
			if(nBinidxF != nMatidxF)
			{
				DestroyDoubleMatrix(pLogBG0);
				DestroyDoubleMatrix(pLogBGF);
				
				nMatidxF = nBinidxF;
				if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
				{
					sprintf(strFileName, "%s0\\%s\\%d.txt", strBGPath, strChr, nMatidxF);
					sprintf(strFileName1, "%s%d\\%s\\%d_f.txt", strBGPath, nBGOrder, strChr, nMatidxF);
				}
				else
				{
					sprintf(strFileName, "%s0/%s/%d.txt", strBGPath, strChr, nMatidxF);
					sprintf(strFileName1, "%s%d/%s/%d_f.txt", strBGPath, nBGOrder, strChr, nMatidxF);
				}
				pLogBG0 = NULL;
				pLogBG0 = DMLOAD(strFileName);
				DMLOGTS(pLogBG0);
				pLogBGF = NULL;
				pLogBGF = DMLOAD(strFileName1);
				DMLOGTS(pLogBGF);
			}

			nW = (int)pBase[ni];
			if(nW < nBaseTypeNum)
			{
				pBLikeF[ni] = DMGETAT(pLogBG0, 0, nW);
				nWordId = nWordId*nBaseTypeNum+nW;
			}
			else
			{
				pBLikeF[ni] = 0.0;
				nWordId = nWordId*nBaseTypeNum;
				nWordBadLen++;
			}
		}

		/* continuing words */
		for(; ni<nLen; ni++)
		{
			nBinidxF = (nChrStart+ni)/nBGStepSize;
			if(nBinidxF != nMatidxF)
			{
				DestroyDoubleMatrix(pLogBG0);
				DestroyDoubleMatrix(pLogBGF);
				
				nMatidxF = nBinidxF;
				if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
				{
					sprintf(strFileName, "%s0\\%s\\%d.txt", strBGPath, strChr, nMatidxF);
					sprintf(strFileName1, "%s%d\\%s\\%d_f.txt", strBGPath, nBGOrder, strChr, nMatidxF);
				}
				else
				{
					sprintf(strFileName, "%s0/%s/%d.txt", strBGPath, strChr, nMatidxF);
					sprintf(strFileName1, "%s%d/%s/%d_f.txt", strBGPath, nBGOrder, strChr, nMatidxF);
				}
				pLogBG0 = NULL;
				pLogBG0 = DMLOAD(strFileName);
				DMLOGTS(pLogBG0);
				pLogBGF = NULL;
				pLogBGF = DMLOAD(strFileName1);
				DMLOGTS(pLogBGF);
			}

			/* get forward likelihood */
			nW = (int)pBase[ni];
			if(nW < nBaseTypeNum)
			{
				if(nWordBadLen == 0)
					pBLikeF[ni] = DMGETAT(pLogBGF, nWordId, nW);
				else
					pBLikeF[ni] = DMGETAT(pLogBG0, 0, nW);
			}
			else
			{
				pBLikeF[ni] = 0.0;
			}
			
			/* update word id */
			/* minus old letter */
			nj = ni-nBGOrder;
			nWB = (int)(pBase[nj]);
			if(nWB < nBaseTypeNum)
			{
				nWordId -= nWB*nScale;
			}
			else
			{
				nWordBadLen--;
			}
						
			/* add new letter */
			nW = (int)pBase[ni];
			if(nW < nBaseTypeNum)
			{
				nWordId = nWordId*nBaseTypeNum+nW;
			}
			else
			{
				nWordId = nWordId*nBaseTypeNum;
				nWordBadLen++;
			}

			/* get backward likelihood */
			nBinidxB = (nChrStart+nj)/nBGStepSize;
			if(nBinidxB != nMatidxB)
			{
				DestroyDoubleMatrix(pLogBG0B);
				DestroyDoubleMatrix(pLogBGB);
				
				nMatidxB = nBinidxB;
				if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
				{
					sprintf(strFileName3, "%s0\\%s\\%d.txt", strBGPath, strChr, nMatidxB);
					sprintf(strFileName2, "%s%d\\%s\\%d_b.txt", strBGPath, nBGOrder, strChr, nMatidxB);
				}
				else
				{
					sprintf(strFileName3, "%s0/%s/%d.txt", strBGPath, strChr, nMatidxB);
					sprintf(strFileName2, "%s%d/%s/%d_b.txt", strBGPath, nBGOrder, strChr, nMatidxB);
				}
				pLogBG0B = NULL;
				pLogBG0B = DMLOAD(strFileName3);
				DMLOGTS(pLogBG0B);
				pLogBGB = NULL;
				pLogBGB = DMLOAD(strFileName2);
				DMLOGTS(pLogBGB);
			}

			if(nWB < nBaseTypeNum)
			{
				if(nWordBadLen == 0)
					pBLikeB[nj] = DMGETAT(pLogBGB, nWordId, nWB);
				else
					pBLikeB[nj] = DMGETAT(pLogBG0B, 0, nWB);
			}
			else
			{
				pBLikeB[nj] = 0.0;
			}
		}

		nj++;
		for(; nj<nLen; nj++)
		{
			nBinidxB = (nChrStart+nj)/nBGStepSize;
			if(nBinidxB != nMatidxB)
			{
				DestroyDoubleMatrix(pLogBG0B);
				DestroyDoubleMatrix(pLogBGB);
				
				nMatidxB = nBinidxB;
				if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
				{
					sprintf(strFileName3, "%s0\\%s\\%d.txt", strBGPath, strChr, nMatidxB);
					sprintf(strFileName2, "%s%d\\%s\\%d_b.txt", strBGPath, nBGOrder, strChr, nMatidxB);
				}
				else
				{
					sprintf(strFileName3, "%s0/%s/%d.txt", strBGPath, strChr, nMatidxB);
					sprintf(strFileName2, "%s%d/%s/%d_b.txt", strBGPath, nBGOrder, strChr, nMatidxB);
				}
				pLogBG0B = NULL;
				pLogBG0B = DMLOAD(strFileName3);
				DMLOGTS(pLogBG0B);
				pLogBGB = NULL;
				pLogBGB = DMLOAD(strFileName2);
				DMLOGTS(pLogBGB);
			}

			nWB = (int)pBase[nj];
			if(nWB < nBaseTypeNum)
			{
				pBLikeB[nj] = DMGETAT(pLogBG0B, 0, nWB);
			}
			else
			{
				pBLikeB[nj] = 0.0;
			}
		}

		DestroyDoubleMatrix(pLogBG0);
		DestroyDoubleMatrix(pLogBG0B);
		DestroyDoubleMatrix(pLogBGF);
		DestroyDoubleMatrix(pLogBGB);

		if( ((nChrStart+nj-1) != nChrEnd) || ((nChrStart+ni-1) != nChrEnd) )
		{
			printf("Error: MotifMap_ScanMatrix_InitBGLogLike_Genome, coordinates do not match!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_InitBGLogLike: Initialize background loglikelihood */
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanMatrix_InitBGLogLike(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,
			int nBGOrder, struct DOUBLEMATRIX *pLogBGF, struct DOUBLEMATRIX *pLogBGB,
			struct DOUBLEMATRIX *pLogBG0)
{
	/* define */
	int nBaseTypeNum,nScale;
	int nWordId,nWordBadLen;
	int nW,nWB;
	int nLen;
	unsigned char *pBase;
	/* unsigned char *pCS; */
	double *pBLikeF,*pBLikeB;
	int ni,nj;

	/* init */
	nLen = pSeqMtf->vSeq[nIndex]->nWidth;
	pBase = pSeqMtf->vSeq[nIndex]->pMatElement;
	pBLikeF = pSeqMtf->vMonitor[0]->pMatElement;
	pBLikeB = pSeqMtf->vMonitor[1]->pMatElement;

	nBaseTypeNum = pLogBG0->nWidth;
	if( nBGOrder == 0)
	{
		nScale = 0;
	}
	else
	{
		nScale = (int)(pow((double)nBaseTypeNum, (double)(nBGOrder-1)));
	}
	
	/* if zero order Markov chain as background */
	if((nBGOrder == 0) || (nLen <= nBGOrder))
	{
		for(ni=0; ni<nLen; ni++)
		{
			nW = (int)pBase[ni];
			if(nW < nBaseTypeNum)
			{
				pBLikeF[ni] = DMGETAT(pLogBG0, 0, nW);
				pBLikeB[ni] = pBLikeF[ni];
			}
			else
			{
				pBLikeF[ni] = 0.0;
				pBLikeB[ni] = 0.0;
			}
		}
	}

	/* if higher order Markov chain as background */
	else
	{
		if( (pLogBGF->nWidth != pLogBGB->nWidth) || (pLogBGF->nHeight != pLogBGB->nHeight) )
		{
			printf("Error: MotifMap_ScanMatrix_InitBGLogLike, background matrix size not match!\n");
			exit(EXIT_FAILURE);
		}

		nWordId = 0;
		nWordBadLen = 0;
		
		/* initial words */
		for(ni=0; ni<nBGOrder; ni++)
		{
			nW = (int)pBase[ni];
			if(nW < nBaseTypeNum)
			{
				pBLikeF[ni] = DMGETAT(pLogBG0, 0, nW);
				nWordId = nWordId*nBaseTypeNum+nW;
			}
			else
			{
				pBLikeF[ni] = 0.0;
				nWordId = nWordId*nBaseTypeNum;
				nWordBadLen++;
			}
		}

		/* continuing words */
		for(; ni<nLen; ni++)
		{
			/* get forward likelihood */
			nW = (int)pBase[ni];
			if(nW < nBaseTypeNum)
			{
				if(nWordBadLen == 0)
					pBLikeF[ni] = DMGETAT(pLogBGF, nWordId, nW);
				else
					pBLikeF[ni] = DMGETAT(pLogBG0, 0, nW);
			}
			else
			{
				pBLikeF[ni] = 0.0;
			}
			
			/* update word id */
			/* minus old letter */
			nj = ni-nBGOrder;
			nWB = (int)(pBase[nj]);
			if(nWB < nBaseTypeNum)
			{
				nWordId -= nWB*nScale;
			}
			else
			{
				nWordBadLen--;
			}
						
			/* add new letter */
			nW = (int)pBase[ni];
			if(nW < nBaseTypeNum)
			{
				nWordId = nWordId*nBaseTypeNum+nW;
			}
			else
			{
				nWordId = nWordId*nBaseTypeNum;
				nWordBadLen++;
			}

			/* get backward likelihood */
			if(nWB < nBaseTypeNum)
			{
				if(nWordBadLen == 0)
					pBLikeB[nj] = DMGETAT(pLogBGB, nWordId, nWB);
				else
					pBLikeB[nj] = DMGETAT(pLogBG0, 0, nWB);
			}
			else
			{
				pBLikeB[nj] = 0.0;
			}
		}

		nj++;
		for(; nj<nLen; nj++)
		{
			nWB = (int)pBase[nj];
			if(nWB < nBaseTypeNum)
			{
				pBLikeB[nj] = DMGETAT(pLogBG0, 0, nWB);
			}
			else
			{
				pBLikeB[nj] = 0.0;
			}
		}
		
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_In_FlexSeqMtf: call motif sites.                   */
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanMatrix_In_FlexSeqMtf(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,
			int nFrom, int nTo, 
			struct DOUBLEMATRIX *pLogPWM, double dR, 
			int nUseCS, double dC, 
			int *pEffecLen, int *pConsEffecLen, int *pSiteNum)
{
	/* define */
	unsigned char *pBase;
	unsigned char *pCS;
	int ni,nj,nLen;
	int nBaseTypeNum;
	struct FLEXMOTIFSITE *pPrev,*pSite;
	double dLike[2];
	int nMask[2];
	double dMaxLike;
	int nStrand;
	int nMotifLen;
	int nActualFrom,nActualTo;

	double dTotConserve;
	double dAveConserve;
	
	/* init */
	if(nFrom > nTo)
	{
		printf("Error: MotifMap_ScanMatrix_In_FlexSeqMtf, coordinates wrong, region start > region end!\n");
		exit(EXIT_FAILURE);
	}

	nBaseTypeNum = pLogPWM->nWidth;
	nMotifLen = pLogPWM->nHeight;
	nLen = pSeqMtf->vSeq[nIndex]->nWidth;

	if( (nLen<nMotifLen) || (nFrom+nMotifLen) > nLen)
	{
		printf("Warning: MotifMap_ScanMatrix_In_FlexSeqMtf, seqence %d length < motif length!\n", pSeqMtf->nId);
		return PROC_SUCCESS;
	}
	else
	{
		nActualFrom = nFrom;
	}

	if( (nTo+nMotifLen) > nLen)
	{
		nActualTo = nLen-nMotifLen;
		if(nActualTo < nActualFrom)
		{
			printf("Error: MotifMap_ScanMatrix_In_FlexSeqMtf, coordinates wrong, region start > region end!\n");
			exit(EXIT_FAILURE);
		}
		else
		{
			/* printf("Warning: MotifMap_ScanMatrix_In_FlexSeqMtf, sequence %d actual to %d!\n", pSeqMtf->nId, nActualTo); */
		}
	}
	else
	{
		nActualTo = nTo;
	}
	
	*pEffecLen = 0;
	*pConsEffecLen = 0;
	*pSiteNum = 0;
	
	pBase = pSeqMtf->vSeq[nIndex]->pMatElement;
	pPrev = NULL;
	
	dTotConserve = 0.0;
	if(nUseCS == 1)
	{
		pCS = pSeqMtf->vScore[0]->pMatElement;
		for(nj=0; nj<(nMotifLen-1); nj++)
		{
			dTotConserve += pCS[nActualFrom+nj];
		}
	}


	for(ni=nActualFrom; ni<=nActualTo; ni++)
	{
		if(nUseCS == 1)
		{
			dTotConserve += pCS[ni+nMotifLen-1];
		}

		/* if repeat or 'N', skip it */
		if(pBase[ni] >= nBaseTypeNum)
		{
			if(nUseCS == 1)
			{
				dTotConserve -= pCS[ni];
			}
			continue;
		}

		/* add effect length */
		*pEffecLen += 1;
		
		/* compute conservation */
		if(nUseCS == 1)
		{
			dAveConserve = dTotConserve/(double)nMotifLen;
			if(dAveConserve < dC)
			{
				dTotConserve -= pCS[ni];
				continue;
			}
			else
			{
				*pConsEffecLen += 1;
			}
		}
		else
		{
			*pConsEffecLen += 1;
		}


		/* get motif likelihood and set posterior probability */
		dLike[0] = MotifMap_ScanMatrix_MotifLikeRatio(pSeqMtf, nIndex, ni, 
					pLogPWM, '+', (nMask+0));
		dLike[1] = MotifMap_ScanMatrix_MotifLikeRatio(pSeqMtf, nIndex, ni, 
					pLogPWM, '-', (nMask+1));
		
		if( (nMask[0]==1) && (nMask[1]==1) )
		{
			if(nUseCS == 1)
			{
				dTotConserve -= pCS[ni];
			}
			continue;
		}
		else if(nMask[0]==1)
		{
			dMaxLike = dLike[1];
			nStrand = 1;
		}
		else if(nMask[1]==1)
		{
			dMaxLike = dLike[0];
			nStrand = 0;
		}
		else
		{
			if(dLike[0] > dLike[1])
			{
				dMaxLike = dLike[0];
				nStrand = 0;
			}
			else
			{
				dMaxLike = dLike[1];
				nStrand = 1;
			}
		}

		dMaxLike = exp(dMaxLike);
		
		/* create new site if above threshold */
		if(dMaxLike >= dR)
		{
			pSite = NULL;
			pSite = FLEXMTFSCREATE();
			if(pSite == NULL)
			{
				printf("Error: MotifMap_ScanMatrix_In_FlexSeqMtf, cannot create new site!\n");
				exit(EXIT_FAILURE);
			}

			pSite->nSeqId = pSeqMtf->nId;
			pSite->nStartPos = ni;
			pSite->nStrand = nStrand;
			pSite->pNext = NULL;
			pSite->pPrev = NULL;
			pSite->dScore = log10(dMaxLike);

			if(pPrev == NULL)
			{
				pSeqMtf->pMotifList = pSite;
				pPrev = pSite;
			}
			else
			{
				pPrev->pNext = pSite;
				pSite->pPrev = pPrev;
				pPrev = pSite;
			}

			pSeqMtf->nSiteNum += 1;
			*pSiteNum += 1;
		}

		if(nUseCS == 1)
		{
			dTotConserve -= pCS[ni];
		}
	}


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_MotifLikeRatio: get loglikelihood ratio for motif. */
/* ----------------------------------------------------------------------- */ 
double MotifMap_ScanMatrix_MotifLikeRatio(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, int nPos, 
					struct DOUBLEMATRIX *pLogPWM, char chStrand,
					int *pMask)
{
	/* define */
	double dL,dLF,dLB,dTemp;
	int ni,nMotifLen,nW;
	int nBaseTypeNum;
	int nSeqLen;
	unsigned char *pBase;
	double *pBLF,*pBLB;
	
	/* init */
	*pMask = 0;
	dL = 0.0;
	dLF = 0.0;
	dLB = 0.0;
	
	/* check */
	if(pLogPWM == NULL)
	{
		printf("Error: MotifMap_ScanMatrix_MotifLikeRatio, no PWM exist for likelihood calculation!\n");
		exit(EXIT_FAILURE);
	}
	nBaseTypeNum = pLogPWM->nWidth;
	
	/* get likelihood */
	nMotifLen = pLogPWM->nHeight;
	nSeqLen = pSeqMtf->vSeq[nIndex]->nWidth;
	if((nPos+nMotifLen) > nSeqLen)
	{
		*pMask = 1;
		return dL;
	}

	pBase = pSeqMtf->vSeq[nIndex]->pMatElement+nPos;
	pBLF = pSeqMtf->vMonitor[0]->pMatElement+nPos;
	pBLB = pSeqMtf->vMonitor[1]->pMatElement+nPos;

	if(chStrand == '+')
	{
		for(ni=0; ni<nMotifLen; ni++)
		{
			nW = (int)(pBase[ni]);
			
			if(nW >= nBaseTypeNum)
			{
				*pMask = 1;
				dL = 0.0;
				return dL;
			}

			dTemp = DMGETAT(pLogPWM, ni, nW);
			dLF += dTemp;
			dLF -= pBLF[ni];
			dLB += dTemp;
			dLB -= pBLB[ni];
		}

		if(dLF < dLB)
			dL = dLF;
		else
			dL = dLB;
	}
	else if(chStrand == '-')
	{
		for(ni=0; ni<nMotifLen; ni++)
		{
			nW = (int)(pBase[ni]);
			
			if(nW >= nBaseTypeNum)
			{
				*pMask = 1;
				dL = 0.0;
				return dL;
			}

			nW = nBaseTypeNum-1-nW;
			dTemp = DMGETAT(pLogPWM, (nMotifLen-1-ni), nW);
			dLF += dTemp;
			dLF -= pBLF[ni];
			dLB += dTemp;
			dLB -= pBLB[ni];
		}
		
		if(dLF < dLB)
			dL = dLF;
		else
			dL = dLB;
	}
	else
	{
		printf("Error: MotifMap_ScanMatrix_MotifLikeRatio, strand information wrong!\n");
		exit(EXIT_FAILURE);
	}

	/* return */
	return dL;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_Export_Single: save motif sites to a file.         */
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanMatrix_Export_Single(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,
						int nMotifLen, char strSeqAlias[], FILE *fpOut)
{
	/* define */
	int nj;
	struct FLEXMOTIFSITE *pSite;
	char chStrand;
	char chBase;
	unsigned char *pBase;

	/* export */
	pSite = pSeqMtf->pMotifList;
	while(pSite != NULL)
	{
		if(pSite->nStrand == 0)
			chStrand = '+';
		else
			chStrand = '-';

		/* fprintf(fpOut, ">%d:%d-%d(%c)\n", pSeqMtf->nId, pSite->nStartPos, pSite->nStartPos+nMotifLen-1,
			chStrand); */

		fprintf(fpOut, "%s\t%d\t%d\t%d\t%c\t%f\t", strSeqAlias, pSeqMtf->nId, pSite->nStartPos, pSite->nStartPos+nMotifLen-1,
			chStrand, pSite->dScore);

		if(pSite->nStrand == 0)
		{
			pBase = pSeqMtf->vSeq[nIndex]->pMatElement+pSite->nStartPos;
			for(nj=0; nj<nMotifLen; nj++)
			{
				switch(*pBase)
				{
					case 0: chBase = 'A';
						break;
					case 1: chBase = 'C';
						break;
					case 2: chBase = 'G';
						break;
					case 3: chBase = 'T';
						break;
					case 10: chBase = 'a';
						break;
					case 11: chBase = 'c';
						break;
					case 12: chBase = 'g';
						break;
					case 13: chBase = 't';
						break;
					default: chBase = 'N';
				}

				fprintf(fpOut, "%c", chBase);
				pBase++;
			}
			fprintf(fpOut, "\n");
		}
		else
		{
			pBase = pSeqMtf->vSeq[nIndex]->pMatElement+pSite->nStartPos+nMotifLen-1;
			for(nj=0; nj<nMotifLen; nj++)
			{
				switch(*pBase)
				{
					case 0: chBase = 'T';
						break;
					case 1: chBase = 'G';
						break;
					case 2: chBase = 'C';
						break;
					case 3: chBase = 'A';
						break;
					case 10: chBase = 't';
						break;
					case 11: chBase = 'g';
						break;
					case 12: chBase = 'c';
						break;
					case 13: chBase = 'a';
						break;
					default: chBase = 'N';
				}

				fprintf(fpOut, "%c", chBase);
				pBase--;
			}
			fprintf(fpOut, "\n");
		}

		pSite = pSite->pNext;
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_Export: export matrix mapping result.              */
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanMatrix_Export(int nSeqCount, struct FLEXSEQMOTIF **vSeqMtf, int nIndex,
			int nMotifLen, double dTotLen, double dTotSite, char strOutputPath[])
{
	/* define */
	FILE *fpOut;
	int ni,nj;
	struct FLEXMOTIFSITE *pSite;
	char chStrand;
	char chBase;
	unsigned char *pBase;

	/* open file */
	fpOut = NULL;
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: MotifMap_ScanMatrix_Export, cannot open file!\n");
		exit(EXIT_FAILURE);
	}

	/* export */
	fprintf(fpOut, "EffecLen= %d\n", (int)dTotLen);
	fprintf(fpOut, "TotalSite= %d\n", (int)dTotSite);
	fprintf(fpOut, "\n");
	for(ni=0; ni<nSeqCount; ni++)
	{
		pSite = vSeqMtf[ni]->pMotifList;
		while(pSite != NULL)
		{
			if(pSite->nStrand == 0)
				chStrand = '+';
			else
				chStrand = '-';

			/* fprintf(fpOut, ">%d:%d-%d(%c)\n", ni, pSite->nStartPos, pSite->nStartPos+nMotifLen-1,
				chStrand); */

			fprintf(fpOut, "%d\t%d\t%d\t%c\t", ni, pSite->nStartPos, pSite->nStartPos+nMotifLen-1,
				chStrand);

			if(pSite->nStrand == 0)
			{
				pBase = vSeqMtf[ni]->vSeq[nIndex]->pMatElement+pSite->nStartPos;
				for(nj=0; nj<nMotifLen; nj++)
				{
					switch(*pBase)
					{
						case 0: chBase = 'A';
							break;
						case 1: chBase = 'C';
							break;
						case 2: chBase = 'G';
							break;
						case 3: chBase = 'T';
							break;
						case 10: chBase = 'a';
							break;
						case 11: chBase = 'c';
							break;
						case 12: chBase = 'g';
							break;
						case 13: chBase = 't';
							break;
						default: chBase = 'N';
					}

					fprintf(fpOut, "%c", chBase);
					pBase++;
				}
				fprintf(fpOut, "\n");
			}
			else
			{
				pBase = vSeqMtf[ni]->vSeq[nIndex]->pMatElement+pSite->nStartPos+nMotifLen-1;
				for(nj=0; nj<nMotifLen; nj++)
				{
					switch(*pBase)
					{
						case 0: chBase = 'T';
							break;
						case 1: chBase = 'G';
							break;
						case 2: chBase = 'C';
							break;
						case 3: chBase = 'A';
							break;
						case 10: chBase = 't';
							break;
						case 11: chBase = 'g';
							break;
						case 12: chBase = 'c';
							break;
						case 13: chBase = 'a';
							break;
						default: chBase = 'N';
					}

					fprintf(fpOut, "%c", chBase);
					pBase--;
				}
				fprintf(fpOut, "\n");
			}

			pSite = pSite->pNext;
		}
	}

	/* close file */
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_Genome_GetCutoff_Main: compute cutoff to calibrate */
/*  motif occurence rate.                                                  */
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanMatrix_Genome_GetCutoff_Main(char strMotifListPath[], char strGenomePath[], 
					char strCodPath[], char strOutputPath[],
					int nBGOrder, char strBGType[], char strBGPath[], int nBGStepSize,
					int nUseCS, double dC, char strCSPath[], int nIncludeRepeat)
{
	/* define */
	FILE *fpIn = NULL;
	FILE *fpStat = NULL;
	FILE *fpOut = NULL;
	
	double dUpper,dLower;
	double dRate;
	double dLR0 = 0.0;
	double dLR = 500.0;

	int nMotifNum = 0;
	char strMotifFile[MED_LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];

	char strTempOutFile[MED_LINE_LENGTH];
	char strTempStatFile[MED_LINE_LENGTH];
	char strTemp[LINE_LENGTH];
	
	int nN[2];
	double dRatio;

	/* process one by one */
	fpIn = NULL;
	fpIn = fopen(strMotifListPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: MotifMap_ScanMatrix_Genome_GetCutoff_Main, cannot open motiflist file!\n");
		exit(EXIT_FAILURE);
	}
	fpOut = NULL;
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: MotifMap_ScanMatrix_Genome_GetCutoff_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	
	nMotifNum = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		sscanf(strLine, "%s %lf", strMotifFile, &dRate);
		printf("Processing motif %d, %s ...\n", nMotifNum, strMotifFile);
		if(dRate <= 0.0)
		{
			printf("Warning: cannot process rate = %f\n", dRate);
			continue;
		}
		
		/* get initial rate */
		dLR0 = 0.0;
		dLR = 500.0;
		sprintf(strTempOutFile, "%s.%d.map", strOutputPath, nMotifNum);
		MotifMap_ScanMatrix_Genome_Main(strMotifFile, strGenomePath, 
					strCodPath, strTempOutFile, dLR, 
					nBGOrder, strBGType, 
					strBGPath, nBGStepSize,
					nUseCS, dC, strCSPath, nIncludeRepeat);
		
		sprintf(strTempStatFile, "%s.stat", strTempOutFile);
		fpStat = NULL;
		fpStat = fopen(strTempStatFile, "r");
		if(fpStat == NULL)
		{
			printf("Error: MotifMap_ScanMatrix_Genome_GetCutoff_Main, cannot open temporary output file!\n");
			exit(EXIT_FAILURE);
		}

		fgets(strLine, LONG_LINE_LENGTH, fpStat);
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		sscanf(strLine, "%s %d", strTemp, nN+1);

		fgets(strLine, LONG_LINE_LENGTH, fpStat);
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		sscanf(strLine, "%s %d", strTemp, nN+0);
		fclose(fpStat);
		
		dRatio = (double)(nN[0])/((double)(nN[1])+1e-6);

		if( fabs(dRatio-dRate) < (0.01*dRate) )
		{
			fprintf(fpOut, "%s\t%d\n", strMotifFile, (int)dLR);
			nMotifNum++;
			continue;
		}

		/* get first upper bound and lower bound */
		if( dRatio>dRate )
		{
			while(dRatio>dRate)
			{
				dLower = dLR;
				dLR = 2.0*dLR;
				sprintf(strTempOutFile, "%s.%d.map", strOutputPath, nMotifNum);
				MotifMap_ScanMatrix_Genome_Main(strMotifFile, strGenomePath, 
							strCodPath, strTempOutFile, dLR, 
							nBGOrder, strBGType, 
							strBGPath, nBGStepSize,
							nUseCS, dC, strCSPath, nIncludeRepeat);
				
				sprintf(strTempStatFile, "%s.stat", strTempOutFile);
				fpStat = NULL;
				fpStat = fopen(strTempStatFile, "r");
				if(fpStat == NULL)
				{
					printf("Error: MotifMap_ScanMatrix_Genome_GetCutoff_Main, cannot open temporary output file!\n");
					exit(EXIT_FAILURE);
				}

				fgets(strLine, LONG_LINE_LENGTH, fpStat);
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				sscanf(strLine, "%s %d", strTemp, nN+1);

				fgets(strLine, LONG_LINE_LENGTH, fpStat);
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				sscanf(strLine, "%s %d", strTemp, nN+0);
				fclose(fpStat);
				
				dRatio = (double)(nN[0])/((double)(nN[1])+1e-6);
			}
			dUpper = dLR;
		}
		else 
		{
			while(dRatio<dRate)
			{
				dUpper = dLR;
				dLR = dLR/2.0;
				MotifMap_ScanMatrix_Genome_Main(strMotifFile, strGenomePath, 
							strCodPath, strTempOutFile, dLR, 
							nBGOrder, strBGType, 
							strBGPath, nBGStepSize,
							nUseCS, dC, strCSPath, nIncludeRepeat);
				
				sprintf(strTempStatFile, "%s.stat", strTempOutFile);
				fpStat = NULL;
				fpStat = fopen(strTempStatFile, "r");
				if(fpStat == NULL)
				{
					printf("Error: MotifMap_ScanMatrix_Genome_GetCutoff_Main, cannot open temporary output file!\n");
					exit(EXIT_FAILURE);
				}

				fgets(strLine, LONG_LINE_LENGTH, fpStat);
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				sscanf(strLine, "%s %d", strTemp, nN+1);

				fgets(strLine, LONG_LINE_LENGTH, fpStat);
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				sscanf(strLine, "%s %d", strTemp, nN+0);
				fclose(fpStat);
				
				dRatio = (double)(nN[0])/((double)(nN[1])+1e-6);
			}
			dLower = dLR;
		}

		/* get desired cutoff */
		dLR = (dLower+dUpper)/2.0;
		while( fabs(dUpper-dLower) > 1.0 )
		{
			MotifMap_ScanMatrix_Genome_Main(strMotifFile, strGenomePath, 
						strCodPath, strTempOutFile, dLR, 
						nBGOrder, strBGType, 
						strBGPath, nBGStepSize,
						nUseCS, dC, strCSPath, nIncludeRepeat);
			
			sprintf(strTempStatFile, "%s.stat", strTempOutFile);
			fpStat = NULL;
			fpStat = fopen(strTempStatFile, "r");
			if(fpStat == NULL)
			{
				printf("Error: MotifMap_ScanMatrix_Genome_GetCutoff_Main, cannot open temporary output file!\n");
				exit(EXIT_FAILURE);
			}

			fgets(strLine, LONG_LINE_LENGTH, fpStat);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			sscanf(strLine, "%s %d", strTemp, nN+1);

			fgets(strLine, LONG_LINE_LENGTH, fpStat);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			sscanf(strLine, "%s %d", strTemp, nN+0);
			fclose(fpStat);
			
			dRatio = (double)(nN[0])/((double)(nN[1])+1e-6);

			if(fabs(dRatio-dRate) < (0.01*dRate))
			{
				break;
			}

			if(dRatio>dRate)
			{
				dLower = dLR;
			}
			else
			{
				dUpper = dLR;
			}

			dLR = (dLower+dUpper)/2.0;
		}
	
		fprintf(fpOut, "%s\t%d\n", strMotifFile, (int)dLR);

		RemoveFiles(strTempOutFile);
		RemoveFiles(strTempStatFile);
	
		nMotifNum++;
	}

	/* close files */
	fclose(fpIn);
	fclose(fpOut);

	
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMatrix_Genome_Summary_Main: compute motif enrichment      */
/*  for targeted genomic regions. All motifs in the motif list will be     */
/*  examined.                                                              */
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanMatrix_Genome_Summary_Main(char strMotifListPath[], char strGenomePath[], 
					char strCodPath[], char strNegCodPath[], char strOutputPath[],
					int nBGOrder, char strBGType[], char strBGPath[], int nBGStepSize,
					int nUseCS, double dC, char strCSPath[], int nIncludeRepeat)
{
	/* define */
	FILE *fpIn = NULL;
	FILE *fpStat = NULL;
	FILE *fpOut = NULL;
	FILE *fpOutNR = NULL;
	double dR;
	int nMotifNum = 0;
	char strMotifFile[MED_LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	char strNROutputPath[MED_LINE_LENGTH];

	char strTempOutFile[MED_LINE_LENGTH];
	char strTempStatFile[MED_LINE_LENGTH];
	char strTemp[LINE_LENGTH];
	int ni,nj;
	
	int nN[6];
	int nNx;
	int nC[6];
	int nCx;
	double dD[6];
	double dE[6];
	double dRatio[6];

	/* process one by one */
	fpIn = NULL;
	fpIn = fopen(strMotifListPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: MotifMap_ScanMatrix_Genome_Summary_Main, cannot open motiflist file!\n");
		exit(EXIT_FAILURE);
	}
	fpOut = NULL;
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: MotifMap_ScanMatrix_Genome_Summary_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	sprintf(strNROutputPath, "%s.nr", strOutputPath);
	fpOutNR = NULL;
	fpOutNR = fopen(strNROutputPath, "w");
	if(fpOutNR == NULL)
	{
		printf("Error: MotifMap_ScanMatrix_Genome_Summary_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "#motif_id\tmotif_file\tn1B\tn2B\tn1C\tn2C\tr1\tn3B\tn4B\tn3C\tn4C\tr2\tr3\n");

	fprintf(fpOutNR, "#motif_id\tmotif_file\tn1B\tn2B\tn1C\tn2C\tr1\tn3B\tn4B\tn3C\tn4C\tr2\tr3\n");

	nMotifNum = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		sscanf(strLine, "%s %lf", strMotifFile, &dR);
		printf("Processing motif %d, %s ...\n", nMotifNum, strMotifFile);
		for(ni=0; ni<6; ni++)
		{
			nN[ni] = 0;
			nC[ni] = 0;
		}
		for(ni=0; ni<6; ni++)
		{
			dD[ni] = 0.0;
			dE[ni] = 0.0;
			dRatio[ni] = 0.0;
		}
		
		/* map nonconserved */
		sprintf(strTempOutFile, "%s.%d_target.map", strOutputPath, nMotifNum);
		MotifMap_ScanMatrix_Genome_Main(strMotifFile, strGenomePath, 
					strCodPath, strTempOutFile, dR, 
					nBGOrder, strBGType, 
					strBGPath, nBGStepSize,
					0, 0.0, "", nIncludeRepeat);
		
		sprintf(strNROutputPath, "%s.%d_target.mapnr", strOutputPath, nMotifNum);
		nN[4] = MotifMap_FilterOverlappingSite(strTempOutFile, strNROutputPath);

		sprintf(strTempStatFile, "%s.stat", strTempOutFile);
		fpStat = NULL;
		fpStat = fopen(strTempStatFile, "r");
		if(fpStat == NULL)
		{
			printf("Error: MotifMap_ScanMatrix_Genome_Summary_Main, cannot open temporary output file!\n");
			exit(EXIT_FAILURE);
		}

		fgets(strLine, LONG_LINE_LENGTH, fpStat);
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		sscanf(strLine, "%s %d", strTemp, nN+1);

		fgets(strLine, LONG_LINE_LENGTH, fpStat);
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		sscanf(strLine, "%s %d", strTemp, nN+0);
		fclose(fpStat);
		

		sprintf(strTempOutFile, "%s.%d_control.map", strOutputPath, nMotifNum);
		MotifMap_ScanMatrix_Genome_Main(strMotifFile, strGenomePath, 
					strNegCodPath, strTempOutFile, dR, 
					nBGOrder, strBGType, 
					strBGPath, nBGStepSize,
					0, 0.0, "", nIncludeRepeat);
	
		sprintf(strNROutputPath, "%s.%d_control.mapnr", strOutputPath, nMotifNum);
		nC[4] = MotifMap_FilterOverlappingSite(strTempOutFile, strNROutputPath);


		sprintf(strTempStatFile, "%s.stat", strTempOutFile);
		fpStat = NULL;
		fpStat = fopen(strTempStatFile, "r");
		if(fpStat == NULL)
		{
			printf("Error: MotifMap_ScanMatrix_Genome_Summary_Main, cannot open temporary output file!\n");
			exit(EXIT_FAILURE);
		}

		fgets(strLine, LONG_LINE_LENGTH, fpStat);
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		sscanf(strLine, "%s %d", strTemp, nC+1);

		fgets(strLine, LONG_LINE_LENGTH, fpStat);
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		sscanf(strLine, "%s %d", strTemp, nC+0);
		fclose(fpStat);


		if(nUseCS == 1)
		{
			sprintf(strTempOutFile, "%s.%d_target_c%d.map", strOutputPath, nMotifNum, (int)dC);
			MotifMap_ScanMatrix_Genome_Main(strMotifFile, strGenomePath, 
					strCodPath, strTempOutFile, dR, 
					nBGOrder, strBGType, 
					strBGPath, nBGStepSize,
					nUseCS, dC, strCSPath, nIncludeRepeat);

			sprintf(strNROutputPath, "%s.%d_target_c%d.mapnr", strOutputPath, nMotifNum, (int)dC);
			nN[5] = MotifMap_FilterOverlappingSite(strTempOutFile, strNROutputPath);

			sprintf(strTempStatFile, "%s.stat", strTempOutFile);
			fpStat = NULL;
			fpStat = fopen(strTempStatFile, "r");
			if(fpStat == NULL)
			{
				printf("Error: MotifMap_ScanMatrix_Genome_Summary_Main, cannot open temporary output file!\n");
				exit(EXIT_FAILURE);
			}
			fgets(strLine, LONG_LINE_LENGTH, fpStat);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			sscanf(strLine, "%s %d", strTemp, &nNx);
			if(nNx != nN[1])
			{
				printf("%d=%d\n", nNx, nN[1]);
				printf("Error: MotifMap_ScanMatrix_Genome_Enrich_Main, wrong count!\n");
				exit(EXIT_FAILURE);
			}

			fgets(strLine, LONG_LINE_LENGTH, fpStat);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			sscanf(strLine, "%s %d", strTemp, nN+2);

			fgets(strLine, LONG_LINE_LENGTH, fpStat);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			sscanf(strLine, "%s %d", strTemp, nN+3);
			fclose(fpStat);

			sprintf(strTempOutFile, "%s.%d_control_c%d.map", strOutputPath, nMotifNum, (int)dC);
			MotifMap_ScanMatrix_Genome_Main(strMotifFile, strGenomePath, 
					strNegCodPath, strTempOutFile, dR, 
					nBGOrder, strBGType, 
					strBGPath, nBGStepSize,
					nUseCS, dC, strCSPath, nIncludeRepeat);

			sprintf(strNROutputPath, "%s.%d_control_c%d.mapnr", strOutputPath, nMotifNum, (int)dC);
			nC[5] = MotifMap_FilterOverlappingSite(strTempOutFile, strNROutputPath);

			sprintf(strTempStatFile, "%s.stat", strTempOutFile);
			fpStat = NULL;
			fpStat = fopen(strTempStatFile, "r");
			if(fpStat == NULL)
			{
				printf("Error: MotifMap_ScanMatrix_Genome_Summary_Main, cannot open temporary output file!\n");
				exit(EXIT_FAILURE);
			}
			fgets(strLine, LONG_LINE_LENGTH, fpStat);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			sscanf(strLine, "%s %d", strTemp, &nCx);
			if(nCx != nC[1])
			{
				printf("%d=%d\n", nCx, nC[1]);
				printf("Error: MotifMap_ScanMatrix_Genome_Enrich_Main, wrong count!\n");
				exit(EXIT_FAILURE);
			}

			fgets(strLine, LONG_LINE_LENGTH, fpStat);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			sscanf(strLine, "%s %d", strTemp, nC+2);

			fgets(strLine, LONG_LINE_LENGTH, fpStat);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			sscanf(strLine, "%s %d", strTemp, nC+3);
			fclose(fpStat);
		}
		

		dD[0] = (double)(nN[0])/((double)(nN[1])+1e-6);
		dE[0] = (double)(nC[0])/((double)(nC[1])+1e-6);
		dD[3] = (double)(nN[4])/((double)(nN[1])+1e-6);
		dE[3] = (double)(nC[4])/((double)(nC[1])+1e-6);

		if(dE[0] > 0.0)
		{
			dRatio[0] = dD[0]/dE[0];
		}
		else
		{
			dRatio[0] = dD[0]/(dE[0]+1e-6);
		}
		if(dE[3] > 0.0)
		{
			dRatio[3] = dD[3]/dE[3];
		}
		else
		{
			dRatio[3] = dD[3]/(dE[3]+1e-6);
		}

		if(nUseCS == 1)
		{
			dD[1] = (double)(nN[2])/((double)(nN[3])+1e-6);
			dE[1] = (double)(nC[2])/((double)(nC[3])+1e-6);

			dD[2] = (double)(nN[2])/((double)(nN[1])+1e-6);
			dE[2] = (double)(nC[2])/((double)(nC[1])+1e-6);

			dD[4] = (double)(nN[5])/((double)(nN[3])+1e-6);
			dE[4] = (double)(nC[5])/((double)(nC[3])+1e-6);

			dD[5] = (double)(nN[5])/((double)(nN[1])+1e-6);
			dE[5] = (double)(nC[5])/((double)(nC[1])+1e-6);


			for(nj=1; nj<3; nj++)
			{
				if(dE[nj] > 0.0)
				{
					dRatio[nj] = dD[nj]/dE[nj];
				}
				else
				{
					dRatio[nj] = dD[nj]/(dE[nj]+1e-6);
				}
			}

			for(nj=4; nj<6; nj++)
			{
				if(dE[nj] > 0.0)
				{
					dRatio[nj] = dD[nj]/dE[nj];
				}
				else
				{
					dRatio[nj] = dD[nj]/(dE[nj]+1e-6);
				}
			}
		}
	
	
		fprintf(fpOut, "%d\t%s\t%d\t%d\t%d\t%d\t%f\t%d\t%d\t%d\t%d\t%f\t%f\n",
			nMotifNum, strMotifFile, nN[0], nN[1], nC[0], nC[1], dRatio[0],
			nN[2], nN[3], nC[2], nC[3], dRatio[1], dRatio[2]);

		fprintf(fpOutNR, "%d\t%s\t%d\t%d\t%d\t%d\t%f\t%d\t%d\t%d\t%d\t%f\t%f\n",
			nMotifNum, strMotifFile, nN[4], nN[1], nC[4], nC[1], dRatio[3],
			nN[5], nN[3], nC[5], nC[3], dRatio[4], dRatio[5]);
	
		nMotifNum++;
	}

	/* close files */
	fclose(fpIn);
	fclose(fpOut);
	fclose(fpOutNR);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanConsensus_Genome_Summary_Main: compute motif enrichment   */
/*  for targeted genomic regions. All motifs in the motif list will be     */
/*  examined.                                                              */
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanConsensus_Genome_Summary_Main(char strMotifListPath[], char strGenomePath[], 
					char strCodPath[], char strNegCodPath[], char strOutputPath[],
					int nUseCS, double dC, char strCSPath[], int nIncludeRepeat)
{
	/* define */
	FILE *fpIn = NULL;
	FILE *fpStat = NULL;
	FILE *fpOut = NULL;
	FILE *fpOutNR = NULL;
	int nMC,nMD;
	int nMotifNum = 0;
	char strMotifFile[MED_LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	char strNROutputPath[MED_LINE_LENGTH];

	char strTempOutFile[MED_LINE_LENGTH];
	char strTempStatFile[MED_LINE_LENGTH];
	char strTemp[LINE_LENGTH];
	int ni,nj;
	
	int nN[6];
	int nNx;
	int nC[6];
	int nCx;
	double dD[6];
	double dE[6];
	double dRatio[6];

	/* process one by one */
	fpIn = NULL;
	fpIn = fopen(strMotifListPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: MotifMap_ScanConsensus_Genome_Summary_Main, cannot open motiflist file!\n");
		exit(EXIT_FAILURE);
	}
	fpOut = NULL;
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: MotifMap_ScanConsensus_Genome_Summary_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	sprintf(strNROutputPath, "%s.nr", strOutputPath);
	fpOutNR = NULL;
	fpOutNR = fopen(strNROutputPath, "w");
	if(fpOutNR == NULL)
	{
		printf("Error: MotifMap_ScanConsensus_Genome_Summary_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "#motif_id\tmotif_file\tn1B\tn2B\tn1C\tn2C\tr1\tn3B\tn4B\tn3C\tn4C\tr2\tr3\n");

	fprintf(fpOutNR, "#motif_id\tmotif_file\tn1B\tn2B\tn1C\tn2C\tr1\tn3B\tn4B\tn3C\tn4C\tr2\tr3\n");

	nMotifNum = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		sscanf(strLine, "%s %d %d", strMotifFile, &nMC, &nMD);
		printf("Processing motif %d, %s ...\n", nMotifNum, strMotifFile);
		for(ni=0; ni<6; ni++)
		{
			nN[ni] = 0;
			nC[ni] = 0;
		}
		for(ni=0; ni<6; ni++)
		{
			dD[ni] = 0.0;
			dE[ni] = 0.0;
			dRatio[ni] = 0.0;
		}
		
		/* map nonconserved */
		sprintf(strTempOutFile, "%s.%d_target.map", strOutputPath, nMotifNum);
		MotifMap_ScanConsensus_Genome_Main(strMotifFile, strGenomePath, 
					strCodPath, strTempOutFile, nMC, nMD, 
					0, 0.0, "", nIncludeRepeat);
		
		sprintf(strNROutputPath, "%s.%d_target.mapnr", strOutputPath, nMotifNum);
		nN[4] = MotifMap_FilterOverlappingSite(strTempOutFile, strNROutputPath);

		sprintf(strTempStatFile, "%s.stat", strTempOutFile);
		fpStat = NULL;
		fpStat = fopen(strTempStatFile, "r");
		if(fpStat == NULL)
		{
			printf("Error: MotifMap_ScanConsensus_Genome_Summary_Main, cannot open temporary output file!\n");
			exit(EXIT_FAILURE);
		}

		fgets(strLine, LONG_LINE_LENGTH, fpStat);
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		sscanf(strLine, "%s %d", strTemp, nN+1);

		fgets(strLine, LONG_LINE_LENGTH, fpStat);
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		sscanf(strLine, "%s %d", strTemp, nN+0);
		fclose(fpStat);
		

		sprintf(strTempOutFile, "%s.%d_control.map", strOutputPath, nMotifNum);
		MotifMap_ScanConsensus_Genome_Main(strMotifFile, strGenomePath, 
					strNegCodPath, strTempOutFile, nMC, nMD, 
					0, 0.0, "", nIncludeRepeat);
	
		sprintf(strNROutputPath, "%s.%d_control.mapnr", strOutputPath, nMotifNum);
		nC[4] = MotifMap_FilterOverlappingSite(strTempOutFile, strNROutputPath);


		sprintf(strTempStatFile, "%s.stat", strTempOutFile);
		fpStat = NULL;
		fpStat = fopen(strTempStatFile, "r");
		if(fpStat == NULL)
		{
			printf("Error: MotifMap_ScanConsensus_Genome_Summary_Main, cannot open temporary output file!\n");
			exit(EXIT_FAILURE);
		}

		fgets(strLine, LONG_LINE_LENGTH, fpStat);
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		sscanf(strLine, "%s %d", strTemp, nC+1);

		fgets(strLine, LONG_LINE_LENGTH, fpStat);
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		sscanf(strLine, "%s %d", strTemp, nC+0);
		fclose(fpStat);


		if(nUseCS == 1)
		{
			sprintf(strTempOutFile, "%s.%d_target_c%d.map", strOutputPath, nMotifNum, (int)dC);
			MotifMap_ScanConsensus_Genome_Main(strMotifFile, strGenomePath, 
					strCodPath, strTempOutFile, nMC, nMD, 
					nUseCS, dC, strCSPath, nIncludeRepeat);

			sprintf(strNROutputPath, "%s.%d_target_c%d.mapnr", strOutputPath, nMotifNum, (int)dC);
			nN[5] = MotifMap_FilterOverlappingSite(strTempOutFile, strNROutputPath);

			sprintf(strTempStatFile, "%s.stat", strTempOutFile);
			fpStat = NULL;
			fpStat = fopen(strTempStatFile, "r");
			if(fpStat == NULL)
			{
				printf("Error: MotifMap_ScanConsensus_Genome_Summary_Main, cannot open temporary output file!\n");
				exit(EXIT_FAILURE);
			}
			fgets(strLine, LONG_LINE_LENGTH, fpStat);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			sscanf(strLine, "%s %d", strTemp, &nNx);
			if(nNx != nN[1])
			{
				printf("%d=%d\n", nNx, nN[1]);
				printf("Error: MotifMap_ScanConsensus_Genome_Enrich_Main, wrong count!\n");
				exit(EXIT_FAILURE);
			}

			fgets(strLine, LONG_LINE_LENGTH, fpStat);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			sscanf(strLine, "%s %d", strTemp, nN+2);

			fgets(strLine, LONG_LINE_LENGTH, fpStat);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			sscanf(strLine, "%s %d", strTemp, nN+3);
			fclose(fpStat);

			sprintf(strTempOutFile, "%s.%d_control_c%d.map", strOutputPath, nMotifNum, (int)dC);
			MotifMap_ScanConsensus_Genome_Main(strMotifFile, strGenomePath, 
					strNegCodPath, strTempOutFile, nMC, nMD, 
					nUseCS, dC, strCSPath, nIncludeRepeat);

			sprintf(strNROutputPath, "%s.%d_control_c%d.mapnr", strOutputPath, nMotifNum, (int)dC);
			nC[5] = MotifMap_FilterOverlappingSite(strTempOutFile, strNROutputPath);

			sprintf(strTempStatFile, "%s.stat", strTempOutFile);
			fpStat = NULL;
			fpStat = fopen(strTempStatFile, "r");
			if(fpStat == NULL)
			{
				printf("Error: MotifMap_ScanConsensus_Genome_Summary_Main, cannot open temporary output file!\n");
				exit(EXIT_FAILURE);
			}
			fgets(strLine, LONG_LINE_LENGTH, fpStat);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			sscanf(strLine, "%s %d", strTemp, &nCx);
			if(nCx != nC[1])
			{
				printf("%d=%d\n", nCx, nC[1]);
				printf("Error: MotifMap_ScanConsensus_Genome_Enrich_Main, wrong count!\n");
				exit(EXIT_FAILURE);
			}

			fgets(strLine, LONG_LINE_LENGTH, fpStat);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			sscanf(strLine, "%s %d", strTemp, nC+2);

			fgets(strLine, LONG_LINE_LENGTH, fpStat);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			sscanf(strLine, "%s %d", strTemp, nC+3);
			fclose(fpStat);
		}
		

		dD[0] = (double)(nN[0])/((double)(nN[1])+1e-6);
		dE[0] = (double)(nC[0])/((double)(nC[1])+1e-6);
		dD[3] = (double)(nN[4])/((double)(nN[1])+1e-6);
		dE[3] = (double)(nC[4])/((double)(nC[1])+1e-6);

		if(dE[0] > 0.0)
		{
			dRatio[0] = dD[0]/dE[0];
		}
		else
		{
			dRatio[0] = dD[0]/(dE[0]+1e-6);
		}
		if(dE[3] > 0.0)
		{
			dRatio[3] = dD[3]/dE[3];
		}
		else
		{
			dRatio[3] = dD[3]/(dE[3]+1e-6);
		}

		if(nUseCS == 1)
		{
			dD[1] = (double)(nN[2])/((double)(nN[3])+1e-6);
			dE[1] = (double)(nC[2])/((double)(nC[3])+1e-6);

			dD[2] = (double)(nN[2])/((double)(nN[1])+1e-6);
			dE[2] = (double)(nC[2])/((double)(nC[1])+1e-6);

			dD[4] = (double)(nN[5])/((double)(nN[3])+1e-6);
			dE[4] = (double)(nC[5])/((double)(nC[3])+1e-6);

			dD[5] = (double)(nN[5])/((double)(nN[1])+1e-6);
			dE[5] = (double)(nC[5])/((double)(nC[1])+1e-6);


			for(nj=1; nj<3; nj++)
			{
				if(dE[nj] > 0.0)
				{
					dRatio[nj] = dD[nj]/dE[nj];
				}
				else
				{
					dRatio[nj] = dD[nj]/(dE[nj]+1e-6);
				}
			}

			for(nj=4; nj<6; nj++)
			{
				if(dE[nj] > 0.0)
				{
					dRatio[nj] = dD[nj]/dE[nj];
				}
				else
				{
					dRatio[nj] = dD[nj]/(dE[nj]+1e-6);
				}
			}
		}
	
	
		fprintf(fpOut, "%d\t%s\t%d\t%d\t%d\t%d\t%f\t%d\t%d\t%d\t%d\t%f\t%f\n",
			nMotifNum, strMotifFile, nN[0], nN[1], nC[0], nC[1], dRatio[0],
			nN[2], nN[3], nC[2], nC[3], dRatio[1], dRatio[2]);

		fprintf(fpOutNR, "%d\t%s\t%d\t%d\t%d\t%d\t%f\t%d\t%d\t%d\t%d\t%f\t%f\n",
			nMotifNum, strMotifFile, nN[4], nN[1], nC[4], nC[1], dRatio[3],
			nN[5], nN[3], nC[5], nC[3], dRatio[4], dRatio[5]);
	
		nMotifNum++;
	}

	/* close files */
	fclose(fpIn);
	fclose(fpOut);
	fclose(fpOutNR);

	/* return */
	return PROC_SUCCESS;
}

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
					int nUseCS, double dC, char strCSPath[], int nIncludeRepeat)
{
	/* define */
	int nRegionNum = 0;
	int nTierNum = 0;
	int nResidual = 0;
	FILE *fpIn = NULL;
	FILE *fpOut = NULL;
	char strLine[LONG_LINE_LENGTH];
	struct INTMATRIX *pN = NULL;
	struct DOUBLEMATRIX *pD = NULL;
	struct DOUBLEMATRIX *pR = NULL;
	char strTempFile[MED_LINE_LENGTH];
	char strTempOutFile[MED_LINE_LENGTH];
	char strTempStatFile[MED_LINE_LENGTH];
	char strTemp[LINE_LENGTH];
	int ni,nj;
	int nN1,nN2,nN3,nN4,nNx;
	double dRatio;
	double dD;
	int nStart, nEnd;

	/* init check */
	if(nTierSize < 1)
	{
		printf("Error: MotifMap_ScanMatrix_Genome_Enrich_Main, tier size must be >=1!\n");
		exit(EXIT_FAILURE);
	}
	
	/* get region number */
	fpIn = NULL;
	fpIn = fopen(strCodPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: MotifMap_ScanMatrix_Genome_Enrich_Main, cannot open target coordinate file!\n");
		exit(EXIT_FAILURE);
	}

	nRegionNum = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		nRegionNum++;
	}

	fclose(fpIn);

	if(nRegionNum == 0)
	{
		printf("Warning: MotifMap_ScanMatrix_Genome_Enrich_Main, there are no regions for mapping!\n");
		return PROC_SUCCESS;
	}

	/* allocate memory */
	nTierNum = nRegionNum/nTierSize;
	nResidual = nRegionNum%nTierSize;
	if(nResidual > 0)
		nTierNum = nTierNum+1;

	pN = NULL;
	pN = CreateIntMatrix((nTierNum+1), 4);
	if(pN == NULL)
	{
		printf("Error: MotifMap_ScanMatrix_Genome_Enrich_Main, cannot create pN!\n");
		exit(EXIT_FAILURE);
	}
	pD = NULL;
	pD = CreateDoubleMatrix((nTierNum+1), 3);
	if(pD == NULL)
	{
		printf("Error: MotifMap_ScanMatrix_Genome_Enrich_Main, cannot create pD!\n");
		exit(EXIT_FAILURE);
	}
	pR = NULL;
	pR = CreateDoubleMatrix((nTierNum+1), 3);
	if(pR == NULL)
	{
		printf("Error: MotifMap_ScanMatrix_Genome_Enrich_Main, cannot create pR!\n");
		exit(EXIT_FAILURE);
	}

	/* process tier by tier */
	fpIn = NULL;
	fpIn = fopen(strCodPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: MotifMap_ScanMatrix_Genome_Enrich_Main, cannot open target coordinate file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	nj = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		nj = ni/nTierSize;
		if(ni%nTierSize == 0)
		{
			if(fpOut != NULL)
			{
				fclose(fpOut);
			}

			sprintf(strTempFile, "%s.tier%d", strCodPath, nj);
			fpOut = NULL;
			fpOut = fopen(strTempFile, "w");
			if(fpOut == NULL)
			{
				printf("Error: MotifMap_ScanMatrix_Genome_Enrich_Main, cannot open temporary tier coordinate file!\n");
				exit(EXIT_FAILURE);
			}
		}

		fprintf(fpOut, "%s\n", strLine);

		ni++;
	}

	if(fpOut != NULL)
	{
		fclose(fpOut);
	}

	fclose(fpIn);

	if(ni != nRegionNum)
	{
		printf("Error: MotifMap_ScanMatrix_Genome_Enrich_Main, region numbers do not match!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nTierNum; ni++)
	{
		/* if(ni == 32)
			ni = 32; */

		printf("processing tier %d ...\n", ni);
		sprintf(strTempFile, "%s.tier%d", strCodPath, ni);
		sprintf(strTempOutFile, "%s.tmp", strOutputPath);
		/* map nonconserved */
		MotifMap_ScanMatrix_Genome_Main(strMotifPath, strGenomePath, 
					strTempFile, strTempOutFile, dR, 
					nBGOrder, strBGType, 
					strBGPath, nBGStepSize,
					0, 0.0, "", nIncludeRepeat);
	
		sprintf(strTempStatFile, "%s.stat", strTempOutFile);
		fpIn = NULL;
		fpIn = fopen(strTempStatFile, "r");
		if(fpIn == NULL)
		{
			printf("Error: MotifMap_ScanMatrix_Genome_Enrich_Main, cannot open temporary output file!\n");
			exit(EXIT_FAILURE);
		}

		fgets(strLine, LONG_LINE_LENGTH, fpIn);
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		sscanf(strLine, "%s %d", strTemp, &nN2);

		fgets(strLine, LONG_LINE_LENGTH, fpIn);
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		sscanf(strLine, "%s %d", strTemp, &nN1);
		fclose(fpIn);

		IMSETAT(pN, ni, 0, nN1);
		IMSETAT(pN, ni, 1, nN2);

		if(nUseCS == 1)
		{
			MotifMap_ScanMatrix_Genome_Main(strMotifPath, strGenomePath, 
					strTempFile, strTempOutFile, dR, 
					nBGOrder, strBGType, 
					strBGPath, nBGStepSize,
					nUseCS, dC, strCSPath, nIncludeRepeat);

			sprintf(strTempStatFile, "%s.stat", strTempOutFile);
			fpIn = NULL;
			fpIn = fopen(strTempStatFile, "r");
			if(fpIn == NULL)
			{
				printf("Error: MotifMap_ScanMatrix_Genome_Enrich_Main, cannot open temporary output file!\n");
				exit(EXIT_FAILURE);
			}
			fgets(strLine, LONG_LINE_LENGTH, fpIn);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			sscanf(strLine, "%s %d", strTemp, &nNx);
			if(nNx != nN2)
			{
				printf("%d=%d\n", nNx, nN2);
				printf("Error: MotifMap_ScanMatrix_Genome_Enrich_Main, wrong count!\n");
				exit(EXIT_FAILURE);
			}

			fgets(strLine, LONG_LINE_LENGTH, fpIn);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			sscanf(strLine, "%s %d", strTemp, &nN3);

			fgets(strLine, LONG_LINE_LENGTH, fpIn);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			sscanf(strLine, "%s %d", strTemp, &nN4);
			fclose(fpIn);

			IMSETAT(pN, ni, 2, nN3);
			IMSETAT(pN, ni, 3, nN4);
		}

		RemoveFiles(strTempFile);
	}

	/* process control regions */
	printf("processing control regions ...\n", ni);
	sprintf(strTempOutFile, "%s.tmp", strOutputPath);
	/* map nonconserved */
	MotifMap_ScanMatrix_Genome_Main(strMotifPath, strGenomePath, 
				strNegCodPath, strTempOutFile, dR, 
				nBGOrder, strBGType, 
				strBGPath, nBGStepSize,
				0, 0.0, "", nIncludeRepeat);

	sprintf(strTempStatFile, "%s.stat", strTempOutFile);
	fpIn = NULL;
	fpIn = fopen(strTempStatFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: MotifMap_ScanMatrix_Genome_Enrich_Main, cannot open temporary output file!\n");
		exit(EXIT_FAILURE);
	}

	fgets(strLine, LONG_LINE_LENGTH, fpIn);
	StrTrimLeft(strLine);
	StrTrimRight(strLine);
	sscanf(strLine, "%s %d", strTemp, &nN2);

	fgets(strLine, LONG_LINE_LENGTH, fpIn);
	StrTrimLeft(strLine);
	StrTrimRight(strLine);
	sscanf(strLine, "%s %d", strTemp, &nN1);
	fclose(fpIn);

	IMSETAT(pN, ni, 0, nN1);
	IMSETAT(pN, ni, 1, nN2);

	if(nUseCS == 1)
	{
		MotifMap_ScanMatrix_Genome_Main(strMotifPath, strGenomePath, 
				strNegCodPath, strTempOutFile, dR, 
				nBGOrder, strBGType, 
				strBGPath, nBGStepSize,
				nUseCS, dC, strCSPath, nIncludeRepeat);

		sprintf(strTempStatFile, "%s.stat", strTempOutFile);
		fpIn = NULL;
		fpIn = fopen(strTempStatFile, "r");
		if(fpIn == NULL)
		{
			printf("Error: MotifMap_ScanMatrix_Genome_Enrich_Main, cannot open temporary output file!\n");
			exit(EXIT_FAILURE);
		}
		fgets(strLine, LONG_LINE_LENGTH, fpIn);
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		sscanf(strLine, "%s %d", strTemp, &nNx);
		if(nNx != nN2)
		{
			printf("%d=%d\n", nNx, nN2);
			printf("Error: MotifMap_ScanMatrix_Genome_Enrich_Main, wrong count!\n");
			exit(EXIT_FAILURE);
		}

		fgets(strLine, LONG_LINE_LENGTH, fpIn);
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		sscanf(strLine, "%s %d", strTemp, &nN3);

		fgets(strLine, LONG_LINE_LENGTH, fpIn);
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		sscanf(strLine, "%s %d", strTemp, &nN4);
		fclose(fpIn);

		IMSETAT(pN, ni, 2, nN3);
		IMSETAT(pN, ni, 3, nN4);
	}

	
	/* compute enrichment level */
	for(ni=0; ni<=nTierNum; ni++)
	{
		dRatio = (double)(IMGETAT(pN, ni, 0))/((double)(IMGETAT(pN, ni, 1))+1e-6);
		DMSETAT(pD, ni, 0, dRatio);
		if(nUseCS == 1)
		{
			dRatio = (double)(IMGETAT(pN, ni, 2))/((double)(IMGETAT(pN, ni, 3))+1e-6);
			DMSETAT(pD, ni, 1, dRatio);

			dRatio = (double)(IMGETAT(pN, ni, 2))/((double)(IMGETAT(pN, ni, 1))+1e-6);
			DMSETAT(pD, ni, 2, dRatio);
		}
	}

	for(nj=0; nj<1; nj++)
	{
		dD = DMGETAT(pD, (pR->nHeight-1), nj);
		if(dD > 0.0)
		{
			for(ni=0; ni<pR->nHeight; ni++)
			{
				dRatio = DMGETAT(pD, ni, nj)/dD;
				DMSETAT(pR, ni, nj, dRatio);
			}
		}
		else
		{
			for(ni=0; ni<pR->nHeight; ni++)
			{
				dRatio = DMGETAT(pD, ni, nj)/(dD+1e-6);
				DMSETAT(pR, ni, nj, dRatio);
			}
		}
	}

	if(nUseCS == 1)
	{
		for(nj=1; nj<pR->nWidth; nj++)
		{
			dD = DMGETAT(pD, (pR->nHeight-1), nj);
			if(dD > 0.0)
			{
				for(ni=0; ni<pR->nHeight; ni++)
				{
					dRatio = DMGETAT(pD, ni, nj)/dD;
					DMSETAT(pR, ni, nj, dRatio);
				}
			}
			else
			{
				for(ni=0; ni<pR->nHeight; ni++)
				{
					dRatio = DMGETAT(pD, ni, nj)/(dD+1e-6);
					DMSETAT(pR, ni, nj, dRatio);
				}
			}
		}
	}

	/* export */
	fpOut = NULL;
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: MotifMap_ScanMatrix_Genome_Enrich_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "#peak_rank\tn1\tn2\tr1\tn3\tn4\tr2\tr3\n");

	for(ni=0; ni<nTierNum; ni++)
	{
		nStart = ni*nTierSize+1;
		nEnd = nStart+nTierSize-1;
		if(nEnd > nRegionNum)
			nEnd = nRegionNum;

		fprintf(fpOut, "%d~%d\t%d\t%d\t%f\t%d\t%d\t%f\t%f\n",
			nStart, nEnd, IMGETAT(pN, ni, 0), IMGETAT(pN, ni, 1), DMGETAT(pR, ni, 0), 
			IMGETAT(pN, ni, 2), IMGETAT(pN, ni, 3), DMGETAT(pR, ni, 1), DMGETAT(pR, ni, 2));
	}

	fprintf(fpOut, "control\t%d\t%d\t%f\t%d\t%d\t%f\t%f\n",
			IMGETAT(pN, ni, 0), IMGETAT(pN, ni, 1), DMGETAT(pR, ni, 0), 
			IMGETAT(pN, ni, 2), IMGETAT(pN, ni, 3), DMGETAT(pR, ni, 1), DMGETAT(pR, ni, 2));

	fclose(fpOut);

	/* release memory */
	DestroyIntMatrix(pN);
	DestroyDoubleMatrix(pD);
	DestroyDoubleMatrix(pR);

	/* remove files */
	sprintf(strTempOutFile, "%s.tmp", strOutputPath);
	RemoveFiles(strTempOutFile);
	sprintf(strTempStatFile, "%s.stat", strTempOutFile);
	RemoveFiles(strTempStatFile);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanConsensus_Genome_Enrich_Main: compute motif enrichment    */
/*  for targeted and ordered genomic regions. Target regions will be       */
/*  grouped into several tiers according to nTierSize. The enrichment      */
/*  level of each tier will be compared to negative control regions.       */
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanConsensus_Genome_Enrich_Main(char strMotifPath[], char strGenomePath[], 
					char strCodPath[], char strNegCodPath[], int nTierSize, 
					char strOutputPath[], int nMC, int nMD, 
					int nUseCS, double dC, char strCSPath[], int nIncludeRepeat)
{
	/* define */
	int nRegionNum = 0;
	int nTierNum = 0;
	int nResidual = 0;
	FILE *fpIn = NULL;
	FILE *fpOut = NULL;
	char strLine[LONG_LINE_LENGTH];
	struct INTMATRIX *pN = NULL;
	struct DOUBLEMATRIX *pD = NULL;
	struct DOUBLEMATRIX *pR = NULL;
	char strTempFile[MED_LINE_LENGTH];
	char strTempOutFile[MED_LINE_LENGTH];
	char strTempStatFile[MED_LINE_LENGTH];
	char strTemp[LINE_LENGTH];
	int ni,nj;
	int nN1,nN2,nN3,nN4,nNx;
	double dRatio;
	double dD;
	int nStart, nEnd;

	/* init check */
	if(nTierSize < 1)
	{
		printf("Error: MotifMap_ScanConsensus_Genome_Enrich_Main, tier size must be >=1!\n");
		exit(EXIT_FAILURE);
	}
	
	/* get region number */
	fpIn = NULL;
	fpIn = fopen(strCodPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: MotifMap_ScanConsensus_Genome_Enrich_Main, cannot open target coordinate file!\n");
		exit(EXIT_FAILURE);
	}

	nRegionNum = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		nRegionNum++;
	}

	fclose(fpIn);

	if(nRegionNum == 0)
	{
		printf("Warning: MotifMap_ScanConsensus_Genome_Enrich_Main, there are no regions for mapping!\n");
		return PROC_SUCCESS;
	}

	/* allocate memory */
	nTierNum = nRegionNum/nTierSize;
	nResidual = nRegionNum%nTierSize;
	if(nResidual > 0)
		nTierNum = nTierNum+1;

	pN = NULL;
	pN = CreateIntMatrix((nTierNum+1), 4);
	if(pN == NULL)
	{
		printf("Error: MotifMap_ScanConsensus_Genome_Enrich_Main, cannot create pN!\n");
		exit(EXIT_FAILURE);
	}
	pD = NULL;
	pD = CreateDoubleMatrix((nTierNum+1), 3);
	if(pD == NULL)
	{
		printf("Error: MotifMap_ScanConsensus_Genome_Enrich_Main, cannot create pD!\n");
		exit(EXIT_FAILURE);
	}
	pR = NULL;
	pR = CreateDoubleMatrix((nTierNum+1), 3);
	if(pR == NULL)
	{
		printf("Error: MotifMap_ScanConsensus_Genome_Enrich_Main, cannot create pR!\n");
		exit(EXIT_FAILURE);
	}

	/* process tier by tier */
	fpIn = NULL;
	fpIn = fopen(strCodPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: MotifMap_ScanConsensus_Genome_Enrich_Main, cannot open target coordinate file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	nj = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		nj = ni/nTierSize;
		if(ni%nTierSize == 0)
		{
			if(fpOut != NULL)
			{
				fclose(fpOut);
			}

			sprintf(strTempFile, "%s.tier%d", strCodPath, nj);
			fpOut = NULL;
			fpOut = fopen(strTempFile, "w");
			if(fpOut == NULL)
			{
				printf("Error: MotifMap_ScanConsensus_Genome_Enrich_Main, cannot open temporary tier coordinate file!\n");
				exit(EXIT_FAILURE);
			}
		}

		fprintf(fpOut, "%s\n", strLine);

		ni++;
	}

	if(fpOut != NULL)
	{
		fclose(fpOut);
	}

	fclose(fpIn);

	if(ni != nRegionNum)
	{
		printf("Error: MotifMap_ScanConsensus_Genome_Enrich_Main, region numbers do not match!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nTierNum; ni++)
	{
		printf("processing tier %d ...\n", ni);
		sprintf(strTempFile, "%s.tier%d", strCodPath, ni);
		sprintf(strTempOutFile, "%s.tmp", strOutputPath);
		/* map nonconserved */
		MotifMap_ScanConsensus_Genome_Main(strMotifPath, strGenomePath, 
					strTempFile, strTempOutFile, nMC, nMD, 
					0, 0.0, "", nIncludeRepeat);
	
		sprintf(strTempStatFile, "%s.stat", strTempOutFile);
		fpIn = NULL;
		fpIn = fopen(strTempStatFile, "r");
		if(fpIn == NULL)
		{
			printf("Error: MotifMap_ScanConsensus_Genome_Enrich_Main, cannot open temporary output file!\n");
			exit(EXIT_FAILURE);
		}

		fgets(strLine, LONG_LINE_LENGTH, fpIn);
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		sscanf(strLine, "%s %d", strTemp, &nN2);

		fgets(strLine, LONG_LINE_LENGTH, fpIn);
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		sscanf(strLine, "%s %d", strTemp, &nN1);
		fclose(fpIn);

		IMSETAT(pN, ni, 0, nN1);
		IMSETAT(pN, ni, 1, nN2);

		if(nUseCS == 1)
		{
			MotifMap_ScanConsensus_Genome_Main(strMotifPath, strGenomePath, 
				strTempFile, strTempOutFile, nMC, nMD, 
				nUseCS, dC, strCSPath, nIncludeRepeat);

			sprintf(strTempStatFile, "%s.stat", strTempOutFile);
			fpIn = NULL;
			fpIn = fopen(strTempStatFile, "r");
			if(fpIn == NULL)
			{
				printf("Error: MotifMap_ScanConsensus_Genome_Enrich_Main, cannot open temporary output file!\n");
				exit(EXIT_FAILURE);
			}
			fgets(strLine, LONG_LINE_LENGTH, fpIn);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			sscanf(strLine, "%s %d", strTemp, &nNx);
			if(nNx != nN2)
			{
				printf("%d=%d\n", nNx, nN2);
				printf("Error: MotifMap_ScanConsensus_Genome_Enrich_Main, wrong count!\n");
				exit(EXIT_FAILURE);
			}

			fgets(strLine, LONG_LINE_LENGTH, fpIn);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			sscanf(strLine, "%s %d", strTemp, &nN3);

			fgets(strLine, LONG_LINE_LENGTH, fpIn);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			sscanf(strLine, "%s %d", strTemp, &nN4);
			fclose(fpIn);

			IMSETAT(pN, ni, 2, nN3);
			IMSETAT(pN, ni, 3, nN4);
		}

		RemoveFiles(strTempFile);
	}

	/* process control regions */
	printf("processing control regions ...\n", ni);
	sprintf(strTempOutFile, "%s.tmp", strOutputPath);
	/* map nonconserved */
	MotifMap_ScanConsensus_Genome_Main(strMotifPath, strGenomePath, 
			strNegCodPath, strTempOutFile, nMC, nMD, 
			0, 0.0, "", nIncludeRepeat);

	sprintf(strTempStatFile, "%s.stat", strTempOutFile);
	fpIn = NULL;
	fpIn = fopen(strTempStatFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: MotifMap_ScanConsensus_Genome_Enrich_Main, cannot open temporary output file!\n");
		exit(EXIT_FAILURE);
	}

	fgets(strLine, LONG_LINE_LENGTH, fpIn);
	StrTrimLeft(strLine);
	StrTrimRight(strLine);
	sscanf(strLine, "%s %d", strTemp, &nN2);

	fgets(strLine, LONG_LINE_LENGTH, fpIn);
	StrTrimLeft(strLine);
	StrTrimRight(strLine);
	sscanf(strLine, "%s %d", strTemp, &nN1);
	fclose(fpIn);

	IMSETAT(pN, ni, 0, nN1);
	IMSETAT(pN, ni, 1, nN2);

	if(nUseCS == 1)
	{
		MotifMap_ScanConsensus_Genome_Main(strMotifPath, strGenomePath, 
				strNegCodPath, strTempOutFile, nMC, nMD, 
				nUseCS, dC, strCSPath, nIncludeRepeat);

		sprintf(strTempStatFile, "%s.stat", strTempOutFile);
		fpIn = NULL;
		fpIn = fopen(strTempStatFile, "r");
		if(fpIn == NULL)
		{
			printf("Error: MotifMap_ScanConsensus_Genome_Enrich_Main, cannot open temporary output file!\n");
			exit(EXIT_FAILURE);
		}
		fgets(strLine, LONG_LINE_LENGTH, fpIn);
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		sscanf(strLine, "%s %d", strTemp, &nNx);
		if(nNx != nN2)
		{
			printf("%d=%d\n", nNx, nN2);
			printf("Error: MotifMap_ScanConsensus_Genome_Enrich_Main, wrong count!\n");
			exit(EXIT_FAILURE);
		}

		fgets(strLine, LONG_LINE_LENGTH, fpIn);
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		sscanf(strLine, "%s %d", strTemp, &nN3);

		fgets(strLine, LONG_LINE_LENGTH, fpIn);
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		sscanf(strLine, "%s %d", strTemp, &nN4);
		fclose(fpIn);

		IMSETAT(pN, ni, 2, nN3);
		IMSETAT(pN, ni, 3, nN4);
	}

	
	/* compute enrichment level */
	for(ni=0; ni<=nTierNum; ni++)
	{
		dRatio = (double)(IMGETAT(pN, ni, 0))/((double)(IMGETAT(pN, ni, 1))+1e-6);
		DMSETAT(pD, ni, 0, dRatio);
		if(nUseCS == 1)
		{
			dRatio = (double)(IMGETAT(pN, ni, 2))/((double)(IMGETAT(pN, ni, 3))+1e-6);
			DMSETAT(pD, ni, 1, dRatio);

			dRatio = (double)(IMGETAT(pN, ni, 2))/((double)(IMGETAT(pN, ni, 1))+1e-6);
			DMSETAT(pD, ni, 2, dRatio);
		}
	}

	for(nj=0; nj<1; nj++)
	{
		dD = DMGETAT(pD, (pR->nHeight-1), nj);
		if(dD > 0.0)
		{
			for(ni=0; ni<pR->nHeight; ni++)
			{
				dRatio = DMGETAT(pD, ni, nj)/dD;
				DMSETAT(pR, ni, nj, dRatio);
			}
		}
		else
		{
			for(ni=0; ni<pR->nHeight; ni++)
			{
				dRatio = DMGETAT(pD, ni, nj)/(dD+1e-6);
				DMSETAT(pR, ni, nj, dRatio);
			}
		}
	}

	if(nUseCS == 1)
	{
		for(nj=1; nj<pR->nWidth; nj++)
		{
			dD = DMGETAT(pD, (pR->nHeight-1), nj);
			if(dD > 0.0)
			{
				for(ni=0; ni<pR->nHeight; ni++)
				{
					dRatio = DMGETAT(pD, ni, nj)/dD;
					DMSETAT(pR, ni, nj, dRatio);
				}
			}
			else
			{
				for(ni=0; ni<pR->nHeight; ni++)
				{
					dRatio = DMGETAT(pD, ni, nj)/(dD+1e-6);
					DMSETAT(pR, ni, nj, dRatio);
				}
			}
		}
	}

	/* export */
	fpOut = NULL;
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: MotifMap_ScanConsensus_Genome_Enrich_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "#peak_rank\tn1\tn2\tr1\tn3\tn4\tr2\tr3\n");

	for(ni=0; ni<nTierNum; ni++)
	{
		nStart = ni*nTierSize+1;
		nEnd = nStart+nTierSize-1;
		if(nEnd > nRegionNum)
			nEnd = nRegionNum;

		fprintf(fpOut, "%d~%d\t%d\t%d\t%f\t%d\t%d\t%f\t%f\n",
			nStart, nEnd, IMGETAT(pN, ni, 0), IMGETAT(pN, ni, 1), DMGETAT(pR, ni, 0), 
			IMGETAT(pN, ni, 2), IMGETAT(pN, ni, 3), DMGETAT(pR, ni, 1), DMGETAT(pR, ni, 2));
	}

	fprintf(fpOut, "control\t%d\t%d\t%f\t%d\t%d\t%f\t%f\n",
			IMGETAT(pN, ni, 0), IMGETAT(pN, ni, 1), DMGETAT(pR, ni, 0), 
			IMGETAT(pN, ni, 2), IMGETAT(pN, ni, 3), DMGETAT(pR, ni, 1), DMGETAT(pR, ni, 2));

	fclose(fpOut);

	/* release memory */
	DestroyIntMatrix(pN);
	DestroyDoubleMatrix(pD);
	DestroyDoubleMatrix(pR);

	/* remove files */
	sprintf(strTempOutFile, "%s.tmp", strOutputPath);
	RemoveFiles(strTempOutFile);
	sprintf(strTempStatFile, "%s.stat", strTempOutFile);
	RemoveFiles(strTempStatFile);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_GetSiteAround_Main: get regions surrounding sites derived     */
/*  from motif mapping                                                     */
/* ----------------------------------------------------------------------- */ 
int MotifMap_GetSiteAround_Main(char strInputPath[], char strGenomePath[], 
			char strOutputPath[], char strSpecies[], int nW, int nCN, int nA)
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	struct INTMATRIX *pChrLen;
	int nChr,nStart,nEnd;
	char chStrand;
	char strChr[LINE_LENGTH];
	char strAlias[LINE_LENGTH];
	double dScore;
	char strSite[LINE_LENGTH];

	char strLine[MED_LINE_LENGTH];
	char strFileName[MED_LINE_LENGTH];
	int nlen;
	int ni;

	/* init */
	StrMakeUpper(strSpecies);
	if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
	{
		nlen = strlen(strGenomePath);
		if(nlen == 0)
		{
			sprintf(strGenomePath, ".\\"); 
		}
		else if(strGenomePath[nlen-1] != '\\')
		{
			strGenomePath[nlen] = '\\';
			strGenomePath[nlen+1] = '\0';
		}
	}
	else
	{
		nlen = strlen(strGenomePath);
		if(nlen == 0)
		{
			sprintf(strGenomePath, "./"); 
		}
		else if(strGenomePath[nlen-1] != '/')
		{
			strGenomePath[nlen] = '/';
			strGenomePath[nlen+1] = '\0';
		}
	}

	/* get chromosome length */
	sprintf(strFileName, "%schrlen.txt", strGenomePath);
	pChrLen = NULL;
	pChrLen = IMLOAD(strFileName);
	if(pChrLen == NULL)
	{
		printf("Error: MotifMap_GetSiteAround_Main, cannot load chromosome length!\n");
		exit(EXIT_FAILURE);
	}

	/* load file */
	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: MotifMap_GetSiteAround_Main, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: MotifMap_GetSiteAround_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		sscanf(strLine, "%s %s %d %d %c %lf %s", strAlias, strChr, &nStart, &nEnd,
			&chStrand, &dScore, strSite);

		nChr = Genome_ChromosomeName_To_Index(strChr, strSpecies);
		if(nChr <= 0)
			continue;

		nStart -= nW;
		if(nStart < 0)
			nStart = 0;
		nEnd += nW;
		if(nEnd >= pChrLen->pMatElement[nChr-1])
			nEnd = pChrLen->pMatElement[nChr-1]-1;

		if(nA == 0)
		{
			if(nCN == 1)
			{
				fprintf(fpOut, "%s\t%d\t%d\t%d\t%c\t%f\t%s\n", strAlias, nChr, nStart, nEnd,
					chStrand, dScore, strSite);
			}
			else
			{
				fprintf(fpOut, "%s\t%s\t%d\t%d\t%c\t%f\t%s\n", strAlias, strChr, nStart, nEnd,
					chStrand, dScore, strSite);
			}
		}
		else
		{
			if(nCN == 1)
			{
				fprintf(fpOut, "%d\t%d\t%d\t%d\t%c\t%f\t%s\n", ni, nChr, nStart, nEnd,
					chStrand, dScore, strSite);
			}
			else
			{
				fprintf(fpOut, "%d\t%s\t%d\t%d\t%c\t%f\t%s\n", ni, strChr, nStart, nEnd,
					chStrand, dScore, strSite);
			}
		}

		ni++;
	}

	fclose(fpIn);
	fclose(fpOut);

	/* release memory */
	DestroyIntMatrix(pChrLen);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_getcluster_Main: get clustered sites.                         */
/*  If distance between 2 sites <= nW, they will be report as clustered.   */
/*  nInputType is the input file type, 0=cod, 1=bed, 2=codp, 3=bedp.       */
/* ----------------------------------------------------------------------- */ 
int MotifMap_getcluster_Main(char strInputPath[], char strOutputPath[],
							 int nW, int nInputType)
{
	/* define */
	FILE *fpIn = NULL;
	FILE *fpOut = NULL;
	FILE *fpClust = NULL;
	char strLine[LONG_LINE_LENGTH];
	char strChr[MED_LINE_LENGTH];
	char strSeqName[MED_LINE_LENGTH];
	int nStart,nEnd,nPos;
	int nClusterNum = 0;
	int nP1,nP2;
	char chStrand;
	char strLastChr[MED_LINE_LENGTH];
	char strLastLine[LONG_LINE_LENGTH];
	int nLastLineWritten;
	int nLastPos = -1;
	int nGoodModule = 0;

	/* open files */
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: MotifMap_getcluster_Main, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: MotifMap_getcluster_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	sprintf(strLine, "%s.cluster.cod", strOutputPath);
	fpClust = fopen(strLine, "w");
	if(fpClust == NULL)
	{
		printf("Error: MotifMap_getcluster_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	
	/* load */
	strcpy(strLastChr, "NULL");
	nLastLineWritten = 0;
	nP1 = 0;
	nP2 = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		if(nInputType == 1)
		{
			sscanf(strLine, "%s %d %d", strChr, &nStart, &nEnd);
			nPos = (nStart+nEnd)/2;
		}
		else if(nInputType == 2)
		{
			sscanf(strLine, "%s %s %d", strSeqName, strChr, &nPos);
			nStart = nPos;
			nEnd = nPos;
		}
		else if(nInputType == 3)
		{
			sscanf(strLine, "%s %d", strChr, &nPos);
			nStart = nPos;
			nEnd = nPos;
		}
		else
		{
			sscanf(strLine, "%s %s %d %d %c", strSeqName, strChr, &nStart, &nEnd, &chStrand);
			nPos = (nStart+nEnd)/2;
		}

		if(strcmp(strChr, strLastChr) != 0)
		{
			if(nGoodModule == 1)
			{
				nClusterNum++;
				fprintf(fpClust, "%d\t%s\t%d\t%d\t+\n", nClusterNum, strLastChr, nP1, nP2);
			}

			strcpy(strLastLine, strLine);
			strcpy(strLastChr, strChr);
			nLastPos = nPos;
			nLastLineWritten = 0;
			nP1 = nStart;
			nP2 = nEnd;
			nGoodModule = 0;
		}
		else
		{
			if( (nPos-nLastPos) <= nW )
			{
				if(nLastLineWritten == 0)
				{
					fprintf(fpOut, "%s\n", strLastLine);
				}
				fprintf(fpOut, "%s\n", strLine);
				strcpy(strLastLine, strLine);
				nLastPos = nPos;
				nLastLineWritten = 1;
				if(nP2 < nEnd)
					nP2 = nEnd;
				nGoodModule = 1;
			}
			else
			{
				if(nGoodModule == 1)
				{
					nClusterNum++;
					fprintf(fpClust, "%d\t%s\t%d\t%d\t+\n", nClusterNum, strLastChr, nP1, nP2);
				}

				strcpy(strLastLine, strLine);
				nLastPos = nPos;
				nLastLineWritten = 0;
				nP1 = nStart;
				nP2 = nEnd;
				nGoodModule = 0;
			}
		}
	}

	if(nGoodModule == 1)
	{
		nClusterNum++;
		fprintf(fpClust, "%d\t%s\t%d\t%d\t+\n", nClusterNum, strLastChr, nP1, nP2);
	}


	/* close file */
	fclose(fpIn);
	fclose(fpOut);
	fclose(fpClust);


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_GetSiteAroundCS_Genome_Main: get mean conservation score      */
/*  TFBS sites.                                                            */
/* ----------------------------------------------------------------------- */ 
int MotifMap_GetSiteAroundCS_Genome_Main(char strInputPath[], char strOutputPath[], 
		int nMotifLen, int nW, char strSpecies[], char strGenomePath[], char strCSPath[])
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	FILE *fpScore;
	struct DOUBLEMATRIX *pStat;
	struct INTMATRIX *pChrLen;
	char strLine[LONG_LINE_LENGTH];
	
	char strAlias[LINE_LENGTH];
	char strChr[LINE_LENGTH];
	int nChr;
	int nStart,nEnd;
	char chStrand;
	int nSiteNum;
	int nLen;
	
	struct BYTEMATRIX *pScore;
	int ni;
	char strFileName[MED_LINE_LENGTH];
	int numread;

	/* init */
	if( nW < 0 )
	{
		printf("Error: MotifMap_GetSiteAroundCS_Genome_Main, window size should be >=0!\n");
		exit(EXIT_FAILURE);
	}

	if(nMotifLen < 0)
	{
		fpIn = NULL;
		fpIn = fopen(strInputPath, "r");
		if(fpIn == NULL)
		{
			printf("Error: MotifMap_GetSiteAroundCS_Genome_Main, cannot open input file!\n");
			exit(EXIT_FAILURE);
		}

		/* process line by line */
		while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;
			if(strLine[0] == '#')
				continue;

			sscanf(strLine, "%s %s %d %d %c", strAlias, strChr, &nStart, &nEnd, &chStrand);
			nMotifLen = nEnd-nStart+1;
			break;
		}

		fclose(fpIn);
	}

	if( nMotifLen <= 0 )
	{
		printf("Error: MotifMap_GetSiteAroundCS_Genome_Main, motif length should be >0!\n");
		exit(EXIT_FAILURE);
	}

	AdjustDirectoryPath(strGenomePath);
	AdjustDirectoryPath(strCSPath);
	sprintf(strFileName, "%schrlen.txt", strGenomePath);
	pChrLen = NULL;
	pChrLen = IMLOAD(strFileName);
	if(pChrLen == NULL)
	{
		printf("Error: MotifMap_GetSiteAroundCS_Genome_Main, cannot open chromosome length!\n");
		exit(EXIT_FAILURE);
	}

	nLen = nMotifLen+2*nW;
	pStat = NULL;
	pStat = CreateDoubleMatrix(nLen,1);
	if(pStat == NULL)
	{
		printf("Error: MotifMap_GetSiteAroundCS_Genome_Main, cannot create statistic matrix!\n");
		exit(EXIT_FAILURE);
	}

	pScore = NULL;
	pScore = CreateByteMatrix(nLen,1);
	if(pScore == NULL)
	{
		printf("Error: MotifMap_GetSiteAroundCS_Genome_Main, cannot create statistic matrix!\n");
		exit(EXIT_FAILURE);
	}

	/* open files */
	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: MotifMap_GetSiteAroundCS_Genome_Main, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}

	/* process line by line */
	nSiteNum = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		sscanf(strLine, "%s %s %d %d %c", strAlias, strChr, &nStart, &nEnd, &chStrand);
		if(nMotifLen != (nEnd-nStart+1) )
		{
			printf("Warning: MotifMap_GetSiteAroundCS_Genome_Main, motif length not consistent. %s\n", strLine);
			continue;
		}

		nChr = Genome_ChromosomeName_To_Index(strChr, strSpecies);
		nChr = nChr-1;
		if( (nChr<0) || (nChr>=pChrLen->nHeight) )
		{
			printf("Warning: MotifMap_GetSiteAroundCS_Genome_Main, chromosome out of range. %s\n", strLine);
			continue;
		}

		nStart = nStart-nW;
		nEnd = nEnd+nW;
		if( (nStart<0) || (nEnd >= pChrLen->pMatElement[nChr]) )
		{
			printf("Warning: MotifMap_GetSiteAroundCS_Genome_Main, chromosome out of range. %s\n", strLine);
			continue;
		}

		/* load */
		sprintf(strFileName, "%s%s.cs", strCSPath, strChr);
		fpScore = NULL;
		fpScore = fopen(strFileName, "rb");
		if(fpScore == NULL)
		{
			printf("Warning: MotifMap_GetSiteAroundCS_Genome_Main, cannot open conservation score file %s!\n", strFileName);
			continue;
		}
		else
		{
			/* load score */
			if( fseek( fpScore, nStart, SEEK_SET ) != 0 )
			{
				printf("Error: MotifMap_Filter_Genome_Main, cannot locate the required sequence %s!\n", strLine);
				exit(EXIT_FAILURE);
			}

			/* read */
			numread = fread(pScore->pMatElement, sizeof(unsigned char), nLen, fpScore);
			if(numread != nLen)
			{
				printf("Error: MotifMap_Filter_Genome_Main, loading error %s!\n", strLine);
				exit(EXIT_FAILURE);
			}

			fclose(fpScore);
		

			/* check */
			if(chStrand == '-')
			{
				for(ni=0; ni<nLen; ni++)
				{
					pStat->pMatElement[ni] += (double)(pScore->pMatElement[nLen-1-ni]);
				}
			}
			else
			{
				for(ni=0; ni<nLen; ni++)
				{
					pStat->pMatElement[ni] += (double)(pScore->pMatElement[ni]);
				}
			}

			nSiteNum++;
		}
	}

	/* close files */
	fclose(fpIn);

	if(nSiteNum > 0)
	{
		DMPDIVTS(pStat, (double)nSiteNum);
	}

	fpOut = NULL;
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: MotifMap_GetSiteAroundCS_Genome_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	/* fprintf(fpOut, "range 0 255\n"); */
	/* fprintf(fpOut, "colors 0x0000ff\n"); */
	fprintf(fpOut, "labels \"position\" \"conservation score\"\n");
	fprintf(fpOut, "data\n");

	for(ni=0; ni<nLen; ni++)
	{
		fprintf(fpOut, "%d\t%f\n", ni-nW, pStat->pMatElement[ni]);
	}

	fclose(fpOut);

	/* release memory */
	DestroyIntMatrix(pChrLen);
	DestroyDoubleMatrix(pStat);
	DestroyByteMatrix(pScore);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_Filter_Genome_Main: fiter TFBS by footprint conservation      */
/*  and mask TFBS in coding sequences.                                     */
/* ----------------------------------------------------------------------- */ 
int MotifMap_Filter_Genome_Main(char strInputPath[], char strOutputPath[], 
		int nUseCS, double dC, char strMaskPath[], char strCSPath[], 
		int nUseCds, double dCds, char strCdsPath[])
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	FILE *fpScore;
	struct INTMATRIX *pMask;
	char strLine[LONG_LINE_LENGTH];
	char strAlias[LINE_LENGTH];
	char strChr[LINE_LENGTH];
	int nStart,nEnd;
	char chStrand;
	int nMotifLen;
	double dCdsCount;
	int nMaskPass;
	int nOK;
	int nOriginSites,nFilterSites;
	struct BYTEMATRIX *pScore;
	int ni;
	char strFileName[MED_LINE_LENGTH];
	int numread;

	/* init */
	if(nUseCS == 1)
	{
		pMask = NULL;
		pMask = IMLOAD(strMaskPath);
		if(pMask == NULL)
		{
			printf("Error: MotifMap_Filter_Genome_Main, cannot open mask file!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* open files */
	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: MotifMap_Filter_Genome_Main, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}
	fpOut = NULL;
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: MotifMap_Filter_Genome_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	/* process line by line */
	nOriginSites = 0;
	nFilterSites = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
		{
			fprintf(fpOut, "%s\n", strLine);
			continue;
		}

		sscanf(strLine, "%s %s %d %d %c", strAlias, strChr, &nStart, &nEnd, &chStrand);
		nMotifLen = nEnd-nStart+1;
		if(nMotifLen <= 0)
		{
			printf("Error: MotifMap_Filter_Genome_Main, motif length should be > 0!\n%s\n", strLine);
			exit(EXIT_FAILURE);
		}

		nOriginSites++;
		nOK = 1;
		pScore = NULL;
		pScore = CreateByteMatrix(1, nMotifLen);
		if(pScore == NULL)
		{
			printf("Error: MotifMap_Filter_Genome_Main, cannot allocate memory to load scores!\n");
			exit(EXIT_FAILURE);
		}

		/* footprint filter */
		if( nUseCS == 1)
		{
			if( (nMotifLen != pMask->nHeight) && (nMotifLen != pMask->nWidth) )
			{
				printf("Error: MotifMap_Filter_Genome_Main, mask length and site length do not match!\n");
				exit(EXIT_FAILURE);
			}

			/* load */
			AdjustDirectoryPath(strCSPath);
			sprintf(strFileName, "%s%s.cs", strCSPath, strChr);
			fpScore = NULL;
			fpScore = fopen(strFileName, "rb");
			if(fpScore == NULL)
			{
				printf("Warning: MotifMap_Filter_Genome_Main, cannot open conservation score file %s!\n", strFileName);
				nOK = 0;
			}
			else
			{
				/* load score */
				if( fseek( fpScore, nStart, SEEK_SET ) != 0 )
				{
					printf("Error: MotifMap_Filter_Genome_Main, cannot locate the required sequence %s!\n", strLine);
					exit(EXIT_FAILURE);
				}

				/* read */
				numread = fread(pScore->pMatElement, sizeof(unsigned char), nMotifLen, fpScore);
				if(numread != nMotifLen)
				{
					printf("Error: MotifMap_Filter_Genome_Main, loading error %s!\n", strLine);
					exit(EXIT_FAILURE);
				}

				fclose(fpScore);
			

				/* check */
				nMaskPass = 1;
				if(chStrand == '-')
				{
					for(ni=0; ni<nMotifLen; ni++)
					{
						if(pMask->pMatElement[ni] > 0)
						{
							if(pScore->pMatElement[nMotifLen-1-ni] < dC)
							{
								nMaskPass = 0;
								break;
							}
						}
					}
				}
				else
				{
					for(ni=0; ni<nMotifLen; ni++)
					{
						if(pMask->pMatElement[ni] > 0)
						{
							if(pScore->pMatElement[ni] < dC)
							{
								nMaskPass = 0;
								break;
							}
						}
					}
				}

				if(nMaskPass == 0)
					nOK = 0;
			}
		}

		/* cds filter */
		if( (nOK == 1) && (nUseCds == 1) )
		{
			dCdsCount = 0.0;

			/* load */
			AdjustDirectoryPath(strCdsPath);
			sprintf(strFileName, "%s%s.cds", strCdsPath, strChr);
			fpScore = NULL;
			fpScore = fopen(strFileName, "rb");
			if(fpScore == NULL)
			{
				printf("Warning: MotifMap_Filter_Genome_Main, cannot open CDS score file %s!\n", strFileName);
				nOK = 0;
			}
			else
			{
				/* load score */
				if( fseek( fpScore, nStart, SEEK_SET ) != 0 )
				{
					printf("Error: MotifMap_Filter_Genome_Main, cannot locate the required sequence %s!\n", strLine);
					exit(EXIT_FAILURE);
				}

				/* read */
				numread = fread(pScore->pMatElement, sizeof(unsigned char), nMotifLen, fpScore);
				if(numread != nMotifLen)
				{
					printf("Error: MotifMap_Filter_Genome_Main, loading error %s!\n", strLine);
					exit(EXIT_FAILURE);
				}

				fclose(fpScore);

				/* check */
				for(ni=0; ni<nMotifLen; ni++)
				{
					if(pScore->pMatElement[ni] == 0)
					{
						dCdsCount += 1.0;
					}
				}

				dCdsCount /= (double)nMotifLen;
				if(dCdsCount < dCds)
					nOK = 0;
			}
		}

		/* destroy scores */
		DestroyByteMatrix(pScore);

		/* export filtered sites */
		if(nOK == 1)
		{
			fprintf(fpOut, "%s\n", strLine);
			nFilterSites++;
		}
	}

	/* close files */
	fclose(fpIn);
	fclose(fpOut);

	if(nUseCS == 1)
	{
		DestroyIntMatrix(pMask);
	}

	/* save statistics */
	sprintf(strFileName, "%s.stat", strOutputPath);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: MotifMap_Filter_Genome_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "OriginSites= %d\n", nOriginSites);
	fprintf(fpOut, "FilteredSites= %d\n", nFilterSites);

	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_SimuScoreDistribution_typeI: get motif mapping score          */
/*  distribution through simulation. This function will tell you if a      */
/*  sequence looks more like a motif than the background or not.           */ 
/* ----------------------------------------------------------------------- */ 
int MotifMap_SimuScoreDistribution_typeI(char strMotifName[], struct DOUBLEMATRIX *pMotifPWM, 
								   int nBGOrder, struct DOUBLEMATRIX *pBG0, 
								   struct DOUBLEMATRIX *pBGF, struct DOUBLEMATRIX *pBGR, 
								   int nSimuNum, struct DOUBLEMATRIX *pQ, 
								   char strCutoffPath[])
{
	/* define */
	struct DOUBLEMATRIX *pScore;
	struct DOUBLEMATRIX *pSortedScore;
	struct DOUBLEMATRIX *pStat;
	int ni,nj,nHeight;
	char strLine[LINE_LENGTH];
	double *pEle,*pEle2;
	double dTotal;
	struct DOUBLEMATRIX *pLogBGF;
	struct DOUBLEMATRIX *pLogBGR;
	struct DOUBLEMATRIX *pLogPWM;
	struct BYTEMATRIX *pSeq;
	int nLen;
	int nBurninRound;

	
	/* check parameter */
	if((pMotifPWM == NULL) || (pBG0 == NULL) || (pBGF == NULL) || (pBGR == NULL) )
	{
		printf("Warning: MotifMap_SimuScoreDistribution, null parameters!\n");
		return PROC_FAILURE;
	}
	nHeight = (int)pow((double)(pBGF->nWidth), (double)(nBGOrder));
	if( (nHeight != pBGF->nHeight) || (nHeight != pBGR->nHeight) )
	{
		printf("Warning: MotifMap_SimuScoreDistribution, dimension not match!\n");
		return PROC_FAILURE;
	}
	if(pQ == NULL)
	{
		printf("Warning: MotifMap_SimuScoreDistribution, no simulation was done since no quantile parameter!\n");
		return PROC_FAILURE;
	}
	if(nSimuNum <= 0)
	{
		printf("Warning: MotifMap_SimuScoreDistribution, no simulation was done since simunum=0!\n");
		return PROC_FAILURE;
	}

	/* init */
	nBurninRound = 10;
	nLen = pMotifPWM->nHeight+(2+nBurninRound)*nBGOrder;
	pSeq = NULL;
	pSeq = CreateByteMatrix(1, nLen);
	if(pSeq == NULL)
	{
		printf("Error: MotifMap_SimuScoreDistribution, cannot create memory for simulating sequence!\n");
		exit(EXIT_FAILURE);
	}

	pLogBGF = NULL;
	pLogBGF = DMCLONE(pBGF);
	pLogBGR = NULL;
	pLogBGR = DMCLONE(pBGR);
	pLogPWM = NULL;
	pLogPWM = DMCLONE(pMotifPWM);
	if( (pLogBGF == NULL) || (pLogBGR == NULL) || (pLogPWM == NULL) )
	{
		printf("Error: MotifMap_SimuScoreDistribution, cannot create log transformed matrix!\n");
		exit(EXIT_FAILURE);
	}

	pEle = pLogPWM->pMatElement;
	for(ni=0; ni<pLogPWM->nHeight; ni++)
	{
		pEle2 = pEle;
		dTotal = 0.0;
		for(nj=0; nj<pLogPWM->nWidth; nj++)
		{
			*pEle2 += 0.5;
			dTotal += (*pEle2);
			pEle2++;
		}
		for(nj=0; nj<pLogPWM->nWidth; nj++)
		{
			*pEle = log(*pEle/dTotal);
			pEle++;
		}
	}

	pEle = pLogBGF->pMatElement;
	for(ni=0; ni<pLogBGF->nHeight; ni++)
	{
		for(nj=0; nj<pLogBGF->nWidth; nj++)
		{
			*pEle = log(*pEle);
			pEle++;
		}
	}

	pEle = pLogBGR->pMatElement;
	for(ni=0; ni<pLogBGR->nHeight; ni++)
	{
		for(nj=0; nj<pLogBGR->nWidth; nj++)
		{
			*pEle = log(*pEle);
			pEle++;
		}
	}

	pScore = NULL;
	pScore = CreateDoubleMatrix(1, nSimuNum);
	if(pScore == NULL)
	{
		printf("Error: MotifMap_SimuScoreDistribution, cannot create score matrix!\n");
		exit(EXIT_FAILURE);
	}

	/* simulation */
	for(ni=0; ni<nSimuNum; ni++)
	{
		pScore->pMatElement[ni] = MotifMap_SimuScore_typeI(pSeq, nBurninRound,
			pLogPWM, nBGOrder, pBG0, pBGF, pBGR, pLogBGF, pLogBGR);
	}


	/* sort score */
	pSortedScore = NULL;
	DMSORTMERGEA_0(pScore, &pSortedScore, NULL);
	pStat = NULL;
	pStat = CreateDoubleMatrix(pQ->nHeight, pQ->nWidth);
	if(pStat == NULL)
	{
		printf("Error: MotifMap_SimuScoreDistribution, cannot get quantile!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<pQ->nWidth; ni++)
	{
		nj = (int)(nSimuNum*(1.0-pQ->pMatElement[ni]));
		if(nj < 0)
			nj = 0;
		if(nj >= nSimuNum)
			nj = nSimuNum-1;
		pStat->pMatElement[ni] = pSortedScore->pMatElement[nj];
	}

	/* output */
	sprintf(strLine, "%s%s.stat", strCutoffPath, strMotifName);
	DMSAVE(pStat, strLine);

	/* destroy */
	DestroyDoubleMatrix(pScore);
	DestroyDoubleMatrix(pSortedScore);
	DestroyDoubleMatrix(pStat);
	DestroyDoubleMatrix(pLogPWM);
	DestroyDoubleMatrix(pLogBGF);
	DestroyDoubleMatrix(pLogBGR);
	DestroyByteMatrix(pSeq);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_SimuScore_typeI: simulate a background sequence and calculate */
/*  its motif mapping score.                                               */
/* ----------------------------------------------------------------------- */ 
double MotifMap_SimuScore_typeI(struct BYTEMATRIX *pSeq, int nBurninRound, 
						  struct DOUBLEMATRIX *pLogPWM, int nBGOrder, 
						  struct DOUBLEMATRIX *pBG0, 
						  struct DOUBLEMATRIX *pBGF, struct DOUBLEMATRIX *pBGR,
						  struct DOUBLEMATRIX *pLogBGF, struct DOUBLEMATRIX *pLogBGR)
{
	/* define */
	double dScore;
	double dMF,dMR,dML;
	double dBF,dBR,dBL;
	int ni,nj,nk,nx;
	unsigned char *pBase;
	double dRand;
	int nMotifLen;
	int nTypeNum;

	int nWordId;
	int nScale;
	
	/* simulation */
	if(nBGOrder == 0)
	{
		pBase = pSeq->pMatElement;
		dMF = 0.0;
		dMR = 0.0;
		dBF = 0.0;
		nMotifLen = pLogPWM->nHeight;
		for(ni=0; ni<nMotifLen; ni++)
		{
			dRand = rand_u();
			for(nj=0; nj<pBG0->nWidth; nj++)
			{
				dRand -= pBG0->pMatElement[nj];
				if(dRand <= 0.0)
					break;
			}
			if(nj == pBG0->nWidth)
				nj--;

			pBase[ni] = (unsigned char)nj;
			dBF += pLogBGF->pMatElement[nj];
			dMF += DMGETAT(pLogPWM, ni, nj);
			dMR += DMGETAT(pLogPWM, (nMotifLen-1-ni), (3-nj));
		}

		if(dMF > dMR)
			dML = dMF;
		else
			dML = dMR;

		dScore = exp(dML-dBF);
	}
	else if(nBGOrder > 0)
	{
		nMotifLen = pLogPWM->nHeight;
		nTypeNum = pBGF->nWidth;
		nScale = (int)pow((double)(nTypeNum), (double)(nBGOrder-1));
		nWordId = 0;
		pBase = pSeq->pMatElement;
		dMF = 0.0;
		dMR = 0.0;
		dBF = 0.0;
		dBR = 0.0;

		/* first burn in */
		nk = 0;
		for(ni=0; ni<nBGOrder; ni++)
		{
			dRand = rand_u();
			for(nj=0; nj<pBG0->nWidth; nj++)
			{
				dRand -= pBG0->pMatElement[nj];
				if(dRand <= 0.0)
					break;
			}
			if(nj == pBG0->nWidth)
				nj--;

			pBase[nk] = (unsigned char)nj;
			nWordId = nWordId*nTypeNum+nj;
			nk++;
		}

		for(nx=0; nx<nBurninRound; nx++)
		{
			for(ni=0; ni<nBGOrder; ni++)
			{
				dRand = rand_u();
				for(nj=0; nj<nTypeNum; nj++)
				{
					dRand -= DMGETAT(pBGF, nWordId, nj);
					if(dRand <= 0.0)
						break;
				}
				if(nj == nTypeNum)
					nj--;

				pBase[nk] = (unsigned char)nj;
				nWordId -= pBase[nk-nBGOrder]*nScale;
				nWordId = nWordId*nTypeNum+nj;
				nk++;
			}
		}

		/* calculate */
		for(ni=0; ni<nBGOrder; ni++)
		{
			dRand = rand_u();
			for(nj=0; nj<nTypeNum; nj++)
			{
				dRand -= DMGETAT(pBGF, nWordId, nj);
				if(dRand <= 0.0)
					break;
			}
			if(nj == nTypeNum)
				nj--;

			pBase[nk] = (unsigned char)nj;
			dBF += DMGETAT(pLogBGF, nWordId, nj);
			dMF += DMGETAT(pLogPWM, ni, nj);
			dMR += DMGETAT(pLogPWM, (nMotifLen-1-ni), (3-nj));

			nWordId -= pBase[nk-nBGOrder]*nScale;
			nWordId = nWordId*nTypeNum+nj;
			nk++;
		}

		for(; ni<nMotifLen; ni++)
		{
			dRand = rand_u();
			for(nj=0; nj<nTypeNum; nj++)
			{
				dRand -= DMGETAT(pBGF, nWordId, nj);
				if(dRand <= 0.0)
					break;
			}
			if(nj == nTypeNum)
				nj--;

			pBase[nk] = (unsigned char)nj;
			dBF += DMGETAT(pLogBGF, nWordId, nj);
			dMF += DMGETAT(pLogPWM, ni, nj);
			dMR += DMGETAT(pLogPWM, (nMotifLen-1-ni), (3-nj));

			nWordId -= pBase[nk-nBGOrder]*nScale;
			nWordId = nWordId*nTypeNum+nj;

			dBR += DMGETAT(pLogBGR, nWordId, pBase[nk-nBGOrder]);

			nk++;
		}

		for(ni=0; ni<nBGOrder; ni++)
		{
			dRand = rand_u();
			for(nj=0; nj<nTypeNum; nj++)
			{
				dRand -= DMGETAT(pBGF, nWordId, nj);
				if(dRand <= 0.0)
					break;
			}
			if(nj == nTypeNum)
				nj--;

			pBase[nk] = (unsigned char)nj;
			nWordId -= pBase[nk-nBGOrder]*nScale;
			nWordId = nWordId*nTypeNum+nj;
			dBR += DMGETAT(pLogBGR, nWordId, pBase[nk-nBGOrder]);

			nk++;
		}

		if(dMF > dMR)
			dML = dMF;
		else
			dML = dMR;

		if(dBF > dBR)
			dBL = dBF;
		else
			dBL = dBR;

		dScore = exp(dML-dBL);

	}
	else
	{
		printf("Error: MotifMap_SimuScore, nBGOrder<0!\n");
	}

	/* return */
	return dScore;
}


/* ----------------------------------------------------------------------- */ 
/*  MotifMap_SimuScoreDistribution_typeII: get motif score distribution    */
/*  through simulation. Motifs will be simulated from PWM and quantiles    */
/*  for log(Scores) will be estimated. This function will give you a sense */
/*  of to what extent a sequence looks like a motif instead of telling     */
/*  if that sequence looks more like a motif than like the background.     */  
/* ----------------------------------------------------------------------- */ 
int MotifMap_SimuScoreDistribution_typeII(char strMotifName[], 
				struct DOUBLEMATRIX *pMotifPWM, int nSimuNum, 
				struct DOUBLEMATRIX *pQ, char strCutoffPath[])
{
	/* define */
	struct DOUBLEMATRIX *pScore;
	struct DOUBLEMATRIX *pSortedScore;
	struct DOUBLEMATRIX *pStat;
	int ni,nj;
	char strLine[LINE_LENGTH];
	double *pEle,*pEle1,*pEle2;
	double dTotal;
	struct DOUBLEMATRIX *pPWM;
	struct DOUBLEMATRIX *pLogPWM;
	struct BYTEMATRIX *pSeq;
	int nLen;
	
	
	/* check parameter */
	if(pMotifPWM == NULL)
	{
		printf("Warning: MotifMap_SimuScoreDistribution_typeII, null parameters!\n");
		return PROC_FAILURE;
	}
	if(pQ == NULL)
	{
		printf("Warning: MotifMap_SimuScoreDistribution_typeII, no simulation was done since no quantile parameter!\n");
		return PROC_FAILURE;
	}
	if(nSimuNum <= 0)
	{
		printf("Warning: MotifMap_SimuScoreDistribution_typeII, no simulation was done since simunum=0!\n");
		return PROC_FAILURE;
	}

	/* init */
	nLen = pMotifPWM->nHeight;
	pSeq = NULL;
	pSeq = CreateByteMatrix(1, nLen);
	if(pSeq == NULL)
	{
		printf("Error: MotifMap_SimuScoreDistribution_typeII, cannot create memory for simulating sequence!\n");
		exit(EXIT_FAILURE);
	}

	pPWM = NULL;
	pPWM = DMCLONE(pMotifPWM);
	if(pPWM == NULL)
	{
		printf("Error: MotifMap_SimuScoreDistribution_typeII, cannot create normalized matrix!\n");
		exit(EXIT_FAILURE);
	}

	pLogPWM = NULL;
	pLogPWM = DMCLONE(pMotifPWM);
	if(pLogPWM == NULL)
	{
		printf("Error: MotifMap_SimuScoreDistribution_typeII, cannot create log transformed matrix!\n");
		exit(EXIT_FAILURE);
	}

	pEle = pLogPWM->pMatElement;
	pEle1 = pPWM->pMatElement;
	for(ni=0; ni<pLogPWM->nHeight; ni++)
	{
		pEle2 = pEle;
		dTotal = 0.0;
		for(nj=0; nj<pLogPWM->nWidth; nj++)
		{
			*pEle2 += 0.5;
			dTotal += (*pEle2);
			pEle2++;
		}
		for(nj=0; nj<pLogPWM->nWidth; nj++)
		{
			*pEle1 = *pEle/dTotal;
			*pEle = log(*pEle/dTotal);
			pEle1++;
			pEle++;
		}
	}

	pScore = NULL;
	pScore = CreateDoubleMatrix(1, nSimuNum);
	if(pScore == NULL)
	{
		printf("Error: MotifMap_SimuScoreDistribution_typeII, cannot create score matrix!\n");
		exit(EXIT_FAILURE);
	}

	/* simulation */
	for(ni=0; ni<nSimuNum; ni++)
	{
		pScore->pMatElement[ni] = MotifMap_SimuScore_typeII(pSeq, pPWM, pLogPWM);
	}


	/* sort score */
	pSortedScore = NULL;
	DMSORTMERGEA_0(pScore, &pSortedScore, NULL);
	pStat = NULL;
	pStat = CreateDoubleMatrix(pQ->nHeight, pQ->nWidth);
	if(pStat == NULL)
	{
		printf("Error: MotifMap_SimuScoreDistribution_typeII, cannot get quantile!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<pQ->nWidth; ni++)
	{
		nj = (int)(nSimuNum*(pQ->pMatElement[ni])-1);
		if(nj < 0)
			nj = 0;
		if(nj >= nSimuNum)
			nj = nSimuNum-1;
		pStat->pMatElement[ni] = pSortedScore->pMatElement[nj];
	}

	/* output */
	sprintf(strLine, "%s%s.true", strCutoffPath, strMotifName);
	DMSAVE(pStat, strLine);

	/* destroy */
	DestroyDoubleMatrix(pScore);
	DestroyDoubleMatrix(pSortedScore);
	DestroyDoubleMatrix(pStat);
	DestroyDoubleMatrix(pPWM);
	DestroyDoubleMatrix(pLogPWM);
	DestroyByteMatrix(pSeq);

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  MotifMap_SimuScore_typeII: simulate a background sequence and          */
/*  calculate its motif mapping score.                                     */
/* ----------------------------------------------------------------------- */ 
double MotifMap_SimuScore_typeII(struct BYTEMATRIX *pSeq, struct DOUBLEMATRIX *pPWM, 
						  struct DOUBLEMATRIX *pLogPWM)
{
	/* define */
	double dScore;
	int ni,nj;
	unsigned char *pBase;
	double dRand;
	int nMotifLen;
	
	/* simulation */
	dScore = 0.0;
	nMotifLen = pPWM->nHeight;
	pBase = pSeq->pMatElement;
	for(ni=0; ni<nMotifLen; ni++)
	{
		dRand = rand_u();
		for(nj=0; nj<pPWM->nWidth; nj++)
		{
			dRand -= DMGETAT(pPWM, ni, nj);
			if(dRand <= 0.0)
				break;
		}
		if(nj == pPWM->nWidth)
			nj--;

		pBase[ni] = (unsigned char)nj;
		dScore += DMGETAT(pLogPWM, ni, nj);
	}

	/* return */
	return dScore;
}


/* ----------------------------------------------------------------------- */ 
/*  MotifMap_SimuScoreDistribution_typeIII: get motif mapping score        */
/*  distribution through simulation. This function will tell you the distn */
/*  of likelihood ratio under the motif model.                             */ 
/* ----------------------------------------------------------------------- */ 
int MotifMap_SimuScoreDistribution_typeIII(char strMotifName[], struct DOUBLEMATRIX *pMotifPWM, 
								   int nBGOrder, struct DOUBLEMATRIX *pBG0, 
								   struct DOUBLEMATRIX *pBGF, struct DOUBLEMATRIX *pBGR, 
								   int nSimuNum, struct DOUBLEMATRIX *pQ, 
								   char strCutoffPath[])
{
	/* define */
	struct DOUBLEMATRIX *pScore;
	struct DOUBLEMATRIX *pSortedScore;
	struct DOUBLEMATRIX *pStat;
	int ni,nj,nHeight;
	char strLine[LINE_LENGTH];
	double *pEle,*pEle2;
	double dTotal;
	struct DOUBLEMATRIX *pLogBGF;
	struct DOUBLEMATRIX *pLogBGR;
	struct DOUBLEMATRIX *pLogPWM;
	struct DOUBLEMATRIX *pPWM;
	struct BYTEMATRIX *pSeq;
	int nLen;
	int nBurninRound;

	
	/* check parameter */
	if((pMotifPWM == NULL) || (pBG0 == NULL) || (pBGF == NULL) || (pBGR == NULL) )
	{
		printf("Warning: MotifMap_SimuScoreDistribution, null parameters!\n");
		return PROC_FAILURE;
	}
	nHeight = (int)pow((double)(pBGF->nWidth), (double)(nBGOrder));
	if( (nHeight != pBGF->nHeight) || (nHeight != pBGR->nHeight) )
	{
		printf("Warning: MotifMap_SimuScoreDistribution, dimension not match!\n");
		return PROC_FAILURE;
	}
	if(pQ == NULL)
	{
		printf("Warning: MotifMap_SimuScoreDistribution, no simulation was done since no quantile parameter!\n");
		return PROC_FAILURE;
	}
	if(nSimuNum <= 0)
	{
		printf("Warning: MotifMap_SimuScoreDistribution, no simulation was done since simunum=0!\n");
		return PROC_FAILURE;
	}

	/* init */
	nBurninRound = 10;
	nLen = pMotifPWM->nHeight+(2+nBurninRound)*nBGOrder;
	pSeq = NULL;
	pSeq = CreateByteMatrix(1, nLen);
	if(pSeq == NULL)
	{
		printf("Error: MotifMap_SimuScoreDistribution, cannot create memory for simulating sequence!\n");
		exit(EXIT_FAILURE);
	}

	pLogBGF = NULL;
	pLogBGF = DMCLONE(pBGF);
	pLogBGR = NULL;
	pLogBGR = DMCLONE(pBGR);
	pLogPWM = NULL;
	pLogPWM = DMCLONE(pMotifPWM);
	pPWM = NULL;
	pPWM = DMCLONE(pMotifPWM);
	if( (pLogBGF == NULL) || (pLogBGR == NULL) || (pLogPWM == NULL) || (pPWM == NULL))
	{
		printf("Error: MotifMap_SimuScoreDistribution, cannot create log transformed matrix!\n");
		exit(EXIT_FAILURE);
	}

	pEle = pLogPWM->pMatElement;
	for(ni=0; ni<pLogPWM->nHeight; ni++)
	{
		pEle2 = pEle;
		dTotal = 0.0;
		for(nj=0; nj<pLogPWM->nWidth; nj++)
		{
			*pEle2 += 0.5;
			dTotal += (*pEle2);
			pEle2++;
		}
		for(nj=0; nj<pLogPWM->nWidth; nj++)
		{
			*pEle = log(*pEle/dTotal);
			pEle++;
		}
	}

	pEle = pPWM->pMatElement;
	for(ni=0; ni<pPWM->nHeight; ni++)
	{
		pEle2 = pEle;
		dTotal = 0.0;
		for(nj=0; nj<pPWM->nWidth; nj++)
		{
			*pEle2 += 0.5;
			dTotal += (*pEle2);
			pEle2++;
		}
		for(nj=0; nj<pPWM->nWidth; nj++)
		{
			*pEle = (*pEle/dTotal);
			pEle++;
		}
	}

	pEle = pLogBGF->pMatElement;
	for(ni=0; ni<pLogBGF->nHeight; ni++)
	{
		for(nj=0; nj<pLogBGF->nWidth; nj++)
		{
			*pEle = log(*pEle);
			pEle++;
		}
	}

	pEle = pLogBGR->pMatElement;
	for(ni=0; ni<pLogBGR->nHeight; ni++)
	{
		for(nj=0; nj<pLogBGR->nWidth; nj++)
		{
			*pEle = log(*pEle);
			pEle++;
		}
	}

	pScore = NULL;
	pScore = CreateDoubleMatrix(1, nSimuNum);
	if(pScore == NULL)
	{
		printf("Error: MotifMap_SimuScoreDistribution, cannot create score matrix!\n");
		exit(EXIT_FAILURE);
	}

	/* simulation */
	for(ni=0; ni<nSimuNum; ni++)
	{
		pScore->pMatElement[ni] = MotifMap_SimuScore_typeIII(pSeq, nBurninRound,
			pLogPWM, nBGOrder, pPWM, pBG0, pBGF, pBGR, pLogBGF, pLogBGR);
	}


	/* sort score */
	pSortedScore = NULL;
	DMSORTMERGEA_0(pScore, &pSortedScore, NULL);
	pStat = NULL;
	pStat = CreateDoubleMatrix(pQ->nHeight, pQ->nWidth);
	if(pStat == NULL)
	{
		printf("Error: MotifMap_SimuScoreDistribution, cannot get quantile!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<pQ->nWidth; ni++)
	{
		nj = (int)(nSimuNum*(pQ->pMatElement[ni])-1);
		if(nj < 0)
			nj = 0;
		if(nj >= nSimuNum)
			nj = nSimuNum-1;
		pStat->pMatElement[ni] = pSortedScore->pMatElement[nj];
	}

	/* output */
	sprintf(strLine, "%s%s.pows", strCutoffPath, strMotifName);
	DMSAVE(pStat, strLine);

	/* destroy */
	DestroyDoubleMatrix(pScore);
	DestroyDoubleMatrix(pSortedScore);
	DestroyDoubleMatrix(pStat);
	DestroyDoubleMatrix(pLogPWM);
	DestroyDoubleMatrix(pPWM);
	DestroyDoubleMatrix(pLogBGF);
	DestroyDoubleMatrix(pLogBGR);
	DestroyByteMatrix(pSeq);

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  MotifMap_SimuScore_typeIII: simulate a motif sequence and calculate    */
/*  its motif mapping score.                                               */
/* ----------------------------------------------------------------------- */ 
double MotifMap_SimuScore_typeIII(struct BYTEMATRIX *pSeq, int nBurninRound, 
						  struct DOUBLEMATRIX *pLogPWM, int nBGOrder, 
						  struct DOUBLEMATRIX *pPWM, struct DOUBLEMATRIX *pBG0,
						  struct DOUBLEMATRIX *pBGF, struct DOUBLEMATRIX *pBGR,
						  struct DOUBLEMATRIX *pLogBGF, struct DOUBLEMATRIX *pLogBGR)
{
	/* define */
	double dScore;
	double dMF,dMR,dML;
	double dBF,dBR,dBL;
	int ni,nj,nk,nx;
	unsigned char *pBase;
	double dRand;
	int nMotifLen;
	int nTypeNum;

	int nWordId;
	int nScale;
	
	/* simulation */
	if(nBGOrder == 0)
	{
		pBase = pSeq->pMatElement;
		dMF = 0.0;
		dMR = 0.0;
		dBF = 0.0;
		nMotifLen = pLogPWM->nHeight;
		for(ni=0; ni<nMotifLen; ni++)
		{
			dRand = rand_u();
			for(nj=0; nj<pPWM->nWidth; nj++)
			{
				dRand -= DMGETAT(pPWM, ni, nj);
				if(dRand <= 0.0)
					break;
			}
			if(nj == pPWM->nWidth)
				nj--;

			pBase[ni] = (unsigned char)nj;
			dBF += pLogBGF->pMatElement[nj];
			dMF += DMGETAT(pLogPWM, ni, nj);
			dMR += DMGETAT(pLogPWM, (nMotifLen-1-ni), (3-nj));
		}

		if(dMF > dMR)
			dML = dMF;
		else
			dML = dMR;

		dScore = exp(dML-dBF);
	}
	else if(nBGOrder > 0)
	{
		nMotifLen = pLogPWM->nHeight;
		nTypeNum = pBGF->nWidth;
		nScale = (int)pow((double)(nTypeNum), (double)(nBGOrder-1));
		nWordId = 0;
		pBase = pSeq->pMatElement;
		dMF = 0.0;
		dMR = 0.0;
		dBF = 0.0;
		dBR = 0.0;

		/* first burn in */
		nk = 0;
		for(ni=0; ni<nBGOrder; ni++)
		{
			dRand = rand_u();
			for(nj=0; nj<pBG0->nWidth; nj++)
			{
				dRand -= pBG0->pMatElement[nj];
				if(dRand <= 0.0)
					break;
			}
			if(nj == pBG0->nWidth)
				nj--;

			pBase[nk] = (unsigned char)nj;
			nWordId = nWordId*nTypeNum+nj;
			nk++;
		}

		for(nx=0; nx<nBurninRound; nx++)
		{
			for(ni=0; ni<nBGOrder; ni++)
			{
				dRand = rand_u();
				for(nj=0; nj<nTypeNum; nj++)
				{
					dRand -= DMGETAT(pBGF, nWordId, nj);
					if(dRand <= 0.0)
						break;
				}
				if(nj == nTypeNum)
					nj--;

				pBase[nk] = (unsigned char)nj;
				nWordId -= pBase[nk-nBGOrder]*nScale;
				nWordId = nWordId*nTypeNum+nj;
				nk++;
			}
		}

		/* calculate */
		for(ni=0; ni<nBGOrder; ni++)
		{
			dRand = rand_u();
			for(nj=0; nj<nTypeNum; nj++)
			{
				dRand -= DMGETAT(pPWM, ni, nj);
				if(dRand <= 0.0)
					break;
			}
			if(nj == nTypeNum)
				nj--;

			pBase[nk] = (unsigned char)nj;
			dBF += DMGETAT(pLogBGF, nWordId, nj);
			dMF += DMGETAT(pLogPWM, ni, nj);
			dMR += DMGETAT(pLogPWM, (nMotifLen-1-ni), (3-nj));

			nWordId -= pBase[nk-nBGOrder]*nScale;
			nWordId = nWordId*nTypeNum+nj;
			nk++;
		}

		for(; ni<nMotifLen; ni++)
		{
			dRand = rand_u();
			for(nj=0; nj<nTypeNum; nj++)
			{
				dRand -= DMGETAT(pPWM, ni, nj);
				if(dRand <= 0.0)
					break;
			}
			if(nj == nTypeNum)
				nj--;

			pBase[nk] = (unsigned char)nj;
			dBF += DMGETAT(pLogBGF, nWordId, nj);
			dMF += DMGETAT(pLogPWM, ni, nj);
			dMR += DMGETAT(pLogPWM, (nMotifLen-1-ni), (3-nj));

			nWordId -= pBase[nk-nBGOrder]*nScale;
			nWordId = nWordId*nTypeNum+nj;

			dBR += DMGETAT(pLogBGR, nWordId, pBase[nk-nBGOrder]);

			nk++;
		}

		for(ni=0; ni<nBGOrder; ni++)
		{
			dRand = rand_u();
			for(nj=0; nj<nTypeNum; nj++)
			{
				dRand -= DMGETAT(pBGF, nWordId, nj);
				if(dRand <= 0.0)
					break;
			}
			if(nj == nTypeNum)
				nj--;

			pBase[nk] = (unsigned char)nj;
			nWordId -= pBase[nk-nBGOrder]*nScale;
			nWordId = nWordId*nTypeNum+nj;
			dBR += DMGETAT(pLogBGR, nWordId, pBase[nk-nBGOrder]);

			nk++;
		}

		if(dMF > dMR)
			dML = dMF;
		else
			dML = dMR;

		if(dBF > dBR)
			dBL = dBF;
		else
			dBL = dBR;

		dScore = exp(dML-dBL);

	}
	else
	{
		printf("Error: MotifMap_SimuScore, nBGOrder<0!\n");
	}

	/* return */
	return dScore;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScoreRefGene: Motif mapping based on already mapped sites     */
/*  This function will score refgene based on known genomic motif mapping. */ 
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScoreRefGene(char strSpecies[], char strVersion[], 
		char strGenomePath[], char strCScorePath[], 
		char strTargetMotifList[], char strMotifMapPath[], double dCSCutoff, int nRepeatMask,
		char strRefGeneFile[], int nTSSUP, int nTSSDOWN, int nTESUP, int nTESDOWN, 
		int nIncludeIntron, int nIncludeExon,
		char strOutPath[])
{
	/* define */

	/* for indexing */
	char strChrLenFile[LINE_LENGTH];
	struct INTMATRIX *pChrLen;
	struct INTMATRIX *pChrContigNum;
	int nContigNum;
	FILE *fpMotifList;
	int nMotifNum;
	struct MOTIFSITEGROUP **vMotifContig;
	char strKnownMapFile[LINE_LENGTH];
	FILE *fpKnownMap;
	struct MOTIFSITE *pNewSite;
	int nLastContigId;
	struct MOTIFSITE *pLastSite;
	
	/* for mapping */
	FILE *fpRefGene;
	FILE *fpAnnotOut;
	FILE *fpScoreOut;
	FILE *fpCS;
	int numread;
	struct BYTEMATRIX *pCS;
	struct BYTEMATRIX *pCodeVec;
	int nCodeStart,nCodeEnd;
	char strConsFile[LINE_LENGTH];
	char strSeqCodeFile[LINE_LENGTH];
	struct tagGenomicRegion *pMapRegion;
	struct tagGenomicRegion *pRegionEle;
	struct DOUBLEMATRIX *pMapScore;
	char strRefLine[LONG_LINE_LENGTH];
	struct tagRefGene *pRefGene;
	struct MOTIFSITE *pMtfSite;
	int nFinished;
	int nEffectiveLen;
	int nRemainLen;

	/* general */
	char strChrName[LINE_LENGTH];
	int nChrId;
	int nContigId;
	int nStart,nEnd;
	char strAlias[LINE_LENGTH];
	double dScore;
	char chStrand;
	char strLine[LINE_LENGTH];
	int ni,nj,nk;
	int nMotifId;
	int nTemp;
	char *chp1,*chp2;

	/* for debug */
	FILE *fpDebug;
	fpDebug = NULL;
	sprintf(strLine, "%s_debug.txt", strOutPath);
	fpDebug = fopen(strLine, "wt");


	/* ---------------------------------------------- */
	/* Indexing the known motifs                      */
	/* ---------------------------------------------- */
	
	/* know genome size */
	sprintf(strChrLenFile, "%schrlen.txt", strGenomePath);
	pChrLen = NULL;
	pChrLen = IMLOAD(strChrLenFile);
	if(pChrLen == NULL)
	{
		printf("Error: MotifMap_ScoreRefGene, cannot load genomesize!\n");
		exit(EXIT_FAILURE);
	}

	pChrContigNum = NULL;
	pChrContigNum = CreateIntMatrix(pChrLen->nHeight, pChrLen->nWidth);
	if(pChrContigNum == NULL)
	{
		printf("Error: MotifMap_ScoreRefGene, cannot index the genome!\n");
		exit(EXIT_FAILURE);
	}

	nContigNum = 0;
	for(ni=0; ni<pChrLen->nHeight; ni++)
	{
		nTemp = pChrLen->pMatElement[ni]/MOTIFMAP_GENOMECONTIG_SIZE;
		if(pChrLen->pMatElement[ni]%MOTIFMAP_GENOMECONTIG_SIZE != 0)
		{
			nTemp += 1;
		}

		nContigNum += nTemp;
		IMSETAT(pChrContigNum, ni, 0, nContigNum);
	}

	/* know motif number */
	nMotifNum = 0;
	fpMotifList = NULL;
	fpMotifList = fopen(strTargetMotifList, "rt");
	if(fpMotifList == NULL)
	{
		printf("Error: MotifMap_ScoreRefGene, cannot open motif list!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LINE_LENGTH, fpMotifList) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		nMotifNum++;
	}

	fclose(fpMotifList);

	if(nMotifNum <= 0)
	{
		printf("Warning: MotifMap_ScoreRefGene, no motifs!\n");
		return PROC_FAILURE;
	}

	/* create index */
	vMotifContig = NULL;
	vMotifContig = (struct MOTIFSITEGROUP **)calloc(nContigNum, sizeof(struct MOTIFSITEGROUP *));
	if(vMotifContig == NULL)
	{
		printf("Error: MotifMap_ScoreRefGene, cannot create indexing system!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nContigNum; ni++)
	{
		vMotifContig[ni] = MTFSGRPCREATE(nMotifNum);
		if(vMotifContig[ni] == NULL)
		{
			printf("Error: MotifMap_ScoreRefGene, cannot create indexing system!\n");
			exit(EXIT_FAILURE);
		}
		vMotifContig[ni]->nGroupId = ni;
	}

	/* initialize the indexing system */
	nMotifId = 0;
	fpMotifList = NULL;
	fpMotifList = fopen(strTargetMotifList, "rt");
	if(fpMotifList == NULL)
	{
		printf("Error: MotifMap_ScoreRefGene, cannot open motif list!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LINE_LENGTH, fpMotifList) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		
		/* open motif mapping file */
		sprintf(strKnownMapFile, "%s%s_%s_%s_map.bed", strMotifMapPath, strLine, strSpecies, strVersion);
		fpKnownMap = NULL;
		fpKnownMap = fopen(strKnownMapFile, "rt");
		if(fpKnownMap == NULL)
		{
			printf("Error: MotifMap_ScoreRefGene, cannot open known motif map!\n");
			exit(EXIT_FAILURE);
		}

		/* load head info */
		if(fgets(strLine, LINE_LENGTH, fpKnownMap) == NULL)
		{
			fclose(fpKnownMap);
			continue;
		}
		if(fgets(strLine, LINE_LENGTH, fpKnownMap) == NULL)
		{
			fclose(fpKnownMap);
			continue;
		}

		/* load motif sites */
		nLastContigId = -1;
		pLastSite = NULL;
	
		while(fgets(strLine, LINE_LENGTH, fpKnownMap) != NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;

			sscanf(strLine, "%s %d %d %s %lf %c", 
				strChrName, &nStart, &nEnd, strAlias, &dScore, &chStrand);

			nChrId = Genome_ChromosomeName_To_Index(strChrName, strSpecies)-1;
			nContigId = MotifMap_GetContigIndex(nChrId, nStart, pChrContigNum, MOTIFMAP_GENOMECONTIG_SIZE);

			pNewSite = NULL;
			pNewSite = MTFSCREATE();
			if(pNewSite == NULL)
			{
				printf("Error: MotifMap_ScoreRefGene, cannot create new motif site!\n");
				exit(EXIT_FAILURE);
			}
			pNewSite->dScore = dScore;
			pNewSite->nMotifType = nMotifId;
			pNewSite->nStartPos = nStart;
			pNewSite->nMotifLen = nEnd-nStart+1;
			if(chStrand == '-')
				pNewSite->nStrand = 1;
			else
				pNewSite->nStrand = 0;

			/* add to indexing system */
			if(nContigId != nLastContigId)
			{
				if(vMotifContig[nContigId]->pSites[nMotifId] != NULL)
				{
					printf("Error: MotifMap_ScoreRefGene, indexing system not initialized!\n");
					exit(EXIT_FAILURE);
				}
				vMotifContig[nContigId]->pSites[nMotifId] = pNewSite;
				nLastContigId = nContigId;
				pLastSite = pNewSite;
			}
			else
			{
				pLastSite->pNext = pNewSite;
				pLastSite = pNewSite;
			}
		}

		/* close file */
		fclose(fpKnownMap);

		nMotifId++;
	}

	fclose(fpMotifList);

	if(nMotifId != nMotifNum)
	{
		printf("Error: MotifMap_ScoreRefGene, motif number not match!\n");
		exit(EXIT_FAILURE);
	}


	/* ---------------------------------------------- */
	/* Load refgenes and score them                   */
	/* ---------------------------------------------- */
	
	pCS = NULL;
	pCS = CreateByteMatrix(1, MOTIFMAP_GENOMECONTIG_SIZE);
	if(pCS == NULL)
	{
		printf("Error: MotifMap_ScoreRefGene, cannot create conservation score space!\n");
		exit(EXIT_FAILURE);
	}

	pCodeVec = NULL;
	if(nRepeatMask == 1)
	{
		pCodeVec = CreateByteMatrix(1, MOTIFMAP_GENOMECONTIG_SIZE);
		if(pCodeVec == NULL)
		{
			printf("Error: MotifMap_ScoreRefGene, cannot create sequence space!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* get map region */
	fpRefGene = NULL;
	fpRefGene = fopen(strRefGeneFile, "rt");
	if(fpRefGene == NULL)
	{
		printf("Error: MotifMap_ScoreRefGene, cannot open refgene list!\n");
		exit(EXIT_FAILURE);
	}

	sprintf(strLine, "%s.txt", strOutPath);
	fpScoreOut = NULL;
	fpScoreOut = fopen(strLine, "wt");
	if(fpScoreOut == NULL)
	{
		printf("Error: MotifMap_ScoreRefGene, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	
	sprintf(strLine, "%s.info", strOutPath);
	fpAnnotOut = NULL;
	fpAnnotOut = fopen(strLine, "wt");
	if(fpAnnotOut == NULL)
	{
		printf("Error: MotifMap_ScoreRefGene, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}


	/* for every motif, map and add up scores */
	pMapScore = NULL;
	pMapScore = CreateDoubleMatrix(1, nMotifNum);
	if(pMapScore == NULL)
	{
		printf("Error: MotifMap_ScoreRefGene, cannot create map score matrix!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpScoreOut, "RefGene\tRefLine\tProbeset\tProbeLine\tClass\tScore\tEffecLen");
	fpMotifList = NULL;
	fpMotifList = fopen(strTargetMotifList, "rt");
	if(fpMotifList == NULL)
	{
		printf("Error: MotifMap_ScoreRefGene, cannot open motif list!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LINE_LENGTH, fpMotifList) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		fprintf(fpScoreOut, "\t%s", strLine);
	}
	fclose(fpMotifList);
	fprintf(fpScoreOut, "\n");

	while(fgets(strRefLine, LONG_LINE_LENGTH, fpRefGene) != NULL)
	{
		StrTrimLeft(strRefLine);
		StrTrimRight(strRefLine);
		if(strRefLine[0] == '\0')
			continue;

		/* load head info */
		chp1 = strRefLine;
		for(nk=0; nk<6; nk++)
		{
			chp2 = strchr(chp1, '\t');
			chp1 = chp2+1;
		}

		*chp2 = '\0';
		fprintf(fpScoreOut, "%s\t", strRefLine);
		fprintf(fpAnnotOut, "%s\t", strRefLine);


		/* load refgene */
		pRefGene = NULL;
		pRefGene = RefGeneCreate();
		if(pRefGene == NULL)
		{
			printf("Error: MotifMap_ScoreRefGene, cannot create refgene!\n");
			exit(EXIT_FAILURE);
		}

		RefGeneInit_FromGenomeLabFormat(pRefGene, chp1, strSpecies);
		
		/* get map regions */
		pMapRegion = NULL;
		pMapRegion = RefGeneSelectRegion(pRefGene, nTSSUP, nTSSDOWN, nTESUP, nTESDOWN, 
			nIncludeIntron, nIncludeExon, pChrLen);

		if(pMapRegion == NULL)
		{
			RefGeneDestroy(pRefGene);
			continue;
		}

		/* get contig index */
		nContigId = MotifMap_GetContigIndex((pMapRegion->nChrom-1), pMapRegion->nStart, 
			pChrContigNum, MOTIFMAP_GENOMECONTIG_SIZE);
		
		/* TODO: get effective length */
		nEffectiveLen = 0;
		pRegionEle = pMapRegion;
		while(pRegionEle != NULL)
		{
			/* write debug file */
			fprintf(fpDebug, "%s\t%s\t%d\t%d\t%c\n", 
				pRefGene->strName, pRegionEle->strChrom, pRegionEle->nStart, pRegionEle->nEnd, pRegionEle->chStrand);
			
			/* get effective len */
			sprintf(strConsFile, "%s%s.cs", strCScorePath, pRegionEle->strChrom);
			sprintf(strSeqCodeFile, "%s%s.sq", strGenomePath, pRegionEle->strChrom);
			fpCS = NULL;
			fpCS = fopen(strConsFile, "rb");
			if(fpCS == NULL)
			{
				printf("Error: MotifMap_ScoreRefGene, cannot open conservation score file!\n");
				exit(EXIT_FAILURE);
			}

			/* load score */
			if( fseek( fpCS, pRegionEle->nStart, SEEK_SET ) != 0 )
			{
				printf("Error: MotifMap_ScoreRefGene, cannot locate the required sequence!\n");
				exit(EXIT_FAILURE);
			}

			/* read */
			nCodeStart = pRegionEle->nStart;
			nRemainLen = pRegionEle->nEnd-pRegionEle->nStart+1;
			while(nRemainLen > 0)
			{
				if(nRemainLen >= MOTIFMAP_GENOMECONTIG_SIZE)
				{
					numread = fread(pCS->pMatElement, sizeof(unsigned char), MOTIFMAP_GENOMECONTIG_SIZE, fpCS);
					if(numread != MOTIFMAP_GENOMECONTIG_SIZE)
					{
						printf("Error: MotifMap_ScoreRefGene, loading error!\n");
						exit(EXIT_FAILURE);
					}

					nCodeEnd = nCodeStart+MOTIFMAP_GENOMECONTIG_SIZE-1;
					if(nRepeatMask == 1)
					{
						numread = Genome_Code_4bit_GetUncompressed(pCodeVec, strSeqCodeFile, nCodeStart, nCodeEnd);
						if(numread != MOTIFMAP_GENOMECONTIG_SIZE)
						{
							printf("Error: MotifMap_ScoreRefGene, loading error!\n");
							exit(EXIT_FAILURE);
						}
					}

					nRemainLen -= MOTIFMAP_GENOMECONTIG_SIZE;
				}
				else
				{
					numread = fread(pCS->pMatElement, sizeof(unsigned char), nRemainLen, fpCS);
					if(numread != nRemainLen)
					{
						printf("Error: MotifMap_ScoreRefGene, loading error!\n");
						exit(EXIT_FAILURE);
					}

					nCodeEnd = nCodeStart+nRemainLen-1;
					if(nRepeatMask == 1)
					{
						numread = Genome_Code_4bit_GetUncompressed(pCodeVec, strSeqCodeFile, nCodeStart, nCodeEnd);
						if(numread != nRemainLen)
						{
							printf("Error: MotifMap_ScoreRefGene, loading error!\n");
							exit(EXIT_FAILURE);
						}
					}

					nRemainLen -= nRemainLen;
				}

				nCodeStart = nCodeEnd+1;

				if(nRepeatMask == 1)
				{
					for(nk=0; nk<numread; nk++)
					{
						if((pCS->pMatElement[nk] >= dCSCutoff) && (pCodeVec->pMatElement[nk] < 4))
							nEffectiveLen++;
					}
				}
				else
				{
					for(nk=0; nk<numread; nk++)
					{
						if(pCS->pMatElement[nk] >= dCSCutoff)
							nEffectiveLen++;
					}
				}
			}

			fclose(fpCS);

			/* get next */
			pRegionEle = pRegionEle->pNext;
		}

		fprintf(fpScoreOut, "%d\t", nEffectiveLen);
		fprintf(fpAnnotOut, "%d\n", nEffectiveLen);



		/* scoring motifs */
		for(ni=0; ni<nMotifNum; ni++)
		{
			dScore = 0.0;
			nFinished = 0;
			pRegionEle = pMapRegion;
			nj = nContigId;
			
			while(nj < pChrContigNum->pMatElement[pRefGene->nChrom-1])
			{
				pMtfSite = vMotifContig[nj]->pSites[ni];
				while(pMtfSite != NULL)
				{
					if(pRegionEle == NULL)
					{
						nFinished = 1;
						break;
					}

					if(pMtfSite->nStartPos > pRegionEle->nEnd)
					{
						pRegionEle = pRegionEle->pNext;
						continue;
					}

					if(pMtfSite->nStartPos >= pRegionEle->nStart)
					{
						/* dScore += pMtfSite->dScore; */
						dScore += 1.0;
					}

					/* get next motif site */
					pMtfSite = pMtfSite->pNext;
				}

				if(nFinished == 1)
					break;

				/* get next contig */
				nj++;
			}

			/* dScore = log10(dScore+1.0); */
			fprintf(fpScoreOut, "% 9.7e ", dScore);
		}

		fprintf(fpScoreOut, "\n");

		/* clear */
		RefGeneDestroy(pRefGene);
		while(pMapRegion != NULL)
		{
			pRegionEle = pMapRegion;
			pMapRegion = pMapRegion->pNext;
			GenomicRegionDestroy(pRegionEle);
		}
	}

	/* close file */
	DestroyByteMatrix(pCS);
	if(nRepeatMask == 1)
	{
		DestroyByteMatrix(pCodeVec);
	}
	DestroyDoubleMatrix(pMapScore);
	fclose(fpRefGene);
	fclose(fpScoreOut);
	fclose(fpAnnotOut);
	fclose(fpDebug);

	/* ---------------------------------------------- */
	/* Release memories                               */
	/* ---------------------------------------------- */

	/* destroy indexing system */
	DestroyIntMatrix(pChrLen);
	DestroyIntMatrix(pChrContigNum);
	for(ni=0; ni<nContigNum; ni++)
	{
		MTFSGRPDESTROY(vMotifContig[ni]);
		vMotifContig[ni] = NULL;
	}
	free(vMotifContig);


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_GetContigIndex: get contig index for a position.              */ 
/* ----------------------------------------------------------------------- */ 
int MotifMap_GetContigIndex(int nChrId, int nStart, 
							struct INTMATRIX *pChrContigNum, int nContigSize)
{
	/* define */
	int nIndex;
	int nTemp;

	/* process */
	if(nChrId >= pChrContigNum->nHeight)
	{
		printf("Error: MotifMap_GetContigIndex, chromosome index out of range!\n");
		exit(EXIT_FAILURE);
	}

	nTemp = nStart/nContigSize;
	if(nChrId == 0)
	{
		if(nTemp >= pChrContigNum->pMatElement[0])
		{
			printf("Warning: MotifMap_GetContigIndex, indexing out of range!\n");
			nIndex = pChrContigNum->pMatElement[0]-1;
		}
		else
		{
			nIndex = nTemp;
		}
	}
	else
	{
		if(nTemp >= (pChrContigNum->pMatElement[nChrId]-pChrContigNum->pMatElement[nChrId-1]))
		{
			printf("Warning: MotifMap_GetContigIndex, indexing out of range!\n");
			nIndex = pChrContigNum->pMatElement[nChrId]-1;
		}
		else
		{
			nIndex = pChrContigNum->pMatElement[nChrId-1]+nTemp;
		}
	}

	/* return */
	return nIndex;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_GetStatistics: read summary statistics from motifmap.         */ 
/* ----------------------------------------------------------------------- */ 
int MotifMap_GetStatistics(char strFileName[], double *pEffecLen, double *pSiteNum)
{
	/* define */
	FILE *fpIn;
	char strLine[MED_LINE_LENGTH];
	char *chSep;

	/* open file */
	fpIn = NULL;
	fpIn = fopen(strFileName, "r");
	if(fpIn == NULL)
	{
		printf("Error: MotifMap_GetStatistics, cannot open summary file!\n");
		exit(EXIT_FAILURE);
	}

	/* read file */
	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strstr(strLine, "EffecLen") == strLine)
		{
			chSep = strchr(strLine, '=');
			chSep++;
			StrTrimLeft(chSep);
			*pEffecLen = atof(chSep);
		}
		else if(strstr(strLine, "TotalSite") == strLine)
		{
			chSep = strchr(strLine, '=');
			chSep++;
			StrTrimLeft(chSep);
			*pSiteNum = atof(chSep);
		}
		else if(strstr(strLine, "ConsLen") == strLine)
		{
		}
		else
		{
			printf("Error: MotifMap_GetStatistics, unknown summary statistics!\n");
			exit(EXIT_FAILURE);
		}
	}


	/* close file */
	fclose(fpIn);

	/* return */
	return PROC_SUCCESS;
}

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
				int nUseCS, char strCSPath[], int nIncludeRepeat, int nW)
{
	/* define */
	/* background type, 0:region, 1:genome */
	int nBGType = 0;

	char strPosFile[LINE_LENGTH];
	char strNegFile[LINE_LENGTH];
	int nPosOnly = 0;
	
	int nMotifNum = 0;
	int nTotSiteNum;
	struct tagString **vMotifName = NULL;
	struct tagString **vMotifFile = NULL;
	struct DOUBLEMATRIX *pMotifR = NULL;
	struct DOUBLEMATRIX *pMotifC = NULL;
	char strMotifName[LINE_LENGTH];
	char strMotifFile[LINE_LENGTH];
	double dR;
	double dC;
	double dELen;
	double dSNum;

	/* summary statistics */
	struct DOUBLEMATRIX *pEffecLenPos;
	struct DOUBLEMATRIX *pSiteNumPos;
	struct DOUBLEMATRIX *pNeighborNumPos;
	struct DOUBLEMATRIX *pNeighborLenPos;
	struct DOUBLEMATRIX *pNeighborSiteNumPos;

	struct DOUBLEMATRIX *pEffecLenNeg;
	struct DOUBLEMATRIX *pSiteNumNeg;
	struct DOUBLEMATRIX *pNeighborNumNeg;

	/* Merged file */
	FILE *fpOut;

	/* others */
	char strFileName[MED_LINE_LENGTH];
	char strMotifPath[MED_LINE_LENGTH];
	char strCodPath[MED_LINE_LENGTH];
	char strMotifMapPath[MED_LINE_LENGTH];
	char strMotifMapNPath[MED_LINE_LENGTH];
	char strCombinedPath[MED_LINE_LENGTH];
	char strSortedPath[MED_LINE_LENGTH];
	char strCommand[MED_LINE_LENGTH];
	

	int nSystemCode;
	int nError = 0;
	FILE *fpIn;
	char strLine[MED_LINE_LENGTH];
	int ni,nj;
	char *chSep;
	
	/* init */
	if(strcmp(strBGType, "GENOME") == 0)
	{
		nBGType = 1;
	}
	else
	{
		nBGType = 0;
	}

	AdjustDirectoryPath(strGenomePath);
	AdjustDirectoryPath(strWorkPath);
	if(nUseCS == 1)
	{
		AdjustDirectoryPath(strCSPath);
	}
	if(nBGType == 1)
	{
		AdjustDirectoryPath(strBGPath);
	}

	/* load parameters */
	sprintf(strFileName, "%s%s", strWorkPath, strInputPath);
	fpIn = NULL;
	fpIn = fopen(strFileName, "r");
	if(fpIn == NULL)
	{
		printf("Error: MotifMap_ScanMultipleMatrix_Genome_Main, cannot open data file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, MED_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] =='\0')
			continue;

		if(strstr(strLine, "[Positive Region]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			strcpy(strPosFile, chSep);
		}
		else if(strstr(strLine, "[Negative Region]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			strcpy(strNegFile, chSep);
		}
		else if(strstr(strLine, "[Motif Number]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			nMotifNum = atoi(chSep);
			if(nMotifNum <=0)
			{
				printf("Error: motif number <= 0\n");
				nError = 1;
				break;
			}

			/* prepare space */
			vMotifName = NULL;
			vMotifName = (struct tagString **)calloc(nMotifNum, sizeof(struct tagString *));
			if(vMotifName == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}

			vMotifFile = NULL;
			vMotifFile = (struct tagString **)calloc(nMotifNum, sizeof(struct tagString *));
			if(vMotifFile == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}

			pMotifR = NULL;
			pMotifR = CreateDoubleMatrix(1, nMotifNum);
			if(pMotifR == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}

			pMotifC = NULL;
			pMotifC = CreateDoubleMatrix(1, nMotifNum);
			if(pMotifC == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
		}

		else if(strstr(strLine, "[Motif]") == strLine)
		{
			ni = 0;
			while(fgets(strLine, MED_LINE_LENGTH, fpIn)!=NULL)
			{
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] =='\0')
					continue;

				if(ni >= nMotifNum)
				{
					printf("Error: MotifMap_ScanMultipleMatrix_Genome_Main, motif number not match!\n");
					exit(EXIT_FAILURE);
				}

				sscanf(strLine, "%s %s %lf %lf", strMotifName, strMotifFile, &dR, &dC);
				StringAddTail(vMotifName+ni, strMotifName);
				StringAddTail(vMotifFile+ni, strMotifFile);
				pMotifR->pMatElement[ni] = dR;
				pMotifC->pMatElement[ni] = dC;

				ni++;
			}

			if(ni != nMotifNum)
			{
				printf("Error: MotifMap_ScanMultipleMatrix_Genome_Main, motif number not match!\n");
				exit(EXIT_FAILURE);
			}
		}

		else
		{
			printf("%s\n", strLine);
			printf("Error: MotifMap_ScanMultipleMatrix_Genome_Main, cannot load parameter file correctly!\n");
			exit(EXIT_FAILURE);
		}

	}

	fclose(fpIn);

	if(nError == 1)
	{
		printf("Error: MotifMap_ScanMultipleMatrix_Genome_Main, cannot load parameter file correctly!\n");
		exit(EXIT_FAILURE);
	}

	if(strcmp(strNegFile, "NULL") == 0)
	{
		nPosOnly = 1;
	}
	else
	{
		nPosOnly = 0;
	}

	/* initialize summary statistics matrix */
	pEffecLenPos = NULL;
	pEffecLenPos = CreateDoubleMatrix(1, nMotifNum);
	if(pEffecLenPos == NULL)
	{
		printf("Error: MotifMap_ScanMultipleMatrix_Genome_Main, cannot create matrix for summary statistics!\n");
		exit(EXIT_FAILURE);
	}

	pSiteNumPos = NULL;
	pSiteNumPos = CreateDoubleMatrix(1, nMotifNum);
	if(pSiteNumPos == NULL)
	{
		printf("Error: MotifMap_ScanMultipleMatrix_Genome_Main, cannot create matrix for summary statistics!\n");
		exit(EXIT_FAILURE);
	}

	pNeighborNumPos = NULL;
	pNeighborNumPos = CreateDoubleMatrix(nMotifNum, nMotifNum);
	if(pNeighborNumPos == NULL)
	{
		printf("Error: MotifMap_ScanMultipleMatrix_Genome_Main, cannot create matrix for summary statistics!\n");
		exit(EXIT_FAILURE);
	}

	pNeighborLenPos = NULL;
	pNeighborLenPos = CreateDoubleMatrix(nMotifNum, nMotifNum);
	if(pNeighborLenPos == NULL)
	{
		printf("Error: MotifMap_ScanMultipleMatrix_Genome_Main, cannot create matrix for summary statistics!\n");
		exit(EXIT_FAILURE);
	}

	pNeighborSiteNumPos = NULL;
	pNeighborSiteNumPos = CreateDoubleMatrix(nMotifNum, nMotifNum);
	if(pNeighborSiteNumPos == NULL)
	{
		printf("Error: MotifMap_ScanMultipleMatrix_Genome_Main, cannot create matrix for summary statistics!\n");
		exit(EXIT_FAILURE);
	}

	pEffecLenNeg = NULL;
	pEffecLenNeg = CreateDoubleMatrix(1, nMotifNum);
	if(pEffecLenNeg == NULL)
	{
		printf("Error: MotifMap_ScanMultipleMatrix_Genome_Main, cannot create matrix for summary statistics!\n");
		exit(EXIT_FAILURE);
	}

	pSiteNumNeg = NULL;
	pSiteNumNeg = CreateDoubleMatrix(1, nMotifNum);
	if(pSiteNumNeg == NULL)
	{
		printf("Error: MotifMap_ScanMultipleMatrix_Genome_Main, cannot create matrix for summary statistics!\n");
		exit(EXIT_FAILURE);
	}

	pNeighborNumNeg = NULL;
	pNeighborNumNeg = CreateDoubleMatrix(nMotifNum, nMotifNum);
	if(pNeighborNumNeg == NULL)
	{
		printf("Error: MotifMap_ScanMultipleMatrix_Genome_Main, cannot create matrix for summary statistics!\n");
		exit(EXIT_FAILURE);
	}

	/* map individual motifs to positive sequences */
	sprintf(strCombinedPath, "%s%s_pos_n.map", strWorkPath, strOutputPath);
	sprintf(strSortedPath, "%s%s_pos_n_sorted.map", strWorkPath, strOutputPath);
	fpOut = NULL;
	fpOut = fopen(strCombinedPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: MotifMap_ScanMultipleMatrix_Genome_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	printf("Positive Mapping:\n");
	nTotSiteNum = 0;
	for(ni=0; ni<nMotifNum; ni++)
	{
		printf("mapping %s...\n", vMotifName[ni]->m_pString);
		sprintf(strMotifPath, "%s%s", strWorkPath, vMotifFile[ni]->m_pString);
		sprintf(strCodPath, "%s%s", strWorkPath, strPosFile);
		sprintf(strMotifMapPath, "%s%s_pos.map", strWorkPath, vMotifName[ni]->m_pString);
		sprintf(strMotifMapNPath, "%s%s_pos_n.map", strWorkPath, vMotifName[ni]->m_pString);
		dR = pMotifR->pMatElement[ni];
		dC = pMotifC->pMatElement[ni];

		/* map motif */
		MotifMap_ScanMatrix_Genome_Main(strMotifPath, strGenomePath, 
					strCodPath, strMotifMapPath, dR, 
					nBGOrder, strBGType, 
					strBGPath, nBGStepSize,
					nUseCS, dC, strCSPath, nIncludeRepeat);	
		
		MotifMap_GetSiteAround_Main(strMotifMapPath, strGenomePath, 
			strMotifMapNPath, strSpecies, 0, 1, 0);
		
		/* get stat */
		sprintf(strFileName, "%s.stat", strMotifMapPath);
		MotifMap_GetStatistics(strFileName, pEffecLenPos->pMatElement+ni, pSiteNumPos->pMatElement+ni);
		nTotSiteNum += (int)(pSiteNumPos->pMatElement[ni]);
		
		/* merge files */
		fpIn = NULL;
		fpIn = fopen(strMotifMapNPath, "r");
		if(fpIn == NULL)
		{
			printf("Error: MotifMap_ScanMultipleMatrix_Genome_Main, cannot open motif mapping file!\n");
			exit(EXIT_FAILURE);
		}
		while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;
			if(strLine[0] == '#')
				continue;

			fprintf(fpOut, "%s\t%d\n", strLine, ni);
		}

		fclose(fpIn);
		
		/* get neighboring window */
		sprintf(strCodPath, "%s%s_pos_%d.cod", strWorkPath, vMotifName[ni]->m_pString, nW);
		MotifMap_GetSiteAround_Main(strMotifMapPath, strGenomePath, 
			strCodPath, strSpecies, nW, 0, 0);

		/* map motif to neighboring window */
		for(nj=0; nj<nMotifNum; nj++)
		{
			sprintf(strMotifPath, "%s%s", strWorkPath, vMotifFile[nj]->m_pString);
			sprintf(strMotifMapPath, "%s%s-in-%s_pos_%d.map", strWorkPath, 
				vMotifName[nj]->m_pString, vMotifName[ni]->m_pString, nW);
			dR = pMotifR->pMatElement[nj];
			dC = pMotifC->pMatElement[nj];

			MotifMap_ScanMatrix_Genome_Main(strMotifPath, strGenomePath, 
						strCodPath, strMotifMapPath, dR, 
						nBGOrder, strBGType, 
						strBGPath, nBGStepSize,
						nUseCS, dC, strCSPath, nIncludeRepeat);
			
			sprintf(strFileName, "%s.stat", strMotifMapPath);
			MotifMap_GetStatistics(strFileName, &dELen, &dSNum);
			DMSETAT(pNeighborLenPos, ni, nj, dELen);
			DMSETAT(pNeighborSiteNumPos, ni, nj, dSNum);
			
			if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
			{
				sprintf(strCommand, "del %s*", strMotifMapPath);
				nSystemCode = system(strCommand);
			}
			else
			{
				sprintf(strCommand, "rm %s*", strMotifMapPath);
				nSystemCode = system(strCommand);
			}
		}

		if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
		{
			sprintf(strCommand, "del %s", strCodPath);
			nSystemCode = system(strCommand);
		}
		else
		{
			sprintf(strCommand, "rm %s", strCodPath);
			nSystemCode = system(strCommand);
		}
	}

	fclose(fpOut);

	/* get pair wise statistics for positive mapping */
	sprintf(strCommand, "sort -o %s -n +1 +2 +3 %s", strSortedPath, strCombinedPath);
	nSystemCode = system(strCommand);
	MotifMap_ScanMultipleMatrix_ClusterTest(nMotifNum, vMotifName, 
		strSortedPath, nTotSiteNum, nW, pNeighborNumPos);
	

	/* map individual motifs to negative sequences */
	if(nPosOnly == 0)
	{
		sprintf(strCombinedPath, "%s%s_neg_n.map", strWorkPath, strOutputPath);
		sprintf(strSortedPath, "%s%s_neg_n_sorted.map", strWorkPath, strOutputPath);
		fpOut = NULL;
		fpOut = fopen(strCombinedPath, "w");
		if(fpOut == NULL)
		{
			printf("Error: MotifMap_ScanMultipleMatrix_Genome_Main, cannot open output file!\n");
			exit(EXIT_FAILURE);
		}

		printf("Negative Mapping:\n");
		nTotSiteNum = 0;
		for(ni=0; ni<nMotifNum; ni++)
		{
			printf("mapping %s...\n", vMotifName[ni]->m_pString);
			sprintf(strMotifPath, "%s%s", strWorkPath, vMotifFile[ni]->m_pString);
			sprintf(strCodPath, "%s%s", strWorkPath, strNegFile);
			sprintf(strMotifMapPath, "%s%s_neg.map", strWorkPath, vMotifName[ni]->m_pString);
			sprintf(strMotifMapNPath, "%s%s_neg_n.map", strWorkPath, vMotifName[ni]->m_pString);
			dR = pMotifR->pMatElement[ni];
			dC = pMotifC->pMatElement[ni];

			/* map motif */
			MotifMap_ScanMatrix_Genome_Main(strMotifPath, strGenomePath, 
						strCodPath, strMotifMapPath, dR, 
						nBGOrder, strBGType, 
						strBGPath, nBGStepSize,
						nUseCS, dC, strCSPath, nIncludeRepeat);
			
			MotifMap_GetSiteAround_Main(strMotifMapPath, strGenomePath, 
				strMotifMapNPath, strSpecies, 0, 1, 0);
			
			/* get stat */
			sprintf(strFileName, "%s.stat", strMotifMapPath);
			MotifMap_GetStatistics(strFileName, pEffecLenNeg->pMatElement+ni, pSiteNumNeg->pMatElement+ni);
			nTotSiteNum += (int)(pSiteNumNeg->pMatElement[ni]);
			
			/* merge files */
			fpIn = NULL;
			fpIn = fopen(strMotifMapNPath, "r");
			if(fpIn == NULL)
			{
				printf("Error: MotifMap_ScanMultipleMatrix_Genome_Main, cannot open motif mapping file!\n");
				exit(EXIT_FAILURE);
			}
			while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
			{
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] == '\0')
					continue;
				if(strLine[0] == '#')
					continue;

				fprintf(fpOut, "%s\t%d\n", strLine, ni);
			}

			fclose(fpIn);
		}

		fclose(fpOut);

		/* get pair wise statistics for negative mapping */
		sprintf(strCommand, "sort -o %s -n +1 +2 +3 %s", strSortedPath, strCombinedPath);
		nSystemCode = system(strCommand);
		MotifMap_ScanMultipleMatrix_ClusterTest(nMotifNum, vMotifName, 
			strSortedPath, nTotSiteNum, nW, pNeighborNumNeg);
	}


	/* export results */
	sprintf(strFileName, "%s%s.stat", strWorkPath, strOutputPath);
	MotifMap_ScanMultipleMatrix_Export(strFileName, nMotifNum, 
		vMotifName, nPosOnly,
		pEffecLenPos, pSiteNumPos, pNeighborNumPos,
		pNeighborLenPos, pNeighborSiteNumPos,
		pEffecLenNeg, pSiteNumNeg, pNeighborNumNeg);
	

	/* release memory */
	for(ni=0; ni<nMotifNum; ni++)
	{
		DeleteString(vMotifName[ni]);
		DeleteString(vMotifFile[ni]);
	}
	free(vMotifName);
	free(vMotifFile);
	DestroyDoubleMatrix(pMotifR);
	DestroyDoubleMatrix(pMotifC);

	DestroyDoubleMatrix(pEffecLenPos);
	DestroyDoubleMatrix(pSiteNumPos);
	DestroyDoubleMatrix(pNeighborNumPos);
	DestroyDoubleMatrix(pNeighborLenPos);
	DestroyDoubleMatrix(pNeighborSiteNumPos);

	DestroyDoubleMatrix(pEffecLenNeg);
	DestroyDoubleMatrix(pSiteNumNeg);
	DestroyDoubleMatrix(pNeighborNumNeg);


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMultipleMatrix_ClusterTest: test clustering of motif site */ 
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanMultipleMatrix_ClusterTest(int nMotifNum, struct tagString **vMotifName, 
			char strSortedSitePath[], int nTotSiteNum, int nW, struct DOUBLEMATRIX *pNeighborNum)
{
	/* define */
	FILE *fpIn;
	struct FLEXMOTIFSITE **vSites;
	int ni,nj,nk;
	char strLine[MED_LINE_LENGTH];
	char strSeqId[LINE_LENGTH];
	int nChr,nStart,nEnd;
	char chStrand;
	double dScore;
	char strSiteSeq[LINE_LENGTH];
	int nSiteType;
	struct INTMATRIX *pNeighborCount;
	double dTemp;

	/* init site database */
	if(nTotSiteNum <= 0)
	{
		printf("Warning: MotifMap_ScanMultipleMatrix_ClusterTest, there are %d (<=0) motif sites!\n", nTotSiteNum);
		return PROC_SUCCESS;
	}
	vSites = NULL;
	vSites = (struct FLEXMOTIFSITE **)calloc(nTotSiteNum, sizeof(struct FLEXMOTIFSITE *));
	if(vSites == NULL)
	{
		printf("Error: MotifMap_ScanMultipleMatrix_ClusterTest, cannot allocate memory for loading site mapping!\n");
		exit(EXIT_FAILURE);
	}

	/* load file */
	fpIn = NULL;
	fpIn = fopen(strSortedSitePath, "r");
	if(fpIn == NULL)
	{
		printf("Error: MotifMap_ScanMultipleMatrix_ClusterTest, cannot open file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(ni >= nTotSiteNum)
		{
			printf("Error: MotifMap_ScanMultipleMatrix_ClusterTest, site number not match!\n");
			exit(EXIT_FAILURE);
		}

		sscanf(strLine, "%s %d %d %d %c %lf %s %d", strSeqId, 
			&nChr, &nStart, &nEnd, &chStrand, &dScore, strSiteSeq, &nSiteType);
		vSites[ni] = FLEXMTFSCREATE();
		vSites[ni]->dScore = dScore;
		vSites[ni]->nMotifType = nSiteType;
		vSites[ni]->nSeqId = nChr;
		vSites[ni]->nStartPos = (nStart+nEnd)/2;
		if(chStrand == '-')
			vSites[ni]->nStrand = 1;
		else
			vSites[ni]->nStrand = 0;

		ni++;
	}

	fclose(fpIn);

	if(ni != nTotSiteNum)
	{
		printf("Error: MotifMap_ScanMultipleMatrix_ClusterTest, site number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* count neighbors */
	for(ni=0; ni<nTotSiteNum; ni++)
	{
		pNeighborCount = NULL;
		pNeighborCount = CreateIntMatrix(1, nMotifNum);
		if(pNeighborCount == NULL)
		{
			printf("Error: MotifMap_ScanMultipleMatrix_ClusterTest, cannot create counting matrix!\n");
			exit(EXIT_FAILURE);
		}

		/* search forward */
		for(nj=ni-1; nj>=0; nj--)
		{
			if(vSites[nj]->nSeqId != vSites[ni]->nSeqId)
				break;
			if(fabs(vSites[nj]->nStartPos-vSites[ni]->nStartPos) > nW)
				break;
			nk = vSites[nj]->nMotifType;
			pNeighborCount->pMatElement[nk] += 1;
		}

		/* search backward */
		for(nj=ni+1; nj<nTotSiteNum; nj++)
		{
			if(vSites[nj]->nSeqId != vSites[ni]->nSeqId)
				break;
			if(fabs(vSites[nj]->nStartPos-vSites[ni]->nStartPos) > nW)
				break;
			nk = vSites[nj]->nMotifType;
			pNeighborCount->pMatElement[nk] += 1;
		}

		/* update summary stat */
		nk = vSites[ni]->nMotifType;
		for(nj=0; nj<nMotifNum; nj++)
		{
			if(pNeighborCount->pMatElement[nj] > 0)
			{
				dTemp = DMGETAT(pNeighborNum, nk, nj)+1.0;
				DMSETAT(pNeighborNum, nk, nj, dTemp);
			}
		}

		DestroyIntMatrix(pNeighborCount);
	}


	/* free memory */
	for(ni=0; ni<nTotSiteNum; ni++)
	{
		FLEXMTFSDESTROY(vSites[ni]);	
	}
	free(vSites);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_ScanMultipleMatrix_Export: save summary statistics.           */ 
/* ----------------------------------------------------------------------- */ 
int MotifMap_ScanMultipleMatrix_Export(char strFileName[], int nMotifNum, 
		struct tagString **vMotifName, int nPosOnly,
		struct DOUBLEMATRIX *pEffecLenPos, struct DOUBLEMATRIX *pSiteNumPos, 
		struct DOUBLEMATRIX *pNeighborNumPos, struct DOUBLEMATRIX *pNeighborLenPos, 
		struct DOUBLEMATRIX *pNeighborSiteNumPos, struct DOUBLEMATRIX *pEffecLenNeg, 
		struct DOUBLEMATRIX *pSiteNumNeg, struct DOUBLEMATRIX *pNeighborNumNeg)
{
	/* define */
	FILE *fpOut;
	int ni,nj;
	int nEffecLenPos,nSiteNumPos,nEffecLenNeg,nSiteNumNeg,nNeighborNumPos,nNeighborNumNeg;
	double dEnrich;

	/* open file */
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: MotifMap_ScanMultipleMatrix_Export, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	/* export motif statistics */
	fprintf(fpOut, "Motif");
	for(ni=0; ni<nMotifNum; ni++)
	{
		fprintf(fpOut, "\t%s", vMotifName[ni]->m_pString);
	}
	fprintf(fpOut, "\n");

	fprintf(fpOut, "Summary");
	for(ni=0; ni<nMotifNum; ni++)
	{
		nEffecLenPos = (int)(pEffecLenPos->pMatElement[ni]);
		nSiteNumPos = (int)(pSiteNumPos->pMatElement[ni]);
		if(nPosOnly == 0)
		{
			nEffecLenNeg = (int)(pEffecLenNeg->pMatElement[ni]);
			nSiteNumNeg = (int)(pSiteNumNeg->pMatElement[ni]);
			
			if( (nEffecLenPos == 0) || (nEffecLenNeg == 0) )
			{
				dEnrich = -1.0;
			}
			else if(nSiteNumNeg == 0)
			{
				dEnrich = (double)(nSiteNumPos*nEffecLenNeg)/((double)(nEffecLenPos)*((double)(nSiteNumNeg)+0.1));
			}
			else
			{
				dEnrich = (double)(nSiteNumPos*nEffecLenNeg)/(double)(nEffecLenPos*nSiteNumNeg);
			}
		}
		else
		{
			nEffecLenNeg = -1;
			nSiteNumNeg = -1;
			dEnrich = -1.0;
		}

		fprintf(fpOut, "\t%f/%d/%d/%d/%d", dEnrich, nSiteNumPos, nEffecLenPos, nSiteNumNeg, nEffecLenNeg);
	}
	fprintf(fpOut, "\n");
	fprintf(fpOut, "\n\n\n");

	/* export clustering statistics */
	fprintf(fpOut, "Motif");
	for(ni=0; ni<nMotifNum; ni++)
	{
		fprintf(fpOut, "\t%s", vMotifName[ni]->m_pString);
	}
	fprintf(fpOut, "\n\n");

	for(ni=0; ni<nMotifNum; ni++)
	{
		fprintf(fpOut, "%s", vMotifName[ni]->m_pString);
		for(nj=0; nj<nMotifNum; nj++)
		{
			nSiteNumPos = (int)(pSiteNumPos->pMatElement[ni]);
			nNeighborNumPos = (int)(DMGETAT(pNeighborNumPos, ni, nj));
			if(nPosOnly == 0)
			{
				nSiteNumNeg = (int)(pSiteNumNeg->pMatElement[ni]);
				nNeighborNumNeg = (int)(DMGETAT(pNeighborNumNeg, ni, nj));
				if( (nSiteNumPos == 0) || (nSiteNumNeg == 0) )
				{
					dEnrich = -1.0;
				}
				else if(nNeighborNumNeg == 0)
				{
					dEnrich = (double)(nNeighborNumPos*nSiteNumNeg)/((double)(nSiteNumPos)*((double)nNeighborNumNeg+0.1));
				}
				else
				{
					dEnrich = (double)(nNeighborNumPos*nSiteNumNeg)/(double)(nSiteNumPos*nNeighborNumNeg);
				}
			}
			else
			{
				nSiteNumNeg = -1;
				nNeighborNumNeg = -1;
				dEnrich = -1.0;
			}

			fprintf(fpOut, "\t%f/%d/%d/%d/%d", dEnrich, nSiteNumPos, nNeighborNumPos,  
				nSiteNumNeg, nNeighborNumNeg);
		}
		fprintf(fpOut, "\n");

		fprintf(fpOut, "%s", vMotifName[ni]->m_pString);
		for(nj=0; nj<nMotifNum; nj++)
		{
			nEffecLenNeg = (int)(pEffecLenPos->pMatElement[nj]);
			nSiteNumNeg = (int)(pSiteNumPos->pMatElement[nj]);
			nEffecLenPos = (int)(DMGETAT(pNeighborLenPos, ni, nj));
			nSiteNumPos = (int)(DMGETAT(pNeighborSiteNumPos, ni, nj));
				
			if( (nEffecLenPos == 0) || (nEffecLenNeg == 0) )
			{
				dEnrich = -1.0;
			}
			else if(nSiteNumNeg == 0)
			{
				dEnrich = (double)(nSiteNumPos*nEffecLenNeg)/((double)(nEffecLenPos)*((double)(nSiteNumNeg)+0.1));
			}
			else
			{
				dEnrich = (double)(nSiteNumPos*nEffecLenNeg)/(double)(nEffecLenPos*nSiteNumNeg);
			}
			
			fprintf(fpOut, "\t%f/%d/%d/%d/%d", dEnrich, nSiteNumPos, nEffecLenPos,  
				nSiteNumNeg, nEffecLenNeg);
		}
		fprintf(fpOut, "\n\n");
	}


	/* close file */
	fclose(fpOut);


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_CountKmer_Sequential_Main: count k-mers in a set of           */
/*  FASTA sequences.                                                       */
/*  dC  = the cutoff for average conservation scores.                      */
/* ----------------------------------------------------------------------- */ 
int MotifMap_CountKmer_Sequential_Main(char strInputPath[], char strSeqFile[],  
					char strOutputPath[], int nK, int nUseCS, 
					double dC, char strCSPrefix[])
{
	/* sequences and conservation */
	/* sequences */
	struct FLEXSEQMOTIF *pSeqMtf = NULL;
	/* sequence number */
	int nSeqCount = 0;
	/* site number */
	struct DOUBLEMATRIX *pKmerCount = NULL;
	int nKmerNum;
	int nBaseTypeNum = 4;
	/* total sequence length and total site number */
	double dTotal = 0.0;
	char strWord[LINE_LENGTH];
	int nTemp;
	int nWordId;
	
	/* motif consensus */
	FILE *fpIn;
	FILE *fpOut;
	char strLine[MED_LINE_LENGTH];
	
	
	/* load sequence */
	char strSeqLine[LONG_LINE_LENGTH];
	struct tagSequence *pNewSeq;
	char strSeqAlias[LINE_LENGTH];
	
	/* other variables */
	int ni,nj;
	
	/* #################################### */
	/* initialize                           */
	/* #################################### */
	if(nK <= 0)
		return PROC_SUCCESS;
	nKmerNum = (int)pow((double)nBaseTypeNum, (double)nK);
	
	/* #################################### */
	/* load initial parameters              */
	/* #################################### */
	AdjustDirectoryPath(strInputPath);

	pKmerCount = NULL;
	pKmerCount = CreateDoubleMatrix(1, nKmerNum);
	if(pKmerCount == NULL)
	{
		printf("Error: MotifMap_CountKmer_Main, cannot allocate memory!\n");
		exit(EXIT_FAILURE);
	}

	/* #################################### */
	/* load sequences and conservation      */
	/* #################################### */
	sprintf(strLine, "%s%s", strInputPath, strSeqFile);
	fpIn = NULL;
	fpIn = fopen(strLine, "r");
	if(fpIn == NULL)
	{
        printf("Error: MotifMap_CountKmer_Main, cannot open sequence file!\n");
		exit(EXIT_FAILURE);
	}

	/* read */
	pNewSeq = NULL;
	while(fgets(strSeqLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strSeqLine);
		StrTrimRight(strSeqLine);
		if(strSeqLine[0] == '\0')
			continue;

		if(strSeqLine[0] == '>')
		{
			if(pNewSeq != NULL)
			{
				/* #################################### */
				/* load sequence and score              */
				/* #################################### */
				pSeqMtf = NULL;
				pSeqMtf = FLEXSEQMTFLOADSEQ_FOR_SCANCONSENSUS_SINGLE(pNewSeq,
					strInputPath, nUseCS, strCSPrefix);
				SequenceDelete(pNewSeq);

				/* #################################### */
				/* count kmer                           */
				/* #################################### */
				MotifMap_CountKmer_In_FlexSeqMtf(pSeqMtf, 0,
						nUseCS, dC, nK, pKmerCount);
				
				/* #################################### */
				/* release memory                       */
				/* #################################### */
				FLEXSEQMTFDESTROY(pSeqMtf);
				nSeqCount++;
			}

			
			pNewSeq = NULL;
			pNewSeq = SequenceCreate();
			if(pNewSeq == NULL)
			{
				fclose(fpIn);
				printf("Error: MotifMap_CountKmer_Main, cannot allocate memory for loading sequences!\n");
				exit(EXIT_FAILURE);
			}
			pNewSeq->m_nIndex = nSeqCount;
			strcpy(pNewSeq->m_strAlias, (strSeqLine+1));
			strcpy(strSeqAlias, (strSeqLine+1));
		}
		else
		{
			if(SequenceAddTail(pNewSeq, strSeqLine) == PROC_FAILURE)
			{
				fclose(fpIn);
				printf("Error: MotifMap_CountKmer_Main, cannot allocate memory for loading sequences!\n");
				exit(EXIT_FAILURE);
			}			
		}
	}

	/* the last sequence */
	if(pNewSeq != NULL)
	{
		/* #################################### */
		/* load sequence and score              */
		/* #################################### */
		pSeqMtf = NULL;
		pSeqMtf = FLEXSEQMTFLOADSEQ_FOR_SCANCONSENSUS_SINGLE(pNewSeq,
			strInputPath, nUseCS, strCSPrefix);
		SequenceDelete(pNewSeq);

		/* #################################### */
		/* count kmer                           */
		/* #################################### */
		MotifMap_CountKmer_In_FlexSeqMtf(pSeqMtf, 0,
				nUseCS, dC, nK, pKmerCount);

		/* #################################### */
		/* release memory                       */
		/* #################################### */
		FLEXSEQMTFDESTROY(pSeqMtf);
		nSeqCount++;
	}

	
	/* close files */
	fclose(fpIn);
	
	/* #################################### */
	/* Export results                       */
	/* #################################### */
	fpOut = NULL;
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
        printf("Error: MotifMap_CountKmer_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	dTotal = 0.0;
	for(ni=0; ni<nKmerNum; ni++)
	{
		dTotal += pKmerCount->pMatElement[ni];

		nWordId = ni;
		for(nj=1; nj<nK; nj++)
		{
			nTemp = nWordId%nBaseTypeNum;
			switch(nTemp)
			{
				case 0: strWord[nK-nj] = 'A';
					break;
				case 1: strWord[nK-nj] = 'C';
					break;
				case 2: strWord[nK-nj] = 'G';
					break;
				case 3: strWord[nK-nj] = 'T';
					break;
			}
			nWordId -= nTemp;
			nWordId /= nBaseTypeNum;
		}

		nTemp = nWordId;
		switch(nTemp)
		{
			case 0: strWord[nK-nj] = 'A';
				break;
			case 1: strWord[nK-nj] = 'C';
				break;
			case 2: strWord[nK-nj] = 'G';
				break;
			case 3: strWord[nK-nj] = 'T';
				break;
		}

		strWord[nK] = '\0';

		fprintf(fpOut, "%s\t%d\n", strWord, (int)(pKmerCount->pMatElement[ni]));
	}

	fprintf(fpOut, "Total\t%d\n", (int)dTotal);

	fclose(fpOut);

	
	/* #################################### */
	/* Release Memory                       */
	/* #################################### */
	DestroyDoubleMatrix(pKmerCount);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_CountKmer_In_FlexSeqMtf: count kmer                           */
/* ----------------------------------------------------------------------- */ 
int MotifMap_CountKmer_In_FlexSeqMtf(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,
				int nUseCS, double dC, int nK, struct DOUBLEMATRIX *pKmerCount)
{
	/* define */
	int nBaseTypeNum = 4;
	int ni,nj,nLen;
	unsigned char *pBase;
	unsigned char *pCS;
	int nWordId = 0;
	int nBadNum = 0;
	int nScale = (int)pow((double)nBaseTypeNum, (double)(nK-1));
	int nOK;
	double dTotConserve,dAveConserve;

	/* init */
	nLen = pSeqMtf->vSeq[nIndex]->nWidth;
	if(nLen <= 0)
	{
		return PROC_SUCCESS;
	}
	pBase = pSeqMtf->vSeq[nIndex]->pMatElement;
	dTotConserve = 0.0;
	if(nUseCS == 1)
	{
		pCS = pSeqMtf->vScore[0]->pMatElement;
	}

	/* scan the first word */
	for(ni=0; ni<(nK-1); ni++)
	{
		if(pBase[ni] < nBaseTypeNum)
		{
			nWordId = nWordId*nBaseTypeNum+pBase[ni];
		}
		else
		{
			nWordId = nWordId*nBaseTypeNum;
			nBadNum += 1;
		}
		
		if(nUseCS == 1)
		{
			dTotConserve += pCS[ni];
		}
	}

	for(; ni<nLen; ni++)
	{
		if(pBase[ni] < nBaseTypeNum)
		{
			nWordId = nWordId*nBaseTypeNum+pBase[ni];
		}
		else
		{
			nWordId = nWordId*nBaseTypeNum;
			nBadNum += 1;
		}

		if(nUseCS == 1)
		{
			dTotConserve += pCS[ni];
		}

		if(nBadNum == 0)
		{
			nOK = 1;
			/* compute conservation */
			if(nUseCS == 1)
			{
				dAveConserve = dTotConserve/(double)nK;
				if(dAveConserve < dC)
				{
					nOK = 0;
				}
			}
			
			if(nOK == 1)
			{
				pKmerCount->pMatElement[nWordId] += 1.0;
			}
		}
		
		nj = ni-nK+1;
		if(pBase[nj] < nBaseTypeNum)
		{
			nWordId = nWordId-pBase[nj]*nScale;
		}
		else
		{
			nWordId = nWordId;
			nBadNum -= 1;
		}

		if(nUseCS == 1)
		{
			dTotConserve -= pCS[nj];
		}
	}


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_Rescore_Main: rescore the TFBS using the background model     */
/*  learnt from current sites.                                             */
/* ----------------------------------------------------------------------- */ 
int MotifMap_Rescore_Main(char strMotifPath[], char strInputPath[], char strOutputPath[])
{
	/* motif matrix */
	struct DOUBLEMATRIX *pMotif = NULL;
	struct DOUBLEMATRIX *pBoostBG = NULL;
		
	/* number of base types */
	int nBaseTypeNum = 4;
	char strSeqName[MED_LINE_LENGTH];
	char strChr[LINE_LENGTH];
	int nStart,nEnd;
	char chStrand;
	double dScore;
	double dNewScore;
	double dAffinity;
	char strSite[MED_LINE_LENGTH];
	char strOriSite[MED_LINE_LENGTH];
	double dTemp;
	
	/* motif consensus */
	FILE *fpIn;
	FILE *fpOut;
	char strLine[LONG_LINE_LENGTH];
	
	/* other variables */
	int ni,nj;
	double *pEle1,*pEle2;
	double dSum;
	
	/* #################################### */
	/* initialize                           */
	/* #################################### */
	
	/* motif PWM */
	pMotif = NULL;
	pMotif = DMLOAD(strMotifPath);
	if( pMotif == NULL)
	{
		printf("Error: MotifMap_Rescore_Main, cannot load init motif pseudocount!\n");
		exit(EXIT_FAILURE);
	}
	
	pBoostBG = NULL;
	pBoostBG = CreateDoubleMatrix(pMotif->nHeight, pMotif->nWidth);
	if(pBoostBG == NULL)
	{
		printf("Error: MotifMap_Rescore_Main, cannot create new background matrix!\n");
		exit(EXIT_FAILURE);
	}

	/* #################################### */
	/* load sites and train background      */
	/* #################################### */
	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: MotifMap_Rescore_Main, cannot open input files!\n");
		exit(EXIT_FAILURE);
	}
	
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;
		
		sscanf(strLine, "%s %s %d %d %c %lf %s", strSeqName,
			strChr, &nStart, &nEnd, &chStrand, &dScore, strSite);

		if((nEnd-nStart+1) != pMotif->nHeight)
		{
			printf("Error: MotifMap_Rescore_Main, site does not match motif length!\n");
			exit(EXIT_FAILURE);
		}

		StrMakeUpper(strSite);
		for(ni=0; ni<pBoostBG->nHeight; ni++)
		{
			switch(strSite[ni])
			{
				case 'A': dTemp = DMGETAT(pBoostBG, ni, 0)+1.0;
					DMSETAT(pBoostBG, ni, 0, dTemp);
					break;
				case 'C': dTemp = DMGETAT(pBoostBG, ni, 1)+1.0;
					DMSETAT(pBoostBG, ni, 1, dTemp);
					break;
				case 'G': dTemp = DMGETAT(pBoostBG, ni, 2)+1.0;
					DMSETAT(pBoostBG, ni, 2, dTemp);
					break;
				case 'T': dTemp = DMGETAT(pBoostBG, ni, 3)+1.0;
					DMSETAT(pBoostBG, ni, 3, dTemp);
					break;
			}
		}
	}
	
	fclose(fpIn);
	
	/* #################################### */
	/* save new background                  */
	/* #################################### */
	sprintf(strLine, "%s.bst", strMotifPath);
	DMSAVE(pBoostBG, strLine);
	
	pEle1 = pMotif->pMatElement;
	for(ni=0; ni<pMotif->nHeight; ni++)
	{
		pEle2 = pEle1;
		dSum = 0.0;
		for(nj=0; nj<pMotif->nWidth; nj++)
		{
			*pEle2 += 1e-3;
			dSum += (*pEle2);
			pEle2++;
		}

		for(nj=0; nj<pMotif->nWidth; nj++)
		{
			*pEle1 = (*pEle1/dSum);
			pEle1++;
		}
	}

	pEle1 = pBoostBG->pMatElement;
	for(ni=0; ni<pBoostBG->nHeight; ni++)
	{
		pEle2 = pEle1;
		dSum = 0.0;
		for(nj=0; nj<pBoostBG->nWidth; nj++)
		{
			*pEle2 += 1e-1;
			dSum += (*pEle2);
			pEle2++;
		}

		for(nj=0; nj<pBoostBG->nWidth; nj++)
		{
			*pEle1 = (*pEle1/dSum);
			pEle1++;
		}
	}

	sprintf(strLine, "%s.tnorm", strMotifPath);
	DMSAVE(pMotif, strLine);
	sprintf(strLine, "%s.bnorm", strMotifPath);
	DMSAVE(pBoostBG, strLine);

	DMLOGTS(pMotif);
	DMLOGTS(pBoostBG);

	/* #################################### */
	/* load sites and rescore them          */
	/* #################################### */
	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: MotifMap_Rescore_Main, cannot open input files!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: MotifMap_Rescore_Main, cannot open output files!\n");
		exit(EXIT_FAILURE);
	}
	
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;
		
		sscanf(strLine, "%s %s %d %d %c %lf %s", strSeqName,
			strChr, &nStart, &nEnd, &chStrand, &dScore, strSite);

		if((nEnd-nStart+1) != pMotif->nHeight)
		{
			printf("Error: MotifMap_Rescore_Main, site does not match motif length!\n");
			exit(EXIT_FAILURE);
		}

		strcpy(strOriSite, strSite);
		StrMakeUpper(strSite);
		dNewScore = 0.0;
		dAffinity = 0.0;
		for(ni=0; ni<pBoostBG->nHeight; ni++)
		{
			switch(strSite[ni])
			{
				case 'A': dNewScore += (DMGETAT(pMotif, ni, 0)-DMGETAT(pBoostBG, ni, 0));
					dAffinity += DMGETAT(pMotif, ni, 0);
					break;
				case 'C': dNewScore += (DMGETAT(pMotif, ni, 1)-DMGETAT(pBoostBG, ni, 1));
					dAffinity += DMGETAT(pMotif, ni, 1);
					break;
				case 'G': dNewScore += (DMGETAT(pMotif, ni, 2)-DMGETAT(pBoostBG, ni, 2));
					dAffinity += DMGETAT(pMotif, ni, 2);
					break;
				case 'T': dNewScore += (DMGETAT(pMotif, ni, 3)-DMGETAT(pBoostBG, ni, 3));
					dAffinity += DMGETAT(pMotif, ni, 3);
					break;
			}
		}

		dNewScore /= log(10.0);
		dAffinity /= log(10.0);

		fprintf(fpOut, "%s\t%s\t%d\t%d\t%c\t%f\t%s\t%f\t%f\n", strSeqName, 
			strChr, nStart, nEnd, chStrand, dNewScore, strOriSite, dScore, dAffinity);
	}
	
	fclose(fpIn);
	fclose(fpOut);

	/* #################################### */
	/* release memory                       */
	/* #################################### */
	DestroyDoubleMatrix(pMotif);
	DestroyDoubleMatrix(pBoostBG);
	
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_FilterOverlappingSite:                                        */
/*  Filter overlapping sites. Return the number of non-overlapping sites.  */
/* ----------------------------------------------------------------------- */ 
int MotifMap_FilterOverlappingSite(char strInFile[], char strOutFile[])
{
	/* define */
	char strLine[LONG_LINE_LENGTH];
	char strName[MED_LINE_LENGTH];
	char strChr[MED_LINE_LENGTH];
	char strChr0[MED_LINE_LENGTH];
	int nStart,nStart0;
	int nEnd,nEnd0;
	FILE *fpIn;
	FILE *fpOut;
	int nNRSiteNum = 0;
	int nMaxStart,nMinEnd;

	/* open files */
	fpIn = NULL;
	fpIn = fopen(strInFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: MotifMap_FilterOverlappingSite, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}
	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: MotifMap_FilterOverlappingSite, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	/* process line by line */
	strcpy(strChr0, "");
	nStart0 = -1;
	nEnd0 = -1;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
		{
			fprintf(fpOut, "%s\n", strLine);
			continue;
		}

		sscanf(strLine, "%s %s %d %d", strName, strChr, &nStart, &nEnd);
		if(nNRSiteNum == 0)
		{
			strcpy(strChr0, strChr);
			nStart0 = nStart;
			nEnd0 = nEnd;
			nNRSiteNum += 1;
			fprintf(fpOut, "%s\n", strLine);
		}
		else
		{
			if(strcmp(strChr0, strChr)==0)
			{
				if(nStart > nStart0)
					nMaxStart = nStart;
				else
					nMaxStart = nStart0;

				if(nEnd < nEnd0)
					nMinEnd = nEnd;
				else
					nMinEnd = nEnd0;

				if(nMaxStart <= nMinEnd)
				{
				}
				else
				{
					nStart0 = nStart;
					nEnd0 = nEnd;
					nNRSiteNum += 1;
					fprintf(fpOut, "%s\n", strLine);
				}

			}
			else
			{
				strcpy(strChr0, strChr);
				nStart0 = nStart;
				nEnd0 = nEnd;
				nNRSiteNum += 1;
				fprintf(fpOut, "%s\n", strLine);
			}
		}
	}

	/* close files */
	fclose(fpOut);
	fclose(fpIn);

	/* return */
	return nNRSiteNum;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifSampler_Quality: Motif sampler based on two stage A/B quality     */
/*  control.                                                               */
/* ----------------------------------------------------------------------- */ 
int MotifSampler_Quality(char strInFile[], char strOutFile[], int nIteration, 
						   int nBGOrder, struct DOUBLEMATRIX *pAFreqPrior,
						   int nMotifNum, 
						   struct DOUBLEMATRIX *pMatrixPrior[],
						   int nSeedNum,
						   struct DOUBLEMATRIX *pSeedMotif[])
{
	/* define */
	/* all sequences */
	struct SEQMOTIF **pSeqMtf;
	/* background */
	struct DOUBLEMATRIX *pBG;
	/* log background */
	struct DOUBLEMATRIX *pLogBG;
	/* background0 */
	struct DOUBLEMATRIX *pBG0;
	/* log background0 */
	struct DOUBLEMATRIX *pLogBG0;
	/* sequence number */
	int nCount;
	/* motif matrix */
	struct MOTIFMATRIX *vMotif[MAX_MOTIF_NUM];
	/* motif length */
	int vMotifLen[MAX_MOTIF_NUM];
	int nMaxMotifLen;

	/* background base transition */
	int nScale;

	/* frequency count */
	double vFreqCount[MAX_TYPE_NUM][2];

	/* number of base types */
	int nBaseTypeNum;
	
	/* others */
	int ni,nj,niter;
	double *pElement;
	/*double *pScoreEle;*/
	/*FILE *fpSeg;*/
	int nRecordNum,nIsRecording;
	/* char strMotifFile[255]; */

	/* initial check */
	nBaseTypeNum = 4;
	nScale = (int)(pow((double)nBaseTypeNum, (double)(nBGOrder-1)));
	
	if(pAFreqPrior == NULL)
	{
		printf("Error: MotifSampler_Quality, prior frequency of observing motif matches must be specified!\n");
		exit(EXIT_FAILURE);
	}
	if( (pAFreqPrior->nHeight != (nMotifNum+1)) || (pAFreqPrior->nWidth != 2) )
	{
		printf("Error: MotifSampler_Quality, prior frequency of observing motif and motif number not match!\n");
		exit(EXIT_FAILURE);
	}
	
	for(ni=0; ni<nMotifNum; ni++)
	{
		if(pMatrixPrior[ni] == NULL)
		{
			printf("Error: MotifSampler_Quality, prior motif matrix must be specified!\n");
			exit(EXIT_FAILURE);
		}
		if(pMatrixPrior[ni]->nWidth != nBaseTypeNum)
		{
			printf("Error: MotifSampler_Quality, priors not properly specified!\n");
			exit(EXIT_FAILURE);
		}
	}
	if(nSeedNum>0)
	{
		for(ni=0; ni<nSeedNum; ni++)
		{
			if(pSeedMotif[ni] == NULL)
			{
				printf("Error: MotifSampler_Quality, prior seed motif matrix must be specified!\n");
				exit(EXIT_FAILURE);
			}
			
			if(pSeedMotif[ni]->nWidth != nBaseTypeNum)
			{
				printf("Error: MotifSampler_Quality, priors not properly specified!\n");
				exit(EXIT_FAILURE);
			}
		}
	}

	/* init prior probability */
	for(ni=0; ni<MAX_TYPE_NUM; ni++)
	{
		for(nj=0; nj<2; nj++)
		{
			vFreqCount[ni][nj] = 0.0;
		}
	}
	for(ni=0; ni<pAFreqPrior->nHeight; ni++)
	{
		for(nj=0; nj<pAFreqPrior->nWidth; nj++)
		{
			vFreqCount[ni][nj] = DMGETAT(pAFreqPrior, ni, nj);
		}
	}
	
	/* init motif */
	for(ni=0; ni<MAX_MOTIF_NUM; ni++)
	{
		vMotif[ni] = NULL;
	}

	/* prepare motif matrix */
	nMaxMotifLen = 0;
	for(ni=0; ni<nMotifNum; ni++)
	{
		/* create */
		vMotif[ni] = MTFMCREATE();
		if(vMotif[ni] == NULL)
		{
			printf("Error: MotifSampler_Quality, cannot create motif matrix!\n");
			exit(EXIT_FAILURE);
		}

		/* init */
		vMotif[ni]->nMotifType = ni;
		vMotif[ni]->pCOUNTp = DMCLONE(pMatrixPrior[ni]);
		vMotifLen[ni] = vMotif[ni]->pCOUNTp->nHeight;
		if(vMotifLen[ni] > nMaxMotifLen)
			nMaxMotifLen = vMotifLen[ni];
	}

	/* load sequences */
	pSeqMtf = NULL;
	pSeqMtf = SEQMTFLOADPRIMARYSEQFROMFASTA_FORQUALITY(&nCount, strInFile);

	/* calculate background */
	pBG = NULL;
	pBG = SEQMTFESTIMATENUCLEICBGMC(pSeqMtf, nCount, 0, nBGOrder);
	pLogBG = NULL;
	pLogBG = DMCLONE(pBG);
	pElement = pLogBG->pMatElement;
	for(ni=0; ni<pLogBG->nHeight; ni++)
	{
		for(nj=0; nj<pLogBG->nWidth; nj++)
		{
			*pElement = log(*pElement);
			pElement++;
		}
	}
	/* DMSAVE(pBG, "bgestimate.txt");*/
	pBG0 = NULL;
	pBG0 = SEQMTFESTIMATENUCLEICBGMC(pSeqMtf, nCount, 0, 0);
	pLogBG0 = NULL;
	pLogBG0 = DMCLONE(pBG0);
	pElement = pLogBG0->pMatElement;
	for(ni=0; ni<pLogBG0->nHeight; ni++)
	{
		for(nj=0; nj<pLogBG0->nWidth; nj++)
		{
			*pElement = log(*pElement);
			pElement++;
		}
	}

	/* get init state */
	for(ni=0; ni<nCount; ni++)
	{
		Quality_RandInit(pSeqMtf[ni], nBGOrder, pBG, nMotifNum, vMotif, vMotifLen, nMaxMotifLen,
					nSeedNum, pSeedMotif, vFreqCount);
		Quality_UpdateNeighborCount(pSeqMtf[ni]);
	}
	/* update motif matrix */
	for(ni=0; ni<nMotifNum; ni++)
	{
		MTFMREFRESH_QUALITY(vMotif[ni]);
		/*sprintf(strMotifFile, "motif%d.txt", ni);
		DMSAVE(vMotif[ni]->pCOUNTp, strMotifFile);*/
	}

	/* quality sampler iterations */
	nRecordNum = MODULE_SAMPLE_NUM;
	if((int)(nIteration/2) < nRecordNum)
		nRecordNum = (int)(nIteration/2);
	for(niter=0; niter<nIteration; niter++)
	{
		if(niter%50 == 0)
		{
			printf("iter=%d...\n", niter);
		}

		if((niter+nRecordNum) < nIteration)
			nIsRecording = 0;
		else
			nIsRecording = 1;

		/* sequence-wise quality sampling */
		for(ni=0; ni<nCount; ni++)
		{
			Quality_Sampler(pSeqMtf[ni], vFreqCount,
						 nBGOrder, pLogBG, pLogBG0,
						 nMotifNum, vMotifLen, nMaxMotifLen, 
						 vMotif, niter, nIteration);
		}

		/* TODO: local mode shifting */

		/* TODO: change motif length */
	}

	
	/* TODO: report result */
	for(ni=0; ni<nCount; ni++)
	{
		Quality_CallMotifSites(pSeqMtf[ni], vMotifLen);
	}
	WriteMotifToFile(strOutFile, nMotifNum, vMotif, pBG0, nCount, pSeqMtf);
	
	/* destroy */
	DestroyDoubleMatrix(pBG);
	DestroyDoubleMatrix(pLogBG);
	DestroyDoubleMatrix(pBG0);
	DestroyDoubleMatrix(pLogBG0);
	for(ni=0; ni<nCount; ni++)
	{
		SEQMTFDESTROY(pSeqMtf[ni]);
		pSeqMtf[ni] = NULL;
	}
	free(pSeqMtf);
	for(ni=0; ni<nMotifNum; ni++)
	{
		MTFMDESTROY(vMotif[ni]);
		vMotif[ni] = NULL;
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Quality_RandInit: Random init.                                         */
/* ----------------------------------------------------------------------- */ 
int Quality_RandInit(struct SEQMOTIF *pSeqMtf, int nBGOrder, 
					 struct DOUBLEMATRIX *pBG, 
					 int nMotifNum, struct MOTIFMATRIX *vMotif[], 
					 int vMotifLen[], int nMaxMotifLen, int nSeedNum, 
					 struct DOUBLEMATRIX *pSeedMotif[], 
					 double vFreqCount[MAX_TYPE_NUM][2])
{
	/* define */
	double vFreqPrior[MAX_TYPE_NUM][2];
	double vFreqPost[MAX_TYPE_NUM];
	int ni,nj,nk,nx,nLen;
	unsigned char *pAi,*pBi,*pBase;
	double dLike;
	int nOK;
	int nPreNum;
	double dTotal;
	double dRand;
	double dPriorTotal;
	int nMask;

	/* check */
	if(nSeedNum > nMotifNum)
	{
		printf("Error: motif seeds should be <= motif number!\n");
		exit(EXIT_FAILURE);
	}

	dPriorTotal = vFreqCount[0][0];
	for(ni=1; ni<=nMotifNum; ni++)
	{
		dPriorTotal += vFreqCount[ni][0]+vFreqCount[ni][1];
	}
	vFreqPrior[0][0] = vFreqCount[0][0]/dPriorTotal;
	for(ni=1; ni<=nMotifNum; ni++)
	{
		for(nj=0; nj<2; nj++)
		{
			vFreqPrior[ni][nj] = vFreqCount[ni][nj]/dPriorTotal;
		}
	}

	/* init */
	nLen = pSeqMtf->ppSeq[0]->nWidth;
	pBase = pSeqMtf->ppSeq[0]->pMatElement;
	for(ni=0; ni<nLen; ni++)
	{
		if(ni < nBGOrder)
			nPreNum = ni;
		else
			nPreNum = nBGOrder;
		

		if(pBase[ni] >= pBG->nWidth)
			continue;

		nMask = 0;
		for(nj=0; nj<nMaxMotifLen; nj++)
		{
			nk = ni+nj;
			if(nk >= nLen)
				break;
			if(pBase[nk] >= pBG->nWidth)
			{
				nMask = 1;
				break;
			}
		}

		if(nMask == 1)
		{
			vFreqCount[0][0] += 1.0;
			continue;
		}


		nOK = 0;
		/* sample with motifs with information */
		for(nj=0; nj<nSeedNum; nj++)
		{
			if(ni+pSeedMotif[nj]->nHeight-1 >= nLen)
			{
				vFreqPost[0] = 1.0;
				vFreqPost[1] = 0.0;
				vFreqPost[2] = 0.0;
			}
			else
			{
				dTotal = 0.0;
				dLike = BaseLikelihood_Background((pBase+ni), pSeedMotif[nj]->nHeight, pBG, nBGOrder, nPreNum, '+');
				vFreqPost[0] = vFreqPrior[0][0]*dLike;
				dTotal += vFreqPost[0];
				dLike = BaseLikelihood_Motif((pBase+ni), pSeedMotif[nj]->nHeight, pSeedMotif[nj], '+');
				vFreqPost[1] = vFreqPrior[nj+1][0]*dLike;
				dTotal += vFreqPost[1];
				dLike = BaseLikelihood_Motif((pBase+ni), pSeedMotif[nj]->nHeight, pSeedMotif[nj], '-');
				vFreqPost[2] = vFreqPrior[nj+1][1]*dLike;
				dTotal += vFreqPost[2];

				for(nx=0; nx<3; nx++)
				{
					vFreqPost[nx] = vFreqPost[nx]/dTotal;
				}
			}

			dRand = rand_u();
			dTotal = 0.0;
			for(nx=0; nx<3; nx++)
			{
				dTotal += vFreqPost[nx];
				if(dRand <= dTotal)
					break;
			}
			if(nx>=3) nx--;
			if(nx>0)
			{
				nOK = 1;
				pAi = pSeqMtf->ppSamplePath[0]->pMatElement+ni;
				pBi = pSeqMtf->ppSamplePath[1]->pMatElement+ni;
				*pAi = (nj+1)*2+nx-1;
				*pBi = 1;
				if(nx == 1)
				{
					AddNucleicMatrixCount(vMotif[nj]->pCOUNTp, '+', (pBase+ni), vMotifLen[nj]);
					vFreqCount[nj+1][0] += 1.0;
				}
				else
				{
					AddNucleicMatrixCount(vMotif[nj]->pCOUNTp, '-', (pBase+ni), vMotifLen[nj]);
					vFreqCount[nj+1][1] += 1.0;
				}
				break;
			}
		}

		if(nOK == 1)
		{
			continue;
		}

		/* sample with motifs without information */
		for(; nj<nMotifNum; nj++)
		{
			if(ni+vMotifLen[nj]-1 >= nLen)
			{
				vFreqPost[0] = 1.0;
				vFreqPost[1] = 0.0;
				vFreqPost[2] = 0.0;
			}
			else
			{
				vFreqPost[0] = (1.0-vFreqPrior[nj+1][0]-vFreqPrior[nj+1][1]);
				vFreqPost[1] = vFreqPrior[nj+1][0];
				vFreqPost[2] = vFreqPrior[nj+1][1];
			}

			dRand = rand_u();
			dTotal = 0.0;
			for(nx=0; nx<3; nx++)
			{
				dTotal += vFreqPost[nx];
				if(dRand <= dTotal)
					break;
			}
			if(nx>=3) nx--;
			if(nx>0)
			{
				nOK = 1;
				pAi = pSeqMtf->ppSamplePath[0]->pMatElement+ni;
				pBi = pSeqMtf->ppSamplePath[1]->pMatElement+ni;
				*pAi = (nj+1)*2+nx-1;
				*pBi = 1;
				if(nx == 1)
				{
					AddNucleicMatrixCount(vMotif[nj]->pCOUNTp, '+', (pBase+ni), vMotifLen[nj]);
					vFreqCount[nj+1][0] += 1.0;
				}
				else
				{
					AddNucleicMatrixCount(vMotif[nj]->pCOUNTp, '-', (pBase+ni), vMotifLen[nj]);
					vFreqCount[nj+1][1] += 1.0;
				}
				break;
			}
		}

		if(nOK == 0)
			vFreqCount[0][0] += 1.0;
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Quality_UpdateNeighborCount: Update neighboring count                  */
/* ----------------------------------------------------------------------- */ 
int Quality_UpdateNeighborCount(struct SEQMOTIF *pSeqMtf)
{
	/* define */
	int ni,nj,nLen,nCount;
	unsigned char *pCi,*pAi,*pSi;
	int nStart,nEnd;

	/* process */
	nLen = pSeqMtf->ppSamplePath[2]->nWidth;
	pAi = pSeqMtf->ppSamplePath[0]->pMatElement;
	pSi = pSeqMtf->ppSamplePath[0]->pMatElement;
	for(ni=0; ni<nLen; ni++)
	{
		if(*pAi > 0)
		{
			nStart = ni-QUALITY_HALFWINDOW;
			nEnd = ni+QUALITY_HALFWINDOW;
			if(nStart < 0)
				nStart = 0;
			if(nEnd >= nLen)
				nEnd = nLen-1;

			nCount = 0;
			for(nj=nStart; nj<=nEnd; nj++)
			{
				if(pSi[nj] > 0)
					nCount++;
			}
			pCi = pSeqMtf->ppSamplePath[2]->pMatElement+ni;
			*pCi = (unsigned char)nCount;
		}

		/* get next */
		pAi++;
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Quality_Sampler: Quality sampler for one sequence.                     */
/* ----------------------------------------------------------------------- */ 
int Quality_Sampler(struct SEQMOTIF *pSeqMtf, double vFreqCount[MAX_TYPE_NUM][2],
					int nBGOrder, struct DOUBLEMATRIX *pLogBG, struct DOUBLEMATRIX *pLogBG0,
					int nMotifNum, int vMotifLen[], int nMaxMotifLen,
					struct MOTIFMATRIX *vMotif[],
					int niter, int nIteration)
{
	/* define */
	int nLen;
	unsigned char nBaseTypeNum;
	unsigned char *pAi,*pBi,*pCi,*pDi,*pBase;
	int ni,nj,nk,nx;
	int nPreNum;
	double dRand;

	/* background/motif likelihood */
	double vLike[MAX_TYPE_NUM][2];
	/* posteror probability */
	double vPost[MAX_TYPE_NUM][2];
	/* probability of true positive */
	double dPBi,dPBj;
	/* likelihood for neighbors Ai */
	double dLn0,dLn1;
	/* # of neighbors */
	int nNeighborCount;
	/* mask: if the site should be ignored for likelihood calculation */
	int nMask,nSiteMask;
	double dMax,dMin,dTotal;
	double dCol[2];

	int nMtfId,nMtfSt;
	int nEndMaxMotifLen;


	/* neighboring site list */
	struct SITEPOS *pNeighbors;
	struct SITEPOS *pNH,*pNT;
	struct SITEPOS *pNPrev,*pNNext;
	int nNstart,nNend;
	double dTemp;

	/* init */
	nBaseTypeNum = (unsigned char)(pLogBG->nWidth);
	nLen = pSeqMtf->ppSeq[0]->nWidth;
	pNeighbors = NULL;
	pNH = NULL;
	pNT = NULL;
	
	/* prepare neighbor list */
	pAi = pSeqMtf->ppSamplePath[0]->pMatElement;
	for(ni=0; ni<QUALITY_HALFWINDOW; ni++)
	{
		if(*pAi > 0)
		{
			/* create new site */
			pNH = NULL;
			pNH = SITEPOSCREATE();
			if(pNH == NULL)
			{
				printf("Error: cannot create neighboring site!\n");
				exit(EXIT_FAILURE);
			}
			pNH->nPos = ni;

			/* add to list */
			if(pNeighbors == NULL)
			{
				pNeighbors = pNH;
				pNT = pNH;
			}
			else
			{
				pNT->pNext = pNH;
				pNT = pNT->pNext;
			}
		}

		/* get next */
		pAi++;
	}

	/* sample base by base */
	pAi = pSeqMtf->ppSamplePath[0]->pMatElement;
	pBase = pSeqMtf->ppSeq[0]->pMatElement;
	for(ni=0; ni<nLen; ni++)
	{
		if(ni == 406 || ni == 454 || ni == 486)
		{
			dTemp = ni;
		}
		if(ni%50 == 0)
		{
			dTemp = (double)ni/50;
		}

		/* if 'N', skip */
		if((*pBase) >= nBaseTypeNum)
		{
			/* get next pos */
			pAi++;
			pBase++;
			continue;
		}

		
		/* -------------------- */
		/* get surrounding Ai's */
		/* -------------------- */

		/* remove head */
		nNstart = ni-QUALITY_HALFWINDOW;
		while(pNeighbors != NULL)
		{
			if(pNeighbors->nPos >= nNstart)
				break;

			pNH = pNeighbors;
			pNeighbors = pNeighbors->pNext;
			SITEPOSDESTROY(pNH);
		}
		if(pNeighbors == NULL)
			pNT = NULL;
		else
		{
			pNT = pNeighbors;
			while(pNT->pNext != NULL)
				pNT = pNT->pNext;
		}
		/* add tail */
		nNend = ni+QUALITY_HALFWINDOW;
		if(nNend < nLen)
		{
			if( (*(pSeqMtf->ppSamplePath[0]->pMatElement+nNend)) > 0 )
			{
				/* create new site */
				pNH = NULL;
				pNH = SITEPOSCREATE();
				if(pNH == NULL)
				{
					printf("Error: cannot create neighboring site!\n");
					exit(EXIT_FAILURE);
				}
				pNH->nPos = nNend;

				/* add to list */
				if(pNeighbors == NULL)
				{
					pNeighbors = pNH;
					pNT = pNH;
				}
				else
				{
					pNT->pNext = pNH;
					pNT = pNT->pNext;
				}
			}
		}

		/* -------------------- */
		/* clear current status */
		/* -------------------- */
		nk = (int)(*pAi);
		if(nk == 0)
		{
			vFreqCount[nk][0] -= 1.0;
		}
		else
		{
			/* update Ai parameters */
			nMtfId = (int)(nk/2);
			nMtfSt = nk%2;
			vFreqCount[nMtfId][nMtfSt] -= 1.0;
			nMtfId -= 1;

			/* if Bi=1, update motif parameters */
			pBi = pSeqMtf->ppSamplePath[1]->pMatElement+ni;
			if( (*pBi) == 1 )
			{
				if(nMtfSt == 0)
				{
					SubtractNucleicMatrixCount(vMotif[nMtfId]->pCOUNTp, '+', pBase, vMotifLen[nMtfId]);
					MTFMREFRESH_QUALITY(vMotif[nMtfId]);
				}
				else
				{
					SubtractNucleicMatrixCount(vMotif[nMtfId]->pCOUNTp, '-', pBase, vMotifLen[nMtfId]);
					MTFMREFRESH_QUALITY(vMotif[nMtfId]);
				}
			}

			/* update surrounding Ai's */
			pNPrev = NULL;
			pNNext = pNeighbors;
			while(pNNext != NULL)
			{
				nx = pNNext->nPos;
				if(nx != ni)
				{
					pCi = pSeqMtf->ppSamplePath[2]->pMatElement+nx;
					*pCi = (*pCi)-1;
					/* get next */
					pNPrev = pNNext;
					pNNext = pNNext->pNext;
				}
				else
				{
					if(pNPrev == NULL)
					{
						pNPrev = pNeighbors;
						pNeighbors = pNeighbors->pNext;
						SITEPOSDESTROY(pNPrev);
						pNPrev = NULL;
						pNNext = pNeighbors;
					}
					else
					{
						pNPrev->pNext = pNNext->pNext;
						SITEPOSDESTROY(pNNext);
						pNNext = pNPrev->pNext;
					}
					pCi = pSeqMtf->ppSamplePath[2]->pMatElement+ni;
					*pCi = (*pCi)-1;
				}
			}
		}

		
		/* -------------------------- */
		/* get predictive probability */
		/* -------------------------- */

		*pAi = 0;
		pBi = pSeqMtf->ppSamplePath[1]->pMatElement+ni;
		*pBi = 0;
		nNeighborCount = SITEPOSLISTGETCOUNT(pNeighbors);
		pCi = pSeqMtf->ppSamplePath[2]->pMatElement+ni;
		if(nk > 0)
		{
			if((int)(*pCi) != nNeighborCount)
			{
				printf("Error: neighbor count not match!\n");
				exit(EXIT_FAILURE);
			}
			*pCi = 0;
		}
		

		/* likelihood */
		if(ni<nBGOrder)
			nPreNum = ni;
		else
			nPreNum = nBGOrder;

		if((ni+nMaxMotifLen-1) < nLen)
		{
			vLike[0][0] = Quality_BaseLogLikelihood_Background(pBase, nMaxMotifLen, pLogBG, nBGOrder, pLogBG0, nPreNum, '+', &nMask);
			nSiteMask = nMask;
			for(nk=0; nk<nMotifNum; nk++)
			{
				vLike[nk+1][0] = Quality_BaseLogLikelihood_Motif(pBase, nMaxMotifLen, vMotifLen[nk], vMotif[nk]->pPWMp, pLogBG, nBGOrder, pLogBG0, '+', &nMask);
				vLike[nk+1][1] = Quality_BaseLogLikelihood_Motif(pBase, nMaxMotifLen, vMotifLen[nk], vMotif[nk]->pPWMp, pLogBG, nBGOrder, pLogBG0, '-', &nMask);
			}
		}
		else
		{
			nEndMaxMotifLen = nLen-ni;
			vLike[0][0] = Quality_BaseLogLikelihood_Background(pBase, nEndMaxMotifLen, pLogBG, nBGOrder, pLogBG0, nPreNum, '+', &nMask);
			nSiteMask = nMask;
			for(nk=0; nk<nMotifNum; nk++)
			{
				if(vMotifLen[nk] <= nEndMaxMotifLen)
				{
					vLike[nk+1][0] = Quality_BaseLogLikelihood_Motif(pBase, nEndMaxMotifLen, vMotifLen[nk], vMotif[nk]->pPWMp, pLogBG, nBGOrder, pLogBG0, '+', &nMask);
					vLike[nk+1][1] = Quality_BaseLogLikelihood_Motif(pBase, nEndMaxMotifLen, vMotifLen[nk], vMotif[nk]->pPWMp, pLogBG, nBGOrder, pLogBG0, '-', &nMask);			
				}
				else
				{
					vLike[nk+1][0] = vLike[0][0]-QUALITY_OVERLAP_PENALTY;
					vLike[nk+1][1] = vLike[0][0]-QUALITY_OVERLAP_PENALTY;
				}
			}
		}
		
		/* neighbor likelihood */
		Quality_Neighbor_Likelihood(pNeighbors, pSeqMtf->ppSamplePath[2]->pMatElement, 
			pSeqMtf->ppSamplePath[1]->pMatElement, &dLn0, &dLn1, niter, nIteration);
		
		/* combine */
		dPBi = Quality_EstimateTruePosProb(nNeighborCount+1, niter, nIteration);
		dPBj = log(1.0-dPBi);
		dPBi = log(dPBi);
		vPost[0][0] = vLike[0][0]+log(vFreqCount[0][0])+dLn0;
		for(nk=0; nk<nMotifNum; nk++)
		{
			for(nx=0; nx<2; nx++)
			{
				dCol[0] = dPBi+vLike[nk+1][nx];
				dCol[1] = dPBj+vLike[0][0];
				if(dCol[0] > dCol[1])
					dMax = dCol[0];
				else
					dMax = dCol[1];
				vPost[nk+1][nx] = dMax+log(exp(dCol[0]-dMax)+exp(dCol[1]-dMax));
				vPost[nk+1][nx] += log(vFreqCount[nk+1][nx])+dLn1;
			}
		}

		/* if there are any 'N', mask the site */
		if(nSiteMask == 1)
		{
			for(nk=0; nk<=nMotifNum; nk++)
			{
				for(nx=0; nx<2; nx++)
				{
					vPost[nk][nx] = 0.0;
				}
			}
			vPost[0][0] = 1.0;
		}
		else
		{
			/* check overlap */
			for(nk=0; nk<nMotifNum; nk++)
			{
				nMask = 0;
				for(nx=0; nx<vMotifLen[nk]; nx++)
				{
					if((ni-nx-1) >= 0)
					{
						if(*(pAi-nx-1) > 0)
						{
							nMask = 1;
							break;
						}
					}
				}
				if(nMask == 1)
					break;
			}

			if(nMask == 1)
			{
				for(nk=0; nk<=nMotifNum; nk++)
				{
					for(nx=0; nx<2; nx++)
					{
						vPost[nk][nx] = 0.0;
					}
				}
				vPost[0][0] = 1.0;
			}
	
			/* other wise get the posterior probability */
			else
			{
				dMin = vPost[0][0];
				for(nk=0; nk<nMotifNum; nk++)
				{
					for(nx=0; nx<2; nx++)
					{
						if(vPost[nk+1][nx] < dMin)
							dMin = vPost[nk+1][nx];
					}
				}

				/* see if there are some overlapping */
				for(nk=0; nk<nMotifNum; nk++)
				{
					nMask = 0;
					pDi = pAi;
					for(nx=0; nx<vMotifLen[nk]; nx++)
					{
						if((ni+nx) >= nLen)
						{
							nMask = 1;
							break;
						}

						if(*pDi > 0)
						{
							nMask = 1;
							break;
						}
						pDi++;
					}

					if(nMask == 1)
					{
						vPost[nk+1][0] = dMin-QUALITY_OVERLAP_PENALTY;
						vPost[nk+1][1] = dMin-QUALITY_OVERLAP_PENALTY;
					}
				}

				dMax = vPost[0][0]; 
				for(nk=0; nk<nMotifNum; nk++)
				{
					for(nx=0; nx<2; nx++)
					{
						if(vPost[nk+1][nx] > dMax)
							dMax = vPost[nk+1][nx];
					}
				}

				vPost[0][0] -= dMax;
				dTotal = exp(vPost[0][0]);
				for(nk=0; nk<nMotifNum; nk++)
				{
					for(nx=0; nx<2; nx++)
					{
						vPost[nk+1][nx] -= dMax;
						dTotal += exp(vPost[nk+1][nx]);
					}
				}
				dTotal = log(dTotal);

				vPost[0][0] = exp(vPost[0][0]-dTotal);
				vPost[0][1] = 0.0;
				for(nk=0; nk<nMotifNum; nk++)
				{
					for(nx=0; nx<2; nx++)
					{
						vPost[nk+1][nx] = exp(vPost[nk+1][nx]-dTotal);
					}
				}
			}
		}

		/* ------------- */
		/* sample status */
		/* ------------- */
		dRand = rand_u();
		dTotal = 0.0;
		for(nk=0; nk<=nMotifNum; nk++)
		{
			for(nx=0; nx<2; nx++)
			{
				dTotal += vPost[nk][nx];
				if(dRand <= dTotal)
					break;
			}

			if(dRand <= dTotal)
				break;
		}
				
		if(nk>nMotifNum || nx>=2)
		{
			printf("Error: floating point error!\n");
			exit(EXIT_FAILURE);
		}

		/* ----------------- */
		/* update parameters */
		/* ----------------- */
		vFreqCount[nk][nx] += 1.0;
		pBi = pSeqMtf->ppSamplePath[1]->pMatElement+ni;
		pCi = pSeqMtf->ppSamplePath[2]->pMatElement+ni;
		if(nk == 0)
		{
			*pAi = 0;
			*pBi = 0;
			*pCi = 0;		
		}
		else
		{
			/* modify current base status */
			*pAi = (unsigned char)(2*nk+nx);
			dPBi = Quality_EstimateTruePosProb(nNeighborCount+1, niter, nIteration);

			dCol[0] = log(dPBi)+vLike[nk][nx];
			dCol[1] = log(1.0-dPBi)+vLike[0][0];

			if(dCol[0] > dCol[1])
				dMax = dCol[0];
			else
				dMax = dCol[1];

			dTotal = log(exp(dCol[0]-dMax)+exp(dCol[1]-dMax));
			dCol[0] = dCol[0]-dMax-dTotal;
			dCol[0] = exp(dCol[0]);
			dCol[1] = 1.0-dCol[0];

			dRand = rand_u();
			if(dRand <= dCol[0])
			{
				*pBi = 1;
				/* add initial TFBS matrix pseudocount */
				switch(nx)
				{
					case 0: AddNucleicMatrixCount(vMotif[nk-1]->pCOUNTp, '+', pBase, vMotifLen[nk-1]);
						MTFMREFRESH_QUALITY(vMotif[nk-1]);
						break;
					case 1:	AddNucleicMatrixCount(vMotif[nk-1]->pCOUNTp, '-', pBase, vMotifLen[nk-1]);
						MTFMREFRESH_QUALITY(vMotif[nk-1]);
						break;
					default: 
						break;
				}
			}
			else
			{
				*pBi = 0;
			}

			/* add sitepos */
			pNH = NULL;
			pNH = SITEPOSCREATE();
			if(pNH == NULL)
			{
				printf("Error: cannot create neighboring site!\n");
				exit(EXIT_FAILURE);
			}
			pNH->nPos = ni;

			if(pNeighbors == NULL)
			{
				pNeighbors = pNH;
			}
			else
			{
				pNPrev = NULL;
				pNNext = pNeighbors;
				while(pNNext != NULL)
				{
					if(pNNext->nPos > ni)
					{
						break;
					}
					/* get next */
					pNPrev = pNNext;
					pNNext = pNNext->pNext;
				}
				if(pNPrev == NULL)
				{
					pNH->pNext = pNeighbors;
					pNeighbors = pNH;
				}
				else
				{
					pNH->pNext = pNNext;
					pNPrev->pNext = pNH;
				}
			}	
		
			/* modify neighbors' counts */
			pNH = pNeighbors;
			while(pNH != NULL)
			{
				nj = pNH->nPos;
				pCi = pSeqMtf->ppSamplePath[2]->pMatElement+nj;
				if(nj == ni)
				{
					*pCi = nNeighborCount+1;
				}
				else
				{
					*pCi = (*pCi)+1;
				}
				pNH = pNH->pNext;
			}
		}

		/* get next pos */
		pAi++;
		pBase++;
	}

	/* clear neighboring position list */
	SITEPOSDESTROYLIST(&pNeighbors);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Quality_EstimateTruePosProb: Get the true positive probability.        */
/* ----------------------------------------------------------------------- */ 
double Quality_EstimateTruePosProb(int nNeighborCount, int niter, int nIteration)
{
	/* define */
	double dProb;

	double d1,d2,d3;
	d1 = 0.04;
	d2 = 0.4;
	d3 = 0.95;
	
	/* get pdf */
	if(nNeighborCount <= 1)
	{
		dProb = 0.999-(1.0-d1)*niter*2.0/nIteration;
		if(dProb < d1)
			dProb = d1;
	}
	else if(nNeighborCount == 2)
	{
		dProb = 0.999-(1.0-d2)*niter*2.0/nIteration;
		if(dProb < d2)
			dProb = d2;
	}
	else if(nNeighborCount == 3)
	{
		dProb = 0.999-(1.0-d3)*niter*2.0/nIteration;
		if(dProb < d3)
			dProb = d3;
	}
	else if(nNeighborCount >= 4)
		dProb = 0.999;

	/* return */
	return dProb;
}

/* ----------------------------------------------------------------------- */ 
/*  Quality_Neighbor_Likelihood: Get the likelihood for neighbors.         */
/* ----------------------------------------------------------------------- */ 
int Quality_Neighbor_Likelihood(struct SITEPOS *pNeighbors, unsigned char *pCount,
								unsigned char *pB, double *pLn0, double *pLn1,
								int niter, int nIteration)
{
	/* define */
	struct SITEPOS *pS;
	int ni,nj;

	/* init */
	*pLn0 = 0.0;
	*pLn1 = 0.0;

	/* calculate */
	pS = pNeighbors;
	while(pS != NULL)
	{
		ni = pS->nPos;
		nj = (int)(pCount[ni]);
		if(pB[ni] == 1)
		{
			*pLn0 = (*pLn0) + log(Quality_EstimateTruePosProb(nj, niter, nIteration));
			*pLn1 = (*pLn1) + log(Quality_EstimateTruePosProb(nj+1, niter, nIteration));
		}
		else
		{
			*pLn0 = (*pLn0) + log(1.0-Quality_EstimateTruePosProb(nj, niter, nIteration));
			*pLn1 = (*pLn1) + log(1.0-Quality_EstimateTruePosProb(nj+1, niter, nIteration));
		}

		/* get next */
		pS = pS->pNext;
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Quality_BaseLogLikelihood_Background: get loglikelihood for background */
/* ----------------------------------------------------------------------- */ 
double Quality_BaseLogLikelihood_Background(unsigned char *pSite, int nMaxMotifLen, 
											struct DOUBLEMATRIX *pLogBG, int nBGOrder, 
											struct DOUBLEMATRIX *pLogBG0, int nPreNum, 
											char chStrand, int *pMask)
{
	/* define */
	double dL;
	int ni,nNnum,nk;
	int nWordId;
	int nBaseTypeNum;
	int nScale;
	int nTruePreNum;

	/* init */
	*pMask = 0;

	dL = 0.0;
	nNnum = 0;
	nBaseTypeNum = pLogBG->nWidth;
	nScale = (int)(pow((double)nBaseTypeNum, (double)(nBGOrder-1)));
	nWordId = 0;

	/* get likelihood */
	if(nBGOrder == 0)
	{
		for(ni=0; ni<nMaxMotifLen; ni++)
		{
			nk = (int)pSite[ni];
			if(nk >= nBaseTypeNum)
			{
				nNnum++;
				break;
			}

			dL += DMGETAT(pLogBG, 0, nk);
		}
	}
	else if(nBGOrder > 0)
	{
		for(ni=nPreNum; ni>0; ni--)
		{
			nk = (int)(*(pSite-ni));
			if(nk >= pLogBG->nWidth)
			{
				nNnum++;
				break;
			}
			nWordId = nWordId*nBaseTypeNum+nk;
		}
		
		if(nNnum > 0)
		{
			nWordId = 0;
			nTruePreNum = 0;
			nNnum = 0;
		}
		else
		{
			nTruePreNum = nPreNum;
		}

		for(ni=0; ni<nBGOrder-nTruePreNum; ni++)
		{
			if(ni == nMaxMotifLen)
				break;
			nk = (int)pSite[ni];
			if(nk >= nBaseTypeNum)
			{
				nNnum++;
				break;
			}
			dL += DMGETAT(pLogBG0, 0, nk);
			nWordId = nWordId*nBaseTypeNum+nk;
		}

		if((nNnum == 0) && (ni<nMaxMotifLen))
		{
			for(; ni<nMaxMotifLen; ni++)
			{
				nk = (int)pSite[ni];
				if(nk >= nBaseTypeNum)
				{
					nNnum++;
					break;
				}

				/* get mc likelihood */
				dL += DMGETAT(pLogBG, nWordId, nk);
				nWordId -= nScale*(int)(pSite[ni-nBGOrder]);
				nWordId = nWordId*nBaseTypeNum+nk;
			}
		}
	}
	else
	{
		/* report wrong */
		printf("Error: BaseLikelihood_Background, nBGOrder should be >=0!\n");
		exit(EXIT_FAILURE);
	}


	if(nNnum > 0)
	{
		*pMask = 1;
		dL = 0.0;
	}

	/* return */
	return dL;
}

/* ----------------------------------------------------------------------- */ 
/*  Quality_BaseLogLikelihood_Motif: get loglikelihood for motif.          */
/* ----------------------------------------------------------------------- */ 
double Quality_BaseLogLikelihood_Motif(unsigned char *pSite, int nMaxMotifLen, 
									   int nMotifLen, struct DOUBLEMATRIX *pLogPWMp, 
									   struct DOUBLEMATRIX *pLogBG, int nBGOrder, 
									   struct DOUBLEMATRIX *pLogBG0, 
									   char chStrand, int *pMask)
{
	/* define */
	double dL;
	int ni,nNnum,nk;
	int nWordId;
	int nBaseTypeNum;
	int nScale;
	int nPreNum;
	unsigned char *pB;

	/* check */
	*pMask = 0;
	if(pLogPWMp == NULL)
	{
		printf("Error: Quality_BaseLogLikelihood_Motif, no PWM matrix exist for likelihood calculation!\n");
		exit(EXIT_FAILURE);
	}
	if(pLogPWMp->nHeight != nMotifLen)
	{
		printf("Error: Quality_BaseLogLikelihood_Motif, motif length not match!\n");
		exit(EXIT_FAILURE);
	}

	/* get likelihood */
	dL = 0.0;
	nNnum = 0;
	if(chStrand == '+')
	{
		for(ni=0; ni<nMotifLen; ni++)
		{
			if(pSite[ni] >= pLogPWMp->nWidth)
			{
				nNnum++;
				break;
			}

			nk = (int)(pSite[ni]);
			dL += DMGETAT(pLogPWMp, ni, nk);
		}
	}
	else if(chStrand == '-')
	{
		for(ni=0; ni<nMotifLen; ni++)
		{
			if(pSite[ni] >= pLogPWMp->nWidth)
			{
				nNnum++;
				break;
			}

			nk = 3-(int)(pSite[ni]);
			dL += DMGETAT(pLogPWMp, (nMotifLen-1-ni), nk);
		}
	}
	else
	{
		printf("Error: BaseLikelihood_Motif, strand information wrong!\n");
		exit(EXIT_FAILURE);
	}

	if(nNnum > 0)
	{
		*pMask = 1;
		dL = 0.0;
	}
	else if(nMotifLen == nMaxMotifLen)
	{
		*pMask = 0;
	}
	else
	{
		/* init */
		nBaseTypeNum = pLogBG->nWidth;
		nScale = (int)(pow((double)nBaseTypeNum, (double)(nBGOrder-1)));
		nWordId = 0;
		
		pB = pSite+nMotifLen;
		/* get nPrenum */
		if(nBGOrder <= nMotifLen)
			nPreNum = nBGOrder;
		else
			nPreNum = nMotifLen;

		for(ni=nPreNum; ni>0; ni--)
		{
			nk = (int)(*(pB-ni));
			nWordId = nWordId*nBaseTypeNum+nk;
		}
		
		for(ni=0; ni<nBGOrder-nPreNum; ni++)
		{
			if((ni+nMotifLen) == nMaxMotifLen)
				break;
			nk = (int)pB[ni];
			if(nk >= nBaseTypeNum)
			{
				nNnum++;
				break;
			}
			dL += DMGETAT(pLogBG0, 0, nk);
			nWordId = nWordId*nBaseTypeNum+nk;
		}

		if((nNnum == 0) && ((ni+nMotifLen)<nMaxMotifLen))
		{
			for(; (ni+nMotifLen)<nMaxMotifLen; ni++)
			{
				nk = (int)pB[ni];
				if(nk >= nBaseTypeNum)
				{
					nNnum++;
					break;
				}

				/* get mc likelihood */
				dL += DMGETAT(pLogBG, nWordId, nk);
				nWordId -= nScale*(int)(pB[ni-nBGOrder]);
				nWordId = nWordId*nBaseTypeNum+nk;
			}
		}

		if(nNnum > 0)
		{
			*pMask = 1;
			dL = 0.0;
		}
	}

	/* return */
	return dL;
}

/* ----------------------------------------------------------------------- */ 
/*  Quality_CallMotifSites: Call motif sites.                              */
/* ----------------------------------------------------------------------- */ 
int Quality_CallMotifSites(struct SEQMOTIF *pSeqMtf, int vMotifLen[])
{
	/* define */
	int ni,nLen;
	unsigned char *pBi,*pAi;
	struct MOTIFSITE *pMtfSite,*pPrevSite;
	int nMtfId,nMtfSt;

	/* init */
	if( pSeqMtf->pMotifSite != NULL)
	{
		MTFSDESTROYLIST(&(pSeqMtf->pMotifSite));
		pSeqMtf->pMotifSite = NULL;
	}
	pMtfSite = NULL;
	pPrevSite = NULL;

	nLen = pSeqMtf->ppSamplePath[1]->nWidth;
	pBi = pSeqMtf->ppSamplePath[1]->pMatElement;
	pAi = pSeqMtf->ppSamplePath[0]->pMatElement;
	for(ni=0; ni<nLen; ni++)
	{
		if(*pBi > 0)
		{
			pMtfSite = NULL;
			pMtfSite = MTFSCREATE();
			
			nMtfId = (*pAi)/2-1;
			nMtfSt = (*pAi)%2;

			pMtfSite->nMotifType = nMtfId;
			pMtfSite->nMotifLen = vMotifLen[nMtfId];
			pMtfSite->nSeqId = 0;
			pMtfSite->nStartPos = ni;
			pMtfSite->nStrand = nMtfSt;
			
			if(pPrevSite == NULL)
			{
				pSeqMtf->pMotifSite = pMtfSite;
				pPrevSite = pMtfSite;
			}
			else
			{
				pPrevSite->pNext = pMtfSite;
				pPrevSite = pMtfSite;
			}
		}

		/* get next */
		pBi++;
		pAi++;
	}

	/* return */
	return PROC_SUCCESS;
}

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
						   struct DOUBLEMATRIX *pSeedMotif[])
{
	/* define */
	/* all sequences */
	struct SEQMOTIF **pSeqMtf;
	/* background */
	struct DOUBLEMATRIX *pBG;
	/* background0 */
	struct DOUBLEMATRIX *pBG0;
	/* sequence number */
	int nCount;
	/* motif matrix */
	struct MOTIFMATRIX *vMotif[MAX_MOTIF_NUM];
	/* motif length */
	int vMotifLen[MAX_MOTIF_NUM];

	/* HMM status num */
	int nHMMStatusNum;
	int nHMMHalfNum;
	int nRecordNum;
	int nIsRecording;

	/* background base transition */
	int nScale;

	/* prior frequency */
	double vFreqPrior[2][MAX_HMMSTATUS_NUM];
	/* frequency count */
	double vFreqCount[2][MAX_HMMSTATUS_NUM];
	double dMotifPTotal;
	
	/* number of base types */
	int nBaseTypeNum;
	/* others */
	int ni,nj,nk,nLen,niter;
	double dTotal;
	double *pElement,*pElement2;
	double *pScoreEle;
	FILE *fpSeg;

	/* initial check */
	nBaseTypeNum = 4;
	nScale = (int)(pow((double)nBaseTypeNum, (double)(nBGOrder-1)));
	nHMMHalfNum = (2*nMotifNum+1);
	nHMMStatusNum = nBGTypeNum*nHMMHalfNum;

	for(ni=0; ni<nBGTypeNum; ni++)
	{
		if(pFreqPrior[ni] == NULL)
		{
			printf("Error: MotifSampler_HMM, prior frequency of observing motif must be specified!\n");
			exit(EXIT_FAILURE);
		}
		if( (pFreqPrior[ni]->nHeight != (nMotifNum+1)) || (pFreqPrior[ni]->nWidth != 2) )
		{
			printf("Error: MotifSampler_HMM, prior frequency of observing motif and motif number not match!\n");
			exit(EXIT_FAILURE);
		}
	}
	if(pModuleTransP == NULL)
	{
		printf("Error: MotifSampler_HMM, module transition probability matrix must be specified!\n");
		exit(EXIT_FAILURE);
	}
	pElement = pModuleTransP->pMatElement;
	for(ni=0; ni<pModuleTransP->nHeight; ni++)
	{
		for(nj=0; nj<pModuleTransP->nWidth; nj++)
		{
			*pElement = log(*pElement);
			pElement++;
		}
	}
	for(ni=0; ni<nMotifNum; ni++)
	{
		if(pMatrixPrior[ni] == NULL)
		{
			printf("Error: MotifSampler_HMM, prior motif matrix must be specified!\n");
			exit(EXIT_FAILURE);
		}
		if(pMatrixPrior[ni]->nWidth != nBaseTypeNum)
		{
			printf("Error: MotifSampler_Collapsed, priors not properly specified!\n");
			exit(EXIT_FAILURE);
		}
	}
	if(nSeedNum>0)
	{
		for(ni=0; ni<nSeedNum; ni++)
		{
			if(pSeedMotif[ni] == NULL)
			{
				printf("Error: MotifSampler_HMM, prior seed motif matrix must be specified!\n");
				exit(EXIT_FAILURE);
			}
			
			if(pSeedMotif[ni]->nWidth != nBaseTypeNum)
			{
				printf("Error: MotifSampler_Collapsed, priors not properly specified!\n");
				exit(EXIT_FAILURE);
			}
		}
	}

	/* init prior probability */
	for(ni=0; ni<nBGTypeNum; ni++)
	{
		nk = 0;
		dTotal = 0.0;
		vFreqCount[ni][nk] = DMGETAT(pFreqPrior[ni], 0, 0);
		dTotal += vFreqCount[ni][nk];
		nk++;
		for(nj=1; nj<pFreqPrior[ni]->nHeight; nj++)
		{
			vFreqCount[ni][nk] = DMGETAT(pFreqPrior[ni], nj, 0);
			dTotal += vFreqCount[ni][nk];
			nk++;
			vFreqCount[ni][nk] = DMGETAT(pFreqPrior[ni], nj, 1);
			dTotal += vFreqCount[ni][nk];
			nk++;
		}

		for(nj=0; nj<nHMMHalfNum; nj++)
		{
			vFreqPrior[ni][nj] = vFreqCount[ni][nj]/dTotal;
		}
	}
	
	/* init motif */
	for(ni=0; ni<MAX_MOTIF_NUM; ni++)
	{
		vMotif[ni] = NULL;
	}
	/* prepare motif matrix */
	for(ni=0; ni<nMotifNum; ni++)
	{
		/* create */
		vMotif[ni] = MTFMCREATE();
		if(vMotif[ni] == NULL)
		{
			printf("Error: MotifSampler_HMM, cannot create motif matrix!\n");
			exit(EXIT_FAILURE);
		}

		/* init */
		vMotif[ni]->nMotifType = ni;
		vMotif[ni]->pCOUNTp = DMCLONE(pMatrixPrior[ni]);
		vMotifLen[ni] = vMotif[ni]->pCOUNTp->nHeight;
	}

	/* load sequences */
	pSeqMtf = NULL;
	pSeqMtf = SEQMTFLOADPRIMARYSEQFROMFASTA(&nCount, strInFile);

	/* calculate background */
	pBG = NULL;
	pBG = SEQMTFESTIMATENUCLEICBGMC(pSeqMtf, nCount, 0, nBGOrder);
	/* DMSAVE(pBG, "bgestimate.txt");*/
	pBG0 = NULL;
	pBG0 = SEQMTFESTIMATENUCLEICBGMC(pSeqMtf, nCount, 0, 0);

	/* get init state */
	for(ni=0; ni<nCount; ni++)
	{
		HMM_RandInit(pSeqMtf[ni], nBGOrder, pBG, nMotifNum, vMotif, vMotifLen,
					nSeedNum, pSeedMotif, vFreqPrior[1]);
	}
	/* update motif matrix */
	for(ni=0; ni<nMotifNum; ni++)
	{
		MTFMSAMPLEPWM(vMotif[ni]);
		/* set pPWMn = log(pPWMp) */
		pElement = vMotif[ni]->pPWMp->pMatElement;
		pElement2 = vMotif[ni]->pPWMn->pMatElement;
		for(nj=0; nj<vMotif[ni]->pPWMp->nHeight; nj++)
		{
			for(nk=0; nk<vMotif[ni]->pPWMp->nWidth; nk++)
			{
				*pElement2 = log(*pElement);
				pElement++;
				pElement2++;
			}
		}
	}

	/* take log of prior */
	for(ni=0; ni<nBGTypeNum; ni++)
	{
		for(nj=0; nj<nHMMHalfNum; nj++)
		{
			vFreqPrior[ni][nj] = log(vFreqPrior[ni][nj]);
		}
	}


	/* hmm calculation iterations */
	nRecordNum = MODULE_SAMPLE_NUM;
	if((int)(nIteration/2) < nRecordNum)
		nRecordNum = (int)(nIteration/2);
	for(niter=0; niter<nIteration; niter++)
	{
		if(niter%50 == 0)
		{
			printf("iter=%d...\n", niter);
		}

		if((niter+nRecordNum) < nIteration)
			nIsRecording = 0;
		else
			nIsRecording = 1;

		/* prepare counting spaces */
		for(ni=0; ni<nMotifNum; ni++)
		{
			DestroyDoubleMatrix(vMotif[ni]->pCOUNTp);
			vMotif[ni]->pCOUNTp = NULL;
			vMotif[ni]->pCOUNTp = DMCLONE(pMatrixPrior[ni]);
		}

		/* initialize the count vectors */
		for(ni=0; ni<nBGTypeNum; ni++)
		{
			nk = 0;
			vFreqCount[ni][nk] = DMGETAT(pFreqPrior[ni], 0, 0);
			nk++;
			for(nj=1; nj<pFreqPrior[ni]->nHeight; nj++)
			{
				vFreqCount[ni][nk] = DMGETAT(pFreqPrior[ni], nj, 0);
				nk++;
				vFreqCount[ni][nk] = DMGETAT(pFreqPrior[ni], nj, 1);
				nk++;
			}
		}

		/* sequence-wise hmm processing */
		for(ni=0; ni<nCount; ni++)
		{
			/* create HMM memory */			
			nLen = pSeqMtf[ni]->ppSamplePath[0]->nWidth;
			SEQMTFCREATEHMM(pSeqMtf[ni], nHMMStatusNum, nLen);
			
			/* forward summation */
			HMM_ForwardSummation(pSeqMtf[ni], nHMMStatusNum,
						 pModuleTransP, vFreqPrior,
						 nBGOrder, pBG,
						 nMotifNum, vMotifLen, 
						 vMotif);
			
			/* backward sampling */
			HMM_BackwardSampling(pSeqMtf[ni], nHMMStatusNum,
						 pModuleTransP, vFreqCount,
						 nMotifNum, vMotifLen, 
						 vMotif, nIsRecording);
			
			/* destroy HMM memory */
			SEQMTFDESTROYHMM(pSeqMtf[ni]);
		}

		/* update motif matrix */
		for(ni=0; ni<nMotifNum; ni++)
		{
			MTFMSAMPLEPWM(vMotif[ni]);
			/* set pPWMn = log(pPWMp) */
			pElement = vMotif[ni]->pPWMp->pMatElement;
			pElement2 = vMotif[ni]->pPWMn->pMatElement;
			for(nj=0; nj<vMotif[ni]->pPWMp->nHeight; nj++)
			{
				for(nk=0; nk<vMotif[ni]->pPWMp->nWidth; nk++)
				{
					*pElement2 = log(*pElement);
					pElement++;
					pElement2++;
				}
			}
		}

		/* update prior frequencies */
		for(ni=0; ni<nBGTypeNum; ni++)
		{
			dTotal = 0.0;
			for(nj=0; nj<nHMMHalfNum; nj++)
			{
				dTotal += vFreqCount[ni][nj];
			}

			for(nj=0; nj<nHMMHalfNum; nj++)
			{
				vFreqPrior[ni][nj] = vFreqCount[ni][nj]/dTotal;
			}

			/* cut at the lower bound and take log */
			if(ni == 0)
			{
				dMotifPTotal = 1.0-vFreqPrior[ni][0];
				if( dMotifPTotal > BG_MOTIF_DENSITY_UPPERBOUND )
				{
					vFreqPrior[ni][0] = 1.0;
					for(nj=1; nj<nHMMHalfNum; nj++)
					{
						vFreqPrior[ni][nj] *= BG_MOTIF_DENSITY_UPPERBOUND/dMotifPTotal;
						vFreqPrior[ni][0] -= vFreqPrior[ni][nj];
					}
				}
			}
			else
			{
				dMotifPTotal = 1.0-vFreqPrior[ni][0];
				if( dMotifPTotal < MODULE_MOTIF_DENSITY_LOWERBOUND )
				{
					vFreqPrior[ni][0] = 1.0;
					for(nj=1; nj<nHMMHalfNum; nj++)
					{
						vFreqPrior[ni][nj] *= MODULE_MOTIF_DENSITY_LOWERBOUND/dMotifPTotal;
						vFreqPrior[ni][0] -= vFreqPrior[ni][nj];
					}
				}
			}

			/* take log */
			for(nj=0; nj<nHMMHalfNum; nj++)
			{
				vFreqPrior[ni][nj] = log(vFreqPrior[ni][nj]);
			}
		}

		/* TODO: local mode shifting */

		/* TODO: change motif length */
	}

	/* get module probability for every base */
	for(ni=0; ni<nCount; ni++)
	{
		pScoreEle = pSeqMtf[ni]->ppScore[0]->pMatElement;
		nLen = pSeqMtf[ni]->ppScore[0]->nWidth;
		for(nj=0; nj<nLen; nj++)
		{
			*pScoreEle = (*pScoreEle)/(double)nRecordNum;
			pScoreEle++;
		}
	}

	/* call motif based on the posterior probability */
	for(ni=0; ni<nBGTypeNum; ni++)
	{
		for(nj=0; nj<nHMMHalfNum; nj++)
		{
			vFreqPrior[ni][nj] = exp(vFreqPrior[ni][nj]);
		}
	}
	for(ni=0; ni<nCount; ni++)
	{
		HMM_CallMotifSites(pSeqMtf[ni], nBGOrder, pBG, nMotifNum, vMotif, vMotifLen,
					vFreqPrior);
	}

	/* report result */
	WriteMotifToFile(strOutFile, nMotifNum, vMotif, pBG0, nCount, pSeqMtf);
	fpSeg = fopen("segout.txt", "wt");
	for(ni=0; ni<nCount; ni++)
	{
		SEQMTFWRITESAMPPATHTOFASTA(fpSeg, pSeqMtf[ni]);
	}
	fclose(fpSeg);

	/* destroy */
	DestroyDoubleMatrix(pBG);
	DestroyDoubleMatrix(pBG0);
	for(ni=0; ni<nCount; ni++)
	{
		SEQMTFDESTROY(pSeqMtf[ni]);
		pSeqMtf[ni] = NULL;
	}
	free(pSeqMtf);
	for(ni=0; ni<nMotifNum; ni++)
	{
		MTFMDESTROY(vMotif[ni]);
		vMotif[ni] = NULL;
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HMM_RandInit: Randomly initialize the motifs.                          */
/* ----------------------------------------------------------------------- */ 
int HMM_RandInit(struct SEQMOTIF *pSeqMtf, 
				 int nBGOrder, struct DOUBLEMATRIX *pBG, 
				 int nMotifNum, struct MOTIFMATRIX *vMotif[], int vMotifLen[],
				 int nSeedNum, struct DOUBLEMATRIX *pSeedMotif[], 
				 double vFreqPrior[])
{
	/* define */
	double vFreqPost[MAX_HMMSTATUS_NUM];
	int ni,nj,nk,nx,nLen;
	unsigned char *pBi;
	double dLike;
	int nOK;
	int nPreNum;
	double dTotal;
	double dRand;

	/* check */
	if(nSeedNum > nMotifNum)
	{
		printf("Error: motif seeds should be <= motif number!\n");
		exit(EXIT_FAILURE);
	}

	/* init */
	nLen = pSeqMtf->ppSeq[0]->nWidth;
	pBi = pSeqMtf->ppSeq[0]->pMatElement;
	for(ni=0; ni<nLen; ni++)
	{
		if(ni < nBGOrder)
			nPreNum = ni;
		else
			nPreNum = nBGOrder;
		nOK = 0;
		/* sample with motifs with information */
		for(nj=0; nj<nSeedNum; nj++)
		{
			if(ni+pSeedMotif[nj]->nHeight-1 >= nLen)
			{
				vFreqPost[0] = 1.0;
				vFreqPost[1] = 0.0;
				vFreqPost[2] = 0.0;
			}
			else
			{
				dTotal = 0.0;
				nk = 2*nj+1;
				dLike = BaseLikelihood_Background(pBi, pSeedMotif[nj]->nHeight, pBG, nBGOrder, nPreNum, '+');
				vFreqPost[0] = (1.0-vFreqPrior[nk]-vFreqPrior[nk+1])*dLike;
				dTotal += vFreqPost[0];
				dLike = BaseLikelihood_Motif(pBi, pSeedMotif[nj]->nHeight, pSeedMotif[nj], '+');
				vFreqPost[1] = vFreqPrior[nk]*dLike;
				dTotal += vFreqPost[1];
				dLike = BaseLikelihood_Motif(pBi, pSeedMotif[nj]->nHeight, pSeedMotif[nj], '-');
				vFreqPost[2] = vFreqPrior[nk+1]*dLike;
				dTotal += vFreqPost[2];

				for(nx=0; nx<3; nx++)
				{
					vFreqPost[nx] = vFreqPost[nx]/dTotal;
				}
			}

			dRand = rand_u();
			dTotal = 0.0;
			for(nx=0; nx<3; nx++)
			{
				dTotal += vFreqPost[nx];
				if(dRand <= dTotal)
					break;
			}
			if(nx>=3) nx--;
			if(nx>0)
			{
				nOK = 1;
				if(nx == 1)
					AddNucleicMatrixCount(vMotif[nj]->pCOUNTp, '+', pBi, vMotifLen[nj]);
				else
					AddNucleicMatrixCount(vMotif[nj]->pCOUNTp, '-', pBi, vMotifLen[nj]);
				break;
			}
		}

		if(nOK == 1)
		{
			pBi++;
			continue;
		}

		/* sample with motifs without information */
		for(; nj<nMotifNum; nj++)
		{
			if(ni+vMotifLen[nj]-1 >= nLen)
			{
				vFreqPost[0] = 1.0;
				vFreqPost[1] = 0.0;
				vFreqPost[2] = 0.0;
			}
			else
			{
				nk = 2*nj+1;
				vFreqPost[0] = (1.0-vFreqPrior[nk]-vFreqPrior[nk+1]);
				vFreqPost[1] = vFreqPrior[nk];
				vFreqPost[2] = vFreqPrior[nk+1];
			}

			dRand = rand_u();
			dTotal = 0.0;
			for(nx=0; nx<3; nx++)
			{
				dTotal += vFreqPost[nx];
				if(dRand <= dTotal)
					break;
			}
			if(nx>=3) nx--;
			if(nx>0)
			{
				nOK = 1;
				if(nx == 1)
					AddNucleicMatrixCount(vMotif[nj]->pCOUNTp, '+', pBi, vMotifLen[nj]);
				else
					AddNucleicMatrixCount(vMotif[nj]->pCOUNTp, '-', pBi, vMotifLen[nj]);
				break;
			}
		}

		/* get next */
		pBi++;
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HMM_CallMotifSites: Call motif sites.                                  */
/* ----------------------------------------------------------------------- */ 
int HMM_CallMotifSites(struct SEQMOTIF *pSeqMtf,
					   int nBGOrder, struct DOUBLEMATRIX *pBG,
					   int nMotifNum, struct MOTIFMATRIX *vMotif[], int vMotifLen[],
					   double vFreqPrior[2][MAX_HMMSTATUS_NUM])
{
	/* define */
	double vFreqPost[MAX_HMMSTATUS_NUM];
	int ni,nj,nk,nx,nLen;
	unsigned char *pBi;
	double *pSi;
	double dLike;
	int nPreNum;
	double dTotal;
	double dP0,dP1,dMP;
	struct MOTIFSITE *pMtfSite,*pPrevSite;

	/* init */
	if( pSeqMtf->pMotifSite != NULL)
	{
		MTFSDESTROYLIST(&(pSeqMtf->pMotifSite));
		pSeqMtf->pMotifSite = NULL;
	}
	pMtfSite = NULL;
	pPrevSite = NULL;


	nLen = pSeqMtf->ppSeq[0]->nWidth;
	pBi = pSeqMtf->ppSeq[0]->pMatElement;
	pSi = pSeqMtf->ppScore[0]->pMatElement;
	for(ni=0; ni<nLen; ni++)
	{
		if(ni < nBGOrder)
			nPreNum = ni;
		else
			nPreNum = nBGOrder;

		/* process motif by motif */
		for(nj=0; nj<nMotifNum; nj++)
		{
			if(ni+vMotifLen[nj]-1 >= nLen)
			{
				vFreqPost[0] = 1.0;
				vFreqPost[1] = 0.0;
				vFreqPost[2] = 0.0;
			}
			else
			{
				nk = 2*nj+1;
				dMP = *pSi;
				dP0 = (1.0-dMP)*vFreqPrior[0][nk]+dMP*vFreqPrior[1][nk];
				dP1 = (1.0-dMP)*vFreqPrior[0][nk+1]+dMP*vFreqPrior[1][nk+1];

				dTotal = 0.0;
				dLike = BaseLikelihood_Background(pBi, vMotifLen[nj], pBG, nBGOrder, nPreNum, '+');
				vFreqPost[0] = (1.0-dP0-dP1)*dLike;
				dTotal += vFreqPost[0];
				dLike = BaseLikelihood_Motif(pBi, vMotifLen[nj], vMotif[nj]->pCOUNTp, '+');
				vFreqPost[1] = dP0*dLike;
				dTotal += vFreqPost[1];
				dLike = BaseLikelihood_Motif(pBi, vMotifLen[nj], vMotif[nj]->pCOUNTp, '-');
				vFreqPost[2] = dP1*dLike;
				dTotal += vFreqPost[2];

				for(nx=0; nx<3; nx++)
				{
					vFreqPost[nx] = vFreqPost[nx]/dTotal;
				}
			}

			if(vFreqPost[0] < MOTIF_CALL_CUTOFF)
			{
				pMtfSite = NULL;
				pMtfSite = MTFSCREATE();
				pMtfSite->nMotifType = nj;
				pMtfSite->nMotifLen = vMotifLen[nj];
				pMtfSite->nSeqId = 0;
				pMtfSite->nStartPos = ni;
				if(vFreqPost[1] >= vFreqPost[2])
					pMtfSite->nStrand = 0;
				else
					pMtfSite->nStrand = 1;

				if(pPrevSite == NULL)
				{
					pSeqMtf->pMotifSite = pMtfSite;
					pPrevSite = pMtfSite;
				}
				else
				{
					pPrevSite->pNext = pMtfSite;
					pPrevSite = pMtfSite;
				}
			}
		}

		/* get next */
		pBi++;
		pSi++;
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HMM_ForwardSummation: Forward integration procedure of HMM.            */
/* ----------------------------------------------------------------------- */ 
int HMM_ForwardSummation(struct SEQMOTIF *pSeqMtf, int nHMMStatusNum,
						 struct DOUBLEMATRIX *pModuleTransP,
						 double vFreqPrior[2][MAX_HMMSTATUS_NUM],
						 int nBGOrder, struct DOUBLEMATRIX *pBG,
						 int nMotifNum, int vMotifLen[], 
						 struct MOTIFMATRIX *vMotif[])
{
	/* DEFINE */
	int nHMMHalfNum;
	int nj,nk,nx;
	struct DOUBLEMATRIX *pHMM;
	unsigned char *pBi;
	double dModuleInitProb;
	int nBaseTypeNum;
	int nWordId,nScale;
	double dBaseTransP;
	int nLen;

	/* current module status & past status */
	int nCurrentMS,nCurrentCellId,nCurrentMotifId,nCurrentMtfStrand;
	int nLastMS;
	double vHMMCol[MAX_HMMSTATUS_NUM];
	double vHMMInit[MAX_HMMSTATUS_NUM];
	double dHMMMax;
	double dMotifLike;
	double dTotal;
	int nMotifStart;

	/* CHECK */
	if(nHMMStatusNum != pSeqMtf->pHMM->nHeight)
	{
		printf("Error: HMM status number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* INIT */
	dModuleInitProb = 0.5;
	nBaseTypeNum = pBG->nWidth;
	nScale = (int)(pow((double)nBaseTypeNum, (double)(nBGOrder-1)));

	pHMM = pSeqMtf->pHMM;
	nHMMHalfNum = (int)(nHMMStatusNum/2);
		
	/* SET INIT STATUS */
	dTotal = 1.0+MARGIN_MOTIF_PROB*(nHMMStatusNum-2);
	for(nk=0; nk<nHMMStatusNum; nk++)
	{
		if(nk==0)
			vHMMInit[nk] = log((1.0-dModuleInitProb)/dTotal);
		else if(nk==nHMMHalfNum)
			vHMMInit[nk] = log(dModuleInitProb/dTotal);
		else
			vHMMInit[nk] = log(MARGIN_MOTIF_PROB/dTotal);
	}

	/* INTEGRATION FOR nj>=0 */
	nLen = pHMM->nWidth;
	pBi = pSeqMtf->ppSeq[0]->pMatElement;
	for(nj=0; nj<nLen; nj++)
	{
		if(*pBi >= nBaseTypeNum)
		{
			printf("Error: base type N not allowed!\n");
			exit(EXIT_FAILURE);
		}

		/* get background base transition prob */
		if(nj < nBGOrder)
			dBaseTransP = log(1.0/((double)nBaseTypeNum));
		else
		{
			nWordId = 0;
			for(nk=nBGOrder; nk>=1; nk--)
			{
				nWordId = nWordId*nBaseTypeNum+(int)(*(pBi-nk));
			}
			dBaseTransP = DMGETAT(pBG, nWordId, (int)(*pBi));
			dBaseTransP = log(dBaseTransP);
		}

		/* process status by status */
		for(nk=0; nk<nHMMStatusNum; nk++)
		{
			/* get current module status */
			if(nk < nHMMHalfNum)
			{
				nCurrentMS = 0;
				nCurrentCellId = nk;
			}
			else
			{
				nCurrentMS = 1;
				nCurrentCellId = nk-nHMMHalfNum;
			}

			/* if background end */
			if(nCurrentCellId == 0)
			{
				/* if background out of range */
				if(nj == 0)
				{
					for(nx=0; nx<nHMMStatusNum; nx++)
					{
						/* get last module status */
						if(nx < nHMMHalfNum)
						{
							nLastMS = 0;
						}
						else
						{
							nLastMS = 1;
						}

						/* module transition likelihood */
						vHMMCol[nx] = vHMMInit[nx] + DMGETAT(pModuleTransP, nLastMS, nCurrentMS)
										+ vFreqPrior[nCurrentMS][0] + dBaseTransP;
						
						/* get max */
						if(nx == 0)
							dHMMMax = vHMMCol[nx];
						else
						{
							if(vHMMCol[nx] > dHMMMax)
								dHMMMax = vHMMCol[nx];
						}
					}
				}
				/* normal background */
				else
				{
					/* process previous status one by one */
					for(nx=0; nx<nHMMStatusNum; nx++)
					{
						/* get last module status */
						if(nx < nHMMHalfNum)
						{
							nLastMS = 0;
						}
						else
						{
							nLastMS = 1;
						}

						/* module transition likelihood */
						vHMMCol[nx] = DMGETAT(pHMM, nx, nj-1) + DMGETAT(pModuleTransP, nLastMS, nCurrentMS)
										+ vFreqPrior[nCurrentMS][0] + dBaseTransP;
						
						/* get max */
						if(nx == 0)
							dHMMMax = vHMMCol[nx];
						else
						{
							if(vHMMCol[nx] > dHMMMax)
								dHMMMax = vHMMCol[nx];
						}
					}
				}

				/* integration */
				dTotal = 0.0;
				for(nx=0; nx<nHMMStatusNum; nx++)
				{
					dTotal += exp(vHMMCol[nx]-dHMMMax);
				}
				dTotal = log(dTotal)+dHMMMax;

				/* set new values */
				DMSETAT(pHMM, nk, nj, dTotal);
			}

			/* if motif end */
			else
			{
				/* get motif id */
				nCurrentMotifId = (int)((nCurrentCellId-1)/2);
				nCurrentMtfStrand = (nCurrentCellId-1)%2;
				
				nMotifStart = (nj-vMotifLen[nCurrentMotifId]+1);

				/* if motif out of range */
				if(nMotifStart < 0)
				{
					/* process previous status one by one */
					for(nx=0; nx<nHMMStatusNum; nx++)
					{
						/* get last module status */
						if(nx < nHMMHalfNum)
						{
							nLastMS = 0;
						}
						else
						{
							nLastMS = 1;
						}

						/* module transition likelihood */
						vHMMCol[nx] = vHMMInit[nx] + DMGETAT(pModuleTransP, nLastMS, nCurrentMS)
							+ DMGETAT(pModuleTransP, nCurrentMS, nCurrentMS)*nj
							+ vFreqPrior[nCurrentMS][nCurrentCellId]
							+ vFreqPrior[nCurrentMS][0]*nj
							+ log(MARGIN_MOTIF_PROB)*(nj+1);
						
						/* get max */
						if(nx == 0)
							dHMMMax = vHMMCol[nx];
						else
						{
							if(vHMMCol[nx] > dHMMMax)
								dHMMMax = vHMMCol[nx];
						}
					}
				}
				/* if motif starts from base 0 */
				else if(nMotifStart == 0)
				{
					if(nCurrentMtfStrand == 0)
					{
						dMotifLike = BaseLogLikelihood_MotifPWM((pBi-(vMotifLen[nCurrentMotifId]-1)),
							vMotifLen[nCurrentMotifId], 
							vMotif[nCurrentMotifId]->pPWMn, '+');
					}
					else
					{
						dMotifLike = BaseLogLikelihood_MotifPWM((pBi-(vMotifLen[nCurrentMotifId]-1)), 
							vMotifLen[nCurrentMotifId], 
							vMotif[nCurrentMotifId]->pPWMn, '-');
					}


					/* process previous status one by one */
					for(nx=0; nx<nHMMStatusNum; nx++)
					{
						/* get last module status */
						if(nx < nHMMHalfNum)
						{
							nLastMS = 0;
						}
						else
						{
							nLastMS = 1;
						}

						/* module transition likelihood */
						vHMMCol[nx] = vHMMInit[nx] + DMGETAT(pModuleTransP, nLastMS, nCurrentMS)
							+ DMGETAT(pModuleTransP, nCurrentMS, nCurrentMS)*nj
							+ vFreqPrior[nCurrentMS][nCurrentCellId]
							+ vFreqPrior[nCurrentMS][0]*nj
							+ dMotifLike;
						
						/* get max */
						if(nx == 0)
							dHMMMax = vHMMCol[nx];
						else
						{
							if(vHMMCol[nx] > dHMMMax)
								dHMMMax = vHMMCol[nx];
						}
					}

				}
				/* normal motifs */
				else
				{
					if(nCurrentMtfStrand == 0)
					{
						dMotifLike = BaseLogLikelihood_MotifPWM((pBi-(vMotifLen[nCurrentMotifId]-1)),
							vMotifLen[nCurrentMotifId], 
							vMotif[nCurrentMotifId]->pPWMn, '+');
					}
					else
					{
						dMotifLike = BaseLogLikelihood_MotifPWM((pBi-(vMotifLen[nCurrentMotifId]-1)), 
							vMotifLen[nCurrentMotifId], 
							vMotif[nCurrentMotifId]->pPWMn, '-');
					}

					for(nx=0; nx<nHMMStatusNum; nx++)
					{
						/* get last module status */
						if(nx < nHMMHalfNum)
						{
							nLastMS = 0;
						}
						else
						{
							nLastMS = 1;
						}

						/* module transition likelihood */
						vHMMCol[nx] = DMGETAT(pHMM, nx, nMotifStart-1) 
									+ DMGETAT(pModuleTransP, nLastMS, nCurrentMS)
									+ DMGETAT(pModuleTransP, nCurrentMS, nCurrentMS)*(vMotifLen[nCurrentMotifId]-1)
									+ vFreqPrior[nCurrentMS][nCurrentCellId] 
									+ vFreqPrior[nCurrentMS][0]*(vMotifLen[nCurrentMotifId]-1)
									+ dMotifLike;
						
						/* get max */
						if(nx == 0)
							dHMMMax = vHMMCol[nx];
						else
						{
							if(vHMMCol[nx] > dHMMMax)
								dHMMMax = vHMMCol[nx];
						}
					}
				}

				/* integration */
				dTotal = 0.0;
				for(nx=0; nx<nHMMStatusNum; nx++)
				{
					dTotal += exp(vHMMCol[nx]-dHMMMax);
				}
				dTotal = log(dTotal)+dHMMMax;

				/* set new values */
				DMSETAT(pHMM, nk, nj, dTotal);
			}
		}

		/* get next base */
		pBi++;
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HMM_BackwardSampling: Backward sampling procedure of HMM.              */
/* ----------------------------------------------------------------------- */ 
int HMM_BackwardSampling(struct SEQMOTIF *pSeqMtf, int nHMMStatusNum,
						 struct DOUBLEMATRIX *pModuleTransP, 
						 double vFreqCount[2][MAX_HMMSTATUS_NUM],
						 int nMotifNum, int vMotifLen[], 
						 struct MOTIFMATRIX *vMotif[],
						 int nRecordModuleOn)
{
	/* define */
	int nj,nx,ny,nLen;
	struct DOUBLEMATRIX *pHMM;
	double dTotal;
	double vHMMCol[MAX_HMMSTATUS_NUM];
	double dHMMMax;
	int nHMMHalfNum;
	double dRand;
	int nCurrentMS,nCurrentCellId,nCurrentMotifId,nCurrentMtfStrand;
	int nNextMS;
	double *pScoreEle;
	unsigned char *pPathEle;

	/* init */
	nHMMHalfNum = (int)(nHMMStatusNum/2);
	pHMM = pSeqMtf->pHMM;
	nLen = pHMM->nWidth;

	/* sample the last base */
	nj=(nLen-1);
	for(nx=0; nx<nHMMStatusNum; nx++)
	{
		vHMMCol[nx] = DMGETAT(pHMM, nx, nj);
		if(nx == 0)
			dHMMMax = vHMMCol[nx];
		else
		{
			if(vHMMCol[nx] > dHMMMax)
				dHMMMax = vHMMCol[nx];
		}
	}
	dTotal = 0.0;
	for(nx=0; nx<nHMMStatusNum; nx++)
	{
		vHMMCol[nx] -= dHMMMax;
		dTotal += exp(vHMMCol[nx]);
	}
	dTotal = log(dTotal);
	for(nx=0; nx<nHMMStatusNum; nx++)
	{
		vHMMCol[nx] = exp(vHMMCol[nx]-dTotal);
	}
	dRand = rand_u();
	dTotal = 0.0;
	for(nx=0; nx<nHMMStatusNum; nx++)
	{
		dTotal += vHMMCol[nx];
		if(dRand <= dTotal)
		{
			break;
		}
	}
	if(nx == nHMMStatusNum)
		nx--;

	/* get current module status */
	if(nx < nHMMHalfNum)
	{
		nCurrentMS = 0;
		nCurrentCellId = nx;
	}
	else
	{
		nCurrentMS = 1;
		nCurrentCellId = nx-nHMMHalfNum;
	}

	/* update frequency counts */
	vFreqCount[nCurrentMS][nCurrentCellId] += 1.0; 
	/* update sample path */
	if(nCurrentCellId == 0)
	{
		BMSETAT(pSeqMtf->ppSamplePath[0], 0, nj, (unsigned char)nx);
		nj--;
	}
	else
	{
		nCurrentMotifId = (int)((nCurrentCellId-1)/2);
		nCurrentMtfStrand = (nCurrentCellId-1)%2;
		if(nCurrentMS == 0)
		{
			for(ny=0; ny<(vMotifLen[nCurrentMotifId]-1); ny++)
			{
				BMSETAT(pSeqMtf->ppSamplePath[0], 0, nj, 0);
				nj--;
			}
		}
		else
		{
			for(ny=0; ny<(vMotifLen[nCurrentMotifId]-1); ny++)
			{
				BMSETAT(pSeqMtf->ppSamplePath[0], 0, nj, (unsigned char)nHMMHalfNum);
				nj--;
			}
		}

		BMSETAT(pSeqMtf->ppSamplePath[0], 0, nj, (unsigned char)nx);
		nj--;

		/* modify motif matrix */
		if(nCurrentMtfStrand == 0)
		{
			AddNucleicMatrixCount(vMotif[nCurrentMotifId]->pCOUNTp, 
				'+', (pSeqMtf->ppSeq[0]->pMatElement+nj+1), 
				vMotifLen[nCurrentMotifId]);
		}
		else
		{
			AddNucleicMatrixCount(vMotif[nCurrentMotifId]->pCOUNTp, 
				'-', (pSeqMtf->ppSeq[0]->pMatElement+nj+1), 
				vMotifLen[nCurrentMotifId]);
		}
	}

	/* prepare for next step */
	nNextMS = nCurrentMS;

	/* sample other bases */
	while(nj>=0)
	{
		for(nx=0; nx<nHMMStatusNum; nx++)
		{
			if(nx < nHMMHalfNum)
			{
				nCurrentMS = 0;
			}
			else
			{
				nCurrentMS = 1;
			}
			vHMMCol[nx] = DMGETAT(pHMM, nx, nj)+DMGETAT(pModuleTransP, nCurrentMS, nNextMS);
			if(nx == 0)
				dHMMMax = vHMMCol[nx];
			else
			{
				if(vHMMCol[nx] > dHMMMax)
					dHMMMax = vHMMCol[nx];
			}
		}
		dTotal = 0.0;
		for(nx=0; nx<nHMMStatusNum; nx++)
		{
			vHMMCol[nx] -= dHMMMax;
			dTotal += exp(vHMMCol[nx]);
		}
		dTotal = log(dTotal);
		for(nx=0; nx<nHMMStatusNum; nx++)
		{
			vHMMCol[nx] = exp(vHMMCol[nx]-dTotal);
		}
		dRand = rand_u();
		dTotal = 0.0;
		for(nx=0; nx<nHMMStatusNum; nx++)
		{
			dTotal += vHMMCol[nx];
			if(dRand <= dTotal)
			{
				break;
			}
		}
		if(nx == nHMMStatusNum)
			nx--;

		/* get current module status */
		if(nx < nHMMHalfNum)
		{
			nCurrentMS = 0;
			nCurrentCellId = nx;
		}
		else
		{
			nCurrentMS = 1;
			nCurrentCellId = nx-nHMMHalfNum;
		}

		/* update frequency counts */
		vFreqCount[nCurrentMS][nCurrentCellId] += 1.0; 
		/* update sample path */
		if(nCurrentCellId == 0)
		{
			BMSETAT(pSeqMtf->ppSamplePath[0], 0, nj, (unsigned char)nx);
			nj--;
		}
		else
		{
			nCurrentMotifId = (int)((nCurrentCellId-1)/2);
			nCurrentMtfStrand = (nCurrentCellId-1)%2;
			if(nCurrentMS == 0)
			{
				for(ny=0; ny<(vMotifLen[nCurrentMotifId]-1); ny++)
				{
					BMSETAT(pSeqMtf->ppSamplePath[0], 0, nj, 0);
					nj--;
				}
			}
			else
			{
				for(ny=0; ny<(vMotifLen[nCurrentMotifId]-1); ny++)
				{
					BMSETAT(pSeqMtf->ppSamplePath[0], 0, nj, (unsigned char)nHMMHalfNum);
					nj--;
				}
			}

			BMSETAT(pSeqMtf->ppSamplePath[0], 0, nj, (unsigned char)nx);
			nj--;

			/* modify motif matrix */
			if(nCurrentMtfStrand == 0)
			{
				AddNucleicMatrixCount(vMotif[nCurrentMotifId]->pCOUNTp, 
					'+', (pSeqMtf->ppSeq[0]->pMatElement+nj+1), 
					vMotifLen[nCurrentMotifId]);
			}
			else
			{
				AddNucleicMatrixCount(vMotif[nCurrentMotifId]->pCOUNTp, 
					'-', (pSeqMtf->ppSeq[0]->pMatElement+nj+1), 
					vMotifLen[nCurrentMotifId]);
			}
		}

		/* prepare for next step */
		nNextMS = nCurrentMS;
	}

	/* record module status */
	if(nRecordModuleOn == 1)
	{
		pScoreEle = pSeqMtf->ppScore[0]->pMatElement;
		pPathEle = pSeqMtf->ppSamplePath[0]->pMatElement;
		for(nj=0; nj<nLen; nj++)
		{
			if(*pPathEle >= nHMMHalfNum)
				*pScoreEle = (*pScoreEle)+1.0;
			pScoreEle++;
			pPathEle++;
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  WriteMotifToFile: Report motif discovery results.                      */
/* ----------------------------------------------------------------------- */ 
int WriteMotifToFile(char strOutFile[], 
					 int nMotifNum, struct MOTIFMATRIX *vMotif[],
					 struct DOUBLEMATRIX *pBG0,
					 int nCount, struct SEQMOTIF *vSeqMtf[])
{
	/* define */
	FILE *fpOut;
	int ni,nj,nk,nSeqId;
	double *pElement;
	char strConsensus[LINE_LENGTH];
	struct MOTIFSITE *pSite;
	unsigned char *pBi;
	double dScore;
	int nMax;
	double dMax;
	int nBaseTypeNum;

	/* report results */
	nBaseTypeNum = pBG0->nWidth;

	fpOut = NULL;
	fpOut = fopen(strOutFile, "wt");
	if(fpOut == NULL)
	{
		printf("Error: cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	/* write motif one by one */
	for(ni=0; ni<nMotifNum; ni++)
	{
		/* get motif score */
		dScore = MTFMSCORE(vMotif[ni], pBG0);

		/* write head */
		fprintf(fpOut, "****** Motif%d ******\n", ni);
		fprintf(fpOut, "Motif Score: %f\n", dScore);
		fprintf(fpOut, "Motif Matrix: \n");
		pElement = vMotif[ni]->pCOUNTp->pMatElement;
		for(nj=0; nj<vMotif[ni]->pCOUNTp->nHeight; nj++)
		{
			dMax = 0.0;
			nMax = nBaseTypeNum;
			for(nk=0; nk<vMotif[ni]->pCOUNTp->nWidth; nk++)
			{
				fprintf(fpOut, "% 9.7e ", *pElement);
				if(*pElement > dMax)
				{
					dMax = *pElement;
					nMax = nk;
				}
				pElement++;
			}
			fprintf(fpOut, "\n");
			switch(nMax)
			{
				case 0: strConsensus[nj] = 'A';
					break;
				case 1: strConsensus[nj] = 'C';
					break;
				case 2: strConsensus[nj] = 'G';
					break;
				case 3: strConsensus[nj] = 'T';
					break;
				default: strConsensus[nj] = 'N';
			}
		}
		strConsensus[nj] = '\0';

		fprintf(fpOut, "\n");
		fprintf(fpOut, "Consensus:\n");
		fprintf(fpOut, "%s\n", strConsensus);

		fprintf(fpOut, "\nMotif Sites:\n");
		for(nk=0; nk<nCount; nk++)
		{
			pSite = vSeqMtf[nk]->pMotifSite;
			while(pSite != NULL)
			{
				/* if match, then write */
				if(pSite->nMotifType != ni)
				{
					pSite = pSite->pNext;
					continue;
				}

				/* head information */
				nSeqId = pSite->nSeqId;
				fprintf(fpOut, ">%d|%d|(%d:%d)|", vSeqMtf[nk]->nId, vSeqMtf[nk]->pSeqId[nSeqId],
						pSite->nStartPos, (pSite->nStartPos+pSite->nMotifLen-1));
				if(pSite->nStrand == 0)
				{
					fprintf(fpOut, "+\n");
				}
				else if(pSite->nStrand == 1)
				{
					fprintf(fpOut, "-\n");
				}
				else
				{
					fprintf(fpOut, "?\n");
				}

				/* write seq */
				if(pSite->nStrand == 0)
				{
					pBi = vSeqMtf[nk]->ppSeq[nSeqId]->pMatElement + pSite->nStartPos;
					for(nj=0; nj<pSite->nMotifLen; nj++)
					{
						switch((int)(*pBi))
						{
							case 0: fprintf(fpOut, "A");
								break;
							case 1: fprintf(fpOut, "C");
								break;
							case 2: fprintf(fpOut, "G");
								break;
							case 3: fprintf(fpOut, "T");
								break;
							default: fprintf(fpOut, "N");
						}
						pBi++;
					}
				}
				else if(pSite->nStrand == 1)
				{
					pBi = vSeqMtf[nk]->ppSeq[nSeqId]->pMatElement + pSite->nStartPos + pSite->nMotifLen - 1;
					for(nj=0; nj<pSite->nMotifLen; nj++)
					{
						switch((int)(*pBi))
						{
							case 0: fprintf(fpOut, "T");
								break;
							case 1: fprintf(fpOut, "G");
								break;
							case 2: fprintf(fpOut, "C");
								break;
							case 3: fprintf(fpOut, "A");
								break;
							default: fprintf(fpOut, "N");
						}
						pBi--;
					}
				}
				fprintf(fpOut, "\n");

				/* get next */
				pSite = pSite->pNext;
			}
		}

		fprintf(fpOut, "\n\n");
	}

	/* close file */
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifSampler_Collapsed: Motif sampler based on collapsed gibbs sampler */
/*  Find one motif per call                                                */
/* ----------------------------------------------------------------------- */ 
int MotifSampler_Collapsed(char strInFile[], char strOutFile[], int nIteration, 
						   int nBGOrder, struct DOUBLEMATRIX *pFreqPrior,
						   struct DOUBLEMATRIX *pMatrixPrior)
{
	/* define */
	/* all sequences */
	struct SEQMOTIF **pSeqMtf;
	/* background */
	struct DOUBLEMATRIX *pBG;
	/* background0 */
	struct DOUBLEMATRIX *pBG0;
	/* sequence number */
	int nCount;
	/* motif matrix */
	struct MOTIFMATRIX *pMotif;
	/* motif length */
	int nMotifLen;
	/* prior frequency */
	double dFreqPrior[3];
	/* frequency count */
	double dFreqCount[3];
	/* posterior frequency */
	double dFreqPost[3];
	/* number of bases could be used for background calculation */
	int nPreNum;
	/* number of base types */
	int nBaseTypeNum;
	/* output file */
	FILE *fpOut;
	/* motif score */
	double dScore;
	/* others */
	int ni,nj,nk,nx,ny,nMask,nLen,niter;
	unsigned char *pAi,*pBi,*pCi;
	double dRand,dTotal,dMax;
	int nMax;
	double *pElement;

	/* initial check */
	if((pFreqPrior == NULL) || (pMatrixPrior == NULL))
	{
		printf("Error: MotifSampler_Collapsed, prior frequency of observing motif and prior motif matrix (both in the form of pseudocounts) must be specified!\n");
		exit(EXIT_FAILURE);
	}
	if(pFreqPrior->nWidth != 3)
	{
		printf("Error: MotifSampler_Collapsed, priors not properly specified!\n");
		exit(EXIT_FAILURE);
	}

	/* set base type number */
	nBaseTypeNum = 4;

	/* load sequences */
	pSeqMtf = NULL;
	pSeqMtf = SEQMTFLOADPRIMARYSEQFROMFASTA(&nCount, strInFile);

	/* calculate background */
	pBG = NULL;
	pBG = SEQMTFESTIMATENUCLEICBGMC(pSeqMtf, nCount, 0, nBGOrder);
	/* DMSAVE(pBG, "bgestimate.txt");*/
	pBG0 = NULL;
	pBG0 = SEQMTFESTIMATENUCLEICBGMC(pSeqMtf, nCount, 0, 0);
	
	/* prepare motif matrix */
	pMotif = NULL;
	pMotif = MTFMCREATE();
	if(pMotif == NULL)
	{
		printf("Error: MotifSampler_Collapsed, cannot create motif matrix!\n");
		exit(EXIT_FAILURE);
	}
	dTotal = 0.0;
	for(ni=0; ni<3; ni++)
	{
		dFreqPrior[ni] = DMGETAT(pFreqPrior, 0, ni);
		dFreqCount[ni] = dFreqPrior[ni];
		dTotal += dFreqPrior[ni];
	}
	if(dTotal <= 0.0)
	{
		printf("Error: MotifSampler_Collapsed, prior frequency (in the form of pseudocounts) of observing motif not properly specified!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=1; ni<3; ni++)
	{
		dFreqPrior[ni] += dFreqPrior[ni-1];
	}
	for(ni=0; ni<3; ni++)
	{
		dFreqPrior[ni] /= dTotal;
	}
	
	/* set initial points */
	nMotifLen = pMatrixPrior->nHeight;
	pMotif->pCOUNTp = DMCLONE(pMatrixPrior);
	for(ni=0; ni<nCount; ni++)
	{
		nLen = pSeqMtf[ni]->ppSamplePath[0]->nWidth;
		pAi = pSeqMtf[ni]->ppSamplePath[0]->pMatElement;
		pBi = pSeqMtf[ni]->ppSeq[0]->pMatElement;
		for(nj=0; nj<(nLen-nMotifLen+1); nj++)
		{
			if((int)(*pBi) >= nBaseTypeNum)
			{
				*pAi = 0;
			}
			else
			{
				dRand = rand_u();
				dTotal = 0.0;
				for(nk=0; nk<3; nk++)
				{
					if(dRand <= dFreqPrior[nk])
						break;
				}
				if(nk>=3)
				{
					printf("Error: floating point error!\n");
					exit(EXIT_FAILURE);
				}
				*pAi = (unsigned char)nk;
				dFreqCount[nk] += 1.0;

				/* add initial TFBS matrix pseudocount */
				switch(nk)
				{
					case 1: AddNucleicMatrixCount(pMotif->pCOUNTp, '+', (pSeqMtf[ni]->ppSeq[0]->pMatElement+nj), nMotifLen);
						break;
					case 2:	AddNucleicMatrixCount(pMotif->pCOUNTp, '-', (pSeqMtf[ni]->ppSeq[0]->pMatElement+nj), nMotifLen);
						break;
					default: 
						break;
				}
			}
			pAi++;
			pBi++;
		}
		for(; nj<nLen; nj++)
		{
			if((int)(*pBi) >= nBaseTypeNum)
			{
				*pAi = 0;
			}
			else
			{
				dFreqCount[0] += 1.0;
			}
			pAi++;
			pBi++;
		}
	}

	/* collapsed sampler */
	for(niter=0; niter<nIteration; niter++)
	{
		if(niter%50 == 0)
		{
			printf("iter=%d...\n", niter);
		}
		for(ni=0; ni<nCount; ni++)
		{
			/* sample */
			nLen = pSeqMtf[ni]->ppSamplePath[0]->nWidth;
			pAi = pSeqMtf[ni]->ppSamplePath[0]->pMatElement;
			pBi = pSeqMtf[ni]->ppSeq[0]->pMatElement;
			for(nj=0; nj<(nLen-nMotifLen+1); nj++)
			{
				/* clear current status */
				if((int)(*pBi) >= nBaseTypeNum)
				{
					/* get next pos */
					pAi++;
					pBi++;
					continue;
				}

				nk = (int)(*pAi);
				dFreqCount[nk] -= 1.0;
				switch(nk)
				{
					case 1: SubtractNucleicMatrixCount(pMotif->pCOUNTp, '+', (pSeqMtf[ni]->ppSeq[0]->pMatElement+nj), nMotifLen);
						break;
					case 2:	SubtractNucleicMatrixCount(pMotif->pCOUNTp, '-', (pSeqMtf[ni]->ppSeq[0]->pMatElement+nj), nMotifLen);
						break;
					default: 
						break;
				}

				*pAi = 0;
				nMask = 0;
				ny = nj-nMotifLen+1;
				if(ny<0)
					ny = 0;
				pCi = pSeqMtf[ni]->ppSamplePath[0]->pMatElement+ny;
				for(; ny<nj+nMotifLen-1; ny++)
				{
					if(*pCi > 0)
					{
						nMask = 1;
						break;
					}
					pCi++;
				}

				/* get posterior */
				if(nMask == 1)
				{
					dFreqPost[0] = 1.0;
					dFreqPost[1] = 0.0;
					dFreqPost[2] = 0.0;
				}
				else
				{
				
					if(nj < nBGOrder)
					{
						nPreNum = nj;
					}
					else
					{
						nPreNum = nBGOrder;
					}

					
					dFreqPost[0] = dFreqCount[0]*BaseLikelihood_Background((pSeqMtf[ni]->ppSeq[0]->pMatElement+nj), nMotifLen, pBG, nBGOrder, nPreNum, '+');
					dFreqPost[1] = dFreqCount[1]*BaseLikelihood_Motif((pSeqMtf[ni]->ppSeq[0]->pMatElement+nj), nMotifLen, pMotif->pCOUNTp, '+');
					dFreqPost[2] = dFreqCount[2]*BaseLikelihood_Motif((pSeqMtf[ni]->ppSeq[0]->pMatElement+nj), nMotifLen, pMotif->pCOUNTp, '-');
					
					/*dFreqPost[0] = dFreqPrior[0]*BaseLikelihood_Background((pSeqMtf[ni]->ppSeq[0]->pMatElement+nj), nMotifLen, pBG, nBGOrder, nPreNum, '+');
					dFreqPost[1] = (dFreqPrior[1]-dFreqPrior[0])*BaseLikelihood_Motif((pSeqMtf[ni]->ppSeq[0]->pMatElement+nj), nMotifLen, pMotif->pCOUNTp, '+');
					dFreqPost[2] = (dFreqPrior[2]-dFreqPrior[1])*BaseLikelihood_Motif((pSeqMtf[ni]->ppSeq[0]->pMatElement+nj), nMotifLen, pMotif->pCOUNTp, '-');
					*/

					dTotal = 0.0;
					for(nx=0; nx<3; nx++)
					{
						dTotal += dFreqPost[nx];
					}
					if(dFreqPost[0] <= 0.0)
						dTotal = 0.0;

					if(dTotal > 0.0)
					{
						for(nx=0; nx<3; nx++)
						{
							dFreqPost[nx] /= dTotal;
						}
					}
					else
					{
						for(nx=0; nx<3; nx++)
						{
							dFreqPost[nx] = dFreqPrior[nx];
						}
						for(nx=2; nx>0; nx--)
						{
							dFreqPost[nx] -= dFreqPost[nx-1];
						}
					}
				}

				/* sample status */
				dRand = rand_u();
				dTotal = 0.0;
				for(nx=0; nx<3; nx++)
				{
					dTotal += dFreqPost[nx];
					if(dRand <= dTotal)
						break;
				}
				nk = nx;
				if(nk>=3)
				{
					printf("Error: floating point error!\n");
					exit(EXIT_FAILURE);
				}

				/* update parameters */
				*pAi = (unsigned char)nk;
				dFreqCount[nk] += 1.0;

				/* add initial TFBS matrix pseudocount */
				switch(nk)
				{
					case 1: AddNucleicMatrixCount(pMotif->pCOUNTp, '+', (pSeqMtf[ni]->ppSeq[0]->pMatElement+nj), nMotifLen);
						break;
					case 2:	AddNucleicMatrixCount(pMotif->pCOUNTp, '-', (pSeqMtf[ni]->ppSeq[0]->pMatElement+nj), nMotifLen);
						break;
					default: 
						break;
				}

				/* get next pos */
				pAi++;
				pBi++;
			}
		}
	}
	
	/* get motif score */
	dScore = MTFMSCORE(pMotif, pBG0);

	/* report results */
	fpOut = NULL;
	fpOut = fopen(strOutFile, "wt");
	if(fpOut == NULL)
	{
		printf("Error: cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	/* write head */
	fprintf(fpOut, "Motif Score: %f\n", dScore);
	pElement = pMotif->pCOUNTp->pMatElement;
	for(nj=0; nj<pMotif->pCOUNTp->nHeight; nj++)
	{
		dMax = 0.0;
		nMax = nBaseTypeNum;
		for(nk=0; nk<pMotif->pCOUNTp->nWidth; nk++)
		{
			fprintf(fpOut, "% 9.7e ", *pElement);
			if(*pElement > dMax)
			{
				dMax = *pElement;
				nMax = nk;
			}

			pElement++;
		}
		switch(nMax)
		{
			case 0: fprintf(fpOut, "A");
				break;
			case 1: fprintf(fpOut, "C");
				break;
			case 2: fprintf(fpOut, "G");
				break;
			case 3: fprintf(fpOut, "T");
				break;
			default: fprintf(fpOut, "N");
		}

		fprintf(fpOut, "\n");
	}

	fprintf(fpOut, "\n");

	for(ni=0; ni<nCount; ni++)
	{
		/* write seq */
		nLen = pSeqMtf[ni]->ppSamplePath[0]->nWidth;
		pAi = pSeqMtf[ni]->ppSamplePath[0]->pMatElement;
		for(nj=0; nj<(nLen-nMotifLen+1); nj++)
		{
			if(*pAi > 0)
			{
				fprintf(fpOut, ">%d|%d|(%d:%d)|", pSeqMtf[ni]->nId, pSeqMtf[ni]->pSeqId[0],
					nj, (nj+nMotifLen-1));
				if(*pAi == 1)
				{
					fprintf(fpOut, "+\n");
				}
				else if(*pAi == 2)
				{
					fprintf(fpOut, "-\n");
				}
				else
				{
					fprintf(fpOut, "?\n");
				}
				if(*pAi == 1)
				{
					pBi = pSeqMtf[ni]->ppSeq[0]->pMatElement+nj;
					for(nk=0; nk<nMotifLen; nk++)
					{
						switch((int)(*pBi))
						{
							case 0: fprintf(fpOut, "A");
								break;
							case 1: fprintf(fpOut, "C");
								break;
							case 2: fprintf(fpOut, "G");
								break;
							case 3: fprintf(fpOut, "T");
								break;
							default: fprintf(fpOut, "N");
						}
						pBi++;
					}
				}
				else if(*pAi == 2)
				{
					pBi = pSeqMtf[ni]->ppSeq[0]->pMatElement+nj+nMotifLen-1;
					for(nk=0; nk<nMotifLen; nk++)
					{
						switch((int)(*pBi))
						{
							case 0: fprintf(fpOut, "T");
								break;
							case 1: fprintf(fpOut, "G");
								break;
							case 2: fprintf(fpOut, "C");
								break;
							case 3: fprintf(fpOut, "A");
								break;
							default: fprintf(fpOut, "N");
						}
						pBi--;
					}
				}
				fprintf(fpOut, "\n");
			}

			/* get next */
			pAi++;
		}
	}
	fclose(fpOut);


	/* destroy */
	DestroyDoubleMatrix(pBG);
	DestroyDoubleMatrix(pBG0);
	for(ni=0; ni<nCount; ni++)
	{
		SEQMTFDESTROY(pSeqMtf[ni]);
		pSeqMtf[ni] = NULL;
	}
	free(pSeqMtf);
	MTFMDESTROY(pMotif);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  AddNucleicMatrixCount: Add count to motif matrix, for nucleic acid.    */
/* ----------------------------------------------------------------------- */ 
int AddNucleicMatrixCount(struct DOUBLEMATRIX *pCountMat, char chStrand, 
						  unsigned char pSite[], int nMotifLen)
{
	/* define */
	int ni,nj,nk;
	double dTemp;

	/* check */
	if(pCountMat == NULL)
	{
		printf("Error: AddNucleicMatrixCount, no count matrix exist for counting!\n");
		exit(EXIT_FAILURE);
	}
	if(pCountMat->nHeight != nMotifLen)
	{
		printf("Error: AddNucleicMatrixCount, motif lengths not match!\n");
		exit(EXIT_FAILURE);
	}

	/* add count from a '+' strand site */
	if(chStrand == '+')
	{
		for(ni=0; ni<nMotifLen; ni++)
		{
			if(pSite[ni] >= 4)
				continue;
			nj = (int)pSite[ni];
			dTemp = DMGETAT(pCountMat, ni, nj) + 1.0;
			DMSETAT(pCountMat, ni, nj, dTemp);
		}
	}
	/* add count from a '-' strand site */
	else if(chStrand == '-')
	{
		for(ni=0; ni<nMotifLen; ni++)
		{
			if(pSite[ni] >= 4)
				continue;
			nj = 3-(int)pSite[ni];
			nk = nMotifLen-1-ni;
			dTemp = DMGETAT(pCountMat, nk, nj) + 1.0;
			DMSETAT(pCountMat, nk, nj, dTemp);
		}
	}
	/* do nothing */
	else
	{
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SubtractNucleicMatrixCount: Add count to motif matrix, nucleic acid.   */
/* ----------------------------------------------------------------------- */ 
int SubtractNucleicMatrixCount(struct DOUBLEMATRIX *pCountMat, char chStrand, 
						  unsigned char pSite[], int nMotifLen)
{
	/* define */
	int ni,nj,nk;
	double dTemp;

	/* check */
	if(pCountMat == NULL)
	{
		printf("Error: SubtractNucleicMatrixCount, no count matrix exist for counting!\n");
		exit(EXIT_FAILURE);
	}
	if(pCountMat->nHeight != nMotifLen)
	{
		printf("Error: SubtractNucleicMatrixCount, motif lengths not match!\n");
		exit(EXIT_FAILURE);
	}

	/* subtract count from a '+' strand site */
	if(chStrand == '+')
	{
		for(ni=0; ni<nMotifLen; ni++)
		{
			if(pSite[ni] >= 4)
				continue;
			nj = (int)pSite[ni];
			dTemp = DMGETAT(pCountMat, ni, nj) - 1.0;
			DMSETAT(pCountMat, ni, nj, dTemp);
		}
	}
	/* subtract count from a '-' strand site */
	else if(chStrand == '-')
	{
		for(ni=0; ni<nMotifLen; ni++)
		{
			if(pSite[ni] >= 4)
				continue;
			nj = 3-(int)pSite[ni];
			nk = nMotifLen-1-ni;
			dTemp = DMGETAT(pCountMat, nk, nj) - 1.0;
			DMSETAT(pCountMat, nk, nj, dTemp);
		}
	}
	/* do nothing */
	else
	{
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  BaseLikelihood_Background: get likelihood for background.              */
/* ----------------------------------------------------------------------- */ 
double BaseLikelihood_Background(unsigned char pSite[], int nMotifLen, struct DOUBLEMATRIX *pBG, int nBGOrder, int nPreNum, char chStrand)
{
	/* define */
	double dL;
	int ni,nNnum,nk;
	int nWordId;
	int nBaseTypeNum;
	int nScale;

	/* init */
	dL = 1.0;
	nNnum = 0;
	nBaseTypeNum = pBG->nWidth;
	nScale = (int)(pow((double)nBaseTypeNum, (double)(nBGOrder-1)));
	nWordId = 0;

	/* get likelihood */
	if(nBGOrder == 0)
	{
		for(ni=0; ni<nMotifLen; ni++)
		{
			nk = (int)pSite[ni];
			if(nk >= pBG->nWidth)
			{
				nNnum++;
				break;
			}

			dL = dL*DMGETAT(pBG, 0, nk);
		}
	}
	else if(nBGOrder > 0)
	{
		for(ni=nPreNum; ni>0; ni--)
		{
			nk = (int)(*(pSite-ni));
			if(nk >= pBG->nWidth)
			{
				nNnum++;
				break;
			}
			nWordId = nWordId*nBaseTypeNum+nk;
		}
		for(ni=0; ni<nBGOrder-nPreNum; ni++)
		{
			if(ni == nMotifLen)
				break;
			nk = (int)pSite[ni];
			if(nk >= pBG->nWidth)
			{
				nNnum++;
				break;
			}
			dL = dL*1.0/(double)nBaseTypeNum;
			nWordId = nWordId*nBaseTypeNum+nk;
		}

		if((nNnum == 0) && (ni<nMotifLen))
		{
			for(; ni<nMotifLen; ni++)
			{
				nk = (int)pSite[ni];
				if(nk >= pBG->nWidth)
				{
					nNnum++;
					break;
				}

				/* get mc likelihood */
				dL = dL*DMGETAT(pBG, nWordId, nk);
				nWordId -= nScale*(int)(pSite[ni-nBGOrder]);
				nWordId = nWordId*nBaseTypeNum+nk;
			}
		}
	}
	else
	{
		/* report wrong */
		printf("Error: BaseLikelihood_Background, nBGOrder should be >=0!\n");
		exit(EXIT_FAILURE);
	}


	if((double)nNnum > 0)
		dL = 0.0;

	/* return */
	return dL;
}

/* ----------------------------------------------------------------------- */ 
/*  BaseLikelihood_Motif: get collapsed likelihood for motif.              */
/* ----------------------------------------------------------------------- */ 
double BaseLikelihood_Motif(unsigned char pSite[], int nMotifLen, struct DOUBLEMATRIX *pCountMat, char chStrand)
{
	/* define */
	double dL;
	int ni,nj,nNnum,nk;
	double dTotal;

	/* check */
	if(pCountMat == NULL)
	{
		printf("Error: BaseLikelihood_Motif, no count matrix exist for likelihood calculation!\n");
		exit(EXIT_FAILURE);
	}
	if(pCountMat->nHeight != nMotifLen)
	{
		printf("Error: BaseLikelihood_Motif, motif length not match!\n");
		exit(EXIT_FAILURE);
	}

	/* get likelihood */
	dL = 1.0;
	nNnum = 0;
	if(chStrand == '+')
	{
		for(ni=0; ni<nMotifLen; ni++)
		{
			if(pSite[ni] >= pCountMat->nWidth)
			{
				nNnum++;
				break;
			}

			dTotal = 0.0;
			for(nj=0; nj<pCountMat->nWidth; nj++)
			{
				dTotal += DMGETAT(pCountMat, ni, nj);
			}
			if(dTotal <= 0.0)
			{
				nNnum++;
				break;
			}
			nk = (int)(pSite[ni]);
			dL = dL*DMGETAT(pCountMat, ni, nk)/dTotal;
		}
	}
	else if(chStrand == '-')
	{
		for(ni=0; ni<nMotifLen; ni++)
		{
			if(pSite[ni] >= pCountMat->nWidth)
			{
				nNnum++;
				break;
			}

			dTotal = 0.0;
			nk = 3-(int)(pSite[ni]);
			for(nj=0; nj<pCountMat->nWidth; nj++)
			{
				dTotal += DMGETAT(pCountMat, (nMotifLen-1-ni), nj);
			}
			if(dTotal <= 0.0)
			{
				nNnum++;
				break;
			}
			dL = dL*DMGETAT(pCountMat, (nMotifLen-1-ni), nk)/dTotal;
		}
	}
	else
	{
		printf("Error: BaseLikelihood_Motif, strand information wrong!\n");
		exit(EXIT_FAILURE);
	}

	if(nNnum > 0)
		dL = 0.0;

	/* return */
	return dL;
}


/* ----------------------------------------------------------------------- */ 
/*  BaseLogLikelihood_MotifPWM: get loglikelihood for motif.               */
/* ----------------------------------------------------------------------- */ 
double BaseLogLikelihood_MotifPWM(unsigned char pSite[], int nMotifLen, struct DOUBLEMATRIX *pLogPWMMat, char chStrand)
{
	/* define */
	double dL;
	int ni,nNnum,nk;

	/* check */
	if(pLogPWMMat == NULL)
	{
		printf("Error: BaseLogLikelihood_MotifPWM, no PWM matrix exist for likelihood calculation!\n");
		exit(EXIT_FAILURE);
	}
	if(pLogPWMMat->nHeight != nMotifLen)
	{
		printf("Error: BaseLogLikelihood_MotifPWM, motif length not match!\n");
		exit(EXIT_FAILURE);
	}

	/* get likelihood */
	dL = 0.0;
	nNnum = 0;
	if(chStrand == '+')
	{
		for(ni=0; ni<nMotifLen; ni++)
		{
			if(pSite[ni] >= pLogPWMMat->nWidth)
			{
				nNnum++;
				break;
			}

			nk = (int)(pSite[ni]);
			dL += DMGETAT(pLogPWMMat, ni, nk);
		}
	}
	else if(chStrand == '-')
	{
		for(ni=0; ni<nMotifLen; ni++)
		{
			if(pSite[ni] >= pLogPWMMat->nWidth)
			{
				nNnum++;
				break;
			}

			nk = 3-(int)(pSite[ni]);
			dL += DMGETAT(pLogPWMMat, (nMotifLen-1-ni), nk);
		}
	}
	else
	{
		printf("Error: BaseLogLikelihood_Motif, strand information wrong!\n");
		exit(EXIT_FAILURE);
	}

	if(nNnum > 0)
	{
		printf("Error: BaseLogLikelihood_Motif, base type N not allowed!\n");
		exit(EXIT_FAILURE);
	}

	/* return */
	return dL;
}

/* ----------------------------------------------------------------------- */ 
/*  BaseLikelihood_MotifPWM: get likelihood for motif.                     */
/* ----------------------------------------------------------------------- */ 
double BaseLikelihood_MotifPWM(unsigned char pSite[], int nMotifLen, struct DOUBLEMATRIX *pPWMMat, char chStrand)
{
	/* define */
	double dL;
	int ni,nNnum,nk;

	/* check */
	if(pPWMMat == NULL)
	{
		printf("Error: BaseLikelihood_MotifPWM, no PWM matrix exist for likelihood calculation!\n");
		exit(EXIT_FAILURE);
	}
	if(pPWMMat->nHeight != nMotifLen)
	{
		printf("Error: BaseLikelihood_MotifPWM, motif length not match!\n");
		exit(EXIT_FAILURE);
	}

	/* get likelihood */
	dL = 1.0;
	nNnum = 0;
	if(chStrand == '+')
	{
		for(ni=0; ni<nMotifLen; ni++)
		{
			if(pSite[ni] >= pPWMMat->nWidth)
			{
				nNnum++;
				break;
			}

			nk = (int)(pSite[ni]);
			dL *= DMGETAT(pPWMMat, ni, nk);
		}
	}
	else if(chStrand == '-')
	{
		for(ni=0; ni<nMotifLen; ni++)
		{
			if(pSite[ni] >= pPWMMat->nWidth)
			{
				nNnum++;
				break;
			}

			nk = 3-(int)(pSite[ni]);
			dL *= DMGETAT(pPWMMat, (nMotifLen-1-ni), nk);
		}
	}
	else
	{
		printf("Error: BaseLikelihood_Motif, strand information wrong!\n");
		exit(EXIT_FAILURE);
	}

	if(nNnum > 0)
	{
		printf("Error: BaseLikelihood_Motif, base type N not allowed!\n");
		exit(EXIT_FAILURE);
	}

	/* return */
	return dL;
}

/* ----------------------------------------------------------------------- */ 
/*                     MTFMCREATE: create PWM of motifs                    */
/* ----------------------------------------------------------------------- */ 
struct MOTIFMATRIX *MTFMCREATE()
{
	/* new motif matrix */
	struct MOTIFMATRIX *pM;

	/* create motif matrix */
	pM = NULL;
	pM = (struct MOTIFMATRIX*)calloc(1, sizeof(struct MOTIFMATRIX));
	if (pM == NULL)
	{
		printf("Error: cannot create motif matrix.\n");
		exit(EXIT_FAILURE);
	}
	
	/* init */
	pM->pCOUNTn = NULL;
	pM->pCOUNTp = NULL;
	pM->pPWMn = NULL;
	pM->pPWMp = NULL;
	pM->pNext = NULL;

	/* return */
	return pM;
}

/* ----------------------------------------------------------------------- */ 
/*                    MTFMDESTROY: destroy PWM of motifs                   */
/* ----------------------------------------------------------------------- */ 
int MTFMDESTROY(struct MOTIFMATRIX *pM)
{
	if(pM == NULL)
		return PROC_SUCCESS;

	/* destroy */
	pM->pNext = NULL;
	if(pM->pCOUNTn != NULL)
	{
		DestroyDoubleMatrix(pM->pCOUNTn);
		pM->pCOUNTn = NULL;
	}
	if(pM->pCOUNTp != NULL)
	{
		DestroyDoubleMatrix(pM->pCOUNTp);
		pM->pCOUNTp = NULL;
	}
	if(pM->pPWMn != NULL)
	{
		DestroyDoubleMatrix(pM->pPWMn);
		pM->pPWMn = NULL;
	}
	if(pM->pPWMp != NULL)
	{
		DestroyDoubleMatrix(pM->pPWMp);
		pM->pPWMp = NULL;
	}

	free(pM);
	pM = NULL;

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*                  MTFMDESTROYLIST: destroy PWM list of motifs            */
/* ----------------------------------------------------------------------- */ 
int MTFMDESTROYLIST(struct MOTIFMATRIX **pMList)
{
	/* new motif matrix */
	struct MOTIFMATRIX *pM;

	if(pMList==NULL)
		return PROC_SUCCESS;

	/* destroy */
	while(*pMList != NULL)
	{
		pM = *pMList;
		*pMList = pM->pNext;
		MTFMDESTROY(pM);
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*                  MTFMSCORE: get motif scores                            */
/* ----------------------------------------------------------------------- */ 
double MTFMSCORE(struct MOTIFMATRIX *pM, struct DOUBLEMATRIX *pBG)
{
	/* define */
	double dScore,dTemp1,dTemp2,dStemp,dTotal;
	int ni,nj;

	/* check */
	if((pM == NULL) || (pBG == NULL))
		return 0.0;

	/* init */
	dScore = 0.0;

	/* get score */
	for(ni=0; ni<pM->pCOUNTp->nHeight; ni++)
	{
		dStemp = 0.0;
		dTotal = 0.0;
		for(nj=0; nj<pM->pCOUNTp->nWidth; nj++)
		{
			dTotal += DMGETAT(pM->pCOUNTp, ni, nj);
		}
		for(nj=0; nj<pM->pCOUNTp->nWidth; nj++)
		{
			dTemp1 = DMGETAT(pM->pCOUNTp, ni, nj)/dTotal;
			dTemp2 = DMGETAT(pBG, 0, nj);
			dStemp += dTemp1*log(dTemp1/dTemp2);
		}
		dScore += log(dTotal)*dStemp;
	}

	dScore /= pM->pCOUNTp->nHeight;

	/* return */
	return dScore;
}

/* ----------------------------------------------------------------------- */ 
/*           MTFMREFRESH_QUALITY: set PWMp = log(COUNTp/sum(COUNTp)        */
/* ----------------------------------------------------------------------- */ 
int MTFMREFRESH_QUALITY(struct MOTIFMATRIX *pM)
{
	/* define */
	double dTotal;
	double *pEle1,*pEle2,*pEle3;
	int ni,nj;

	/* check */
	if(pM == NULL)
		return PROC_FAILURE;
	if(pM->pPWMp != NULL)
	{
		DestroyDoubleMatrix(pM->pPWMp);
		pM->pPWMp = NULL;
	}

	/* process */
	pM->pPWMp = CreateDoubleMatrix(pM->pCOUNTp->nHeight, pM->pCOUNTp->nWidth);
	pEle1 = pM->pCOUNTp->pMatElement;
	pEle2 = pM->pPWMp->pMatElement;
	for(ni=0; ni<pM->pCOUNTp->nHeight; ni++)
	{
		dTotal = 0.0;
		pEle3 = pEle1;
		for(nj=0; nj<pM->pCOUNTp->nWidth; nj++)
		{
			dTotal += (*pEle3);
			pEle3++;
		}

		for(nj=0; nj<pM->pCOUNTp->nWidth; nj++)
		{
			*pEle2 = log(*pEle1/dTotal);
			pEle1++;
			pEle2++;
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*                     MTFSCREATE: create site of motifs                   */
/* ----------------------------------------------------------------------- */ 
struct MOTIFSITE *MTFSCREATE()
{
	/* new motif site */
	struct MOTIFSITE *pS;

	/* create motif site */
	pS = NULL;
	pS = (struct MOTIFSITE*)calloc(1, sizeof(struct MOTIFSITE));
	if (pS == NULL)
	{
		printf("Error: cannot create motif site.\n");
		exit(EXIT_FAILURE);
	}

	pS->pNext = NULL;

	/* return */
	return pS;
}

/* ----------------------------------------------------------------------- */ 
/*                    MTFSDESTROY: destroy motif site                      */
/* ----------------------------------------------------------------------- */ 
int MTFSDESTROY(struct MOTIFSITE *pS)
{
	if(pS == NULL)
		return PROC_SUCCESS;

	/* destroy */
	pS->pNext = NULL;
	
	free(pS);
	pS = NULL;

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*                  MTFSDESTROYLIST: destroy list of motif sites           */
/* ----------------------------------------------------------------------- */ 
int MTFSDESTROYLIST(struct MOTIFSITE **pSList)
{
	/* new motif site */
	struct MOTIFSITE *pS;

	if(pSList==NULL)
		return PROC_SUCCESS;

	/* destroy */
	while(*pSList != NULL)
	{
		pS = *pSList;
		*pSList = pS->pNext;
		MTFSDESTROY(pS);
	}

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*                  MTFMSAMPLEPWM: sample PWM from counts                  */
/* ----------------------------------------------------------------------- */ 
int MTFMSAMPLEPWM(struct MOTIFMATRIX *pM)
{
	int ni,nj;
	double *pEle1,*pEle2;

	if(pM == NULL)
		return PROC_FAILURE;

	/* sample */
	DestroyDoubleMatrix(pM->pPWMp);
	pM->pPWMp = DMRANDPRODDIRICHLET(pM->pCOUNTp);
	DestroyDoubleMatrix(pM->pPWMn);
	pM->pPWMn = CreateDoubleMatrix(pM->pPWMp->nHeight, pM->pPWMp->nWidth);

	pEle1 = pM->pPWMp->pMatElement;
	pEle2 = pM->pPWMn->pMatElement + pM->pPWMn->nHeight*pM->pPWMn->nWidth - 1;
	for(ni=0; ni<pM->pPWMp->nHeight; ni++)
	{
		for(nj=0; nj<pM->pPWMp->nWidth; nj++)
		{
			*pEle2 = *pEle1;
			pEle1++;
			pEle2--;
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*                  SEQMTFCREATE: create seq/motif complex                 */
/* ----------------------------------------------------------------------- */ 
struct SEQMOTIF *SEQMTFCREATE(int nId, int nSeqNum, int nScoreNum, int nPathNum)
{
	/* new complex */
	struct SEQMOTIF *pSeqMtf;

	/* init */
	if((nSeqNum<0) || (nScoreNum<0) || (nPathNum<0))
	{
		printf("Error: seqmtfcreate, parameters must be no less than 0!\n");
		exit(EXIT_FAILURE);
	}

	/* create motif matrix */
	pSeqMtf = NULL;
	pSeqMtf = (struct SEQMOTIF *)calloc(1, sizeof(struct SEQMOTIF));
	if (pSeqMtf == NULL)
	{
		printf("Error: cannot create seq/motif complex.\n");
		exit(EXIT_FAILURE);
	}
	
	/* init */
	pSeqMtf->nId = nId;
	pSeqMtf->nSeqNum = nSeqNum;
	if(nSeqNum > 0)
	{
		pSeqMtf->pSeqId = (int *)calloc(nSeqNum, sizeof(int));
		pSeqMtf->ppSeq = (struct BYTEMATRIX **)calloc(nSeqNum, sizeof(struct BYTEMATRIX *));
	}

	pSeqMtf->nScoreNum = nScoreNum;
	if(nScoreNum > 0)
	{
		pSeqMtf->pScoreId = (int *)calloc(nScoreNum, sizeof(int));
		pSeqMtf->ppScore = (struct DOUBLEMATRIX **)calloc(nScoreNum, sizeof(struct DOUBLEMATRIX *));
	}

	pSeqMtf->nPathNum = nPathNum;
	if(nPathNum > 0)
	{
		pSeqMtf->pPathId = (int *)calloc(nPathNum, sizeof(int));
		pSeqMtf->ppSamplePath = (struct BYTEMATRIX **)calloc(nPathNum, sizeof(struct BYTEMATRIX *));
	}

	pSeqMtf->pHMM = NULL;
	pSeqMtf->pMotifSite = NULL;

	/* return */
	return pSeqMtf;
}

/* ----------------------------------------------------------------------- */ 
/*                  SEQMTFDESTROY: destroy seq/motif complex               */
/* ----------------------------------------------------------------------- */ 
int SEQMTFDESTROY(struct SEQMOTIF *pSeqMtf)
{
	int ni;

	/* initial check */
	if(pSeqMtf == NULL)
		return PROC_SUCCESS;

	/* destroy sequence */
	for(ni=0; ni<pSeqMtf->nSeqNum; ni++)
	{
		DestroyByteMatrix((pSeqMtf->ppSeq)[ni]);
		(pSeqMtf->ppSeq)[ni] = NULL;
	}
	free(pSeqMtf->ppSeq);
	pSeqMtf->ppSeq = NULL;
	free(pSeqMtf->pSeqId);
	pSeqMtf->pSeqId = NULL;
	pSeqMtf->nSeqNum = 0;

	/* destroy score */
	for(ni=0; ni<pSeqMtf->nScoreNum; ni++)
	{
		DestroyDoubleMatrix((pSeqMtf->ppScore)[ni]);
		(pSeqMtf->ppScore)[ni] = NULL;
	}
	free(pSeqMtf->ppScore);
	pSeqMtf->ppScore = NULL;
	free(pSeqMtf->pScoreId);
	pSeqMtf->pScoreId = NULL;
	pSeqMtf->nScoreNum = 0;

	/* destroy path */
	for(ni=0; ni<pSeqMtf->nPathNum; ni++)
	{
		DestroyByteMatrix((pSeqMtf->ppSamplePath)[ni]);
		(pSeqMtf->ppSamplePath)[ni] = NULL;
	}
	free(pSeqMtf->ppSamplePath);
	pSeqMtf->ppSamplePath = NULL;
	free(pSeqMtf->pPathId);
	pSeqMtf->pPathId = NULL;
	pSeqMtf->nPathNum = 0;

	/* destroy motif sites */
	MTFSDESTROYLIST(&(pSeqMtf->pMotifSite));
	pSeqMtf->nSiteNum = 0;

	/* destroy HMM */
	DestroyDoubleMatrix(pSeqMtf->pHMM);
	pSeqMtf->pHMM = NULL;

	/* free the parent structure */
	free(pSeqMtf);
	pSeqMtf = NULL;

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*              SEQMTFWRITESEQTOFASTA: write sequence to fasta file        */
/* ----------------------------------------------------------------------- */ 
int SEQMTFWRITESEQTOFASTA(FILE *fpOut, struct SEQMOTIF *pSeqMtf)
{
	int ni,nj,nk;
	unsigned char *pBase;

	/* init check */
	if((pSeqMtf == NULL) || (fpOut == NULL))
	{
		return PROC_FAILURE;
	}

	/* write */
	for(ni=0; ni<pSeqMtf->nSeqNum; ni++)
	{
		fprintf(fpOut, ">%d|%d\n", pSeqMtf->nId, pSeqMtf->pSeqId[ni]);
		if(pSeqMtf->ppSeq[ni] == NULL)
			continue;
		nk = 0;
		pBase = pSeqMtf->ppSeq[ni]->pMatElement;
		for(nj=0; nj<pSeqMtf->ppSeq[ni]->nWidth; nj++)
		{
			switch(*pBase)
			{
				case 0: fprintf(fpOut, "A");
					break;
				case 1: fprintf(fpOut, "C");
					break;
				case 2: fprintf(fpOut, "G");
					break;
				case 3: fprintf(fpOut, "T");
					break;
				default: fprintf(fpOut, "N");
			}

			nk++;
			pBase++;

			if(nk == FASTA_LINE_LEN)
			{
				fprintf(fpOut, "\n");
				nk = 0;
			}
		}
		fprintf(fpOut, "\n");
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*          SEQMTFWRITESAMPPATHTOFASTA: write sample path to fasta file    */
/* ----------------------------------------------------------------------- */ 
int SEQMTFWRITESAMPPATHTOFASTA(FILE *fpOut, struct SEQMOTIF *pSeqMtf)
{
	int ni,nj,nk;
	unsigned char *pBase;

	/* init check */
	if((pSeqMtf == NULL) || (fpOut == NULL))
	{
		return PROC_FAILURE;
	}

	/* write */
	for(ni=0; ni<pSeqMtf->nPathNum; ni++)
	{
		fprintf(fpOut, ">seqmtf%d|path%d\n", pSeqMtf->nId, pSeqMtf->pPathId[ni]);
		if(pSeqMtf->ppSamplePath[ni] == NULL)
			continue;
		nk = 0;
		pBase = pSeqMtf->ppSamplePath[ni]->pMatElement;
		for(nj=0; nj<pSeqMtf->ppSamplePath[ni]->nWidth; nj++)
		{
			fprintf(fpOut, "%d", (int)(*pBase));
			nk++;
			pBase++;

			if(nk == FASTA_LINE_LEN)
			{
				fprintf(fpOut, "\n");
				nk = 0;
			}
		}
		fprintf(fpOut, "\n");
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*          SEQMTFCREATESEQ: allocate memory required for a sequence       */
/* ----------------------------------------------------------------------- */ 
int SEQMTFCREATESEQ(struct SEQMOTIF *pSeqMtf, int nIndex, int nLen)
{
	/* check */
	if(pSeqMtf == NULL)
		return PROC_FAILURE;
	if( (nIndex >= pSeqMtf->nSeqNum) || (nIndex < 0) )
	{
		printf("Error: SEQMTFCREATESEQ, nIndex out of range!\n");
		exit(EXIT_FAILURE);
	}
	if(pSeqMtf->ppSeq[nIndex] != NULL)
	{
		DestroyByteMatrix(pSeqMtf->ppSeq[nIndex]);
		pSeqMtf->ppSeq[nIndex] = NULL;
	}

	/* create */
	pSeqMtf->ppSeq[nIndex] = CreateByteMatrix(1, nLen);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*          SEQMTFDESTROYSEQ: release memory required for a sequence       */
/* ----------------------------------------------------------------------- */ 
int SEQMTFDESTROYSEQ(struct SEQMOTIF *pSeqMtf, int nIndex)
{
	/* check */
	if(pSeqMtf == NULL)
		return PROC_FAILURE;
	if( (nIndex >= pSeqMtf->nSeqNum) || (nIndex < 0) )
	{
		printf("Error: SEQMTFDESTROYSEQ, nIndex out of range!\n");
		exit(EXIT_FAILURE);
	}
	if(pSeqMtf->ppSeq[nIndex] != NULL)
	{
		DestroyByteMatrix(pSeqMtf->ppSeq[nIndex]);
		pSeqMtf->ppSeq[nIndex] = NULL;
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*          SEQMTFCREATESCORE: allocate memory required for a score        */
/* ----------------------------------------------------------------------- */ 
int SEQMTFCREATESCORE(struct SEQMOTIF *pSeqMtf, int nIndex, int nLen)
{
	/* check */
	if(pSeqMtf == NULL)
		return PROC_FAILURE;
	if( (nIndex >= pSeqMtf->nScoreNum) || (nIndex < 0) )
	{
		printf("Error: SEQMTFCREATESCORE, nIndex out of range!\n");
		exit(EXIT_FAILURE);
	}
	if(pSeqMtf->ppScore[nIndex] != NULL)
	{
		DestroyDoubleMatrix(pSeqMtf->ppScore[nIndex]);
		pSeqMtf->ppScore[nIndex] = NULL;
	}

	/* create */
	pSeqMtf->ppScore[nIndex] = CreateDoubleMatrix(1, nLen);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*          SEQMTFDESTROYSCORE: release memory required for a score        */
/* ----------------------------------------------------------------------- */ 
int SEQMTFDESTROYSCORE(struct SEQMOTIF *pSeqMtf, int nIndex)
{
	/* check */
	if(pSeqMtf == NULL)
		return PROC_FAILURE;
	if( (nIndex >= pSeqMtf->nScoreNum) || (nIndex < 0) )
	{
		printf("Error: SEQMTFDESTROYSCORE, nIndex out of range!\n");
		exit(EXIT_FAILURE);
	}
	if(pSeqMtf->ppScore[nIndex] != NULL)
	{
		DestroyDoubleMatrix(pSeqMtf->ppScore[nIndex]);
		pSeqMtf->ppScore[nIndex] = NULL;
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*          SEQMTFCREATEPATH: allocate memory required for a sample path   */
/* ----------------------------------------------------------------------- */ 
int SEQMTFCREATEPATH(struct SEQMOTIF *pSeqMtf, int nIndex, int nLen)
{
	/* check */
	if(pSeqMtf == NULL)
		return PROC_FAILURE;
	if( (nIndex >= pSeqMtf->nPathNum) || (nIndex < 0) )
	{
		printf("Error: SEQMTFCREATEPATH, nIndex out of range!\n");
		exit(EXIT_FAILURE);
	}
	if(pSeqMtf->ppSamplePath[nIndex] != NULL)
	{
		DestroyByteMatrix(pSeqMtf->ppSamplePath[nIndex]);
		pSeqMtf->ppSamplePath[nIndex] = NULL;
	}

	/* create */
	pSeqMtf->ppSamplePath[nIndex] = CreateByteMatrix(1, nLen);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*          SEQMTFDESTROYPATH: release memory required for a sample path   */
/* ----------------------------------------------------------------------- */ 
int SEQMTFDESTROYPATH(struct SEQMOTIF *pSeqMtf, int nIndex)
{
	/* check */
	if(pSeqMtf == NULL)
		return PROC_FAILURE;
	if( (nIndex >= pSeqMtf->nPathNum) || (nIndex < 0) )
	{
		printf("Error: SEQMTFDESTROYPATH, nIndex out of range!\n");
		exit(EXIT_FAILURE);
	}
	if(pSeqMtf->ppSamplePath[nIndex] != NULL)
	{
		DestroyByteMatrix(pSeqMtf->ppSamplePath[nIndex]);
		pSeqMtf->ppSamplePath[nIndex] = NULL;
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*          SEQMTFCREATEHMM: allocate memory required for HMM.             */
/* ----------------------------------------------------------------------- */ 
int SEQMTFCREATEHMM(struct SEQMOTIF *pSeqMtf, int nStatusNum, int nLen)
{
	/* check */
	if(pSeqMtf == NULL)
		return PROC_FAILURE;
	if( (nStatusNum <= 0) || (nLen <= 0) )
	{
		printf("Error: SEQMTFCREATEHMM, nStatusNum and nLen must be greater than 0!\n");
		exit(EXIT_FAILURE);
	}
	if(pSeqMtf->pHMM != NULL)
	{
		DestroyDoubleMatrix(pSeqMtf->pHMM);
		pSeqMtf->pHMM = NULL;
	}

	/* create */
	pSeqMtf->pHMM = CreateDoubleMatrix(nStatusNum, nLen);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*          SEQMTFDESTROYHMM: release memory required for HMM.             */
/* ----------------------------------------------------------------------- */ 
int SEQMTFDESTROYHMM(struct SEQMOTIF *pSeqMtf)
{
	/* check */
	if(pSeqMtf == NULL)
		return PROC_FAILURE;
	
	if(pSeqMtf->pHMM != NULL)
	{
		DestroyDoubleMatrix(pSeqMtf->pHMM);
		pSeqMtf->pHMM = NULL;
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*          SEQMTFCODENUCLEICSEQ: code nucleic acid sequences.             */
/* ----------------------------------------------------------------------- */ 
int SEQMTFCODENUCLEICSEQ(struct SEQMOTIF *pSeqMtf, int nIndex, int nLen, char strSeq[])
{
	int ni;
	unsigned char *pEle;

	/* check */
	if(pSeqMtf == NULL)
		return PROC_FAILURE;
	if( (nIndex >= pSeqMtf->nSeqNum) || (nIndex < 0) )
	{
		printf("Error: SEQMTFCODENUCLEICSEQ, nIndex out of range!\n");
		exit(EXIT_FAILURE);
	}

	if( ((int)(strlen(strSeq)) != nLen) || (pSeqMtf->ppSeq[nIndex]->nWidth != nLen) )
	{
		printf("Error: SEQMTFCODENUCLEICSEQ, lengths not match!\n");
		exit(EXIT_FAILURE);
	}

	/* code */
	pEle = pSeqMtf->ppSeq[nIndex]->pMatElement;
	for(ni=0; ni<nLen; ni++)
	{
		switch(strSeq[ni])
		{
			case 'a': *pEle = 0;
				break;
			case 'A': *pEle = 0;
				break;
			case 'c': *pEle = 1;
				break;
			case 'C': *pEle = 1;
				break;
			case 'g': *pEle = 2;
				break;
			case 'G': *pEle = 2;
				break;
			case 't': *pEle = 3;
				break;
			case 'T': *pEle = 3;
				break;
			default: *pEle = 4;
		}
		pEle++;
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*       SEQMTFLOADPRIMARYSEQFROMFASTA: load sequence from fasta file.     */
/* Return a vector of pointers to newly created motif/seq complexes, and   */
/* nCount, the number of loaded complexes.                                 */
/* ----------------------------------------------------------------------- */
struct SEQMOTIF **SEQMTFLOADPRIMARYSEQFROMFASTA(int *pCount, char strFilePath[])
{
	/* define */
	int nCount;
	struct tagSequence *pSeqList;
	struct tagSequence *pCurrentSeq;
	struct SEQMOTIF **pSeqMtf;
	int ni,nLen;

	/* init */
	pSeqMtf = NULL;

	/* load sequences from fasta file */
	pSeqList = NULL;
	nCount = LoadFullSequenceList(strFilePath, &pSeqList);
	if(nCount <= 0)
		return NULL;

	/* create complexes */
	pSeqMtf = (struct SEQMOTIF **)calloc(nCount, sizeof(struct SEQMOTIF *));
	if(pSeqMtf == NULL)
		return NULL;
	/* code sequence */
	pCurrentSeq = pSeqList;
	for(ni=0; ni<nCount; ni++)
	{
		pSeqMtf[ni] = NULL;
		pSeqMtf[ni] = SEQMTFCREATE(ni, 1, 1, 1);
		if(pSeqMtf[ni] == NULL)
		{
			printf("Error: SEQMTFLOADPRIMARYSEQFROMFASTA, cannot create sequence/motif complex!\n");
			exit(EXIT_FAILURE);
		}
		nLen = pCurrentSeq->m_nLength;
		SEQMTFCREATESEQ(pSeqMtf[ni], 0, nLen);
		SEQMTFCREATESCORE(pSeqMtf[ni], 0, nLen);
		SEQMTFCREATEPATH(pSeqMtf[ni], 0, nLen);
		SEQMTFCODENUCLEICSEQ(pSeqMtf[ni], 0, nLen, pCurrentSeq->m_pSequence->m_pString);

		/* get next seq */
		pCurrentSeq = pCurrentSeq->m_pNext;
	}

	/* destroy fasta sequences */
	SequenceListClear(&pSeqList);

	/* return */
	*pCount = nCount;
	return pSeqMtf;
}

/* ----------------------------------------------------------------------- */ 
/* SEQMTFLOADPRIMARYSEQFROMFASTA_FORQUALITY: load sequence from fasta file */
/* Return a vector of pointers to newly created motif/seq complexes, and   */
/* nCount, the number of loaded complexes.                                 */
/* ----------------------------------------------------------------------- */
struct SEQMOTIF **SEQMTFLOADPRIMARYSEQFROMFASTA_FORQUALITY(int *pCount, char strFilePath[])
{
	/* define */
	int nCount;
	struct tagSequence *pSeqList;
	struct tagSequence *pCurrentSeq;
	struct SEQMOTIF **pSeqMtf;
	int ni,nj,nLen;

	/* init */
	pSeqMtf = NULL;

	/* load sequences from fasta file */
	pSeqList = NULL;
	nCount = LoadFullSequenceList(strFilePath, &pSeqList);
	if(nCount <= 0)
		return NULL;

	/* create complexes */
	pSeqMtf = (struct SEQMOTIF **)calloc(nCount, sizeof(struct SEQMOTIF *));
	if(pSeqMtf == NULL)
		return NULL;
	/* code sequence */
	pCurrentSeq = pSeqList;
	for(ni=0; ni<nCount; ni++)
	{
		pSeqMtf[ni] = NULL;
		pSeqMtf[ni] = SEQMTFCREATE(ni, 1, 0, 3);
		if(pSeqMtf[ni] == NULL)
		{
			printf("Error: SEQMTFLOADPRIMARYSEQFROMFASTA_FORQUALITY, cannot create sequence/motif complex!\n");
			exit(EXIT_FAILURE);
		}
		nLen = pCurrentSeq->m_nLength;
		SEQMTFCREATESEQ(pSeqMtf[ni], 0, nLen);
		for(nj=0; nj<3; nj++)
		{
			pSeqMtf[ni]->pPathId[nj] = nj;
			SEQMTFCREATEPATH(pSeqMtf[ni], nj, nLen);
		}
		SEQMTFCODENUCLEICSEQ(pSeqMtf[ni], 0, nLen, pCurrentSeq->m_pSequence->m_pString);

		/* get next seq */
		pCurrentSeq = pCurrentSeq->m_pNext;
	}

	/* destroy fasta sequences */
	SequenceListClear(&pSeqList);

	/* return */
	*pCount = nCount;
	return pSeqMtf;
}


/* ----------------------------------------------------------------------- */ 
/* SEQMTFLOADPRIMARYSEQFROMFASTA_FORREGRESSOR: load sequence from fasta    */
/* file and corresponding conservation scores.                             */
/* Return a vector of pointers to newly created motif/seq complexes, and   */
/* nCount, the number of loaded complexes.                                 */
/* ----------------------------------------------------------------------- */
struct SEQMOTIF **SEQMTFLOADPRIMARYSEQFROMFASTA_FORREGRESSOR(int *pCount, char strFilePath[], char strFastaName[])
{
	/* define */
	int nCount;
	struct tagSequence *pSeqList;
	struct tagSequence *pCurrentSeq;
	struct SEQMOTIF **pSeqMtf;
	int ni,nLen;
	char strConsName[LINE_LENGTH];
	FILE *fpCons;
	int nReadCount;
	char strPath[LINE_LENGTH];

	/* init */
	pSeqMtf = NULL;
	sprintf(strPath, strFilePath);
	StrTrimRight(strPath);
	ni = strlen(strPath);
	if(strPath[ni-1] == '\\')
		strPath[ni-1] = '\0';

	sprintf(strConsName, "%s\\%s", strPath, strFastaName);

	/* load sequences from fasta file */
	pSeqList = NULL;
	nCount = LoadFullSequenceList(strConsName, &pSeqList);
	if(nCount <= 0)
		return NULL;

	/* create complexes */
	pSeqMtf = (struct SEQMOTIF **)calloc(nCount, sizeof(struct SEQMOTIF *));
	if(pSeqMtf == NULL)
		return NULL;
	/* code sequence */
	pCurrentSeq = pSeqList;
	for(ni=0; ni<nCount; ni++)
	{
		pSeqMtf[ni] = NULL;
		pSeqMtf[ni] = SEQMTFCREATE(ni, 1, 0, 1);
		if(pSeqMtf[ni] == NULL)
		{
			printf("Error: SEQMTFLOADPRIMARYSEQFROMFASTA_FORREGRESSOR, cannot create sequence/motif complex!\n");
			exit(EXIT_FAILURE);
		}
		nLen = pCurrentSeq->m_nLength;
		SEQMTFCREATESEQ(pSeqMtf[ni], 0, nLen);
		SEQMTFCREATEPATH(pSeqMtf[ni], 0, nLen);
		SEQMTFCODENUCLEICSEQ(pSeqMtf[ni], 0, nLen, pCurrentSeq->m_pSequence->m_pString);
		sprintf(strConsName, "%s\\%s.cs", strPath, pCurrentSeq->m_strAlias); 
		if( (fpCons = fopen(strConsName, "rb" )) != NULL )
		{
			  /* load conservation */
			  nReadCount = fread( pSeqMtf[ni]->ppSamplePath[0]->pMatElement, sizeof( unsigned char ), nLen, fpCons );
      
			  if(nReadCount != nLen)
			  {
				  printf("Error: load conservation score failure!\n");
				  fclose( fpCons );
				  exit(EXIT_FAILURE);
			  }

			  /*for(nj=0; nj<25; nj++)
			  {
				  printf("%d\n", (int)pSeqMtf[ni]->ppSamplePath[0]->pMatElement[nj]);
			  }*/

			  fclose( fpCons );
		}
		else
		{
			printf("Error: load conservation score failure!\n");
			exit(EXIT_FAILURE);
		}

		/* get next seq */
		pCurrentSeq = pCurrentSeq->m_pNext;
	}

	/* destroy fasta sequences */
	SequenceListClear(&pSeqList);

	/* return */
	*pCount = nCount;
	return pSeqMtf;
}


/* ----------------------------------------------------------------------- */ 
/*       SEQMTFESTIMATENUCLEICBGMC: estimate markov background.            */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *SEQMTFESTIMATENUCLEICBGMC(struct SEQMOTIF **pSeqMtf, int nCount, int nIndex, int nOrder)
{
	/* define */
	struct DOUBLEMATRIX *pBG;
	int nHeight,nWidth,ni,nj,nk,nz,nLen;
	int nWordId,nBadLen;
	unsigned char *pEle;
	int nSLabel[BGMC_MAX_ORDER];
	double dTemp;
	double *pBaseCount,*pBaseStart;
	double dTotalCount;
	int nScale;

	/* check */
	if((pSeqMtf == NULL) || (nCount <= 0) || (nOrder < 0) || (nIndex < 0) || (nOrder > BGMC_MAX_ORDER))
	{
		printf("Error: SEQMTFESTIMATENUCLEICBGMC, parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	
	/* init */
	pBG = NULL;
	nWidth = 4;
	nHeight = (int)(pow((double)nWidth, (double)nOrder));
	nScale = (int)(pow((double)nWidth, (double)(nOrder-1)));
	pBG = CreateDoubleMatrix(nHeight, nWidth);
	if(pBG == NULL)
	{
		printf("Error: SEQMTFESTIMATENUCLEICBGMC, can't create background matrix!\n");
		exit(EXIT_FAILURE);
	}
	pBaseStart = pBG->pMatElement;
	for(ni=0; ni<nHeight; ni++)
	{
		for(nj=0; nj<nWidth; nj++)
		{
			*pBaseStart = 1.0;
			pBaseStart++;
		}
	}

	/* create */
	for(ni=0; ni<nCount; ni++)
	{
		if(pSeqMtf[ni] == NULL)
			continue;
		if(nIndex >= pSeqMtf[ni]->nSeqNum)
			continue;
		if(pSeqMtf[ni]->ppSeq[nIndex] == NULL)
			continue;

		for(nj=0; nj<BGMC_MAX_ORDER; nj++)
			nSLabel[nj] = 0;
		nLen = pSeqMtf[ni]->ppSeq[nIndex]->nWidth;
		
		if(nOrder == 0)
		{
			pEle = pSeqMtf[ni]->ppSeq[nIndex]->pMatElement;
			for(nj=0; nj<nLen; nj++)
			{
				nWordId = (int)(*pEle);
				if(nWordId <nWidth)
				{
					dTemp = DMGETAT(pBG, 0, nWordId)+1.0;
					DMSETAT(pBG, 0, nWordId, dTemp);
				}
				pEle++;
			}
		}
		else
		{
			nWordId = 0;
			nBadLen = 0;
			pEle = pSeqMtf[ni]->ppSeq[nIndex]->pMatElement;
			for(nj=0; nj<nOrder; nj++)
			{
				nSLabel[nj] = (int)(*pEle);
				nk = (int)(*pEle);
				if(nk < nWidth)
				{
					nWordId = nWidth*nWordId+nk;
				}
				else
				{
					nWordId = nWidth*nWordId;
					nBadLen++;
				}
				
				pEle++;
			}
			for(; nj<nLen; nj++)
			{
				nk = (int)(*pEle);
				if( (nk < nWidth) && (nBadLen == 0) )
				{
					dTemp = DMGETAT(pBG, nWordId, nk)+1.0;
					DMSETAT(pBG, nWordId, nk, dTemp);
				}

				if(nSLabel[0] < nWidth)
				{
					nWordId -= nSLabel[0]*nScale;
				}
				else
				{
					nBadLen--;
				}
				for(nz=0; nz<(nOrder-1); nz++)
				{
					nSLabel[nz] = nSLabel[nz+1];
				}
				nSLabel[nz] = nk;
				
				if(nk<nWidth)
				{
					nWordId = nWordId*nWidth+nk;
				}
				else
				{
					nWordId = nWordId*nWidth;
					nBadLen++;
				}

				/* get next */
				pEle++;
			}
		}
	}

	/* normalize */
	pBaseStart = pBG->pMatElement;
	for(ni=0; ni<nHeight; ni++)
	{
		pBaseCount = pBaseStart;
		dTotalCount = 0.0;
		for(nj=0; nj<nWidth; nj++)
		{
			dTotalCount += (*pBaseCount);
			pBaseCount++;
		}

		if(dTotalCount>0.0)
		{
			for(nj=0; nj<nWidth; nj++)
			{
				*pBaseStart = *pBaseStart/dTotalCount;
				pBaseStart++;
			}
		}
		else
		{
			for(nj=0; nj<nWidth; nj++)
			{
				*pBaseStart = 1.0/(double)nWidth;
				pBaseStart++;
			}
		}
	}

	/* return */
	return pBG;
}

/* ----------------------------------------------------------------------- */ 
/*       SEQMTFESTIMATENUCLEICBGMC: estimate markov background.            */
/* ----------------------------------------------------------------------- */ 
double SEQMTF_MOTIFMAP_SCORE(struct SEQMOTIF *pSeqMtf, int nIndex, 
							 int nBGOrder, struct DOUBLEMATRIX *pBG,
							 struct DOUBLEMATRIX *pMotif,
							 double dRatioCutoff)
{
	/* define */
	double dScore,dLocal,dCons;
	int nLen,ni,nj,nk,nz;
	int nWordId,nBadLen;
	int nLocalWordId,nLocalBadLen;
	int nTotalCount;
	unsigned char *pBase;
	unsigned char *pCons;
	int nWidth,nScale;
	int nMotifLen;
	int nSkip;
	double dRatio;

	/* init */
	dScore = 0.0;
	nWordId = 0;
	nBadLen = 0;
	nTotalCount = 0;
	nWidth = pBG->nWidth;
	nMotifLen = pMotif->nHeight;
	if(nBGOrder == 0)
		nScale = 0;
	else
		nScale = (int)pow((double)nWidth, (double)(nBGOrder-1));
	
	/* get motif score*/
	nLen = pSeqMtf->ppSeq[nIndex]->nWidth;
	if(nBGOrder == 0)
	{
		pBase = pSeqMtf->ppSeq[nIndex]->pMatElement;
		pCons = pSeqMtf->ppSamplePath[0]->pMatElement;
		dLocal = 0.0;
		for(ni=0; ni<nLen-nMotifLen+1; ni++)
		{
			/* get conservation score */
			dCons = 0.0;
			for(nj=0; nj<nMotifLen; nj++)
			{
				dCons += (double)pCons[nj];
			}
			dCons /= (double)nMotifLen;

			/* if not conserved, skip */
			if(dCons < MOTIFMAP_CONSERVE_CUTOFF)
			{
				/* get next */
				pBase++;
				pCons++;
				continue;
			}

			/* otherwise process it */
			dLocal = 0.0;
			nBadLen = 0;
			for(nj=0; nj<nMotifLen; nj++)
			{
				nk = (int)(pBase[nj]);
				if(nk >= nWidth)
				{
					nBadLen++;
					break;
				}
				dLocal += DMGETAT(pMotif, nj, nk)-DMGETAT(pBG, 0, nk);
			}

			if(nBadLen > 0)
			{
				/* get next */
				pBase++;
				pCons++;
				continue;
			}

			/* add up to total score */
			dRatio = exp(dLocal);
			if(dRatio >= dRatioCutoff)
			{
				dScore += dRatio;
				nTotalCount++;
			}

			/* get next */
			pBase++;
			pCons++;
		}

		/* get average */
		dScore = log(dScore);
	}

	/* if using markov background */
	else
	{
		pBase = pSeqMtf->ppSeq[nIndex]->pMatElement;
		pCons = pSeqMtf->ppSamplePath[0]->pMatElement;
	
		nWordId = 0;
		nBadLen = 0;
		for(ni=0; ni<nBGOrder; ni++)
		{
			nk = (int)(*pBase);
			if(nk < nWidth)
			{
				nWordId = nWidth*nWordId+nk;
			}
			else
			{
				nWordId = nWidth*nWordId;
				nBadLen++;
			}
			
			/* get next */
			pBase++;
			pCons++;
		}

		for(; ni<nLen-nMotifLen+1; ni++)
		{
			dLocal = 0.0;
			nSkip = 0;
			nLocalWordId = nWordId;
			nLocalBadLen = nBadLen;
			for(nj=0; nj<nMotifLen; nj++)
			{
				nk = (int)(pBase[nj]);
				if((nk >= nWidth) || (nLocalBadLen > 0))
				{
					nSkip = 1;
					break;
				}
				dLocal += DMGETAT(pMotif, nj, nk)-DMGETAT(pBG, nLocalWordId, nk);
				
				nz = (int)pBase[nj-nBGOrder];
				if(nz < nWidth)
				{
					nLocalWordId -= nz*nScale;
				}
				else
				{
					nLocalBadLen--;
				}

				if(nk < nWidth)
				{
					nLocalWordId = nWidth*nLocalWordId+nk;
				}
				else
				{
					nLocalWordId = nWidth*nLocalWordId;
					nLocalBadLen++;
				}
			}

			/* add up to total score */
			if(nSkip == 0)
			{
				dRatio = exp(dLocal);
				if(dRatio >= dRatioCutoff)
				{
					dScore += dRatio;
					nTotalCount++;
				}
			}

			nk = *pBase;
			nz = *(pBase-nBGOrder);
			if(nz < nWidth)
			{
				nWordId -= nz*nScale;
			}
			else
			{
				nBadLen--;
			}
			if(nk < nWidth)
			{
				nWordId = nWidth*nWordId+nk;
			}
			else
			{
				nWordId = nWidth*nWordId;
				nBadLen++;
			}
			
			/* get next */
			pBase++;
			pCons++;
		}
	}

	/* return */
	return dScore;
}

/* ----------------------------------------------------------------------- */ 
/*      SEQMTFADDCOUNTTONUCLEICBGMC: add count to markov background.       */
/* ----------------------------------------------------------------------- */ 
int SEQMTFADDCOUNTTONUCLEICBGMC(struct DOUBLEMATRIX *pBG, struct SEQMOTIF *pSeqMtf, int nIndex, int nOrder)
{
	/* define */
	int nHeight,nWidth,nScale;
	int nSLabel[BGMC_MAX_ORDER];
	int nz,nj,nk,nLen,nBadLen;
	unsigned char *pEle;
	int nWordId;
	double dTemp;

	/* check */
	if(pBG == NULL)
	{
		printf("Error: SEQMTFADDCOUNTTONUCLEICBGMC, no background matrix!\n");
		exit(EXIT_FAILURE);
	}
	if(pSeqMtf == NULL)
	{
		printf("Error: SEQMTFADDCOUNTTONUCLEICBGMC, no sequence!\n");
		exit(EXIT_FAILURE);
	}
	if(nIndex >= pSeqMtf->nSeqNum)
	{
		printf("Error: SEQMTFADDCOUNTTONUCLEICBGMC, index out of range!\n");
		exit(EXIT_FAILURE);
	}
	if(pSeqMtf->ppSeq[nIndex] == NULL)
	{
		printf("Error: SEQMTFADDCOUNTTONUCLEICBGMC, no sequence!\n");
		exit(EXIT_FAILURE);
	}
	
	
	/* init */
	nWidth = pBG->nWidth;
	nHeight = (int)(pow((double)nWidth, (double)nOrder));
	if(nHeight != pBG->nHeight)
	{
		printf("Error: SEQMTFADDCOUNTTONUCLEICBGMC, dimension not match!\n");
		exit(EXIT_FAILURE);
	}

	nScale = (int)(pow((double)nWidth, (double)(nOrder-1)));

	/* create */
	for(nj=0; nj<BGMC_MAX_ORDER; nj++)
		nSLabel[nj] = 0;
	nLen = pSeqMtf->ppSeq[nIndex]->nWidth;
	
	if(nOrder == 0)
	{
		pEle = pSeqMtf->ppSeq[nIndex]->pMatElement;
		for(nj=0; nj<nLen; nj++)
		{
			nWordId = (int)(*pEle);
			if(nWordId < nWidth)
			{
				dTemp = DMGETAT(pBG, 0, nWordId)+1.0;
				DMSETAT(pBG, 0, nWordId, dTemp);
			}
			pEle++;
		}
	}
	else
	{
		nWordId = 0;
		nBadLen = 0;
		pEle = pSeqMtf->ppSeq[nIndex]->pMatElement;
		for(nj=0; nj<nOrder; nj++)
		{
			nSLabel[nj] = (int)(*pEle);
			nk = (int)(*pEle);
			if(nk < nWidth)
			{
				nWordId = nWidth*nWordId+nk;
			}
			else
			{
				nWordId = nWidth*nWordId;
				nBadLen++;
			}
			
			pEle++;
		}
		for(; nj<nLen; nj++)
		{
			nk = (int)(*pEle);
			if( (nk < nWidth) && (nBadLen == 0) )
			{
				dTemp = DMGETAT(pBG, nWordId, nk)+1.0;
				DMSETAT(pBG, nWordId, nk, dTemp);
			}

			if(nSLabel[0] < nWidth)
			{
				nWordId -= nSLabel[0]*nScale;
			}
			else
			{
				nBadLen--;
			}
			for(nz=0; nz<(nOrder-1); nz++)
			{
				nSLabel[nz] = nSLabel[nz+1];
			}
			nSLabel[nz] = nk;
			
			if(nk<nWidth)
			{
				nWordId = nWordId*nWidth+nk;
			}
			else
			{
				nWordId = nWordId*nWidth;
				nBadLen++;
			}

			/* get next */
			pEle++;
		}
	}
	

	
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*      SEQMTF_LOADSEQFROMGENOME: load sequences from genome.              */
/* This function will load sequences from fpSeq (from nStart base to nEnd  */
/* base), and write it to the nSeqId th sequence of pSeqMtf, starting from */
/* nOffset th bp.                                                          */
/* ----------------------------------------------------------------------- */ 
int SEQMTF_LOADSEQFROMGENOME(struct SEQMOTIF *pSeqMtf, int nSeqId, int nOffset,
							 FILE *fpSeq, int nStart, int nEnd, char chStrand)
{
	/* define */
	unsigned char *pBase;
	int nLen;
	int nP1,nP2,nR1,nR2;
	int numread;
	int ni,nk;
	unsigned char bBase,bChar;

	/* init */
	if( (pSeqMtf == NULL) || (fpSeq == NULL) )
	{
		printf("Warning: SEQMTF_LOADSEQFROMGENOME, null seqmotif or seqfile pointer!\n");
		return PROC_FAILURE;
	}
	if( (nSeqId<0) || (nSeqId >= pSeqMtf->nSeqNum) )
	{
		printf("Warning: SEQMTF_LOADSEQFROMGENOME, sequence id out of range!\n");
		return PROC_FAILURE;
	}
	if( (nOffset<0) || (nOffset >= pSeqMtf->ppSeq[nSeqId]->nWidth) )
	{
		printf("Warning: SEQMTF_LOADSEQFROMGENOME, offset out of range!\n");
		return PROC_FAILURE;
	}
	nLen = (nEnd-nStart+1);
	if( (nLen<=0) || (nLen >  (pSeqMtf->ppSeq[nSeqId]->nWidth-nOffset) ) )
	{
		printf("Warning: SEQMTF_LOADSEQFROMGENOME, the length of the sequence to be retrieved is larger than the length allowed by seqmotif!\n");
		return PROC_FAILURE;
	}
	if((chStrand != '+') && (chStrand != '-'))
	{
		printf("Warning: SEQMTF_LOADSEQFROMGENOME, strand information not specified!\n");
	}


	/* get seq */
	pBase = pSeqMtf->ppSeq[nSeqId]->pMatElement+nOffset;
	nP1 = (int)(nStart/2);
	nR1 = nStart%2;
	nP2 = (int)(nEnd/2);
	nR2 = nEnd%2;
	nk = 0;

	if( fseek( fpSeq, nP1, SEEK_SET ) != 0 )
	{
		printf("Error: SEQMTF_LOADSEQFROMGENOME, cannot locate the required sequence!\n");
		exit(EXIT_FAILURE);
	}


	/* first base */
	numread = fread( &bChar, sizeof( unsigned char ), 1, fpSeq );
	if(nR1 == 0)
	{
		bBase = bChar >> 4;
		switch(bBase)
		{
			case 0: *pBase = 0;
				break;
			case 1: *pBase = 1;
				break;
			case 2: *pBase = 2;
				break;
			case 3: *pBase = 3;
				break;
			case 4: *pBase = 0;
				break;
			case 5: *pBase = 1;
				break;
			case 6: *pBase = 2;
				break;
			case 7: *pBase = 3;
				break;
			default: *pBase = 4;
		}
		pBase++;
		nk++;
	}

	if(nk < nLen)
	{
		bBase = bChar & 0x0F;
		switch(bBase)
		{
			case 0: *pBase = 0;
				break;
			case 1: *pBase = 1;
				break;
			case 2: *pBase = 2;
				break;
			case 3: *pBase = 3;
				break;
			case 4: *pBase = 0;
				break;
			case 5: *pBase = 1;
				break;
			case 6: *pBase = 2;
				break;
			case 7: *pBase = 3;
				break;
			default: *pBase = 4;
		}
		pBase++;
		nk++;
	}

	/* middle bases */
	if(nk < nLen)
	{
		for(ni=(nP1+1); ni<nP2; ni++)
		{
			numread = fread( &bChar, sizeof( unsigned char ), 1, fpSeq );
			bBase = bChar >> 4;
			switch(bBase)
			{
				case 0: *pBase = 0;
					break;
				case 1: *pBase = 1;
					break;
				case 2: *pBase = 2;
					break;
				case 3: *pBase = 3;
					break;
				case 4: *pBase = 0;
					break;
				case 5: *pBase = 1;
					break;
				case 6: *pBase = 2;
					break;
				case 7: *pBase = 3;
					break;
				default: *pBase = 4;
			}
			pBase++;
			nk++;

			bBase = bChar & 0x0F;
			switch(bBase)
			{
				case 0: *pBase = 0;
					break;
				case 1: *pBase = 1;
					break;
				case 2: *pBase = 2;
					break;
				case 3: *pBase = 3;
					break;
				case 4: *pBase = 0;
					break;
				case 5: *pBase = 1;
					break;
				case 6: *pBase = 2;
					break;
				case 7: *pBase = 3;
					break;
				default: *pBase = 4;
			}
			pBase++;
			nk++;
		}
	}
	
	/* last base */
	if(nk < nLen)
	{
		numread = fread( &bChar, sizeof( unsigned char ), 1, fpSeq );
		bBase = bChar >> 4;
		switch(bBase)
		{
			case 0: *pBase = 0;
				break;
			case 1: *pBase = 1;
				break;
			case 2: *pBase = 2;
				break;
			case 3: *pBase = 3;
				break;
			case 4: *pBase = 0;
				break;
			case 5: *pBase = 1;
				break;
			case 6: *pBase = 2;
				break;
			case 7: *pBase = 3;
				break;
			default: *pBase = 4;
		}
		pBase++;
		nk++;
	}

	if(nk < nLen)
	{
		if(nR2 == 1)
		{
			bBase = bChar & 0x0F;
			switch(bBase)
			{
				case 0: *pBase = 0;
					break;
				case 1: *pBase = 1;
					break;
				case 2: *pBase = 2;
					break;
				case 3: *pBase = 3;
					break;
				case 4: *pBase = 0;
					break;
				case 5: *pBase = 1;
					break;
				case 6: *pBase = 2;
					break;
				case 7: *pBase = 3;
					break;
				default: *pBase = 4;
			}
			pBase++;
			nk++;
		}
	}

	if(nk != nLen)
	{
		printf("Error: SEQMTF_LOADSEQFROMGENOME, sequence length not match!\n");
		exit(EXIT_FAILURE);
	}

	/* get reverse complement if necessary */
	if(chStrand == '-')
	{
		pBase = pSeqMtf->ppSeq[nSeqId]->pMatElement+nOffset;
		ni = 0;
		nk = nLen-1;
		while(nk>=ni)
		{
			bBase = pBase[ni];
			if(pBase[nk]<4)
			{
				pBase[ni] = 3-pBase[nk];
			}
			else
			{
				pBase[ni] = 4;
			}

			if(bBase<4)
			{
				pBase[nk] = 3-bBase;
			}
			else
			{
				pBase[nk] = 4;
			}

			ni++;
			nk--;
		}
	}

	/* return */
	return PROC_SUCCESS;
}


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
											int nFromPos, int nToPos, int nOrder)
{
	/* define */
	int nHeight,nWidth,nScale;
	int nFLabel[BGMC_MAX_ORDER];
	int nRLabel[BGMC_MAX_ORDER];

	int nz,nj,nk,nLen,nFBadLen,nRBadLen;
	unsigned char *pEle;
	int nFWordId,nRWordId,nWordId;
	double dTemp;

	/* check */
	if((pBGF == NULL) || (pBGR == NULL))
	{
		printf("Error: SEQMTF_ADDCOUNTTONUCLEICBGMC_BOTHSTRAND, no background matrix!\n");
		exit(EXIT_FAILURE);
	}
	if(pSeqMtf == NULL)
	{
		printf("Error: SEQMTF_ADDCOUNTTONUCLEICBGMC_BOTHSTRAND, no sequence!\n");
		exit(EXIT_FAILURE);
	}
	if((nIndex < 0) || (nIndex >= pSeqMtf->nSeqNum))
	{
		printf("Error: SEQMTF_ADDCOUNTTONUCLEICBGMC_BOTHSTRAND, index out of range!\n");
		exit(EXIT_FAILURE);
	}
	if(pSeqMtf->ppSeq[nIndex] == NULL)
	{
		printf("Error: SEQMTF_ADDCOUNTTONUCLEICBGMC_BOTHSTRAND, no sequence!\n");
		exit(EXIT_FAILURE);
	}
	if( (nEffectiveLen > pSeqMtf->ppSeq[nIndex]->nWidth) ||
		(nFromPos<nOrder) || (nFromPos>=nEffectiveLen) ||
		(nToPos<0) || (nToPos>=(nEffectiveLen-nOrder)) || (nFromPos>nToPos))
	{
		printf("Error: SEQMTF_ADDCOUNTTONUCLEICBGMC_BOTHSTRAND, sequence out of range!\n");
		exit(EXIT_FAILURE);
	}
	if(pBGF->nWidth != pBGR->nWidth)
	{
		printf("Error: SEQMTF_ADDCOUNTTONUCLEICBGMC_BOTHSTRAND, matrix dimension not match!\n");
		exit(EXIT_FAILURE);
	}

	
	/* init */
	nWidth = pBGF->nWidth;
	nHeight = (int)(pow((double)nWidth, (double)nOrder));
	if( (nHeight != pBGF->nHeight) || (nHeight != pBGR->nHeight) ) 
	{
		printf("Error: SEQMTF_ADDCOUNTTONUCLEICBGMC_BOTHSTRAND, matrix dimension not match!\n");
		exit(EXIT_FAILURE);
	}

	nScale = (int)(pow((double)nWidth, (double)(nOrder-1)));

	/* create */
	for(nj=0; nj<BGMC_MAX_ORDER; nj++)
	{
		nFLabel[nj] = 0;
		nRLabel[nj] = 0;
	}

	nLen = pSeqMtf->ppSeq[nIndex]->nWidth;
	
	if(nOrder == 0)
	{
		pEle = pSeqMtf->ppSeq[nIndex]->pMatElement;
		for(nj=nFromPos; nj<=nToPos; nj++)
		{
			nWordId = (int)(*pEle);
			if(nWordId < nWidth)
			{
				dTemp = DMGETAT(pBGF, 0, nWordId)+1.0;
				DMSETAT(pBGF, 0, nWordId, dTemp);
				dTemp = DMGETAT(pBGR, 0, nWordId)+1.0;
				DMSETAT(pBGR, 0, nWordId, dTemp);
			}
			pEle++;
		}
	}
	else
	{
		nFWordId = 0;
		nFBadLen = 0;
		nRWordId = 0;
		nRBadLen = 0;
		pEle = pSeqMtf->ppSeq[nIndex]->pMatElement;
		
		for(nj=0; nj<nOrder; nj++)
		{
			nFLabel[nj] = (int)(pEle[nFromPos-nOrder+nj]);
			nRLabel[nj] = (int)(pEle[nFromPos+nj]);

			if(nFLabel[nj] < nWidth)
			{
				nFWordId = nWidth*nFWordId+nFLabel[nj];
			}
			else
			{
				nFWordId = nWidth*nFWordId;
				nFBadLen++;
			}
			
			if(nRLabel[nj] < nWidth)
			{
				nRWordId = nWidth*nRWordId+nRLabel[nj];
			}
			else
			{
				nRWordId = nWidth*nRWordId;
				nRBadLen++;
			}
		}

		for(nj=nFromPos; nj<=nToPos; nj++)
		{
			if(nRLabel[0] < nWidth)
			{
				nRWordId -= nRLabel[0]*nScale;
			}
			else
			{
				nRBadLen--;
			}
			for(nz=0; nz<(nOrder-1); nz++)
			{
				nRLabel[nz] = nRLabel[nz+1];
			}
			nk = (int)(pEle[nj+nOrder]);
			nRLabel[nz] = nk;
			if(nk<nWidth)
			{
				nRWordId = nRWordId*nWidth+nk;
			}
			else
			{
				nRWordId = nRWordId*nWidth;
				nRBadLen++;
			}

			nk = (int)(pEle[nj]);
			if(nk < nWidth)
			{
				if(nFBadLen == 0)
				{
					dTemp = DMGETAT(pBGF, nFWordId, nk)+1.0;
					DMSETAT(pBGF, nFWordId, nk, dTemp);
				}
				if(nRBadLen == 0)
				{
					dTemp = DMGETAT(pBGR, nRWordId, nk)+1.0;
					DMSETAT(pBGR, nRWordId, nk, dTemp);
				}
			}

			if(nFLabel[0] < nWidth)
			{
				nFWordId -= nFLabel[0]*nScale;
			}
			else
			{
				nFBadLen--;
			}
			for(nz=0; nz<(nOrder-1); nz++)
			{
				nFLabel[nz] = nFLabel[nz+1];
			}
			nFLabel[nz] = nk;
			if(nk<nWidth)
			{
				nFWordId = nFWordId*nWidth+nk;
			}
			else
			{
				nFWordId = nFWordId*nWidth;
				nFBadLen++;
			}			
		}
	}
	
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/* SEQMTF_SETBGMCSCORE_BOTHSTRAND: get background log-likelihood for       */
/* every base in the sequences.                                            */
/* ----------------------------------------------------------------------- */ 
int SEQMTF_SETBGMCSCORE_BOTHSTRAND(struct SEQMOTIF *pSeqMtf, int nIndex, int nEffectiveLen, 
								   int nFromPos, int nToPos, 
								   struct DOUBLEMATRIX *pBGF, struct DOUBLEMATRIX *pBGR, 
								   int nOrder)
{
	/* define */
	int nHeight,nWidth,nScale;
	int nFLabel[BGMC_MAX_ORDER];
	int nRLabel[BGMC_MAX_ORDER];

	int nz,nj,nk,nLen,nFBadLen,nRBadLen;
	unsigned char *pEle;
	double *pBGFEle;
	double *pBGREle;
	int nFWordId,nRWordId,nWordId;
	double dTemp;

	/* check */
	if((pBGF == NULL) || (pBGR == NULL))
	{
		printf("Error: SEQMTF_SETBGMCSCORE_BOTHSTRAND, no background matrix!\n");
		exit(EXIT_FAILURE);
	}
	if(pSeqMtf == NULL)
	{
		printf("Error: SEQMTF_SETBGMCSCORE_BOTHSTRAND, no sequence!\n");
		exit(EXIT_FAILURE);
	}
	if((nIndex < 0) || (nIndex >= pSeqMtf->nSeqNum))
	{
		printf("Error: SEQMTF_SETBGMCSCORE_BOTHSTRAND, index out of range!\n");
		exit(EXIT_FAILURE);
	}
	if(pSeqMtf->ppSeq[nIndex] == NULL)
	{
		printf("Error: SEQMTF_SETBGMCSCORE_BOTHSTRAND, no sequence!\n");
		exit(EXIT_FAILURE);
	}
	if(pSeqMtf->nScoreNum < 2)
	{
		printf("Error: SEQMTF_SETBGMCSCORE_BOTHSTRAND, no place for storing scores!\n");
		exit(EXIT_FAILURE);
	}
	if((pSeqMtf->ppScore[0] == NULL) || (pSeqMtf->ppScore[1] == NULL))
	{
		printf("Error: SEQMTF_SETBGMCSCORE_BOTHSTRAND, no place for storing scores!\n");
		exit(EXIT_FAILURE);
	}
	if( (nEffectiveLen > pSeqMtf->ppSeq[nIndex]->nWidth) ||
		(nFromPos<nOrder) || (nFromPos>=nEffectiveLen) ||
		(nToPos<0) || (nToPos>=(nEffectiveLen-nOrder)) || (nFromPos>nToPos))
	{
		printf("Error: SEQMTF_SETBGMCSCORE_BOTHSTRAND, sequence out of range!\n");
		exit(EXIT_FAILURE);
	}
	if(pBGF->nWidth != pBGR->nWidth)
	{
		printf("Error: SEQMTF_SETBGMCSCORE_BOTHSTRAND, matrix dimension not match!\n");
		exit(EXIT_FAILURE);
	}

	
	/* init */
	nWidth = pBGF->nWidth;
	nHeight = (int)(pow((double)nWidth, (double)nOrder));
	if( (nHeight != pBGF->nHeight) || (nHeight != pBGR->nHeight) ) 
	{
		printf("Error: SEQMTF_SETBGMCSCORE_BOTHSTRAND, matrix dimension not match!\n");
		exit(EXIT_FAILURE);
	}

	nScale = (int)(pow((double)nWidth, (double)(nOrder-1)));

	/* create */
	for(nj=0; nj<BGMC_MAX_ORDER; nj++)
	{
		nFLabel[nj] = 0;
		nRLabel[nj] = 0;
	}

	nLen = pSeqMtf->ppSeq[nIndex]->nWidth;
	
	if(nOrder == 0)
	{
		pEle = pSeqMtf->ppSeq[nIndex]->pMatElement;
		pBGFEle = pSeqMtf->ppScore[0]->pMatElement;
		pBGREle = pSeqMtf->ppScore[1]->pMatElement;

		for(nj=nFromPos; nj<=nToPos; nj++)
		{
			nWordId = (int)(pEle[nj]);
			if(nWordId < nWidth)
			{
				dTemp = DMGETAT(pBGF, 0, nWordId);
				pBGFEle[nj] = dTemp;
				dTemp = DMGETAT(pBGR, 0, nWordId);
				pBGREle[nj] = dTemp;
			}
			else
			{
				pBGFEle[nj] = 1.0;
				pBGREle[nj] = 1.0;
			}
		}
	}
	else
	{
		nFWordId = 0;
		nFBadLen = 0;
		nRWordId = 0;
		nRBadLen = 0;
		pEle = pSeqMtf->ppSeq[nIndex]->pMatElement;
		pBGFEle = pSeqMtf->ppScore[0]->pMatElement;
		pBGREle = pSeqMtf->ppScore[1]->pMatElement;
		
		for(nj=0; nj<nOrder; nj++)
		{
			nFLabel[nj] = (int)(pEle[nFromPos-nOrder+nj]);
			nRLabel[nj] = (int)(pEle[nFromPos+nj]);

			if(nFLabel[nj] < nWidth)
			{
				nFWordId = nWidth*nFWordId+nFLabel[nj];
			}
			else
			{
				nFWordId = nWidth*nFWordId;
				nFBadLen++;
			}
			
			if(nRLabel[nj] < nWidth)
			{
				nRWordId = nWidth*nRWordId+nRLabel[nj];
			}
			else
			{
				nRWordId = nWidth*nRWordId;
				nRBadLen++;
			}
		}

		for(nj=nFromPos; nj<=nToPos; nj++)
		{
			if(nRLabel[0] < nWidth)
			{
				nRWordId -= nRLabel[0]*nScale;
			}
			else
			{
				nRBadLen--;
			}
			for(nz=0; nz<(nOrder-1); nz++)
			{
				nRLabel[nz] = nRLabel[nz+1];
			}
			nk = (int)(pEle[nj+nOrder]);
			nRLabel[nz] = nk;
			if(nk<nWidth)
			{
				nRWordId = nRWordId*nWidth+nk;
			}
			else
			{
				nRWordId = nRWordId*nWidth;
				nRBadLen++;
			}

			nk = (int)(pEle[nj]);
			if(nk < nWidth)
			{
				if(nFBadLen == 0)
				{
					dTemp = DMGETAT(pBGF, nFWordId, nk);
					pBGFEle[nj] = dTemp;
				}
				else
				{
					pBGFEle[nj] = 1.0;
				}

				if(nRBadLen == 0)
				{
					dTemp = DMGETAT(pBGR, nRWordId, nk);
					pBGREle[nj] = dTemp;
				}
				else
				{
					pBGREle[nj] = 1.0;
				}
			}
			else
			{
				pBGFEle[nj] = 1.0;
				pBGREle[nj] = 1.0;
			}

			if(nFLabel[0] < nWidth)
			{
				nFWordId -= nFLabel[0]*nScale;
			}
			else
			{
				nFBadLen--;
			}
			for(nz=0; nz<(nOrder-1); nz++)
			{
				nFLabel[nz] = nFLabel[nz+1];
			}
			nFLabel[nz] = nk;
			if(nk<nWidth)
			{
				nFWordId = nFWordId*nWidth+nk;
			}
			else
			{
				nFWordId = nFWordId*nWidth;
				nFBadLen++;
			}			
		}
	}
	
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/* SEQMTF_MAPMOTIF_BOTHSTRANDASBG: map a motif to a genomic sequence.      */
/* the map locations will be recorded in motif site list.                  */
/* ----------------------------------------------------------------------- */ 
int SEQMTF_MAPMOTIF_BOTHSTRANDASBG(struct DOUBLEMATRIX *pMotifPWM, int nMotifId, double dCutoff, 
								   struct SEQMOTIF *pSeqMtf, int nIndex, int nEffectiveLen, 
								   int nFromPos, int nToPos)
{
	/* define */
	unsigned char *pBase;
	unsigned char bBase;
	double *pBGFEle,*pBGREle;
	int ni,nj,nk;
	int nMotifLen;
	int nBaseBadLen, nBGFBadLen, nBGRBadLen;
	int nWidth;

	double dRatio;
	double dBGF,dBGR,dMF,dMR;
	double dML,dBGL;

	/* pointer for adding motif sites */
	struct MOTIFSITE *pLastSite;
	struct MOTIFSITE *pNewSite;
	int nStrand;

	/* check parameters */
	if(pMotifPWM == NULL)
	{
		printf("Warning: SEQMTF_MAPMOTIF_BOTHSTRANDASBG, no PWM for motif mapping!\n");
		return PROC_FAILURE;
	}
	
	nMotifLen = pMotifPWM->nHeight;

	if(pSeqMtf == NULL)
	{
		printf("Error: SEQMTF_MAPMOTIF_BOTHSTRANDASBG, no sequence!\n");
		exit(EXIT_FAILURE);
	}
	if((nIndex < 0) || (nIndex >= pSeqMtf->nSeqNum))
	{
		printf("Error: SEQMTF_MAPMOTIF_BOTHSTRANDASBG, index out of range!\n");
		exit(EXIT_FAILURE);
	}
	if(pSeqMtf->ppSeq[nIndex] == NULL)
	{
		printf("Error: SEQMTF_MAPMOTIF_BOTHSTRANDASBG, no sequence!\n");
		exit(EXIT_FAILURE);
	}
	if(pSeqMtf->nScoreNum < 2)
	{
		printf("Error: SEQMTF_MAPMOTIF_BOTHSTRANDASBG, no place for storing scores!\n");
		exit(EXIT_FAILURE);
	}
	if((pSeqMtf->ppScore[0] == NULL) || (pSeqMtf->ppScore[1] == NULL))
	{
		printf("Error: SEQMTF_MAPMOTIF_BOTHSTRANDASBG, no place for storing scores!\n");
		exit(EXIT_FAILURE);
	}
	if( (nEffectiveLen > pSeqMtf->ppSeq[nIndex]->nWidth) ||
		(nFromPos<0) || (nFromPos>(nEffectiveLen-nMotifLen)) ||
		(nToPos<0) || (nToPos>(nEffectiveLen-nMotifLen)) || (nFromPos>nToPos))
	{
		printf("Error: SEQMTF_MAPMOTIF_BOTHSTRANDASBG, sequence out of range!\n");
		exit(EXIT_FAILURE);
	}

	/* init */
	nWidth = pMotifPWM->nWidth;
	
	if(pSeqMtf->pMotifSite == NULL)
	{
		if(pSeqMtf->nSiteNum != 0)
		{
			printf("Error: SEQMTF_MAPMOTIF_BOTHSTRANDASBG, nMotifSite number and recorded motif sites not match\n");
			exit(EXIT_FAILURE);
		}
		pLastSite = NULL;
	}
	else
	{
		pLastSite = pSeqMtf->pMotifSite;
		ni = 1;
		while(pLastSite->pNext != NULL)
		{
			pLastSite = pLastSite->pNext;
			ni++;
		}
		if(ni != pSeqMtf->nSiteNum)
		{
			printf("Error: SEQMTF_MAPMOTIF_BOTHSTRANDASBG, nMotifSite number and recorded motif sites not match\n");
			exit(EXIT_FAILURE);
		}
	}

		
	/* get motif score*/
	pBase = pSeqMtf->ppSeq[nIndex]->pMatElement;
	pBGFEle = pSeqMtf->ppScore[0]->pMatElement;
	pBGREle = pSeqMtf->ppScore[1]->pMatElement;
	
	nBaseBadLen = 0;
	nBGFBadLen = 0;
	nBGRBadLen = 0;
	dBGF = 0.0;
	dBGR = 0.0;
	dMF = 0.0;
	dMR = 0.0;
	
	/* prepare background */
	for(ni=nFromPos; ni<(nFromPos+nMotifLen); ni++)
	{
		bBase = (int)(pBase[ni]);
		
		if(bBase >= nWidth)
		{
			nBaseBadLen++;
		}

		if(pBGFEle[ni] >= 0.0)
		{
			nBGFBadLen++;
		}
		else
		{
			dBGF += pBGFEle[ni];
		}

		if(pBGREle[ni] >= 0.0)
		{
			nBGRBadLen++;
		}
		else
		{
			dBGR += pBGREle[ni];				
		}
	}

	/* first base */
	if( (nBaseBadLen == 0) && (nBGFBadLen == 0) && (nBGRBadLen == 0) )
	{
		dMF = 0.0;
		dMR = 0.0;
		for(nj=0; nj<nMotifLen; nj++)
		{
			bBase = (int)(pBase[nFromPos+nj]);
			dMF += DMGETAT(pMotifPWM, nj, bBase);
			dMR += DMGETAT(pMotifPWM, (nMotifLen-1-nj), (3-bBase));
		}

		if(dMF>dMR)
		{
			dML = dMF;
			nStrand = 0;
		}
		else
		{
			dML = dMR;
			nStrand = 1;
		}

		if(dBGF>dBGR)
			dBGL = dBGF;
		else
			dBGL = dBGR;

		/* likelihood ratio */
		dRatio = exp(dML-dBGL);
		if(dRatio >= dCutoff)
		{
			/* add new motif site */
			pNewSite = NULL;
			pNewSite = MTFSCREATE();
			if(pNewSite == NULL)
			{
				printf("Error: SEQMTF_MAPMOTIF_BOTHSTRANDASBG, cannot create motif site!\n");
				exit(EXIT_FAILURE);
			}
			pNewSite->dScore = dRatio;
			pNewSite->nMotifLen = nMotifLen;
			pNewSite->nMotifType = nMotifId;
			pNewSite->nSeqId = nIndex;
			pNewSite->nStartPos = nFromPos;
			pNewSite->nStrand = nStrand;
			pNewSite->pNext = NULL;

			if(pLastSite == NULL)
			{
				pSeqMtf->pMotifSite = pNewSite;
			}
			else
			{
				pLastSite->pNext = pNewSite;
			}
			pLastSite = pNewSite;
			pSeqMtf->nSiteNum = pSeqMtf->nSiteNum+1;
		}
	}


	/* motif map */
	for(ni=(nFromPos+1); ni<=nToPos; ni++)
	{
		/* remove the oldest base */
		nk = ni-1;
		bBase = (int)(pBase[nk]);
		
		if(bBase >= nWidth)
		{
			nBaseBadLen--;
		}

		if(pBGFEle[nk] >= 0.0)
		{
			nBGFBadLen--;
		}
		else
		{
			dBGF -= pBGFEle[nk];
		}

		if(pBGREle[nk] >= 0.0)
		{
			nBGRBadLen--;
		}
		else
		{
			dBGR -= pBGREle[nk];				
		}

		/* add the newest base */
		nk = ni+nMotifLen-1;
		bBase = (int)(pBase[nk]);
		
		if(bBase >= nWidth)
		{
			nBaseBadLen++;
		}

		if(pBGFEle[nk] >= 0.0)
		{
			nBGFBadLen++;
		}
		else
		{
			dBGF += pBGFEle[nk];
		}

		if(pBGREle[nk] >= 0.0)
		{
			nBGRBadLen++;
		}
		else
		{
			dBGR += pBGREle[nk];				
		}

		/* calculate motif score */
		if( (nBaseBadLen == 0) && (nBGFBadLen == 0) && (nBGRBadLen == 0) )
		{
			dMF = 0.0;
			dMR = 0.0;
			for(nj=0; nj<nMotifLen; nj++)
			{
				bBase = (int)(pBase[ni+nj]);
				dMF += DMGETAT(pMotifPWM, nj, bBase);
				dMR += DMGETAT(pMotifPWM, (nMotifLen-1-nj), (3-bBase));
			}

			if(dMF>dMR)
			{
				dML = dMF;
				nStrand = 0;
			}
			else
			{
				dML = dMR;
				nStrand = 1;
			}

			if(dBGF>dBGR)
				dBGL = dBGF;
			else
				dBGL = dBGR;

			/* likelihood ratio */
			dRatio = exp(dML-dBGL);
			if(dRatio >= dCutoff)
			{
				/* add new motif site */
				pNewSite = NULL;
				pNewSite = MTFSCREATE();
				if(pNewSite == NULL)
				{
					printf("Error: SEQMTF_MAPMOTIF_BOTHSTRANDASBG, cannot create motif site!\n");
					exit(EXIT_FAILURE);
				}
				pNewSite->dScore = dRatio;
				pNewSite->nMotifLen = nMotifLen;
				pNewSite->nMotifType = nMotifId;
				pNewSite->nSeqId = nIndex;
				pNewSite->nStartPos = ni;
				pNewSite->nStrand = nStrand;
				pNewSite->pNext = NULL;

				if(pLastSite == NULL)
				{
					pSeqMtf->pMotifSite = pNewSite;
				}
				else
				{
					pLastSite->pNext = pNewSite;
				}
				pLastSite = pNewSite;
				pSeqMtf->nSiteNum = pSeqMtf->nSiteNum+1;
			}
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/* SEQMTF_WRITEMOTIFSITETOBED_THENCLEARSITES: write motifs to bed files.   */
/* ----------------------------------------------------------------------- */ 
int SEQMTF_WRITEMOTIFSITETOBED_THENCLEARSITES(struct SEQMOTIF *pSeqMtf, int nMaxMotifNum, 
											  struct tagString **vMotifName, 
											  char strChrName[], int nOffset,
											  FILE **vfpMapOut)
{
	/* define */
	struct MOTIFSITE *pSite;
	int ni;
	char chStrand;
	
	/* check parameters */
	if(pSeqMtf == NULL)
	{
		printf("Warning: SEQMTF_WRITEMOTIFSITETOBED_THENCLEARSITES, no sequence!\n");
		return PROC_FAILURE;
	}
	if((vMotifName == NULL) || (vfpMapOut == NULL)) 
	{
		printf("Warning: SEQMTF_WRITEMOTIFSITETOBED_THENCLEARSITES, lack output information!\n");
		return PROC_FAILURE;
	}

	/* save */
	while(pSeqMtf->pMotifSite != NULL)
	{
		pSite = pSeqMtf->pMotifSite;
		pSeqMtf->pMotifSite = pSite->pNext;
		pSeqMtf->nSiteNum = pSeqMtf->nSiteNum-1;

		ni = pSite->nMotifType;
		if(ni >= nMaxMotifNum)
		{
			printf("Warning: SEQMTF_WRITEMOTIFSITETOBED_THENCLEARSITES, motif index out of range!\n");
			MTFSDESTROY(pSite);
			continue;
		}

		if(vfpMapOut[ni] == NULL || vMotifName[ni] == NULL)
		{
			printf("Warning: SEQMTF_WRITEMOTIFSITETOBED_THENCLEARSITES, lack output information!\n");
			MTFSDESTROY(pSite);
		}
		else
		{
			if(pSite->nStrand == 1)
			{
				chStrand = '-';
			}
			else
			{
				chStrand = '+';
			}
			fprintf(vfpMapOut[ni], "%s\t%d\t%d\t%s\t%f\t%c\n", strChrName, 
				(nOffset+pSite->nStartPos), (nOffset+pSite->nStartPos+pSite->nMotifLen-1),
				vMotifName[ni]->m_pString, pSite->dScore, chStrand);
			MTFSDESTROY(pSite);
		}
	}

	if(pSeqMtf->nSiteNum != 0)
	{
		printf("Error: SEQMTF_WRITEMOTIFSITETOBED_THENCLEARSITES, motif number not correct!\n");
		exit(EXIT_FAILURE);
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/* SEQMTF_LOADSEQCSFROMGENOME_FORMOTIFFILTER: load sequence and            */
/* conservation score.                                                     */
/* ----------------------------------------------------------------------- */ 
int SEQMTF_LOADSEQCSFROMGENOME_FORMOTIFFILTER(struct SEQMOTIF *pSeqMtf, int nSeqId, 
											  char strChrName[], int nStart, int nEnd, 
											  char strGenomePath[], char strCScorePath[])
{
	/* define */
	char strLine[LINE_LENGTH];
	FILE *fpSeq;
	FILE *fpCS;

	unsigned char *pBase;
	unsigned char *pRepeat;
	unsigned char *pCScore;
	int nLen;
	int nP1,nP2,nR1,nR2;
	int numread;
	int ni,nk;
	unsigned char bBase,bChar;

	/* init */
	if(pSeqMtf == NULL)
	{
		printf("Error: SEQMTF_LOADSEQCSFROMGENOME_FORMOTIFFILTER, null seqmotif!\n");
		exit(EXIT_FAILURE);
	}
	if( (nSeqId<0) || (nSeqId >= pSeqMtf->nSeqNum) )
	{
		printf("Error: SEQMTF_LOADSEQCSFROMGENOME_FORMOTIFFILTER, sequence id out of range!\n");
		exit(EXIT_FAILURE);
	}
	nLen = (nEnd-nStart+1);
	if( (nLen<=0) || (nLen > pSeqMtf->ppSeq[nSeqId]->nWidth ) )
	{
		printf("Error: SEQMTF_LOADSEQCSFROMGENOME_FORMOTIFFILTER, the length of the sequence to be retrieved is larger than the length allowed by seqmotif!\n");
		exit(EXIT_FAILURE);
	}
	if(pSeqMtf->nPathNum < 2)
	{
		printf("Error: SEQMTF_LOADSEQCSFROMGENOME_FORMOTIFFILTER, no space for conservation score and repeatmasker!\n");
		exit(EXIT_FAILURE);
	}


	/* get sequences */
	sprintf(strLine, "%s%s.sq", strGenomePath, strChrName);
	fpSeq = NULL;
	fpSeq = fopen(strLine, "rb");
	if(fpSeq == NULL)
	{
		printf("Error: SEQMTF_LOADSEQCSFROMGENOME_FORMOTIFFILTER, cannot open *.sq file!\n");
		exit(EXIT_FAILURE);
	}

	/* get seq */
	pBase = pSeqMtf->ppSeq[nSeqId]->pMatElement;
	pRepeat = pSeqMtf->ppSamplePath[1]->pMatElement;
	nP1 = (int)(nStart/2);
	nR1 = nStart%2;
	nP2 = (int)(nEnd/2);
	nR2 = nEnd%2;
	nk = 0;

	if( fseek( fpSeq, nP1, SEEK_SET ) != 0 )
	{
		printf("Error: SEQMTF_LOADSEQCSFROMGENOME_FORMOTIFFILTER, cannot locate the required sequence!\n");
		exit(EXIT_FAILURE);
	}


	/* first base */
	numread = fread( &bChar, sizeof( unsigned char ), 1, fpSeq );
	if(nR1 == 0)
	{
		bBase = bChar >> 4;
		switch(bBase)
		{
			case 0: *pBase = 0;
				*pRepeat = 0;
				break;
			case 1: *pBase = 1;
				*pRepeat = 0;
				break;
			case 2: *pBase = 2;
				*pRepeat = 0;
				break;
			case 3: *pBase = 3;
				*pRepeat = 0;
				break;
			case 4: *pBase = 0;
				*pRepeat = 1;
				break;
			case 5: *pBase = 1;
				*pRepeat = 1;
				break;
			case 6: *pBase = 2;
				*pRepeat = 1;
				break;
			case 7: *pBase = 3;
				*pRepeat = 1;
				break;
			default: *pBase = 4;
				*pRepeat = 1;
		}
		pBase++;
		pRepeat++;
		nk++;
	}

	if(nk < nLen)
	{
		bBase = bChar & 0x0F;
		switch(bBase)
		{
			case 0: *pBase = 0;
				*pRepeat = 0;
				break;
			case 1: *pBase = 1;
				*pRepeat = 0;
				break;
			case 2: *pBase = 2;
				*pRepeat = 0;
				break;
			case 3: *pBase = 3;
				*pRepeat = 0;
				break;
			case 4: *pBase = 0;
				*pRepeat = 1;
				break;
			case 5: *pBase = 1;
				*pRepeat = 1;
				break;
			case 6: *pBase = 2;
				*pRepeat = 1;
				break;
			case 7: *pBase = 3;
				*pRepeat = 1;
				break;
			default: *pBase = 4;
				*pRepeat = 1;
		}
		pBase++;
		pRepeat++;
		nk++;
	}

	/* middle bases */
	if(nk < nLen)
	{
		for(ni=(nP1+1); ni<nP2; ni++)
		{
			numread = fread( &bChar, sizeof( unsigned char ), 1, fpSeq );
			bBase = bChar >> 4;
			switch(bBase)
			{
				case 0: *pBase = 0;
					*pRepeat = 0;
					break;
				case 1: *pBase = 1;
					*pRepeat = 0;
					break;
				case 2: *pBase = 2;
					*pRepeat = 0;
					break;
				case 3: *pBase = 3;
					*pRepeat = 0;
					break;
				case 4: *pBase = 0;
					*pRepeat = 1;
					break;
				case 5: *pBase = 1;
					*pRepeat = 1;
					break;
				case 6: *pBase = 2;
					*pRepeat = 1;
					break;
				case 7: *pBase = 3;
					*pRepeat = 1;
					break;
				default: *pBase = 4;
					*pRepeat = 1;
			}
			pBase++;
			pRepeat++;
			nk++;

			bBase = bChar & 0x0F;
			switch(bBase)
			{
				case 0: *pBase = 0;
					*pRepeat = 0;
					break;
				case 1: *pBase = 1;
					*pRepeat = 0;
					break;
				case 2: *pBase = 2;
					*pRepeat = 0;
					break;
				case 3: *pBase = 3;
					*pRepeat = 0;
					break;
				case 4: *pBase = 0;
					*pRepeat = 1;
					break;
				case 5: *pBase = 1;
					*pRepeat = 1;
					break;
				case 6: *pBase = 2;
					*pRepeat = 1;
					break;
				case 7: *pBase = 3;
					*pRepeat = 1;
					break;
				default: *pBase = 4;
					*pRepeat = 1;
			}
			pBase++;
			pRepeat++;
			nk++;
		}
	}
	
	/* last base */
	if(nk < nLen)
	{
		numread = fread( &bChar, sizeof( unsigned char ), 1, fpSeq );
		bBase = bChar >> 4;
		switch(bBase)
		{
			case 0: *pBase = 0;
				*pRepeat = 0;
				break;
			case 1: *pBase = 1;
				*pRepeat = 0;
				break;
			case 2: *pBase = 2;
				*pRepeat = 0;
				break;
			case 3: *pBase = 3;
				*pRepeat = 0;
				break;
			case 4: *pBase = 0;
				*pRepeat = 1;
				break;
			case 5: *pBase = 1;
				*pRepeat = 1;
				break;
			case 6: *pBase = 2;
				*pRepeat = 1;
				break;
			case 7: *pBase = 3;
				*pRepeat = 1;
				break;
			default: *pBase = 4;
				*pRepeat = 1;
		}
		pBase++;
		pRepeat++;
		nk++;
	}

	if(nk < nLen)
	{
		if(nR2 == 1)
		{
			bBase = bChar & 0x0F;
			switch(bBase)
			{
				case 0: *pBase = 0;
					*pRepeat = 0;
					break;
				case 1: *pBase = 1;
					*pRepeat = 0;
					break;
				case 2: *pBase = 2;
					*pRepeat = 0;
					break;
				case 3: *pBase = 3;
					*pRepeat = 0;
					break;
				case 4: *pBase = 0;
					*pRepeat = 1;
					break;
				case 5: *pBase = 1;
					*pRepeat = 1;
					break;
				case 6: *pBase = 2;
					*pRepeat = 1;
					break;
				case 7: *pBase = 3;
					*pRepeat = 1;
					break;
				default: *pBase = 4;
					*pRepeat = 1;
			}
			pBase++;
			pRepeat++;
			nk++;
		}
	}

	if(nk != nLen)
	{
		printf("Error: SEQMTF_LOADSEQCSFROMGENOME_FORMOTIFFILTER, sequence length not match!\n");
		exit(EXIT_FAILURE);
	}

	/* close seq file */
	fclose(fpSeq);


	/* load conservation score */
	sprintf(strLine, "%s%s.cs", strCScorePath, strChrName);
	fpCS = NULL;
	fpCS = fopen(strLine, "rb");
	if(fpCS == NULL)
	{
		printf("Error: SEQMTF_LOADSEQCSFROMGENOME_FORMOTIFFILTER, cannot open *.cs file!\n");
		exit(EXIT_FAILURE);
	}

	/* load */
	pCScore = pSeqMtf->ppSamplePath[0]->pMatElement;

	if( fseek( fpCS, nStart, SEEK_SET ) != 0 )
	{
		printf("Error: SEQMTF_LOADSEQCSFROMGENOME_FORMOTIFFILTER, cannot locate the required sequence!\n");
		exit(EXIT_FAILURE);
	}

	/* read */
	numread = fread(pCScore, sizeof(unsigned char), nLen, fpCS);
	if(numread != nLen)
	{
		printf("Error: SEQMTF_LOADSEQCSFROMGENOME_FORMOTIFFILTER, loading error!\n");
		exit(EXIT_FAILURE);
	}

	/* close cs file */
	fclose(fpCS);
	
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*             MCMOVECREATE: create markov chain move elements             */
/* ----------------------------------------------------------------------- */ 
struct MARKOVCHAINMOVE *MCMOVECREATE()
{
	/* new complex */
	struct MARKOVCHAINMOVE *pM;

	/* create */
	pM = NULL;
	pM = (struct MARKOVCHAINMOVE *)calloc(1, sizeof(struct MARKOVCHAINMOVE));
	if (pM == NULL)
	{
		printf("Error: cannot create MCMOVE complex.\n");
		exit(EXIT_FAILURE);
	}
	
	/* init */
	
	/* return */
	return pM;
}

/* ----------------------------------------------------------------------- */ 
/*             MCMOVEDESTROY: destroy markov chain move elements           */
/* ----------------------------------------------------------------------- */ 
int MCMOVEDESTROY(struct MARKOVCHAINMOVE *pM)
{
	if(pM == NULL)
		return PROC_SUCCESS;

	/* destroy */
	pM->pNext = NULL;
	
	free(pM);
	pM = NULL;

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*             MCMOVEDESTROYCIRCLE: destroy markov chain circle            */
/* ----------------------------------------------------------------------- */ 
int MCMOVEDESTROYCIRCLE(struct MARKOVCHAINMOVE **pMList)
{
	struct MARKOVCHAINMOVE *pM;

	if(pMList==NULL)
		return PROC_SUCCESS;

	/* destroy */
	pM = *pMList;
	*pMList = pM->pNext;
	if(pM == *pMList)
	{
		pM->pNext = NULL;
		MCMOVEDESTROY(pM);
		*pMList = NULL;
	}
	else
	{
		pM->pNext = NULL;
		while(*pMList != NULL)
		{
			pM = *pMList;
			*pMList = pM->pNext;
			MCMOVEDESTROY(pM);
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*             SITEPOSCREATE: create site pos elements                     */
/* ----------------------------------------------------------------------- */ 
struct SITEPOS *SITEPOSCREATE()
{
	/* new complex */
	struct SITEPOS *pS;

	/* create */
	pS = NULL;
	pS = (struct SITEPOS *)calloc(1, sizeof(struct SITEPOS));
	if (pS == NULL)
	{
		printf("Error: cannot create SITEPOS element.\n");
		exit(EXIT_FAILURE);
	}
	
	/* init */
	
	/* return */
	return pS;
}

/* ----------------------------------------------------------------------- */ 
/*             SITEPOSDESTROY: destroy site pos elements                   */
/* ----------------------------------------------------------------------- */ 
int SITEPOSDESTROY(struct SITEPOS *pS)
{
	if(pS == NULL)
		return PROC_SUCCESS;

	/* destroy */
	pS->pNext = NULL;
	
	free(pS);
	pS = NULL;

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*             SITEPOSDESTROYLIST: destroy position list                   */
/* ----------------------------------------------------------------------- */ 
int SITEPOSDESTROYLIST(struct SITEPOS **pSList)
{
	/* new motif site */
	struct SITEPOS *pS;

	if(pSList==NULL)
		return PROC_SUCCESS;

	/* destroy */
	while(*pSList != NULL)
	{
		pS = *pSList;
		*pSList = pS->pNext;
		SITEPOSDESTROY(pS);
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*             SITEPOSLISTGETCOUNT: get # of site pos elements             */
/* ----------------------------------------------------------------------- */ 
int SITEPOSLISTGETCOUNT(struct SITEPOS *pSList)
{
	/* define */
	int nCount;
	struct SITEPOS *pS;

	/* count */
	if(pSList == NULL)
		return 0;
	else
	{
		nCount = 0;
		pS = pSList;
		while(pS != NULL)
		{
			nCount++;
			pS = pS->pNext;
		}
	}

	/* return */
	return nCount;
}

/* ----------------------------------------------------------------------- */ 
/*              MTFSGRPCREATE: create motif site group                     */
/* ----------------------------------------------------------------------- */ 
struct MOTIFSITEGROUP *MTFSGRPCREATE(int nMotifNum)
{
	/* new motif site */
	struct MOTIFSITEGROUP *pG;

	/* create motif site */
	if(nMotifNum <= 0)
		return NULL;

	pG = NULL;
	pG = (struct MOTIFSITEGROUP *)calloc(1, sizeof(struct MOTIFSITEGROUP));
	if (pG == NULL)
	{
		printf("Error: cannot create motif site group!\n");
		exit(EXIT_FAILURE);
	}

	pG->nMotifNum = nMotifNum;
	pG->pSites = (struct MOTIFSITE **)calloc(nMotifNum, sizeof(struct MOTIFSITE *));
	if(pG->pSites == NULL)
	{
		printf("Error: cannot create motif site group!\n");
		exit(EXIT_FAILURE);
	}

	/* return */
	return pG;
}

/* ----------------------------------------------------------------------- */ 
/*             MTFSGRPDESTROY: destroy motif site group                    */
/* ----------------------------------------------------------------------- */ 
int MTFSGRPDESTROY(struct MOTIFSITEGROUP *pG)
{
	/* define */
	int ni;

	/* init check */
	if(pG == NULL)
		return PROC_SUCCESS;

	/* destroy */
	for(ni=0; ni<pG->nMotifNum; ni++)
	{
		MTFSDESTROYLIST(pG->pSites+ni);
		pG->pSites[ni] = NULL;
	}
	if(pG->pSites != NULL)
	{
		free(pG->pSites);
	}
	
	/* free all */
	free(pG);
	pG = NULL;

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*             FLEXMTFSCREATE: create site of flexmotifs                   */
/* ----------------------------------------------------------------------- */ 
struct FLEXMOTIFSITE *FLEXMTFSCREATE()
{
	/* new motif site */
	struct FLEXMOTIFSITE *pS;

	/* create motif site */
	pS = NULL;
	pS = (struct FLEXMOTIFSITE*)calloc(1, sizeof(struct FLEXMOTIFSITE));
	if (pS == NULL)
	{
		printf("Error: cannot create flex motif site.\n");
		exit(EXIT_FAILURE);
	}

	pS->pContextPos = NULL;
	pS->pContextMotif = NULL;
	pS->pContextStrand = NULL;
	pS->pPrev = NULL;
	pS->pNext = NULL;

	/* return */
	return pS;
}

/* ----------------------------------------------------------------------- */ 
/*             FLEXMTFSDESTROY: destroy flexmotif site                     */
/* ----------------------------------------------------------------------- */ 
int FLEXMTFSDESTROY(struct FLEXMOTIFSITE *pS)
{
	if(pS == NULL)
		return PROC_SUCCESS;

	/* destroy */
	DestroyIntMatrix(pS->pContextPos);
	DestroyIntMatrix(pS->pContextMotif);
	DestroyIntMatrix(pS->pContextStrand);
	pS->pPrev = NULL;
	pS->pNext = NULL;
	
	free(pS);
	pS = NULL;

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*             FLEXMTFSDESTROYLIST: destroy list of flex motif sites       */
/* ----------------------------------------------------------------------- */ 
int FLEXMTFSDESTROYLIST(struct FLEXMOTIFSITE **pSList)
{
	/* new motif site */
	struct FLEXMOTIFSITE *pS;

	if(pSList==NULL)
		return PROC_SUCCESS;

	/* destroy */
	while(*pSList != NULL)
	{
		pS = *pSList;
		*pSList = pS->pNext;
		FLEXMTFSDESTROY(pS);
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/* FLEXMTFSADDTOLIST_SORTBYPROB: add flex motif sites to list, and sort    */
/* them according to the posterior probability.                            */
/* ----------------------------------------------------------------------- */ 
int FLEXMTFSADDTOLIST_SORTBYPROB(struct FLEXMOTIFSITE **pSList, struct FLEXMOTIFSITE *pSite)
{
	/* define */
	struct FLEXMOTIFSITE *pPrev,*pNext;

	/* check */
	if(pSList == NULL)
	{
		printf("Error: FLEXMTFSADDTOLIST_SORTBYPROB, null flex motif site list!\n");
		exit(EXIT_FAILURE);
	}

	/* add */
	if(*pSList == NULL)
	{
		*pSList = pSite;
	}
	else
	{
		pPrev = NULL;
		pNext = *pSList;
		while(pNext != NULL)
		{
			if(fabs(pNext->dProb-pSite->dProb) < 1e-6)
			{
				if( pNext->dScore < pSite->dScore )
				{
					break;
				}
				else
				{
					pPrev = pNext;
					pNext = pNext->pNext;
				}
			}
			else if(pNext->dProb < pSite->dProb)
			{
				break;
			}
			else
			{
				pPrev = pNext;
				pNext = pNext->pNext;
			}
		}

		if(pNext != NULL)
		{
			if(pPrev != NULL)
			{
				pPrev->pNext = pSite;
				pSite->pPrev = pPrev;
				pSite->pNext = pNext;
				pNext->pPrev = pSite;
			}
			else
			{
				pSite->pNext = pNext;
				pNext->pPrev = pSite;
				*pSList = pSite;
			}
		}
		else
		{
			if(pPrev != NULL)
			{
				pPrev->pNext = pSite;
				pSite->pPrev = pPrev;
			}
			else
			{
				printf("Error: FLEXMTFSADDTOLIST_SORTBYPROB, search list wrong!\n");
				exit(EXIT_FAILURE);
			}
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/* FLEXMTFSADDTOLIST_SORTBYPOS: add flex motif sites to list, and sort     */
/* them according to their position.                                       */
/* ----------------------------------------------------------------------- */ 
int FLEXMTFSADDTOLIST_SORTBYPOS(struct FLEXMOTIFSITE **pSList, struct FLEXMOTIFSITE *pSite)
{
	/* define */
	struct FLEXMOTIFSITE *pPrev,*pNext;

	/* check */
	if(pSList == NULL)
	{
		printf("Error: FLEXMTFSADDTOLIST_SORTBYPOS, null flex motif site list!\n");
		exit(EXIT_FAILURE);
	}

	/* add */
	if(*pSList == NULL)
	{
		*pSList = pSite;
	}
	else
	{
		pPrev = NULL;
		pNext = *pSList;
		while(pNext != NULL)
		{
			if(pNext->nSeqId > pSite->nSeqId)
			{
				break;
			}
			else if( (pNext->nSeqId == pSite->nSeqId) && (pNext->nStartPos > pSite->nStartPos) )
			{
				break;
			}
			else
			{
				pPrev = pNext;
				pNext = pNext->pNext;
			}
		}

		if(pNext != NULL)
		{
			if(pPrev != NULL)
			{
				pPrev->pNext = pSite;
				pSite->pPrev = pPrev;
				pSite->pNext = pNext;
				pNext->pPrev = pSite;
			}
			else
			{
				pSite->pNext = pNext;
				pNext->pPrev = pSite;
				*pSList = pSite;
			}
		}
		else
		{
			if(pPrev != NULL)
			{
				pPrev->pNext = pSite;
				pSite->pPrev = pPrev;
			}
			else
			{
				printf("Error: FLEXMTFSADDTOLIST_SORTBYPOS, search list wrong!\n");
				exit(EXIT_FAILURE);
			}
		}
	}
	
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*             FLEXMTFSREMOVESITE: remove one site from the matrix         */
/* ----------------------------------------------------------------------- */ 
int FLEXMTFSREMOVESITE(struct FLEXMOTIFSITE *pS, int nPos)
{
	/* define */
	int nNum,nNewNum,ni,nj;
	struct INTMATRIX *pPos;
	struct INTMATRIX *pMotif;
	struct INTMATRIX *pStrand;


	/* check */
	if(pS == NULL)
		return 0;
	if(pS->pContextPos == NULL)
	{
		printf("Warning: FLEXMTFSREMOVESITE, no site to remove!\n");
		return 0;
	}

	/* init */
	nNum = pS->pContextPos->nWidth;
	nNewNum = nNum-1;
	
	/* remove site */
	if(nNewNum == 0)
	{
		if(nPos != pS->pContextPos->pMatElement[0])
		{
			printf("Error: FLEXMTFSREMOVESITE, site coordinates wrong!\n");
			exit(EXIT_FAILURE);
		}
		DestroyIntMatrix(pS->pContextPos);
		pS->pContextPos = NULL;
		DestroyIntMatrix(pS->pContextMotif);
		pS->pContextMotif = NULL;
		DestroyIntMatrix(pS->pContextStrand);
		pS->pContextStrand = NULL;
	}
	else
	{
		pPos = NULL;
		pPos = CreateIntMatrix(1, nNewNum);
		if(pPos == NULL)
		{
			printf("Error: FLEXMTFSREMOVESITE, cannot reduce context sites!\n");
			exit(EXIT_FAILURE);
		}

		pMotif = NULL;
		pMotif = CreateIntMatrix(1, nNewNum);
		if(pMotif == NULL)
		{
			printf("Error: FLEXMTFSREMOVESITE, cannot reduce context sites!\n");
			exit(EXIT_FAILURE);
		}


		pStrand = NULL;
		pStrand = CreateIntMatrix(1, nNewNum);
		if(pStrand == NULL)
		{
			printf("Error: FLEXMTFSREMOVESITE, cannot reduce context sites!\n");
			exit(EXIT_FAILURE);
		}

		nj = 0;
		for(ni=0; ni<nNum; ni++)
		{
			if(pS->pContextPos->pMatElement[ni] != nPos)
			{
				pPos->pMatElement[nj] = pS->pContextPos->pMatElement[ni];
				pMotif->pMatElement[nj] = pS->pContextMotif->pMatElement[ni];
				pStrand->pMatElement[nj] = pS->pContextStrand->pMatElement[ni];
				nj++;
			}
		}

		DestroyIntMatrix(pS->pContextPos);
		pS->pContextPos = pPos;
		DestroyIntMatrix(pS->pContextMotif);
		pS->pContextMotif = pMotif;
		DestroyIntMatrix(pS->pContextStrand);
		pS->pContextStrand = pStrand;
	}

	/* return */
	return nNewNum;
}

/* ----------------------------------------------------------------------- */ 
/*            FLEXMTFSADDSITE: add one site to the context matrix.         */
/* ----------------------------------------------------------------------- */ 
int FLEXMTFSADDSITE(struct FLEXMOTIFSITE *pS, int nPos, int nMotifId, int nStrand)
{
	/* define */
	int nNum,nNewNum,ni,nj;
	struct INTMATRIX *pPos;
	struct INTMATRIX *pMotif;
	struct INTMATRIX *pStrand;
	int nAdd;

	/* check */
	if(pS == NULL)
	{
		printf("Error: FLEXMTFSADDSITE, null motif site structure!\n");
		exit(EXIT_FAILURE);
	}
	if(pS->pContextPos == NULL)
	{
		nNum = 0;
	}
	else
	{
		nNum = pS->pContextPos->nWidth;
	}

	/* init */
	nNewNum = nNum+1;

	pPos = NULL;
	pPos = CreateIntMatrix(1, nNewNum);
	if(pPos == NULL)
	{
		printf("Error: FLEXMTFSADDSITE, cannot add context sites!\n");
		exit(EXIT_FAILURE);
	}

	pMotif = NULL;
	pMotif = CreateIntMatrix(1, nNewNum);
	if(pMotif == NULL)
	{
		printf("Error: FLEXMTFSADDSITE, cannot add context sites!\n");
		exit(EXIT_FAILURE);
	}

	pStrand = NULL;
	pStrand = CreateIntMatrix(1, nNewNum);
	if(pStrand == NULL)
	{
		printf("Error: FLEXMTFSADDSITE, cannot add context sites!\n");
		exit(EXIT_FAILURE);
	}

	nAdd = 0;
	nj = 0;
	for(ni=0; ni<nNum; ni++)
	{
		if( (pS->pContextPos->pMatElement[ni] > nPos) && (nAdd == 0))
		{
			pPos->pMatElement[nj] = nPos;
			pMotif->pMatElement[nj] = nMotifId;
			pStrand->pMatElement[nj] = nStrand;
			nj++;
			nAdd = 1;
		}
		pPos->pMatElement[nj] = pS->pContextPos->pMatElement[ni];
		pMotif->pMatElement[nj] = pS->pContextMotif->pMatElement[ni];
		pStrand->pMatElement[nj] = pS->pContextStrand->pMatElement[ni];
		nj++;
	}
	if(nAdd == 0)
	{
		pPos->pMatElement[nj] = nPos;
		pMotif->pMatElement[nj] = nMotifId;
		pStrand->pMatElement[nj] = nStrand;
	}

	DestroyIntMatrix(pS->pContextPos);
	pS->pContextPos = pPos;
	DestroyIntMatrix(pS->pContextMotif);
	pS->pContextMotif = pMotif;
	DestroyIntMatrix(pS->pContextStrand);
	pS->pContextStrand = pStrand;

	/* return */
	return nNewNum;
}


/* ----------------------------------------------------------------------- */ 
/*             FLEXSEQMTFCREATE: create flex seq/motif complex             */
/* ----------------------------------------------------------------------- */ 
struct FLEXSEQMOTIF *FLEXSEQMTFCREATE(int nId, int nSeqNum, int nScoreNum, 
									  int nStatusNum, int nMonitorNum)
{
	/* new complex */
	struct FLEXSEQMOTIF *pSeqMtf;

	/* init */
	if((nSeqNum<0) || (nScoreNum<0) || (nStatusNum<0) || (nMonitorNum<0))
	{
		printf("Error: flexseqmtfcreate, parameters must be no less than 0!\n");
		exit(EXIT_FAILURE);
	}

	/* create motif matrix */
	pSeqMtf = NULL;
	pSeqMtf = (struct FLEXSEQMOTIF *)calloc(1, sizeof(struct FLEXSEQMOTIF));
	if (pSeqMtf == NULL)
	{
		printf("Error: cannot create flex seq/motif complex.\n");
		exit(EXIT_FAILURE);
	}
	
	/* init */
	pSeqMtf->nId = nId;
	strcpy(pSeqMtf->strAlias, "");
	pSeqMtf->nSeqLen = 0;

	pSeqMtf->nSeqNum = nSeqNum;
	if(nSeqNum > 0)
	{
		pSeqMtf->pSeqId = (int *)calloc(nSeqNum, sizeof(int));
		pSeqMtf->vSeq = (struct BYTEMATRIX **)calloc(nSeqNum, sizeof(struct BYTEMATRIX *));
	}

	pSeqMtf->nScoreNum = nScoreNum;
	if(nScoreNum > 0)
	{
		pSeqMtf->pScoreId = (int *)calloc(nScoreNum, sizeof(int));
		pSeqMtf->vScore = (struct BYTEMATRIX **)calloc(nScoreNum, sizeof(struct BYTEMATRIX *));
	}

	pSeqMtf->nStatusNum = nStatusNum;
	if(nStatusNum > 0)
	{
		pSeqMtf->pStatusId = (int *)calloc(nStatusNum, sizeof(int));
		pSeqMtf->vStatus = (struct BYTEMATRIX **)calloc(nStatusNum, sizeof(struct BYTEMATRIX *));
	}

	pSeqMtf->nMonitorNum = nMonitorNum;
	if(nMonitorNum > 0)
	{
		pSeqMtf->pMonitorId = (int *)calloc(nMonitorNum, sizeof(int));
		pSeqMtf->vMonitor = (struct DOUBLEMATRIX **)calloc(nMonitorNum, sizeof(struct DOUBLEMATRIX *));
	}

	pSeqMtf->vMotif = NULL;
	pSeqMtf->nSiteNum = 0;
	pSeqMtf->pMotifList = NULL;
	pSeqMtf->pTrialMotifList = NULL;

	/* return */
	return pSeqMtf;
}

/* ----------------------------------------------------------------------- */ 
/*             FLEXSEQMTFDESTROY: destroy flex seq/motif complex           */
/* ----------------------------------------------------------------------- */ 
int FLEXSEQMTFDESTROY(struct FLEXSEQMOTIF *pSeqMtf)
{
	int ni;

	/* initial check */
	if(pSeqMtf == NULL)
		return PROC_SUCCESS;

	/* destroy sequence */
	for(ni=0; ni<pSeqMtf->nSeqNum; ni++)
	{
		DestroyByteMatrix((pSeqMtf->vSeq)[ni]);
		(pSeqMtf->vSeq)[ni] = NULL;
	}
	free(pSeqMtf->vSeq);
	pSeqMtf->vSeq = NULL;
	free(pSeqMtf->pSeqId);
	pSeqMtf->pSeqId = NULL;
	pSeqMtf->nSeqNum = 0;

	/* destroy score */
	for(ni=0; ni<pSeqMtf->nScoreNum; ni++)
	{
		DestroyByteMatrix((pSeqMtf->vScore)[ni]);
		(pSeqMtf->vScore)[ni] = NULL;
	}
	free(pSeqMtf->vScore);
	pSeqMtf->vScore = NULL;
	free(pSeqMtf->pScoreId);
	pSeqMtf->pScoreId = NULL;
	pSeqMtf->nScoreNum = 0;

	/* destroy status */
	for(ni=0; ni<pSeqMtf->nStatusNum; ni++)
	{
		DestroyByteMatrix((pSeqMtf->vStatus)[ni]);
		(pSeqMtf->vStatus)[ni] = NULL;
	}
	free(pSeqMtf->vStatus);
	pSeqMtf->vStatus = NULL;
	free(pSeqMtf->pStatusId);
	pSeqMtf->pStatusId = NULL;
	pSeqMtf->nStatusNum = 0;

	/* destroy monitor */
	for(ni=0; ni<pSeqMtf->nMonitorNum; ni++)
	{
		DestroyDoubleMatrix((pSeqMtf->vMonitor)[ni]);
		(pSeqMtf->vMonitor)[ni] = NULL;
	}
	free(pSeqMtf->vMonitor);
	pSeqMtf->vMonitor = NULL;
	free(pSeqMtf->pMonitorId);
	pSeqMtf->pMonitorId = NULL;
	pSeqMtf->nMonitorNum = 0;

	/* destroy motif sites */
	if(pSeqMtf->vMotif != NULL)
	{
		free(pSeqMtf->vMotif);
	}
	FLEXMTFSDESTROYLIST(&(pSeqMtf->pMotifList));
	FLEXMTFSDESTROYLIST(&(pSeqMtf->pTrialMotifList));
	pSeqMtf->nSiteNum = 0;

	/* free the parent structure */
	free(pSeqMtf);
	pSeqMtf = NULL;

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*      FLEXSEQMTFCREATESEQ: allocate memory required for a sequence       */
/* ----------------------------------------------------------------------- */ 
int FLEXSEQMTFCREATESEQ(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, int nLen)
{
	/* check */
	if(pSeqMtf == NULL)
		return PROC_FAILURE;
	if( (nIndex >= pSeqMtf->nSeqNum) || (nIndex < 0) )
	{
		printf("Error: FLEXSEQMTFCREATESEQ, nIndex out of range!\n");
		exit(EXIT_FAILURE);
	}
	if(pSeqMtf->vSeq[nIndex] != NULL)
	{
		DestroyByteMatrix(pSeqMtf->vSeq[nIndex]);
		pSeqMtf->vSeq[nIndex] = NULL;
	}

	/* create */
	pSeqMtf->vSeq[nIndex] = CreateByteMatrix(1, nLen);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*      FLEXSEQMTFDESTROYSEQ: release memory required for a sequence       */
/* ----------------------------------------------------------------------- */ 
int FLEXSEQMTFDESTROYSEQ(struct FLEXSEQMOTIF *pSeqMtf, int nIndex)
{
	/* check */
	if(pSeqMtf == NULL)
		return PROC_FAILURE;
	if( (nIndex >= pSeqMtf->nSeqNum) || (nIndex < 0) )
	{
		printf("Error: FLEXSEQMTFDESTROYSEQ, nIndex out of range!\n");
		exit(EXIT_FAILURE);
	}
	if(pSeqMtf->vSeq[nIndex] != NULL)
	{
		DestroyByteMatrix(pSeqMtf->vSeq[nIndex]);
		pSeqMtf->vSeq[nIndex] = NULL;
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*      FLEXSEQMTFCREATESCORE: allocate memory required for a score        */
/* ----------------------------------------------------------------------- */ 
int FLEXSEQMTFCREATESCORE(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, int nLen)
{
	/* check */
	if(pSeqMtf == NULL)
		return PROC_FAILURE;
	if( (nIndex >= pSeqMtf->nScoreNum) || (nIndex < 0) )
	{
		printf("Error: FLEXSEQMTFCREATESCORE, nIndex out of range!\n");
		exit(EXIT_FAILURE);
	}
	if(pSeqMtf->vScore[nIndex] != NULL)
	{
		DestroyByteMatrix(pSeqMtf->vScore[nIndex]);
		pSeqMtf->vScore[nIndex] = NULL;
	}

	/* create */
	pSeqMtf->vScore[nIndex] = CreateByteMatrix(1, nLen);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*      FLEXSEQMTFDESTROYSCORE: release memory required for a score        */
/* ----------------------------------------------------------------------- */ 
int FLEXSEQMTFDESTROYSCORE(struct FLEXSEQMOTIF *pSeqMtf, int nIndex)
{
	/* check */
	if(pSeqMtf == NULL)
		return PROC_FAILURE;
	if( (nIndex >= pSeqMtf->nScoreNum) || (nIndex < 0) )
	{
		printf("Error: FLEXSEQMTFDESTROYSCORE, nIndex out of range!\n");
		exit(EXIT_FAILURE);
	}
	if(pSeqMtf->vScore[nIndex] != NULL)
	{
		DestroyByteMatrix(pSeqMtf->vScore[nIndex]);
		pSeqMtf->vScore[nIndex] = NULL;
	}

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*      FLEXSEQMTFCREATESTATUS: allocate memory required for a status      */
/* ----------------------------------------------------------------------- */ 
int FLEXSEQMTFCREATESTATUS(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, int nLen)
{
	/* check */
	if(pSeqMtf == NULL)
		return PROC_FAILURE;
	if( (nIndex >= pSeqMtf->nStatusNum) || (nIndex < 0) )
	{
		printf("Error: FLEXSEQMTFCREATESTATUS, nIndex out of range!\n");
		exit(EXIT_FAILURE);
	}
	if(pSeqMtf->vStatus[nIndex] != NULL)
	{
		DestroyByteMatrix(pSeqMtf->vStatus[nIndex]);
		pSeqMtf->vStatus[nIndex] = NULL;
	}

	/* create */
	pSeqMtf->vStatus[nIndex] = CreateByteMatrix(1, nLen);

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*      FLEXSEQMTFDESTROYSTATUS: release memory required for a status      */
/* ----------------------------------------------------------------------- */ 
int FLEXSEQMTFDESTROYSTATUS(struct FLEXSEQMOTIF *pSeqMtf, int nIndex)
{
	/* check */
	if(pSeqMtf == NULL)
		return PROC_FAILURE;
	if( (nIndex >= pSeqMtf->nStatusNum) || (nIndex < 0) )
	{
		printf("Error: FLEXSEQMTFDESTROYSTATUS, nIndex out of range!\n");
		exit(EXIT_FAILURE);
	}
	if(pSeqMtf->vStatus[nIndex] != NULL)
	{
		DestroyByteMatrix(pSeqMtf->vStatus[nIndex]);
		pSeqMtf->vStatus[nIndex] = NULL;
	}

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*      FLEXSEQMTFCREATEMONITOR: allocate memory required for a monitor    */
/* ----------------------------------------------------------------------- */ 
int FLEXSEQMTFCREATEMONITOR(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, int nLen)
{
	/* check */
	if(pSeqMtf == NULL)
		return PROC_FAILURE;
	if( (nIndex >= pSeqMtf->nMonitorNum) || (nIndex < 0) )
	{
		printf("Error: FLEXSEQMTFCREATEMONITOR, nIndex out of range!\n");
		exit(EXIT_FAILURE);
	}
	if(pSeqMtf->vMonitor[nIndex] != NULL)
	{
		DestroyDoubleMatrix(pSeqMtf->vMonitor[nIndex]);
		pSeqMtf->vMonitor[nIndex] = NULL;
	}

	/* create */
	pSeqMtf->vMonitor[nIndex] = CreateDoubleMatrix(1, nLen);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*      FLEXSEQMTFDESTROYMONITOR: release memory required for a monitor    */
/* ----------------------------------------------------------------------- */ 
int FLEXSEQMTFDESTROYMONITOR(struct FLEXSEQMOTIF *pSeqMtf, int nIndex)
{
	/* check */
	if(pSeqMtf == NULL)
		return PROC_FAILURE;
	if( (nIndex >= pSeqMtf->nMonitorNum) || (nIndex < 0) )
	{
		printf("Error: FLEXSEQMTFDESTROYMONITOR, nIndex out of range!\n");
		exit(EXIT_FAILURE);
	}
	if(pSeqMtf->vMonitor[nIndex] != NULL)
	{
		DestroyDoubleMatrix(pSeqMtf->vMonitor[nIndex]);
		pSeqMtf->vMonitor[nIndex] = NULL;
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*      FLEXSEQMTFCREATEMOTIF: allocate memory required for motifs         */
/* ----------------------------------------------------------------------- */ 
int FLEXSEQMTFCREATEMOTIF(struct FLEXSEQMOTIF *pSeqMtf, int nLen)
{
	/* check */
	if(pSeqMtf == NULL)
		return PROC_FAILURE;
	
	if(pSeqMtf->vMotif != NULL)
	{
		free(pSeqMtf->vMotif);
		pSeqMtf->vMotif = NULL;
	}

	/* create */
	pSeqMtf->vMotif = (struct FLEXMOTIFSITE **)calloc(nLen, sizeof(struct FLEXMOTIFSITE *));
	if(pSeqMtf->vMotif == NULL)
	{
		printf("Error: FLEXSEQMTFCREATEMOTIF, cannot create memory for tracking motifs!\n");
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*      FLEXSEQMTFDESTROYMOTIF: release memory required for motifs         */
/* ----------------------------------------------------------------------- */ 
int FLEXSEQMTFDESTROYMOTIF(struct FLEXSEQMOTIF *pSeqMtf)
{
	/* check */
	if(pSeqMtf == NULL)
		return PROC_FAILURE;

	if(pSeqMtf->vMotif != NULL)
	{
		free(pSeqMtf->vMotif);
		pSeqMtf->vMotif = NULL;
	}

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/* FLEXSEQMTFLOADSEQ: load sequence and scores from input files            */
/* Return a vector of pointers to newly created flexmotif/seq complexes,   */
/* and nSeqCount, the number of loaded complexes.                          */
/* ----------------------------------------------------------------------- */
struct FLEXSEQMOTIF **FLEXSEQMTFLOADSEQ(int *nSeqCount, char strWorkPath[], 
										char strSeqFile[],
										int nUseCS, char strCSPrefix[],
										int nMotifNum)
{
	/* define */
	int nCount;
	struct tagSequence *pSeqList;
	struct tagSequence *pCurrentSeq;
	struct FLEXSEQMOTIF **vSeqMtf;
	int ni,nj,nLen;
	int nSeqNum = 1;
	int nScoreNum = 0;
	/* status[0]: Ai, true motif; status[1]: Si, site occupied */
	int nStatusNum = 2;
	/* monitor[0]: background likelihood; monitor[1,2...]: posterior of motifs */
	int nMonitorNum = 1;
	char strFilePath[MED_LINE_LENGTH];

	/* init */
	vSeqMtf = NULL;
	nMonitorNum = 1+nMotifNum;
	if(nUseCS == 1)
		nScoreNum = 1;
	else
		nScoreNum = 0;

	/* load sequences from fasta file */
	sprintf(strFilePath, "%s%s", strWorkPath, strSeqFile);
	pSeqList = NULL;
	nCount = LoadFullSequenceList(strFilePath, &pSeqList);
	if(nCount <= 0)
		return NULL;

	/* create complexes */
	vSeqMtf = (struct FLEXSEQMOTIF **)calloc(nCount, sizeof(struct FLEXSEQMOTIF *));
	if(vSeqMtf == NULL)
		return NULL;
	
	/* code sequence */
	pCurrentSeq = pSeqList;
	for(ni=0; ni<nCount; ni++)
	{
		vSeqMtf[ni] = NULL;
		vSeqMtf[ni] = FLEXSEQMTFCREATE(ni, nSeqNum, nScoreNum, nStatusNum, nMonitorNum);
		if(vSeqMtf[ni] == NULL)
		{
			printf("Error: FLEXSEQMTFLOADSEQ, cannot create flex sequence/motif complex!\n");
			exit(EXIT_FAILURE);
		}
		nLen = pCurrentSeq->m_nLength;

		/* init length */
		strcpy(vSeqMtf[ni]->strAlias, pCurrentSeq->m_strAlias);
		vSeqMtf[ni]->nSiteNum = 0;
		vSeqMtf[ni]->nSeqLen = nLen;

		/* create space for sequence */
		for(nj=0; nj<nSeqNum; nj++)
		{
			vSeqMtf[ni]->pSeqId[nj] = nj;
			FLEXSEQMTFCREATESEQ(vSeqMtf[ni], nj, nLen);
		}

		/* create space for conservation score etc. */
		for(nj=0; nj<nScoreNum; nj++)
		{
			vSeqMtf[ni]->pScoreId[nj] = nj;
			FLEXSEQMTFCREATESCORE(vSeqMtf[ni], nj, nLen);
		}

		/* create space for sampler status */
		for(nj=0; nj<nStatusNum; nj++)
		{
			vSeqMtf[ni]->pStatusId[nj] = nj;
			FLEXSEQMTFCREATESTATUS(vSeqMtf[ni], nj, nLen);
		}

		/* create space for sampler monitors */
		for(nj=0; nj<nMonitorNum; nj++)
		{
			vSeqMtf[ni]->pMonitorId[nj] = nj;
			FLEXSEQMTFCREATEMONITOR(vSeqMtf[ni], nj, nLen);
		}

		/* create space for motif status */
		FLEXSEQMTFCREATEMOTIF(vSeqMtf[ni], nLen);

		/* Code sequence */
		FLEXSEQMTFCODENUCLEICSEQ(vSeqMtf[ni], 0, nLen, pCurrentSeq->m_pSequence->m_pString);

		/* load conservation score */
		if(nUseCS == 1)
		{
			sprintf(strFilePath, "%s%s%d_%s.cs", strWorkPath, strCSPrefix, ni, pCurrentSeq->m_strAlias);
			FLEXSEQMTFLOADCS(vSeqMtf[ni], 0, nLen, strFilePath);
		}

		/* get next seq */
		pCurrentSeq = pCurrentSeq->m_pNext;
	}

	/* destroy fasta sequences */
	SequenceListClear(&pSeqList);

	/* return */
	*nSeqCount = nCount;
	return vSeqMtf;
}

/* ----------------------------------------------------------------------- */ 
/* FLEXSEQMTFLOADSEQ_FOR_BGFIT: load sequence and scores from input files  */
/* This loading is for fitting background likelihood.                      */
/* Return a vector of pointers to newly created flexmotif/seq complexes,   */
/* and nSeqCount, the number of loaded complexes.                          */
/* ----------------------------------------------------------------------- */
struct FLEXSEQMOTIF **FLEXSEQMTFLOADSEQ_FOR_BGFIT(int *nSeqCount, char strWorkPath[], 
										char strSeqFile[],
										int nUseCS, char strCSPrefix[],
										int nMotifNum)
{
	/* define */
	int nCount;
	struct tagSequence *pSeqList;
	struct tagSequence *pCurrentSeq;
	struct FLEXSEQMOTIF **vSeqMtf;
	int ni,nj,nLen;
	int nSeqNum = 1;
	int nScoreNum = 0;
	/* status[0]: Ai, true motif; status[1]: Si, site occupied */
	int nStatusNum = 2;
	/* monitor[0]: background likelihood; monitor[1,2...]: posterior of motifs */
	int nMonitorNum = 5;
	char strFilePath[MED_LINE_LENGTH];

	/* init */
	vSeqMtf = NULL;
	nMonitorNum += 2*nMotifNum;
	nStatusNum += 2*nMotifNum;
	if(nUseCS == 1)
		nScoreNum = 1;
	else
		nScoreNum = 0;

	/* load sequences from fasta file */
	sprintf(strFilePath, "%s%s", strWorkPath, strSeqFile);
	pSeqList = NULL;
	nCount = LoadFullSequenceList(strFilePath, &pSeqList);
	if(nCount <= 0)
		return NULL;

	/* create complexes */
	vSeqMtf = (struct FLEXSEQMOTIF **)calloc(nCount, sizeof(struct FLEXSEQMOTIF *));
	if(vSeqMtf == NULL)
		return NULL;
	
	/* code sequence */
	pCurrentSeq = pSeqList;
	for(ni=0; ni<nCount; ni++)
	{
		vSeqMtf[ni] = NULL;
		vSeqMtf[ni] = FLEXSEQMTFCREATE(ni, nSeqNum, nScoreNum, nStatusNum, nMonitorNum);
		if(vSeqMtf[ni] == NULL)
		{
			printf("Error: FLEXSEQMTFLOADSEQ_FOR_BGFIT, cannot create flex sequence/motif complex!\n");
			exit(EXIT_FAILURE);
		}
		nLen = pCurrentSeq->m_nLength;

		/* init length */
		strcpy(vSeqMtf[ni]->strAlias, pCurrentSeq->m_strAlias);
		vSeqMtf[ni]->nSiteNum = 0;
		vSeqMtf[ni]->nSeqLen = nLen;

		/* create space for sequence */
		for(nj=0; nj<nSeqNum; nj++)
		{
			vSeqMtf[ni]->pSeqId[nj] = nj;
			FLEXSEQMTFCREATESEQ(vSeqMtf[ni], nj, nLen);
		}

		/* create space for conservation score etc. */
		for(nj=0; nj<nScoreNum; nj++)
		{
			vSeqMtf[ni]->pScoreId[nj] = nj;
			FLEXSEQMTFCREATESCORE(vSeqMtf[ni], nj, nLen);
		}

		/* create space for sampler status */
		for(nj=0; nj<nStatusNum; nj++)
		{
			vSeqMtf[ni]->pStatusId[nj] = nj;
			FLEXSEQMTFCREATESTATUS(vSeqMtf[ni], nj, nLen);
		}

		/* create space for sampler monitors */
		for(nj=0; nj<nMonitorNum; nj++)
		{
			vSeqMtf[ni]->pMonitorId[nj] = nj;
			FLEXSEQMTFCREATEMONITOR(vSeqMtf[ni], nj, nLen);
		}

		/* create space for motif status */
		FLEXSEQMTFCREATEMOTIF(vSeqMtf[ni], nLen);

		/* Code sequence */
		FLEXSEQMTFCODENUCLEICSEQ(vSeqMtf[ni], 0, nLen, pCurrentSeq->m_pSequence->m_pString);

		/* load conservation score */
		if(nUseCS == 1)
		{
			sprintf(strFilePath, "%s%s%d_%s.cs", strWorkPath, strCSPrefix, ni, pCurrentSeq->m_strAlias);
			FLEXSEQMTFLOADCS(vSeqMtf[ni], 0, nLen, strFilePath);
		}

		/* get next seq */
		pCurrentSeq = pCurrentSeq->m_pNext;
	}

	/* destroy fasta sequences */
	SequenceListClear(&pSeqList);

	/* return */
	*nSeqCount = nCount;
	return vSeqMtf;
}

/* ----------------------------------------------------------------------- */ 
/* FLEXSEQMTFLOADSEQ_FOR_SCANCONSENSUS_GENOME: load single sequence and    */
/* scores from genome. Return a pointer to newly created                   */
/* flexmotif/seq complexes.                                                */
/* ----------------------------------------------------------------------- */
struct FLEXSEQMOTIF *FLEXSEQMTFLOADSEQ_FOR_SCANCONSENSUS_GENOME(int nSeqId,
			char strGenomePath[], char strChr[], int nStart, int nEnd, 
			int nUseCS, char strCSPath[])
{
	/* define */
	struct FLEXSEQMOTIF *pSeqMtf;
	unsigned char *pBase;
	int nj,nLen;
	int nSeqNum = 1;
	int nScoreNum = 0;
	/* status[0]: Ai, true motif; status[1]: Si, site occupied */
	int nStatusNum = 0;
	/* monitor[0]: background likelihood; monitor[1,2...]: posterior of motifs */
	int nMonitorNum = 0;
	char strFilePath[MED_LINE_LENGTH];
	FILE *fpIn;
	/* for loading sequence */
	int nP1,nP2,nR1,nR2;
	int numread;
	int ni,nk;
	unsigned char bBase,bChar;


	/* init */
	pSeqMtf = NULL;
	nMonitorNum = 0;
	if(nUseCS == 1)
		nScoreNum = 1;
	else
		nScoreNum = 0;

	/* load sequences */
	pSeqMtf = NULL;
	pSeqMtf = FLEXSEQMTFCREATE(nSeqId, nSeqNum, nScoreNum, nStatusNum, nMonitorNum);
	if(pSeqMtf == NULL)
	{
		printf("Error: FLEXSEQMTFLOADSEQ_FOR_SCANCONSENSUS_GENOME, cannot create flex sequence/motif complex!\n");
		exit(EXIT_FAILURE);
	}
	nLen = nEnd-nStart+1;
	if(nLen <= 0)
	{
		printf("Error: FLEXSEQMTFLOADSEQ_FOR_SCANCONSENSUS_GENOME, sequence length <= 0!\n");
		exit(EXIT_FAILURE);
	}

	/* init length */
	pSeqMtf->nSiteNum = 0;
	pSeqMtf->nSeqLen = nLen;

	/* create space for sequence */
	for(nj=0; nj<nSeqNum; nj++)
	{
		pSeqMtf->pSeqId[nj] = nj;
		FLEXSEQMTFCREATESEQ(pSeqMtf, nj, nLen);
	}

	/* create space for conservation score etc. */
	for(nj=0; nj<nScoreNum; nj++)
	{
		pSeqMtf->pScoreId[nj] = nj;
		FLEXSEQMTFCREATESCORE(pSeqMtf, nj, nLen);
	}

	/* create space for sampler status */
	for(nj=0; nj<nStatusNum; nj++)
	{
		pSeqMtf->pStatusId[nj] = nj;
		FLEXSEQMTFCREATESTATUS(pSeqMtf, nj, nLen);
	}

	/* create space for sampler monitors */
	for(nj=0; nj<nMonitorNum; nj++)
	{
		pSeqMtf->pMonitorId[nj] = nj;
		FLEXSEQMTFCREATEMONITOR(pSeqMtf, nj, nLen);
	}

	/* Load sequence */
	/* init */
	pBase = pSeqMtf->vSeq[0]->pMatElement;
	sprintf(strFilePath, "%s%s.sq", strGenomePath, strChr);
	fpIn = NULL;
	fpIn = fopen(strFilePath, "rb");
	if(fpIn == NULL)
	{
		printf("Warning: FLEXSEQMTFLOADSEQ_FOR_SCANCONSENSUS_GENOME, cannot open sequence file %s!\n", strFilePath);
		FLEXSEQMTFDESTROY(pSeqMtf);
		return NULL;
	}

	/* load seq */
	nP1 = (int)(nStart/2);
	nR1 = nStart%2;
	nP2 = (int)(nEnd/2);
	nR2 = nEnd%2;

	if( fseek( fpIn, nP1, SEEK_SET ) != 0 )
	{
		printf("Error: FLEXSEQMTFLOADSEQ_FOR_SCANCONSENSUS_GENOME, cannot locate the required sequence!\n");
		exit(EXIT_FAILURE);
	}

	nk = 0;

	/* first base */
	numread = fread( &bChar, sizeof( unsigned char ), 1, fpIn );
	if(nR1 == 0)
	{
		bBase = bChar >> 4;
		switch(bBase)
		{
			case 0: pBase[nk] = 0;
				break;
			case 1: pBase[nk] = 1;
				break;
			case 2: pBase[nk] = 2;
				break;
			case 3: pBase[nk] = 3;
				break;
			case 4: pBase[nk] = 10;
				break;
			case 5: pBase[nk] = 11;
				break;
			case 6: pBase[nk] = 12;
				break;
			case 7: pBase[nk] = 13;
				break;
			default: pBase[nk] = 4;
		}
		nk++;
	}

	if(nk < nLen)
	{
		bBase = bChar & 0x0F;
		switch(bBase)
		{
			case 0: pBase[nk] = 0;
				break;
			case 1: pBase[nk] = 1;
				break;
			case 2: pBase[nk] = 2;
				break;
			case 3: pBase[nk] = 3;
				break;
			case 4: pBase[nk] = 10;
				break;
			case 5: pBase[nk] = 11;
				break;
			case 6: pBase[nk] = 12;
				break;
			case 7: pBase[nk] = 13;
				break;
			default: pBase[nk] = 4;
		}
		nk++;
	}

	/* middle bases */
	if(nk < nLen)
	{
		for(ni=(nP1+1); ni<nP2; ni++)
		{
			numread = fread( &bChar, sizeof( unsigned char ), 1, fpIn );
			bBase = bChar >> 4;
			switch(bBase)
			{
				case 0: pBase[nk] = 0;
					break;
				case 1: pBase[nk] = 1;
					break;
				case 2: pBase[nk] = 2;
					break;
				case 3: pBase[nk] = 3;
					break;
				case 4: pBase[nk] = 10;
					break;
				case 5: pBase[nk] = 11;
					break;
				case 6: pBase[nk] = 12;
					break;
				case 7: pBase[nk] = 13;
					break;
				default: pBase[nk] = 4;
			}
			nk++;

			bBase = bChar & 0x0F;
			switch(bBase)
			{
				case 0: pBase[nk] = 0;
					break;
				case 1: pBase[nk] = 1;
					break;
				case 2: pBase[nk] = 2;
					break;
				case 3: pBase[nk] = 3;
					break;
				case 4: pBase[nk] = 10;
					break;
				case 5: pBase[nk] = 11;
					break;
				case 6: pBase[nk] = 12;
					break;
				case 7: pBase[nk] = 13;
					break;
				default: pBase[nk] = 4;
			}
			nk++;
		}
	}
	
	/* last base */
	if(nk < nLen)
	{
		numread = fread( &bChar, sizeof( unsigned char ), 1, fpIn );
		bBase = bChar >> 4;
		switch(bBase)
		{
			case 0: pBase[nk] = 0;
				break;
			case 1: pBase[nk] = 1;
				break;
			case 2: pBase[nk] = 2;
				break;
			case 3: pBase[nk] = 3;
				break;
			case 4: pBase[nk] = 10;
				break;
			case 5: pBase[nk] = 11;
				break;
			case 6: pBase[nk] = 12;
				break;
			case 7: pBase[nk] = 13;
				break;
			default: pBase[nk] = 4;
		}
		nk++;
	}

	if(nk < nLen)
	{
		if(nR2 == 1)
		{
			bBase = bChar & 0x0F;
			switch(bBase)
			{
				case 0: pBase[nk] = 0;
					break;
				case 1: pBase[nk] = 1;
					break;
				case 2: pBase[nk] = 2;
					break;
				case 3: pBase[nk] = 3;
					break;
				case 4: pBase[nk] = 10;
					break;
				case 5: pBase[nk] = 11;
					break;
				case 6: pBase[nk] = 12;
					break;
				case 7: pBase[nk] = 13;
					break;
				default: pBase[nk] = 4;
			}
			nk++;
		}
	}

	/* close file */
	fclose(fpIn);

	if(nk != nLen)
	{
		printf("Error: FLEXSEQMTFLOADSEQ_FOR_SCANCONSENSUS_GENOME, incorrect sequence length!\n");
		printf("%s:%d-%d(%d=%d)\n", strChr, nStart, nEnd, nLen, nk);
		exit(EXIT_FAILURE);
	}


	/* load conservation score */
	if(nUseCS == 1)
	{
		sprintf(strFilePath, "%s%s.cs", strCSPath, strChr);
		fpIn = NULL;
		fpIn = fopen(strFilePath, "rb");
		if(fpIn == NULL)
		{
			printf("Error: FLEXSEQMTFLOADSEQ_FOR_SCANCONSENSUS_GENOME, cannot open conservation score file %s!\n", strFilePath);
			FLEXSEQMTFDESTROY(pSeqMtf);
			return NULL;
		}

		/* load score */
		if( fseek( fpIn, nStart, SEEK_SET ) != 0 )
		{
			printf("Error: FLEXSEQMTFLOADSEQ_FOR_SCANCONSENSUS_GENOME, cannot locate the required sequence!\n");
			exit(EXIT_FAILURE);
		}

		numread = fread((pSeqMtf->vScore[0])->pMatElement, sizeof(unsigned char), nLen, fpIn);
		if(numread != nLen)
		{
			printf("Error: FLEXSEQMTFLOADSEQ_FOR_SCANCONSENSUS_GENOME, loading error!\n");
			exit(EXIT_FAILURE);
		}

		fclose(fpIn);
	}

	/* return */
	return pSeqMtf;
}


/* ----------------------------------------------------------------------- */ 
/* FLEXSEQMTFLOADSEQ_FOR_SCANCONSENSUS_SINGLE: load single sequence and    */
/* scores from input sequences. Return a pointer to newly created          */
/* flexmotif/seq complexes.                                                */
/* ----------------------------------------------------------------------- */
struct FLEXSEQMOTIF *FLEXSEQMTFLOADSEQ_FOR_SCANCONSENSUS_SINGLE(struct tagSequence *pSeq,
					char strWorkPath[], int nUseCS, char strCSPrefix[])
{
	/* define */
	struct FLEXSEQMOTIF *pSeqMtf;
	int nj,nLen;
	int nSeqNum = 1;
	int nScoreNum = 0;
	/* status[0]: Ai, true motif; status[1]: Si, site occupied */
	int nStatusNum = 0;
	/* monitor[0]: background likelihood; monitor[1,2...]: posterior of motifs */
	int nMonitorNum = 0;
	char strFilePath[MED_LINE_LENGTH];

	/* init */
	if(pSeq == NULL)	
	{
		printf("Warning: FLEXSEQMTFLOADSEQ_FOR_SCANCONSENSUS_SINGLE, empty sequences!\n");
		return NULL;
	}

	pSeqMtf = NULL;
	nMonitorNum = 0;
	if(nUseCS == 1)
		nScoreNum = 1;
	else
		nScoreNum = 0;

	/* load sequences from fasta file */
	pSeqMtf = NULL;
	pSeqMtf = FLEXSEQMTFCREATE(pSeq->m_nIndex, nSeqNum, nScoreNum, nStatusNum, nMonitorNum);
	if(pSeqMtf == NULL)
	{
		printf("Error: FLEXSEQMTFLOADSEQ_FOR_SCANCONSENSUS, cannot create flex sequence/motif complex!\n");
		exit(EXIT_FAILURE);
	}
	nLen = pSeq->m_nLength;

	/* init length */
	pSeqMtf->nSiteNum = 0;
	pSeqMtf->nSeqLen = nLen;

	/* create space for sequence */
	for(nj=0; nj<nSeqNum; nj++)
	{
		pSeqMtf->pSeqId[nj] = nj;
		FLEXSEQMTFCREATESEQ(pSeqMtf, nj, nLen);
	}

	/* create space for conservation score etc. */
	for(nj=0; nj<nScoreNum; nj++)
	{
		pSeqMtf->pScoreId[nj] = nj;
		FLEXSEQMTFCREATESCORE(pSeqMtf, nj, nLen);
	}

	/* create space for sampler status */
	for(nj=0; nj<nStatusNum; nj++)
	{
		pSeqMtf->pStatusId[nj] = nj;
		FLEXSEQMTFCREATESTATUS(pSeqMtf, nj, nLen);
	}

	/* create space for sampler monitors */
	for(nj=0; nj<nMonitorNum; nj++)
	{
		pSeqMtf->pMonitorId[nj] = nj;
		FLEXSEQMTFCREATEMONITOR(pSeqMtf, nj, nLen);
	}

	/* Code sequence */
	FLEXSEQMTFCODENUCLEICSEQ(pSeqMtf, 0, nLen, pSeq->m_pSequence->m_pString);

	/* load conservation score */
	if(nUseCS == 1)
	{
		sprintf(strFilePath, "%s%s%d_%s.cs", strWorkPath, strCSPrefix, pSeq->m_nIndex, pSeq->m_strAlias);
		FLEXSEQMTFLOADCS(pSeqMtf, 0, nLen, strFilePath);
	}

	/* return */
	return pSeqMtf;
}

/* ----------------------------------------------------------------------- */ 
/* FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_SINGLE: load single sequence and       */
/* scores from input sequences. Return a pointer to newly created          */
/* flexmotif/seq complexes.                                                */
/* ----------------------------------------------------------------------- */
struct FLEXSEQMOTIF *FLEXSEQMTFLOADSEQ_FOR_COUNTMCBACKGROUND_GENOME(int nSeqId,
				char strGenomePath[], char strChr[], int nStart, int nEnd,
				int nUseCS, char strCSPath[])
{
	/* define */
	struct FLEXSEQMOTIF *pSeqMtf;
	unsigned char *pBase;
	int nj,nLen;
	int nSeqNum = 1;
	int nScoreNum = 0;
	/* status[0]: Ai, true motif; status[1]: Si, site occupied */
	int nStatusNum = 0;
	/* monitor[0]: background likelihood; monitor[1,2...]: posterior of motifs */
	int nMonitorNum = 0;
	char strFilePath[MED_LINE_LENGTH];
	FILE *fpIn;
	/* for loading sequence */
	int nP1,nP2,nR1,nR2;
	int numread;
	int ni,nk;
	unsigned char bBase,bChar;

	/* init */
	if(nUseCS == 1)
		nScoreNum = 1;
	else
		nScoreNum = 0;

	/* load sequences */
	pSeqMtf = NULL;
	pSeqMtf = FLEXSEQMTFCREATE(nSeqId, nSeqNum, nScoreNum, nStatusNum, nMonitorNum);
	if(pSeqMtf == NULL)
	{
		printf("Error: FLEXSEQMTFLOADSEQ_FOR_COUNTMCBACKGROUND_GENOME, cannot create flex sequence/motif complex!\n");
		exit(EXIT_FAILURE);
	}
	nLen = nEnd-nStart+1;
	if(nLen <= 0)
	{
		printf("Error: FLEXSEQMTFLOADSEQ_FOR_COUNTMCBACKGROUND_GENOME, sequence length <= 0!\n");
		exit(EXIT_FAILURE);
	}

	/* init length */
	pSeqMtf->nSiteNum = 0;
	pSeqMtf->nSeqLen = nLen;

	/* create space for sequence */
	for(nj=0; nj<nSeqNum; nj++)
	{
		pSeqMtf->pSeqId[nj] = nj;
		FLEXSEQMTFCREATESEQ(pSeqMtf, nj, nLen);
	}

	/* create space for conservation score etc. */
	for(nj=0; nj<nScoreNum; nj++)
	{
		pSeqMtf->pScoreId[nj] = nj;
		FLEXSEQMTFCREATESCORE(pSeqMtf, nj, nLen);
	}

	/* create space for sampler status */
	for(nj=0; nj<nStatusNum; nj++)
	{
		pSeqMtf->pStatusId[nj] = nj;
		FLEXSEQMTFCREATESTATUS(pSeqMtf, nj, nLen);
	}

	/* create space for sampler monitors */
	for(nj=0; nj<nMonitorNum; nj++)
	{
		pSeqMtf->pMonitorId[nj] = nj;
		FLEXSEQMTFCREATEMONITOR(pSeqMtf, nj, nLen);
	}

	/* Load sequence */
	/* init */
	pBase = pSeqMtf->vSeq[0]->pMatElement;
	sprintf(strFilePath, "%s%s.sq", strGenomePath, strChr);
	fpIn = NULL;
	fpIn = fopen(strFilePath, "rb");
	if(fpIn == NULL)
	{
		printf("Warning: FLEXSEQMTFLOADSEQ_FOR_COUNTMCBACKGROUND_GENOME, cannot open sequence file %s!\n", strFilePath);
		FLEXSEQMTFDESTROY(pSeqMtf);
		return NULL;
	}

	/* load seq */
	nP1 = (int)(nStart/2);
	nR1 = nStart%2;
	nP2 = (int)(nEnd/2);
	nR2 = nEnd%2;

	if( fseek( fpIn, nP1, SEEK_SET ) != 0 )
	{
		printf("Error: FLEXSEQMTFLOADSEQ_FOR_COUNTMCBACKGROUND_GENOME, cannot locate the required sequence!\n");
		exit(EXIT_FAILURE);
	}

	nk = 0;

	/* first base */
	numread = fread( &bChar, sizeof( unsigned char ), 1, fpIn );
	if(nR1 == 0)
	{
		bBase = bChar >> 4;
		switch(bBase)
		{
			case 0: pBase[nk] = 0;
				break;
			case 1: pBase[nk] = 1;
				break;
			case 2: pBase[nk] = 2;
				break;
			case 3: pBase[nk] = 3;
				break;
			case 4: pBase[nk] = 10;
				break;
			case 5: pBase[nk] = 11;
				break;
			case 6: pBase[nk] = 12;
				break;
			case 7: pBase[nk] = 13;
				break;
			default: pBase[nk] = 4;
		}
		nk++;
	}

	if(nk < nLen)
	{
		bBase = bChar & 0x0F;
		switch(bBase)
		{
			case 0: pBase[nk] = 0;
				break;
			case 1: pBase[nk] = 1;
				break;
			case 2: pBase[nk] = 2;
				break;
			case 3: pBase[nk] = 3;
				break;
			case 4: pBase[nk] = 10;
				break;
			case 5: pBase[nk] = 11;
				break;
			case 6: pBase[nk] = 12;
				break;
			case 7: pBase[nk] = 13;
				break;
			default: pBase[nk] = 4;
		}
		nk++;
	}

	/* middle bases */
	if(nk < nLen)
	{
		for(ni=(nP1+1); ni<nP2; ni++)
		{
			numread = fread( &bChar, sizeof( unsigned char ), 1, fpIn );
			bBase = bChar >> 4;
			switch(bBase)
			{
				case 0: pBase[nk] = 0;
					break;
				case 1: pBase[nk] = 1;
					break;
				case 2: pBase[nk] = 2;
					break;
				case 3: pBase[nk] = 3;
					break;
				case 4: pBase[nk] = 10;
					break;
				case 5: pBase[nk] = 11;
					break;
				case 6: pBase[nk] = 12;
					break;
				case 7: pBase[nk] = 13;
					break;
				default: pBase[nk] = 4;
			}
			nk++;

			bBase = bChar & 0x0F;
			switch(bBase)
			{
				case 0: pBase[nk] = 0;
					break;
				case 1: pBase[nk] = 1;
					break;
				case 2: pBase[nk] = 2;
					break;
				case 3: pBase[nk] = 3;
					break;
				case 4: pBase[nk] = 10;
					break;
				case 5: pBase[nk] = 11;
					break;
				case 6: pBase[nk] = 12;
					break;
				case 7: pBase[nk] = 13;
					break;
				default: pBase[nk] = 4;
			}
			nk++;
		}
	}
	
	/* last base */
	if(nk < nLen)
	{
		numread = fread( &bChar, sizeof( unsigned char ), 1, fpIn );
		bBase = bChar >> 4;
		switch(bBase)
		{
			case 0: pBase[nk] = 0;
				break;
			case 1: pBase[nk] = 1;
				break;
			case 2: pBase[nk] = 2;
				break;
			case 3: pBase[nk] = 3;
				break;
			case 4: pBase[nk] = 10;
				break;
			case 5: pBase[nk] = 11;
				break;
			case 6: pBase[nk] = 12;
				break;
			case 7: pBase[nk] = 13;
				break;
			default: pBase[nk] = 4;
		}
		nk++;
	}

	if(nk < nLen)
	{
		if(nR2 == 1)
		{
			bBase = bChar & 0x0F;
			switch(bBase)
			{
				case 0: pBase[nk] = 0;
					break;
				case 1: pBase[nk] = 1;
					break;
				case 2: pBase[nk] = 2;
					break;
				case 3: pBase[nk] = 3;
					break;
				case 4: pBase[nk] = 10;
					break;
				case 5: pBase[nk] = 11;
					break;
				case 6: pBase[nk] = 12;
					break;
				case 7: pBase[nk] = 13;
					break;
				default: pBase[nk] = 4;
			}
			nk++;
		}
	}

	/* close file */
	fclose(fpIn);

	if(nk != nLen)
	{
		printf("Error: FLEXSEQMTFLOADSEQ_FOR_COUNTMCBACKGROUND_GENOME, incorrect sequence length!\n");
		exit(EXIT_FAILURE);
	}


	/* load conservation score */
	if(nUseCS == 1)
	{
		sprintf(strFilePath, "%s%s.cs", strCSPath, strChr);
		fpIn = NULL;
		fpIn = fopen(strFilePath, "rb");
		if(fpIn == NULL)
		{
			printf("Error: FLEXSEQMTFLOADSEQ_FOR_COUNTMCBACKGROUND_GENOME, cannot open conservation score file %s!\n", strFilePath);
			FLEXSEQMTFDESTROY(pSeqMtf);
			return NULL;
		}

		/* load score */
		if( fseek( fpIn, nStart, SEEK_SET ) != 0 )
		{
			printf("Error: FLEXSEQMTFLOADSEQ_FOR_COUNTMCBACKGROUND_GENOME, cannot locate the required sequence!\n");
			exit(EXIT_FAILURE);
		}

		numread = fread((pSeqMtf->vScore[0])->pMatElement, sizeof(unsigned char), nLen, fpIn);
		if(numread != nLen)
		{
			printf("Error: FLEXSEQMTFLOADSEQ_FOR_COUNTMCBACKGROUND_GENOME, loading error!\n");
			exit(EXIT_FAILURE);
		}

		fclose(fpIn);
	}

	/* return */
	return pSeqMtf;
}

/* ----------------------------------------------------------------------- */ 
/* FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_GENOME: load single sequence and       */
/* scores from input sequences. Return a pointer to newly created          */
/* flexmotif/seq complexes.                                                */
/* ----------------------------------------------------------------------- */
struct FLEXSEQMOTIF *FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_GENOME(int nSeqId,
				char strGenomePath[], char strChr[], int nStart, int nEnd,
				int nUseCS, char strCSPath[])
{
	/* define */
	struct FLEXSEQMOTIF *pSeqMtf;
	unsigned char *pBase;
	int nj,nLen;
	int nSeqNum = 1;
	int nScoreNum = 0;
	/* status[0]: Ai, true motif; status[1]: Si, site occupied */
	int nStatusNum = 0;
	/* monitor[0]: background likelihood; monitor[1,2...]: posterior of motifs */
	int nMonitorNum = 2;
	char strFilePath[MED_LINE_LENGTH];
	FILE *fpIn;
	/* for loading sequence */
	int nP1,nP2,nR1,nR2;
	int numread;
	int ni,nk;
	unsigned char bBase,bChar;

	/* init */
	if(nUseCS == 1)
		nScoreNum = 1;
	else
		nScoreNum = 0;

	/* load sequences */
	pSeqMtf = NULL;
	pSeqMtf = FLEXSEQMTFCREATE(nSeqId, nSeqNum, nScoreNum, nStatusNum, nMonitorNum);
	if(pSeqMtf == NULL)
	{
		printf("Error: FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_GENOME, cannot create flex sequence/motif complex!\n");
		exit(EXIT_FAILURE);
	}
	nLen = nEnd-nStart+1;
	if(nLen <= 0)
	{
		printf("Error: FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_GENOME, sequence length <= 0!\n");
		exit(EXIT_FAILURE);
	}

	/* init length */
	pSeqMtf->nSiteNum = 0;
	pSeqMtf->nSeqLen = nLen;

	/* create space for sequence */
	for(nj=0; nj<nSeqNum; nj++)
	{
		pSeqMtf->pSeqId[nj] = nj;
		FLEXSEQMTFCREATESEQ(pSeqMtf, nj, nLen);
	}

	/* create space for conservation score etc. */
	for(nj=0; nj<nScoreNum; nj++)
	{
		pSeqMtf->pScoreId[nj] = nj;
		FLEXSEQMTFCREATESCORE(pSeqMtf, nj, nLen);
	}

	/* create space for sampler status */
	for(nj=0; nj<nStatusNum; nj++)
	{
		pSeqMtf->pStatusId[nj] = nj;
		FLEXSEQMTFCREATESTATUS(pSeqMtf, nj, nLen);
	}

	/* create space for sampler monitors */
	for(nj=0; nj<nMonitorNum; nj++)
	{
		pSeqMtf->pMonitorId[nj] = nj;
		FLEXSEQMTFCREATEMONITOR(pSeqMtf, nj, nLen);
	}

	/* Load sequence */
	/* init */
	pBase = pSeqMtf->vSeq[0]->pMatElement;
	sprintf(strFilePath, "%s%s.sq", strGenomePath, strChr);
	fpIn = NULL;
	fpIn = fopen(strFilePath, "rb");
	if(fpIn == NULL)
	{
		printf("Warning: FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_GENOME, cannot open sequence file %s!\n", strFilePath);
		FLEXSEQMTFDESTROY(pSeqMtf);
		return NULL;
	}

	/* load seq */
	nP1 = (int)(nStart/2);
	nR1 = nStart%2;
	nP2 = (int)(nEnd/2);
	nR2 = nEnd%2;

	if( fseek( fpIn, nP1, SEEK_SET ) != 0 )
	{
		printf("Error: FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_GENOME, cannot locate the required sequence!\n");
		exit(EXIT_FAILURE);
	}

	nk = 0;

	/* first base */
	numread = fread( &bChar, sizeof( unsigned char ), 1, fpIn );
	if(nR1 == 0)
	{
		bBase = bChar >> 4;
		switch(bBase)
		{
			case 0: pBase[nk] = 0;
				break;
			case 1: pBase[nk] = 1;
				break;
			case 2: pBase[nk] = 2;
				break;
			case 3: pBase[nk] = 3;
				break;
			case 4: pBase[nk] = 10;
				break;
			case 5: pBase[nk] = 11;
				break;
			case 6: pBase[nk] = 12;
				break;
			case 7: pBase[nk] = 13;
				break;
			default: pBase[nk] = 4;
		}
		nk++;
	}

	if(nk < nLen)
	{
		bBase = bChar & 0x0F;
		switch(bBase)
		{
			case 0: pBase[nk] = 0;
				break;
			case 1: pBase[nk] = 1;
				break;
			case 2: pBase[nk] = 2;
				break;
			case 3: pBase[nk] = 3;
				break;
			case 4: pBase[nk] = 10;
				break;
			case 5: pBase[nk] = 11;
				break;
			case 6: pBase[nk] = 12;
				break;
			case 7: pBase[nk] = 13;
				break;
			default: pBase[nk] = 4;
		}
		nk++;
	}

	/* middle bases */
	if(nk < nLen)
	{
		for(ni=(nP1+1); ni<nP2; ni++)
		{
			numread = fread( &bChar, sizeof( unsigned char ), 1, fpIn );
			bBase = bChar >> 4;
			switch(bBase)
			{
				case 0: pBase[nk] = 0;
					break;
				case 1: pBase[nk] = 1;
					break;
				case 2: pBase[nk] = 2;
					break;
				case 3: pBase[nk] = 3;
					break;
				case 4: pBase[nk] = 10;
					break;
				case 5: pBase[nk] = 11;
					break;
				case 6: pBase[nk] = 12;
					break;
				case 7: pBase[nk] = 13;
					break;
				default: pBase[nk] = 4;
			}
			nk++;

			bBase = bChar & 0x0F;
			switch(bBase)
			{
				case 0: pBase[nk] = 0;
					break;
				case 1: pBase[nk] = 1;
					break;
				case 2: pBase[nk] = 2;
					break;
				case 3: pBase[nk] = 3;
					break;
				case 4: pBase[nk] = 10;
					break;
				case 5: pBase[nk] = 11;
					break;
				case 6: pBase[nk] = 12;
					break;
				case 7: pBase[nk] = 13;
					break;
				default: pBase[nk] = 4;
			}
			nk++;
		}
	}
	
	/* last base */
	if(nk < nLen)
	{
		numread = fread( &bChar, sizeof( unsigned char ), 1, fpIn );
		bBase = bChar >> 4;
		switch(bBase)
		{
			case 0: pBase[nk] = 0;
				break;
			case 1: pBase[nk] = 1;
				break;
			case 2: pBase[nk] = 2;
				break;
			case 3: pBase[nk] = 3;
				break;
			case 4: pBase[nk] = 10;
				break;
			case 5: pBase[nk] = 11;
				break;
			case 6: pBase[nk] = 12;
				break;
			case 7: pBase[nk] = 13;
				break;
			default: pBase[nk] = 4;
		}
		nk++;
	}

	if(nk < nLen)
	{
		if(nR2 == 1)
		{
			bBase = bChar & 0x0F;
			switch(bBase)
			{
				case 0: pBase[nk] = 0;
					break;
				case 1: pBase[nk] = 1;
					break;
				case 2: pBase[nk] = 2;
					break;
				case 3: pBase[nk] = 3;
					break;
				case 4: pBase[nk] = 10;
					break;
				case 5: pBase[nk] = 11;
					break;
				case 6: pBase[nk] = 12;
					break;
				case 7: pBase[nk] = 13;
					break;
				default: pBase[nk] = 4;
			}
			nk++;
		}
	}

	/* close file */
	fclose(fpIn);

	if(nk != nLen)
	{
		printf("Error: FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_GENOME, incorrect sequence length!\n");
		exit(EXIT_FAILURE);
	}


	/* load conservation score */
	if(nUseCS == 1)
	{
		sprintf(strFilePath, "%s%s.cs", strCSPath, strChr);
		fpIn = NULL;
		fpIn = fopen(strFilePath, "rb");
		if(fpIn == NULL)
		{
			printf("Error: FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_GENOME, cannot open conservation score file %s!\n", strFilePath);
			FLEXSEQMTFDESTROY(pSeqMtf);
			return NULL;
		}

		/* load score */
		if( fseek( fpIn, nStart, SEEK_SET ) != 0 )
		{
			printf("Error: FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_GENOME, cannot locate the required sequence!\n");
			exit(EXIT_FAILURE);
		}

		numread = fread((pSeqMtf->vScore[0])->pMatElement, sizeof(unsigned char), nLen, fpIn);
		if(numread != nLen)
		{
			printf("Error: FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_GENOME, loading error!\n");
			exit(EXIT_FAILURE);
		}

		fclose(fpIn);
	}

	/* return */
	return pSeqMtf;
}

/* ----------------------------------------------------------------------- */ 
/* FLEXSEQMTFLOADSEQ_FOR_INTEGRATEMATRIX_GENOME: load single sequence and  */
/* scores from input sequences. Return a pointer to newly created          */
/* flexmotif/seq complexes.                                                */
/* ----------------------------------------------------------------------- */
struct FLEXSEQMOTIF *FLEXSEQMTFLOADSEQ_FOR_INTEGRATEMATRIX_GENOME(int nSeqId,
				char strGenomePath[], char strChr[], int nStart, int nEnd,
				int nUseCS, char strCSPath[])
{
	/* define */
	struct FLEXSEQMOTIF *pSeqMtf;
	unsigned char *pBase;
	int nj,nLen;
	int nSeqNum = 1;
	int nScoreNum = 0;
	/* status[0]: Ai, true motif; status[1]: Si, site occupied */
	int nStatusNum = 0;
	/* monitor[0]: background likelihood; monitor[1,2...]: posterior of motifs */
	int nMonitorNum = 5;
	char strFilePath[MED_LINE_LENGTH];
	FILE *fpIn;
	/* for loading sequence */
	int nP1,nP2,nR1,nR2;
	int numread;
	int ni,nk;
	unsigned char bBase,bChar;

	/* init */
	if(nUseCS == 1)
		nScoreNum = 1;
	else
		nScoreNum = 0;

	/* load sequences */
	pSeqMtf = NULL;
	pSeqMtf = FLEXSEQMTFCREATE(nSeqId, nSeqNum, nScoreNum, nStatusNum, nMonitorNum);
	if(pSeqMtf == NULL)
	{
		printf("Error: FLEXSEQMTFLOADSEQ_FOR_INTEGRATEMATRIX_GENOME, cannot create flex sequence/motif complex!\n");
		exit(EXIT_FAILURE);
	}
	nLen = nEnd-nStart+1;
	if(nLen <= 0)
	{
		printf("Error: FLEXSEQMTFLOADSEQ_FOR_INTEGRATEMATRIX_GENOME, sequence length <= 0!\n");
		exit(EXIT_FAILURE);
	}

	/* init length */
	pSeqMtf->nSiteNum = 0;
	pSeqMtf->nSeqLen = nLen;

	/* create space for sequence */
	for(nj=0; nj<nSeqNum; nj++)
	{
		pSeqMtf->pSeqId[nj] = nj;
		FLEXSEQMTFCREATESEQ(pSeqMtf, nj, nLen);
	}

	/* create space for conservation score etc. */
	for(nj=0; nj<nScoreNum; nj++)
	{
		pSeqMtf->pScoreId[nj] = nj;
		FLEXSEQMTFCREATESCORE(pSeqMtf, nj, nLen);
	}

	/* create space for sampler status */
	for(nj=0; nj<nStatusNum; nj++)
	{
		pSeqMtf->pStatusId[nj] = nj;
		FLEXSEQMTFCREATESTATUS(pSeqMtf, nj, nLen);
	}

	/* create space for sampler monitors */
	for(nj=0; nj<nMonitorNum; nj++)
	{
		pSeqMtf->pMonitorId[nj] = nj;
		FLEXSEQMTFCREATEMONITOR(pSeqMtf, nj, nLen);
	}

	/* Load sequence */
	/* init */
	pBase = pSeqMtf->vSeq[0]->pMatElement;
	sprintf(strFilePath, "%s%s.sq", strGenomePath, strChr);
	fpIn = NULL;
	fpIn = fopen(strFilePath, "rb");
	if(fpIn == NULL)
	{
		printf("Warning: FLEXSEQMTFLOADSEQ_FOR_INTEGRATEMATRIX_GENOME, cannot open sequence file %s!\n", strFilePath);
		FLEXSEQMTFDESTROY(pSeqMtf);
		return NULL;
	}

	/* load seq */
	nP1 = (int)(nStart/2);
	nR1 = nStart%2;
	nP2 = (int)(nEnd/2);
	nR2 = nEnd%2;

	if( fseek( fpIn, nP1, SEEK_SET ) != 0 )
	{
		printf("Error: FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_GENOME, cannot locate the required sequence!\n");
		exit(EXIT_FAILURE);
	}

	nk = 0;

	/* first base */
	numread = fread( &bChar, sizeof( unsigned char ), 1, fpIn );
	if(nR1 == 0)
	{
		bBase = bChar >> 4;
		switch(bBase)
		{
			case 0: pBase[nk] = 0;
				break;
			case 1: pBase[nk] = 1;
				break;
			case 2: pBase[nk] = 2;
				break;
			case 3: pBase[nk] = 3;
				break;
			case 4: pBase[nk] = 10;
				break;
			case 5: pBase[nk] = 11;
				break;
			case 6: pBase[nk] = 12;
				break;
			case 7: pBase[nk] = 13;
				break;
			default: pBase[nk] = 4;
		}
		nk++;
	}

	if(nk < nLen)
	{
		bBase = bChar & 0x0F;
		switch(bBase)
		{
			case 0: pBase[nk] = 0;
				break;
			case 1: pBase[nk] = 1;
				break;
			case 2: pBase[nk] = 2;
				break;
			case 3: pBase[nk] = 3;
				break;
			case 4: pBase[nk] = 10;
				break;
			case 5: pBase[nk] = 11;
				break;
			case 6: pBase[nk] = 12;
				break;
			case 7: pBase[nk] = 13;
				break;
			default: pBase[nk] = 4;
		}
		nk++;
	}

	/* middle bases */
	if(nk < nLen)
	{
		for(ni=(nP1+1); ni<nP2; ni++)
		{
			numread = fread( &bChar, sizeof( unsigned char ), 1, fpIn );
			bBase = bChar >> 4;
			switch(bBase)
			{
				case 0: pBase[nk] = 0;
					break;
				case 1: pBase[nk] = 1;
					break;
				case 2: pBase[nk] = 2;
					break;
				case 3: pBase[nk] = 3;
					break;
				case 4: pBase[nk] = 10;
					break;
				case 5: pBase[nk] = 11;
					break;
				case 6: pBase[nk] = 12;
					break;
				case 7: pBase[nk] = 13;
					break;
				default: pBase[nk] = 4;
			}
			nk++;

			bBase = bChar & 0x0F;
			switch(bBase)
			{
				case 0: pBase[nk] = 0;
					break;
				case 1: pBase[nk] = 1;
					break;
				case 2: pBase[nk] = 2;
					break;
				case 3: pBase[nk] = 3;
					break;
				case 4: pBase[nk] = 10;
					break;
				case 5: pBase[nk] = 11;
					break;
				case 6: pBase[nk] = 12;
					break;
				case 7: pBase[nk] = 13;
					break;
				default: pBase[nk] = 4;
			}
			nk++;
		}
	}
	
	/* last base */
	if(nk < nLen)
	{
		numread = fread( &bChar, sizeof( unsigned char ), 1, fpIn );
		bBase = bChar >> 4;
		switch(bBase)
		{
			case 0: pBase[nk] = 0;
				break;
			case 1: pBase[nk] = 1;
				break;
			case 2: pBase[nk] = 2;
				break;
			case 3: pBase[nk] = 3;
				break;
			case 4: pBase[nk] = 10;
				break;
			case 5: pBase[nk] = 11;
				break;
			case 6: pBase[nk] = 12;
				break;
			case 7: pBase[nk] = 13;
				break;
			default: pBase[nk] = 4;
		}
		nk++;
	}

	if(nk < nLen)
	{
		if(nR2 == 1)
		{
			bBase = bChar & 0x0F;
			switch(bBase)
			{
				case 0: pBase[nk] = 0;
					break;
				case 1: pBase[nk] = 1;
					break;
				case 2: pBase[nk] = 2;
					break;
				case 3: pBase[nk] = 3;
					break;
				case 4: pBase[nk] = 10;
					break;
				case 5: pBase[nk] = 11;
					break;
				case 6: pBase[nk] = 12;
					break;
				case 7: pBase[nk] = 13;
					break;
				default: pBase[nk] = 4;
			}
			nk++;
		}
	}

	/* close file */
	fclose(fpIn);

	if(nk != nLen)
	{
		printf("Error: FLEXSEQMTFLOADSEQ_FOR_INTEGRATEMATRIX_GENOME, incorrect sequence length!\n");
		exit(EXIT_FAILURE);
	}


	/* load conservation score */
	if(nUseCS == 1)
	{
		sprintf(strFilePath, "%s%s.cs", strCSPath, strChr);
		fpIn = NULL;
		fpIn = fopen(strFilePath, "rb");
		if(fpIn == NULL)
		{
			printf("Error: FLEXSEQMTFLOADSEQ_FOR_INTEGRATEMATRIX_GENOME, cannot open conservation score file %s!\n", strFilePath);
			FLEXSEQMTFDESTROY(pSeqMtf);
			return NULL;
		}

		/* load score */
		if( fseek( fpIn, nStart, SEEK_SET ) != 0 )
		{
			printf("Error: FLEXSEQMTFLOADSEQ_FOR_INTEGRATEMATRIX_GENOME, cannot locate the required sequence!\n");
			exit(EXIT_FAILURE);
		}

		numread = fread((pSeqMtf->vScore[0])->pMatElement, sizeof(unsigned char), nLen, fpIn);
		if(numread != nLen)
		{
			printf("Error: FLEXSEQMTFLOADSEQ_FOR_INTEGRATEMATRIX_GENOME, loading error!\n");
			exit(EXIT_FAILURE);
		}

		fclose(fpIn);
	}

	/* return */
	return pSeqMtf;
}

/* ----------------------------------------------------------------------- */ 
/* FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_SINGLE: load single sequence and       */
/* scores from input sequences. Return a pointer to newly created          */
/* flexmotif/seq complexes.                                                */
/* ----------------------------------------------------------------------- */
struct FLEXSEQMOTIF *FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_SINGLE(struct tagSequence *pSeq,
					char strWorkPath[], int nUseCS, char strCSPrefix[])
{
	/* define */
	struct FLEXSEQMOTIF *pSeqMtf;
	int nj,nLen;
	int nSeqNum = 1;
	int nScoreNum = 0;
	/* status[0]: Ai, true motif; status[1]: Si, site occupied */
	int nStatusNum = 0;
	/* monitor[0] and monitor[1]: background likelihood; monitor[2,3...]: posterior of motifs */
	int nMonitorNum = 2;
	char strFilePath[MED_LINE_LENGTH];

	/* init */
	if(pSeq == NULL)	
	{
		printf("Warning: FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_SINGLE, empty sequences!\n");
		return NULL;
	}

	pSeqMtf = NULL;
	nMonitorNum = 2;
	if(nUseCS == 1)
		nScoreNum = 1;
	else
		nScoreNum = 0;

	/* load sequences from fasta file */
	pSeqMtf = NULL;
	pSeqMtf = FLEXSEQMTFCREATE(pSeq->m_nIndex, nSeqNum, nScoreNum, nStatusNum, nMonitorNum);
	if(pSeqMtf == NULL)
	{
		printf("Error: FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_SINGLE, cannot create flex sequence/motif complex!\n");
		exit(EXIT_FAILURE);
	}
	nLen = pSeq->m_nLength;

	/* init length */
	pSeqMtf->nSiteNum = 0;
	pSeqMtf->nSeqLen = nLen;

	/* create space for sequence */
	for(nj=0; nj<nSeqNum; nj++)
	{
		pSeqMtf->pSeqId[nj] = nj;
		FLEXSEQMTFCREATESEQ(pSeqMtf, nj, nLen);
	}

	/* create space for conservation score etc. */
	for(nj=0; nj<nScoreNum; nj++)
	{
		pSeqMtf->pScoreId[nj] = nj;
		FLEXSEQMTFCREATESCORE(pSeqMtf, nj, nLen);
	}

	/* create space for sampler status */
	for(nj=0; nj<nStatusNum; nj++)
	{
		pSeqMtf->pStatusId[nj] = nj;
		FLEXSEQMTFCREATESTATUS(pSeqMtf, nj, nLen);
	}

	/* create space for sampler monitors */
	for(nj=0; nj<nMonitorNum; nj++)
	{
		pSeqMtf->pMonitorId[nj] = nj;
		FLEXSEQMTFCREATEMONITOR(pSeqMtf, nj, nLen);
	}

	/* Code sequence */
	FLEXSEQMTFCODENUCLEICSEQ(pSeqMtf, 0, nLen, pSeq->m_pSequence->m_pString);

	/* load conservation score */
	if(nUseCS == 1)
	{
		sprintf(strFilePath, "%s%s%d_%s.cs", strWorkPath, strCSPrefix, pSeq->m_nIndex, pSeq->m_strAlias);
		FLEXSEQMTFLOADCS(pSeqMtf, 0, nLen, strFilePath);
	}

	/* return */
	return pSeqMtf;
}

/* ----------------------------------------------------------------------- */ 
/* FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_GROUP: load sequence&scores from input */
/* files. Return a vector of pointers to newly created flexmotif/seq       */
/* complexes and nSeqCount, the number of loaded complexes.                */
/* ----------------------------------------------------------------------- */
struct FLEXSEQMOTIF **FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_GROUP(int *nSeqCount, 
					char strWorkPath[], char strSeqFile[], int nUseCS, char strCSPrefix[])
{
	/* define */
	int nCount = 0;
	struct tagSequence *pSeqList;
	struct tagSequence *pCurrentSeq;
	struct FLEXSEQMOTIF **vSeqMtf;
	int ni,nj,nLen;
	int nSeqNum = 1;
	int nScoreNum = 0;
	/* status[0]: Ai, true motif; status[1]: Si, site occupied */
	int nStatusNum = 0;
	/* monitor[0]: background likelihood; monitor[1,2...]: posterior of motifs */
	int nMonitorNum = 2;
	char strFilePath[MED_LINE_LENGTH];

	/* init */
	vSeqMtf = NULL;
	nMonitorNum = 2;
	if(nUseCS == 1)
		nScoreNum = 1;
	else
		nScoreNum = 0;

	/* load sequences from fasta file */
	sprintf(strFilePath, "%s%s", strWorkPath, strSeqFile);
	pSeqList = NULL;
	nCount = LoadFullSequenceList(strFilePath, &pSeqList);
	if(nCount <= 0)
		return NULL;

	/* create complexes */
	vSeqMtf = (struct FLEXSEQMOTIF **)calloc(nCount, sizeof(struct FLEXSEQMOTIF *));
	if(vSeqMtf == NULL)
		return NULL;
	
	/* code sequence */
	pCurrentSeq = pSeqList;
	for(ni=0; ni<nCount; ni++)
	{
		vSeqMtf[ni] = NULL;
		vSeqMtf[ni] = FLEXSEQMTFCREATE(ni, nSeqNum, nScoreNum, nStatusNum, nMonitorNum);
		if(vSeqMtf[ni] == NULL)
		{
			printf("Error: FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX, cannot create flex sequence/motif complex!\n");
			exit(EXIT_FAILURE);
		}
		nLen = pCurrentSeq->m_nLength;
		
		/* init length */
		strcpy(vSeqMtf[ni]->strAlias, pCurrentSeq->m_strAlias);
		vSeqMtf[ni]->nSiteNum = 0;
		vSeqMtf[ni]->nSeqLen = nLen;

		/* create space for sequence */
		for(nj=0; nj<nSeqNum; nj++)
		{
			vSeqMtf[ni]->pSeqId[nj] = nj;
			FLEXSEQMTFCREATESEQ(vSeqMtf[ni], nj, nLen);
		}

		/* create space for conservation score etc. */
		for(nj=0; nj<nScoreNum; nj++)
		{
			vSeqMtf[ni]->pScoreId[nj] = nj;
			FLEXSEQMTFCREATESCORE(vSeqMtf[ni], nj, nLen);
		}

		/* create space for sampler status */
		for(nj=0; nj<nStatusNum; nj++)
		{
			vSeqMtf[ni]->pStatusId[nj] = nj;
			FLEXSEQMTFCREATESTATUS(vSeqMtf[ni], nj, nLen);
		}

		/* create space for sampler monitors */
		for(nj=0; nj<nMonitorNum; nj++)
		{
			vSeqMtf[ni]->pMonitorId[nj] = nj;
			FLEXSEQMTFCREATEMONITOR(vSeqMtf[ni], nj, nLen);
		}

		/* create space for motif status */
		/* FLEXSEQMTFCREATEMOTIF(vSeqMtf[ni], nLen); */

		/* Code sequence */
		FLEXSEQMTFCODENUCLEICSEQ(vSeqMtf[ni], 0, nLen, pCurrentSeq->m_pSequence->m_pString);

		/* load conservation score */
		if(nUseCS == 1)
		{
			sprintf(strFilePath, "%s%s%d_%s.cs", strWorkPath, strCSPrefix, ni, pCurrentSeq->m_strAlias);
			FLEXSEQMTFLOADCS(vSeqMtf[ni], 0, nLen, strFilePath);
		}

		/* get next seq */
		pCurrentSeq = pCurrentSeq->m_pNext;
	}

	/* destroy fasta sequences */
	SequenceListClear(&pSeqList);

	/* return */
	*nSeqCount = nCount;
	return vSeqMtf;
}

/* ----------------------------------------------------------------------- */ 
/* FLEXSEQMTFLOADSEQ_FOR_HASHGENOME: load single sequence and              */
/* scores from input sequences. Return a pointer to newly created          */
/* flexmotif/seq complexes.                                                */
/* ----------------------------------------------------------------------- */
struct FLEXSEQMOTIF *FLEXSEQMTFLOADSEQ_FOR_HASHGENOME(int nSeqId,
				char strGenomePath[], char strChr[], int nStart, int nEnd,
				int nUseCS, char strCSPath[])
{
	/* define */
	struct FLEXSEQMOTIF *pSeqMtf;
	unsigned char *pBase;
	int nj,nLen;
	int nSeqNum = 1;
	int nScoreNum = 0;
	/* status[0]: Ai, true motif; status[1]: Si, site occupied */
	int nStatusNum = 0;
	/* monitor[0]: background likelihood; monitor[1,2...]: posterior of motifs */
	int nMonitorNum = 0;
	char strFilePath[MED_LINE_LENGTH];
	FILE *fpIn;
	/* for loading sequence */
	int nP1,nP2,nR1,nR2;
	int numread;
	int ni,nk;
	unsigned char bBase,bChar;

	/* init */
	if(nUseCS == 1)
		nScoreNum = 1;
	else
		nScoreNum = 0;

	/* load sequences */
	pSeqMtf = NULL;
	pSeqMtf = FLEXSEQMTFCREATE(nSeqId, nSeqNum, nScoreNum, nStatusNum, nMonitorNum);
	if(pSeqMtf == NULL)
	{
		printf("Error: FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_GENOME, cannot create flex sequence/motif complex!\n");
		exit(EXIT_FAILURE);
	}
	nLen = nEnd-nStart+1;
	if(nLen <= 0)
	{
		printf("Error: FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_GENOME, sequence length <= 0!\n");
		exit(EXIT_FAILURE);
	}

	/* init length */
	pSeqMtf->nSiteNum = 0;
	pSeqMtf->nSeqLen = nLen;

	/* create space for sequence */
	for(nj=0; nj<nSeqNum; nj++)
	{
		pSeqMtf->pSeqId[nj] = nj;
		FLEXSEQMTFCREATESEQ(pSeqMtf, nj, nLen);
	}

	/* create space for conservation score etc. */
	for(nj=0; nj<nScoreNum; nj++)
	{
		pSeqMtf->pScoreId[nj] = nj;
		FLEXSEQMTFCREATESCORE(pSeqMtf, nj, nLen);
	}

	/* create space for sampler status */
	for(nj=0; nj<nStatusNum; nj++)
	{
		pSeqMtf->pStatusId[nj] = nj;
		FLEXSEQMTFCREATESTATUS(pSeqMtf, nj, nLen);
	}

	/* create space for sampler monitors */
	for(nj=0; nj<nMonitorNum; nj++)
	{
		pSeqMtf->pMonitorId[nj] = nj;
		FLEXSEQMTFCREATEMONITOR(pSeqMtf, nj, nLen);
	}

	/* Load sequence */
	/* init */
	pBase = pSeqMtf->vSeq[0]->pMatElement;
	sprintf(strFilePath, "%s%s.sq", strGenomePath, strChr);
	fpIn = NULL;
	fpIn = fopen(strFilePath, "rb");
	if(fpIn == NULL)
	{
		printf("Warning: FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_GENOME, cannot open sequence file %s!\n", strFilePath);
		FLEXSEQMTFDESTROY(pSeqMtf);
		return NULL;
	}

	/* load seq */
	nP1 = (int)(nStart/2);
	nR1 = nStart%2;
	nP2 = (int)(nEnd/2);
	nR2 = nEnd%2;

	if( fseek( fpIn, nP1, SEEK_SET ) != 0 )
	{
		printf("Error: FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_GENOME, cannot locate the required sequence!\n");
		exit(EXIT_FAILURE);
	}

	nk = 0;

	/* first base */
	numread = fread( &bChar, sizeof( unsigned char ), 1, fpIn );
	if(nR1 == 0)
	{
		bBase = bChar >> 4;
		switch(bBase)
		{
			case 0: pBase[nk] = 0;
				break;
			case 1: pBase[nk] = 1;
				break;
			case 2: pBase[nk] = 2;
				break;
			case 3: pBase[nk] = 3;
				break;
			case 4: pBase[nk] = 10;
				break;
			case 5: pBase[nk] = 11;
				break;
			case 6: pBase[nk] = 12;
				break;
			case 7: pBase[nk] = 13;
				break;
			default: pBase[nk] = 4;
		}
		nk++;
	}

	if(nk < nLen)
	{
		bBase = bChar & 0x0F;
		switch(bBase)
		{
			case 0: pBase[nk] = 0;
				break;
			case 1: pBase[nk] = 1;
				break;
			case 2: pBase[nk] = 2;
				break;
			case 3: pBase[nk] = 3;
				break;
			case 4: pBase[nk] = 10;
				break;
			case 5: pBase[nk] = 11;
				break;
			case 6: pBase[nk] = 12;
				break;
			case 7: pBase[nk] = 13;
				break;
			default: pBase[nk] = 4;
		}
		nk++;
	}

	/* middle bases */
	if(nk < nLen)
	{
		for(ni=(nP1+1); ni<nP2; ni++)
		{
			numread = fread( &bChar, sizeof( unsigned char ), 1, fpIn );
			bBase = bChar >> 4;
			switch(bBase)
			{
				case 0: pBase[nk] = 0;
					break;
				case 1: pBase[nk] = 1;
					break;
				case 2: pBase[nk] = 2;
					break;
				case 3: pBase[nk] = 3;
					break;
				case 4: pBase[nk] = 10;
					break;
				case 5: pBase[nk] = 11;
					break;
				case 6: pBase[nk] = 12;
					break;
				case 7: pBase[nk] = 13;
					break;
				default: pBase[nk] = 4;
			}
			nk++;

			bBase = bChar & 0x0F;
			switch(bBase)
			{
				case 0: pBase[nk] = 0;
					break;
				case 1: pBase[nk] = 1;
					break;
				case 2: pBase[nk] = 2;
					break;
				case 3: pBase[nk] = 3;
					break;
				case 4: pBase[nk] = 10;
					break;
				case 5: pBase[nk] = 11;
					break;
				case 6: pBase[nk] = 12;
					break;
				case 7: pBase[nk] = 13;
					break;
				default: pBase[nk] = 4;
			}
			nk++;
		}
	}
	
	/* last base */
	if(nk < nLen)
	{
		numread = fread( &bChar, sizeof( unsigned char ), 1, fpIn );
		bBase = bChar >> 4;
		switch(bBase)
		{
			case 0: pBase[nk] = 0;
				break;
			case 1: pBase[nk] = 1;
				break;
			case 2: pBase[nk] = 2;
				break;
			case 3: pBase[nk] = 3;
				break;
			case 4: pBase[nk] = 10;
				break;
			case 5: pBase[nk] = 11;
				break;
			case 6: pBase[nk] = 12;
				break;
			case 7: pBase[nk] = 13;
				break;
			default: pBase[nk] = 4;
		}
		nk++;
	}

	if(nk < nLen)
	{
		if(nR2 == 1)
		{
			bBase = bChar & 0x0F;
			switch(bBase)
			{
				case 0: pBase[nk] = 0;
					break;
				case 1: pBase[nk] = 1;
					break;
				case 2: pBase[nk] = 2;
					break;
				case 3: pBase[nk] = 3;
					break;
				case 4: pBase[nk] = 10;
					break;
				case 5: pBase[nk] = 11;
					break;
				case 6: pBase[nk] = 12;
					break;
				case 7: pBase[nk] = 13;
					break;
				default: pBase[nk] = 4;
			}
			nk++;
		}
	}

	/* close file */
	fclose(fpIn);

	if(nk != nLen)
	{
		printf("Error: FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_GENOME, incorrect sequence length!\n");
		exit(EXIT_FAILURE);
	}


	/* load conservation score */
	if(nUseCS == 1)
	{
		sprintf(strFilePath, "%s%s.cs", strCSPath, strChr);
		fpIn = NULL;
		fpIn = fopen(strFilePath, "rb");
		if(fpIn == NULL)
		{
			printf("Error: FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_GENOME, cannot open conservation score file %s!\n", strFilePath);
			FLEXSEQMTFDESTROY(pSeqMtf);
			return NULL;
		}

		/* load score */
		if( fseek( fpIn, nStart, SEEK_SET ) != 0 )
		{
			printf("Error: FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_GENOME, cannot locate the required sequence!\n");
			exit(EXIT_FAILURE);
		}

		numread = fread((pSeqMtf->vScore[0])->pMatElement, sizeof(unsigned char), nLen, fpIn);
		if(numread != nLen)
		{
			printf("Error: FLEXSEQMTFLOADSEQ_FOR_SCANMATRIX_GENOME, loading error!\n");
			exit(EXIT_FAILURE);
		}

		fclose(fpIn);
	}

	/* return */
	return pSeqMtf;
}

/* ----------------------------------------------------------------------- */ 
/*      FLEXSEQMTFCODENUCLEICSEQ: code nucleic acid sequences.             */
/* ----------------------------------------------------------------------- */ 
int FLEXSEQMTFCODENUCLEICSEQ(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, int nLen, char strSeq[])
{
	int ni;
	unsigned char *pEle;

	/* check */
	if(pSeqMtf == NULL)
		return PROC_FAILURE;

	if( (nIndex >= pSeqMtf->nSeqNum) || (nIndex < 0) )
	{
		printf("Error: FLEXSEQMTFCODENUCLEICSEQ, nIndex out of range!\n");
		exit(EXIT_FAILURE);
	}

	if( ((int)(strlen(strSeq)) != nLen) || (pSeqMtf->vSeq[nIndex]->nWidth != nLen) )
	{
		printf("Error: FLEXSEQMTFCODENUCLEICSEQ, lengths not match!\n");
		exit(EXIT_FAILURE);
	}

	/* code */
	pEle = pSeqMtf->vSeq[nIndex]->pMatElement;
	for(ni=0; ni<nLen; ni++)
	{
		switch(strSeq[ni])
		{
			case 'a': *pEle = 10;
				break;
			case 'A': *pEle = 0;
				break;
			case 'c': *pEle = 11;
				break;
			case 'C': *pEle = 1;
				break;
			case 'g': *pEle = 12;
				break;
			case 'G': *pEle = 2;
				break;
			case 't': *pEle = 13;
				break;
			case 'T': *pEle = 3;
				break;
			default: *pEle = 4;
		}
		pEle++;
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*      FLEXSEQMTFLOADCS: load conservation score.                         */
/* ----------------------------------------------------------------------- */ 
int FLEXSEQMTFLOADCS(struct FLEXSEQMOTIF *pSeqMtf, int nIndex, int nSeqLen, char strConsFile[])
{
	/* define */
	FILE *fpIn;
	int numread;
	int nLen;

	/* init */
	nLen = nSeqLen;
	if(nLen <= 0)
		return PROC_FAILURE;

	fpIn = NULL;
	fpIn = fopen(strConsFile, "rb");
	if(fpIn == NULL)
	{
		printf("Error: FLEXSEQMTFLOADCS, cannot open conservation score file!\n");
		exit(EXIT_FAILURE);
	}

	/* load score */
	numread = fread((pSeqMtf->vScore[nIndex])->pMatElement, sizeof(unsigned char), nLen, fpIn);
	if(numread != nLen)
	{
		printf("Error: FLEXSEQMTFLOADCS, loading error!\n");
		exit(EXIT_FAILURE);
	}

	fclose(fpIn);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*    FLEXSEQMTFESTIMATENUCLEICBGMC: estimate markov background.           */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *FLEXSEQMTFESTIMATENUCLEICBGMC(struct FLEXSEQMOTIF **vSeqMtf, int nCount, 
								int nIndex, int nOrder, int nBaseTypeNum)
{
	/* define */
	struct DOUBLEMATRIX *pBG;
	int nHeight,nWidth,ni,nj,nk,nz,nLen;
	int nWordId,nBadLen;
	unsigned char *pEle;
	int nSLabel[BGMC_MAX_ORDER];
	double dTemp;
	double *pBaseCount,*pBaseStart;
	double dTotalCount;
	int nScale;

	/* check */
	if((vSeqMtf == NULL) || (nCount <= 0) || (nOrder < 0) || (nIndex < 0) || (nOrder > BGMC_MAX_ORDER))
	{
		printf("Error: FLEXSEQMTFESTIMATENUCLEICBGMC, parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	
	/* init */
	pBG = NULL;
	nWidth = nBaseTypeNum;
	nHeight = (int)(pow((double)nWidth, (double)nOrder));
	nScale = (int)(pow((double)nWidth, (double)(nOrder-1)));
	pBG = CreateDoubleMatrix(nHeight, nWidth);
	if(pBG == NULL)
	{
		printf("Error: FLEXSEQMTFESTIMATENUCLEICBGMC, can't create background matrix!\n");
		exit(EXIT_FAILURE);
	}
	pBaseStart = pBG->pMatElement;
	for(ni=0; ni<nHeight; ni++)
	{
		for(nj=0; nj<nWidth; nj++)
		{
			*pBaseStart = 1.0;
			pBaseStart++;
		}
	}

	/* create */
	for(ni=0; ni<nCount; ni++)
	{
		if(vSeqMtf[ni] == NULL)
			continue;
		if(nIndex >= vSeqMtf[ni]->nSeqNum)
			continue;
		if(vSeqMtf[ni]->vSeq[nIndex] == NULL)
			continue;

		for(nj=0; nj<BGMC_MAX_ORDER; nj++)
			nSLabel[nj] = 0;
		nLen = vSeqMtf[ni]->vSeq[nIndex]->nWidth;
		
		if(nOrder == 0)
		{
			pEle = vSeqMtf[ni]->vSeq[nIndex]->pMatElement;
			for(nj=0; nj<nLen; nj++)
			{
				nWordId = (int)(*pEle);
				if(nWordId <nWidth)
				{
					dTemp = DMGETAT(pBG, 0, nWordId)+1.0;
					DMSETAT(pBG, 0, nWordId, dTemp);
				}
				pEle++;
			}
		}
		else
		{
			nWordId = 0;
			nBadLen = 0;
			pEle = vSeqMtf[ni]->vSeq[nIndex]->pMatElement;
			for(nj=0; nj<nOrder; nj++)
			{
				nSLabel[nj] = (int)(*pEle);
				nk = (int)(*pEle);
				if(nk < nWidth)
				{
					nWordId = nWidth*nWordId+nk;
				}
				else
				{
					nWordId = nWidth*nWordId;
					nBadLen++;
				}
				
				pEle++;
			}
			for(; nj<nLen; nj++)
			{
				nk = (int)(*pEle);
				if( (nk < nWidth) && (nBadLen == 0) )
				{
					dTemp = DMGETAT(pBG, nWordId, nk)+1.0;
					DMSETAT(pBG, nWordId, nk, dTemp);
				}

				if(nSLabel[0] < nWidth)
				{
					nWordId -= nSLabel[0]*nScale;
				}
				else
				{
					nBadLen--;
				}
				for(nz=0; nz<(nOrder-1); nz++)
				{
					nSLabel[nz] = nSLabel[nz+1];
				}
				nSLabel[nz] = nk;
				
				if(nk<nWidth)
				{
					nWordId = nWordId*nWidth+nk;
				}
				else
				{
					nWordId = nWordId*nWidth;
					nBadLen++;
				}

				/* get next */
				pEle++;
			}
		}
	}

	/* normalize */
	pBaseStart = pBG->pMatElement;
	for(ni=0; ni<nHeight; ni++)
	{
		pBaseCount = pBaseStart;
		dTotalCount = 0.0;
		for(nj=0; nj<nWidth; nj++)
		{
			dTotalCount += (*pBaseCount);
			pBaseCount++;
		}

		if(dTotalCount>0.0)
		{
			for(nj=0; nj<nWidth; nj++)
			{
				*pBaseStart = *pBaseStart/dTotalCount;
				pBaseStart++;
			}
		}
		else
		{
			for(nj=0; nj<nWidth; nj++)
			{
				*pBaseStart = 1.0/(double)nWidth;
				pBaseStart++;
			}
		}
	}

	/* return */
	return pBG;
}

/* ----------------------------------------------------------------------- */ 
/*      FLEXMTFMCREATE: create PWM of flex motifs                          */
/* ----------------------------------------------------------------------- */ 
struct FLEXMOTIFMATRIX *FLEXMTFMCREATE()
{
	/* new motif matrix */
	struct FLEXMOTIFMATRIX *pM;

	/* create motif matrix */
	pM = NULL;
	pM = (struct FLEXMOTIFMATRIX *)calloc(1, sizeof(struct FLEXMOTIFMATRIX));
	if (pM == NULL)
	{
		printf("Error: FLEXMOTIFMATRIX, cannot create flexmotif matrix.\n");
		exit(EXIT_FAILURE);
	}
	
	/* init */
	pM->dCm = 0.0;
	pM->dCp = 0.0;
	pM->dPm = 0.0;
	pM->dPp = 0.0;
	pM->pPriorCount = NULL;
	pM->pSampleCount = NULL;
	pM->pPWM = NULL;

	pM->nHeadExtraNum = 0;
	pM->nTailExtraNum = 0;
	pM->pPriorHeadExtra = NULL;
	pM->pPriorTailExtra = NULL;

	/* return */
	return pM;
}

/* ----------------------------------------------------------------------- */ 
/*      FLEXMTFMDESTROY: destroy PWM of flex motifs                        */
/* ----------------------------------------------------------------------- */ 
int FLEXMTFMDESTROY(struct FLEXMOTIFMATRIX *pM)
{
	if(pM == NULL)
		return PROC_SUCCESS;

	/* destroy */
	if(pM->pPriorCount != NULL)
	{
		DestroyDoubleMatrix(pM->pPriorCount);
		pM->pPriorCount = NULL;
	}
	if(pM->pSampleCount != NULL)
	{
		DestroyDoubleMatrix(pM->pSampleCount);
		pM->pSampleCount = NULL;
	}
	if(pM->pPWM != NULL)
	{
		DestroyDoubleMatrix(pM->pPWM);
		pM->pPWM = NULL;
	}
	if(pM->pPriorHeadExtra != NULL)
	{
		DestroyDoubleMatrix(pM->pPriorHeadExtra);
		pM->pPriorHeadExtra = NULL;
	}
	if(pM->pPriorTailExtra != NULL)
	{
		DestroyDoubleMatrix(pM->pPriorTailExtra);
		pM->pPriorTailExtra = NULL;
	}

	free(pM);
	pM = NULL;

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  FLEXMTFMREFRESH: update PWM based on current count matrix.             */
/* ----------------------------------------------------------------------- */ 
int FLEXMTFMREFRESH(struct FLEXMOTIFMATRIX *pM)
{
	/* define */
	double dTotal;
	double *pEle1,*pEle2,*pEle3,*pEle4;
	int ni,nj;

	/* check */
	if(pM == NULL)
		return PROC_FAILURE;

	/* update PWM */
	pEle1 = pM->pPriorCount->pMatElement;
	pEle2 = pM->pSampleCount->pMatElement;
	pEle3 = pM->pPWM->pMatElement;

	for(ni=0; ni<pM->pPriorCount->nHeight; ni++)
	{
		dTotal = 0.0;
		pEle4 = pEle3;
		for(nj=0; nj<pM->pPriorCount->nWidth; nj++)
		{
			*pEle4 = (*pEle1)+(*pEle2);
			dTotal += (*pEle4);
			pEle1++;
			pEle2++;
			pEle4++;
		}

		for(nj=0; nj<pM->pPriorCount->nWidth; nj++)
		{
			*pEle3 = log(*pEle3/dTotal);
			/* if(*pEle3 < -1e10)
			{
				*pEle3 = *pEle3;
			} */
			pEle3++;
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  FLEXMTFMSCORE: get motif score.                                        */
/* ----------------------------------------------------------------------- */ 
double FLEXMTFMSCORE(struct FLEXMOTIFMATRIX *pM, struct DOUBLEMATRIX *pBG)
{
	/* define */
	double dScore,dTemp1,dTemp2,dStemp,dTotal;
	int ni,nj;

	/* check */
	if((pM == NULL) || (pBG == NULL))
		return 0.0;

	/* init */
	dScore = 0.0;

	/* get score */
	for(ni=0; ni<pM->pPriorCount->nHeight; ni++)
	{
		dStemp = 0.0;
		dTotal = 0.0;
		for(nj=0; nj<pM->pPriorCount->nWidth; nj++)
		{
			dTotal += DMGETAT(pM->pPriorCount, ni, nj)+DMGETAT(pM->pSampleCount, ni, nj);
		}
		for(nj=0; nj<pM->pPriorCount->nWidth; nj++)
		{
			dTemp1 = (DMGETAT(pM->pPriorCount, ni, nj)+DMGETAT(pM->pSampleCount, ni, nj))/dTotal;
			dTemp2 = DMGETAT(pBG, 0, nj);
			dStemp += dTemp1*log(dTemp1/dTemp2);
		}

		dScore += log(dTotal)*dStemp;
	}

	dScore /= pM->pPriorCount->nHeight;

	/* return */
	return dScore;
}

/* ----------------------------------------------------------------------- */ 
/*  FLEXMTFMSHIFTPRIORCOUNT: shift prior count.                            */
/* ----------------------------------------------------------------------- */
int FLEXMTFMSHIFTPRIORCOUNT(struct FLEXMOTIFMATRIX *pM, int nMotifLen, 
							int nOffset, double dDefaultPrior)
{
	/* define */
	struct DOUBLEMATRIX *pNewPriorCount;
	int nDHead,nEHead,nBS,nBE,nETail,nDTail;
	double *pEle1,*pEle2;
	int ni,nj;
	int nBaseTypeNum;

	/* init */
	nBaseTypeNum = pM->pPriorCount->nWidth;
	pNewPriorCount = NULL;
	pNewPriorCount = CreateDoubleMatrix(nMotifLen, nBaseTypeNum);
	if(pNewPriorCount == NULL)
	{
		printf("Error: FLEXMTFMSHIFTPRIORCOUNT, cannot create new prior count matrix!\n");
		exit(EXIT_FAILURE);
	}

	/* add head */
	if(nOffset < 0)
	{
		if( pM->nHeadExtraNum <= (-nOffset) )
		{
			if(pM->nHeadExtraNum >= 0)
			{
				nEHead = pM->nHeadExtraNum;
				nDHead = (-nOffset)-pM->nHeadExtraNum;
			}
			else
			{
				nEHead = 0;
				nDHead = -nOffset;
			}
		}
		else
		{
			nEHead = (-nOffset);
			nDHead = 0;
		}
	}
	else
	{
		nEHead = 0;
		nDHead = 0;
	}
	
	if((nMotifLen-nEHead-nDHead) > 0)
	{
		if(nOffset <= 0)
		{
			nBS = 0;
			nBE = pM->pPriorCount->nHeight;
		}
		else
		{
			nBS = nOffset;
			nBE = pM->pPriorCount->nHeight;
		}
	}
	else
	{
		nBS = 0;
		nBE = 0;
		nETail = 0;
		nDTail = 0;
	}

	if( (nBE-nBS+nEHead+nDHead) >= nMotifLen )
	{
		nBE = nMotifLen-nEHead-nDHead+nBS;
		nETail = 0;
		nDTail = 0;
	}
	else if( (nBE-nBS+nEHead+nDHead) < nMotifLen )
	{
		nETail = nMotifLen-(nBE-nBS+nEHead+nDHead);
		if( nETail > pM->nTailExtraNum )
		{
			if(pM->nTailExtraNum >= 0)
			{
				nDTail = nETail-pM->nTailExtraNum;
				nETail = pM->nTailExtraNum;
			}
			else
			{
				nDTail = nETail;
				nETail = 0;
			}
		}
		else
		{
			nDTail = 0;
		}
	}

	/* get shifted prior matrix */
	pEle1 = pNewPriorCount->pMatElement;
	for(ni=0; ni<nDHead; ni++)
	{
		for(nj=0; nj<nBaseTypeNum; nj++)
		{
			*pEle1 = dDefaultPrior;
			pEle1++;
		}
	}

	if(nEHead > 0)
	{
		pEle2 = pM->pPriorHeadExtra->pMatElement+(pM->nHeadExtraNum-nEHead)*nBaseTypeNum;
	}
	for(ni=0; ni<nEHead; ni++)
	{
		for(nj=0; nj<nBaseTypeNum; nj++)
		{
			*pEle1 = *pEle2;
			pEle1++;
			pEle2++;
		}
	}

	pEle2 = pM->pPriorCount->pMatElement+nBS*nBaseTypeNum;
	for(ni=nBS; ni<nBE; ni++)
	{
		for(nj=0; nj<nBaseTypeNum; nj++)
		{
			*pEle1 = *pEle2;
			pEle1++;
			pEle2++;
		}
	}

	for(ni=0; ni<nETail; ni++)
	{
		pEle2 = pM->pPriorTailExtra->pMatElement+(pM->nTailExtraNum-ni-1)*nBaseTypeNum;
		for(nj=0; nj<nBaseTypeNum; nj++)
		{
			*pEle1 = *pEle2;
			pEle1++;
			pEle2++;
		}
	}

	for(ni=0; ni<nDTail; ni++)
	{
		for(nj=0; nj<nBaseTypeNum; nj++)
		{
			*pEle1 = dDefaultPrior;
			pEle1++;
		}
	}

	/* update */
	DestroyDoubleMatrix(pM->pPriorCount);
	pM->pPriorCount = pNewPriorCount;

	/* return */
	return PROC_SUCCESS;
}


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
					int nUseCS, double dC, char strCSPath[])
{
	/* background type: 0=region; 1=genome */
	int nBGType = 0;

	/* sequences and conservation */
	/* chromosome length */
	struct INTMATRIX *pChrLen = NULL;
	int nChr;
	int nCurrentChr, nCurrentPos;
	int nNewSampLine;
	int nSampEnd;
	
	/* all sequences */
	struct FLEXSEQMOTIF *pSeqMtf = NULL;
	/* sequence number */
	int nSeqCount = 0;
	/* site number */
	int nEffecLen = 0;
	int nConsEffecLen = 0;
	int nSiteNum = 0;
	/* total sequence length and total site number */
	double dTotLen = 0.0;
	double dTotConsLen = 0.0;
	double dTotSite = 0.0;
	
	/* background */
	struct DOUBLEMATRIX *pBGF = NULL;
	struct DOUBLEMATRIX *pBGB = NULL;
	/* log background */
	struct DOUBLEMATRIX *pLogBGF = NULL;
	struct DOUBLEMATRIX *pLogBGB = NULL;
	/* background0 */
	struct DOUBLEMATRIX *pBG0 = NULL;
	/* log background0 */
	struct DOUBLEMATRIX *pLogBG0 = NULL;
	
	/* motif matrix */
	struct DOUBLEMATRIX *pMotif = NULL;
	
	/* number of base types */
	int nBaseTypeNum = 4;
	int nScale = 0;
	int nBGHeight;

	/* motif consensus */
	FILE *fpIn;
	FILE *fpSamp;
	FILE *fpOut;
	FILE *fpOut2;
	char strLine[LONG_LINE_LENGTH];
	char strSeqAlias[LINE_LENGTH];

	/* coordinates */
	int nMotifLen;
	char strChr[LINE_LENGTH];
	int nStart,nEnd;
	int nM1,nM2;
	int nActualStart,nActualEnd;
	int nFrom,nTo;
	
	/* other variables */
	int ni,nj;
	double *pEle1,*pEle2;
	double dSum;
	
	/* #################################### */
	/* initialize                           */
	/* #################################### */
	nBaseTypeNum = 4;
	nScale = (int)(pow((double)nBaseTypeNum, (double)(nBGOrder-1)));
	nBGHeight = (int)(pow((double)nBaseTypeNum, (double)nBGOrder));
	
	if(strcmp(strBGType, "GENOME") == 0)
	{
		nBGType = 1;
	}
	else
	{
		nBGType = 0;
	}
	
	/* #################################### */
	/* load initial parameters              */
	/* #################################### */
	AdjustDirectoryPath(strGenomePath);
	sprintf(strLine, "%schrlen.txt", strGenomePath);
	pChrLen = IMLOAD(strLine);
	if(pChrLen == NULL)
	{
		printf("Error: MotifMap_ScanMatrix_Integrate_Genome_Main, cannot load chromosome length!\n");
		exit(EXIT_FAILURE);
	}
	if(nUseCS == 1)
	{
		AdjustDirectoryPath(strCSPath);
	}
	if(nBGType == 1)
	{
		AdjustDirectoryPath(strBGPath);
	}

	/* #################################### */
	/* load motif                           */
	/* #################################### */
	pMotif = NULL;
	pMotif = DMLOAD(strMotifPath);
	if( pMotif == NULL)
	{
		printf("Error: MotifMap_ScanMatrix_Integrate_Genome_Main, cannot load init motif pseudocount!\n");
		exit(EXIT_FAILURE);
	}

	nMotifLen = pMotif->nHeight;

	pEle1 = pMotif->pMatElement;
	for(ni=0; ni<pMotif->nHeight; ni++)
	{
		pEle2 = pEle1;
		dSum = 0.0;
		for(nj=0; nj<pMotif->nWidth; nj++)
		{
			*pEle2 += 1e-3;
			dSum += (*pEle2);
			pEle2++;
		}

		for(nj=0; nj<pMotif->nWidth; nj++)
		{
			*pEle1 = log(*pEle1/dSum);
			pEle1++;
		}
	}


	/* #################################### */
	/* compute background if needed         */
	/* #################################### */
	
	/* genomic local background */
	if(nBGType == 1)
	{
	}
	/* region based background */
	else
	{
		pBGF = NULL;
		pBGF = CreateDoubleMatrix(nBGHeight, nBaseTypeNum);
		pBGB = NULL;
		pBGB = CreateDoubleMatrix(nBGHeight, nBaseTypeNum);
		pBG0 = NULL;
		pBG0 = CreateDoubleMatrix(1, nBaseTypeNum);
		for(ni=0; ni<nBGHeight; ni++)
		{
			for(nj=0; nj<nBaseTypeNum; nj++)
			{
				DMSETAT(pBGF, ni, nj, 1.0);
				DMSETAT(pBGB, ni, nj, 1.0);
			}
		}
		for(nj=0; nj<nBaseTypeNum; nj++)
		{
			DMSETAT(pBG0, 0, nj, 1.0);
		}
		
		if( (pBGF == NULL) || (pBGB == NULL) || (pBG0 == NULL) )
		{
			printf("Error: MotifMap_ScanMatrix_Integrate_Genome_Main, failure to create background markov matrix!\n");
			exit(EXIT_FAILURE);
		}

		/* compute background */
		MotifMap_ScanMatrix_Genome_ComputeBackground_Region(strGenomePath,
					strCodPath, nBGOrder, pBGF, pBGB, pBG0);
		
		pLogBGF = NULL;
		pLogBGF = DMCLONE(pBGF);
		DMLOGTS(pLogBGF);

		pLogBGB = NULL;
		pLogBGB = DMCLONE(pBGB);
		DMLOGTS(pLogBGB);

		pLogBG0 = NULL;
		pLogBG0 = DMCLONE(pBG0);
		DMLOGTS(pLogBG0);

		sprintf(strLine, "%s.bg0", strOutputPath);
		DMSAVE(pBG0, strLine);
		sprintf(strLine, "%s.bgf", strOutputPath);
		DMSAVE(pBGF, strLine);
		sprintf(strLine, "%s.bgb", strOutputPath);
		DMSAVE(pBGB, strLine);
	}

	/* #################################### */
	/* load sequences and conservation      */
	/* #################################### */
	fpIn = NULL;
	fpIn = fopen(strCodPath, "r");
	if(fpIn == NULL)
	{
        printf("Error: MotifMap_ScanMatrix_Integrate_Genome_Main, cannot open coordinates file!\n");
		exit(EXIT_FAILURE);
	}
	
	fpSamp = NULL;
	fpSamp = fopen(strSampPath, "r");
	if(fpSamp == NULL)
	{
        printf("Error: MotifMap_ScanMatrix_Integrate_Genome_Main, cannot open coordinates file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
        printf("Error: MotifMap_ScanMatrix_Integrate_Genome_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	
	sprintf(strLine, "%s.stat", strOutputPath);
	fpOut2 = NULL;
	fpOut2 = fopen(strLine, "w");
	if(fpOut2 == NULL)
	{
        printf("Error: MotifMap_ScanMatrix_Integrate_Genome_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	/* read */
	nSampEnd = 0;
	nNewSampLine = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpSamp) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] != '\0')
		{
			sscanf(strLine, "%d %d", &nCurrentChr, &nCurrentPos);
			nNewSampLine = 1;
			break;
		}
	}
	if(nNewSampLine == 0)
	{
		nSampEnd = 1;
		printf("MotifMap_ScanMatrix_Integrate_Genome_Main, end of sampling file!\n");
	
		/* close files */
		fclose(fpIn);
		fclose(fpSamp);
		fclose(fpOut);
		fclose(fpOut2);

		/* #################################### */
		/* release memory                       */
		/* #################################### */
		DestroyDoubleMatrix(pMotif);
		DestroyIntMatrix(pChrLen);
		if(nBGType == 0)
		{
			DestroyDoubleMatrix(pBGF);
			DestroyDoubleMatrix(pBGB);
			DestroyDoubleMatrix(pLogBGF);
			DestroyDoubleMatrix(pLogBGB);
			DestroyDoubleMatrix(pBG0);
			DestroyDoubleMatrix(pLogBG0);
		}
	}
	
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		sscanf(strLine, "%s %s %d %d", strSeqAlias, strChr,
			&nStart, &nEnd);
        if(nStart > nEnd)
		{
			printf("Warning: MotifMap_ScanMatrix_Integrate_Genome_Main, start > end: %s!\n", strLine);
			continue;
		}
		nChr = Genome_ChromosomeName_To_Index(strChr, strSpecies);
		if(nStart < 0)
			nStart = 0;
		if(nEnd < 0)
			nEnd = 0;
		if(nStart >= pChrLen->pMatElement[nChr-1])
			nStart = pChrLen->pMatElement[nChr-1]-1;
		if(nEnd >= pChrLen->pMatElement[nChr-1])
			nEnd = pChrLen->pMatElement[nChr-1]-1;

		/* process segment by segment */
		nM1 = nStart;
		while(nM1 <= nEnd)
		{
			nM2 = nM1+GENOME_CONTIG_LEN-1;
			if(nM2 > nEnd)
				nM2 = nEnd;

			nActualStart = nM1-nW-nBGOrder;
			if(nActualStart < 0)
			{
				nActualStart = 0;
			}
			nFrom = nM1-nActualStart;
			nTo = nM2-nActualStart;
			nActualEnd = nM2+nW+nBGOrder+nMotifLen-1;
			if(nActualEnd >= pChrLen->pMatElement[nChr-1])
			{
				nActualEnd = pChrLen->pMatElement[nChr-1]-1;
			}

			/* #################################### */ 
			/* load sequence and score              */
			/* #################################### */
			pSeqMtf = NULL;
			pSeqMtf = FLEXSEQMTFLOADSEQ_FOR_INTEGRATEMATRIX_GENOME(nSeqCount,
				strGenomePath, strChr, nActualStart, nActualEnd, 
				nUseCS, strCSPath);
			if(pSeqMtf == NULL)
			{
				break;
			}
			if(nRepeat == 1)
			{
				FLEXSEQMTF_REPEATUNMASK(pSeqMtf, 0);
			}

			/* #################################### */
			/* set background                       */
			/* #################################### */
			/* genomic local background */
			if(nBGType == 1)
			{
				MotifMap_ScanMatrix_InitBGLogLike_Genome(pSeqMtf, 0,
					strChr, nActualStart, nActualEnd, 
					nBGOrder, strBGPath, nBGStepSize);
			}
			/* region based background */
			else
			{
				MotifMap_ScanMatrix_InitBGLogLike(pSeqMtf, 0,
					nBGOrder, pLogBGF, pLogBGB, pLogBG0);
			}
			
			/* #################################### */
			/* scan motif                           */
			/* #################################### */
			MotifMap_IntegrateMatrix_In_FlexSeqMtf(pSeqMtf, 0,
				nW, pMotif, dR, 
				nUseCS, dC, &nEffecLen, &nConsEffecLen, &nSiteNum);
			dTotLen += nEffecLen;
			dTotConsLen += nConsEffecLen;
			dTotSite += nSiteNum;

			/* #################################### */
			/* export results                       */
			/* #################################### */
			MotifMap_IntegrateMatrix_Export_Genome(pSeqMtf, 0,
				nChr, nActualStart, nFrom, nTo, 
				&nCurrentChr, &nCurrentPos, fpSamp, fpOut, &nSampEnd);

			/* #################################### */
			/* release memory                       */
			/* #################################### */
			FLEXSEQMTFDESTROY(pSeqMtf);
			nM1 = nM2+1;

			if(nSampEnd == 1)
			{
				break;
			}
		}
		nSeqCount++;
		if(nSampEnd == 1)
		{
			break;
		}
	}


	/* save statistics */
	fprintf(fpOut2, "EffecLen= %d\n", (int)dTotLen);
	fprintf(fpOut2, "TotalSite= %d\n", (int)dTotSite);
	fprintf(fpOut2, "ConsLen= %d\n", (int)dTotConsLen);
	
	/* close files */
	fclose(fpIn);
	fclose(fpSamp);
	fclose(fpOut);
	fclose(fpOut2);

	/* #################################### */
	/* release memory                       */
	/* #################################### */
	DestroyDoubleMatrix(pMotif);
	DestroyIntMatrix(pChrLen);
	if(nBGType == 0)
	{
		DestroyDoubleMatrix(pBGF);
		DestroyDoubleMatrix(pBGB);
		DestroyDoubleMatrix(pLogBGF);
		DestroyDoubleMatrix(pLogBGB);
		DestroyDoubleMatrix(pBG0);
		DestroyDoubleMatrix(pLogBG0);
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_IntegrateMatrix_In_FlexSeqMtf: integrate motif signals.       */
/* ----------------------------------------------------------------------- */ 
int MotifMap_IntegrateMatrix_In_FlexSeqMtf(struct FLEXSEQMOTIF *pSeqMtf, 
			int nIndex, int nW,
			struct DOUBLEMATRIX *pLogPWM, double dR, 
			int nUseCS, double dC, 
			int *pEffecLen, int *pConsEffecLen, int *pSiteNum)
{
	/* define */
	unsigned char *pBase;
	unsigned char *pCS;
	double *pLF,*pLB,*pLW;
	int ni,nj,nLen;
	int nBaseTypeNum;
	struct FLEXMOTIFSITE *pPrev;
	double dLike[2];
	int nMask[2];
	int nMotifLen;
	int nActualTo;
	double dTotConserve;
	double dAveConserve;
	double dTotal;
	int n2W = 2*nW+1;
	
	/* init */
	nBaseTypeNum = pLogPWM->nWidth;
	nMotifLen = pLogPWM->nHeight;
	nLen = pSeqMtf->vSeq[nIndex]->nWidth;

	if(nLen<nMotifLen)
	{
		printf("Warning: MotifMap_IntegrateMatrix_In_FlexSeqMtf, seqence %d length < motif length!\n", pSeqMtf->nId);
		return PROC_SUCCESS;
	}
	nActualTo = nLen-nMotifLen;
	
	*pEffecLen = 0;
	*pConsEffecLen = 0;
	*pSiteNum = 0;
	
	pBase = pSeqMtf->vSeq[nIndex]->pMatElement;
	pLF = pSeqMtf->vMonitor[2]->pMatElement;
	pLB = pSeqMtf->vMonitor[3]->pMatElement;
	pLW = pSeqMtf->vMonitor[4]->pMatElement;

	pPrev = NULL;
	
	dTotConserve = 0.0;
	if(nUseCS == 1)
	{
		pCS = pSeqMtf->vScore[0]->pMatElement;
		for(nj=0; nj<(nMotifLen-1); nj++)
		{
			dTotConserve += pCS[nj];
		}
	}

	for(ni=0; ni<=nActualTo; ni++)
	{
		if(nUseCS == 1)
		{
			dTotConserve += pCS[ni+nMotifLen-1];
		}

		/* if repeat or 'N', skip it */
		if(pBase[ni] >= nBaseTypeNum)
		{
			if(nUseCS == 1)
			{
				dTotConserve -= pCS[ni];
			}
			continue;
		}

		/* add effect length */
		*pEffecLen += 1;
		
		/* compute conservation */
		if(nUseCS == 1)
		{
			dAveConserve = dTotConserve/(double)nMotifLen;
			if(dAveConserve < dC)
			{
				dTotConserve -= pCS[ni];
				continue;
			}
			else
			{
				*pConsEffecLen += 1;
			}
		}
		else
		{
			*pConsEffecLen += 1;
		}


		/* get motif likelihood and set posterior probability */
		dLike[0] = MotifMap_ScanMatrix_MotifLikeRatio(pSeqMtf, nIndex, ni, 
					pLogPWM, '+', (nMask+0));
		dLike[1] = MotifMap_ScanMatrix_MotifLikeRatio(pSeqMtf, nIndex, ni, 
					pLogPWM, '-', (nMask+1));
		
		if(nMask[0] == 0)
		{
			dLike[0] = exp(dLike[0]);
			if(dLike[0] >= dR)
				pLF[ni] = dLike[0];
		}

		if(nMask[1] == 0)
		{
			dLike[1] = exp(dLike[1]);
			if(dLike[1] >= dR)
				pLB[ni] = dLike[1];
		}

		if(nUseCS == 1)
		{
			dTotConserve -= pCS[ni];
		}
	}


	dTotal = 0.0;
	for(ni=0; ni<nW; ni++)
	{
		if(ni >= nLen)
			break;
		dTotal = dTotal+pLF[ni]+pLB[ni];
	}
	for(ni=0; ni<=nW; ni++)
	{
		if(ni >= nLen)
			break;

		nj = ni+nW;
		if(nj<nLen)
			dTotal = dTotal+pLF[nj]+pLB[nj];
		
		pLW[ni] = dTotal;
	}
	for(; ni<nLen; ni++)
	{
		nj = ni-nW-1;
		if(nj>=0)
			dTotal = dTotal-pLF[nj]-pLB[nj];
		nj = ni+nW;
		if(nj<nLen)
			dTotal = dTotal+pLF[nj]+pLB[nj];
		pLW[ni] = dTotal;
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MotifMap_IntegrateMatrix_Export_Genome: export integrated singal.      */
/* ----------------------------------------------------------------------- */ 
int MotifMap_IntegrateMatrix_Export_Genome(struct FLEXSEQMOTIF *pSeqMtf, int nIndex,
				int nChr, int nActualStart, int nFrom, int nTo, 
				int *pnCurrentChr, int *pnCurrentPos, FILE *fpSamp, FILE *fpOut,
				int *pnSampEnd)
{
	/* define */
	char strLine[LONG_LINE_LENGTH];
	int nNewSampLine;
	int nj;

	/* output */
	if(nChr < (*pnCurrentChr))
	{
		return PROC_SUCCESS;
	}

	if( (nChr == (*pnCurrentChr)) && ((nActualStart+nTo) < (*pnCurrentPos)) )
	{
		return PROC_SUCCESS;
	}

	while(nChr > (*pnCurrentChr))
	{
		nNewSampLine = 0;
		while(fgets(strLine, LONG_LINE_LENGTH, fpSamp) != NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] != '\0')
			{
				sscanf(strLine, "%d %d", pnCurrentChr, pnCurrentPos);
				nNewSampLine = 1;
				break;
			}
		}

		if(nNewSampLine == 0)
		{
			*pnSampEnd = 1;
			return PROC_SUCCESS;
		}
	}

	if(nChr < (*pnCurrentChr))
	{
		return PROC_SUCCESS;
	}

	while( (nChr == (*pnCurrentChr)) && ((nActualStart+nFrom) > (*pnCurrentPos)) )
	{
		nNewSampLine = 0;
		while(fgets(strLine, LONG_LINE_LENGTH, fpSamp) != NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] != '\0')
			{
				sscanf(strLine, "%d %d", pnCurrentChr, pnCurrentPos);
				nNewSampLine = 1;
				break;
			}
		}

		if(nNewSampLine == 0)
		{
			*pnSampEnd = 1;
			return PROC_SUCCESS;
		}
	}

	if(nChr < (*pnCurrentChr))
	{
		return PROC_SUCCESS;
	}

	if((nActualStart+nTo) < (*pnCurrentPos))
	{
		return PROC_SUCCESS;
	}

	while( (nChr == (*pnCurrentChr)) && ((*pnCurrentPos) >= (nActualStart+nFrom)) 
		&& ((*pnCurrentPos) <= (nActualStart+nTo)) )
	{
		nj = *pnCurrentPos-nActualStart;
		fprintf(fpOut, "%d\t%d\t%9.7e\n", *pnCurrentChr, *pnCurrentPos, 
			pSeqMtf->vMonitor[4]->pMatElement[nj]);

		nNewSampLine = 0;
		while(fgets(strLine, LONG_LINE_LENGTH, fpSamp) != NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] != '\0')
			{
				sscanf(strLine, "%d %d", pnCurrentChr, pnCurrentPos);
				nNewSampLine = 1;
				break;
			}
		}

		if(nNewSampLine == 0)
		{
			*pnSampEnd = 1;
			return PROC_SUCCESS;
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  FLEXSEQMTF_REPEATUNMASK: unmask repeat sequences.                      */
/* ----------------------------------------------------------------------- */ 
int FLEXSEQMTF_REPEATUNMASK(struct FLEXSEQMOTIF *pSeqMtf, int nIndex)
{
	/* define */
	unsigned char *vSeq;
	int ni;

	/* init */
	if(pSeqMtf == NULL)
		return PROC_SUCCESS;

	vSeq = pSeqMtf->vSeq[nIndex]->pMatElement;
	for(ni=0; ni<pSeqMtf->nSeqLen; ni++)
	{
		if( (vSeq[ni] >= 10) && (vSeq[ni] <= 13) )
		{
			vSeq[ni] -= 10;
		}
	}

	/* return */
	return PROC_SUCCESS;
}