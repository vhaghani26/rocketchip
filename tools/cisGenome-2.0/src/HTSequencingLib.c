/* ----------------------------------------------------------------------- */
/*  HTSequencingLib.c : implementation of the high throughput sequencing   */
/*  library                                                                */
/*  Author : Ji HongKai ; Time: 2007.10                                    */
/* ----------------------------------------------------------------------- */

#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
#include "limits.h"

#include "MathLib.h"
#include "MatrixLib.h"
#include "RandomLib.h"
#include "StringLib.h"
#include "GenomeLib.h"
#include "AffyLib.h"
#include "TilingArrayLib.h"
#include "HTSequencingLib.h"

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Unique()                                                       */
/*  Take a sorted aln file as input, remove redundant reads.               */
/*  At each position, keep at most nMax reads.                             */
/* ----------------------------------------------------------------------- */ 
int HTS_Aln2Unique(char strInputPath[], char strOutputPath[], int nMax)
{	
	/* define */
	int ni;
	char strLastChr[MED_LINE_LENGTH];
	char strChr[MED_LINE_LENGTH];
	int nLastPos,nPos;
	int nStrand;
	int *vLastStrand;
	int nLastC;
	char strLine[MED_LINE_LENGTH];

	FILE *fpIn;
	FILE *fpOut;

	/* fopen files */
	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: HTS_Aln2Unique, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: HTS_Aln2Unique, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	/* create strand */
	vLastStrand = NULL;
	vLastStrand = (int *)calloc(nMax, sizeof(int));

	/* read file */
	strcpy(strLastChr, "");
	nLastPos = -1;
	nLastC = 0;
	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		SeqClust_ReadPos_From_Aln(strLine, strChr, &nPos, &nStrand);

		if( (strcmp(strChr, strLastChr) == 0) && (nPos == nLastPos) )
		{
			if(nLastC < nMax)
				vLastStrand[nLastC] = nStrand;
			nLastC++;
		}
		else
		{
			if(nLastC > 0)
			{
				if(nLastC > nMax)
					nLastC = nMax;

				for(ni=0; ni<nLastC; ni++)
				{
					if(vLastStrand[ni] == 0)
						fprintf(fpOut, "%s\t%d\t+\n", strLastChr, nLastPos);
					else
						fprintf(fpOut, "%s\t%d\t-\n", strLastChr, nLastPos);
				}
			}

			strcpy(strLastChr, strChr);
			nLastPos = nPos;
			vLastStrand[0] = nStrand;
			nLastC = 1;
		}
	}

	if(nLastC > 0)
	{
		if(nLastC > nMax)
			nLastC = nMax;
		
		for(ni=0; ni<nLastC; ni++)
		{
			if(vLastStrand[ni] == 0)
				fprintf(fpOut, "%s\t%d\t+\n", strLastChr, nLastPos);
			else
				fprintf(fpOut, "%s\t%d\t-\n", strLastChr, nLastPos);
		}
	}

	/* close files */
	fclose(fpIn);
	fclose(fpOut);
	free(vLastStrand);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Window_Main()                                                  */
/*  Convert high throughput sequencing alignment to windowed bar tiling    */
/*  data. Each sample will be scaled to have the same number of aligned    */
/*  reads if specified.                                                    */
/* ----------------------------------------------------------------------- */ 
int HTS_Aln2Window_Main(char strParamPath[])
{	
	/* Sample name */
	struct tagHTSAln2WinParam *pParam = NULL;
	struct INTMATRIX *pChrLen = NULL;
	struct tagBARData *pBARData = NULL;
	int ni;
	double dCount;
	char strFileName[LONG_LINE_LENGTH];
	struct INTMATRIX *pCol;
	int nx,ny;
	FILE *fpOut;

	/* Load Parameters */
	pParam = HTS_Aln2Window_LoadParamter(strParamPath);
	if(pParam == NULL)
	{
		printf("Error: HTS_Aln2Window_Main, cannot allocate memory for loading parameters!\n");
		exit(EXIT_FAILURE);
	}

	/* Load Chromosome */
	pChrLen = IMLOAD(pParam->strChrLen);
	if(pChrLen == NULL)
	{
		printf("Error: HTS_Aln2Window_Main, cannot load chromosome length!\n");
		exit(EXIT_FAILURE);
	}

	/* prepare BAR object */
	pBARData = HTS_Aln2Window_PrepareBARData(pParam, pChrLen);
	if(pBARData == NULL)
	{
		printf("Error: HTS_Aln2Window_Main, cannot create bar data object!\n");
		exit(EXIT_FAILURE);
	}

	/* process data one by one */
	pCol = NULL;
	pCol = CreateIntMatrix(1,pBARData->nColNum);
	if(pCol == NULL)
	{
		printf("Error: HTS_Aln2Window_Main, cannot create output column information!\n");
		exit(EXIT_FAILURE);
	}
	pCol->pMatElement[0] = 1;
	pCol->pMatElement[1] = 1;
	for(ni=0; ni<pParam->nSampleNum; ni++)
	{
		/* count */
		sprintf(strFileName, "%s%s", pParam->strDataFolder, pParam->vSampleFile[ni]->m_pString);
		dCount = 0;
		dCount = HTS_Aln2Window_CountHits(pBARData, strFileName, pParam->nW);

		/* scaling */
		if(dCount > 0.0)
		{
			if(pParam->dScaling > 0.0)
			{
				/* HTS_Aln2Window_Scaling(pBARData, pParam->dScaling/dCount, pParam->nW); */
				/* HTS_Aln2Window_Scaling(pBARData, pParam->nW); */
			}
		}

		/* output */
		sprintf(strFileName, "%s%s.bar", pParam->strExportFolder, pParam->vSampleFile[ni]->m_pString);
		Affy_SaveBAR_Columns_Fast(strFileName, pBARData, pCol);
		sprintf(strFileName, "%s%s.bar.txt", pParam->strExportFolder, pParam->vSampleFile[ni]->m_pString);
		fpOut = NULL;
		fpOut = fopen(strFileName, "w");
		for(nx=0; nx<3; nx++)
		{
			for(ny=0; ny<pBARData->vSeqData[nx]->nDataNum; ny++)
			{
				fprintf(fpOut, "%f\n", pBARData->vSeqData[nx]->vData[1]->pMatElement[ny]);
			}
		}
		fclose(fpOut);
	}
	DestroyIntMatrix(pCol);

	/* write cgw file */
	HTS_Aln2Window_WriteCGW(pParam);
	
	/* write cgb file */
	HTS_Aln2Window_WriteCGB(pParam);

	/* Destroy Paramter */
	HTSAln2WinParam_Destroy(&pParam);
	DestroyIntMatrix(pChrLen);
	pChrLen = NULL;
	Affy_BARData_Destroy(&pBARData);

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  HTSAln2WinParam_Create()                                               */
/*  Create hts_aln2window paramter.                                        */
/* ----------------------------------------------------------------------- */ 
struct tagHTSAln2WinParam *HTSAln2WinParam_Create()
{
	/* define */
	struct tagHTSAln2WinParam *pParam = NULL;

	/* create */
	pParam = (struct tagHTSAln2WinParam *)calloc(1, sizeof(struct tagHTSAln2WinParam));
	if(pParam == NULL)
	{
		printf("Error: HTSAln2WinParam_Create, cannot allocate memory for creating HTSAln2WinParam object!\n");
		exit(EXIT_FAILURE);
	}

	/* return */
	return pParam;
}

/* ----------------------------------------------------------------------- */ 
/*  HTSAln2WinParam_Destroy()                                              */
/*  Destroy hts_aln2window paramter.                                       */
/* ----------------------------------------------------------------------- */ 
void HTSAln2WinParam_Destroy(struct tagHTSAln2WinParam **pParam)
{
	int ni;

	if(pParam == NULL)
		return;

	if(*pParam == NULL)
		return;

	DestroyIntMatrix((*pParam)->vGroupLabel);
	(*pParam)->vGroupLabel = NULL;

	for(ni=0; ni<(*pParam)->nSampleNum; ni++)
	{
		DeleteString((*pParam)->vSampleAlias[ni]);
		(*pParam)->vSampleAlias[ni] = NULL;
		DeleteString((*pParam)->vSampleFile[ni]);
		(*pParam)->vSampleFile[ni] = NULL;
	}
	free((*pParam)->vSampleAlias);
	(*pParam)->vSampleAlias = NULL;
	free((*pParam)->vSampleFile);
	(*pParam)->vSampleFile = NULL;

	free(*pParam);
	*pParam = NULL;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Window_LoadParamter()                                          */
/*  Load hts_aln2window paramter.                                          */
/* ----------------------------------------------------------------------- */ 
struct tagHTSAln2WinParam *HTS_Aln2Window_LoadParamter(char strParamPath[])
{
	/* define */
	struct tagHTSAln2WinParam *pParam = NULL;
	FILE *fpIn = NULL;
	char strLine[LONG_LINE_LENGTH];
	char strKey[MED_LINE_LENGTH];
	char strValue[LONG_LINE_LENGTH];
	char *chp;
	int ni;

	/* load */
	pParam = HTSAln2WinParam_Create();
	if(pParam == NULL)
	{
		printf("Error: HTS_Aln2Window_LoadParamter, cannot allocate memory for loading parameter!\n");
		exit(EXIT_FAILURE);
	}

	fpIn = fopen(strParamPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: HTS_Aln2Window_LoadParamter, cannot open input parameter file!\n");
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

		chp = strchr(strLine, '=');
		if(chp == NULL)
		{
			printf("Error: HTS_Aln2Window_LoadParamter, wrong parameter file!\n");
			exit(EXIT_FAILURE);
		}

		*chp = '\0';
		strcpy(strKey, strLine);
		StrTrimRight(strKey);

		chp++;
		strcpy(strValue, chp);
		StrTrimLeft(strValue);

		if(strcmp(strKey, "project_name") == 0)
		{
			strcpy(pParam->strProjectName, strValue);
		}
		else if(strcmp(strKey, "chrlist") == 0)
		{
			strcpy(pParam->strChrList, strValue);
		}
		else if(strcmp(strKey, "chrlen") == 0)
		{
			strcpy(pParam->strChrLen, strValue);
		}
		else if(strcmp(strKey, "data_folder") == 0)
		{
			strcpy(pParam->strDataFolder, strValue);
			AdjustDirectoryPath(pParam->strDataFolder);
		}
		else if(strcmp(strKey, "sample_num") == 0)
		{
			pParam->nSampleNum = atoi(strValue);
			if( pParam->nSampleNum <= 0 )
			{
				printf("Error: HTS_Aln2Window_LoadParamter, no sample available!\n");
				exit(EXIT_FAILURE);
			}

			pParam->vSampleAlias = (struct tagString **)calloc(pParam->nSampleNum, sizeof(struct tagString *));
			pParam->vSampleFile = (struct tagString **)calloc(pParam->nSampleNum, sizeof(struct tagString *));
			pParam->vGroupLabel = CreateIntMatrix(1, pParam->nSampleNum);

			if( (pParam->vSampleAlias == NULL) || (pParam->vSampleFile == NULL) || (pParam->vGroupLabel == NULL))
			{
				printf("Error: HTS_Aln2Window_LoadParamter, cannot create memory for loading sample information!\n");
				exit(EXIT_FAILURE);
			}

		}
		else if(strstr(strKey, "sample_alias") == strKey)
		{
			chp = strrchr(strKey, '_');
			if(chp == NULL)
			{
				printf("Error: HTS_Aln2Window_LoadParamter, wrong parameter file!\n");
				exit(EXIT_FAILURE);
			}
			chp++;
			ni = atoi(chp)-1;
			StringAddTail(pParam->vSampleAlias+ni, strValue);
		}
		else if(strstr(strKey, "sample_group") == strKey)
		{
			chp = strrchr(strKey, '_');
			if(chp == NULL)
			{
				printf("Error: HTS_Aln2Window_LoadParamter, wrong parameter file!\n");
				exit(EXIT_FAILURE);
			}
			chp++;
			ni = atoi(chp)-1;
			pParam->vGroupLabel->pMatElement[ni] = atoi(strValue);
		}
		else if(strstr(strKey, "sample_file") == strKey)
		{
			
			chp = strrchr(strKey, '_');
			if(chp == NULL)
			{
				printf("Error: HTS_Aln2Window_LoadParamter, wrong parameter file!\n");
				exit(EXIT_FAILURE);
			}
			chp++;
			ni = atoi(chp)-1;
			StringAddTail(pParam->vSampleFile+ni, strValue);
		}
		else if(strcmp(strKey, "group_num") == 0)
		{
			pParam->nGroupNum = atoi(strValue);
		}
		else if(strcmp(strKey, "scaling") == 0)
		{
			pParam->dScaling = atof(strValue);
		}
		else if(strcmp(strKey, "window_size") == 0)
		{
			pParam->nW = atoi(strValue);
		}
		else if(strcmp(strKey, "step_size") == 0)
		{
			pParam->nS = atoi(strValue);
		}
		else if(strcmp(strKey, "export_folder") == 0)
		{
			strcpy(pParam->strExportFolder, strValue);
			AdjustDirectoryPath(pParam->strExportFolder);
		}
	}

	fclose(fpIn);

	/* return */
	return pParam;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Window_PrepareBARData()                                        */
/*  Prepare the bar data object.                                           */
/* ----------------------------------------------------------------------- */ 
struct tagBARData * HTS_Aln2Window_PrepareBARData(struct tagHTSAln2WinParam *pParam, 
				struct INTMATRIX *pChrLen)
{
	/* define */
	struct tagBARData *pBARData;
	int nSeqNum;
	int nColNum;
	int nDataNum;
	FILE *fpIn;
	char strLine[LONG_LINE_LENGTH];
	int ni,nj,nk;
	
	/* init */
	nSeqNum = pChrLen->nHeight;
	nColNum = 2;

	/* create */
	pBARData = Affy_BARData_Create();
	if(pBARData == NULL)
	{
		printf("Error: HTS_Aln2Window_PrepareBARData, BARData object was not created!\n");
		exit(EXIT_FAILURE);
	}
	strcpy(pBARData->strMagicnumber, "barr\r\n\032\n");
	pBARData->fVersionnumber = 2.0;
    pBARData->nSeqNum = nSeqNum;
	pBARData->nColNum = nColNum;
	pBARData->nParamNum = 0;
	pBARData->vParamName = NULL;
	pBARData->vParamValue = NULL;
	pBARData->pFieldType = CreateIntMatrix(1, pBARData->nColNum);
	if(pBARData->pFieldType == NULL)
	{
		printf("Error: HTS_Aln2Window_PrepareBARData, cannot allocate memory for field type.\n");
		exit(EXIT_FAILURE);
	}
	pBARData->pFieldType->pMatElement[0] = 2;
	for(ni=1; ni<pBARData->nColNum; ni++)
		pBARData->pFieldType->pMatElement[ni] = 1;

	pBARData->vSeqData = (struct tagBARSeq **)calloc(pBARData->nSeqNum, sizeof(struct tagBARSeq *));
	if(pBARData->vSeqData == NULL)
	{
		printf("Error: HTS_Aln2Window_PrepareBARData, cannot allocate memory for loading sequence data.\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<pBARData->nSeqNum; ni++)
	{
		/* create BARSeq object */
		pBARData->vSeqData[ni] = Affy_BARSeq_Create();
		if(pBARData->vSeqData[ni] == NULL)
		{
			printf("Error: HTS_Aln2Window_PrepareBARData, cannot create BARSeq object.\n");
			exit(EXIT_FAILURE);
		}

		pBARData->vSeqData[ni]->pSeqGroupName = NULL;
		pBARData->vSeqData[ni]->pSeqVersion = NULL;
		pBARData->vSeqData[ni]->nParamNum = 0;
		pBARData->vSeqData[ni]->vParamName = NULL;
		pBARData->vSeqData[ni]->vParamValue = NULL;

		pBARData->vSeqData[ni]->nColNum = pBARData->nColNum;

		nDataNum = pChrLen->pMatElement[ni]/pParam->nW;
		if(pChrLen->pMatElement[ni]%pParam->nW != 0)
			nDataNum++;

		pBARData->vSeqData[ni]->nDataNum = nDataNum;
		pBARData->vSeqData[ni]->vData = (struct DOUBLEMATRIX **)calloc(pBARData->vSeqData[ni]->nColNum, sizeof(struct DOUBLEMATRIX *));
		if(pBARData->vSeqData[ni]->vData == NULL)
		{
			printf("Error: HTS_Aln2Window_PrepareBARData, cannot allocate memory for loading bpmap coordinate data!\n");
			exit(EXIT_FAILURE);
		}
		if(pBARData->vSeqData[ni]->nDataNum > 0)
		{
			for(nj=0; nj<pBARData->vSeqData[ni]->nColNum; nj++)
			{
				pBARData->vSeqData[ni]->vData[nj] = CreateDoubleMatrix(1, pBARData->vSeqData[ni]->nDataNum );
				if(pBARData->vSeqData[ni]->vData[nj] == NULL)
				{
					printf("Error: HTS_Aln2Window_PrepareBARData, cannot allocate memory for loading bpmap coordinate data!\n");
					exit(EXIT_FAILURE);
				}
			}
		}

		nk = pParam->nW/2;
		for(nj=0; nj<nDataNum; nj++)
		{
			pBARData->vSeqData[ni]->vData[0]->pMatElement[nj] = nk;
			nk += pParam->nW;
		}
	}

	fpIn = NULL;
	fpIn = fopen(pParam->strChrList, "r");
	if(fpIn == NULL)
	{
		printf("Error: HTS_Aln2Window_PrepareBARData, cannot open chrlist file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		if(ni >= nSeqNum)
		{
			printf("Error: HTS_Aln2Window_PrepareBARData, chromosome number wrong!\n");
			exit(EXIT_FAILURE);
		}

		StringAddTail(&(pBARData->vSeqData[ni]->pSeqName), strLine);
		ni++;
	}
	
	fclose(fpIn);

	if(ni != nSeqNum)
	{
		printf("Error: HTS_Aln2Window_PrepareBARData, chromosome number wrong!\n");
		exit(EXIT_FAILURE);
	}

	/* return */
	return pBARData;
}


/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Window_CountHits()                                             */
/*  Count hits in each window.                                             */
/* ----------------------------------------------------------------------- */ 
double HTS_Aln2Window_CountHits(struct tagBARData *pBARData, char strInFile[], int nW)
{
	/* define */
	FILE *fpIn = NULL;
	char strLine[LONG_LINE_LENGTH];
	char strChr[MED_LINE_LENGTH];
	int nPos,nId;
	int ni;
	double dCount = 0.0;

	/* reset matrix */
	for(ni=0; ni<pBARData->nSeqNum; ni++)
	{
		DestroyDoubleMatrix(pBARData->vSeqData[ni]->vData[1]);
		pBARData->vSeqData[ni]->vData[1] = NULL;
		pBARData->vSeqData[ni]->vData[1] = CreateDoubleMatrix(1, pBARData->vSeqData[ni]->nDataNum);
		if(pBARData->vSeqData[ni]->vData[1] == NULL)
		{
			printf("Error: HTS_Aln2Window_PrepareBARData, cannot create memory for counting!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* open file */
	fpIn = fopen(strInFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: HTS_Aln2Window_PrepareBARData, cannot open input file!\n");
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

		sscanf(strLine, "%s %d", strChr, &nPos);
		for(ni=0; ni<pBARData->nSeqNum; ni++)
		{
			if(strcmp(strChr, pBARData->vSeqData[ni]->pSeqName->m_pString) == 0)
			{
				nId = nPos/nW;
				if(nId < pBARData->vSeqData[ni]->nDataNum)
				{
					pBARData->vSeqData[ni]->vData[1]->pMatElement[nId] += 1;
					dCount += 1.0;
				}
				break;
			}
		}
	}


	/* close file */
	fclose(fpIn);

	/* return */ 
	return dCount;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Window_Scaling()                                               */
/*  Scaling HTS data.                                                      */
/* ----------------------------------------------------------------------- */ 
int HTS_Aln2Window_Scaling(struct tagBARData *pBARData, double dScale)
{
	int ni,nj;

	/* 
	for(ni=0; ni<pBARData->nSeqNum; ni++)
	{
		for(nj=0; nj<pBARData->vSeqData[ni]->nDataNum; nj++)
		{
			pBARData->vSeqData[ni]->vData[1]->pMatElement[nj] *= dScale;
		}
	} */

	double dM = 0.0;
	double dV = 0.0;
	double dTotalC = 0.0;
	double dA,dB;
	double dSum = 0.0;

	for(ni=0; ni<pBARData->nSeqNum; ni++)
	{
		dTotalC += pBARData->vSeqData[ni]->nDataNum;
		for(nj=0; nj<pBARData->vSeqData[ni]->nDataNum; nj++)
		{
			dM += pBARData->vSeqData[ni]->vData[1]->pMatElement[nj];
		}
	}
	dM /= dTotalC;

	for(ni=0; ni<pBARData->nSeqNum; ni++)
	{
		for(nj=0; nj<pBARData->vSeqData[ni]->nDataNum; nj++)
		{
			dV += (pBARData->vSeqData[ni]->vData[1]->pMatElement[nj]-dM)*(pBARData->vSeqData[ni]->vData[1]->pMatElement[nj]-dM);
		}
	}
	dV /= (dTotalC-1.0);

	dB = dM/(dV-dM);
	dA = dB*dM;


	for(ni=0; ni<pBARData->nSeqNum; ni++)
	{
		for(nj=0; nj<pBARData->vSeqData[ni]->nDataNum; nj++)
		{
			pBARData->vSeqData[ni]->vData[1]->pMatElement[nj] = (pBARData->vSeqData[ni]->vData[1]->pMatElement[nj]+dA)/(1.0+dB);
			dSum += pBARData->vSeqData[ni]->vData[1]->pMatElement[nj];
		}
	}

	for(ni=0; ni<pBARData->nSeqNum; ni++)
	{
		for(nj=0; nj<pBARData->vSeqData[ni]->nDataNum; nj++)
		{
			pBARData->vSeqData[ni]->vData[1]->pMatElement[nj] /= dSum;
		}
	}

	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Window_WriteCGW()                                              */
/*  Write to CGW file.                                                     */
/* ----------------------------------------------------------------------- */ 
int HTS_Aln2Window_WriteCGW(struct tagHTSAln2WinParam *pParam)
{
	/* define */
	FILE *fpOut;
	char strFileName[LONG_LINE_LENGTH];
	int ni;

	/* write to CGW file */
	sprintf(strFileName, "%s%s.cgw", pParam->strExportFolder, pParam->strProjectName);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: HTS_Aln2Window_WriteCGW, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "[item1]\n");
	fprintf(fpOut, "type=bartilingexp\n");
	fprintf(fpOut, "item_name=%s\n", pParam->strProjectName);
	fprintf(fpOut, "bar_folder=%s\n", pParam->strExportFolder);
	fprintf(fpOut, "lib_num=1\n");
	fprintf(fpOut, "sample_num=%d\n", pParam->nSampleNum);
	fprintf(fpOut, "group_num=%d\n", pParam->nGroupNum);
	for(ni=0; ni<pParam->nSampleNum; ni++)
	{
		fprintf(fpOut, "sample_alias_%d=%s\n", ni+1, pParam->vSampleAlias[ni]->m_pString);
		fprintf(fpOut, "sample_group_%d=%d\n", ni+1, pParam->vGroupLabel->pMatElement[ni]);
		fprintf(fpOut, "bar_file_%d_1=%s.bar\n", ni+1, pParam->vSampleFile[ni]->m_pString);
	}
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Window_WriteCGB()                                              */
/*  Write to CGB file.                                                     */
/* ----------------------------------------------------------------------- */
int HTS_Aln2Window_WriteCGB(struct tagHTSAln2WinParam *pParam)
{
	/* define */
	FILE *fpOut;
	char strFileName[LONG_LINE_LENGTH];
	char strCommand[LONG_LINE_LENGTH];
	int ni,nj;
	int nTrackNum;
	char strInName[MED_LINE_LENGTH];
	char strOutName[MED_LINE_LENGTH];

	/* write to CGW file */
	nTrackNum = 2*pParam->nSampleNum;

	sprintf(strFileName, "%s%s.ini", pParam->strExportFolder, pParam->strProjectName);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: HTS_Aln2Window_WriteCGB, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "[session]\n");
	fprintf(fpOut, "type=genome\n");
	fprintf(fpOut, "[genome]\n");
	fprintf(fpOut, "num_tracks=%d\n", nTrackNum);
	fprintf(fpOut, "region=chr1:1-5000\n");

	nj = 1;
	for(ni=0; ni<pParam->nSampleNum; ni++)
	{
		sprintf(strInName, "%s%s", pParam->strDataFolder, pParam->vSampleFile[ni]->m_pString);
		sprintf(strOutName, "%s%s.sort", pParam->strExportFolder, pParam->vSampleFile[ni]->m_pString);
		sprintf(strCommand, "tablesorter %s %s", strInName, strOutName);
		system(strCommand);

		strcpy(strInName, strOutName);
		sprintf(strOutName, "%s.raw", pParam->vSampleFile[ni]->m_pString);
		TileMapv2_TXT2BAR(strInName, pParam->strExportFolder, strOutName);

		fprintf(fpOut, "[track%d]\n", nj);
		fprintf(fpOut, "title=%s_raw\n", pParam->vSampleAlias[ni]->m_pString);
		fprintf(fpOut, "type=signal\n");
		fprintf(fpOut, "pic_height=50\n");
		fprintf(fpOut, "src_filename=%s%s.raw.bar\n", pParam->strExportFolder, pParam->vSampleFile[ni]->m_pString);
		nj++;

		fprintf(fpOut, "[track%d]\n", nj);
		fprintf(fpOut, "title=%s_win\n", pParam->vSampleAlias[ni]->m_pString);
		fprintf(fpOut, "type=signal\n");
		fprintf(fpOut, "src_filename=%s%s.bar\n", pParam->strExportFolder, pParam->vSampleFile[ni]->m_pString);
		nj++;
	}
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Diff_Main()                                                    */
/*  Detect differentially aligned regions.                                 */
/* ----------------------------------------------------------------------- */ 
int HTS_Aln2Diff_Main(char strParamPath[])
{	
	/* Sample name */
	struct tagHTSAln2WinParam *pParam = NULL;
	struct INTMATRIX *pChrLen = NULL;
	struct tagBARData *pBARData = NULL;
	int ni;
	double dCount;
	char strFileName[LONG_LINE_LENGTH];
	char strInName[MED_LINE_LENGTH];
	char strOutName[MED_LINE_LENGTH];
	char strCommand[MED_LINE_LENGTH];
	struct INTMATRIX *pCol;
	double *vP;
	struct DOUBLEMATRIX *pCount = NULL;

	/* Load Parameters */
	pParam = HTS_Aln2Window_LoadParamter(strParamPath);
	if(pParam == NULL)
	{
		printf("Error: HTS_Aln2Diff_Main, cannot allocate memory for loading parameters!\n");
		exit(EXIT_FAILURE);
	}

	/* Load Chromosome */
	pChrLen = IMLOAD(pParam->strChrLen);
	if(pChrLen == NULL)
	{
		printf("Error: HTS_Aln2Diff_Main, cannot load chromosome length!\n");
		exit(EXIT_FAILURE);
	}

	/* prepare BAR object */
	pBARData = HTS_Aln2Diff_PrepareBARData(pParam, pChrLen);
	if(pBARData == NULL)
	{
		printf("Error: HTS_Aln2Diff_Main, cannot create bar data object!\n");
		exit(EXIT_FAILURE);
	}

	/* process data one by one */
	pCol = NULL;
	pCol = CreateIntMatrix(1,pBARData->nColNum);
	if(pCol == NULL)
	{
		printf("Error: HTS_Aln2Diff_Main, cannot create output column information!\n");
		exit(EXIT_FAILURE);
	}
	pCol->pMatElement[0] = 1;
	for(ni=0; ni<pParam->nSampleNum; ni++)
	{
		/* count */
		sprintf(strFileName, "%s%s", pParam->strDataFolder, pParam->vSampleFile[ni]->m_pString);
		dCount = 0;
		dCount = HTS_Aln2Diff_CountHits(pBARData, strFileName, pParam->nW, ni);

		/* output */
		pCol->pMatElement[1+ni] = 1;
		sprintf(strFileName, "%s%s.bar", pParam->strExportFolder, pParam->vSampleFile[ni]->m_pString);
		Affy_SaveBAR_Columns_Fast(strFileName, pBARData, pCol);
		pCol->pMatElement[1+ni] = 0;
	}
	DestroyIntMatrix(pCol);

	/* pCount = CreateDoubleMatrix(1001,1001);
	sprintf(strTemp, "%sdave_stat.txt", pParam->strExportFolder);
	for(ni=0; ni<pBARData->nSeqNum; ni++)
	{
		for(nj=0; nj<pBARData->vSeqData[ni]->nDataNum; nj++)
		{
			nx = (int)pBARData->vSeqData[ni]->vData[1]->pMatElement[nj];
			nk = nx+(int)pBARData->vSeqData[ni]->vData[2]->pMatElement[nj];

			if(nk>1000)
				nk = 1000;
			if(nx>1000)
				nx = 1000;
			nt = (int)DMGETAT(pCount, nk, nx)+1;
			DMSETAT(pCount, nk, nx, nt);
		}
	}
	DMSAVE(pCount, strTemp);
	DestroyDoubleMatrix(pCount);	
	*/

	/* estimate paramters */
	vP = NULL;
	/* vP = (double *)calloc(pParam->nSampleNum+2, sizeof(double)); */
	/* vP = (double *)calloc(2, sizeof(double)); */
	vP = (double *)calloc(9, sizeof(double)); 
	if(vP == NULL)
	{
		printf("Error: HTS_Aln2Diff_Main, cannot create vP!\n");
		exit(EXIT_FAILURE);
	}
	
	/* HTS_Aln2Diff_EstimateP(pBARData, vP, pParam->nSampleNum); */

	/* estimate prior probability */
	/* dPriorP = 0.01;
	HTS_Aln2Diff_EstimatePrior2Col(pBARData, pParam->nSampleNum, vP, &dPriorP); */
	/* HTS_Aln2Diff_EstimatePBinom(pBARData, pParam->nSampleNum, vP); */
	
	HTS_Aln2Diff_EstimateP3(pBARData, pParam->nSampleNum, vP);

	/* output */
	for(ni=0; ni<pParam->nSampleNum; ni++)
	{
		sprintf(strInName, "%s%s", pParam->strDataFolder, pParam->vSampleFile[ni]->m_pString);
		sprintf(strOutName, "%s%s.sort", pParam->strExportFolder, pParam->vSampleFile[ni]->m_pString);
		sprintf(strCommand, "tablesorter %s %s", strInName, strOutName);
		system(strCommand);
	}
	
	HTS_Aln2Diff_CallRegion(pParam, vP);
	
	/* Destroy Paramter */
	free(vP);
	HTSAln2WinParam_Destroy(&pParam);
	DestroyIntMatrix(pChrLen);
	pChrLen = NULL;
	Affy_BARData_Destroy(&pBARData);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Diff_PrepareBARData()                                          */
/*  Prepare the bar data object.                                           */
/* ----------------------------------------------------------------------- */ 
struct tagBARData * HTS_Aln2Diff_PrepareBARData(struct tagHTSAln2WinParam *pParam, 
				struct INTMATRIX *pChrLen)
{
	/* define */
	struct tagBARData *pBARData;
	int nSeqNum;
	int nColNum;
	int nDataNum;
	FILE *fpIn;
	char strLine[LONG_LINE_LENGTH];
	int ni,nj,nk;
	
	/* init */
	nSeqNum = pChrLen->nHeight;
	nColNum = pParam->nSampleNum+1;

	/* create */
	pBARData = Affy_BARData_Create();
	if(pBARData == NULL)
	{
		printf("Error: HTS_Aln2Diff_PrepareBARData, BARData object was not created!\n");
		exit(EXIT_FAILURE);
	}
	strcpy(pBARData->strMagicnumber, "barr\r\n\032\n");
	pBARData->fVersionnumber = 2.0;
    pBARData->nSeqNum = nSeqNum;
	pBARData->nColNum = nColNum;
	pBARData->nParamNum = 0;
	pBARData->vParamName = NULL;
	pBARData->vParamValue = NULL;
	pBARData->pFieldType = CreateIntMatrix(1, pBARData->nColNum);
	if(pBARData->pFieldType == NULL)
	{
		printf("Error: HTS_Aln2Diff_PrepareBARData, cannot allocate memory for field type.\n");
		exit(EXIT_FAILURE);
	}
	pBARData->pFieldType->pMatElement[0] = 2;
	for(ni=1; ni<pBARData->nColNum; ni++)
		pBARData->pFieldType->pMatElement[ni] = 1;

	pBARData->vSeqData = (struct tagBARSeq **)calloc(pBARData->nSeqNum, sizeof(struct tagBARSeq *));
	if(pBARData->vSeqData == NULL)
	{
		printf("Error: HTS_Aln2Diff_PrepareBARData, cannot allocate memory for loading sequence data.\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<pBARData->nSeqNum; ni++)
	{
		/* create BARSeq object */
		pBARData->vSeqData[ni] = Affy_BARSeq_Create();
		if(pBARData->vSeqData[ni] == NULL)
		{
			printf("Error: HTS_Aln2Diff_PrepareBARData, cannot create BARSeq object.\n");
			exit(EXIT_FAILURE);
		}

		pBARData->vSeqData[ni]->pSeqGroupName = NULL;
		pBARData->vSeqData[ni]->pSeqVersion = NULL;
		pBARData->vSeqData[ni]->nParamNum = 0;
		pBARData->vSeqData[ni]->vParamName = NULL;
		pBARData->vSeqData[ni]->vParamValue = NULL;

		pBARData->vSeqData[ni]->nColNum = pBARData->nColNum;

		nDataNum = pChrLen->pMatElement[ni]/pParam->nW;
		if(pChrLen->pMatElement[ni]%pParam->nW != 0)
			nDataNum++;

		pBARData->vSeqData[ni]->nDataNum = nDataNum;
		pBARData->vSeqData[ni]->vData = (struct DOUBLEMATRIX **)calloc(pBARData->vSeqData[ni]->nColNum, sizeof(struct DOUBLEMATRIX *));
		if(pBARData->vSeqData[ni]->vData == NULL)
		{
			printf("Error: HTS_Aln2Diff_PrepareBARData, cannot allocate memory for loading bpmap coordinate data!\n");
			exit(EXIT_FAILURE);
		}
		if(pBARData->vSeqData[ni]->nDataNum > 0)
		{
			for(nj=0; nj<pBARData->vSeqData[ni]->nColNum; nj++)
			{
				pBARData->vSeqData[ni]->vData[nj] = CreateDoubleMatrix(1, pBARData->vSeqData[ni]->nDataNum );
				if(pBARData->vSeqData[ni]->vData[nj] == NULL)
				{
					printf("Error: HTS_Aln2Diff_PrepareBARData, cannot allocate memory for loading bpmap coordinate data!\n");
					exit(EXIT_FAILURE);
				}
			}
		}

		nk = pParam->nW/2;
		for(nj=0; nj<nDataNum; nj++)
		{
			pBARData->vSeqData[ni]->vData[0]->pMatElement[nj] = nk;
			nk += pParam->nW;
		}
	}

	fpIn = NULL;
	fpIn = fopen(pParam->strChrList, "r");
	if(fpIn == NULL)
	{
		printf("Error: HTS_Aln2Diff_PrepareBARData, cannot open chrlist file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		if(ni >= nSeqNum)
		{
			printf("Error: HTS_Aln2Diff_PrepareBARData, chromosome number wrong!\n");
			exit(EXIT_FAILURE);
		}

		StringAddTail(&(pBARData->vSeqData[ni]->pSeqName), strLine);
		ni++;
	}
	
	fclose(fpIn);

	if(ni != nSeqNum)
	{
		printf("Error: HTS_Aln2Diff_PrepareBARData, chromosome number wrong!\n");
		exit(EXIT_FAILURE);
	}

	/* return */
	return pBARData;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Diff_PrepareRecountingData()                                   */
/*  Prepare the bar data object.                                           */
/* ----------------------------------------------------------------------- */ 
struct tagBARData * HTS_Aln2Diff_PrepareRecountingData(struct tagHTSAln2WinParam *pParam, 
				struct INTMATRIX *pChrLen)
{
	/* define */
	struct tagBARData *pBARData;
	int nSeqNum;
	int nColNum;
	int nDataNum;
	FILE *fpIn;
	char strLine[LONG_LINE_LENGTH];
	int ni,nj,nk;
	
	/* init */
	nSeqNum = pChrLen->nHeight;
	nColNum = pParam->nSampleNum+2;

	/* create */
	pBARData = Affy_BARData_Create();
	if(pBARData == NULL)
	{
		printf("Error: HTS_Aln2Diff_PrepareRecountingData, BARData object was not created!\n");
		exit(EXIT_FAILURE);
	}
	strcpy(pBARData->strMagicnumber, "barr\r\n\032\n");
	pBARData->fVersionnumber = 2.0;
    pBARData->nSeqNum = nSeqNum;
	pBARData->nColNum = nColNum;
	pBARData->nParamNum = 0;
	pBARData->vParamName = NULL;
	pBARData->vParamValue = NULL;
	pBARData->pFieldType = CreateIntMatrix(1, pBARData->nColNum);
	if(pBARData->pFieldType == NULL)
	{
		printf("Error: HTS_Aln2Diff_PrepareRecountingData, cannot allocate memory for field type.\n");
		exit(EXIT_FAILURE);
	}
	pBARData->pFieldType->pMatElement[0] = 2;
	for(ni=1; ni<pBARData->nColNum; ni++)
		pBARData->pFieldType->pMatElement[ni] = 1;

	pBARData->vSeqData = (struct tagBARSeq **)calloc(pBARData->nSeqNum, sizeof(struct tagBARSeq *));
	if(pBARData->vSeqData == NULL)
	{
		printf("Error: HTS_Aln2Diff_PrepareRecountingData, cannot allocate memory for loading sequence data.\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<pBARData->nSeqNum; ni++)
	{
		/* create BARSeq object */
		pBARData->vSeqData[ni] = Affy_BARSeq_Create();
		if(pBARData->vSeqData[ni] == NULL)
		{
			printf("Error: HTS_Aln2Diff_PrepareRecountingData, cannot create BARSeq object.\n");
			exit(EXIT_FAILURE);
		}

		pBARData->vSeqData[ni]->pSeqGroupName = NULL;
		pBARData->vSeqData[ni]->pSeqVersion = NULL;
		pBARData->vSeqData[ni]->nParamNum = 0;
		pBARData->vSeqData[ni]->vParamName = NULL;
		pBARData->vSeqData[ni]->vParamValue = NULL;

		pBARData->vSeqData[ni]->nColNum = pBARData->nColNum;

		nDataNum = pChrLen->pMatElement[ni]/pParam->nS;
		if(pChrLen->pMatElement[ni]%pParam->nS != 0)
			nDataNum++;

		pBARData->vSeqData[ni]->nDataNum = nDataNum;
		pBARData->vSeqData[ni]->vData = (struct DOUBLEMATRIX **)calloc(pBARData->vSeqData[ni]->nColNum, sizeof(struct DOUBLEMATRIX *));
		if(pBARData->vSeqData[ni]->vData == NULL)
		{
			printf("Error: HTS_Aln2Diff_PrepareRecountingData, cannot allocate memory for loading bpmap coordinate data!\n");
			exit(EXIT_FAILURE);
		}
		if(pBARData->vSeqData[ni]->nDataNum > 0)
		{
			for(nj=0; nj<pBARData->vSeqData[ni]->nColNum; nj++)
			{
				pBARData->vSeqData[ni]->vData[nj] = CreateDoubleMatrix(1, pBARData->vSeqData[ni]->nDataNum );
				if(pBARData->vSeqData[ni]->vData[nj] == NULL)
				{
					printf("Error: HTS_Aln2Diff_PrepareRecountingData, cannot allocate memory for loading bpmap coordinate data!\n");
					exit(EXIT_FAILURE);
				}
			}
		}

		nk = pParam->nS/2;
		for(nj=0; nj<nDataNum; nj++)
		{
			pBARData->vSeqData[ni]->vData[0]->pMatElement[nj] = nk;
			nk += pParam->nS;
		}
	}

	fpIn = NULL;
	fpIn = fopen(pParam->strChrList, "r");
	if(fpIn == NULL)
	{
		printf("Error: HTS_Aln2Diff_PrepareRecountingData, cannot open chrlist file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		if(ni >= nSeqNum)
		{
			printf("Error: HTS_Aln2Diff_PrepareRecountingData, chromosome number wrong!\n");
			exit(EXIT_FAILURE);
		}

		StringAddTail(&(pBARData->vSeqData[ni]->pSeqName), strLine);
		ni++;
	}
	
	fclose(fpIn);

	if(ni != nSeqNum)
	{
		printf("Error: HTS_Aln2Diff_PrepareRecountingData, chromosome number wrong!\n");
		exit(EXIT_FAILURE);
	}

	/* return */
	return pBARData;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Diff_CountHits()                                               */
/*  Count hits in each window.                                             */
/* ----------------------------------------------------------------------- */ 
double HTS_Aln2Diff_CountHits(struct tagBARData *pBARData, char strInFile[], int nW, int nSampleID)
{
	/* define */
	FILE *fpIn = NULL;
	char strLine[LONG_LINE_LENGTH];
	char strChr[MED_LINE_LENGTH];
	int nPos,nId;
	int ni;
	double dCount = 0.0;
	int nCol;

	/* open file */
	nCol = nSampleID+1;

	fpIn = fopen(strInFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: HTS_Aln2Diff_CountHits, cannot open input file!\n");
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

		sscanf(strLine, "%s %d", strChr, &nPos);
		for(ni=0; ni<pBARData->nSeqNum; ni++)
		{
			if(strcmp(strChr, pBARData->vSeqData[ni]->pSeqName->m_pString) == 0)
			{
				nId = nPos/nW;
				if(nId < pBARData->vSeqData[ni]->nDataNum)
				{
					pBARData->vSeqData[ni]->vData[nCol]->pMatElement[nId] += 1;
					dCount += 1.0;
				}
				break;
			}
		}
	}


	/* close file */
	fclose(fpIn);

	/* return */ 
	return dCount;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Diff_EstimateP()                                               */
/*  Estimate hyperparameter.                                               */
/* ----------------------------------------------------------------------- */ 
int HTS_Aln2Diff_EstimateP(struct tagBARData *pBARData, double *vP, int nSampleNum)
{
	int ni,nj,nk;
	double *vM = NULL;
	double *vV = NULL;

	double dM = 0.0;
	double dV = 0.0;
	double dTotalC = 0.0;
	double dA,dB;
	double dSum = 0.0;

	vM = (double *)calloc(nSampleNum, sizeof(double));
	vV = (double *)calloc(nSampleNum, sizeof(double));

	if(vM == NULL || vV == NULL)
	{
		printf("Error: HTS_Aln2Diff_EstimateP, cannot estimate parameter!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<pBARData->nSeqNum; ni++)
	{
		dTotalC += pBARData->vSeqData[ni]->nDataNum;
		for(nk=1; nk<=nSampleNum; nk++)
		{
			for(nj=0; nj<pBARData->vSeqData[ni]->nDataNum; nj++)
			{
				vM[nk-1] += pBARData->vSeqData[ni]->vData[nk]->pMatElement[nj];
			}
		}
	}
	for(nk=0; nk<nSampleNum; nk++)
	{
		vM[nk] /= dTotalC;
	}

	for(ni=0; ni<pBARData->nSeqNum; ni++)
	{
		for(nk=1; nk<=nSampleNum; nk++)
		{
			for(nj=0; nj<pBARData->vSeqData[ni]->nDataNum; nj++)
			{
				vV[nk-1] += (pBARData->vSeqData[ni]->vData[nk]->pMatElement[nj]-vM[nk-1])*(pBARData->vSeqData[ni]->vData[nk]->pMatElement[nj]-vM[nk-1]);
			}
		}
	}
	for(nk=0; nk<nSampleNum; nk++)
	{
		vV[nk] /= (dTotalC-1.0);
	}

	dA = 0.0;
	for(nk=0; nk<nSampleNum; nk++)
	{
		vP[nk] = vM[nk]*vM[nk]/(vV[nk]-vM[nk]);
		dA += vP[nk];
	}
	dA /= nSampleNum;

	dB = 0.0;
	for(nk=0; nk<nSampleNum; nk++)
	{
		dB += (vP[nk]-dA)*(vP[nk]-dA);
	}
	dB = sqrt(dB/(nSampleNum-1.0));

	printf("NegBin CV(alpha) = %f\n", dB/dA);
	if(dB/dA > 0.1)
	{
		printf("Warning: when NegBin CV(alpha) > 0.1, the model assumptions used here may not hold!\n");
	}

	vP[nSampleNum] = dA;
	for(nk=0; nk<nSampleNum; nk++)
	{
		vP[nk] = vP[nSampleNum]/vM[nk];
	}
	vP[nSampleNum+1] = vP[0];
	for(nk=0; nk<nSampleNum; nk++)
	{
		vP[nk] = vP[nSampleNum+1]/vP[nk];
	}

	free(vM);
	free(vV);

	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Diff_EstimatePrior2Col()                                       */
/*  Estimate prior probability of difference                               */
/* ----------------------------------------------------------------------- */ 
int HTS_Aln2Diff_EstimatePrior2Col(struct tagBARData *pBARData, int nSampleNum, 
								   double *vP, double *dPriorP)
{
	/* define */
	double dL0,dL1;
	int ni,nj;
	double vC[2];
	double dX,dY,dT,dT1,dT2;
	double dA,dB;
	double dTemp = *dPriorP;

	*dPriorP = dTemp/2;

	while( fabs(dTemp-(*dPriorP))/(*dPriorP) > 0.01 )
	{
		*dPriorP = dTemp;

		for(ni=0; ni<2; ni++)
			vC[ni] = 0.0;

		for(ni=0; ni<pBARData->nSeqNum; ni++)
		{
			for(nj=0; nj<pBARData->vSeqData[ni]->nDataNum; nj++)
			{
				dA = vP[nSampleNum];
				dB = vP[nSampleNum+1];

				dX = pBARData->vSeqData[ni]->vData[1]->pMatElement[nj];
				dY = pBARData->vSeqData[ni]->vData[2]->pMatElement[nj];

				dT = vP[0]+vP[1]+dB;

				dL0 = gammaln(dX+dY+dA)-gammaln(dA)-gammaln(dX+1.0)-gammaln(dY+1.0)
					+dA*log(dB/dT)+dX*log(vP[0]/dT)+dY*log(vP[1]/dT);

				dT1 = vP[0]+dB;
				dT2 = vP[1]+dB;
				dL1 = gammaln(dX+dA)-gammaln(dA)-gammaln(dX+1.0)+dA*log(dB/dT1)+dX*log(vP[0]/dT1)
					+gammaln(dY+dA)-gammaln(dA)-gammaln(dY+1.0)+dA*log(dB/dT2)+dY*log(vP[1]/dT2);

				dTemp = (1-(*dPriorP))*exp(dL0-dL1)/(*dPriorP);
				dTemp = 1.0/(1.0+dTemp);
				vC[1] += dTemp;
				vC[0] += 1.0-dTemp;
			}
		}

		dTemp = vC[1]/(vC[1]+vC[0]);
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Diff_EstimatePBinom()                                          */
/*  Estimate hyperparamter for mixture binominal                           */
/* ----------------------------------------------------------------------- */ 
int HTS_Aln2Diff_EstimatePBinom(struct tagBARData *pBARData, int nSampleNum, 
								   double *vP)
{
	/* define */
	double dL0,dL1;
	int ni,nj;
	double vC[2];
	double dX,dY;
	double dA,dB;
	double dTemp;
	double dTotal;

	vC[0] = 0.01;
	vC[1] = 0.0;
	
	dA = 0.0;
	dB = 0.0;
	dTotal = 0.0;
	for(ni=0; ni<pBARData->nSeqNum; ni++)
	{
		dTotal += pBARData->vSeqData[ni]->nDataNum;
		for(nj=0; nj<pBARData->vSeqData[ni]->nDataNum; nj++)
		{
			dA += pBARData->vSeqData[ni]->vData[1]->pMatElement[nj];
			dB += pBARData->vSeqData[ni]->vData[2]->pMatElement[nj];		
		}
	}

	vC[1] = dA/(dA+dB);
	vP[0] = vC[0]/2.0;
	vP[1] = vC[1]/2.0;

	while( (fabs(vC[0]-vP[0])/vP[0]>0.01) || (fabs(vC[1]-vP[1])/vP[1]>0.01) )
	{
		for(ni=0; ni<2; ni++)
		{
			vP[ni] = vC[ni];
			vC[ni] = 0.0;
		}

			
		dA = 0.0;
		dB = 0.0;
		for(ni=0; ni<pBARData->nSeqNum; ni++)
		{
			for(nj=0; nj<pBARData->vSeqData[ni]->nDataNum; nj++)
			{
				dX = pBARData->vSeqData[ni]->vData[1]->pMatElement[nj];
				dY = pBARData->vSeqData[ni]->vData[2]->pMatElement[nj];

				dL0 = log(1.0-vP[0])+gammaln(dX+dY+1.0)-gammaln(dX+1.0)-gammaln(dY+1.0)+dX*log(vP[1])+dY*log(1.0-vP[1]);

				dL1 = log(vP[0])-log(dX+dY+1.0);

				dTemp = exp(dL0-dL1);
				dTemp = 1.0/(1.0+dTemp);
				vC[0] += dTemp;
				dA += (1.0-dTemp)*dX;
				dB += (1.0-dTemp)*dY;
			}
		}

		vC[0] /= dTotal;
		vC[1] = dA/(dA+dB);
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Diff_EstimateP3()                                              */
/*  Estimate hyperparamter for mixture binominal                           */
/* ----------------------------------------------------------------------- */ 
int HTS_Aln2Diff_EstimateP3(struct tagBARData *pBARData, int nSampleNum, 
								   double *vP)
{
	/* define */
	int nMinC = 2;
	int nParamNum = 9;
	double vC[9];
	double vBackUp[9];
	double dQ0,dQ1,dQ2;
	double dLam0,dLam1,dPro0,dPro1;
	double dBetaP;
	double dL0,dL1,dL2;
	int ni,nj,nk;
	double dX,dY,dN;
	double dA,dB,dC,dD;
	double dTemp,dTemp2;
	double dTotal;
	double dError;
	double logBetaP, logBetaPC;
	double logPi0, logPi1, logPi2, logP0, logP0C, logLambda, logGab, logGn;
	double logGam;
	double dLMax;
	int nMaxIter = 20;
	int nIter;
	double dLike0,dLike1;
	double logGalpha0,logGalpha;
	double dM1Lam,dM2Lam,dM1Beta,dM2Beta,dMTBeta;
	struct DOUBLEMATRIX *pDQ0,*pDQ1,*pDQ2,*pDLike;
	int nDQSIZE;
	double vSumDist[256];
	int nTestParameter = 0;
	int nSuppressDisplay = 0;
	char strFileName[255];
	FILE *fpOut;


	vC[0] = 0.98; /* pi_0 */
	vC[1] = 0.01; /* pi_1 */
	vC[2] = 0.01; /* pi_2 */
	vC[3] = 0.05; /* lambda_0 */
	vC[4] = 0.5; /* p_0 */
	vC[5] = 0.05; /* beta */
	vC[6] = 1.0; /* alpha */
	vC[7] = 1.0; /* a */
	vC[8] = 1.0; /* b */

	nDQSIZE = 21;
	pDQ0 = NULL;
	pDQ0 = CreateDoubleMatrix(nDQSIZE,nDQSIZE);
	pDQ1 = NULL;
	pDQ1 = CreateDoubleMatrix(nDQSIZE,nDQSIZE);
	pDQ2 = NULL;
	pDQ2 = CreateDoubleMatrix(nDQSIZE,nDQSIZE);
	pDLike = NULL;
	pDLike = CreateDoubleMatrix(nDQSIZE,nDQSIZE);
	if( (pDQ0 == NULL) || (pDQ1 == NULL) || (pDQ2 == NULL) || (pDLike == NULL) )
	{
		printf("Error: cannot create Q database!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<256; ni++)
	{
		vSumDist[ni] = 0.0;
	}

	dA = 0.0;
	dB = 0.0;
	dC = 0.0;
	dD = 0.0;
	dTotal = 0.0;
	for(ni=0; ni<pBARData->nSeqNum; ni++)
	{
		dTotal += pBARData->vSeqData[ni]->nDataNum;
		for(nj=0; nj<pBARData->vSeqData[ni]->nDataNum; nj++)
		{
			dA += pBARData->vSeqData[ni]->vData[1]->pMatElement[nj];
			dB += pBARData->vSeqData[ni]->vData[2]->pMatElement[nj];
			nk = (int)(pBARData->vSeqData[ni]->vData[1]->pMatElement[nj]+pBARData->vSeqData[ni]->vData[2]->pMatElement[nj]);
			if(nk > 255)
				nk = 255;
			vSumDist[nk] += 1;
			if(nk >= nMinC)
			{
				dC += nk;
				dD += 1.0;
			}
		}
	}

	vC[3] = (dA+dB)/dTotal;
	if((double)nMinC < vC[3])
	{
		nMinC = (int)(vC[3])+2;
	}
	vC[4] = dA/(dA+dB);
	dC /= dD;
	if(dC < nMinC)
		dC = nMinC+1;
	vC[5] = 1/(2.0*dC);
	dD = vC[3];

	dTemp = 0.0;
	for(nk=0; nk<256; nk++)
	{
		dTemp += vSumDist[nk];
	}
	for(nk=0; nk<256; nk++)
	{
		vSumDist[nk] /= dTemp;
	}
	strcpy(strFileName, "wincount_summary.txt");
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	for(nk=0; nk<256; nk++)
	{
		fprintf(fpOut, "%f\n", vSumDist[nk]);
	}
	fclose(fpOut);

	for(ni=0; ni<nParamNum; ni++)
	{
		vP[ni] = vC[ni]/2.0;
	}
	dError = 1.0;

	dLike0 = -1e20;
	dLike1 = -1e20;
	nIter = 0;
	nTestParameter = 0;
	while( dError > 0.01 )
	{
		if(nIter%10 == 0)
		{
			if(nSuppressDisplay == 1)
			{
				nSuppressDisplay = 0;
			}
			else
			{
				printf("EM iteration %d ...\n", nIter);
			}
		}

		if(nTestParameter == 1)
		{
			for(ni=0; ni<5; ni++)
			{
				vBackUp[ni] = vC[ni];
				vP[ni] = vC[ni];
				vC[ni] = 0.0;
			}
			for(; ni<nParamNum; ni++)
			{
				vBackUp[ni] = vP[ni];
				vP[ni] = vC[ni];
				vC[ni] = 0.0;
			}
		}
		else
		{
			for(ni=0; ni<nParamNum; ni++)
			{
				vP[ni] = vC[ni];
				vC[ni] = 0.0;
			}
		}

		dA = 0.0;
		dB = 0.0;
		dLam0 = 0.0;
		dLam1 = 0.0;
		dPro0 = 0.0;
		dPro1 = 0.0;
		dLike0 = dLike1;
		dLike1 = 0.0;
		dM1Lam = 0.0;
		dM2Lam = 0.0;
		dM1Beta = 0.0;
		dM2Beta = 0.0;
		dMTBeta = 0.0;

		logPi0 = log(vP[0]);
		logPi1 = log(vP[1]);
		logPi2 = log(vP[2]);

		logLambda = log(vP[3]);

		logP0 = log(vP[4]);
		logP0C = log(1.0-vP[4]);
		
		dBetaP = 1.0/(1.0+vP[5]);
		logBetaP = log(dBetaP);
		logBetaPC = log(1.0-dBetaP);

		logGab = gammaln(vP[7]+vP[8])-gammaln(vP[7])-gammaln(vP[8]);
		logGalpha0 = gammaln(vP[6]);
			

		for(ni=0; ni<nDQSIZE; ni++)
		{
			for(nj=0; nj<nDQSIZE; nj++)
			{
				dX = ni;
				dY = nj;
				dN = dX+dY;
			
				if((int)dN < nMinC)
				{
					logGam = gammaln(dX+1.0)+gammaln(dY+1.0);
					dL0 = logPi0+dN*logLambda-vP[3]+dX*logP0+dY*logP0C-logGam;
					dQ0 = 1.0;
					dQ1 = 0.0;
					dQ2 = 0.0;
					DMSETAT(pDQ0, ni, nj, dQ0);
					DMSETAT(pDQ1, ni, nj, dQ1);
					DMSETAT(pDQ2, ni, nj, dQ2);
					DMSETAT(pDLike, ni, nj, dL0);
				}
				else
				{
					logGam = gammaln(dX+1.0)+gammaln(dY+1.0);
					logGalpha = gammaln(dN+vP[6]-nMinC);
					logGn = 0.0;
					for(nk=0; nk<nMinC; nk++)
						logGn += log(dN-nk);

					dL0 = logPi0+dN*logLambda-vP[3]+dX*logP0+dY*logP0C-logGam;
					dLMax = dL0;

					dL1 = logPi1+logGalpha-logGalpha0+vP[6]*logBetaPC+(dN-nMinC)*logBetaP+dX*logP0+dY*logP0C-logGam+logGn;
					if(dL1 > dLMax)
						dLMax = dL1;

					dL2 = logPi2+logGalpha-logGalpha0+vP[6]*logBetaPC+(dN-nMinC)*logBetaP+logGab-logGam+logGn+gammaln(dX+vP[7])+gammaln(dY+vP[8])-gammaln(dN+vP[7]+vP[8]);
					if(dL2 > dLMax)
						dLMax = dL2;

					dL0 -= dLMax;
					dL1 -= dLMax;
					dL2 -= dLMax;

					dL0 = exp(dL0);
					dL1 = exp(dL1);
					dL2 = exp(dL2);
					dTemp = dL0+dL1+dL2;
					dQ0 = dL0/dTemp;
					dQ1 = dL1/dTemp;
					dQ2 = dL2/dTemp;

					DMSETAT(pDQ0, ni, nj, dQ0);
					DMSETAT(pDQ1, ni, nj, dQ1);
					DMSETAT(pDQ2, ni, nj, dQ2);
					dTemp = dLMax+log(dTemp);
					DMSETAT(pDLike, ni, nj, dTemp);
				}
			}
		}

		for(ni=0; ni<pBARData->nSeqNum; ni++)
		{
			for(nj=0; nj<pBARData->vSeqData[ni]->nDataNum; nj++)
			{
				dX = pBARData->vSeqData[ni]->vData[1]->pMatElement[nj];
				dY = pBARData->vSeqData[ni]->vData[2]->pMatElement[nj];
				dN = dX+dY;

				if( ((int)dX<nDQSIZE) && ((int)dY<nDQSIZE) )
				{
					dQ0 = DMGETAT(pDQ0, (int)dX, (int)dY);
					dQ1 = DMGETAT(pDQ1, (int)dX, (int)dY);
					dQ2 = DMGETAT(pDQ2, (int)dX, (int)dY);
					dLike1 += DMGETAT(pDLike, (int)dX, (int)dY);
				}
				else
				{
					logGam = gammaln(dX+1.0)+gammaln(dY+1.0);
					logGalpha = gammaln(dN+vP[6]-nMinC);
					logGn = 0.0;
					for(nk=0; nk<nMinC; nk++)
						logGn += log(dN-nk);

					dL0 = logPi0+dN*logLambda-vP[3]+dX*logP0+dY*logP0C-logGam;
					dLMax = dL0;

					dL1 = logPi1+logGalpha-logGalpha0+vP[6]*logBetaPC+(dN-nMinC)*logBetaP+dX*logP0+dY*logP0C-logGam+logGn;
					if(dL1 > dLMax)
						dLMax = dL1;

					dL2 = logPi2+logGalpha-logGalpha0+vP[6]*logBetaPC+(dN-nMinC)*logBetaP+logGab-logGam+logGn+gammaln(dX+vP[7])+gammaln(dY+vP[8])-gammaln(dN+vP[7]+vP[8]);
					if(dL2 > dLMax)
						dLMax = dL2;

					dL0 -= dLMax;
					dL1 -= dLMax;
					dL2 -= dLMax;

					dL0 = exp(dL0);
					dL1 = exp(dL1);
					dL2 = exp(dL2);
					dTemp = dL0+dL1+dL2;
					dQ0 = dL0/dTemp;
					dQ1 = dL1/dTemp;
					dQ2 = dL2/dTemp;

					dLike1 += dLMax+log(dTemp);
				}

				vC[0] += dQ0;
				vC[1] += dQ1;
				vC[2] += dQ2;

				dLam0 += dQ0;
				dLam1 += dQ0*dN;

				dPro0 += (dQ0+dQ1)*dN;
				dPro1 += (dQ0+dQ1)*dX;

				dA += (dQ1+dQ2);
				dB += (dQ1+dQ2)*dN;

				dM1Lam += (dQ1+dQ2)*dN;
				dM2Lam += (dQ1+dQ2)*dN*dN;
				/* if(dN > 0.5)
				{ */
					dTemp = (dX+vP[7])/(dN+vP[7]+vP[8]);
					dM1Beta += dQ2*dTemp;
					dM2Beta += dQ2*dTemp*dTemp;
					dMTBeta += dQ2;
				/* } */
			}
		}

		dM1Lam /= dA;
		dM2Lam /= dA;
		dM1Beta /= dMTBeta;
		dM2Beta /= dMTBeta;

		if(nTestParameter == 1)
		{
			nTestParameter = 0;
			if(dLike1 < dLike0)
			{
				dLike1 = dLike0;
				for(nk=0; nk<nParamNum; nk++)
				{
					vC[nk] = vBackUp[nk];
				}
				dError = 1.0;
				nSuppressDisplay = 1;
				continue;
			}
		}

		dTotal = vC[0]+vC[1]+vC[2];
		vC[0] /= dTotal;
		vC[1] /= dTotal;
		vC[2] /= dTotal;
		vC[3] = dLam1/dLam0;
		if(vC[3] > dD)
			vC[3] = dD;

		/* vC[4] = vP[4]; */
		vC[4] = dPro1/dPro0;
		/* vC[5] = dA/dB;
		if(vC[5] > 1.0)
			vC[5] = 1.0; */
		vC[5] = vP[5];
		vC[6] = vP[6];
		vC[7] = vP[7];
		vC[8] = vP[8];

		if(dLike1 < dLike0)
		{
			printf("Warning: likelihood decreasing!\n");
			/* break; */
		}

		dError = 0.0;
		for(ni=0; ni<nParamNum; ni++)
		{
			dTemp = fabs(vC[ni]-vP[ni])/vP[ni];
			if(dTemp > dError)
				dError = dTemp;
		}

		
		nIter++;

		if(nIter%5 == 0)
		{
			if(nIter%10 == 0)
			{
				dTemp2 = dM1Beta*dM1Beta*(1.0-dM1Beta)/(dM2Beta-dM1Beta*dM1Beta)-dM1Beta;
				if( dTemp2 > 0.0 )
				{
					vC[7] = dTemp2;
					vC[8] = (1.0-dM1Beta)*vC[7]/dM1Beta;
				}
			}
			else
			{
				dTemp = dM1Lam/(dM2Lam-dM1Lam*dM1Lam-dM1Lam);
				
				if(dM1Lam >= dC)
				{
					if( dTemp > 0.0 )
					{
						vC[5] = dTemp;
						vC[6] = dM1Lam*vC[5];
					}
				}
			}

			nTestParameter = 1;
		}

		if(nIter >= nMaxIter)
			break;
	}


	DestroyDoubleMatrix(pDQ0);
	DestroyDoubleMatrix(pDQ1);
	DestroyDoubleMatrix(pDQ2);
	DestroyDoubleMatrix(pDLike);

	printf("Estimated Parameters [Relative Error]:\n");
	printf("pi_0 = %f [%f%%]\n", vP[0], 100.0*fabs(vC[0]-vP[0])/vP[0]);
	printf("pi_1 = %f [%f%%]\n", vP[1], 100.0*fabs(vC[1]-vP[1])/vP[1]);
	printf("pi_2 = %f [%f%%]\n", vP[2], 100.0*fabs(vC[2]-vP[2])/vP[2]);
	printf("lambda_0 = %f [%f%%]\n", vP[3], 100.0*fabs(vC[3]-vP[3])/vP[3]);
	printf("p_0 = %f [%f%%]\n", vP[4], 100.0*fabs(vC[4]-vP[4])/vP[4]);
	printf("alpha = %f [%f%%]\n", vP[6], 100.0*fabs(vC[6]-vP[6])/vP[6]);
	printf("beta = %f [%f%%]\n", vP[5], 100.0*fabs(vC[5]-vP[5])/vP[5]);
	printf("a = %f [%f%%]\n", vP[7], 100.0*fabs(vC[7]-vP[7])/vP[7]);
	printf("b = %f [%f%%]\n", vP[8], 100.0*fabs(vC[8]-vP[8])/vP[8]);

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Diff_CallRegion()                                              */
/*  Call differentially expressed regions.                                 */
/* ----------------------------------------------------------------------- */ 
int HTS_Aln2Diff_CallRegion(struct tagHTSAln2WinParam *pParam, double *vP)
{
	/* define */
	FILE *vfpIn[2];
	FILE *fpOut;
	int ni,nj,nk;
	char strFileName[MED_LINE_LENGTH];
	char strOutFileName[MED_LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	char strChr1[LINE_LENGTH];
	char strChr2[LINE_LENGTH];
	char strChr[LINE_LENGTH];
	int nW1,nW2;
	int nR1,nR2;
	int nPos1,nPos2;
	int nX,nY;
	int nl,nu;
	int nStrCmpR;
	int nFile1End,nFile2End;
	struct tagBARData *pBARData;
	int nWrite1;
	int nStart,nEnd;

	int nMinC = 2;
	int nParamNum = 9;
	double dQ0,dQ1,dQ2;
	double dLam0,dLam1,dPro0,dPro1;
	double dBetaP;
	double dL0,dL1,dL2;
	double dX,dY,dN;
	double dA,dB,dC,dD;
	double dTemp,dTemp2;
	double dTotal;
	double logBetaP, logBetaPC;
	double logPi0, logPi1, logPi2, logP0, logP0C, logLambda, logGab, logGn;
	double logGam;
	double dLMax;
	double dLike0,dLike1;
	double logGalpha0,logGalpha;
	struct DOUBLEMATRIX *pDQ0,*pDQ1,*pDQ2;
	int nDQSIZE = 21;
	double dMaxQ = 0.0;
	double dMaxR = 0.0;
	int nMaxX = 0;
	int nMaxY = 0;
	int nRID = 0;

	/* prepare display file */
	for(ni=0; ni<2; ni++)
	{
		vfpIn[ni] = NULL;
		sprintf(strFileName, "%s%s.sort", pParam->strExportFolder, pParam->vSampleFile[ni]->m_pString);
		vfpIn[ni] = fopen(strFileName, "r");
		if(vfpIn[ni] == NULL)
		{
			printf("Error: HTS_Aln2Diff_CallRegion, cannot open input file!\n");
			exit(EXIT_FAILURE);
		}

		fpOut = NULL;
		sprintf(strFileName, "%s%s.sortm", pParam->strExportFolder, pParam->vSampleFile[ni]->m_pString);
		fpOut = fopen(strFileName, "w");
		if(fpOut == NULL)
		{
			printf("Error: HTS_Aln2Diff_CallRegion, cannot open output file!\n");
			exit(EXIT_FAILURE);
		}

		fprintf(fpOut, "#chr\tpos\t%s.raw\n", pParam->vSampleFile[ni]->m_pString);
		fprintf(fpOut, "#chr\tpos\t1\n");

		strcpy(strChr2, "");
		nPos2 = -1;
		dX = 0.0;
		while(fgets(strLine, LONG_LINE_LENGTH, vfpIn[ni]) != NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;
			if(strLine[0] == '#')
				continue;

			sscanf(strLine, "%s %d", strChr1, &nPos1);
			if( (strcmp(strChr1, strChr2) == 0) && (nPos1 == nPos2) )
			{
				dX += 1.0;
			}
			else
			{
				if(strcmp(strChr2, "") != 0)
				{
					fprintf(fpOut, "%s\t%d\t%f\n", strChr2, nPos2, dX);
				}
				strcpy(strChr2, strChr1);
				nPos2 = nPos1;
				dX = 1.0;
			}
		}

		if(strcmp(strChr2, "") != 0)
		{
			fprintf(fpOut, "%s\t%d\t%f\n", strChr2, nPos2, dX);
		}

		fclose(vfpIn[ni]);
		fclose(fpOut);

		sprintf(strOutFileName, "%s_raw.cgw", pParam->vSampleFile[ni]->m_pString);
		TileMapv2_TXT2BAR(strFileName, pParam->strExportFolder, strOutFileName);
	}

	/* open files */
	for(ni=0; ni<2; ni++)
	{
		vfpIn[ni] = NULL;
		sprintf(strFileName, "%s%s.sort", pParam->strExportFolder, pParam->vSampleFile[ni]->m_pString);
		vfpIn[ni] = fopen(strFileName, "r");
		if(vfpIn[ni] == NULL)
		{
			printf("Error: HTS_Aln2Diff_CallRegion, cannot open input file!\n");
			exit(EXIT_FAILURE);
		}
	}

	fpOut = NULL;
	sprintf(strFileName, "%s%s.cod", pParam->strExportFolder, pParam->strProjectName);
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: HTS_Aln2Diff_CallRegion, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "#chr\tpos\t%s_raw\n", pParam->strProjectName);
	fprintf(fpOut, "#chr\tpos\t1\n");

	/* load file */
	nFile1End = 1;
	while(fgets(strLine, LONG_LINE_LENGTH, vfpIn[0]) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		sscanf(strLine, "%s %d", strChr1, &nPos1);
		nFile1End = 0;
		break;
	}

	nFile2End = 1;
	while(fgets(strLine, LONG_LINE_LENGTH, vfpIn[1]) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		sscanf(strLine, "%s %d", strChr2, &nPos2);
		nFile2End = 0;
		break;
	}

	while( (nFile1End == 0) || (nFile2End == 0) )
	{
		if(nFile1End == 1)
		{
			fprintf(fpOut, "%s\t%d\t-1\n", strChr2, nPos2);
			while(fgets(strLine, LONG_LINE_LENGTH, vfpIn[1]) != NULL)
			{
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] == '\0')
					continue;
				if(strLine[0] == '#')
					continue;

				sscanf(strLine, "%s %d", strChr2, &nPos2);
				fprintf(fpOut, "%s\t%d\t-1\n", strChr2, nPos2);
			}
			nFile2End = 1;
		}
		else if(nFile2End == 1)
		{
			fprintf(fpOut, "%s\t%d\t1\n", strChr1, nPos1);
			while(fgets(strLine, LONG_LINE_LENGTH, vfpIn[0]) != NULL)
			{
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] == '\0')
					continue;
				if(strLine[0] == '#')
					continue;

				sscanf(strLine, "%s %d", strChr1, &nPos1);
				fprintf(fpOut, "%s\t%d\t1\n", strChr1, nPos1);
			}
			nFile1End = 1;
		}
		else
		{
			nStrCmpR = strcmp(strChr1, strChr2);
			if(nStrCmpR < 0)
			{
				nWrite1 = 1;
			}
			else if(nStrCmpR == 0 && nPos1 <= nPos2)
			{
				nWrite1 = 1;
			}
			else
			{
				nWrite1 = 0;
			}

			if(nWrite1 == 1)
			{
				fprintf(fpOut, "%s\t%d\t1\n", strChr1, nPos1);
				nFile1End = 1;
				while(fgets(strLine, LONG_LINE_LENGTH, vfpIn[0]) != NULL)
				{
					StrTrimLeft(strLine);
					StrTrimRight(strLine);
					if(strLine[0] == '\0')
						continue;
					if(strLine[0] == '#')
						continue;

					sscanf(strLine, "%s %d", strChr1, &nPos1);
					nFile1End = 0;
					break;
				}
			}
			else
			{
				fprintf(fpOut, "%s\t%d\t-1\n", strChr2, nPos2);
				nFile2End = 1;
				while(fgets(strLine, LONG_LINE_LENGTH, vfpIn[1]) != NULL)
				{
					StrTrimLeft(strLine);
					StrTrimRight(strLine);
					if(strLine[0] == '\0')
						continue;
					if(strLine[0] == '#')
						continue;

					sscanf(strLine, "%s %d", strChr2, &nPos2);
					nFile2End = 0;
					break;
				}
			}

		}
	}

	/* close files */
	for(ni=0; ni<2; ni++)
		fclose(vfpIn[ni]);
	fclose(fpOut);


	/* convert txt to bar */
	sprintf(strOutFileName, "%s_raw.cgw", pParam->strProjectName);
	TileMapv2_TXT2BAR(strFileName, pParam->strExportFolder, strOutFileName);
	sprintf(strOutFileName, "%s%s_raw.bar", pParam->strExportFolder, pParam->strProjectName);
	pBARData = NULL;
	pBARData = Affy_LoadBAR_Fast(strOutFileName);
	if(pBARData == NULL)
	{
		printf("Error: cannot load bar data!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	sprintf(strFileName, "%s%s.cod", pParam->strExportFolder, pParam->strProjectName);
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: HTS_Aln2Diff_CallRegion, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	pDQ0 = NULL;
	pDQ0 = CreateDoubleMatrix(nDQSIZE,nDQSIZE);
	pDQ1 = NULL;
	pDQ1 = CreateDoubleMatrix(nDQSIZE,nDQSIZE);
	pDQ2 = NULL;
	pDQ2 = CreateDoubleMatrix(nDQSIZE,nDQSIZE);
	if( (pDQ0 == NULL) || (pDQ1 == NULL) || (pDQ2 == NULL) )
	{
		printf("Error: cannot create Q database!\n");
		exit(EXIT_FAILURE);
	}

	logPi0 = log(vP[0]);
	logPi1 = log(vP[1]);
	logPi2 = log(vP[2]);

	logLambda = log(vP[3]);

	logP0 = log(vP[4]);
	logP0C = log(1.0-vP[4]);
	
	dBetaP = 1.0/(1.0+vP[5]);
	logBetaP = log(dBetaP);
	logBetaPC = log(1.0-dBetaP);

	logGab = gammaln(vP[7]+vP[8])-gammaln(vP[7])-gammaln(vP[8]);
	logGalpha0 = gammaln(vP[6]);
		

	for(ni=0; ni<nDQSIZE; ni++)
	{
		for(nj=0; nj<nDQSIZE; nj++)
		{
			dX = ni;
			dY = nj;
			dN = dX+dY;
		
			if((int)dN < nMinC)
			{
				logGam = gammaln(dX+1.0)+gammaln(dY+1.0);
				dL0 = logPi0+dN*logLambda-vP[3]+dX*logP0+dY*logP0C-logGam;
				dQ0 = 1.0;
				dQ1 = 0.0;
				dQ2 = 0.0;
				DMSETAT(pDQ0, ni, nj, dQ0);
				DMSETAT(pDQ1, ni, nj, dQ1);
				DMSETAT(pDQ2, ni, nj, dQ2);
			}
			else
			{
				logGam = gammaln(dX+1.0)+gammaln(dY+1.0);
				logGalpha = gammaln(dN+vP[6]-nMinC);
				logGn = 0.0;
				for(nk=0; nk<nMinC; nk++)
					logGn += log(dN-nk);

				dL0 = logPi0+dN*logLambda-vP[3]+dX*logP0+dY*logP0C-logGam;
				dLMax = dL0;

				dL1 = logPi1+logGalpha-logGalpha0+vP[6]*logBetaPC+(dN-nMinC)*logBetaP+dX*logP0+dY*logP0C-logGam+logGn;
				if(dL1 > dLMax)
					dLMax = dL1;

				dL2 = logPi2+logGalpha-logGalpha0+vP[6]*logBetaPC+(dN-nMinC)*logBetaP+logGab-logGam+logGn+gammaln(dX+vP[7])+gammaln(dY+vP[8])-gammaln(dN+vP[7]+vP[8]);
				if(dL2 > dLMax)
					dLMax = dL2;

				dL0 -= dLMax;
				dL1 -= dLMax;
				dL2 -= dLMax;

				dL0 = exp(dL0);
				dL1 = exp(dL1);
				dL2 = exp(dL2);
				dTemp = dL0+dL1+dL2;
				dQ0 = dL0/dTemp;
				dQ1 = dL1/dTemp;
				dQ2 = dL2/dTemp;

				DMSETAT(pDQ0, ni, nj, dQ0);
				DMSETAT(pDQ1, ni, nj, dQ1);
				DMSETAT(pDQ2, ni, nj, dQ2);
			}
		}
	}

	/* check */
	nRID = 0;
	for(ni=0; ni<pBARData->nSeqNum; ni++)
	{
		if(pBARData->vSeqData[ni]->nDataNum <= 0)
			continue;

		nR1 = -1;
		nR2 = -1;

		nj = pBARData->vSeqData[ni]->nDataNum-1;
		nStart = (int)(pBARData->vSeqData[ni]->vData[0]->pMatElement[0]-pParam->nW+1);
		nEnd = (int)(pBARData->vSeqData[ni]->vData[0]->pMatElement[nj]+1);
		nX = 0;
		nY = 0;
		nl = 0;
		nu = 0;
		for(nW1=nStart; nW1<nEnd; nW1++)
		{
			if(nl < pBARData->vSeqData[ni]->nDataNum)
			{
				while( pBARData->vSeqData[ni]->vData[0]->pMatElement[nl] < nW1 )
				{
					if(pBARData->vSeqData[ni]->vData[1]->pMatElement[nl] > 0)
						nX--;
					else
						nY--;
					nl++;
					if(nl >= pBARData->vSeqData[ni]->nDataNum)
						break;
				}
			}
			
			if(nu < pBARData->vSeqData[ni]->nDataNum)
			{

				while( pBARData->vSeqData[ni]->vData[0]->pMatElement[nu] < (nW1+pParam->nW) )
				{
					if(pBARData->vSeqData[ni]->vData[1]->pMatElement[nu] > 0)
						nX++;
					else
						nY++;
					nu++;

					if(nu >= pBARData->vSeqData[ni]->nDataNum)
						break;

				}
			}

			/* compute */
			if( (nX < nDQSIZE) && (nY < nDQSIZE) )
			{
				dQ1 = DMGETAT(pDQ1, nX, nY);
				dQ2 = DMGETAT(pDQ2, nX, nY);
			}
			else
			{
				dX = nX;
				dY = nY;
				dN = dX+dY;

				logGam = gammaln(dX+1.0)+gammaln(dY+1.0);
				logGalpha = gammaln(dN+vP[6]-nMinC);
				logGn = 0.0;
				for(nk=0; nk<nMinC; nk++)
					logGn += log(dN-nk);

				dL0 = logPi0+dN*logLambda-vP[3]+dX*logP0+dY*logP0C-logGam;
				dLMax = dL0;

				dL1 = logPi1+logGalpha-logGalpha0+vP[6]*logBetaPC+(dN-nMinC)*logBetaP+dX*logP0+dY*logP0C-logGam+logGn;
				if(dL1 > dLMax)
					dLMax = dL1;

				dL2 = logPi2+logGalpha-logGalpha0+vP[6]*logBetaPC+(dN-nMinC)*logBetaP+logGab-logGam+logGn+gammaln(dX+vP[7])+gammaln(dY+vP[8])-gammaln(dN+vP[7]+vP[8]);
				if(dL2 > dLMax)
					dLMax = dL2;

				dL0 -= dLMax;
				dL1 -= dLMax;
				dL2 -= dLMax;

				dL0 = exp(dL0);
				dL1 = exp(dL1);
				dL2 = exp(dL2);
				dTemp = dL0+dL1+dL2;
				dQ0 = dL0/dTemp;
				dQ1 = dL1/dTemp;
				dQ2 = dL2/dTemp;
			}
			
			if(dQ2 > 0.5)
			{
				if(nR1 < 0)
				{
					nR1 = nW1;
					nR2 = nW1;
				}
				else
				{
					nR2 = nW1;
				}

				if(dQ2 > dMaxQ)
				{
					dMaxQ = dQ2;
					dMaxR = ((double)nX/(double)(nX+nY) - vP[4])/vP[4];
					nMaxX = nX;
					nMaxY = nY;
				}
			}
			else
			{
				if(nR1 >= 0)
				{
					/* save to file */
					nR2 += pParam->nW-1;
					if( (nMaxX+nMaxY) > 10)
					{
						fprintf(fpOut, "%d\t%s\t%d\t%d\t+\t%d\t%d\t%f\t%f\n", nRID, pBARData->vSeqData[ni]->pSeqName->m_pString,
							nR1, nR2, nMaxX, nMaxY, dMaxQ, dMaxR);
						nRID++;
					}
				}

				nR1 = -1;
				nR2 = -1;
				dMaxQ = 0.0;
				dMaxR = 0.0;
				nMaxX = 0;
				nMaxY = 0;
			}
		}

		if(nR1 >= 0)
		{
			/* save to file */
			nR2 += pParam->nW-1;
			if( (nMaxX+nMaxY) > 10)
			{
				fprintf(fpOut, "%d\t%s\t%d\t%d\t+\t%d\t%d\t%f\t%f\n", nRID, pBARData->vSeqData[ni]->pSeqName->m_pString,
					nR1, nR2, nMaxX, nMaxY, dMaxQ, dMaxR);
				nRID++;
			}
		}
	}

	fclose(fpOut);

	Affy_BARData_Destroy(&pBARData);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Enrich_Main()                                                  */
/*  Detect enriched alignment regions.                                     */
/* ----------------------------------------------------------------------- */ 
int HTS_Aln2Enrich_Main(char strParamPath[])
{	
	/* Sample name */
	struct tagHTSAln2WinParam *pParam = NULL;
	struct INTMATRIX *pChrLen = NULL;
	struct tagBARData *pBARData = NULL;
	int ni;
	double dCount;
	char strFileName[LONG_LINE_LENGTH];
	struct INTMATRIX *pCol;
	double *vP;

	/* Load Parameters */
	pParam = HTS_Aln2Window_LoadParamter(strParamPath);
	if(pParam == NULL)
	{
		printf("Error: HTS_Aln2Enrich_Main, cannot allocate memory for loading parameters!\n");
		exit(EXIT_FAILURE);
	}

	/* Load Chromosome */
	pChrLen = IMLOAD(pParam->strChrLen);
	if(pChrLen == NULL)
	{
		printf("Error: HTS_Aln2Enrich_Main, cannot load chromosome length!\n");
		exit(EXIT_FAILURE);
	}

	/* prepare BAR object */
	pBARData = HTS_Aln2Diff_PrepareBARData(pParam, pChrLen);
	if(pBARData == NULL)
	{
		printf("Error: HTS_Aln2Enrich_Main, cannot create bar data object!\n");
		exit(EXIT_FAILURE);
	}

	/* process data one by one */
	pCol = NULL;
	pCol = CreateIntMatrix(1,pBARData->nColNum);
	if(pCol == NULL)
	{
		printf("Error: HTS_Aln2Enrich_Main, cannot create output column information!\n");
		exit(EXIT_FAILURE);
	}
	pCol->pMatElement[0] = 2;
	for(ni=0; ni<pParam->nSampleNum; ni++)
	{
		/* count */
		sprintf(strFileName, "%s%s", pParam->strDataFolder, pParam->vSampleFile[ni]->m_pString);
		dCount = 0;
		dCount = HTS_Aln2Diff_CountHits(pBARData, strFileName, pParam->nW, ni);

		/* output */
		pCol->pMatElement[1+ni] = 1;
		sprintf(strFileName, "%s%s.bar", pParam->strExportFolder, pParam->vSampleFile[ni]->m_pString);
		Affy_SaveBAR_Columns_Fast(strFileName, pBARData, pCol);
		pCol->pMatElement[1+ni] = 0;
	}
	DestroyIntMatrix(pCol);
	
	/* estimate paramters */
	vP = NULL;
	vP = (double *)calloc(6, sizeof(double)); 
	if(vP == NULL)
	{
		printf("Error: HTS_Aln2Enrich_Main, cannot create vP!\n");
		exit(EXIT_FAILURE);
	}
	
	/* estimate prior probability */
	HTS_Aln2Enrich_EstimateP(pBARData, pParam->nSampleNum, vP);

	/* Destroy Paramter */
	free(vP);
	HTSAln2WinParam_Destroy(&pParam);
	DestroyIntMatrix(pChrLen);
	pChrLen = NULL;
	Affy_BARData_Destroy(&pBARData);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2Enrich_EstimateP()                                             */
/*  Estimate hyperparamter for mixture binominal                           */
/* ----------------------------------------------------------------------- */ 
int HTS_Aln2Enrich_EstimateP(struct tagBARData *pBARData, int nSampleNum, 
								   double *vP)
{
	/* define */
	int nMinC = 2;
	int nParamNum = 6;
	double vC[6];
	double vBackUp[6];
	double dQ0,dQ1;
	double dLam0,dLam1;
	double dBetaP;
	double dL0,dL1;
	int ni,nj,nk;
	double dN;
	double dA,dB,dC,dD;
	double dTemp;
	double dTotal;
	double dError;
	double logBetaP, logBetaPC;
	double logPi0, logPi1, logLambda, logGn;
	double logGam;
	double dLMax;
	int nMaxIter = 50;
	int nIter;
	double dLike0,dLike1;
	double logGalpha0,logGalpha;
	double dM1Lam,dM2Lam;
	struct DOUBLEMATRIX *pDQ0,*pDQ1,*pDLike;
	int nDQSIZE;
	double vSumDist[256];
	int nTestParameter = 0;
	int nSuppressDisplay = 0;
	char strFileName[255];
	FILE *fpOut;


	vC[0] = 0.99; /* pi_0 */
	vC[1] = 0.01; /* pi_1 */
	vC[2] = 0.05; /* lambda_0 */
	vC[3] = 0.5; /* p_0 */
	vC[4] = 0.05; /* beta */
	vC[5] = 1.0; /* alpha */
	
	nDQSIZE = 21;
	pDQ0 = NULL;
	pDQ0 = CreateDoubleMatrix(1,nDQSIZE);
	pDQ1 = NULL;
	pDQ1 = CreateDoubleMatrix(1,nDQSIZE);
	pDLike = NULL;
	pDLike = CreateDoubleMatrix(1,nDQSIZE);
	if( (pDQ0 == NULL) || (pDQ1 == NULL) || (pDLike == NULL) )
	{
		printf("Error: cannot create Q database!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<256; ni++)
	{
		vSumDist[ni] = 0.0;
	}

	dA = 0.0;
	dB = 0.0;
	dC = 0.0;
	dD = 0.0;
	dTotal = 0.0;
	for(ni=0; ni<pBARData->nSeqNum; ni++)
	{
		dTotal += pBARData->vSeqData[ni]->nDataNum;
		for(nj=0; nj<pBARData->vSeqData[ni]->nDataNum; nj++)
		{
			dA += pBARData->vSeqData[ni]->vData[1]->pMatElement[nj];
			nk = (int)(pBARData->vSeqData[ni]->vData[1]->pMatElement[nj]);
			if(nk > 255)
				nk = 255;
			vSumDist[nk] += 1;
			if(nk >= nMinC)
			{
				dC += nk;
				dD += 1.0;
			}
		}
	}

	vC[2] = dA/dTotal;
	if((double)nMinC < vC[2])
	{
		nMinC = (int)(vC[2])+2;
	}
	
	dC /= dD;
	if(dC < nMinC)
		dC = nMinC+1;
	vC[4] = 1/(2.0*dC);
	dD = vC[2];

	dTemp = 0.0;
	for(nk=0; nk<256; nk++)
	{
		dTemp += vSumDist[nk];
	}
	for(nk=0; nk<256; nk++)
	{
		vSumDist[nk] /= dTemp;
	}
	strcpy(strFileName, "wincount_summary.txt");
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	for(nk=0; nk<256; nk++)
	{
		fprintf(fpOut, "%e\n", vSumDist[nk]);
	}
	fclose(fpOut);

	for(ni=0; ni<nParamNum; ni++)
	{
		vP[ni] = vC[ni]/2.0;
	}
	dError = 1.0;

	dLike0 = -1e20;
	dLike1 = -1e20;
	nIter = 0;
	nTestParameter = 0;
	while( dError > 0.01 )
	{
		if(nIter%10 == 0)
		{
			if(nSuppressDisplay == 1)
			{
				nSuppressDisplay = 0;
			}
			else
			{
				printf("EM iteration %d ...\n", nIter);
			}
		}

		if(nTestParameter == 1)
		{
			for(ni=0; ni<4; ni++)
			{
				vBackUp[ni] = vC[ni];
				vP[ni] = vC[ni];
				vC[ni] = 0.0;
			}
			for(; ni<nParamNum; ni++)
			{
				vBackUp[ni] = vP[ni];
				vP[ni] = vC[ni];
				vC[ni] = 0.0;
			}
		}
		else
		{
			for(ni=0; ni<nParamNum; ni++)
			{
				vP[ni] = vC[ni];
				vC[ni] = 0.0;
			}
		}

		dA = 0.0;
		dB = 0.0;
		dLam0 = 0.0;
		dLam1 = 0.0;
		dLike0 = dLike1;
		dLike1 = 0.0;
		dM1Lam = 0.0;
		dM2Lam = 0.0;

		logPi0 = log(vP[0]);
		logPi1 = log(vP[1]);
		logLambda = log(vP[2]);		
		dBetaP = 1.0/(1.0+vP[4]);
		logBetaP = log(dBetaP);
		logBetaPC = log(1.0-dBetaP);
		logGalpha0 = gammaln(vP[5]);
	
		for(ni=0; ni<nDQSIZE; ni++)
		{
			dN = ni;
			if((int)dN < nMinC)
			{
				logGn = gammaln(dN+1.0);
				dL0 = logPi0+dN*logLambda-vP[2]-logGn;
				dQ0 = 1.0;
				dQ1 = 0.0;
				DMSETAT(pDQ0, 0, nj, dQ0);
				DMSETAT(pDQ1, 0, nj, dQ1);
				DMSETAT(pDLike, 0, nj, dL0);
			}
			else
			{
				logGn = gammaln(dN+1.0);
				logGam = gammaln(dN-nMinC);
				logGalpha = gammaln(dN+vP[5]-nMinC);
				
				dL0 = logPi0+dN*logLambda-vP[2]-logGn;
				dLMax = dL0;

				dL1 = logPi1+logGalpha-logGalpha0+vP[5]*logBetaPC+(dN-nMinC)*logBetaP-logGam;
				if(dL1 > dLMax)
					dLMax = dL1;

				dL0 -= dLMax;
				dL1 -= dLMax;

				dL0 = exp(dL0);
				dL1 = exp(dL1);
				dTemp = dL0+dL1;
				dQ0 = dL0/dTemp;
				dQ1 = dL1/dTemp;

				DMSETAT(pDQ0, 0, nj, dQ0);
				DMSETAT(pDQ1, 0, nj, dQ1);
				dTemp = dLMax+log(dTemp);
				DMSETAT(pDLike, 0, nj, dTemp);
			}
		}

		for(ni=0; ni<pBARData->nSeqNum; ni++)
		{
			for(nj=0; nj<pBARData->vSeqData[ni]->nDataNum; nj++)
			{
				dN = pBARData->vSeqData[ni]->vData[1]->pMatElement[nj];
				
				if((int)dN<nDQSIZE)
				{
					dQ0 = DMGETAT(pDQ0, 0, (int)dN);
					dQ1 = DMGETAT(pDQ1, 0, (int)dN);
					dLike1 += DMGETAT(pDLike, 0, (int)dN);
				}
				else
				{
					logGn = gammaln(dN+1.0);
					logGam = gammaln(dN-nMinC);
					logGalpha = gammaln(dN+vP[5]-nMinC);
					
					dL0 = logPi0+dN*logLambda-vP[2]-logGn;
					dLMax = dL0;

					dL1 = logPi1+logGalpha-logGalpha0+vP[5]*logBetaPC+(dN-nMinC)*logBetaP-logGam;
					if(dL1 > dLMax)
						dLMax = dL1;

					dL0 -= dLMax;
					dL1 -= dLMax;

					dL0 = exp(dL0);
					dL1 = exp(dL1);
					dTemp = dL0+dL1;
					dQ0 = dL0/dTemp;
					dQ1 = dL1/dTemp;
					dLike1 += dLMax+log(dTemp);
				}

				vC[0] += dQ0;
				vC[1] += dQ1;

				dLam0 += dQ0;
				dLam1 += dQ0*dN;

				dA += dQ1;
				dB += dQ1*dN;

				dM1Lam += dQ1*dN;
				dM2Lam += dQ1*dN*dN;
			}
		}

		dM1Lam /= dA;
		dM2Lam /= dA;

		if(nTestParameter == 1)
		{
			nTestParameter = 0;
			if(dLike1 < dLike0)
			{
				dLike1 = dLike0;
				for(nk=0; nk<nParamNum; nk++)
				{
					vC[nk] = vBackUp[nk];
				}
				dError = 1.0;
				nSuppressDisplay = 1;
				continue;
			}
		}

		dTotal = vC[0]+vC[1];
		vC[0] /= dTotal;
		vC[1] /= dTotal;
		vC[2] = dLam1/dLam0;
		if(vC[2] > dD)
			vC[2] = dD;

		vC[4] = vP[4];
		vC[5] = vP[5];
		
		if(dLike1 < dLike0)
		{
			printf("Warning: likelihood decreasing!\n");
			/* break; */
		}

		dError = 0.0;
		for(ni=0; ni<nParamNum; ni++)
		{
			dTemp = fabs(vC[ni]-vP[ni])/vP[ni];
			if(dTemp > dError)
				dError = dTemp;
		}

		
		nIter++;

		if(nIter%5 == 0)
		{
			dTemp = dM1Lam/(dM2Lam-dM1Lam*dM1Lam-dM1Lam);
				
			if(dM1Lam >= dC)
			{
				if( dTemp > 0.0 )
				{
					vC[4] = dTemp;
					vC[5] = dM1Lam*vC[4];
				}
			}

			nTestParameter = 1;
		}

		if(nIter >= nMaxIter)
			break;
	}


	DestroyDoubleMatrix(pDQ0);
	DestroyDoubleMatrix(pDQ1);
	DestroyDoubleMatrix(pDLike);

	printf("Estimated Parameters [Relative Error]:\n");
	printf("pi_0 = %f [%f%%]\n", vP[0], 100.0*fabs(vC[0]-vP[0])/vP[0]);
	printf("pi_1 = %f [%f%%]\n", vP[1], 100.0*fabs(vC[1]-vP[1])/vP[1]);
	printf("lambda_0 = %f [%f%%]\n", vP[2], 100.0*fabs(vC[2]-vP[2])/vP[2]);
	printf("alpha = %f [%f%%]\n", vP[5], 100.0*fabs(vC[5]-vP[5])/vP[5]);
	printf("beta = %f [%f%%]\n", vP[4], 100.0*fabs(vC[4]-vP[4])/vP[4]);
	
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2BAR()                                                          */
/*  Convert high throughput sequencing alignment to bar file.              */
/* ----------------------------------------------------------------------- */ 
int HTS_Aln2BAR(char strTXTFile[], char strBARFile[])
{
	/* define */
	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));

	struct tagBARData *pBARData = NULL;
	FILE *fpIn;
	FILE *fpOut;

	/* variables */
	int nFieldNum = 0;
	int nSeqNum = 0;
	int nColNum = 0;
	struct tagString *pSeqInfo = NULL;
	int nProbeNum = 0;
	int nTotalProbeNum = 0;
	char strLine[LONG_LINE_LENGTH];
	char strTemp[LONG_LINE_LENGTH];
	char *chp1,*chp2;
	char strChr[MED_LINE_LENGTH];
	char strLastChr[MED_LINE_LENGTH];
	int nLastPos;
	int nPos,nAlnCount;
	char strBARTmpFileName[MED_LINE_LENGTH];
	struct INTMATRIX *pCol = NULL;
	
	/* count */
	int ni,nj,nk;

	/* load */
	fpIn = NULL;
	fpIn = fopen(strTXTFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: HTS_Aln2BAR, cannot open source file!\n");
		exit(EXIT_FAILURE);
	}

	sprintf(strBARTmpFileName, "%s.tmp", strBARFile);
	fpOut = NULL;
	fpOut = fopen(strBARTmpFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: HTS_Aln2BAR, cannot open temporary output file!\n");
		exit(EXIT_FAILURE);
	}

	/* load input */
	pSeqInfo = NULL;
	strcpy(strLastChr, "");
	nLastPos = -1;
	nProbeNum = 0;
	nColNum = 1;
	nSeqNum = 0;
	nFieldNum = 3;

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '0')
			continue;
		if(strLine[0] == '#')
		{
			/* fprintf(fpOut, "%s\n", strLine); */
			continue;
		}

		sscanf(strLine, "%s %d", strChr, &nPos);
		
		if(strcmp(strChr, strLastChr) != 0)
		{
			if(strcmp(strLastChr, "") != 0)
			{
				sprintf(strTemp, "%s\t%d\t", strLastChr, nProbeNum);
				StringAddTail(&pSeqInfo, strTemp);
				fprintf(fpOut, "%s\t%d\t%d\n", strLastChr, nLastPos, nAlnCount);
			}

			nSeqNum++;
			strcpy(strLastChr, strChr);
			nLastPos = nPos;
			nProbeNum = 1;
			nAlnCount = 1;
		}
		else
		{
			if(nPos != nLastPos)
			{
				fprintf(fpOut, "%s\t%d\t%d\n", strLastChr, nLastPos, nAlnCount);
				nLastPos = nPos;
				nAlnCount = 1;
				nProbeNum++;
			}
			else
			{
				nAlnCount += 1;
			}
		}
	}

	if(strcmp(strLastChr, "") != 0)
	{
		sprintf(strTemp, "%s\t%d\t", strLastChr, nProbeNum);
		StringAddTail(&pSeqInfo, strTemp);
		fprintf(fpOut, "%s\t%d\t%d\n", strLastChr, nLastPos, nAlnCount);
	}

	fclose(fpIn);
	fclose(fpOut);


	/* load head */
	if(nSeqNum <= 0)
	{
		printf("Warning: empty sequences!\n");
		return PROC_FAILURE;
	}

	/* create BAR object */
	pBARData = Affy_BARData_Create();
	if(pBARData == NULL)
	{
		printf("Error: HTS_Aln2BAR, BARData object was not created!\n");
		exit(EXIT_FAILURE);
	}
	strcpy(pBARData->strMagicnumber, "barr\r\n\032\n");
	pBARData->fVersionnumber = 2.0;
    pBARData->nSeqNum = nSeqNum;
	pBARData->nColNum = nColNum+1;
	pBARData->nParamNum = 0;
	pBARData->vParamName = NULL;
	pBARData->vParamValue = NULL;
	pBARData->pFieldType = CreateIntMatrix(1, pBARData->nColNum);
	if(pBARData->pFieldType == NULL)
	{
		printf("Error: HTS_Aln2BAR, cannot allocate memory for field type.\n");
		exit(EXIT_FAILURE);
	}
	pBARData->pFieldType->pMatElement[0] = 2;
	for(ni=1; ni<pBARData->nColNum; ni++)
		pBARData->pFieldType->pMatElement[ni] = 1;

	pBARData->vSeqData = (struct tagBARSeq **)calloc(pBARData->nSeqNum, sizeof(struct tagBARSeq *));
	if(pBARData->vSeqData == NULL)
	{
		printf("Error: TileMapv2_TXT2BAR, cannot allocate memory for loading sequence data.\n");
		exit(EXIT_FAILURE);
	}

	chp1 = pSeqInfo->m_pString;
	for(ni=0; ni<pBARData->nSeqNum; ni++)
	{
		/* create BARSeq object */
		pBARData->vSeqData[ni] = Affy_BARSeq_Create();
		if(pBARData->vSeqData[ni] == NULL)
		{
			printf("Error: TileMapv2_TXT2BAR, cannot create BARSeq object.\n");
			exit(EXIT_FAILURE);
		}

		pBARData->vSeqData[ni]->pSeqGroupName = NULL;
		pBARData->vSeqData[ni]->pSeqVersion = NULL;
		pBARData->vSeqData[ni]->nParamNum = 0;
		pBARData->vSeqData[ni]->vParamName = NULL;
		pBARData->vSeqData[ni]->vParamValue = NULL;

		pBARData->vSeqData[ni]->nColNum = pBARData->nColNum;

		chp2 = strchr(chp1, '\t');
		*chp2 = '\0';
		strcpy(strTemp, chp1);
		StringAddTail(&(pBARData->vSeqData[ni]->pSeqName), strTemp);
		chp1 = chp2+1;

		chp2 = strchr(chp1, '\t');
		*chp2 = '\0';
		strcpy(strTemp, chp1);
		pBARData->vSeqData[ni]->nDataNum = atoi(strTemp);
		chp1 = chp2+1;


		pBARData->vSeqData[ni]->vData = (struct DOUBLEMATRIX **)calloc(pBARData->vSeqData[ni]->nColNum, sizeof(struct DOUBLEMATRIX *));
		if(pBARData->vSeqData[ni]->vData == NULL)
		{
			printf("Error: HTS_Aln2BAR, cannot allocate memory for loading bpmap coordinate data!\n");
			exit(EXIT_FAILURE);
		}
		if(pBARData->vSeqData[ni]->nDataNum > 0)
		{
			for(nj=0; nj<pBARData->vSeqData[ni]->nColNum; nj++)
			{
				pBARData->vSeqData[ni]->vData[nj] = CreateDoubleMatrix(1, pBARData->vSeqData[ni]->nDataNum );
				if(pBARData->vSeqData[ni]->vData[nj] == NULL)
				{
					printf("Error: HTS_Aln2BAR, cannot allocate memory for loading bpmap coordinate data!\n");
					exit(EXIT_FAILURE);
				}
			}
		}
	}
	
	/* load data */
	fpIn = NULL;
	fpIn = fopen(strBARTmpFileName, "r");
	if(fpIn == NULL)
	{
		printf("Error: HTS_Aln2BAR, cannot open source file!\n");
		exit(EXIT_FAILURE);
	}

	strcpy(strLastChr, "");
	nSeqNum = -1;
	nProbeNum = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '0')
			continue;
		if(strLine[0] == '#')
			continue;

		chp1 = strLine;
		chp2 = strchr(chp1, '\t');
		if(chp2 == NULL)
		{
			printf("Error: HTS_Aln2BAR, wrong input file format!\n");
			exit(EXIT_FAILURE);
		}
		*chp2 = '\0';
		StrTrimRight(chp1);
		strcpy(strChr, chp1);

		if(strcmp(strChr, strLastChr) != 0)
		{
			if(strcmp(strLastChr, "") != 0)
			{
				if(nProbeNum != (pBARData->vSeqData[nSeqNum]->nDataNum-1))
				{
					printf("Error: HTS_Aln2BAR, cannot read input file correctly!");
					exit(EXIT_FAILURE);
				}
			}

			nSeqNum++;
			strcpy(strLastChr, strChr);
			nProbeNum = 0;
		}
		else
		{
			nProbeNum++;
		}

		chp1 = chp2+1;
		chp2 = strchr(chp1, '\t');
		nk = 0;
		while(chp2 != NULL)
		{
			if(nk >= pBARData->vSeqData[nSeqNum]->nColNum)
			{
				printf("Error: HTS_Aln2BAR, input file format error, column number inconsistent!");
				exit(EXIT_FAILURE);
			}

			*chp2 = '\0';
			StrTrimLeft(chp1);
			StrTrimRight(chp1);
			pBARData->vSeqData[nSeqNum]->vData[nk]->pMatElement[nProbeNum] = atof(chp1);
			nk++;

			chp1 = chp2+1;
			chp2 = strchr(chp1, '\t');
		}

		if(nk >= pBARData->vSeqData[nSeqNum]->nColNum)
		{
			printf("Error: HTS_Aln2BAR, input file format error, column number inconsistent!");
			exit(EXIT_FAILURE);
		}
		
		StrTrimLeft(chp1);
		StrTrimRight(chp1);
		pBARData->vSeqData[nSeqNum]->vData[nk]->pMatElement[nProbeNum] = atof(chp1);
		nk++;

		if(nk != pBARData->vSeqData[nSeqNum]->nColNum)
		{
			printf("Error: HTS_Aln2BAR, input file format error, column number inconsistent!");
			exit(EXIT_FAILURE);
		}
	}

	
	/* close file */
	fclose(fpIn);

	if(nProbeNum != (pBARData->vSeqData[nSeqNum]->nDataNum-1))
	{
		printf("Error: HTS_Aln2BAR, cannot read input file correctly!");
		exit(EXIT_FAILURE);
	}

	if(nSeqNum != pBARData->nSeqNum-1)
	{
		printf("Error: HTS_Aln2BAR, cannot read input file correctly!");
		exit(EXIT_FAILURE);
	}
	
	/* write to BAR file */
	pCol = NULL;
	pCol = CreateIntMatrix(1,pBARData->nColNum);
	if(pCol == NULL)
	{
		printf("Error: HTS_Aln2BAR, cannot create output column information!\n");
		exit(EXIT_FAILURE);
	}
	
	pCol->pMatElement[0] = 1;
	for(ni=0; ni<nColNum; ni++)
	{
		pCol->pMatElement[1+ni] = 1;
		Affy_SaveBAR_Columns_Fast(strBARFile, pBARData, pCol);
		pCol->pMatElement[1+ni] = 0;
	}
	DestroyIntMatrix(pCol);

	
	/* release memory */
	DeleteString(pSeqInfo);
	pSeqInfo = NULL;

	Affy_BARData_Destroy(&pBARData);

	RemoveFiles(strBARTmpFileName);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2BARv2()                                                        */
/*  Convert high throughput sequencing alignment to bar file.              */
/* ----------------------------------------------------------------------- */ 
int HTS_Aln2BARv2(char strTXTFile[], char strBARFile[])
{
	/* define */
	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));

	struct tagBARData *pBARData = NULL;
	FILE *fpIn;
	FILE *fpOut;

	/* variables */
	int nFieldNum = 0;
	int nSeqNum = 0;
	int nColNum = 0;
	struct tagString *pSeqInfo = NULL;
	int nProbeNum = 0;
	char chStrand;
	int nTotalProbeNum = 0;
	char strLine[LONG_LINE_LENGTH];
	char strTemp[LONG_LINE_LENGTH];
	char *chp1,*chp2;
	char strChr[MED_LINE_LENGTH];
	char strLastChr[MED_LINE_LENGTH];
	int nLastPos;
	int nPos,nAlnCount;
	int nAlnCountF = 0;
	int nAlnCountR = 0;
	char strBARTmpFileName[MED_LINE_LENGTH];
	struct INTMATRIX *pCol = NULL;
	
	/* count */
	int ni,nj,nk;

	/* load */
	fpIn = NULL;
	fpIn = fopen(strTXTFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: HTS_Aln2BAR, cannot open source file!\n");
		exit(EXIT_FAILURE);
	}

	sprintf(strBARTmpFileName, "%s.tmp", strBARFile);
	fpOut = NULL;
	fpOut = fopen(strBARTmpFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: HTS_Aln2BAR, cannot open temporary output file!\n");
		exit(EXIT_FAILURE);
	}

	/* load input */
	pSeqInfo = NULL;
	strcpy(strLastChr, "");
	nLastPos = -1;
	nProbeNum = 0;
	nColNum = 3;
	nSeqNum = 0;
	nFieldNum = 3;

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '0')
			continue;
		if(strLine[0] == '#')
		{
			/* fprintf(fpOut, "%s\n", strLine); */
			continue;
		}

		sscanf(strLine, "%s %d %c", strChr, &nPos, &chStrand);
		
		if(strcmp(strChr, strLastChr) != 0)
		{
			if(strcmp(strLastChr, "") != 0)
			{
				sprintf(strTemp, "%s\t%d\t", strLastChr, nProbeNum);
				StringAddTail(&pSeqInfo, strTemp);
				fprintf(fpOut, "%s\t%d\t%d\t%d\t%d\n", strLastChr, nLastPos, nAlnCount, nAlnCountF, nAlnCountR);
			}

			nSeqNum++;
			strcpy(strLastChr, strChr);
			nLastPos = nPos;
			nProbeNum = 1;
			nAlnCount = 1;
			if( (chStrand == 'R') || (chStrand == '-') || (chStrand == '1') || (chStrand == 'r') )
			{
				nAlnCountF = 0;
				nAlnCountR = 1;
			}
			else
			{
				nAlnCountF = 1;
				nAlnCountR = 0;
			}
		}
		else
		{
			if(nPos != nLastPos)
			{
				fprintf(fpOut, "%s\t%d\t%d\t%d\t%d\n", strLastChr, nLastPos, nAlnCount, nAlnCountF, nAlnCountR);
				nLastPos = nPos;
				nAlnCount = 1;
				if( (chStrand == 'R') || (chStrand == '-') || (chStrand == '1') || (chStrand == 'r') )
				{
					nAlnCountF = 0;
					nAlnCountR = 1;
				}
				else
				{
					nAlnCountF = 1;
					nAlnCountR = 0;
				}
				nProbeNum++;
			}
			else
			{
				nAlnCount += 1;
				if( (chStrand == 'R') || (chStrand == '-') || (chStrand == '1') || (chStrand == 'r') )
				{
					nAlnCountR += 1;
				}
				else
				{
					nAlnCountF += 1;
				}
			}
		}
	}

	if(strcmp(strLastChr, "") != 0)
	{
		sprintf(strTemp, "%s\t%d\t", strLastChr, nProbeNum);
		StringAddTail(&pSeqInfo, strTemp);
		fprintf(fpOut, "%s\t%d\t%d\t%d\t%d\n", strLastChr, nLastPos, nAlnCount, nAlnCountF, nAlnCountR);
	}

	fclose(fpIn);
	fclose(fpOut);


	/* load head */
	if(nSeqNum <= 0)
	{
		printf("Warning: empty sequences!\n");
		return PROC_FAILURE;
	}

	/* create BAR object */
	pBARData = Affy_BARData_Create();
	if(pBARData == NULL)
	{
		printf("Error: HTS_Aln2BAR, BARData object was not created!\n");
		exit(EXIT_FAILURE);
	}
	strcpy(pBARData->strMagicnumber, "barr\r\n\032\n");
	pBARData->fVersionnumber = 2.0;
    pBARData->nSeqNum = nSeqNum;
	pBARData->nColNum = nColNum+1;
	pBARData->nParamNum = 0;
	pBARData->vParamName = NULL;
	pBARData->vParamValue = NULL;
	pBARData->pFieldType = CreateIntMatrix(1, pBARData->nColNum);
	if(pBARData->pFieldType == NULL)
	{
		printf("Error: HTS_Aln2BAR, cannot allocate memory for field type.\n");
		exit(EXIT_FAILURE);
	}
	pBARData->pFieldType->pMatElement[0] = 2;
	for(ni=1; ni<pBARData->nColNum; ni++)
		pBARData->pFieldType->pMatElement[ni] = 1;

	pBARData->vSeqData = (struct tagBARSeq **)calloc(pBARData->nSeqNum, sizeof(struct tagBARSeq *));
	if(pBARData->vSeqData == NULL)
	{
		printf("Error: TileMapv2_TXT2BAR, cannot allocate memory for loading sequence data.\n");
		exit(EXIT_FAILURE);
	}

	chp1 = pSeqInfo->m_pString;
	for(ni=0; ni<pBARData->nSeqNum; ni++)
	{
		/* create BARSeq object */
		pBARData->vSeqData[ni] = Affy_BARSeq_Create();
		if(pBARData->vSeqData[ni] == NULL)
		{
			printf("Error: TileMapv2_TXT2BAR, cannot create BARSeq object.\n");
			exit(EXIT_FAILURE);
		}

		pBARData->vSeqData[ni]->pSeqGroupName = NULL;
		pBARData->vSeqData[ni]->pSeqVersion = NULL;
		pBARData->vSeqData[ni]->nParamNum = 0;
		pBARData->vSeqData[ni]->vParamName = NULL;
		pBARData->vSeqData[ni]->vParamValue = NULL;

		pBARData->vSeqData[ni]->nColNum = pBARData->nColNum;

		chp2 = strchr(chp1, '\t');
		*chp2 = '\0';
		strcpy(strTemp, chp1);
		StringAddTail(&(pBARData->vSeqData[ni]->pSeqName), strTemp);
		chp1 = chp2+1;

		chp2 = strchr(chp1, '\t');
		*chp2 = '\0';
		strcpy(strTemp, chp1);
		pBARData->vSeqData[ni]->nDataNum = atoi(strTemp);
		chp1 = chp2+1;


		pBARData->vSeqData[ni]->vData = (struct DOUBLEMATRIX **)calloc(pBARData->vSeqData[ni]->nColNum, sizeof(struct DOUBLEMATRIX *));
		if(pBARData->vSeqData[ni]->vData == NULL)
		{
			printf("Error: HTS_Aln2BAR, cannot allocate memory for loading bpmap coordinate data!\n");
			exit(EXIT_FAILURE);
		}
		if(pBARData->vSeqData[ni]->nDataNum > 0)
		{
			for(nj=0; nj<pBARData->vSeqData[ni]->nColNum; nj++)
			{
				pBARData->vSeqData[ni]->vData[nj] = CreateDoubleMatrix(1, pBARData->vSeqData[ni]->nDataNum );
				if(pBARData->vSeqData[ni]->vData[nj] == NULL)
				{
					printf("Error: HTS_Aln2BAR, cannot allocate memory for loading bpmap coordinate data!\n");
					exit(EXIT_FAILURE);
				}
			}
		}
	}
	
	/* load data */
	fpIn = NULL;
	fpIn = fopen(strBARTmpFileName, "r");
	if(fpIn == NULL)
	{
		printf("Error: HTS_Aln2BAR, cannot open source file!\n");
		exit(EXIT_FAILURE);
	}

	strcpy(strLastChr, "");
	nSeqNum = -1;
	nProbeNum = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '0')
			continue;
		if(strLine[0] == '#')
			continue;

		chp1 = strLine;
		chp2 = strchr(chp1, '\t');
		if(chp2 == NULL)
		{
			printf("Error: HTS_Aln2BAR, wrong input file format!\n");
			exit(EXIT_FAILURE);
		}
		*chp2 = '\0';
		StrTrimRight(chp1);
		strcpy(strChr, chp1);

		if(strcmp(strChr, strLastChr) != 0)
		{
			if(strcmp(strLastChr, "") != 0)
			{
				if(nProbeNum != (pBARData->vSeqData[nSeqNum]->nDataNum-1))
				{
					printf("Error: HTS_Aln2BAR, cannot read input file correctly!");
					exit(EXIT_FAILURE);
				}
			}

			nSeqNum++;
			strcpy(strLastChr, strChr);
			nProbeNum = 0;
		}
		else
		{
			nProbeNum++;
		}

		chp1 = chp2+1;
		chp2 = strchr(chp1, '\t');
		nk = 0;
		while(chp2 != NULL)
		{
			if(nk >= pBARData->vSeqData[nSeqNum]->nColNum)
			{
				printf("Error: HTS_Aln2BAR, input file format error, column number inconsistent!");
				exit(EXIT_FAILURE);
			}

			*chp2 = '\0';
			StrTrimLeft(chp1);
			StrTrimRight(chp1);
			pBARData->vSeqData[nSeqNum]->vData[nk]->pMatElement[nProbeNum] = atof(chp1);
			nk++;

			chp1 = chp2+1;
			chp2 = strchr(chp1, '\t');
		}

		if(nk >= pBARData->vSeqData[nSeqNum]->nColNum)
		{
			printf("Error: HTS_Aln2BAR, input file format error, column number inconsistent!");
			exit(EXIT_FAILURE);
		}
		
		StrTrimLeft(chp1);
		StrTrimRight(chp1);
		pBARData->vSeqData[nSeqNum]->vData[nk]->pMatElement[nProbeNum] = atof(chp1);
		nk++;

		if(nk != pBARData->vSeqData[nSeqNum]->nColNum)
		{
			printf("Error: HTS_Aln2BAR, input file format error, column number inconsistent!");
			exit(EXIT_FAILURE);
		}
	}

	
	/* close file */
	fclose(fpIn);

	if(nProbeNum != (pBARData->vSeqData[nSeqNum]->nDataNum-1))
	{
		printf("Error: HTS_Aln2BAR, cannot read input file correctly!");
		exit(EXIT_FAILURE);
	}

	if(nSeqNum != pBARData->nSeqNum-1)
	{
		printf("Error: HTS_Aln2BAR, cannot read input file correctly!");
		exit(EXIT_FAILURE);
	}
	
	/* write to BAR file */
	pCol = NULL;
	pCol = CreateIntMatrix(1,pBARData->nColNum);
	if(pCol == NULL)
	{
		printf("Error: HTS_Aln2BAR, cannot create output column information!\n");
		exit(EXIT_FAILURE);
	}
	
	pCol->pMatElement[0] = 1;
	pCol->pMatElement[1] = 1;
	Affy_SaveBAR_Columns_Fast(strBARFile, pBARData, pCol);
	pCol->pMatElement[1] = 0;
	pCol->pMatElement[2] = 1;
	sprintf(strTemp, "%s_F.bar", strBARFile);
	Affy_SaveBAR_Columns_Fast(strTemp, pBARData, pCol);
	pCol->pMatElement[2] = 0;
	pCol->pMatElement[3] = 1;
	sprintf(strTemp, "%s_R.bar", strBARFile);
	Affy_SaveBAR_Columns_Fast(strTemp, pBARData, pCol);
	DestroyIntMatrix(pCol);

	
	/* release memory */
	DeleteString(pSeqInfo);
	pSeqInfo = NULL;

	Affy_BARData_Destroy(&pBARData);

	RemoveFiles(strBARTmpFileName);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_WindowSummary()                                                    */
/*  Summarize window counts.                                               */
/* ----------------------------------------------------------------------- */ 
int HTS_WindowSummary(char strBARFile[], char strChrList[], char strChrLen[], 
					  int nW, char strOutFile[], struct tagBARData *pRepeatData)
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;

	struct tagBARData *pBARData = NULL;
	struct INTMATRIX *pChrLen = NULL;
	struct DOUBLEMATRIX *pCount = NULL;
	int nMaxC = 0;
	int nUpperC = 4096;
	int *vC;
	double dTotal;
	
	char strLine[LONG_LINE_LENGTH];
	int ni,nj,nk,nl,nx,nLen;

	double dR0,dR1,dR2;
	double dlambda,dpoisp,dalpha,dbeta,dnegp;
	double dP1,dP2,dTemp;

	/* initial check */
	if(nW <= 0)
	{
		nW = 100;
		printf("Warning: Window size<=0! Proceed with default window size = 100\n");
	}

	/* load bar data */
	pBARData = NULL;
	pBARData = Affy_LoadBAR_Fast(strBARFile);
	if(pBARData == NULL)
	{
		printf("Error: HTS_WindowSummary, cannot load raw data!\n");
		exit(EXIT_FAILURE);
	}

	/* load chromosome length */
	pChrLen = IMLOAD(strChrLen);
	if(pChrLen == NULL)
	{
		printf("Error: HTS_WindowSummary, cannot load chromosome length!\n");
		exit(EXIT_FAILURE);
	}

	/* create count matrix */
	pCount = CreateDoubleMatrix(nUpperC, 1);
	if(pCount == NULL)
	{
		printf("Error: HTS_WindowSummary, cannot create count matrix!\n");
		exit(EXIT_FAILURE);
	}

	/* process chromosome one by one */
	fpIn = NULL;
	fpIn = fopen(strChrList, "r");
	if(fpIn == NULL)
	{
		printf("Error: HTS_WindowSummary, cannot open chrlist file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	nMaxC = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		nLen = pChrLen->pMatElement[ni]/nW;
		if(pChrLen->pMatElement[ni]%nW != 0)
			nLen += 1;

		/* prepare memory */
		vC = NULL;
		vC = (int *)calloc(nLen, sizeof(int));
		if(vC == NULL)
		{
			printf("Error: HTS_WindowSummary, cannot create memory for counting!\n");
			exit(EXIT_FAILURE);
		}

		/* find matching chromosome */
		if(pRepeatData == NULL)
		{
			nx = -1;
		}
		else
		{
			for(nx=0; nx<pRepeatData->nSeqNum; nx++)
			{
				if(strcmp(strLine, pRepeatData->vSeqData[nx]->pSeqName->m_pString) == 0)
					break;
			}

			if(nx >= pRepeatData->nSeqNum)
			{
				printf("Error: repeat masking file does not provide complete genome coverage!\n");
				exit(EXIT_FAILURE);
			}

			if(nLen != pRepeatData->vSeqData[nx]->nDataNum)
			{
				printf("Error: the length of repeat masking file does not match the genome length!\n");
				exit(EXIT_FAILURE);
			}
		}

		for(nj=0; nj<pBARData->nSeqNum; nj++)
		{
			if(strcmp(strLine, pBARData->vSeqData[nj]->pSeqName->m_pString) != 0)
				continue;

			for(nk=0; nk<pBARData->vSeqData[nj]->nDataNum; nk++)
			{
				nl = (int)(pBARData->vSeqData[nj]->vData[0]->pMatElement[nk])/nW;
				if(nl == nLen)
					nl = nLen-1;
				else if(nl > nLen)
				{
					printf("Error: HTS_WindowSummary, index out of range! Please check if the genomes are matching\n");
					exit(EXIT_FAILURE);
				}

				vC[nl] += (int)(pBARData->vSeqData[nj]->vData[1]->pMatElement[nk]);
			}
		}


		/* count */
		for(nj=0; nj<nLen; nj++)
		{
			nk = vC[nj];
			if(pRepeatData != NULL)
			{
				if(pRepeatData->vSeqData[nx]->vData[1]->pMatElement[nj] == 1)
					continue;
			}


			if(nk >= nUpperC)
				nk = nUpperC-1;

			if(nk>nMaxC)
				nMaxC = nk;

			pCount->pMatElement[nk] += 1;
		}

		/* free memory */
		free(vC);

		ni++;
	}
	
	fclose(fpIn);


	/* output summary statistics */
	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: HTS_WindowSummary, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	dTotal = 0.0;
	for(ni=0; ni<=nMaxC; ni++)
	{
		dTotal += pCount->pMatElement[ni];
	}

	/* estimate poisson and neg-binonmial */
	dR0 = pCount->pMatElement[0]/dTotal;
	dR1 = pCount->pMatElement[1]/dTotal;
	dR2 = pCount->pMatElement[2]/dTotal;
	dR2 = dR2/dR1;
	dR1 = dR1/dR0;

	dlambda = dR1;
	dpoisp = dR0/exp(-dlambda);
	if(dpoisp > 1.0)
		dpoisp = 1.0;

	dalpha = dR1/(2.0*dR2-dR1);
	dbeta = 1.0/(2.0*dR2-dR1)-1.0;
	dnegp = dR0/pow( (dbeta/(dbeta+1.0)), dalpha);
	if(dnegp > 1.0)
		dnegp = 1.0;


	fprintf(fpOut,"# Window_Size=%d\tPoisson_Lambda=%f\tPoisson_p=%f\tNegBinomial_Alpha=%f\tNegBinomial_Beta=%f\tNegBinomial_p=%f\n", nW, dlambda, dpoisp, dalpha, dbeta, dnegp);
	fprintf(fpOut, "#No_of_reads/window\tNo_of_window\tpercentage\tpoisson_expected\tpoisson_exp/obs\tnegbinomial_expected\tnegbinomial_exp/obs\n");
	for(ni=0; ni<=nMaxC; ni++)
	{
		dP1 = 0.0;
		dP1 = ni*log(dlambda)-dlambda-gammaln(ni+1.0)+log(dpoisp);
		dP1 = exp(dP1);

		dP2 = 0.0;
		dP2 = gammaln(ni+dalpha)-gammaln(ni+1.0)-gammaln(dalpha)+dalpha*log(dbeta/(dbeta+1.0))-ni*log(dbeta+1.0)+log(dnegp);
		dP2 = exp(dP2);

		dTemp = pCount->pMatElement[ni]/dTotal;
		
		fprintf(fpOut, "%d\t%d\t%f\t%f\t%f\t%f\t%f\n", ni, (int)(pCount->pMatElement[ni]), dTemp, dP1, dP1/(dTemp+1e-20), dP2, dP2/(dTemp+1e-20));
	}

	fclose(fpOut);

	/* destroy memory */
	Affy_BARData_Destroy(&pBARData);
	DestroyIntMatrix(pChrLen);
	DestroyDoubleMatrix(pCount);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_WindowSummaryv2()                                                  */
/*  Summarize window counts.                                               */
/* ----------------------------------------------------------------------- */ 
int HTS_WindowSummaryv2(char strBARFile[], char strChrList[], char strChrLen[], 
					  int nW, char strOutFile[], int nCombineShift, char strRepeatFile[])
{
	/* define */
	char strFileName[MED_LINE_LENGTH];
	char strFileName1[MED_LINE_LENGTH];
	char strFileName2[MED_LINE_LENGTH];
	char strFileName3[MED_LINE_LENGTH];
	FILE *fpIn;
	FILE *fpOut;
	char strLine[LONG_LINE_LENGTH];
	struct tagBARData *pRepeatData = NULL;

	/* run one by one */
	sprintf(strFileName1, "%s.tmp1", strOutFile);
	sprintf(strFileName2, "%s.tmp2", strOutFile);
	sprintf(strFileName3, "%s.tmp3", strOutFile);
	if(strcmp(strRepeatFile, "") != 0)
	{
		pRepeatData = Affy_LoadBAR_Fast(strRepeatFile);
		if(pRepeatData == NULL)
		{
			printf("Error: cannot load repeat masking file!\n");
			exit(EXIT_FAILURE);
		}
	}

	if(nCombineShift == 1)
	{
		sprintf(strFileName, "%s_C.bar", strBARFile);
		HTS_WindowSummary(strFileName, strChrList, strChrLen, nW, strFileName1, pRepeatData);
	}
	else
	{
		HTS_WindowSummary(strBARFile, strChrList, strChrLen, nW, strFileName1, pRepeatData);
	}
	
	sprintf(strFileName, "%s_F.bar", strBARFile);
	HTS_WindowSummary(strFileName, strChrList, strChrLen, nW, strFileName2, pRepeatData);
	sprintf(strFileName, "%s_R.bar", strBARFile);
	HTS_WindowSummary(strFileName, strChrList, strChrLen, nW, strFileName3, pRepeatData);

	Affy_BARData_Destroy(&pRepeatData);

	/* summarize results */
	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: HTS_WindowSummaryv2, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "#Forward+Reverse_Combined\n");
	fpIn = NULL;
	fpIn = fopen(strFileName1, "r");
	if(fpIn == NULL)
	{
		printf("Error: HTS_WindowSummaryv2, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);

		fprintf(fpOut, "%s\n", strLine);
	}

	fclose(fpIn);

	fprintf(fpOut, "\n#Forward_Only\n");
	fpIn = NULL;
	fpIn = fopen(strFileName2, "r");
	if(fpIn == NULL)
	{
		printf("Error: HTS_WindowSummaryv2, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);

		fprintf(fpOut, "%s\n", strLine);
	}

	fclose(fpIn);

	fprintf(fpOut, "\n#Reverse_Only\n");
	fpIn = NULL;
	fpIn = fopen(strFileName3, "r");
	if(fpIn == NULL)
	{
		printf("Error: HTS_WindowSummaryv2, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);

		fprintf(fpOut, "%s\n", strLine);
	}

	fclose(fpIn);

	fclose(fpOut);

	/* clear temp files */
	RemoveFiles(strFileName1);
	RemoveFiles(strFileName2);
	RemoveFiles(strFileName3);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_OneSample_Main()                                            */
/*  Find enriched region from a sequencing data set.                       */
/* ----------------------------------------------------------------------- */ 
int HTS_Enrich_OneSample_Main(char strBARFile[], int nW, int nS, int nCutoff,
					  int nMinLen, int nMaxGap,
					  char strExportFolder[], char strOutFileTitle[])
{
	/* define */
	struct tagBARData *pBARData = NULL;
	struct DOUBLEMATRIX *pRegion0;
	struct DOUBLEMATRIX *pRegion;
	char strFileName[MED_LINE_LENGTH];
	struct INTMATRIX *pType = NULL;
	struct INTMATRIX *pPriority = NULL;
	struct DOUBLEMATRIX *pRegionSort = NULL;
	struct LONGMATRIX *pRegionSid = NULL;

	/* initial check */
	AdjustDirectoryPath(strExportFolder);
	
	/* load bar data */
	pBARData = NULL;
	pBARData = Affy_LoadBAR_Fast(strBARFile);
	if(pBARData == NULL)
	{
		printf("Error: HTS_Enrich_OneSample_Main, cannot load raw data!\n");
		exit(EXIT_FAILURE);
	}

	/* initial region call */
	pRegion0 = NULL;
	pRegion0 = HTS_Enrich_OneSample_CallRegion_Initial(pBARData, nW, nS, nCutoff,
					  strExportFolder, strOutFileTitle);

	/* merge and filter regions */
	pRegion = HTS_Enrich_OneSample_MergeRegion(pRegion0, nMinLen, nMaxGap);
	DestroyDoubleMatrix(pRegion0);
	

	/* collect region information */
	HTS_Enrich_OneSample_RegionCollectInfo(pRegion, pBARData);
	
	/* sort regions */
	pType = NULL;
	pType = CreateIntMatrix(1, pRegion->nWidth);
	if(pType == NULL)
	{
		printf("Error: HTS_Enrich_OneSample_Main, cannot allocate memory for summarizing results!\n");
		exit(EXIT_FAILURE);
	}
	pType->pMatElement[0] = 2;
	pType->pMatElement[1] = 2;
	pType->pMatElement[2] = 2;
	pType->pMatElement[3] = 2;
	pType->pMatElement[4] = 2;
	pType->pMatElement[5] = 2;
	pType->pMatElement[6] = 2;

	pPriority = NULL;
	pPriority = CreateIntMatrix(1, 3);
	if(pPriority == NULL)
	{
		printf("Error: HTS_Enrich_OneSample_Main, cannot allocate memory for summarizing results!\n");
		exit(EXIT_FAILURE);
	}
	pPriority->pMatElement[0] = 4;
	pPriority->pMatElement[1] = 6;
	pPriority->pMatElement[2] = 3;

	pRegionSort = NULL;
	pRegionSid = NULL;
	DMSORTROWS(pRegion, pType, pPriority, &pRegionSort, &pRegionSid);
	
	DestroyIntMatrix(pType);
	DestroyIntMatrix(pPriority);
	DestroyDoubleMatrix(pRegion);
	
	/* export result */
	/* rank, chr, start, end, strand, length, biggest n, pos of biggest n, total reads */
	sprintf(strFileName, "%s%s.cod", strExportFolder, strOutFileTitle);
	HTS_Enrich_OneSample_ExportResults(pRegionSort, pBARData, strFileName);
	DestroyDoubleMatrix(pRegionSort);
	DestroyLongMatrix(pRegionSid);


	/* destroy memory */
	Affy_BARData_Destroy(&pBARData);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_OneSamplev2_Main()                                          */
/*  Find enriched region from a sequencing data set.                       */
/* ----------------------------------------------------------------------- */ 
int HTS_Enrich_OneSamplev2_Main_Old(char strBARFile[], int nW, int nS, int nCutoff,
					  int nCutoffF, int nCutoffR, int nMinLen, int nMaxGap,
					  char strExportFolder[], char strOutFileTitle[])
{
	/* define */
	struct tagBARData *pBARData = NULL;
	struct tagBARData *pBARDataF = NULL;
	struct tagBARData *pBARDataR = NULL;
	struct DOUBLEMATRIX *pRegion0;
	struct DOUBLEMATRIX *pRegion;
	struct DOUBLEMATRIX *pNewRegion;
	struct DOUBLEMATRIX *pRegionF;
	struct DOUBLEMATRIX *pRegionR;
	char strFileName[MED_LINE_LENGTH];
	char strOutFileName[MED_LINE_LENGTH];
	struct INTMATRIX *pType = NULL;
	struct INTMATRIX *pPriority = NULL;
	struct DOUBLEMATRIX *pRegionSort = NULL;
	struct LONGMATRIX *pRegionSid = NULL;
	int ni;

	/* initial check */
	AdjustDirectoryPath(strExportFolder);

	
	/* load bar data */
	pBARData = NULL;
	pBARData = Affy_LoadBAR_Fast(strBARFile);
	if(pBARData == NULL)
	{
		printf("Error: HTS_Enrich_OneSamplev2_Main, cannot load raw data!\n");
		exit(EXIT_FAILURE);
	}

	/* initial region call */
	pRegion0 = NULL;
	pRegion0 = HTS_Enrich_OneSample_CallRegion_Initial(pBARData, nW, nS, nCutoff,
					  strExportFolder, strOutFileTitle);

	/* merge and filter regions */
	pRegion = HTS_Enrich_OneSample_MergeRegion(pRegion0, nMinLen, nMaxGap);
	DestroyDoubleMatrix(pRegion0);

	/* load forward */
	sprintf(strFileName, "%s_F.bar", strBARFile);
	sprintf(strOutFileName, "%s_F", strOutFileTitle);
	pBARDataF = NULL;
	pBARDataF = Affy_LoadBAR_Fast(strFileName);
	if(pBARDataF == NULL)
	{
		printf("Error: HTS_Enrich_OneSamplev2_Main, cannot load raw data!\n");
		exit(EXIT_FAILURE);
	}

	/* initial region call */
	pRegion0 = NULL;
	pRegion0 = HTS_Enrich_OneSample_CallRegion_Initial(pBARDataF, nW, nS, nCutoffF,
					  strExportFolder, strOutFileName);

	/* merge and filter regions */
	pRegionF = HTS_Enrich_OneSample_MergeRegion(pRegion0, nMinLen, nMaxGap);
	DestroyDoubleMatrix(pRegion0);

	/* load reverse */
	sprintf(strFileName, "%s_R.bar", strBARFile);
	sprintf(strOutFileName, "%s_R", strOutFileTitle);
	pBARDataR = NULL;
	pBARDataR = Affy_LoadBAR_Fast(strFileName);
	if(pBARDataR == NULL)
	{
		printf("Error: HTS_Enrich_OneSamplev2_Main, cannot load raw data!\n");
		exit(EXIT_FAILURE);
	}

	/* initial region call */
	pRegion0 = NULL;
	pRegion0 = HTS_Enrich_OneSample_CallRegion_Initial(pBARDataR, nW, nS, nCutoffR,
					  strExportFolder, strOutFileName);

	/* merge and filter regions */
	pRegionR = HTS_Enrich_OneSample_MergeRegion(pRegion0, nMinLen, nMaxGap);
	DestroyDoubleMatrix(pRegion0);

	/* collect region information */
	HTS_Enrich_OneSample_RegionCollectInfo(pRegion, pBARData);
	HTS_Enrich_OneSample_RegionCollectInfo(pRegionF, pBARDataF);
	HTS_Enrich_OneSample_RegionCollectInfo(pRegionR, pBARDataR);
	
	/* match forward and reverse regions */
	pNewRegion = pRegion;
	pRegion = NULL;
	pRegion = HTS_Enrich_OneSamplev2_MatchFRRegions(pNewRegion, pRegionF, pRegionR);
	DestroyDoubleMatrix(pNewRegion);
	
	/* sort regions */
	pType = NULL;
	pType = CreateIntMatrix(1, pRegion->nWidth);
	if(pType == NULL)
	{
		printf("Error: HTS_Enrich_OneSample_Main, cannot allocate memory for summarizing results!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<pType->nWidth; ni++)
	{
		IMSETAT(pType, 0, ni, 2);
	}

	pPriority = NULL;
	pPriority = CreateIntMatrix(1, 3);
	if(pPriority == NULL)
	{
		printf("Error: HTS_Enrich_OneSample_Main, cannot allocate memory for summarizing results!\n");
		exit(EXIT_FAILURE);
	}
	pPriority->pMatElement[0] = 4;
	pPriority->pMatElement[1] = 6;
	pPriority->pMatElement[2] = 3;

	pRegionSort = NULL;
	pRegionSid = NULL;
	DMSORTROWS(pRegion, pType, pPriority, &pRegionSort, &pRegionSid);
	
	DestroyIntMatrix(pType);
	DestroyIntMatrix(pPriority);
	DestroyDoubleMatrix(pRegion);
	DestroyDoubleMatrix(pRegionF);
	DestroyDoubleMatrix(pRegionR);
	
	/* export result */
	/* rank, chr, start, end, strand, length, biggest n, pos of biggest n, total reads */
	sprintf(strFileName, "%s%s.cod", strExportFolder, strOutFileTitle);
	HTS_Enrich_OneSamplev2_ExportResults(pRegionSort, pBARData, strFileName, 0, 0, 0, 0, 0);
	DestroyDoubleMatrix(pRegionSort);
	DestroyLongMatrix(pRegionSid);


	/* destroy memory */
	Affy_BARData_Destroy(&pBARData);
	Affy_BARData_Destroy(&pBARDataF);
	Affy_BARData_Destroy(&pBARDataR);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_OneSamplev2_Main()                                          */
/*  Find enriched region from a sequencing data set.                       */
/* ----------------------------------------------------------------------- */ 
int HTS_Enrich_OneSamplev2_Main(char strBARFile[], int nW, int nS, int nCutoff, 
					  int nCutoffF, int nCutoffR, int nMinLen, int nMaxGap, 
					  char strExportFolder[], char strOutFileTitle[],
					  int nBR, int nBRL, int nSSF, int nCombineShift)
{
	/* define */
	struct tagBARData *pBARData = NULL;
	struct tagBARData *pBARDataF = NULL;
	struct tagBARData *pBARDataR = NULL;
	struct DOUBLEMATRIX *pRegion0;
	struct DOUBLEMATRIX *pRegion;
	struct DOUBLEMATRIX *pNewRegion;
	char strFileName[MED_LINE_LENGTH];
	char strOutFileName[MED_LINE_LENGTH];
	struct INTMATRIX *pType = NULL;
	struct INTMATRIX *pPriority = NULL;
	struct DOUBLEMATRIX *pRegionSort = NULL;
	struct LONGMATRIX *pRegionSid = NULL;
	int ni;

	/* initial check */
	AdjustDirectoryPath(strExportFolder);

	/* load bar data */
	pBARData = NULL;
	if(nCombineShift == 1)
	{
		sprintf(strFileName, "%s_C.bar", strBARFile);
		pBARData = Affy_LoadBAR_Fast(strFileName);
	}
	else
	{
		pBARData = Affy_LoadBAR_Fast(strBARFile);
	}
	if(pBARData == NULL)
	{
		printf("Error: HTS_Enrich_OneSamplev2_Main, cannot load raw data!\n");
		exit(EXIT_FAILURE);
	}

	/* initial region call */
	pRegion0 = NULL;
	pRegion0 = HTS_Enrich_OneSample_CallRegion_Initial(pBARData, nW, nS, nCutoff,
					  strExportFolder, strOutFileTitle);

	/* merge and filter regions */
	pRegion = HTS_Enrich_OneSample_MergeRegion(pRegion0, nMinLen, nMaxGap);
	DestroyDoubleMatrix(pRegion0);

	/* load forward */
	sprintf(strFileName, "%s_F.bar", strBARFile);
	sprintf(strOutFileName, "%s_F", strOutFileTitle);
	pBARDataF = NULL;
	pBARDataF = Affy_LoadBAR_Fast(strFileName);
	if(pBARDataF == NULL)
	{
		printf("Error: HTS_Enrich_OneSamplev2_Main, cannot load raw data!\n");
		exit(EXIT_FAILURE);
	}

	/* initial region call */
	pRegion0 = NULL;
	pRegion0 = HTS_Enrich_OneSample_CallRegion_Initial(pBARDataF, nW, nS, nCutoffF,
					  strExportFolder, strOutFileName);
	DestroyDoubleMatrix(pRegion0);


	/* load reverse */
	sprintf(strFileName, "%s_R.bar", strBARFile);
	sprintf(strOutFileName, "%s_R", strOutFileTitle);
	pBARDataR = NULL;
	pBARDataR = Affy_LoadBAR_Fast(strFileName);
	if(pBARDataR == NULL)
	{
		printf("Error: HTS_Enrich_OneSamplev2_Main, cannot load raw data!\n");
		exit(EXIT_FAILURE);
	}

	/* initial region call */
	pRegion0 = NULL;
	pRegion0 = HTS_Enrich_OneSample_CallRegion_Initial(pBARDataR, nW, nS, nCutoffR,
					  strExportFolder, strOutFileName);
	DestroyDoubleMatrix(pRegion0);

	/* collect region information */
	HTS_Enrich_OneSample_RegionCollectInfo(pRegion, pBARData);
	pNewRegion = pRegion;
	pRegion = NULL;
	pRegion = HTS_Enrich_OneSamplev2_RegionCollectInfo(pNewRegion, pBARDataF, nW);
	DestroyDoubleMatrix(pNewRegion);
	pNewRegion = NULL;
	pNewRegion = pRegion;
	pRegion = HTS_Enrich_OneSamplev2_RegionCollectInfo(pNewRegion, pBARDataR, nW);
	DestroyDoubleMatrix(pNewRegion);

	/* sort regions */
	pType = NULL;
	pType = CreateIntMatrix(1, pRegion->nWidth);
	if(pType == NULL)
	{
		printf("Error: HTS_Enrich_OneSample_Main, cannot allocate memory for summarizing results!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<pType->nWidth; ni++)
	{
		IMSETAT(pType, 0, ni, 2);
	}

	pPriority = NULL;
	pPriority = CreateIntMatrix(1, 3);
	if(pPriority == NULL)
	{
		printf("Error: HTS_Enrich_OneSample_Main, cannot allocate memory for summarizing results!\n");
		exit(EXIT_FAILURE);
	}
	pPriority->pMatElement[0] = 4;
	pPriority->pMatElement[1] = 6;
	pPriority->pMatElement[2] = 3;

	pRegionSort = NULL;
	pRegionSid = NULL;
	DMSORTROWS(pRegion, pType, pPriority, &pRegionSort, &pRegionSid);
	
	DestroyIntMatrix(pType);
	DestroyIntMatrix(pPriority);
	DestroyDoubleMatrix(pRegion);
	
	/* export result */
	/* rank, chr, start, end, strand, length, biggest n, pos of biggest n, total reads */
	sprintf(strFileName, "%s%s.cod", strExportFolder, strOutFileTitle);
	HTS_Enrich_OneSamplev2_ExportResults(pRegionSort, pBARData, strFileName, 
		nBR, nBRL, nSSF, nCutoffF, nCutoffR);
	DestroyDoubleMatrix(pRegionSort);
	DestroyLongMatrix(pRegionSid);


	/* destroy memory */
	Affy_BARData_Destroy(&pBARData);
	Affy_BARData_Destroy(&pBARDataF);
	Affy_BARData_Destroy(&pBARDataR);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_OneSample_CallRegion_Initial()                              */
/*  Search for enriched windows.                                           */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *HTS_Enrich_OneSample_CallRegion_Initial(struct tagBARData *pBARData,
			int nW, int nS, int nCutoff, char strExportFolder[], char strOutFileTitle[])
{
	/* define */
	FILE *fpOut;
	FILE *fpReg;
	char strWinBarTmpFile[MED_LINE_LENGTH];
	char strRegTmpFile[MED_LINE_LENGTH];
	char strFileName[MED_LINE_LENGTH];
	struct DOUBLEMATRIX *pRegion = NULL;
	int ni,nj,nk,nl,nx;
	int nP1,nP2;
	int nN;
	int nW2 = nW/2;
	int nStart,nEnd,nMaxN,nMaxNPos,nMaxNPos2;

	/* init */
	if(pBARData == NULL)
		return NULL;

	/* process one by one */
	sprintf(strWinBarTmpFile, "%s%s.bar.tmp", strExportFolder, strOutFileTitle);
	fpOut = NULL;
	fpOut = fopen(strWinBarTmpFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: HTS_Enrich_OneSample_CallRegion_Initial, cannot open temporary file to write display information!\n");
		exit(EXIT_FAILURE);
	}
	fprintf(fpOut, "#chr\tpos\t%s\n", strOutFileTitle);
	fprintf(fpOut, "#chr\tpos\t1\n");

	sprintf(strRegTmpFile, "%s%s.regtmp", strExportFolder, strOutFileTitle);
	fpReg = NULL;
	fpReg = fopen(strRegTmpFile, "w");
	if(fpReg == NULL)
	{
		printf("Error: HTS_Enrich_OneSample_CallRegion_Initial, cannot open temporary file to write region information!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<pBARData->nSeqNum; ni++)
	{
		if(pBARData->vSeqData[ni]->nDataNum <= 0)
			continue;

		nx = pBARData->vSeqData[ni]->nDataNum-1;
		nP1 = (int)(pBARData->vSeqData[ni]->vData[0]->pMatElement[0])-nW;
		nP2 = (int)(pBARData->vSeqData[ni]->vData[0]->pMatElement[nx]);
		nStart = -1;
		nEnd = -1;
		nMaxN = 0;
		nMaxNPos = -1;
		nMaxNPos2 = -1;

		nN = 0;
		nk = 0;
		nl = 0;
		for(nj=nP1; nj<=nP2; nj++)
		{
			for(; nl<pBARData->vSeqData[ni]->nDataNum; nl++)
			{
				if((int)(pBARData->vSeqData[ni]->vData[0]->pMatElement[nl]) >= nj)
					break;
				
				nN -= (int)(pBARData->vSeqData[ni]->vData[1]->pMatElement[nl]);
			}

			for(; nk<pBARData->vSeqData[ni]->nDataNum; nk++)
			{
				if((int)(pBARData->vSeqData[ni]->vData[0]->pMatElement[nk]) >= (nj+nW))
					break;

				nN += (int)(pBARData->vSeqData[ni]->vData[1]->pMatElement[nk]);
			}

			if( (nj+nW2)%nS == 0)
			{
				if( ((nj+nW2)>=0) && (nN>0) )
					fprintf(fpOut, "%s\t%d\t%d\n", pBARData->vSeqData[ni]->pSeqName->m_pString, (int)(nj+nW2), nN);
			}

			if( nN>= nCutoff )
			{
				if(nStart < 0)
				{
					nStart = nj;
					nEnd = nj+nW-1;
					nMaxN = nN;
					nMaxNPos = nj+nW2;
					nMaxNPos2 = nMaxNPos;
				}
				else
				{
					nEnd = nj+nW-1;
					if(nN > nMaxN)
					{
						nMaxN = nN;
						nMaxNPos = nj+nW2;
						nMaxNPos2 = nMaxNPos;
					}
					else if(nN == nMaxN)
					{
						nMaxNPos2 = nj+nW2;
					}
				}
			}
			else
			{
				if(nStart >= 0)
				{
					fprintf(fpReg, "%d\t%d\t%d\t%d\t%d\n", ni, nStart, nEnd, nMaxN, (nMaxNPos+nMaxNPos2)/2);
					nStart = -1;
					nEnd = -1;
					nMaxN = 0;
					nMaxNPos = -1;
					nMaxNPos2 = -1;
				}
			}
		}

		if(nStart >= 0)
		{
			fprintf(fpReg, "%d\t%d\t%d\t%d\t%d\n", ni, nStart, nEnd, nMaxN, (nMaxNPos+nMaxNPos2)/2);
			nStart = -1;
			nEnd = -1;
			nMaxN = 0;
			nMaxNPos = -1;
			nMaxNPos2 = -1;
		}
	}

	fclose(fpOut);
	fclose(fpReg);

	/* load region information */
	pRegion = DMLOAD(strRegTmpFile);

	/* convert temp file to bar file */
	sprintf(strFileName, "%s.cgw", strOutFileTitle);
	TileMapv2_TXT2BAR(strWinBarTmpFile, strExportFolder, strFileName);
	
	/* remove temp file */
	RemoveFiles(strWinBarTmpFile);
	RemoveFiles(strRegTmpFile);

	/* return */
	return pRegion;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_OneSample_MergeRegion()                                     */
/*  Merge and filter regions.                                              */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *HTS_Enrich_OneSample_MergeRegion(struct DOUBLEMATRIX *pRegion0, 
	int nMinLen, int nMaxGap)
{
	/* define */
	struct DOUBLEMATRIX *pMReg = NULL;
	struct DOUBLEMATRIX *pRegion = NULL;
	int nRegNum;
	int ni;
	int nChr0,nChr;
	int nStart0,nStart;
	int nEnd0,nEnd;
	int nMaxN0,nMaxN;
	int nMaxNPos0,nMaxNPos;
	int nIsNew;
	int nLen;

	/* init check */
	if(pRegion0 == NULL)
		return NULL;
	if(pRegion0->nHeight <= 0)
		return NULL;

	/* merge */
	pMReg = CreateDoubleMatrix(pRegion0->nHeight, pRegion0->nWidth+2);
	if(pMReg == NULL)
	{
		printf("Error: HTS_Enrich_OneSample_MergeRegion, cannot create memory for merging regions!\n");
		exit(EXIT_FAILURE);
	}

	nRegNum = 0;
	nChr0 = (int)(DMGETAT(pRegion0, 0, 0));
	nStart0 = (int)(DMGETAT(pRegion0, 0, 1));
	nEnd0 = (int)(DMGETAT(pRegion0, 0, 2));
	nMaxN0 = (int)(DMGETAT(pRegion0, 0, 3));
	nMaxNPos0 = (int)(DMGETAT(pRegion0, 0, 4));

	for(ni=1; ni<pRegion0->nHeight; ni++)
	{
		nIsNew = 0;
		nChr = (int)(DMGETAT(pRegion0, ni, 0));
		nStart = (int)(DMGETAT(pRegion0, ni, 1));
		nEnd = (int)(DMGETAT(pRegion0, ni, 2));
		nMaxN = (int)(DMGETAT(pRegion0, ni, 3));
		nMaxNPos = (int)(DMGETAT(pRegion0, ni, 4));

		if(nChr != nChr0)
		{
			nIsNew = 1;
		}
		else
		{
			if( (nStart-nEnd0) > nMaxGap)
				nIsNew = 1;
		}

		if(nIsNew == 1)
		{
			nLen = nEnd0-nStart0+1;
			if(nLen >= nMinLen)
			{
				DMSETAT(pMReg, nRegNum, 0, (double)nChr0);
				DMSETAT(pMReg, nRegNum, 1, (double)nStart0);
				DMSETAT(pMReg, nRegNum, 2, (double)nEnd0);
				DMSETAT(pMReg, nRegNum, 3, (double)nLen);
				DMSETAT(pMReg, nRegNum, 4, (double)nMaxN0);
				DMSETAT(pMReg, nRegNum, 5, (double)nMaxNPos0);
				nRegNum++;
			}

			nChr0 = nChr;
			nStart0 = nStart;
			nEnd0 = nEnd;
			nMaxN0 = nMaxN;
			nMaxNPos0 = nMaxNPos;
		}
		else
		{
			if(nEnd > nEnd0)
				nEnd0 = nEnd;
			if(nMaxN > nMaxN0)
			{
				nMaxN0 = nMaxN;
				nMaxNPos0 = nMaxNPos;
			}
		}
	}

	nLen = nEnd0-nStart0+1;
	if(nLen >= nMinLen)
	{
		DMSETAT(pMReg, nRegNum, 0, (double)nChr0);
		DMSETAT(pMReg, nRegNum, 1, (double)nStart0);
		DMSETAT(pMReg, nRegNum, 2, (double)nEnd0);
		DMSETAT(pMReg, nRegNum, 3, (double)nLen);
		DMSETAT(pMReg, nRegNum, 4, (double)nMaxN0);
		DMSETAT(pMReg, nRegNum, 5, (double)nMaxNPos0);
		nRegNum++;
	}


	/* prepare the regions */
	pRegion = CreateDoubleMatrix(nRegNum, pMReg->nWidth);
	if(pRegion == NULL)
	{
		printf("Error: HTS_Enrich_OneSample_MergeRegion, cannot create memory for merging regions!\n");
		exit(EXIT_FAILURE);
	}
	memcpy(pRegion->pMatElement, pMReg->pMatElement, (sizeof(double)*nRegNum*pMReg->nWidth));

	DestroyDoubleMatrix(pMReg);

	/* return */
	return pRegion;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_OneSample_MergeRegion()                                     */
/*  Merge and filter regions.                                              */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *HTS_Enrich_OneSamplev2_MatchFRRegions(struct DOUBLEMATRIX *pRegion, 
		struct DOUBLEMATRIX *pRegionF, struct DOUBLEMATRIX *pRegionR)
{
	/* define */
	struct DOUBLEMATRIX *pNewRegion = NULL;
	int ni,nj,nk;
	double dTemp;
	int nChr,nStart,nEnd;

	/* init check */
	if( (pRegion == NULL) || (pRegionF == NULL) || (pRegionR == NULL) )
	{
		printf("Warning: HTS_Enrich_OneSamplev2_MatchFRRegions, empty region!\n");
		return NULL;
	}


	/* create */
	pNewRegion = CreateDoubleMatrix(pRegion->nHeight, pRegion->nWidth+10);
	if(pNewRegion == NULL)
	{
		printf("Error: HTS_Enrich_OneSamplev2_MatchFRRegions, cannot create matrix for linking regions!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<pRegion->nHeight; ni++)
	{
		nChr = (int)(DMGETAT(pRegion, ni, 0));
		nStart = (int)(DMGETAT(pRegion, ni, 1));
		nEnd = (int)(DMGETAT(pRegion, ni, 2));

		for(nj=0; nj<pRegion->nWidth; nj++)
		{
			dTemp = DMGETAT(pRegion, ni, nj);
			DMSETAT(pNewRegion, ni, nj, dTemp);
		}

		nk = HTS_Enrich_OneSamplev2_FindOverlap(pRegionF, nChr, nStart, nEnd);
		if(nk < 0)
		{
			dTemp = -1.0;
			for(; nj<(5+pRegion->nWidth); nj++)
			{
				DMSETAT(pNewRegion, ni, nj, dTemp);
			}
		}
		else
		{
			dTemp = DMGETAT(pRegionF, nk, 1);
			DMSETAT(pNewRegion, ni, nj, dTemp);
			nj++;
			dTemp = DMGETAT(pRegionF, nk, 2);
			DMSETAT(pNewRegion, ni, nj, dTemp);
			nj++;
			dTemp = DMGETAT(pRegionF, nk, 4);
			DMSETAT(pNewRegion, ni, nj, dTemp);
			nj++;
			dTemp = DMGETAT(pRegionF, nk, 5);
			DMSETAT(pNewRegion, ni, nj, dTemp);
			nj++;
			dTemp = DMGETAT(pRegionF, nk, 6);
			DMSETAT(pNewRegion, ni, nj, dTemp);
			nj++;
		}

		nk = HTS_Enrich_OneSamplev2_FindOverlap(pRegionR, nChr, nStart, nEnd);
		if(nk < 0)
		{
			dTemp = -1.0;
			for(; nj<(10+pRegion->nWidth); nj++)
			{
				DMSETAT(pNewRegion, ni, nj, dTemp);
			}
		}
		else
		{
			dTemp = DMGETAT(pRegionR, nk, 1);
			DMSETAT(pNewRegion, ni, nj, dTemp);
			nj++;
			dTemp = DMGETAT(pRegionR, nk, 2);
			DMSETAT(pNewRegion, ni, nj, dTemp);
			nj++;
			dTemp = DMGETAT(pRegionR, nk, 4);
			DMSETAT(pNewRegion, ni, nj, dTemp);
			nj++;
			dTemp = DMGETAT(pRegionR, nk, 5);
			DMSETAT(pNewRegion, ni, nj, dTemp);
			nj++;
			dTemp = DMGETAT(pRegionR, nk, 6);
			DMSETAT(pNewRegion, ni, nj, dTemp);
			nj++;
		}
	}

	/* return */
	return pNewRegion;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_OneSamplev2_FindOverlap()                                   */
/*  Find overlap region.                                                   */
/* ----------------------------------------------------------------------- */ 
int HTS_Enrich_OneSamplev2_FindOverlap(struct DOUBLEMATRIX *pRegion, 
					int nChr, int nStart, int nEnd)
{
	/* define */
	int ni,nj,nk;
	int nIdx = -1;
	int nChr2,nStart2,nEnd2;
	int nMaxN,nN;

	/* init */
	if(pRegion == NULL)
		return -1;

	ni = 0;
	nj = pRegion->nHeight-1;
	while((nj-ni) > 1)
	{
		nk = (nj+ni)/2;
		nChr2 = (int)(DMGETAT(pRegion, nk, 0));
		nStart2 = (int)(DMGETAT(pRegion, nk, 1));
		nEnd2 = (int)(DMGETAT(pRegion, nk, 2));

		if(nChr < nChr2)
		{
			nj = nk;
		}
		else if(nChr > nChr2)
		{
			ni = nk;
		}
		else if(nEnd < nStart2)
		{
			nj = nk;
		}
		else if(nStart > nEnd2)
		{
			ni = nk;
		}
		else
		{
			nIdx = nk;
			nMaxN = (int)(DMGETAT(pRegion, ni, 4));

			ni = nk-1;
			while(ni >= 0)
			{
				nChr2 = (int)(DMGETAT(pRegion, ni, 0));
				nStart2 = (int)(DMGETAT(pRegion, ni, 1));
				nEnd2 = (int)(DMGETAT(pRegion, ni, 2));
				nN = (int)(DMGETAT(pRegion, ni, 4));

				if( (nChr == nChr2) && (nStart <= nEnd2) && ( nEnd>= nStart2) )
				{
					if(nN > nMaxN)
					{
						nMaxN = nN;
						nIdx = ni;
					}
				}
				else
				{
					break;
				}

				ni--;
			}

			ni = nk+1;
			while(ni < pRegion->nHeight)
			{
				nChr2 = (int)(DMGETAT(pRegion, ni, 0));
				nStart2 = (int)(DMGETAT(pRegion, ni, 1));
				nEnd2 = (int)(DMGETAT(pRegion, ni, 2));
				nN = (int)(DMGETAT(pRegion, ni, 4));

				if( (nChr == nChr2) && (nStart <= nEnd2) && ( nEnd>= nStart2) )
				{
					if(nN > nMaxN)
					{
						nMaxN = nN;
						nIdx = ni;
					}
				}
				else
				{
					break;
				}

				ni++;
			}

			break;
		}
	}

	if(nIdx < 0)
	{
		nChr2 = (int)(DMGETAT(pRegion, ni, 0));
		nStart2 = (int)(DMGETAT(pRegion, ni, 1));
		nEnd2 = (int)(DMGETAT(pRegion, ni, 2));

		if( (nChr == nChr2) && (nStart <= nEnd2) && ( nEnd>= nStart2) )
		{
			nIdx = ni;
		}
	}

	if(nIdx < 0)
	{
		nChr2 = (int)(DMGETAT(pRegion, nj, 0));
		nStart2 = (int)(DMGETAT(pRegion, nj, 1));
		nEnd2 = (int)(DMGETAT(pRegion, nj, 2));

		if( (nChr == nChr2) && (nStart <= nEnd2) && ( nEnd>= nStart2) )
		{
			nIdx = nj;
		}
	}

	/* return */
	return nIdx;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_OneSample_ExportResults()                                   */
/*  Export results.                                                        */
/* ----------------------------------------------------------------------- */ 
int	HTS_Enrich_OneSample_ExportResults(struct DOUBLEMATRIX *pRegion, 
				struct tagBARData *pBARData, char strFileName[])
{
	/* define */
	FILE *fpOut;
	int ni,nj,nk;

	/* open file */
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: HTS_Enrich_OneSample_ExportResults, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "#rank\tchr\tstart\tend\tstrand\tlength\tmaxN\tmaxN_pos\ttotal_reads\n");

	/* write */
	for(nj=0; nj<pRegion->nHeight; nj++)
	{
		ni = pRegion->nHeight-1-nj; 
		nk = (int)(DMGETAT(pRegion, ni,0));
		fprintf(fpOut, "%d\t%s\t%d\t%d\t+\t%d\t%d\t%d\t%d\n", nj+1, pBARData->vSeqData[nk]->pSeqName->m_pString,
			(int)(DMGETAT(pRegion, ni,1)), (int)(DMGETAT(pRegion, ni,2)), (int)(DMGETAT(pRegion, ni,3)),
			(int)(DMGETAT(pRegion, ni,4)), (int)(DMGETAT(pRegion, ni,5)), (int)(DMGETAT(pRegion, ni,6)));
	}

	/* close file */
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_OneSamplev2_ExportResults()                                 */
/*  Export results.                                                        */
/* ----------------------------------------------------------------------- */ 
int	HTS_Enrich_OneSamplev2_ExportResults_Old(struct DOUBLEMATRIX *pRegion, 
				struct tagBARData *pBARData, char strFileName[])
{
	/* define */
	FILE *fpOut;
	int ni,nj,nk;

	/* open file */
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: HTS_Enrich_OneSample_ExportResults, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "#rank\tchr\tstart\tend\tstrand\tlength\tmaxN\tmaxN_pos\ttotal_reads\tFstart\tFend\tFmaxN\tFmaxN_pos\tFtot_reads\tRstart\tRend\tRmaxN\tRmaxN_pos\tRtot_reads\n");

	/* write */
	for(nj=0; nj<pRegion->nHeight; nj++)
	{
		ni = pRegion->nHeight-1-nj; 
		nk = (int)(DMGETAT(pRegion, ni,0));
		fprintf(fpOut, "%d\t%s\t%d\t%d\t+\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", 
			nj+1, pBARData->vSeqData[nk]->pSeqName->m_pString,
			(int)(DMGETAT(pRegion, ni,1)), (int)(DMGETAT(pRegion, ni,2)), (int)(DMGETAT(pRegion, ni,3)),
			(int)(DMGETAT(pRegion, ni,4)), (int)(DMGETAT(pRegion, ni,5)), (int)(DMGETAT(pRegion, ni,6)),
			(int)(DMGETAT(pRegion, ni,7)), (int)(DMGETAT(pRegion, ni,8)), (int)(DMGETAT(pRegion, ni,9)),
			(int)(DMGETAT(pRegion, ni,10)), (int)(DMGETAT(pRegion, ni,11)), (int)(DMGETAT(pRegion, ni,12)),
			(int)(DMGETAT(pRegion, ni,13)), (int)(DMGETAT(pRegion, ni,14)), (int)(DMGETAT(pRegion, ni,15)),
			(int)(DMGETAT(pRegion, ni,16)));
	}

	/* close file */
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_OneSamplev2_ExportResults()                                 */
/*  Export results.                                                        */
/* ----------------------------------------------------------------------- */ 
int	HTS_Enrich_OneSamplev2_ExportResults(struct DOUBLEMATRIX *pRegion, 
				struct tagBARData *pBARData, char strFileName[], int nBR, int nBRL,
				int nSSF, int nSSFF, int nSSFR)
{
	/* define */
	FILE *fpOut;
	int ni,nj,nk,nz;
	double dMaxN,dMinN,dTemp;
	int nRegLen1,nRegLen2,nMR,nMF,nB1,nB2,nBRW,nP1,nP2;

	/* open file */
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: HTS_Enrich_OneSample_ExportResults, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "#rank\tchr\tstart\tend\tstrand\tlength\tmaxN\tmaxN_pos\ttotal_reads\tFtot_reads\tFmaxN\tFmaxN_pos\tRtot_reads\tRmaxN\tRmaxN_pos\tRmaxpos-Fmaxpos\tDelta\n");

	/* write */
	nz = 0;
	for(nj=0; nj<pRegion->nHeight; nj++)
	{
		ni = pRegion->nHeight-1-nj; 
		nk = (int)(DMGETAT(pRegion, ni, 0));
		nP1 = (int)(DMGETAT(pRegion, ni, 1));
		nP2 = (int)(DMGETAT(pRegion, ni, 2));
		nRegLen1 = (int)(DMGETAT(pRegion, ni,3));
		nB1 = (int)(DMGETAT(pRegion, ni,9));
		nB2 = (int)(DMGETAT(pRegion, ni,12));
		nRegLen2 = nB2-nB1+1;
		nMF = (int)(DMGETAT(pRegion, ni,8));
		nMR = (int)(DMGETAT(pRegion, ni,11));

		if(nSSF == 1)
		{
			if( (nMF < nSSFF) || (nMR < nSSFR) )
				continue;
		}
		if(nBR == 1)
		{
			if(nRegLen2 >= 0)
			{
				if(nRegLen2 < nBRL)
				{
					nBRW = (nBRL-nRegLen2)/2;
					nP1 = nB1-nBRW;
					nP2 = nB2+nBRW;
				}
				else
				{
					nP1 = nB1;
					nP2 = nB2;
				}
				nRegLen1 = nP2-nP1+1;
			}
		}

		fprintf(fpOut, "%d\t%s\t%d\t%d\t+\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t", 
			nz+1, pBARData->vSeqData[nk]->pSeqName->m_pString,
			nP1, nP2, nRegLen1,
			(int)(DMGETAT(pRegion, ni,4)), (int)(DMGETAT(pRegion, ni,5)), (int)(DMGETAT(pRegion, ni,6)),
			(int)(DMGETAT(pRegion, ni,7)), nMF, nB1,
			(int)(DMGETAT(pRegion, ni,10)), nMR, nB2,
			nRegLen2);

		dMaxN = DMGETAT(pRegion, ni,8);
		dMinN = DMGETAT(pRegion, ni,11);
		if(dMaxN < dMinN)
		{
			dTemp = dMinN;
			dMinN = dMaxN;
			dMaxN = dTemp;
		}
		dTemp = (dMinN+1.0)/(dMaxN+1.0);

		fprintf(fpOut, "%f\n", dTemp);
		nz++;
	}

	/* close file */
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_OneSample_RegionCollectInfo()                               */
/*  Collect reads number.                                                  */
/* ----------------------------------------------------------------------- */ 
int HTS_Enrich_OneSample_RegionCollectInfo(struct DOUBLEMATRIX *pRegion, struct tagBARData *pBARData)
{
	/* define */
	int ni,nj,nk,nl;
	int nChr,nStart,nEnd;
	int nP1,nP2;

	/* init */
	if( (pRegion == NULL) || (pBARData == NULL) )
		return PROC_SUCCESS;

	/* collect */
	for(ni=0; ni<pRegion->nHeight; ni++)
	{
		nChr = (int)(DMGETAT(pRegion, ni, 0));
		nStart = (int)(DMGETAT(pRegion, ni, 1));
		nEnd = (int)(DMGETAT(pRegion, ni, 2));

		if(pBARData->vSeqData[nChr]->nDataNum <= 0)
			continue;

		nj = 0;
		nk = pBARData->vSeqData[nChr]->nDataNum-1;

		if(nStart > pBARData->vSeqData[nChr]->vData[0]->pMatElement[nk])
		{
			nP1 = nk+1;
		}
		else if(nStart <= pBARData->vSeqData[nChr]->vData[0]->pMatElement[0])
		{
			nP1 = 0;
		}
		else
		{
			while( (nk-nj) > 1)
			{
				nl = (nk+nj)/2;
				if( pBARData->vSeqData[nChr]->vData[0]->pMatElement[nl] >= nStart)
				{
					nk = nl;
				}
				else
				{
					nj = nl;
				}
			}
			nP1 = nk;
		}

		nj = 0;
		nk = pBARData->vSeqData[nChr]->nDataNum-1;

		if(nEnd >= pBARData->vSeqData[nChr]->vData[0]->pMatElement[nk])
		{
			nP2 = nk;
		}
		else if(nEnd < pBARData->vSeqData[nChr]->vData[0]->pMatElement[0])
		{
			nP2 = -1;
		}
		else
		{
			while( (nk-nj) > 1)
			{
				nl = (nk+nj)/2;
				if( pBARData->vSeqData[nChr]->vData[0]->pMatElement[nl] > nEnd)
				{
					nk = nl;
				}
				else
				{
					nj = nl;
				}
			}
			nP2 = nj;
		}

		nk = 0;
		for(nj=nP1; nj<=nP2; nj++)
		{
			nk += (int)(pBARData->vSeqData[nChr]->vData[1]->pMatElement[nj]);
		}

		DMSETAT(pRegion, ni, 6, (double)nk);
	}


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_OneSamplev2_RegionCollectInfo()                             */
/*  Collect forward/reverse reads info                                     */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *HTS_Enrich_OneSamplev2_RegionCollectInfo(struct DOUBLEMATRIX *pRegion, 
				struct tagBARData *pBARData, int nW)
{
	/* define */
	struct DOUBLEMATRIX *pNewRegion = NULL;

	/* define */
	int ni,nj,nk,nl;
	int nChr,nStart,nEnd;
	int nP1,nP2,nQ1,nQ2;
	double dTemp;
	int nN = 0;
	int nMaxN = 0;
	int nMaxNPos,nMaxNPos2;
	int nW2 = nW/2;

	
	/* init */
	if( (pRegion == NULL) || (pBARData == NULL) )
		return NULL;

	pNewRegion = CreateDoubleMatrix(pRegion->nHeight, pRegion->nWidth+3);
	if(pNewRegion == NULL)
	{
		printf("Error: cannot create memory for collecting forward/reverse read information!\n");
		exit(EXIT_FAILURE);
	}

	/* collect */
	for(ni=0; ni<pRegion->nHeight; ni++)
	{
		for(nj=0; nj<pRegion->nWidth; nj++)
		{
			dTemp = DMGETAT(pRegion, ni, nj);
			DMSETAT(pNewRegion, ni, nj, dTemp);
		}

		nChr = (int)(DMGETAT(pRegion, ni, 0));
		nStart = (int)(DMGETAT(pRegion, ni, 1));
		nEnd = (int)(DMGETAT(pRegion, ni, 2));

		if(pBARData->vSeqData[nChr]->nDataNum <= 0)
			continue;

		nj = 0;
		nk = pBARData->vSeqData[nChr]->nDataNum-1;

		if(nStart > pBARData->vSeqData[nChr]->vData[0]->pMatElement[nk])
		{
			nP1 = nk+1;
		}
		else if(nStart <= pBARData->vSeqData[nChr]->vData[0]->pMatElement[0])
		{
			nP1 = 0;
		}
		else
		{
			while( (nk-nj) > 1)
			{
				nl = (nk+nj)/2;
				if( pBARData->vSeqData[nChr]->vData[0]->pMatElement[nl] >= nStart)
				{
					nk = nl;
				}
				else
				{
					nj = nl;
				}
			}
			nP1 = nk;
		}

		nj = 0;
		nk = pBARData->vSeqData[nChr]->nDataNum-1;

		if(nEnd >= pBARData->vSeqData[nChr]->vData[0]->pMatElement[nk])
		{
			nP2 = nk;
		}
		else if(nEnd < pBARData->vSeqData[nChr]->vData[0]->pMatElement[0])
		{
			nP2 = -1;
		}
		else
		{
			while( (nk-nj) > 1)
			{
				nl = (nk+nj)/2;
				if( pBARData->vSeqData[nChr]->vData[0]->pMatElement[nl] > nEnd)
				{
					nk = nl;
				}
				else
				{
					nj = nl;
				}
			}
			nP2 = nj;
		}

		nk = 0;
		for(nj=nP1; nj<=nP2; nj++)
		{
			nk += (int)(pBARData->vSeqData[nChr]->vData[1]->pMatElement[nj]);
		}

		DMSETAT(pNewRegion, ni, pRegion->nWidth, (double)nk);


		/* fine peak */
		if(nP1 > nP2)
			continue;

		nQ1 = nStart;
		nQ2 = nEnd-nW+1;

		nMaxN = -1;
		nMaxNPos = -1;
		nMaxNPos2 = -1;

		nN = 0;
		nk = nP1;
		nl = nP1;
		for(nj=nQ1; nj<=nQ2; nj++)
		{
			for(; nl<=nP2; nl++)
			{
				if((int)(pBARData->vSeqData[nChr]->vData[0]->pMatElement[nl]) >= nj)
					break;
				
				nN -= (int)(pBARData->vSeqData[nChr]->vData[1]->pMatElement[nl]);
			}

			for(; nk<=nP2; nk++)
			{
				if((int)(pBARData->vSeqData[nChr]->vData[0]->pMatElement[nk]) >= (nj+nW))
					break;

				nN += (int)(pBARData->vSeqData[nChr]->vData[1]->pMatElement[nk]);
			}

			if(nN > nMaxN)
			{
				nMaxN = nN;
				nMaxNPos = nj+nW2;
				nMaxNPos2 = nMaxNPos;
			}
			else if(nN == nMaxN)
			{
				nMaxNPos2 = nj+nW2;
			}
		}

		nMaxNPos = (nMaxNPos+nMaxNPos2)/2;
		DMSETAT(pNewRegion, ni, (pRegion->nWidth+1), (double)nMaxN);
		DMSETAT(pNewRegion, ni, (pRegion->nWidth+2), (double)nMaxNPos);
	}

	/* return */
	return pNewRegion;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_TwoSample_WindowSummary()                                          */
/*  Summarize window counts.                                               */
/* ----------------------------------------------------------------------- */ 
int HTS_TwoSample_WindowSummary(char strPosBARFile[], char strNegBARFile[],
					char strChrList[], char strChrLen[], 
					int nW, char strOutFile[])
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;

	struct tagBARData *pBARDataPos = NULL;
	struct tagBARData *pBARDataNeg = NULL;
	struct INTMATRIX *pChrLen = NULL;
	struct DOUBLEMATRIX *pCount = NULL;
	struct DOUBLEMATRIX *pCount2D = NULL;
	
	struct DOUBLEMATRIX *pFDR = NULL;
	double dPi0 = 0.0;
	
	int nMaxC = 0;
	int nMaxCPos = 0;
	int nMaxCNeg = 0;
	int nUpperC = 4096;
	int *vCPos;
	int *vCNeg;

	double dTotal;
	
	char strLine[LONG_LINE_LENGTH];
	int ni,nj,nk,nl,nLen;

	double dR0,dR1,dR2;
	double dlambda,dpoisp,dalpha,dbeta,dnegp;
	double dP1,dP2,dTemp;

	/* initial check */
	if(nW <= 0)
	{
		nW = 100;
		printf("Warning: Window size<=0! Proceed with default window size = 100\n");
	}

	/* load bar data */
	pBARDataPos = NULL;
	pBARDataPos = Affy_LoadBAR_Fast(strPosBARFile);
	if(pBARDataPos == NULL)
	{
		printf("Error: HTS_TwoSample_WindowSummary, cannot load raw data!\n");
		exit(EXIT_FAILURE);
	}

	pBARDataNeg = NULL;
	pBARDataNeg = Affy_LoadBAR_Fast(strNegBARFile);
	if(pBARDataNeg == NULL)
	{
		printf("Error: HTS_TwoSample_WindowSummary, cannot load raw data!\n");
		exit(EXIT_FAILURE);
	}

	/* load chromosome length */
	pChrLen = IMLOAD(strChrLen);
	if(pChrLen == NULL)
	{
		printf("Error: HTS_TwoSample_WindowSummary, cannot load chromosome length!\n");
		exit(EXIT_FAILURE);
	}

	/* create count matrix */
	pCount = CreateDoubleMatrix(nUpperC, 1);
	if(pCount == NULL)
	{
		printf("Error: HTS_TwoSample_WindowSummary, cannot create count matrix!\n");
		exit(EXIT_FAILURE);
	}

	pCount2D = CreateDoubleMatrix(nUpperC, nUpperC);
	if(pCount2D == NULL)
	{
		printf("Error: HTS_TwoSample_WindowSummary, cannot create count matrix!\n");
		exit(EXIT_FAILURE);
	}

	/* process chromosome one by one */
	fpIn = NULL;
	fpIn = fopen(strChrList, "r");
	if(fpIn == NULL)
	{
		printf("Error: HTS_TwoSample_WindowSummary, cannot open chrlist file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	nMaxCPos = 0;
	nMaxCNeg = 0;
	nMaxC = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		nLen = pChrLen->pMatElement[ni]/nW;
		if(pChrLen->pMatElement[ni]%nW != 0)
			nLen += 1;

		/* prepare memory */
		vCPos = NULL;
		vCPos = (int *)calloc(nLen, sizeof(int));
		if(vCPos == NULL)
		{
			printf("Error: HTS_TwoSample_WindowSummary, cannot create memory for counting!\n");
			exit(EXIT_FAILURE);
		}

		vCNeg = NULL;
		vCNeg = (int *)calloc(nLen, sizeof(int));
		if(vCNeg == NULL)
		{
			printf("Error: HTS_TwoSample_WindowSummary, cannot create memory for counting!\n");
			exit(EXIT_FAILURE);
		}

		/* find matching chromosome */
		for(nj=0; nj<pBARDataPos->nSeqNum; nj++)
		{
			if(strcmp(strLine, pBARDataPos->vSeqData[nj]->pSeqName->m_pString) != 0)
				continue;

			for(nk=0; nk<pBARDataPos->vSeqData[nj]->nDataNum; nk++)
			{
				nl = (int)(pBARDataPos->vSeqData[nj]->vData[0]->pMatElement[nk])/nW;
				if(nl == nLen)
					nl = nLen-1;
				else if(nl > nLen)
				{
					printf("Error: HTS_TwoSample_WindowSummary, index out of range! Please check if the genomes are matching\n");
					exit(EXIT_FAILURE);
				}

				vCPos[nl] += (int)(pBARDataPos->vSeqData[nj]->vData[1]->pMatElement[nk]);
			}
		}

		for(nj=0; nj<pBARDataNeg->nSeqNum; nj++)
		{
			if(strcmp(strLine, pBARDataNeg->vSeqData[nj]->pSeqName->m_pString) != 0)
				continue;

			for(nk=0; nk<pBARDataNeg->vSeqData[nj]->nDataNum; nk++)
			{
				nl = (int)(pBARDataNeg->vSeqData[nj]->vData[0]->pMatElement[nk])/nW;
				if(nl == nLen)
					nl = nLen-1;
				else if(nl > nLen)
				{
					printf("Error: HTS_TwoSample_WindowSummary, index out of range! Please check if the genomes are matching\n");
					exit(EXIT_FAILURE);
				}

				vCNeg[nl] += (int)(pBARDataNeg->vSeqData[nj]->vData[1]->pMatElement[nk]);
			}
		}

		/* count */
		for(nj=0; nj<nLen; nj++)
		{
			nk = vCPos[nj]+vCNeg[nj];
			if(nk >= nUpperC)
				nk = nUpperC-1;

			if(nk>nMaxC)
				nMaxC = nk;

			pCount->pMatElement[nk] += 1;

			nl = vCPos[nj];
			if(nl >= nUpperC)
				nl = nUpperC-1;

			if(vCPos[nj]>nMaxCPos)
				nMaxCPos = vCPos[nj];
			if(vCNeg[nj]>nMaxCNeg)
				nMaxCNeg = vCNeg[nj];


			dTemp = DMGETAT(pCount2D, nk, nl)+1.0;
			DMSETAT(pCount2D, nk, nl, dTemp);
		}

		/* free memory */
		free(vCPos);
		free(vCNeg);

		ni++;
	}
	
	fclose(fpIn);

	if(nMaxCPos >= nUpperC)
		nMaxCPos = nUpperC-1;
	if(nMaxCNeg >= nUpperC)
		nMaxCNeg = nUpperC-1;

	/* TODO: compute FDR */
	HTS_TwoSample_FDR(pCount, pCount2D, nMaxC, nMaxCPos, &pFDR, &dPi0);
	sprintf(strLine, "%s.fdr", strOutFile);
	DMSAVE(pFDR, strLine);
	DestroyDoubleMatrix(pFDR);

	/* output summary statistics */
	dTotal = 0.0;
	for(ni=0; ni<=nMaxC; ni++)
	{
		dTotal += pCount->pMatElement[ni];
	}

	/* estimate poisson and neg-binonmial */
	dR0 = pCount->pMatElement[0]/dTotal;
	dR1 = pCount->pMatElement[1]/dTotal;
	dR2 = pCount->pMatElement[2]/dTotal;
	dR2 = dR2/dR1;
	dR1 = dR1/dR0;

	dlambda = dR1;
	dpoisp = dR0/exp(-dlambda);
	if(dpoisp > 1.0)
		dpoisp = 1.0;

	dalpha = dR1/(2.0*dR2-dR1);
	dbeta = 1.0/(2.0*dR2-dR1)-1.0;
	dnegp = dR0/pow( (dbeta/(dbeta+1.0)), dalpha);
	if(dnegp > 1.0)
		dnegp = 1.0;


	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: HTS_WindowSummary, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "#Window_Size=%d\tdP0_hat=%f\n", nW, dPi0);
	fprintf(fpOut, "#Poisson_Lambda=%f\tPoisson_p=%f\n", dlambda, dpoisp);
	fprintf(fpOut, "#NegBinomial_Alpha=%f\tNegBinomial_Beta=%f\tNegBinomial_p=%f\n", dalpha, dbeta, dnegp);
	fprintf(fpOut, "#No_of_reads/window\tNo_of_window\tpercentage\tpoisson_expected\tpoisson_exp/obs\tnegbinomial_expected\tnegbinomial_exp/obs\n");
	for(ni=0; ni<=nMaxC; ni++)
	{
		dP1 = 0.0;
		dP1 = ni*log(dlambda)-dlambda-gammaln(ni+1.0)+log(dpoisp);
		dP1 = exp(dP1);

		dP2 = 0.0;
		dP2 = gammaln(ni+dalpha)-gammaln(ni+1.0)-gammaln(dalpha)+dalpha*log(dbeta/(dbeta+1.0))-ni*log(dbeta+1.0)+log(dnegp);
		dP2 = exp(dP2);

		dTemp = pCount->pMatElement[ni]/dTotal;
		
		fprintf(fpOut, "%d\t%d\t%f\t%f\t%f\t%f\t%f\n", ni, (int)(pCount->pMatElement[ni]), dTemp, dP1, dP1/(dTemp+1e-20), dP2, dP2/(dTemp+1e-20));
	}

	fprintf(fpOut, "\n\n");
	fprintf(fpOut, "#total_reads/window\tpos_reads/window\n");
	fprintf(fpOut, "#");
	for(ni=0; ni<=nMaxCPos; ni++)
	{
		fprintf(fpOut, "\t%d", ni);
	}
	fprintf(fpOut, "\n");

	for(ni=0; ni<=nMaxC; ni++)
	{
		fprintf(fpOut, "%d", ni);
		for(nj=0; nj<=nMaxCPos; nj++)
		{
			dTemp = DMGETAT(pCount2D, ni, nj);
			fprintf(fpOut, "\t%d", (int)dTemp);
		}
		fprintf(fpOut, "\n");
	}

	fclose(fpOut);

	/* destroy memory */
	Affy_BARData_Destroy(&pBARDataPos);
	Affy_BARData_Destroy(&pBARDataNeg);
	DestroyIntMatrix(pChrLen);
	DestroyDoubleMatrix(pCount);
	DestroyDoubleMatrix(pCount2D);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_TwoSample_WindowSummaryv2()                                        */
/*  Summarize window counts.                                               */
/* ----------------------------------------------------------------------- */ 
int HTS_TwoSample_WindowSummaryv2(char strPosBARFile[], char strNegBARFile[],
					char strChrList[], char strChrLen[], 
					int nW, char strOutFile[], int nCombineShift)
{
	/* define */
	char strFileNamePos[MED_LINE_LENGTH];
	char strFileNameNeg[MED_LINE_LENGTH];
	char strFileName1[MED_LINE_LENGTH];
	char strFileName2[MED_LINE_LENGTH];
	FILE *fpIn;
	FILE *fpOut;
	char strLine[LONG_LINE_LENGTH];

	/* run one by one */
	sprintf(strFileName1, "%s.pos", strOutFile);
	sprintf(strFileName2, "%s.neg", strOutFile);

	if(nCombineShift == 1)
	{
		sprintf(strFileNamePos, "%s_C.bar", strPosBARFile);
		sprintf(strFileNameNeg, "%s_C.bar", strNegBARFile);
		HTS_TwoSample_WindowSummary(strFileNamePos, strFileNameNeg,
					strChrList, strChrLen, 
					nW, strOutFile);
	}
	else
	{
		HTS_TwoSample_WindowSummary(strPosBARFile, strNegBARFile,
					strChrList, strChrLen, 
					nW, strOutFile);
	}

	sprintf(strFileNamePos, "%s_F.bar", strPosBARFile);
	sprintf(strFileNameNeg, "%s_F.bar", strNegBARFile);
	HTS_TwoSample_WindowSummary(strFileNamePos, strFileNameNeg,
					strChrList, strChrLen, 
					nW, strFileName1);

	sprintf(strFileNamePos, "%s_R.bar", strPosBARFile);
	sprintf(strFileNameNeg, "%s_R.bar", strNegBARFile);
	HTS_TwoSample_WindowSummary(strFileNamePos, strFileNameNeg,
					strChrList, strChrLen, 
					nW, strFileName2);

	/* summarize results */
	fpOut = NULL;
	fpOut = fopen(strOutFile, "at");
	if(fpOut == NULL)
	{
		printf("Error: HTS_TwoSample_WindowSummaryv2, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "\n#Forward_Only\n");
	fpIn = NULL;
	fpIn = fopen(strFileName1, "r");
	if(fpIn == NULL)
	{
		printf("Error: HTS_WindowSummaryv2, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);

		fprintf(fpOut, "%s\n", strLine);
	}

	fclose(fpIn);

	fprintf(fpOut, "\n#Reverse_Only\n");
	fpIn = NULL;
	fpIn = fopen(strFileName2, "r");
	if(fpIn == NULL)
	{
		printf("Error: HTS_WindowSummaryv2, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);

		fprintf(fpOut, "%s\n", strLine);
	}

	fclose(fpIn);

	fclose(fpOut);

	/* clear temp files */
	RemoveFiles(strFileName1);
	RemoveFiles(strFileName2);

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  HTS_TwoSample_FDR()                                                    */
/*  Compute FDR for two sample comparison.                                 */
/* ----------------------------------------------------------------------- */ 
int HTS_TwoSample_FDR(struct DOUBLEMATRIX *pCount, struct DOUBLEMATRIX *pCount2D, 
					  int nMaxC, int nMaxCPos, struct DOUBLEMATRIX **pFDR, 
					  double *dP0)
{
	/* define */
	int nR2I = 1;
	int nR2J;
	int ni,nj;
	double dR;
	double dRCut = 2.0;
	double dPi0 = 0.5;
	double dL0,dL1;
	double dT,dE,dO,dF;
	double dTemp1,dTemp2;
	double dMaxE;
	int nMaxj;

	/* init */
	if( (nMaxC >= pCount2D->nHeight) || (nMaxCPos >= pCount2D->nWidth) )
	{
		printf("Error: HTS_TwoSample_FDR, max count inconsistent!\n");
		exit(EXIT_FAILURE);
	}

	if( (nMaxC < 1) || (nMaxCPos < 1) )
	{
		printf("Error: HTS_TwoSample_FDR, too few data to estimate FDR!\n");
		return PROC_FAILURE;
	}

	/* get p0 */
	dPi0 = (DMGETAT(pCount2D, 1, 0)+1.0)/(DMGETAT(pCount2D, 1, 1)+1.0);
	dPi0 = 1.0/(dPi0+1.0);
	(*dP0) = dPi0;
	dL0 = log(1.0-dPi0);
	dL1 = log(dPi0);

	/* get r2 cutoff */
	for(ni=1; ni<=nMaxC; ni++)
	{
		dR = pCount->pMatElement[ni]/(double)ni;
		if(dR < dRCut)
		{
			break;
		}
	}

	nR2I = ni;
	nR2J = nR2I;
	if((nMaxCPos+1) < nR2J)
		nR2J = nMaxCPos+1;

	/* get proportion, expectation, and fdr */
	(*pFDR) = NULL;
	(*pFDR) = CreateDoubleMatrix(nR2I, nR2J);
	if( *pFDR == NULL )
	{
		printf("Error: HTS_TwoSample_FDR, cannot create FDR matrix!\n");
	}

	DMSETAT((*pFDR), 0, 0, 1.0);
	for(ni=1; ni<nR2I; ni++)
	{
		dT = log(pCount->pMatElement[ni]+ni+1);
		nMaxj = -1;

		/* compute */
		for(nj=0; nj<=ni; nj++)
		{
			if(nj >= nR2J)
				break;
			
			dE = dT+gammaln((double)(ni+1))-gammaln((double)(nj+1))-gammaln((double)(ni-nj+1))+(ni-nj)*dL0+nj*dL1;
			
			if(nMaxj == -1)
			{
				nMaxj = nj;
				dMaxE = dE;
			}
			else
			{
				if(dE > dMaxE)
				{
					dMaxE = dE;
					nMaxj = nj;
				}
			}

			dE = exp(dE);

			dO = DMGETAT(pCount2D, ni, nj)+1.0;
			dF = dE/dO;
			if( dF > 1.0 )
				dF = 1.0;

			DMSETAT( (*pFDR), ni, nj, dF);
		}

		/* adjust monotonicity */
		for(nj = nMaxj-1; nj>=0; nj--)
		{
			dTemp1 = DMGETAT((*pFDR), ni, nj);
			dTemp2 = DMGETAT((*pFDR), ni, (nj+1));
			if(dTemp1 > dTemp2)
				dTemp1 = dTemp2;
			DMSETAT( (*pFDR), ni, nj, dTemp1);
		}

		for(nj = nMaxj+1; nj<=ni; nj++)
		{
			if(nj >= nR2J)
				break;

			dTemp1 = DMGETAT((*pFDR), ni, nj);
			dTemp2 = DMGETAT((*pFDR), ni, (nj-1));
			if(dTemp1 > dTemp2)
				dTemp1 = dTemp2;
			DMSETAT( (*pFDR), ni, nj, dTemp1);
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_TwoSample_Main()                                            */
/*  Find enriched region from a sequencing data set.                       */
/* ----------------------------------------------------------------------- */ 
int HTS_Enrich_TwoSample_Main(char strPosBARFile[], char strNegBARFile[], 
					  int nOneSide, int nW, int nS, int nTCut, char strFDRFile[], double dFDRCut,
					  int nMinLen, int nMaxGap, double dP0,
					  char strExportFolder[], char strOutFileTitle[])
{
	/* define */
	struct tagBARData *pBARDataPos = NULL;
	struct tagBARData *pBARDataNeg = NULL;
	struct DOUBLEMATRIX *pFDR = NULL;

	struct DOUBLEMATRIX *pRegion0;
	struct DOUBLEMATRIX *pRegion;
	char strFileName[MED_LINE_LENGTH];
	struct INTMATRIX *pType = NULL;
	struct INTMATRIX *pPriority = NULL;
	struct DOUBLEMATRIX *pRegionSort = NULL;
	struct LONGMATRIX *pRegionSid = NULL;

	/* initial check */
	AdjustDirectoryPath(strExportFolder);
	
	/* load bar data */
	pBARDataPos = NULL;
	pBARDataPos = Affy_LoadBAR_Fast(strPosBARFile);
	if(pBARDataPos == NULL)
	{
		printf("Error: HTS_Enrich_TwoSample_Main, cannot load raw data!\n");
		exit(EXIT_FAILURE);
	}

	pBARDataNeg = NULL;
	pBARDataNeg = Affy_LoadBAR_Fast(strNegBARFile);
	if(pBARDataNeg == NULL)
	{
		printf("Error: HTS_Enrich_TwoSample_Main, cannot load raw data!\n");
		exit(EXIT_FAILURE);
	}

	pFDR = NULL;
	pFDR = DMLOAD(strFDRFile);
	if(pFDR == NULL)
	{
		printf("Error: HTS_Enrich_TwoSample_Main, cannot load FDR!\n");
		exit(EXIT_FAILURE);
	}


	/* initial region call */
	pRegion0 = NULL;
	pRegion0 = HTS_Enrich_TwoSample_CallRegion_Initial(pBARDataPos, pBARDataNeg,
					  nW, nS, nTCut, pFDR, dFDRCut, dP0, nOneSide,
					  strExportFolder, strOutFileTitle);

	/* merge and filter regions */
	pRegion = HTS_Enrich_TwoSample_MergeRegion(pRegion0, nMinLen, nMaxGap);
	DestroyDoubleMatrix(pRegion0);
	

	/* collect region information */
	HTS_Enrich_TwoSample_RegionCollectInfo(pRegion, pBARDataPos, pBARDataNeg);
	
	/* sort regions */
	pType = NULL;
	pType = CreateIntMatrix(1, pRegion->nWidth);
	if(pType == NULL)
	{
		printf("Error: HTS_Enrich_TwoSample_Main, cannot allocate memory for summarizing results!\n");
		exit(EXIT_FAILURE);
	}
	pType->pMatElement[0] = 2;
	pType->pMatElement[1] = 2;
	pType->pMatElement[2] = 2;
	pType->pMatElement[3] = 2;
	pType->pMatElement[4] = 2;
	pType->pMatElement[5] = 1;
	pType->pMatElement[6] = 2;
	pType->pMatElement[7] = 1;
	pType->pMatElement[8] = 2;
	pType->pMatElement[9] = 2;
	pType->pMatElement[10] = 2;

	pPriority = NULL;
	pPriority = CreateIntMatrix(1, 3);
	if(pPriority == NULL)
	{
		printf("Error: HTS_Enrich_OneSample_Main, cannot allocate memory for summarizing results!\n");
		exit(EXIT_FAILURE);
	}
	pPriority->pMatElement[0] = 5;
	pPriority->pMatElement[1] = 7;
	pPriority->pMatElement[2] = 4;

	pRegionSort = NULL;
	pRegionSid = NULL;
	DMSORTROWS(pRegion, pType, pPriority, &pRegionSort, &pRegionSid);
	
	DestroyIntMatrix(pType);
	DestroyIntMatrix(pPriority);
	DestroyDoubleMatrix(pRegion);
	
	/* export result */
	/* rank, chr, start, end, strand, length, biggest n, pos of biggest n, total reads */
	sprintf(strFileName, "%s%s.cod", strExportFolder, strOutFileTitle);
	HTS_Enrich_TwoSample_ExportResults(pRegionSort, pBARDataPos, pBARDataNeg, strFileName);
	DestroyDoubleMatrix(pRegionSort);
	DestroyLongMatrix(pRegionSid);


	/* destroy memory */
	Affy_BARData_Destroy(&pBARDataPos);
	Affy_BARData_Destroy(&pBARDataNeg);
	DestroyDoubleMatrix(pFDR);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_TwoSamplev2_Main()                                          */
/*  Find enriched region from a sequencing data set.                       */
/* ----------------------------------------------------------------------- */ 
int HTS_Enrich_TwoSamplev2_Main(char strPosBARFile[], char strNegBARFile[], 
					  int nOneSide, int nW, int nS, int nTCut, char strFDRFile[], double dFDRCut,
					  int nMinLen, int nMaxGap, double dP0,
					  char strExportFolder[], char strOutFileTitle[],
					  int nBR, int nBRL, int nSSF, int nSSFF, int nSSFR, int nCombineShift,
					  double dFC, double dTFC)
{
	/* define */
	struct tagBARData *pBARDataPos = NULL;
	struct tagBARData *pBARDataNeg = NULL;
	struct DOUBLEMATRIX *pFDR = NULL;
	struct DOUBLEMATRIX *pFDRF = NULL;
	struct DOUBLEMATRIX *pFDRR = NULL;

	struct tagBARData *pBARDataPosF = NULL;
	struct tagBARData *pBARDataNegF = NULL;
	struct tagBARData *pBARDataPosR = NULL;
	struct tagBARData *pBARDataNegR = NULL;
	

	struct DOUBLEMATRIX *pRegion0;
	struct DOUBLEMATRIX *pRegion;
	struct DOUBLEMATRIX *pNewRegion;
	char strFileName[MED_LINE_LENGTH];
	struct INTMATRIX *pType = NULL;
	struct INTMATRIX *pPriority = NULL;
	struct DOUBLEMATRIX *pRegionSort = NULL;
	struct LONGMATRIX *pRegionSid = NULL;
	int ni;
	char *chp;

	/* initial check */
	AdjustDirectoryPath(strExportFolder);
	
	/* load bar data */
	pBARDataPos = NULL;
	if(nCombineShift == 1)
	{
		sprintf(strFileName, "%s_C.bar", strPosBARFile);
		pBARDataPos = Affy_LoadBAR_Fast(strFileName);
	}
	else
	{
		pBARDataPos = Affy_LoadBAR_Fast(strPosBARFile);
	}
	if(pBARDataPos == NULL)
	{
		printf("Error: HTS_Enrich_TwoSample_Main, cannot load raw data!\n");
		exit(EXIT_FAILURE);
	}

	pBARDataNeg = NULL;
	if(nCombineShift == 1)
	{
		sprintf(strFileName, "%s_C.bar", strNegBARFile);
		pBARDataNeg = Affy_LoadBAR_Fast(strFileName);
	}
	else
	{
		pBARDataNeg = Affy_LoadBAR_Fast(strNegBARFile);
	}
	if(pBARDataNeg == NULL)
	{
		printf("Error: HTS_Enrich_TwoSample_Main, cannot load raw data!\n");
		exit(EXIT_FAILURE);
	}

	pFDR = NULL;
	pFDR = DMLOAD(strFDRFile);
	if(pFDR == NULL)
	{
		printf("Error: HTS_Enrich_TwoSample_Main, cannot load FDR!\n");
		exit(EXIT_FAILURE);
	}

	chp = strstr(strFDRFile, ".fdr");
	if(chp == NULL)
	{
		printf("Error: HTS_Enrich_TwoSample_Main, FDR file needs to be ended with .fdr!\n");
		exit(EXIT_FAILURE);
	}
	*chp = '\0';
	sprintf(strFileName, "%s.pos.fdr", strFDRFile);
	pFDRF = NULL;
	pFDRF = DMLOAD(strFileName);
	if(pFDRF == NULL)
	{
		printf("Error: HTS_Enrich_TwoSample_Main, cannot load FDR!\n");
		exit(EXIT_FAILURE);
	}
	sprintf(strFileName, "%s.neg.fdr", strFDRFile);
	pFDRR = NULL;
	pFDRR = DMLOAD(strFileName);
	if(pFDRR == NULL)
	{
		printf("Error: HTS_Enrich_TwoSample_Main, cannot load FDR!\n");
		exit(EXIT_FAILURE);
	}


	/* initial region call */
	pRegion0 = NULL;
	pRegion0 = HTS_Enrich_TwoSample_CallRegion_Initial(pBARDataPos, pBARDataNeg,
					  nW, nS, nTCut, pFDR, dFDRCut, dP0, nOneSide,
					  strExportFolder, strOutFileTitle);

	/* merge and filter regions */
	pRegion = HTS_Enrich_TwoSample_MergeRegion(pRegion0, nMinLen, nMaxGap);
	DestroyDoubleMatrix(pRegion0);
	

	/* collect region information */
	HTS_Enrich_TwoSample_RegionCollectInfo(pRegion, pBARDataPos, pBARDataNeg);

	/* process forward */
	/* load bar data */
	sprintf(strFileName, "%s_F.bar", strPosBARFile);
	pBARDataPosF = NULL;
	pBARDataPosF = Affy_LoadBAR_Fast(strFileName);
	if(pBARDataPosF == NULL)
	{
		printf("Error: HTS_Enrich_TwoSample_Main, cannot load raw data!\n");
		exit(EXIT_FAILURE);
	}

	sprintf(strFileName, "%s_F.bar", strNegBARFile);
	pBARDataNegF = NULL;
	pBARDataNegF = Affy_LoadBAR_Fast(strFileName);
	if(pBARDataNegF == NULL)
	{
		printf("Error: HTS_Enrich_TwoSample_Main, cannot load raw data!\n");
		exit(EXIT_FAILURE);
	}

	/* initial region call */
	sprintf(strFileName, "%s_F", strOutFileTitle);
	pRegion0 = NULL;
	pRegion0 = HTS_Enrich_TwoSample_CallRegion_Initial(pBARDataPosF, pBARDataNegF,
					  nW, nS, nTCut, pFDR, dFDRCut, dP0, nOneSide,
					  strExportFolder, strFileName);
	DestroyDoubleMatrix(pRegion0);

	/* process reverse */
	/* load bar data */
	sprintf(strFileName, "%s_R.bar", strPosBARFile);
	pBARDataPosR = NULL;
	pBARDataPosR = Affy_LoadBAR_Fast(strFileName);
	if(pBARDataPosR == NULL)
	{
		printf("Error: HTS_Enrich_TwoSample_Main, cannot load raw data!\n");
		exit(EXIT_FAILURE);
	}

	sprintf(strFileName, "%s_R.bar", strNegBARFile);
	pBARDataNegR = NULL;
	pBARDataNegR = Affy_LoadBAR_Fast(strFileName);
	if(pBARDataNegR == NULL)
	{
		printf("Error: HTS_Enrich_TwoSample_Main, cannot load raw data!\n");
		exit(EXIT_FAILURE);
	}

	/* initial region call */
	sprintf(strFileName, "%s_R", strOutFileTitle);
	pRegion0 = NULL;
	pRegion0 = HTS_Enrich_TwoSample_CallRegion_Initial(pBARDataPosR, pBARDataNegR,
					  nW, nS, nTCut, pFDR, dFDRCut, dP0, nOneSide,
					  strExportFolder, strFileName);
	DestroyDoubleMatrix(pRegion0);
	

	pNewRegion = pRegion;
	pRegion = NULL;
	pRegion = HTS_Enrich_TwoSamplev2_RegionCollectInfo(pNewRegion, pBARDataPosF, pBARDataNegF, nW, pFDRF, dP0, nOneSide);
	DestroyDoubleMatrix(pNewRegion);
	pNewRegion = NULL;
	pNewRegion = pRegion;
	pRegion = HTS_Enrich_TwoSamplev2_RegionCollectInfo(pNewRegion, pBARDataPosR, pBARDataNegR, nW, pFDRR, dP0, nOneSide);
	DestroyDoubleMatrix(pNewRegion);

	/* sort regions */
	pType = NULL;
	pType = CreateIntMatrix(1, pRegion->nWidth);
	if(pType == NULL)
	{
		printf("Error: HTS_Enrich_TwoSample_Main, cannot allocate memory for summarizing results!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<pType->nWidth; ni++)
		pType->pMatElement[ni] = 2;
	pType->pMatElement[5] = 1;
	pType->pMatElement[7] = 1;
	
	pPriority = NULL;
	pPriority = CreateIntMatrix(1, 3);
	if(pPriority == NULL)
	{
		printf("Error: HTS_Enrich_TwoSample_Main, cannot allocate memory for summarizing results!\n");
		exit(EXIT_FAILURE);
	}
	pPriority->pMatElement[0] = 5;
	pPriority->pMatElement[1] = 7;
	pPriority->pMatElement[2] = 4;

	pRegionSort = NULL;
	pRegionSid = NULL;
	DMSORTROWS(pRegion, pType, pPriority, &pRegionSort, &pRegionSid);
	
	DestroyIntMatrix(pType);
	DestroyIntMatrix(pPriority);
	DestroyDoubleMatrix(pRegion);
	
	/* export result */
	/* rank, chr, start, end, strand, length, biggest n, pos of biggest n, total reads */
	sprintf(strFileName, "%s%s.cod", strExportFolder, strOutFileTitle);
	HTS_Enrich_TwoSamplev2_ExportResults(pRegionSort, pBARDataPos, pBARDataNeg, strFileName, 
		nBR, nBRL, nSSF, nSSFF, nSSFR, dFC, dTFC);
	DestroyDoubleMatrix(pRegionSort);
	DestroyLongMatrix(pRegionSid);


	/* destroy memory */
	Affy_BARData_Destroy(&pBARDataPos);
	Affy_BARData_Destroy(&pBARDataNeg);
	Affy_BARData_Destroy(&pBARDataPosF);
	Affy_BARData_Destroy(&pBARDataNegF);
	Affy_BARData_Destroy(&pBARDataPosR);
	Affy_BARData_Destroy(&pBARDataNegR);
	DestroyDoubleMatrix(pFDR);
	DestroyDoubleMatrix(pFDRF);
	DestroyDoubleMatrix(pFDRR);

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_TwoSample_CallRegion_Initial()                              */
/*  Search for differentially expressed windows.                           */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *HTS_Enrich_TwoSample_CallRegion_Initial(struct tagBARData *pBARDataPos,
			struct tagBARData *pBARDataNeg, int nW, int nS, int nTCut,
			struct DOUBLEMATRIX *pFDR, double dFDRCut, double dP0, int nOneSide,
			char strExportFolder[], char strOutFileTitle[])
{
	/* define */
	FILE *fpOut;
	FILE *fpReg;

	char strWinBarTmpFile[MED_LINE_LENGTH];
	int nCmpResult;
	double dR0 = dP0/(1.0-dP0);
	double dLog2 = log(2.0);
	
	char strRegTmpFile[MED_LINE_LENGTH];
	char strFileName[MED_LINE_LENGTH];
	struct DOUBLEMATRIX *pRegion = NULL;
	int ni1,ni2,nj,nk,nl,nk2,nl2,nx;
	int nP1,nP2;
	int nN1,nN2,nN;
	int nN1Sub,nNSub;
	double dFC,dTemp;
	int nW2 = nW/2;
	int nStart,nEnd;
	double dMinFDR;
	int nMinFDRPos1;
	int nMinFDRPos2;
	double dMaxFC;
	int nMaxFCPos1;
	int nMaxFCPos2;
	double dMinFC;
	int nMinFCPos1;
	int nMinFCPos2;
	int nPosSide;
	int nMaxN1,nMaxN2;
	int nMaxN1Pos1,nMaxN1Pos2;
	int nMaxN2Pos1,nMaxN2Pos2;



	/* init */
	if((pBARDataPos == NULL) || (pBARDataNeg == NULL))
		return NULL;

	/* process one by one */
	sprintf(strWinBarTmpFile, "%s%s.bar.tmp", strExportFolder, strOutFileTitle);
	
	fpOut = NULL;
	fpOut = fopen(strWinBarTmpFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: HTS_Enrich_TwoSample_CallRegion_Initial, cannot open temporary file to write display information!\n");
		exit(EXIT_FAILURE);
	}
	fprintf(fpOut, "#chr\tpos\t%s.pos\t%s.neg\t%s.log2fc\n", strOutFileTitle, strOutFileTitle, strOutFileTitle);
	fprintf(fpOut, "#chr\tpos\t1\t1\t1\n");

	sprintf(strRegTmpFile, "%s%s.regtmp", strExportFolder, strOutFileTitle);
	fpReg = NULL;
	fpReg = fopen(strRegTmpFile, "w");
	if(fpReg == NULL)
	{
		printf("Error: HTS_Enrich_TwoSample_CallRegion_Initial, cannot open temporary file to write region information!\n");
		exit(EXIT_FAILURE);
	}

	ni1 = 0;
	ni2 = 0;
	while( (ni1<pBARDataPos->nSeqNum) || (ni2<pBARDataNeg->nSeqNum) )
	{
		if( (ni1<pBARDataPos->nSeqNum) && (ni2<pBARDataNeg->nSeqNum) )
			nCmpResult = strcmp(pBARDataPos->vSeqData[ni1]->pSeqName->m_pString, pBARDataNeg->vSeqData[ni2]->pSeqName->m_pString);
		else if( ni1<pBARDataPos->nSeqNum )
			nCmpResult = -1;
		else
			nCmpResult = 1;

		/* process both bar_pos and bar_neg only */
		if(nCmpResult == 0)
		{
			if( (pBARDataPos->vSeqData[ni1]->nDataNum <= 0) && (pBARDataNeg->vSeqData[ni2]->nDataNum <= 0) )
			{
				ni1++;
				ni2++;
				continue;
			}

			if( pBARDataPos->vSeqData[ni1]->nDataNum > 0)
			{
				nx = pBARDataPos->vSeqData[ni1]->nDataNum-1;
				nP1 = (int)(pBARDataPos->vSeqData[ni1]->vData[0]->pMatElement[0])-nW;
				nP2 = (int)(pBARDataPos->vSeqData[ni1]->vData[0]->pMatElement[nx]);
			
				if( pBARDataNeg->vSeqData[ni2]->nDataNum > 0)
				{
					nx = pBARDataNeg->vSeqData[ni2]->nDataNum-1;
					if( (int)(pBARDataNeg->vSeqData[ni2]->vData[0]->pMatElement[0])-nW < nP1)
						nP1 = (int)(pBARDataNeg->vSeqData[ni2]->vData[0]->pMatElement[0])-nW;
					if( (int)(pBARDataNeg->vSeqData[ni2]->vData[0]->pMatElement[nx]) > nP2)
						nP2 = (int)(pBARDataNeg->vSeqData[ni2]->vData[0]->pMatElement[nx]);
				}
			}
			else
			{
				nx = pBARDataNeg->vSeqData[ni2]->nDataNum-1;
				nP1 = (int)(pBARDataNeg->vSeqData[ni2]->vData[0]->pMatElement[0])-nW;
				nP2 = (int)(pBARDataNeg->vSeqData[ni2]->vData[0]->pMatElement[nx]);
			}
			
			nStart = -1;
			nEnd = -1;

			dMinFDR = 1.0;
			nMinFDRPos1 = -1;
			nMinFDRPos2 = -1;

			dMaxFC = 0.0;
			nMaxFCPos1 = -1;
			nMaxFCPos2 = -1;

			dMinFC = 0.0;
			nMinFCPos1 = -1;
			nMinFCPos2 = -1;

			nMaxN1 = 0;
			nMaxN1Pos1 = -1;
			nMaxN1Pos2 = -1;

			nMaxN2 = 0;
			nMaxN2Pos1 = -1;
			nMaxN2Pos2 = -1;

			nN1 = 0;
			nN2 = 0;
			dFC = 0.0;
			nk = 0;
			nl = 0;
			nk2 = 0;
			nl2 = 0;
			for(nj=nP1; nj<=nP2; nj++)
			{
				for(; nl<pBARDataPos->vSeqData[ni1]->nDataNum; nl++)
				{
					if((int)(pBARDataPos->vSeqData[ni1]->vData[0]->pMatElement[nl]) >= nj)
						break;
					
					nN1 -= (int)(pBARDataPos->vSeqData[ni1]->vData[1]->pMatElement[nl]);
				}

				for(; nk<pBARDataPos->vSeqData[ni1]->nDataNum; nk++)
				{
					if((int)(pBARDataPos->vSeqData[ni1]->vData[0]->pMatElement[nk]) >= (nj+nW))
						break;

					nN1 += (int)(pBARDataPos->vSeqData[ni1]->vData[1]->pMatElement[nk]);
				}

				for(; nl2<pBARDataNeg->vSeqData[ni2]->nDataNum; nl2++)
				{
					if((int)(pBARDataNeg->vSeqData[ni2]->vData[0]->pMatElement[nl2]) >= nj)
						break;
					
					nN2 -= (int)(pBARDataNeg->vSeqData[ni2]->vData[1]->pMatElement[nl2]);
				}

				for(; nk2<pBARDataNeg->vSeqData[ni2]->nDataNum; nk2++)
				{
					if((int)(pBARDataNeg->vSeqData[ni2]->vData[0]->pMatElement[nk2]) >= (nj+nW))
						break;

					nN2 += (int)(pBARDataNeg->vSeqData[ni2]->vData[1]->pMatElement[nk2]);
				}

				dTemp = nN2*dR0;
				if(dTemp < (double)nN1)
				{
					dFC = log((double)(nN1+1)/(dTemp+1.0))/dLog2;
				}
				else
				{
					dFC = -log((dTemp+1.0)/(double)(nN1+1))/dLog2;
				}

				nN = nN1+nN2;
				
				if( (nj+nW2)%nS == 0)
				{
					if( ((nj+nW2)>=0) && (nN>0) )
						fprintf(fpOut, "%s\t%d\t%d\t%d\t%f\n", pBARDataPos->vSeqData[ni1]->pSeqName->m_pString, (int)(nj+nW2), nN1, nN2, dFC);
				}

				/* get fdr */
				nN1Sub = nN1;
				if(nN < pFDR->nHeight)
				{
					if(nN1 >= pFDR->nWidth)
						nN1Sub = pFDR->nWidth-1;

					dTemp = DMGETAT(pFDR, nN, nN1Sub);
				}
				else
				{
					nNSub = pFDR->nHeight-1;
					dTemp = (double)(nNSub*nN1)/(double)nN;
					nN1Sub = (int)dTemp;
					if(dTemp-nN1Sub > 0.5)
						nN1Sub += 1;
					
					if(nN1Sub >= pFDR->nWidth)
						nN1Sub = pFDR->nWidth-1;

					dTemp = DMGETAT(pFDR, nNSub, nN1Sub);
				}

				if( nOneSide == 1 )
				{
					if( dFC < 0.0 )
						dTemp = 1.0;
				}

				if( (nN >= nTCut) && (dTemp <= dFDRCut) )
				{
					if(nStart < 0)
					{
						nStart = nj;
						nEnd = nj+nW-1;

						dMinFDR = dTemp;
						nMinFDRPos1 = nj+nW2;
						nMinFDRPos2 = nMinFDRPos1;
						if(dFC >= 0.0)
							nPosSide = 1;
						else
							nPosSide = 0;

						dMaxFC = dFC;
						nMaxFCPos1 = nj+nW2;
						nMaxFCPos2 = nMaxFCPos1;

						dMinFC = dFC;
						nMinFCPos1 = nj+nW2;
						nMinFCPos2 = nMinFCPos1;

						nMaxN1 = nN1;
						nMaxN1Pos1 = nj+nW2;
						nMaxN1Pos2 = nMaxN1Pos1;

						
						nMaxN2 = nN2;
						nMaxN2Pos1 = nj+nW2;
						nMaxN2Pos2 = nMaxN2Pos1;
					}
					else
					{
						nEnd = nj+nW-1;

						if( fabs(dTemp-dMinFDR) < 1e-10 )
						{
							nMinFDRPos2 = nj+nW2;
						}
						else if(dTemp < dMinFDR)
						{
							dMinFDR = dTemp;
							nMinFDRPos1 = nj+nW2;
							nMinFDRPos2 = nMinFDRPos1;

							if(dFC >= 0.0)
								nPosSide = 1;
							else
								nPosSide = 0;
						}

						if( fabs(dFC-dMaxFC) < 1e-10 )
						{
							nMaxFCPos2 = nj+nW2;
						}
						else if(dFC > dMaxFC)
						{
							dMaxFC = dFC;
							nMaxFCPos1 = nj+nW2;
							nMaxFCPos2 = nMaxFCPos1;
						}

						if( fabs(dFC-dMinFC) < 1e-10 )
						{
							nMinFCPos2 = nj+nW2;
						}
						else if(dFC < dMinFC)
						{
							dMinFC = dFC;
							nMinFCPos1 = nj+nW2;
							nMinFCPos2 = nMinFCPos1;
						}

						if(nN1 == nMaxN1)
						{
							nMaxN1Pos2 = nj+nW2;
						}
						else if(nN1 > nMaxN1)
						{
							nMaxN1 = nN1;
							nMaxN1Pos1 = nj+nW2;
							nMaxN1Pos2 = nMaxN1Pos1;
						}

						if(nN2 == nMaxN2)
						{
							nMaxN2Pos2 = nj+nW2;
						}
						else if(nN2 > nMaxN2)
						{
							nMaxN2 = nN2;
							nMaxN2Pos1 = nj+nW2;
							nMaxN2Pos2 = nMaxN2Pos1;
						}
					}
				}
				else
				{
					if(nStart >= 0)
					{
						if(nPosSide == 1)
						{
							fprintf(fpReg, "0\t%d\t%d\t%d\t%f\t%d\t%f\t%d\t%d\t%d\t%d\t%d\n", ni1, nStart, nEnd, 
								dMinFDR, (nMinFDRPos1+nMinFDRPos2)/2,
								dMaxFC, (nMaxFCPos1+nMaxFCPos2)/2,
								nMaxN1, (nMaxN1Pos1+nMaxN1Pos2)/2,
								nMaxN2, (nMaxN2Pos1+nMaxN2Pos2)/2);
						}
						else
						{
							fprintf(fpReg, "0\t%d\t%d\t%d\t%f\t%d\t%f\t%d\t%d\t%d\t%d\t%d\n", ni1, nStart, nEnd, 
								dMinFDR, (nMinFDRPos1+nMinFDRPos2)/2,
								dMinFC, (nMinFCPos1+nMinFCPos2)/2,
								nMaxN1, (nMaxN1Pos1+nMaxN1Pos2)/2,
								nMaxN2, (nMaxN2Pos1+nMaxN2Pos2)/2);
						}
						nStart = -1;
						nEnd = -1;
						
						dMinFDR = 1.0;
						nMinFDRPos1 = -1;
						nMinFDRPos2 = -1;

						dMaxFC = 0.0;
						nMaxFCPos1 = -1;
						nMaxFCPos2 = -1;

						dMinFC = 0.0;
						nMinFCPos1 = -1;
						nMinFCPos2 = -1;

						nMaxN1 = 0;
						nMaxN1Pos1 = -1;
						nMaxN1Pos2 = -1;
						nMaxN2 = 0;
						nMaxN2Pos1 = -1;
						nMaxN2Pos2 = -1;
					}
				}
			}

			if(nStart >= 0)
			{
				if(nPosSide == 1)
				{
					fprintf(fpReg, "0\t%d\t%d\t%d\t%f\t%d\t%f\t%d\t%d\t%d\t%d\t%d\n", ni1, nStart, nEnd, 
						dMinFDR, (nMinFDRPos1+nMinFDRPos2)/2,
						dMaxFC, (nMaxFCPos1+nMaxFCPos2)/2,
						nMaxN1, (nMaxN1Pos1+nMaxN1Pos2)/2,
						nMaxN2, (nMaxN2Pos1+nMaxN2Pos2)/2);
				}
				else
				{
					fprintf(fpReg, "0\t%d\t%d\t%d\t%f\t%d\t%f\t%d\t%d\t%d\t%d\t%d\n", ni1, nStart, nEnd, 
						dMinFDR, (nMinFDRPos1+nMinFDRPos2)/2,
						dMinFC, (nMinFCPos1+nMinFCPos2)/2,
						nMaxN1, (nMaxN1Pos1+nMaxN1Pos2)/2,
						nMaxN2, (nMaxN2Pos1+nMaxN2Pos2)/2);
				}
				nStart = -1;
				nEnd = -1;
				dMinFDR = 1.0;
				nMinFDRPos1 = -1;
				nMinFDRPos2 = -1;

				dMaxFC = 0.0;
				nMaxFCPos1 = -1;
				nMaxFCPos2 = -1;

				dMinFC = 0.0;
				nMinFCPos1 = -1;
				nMinFCPos2 = -1;

				nMaxN1 = 0;
				nMaxN1Pos1 = -1;
				nMaxN1Pos2 = -1;
				nMaxN2 = 0;
				nMaxN2Pos1 = -1;
				nMaxN2Pos2 = -1;
			}

			ni1++;
			ni2++;
		}

		/* process bar_pos only */
		else if(nCmpResult < 0)
		{
			if(pBARDataPos->vSeqData[ni1]->nDataNum <= 0)
			{
				ni1++;
				continue;
			}

			nx = pBARDataPos->vSeqData[ni1]->nDataNum-1;
			nP1 = (int)(pBARDataPos->vSeqData[ni1]->vData[0]->pMatElement[0])-nW;
			nP2 = (int)(pBARDataPos->vSeqData[ni1]->vData[0]->pMatElement[nx]);
			nStart = -1;
			nEnd = -1;

			dMinFDR = 1.0;
			nMinFDRPos1 = -1;
			nMinFDRPos2 = -1;

			dMaxFC = 0.0;
			nMaxFCPos1 = -1;
			nMaxFCPos2 = -1;

			dMinFC = 0.0;
			nMinFCPos1 = -1;
			nMinFCPos2 = -1;

			nMaxN1 = 0;
			nMaxN1Pos1 = -1;
			nMaxN1Pos2 = -1;
			nMaxN2 = 0;
			nMaxN2Pos1 = -1;
			nMaxN2Pos2 = -1;


			nN1 = 0;
			nN2 = 0;
			dFC = 0.0;
			nk = 0;
			nl = 0;
			for(nj=nP1; nj<=nP2; nj++)
			{
				for(; nl<pBARDataPos->vSeqData[ni1]->nDataNum; nl++)
				{
					if((int)(pBARDataPos->vSeqData[ni1]->vData[0]->pMatElement[nl]) >= nj)
						break;
					
					nN1 -= (int)(pBARDataPos->vSeqData[ni1]->vData[1]->pMatElement[nl]);
				}

				for(; nk<pBARDataPos->vSeqData[ni1]->nDataNum; nk++)
				{
					if((int)(pBARDataPos->vSeqData[ni1]->vData[0]->pMatElement[nk]) >= (nj+nW))
						break;

					nN1 += (int)(pBARDataPos->vSeqData[ni1]->vData[1]->pMatElement[nk]);
				}

				dTemp = nN2*dR0;
				if(dTemp < (double)nN1)
				{
					dFC = log((double)(nN1+1)/(dTemp+1.0))/dLog2;
				}
				else
				{
					dFC = -log((dTemp+1.0)/(double)(nN1+1))/dLog2;
				}

				
				nN = nN1+nN2;

				if( (nj+nW2)%nS == 0)
				{
					if( ((nj+nW2)>=0) && (nN>0) )
						fprintf(fpOut, "%s\t%d\t%d\t%d\t%f\n", pBARDataPos->vSeqData[ni1]->pSeqName->m_pString, (int)(nj+nW2), nN1, nN2, dFC);
				}

				/* get fdr */
				nN1Sub = nN1;
				if(nN < pFDR->nHeight)
				{
					if(nN1 >= pFDR->nWidth)
						nN1Sub = pFDR->nWidth-1;

					dTemp = DMGETAT(pFDR, nN, nN1Sub);
				}
				else
				{
					nNSub = pFDR->nHeight-1;
					dTemp = (double)(nNSub*nN1)/(double)nN;
					nN1Sub = (int)dTemp;
					if(dTemp-nN1Sub > 0.5)
						nN1Sub += 1;
					
					if(nN1Sub >= pFDR->nWidth)
						nN1Sub = pFDR->nWidth-1;

					dTemp = DMGETAT(pFDR, nNSub, nN1Sub);
				}

				if( nOneSide == 1 )
				{
					if( dFC < 0.0 )
						dTemp = 1.0;
				}

				if( (nN >= nTCut) && (dTemp <= dFDRCut) )
				{
					if(nStart < 0)
					{
						nStart = nj;
						nEnd = nj+nW-1;

						dMinFDR = dTemp;
						nMinFDRPos1 = nj+nW2;
						nMinFDRPos2 = nMinFDRPos1;
						if(dFC >= 0.0)
							nPosSide = 1;
						else
							nPosSide = 0;

						dMaxFC = dFC;
						nMaxFCPos1 = nj+nW2;
						nMaxFCPos2 = nMaxFCPos1;

						dMinFC = dFC;
						nMinFCPos1 = nj+nW2;
						nMinFCPos2 = nMinFCPos1;

						nMaxN1 = nN1;
						nMaxN1Pos1 = nj+nW2;
						nMaxN1Pos2 = nMaxN1Pos1;

						nMaxN2 = nN2;
						nMaxN2Pos1 = nj+nW2;
						nMaxN2Pos2 = nMaxN2Pos1;
					}
					else
					{
						nEnd = nj+nW-1;

						if( fabs(dTemp-dMinFDR) < 1e-10 )
						{
							nMinFDRPos2 = nj+nW2;
						}
						else if(dTemp < dMinFDR)
						{
							dMinFDR = dTemp;
							nMinFDRPos1 = nj+nW2;
							nMinFDRPos2 = nMinFDRPos1;

							if(dFC >= 0.0)
								nPosSide = 1;
							else
								nPosSide = 0;
						}

						if( fabs(dFC-dMaxFC) < 1e-10 )
						{
							nMaxFCPos2 = nj+nW2;
						}
						else if(dFC > dMaxFC)
						{
							dMaxFC = dFC;
							nMaxFCPos1 = nj+nW2;
							nMaxFCPos2 = nMaxFCPos1;
						}

						if( fabs(dFC-dMinFC) < 1e-10 )
						{
							nMinFCPos2 = nj+nW2;
						}
						else if(dFC < dMinFC)
						{
							dMinFC = dFC;
							nMinFCPos1 = nj+nW2;
							nMinFCPos2 = nMinFCPos1;
						}

						if(nN1 == nMaxN1)
						{
							nMaxN1Pos2 = nj+nW2;
						}
						else if(nN1 > nMaxN1)
						{
							nMaxN1 = nN1;
							nMaxN1Pos1 = nj+nW2;
							nMaxN1Pos2 = nMaxN1Pos1;
						}

						if(nN2 == nMaxN2)
						{
							nMaxN2Pos2 = nj+nW2;
						}
						else if(nN2 > nMaxN2)
						{
							nMaxN2 = nN2;
							nMaxN2Pos1 = nj+nW2;
							nMaxN2Pos2 = nMaxN2Pos1;
						}
					}
				}
				else
				{
					if(nStart >= 0)
					{
						if(nPosSide == 1)
						{
							fprintf(fpReg, "1\t%d\t%d\t%d\t%f\t%d\t%f\t%d\t%d\t%d\t%d\t%d\n", ni1, nStart, nEnd, 
								dMinFDR, (nMinFDRPos1+nMinFDRPos2)/2,
								dMaxFC, (nMaxFCPos1+nMaxFCPos2)/2,
								nMaxN1, (nMaxN1Pos1+nMaxN1Pos2)/2,
								nMaxN2, (nMaxN2Pos1+nMaxN2Pos2)/2);
						}
						else
						{
							fprintf(fpReg, "1\t%d\t%d\t%d\t%f\t%d\t%f\t%d\t%d\t%d\t%d\t%d\n", ni1, nStart, nEnd, 
								dMinFDR, (nMinFDRPos1+nMinFDRPos2)/2,
								dMinFC, (nMinFCPos1+nMinFCPos2)/2,
								nMaxN1, (nMaxN1Pos1+nMaxN1Pos2)/2,
								nMaxN2, (nMaxN2Pos1+nMaxN2Pos2)/2);
						}
						nStart = -1;
						nEnd = -1;
						
						dMinFDR = 1.0;
						nMinFDRPos1 = -1;
						nMinFDRPos2 = -1;

						dMaxFC = 0.0;
						nMaxFCPos1 = -1;
						nMaxFCPos2 = -1;

						dMinFC = 0.0;
						nMinFCPos1 = -1;
						nMinFCPos2 = -1;

						nMaxN1 = 0;
						nMaxN1Pos1 = -1;
						nMaxN1Pos2 = -1;
						nMaxN2 = 0;
						nMaxN2Pos1 = -1;
						nMaxN2Pos2 = -1;
					}
				}
			}

			if(nStart >= 0)
			{
				if(nPosSide == 1)
				{
					fprintf(fpReg, "1\t%d\t%d\t%d\t%f\t%d\t%f\t%d\t%d\t%d\t%d\t%d\n", ni1, nStart, nEnd, 
						dMinFDR, (nMinFDRPos1+nMinFDRPos2)/2,
						dMaxFC, (nMaxFCPos1+nMaxFCPos2)/2,
						nMaxN1, (nMaxN1Pos1+nMaxN1Pos2)/2,
						nMaxN2, (nMaxN2Pos1+nMaxN2Pos2)/2);
				}
				else
				{
					fprintf(fpReg, "1\t%d\t%d\t%d\t%f\t%d\t%f\t%d\t%d\t%d\t%d\t%d\n", ni1, nStart, nEnd, 
						dMinFDR, (nMinFDRPos1+nMinFDRPos2)/2,
						dMinFC, (nMinFCPos1+nMinFCPos2)/2,
						nMaxN1, (nMaxN1Pos1+nMaxN1Pos2)/2,
						nMaxN2, (nMaxN2Pos1+nMaxN2Pos2)/2);
				}
				nStart = -1;
				nEnd = -1;
				dMinFDR = 1.0;
				nMinFDRPos1 = -1;
				nMinFDRPos2 = -1;

				dMaxFC = 0.0;
				nMaxFCPos1 = -1;
				nMaxFCPos2 = -1;

				dMinFC = 0.0;
				nMinFCPos1 = -1;
				nMinFCPos2 = -1;

				nMaxN1 = 0;
				nMaxN1Pos1 = -1;
				nMaxN1Pos2 = -1;
				nMaxN2 = 0;
				nMaxN2Pos1 = -1;
				nMaxN2Pos2 = -1;
			}

			/* add seqnum */
			ni1++;
		}

		/* process bar_neg only */
		else
		{
			if(pBARDataNeg->vSeqData[ni2]->nDataNum <= 0)
			{
				ni2++;
				continue;
			}

			nx = pBARDataNeg->vSeqData[ni2]->nDataNum-1;
			nP1 = (int)(pBARDataNeg->vSeqData[ni2]->vData[0]->pMatElement[0])-nW;
			nP2 = (int)(pBARDataNeg->vSeqData[ni2]->vData[0]->pMatElement[nx]);
			nStart = -1;
			nEnd = -1;

			dMinFDR = 1.0;
			nMinFDRPos1 = -1;
			nMinFDRPos2 = -1;

			dMaxFC = 0.0;
			nMaxFCPos1 = -1;
			nMaxFCPos2 = -1;

			dMinFC = 0.0;
			nMinFCPos1 = -1;
			nMinFCPos2 = -1;

			nMaxN1 = 0;
			nMaxN1Pos1 = -1;
			nMaxN1Pos2 = -1;
			nMaxN2 = 0;
			nMaxN2Pos1 = -1;
			nMaxN2Pos2 = -1;

			nN1 = 0;
			nN2 = 0;
			dFC = 0.0;
			nk = 0;
			nl = 0;
			for(nj=nP1; nj<=nP2; nj++)
			{
				for(; nl<pBARDataNeg->vSeqData[ni2]->nDataNum; nl++)
				{
					if((int)(pBARDataNeg->vSeqData[ni2]->vData[0]->pMatElement[nl]) >= nj)
						break;
					
					nN2 -= (int)(pBARDataNeg->vSeqData[ni2]->vData[1]->pMatElement[nl]);
				}

				for(; nk<pBARDataNeg->vSeqData[ni2]->nDataNum; nk++)
				{
					if((int)(pBARDataNeg->vSeqData[ni2]->vData[0]->pMatElement[nk]) >= (nj+nW))
						break;

					nN2 += (int)(pBARDataNeg->vSeqData[ni2]->vData[1]->pMatElement[nk]);
				}

				dTemp = nN2*dR0;
				if(dTemp < (double)nN1)
				{
					dFC = log((double)(nN1+1)/(dTemp+1.0))/dLog2;
				}
				else
				{
					dFC = -log((dTemp+1.0)/(double)(nN1+1))/dLog2;
				}

				nN = nN1+nN2;

				if( (nj+nW2)%nS == 0)
				{
					if( ((nj+nW2)>=0) && (nN>0) )
						fprintf(fpOut, "%s\t%d\t%d\t%d\t%f\n", pBARDataNeg->vSeqData[ni2]->pSeqName->m_pString, (int)(nj+nW2), nN1, nN2, dFC);
				}

				/* get fdr */
				nN1Sub = nN1;
				if(nN < pFDR->nHeight)
				{
					if(nN1 >= pFDR->nWidth)
						nN1Sub = pFDR->nWidth-1;

					dTemp = DMGETAT(pFDR, nN, nN1Sub);
				}
				else
				{
					nNSub = pFDR->nHeight-1;
					dTemp = (double)(nNSub*nN1)/(double)nN;
					nN1Sub = (int)dTemp;
					if(dTemp-nN1Sub > 0.5)
						nN1Sub += 1;
					
					if(nN1Sub >= pFDR->nWidth)
						nN1Sub = pFDR->nWidth-1;

					dTemp = DMGETAT(pFDR, nNSub, nN1Sub);
				}

				if( nOneSide == 1 )
				{
					if( dFC < 0.0 )
						dTemp = 1.0;
				}

				if( (nN >= nTCut) && (dTemp <= dFDRCut) )
				{
					if(nStart < 0)
					{
						nStart = nj;
						nEnd = nj+nW-1;

						dMinFDR = dTemp;
						nMinFDRPos1 = nj+nW2;
						nMinFDRPos2 = nMinFDRPos1;
						if(dFC >= 0.0)
							nPosSide = 1;
						else
							nPosSide = 0;

						dMaxFC = dFC;
						nMaxFCPos1 = nj+nW2;
						nMaxFCPos2 = nMaxFCPos1;

						dMinFC = dFC;
						nMinFCPos1 = nj+nW2;
						nMinFCPos2 = nMinFCPos1;

						nMaxN1 = nN1;
						nMaxN1Pos1 = nj+nW2;
						nMaxN1Pos2 = nMaxN1Pos1;
						nMaxN2 = nN2;
						nMaxN2Pos1 = nj+nW2;
						nMaxN2Pos2 = nMaxN2Pos1;
					}
					else
					{
						nEnd = nj+nW-1;

						if( fabs(dTemp-dMinFDR) < 1e-10 )
						{
							nMinFDRPos2 = nj+nW2;
						}
						else if(dTemp < dMinFDR)
						{
							dMinFDR = dTemp;
							nMinFDRPos1 = nj+nW2;
							nMinFDRPos2 = nMinFDRPos1;

							if(dFC >= 0.0)
								nPosSide = 1;
							else
								nPosSide = 0;
						}

						if( fabs(dFC-dMaxFC) < 1e-10 )
						{
							nMaxFCPos2 = nj+nW2;
						}
						else if(dFC > dMaxFC)
						{
							dMaxFC = dFC;
							nMaxFCPos1 = nj+nW2;
							nMaxFCPos2 = nMaxFCPos1;
						}

						if( fabs(dFC-dMinFC) < 1e-10 )
						{
							nMinFCPos2 = nj+nW2;
						}
						else if(dFC < dMinFC)
						{
							dMinFC = dFC;
							nMinFCPos1 = nj+nW2;
							nMinFCPos2 = nMinFCPos1;
						}

						if(nN1 == nMaxN1)
						{
							nMaxN1Pos2 = nj+nW2;
						}
						else if(nN1 > nMaxN1)
						{
							nMaxN1 = nN1;
							nMaxN1Pos1 = nj+nW2;
							nMaxN1Pos2 = nMaxN1Pos1;
						}

						if(nN2 == nMaxN2)
						{
							nMaxN2Pos2 = nj+nW2;
						}
						else if(nN2 > nMaxN2)
						{
							nMaxN2 = nN2;
							nMaxN2Pos1 = nj+nW2;
							nMaxN2Pos2 = nMaxN2Pos1;
						}
					}
				}
				else
				{
					if(nStart >= 0)
					{
						if(nPosSide == 1)
						{
							fprintf(fpReg, "2\t%d\t%d\t%d\t%f\t%d\t%f\t%d\t%d\t%d\t%d\t%d\n", ni2, nStart, nEnd, 
								dMinFDR, (nMinFDRPos1+nMinFDRPos2)/2,
								dMaxFC, (nMaxFCPos1+nMaxFCPos2)/2,
								nMaxN1, (nMaxN1Pos1+nMaxN1Pos2)/2,
								nMaxN2, (nMaxN2Pos1+nMaxN2Pos2)/2);
						}
						else
						{
							fprintf(fpReg, "2\t%d\t%d\t%d\t%f\t%d\t%f\t%d\t%d\t%d\t%d\t%d\n", ni2, nStart, nEnd, 
								dMinFDR, (nMinFDRPos1+nMinFDRPos2)/2,
								dMinFC, (nMinFCPos1+nMinFCPos2)/2,
								nMaxN1, (nMaxN1Pos1+nMaxN1Pos2)/2,
								nMaxN2, (nMaxN2Pos1+nMaxN2Pos2)/2);
						}
						nStart = -1;
						nEnd = -1;
						
						dMinFDR = 1.0;
						nMinFDRPos1 = -1;
						nMinFDRPos2 = -1;

						dMaxFC = 0.0;
						nMaxFCPos1 = -1;
						nMaxFCPos2 = -1;

						dMinFC = 0.0;
						nMinFCPos1 = -1;
						nMinFCPos2 = -1;

						nMaxN1 = 0;
						nMaxN1Pos1 = -1;
						nMaxN1Pos2 = -1;
						nMaxN2 = 0;
						nMaxN2Pos1 = -1;
						nMaxN2Pos2 = -1;
					}
				}
			}

			if(nStart >= 0)
			{
				if(nPosSide == 1)
				{
					fprintf(fpReg, "2\t%d\t%d\t%d\t%f\t%d\t%f\t%d\t%d\t%d\t%d\t%d\n", ni2, nStart, nEnd, 
						dMinFDR, (nMinFDRPos1+nMinFDRPos2)/2,
						dMaxFC, (nMaxFCPos1+nMaxFCPos2)/2,
						nMaxN1, (nMaxN1Pos1+nMaxN1Pos2)/2,
						nMaxN2, (nMaxN2Pos1+nMaxN2Pos2)/2);
				}
				else
				{
					fprintf(fpReg, "2\t%d\t%d\t%d\t%f\t%d\t%f\t%d\t%d\t%d\t%d\t%d\n", ni2, nStart, nEnd, 
						dMinFDR, (nMinFDRPos1+nMinFDRPos2)/2,
						dMinFC, (nMinFCPos1+nMinFCPos2)/2,
						nMaxN1, (nMaxN1Pos1+nMaxN1Pos2)/2,
						nMaxN2, (nMaxN2Pos1+nMaxN2Pos2)/2);
				}
				nStart = -1;
				nEnd = -1;
				dMinFDR = 1.0;
				nMinFDRPos1 = -1;
				nMinFDRPos2 = -1;

				dMaxFC = 0.0;
				nMaxFCPos1 = -1;
				nMaxFCPos2 = -1;

				dMinFC = 0.0;
				nMinFCPos1 = -1;
				nMinFCPos2 = -1;

				nMaxN1 = 0;
				nMaxN1Pos1 = -1;
				nMaxN1Pos2 = -1;
				nMaxN2 = 0;
				nMaxN2Pos1 = -1;
				nMaxN2Pos2 = -1;
			}

			ni2++;
		}
	}

	fclose(fpOut);
	fclose(fpReg);

	/* load region information */
	pRegion = DMLOAD(strRegTmpFile);

	/* convert temp file to bar file */
	sprintf(strFileName, "%s.cgw", strOutFileTitle);
	TileMapv2_TXT2BAR(strWinBarTmpFile, strExportFolder, strFileName);
	
	/* remove temp file */
	RemoveFiles(strWinBarTmpFile);
	RemoveFiles(strRegTmpFile);

	/* return */
	return pRegion;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_TwoSample_MergeRegion()                                     */
/*  Merge and filter regions.                                              */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *HTS_Enrich_TwoSample_MergeRegion(struct DOUBLEMATRIX *pRegion0, 
	int nMinLen, int nMaxGap)
{
	/* define */
	struct DOUBLEMATRIX *pMReg = NULL;
	struct DOUBLEMATRIX *pRegion = NULL;
	int nRegNum;
	int ni;
	int nBar0,nBar;
	int nChr0,nChr;
	int nStart0,nStart;
	int nEnd0,nEnd;
	double dMinFDR0,dMinFDR;
	int nMinFDRPos0,nMinFDRPos;
	double dMFC0,dMFC;
	int nMFCPos0,nMFCPos;
	int nMaxN10,nMaxN1;
	int nMaxN1Pos0,nMaxN1Pos;
	int nMaxN20,nMaxN2;
	int nMaxN2Pos0,nMaxN2Pos;
	int nPSide0,nPSide;
	int nIsNew;
	int nLen;

	/* init check */
	if(pRegion0 == NULL)
		return NULL;
	if(pRegion0->nHeight <= 0)
		return NULL;

	/* merge */
	pMReg = CreateDoubleMatrix(pRegion0->nHeight, pRegion0->nWidth+3);
	if(pMReg == NULL)
	{
		printf("Error: HTS_Enrich_TwoSample_MergeRegion, cannot create memory for merging regions!\n");
		exit(EXIT_FAILURE);
	}

	nRegNum = 0;
	nBar0 = (int)(DMGETAT(pRegion0, 0, 0));
	nChr0 = (int)(DMGETAT(pRegion0, 0, 1));
	nStart0 = (int)(DMGETAT(pRegion0, 0, 2));
	nEnd0 = (int)(DMGETAT(pRegion0, 0, 3));
	dMinFDR0 = DMGETAT(pRegion0, 0, 4);
	nMinFDRPos0 = (int)(DMGETAT(pRegion0, 0, 5));
	dMFC0 = DMGETAT(pRegion0, 0, 6);
	nMFCPos0 = (int)(DMGETAT(pRegion0, 0, 7));
	nMaxN10 = (int)(DMGETAT(pRegion0, 0, 8));
	nMaxN1Pos0 = (int)(DMGETAT(pRegion0, 0, 9));
	nMaxN20 = (int)(DMGETAT(pRegion0, 0, 10));
	nMaxN2Pos0 = (int)(DMGETAT(pRegion0, 0, 11));

	if(dMFC0 > 0.0)
		nPSide0 = 1;
	else
		nPSide0 = -1;

	for(ni=1; ni<pRegion0->nHeight; ni++)
	{
		nIsNew = 0;

		nBar = (int)(DMGETAT(pRegion0, ni, 0));
		nChr = (int)(DMGETAT(pRegion0, ni, 1));
		nStart = (int)(DMGETAT(pRegion0, ni, 2));
		nEnd = (int)(DMGETAT(pRegion0, ni, 3));
		dMinFDR = DMGETAT(pRegion0, ni, 4);
		nMinFDRPos = (int)(DMGETAT(pRegion0, ni, 5));
		dMFC = DMGETAT(pRegion0, ni, 6);
		nMFCPos = (int)(DMGETAT(pRegion0, ni, 7));
		nMaxN1 = (int)(DMGETAT(pRegion0, ni, 8));
		nMaxN1Pos = (int)(DMGETAT(pRegion0, ni, 9));
		nMaxN2 = (int)(DMGETAT(pRegion0, ni, 10));
		nMaxN2Pos = (int)(DMGETAT(pRegion0, ni, 11));

		if(dMFC > 0.0)
			nPSide = 1;
		else
			nPSide = -1;


		if( (nBar != nBar0) || (nChr != nChr0) )
		{
			nIsNew = 1;
		}
		else
		{
			if( (nStart-nEnd0) > nMaxGap)
				nIsNew = 1;
			else if( nPSide0 != nPSide )
				nIsNew = 1;
		}

		if(nIsNew == 1)
		{
			nLen = nEnd0-nStart0+1;
			if(nLen >= nMinLen)
			{
				DMSETAT(pMReg, nRegNum, 0, (double)nBar0);
				DMSETAT(pMReg, nRegNum, 1, (double)nChr0);
				DMSETAT(pMReg, nRegNum, 2, (double)nStart0);
				DMSETAT(pMReg, nRegNum, 3, (double)nEnd0);
				DMSETAT(pMReg, nRegNum, 4, (double)nLen);
				DMSETAT(pMReg, nRegNum, 5, -dMinFDR0);
				DMSETAT(pMReg, nRegNum, 6, (double)nMinFDRPos0);
				DMSETAT(pMReg, nRegNum, 7, dMFC0);
				DMSETAT(pMReg, nRegNum, 8, (double)nMFCPos0);
				DMSETAT(pMReg, nRegNum, 9, (double)nMaxN10);
				DMSETAT(pMReg, nRegNum, 10, (double)nMaxN1Pos0);
				DMSETAT(pMReg, nRegNum, 11, (double)nMaxN20);
				DMSETAT(pMReg, nRegNum, 12, (double)nMaxN2Pos0);
				nRegNum++;
			}

			nBar0 = nBar;
			nChr0 = nChr;
			nStart0 = nStart;
			nEnd0 = nEnd;
			dMinFDR0 = dMinFDR;
			nMinFDRPos0 = nMinFDRPos;
			dMFC0 = dMFC;
			nMFCPos0 = nMFCPos;
			nMaxN10 = nMaxN1;
			nMaxN1Pos0 = nMaxN1Pos;
			nMaxN20 = nMaxN2;
			nMaxN2Pos0 = nMaxN2Pos;
			nPSide0 = nPSide;
		}
		else
		{
			if(nEnd > nEnd0)
				nEnd0 = nEnd;

			if(dMinFDR < dMinFDR0)
			{
				dMinFDR0 = dMinFDR;
				nMinFDRPos0 = nMinFDRPos;
			}

			if(nPSide0 == 1)
			{
				if(dMFC > dMFC0)
				{
					dMFC0 = dMFC;
					nMFCPos0 = nMFCPos;
				}
			}
			else
			{
				if(dMFC < dMFC0)
				{
					dMFC0 = dMFC;
					nMFCPos0 = nMFCPos;
				}
			}

			if(nMaxN1 > nMaxN10)
			{
				nMaxN10 = nMaxN1;
				nMaxN1Pos0 = nMaxN1Pos;
			}

			if(nMaxN2 > nMaxN20)
			{
				nMaxN20 = nMaxN2;
				nMaxN2Pos0 = nMaxN2Pos;
			}
		}
	}

	nLen = nEnd0-nStart0+1;
	if(nLen >= nMinLen)
	{
		DMSETAT(pMReg, nRegNum, 0, (double)nBar0);
		DMSETAT(pMReg, nRegNum, 1, (double)nChr0);
		DMSETAT(pMReg, nRegNum, 2, (double)nStart0);
		DMSETAT(pMReg, nRegNum, 3, (double)nEnd0);
		DMSETAT(pMReg, nRegNum, 4, (double)nLen);
		DMSETAT(pMReg, nRegNum, 5, -dMinFDR0);
		DMSETAT(pMReg, nRegNum, 6, (double)nMinFDRPos0);
		DMSETAT(pMReg, nRegNum, 7, dMFC0);
		DMSETAT(pMReg, nRegNum, 8, (double)nMFCPos0);
		DMSETAT(pMReg, nRegNum, 9, (double)nMaxN10);
		DMSETAT(pMReg, nRegNum, 10, (double)nMaxN1Pos0);
		DMSETAT(pMReg, nRegNum, 11, (double)nMaxN20);
		DMSETAT(pMReg, nRegNum, 12, (double)nMaxN2Pos0);
		nRegNum++;
	}


	/* prepare the regions */
	pRegion = CreateDoubleMatrix(nRegNum, pMReg->nWidth);
	if(pRegion == NULL)
	{
		printf("Error: HTS_Enrich_TwoSample_MergeRegion, cannot create memory for merging regions!\n");
		exit(EXIT_FAILURE);
	}
	memcpy(pRegion->pMatElement, pMReg->pMatElement, (sizeof(double)*nRegNum*pMReg->nWidth));

	DestroyDoubleMatrix(pMReg);

	/* return */
	return pRegion;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_TwoSample_ExportResults()                                   */
/*  Export results.                                                        */
/* ----------------------------------------------------------------------- */ 
int	HTS_Enrich_TwoSample_ExportResults(struct DOUBLEMATRIX *pRegion, 
			struct tagBARData *pBARDataPos, struct tagBARData *pBARDataNeg,
			char strFileName[])
{
	/* define */
	FILE *fpOut;
	int ni,nj,nk,nl;

	/* open file */
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: HTS_Enrich_TwoSample_ExportResults, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "#rank\tchr\tstart\tend\tstrand\tlength\tminFDR\tminFDR_pos\tmax|FC|\tmax|FC|_pos\tpos_read_num\tneg_read_num\n");

	/* write */
	for(nj=0; nj<pRegion->nHeight; nj++)
	{
		ni = pRegion->nHeight-1-nj; 
		nk = (int)(DMGETAT(pRegion, ni,0));
		nl = (int)(DMGETAT(pRegion, ni,1));
		if(nk <= 1)
		{
			fprintf(fpOut, "%d\t%s\t%d\t%d\t+\t%d\t%f\t%d\t%f\t%d\t%d\t%d\n", nj+1, pBARDataPos->vSeqData[nl]->pSeqName->m_pString,
				(int)(DMGETAT(pRegion, ni,2)), (int)(DMGETAT(pRegion, ni,3)), (int)(DMGETAT(pRegion, ni,4)),
				-DMGETAT(pRegion, ni,5), (int)(DMGETAT(pRegion, ni,6)), DMGETAT(pRegion, ni,7), (int)(DMGETAT(pRegion, ni,8)),
				(int)DMGETAT(pRegion, ni,9), (int)(DMGETAT(pRegion, ni,10)));
		}
		else
		{
			fprintf(fpOut, "%d\t%s\t%d\t%d\t+\t%d\t%f\t%d\t%f\t%d\t%d\t%d\n", nj+1, pBARDataNeg->vSeqData[nl]->pSeqName->m_pString,
				(int)(DMGETAT(pRegion, ni,2)), (int)(DMGETAT(pRegion, ni,3)), (int)(DMGETAT(pRegion, ni,4)),
				-DMGETAT(pRegion, ni,5), (int)(DMGETAT(pRegion, ni,6)), DMGETAT(pRegion, ni,7), (int)(DMGETAT(pRegion, ni,8)),
				(int)DMGETAT(pRegion, ni,9), (int)(DMGETAT(pRegion, ni,10)));
		}
	}

	/* close file */
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_TwoSamplev2_ExportResults()                                 */
/*  Export results.                                                        */
/* ----------------------------------------------------------------------- */ 
int	HTS_Enrich_TwoSamplev2_ExportResults(struct DOUBLEMATRIX *pRegion, 
			struct tagBARData *pBARDataPos, struct tagBARData *pBARDataNeg,
			char strFileName[], int nBR, int nBRL, int nSSF, int nSSFF, int nSSFR,
			double dFC, double dTFC)
{
	/* define */
	FILE *fpOut;
	int ni,nj,nk,nl,nx,ny,nz;
	double dMin,dMax,dTemp;
	int nP1,nP2,nB1,nB2,nMF,nMR,nBRW;
	int nRegLen1,nRegLen2;
	double df;

	/* open file */
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: HTS_Enrich_TwoSample_ExportResults, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "#rank\tchr\tstart\tend\tstrand\tlength\tminFDR\tminFDR_pos\tmax|FC|\tmax|FC|_pos\tpos_read_num\tneg_read_num\tFminFDR\tFminFDR_Pos\tFmax|FC|\tFmax|FC|_pos\tFpos_readnum\tFneg_readnum\tFmaxpos_readnum\tFmaxpos_readnum_pos\tFmaxneg_readnum\tFmaxneg_readnum_pos\tRminFDR\tRminFDR_Pos\tRmax|FC|\tRmax|FC|_pos\tRpos_readnum\tRneg_readnum\tRmaxpos_readnum\tRmaxpos_readnum_pos\tRmaxneg_readnum\tRmaxneg_readnum_pos\tRmode-Fmode\tDelta\tmaxPosN\tmaxPosN_pos\tmaxNegN\tmaxNegN_pos\n");

	/* write */
	nz = 0;
	for(nj=0; nj<pRegion->nHeight; nj++)
	{
		ni = pRegion->nHeight-1-nj; 
		nk = (int)(DMGETAT(pRegion, ni,0));
		nl = (int)(DMGETAT(pRegion, ni,1));
		nP1 = (int)(DMGETAT(pRegion, ni,2));
		nP2 = (int)(DMGETAT(pRegion, ni,3));
		nRegLen1 = (int)(DMGETAT(pRegion, ni,4));
		nB1 = (int)(DMGETAT(pRegion, ni, 18));
		nB2 = (int)(DMGETAT(pRegion, ni, 28));
		nRegLen2 = nB2-nB1+1;
		nMF = (int)(DMGETAT(pRegion, ni, 21));
		nMR = (int)(DMGETAT(pRegion, ni, 31));
		df = DMGETAT(pRegion, ni,7);

		if(df < dFC)
			continue;
		if(DMGETAT(pRegion, ni,13) < dTFC*DMGETAT(pRegion, ni,14))
			continue;

		if(nSSF == 1)
		{
			if( (nMF < nSSFF) || (nMR < nSSFR) )
				continue;
		}
		if(nBR == 1)
		{
			if(nRegLen2 >= 0)
			{
				if(nRegLen2 < nBRL)
				{
					nBRW = (nBRL-nRegLen2)/2;
					nP1 = nB1-nBRW;
					nP2 = nB2+nBRW;
				}
				else
				{
					nP1 = nB1;
					nP2 = nB2;
				}
				nRegLen1 = nP2-nP1+1;
			}
		}

		if(nk <= 1)
		{
			fprintf(fpOut, "%d\t%s\t%d\t%d\t+\t%d\t%f\t%d\t%f\t%d\t%d\t%d\t", nz+1, pBARDataPos->vSeqData[nl]->pSeqName->m_pString,
				nP1, nP2, nRegLen1,
				-DMGETAT(pRegion, ni,5), (int)(DMGETAT(pRegion, ni,6)), DMGETAT(pRegion, ni,7), (int)(DMGETAT(pRegion, ni,8)),
				(int)DMGETAT(pRegion, ni,13), (int)(DMGETAT(pRegion, ni,14)));
		}
		else
		{
			fprintf(fpOut, "%d\t%s\t%d\t%d\t+\t%d\t%f\t%d\t%f\t%d\t%d\t%d\t", nz+1, pBARDataNeg->vSeqData[nl]->pSeqName->m_pString,
				nP1, nP2, nRegLen1,
				-DMGETAT(pRegion, ni,5), (int)(DMGETAT(pRegion, ni,6)), DMGETAT(pRegion, ni,7), (int)(DMGETAT(pRegion, ni,8)),
				(int)DMGETAT(pRegion, ni,13), (int)(DMGETAT(pRegion, ni,14)));
		}

		ny = 15;
		for(nx=0; nx<2; nx++)
		{
			fprintf(fpOut, "%f\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t", (DMGETAT(pRegion, ni, ny+4)), (int)(DMGETAT(pRegion, ni, ny+5)),
				(DMGETAT(pRegion, ni, ny+2)), (int)(DMGETAT(pRegion, ni, ny+3)), 
				(int)(DMGETAT(pRegion, ni, ny)), (int)(DMGETAT(pRegion, ni, ny+1)),
				(int)(DMGETAT(pRegion, ni, ny+6)), (int)(DMGETAT(pRegion, ni, ny+7)),
				(int)(DMGETAT(pRegion, ni, ny+8)), (int)(DMGETAT(pRegion, ni, ny+9)) );
			ny += 10;
		}

		fprintf(fpOut, "%d\t", nRegLen2);

		dMax = DMGETAT(pRegion, ni, 17);
		dMin = DMGETAT(pRegion, ni, 27);
		if(dMin > dMax)
		{
			dTemp = dMax;
			dMax = dMin;
			dMin = dTemp;
		}

		dTemp = pow(2.0, dMin)/pow(2.0, dMax);

		fprintf(fpOut, "%f\t", dTemp);

		fprintf(fpOut, "%d\t%d\t%d\t%d\n", (int)(DMGETAT(pRegion, ni,9)), (int)(DMGETAT(pRegion, ni,10)),
				(int)DMGETAT(pRegion, ni,11), (int)(DMGETAT(pRegion, ni,12)));

		nz++;
	}

	/* close file */
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_TwoSample_RegionCollectInfo()                               */
/*  Collect reads number.                                                  */
/* ----------------------------------------------------------------------- */ 
int HTS_Enrich_TwoSample_RegionCollectInfo(struct DOUBLEMATRIX *pRegion, 
		struct tagBARData *pBARDataPos, struct tagBARData *pBARDataNeg)
{
	/* define */
	int ni,nj,nk,nl;
	int nBar,nChr,nStart,nEnd;
	int nP1,nP2;
	char strChr[MED_LINE_LENGTH];
	int nPosChr,nNegChr;

	/* init */
	if( (pRegion == NULL) || (pBARDataPos == NULL) || (pBARDataNeg == NULL) )
		return PROC_SUCCESS;

	/* collect */
	for(ni=0; ni<pRegion->nHeight; ni++)
	{
		nBar = (int)(DMGETAT(pRegion, ni, 0));
		nChr = (int)(DMGETAT(pRegion, ni, 1));
		nStart = (int)(DMGETAT(pRegion, ni, 2));
		nEnd = (int)(DMGETAT(pRegion, ni, 3));

		if(nBar == 0)
		{
			nPosChr = nChr;
			strcpy(strChr, pBARDataPos->vSeqData[nChr]->pSeqName->m_pString);
			for(nj=0; nj<pBARDataNeg->nSeqNum; nj++)
			{
				if(strcmp(strChr, pBARDataNeg->vSeqData[nj]->pSeqName->m_pString) == 0)
				{
					nNegChr = nj;
					break;
				}
			}
		}
		else if(nBar == 1)
		{
			nPosChr = nChr;
			nNegChr = -1; 
		}
		else
		{
			nPosChr = -1;
			nNegChr = nChr;
		}


		if(nPosChr >= 0)
		{
			if(pBARDataPos->vSeqData[nPosChr]->nDataNum > 0)
			{
		
				nj = 0;
				nk = pBARDataPos->vSeqData[nPosChr]->nDataNum-1;

				if(nStart > pBARDataPos->vSeqData[nPosChr]->vData[0]->pMatElement[nk])
				{
					nP1 = nk+1;
				}
				else if(nStart <= pBARDataPos->vSeqData[nPosChr]->vData[0]->pMatElement[0])
				{
					nP1 = 0;
				}
				else
				{
					while( (nk-nj) > 1)
					{
						nl = (nk+nj)/2;
						if( pBARDataPos->vSeqData[nPosChr]->vData[0]->pMatElement[nl] >= nStart)
						{
							nk = nl;
						}
						else
						{
							nj = nl;
						}
					}
					nP1 = nk;
				}

				nj = 0;
				nk = pBARDataPos->vSeqData[nPosChr]->nDataNum-1;

				if(nEnd >= pBARDataPos->vSeqData[nPosChr]->vData[0]->pMatElement[nk])
				{
					nP2 = nk;
				}
				else if(nEnd < pBARDataPos->vSeqData[nPosChr]->vData[0]->pMatElement[0])
				{
					nP2 = -1;
				}
				else
				{
					while( (nk-nj) > 1)
					{
						nl = (nk+nj)/2;
						if( pBARDataPos->vSeqData[nPosChr]->vData[0]->pMatElement[nl] > nEnd)
						{
							nk = nl;
						}
						else
						{
							nj = nl;
						}
					}
					nP2 = nj;
				}

				nk = 0;
				for(nj=nP1; nj<=nP2; nj++)
				{
					nk += (int)(pBARDataPos->vSeqData[nPosChr]->vData[1]->pMatElement[nj]);
				}

				DMSETAT(pRegion, ni, 13, (double)nk);
			}
		}

		if(nNegChr >= 0)
		{
			if(pBARDataNeg->vSeqData[nNegChr]->nDataNum > 0)
			{
		
				nj = 0;
				nk = pBARDataNeg->vSeqData[nNegChr]->nDataNum-1;

				if(nStart > pBARDataNeg->vSeqData[nNegChr]->vData[0]->pMatElement[nk])
				{
					nP1 = nk+1;
				}
				else if(nStart <= pBARDataNeg->vSeqData[nNegChr]->vData[0]->pMatElement[0])
				{
					nP1 = 0;
				}
				else
				{
					while( (nk-nj) > 1)
					{
						nl = (nk+nj)/2;
						if( pBARDataNeg->vSeqData[nNegChr]->vData[0]->pMatElement[nl] >= nStart)
						{
							nk = nl;
						}
						else
						{
							nj = nl;
						}
					}
					nP1 = nk;
				}

				nj = 0;
				nk = pBARDataNeg->vSeqData[nNegChr]->nDataNum-1;

				if(nEnd >= pBARDataNeg->vSeqData[nNegChr]->vData[0]->pMatElement[nk])
				{
					nP2 = nk;
				}
				else if(nEnd < pBARDataNeg->vSeqData[nNegChr]->vData[0]->pMatElement[0])
				{
					nP2 = -1;
				}
				else
				{
					while( (nk-nj) > 1)
					{
						nl = (nk+nj)/2;
						if( pBARDataNeg->vSeqData[nNegChr]->vData[0]->pMatElement[nl] > nEnd)
						{
							nk = nl;
						}
						else
						{
							nj = nl;
						}
					}
					nP2 = nj;
				}

				nk = 0;
				for(nj=nP1; nj<=nP2; nj++)
				{
					nk += (int)(pBARDataNeg->vSeqData[nNegChr]->vData[1]->pMatElement[nj]);
				}

				DMSETAT(pRegion, ni, 14, (double)nk);
			}
		}
	}


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Enrich_TwoSamplev2_RegionCollectInfo()                             */
/*  Collect reads number.                                                  */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *HTS_Enrich_TwoSamplev2_RegionCollectInfo(struct DOUBLEMATRIX *pRegion, 
		struct tagBARData *pBARDataPos, struct tagBARData *pBARDataNeg, int nW, 
		struct DOUBLEMATRIX *pFDR, double dP0, int nOneSide)
{
	/* define */
	struct DOUBLEMATRIX *pNewRegion = NULL;
	int ni,nj,nk,nl,nk1,nl1,nk2,nl2,nu1,nu2;
	int nBar,nChr,nStart,nEnd;
	int nP1,nP2;
	char strChr[MED_LINE_LENGTH];
	int nPosChr,nNegChr;
	double dTemp;
	int nQ1,nQ2;
	double dR0 = dP0/(1.0-dP0);
	double dLog2 = log(2.0);
	double dMaxFC,dMinFC,dFC,dMinFDR,dFDR;
	int nMaxFCPos1,nMaxFCPos2,nMinFCPos1,nMinFCPos2,nMinFDRPos1,nMinFDRPos2;
	int nPMaxPos1,nPMaxPos2,nNMaxPos1,nNMaxPos2;
	int nPMax,nNMax;
	int nN1,nN2,nN,nNSub,nN1Sub;
	int nGlobalFCPos;
	int nW2 = nW/2;

	/* init */
	if( (pRegion == NULL) || (pBARDataPos == NULL) || (pBARDataNeg == NULL) )
		return NULL;

	pNewRegion = CreateDoubleMatrix(pRegion->nHeight, pRegion->nWidth+10);
	if(pNewRegion == NULL)
	{
		printf("Error: HTS_Enrich_TwoSamplev2_RegionCollectInfo, cannot create new region matrix.\n");
		exit(EXIT_FAILURE);
	}

	/* collect */
	for(ni=0; ni<pRegion->nHeight; ni++)
	{
		for(nj=0; nj<pRegion->nWidth; nj++)
		{
			dTemp = DMGETAT(pRegion, ni, nj);
			DMSETAT(pNewRegion, ni, nj, dTemp);
		}

		nBar = (int)(DMGETAT(pRegion, ni, 0));
		nChr = (int)(DMGETAT(pRegion, ni, 1));
		nStart = (int)(DMGETAT(pRegion, ni, 2));
		nEnd = (int)(DMGETAT(pRegion, ni, 3));

		if(DMGETAT(pRegion, ni, 7) >= 0.0)
			nGlobalFCPos = 1;
		else
			nGlobalFCPos = 0;

		if(nBar == 0)
		{
			nPosChr = nChr;
			strcpy(strChr, pBARDataPos->vSeqData[nChr]->pSeqName->m_pString);
			for(nj=0; nj<pBARDataNeg->nSeqNum; nj++)
			{
				if(strcmp(strChr, pBARDataNeg->vSeqData[nj]->pSeqName->m_pString) == 0)
				{
					nNegChr = nj;
					break;
				}
			}
		}
		else if(nBar == 1)
		{
			nPosChr = nChr;
			nNegChr = -1; 
		}
		else
		{
			nPosChr = -1;
			nNegChr = nChr;
		}

		nk1 = 0;
		nl1 = 0;
		nu1 = -1;
		nk2 = 0;
		nl2 = 0;
		nu2 = -1;

		if(nPosChr >= 0)
		{
			if(pBARDataPos->vSeqData[nPosChr]->nDataNum > 0)
			{
		
				nj = 0;
				nk = pBARDataPos->vSeqData[nPosChr]->nDataNum-1;

				if(nStart > pBARDataPos->vSeqData[nPosChr]->vData[0]->pMatElement[nk])
				{
					nP1 = nk+1;
				}
				else if(nStart <= pBARDataPos->vSeqData[nPosChr]->vData[0]->pMatElement[0])
				{
					nP1 = 0;
				}
				else
				{
					while( (nk-nj) > 1)
					{
						nl = (nk+nj)/2;
						if( pBARDataPos->vSeqData[nPosChr]->vData[0]->pMatElement[nl] >= nStart)
						{
							nk = nl;
						}
						else
						{
							nj = nl;
						}
					}
					nP1 = nk;
				}

				nj = 0;
				nk = pBARDataPos->vSeqData[nPosChr]->nDataNum-1;

				if(nEnd >= pBARDataPos->vSeqData[nPosChr]->vData[0]->pMatElement[nk])
				{
					nP2 = nk;
				}
				else if(nEnd < pBARDataPos->vSeqData[nPosChr]->vData[0]->pMatElement[0])
				{
					nP2 = -1;
				}
				else
				{
					while( (nk-nj) > 1)
					{
						nl = (nk+nj)/2;
						if( pBARDataPos->vSeqData[nPosChr]->vData[0]->pMatElement[nl] > nEnd)
						{
							nk = nl;
						}
						else
						{
							nj = nl;
						}
					}
					nP2 = nj;
				}

				nk = 0;
				for(nj=nP1; nj<=nP2; nj++)
				{
					nk += (int)(pBARDataPos->vSeqData[nPosChr]->vData[1]->pMatElement[nj]);
				}

				DMSETAT(pNewRegion, ni, pRegion->nWidth, (double)nk);
				nk1 = nP1;
				nl1 = nP1;
				nu1 = nP2;
			}
		}

		if(nNegChr >= 0)
		{
			if(pBARDataNeg->vSeqData[nNegChr]->nDataNum > 0)
			{
		
				nj = 0;
				nk = pBARDataNeg->vSeqData[nNegChr]->nDataNum-1;

				if(nStart > pBARDataNeg->vSeqData[nNegChr]->vData[0]->pMatElement[nk])
				{
					nP1 = nk+1;
				}
				else if(nStart <= pBARDataNeg->vSeqData[nNegChr]->vData[0]->pMatElement[0])
				{
					nP1 = 0;
				}
				else
				{
					while( (nk-nj) > 1)
					{
						nl = (nk+nj)/2;
						if( pBARDataNeg->vSeqData[nNegChr]->vData[0]->pMatElement[nl] >= nStart)
						{
							nk = nl;
						}
						else
						{
							nj = nl;
						}
					}
					nP1 = nk;
				}

				nj = 0;
				nk = pBARDataNeg->vSeqData[nNegChr]->nDataNum-1;

				if(nEnd >= pBARDataNeg->vSeqData[nNegChr]->vData[0]->pMatElement[nk])
				{
					nP2 = nk;
				}
				else if(nEnd < pBARDataNeg->vSeqData[nNegChr]->vData[0]->pMatElement[0])
				{
					nP2 = -1;
				}
				else
				{
					while( (nk-nj) > 1)
					{
						nl = (nk+nj)/2;
						if( pBARDataNeg->vSeqData[nNegChr]->vData[0]->pMatElement[nl] > nEnd)
						{
							nk = nl;
						}
						else
						{
							nj = nl;
						}
					}
					nP2 = nj;
				}

				nk = 0;
				for(nj=nP1; nj<=nP2; nj++)
				{
					nk += (int)(pBARDataNeg->vSeqData[nNegChr]->vData[1]->pMatElement[nj]);
				}

				DMSETAT(pNewRegion, ni, (pRegion->nWidth+1), (double)nk);
				nk2 = nP1;
				nl2 = nP1;
				nu2 = nP2;
			}
		}

		nQ1 = nStart;
		nQ2 = nEnd-nW+1;

		dMinFDR = 10.0;
		nMinFDRPos1 = -1;
		nMinFDRPos2 = -1;

		dMaxFC = -1e20;
		nMaxFCPos1 = -1;
		nMaxFCPos2 = -1;

		dMinFC = 1e20;
		nMinFCPos1 = -1;
		nMinFCPos2 = -1;

		nPMax = -1;
		nNMax = -1;
		nPMaxPos1 = -1;
		nPMaxPos2 = -1;
		nNMaxPos1 = -1;
		nNMaxPos2 = -1;


		nN1 = 0;
		nN2 = 0;
		dFC = 0.0;

		for(nj=nQ1; nj<=nQ2; nj++)
		{
			if(nPosChr >= 0)
			{
				for(; nl1<=nu1; nl1++)
				{
					if((int)(pBARDataPos->vSeqData[nPosChr]->vData[0]->pMatElement[nl1]) >= nj)
						break;
					
					nN1 -= (int)(pBARDataPos->vSeqData[nPosChr]->vData[1]->pMatElement[nl1]);
				}

				for(; nk1<=nu1; nk1++)
				{
					if((int)(pBARDataPos->vSeqData[nPosChr]->vData[0]->pMatElement[nk1]) >= (nj+nW))
						break;

					nN1 += (int)(pBARDataPos->vSeqData[nPosChr]->vData[1]->pMatElement[nk1]);
				}
			}

			if(nNegChr >= 0)
			{
				for(; nl2<=nu2; nl2++)
				{
					if((int)(pBARDataNeg->vSeqData[nNegChr]->vData[0]->pMatElement[nl2]) >= nj)
						break;
					
					nN2 -= (int)(pBARDataNeg->vSeqData[nNegChr]->vData[1]->pMatElement[nl2]);
				}

				for(; nk2<=nu2; nk2++)
				{
					if((int)(pBARDataNeg->vSeqData[nNegChr]->vData[0]->pMatElement[nk2]) >= (nj+nW))
						break;

					nN2 += (int)(pBARDataNeg->vSeqData[nNegChr]->vData[1]->pMatElement[nk2]);
				}
			}

			if( nPMax == nN1 )
			{
				nPMaxPos2 = nj+nW2;
			}
			else if( nN1 > nPMax )
			{
				nPMax = nN1;
				nPMaxPos1 = nj+nW2;
				nPMaxPos2 = nPMaxPos1;
			}

			if( nNMax == nN2 )
			{
				nNMaxPos2 = nj+nW2;
			}
			else if( nN2 > nNMax )
			{
				nNMax = nN2;
				nNMaxPos1 = nj+nW2;
				nNMaxPos2 = nPMaxPos1;
			}

			dTemp = nN2*dR0;
			dFC = log((double)(nN1+1)/(dTemp+1.0))/dLog2;
			
			if( fabs(dFC-dMaxFC) < 1e-10 )
			{
				nMaxFCPos2 = nj+nW2;
			}
			else if(dFC > dMaxFC)
			{
				dMaxFC = dFC;
				nMaxFCPos1 = nj+nW2;
				nMaxFCPos2 = nMaxFCPos1;
			}

			if( fabs(dFC-dMinFC) < 1e-10 )
			{
				nMinFCPos2 = nj+nW2;
			}
			else if(dFC < dMinFC)
			{
				dMinFC = dFC;
				nMinFCPos1 = nj+nW2;
				nMinFCPos2 = nMinFCPos1;
			}


			/* get fdr */
			if( ((nGlobalFCPos == 1) && (dFC>=0.0)) || ((nGlobalFCPos == 0) && (dFC<0.0)) )
			{
				nN = nN1+nN2;
				nN1Sub = nN1;
				if(nN < pFDR->nHeight)
				{
					if(nN1 >= pFDR->nWidth)
						nN1Sub = pFDR->nWidth-1;

					dFDR = DMGETAT(pFDR, nN, nN1Sub);
				}
				else
				{
					nNSub = pFDR->nHeight-1;
					dTemp = (double)(nNSub*nN1)/(double)nN;
					nN1Sub = (int)dTemp;
					if(dTemp-nN1Sub > 0.5)
						nN1Sub += 1;
					
					if(nN1Sub >= pFDR->nWidth)
						nN1Sub = pFDR->nWidth-1;

					dFDR = DMGETAT(pFDR, nNSub, nN1Sub);
				}
			}
			else
			{
				dFDR = 1.0;
			}

			if( nOneSide == 1 )
			{
				if( dFC < 0.0 )
					dFDR = 1.0;
			}
		
			if( fabs(dFDR-dMinFDR) < 1e-10 )
			{
				nMinFDRPos2 = nj+nW2;
			}
			else if(dFDR < dMinFDR)
			{
				dMinFDR = dFDR;
				nMinFDRPos1 = nj+nW2;
				nMinFDRPos2 = nMinFDRPos1;
			}
		}
	
		if(nPMax < 0)
		{
			nPMax = 0;
			nPMaxPos1 = -1;
			nPMaxPos2 = -1;
		}

		if(nNMax < 0)
		{
			nNMax = 0;
			nNMaxPos1 = -1;
			nNMaxPos2 = -1;
		}

		if(dMaxFC < (-1e20+1.0))
		{
			dMaxFC = 0.0;
			nMaxFCPos1 = -1;
			nMaxFCPos2 = -1;
		}

		if(dMinFC > (1e20-1.0))
		{
			dMinFC = 0.0;
			nMinFCPos1 = -1;
			nMinFCPos2 = -1;
		}

		if(dMinFDR > (1.01))
		{
			dMinFDR = 1.0;
			nMinFDRPos1 = -1;
			nMinFDRPos2 = -1;
		}

		if( nOneSide == 1 )
		{
			DMSETAT(pNewRegion, ni, (pRegion->nWidth+2), dMaxFC);
			nk = (nMaxFCPos1+nMaxFCPos2)/2;
			DMSETAT(pNewRegion, ni, (pRegion->nWidth+3), (double)nk);
		}
		else
		{
			if(DMGETAT(pRegion, ni, 7) >= 0.0)
			{
				DMSETAT(pNewRegion, ni, (pRegion->nWidth+2), dMaxFC);
				nk = (nMaxFCPos1+nMaxFCPos2)/2;
				DMSETAT(pNewRegion, ni, (pRegion->nWidth+3), (double)nk);
			}
			else
			{
				DMSETAT(pNewRegion, ni, (pRegion->nWidth+2), dMinFC);
				nk = (nMinFCPos1+nMinFCPos2)/2;
				DMSETAT(pNewRegion, ni, (pRegion->nWidth+3), (double)nk);
			}
		}
		DMSETAT(pNewRegion, ni, (pRegion->nWidth+4), dMinFDR);
		nk = (nMinFDRPos1+nMinFDRPos2)/2;
		DMSETAT(pNewRegion, ni, (pRegion->nWidth+5), (double)nk);

		DMSETAT(pNewRegion, ni, (pRegion->nWidth+6), nPMax);
		nk = (nPMaxPos1+nPMaxPos2)/2;
		DMSETAT(pNewRegion, ni, (pRegion->nWidth+7), (double)nk);

		DMSETAT(pNewRegion, ni, (pRegion->nWidth+8), nNMax);
		nk = (nNMaxPos1+nNMaxPos2)/2;
		DMSETAT(pNewRegion, ni, (pRegion->nWidth+9), (double)nk);
	}

	/* return */
	return pNewRegion;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_AlnShift2BAR()                                                     */
/*  shift base pairs.                                                      */
/* ----------------------------------------------------------------------- */ 
int HTS_AlnShift2BAR(char strInputPath[], char strOutputPath[], int nS)
{
	/* define */
	char strFileName[MED_LINE_LENGTH];
	struct tagBARData *pBARDataF,*pBARDataR,*pBARData;
	char strBARTmpFileName[MED_LINE_LENGTH];
	FILE *fpIn,*fpOut;
	int nSeqNum = 0;
	int ni,nfi,nri,nj;
	int nk1,nk2;
	int nProbeNum = 0;
	int nCPos1,nCPos2;

	int nFieldNum = 0;
	int nColNum = 0;
	struct tagString *pSeqInfo = NULL;
	
	char strLine[LONG_LINE_LENGTH];
	char strTemp[LONG_LINE_LENGTH];
	char *chp1,*chp2;
	char strChr[MED_LINE_LENGTH];
	char strLastChr[MED_LINE_LENGTH];
	int nLastPos;
	int nMaxPos,nPos,nAlnCount;
	struct INTMATRIX *pCol = NULL;
	int nk;


	/* load all */
	pBARData = NULL;
	pBARData = Affy_LoadBAR_Fast(strInputPath);
	if(pBARData == NULL)
	{
		printf("Error: HTS_AlnShift2BAR, cannot load raw data!\n");
		exit(EXIT_FAILURE);
	}

	/* load forward */
	sprintf(strFileName, "%s_F.bar", strInputPath);
	pBARDataF = NULL;
	pBARDataF = Affy_LoadBAR_Fast(strFileName);
	if(pBARDataF == NULL)
	{
		printf("Error: HTS_AlnShift2BAR, cannot load raw data!\n");
		exit(EXIT_FAILURE);
	}

	/* load reverse */
	sprintf(strFileName, "%s_R.bar", strInputPath);
	pBARDataR = NULL;
	pBARDataR = Affy_LoadBAR_Fast(strFileName);
	if(pBARDataR == NULL)
	{
		printf("Error: HTS_AlnShift2BAR, cannot load raw data!\n");
		exit(EXIT_FAILURE);
	}

	sprintf(strBARTmpFileName, "%s.tmp", strOutputPath);
	fpOut = NULL;
	fpOut = fopen(strBARTmpFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: HTS_AlnShift2BAR, cannot open temporary output file!\n");
		exit(EXIT_FAILURE);
	}

	nSeqNum = pBARData->nSeqNum;
	for(ni=0; ni<nSeqNum; ni++)
	{
		nj = pBARData->vSeqData[ni]->nDataNum-1;
		if(nj >= 0)
			nMaxPos = (int)(pBARData->vSeqData[ni]->vData[0]->pMatElement[nj]);
		else
			nMaxPos = 0;

		for(nfi=0; nfi<pBARDataF->nSeqNum; nfi++)
		{
			if(strcmp(pBARData->vSeqData[ni]->pSeqName->m_pString, pBARDataF->vSeqData[nfi]->pSeqName->m_pString) == 0)
			{
				break;
			}
		}

		for(nri=0; nri<pBARDataR->nSeqNum; nri++)
		{
			if(strcmp(pBARData->vSeqData[ni]->pSeqName->m_pString, pBARDataR->vSeqData[nri]->pSeqName->m_pString) == 0)
			{
				break;
			}
		}

		nProbeNum = 0;
		if( (nfi>=pBARDataF->nSeqNum) && (nri>=pBARDataR->nSeqNum) )
		{
			printf("Warning: empty sequence %s\n", pBARData->vSeqData[ni]->pSeqName->m_pString);
		}
		else if(nfi>=pBARDataF->nSeqNum)
		{
			nAlnCount = 0;
			for(nj=0; nj<pBARDataR->vSeqData[nri]->nDataNum; nj++)
			{
				if(pBARDataR->vSeqData[nri]->vData[0]->pMatElement[nj] > nS)
					break;

				nAlnCount += (int)(pBARDataR->vSeqData[nri]->vData[1]->pMatElement[nj]);
			}

			if(nAlnCount > 0)
			{
				fprintf(fpOut, "%s\t0\t%d\n", pBARData->vSeqData[ni]->pSeqName->m_pString, nAlnCount);
				nProbeNum = 1;
			}

			for(; nj<pBARDataR->vSeqData[nri]->nDataNum; nj++)
			{
				nPos = (int)(pBARDataR->vSeqData[nri]->vData[0]->pMatElement[nj])-nS;
				fprintf(fpOut, "%s\t%d\t%d\n", pBARData->vSeqData[ni]->pSeqName->m_pString, nPos, pBARDataR->vSeqData[nri]->vData[1]->pMatElement[nj]);
				nProbeNum++;
			}
		}
		else if(nri>=pBARDataR->nSeqNum)
		{
			for(nj=0; nj<pBARDataF->vSeqData[nfi]->nDataNum; nj++)
			{
				nPos = (int)(pBARDataF->vSeqData[nfi]->vData[0]->pMatElement[nj])+nS;
				if(nPos >= nMaxPos)
				{
					break;
				}

				fprintf(fpOut, "%s\t%d\t%d\n", pBARData->vSeqData[ni]->pSeqName->m_pString, nPos, pBARDataF->vSeqData[nfi]->vData[1]->pMatElement[nj]);
				nProbeNum++;
			}

			nAlnCount = 0;
			for(; nj<pBARDataF->vSeqData[nfi]->nDataNum; nj++)
			{
				nAlnCount += (int)(pBARDataF->vSeqData[nfi]->vData[1]->pMatElement[nj]);
			}
			if(nAlnCount > 0)
			{
				fprintf(fpOut, "%s\t%d\t%d\n", pBARData->vSeqData[ni]->pSeqName->m_pString, nMaxPos, nAlnCount);
				nProbeNum++;
			}
		}
		else
		{
			nk1 = 0;
			nk2 = 0;
			
			nAlnCount = 0;
			for(nk2=0; nk2<pBARDataR->vSeqData[nri]->nDataNum; nk2++)
			{
				if(pBARDataR->vSeqData[nri]->vData[0]->pMatElement[nk2] >= nS)
					break;

				nAlnCount += (int)(pBARDataR->vSeqData[nri]->vData[1]->pMatElement[nk2]);
			}

			nLastPos = 0;
			while( (nk1<pBARDataF->vSeqData[nfi]->nDataNum) || (nk2<pBARDataR->vSeqData[nri]->nDataNum) )
			{
				for(; nk1<pBARDataF->vSeqData[nfi]->nDataNum; nk1++)
				{
					nCPos1 = (int)(pBARDataF->vSeqData[nfi]->vData[0]->pMatElement[nk1])+nS;
					if(nCPos1 > nMaxPos)
						nCPos1 = nMaxPos;

					if(nCPos1 > nLastPos)
						break;

					nAlnCount += (int)(pBARDataF->vSeqData[nfi]->vData[1]->pMatElement[nk1]);
				}
				if(nk1 == pBARDataF->vSeqData[nfi]->nDataNum)
					nCPos1 = nMaxPos+1;
				

				for(; nk2<pBARDataR->vSeqData[nri]->nDataNum; nk2++)
				{
					nCPos2 = (int)(pBARDataR->vSeqData[nri]->vData[0]->pMatElement[nk2])-nS;
					if(nCPos2 > nMaxPos)
						nCPos2 = nMaxPos;

					if(nCPos2 > nLastPos)
						break;

					nAlnCount += (int)(pBARDataR->vSeqData[nri]->vData[1]->pMatElement[nk2]);
				}
				if(nk2 == pBARDataR->vSeqData[nri]->nDataNum)
					nCPos2 = nMaxPos+1;

				if(nAlnCount > 0)
				{
					fprintf(fpOut, "%s\t%d\t%d\n", pBARData->vSeqData[ni]->pSeqName->m_pString, nLastPos, nAlnCount);
					nProbeNum++;
				}

				nAlnCount = 0;
				if(nCPos1 < nCPos2)
					nLastPos = nCPos1;
				else
					nLastPos = nCPos2;
			}
		}

		sprintf(strTemp, "%s\t%d\t", pBARData->vSeqData[ni]->pSeqName->m_pString, nProbeNum);
		StringAddTail(&pSeqInfo, strTemp);
	}


	fclose(fpOut);

	Affy_BARData_Destroy(&pBARData);
	Affy_BARData_Destroy(&pBARDataF);
	Affy_BARData_Destroy(&pBARDataR);


	/* load input */
	nColNum = 1;
	nFieldNum = 1;

	/* load head */
	if(nSeqNum <= 0)
	{
		printf("Warning: empty sequences!\n");
		return PROC_FAILURE;
	}

	/* create BAR object */
	pBARData = NULL;
	pBARData = Affy_BARData_Create();
	if(pBARData == NULL)
	{
		printf("Error: HTS_AlnShift2BAR, BARData object was not created!\n");
		exit(EXIT_FAILURE);
	}
	strcpy(pBARData->strMagicnumber, "barr\r\n\032\n");
	pBARData->fVersionnumber = 2.0;
    pBARData->nSeqNum = nSeqNum;
	pBARData->nColNum = nColNum+1;
	pBARData->nParamNum = 0;
	pBARData->vParamName = NULL;
	pBARData->vParamValue = NULL;
	pBARData->pFieldType = CreateIntMatrix(1, pBARData->nColNum);
	if(pBARData->pFieldType == NULL)
	{
		printf("Error: HTS_AlnShift2BAR, cannot allocate memory for field type.\n");
		exit(EXIT_FAILURE);
	}
	pBARData->pFieldType->pMatElement[0] = 2;
	for(ni=1; ni<pBARData->nColNum; ni++)
		pBARData->pFieldType->pMatElement[ni] = 1;

	pBARData->vSeqData = (struct tagBARSeq **)calloc(pBARData->nSeqNum, sizeof(struct tagBARSeq *));
	if(pBARData->vSeqData == NULL)
	{
		printf("Error: HTS_AlnShift2BAR, cannot allocate memory for loading sequence data.\n");
		exit(EXIT_FAILURE);
	}

	chp1 = pSeqInfo->m_pString;
	for(ni=0; ni<pBARData->nSeqNum; ni++)
	{
		/* create BARSeq object */
		pBARData->vSeqData[ni] = Affy_BARSeq_Create();
		if(pBARData->vSeqData[ni] == NULL)
		{
			printf("Error: HTS_AlnShift2BAR, cannot create BARSeq object.\n");
			exit(EXIT_FAILURE);
		}

		pBARData->vSeqData[ni]->pSeqGroupName = NULL;
		pBARData->vSeqData[ni]->pSeqVersion = NULL;
		pBARData->vSeqData[ni]->nParamNum = 0;
		pBARData->vSeqData[ni]->vParamName = NULL;
		pBARData->vSeqData[ni]->vParamValue = NULL;

		pBARData->vSeqData[ni]->nColNum = pBARData->nColNum;

		chp2 = strchr(chp1, '\t');
		*chp2 = '\0';
		strcpy(strTemp, chp1);
		StringAddTail(&(pBARData->vSeqData[ni]->pSeqName), strTemp);
		chp1 = chp2+1;

		chp2 = strchr(chp1, '\t');
		*chp2 = '\0';
		strcpy(strTemp, chp1);
		pBARData->vSeqData[ni]->nDataNum = atoi(strTemp);
		chp1 = chp2+1;


		pBARData->vSeqData[ni]->vData = (struct DOUBLEMATRIX **)calloc(pBARData->vSeqData[ni]->nColNum, sizeof(struct DOUBLEMATRIX *));
		if(pBARData->vSeqData[ni]->vData == NULL)
		{
			printf("Error: HTS_AlnShift2BAR, cannot allocate memory for loading bpmap coordinate data!\n");
			exit(EXIT_FAILURE);
		}
		if(pBARData->vSeqData[ni]->nDataNum > 0)
		{
			for(nj=0; nj<pBARData->vSeqData[ni]->nColNum; nj++)
			{
				pBARData->vSeqData[ni]->vData[nj] = CreateDoubleMatrix(1, pBARData->vSeqData[ni]->nDataNum );
				if(pBARData->vSeqData[ni]->vData[nj] == NULL)
				{
					printf("Error: HTS_AlnShift2BAR, cannot allocate memory for loading bpmap coordinate data!\n");
					exit(EXIT_FAILURE);
				}
			}
		}
	}
	
	/* load data */
	fpIn = NULL;
	fpIn = fopen(strBARTmpFileName, "r");
	if(fpIn == NULL)
	{
		printf("Error: HTS_AlnShift2BAR, cannot open source file!\n");
		exit(EXIT_FAILURE);
	}

	strcpy(strLastChr, "");
	nSeqNum = -1;
	nProbeNum = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '0')
			continue;
		if(strLine[0] == '#')
			continue;

		chp1 = strLine;
		chp2 = strchr(chp1, '\t');
		if(chp2 == NULL)
		{
			printf("Error: HTS_AlnShift2BAR, wrong input file format!\n");
			exit(EXIT_FAILURE);
		}
		*chp2 = '\0';
		StrTrimRight(chp1);
		strcpy(strChr, chp1);

		if(strcmp(strChr, strLastChr) != 0)
		{
			if(strcmp(strLastChr, "") != 0)
			{
				if(nProbeNum != (pBARData->vSeqData[nSeqNum]->nDataNum-1))
				{
					printf("Error: HTS_Aln2BAR, cannot read input file correctly!");
					exit(EXIT_FAILURE);
				}
			}

			nSeqNum++;
			strcpy(strLastChr, strChr);
			nProbeNum = 0;
		}
		else
		{
			nProbeNum++;
		}

		chp1 = chp2+1;
		chp2 = strchr(chp1, '\t');
		nk = 0;
		while(chp2 != NULL)
		{
			if(nk >= pBARData->vSeqData[nSeqNum]->nColNum)
			{
				printf("Error: HTS_AlnShift2BAR, input file format error, column number inconsistent!");
				exit(EXIT_FAILURE);
			}

			*chp2 = '\0';
			StrTrimLeft(chp1);
			StrTrimRight(chp1);
			pBARData->vSeqData[nSeqNum]->vData[nk]->pMatElement[nProbeNum] = atof(chp1);
			nk++;

			chp1 = chp2+1;
			chp2 = strchr(chp1, '\t');
		}

		if(nk >= pBARData->vSeqData[nSeqNum]->nColNum)
		{
			printf("Error: HTS_AlnShift2BAR, input file format error, column number inconsistent!");
			exit(EXIT_FAILURE);
		}
		
		StrTrimLeft(chp1);
		StrTrimRight(chp1);
		pBARData->vSeqData[nSeqNum]->vData[nk]->pMatElement[nProbeNum] = atof(chp1);
		nk++;

		if(nk != pBARData->vSeqData[nSeqNum]->nColNum)
		{
			printf("Error: HTS_AlnShift2BAR, input file format error, column number inconsistent!");
			exit(EXIT_FAILURE);
		}
	}

	
	/* close file */
	fclose(fpIn);

	if(nProbeNum != (pBARData->vSeqData[nSeqNum]->nDataNum-1))
	{
		printf("Error: HTS_AlnShift2BAR, cannot read input file correctly!");
		exit(EXIT_FAILURE);
	}

	if(nSeqNum != pBARData->nSeqNum-1)
	{
		printf("Error: HTS_AlnShift2BAR, cannot read input file correctly!");
		exit(EXIT_FAILURE);
	}
	
	/* write to BAR file */
	pCol = NULL;
	pCol = CreateIntMatrix(1,pBARData->nColNum);
	if(pCol == NULL)
	{
		printf("Error: HTS_AlnShift2BAR, cannot create output column information!\n");
		exit(EXIT_FAILURE);
	}
	
	pCol->pMatElement[0] = 1;
	pCol->pMatElement[1] = 1;
	Affy_SaveBAR_Columns_Fast(strOutputPath, pBARData, pCol);
	DestroyIntMatrix(pCol);

	
	/* release memory */
	DeleteString(pSeqInfo);
	pSeqInfo = NULL;

	Affy_BARData_Destroy(&pBARData);

	RemoveFiles(strBARTmpFileName);
	
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_WindowSummaryPaper()                                               */
/*  Summarize window counts.                                               */
/* ----------------------------------------------------------------------- */ 
int HTS_WindowSummaryPaper(char strBARFile[], char strChrList[], char strChrLen[], 
					  int nW, char strOutFile[], int nCombineShift)
{
	/* define */
	char strFileName[MED_LINE_LENGTH];
	FILE *fpIn;
	FILE *fpOut;
	char strLine[LONG_LINE_LENGTH];
	struct tagBARData *pBARData = NULL;
	struct INTMATRIX *pChrLen = NULL;
	struct DOUBLEMATRIX *pCount = NULL;
	int *vC;
	int ni,nj,nk,nl,nLen;


	/* initial check */
	if(nW <= 0)
	{
		nW = 100;
		printf("Warning: Window size<=0! Proceed with default window size = 100\n");
	}

	/* load chromosome length */
	pChrLen = IMLOAD(strChrLen);
	if(pChrLen == NULL)
	{
		printf("Error: HTS_WindowSummary, cannot load chromosome length!\n");
		exit(EXIT_FAILURE);
	}


	/* load bar data */
	pBARData = NULL;
	if(nCombineShift == 1)
	{
		sprintf(strFileName, "%s_C.bar", strBARFile);
		pBARData = Affy_LoadBAR_Fast(strFileName);
	}
	else
	{
		pBARData = Affy_LoadBAR_Fast(strBARFile);
	}
	if(pBARData == NULL)
	{
		printf("Error: HTS_WindowSummary, cannot load raw data!\n");
		exit(EXIT_FAILURE);
	}


	/* summarize results */
	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: HTS_WindowSummaryv2, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	
	/* process chromosome one by one */
	fpIn = NULL;
	fpIn = fopen(strChrList, "r");
	if(fpIn == NULL)
	{
		printf("Error: HTS_WindowSummary, cannot open chrlist file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		nLen = pChrLen->pMatElement[ni]/nW;
		if(pChrLen->pMatElement[ni]%nW != 0)
			nLen += 1;

		/* prepare memory */
		vC = NULL;
		vC = (int *)calloc(nLen, sizeof(int));
		if(vC == NULL)
		{
			printf("Error: HTS_WindowSummary, cannot create memory for counting!\n");
			exit(EXIT_FAILURE);
		}

		/* find matching chromosome */
		for(nj=0; nj<pBARData->nSeqNum; nj++)
		{
			if(strcmp(strLine, pBARData->vSeqData[nj]->pSeqName->m_pString) != 0)
				continue;

			for(nk=0; nk<pBARData->vSeqData[nj]->nDataNum; nk++)
			{
				nl = (int)(pBARData->vSeqData[nj]->vData[0]->pMatElement[nk])/nW;
				if(nl == nLen)
					nl = nLen-1;
				else if(nl > nLen)
				{
					printf("Error: HTS_WindowSummary, index out of range! Please check if the genomes are matching\n");
					exit(EXIT_FAILURE);
				}

				vC[nl] += (int)(pBARData->vSeqData[nj]->vData[1]->pMatElement[nk]);
			}
		}


		/* count */
		for(nj=0; nj<nLen; nj++)
		{
			fprintf(fpOut, "%d\n", vC[nj]);
		}

		/* free memory */
		free(vC);

		ni++;
	}
	
	fclose(fpIn);
	fclose(fpOut);

	/* destroy memory */
	Affy_BARData_Destroy(&pBARData);
	DestroyIntMatrix(pChrLen);
	DestroyDoubleMatrix(pCount);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_WindowSummaryv2()                                                  */
/*  Summarize window counts.                                               */
/* ----------------------------------------------------------------------- */ 
int HTS_CreateRepeatFilter(char strGenomePath[], char strChrList[], char strChrLen[], 
					 int nW, double dR, char strExportFolder[], char strOutFileTitle[])
{
	/* define */
	FILE *fpIn;
	FILE *fpSeq;
	FILE *fpOut;
	struct INTMATRIX *pChrLen;
	char strLine[MED_LINE_LENGTH];
	char strSeqFile[MED_LINE_LENGTH];
	int nFileLen;
	int nChr;
	int ni,nj,nk;
	int nWinNum = 0;
	int nM = 1000000;
	int nWN = 0;
	unsigned char *vB;
	unsigned char nMask1 = 192;
	unsigned char nMask2 = 12;
	unsigned char nTemp;
	int nRepeatN;
	double dRepeatR;
	char strOutFile[MED_LINE_LENGTH];
	int numread;

	/* init */
	AdjustDirectoryPath(strGenomePath);
	AdjustDirectoryPath(strExportFolder);
	nWN = nM/nW;
	if(nWN == 0)
		nWN = 1;
	nM = (int)(nWN*nW);
	vB = NULL;
	vB = (unsigned char *)calloc(nM, sizeof(unsigned char));
	if(vB == NULL)
	{
		printf("Error: cannot allocate memory!\n");
		exit(EXIT_FAILURE);
	}

	/* load chromosome length */
	pChrLen = NULL;
	pChrLen = IMLOAD(strChrLen);
	if(pChrLen == NULL)
	{
		printf("Error: cannot load chromosome length!\n");
		exit(EXIT_FAILURE);
	}

	/* open file */
	sprintf(strOutFile, "%s%s.bar.tmp", strExportFolder, strOutFileTitle);
	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	fprintf(fpOut, "#chr\tpos\t%s\n", strOutFileTitle);
	fprintf(fpOut, "#chr\tpos\t1\n");

	fpIn = NULL;
	fpIn = fopen(strChrList, "r");
	if(fpIn == NULL)
	{
		printf("Error: cannot open input file!\n");
		exit(EXIT_FAILURE);
	}

	/* process chromosome by chromosome */
	nChr = 0;
	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		sprintf(strSeqFile, "%s%s.sq", strGenomePath, strLine);
		
		fpSeq = NULL;
		fpSeq = fopen(strSeqFile, "rb");
		if(fpSeq == NULL)
		{
			printf("Error: cannot open sequence file!\n");
			exit(EXIT_FAILURE);
		}

		if(nChr >= pChrLen->nHeight)
		{
			printf("Error: chromosome number not match!\n");
			exit(EXIT_FAILURE);
		}
		nFileLen = pChrLen->pMatElement[nChr]/2;
		if(pChrLen->pMatElement[nChr]%2 != 0)
			nFileLen += 1;
		nWinNum = 0;

		ni = 0;
		while(ni<nFileLen)
		{
			nj = ni+nM-1;
			if(nj >= nFileLen)
				nj = nFileLen-1;

			numread = (int)(fread(vB, sizeof(unsigned char), (nj-ni+1), fpSeq));
			if(numread != (nj-ni+1))
			{
				printf("Error: File reading error!\n");
				exit(EXIT_FAILURE);
			}

			nk = 0;
			nRepeatN = 0;
			dRepeatR = 0.0;
			for(nj=0; nj<numread; nj++)
			{
				nTemp = vB[nj]&nMask1;
				if(nTemp>0)
					nRepeatN++;

				nk++;
				if(nk == nW)
				{
					dRepeatR = (double)nRepeatN/(double)nW;
					if(dRepeatR >= dR)
						fprintf(fpOut, "%s\t%d\t1\n", strLine, (int)(nWinNum*nW));
					else
						fprintf(fpOut, "%s\t%d\t0\n", strLine, (int)(nWinNum*nW));
					nk = 0;
					nRepeatN = 0;
					dRepeatR = 0.0;
					nWinNum++;
				}

				nTemp = vB[nj]&nMask2;
				if(nTemp>0)
					nRepeatN++;

				nk++;
				if(nk == nW)
				{
					dRepeatR = (double)nRepeatN/(double)nW;
					if(dRepeatR >= dR)
						fprintf(fpOut, "%s\t%d\t1\n", strLine, (int)(nWinNum*nW));
					else
						fprintf(fpOut, "%s\t%d\t0\n", strLine, (int)(nWinNum*nW));
					nk = 0;
					nRepeatN = 0;
					dRepeatR = 0.0;
					nWinNum++;
				}
			}

			if(nk != 0)
			{
				dRepeatR = (double)nRepeatN/(double)nW;
				if(dRepeatR >= dR)
					fprintf(fpOut, "%s\t%d\t1\n", strLine, (int)(nWinNum*nW));
				else
					fprintf(fpOut, "%s\t%d\t0\n", strLine, (int)(nWinNum*nW));
				nk = 0;
				nRepeatN = 0;
				dRepeatR = 0.0;
				nWinNum++;
			}

			ni += nM;
		}

		fclose(fpSeq);
		nChr++;
	}

	fclose(fpIn);
	fclose(fpOut);
	free(vB);
	DestroyIntMatrix(pChrLen);
	
	sprintf(strLine, "%s.cgw", strOutFileTitle);
	TileMapv2_TXT2BAR(strOutFile, strExportFolder, strLine);
	RemoveFiles(strOutFile);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_FilterRepeatReads()                                                */
/*  Filter repeat reads.                                                   */
/* ----------------------------------------------------------------------- */ 
int HTS_FilterRepeatReads(char strInputFile[], char strOutputFile[], char strGenomePath[])
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	FILE *fpSeq;
	char strLine[LONG_LINE_LENGTH];
	char strChr[MED_LINE_LENGTH];
	int nPos;
	char strSeqFile[MED_LINE_LENGTH];
	int ni,nj;
	unsigned char vB[1];
	unsigned char nMask1 = 192;
	unsigned char nMask2 = 12;
	unsigned char nTemp;

	/* init */
	AdjustDirectoryPath(strGenomePath);

	/* open files */
	fpIn = NULL;
	fpIn = fopen(strInputFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: HTS_FilterRepeatReads, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutputFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: HTS_FilterRepeatReads, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}

	/* read file */
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		sscanf(strLine, "%s %d", strChr, &nPos);
		sprintf(strSeqFile, "%s%s.sq", strGenomePath, strChr);

		fpSeq = NULL;
		fpSeq = fopen(strSeqFile, "rb");
		if(fpSeq == NULL)
			continue;

		ni = nPos/2;
		nj = nPos%2;

		if(fseek(fpSeq, ni, SEEK_SET) != 0)
		{
			fclose(fpSeq);
			continue;
		}

		if(fread(vB, sizeof(unsigned char), 1, fpSeq) != 1)
		{
			fclose(fpSeq);
			continue;
		}

		if(nj == 0)
		{
			nTemp = vB[0]&nMask1;
		}
		else
		{
			nTemp = vB[0]&nMask2;
		}

		if(nTemp == 0)
		{
			fprintf(fpOut, "%s\n", strLine);
		}

		fclose(fpSeq);
	}


	/* close files */
	fclose(fpIn);
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  HTS_CollectReads()                                                     */
/*  Collecte reads in specific regions.                                    */
/* ----------------------------------------------------------------------- */ 
int HTS_CollectReads(char strInputFile[], char strOutputFile[], char strBARFile[])
{
	/* define */
	struct tagBARData *pBARData;
	FILE *fpIn;
	FILE *fpOut;
	FILE *fpOut2;
	char strLine[LONG_LINE_LENGTH];
	char strID[MED_LINE_LENGTH];
	char strChr[MED_LINE_LENGTH];
	int nStart,nEnd;
	int ni,nj,nk,nx,nz;
	int nChr;
	int nP1,nP2,nP0;
	double dTotalCount = 0.0;
	double dRegionCount = 0.0;

	/* load bar data */
	pBARData = NULL;
	pBARData = Affy_LoadBAR_Fast(strBARFile);
	if(pBARData == NULL)
	{
		printf("Error: HTS_CollectReads, cannot load bar file!\n");
		exit(EXIT_FAILURE);
	}

	/* open file */
	fpIn = NULL;
	fpIn = fopen(strInputFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: HTS_CollectReads, cannot open input file!\n");
		exit(EXIT_FAILURE);	
	}

	fpOut = NULL;
	fpOut = fopen(strOutputFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: HTS_CollectReads, cannot open output file!\n");
		exit(EXIT_FAILURE);	
	}

	sprintf(strLine, "%s.stat", strOutputFile);
	fpOut2 = NULL;
	fpOut2 = fopen(strLine, "w");
	if(fpOut2 == NULL)
	{
		printf("Error: HTS_CollectReads, cannot open output file!\n");
		exit(EXIT_FAILURE);	
	}

	/* load data */
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		sscanf(strLine, "%s %s %d %d", strID, strChr, &nStart, &nEnd);
		dRegionCount = 0.0;

		nChr = -1;
		for(ni=0; ni<pBARData->nSeqNum; ni++)
		{
			if(strcmp(strChr, pBARData->vSeqData[ni]->pSeqName->m_pString) == 0)
			{
				nChr = ni;
				if(pBARData->vSeqData[nChr]->nDataNum > 0)
				{
					nP1 = 0;
					nP2 = pBARData->vSeqData[nChr]->nDataNum-1;
					
					if((int)(pBARData->vSeqData[nChr]->vData[0]->pMatElement[nP2]) < nStart)
						continue;
					if((int)(pBARData->vSeqData[nChr]->vData[0]->pMatElement[nP1]) > nEnd)
						continue;

					while(nP2-nP1 > 1)
					{
						nP0 = (nP1+nP2)/2;
						if(nStart <= (int)(pBARData->vSeqData[nChr]->vData[0]->pMatElement[nP0]))
							nP2 = nP0;
						else
							nP1 = nP0;
					}

					if(pBARData->vSeqData[nChr]->vData[0]->pMatElement[nP1] >= nStart)
						nj = nP1;
					else if(pBARData->vSeqData[nChr]->vData[0]->pMatElement[nP2] >= nStart)
						nj = nP2;

					nP1 = 0;
					nP2 = pBARData->vSeqData[nChr]->nDataNum-1;
					while(nP2-nP1 > 1)
					{
						nP0 = (nP1+nP2)/2;
						if(nEnd >= (int)(pBARData->vSeqData[nChr]->vData[0]->pMatElement[nP0]))
							nP1 = nP0;
						else
							nP2 = nP0;
					}

					if(pBARData->vSeqData[nChr]->vData[0]->pMatElement[nP2] <= nEnd)
						nk = nP2;
					else if(pBARData->vSeqData[nChr]->vData[0]->pMatElement[nP1] <= nEnd)
						nk = nP1;


					for(nz=nj; nz<=nk; nz++)
					{
						dRegionCount += pBARData->vSeqData[nChr]->vData[1]->pMatElement[nz];
						dTotalCount += pBARData->vSeqData[nChr]->vData[1]->pMatElement[nz];
						for(nx=0; nx<pBARData->vSeqData[nChr]->vData[1]->pMatElement[nz]; nx++)
						{
							fprintf(fpOut, "%s\t%d\tF\n", strChr, (int)(pBARData->vSeqData[nChr]->vData[0]->pMatElement[nz]));
						}
					}
				}
			}
		}
		fprintf(fpOut2, "%s\t%s\t%d\t%d\t+\t%d\n", strID, strChr, nStart, nEnd, (int)dRegionCount);
	}

	/* clear work space */
	fclose(fpIn);
	fclose(fpOut);
	fclose(fpOut2);
	Affy_BARData_Destroy(&pBARData);

	printf("Total Reads in the Target Regions = %d\n", (int)dTotalCount);

	/* return */
	return PROC_SUCCESS;

}


/* ----------------------------------------------------------------------- */ 
/*  HTS_CollectProfile_Main()                                                   */
/*  Collecte read profile in specific regions.                             */
/* ----------------------------------------------------------------------- */ 
int HTS_CollectProfile_Main(char strInputFile[], char strOutputFile[], char strBARFile[], int nW, int nS)
{
	/* define */
	char strFileName[MED_LINE_LENGTH];
	struct tagBARData *pBARDataF;
	struct tagBARData *pBARDataR;
	FILE *fpIn;
	FILE *fpOut;
	char strLine[LONG_LINE_LENGTH];
	char strAlias[MED_LINE_LENGTH];
	char strChr[MED_LINE_LENGTH];
	int nPos;
	struct INTMATRIX *pCF,*pCR;
	int ni;

	
	/* load bar data */
	sprintf(strFileName, "%s_F.bar", strBARFile);
	pBARDataF = NULL;
	pBARDataF = Affy_LoadBAR_Fast(strFileName);
	if(pBARDataF == NULL)
	{
		printf("Error: HTS_CollectProfile, cannot load bar file!\n");
		exit(EXIT_FAILURE);
	}

	sprintf(strFileName, "%s_R.bar", strBARFile);
	pBARDataR = NULL;
	pBARDataR = Affy_LoadBAR_Fast(strFileName);
	if(pBARDataR == NULL)
	{
		printf("Error: HTS_CollectProfile, cannot load bar file!\n");
		exit(EXIT_FAILURE);
	}

	/* open file */
	fpIn = NULL;
	fpIn = fopen(strInputFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: HTS_CollectProfile, cannot open input file!\n");
		exit(EXIT_FAILURE);	
	}

	fpOut = NULL;
	fpOut = fopen(strOutputFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: HTS_CollectProfile, cannot open output file!\n");
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

		sscanf(strLine, "%s %s %d", strAlias, strChr, &nPos);

		pCF = NULL;
		pCF = HTS_CollectProfile(strChr, nPos, nW, nS, pBARDataF);
		if(pCF == NULL)
		{
			printf("Error: cannot collect profile!\n");
			exit(EXIT_FAILURE);
		}

		pCR = NULL;
		pCR = HTS_CollectProfile(strChr, nPos, nW, nS, pBARDataR);
		if(pCR == NULL)
		{
			printf("Error: cannot collect profile!\n");
			exit(EXIT_FAILURE);
		}

		fprintf(fpOut, "%s\t%s\t%d", strAlias, strChr, nPos);
		for(ni=0; ni<pCF->nWidth; ni++)
		{
			fprintf(fpOut, "\t%d", pCF->pMatElement[ni]);
		}
		for(ni=0; ni<pCR->nWidth; ni++)
		{
			fprintf(fpOut, "\t%d", pCR->pMatElement[ni]);
		}
		fprintf(fpOut, "\n");

		DestroyIntMatrix(pCF);
		DestroyIntMatrix(pCR);
		pCF = NULL;
		pCR = NULL;

	}

	/* close file */
	fclose(fpIn);
	fclose(fpOut);

	/* release memory */
	Affy_BARData_Destroy(&pBARDataF);
	Affy_BARData_Destroy(&pBARDataR);

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  HTS_CollectProfile()                                                   */
/*  Collecte read profile in specific regions.                             */
/* ----------------------------------------------------------------------- */ 
struct INTMATRIX* HTS_CollectProfile(char strChr[], int nPos, int nW, int nS, struct tagBARData *pBARData)
{
	/* define */
	int ni,nj,nk,nz,nx;
	struct INTMATRIX *pC;
	int nL,nE;
	int nChr,nStart,nEnd;
	int nP0,nP1,nP2;
	double dDist;

	/* init */
	nL = nW/nS;
	if(nW%nS != 0)
		nL++;
	nE = (int)(nL*nS);
	nL = nL;
	pC = NULL;
	pC = CreateIntMatrix(1, (int)(2*nL));
	if(pC == NULL)
	{
		printf("Error: HTS_CollectProfile, cannot create counting matrix!\n");
		exit(EXIT_FAILURE);
	}
	nStart = nPos-nE;
	nEnd = nPos+nE;

	/* find boundary */
	for(ni=0; ni<pBARData->nSeqNum; ni++)
	{
		if(strcmp(strChr, pBARData->vSeqData[ni]->pSeqName->m_pString) == 0)
		{
			nChr = ni;
			if(pBARData->vSeqData[nChr]->nDataNum > 0)
			{
				nP1 = 0;
				nP2 = pBARData->vSeqData[nChr]->nDataNum-1;
				
				if((int)(pBARData->vSeqData[nChr]->vData[0]->pMatElement[nP2]) < nStart)
					continue;
				if((int)(pBARData->vSeqData[nChr]->vData[0]->pMatElement[nP1]) > nEnd)
					continue;

				while(nP2-nP1 > 1)
				{
					nP0 = (nP1+nP2)/2;
					if(nStart <= (int)(pBARData->vSeqData[nChr]->vData[0]->pMatElement[nP0]))
						nP2 = nP0;
					else
						nP1 = nP0;
				}

				if(pBARData->vSeqData[nChr]->vData[0]->pMatElement[nP1] >= nStart)
					nj = nP1;
				else if(pBARData->vSeqData[nChr]->vData[0]->pMatElement[nP2] >= nStart)
					nj = nP2;

				nP1 = 0;
				nP2 = pBARData->vSeqData[nChr]->nDataNum-1;
				while(nP2-nP1 > 1)
				{
					nP0 = (nP1+nP2)/2;
					if(nEnd >= (int)(pBARData->vSeqData[nChr]->vData[0]->pMatElement[nP0]))
						nP1 = nP0;
					else
						nP2 = nP0;
				}

				if(pBARData->vSeqData[nChr]->vData[0]->pMatElement[nP2] <= nEnd)
					nk = nP2;
				else if(pBARData->vSeqData[nChr]->vData[0]->pMatElement[nP1] <= nEnd)
					nk = nP1;


				for(nz=nj; nz<=nk; nz++)
				{
					dDist = pBARData->vSeqData[nChr]->vData[0]->pMatElement[nz]-nPos;
					if(dDist <= 0)
					{
						nx = (int)(-dDist)/nS;
						if(nx < nL)
						{
							pC->pMatElement[nL-1-nx] += (int)(pBARData->vSeqData[nChr]->vData[1]->pMatElement[nz]);
						}
					}
					else
					{
						nx = (int)(dDist)/nS;
						if(nx < nL)
						{
							pC->pMatElement[nL+nx] += (int)(pBARData->vSeqData[nChr]->vData[1]->pMatElement[nz]);
						}
					}					
				}
			}
		}
	}

	/* return */
	return pC;
}

/* ----------------------------------------------------------------------- */ 
/*  CNV_Aln2Window_Main()                                                  */
/*  Convert copy number sequencing alignment to windowed bar tiling        */
/*  data.                                                                  */
/* ----------------------------------------------------------------------- */ 
int CNV_Aln2Window_Main(char strInFile[], char strOutFile[], 
						char strChrListFile[], char strChrLenFile[], 
						int nW, int nL, int nN)
{	
	/* Sample name */
	struct INTMATRIX *pChrLen = NULL;
	struct tagBARData **vBARData = NULL;
	int ni,nj,nx,ny,nId,nBARNum;
	char strFileName[LONG_LINE_LENGTH];
	
	struct INTMATRIX *pCol;
	int *vN;
	int nCW;
	
	FILE *fpIn;
	char strLine[LONG_LINE_LENGTH];
	char strChr[LINE_LENGTH];
	int nPos;
	char strStrand[LINE_LENGTH];
	

	/* Load Chromosome */
	pChrLen = IMLOAD(strChrLenFile);
	if(pChrLen == NULL)
	{
		printf("Error: CNV_Aln2Window_Main, cannot load chromosome length!\n");
		exit(EXIT_FAILURE);
	}

	/* prepare BAR object */
	vBARData = NULL;
	nBARNum = (int)(nN*2);
	vBARData = (struct tagBARData **)calloc(nBARNum, sizeof(struct tagBARData *));
	if(vBARData == NULL)
	{
		printf("Error: CNV_Aln2Window_Main, cannot create bar data vector!\n");
		exit(EXIT_FAILURE);
	}

	vN = (int *)calloc(nN, sizeof(int));
	if(vN == NULL)
	{
		printf("Error: CNV_Aln2Window_Main, cannot create window size vector!\n");
		exit(EXIT_FAILURE);
	}

	nj = 0;
	nCW = nW;
	for(ni=0; ni<nN; ni++)
	{
		vBARData[nj] = CNV_Aln2Window_PrepareBARData(nCW, pChrLen, strChrListFile);
		if(vBARData[nj] == NULL)
		{
			printf("Error: CNV_Aln2Window_Main, cannot create bar data object!\n");
			exit(EXIT_FAILURE);
		}
		nj++;

		vBARData[nj] = CNV_Aln2Window_PrepareBARData(nCW, pChrLen, strChrListFile);
		if(vBARData[nj] == NULL)
		{
			printf("Error: CNV_Aln2Window_Main, cannot create bar data object!\n");
			exit(EXIT_FAILURE);
		}
		nj++;

		vN[ni] = nCW;
		nCW = (int)(nCW*nL);
	}

	/* process data one by one */
	fpIn = fopen(strInFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: CNV_Aln2Window_Main, cannot open input file!\n");
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

		sscanf(strLine, "%s %d %s", strChr, &nPos, strStrand);
		
		for(nx=0; nx<nN; nx++)
		{
			if( (strcmp(strStrand, "-") == 0) || (strcmp(strStrand, "R") == 0) )
				ny = (int)(2*nx+1);
			else
				ny = (int)(2*nx);

			for(ni=0; ni<vBARData[ny]->nSeqNum; ni++)
			{
				if(strcmp(strChr, vBARData[ny]->vSeqData[ni]->pSeqName->m_pString) == 0)
				{
					nId = nPos/vN[nx];
					if(nId < vBARData[ny]->vSeqData[ni]->nDataNum)
					{
						vBARData[ny]->vSeqData[ni]->vData[1]->pMatElement[nId] += 1;
					}
					break;
				}
			}
		}
	}

	fclose(fpIn);

	/* output */
	pCol = NULL;
	pCol = CreateIntMatrix(1,vBARData[0]->nColNum);
	if(pCol == NULL)
	{
		printf("Error: CNV_Aln2Window_Main, cannot create output column information!\n");
		exit(EXIT_FAILURE);
	}
	pCol->pMatElement[0] = 1;
	pCol->pMatElement[1] = 1;

	nj = 0;
	for(ni=0; ni<nN; ni++)
	{
		/* output */
		sprintf(strFileName, "%s_w%d_F.bar", strOutFile, vN[ni]);
		Affy_SaveBAR_Columns_Fast(strFileName, vBARData[nj], pCol);
		nj++;

		sprintf(strFileName, "%s_w%d_R.bar", strOutFile, vN[ni]);
		Affy_SaveBAR_Columns_Fast(strFileName, vBARData[nj], pCol);
		nj++;
	}
	DestroyIntMatrix(pCol);
		
	/* Destroy Paramter */
	DestroyIntMatrix(pChrLen);
	pChrLen = NULL;
	free(vN);

	nj = 0;
	for(ni=0; ni<nN; ni++)
	{
		Affy_BARData_Destroy(vBARData+nj);
		nj++;
		Affy_BARData_Destroy(vBARData+nj);
		nj++;
	}
	free(vBARData);
	
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  CNV_Aln2WindowC_Main()                                                 */
/*  Convert copy number sequencing alignment to windowed bar tiling        */
/*  data.                                                                  */
/* ----------------------------------------------------------------------- */ 
int CNV_Aln2WindowC_Main(char strInFile[], char strOutFile[], 
						char strChrListFile[], char strChrLenFile[], 
						int nW, int nL, int nN)
{	
	/* Sample name */
	struct INTMATRIX *pChrLen = NULL;
	struct tagBARData **vBARData = NULL;
	int ni,nx,nId,nBARNum;
	char strFileName[LONG_LINE_LENGTH];
	
	struct INTMATRIX *pCol;
	int *vN;
	int nCW;
	
	FILE *fpIn;
	char strLine[LONG_LINE_LENGTH];
	char strChr[LINE_LENGTH];
	int nPos;
	char strStrand[LINE_LENGTH];
	

	/* Load Chromosome */
	pChrLen = IMLOAD(strChrLenFile);
	if(pChrLen == NULL)
	{
		printf("Error: CNV_Aln2WindowC_Main, cannot load chromosome length!\n");
		exit(EXIT_FAILURE);
	}

	/* prepare BAR object */
	vBARData = NULL;
	nBARNum = nN;
	vBARData = (struct tagBARData **)calloc(nBARNum, sizeof(struct tagBARData *));
	if(vBARData == NULL)
	{
		printf("Error: CNV_Aln2WindowC_Main, cannot create bar data vector!\n");
		exit(EXIT_FAILURE);
	}

	vN = (int *)calloc(nN, sizeof(int));
	if(vN == NULL)
	{
		printf("Error: CNV_Aln2WindowC_Main, cannot create window size vector!\n");
		exit(EXIT_FAILURE);
	}

	nCW = nW;
	for(ni=0; ni<nN; ni++)
	{
		vBARData[ni] = CNV_Aln2Window_PrepareBARData(nCW, pChrLen, strChrListFile);
		if(vBARData[ni] == NULL)
		{
			printf("Error: CNV_Aln2WindowC_Main, cannot create bar data object!\n");
			exit(EXIT_FAILURE);
		}

		vN[ni] = nCW;
		nCW = (int)(nCW*nL);
	}

	/* process data one by one */
	fpIn = fopen(strInFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: CNV_Aln2WindowC_Main, cannot open input file!\n");
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

		sscanf(strLine, "%s %d %s", strChr, &nPos, strStrand);
		
		for(nx=0; nx<nN; nx++)
		{
			for(ni=0; ni<vBARData[nx]->nSeqNum; ni++)
			{
				if(strcmp(strChr, vBARData[nx]->vSeqData[ni]->pSeqName->m_pString) == 0)
				{
					nId = nPos/vN[nx];
					if(nId < vBARData[nx]->vSeqData[ni]->nDataNum)
					{
						vBARData[nx]->vSeqData[ni]->vData[1]->pMatElement[nId] += 1;
					}
					break;
				}
			}
		}
	}

	fclose(fpIn);

	/* output */
	pCol = NULL;
	pCol = CreateIntMatrix(1,vBARData[0]->nColNum);
	if(pCol == NULL)
	{
		printf("Error: CNV_Aln2WindowC_Main, cannot create output column information!\n");
		exit(EXIT_FAILURE);
	}
	pCol->pMatElement[0] = 1;
	pCol->pMatElement[1] = 1;

	for(ni=0; ni<nN; ni++)
	{
		/* output */
		sprintf(strFileName, "%s_w%d.bar", strOutFile, vN[ni]);
		Affy_SaveBAR_Columns_Fast(strFileName, vBARData[ni], pCol);
	}
	DestroyIntMatrix(pCol);
		
	/* Destroy Paramter */
	DestroyIntMatrix(pChrLen);
	pChrLen = NULL;
	free(vN);

	for(ni=0; ni<nN; ni++)
	{
		Affy_BARData_Destroy(vBARData+ni);
	}
	free(vBARData);
	
	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  CNV_Aln2Window_PrepareBARData()                                        */
/*  Prepare the bar data object.                                           */
/* ----------------------------------------------------------------------- */ 
struct tagBARData * CNV_Aln2Window_PrepareBARData(int nW, struct INTMATRIX *pChrLen, 
				char strChrListFile[])
{
	/* define */
	struct tagBARData *pBARData;
	int nSeqNum;
	int nColNum;
	int nDataNum;
	FILE *fpIn;
	char strLine[LONG_LINE_LENGTH];
	int ni,nj,nk;
	
	/* init */
	nSeqNum = pChrLen->nHeight;
	nColNum = 2;

	/* create */
	pBARData = Affy_BARData_Create();
	if(pBARData == NULL)
	{
		printf("Error: CNV_Aln2Window_PrepareBARData, BARData object was not created!\n");
		exit(EXIT_FAILURE);
	}
	strcpy(pBARData->strMagicnumber, "barr\r\n\032\n");
	pBARData->fVersionnumber = 2.0;
    pBARData->nSeqNum = nSeqNum;
	pBARData->nColNum = nColNum;
	pBARData->nParamNum = 0;
	pBARData->vParamName = NULL;
	pBARData->vParamValue = NULL;
	pBARData->pFieldType = CreateIntMatrix(1, pBARData->nColNum);
	if(pBARData->pFieldType == NULL)
	{
		printf("Error: CNV_Aln2Window_PrepareBARData, cannot allocate memory for field type.\n");
		exit(EXIT_FAILURE);
	}
	pBARData->pFieldType->pMatElement[0] = 2;
	for(ni=1; ni<pBARData->nColNum; ni++)
		pBARData->pFieldType->pMatElement[ni] = 1;

	pBARData->vSeqData = (struct tagBARSeq **)calloc(pBARData->nSeqNum, sizeof(struct tagBARSeq *));
	if(pBARData->vSeqData == NULL)
	{
		printf("Error: CNV_Aln2Window_PrepareBARData, cannot allocate memory for loading sequence data.\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<pBARData->nSeqNum; ni++)
	{
		/* create BARSeq object */
		pBARData->vSeqData[ni] = Affy_BARSeq_Create();
		if(pBARData->vSeqData[ni] == NULL)
		{
			printf("Error: CNV_Aln2Window_PrepareBARData, cannot create BARSeq object.\n");
			exit(EXIT_FAILURE);
		}

		pBARData->vSeqData[ni]->pSeqGroupName = NULL;
		pBARData->vSeqData[ni]->pSeqVersion = NULL;
		pBARData->vSeqData[ni]->nParamNum = 0;
		pBARData->vSeqData[ni]->vParamName = NULL;
		pBARData->vSeqData[ni]->vParamValue = NULL;

		pBARData->vSeqData[ni]->nColNum = pBARData->nColNum;

		nDataNum = pChrLen->pMatElement[ni]/nW;
		if(pChrLen->pMatElement[ni]%nW != 0)
			nDataNum++;

		pBARData->vSeqData[ni]->nDataNum = nDataNum;
		pBARData->vSeqData[ni]->vData = (struct DOUBLEMATRIX **)calloc(pBARData->vSeqData[ni]->nColNum, sizeof(struct DOUBLEMATRIX *));
		if(pBARData->vSeqData[ni]->vData == NULL)
		{
			printf("Error: CNV_Aln2Window_PrepareBARData, cannot allocate memory for loading bpmap coordinate data!\n");
			exit(EXIT_FAILURE);
		}
		if(pBARData->vSeqData[ni]->nDataNum > 0)
		{
			for(nj=0; nj<pBARData->vSeqData[ni]->nColNum; nj++)
			{
				pBARData->vSeqData[ni]->vData[nj] = CreateDoubleMatrix(1, pBARData->vSeqData[ni]->nDataNum );
				if(pBARData->vSeqData[ni]->vData[nj] == NULL)
				{
					printf("Error: CNV_Aln2Window_PrepareBARData, cannot allocate memory for loading bpmap coordinate data!\n");
					exit(EXIT_FAILURE);
				}
			}
		}

		nk = nW/2;
		for(nj=0; nj<nDataNum; nj++)
		{
			pBARData->vSeqData[ni]->vData[0]->pMatElement[nj] = nk;
			nk += nW;
		}
	}

	fpIn = NULL;
	fpIn = fopen(strChrListFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: CNV_Aln2Window_PrepareBARData, cannot open chrlist file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		if(ni >= nSeqNum)
		{
			printf("Error: CNV_Aln2Window_PrepareBARData, chromosome number wrong!\n");
			exit(EXIT_FAILURE);
		}

		StringAddTail(&(pBARData->vSeqData[ni]->pSeqName), strLine);
		ni++;
	}
	
	fclose(fpIn);

	if(ni != nSeqNum)
	{
		printf("Error: CNV_Aln2Window_PrepareBARData, chromosome number wrong!\n");
		exit(EXIT_FAILURE);
	}

	/* return */
	return pBARData;
}

/* ----------------------------------------------------------------------- */ 
/*  CNV_Repeat2Window_Main()                                               */
/*  Convert repeat percentage to windowed bar tiling data.                 */
/* ----------------------------------------------------------------------- */ 
int CNV_Repeat2Window_Main(char strGenomePath[], char strOutFile[], 
						char strChrListFile[], char strChrLenFile[], 
						int nW, int nL, int nN)
{	
	/* Sample name */
	struct INTMATRIX *pChrLen = NULL;
	struct tagBARData **vBARData = NULL;
	int ni,nj,nx,ny,nz,nw,nId;
	char strFileName[LONG_LINE_LENGTH];
	int nContigLen = 1000000;
	unsigned char *vBase;
	int nHChrLen,nlen,nOdd;
	
	struct INTMATRIX *pCol;
	int *vN;
	int *vI;
	int nCW;
	
	FILE *fpIn;
	FILE *fpChr;
	char strChr[LINE_LENGTH];
	int nChr;
	int nPos;
	unsigned char bBase,bChar;
	
	/* init */
	AdjustDirectoryPath(strGenomePath);

	/* Load Chromosome */
	pChrLen = IMLOAD(strChrLenFile);
	if(pChrLen == NULL)
	{
		printf("Error: CNV_Repeat2Window_Main, cannot load chromosome length!\n");
		exit(EXIT_FAILURE);
	}

	/* prepare BAR object */
	vBARData = NULL;
	vBARData = (struct tagBARData **)calloc(nN, sizeof(struct tagBARData *));
	if(vBARData == NULL)
	{
		printf("Error: CNV_Repeat2Window_Main, cannot create bar data vector!\n");
		exit(EXIT_FAILURE);
	}

	vN = (int *)calloc(nN, sizeof(int));
	if(vN == NULL)
	{
		printf("Error: CNV_Repeat2Window_Main, cannot create window size vector!\n");
		exit(EXIT_FAILURE);
	}

	vI = (int *)calloc(nN, sizeof(int));
	if(vI == NULL)
	{
		printf("Error: CNV_Repeat2Window_Main, cannot create bar file indicator vector!\n");
		exit(EXIT_FAILURE);
	}

	vBase = (unsigned char *)calloc(nContigLen, sizeof(unsigned char));
	if(vBase == NULL)
	{
		printf("Error: CNV_Repeat2Window_Main, cannot create memory for loading sequences!\n");
		exit(EXIT_FAILURE);
	}

	nCW = nW;
	for(ni=0; ni<nN; ni++)
	{
		vBARData[ni] = CNV_Aln2Window_PrepareBARData(nCW, pChrLen, strChrListFile);
		if(vBARData[ni] == NULL)
		{
			printf("Error: CNV_Repeat2Window_Main, cannot create bar data object!\n");
			exit(EXIT_FAILURE);
		}
		vN[ni] = nCW;
		nCW = (int)(nCW*nL);
	}

	/* process data one by one */
	fpIn = fopen(strChrListFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: CNV_Repeat2Window_Main, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}

	nChr = 0;
	while(fgets(strChr, LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strChr);
		StrTrimRight(strChr);
		if(strChr[0] == '\0')
			continue;
		if(strChr[0] == '#')
			continue;

		for(ni=0; ni<nN; ni++)
			vI[ni] = -1;

		for(ni=0; ni<nN; ni++)
		{
			for(nj=0; nj<vBARData[ni]->nSeqNum; nj++)
			{
				if(strcmp(strChr, vBARData[ni]->vSeqData[nj]->pSeqName->m_pString) == 0)
				{
					vI[ni] = nj;
					break;
				}
			}
		}

		sprintf(strFileName, "%s%s.sq", strGenomePath, strChr);
		fpChr = NULL;
		fpChr = fopen(strFileName, "rb");
		if(fpChr == NULL)
		{
			printf("Error: CNV_Repeat2Window_Main, cannot open sequence file!\n");
			exit(EXIT_FAILURE);
		}

		nHChrLen = (int)(pChrLen->pMatElement[nChr]/2);
		if(pChrLen->pMatElement[nChr]%2 != 0)
		{
			nHChrLen++;
			nOdd = 1;
		}
		else
		{
			nOdd = 0;
		}

		nPos = 0;
		nx = 0;
		while(nx < nHChrLen)
		{
			ny = nx+nContigLen-1;
			if(ny >= nHChrLen)
				ny = nHChrLen-1;
			nlen = ny-nx+1;

			if( fread(vBase, sizeof(unsigned char), nlen, fpChr) != nlen)
			{
				printf("Error: CNV_Repeat2Window_Main, cannot load sequence correctly!\n");
				exit(EXIT_FAILURE);
			}

			for(nz=0; nz<nlen; nz++)
			{
				bChar = vBase[nz] >> 4; 
				bBase = bChar & 0x0F;
				if((bBase > 3) && (nPos<pChrLen->pMatElement[nChr]))
				{
					for(ni=0; ni<nN; ni++)
					{
						nw = vI[ni];
						if(nw >= 0)
						{
							nId = nPos/vN[ni];
							if(nId < vBARData[ni]->vSeqData[nw]->nDataNum)
							{
								vBARData[ni]->vSeqData[nw]->vData[1]->pMatElement[nId] += 1;
							}
						}
					}
				}
				nPos++;

				bBase = vBase[nz] & 0x0F;
				if((bBase > 3) && (nPos<pChrLen->pMatElement[nChr]))
				{
					for(ni=0; ni<nN; ni++)
					{
						nw = vI[ni];
						if(nw >= 0)
						{
							nId = nPos/vN[ni];
							if(nId < vBARData[ni]->vSeqData[nw]->nDataNum)
							{
								vBARData[ni]->vSeqData[nw]->vData[1]->pMatElement[nId] += 1;
							}
						}
					}
				}
				nPos++;
			}

			nx = ny+1;
		}

		fclose(fpChr);

		if(nOdd == 0)
		{
			if(nPos != pChrLen->pMatElement[nChr])
			{
				printf("Error: CNV_Repeat2Window_Main, loaded length of sequence does not match chromosome length!\n");
				exit(EXIT_FAILURE);
			}
		}
		else
		{
			if((nPos-1) != pChrLen->pMatElement[nChr])
			{
				printf("Error: CNV_Repeat2Window_Main, loaded length of sequence does not match chromosome length!\n");
				exit(EXIT_FAILURE);
			}
		}


		nChr++;
	}

	fclose(fpIn);

	/* normalize */
	for(ni=0; ni<nN; ni++)
	{
		for(nj=0; nj<vBARData[ni]->nSeqNum; nj++)
		{
			nw = pChrLen->pMatElement[nj]%vN[ni];

			for(nx=0; nx<vBARData[ni]->vSeqData[nj]->nDataNum; nx++)
			{
				vBARData[ni]->vSeqData[nj]->vData[1]->pMatElement[nx] /= vN[ni];
			}
			if(nw > 0)
			{
				nx--;
				vBARData[ni]->vSeqData[nj]->vData[1]->pMatElement[nx] = vBARData[ni]->vSeqData[nj]->vData[1]->pMatElement[nx]*vN[ni]/nw;
			}
		}
	}

	/* output */
	pCol = NULL;
	pCol = CreateIntMatrix(1,vBARData[0]->nColNum);
	if(pCol == NULL)
	{
		printf("Error: CNV_Repeat2Window_Main, cannot create output column information!\n");
		exit(EXIT_FAILURE);
	}
	pCol->pMatElement[0] = 1;
	pCol->pMatElement[1] = 1;

	for(ni=0; ni<nN; ni++)
	{
		/* output */
		sprintf(strFileName, "%s_w%d.bar", strOutFile, vN[ni]);
		Affy_SaveBAR_Columns_Fast(strFileName, vBARData[ni], pCol);
	}
	DestroyIntMatrix(pCol);
		
	/* Destroy Paramter */
	DestroyIntMatrix(pChrLen);
	pChrLen = NULL;
	free(vN);
	free(vI);
	free(vBase);

	for(ni=0; ni<nN; ni++)
	{
		Affy_BARData_Destroy(vBARData+ni);
	}
	free(vBARData);
	
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RNASEQ_CountReadPerTranscript_Main()                                   */
/*  Count reads for all transcripts.                                       */
/*  nInputType = 0: input is a single file specified by strInFile.         */
/*  nInputType = 1: input is a list of files specified by strInFile.       */
/* ----------------------------------------------------------------------- */ 
int RNASEQ_CountReadPerTranscript_Main(char strInFile[], int nInputType,
						char strOutFile[], char strDatabaseFile[], 
						int nDatabaseType, char strSpecies[],
						int nStandardize)
{
	/* out */
	struct tagRefGene **vSourceRefGene = NULL;
	int nSourceRefNum = 0;
	int nSampleNum = 0;
	FILE *fpIn;
	FILE *fpOut;
	char strLine[LONG_LINE_LENGTH];
	struct DOUBLEMATRIX **vE;
	int nUp = 100;
	int ni,nj,nx,ny,nz,nu;
	struct tagString **vSampleName = NULL;
	double dT;
	char strChr[MED_LINE_LENGTH];
	int nPos;
	char strStrand[LINE_LENGTH];
	int nChr;
	int nGeneLen;
	struct INTMATRIX *pGeneLen;
	
	/* load source refgene */
	vSourceRefGene = NULL;
	vSourceRefGene = RefGene_LoadDatabase(strDatabaseFile, nDatabaseType, 
		strSpecies, &nSourceRefNum);
	if(vSourceRefGene == NULL)
	{
		printf("Warning: RNASEQ_CountReadPerTranscript_Main, null database!\n");
		return PROC_SUCCESS;
	}

	/* count sample number */
	if(nInputType == 0)
	{
		nSampleNum = 1;
	}
	else
	{
		nSampleNum = 0;
		fpIn = NULL;
		fpIn = fopen(strInFile, "r");
		if(fpIn == NULL)
		{
			printf("Error: RNASEQ_CountReadPerTranscript_Main, cannot open input file list!\n");
			exit(EXIT_FAILURE);
		}

		while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;
			if(strLine[0] == '#')
				continue;

			nSampleNum++;
		}

		fclose(fpIn);
	}

	/* prepare sample names */
	vSampleName = NULL;
	vSampleName = (struct tagString **)calloc(nSampleNum, sizeof(struct tagString *));
	if(vSampleName == NULL)
	{
		printf("Error: RNASEQ_CountReadPerTranscript_Main, cannot allocate memory for storing sample names!\n");
		exit(EXIT_FAILURE);
	}
	if(nInputType == 0)
	{
		StringAddTail(vSampleName+0, strInFile);
	}
	else
	{
		fpIn = NULL;
		fpIn = fopen(strInFile, "r");
		if(fpIn == NULL)
		{
			printf("Error: RNASEQ_CountReadPerTranscript_Main, cannot open input file list!\n");
			exit(EXIT_FAILURE);
		}

		ni = 0;
		while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;
			if(strLine[0] == '#')
				continue;

			StringAddTail(vSampleName+ni, strLine);
			ni++;
		}

		fclose(fpIn);

		if(ni != nSampleNum)
		{
			printf("Error: RNASEQ_CountReadPerTranscript_Main, sample number inconsistent!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* create expression matrix */
	vE = NULL;
	vE = (struct DOUBLEMATRIX **)calloc(nSampleNum, sizeof(struct DOUBLEMATRIX *));
	if(vE == NULL)
	{
		printf("Error: RNASEQ_CountReadPerTranscript_Main, cannot allocate memory for master expression vector!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nSampleNum; ni++)
	{
		vE[ni] = CreateDoubleMatrix(1, nSourceRefNum);
		if(vE[ni] == NULL)
		{
			printf("Error: RNASEQ_CountReadPerTranscript_Main, cannot allocate memory for storing expression values!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* process one by one */
	for(ni=0; ni<nSampleNum; ni++)
	{
		/* open file */
		fpIn = NULL;
		fpIn = fopen(vSampleName[ni]->m_pString, "r");
		if(fpIn == NULL)
		{
			printf("Error: RNASEQ_CountReadPerTranscript_Main, cannot open %s!\n", vSampleName[ni]->m_pString);
			exit(EXIT_FAILURE);
		}

		dT = 0.0;

		/* count */
		while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;
			if(strLine[0] == '#')
				continue;

			sscanf(strLine, "%s %d %s", strChr, &nPos, strStrand);
			nChr = Genome_ChromosomeName_To_Index(strChr, strSpecies);
			if(nChr < 0)
				continue;

			dT += 1;

			nx = 0;
			ny = nSourceRefNum-1;
			while(ny-nx > 1)
			{
				nz = (int)((nx+ny)/2);

				if(nChr < vSourceRefGene[nz]->nChrom)
				{
					ny = nz;
				}
				else if(nChr > vSourceRefGene[nz]->nChrom)
				{
					nx = nz;
				}
				else
				{
					if(nPos < vSourceRefGene[nz]->nTxStart)
						ny = nz;
					else if(nPos >= vSourceRefGene[nz]->nTxStart)
						nx = nz;
				}
			}

			
			if( (nChr == vSourceRefGene[ny]->nChrom) && (nPos < vSourceRefGene[ny]->nTxStart) )
			{
				nz = ny;
			}
			else if( (nChr == vSourceRefGene[nx]->nChrom) && (nPos < vSourceRefGene[nx]->nTxStart) )
			{
				nz = nx;
			}
			else if( (nChr < vSourceRefGene[ny]->nChrom) && (nChr == vSourceRefGene[nx]->nChrom) )
			{
				nz = nx;
			}
			else
			{
				continue;
			}

			for(nu=0; nu<nUp; nu++)
			{
				if(nChr != vSourceRefGene[nz]->nChrom)
					break;

				if((nPos >= vSourceRefGene[nz]->nTxStart) && (nPos <= vSourceRefGene[nz]->nTxEnd))
				{
					if( RNASEQ_HitTest(nChr, nPos, vSourceRefGene[nz]) == 1 )
					{
						vE[ni]->pMatElement[nz] += 1;
					}
				}

				nz--;
				if(nz < 0)
					break;
			}
		}

		/* close file */
		fclose(fpIn);

		if(nStandardize == 1)
		{
			if(dT > 0.5)
			{
				for(nz=0; nz<nSourceRefNum; nz++)
				{
					vE[ni]->pMatElement[nz] = vE[ni]->pMatElement[nz]*1e9/dT;
				}
			}
		}
	}

	/* standardize by length and read count */
	pGeneLen = NULL;
	pGeneLen = CreateIntMatrix(1, nSourceRefNum);
	if(pGeneLen == NULL)
	{
		printf("Error: RNASEQ_CountReadPerTranscript_Main, cannot allocate memory for storing gene length!\n");
		exit(EXIT_FAILURE);
	}
	for(nz=0; nz<nSourceRefNum; nz++)
	{
		nGeneLen = 0;
		for(ni=0; ni<vSourceRefGene[nz]->nExonCount; ni++)
		{
			nGeneLen += IMGETAT(vSourceRefGene[nz]->pmatExonStartsEnds, ni, 1)-IMGETAT(vSourceRefGene[nz]->pmatExonStartsEnds, ni, 0)+1;
		}

		if(nStandardize == 1)
		{
			if(nGeneLen > 0)
			{
				for(nj=0; nj<nSampleNum; nj++)
				{
					vE[nj]->pMatElement[nz] /= nGeneLen;
				}
			}
		}

		pGeneLen->pMatElement[nz] = nGeneLen;
	}

	/* export */
	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: RNASEQ_CountReadPerTranscript_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	fprintf(fpOut, "#Gene\tRefSeq\tEntrezID\tChr\tStart\tEnd\tStrand\tExonLen");
	for(nj=0; nj<nSampleNum; nj++)
	{
		GetFileName(vSampleName[nj]->m_pString, strLine);
		fprintf(fpOut, "\t%s", strLine);
	}
	fprintf(fpOut, "\n");
	for(ni=0; ni<nSourceRefNum; ni++)
	{
		fprintf(fpOut, "%s\t%s\t%d\t%s\t%d\t%d\t%c\t%d", vSourceRefGene[ni]->strGene, vSourceRefGene[ni]->strName, vSourceRefGene[ni]->nGeneID, 
			vSourceRefGene[ni]->strChrom, vSourceRefGene[ni]->nTxStart, vSourceRefGene[ni]->nTxEnd, vSourceRefGene[ni]->chStrand, pGeneLen->pMatElement[ni]);
		for(nj=0; nj<nSampleNum; nj++)
		{
			fprintf(fpOut, "\t%e", vE[nj]->pMatElement[ni]);
		}
		fprintf(fpOut, "\n");
	}
	fclose(fpOut);
	
	/* clear database */
	RefGene_ClearDatabase(&vSourceRefGene, nSourceRefNum);
	for(ni=0; ni<nSampleNum; ni++)
	{
		DestroyDoubleMatrix(vE[ni]);
		vE[ni] = NULL;
		DeleteString(vSampleName[ni]);
		vSampleName[ni] = NULL;
	}
	free(vE);
	free(vSampleName);
	DestroyIntMatrix(pGeneLen);
	
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RNASEQ_HitTest()                                                       */
/*  Test if a read is aligned to an exon of a gene.                        */
/* ----------------------------------------------------------------------- */ 
int RNASEQ_HitTest(int nChr, int nPos, struct tagRefGene *pRefGene)
{
	/* define */
	int nHit = 0;
	int ni;

	/* init check */
	if(pRefGene == NULL)
		return nHit;
	if(pRefGene->nChrom != nChr)
		return nHit;

	/* test */
	for(ni=0; ni<pRefGene->nExonCount; ni++)
	{
		if( (nPos >= IMGETAT(pRefGene->pmatExonStartsEnds, ni, 0)) && (nPos <= IMGETAT(pRefGene->pmatExonStartsEnds, ni, 1)) )
		{
			nHit = 1;
			break;
		}
	}

	/* return */
	return nHit;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_CountReads4RefGene_Main()                                          */
/*  Count reads for all transcripts.                                       */
/* ----------------------------------------------------------------------- */ 
int HTS_CountReads4RefGene_Main(char strDatabasePath[], int nDatabaseType,
			char strSpecies[], char strInputPath[], int nInputType, char strOutputPath[],
			int nRefType, int nUP, int nDOWN, int nStandardizebyT, int nStandardizebyL)
{
	/* out */
	struct tagRefGene **vSourceRefGene = NULL;
	int nSourceRefNum = 0;
	int nSampleNum = 0;
	FILE *fpIn;
	FILE *fpOut;
	char strLine[LONG_LINE_LENGTH];
	struct DOUBLEMATRIX **vE;
	int ni,nj,nk,nz;
	struct tagString **vSampleName = NULL;
	double dT;
	int nStart,nEnd;
	struct INTMATRIX *pGeneLen;
	struct tagBARData *pBARData;
	
	/* load source refgene */
	vSourceRefGene = NULL;
	vSourceRefGene = RefGene_LoadDatabase(strDatabasePath, nDatabaseType, 
		strSpecies, &nSourceRefNum);
	if(vSourceRefGene == NULL)
	{
		printf("Warning: HTS_CountReads4RefGene_Main, null database!\n");
		return PROC_SUCCESS;
	}
	
	/* compute gene length */
	pGeneLen = NULL;
	pGeneLen = CreateIntMatrix(1, nSourceRefNum);
	if(pGeneLen == NULL)
	{
		printf("Error: HTS_CountReads4RefGene_Main, cannot allocate memory for storing gene length!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nSourceRefNum; ni++)
	{
		if(nRefType == 0)
		{
			if(vSourceRefGene[ni]->chStrand == '-')
			{
				nStart = vSourceRefGene[ni]->nTxStart-nDOWN;
				nEnd = vSourceRefGene[ni]->nTxEnd+nUP;
			}
			else
			{
				nStart = vSourceRefGene[ni]->nTxStart-nUP;
				nEnd = vSourceRefGene[ni]->nTxEnd+nDOWN;
			}
		}
		else if(nRefType == 1)
		{
			if(vSourceRefGene[ni]->chStrand == '-')
			{
				nStart = vSourceRefGene[ni]->nTxEnd-nDOWN;
				nEnd = vSourceRefGene[ni]->nTxEnd+nUP;
			}
			else
			{
				nStart = vSourceRefGene[ni]->nTxStart-nUP;
				nEnd = vSourceRefGene[ni]->nTxStart+nDOWN;
			}
		}
		else if(nRefType == 2)
		{
			if(vSourceRefGene[ni]->chStrand == '-')
			{
				nStart = vSourceRefGene[ni]->nTxStart-nDOWN;
				nEnd = vSourceRefGene[ni]->nTxStart+nUP;
			}
			else
			{
				nStart = vSourceRefGene[ni]->nTxEnd-nUP;
				nEnd = vSourceRefGene[ni]->nTxEnd+nDOWN;
			}
		}
		else if(nRefType == 3)
		{
			if(vSourceRefGene[ni]->chStrand == '-')
			{
				nStart = vSourceRefGene[ni]->nCdsStart-nDOWN;
				nEnd = vSourceRefGene[ni]->nCdsEnd+nUP;
			}
			else
			{
				nStart = vSourceRefGene[ni]->nCdsStart-nUP;
				nEnd = vSourceRefGene[ni]->nCdsEnd+nDOWN;
			}
		}
		else if(nRefType == 4)
		{
			if(vSourceRefGene[ni]->chStrand == '-')
			{
				nStart = vSourceRefGene[ni]->nCdsEnd-nDOWN;
				nEnd = vSourceRefGene[ni]->nCdsEnd+nUP;
			}
			else
			{
				nStart = vSourceRefGene[ni]->nCdsStart-nUP;
				nEnd = vSourceRefGene[ni]->nCdsStart+nDOWN;
			}
		}
		else if(nRefType == 5)
		{
			if(vSourceRefGene[ni]->chStrand == '-')
			{
				nStart = vSourceRefGene[ni]->nCdsStart-nDOWN;
				nEnd = vSourceRefGene[ni]->nCdsStart+nUP;
			}
			else
			{
				nStart = vSourceRefGene[ni]->nCdsEnd-nUP;
				nEnd = vSourceRefGene[ni]->nCdsEnd+nDOWN;
			}
		}

		pGeneLen->pMatElement[ni] = nEnd-nStart+1;
		vSourceRefGene[ni]->nTxStart = nStart;
		vSourceRefGene[ni]->nTxEnd = nEnd;
	}


	/* count sample number */
	if(nInputType == 0)
	{
		nSampleNum = 1;
	}
	else
	{
		nSampleNum = 0;
		fpIn = NULL;
		fpIn = fopen(strInputPath, "r");
		if(fpIn == NULL)
		{
			printf("Error: HTS_CountReads4RefGene_Main, cannot open input file list!\n");
			exit(EXIT_FAILURE);
		}

		while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;
			if(strLine[0] == '#')
				continue;

			nSampleNum++;
		}

		fclose(fpIn);
	}

	/* prepare sample names */
	vSampleName = NULL;
	vSampleName = (struct tagString **)calloc(nSampleNum, sizeof(struct tagString *));
	if(vSampleName == NULL)
	{
		printf("Error: HTS_CountReads4RefGene_Main, cannot allocate memory for storing sample names!\n");
		exit(EXIT_FAILURE);
	}
	if(nInputType == 0)
	{
		StringAddTail(vSampleName+0, strInputPath);
	}
	else
	{
		fpIn = NULL;
		fpIn = fopen(strInputPath, "r");
		if(fpIn == NULL)
		{
			printf("Error: HTS_CountReads4RefGene_Main, cannot open input file list!\n");
			exit(EXIT_FAILURE);
		}

		ni = 0;
		while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;
			if(strLine[0] == '#')
				continue;

			StringAddTail(vSampleName+ni, strLine);
			ni++;
		}

		fclose(fpIn);

		if(ni != nSampleNum)
		{
			printf("Error: HTS_CountReads4RefGene_Main, sample number inconsistent!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* create expression matrix */
	vE = NULL;
	vE = (struct DOUBLEMATRIX **)calloc(nSampleNum, sizeof(struct DOUBLEMATRIX *));
	if(vE == NULL)
	{
		printf("Error: HTS_CountReads4RefGene_Main, cannot allocate memory for master count vector!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nSampleNum; ni++)
	{
		vE[ni] = CreateDoubleMatrix(1, nSourceRefNum);
		if(vE[ni] == NULL)
		{
			printf("Error: HTS_CountReads4RefGene_Main, cannot allocate memory for storing counts!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* process one sample by one sample */
	for(ni=0; ni<nSampleNum; ni++)
	{
		/* read data */
		pBARData = NULL;
		pBARData = Affy_LoadBAR_Fast(vSampleName[ni]->m_pString);
		if(pBARData == NULL)
		{
			printf("Warning: cannot open the BAR file %s\n", vSampleName[ni]->m_pString);
			continue;
		}

		/* count total reads */
		dT = 0.0;
		for(nj=0; nj<pBARData->nSeqNum; nj++)
		{
			for(nk=0; nk<pBARData->vSeqData[nj]->nDataNum; nk++)
			{
				if(pBARData->vSeqData[nj]->vData[1]->pMatElement[nk] < 0.0)
				{
					printf("Warning: negative read count!\n");
				}

				dT += pBARData->vSeqData[nj]->vData[1]->pMatElement[nk];
			}
		}


		/* process genes */
		for(nj=0; nj<nSourceRefNum; nj++)
		{
			nStart = vSourceRefGene[nj]->nTxStart;
			nEnd = vSourceRefGene[nj]->nTxEnd;
			vE[ni]->pMatElement[nj] = HTS_CountReads4Region(vSourceRefGene[nj]->strChrom, nStart, nEnd, pBARData);
		}

		/* standardize by total count */
		if(nStandardizebyT == 1)
		{
			if(dT > 0.5)
			{
				for(nz=0; nz<nSourceRefNum; nz++)
				{
					vE[ni]->pMatElement[nz] = vE[ni]->pMatElement[nz]*1e6/dT;
				}
			}
		}

		/* standardize by gene length */
		if(nStandardizebyL == 1)
		{
			for(nz=0; nz<nSourceRefNum; nz++)
			{
				if(pGeneLen->pMatElement[nz] > 0)
				{
					vE[ni]->pMatElement[nz] = vE[ni]->pMatElement[nz]*1e3/(double)(pGeneLen->pMatElement[nz]);
				}
				else
				{
					vE[ni]->pMatElement[nz] = 0.0;
				}
			}
		}

		Affy_BARData_Destroy(&pBARData);
	}

	/* export */
	fpOut = NULL;
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: HTS_CountReads4RefGene_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	fprintf(fpOut, "#Gene\tRefSeq\tEntrezID\tChr\tStart\tEnd\tStrand\tGeneLen");
	for(nj=0; nj<nSampleNum; nj++)
	{
		GetFileName(vSampleName[nj]->m_pString, strLine);
		fprintf(fpOut, "\t%s", strLine);
	}
	fprintf(fpOut, "\n");
	for(ni=0; ni<nSourceRefNum; ni++)
	{
		fprintf(fpOut, "%s\t%s\t%d\t%s\t%d\t%d\t%c\t%d", vSourceRefGene[ni]->strGene, vSourceRefGene[ni]->strName, vSourceRefGene[ni]->nGeneID, 
			vSourceRefGene[ni]->strChrom, vSourceRefGene[ni]->nTxStart, vSourceRefGene[ni]->nTxEnd, vSourceRefGene[ni]->chStrand, pGeneLen->pMatElement[ni]);
		for(nj=0; nj<nSampleNum; nj++)
		{
			fprintf(fpOut, "\t%e", vE[nj]->pMatElement[ni]);
		}
		fprintf(fpOut, "\n");
	}
	fclose(fpOut);
	
	/* clear database */
	RefGene_ClearDatabase(&vSourceRefGene, nSourceRefNum);
	for(ni=0; ni<nSampleNum; ni++)
	{
		DestroyDoubleMatrix(vE[ni]);
		vE[ni] = NULL;
		DeleteString(vSampleName[ni]);
		vSampleName[ni] = NULL;
	}
	free(vE);
	free(vSampleName);
	DestroyIntMatrix(pGeneLen);
	
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_CountReads4Region()                                                */
/*  Count reads for a region.                                              */
/* ----------------------------------------------------------------------- */ 
double HTS_CountReads4Region(char strChr[], int nStart, int nEnd, struct tagBARData *pBARData)
{
	/* define */
	double dT = 0.0;
	int ni,nx,ny,nz;
	int nP;

	/* init */
	nP = (int)((nStart+nEnd)/2);
	if(pBARData == NULL)
	{
		return dT;
	}

	/* find chromosome */
	for(ni=0; ni<pBARData->nSeqNum; ni++)
	{
		if(strcmp(strChr, pBARData->vSeqData[ni]->pSeqName->m_pString) != 0)
		{
			continue;
		}

		nx = 0;
		ny = pBARData->vSeqData[ni]->nDataNum-1;
		
		if(pBARData->vSeqData[ni]->vData[0]->pMatElement[nx] > nEnd)
		{
			continue;
		}

		if(pBARData->vSeqData[ni]->vData[0]->pMatElement[ny] < nStart)
		{
			continue;
		}


		while(ny-nx > 1)
		{
			nz = (int)((nx+ny)/2);

			if( (pBARData->vSeqData[ni]->vData[0]->pMatElement[nz] >= nStart) && (pBARData->vSeqData[ni]->vData[0]->pMatElement[nz] <= nEnd) )
			{
				break;
			}

			if(pBARData->vSeqData[ni]->vData[0]->pMatElement[nz] > nEnd)
			{
				ny = nz;
			}
			else
			{
				nx = nz;
			}
		}

		if( (pBARData->vSeqData[ni]->vData[0]->pMatElement[nz] < nStart) || (pBARData->vSeqData[ni]->vData[0]->pMatElement[nz] > nEnd) )
		{
			continue;
		}
		else
		{
			dT += pBARData->vSeqData[ni]->vData[1]->pMatElement[nz];

			nx = nz-1;
			while(nx >= 0)
			{
				if(pBARData->vSeqData[ni]->vData[0]->pMatElement[nx] < nStart)
					break;
				dT += pBARData->vSeqData[ni]->vData[1]->pMatElement[nx];
				nx--;
			}

			ny = nz+1;
			while(ny < pBARData->vSeqData[ni]->nDataNum)
			{
				if(pBARData->vSeqData[ni]->vData[0]->pMatElement[ny] > nEnd)
					break;
				dT += pBARData->vSeqData[ni]->vData[1]->pMatElement[ny];
				ny++;
			}
		}
	}

	/* return */
	return dT;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_Seg_Main()                                                    */
/*  Clustering analysis of ChIP-seq.                                       */
/* ----------------------------------------------------------------------- */ 
int SeqClust_Seg_Main(char strInputPath[], char strOutputPath[], char strOutputFile[],
				  int nBinSize, int nKernelType, int nKernelStep, 
				  int nKernelLen, int nKernelBand,
				  int nSegType, char strSegFile[], int nDistType, 
				  int nUp, int nDown, int nDatabaseType, 
				  int nMemBlock, int nCorrBlock, int nGridNum,
				  int nCutType, double dCutL, double dCutH,
				  int nExportBAR)
{
	/* define */
	int nResult;
	int ni;
	char strFileName[MED_LINE_LENGTH];

	/* For Step 1*/
	int nFileNum = 0;
	int nChrNum = 0;
	struct tagString **vChrName = NULL;
	struct tagString **vFileName = NULL;
	struct INTMATRIX *pChrLen = NULL;
	struct DOUBLEMATRIX *pFileReadCount = NULL;
	struct BYTEMATRIX *pIncludeChr = NULL;
	struct DOUBLEMATRIX *pFileNormFactor = NULL;
	struct DOUBLEMATRIX *pSmoothK = NULL;
	int nSmoothStart = 0;
	int nSmoothNum = 0;
	double dSmoothFactor = 0.0;
	int nBaseLineFileId = 0;
	
	/* For Step 3 */
	double dBinMax = 0.0;
	double dBinMin = 0.0;
	double dCL = 0.0;
	double dCH = 0.0;



	/* ---------------------------------------------------------*/
	/* STEP1: First scan, find chromosome number and length     */
	/*        Find read number for each sample                  */
	/* ---------------------------------------------------------*/
	printf("/* ------------------------------------------------- */\n");
	printf("/* STEP 1: Initialize                                */\n");
	printf("/* ------------------------------------------------- */\n");
	rand_u_init(13);
	AdjustDirectoryPath(strOutputPath);
	dSmoothFactor = 1.0;
	if(nKernelType > 0)
	{
		nResult = SeqClust_InitializeKernel(&pSmoothK, nKernelType, 
			nKernelStep, nKernelLen, nKernelBand, &nSmoothStart, &nSmoothNum);
		dSmoothFactor = 0.0;
		for(ni=0; ni<pSmoothK->nWidth; ni++)
		{
			dSmoothFactor += pSmoothK->pMatElement[ni];
		}
	}
	printf("Smooth Factor = %f\n", dSmoothFactor); 
	
	nResult = SeqClust_Initialize(strInputPath, &nFileNum, &vFileName, 
			&pFileReadCount, &pFileNormFactor, &nBaseLineFileId, &nChrNum, &vChrName, &pChrLen);
	for(ni=0; ni<pFileReadCount->nWidth; ni++)
	{
		pFileReadCount->pMatElement[ni] *= dSmoothFactor;
	}


	/* ---------------------------------------------------------*/
	/* STEP2: Count reads for genomic bins                      */
	/* ---------------------------------------------------------*/
	printf("\n");
	printf("/* ------------------------------------------------- */\n");
	printf("/* STEP 2: Count reads for genomic bins              */\n");
	printf("/* ------------------------------------------------- */\n");
	pIncludeChr = NULL;
	pIncludeChr = CreateByteMatrix(1, nChrNum);
	if(pIncludeChr == NULL)
	{
		printf("Error: SeqClust_Seg_Main, cannot create chromosome inclusion vector!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nChrNum; ni++)
		pIncludeChr->pMatElement[ni] = 1;
	nResult = SeqClust_CountBinReads(strInputPath, strOutputPath,
		strOutputFile, nFileNum, nChrNum, vChrName, pChrLen, nBinSize,
		nKernelType, pSmoothK, nSmoothStart, nKernelStep, nSmoothNum, pIncludeChr, nExportBAR);
	
	/* ---------------------------------------------------------*/
	/* STEP3 - STEP 5                                           */
	/* ---------------------------------------------------------*/
	if(nSegType != 0)
	{
		printf("\n");
		printf("/* ------------------------------------------------- */\n");
		printf("/* SKIP STEP 3 & 4. Use user supplied intervals.     */\n");
		printf("/* ------------------------------------------------- */\n");
	}
	/* automatic segmentation */
	else
	{
		/* ---------------------------------------------------------*/
		/* STEP3: Compute summary statistic for genomic bins        */
		/* ---------------------------------------------------------*/
		printf("\n");
		printf("/* ------------------------------------------------- */\n");
		printf("/* STEP 3: Compute summary statistics for bins       */\n");
		printf("/* ------------------------------------------------- */\n");
		nResult = SeqClust_ComputeBinStats(strOutputPath, strOutputFile, 
			nFileNum, vFileName, pFileNormFactor, nChrNum, vChrName, pChrLen,
			nBinSize, nMemBlock, nCorrBlock, &dBinMax, &dBinMin);

		/* cutoff by Percentile */
		if(nCutType == 0)
		{
			/* dBinMax = 34.0;
			dBinMin = -1.9; */
			SeqClust_QuantileCutoff(strOutputPath, strOutputFile, 
				nChrNum, vChrName, pChrLen, nBinSize, dBinMax, dBinMin, 
				nGridNum, dCutL, dCutH, &dCL, &dCH, pIncludeChr);
		}
		/* cutoff by FDR */
		else if(nCutType == 1)
		{
			SeqClust_FDRCutoff(strOutputPath, strOutputFile, 
				nChrNum, vChrName, pChrLen, nBinSize, dBinMax, dBinMin, 
				nGridNum, dCutL, dCutH, &dCL, &dCH, pIncludeChr);
		}
		/* cutoff supplied by users */
		else
		{
			dCL = dCutL;
			dCH = dCutH;
		}
		printf("Max. Bin Stat. = %f\n", dBinMax);
		printf("Min. Bin Stat. = %f\n", dBinMin);
		printf("Higher Cutoff = %f\n", dCH);
		printf("Lower Cutoff = %f\n", dCL);

		/* ---------------------------------------------------------*/
		/* STEP4: Genome segmentation                               */
		/* ---------------------------------------------------------*/
		printf("\n");
		printf("/* ------------------------------------------------- */\n");
		printf("/* STEP 4: Find correlation blocks                   */\n");
		printf("/* ------------------------------------------------- */\n");
		SeqClust_Segmentation(strOutputPath, strOutputFile, 
				nChrNum, vChrName, pChrLen, nBinSize, dCL, dCH);

		/* ---------------------------------------------------------*/
		/* STEP5: Collect information for genomic intervals         */
		/* ---------------------------------------------------------*/
		printf("\n");
		printf("/* ------------------------------------------------- */\n");
		printf("/* STEP 5: Collect data for clustering               */\n");
		printf("/* ------------------------------------------------- */\n");

		SeqClust_CollectData_RefineFileNormFactor(strOutputPath, strOutputFile, 
				nFileNum, vFileName, nBaseLineFileId,
				pFileReadCount, pFileNormFactor, dSmoothFactor,
				nChrNum, vChrName, pChrLen, nBinSize);

		SeqClust_CollectData_ForAutoSeg(strOutputPath, strOutputFile, 
				nFileNum, vFileName, pFileNormFactor, dSmoothFactor,
				nChrNum, vChrName, pChrLen, nBinSize, 
				pIncludeChr);
	}

	/* ---------------------------------------------------------*/
	/* STEP6: K-means clustering                                */
	/*        Perform from k_min to k_max. Choose the best      */
	/*        using BIC.                                        */
	/* ---------------------------------------------------------*/

	/* ---------------------------------------------------------*/
	/* STEP7: Release memory                                    */
	/* ---------------------------------------------------------*/
	DestroyDoubleMatrix(pSmoothK);
	DestroyIntMatrix(pChrLen);
	DestroyDoubleMatrix(pFileReadCount);
	DestroyDoubleMatrix(pFileNormFactor);
	DestroyByteMatrix(pIncludeChr);
	for(ni=0; ni<nChrNum; ni++)
	{
		DeleteString(vChrName[ni]);
		vChrName[ni] = NULL;
	}
	free(vChrName);
	for(ni=0; ni<nFileNum; ni++)
	{
		DeleteString(vFileName[ni]);
		vFileName[ni] = NULL;
	}
	free(vFileName);
	sprintf(strFileName, "%s%s_*.bincount", strOutputPath, strOutputFile);
	RemoveFiles(strFileName);
	sprintf(strFileName, "%s%s_*.tmpreg", strOutputPath, strOutputFile);
	RemoveFiles(strFileName);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_InitializeKernel()                                            */
/*  Clustering analysis of ChIP-seq: Initialize smoothing kernel vector    */ 
/*  Kernel Type                                                             */
/*    0: no kernel, use raw read count                                     */
/*    1(default): one-sided exponential                                    */
/*    2: two-sided exponential                                             */
/*    3: one-sided gaussian                                                */
/*    4: two-sided gaussian                                                */
/* ----------------------------------------------------------------------- */
int SeqClust_InitializeKernel(struct DOUBLEMATRIX **ppSmoothK, int nKernelType, 
			int nKernelStep, int nKernelLen, int nKernelBand, int *pKStart, int *pKNum)
{
	/* define */
	int nW,nH;
	int ni;
	double dx;
	double dSum; 
	double dMax;

	/* init */
	if((nKernelType > 4) || (nKernelType < 0))
	{
		printf("Error: SeqClust_InitializeKernel, unsupported kernel type!\n");
		exit(EXIT_FAILURE);
	}
	if(nKernelType == 0)
		return PROC_SUCCESS;

	if(nKernelBand < 0)
	{
		printf("Error: SeqClust_InitializeKernel, kernel bandwidth must be a positive integer!\n");
		exit(EXIT_FAILURE);
	}

	nH = nKernelLen/nKernelStep;
	if(nKernelType % 2 == 0)
	{
		nW = nH*2+1;
		*pKStart = -nH*nKernelStep;
		*pKNum = nW;
	}
	else
	{
		nW = nH+1;
		*pKStart = 0;
		*pKNum = nW;
	}

	*ppSmoothK = NULL;
	*ppSmoothK = CreateDoubleMatrix(1, nW);
	if(*ppSmoothK == NULL)
	{
		printf("Error: SeqClust_InitializeKernel, cannot create vector for smoothing kernel!\n");
		exit(EXIT_FAILURE);
	}

	/* compute two-sided kernel */
	dx = *pKStart;
	dSum = 0.0;
	dMax = 0.0;
	for(ni=0; ni<nW; ni++)
	{
		if( (nKernelType == 1) || (nKernelType == 2) )
		{
			(*ppSmoothK)->pMatElement[ni] = exp(-dx/nKernelBand)/nKernelBand;
		}
		else if( (nKernelType == 3) || (nKernelType == 4) )
		{
			(*ppSmoothK)->pMatElement[ni] = dnorm(dx, 0.0, nKernelBand);
		}

		if((*ppSmoothK)->pMatElement[ni] > dMax)
			dMax = (*ppSmoothK)->pMatElement[ni];

		dSum += (*ppSmoothK)->pMatElement[ni];
		dx += nKernelStep;
	}

	/* normalize */
	for(ni=0; ni<nW; ni++)
	{
		(*ppSmoothK)->pMatElement[ni] = (*ppSmoothK)->pMatElement[ni]/dMax;
	}
	
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_Initialize()                                                  */
/*  Clustering analysis of ChIP-seq: Initialize                            */
/* ----------------------------------------------------------------------- */ 
int SeqClust_Initialize(char strInputPath[], int *pFileNum, struct tagString ***vFileName,
						struct DOUBLEMATRIX **ppFileReadCount, 
						struct DOUBLEMATRIX **ppFileNormFactor, int *pBaseLineFileId, 
						int *pChrNum, 
						struct tagString ***vvChrName, struct INTMATRIX **ppChrLen)
{
	/* define */
	int nMaxChrNum = 65535;	
	int ni,nId,nk;
	int nChrNum = 0;
	int nFileNum = 0;
	int nReadCount;
	double dBaseCount;
	struct tagProbeGenomeInfo **vChrList = NULL;
	
	FILE *fpIn;
	FILE *fpRead;
	char strLine[LONG_LINE_LENGTH];
	struct DOUBLEMATRIX *pSampleReadCount;
	struct DOUBLEMATRIX *pSampleReadCountSort;
	struct LONGMATRIX *pSampleReadCountIdx;
	char strChr[LINE_LENGTH];
	char strFileName[MED_LINE_LENGTH];
	int nPos;
	int nStrand;

	/* Count file number */
	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: SeqClust_Initialize, cannot open input file list!\n");
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

		nFileNum++;
	}
	fclose(fpIn);
	printf("No. of Samples (Files) = %d\n", nFileNum); 

	/* prepare memory for file information */
	*pFileNum = nFileNum;

	*vFileName = NULL;
	*vFileName = (struct tagString **)calloc(nFileNum, sizeof(struct tagString *));
	if(*vFileName == NULL)
	{
		printf("Error: SeqClust_Initialize, cannot create vector for file names!\n");
		exit(EXIT_FAILURE);
	}

	*ppFileReadCount = NULL;
	*ppFileReadCount = CreateDoubleMatrix(1, nFileNum);
	if(*ppFileReadCount == NULL)
	{
		printf("Error: SeqClust_Initialize, cannot create vector for file read counts!\n");
		exit(EXIT_FAILURE);
	}

	pSampleReadCount = NULL;
	pSampleReadCount = CreateDoubleMatrix(1, nFileNum);
	if(pSampleReadCount == NULL)
	{
		printf("Error: SeqClust_Initialize, cannot create vector for file read counts!\n");
		exit(EXIT_FAILURE);
	}

	*ppFileNormFactor = NULL;
	*ppFileNormFactor = CreateDoubleMatrix(1, nFileNum);
	if(*ppFileNormFactor == NULL)
	{
		printf("Error: SeqClust_Initialize, cannot create vector for file normalizing factor!\n");
		exit(EXIT_FAILURE);
	}

	/* process files one by one */
	vChrList = NULL;
	vChrList = (struct tagProbeGenomeInfo **)calloc(nMaxChrNum, sizeof(struct tagProbeGenomeInfo *));
	if(vChrList == NULL)
	{
		printf("Error: SeqClust_Initialize, cannot create vector for chromosomes!\n");
		exit(EXIT_FAILURE);
	}
	nChrNum = 0;

	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: SeqClust_Initialize, cannot open input file list!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		/* read alignment file */
		GetFileName(strLine, strFileName);
		StringAddTail((*vFileName)+ni, strFileName);

		printf("  Processing %s ...\n", strLine);
		nReadCount = 0;
		fpRead = NULL;
		fpRead = fopen(strLine, "r");
		if(fpRead == NULL)
		{
			printf("Error: SeqClust_Initialize, cannot open read alignment file!\n");
			exit(EXIT_FAILURE);
		}
		while(fgets(strLine, LONG_LINE_LENGTH, fpRead) != NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;
			if(strLine[0] == '#')
				continue;

			SeqClust_ReadPos_From_Aln(strLine, strChr, &nPos, &nStrand);

			/* find chromosome and update chromosome length */
			nId = SeqClust_FindChrInfo(vChrList, nMaxChrNum, &nChrNum, strChr);
			if( (nId < 0) || (nId >= nMaxChrNum) )
			{
				printf("Error: SeqClust_Initialize, cannot find the matching chromosome!\n");
				exit(EXIT_FAILURE);
			}
			if(nPos > vChrList[nId]->nPos)
			{
				vChrList[nId]->nPos = nPos;
			}
			
			nReadCount++;
		}	

		/* close alignment file */
		fclose(fpRead);

		/* update file information */
		pSampleReadCount->pMatElement[ni] = nReadCount;
		(*ppFileReadCount)->pMatElement[ni] = nReadCount;

		ni++;
	}
	
	fclose(fpIn);
	
	if(ni != nFileNum)
	{
		printf("Error: SeqClust_Initialize, inconsistent file number!\n");
		exit(EXIT_FAILURE);
	}

	/* return chromosome name and length */
	*pChrNum = nChrNum;
	*vvChrName = (struct tagString **)calloc(nChrNum, sizeof(struct tagString *));
	if(*vvChrName == NULL)
	{
		printf("Error: SeqClust_Initialize, cannot create memory for chromosome name!\n");
		exit(EXIT_FAILURE);
	}

	*ppChrLen = NULL;
	*ppChrLen = CreateIntMatrix(1, nChrNum);
	if(*ppChrLen == NULL)
	{
		printf("Error: SeqClust_Initialize, cannot create memory for chromosome length!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nChrNum; ni++)
	{
		(*vvChrName)[ni] = vChrList[ni]->pProbe;
		vChrList[ni]->pProbe = NULL;
		(*ppChrLen)->pMatElement[ni] = vChrList[ni]->nPos + 1;

		ProbeGenomeInfoDestroy(vChrList+ni);
	}
	free(vChrList);

	if(ni != nChrNum)
	{
		printf("Error: SeqClust_Initialize, inconsistent chromosome number!\n");
		exit(EXIT_FAILURE);
	}

	/* compute normalizing factor */
	pSampleReadCountSort = NULL;
	pSampleReadCountIdx = NULL;
	DMSORTMERGEA_0(pSampleReadCount, &pSampleReadCountSort, &pSampleReadCountIdx);
	if(nFileNum % 2 == 0)
	{
		nk = nFileNum/2 - 1;
	}
	else
	{
		nk = nFileNum/2;
	}
	dBaseCount = pSampleReadCountSort->pMatElement[nk];
	*pBaseLineFileId = pSampleReadCountIdx->pMatElement[nk];
	for(ni=0; ni<nFileNum; ni++)
	{
		if( (*ppFileReadCount)->pMatElement[ni] <= 0)
		{
			printf("Error: SeqClust_Initialize, sample read count for the %d-th file is zero, cannot perform normalization!\n", ni+1);
			exit(EXIT_FAILURE);
		}
		(*ppFileNormFactor)->pMatElement[ni] = dBaseCount/(*ppFileReadCount)->pMatElement[ni];
	}


	/* print information */
	printf("\nFile\tNo_of_Reads\tScaling_Factor\n");
	for(ni=0; ni<(*pFileNum); ni++)
	{
		printf("%s\t%d\t%f\n", (*vFileName)[ni]->m_pString, (int)((*ppFileReadCount)->pMatElement[ni]), (*ppFileNormFactor)->pMatElement[ni]);
	}

	printf("\nNo. of chromosomes = %d\n", *pChrNum);
	printf("\tChromosome\tMax_Coordinate\n");
	for(ni=0; ni<(*pChrNum); ni++)
	{
		printf("\t%s\t%d\n", (*vvChrName)[ni]->m_pString, (*ppChrLen)->pMatElement[ni]);
	}

	/* release memory */
	DestroyDoubleMatrix(pSampleReadCount);
	DestroyDoubleMatrix(pSampleReadCountSort);
	DestroyLongMatrix(pSampleReadCountIdx);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_ReadPos_From_Aln()                                            */
/*  Clustering analysis of ChIP-seq: Parse coordinates from aln file       */
/* ----------------------------------------------------------------------- */ 
int SeqClust_ReadPos_From_Aln(char strLine[], char strChr[], int *pPos, int *pStrand)
{
	/* define */
	char chStrand;

	sscanf(strLine, "%s %d %c", strChr, pPos, &chStrand);
	if((chStrand == '-') || (chStrand == 'r') || (chStrand == 'R'))
		*pStrand = 1;
	else
		*pStrand = 0;

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_ReadPos_From_Cod()                                            */
/*  Clustering analysis of ChIP-seq: Parse coordinates from cod file       */
/* ----------------------------------------------------------------------- */ 
int SeqClust_ReadPos_From_Cod(char strLine[], char strAlias[], char strChr[], 
							  int *pPos1, int *pPos2, int *pStrand)
{
	/* define */
	char chStrand;

	sscanf(strLine, "%s %s %d %d %c", strAlias, strChr, pPos1, pPos2, &chStrand);
	if((chStrand == '-') || (chStrand == 'r') || (chStrand == 'R'))
		*pStrand = 1;
	else
		*pStrand = 0;

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_FindChr()                                                     */
/*  Clustering analysis of ChIP-seq: Find matching chromosome              */
/*  If no matching chromosome is found, the chromosome will be added to    */
/*  the chromosome list. If there are more than nMaxChrNum, return -1.     */
/*  otherwise return the index of the matching chromosome.                 */
/* ----------------------------------------------------------------------- */ 
int SeqClust_FindChrInfo(struct tagProbeGenomeInfo **vChrList, int nMaxChrNum, 
						 int *pChrNum, char strChr[])
{
	/* define */
	int nId  = -1;
	int ni,nx,ny,nz;
	struct tagProbeGenomeInfo *pChrNode;
	int ncomp;

	/* add new chromosome to empty list */
	if(*pChrNum == 0)
	{
		if(*pChrNum == nMaxChrNum)
		{
			printf("Error: SeqClust_FindChrInfo, chromosome number exceeds %d!\n", nMaxChrNum);
			nId = -1;
		}
		else
		{
			nId = 0;
			pChrNode = NULL;
			pChrNode = ProbeGenomeInfoCreate();
			if(pChrNode == NULL)
			{
				printf("Error: SeqClust_FindChrInfo, cannot create new chromosome node!\n");
				exit(EXIT_FAILURE);
			}
			pChrNode->nPos = 0;
			StringAddTail(&(pChrNode->pProbe), strChr);
			vChrList[nId] = pChrNode;
			*pChrNum = (*pChrNum) + 1;
		}
	}
	/* find chromosome */
	else
	{
		nx = 0;
		ny = (*pChrNum)-1;
		
		/* if smaller than the first chromosome */
		if(strcmp(strChr, vChrList[nx]->pProbe->m_pString) < 0)
		{
			if(*pChrNum == nMaxChrNum)
			{
				printf("Error: SeqClust_FindChrInfo, chromosome number exceeds %d!\n", nMaxChrNum);
				nId = -1;
			}
			else
			{
				for(ni=(*pChrNum); ni>0; ni--)
				{
					vChrList[ni] = vChrList[ni-1];
				}

				nId = 0;
				pChrNode = NULL;
				pChrNode = ProbeGenomeInfoCreate();
				if(pChrNode == NULL)
				{
					printf("Error: SeqClust_FindChrInfo, cannot create new chromosome node!\n");
					exit(EXIT_FAILURE);
				}
				pChrNode->nPos = 0;
				StringAddTail(&(pChrNode->pProbe), strChr);
				vChrList[nId] = pChrNode;
				*pChrNum = (*pChrNum) + 1;
			}

			return nId;
		}
		
		/* if bigger than the last chromosome */
		if(strcmp(strChr, vChrList[ny]->pProbe->m_pString) > 0)
		{
			if(*pChrNum == nMaxChrNum)
			{
				printf("Error: SeqClust_FindChrInfo, chromosome number exceeds %d!\n", nMaxChrNum);
				nId = -1;
			}
			else
			{
				nId = *pChrNum;
				pChrNode = NULL;
				pChrNode = ProbeGenomeInfoCreate();
				if(pChrNode == NULL)
				{
					printf("Error: SeqClust_FindChrInfo, cannot create new chromosome node!\n");
					exit(EXIT_FAILURE);
				}
				pChrNode->nPos = 0;
				StringAddTail(&(pChrNode->pProbe), strChr);
				vChrList[nId] = pChrNode;
				*pChrNum = (*pChrNum) + 1;
			}

			return nId;
		}

		/* if in between */
		nId = -1;
		while(ny-nx > 1)
		{
			nz = (ny+nx)/2;
			ncomp = strcmp(strChr, vChrList[nz]->pProbe->m_pString);
			if(ncomp == 0)
			{
				nId = nz;
				break;
			}
			else if(ncomp < 0)
			{
				ny = nz;
			}
			else
			{
				nx = nz;
			}
		}

		if(nId < 0)
		{
			if(strcmp(strChr, vChrList[nx]->pProbe->m_pString) == 0)
			{
				nId = nx;
			}
			else if(strcmp(strChr, vChrList[ny]->pProbe->m_pString) == 0)
			{
				nId = ny;
			}
			else
			{
				if(*pChrNum == nMaxChrNum)
				{
					printf("Error: SeqClust_FindChrInfo, chromosome number exceeds %d!\n", nMaxChrNum);
					nId = -1;
				}
				else
				{
					for(ni=(*pChrNum); ni>ny; ni--)
					{
						vChrList[ni] = vChrList[ni-1];
					}

					nId = ny;
					pChrNode = NULL;
					pChrNode = ProbeGenomeInfoCreate();
					if(pChrNode == NULL)
					{
						printf("Error: SeqClust_FindChrInfo, cannot create new chromosome node!\n");
						exit(EXIT_FAILURE);
					}
					pChrNode->nPos = 0;
					StringAddTail(&(pChrNode->pProbe), strChr);
					vChrList[nId] = pChrNode;
					*pChrNum = (*pChrNum) + 1;
				}
			}
		}
	}

	/* return */
	return nId;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_CountBinReads()                                               */
/*  Count reads for genomic bins                                           */
/* ----------------------------------------------------------------------- */ 
int SeqClust_CountBinReads(char strInputPath[], char strOutputPath[],
		char strOutputFile[], int nFileNum, int nChrNum, 
		struct tagString **vChrName, struct INTMATRIX *pChrLen, int nBinSize,
		int nKernelType, struct DOUBLEMATRIX *pSmoothK, int nSmoothStart, 
		int nKernelStep, int nSmoothNum, struct BYTEMATRIX *pIncludeChr,
		int nExportBAR)
{
	/* define */
	int ni;
	int nResult;
	FILE *fpIn;
	char strLine[LONG_LINE_LENGTH];

	/* open file */
	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: SeqClust_CountBinReads, cannot open input file list!\n");
		exit(EXIT_FAILURE);
	}

	/* process one by one */
	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		/* read alignment file */
		printf("  Processing %s ...\n", strLine);
		nResult = SeqClust_CountBinReads_SingleFile(strLine, strOutputPath,
				strOutputFile, nChrNum, vChrName, pChrLen, nBinSize,
				nKernelType, pSmoothK, nSmoothStart, nKernelStep, nSmoothNum,
				pIncludeChr, nExportBAR);		
	
		ni++;
	}
	
	/* close file */
	fclose(fpIn);
	
	if(ni != nFileNum)
	{
		printf("Error: SeqClust_CountBinReads, inconsistent file number!\n");
		exit(EXIT_FAILURE);
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_CountBinReads_SingleFile()                                    */
/*  Count reads for genomic bins for a single sample file                  */
/* ----------------------------------------------------------------------- */ 
int SeqClust_CountBinReads_SingleFile(char strInputFile[], char strOutputPath[],
				char strOutputFile[], int nChrNum, struct tagString **vChrName, 
				struct INTMATRIX *pChrLen, int nBinSize, 		
				int nKernelType, struct DOUBLEMATRIX *pSmoothK, int nSmoothStart, 
				int nKernelStep, int nSmoothNum, struct BYTEMATRIX *pIncludeChr,
				int nExportBAR)
{
	/* define */
	float **vBinC = NULL;
	int ni,nj,nx;
	int *vBinNum;
	int nBinNum;
	char strFileName[LINE_LENGTH];
	char strOutFileName[MED_LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	FILE *fpRead;
	char strChr[LINE_LENGTH];
	int nPos;
	int nStrand;
	int nChrId,nBinId;
	char strCGWFileName[MED_LINE_LENGTH];

	/* init */
	GetFileName(strInputFile, strFileName);
	vBinC = NULL;
	vBinC = (float **)calloc(nChrNum, sizeof(float *));
	if(vBinC == NULL)
	{
		printf("Error: SeqClust_CountBinReads_SingleFile, cannot create vector for genomic bins!\n");
		exit(EXIT_FAILURE);
	}

	vBinNum = NULL;
	vBinNum = (int *)calloc(nChrNum, sizeof(int));
	if(vBinNum == NULL)
	{
		printf("Error: SeqClust_CountBinReads_SingleFile, cannot create vector for genomic bin size!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nChrNum; ni++)
	{
		nBinNum = pChrLen->pMatElement[ni]/nBinSize;
		if(pChrLen->pMatElement[ni] % nBinSize != 0)
			nBinNum++;

		vBinC[ni] = (float *)calloc(nBinNum, sizeof(float));
		if(vBinC[ni] == NULL)
		{
			printf("Error: SeqClust_CountBinReads_SingleFile, insufficient memory for creating genomic bins, try a larger bin size!\n");
			exit(EXIT_FAILURE);
		}		

		vBinNum[ni] = nBinNum;
	}

	/* count */
	fpRead = NULL;
	fpRead = fopen(strInputFile, "r");
	if(fpRead == NULL)
	{
		printf("Error: SeqClust_CountBinReads_SingleFile, cannot open read alignment file!\n");
		exit(EXIT_FAILURE);
	}
	while(fgets(strLine, LONG_LINE_LENGTH, fpRead) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		SeqClust_ReadPos_From_Aln(strLine, strChr, &nPos, &nStrand);

		/* find chromosome and update chromosome length */
		nChrId = SeqClust_FindChr(vChrName, nChrNum, strChr);
		if( (nChrId < 0) || (nChrId >= nChrNum) )
		{
			printf("Warning: SeqClust_CountBinReads_SingleFile, cannot find the matching chromosome for %s:%d!\n",strChr, nPos);
			continue;
		}

		/* update bin count */
		if(nKernelType == 0)
		{
			nBinId = nPos/nBinSize;
			if(nBinId >= vBinNum[nChrId])
			{
				printf("Error: SeqClust_CountBinReads_SingleFile, inconsistent bin number!\n");
				exit(EXIT_FAILURE);
			}
			vBinC[nChrId][nBinId] += 1;
		}
		else
		{
			nBinId = nPos/nBinSize;
			if(nBinId >= vBinNum[nChrId])
			{
				printf("Error: SeqClust_CountBinReads_SingleFile, inconsistent bin number!\n");
				exit(EXIT_FAILURE);
			}

			/* '-' strand */
			if(nStrand == 1)
			{
				nx = nPos-nSmoothStart;
				for(ni=0; ni<nSmoothNum; ni++)
				{
					nBinId = nx/nBinSize;
					if( (nBinId >= 0) && (nBinId < vBinNum[nChrId]) )
					{
						vBinC[nChrId][nBinId] += (float)(pSmoothK->pMatElement[ni]);
					}
					nx -= nKernelStep;
				}
			}
			/* '+' strand */
			else
			{
				nx = nPos+nSmoothStart;
				for(ni=0; ni<nSmoothNum; ni++)
				{
					nBinId = nx/nBinSize;
					if( (nBinId >= 0) && (nBinId < vBinNum[nChrId]) )
					{
						vBinC[nChrId][nBinId] += (float)(pSmoothK->pMatElement[ni]);
					}
					nx += nKernelStep;
				}
			}
		}
	}	

	/* close alignment file */
	fclose(fpRead);

	/* create chromosome indicator */
	for(ni=0; ni<nChrNum; ni++)
	{
		nx = 0;
		for(nj=0; nj<vBinNum[ni]; nj++)
		{
			if(vBinC[ni][nj] > 1e-6)
			{
				nx = 1;
				break;
			}
		}

		if(nx == 0)
		{
			pIncludeChr->pMatElement[ni] = 0;
		}
	}


	/* save & release memory */
	if(nExportBAR == 1)
	{
		sprintf(strOutFileName, "%s%s_%s.bar.txt", strOutputPath, strOutputFile, strFileName);
		fpRead = NULL;
		fpRead = fopen(strOutFileName, "w");
		if(fpRead == NULL)
		{
			printf("Error: cannot open output file!\n");
			exit(EXIT_FAILURE);
		}

		fprintf(fpRead, "#chr\tpos\t%s_binc\n", strFileName);
		fprintf(fpRead, "#chr\tpos\t1\n");
		for(ni=0; ni<nChrNum; ni++)
		{
			nPos = 12;
			for(nx=0; nx<vBinNum[ni]; nx++)
			{
				fprintf(fpRead, "%s\t%d\t%f\n", vChrName[ni]->m_pString, nPos, vBinC[ni][nx]);
				nPos += nBinSize;
			}
		}

		fclose(fpRead);

		sprintf(strCGWFileName, "%s.cgw", strOutputFile);
		TileMapv2_TXT2BAR(strOutFileName, strOutputPath, strCGWFileName);
	}
	
	for(ni=0; ni<nChrNum; ni++)
	{
		sprintf(strOutFileName, "%s%s_%s_%s.bincount", strOutputPath, strOutputFile, strFileName, vChrName[ni]->m_pString);
		TileMapv2_SaveToBinaryFile((void *)(vBinC[ni]), sizeof(float), vBinNum[ni], strOutFileName);

		free(vBinC[ni]);
		vBinC[ni] = NULL;
	}
	free(vBinC);
	free(vBinNum);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_FindChr()                                                     */
/*  Find matching chromosome                                               */
/*  If no matching chromosome is found, return -1.                         */
/*  otherwise return the index of the matching chromosome.                 */
/* ----------------------------------------------------------------------- */ 
int SeqClust_FindChr(struct tagString **vChrName, int nChrNum, char strChr[])
{
	/* define */
	int nId = -1;
	int nx,ny,nz;
	int ncomp;

	/* init */
	if(nChrNum <= 0)
		return nId;

	nx = 0;
	ny = nChrNum-1;
	
	/* if smaller than the first chromosome */
	ncomp = strcmp(strChr, vChrName[nx]->m_pString);
	if(ncomp < 0)
	{
		return nId;
	}
	else if(ncomp == 0)
	{
		return nx;
	}
	
	/* if bigger than the last chromosome */
	ncomp = strcmp(strChr, vChrName[ny]->m_pString);
	if(ncomp > 0)
	{
		return nId;
	}
	else if(ncomp == 0)
	{
		return ny;
	}

	/* if in between */
	nId = -1;
	while(ny-nx > 1)
	{
		nz = (ny+nx)/2;
		ncomp = strcmp(strChr, vChrName[nz]->m_pString);
		if(ncomp == 0)
		{
			nId = nz;
			break;
		}
		else if(ncomp < 0)
		{
			ny = nz;
		}
		else
		{
			nx = nz;
		}
	}

	if(nId < 0)
	{
		if(strcmp(strChr, vChrName[nx]->m_pString) == 0)
		{
			nId = nx;
		}
		else if(strcmp(strChr, vChrName[ny]->m_pString) == 0)
		{
			nId = ny;
		}
	}


	/* return */
	return nId;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_ComputeBinStats()                                             */
/*  Compute summary statistics for genomic bins.                           */
/* ----------------------------------------------------------------------- */ 
int SeqClust_ComputeBinStats(char strOutputPath[], char strOutputFile[], 
			int nFileNum, struct tagString **vFileName, struct DOUBLEMATRIX *pFileNormFactor, 
			int nChrNum, struct tagString **vChrName, struct INTMATRIX *pChrLen,
			int nBinSize, int nMemBlock, int nCorrBlock, double *pBinMax, double *pBinMin)
{
	/* define */
	int ni;
	int nBinNum;

	/* process chromosome by chromosome */
	for(ni=0; ni<nChrNum; ni++)
	{
		nBinNum = pChrLen->pMatElement[ni]/nBinSize;
		if(pChrLen->pMatElement[ni] % nBinSize != 0)
			nBinNum++;

		SeqClust_ComputeBinStats_Chr(strOutputPath, strOutputFile, 
			nFileNum, vFileName, pFileNormFactor, vChrName[ni]->m_pString, 
			nBinNum, nBinSize, nMemBlock, nCorrBlock, pBinMax, pBinMin);
	}


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_ComputeBinStats_Chr()                                         */
/*  Compute summary statistics for genomic bins for a single chromosome    */
/* ----------------------------------------------------------------------- */ 
int SeqClust_ComputeBinStats_Chr(char strOutputPath[], char strOutputFile[], 
			int nFileNum, struct tagString **vFileName, struct DOUBLEMATRIX *pFileNormFactor, 
			char strChr[], int nBinNum, int nBinSize, int nMemBlock, int nCorrBlock,
			double *pBinMax, double *pBinMin)
{
	/* define */
	float *vV = NULL;
	float *vM = NULL;
	float *vC = NULL;
	float *vC0 = NULL;
	float *vCorr = NULL;
	float *vCov = NULL;
	float *vMAV = NULL;
	float *vMAM = NULL;
	char strOutFileName[MED_LINE_LENGTH];
	char strInFileName[MED_LINE_LENGTH];
	float **vD = NULL;
	/* struct DOUBLEMATRIX *pLogNormFactor; */
	double dLog2 = log(2.0);

	FILE **vfpIn;
	int ni,nj,nk,nw1,nw2,nu,nz;
	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));

	int nP0,nP1;
	int nQ0,nQ1;
	int nCount,nQCount;
	double dMu,dVar,dCorr,dTemp,dCov,dC,dC0,dMAV,dMAM,dRand;
	float *pV,*pC,*pC0,*pM,*pCorr,*pCov,*pMAV,*pMAM;
	int nRandi;

	FILE *fpRead;
	char strLine[MED_LINE_LENGTH];
	char strTileMapFolder[MED_LINE_LENGTH];

	/* init */
	/* pLogNormFactor = NULL;
	pLogNormFactor = CreateDoubleMatrix(1, nFileNum);
	if(pLogNormFactor == NULL)
	{
		printf("Error: SeqClust_ComputeBinStats_Chr, cannot allocate memory for log norm factor!\n");
		exit(EXIT_FAILURE);
	} */

	vV = (float *)calloc(nBinNum, sizeof(float));
	if(vV == NULL)
	{
		printf("Error: SeqClust_ComputeBinStats_Chr, cannot allocate memory for variance vector!\n");
		exit(EXIT_FAILURE);
	}
	vC = (float *)calloc(nBinNum, sizeof(float));
	if(vC == NULL)
	{
		printf("Error: SeqClust_ComputeBinStats_Chr, cannot allocate memory for correlation vector!\n");
		exit(EXIT_FAILURE);
	}
	vC0 = (float *)calloc(nBinNum, sizeof(float));
	if(vC0 == NULL)
	{
		printf("Error: SeqClust_ComputeBinStats_Chr, cannot allocate memory for correlation vector!\n");
		exit(EXIT_FAILURE);
	} 
	vM = (float *)calloc(nBinNum, sizeof(float));
	if(vM == NULL)
	{
		printf("Error: SeqClust_ComputeBinStats_Chr, cannot allocate memory for read count vector!\n");
		exit(EXIT_FAILURE);
	}
	vCorr = (float *)calloc(nBinNum, sizeof(float));
	if(vCorr == NULL)
	{
		printf("Error: SeqClust_ComputeBinStats_Chr, cannot allocate memory for read count vector!\n");
		exit(EXIT_FAILURE);
	}
	vCov = (float *)calloc(nBinNum, sizeof(float));
	if(vCov == NULL)
	{
		printf("Error: SeqClust_ComputeBinStats_Chr, cannot allocate memory for read count vector!\n");
		exit(EXIT_FAILURE);
	}
	vMAV = (float *)calloc(nBinNum, sizeof(float));
	if(vMAV == NULL)
	{
		printf("Error: SeqClust_ComputeBinStats_Chr, cannot allocate memory for read count vector!\n");
		exit(EXIT_FAILURE);
	}
	vMAM = (float *)calloc(nBinNum, sizeof(float));
	if(vMAM == NULL)
	{
		printf("Error: SeqClust_ComputeBinStats_Chr, cannot allocate memory for read count vector!\n");
		exit(EXIT_FAILURE);
	}

	vfpIn = NULL;
	vfpIn = (FILE **)calloc(nFileNum, sizeof(FILE *));
	if(vfpIn == NULL)
	{
		printf("Error: SeqClust_ComputeBinStats_Chr, cannot create vector for file pointers!\n");
		exit(EXIT_FAILURE);
	}

	
	for(ni=0; ni<nFileNum; ni++)
	{
		sprintf(strInFileName, "%s%s_%s_%s.bincount", strOutputPath, strOutputFile, vFileName[ni]->m_pString, strChr);
		vfpIn[ni] = fopen(strInFileName, "rb");
		if(vfpIn[ni] == NULL)
		{
			printf("Error: SeqClust_ComputeBinStats_Chr, cannot open file %s!\n", strInFileName);
			exit(EXIT_FAILURE);
		}
	}

	vD = NULL;
	vD = (float **)calloc(nFileNum, sizeof(float *));
	if(vD == NULL)
	{
		printf("Error: SeqClust_ComputeBinStats_Chr, cannot create vector for summary statistics!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nFileNum; ni++)
	{
		vD[ni] = (float *)calloc((nMemBlock+nCorrBlock+nCorrBlock), sizeof(float));
		if(vD[ni] == NULL)
		{
			printf("Error: SeqClust_ComputeBinStats_Chr, cannot create memory block for computing summary statistics!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* processing */
	nP0 = 0;
	while(nP0 < nBinNum)
	{
		/* determine memory block */
		nP1 = nP0 + nMemBlock-1;
		if(nP1 >= nBinNum)
			nP1 = nBinNum-1;
		nCount = nP1-nP0+1;
		pV = vV+nP0;
		pC = vC+nP0;
		pC0 = vC0+nP0;
		pM = vM+nP0;
		pCorr = vCorr+nP0;
		pCov = vCov+nP0;
		pMAV = vMAV+nP0;
		pMAM = vMAM+nP0;

		nQ0 = nP0-nCorrBlock;
		nQ1 = nP1+nCorrBlock;
		if(nQ0 < 0)
			nQ0 = 0;
		if(nQ1 >= nBinNum)
			nQ1 = nBinNum-1;
		nQCount = nQ1-nQ0+1;

		/* read data */
		for(ni=0; ni<nFileNum; ni++)
		{
			fseek(vfpIn[ni], nQ0*sizeof(float), SEEK_SET);
			if( little_endian_fread(vD[ni], sizeof(float), nQCount, vfpIn[ni], little_endian_machine) != nQCount)
			{
				printf("Error: SeqClust_ComputeBinStats_Chr, incorrect loading, number of bins inconsistent!\n");
				exit(EXIT_FAILURE);
			}
		}

		/* process */
		nk = nQ0-nP0;
		for(nj=0; nj<nQCount; nj++)
		{
			dMu = 0.0;
			dVar = 0.0;
			for(ni=0; ni<nFileNum; ni++)
			{
				dTemp = vD[ni][nj]*(pFileNormFactor->pMatElement[ni]);
				if(dTemp < 1.0)
					dTemp = 1.0;
				vD[ni][nj] = (float)(log(dTemp)/dLog2);
				dMu += vD[ni][nj];
			}

			dMu /= nFileNum;

			for(ni=0; ni<nFileNum; ni++)
			{
				dTemp = vD[ni][nj]-dMu;
				dVar += (dTemp*dTemp);
				vD[ni][nj] = (float)dTemp;
			}

			if(nFileNum > 1)
				dVar /= (nFileNum-1);
			else
				dVar = 0.0;

			dVar = sqrt(dVar);
			
			if(nk < 0)
			{
				if( (fabs(pM[nk] - dMu) > 1e-6) || (fabs(pV[nk] - dVar) > 1e-6) )
				{
					printf("Error: SeqClust_ComputeBinStats_Chr, inconsistent computation results!\n");
					exit(EXIT_FAILURE);
				}
			}

			*(pM+nk) = (float)dMu;
			*(pV+nk) = (float)dVar;

			nk++;
		}

		nk = nQ0-nP0;
		for(nj=0; nj< nQCount; nj++)
		{
			if( (nk >= 0) && (nk < nCount) )
			{
				dC = 0.0;
				dC0 = 0.0;
				dCorr = 0.0;
				dCov = 0.0;
				dMAV = 0.0;
				dMAM = 0.0;
				nz = 0;
				dRand = rand_u();
				nRandi = (int)(nCount*dRand);
				while( ((nRandi-nCorrBlock) < 0) || ( (nRandi+nCorrBlock) >= nCount) )
				{
					dRand = rand_u();
					nRandi = (int)(nQCount*dRand);
					nz++;
					if(nz >= 10000)
					{
						printf("Warning: difficult to generate random permutation!\n");
						nRandi = nj;
						break;
					}
				}

				for(nu=1; nu<=nCorrBlock; nu++)
				{
					nw1 = nj-nu;
					nw2 = nj+nu;
					if( (nw1 < 0) || (nw2 >=nQCount) )
					{
						continue;
					}

					dTemp = 0.0;
					for(ni=0; ni<nFileNum; ni++)
					{
						dTemp += vD[ni][nw1]*vD[ni][nw2];
					}

					if(nFileNum > 1)
					{
						dTemp = dTemp/(nFileNum-1);
						dCov += dTemp;
						dTemp = dTemp/(*(pV+nk-nu)+1e-6)/(*(pV+nk+nu)+1e-6);
						dCorr += dTemp;
						dC0 += dTemp*(*(pV+nRandi-nu)+1e-6)*(*(pV+nRandi+nu)+1e-6);
					}

					dMAM += (*(pM+nk-nu)) + (*(pM+nk+nu));
					dMAV += (*(pV+nk-nu))*(*(pV+nk-nu)) + (*(pV+nk+nu))*(*(pV+nk+nu));
				}
				
				dTemp = *(pV+nRandi);
				dTemp = dTemp*dTemp;
				dC0 += dTemp;

				dTemp = *(pV+nk);
				dTemp = dTemp*dTemp;
				dC = dCov+dTemp;

				dMAV += dTemp;
				dMAM += (*(pM+nk));

				dC /= (nCorrBlock+1);
				dC0 /= (nCorrBlock+1);
				dCorr /= nCorrBlock;
				dCov /= nCorrBlock;
				dMAV /= (nCorrBlock+nCorrBlock+1);
				dMAM /= (nCorrBlock+nCorrBlock+1);

				*(pC+nk) = (float)dC;
				if(dC > (*pBinMax))
					*pBinMax = dC;
				if(dC < (*pBinMin))
					*pBinMin = dC;
				*(pC0+nk) = (float)dC0;
				*(pCorr+nk) = (float)dCorr;
				*(pCov+nk) = (float)dCov;
				*(pMAV+nk) = (float)dMAV;
				*(pMAM+nk) = (float)dMAM;
			}

			nk++;
		}

		/* get next memory block */
		nP0 = nP1+1;
	}

	
	/* close files and free memory blocks */
	for(ni=0; ni<nFileNum; ni++)
	{
		fclose(vfpIn[ni]);
		vfpIn[ni] = NULL;
		free(vD[ni]);
	}
	free(vfpIn);
	free(vD);

	/* save data to files */
	sprintf(strOutFileName, "%s%s_%s.cma", strOutputPath, strOutputFile, strChr);
	TileMapv2_SaveToBinaryFile((void *)vC, sizeof(float), nBinNum, strOutFileName);
	sprintf(strOutFileName, "%s%s_%s.cmo", strOutputPath, strOutputFile, strChr);
	TileMapv2_SaveToBinaryFile((void *)vC0, sizeof(float), nBinNum, strOutFileName);
	/* sprintf(strOutFileName, "%s%s_%s.ave", strOutputPath, strOutputFile, strChr);
	TileMapv2_SaveToBinaryFile((void *)vM, sizeof(float), nBinNum, strOutFileName); */

	/* FOR DEBUG & EXPLORE */
	/* if(strcmp(strChr, "chr10") == 0) */
	if(1 == 0)
	{
		strcpy(strTileMapFolder, "/home/bst/faculty/hji/projects/cisgenome_project/bin/");
		/* strcpy(strTileMapFolder, "D:\\Projects\\cisgenome_project\\bin\\"); */

		sprintf(strOutFileName, "%s%s_%s.var.bar.txt", strOutputPath, strOutputFile, strChr);
		fpRead = NULL;
		fpRead = fopen(strOutFileName, "w");
		if(fpRead == NULL)
		{
			printf("Error: cannot open output file!\n");
			exit(EXIT_FAILURE);
		}

		fprintf(fpRead, "#chr\tpos\t%s_var_%s\n", strOutputFile,strChr);
		fprintf(fpRead, "#chr\tpos\t1\n");
		nP0 = 12;
		for(ni=0; ni<nBinNum; ni++)
		{
			fprintf(fpRead, "%s\t%d\t%f\n", strChr, nP0, vV[ni]);
			nP0 += nBinSize;
		}

		fclose(fpRead);
		
		
		sprintf(strLine, "%stilemapv2_txt2bar -i %s -d %s -o %s.cgw", strTileMapFolder, strOutFileName, strOutputPath, strOutputFile);

		system(strLine);


		sprintf(strOutFileName, "%s%s_%s.ave.bar.txt", strOutputPath, strOutputFile, strChr);
		fpRead = NULL;
		fpRead = fopen(strOutFileName, "w");
		if(fpRead == NULL)
		{
			printf("Error: cannot open output file!\n");
			exit(EXIT_FAILURE);
		}

		fprintf(fpRead, "#chr\tpos\t%s_ave_%s\n", strOutputFile,strChr);
		fprintf(fpRead, "#chr\tpos\t1\n");
		nP0 = 12;
		for(ni=0; ni<nBinNum; ni++)
		{
			fprintf(fpRead, "%s\t%d\t%f\n", strChr, nP0, vM[ni]);
			nP0 += nBinSize;
		}

		fclose(fpRead);
		
		sprintf(strLine, "%stilemapv2_txt2bar -i %s -d %s -o %s.cgw", strTileMapFolder, strOutFileName, strOutputPath, strOutputFile);
		system(strLine);

		sprintf(strOutFileName, "%s%s_%s.cor.bar.txt", strOutputPath, strOutputFile, strChr);
		fpRead = NULL;
		fpRead = fopen(strOutFileName, "w");
		if(fpRead == NULL)
		{
			printf("Error: cannot open output file!\n");
			exit(EXIT_FAILURE);
		}

		fprintf(fpRead, "#chr\tpos\t%s_cor_%s\n", strOutputFile,strChr);
		fprintf(fpRead, "#chr\tpos\t1\n");
		nP0 = 12;
		for(ni=0; ni<nBinNum; ni++)
		{
			fprintf(fpRead, "%s\t%d\t%f\n", strChr, nP0, vCorr[ni]);
			nP0 += nBinSize;
		}

		fclose(fpRead);
		
		sprintf(strLine, "%stilemapv2_txt2bar -i %s -d %s -o %s.cgw", strTileMapFolder, strOutFileName, strOutputPath, strOutputFile);
		system(strLine);

		sprintf(strOutFileName, "%s%s_%s.cov.bar.txt", strOutputPath, strOutputFile, strChr);
		fpRead = NULL;
		fpRead = fopen(strOutFileName, "w");
		if(fpRead == NULL)
		{
			printf("Error: cannot open output file!\n");
			exit(EXIT_FAILURE);
		}

		fprintf(fpRead, "#chr\tpos\t%s_cov_%s\n", strOutputFile,strChr);
		fprintf(fpRead, "#chr\tpos\t1\n");
		nP0 = 12;
		for(ni=0; ni<nBinNum; ni++)
		{
			fprintf(fpRead, "%s\t%d\t%f\n", strChr, nP0, vCov[ni]);
			nP0 += nBinSize;
		}

		fclose(fpRead);
		
		sprintf(strLine, "%stilemapv2_txt2bar -i %s -d %s -o %s.cgw", strTileMapFolder, strOutFileName, strOutputPath, strOutputFile);
		system(strLine);

		sprintf(strOutFileName, "%s%s_%s.ccc.bar.txt", strOutputPath, strOutputFile, strChr);
		fpRead = NULL;
		fpRead = fopen(strOutFileName, "w");
		if(fpRead == NULL)
		{
			printf("Error: cannot open output file!\n");
			exit(EXIT_FAILURE);
		}

		fprintf(fpRead, "#chr\tpos\t%s_ccc_%s\n", strOutputFile,strChr);
		fprintf(fpRead, "#chr\tpos\t1\n");
		nP0 = 12;
		for(ni=0; ni<nBinNum; ni++)
		{
			fprintf(fpRead, "%s\t%d\t%f\n", strChr, nP0, vC[ni]);
			nP0 += nBinSize;
		}

		fclose(fpRead);

		sprintf(strLine, "%stilemapv2_txt2bar -i %s -d %s -o %s.cgw", strTileMapFolder, strOutFileName, strOutputPath, strOutputFile);
		system(strLine);

		sprintf(strOutFileName, "%s%s_%s.cco.bar.txt", strOutputPath, strOutputFile, strChr);
		fpRead = NULL;
		fpRead = fopen(strOutFileName, "w");
		if(fpRead == NULL)
		{
			printf("Error: cannot open output file!\n");
			exit(EXIT_FAILURE);
		}

		fprintf(fpRead, "#chr\tpos\t%s_cco_%s\n", strOutputFile,strChr);
		fprintf(fpRead, "#chr\tpos\t1\n");
		nP0 = 12;
		for(ni=0; ni<nBinNum; ni++)
		{
			fprintf(fpRead, "%s\t%d\t%f\n", strChr, nP0, vC0[ni]);
			nP0 += nBinSize;
		}

		fclose(fpRead);

		sprintf(strLine, "%stilemapv2_txt2bar -i %s -d %s -o %s.cgw", strTileMapFolder, strOutFileName, strOutputPath, strOutputFile);
		system(strLine);

		sprintf(strOutFileName, "%s%s_%s.mav.bar.txt", strOutputPath, strOutputFile, strChr);
		fpRead = NULL;
		fpRead = fopen(strOutFileName, "w");
		if(fpRead == NULL)
		{
			printf("Error: cannot open output file!\n");
			exit(EXIT_FAILURE);
		}

		fprintf(fpRead, "#chr\tpos\t%s_mav_%s\n", strOutputFile,strChr);
		fprintf(fpRead, "#chr\tpos\t1\n");
		nP0 = 12;
		for(ni=0; ni<nBinNum; ni++)
		{
			fprintf(fpRead, "%s\t%d\t%f\n", strChr, nP0, vMAV[ni]);
			nP0 += nBinSize;
		}

		fclose(fpRead);
		
		sprintf(strLine, "%stilemapv2_txt2bar -i %s -d %s -o %s.cgw", strTileMapFolder, strOutFileName, strOutputPath, strOutputFile);
		system(strLine);

		sprintf(strOutFileName, "%s%s_%s.mam.bar.txt", strOutputPath, strOutputFile, strChr);
		fpRead = NULL;
		fpRead = fopen(strOutFileName, "w");
		if(fpRead == NULL)
		{
			printf("Error: cannot open output file!\n");
			exit(EXIT_FAILURE);
		}

		fprintf(fpRead, "#chr\tpos\t%s_mam_%s\n", strOutputFile,strChr);
		fprintf(fpRead, "#chr\tpos\t1\n");
		nP0 = 12;
		for(ni=0; ni<nBinNum; ni++)
		{
			fprintf(fpRead, "%s\t%d\t%f\n", strChr, nP0, vMAM[ni]);
			nP0 += nBinSize;
		}

		fclose(fpRead);
		
		sprintf(strLine, "%stilemapv2_txt2bar -i %s -d %s -o %s.cgw", strTileMapFolder, strOutFileName, strOutputPath, strOutputFile);
		system(strLine);
	}


	/* release memory */
	free(vV);
	free(vC);
	free(vC0);
	free(vM);
	free(vCorr);
	free(vCov);
	free(vMAV);
	free(vMAM);
	
	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  SeqClust_QuantileCutoff()                                              */
/*  Find quantile cutoffs, save to dCL and dCH.                            */
/* ----------------------------------------------------------------------- */ 
int SeqClust_QuantileCutoff(char strOutputPath[], char strOutputFile[], 
		int nChrNum, struct tagString **vChrName, struct INTMATRIX *pChrLen, 
		int nBinSize, double dBinMax, double dBinMin, 
		int nGridNum, double dCutL, double dCutH, double *pCL, double *pCH,
		struct BYTEMATRIX *pIncludeChr)
{
	/* define */
	int ni,nj,nk;
	int *vGrid = NULL;
	int nBinNum;
	float *vC;
	char strFileName[MED_LINE_LENGTH];
	double dStep;
	int nTotal = 0;
	int nSum;
	double dPrc;
	int nCLOK = 0;
	int nCHOK = 0;

	/* init */
	dStep = (dBinMax-dBinMin)/nGridNum;

	vGrid = (int *)calloc(nGridNum, sizeof(int));
	if(vGrid == NULL)
	{
		printf("Error: SeqClust_QuantileCutoff, cannot allocate memory for grid computation!\n");
		exit(EXIT_FAILURE);
	}

	/* process chromosome by chromosome */
	for(ni=0; ni<nChrNum; ni++)
	{
		/* skip chromosome if necessary */
		if(pIncludeChr->pMatElement[ni] == 0)
			continue;

		/* load data */
		nBinNum = pChrLen->pMatElement[ni]/nBinSize;
		if(pChrLen->pMatElement[ni] % nBinSize != 0)
			nBinNum++;

		nTotal += nBinNum;

		vC = (float *)calloc(nBinNum, sizeof(float));
		if(vC == NULL)
		{
			printf("Error: SeqClust_QuantileCutoff, cannot allocate memory for correlation vector!\n");
			exit(EXIT_FAILURE);
		}

		sprintf(strFileName, "%s%s_%s.cma", strOutputPath, strOutputFile, vChrName[ni]->m_pString);

		if( TileMapv2_LoadFromBinaryFile(vC, sizeof(float), nBinNum, strFileName) != PROC_SUCCESS)
		{
			printf("Error: SeqClust_QuantileCutoff, cannot read bin statistics correctly!\n");
			exit(EXIT_FAILURE);
		}

		/* count */
		for(nj=0; nj<nBinNum; nj++)
		{
			nk = (int)((vC[nj]-dBinMin)/dStep);
			if(nk < 0)
				nk = 0;
			if(nk >= nGridNum)
				nk = nGridNum-1;

			vGrid[nk] += 1;
		}

		/* free memory */
		free(vC);
	}

	/* find cutoff */
	nSum = 0;
	for(nj=nGridNum-1; nj>=0; nj--)
	{
		nSum += vGrid[nj];
		dPrc = (double)nSum/(double)nTotal;
		if(dPrc > (1.0 - dCutH))
		{
			if(nCHOK == 0)
			{
				*pCH = dBinMin + dStep*(nj+1);
				nCHOK = 1;
			}
		}

		if(dPrc > (1.0 - dCutL))
		{
			if(nCLOK == 0)
			{
				*pCL = dBinMin + dStep*(nj+1);
				nCLOK = 1;
			}
		}

		if((nCHOK == 1) && (nCLOK ==1))
			break;
	}

	/* free memory */
	free(vGrid);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_FDRCutoff()                                                   */
/*  Find FDR cutoffs, save to dCL and dCH.                                 */
/* ----------------------------------------------------------------------- */ 
int SeqClust_FDRCutoff(char strOutputPath[], char strOutputFile[], 
		int nChrNum, struct tagString **vChrName, struct INTMATRIX *pChrLen, 
		int nBinSize, double dBinMax, double dBinMin, 
		int nGridNum, double dCutL, double dCutH, double *pCL, double *pCH,
		struct BYTEMATRIX *pIncludeChr)
{
	/* define */
	int ni,nj,nk;
	int *vGrid = NULL;
	int *vGrid0 = NULL;
	double *vFDR = NULL;
	int nBinNum;
	float *vC;
	char strFileName[MED_LINE_LENGTH];
	double dStep;
	int nTotal = 0;
	int nSum,nSum0;
	int nCLOK = 0;
	int nCHOK = 0;

	/* init */
	dStep = (dBinMax-dBinMin)/nGridNum;

	vGrid = (int *)calloc(nGridNum, sizeof(int));
	if(vGrid == NULL)
	{
		printf("Error: SeqClust_FDRCutoff, cannot allocate memory for grid computation!\n");
		exit(EXIT_FAILURE);
	}
	vGrid0 = (int *)calloc(nGridNum, sizeof(int));
	if(vGrid0 == NULL)
	{
		printf("Error: SeqClust_FDRCutoff, cannot allocate memory for grid computation!\n");
		exit(EXIT_FAILURE);
	}
	vFDR = (double *)calloc(nGridNum, sizeof(double));
	if(vFDR == NULL)
	{
		printf("Error: SeqClust_FDRCutoff, cannot allocate memory for grid computation!\n");
		exit(EXIT_FAILURE);
	}

	/* process chromosome by chromosome */
	for(ni=0; ni<nChrNum; ni++)
	{
		/* skip chromosome if necessary */
		if(pIncludeChr->pMatElement[ni] == 0)
			continue;

		/* load data */
		nBinNum = pChrLen->pMatElement[ni]/nBinSize;
		if(pChrLen->pMatElement[ni] % nBinSize != 0)
			nBinNum++;

		nTotal += nBinNum;

		vC = (float *)calloc(nBinNum, sizeof(float));
		if(vC == NULL)
		{
			printf("Error: SeqClust_FDRCutoff, cannot allocate memory for correlation vector!\n");
			exit(EXIT_FAILURE);
		}

		/* count statistics */
		sprintf(strFileName, "%s%s_%s.cma", strOutputPath, strOutputFile, vChrName[ni]->m_pString);

		if( TileMapv2_LoadFromBinaryFile(vC, sizeof(float), nBinNum, strFileName) != PROC_SUCCESS)
		{
			printf("Error: SeqClust_FDRCutoff, cannot read bin statistics correctly!\n");
			exit(EXIT_FAILURE);
		}

		for(nj=0; nj<nBinNum; nj++)
		{
			nk = (int)((vC[nj]-dBinMin)/dStep);
			if(nk < 0)
				nk = 0;
			if(nk >= nGridNum)
				nk = nGridNum-1;

			vGrid[nk] += 1;
		}

		/* count null statistics */
		sprintf(strFileName, "%s%s_%s.cmo", strOutputPath, strOutputFile, vChrName[ni]->m_pString);

		if( TileMapv2_LoadFromBinaryFile(vC, sizeof(float), nBinNum, strFileName) != PROC_SUCCESS)
		{
			printf("Error: SeqClust_FDRCutoff, cannot read null statistics correctly!\n");
			exit(EXIT_FAILURE);
		}

		for(nj=0; nj<nBinNum; nj++)
		{
			nk = (int)((vC[nj]-dBinMin)/dStep);
			if(nk < 0)
				nk = 0;
			if(nk >= nGridNum)
				nk = nGridNum-1;

			vGrid0[nk] += 1;
		}

		/* free memory */
		free(vC);
	}

	/* find cutoff */
	nSum = 0;
	nSum0 = 0;
	for(nj=nGridNum-1; nj>=0; nj--)
	{
		nSum += vGrid[nj];
		nSum0 += vGrid0[nj];
		vFDR[nj] = (double)nSum0/(double)(nSum+1e-6);
	}
	for(nj=1; nj<nGridNum; nj++)
	{
		if(vFDR[nj] > vFDR[nj-1])
			vFDR[nj] = vFDR[nj-1];
	}
	for(nj=nGridNum-1; nj>=0; nj--)
	{
		if(vFDR[nj] > dCutH)
		{
			if(nCHOK == 0)
			{
				*pCH = dBinMin + dStep*(nj+1);
				nCHOK = 1;
			}
		}

		if(vFDR[nj] > dCutL)
		{
			if(nCLOK == 0)
			{
				*pCL = dBinMin + dStep*(nj+1);
				nCLOK = 1;
			}
		}

		if((nCHOK == 1) && (nCLOK ==1))
			break;
	}

	/* free memory */
	free(vGrid);
	free(vGrid0);
	free(vFDR);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_Segmentation()                                                */
/*  Genome segmentation. Divide genome into correlation blocks.            */
/* ----------------------------------------------------------------------- */ 
int SeqClust_Segmentation(char strOutputPath[], char strOutputFile[], 
				int nChrNum, struct tagString **vChrName, struct INTMATRIX *pChrLen, 
				int nBinSize, double dCL, double dCH)
{
	/* define */
	int ni,nj;
	int nBinNum;
	float *vC;
	char strFileName[MED_LINE_LENGTH];
	FILE *fpOut;
	char strOutFileName[MED_LINE_LENGTH];
	int nP1,nP2;
	float fMax;
	
	/* process chromosome by chromosome */
	for(ni=0; ni<nChrNum; ni++)
	{
		/* load data */
		nBinNum = pChrLen->pMatElement[ni]/nBinSize;
		if(pChrLen->pMatElement[ni] % nBinSize != 0)
			nBinNum++;

		vC = (float *)calloc(nBinNum, sizeof(float));
		if(vC == NULL)
		{
			printf("Error: SeqClust_QuantileCutoff, cannot allocate memory for correlation vector!\n");
			exit(EXIT_FAILURE);
		}

		sprintf(strFileName, "%s%s_%s.cma", strOutputPath, strOutputFile, vChrName[ni]->m_pString);

		if( TileMapv2_LoadFromBinaryFile(vC, sizeof(float), nBinNum, strFileName) != PROC_SUCCESS)
		{
			printf("Error: SeqClust_QuantileCutoff, cannot read bin statistics correctly!\n");
			exit(EXIT_FAILURE);
		}

		/* open output file */
		sprintf(strOutFileName, "%s%s_%s.tmpreg", strOutputPath, strOutputFile, vChrName[ni]->m_pString);
		fpOut = NULL;
		fpOut = fopen(strOutFileName, "w");
		if(fpOut == NULL)
		{
			printf("Error: SeqClust_QuantileCutoff, cannot open file %s to export data!\n", strOutFileName);
			exit(EXIT_FAILURE);
		}

		/* find peaks */
		nP1 = -1;
		nP2 = -1;
		fMax = -1000000.0;
		for(nj=0; nj<nBinNum; nj++)
		{
			if(vC[nj] >= dCL)
			{
				if(nP1 < 0)
				{
					nP1 = nj;
					nP2 = nj;
					fMax = vC[nj];
				}
				else
				{
					nP2 = nj;
					if(vC[nj] > fMax)
						fMax = vC[nj];
				}
			}
			else
			{
				if(nP1 >= 0)
				{
					if(fMax > dCH)
					{
						fprintf(fpOut, "%d\t%d\t%f\n", nP1, nP2, fMax);
					}
					nP1 = -1;
					nP2 = -1;
					fMax = -1000000.0;
				}
			}
		}

		if(nP1 >= 0)
		{
			if(fMax > dCH)
			{
				fprintf(fpOut, "%d\t%d\t%f\n", nP1, nP2, fMax);
			}
			nP1 = -1;
			nP2 = -1;
			fMax = -1000000.0;
		}

		/* free memory */
		free(vC);

		/* close file */
		fclose(fpOut);
	}
	
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_CollectData_ForAutoSeg()                                      */
/*  collect data from automatically determined blocks for clustering.      */
/* ----------------------------------------------------------------------- */ 
int SeqClust_CollectData_ForAutoSeg(char strOutputPath[], char strOutputFile[], 
				int nFileNum, struct tagString **vFileName, 
				struct DOUBLEMATRIX *pFileNormFactor, double dSmoothFactor,
				int nChrNum, struct tagString **vChrName, 
				struct INTMATRIX *pChrLen, int nBinSize, 
				struct BYTEMATRIX *pIncludeChr)
{
	/* define */
	int ni,nj,nk;
	FILE *fpOut;
	char strOutFileName[MED_LINE_LENGTH];
	char strFileName[MED_LINE_LENGTH];
	int nBinNum;
	struct DOUBLEMATRIX *pReg;
	double **vData;
	int nP1,nP2;
	int nBinHalf;
	double dMax;

	/* init */
	nBinHalf = nBinSize/2;

	/* open file */
	sprintf(strOutFileName, "%s%s_data.txt", strOutputPath, strOutputFile);
	fpOut = NULL;
	fpOut = fopen(strOutFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: SeqClust_CollectData_ForAutoSeg, cannot open file to export data for clustering!\n");
		exit(EXIT_FAILURE);
	}
	fprintf(fpOut, "#chr\tstart\tend\tlength\tmax_stat");
	for(nj=0; nj<nFileNum; nj++)
	{
		fprintf(fpOut, "\t%s", vFileName[nj]->m_pString);
	}
	fprintf(fpOut, "\n");

	/* process chromosome by chromosome */
	for(ni=0; ni<nChrNum; ni++)
	{
		/* skip chromosome if necessary */
		if(pIncludeChr->pMatElement[ni] == 0)
			continue;

		/* init */
		nBinNum = pChrLen->pMatElement[ni]/nBinSize;
		if(pChrLen->pMatElement[ni] % nBinSize != 0)
			nBinNum++;

		/* load data */
		sprintf(strFileName, "%s%s_%s.tmpreg", strOutputPath, strOutputFile, vChrName[ni]->m_pString);
		pReg = NULL;
		pReg = DMLOAD(strFileName);
		if(pReg == NULL)
		{
			printf("  No regions loaded for %s.\n", vChrName[ni]->m_pString);
			continue;
		}
		if(pReg->nHeight <= 0)
		{
			printf("  No regions loaded for %s.\n", vChrName[ni]->m_pString);
			DestroyDoubleMatrix(pReg);
			continue;
		}

		vData = NULL;
		vData = (double **)calloc(nFileNum, sizeof(double *));
		if(vData == NULL)
		{
			printf("Error: SeqClust_CollectData_ForAutoSeg, cannot create the vector for collecting clustering data!\n");
			exit(EXIT_FAILURE);
		}
		for(nj=0; nj<nFileNum; nj++)
		{
			vData[nj] = (double *)calloc(pReg->nHeight, sizeof(double));
			if(vData[nj] == NULL)
			{
				printf("Error: SeqClust_CollectData_ForAutoSeg, cannot create memory for collecting clustering data!\n");
				exit(EXIT_FAILURE);
			}

			sprintf(strFileName, "%s%s_%s_%s.bincount", strOutputPath, strOutputFile, 
				vFileName[nj]->m_pString, vChrName[ni]->m_pString);

			/* process file by file */
			SeqClust_CollectData_ForAutoSeg_SingleSample(pReg, vData[nj], strFileName, nBinNum);
		}
	
		/* export to files */
		for(nk = 0; nk<pReg->nHeight; nk++)
		{
			nP1 = (int)(DMGETAT(pReg, nk, 0)*nBinSize)+nBinHalf;
			nP2 = (int)(DMGETAT(pReg, nk, 1)*nBinSize)+nBinHalf;
			dMax = DMGETAT(pReg, nk, 2);

			fprintf(fpOut, "%s\t%d\t%d\t%d\t%f", vChrName[ni]->m_pString, nP1, nP2, nP2-nP1+1, dMax);

			for(nj=0; nj<nFileNum; nj++)
			{
				fprintf(fpOut, "\t%f", vData[nj][nk]*pFileNormFactor->pMatElement[nj]/(dSmoothFactor+1e-6));
			}

			fprintf(fpOut, "\n");
		}

		/* release memory */
		DestroyDoubleMatrix(pReg);
		for(nj=0; nj<nFileNum; nj++)
		{
			free(vData[nj]);
			vData[nj] = NULL;
		}
		free(vData);
	}

	/* close file */
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_CollectData_ForAutoSeg_SingleSample()                         */
/*  collect data from a single sample (file) for a region list.            */
/* ----------------------------------------------------------------------- */ 
int SeqClust_CollectData_ForAutoSeg_SingleSample(struct DOUBLEMATRIX *pReg, 
					double *pData, char strFileName[], int nBinNum)
{
	/* define */
	float *vC;
	int ni,nj;
	int nP1,nP2;
	double dSum;

	/* init */
	vC = (float *)calloc(nBinNum, sizeof(float));
	if(vC == NULL)
	{
		printf("Error: SeqClust_CollectData_ForAutoSeg_SingleSample, cannot allocate memory for correlation vector!\n");
		exit(EXIT_FAILURE);
	}

	if( TileMapv2_LoadFromBinaryFile(vC, sizeof(float), nBinNum, strFileName) != PROC_SUCCESS)
	{
		printf("Error: SeqClust_CollectData_ForAutoSeg_SingleSample, cannot read bin statistics correctly!\n");
		exit(EXIT_FAILURE);
	}

	/* collect data */
	for(ni=0; ni<pReg->nHeight; ni++)
	{
		nP1 = (int)(DMGETAT(pReg, ni, 0));
		nP2 = (int)(DMGETAT(pReg, ni, 1));
		dSum = 0.0;

		for(nj=nP1; nj<=nP2; nj++)
		{
			dSum += vC[nj];
		}

		pData[ni] = dSum;
	}

	/* release memory */
	free(vC);

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  SeqClust_CollectData_RefineFileNormFactor()                            */
/*  Recompute file normalizing factor after excluding variable blocks.     */
/* ----------------------------------------------------------------------- */ 
int SeqClust_CollectData_RefineFileNormFactor(char strOutputPath[], char strOutputFile[], 
				int nFileNum, struct tagString **vFileName, int nBaseFileId,
				struct DOUBLEMATRIX *pFileReadCount,
				struct DOUBLEMATRIX *pFileNormFactor, double dSmoothFactor,
				int nChrNum, struct tagString **vChrName, 
				struct INTMATRIX *pChrLen, int nBinSize)
{
	/* define */
	int ni,nj,nk;
	char strFileName[MED_LINE_LENGTH];
	int nBinNum;
	struct DOUBLEMATRIX *pReg;
	double **vData;

	/* process chromosome by chromosome */
	for(ni=0; ni<nChrNum; ni++)
	{
		/* init */
		nBinNum = pChrLen->pMatElement[ni]/nBinSize;
		if(pChrLen->pMatElement[ni] % nBinSize != 0)
			nBinNum++;

		/* load data */
		sprintf(strFileName, "%s%s_%s.tmpreg", strOutputPath, strOutputFile, vChrName[ni]->m_pString);
		pReg = NULL;
		pReg = DMLOAD(strFileName);
		if(pReg == NULL)
		{
			printf("  No regions loaded for %s.\n", vChrName[ni]->m_pString);
			continue;
		}
		if(pReg->nHeight <= 0)
		{
			printf("  No regions loaded for %s.\n", vChrName[ni]->m_pString);
			DestroyDoubleMatrix(pReg);
			continue;
		}

		vData = NULL;
		vData = (double **)calloc(nFileNum, sizeof(double *));
		if(vData == NULL)
		{
			printf("Error: SeqClust_CollectData_RefineFileNormFactor, cannot create the vector for collecting clustering data!\n");
			exit(EXIT_FAILURE);
		}
		for(nj=0; nj<nFileNum; nj++)
		{
			vData[nj] = (double *)calloc(pReg->nHeight, sizeof(double));
			if(vData[nj] == NULL)
			{
				printf("Error: SeqClust_CollectData_RefineFileNormFactor, cannot create memory for collecting clustering data!\n");
				exit(EXIT_FAILURE);
			}

			sprintf(strFileName, "%s%s_%s_%s.bincount", strOutputPath, strOutputFile, 
				vFileName[nj]->m_pString, vChrName[ni]->m_pString);

			/* process file by file */
			SeqClust_CollectData_ForAutoSeg_SingleSample(pReg, vData[nj], strFileName, nBinNum);
		}
	
		/* export to files */
		for(nk = 0; nk<pReg->nHeight; nk++)
		{
			for(nj=0; nj<nFileNum; nj++)
			{
				pFileReadCount->pMatElement[nj] -= vData[nj][nk];
			}
		}

		/* release memory */
		DestroyDoubleMatrix(pReg);
		for(nj=0; nj<nFileNum; nj++)
		{
			free(vData[nj]);
			vData[nj] = NULL;
		}
		free(vData);
	}

	/* recompute normalizing factor */
	printf("Adjusted File Normalizing Constant:\n");
	printf("File\tEffec_Read_No\tOld_NormFactor\tNew_NormFactor\n");
	for(ni=0; ni<nFileNum; ni++)
	{
		printf("%s\t%d\t%f\t", vFileName[ni]->m_pString, (int)(pFileReadCount->pMatElement[ni]/(dSmoothFactor+1e-6)), pFileNormFactor->pMatElement[ni]);
		pFileNormFactor->pMatElement[ni] = pFileReadCount->pMatElement[nBaseFileId]/pFileReadCount->pMatElement[ni];
		printf("%f\n", pFileNormFactor->pMatElement[ni]);
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_Main()                                                        */
/*  Model Based Clustering of ChIP-seq.                                    */
/* ----------------------------------------------------------------------- */ 
int SeqClust_Main(char strInputPath[], char strOutputPath[], char strOutputFile[],
			int nSkipCol, int nKmin, int nKmax, int nKr, int nMethod, 
			int nBmin, int nBmax, int nSeed,
			int nTransform, double dTL, int nRowStandardize, double dCut,
			int nMaxIter, double dTol)
{
	/* define */
	int nGeneNum = 0;
	int nSampleNum = 0;
	struct tagString **vInfo = NULL;
	double **vData = NULL;
	int ni;
	int nBestB = 0;
	int nBestK = 0;
	int nBestId = 0;
	int nGroupNum = 0;
	int *vGroupSize = NULL;
	int **vGroupId = NULL;


	/* ---------------------------------------------------------*/
	/* STEP1: Preprocessing                                     */
	/* ---------------------------------------------------------*/
	printf("/* ------------------------------------------------- */\n");
	printf("/* STEP 1: Preprocessing                             */\n");
	printf("/* ------------------------------------------------- */\n");
	rand_u_init(nSeed);
	AdjustDirectoryPath(strOutputPath);
	SeqClust_LoadData_Kmeans(strInputPath, nSkipCol, 
		&nGeneNum, &nSampleNum, &vInfo, &vData, 
		nTransform, dTL, nRowStandardize,
		&nGroupNum, &vGroupSize, &vGroupId);
	printf("Region number = %d\n", nGeneNum);
	printf("Sample number = %d\n", nSampleNum);


	/* ---------------------------------------------------------*/
	/* STEP2: Clustering                                        */
	/* ---------------------------------------------------------*/
	printf("\n");
	printf("/* ------------------------------------------------- */\n");
	printf("/* STEP 2: Clustering                                */\n");
	printf("/* ------------------------------------------------- */\n");

	if(nMethod <= 2)
	{
		SeqClust_FixedBlock(nGeneNum, nSampleNum, vData, 
			nKmin, nKmax, nKr, nMethod, 
			nMaxIter, dTol, 
			strOutputPath, strOutputFile,	
			nGroupNum, vGroupSize, vGroupId,
			&nBestK, &nBestId);
	}
	else if(nMethod == 3)
	{
		SeqClust_AdaptiveBlock(nGeneNum, nSampleNum, vData, 
			nKmin, nKmax, nKr, nMethod, nBmin, nBmax, 
			nMaxIter, dTol, 
			strOutputPath, strOutputFile,	
			&nBestB, &nBestK, &nBestId);
	}
	else if(nMethod == 4)
	{
		if( (nBmin < 1) || (nBmax > (nSampleNum-1)/2) )
		{
			printf("Error: Number of factors should be between 1 and %d!\n", (nSampleNum-1)/2);
			exit(EXIT_FAILURE);
		}

		SeqClust_AdaptiveFactor(nGeneNum, nSampleNum, vData, 
			nKmin, nKmax, nKr, nMethod, nBmin, nBmax, 
			nMaxIter, dTol, 
			strOutputPath, strOutputFile,	
			&nBestB, &nBestK, &nBestId);
	}

	/* ---------------------------------------------------------*/
	/* STEP4: Hierarchical clustering of samples and clusters   */
	/* ---------------------------------------------------------*/
	printf("\n");
	printf("/* ------------------------------------------------- */\n");
	printf("/* STEP 4: Cluster samples and clusters              */\n");
	printf("/* ------------------------------------------------- */\n");

	/* ---------------------------------------------------------*/
	/* STEP5: Export results                                    */
	/* ---------------------------------------------------------*/
	printf("\n");
	printf("/* ------------------------------------------------- */\n");
	printf("/* STEP 5: Export results                            */\n");
	printf("/* ------------------------------------------------- */\n");

	/* ---------------------------------------------------------*/
	/* STEP6: Clear memory                                      */
	/* ---------------------------------------------------------*/
	for(ni=0; ni<nGeneNum; ni++)
	{
		free(vData[ni]);
		DeleteString(vInfo[ni]);
	}
	free(vData);
	free(vInfo);
	if(nGroupNum > 0)
	{
		for(ni=0; ni<nGroupNum; ni++)
		{
			free(vGroupId[ni]);
			vGroupId[ni] = NULL;
		}
		free(vGroupId);
		free(vGroupSize);
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_FixedBlock()                                                  */
/*  Fixed covariance block clustering.                                     */
/* ----------------------------------------------------------------------- */ 
int SeqClust_FixedBlock(int nGeneNum, int nSampleNum, double **vData, 
			int nKmin, int nKmax, int nKr, int nMethod, 
			int nMaxIter, double dTol, 
			char strOutputPath[], char strOutputFile[],	
			int nGroupNum, int *vGroupSize, int **vGroupId,
			int *pBestK, int *pBestId)
{
	/* define */
	int ni,nj,nk;
	double *vBIC = NULL;
	int *vBID = NULL;
	double dResult = 0.0;
	double dBestBIC = 0.0;
	int nBestK = 0;
	int nBestId = 0;
	char strFileName[MED_LINE_LENGTH];
	
	vBIC = (double *)calloc(nKmax-nKmin+1, sizeof(double));
	vBID = (int *)calloc(nKmax-nKmin+1, sizeof(int));
	if( (vBIC == NULL) || (vBID == NULL) )
	{
		printf("Error: SeqClust_Main, cannot create memory for monitering BIC.\n");
		exit(EXIT_FAILURE);
	}
	for(nk=nKmin; nk<=nKmax; nk++)
		vBIC[nk-nKmin] = 1e6;

	for(nk=nKmin; nk<=nKmax; nk++)
	{
		nj = nk-nKmin;
		for(ni=0; ni<nKr; ni++)
		{
			sprintf(strFileName, "%s_K%d_T%d", strOutputFile, nk, ni);

			printf(" K=%d, Try=%d ... \n", nk, ni);
			if(nMethod == 1)
			{
				dResult = SeqClust_Kpolya(nGeneNum, nSampleNum, vData, nk, 
					nMaxIter, dTol, strOutputPath, strFileName);
			}
			else if(nMethod == 2)
			{
				dResult = SeqClust_KMVnorm(nGeneNum, nSampleNum, vData, nk, 
					nMaxIter, dTol, strOutputPath, strFileName, 
					nGroupNum, vGroupSize, vGroupId);
			}
			else
			{
				dResult = SeqClust_Knorm(nGeneNum, nSampleNum, vData, nk, 
					nMaxIter, dTol, strOutputPath, strFileName);
			}

			if(ni == 0)
			{
				vBIC[nj] = dResult;
				vBID[nj] = ni;
			}
			else if(dResult < vBIC[nj])
			{
				vBIC[nj] = dResult;
				vBID[nj] = ni;
			}
		}

		if(nk==nKmin)
		{
			dBestBIC = vBIC[nj];
			nBestK = nk;
			nBestId = vBID[nj];
		}
		else if(vBIC[nj] < dBestBIC)
		{
			dBestBIC = vBIC[nj];
			nBestK = nk;
			nBestId = vBID[nj];
		}
		else
		{
			if((nk-nBestK) >= 3)
				break;
		}
	}

	printf("\n");
	printf("/* ------------------------------------------------- */\n");
	printf("/* STEP 3: Find the best cluster result              */\n");
	printf("/* ------------------------------------------------- */\n");
	
	printf("Clust_Num\tBIC\tBest_Trial_Id\n");
	for(nk=nKmin; nk<=nKmax; nk++)
	{
		nj = nk-nKmin;
		printf("%d\t%f\t%d\n", nk, vBIC[nj], vBID[nj]);
	}
	printf("Optimal Cluster Number = %d\n", nBestK);
	printf("Optimal Cluster Trial = %d\n", nBestId);

	*pBestK = nBestK;
	*pBestId = nBestId;

	/* release memory */
	free(vBIC);
	free(vBID);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_AdaptiveBlock()                                               */
/*  Adaptive covariance block clustering.                                  */
/* ----------------------------------------------------------------------- */ 
int SeqClust_AdaptiveBlock(int nGeneNum, int nSampleNum, double **vData, 
			int nKmin, int nKmax, int nKr, int nMethod, int nBmin, int nBmax,
			int nMaxIter, double dTol, 
			char strOutputPath[], char strOutputFile[],	
			int *pBestB, int *pBestK, int *pBestId)
{
	/* define */
	int ni,nk,nx,ny;
	double *vBIC = NULL;
	int *vBID = NULL;
	double dResult = 0.0;
	double dBestBIC = 0.0;
	int nBestB = 0;
	int nBestK = 0;
	int nBestId = 0;
	char strFileName[MED_LINE_LENGTH];
	int nGroupNum = 0;
	int *vGroupSize = NULL;
	int **vGroupId = NULL;
	int nDistanceType = 1; /* 0=correlation, 1=absolute correlation*/
	int nMergeType = 0; /* average linkage */
	int nExportCluster = 1;

	/* step1: cluster samples */
	printf("Cluster samples ... \n");
	HierarchicalCluster_Column(nGeneNum, nSampleNum, vData, nDistanceType, nMergeType, 
		nExportCluster, strOutputPath, strOutputFile, nBmin, nBmax);

	/* step2: cluster for each block configuration */
	vBIC = (double *)calloc((nKmax-nKmin+1)*(nBmax-nBmin+1), sizeof(double));
	vBID = (int *)calloc((nKmax-nKmin+1)*(nBmax-nBmin+1), sizeof(int));
	if( (vBIC == NULL) || (vBID == NULL) )
	{
		printf("Error: SeqClust_Main, cannot create memory for monitering BIC.\n");
		exit(EXIT_FAILURE);
	}

	ny = 0;
	for(ni=nBmin; ni<=nBmax; ni++)
	{
		for(nk=nKmin; nk<=nKmax; nk++)
		{
			vBIC[ny] = 1e6;
			ny++;
		}
	}

	ny = 0;
	for(ni=nBmin; ni<=nBmax; ni++)
	{
		sprintf(strFileName, "%s%s_sampleblock_B%d.txt", strOutputPath, strOutputFile, ni);
		SeqClust_ReadCovBlock(strFileName, &nGroupNum, &vGroupSize, &vGroupId);

		for(nk=nKmin; nk<=nKmax; nk++)
		{
			for(nx=0; nx<nKr; nx++)
			{
				sprintf(strFileName, "%s_B%d_K%d_T%d", strOutputFile, ni, nk, nx);

				printf(" B=%d, K=%d, Try=%d ... \n", ni, nk, nx);
				dResult = SeqClust_KMVnorm(nGeneNum, nSampleNum, vData, nk, 
					nMaxIter, dTol, strOutputPath, strFileName, 
					nGroupNum, vGroupSize, vGroupId);
				
				if(nx == 0)
				{
					vBIC[ny] = dResult;
					vBID[ny] = nx;
				}
				else if(dResult < vBIC[ny])
				{
					vBIC[ny] = dResult;
					vBID[ny] = nx;
				}
			}

			if(ny==0)
			{
				dBestBIC = vBIC[ny];
				nBestB = ni;
				nBestK = nk;
				nBestId = vBID[ny];
			}
			else if(vBIC[ny] < dBestBIC)
			{
				dBestBIC = vBIC[ny];
				nBestB = ni;
				nBestK = nk;
				nBestId = vBID[ny];
			}

			ny++;
		}

		for(nk=0; nk<nGroupNum; nk++)
		{
			free(vGroupId[nk]);
			vGroupId[nk] = NULL;
		}
		free(vGroupId);
		free(vGroupSize);
		vGroupId = NULL;
		vGroupSize = NULL;
	}

	printf("\n");
	printf("/* ------------------------------------------------- */\n");
	printf("/* STEP 3: Find the best cluster result              */\n");
	printf("/* ------------------------------------------------- */\n");
	
	printf("Block_Num\tClust_Num\tBIC\tBest_Trial_Id\n");
	ny = 0;
	for(ni=nBmin; ni<=nBmax; ni++)
	{
		for(nk=nKmin; nk<=nKmax; nk++)
		{
			printf("%d\t%d\t%f\t%d\n", ni, nk, vBIC[ny], vBID[ny]);
			ny++;
		}
	}
	printf("Optimal Block Number = %d\n", nBestB);
	printf("Optimal Cluster Number = %d\n", nBestK);
	printf("Optimal Cluster Trial = %d\n", nBestId);

	*pBestB = nBestB;
	*pBestK = nBestK;
	*pBestId = nBestId;

	/* release memory */
	free(vBIC);
	free(vBID);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_AdaptiveFactor()                                              */
/*  Adaptive mixture factor model based clustering.                        */
/* ----------------------------------------------------------------------- */ 
int SeqClust_AdaptiveFactor(int nGeneNum, int nSampleNum, double **vData, 
			int nKmin, int nKmax, int nKr, int nMethod, int nBmin, int nBmax,
			int nMaxIter, double dTol, 
			char strOutputPath[], char strOutputFile[],	
			int *pBestB, int *pBestK, int *pBestId)
{
	/* define */
	int ni,nk,nx,ny;
	double *vBIC = NULL;
	int *vBID = NULL;
	double dResult = 0.0;
	double dBestBIC = 0.0;
	int nBestB = 0;
	int nBestK = 0;
	int nBestId = 0;
	char strFileName[MED_LINE_LENGTH];
	int nExportCluster = 1;

	/* step1: cluster for each block configuration */
	vBIC = (double *)calloc((nKmax-nKmin+1)*(nBmax-nBmin+1), sizeof(double));
	vBID = (int *)calloc((nKmax-nKmin+1)*(nBmax-nBmin+1), sizeof(int));
	if( (vBIC == NULL) || (vBID == NULL) )
	{
		printf("Error: SeqClust_Main, cannot create memory for monitering BIC.\n");
		exit(EXIT_FAILURE);
	}

	ny = 0;
	for(ni=nBmin; ni<=nBmax; ni++)
	{
		for(nk=nKmin; nk<=nKmax; nk++)
		{
			vBIC[ny] = 1e6;
			ny++;
		}
	}

	ny = 0;
	for(ni=nBmin; ni<=nBmax; ni++)
	{
		for(nk=nKmin; nk<=nKmax; nk++)
		{
			for(nx=0; nx<nKr; nx++)
			{
				sprintf(strFileName, "%s_F%d_K%d_T%d", strOutputFile, ni, nk, nx);

				printf(" F=%d, K=%d, Try=%d ... \n", ni, nk, nx);
				dResult = SeqClust_KMFnorm(nGeneNum, nSampleNum, vData, ni, nk, 
					nMaxIter, dTol, strOutputPath, strFileName);
				
				if(nx == 0)
				{
					vBIC[ny] = dResult;
					vBID[ny] = nx;
				}
				else if(dResult < vBIC[ny])
				{
					vBIC[ny] = dResult;
					vBID[ny] = nx;
				}
			}

			if(ny==0)
			{
				dBestBIC = vBIC[ny];
				nBestB = ni;
				nBestK = nk;
				nBestId = vBID[ny];
			}
			else if(vBIC[ny] < dBestBIC)
			{
				dBestBIC = vBIC[ny];
				nBestB = ni;
				nBestK = nk;
				nBestId = vBID[ny];
			}

			ny++;
		}
	}

	printf("\n");
	printf("/* ------------------------------------------------- */\n");
	printf("/* STEP 3: Find the best cluster result              */\n");
	printf("/* ------------------------------------------------- */\n");
	
	printf("Factor_Num\tClust_Num\tBIC\tBest_Trial_Id\n");
	ny = 0;
	for(ni=nBmin; ni<=nBmax; ni++)
	{
		for(nk=nKmin; nk<=nKmax; nk++)
		{
			printf("%d\t%d\t%f\t%d\n", ni, nk, vBIC[ny], vBID[ny]);
			ny++;
		}
	}
	printf("Optimal Factor Number = %d\n", nBestB);
	printf("Optimal Cluster Number = %d\n", nBestK);
	printf("Optimal Cluster Trial = %d\n", nBestId);

	*pBestB = nBestB;
	*pBestK = nBestK;
	*pBestId = nBestId;

	/* release memory */
	free(vBIC);
	free(vBID);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_ReadCovBlock()                                                */
/*  Load covariance blocks for clustering analysis.                        */
/* ----------------------------------------------------------------------- */ 
int SeqClust_ReadCovBlock(char strFileName[], int *pGroupNum, 
						  int **vGroupSize, int ***vGroupId)
{
	/* define */
	FILE *fpIn;
	char strLine[LONG_LINE_LENGTH];
	int ni,nj;
	
	/* open files */
	fpIn = NULL;
	fpIn = fopen(strFileName, "r");
	if(fpIn == NULL)
	{
		printf("Error: SeqClust_ReadCovBlock, cannot open covariance block file!\n");
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

		*pGroupNum = atoi(strLine);
		
		*vGroupSize = NULL;
		*vGroupSize = (int *)calloc(*pGroupNum, sizeof(int));
		if(*vGroupSize == NULL)
		{
			printf("Error: SeqClust_ReadCovBlock, cannot create group size vector\n");
			exit(EXIT_FAILURE);
		}

		*vGroupId = NULL;
		*vGroupId = (int **)calloc(*pGroupNum, sizeof(int *));
		if(*vGroupId == NULL)
		{
			printf("Error: SeqClust_ReadCovBlock, cannot create group id vector\n");
			exit(EXIT_FAILURE);
		}

		for(ni=0; ni<(*pGroupNum); ni++)
		{
			if(fgets(strLine, LONG_LINE_LENGTH, fpIn) == NULL)
			{
				printf("Error: SeqClust_ReadCovBlock, corrupt file\n");
				exit(EXIT_FAILURE);
			}
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] != '>')
			{
				printf("Error: SeqClust_ReadCovBlock, corrupt file, no cluster indicator\n");
				exit(EXIT_FAILURE);
			}
			(*vGroupSize)[ni] = atoi(strLine+1);

			(*vGroupId)[ni] = (int *)calloc((*vGroupSize)[ni], sizeof(int));
			if((*vGroupId)[ni] == NULL)
			{
				printf("Error: SeqClust_ReadCovBlock, cannot create group id vector\n");
				exit(EXIT_FAILURE);
			}

			for(nj=0; nj<(*vGroupSize)[ni]; nj++)
			{
				if(fgets(strLine, LONG_LINE_LENGTH, fpIn) == NULL)
				{
					printf("Error: SeqClust_ReadCovBlock, corrupt file, cannot find sample id\n");
					exit(EXIT_FAILURE);
				}
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				(*vGroupId)[ni][nj] = atoi(strLine);
			}
		}
	}

	/* close file */
	fclose(fpIn);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_LoadData_Kmeans()                                             */
/*  Load data for K-means clustering.                                      */
/* ----------------------------------------------------------------------- */ 
int SeqClust_LoadData_Kmeans(char strInputPath[], int nSkipCol, 
		int *pGeneNum, int *pSampleNum, 
		struct tagString ***vInfo, double ***vData, 
		int nTransform, double dTL, int nRowStandardize,
		int *pGroupNum, int **vGroupSize, int ***vGroupId)
{
	/* define */
	FILE *fpIn;
	char strLine[LONG_LINE_LENGTH];
	int nGeneNum = 0;
	int nSampleNum = 0;
	int nFirstLine = 1;
	char *chp1,*chp2;
	int ni,nj,nk;
	double dLog2 = log(2.0);
	double dRowM,dRowSD,dTemp;
	struct DOUBLEMATRIX *pGID;
	struct DOUBLEMATRIX *pGIDSort;
	struct LONGMATRIX *pGIDSid;

	/* obtain gene number and sample number */
	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: SeqClust_LoadData_Kmeans, cannot open the input file! \n");
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
		
		if(nFirstLine == 1)
		{
			chp1 = strLine;
			chp2 = strchr(chp1, '\t');
			while(chp2 != NULL)
			{
				nSampleNum++;
				chp1 = chp2+1;
				chp2 = strchr(chp1, '\t');
			}
			nSampleNum = nSampleNum+1-nSkipCol;
			nFirstLine = 0;
		}

		if(strstr(strLine, "GROUPID") == strLine)
		{
			pGIDSort = NULL;
			pGIDSid = NULL;
			pGID = NULL;
			pGID = CreateDoubleMatrix(1, nSampleNum);
			if(pGID == NULL)
			{
				printf("Error: SeqClust_LoadData_Kmeans, cannot create memory for loading group ids\n");
				exit(EXIT_FAILURE);
			}

			chp1 = strLine;
			for(nj=0; nj<nSkipCol; nj++)
			{
				chp2 = strchr(chp1, '\t');
				chp1 = chp2+1;
			}

			nj = 0;
			chp2 = strchr(chp1, '\t');
			while(chp2 != NULL)
			{
				*chp2 = '\0';
				pGID->pMatElement[nj] = atof(chp1);

				chp1 = chp2+1;
				chp2 = strchr(chp1, '\t');
				nj++;
			}

			pGID->pMatElement[nj] = atof(chp1);
			nj++;
			if(nj != nSampleNum)
			{
				printf("Error: SeqClust_LoadData_Kmeans, sample number does not match\n");
				exit(EXIT_FAILURE);
			}

			DMSORTMERGEA_0(pGID, &pGIDSort, &pGIDSid);

			*pGroupNum = 1;
			for(nj=1; nj<nSampleNum; nj++)
			{
				if((int)(pGIDSort->pMatElement[nj]) != (int)(pGIDSort->pMatElement[nj-1]))
					*pGroupNum += 1;
			}
			
			*vGroupSize = NULL;
			*vGroupSize = (int *)calloc(*pGroupNum, sizeof(int));
			if(*vGroupSize == NULL)
			{
				printf("Error: SeqClust_LoadData_Kmeans, cannot create group size vector\n");
				exit(EXIT_FAILURE);
			}

			ni = 0;
			nk = 1;
			for(nj=1; nj<nSampleNum; nj++)
			{
				if((int)(pGIDSort->pMatElement[nj]) == (int)(pGIDSort->pMatElement[nj-1]))
				{
					nk++;
				}
				else
				{
					(*vGroupSize)[ni] = nk;
					ni++;
					nk = 1;
				}
			}

			(*vGroupSize)[ni] = nk;
			ni++;
			if(ni != (*pGroupNum))
			{
				printf("Error: SeqClust_LoadData_Kmeans, group number does not match\n");
				exit(EXIT_FAILURE);
			}

			*vGroupId = NULL;
			*vGroupId = (int **)calloc(*pGroupNum, sizeof(int *));
			if(*vGroupId == NULL)
			{
				printf("Error: SeqClust_LoadData_Kmeans, cannot create group id vector\n");
				exit(EXIT_FAILURE);
			}
			for(ni = 0; ni<(*pGroupNum); ni++)
			{
				(*vGroupId)[ni] = (int *)calloc((*vGroupSize)[ni], sizeof(int));
				if((*vGroupId)[ni] == NULL)
				{
					printf("Error: SeqClust_LoadData_Kmeans, cannot create group id vector\n");
					exit(EXIT_FAILURE);
				}
			}

			ni = 0;
			(*vGroupId)[ni][0] = pGIDSid->pMatElement[0];
			nk = 1;

			for(nj=1; nj<nSampleNum; nj++)
			{
				if((int)(pGIDSort->pMatElement[nj]) == (int)(pGIDSort->pMatElement[nj-1]))
				{
					(*vGroupId)[ni][nk] = pGIDSid->pMatElement[nj];
					nk++;
				}
				else
				{
					ni++;
					(*vGroupId)[ni][0] = pGIDSid->pMatElement[nj];
					nk = 1;
				}
			}
			
			DestroyDoubleMatrix(pGID);
			DestroyDoubleMatrix(pGIDSort);
			DestroyLongMatrix(pGIDSid);
			continue;
		}

		nGeneNum++;
	}

	fclose(fpIn);

	*pGeneNum = nGeneNum;
	*pSampleNum = nSampleNum;

	/* create memory */
	*vInfo = NULL;
	*vInfo = (struct tagString **)calloc(nGeneNum, sizeof(struct tagString *));
	if(*vInfo == NULL)
	{
		printf("Error: SeqClust_LoadData_Kmeans, cannot create memory for storing information! \n");
		exit(EXIT_FAILURE);
	}

	*vData = NULL;
	*vData = (double **)calloc(nGeneNum, sizeof(double *));
	if(*vData == NULL)
	{
		printf("Error: SeqClust_LoadData_Kmeans, cannot create memory for storing data! \n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nGeneNum; ni++)
	{
		(*vData)[ni] = (double *)calloc(nSampleNum, sizeof(double));
		if((*vData)[ni] == NULL)
		{
			printf("Error: SeqClust_LoadData_Kmeans, cannot create memory for storing data! \n");
			exit(EXIT_FAILURE);
		}
	}


	/* load data */
	ni = 0;

	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: SeqClust_LoadData_Kmeans, cannot open the input file! \n");
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
		if(strstr(strLine, "GROUPID") == strLine)
			continue;

		chp1 = strLine;
		for(nj=0; nj<nSkipCol; nj++)
		{
			chp2 = strchr(chp1, '\t');
			chp1 = chp2+1;
		}

		nj = 0;
		chp2 = strchr(chp1, '\t');
		while(chp2 != NULL)
		{
			*chp2 = '\0';
			(*vData)[ni][nj] = atof(chp1);

			chp1 = chp2+1;
			chp2 = strchr(chp1, '\t');
			nj++;
		}

		(*vData)[ni][nj] = atof(chp1);
		nj++;

		if(nj != nSampleNum)
		{
			printf("Error: SeqClust_LoadData_Kmeans, sample number not consistent! \n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	fclose(fpIn);
	
	if(ni != nGeneNum)
	{
		printf("Error: SeqClust_LoadData_Kmeans, gene number not consistent! \n");
		exit(EXIT_FAILURE);
	}

	/* preprocessing */
	if(nTransform == 1)
	{
		for(ni=0; ni<nGeneNum; ni++)
		{
			for(nj=0; nj<nSampleNum; nj++)
			{
				if((*vData)[ni][nj] < dTL)
					(*vData)[ni][nj] = dTL;
				(*vData)[ni][nj] = log((*vData)[ni][nj])/dLog2;
			}
		}
	}
	
	if(nRowStandardize == 1)
	{
		for(ni=0; ni<nGeneNum; ni++)
		{
			dRowM = 0.0;
			for(nj=0; nj<nSampleNum; nj++)
			{
				dRowM += (*vData)[ni][nj];
			}
			dRowM /= nSampleNum;

			dRowSD = 0.0;
			for(nj=0; nj<nSampleNum; nj++)
			{
				dTemp = (*vData)[ni][nj]-dRowM;
				dRowSD += dTemp*dTemp;
			}
			if(nSampleNum > 1)
				dRowSD /= (nSampleNum-1);
			else
				dRowSD = 0;

			dRowSD = sqrt(dRowSD)+1e-6;

			for(nj=0; nj<nSampleNum; nj++)
			{
				(*vData)[ni][nj] = ((*vData)[ni][nj]-dRowM)/dRowSD;
			}
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_Knorm()                                                       */
/*  Model based clustering (K normal distributions with diagonal COV       */
/* ----------------------------------------------------------------------- */ 
double SeqClust_Knorm(int nGeneNum, int nSampleNum, double **vData, 
				int nK, int nMaxIter, double dTol, 
				char strOutputPath[], char strOutputFile[])
{
	/* define */
	double dBIC = 0.0;

	double *vCP = NULL;
	double **vCMean = NULL;
	double **vCVar = NULL;
	
	double *vCPNew = NULL;
	double **vCMeanNew = NULL;
	double **vCVarNew = NULL;

	double **vLike = NULL;
	double dLike = 0.0;
	double dLikeNew = 0.0;

	double dErr = 1e6;
	int nIter = 0;
	int nFinal = 0;

	int ni,nj,nk;
	double dRand;
	double dpi = 4.0*atan(1.0);
	double dnormfactor = -0.5*log(2*dpi);

	double dLMax,dLSum,dRErr;

	FILE *fpOut;
	char strFileName[MED_LINE_LENGTH];

	/* Initialize */
	vCP = (double *)calloc(nK, sizeof(double));
	vCPNew = (double *)calloc(nK, sizeof(double));
	if( (vCP == NULL) || (vCPNew == NULL) )
	{
		printf("Error: SeqClust_Knorm, cannot create memory for prior probabilities of clusters\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nK; ni++)
	{
		vCP[ni] = 1.0/(double)nK;
		vCPNew[ni] = vCP[ni];
	}

	vCMean = (double **)calloc(nK, sizeof(double *));
	vCMeanNew = (double **)calloc(nK, sizeof(double *));
	vCVar = (double **)calloc(nK, sizeof(double *));
	vCVarNew = (double **)calloc(nK, sizeof(double *));
	if( (vCMean == NULL) || (vCMeanNew == NULL) || (vCVar == NULL) || (vCVarNew == NULL) )
	{
		printf("Error: SeqClust_Knorm, cannot create memory for cluster mean and variance\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nK; ni++)
	{
		vCMean[ni] = (double *)calloc(nSampleNum, sizeof(double));
		vCMeanNew[ni] = (double *)calloc(nSampleNum, sizeof(double));
		vCVar[ni] = (double *)calloc(nSampleNum, sizeof(double));
		vCVarNew[ni] = (double *)calloc(nSampleNum, sizeof(double));
		if( (vCMean[ni] == NULL) || (vCMeanNew[ni] == NULL) || (vCVar[ni] == NULL) || (vCVarNew[ni] == NULL) )
		{
			printf("Error: SeqClust_Knorm, cannot create memory for cluster mean and variance\n");
			exit(EXIT_FAILURE);
		}
		for(nj=0; nj<nSampleNum; nj++)
		{
			vCVar[ni][nj] = 1.0;
			vCVarNew[ni][nj] = 1.0;
		}

		dRand = rand_u();
		nk = (int)(dRand*nGeneNum);
		if(nk >= nGeneNum)
			nk = nGeneNum-1;
		if(nk < 0)
			nk = 0;
		for(nj=0; nj<nSampleNum; nj++)
		{
			vCMean[ni][nj] = vData[nk][nj];
			vCMeanNew[ni][nj] = vData[nk][nj];
		}
	}

	vLike = (double **)calloc(nGeneNum, sizeof(double *));
	if( vLike == NULL )
	{
		printf("Error: SeqClust_Knorm, cannot create memory for computing likelihood\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nGeneNum; ni++)
	{
		vLike[ni] = (double *)calloc(nK, sizeof(double));
		if( vLike[ni] == NULL )
		{
			printf("Error: SeqClust_Knorm, cannot create memory for computing likelihood\n");
			exit(EXIT_FAILURE);
		}
	}

	/* EM */
	for(nIter=0; nIter<nMaxIter; nIter++)
	{
		/* prepare parameter */
		for(nk=0; nk<nK; nk++)
			vCPNew[nk] = log(vCP[nk]);
		for(nk=0; nk<nK; nk++)
		{
			for(nj=0; nj<nSampleNum; nj++)
			{
				vCVarNew[nk][nj] = dnormfactor-0.5*log(vCVar[nk][nj]);
			}
		}

		/* compute likelihood & posterior */
		dLikeNew = 0.0;
		for(ni=0; ni<nGeneNum; ni++)
		{
			for(nk=0; nk<nK; nk++)
			{
				vLike[ni][nk] = vCPNew[nk];
				for(nj=0; nj<nSampleNum; nj++)
				{
					vLike[ni][nk] += vCVarNew[nk][nj]-(vData[ni][nj]-vCMean[nk][nj])*(vData[ni][nj]-vCMean[nk][nj])/2/vCVar[nk][nj];
				}

				if(nk == 0)
					dLMax = vLike[ni][nk];
				else if(vLike[ni][nk] > dLMax)
					dLMax = vLike[ni][nk];
			}

			dLSum = 0.0;
			for(nk=0; nk<nK; nk++)
			{
				vLike[ni][nk] = exp(vLike[ni][nk]-dLMax);
				dLSum += vLike[ni][nk];
			}

			dLikeNew += dLMax+log(dLSum);

			for(nk=0; nk<nK; nk++)
			{
				vLike[ni][nk] = vLike[ni][nk]/dLSum;
			}
		}

		if(nFinal == 1)
		{
			dLike = dLikeNew;
			break;
		}

		/* update parameters */
		for(nk=0; nk<nK; nk++)
		{
			vCPNew[nk] = 0.0;
			for(nj=0; nj<nSampleNum; nj++)
			{
				vCMeanNew[nk][nj] = 0.0;
				vCVarNew[nk][nj] = 0.0;
			}
		}
		for(ni=0; ni<nGeneNum; ni++)
		{
			for(nk=0; nk<nK; nk++)
			{
				vCPNew[nk] += vLike[ni][nk];
				for(nj=0; nj<nSampleNum; nj++)
				{
					vCMeanNew[nk][nj] += vLike[ni][nk]*vData[ni][nj];
					vCVarNew[nk][nj] += vLike[ni][nk]*vData[ni][nj]*vData[ni][nj];
				}
			}
		}

		dErr = 0.0;
		for(nk=0; nk<nK; nk++)
		{
			for(nj=0; nj<nSampleNum; nj++)
			{
				vCMeanNew[nk][nj]  /= vCPNew[nk];
				dRErr = fabs( (vCMeanNew[nk][nj]-vCMean[nk][nj])/(vCMean[nk][nj]+1e-6));
				if(dRErr > dErr)
					dErr = dRErr;

				vCVarNew[nk][nj]  = vCVarNew[nk][nj]/vCPNew[nk]-vCMeanNew[nk][nj]*vCMeanNew[nk][nj];
				dRErr = fabs( (vCVarNew[nk][nj]-vCVar[nk][nj])/(vCVar[nk][nj]+1e-6));
				if(dRErr > dErr)
					dErr = dRErr;
			}

			vCPNew[nk] /= nGeneNum;
			dRErr = fabs( (vCPNew[nk]-vCP[nk])/(vCP[nk]+1e-6) );
			if(dRErr > dErr)
				dErr = dRErr;
		}

				
		/* evaluate stopping criteria */
		if(nIter > 0)
		{
			if(dLikeNew < dLike)
			{
				if( fabs(dLikeNew-dLike)/(dLike+1e-6) > 1e-6)
					printf("Warning: decreasing likelihood!\n");
				dLike = dLikeNew;
				break;
			}
		}

		if(dErr < dTol)
		{
			nFinal = 1;
			nIter--;
		}

		for(nk=0; nk<nK; nk++)
		{
			for(nj=0; nj<nSampleNum; nj++)
			{
				vCMean[nk][nj] = vCMeanNew[nk][nj];
				vCVar[nk][nj]  = vCVarNew[nk][nj];
			}
			vCP[nk] = vCPNew[nk];
			dLike = dLikeNew;
		}
	}

	/* Export results */
	sprintf(strFileName, "%s%s.pp", strOutputPath, strOutputFile);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: SeqClust_Knorm, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nGeneNum; ni++)
	{
		fprintf(fpOut, "%d", ni);
		for(nk=0; nk<nK; nk++)
		{
			fprintf(fpOut, "\t%f", vLike[ni][nk]); 
		}
		fprintf(fpOut, "\n");
	}
	fclose(fpOut);

	sprintf(strFileName, "%s%s.prior", strOutputPath, strOutputFile);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: SeqClust_Knorm, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	for(nk=0; nk<nK; nk++)
	{
		if(nk != 0)
			fprintf(fpOut, "\t");

		fprintf(fpOut, "%f", vCP[nk]); 
	}
	fprintf(fpOut, "\n");
	fclose(fpOut);

	sprintf(strFileName, "%s%s.mu", strOutputPath, strOutputFile);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: SeqClust_Knorm, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	for(nk=0; nk<nK; nk++)
	{
		for(nj=0; nj<nSampleNum; nj++)
		{
			if(nj != 0)
				fprintf(fpOut, "\t");
			fprintf(fpOut, "%f", vCMean[nk][nj]); 
		}
		fprintf(fpOut, "\n");
	}
	fclose(fpOut);

	sprintf(strFileName, "%s%s.va", strOutputPath, strOutputFile);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: SeqClust_Knorm, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	for(nk=0; nk<nK; nk++)
	{
		for(nj=0; nj<nSampleNum; nj++)
		{
			if(nj != 0)
				fprintf(fpOut, "\t");
			fprintf(fpOut, "%f", vCVar[nk][nj]); 
		}
		fprintf(fpOut, "\n");
	}
	fclose(fpOut);

	/* Release memory */
	for(ni=0; ni<nGeneNum; ni++)
	{
		free(vLike[ni]);
	}
	free(vLike);
	free(vCP);
	free(vCPNew);
	for(ni=0; ni<nK; ni++)
	{
		free(vCMean[ni]);
		free(vCMeanNew[ni]);
		free(vCVar[ni]);
		free(vCVarNew[ni]);
	}
	free(vCMean);
	free(vCMeanNew);
	free(vCVar);
	free(vCVarNew);

	/* return */
	dBIC = -2.0*dLike+(nK-1+nK*nSampleNum*2)*log((double)(nGeneNum));
	printf(" iteration number = %d;  relative err = %f; BIC = %f \n", nIter, dErr, dBIC);

	return dBIC;
}


/* ----------------------------------------------------------------------- */ 
/*  SeqClust_KMVnorm()                                                     */
/*  Model based clustering (K normal distributions with block covariance   */
/*  matrix.                                                                */
/* ----------------------------------------------------------------------- */ 
double SeqClust_KMVnorm(int nGeneNum, int nSampleNum, double **vData, 
				int nK, int nMaxIter, double dTol, 
				char strOutputPath[], char strOutputFile[],
				int nGroupNum, int *vGroupSize, int **vGroupId)
{
	/* define */
	double dBIC = 0.0;

	double *vCP = NULL;
	double *vCPNew = NULL;

	double **vCMean = NULL;
	double **vCMeanNew = NULL;
	
	struct DOUBLEMATRIX ***vCVar = NULL;
	struct DOUBLEMATRIX ***vCVarNew = NULL;

	struct DOUBLEMATRIX *pD;
	struct DOUBLEMATRIX *pV;
	double *vTemp1;

	double **vLike = NULL;
	double dLike = 0.0;
	double dLikeNew = 0.0;

	int nSingularExit = 0;

	double dErr = 1e6;
	int nIter = 0;
	int nFinal = 0;

	double *pinvV,*pinvV2;
	int ni,nj,nk,nx,ny,nu,nw;
	double dRand;
	double dpi = 4.0*atan(1.0);
	double dnormfactor;

	double dLMax,dLSum,dRErr;
	int nParamNum;

	FILE *fpOut;
	char strFileName[MED_LINE_LENGTH];

	/* Initialize */
	dnormfactor = -nSampleNum*log(2*dpi)/2.0;
	vCP = (double *)calloc(nK, sizeof(double));
	vCPNew = (double *)calloc(nK, sizeof(double));
	if( (vCP == NULL) || (vCPNew == NULL) )
	{
		printf("Error: SeqClust_KMVnorm, cannot create memory for prior probabilities of clusters\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nK; ni++)
	{
		vCP[ni] = 1.0/(double)nK;
		vCPNew[ni] = vCP[ni];
	}

	vCMean = (double **)calloc(nK, sizeof(double *));
	vCMeanNew = (double **)calloc(nK, sizeof(double *));
	vCVar = (struct DOUBLEMATRIX ***)calloc(nK, sizeof(struct DOUBLEMATRIX **));
	vCVarNew = (struct DOUBLEMATRIX ***)calloc(nK, sizeof(struct DOUBLEMATRIX **));
	if( (vCMean == NULL) || (vCMeanNew == NULL) || (vCVar == NULL) || (vCVarNew == NULL) )
	{
		printf("Error: SeqClust_KMVnorm, cannot create memory for cluster mean and variance\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nK; ni++)
	{
		vCMean[ni] = (double *)calloc(nSampleNum, sizeof(double));
		vCMeanNew[ni] = (double *)calloc(nSampleNum, sizeof(double));
		vCVar[ni] = (struct DOUBLEMATRIX **)calloc(nGroupNum, sizeof(struct DOUBLEMATRIX *));
		vCVarNew[ni] = (struct DOUBLEMATRIX **)calloc(nGroupNum, sizeof(struct DOUBLEMATRIX *));
		if( (vCMean[ni] == NULL) || (vCMeanNew[ni] == NULL) || (vCVar[ni] == NULL) || (vCVarNew[ni] == NULL) )
		{
			printf("Error: SeqClust_KMVnorm, cannot create memory for cluster mean and variance\n");
			exit(EXIT_FAILURE);
		}
		for(nj=0; nj<nGroupNum; nj++)
		{
			vCVar[ni][nj] = CreateDoubleMatrix(vGroupSize[nj], vGroupSize[nj]);
			if( vCVar[ni][nj] == NULL )
			{
				printf("Error: SeqClust_KMVnorm, cannot create memory for cluster variance\n");
				exit(EXIT_FAILURE);
			}
			for(nk=0; nk<vGroupSize[nj]; nk++)
			{
				DMSETAT(vCVar[ni][nj], nk, nk, 1.0);
			}
		}

		dRand = rand_u();
		nk = (int)(dRand*nGeneNum);
		if(nk >= nGeneNum)
			nk = nGeneNum-1;
		if(nk < 0)
			nk = 0;
		for(nj=0; nj<nSampleNum; nj++)
		{
			vCMean[ni][nj] = vData[nk][nj];
			vCMeanNew[ni][nj] = vData[nk][nj];
		}
	}

	vLike = (double **)calloc(nGeneNum, sizeof(double *));
	if( vLike == NULL )
	{
		printf("Error: SeqClust_KMVnorm, cannot create memory for computing likelihood\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nGeneNum; ni++)
	{
		vLike[ni] = (double *)calloc(nK, sizeof(double));
		if( vLike[ni] == NULL )
		{
			printf("Error: SeqClust_KMVnorm, cannot create memory for computing likelihood\n");
			exit(EXIT_FAILURE);
		}
	}

	vTemp1 = NULL;
	vTemp1 = (double *)calloc(nK, sizeof(double));
	if(vTemp1 == NULL)
	{
		printf("Error: SeqClust_KMVnorm, cannot create memory for intermediate results\n");
		exit(EXIT_FAILURE);
	}

	/* EM */
	for(nIter=0; nIter<nMaxIter; nIter++)
	{
		/* prepare parameter */
		for(nk=0; nk<nK; nk++)
			vCPNew[nk] = log(vCP[nk]);

		/* decomposition */
		for(nk=0; nk<nK; nk++)
		{
			vTemp1[nk] = dnormfactor;

			for(nj=0; nj<nGroupNum; nj++)
			{
				vCVarNew[nk][nj] = NULL;

				/* find determinant */
				pD = NULL;
				pV = NULL;
				DMSYMEIGEN(vCVar[nk][nj], &pD, &pV);
				for(ni=0; ni<vGroupSize[nj]; ni++)
				{
					if(pD->pMatElement[ni] < 1e-6)
					{
						printf("Warning: SeqClust_KMVnorm, covariance matrix nearly singular at K = %d! Skip ...\n", nK);
						nSingularExit = 1;
						DestroyDoubleMatrix(pD);
						DestroyDoubleMatrix(pV);
						break;
					}

					vTemp1[nk] -= 0.5*log(pD->pMatElement[ni]);
				}

				if(nSingularExit == 1)
					break;
				
				/* find inverse */
				DMSYMINVEIGEN(pD, pV, (vCVarNew[nk]+nj));

				/* destroy matrix */
				DestroyDoubleMatrix(pD);
				DestroyDoubleMatrix(pV);
			}

			if(nSingularExit == 1)
				break;
		}

		if(nSingularExit == 1)
		{
			for(ni=0; ni<nk; ni++)
			{
				for(nx=0; nx<nGroupNum; nx++)
				{
					DestroyDoubleMatrix(vCVarNew[ni][nx]);
					vCVarNew[ni][nx] = NULL;
				}
			}
			for(nx=0; nx<nj; nx++)
			{
				DestroyDoubleMatrix(vCVarNew[nk][nx]);
				vCVarNew[nk][nx] = NULL;
			}
			break;
		}

		/* compute likelihood & posterior */
		dLikeNew = 0.0;
		for(ni=0; ni<nGeneNum; ni++)
		{
			for(nk=0; nk<nK; nk++)
			{
				vLike[ni][nk] = vCPNew[nk]+vTemp1[nk];
				for(nj=0; nj<nGroupNum; nj++)
				{
					vLike[ni][nk] += SeqClust_KMVnorm_Quadrature(vData[ni], vCMean[nk], vCVarNew[nk][nj], vGroupId[nj], vGroupSize[nj]);
				}

				if(nk == 0)
					dLMax = vLike[ni][nk];
				else if(vLike[ni][nk] > dLMax)
					dLMax = vLike[ni][nk];
			}

			dLSum = 0.0;
			for(nk=0; nk<nK; nk++)
			{
				vLike[ni][nk] = exp(vLike[ni][nk]-dLMax);
				dLSum += vLike[ni][nk];
			}

			dLikeNew += dLMax+log(dLSum);

			for(nk=0; nk<nK; nk++)
			{
				vLike[ni][nk] = vLike[ni][nk]/dLSum;
			}
		}

		if(nFinal == 1)
		{
			dLike = dLikeNew;
			for(nk=0; nk<nK; nk++)
			{
				for(nj=0; nj<nGroupNum; nj++)
				{
					DestroyDoubleMatrix(vCVarNew[nk][nj]);
					vCVarNew[nk][nj] = NULL;
				}
			}
			break;
		}

		/* update parameters */
		for(nk=0; nk<nK; nk++)
		{
			vCPNew[nk] = 0.0;
			for(nj=0; nj<nSampleNum; nj++)
			{
				vCMeanNew[nk][nj] = 0.0;
			}

			for(nj=0; nj<nGroupNum; nj++)
			{			
				DestroyDoubleMatrix(vCVarNew[nk][nj]);
				vCVarNew[nk][nj] = NULL;
				vCVarNew[nk][nj] = CreateDoubleMatrix(vGroupSize[nj], vGroupSize[nj]);
				if(vCVarNew[nk][nj] == NULL)
				{
					printf("Error: SeqClust_KMVnorm, cannot create memory for new covariance matrix\n");
					exit(EXIT_FAILURE);
				}
			}
		}
		for(ni=0; ni<nGeneNum; ni++)
		{
			for(nk=0; nk<nK; nk++)
			{
				vCPNew[nk] += vLike[ni][nk];
				for(nj=0; nj<nSampleNum; nj++)
				{
					vCMeanNew[nk][nj] += vLike[ni][nk]*vData[ni][nj];
				}

				for(nj=0; nj<nGroupNum; nj++)
				{
					pinvV = vCVarNew[nk][nj]->pMatElement;
					for(nx=0; nx<vGroupSize[nj]; nx++)
					{
						nu = vGroupId[nj][nx];
						for(ny=0; ny<vGroupSize[nj]; ny++)
						{
							nw = vGroupId[nj][ny];
							(*pinvV) += vLike[ni][nk]*vData[ni][nu]*vData[ni][nw];
							pinvV++;
						}
					}
				}
			}
		}

		dErr = 0.0;
		for(nk=0; nk<nK; nk++)
		{
			for(nj=0; nj<nSampleNum; nj++)
			{
				vCMeanNew[nk][nj]  /= vCPNew[nk];
				dRErr = fabs( (vCMeanNew[nk][nj]-vCMean[nk][nj])/(vCMean[nk][nj]+1e-6));
				if(dRErr > dErr)
					dErr = dRErr;
			}

			for(nj=0; nj<nGroupNum; nj++)
			{
				pinvV = vCVarNew[nk][nj]->pMatElement;
				pinvV2 = vCVar[nk][nj]->pMatElement;
				for(nx=0; nx<vGroupSize[nj]; nx++)
				{
					nu = vGroupId[nj][nx];
					for(ny=0; ny<vGroupSize[nj]; ny++)
					{
						nw = vGroupId[nj][ny];
						(*pinvV) = (*pinvV)/vCPNew[nk]-vCMeanNew[nk][nu]*vCMeanNew[nk][nw];

						dRErr = fabs( ((*pinvV)-(*pinvV2))/((*pinvV2)+1e-6) );
						if(dRErr > dErr)
							dErr = dRErr;

						pinvV++;
						pinvV2++;
					}
				}
			}

			vCPNew[nk] /= nGeneNum;
			dRErr = fabs( (vCPNew[nk]-vCP[nk])/(vCP[nk]+1e-6) );
			if(dRErr > dErr)
				dErr = dRErr;
		}

		/* evaluate stopping criteria */
		if(nIter > 0)
		{
			if(dLikeNew < dLike)
			{
				if( fabs(dLikeNew-dLike)/(dLike+1e-6) > 1e-6)
					printf("Warning: decreasing likelihood!\n");
				dLike = dLikeNew;
				for(nk=0; nk<nK; nk++)
				{
					for(nj=0; nj<nGroupNum; nj++)
					{
						DestroyDoubleMatrix(vCVarNew[nk][nj]);
						vCVarNew[nk][nj] = NULL;
					}
				}
				break;
			}
		}

		if(dErr < dTol)
		{
			nFinal = 1;
			nIter--;
		}

		for(nk=0; nk<nK; nk++)
		{
			for(nj=0; nj<nSampleNum; nj++)
			{
				vCMean[nk][nj] = vCMeanNew[nk][nj];
			}
			for(nj=0; nj<nGroupNum; nj++)
			{
				DestroyDoubleMatrix(vCVar[nk][nj]);
				vCVar[nk][nj]  = vCVarNew[nk][nj];
				vCVarNew[nk][nj] = NULL;
			}
			vCP[nk] = vCPNew[nk];
		}
		dLike = dLikeNew;
	}

	/* Export results */
	sprintf(strFileName, "%s%s.pp", strOutputPath, strOutputFile);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: SeqClust_KMVnorm, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nGeneNum; ni++)
	{
		fprintf(fpOut, "%d", ni);
		for(nk=0; nk<nK; nk++)
		{
			fprintf(fpOut, "\t%f", vLike[ni][nk]); 
		}
		fprintf(fpOut, "\n");
	}
	fclose(fpOut);

	sprintf(strFileName, "%s%s.prior", strOutputPath, strOutputFile);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: SeqClust_KMVnorm, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	for(nk=0; nk<nK; nk++)
	{
		if(nk != 0)
			fprintf(fpOut, "\t");

		fprintf(fpOut, "%f", vCP[nk]); 
	}
	fprintf(fpOut, "\n");
	fclose(fpOut);

	sprintf(strFileName, "%s%s.mu", strOutputPath, strOutputFile);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: SeqClust_KMVnorm, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	for(nk=0; nk<nK; nk++)
	{
		for(nj=0; nj<nSampleNum; nj++)
		{
			if(nj != 0)
				fprintf(fpOut, "\t");
			fprintf(fpOut, "%f", vCMean[nk][nj]); 
		}
		fprintf(fpOut, "\n");
	}
	fclose(fpOut);

	sprintf(strFileName, "%s%s.va", strOutputPath, strOutputFile);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: SeqClust_KMVnorm, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	for(nk=0; nk<nK; nk++)
	{
		fprintf(fpOut, "==cluster_%d==\n\n", nk);
		for(nj=0; nj<nGroupNum; nj++)
		{
			fprintf(fpOut, ">group_%d\n", nj);
			ny = 0;
			for(ni=0; ni<vGroupSize[nj]; ni++)
			{
				for(nx=0; nx<vGroupSize[nj]; nx++)
				{
					if(nx != 0)
						fprintf(fpOut, "\t");
					fprintf(fpOut, "%f", vCVar[nk][nj]->pMatElement[ny]);
					ny++;
				}
				fprintf(fpOut, "\n");
			}
			fprintf(fpOut, "\n"); 
		}
	}
	fclose(fpOut);

	/* Release memory */
	for(ni=0; ni<nGeneNum; ni++)
	{
		free(vLike[ni]);
	}
	free(vLike);
	free(vCP);
	free(vCPNew);
	for(ni=0; ni<nK; ni++)
	{
		free(vCMean[ni]);
		free(vCMeanNew[ni]);
		for(nj=0; nj<nGroupNum; nj++)
		{
			DestroyDoubleMatrix(vCVar[ni][nj]);
			vCVar[ni][nj] = NULL;
		}
		free(vCVar[ni]);
		free(vCVarNew[ni]);
	}
	free(vCMean);
	free(vCMeanNew);
	free(vCVar);
	free(vCVarNew);
	free(vTemp1);

	/* return */
	if(nSingularExit == 1)
		dBIC = 1e30;
	else
	{
		nParamNum = 0;
		for(nj=0; nj<nGroupNum; nj++)
		{
			nParamNum += (vGroupSize[nj]+1)*vGroupSize[nj]/2;
		}
		nParamNum = nK-1+nK*(nSampleNum+nParamNum);
		dBIC = -2.0*dLike+nParamNum*log((double)(nGeneNum));
	}
	printf(" iteration number = %d;  relative err = %f; logLike = %f; parameter number = %d; BIC = %f \n", nIter, dErr, dLike, nParamNum, dBIC);

	return dBIC;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_KMVnorm_Quadrature()                                          */
/*  compute normal quadrature -(x-mu)'invV(x-mu)/2                         */
/* ----------------------------------------------------------------------- */ 
double SeqClust_KMVnorm_Quadrature(double *vX, double *vMu, 
					struct DOUBLEMATRIX *pinvV, int *vGroupId, int nGroupSize)
{
	/* define */
	double dR = 0.0;
	double *vTemp1,*vTemp2;
	int ni,nj,nk;
	double *pE;

	/* compute */
	vTemp1 = NULL;
	vTemp2 = NULL;
	vTemp1 = (double *)calloc(nGroupSize, sizeof(double));
	vTemp2 = (double *)calloc(nGroupSize, sizeof(double));

	for(ni=0; ni<nGroupSize; ni++)
	{
		nk = vGroupId[ni];
		vTemp1[ni] = vX[nk]-vMu[nk];
		vTemp2[ni] = 0.0;
	}

	pE = pinvV->pMatElement;
	for(ni=0; ni<nGroupSize; ni++)
	{
		for(nj=0; nj<nGroupSize; nj++)
		{
			vTemp2[ni] += (*pE)*vTemp1[nj];
			pE++;
		}
	}

	for(ni=0; ni<nGroupSize; ni++)
		dR += vTemp1[ni]*vTemp2[ni];

	dR = -dR/2.0;

	/* release memory */
	free(vTemp1);
	free(vTemp2);

	/* return */
	return dR;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_KMFnorm_Quadrature()                                          */
/*  compute normal quadrature -(x-mu)'invV(x-mu)/2                         */
/* ----------------------------------------------------------------------- */ 
double SeqClust_KMFnorm_Quadrature(double *vX, double *vMu, struct DOUBLEMATRIX *pinvV)
{
	/* define */
	double dR = 0.0;
	double *vTemp1,*vTemp2;
	int ni,nj;
	double *pE;

	/* compute */
	vTemp1 = NULL;
	vTemp2 = NULL;
	vTemp1 = (double *)calloc(pinvV->nHeight, sizeof(double));
	vTemp2 = (double *)calloc(pinvV->nHeight, sizeof(double));

	for(ni=0; ni<pinvV->nHeight; ni++)
	{
		vTemp1[ni] = vX[ni]-vMu[ni];
		vTemp2[ni] = 0.0;
	}

	pE = pinvV->pMatElement;
	for(ni=0; ni<pinvV->nHeight; ni++)
	{
		for(nj=0; nj<pinvV->nHeight; nj++)
		{
			vTemp2[ni] += (*pE)*vTemp1[nj];
			pE++;
		}
	}

	for(ni=0; ni<pinvV->nHeight; ni++)
		dR += vTemp1[ni]*vTemp2[ni];

	dR = -dR/2.0;

	/* release memory */
	free(vTemp1);
	free(vTemp2);

	/* return */
	return dR;
}


/* ----------------------------------------------------------------------- */ 
/*  SeqClust_KMFnorm()                                                     */
/*  Mixture factor model based clustering (K clusters with covariance      */
/*  matrix modeled by nF factors).                                         */
/* ----------------------------------------------------------------------- */ 
double SeqClust_KMFnorm(int nGeneNum, int nSampleNum, double **vData, 
				int nF, int nK, int nMaxIter, double dTol, 
				char strOutputPath[], char strOutputFile[])
{
	/* define */
	double dBIC = 0.0;

	double *vCP = NULL;
	double *vCPNew = NULL;

	double **vCMean = NULL;
	double **vCMeanNew = NULL;
	
	struct DOUBLEMATRIX **vCVar = NULL;
	struct DOUBLEMATRIX **vCVarInv = NULL;
	
	struct DOUBLEMATRIX **vCFac = NULL;
	struct DOUBLEMATRIX **vCFacNew = NULL;

	double **vCPhi = NULL;
	double **vCPhiNew = NULL;
	
	double **vEZ = NULL;

	double *vTemp1;

	struct DOUBLEMATRIX *pD;
	struct DOUBLEMATRIX *pV;
	
	int nSingularExit = 0;

	double dTemp;
	double dRand;
	double dpi = 4.0*atan(1.0);
	double dnormfactor;
	int ni,nj,nk,nl,nx,ny;
	double *pEle,*pEle2,*pEle3;
	
	double **vLike = NULL;
	double dLike = 0.0;
	double dLikeNew = 0.0;
	double dLMax,dLSum,dRErr;

	struct DOUBLEMATRIX *pBeta;
	struct DOUBLEMATRIX *pZZ0;

	struct DOUBLEMATRIX *pZZSum;
	struct DOUBLEMATRIX *pZSum;

	struct DOUBLEMATRIX *pTempMat1;
	struct DOUBLEMATRIX *pTempMat2;
	
	double dErr = 1e6;
	int nIter = 0;
	int nFinal = 0;

	int nParamNum = 0;
	FILE *fpOut;
	char strFileName[MED_LINE_LENGTH];

	/* Initialize */
	dnormfactor = -nSampleNum*log(2*dpi)/2.0;

	/* prior abundance */
	vCP = (double *)calloc(nK, sizeof(double));
	vCPNew = (double *)calloc(nK, sizeof(double));
	if( (vCP == NULL) || (vCPNew == NULL) )
	{
		printf("Error: SeqClust_KMFnorm, cannot create memory for prior probabilities of clusters\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nK; ni++)
	{
		vCP[ni] = 1.0/(double)nK;
		vCPNew[ni] = vCP[ni];
	}

	/* prior mean */
	vCMean = (double **)calloc(nK, sizeof(double *));
	vCMeanNew = (double **)calloc(nK, sizeof(double *));
	/* prior variance */
	vCVar = (struct DOUBLEMATRIX **)calloc(nK, sizeof(struct DOUBLEMATRIX *));
	vCVarInv = (struct DOUBLEMATRIX **)calloc(nK, sizeof(struct DOUBLEMATRIX *));
	/* prior factor */
	vCFac = (struct DOUBLEMATRIX **)calloc(nK, sizeof(struct DOUBLEMATRIX *));
	vCFacNew = (struct DOUBLEMATRIX **)calloc(nK, sizeof(struct DOUBLEMATRIX *));
	/* prior residual var */
	vCPhi = (double **)calloc(nK, sizeof(double *));
	vCPhiNew = (double **)calloc(nK, sizeof(double *));

	if( (vCMean == NULL) || (vCMeanNew == NULL) || (vCVar == NULL) || (vCVarInv == NULL)
		|| (vCFac == NULL) || (vCFacNew == NULL) || (vCPhi == NULL) || (vCPhiNew == NULL) )
	{
		printf("Error: SeqClust_KMFnorm, cannot create memory for cluster mean and variance\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nK; ni++)
	{
		vCMean[ni] = (double *)calloc(nSampleNum, sizeof(double));
		vCMeanNew[ni] = (double *)calloc(nSampleNum, sizeof(double));
		vCVar[ni] = NULL;
		vCVarInv[ni] = NULL;
		vCFac[ni] = CreateDoubleMatrix(nSampleNum, nF);
		vCFacNew[ni] = NULL;
		vCPhi[ni] = (double *)calloc(nSampleNum, sizeof(double));
		vCPhiNew[ni] = (double *)calloc(nSampleNum, sizeof(double));
		if( (vCMean[ni] == NULL) || (vCMeanNew[ni] == NULL) || 
			(vCFac[ni] == NULL) || (vCPhi[ni] == NULL) || (vCPhiNew[ni] == NULL) )
		{
			printf("Error: SeqClust_KMFnorm, cannot create memory for cluster mean and variance\n");
			exit(EXIT_FAILURE);
		}

		for(nj=0; nj<nSampleNum; nj++)
		{
			vCPhi[ni][nj] = 1.0;
			vCPhiNew[ni][nj] = 1.0;
		}

		pEle = vCFac[ni]->pMatElement;
		for(nj=0; nj<nSampleNum; nj++)
		{
			for(nk=0; nk<nF; nk++)
			{
				*pEle = 0.01*normrnd(0,1)/sqrt((double)nF);
				pEle++;
			}
		}

		vCFacNew[ni] = NULL;
		vCFacNew[ni] = DM_T(vCFac[ni]);
		if(vCFacNew[ni] == NULL)
		{
			printf("Error: SeqClust_KMVnorm, cannot transpose the factor matrix!\n");
			exit(EXIT_FAILURE);
		}

		vCVar[ni] = DMMUL(vCFac[ni], vCFacNew[ni]);
		if(vCVar[ni] == NULL)
		{
			printf("Error: SeqClust_KMVnorm, cannot initialize the variance matrix!\n");
			exit(EXIT_FAILURE);
		}
		for(nj=0; nj<nSampleNum; nj++)
		{
			dTemp = DMGETAT(vCVar[ni], nj, nj)+vCPhi[ni][nj];
			DMSETAT(vCVar[ni], nj, nj, dTemp);
		}

		dRand = rand_u();
		nk = (int)(dRand*nGeneNum);
		if(nk >= nGeneNum)
			nk = nGeneNum-1;
		if(nk < 0)
			nk = 0;
		for(nj=0; nj<nSampleNum; nj++)
		{
			vCMean[ni][nj] = vData[nk][nj];
			vCMeanNew[ni][nj] = vData[nk][nj];
		}
	}

	/* posterior probability */
	/* stroage for EZ, reused by different clusters */
	vLike = (double **)calloc(nGeneNum, sizeof(double *));
	vEZ = (double **)calloc(nGeneNum, sizeof(double *));
	if( (vLike == NULL) || (vEZ == NULL) )
	{
		printf("Error: SeqClust_KMFnorm, cannot create memory for computing likelihood or EZ.\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nGeneNum; ni++)
	{
		vLike[ni] = (double *)calloc(nK, sizeof(double));
		if( vLike[ni] == NULL )
		{
			printf("Error: SeqClust_KMFnorm, cannot create memory for computing likelihood\n");
			exit(EXIT_FAILURE);
		}

		vEZ[ni] = (double *)calloc(nF, sizeof(double));
		if(vEZ[ni] == NULL)
		{
			printf("Error: SeqClust_KMFnorm, cannot create memory for EZ.\n");
			exit(EXIT_FAILURE);
		}
	}

	/* auxiliary variables */

	/* for computing likelihood */
	vTemp1 = NULL;
	vTemp1 = (double *)calloc(nK, sizeof(double));
	if(vTemp1 == NULL)
	{
		printf("Error: SeqClust_KMFnorm, cannot create memory for intermediate results\n");
		exit(EXIT_FAILURE);
	}

	/* EM */
	for(nIter=0; nIter<nMaxIter; nIter++)
	{
		/* prepare parameter */
		for(nk=0; nk<nK; nk++)
		{
			/* prior & normal factors */
			vCPNew[nk] = log(vCP[nk]);
			vTemp1[nk] = dnormfactor;

			/* find determinant */
			pD = NULL;
			pV = NULL;
			DMSYMEIGEN(vCVar[nk], &pD, &pV);
			for(ni=0; ni<nSampleNum; ni++)
			{
				if(pD->pMatElement[ni] < 1e-6)
				{
					printf("Warning: SeqClust_KMFnorm, covariance matrix nearly singular at K = %d! Skip ...\n", nK);
					nSingularExit = 1;
					DestroyDoubleMatrix(pD);
					DestroyDoubleMatrix(pV);
					break;
				}

				vTemp1[nk] -= 0.5*log(pD->pMatElement[ni]);
			}

			if(nSingularExit == 1)
				break;
				
			/* find inverse */
			vCVarInv[nk] = NULL;
			DMSYMINVEIGEN(pD, pV, (vCVarInv+nk));
			if(vCVarInv[nk] == NULL)
			{
				printf("Error: SeqClust_KMFnorm, cannot compute inverse covariance!\n");
				exit(EXIT_FAILURE);
			}

			/* destroy matrix */
			DestroyDoubleMatrix(pD);
			DestroyDoubleMatrix(pV);
		}

		if(nSingularExit == 1)
		{
			for(ni=0; ni<nk; ni++)
			{
				DestroyDoubleMatrix(vCVarInv[nk]);
				vCVarInv[nk] = NULL;
			}

			break;
		}

		/* compute likelihood & posterior */
		dLikeNew = 0.0;
		for(ni=0; ni<nGeneNum; ni++)
		{
			for(nk=0; nk<nK; nk++)
			{
				vLike[ni][nk] = vCPNew[nk]+vTemp1[nk]+SeqClust_KMFnorm_Quadrature(vData[ni], vCMean[nk], vCVarInv[nk]);

				if(nk == 0)
					dLMax = vLike[ni][nk];
				else if(vLike[ni][nk] > dLMax)
					dLMax = vLike[ni][nk];
			}

			dLSum = 0.0;
			for(nk=0; nk<nK; nk++)
			{
				vLike[ni][nk] = exp(vLike[ni][nk]-dLMax);
				dLSum += vLike[ni][nk];
			}

			dLikeNew += dLMax+log(dLSum);

			for(nk=0; nk<nK; nk++)
			{
				vLike[ni][nk] = vLike[ni][nk]/dLSum;
			}
		}

		if(nFinal == 1)
		{
			dLike = dLikeNew;
			for(nk=0; nk<nK; nk++)
			{
				DestroyDoubleMatrix(vCVarInv[nk]);
				vCVarInv[nk] = NULL;
			}
		
			break;
		}

		/* for each cluster, compute EZ, EZZ' and estimate Factor & Phi */
		dErr = 0.0;
		for(nk=0; nk<nK; nk++)
		{
			/* first compute Beta */
			pBeta = NULL;
			pBeta = DMMUL(vCFacNew[nk], vCVarInv[nk]);
			if(pBeta == NULL)
			{
				printf("Error: SeqClust_KMFnorm, cannot compute beta.\n");
				exit(EXIT_FAILURE);
			}

			/* then compute I-Beta*Factor */
			pZZ0 = NULL;
			pZZ0 = DMMUL(pBeta, vCFac[nk]);
			if(pZZ0 == NULL)
			{
				printf("Error: SeqClust_KMFnorm, cannot compute beta*factor.\n");
				exit(EXIT_FAILURE);
			}
			pEle = pZZ0->pMatElement;
			for(nj=0; nj<nF; nj++)
			{
				for(nl=0; nl<nF; nl++)
				{
					if(nl == nj)
						*pEle = 1.0-(*pEle);
					else
						*pEle = -(*pEle);

					pEle++;
				}
			}

			/* next prepare ZZSum and ZSum */
			pZZSum = NULL;
			pZZSum = CreateDoubleMatrix(nF+1, nF+1);
			pZSum = NULL;
			pZSum = CreateDoubleMatrix(nSampleNum, nF+1);
			if( (pZZSum == NULL) || (pZSum == NULL) )
			{
				printf("Error: SeqClust_KMFnorm, cannot create memory for computing Sum(Ew * EZZ') or Sum(Ew * EZ) \n");
				exit(EXIT_FAILURE);
			}

			/* scan genes to compute EZ, EZZ' */
			for(ni=0; ni<nGeneNum; ni++)
			{
				/* EZ */
				pEle = pBeta->pMatElement;
				for(nl=0; nl<nF; nl++)
				{
					vEZ[ni][nl] = 0.0;
					for(nj=0; nj<nSampleNum; nj++)
					{
						vEZ[ni][nl] += (*pEle)*(vData[ni][nj]-vCMean[nk][nj]);
						pEle++;
					}
				}

				/* EZSum */
				pEle = pZSum->pMatElement;
				for(nj=0; nj<nSampleNum; nj++)
				{
					for(nl=0; nl<nF; nl++)
					{
						(*pEle) += vLike[ni][nk]*vData[ni][nj]*vEZ[ni][nl];
						pEle++;
					}

					(*pEle) += vLike[ni][nk]*vData[ni][nj];
					pEle++;
				}

				/* EZZSum */
				pEle = pZZSum->pMatElement;
				pEle2 = pZZ0->pMatElement;
				for(nj=0; nj<nF; nj++)
				{
					for(nl=0; nl<nF; nl++)
					{
						(*pEle) += vLike[ni][nk]*((*pEle2)+vEZ[ni][nj]*vEZ[ni][nl]);

						pEle++;
						pEle2++;
					}

					(*pEle) += vLike[ni][nk]*vEZ[ni][nj];
					pEle++;
				}
				for(nl=0; nl<nF; nl++)
				{
					(*pEle) += vLike[ni][nk]*vEZ[ni][nl];
					pEle++;
				}
				(*pEle) += vLike[ni][nk];
			}

			/* compute new mu, new factor */
			pD = NULL;
			pV = NULL;
			DMSYMEIGEN(pZZSum, &pD, &pV);
			for(ni=0; ni<=nF; ni++)
			{
				if(pD->pMatElement[ni] < 1e-6)
				{
					printf("Warning: SeqClust_KMFnorm, EZZSum matrix nearly singular!\n", nK);
				}
			}

			/* find inverse */
			DestroyDoubleMatrix(vCFacNew[nk]);
			vCFacNew[nk] = NULL;
			pTempMat1 = NULL;
			pTempMat2 = NULL;
			DMSYMINVEIGEN(pD, pV, &pTempMat1);
			if(pTempMat1 == NULL)
			{
				printf("Error: SeqClust_KMFnorm, cannot compute inverse of Sum(Ew * EZZ)!\n");
				exit(EXIT_FAILURE);
			}
			pTempMat2 = DMMUL(pZSum, pTempMat1);
			if(pTempMat2 == NULL)
			{
				printf("Error: SeqClust_KMFnorm, cannot compute new factor and mean!\n");
				exit(EXIT_FAILURE);
			}
			vCFacNew[nk] = CreateDoubleMatrix(nSampleNum, nF);
			if(vCFacNew[nk] == NULL)
			{
				printf("Error: SeqClust_KMFnorm, cannot allocate memory for new factor matrix!\n");
				exit(EXIT_FAILURE);
			}

			pEle = pTempMat2->pMatElement;
			pEle2 = vCFacNew[nk]->pMatElement;
			pEle3 = vCFac[nk]->pMatElement;
			for(nj=0; nj<nSampleNum; nj++)
			{
				for(nl=0; nl<nF; nl++)
				{
					(*pEle2) = (*pEle);
					
					dRErr = fabs( ((*pEle2)-(*pEle3))/((*pEle3)+1e-6) );
					if(dRErr > dErr)
						dErr = dRErr;

					pEle++;
					pEle2++;
					pEle3++;
				}

				vCMeanNew[nk][nj] = (*pEle);
				dRErr = fabs( (vCMeanNew[nk][nj]-vCMean[nk][nj])/(vCMean[nk][nj]+1e-6) );
				if(dRErr > dErr)
					dErr = dRErr;

				pEle++;
			}

			DestroyDoubleMatrix(pD);
			DestroyDoubleMatrix(pV);
			DestroyDoubleMatrix(pTempMat1);

			/* scan genes to compute Phi' */
			vCPNew[nk] = 0.0;
			for(nj=0; nj<nSampleNum; nj++)
			{
				vCPhiNew[nk][nj] = 0.0;
			}

			for(ni=0; ni<nGeneNum; ni++)
			{
				/* post p */
				vCPNew[nk] += vLike[ni][nk];

				/* Phi */
				pEle = pTempMat2->pMatElement;
				for(nj=0; nj<nSampleNum; nj++)
				{
					dTemp = 0.0;
					for(nl=0; nl<nF; nl++)
					{
						dTemp += (*pEle)*vEZ[ni][nl];
						pEle++;
					}
					dTemp += (*pEle);
					pEle++;

					vCPhiNew[nk][nj] += vLike[ni][nk]*(vData[ni][nj]-dTemp)*vData[ni][nj];
				}
			}

			for(nj=0; nj<nSampleNum; nj++)
			{
				vCPhiNew[nk][nj] /= vCPNew[nk];
				dRErr = fabs( (vCPhiNew[nk][nj]-vCPhi[nk][nj])/(vCPhi[nk][nj]+1e-6) );
				if(dRErr > dErr)
					dErr = dRErr;
			}

			DestroyDoubleMatrix(pTempMat2);

			/* compute prior */
			vCPNew[nk] /= nGeneNum;
			dRErr = fabs( (vCPNew[nk]-vCP[nk])/(vCP[nk]+1e-6) );
			if(dRErr > dErr)
				dErr = dRErr;

			/* clear memory */
			DestroyDoubleMatrix(pBeta);
			DestroyDoubleMatrix(pZZ0);
			DestroyDoubleMatrix(pZZSum);
			DestroyDoubleMatrix(pZSum);
		}

		
		/* evaluate stopping criteria */
		if(nIter > 0)
		{
			if(dLikeNew < dLike)
			{
				if( fabs(dLikeNew-dLike)/(dLike+1e-6) > 1e-6)
					printf("Warning: decreasing likelihood!\n");
				dLike = dLikeNew;
				
				break;
			}
		}

		if(dErr < dTol)
		{
			nFinal = 1;
			nIter--;
		}

		/* update parameter */
		for(nk=0; nk<nK; nk++)
		{
			vCP[nk] = vCPNew[nk];

			for(nj=0; nj<nSampleNum; nj++)
			{
				vCMean[nk][nj] = vCMeanNew[nk][nj];
				vCPhi[nk][nj] = vCPhiNew[nk][nj];
			}

			DestroyDoubleMatrix(vCFac[nk]);
			vCFac[nk]  = vCFacNew[nk];
			vCFacNew[nk] = NULL;

			DestroyDoubleMatrix(vCVar[nk]);
			DestroyDoubleMatrix(vCVarInv[nk]);
			vCVar[nk] = NULL;
			vCVarInv[nk] = NULL;

			vCFacNew[nk] = DM_T(vCFac[nk]);
			vCVar[nk] = DMMUL(vCFac[nk], vCFacNew[nk]);
			if(vCVar[nk] == NULL)
			{
				printf("Error: SeqClust_KMVnorm, cannot compute new variance matrix!\n");
				exit(EXIT_FAILURE);
			}
			for(nj=0; nj<nSampleNum; nj++)
			{
				dTemp = DMGETAT(vCVar[nk], nj, nj)+vCPhiNew[nk][nj];
				DMSETAT(vCVar[nk], nj, nj, dTemp);
			}			
		}
		dLike = dLikeNew;
	}

	/* Export results */
	sprintf(strFileName, "%s%s.pp", strOutputPath, strOutputFile);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: SeqClust_KMFnorm, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nGeneNum; ni++)
	{
		fprintf(fpOut, "%d", ni);
		for(nk=0; nk<nK; nk++)
		{
			fprintf(fpOut, "\t%f", vLike[ni][nk]); 
		}
		fprintf(fpOut, "\n");
	}
	fclose(fpOut);

	sprintf(strFileName, "%s%s.prior", strOutputPath, strOutputFile);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: SeqClust_KMFnorm, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	for(nk=0; nk<nK; nk++)
	{
		if(nk != 0)
			fprintf(fpOut, "\t");

		fprintf(fpOut, "%f", vCP[nk]); 
	}
	fprintf(fpOut, "\n");
	fclose(fpOut);

	sprintf(strFileName, "%s%s.mu", strOutputPath, strOutputFile);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: SeqClust_KMFnorm, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	for(nk=0; nk<nK; nk++)
	{
		for(nj=0; nj<nSampleNum; nj++)
		{
			if(nj != 0)
				fprintf(fpOut, "\t");
			fprintf(fpOut, "%f", vCMean[nk][nj]); 
		}
		fprintf(fpOut, "\n");
	}
	fclose(fpOut);

	sprintf(strFileName, "%s%s.va", strOutputPath, strOutputFile);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: SeqClust_KMFnorm, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	for(nk=0; nk<nK; nk++)
	{
		ny = 0;
		fprintf(fpOut, "==cluster_%d==\n\n", nk);
		for(ni=0; ni<nSampleNum; ni++)
		{
			for(nx=0; nx<nSampleNum; nx++)
			{
				if(nx != 0)
					fprintf(fpOut, "\t");
				fprintf(fpOut, "%f", vCVar[nk]->pMatElement[ny]);
				ny++;
			}
			fprintf(fpOut, "\n");
		}
		fprintf(fpOut, "\n"); 
	}
	fclose(fpOut);

	/* Release memory */
	for(ni=0; ni<nGeneNum; ni++)
	{
		free(vLike[ni]);
		free(vEZ[ni]);
	}
	free(vLike);
	free(vEZ);

	free(vCP);
	free(vCPNew);
	for(ni=0; ni<nK; ni++)
	{
		free(vCMean[ni]);
		free(vCMeanNew[ni]);
		DestroyDoubleMatrix(vCVar[ni]);
		vCVar[ni] = NULL;
		DestroyDoubleMatrix(vCFac[ni]);
		vCFac[ni] = NULL;
		DestroyDoubleMatrix(vCFacNew[ni]);
		vCFacNew[ni] = NULL;
		free(vCPhi[ni]);
		free(vCPhiNew[ni]);
	}
	free(vCMean);
	free(vCMeanNew);
	free(vCVar);
	free(vCVarInv);
	free(vCFac);
	free(vCFacNew);
	free(vCPhi);
	free(vCPhiNew);
	free(vTemp1);

	/* return */
	if(nSingularExit == 1)
		dBIC = 1e30;
	else
	{
		nParamNum = nK*nSampleNum*(nF+2)+nK-1;
		dBIC = -2.0*dLike+nParamNum*log((double)(nGeneNum));
	}
	printf(" iteration number = %d;  relative err = %f; logLike = %f; parameter number = %d; BIC = %f \n", nIter, dErr, dLike, nParamNum, dBIC);

	return dBIC;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_Kpolya()                                                      */
/*  Model based clustering (K polya distributions)                         */
/* ----------------------------------------------------------------------- */ 
double SeqClust_Kpolya(int nGeneNum, int nSampleNum, double **vData, 
				int nK, int nMaxIter, double dTol, 
				char strOutputPath[], char strOutputFile[])
{
	/* define */
	double dBIC = 0.0;

	double *vM = NULL;
	double *vMHist = NULL;
	double **vXHist = NULL;
	int nMMax = 0;
	int nXMax = 0;

	double *vCP = NULL;
	double **vCA = NULL;
	
	double *vCPNew = NULL;
	double **vCANew = NULL;
	
	double **vLike = NULL;
	double dLike = 0.0;
	double dLikeNew = 0.0;

	double dErr = 1e6;
	int nIter = 0;
	int nFinal = 0;

	int ni,nj,nk,nl;
	double dRand;

	double dLMax,dLSum,dRErr;

	FILE *fpOut;
	char strFileName[MED_LINE_LENGTH];

	/* Initialize */
	vM = (double *)calloc(nGeneNum, sizeof(double));
	if(vM == NULL)
	{
		printf("Error: SeqClust_Kpolya, cannot create memory for gene's total read count\n");
		exit(EXIT_FAILURE);
	}
	nMMax = 0;
	nXMax = 0;
	for(ni=0; ni<nGeneNum; ni++)
	{
		dLSum = 0.0;
		for(nj=0; nj<nSampleNum; nj++)
		{
			dLSum += vData[ni][nj];
			if((int)(vData[ni][nj]) > nXMax)
				nXMax = (int)(vData[ni][nj]);
		}
		if( (int)dLSum > nMMax)
			nMMax = (int)dLSum;
		vM[ni] = dLSum;
	}
	nMMax += 2;
	nXMax += 2;

	vMHist = (double *)calloc(nMMax, sizeof(double));
	if(vMHist == NULL)
	{
		printf("Error: SeqClust_Kpolya, cannot create histogram for gene's total read count\n");
		exit(EXIT_FAILURE);
	}
	vXHist = (double **)calloc(nXMax, sizeof(double *));
	if(vXHist == NULL)
	{
		printf("Error: SeqClust_Kpolya, cannot create histogram for gene's sample read count\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nXMax; ni++)
	{
		vXHist[ni] = (double *)calloc(nSampleNum, sizeof(double));
		if(vXHist[ni] == NULL)
		{
			printf("Error: SeqClust_Kpolya, cannot create histogram for gene's sample read count\n");
			exit(EXIT_FAILURE);
		}
	}

	vCP = (double *)calloc(nK, sizeof(double));
	vCPNew = (double *)calloc(nK, sizeof(double));
	if( (vCP == NULL) || (vCPNew == NULL) )
	{
		printf("Error: SeqClust_Kpolya, cannot create memory for prior probabilities of clusters\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nK; ni++)
	{
		vCP[ni] = 1.0/(double)nK;
		vCPNew[ni] = vCP[ni];
	}

	vCA = (double **)calloc(nK, sizeof(double *));
	vCANew = (double **)calloc(nK, sizeof(double *));
	if( (vCA == NULL) || (vCANew == NULL) )
	{
		printf("Error: SeqClust_Kpolya, cannot create memory for cluster mean and variance\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nK; ni++)
	{
		vCA[ni] = (double *)calloc(nSampleNum, sizeof(double));
		vCANew[ni] = (double *)calloc(nSampleNum, sizeof(double));
		if( (vCA[ni] == NULL) || (vCANew[ni] == NULL) )
		{
			printf("Error: SeqClust_Kpolya, cannot create memory for cluster mean and variance\n");
			exit(EXIT_FAILURE);
		}

		dRand = rand_u();
		nk = (int)(dRand*nGeneNum);
		if(nk >= nGeneNum)
			nk = nGeneNum-1;
		if(nk < 0)
			nk = 0;
		
		dLSum = 0.0;
		for(nj=0; nj<nSampleNum; nj++)
		{
			dLSum += vData[nk][nj]+1;
		}

		for(nj=0; nj<nSampleNum; nj++)
		{
			vCA[ni][nj] = nSampleNum*(vData[nk][nj]+1)/dLSum;
			vCANew[ni][nj] = vCA[ni][nj];
		}
	}

	vLike = (double **)calloc(nGeneNum, sizeof(double *));
	if( vLike == NULL )
	{
		printf("Error: SeqClust_Kpolya, cannot create memory for computing likelihood\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nGeneNum; ni++)
	{
		vLike[ni] = (double *)calloc(nK, sizeof(double));
		if( vLike[ni] == NULL )
		{
			printf("Error: SeqClust_Kpolya, cannot create memory for computing likelihood\n");
			exit(EXIT_FAILURE);
		}
	}

	/* EM */
	for(nIter=0; nIter<nMaxIter; nIter++)
	{
		/* prepare parameter */
		for(nk=0; nk<nK; nk++)
		{
			vCPNew[nk] = log(vCP[nk]);
			dLSum = 0.0;
			for(nj=0; nj<nSampleNum; nj++)
			{
				dLSum += vCA[nk][nj];
				vCPNew[nk] -= gammaln(vCA[nk][nj]);
			}
			vCPNew[nk] = vCPNew[nk]+gammaln(dLSum);
			vCANew[nk][0] = dLSum;
		}

		/* compute likelihood & posterior */
		dLikeNew = 0.0;
		for(ni=0; ni<nGeneNum; ni++)
		{
			for(nk=0; nk<nK; nk++)
			{
				vLike[ni][nk] = vCPNew[nk];
				for(nj=0; nj<nSampleNum; nj++)
				{
					vLike[ni][nk] += gammaln(vData[ni][nj]+vCA[nk][nj]);
				}
				vLike[ni][nk] -= gammaln(vCANew[nk][0]+vM[ni]);

				if(nk == 0)
					dLMax = vLike[ni][nk];
				else if(vLike[ni][nk] > dLMax)
					dLMax = vLike[ni][nk];
			}

			dLSum = 0.0;
			for(nk=0; nk<nK; nk++)
			{
				vLike[ni][nk] = exp(vLike[ni][nk]-dLMax);
				dLSum += vLike[ni][nk];
			}

			dLikeNew += dLMax+log(dLSum);

			for(nk=0; nk<nK; nk++)
			{
				vLike[ni][nk] = vLike[ni][nk]/dLSum;
			}
		}

		if(nFinal == 1)
		{
			dLike = dLikeNew;
			break;
		}

		/* update parameters */
		for(nk=0; nk<nK; nk++)
		{
			vCPNew[nk] = 0.0;
			for(nj=0; nj<nSampleNum; nj++)
			{
				vCANew[nk][nj] = 0.0;
			}
		}

		/* update alpha */
		for(nk=0; nk<nK; nk++)
		{
			for(ni=0; ni<nMMax; ni++)
				vMHist[ni] = 0.0;
			for(ni=0; ni<nXMax; ni++)
				for(nj=0; nj<nSampleNum; nj++)
					vXHist[ni][nj] = 0.0;

			for(ni=0; ni<nGeneNum; ni++)
			{
				for(nl=0; nl<=(vM[ni]-0.9); nl++)
				{
					vMHist[nl] += vLike[ni][nk];
				}

				for(nj=0; nj<nSampleNum; nj++)
				{
					for(nl=0; nl<=(vData[ni][nj]-0.9); nl++)
					{
						vXHist[nl][nj] += vLike[ni][nk];
					}
				}
			}

			SeqClust_Kpolya_MM(nSampleNum, nMMax, nXMax, vMHist, vXHist, vCA[nk], vCANew[nk]);
		}

		/* update p */
		for(ni=0; ni<nGeneNum; ni++)
		{
			for(nk=0; nk<nK; nk++)
			{
				vCPNew[nk] += vLike[ni][nk];
			}
		}

		/* compute error */
		dErr = 0.0;
		for(nk=0; nk<nK; nk++)
		{
			vCPNew[nk] /= nGeneNum;
			dRErr = fabs( (vCPNew[nk]-vCP[nk])/(vCP[nk]+1e-6) );
			if(dRErr > dErr)
				dErr = dRErr;

			for(nj=0; nj<nSampleNum; nj++)
			{
				dRErr = fabs( (vCANew[nk][nj]-vCA[nk][nj])/(vCA[nk][nj]+1e-6) );
				if(dRErr > dErr)
					dErr = dRErr;
			}
		}

				
		/* evaluate stopping criteria */
		if(nIter > 0)
		{
			if(dLikeNew < dLike)
			{
				if( fabs(dLikeNew-dLike)/(dLike+1e-6) > 1e-6)
					printf("Warning: decreasing likelihood!\n");
				dLike = dLikeNew;
				break;
			}
		}

		if(dErr < dTol)
		{
			nFinal = 1;
			nIter--;
		}

		for(nk=0; nk<nK; nk++)
		{
			for(nj=0; nj<nSampleNum; nj++)
			{
				vCA[nk][nj] = vCANew[nk][nj];
			}
			vCP[nk] = vCPNew[nk];
			dLike = dLikeNew;
		}
	}

	/* Export results */
	sprintf(strFileName, "%s%s.pp", strOutputPath, strOutputFile);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: SeqClust_Kpolya, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nGeneNum; ni++)
	{
		fprintf(fpOut, "%d", ni);
		for(nk=0; nk<nK; nk++)
		{
			fprintf(fpOut, "\t%f", vLike[ni][nk]); 
		}
		fprintf(fpOut, "\n");
	}
	fclose(fpOut);

	sprintf(strFileName, "%s%s.prior", strOutputPath, strOutputFile);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: SeqClust_Kpolya, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	for(nk=0; nk<nK; nk++)
	{
		if(nk != 0)
			fprintf(fpOut, "\t");

		fprintf(fpOut, "%f", vCP[nk]); 
	}
	fprintf(fpOut, "\n");
	fclose(fpOut);

	sprintf(strFileName, "%s%s.alpha", strOutputPath, strOutputFile);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: SeqClust_Kpolya, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	for(nk=0; nk<nK; nk++)
	{
		for(nj=0; nj<nSampleNum; nj++)
		{
			if(nj != 0)
				fprintf(fpOut, "\t");
			fprintf(fpOut, "%f", vCA[nk][nj]); 
		}
		fprintf(fpOut, "\n");
	}
	fclose(fpOut);

	/* Release memory */
	free(vM);
	free(vMHist);
	for(ni=0; ni<nXMax; ni++)
	{
		free(vXHist[ni]);
	}
	free(vXHist);

	for(ni=0; ni<nGeneNum; ni++)
	{
		free(vLike[ni]);
	}
	free(vLike);
	free(vCP);
	free(vCPNew);
	for(ni=0; ni<nK; ni++)
	{
		free(vCA[ni]);
		free(vCANew[ni]);
	}
	free(vCA);
	free(vCANew);
	
	/* return */
	dBIC = -2.0*dLike+(nK-1+nK*nSampleNum)*log((double)(nGeneNum));
	printf(" iteration number = %d;  relative err = %f; BIC = %f \n", nIter, dErr, dBIC);

	return dBIC;
}


/* ----------------------------------------------------------------------- */ 
/*  SeqClust_Kpolya_MM()                                                   */
/*  MM algorithm for fitting polya distribution.                           */
/* ----------------------------------------------------------------------- */ 
int SeqClust_Kpolya_MM(int nSampleNum, int nMMax, int nXMax, 
					   double *vMHist, double **vXHist, 
					   double *vAold, double *vAnew)
{
	/* define */
	double *vA0 = NULL;
	double *vAmid1 = NULL;
	double *vAmid2 = NULL;
	double *vA1 = NULL;
	int nMaxIter = 1000;
	int ni,nj,nk;
	double dErr;
	double dS1,dS2;
	double dASum;
	double dI1,dI2;
	double dLike0,dLike1;
	int nAPositive;

	/* initialize */
	vA0 = (double *)calloc(nSampleNum, sizeof(double));
	vAmid1 = (double *)calloc(nSampleNum, sizeof(double));
	vAmid2 = (double *)calloc(nSampleNum, sizeof(double));
	vA1 = (double *)calloc(nSampleNum, sizeof(double));
	if( (vA0 == NULL) || (vAmid1 == NULL) || (vAmid2 == NULL) || (vA1 == NULL) )
	{
		printf("Error: SeqClust_Kpolya_MM, cannot create memory for updating alpha!\n");
		exit(EXIT_FAILURE);
	}
	for(nj=0; nj<nSampleNum; nj++)
		vA0[nj] = vAold[nj];

	dLike0 = SeqClust_Kpolya_MM_LogLike(nSampleNum, nMMax, nXMax, vMHist, vXHist, vA0);
	
	/* MM */
	for(ni=0; ni<nMaxIter; ni++)
	{
		/* map 1*/
		dASum = 0.0;
		for(nj=0; nj<nSampleNum; nj++)
			dASum += vA0[nj];

		dS2 = 0.0;
		for(nk=0; nk<nMMax; nk++)
		{
			dS2 += vMHist[nk]/(dASum+nk);
		}

		for(nj=0; nj<nSampleNum; nj++)
		{
			dS1 = vXHist[0][nj];
			for(nk=1; nk<nXMax; nk++)
			{
				dS1 += vXHist[nk][nj]*vA0[nj]/(vA0[nj]+nk);
			}

			vAmid1[nj] = dS1/dS2;
		}

		/* map 2 */
		dASum = 0.0;
		for(nj=0; nj<nSampleNum; nj++)
			dASum += vAmid1[nj];

		dS2 = 0.0;
		for(nk=0; nk<nMMax; nk++)
		{
			dS2 += vMHist[nk]/(dASum+nk);
		}

		for(nj=0; nj<nSampleNum; nj++)
		{
			dS1 = vXHist[0][nj];
			for(nk=1; nk<nXMax; nk++)
			{
				dS1 += vXHist[nk][nj]*vAmid1[nj]/(vAmid1[nj]+nk);
			}
		
			vAmid2[nj] = dS1/dS2;
		}

		/* accelarator */
		dI1 = 0.0;
		dI2 = 0.0;
		for(nj=0; nj<nSampleNum; nj++)
		{
			vA1[nj] = vAmid1[nj]-vA0[nj];
			vAmid1[nj] = vAmid2[nj]-vAmid1[nj]-vA1[nj];
			dI1 += vA1[nj]*vA1[nj];
			dI2 += vA1[nj]*vAmid1[nj];
		}
		dI1 /= dI2;

		nAPositive = 1;
		for(nj=0; nj<nSampleNum; nj++)
		{
			vA1[nj] = vA0[nj]-2.0*dI1*vA1[nj]+dI1*dI1*vAmid1[nj];
			if(vA1[nj] < 0.0)
			{
				nAPositive = 0;
				break;
			}
		}
		
		if(nAPositive == 1)
		{
			dLike1 = SeqClust_Kpolya_MM_LogLike(nSampleNum, nMMax, nXMax, vMHist, vXHist, vA1);
		}
		else
		{
			dLike1 = dLike0 - 1.0;
		}

		if(dLike1 < dLike0)
		{
			for(nj=0; nj<nSampleNum; nj++)
				vA1[nj] = vAmid2[nj];

			dLike1 = SeqClust_Kpolya_MM_LogLike(nSampleNum, nMMax, nXMax, vMHist, vXHist, vA1);
		}

		dErr = fabs(dLike1-dLike0)/(fabs(dLike0)+1.0);

		if(dLike1 < dLike0)
		{
			if(dErr > 1e-6)
			{
				printf("Warning: SeqClust_Kpolya_MM, decreasing likelihood!\n");
			}
			break;
		}

		if(dErr < 1e-9)
			break;

		for(nj=0; nj<nSampleNum; nj++)
			vA0[nj] = vA1[nj];
		dLike0 = dLike1;
	}

	for(nj=0; nj<nSampleNum; nj++)
		vAnew[nj] = vA1[nj];

	/* release memory */
	free(vA0);
	free(vAmid1);
	free(vAmid2);
	free(vA1);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_Kpolya_MM_LogLike()                                           */
/*  Compute Polya loglikelihood for MM algorithm                           */
/* ----------------------------------------------------------------------- */ 
double SeqClust_Kpolya_MM_LogLike(int nSampleNum, int nMMax, int nXMax, 
					   double *vMHist, double **vXHist, double *vA)
{
	/* define */
	double dLike = 0.0;
	int ni,nj;
	double dASum;

	/* compute */
	dASum = 0.0;
	for(nj=0; nj<nSampleNum; nj++)
		dASum += vA[nj];

	for(ni=0; ni<nMMax; ni++)
		dLike -= vMHist[ni]*log(dASum+ni);
	for(ni=0; ni<nXMax; ni++)
	{
		for(nj=0; nj<nSampleNum; nj++)
		{
			dLike += vXHist[ni][nj]*log(vA[nj]+ni);
		}
	}

	/* return */
	return dLike;
}


/* ----------------------------------------------------------------------- */ 
/*  SeqClust_Count_Main()                                                  */
/*  Get read count for genomic regions.                                    */
/* ----------------------------------------------------------------------- */ 
int SeqClust_Count_Main(char strInputPath[], char strDataPath[], char strOutputFile[],
						int nExtendLen)
{
	/* define */
	struct DOUBLEMATRIX *pRegion = NULL;
	struct tagString **vRegionName = NULL;
	struct DOUBLEMATRIX *pSortRegion = NULL;
	struct LONGMATRIX *pSortId = NULL;
	struct tagString **vChrName = NULL;
	struct tagString **vFileName = NULL;
	int nChrNum = 0;
	int ni;
	long nk;
	int nFileNum = 0;
	FILE *fpData;
	char strLine[LONG_LINE_LENGTH];
	char strFileName[MED_LINE_LENGTH];
	double **vData = NULL;
	double **vDataS = NULL;

	/* STEP 1: LOAD AND SORT REGIONS */
	printf("Loading regions ...\n");
	SeqClust_Count_LoadRegion(strInputPath, &pRegion, &vRegionName, &pSortRegion, &pSortId, &vChrName, &nChrNum);

	/* STEP 2: COUNT */
	/* count file number */
	printf("Counting ...\n");

	fpData = NULL;
	fpData = fopen(strDataPath, "r");
	if(fpData == NULL)
	{
		printf("Error: SeqClust_Count_Main, cannot open sample list file!\n");
		exit(EXIT_FAILURE);
	}
	while(fgets(strLine, LONG_LINE_LENGTH, fpData)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;
		
		nFileNum++;
	}
	fclose(fpData);

	/* prepare memory */ 
	vData = (double **)calloc(pRegion->nHeight, sizeof(double *));
	vDataS = (double **)calloc(pRegion->nHeight, sizeof(double *));
	if( (vData == NULL) || (vDataS == NULL) )
	{
		printf("Error: SeqClust_Count_Main, cannot create memory for count data!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<pRegion->nHeight; ni++)
	{
		vData[ni] = (double *)calloc(nFileNum, sizeof(double));
		if(vData[ni] == NULL)
		{
			printf("Error: SeqClust_Count_Main, cannot create memory for count data!\n");
			exit(EXIT_FAILURE);
		}
	}

	vFileName = (struct tagString **)calloc(nFileNum, sizeof(struct tagString *));
	if(vFileName == NULL)
	{
		printf("Error: SeqClust_Count_Main, cannot create memory for file names!\n");
		exit(EXIT_FAILURE);
	}

	/* count reads */
	fpData = NULL;
	fpData = fopen(strDataPath, "r");
	if(fpData == NULL)
	{
		printf("Error: SeqClust_Count_Main, cannot open sample list file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpData)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;
		
		printf(" %s ...\n", strLine);
		GetFileName(strLine, strFileName);
		StringAddTail(vFileName+ni, strFileName);
		SeqClust_Count_SingleSample(strLine, pSortRegion, vChrName, nChrNum, vData, ni, nExtendLen);

		ni++;
	}
	fclose(fpData);
	
	if(ni != nFileNum)
	{
		printf("Error: SeqClust_Count_Main, inconsistent file number!\n");
		exit(EXIT_FAILURE);
	}


	/* STEP 3: RESORT REGIONS */
	for(ni=0; ni<pRegion->nHeight; ni++)
	{
		nk = pSortId->pMatElement[ni];
		vDataS[nk] = vData[ni];
		vData[ni] = NULL;
	}

	/* STEP 4: EXPORT */
	printf("Export data ...\n");
	SeqClust_Count_ExportCod(strOutputFile, pRegion, vRegionName, vChrName, nChrNum, vDataS, 
		nFileNum, vFileName);

	/* STEP 5: CLEAR MEMORY */
	for(ni=0; ni<nFileNum; ni++)
	{
		DeleteString(vFileName[ni]);
		vFileName[ni] = NULL;
	}
	free(vFileName);
	for(ni=0; ni<pRegion->nHeight; ni++)
	{
		free(vDataS[ni]);
	}
	free(vData);
	free(vDataS);
	for(ni=0; ni<nChrNum; ni++)
	{
		DeleteString(vChrName[ni]);
		vChrName[ni] = NULL;
	}
	free(vChrName);
	for(ni=0; ni<pRegion->nHeight; ni++)
	{
		DeleteString(vRegionName[ni]);
		vRegionName[ni] = NULL;
	}
	free(vRegionName);

	DestroyDoubleMatrix(pRegion);
	DestroyDoubleMatrix(pSortRegion);
	DestroyLongMatrix(pSortId);
	
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_Count_LoadRegion()                                            */
/*  Load input coordinates.                                                */
/* ----------------------------------------------------------------------- */ 
int SeqClust_Count_LoadRegion(char strInputPath[], struct DOUBLEMATRIX **pRegion,
				struct tagString ***vvRegionName,
				struct DOUBLEMATRIX **pSortRegion, struct LONGMATRIX **pSortId,
				struct tagString ***vvChrName, int *pChrNum)
{
	/* define */
	struct INTMATRIX *pType = NULL;
	struct INTMATRIX *pPriority = NULL;
	int nMaxChrNum = 65535;	
	int nChrNum = 0;
	struct tagProbeGenomeInfo **vChrList;
	int nRegionNum = 0;

	int ni,nId;
	FILE *fpIn;
	char strLine[LONG_LINE_LENGTH];
	
	char strAlias[LINE_LENGTH];
	char strChr[LINE_LENGTH];
	int nPos1,nPos2;
	int nStrand;

	/* STEP 1: init chromosome index */
	vChrList = NULL;
	vChrList = (struct tagProbeGenomeInfo **)calloc(nMaxChrNum, sizeof(struct tagProbeGenomeInfo *));
	if(vChrList == NULL)
	{
		printf("Error: SeqClust_Count_LoadRegion, cannot create vector for chromosomes!\n");
		exit(EXIT_FAILURE);
	}
	nChrNum = 0;

	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: SeqClust_Count_LoadRegion, cannot open input file list!\n");
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

		SeqClust_ReadPos_From_Cod(strLine, strAlias, strChr, &nPos1, &nPos2, &nStrand);
		nId = SeqClust_FindChrInfo(vChrList, nMaxChrNum, &nChrNum, strChr);
		if( (nId < 0) || (nId >= nMaxChrNum) )
		{
			printf("Error: SeqClust_Count_LoadRegion, cannot find the matching chromosome!\n");
			exit(EXIT_FAILURE);
		}
		
		nRegionNum++;
	}
	
	fclose(fpIn);
	

	/* return chromosome name and length */
	*pChrNum = nChrNum;
	*vvChrName = (struct tagString **)calloc(nChrNum, sizeof(struct tagString *));
	if(*vvChrName == NULL)
	{
		printf("Error: SeqClust_Count_LoadRegion, cannot create memory for chromosome name!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nChrNum; ni++)
	{
		(*vvChrName)[ni] = vChrList[ni]->pProbe;
		vChrList[ni]->pProbe = NULL;
		ProbeGenomeInfoDestroy(vChrList+ni);
	}
	free(vChrList);

	if(ni != nChrNum)
	{
		printf("Error: SeqClust_Count_LoadRegion, inconsistent chromosome number!\n");
		exit(EXIT_FAILURE);
	}

	/* Step 2: load regions */
	*pRegion = NULL;
	*pRegion = CreateDoubleMatrix(nRegionNum, 4);
	if(*pRegion == NULL)
	{
		printf("Error: SeqClust_Count_LoadRegion, cannot allocate memory for regions!\n");
		exit(EXIT_FAILURE);
	}

	*vvRegionName = NULL;
	*vvRegionName = (struct tagString **)calloc(nRegionNum, sizeof(struct tagString *));
	if(*vvRegionName == NULL)
	{
		printf("Error: SeqClust_Count_LoadRegion, cannot allocate memory for region names!\n");
		exit(EXIT_FAILURE);
	}

	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: SeqClust_Count_LoadRegion, cannot open input file list!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		SeqClust_ReadPos_From_Cod(strLine, strAlias, strChr, &nPos1, &nPos2, &nStrand);

		nId = SeqClust_FindChr((*vvChrName), nChrNum, strChr);
		if( (nId < 0) || (nId >= nChrNum) )
		{
			printf("Error: SeqClust_Count_LoadRegion, cannot find the matching chromosome!\n");
			exit(EXIT_FAILURE);
		}

		DMSETAT(*pRegion, ni, 0, nId);
		DMSETAT(*pRegion, ni, 1, nPos1);
		DMSETAT(*pRegion, ni, 2, nPos2);
		DMSETAT(*pRegion, ni, 3, nStrand);
		StringAddTail((*vvRegionName)+ni, strAlias);

		ni++;
	}
	
	fclose(fpIn);
	
	if(ni != nRegionNum)
	{
		printf("Error: SeqClust_Count_LoadRegion, inconsistent region number!\n");
		exit(EXIT_FAILURE);
	}


	/* Step 3: sort */
	pType = CreateIntMatrix(1, 4);
	if(pType == NULL)
	{
		printf("Error: SeqClust_Count_LoadRegion, cannot create field type for the input region list!\n");
		exit(EXIT_FAILURE);
	}
	pType->pMatElement[0] = 2;
	pType->pMatElement[1] = 2;
	pType->pMatElement[2] = 2;
	pType->pMatElement[3] = 2;

	pPriority = CreateIntMatrix(1, 3);
	if(pPriority == NULL)
	{
		printf("Error: SeqClust_Count_LoadRegion, cannot create priority level for the input region list!\n");
		exit(EXIT_FAILURE);
	}
	pPriority->pMatElement[0] = 0;
	pPriority->pMatElement[1] = 1;
	pPriority->pMatElement[2] = 2;

	DMSORTROWS( (*pRegion), pType, pPriority, pSortRegion, pSortId);

	/* release memory */
	DestroyIntMatrix(pType);
	DestroyIntMatrix(pPriority);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_Count_SingleSample()                                          */
/*  Count read for a single file.                                          */
/* ----------------------------------------------------------------------- */ 
int SeqClust_Count_SingleSample(char strFileName[], struct DOUBLEMATRIX *pRegion, 
							struct tagString **vChrName, int nChrNum, 
							double **vData, int nDataCol, int nExtendLen)
{
	/* define */
	FILE *fpIn = NULL;
	char strLine[LONG_LINE_LENGTH];
	char strChr[LINE_LENGTH];
	int nPos;
	int nStrand;
	int nPS,nPE;
	int nChrId;
	int nx,ny,nz;
	int nCmp;
	int nHit;


	/* read */
	fpIn = fopen(strFileName, "r");
	if(fpIn == NULL)
	{
		printf("Error: SeqClust_Count_SingleSample, cannot open file %s\n", strFileName);
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

		SeqClust_ReadPos_From_Aln(strLine, strChr, &nPos, &nStrand);
		nChrId = SeqClust_FindChr(vChrName, nChrNum, strChr);
		if( (nChrId < 0) || (nChrId >= nChrNum) )
		{
			continue;
		}

		if(nStrand == 1)
		{
			nPS = nPos-nExtendLen;
			nPE = nPos;
		}
		else
		{
			nPS = nPos;
			nPE = nPos+nExtendLen;
		}

		/* find matching regions */
		nx = 0;
		ny = pRegion->nHeight-1;

		nHit = 0;

		nCmp = SeqClust_Count_RegCompare(nChrId, nPS, nPE, pRegion, nx);
		if( nCmp < 0 )
		{
			continue;
		}
		else if( (nCmp == 0) || (nCmp == 1) )
		{
			nHit = 1;
			nz = nx;
		}
		
		if(nHit == 0)
		{
			nCmp = SeqClust_Count_RegCompare(nChrId, nPS, nPE, pRegion, ny);
			if( nCmp > 1 )
			{
				continue;
			}
			else if( (nCmp == 0) || (nCmp == 1) )
			{
				nHit = 1;
				nz = ny;
			}
		}

		if(nHit == 0)
		{
			while(ny-nx > 1)
			{
				nz = (int)((nx+ny)/2);
				nCmp = SeqClust_Count_RegCompare(nChrId, nPS, nPE, pRegion, nz);

				if(nCmp < 0)
				{
					ny = nz;
				}
				else if(nCmp > 1)
				{
					nx = nz;
				}
				else
				{
					nHit = 1;
					break;
				}
			}
		}

		if(nHit == 0)
		{
			nCmp = SeqClust_Count_RegCompare(nChrId, nPS, nPE, pRegion, nx);
			if( (nCmp == 0) || (nCmp == 1) )
			{
				nz = nx;
				nHit = 1;
			}

			nCmp = SeqClust_Count_RegCompare(nChrId, nPS, nPE, pRegion, ny);
			if( (nCmp == 0) || (nCmp == 1) )
			{
				nz = ny;
				nHit = 1;
			}
		}

		if(nHit == 1)
		{

			vData[nz][nDataCol] += 1;

			for(nx=nz-1; nx>=0; nx--)
			{
				nCmp = SeqClust_Count_RegCompare(nChrId, nPS, nPE, pRegion, nx);
				if(nCmp > 1)
					break;

				if((nCmp == 0) || (nCmp == 1))
					vData[nx][nDataCol] += 1;
			}

			for(ny=nz+1; ny<pRegion->nHeight; ny++)
			{
				nCmp = SeqClust_Count_RegCompare(nChrId, nPS, nPE, pRegion, ny);
				if(nCmp < 0)
					break;

				if((nCmp == 0) || (nCmp == 1))
					vData[ny][nDataCol] += 1;
			}
		}
	}

	fclose(fpIn);


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_Count_ExportCod()                                             */
/*  Save count data to output files.                                       */
/* ----------------------------------------------------------------------- */ 
int SeqClust_Count_ExportCod(char strOutputFile[], struct DOUBLEMATRIX *pRegion, 
							 struct tagString **vRegionName, 
							 struct tagString **vChrName, int nChrNum, 
							 double **vDataS, int nFileNum,
							 struct tagString **vFileName)
{
	/* define */
	FILE *fpOut = NULL;
	int ni,nj;
	int nChr;
	int nPos1,nPos2,nStrand;
	char chStrand;

	/* write */
	fpOut = NULL;
	fpOut = fopen(strOutputFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: SeqClust_Count_ExportCod, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "#name\tchr\tstart\tend\tstrand");
	for(ni=0; ni<nFileNum; ni++)
	{
		fprintf(fpOut, "\t%s", vFileName[ni]->m_pString);
	}
	fprintf(fpOut, "\n");

	for(ni=0; ni<pRegion->nHeight; ni++)
	{
		nChr = (int)(DMGETAT(pRegion, ni, 0));
		nPos1 = (int)(DMGETAT(pRegion, ni, 1));
		nPos2 = (int)(DMGETAT(pRegion, ni, 2));
		nStrand = (int)(DMGETAT(pRegion, ni, 3));
		if(nStrand == 1)
			chStrand = '-';
		else
			chStrand = '+';

		fprintf(fpOut, "%s\t%s\t%d\t%d\t%c", vRegionName[ni]->m_pString, 
			vChrName[nChr]->m_pString, nPos1, nPos2, chStrand);
		
		for(nj=0; nj<nFileNum; nj++)
		{
			fprintf(fpOut, "\t%d", (int)(vDataS[ni][nj]));
		}

		fprintf(fpOut, "\n");
	}

	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClust_Count_RegCompare()                                            */
/*  Compare two regions.                                                   */
/* ----------------------------------------------------------------------- */ 
int SeqClust_Count_RegCompare(int nChrId, int nStart, int nEnd, 
							  struct DOUBLEMATRIX *pRegion, int nRegionId)
{
	/* define */
	int nCmp = 0;
	int nChr;
	int nP1,nP2;

	nChr = (int)(DMGETAT(pRegion, nRegionId, 0)+0.01);
	nP1 = (int)(DMGETAT(pRegion, nRegionId, 1));
	nP2 = (int)(DMGETAT(pRegion, nRegionId, 2));

	if(nChrId < nChr)
	{
		nCmp = -1;
	}
	else if(nChrId > nChr)
	{
		nCmp = 2;
	}
	else
	{
		if(nEnd < nP1)
		{
			nCmp = -1;
		}
		else if(nStart > nP2)
		{
			nCmp = 2;
		}
		else if(nP1 >= nStart)
		{
			nCmp = 0;
		}
		else
		{
			nCmp = 1;
		}
	}

	/* return */
	return nCmp;
}

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
		int nMinClustNum, int nMaxClustNum)
{
	/* define */
	double **vDist = NULL;
	double **vTempDist = NULL;
	double **vTempDistNew = NULL;
	int ni,nj,nHeight,nx,ny;
	struct tagHCNode **vHCNode = NULL;
	struct tagHCNode *pNewNode = NULL;
	int nRootNum = nColNum;
	char strFileName[MED_LINE_LENGTH];
	int nM1,nM2,nLeafNum;
	double dMinDist;
	int nTotalNode = 0;
	double dNewDist;

	/* ----------------------- */
	/* step1: compute distance */
	/* ----------------------- */
	vDist = (double **)calloc(nColNum, sizeof(double *));
	if(vDist == NULL)
	{
		printf("Error: HierarchicalCluster_Column, cannot create memory for computing distance!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nColNum; ni++)
	{
		vDist[ni] = (double *)calloc(nColNum, sizeof(double));
		if(vDist[ni] == NULL)
		{
			printf("Error: HierarchicalCluster_Column, cannot create memory for computing distance!\n");
			exit(EXIT_FAILURE);
		}
	}

	if(nDistanceType == 0)
	{
		HierarchicalCluster_Column_Dist_Corr(nRowNum, nColNum, vData, vDist);
	}
	else if(nDistanceType == 1)
	{
		HierarchicalCluster_Column_Dist_AbsCorr(nRowNum, nColNum, vData, vDist);
	}

	/* ----------------------- */
	/* step2: clustering       */
	/* ----------------------- */
	vHCNode = (struct tagHCNode **)calloc(nColNum, sizeof(struct tagHCNode *));
	if(vHCNode == NULL)
	{
		printf("Error: HierarchicalCluster_Column, cannot create memory for cluster node vector!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nColNum; ni++)
	{
		vHCNode[ni] = HCNode_Create(nTotalNode, 0, 0, 1);
		if(vHCNode[ni] == NULL)
		{
			printf("Error: HierarchicalCluster_Column, cannot create memory for cluster node!\n");
			exit(EXIT_FAILURE);
		}
		vHCNode[ni]->vLeafId[0] = ni;
		nTotalNode++;
	}
	nHeight = 1;

	vTempDist = NULL;
	vTempDist = (double **)calloc(nColNum, sizeof(double *));
	if(vTempDist == NULL)
	{
		printf("Error: HierarchicalCluster_Column, cannot create memory for temporary distance!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nColNum; ni++)
	{
		vTempDist[ni] = (double *)calloc(nColNum, sizeof(double));
		if(vTempDist[ni] == NULL)
		{
			printf("Error: HierarchicalCluster_Column, cannot create memory for temporary distance!\n");
			exit(EXIT_FAILURE);
		}
		memcpy(vTempDist[ni], vDist[ni], nColNum*sizeof(double));
	}

	/* export clusters */
	if( (nExportCluster == 1) && (nRootNum >= nMinClustNum) && (nRootNum <= nMaxClustNum) )
	{
		sprintf(strFileName, "%s%s_sampleblock_B%d.txt", strOutputPath, strOutputFile, nRootNum);
		HierarchicalCluster_ExportCluster(vHCNode, nRootNum, strFileName);
	}

	while(nRootNum > 1)
	{
		/* find closest node pair */
		dMinDist = vTempDist[0][1];
		nM1 = 0;
		nM2 = 1;
		for(ni=0; ni<nRootNum; ni++)
		{
			for(nj=ni+1; nj<nRootNum; nj++)
			{
				if(vTempDist[ni][nj] < dMinDist)
				{
					nM1 = ni;
					nM2 = nj;
					dMinDist = vTempDist[ni][nj];
				}
			}
		}

		/* merge node */
		nLeafNum = vHCNode[nM1]->nLeafNum + vHCNode[nM2]->nLeafNum;
		pNewNode = NULL;
		pNewNode = HCNode_Create(nTotalNode, nHeight, 2, nLeafNum);
		if(pNewNode == NULL)
		{
			printf("Error: HierarchicalCluster_Column, cannot create memory for new node!\n");
			exit(EXIT_FAILURE);
		}
		memcpy(pNewNode->vLeafId, vHCNode[nM1]->vLeafId, vHCNode[nM1]->nLeafNum*sizeof(int));
		memcpy(pNewNode->vLeafId+vHCNode[nM1]->nLeafNum, vHCNode[nM2]->vLeafId, vHCNode[nM2]->nLeafNum*sizeof(int));
		free(vHCNode[nM1]->vLeafId);
		vHCNode[nM1]->vLeafId = NULL;
		vHCNode[nM1]->nLeafNum = 0;

		free(vHCNode[nM2]->vLeafId);
		vHCNode[nM2]->vLeafId = NULL;
		vHCNode[nM2]->nLeafNum = 0;

		pNewNode->vChildNode[0] = vHCNode[nM1];
		pNewNode->vChildNode[1] = vHCNode[nM2];
		
		vHCNode[nM1] = NULL;
		vHCNode[nM2] = NULL;

		nTotalNode++;

		/* recompute distance */
		vTempDistNew = NULL;
		vTempDistNew = (double **)calloc(nRootNum-1, sizeof(double *));
		if(vTempDistNew == NULL)
		{
			printf("Error: HierarchicalCluster_Column, cannot create memory for temporary distance!\n");
			exit(EXIT_FAILURE);
		}
		for(ni=0; ni<(nRootNum-1); ni++)
		{
			vTempDistNew[ni] = (double *)calloc((nRootNum-1), sizeof(double));
			if(vTempDistNew[ni] == NULL)
			{
				printf("Error: HierarchicalCluster_Column, cannot create memory for temporary distance!\n");
				exit(EXIT_FAILURE);
			}
		}

		for(ni=0; ni<nRootNum; ni++)
		{
			if( (ni == nM1) || (ni == nM2) )
				continue;

			if( (ni > nM1) && (ni < nM2) )
			{
				nx = ni-1;
				vHCNode[nx] = vHCNode[ni];
				vHCNode[ni] = NULL;
			}
			else if(ni > nM2)
			{
				nx = ni-2;
				vHCNode[nx] = vHCNode[ni];
				vHCNode[ni] = NULL;
			}
			else
			{
				nx = ni;
			}

			for(nj=ni+1; nj<nRootNum; nj++)
			{
				if( (nj == nM1) || (nj == nM2) )
					continue;

				if( (nj > nM1) && (nj < nM2) )
					ny = nj-1;
				else if(nj > nM2)
					ny = nj-2;
				else
					ny = nj;

				vTempDistNew[nx][ny] = vTempDist[ni][nj];
				vTempDistNew[ny][nx] = vTempDist[ni][nj];
			}
		}

		ni = nRootNum-2;
		vHCNode[ni] = pNewNode;
		for(nj=0; nj<(nRootNum-2); nj++)
		{
			/* average linkage */
			if(nMergeType == 0)
			{
				dNewDist = HierarchicalCluster_UpdateDist_AverageLinkage(pNewNode, vHCNode[nj], vDist);
			}
			/* complete linkage */
			/* single linkage */

			vTempDistNew[ni][nj] = dNewDist;
			vTempDistNew[nj][ni] = dNewDist;
		}

		for(ni=0; ni<nRootNum; ni++)
		{
			free(vTempDist[ni]);
		}
		
		free(vTempDist);
		vTempDist = vTempDistNew;
		vTempDistNew = NULL;

		/* update height */
		nHeight++;
		nRootNum--;

		/* export clusters */
		if( (nExportCluster == 1) && (nRootNum >= nMinClustNum) && (nRootNum <= nMaxClustNum) )
		{
			sprintf(strFileName, "%s%s_sampleblock_B%d.txt", strOutputPath, strOutputFile, nRootNum);
			HierarchicalCluster_ExportCluster(vHCNode, nRootNum, strFileName);
		}
	}

	/* free memory */
	for(ni=0; ni<nRootNum; ni++)
	{
		free(vTempDist[ni]);
	}
	free(vTempDist);
	
	/* CLEAR ALL NODES IN THE TREE */
	HCNode_ClearTree(vHCNode);
	free(vHCNode);

	for(ni=0; ni<nColNum; ni++)
	{
		free(vDist[ni]);
	}
	free(vDist);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HierarchicalCluster_Column_Dist_Corr()                                 */
/*  Compute correlation distance.                                          */
/* ----------------------------------------------------------------------- */ 
int HierarchicalCluster_Column_Dist_Corr(int nRowNum, int nColNum, 
					double **vData, double **vDist)
{
	/* define */
	double *vM;
	double *vS;
	double *vT;
	int ni,nj,nk;

	/* init */
	vM = (double *)calloc(nColNum, sizeof(double));
	vS = (double *)calloc(nColNum, sizeof(double));
	vT = (double *)calloc(nColNum, sizeof(double));
	if( (vM == NULL) || (vS == NULL) || (vT == NULL))
	{
		printf("Error: HierarchicalCluster_Column_Dist_Corr, cannot allocate memory for computing correlation distance!\n");
		exit(EXIT_FAILURE);
	}

	/* compute mean */
	for(ni=0; ni<nRowNum; ni++)
	{
		for(nj=0; nj<nColNum; nj++)
		{
			vM[nj] += vData[ni][nj];
		}
	}

	for(nj=0; nj<nColNum; nj++)
		vM[nj] /= nRowNum;

	for(ni=0; ni<nRowNum; ni++)
	{
		for(nj=0; nj<nColNum; nj++)
		{
			vT[nj] = vData[ni][nj]-vM[nj];
		}

		for(nj=0; nj<nColNum; nj++)
		{
			vDist[nj][nj] += vT[nj]*vT[nj];
			for(nk=nj+1; nk<nColNum; nk++)
			{
				vDist[nj][nk] += vT[nj]*vT[nk];
			}
		}
	}

	for(nj=0; nj<nColNum; nj++)
	{
		vS[nj] = sqrt(vDist[nj][nj]);
	}

	for(nj=0; nj<nColNum; nj++)
	{
		vDist[nj][nj] = 0.0;
		for(nk=nj+1; nk<nColNum; nk++)
		{
			vDist[nj][nk] = vDist[nj][nk]/(vS[nj]*vS[nk]+1e-6);
			vDist[nj][nk] = (1.0-vDist[nj][nk])/2.0;
			vDist[nk][nj] = vDist[nj][nk];
		}
	}

	/* release memory */
	free(vM);
	free(vS);
	free(vT);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HierarchicalCluster_Column_Dist_AbsCorr()                              */
/*  Compute absolute correlation distance.                                 */
/* ----------------------------------------------------------------------- */ 
int HierarchicalCluster_Column_Dist_AbsCorr(int nRowNum, int nColNum, 
					double **vData, double **vDist)
{
	/* define */
	int nj,nk;

	/* init */
	HierarchicalCluster_Column_Dist_Corr(nRowNum, nColNum, vData, vDist);

	for(nj=0; nj<nColNum; nj++)
	{
		for(nk=nj+1; nk<nColNum; nk++)
		{
			vDist[nj][nk] = 1.0-vDist[nj][nk]*2.0;
			vDist[nj][nk] = 1.0-fabs(vDist[nj][nk]);
			vDist[nk][nj] = vDist[nj][nk];
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HCNode_Create()                                                        */
/*  Create node for hierarchical clustering.                               */
/* ----------------------------------------------------------------------- */ 
struct tagHCNode *HCNode_Create(int nNodeId, int nHeight, int nChildNum, int nLeafNum)
{
	/* define */
	struct tagHCNode *pNode = NULL;

	pNode = (struct tagHCNode *)calloc(1, sizeof(struct tagHCNode));
	if(pNode == NULL)
	{
		printf("Warning: HCNode_Create, cannot create HCNode!\n");
		return NULL;
	}

	pNode->nNodeId = nNodeId;
	pNode->nHeight = nHeight;
	pNode->nChildNum = nChildNum;
	pNode->nLeafNum = nLeafNum;

	pNode->vChildNode = NULL;
	if(nChildNum > 0)
	{
		pNode->vChildNode = (struct tagHCNode **)calloc(nChildNum, sizeof(struct tagHCNode *));
		if(pNode->vChildNode == NULL)
		{
			printf("Error: HCNode_Create, cannot create memory for child nodes!\n");
			return NULL;
		}
	}

	pNode->vLeafId = NULL;
	if(nLeafNum > 0)
	{
		pNode->vLeafId = (int *)calloc(nLeafNum, sizeof(int));
		if(pNode->vLeafId == NULL)
		{
			printf("Error: HCNode_Create, cannot create memory for leaf nodes!\n");
			return NULL;
		}
	}

	/* return */
	return pNode;
}

/* ----------------------------------------------------------------------- */ 
/*  HCNode_Delete()                                                        */
/*  delete a node for hierarchical clustering.                             */
/* ----------------------------------------------------------------------- */ 
int HCNode_Delete(struct tagHCNode **pNode)
{
	/* delete */

	if(pNode == NULL)
		return PROC_SUCCESS;

	if(*pNode == NULL)
		return PROC_SUCCESS;

	if( (*pNode)->vChildNode != NULL)
	{
		free((*pNode)->vChildNode);
		(*pNode)->nChildNum = 0;
	}

	if( (*pNode)->vLeafId != NULL)
	{
		free((*pNode)->vLeafId);
		(*pNode)->nLeafNum = 0;
	}

	free(*pNode);
	*pNode = NULL;

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  HCNode_ClearTree()                                                     */
/*  delete all nodes in a tree.                                            */
/* ----------------------------------------------------------------------- */ 
int HCNode_ClearTree(struct tagHCNode **pNode)
{
	/* define */
	int ni;

	/* delete */
	if(pNode == NULL)
		return PROC_SUCCESS;

	if(*pNode == NULL)
		return PROC_SUCCESS;

	for(ni=0; ni<(*pNode)->nChildNum; ni++)
	{
		HCNode_ClearTree((*pNode)->vChildNode+ni);
	}

	if( (*pNode)->vChildNode != NULL)
	{
		free((*pNode)->vChildNode);
		(*pNode)->nChildNum = 0;
	}

	if( (*pNode)->vLeafId != NULL)
	{
		free((*pNode)->vLeafId);
		(*pNode)->nLeafNum = 0;
	}

	free(*pNode);
	*pNode = NULL;

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  HierarchicalCluster_UpdateDist_AverageLinkage()                        */
/*  compute average linkage distance.        .                             */
/* ----------------------------------------------------------------------- */ 
double HierarchicalCluster_UpdateDist_AverageLinkage(struct tagHCNode *pNode1,
				struct tagHCNode *pNode2, double **vDist)
{
	/* define */
	double dD = 0.0;
	int ni,nj,nx,ny;

	if( (pNode1->nLeafNum <= 0) || (pNode2->nLeafNum <= 0) )
		return dD;

	for(ni=0; ni<pNode1->nLeafNum; ni++)
	{
		nx = pNode1->vLeafId[ni];
		for(nj=0; nj<pNode2->nLeafNum; nj++)
		{
			ny = pNode2->vLeafId[nj];
			dD += vDist[nx][ny];
		}
	}

	dD = dD/pNode1->nLeafNum/pNode2->nLeafNum;

	/* return */
	return dD;
}


/* ----------------------------------------------------------------------- */ 
/*  HierarchicalCluster_UpdateDist_AverageLinkage()                        */
/*  compute average linkage distance.        .                             */
/* ----------------------------------------------------------------------- */ 
int HierarchicalCluster_ExportCluster(struct tagHCNode **vHCNode, 
					int nNodeNum, char strOutFile[])
{
	/* define */
	FILE *fpOut; 
	int ni,nj;

	/* export */
	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: HierarchicalCluster_ExportCluster, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	
	fprintf(fpOut, "%d\n", nNodeNum);
	for(ni=0; ni<nNodeNum; ni++)
	{
		fprintf(fpOut, ">%d\n", vHCNode[ni]->nLeafNum);
		for(nj=0; nj<vHCNode[ni]->nLeafNum; nj++)
		{
			fprintf(fpOut, "%d\n", vHCNode[ni]->vLeafId[nj]);
		}
	}

	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}

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
				  int nBlockLenCut, int nExportBAR)
{
	/* define */
	int nResult;
	int ni;
	char strFileName[MED_LINE_LENGTH];

	/* For Step 1*/
	int nFileNum = 0;
	int nChrNum = 0;
	struct tagString **vChrName = NULL;
	struct tagString **vFileName = NULL;
	struct INTMATRIX *pChrLen = NULL;
	struct DOUBLEMATRIX *pFileReadCount = NULL;
	struct BYTEMATRIX *pIncludeChr = NULL;
	struct DOUBLEMATRIX *pFileNormFactor = NULL;
	struct DOUBLEMATRIX *pSmoothK = NULL;
	int nSmoothStart = 0;
	int nSmoothNum = 0;
	double dSmoothFactor = 0.0;
	int nBaseLineFileId = 0;
	
	/* For Step 3 */
	double dBinMax = 0.0;
	double dBinMin = 1000000.0;
	double dCL = 0.0;
	double dCH = 0.0;



	/* ---------------------------------------------------------*/
	/* STEP1: First scan, find chromosome number and length     */
	/*        Find read number for each sample                  */
	/* ---------------------------------------------------------*/
	printf("/* ------------------------------------------------- */\n");
	printf("/* STEP 1: Initialize                                */\n");
	printf("/* ------------------------------------------------- */\n");
	rand_u_init(13);
	AdjustDirectoryPath(strOutputPath);
	dSmoothFactor = 1.0;
	if(nKernelType > 0)
	{
		nResult = SeqClust_InitializeKernel(&pSmoothK, nKernelType, 
			nKernelStep, nKernelLen, nKernelBand, &nSmoothStart, &nSmoothNum);
		dSmoothFactor = 0.0;
		for(ni=0; ni<pSmoothK->nWidth; ni++)
		{
			dSmoothFactor += pSmoothK->pMatElement[ni];
		}
	}
	printf("Smooth Factor = %f\n", dSmoothFactor); 
	
	nResult = SeqClust_Initialize(strInputPath, &nFileNum, &vFileName, 
			&pFileReadCount, &pFileNormFactor, &nBaseLineFileId, &nChrNum, &vChrName, &pChrLen);
	for(ni=0; ni<pFileReadCount->nWidth; ni++)
	{
		pFileReadCount->pMatElement[ni] *= dSmoothFactor;
	}


	/* ---------------------------------------------------------*/
	/* STEP2: Count reads for genomic bins                      */
	/* ---------------------------------------------------------*/
	printf("\n");
	printf("/* ------------------------------------------------- */\n");
	printf("/* STEP 2: Count reads for genomic bins              */\n");
	printf("/* ------------------------------------------------- */\n");
	pIncludeChr = NULL;
	pIncludeChr = CreateByteMatrix(1, nChrNum);
	if(pIncludeChr == NULL)
	{
		printf("Error: SeqClustDP_Seg_Main, cannot create chromosome inclusion vector!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nChrNum; ni++)
		pIncludeChr->pMatElement[ni] = 1;
	nResult = SeqClust_CountBinReads(strInputPath, strOutputPath,
		strOutputFile, nFileNum, nChrNum, vChrName, pChrLen, nBinSize,
		nKernelType, pSmoothK, nSmoothStart, nKernelStep, nSmoothNum, pIncludeChr,
		nExportBAR);
	

	/* ---------------------------------------------------------*/
	/* STEP3 - STEP 5                                           */
	/* ---------------------------------------------------------*/
	if(nSegType != 0)
	{
		printf("\n");
		printf("/* ------------------------------------------------- */\n");
		printf("/* SKIP STEP 3 & 4. Use user supplied intervals.     */\n");
		printf("/* ------------------------------------------------- */\n");
	}
	/* automatic segmentation */
	else
	{
		/* ---------------------------------------------------------*/
		/* STEP3: Compute summary statistic for genomic bins        */
		/* ---------------------------------------------------------*/
		printf("\n");
		printf("/* ------------------------------------------------- */\n");
		printf("/* STEP 3: Compute summary statistics for bins       */\n");
		printf("/* ------------------------------------------------- */\n");
		nResult = SeqClustDP_ComputeBinStats(strOutputPath, strOutputFile, 
			nFileNum, vFileName, pFileNormFactor, nChrNum, vChrName, pChrLen,
			nBinSize, nMemBlock, nCorrBlock, &dBinMax, &dBinMin, nExportBAR);
	

		/* cutoff by Percentile */
		if(nCutType == 0)
		{
			/* dBinMax = 34.0;
			dBinMin = -1.9; */
			SeqClustDP_QuantileCutoff(strOutputPath, strOutputFile, 
				nChrNum, vChrName, pChrLen, nBinSize, dBinMax, dBinMin, 
				nGridNum, dCutL, dCutH, &dCL, &dCH, pIncludeChr);
		}
		/* cutoff by FDR */
		/* else if(nCutType == 1)
		{
			SeqClust_FDRCutoff(strOutputPath, strOutputFile, 
				nChrNum, vChrName, pChrLen, nBinSize, dBinMax, dBinMin, 
				nGridNum, dCutL, dCutH, &dCL, &dCH, pIncludeChr);
		} */
		/* cutoff supplied by users */
		else
		{
			dCL = dCutL;
			dCH = dCutH;
		}
		printf("Max. Bin Stat. = %f\n", dBinMax);
		printf("Min. Bin Stat. = %f\n", dBinMin);
		printf("Higher Cutoff = %f\n", dCH);
		printf("Lower Cutoff = %f\n", dCL);
	

		/* ---------------------------------------------------------*/
		/* STEP4: Genome segmentation                               */
		/* ---------------------------------------------------------*/
		printf("\n");
		printf("/* ------------------------------------------------- */\n");
		printf("/* STEP 4: Find correlation blocks                   */\n");
		printf("/* ------------------------------------------------- */\n");
		/* for debug */
		SeqClustDP_Segmentation(strOutputPath, strOutputFile, 
				nChrNum, vChrName, pChrLen, nBinSize, dCL, dCH);

		/* ---------------------------------------------------------*/
		/* STEP5: Update normalizing factors                        */
		/* ---------------------------------------------------------*/
		printf("\n");
		printf("/* ------------------------------------------------- */\n");
		printf("/* STEP 5: Update normalizing factors                */\n");
		printf("/* ------------------------------------------------- */\n");

		SeqClust_CollectData_RefineFileNormFactor(strOutputPath, strOutputFile, 
				nFileNum, vFileName, nBaseLineFileId,
				pFileReadCount, pFileNormFactor, dSmoothFactor,
				nChrNum, vChrName, pChrLen, nBinSize);

		/* ---------------------------------------------------------*/
		/* STEP6: Detect correlation blocks by DP                   */
		/* ---------------------------------------------------------*/
		printf("\n");
		printf("/* ------------------------------------------------- */\n");
		printf("/* STEP 6: Detect correlation blocks by DP           */\n");
		printf("/* ------------------------------------------------- */\n");
		SeqClustDP_FindBlock(strOutputPath, strOutputFile, 
				nFileNum, vFileName, pFileNormFactor, dSmoothFactor,
				nChrNum, vChrName, pChrLen, nBinSize, pIncludeChr,
				nBlockLenCut);

		/* ---------------------------------------------------------*/
		/* STEP7: Collect data for pattern analysis                 */
		/* ---------------------------------------------------------*/
		printf("\n");
		printf("/* ------------------------------------------------- */\n");
		printf("/* STEP 7: Collect data for pattern analysis         */\n");
		printf("/* ------------------------------------------------- */\n");

		SeqClust_CollectData_ForAutoSeg(strOutputPath, strOutputFile, 
				nFileNum, vFileName, pFileNormFactor, dSmoothFactor,
				nChrNum, vChrName, pChrLen, nBinSize, 
				pIncludeChr);
	}

	/* ---------------------------------------------------------*/
	/* STEP6: K-means clustering                                */
	/*        Perform from k_min to k_max. Choose the best      */
	/*        using BIC.                                        */
	/* ---------------------------------------------------------*/

	/* ---------------------------------------------------------*/
	/* STEP7: Release memory                                    */
	/* ---------------------------------------------------------*/
	DestroyDoubleMatrix(pSmoothK);
	DestroyIntMatrix(pChrLen);
	DestroyDoubleMatrix(pFileReadCount);
	DestroyDoubleMatrix(pFileNormFactor);
	DestroyByteMatrix(pIncludeChr);
	for(ni=0; ni<nChrNum; ni++)
	{
		DeleteString(vChrName[ni]);
		vChrName[ni] = NULL;
	}
	free(vChrName);
	for(ni=0; ni<nFileNum; ni++)
	{
		DeleteString(vFileName[ni]);
		vFileName[ni] = NULL;
	}
	free(vFileName);
	sprintf(strFileName, "%s%s_*.bincount", strOutputPath, strOutputFile);
	RemoveFiles(strFileName);
	sprintf(strFileName, "%s%s_*.tmpreg", strOutputPath, strOutputFile);
	RemoveFiles(strFileName);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClustDP_ComputeBinStats()                                           */
/*  Compute summary statistics for genomic bins.                           */
/* ----------------------------------------------------------------------- */ 
int SeqClustDP_ComputeBinStats(char strOutputPath[], char strOutputFile[], 
			int nFileNum, struct tagString **vFileName, struct DOUBLEMATRIX *pFileNormFactor, 
			int nChrNum, struct tagString **vChrName, struct INTMATRIX *pChrLen,
			int nBinSize, int nMemBlock, int nCorrBlock, double *pBinMax, double *pBinMin,
			int nExportBAR)
{
	/* define */
	int ni;
	int nBinNum;

	/* process chromosome by chromosome */
	for(ni=0; ni<nChrNum; ni++)
	{
		printf(" %s ...", vChrName[ni]->m_pString);
		nBinNum = pChrLen->pMatElement[ni]/nBinSize;
		if(pChrLen->pMatElement[ni] % nBinSize != 0)
			nBinNum++;

		SeqClustDP_ComputeBinStats_Chr(strOutputPath, strOutputFile, 
			nFileNum, vFileName, pFileNormFactor, vChrName[ni]->m_pString, 
			nBinNum, nBinSize, nMemBlock, nCorrBlock, pBinMax, pBinMin, nExportBAR);
	}
	printf("\n");


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClustDP_ComputeBinStats_Chr()                                       */
/*  Compute summary statistics for genomic bins for a single chromosome    */
/* ----------------------------------------------------------------------- */ 
int SeqClustDP_ComputeBinStats_Chr(char strOutputPath[], char strOutputFile[], 
			int nFileNum, struct tagString **vFileName, struct DOUBLEMATRIX *pFileNormFactor, 
			char strChr[], int nBinNum, int nBinSize, int nMemBlock, int nCorrBlock,
			double *pBinMax, double *pBinMin, int nExportBAR)
{
	/* define */
	float *vM = NULL;
	float *vMAM = NULL;
	char strOutFileName[MED_LINE_LENGTH];
	char strInFileName[MED_LINE_LENGTH];
	float **vD = NULL;
	/* struct DOUBLEMATRIX *pLogNormFactor; */
	double dLog2 = log(2.0);

	FILE **vfpIn;
	int ni,nj,nk,nw1,nw2,nu;
	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));

	int nP0,nP1;
	int nQ0,nQ1;
	int nCount,nQCount;
	double dMu,dMAM;
	float *pM,*pMAM;

	FILE *fpRead;
	char strLine[MED_LINE_LENGTH];

	/* init */
	vM = (float *)calloc(nBinNum, sizeof(float));
	if(vM == NULL)
	{
		printf("Error: SeqClustDP_ComputeBinStats_Chr, cannot allocate memory for read count vector!\n");
		exit(EXIT_FAILURE);
	}
	vMAM = (float *)calloc(nBinNum, sizeof(float));
	if(vMAM == NULL)
	{
		printf("Error: SeqClustDP_ComputeBinStats_Chr, cannot allocate memory for read count vector!\n");
		exit(EXIT_FAILURE);
	}

	vfpIn = NULL;
	vfpIn = (FILE **)calloc(nFileNum, sizeof(FILE *));
	if(vfpIn == NULL)
	{
		printf("Error: SeqClustDP_ComputeBinStats_Chr, cannot create vector for file pointers!\n");
		exit(EXIT_FAILURE);
	}

	
	for(ni=0; ni<nFileNum; ni++)
	{
		sprintf(strInFileName, "%s%s_%s_%s.bincount", strOutputPath, strOutputFile, vFileName[ni]->m_pString, strChr);
		vfpIn[ni] = fopen(strInFileName, "rb");
		if(vfpIn[ni] == NULL)
		{
			printf("Error: SeqClustDP_ComputeBinStats_Chr, cannot open file %s!\n", strInFileName);
			exit(EXIT_FAILURE);
		}
	}

	vD = NULL;
	vD = (float **)calloc(nFileNum, sizeof(float *));
	if(vD == NULL)
	{
		printf("Error: SeqClustDP_ComputeBinStats_Chr, cannot create vector for summary statistics!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nFileNum; ni++)
	{
		vD[ni] = (float *)calloc((nMemBlock+nCorrBlock+nCorrBlock), sizeof(float));
		if(vD[ni] == NULL)
		{
			printf("Error: SeqClustDP_ComputeBinStats_Chr, cannot create memory block for computing summary statistics!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* processing */
	nP0 = 0;
	while(nP0 < nBinNum)
	{
		/* determine memory block */
		nP1 = nP0 + nMemBlock - 1;
		if(nP1 >= nBinNum)
			nP1 = nBinNum-1;
		nCount = nP1-nP0+1;
		pM = vM+nP0;
		pMAM = vMAM+nP0;

		nQ0 = nP0-nCorrBlock;
		nQ1 = nP1+nCorrBlock;
		if(nQ0 < 0)
			nQ0 = 0;
		if(nQ1 >= nBinNum)
			nQ1 = nBinNum-1;
		nQCount = nQ1-nQ0+1;

		/* read data */
		for(ni=0; ni<nFileNum; ni++)
		{
			fseek(vfpIn[ni], nQ0*sizeof(float), SEEK_SET);
			if( little_endian_fread(vD[ni], sizeof(float), nQCount, vfpIn[ni], little_endian_machine) != nQCount)
			{
				printf("Error: SeqClustDP_ComputeBinStats_Chr, incorrect loading, number of bins inconsistent!\n");
				exit(EXIT_FAILURE);
			}
		}

		/* process */
		nk = nQ0-nP0;
		for(nj=0; nj<nQCount; nj++)
		{
			dMu = 0.0;
			for(ni=0; ni<nFileNum; ni++)
			{
				vD[ni][nj] = (float)(vD[ni][nj]*(pFileNormFactor->pMatElement[ni]));
		
				if(vD[ni][nj] > dMu)
					dMu = vD[ni][nj];
			}

			*(pM+nk) = (float)dMu;

			nk++;
		}

		nk = nQ0-nP0;
		for(nj=0; nj< nQCount; nj++)
		{
			if( (nk >= 0) && (nk < nCount) )
			{
				dMAM = 0.0;
				
				for(nu=1; nu<=nCorrBlock; nu++)
				{
					nw1 = nj-nu;
					nw2 = nj+nu;
					if( (nw1 < 0) || (nw2 >=nQCount) )
					{
						continue;
					}

					dMAM += (*(pM+nk-nu)) + (*(pM+nk+nu));
				}
				
				dMAM += (*(pM+nk));

				dMAM /= (nCorrBlock+nCorrBlock+1);

				*(pMAM+nk) = (float)dMAM;

				if( dMAM < *pBinMin )
					*pBinMin = dMAM;

				if( dMAM > *pBinMax )
					*pBinMax = dMAM;
			}

			nk++;
		}

		/* get next memory block */
		nP0 = nP1+1;
	}

	
	/* close files and free memory blocks */
	for(ni=0; ni<nFileNum; ni++)
	{
		fclose(vfpIn[ni]);
		vfpIn[ni] = NULL;
		free(vD[ni]);
	}
	free(vfpIn);
	free(vD);

	/* save data to files */
	sprintf(strOutFileName, "%s%s_%s.mam", strOutputPath, strOutputFile, strChr);
	TileMapv2_SaveToBinaryFile((void *)vMAM, sizeof(float), nBinNum, strOutFileName);

	/* FOR DEBUG & EXPLORE */
	if(nExportBAR == 1)
	{
		sprintf(strOutFileName, "%s%s_%s.mam.bar.txt", strOutputPath, strOutputFile, strChr);
		fpRead = NULL;
		fpRead = fopen(strOutFileName, "w");
		if(fpRead == NULL)
		{
			printf("Error: cannot open output file!\n");
			exit(EXIT_FAILURE);
		}

		fprintf(fpRead, "#chr\tpos\t%s_mam_%s\n", strOutputFile,strChr);
		fprintf(fpRead, "#chr\tpos\t1\n");
		nP0 = 12;
		for(ni=0; ni<nBinNum; ni++)
		{
			fprintf(fpRead, "%s\t%d\t%f\n", strChr, nP0, vMAM[ni]);
			nP0 += nBinSize;
		}

		fclose(fpRead);
		
		sprintf(strLine, "%s.cgw", strOutputFile);
		TileMapv2_TXT2BAR(strOutFileName, strOutputPath, strLine);
	}


	/* release memory */
	free(vM);
	free(vMAM);
	
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClustDP_QuantileCutoff()                                            */
/*  Find quantile cutoffs, save to dCL and dCH.                            */
/* ----------------------------------------------------------------------- */ 
int SeqClustDP_QuantileCutoff(char strOutputPath[], char strOutputFile[], 
		int nChrNum, struct tagString **vChrName, struct INTMATRIX *pChrLen, 
		int nBinSize, double dBinMax, double dBinMin, 
		int nGridNum, double dCutL, double dCutH, double *pCL, double *pCH,
		struct BYTEMATRIX *pIncludeChr)
{
	/* define */
	int ni,nj,nk;
	int *vGrid = NULL;
	int nBinNum;
	float *vC;
	char strFileName[MED_LINE_LENGTH];
	double dStep;
	int nTotal = 0;
	int nSum;
	double dPrc;
	int nCLOK = 0;
	int nCHOK = 0;

	/* init */
	dStep = (dBinMax-dBinMin)/nGridNum;

	vGrid = (int *)calloc(nGridNum, sizeof(int));
	if(vGrid == NULL)
	{
		printf("Error: SeqClustDP_QuantileCutoff, cannot allocate memory for grid computation!\n");
		exit(EXIT_FAILURE);
	}

	/* process chromosome by chromosome */
	for(ni=0; ni<nChrNum; ni++)
	{
		/* skip chromosome if necessary */
		if(pIncludeChr->pMatElement[ni] == 0)
			continue;

		/* load data */
		nBinNum = pChrLen->pMatElement[ni]/nBinSize;
		if(pChrLen->pMatElement[ni] % nBinSize != 0)
			nBinNum++;

		nTotal += nBinNum;

		vC = (float *)calloc(nBinNum, sizeof(float));
		if(vC == NULL)
		{
			printf("Error: SeqClustDP_QuantileCutoff, cannot allocate memory for correlation vector!\n");
			exit(EXIT_FAILURE);
		}

		sprintf(strFileName, "%s%s_%s.mam", strOutputPath, strOutputFile, vChrName[ni]->m_pString);

		if( TileMapv2_LoadFromBinaryFile(vC, sizeof(float), nBinNum, strFileName) != PROC_SUCCESS)
		{
			printf("Error: SeqClust_QuantileCutoff, cannot read bin statistics correctly!\n");
			exit(EXIT_FAILURE);
		}

		/* count */
		for(nj=0; nj<nBinNum; nj++)
		{
			nk = (int)((vC[nj]-dBinMin)/dStep);
			if(nk < 0)
				nk = 0;
			if(nk >= nGridNum)
				nk = nGridNum-1;

			vGrid[nk] += 1;
		}

		/* free memory */
		free(vC);
	}

	/* find cutoff */
	nSum = 0;
	for(nj=nGridNum-1; nj>=0; nj--)
	{
		nSum += vGrid[nj];
		dPrc = (double)nSum/(double)nTotal;
		if(dPrc > (1.0 - dCutH))
		{
			if(nCHOK == 0)
			{
				*pCH = dBinMin + dStep*(nj+1);
				nCHOK = 1;
			}
		}

		if(dPrc > (1.0 - dCutL))
		{
			if(nCLOK == 0)
			{
				*pCL = dBinMin + dStep*(nj+1);
				nCLOK = 1;
			}
		}

		if((nCHOK == 1) && (nCLOK ==1))
			break;
	}

	/* free memory */
	free(vGrid);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClustDP_Segmentation()                                              */
/*  Genome segmentation. Divide genome into correlation blocks.            */
/* ----------------------------------------------------------------------- */ 
int SeqClustDP_Segmentation(char strOutputPath[], char strOutputFile[], 
				int nChrNum, struct tagString **vChrName, struct INTMATRIX *pChrLen, 
				int nBinSize, double dCL, double dCH)
{
	/* define */
	int ni,nj;
	int nBinNum;
	float *vC;
	char strFileName[MED_LINE_LENGTH];
	FILE *fpOut;
	char strOutFileName[MED_LINE_LENGTH];
	int nP1,nP2;
	float fMax;
	
	/* process chromosome by chromosome */
	for(ni=0; ni<nChrNum; ni++)
	{
		printf(" %s ... ", vChrName[ni]->m_pString);
		/* load data */
		nBinNum = pChrLen->pMatElement[ni]/nBinSize;
		if(pChrLen->pMatElement[ni] % nBinSize != 0)
			nBinNum++;

		vC = (float *)calloc(nBinNum, sizeof(float));
		if(vC == NULL)
		{
			printf("Error: SeqClustDP_QuantileCutoff, cannot allocate memory for correlation vector!\n");
			exit(EXIT_FAILURE);
		}

		sprintf(strFileName, "%s%s_%s.mam", strOutputPath, strOutputFile, vChrName[ni]->m_pString);

		if( TileMapv2_LoadFromBinaryFile(vC, sizeof(float), nBinNum, strFileName) != PROC_SUCCESS)
		{
			printf("Error: SeqClustDP_QuantileCutoff, cannot read bin statistics correctly!\n");
			exit(EXIT_FAILURE);
		}

		/* open output file */
		sprintf(strOutFileName, "%s%s_%s.tmpreg", strOutputPath, strOutputFile, vChrName[ni]->m_pString);
		fpOut = NULL;
		fpOut = fopen(strOutFileName, "w");
		if(fpOut == NULL)
		{
			printf("Error: SeqClustDP_QuantileCutoff, cannot open file %s to export data!\n", strOutFileName);
			exit(EXIT_FAILURE);
		}

		/* find peaks */
		nP1 = -1;
		nP2 = -1;
		fMax = -1000000.0;
		for(nj=0; nj<nBinNum; nj++)
		{
			if(vC[nj] >= dCL)
			{
				if(nP1 < 0)
				{
					nP1 = nj;
					nP2 = nj;
					fMax = vC[nj];
				}
				else
				{
					nP2 = nj;
					if(vC[nj] > fMax)
						fMax = vC[nj];
				}
			}
			else
			{
				if(nP1 >= 0)
				{
					if(fMax > dCH)
					{
						fprintf(fpOut, "%d\t%d\t%f\n", nP1, nP2, fMax);
					}
					nP1 = -1;
					nP2 = -1;
					fMax = -1000000.0;
				}
			}
		}

		if(nP1 >= 0)
		{
			if(fMax > dCH)
			{
				fprintf(fpOut, "%d\t%d\t%f\n", nP1, nP2, fMax);
			}
			nP1 = -1;
			nP2 = -1;
			fMax = -1000000.0;
		}

		/* free memory */
		free(vC);

		/* close file */
		fclose(fpOut);
	}

	printf("\n");
	
	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  SeqClustDP_FindBlock()                                                 */
/*  Find correlation blocks.                                               */
/* ----------------------------------------------------------------------- */ 
int SeqClustDP_FindBlock(char strOutputPath[], char strOutputFile[], 
				int nFileNum, struct tagString **vFileName, 
				struct DOUBLEMATRIX *pFileNormFactor, double dSmoothFactor,
				int nChrNum, struct tagString **vChrName, 
				struct INTMATRIX *pChrLen, int nBinSize, 
				struct BYTEMATRIX *pIncludeChr, int nBlockLenCut)
{
	/* define */
	int ni,nj;
	FILE *fpOut;
	char strFileName[MED_LINE_LENGTH];
	struct DOUBLEMATRIX *pReg;
	int nP1,nP2;
	int nBCut;

	/* init */
	nBCut = (int)(nBlockLenCut/nBinSize+1);

	/* process chromosome by chromosome */
	for(ni=0; ni<nChrNum; ni++)
	{
		/* skip chromosome if necessary */
		if(pIncludeChr->pMatElement[ni] == 0)
			continue;

		printf(" %s ... ", vChrName[ni]->m_pString);
		/* load data */
		sprintf(strFileName, "%s%s_%s.tmpreg", strOutputPath, strOutputFile, vChrName[ni]->m_pString);
		pReg = NULL;
		pReg = DMLOAD(strFileName);
		if(pReg == NULL)
		{
			printf("  No regions loaded for %s.\n", vChrName[ni]->m_pString);
			continue;
		}
		if(pReg->nHeight <= 0)
		{
			printf("  No regions loaded for %s.\n", vChrName[ni]->m_pString);
			DestroyDoubleMatrix(pReg);
			continue;
		}

		fpOut = NULL;
		fpOut = fopen(strFileName, "w");
		if(fpOut == NULL)
		{
			printf("Error: SeqClustDP_FindBlock, cannot open file to export DP segmentation results.\n");
			exit(EXIT_FAILURE);
		}

		for(nj=0; nj<pReg->nHeight; nj++)
		{
			nP1 = (int)(DMGETAT(pReg, nj, 0));
			nP2 = (int)(DMGETAT(pReg, nj, 1));

			SeqClustDP_FindBlock_SingleRegion_BinStandardize(vChrName[ni]->m_pString, nP1, nP2, nFileNum, 
				strOutputPath, strOutputFile, vFileName, pFileNormFactor, dSmoothFactor, fpOut,
				nBCut);
		}

		/* close file */
		fclose(fpOut);
	}

	printf("\n");

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClustDP_FindBlock_SingleRegion()                                    */
/*  Using dynamic programming to find blocks in a single region.           */
/* ----------------------------------------------------------------------- */ 
int SeqClustDP_FindBlock_SingleRegion(char strChr[], int nStart, int nEnd, int nFileNum, 
				char strOutputPath[], char strOutputFile[], 
				struct tagString **vFileName, struct DOUBLEMATRIX *pFileNormFactor, 
				double dSmoothFactor, FILE *fpOut, int nBlockLenCut)
{
	/* define */
	float **vData;
	int nLen;
	int ni,nj,nk;
	char strInFileName[MED_LINE_LENGTH];
	FILE *fpIn;
	double dLog2 = log(2.0);
	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));
	/* optimal score vector */
	double *vS;
	/* indicator of the optimal path */
	int *vI;
	/* temp vector used for DP */
	double *vT;
	double dScore,dBestScore;
	int nBestI,nCurrentI;
	/* Block number */
	int nRegionNum;
	int *vP1;
	int *vP2;
	double *vPS;

	double *vMu;
	double *vVar;

	/* Load Data */
	nLen = nEnd-nStart+1;

	vData = NULL;
	vData = (float **)calloc(nFileNum, sizeof(float *));
	if(vData == NULL)
	{
		printf("Error: SeqClustDP_FindBlock_SingleRegion, cannot allocate memory for loading data\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nFileNum; ni++)
	{
		vData[ni] = (float *)calloc(nLen, sizeof(float));
		if(vData[ni] == NULL)
		{
			printf("Error: SeqClustDP_FindBlock_SingleRegion, cannot allocate memory for loading data\n");
			exit(EXIT_FAILURE);
		}

		sprintf(strInFileName, "%s%s_%s_%s.bincount", strOutputPath, strOutputFile, vFileName[ni]->m_pString, strChr);
		fpIn = NULL;
		fpIn = fopen(strInFileName, "rb");
		if(fpIn == NULL)
		{
			printf("Error: SeqClustDP_FindBlock_SingleRegion, cannot open file %s!\n", strInFileName);
			exit(EXIT_FAILURE);
		}

		fseek(fpIn, nStart*sizeof(float), SEEK_SET);
		if( little_endian_fread(vData[ni], sizeof(float), nLen, fpIn, little_endian_machine) != nLen)
		{
			printf("Error: SeqClustDP_FindBlock_SingleRegion, incorrect loading, number of bins inconsistent!\n");
			exit(EXIT_FAILURE);
		}

		fclose(fpIn);

		for(nj=0; nj<nLen; nj++)
		{
			vData[ni][nj] = (float)(vData[ni][nj]*(pFileNormFactor->pMatElement[ni])/(dSmoothFactor+1e-6));
			vData[ni][nj] = (float)(log(vData[ni][nj]+1.0)/dLog2);
		}
	}

	/* DP */
	vS = NULL;
	vS = (double *)calloc(nLen, sizeof(double));
	vI = NULL;
	vI = (int *)calloc(nLen, sizeof(int));
	vT = NULL;
	vT = (double *)calloc(nLen, sizeof(double));
	vMu = NULL;
	vMu = (double *)calloc(nFileNum, sizeof(double));
	vVar = NULL;
	vVar = (double *)calloc(nFileNum, sizeof(double));

	if( (vS == NULL) || (vI == NULL) || (vT == NULL) || (vMu == NULL) || (vVar == NULL))
	{
		printf("Error: SeqClustDP_FindBlock_SingleRegion, cannot allocate memory for dynamic programming!\n");
		exit(EXIT_FAILURE);
	}

	vS[0] = 1e10;
	vI[0] = 0;

	for(ni=1; ni<nLen; ni++)
	{
		/* compute score */
		for(nk=0; nk<nFileNum; nk++)
		{
			vMu[nk] = vData[nk][ni];
			vVar[nk] = 0.0;
		}
		
		for(nj=ni-1; nj>=0; nj--)
		{
			dScore = SeqClustDP_BlockBIC_Normal(vData, nFileNum, nj, ni, nLen, vMu, vVar);
			
			if(nj == 0)
				vT[nj] = dScore;
			else
				vT[nj] = vS[nj-1]+dScore;
		}

		/* find optimal score */
		dBestScore = vT[0];
		nBestI = 0;
		for(nj=1; nj<ni; nj++)
		{
			if(vT[nj] < dBestScore)
			{
				dBestScore = vT[nj];
				nBestI = nj;
			}
		}

		vS[ni] = dBestScore;
		vI[ni] = nBestI;
	}

	/* Find Region Number */
	nRegionNum = 0;
	nBestI = vI[nLen-1];
	while(nBestI > 0)
	{
		nRegionNum++;
		nBestI = vI[nBestI-1];
	}
	nRegionNum++;
	
	vP1 = NULL;
	vP2 = NULL;
	vPS = NULL;
	vP1 = (int *)calloc(nRegionNum, sizeof(int));
	vP2 = (int *)calloc(nRegionNum, sizeof(int));
	vPS = (double *)calloc(nRegionNum, sizeof(double));
	if( (vP1 == NULL) || (vP2 == NULL) || (vPS == NULL) )
	{
		printf("Error: SeqClustDP_FindBlock_SingleRegion, cannot allocate memory for back tracing blocks!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	nCurrentI = nLen-1;
	nBestI = vI[nCurrentI];
	while(nBestI > 0)
	{
		vP1[ni] = nBestI;
		vP2[ni] = nCurrentI;
		vPS[ni] = vS[nCurrentI]-vS[nBestI-1];

		ni++;
		nCurrentI = nBestI-1;
		nBestI = vI[nCurrentI];
	}

	vP1[ni] = nBestI;
	vP2[ni] = nCurrentI;
	vPS[ni] = vS[nCurrentI];
	ni++;

	if(ni != nRegionNum)
	{
		printf("Error: SeqClustDP_FindBlock_SingleRegion, region number not match!\n");
		exit(EXIT_FAILURE);
	}
 
	/* Output Regions */
	for(ni=nRegionNum-1; ni>=0; ni--)
	{
		if((vP2[ni]-vP1[ni]+1) >= nBlockLenCut)
			fprintf(fpOut, "%d\t%d\t%f\n", vP1[ni]+nStart, vP2[ni]+nStart, vPS[ni]);
	}

	/* release memory */
	free(vMu);
	free(vVar);
	free(vP1);
	free(vP2);
	free(vPS);
	free(vS);
	free(vI);
	free(vT);
	for(ni=0; ni<nFileNum; ni++)
	{
		free(vData[ni]);
	}
	free(vData);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClustDP_FindBlock_SingleRegion_BinStandardize()                     */
/*  Using dynamic programming to find blocks in a single region.           */
/*  Standardize bin first.                                                 */
/* ----------------------------------------------------------------------- */ 
int SeqClustDP_FindBlock_SingleRegion_BinStandardize(char strChr[], int nStart, int nEnd, int nFileNum, 
				char strOutputPath[], char strOutputFile[], 
				struct tagString **vFileName, struct DOUBLEMATRIX *pFileNormFactor, 
				double dSmoothFactor, FILE *fpOut, int nBlockLenCut)
{
	/* define */
	float **vData;
	int nLen;
	int ni,nj,nk;
	char strInFileName[MED_LINE_LENGTH];
	FILE *fpIn;
	double dLog2 = log(2.0);
	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));
	/* optimal score vector */
	double *vS;
	/* indicator of the optimal path */
	int *vI;
	/* temp vector used for DP */
	double *vT;
	double dScore,dBestScore;
	int nBestI,nCurrentI;
	/* Block number */
	int nRegionNum;
	int *vP1;
	int *vP2;
	double *vPS;

	double *vMu;
	double *vVar;
	double dTemp,dTemp1;

	/* Load Data */
	if(nFileNum <= 1)
	{
		printf("Error: SeqClustDP_FindBlock_SingleRegion, file number must be >1 to use the correlation based segmentation!\n");
		exit(EXIT_FAILURE);
	}


	nLen = nEnd-nStart+1;

	vData = NULL;
	vData = (float **)calloc(nFileNum, sizeof(float *));
	if(vData == NULL)
	{
		printf("Error: SeqClustDP_FindBlock_SingleRegion, cannot allocate memory for loading data\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nFileNum; ni++)
	{
		vData[ni] = (float *)calloc(nLen, sizeof(float));
		if(vData[ni] == NULL)
		{
			printf("Error: SeqClustDP_FindBlock_SingleRegion, cannot allocate memory for loading data\n");
			exit(EXIT_FAILURE);
		}

		sprintf(strInFileName, "%s%s_%s_%s.bincount", strOutputPath, strOutputFile, vFileName[ni]->m_pString, strChr);
		fpIn = NULL;
		fpIn = fopen(strInFileName, "rb");
		if(fpIn == NULL)
		{
			printf("Error: SeqClustDP_FindBlock_SingleRegion, cannot open file %s!\n", strInFileName);
			exit(EXIT_FAILURE);
		}

		fseek(fpIn, nStart*sizeof(float), SEEK_SET);
		if( little_endian_fread(vData[ni], sizeof(float), nLen, fpIn, little_endian_machine) != nLen)
		{
			printf("Error: SeqClustDP_FindBlock_SingleRegion, incorrect loading, number of bins inconsistent!\n");
			exit(EXIT_FAILURE);
		}

		fclose(fpIn);

		for(nj=0; nj<nLen; nj++)
		{
			vData[ni][nj] = (float)(vData[ni][nj]*(pFileNormFactor->pMatElement[ni])/(dSmoothFactor+1e-6));
			vData[ni][nj] = (float)(log(vData[ni][nj]+1.0)/dLog2);
		}
	}

	for(nj=0; nj<nLen; nj++)
	{
		dTemp = 0.0;
		dTemp1 = 0.0;
		for(ni=0; ni<nFileNum; ni++)
		{
			dTemp += vData[ni][nj];
			dTemp1 += vData[ni][nj]*vData[ni][nj];
		}
		dTemp /= nFileNum;
		dTemp1 = dTemp1-nFileNum*dTemp*dTemp;

		dTemp1 /= (nFileNum-1);
		dTemp1 = sqrt(dTemp1)+1e-6;
		
		for(ni=0; ni<nFileNum; ni++)
		{
			vData[ni][nj] = (float)((vData[ni][nj]-dTemp)/dTemp1);
		}
	}

	/* DP */
	vS = NULL;
	vS = (double *)calloc(nLen, sizeof(double));
	vI = NULL;
	vI = (int *)calloc(nLen, sizeof(int));
	vT = NULL;
	vT = (double *)calloc(nLen, sizeof(double));
	vMu = NULL;
	vMu = (double *)calloc(nFileNum-1, sizeof(double));
	vVar = NULL;
	vVar = (double *)calloc(nFileNum-1, sizeof(double));

	if( (vS == NULL) || (vI == NULL) || (vT == NULL) || (vMu == NULL) || (vVar == NULL))
	{
		printf("Error: SeqClustDP_FindBlock_SingleRegion, cannot allocate memory for dynamic programming!\n");
		exit(EXIT_FAILURE);
	}

	vS[0] = 1e10;
	vI[0] = 0;

	for(ni=1; ni<nLen; ni++)
	{
		/* compute score */
		for(nk=0; nk<nFileNum-1; nk++)
		{
			vMu[nk] = vData[nk][ni];
			vVar[nk] = 0.0;
		}
		
		for(nj=ni-1; nj>=0; nj--)
		{
			dScore = SeqClustDP_BlockBIC_Normal(vData, nFileNum-1, nj, ni, nLen, vMu, vVar);
			
			if(nj == 0)
				vT[nj] = dScore;
			else
				vT[nj] = vS[nj-1]+dScore;
		}

		/* find optimal score */
		dBestScore = vT[0];
		nBestI = 0;
		for(nj=1; nj<ni; nj++)
		{
			if(vT[nj] < dBestScore)
			{
				dBestScore = vT[nj];
				nBestI = nj;
			}
		}

		vS[ni] = dBestScore;
		vI[ni] = nBestI;
	}

	/* Find Region Number */
	nRegionNum = 0;
	nBestI = vI[nLen-1];
	while(nBestI > 0)
	{
		nRegionNum++;
		nBestI = vI[nBestI-1];
	}
	nRegionNum++;
	
	vP1 = NULL;
	vP2 = NULL;
	vPS = NULL;
	vP1 = (int *)calloc(nRegionNum, sizeof(int));
	vP2 = (int *)calloc(nRegionNum, sizeof(int));
	vPS = (double *)calloc(nRegionNum, sizeof(double));
	if( (vP1 == NULL) || (vP2 == NULL) || (vPS == NULL) )
	{
		printf("Error: SeqClustDP_FindBlock_SingleRegion, cannot allocate memory for back tracing blocks!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	nCurrentI = nLen-1;
	nBestI = vI[nCurrentI];
	while(nBestI > 0)
	{
		vP1[ni] = nBestI;
		vP2[ni] = nCurrentI;
		vPS[ni] = vS[nCurrentI]-vS[nBestI-1];

		ni++;
		nCurrentI = nBestI-1;
		nBestI = vI[nCurrentI];
	}

	vP1[ni] = nBestI;
	vP2[ni] = nCurrentI;
	vPS[ni] = vS[nCurrentI];
	ni++;

	if(ni != nRegionNum)
	{
		printf("Error: SeqClustDP_FindBlock_SingleRegion, region number not match!\n");
		exit(EXIT_FAILURE);
	}
 
	/* Output Regions */
	for(ni=nRegionNum-1; ni>=0; ni--)
	{
		if((vP2[ni]-vP1[ni]+1) >= nBlockLenCut)
			fprintf(fpOut, "%d\t%d\t%f\n", vP1[ni]+nStart, vP2[ni]+nStart, vPS[ni]);
	}

	/* release memory */
	free(vMu);
	free(vVar);
	free(vP1);
	free(vP2);
	free(vPS);
	free(vS);
	free(vI);
	free(vT);
	for(ni=0; ni<nFileNum; ni++)
	{
		free(vData[ni]);
	}
	free(vData);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqClustDP_BlockBIC_Normal()                                           */
/*  Compute BIC for a block using independent normal.                      */
/* ----------------------------------------------------------------------- */ 
double SeqClustDP_BlockBIC_Normal(float **vData, int nFileNum, int nStart, int nEnd, 
				int nRegionNum, double *vMu, double *vVar)
{
	/* define */
	double dScore = 0.0;
	int ni,nj;
	int n1 = 1;
	int n2 = nEnd-nStart;
	double *vNewMu;
	double dTemp1,dTemp2,dVarShrink;
	double dpi = 4.0*atan(1.0);
	double dNormConst = -0.5*log(2.0*dpi);

	/* update parameters */
	vNewMu = NULL;
	vNewMu = (double *)calloc(nFileNum, sizeof(double));
	if(vNewMu == NULL)
	{
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nFileNum; ni++)
	{
		vNewMu[ni] = (vMu[ni]*n2+vData[ni][nStart])/(n2+1);
	}
	for(ni=0; ni<nFileNum; ni++)
	{
		dTemp1 = vData[ni][nStart]-vNewMu[ni];
		dTemp2 = vMu[ni]-vNewMu[ni];
		vVar[ni] = ((n2-1)*vVar[ni]+dTemp1*dTemp1+n2*dTemp2*dTemp2)/(n1+n2-1);
		dVarShrink = (4*1+(n1+n2-1)*vVar[ni])/(n1+n2+4+2);
		vMu[ni] = vNewMu[ni];
		/* vNewMu[ni] = -0.5*log(vVar[ni]); */
		vNewMu[ni] = -0.5*log(dVarShrink);
	}
	
	/* compute -2*loglike */
	for(ni=nStart; ni<=nEnd; ni++)
	{
		for(nj=0; nj<nFileNum; nj++)
		{
			dTemp1 = vData[nj][ni]-vMu[nj];
			/* dScore += dNormConst + vNewMu[nj] - 0.5*dTemp1*dTemp1/vVar[nj]; */
			dScore += dNormConst + vNewMu[nj] - 0.5*dTemp1*dTemp1/dVarShrink;
		}
	}
	dScore = -2.0*dScore;

	/* add penalty */
	dScore += 2*nFileNum*log((double)nRegionNum); 

	/* release memory */
	free(vNewMu);

	/* return */
	return dScore;
}

/* ----------------------------------------------------------------------- */ 
/*  HTS_Aln2WinBAR()                                                       */
/*  Convert alignment to bar file after read extension.                    */
/* ----------------------------------------------------------------------- */ 
int HTS_Aln2WinBAR(char strTXTFile[], char strOutputFolder[], char strBARHeader[], 
				   int nExtLen, int nBinSize, int nExportStrandProfile,
				   char strSpecies[], char strChrLenFile[])
{
	/* define */
	int nChrNum;
	struct INTMATRIX *pChrLen = NULL;
	char strHeader[MED_LINE_LENGTH];
	int nExportType;
	
	/* get chrnum, chrlen */
	AdjustDirectoryPath(strOutputFolder);
	pChrLen = IMLOAD(strChrLenFile);
	if(pChrLen == NULL)
	{
		printf("Error: HTS_Aln2WinBAR, cannot load chromosome length!\n");
		exit(EXIT_FAILURE);
	}
	nChrNum = pChrLen->nHeight;
	if(nChrNum == 0)
	{
		printf("Warning: HTS_Aln2WinBAR, no chromosome available!\n");
		return PROC_SUCCESS;
	}
	
	/* export combined */
	nExportType = 2;
	HTS_Aln2WinBAR_Export(strTXTFile, strOutputFolder, strBARHeader, 
				   nExtLen, nBinSize, strSpecies, nChrNum, pChrLen, nExportType);
	
	if(nExportStrandProfile == 1)
	{
		/* export forward */
		nExportType = 0;
		sprintf(strHeader, "%s_F", strBARHeader);
		HTS_Aln2WinBAR_Export(strTXTFile, strOutputFolder, strHeader, 
				   nExtLen, nBinSize, strSpecies, nChrNum, pChrLen, nExportType);

		/* export reverse */
		nExportType = 1;
		sprintf(strHeader, "%s_R", strBARHeader);
		HTS_Aln2WinBAR_Export(strTXTFile, strOutputFolder, strHeader, 
				   nExtLen, nBinSize, strSpecies, nChrNum, pChrLen, nExportType);
	}

	/* release memory */
	DestroyIntMatrix(pChrLen);

	/* return */
	return PROC_SUCCESS;
}

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
				 int nExportType)
{
	/* define */
	float **vBinC = NULL;
	int ni,nx;
	int *vBinNum;
	int nBinNum;
	int nHalfExt;

	char strOutFileName[MED_LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	FILE *fpRead;
	char strChr[LINE_LENGTH];
	int nPos;
	int nStrand;
	int nChrId,nBinId,nEndBinId;
	char strCGWFileName[MED_LINE_LENGTH];

	/* init */
	nHalfExt = nExtLen/2;

	vBinC = NULL;
	vBinC = (float **)calloc(nChrNum, sizeof(float *));
	if(vBinC == NULL)
	{
		printf("Error: HTS_Aln2WinBAR_Export, cannot create vector for genomic bins!\n");
		exit(EXIT_FAILURE);
	}

	vBinNum = NULL;
	vBinNum = (int *)calloc(nChrNum, sizeof(int));
	if(vBinNum == NULL)
	{
		printf("Error: HTS_Aln2WinBAR_Export, cannot create vector for genomic bin size!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nChrNum; ni++)
	{
		nBinNum = pChrLen->pMatElement[ni]/nBinSize;
		if(pChrLen->pMatElement[ni] % nBinSize != 0)
			nBinNum++;

		vBinC[ni] = (float *)calloc(nBinNum, sizeof(float));
		if(vBinC[ni] == NULL)
		{
			printf("Error: HTS_Aln2WinBAR_Export, insufficient memory for creating genomic bins, try a larger bin size!\n");
			exit(EXIT_FAILURE);
		}		

		vBinNum[ni] = nBinNum;
	}

	/* count */
	fpRead = NULL;
	fpRead = fopen(strTXTFile, "r");
	if(fpRead == NULL)
	{
		printf("Error: HTS_Aln2WinBAR_Export, cannot open read alignment file!\n");
		exit(EXIT_FAILURE);
	}
	while(fgets(strLine, LONG_LINE_LENGTH, fpRead) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		SeqClust_ReadPos_From_Aln(strLine, strChr, &nPos, &nStrand);

		/* select which read to export */
		if( (nExportType == 0) && (nStrand == 1) )
			continue;
		if( (nExportType == 1) && (nStrand == 0) )
			continue;

		/* find chromosome and update chromosome length */
		nChrId = Genome_ChromosomeName_To_Index(strChr, strSpecies)-1;
		if( (nChrId < 0) || (nChrId >= nChrNum) )
		{
			printf("Warning: HTS_Aln2WinBAR_Export, cannot find the matching chromosome for %s:%d!\n",strChr, nPos);
			continue;
		}

		/* update combined bin count */
		nBinId = nPos/nBinSize;
		if( (nBinId < 0 ) || (nBinId >= vBinNum[nChrId]) )
		{
			printf("Error: HTS_Aln2WinBAR_Export, inconsistent bin number!\n");
			exit(EXIT_FAILURE);
		}

		if(nExportType != 2)
		{
			nBinId = (nPos-nHalfExt)/nBinSize;
			nEndBinId = (nPos+nHalfExt)/nBinSize;
		}
		else
		{
			if(nStrand == 1)
			{
				nBinId = (nPos-nExtLen)/nBinSize;
				nEndBinId = nPos/nBinSize;
			}
			else
			{
				nBinId = nPos/nBinSize;
				nEndBinId = (nPos+nExtLen)/nBinSize;
			}
		}

		if(nBinId < 0)
			nBinId = 0;
		if(nBinId >= vBinNum[nChrId])
			nBinId = vBinNum[nChrId]-1;

		if(nEndBinId < 0)
			nEndBinId = 0;
		if(nEndBinId >= vBinNum[nChrId])
			nEndBinId = vBinNum[nChrId]-1;
		
		while(nBinId <= nEndBinId)
		{
			vBinC[nChrId][nBinId] += 1.0;
			nBinId++;
		}
	}	

	/* close alignment file */
	fclose(fpRead);

	/* save & release memory */
	sprintf(strOutFileName, "%s%s.bar.txt", strOutputFolder, strBARHeader);
	fpRead = NULL;
	fpRead = fopen(strOutFileName, "w");
	if(fpRead == NULL)
	{
		printf("Error: cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpRead, "#chr\tpos\t%s\n", strBARHeader);
	fprintf(fpRead, "#chr\tpos\t1\n");
	for(ni=0; ni<nChrNum; ni++)
	{
		Genome_Index_To_ChromosomeName(strChr, strSpecies, ni+1);
		nPos = (int)(nBinSize/2);
		for(nx=0; nx<vBinNum[ni]; nx++)
		{
			if(vBinC[ni][nx] > 0.5)
				fprintf(fpRead, "%s\t%d\t%f\n", strChr, nPos, vBinC[ni][nx]);
			nPos += nBinSize;
		}
	}

	fclose(fpRead);

	sprintf(strCGWFileName, "%s.cgw", strBARHeader);
	TileMapv2_TXT2BAR(strOutFileName, strOutputFolder, strCGWFileName);
	RemoveFiles(strOutFileName);
	sprintf(strLine, "%s%s", strOutputFolder, strCGWFileName);
	RemoveFiles(strLine);

	/* release memory */
	for(ni=0; ni<nChrNum; ni++)
	{
		free(vBinC[ni]);
		vBinC[ni] = NULL;
	}
	free(vBinC);
	free(vBinNum);

	/* return */
	return PROC_SUCCESS;
}

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
				  int nLFCAdj, int nLFCWin)
{
	/* define */
	int nResult;
	int ni;
	/* char strFileName[MED_LINE_LENGTH]; */

	/* For Peak Calling */
	int nIPNum = 0;
	int nCTNum = 0;
	int nFileNum = 0;
	int nChrNum = 0;
	struct tagString **vChrName = NULL;
	struct tagString **vFileName = NULL;
	int *vGroupId = NULL;
	struct INTMATRIX *pChrLen = NULL;
	struct DOUBLEMATRIX *pFileReadCount = NULL;
	struct BYTEMATRIX *pIncludeChr = NULL;
	struct DOUBLEMATRIX *pFileNormFactor = NULL;
	struct INTMATRIX *pPeakBoundary = NULL;
	struct DOUBLEMATRIX *pPeakRawData = NULL;
	int nSmoothStart = 0;
	int nSmoothNum = 0;
	double dSmoothFactor = 0.0;
	int nBaseLineFileId = 0;
	
	/* Peaks */
	struct DOUBLEMATRIX *pRegion = NULL;
	struct DOUBLEMATRIX *pRegionSort = NULL;
	struct LONGMATRIX *pRegionSid = NULL;
	struct DOUBLEMATRIX *pRegionFDR = NULL;


	/* ---------------------------------------------------------*/
	/* STEP1: First scan, find chromosome number and length     */
	/*        Find read number for each sample                  */
	/* ---------------------------------------------------------*/
	printf("/* ------------------------------------------------- */\n");
	printf("/* STEP 1: Initialize                                */\n");
	printf("/* ------------------------------------------------- */\n");
	rand_u_init(13);
	AdjustDirectoryPath(strOutputPath);
	dSmoothFactor = 1.0+(double)nExtLen/(double)nBinSize;
	
	nResult = SeqPeak_Initialize(strInputPath, &nIPNum, &nCTNum, &nFileNum, &vFileName, &vGroupId, 
			&pFileReadCount, &pFileNormFactor, &nBaseLineFileId, &nChrNum, &vChrName, &pChrLen);

	if( nIPNum <= 0 )
	{
		printf("Error: SeqPeak_Main, there is no IP samples!\n");
		exit(EXIT_FAILURE);
	}

	/* ---------------------------------------------------------*/
	/* STEP2: Count reads for genomic bins                      */
	/* ---------------------------------------------------------*/
	printf("\n");
	printf("/* ------------------------------------------------- */\n");
	printf("/* STEP 2: Count reads for genomic bins              */\n");
	printf("/* ------------------------------------------------- */\n");
	pIncludeChr = NULL;
	pIncludeChr = CreateByteMatrix(1, nChrNum);
	if(pIncludeChr == NULL)
	{
		printf("Error: SeqPeak_Main, cannot create chromosome inclusion vector!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nChrNum; ni++)
		pIncludeChr->pMatElement[ni] = 1;
	nResult = SeqPeak_CountBinReads(strInputPath, strOutputPath,
		strOutputFile, nFileNum, nChrNum, vChrName, pChrLen, nBinSize,
		nExtLen, pIncludeChr, nExportBAR);


	/* ---------------------------------------------------------*/
	/* STEP3: Compute summary statistic for genomic bins        */
	/* ---------------------------------------------------------*/
	printf("\n");
	printf("/* ------------------------------------------------- */\n");
	printf("/* STEP 3: Compute summary statistics for bins       */\n");
	printf("/* ------------------------------------------------- */\n");
	nResult = SeqPeak_ComputeBinStats(strOutputPath, strOutputFile, 
		nFileNum, nIPNum, nCTNum, vFileName, vGroupId, pFileNormFactor, nChrNum, vChrName, pChrLen,
		nBinSize, nWinSize, dSmoothFactor, pFileReadCount, nTStandardize, nPoisFilter, nPoisWin,
		nLFCAdj, nLFCWin);
	
	/* ---------------------------------------------------------*/
	/* STEP4: Peak Calling                                      */
	/* ---------------------------------------------------------*/
	printf("\n");
	printf("/* ------------------------------------------------- */\n");
	printf("/* STEP 4: Find peaks                                */\n");
	printf("/* ------------------------------------------------- */\n");
	
	/* find positive peaks */
	SeqPeak_Segmentation(strOutputPath, strOutputFile, 
			nChrNum, vChrName, pChrLen, nBinSize, dCutoff, nMaxGap, nMinLen, 0, nPoisFilter, 1, dPoisCut);

	if( nCTNum > 0 )
	{
		/* find negative peaks */
		SeqPeak_Segmentation(strOutputPath, strOutputFile, 
				nChrNum, vChrName, pChrLen, nBinSize, dCutoff, nMaxGap, nMinLen, 1, nPoisFilter, 0, dPoisCut);
		
		/* compute peak FDR */
		SeqPeak_FDR_Flip(strOutputPath, strOutputFile, &pRegion, &pRegionSort, 
			&pRegionSid, &pRegionFDR);
	}
	else if( nLFCAdj == 1 )
	{
		/* find negative peaks */
		SeqPeak_Segmentation(strOutputPath, strOutputFile, 
			nChrNum, vChrName, pChrLen, nBinSize, dCutoff, nMaxGap, nMinLen, 1, 0, 1, dPoisCut);

		/* compute peak FDR */
		SeqPeak_FDR_Flip(strOutputPath, strOutputFile, &pRegion, &pRegionSort, 
			&pRegionSid, &pRegionFDR);
	}
	else
	{
		/* compute peak FDR */
		SeqPeak_FDR_Flip(strOutputPath, strOutputFile, &pRegion, &pRegionSort, 
			&pRegionSid, &pRegionFDR);

		for(ni=0; ni<pRegionFDR->nWidth; ni++)
			pRegionFDR->pMatElement[ni] = -1.0;
	}

	/* ---------------------------------------------------------*/
	/* STEP5: Refine Peak Boundaries                            */
	/* ---------------------------------------------------------*/
	printf("\n");
	printf("/* ------------------------------------------------- */\n");
	printf("/* STEP 5: Refine peak boundaries                    */\n");
	printf("/* ------------------------------------------------- */\n");
	if(nBoundaryRefine == 1)
	{
		nResult = SeqPeak_RefineBoundary(pRegion, nBinSize, nBRWin,
				strInputPath, strOutputPath, strOutputFile, 
				nIPNum, nCTNum, nFileNum, vFileName, vGroupId,
				pFileNormFactor, nChrNum, vChrName, pChrLen, nExtLen,
				&pPeakBoundary);
	}

	/* ---------------------------------------------------------*/
	/* STEP6: Collect data for pattern analysis                 */
	/* ---------------------------------------------------------*/
	printf("\n");
	printf("/* ------------------------------------------------- */\n");
	printf("/* STEP 6: Collect data for pattern analysis         */\n");
	printf("/* ------------------------------------------------- */\n");
	if(nCollectRawData == 1)
	{
		nResult = SeqPeak_CollectRawData(pRegion, strOutputPath, strOutputFile, 
				nFileNum, vFileName, pFileNormFactor, dSmoothFactor,
				nChrNum, vChrName, pChrLen, nBinSize, &pPeakRawData);
	}	

	/* ---------------------------------------------------------*/
	/* STEP7: Convert BAR files                                 */
	/* ---------------------------------------------------------*/
	printf("\n");
	printf("/* ------------------------------------------------- */\n");
	printf("/* STEP 7: Export results                            */\n");
	printf("/* ------------------------------------------------- */\n");
	nResult = SeqPeak_ExportPeaks(strOutputPath, strOutputFile, 
		pRegionSort, pRegionSid, pRegionFDR, nBinSize,
		nChrNum, vChrName, nBoundaryRefine, pPeakBoundary,
		nFileNum, vFileName, nCollectRawData, pPeakRawData,
		nPoisFilter, dPoisCut);

	nResult = SeqPeak_CleanFiles(strOutputPath, strOutputFile, 
		nFileNum, vFileName, 
		nChrNum, vChrName, pChrLen, nBinSize,
		nExportBAR, nKeepTempFiles, nPoisFilter);

	/* ---------------------------------------------------------*/
	/* STEP8: Release memory                                    */
	/* ---------------------------------------------------------*/
	DestroyDoubleMatrix(pRegion);
	DestroyDoubleMatrix(pRegionSort);
	DestroyLongMatrix(pRegionSid);
	DestroyDoubleMatrix(pRegionFDR);
	DestroyIntMatrix(pPeakBoundary);
	DestroyDoubleMatrix(pPeakRawData);

	DestroyIntMatrix(pChrLen);
	DestroyDoubleMatrix(pFileReadCount);
	DestroyDoubleMatrix(pFileNormFactor);
	DestroyByteMatrix(pIncludeChr);
	for(ni=0; ni<nChrNum; ni++)
	{
		DeleteString(vChrName[ni]);
		vChrName[ni] = NULL;
	}
	free(vChrName);
	for(ni=0; ni<nFileNum; ni++)
	{
		DeleteString(vFileName[ni]);
		vFileName[ni] = NULL;
	}
	free(vFileName);
	free(vGroupId);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqPeak_Initialize()                                                   */
/*  SeqPeak peak detection: Initialize                                     */
/* ----------------------------------------------------------------------- */ 
int SeqPeak_Initialize(char strInputPath[], int *pIPNum, int *pCTNum, int *pFileNum, 
						struct tagString ***vFileName, int **vGroupId,
						struct DOUBLEMATRIX **ppFileReadCount, 
						struct DOUBLEMATRIX **ppFileNormFactor, int *pBaseLineFileId, 
						int *pChrNum, 
						struct tagString ***vvChrName, struct INTMATRIX **ppChrLen)
{
	/* define */
	int nMaxChrNum = 65535;	
	int ni,nId,nk;
	int nChrNum = 0;
	int nIPNum = 0;
	int nCTNum = 0;
	int nFileNum = 0;
	int nReadCount;
	double dBaseCount;
	struct tagProbeGenomeInfo **vChrList = NULL;
	
	FILE *fpIn;
	FILE *fpRead;
	char strLine[LONG_LINE_LENGTH];
	struct DOUBLEMATRIX *pSampleReadCount;
	struct DOUBLEMATRIX *pSampleReadCountSort;
	struct LONGMATRIX *pSampleReadCountIdx;
	char strChr[LINE_LENGTH];
	char strFileName[MED_LINE_LENGTH];
	int nPos;
	int nStrand;
	char *chp;

	/* Count file number */
	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: SeqPeak_Initialize, cannot open input file list!\n");
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

		nFileNum++;

		chp = strchr(strLine, '\t');
		if(chp == NULL)
		{
			printf("Error: SeqPeak_Initialize, you need to provide group ids for each sample (1 for IP, 0 for control)!\n");
			exit(EXIT_FAILURE);
		}

		chp++;
		StrTrimLeft(chp);
		if(atoi(chp) == 1)
			nIPNum++;
		else
			nCTNum++;
	}
	fclose(fpIn);
	printf("No. of IP Samples (Files) = %d\n", nIPNum);
	printf("No. of CT Samples (Files) = %d\n", nCTNum);
	printf("No. of Total Samples (Files) = %d\n", nFileNum);

	/* prepare memory for file information */
	*pIPNum = nIPNum;
	*pCTNum = nCTNum;
	*pFileNum = nFileNum;


	*vFileName = NULL;
	*vFileName = (struct tagString **)calloc(nFileNum, sizeof(struct tagString *));
	if(*vFileName == NULL)
	{
		printf("Error: SeqPeak_Initialize, cannot create vector for file names!\n");
		exit(EXIT_FAILURE);
	}

	*vGroupId = NULL;
	*vGroupId = (int *)calloc(nFileNum, sizeof(int));
	if(*vGroupId == NULL)
	{
		printf("Error: SeqPeak_Initialize, cannot create vector for file group ids!\n");
		exit(EXIT_FAILURE);
	}

	*ppFileReadCount = NULL;
	*ppFileReadCount = CreateDoubleMatrix(1, nFileNum);
	if(*ppFileReadCount == NULL)
	{
		printf("Error: SeqPeak_Initialize, cannot create vector for file read counts!\n");
		exit(EXIT_FAILURE);
	}

	pSampleReadCount = NULL;
	pSampleReadCount = CreateDoubleMatrix(1, nFileNum);
	if(pSampleReadCount == NULL)
	{
		printf("Error: SeqPeak_Initialize, cannot create vector for file read counts!\n");
		exit(EXIT_FAILURE);
	}

	*ppFileNormFactor = NULL;
	*ppFileNormFactor = CreateDoubleMatrix(1, nFileNum);
	if(*ppFileNormFactor == NULL)
	{
		printf("Error: SeqPeak_Initialize, cannot create vector for file normalizing factor!\n");
		exit(EXIT_FAILURE);
	}

	/* process files one by one */
	vChrList = NULL;
	vChrList = (struct tagProbeGenomeInfo **)calloc(nMaxChrNum, sizeof(struct tagProbeGenomeInfo *));
	if(vChrList == NULL)
	{
		printf("Error: SeqPeak_Initialize, cannot create vector for chromosomes!\n");
		exit(EXIT_FAILURE);
	}
	nChrNum = 0;

	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: SeqPeak_Initialize, cannot open input file list!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		/* read alignment file */
		chp = strchr(strLine, '\t');
		*chp = '\0';
		chp++;
		StrTrimLeft(chp);
		(*vGroupId)[ni] = atoi(chp);

		GetFileName(strLine, strFileName);
		StringAddTail((*vFileName)+ni, strFileName);

		printf("  Processing %s ...\n", strLine);
		nReadCount = 0;
		fpRead = NULL;
		fpRead = fopen(strLine, "r");
		if(fpRead == NULL)
		{
			printf("Error: SeqPeak_Initialize, cannot open read alignment file!\n");
			exit(EXIT_FAILURE);
		}
		while(fgets(strLine, LONG_LINE_LENGTH, fpRead) != NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;
			if(strLine[0] == '#')
				continue;

			SeqClust_ReadPos_From_Aln(strLine, strChr, &nPos, &nStrand);

			/* find chromosome and update chromosome length */
			nId = SeqClust_FindChrInfo(vChrList, nMaxChrNum, &nChrNum, strChr);
			if( (nId < 0) || (nId >= nMaxChrNum) )
			{
				printf("Error: SeqPeak_Initialize, cannot find the matching chromosome!\n");
				exit(EXIT_FAILURE);
			}
			if(nPos > vChrList[nId]->nPos)
			{
				vChrList[nId]->nPos = nPos;
			}
			
			nReadCount++;
		}	

		/* close alignment file */
		fclose(fpRead);

		/* update file information */
		pSampleReadCount->pMatElement[ni] = nReadCount;
		(*ppFileReadCount)->pMatElement[ni] = nReadCount;

		ni++;
	}
	
	fclose(fpIn);
	
	if(ni != nFileNum)
	{
		printf("Error: SeqPeak_Initialize, inconsistent file number!\n");
		exit(EXIT_FAILURE);
	}

	/* return chromosome name and length */
	*pChrNum = nChrNum;
	*vvChrName = (struct tagString **)calloc(nChrNum, sizeof(struct tagString *));
	if(*vvChrName == NULL)
	{
		printf("Error: SeqPeak_Initialize, cannot create memory for chromosome name!\n");
		exit(EXIT_FAILURE);
	}

	*ppChrLen = NULL;
	*ppChrLen = CreateIntMatrix(1, nChrNum);
	if(*ppChrLen == NULL)
	{
		printf("Error: SeqPeak_Initialize, cannot create memory for chromosome length!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nChrNum; ni++)
	{
		(*vvChrName)[ni] = vChrList[ni]->pProbe;
		vChrList[ni]->pProbe = NULL;
		(*ppChrLen)->pMatElement[ni] = vChrList[ni]->nPos + 1;

		ProbeGenomeInfoDestroy(vChrList+ni);
	}
	free(vChrList);

	if(ni != nChrNum)
	{
		printf("Error: SeqPeak_Initialize, inconsistent chromosome number!\n");
		exit(EXIT_FAILURE);
	}

	/* compute normalizing factor */
	pSampleReadCountSort = NULL;
	pSampleReadCountIdx = NULL;
	DMSORTMERGEA_0(pSampleReadCount, &pSampleReadCountSort, &pSampleReadCountIdx);
	/* if(nFileNum % 2 == 0)
	{
		nk = nFileNum/2 - 1;
	}
	else
	{
		nk = nFileNum/2;
	} */
	nk = 0;
	dBaseCount = pSampleReadCountSort->pMatElement[nk];
	*pBaseLineFileId = pSampleReadCountIdx->pMatElement[nk];
	for(ni=0; ni<nFileNum; ni++)
	{
		if( (*ppFileReadCount)->pMatElement[ni] <= 0)
		{
			printf("Error: SeqPeak_Initialize, sample read count for the %d-th file is zero, cannot perform normalization!\n", ni+1);
			exit(EXIT_FAILURE);
		}
		(*ppFileNormFactor)->pMatElement[ni] = dBaseCount/(*ppFileReadCount)->pMatElement[ni];
	}


	/* print information */
	printf("\nFile\tNo_of_Reads\tScaling_Factor\n");
	for(ni=0; ni<(*pFileNum); ni++)
	{
		printf("%s\t%d\t%f\n", (*vFileName)[ni]->m_pString, (int)((*ppFileReadCount)->pMatElement[ni]), (*ppFileNormFactor)->pMatElement[ni]);
	}

	printf("\nNo. of chromosomes = %d\n", *pChrNum);
	printf("\tChromosome\tMax_Coordinate\n");
	for(ni=0; ni<(*pChrNum); ni++)
	{
		printf("\t%s\t%d\n", (*vvChrName)[ni]->m_pString, (*ppChrLen)->pMatElement[ni]);
	}

	/* release memory */
	DestroyDoubleMatrix(pSampleReadCount);
	DestroyDoubleMatrix(pSampleReadCountSort);
	DestroyLongMatrix(pSampleReadCountIdx);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqPeak_CountBinReads()                                                */
/*  Count reads for genomic bins                                           */
/* ----------------------------------------------------------------------- */ 
int SeqPeak_CountBinReads(char strInputPath[], char strOutputPath[],
		char strOutputFile[], int nFileNum, int nChrNum, 
		struct tagString **vChrName, struct INTMATRIX *pChrLen, int nBinSize,
		int nExtLen, struct BYTEMATRIX *pIncludeChr,
		int nExportBAR)
{
	/* define */
	int ni;
	int nResult;
	FILE *fpIn;
	char strLine[LONG_LINE_LENGTH];
	char *chp;

	/* open file */
	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: SeqPeak_CountBinReads, cannot open input file list!\n");
		exit(EXIT_FAILURE);
	}

	/* process one by one */
	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		chp = strchr(strLine, '\t');
		*chp = '\0';
		StrTrimRight(strLine);

		/* read alignment file */
		printf("  Processing %s ...\n", strLine);
		nResult = SeqPeak_CountBinReads_SingleFile(strLine, strOutputPath,
				strOutputFile, nChrNum, vChrName, pChrLen, nBinSize,
				nExtLen, pIncludeChr, nExportBAR);		
	
		ni++;
	}
	
	/* close file */
	fclose(fpIn);
	
	if(ni != nFileNum)
	{
		printf("Error: SeqPeak_CountBinReads, inconsistent file number!\n");
		exit(EXIT_FAILURE);
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqPeak_CountBinReads_SingleFile()                                     */
/*  Count reads for genomic bins for a single sample file                  */
/* ----------------------------------------------------------------------- */ 
int SeqPeak_CountBinReads_SingleFile(char strInputFile[], char strOutputPath[],
				char strOutputFile[], int nChrNum, struct tagString **vChrName, 
				struct INTMATRIX *pChrLen, int nBinSize, 		
				int nExtLen, struct BYTEMATRIX *pIncludeChr,
				int nExportBAR)
{
	/* define */
	float **vBinC = NULL;
	int ni,nj,nx;
	int *vBinNum;
	int nBinNum;
	char strFileName[LINE_LENGTH];
	char strOutFileName[MED_LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	FILE *fpRead;
	char strChr[LINE_LENGTH];
	int nPos,nEndPos;
	int nStrand;
	int nChrId,nBinId,nEndBinId;

	/* init */
	GetFileName(strInputFile, strFileName);
	vBinC = NULL;
	vBinC = (float **)calloc(nChrNum, sizeof(float *));
	if(vBinC == NULL)
	{
		printf("Error: SeqPeak_CountBinReads_SingleFile, cannot create vector for genomic bins!\n");
		exit(EXIT_FAILURE);
	}

	vBinNum = NULL;
	vBinNum = (int *)calloc(nChrNum, sizeof(int));
	if(vBinNum == NULL)
	{
		printf("Error: SeqPeak_CountBinReads_SingleFile, cannot create vector for genomic bin size!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nChrNum; ni++)
	{
		nBinNum = pChrLen->pMatElement[ni]/nBinSize;
		if(pChrLen->pMatElement[ni] % nBinSize != 0)
			nBinNum++;

		vBinC[ni] = (float *)calloc(nBinNum, sizeof(float));
		if(vBinC[ni] == NULL)
		{
			printf("Error: SeqPeak_CountBinReads_SingleFile, insufficient memory for creating genomic bins, try a larger bin size!\n");
			exit(EXIT_FAILURE);
		}		

		vBinNum[ni] = nBinNum;
	}

	/* count */
	fpRead = NULL;
	fpRead = fopen(strInputFile, "r");
	if(fpRead == NULL)
	{
		printf("Error: SeqPeak_CountBinReads_SingleFile, cannot open read alignment file!\n");
		exit(EXIT_FAILURE);
	}
	while(fgets(strLine, LONG_LINE_LENGTH, fpRead) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		SeqClust_ReadPos_From_Aln(strLine, strChr, &nPos, &nStrand);

		/* find chromosome and update chromosome length */
		nChrId = SeqClust_FindChr(vChrName, nChrNum, strChr);
		if( (nChrId < 0) || (nChrId >= nChrNum) )
		{
			printf("Warning: SeqPeak_CountBinReads_SingleFile, cannot find the matching chromosome for %s:%d!\n",strChr, nPos);
			continue;
		}

		/* update bin count */
		nBinId = nPos/nBinSize;
		if( (nBinId < 0 ) || (nBinId >= vBinNum[nChrId]) )
		{
			printf("Error: SeqPeak_CountBinReads_SingleFile, inconsistent bin number!\n");
			exit(EXIT_FAILURE);
		}

		/* '-' strand */
		if(nStrand == 1)
		{
			nEndPos = nPos-nExtLen;
			nEndBinId = nEndPos/nBinSize;
			if(nEndBinId < 0)
				nEndBinId = 0;
			if(nEndBinId >= vBinNum[nChrId])
				nEndBinId = vBinNum[nChrId]-1;

			while(nBinId >= nEndBinId)
			{
				vBinC[nChrId][nBinId] += 1.0;
				nBinId--;
			}
		}
		/* '+' strand */
		else
		{
			nEndPos = nPos+nExtLen;
			nEndBinId = nEndPos/nBinSize;
			if(nEndBinId < 0)
				nEndBinId = 0;
			if(nEndBinId >= vBinNum[nChrId])
				nEndBinId = vBinNum[nChrId]-1;

			while(nBinId <= nEndBinId)
			{
				vBinC[nChrId][nBinId] += 1.0;
				nBinId++;
			}
		}
	}	

	/* close alignment file */
	fclose(fpRead);

	/* create chromosome indicator */
	for(ni=0; ni<nChrNum; ni++)
	{
		nx = 0;
		for(nj=0; nj<vBinNum[ni]; nj++)
		{
			if(vBinC[ni][nj] > 1e-6)
			{
				nx = 1;
				break;
			}
		}

		if(nx == 0)
		{
			pIncludeChr->pMatElement[ni] = 0;
		}
	}

	/* save & release memory */
	if(nExportBAR == 1)
	{
		sprintf(strOutFileName, "%s%s_%s.bar.txt", strOutputPath, strOutputFile, strFileName);
		fpRead = NULL;
		fpRead = fopen(strOutFileName, "w");
		if(fpRead == NULL)
		{
			printf("Error: cannot open output file!\n");
			exit(EXIT_FAILURE);
		}

		fprintf(fpRead, "#chr\tpos\t%s_binc\n", strFileName);
		fprintf(fpRead, "#chr\tpos\t1\n");

		/* process chromosome by chromosome */
		for(ni=0; ni<nChrNum; ni++)
		{
			nPos = nBinSize/2;
			for(nx=0; nx<vBinNum[ni]; nx++)
			{
				if(vBinC[ni][nx] > 1e-6)
					fprintf(fpRead, "%s\t%d\t%f\n", vChrName[ni]->m_pString, nPos, vBinC[ni][nx]);
				nPos += nBinSize;
			}
		}

		fclose(fpRead);
	}

	for(ni=0; ni<nChrNum; ni++)
	{
		sprintf(strOutFileName, "%s%s_%s_%s.bincount", strOutputPath, strOutputFile, strFileName, vChrName[ni]->m_pString);
		TileMapv2_SaveToBinaryFile((void *)(vBinC[ni]), sizeof(float), vBinNum[ni], strOutFileName);

		free(vBinC[ni]);
		vBinC[ni] = NULL;
	}
	free(vBinC);
	free(vBinNum);

	/* return */
	return PROC_SUCCESS;
}


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
			int nLFCAdj, int nLFCWin)
{
	/* define */
	int ni;
	int nBinNum;
	double dS2M = 0.0;
	double dS2S = 0.0;
	double dTM = 0.0;
	double dTS = 0.0;
	int nTotalBinNum = 0;
	int nEffectBinNum = 0;
	int nDf;
	double dB = 0.0;

	/* process chromosome by chromosome to get initial statistics */
	for(ni=0; ni<nChrNum; ni++)
	{
		nBinNum = pChrLen->pMatElement[ni]/nBinSize;
		if(pChrLen->pMatElement[ni] % nBinSize != 0)
			nBinNum++;

		SeqPeak_ComputeBinStats_Chr_Initial(strOutputPath, strOutputFile, 
			nFileNum, nIPNum, nCTNum, vFileName, vGroupId, pFileNormFactor, vChrName[ni]->m_pString, 
			nBinNum, nBinSize, nWinSize, dSmoothFactor, 
			pFileReadCount, nPoisFilter, nPoisWin, 
			nLFCAdj, nLFCWin,
			&dS2M, &dS2S, &nEffectBinNum);

		nTotalBinNum += nBinNum;
	}

	/* compute shrinkage factor */
	if( (nIPNum > 1) || (nCTNum > 1) )
	{
		dS2M /= nEffectBinNum;
		dS2S -= nEffectBinNum*dS2M*dS2M;
		nDf = 0;
		if(nCTNum > 0)
			nDf += nCTNum-1;
		if(nIPNum > 0)
			nDf += nIPNum-1;

		dB = 2.0*(nEffectBinNum-1)/(nDf+2.0)/nEffectBinNum + 2.0*dS2M*dS2M*(nEffectBinNum-1)/(nDf+2.0)/dS2S;
		if(dB < 0.0)
			dB = 0.0;
		if(dB > 1.0)
			dB = 1.0;
	}
	else
	{
		dS2M = 0.0;
		dS2S = 0.0;
		nDf = 0;
		dB = 0.0;
	}

	/* process chromosome by chromosome to get log2 fc and t-statistics with variance shrinkage */
	for(ni=0; ni<nChrNum; ni++)
	{
		nBinNum = pChrLen->pMatElement[ni]/nBinSize;
		if(pChrLen->pMatElement[ni] % nBinSize != 0)
			nBinNum++;
		SeqPeak_ComputeBinStats_Chr_Tstat(strOutputPath, strOutputFile, 
			vChrName[ni]->m_pString, nBinNum, dB, dS2M, nIPNum, nCTNum, &dTM, &dTS);
	}

	/* standardize t-statistics */
	dTM /= nTotalBinNum;
	dTS = (dTS-nTotalBinNum*dTM*dTM)/nTotalBinNum;
	dTS = sqrt(dTS)+1e-6;
	if(nTStandardize == 1)
	{
		for(ni=0; ni<nChrNum; ni++)
		{
			nBinNum = pChrLen->pMatElement[ni]/nBinSize;
			if(pChrLen->pMatElement[ni] % nBinSize != 0)
				nBinNum++;

			SeqPeak_ComputeBinStats_Chr_Tstandard(strOutputPath, strOutputFile, 
				vChrName[ni]->m_pString, nBinNum, dTM, dTS);
		}
		
	}

	/* return */
	return PROC_SUCCESS;
}

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
			double *pS2M, double *pS2S, int *pEffectBinNum)
{
	/* define */
	/* data */
	float **vD = NULL;
	/* window average */
	float **vW = NULL;
	/* group mean */
	float *vM0 = NULL;
	float *vM1 = NULL;
	/* sample variance */
	float *vV = NULL;
	/* log2 fc */
	float *vFC = NULL;
	
	/* total CT read count */
	float *vT0 = NULL;
	/* total IP read count */
	float *vT1 = NULL;
	/* CT poisson p-value */
	float *vP0 = NULL;
	/* IP poisson p-value */
	float *vP1 = NULL;

	char strOutFileName[MED_LINE_LENGTH];
	char strInFileName[MED_LINE_LENGTH];
	int nWinNum;
	float dSum;
	int nPoisWinNum = 0;
	int nPoisWinSize = 0;
	float dSumT = 0.0;
	int nLFCWinNum = 0;
	int nLFCWinSize = 0;
	float dSumL = 0.0;
	double dIPCount,dCTCount; 
	float dTemp;
	int nDf;
	double dLambda1,dLambda2;
	double dTempP,dTempL1,dTempL0;

	/* struct DOUBLEMATRIX *pLogNormFactor; */
	double dLog2 = log(2.0);
	double dLog10 = log(10.0);

	FILE **vfpIn;
	int ni,nj,nk,nl,nx,ny;
	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));

	/* init */
	nWinNum = 2*nWinSize+1;
	if(nPoisFilter == 1)
	{
		nPoisWinSize = nPoisWin/nBinSize/2;
		if(nPoisWinSize < 1)
			nPoisWinSize = 1;
		nPoisWinNum = 2*nPoisWinSize+1;
	}
	if(nLFCAdj == 1)
	{
		nLFCWinSize = nLFCWin/nBinSize/2;
		if(nLFCWinSize <= nWinSize)
			nLFCWinSize = nWinSize+1;
		nLFCWinNum = 2*nLFCWinSize+1;
	}

	vD = NULL;
	vD = (float **)calloc(nFileNum, sizeof(float *));
	if(vD == NULL)
	{
		printf("Error: SeqPeak_ComputeBinStats_Chr_Initial, cannot create vector for data!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nFileNum; ni++)
	{
		vD[ni] = (float *)calloc(nBinNum, sizeof(float));
		if(vD[ni] == NULL)
		{
			printf("Error: SeqPeak_ComputeBinStats_Chr_Initial, cannot create memory block for data!\n");
			exit(EXIT_FAILURE);
		}
	}

	vW = NULL;
	vW = (float **)calloc(nFileNum, sizeof(float *));
	if(vW == NULL)
	{
		printf("Error: SeqPeak_ComputeBinStats_Chr_Initial, cannot create vector for window statistics!\n");
		exit(EXIT_FAILURE);
	}

	vM0 = (float *)calloc(nBinNum, sizeof(float));
	if(vM0 == NULL)
	{
		printf("Error: SeqPeak_ComputeBinStats_Chr_Initial, cannot allocate memory for M0!\n");
		exit(EXIT_FAILURE);
	}

	vM1 = (float *)calloc(nBinNum, sizeof(float));
	if(vM1 == NULL)
	{
		printf("Error: SeqPeak_ComputeBinStats_Chr_Initial, cannot allocate memory for M1!\n");
		exit(EXIT_FAILURE);
	}

	vV = (float *)calloc(nBinNum, sizeof(float));
	if(vV == NULL)
	{
		printf("Error: SeqPeak_ComputeBinStats_Chr_Initial, cannot allocate memory for V!\n");
		exit(EXIT_FAILURE);
	}

	vFC = (float *)calloc(nBinNum, sizeof(float));
	if(vFC == NULL)
	{
		printf("Error: SeqPeak_ComputeBinStats_Chr_Initial, cannot allocate memory for FC!\n");
		exit(EXIT_FAILURE);
	}

	if(nPoisFilter == 1)
	{
		vT0 = (float *)calloc(nBinNum, sizeof(float));
		if(vT0 == NULL)
		{
			printf("Error: SeqPeak_ComputeBinStats_Chr_Initial, cannot allocate memory for T0!\n");
			exit(EXIT_FAILURE);
		}

		vT1 = (float *)calloc(nBinNum, sizeof(float));
		if(vT1 == NULL)
		{
			printf("Error: SeqPeak_ComputeBinStats_Chr_Initial, cannot allocate memory for T1!\n");
			exit(EXIT_FAILURE);
		}

		vP0 = (float *)calloc(nBinNum, sizeof(float));
		if(vP0 == NULL)
		{
			printf("Error: SeqPeak_ComputeBinStats_Chr_Initial, cannot allocate memory for P0!\n");
			exit(EXIT_FAILURE);
		}

		vP1 = (float *)calloc(nBinNum, sizeof(float));
		if(vP1 == NULL)
		{
			printf("Error: SeqPeak_ComputeBinStats_Chr_Initial, cannot allocate memory for P1!\n");
			exit(EXIT_FAILURE);
		}
	}

	vfpIn = NULL;
	vfpIn = (FILE **)calloc(nFileNum, sizeof(FILE *));
	if(vfpIn == NULL)
	{
		printf("Error: SeqPeak_ComputeBinStats_Chr_Initial, cannot create vector for file pointers!\n");
		exit(EXIT_FAILURE);
	}

	/* load data */
	for(ni=0; ni<nFileNum; ni++)
	{
		sprintf(strInFileName, "%s%s_%s_%s.bincount", strOutputPath, strOutputFile, vFileName[ni]->m_pString, strChr);
		vfpIn[ni] = fopen(strInFileName, "rb");
		if(vfpIn[ni] == NULL)
		{
			printf("Error: SeqPeak_ComputeBinStats_Chr_Initial, cannot open file %s!\n", strInFileName);
			exit(EXIT_FAILURE);
		}
		
		if( little_endian_fread(vD[ni], sizeof(float), nBinNum, vfpIn[ni], little_endian_machine) != nBinNum)
		{
			printf("Error: SeqPeak_ComputeBinStats_Chr_Initial, incorrect loading, number of bins inconsistent!\n");
			exit(EXIT_FAILURE);
		}

		fclose(vfpIn[ni]);
	}
	free(vfpIn);

	/* compute window average */
	dIPCount = 0.0;
	dCTCount = 0.0;
	for(ni=0; ni<nFileNum; ni++)
	{
		vW[ni] = (float *)calloc(nBinNum, sizeof(float));
		if(vW[ni] == NULL)
		{
			printf("Error: SeqPeak_ComputeBinStats_Chr_Initial, cannot create memory block for window statistics!\n");
			exit(EXIT_FAILURE);
		}

		if(nBinNum <= nWinSize)
			continue;

		if(nPoisFilter == 1)
		{
			if(nBinNum <= nPoisWinSize)
				continue;
		}
		if(nLFCAdj == 1)
		{
			if(nBinNum <= nLFCWinSize)
				continue;
		}
		
		if(vGroupId[ni] == 0)
			dCTCount += pFileReadCount->pMatElement[ni];
		else
			dIPCount += pFileReadCount->pMatElement[ni];

		dSum = 0.0;
		for(nj=0; nj<nWinSize; nj++)
		{
			dSum += vD[ni][nj];
		}

		if(nPoisFilter == 1)
		{
			dSumT = 0.0;
			for(nj=0; nj<nPoisWinSize; nj++)
			{
				dSumT += vD[ni][nj];
			}
		}

		if(nLFCAdj == 1)
		{
			dSumL = 0.0;
			for(nj=0; nj<nLFCWinSize; nj++)
			{
				dSumL += vD[ni][nj];
			}
		}

		for(nj=0; nj<nBinNum; nj++)
		{
			nk = nj+nWinSize;
			if(nk < nBinNum)
				dSum = dSum + vD[ni][nk];
			else
				nk = nBinNum-1;

			nl = nj-nWinSize-1;
			if(nl >= 0)
				dSum = dSum - vD[ni][nl];
			else
				nl = -1;
			nx = nk-nl;

			if(nLFCAdj == 0)
			{
				vW[ni][nj] = (float)(log(1.0+dSum*pFileNormFactor->pMatElement[ni]/dSmoothFactor)/dLog2);
			}
			else
			{
				nk = nj+nLFCWinSize;
				if(nk < nBinNum)
					dSumL = dSumL + vD[ni][nk];
				else
					nk = nBinNum-1;

				nl = nj-nLFCWinSize-1;
				if(nl >= 0)
					dSumL = dSumL - vD[ni][nl];
				else
					nl = -1;

				ny = nk-nl;

				dTempL1 = dSum*pFileNormFactor->pMatElement[ni]/dSmoothFactor;
				dTempL0 = (dSumL-dSum)*pFileNormFactor->pMatElement[ni]*nx/dSmoothFactor/(ny-nx);
				
				vW[ni][nj] = (float)(log( (1.0+dTempL1)/(1.0+dTempL0) )/dLog2);
			}

			if(vGroupId[ni] == 0)
				vM0[nj] += vW[ni][nj];
			else
				vM1[nj] += vW[ni][nj];


			/* prepare for poisson filter */
			if(nPoisFilter == 1)
			{
				nk = nj+nPoisWinSize;
				if(nk < nBinNum)
					dSumT = dSumT + vD[ni][nk];
				nk = nj-nPoisWinSize-1;
				if(nk >= 0)
					dSumT = dSumT - vD[ni][nk];

				if(vGroupId[ni] == 0)
				{
					vT0[nj] += vD[ni][nj];
					vP0[nj] += dSumT;
				}
				else
				{
					vT1[nj] += vD[ni][nj];
					vP1[nj] += dSumT;
				}
			}
		}
	}

	/* compute mean and variance */
	nDf = 0;
	if(nCTNum > 0)
		nDf += nCTNum-1;
	if(nIPNum > 0)
		nDf += nIPNum-1;

	for(nj=0; nj<nBinNum; nj++)
	{
		if(nCTNum > 0)
			vM0[nj] /= nCTNum;

		if(nIPNum > 0)
			vM1[nj] /= nIPNum;

		vFC[nj] = vM1[nj]-vM0[nj];

		if(nDf > 0)
		{
			for(ni=0; ni<nFileNum; ni++)
			{
				if(vGroupId[ni] == 0)
				{
					dTemp = vW[ni][nj] - vM0[nj];
				}
				else
				{
					dTemp = vW[ni][nj] - vM1[nj];
				}

				vV[nj] += dTemp*dTemp;
			}

			vV[nj] /= nDf;

			if(vV[nj] > 1e-6)
			{
				*pS2M += vV[nj];
				*pS2S += (vV[nj]*vV[nj]);
				*pEffectBinNum = (*pEffectBinNum)+1;
			}
		}

		if(nPoisFilter == 1)
		{
			nk = nj+nPoisWinSize;
			if(nk >= nBinNum)
				nk = nBinNum-1;

			nl = nj-nPoisWinSize;
			if(nl < 0)
				nl = 0;

			nk = nk-nl+1;

			if(nCTNum > 0)
			{
				dLambda1 = vT0[nj];
				dLambda2 = vP0[nj]/nk;
				if(dLambda1 < dLambda2)
					dLambda1 = dLambda2;
				dLambda1 = dLambda1*dIPCount/dCTCount;
			}
			else
			{
				dLambda1 = (vP1[nj]-vT1[nj])/(nk-1);
			}
			if(vT1[nj] < 1e-6)
				vP1[nj] = 0.0;
			else
			{
				dTempP = gammp(vT1[nj], dLambda1);
				if(dTempP < 1e-100)
					dTempP = 1e-100;
				vP1[nj] = (float)(-log(dTempP)/dLog10);
			}

			if(nCTNum > 0)
			{
				dLambda1 = vT1[nj];
				dLambda2 = vP1[nj]/nk;
				if(dLambda1 < dLambda2)
					dLambda1 = dLambda2;
				dLambda1 = dLambda1*dCTCount/dIPCount;

				if(vT0[nj] < 1e-6)
					vP0[nj] = 0.0;
				else
				{
					dTempP = gammp(vT0[nj], dLambda1);
					if(dTempP < 1e-100)
						dTempP = 1e-100;
					vP0[nj] = (float)(-log(dTempP)/dLog10);
				}
			}
		}
	}

	/* save data to files */
	sprintf(strOutFileName, "%s%s_%s.log2fc", strOutputPath, strOutputFile, strChr);
	TileMapv2_SaveToBinaryFile((void *)vFC, sizeof(float), nBinNum, strOutFileName);

	sprintf(strOutFileName, "%s%s_%s.v", strOutputPath, strOutputFile, strChr);
	TileMapv2_SaveToBinaryFile((void *)vV, sizeof(float), nBinNum, strOutFileName);
	
	if(nPoisFilter == 1)
	{
		sprintf(strOutFileName, "%s%s_%s.p1", strOutputPath, strOutputFile, strChr);
		TileMapv2_SaveToBinaryFile((void *)vP1, sizeof(float), nBinNum, strOutFileName);

		if(nCTNum > 0)
		{
			sprintf(strOutFileName, "%s%s_%s.p0", strOutputPath, strOutputFile, strChr);
			TileMapv2_SaveToBinaryFile((void *)vP0, sizeof(float), nBinNum, strOutFileName);
		}
	}
	
	/* release memory */
	for(ni=0; ni<nFileNum; ni++)
	{
		free(vD[ni]);
		vD[ni] = NULL;

		free(vW[ni]);
		vW[ni] = NULL;
	}
	free(vD);
	free(vW);
	free(vM0);
	free(vM1);
	free(vV);
	free(vFC);

	if(nPoisFilter == 1)
	{
		free(vT0);
		free(vT1);
		free(vP0);
		free(vP1);
	}
	
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqPeak_ComputeBinStats_Chr_Tstat()                                    */
/*  Compute t-statistics.                                                  */
/* ----------------------------------------------------------------------- */ 
int SeqPeak_ComputeBinStats_Chr_Tstat(char strOutputPath[], char strOutputFile[], 
			char strChr[], int nBinNum, double dB, double dVM, int nIPNum, int nCTNum, 
			double *pTM, double *pTS)
{
	/* define */
	/* sample variance */
	float *vV = NULL;
	/* log2 fc */
	float *vM = NULL;
	/* t-stat */
	float *vT = NULL;
	/* file */
	FILE *fpIn;
	/* use var */
	int nUseVar = 0;

	char strOutFileName[MED_LINE_LENGTH];
	char strInFileName[MED_LINE_LENGTH];
	double dScale;
	double dTemp;

	int ni;
	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));

	/* init */
	dScale = 0.0;
	if(nIPNum > 0)
		dScale += 1.0/(float)nIPNum;
	if(nCTNum > 0)
		dScale += 1.0/(float)nCTNum;
	if( (nIPNum > 1) || (nCTNum > 1) )
		nUseVar = 1;

	vM = (float *)calloc(nBinNum, sizeof(float));
	if(vM == NULL)
	{
		printf("Error: SeqPeak_ComputeBinStats_Chr_Tstat, cannot allocate memory for M!\n");
		exit(EXIT_FAILURE);
	}

	vV = (float *)calloc(nBinNum, sizeof(float));
	if(vV == NULL)
	{
		printf("Error: SeqPeak_ComputeBinStats_Chr_Tstat, cannot allocate memory for V!\n");
		exit(EXIT_FAILURE);
	}

	vT = (float *)calloc(nBinNum, sizeof(float));
	if(vT == NULL)
	{
		printf("Error: SeqPeak_ComputeBinStats_Chr_Tstat, cannot allocate memory for T!\n");
		exit(EXIT_FAILURE);
	}

	/* load D */
	sprintf(strInFileName, "%s%s_%s.log2fc", strOutputPath, strOutputFile, strChr);
	fpIn = NULL;
	fpIn = fopen(strInFileName, "rb");
	if(fpIn == NULL)
	{
		printf("Error: SeqPeak_ComputeBinStats_Chr_Tstat, cannot open file %s!\n", strInFileName);
		exit(EXIT_FAILURE);
	}
	
	if( little_endian_fread(vM, sizeof(float), nBinNum, fpIn, little_endian_machine) != nBinNum)
	{
		printf("Error: SeqPeak_ComputeBinStats_Chr_Tstat, incorrect loading, number of bins inconsistent!\n");
		exit(EXIT_FAILURE);
	}

	fclose(fpIn);

	/* load V */
	sprintf(strInFileName, "%s%s_%s.v", strOutputPath, strOutputFile, strChr);
	fpIn = NULL;
	fpIn = fopen(strInFileName, "rb");
	if(fpIn == NULL)
	{
		printf("Error: SeqPeak_ComputeBinStats_Chr_Tstat, cannot open file %s!\n", strInFileName);
		exit(EXIT_FAILURE);
	}
	
	if( little_endian_fread(vV, sizeof(float), nBinNum, fpIn, little_endian_machine) != nBinNum)
	{
		printf("Error: SeqPeak_ComputeBinStats_Chr_Tstat, incorrect loading, number of bins inconsistent!\n");
		exit(EXIT_FAILURE);
	}

	fclose(fpIn);

	/* compute t-stat */
	for(ni=0; ni<nBinNum; ni++)
	{
		if(nUseVar == 1)
		{
			dTemp = (1.0-dB)*vV[ni]+dB*dVM;
			dTemp = sqrt(dScale*dTemp)+1e-6;
			vT[ni] = (float)(vM[ni]/dTemp);
		}
		else
		{
			vT[ni] = vM[ni];
		}

		*pTM = (*pTM) + vT[ni];
		*pTS = (*pTS) + vT[ni]*vT[ni];
	}

	/* save t-stat */
	sprintf(strOutFileName, "%s%s_%s.t", strOutputPath, strOutputFile, strChr);
	TileMapv2_SaveToBinaryFile((void *)vT, sizeof(float), nBinNum, strOutFileName);

	/* release memory */
	free(vM);
	free(vV);
	free(vT);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqPeak_ComputeBinStats_Chr_Tstandard()                                */
/*  Standardize t-statistics.                                              */
/* ----------------------------------------------------------------------- */
int SeqPeak_ComputeBinStats_Chr_Tstandard(char strOutputPath[], char strOutputFile[],
			char strChr[], int nBinNum, double dTM, double dTS)
{
	/* define */
	/* t-stat */
	float *vT = NULL;
	/* file */
	FILE *fpIn;

	char strOutFileName[MED_LINE_LENGTH];
	char strInFileName[MED_LINE_LENGTH];

	int ni;
	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));

	/* init */
	vT = (float *)calloc(nBinNum, sizeof(float));
	if(vT == NULL)
	{
		printf("Error: SeqPeak_ComputeBinStats_Chr_Tstandard, cannot allocate memory for T!\n");
		exit(EXIT_FAILURE);
	}

	/* load T */
	sprintf(strInFileName, "%s%s_%s.t", strOutputPath, strOutputFile, strChr);
	fpIn = NULL;
	fpIn = fopen(strInFileName, "rb");
	if(fpIn == NULL)
	{
		printf("Error: SeqPeak_ComputeBinStats_Chr_Tstandard, cannot open file %s!\n", strInFileName);
		exit(EXIT_FAILURE);
	}
	
	if( little_endian_fread(vT, sizeof(float), nBinNum, fpIn, little_endian_machine) != nBinNum)
	{
		printf("Error: SeqPeak_ComputeBinStats_Chr_Tstandard, incorrect loading, number of bins inconsistent!\n");
		exit(EXIT_FAILURE);
	}

	fclose(fpIn);

	/* compute t-stat */
	for(ni=0; ni<nBinNum; ni++)
	{
		vT[ni] = (float)((vT[ni]-dTM)/dTS);
	}

	/* save t-stat */
	sprintf(strOutFileName, "%s%s_%s.t", strOutputPath, strOutputFile, strChr);
	TileMapv2_SaveToBinaryFile((void *)vT, sizeof(float), nBinNum, strOutFileName);

	/* release memory */
	free(vT);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqPeak_CleanFiles()                                                   */
/*  Export results to BAR files and clean intermediate results.            */
/* ----------------------------------------------------------------------- */
int SeqPeak_CleanFiles(char strOutputPath[], char strOutputFile[],
		int nFileNum, struct tagString **vFileName,
		int nChrNum, struct tagString **vChrName, struct INTMATRIX * pChrLen, int nBinSize,
		int nExportBAR, int nKeepTempFiles, int nPoisFilter)
{
	float *vT = NULL;
	FILE *fpOut;
	FILE *fpIn;
	int ni,nj,nk;
	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));
	int nBinNum;
	int nPos;

	char strOutFileName[MED_LINE_LENGTH];
	char strFileName[LONG_LINE_LENGTH];
	char strCGWFileName[MED_LINE_LENGTH];


if(1) {
	/* save log2 fc file */
	printf("Export final statistics to BAR files ...\n");
	sprintf(strOutFileName, "%s%s_log2fc.bar.txt", strOutputPath, strOutputFile);

	fpOut = NULL;
	fpOut = fopen(strOutFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "#chr\tpos\t%s_log2fc\n", strOutputFile);
	fprintf(fpOut, "#chr\tpos\t1\n");

	for(ni=0; ni<nChrNum; ni++)
	{
		nBinNum = pChrLen->pMatElement[ni]/nBinSize;
		if(pChrLen->pMatElement[ni] % nBinSize != 0)
			nBinNum++;

		/* create memory */
		vT = NULL;
		vT = (float *)calloc(nBinNum, sizeof(float));
		if(vT == NULL)
		{
			printf("Error: SeqPeak_CleanFiles, cannot allocate memory for T!\n");
			exit(EXIT_FAILURE);
		}

		/* load data */
		sprintf(strFileName, "%s%s_%s.log2fc", strOutputPath, strOutputFile, vChrName[ni]->m_pString);
		fpIn = NULL;
		fpIn = fopen(strFileName, "rb");
		if(fpIn == NULL)
		{
			printf("Error: SeqPeak_CleanFiles, cannot open file %s!\n", strFileName);
			exit(EXIT_FAILURE);
		}
	
		if( little_endian_fread(vT, sizeof(float), nBinNum, fpIn, little_endian_machine) != nBinNum)
		{
			printf("Error: SeqPeak_CleanFiles, incorrect loading, number of bins inconsistent!\n");
			exit(EXIT_FAILURE);
		}

		fclose(fpIn);

		/* write to file */
		nPos = nBinSize/2;
		for(nj=0; nj<nBinNum; nj++)
		{
			if( (vT[nj] > 0.5) || (vT[nj] < -0.5) ) 
				fprintf(fpOut, "%s\t%d\t%f\n", vChrName[ni]->m_pString, nPos, vT[nj]);
			nPos += nBinSize;
		}

		/* free memory */
		free(vT);
	}
	fclose(fpOut);
	sprintf(strCGWFileName, "%s.cgw", strOutputFile);
	TileMapv2_TXT2BAR(strOutFileName, strOutputPath, strCGWFileName);
	RemoveFiles(strOutFileName);

	/* save t file */
	sprintf(strOutFileName, "%s%s_t.bar.txt", strOutputPath, strOutputFile);

	fpOut = NULL;
	fpOut = fopen(strOutFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "#chr\tpos\t%s_t\n", strOutputFile);
	fprintf(fpOut, "#chr\tpos\t1\n");

	for(ni=0; ni<nChrNum; ni++)
	{
		nBinNum = pChrLen->pMatElement[ni]/nBinSize;
		if(pChrLen->pMatElement[ni] % nBinSize != 0)
			nBinNum++;

		/* create memory */
		vT = NULL;
		vT = (float *)calloc(nBinNum, sizeof(float));
		if(vT == NULL)
		{
			printf("Error: SeqPeak_CleanFiles, cannot allocate memory for T!\n");
			exit(EXIT_FAILURE);
		}

		/* load data */
		sprintf(strFileName, "%s%s_%s.t", strOutputPath, strOutputFile, vChrName[ni]->m_pString);
		fpIn = NULL;
		fpIn = fopen(strFileName, "rb");
		if(fpIn == NULL)
		{
			printf("Error: SeqPeak_CleanFiles, cannot open file %s!\n", strFileName);
			exit(EXIT_FAILURE);
		}
	
		if( little_endian_fread(vT, sizeof(float), nBinNum, fpIn, little_endian_machine) != nBinNum)
		{
			printf("Error: SeqPeak_CleanFiles, incorrect loading, number of bins inconsistent!\n");
			exit(EXIT_FAILURE);
		}

		fclose(fpIn);

		/* write to file */
		nPos = nBinSize/2;
		for(nj=0; nj<nBinNum; nj++)
		{
			if( (vT[nj] > 1.0) || (vT[nj] < -1.0) ) 
				fprintf(fpOut, "%s\t%d\t%f\n", vChrName[ni]->m_pString, nPos, vT[nj]);
			nPos += nBinSize;
		}

		/* free memory */
		free(vT);
	}
	fclose(fpOut);
	sprintf(strCGWFileName, "%s.cgw", strOutputFile);
	TileMapv2_TXT2BAR(strOutFileName, strOutputPath, strCGWFileName);
	RemoveFiles(strOutFileName);
	RemoveFiles(strCGWFileName);
}

	/* remove redundant files */
	if(nKeepTempFiles == 0)
	{
		printf("Clean intermediate files ...\n");
		sprintf(strFileName, "%s%s_*.t", strOutputPath, strOutputFile);
		RemoveFiles(strFileName);
		sprintf(strFileName, "%s%s_*.v", strOutputPath, strOutputFile);
		RemoveFiles(strFileName);
		sprintf(strFileName, "%s%s_*.log2fc", strOutputPath, strOutputFile);
		RemoveFiles(strFileName);
		sprintf(strFileName, "%s%s_*.bincount", strOutputPath, strOutputFile);
		RemoveFiles(strFileName);
		sprintf(strFileName, "%s%s.pos.reg", strOutputPath, strOutputFile);
		RemoveFiles(strFileName);
		sprintf(strFileName, "%s%s.neg.reg", strOutputPath, strOutputFile);
		RemoveFiles(strFileName);

		if(nPoisFilter == 1)
		{
			sprintf(strFileName, "%s%s_*.p1", strOutputPath, strOutputFile);
			RemoveFiles(strFileName);
			sprintf(strFileName, "%s%s_*.p0", strOutputPath, strOutputFile);
			RemoveFiles(strFileName);
		}
	}

	/* save file bin stat */
	if(nExportBAR == 1)
	{
		for(nk=0; nk<nFileNum; nk++)
		{
			printf("Export bin statistics of %s to BAR file...\n", vFileName[nk]->m_pString);
			sprintf(strOutFileName, "%s%s_%s.bar.txt", strOutputPath, strOutputFile, vFileName[nk]->m_pString);
			sprintf(strCGWFileName, "%s.cgw", strOutputFile);
			TileMapv2_TXT2BAR(strOutFileName, strOutputPath, strCGWFileName);
			
			RemoveFiles(strOutFileName);
			sprintf(strFileName, "%s%s", strOutputPath, strCGWFileName);
			RemoveFiles(strCGWFileName);
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqPeak_Segmentation()                                                 */
/*  Genome segmentation. Find peaks.                                       */
/* ----------------------------------------------------------------------- */ 
int SeqPeak_Segmentation(char strOutputPath[], char strOutputFile[], 
				int nChrNum, struct tagString **vChrName, struct INTMATRIX *pChrLen, 
				int nBinSize, double dCutoff, int nMaxGap, int nMinLen, 
				int nFlipGroupLabel, int nPoisFilter, 
				int nPoisFilterLabel, double dPoisCut)
{
	/* define */
	int ni,nj,nk;
	int nBinNum;
	float *vT;
	float *vF;
	float *vP;
	char strLine[MED_LINE_LENGTH];
	char strFileName[MED_LINE_LENGTH];
	FILE *fpOut;
	FILE *fpReg;
	char strOutFileName[MED_LINE_LENGTH];
	int nP1,nP2,nTMax,nFMax,nPvMin;
	float tMax;
	float fcMax;
	float pvMin;
	int nRegionNum = 0;
	int nNrRegionNum = 0;
	int nLastP1,nLastP2,nLastTMax,nLastFMax,nLastPvMin;
	float dLastT,dLastF,dLastPv;
	int nPassCut;
		
	/* open the region file */
	if(nFlipGroupLabel == 0)
		sprintf(strOutFileName, "%s%s.pos.reg", strOutputPath, strOutputFile);
	else
		sprintf(strOutFileName, "%s%s.neg.reg", strOutputPath, strOutputFile);

	fpReg = NULL;
	fpReg = fopen(strOutFileName, "w");
	if(fpReg == NULL)
	{
		printf("Error: SeqPeak_Segmentation, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	
	/* process chromosome by chromosome */
	for(ni=0; ni<nChrNum; ni++)
	{
		/* load data */
		nBinNum = pChrLen->pMatElement[ni]/nBinSize;
		if(pChrLen->pMatElement[ni] % nBinSize != 0)
			nBinNum++;

		vT = (float *)calloc(nBinNum, sizeof(float));
		if(vT == NULL)
		{
			printf("Error: SeqPeak_Segmentation, cannot allocate memory for bin stat!\n");
			exit(EXIT_FAILURE);
		}

		vF = (float *)calloc(nBinNum, sizeof(float));
		if(vF == NULL)
		{
			printf("Error: SeqPeak_Segmentation, cannot allocate memory for bin stat!\n");
			exit(EXIT_FAILURE);
		}

		sprintf(strFileName, "%s%s_%s.t", strOutputPath, strOutputFile, vChrName[ni]->m_pString);
		if( TileMapv2_LoadFromBinaryFile(vT, sizeof(float), nBinNum, strFileName) != PROC_SUCCESS)
		{
			printf("Error: SeqPeak_Segmentation, cannot read bin statistics correctly!\n");
			exit(EXIT_FAILURE);
		}

		sprintf(strFileName, "%s%s_%s.log2fc", strOutputPath, strOutputFile, vChrName[ni]->m_pString);
		if( TileMapv2_LoadFromBinaryFile(vF, sizeof(float), nBinNum, strFileName) != PROC_SUCCESS)
		{
			printf("Error: SeqPeak_Segmentation, cannot read bin statistics correctly!\n");
			exit(EXIT_FAILURE);
		}

		
		if(nPoisFilter == 1)
		{
			vP = (float *)calloc(nBinNum, sizeof(float));
			if(vP == NULL)
			{
				printf("Error: SeqPeak_Segmentation, cannot allocate memory for bin stat!\n");
				exit(EXIT_FAILURE);
			}

			if(nPoisFilterLabel == 1)		
				sprintf(strFileName, "%s%s_%s.p1", strOutputPath, strOutputFile, vChrName[ni]->m_pString);
			else
				sprintf(strFileName, "%s%s_%s.p0", strOutputPath, strOutputFile, vChrName[ni]->m_pString);

			if( TileMapv2_LoadFromBinaryFile(vP, sizeof(float), nBinNum, strFileName) != PROC_SUCCESS)
			{
				printf("Error: SeqPeak_Segmentation, cannot read bin statistics correctly!\n");
				exit(EXIT_FAILURE);
			}
		}

		/* open output file */
		sprintf(strOutFileName, "%s%s_%s.tmpreg", strOutputPath, strOutputFile, vChrName[ni]->m_pString);
		fpOut = NULL;
		fpOut = fopen(strOutFileName, "w");
		if(fpOut == NULL)
		{
			printf("Error: SeqPeak_Segmentation, cannot open file %s to export data!\n", strOutFileName);
			exit(EXIT_FAILURE);
		}

		/* find peaks */
		nP1 = -1;
		nP2 = -1;
		nTMax = -1;
		tMax = -1000000.0;
		nFMax = -1;
		fcMax = -1000000.0;
		pvMin = -1.0;
		nPvMin = -1;
		nRegionNum = 0;
		for(nj=0; nj<nBinNum; nj++)
		{
			if(nFlipGroupLabel == 1)
			{
				vT[nj] = -vT[nj];
				vF[nj] = -vF[nj];
			}

			nPassCut = 0;
			if(vT[nj] >= dCutoff)
			{
				nPassCut = 1;
			}

			if(nPassCut == 1)
			{
				if(nP1 < 0)
				{
					nP1 = nj;
					nP2 = nj;
					nTMax = nj;
					nFMax = nj;
					tMax = vT[nj];
					fcMax = vF[nj];

					if(nPoisFilter == 1)
					{
						nPvMin = nj;
						pvMin = vP[nj];						
					}
				}
				else
				{
					nP2 = nj;
					if(vT[nj] > tMax)
					{
						tMax = vT[nj];
						nTMax = nj;
					}
					if(vF[nj] > fcMax)
					{
						fcMax = vF[nj];
						nFMax = nj;
					}
					if(nPoisFilter == 1) 
					{
						if(vP[nj] > pvMin)
						{
							pvMin = vP[nj];
							nPvMin = nj;
						}
					}
				}
			}
			else
			{
				if(nP1 >= 0)
				{
					fprintf(fpOut, "%d\t%d\t%f\t%d\t%f\t%d\t%f\t%d\n", nP1, nP2, tMax, nTMax, fcMax, nFMax, pvMin, nPvMin);
					nRegionNum++;
					nP1 = -1;
					nP2 = -1;
					nTMax = -1;
					nFMax = -1;
					nPvMin = -1;
					tMax = -1000000.0;
					fcMax = -1000000.0;
					pvMin = -1.0;
				}
			}
		}

		if(nP1 >= 0)
		{
			fprintf(fpOut, "%d\t%d\t%f\t%d\t%f\t%d\t%f\t%d\n", nP1, nP2, tMax, nTMax, fcMax, nFMax, pvMin, nPvMin);
			nRegionNum++;
			nP1 = -1;
			nP2 = -1;
			nTMax = -1;
			nFMax = -1;
			nPvMin = -1;
			tMax = -1000000.0;
			fcMax = -1000000.0;
			pvMin = -1.0;
		}

		/* free memory */
		free(vT);
		free(vF);
		if(nPoisFilter == 1)
			free(vP);

		/* close file */
		fclose(fpOut);

		if(nRegionNum > 0)
		{
			/* reload regions */
			fpOut = NULL;
			fpOut = fopen(strOutFileName, "r");
			if(fpOut == NULL)
			{
				printf("Error: SeqPeak_Segmentation, cannot open file %s to load data!\n", strOutFileName);
				exit(EXIT_FAILURE);
			}

			fgets(strLine, MED_LINE_LENGTH, fpOut);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			sscanf(strLine, "%d %d %f %d %f %d %f %d", &nLastP1, &nLastP2, &dLastT, &nLastTMax, &dLastF, &nLastFMax, &dLastPv, &nLastPvMin);

			nk = 1;
			while(fgets(strLine, MED_LINE_LENGTH, fpOut) != NULL)
			{
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] == '\0')
					continue;

				if(nk >= nRegionNum)
				{
					printf("Error: SeqPeak_Segmentation, inconsistent region number!\n");
					exit(EXIT_FAILURE);
				}

				sscanf(strLine, "%d %d %f %d %f %d %f %d", &nP1, &nP2, &tMax, &nTMax, &fcMax, &nFMax, &pvMin, &nPvMin);

				if( (nP1-nLastP2)*nBinSize <= nMaxGap )
				{
					nLastP2 = nP2;
					if(tMax > dLastT)
					{
						dLastT = tMax;
						nLastTMax = nTMax;
					}
					if(fcMax > dLastF)
					{
						dLastF = fcMax;
						nLastFMax = nFMax;
					}
					if(nPoisFilter == 1)
					{
						if(pvMin > dLastPv)
						{
							dLastPv = pvMin;
							nLastPvMin = nPvMin;
						}
					}
				}
				else 
				{
					/* save old region */
					if( (nLastP2-nLastP1+1)*nBinSize >= nMinLen )
					{
						fprintf(fpReg, "%d\t%d\t%d\t%f\t%d\t%f\t%d\t%f\t%d\n", ni, nLastP1, nLastP2, dLastT, nLastTMax, dLastF, nLastFMax, dLastPv, nLastPvMin);
						nNrRegionNum++;
					}

					/* create new region */
					nLastP1 = nP1;
					nLastP2 = nP2;
					nLastTMax = nTMax;
					nLastFMax = nFMax;
					nLastPvMin = nPvMin;
					dLastT = tMax;
					dLastF = fcMax;
					dLastPv = pvMin;
				}

				nk++;
			}

			if( (nLastP2-nLastP1+1)*nBinSize >= nMinLen )
			{
				fprintf(fpReg, "%d\t%d\t%d\t%f\t%d\t%f\t%d\t%f\t%d\n", ni, nLastP1, nLastP2, dLastT, nLastTMax, dLastF, nLastFMax, dLastPv, nLastPvMin);
				nNrRegionNum++;
			}

			fclose(fpOut);

			if(nk != nRegionNum)
			{
				printf("Error: SeqPeak_Segmentation, inconsistent region number!\n");
				exit(EXIT_FAILURE);
			}
		}
	}

	/* close file */
	fclose(fpReg);

	/* remove temporary files */
	sprintf(strOutFileName, "%s%s_*.tmpreg", strOutputPath, strOutputFile);
	RemoveFiles(strOutFileName);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqPeak_FDR_Flip()                                                     */
/*  Compute FDR.                                                           */
/* ----------------------------------------------------------------------- */ 
int SeqPeak_FDR_Flip(char strOutputPath[], char strOutputFile[], 
			struct DOUBLEMATRIX **ppRegion,
			struct DOUBLEMATRIX **ppRegionSort, struct LONGMATRIX **ppRegionSid,
			struct DOUBLEMATRIX **ppRegionFDR)
{
	/* define */
	char strFileName[MED_LINE_LENGTH];
	struct DOUBLEMATRIX *pRegion = NULL;
	struct DOUBLEMATRIX *pRegionCT = NULL;
	struct INTMATRIX *pPriority = NULL;
	struct INTMATRIX *pType = NULL;
	struct DOUBLEMATRIX *pRegionMaxSort = NULL;
	struct LONGMATRIX *pRegionMaxSid = NULL;
	struct DOUBLEMATRIX *pRegionCTMaxSort = NULL;
	struct LONGMATRIX *pRegionCTMaxSid = NULL;
	struct DOUBLEMATRIX *pRegionMaxFDR = NULL;
	int nCol;

	/* load regions */
	sprintf(strFileName, "%s%s.pos.reg", strOutputPath, strOutputFile);
	pRegion = DMLOAD(strFileName);
	if(pRegion == NULL)
	{
		return PROC_SUCCESS;
	}

	/* data type */
	pType = CreateIntMatrix(1, pRegion->nWidth);
	if(pType == NULL)
	{
		printf("Error: SeqPeak_FDR_Flip, cannot create column type!\n");
		exit(EXIT_FAILURE);
	}
	pType->pMatElement[0] = 2;
	pType->pMatElement[1] = 2;
	pType->pMatElement[2] = 2;
	pType->pMatElement[3] = 1;
	pType->pMatElement[4] = 2;
	pType->pMatElement[5] = 1;
	pType->pMatElement[6] = 2;
	pType->pMatElement[7] = 1;
	pType->pMatElement[8] = 2;

	/* sorting index */
	pPriority = NULL;
	pPriority = CreateIntMatrix(1, 1);
	if(pPriority == NULL)
	{
		printf("Error: SeqPeak_FDR_Flip, cannot create sorting key!\n");
		exit(EXIT_FAILURE);
	}
	pPriority->pMatElement[0] = 3;

	/* sort peak regions */
	DMSORTROWS(pRegion, pType, pPriority, &pRegionMaxSort, &pRegionMaxSid);

	pRegionMaxFDR = CreateDoubleMatrix(1, pRegion->nHeight);
	if(pRegionMaxFDR == NULL)
	{
		printf("Error: SeqPeak_FDR_Flip, cannot create matrix for region FDR!\n");
		exit(EXIT_FAILURE);
	}

	/* load control region */
	sprintf(strFileName, "%s%s.neg.reg", strOutputPath, strOutputFile);
	pRegionCT = DMLOAD(strFileName);
	
	if(pRegionCT != NULL)
	{
		/* sort region, regionCT */
		pRegionCTMaxSort = NULL;
		pRegionCTMaxSid = NULL;
		DMSORTROWS(pRegionCT, pType, pPriority, &pRegionCTMaxSort, &pRegionCTMaxSid);

		/* count */
		nCol = 3;
		TileMapv2_RegionFDR_Count(pRegionMaxSort, pRegionCTMaxSort, pRegionMaxFDR, nCol);

		TileMapv2_RegionFDR_Compute(pRegionMaxFDR);

		/* -------------- */
		/* release memory */
		/* -------------- */
		DestroyDoubleMatrix(pRegionCTMaxSort);
		DestroyLongMatrix(pRegionCTMaxSid);
		DestroyDoubleMatrix(pRegionCT);
	}
	else
	{
		/* FDR == 0*/
	}

	/* release memory */
	DestroyIntMatrix(pPriority);
	DestroyIntMatrix(pType);

	*ppRegion = pRegion;
	pRegion = NULL;
	*ppRegionSort = pRegionMaxSort;
	pRegionMaxSort = NULL;
	*ppRegionSid = pRegionMaxSid;
	pRegionMaxSid = NULL;
	*ppRegionFDR = pRegionMaxFDR;
	pRegionMaxFDR = NULL;

	/* return */
	return PROC_SUCCESS;
}

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
		int nPoisFilter, double dPoisCut)
{
	/* define */
	FILE *fpOut = NULL;
	char strFileName[MED_LINE_LENGTH];
	int ni,nRi,nRank;
	long nj,nk;
	int nChrId;
	int nP1,nP2,nTMax,nFMax,nPvMin;
	float tMax,fcMax,pvMin;
	int nLeftB,nRightB,nModeB,nMidB;
	double dLog10PoisCut = -log(dPoisCut)/log(10.0);

	/* open file */
	sprintf(strFileName, "%s%s_peak.cod", strOutputPath, strOutputFile);
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: SeqPeak_ExportPeaks, cannot open the output file to export peaks!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "#rank\tchromosome\tstart\tend\tstrand\tpeak_length\tFDR");
	if(nBoundaryRefine == 1)
		fprintf(fpOut, "\tleft_peakboundary\tright_peakboundary\tpeak_summit\tbound_center\tbound_width");
	fprintf(fpOut, "\tmaxT\tmaxT_pos\tmax_log2FC\tmaxFC_pos\tminuslog10_minPoisP\tminPoisP_pos");
	if(nCollectRawData == 1)
	{
		for(ni=0; ni<nFileNum; ni++)
			fprintf(fpOut, "\t%s", vFileName[ni]->m_pString);
	}
	fprintf(fpOut, "\n");

	/* write results */
	if(pRegionSort != NULL)
	{
		nRank = 0;
		for(ni=0; ni<pRegionSort->nHeight; ni++)
		{
			nRi = pRegionSort->nHeight-ni-1;
			nChrId = (int)(DMGETAT(pRegionSort, nRi, 0));
			if( (nChrId < 0) || (nChrId >= nChrNum) )
			{
				printf("Error: SeqPeak_ExportPeaks, chromosome index out of range!\n");
				exit(EXIT_FAILURE);
			}

			nP1 = (int)( nBinSize * DMGETAT(pRegionSort, nRi, 1) );
			nP2 = (int)( nBinSize * DMGETAT(pRegionSort, nRi, 2) ) + nBinSize - 1;
			tMax = (float)(DMGETAT(pRegionSort, nRi, 3));
			nTMax = (int)( nBinSize * DMGETAT(pRegionSort, nRi, 4) ) + nBinSize/2;
			fcMax = (float)(DMGETAT(pRegionSort, nRi, 5));
			nFMax = (int)( nBinSize * DMGETAT(pRegionSort, nRi, 6) ) + nBinSize/2;
			pvMin = (float)(DMGETAT(pRegionSort, nRi, 7));
			nPvMin = (int)( nBinSize * DMGETAT(pRegionSort, nRi, 8) ) + nBinSize/2;

			if(nPoisFilter == 1)
			{
				if(pvMin < dLog10PoisCut)
					continue;
			}

			nRank++;

			fprintf(fpOut, "%d\t%s\t%d\t%d\t+\t%d\t%f", nRank, vChrName[nChrId]->m_pString, nP1, nP2, nP2-nP1+1,
				pRegionFDR->pMatElement[nRi]);

			nj = pRegionSid->pMatElement[nRi];
			if(nBoundaryRefine == 1)
			{
				nLeftB = IMGETAT(pPeakBoundary, nj, 0);
				nRightB = IMGETAT(pPeakBoundary, nj, 1);
				nModeB = IMGETAT(pPeakBoundary, nj, 2);
				nMidB = (nLeftB+nRightB)/2;
				fprintf(fpOut, "\t%d\t%d\t%d\t%d\t%d", nLeftB, nRightB, nModeB, nMidB, nRightB-nLeftB+1);
			}

			fprintf(fpOut, "\t%f\t%d\t%f\t%d\t%f\t%d", tMax, nTMax, fcMax, nFMax, pvMin, nPvMin);

			if(nCollectRawData == 1)
			{
				for(nk=0; nk<nFileNum; nk++)
					fprintf(fpOut, "\t%f", (float)(DMGETAT(pPeakRawData, nj, nk)));
			}

			fprintf(fpOut, "\n");
		}
	}

	/* close file */
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}

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
			struct INTMATRIX **ppPeakBoundary)
{
	/* define */
	int nRegionNum;
	int *vChrId = NULL;
	int *vStart = NULL;
	int *vEnd = NULL;
	int *vBRWinNum = NULL;
	float **vPosCount = NULL;
	float **vNegCount = NULL;
	float **vTotCount = NULL;
	int ni,nRegLen;


	/* STEP 0: initialize */
	if(pRegion == NULL)
		return PROC_SUCCESS;
	nRegionNum = pRegion->nHeight;
	if(nRegionNum == 0)
		return PROC_SUCCESS;

	/* STEP 1: Prepare Regions & Space */
	vChrId = (int *)calloc(nRegionNum, sizeof(int));
	if(vChrId == NULL)
	{
		printf("Error: SeqPeak_RefineBoundary, cannot create memory for vChrId!\n");
		exit(EXIT_FAILURE);
	}

	vStart = (int *)calloc(nRegionNum, sizeof(int));
	if(vStart == NULL)
	{
		printf("Error: SeqPeak_RefineBoundary, cannot create memory for vStart!\n");
		exit(EXIT_FAILURE);
	}

	vEnd = (int *)calloc(nRegionNum, sizeof(int));
	if(vEnd == NULL)
	{
		printf("Error: SeqPeak_RefineBoundary, cannot create memory for vEnd!\n");
		exit(EXIT_FAILURE);
	}

	vBRWinNum = (int *)calloc(nRegionNum, sizeof(int));
	if(vBRWinNum == NULL)
	{
		printf("Error: SeqPeak_RefineBoundary, cannot create memory for vBRWinNum!\n");
		exit(EXIT_FAILURE);
	}

	vPosCount = (float **)calloc(nRegionNum, sizeof(float *));
	if(vPosCount == NULL)
	{
		printf("Error: SeqPeak_RefineBoundary, cannot create memory for vPosCount!\n");
		exit(EXIT_FAILURE);
	}

	vNegCount = (float **)calloc(nRegionNum, sizeof(float *));
	if(vNegCount == NULL)
	{
		printf("Error: SeqPeak_RefineBoundary, cannot create memory for vNegCount!\n");
		exit(EXIT_FAILURE);
	}

	vTotCount = (float **)calloc(nRegionNum, sizeof(float *));
	if(vTotCount == NULL)
	{
		printf("Error: SeqPeak_RefineBoundary, cannot create memory for vNegCount!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nRegionNum; ni++)
	{
		vChrId[ni] = (int)(DMGETAT(pRegion, ni, 0));
		vStart[ni] = (int)(nBinSize * DMGETAT(pRegion, ni, 1));
		vEnd[ni] = (int)(nBinSize * DMGETAT(pRegion, ni, 2)) + nBinSize - 1;
		nRegLen = vEnd[ni]-vStart[ni]+1;
		vBRWinNum[ni] = nRegLen/nBRWin;
		if( nRegLen%nBRWin != 0 )
			vBRWinNum[ni] += 1;

		vPosCount[ni] = (float *)calloc(vBRWinNum[ni], sizeof(float));
		if(vPosCount[ni] == NULL)
		{
			printf("Error: SeqPeak_RefineBoundary, cannot create memory for vPosCount elements!\n");
			exit(EXIT_FAILURE);
		}

		vNegCount[ni] = (float *)calloc(vBRWinNum[ni], sizeof(float));
		if(vNegCount[ni] == NULL)
		{
			printf("Error: SeqPeak_RefineBoundary, cannot create memory for vNegCount elements!\n");
			exit(EXIT_FAILURE);
		}

		vTotCount[ni] = (float *)calloc(vBRWinNum[ni], sizeof(float));
		if(vTotCount[ni] == NULL)
		{
			printf("Error: SeqPeak_RefineBoundary, cannot create memory for vTotCount elements!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* STEP 2: Bring in reads */
	SeqPeak_RefineBoundary_RegionProfile(nRegionNum, 
		vChrId, vStart, vEnd, vBRWinNum, 
		vPosCount, vNegCount, vTotCount,
		nBRWin,	strInputPath, nIPNum, nCTNum, 
		nFileNum, vGroupId, pFileNormFactor, 
		nChrNum, vChrName, pChrLen, nExtLen);

	/* STEP 3: Find peak and boundaries for each region */
	*ppPeakBoundary = NULL;
	*ppPeakBoundary = CreateIntMatrix(nRegionNum, 3);
	if(*ppPeakBoundary == NULL)
	{
		printf("Error: SeqPeak_RefineBoundary, cannot create memory for boundary coordinates!\n");
		exit(EXIT_FAILURE);
	}
	SeqPeak_RefineBoundary_DetectBoundary(nRegionNum, 
		vChrId, vStart, vEnd, vBRWinNum, 
		vPosCount, vNegCount, vTotCount,
		nBRWin,	nExtLen, *ppPeakBoundary);


	/* release memory */
	free(vChrId);
	free(vStart);
	free(vEnd);
	free(vBRWinNum);
	for(ni=0; ni<nRegionNum; ni++)
	{
		free(vPosCount[ni]);
		vPosCount[ni] = NULL;
		free(vNegCount[ni]);
		vNegCount[ni] = NULL;
		free(vTotCount[ni]);
		vTotCount[ni] = NULL;
	}
	free(vPosCount);
	free(vNegCount);
	free(vTotCount);

	/* return */
	return PROC_SUCCESS;
}

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
		struct INTMATRIX *pChrLen, int nExtLen)
{
	/* define */
	FILE *fpIn;
	FILE *fpRead;
	char strLine[LONG_LINE_LENGTH];
	char *chp;
	char strChr[LINE_LENGTH];
	int nPos,nStrand;
	int nId;
	int nO1,nO2;
	int nP1,nP2;
	int nQ1,nQ2;
	int nx,ny;
	int ni,nIPi,nj,nk;

	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: SeqPeak_RefineBoundary_RegionProfile, cannot open input file list!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	nIPi = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		/* read alignment file */
		chp = strchr(strLine, '\t');
		*chp = '\0';
		
		if(vGroupId[ni] == 0)
		{
			ni++;
			continue;
		}

		nIPi++;
		printf("  Processing %s ...\n", strLine);
		
		fpRead = NULL;
		fpRead = fopen(strLine, "r");
		if(fpRead == NULL)
		{
			printf("Error: SeqPeak_RefineBoundary_RegionProfile, cannot open read alignment file!\n");
			exit(EXIT_FAILURE);
		}
		while(fgets(strLine, LONG_LINE_LENGTH, fpRead) != NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;
			if(strLine[0] == '#')
				continue;

			SeqClust_ReadPos_From_Aln(strLine, strChr, &nPos, &nStrand);

			/* find chromosome and update chromosome length */
			nId = SeqClust_FindChr(vChrName, nChrNum, strChr);
			if( (nId < 0) || (nId >= nChrNum) )
			{
				printf("Warning: SeqPeak_RefineBoundary_RegionProfile, cannot find the matching chromosome for %s:%d!\n",strChr, nPos);
				continue;
			}

			/* add pos and neg read to region */
			nx = 0;
			ny = -1;
			nO1 = nPos-nExtLen/2;
			nO2 = nPos+nExtLen/2;
			if(nStrand == 0)
			{
				nP1 = nPos;
				nP2 = nPos+nExtLen;
				/* find hit id nx, ny */
				SeqPeak_RefineBoundary_FindOverlapRegions(nId, nO1, nP2, 
					nRegionNum, vChrId, vStart, vEnd, &nx, &ny);
			}
			else
			{
				nP1 = nPos-nExtLen;
				nP2 = nPos;
				/* find hit id nx, ny */
				SeqPeak_RefineBoundary_FindOverlapRegions(nId, nP1, nO2, 
					nRegionNum, vChrId, vStart, vEnd, &nx, &ny);
			}
			
			/* update profile */
			for(nk=nx; nk<=ny; nk++)
			{
				/* update total profile */
				nQ1 = nP1;
				nQ2 = nP2;
				if(nQ1 < vStart[nk])
					nQ1 = vStart[nk];
				if(nQ2 > vEnd[nk])
					nQ2 = vEnd[nk];
				
				if(nQ1 <= nQ2)
				{
					nQ1 = (nQ1-vStart[nk])/nBRWin;
					nQ2 = (nQ2-vStart[nk])/nBRWin;
					
					if( (nQ1 > vBRWinNum[nk]) || (nQ2 > vBRWinNum[nk]) )
					{
						printf("Error: SeqPeak_RefineBoundary_RegionProfile, window index out of range!\n");
						exit(EXIT_FAILURE);
					}

					for(nj=nQ1; nj<=nQ2; nj++)
						vTotCount[nk][nj] += (float)(pFileNormFactor->pMatElement[ni]);
				}

				/* update strand profile */
				nQ1 = nO1;
				nQ2 = nO2;
				if(nQ1 < vStart[nk])
					nQ1 = vStart[nk];
				if(nQ2 > vEnd[nk])
					nQ2 = vEnd[nk];
				
				if(nQ1 <= nQ2)
				{
					nQ1 = (nQ1-vStart[nk])/nBRWin;
					nQ2 = (nQ2-vStart[nk])/nBRWin;
					
					if( (nQ1 > vBRWinNum[nk]) || (nQ2 > vBRWinNum[nk]) )
					{
						printf("Error: SeqPeak_RefineBoundary_RegionProfile, window index out of range!\n");
						exit(EXIT_FAILURE);
					}

					if(nStrand == 0)
					{
						for(nj=nQ1; nj<=nQ2; nj++)
							vPosCount[nk][nj] += (float)(pFileNormFactor->pMatElement[ni]);
					}
					else
					{
						for(nj=nQ1; nj<=nQ2; nj++)
							vNegCount[nk][nj] += (float)(pFileNormFactor->pMatElement[ni]);
					}
				}
			}
		}	

		/* close alignment file */
		fclose(fpRead);
		ni++;
	}
	
	fclose(fpIn);
	
	if(ni != nFileNum)
	{
		printf("Error: SeqPeak_RefineBoundary, inconsistent file number!\n");
		exit(EXIT_FAILURE);
	}

	if(nIPi != nIPNum)
	{
		printf("Error: SeqPeak_RefineBoundary, inconsistent IP sample number!\n");
		exit(EXIT_FAILURE);
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqPeak_RefineBoundary_FindOverlapRegions()                            */
/*  Find regions overlapping with an input region. Return region index.    */
/* ----------------------------------------------------------------------- */ 
int SeqPeak_RefineBoundary_FindOverlapRegions(int nChrId, int nP1, int nP2, 
		int nRegionNum, int *vChrId, int *vStart, int *vEnd, int *px, int *py)
{
	/* define */
	int ni,nj,nk;
	int nResult;
	int nHit;

	/* init */
	ni = 0;
	nj = nRegionNum-1;
	nk = -1;
	nHit = 0;


	/* compare to first region */
	nResult = SeqPeak_RegionCompare(nChrId, nP1, nP2, vChrId[ni], vStart[ni], vEnd[ni]);
	if( nResult < 0 )
	{
		*px = 0;
		*py = -1;
		return PROC_SUCCESS;
	}
	else if( nResult == 0)
	{
		nHit = 1;
		nk = 0;
	}

	/* compare to last region */
	nResult = SeqPeak_RegionCompare(nChrId, nP1, nP2, vChrId[nj], vStart[nj], vEnd[nj]);
	if(  nResult > 0 )
	{
		*px = nRegionNum;
		*py = nRegionNum-1;
		return PROC_SUCCESS;
	}
	else if( nResult == 0)
	{
		nHit = 1;
		nk = nRegionNum-1;
	}

	/* compare middle regions */
	if(nk < 0)
	{
		while( (nj-ni) > 1)
		{
			nk = (ni+nj)/2;
			nResult = SeqPeak_RegionCompare(nChrId, nP1, nP2, vChrId[nk], vStart[nk], vEnd[nk]);

			if(nResult == 0)
			{
				nHit = 1;
				break;
			}
			else if(nResult < 0)
			{
				nj = nk;
			}
			else
			{
				ni = nk;
			}
		}

		if(nHit == 0)
		{
			nResult = SeqPeak_RegionCompare(nChrId, nP1, nP2, vChrId[ni], vStart[ni], vEnd[ni]);
			if(nResult == 0)
			{
				nHit = 1;
				nk = ni;
			}

			nResult = SeqPeak_RegionCompare(nChrId, nP1, nP2, vChrId[nj], vStart[nj], vEnd[nj]);
			if(nResult == 0)
			{
				nHit = 1;
				nk = nj;
			}
		}
	}

	/* extend matching */
	if(nHit == 0)
	{
		*px = 0;
		*py = -1;
		return PROC_SUCCESS;
	}

	for(ni=nk-1; ni>=0; ni--)
	{
		nResult = SeqPeak_RegionCompare(vChrId[ni], vStart[ni], vEnd[ni], nChrId, nP1, nP2);
		if(nResult < 0)
			break;
	}
	*px = ni+1;

	for(nj=nk+1; nj<nRegionNum; nj++)
	{
		nResult = SeqPeak_RegionCompare(nChrId, nP1, nP2, vChrId[nj], vStart[nj], vEnd[nj]);
		if(nResult < 0)
			break;
	}
	*py = nj-1;

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqPeak_RegionCompare()                                                */
/*  Compare two regions.                                                   */
/*  Return  -1 if region1 < region2;                                       */
/*           0 if region1 = region2;                                       */
/*           1 if region1 > region2;                                       */
/* ----------------------------------------------------------------------- */ 
int SeqPeak_RegionCompare(int nChr1, int nStart1, int nEnd1,
						  int nChr2, int nStart2, int nEnd2)
{
	/* define */
	int nResult;

	/* compare */
	if(nChr1 < nChr2)
		nResult = -1;
	else if(nChr1 > nChr2)
		nResult = 1;
	else
	{
		if(nEnd1 < nStart2)
			nResult = -1;
		else if(nStart1 > nEnd2)
			nResult = 1;
		else
			nResult = 0;
	}

	/* return */
	return nResult;
}

/* ----------------------------------------------------------------------- */ 
/*  SeqPeak_RefineBoundary_DetectBoundary()                                */
/*  Find boundaries.                                                       */
/* ----------------------------------------------------------------------- */ 
int SeqPeak_RefineBoundary_DetectBoundary(int nRegionNum, 
		int *vChrId, int *vStart, int *vEnd, int *vBRWinNum, 
		float **vPosCount, float **vNegCount, float **vTotCount,
		int nBRWin, int nExtLen, struct INTMATRIX *pPeakBoundary)
{
	/* define */
	int ni,nj;
	float fMax,rMax,tMax;
	int fPos,rPos,tPos;

	/* process region by region */
	for(ni=0; ni<nRegionNum; ni++)
	{
		fMax = vPosCount[ni][0];
		rMax = vNegCount[ni][0];
		tMax = vTotCount[ni][0];
		fPos = 0;
		rPos = 0;
		tPos = 0;

		for(nj=1; nj<vBRWinNum[ni]; nj++)
		{
			if(vPosCount[ni][nj] > fMax)
			{
				fMax = vPosCount[ni][nj];
				fPos = nj;
			}

			if(vNegCount[ni][nj] > rMax)
			{
				rMax = vNegCount[ni][nj];
				rPos = nj;
			}

			if(vTotCount[ni][nj] > tMax)
			{
				tMax = vTotCount[ni][nj];
				tPos = nj;
			}
		}

		fPos = vStart[ni]+fPos*nBRWin+nBRWin/2;
		rPos = vStart[ni]+rPos*nBRWin+nBRWin/2;
		tPos = vStart[ni]+tPos*nBRWin+nBRWin/2;
		IMSETAT(pPeakBoundary, ni, 0, fPos);
		IMSETAT(pPeakBoundary, ni, 1, rPos);
		IMSETAT(pPeakBoundary, ni, 2, tPos);
	}

	/* return */
	return PROC_SUCCESS;
}

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
		struct DOUBLEMATRIX **ppPeakRawData)
{
	/* define */
	int ni,nRegionNum;
	int *vChrId = NULL;
	int *vStart = NULL;
	int *vEnd = NULL;

	/* init */
	*ppPeakRawData = NULL;

	if(pRegion == NULL)
		return PROC_SUCCESS;
	if(pRegion->nHeight == 0)
		return PROC_SUCCESS;
	if(nFileNum == 0)
		return PROC_SUCCESS;

	nRegionNum = pRegion->nHeight;

	/* create data matrix */
	*ppPeakRawData = CreateDoubleMatrix(nRegionNum, nFileNum);
	if(*ppPeakRawData == NULL)
	{
		printf("Error: SeqPeak_CollectRawData, cannot create matrix for storing normalized raw data!\n");
		exit(EXIT_FAILURE);
	}

	vChrId = (int *)calloc(nRegionNum, sizeof(int));
	if(vChrId == NULL)
	{
		printf("Error: SeqPeak_CollectRawData, cannot create memory for vChrId!\n");
		exit(EXIT_FAILURE);
	}

	vStart = (int *)calloc(nRegionNum, sizeof(int));
	if(vStart == NULL)
	{
		printf("Error: SeqPeak_CollectRawData, cannot create memory for vStart!\n");
		exit(EXIT_FAILURE);
	}

	vEnd = (int *)calloc(nRegionNum, sizeof(int));
	if(vEnd == NULL)
	{
		printf("Error: SeqPeak_CollectRawData, cannot create memory for vEnd!\n");
		exit(EXIT_FAILURE);
	}

	/* get region coordinates */
	for(ni=0; ni<nRegionNum; ni++)
	{
		vChrId[ni] = (int)(DMGETAT(pRegion, ni, 0));
		vStart[ni] = (int)(DMGETAT(pRegion, ni, 1));
		vEnd[ni] = (int)(DMGETAT(pRegion, ni, 2));
	}

	/* process files one by one */
	for(ni=0; ni<nFileNum; ni++)
	{
		SeqPeak_CollectRawData_SingleFile(nRegionNum, 
			vChrId, vStart, vEnd,
			strOutputPath, strOutputFile, 
			vFileName[ni]->m_pString, pFileNormFactor->pMatElement[ni],
			dSmoothFactor, nChrNum, vChrName, pChrLen, 
			nBinSize, *ppPeakRawData, ni);
	}

	/* release memory */
	free(vChrId);
	free(vStart);
	free(vEnd);

	/* return */
	return PROC_SUCCESS;
}

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
		int nCol)
{
	float **vBinC = NULL;
	int ni,nj;
	int *vBinNum;
	int nBinNum;
	char strBinFileName[MED_LINE_LENGTH];
	int nChrId;
	double dSum;

	/* init */
	vBinC = NULL;
	vBinC = (float **)calloc(nChrNum, sizeof(float *));
	if(vBinC == NULL)
	{
		printf("Error: SeqPeak_CollectRawData_SingleFile, cannot create vector for genomic bins!\n");
		exit(EXIT_FAILURE);
	}

	vBinNum = NULL;
	vBinNum = (int *)calloc(nChrNum, sizeof(int));
	if(vBinNum == NULL)
	{
		printf("Error: SeqPeak_CollectRawData_SingleFile, cannot create vector for genomic bin size!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nChrNum; ni++)
	{
		nBinNum = pChrLen->pMatElement[ni]/nBinSize;
		if(pChrLen->pMatElement[ni] % nBinSize != 0)
			nBinNum++;

		vBinC[ni] = (float *)calloc(nBinNum, sizeof(float));
		if(vBinC[ni] == NULL)
		{
			printf("Error: SeqPeak_CollectRawData_SingleFile, insufficient memory for creating genomic bins, try a larger bin size!\n");
			exit(EXIT_FAILURE);
		}		

		vBinNum[ni] = nBinNum;
	}

	/* bring in bin count data */
	for(ni=0; ni<nChrNum; ni++)
	{
		sprintf(strBinFileName, "%s%s_%s_%s.bincount", strOutputPath, strOutputFile, strFileName, vChrName[ni]->m_pString);
		if( TileMapv2_LoadFromBinaryFile(vBinC[ni], sizeof(float), vBinNum[ni], strBinFileName) != PROC_SUCCESS)
		{
			printf("Error: SeqPeak_CollectRawData_SingleFile, cannot read bin statistics correctly!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* collect region data */
	for(ni=0; ni<nRegionNum; ni++)
	{
		dSum = 0.0;
		nChrId = vChrId[ni];
		for(nj=vStart[ni]; nj<=vEnd[ni]; nj++)
		{
			dSum += vBinC[nChrId][nj];
		}
		dSum = dSum*dFileNormFactor/(dSmoothFactor+1e-6);

		DMSETAT(pPeakRawData, ni, nCol, dSum);
	}

	/* release memory */
	for(ni=0; ni<nChrNum; ni++)
	{
		free(vBinC[ni]);
		vBinC[ni] = NULL;
	}
	free(vBinC);
	free(vBinNum);

	/* return */
	return PROC_SUCCESS;
}
