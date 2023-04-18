/* ----------------------------------------------------------------------- */
/*  MicroarrayLib.c : implementation of the microarray library             */
/*  Author : Ji HongKai ; Time: 2004.08                                    */
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
#include "MicroarrayLib.h"
#include "GenomeLib.h"
#include "AffyLib.h"

/* ----------------------------------------------------------------------- */ 
/*  Expression_GetSpecificProbe_Main()                                     */
/*  Get specific probe from the expression data.                           */
/*  The first row of the raw data will be ignored.                         */
/*  Return PROC_SUCCESS if success.                                        */
/* ----------------------------------------------------------------------- */ 
int Expression_GetSpecificProbe_Main(char strDatabasePath[], 
			char strInputPath[], int nColumn, char strOutputPath[])
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	char strHeader[LONG_LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	char strProbeId[LINE_LENGTH];
	char *chp1,*chp2;
	int nj;
	int nPairNum = 0;
	int nFind;
	struct tagStringPair **vDataMap = NULL;

	/* load database */
	strcpy(strHeader, "");
	vDataMap = Affy_LoadDatabase_ExpressionData(strDatabasePath, &nPairNum, strHeader);

	/* init */
	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: Expression_GetSpecificProbe_Main, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: Expression_GetSpecificProbe_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "%s\n", strHeader);

	/* process one by one */
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		/* find refgene id */
		chp1 = strLine;
		for(nj=0; nj<nColumn; nj++)
		{
			chp2 = strchr(chp1, '\t');
			if(chp2 == NULL)
				break;

			chp1 = chp2+1;
		}

		if(nj != nColumn)
		{
			strcpy(strProbeId, "---");	
		}
		else
		{
			sscanf(chp1, "%s", strProbeId);

			if( (strcmp(strProbeId, "---") == 0) || (strcmp(strProbeId, "NA") == 0) )
			{
				strcpy(strProbeId, "---");
			}
			else
			{
			}
		}

		/* search for affy */
		if(strcmp(strProbeId, "---") == 0)
		{
		}
		else
		{
			nFind = 0;
			for(nj=0; nj<nPairNum; nj++)
			{
				if(vDataMap[nj] != NULL)
				{
					if(strcmp(strProbeId, vDataMap[nj]->m_pStr1->m_pString) == 0)
					{
						nFind = 1;
						fprintf(fpOut, "%s\n", vDataMap[nj]->m_pStr2->m_pString);
					}
				}
			}

			if(nFind == 0)
			{
				fprintf(fpOut, "\n");
			}
		}
	}

	/* close files */
	fclose(fpIn);
	fclose(fpOut);

	/* clear memory */
	StringPair_ClearDatabase(&vDataMap, nPairNum);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Expression_GetNonRedundantProbe_Main()                                 */
/*  Get non-redundant probe from the expression data.                      */
/*  The first row of the raw data will be ignored.                         */
/*  Return PROC_SUCCESS if success.                                        */
/* ----------------------------------------------------------------------- */ 
int Expression_GetNonRedundantProbe_Main(char strInputPath[], char strOutputPath[])
{
	/* define */
	FILE *fpOut;
	char strHeader[LONG_LINE_LENGTH];
	int ni,nj;
	int nPairNum = 0;
	struct tagStringPair **vDataMap = NULL;
	struct BYTEMATRIX *pMask;

	/* load database */
	strcpy(strHeader, "");
	vDataMap = Affy_LoadDatabase_ExpressionData(strInputPath, &nPairNum, strHeader);
	if(nPairNum == 0)
	{
		printf("Warning: Expression_GetNonRedundantProbe_Main, null data!\n");
		return PROC_SUCCESS;
	}

	pMask = NULL;
	pMask = CreateByteMatrix(1, nPairNum);
	if(pMask == NULL)
	{
		printf("Error: Expression_GetNonRedundantProbe_Main, cannot create mask matrix!\n");
		exit(EXIT_FAILURE);
	}

	/* trim redundancy */
	for(ni=0; ni<nPairNum; ni++)
	{
		if(vDataMap[ni] == NULL)
		{
			pMask->pMatElement[ni] = 1;
			continue;
		}
		
		for(nj=0; nj<ni; nj++)
		{
			if(pMask->pMatElement[nj] == 0)
			{
				if(strcmp(vDataMap[ni]->m_pStr1->m_pString, vDataMap[nj]->m_pStr1->m_pString) == 0)
				{
					pMask->pMatElement[ni] = 1;
					break;
				}
			}
		}
	}


	/* init */
	fpOut = NULL;
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: Expression_GetNonRedundantProbe_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "%s\n", strHeader);

	/* process one by one */
	for(ni=0; ni<nPairNum; ni++)
	{
		if(pMask->pMatElement[ni] == 0)
		{
			fprintf(fpOut, "%s\n", vDataMap[ni]->m_pStr2->m_pString);
		}
	}
	
	/* close files */
	fclose(fpOut);

	/* clear memory */
	DestroyByteMatrix(pMask);
	StringPair_ClearDatabase(&vDataMap, nPairNum);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneSelection_Main()                                        */
/*  Gene selection based on criteria specified in the criteria file        */
/*  Return PROC_SUCCESS if success.                                        */
/* ----------------------------------------------------------------------- */ 
int Expression_GeneSelection_Main(char strDataFile[], char strGeneInfoFile[],
								  char strCriteriaFile[], char strOutFile[])
{
	/* ------- */
	/* define  */
	/* ------- */

	/* array */
	int nArrayNum;
	int nProbeNum;
	struct DOUBLEMATRIX **vArray;

	/* class id */
	int nClassNum;
	struct INTMATRIX *pClassSize;
	struct INTMATRIX *pClassID;
	struct INTMATRIX *pPermClassID;
	struct INTMATRIX *pDataClassID;

	/* variance group */
	int nVargroupNum;
	struct INTMATRIX *pVargroupSize;
	struct INTMATRIX **vVargroupMap;

	/* permutation group */
	int nCentralize = 1;
	int nPermgroupNum;
	struct INTMATRIX *pPermgroupSize;
	struct INTMATRIX **vPermgroupMap;

	/* files */
	FILE *fpData;
	FILE *fpCriteria;
	FILE *fpProbe;

	/* strings */
	char strFileName[MED_LINE_LENGTH];
	char strLine[MED_LINE_LENGTH];
	char strLongLine[LONG_LINE_LENGTH];
	char strComparisons[MED_LINE_LENGTH];
	int nPairwiseComp;
	
	/* scores */
	struct DOUBLEMATRIX *pScore;
	struct DOUBLEMATRIX *pSortScore;
	struct DOUBLEMATRIX *pPermScore;
	struct DOUBLEMATRIX *pPermSortScore;
	struct LONGMATRIX *pSortID;
	struct DOUBLEMATRIX *pFDR;

	/* pointers */
	double *vP1,*vP2,*vP3;
	int nx,ny;

	/* parameters */
	double dTruncLow;
	int nTakeLog;
	int nFDRPermNum;
	int nCycPermNum;
	int nOutputNum;

	/* error */
	int nError;

	/* count */
	int ni,nj,nk;
	char *chSep,*chSep2;
	double dValue;
	int nTemp;


	/* --------- */
	/* load data */
	/* --------- */
	printf("Expression_GeneSelection_Main:\n");
	printf("loading data...\n");
	
	/* header info */
	nError = 0;
	fpCriteria = NULL;
	fpCriteria = fopen(strCriteriaFile, "rt");
	if(fpCriteria == NULL)
	{
		printf("Error: Expression_GeneSelection_Main, cannot open comparison info file!\n");
		exit(EXIT_FAILURE);
	}
	while(fgets(strLine, MED_LINE_LENGTH, fpCriteria) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		/* load basic info */
		if(strcmp(strLine, "[Basic Info]") == 0)
		{
			/* array number */
			fgets(strLine, MED_LINE_LENGTH, fpCriteria);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			chSep = strstr(strLine, "array");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			nArrayNum = atoi(chSep);

			/* probe number */
			fgets(strLine, MED_LINE_LENGTH, fpCriteria);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			chSep = strstr(strLine, "probeset");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			nProbeNum = atoi(chSep);

			/* class number */
			fgets(strLine, MED_LINE_LENGTH, fpCriteria);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			chSep = strstr(strLine, "group");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			nClassNum = atoi(chSep);
		}

		/* load class info */
		else if(strcmp(strLine, "[Group ID]") == 0)
		{
			fgets(strLongLine, LONG_LINE_LENGTH, fpCriteria);
			pDataClassID = Expression_GeneSelection_LoadGroupId(strLongLine);
			if(pDataClassID == NULL)
			{
				nError = 1;
				break;
			}
		}

		/* load comparisons */
		else if(strcmp(strLine, "[Comparisons]") == 0)
		{
			fgets(strComparisons, MED_LINE_LENGTH, fpCriteria);
			StrTrimLeft(strComparisons);
			StrTrimRight(strComparisons);
			if(strComparisons[0] == '\0')
			{
				nError = 1;
				break;
			}

			nPairwiseComp = 0;
			chSep = strpbrk(strComparisons, "<>" );
			while(chSep != NULL)
			{
				nPairwiseComp++;
				chSep = strpbrk((chSep+1), "<>" );
			}

			if(nPairwiseComp > 1)
				nPairwiseComp = 0;
		}

		/* load preprocessing */
		else if(strcmp(strLine, "[Preprocessing Setup]") == 0)
		{
			/* truncate lower bound */
			fgets(strLine, MED_LINE_LENGTH, fpCriteria);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			chSep = strstr(strLine, "lower");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			dTruncLow = atof(chSep);

			/* take log */
			fgets(strLine, MED_LINE_LENGTH, fpCriteria);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			chSep = strstr(strLine, "log2");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			nTakeLog = atoi(chSep);
		}
		
		/* load output setup */
		else if(strcmp(strLine, "[Output Setup]") == 0)
		{
			/* output number */
			fgets(strLine, MED_LINE_LENGTH, fpCriteria);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			chSep = strstr(strLine, "print");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			nOutputNum = atoi(chSep);
		}

		/* load simulation setup */
		else if(strcmp(strLine, "[Simulation Setup]") == 0)
		{
			/* FDR permutation number */
			fgets(strLine, MED_LINE_LENGTH, fpCriteria);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			chSep = strstr(strLine, "FDR");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			nFDRPermNum = atoi(chSep);

			/* within cycle permutation number */
			fgets(strLine, MED_LINE_LENGTH, fpCriteria);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			chSep = strstr(strLine, "cycle");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			nCycPermNum = atoi(chSep);
		}

		/* permutation group */
		else if(strcmp(strLine, "[Permutation Setup]") == 0)
		{
			/* centralization */
			fgets(strLine, MED_LINE_LENGTH, fpCriteria);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			chSep = strstr(strLine, "centralize");
			if(chSep != NULL)
			{
				chSep = strstr(strLine, "=");
				if(chSep == NULL)
				{
					nError = 1;
					break;
				}
				chSep++;
				StrTrimLeft(chSep);
				nCentralize = atoi(chSep);

				fgets(strLine, MED_LINE_LENGTH, fpCriteria);
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
			}
			else
			{
			}
			
			/* group number */
			chSep = strstr(strLine, "permutation");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			nPermgroupNum = atoi(chSep);

			/* groups */
			pPermgroupSize = NULL;
			pPermgroupSize = CreateIntMatrix(nPermgroupNum, 1);
			if(pPermgroupSize == NULL)
			{
				nError = 1;
				break;
			}
			vPermgroupMap = NULL;
			vPermgroupMap = (struct INTMATRIX **)calloc(nPermgroupNum, sizeof(struct INTMATRIX *));
			if(vPermgroupMap == NULL)
			{
				DestroyIntMatrix(pPermgroupSize);
				nError = 1;
				break;
			}

			for(ni=0; ni<nPermgroupNum; ni++)
			{
				fgets(strLine, MED_LINE_LENGTH, fpCriteria);
				vPermgroupMap[ni] = Expression_GeneSelection_LoadGroupId(strLine);
				if(vPermgroupMap[ni] == NULL)
				{
					nError = 1;
					break;
				}
				IMSETAT(pPermgroupSize, ni, 0, vPermgroupMap[ni]->nWidth);
			}
		}

		/* variance group */
		else if(strcmp(strLine, "[Variance Setup]") == 0)
		{
			/* group number */
			fgets(strLine, MED_LINE_LENGTH, fpCriteria);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			chSep = strstr(strLine, "variance");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			nVargroupNum = atoi(chSep);

			/* groups */
			pVargroupSize = NULL;
			pVargroupSize = CreateIntMatrix(nVargroupNum, 1);
			if(pVargroupSize == NULL)
			{
				nError = 1;
				break;
			}
			vVargroupMap = NULL;
			vVargroupMap = (struct INTMATRIX **)calloc(nVargroupNum, sizeof(struct INTMATRIX *));
			if(vVargroupMap == NULL)
			{
				DestroyIntMatrix(pVargroupSize);
				nError = 1;
				break;
			}

			for(ni=0; ni<nVargroupNum; ni++)
			{
				fgets(strLine, MED_LINE_LENGTH, fpCriteria);
				vVargroupMap[ni] = Expression_GeneSelection_LoadGroupId(strLine);
				if(vVargroupMap[ni] == NULL)
				{
					nError = 1;
					break;
				}
				IMSETAT(pVargroupSize, ni, 0, vVargroupMap[ni]->nWidth);
			}
		}

		/* do nothing */
		else
		{
		}
	}

	fclose(fpCriteria);

	if(nError > 0)
	{
		printf("Error: Expression_GeneSelection_Main, comparison info file format wrong!\n");
		exit(EXIT_FAILURE);
	}

	/* prepare space */
	pClassID = NULL;
	pClassID = CreateIntMatrix(1, nArrayNum);
	pClassSize = NULL;
	pClassSize = CreateIntMatrix(1, nClassNum);
	nj = 0;
	for(ni=0; ni<pDataClassID->nWidth; ni++)
	{
		nk = IMGETAT(pDataClassID, 0, ni);
		if( nk > nClassNum )
		{
			printf("Error: Expression_GeneSelection_Main, group id out of range!\n");
			exit(EXIT_FAILURE);
		}
		if( nk>0 )
		{
			if(nj >= nArrayNum)
			{
				printf("Error: Expression_GeneSelection_Main, group id out of range!\n");
				exit(EXIT_FAILURE);
			}
			IMSETAT(pClassID, 0, nj, nk);
			nTemp = IMGETAT(pClassSize, 0, (nk-1))+1;
			IMSETAT(pClassSize, 0, (nk-1), nTemp);
			nj++;
		}
	}
	if(nj != nArrayNum)
	{
		printf("Error: Expression_GeneSelection_Main, array number not match!\n");
		exit(EXIT_FAILURE);
	}

	vArray = NULL;
	vArray = (struct DOUBLEMATRIX **)calloc(nArrayNum, sizeof(struct DOUBLEMATRIX *));
	if(vArray == NULL)
	{
		printf("Error: Expression_GeneSelection_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nArrayNum; ni++)
	{
		vArray[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vArray[ni] == NULL)
		{
			printf("Error: Expression_GeneSelection_Main, cannot create data space!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* data */
	fpData = NULL;
	fpData = fopen(strDataFile, "rt");
	if(fpData == NULL)
	{
		printf("Error: Expression_GeneSelection_Main, cannot open data file!\n");
		exit(EXIT_FAILURE);
	}
	sprintf(strFileName, "%s.pb", strOutFile);
	fpProbe = NULL;
	fpProbe = fopen(strFileName, "wt");
	if(fpProbe == NULL)
	{
		printf("Error: Expression_GeneSelection_Main, cannot open probe file!\n");
		exit(EXIT_FAILURE);
	}
	fgets(strLongLine, LONG_LINE_LENGTH, fpData);
	nj = 0;
	while(fgets(strLongLine, LONG_LINE_LENGTH, fpData) != NULL)
	{
		StrTrimLeft(strLongLine);
		StrTrimRight(strLongLine);
		if(strLongLine[0] == '\0')
			continue;

		if(nj >= nProbeNum)
		{
			printf("Error: Expression_GeneSelection_Main, probeset number not match!\n");
			exit(EXIT_FAILURE);
		}

		ni = 0;
		nk = 0;
		chSep = strchr(strLongLine, '\t');
		if(chSep != NULL)
			*chSep = '\0';
		fprintf(fpProbe, "%s\n", strLongLine);

		while(chSep != NULL)
		{
			if( nk >= pDataClassID->nWidth )
			{
				printf("Error: Expression_GeneSelection_Main, array number not match!\n");
				exit(EXIT_FAILURE);
			}

			chSep++;
			chSep2 = strchr(chSep, '\t');
		
			if( pDataClassID->pMatElement[nk] > 0 )
			{
				/* middle number */
				if(chSep2 != NULL)
				{
					*chSep2 = '\0';

					if(chSep == chSep2)
					{
						dValue = 0.0;
					}
					else
					{
						dValue = atof(chSep);
					}
				}
				/* last number */
				else
				{
					if(chSep == chSep2)
					{
						dValue = 0.0;
					}
					else
					{
						dValue = atof(chSep);
					}
				}

				if(dValue < dTruncLow)
					dValue = dTruncLow;
				if(nTakeLog == 1)
					dValue = log(dValue)/log(2.0);
			
				DMSETAT(vArray[ni], nj, 0, dValue);
				ni++;
			}

			/* get next */
			nk++;
			chSep = chSep2;
		}

		/* get next line */
		nj++;
	}

	fclose(fpData);
	fclose(fpProbe);
	DestroyIntMatrix(pDataClassID);

	if(nj != nProbeNum)
	{
		printf("Error: Expression_GeneSelection_Main, probeset number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* -------------------- */
	/* get original ranking */
	/* -------------------- */
	printf("rank genes...\n");
	pScore = NULL;
	pSortScore = NULL;
	pSortID = NULL;
	if(nPairwiseComp == 1)
	{
		pScore = Expression_GeneSelection_tTest(nProbeNum, nArrayNum, vArray,
									nClassNum, pClassSize, pClassID,
									nVargroupNum, pVargroupSize, vVargroupMap, 
									strComparisons);
	}
	else
	{
		pScore = Expression_GeneSelection_MonteCarlo(nProbeNum, nArrayNum, vArray,
									nClassNum, pClassSize, pClassID,
									nVargroupNum, pVargroupSize, vVargroupMap, 
									nCycPermNum, strComparisons);
	}

	if( DMSORTMERGEA_0(pScore, &pSortScore, &pSortID) == PROC_FAILURE )
	{
		printf("Error: Expression_GeneSelection_Main, sorting score failure!\n");
		exit(EXIT_FAILURE);
	}

	/* ---------------------------- */
	/* permutations to estimate FDR */
	/* ---------------------------- */
	pFDR = NULL;
	if(nFDRPermNum > 0)
	{
		/* centralization */
		if(nCentralize == 1)
		{
			Expression_GeneSelection_SubtractMean(nProbeNum, 
					nArrayNum, vArray, nClassNum, pClassSize, 
					pClassID);
		}

		printf("estimate FDR...\n");
		pFDR = CreateDoubleMatrix(1, nProbeNum);

		/* cycles */
		for(ni=0; ni<nFDRPermNum; ni++)
		{
			printf("perm %d...\n", ni);
			/* permute the class label */
			pPermClassID = NULL;
			pPermClassID = Expression_GeneSelection_ClassPerm(nClassNum, pClassID, pClassSize,
				nPermgroupNum, pPermgroupSize, vPermgroupMap);

			/* rerun a single cycle */
			pPermScore = NULL;
			pPermSortScore = NULL;
			if(nPairwiseComp == 1)
			{
				pPermScore = Expression_GeneSelection_tTest(nProbeNum, nArrayNum, vArray,
									nClassNum, pClassSize, pPermClassID,
									nVargroupNum, pVargroupSize, vVargroupMap, 
									strComparisons);
			}
			else
			{
				pPermScore = Expression_GeneSelection_MonteCarlo(nProbeNum, nArrayNum, vArray,
									nClassNum, pClassSize, pPermClassID,
									nVargroupNum, pVargroupSize, vVargroupMap, 
									nCycPermNum, strComparisons);
			}

			/* sort scores */
			if( DMSORTMERGEA_0(pPermScore, &pPermSortScore, NULL) == PROC_FAILURE )
			{
				printf("Error: Expression_GeneSelection_Main, sorting score failure!\n");
				exit(EXIT_FAILURE);
			}

			/* add FDR count */
			vP1 = pFDR->pMatElement;
			vP2 = pSortScore->pMatElement;
			vP3 = pPermSortScore->pMatElement;
			ny = 0;
			for(nx=0; nx<nProbeNum; nx++)
			{
				for(; ny<nProbeNum; ny++)
				{
					if(vP3[ny] > vP2[nx])
						break;
				}

				if(ny >= nProbeNum)
				{
					vP1[nx] += (double)nProbeNum;
				}
				else
				{
					vP1[nx] += (double)ny;
				}
			}

			/* release memory */
			DestroyIntMatrix(pPermClassID);
			DestroyDoubleMatrix(pPermScore);
			DestroyDoubleMatrix(pPermSortScore);
		}

		/* normalize FDR */
		vP1 = pFDR->pMatElement;
		for(ni=0; ni<nProbeNum; ni++)
		{
			vP1[ni] = vP1[ni]/(double)(nFDRPermNum*(ni+1));
		}
		for(ni=nProbeNum-1; ni>0; ni--)
		{
			if(vP1[ni-1] > vP1[ni])
				vP1[ni-1] = vP1[ni];
		}
	}

	/* -------------- */
	/* release memory */
	/* -------------- */

	/* destroy arrays */
	for(ni=0; ni<nArrayNum; ni++)
	{
		DestroyDoubleMatrix(vArray[ni]);
		vArray[ni] = NULL;
	}
	free(vArray);

	/* destroy class ids */
	DestroyIntMatrix(pClassSize);
	DestroyIntMatrix(pClassID);

	/* destroy variance group ids */
	for(ni=0; ni<nVargroupNum; ni++)
	{
		DestroyIntMatrix(vVargroupMap[ni]);
		vVargroupMap[ni] = NULL;
	}
	free(vVargroupMap);
	DestroyIntMatrix(pVargroupSize);

	/* destroy permutation group ids */
	for(ni=0; ni<nPermgroupNum; ni++)
	{
		DestroyIntMatrix(vPermgroupMap[ni]);
		vPermgroupMap[ni] = NULL;
	}
	free(vPermgroupMap);
	DestroyIntMatrix(pPermgroupSize);
	
	/* -------------------------------- */
	/* add gene info to the output file */
	/* -------------------------------- */
	printf("output top genes and link gene information...\n");
	sprintf(strFileName, "%s.pb", strOutFile);
	Expression_GeneSelection_Output_WithRandomControl_Fast(nProbeNum, nOutputNum, 
		pScore, pSortID, pSortScore, pFDR,
		strFileName, strGeneInfoFile, strOutFile);

	/* destroy results */
	DestroyDoubleMatrix(pScore);
	DestroyDoubleMatrix(pFDR);
	DestroyDoubleMatrix(pSortScore);
	DestroyLongMatrix(pSortID);

	/* ------ */
	/* return */
	/* ------ */
	printf("Done!\n");
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneSelection_LoadGroupId()                                 */
/*  Load group ids                                                         */
/* ----------------------------------------------------------------------- */ 
struct INTMATRIX *Expression_GeneSelection_LoadGroupId(char strInLine[])
{
	/* define */
	int nElement[IM_LOAD_MAX_WIDTH];
	struct INTMATRIX *pMat;
	int nElementCount;

	/* init */
	pMat = NULL;

	/* load */
	nElementCount = IntLoadRowVector(strInLine, nElement);
	if(nElementCount == 0)
		return NULL;	

	if(IMADDROW(&pMat, nElement, nElementCount) == PROC_FAILURE)
	{
		DestroyIntMatrix(pMat);
		pMat = NULL;
	}

	/* return */
	return pMat;
}

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
					int nCycPermNum, char strComparisons[])
{
	/* ------ */
	/* define */
	/* ------ */
	struct INTMATRIX *pClassSizeCopy;
	struct INTMATRIX *pDf;
	struct DOUBLEMATRIX *pScore;
	struct DOUBLEMATRIX **vMeans;
	struct DOUBLEMATRIX **vVars;
	struct DOUBLEMATRIX **vSDs;
	/* struct DOUBLEMATRIX *pMeanDraws; */
	struct DOUBLEMATRIX **vMeanDraws;
	struct DOUBLEMATRIX *pSigCoef;
	double *vExp,*vSum,*vAve,*vSigma,*vMu;
	unsigned char *vEval;
	int nClustId;
	double dDenom,dTemp;
	double dVarmean,dVarsst;
	double dB,dV2,dK,dN;

	char vLogic[MED_LINE_LENGTH];
	char vTNumber[LINE_LENGTH];
	double vGid[MED_LINE_LENGTH];
	struct DOUBLEMATRIX *vVid[MED_LINE_LENGTH];
	struct BYTEMATRIX *vLid[MED_LINE_LENGTH];
	int nLogicLen;
	int nLeftNum,nRightNum,nTNumLen;
	int nTrue;
	struct BYTEMATRIX *pTrueVec;
	struct DOUBLEMATRIX *pTrueVal;

	int ni,nj,nk,nLen;

	/* ------------- */
	/* expression    */
	/* ------------- */
	nLogicLen = 0;
	StrTrimLeft(strComparisons);
	StrTrimRight(strComparisons);
	if(strComparisons[0] == '\0')
	{
		printf("Error: Expression_GeneSelection_MonteCarlo, no logic expression formula available!\n");
		exit(EXIT_FAILURE);
	}
	nLen = strlen(strComparisons);
	nj = 0;
	nLeftNum = 0;
	nRightNum = 0;
	nTNumLen = 0;
	for(ni=0; ni<nLen; ni++)
	{
		if( strComparisons[ni] == '(' )
		{
			if(nTNumLen != 0)
			{
				printf("Error: Expression_GeneSelection_MonteCarlo, logic expression wrong!\n");
				exit(EXIT_FAILURE);
			}

			vLogic[nj] = strComparisons[ni];
			vGid[nj] = -1.0;
			nj++;
			nLeftNum++;
		}
		else if(strComparisons[ni] == ')')
		{
			if(nTNumLen > 0)
			{
				vTNumber[nTNumLen] = '\0';
				vLogic[nj] = 'G';
				vGid[nj] = atof(vTNumber);
				nTNumLen = 0;
				nj++;
			}
			vLogic[nj] = strComparisons[ni];
			vGid[nj] = -1.0;
			nj++;
			nRightNum++;
		}
		else if( (strComparisons[ni] == '<') || (strComparisons[ni] == '>') 
			|| (strComparisons[ni] == '&') || (strComparisons[ni] == '|') )
		{
			if(nTNumLen > 0)
			{
				vTNumber[nTNumLen] = '\0';
				vLogic[nj] = 'G';
				vGid[nj] = atof(vTNumber);
				nTNumLen = 0;
				nj++;
			}
			vLogic[nj] = strComparisons[ni];
			vGid[nj] = -1.0;
			nj++;
			
		}
		else if( (strComparisons[ni] >= '0') && (strComparisons[ni] <= '9') )
		{
			vTNumber[nTNumLen] = strComparisons[ni];
			nTNumLen++;
		}
		else if( (strComparisons[ni] == ' ') || (strComparisons[ni] == '\t') )
		{
			/* do nothing */
		}
		else
		{
			printf("Error: Expression_GeneSelection_MonteCarlo, %c not supported in logic expressions!\n", strComparisons[ni]);
			exit(EXIT_FAILURE);
		}
	}
	if(nTNumLen > 0)
	{
		vTNumber[nTNumLen] = '\0';
		vLogic[nj] = 'G';
		vGid[nj] = atof(vTNumber);
		nTNumLen = 0;
		nj++;
	}
	nLogicLen = nj;

	if(nLeftNum != nRightNum)
	{
		printf("Error: Expression_GeneSelection_MonteCarlo, lack '(' or ')'!\n");
		exit(EXIT_FAILURE);
	}

	/* for debug purpose*/
	/* for(ni=0; ni<nLogicLen; ni++)
	{
		printf("%c %d\n", vLogic[ni], (int)(vGid[ni]));
	} */
	
	/* ---- */
	/* init */
	/* ---- */
	pScore = NULL;
	pScore = CreateDoubleMatrix(1, nProbeNum);
	pClassSizeCopy = NULL;
	pClassSizeCopy = CreateIntMatrix(1, nClassNum);
	pDf = NULL;
	pDf = CreateIntMatrix(1, nVargroupNum);


	/* ------------- */
	/* create memory */
	/* ------------- */
	vMeans = NULL;
	vMeans = (struct DOUBLEMATRIX **)calloc(nClassNum, sizeof(struct DOUBLEMATRIX *));
	if(vMeans == NULL)
	{
		printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
		exit(EXIT_FAILURE);
	}
	vVars = NULL;
	vVars = (struct DOUBLEMATRIX **)calloc(nClassNum, sizeof(struct DOUBLEMATRIX *));
	if(vVars == NULL)
	{
		printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
		exit(EXIT_FAILURE);
	}
	vSDs = NULL;
	vSDs = (struct DOUBLEMATRIX **)calloc(nVargroupNum, sizeof(struct DOUBLEMATRIX *));
	if(vSDs == NULL)
	{
		printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
		exit(EXIT_FAILURE);
	}
	vMeanDraws = NULL;
	vMeanDraws = (struct DOUBLEMATRIX **)calloc(nClassNum, sizeof(struct DOUBLEMATRIX *));
	if(vMeanDraws == NULL)
	{
		printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
		exit(EXIT_FAILURE);
	}
	
	for(ni=0; ni<nClassNum; ni++)
	{
		vMeans[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vMeans[ni] == NULL)
		{
			printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
			exit(EXIT_FAILURE);
		}

		vVars[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vVars[ni] == NULL)
		{
			printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
			exit(EXIT_FAILURE);
		}

		vMeanDraws[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vMeanDraws[ni] == NULL)
		{
			printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	for(ni=0; ni<nVargroupNum; ni++)
	{
		vSDs[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vSDs[ni] == NULL)
		{
			printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation!\n");
			exit(EXIT_FAILURE);
		}
	}

	/*pMeanDraws = NULL;
	pMeanDraws = CreateDoubleMatrix(1,nClassNum);
	*/

	/* ---------------- */
	/* get mean and var */
	/* ---------------- */

	/* mean */
	for(ni=0; ni<nArrayNum; ni++)
	{
		nClustId = (pClassID->pMatElement)[ni]-1;
		(pClassSizeCopy->pMatElement)[nClustId] += 1;
		vSum = vMeans[nClustId]->pMatElement;
		vExp = vArray[ni]->pMatElement;
		for(nj=0; nj<nProbeNum; nj++)
		{
			vSum[nj] += vExp[nj];
		}
	}
	
	for(ni=0; ni<nClassNum; ni++)
	{
		if((pClassSizeCopy->pMatElement)[ni] != (pClassSize->pMatElement)[ni])
		{
			printf("Error: Expression_GeneSelection_MonteCarlo, class size not match!\n");
			exit(EXIT_FAILURE);
		}
		if((pClassSizeCopy->pMatElement)[ni] > 0)
		{
			dDenom = (double)((pClassSizeCopy->pMatElement)[ni]);
			vSum = vMeans[ni]->pMatElement;
			for(nj=0; nj<nProbeNum; nj++)
			{
				vSum[nj] /= dDenom;
			}
		}
	}

	/* variance: sum of squares */
	for(ni=0; ni<nArrayNum; ni++)
	{
		nClustId = (pClassID->pMatElement)[ni]-1;
		vSum = vVars[nClustId]->pMatElement;
		vAve = vMeans[nClustId]->pMatElement;
		vExp = vArray[ni]->pMatElement;
		for(nj=0; nj<nProbeNum; nj++)
		{
			dTemp = (vExp[nj]-vAve[nj]);		
			vSum[nj] += dTemp*dTemp;
		}
	}

	/* variance: estimates */
	for(ni=0; ni<nVargroupNum; ni++)
	{
		if(vVargroupMap[ni]->nWidth != pVargroupSize->pMatElement[ni])
		{
			printf("Error: Expression_GeneSelection_MonteCarlo, variance group number not match!\n");
			exit(EXIT_FAILURE);
		}

		for(nj=0; nj<vVargroupMap[ni]->nWidth; nj++)
		{
			nClustId = (vVargroupMap[ni]->pMatElement)[nj]-1;
			if((pClassSize->pMatElement)[nClustId] > 0)
				(pDf->pMatElement)[ni] += (pClassSize->pMatElement)[nClustId]-1;
			vSum = vSDs[ni]->pMatElement;
			vExp = vVars[nClustId]->pMatElement;
			for(nk=0; nk<nProbeNum; nk++)
			{
				vSum[nk] += vExp[nk];
			}
		}
	}
	

	/* ------------- */
	/* shrinking var */
	/* ------------- */
	for(ni=0; ni<nVargroupNum; ni++)
	{
		dVarmean = 0.0;
		dVarsst = 0.0;
		dB = 0.0;

		if((pDf->pMatElement)[ni] > 0)
		{
			dDenom = (double)(pDf->pMatElement)[ni];
			vSum = vSDs[ni]->pMatElement;

			/* xbar */
			for(nj=0; nj<nProbeNum; nj++)
			{
				vSum[nj] /= dDenom;
				dVarmean += vSum[nj];
			}
			dVarmean /= (double)nProbeNum;

			/* sst */
			for(nj=0; nj<nProbeNum; nj++)
			{
				dTemp = vSum[nj]-dVarmean;
				dVarsst += dTemp*dTemp;
			}

			/* shrinkage factor */
			dV2 = 2.0/dDenom;
			dK = (double)nProbeNum;
			dN = 1.0;
			if((dK > 0.0) && (dVarsst > 0.0))
				dB = dV2*(dK-1)/((dN+dV2)*dK)+dN*(dK-1)*dV2*dVarmean*dVarmean/((dN+dV2)*dN*dVarsst);
				/* dB = 0.0; */
			else
				dB = 0.0;

			if(dB < 0.0)
				dB = 0.0;
			else if(dB > 1.0)
				dB = 1.0;

			/* shrink variance */
			for(nj=0; nj<nProbeNum; nj++)
			{
				vSum[nj] = (1-dB)*vSum[nj]+dB*dVarmean;
				vSum[nj] = sqrt(vSum[nj])+1e-16;
			}
		}
	}

	for(ni=0; ni<nClassNum; ni++)
	{
		DestroyDoubleMatrix(vVars[ni]);
		vVars[ni] = NULL;
	}

	for(ni=0; ni<nVargroupNum; ni++)
	{
		for(nj=0; nj<vVargroupMap[ni]->nWidth; nj++)
		{
			nClustId = (vVargroupMap[ni]->pMatElement)[nj]-1;
			vVars[nClustId] = vSDs[ni];
		}
	}

	/* ----------- */
	/* simulation  */
	/* ----------- */
	pSigCoef = NULL;
	pSigCoef = CreateDoubleMatrix(1, nClassNum);
	vExp = pSigCoef->pMatElement;
	for(ni=0; ni<nClassNum; ni++)
	{
		if((pClassSize->pMatElement)[ni] > 0)
		{
			dDenom = (double)((pClassSize->pMatElement)[ni]);
			vExp[ni] = sqrt(1.0/dDenom);
		}
		else
		{
			vExp[ni] = 0.0;
		}
	}

	for(ni=0; ni<nLogicLen; ni++)
	{
		if(vGid[ni] > 0.0)
		{
			nk = (int)(vGid[ni]-1.0);
			vVid[ni] = vMeanDraws[nk];
		}
		else
		{
			vVid[ni] = NULL;
		}
		vLid[ni] = NULL;
	}

	for(ni=0; ni<nCycPermNum; ni++)
	{
		if(ni%100 == 0)
		{
			printf("iter %d...\n", ni);
		}

		vExp = pSigCoef->pMatElement;
		
		/* probeset by probeset */
		for(nk=0; nk<nClassNum; nk++)
		{
			vAve = vMeanDraws[nk]->pMatElement;
			vMu = vMeans[nk]->pMatElement;
			vSigma = vVars[nk]->pMatElement;
			dTemp = normrnd(0.0, 1.0);
			for(nj=0; nj<nProbeNum; nj++)
			{
				/* simulate */
				vAve[nj] = vMu[nj]+vExp[nk]*vSigma[nj]*dTemp;
			}
		}

		

		/* evaluate */
		pTrueVec = NULL;
		pTrueVal = NULL;
		nTrue = Expression_GeneSelection_EvaluateVec(nProbeNum, vLogic, vVid, vLid, nLogicLen, 0, &pTrueVec, &pTrueVal);
		if(nTrue != 4)
		{
			printf("Error: logic evaluation wrong!\n");
			exit(EXIT_FAILURE);
		}

		/* nTrue = Expression_GeneSelection_Evaluate(pMeanDraws, vLogic, vGid, nLogicLen, 0); */

		/* add */
		vSum = pScore->pMatElement;
		/* vSum[nj] += (double)(nTrue); */
		
		vEval = pTrueVec->pMatElement;
		for(nj=0; nj<nProbeNum; nj++)
		{
			vSum[nj] += (double)(vEval[nj]);
		}
		DestroyByteMatrix(pTrueVec);
	}
	DestroyDoubleMatrix(pSigCoef);

	/* normalize */
	dDenom = (double)nCycPermNum;
	vSum = pScore->pMatElement;
	for(nj=0; nj<nProbeNum; nj++)
	{
		vSum[nj] = 1.0-vSum[nj]/dDenom;
	}

	/* -------------- */
	/* release memory */
	/* -------------- */
	for(ni=0; ni<nClassNum; ni++)
	{
		DestroyDoubleMatrix(vMeans[ni]);
		vMeans[ni] = NULL;
		
		DestroyDoubleMatrix(vMeanDraws[ni]);
		vMeanDraws[ni] = NULL;
	}

	for(ni=0; ni<nVargroupNum; ni++)
	{
		DestroyDoubleMatrix(vSDs[ni]);
		vSDs[ni] = NULL;		
	}

	free(vMeans);
	free(vVars);
	free(vSDs);
	free(vMeanDraws);

	/* DestroyDoubleMatrix(pMeanDraws); */
	DestroyIntMatrix(pClassSizeCopy);
	DestroyIntMatrix(pDf);

	/* ------ */
	/* return */
	/* ------ */
	return pScore;
}

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
					char strComparisons[])
{
	/* ------ */
	/* define */
	/* ------ */
	struct INTMATRIX *pClassSizeCopy;
	struct INTMATRIX *pDf;
	struct DOUBLEMATRIX *pScore;
	struct DOUBLEMATRIX **vMeans;
	struct DOUBLEMATRIX **vVars;
	struct DOUBLEMATRIX **vSDs;
	
	struct DOUBLEMATRIX *pSigCoef;
	double *vExp,*vSum,*vAve,*vSigma,*vMu,*vSigma2,*vMu2;
	double dExp,dExp2,dVarTemp;
	int nClustId;
	double dDenom,dTemp;
	double dVarmean,dVarsst;
	double dB,dV2,dK,dN;

	char vLogic[MED_LINE_LENGTH];
	char vTNumber[LINE_LENGTH];
	double vGid[MED_LINE_LENGTH];
	
	int nLogicLen;
	int nLeftNum,nRightNum,nTNumLen;

	int ni,nj,nk,nLen;

	int ng1,ng2;

	/* ------------- */
	/* expression    */
	/* ------------- */
	nLogicLen = 0;
	StrTrimLeft(strComparisons);
	StrTrimRight(strComparisons);
	if(strComparisons[0] == '\0')
	{
		printf("Error: Expression_GeneSelection_t-Test, no logic expression formula available!\n");
		exit(EXIT_FAILURE);
	}
	nLen = strlen(strComparisons);
	nj = 0;
	nLeftNum = 0;
	nRightNum = 0;
	nTNumLen = 0;
	for(ni=0; ni<nLen; ni++)
	{
		if( strComparisons[ni] == '(' )
		{
			if(nTNumLen != 0)
			{
				printf("Error: Expression_GeneSelection_t-Test, logic expression wrong!\n");
				exit(EXIT_FAILURE);
			}

			nLeftNum++;
		}
		else if(strComparisons[ni] == ')')
		{
			if(nTNumLen > 0)
			{
				vTNumber[nTNumLen] = '\0';
				vLogic[nj] = 'G';
				vGid[nj] = atof(vTNumber);
				nTNumLen = 0;
				nj++;
			}
			nRightNum++;
		}
		else if( (strComparisons[ni] == '<') || (strComparisons[ni] == '>') 
			|| (strComparisons[ni] == '&') || (strComparisons[ni] == '|') )
		{
			if(nTNumLen > 0)
			{
				vTNumber[nTNumLen] = '\0';
				vLogic[nj] = 'G';
				vGid[nj] = atof(vTNumber);
				nTNumLen = 0;
				nj++;
			}
			vLogic[nj] = strComparisons[ni];
			vGid[nj] = -1.0;
			nj++;	
		}
		else if( (strComparisons[ni] >= '0') && (strComparisons[ni] <= '9') )
		{
			vTNumber[nTNumLen] = strComparisons[ni];
			nTNumLen++;
		}
		else if( (strComparisons[ni] == ' ') || (strComparisons[ni] == '\t') )
		{
			/* do nothing */
		}
		else
		{
			printf("Error: Expression_GeneSelection_t-Test, %c not supported in logic expressions!\n", strComparisons[ni]);
			exit(EXIT_FAILURE);
		}
	}
	if(nTNumLen > 0)
	{
		vTNumber[nTNumLen] = '\0';
		vLogic[nj] = 'G';
		vGid[nj] = atof(vTNumber);
		nTNumLen = 0;
		nj++;
	}
	nLogicLen = nj;

	if(nLeftNum != nRightNum)
	{
		printf("Error: Expression_GeneSelection_t-Test, lack '(' or ')'!\n");
		exit(EXIT_FAILURE);
	}

	if(nLogicLen != 3)
	{
		printf("Error: Expression_GeneSelection_t-Test, logic expression wrong!\n");
		exit(EXIT_FAILURE);
	}
	if((vLogic[0] != 'G') || (vLogic[2] != 'G'))
	{
		printf("Error: Expression_GeneSelection_t-Test, logic expression wrong!\n");
		exit(EXIT_FAILURE);
	}
	if((vLogic[1] != '<') && (vLogic[1] != '>'))
	{
		printf("Error: Expression_GeneSelection_t-Test, logic expression wrong!\n");
		exit(EXIT_FAILURE);
	}

	if(vLogic[1] == '<')
	{
		ng1 = (int)vGid[0]-1;
		ng2 = (int)vGid[2]-1;
	}
	else if(vLogic[1] == '>')
	{
		ng1 = (int)vGid[2]-1;
		ng2 = (int)vGid[0]-1;
	}
	else
	{
		printf("Error: Expression_GeneSelection_t-Test, logic expression wrong!\n");
		exit(EXIT_FAILURE);
	}

	if((ng1<0) || (ng1>=nClassNum) || (ng2<0) || (ng2>=nClassNum))
	{
		printf("Error: Expression_GeneSelection_t-Test, group id out of range!\n");
		exit(EXIT_FAILURE);
	}

	/* for debug purpose*/
	/* for(ni=0; ni<nLogicLen; ni++)
	{
		printf("%c %d\n", vLogic[ni], (int)(vGid[ni]));
	} */
	
	/* ---- */
	/* init */
	/* ---- */
	pScore = NULL;
	pScore = CreateDoubleMatrix(1, nProbeNum);
	pClassSizeCopy = NULL;
	pClassSizeCopy = CreateIntMatrix(1, nClassNum);
	pDf = NULL;
	pDf = CreateIntMatrix(1, nVargroupNum);


	/* ------------- */
	/* create memory */
	/* ------------- */
	vMeans = NULL;
	vMeans = (struct DOUBLEMATRIX **)calloc(nClassNum, sizeof(struct DOUBLEMATRIX *));
	if(vMeans == NULL)
	{
		printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
		exit(EXIT_FAILURE);
	}
	vVars = NULL;
	vVars = (struct DOUBLEMATRIX **)calloc(nClassNum, sizeof(struct DOUBLEMATRIX *));
	if(vVars == NULL)
	{
		printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
		exit(EXIT_FAILURE);
	}
	vSDs = NULL;
	vSDs = (struct DOUBLEMATRIX **)calloc(nVargroupNum, sizeof(struct DOUBLEMATRIX *));
	if(vSDs == NULL)
	{
		printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
		exit(EXIT_FAILURE);
	}
		
	for(ni=0; ni<nClassNum; ni++)
	{
		vMeans[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vMeans[ni] == NULL)
		{
			printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
			exit(EXIT_FAILURE);
		}

		vVars[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vVars[ni] == NULL)
		{
			printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	for(ni=0; ni<nVargroupNum; ni++)
	{
		vSDs[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vSDs[ni] == NULL)
		{
			printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* ---------------- */
	/* get mean and var */
	/* ---------------- */

	/* mean */
	for(ni=0; ni<nArrayNum; ni++)
	{
		nClustId = (pClassID->pMatElement)[ni]-1;
		(pClassSizeCopy->pMatElement)[nClustId] += 1;
		vSum = vMeans[nClustId]->pMatElement;
		vExp = vArray[ni]->pMatElement;
		for(nj=0; nj<nProbeNum; nj++)
		{
			vSum[nj] += vExp[nj];
		}
	}
	
	for(ni=0; ni<nClassNum; ni++)
	{
		if((pClassSizeCopy->pMatElement)[ni] != (pClassSize->pMatElement)[ni])
		{
			printf("Error: Expression_GeneSelection_t-Test, class size not match!\n");
			exit(EXIT_FAILURE);
		}
		if((pClassSizeCopy->pMatElement)[ni] > 0)
		{
			dDenom = (double)((pClassSizeCopy->pMatElement)[ni]);
			vSum = vMeans[ni]->pMatElement;
			for(nj=0; nj<nProbeNum; nj++)
			{
				vSum[nj] /= dDenom;
			}
		}
	}

	/* variance: sum of squares */
	for(ni=0; ni<nArrayNum; ni++)
	{
		nClustId = (pClassID->pMatElement)[ni]-1;
		vSum = vVars[nClustId]->pMatElement;
		vAve = vMeans[nClustId]->pMatElement;
		vExp = vArray[ni]->pMatElement;
		for(nj=0; nj<nProbeNum; nj++)
		{
			dTemp = (vExp[nj]-vAve[nj]);		
			vSum[nj] += dTemp*dTemp;
		}
	}

	/* variance: estimates */
	for(ni=0; ni<nVargroupNum; ni++)
	{
		if(vVargroupMap[ni]->nWidth != pVargroupSize->pMatElement[ni])
		{
			printf("Error: Expression_GeneSelection_t-Test, variance group number not match!\n");
			exit(EXIT_FAILURE);
		}

		for(nj=0; nj<vVargroupMap[ni]->nWidth; nj++)
		{
			nClustId = (vVargroupMap[ni]->pMatElement)[nj]-1;
			if((pClassSize->pMatElement)[nClustId] > 0)
				(pDf->pMatElement)[ni] += (pClassSize->pMatElement)[nClustId]-1;
			vSum = vSDs[ni]->pMatElement;
			vExp = vVars[nClustId]->pMatElement;
			for(nk=0; nk<nProbeNum; nk++)
			{
				vSum[nk] += vExp[nk];
			}
		}
	}
	

	/* ------------- */
	/* shrinking var */
	/* ------------- */
	for(ni=0; ni<nVargroupNum; ni++)
	{
		dVarmean = 0.0;
		dVarsst = 0.0;
		dB = 0.0;

		if((pDf->pMatElement)[ni] > 0)
		{
			dDenom = (double)(pDf->pMatElement)[ni];
			vSum = vSDs[ni]->pMatElement;

			/* xbar */
			for(nj=0; nj<nProbeNum; nj++)
			{
				vSum[nj] /= dDenom;
				dVarmean += vSum[nj];
			}
			dVarmean /= (double)nProbeNum;

			/* sst */
			for(nj=0; nj<nProbeNum; nj++)
			{
				dTemp = vSum[nj]-dVarmean;
				dVarsst += dTemp*dTemp;
			}

			/* shrinkage factor */
			dV2 = 2.0/dDenom;
			dK = (double)nProbeNum;
			dN = 1.0;
			if((dK > 0.0) && (dVarsst > 0.0))
				dB = dV2*(dK-1)/((dN+dV2)*dK)+dN*(dK-1)*dV2*dVarmean*dVarmean/((dN+dV2)*dN*dVarsst);
				/* dB = 0.0; */
			else
				dB = 0.0;

			if(dB < 0.0)
				dB = 0.0;
			else if(dB > 1.0)
				dB = 1.0;

			/* shrink variance */
			for(nj=0; nj<nProbeNum; nj++)
			{
				vSum[nj] = (1-dB)*vSum[nj]+dB*dVarmean+1e-16;
			}
		}
	}

	for(ni=0; ni<nClassNum; ni++)
	{
		DestroyDoubleMatrix(vVars[ni]);
		vVars[ni] = NULL;
	}

	for(ni=0; ni<nVargroupNum; ni++)
	{
		for(nj=0; nj<vVargroupMap[ni]->nWidth; nj++)
		{
			nClustId = (vVargroupMap[ni]->pMatElement)[nj]-1;
			vVars[nClustId] = vSDs[ni];
		}
	}

	/* ----------- */
	/* simulation  */
	/* ----------- */
	pSigCoef = NULL;
	pSigCoef = CreateDoubleMatrix(1, nClassNum);
	vExp = pSigCoef->pMatElement;
	for(ni=0; ni<nClassNum; ni++)
	{
		if((pClassSize->pMatElement)[ni] > 0)
		{
			dDenom = (double)((pClassSize->pMatElement)[ni]);
			vExp[ni] = 1.0/dDenom;
		}
		else
		{
			vExp[ni] = 0.0;
		}
	}


	vMu = vMeans[ng1]->pMatElement;
	vSigma = vVars[ng1]->pMatElement;
	vMu2 = vMeans[ng2]->pMatElement;
	vSigma2 = vVars[ng2]->pMatElement;
	dExp = pSigCoef->pMatElement[ng1];
	dExp2 = pSigCoef->pMatElement[ng2];
	vSum = pScore->pMatElement;

	/* probeset by probeset */
	for(nj=0; nj<nProbeNum; nj++)
	{
		if((vSigma[nj] <= 0.0) || (vSigma2[nj] <= 0.0))
		{
			printf("Warning: Expression_GeneSelection_t-Test, variance=0, may not have enough sample to estimate variance!\n");
		}
		dVarTemp = sqrt(vSigma[nj]*dExp+vSigma2[nj]*dExp2);
		if(dVarTemp > 0.0)
			vSum[nj] = (vMu[nj]-vMu2[nj])/dVarTemp;
		else
			vSum[nj] = 0.0;
	}
	
	DestroyDoubleMatrix(pSigCoef);

	/* -------------- */
	/* release memory */
	/* -------------- */
	for(ni=0; ni<nClassNum; ni++)
	{
		DestroyDoubleMatrix(vMeans[ni]);
		vMeans[ni] = NULL;
	}

	for(ni=0; ni<nVargroupNum; ni++)
	{
		DestroyDoubleMatrix(vSDs[ni]);
		vSDs[ni] = NULL;		
	}

	free(vMeans);
	free(vVars);
	free(vSDs);

	DestroyIntMatrix(pClassSizeCopy);
	DestroyIntMatrix(pDf);

	/* ------ */
	/* return */
	/* ------ */
	return pScore;
}


/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneSelection_SubtractMean()                                */
/*  Get a score for every gene.                                            */
/* ----------------------------------------------------------------------- */ 
int Expression_GeneSelection_SubtractMean(int nProbeNum, 
					int nArrayNum, struct DOUBLEMATRIX **vArray,
					int nClassNum, struct INTMATRIX *pClassSize, 
					struct INTMATRIX *pClassID)
{
	/* ------ */
	/* define */
	/* ------ */
	struct INTMATRIX *pClassSizeCopy;
	struct DOUBLEMATRIX **vMeans;
	int nClustId;
	int ni,nj;
	double *vExp,*vSum;
	double dDenom;
		
	/* ---- */
	/* init */
	/* ---- */
	pClassSizeCopy = NULL;
	pClassSizeCopy = CreateIntMatrix(1, nClassNum);
	
	/* ------------- */
	/* create memory */
	/* ------------- */
	vMeans = NULL;
	vMeans = (struct DOUBLEMATRIX **)calloc(nClassNum, sizeof(struct DOUBLEMATRIX *));
	if(vMeans == NULL)
	{
		printf("Error: Expression_GeneSelection_SubtractMean, cannot create space for simulation.!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nClassNum; ni++)
	{
		vMeans[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vMeans[ni] == NULL)
		{
			printf("Error: Expression_GeneSelection_SubtractMean, cannot create space for simulation.!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* ---------------- */
	/* get mean and var */
	/* ---------------- */

	/* mean */
	for(ni=0; ni<nArrayNum; ni++)
	{
		nClustId = (pClassID->pMatElement)[ni]-1;
		(pClassSizeCopy->pMatElement)[nClustId] += 1;
		vSum = vMeans[nClustId]->pMatElement;
		vExp = vArray[ni]->pMatElement;
		for(nj=0; nj<nProbeNum; nj++)
		{
			vSum[nj] += vExp[nj];
		}
	}
	
	for(ni=0; ni<nClassNum; ni++)
	{
		if((pClassSizeCopy->pMatElement)[ni] != (pClassSize->pMatElement)[ni])
		{
			printf("Error: Expression_GeneSelection_SubtractMean, class size not match!\n");
			exit(EXIT_FAILURE);
		}
		if((pClassSizeCopy->pMatElement)[ni] > 0)
		{
			dDenom = (double)((pClassSizeCopy->pMatElement)[ni]);
			vSum = vMeans[ni]->pMatElement;
			for(nj=0; nj<nProbeNum; nj++)
			{
				vSum[nj] /= dDenom;
			}
		}
	}

	for(ni=0; ni<nArrayNum; ni++)
	{
		nClustId = (pClassID->pMatElement)[ni]-1;
		vSum = vMeans[nClustId]->pMatElement;
		vExp = vArray[ni]->pMatElement;
		for(nj=0; nj<nProbeNum; nj++)
		{
			vExp[nj] -= vSum[nj];
		}
	}

	/* -------------- */
	/* release memory */
	/* -------------- */
	for(ni=0; ni<nClassNum; ni++)
	{
		DestroyDoubleMatrix(vMeans[ni]);
		vMeans[ni] = NULL;
	}
	free(vMeans);
	DestroyIntMatrix(pClassSizeCopy);

	/* ------ */
	/* return */
	/* ------ */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneSelection_Evaluate()                                    */
/*  Evaluate a expression.                                                 */
/*  return 1 if the expression is True; 0 if the expression is false.      */
/* ----------------------------------------------------------------------- */ 
int Expression_GeneSelection_Evaluate(struct DOUBLEMATRIX *pDraws, 
									  char vLogic[], double vGid[], 
									  int nLogicLen, int nSimple)
{
	/* define */
	int nResult;
	char vStack[MED_LINE_LENGTH];
	char vLogic0[MED_LINE_LENGTH];
	double vGid0[MED_LINE_LENGTH];
	double vGidTemp[MED_LINE_LENGTH];
	int ni,nj,nk,nLen,nId;
	int nSubLen;
	double dTemp1,dTemp2;
	int nTemp1,nTemp2;

	/* init */
	nResult = 0;

	/* judge */
	
	/* priority 1: '(', ')' and replace group id by values */
	if(nSimple == 0)
	{
		nj = 0;
		for(ni=0; ni<nLogicLen; ni++)
		{
			if(vLogic[ni] == '(')
			{
				vStack[nj] = vLogic[ni];
				vGidTemp[nj] = vGid[ni];
				nj++;
			}
			else if(vLogic[ni] == ')')
			{
				nk = nj-1;
				nSubLen = 0;
				while(vStack[nk] != '(')
				{
					nSubLen++;
					nk--;
				}
				vGidTemp[nk] = (double)Expression_GeneSelection_Evaluate(pDraws, (vStack+nk+1), (vGidTemp+nk+1), nSubLen, 1);
				vStack[nk] = 'V';
				nj = nk+1;
			}
			else if(vLogic[ni] == 'G')
			{
				vStack[nj] = 'V';
				nId = (int)(vGid[ni]-1.0);
				vGidTemp[nj] = (pDraws->pMatElement)[nId];
				nj++;
			}
			else
			{
				vStack[nj] = vLogic[ni];
				vGidTemp[nj] = vGid[ni];
				nj++;
			}
		}

		nResult = Expression_GeneSelection_Evaluate(pDraws, vStack, vGidTemp, nj, 1);
	}

	else
	{
		/* priority 2: '<', '>' */
		nj = 0;
		for(ni=0; ni<nLogicLen; ni++)
		{
			if(vLogic[ni] == '<')
			{
				if( (vStack[nj-1] != 'V') || (vLogic[ni+1] != 'V') )
				{
					printf("Error: logic expression error!\n");
					exit(EXIT_FAILURE);
				}

				dTemp1 = vGidTemp[nj-1];
				dTemp2 = vGid[ni+1];
				if( dTemp1 < dTemp2 )
				{
					vStack[nj-1] = 'V';
					vGidTemp[nj-1] = 1.0;
				}
				else
				{
					vStack[nj-1] = 'V';
					vGidTemp[nj-1] = 0.0;
				}
				ni++;
			}
			else if(vLogic[ni] == '>')
			{
				if( (vStack[nj-1] != 'V') || (vLogic[ni+1] != 'V') )
				{
					printf("Error: logic expression error!\n");
					exit(EXIT_FAILURE);
				}

				dTemp1 = vGidTemp[nj-1];
				dTemp2 = vGid[ni+1];
				if( dTemp1 > dTemp2 )
				{
					vStack[nj-1] = 'V';
					vGidTemp[nj-1] = 1.0;
				}
				else
				{
					vStack[nj-1] = 'V';
					vGidTemp[nj-1] = 0.0;
				}
				ni++;
			}
			else
			{
				vStack[nj] = vLogic[ni];
				vGidTemp[nj] = vGid[ni];
				nj++;
			}
		}

		/* priority 3: '&' */
		nLen = nj;
		for(nj=0; nj<nLen; nj++)
		{
			vLogic0[nj] = vStack[nj];
			vGid0[nj] = vGidTemp[nj];
		}
		nj = 0;
		for(ni=0; ni<nLen; ni++)
		{
			if(vLogic0[ni] == '&')
			{
				if( (vStack[nj-1] != 'V') || (vLogic0[ni+1] != 'V') )
				{
					printf("Error: logic expression error!\n");
					exit(EXIT_FAILURE);
				}

				nTemp1 = (int)vGidTemp[nj-1];
				nTemp2 = (int)vGid0[ni+1];
				if( (nTemp1 == 0) || (nTemp2 == 0) )
				{
					vStack[nj-1] = 'V';
					vGidTemp[nj-1] = 0.0;
				}
				else
				{
					vStack[nj-1] = 'V';
					vGidTemp[nj-1] = 1.0;
				}
				ni++;
			}
			else
			{
				vStack[nj] = vLogic0[ni];
				vGidTemp[nj] = vGid0[ni];
				nj++;
			}
		}

		/* priority 4: '|' */
		nLen = nj;
		for(nj=0; nj<nLen; nj++)
		{
			vLogic0[nj] = vStack[nj];
			vGid0[nj] = vGidTemp[nj];
		}
		nj = 0;
		for(ni=0; ni<nLen; ni++)
		{
			if(vLogic0[ni] == '|')
			{
				if( (vStack[nj-1] != 'V') || (vLogic0[ni+1] != 'V') )
				{
					printf("Error: logic expression error!\n");
					exit(EXIT_FAILURE);
				}

				nTemp1 = (int)vGidTemp[nj-1];
				nTemp2 = (int)vGid0[ni+1];
				if( (nTemp1 == 0) && (nTemp2 == 0) )
				{
					vStack[nj-1] = 'V';
					vGidTemp[nj-1] = 0.0;
				}
				else
				{
					vStack[nj-1] = 'V';
					vGidTemp[nj-1] = 1.0;
				}
				ni++;
			}
			else
			{
				vStack[nj] = vLogic0[ni];
				vGidTemp[nj] = vGid0[ni];
				nj++;
			}
		}

		/* result */
		nLen = nj;
		if(nLen == 0)
			nResult = 1;
		else if(nLen == 1)
		{
			if(vStack[0] == 'V')
			{
				nResult = ((int)(vGidTemp[0]) != 0);
			}
			else
			{
				printf("Error: logic expression error!\n");
				exit(EXIT_FAILURE);
			}
		}
		else
		{
			printf("Error: logic expression error!\n");
			exit(EXIT_FAILURE);
		}

	}

	/* return */
	return nResult;
}

/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneSelection_EvaluateVec()                                 */
/*  Evaluate a expression.                                                 */
/*  return 0 if any error happened.                                        */
/*  return 1 if the result is a value vector pointer which is not new and  */
/*  shouldn't be deleted.                                                  */
/*  return 2 if the result is a value vector pointer which is new and      */
/*  need to be deleted later.                                              */
/*  return 3 if the result is a logic vector pointer which is not new and  */
/*  shouldn't be deleted.                                                  */
/*  return 4 if the result is a logic vector pointer which is new and      */
/*  need to be deleted later.                                              */
/*  V - double value, deletable                                            */
/*  G - double value, not deletable                                        */
/*  L - logic value, deletable                                             */
/*  T - logic value, not deletable                                         */
/* ----------------------------------------------------------------------- */ 
int Expression_GeneSelection_EvaluateVec(int nProbeNum,
									  char vLogic[], 
									  struct DOUBLEMATRIX *vVid[], 
									  struct BYTEMATRIX *vLid[],
									  int nLogicLen, int nSimple,
									  struct BYTEMATRIX **pResultVec, 
									  struct DOUBLEMATRIX **pResultVal)
{
	/* define */
	int nResultType;
	struct BYTEMATRIX *pResultL;
	struct DOUBLEMATRIX *pResultV;

	char vStack[MED_LINE_LENGTH];
	char vLogic0[MED_LINE_LENGTH];
	struct DOUBLEMATRIX *vVid0[MED_LINE_LENGTH];
	struct DOUBLEMATRIX *vVidTemp[MED_LINE_LENGTH];
	struct BYTEMATRIX *vLid0[MED_LINE_LENGTH];
	struct BYTEMATRIX *vLidTemp[MED_LINE_LENGTH];
	int ni,nj,nk,nx,nLen;
	int nSubLen;
	int nTemp1,nTemp2;
	int nTempType;

	struct BYTEMATRIX *pLogicTemp;
	double *vV1,*vV2;
	unsigned char *vL1,*vL2,*vL3;

	/* init */
	pResultL = NULL;
	pResultV = NULL;

	/* judge */
	
	/* priority 1: '(', ')' and replace group id by values */
	if(nSimple == 0)
	{
		nj = 0;
		for(ni=0; ni<nLogicLen; ni++)
		{
			if(vLogic[ni] == '(')
			{
				vStack[nj] = vLogic[ni];
				vVidTemp[nj] = vVid[ni];
				vLidTemp[nj] = vLid[ni];
				nj++;
			}
			else if(vLogic[ni] == ')')
			{
				nk = nj-1;
				nSubLen = 0;
				while(vStack[nk] != '(')
				{
					nSubLen++;
					nk--;
				}
				pResultL = NULL;
				pResultV = NULL;
				nTempType = Expression_GeneSelection_EvaluateVec(nProbeNum,
									  (vStack+nk+1), 
									  (vVidTemp+nk+1), 
									  (vLidTemp+nk+1),
									  nSubLen, 1,
									  &pResultL, 
									  &pResultV);

				/* destroy old */
				for(nx=nk; nx<nj; nx++)
				{
					if(vStack[nx] == 'V')
					{
						DestroyDoubleMatrix(vVidTemp[nx]);
						vVidTemp[nx] = NULL;
						if(vLidTemp[nx] != NULL)
						{
							printf("Error: pointer wrong!\n");
							exit(EXIT_FAILURE);
						}
					}
					else if(vStack[nx] == 'G')
					{
						vVidTemp[nx] = NULL;
						if(vLidTemp[nx] != NULL)
						{
							printf("Error: pointer wrong!\n");
							exit(EXIT_FAILURE);
						}
					}
					else if(vStack[nx] == 'L')
					{
						DestroyByteMatrix(vLidTemp[nx]);
						vLidTemp[nx] = NULL;
						if(vVidTemp[nx] != NULL)
						{
							printf("Error: pointer wrong!\n");
							exit(EXIT_FAILURE);
						}
					}
					else if(vStack[nx] == 'T')
					{
						vLidTemp[nx] = NULL;
						if(vVidTemp[nx] != NULL)
						{
							printf("Error: pointer wrong!\n");
							exit(EXIT_FAILURE);
						}
					}
					else
					{
						if((vVidTemp[nx] != NULL) || (vLidTemp[nx] != NULL))
						{
							printf("Error: pointer wrong!\n");
							exit(EXIT_FAILURE);
						}
					}
				}


				/* establish new */
				if( nTempType == 1 )
				{
					vStack[nk] = 'G';
					vVidTemp[nk] = pResultV;
					pResultV = NULL;
				}
				else if( nTempType == 2 )
				{
					vStack[nk] = 'V';
					vVidTemp[nk] = pResultV;
					pResultV = NULL;
				}
				else if( nTempType == 3 )
				{
					vStack[nk] = 'T';
					vLidTemp[nk] = pResultL;
					pResultL = NULL;
				}
				else if( nTempType == 4 )
				{
					vStack[nk] = 'L';
					vLidTemp[nk] = pResultL;
					pResultL = NULL;
				}
				else
				{
					printf("Error: error in logic expression evaluation!\n");
					exit(EXIT_FAILURE);
				}
				nj = nk+1;
			}
			else
			{
				vStack[nj] = vLogic[ni];
				vVidTemp[nj] = vVid[ni];
				vLidTemp[nj] = vLid[ni];
				nj++;
			}
		}

		nLen = nj;
		nResultType = Expression_GeneSelection_EvaluateVec(nProbeNum,
									  vStack, vVidTemp, vLidTemp,
									  nLen, 1,
									  &pResultL, 
									  &pResultV);
		
		/* empty stack */
		for(nx=0; nx<nLen; nx++)
		{
			if(vStack[nx] == 'V')
			{
				DestroyDoubleMatrix(vVidTemp[nx]);
				vVidTemp[nx] = NULL;
				if(vLidTemp[nx] != NULL)
				{
					printf("Error: pointer wrong!\n");
					exit(EXIT_FAILURE);
				}
			}
			else if(vStack[nx] == 'G')
			{
				vVidTemp[nx] = NULL;
				if(vLidTemp[nx] != NULL)
				{
					printf("Error: pointer wrong!\n");
					exit(EXIT_FAILURE);
				}
			}
			else if(vStack[nx] == 'L')
			{
				DestroyByteMatrix(vLidTemp[nx]);
				vLidTemp[nx] = NULL;
				if(vVidTemp[nx] != NULL)
				{
					printf("Error: pointer wrong!\n");
					exit(EXIT_FAILURE);
				}
			}
			else if(vStack[nx] == 'T')
			{
				vLidTemp[nx] = NULL;
				if(vVidTemp[nx] != NULL)
				{
					printf("Error: pointer wrong!\n");
					exit(EXIT_FAILURE);
				}
			}
			else
			{
				if((vVidTemp[nx] != NULL) || (vLidTemp[nx] != NULL))
				{
					printf("Error: pointer wrong!\n");
					exit(EXIT_FAILURE);
				}
			}
		}

		/* return result */
		if( (nResultType == 1) ||  (nResultType == 2) )
		{
			if((pResultL != NULL) || (pResultV == NULL))
			{
				printf("Error: error in logic expression evaluation!\n");
				exit(EXIT_FAILURE);
			}
		}
		else if( (nResultType == 3) ||  (nResultType == 4) )
		{
			if((pResultL == NULL) || (pResultV != NULL))
			{
				printf("Error: error in logic expression evaluation!\n");
				exit(EXIT_FAILURE);
			}
		}
		else
		{
			printf("Error: error in logic expression evaluation!\n");
			exit(EXIT_FAILURE);
		}

		*pResultVec = pResultL;
		*pResultVal = pResultV;
	}

	else
	{
		/* priority 2: '<', '>' */
		nLen = nLogicLen;
		if(nLen == 0)
		{
			*pResultVec = NULL;
			*pResultVal = NULL;
			nResultType = 0;
			return nResultType;
		}
		if(nLen == 1)
		{
			if(vLogic[0] == 'G')
			{
				*pResultVal = vVid[0];
				*pResultVec = NULL;
				nResultType = 1;
				return nResultType;
			}
			else if(vLogic[0] == 'V')
			{
				vLogic[0] = 'G';
				*pResultVal = vVid[0];
				*pResultVec = NULL;
				nResultType = 2;
				return nResultType;
			}
			else if(vLogic[0] == 'T')
			{
				*pResultVal = NULL;
				*pResultVec = vLid[0];
				nResultType = 3;
				return nResultType;
			}
			else if(vLogic[0] == 'L')
			{
				vLogic[0] = 'T';
				*pResultVal = NULL;
				*pResultVec = vLid[0];
				nResultType = 4;
				return nResultType;
			}
			else
			{
				*pResultVec = NULL;
				*pResultVal = NULL;
				nResultType = 0;
				return nResultType;
			}
		}
		
		/* nLen >= 2 */
		for(nj=0; nj<nLen; nj++)
		{
			vLogic0[nj] = vLogic[nj];
			vVid0[nj] = vVid[nj];
			vLid0[nj] = vLid[nj];
			if(vLogic0[nj] == 'V')
				vLogic0[nj] = 'G';
			else if(vLogic0[nj] == 'L')
				vLogic0[nj] = 'T';
		}

		nj = 0;
		for(ni=0; ni<nLen; ni++)
		{
			if(vLogic0[ni] == '<')
			{
				if((vStack[nj-1] == 'G') || (vStack[nj-1] == 'V'))
				{
					nTemp1 = 1;
				}
				else
				{
					nTemp1 = 0;
				}

				if((vLogic0[ni+1] == 'G') || (vLogic0[ni+1] == 'V'))
				{
					nTemp2 = 1;
				}
				else
				{
					nTemp2 = 0;
				}

				if( (nTemp1==0) || (nTemp2==0) )
				{
					printf("Error: logic expression error!\n");
					exit(EXIT_FAILURE);
				}

				pLogicTemp = NULL;
				pLogicTemp = CreateByteMatrix(1, nProbeNum);
				if(pLogicTemp == NULL)
				{
					printf("Error: cannot create space for logic vector!\n");
					exit(EXIT_FAILURE);
				}

				vV1 = vVidTemp[nj-1]->pMatElement;
				vV2 = vVid0[ni+1]->pMatElement;
				vL1 = pLogicTemp->pMatElement;

				for(nk=0; nk<nProbeNum; nk++)
				{
					if( vV1[nk] < vV2[nk] )
					{
						vL1[nk] = 1;
					}
					else
					{
						vL1[nk] = 0;
					}
				}

				if(vStack[nj-1] == 'V')
				{
					DestroyDoubleMatrix(vVidTemp[nj-1]);
				}
				if(vLogic0[ni+1] == 'V')
				{
					DestroyDoubleMatrix(vVid0[ni+1]);
				}
				if(vLidTemp[nj-1] != NULL)
				{
					printf("Error: logic expression error!\n");
					exit(EXIT_FAILURE);
				}
				vVidTemp[nj-1] = NULL;
				vLidTemp[nj-1] = NULL;
				vVid0[ni+1] = NULL;

				vStack[nj-1] = 'L';
				vLidTemp[nj-1] = pLogicTemp;
				ni++;
			}
			else if(vLogic[ni] == '>')
			{
				if((vStack[nj-1] == 'G') || (vStack[nj-1] == 'V'))
				{
					nTemp1 = 1;
				}
				else
				{
					nTemp1 = 0;
				}

				if((vLogic0[ni+1] == 'G') || (vLogic0[ni+1] == 'V'))
				{
					nTemp2 = 1;
				}
				else
				{
					nTemp2 = 0;
				}

				if( (nTemp1==0) || (nTemp2==0) )
				{
					printf("Error: logic expression error!\n");
					exit(EXIT_FAILURE);
				}

				pLogicTemp = NULL;
				pLogicTemp = CreateByteMatrix(1, nProbeNum);
				if(pLogicTemp == NULL)
				{
					printf("Error: cannot create space for logic vector!\n");
					exit(EXIT_FAILURE);
				}

				vV1 = vVidTemp[nj-1]->pMatElement;
				vV2 = vVid0[ni+1]->pMatElement;
				vL1 = pLogicTemp->pMatElement;

				for(nk=0; nk<nProbeNum; nk++)
				{
					if( vV1[nk] > vV2[nk] )
					{
						vL1[nk] = 1;
					}
					else
					{
						vL1[nk] = 0;
					}
				}



				if(vStack[nj-1] == 'V')
				{
					DestroyDoubleMatrix(vVidTemp[nj-1]);
				}
				if(vLogic0[ni+1] == 'V')
				{
					DestroyDoubleMatrix(vVid0[ni+1]);
				}
				if(vLidTemp[nj-1] != NULL)
				{
					printf("Error: logic expression error!\n");
					exit(EXIT_FAILURE);
				}
				vVidTemp[nj-1] = NULL;
				vLidTemp[nj-1] = NULL;
				vVid0[ni+1] = NULL;

				vStack[nj-1] = 'L';
				vLidTemp[nj-1] = pLogicTemp;
				ni++;
			}
			else
			{
				vStack[nj] = vLogic0[ni];
				vVidTemp[nj] = vVid0[ni];
				vLidTemp[nj] = vLid0[ni];
				nj++;
			}
		}

		/* priority 3: '&' */
		nLen = nj;
		if(nLen == 0)
		{
			*pResultVec = NULL;
			*pResultVal = NULL;
			nResultType = 0;
			return nResultType;
		}
		if(nLen == 1)
		{
			if(vStack[0] == 'G')
			{
				*pResultVal = vVidTemp[0];
				*pResultVec = NULL;
				nResultType = 1;
				return nResultType;
			}
			else if(vStack[0] == 'V')
			{
				vStack[0] = 'G';
				*pResultVal = vVidTemp[0];
				*pResultVec = NULL;
				nResultType = 2;
				return nResultType;
			}
			else if(vStack[0] == 'T')
			{
				*pResultVal = NULL;
				*pResultVec = vLidTemp[0];
				nResultType = 3;
				return nResultType;
			}
			else if(vStack[0] == 'L')
			{
				vStack[0] = 'T';
				*pResultVal = NULL;
				*pResultVec = vLidTemp[0];
				nResultType = 4;
				return nResultType;
			}
			else
			{
				*pResultVec = NULL;
				*pResultVal = NULL;
				nResultType = 0;
				return nResultType;
			}
		}
		
		
		for(nj=0; nj<nLen; nj++)
		{
			vLogic0[nj] = vStack[nj];
			vVid0[nj] = vVidTemp[nj];
			vLid0[nj] = vLidTemp[nj];
		}
		nj = 0;
		for(ni=0; ni<nLen; ni++)
		{
			if(vLogic0[ni] == '&')
			{
				if((vStack[nj-1] == 'L') || (vStack[nj-1] == 'T'))
				{
					nTemp1 = 1;
				}
				else
				{
					nTemp1 = 0;
				}

				if((vLogic0[ni+1] == 'L') || (vLogic0[ni+1] == 'T'))
				{
					nTemp2 = 1;
				}
				else
				{
					nTemp2 = 0;
				}

				if( (nTemp1==0) || (nTemp2==0) )
				{
					printf("Error: logic expression error!\n");
					exit(EXIT_FAILURE);
				}

				pLogicTemp = NULL;
				pLogicTemp = CreateByteMatrix(1, nProbeNum);
				if(pLogicTemp == NULL)
				{
					printf("Error: cannot create space for logic vector!\n");
					exit(EXIT_FAILURE);
				}

				vL1 = vLidTemp[nj-1]->pMatElement;
				vL2 = vLid0[ni+1]->pMatElement;
				vL3 = pLogicTemp->pMatElement;
				for(nk=0; nk<nProbeNum; nk++)
				{
					if( (vL1[nk] == 0) || (vL2[nk] == 0) )
					{
						vL3[nk] = 0;
					}
					else
					{
						vL3[nk] = 1;
					}
				}

				if(vStack[nj-1] == 'L')
				{
					DestroyByteMatrix(vLidTemp[nj-1]);
				}
				if(vLogic0[ni+1] == 'L')
				{
					DestroyByteMatrix(vLid0[ni+1]);
				}
				if(vVidTemp[nj-1] != NULL)
				{
					printf("Error: logic expression error!\n");
					exit(EXIT_FAILURE);
				}
				vVidTemp[nj-1] = NULL;
				vLidTemp[nj-1] = NULL;
				vLid0[ni+1] = NULL;

				vStack[nj-1] = 'L';
				vLidTemp[nj-1] = pLogicTemp;
				ni++;
			}
			else
			{
				vStack[nj] = vLogic0[ni];
				vVidTemp[nj] = vVid0[ni];
				vLidTemp[nj] = vLid0[ni];
				nj++;
			}
		}

		/* priority 4: '|' */
		nLen = nj;
		if(nLen == 0)
		{
			*pResultVec = NULL;
			*pResultVal = NULL;
			nResultType = 0;
			return nResultType;
		}
		if(nLen == 1)
		{
			if(vStack[0] == 'G')
			{
				*pResultVal = vVidTemp[0];
				*pResultVec = NULL;
				nResultType = 1;
				return nResultType;
			}
			else if(vStack[0] == 'V')
			{
				vStack[0] = 'G';
				*pResultVal = vVidTemp[0];
				*pResultVec = NULL;
				nResultType = 2;
				return nResultType;
			}
			else if(vStack[0] == 'T')
			{
				*pResultVal = NULL;
				*pResultVec = vLidTemp[0];
				nResultType = 3;
				return nResultType;
			}
			else if(vStack[0] == 'L')
			{
				vStack[0] = 'T';
				*pResultVal = NULL;
				*pResultVec = vLidTemp[0];
				nResultType = 4;
				return nResultType;
			}
			else
			{
				*pResultVec = NULL;
				*pResultVal = NULL;
				nResultType = 0;
				return nResultType;
			}
		}
		
		
		for(nj=0; nj<nLen; nj++)
		{
			vLogic0[nj] = vStack[nj];
			vVid0[nj] = vVidTemp[nj];
			vLid0[nj] = vLidTemp[nj];
		}
		nj = 0;
		for(ni=0; ni<nLen; ni++)
		{
			if(vLogic0[ni] == '|')
			{
				if((vStack[nj-1] == 'L') || (vStack[nj-1] == 'T'))
				{
					nTemp1 = 1;
				}
				else
				{
					nTemp1 = 0;
				}

				if((vLogic0[ni+1] == 'L') || (vLogic0[ni+1] == 'T'))
				{
					nTemp2 = 1;
				}
				else
				{
					nTemp2 = 0;
				}

				if( (nTemp1==0) || (nTemp2==0) )
				{
					printf("Error: logic expression error!\n");
					exit(EXIT_FAILURE);
				}

				pLogicTemp = NULL;
				pLogicTemp = CreateByteMatrix(1, nProbeNum);
				if(pLogicTemp == NULL)
				{
					printf("Error: cannot create space for logic vector!\n");
					exit(EXIT_FAILURE);
				}

				vL1 = vLidTemp[nj-1]->pMatElement;
				vL2 = vLid0[ni+1]->pMatElement;
				vL3 = pLogicTemp->pMatElement;
				for(nk=0; nk<nProbeNum; nk++)
				{
					if( (vL1[nk] == 0 ) && (vL2[nk] == 0) )
					{
						vL3[nk] = 0;
					}
					else
					{
						vL3[nk] = 1;
					}
				}

				if(vStack[nj-1] == 'L')
				{
					DestroyByteMatrix(vLidTemp[nj-1]);
				}
				if(vLogic0[ni+1] == 'L')
				{
					DestroyByteMatrix(vLid0[ni+1]);
				}
				if(vVidTemp[nj-1] != NULL)
				{
					printf("Error: logic expression error!\n");
					exit(EXIT_FAILURE);
				}
				vVidTemp[nj-1] = NULL;
				vLidTemp[nj-1] = NULL;
				vLid0[ni+1] = NULL;

				vStack[nj-1] = 'L';
				vLidTemp[nj-1] = pLogicTemp;
				ni++;
			}
			else
			{
				vStack[nj] = vLogic0[ni];
				vVidTemp[nj] = vVid0[ni];
				vLidTemp[nj] = vLid0[ni];
				nj++;
			}
		}

		/* result */
		nLen = nj;
		if(nLen == 0)
		{
			*pResultVec = NULL;
			*pResultVal = NULL;
			nResultType = 0;
		}
		else if(nLen == 1)
		{
			if(vStack[0] == 'G')
			{
				*pResultVal = vVidTemp[0];
				*pResultVec = NULL;
				nResultType = 1;
			}
			else if(vStack[0] == 'V')
			{
				vStack[0] = 'G';
				*pResultVal = vVidTemp[0];
				*pResultVec = NULL;
				nResultType = 2;
			}
			else if(vStack[0] == 'T')
			{
				*pResultVal = NULL;
				*pResultVec = vLidTemp[0];
				nResultType = 3;
			}
			else if(vStack[0] == 'L')
			{
				vStack[0] = 'T';
				*pResultVal = NULL;
				*pResultVec = vLidTemp[0];
				nResultType = 4;
			}
			else
			{
				*pResultVec = NULL;
				*pResultVal = NULL;
				nResultType = 0;
			}
		}

	}

	/* return */
	return nResultType;
}


/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneSelection_ClassPerm()                                   */
/*  permute class labels.                                                  */
/*  return permuted class labels.                                          */
/* ----------------------------------------------------------------------- */ 
struct INTMATRIX *Expression_GeneSelection_ClassPerm(int nClassNum,
				struct INTMATRIX *pClassID, struct INTMATRIX *pClassSize, 
				int nPermgroupNum, struct INTMATRIX *pPermgroupSize, 
				struct INTMATRIX *vPermgroupMap[])
{
	/* define */
	struct INTMATRIX *pPermClassID;
	struct DOUBLEMATRIX *pRandMat;
	struct DOUBLEMATRIX *pSortRand;
	struct LONGMATRIX *pSortID;
	int ni,nj,nk,nx;
	int nLen,nId,nId0;
	int nOK;
	int nSum;

	/* init */
	if(pClassID == NULL)
		return NULL;
	pPermClassID = NULL;
	pPermClassID = CreateIntMatrix(pClassID->nHeight, pClassID->nWidth);

	/* permutation */
	for(ni=0; ni<nPermgroupNum; ni++)
	{
		if(pPermgroupSize->pMatElement[ni] != vPermgroupMap[ni]->nWidth)
		{
			printf("Error: permutation error!\n");
			exit(EXIT_FAILURE);
		}

		nLen = 0;
		for(nj=0; nj<vPermgroupMap[ni]->nWidth; nj++)
		{
			nId = vPermgroupMap[ni]->pMatElement[nj] - 1;
			nLen += pClassSize->pMatElement[nId];
		}

		if(nLen <= 0)
			continue;

		/* rand */
		pRandMat = NULL;
		pRandMat = DMRANDU(1, nLen);
		if(pRandMat == NULL)
		{
			printf("Error: permutation error!\n");
			exit(EXIT_FAILURE);
		}
		pSortRand = NULL;
		pSortID = NULL;
		if( DMSORTMERGEA_0(pRandMat, &pSortRand, &pSortID) == PROC_FAILURE )
		{
			printf("Error: permutation error!\n");
			exit(EXIT_FAILURE);
		}

		nx = 0;
		for(nj=0; nj<pClassID->nWidth; nj++)
		{
			/* match class label */
			nId0 = pClassID->pMatElement[nj];
			nOK = 0;
			for(nk=0; nk<vPermgroupMap[ni]->nWidth; nk++)
			{
				nId = vPermgroupMap[ni]->pMatElement[nk];
				if(nId0 == nId)
				{
					nOK = 1;
					break;
				}
			}

			/* if not match, continue */
			if(nOK == 0)
				continue;

			/* if match, get new label */
			nId = pSortID->pMatElement[nx];
			nSum = 0;
			for(nk=0; nk<vPermgroupMap[ni]->nWidth; nk++)
			{
				nId0 = vPermgroupMap[ni]->pMatElement[nk]-1;
				nSum += pClassSize->pMatElement[nId0];
				if(nId < nSum)
				{
					nId = nId0+1;
					break;
				}
			}


			pPermClassID->pMatElement[nj] = nId;

			/* get next */
			nx++;
		}

		/* release memory */
		DestroyDoubleMatrix(pRandMat);
		DestroyDoubleMatrix(pSortRand);
		DestroyLongMatrix(pSortID);
	}


	/* IMSAVE(pPermClassID, "debugperm.txt"); */

	/* return */
	return pPermClassID;
}

/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneSelection_Output()                                      */
/*  Output results.                                                        */
/* ----------------------------------------------------------------------- */ 
int Expression_GeneSelection_Output(int nProbeNum, int nOutputNum, 
				struct DOUBLEMATRIX *pScore, struct LONGMATRIX *pSortID, 
				struct DOUBLEMATRIX *pSortScore, struct DOUBLEMATRIX *pFDR,
				char strProbeFile[], char strGeneInfoFile[], char strOutFile[])
{
	/* define */
	struct tagString **vProbeName;
	int ni,nj;
	char strOutPath[LINE_LENGTH];
	char strLine[LINE_LENGTH];
	char strLine2[LONG_LINE_LENGTH];
	char strProbe2[LINE_LENGTH];
	char strProbe[LINE_LENGTH];
	FILE *fpIn;
	FILE *fpOut;
	FILE *fpInfo;
	int nIncludeInfo;
	long *vSID;
	double *vSS;
	double *vSF;
	char *chSep;

	/* init */
	vProbeName = NULL;
	vProbeName = (struct tagString **)calloc(nProbeNum, sizeof(struct tagString *));

	/* load probe name */
	fpIn = NULL;
	fpIn = fopen(strProbeFile, "rt");
	if(fpIn == NULL)
	{
		printf("Error: Expression_GeneSelection_Output, cannot open probeset name file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(ni >= nProbeNum)
		{
			printf("Error: Expression_GeneSelection_Output, probeset number not match!\n");
			exit(EXIT_FAILURE);
		}

		StringAddTail((vProbeName+ni), strLine);
		ni++;
	}

	fclose(fpIn);
	if(ni != nProbeNum)
	{
		printf("Error: Expression_GeneSelection_Output, probeset number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* output unsorted */
	sprintf(strOutPath, "%s.ori", strOutFile);
	fpOut = NULL;
	fpOut = fopen(strOutPath, "wt");
	if(fpOut == NULL)
	{
		printf("Error: Expression_GeneSelection_Output, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	vSS = pScore->pMatElement;
	fprintf(fpOut, "probeset_id\tscore\n");
	for(ni=0; ni<nProbeNum; ni++)
	{
		if(vProbeName[ni] != NULL)
		{
			fprintf(fpOut, "%s\t", vProbeName[ni]->m_pString);
		}
		else
		{
			fprintf(fpOut, "NA\t");
		}
		
		fprintf(fpOut, "% 9.7e\n", vSS[ni]);
	}

	fclose(fpOut);

	/* output selected probes */
	sprintf(strOutPath, "%s.txt", strOutFile);
	fpOut = NULL;
	fpOut = fopen(strOutPath, "wt");
	if(fpOut == NULL)
	{
		printf("Error: Expression_GeneSelection_Output, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	
	if(strcmp(strGeneInfoFile, "NULL") == 0)
	{
		nIncludeInfo = 0;
	}
	else
	{
		nIncludeInfo = 1;
	}

	vSID = pSortID->pMatElement;
	vSS = pSortScore->pMatElement;
	if(pFDR != NULL)
	{
		vSF = pFDR->pMatElement;
	}

	fprintf(fpOut, "probeset_id\tline_id\trank\tscore\tFDR\tInfo\n");
	for(ni=0; ni<nOutputNum; ni++)
	{
		if(ni >= nProbeNum)
			break;
		nj = pSortID->pMatElement[ni];

		if(vProbeName[nj] != NULL)
		{
			fprintf(fpOut, "%s\t", vProbeName[nj]->m_pString);
		}
		else
		{
			fprintf(fpOut, "NA\t");
		}
		
		fprintf(fpOut, "%d\t%d\t% 9.7e\t", (nj+1), (ni+1), vSS[ni]);
		
		if(pFDR != NULL)
		{
			fprintf(fpOut, "% 9.7e\t", vSF[ni]);
		}
		else
		{
			fprintf(fpOut, "NA\t");
		}

		if(nIncludeInfo == 1)
		{
			fpInfo = NULL;
			fpInfo = fopen(strGeneInfoFile, "rt");
			if(fpInfo != NULL)
			{
				strcpy(strProbe, vProbeName[nj]->m_pString);
				fgets(strLine2, LONG_LINE_LENGTH, fpInfo);
				while(fgets(strLine2, LONG_LINE_LENGTH, fpInfo) != NULL)
				{
					sscanf(strLine2, "%s ", strProbe2);
					if(strcmp(strProbe, strProbe2) == 0)
					{
						StrTrimLeft(strLine2);
						chSep = strLine2+strlen(strProbe2);
						StrTrimLeft(chSep);
						StrTrimRight(chSep);
						fprintf(fpOut, "%s", chSep);
						break;
					}
				}
				fclose(fpInfo);
			}
		}
		
		fprintf(fpOut, "\n");
	}

	fclose(fpOut);

	/* release memory */
	for(ni=0; ni<nProbeNum; ni++)
	{
		DeleteString(vProbeName[ni]);
		vProbeName[ni] = NULL;
	}
	free(vProbeName);

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneSelection_Output_1()                                    */
/*  Output results.                                                        */
/* ----------------------------------------------------------------------- */ 
int Expression_GeneSelection_Output_1(int nProbeNum, int nOutputNum, 
				struct DOUBLEMATRIX *pScore, struct LONGMATRIX *pSortID, 
				struct DOUBLEMATRIX *pSortScore, struct DOUBLEMATRIX *pFDR,
				char strProbeFile[], char strGeneInfoFile[], char strOutFile[])
{
	/* define */
	struct tagString **vProbeName;
	int ni,nj;
	char strOutPath[LINE_LENGTH];
	char strLine[LINE_LENGTH];
	char strLine2[LONG_LINE_LENGTH];
	char strProbe2[LINE_LENGTH];
	char strProbe[LINE_LENGTH];
	char strGenShort[LINE_LENGTH];
	FILE *fpIn;
	FILE *fpOut;
	FILE *fpInfo;
	int nIncludeInfo;
	long *vSID;
	double *vSS;
	double *vSF;
	char *chSep,*chGen;

	/* init */
	vProbeName = NULL;
	vProbeName = (struct tagString **)calloc(nProbeNum, sizeof(struct tagString *));

	/* load probe name */
	fpIn = NULL;
	fpIn = fopen(strProbeFile, "rt");
	if(fpIn == NULL)
	{
		printf("Error: Expression_GeneSelection_Output, cannot open probeset name file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(ni >= nProbeNum)
		{
			printf("Error: Expression_GeneSelection_Output, probeset number not match!\n");
			exit(EXIT_FAILURE);
		}

		StringAddTail((vProbeName+ni), strLine);
		ni++;
	}

	fclose(fpIn);
	if(ni != nProbeNum)
	{
		printf("Error: Expression_GeneSelection_Output, probeset number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* output unsorted */
	sprintf(strOutPath, "%s.ori", strOutFile);
	fpOut = NULL;
	fpOut = fopen(strOutPath, "wt");
	if(fpOut == NULL)
	{
		printf("Error: Expression_GeneSelection_Output, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	vSS = pScore->pMatElement;
	fprintf(fpOut, "probeset_id\tscore\n");
	for(ni=0; ni<nProbeNum; ni++)
	{
		if(vProbeName[ni] != NULL)
		{
			fprintf(fpOut, "%s\t", vProbeName[ni]->m_pString);
		}
		else
		{
			fprintf(fpOut, "NA\t");
		}
		
		fprintf(fpOut, "% 9.7e\n", vSS[ni]);
	}

	fclose(fpOut);

	/* output */
	sprintf(strOutPath, "%s.txt", strOutFile);
	fpOut = NULL;
	fpOut = fopen(strOutPath, "wt");
	if(fpOut == NULL)
	{
		printf("Error: Expression_GeneSelection_Output, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	
	if(strcmp(strGeneInfoFile, "NULL") == 0)
	{
		nIncludeInfo = 0;
	}
	else
	{
		nIncludeInfo = 1;
	}

	vSID = pSortID->pMatElement;
	vSS = pSortScore->pMatElement;
	if(pFDR != NULL)
	{
		vSF = pFDR->pMatElement;
	}

	fprintf(fpOut, "probeset_id\trank\tscore\tFDR\tGene\tUG_Gene\tInfo\n");
	for(ni=0; ni<nOutputNum; ni++)
	{
		if(ni >= nProbeNum)
			break;
		nj = pSortID->pMatElement[ni];

		if(vProbeName[nj] != NULL)
		{
			fprintf(fpOut, "%s\t", vProbeName[nj]->m_pString);
		}
		else
		{
			fprintf(fpOut, "NA\t");
		}
		
		fprintf(fpOut, "%d\t% 9.7e\t", (ni+1), vSS[ni]);
		
		if(pFDR != NULL)
		{
			fprintf(fpOut, "% 9.7e\t", vSF[ni]);
		}
		else
		{
			fprintf(fpOut, "NA\t");
		}

		if(nIncludeInfo == 1)
		{
			fpInfo = NULL;
			fpInfo = fopen(strGeneInfoFile, "rt");
			if(fpInfo != NULL)
			{
				strcpy(strProbe, vProbeName[nj]->m_pString);
				fgets(strLine2, LONG_LINE_LENGTH, fpInfo);
				while(fgets(strLine2, LONG_LINE_LENGTH, fpInfo) != NULL)
				{
					sscanf(strLine2, "%s ", strProbe2);
					if(strcmp(strProbe, strProbe2) == 0)
					{
						StrTrimLeft(strLine2);
						chSep = strLine2+strlen(strProbe2);
						StrTrimLeft(chSep);
						StrTrimRight(chSep);
						chGen = strstr(chSep, "/GEN=");
						if(chGen != NULL)
						{
							sscanf(chGen, "%s ", strGenShort);
							fprintf(fpOut, "%s\t", (strGenShort+5));
						}
						else
						{
							fprintf(fpOut, "NA\t");
						}
						chGen = strstr(chSep, "/UG_GENE=");
						if(chGen != NULL)
						{
							sscanf(chGen, "%s ", strGenShort);
							fprintf(fpOut, "%s\t", (strGenShort+9));
						}
						else
						{
							fprintf(fpOut, "NA\t");
						}
						fprintf(fpOut, "%s", chSep);
						break;
					}
				}
				fclose(fpInfo);
			}
		}
		
		fprintf(fpOut, "\n");
	}

	fclose(fpOut);

	/* release memory */
	for(ni=0; ni<nProbeNum; ni++)
	{
		DeleteString(vProbeName[ni]);
		vProbeName[ni] = NULL;
	}
	free(vProbeName);

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneSelection_Output_WithRandomControl()                    */
/*  Output results, a group of random control will be selected.            */
/* ----------------------------------------------------------------------- */ 
int Expression_GeneSelection_Output_WithRandomControl(int nProbeNum, int nOutputNum, 
				struct DOUBLEMATRIX *pScore, struct LONGMATRIX *pSortID, 
				struct DOUBLEMATRIX *pSortScore, struct DOUBLEMATRIX *pFDR,
				char strProbeFile[], char strGeneInfoFile[], char strOutFile[])
{
	/* define */
	struct tagString **vProbeName;
	int ni,nj,nk;
	char strOutPath[LINE_LENGTH];
	char strLine[LINE_LENGTH];
	char strLine2[LONG_LINE_LENGTH];
	char strProbe2[LINE_LENGTH];
	char strProbe[LINE_LENGTH];
	FILE *fpIn;
	FILE *fpOut;
	FILE *fpInfo;
	int nIncludeInfo;
	long *vSID;
	double *vSS;
	double *vSF;
	char *chSep;

	struct DOUBLEMATRIX *pRandMat;
	struct DOUBLEMATRIX *pRandSort;
	struct LONGMATRIX *pRandSortId;


	/* init */
	vProbeName = NULL;
	vProbeName = (struct tagString **)calloc(nProbeNum, sizeof(struct tagString *));

	/* load probe name */
	fpIn = NULL;
	fpIn = fopen(strProbeFile, "rt");
	if(fpIn == NULL)
	{
		printf("Error: Expression_GeneSelection_Output, cannot open probeset name file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(ni >= nProbeNum)
		{
			printf("Error: Expression_GeneSelection_Output, probeset number not match!\n");
			exit(EXIT_FAILURE);
		}

		StringAddTail((vProbeName+ni), strLine);
		ni++;
	}

	fclose(fpIn);
	if(ni != nProbeNum)
	{
		printf("Error: Expression_GeneSelection_Output, probeset number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* output unsorted */
	sprintf(strOutPath, "%s.ori", strOutFile);
	fpOut = NULL;
	fpOut = fopen(strOutPath, "wt");
	if(fpOut == NULL)
	{
		printf("Error: Expression_GeneSelection_Output, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	vSS = pScore->pMatElement;
	fprintf(fpOut, "probeset_id\tscore\n");
	for(ni=0; ni<nProbeNum; ni++)
	{
		if(vProbeName[ni] != NULL)
		{
			fprintf(fpOut, "%s\t", vProbeName[ni]->m_pString);
		}
		else
		{
			fprintf(fpOut, "NA\t");
		}
		
		fprintf(fpOut, "% 9.7e\n", vSS[ni]);
	}

	fclose(fpOut);

	/* output selected probes */
	sprintf(strOutPath, "%s.txt", strOutFile);
	fpOut = NULL;
	fpOut = fopen(strOutPath, "wt");
	if(fpOut == NULL)
	{
		printf("Error: Expression_GeneSelection_Output, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	
	if(strcmp(strGeneInfoFile, "NULL") == 0)
	{
		nIncludeInfo = 0;
	}
	else
	{
		nIncludeInfo = 1;
	}

	vSID = pSortID->pMatElement;
	vSS = pSortScore->pMatElement;
	if(pFDR != NULL)
	{
		vSF = pFDR->pMatElement;
	}

	fprintf(fpOut, "probeset_id\tline_id\trank\tscore\tFDR\tInfo\n");
	for(ni=0; ni<nOutputNum; ni++)
	{
		if(ni >= nProbeNum)
			break;
		nj = pSortID->pMatElement[ni];

		if(vProbeName[nj] != NULL)
		{
			fprintf(fpOut, "%s\t", vProbeName[nj]->m_pString);
		}
		else
		{
			fprintf(fpOut, "NA\t");
		}
		
		fprintf(fpOut, "%d\t%d\t% 9.7e\t", (nj+1), (ni+1), vSS[ni]);
		
		if(pFDR != NULL)
		{
			fprintf(fpOut, "% 9.7e\t", vSF[ni]);
		}
		else
		{
			fprintf(fpOut, "NA\t");
		}

		if(nIncludeInfo == 1)
		{
			fpInfo = NULL;
			fpInfo = fopen(strGeneInfoFile, "rt");
			if(fpInfo != NULL)
			{
				strcpy(strProbe, vProbeName[nj]->m_pString);
				fgets(strLine2, LONG_LINE_LENGTH, fpInfo);
				while(fgets(strLine2, LONG_LINE_LENGTH, fpInfo) != NULL)
				{
					sscanf(strLine2, "%s ", strProbe2);
					if(strcmp(strProbe, strProbe2) == 0)
					{
						StrTrimLeft(strLine2);
						chSep = strLine2+strlen(strProbe2);
						StrTrimLeft(chSep);
						StrTrimRight(chSep);
						fprintf(fpOut, "%s", chSep);
						break;
					}
				}
				fclose(fpInfo);
			}
		}
		
		fprintf(fpOut, "\n");
	}

	fclose(fpOut);

	/* output random control probes */
	pRandMat = NULL;
	pRandSort = NULL;
	pRandSortId = NULL;
	pRandMat = DMRANDU(1, nProbeNum);
	DMSORTMERGEA_0(pRandMat, &pRandSort, &pRandSortId);

	sprintf(strOutPath, "%s_ctr.txt", strOutFile);
	fpOut = NULL;
	fpOut = fopen(strOutPath, "wt");
	if(fpOut == NULL)
	{
		printf("Error: Expression_GeneSelection_Output, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	
	if(strcmp(strGeneInfoFile, "NULL") == 0)
	{
		nIncludeInfo = 0;
	}
	else
	{
		nIncludeInfo = 0;
	}

	vSID = pRandSortId->pMatElement;
	vSS = pSortScore->pMatElement;
	if(pFDR != NULL)
	{
		vSF = pFDR->pMatElement;
	}

	fprintf(fpOut, "probeset_id\tline_id\trank\tscore\tFDR\tInfo\n");
	for(ni=0; ni<nOutputNum; ni++)
	{
		if(ni >= nProbeNum)
			break;
		nk = pRandSortId->pMatElement[ni];
		nj = pSortID->pMatElement[nk];

		if(vProbeName[nj] != NULL)
		{
			fprintf(fpOut, "%s\t", vProbeName[nj]->m_pString);
		}
		else
		{
			fprintf(fpOut, "NA\t");
		}
		
		fprintf(fpOut, "%d\t%d\t% 9.7e\t", (nj+1), (nk+1), vSS[nk]);
		
		if(pFDR != NULL)
		{
			fprintf(fpOut, "% 9.7e\t", vSF[nk]);
		}
		else
		{
			fprintf(fpOut, "NA\t");
		}

		if(nIncludeInfo == 1)
		{
			fpInfo = NULL;
			fpInfo = fopen(strGeneInfoFile, "rt");
			if(fpInfo != NULL)
			{
				strcpy(strProbe, vProbeName[nj]->m_pString);
				fgets(strLine2, LONG_LINE_LENGTH, fpInfo);
				while(fgets(strLine2, LONG_LINE_LENGTH, fpInfo) != NULL)
				{
					sscanf(strLine2, "%s ", strProbe2);
					if(strcmp(strProbe, strProbe2) == 0)
					{
						StrTrimLeft(strLine2);
						chSep = strLine2+strlen(strProbe2);
						StrTrimLeft(chSep);
						StrTrimRight(chSep);
						fprintf(fpOut, "%s", chSep);
						break;
					}
				}
				fclose(fpInfo);
			}
		}
		
		fprintf(fpOut, "\n");
	}

	fclose(fpOut);


	/* release memory */
	DestroyDoubleMatrix(pRandMat);
	DestroyDoubleMatrix(pRandSort);
	DestroyLongMatrix(pRandSortId);

	for(ni=0; ni<nProbeNum; ni++)
	{
		DeleteString(vProbeName[ni]);
		vProbeName[ni] = NULL;
	}
	free(vProbeName);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneSelection_Output_WithRandomControl_Fast()               */
/*  Fast version of Expression_GeneSelection_Output_WithRandomControl.     */
/*  Requires annotation ordering to be the same as probeset odering.       */
/*  Output results, a group of random control will be selected.            */
/* ----------------------------------------------------------------------- */ 
int Expression_GeneSelection_Output_WithRandomControl_Fast(int nProbeNum, int nOutputNum, 
				struct DOUBLEMATRIX *pScore, struct LONGMATRIX *pSortID, 
				struct DOUBLEMATRIX *pSortScore, struct DOUBLEMATRIX *pFDR,
				char strProbeFile[], char strGeneInfoFile[], char strOutFile[])
{
	/* define */
	struct tagString **vProbeName;
	struct tagString **vAnnotation;
	int ni,nj,nk;
	char strOutPath[LINE_LENGTH];
	char strLine[LINE_LENGTH];
	char strLine2[LONG_LINE_LENGTH];
	char strProbe[LINE_LENGTH];
	FILE *fpIn;
	FILE *fpOut;
	int nIncludeInfo;
	long *vSID;
	double *vSS;
	double *vSF;
	char *chSep;

	struct DOUBLEMATRIX *pRandMat;
	struct DOUBLEMATRIX *pRandSort;
	struct LONGMATRIX *pRandSortId;


	/* init */
	vProbeName = NULL;
	vProbeName = (struct tagString **)calloc(nProbeNum, sizeof(struct tagString *));

	/* load probe name */
	fpIn = NULL;
	fpIn = fopen(strProbeFile, "rt");
	if(fpIn == NULL)
	{
		printf("Error: Expression_GeneSelection_Output, cannot open probeset name file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(ni >= nProbeNum)
		{
			printf("Error: Expression_GeneSelection_Output, probeset number not match!\n");
			exit(EXIT_FAILURE);
		}

		StringAddTail((vProbeName+ni), strLine);
		ni++;
	}

	fclose(fpIn);
	if(ni != nProbeNum)
	{
		printf("Error: Expression_GeneSelection_Output, probeset number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* output unsorted */
	sprintf(strOutPath, "%s.ori", strOutFile);
	fpOut = NULL;
	fpOut = fopen(strOutPath, "wt");
	if(fpOut == NULL)
	{
		printf("Error: Expression_GeneSelection_Output, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	vSS = pScore->pMatElement;
	fprintf(fpOut, "probeset_id\tscore\n");
	for(ni=0; ni<nProbeNum; ni++)
	{
		if(vProbeName[ni] != NULL)
		{
			fprintf(fpOut, "%s\t", vProbeName[ni]->m_pString);
		}
		else
		{
			fprintf(fpOut, "NA\t");
		}
		
		fprintf(fpOut, "% 9.7e\n", vSS[ni]);
	}

	fclose(fpOut);

	if(strcmp(strGeneInfoFile, "NULL") == 0)
	{
		nIncludeInfo = 0;
	}
	else
	{
		nIncludeInfo = 1;
	}

	if(nIncludeInfo == 1)
	{
		/* load annotation */
		fpIn = NULL;
		fpIn = fopen(strGeneInfoFile, "rt");
		if(fpIn == NULL)
		{
			printf("Error: Expression_GeneSelection_OutputWithRandomControl_Fast, cannot open file to load annotation!\n");
			exit(EXIT_FAILURE);
		}

		vAnnotation = NULL;
		vAnnotation = (struct tagString **)calloc(nProbeNum, sizeof(struct tagString *));
		if(vAnnotation == NULL)
		{
			printf("Error: Expression_GeneSelection_OutputWithRandomControl_Fast, cannot allocate enough memory for loading annotation!\n");
			exit(EXIT_FAILURE);
		}

		fgets(strLine2, LONG_LINE_LENGTH, fpIn);
		ni = 0;
		while(fgets(strLine2, LONG_LINE_LENGTH, fpIn) != NULL)
		{
			StrTrimLeft(strLine2);
			StrTrimRight(strLine2);
			if(strLine2[0] == '\0')
				continue;

			if(ni >= nProbeNum)
			{
				printf("Error: Expression_GeneSelection_OutputWithRandomControl_Fast, probeset number not match!\n");
				exit(EXIT_FAILURE);
			}

			sscanf(strLine2, "%s", strProbe);
			if(strcmp(strProbe, vProbeName[ni]->m_pString) != 0)
			{
				printf("Error: Expression_GeneSelection_OutputWithRandomControl_Fast, annotation order not match original score order!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine2, '\t');
			chSep++;
			StringAddTail((vAnnotation+ni), chSep);
			/* StringAddTail((vAnnotation+ni), strLine2); */
			ni++;
		}

		fclose(fpIn);

		if(ni != nProbeNum)
		{
			printf("Error: Expression_GeneSelection_OutputWithRandomControl_Fast, probeset number not match!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* output selected probes */
	sprintf(strOutPath, "%s.txt", strOutFile);
	fpOut = NULL;
	fpOut = fopen(strOutPath, "wt");
	if(fpOut == NULL)
	{
		printf("Error: Expression_GeneSelection_Output, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	

	vSID = pSortID->pMatElement;
	vSS = pSortScore->pMatElement;
	if(pFDR != NULL)
	{
		vSF = pFDR->pMatElement;
	}

	fprintf(fpOut, "probeset_id\tline_id\trank\tscore\tFDR\tInfo\n");
	for(ni=0; ni<nOutputNum; ni++)
	{
		if(ni >= nProbeNum)
			break;
		nj = pSortID->pMatElement[ni];

		if(vProbeName[nj] != NULL)
		{
			fprintf(fpOut, "%s\t", vProbeName[nj]->m_pString);
		}
		else
		{
			fprintf(fpOut, "NA\t");
		}
		
		fprintf(fpOut, "%d\t%d\t% 9.7e\t", (nj+1), (ni+1), vSS[ni]);
		
		if(pFDR != NULL)
		{
			fprintf(fpOut, "% 9.7e\t", vSF[ni]);
		}
		else
		{
			fprintf(fpOut, "NA\t");
		}

		if(nIncludeInfo == 1)
		{
			fprintf(fpOut, "%s", vAnnotation[nj]->m_pString);
		}
		
		fprintf(fpOut, "\n");
	}

	fclose(fpOut);

	/* output random control probes */
	pRandMat = NULL;
	pRandSort = NULL;
	pRandSortId = NULL;
	pRandMat = DMRANDU(1, nProbeNum);
	DMSORTMERGEA_0(pRandMat, &pRandSort, &pRandSortId);

	sprintf(strOutPath, "%s_ctr.txt", strOutFile);
	fpOut = NULL;
	fpOut = fopen(strOutPath, "wt");
	if(fpOut == NULL)
	{
		printf("Error: Expression_GeneSelection_Output, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	
	if(strcmp(strGeneInfoFile, "NULL") == 0)
	{
		nIncludeInfo = 0;
	}
	else
	{
		nIncludeInfo = 0;
	}

	vSID = pRandSortId->pMatElement;
	vSS = pSortScore->pMatElement;
	if(pFDR != NULL)
	{
		vSF = pFDR->pMatElement;
	}

	fprintf(fpOut, "probeset_id\tline_id\trank\tscore\tFDR\tInfo\n");
	for(ni=0; ni<nOutputNum; ni++)
	{
		if(ni >= nProbeNum)
			break;
		nk = pRandSortId->pMatElement[ni];
		nj = pSortID->pMatElement[nk];

		if(vProbeName[nj] != NULL)
		{
			fprintf(fpOut, "%s\t", vProbeName[nj]->m_pString);
		}
		else
		{
			fprintf(fpOut, "NA\t");
		}
		
		fprintf(fpOut, "%d\t%d\t% 9.7e\t", (nj+1), (nk+1), vSS[nk]);
		
		if(pFDR != NULL)
		{
			fprintf(fpOut, "% 9.7e\t", vSF[nk]);
		}
		else
		{
			fprintf(fpOut, "NA\t");
		}

		if(nIncludeInfo == 1)
		{
			fprintf(fpOut, "%s", vAnnotation[nj]->m_pString);
		}
		
		fprintf(fpOut, "\n");
	}

	fclose(fpOut);


	/* release memory */
	DestroyDoubleMatrix(pRandMat);
	DestroyDoubleMatrix(pRandSort);
	DestroyLongMatrix(pRandSortId);

	for(ni=0; ni<nProbeNum; ni++)
	{
		DeleteString(vProbeName[ni]);
		vProbeName[ni] = NULL;
	}
	free(vProbeName);

	if(nIncludeInfo == 1)
	{
		for(ni=0; ni<nProbeNum; ni++)
		{
			DeleteString(vAnnotation[ni]);
			vAnnotation[ni] = NULL;
		}
		free(vAnnotation);
	}
	
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneRankByLocusLink_Main()                                  */
/*  Gene selection based on combining multiple probesets by Locuslink ID   */
/*  Return PROC_SUCCESS if success.                                        */
/* ----------------------------------------------------------------------- */ 
int Expression_GeneRankByLocusLink_Main(char strParamFile[])
{
	/* define */
	int nProbeNum;
	double dResolution;
	char strScorePath[LINE_LENGTH];
	char strTransform[LINE_LENGTH];
	char strMapPath[LINE_LENGTH];
	char strOutPath[LINE_LENGTH];
	char strLine[MED_LINE_LENGTH];
	char *chSep;
	int nError = 0;
	FILE *fpIn;
	int nMaxIter;
	double dConvth;
	int nOutputNum;

	/* -------------- */
	/* load parameter */
	/* -------------- */
	fpIn = NULL;
	fpIn = fopen(strParamFile, "rt");
	if(fpIn == NULL)
	{
		printf("Error: Expression_GeneRankByLocusLink_Main, cannot open parameter file!\n");
		exit(EXIT_FAILURE);
	}
	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		/* load basic info */
		if(strstr(strLine, "probeset") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			nProbeNum = atoi(chSep);
		}
		else if(strstr(strLine, "score") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strScorePath, chSep);
		}
		else if(strstr(strLine, "transform") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strTransform, chSep);
		}
		else if(strstr(strLine, "resolution") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			dResolution = atof(chSep);
		}
		else if(strstr(strLine, "locuslink") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strMapPath, chSep);
		}
		else if(strstr(strLine, "max") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			nMaxIter = atoi(chSep);
		}
		else if(strstr(strLine, "convergence") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			dConvth = atof(chSep);
		}
		else if((strstr(strLine, "output") == strLine) && (strstr(strLine, "path") != NULL) )
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strOutPath, chSep);
		}
		else if((strstr(strLine, "output") == strLine) && (strstr(strLine, "number") != NULL) )
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			nOutputNum = atoi(chSep);
		}
		
	}

	fclose(fpIn);

	/* analysis */
	Expression_GeneRankByLocusLink_TestPaper(nProbeNum, strScorePath, strTransform,
		dResolution, strMapPath, nMaxIter, dConvth, strOutPath, nOutputNum);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneRankByLocusLink()                                       */
/*  Gene selection based on combining multiple probesets by Locuslink ID   */
/*  Return PROC_SUCCESS if success.                                        */
/* ----------------------------------------------------------------------- */ 
int Expression_GeneRankByLocusLink(int nProbeNum, char strScorePath[],
								   char strTransform[], double dResolution,
								   char strMapPath[], int nMaxIter, 
								   double dPrecision, char strOutPath[],
								   int nOutputNum)
{
	/* define */
	/* score and locuslink map */
	struct DOUBLEMATRIX *pScore;
	struct DOUBLEMATRIX *pWorkScore;
	struct DOUBLEMATRIX *pLocus;
	struct tagString **vProbeName;
	int nLocusNum,nLocusNANum;

	/* sort locus */
	double *vLL;
	long *vLI;
	struct DOUBLEMATRIX *pSortLocus;
	struct LONGMATRIX *pSortIndex;
	struct INTMATRIX *pUniqLocus;

	/* read file */
	FILE *fpIn;
	char strLine[LINE_LENGTH];
	char strLine2[LONG_LINE_LENGTH];
	char strProbe[LINE_LENGTH];
	char strTemp[LINE_LENGTH];
	char *chSep1,*chSep2;

	/* count */
	int ni,nj,nk,nl,nx,ny;
	int nignore;
	double dignoremax;

	/* hierarchical model */
	/* sample mean */
	struct DOUBLEMATRIX *pX;
	/* sample sum of squared residuals */
	struct DOUBLEMATRIX *pS;
	/* sample size */
	struct INTMATRIX *pN;
	/* group map */
	struct INTMATRIX **vLocusMap;

	double *vX,*vS,*vWScore;
	int *vN;
	double dTemp,dSTemp;

	/* final score */
	struct DOUBLEMATRIX *pPostMean;
	struct DOUBLEMATRIX *pSortPostMean;
	struct LONGMATRIX *pSortPostIndex;

	/* ----------- */
	/* init        */
	/* ----------- */
	if(nProbeNum <=0 )
	{
		printf("Warning: Expression_GeneRankByLocusLink, probeset number = 0!\n");
		return PROC_FAILURE;
	}
	
	/* ----------- */
	/* load score  */
	/* ----------- */
	printf("Load original score...\n");
	fpIn = NULL;
	fpIn = fopen(strScorePath, "rt");
	if(fpIn == NULL)
	{
		printf("Error: Expression_GeneRankByLocusLink, cannot open file to load score!\n");
		exit(EXIT_FAILURE);
	}

	pScore = NULL;
	pScore = CreateDoubleMatrix(1, nProbeNum);
	if(pScore == NULL)
	{
		printf("Error: Expression_GeneRankByLocusLink, cannot allocate enough memory for loading score!\n");
		exit(EXIT_FAILURE);
	}

	vProbeName = NULL;
	vProbeName = (struct tagString **)calloc(nProbeNum, sizeof(struct tagString *));
	if(vProbeName == NULL)
	{
		printf("Error: Expression_GeneRankByLocusLink, cannot allocate enough memory for loading probeset id!\n");
		exit(EXIT_FAILURE);
	}

	fgets(strLine, LINE_LENGTH, fpIn);
	ni = 0;
	while(fgets(strLine, LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(ni >= nProbeNum)
		{
			printf("Error: Expression_GeneRankByLocusLink, probeset number not match!\n");
			exit(EXIT_FAILURE);
		}

		sscanf(strLine, "%s %lf", strProbe, (pScore->pMatElement+ni));
		StringAddTail((vProbeName+ni), strProbe);
		ni++;
	}

	fclose(fpIn);

	if(ni != nProbeNum)
	{
		printf("Error: Expression_GeneRankByLocusLink, probeset number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* ----------- */
	/* transform   */
	/* ----------- */
	printf("Transform original score...\n");
	pWorkScore = NULL;
	pWorkScore = Expression_GeneRankByLocusLink_ScoreTransform(pScore, dResolution, strTransform);
	if(pWorkScore == NULL)
	{
		printf("Error: Expression_GeneRankByLocusLink, cannot allocate enough memory for loading working score!\n");
		exit(EXIT_FAILURE);
	}

	/* -------------- */
	/* load locuslink */
	/* -------------- */
	printf("Create locus link map...\n");
	fpIn = NULL;
	fpIn = fopen(strMapPath, "rt");
	if(fpIn == NULL)
	{
		printf("Error: Expression_GeneRankByLocusLink, cannot open file to load locus link!\n");
		exit(EXIT_FAILURE);
	}

	pLocus = NULL;
	pLocus = CreateDoubleMatrix(1, nProbeNum);
	if(pLocus == NULL)
	{
		printf("Error: Expression_GeneRankByLocusLink, cannot allocate enough memory for loading locus link!\n");
		exit(EXIT_FAILURE);
	}

	fgets(strLine2, LONG_LINE_LENGTH, fpIn);
	ni = 0;
	nLocusNANum = 0;
	while(fgets(strLine2, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine2);
		StrTrimRight(strLine2);
		if(strLine2[0] == '\0')
			continue;

		if(ni >= nProbeNum)
		{
			printf("Error: Expression_GeneRankByLocusLink, probeset number not match!\n");
			exit(EXIT_FAILURE);
		}

		chSep1 = strchr(strLine2, '\t');
		*chSep1 = '\0';
		strcpy(strProbe, strLine2);
		StrTrimRight(strProbe);
		if(strcmp(strProbe, vProbeName[ni]->m_pString) != 0)
		{
			printf("Error: locuslink map not match the probeset id in the score file!\n");
			exit(EXIT_FAILURE);
		}

		chSep1++;
		chSep2 = strchr(chSep1, '\t');
		chSep1 = chSep2+1;
		chSep2 = strchr(chSep1, '\t');
		*chSep2 = '\0';
		strcpy(strTemp, chSep1);
		StrTrimRight(strTemp);
		StrTrimLeft(strTemp);
		if(strTemp[0] == '\0')
		{
			nLocusNANum++;
			pLocus->pMatElement[ni] = -nLocusNANum;
		}
		else if(strcmp(strTemp, "---") == 0)
		{
			nLocusNANum++;
			pLocus->pMatElement[ni] = -nLocusNANum;
		}
		else
		{
			pLocus->pMatElement[ni] = atof(strTemp);
		}
		
		ni++;
	}

	fclose(fpIn);

	if(ni != nProbeNum)
	{
		printf("Error: Expression_GeneRankByLocusLink, probeset number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* ------------- */
	/* sort & group  */
	/* ------------- */
	/* sort */
	pSortLocus = NULL;
	pSortIndex = NULL;
	if( DMSORTMERGEA_0(pLocus, &pSortLocus, &pSortIndex) == PROC_FAILURE )
	{
		printf("Error: Expression_GeneRankByLocusLink, sorting locus link failure!\n");
		exit(EXIT_FAILURE);
	}

	/* get locus link group number */
	nLocusNum = 1;
	vLL = pSortLocus->pMatElement;
	for(ni=1; ni<nProbeNum; ni++)
	{
		if(vLL[ni] != vLL[ni-1])
		{
			nLocusNum++;
		}
	}

	/* prepare mean, sst, df, etc. */
	vLI = pSortIndex->pMatElement;
	vWScore = pWorkScore->pMatElement;

	pUniqLocus = CreateIntMatrix(1, nLocusNum);
	pN = CreateIntMatrix(1, nLocusNum);
	pX = CreateDoubleMatrix(1, nLocusNum);
	pS = CreateDoubleMatrix(1, nLocusNum);
	vN = pN->pMatElement;
	vX = pX->pMatElement;
	vS = pS->pMatElement;
	
	vLocusMap = (struct INTMATRIX **)calloc(nLocusNum, sizeof(struct INTMATRIX *));
	
	nj = 0;
	nl = 1;
	ny = 0;
	for(ni=1; ni<nProbeNum; ni++)
	{
		if(vLL[ni] != vLL[ni-1])
		{
			vLocusMap[ny] = CreateIntMatrix(1, nl);
			dTemp = 0.0;
			nignore = 0;
			for(nk=0; nk<nl; nk++)
			{
				nx = (int)(vLI[ni-nl+nk]);
				if(nk==0)
				{
					dignoremax = vWScore[nx];
					nignore = 0;
				}
				else
				{
					if(vWScore[nx] > dignoremax)
					{
						dignoremax = vWScore[nx];
						nignore = nk;
					}
				}
				dTemp += vWScore[nx];
				vLocusMap[ny]->pMatElement[nk] = nx;
			}

			if(nl == 1)
			{
				vX[ny] = dTemp;
				vS[ny] = 0.0;
				vN[ny] = nl;
			}
			else
			{
				dTemp -= dignoremax;
				dTemp /= (double)(nl-1);
				vX[ny] = dTemp;
		
				dSTemp = 0.0;
				for(nk=0; nk<nl; nk++)
				{
					if(nk == nignore)
						continue;

					nx = (int)(vLI[ni-nl+nk]);
					dSTemp += (vWScore[nx]-dTemp)*(vWScore[nx]-dTemp);
				}

				vS[ny] = dSTemp;
				vN[ny] = nl-1;
			}

			pUniqLocus->pMatElement[ny] = (int)(vLL[ni-1]);
			nj = ni;
			nl = 1;
			ny++;
		}
		else
		{
			nl++;
		}
	}

	vLocusMap[ny] = CreateIntMatrix(1, nl);
	dTemp = 0.0;
	nignore = 0;
	for(nk=0; nk<nl; nk++)
	{
		nx = (int)(vLI[ni-nl+nk]);
		if(nk==0)
		{
			dignoremax = vWScore[nx];
			nignore = 0;
		}
		else
		{
			if(vWScore[nx] > dignoremax)
			{
				dignoremax = vWScore[nx];
				nignore = nk;
			}
		}
		dTemp += vWScore[nx];
		vLocusMap[ny]->pMatElement[nk] = nx;
	}

	if(nl == 1)
	{
		vX[ny] = dTemp;
		vS[ny] = 0.0;
		vN[ny] = nl;
	}
	else
	{
		dTemp -= dignoremax;
		dTemp /= (double)(nl-1);
		vX[ny] = dTemp;

		dSTemp = 0.0;
		for(nk=0; nk<nl; nk++)
		{
			if(nk == nignore)
				continue;

			nx = (int)(vLI[ni-nl+nk]);
			dSTemp += (vWScore[nx]-dTemp)*(vWScore[nx]-dTemp);
		}

		vS[ny] = dSTemp;
		vN[ny] = nl-1;
	}

	pUniqLocus->pMatElement[ny] = (int)(vLL[ni-1]);

	ny++;

	if(ny != nLocusNum)
	{
		printf("Error: Expression_GeneRankByLocusLink, locus link group number not match!\n");
		exit(EXIT_FAILURE);
	}
	

	/* -------------- */
	/* posterior var  */
	/* -------------- */
	printf("Calculate new score...\n");
	Expression_GeneRankByLocusLink_EstimateVar(pS, pN);
	for(ni=0; ni<nLocusNum; ni++)
	{
		if(vS[ni] > 0.0)
		{
			vS[ni] = (double)(vN[ni])/vS[ni];
		}
		else
		{
			printf("Error: Expression_GeneRankByLocusLink, try to get precision using zero variance!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* -------------- */
	/* posterior mean */
	/* -------------- */
	pPostMean = Expression_GeneRankByLocusLink_EstimateMean(pX, pS, pN, nMaxIter, dPrecision);

	DestroyDoubleMatrix(pX);
	DestroyDoubleMatrix(pS);
	DestroyIntMatrix(pN);

	Expression_GeneRankByLocusLink_ScoreInverseTransform(pPostMean, strTransform);

	/* -------------- */
	/* rank genes     */
	/* -------------- */
	printf("Rank genes by new score...\n");
	pSortPostMean = NULL;
	pSortPostIndex = NULL;
	if( DMSORTMERGEA_0(pPostMean, &pSortPostMean, &pSortPostIndex) == PROC_FAILURE )
	{
		printf("Error: Expression_GeneRankByLocusLink, sorting final score failure!\n");
		exit(EXIT_FAILURE);
	}

	/* ----------- */
	/* output      */
	/* ----------- */
	printf("Link annotations...\n");

	/* sortpostindex -> sortlocus -> locusmap -> probeset -> annotation */
	Expression_GeneRankByLocusLink_Output(nLocusNum, pSortPostMean, pSortPostIndex, 
		pUniqLocus, vLocusMap, nProbeNum, vProbeName, pScore, strMapPath, strOutPath, nOutputNum);

	/* ----------- */
	/* clear mem   */
	/* ----------- */
	DestroyIntMatrix(pUniqLocus);
	DestroyDoubleMatrix(pScore);
	DestroyDoubleMatrix(pWorkScore);
	DestroyDoubleMatrix(pLocus);

	for(ni=0; ni<nProbeNum; ni++)
	{
		DeleteString(vProbeName[ni]);
		vProbeName[ni] = NULL;
	}
	free(vProbeName);

	DestroyDoubleMatrix(pSortLocus);
	DestroyLongMatrix(pSortIndex);

	for(ni=0; ni<nLocusNum; ni++)
	{
		DestroyIntMatrix(vLocusMap[ni]);
		vLocusMap[ni] = NULL;
	}
	free(vLocusMap);

	DestroyDoubleMatrix(pPostMean);
	DestroyDoubleMatrix(pSortPostMean);
	DestroyLongMatrix(pSortPostIndex);

	printf("Done!\n");

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneRankByLocusLink()                                       */
/*  Gene selection based on combining multiple probesets by Locuslink ID   */
/*  Return PROC_SUCCESS if success.                                        */
/* ----------------------------------------------------------------------- */ 
int Expression_GeneRankByLocusLink_TestPaper(int nProbeNum, char strScorePath[],
								   char strTransform[], double dResolution,
								   char strMapPath[], int nMaxIter, 
								   double dPrecision, char strOutPath[],
								   int nOutputNum)
{
	/* define */
	/* score and locuslink map */
	struct DOUBLEMATRIX *pScore;
	struct DOUBLEMATRIX *pWorkScore;
	struct DOUBLEMATRIX *pLocus;
	struct tagString **vProbeName;
	int nLocusNum,nLocusNANum;

	/* sort locus */
	double *vLL;
	long *vLI;
	struct DOUBLEMATRIX *pSortLocus;
	struct LONGMATRIX *pSortIndex;
	struct INTMATRIX *pUniqLocus;

	/* read file */
	FILE *fpIn;
	char strLine[LINE_LENGTH];
	char strLine2[LONG_LINE_LENGTH];
	char strProbe[LINE_LENGTH];
	char strTemp[LINE_LENGTH];
	char *chSep1,*chSep2;

	/* count */
	int ni,nj,nk,nl,nx,ny;
	int nignore;
	double dignoremax;

	/* hierarchical model */
	/* sample mean */
	struct DOUBLEMATRIX *pX;
	/* sample sum of squared residuals */
	struct DOUBLEMATRIX *pS;
	/* sample size */
	struct INTMATRIX *pN;
	/* group map */
	struct INTMATRIX **vLocusMap;

	double *vX,*vS,*vWScore;
	int *vN;
	double dTemp,dSTemp;

	/* final score */
	struct DOUBLEMATRIX *pPostMean;
	struct DOUBLEMATRIX *pSortPostMean;
	struct LONGMATRIX *pSortPostIndex;

	/* ----------- */
	/* init        */
	/* ----------- */
	if(nProbeNum <=0 )
	{
		printf("Warning: Expression_GeneRankByLocusLink, probeset number = 0!\n");
		return PROC_FAILURE;
	}
	
	/* ----------- */
	/* load score  */
	/* ----------- */
	printf("Load original score...\n");
	fpIn = NULL;
	fpIn = fopen(strScorePath, "rt");
	if(fpIn == NULL)
	{
		printf("Error: Expression_GeneRankByLocusLink, cannot open file to load score!\n");
		exit(EXIT_FAILURE);
	}

	pScore = NULL;
	pScore = CreateDoubleMatrix(1, nProbeNum);
	if(pScore == NULL)
	{
		printf("Error: Expression_GeneRankByLocusLink, cannot allocate enough memory for loading score!\n");
		exit(EXIT_FAILURE);
	}

	vProbeName = NULL;
	vProbeName = (struct tagString **)calloc(nProbeNum, sizeof(struct tagString *));
	if(vProbeName == NULL)
	{
		printf("Error: Expression_GeneRankByLocusLink, cannot allocate enough memory for loading probeset id!\n");
		exit(EXIT_FAILURE);
	}

	/* fgets(strLine, LINE_LENGTH, fpIn); */
	ni = 0;
	while(fgets(strLine, LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(ni >= nProbeNum)
		{
			printf("Error: Expression_GeneRankByLocusLink, probeset number not match!\n");
			exit(EXIT_FAILURE);
		}

		sscanf(strLine, "%lf", (pScore->pMatElement+ni));
		/* StringAddTail((vProbeName+ni), strProbe); */
		ni++;
	}

	fclose(fpIn);

	if(ni != nProbeNum)
	{
		printf("Error: Expression_GeneRankByLocusLink, probeset number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* ----------- */
	/* transform   */
	/* ----------- */
	printf("Transform original score...\n");
	pWorkScore = NULL;
	pWorkScore = Expression_GeneRankByLocusLink_ScoreTransform(pScore, dResolution, strTransform);
	if(pWorkScore == NULL)
	{
		printf("Error: Expression_GeneRankByLocusLink, cannot allocate enough memory for loading working score!\n");
		exit(EXIT_FAILURE);
	}

	/* -------------- */
	/* load locuslink */
	/* -------------- */
	printf("Create locus link map...\n");
	fpIn = NULL;
	fpIn = fopen(strMapPath, "rt");
	if(fpIn == NULL)
	{
		printf("Error: Expression_GeneRankByLocusLink, cannot open file to load locus link!\n");
		exit(EXIT_FAILURE);
	}

	pLocus = NULL;
	pLocus = CreateDoubleMatrix(1, nProbeNum);
	if(pLocus == NULL)
	{
		printf("Error: Expression_GeneRankByLocusLink, cannot allocate enough memory for loading locus link!\n");
		exit(EXIT_FAILURE);
	}

	fgets(strLine2, LONG_LINE_LENGTH, fpIn);
	ni = 0;
	nLocusNANum = 0;
	while(fgets(strLine2, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine2);
		StrTrimRight(strLine2);
		if(strLine2[0] == '\0')
			continue;

		if(ni >= nProbeNum)
		{
			printf("Error: Expression_GeneRankByLocusLink, probeset number not match!\n");
			exit(EXIT_FAILURE);
		}

		chSep1 = strchr(strLine2, '\t');
		*chSep1 = '\0';
		strcpy(strProbe, strLine2);
		StrTrimRight(strProbe);
		/* if(strcmp(strProbe, vProbeName[ni]->m_pString) != 0)
		{
			printf("Error: locuslink map not match the probeset id in the score file!\n");
			exit(EXIT_FAILURE);
		} */
		StringAddTail(vProbeName+ni, strProbe);

		chSep1++;
		chSep2 = strchr(chSep1, '\t');
		chSep1 = chSep2+1;
		chSep2 = strchr(chSep1, '\t');
		*chSep2 = '\0';
		strcpy(strTemp, chSep1);
		StrTrimRight(strTemp);
		StrTrimLeft(strTemp);
		if(strTemp[0] == '\0')
		{
			nLocusNANum++;
			pLocus->pMatElement[ni] = -nLocusNANum;
		}
		else if(strcmp(strTemp, "---") == 0)
		{
			nLocusNANum++;
			pLocus->pMatElement[ni] = -nLocusNANum;
		}
		else
		{
			pLocus->pMatElement[ni] = atof(strTemp);
		}
		
		ni++;
	}

	fclose(fpIn);

	if(ni != nProbeNum)
	{
		printf("Error: Expression_GeneRankByLocusLink, probeset number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* ------------- */
	/* sort & group  */
	/* ------------- */
	/* sort */
	pSortLocus = NULL;
	pSortIndex = NULL;
	if( DMSORTMERGEA_0(pLocus, &pSortLocus, &pSortIndex) == PROC_FAILURE )
	{
		printf("Error: Expression_GeneRankByLocusLink, sorting locus link failure!\n");
		exit(EXIT_FAILURE);
	}

	/* get locus link group number */
	nLocusNum = 1;
	vLL = pSortLocus->pMatElement;
	for(ni=1; ni<nProbeNum; ni++)
	{
		if(vLL[ni] != vLL[ni-1])
		{
			nLocusNum++;
		}
	}

	/* prepare mean, sst, df, etc. */
	vLI = pSortIndex->pMatElement;
	vWScore = pWorkScore->pMatElement;

	pUniqLocus = CreateIntMatrix(1, nLocusNum);
	pN = CreateIntMatrix(1, nLocusNum);
	pX = CreateDoubleMatrix(1, nLocusNum);
	pS = CreateDoubleMatrix(1, nLocusNum);
	vN = pN->pMatElement;
	vX = pX->pMatElement;
	vS = pS->pMatElement;
	
	vLocusMap = (struct INTMATRIX **)calloc(nLocusNum, sizeof(struct INTMATRIX *));
	
	nj = 0;
	nl = 1;
	ny = 0;
	for(ni=1; ni<nProbeNum; ni++)
	{
		if(vLL[ni] != vLL[ni-1])
		{
			vLocusMap[ny] = CreateIntMatrix(1, nl);
			dTemp = 0.0;
			nignore = 0;
			for(nk=0; nk<nl; nk++)
			{
				nx = (int)(vLI[ni-nl+nk]);
				if(nk==0)
				{
					dignoremax = vWScore[nx];
					nignore = 0;
				}
				else
				{
					if(vWScore[nx] > dignoremax)
					{
						dignoremax = vWScore[nx];
						nignore = nk;
					}
				}
				dTemp += vWScore[nx];
				vLocusMap[ny]->pMatElement[nk] = nx;
			}

			if(nl == 1)
			{
				vX[ny] = dTemp;
				vS[ny] = 0.0;
				vN[ny] = nl;
			}
			else
			{
				dTemp -= dignoremax;
				dTemp /= (double)(nl-1);
				vX[ny] = dTemp;
		
				dSTemp = 0.0;
				for(nk=0; nk<nl; nk++)
				{
					if(nk == nignore)
						continue;

					nx = (int)(vLI[ni-nl+nk]);
					dSTemp += (vWScore[nx]-dTemp)*(vWScore[nx]-dTemp);
				}

				vS[ny] = dSTemp;
				vN[ny] = nl-1;
			}

			pUniqLocus->pMatElement[ny] = (int)(vLL[ni-1]);
			nj = ni;
			nl = 1;
			ny++;
		}
		else
		{
			nl++;
		}
	}

	vLocusMap[ny] = CreateIntMatrix(1, nl);
	dTemp = 0.0;
	nignore = 0;
	for(nk=0; nk<nl; nk++)
	{
		nx = (int)(vLI[ni-nl+nk]);
		if(nk==0)
		{
			dignoremax = vWScore[nx];
			nignore = 0;
		}
		else
		{
			if(vWScore[nx] > dignoremax)
			{
				dignoremax = vWScore[nx];
				nignore = nk;
			}
		}
		dTemp += vWScore[nx];
		vLocusMap[ny]->pMatElement[nk] = nx;
	}

	if(nl == 1)
	{
		vX[ny] = dTemp;
		vS[ny] = 0.0;
		vN[ny] = nl;
	}
	else
	{
		dTemp -= dignoremax;
		dTemp /= (double)(nl-1);
		vX[ny] = dTemp;

		dSTemp = 0.0;
		for(nk=0; nk<nl; nk++)
		{
			if(nk == nignore)
				continue;

			nx = (int)(vLI[ni-nl+nk]);
			dSTemp += (vWScore[nx]-dTemp)*(vWScore[nx]-dTemp);
		}

		vS[ny] = dSTemp;
		vN[ny] = nl-1;
	}

	pUniqLocus->pMatElement[ny] = (int)(vLL[ni-1]);

	ny++;

	if(ny != nLocusNum)
	{
		printf("Error: Expression_GeneRankByLocusLink, locus link group number not match!\n");
		exit(EXIT_FAILURE);
	}
	

	/* -------------- */
	/* posterior var  */
	/* -------------- */
	printf("Calculate new score...\n");
	Expression_GeneRankByLocusLink_EstimateVar(pS, pN);
	for(ni=0; ni<nLocusNum; ni++)
	{
		if(vS[ni] > 0.0)
		{
			vS[ni] = (double)(vN[ni])/vS[ni];
		}
		else
		{
			printf("Error: Expression_GeneRankByLocusLink, try to get precision using zero variance!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* -------------- */
	/* posterior mean */
	/* -------------- */
	pPostMean = Expression_GeneRankByLocusLink_EstimateMean(pX, pS, pN, nMaxIter, dPrecision);

	DestroyDoubleMatrix(pX);
	DestroyDoubleMatrix(pS);
	DestroyIntMatrix(pN);

	Expression_GeneRankByLocusLink_ScoreInverseTransform(pPostMean, strTransform);

	/* -------------- */
	/* rank genes     */
	/* -------------- */
	printf("Rank genes by new score...\n");
	pSortPostMean = NULL;
	pSortPostIndex = NULL;
	if( DMSORTMERGEA_0(pPostMean, &pSortPostMean, &pSortPostIndex) == PROC_FAILURE )
	{
		printf("Error: Expression_GeneRankByLocusLink, sorting final score failure!\n");
		exit(EXIT_FAILURE);
	}

	/* ----------- */
	/* output      */
	/* ----------- */
	printf("Link annotations...\n");

	/* sortpostindex -> sortlocus -> locusmap -> probeset -> annotation */
	Expression_GeneRankByLocusLink_Output(nLocusNum, pSortPostMean, pSortPostIndex, 
		pUniqLocus, vLocusMap, nProbeNum, vProbeName, pScore, strMapPath, strOutPath, nOutputNum);

	/* ----------- */
	/* clear mem   */
	/* ----------- */
	DestroyIntMatrix(pUniqLocus);
	DestroyDoubleMatrix(pScore);
	DestroyDoubleMatrix(pWorkScore);
	DestroyDoubleMatrix(pLocus);

	for(ni=0; ni<nProbeNum; ni++)
	{
		DeleteString(vProbeName[ni]);
		vProbeName[ni] = NULL;
	}
	free(vProbeName);

	DestroyDoubleMatrix(pSortLocus);
	DestroyLongMatrix(pSortIndex);

	for(ni=0; ni<nLocusNum; ni++)
	{
		DestroyIntMatrix(vLocusMap[ni]);
		vLocusMap[ni] = NULL;
	}
	free(vLocusMap);

	DestroyDoubleMatrix(pPostMean);
	DestroyDoubleMatrix(pSortPostMean);
	DestroyLongMatrix(pSortPostIndex);

	printf("Done!\n");

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneRankByLocusLink_ScoreTransform()                        */
/*  transform original scores to working scores so that normality is       */
/*  appropriate. Acceptable transformations are logit, rank, identity      */
/*  return the transformed matrix.                                         */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *Expression_GeneRankByLocusLink_ScoreTransform(struct DOUBLEMATRIX *pScore, 
										double dResolution, char strTransform[])
{
	/* define */
	struct DOUBLEMATRIX *pNewScore;
	int ni,nj,nk,nl,nx;
	char strTransformType[LINE_LENGTH];
	double *pEle1,*pEle2;
	double dTemp;
	struct LONGMATRIX *pSortIndex;
	struct DOUBLEMATRIX *pSortScore;
	long *vId;
	

	/* init */
	if(pScore == NULL)
		return NULL;

	/* transform */
	strcpy(strTransformType, strTransform);
	StrMakeUpper(strTransformType);

	/* logit */
	if(strcmp(strTransformType, "LOGIT") == 0)
	{
		pNewScore = NULL;
		pNewScore = CreateDoubleMatrix(pScore->nHeight, pScore->nWidth);
		if(pNewScore == NULL)
		{
			printf("Error: Expression_GeneRankByLocusLink_ScoreTransform, cannot allocate enough memory for storing working score!\n");
			exit(EXIT_FAILURE);
		}

		pEle1 = pNewScore->pMatElement;
		pEle2 = pScore->pMatElement;
		dTemp = dResolution/2.0;
		for(ni=0; ni<pScore->nHeight; ni++)
		{
			for(nj=0; nj<pScore->nWidth; nj++)
			{
				if(*pEle2 < dTemp)
				{
					*pEle1 = log(dTemp/(1.0-dTemp));
				}
				else if(*pEle2 > 1.0-dTemp)
				{
					*pEle1 = log((1.0-dTemp)/dTemp);
				}
				else
				{
					*pEle1 = log((*pEle2)/(1.0-*pEle2));
				}
				pEle1++;
				pEle2++;
			}
		}
	}
	/* rank */
	else if(strcmp(strTransformType, "RANK") == 0)
	{
		pSortIndex = NULL;
		pSortScore = NULL;
		if( DMSORTMERGEA_0(pScore, &pSortScore, &pSortIndex) == PROC_FAILURE )
		{
			printf("Error: Expression_GeneRankByLocusLink_ScoreTransform, sorting score failure!\n");
			exit(EXIT_FAILURE);
		}

		pNewScore = NULL;
		pNewScore = CreateDoubleMatrix(pScore->nHeight, pScore->nWidth);
		if(pNewScore == NULL)
		{
			printf("Error: Expression_GeneRankByLocusLink_ScoreTransform, cannot allocate enough memory for storing working score!\n");
			exit(EXIT_FAILURE);
		}

		ni=0;
		pEle1 = pNewScore->pMatElement;
		pEle2 = pSortScore->pMatElement;
		vId = pSortIndex->pMatElement;
		nl = 1;
		for(nj=1; nj<pScore->nWidth; nj++)
		{
			if(pEle2[nj] > pEle2[nj-1])
			{
				for(nx=0; nx<nl; nx++)
				{
					nk = (int)(vId[nj-1-nx]);
					pEle1[nk] = (double)(ni+nj-1)/2.0;
				}
				
				ni = nj;
				nl = 1;
			}
			else
			{
				nl++;
			}
		}
		for(nx=0; nx<nl; nx++)
		{
			nk = (int)(vId[nj-1-nx]);
			pEle1[nk] = (double)(ni+nj-1)/2.0;
		}

		for(nj=0; nj<pNewScore->nWidth; nj++)
		{
			pEle1[nj] = (pEle1[nj]+0.5)/(double)(pNewScore->nWidth);
			pEle1[nj] = log(pEle1[nj]/(1.0-pEle1[nj]));
		}

		DestroyDoubleMatrix(pSortScore);
		DestroyLongMatrix(pSortIndex);
	}
	/* identity */
	else
	{
		pNewScore = NULL;
		pNewScore = DMCLONE(pScore);
	}

	/* return */
	return pNewScore;
}

/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneRankByLocusLink_ScoreInverseTransform()                 */
/* ----------------------------------------------------------------------- */ 
int Expression_GeneRankByLocusLink_ScoreInverseTransform(struct DOUBLEMATRIX *pScore, 
														 char strTransform[])
{
	/* define */
	char strTransformType[LINE_LENGTH];
	double *vS;
	double dTemp;
	int ni;

	/* init */
	if(pScore == NULL)
		return PROC_FAILURE;

	/* transform */
	strcpy(strTransformType, strTransform);
	StrMakeUpper(strTransformType);

	/* logit */
	if(strcmp(strTransformType, "LOGIT") == 0)
	{
		vS = pScore->pMatElement;
		for(ni=0; ni<pScore->nWidth; ni++)
		{
			dTemp = exp(vS[ni]);
			vS[ni] = dTemp/(1.0+dTemp);
		}
	}
	else
	{
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneRankByLocusLink_EstimateVar()                           */
/*  estimate the variance                                                  */
/* ----------------------------------------------------------------------- */ 
int Expression_GeneRankByLocusLink_EstimateVar(struct DOUBLEMATRIX *pS, 
							struct INTMATRIX *pN)
{
	/* define */
	double dMu;
	double dM2;
	double dDfTot;
	double dDf;
	double *vS;
	int *vN;
	double dTot;
	double dTemp;
	int ni;
	double dB;

	/* init */
	if((pS == NULL) || (pN == NULL))
	{
		printf("Warning: Expression_GeneRankByLocusLink_EstimateVar, no S or N matrix for variance estimation!\n");
		return PROC_FAILURE;
	}

	/* get EMU and EMU^2 */
	dMu = 0.0;
	dM2 = 0.0;
	dTot = 0.0;
	dDfTot = 0.0;
	vS = pS->pMatElement;
	vN = pN->pMatElement;
	for(ni=0; ni<pS->nWidth; ni++)
	{
		if(vN[ni] <= 1)
			continue;
		
		dDf = (double)(vN[ni]-1);
		dDfTot += dDf;
		dMu += vS[ni];

		vS[ni] /= dDf;
		dTemp = 2.0*vS[ni]*vS[ni]/dDf;
		dM2 += dTemp/((2.0/dDf)*(2.0/dDf+1.0));

		dTot += 1.0;
	}
		

	if(dTot > 0.0)
	{
		dMu /= dDfTot;
		dM2 /= dTot;
	}
	else
	{
		printf("Error: Expression_GeneRankByLocusLink_EstimateVar, pooled variance not estimable!\n");
		exit(EXIT_FAILURE);
	}

	/* shrink */
	for(ni=0; ni<pS->nWidth; ni++)
	{
		if(vN[ni] <= 1)
		{
			vS[ni] = dMu;
		}
		else
		{
			dDf = (double)(vN[ni]-1);
			dTemp = 2.0/dDf;
			dB = dTemp*dM2/((dTemp+1.0)*dM2-dMu*dMu);
			if(dB > 1.0) 
				dB = 1.0;
			else if(dB < 1e-2)
				dB = 1e-2;
			vS[ni] = (1.0-dB)*vS[ni]+dB*dMu;
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneRankByLocusLink_EstimateVar0()                          */
/*  estimate the variance                                                  */
/* ----------------------------------------------------------------------- */ 
int Expression_GeneRankByLocusLink_EstimateVar0(struct DOUBLEMATRIX *pS, 
							struct INTMATRIX *pN)
{
	/* define */
	double dSMean;
	double dSVar;
	double dSNum;
	double dDf;
	double dB;
	
	double *vS;
	int *vN;
	int ni;

	/* init */
	if((pS == NULL) || (pN == NULL))
	{
		printf("Warning: Expression_GeneRankByLocusLink_EstimateVar, no S or N matrix for variance estimation!\n");
		return PROC_FAILURE;
	}

	/* get mean of S */
	dSMean = 0.0;
	dSNum = 0.0;
	dDf = 0.0;
	vS = pS->pMatElement;
	vN = pN->pMatElement;
	for(ni=0; ni<pS->nWidth; ni++)
	{
		if(vN[ni] <= 1)
			continue;
		
		dSMean += vS[ni];
		dSNum += 1.0;
		dDf += vN[ni]-1.0;
	}
	if(dDf > 0.0)
	{
		dSMean /= dDf;
	}
	else
	{
		printf("Error: Expression_GeneRankByLocusLink_EstimateVar, pooled variance not estimable!\n");
		exit(EXIT_FAILURE);
	}

	/* get var of S */
	dSVar = 0.0;
	for(ni=0; ni<pS->nWidth; ni++)
	{
		if(vN[ni] <= 1)
			continue;
		
		vS[ni] /= (double)(vN[ni]-1);
		dSVar += (vS[ni]-dSMean)*(vS[ni]-dSMean);
	}
	if(dSNum > 1.0)
	{
		dSVar /= (dSNum-1.0);
	}
	else
	{
		printf("Error: Expression_GeneRankByLocusLink_EstimateVar, var of pooled variance not estimable!\n");
		exit(EXIT_FAILURE);
	}

	/* shrink */
	for(ni=0; ni<pS->nWidth; ni++)
	{
		if(vN[ni] <= 1)
		{
			vS[ni] = dSMean;
		}
		else
		{
			dB = 2*vS[ni]*vS[ni]/(double)(vN[ni]-1)/dSVar;
			if(dB > 1.0) 
				dB = 1.0;
			else if(dB < 1e-2)
				dB = 1e-2;
			vS[ni] = (1.0-dB)*vS[ni]+dB*dSMean;
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneRankByLocusLink_EstimateMean()                          */
/*  estimate the variance                                                  */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *Expression_GeneRankByLocusLink_EstimateMean(struct DOUBLEMATRIX *pX, 
							struct DOUBLEMATRIX *pS, struct INTMATRIX *pN,
							int nMaxIter, double dPrecision)
{
	/* define */
	struct DOUBLEMATRIX *pPostMean;

	double dMu[2];
	double dVar[2];
	
	double *vX,*vY,*vS;

	int ni,niter;
	int nProbeNum;
	double dError;

	double dTemp,dPrec0;

	/* init */
	if( (pX == NULL) || (pS == NULL) || (pN == NULL) )
	{
		printf("Warning: Expression_GeneRankByLocusLink_EstimateMean, parameter matrix does not exist!\n");
		return NULL;
	}
	pPostMean = NULL;
	pPostMean = CreateDoubleMatrix(pX->nHeight, pX->nWidth);
	if(pPostMean == NULL)
	{
		printf("Error: Expression_GeneRankByLocusLink_EstimateMean, cannot create matrix for EM!\n");
		exit(EXIT_FAILURE);
	}

	vX = pX->pMatElement;
	vY = pPostMean->pMatElement;
	vS = pS->pMatElement;
	
	for(ni=0; ni<2; ni++)
	{
		dMu[ni] = 0.0;
		dVar[ni] = 0.0;
	}

	nProbeNum = pX->nWidth;

	/* initial estimate */
	dMu[0] = 0.0;
	for(ni=0; ni<nProbeNum; ni++)
	{
		dMu[0] += vX[ni];
	}
	dMu[0] /= (double)nProbeNum;

	dVar[0] = 0.0;
	for(ni=0; ni<nProbeNum; ni++)
	{
		dVar[0] += (vX[ni]-dMu[0])*(vX[ni]-dMu[0]);
	}
	dVar[0] /= (double)nProbeNum;

	/* EM */
	dError = 1e6;
	niter = 0;
	while(dError >= dPrecision)
	{
		if(niter >= nMaxIter)
			break;

		if(niter%10 == 0)
			printf("iter %d...\n", niter);
		
		/* E-step */
		dMu[1] = 0.0;
		dVar[1] = 0.0;
		dPrec0 = 1.0/dVar[0];
		for(ni=0; ni<nProbeNum; ni++)
		{
			dTemp = vS[ni]+dPrec0;
			vY[ni] = (vS[ni]*vX[ni]+dPrec0*dMu[0])/dTemp;
			dMu[1] += vY[ni];
			dVar[1] += 1.0/dTemp+vY[ni]*vY[ni];
		}
		
		/* M-step */
		dMu[1] /= (double)nProbeNum;
		dVar[1] = dVar[1]/(double)nProbeNum - dMu[1]*dMu[1];

		/* Error */
		dError = fabs(dMu[1]-dMu[0]);
		dTemp = fabs(dVar[1]-dVar[0]);
		if(dTemp > dError)
			dError = dTemp;

		/* next step */
		dMu[0] = dMu[1];
		dVar[0] = dVar[1];
		niter++;
	}

	/* updata posterior mean */
	dPrec0 = 1.0/dVar[0];
	for(ni=0; ni<nProbeNum; ni++)
	{
		dTemp = vS[ni]+dPrec0;
		vY[ni] = (vS[ni]*vX[ni]+dPrec0*dMu[0])/dTemp;
	}

	/* return */
	return pPostMean;
}

/* ----------------------------------------------------------------------- */ 
/*  Expression_GeneRankByLocusLink_Output()                                */
/*  output the selected gene                                               */
/* ----------------------------------------------------------------------- */ 
int Expression_GeneRankByLocusLink_Output(int nLocusNum, 
		struct DOUBLEMATRIX *pSortScore, struct LONGMATRIX *pSortIndex, 
		struct INTMATRIX *pLocus, struct INTMATRIX **vLocusMap, 
		int nProbeNum, struct tagString **vProbeName, 
		struct DOUBLEMATRIX *pScore,
		char strMapPath[], char strOutPath[], int nOutputNum)
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;

	struct tagString **vAnnotation;
	int ni,nj,nLineId,nk;
	char strLine[LONG_LINE_LENGTH];
	char strProbe[LINE_LENGTH];
	char *chSep;
	struct DOUBLEMATRIX *pGroupScore;
	struct DOUBLEMATRIX *pSortGroupScore;
	struct LONGMATRIX *pSortGroupIndex;
	struct INTMATRIX *pGroupLine;

	/* for paper test */
	/* struct INTMATRIX *pProbeRank;
	struct INTMATRIX *pLocusId;
	struct DOUBLEMATRIX *pNewScore;
	pProbeRank = NULL;
	pProbeRank = CreateIntMatrix(nProbeNum, 1);
	if(pProbeRank == NULL)
	{
		printf("Error: Expression_GeneRankByLocusLink_Output, cannot open restore rank at probelevel!\n");
		exit(EXIT_FAILURE);
	}
	pLocusId = NULL;
	pLocusId = CreateIntMatrix(nProbeNum, 1);
	if(pLocusId == NULL)
	{
		printf("Error: Expression_GeneRankByLocusLink_Output, cannot open restore rank at probelevel!\n");
		exit(EXIT_FAILURE);
	}
	pNewScore = NULL;
	pNewScore = CreateDoubleMatrix(nProbeNum, 1);
	if(pNewScore == NULL)
	{
		printf("Error: Expression_GeneRankByLocusLink_Output, cannot open restore rank at probelevel!\n");
		exit(EXIT_FAILURE);
	} */

	/* load annotation */
	fpIn = NULL;
	fpIn = fopen(strMapPath, "rt");
	if(fpIn == NULL)
	{
		printf("Error: Expression_GeneRankByLocusLink_Output, cannot open file to load annotation!\n");
		exit(EXIT_FAILURE);
	}

	vAnnotation = NULL;
	vAnnotation = (struct tagString **)calloc(nProbeNum, sizeof(struct tagString *));
	if(vAnnotation == NULL)
	{
		printf("Error: Expression_GeneRankByLocusLink_Output, cannot allocate enough memory for loading annotation!\n");
		exit(EXIT_FAILURE);
	}

	fgets(strLine, LONG_LINE_LENGTH, fpIn);
	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(ni >= nProbeNum)
		{
			printf("Error: Expression_GeneRankByLocusLink_Output, probeset number not match!\n");
			exit(EXIT_FAILURE);
		}

		sscanf(strLine, "%s", strProbe);
		if(strcmp(strProbe, vProbeName[ni]->m_pString) != 0)
		{
			printf("Error: Expression_GeneRankByLocusLink_Output, annotation order not match original score order!\n");
			exit(EXIT_FAILURE);
		}
		chSep = strchr(strLine, '\t');
		chSep++;
		StringAddTail((vAnnotation+ni), chSep);
		ni++;
	}

	fclose(fpIn);

	if(ni != nProbeNum)
	{
		printf("Error: Expression_GeneRankByLocusLink_Output, probeset number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* write file */
	fpOut = NULL;
	fpOut = fopen(strOutPath, "wt");
	if(fpOut == NULL)
	{
		printf("Error: Expression_GeneRankByLocusLink_Output, cannot open file to write annotation!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "Rank\tProbe_Line\tLocuslink\tProbeset\tNew_Score\tOld_Score\tAnnotation\n");
	for(ni=0; ni<nOutputNum; ni++)
	{
		nLineId = (int)(pSortIndex->pMatElement[ni]);
		pGroupScore = NULL;
		pSortGroupScore = NULL;
		pSortGroupIndex = NULL;
		pGroupLine = NULL;
		pGroupScore = CreateDoubleMatrix(1, vLocusMap[nLineId]->nWidth);
		pGroupLine = CreateIntMatrix(1, vLocusMap[nLineId]->nWidth);

		for(nj=0; nj<vLocusMap[nLineId]->nWidth; nj++)
		{
			nk = vLocusMap[nLineId]->pMatElement[nj];
			pGroupScore->pMatElement[nj] = pScore->pMatElement[nk];
			pGroupLine->pMatElement[nj] = nk;
		}

		if( DMSORTMERGEA_0(pGroupScore, &pSortGroupScore, &pSortGroupIndex) == PROC_FAILURE )
		{
			printf("Error: Expression_GeneRankByLocusLink_Output, sorting score failure!\n");
			exit(EXIT_FAILURE);
		}

		for(nj=0; nj<vLocusMap[nLineId]->nWidth; nj++)
		{
			nk = pSortGroupIndex->pMatElement[nj];
			nk = pGroupLine->pMatElement[nk];
			fprintf(fpOut, "%d\t%d\t%d\t%s\t% 9.7e\t% 9.7e\t%s\n", (ni+1), nk,
				pLocus->pMatElement[nLineId], vProbeName[nk]->m_pString,
				pSortScore->pMatElement[ni], pScore->pMatElement[nk],
				vAnnotation[nk]->m_pString);
		}

		
		DestroyDoubleMatrix(pGroupScore);
		DestroyDoubleMatrix(pSortGroupScore);
		DestroyLongMatrix(pSortGroupIndex);
		DestroyIntMatrix(pGroupLine);
	}

	fclose(fpOut);


	/* for paper test */
	/* nx = 0;
	for(ni=0; ni<pSortIndex->nWidth; ni++)
	{
		nLineId = (int)(pSortIndex->pMatElement[ni]);
		pGroupScore = NULL;
		pSortGroupScore = NULL;
		pSortGroupIndex = NULL;
		pGroupLine = NULL;
		pGroupScore = CreateDoubleMatrix(1, vLocusMap[nLineId]->nWidth);
		pGroupLine = CreateIntMatrix(1, vLocusMap[nLineId]->nWidth);

		for(nj=0; nj<vLocusMap[nLineId]->nWidth; nj++)
		{
			nk = vLocusMap[nLineId]->pMatElement[nj];
			pGroupScore->pMatElement[nj] = pScore->pMatElement[nk];
			pGroupLine->pMatElement[nj] = nk;
		}

		if( DMSORTMERGEA_0(pGroupScore, &pSortGroupScore, &pSortGroupIndex) == PROC_FAILURE )
		{
			printf("Error: Expression_GeneRankByLocusLink_Output, sorting score failure!\n");
			exit(EXIT_FAILURE);
		}

		for(nj=0; nj<vLocusMap[nLineId]->nWidth; nj++)
		{
			nk = pSortGroupIndex->pMatElement[nj];
			nk = pGroupLine->pMatElement[nk];
			pProbeRank->pMatElement[nk] = nx+1;
			pLocusId->pMatElement[nk] = pLocus->pMatElement[nLineId];
			pNewScore->pMatElement[nk] = pSortScore->pMatElement[ni];
			nx++;
		}

		
		DestroyDoubleMatrix(pGroupScore);
		DestroyDoubleMatrix(pSortGroupScore);
		DestroyLongMatrix(pSortGroupIndex);
		DestroyIntMatrix(pGroupLine);
	}

	if(nx!=nProbeNum)
	{
		printf("Error: Expression_GeneRankByLocusLink_Output, probenumber not match!\n");
		exit(EXIT_FAILURE);
	}

	sprintf(strLine, "%s_locrank", strOutPath);
	IMSAVE(pProbeRank, strLine);
	DestroyIntMatrix(pProbeRank);
	sprintf(strLine, "%s_locid", strOutPath);
	IMSAVE(pLocusId, strLine);
	DestroyIntMatrix(pLocusId);
	sprintf(strLine, "%s", strOutPath);
	DMSAVE(pNewScore, strLine);
	DestroyDoubleMatrix(pNewScore);
	/* for paper test */


	/* release memory */
	for(ni=0; ni<nProbeNum; ni++)
	{
		DeleteString(vAnnotation[ni]);
		vAnnotation[ni] = NULL;
	}
	free(vAnnotation);


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Expression_PowExpressLocusLink_Main()                                  */
/*  This is the algorithm used in powerexpress.                            */
/*  Gene selection based on combining multiple probesets by Locuslink ID   */
/*  Return PROC_SUCCESS if success.                                        */
/* ----------------------------------------------------------------------- */ 
int Expression_PowExpressLocusLink_Main(char strParamFile[])
{
	/* define */
	int nProbeNum;
	double dResolution;
	char strScorePath[LINE_LENGTH];
	char strTransform[LINE_LENGTH];
	char strMapPath[LINE_LENGTH];
	char strOutPath[LINE_LENGTH];
	char strLine[MED_LINE_LENGTH];
	char *chSep;
	int nError = 0;
	FILE *fpIn;
	int nOutputNum;

	/* -------------- */
	/* load parameter */
	/* -------------- */
	fpIn = NULL;
	fpIn = fopen(strParamFile, "rt");
	if(fpIn == NULL)
	{
		printf("Error: Expression_GeneRankByLocusLink_Main, cannot open parameter file!\n");
		exit(EXIT_FAILURE);
	}
	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		/* load basic info */
		if(strstr(strLine, "probeset") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			nProbeNum = atoi(chSep);
		}
		else if(strstr(strLine, "score") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strScorePath, chSep);
		}
		else if(strstr(strLine, "transform") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strTransform, chSep);
		}
		else if(strstr(strLine, "resolution") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			dResolution = atof(chSep);
		}
		else if(strstr(strLine, "locuslink") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strMapPath, chSep);
		}
		else if((strstr(strLine, "output") == strLine) && (strstr(strLine, "path") != NULL) )
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strOutPath, chSep);
		}
		else if((strstr(strLine, "output") == strLine) && (strstr(strLine, "number") != NULL) )
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			nOutputNum = atoi(chSep);
		}
		
	}

	fclose(fpIn);

	/* analysis */
	Expression_PowExpressLocusLink(nProbeNum, strScorePath, strTransform,
		dResolution, strMapPath, strOutPath, nOutputNum);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Expression_PowExpressLocusLink()                                       */
/*  This is the algorithm used by PowerExpress.                            */
/*  Gene selection based on combining multiple probesets by Locuslink ID   */
/*  Return PROC_SUCCESS if success.                                        */
/* ----------------------------------------------------------------------- */ 
int Expression_PowExpressLocusLink(int nProbeNum, char strScorePath[],
								   char strTransform[], double dResolution,
								   char strMapPath[], char strOutPath[],
								   int nOutputNum)
{
	/* define */
	/* score and locuslink map */
	struct DOUBLEMATRIX *pScore;
	struct DOUBLEMATRIX *pWorkScore;
	struct DOUBLEMATRIX *pLocus;
	struct tagString **vProbeName;
	int nLocusNum,nLocusNANum;

	/* sort locus */
	double *vLL;
	long *vLI;
	struct DOUBLEMATRIX *pSortLocus;
	struct LONGMATRIX *pSortIndex;
	struct INTMATRIX *pUniqLocus;

	/* read file */
	FILE *fpIn;
	char strLine[LINE_LENGTH];
	char strLine2[LONG_LINE_LENGTH];
	char strProbe[LINE_LENGTH];
	char strTemp[LINE_LENGTH];
	char *chSep1,*chSep2;

	/* count */
	int ni,nj,nk,nl,nx,ny;
	int nignore;
	double dignoremax;

	/* hierarchical model */
	/* sample mean */
	struct DOUBLEMATRIX *pX;
	/* sample sum of squared residuals */
	struct DOUBLEMATRIX *pS;
	/* sample size */
	struct INTMATRIX *pN;
	/* group map */
	struct INTMATRIX **vLocusMap;

	double *vX,*vS,*vWScore,*vPost;
	int *vN;
	double dTemp,dSTemp;

	/* shrinkage */
	double dPhi,dDf,dTm,dTau,dTauN,dB,dBtemp;
	int nNnum;
	struct DOUBLEMATRIX *pGroupMean;
	struct DOUBLEMATRIX *pGroupVar;
	struct DOUBLEMATRIX *pGroupN;


	/* final score */
	struct DOUBLEMATRIX *pPostMean;
	struct DOUBLEMATRIX *pSortPostMean;
	struct LONGMATRIX *pSortPostIndex;

	/* ----------- */
	/* init        */
	/* ----------- */
	if(nProbeNum <=0 )
	{
		printf("Warning: Expression_PowExpressLocusLink, probeset number = 0!\n");
		return PROC_FAILURE;
	}
	
	/* ----------- */
	/* load score  */
	/* ----------- */
	printf("Load original score...\n");
	fpIn = NULL;
	fpIn = fopen(strScorePath, "rt");
	if(fpIn == NULL)
	{
		printf("Error: Expression_PowExpressLocusLink, cannot open file to load score!\n");
		exit(EXIT_FAILURE);
	}

	pScore = NULL;
	pScore = CreateDoubleMatrix(1, nProbeNum);
	if(pScore == NULL)
	{
		printf("Error: Expression_PowExpressLocusLink, cannot allocate enough memory for loading score!\n");
		exit(EXIT_FAILURE);
	}

	vProbeName = NULL;
	vProbeName = (struct tagString **)calloc(nProbeNum, sizeof(struct tagString *));
	if(vProbeName == NULL)
	{
		printf("Error: Expression_PowExpressLocusLink, cannot allocate enough memory for loading probeset id!\n");
		exit(EXIT_FAILURE);
	}

	fgets(strLine, LINE_LENGTH, fpIn);
	ni = 0;
	while(fgets(strLine, LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(ni >= nProbeNum)
		{
			printf("Error: Expression_PowExpressLocusLink, probeset number not match!\n");
			exit(EXIT_FAILURE);
		}

		sscanf(strLine, "%s %lf", strProbe, (pScore->pMatElement+ni));
		StringAddTail((vProbeName+ni), strProbe);
		ni++;
	}

	fclose(fpIn);

	if(ni != nProbeNum)
	{
		printf("Error: Expression_PowExpressLocusLink, probeset number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* ----------- */
	/* transform   */
	/* ----------- */
	printf("Transform original score...\n");
	pWorkScore = NULL;
	pWorkScore = Expression_GeneRankByLocusLink_ScoreTransform(pScore, dResolution, strTransform);
	if(pWorkScore == NULL)
	{
		printf("Error: Expression_PowExpressLocusLink, cannot allocate enough memory for loading working score!\n");
		exit(EXIT_FAILURE);
	}

	/* -------------- */
	/* load locuslink */
	/* -------------- */
	printf("Create locus link map...\n");
	fpIn = NULL;
	fpIn = fopen(strMapPath, "rt");
	if(fpIn == NULL)
	{
		printf("Error: Expression_PowExpressLocusLink, cannot open file to load locus link!\n");
		exit(EXIT_FAILURE);
	}

	pLocus = NULL;
	pLocus = CreateDoubleMatrix(1, nProbeNum);
	if(pLocus == NULL)
	{
		printf("Error: Expression_PowExpressLocusLink, cannot allocate enough memory for loading locus link!\n");
		exit(EXIT_FAILURE);
	}

	fgets(strLine2, LONG_LINE_LENGTH, fpIn);
	ni = 0;
	nLocusNANum = 0;
	while(fgets(strLine2, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine2);
		StrTrimRight(strLine2);
		if(strLine2[0] == '\0')
			continue;

		if(ni >= nProbeNum)
		{
			printf("Error: Expression_PowExpressLocusLink, probeset number not match!\n");
			exit(EXIT_FAILURE);
		}

		chSep1 = strchr(strLine2, '\t');
		*chSep1 = '\0';
		strcpy(strProbe, strLine2);
		StrTrimRight(strProbe);
		if(strcmp(strProbe, vProbeName[ni]->m_pString) != 0)
		{
			printf("Error: locuslink map not match the probeset id in the score file!\n");
			exit(EXIT_FAILURE);
		}

		chSep1++;
		chSep2 = strchr(chSep1, '\t');
		chSep1 = chSep2+1;
		chSep2 = strchr(chSep1, '\t');
		if(chSep2 != NULL)
			*chSep2 = '\0';
		strcpy(strTemp, chSep1);
		StrTrimRight(strTemp);
		StrTrimLeft(strTemp);
		if(strTemp[0] == '\0')
		{
			nLocusNANum++;
			pLocus->pMatElement[ni] = -nLocusNANum;
		}
		else if(strcmp(strTemp, "---") == 0)
		{
			nLocusNANum++;
			pLocus->pMatElement[ni] = -nLocusNANum;
		}
		else
		{
			pLocus->pMatElement[ni] = atof(strTemp);
		}
		
		ni++;
	}

	fclose(fpIn);

	if(ni != nProbeNum)
	{
		printf("Error: Expression_PowExpressLocusLink, probeset number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* ------------- */
	/* sort & group  */
	/* ------------- */
	/* sort */
	pSortLocus = NULL;
	pSortIndex = NULL;
	if( DMSORTMERGEA_0(pLocus, &pSortLocus, &pSortIndex) == PROC_FAILURE )
	{
		printf("Error: Expression_PowExpressLocusLink, sorting locus link failure!\n");
		exit(EXIT_FAILURE);
	}

	/* get locus link group number */
	nLocusNum = 1;
	vLL = pSortLocus->pMatElement;
	for(ni=1; ni<nProbeNum; ni++)
	{
		if(vLL[ni] != vLL[ni-1])
		{
			nLocusNum++;
		}
	}

	/* prepare mean, sst, df, etc. */
	vLI = pSortIndex->pMatElement;
	vWScore = pWorkScore->pMatElement;

	pUniqLocus = CreateIntMatrix(1, nLocusNum);
	pN = CreateIntMatrix(1, nLocusNum);
	pX = CreateDoubleMatrix(1, nLocusNum);
	pS = CreateDoubleMatrix(1, nLocusNum);
	vN = pN->pMatElement;
	vX = pX->pMatElement;
	vS = pS->pMatElement;
	
	vLocusMap = (struct INTMATRIX **)calloc(nLocusNum, sizeof(struct INTMATRIX *));
	
	nj = 0;
	nl = 1;
	ny = 0;
	for(ni=1; ni<nProbeNum; ni++)
	{
		if(vLL[ni] != vLL[ni-1])
		{
			vLocusMap[ny] = CreateIntMatrix(1, nl);
			dTemp = 0.0;
			nignore = 0;
			for(nk=0; nk<nl; nk++)
			{
				nx = (int)(vLI[ni-nl+nk]);
				if(nk==0)
				{
					dignoremax = vWScore[nx];
					nignore = 0;
				}
				else
				{
					if(vWScore[nx] > dignoremax)
					{
						dignoremax = vWScore[nx];
						nignore = nk;
					}
				}
				dTemp += vWScore[nx];
				vLocusMap[ny]->pMatElement[nk] = nx;
			}

			if(nl == 1)
			{
				vX[ny] = dTemp;
				vS[ny] = 0.0;
				vN[ny] = nl;
			}
			else
			{
				/* dTemp -= dignoremax; */
				/* dTemp /= (double)(nl-1); */
				dTemp /= (double)(nl);
				vX[ny] = dTemp;
		
				dSTemp = 0.0;
				for(nk=0; nk<nl; nk++)
				{
					/* if(nk == nignore)
						continue;
					*/

					nx = (int)(vLI[ni-nl+nk]);
					dSTemp += (vWScore[nx]-dTemp)*(vWScore[nx]-dTemp);
				}

				vS[ny] = dSTemp;
				/* vN[ny] = nl-1; */
				vN[ny] = nl;
			}

			pUniqLocus->pMatElement[ny] = (int)(vLL[ni-1]);
			nj = ni;
			nl = 1;
			ny++;
		}
		else
		{
			nl++;
		}
	}

	vLocusMap[ny] = CreateIntMatrix(1, nl);
	dTemp = 0.0;
	nignore = 0;
	for(nk=0; nk<nl; nk++)
	{
		nx = (int)(vLI[ni-nl+nk]);
		if(nk==0)
		{
			dignoremax = vWScore[nx];
			nignore = 0;
		}
		else
		{
			if(vWScore[nx] > dignoremax)
			{
				dignoremax = vWScore[nx];
				nignore = nk;
			}
		}
		dTemp += vWScore[nx];
		vLocusMap[ny]->pMatElement[nk] = nx;
	}

	if(nl == 1)
	{
		vX[ny] = dTemp;
		vS[ny] = 0.0;
		vN[ny] = nl;
	}
	else
	{
		/* dTemp -= dignoremax;
		dTemp /= (double)(nl-1); */
		dTemp /= (double)(nl);
		vX[ny] = dTemp;

		dSTemp = 0.0;
		for(nk=0; nk<nl; nk++)
		{
			/* if(nk == nignore)
				continue;
			*/

			nx = (int)(vLI[ni-nl+nk]);
			dSTemp += (vWScore[nx]-dTemp)*(vWScore[nx]-dTemp);
		}

		vS[ny] = dSTemp;
		/* vN[ny] = nl-1; */
		vN[ny] = nl;
	}

	pUniqLocus->pMatElement[ny] = (int)(vLL[ni-1]);

	ny++;

	if(ny != nLocusNum)
	{
		printf("Error: Expression_PowExpressLocusLink, locus link group number not match!\n");
		exit(EXIT_FAILURE);
	}
	

	/* -------------- */
	/* shrinkage B    */
	/* -------------- */
	printf("Calculate new score...\n");

	dTm = 0.0;
	dPhi = 0.0;
	dDf = 0.0;
	nNnum = 0;
	for(ni=0; ni<nLocusNum; ni++)
	{
		dTm += vX[ni];
		if(vN[ni] > 1)
		{
			dPhi += vS[ni];
			dDf += (vN[ni]-1.0);
		}
		if(vN[ni] > nNnum)
		{
			nNnum = vN[ni];
		}
	}
	if(dDf <= 0.0)
	{
		printf("Error: Expression_PowExpressLocusLink, no degree of freedom to estimate variance!\n");
		exit(EXIT_FAILURE);
	}
	dPhi /= dDf;
	if(nLocusNum <= 0)
	{
		printf("Error: Expression_PowExpressLocusLink, no gene available!\n");
		exit(EXIT_FAILURE);
	}
	dTm /= (double)(nLocusNum);

	
	/* -------------- */
	/* posterior mean */
	/* -------------- */
	pGroupMean = NULL;
	pGroupMean = CreateDoubleMatrix(1, nNnum);
	if(pGroupMean == NULL)
	{
		printf("Error: Expression_PowExpressLocusLink, cannot complete shrinkage estimation!\n");
		exit(EXIT_FAILURE);
	}
	pGroupVar = NULL;
	pGroupVar = CreateDoubleMatrix(1, nNnum);
	if(pGroupVar == NULL)
	{
		printf("Error: Expression_PowExpressLocusLink, cannot complete shrinkage estimation!\n");
		exit(EXIT_FAILURE);
	}
	pGroupN = NULL;
	pGroupN = CreateDoubleMatrix(1, nNnum);
	if(pGroupN == NULL)
	{
		printf("Error: Expression_PowExpressLocusLink, cannot complete shrinkage estimation!\n");
		exit(EXIT_FAILURE);
	}
	
	for(ni=0; ni<nLocusNum; ni++)
	{
		nk = vN[ni]-1;
		pGroupMean->pMatElement[nk] += vX[ni];
		pGroupN->pMatElement[nk] = pGroupN->pMatElement[nk]+1.0;
	}
	for(nk=0; nk<nNnum; nk++)
	{
		if( pGroupN->pMatElement[nk] > 0.0 )
		{
			pGroupMean->pMatElement[nk] /= pGroupN->pMatElement[nk];
		}
		else
		{
			pGroupMean->pMatElement[nk] = 0.0;
		}
	}
	for(ni=0; ni<nLocusNum; ni++)
	{
		nk = vN[ni]-1;
		pGroupVar->pMatElement[nk] += (vX[ni]-pGroupMean->pMatElement[nk])*(vX[ni]-pGroupMean->pMatElement[nk]);
	}

	dTau = 0.0;
	dTauN = 0.0;
	for(nk=0; nk<nNnum; nk++)
	{
		if( pGroupN->pMatElement[nk] > 1.0 )
		{
			pGroupVar->pMatElement[nk] /= (pGroupN->pMatElement[nk]-1.0);
		}
		else
		{
			pGroupVar->pMatElement[nk] = 0.0;
		}

		if( pGroupN->pMatElement[nk] > (double)POWEXPRESSLOC_MINWIN )
		{
			dTemp = pGroupVar->pMatElement[nk]-dPhi/(double)(nk+1);
			if(  dTemp > 0.0)
			{				
				dTau += dTemp;
				dTauN += 1.0;
			}
		}
	}
	if(dTauN <= 0.0)
	{
		printf("Error: Expression_PowExpressLocusLink, cannot complete shrinkage estimation!\n");
		exit(EXIT_FAILURE);
	}
	dTau /= dTauN;

	pPostMean = NULL;
	pPostMean = CreateDoubleMatrix(1, nLocusNum);
	if(pPostMean == NULL)
	{
		printf("Error: Expression_PowExpressLocusLink, cannot complete shrinkage estimation!\n");
		exit(EXIT_FAILURE);
	}

	if((dTau <= 0.0) || (dPhi <= 0.0))
	{
		printf("Error: Expression_PowExpressLocusLink, cannot complete shrinkage estimation!\n");
		exit(EXIT_FAILURE);
	}

	vPost = pPostMean->pMatElement;
	for(ni=0; ni<nLocusNum; ni++)
	{
		dBtemp = dPhi/vN[ni];
		dB = dBtemp/(dBtemp+dTau);
		vPost[ni] = (1.0-dB)*vX[ni]+dB*dTm;
	}

	DestroyDoubleMatrix(pGroupMean);
	DestroyDoubleMatrix(pGroupVar);
	DestroyDoubleMatrix(pGroupN);
	
	DestroyDoubleMatrix(pX);
	DestroyDoubleMatrix(pS);
	DestroyIntMatrix(pN);

	Expression_GeneRankByLocusLink_ScoreInverseTransform(pPostMean, strTransform);

	/* -------------- */
	/* rank genes     */
	/* -------------- */
	printf("Rank genes by new score...\n");
	pSortPostMean = NULL;
	pSortPostIndex = NULL;
	if( DMSORTMERGEA_0(pPostMean, &pSortPostMean, &pSortPostIndex) == PROC_FAILURE )
	{
		printf("Error: Expression_PowExpressLocusLink, sorting final score failure!\n");
		exit(EXIT_FAILURE);
	}

	/* ----------- */
	/* output      */
	/* ----------- */
	printf("Link annotations...\n");

	/* sortpostindex -> sortlocus -> locusmap -> probeset -> annotation */
	Expression_GeneRankByLocusLink_Output(nLocusNum, pSortPostMean, pSortPostIndex, 
		pUniqLocus, vLocusMap, nProbeNum, vProbeName, pScore, strMapPath, strOutPath, nOutputNum);

	/* ----------- */
	/* clear mem   */
	/* ----------- */
	DestroyIntMatrix(pUniqLocus);
	DestroyDoubleMatrix(pScore);
	DestroyDoubleMatrix(pWorkScore);
	DestroyDoubleMatrix(pLocus);

	for(ni=0; ni<nProbeNum; ni++)
	{
		DeleteString(vProbeName[ni]);
		vProbeName[ni] = NULL;
	}
	free(vProbeName);

	DestroyDoubleMatrix(pSortLocus);
	DestroyLongMatrix(pSortIndex);

	for(ni=0; ni<nLocusNum; ni++)
	{
		DestroyIntMatrix(vLocusMap[ni]);
		vLocusMap[ni] = NULL;
	}
	free(vLocusMap);

	DestroyDoubleMatrix(pPostMean);
	DestroyDoubleMatrix(pSortPostMean);
	DestroyLongMatrix(pSortPostIndex);

	printf("Done!\n");

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Expression_Normalization_Quantile_Main()                               */
/*  Quantile normalization.                                                */
/*  Return PROC_SUCCESS if success.                                        */
/* ----------------------------------------------------------------------- */ 
int Expression_Normalization_Quantile_Main(int nArrayNum, int nProbeNum, 
								char strDataFile[], char strOutFile[],
								int nTakeLog, double dTruncLow)
{
	/* define */
	struct tagString **vChr;
	struct INTMATRIX *pPosition;
	struct DOUBLEMATRIX **vArray;
	struct DOUBLEMATRIX **vSortArray;
	struct LONGMATRIX **vSortIndex;
	FILE *fpData;
	FILE *fpOut;
	int ni,nj,nk;
	char strLongLine[LONG_LINE_LENGTH];
	char *chSep,*chSep2;
	double dValue;

	/* check */
	if( (nArrayNum <= 0) || (nProbeNum <= 0) )
	{
		printf("Warning: Normalization_Quantile_Main, array or probe number <= 0!\n");
		return PROC_SUCCESS;
	}

	/* init */
	pPosition = NULL;
	pPosition = CreateIntMatrix(1, nProbeNum);
	if(pPosition == NULL)
	{
		printf("Error: Normalization_Quantile_Main, cannot create position matrix!\n");
		exit(EXIT_FAILURE);
	}

	vChr = NULL;
	vChr = (struct tagString **)calloc(nProbeNum, sizeof(struct tagString *));
	if(vChr == NULL)
	{
		printf("Error: Normalization_Quantile_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}

	vArray = NULL;
	vArray = (struct DOUBLEMATRIX **)calloc(nArrayNum, sizeof(struct DOUBLEMATRIX *));
	if(vArray == NULL)
	{
		printf("Error: Normalization_Quantile_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nArrayNum; ni++)
	{
		vArray[ni] = CreateDoubleMatrix(1, nProbeNum);
		if(vArray[ni] == NULL)
		{
			printf("Error: Normalization_Quantile_Main, cannot create data space!\n");
			exit(EXIT_FAILURE);
		}
	}

	vSortArray = NULL;
	vSortArray = (struct DOUBLEMATRIX **)calloc(nArrayNum, sizeof(struct DOUBLEMATRIX *));
	if(vSortArray == NULL)
	{
		printf("Error: Normalization_Quantile_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}
	vSortIndex = NULL;
	vSortIndex = (struct LONGMATRIX **)calloc(nArrayNum, sizeof(struct LONGMATRIX *));
	if(vSortIndex == NULL)
	{
		printf("Error: Normalization_Quantile_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}

	/* data */
	fpOut = NULL;
	fpOut = fopen(strOutFile, "wt");
	if(fpOut == NULL)
	{
		printf("Error: Normalization_Quantile_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fpData = NULL;
	fpData = fopen(strDataFile, "rt");
	if(fpData == NULL)
	{
		printf("Error: Normalization_Quantile_Main, cannot open data file!\n");
		exit(EXIT_FAILURE);
	}
	
	fgets(strLongLine, LONG_LINE_LENGTH, fpData);
	fprintf(fpOut, "%s", strLongLine);

	nj = 0;
	while(fgets(strLongLine, LONG_LINE_LENGTH, fpData) != NULL)
	{
		StrTrimLeft(strLongLine);
		StrTrimRight(strLongLine);
		if(strLongLine[0] == '\0')
			continue;

		if(nj >= nProbeNum)
		{
			printf("Error: Normalization_Quantile_Main, probeset number not match!\n");
			exit(EXIT_FAILURE);
		}

		ni = 0;
		
		/* chromosome */
		chSep = strchr(strLongLine, '\t');
		if(chSep != NULL)
			*chSep = '\0';
		StringAddTail((vChr+nj), strLongLine);
		
		/* position */
		chSep++;
		chSep2 = strchr(chSep, '\t');
		if(chSep2 != NULL)
			*chSep2 = '\0';
		pPosition->pMatElement[nj] = atoi(chSep);
		chSep = chSep2;
		

		while(chSep != NULL)
		{
			if( ni >= nArrayNum )
			{
				printf("Error: Normalization_Quantile_Main, array number not match!\n");
				exit(EXIT_FAILURE);
			}

			chSep++;
			chSep2 = strchr(chSep, '\t');
		
			/* middle number */
			if(chSep2 != NULL)
			{
				*chSep2 = '\0';

				if(chSep == chSep2)
				{
					dValue = 0.0;
				}
				else
				{
					dValue = atof(chSep);
				}
			}
			/* last number */
			else
			{
				if(chSep == chSep2)
				{
					dValue = 0.0;
				}
				else
				{
					dValue = atof(chSep);
				}
			}

			if(dValue < dTruncLow)
				dValue = dTruncLow;
			if(nTakeLog == 1)
				dValue = log(dValue)/log(2.0);
		
			DMSETAT(vArray[ni], 0, nj, dValue);
			ni++;

			/* get next */
			chSep = chSep2;
		}

		if(ni!=nArrayNum)
		{
			printf("Error: Normalization_Quantile_Main, array number not match!\n");
			exit(EXIT_FAILURE);
		}

		/* get next line */
		nj++;
	}

	fclose(fpData);

	if(nj != nProbeNum)
	{
		printf("Error: Normalization_Quantile_Main, probe number not match!\n");
		exit(EXIT_FAILURE);
	}


	/* normalization */
	for(ni=0; ni<nArrayNum; ni++)
	{
		DMSORTMERGEA_0(vArray[ni], (vSortArray+ni), (vSortIndex+ni));
	}
	for(nj=0; nj<nProbeNum; nj++)
	{
		dValue = 0.0;
		for(ni=0; ni<nArrayNum; ni++)
		{
			dValue += vSortArray[ni]->pMatElement[nj];
		}
		dValue /= (double)nArrayNum;
		for(ni=0; ni<nArrayNum; ni++)
		{
			vSortArray[ni]->pMatElement[nj] = dValue;
		}
	}

	/* save */
	for(nj=0; nj<nProbeNum; nj++)
	{
		for(ni=0; ni<nArrayNum; ni++)
		{
			nk = (int)(vSortIndex[ni]->pMatElement[nj]);
			vArray[ni]->pMatElement[nk] = vSortArray[ni]->pMatElement[nj];
		}
	}

	for(nj=0; nj<nProbeNum; nj++)
	{
		fprintf(fpOut, "%s\t%d", vChr[nj]->m_pString, pPosition->pMatElement[nj]);
		/* fprintf(fpOut, "%s", vChr[nj]->m_pString); */
		for(ni=0; ni<nArrayNum; ni++)
		{
			fprintf(fpOut, "\t%9.7e", vArray[ni]->pMatElement[nj]);
		}
		fprintf(fpOut, "\n");
	}


	/* release memory */
	fclose(fpOut);
	DestroyIntMatrix(pPosition);
	for(ni=0; ni<nArrayNum; ni++)
	{
		DestroyDoubleMatrix(vArray[ni]);
		vArray[ni] = NULL;
		DestroyDoubleMatrix(vSortArray[ni]);
		vSortArray[ni] = NULL;
		DestroyLongMatrix(vSortIndex[ni]);
		vSortIndex[ni] = NULL;
	}
	free(vArray);
	free(vSortArray);
	free(vSortIndex);

	for(ni=0; ni<nProbeNum; ni++)
	{
		DeleteString(vChr[ni]);
		vChr[ni] = NULL;
	}
	free(vChr);

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  Expression_Normalization_Quantile()                                    */
/*  Quantile normalization.                                                */
/*  Return PROC_SUCCESS if success.                                        */
/* ----------------------------------------------------------------------- */ 
int Expression_Normalization_Quantile(int nArrayNum, int nProbeNum, 
								struct DOUBLEMATRIX **vArray)
{
	/* define */
	struct DOUBLEMATRIX **vSortArray;
	struct LONGMATRIX **vSortIndex;
	int ni,nj,nk;
	double dValue;
	
	/* check */
	if( (nArrayNum <= 0) || (nProbeNum <= 0) )
	{
		printf("Warning: Normalization_Quantile, array or probe number <= 0!\n");
		return PROC_SUCCESS;
	}

	/* init */
	vSortArray = NULL;
	vSortArray = (struct DOUBLEMATRIX **)calloc(nArrayNum, sizeof(struct DOUBLEMATRIX *));
	if(vSortArray == NULL)
	{
		printf("Error: Normalization_Quantile_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}
	vSortIndex = NULL;
	vSortIndex = (struct LONGMATRIX **)calloc(nArrayNum, sizeof(struct LONGMATRIX *));
	if(vSortIndex == NULL)
	{
		printf("Error: Normalization_Quantile_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}

	/* normalization */
	for(ni=0; ni<nArrayNum; ni++)
	{
		DMSORTMERGEA_0(vArray[ni], (vSortArray+ni), (vSortIndex+ni));
	}
	for(nj=0; nj<nProbeNum; nj++)
	{
		dValue = 0.0;
		for(ni=0; ni<nArrayNum; ni++)
		{
			dValue += vSortArray[ni]->pMatElement[nj];
		}
		dValue /= (double)nArrayNum;
		for(ni=0; ni<nArrayNum; ni++)
		{
			vSortArray[ni]->pMatElement[nj] = dValue;
		}
	}

	/* save */
	for(nj=0; nj<nProbeNum; nj++)
	{
		for(ni=0; ni<nArrayNum; ni++)
		{
			nk = (int)(vSortIndex[ni]->pMatElement[nj]);
			vArray[ni]->pMatElement[nk] = vSortArray[ni]->pMatElement[nj];
		}
	}

	/* release memory */
	for(ni=0; ni<nArrayNum; ni++)
	{
		DestroyDoubleMatrix(vSortArray[ni]);
		vSortArray[ni] = NULL;
		DestroyLongMatrix(vSortIndex[ni]);
		vSortIndex[ni] = NULL;
	}
	free(vSortArray);
	free(vSortIndex);

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  Tiling_ProbeSelection_Main()                                           */
/*  Probe selection based on criteria specified in the criteria file       */
/*  Return PROC_SUCCESS if success.                                        */
/* ----------------------------------------------------------------------- */ 
int Tiling_ProbeSelection_Main(char strDataFile[], char strGeneInfoFile[],
							   char strCriteriaFile[], char strOutFile[])
{
	/* ------- */
	/* define  */
	/* ------- */

	/* array */
	int nArrayNum;
	int nProbeNum;
	struct DOUBLEMATRIX **vArray;

	/* class id */
	int nClassNum;
	struct INTMATRIX *pClassSize;
	struct INTMATRIX *pClassID;
	struct INTMATRIX *pPermClassID;
	struct INTMATRIX *pDataClassID;

	/* variance group */
	int nVargroupNum;
	struct INTMATRIX *pVargroupSize;
	struct INTMATRIX **vVargroupMap;

	/* permutation group */
	int nPermgroupNum;
	struct INTMATRIX *pPermgroupSize;
	struct INTMATRIX **vPermgroupMap;

	/* files */
	FILE *fpData;
	FILE *fpCriteria;
	FILE *fpProbe;

	/* strings */
	char strFileName[MED_LINE_LENGTH];
	char strLine[MED_LINE_LENGTH];
	char strLongLine[LONG_LINE_LENGTH];
	char strComparisons[MED_LINE_LENGTH];
	int nPairwiseComp;
	
	/* scores */
	struct DOUBLEMATRIX *pScore;
	struct DOUBLEMATRIX *pSortScore;
	struct DOUBLEMATRIX *pPermScore;
	struct DOUBLEMATRIX *pPermSortScore;
	struct LONGMATRIX *pSortID;
	struct DOUBLEMATRIX *pFDR;

	/* pointers */
	double *vP1,*vP2,*vP3;
	int nx,ny;

	/* parameters */
	double dTruncLow;
	int nTakeLog;
	int nFDRPermNum;
	int nCycPermNum;
	int nOutputNum;

	/* error */
	int nError;

	/* count */
	int ni,nj,nk;
	char *chSep,*chSep2;
	double dValue;
	int nTemp;


	/* --------- */
	/* load data */
	/* --------- */
	printf("Tiling_ProbeSelection_Main:\n");
	printf("loading data...\n");
	
	/* header info */
	nError = 0;
	fpCriteria = NULL;
	fpCriteria = fopen(strCriteriaFile, "rt");
	if(fpCriteria == NULL)
	{
		printf("Error: Tiling_ProbeSelection_Main, cannot open comparison info file!\n");
		exit(EXIT_FAILURE);
	}
	while(fgets(strLine, MED_LINE_LENGTH, fpCriteria) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		/* load basic info */
		if(strcmp(strLine, "[Basic Info]") == 0)
		{
			/* array number */
			fgets(strLine, MED_LINE_LENGTH, fpCriteria);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			chSep = strstr(strLine, "array");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			nArrayNum = atoi(chSep);

			/* probe number */
			fgets(strLine, MED_LINE_LENGTH, fpCriteria);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			chSep = strstr(strLine, "probeset");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			nProbeNum = atoi(chSep);

			/* class number */
			fgets(strLine, MED_LINE_LENGTH, fpCriteria);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			chSep = strstr(strLine, "group");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			nClassNum = atoi(chSep);
		}

		/* load class info */
		else if(strcmp(strLine, "[Group ID]") == 0)
		{
			fgets(strLongLine, LONG_LINE_LENGTH, fpCriteria);
			pDataClassID = Expression_GeneSelection_LoadGroupId(strLongLine);
			if(pDataClassID == NULL)
			{
				nError = 1;
				break;
			}
		}

		/* load comparisons */
		else if(strcmp(strLine, "[Comparisons]") == 0)
		{
			fgets(strComparisons, MED_LINE_LENGTH, fpCriteria);
			StrTrimLeft(strComparisons);
			StrTrimRight(strComparisons);
			if(strComparisons[0] == '\0')
			{
				nError = 1;
				break;
			}

			nPairwiseComp = 0;
			chSep = strpbrk(strComparisons, "<>" );
			while(chSep != NULL)
			{
				nPairwiseComp++;
				chSep = strpbrk((chSep+1), "<>" );
			}

			if(nPairwiseComp > 1)
				nPairwiseComp = 0;
		}

		/* load preprocessing */
		else if(strcmp(strLine, "[Preprocessing Setup]") == 0)
		{
			/* truncate lower bound */
			fgets(strLine, MED_LINE_LENGTH, fpCriteria);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			chSep = strstr(strLine, "lower");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			dTruncLow = atof(chSep);

			/* take log */
			fgets(strLine, MED_LINE_LENGTH, fpCriteria);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			chSep = strstr(strLine, "log2");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			nTakeLog = atoi(chSep);
		}
		
		/* load output setup */
		else if(strcmp(strLine, "[Output Setup]") == 0)
		{
			/* output number */
			fgets(strLine, MED_LINE_LENGTH, fpCriteria);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			chSep = strstr(strLine, "print");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			nOutputNum = atoi(chSep);
		}

		/* load simulation setup */
		else if(strcmp(strLine, "[Simulation Setup]") == 0)
		{
			/* FDR permutation number */
			fgets(strLine, MED_LINE_LENGTH, fpCriteria);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			chSep = strstr(strLine, "FDR");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			nFDRPermNum = atoi(chSep);

			/* within cycle permutation number */
			fgets(strLine, MED_LINE_LENGTH, fpCriteria);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			chSep = strstr(strLine, "cycle");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			nCycPermNum = atoi(chSep);
		}

		/* permutation group */
		else if(strcmp(strLine, "[Permutation Setup]") == 0)
		{
			/* group number */
			fgets(strLine, MED_LINE_LENGTH, fpCriteria);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			chSep = strstr(strLine, "permutation");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			nPermgroupNum = atoi(chSep);

			/* groups */
			pPermgroupSize = NULL;
			pPermgroupSize = CreateIntMatrix(nPermgroupNum, 1);
			if(pPermgroupSize == NULL)
			{
				nError = 1;
				break;
			}
			vPermgroupMap = NULL;
			vPermgroupMap = (struct INTMATRIX **)calloc(nPermgroupNum, sizeof(struct INTMATRIX *));
			if(vPermgroupMap == NULL)
			{
				DestroyIntMatrix(pPermgroupSize);
				nError = 1;
				break;
			}

			for(ni=0; ni<nPermgroupNum; ni++)
			{
				fgets(strLine, MED_LINE_LENGTH, fpCriteria);
				vPermgroupMap[ni] = Expression_GeneSelection_LoadGroupId(strLine);
				if(vPermgroupMap[ni] == NULL)
				{
					nError = 1;
					break;
				}
				IMSETAT(pPermgroupSize, ni, 0, vPermgroupMap[ni]->nWidth);
			}
		}

		/* variance group */
		else if(strcmp(strLine, "[Variance Setup]") == 0)
		{
			/* group number */
			fgets(strLine, MED_LINE_LENGTH, fpCriteria);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			chSep = strstr(strLine, "variance");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			nVargroupNum = atoi(chSep);

			/* groups */
			pVargroupSize = NULL;
			pVargroupSize = CreateIntMatrix(nVargroupNum, 1);
			if(pVargroupSize == NULL)
			{
				nError = 1;
				break;
			}
			vVargroupMap = NULL;
			vVargroupMap = (struct INTMATRIX **)calloc(nVargroupNum, sizeof(struct INTMATRIX *));
			if(vVargroupMap == NULL)
			{
				DestroyIntMatrix(pVargroupSize);
				nError = 1;
				break;
			}

			for(ni=0; ni<nVargroupNum; ni++)
			{
				fgets(strLine, MED_LINE_LENGTH, fpCriteria);
				vVargroupMap[ni] = Expression_GeneSelection_LoadGroupId(strLine);
				if(vVargroupMap[ni] == NULL)
				{
					nError = 1;
					break;
				}
				IMSETAT(pVargroupSize, ni, 0, vVargroupMap[ni]->nWidth);
			}
		}

		/* do nothing */
		else
		{
		}
	}

	fclose(fpCriteria);

	if(nError > 0)
	{
		printf("Error: Tiling_ProbeSelection_Main, comparison info file format wrong!\n");
		exit(EXIT_FAILURE);
	}

	/* prepare space */
	pClassID = NULL;
	pClassID = CreateIntMatrix(1, nArrayNum);
	pClassSize = NULL;
	pClassSize = CreateIntMatrix(1, nClassNum);
	nj = 0;
	for(ni=0; ni<pDataClassID->nWidth; ni++)
	{
		nk = IMGETAT(pDataClassID, 0, ni);
		if( nk > nClassNum )
		{
			printf("Error: Tiling_ProbeSelection_Main, group id out of range!\n");
			exit(EXIT_FAILURE);
		}
		if( nk>0 )
		{
			if(nj >= nArrayNum)
			{
				printf("Error: Tiling_ProbeSelection_Main, group id out of range!\n");
				exit(EXIT_FAILURE);
			}
			IMSETAT(pClassID, 0, nj, nk);
			nTemp = IMGETAT(pClassSize, 0, (nk-1))+1;
			IMSETAT(pClassSize, 0, (nk-1), nTemp);
			nj++;
		}
	}
	if(nj != nArrayNum)
	{
		printf("Error: Tiling_ProbeSelection_Main, array number not match!\n");
		exit(EXIT_FAILURE);
	}

	vArray = NULL;
	vArray = (struct DOUBLEMATRIX **)calloc(nArrayNum, sizeof(struct DOUBLEMATRIX *));
	if(vArray == NULL)
	{
		printf("Error: Tiling_ProbeSelection_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nArrayNum; ni++)
	{
		vArray[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vArray[ni] == NULL)
		{
			printf("Error: Tiling_ProbeSelection_Main, cannot create data space!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* data */
	fpData = NULL;
	fpData = fopen(strDataFile, "rt");
	if(fpData == NULL)
	{
		printf("Error: Tiling_ProbeSelection_Main, cannot open data file!\n");
		exit(EXIT_FAILURE);
	}
	sprintf(strFileName, "%s.pb", strOutFile);
	fpProbe = NULL;
	fpProbe = fopen(strFileName, "wt");
	if(fpProbe == NULL)
	{
		printf("Error: Tiling_ProbeSelection_Main, cannot open probe file!\n");
		exit(EXIT_FAILURE);
	}
	fgets(strLongLine, LONG_LINE_LENGTH, fpData);
	nj = 0;
	while(fgets(strLongLine, LONG_LINE_LENGTH, fpData) != NULL)
	{
		StrTrimLeft(strLongLine);
		StrTrimRight(strLongLine);
		if(strLongLine[0] == '\0')
			continue;

		if(nj >= nProbeNum)
		{
			printf("Error: Tiling_ProbeSelection_Main, probeset number not match!\n");
			exit(EXIT_FAILURE);
		}

		ni = 0;
		nk = 0;
		chSep = strchr(strLongLine, '\t');
		if(chSep != NULL)
			*chSep = '\0';
		fprintf(fpProbe, "%s\t", strLongLine);

		chSep++;
		chSep2 = strchr(chSep, '\t');
		if(chSep2 != NULL)
			*chSep2 = '\0';
		fprintf(fpProbe, "%s\n", chSep);
		chSep = chSep2;

		while(chSep != NULL)
		{
			if( nk >= pDataClassID->nWidth )
			{
				printf("Error: Tiling_ProbeSelection_Main, array number not match!\n");
				exit(EXIT_FAILURE);
			}

			chSep++;
			chSep2 = strchr(chSep, '\t');
		
			if( pDataClassID->pMatElement[nk] > 0 )
			{
				/* middle number */
				if(chSep2 != NULL)
				{
					*chSep2 = '\0';

					if(chSep == chSep2)
					{
						dValue = 0.0;
					}
					else
					{
						dValue = atof(chSep);
					}
				}
				/* last number */
				else
				{
					if(chSep == chSep2)
					{
						dValue = 0.0;
					}
					else
					{
						dValue = atof(chSep);
					}
				}

				if(dValue < dTruncLow)
					dValue = dTruncLow;
				if(nTakeLog == 1)
					dValue = log(dValue)/log(2.0);
			
				DMSETAT(vArray[ni], nj, 0, dValue);
				ni++;
			}

			/* get next */
			nk++;
			chSep = chSep2;
		}

		/* get next line */
		nj++;
	}

	fclose(fpData);
	fclose(fpProbe);
	DestroyIntMatrix(pDataClassID);

	if(nj != nProbeNum)
	{
		printf("Error: Tiling_ProbeSelection_Main, probeset number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* -------------------- */
	/* get original ranking */
	/* -------------------- */
	printf("rank genes...\n");
	pScore = NULL;
	pSortScore = NULL;
	pSortID = NULL;
	if(nPairwiseComp == 1)
	{
		pScore = Tiling_ProbeSelection_tTest(nProbeNum, nArrayNum, vArray,
									nClassNum, pClassSize, pClassID,
									nVargroupNum, pVargroupSize, vVargroupMap, 
									strComparisons);
	}
	else
	{
		pScore = Tiling_ProbeSelection_MonteCarlo(nProbeNum, nArrayNum, vArray,
									nClassNum, pClassSize, pClassID,
									nVargroupNum, pVargroupSize, vVargroupMap, 
									nCycPermNum, strComparisons);
	}

	if( DMSORTMERGEA_0(pScore, &pSortScore, &pSortID) == PROC_FAILURE )
	{
		printf("Error: Tiling_ProbeSelection_Main, sorting score failure!\n");
		exit(EXIT_FAILURE);
	}

	/* ---------------------------- */
	/* permutations to estimate FDR */
	/* ---------------------------- */
	pFDR = NULL;
	if(nFDRPermNum > 0)
	{
		printf("estimate FDR...\n");
		pFDR = CreateDoubleMatrix(1, nProbeNum);

		/* cycles */
		for(ni=0; ni<nFDRPermNum; ni++)
		{
			printf("perm %d...\n", ni);
			/* permute the class label */
			pPermClassID = NULL;
			pPermClassID = Expression_GeneSelection_ClassPerm(nClassNum, pClassID, pClassSize,
				nPermgroupNum, pPermgroupSize, vPermgroupMap);

			/* rerun a single cycle */
			pPermScore = NULL;
			pPermSortScore = NULL;
			if(nPairwiseComp == 1)
			{
				pPermScore = Tiling_ProbeSelection_tTest(nProbeNum, nArrayNum, vArray,
									nClassNum, pClassSize, pPermClassID,
									nVargroupNum, pVargroupSize, vVargroupMap, 
									strComparisons);
			}
			else
			{
				pPermScore = Tiling_ProbeSelection_MonteCarlo(nProbeNum, nArrayNum, vArray,
									nClassNum, pClassSize, pPermClassID,
									nVargroupNum, pVargroupSize, vVargroupMap, 
									nCycPermNum, strComparisons);
			}

			/* sort scores */
			if( DMSORTMERGEA_0(pPermScore, &pPermSortScore, NULL) == PROC_FAILURE )
			{
				printf("Error: Tiling_ProbeSelection_Main, sorting score failure!\n");
				exit(EXIT_FAILURE);
			}

			/* add FDR count */
			vP1 = pFDR->pMatElement;
			vP2 = pSortScore->pMatElement;
			vP3 = pPermSortScore->pMatElement;
			ny = 0;
			for(nx=0; nx<nProbeNum; nx++)
			{
				for(; ny<nProbeNum; ny++)
				{
					if(vP3[ny] > vP2[nx])
						break;
				}

				if(ny >= nProbeNum)
				{
					vP1[nx] += (double)nProbeNum;
				}
				else
				{
					vP1[nx] += (double)ny;
				}
			}

			/* release memory */
			DestroyIntMatrix(pPermClassID);
			DestroyDoubleMatrix(pPermScore);
			DestroyDoubleMatrix(pPermSortScore);
		}

		/* normalize FDR */
		vP1 = pFDR->pMatElement;
		for(ni=0; ni<nProbeNum; ni++)
		{
			vP1[ni] = vP1[ni]/(double)(nFDRPermNum*(ni+1));
		}
		for(ni=nProbeNum-1; ni>0; ni--)
		{
			if(vP1[ni-1] > vP1[ni])
				vP1[ni-1] = vP1[ni];
		}
	}

	/* -------------- */
	/* release memory */
	/* -------------- */

	/* destroy arrays */
	for(ni=0; ni<nArrayNum; ni++)
	{
		DestroyDoubleMatrix(vArray[ni]);
		vArray[ni] = NULL;
	}
	free(vArray);

	/* destroy class ids */
	DestroyIntMatrix(pClassSize);
	DestroyIntMatrix(pClassID);

	/* destroy variance group ids */
	for(ni=0; ni<nVargroupNum; ni++)
	{
		DestroyIntMatrix(vVargroupMap[ni]);
		vVargroupMap[ni] = NULL;
	}
	free(vVargroupMap);
	DestroyIntMatrix(pVargroupSize);

	/* destroy permutation group ids */
	for(ni=0; ni<nPermgroupNum; ni++)
	{
		DestroyIntMatrix(vPermgroupMap[ni]);
		vPermgroupMap[ni] = NULL;
	}
	free(vPermgroupMap);
	DestroyIntMatrix(pPermgroupSize);
	
	/* -------------------------------- */
	/* add gene info to the output file */
	/* -------------------------------- */
	printf("output top genes and link gene information...\n");
	sprintf(strFileName, "%s.pb", strOutFile);
	Tiling_ProbeSelection_Output_WithRandomControl_ForTest(nProbeNum, nOutputNum, 
		pScore, pSortID, pSortScore, pFDR,
		strFileName, strGeneInfoFile, strOutFile);

	/* destroy results */
	DestroyDoubleMatrix(pScore);
	DestroyDoubleMatrix(pFDR);
	DestroyDoubleMatrix(pSortScore);
	DestroyLongMatrix(pSortID);

	/* ------ */
	/* return */
	/* ------ */
	printf("Done!\n");
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  Tiling_ProbeSelection_Output_WithRandomControl()                       */
/*  Output results, a group of random control will be selected.            */
/* ----------------------------------------------------------------------- */ 
int Tiling_ProbeSelection_Output_WithRandomControl(int nProbeNum, int nOutputNum, 
				struct DOUBLEMATRIX *pScore, struct LONGMATRIX *pSortID, 
				struct DOUBLEMATRIX *pSortScore, struct DOUBLEMATRIX *pFDR,
				char strProbeFile[], char strGeneInfoFile[], char strOutFile[])
{
	/* define */
	struct tagString **vProbeName;
	int ni,nj,nk;
	char strOutPath[LINE_LENGTH];
	char strLine[LINE_LENGTH];
	char strLine2[LONG_LINE_LENGTH];
	char strProbe2[LINE_LENGTH];
	char strProbe[LINE_LENGTH];
	FILE *fpIn;
	FILE *fpOut;
	FILE *fpInfo;
	int nIncludeInfo;
	long *vSID;
	double *vSS;
	double *vSF;
	int nChrPos;
	char strChr[LINE_LENGTH];

	struct DOUBLEMATRIX *pRandMat;
	struct DOUBLEMATRIX *pRandSort;
	struct LONGMATRIX *pRandSortId;

	/* temporary */
	/* struct LONGMATRIX *pRank; */


	/* init */
	vProbeName = NULL;
	vProbeName = (struct tagString **)calloc(nProbeNum, sizeof(struct tagString *));

	/* load probe name */
	fpIn = NULL;
	fpIn = fopen(strProbeFile, "rt");
	if(fpIn == NULL)
	{
		printf("Error: Tiling_ProbeSelection_Output, cannot open probeset name file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(ni >= nProbeNum)
		{
			printf("Error: Tiling_ProbeSelection_Output, probeset number not match!\n");
			exit(EXIT_FAILURE);
		}

		StringAddTail((vProbeName+ni), strLine);
		ni++;
	}

	fclose(fpIn);
	if(ni != nProbeNum)
	{
		printf("Error: Tiling_ProbeSelection_Output, probeset number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* output unsorted */
	sprintf(strOutPath, "%s.ori", strOutFile);
	fpOut = NULL;
	fpOut = fopen(strOutPath, "wt");
	if(fpOut == NULL)
	{
		printf("Error: Tiling_ProbeSelection_Output, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	vSS = pScore->pMatElement;
	
	
	fprintf(fpOut, "chromosome\tposition\tscore\n");
	for(ni=0; ni<nProbeNum; ni++)
	{
		if(vProbeName[ni] != NULL)
		{
			fprintf(fpOut, "%s\t", vProbeName[ni]->m_pString);
		}
		else
		{
			fprintf(fpOut, "NA\tNA\t");
		}
		
		fprintf(fpOut, "% 9.7e\n", vSS[ni]);
	}

	/* FOR DEBUG: temporary */
	/* sprintf(strOutPath, "%s.rnk", strOutFile);
	pRank = NULL;
	pRank = CreateLongMatrix(pSortID->nWidth,1);
	vSID = pSortID->pMatElement;
	for(ni=0; ni<nProbeNum; ni++)
	{
		nj = vSID[ni];
		pRank->pMatElement[nj] = ni+1;
	}
	LMSAVE(pRank, strOutPath);
	DestroyLongMatrix(pRank);
	*/


	/* for DEBUG: */
	/* for(ni=0; ni<nProbeNum; ni++)
	{
		if(vProbeName[ni] != NULL)
		{
			sscanf(vProbeName[ni]->m_pString, "%s %d", strChr, &nChrPos);
			fprintf(fpOut, "%d\t", nChrPos);
		}
		else
		{
			fprintf(fpOut, "NA\t");
		}
		
		fprintf(fpOut, "% 9.7e\n", vSS[ni]);
	} */

	fclose(fpOut);

	/* output selected probes */
	sprintf(strOutPath, "%s.txt", strOutFile);
	fpOut = NULL;
	fpOut = fopen(strOutPath, "wt");
	if(fpOut == NULL)
	{
		printf("Error: Tiling_ProbeSelection_Output, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	
	if(strcmp(strGeneInfoFile, "NULL") == 0)
	{
		nIncludeInfo = 0;
	}
	else
	{
		nIncludeInfo = 1;
	}

	vSID = pSortID->pMatElement;
	vSS = pSortScore->pMatElement;
	if(pFDR != NULL)
	{
		vSF = pFDR->pMatElement;
	}

	fprintf(fpOut, "chromosome\tposition\tline_id\trank\tscore\tFDR\tInfo\n");
	for(ni=0; ni<nOutputNum; ni++)
	{
		if(ni >= nProbeNum)
			break;
		nj = pSortID->pMatElement[ni];

		if(vProbeName[nj] != NULL)
		{
			fprintf(fpOut, "%s\t", vProbeName[nj]->m_pString);
		}
		else
		{
			fprintf(fpOut, "NA\tNA\t");
		}
		
		fprintf(fpOut, "%d\t%d\t% 9.7e\t", (nj+1), (ni+1), vSS[ni]);
		
		if(pFDR != NULL)
		{
			fprintf(fpOut, "% 9.7e\t", vSF[ni]);
		}
		else
		{
			fprintf(fpOut, "NA\t");
		}

		if(nIncludeInfo == 1)
		{
			/* TODO: link gene information */
		}
		
		fprintf(fpOut, "\n");
	}

	fclose(fpOut);

	/* output random control probes */
	pRandMat = NULL;
	pRandSort = NULL;
	pRandSortId = NULL;
	pRandMat = DMRANDU(1, nProbeNum);
	DMSORTMERGEA_0(pRandMat, &pRandSort, &pRandSortId);

	sprintf(strOutPath, "%s_ctr.txt", strOutFile);
	fpOut = NULL;
	fpOut = fopen(strOutPath, "wt");
	if(fpOut == NULL)
	{
		printf("Error: Tiling_GeneSelection_Output, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	
	if(strcmp(strGeneInfoFile, "NULL") == 0)
	{
		nIncludeInfo = 1;
	}
	else
	{
		nIncludeInfo = 0;
	}

	vSID = pRandSortId->pMatElement;
	vSS = pSortScore->pMatElement;
	if(pFDR != NULL)
	{
		vSF = pFDR->pMatElement;
	}

	fprintf(fpOut, "chromosome\tposition\tline_id\trank\tscore\tFDR\tInfo\n");
	for(ni=0; ni<nOutputNum; ni++)
	{
		if(ni >= nProbeNum)
			break;
		nk = pRandSortId->pMatElement[ni];
		nj = pSortID->pMatElement[nk];

		if(vProbeName[nj] != NULL)
		{
			fprintf(fpOut, "%s\t", vProbeName[nj]->m_pString);
		}
		else
		{
			fprintf(fpOut, "NA\tNA\t");
		}
		
		fprintf(fpOut, "%d\t%d\t% 9.7e\t", (nj+1), (nk+1), vSS[nk]);
		
		if(pFDR != NULL)
		{
			fprintf(fpOut, "% 9.7e\t", vSF[nk]);
		}
		else
		{
			fprintf(fpOut, "NA\t");
		}

		if(nIncludeInfo == 1)
		{
			/* TODO: link information */
		}
		
		fprintf(fpOut, "\n");
	}

	fclose(fpOut);


	/* release memory */
	DestroyDoubleMatrix(pRandMat);
	DestroyDoubleMatrix(pRandSort);
	DestroyLongMatrix(pRandSortId);

	for(ni=0; ni<nProbeNum; ni++)
	{
		DeleteString(vProbeName[ni]);
		vProbeName[ni] = NULL;
	}
	free(vProbeName);

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  Tiling_ProbeSelection_Output_WithRandomControl_ForTest()               */
/*  Output results, a group of random control will be selected.            */
/* ----------------------------------------------------------------------- */ 
int Tiling_ProbeSelection_Output_WithRandomControl_ForTest(int nProbeNum, int nOutputNum, 
				struct DOUBLEMATRIX *pScore, struct LONGMATRIX *pSortID, 
				struct DOUBLEMATRIX *pSortScore, struct DOUBLEMATRIX *pFDR,
				char strProbeFile[], char strGeneInfoFile[], char strOutFile[])
{
	/* define */
	struct tagString **vProbeName;
	int ni,nj,nk;
	char strOutPath[LINE_LENGTH];
	char strLine[LINE_LENGTH];
	char strLine2[LONG_LINE_LENGTH];
	char strProbe2[LINE_LENGTH];
	char strProbe[LINE_LENGTH];
	FILE *fpIn;
	FILE *fpOut;
	FILE *fpInfo;
	int nIncludeInfo;
	long *vSID;
	double *vSS;
	double *vSF;
	char *chSep;
	int nChrPos;
	char strChr[LINE_LENGTH];

	struct DOUBLEMATRIX *pRandMat;
	struct DOUBLEMATRIX *pRandSort;
	struct LONGMATRIX *pRandSortId;
	struct INTMATRIX *pChrPos;
	int nWinStart;
	int nWinCount;
	double dWinMean;
	double dWinMin;
	double dWinS2;
	int *vPosVec;

	/* temporary */
	struct LONGMATRIX *pRank;
	struct DOUBLEMATRIX *pNewScore;


	/* init */
	vProbeName = NULL;
	vProbeName = (struct tagString **)calloc(nProbeNum, sizeof(struct tagString *));
	pChrPos = NULL;
	pChrPos = CreateIntMatrix(1, nProbeNum);


	/* load probe name */
	fpIn = NULL;
	fpIn = fopen(strProbeFile, "rt");
	if(fpIn == NULL)
	{
		printf("Error: Tiling_ProbeSelection_Output, cannot open probeset name file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(ni >= nProbeNum)
		{
			printf("Error: Tiling_ProbeSelection_Output, probeset number not match!\n");
			exit(EXIT_FAILURE);
		}


		StringAddTail((vProbeName+ni), strLine);

		sscanf(strLine, "%s %d", strChr, &nChrPos);
		pChrPos->pMatElement[ni] = nChrPos;

		ni++;
	}

	fclose(fpIn);
	if(ni != nProbeNum)
	{
		printf("Error: Tiling_ProbeSelection_Output, probeset number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* FOR DEBUG: temporary */
	sprintf(strOutPath, "%s.rnk", strOutFile);
	pRank = NULL;
	pRank = CreateLongMatrix(pSortID->nWidth,1);
	vSID = pSortID->pMatElement;
	for(ni=0; ni<nProbeNum; ni++)
	{
		nj = vSID[ni];
		pRank->pMatElement[nj] = ni+1;
	}
	LMSAVE(pRank, strOutPath);
	DestroyLongMatrix(pRank);
	

	/* output unsorted */
	sprintf(strOutPath, "%s.ori", strOutFile);
	fpOut = NULL;
	fpOut = fopen(strOutPath, "wt");
	if(fpOut == NULL)
	{
		printf("Error: Tiling_ProbeSelection_Output, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	vSS = pScore->pMatElement;
	
	
	/* fprintf(fpOut, "chromosome\tposition\tscore\n");
	for(ni=0; ni<nProbeNum; ni++)
	{
		if(vProbeName[ni] != NULL)
		{
			fprintf(fpOut, "%s\t", vProbeName[ni]->m_pString);
		}
		else
		{
			fprintf(fpOut, "NA\tNA\t");
		}
		
		fprintf(fpOut, "% 9.7e\n", vSS[ni]);
	} */

	/* for DEBUG: */
	for(ni=0; ni<nProbeNum; ni++)
	{
		if(vProbeName[ni] != NULL)
		{
			sscanf(vProbeName[ni]->m_pString, "%s %d", strChr, &nChrPos);
			fprintf(fpOut, "%d\t", nChrPos);
		}
		else
		{
			fprintf(fpOut, "NA\t");
		}
		
		fprintf(fpOut, "% 9.7e\n", vSS[ni]);
	}

	fclose(fpOut);

	/* output windows */
	/* pNewScore = Expression_GeneRankByLocusLink_ScoreTransform(pScore, 0.0001, "RANK"); */
	/* pNewScore = Expression_GeneRankByLocusLink_ScoreTransform(pScore, 0.0001, "LOGIT"); */
	pNewScore = Expression_GeneRankByLocusLink_ScoreTransform(pScore, 0.0001, "Identity");

	/* sprintf(strOutPath, "%s.wnd", strOutFile);
	fpOut = NULL;
	fpOut = fopen(strOutPath, "wt");
	if(fpOut == NULL)
	{
		printf("Error: Tiling_ProbeSelection_Output, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	nWinStart = 0;
	ni = 0;
	vPosVec = pChrPos->pMatElement;
	vSS = pNewScore->pMatElement;
	while(ni<nProbeNum)
	{
		if( vPosVec[ni] - vPosVec[nWinStart] >= 100 )
		{
			nWinCount = ni-nWinStart;
			dWinMean = 0.0;
			dWinMin = 1e6;
			dWinS2 = 0.0;
			for(nj=nWinStart; nj<ni; nj++)
			{
				dWinMean += vSS[nj];
				if(vSS[nj] < dWinMin)
					dWinMin = vSS[nj];
			}
			dWinMean /= nWinCount;
			for(nj=nWinStart; nj<ni; nj++)
			{
				dWinS2 += (vSS[nj]-dWinMean)*(vSS[nj]-dWinMean);
			}
			if(nWinCount <= 1)
				dWinS2 = 0.0;
			else
				dWinS2 /= (double)(nWinCount-1);
			fprintf(fpOut, "%d\t%f\t%f\t%f\n", nWinCount, dWinMin, dWinMean, dWinS2);

			nWinStart = ni;
		}
		ni++;
	}

	nWinCount = ni-nWinStart;
	dWinMean = 0.0;
	dWinMin = 1e6;
	dWinS2 = 0.0;
	for(nj=nWinStart; nj<ni; nj++)
	{
		dWinMean += vSS[nj];
		if(vSS[nj] < dWinMin)
			dWinMin = vSS[nj];
	}
	dWinMean /= nWinCount;
	for(nj=nWinStart; nj<ni; nj++)
	{
		dWinS2 += (vSS[nj]-dWinMean)*(vSS[nj]-dWinMean);
	}
	if(nWinCount <= 1)
		dWinS2 = 0.0;
	else
		dWinS2 /= (double)(nWinCount-1);
	fprintf(fpOut, "%d\t%f\t%f\t%f\n", nWinCount, dWinMin, dWinMean, dWinS2);

	nWinStart = ni;

	fclose(fpOut);

	/* output selected probes */
	/* sprintf(strOutPath, "%s.txt", strOutFile);
	fpOut = NULL;
	fpOut = fopen(strOutPath, "wt");
	if(fpOut == NULL)
	{
		printf("Error: Tiling_ProbeSelection_Output, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	
	if(strcmp(strGeneInfoFile, "NULL") == 0)
	{
		nIncludeInfo = 0;
	}
	else
	{
		nIncludeInfo = 1;
	}

	vSID = pSortID->pMatElement;
	vSS = pSortScore->pMatElement;
	if(pFDR != NULL)
	{
		vSF = pFDR->pMatElement;
	}

	fprintf(fpOut, "chromosome\tposition\tline_id\trank\tscore\tFDR\tInfo\n");
	for(ni=0; ni<nOutputNum; ni++)
	{
		if(ni >= nProbeNum)
			break;
		nj = pSortID->pMatElement[ni];

		if(vProbeName[nj] != NULL)
		{
			fprintf(fpOut, "%s\t", vProbeName[nj]->m_pString);
		}
		else
		{
			fprintf(fpOut, "NA\tNA\t");
		}
		
		fprintf(fpOut, "%d\t%d\t% 9.7e\t", (nj+1), (ni+1), vSS[ni]);
		
		if(pFDR != NULL)
		{
			fprintf(fpOut, "% 9.7e\t", vSF[ni]);
		}
		else
		{
			fprintf(fpOut, "NA\t");
		}

		if(nIncludeInfo == 1)
		{
			/* TODO: link gene information */
	/*	}
		
		fprintf(fpOut, "\n");
	}

	fclose(fpOut);

	/* output random control probes */
	/* pRandMat = NULL;
	pRandSort = NULL;
	pRandSortId = NULL;
	pRandMat = DMRANDU(1, nProbeNum);
	DMSORTMERGEA_0(pRandMat, &pRandSort, &pRandSortId);

	sprintf(strOutPath, "%s_ctr.txt", strOutFile);
	fpOut = NULL;
	fpOut = fopen(strOutPath, "wt");
	if(fpOut == NULL)
	{
		printf("Error: Tiling_GeneSelection_Output, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	
	if(strcmp(strGeneInfoFile, "NULL") == 0)
	{
		nIncludeInfo = 1;
	}
	else
	{
		nIncludeInfo = 0;
	}

	vSID = pRandSortId->pMatElement;
	vSS = pSortScore->pMatElement;
	if(pFDR != NULL)
	{
		vSF = pFDR->pMatElement;
	}

	fprintf(fpOut, "chromosome\tposition\tline_id\trank\tscore\tFDR\tInfo\n");
	for(ni=0; ni<nOutputNum; ni++)
	{
		if(ni >= nProbeNum)
			break;
		nk = pRandSortId->pMatElement[ni];
		nj = pSortID->pMatElement[nk];

		if(vProbeName[nj] != NULL)
		{
			fprintf(fpOut, "%s\t", vProbeName[nj]->m_pString);
		}
		else
		{
			fprintf(fpOut, "NA\tNA\t");
		}
		
		fprintf(fpOut, "%d\t%d\t% 9.7e\t", (nj+1), (nk+1), vSS[nk]);
		
		if(pFDR != NULL)
		{
			fprintf(fpOut, "% 9.7e\t", vSF[nk]);
		}
		else
		{
			fprintf(fpOut, "NA\t");
		}

		if(nIncludeInfo == 1)
		{
			/* TODO: link information */
	/*	}
		
		fprintf(fpOut, "\n");
	}

	fclose(fpOut);
	*/

	/* release memory */
	/* DestroyDoubleMatrix(pRandMat);
	DestroyDoubleMatrix(pRandSort);
	DestroyLongMatrix(pRandSortId);
	*/
	
	for(ni=0; ni<nProbeNum; ni++)
	{
		DeleteString(vProbeName[ni]);
		vProbeName[ni] = NULL;
	}
	free(vProbeName);

	DestroyIntMatrix(pChrPos);
	DestroyDoubleMatrix(pNewScore);

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  Tiling_BindingRegionSelection_HMM_Main()                               */
/*  HMM Forward-Backward algorithm for calculating the posterior prob. of  */
/*  binding for each position.                                             */
/* ----------------------------------------------------------------------- */ 
int Tiling_BindingRegionSelection_HMM_Main(int nProbeNum, char strDataPath[],
				char strTransformType[], double dScoreResolution,
				char strTransitionPath[], char strEmissionPath[], 
				double dGapDist, double dPosteriorCutoff,
				char strOutPath[])
{
	/* define */
	struct DOUBLEMATRIX *pScore;
	struct DOUBLEMATRIX *pNewScore;
	struct DOUBLEMATRIX *pPosition;
	struct DOUBLEMATRIX *pTransition;
	struct DOUBLEMATRIX *pEmission;
	struct DOUBLEMATRIX *pStationary;
	int nStateNum;
	struct DOUBLEMATRIX **vPosterior;
	double dTemp;
	int ni,nj;
	double *pEle;

	/* for loading data */
	char strLine[LINE_LENGTH];
	char strChr[LINE_LENGTH];
	int nPos;
	double dPos,dScore;
	FILE *fpIn;
	FILE *fpOut;

	
	/* check */

	/* init */
	pTransition = NULL;
	pTransition = DMLOAD(strTransitionPath);
	if(pTransition == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_Main, cannot load transition probability!\n");
		exit(EXIT_FAILURE);
	}

	pEmission = NULL;
	pEmission = DMLOAD(strEmissionPath);
	if(pEmission == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_Main, cannot load emission probability!\n");
		exit(EXIT_FAILURE);
	}

	if((pTransition->nHeight+1) != pEmission->nHeight)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_Main, state number not match!\n");
		exit(EXIT_FAILURE);
	}

	nStateNum = pTransition->nHeight;
	if(nStateNum != 2)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_Main, this function only support 2 states now!\n");
		exit(EXIT_FAILURE);
	}

	pStationary = NULL;
	pStationary = CreateDoubleMatrix(1, nStateNum);
	dTemp = DMGETAT(pTransition, 1, 0)+DMGETAT(pTransition, 0, 1);
	pStationary->pMatElement[0] = log(DMGETAT(pTransition, 1, 0)/dTemp);
	pStationary->pMatElement[1] = log(DMGETAT(pTransition, 0, 1)/dTemp);

	pEle = pTransition->pMatElement;
	for(ni=0; ni<pTransition->nHeight; ni++)
	{
		for(nj=0; nj<pTransition->nWidth; nj++)
		{
			*pEle = log(*pEle);
			pEle++;
		}
	}

	pEle = pEmission->pMatElement+pEmission->nWidth;
	for(ni=1; ni<pEmission->nHeight; ni++)
	{
		for(nj=0; nj<pEmission->nWidth; nj++)
		{
			*pEle = log(*pEle);
			pEle++;
		}
	}

	

	/* prepare space */
	pScore = NULL;
	pScore = CreateDoubleMatrix(1, nProbeNum);
	if(pScore == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}

	pPosition = NULL;
	pPosition = CreateDoubleMatrix(1, nProbeNum);
	if(pPosition == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}

	vPosterior = NULL;
	vPosterior = (struct DOUBLEMATRIX **)calloc(nStateNum, sizeof(struct DOUBLEMATRIX *));
	if(vPosterior == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nStateNum; ni++)
	{
		vPosterior[ni] = CreateDoubleMatrix(1, nProbeNum);
		if(vPosterior[ni] == NULL)
		{
			printf("Error: Tiling_BindingRegionSelection_HMM_Main, cannot create data space!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* load data */
	fpIn = NULL;
	fpIn = fopen(strDataPath, "rt");
	if(fpIn == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_Main, cannot open data file!\n");
		exit(EXIT_FAILURE);
	}

	/* fgets(strLine, LINE_LENGTH, fpIn); */

	ni = 0;
	while(fgets(strLine, LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		/* sscanf(strLine, "%s %d %lf", strChr, &nPos, &dScore); */
		sscanf(strLine, "%lf %lf", &dPos, &dScore);
		pScore->pMatElement[ni] = dScore;
		pPosition->pMatElement[ni] = dPos;
		ni++;
	}

	fclose(fpIn);

	if(ni != nProbeNum)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_Main, probe number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* transform score */
	pNewScore = NULL;
	pNewScore = Tiling_BindingRegionSelection_ScoreTransform(pScore, dScoreResolution, strTransformType);
	DestroyDoubleMatrix(pScore);

	/* calculate posterior probability */
	Tiling_BindingRegionSelection_HMM(nProbeNum, nStateNum,
				pStationary, pTransition, pEmission, dGapDist,
				pNewScore, pPosition, vPosterior);

	/* save */
	sprintf(strLine, "%s_pl.txt", strOutPath);
	fpOut = NULL;
	fpOut = fopen(strLine, "wt");
	if(fpOut == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nProbeNum; ni++)
	{
		for(nj=0; nj<nStateNum; nj++)
		{
			fprintf(fpOut, "%9.7e\t", vPosterior[nj]->pMatElement[ni]);
		}
		fprintf(fpOut, "\n");
	}

	fclose(fpOut);

	/* save to *.bed file */
	sprintf(strLine, "%s_rg.bed", strOutPath);
	Tiling_BindingRegionSelection_OutputToBed(nProbeNum, pPosition, vPosterior, 1, 
		dPosteriorCutoff, dGapDist, strLine);


	/* release memory */
	DestroyDoubleMatrix(pTransition);
	DestroyDoubleMatrix(pEmission);
	DestroyDoubleMatrix(pStationary);
	DestroyDoubleMatrix(pNewScore);
	DestroyDoubleMatrix(pPosition);
	for(ni=0; ni<nStateNum; ni++)
	{
		DestroyDoubleMatrix(vPosterior[ni]);
		vPosterior[ni] = NULL;
	}
	free(vPosterior);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Tiling_BindingRegionSelection_HMM()                                    */
/*  HMM Forward-Backward algorithm for calculating the posterior prob. of  */
/*  binding for each position.                                             */
/* ----------------------------------------------------------------------- */ 
int Tiling_BindingRegionSelection_HMM(int nProbeNum, int nStateNum,
				struct DOUBLEMATRIX *pStationary, struct DOUBLEMATRIX *pTransition, 
				struct DOUBLEMATRIX *pEmission, double dGapDist, 
				struct DOUBLEMATRIX *pScore, struct DOUBLEMATRIX *pPosition, 
				struct DOUBLEMATRIX **vPosterior)
{
	/* define */
	struct DOUBLEMATRIX **vForwardSum;
	struct DOUBLEMATRIX **vBackwardSum;
	struct DOUBLEMATRIX *pDP;
	double dSum,dSumF,dSumB;
	double dTemp,dMax;
	double dDist;
	int ni,nj,nk;
	double dInitMax;

	int nEqualLenInterval;
	double dIntS,dIntE,dIntStep;

	/* check */
	if( (nProbeNum <= 0) || (nStateNum <= 0) )
	{
		printf("Error: Tiling_BindingRegionSelection_HMM, nProbeNum/nStateNum <=0!\n");
		exit(EXIT_FAILURE);
	}
	if((pScore == NULL) || (pStationary == NULL) || (pTransition == NULL) 
		|| (pEmission == NULL) || (pPosition == NULL) || (vPosterior == NULL))
	{
		printf("Error: Tiling_BindingRegionSelection_HMM, no input data/parameters!\n");
		exit(EXIT_FAILURE);
	}
	if(pScore->nWidth != nProbeNum)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM, probe number not match!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nStateNum; ni++)
	{
		if(vPosterior[ni]->nWidth != nProbeNum)
		{
			printf("Error: Tiling_BindingRegionSelection_HMM, probe number not match!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* init */
	dInitMax = -DM_ACCESS_VIOLATION;
	
	/* by default, use intervals of equal length for likelihood calculation */
	nEqualLenInterval = 1;
	if(pEmission->nWidth <= 2)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM, emission probability need to be set for more than two intervals!\n");
		exit(EXIT_FAILURE);
	}

	dIntS = DMGETAT(pEmission, 0, 0);
	dIntE = DMGETAT(pEmission, 0, (pEmission->nWidth-2));
	
	dIntStep = (dIntE-dIntS)/(double)(pEmission->nWidth-2);
	if(dIntStep < 0.0)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM, negative interval step length!\n");
		exit(EXIT_FAILURE);
	}


	/* prepare space */
	vForwardSum = NULL;
	vForwardSum = (struct DOUBLEMATRIX **)calloc(nStateNum, sizeof(struct DOUBLEMATRIX *));
	if(vForwardSum == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM, cannot create space for HMM calculation!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nStateNum; ni++)
	{
		vForwardSum[ni] = CreateDoubleMatrix(1, nProbeNum);
		if(vForwardSum[ni] == NULL)
		{
			printf("Error: Tiling_BindingRegionSelection_HMM, cannot create space for HMM calculation!\n");
			exit(EXIT_FAILURE);
		}
	}

	vBackwardSum = NULL;
	vBackwardSum = (struct DOUBLEMATRIX **)calloc(nStateNum, sizeof(struct DOUBLEMATRIX *));
	if(vBackwardSum == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM, cannot create space for HMM calculation!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nStateNum; ni++)
	{
		vBackwardSum[ni] = CreateDoubleMatrix(1, nProbeNum);
		if(vBackwardSum[ni] == NULL)
		{
			printf("Error: Tiling_BindingRegionSelection_HMM, cannot create space for HMM calculation!\n");
			exit(EXIT_FAILURE);
		}
	}

	pDP = NULL;
	pDP = CreateDoubleMatrix(nStateNum, 1);
	if(pDP == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM, cannot create space for HMM calculation!\n");
		exit(EXIT_FAILURE);
	}

	/* forward summation */
	for(ni=0; ni<nStateNum; ni++)
	{
		dTemp = pStationary->pMatElement[ni]+Tiling_BindingRegionSelection_HMM_GetEmission(pEmission, ni, pScore->pMatElement[0], nEqualLenInterval, dIntS, dIntE, dIntStep);
		vForwardSum[ni]->pMatElement[0] = dTemp;
	}

	for(nj=1; nj<nProbeNum; nj++)
	{
		dDist = pPosition->pMatElement[nj] - pPosition->pMatElement[nj-1];
		for(ni=0; ni<nStateNum; ni++)
		{
			dMax = dInitMax;
			for(nk=0; nk<nStateNum; nk++)
			{
				dTemp = vForwardSum[nk]->pMatElement[nj-1] + Tiling_BindingRegionSelection_HMM_GetTransition(pStationary, pTransition, nk, ni, dDist, dGapDist); 
				if(dTemp > dMax)
				{
					dMax = dTemp;
				}
				pDP->pMatElement[nk] = dTemp;
			}

			dSum = 0.0;
			for(nk=0; nk<nStateNum; nk++)
			{
				dTemp = pDP->pMatElement[nk]-dMax;
				dSum += exp(dTemp);
			}
			dSum = log(dSum)+dMax+Tiling_BindingRegionSelection_HMM_GetEmission(pEmission, ni, pScore->pMatElement[nj], nEqualLenInterval, dIntS, dIntE, dIntStep);
			vForwardSum[ni]->pMatElement[nj] = dSum;
		}
	}

	
	dMax = dInitMax;
	for(ni=0; ni<nStateNum; ni++)
	{
		dTemp = vForwardSum[ni]->pMatElement[nj-1]; 
		if(dTemp > dMax)
		{
			dMax = dTemp;
		}
		pDP->pMatElement[ni] = dTemp;
	}

	dSum = 0.0;
	for(ni=0; ni<nStateNum; ni++)
	{
		dTemp = pDP->pMatElement[ni]-dMax;
		dSum += exp(dTemp);
	}
	dSumF = log(dSum)+dMax;


	/* backward summation */
	for(ni=0; ni<nStateNum; ni++)
	{
		vBackwardSum[ni]->pMatElement[nProbeNum-1] = 0.0;
	}
	for(nj=nProbeNum-2; nj>=0; nj--)
	{
		dDist = pPosition->pMatElement[nj+1]-pPosition->pMatElement[nj];
		for(ni=0; ni<nStateNum; ni++)
		{
			dMax = dInitMax;
			for(nk=0; nk<nStateNum; nk++)
			{
				dTemp = vBackwardSum[nk]->pMatElement[nj+1] 
					+ Tiling_BindingRegionSelection_HMM_GetTransition(pStationary, pTransition, ni, nk, dDist, dGapDist) 
					+ Tiling_BindingRegionSelection_HMM_GetEmission(pEmission, nk, pScore->pMatElement[nj+1], nEqualLenInterval, dIntS, dIntE, dIntStep); 
				if(dTemp > dMax)
				{
					dMax = dTemp;
				}
				pDP->pMatElement[nk] = dTemp;
			}

			dSum = 0.0;
			for(nk=0; nk<nStateNum; nk++)
			{
				dTemp = pDP->pMatElement[nk]-dMax;
				dSum += exp(dTemp);
			}
			dSum = log(dSum)+dMax;
			vBackwardSum[ni]->pMatElement[nj] = dSum;
		}
	}

	dMax = dInitMax;
	for(ni=0; ni<nStateNum; ni++)
	{
		dTemp = pStationary->pMatElement[ni] 
			+ Tiling_BindingRegionSelection_HMM_GetEmission(pEmission, ni, pScore->pMatElement[0], nEqualLenInterval, dIntS, dIntE, dIntStep)
			+ vBackwardSum[ni]->pMatElement[0]; 
		if(dTemp > dMax)
		{
			dMax = dTemp;
		}
		pDP->pMatElement[ni] = dTemp;
	}

	dSum = 0.0;
	for(ni=0; ni<nStateNum; ni++)
	{
		dTemp = pDP->pMatElement[ni]-dMax;
		dSum += exp(dTemp);
	}
	dSumB = log(dSum)+dMax;

	if(fabs(dSumF-dSumB) > 1e-3)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM, forward and backbward summation not match!\n");
		exit(EXIT_FAILURE);
	}
	dSumF = (dSumF+dSumB)/2.0;

	/* posterior calculation */
	for(nj=0; nj<nProbeNum; nj++)
	{
		dSum = 0.0;
		for(ni=0; ni<nStateNum; ni++)
		{
			dTemp = vForwardSum[ni]->pMatElement[nj]+vBackwardSum[ni]->pMatElement[nj]-dSumF;
			dTemp = exp(dTemp);
			dSum += dTemp;
			vPosterior[ni]->pMatElement[nj] = dTemp;
		}
		for(ni=0; ni<nStateNum; ni++)
		{
			vPosterior[ni]->pMatElement[nj] /= dSum;
		}
	}

	/* destroy space */
	for(ni=0; ni<nStateNum; ni++)
	{
		DestroyDoubleMatrix(vForwardSum[ni]);
		vForwardSum[ni] = NULL;
		DestroyDoubleMatrix(vBackwardSum[ni]);
		vBackwardSum[ni] = NULL;
	}
	free(vForwardSum);
	free(vBackwardSum);

	DestroyDoubleMatrix(pDP);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Tiling_BindingRegionSelection_ScoreTransform()                         */
/*  transform original scores to working scores                            */
/*  appropriate. Acceptable transformations are logit, rank, identity      */
/*  return the transformed matrix.                                         */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *Tiling_BindingRegionSelection_ScoreTransform(struct DOUBLEMATRIX *pScore, 
										double dResolution, char strTransform[])
{
	/* define */
	struct DOUBLEMATRIX *pNewScore;
	int ni,nj,nk,nl,nx;
	char strTransformType[LINE_LENGTH];
	double *pEle1,*pEle2;
	double dTemp;
	struct LONGMATRIX *pSortIndex;
	struct DOUBLEMATRIX *pSortScore;
	long *vId;
	

	/* init */
	if(pScore == NULL)
		return NULL;

	/* transform */
	strcpy(strTransformType, strTransform);
	StrMakeUpper(strTransformType);

	/* logit */
	if(strcmp(strTransformType, "LOGIT") == 0)
	{
		pNewScore = NULL;
		pNewScore = CreateDoubleMatrix(pScore->nHeight, pScore->nWidth);
		if(pNewScore == NULL)
		{
			printf("Error: Tiling_BindingRegionSelection_ScoreTransform, cannot allocate enough memory for storing working score!\n");
			exit(EXIT_FAILURE);
		}

		pEle1 = pNewScore->pMatElement;
		pEle2 = pScore->pMatElement;
		dTemp = dResolution/2.0;
		for(ni=0; ni<pScore->nHeight; ni++)
		{
			for(nj=0; nj<pScore->nWidth; nj++)
			{
				if(*pEle2 < dTemp)
				{
					*pEle1 = log(dTemp/(1.0-dTemp));
				}
				else if(*pEle2 > 1.0-dTemp)
				{
					*pEle1 = log((1.0-dTemp)/dTemp);
				}
				else
				{
					*pEle1 = log((*pEle2)/(1.0-*pEle2));
				}
				pEle1++;
				pEle2++;
			}
		}
	}
	/* inv logit */
	else if(strcmp(strTransformType, "INVLOGIT") == 0)
	{
		pNewScore = NULL;
		pNewScore = CreateDoubleMatrix(pScore->nHeight, pScore->nWidth);
		if(pNewScore == NULL)
		{
			printf("Error: Tiling_BindingRegionSelection_ScoreTransform, cannot allocate enough memory for storing working score!\n");
			exit(EXIT_FAILURE);
		}

		pEle1 = pNewScore->pMatElement;
		pEle2 = pScore->pMatElement;
		for(ni=0; ni<pScore->nHeight; ni++)
		{
			for(nj=0; nj<pScore->nWidth; nj++)
			{
				dTemp = exp(*pEle2);
				*pEle1 = dTemp/(1.0+dTemp);
				pEle1++;
				pEle2++;
			}
		}
	}
	/* rank */
	else if(strcmp(strTransformType, "RANK") == 0)
	{
		pSortIndex = NULL;
		pSortScore = NULL;
		if( DMSORTMERGEA_0(pScore, &pSortScore, &pSortIndex) == PROC_FAILURE )
		{
			printf("Error: Tiling_BindingRegionSelection_ScoreTransform, sorting score failure!\n");
			exit(EXIT_FAILURE);
		}

		pNewScore = NULL;
		pNewScore = CreateDoubleMatrix(pScore->nHeight, pScore->nWidth);
		if(pNewScore == NULL)
		{
			printf("Error: Tiling_BindingRegionSelection_ScoreTransform, cannot allocate enough memory for storing working score!\n");
			exit(EXIT_FAILURE);
		}

		ni=0;
		pEle1 = pNewScore->pMatElement;
		pEle2 = pSortScore->pMatElement;
		vId = pSortIndex->pMatElement;
		nl = 1;
		for(nj=1; nj<pScore->nWidth; nj++)
		{
			if(pEle2[nj] > pEle2[nj-1])
			{
				for(nx=0; nx<nl; nx++)
				{
					nk = (int)(vId[nj-1-nx]);
					pEle1[nk] = (double)(ni+nj-1)/2.0;
				}
				
				ni = nj;
				nl = 1;
			}
			else
			{
				nl++;
			}
		}
		for(nx=0; nx<nl; nx++)
		{
			nk = (int)(vId[nj-1-nx]);
			pEle1[nk] = (double)(ni+nj-1)/2.0;
		}

		for(nj=0; nj<pNewScore->nWidth; nj++)
		{
			pEle1[nj] = (pEle1[nj]+0.5)/(double)(pNewScore->nWidth);
			pEle1[nj] = log(pEle1[nj]/(1.0-pEle1[nj]));
		}

		DestroyDoubleMatrix(pSortScore);
		DestroyLongMatrix(pSortIndex);
	}
	/* identity */
	else
	{
		pNewScore = NULL;
		pNewScore = DMCLONE(pScore);
	}

	/* return */
	return pNewScore;
}

/* ----------------------------------------------------------------------- */ 
/*  Tiling_BindingRegionSelection_HMM_GetTransition()                      */
/*  Get transition probability for HMM                                     */
/* ----------------------------------------------------------------------- */ 
double Tiling_BindingRegionSelection_HMM_GetTransition(struct DOUBLEMATRIX *pStationary,
					struct DOUBLEMATRIX *pTransition, 
					int nFromS, int nToS, double dDist, double dGapDist)
{
	/* define */
	double dTP;

	/* process */
	if(dDist > dGapDist)
	{
		dTP = DMGETAT(pStationary, 0, nToS);
	}
	else
	{
		dTP = DMGETAT(pTransition, nFromS, nToS);
	}

		
	/* return */
	return dTP;
}

/* ----------------------------------------------------------------------- */ 
/*  Tiling_BindingRegionSelection_HMM_GetEmission()                        */
/*  Get emission probability for HMM                                       */
/* ----------------------------------------------------------------------- */ 
double Tiling_BindingRegionSelection_HMM_GetEmission(struct DOUBLEMATRIX *pEmission, 
					int nState, double dScore, int nEqualLenInterval, 
					double dIntS, double dIntE, double dIntStep)
{
	/* define */
	double dEP;
	int nIndex;

	/* process */
	if(nEqualLenInterval == 1)
	{
		if( dScore <= dIntS)
		{
			nIndex = 0;
		}
		else
		{
			nIndex = (int)((dScore-dIntS)/dIntStep)+1;
		}
		if(nIndex >= pEmission->nWidth)
			nIndex--;
	}
	else
	{
		for(nIndex=0; nIndex<pEmission->nWidth; nIndex++)
		{
			if(dScore <= DMGETAT(pEmission, 0, nIndex))
				break;
		}
		if(nIndex == pEmission->nWidth)
			nIndex--;
	}

	dEP = DMGETAT(pEmission, (nState+1), nIndex);

	/* return */
	return dEP;
}

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
					char strComparisons[])
{
	/* ------ */
	/* define */
	/* ------ */
	struct INTMATRIX *pClassSizeCopy;
	struct INTMATRIX *pDf;
	struct DOUBLEMATRIX *pScore;
	struct DOUBLEMATRIX **vMeans;
	struct DOUBLEMATRIX **vVars;
	struct DOUBLEMATRIX **vSDs;
	
	struct DOUBLEMATRIX *pSigCoef;
	double *vExp,*vSum,*vAve,*vSigma,*vMu,*vSigma2,*vMu2;
	double dExp,dExp2,dVarTemp;
	int nClustId;
	double dDenom,dTemp;
	double dVarmean,dVarsst;
	double dB,dV2,dK,dN;

	char vLogic[MED_LINE_LENGTH];
	char vTNumber[LINE_LENGTH];
	double vGid[MED_LINE_LENGTH];
	
	int nLogicLen;
	int nLeftNum,nRightNum,nTNumLen;

	int ni,nj,nk,nLen;

	int ng1,ng2;

	/* ------------- */
	/* expression    */
	/* ------------- */
	nLogicLen = 0;
	StrTrimLeft(strComparisons);
	StrTrimRight(strComparisons);
	if(strComparisons[0] == '\0')
	{
		printf("Error: Expression_GeneSelection_t-Test, no logic expression formula available!\n");
		exit(EXIT_FAILURE);
	}
	nLen = strlen(strComparisons);
	nj = 0;
	nLeftNum = 0;
	nRightNum = 0;
	nTNumLen = 0;
	for(ni=0; ni<nLen; ni++)
	{
		if( strComparisons[ni] == '(' )
		{
			if(nTNumLen != 0)
			{
				printf("Error: Expression_GeneSelection_t-Test, logic expression wrong!\n");
				exit(EXIT_FAILURE);
			}

			nLeftNum++;
		}
		else if(strComparisons[ni] == ')')
		{
			if(nTNumLen > 0)
			{
				vTNumber[nTNumLen] = '\0';
				vLogic[nj] = 'G';
				vGid[nj] = atof(vTNumber);
				nTNumLen = 0;
				nj++;
			}
			nRightNum++;
		}
		else if( (strComparisons[ni] == '<') || (strComparisons[ni] == '>') 
			|| (strComparisons[ni] == '&') || (strComparisons[ni] == '|') )
		{
			if(nTNumLen > 0)
			{
				vTNumber[nTNumLen] = '\0';
				vLogic[nj] = 'G';
				vGid[nj] = atof(vTNumber);
				nTNumLen = 0;
				nj++;
			}
			vLogic[nj] = strComparisons[ni];
			vGid[nj] = -1.0;
			nj++;	
		}
		else if( (strComparisons[ni] >= '0') && (strComparisons[ni] <= '9') )
		{
			vTNumber[nTNumLen] = strComparisons[ni];
			nTNumLen++;
		}
		else if( (strComparisons[ni] == ' ') || (strComparisons[ni] == '\t') )
		{
			/* do nothing */
		}
		else
		{
			printf("Error: Expression_GeneSelection_t-Test, %c not supported in logic expressions!\n", strComparisons[ni]);
			exit(EXIT_FAILURE);
		}
	}
	if(nTNumLen > 0)
	{
		vTNumber[nTNumLen] = '\0';
		vLogic[nj] = 'G';
		vGid[nj] = atof(vTNumber);
		nTNumLen = 0;
		nj++;
	}
	nLogicLen = nj;

	if(nLeftNum != nRightNum)
	{
		printf("Error: Expression_GeneSelection_t-Test, lack '(' or ')'!\n");
		exit(EXIT_FAILURE);
	}

	if(nLogicLen != 3)
	{
		printf("Error: Expression_GeneSelection_t-Test, logic expression wrong!\n");
		exit(EXIT_FAILURE);
	}
	if((vLogic[0] != 'G') || (vLogic[2] != 'G'))
	{
		printf("Error: Expression_GeneSelection_t-Test, logic expression wrong!\n");
		exit(EXIT_FAILURE);
	}
	if((vLogic[1] != '<') && (vLogic[1] != '>'))
	{
		printf("Error: Expression_GeneSelection_t-Test, logic expression wrong!\n");
		exit(EXIT_FAILURE);
	}

	if(vLogic[1] == '<')
	{
		ng1 = (int)vGid[0]-1;
		ng2 = (int)vGid[2]-1;
	}
	else if(vLogic[1] == '>')
	{
		ng1 = (int)vGid[2]-1;
		ng2 = (int)vGid[0]-1;
	}
	else
	{
		printf("Error: Expression_GeneSelection_t-Test, logic expression wrong!\n");
		exit(EXIT_FAILURE);
	}

	if((ng1<0) || (ng1>=nClassNum) || (ng2<0) || (ng2>=nClassNum))
	{
		printf("Error: Expression_GeneSelection_t-Test, group id out of range!\n");
		exit(EXIT_FAILURE);
	}

	/* for debug purpose*/
	/* for(ni=0; ni<nLogicLen; ni++)
	{
		printf("%c %d\n", vLogic[ni], (int)(vGid[ni]));
	} */
	
	/* ---- */
	/* init */
	/* ---- */
	pScore = NULL;
	pScore = CreateDoubleMatrix(1, nProbeNum);
	pClassSizeCopy = NULL;
	pClassSizeCopy = CreateIntMatrix(1, nClassNum);
	pDf = NULL;
	pDf = CreateIntMatrix(1, nVargroupNum);


	/* ------------- */
	/* create memory */
	/* ------------- */
	vMeans = NULL;
	vMeans = (struct DOUBLEMATRIX **)calloc(nClassNum, sizeof(struct DOUBLEMATRIX *));
	if(vMeans == NULL)
	{
		printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
		exit(EXIT_FAILURE);
	}
	vVars = NULL;
	vVars = (struct DOUBLEMATRIX **)calloc(nClassNum, sizeof(struct DOUBLEMATRIX *));
	if(vVars == NULL)
	{
		printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
		exit(EXIT_FAILURE);
	}
	vSDs = NULL;
	vSDs = (struct DOUBLEMATRIX **)calloc(nVargroupNum, sizeof(struct DOUBLEMATRIX *));
	if(vSDs == NULL)
	{
		printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
		exit(EXIT_FAILURE);
	}
		
	for(ni=0; ni<nClassNum; ni++)
	{
		vMeans[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vMeans[ni] == NULL)
		{
			printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
			exit(EXIT_FAILURE);
		}

		vVars[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vVars[ni] == NULL)
		{
			printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	for(ni=0; ni<nVargroupNum; ni++)
	{
		vSDs[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vSDs[ni] == NULL)
		{
			printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* ---------------- */
	/* get mean and var */
	/* ---------------- */

	/* mean */
	for(ni=0; ni<nArrayNum; ni++)
	{
		nClustId = (pClassID->pMatElement)[ni]-1;
		(pClassSizeCopy->pMatElement)[nClustId] += 1;
		vSum = vMeans[nClustId]->pMatElement;
		vExp = vArray[ni]->pMatElement;
		for(nj=0; nj<nProbeNum; nj++)
		{
			vSum[nj] += vExp[nj];
		}
	}
	
	for(ni=0; ni<nClassNum; ni++)
	{
		if((pClassSizeCopy->pMatElement)[ni] != (pClassSize->pMatElement)[ni])
		{
			printf("Error: Expression_GeneSelection_t-Test, class size not match!\n");
			exit(EXIT_FAILURE);
		}
		if((pClassSizeCopy->pMatElement)[ni] > 0)
		{
			dDenom = (double)((pClassSizeCopy->pMatElement)[ni]);
			vSum = vMeans[ni]->pMatElement;
			for(nj=0; nj<nProbeNum; nj++)
			{
				vSum[nj] /= dDenom;
			}
		}
	}

	/* variance: sum of squares */
	for(ni=0; ni<nArrayNum; ni++)
	{
		nClustId = (pClassID->pMatElement)[ni]-1;
		vSum = vVars[nClustId]->pMatElement;
		vAve = vMeans[nClustId]->pMatElement;
		vExp = vArray[ni]->pMatElement;
		for(nj=0; nj<nProbeNum; nj++)
		{
			dTemp = (vExp[nj]-vAve[nj]);		
			vSum[nj] += dTemp*dTemp;
		}
	}

	/* variance: estimates */
	for(ni=0; ni<nVargroupNum; ni++)
	{
		if(vVargroupMap[ni]->nWidth != pVargroupSize->pMatElement[ni])
		{
			printf("Error: Expression_GeneSelection_t-Test, variance group number not match!\n");
			exit(EXIT_FAILURE);
		}

		for(nj=0; nj<vVargroupMap[ni]->nWidth; nj++)
		{
			nClustId = (vVargroupMap[ni]->pMatElement)[nj]-1;
			if((pClassSize->pMatElement)[nClustId] > 0)
				(pDf->pMatElement)[ni] += (pClassSize->pMatElement)[nClustId]-1;
			vSum = vSDs[ni]->pMatElement;
			vExp = vVars[nClustId]->pMatElement;
			for(nk=0; nk<nProbeNum; nk++)
			{
				vSum[nk] += vExp[nk];
			}
		}
	}
	

	/* ------------- */
	/* shrinking var */
	/* ------------- */
	for(ni=0; ni<nVargroupNum; ni++)
	{
		dVarmean = 0.0;
		dVarsst = 0.0;
		dB = 0.0;

		if((pDf->pMatElement)[ni] > 0)
		{
			dDenom = (double)(pDf->pMatElement)[ni];
			vSum = vSDs[ni]->pMatElement;

			/* xbar */
			for(nj=0; nj<nProbeNum; nj++)
			{
				vSum[nj] /= dDenom;
				dVarmean += vSum[nj];
			}
			dVarmean /= (double)nProbeNum;

			/* sst */
			for(nj=0; nj<nProbeNum; nj++)
			{
				dTemp = vSum[nj]-dVarmean;
				dVarsst += dTemp*dTemp;
			}

			/* shrinkage factor */
			dV2 = 2.0/dDenom;
			dK = (double)nProbeNum;
			dN = 1.0;
			if((dK > 0.0) && (dVarsst > 0.0))
				dB = dV2*(dK-1)/((dN+dV2)*dK)+dN*(dK-1)*dV2*dVarmean*dVarmean/((dN+dV2)*dN*dVarsst);
				/* dB = 0.0; */
			else
				dB = 0.0;

			/* shrink variance */
			for(nj=0; nj<nProbeNum; nj++)
			{
				vSum[nj] = (1-dB)*vSum[nj]+dB*dVarmean;
				/* vSum[nj] += 1e-16; */
			}
		}
	}

	for(ni=0; ni<nClassNum; ni++)
	{
		DestroyDoubleMatrix(vVars[ni]);
		vVars[ni] = NULL;
	}

	for(ni=0; ni<nVargroupNum; ni++)
	{
		for(nj=0; nj<vVargroupMap[ni]->nWidth; nj++)
		{
			nClustId = (vVargroupMap[ni]->pMatElement)[nj]-1;
			vVars[nClustId] = vSDs[ni];
		}
	}

	/* ----------- */
	/* simulation  */
	/* ----------- */
	pSigCoef = NULL;
	pSigCoef = CreateDoubleMatrix(1, nClassNum);
	vExp = pSigCoef->pMatElement;
	for(ni=0; ni<nClassNum; ni++)
	{
		if((pClassSize->pMatElement)[ni] > 0)
		{
			dDenom = (double)((pClassSize->pMatElement)[ni]);
			vExp[ni] = 1.0/dDenom;
		}
		else
		{
			vExp[ni] = 0.0;
		}
	}


	vMu = vMeans[ng1]->pMatElement;
	vSigma = vVars[ng1]->pMatElement;
	vMu2 = vMeans[ng2]->pMatElement;
	vSigma2 = vVars[ng2]->pMatElement;
	dExp = pSigCoef->pMatElement[ng1];
	dExp2 = pSigCoef->pMatElement[ng2];
	vSum = pScore->pMatElement;

	/* probeset by probeset */
	for(nj=0; nj<nProbeNum; nj++)
	{
		if((vSigma[nj] <= 0.0) || (vSigma2[nj] <= 0.0))
		{
			printf("Warning: Expression_GeneSelection_t-Test, variance=0, may not have enough sample to estimate variance!\n");
		}
		dVarTemp = sqrt(vSigma[nj]*dExp+vSigma2[nj]*dExp2);
		if(dVarTemp > 0.0)
			vSum[nj] = (vMu[nj]-vMu2[nj])/dVarTemp;
		else
			vSum[nj] = 0.0;
	}
	
	DestroyDoubleMatrix(pSigCoef);

	/* -------------- */
	/* release memory */
	/* -------------- */
	for(ni=0; ni<nClassNum; ni++)
	{
		DestroyDoubleMatrix(vMeans[ni]);
		vMeans[ni] = NULL;
	}

	for(ni=0; ni<nVargroupNum; ni++)
	{
		DestroyDoubleMatrix(vSDs[ni]);
		vSDs[ni] = NULL;		
	}

	free(vMeans);
	free(vVars);
	free(vSDs);

	DestroyIntMatrix(pClassSizeCopy);
	DestroyIntMatrix(pDf);

	/* ------ */
	/* return */
	/* ------ */
	return pScore;
}


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
					int nCycPermNum, char strComparisons[])
{
	/* ------ */
	/* define */
	/* ------ */
	struct INTMATRIX *pClassSizeCopy;
	struct INTMATRIX *pDf;
	struct DOUBLEMATRIX *pScore;
	struct DOUBLEMATRIX **vMeans;
	struct DOUBLEMATRIX **vVars;
	struct DOUBLEMATRIX **vSDs;
	/* struct DOUBLEMATRIX *pMeanDraws; */
	struct DOUBLEMATRIX **vMeanDraws;
	struct DOUBLEMATRIX *pSigCoef;
	double *vExp,*vSum,*vAve,*vSigma,*vMu;
	unsigned char *vEval;
	int nClustId;
	double dDenom,dTemp;
	double dVarmean,dVarsst;
	double dB,dV2,dK,dN;

	char vLogic[MED_LINE_LENGTH];
	char vTNumber[LINE_LENGTH];
	double vGid[MED_LINE_LENGTH];
	struct DOUBLEMATRIX *vVid[MED_LINE_LENGTH];
	struct BYTEMATRIX *vLid[MED_LINE_LENGTH];
	int nLogicLen;
	int nLeftNum,nRightNum,nTNumLen;
	int nTrue;
	struct BYTEMATRIX *pTrueVec;
	struct DOUBLEMATRIX *pTrueVal;

	int ni,nj,nk,nLen;

	/* ------------- */
	/* expression    */
	/* ------------- */
	nLogicLen = 0;
	StrTrimLeft(strComparisons);
	StrTrimRight(strComparisons);
	if(strComparisons[0] == '\0')
	{
		printf("Error: Expression_GeneSelection_MonteCarlo, no logic expression formula available!\n");
		exit(EXIT_FAILURE);
	}
	nLen = strlen(strComparisons);
	nj = 0;
	nLeftNum = 0;
	nRightNum = 0;
	nTNumLen = 0;
	for(ni=0; ni<nLen; ni++)
	{
		if( strComparisons[ni] == '(' )
		{
			if(nTNumLen != 0)
			{
				printf("Error: Expression_GeneSelection_MonteCarlo, logic expression wrong!\n");
				exit(EXIT_FAILURE);
			}

			vLogic[nj] = strComparisons[ni];
			vGid[nj] = -1.0;
			nj++;
			nLeftNum++;
		}
		else if(strComparisons[ni] == ')')
		{
			if(nTNumLen > 0)
			{
				vTNumber[nTNumLen] = '\0';
				vLogic[nj] = 'G';
				vGid[nj] = atof(vTNumber);
				nTNumLen = 0;
				nj++;
			}
			vLogic[nj] = strComparisons[ni];
			vGid[nj] = -1.0;
			nj++;
			nRightNum++;
		}
		else if( (strComparisons[ni] == '<') || (strComparisons[ni] == '>') 
			|| (strComparisons[ni] == '&') || (strComparisons[ni] == '|') )
		{
			if(nTNumLen > 0)
			{
				vTNumber[nTNumLen] = '\0';
				vLogic[nj] = 'G';
				vGid[nj] = atof(vTNumber);
				nTNumLen = 0;
				nj++;
			}
			vLogic[nj] = strComparisons[ni];
			vGid[nj] = -1.0;
			nj++;
			
		}
		else if( (strComparisons[ni] >= '0') && (strComparisons[ni] <= '9') )
		{
			vTNumber[nTNumLen] = strComparisons[ni];
			nTNumLen++;
		}
		else if( (strComparisons[ni] == ' ') || (strComparisons[ni] == '\t') )
		{
			/* do nothing */
		}
		else
		{
			printf("Error: Expression_GeneSelection_MonteCarlo, %c not supported in logic expressions!\n", strComparisons[ni]);
			exit(EXIT_FAILURE);
		}
	}
	if(nTNumLen > 0)
	{
		vTNumber[nTNumLen] = '\0';
		vLogic[nj] = 'G';
		vGid[nj] = atof(vTNumber);
		nTNumLen = 0;
		nj++;
	}
	nLogicLen = nj;

	if(nLeftNum != nRightNum)
	{
		printf("Error: Expression_GeneSelection_MonteCarlo, lack '(' or ')'!\n");
		exit(EXIT_FAILURE);
	}

	/* for debug purpose*/
	/* for(ni=0; ni<nLogicLen; ni++)
	{
		printf("%c %d\n", vLogic[ni], (int)(vGid[ni]));
	} */
	
	/* ---- */
	/* init */
	/* ---- */
	pScore = NULL;
	pScore = CreateDoubleMatrix(1, nProbeNum);
	pClassSizeCopy = NULL;
	pClassSizeCopy = CreateIntMatrix(1, nClassNum);
	pDf = NULL;
	pDf = CreateIntMatrix(1, nVargroupNum);


	/* ------------- */
	/* create memory */
	/* ------------- */
	vMeans = NULL;
	vMeans = (struct DOUBLEMATRIX **)calloc(nClassNum, sizeof(struct DOUBLEMATRIX *));
	if(vMeans == NULL)
	{
		printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
		exit(EXIT_FAILURE);
	}
	vVars = NULL;
	vVars = (struct DOUBLEMATRIX **)calloc(nClassNum, sizeof(struct DOUBLEMATRIX *));
	if(vVars == NULL)
	{
		printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
		exit(EXIT_FAILURE);
	}
	vSDs = NULL;
	vSDs = (struct DOUBLEMATRIX **)calloc(nVargroupNum, sizeof(struct DOUBLEMATRIX *));
	if(vSDs == NULL)
	{
		printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
		exit(EXIT_FAILURE);
	}
	vMeanDraws = NULL;
	vMeanDraws = (struct DOUBLEMATRIX **)calloc(nClassNum, sizeof(struct DOUBLEMATRIX *));
	if(vMeanDraws == NULL)
	{
		printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
		exit(EXIT_FAILURE);
	}
	
	for(ni=0; ni<nClassNum; ni++)
	{
		vMeans[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vMeans[ni] == NULL)
		{
			printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
			exit(EXIT_FAILURE);
		}

		vVars[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vVars[ni] == NULL)
		{
			printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
			exit(EXIT_FAILURE);
		}

		vMeanDraws[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vMeanDraws[ni] == NULL)
		{
			printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	for(ni=0; ni<nVargroupNum; ni++)
	{
		vSDs[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vSDs[ni] == NULL)
		{
			printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation!\n");
			exit(EXIT_FAILURE);
		}
	}

	/*pMeanDraws = NULL;
	pMeanDraws = CreateDoubleMatrix(1,nClassNum);
	*/

	/* ---------------- */
	/* get mean and var */
	/* ---------------- */

	/* mean */
	for(ni=0; ni<nArrayNum; ni++)
	{
		nClustId = (pClassID->pMatElement)[ni]-1;
		(pClassSizeCopy->pMatElement)[nClustId] += 1;
		vSum = vMeans[nClustId]->pMatElement;
		vExp = vArray[ni]->pMatElement;
		for(nj=0; nj<nProbeNum; nj++)
		{
			vSum[nj] += vExp[nj];
		}
	}
	
	for(ni=0; ni<nClassNum; ni++)
	{
		if((pClassSizeCopy->pMatElement)[ni] != (pClassSize->pMatElement)[ni])
		{
			printf("Error: Expression_GeneSelection_MonteCarlo, class size not match!\n");
			exit(EXIT_FAILURE);
		}
		if((pClassSizeCopy->pMatElement)[ni] > 0)
		{
			dDenom = (double)((pClassSizeCopy->pMatElement)[ni]);
			vSum = vMeans[ni]->pMatElement;
			for(nj=0; nj<nProbeNum; nj++)
			{
				vSum[nj] /= dDenom;
			}
		}
	}

	/* variance: sum of squares */
	for(ni=0; ni<nArrayNum; ni++)
	{
		nClustId = (pClassID->pMatElement)[ni]-1;
		vSum = vVars[nClustId]->pMatElement;
		vAve = vMeans[nClustId]->pMatElement;
		vExp = vArray[ni]->pMatElement;
		for(nj=0; nj<nProbeNum; nj++)
		{
			dTemp = (vExp[nj]-vAve[nj]);		
			vSum[nj] += dTemp*dTemp;
		}
	}

	/* variance: estimates */
	for(ni=0; ni<nVargroupNum; ni++)
	{
		if(vVargroupMap[ni]->nWidth != pVargroupSize->pMatElement[ni])
		{
			printf("Error: Expression_GeneSelection_MonteCarlo, variance group number not match!\n");
			exit(EXIT_FAILURE);
		}

		for(nj=0; nj<vVargroupMap[ni]->nWidth; nj++)
		{
			nClustId = (vVargroupMap[ni]->pMatElement)[nj]-1;
			if((pClassSize->pMatElement)[nClustId] > 0)
				(pDf->pMatElement)[ni] += (pClassSize->pMatElement)[nClustId]-1;
			vSum = vSDs[ni]->pMatElement;
			vExp = vVars[nClustId]->pMatElement;
			for(nk=0; nk<nProbeNum; nk++)
			{
				vSum[nk] += vExp[nk];
			}
		}
	}
	

	/* ------------- */
	/* shrinking var */
	/* ------------- */
	for(ni=0; ni<nVargroupNum; ni++)
	{
		dVarmean = 0.0;
		dVarsst = 0.0;
		dB = 0.0;

		if((pDf->pMatElement)[ni] > 0)
		{
			dDenom = (double)(pDf->pMatElement)[ni];
			vSum = vSDs[ni]->pMatElement;

			/* xbar */
			for(nj=0; nj<nProbeNum; nj++)
			{
				vSum[nj] /= dDenom;
				dVarmean += vSum[nj];
			}
			dVarmean /= (double)nProbeNum;

			/* sst */
			for(nj=0; nj<nProbeNum; nj++)
			{
				dTemp = vSum[nj]-dVarmean;
				dVarsst += dTemp*dTemp;
			}

			/* shrinkage factor */
			dV2 = 2.0/dDenom;
			dK = (double)nProbeNum;
			dN = 1.0;
			if((dK > 0.0) && (dVarsst > 0.0))
				dB = dV2*(dK-1)/((dN+dV2)*dK)+dN*(dK-1)*dV2*dVarmean*dVarmean/((dN+dV2)*dN*dVarsst);
				/* dB = 0.0; */
			else
				dB = 0.0;

			/* shrink variance */
			for(nj=0; nj<nProbeNum; nj++)
			{
				vSum[nj] = (1-dB)*vSum[nj]+dB*dVarmean;
				vSum[nj] = sqrt(vSum[nj]);
 				/* vSum[nj] = sqrt(vSum[nj]) + 1e-16; */
			}
		}
	}

	for(ni=0; ni<nClassNum; ni++)
	{
		DestroyDoubleMatrix(vVars[ni]);
		vVars[ni] = NULL;
	}

	for(ni=0; ni<nVargroupNum; ni++)
	{
		for(nj=0; nj<vVargroupMap[ni]->nWidth; nj++)
		{
			nClustId = (vVargroupMap[ni]->pMatElement)[nj]-1;
			vVars[nClustId] = vSDs[ni];
		}
	}

	/* ----------- */
	/* simulation  */
	/* ----------- */
	pSigCoef = NULL;
	pSigCoef = CreateDoubleMatrix(1, nClassNum);
	vExp = pSigCoef->pMatElement;
	for(ni=0; ni<nClassNum; ni++)
	{
		if((pClassSize->pMatElement)[ni] > 0)
		{
			dDenom = (double)((pClassSize->pMatElement)[ni]);
			vExp[ni] = sqrt(1.0/dDenom);
		}
		else
		{
			vExp[ni] = 0.0;
		}
	}

	for(ni=0; ni<nLogicLen; ni++)
	{
		if(vGid[ni] > 0.0)
		{
			nk = (int)(vGid[ni]-1.0);
			vVid[ni] = vMeanDraws[nk];
		}
		else
		{
			vVid[ni] = NULL;
		}
		vLid[ni] = NULL;
	}

	for(ni=0; ni<nCycPermNum; ni++)
	{
		if(ni%100 == 0)
		{
			printf("iter %d...\n", ni);
		}

		vExp = pSigCoef->pMatElement;
		
		/* probeset by probeset */
		for(nk=0; nk<nClassNum; nk++)
		{
			vAve = vMeanDraws[nk]->pMatElement;
			vMu = vMeans[nk]->pMatElement;
			vSigma = vVars[nk]->pMatElement;
			dTemp = normrnd(0.0, 1.0);
			for(nj=0; nj<nProbeNum; nj++)
			{
				/* simulate */
				vAve[nj] = vMu[nj]+vExp[nk]*vSigma[nj]*dTemp;
			}
		}

		

		/* evaluate */
		pTrueVec = NULL;
		pTrueVal = NULL;
		nTrue = Expression_GeneSelection_EvaluateVec(nProbeNum, vLogic, vVid, vLid, nLogicLen, 0, &pTrueVec, &pTrueVal);
		if(nTrue != 4)
		{
			printf("Error: logic evaluation wrong!\n");
			exit(EXIT_FAILURE);
		}

		/* nTrue = Expression_GeneSelection_Evaluate(pMeanDraws, vLogic, vGid, nLogicLen, 0); */

		/* add */
		vSum = pScore->pMatElement;
		/* vSum[nj] += (double)(nTrue); */
		
		vEval = pTrueVec->pMatElement;
		for(nj=0; nj<nProbeNum; nj++)
		{
			vSum[nj] += (double)(vEval[nj]);
		}
		DestroyByteMatrix(pTrueVec);
	}
	DestroyDoubleMatrix(pSigCoef);

	/* normalize */
	dDenom = (double)nCycPermNum;
	vSum = pScore->pMatElement;
	for(nj=0; nj<nProbeNum; nj++)
	{
		vSum[nj] = 1.0-vSum[nj]/dDenom;
	}

	/* -------------- */
	/* release memory */
	/* -------------- */
	for(ni=0; ni<nClassNum; ni++)
	{
		DestroyDoubleMatrix(vMeans[ni]);
		vMeans[ni] = NULL;
		
		DestroyDoubleMatrix(vMeanDraws[ni]);
		vMeanDraws[ni] = NULL;
	}

	for(ni=0; ni<nVargroupNum; ni++)
	{
		DestroyDoubleMatrix(vSDs[ni]);
		vSDs[ni] = NULL;		
	}

	free(vMeans);
	free(vVars);
	free(vSDs);
	free(vMeanDraws);

	/* DestroyDoubleMatrix(pMeanDraws); */
	DestroyIntMatrix(pClassSizeCopy);
	DestroyIntMatrix(pDf);

	/* ------ */
	/* return */
	/* ------ */
	return pScore;
}

/* ----------------------------------------------------------------------- */ 
/*  Tiling_BindingRegionSelection_OutputToBed()                            */
/*  Output binding region to a bed file.                                   */
/* ----------------------------------------------------------------------- */ 
int Tiling_BindingRegionSelection_OutputToBed(int nProbeNum, 
					struct DOUBLEMATRIX *pPosition, struct DOUBLEMATRIX **vPosterior, 
					int nStateId, double dPosteriorCutoff, double dGapDist, 
					char strOutPath[])
{
	/* define */
	FILE *fpOut;
	double dDist;
	int ni;
	int nStart,nEnd;
	int nP1,nP2;
	double dMeanPost;
	int nProbeCount;

	/* init */
	if(pPosition->nWidth != nProbeNum)
	{
		printf("Error: probe number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* output */
	fpOut = NULL;
	fpOut = fopen(strOutPath, "wt");
	if(fpOut == NULL)
	{
		printf("Error: cannot open *.bed file for output!\n");
	}

	dDist = 0.0;
	nProbeCount = 0;
	dMeanPost = 0.0;
	nStart = -1;
	nEnd = -1;
	nP1 = -1;
	nP2 = -1;
	
	if(vPosterior[nStateId]->pMatElement[0] > dPosteriorCutoff)
	{
		nStart = (int)(pPosition->pMatElement[0])-12;
		nEnd = (int)(pPosition->pMatElement[0])+12;
		nP1 = 0;
		nP2 = 0;
		nProbeCount = 1;
		dMeanPost = vPosterior[nStateId]->pMatElement[0];
	}

	for(ni=1; ni<nProbeNum; ni++)
	{
		dDist = pPosition->pMatElement[ni]-pPosition->pMatElement[ni-1];
		if(dDist > dGapDist)
		{
			if(nProbeCount > 0)
			{
				dMeanPost /= (double)nProbeCount;
				fprintf(fpOut, "%d\t%d\t%d\t%d\t%d\t%f\n", nStart, nEnd, nP1, nP2, nProbeCount, dMeanPost);
				nProbeCount = 0;
				dMeanPost = 0.0;
				nStart = -1;
				nEnd = -1;
				nP1 = -1;
				nP2 = -1;
			}
			else
			{
			}

			if(vPosterior[nStateId]->pMatElement[ni] > dPosteriorCutoff)
			{
				nStart = (int)(pPosition->pMatElement[ni])-12;
				nEnd = (int)(pPosition->pMatElement[ni])+12;
				nP1 = ni;
				nP2 = ni;
				nProbeCount = 1;
				dMeanPost = vPosterior[nStateId]->pMatElement[ni];
			}
		}
		else
		{
			if(vPosterior[nStateId]->pMatElement[ni] > dPosteriorCutoff)
			{
				if(nProbeCount > 0)
				{
					nEnd = (int)(pPosition->pMatElement[ni])+12;
					nP2 = ni;
					nProbeCount += 1;
					dMeanPost += vPosterior[nStateId]->pMatElement[ni];
				}
				else
				{
					nStart = (int)(pPosition->pMatElement[ni])-12;
					nEnd = (int)(pPosition->pMatElement[ni])+12;
					nP1 = ni;
					nP2 = ni;
					nProbeCount = 1;
					dMeanPost = vPosterior[nStateId]->pMatElement[ni];
				}
			}
			else
			{
				if(nProbeCount > 0)
				{
					dMeanPost /= (double)nProbeCount;
					fprintf(fpOut, "%d\t%d\t%d\t%d\t%d\t%f\n", nStart, nEnd, nP1, nP2, nProbeCount, dMeanPost);
					nProbeCount = 0;
					dMeanPost = 0.0;
					nStart = -1;
					nEnd = -1;
					nP1 = -1;
					nP2 = -1;
				}
			}
		}
	}

	if(nProbeCount > 0)
	{
		dMeanPost /= (double)nProbeCount;
		fprintf(fpOut, "%d\t%d\t%d\t%d\t%d\t%f\n", nStart, nEnd, nP1, nP2, nProbeCount, dMeanPost);
	}

	/* close file */
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}


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
				char strOutPath[])
{
	/* define */
	struct DOUBLEMATRIX *pScore;
	struct DOUBLEMATRIX *pNewScore;
	struct DOUBLEMATRIX *pPosition;
	struct DOUBLEMATRIX *pTransition;
	struct DOUBLEMATRIX *pEmission;
	struct DOUBLEMATRIX *pStationary;
	struct DOUBLEMATRIX *pStationaryTrue;
	struct DOUBLEMATRIX *pTransitionNew;
	struct DOUBLEMATRIX *pEmissionNew;
	struct DOUBLEMATRIX *pStationaryNew;

	int nStateNum;
	struct DOUBLEMATRIX **vPosterior;
	double dTemp;
	int ni,nj;
	double *pEle;

	/* for loading data */
	char strLine[LINE_LENGTH];
	char strChr[LINE_LENGTH];
	int nPos;
	double dScore;
	FILE *fpIn;
	FILE *fpOut;

	
	/* check */

	/* init */
	pTransition = NULL;
	pTransition = DMLOAD(strTransitionPath);
	if(pTransition == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_Main, cannot load transition probability!\n");
		exit(EXIT_FAILURE);
	}

	pEmission = NULL;
	pEmission = DMLOAD(strEmissionPath);
	if(pEmission == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_Main, cannot load emission probability!\n");
		exit(EXIT_FAILURE);
	}

	if((pTransition->nHeight+1) != pEmission->nHeight)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_Main, state number not match!\n");
		exit(EXIT_FAILURE);
	}

	nStateNum = pTransition->nHeight;
	if(nStateNum != 2)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_Main, this function only support 2 states now!\n");
		exit(EXIT_FAILURE);
	}

	pStationary = NULL;
	pStationary = CreateDoubleMatrix(1, nStateNum);
	dTemp = DMGETAT(pTransition, 1, 0)+DMGETAT(pTransition, 0, 1);
	pStationary->pMatElement[0] = log(DMGETAT(pTransition, 1, 0)/dTemp);
	pStationary->pMatElement[1] = log(DMGETAT(pTransition, 0, 1)/dTemp);

	pStationaryTrue = NULL;
	pStationaryTrue = DMCLONE(pStationary);

	pEle = pTransition->pMatElement;
	for(ni=0; ni<pTransition->nHeight; ni++)
	{
		for(nj=0; nj<pTransition->nWidth; nj++)
		{
			*pEle = log(*pEle);
			pEle++;
		}
	}

	pEle = pEmission->pMatElement+pEmission->nWidth;
	for(ni=1; ni<pEmission->nHeight; ni++)
	{
		for(nj=0; nj<pEmission->nWidth; nj++)
		{
			*pEle = log(*pEle);
			pEle++;
		}
	}

	

	/* prepare space */
	pScore = NULL;
	pScore = CreateDoubleMatrix(1, nProbeNum);
	if(pScore == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}

	pPosition = NULL;
	pPosition = CreateDoubleMatrix(1, nProbeNum);
	if(pPosition == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}

	
	/* load data */
	fpIn = NULL;
	fpIn = fopen(strDataPath, "rt");
	if(fpIn == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_Main, cannot open data file!\n");
		exit(EXIT_FAILURE);
	}

	/* fgets(strLine, LINE_LENGTH, fpIn); */

	ni = 0;
	while(fgets(strLine, LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		/* sscanf(strLine, "%s %d %lf", strChr, &nPos, &dScore); */
		sscanf(strLine, "%d %lf", &nPos, &dScore);
		pScore->pMatElement[ni] = dScore;
		pPosition->pMatElement[ni] = nPos;
		ni++;
	}

	fclose(fpIn);

	if(ni != nProbeNum)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_Main, probe number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* transform score */
	pNewScore = NULL;
	pNewScore = Tiling_BindingRegionSelection_ScoreTransform(pScore, dScoreResolution, strTransformType);
	DestroyDoubleMatrix(pScore);

	/* Baum-Welch to estimate parameters */
	pStationaryNew = NULL;
	pStationaryNew = DMCLONE(pStationary);
	if(pStationaryNew == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_Main, cannot create matrix for Baum-Welch!\n");
		exit(EXIT_FAILURE);
	}
	pTransitionNew = NULL;
	pTransitionNew = DMCLONE(pTransition);
	if(pTransitionNew == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_Main, cannot create matrix for Baum-Welch!\n");
		exit(EXIT_FAILURE);
	}
	pEmissionNew = NULL;
	pEmissionNew = DMCLONE(pEmission);
	if(pEmissionNew == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_Main, cannot create matrix for Baum-Welch!\n");
		exit(EXIT_FAILURE);
	}

	Tiling_BindingRegionSelection_HMM_BaumWelch(nProbeNum, nStateNum,
				pStationary, pTransition, pEmission, dGapDist,
				pNewScore, pPosition,
				pStationaryNew, pTransitionNew, pEmissionNew,
				pStationaryTrue, dStopCut, nMaxIter);

	/* save paramters */
	sprintf(strLine, "%s_transitionp.txt", strOutPath);
	DMSAVE(pTransitionNew, strLine);
	sprintf(strLine, "%s_emissionp.txt", strOutPath);
	DMSAVE(pEmissionNew, strLine);
	sprintf(strLine, "%s_stationaryp.txt", strOutPath);
	DMSAVE(pStationaryNew, strLine);
	sprintf(strLine, "%s_stationarytruep.txt", strOutPath);
	DMSAVE(pStationaryTrue, strLine);
	
	/* transform parameters */
	pEle = pTransitionNew->pMatElement;
	for(ni=0; ni<pTransitionNew->nHeight; ni++)
	{
		for(nj=0; nj<pTransitionNew->nWidth; nj++)
		{
			*pEle = log(*pEle);
			pEle++;
		}
	}

	pEle = pEmissionNew->pMatElement+pEmissionNew->nWidth;
	for(ni=1; ni<pEmissionNew->nHeight; ni++)
	{
		for(nj=0; nj<pEmissionNew->nWidth; nj++)
		{
			*pEle = log(*pEle);
			pEle++;
		}
	}

	pEle = pStationaryNew->pMatElement;
	for(ni=0; ni<pStationaryNew->nHeight; ni++)
	{
		for(nj=0; nj<pStationaryNew->nWidth; nj++)
		{
			*pEle = log(*pEle);
			pEle++;
		}
	}


	/* calculate posterior probability */
	vPosterior = NULL;
	vPosterior = (struct DOUBLEMATRIX **)calloc(nStateNum, sizeof(struct DOUBLEMATRIX *));
	if(vPosterior == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nStateNum; ni++)
	{
		vPosterior[ni] = CreateDoubleMatrix(1, nProbeNum);
		if(vPosterior[ni] == NULL)
		{
			printf("Error: Tiling_BindingRegionSelection_HMM_Main, cannot create data space!\n");
			exit(EXIT_FAILURE);
		}
	}

	Tiling_BindingRegionSelection_HMM(nProbeNum, nStateNum,
				pStationaryNew, pTransitionNew, pEmissionNew, dGapDist,
				pNewScore, pPosition, vPosterior);

	/* save */
	sprintf(strLine, "%s_pl.txt", strOutPath);
	fpOut = NULL;
	fpOut = fopen(strLine, "wt");
	if(fpOut == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nProbeNum; ni++)
	{
		for(nj=0; nj<nStateNum; nj++)
		{
			fprintf(fpOut, "%f\t", vPosterior[nj]->pMatElement[ni]);
		}
		fprintf(fpOut, "\n");
	}

	fclose(fpOut);

	/* save to *.bed file */
	sprintf(strLine, "%s_rg.bed", strOutPath);
	Tiling_BindingRegionSelection_OutputToBed(nProbeNum, pPosition, vPosterior, 1, 
		dPosteriorCutoff, dGapDist, strLine);


	/* release memory */
	DestroyDoubleMatrix(pTransition);
	DestroyDoubleMatrix(pEmission);
	DestroyDoubleMatrix(pStationary);
	DestroyDoubleMatrix(pStationaryTrue);
	DestroyDoubleMatrix(pTransitionNew);
	DestroyDoubleMatrix(pEmissionNew);
	DestroyDoubleMatrix(pStationaryNew);
	DestroyDoubleMatrix(pNewScore);
	DestroyDoubleMatrix(pPosition);
	for(ni=0; ni<nStateNum; ni++)
	{
		DestroyDoubleMatrix(vPosterior[ni]);
		vPosterior[ni] = NULL;
	}
	free(vPosterior);

	/* return */
	return PROC_SUCCESS;
}


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
				double dStopCut, int nMaxIter)
{
	/* define */
	struct DOUBLEMATRIX *pStationaryOld;
	struct DOUBLEMATRIX *pTransitionOld;
	struct DOUBLEMATRIX *pEmissionOld;
	int nIter;
	double dPrevLike,dError;

	struct DOUBLEMATRIX **vForwardSum;
	struct DOUBLEMATRIX **vBackwardSum;
	struct DOUBLEMATRIX *pDP;
	struct DOUBLEMATRIX *pTP;
	double dSum,dSumF,dSumB;
	double dTemp,dMax;
	double dDist;
	int ni,nj,nk,nx;
	double dInitMax;

	int nEqualLenInterval;
	double dIntS,dIntE,dIntStep;
	double *pEle,*pEle2;

	/* check */
	if( (nProbeNum <= 0) || (nStateNum <= 0) )
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_BaumWelch, nProbeNum/nStateNum <=0!\n");
		exit(EXIT_FAILURE);
	}
	if((pScore == NULL) || (pStationary == NULL) || (pTransition == NULL) 
		|| (pEmission == NULL) || (pPosition == NULL) || (pTransitionNew == NULL) 
		|| (pEmissionNew == NULL) || (pStationaryNew == NULL))
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_BaumWelch, no input data/parameters!\n");
		exit(EXIT_FAILURE);
	}
	if(pScore->nWidth != nProbeNum)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_BaumWelch, probe number not match!\n");
		exit(EXIT_FAILURE);
	}
	

	/* init */
	dInitMax = -DM_ACCESS_VIOLATION;
	
	/* by default, use intervals of equal length for likelihood calculation */
	nEqualLenInterval = 1;
	if(pEmission->nWidth <= 2)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_BaumWelch, emission probability need to be set for more than two intervals!\n");
		exit(EXIT_FAILURE);
	}

	dIntS = DMGETAT(pEmission, 0, 0);
	dIntE = DMGETAT(pEmission, 0, (pEmission->nWidth-2));
	
	dIntStep = (dIntE-dIntS)/(double)(pEmission->nWidth-2);
	if(dIntStep < 0.0)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_BaumWelch, negative interval step length!\n");
		exit(EXIT_FAILURE);
	}


	/* prepare space */
	vForwardSum = NULL;
	vForwardSum = (struct DOUBLEMATRIX **)calloc(nStateNum, sizeof(struct DOUBLEMATRIX *));
	if(vForwardSum == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_BaumWelch, cannot create space for HMM calculation!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nStateNum; ni++)
	{
		vForwardSum[ni] = CreateDoubleMatrix(1, nProbeNum);
		if(vForwardSum[ni] == NULL)
		{
			printf("Error: Tiling_BindingRegionSelection_HMM_BaumWelch, cannot create space for HMM calculation!\n");
			exit(EXIT_FAILURE);
		}
	}

	vBackwardSum = NULL;
	vBackwardSum = (struct DOUBLEMATRIX **)calloc(nStateNum, sizeof(struct DOUBLEMATRIX *));
	if(vBackwardSum == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_BaumWelch, cannot create space for HMM calculation!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nStateNum; ni++)
	{
		vBackwardSum[ni] = CreateDoubleMatrix(1, nProbeNum);
		if(vBackwardSum[ni] == NULL)
		{
			printf("Error: Tiling_BindingRegionSelection_HMM_BaumWelch, cannot create space for HMM calculation!\n");
			exit(EXIT_FAILURE);
		}
	}

	pDP = NULL;
	pDP = CreateDoubleMatrix(nStateNum, 1);
	if(pDP == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_BaumWelch, cannot create space for HMM calculation!\n");
		exit(EXIT_FAILURE);
	}

	pTP = NULL;
	pTP = CreateDoubleMatrix(nStateNum, nStateNum);
	if(pTP == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_BaumWelch, cannot create space for HMM calculation!\n");
		exit(EXIT_FAILURE);
	}

	pStationaryOld = NULL;
	pStationaryOld = DMCLONE(pStationary);
	if(pStationaryOld == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_BaumWelch, cannot transfer parameter correctly!\n");
		exit(EXIT_FAILURE);
	}
	pEmissionOld = NULL;
	pEmissionOld = DMCLONE(pEmission);
	if(pEmissionOld == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_BaumWelch, cannot transfer parameter correctly!\n");
		exit(EXIT_FAILURE);
	}
	pTransitionOld = NULL;
	pTransitionOld = DMCLONE(pTransition);
	if(pTransitionOld == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_BaumWelch, cannot transfer parameter correctly!\n");
		exit(EXIT_FAILURE);
	}

	/* Baum-Welch */
	dPrevLike = -DM_ACCESS_VIOLATION;
	for(nIter=0; nIter<nMaxIter; nIter++)
	{
		printf("Iter %d...\n", nIter);

		/* forward summation */
		for(ni=0; ni<nStateNum; ni++)
		{
			dTemp = pStationaryOld->pMatElement[ni]+Tiling_BindingRegionSelection_HMM_GetEmission(pEmissionOld, ni, pScore->pMatElement[0], nEqualLenInterval, dIntS, dIntE, dIntStep);
			vForwardSum[ni]->pMatElement[0] = dTemp;
		}

		for(nj=1; nj<nProbeNum; nj++)
		{
			dDist = pPosition->pMatElement[nj] - pPosition->pMatElement[nj-1];
			for(ni=0; ni<nStateNum; ni++)
			{
				dMax = dInitMax;
				for(nk=0; nk<nStateNum; nk++)
				{
					dTemp = vForwardSum[nk]->pMatElement[nj-1] + Tiling_BindingRegionSelection_HMM_GetTransition(pStationaryOld, pTransitionOld, nk, ni, dDist, dGapDist); 
					if(dTemp > dMax)
					{
						dMax = dTemp;
					}
					pDP->pMatElement[nk] = dTemp;
				}

				dSum = 0.0;
				for(nk=0; nk<nStateNum; nk++)
				{
					dTemp = pDP->pMatElement[nk]-dMax;
					dSum += exp(dTemp);
				}
				dSum = log(dSum)+dMax+Tiling_BindingRegionSelection_HMM_GetEmission(pEmissionOld, ni, pScore->pMatElement[nj], nEqualLenInterval, dIntS, dIntE, dIntStep);
				vForwardSum[ni]->pMatElement[nj] = dSum;
			}
		}

		
		dMax = dInitMax;
		for(ni=0; ni<nStateNum; ni++)
		{
			dTemp = vForwardSum[ni]->pMatElement[nj-1]; 
			if(dTemp > dMax)
			{
				dMax = dTemp;
			}
			pDP->pMatElement[ni] = dTemp;
		}

		dSum = 0.0;
		for(ni=0; ni<nStateNum; ni++)
		{
			dTemp = pDP->pMatElement[ni]-dMax;
			dSum += exp(dTemp);
		}
		dSumF = log(dSum)+dMax;


		/* backward summation */
		for(ni=0; ni<nStateNum; ni++)
		{
			vBackwardSum[ni]->pMatElement[nProbeNum-1] = 0.0;
		}
		for(nj=nProbeNum-2; nj>=0; nj--)
		{
			dDist = pPosition->pMatElement[nj+1]-pPosition->pMatElement[nj];
			for(ni=0; ni<nStateNum; ni++)
			{
				dMax = dInitMax;
				for(nk=0; nk<nStateNum; nk++)
				{
					dTemp = vBackwardSum[nk]->pMatElement[nj+1] 
						+ Tiling_BindingRegionSelection_HMM_GetTransition(pStationaryOld, pTransitionOld, ni, nk, dDist, dGapDist) 
						+ Tiling_BindingRegionSelection_HMM_GetEmission(pEmissionOld, nk, pScore->pMatElement[nj+1], nEqualLenInterval, dIntS, dIntE, dIntStep); 
					if(dTemp > dMax)
					{
						dMax = dTemp;
					}
					pDP->pMatElement[nk] = dTemp;
				}

				dSum = 0.0;
				for(nk=0; nk<nStateNum; nk++)
				{
					dTemp = pDP->pMatElement[nk]-dMax;
					dSum += exp(dTemp);
				}
				dSum = log(dSum)+dMax;
				vBackwardSum[ni]->pMatElement[nj] = dSum;
			}
		}

		dMax = dInitMax;
		for(ni=0; ni<nStateNum; ni++)
		{
			dTemp = pStationaryOld->pMatElement[ni] 
				+ Tiling_BindingRegionSelection_HMM_GetEmission(pEmissionOld, ni, pScore->pMatElement[0], nEqualLenInterval, dIntS, dIntE, dIntStep)
				+ vBackwardSum[ni]->pMatElement[0]; 
			if(dTemp > dMax)
			{
				dMax = dTemp;
			}
			pDP->pMatElement[ni] = dTemp;
		}

		dSum = 0.0;
		for(ni=0; ni<nStateNum; ni++)
		{
			dTemp = pDP->pMatElement[ni]-dMax;
			dSum += exp(dTemp);
		}
		dSumB = log(dSum)+dMax;

		if(fabs(dSumF-dSumB) > 1e-3)
		{
			printf("Error: Tiling_BindingRegionSelection_HMM_BaumWelch, forward and backbward summation not match!\n");
			exit(EXIT_FAILURE);
		}
		dSumF = (dSumF+dSumB)/2.0;

		/* posterior calculation */
		pEle = pTransitionNew->pMatElement;
		for(ni=0; ni<pTransitionNew->nHeight; ni++)
		{
			for(nj=0; nj<pTransitionNew->nWidth; nj++)
			{
				*pEle = 0.1;
				pEle++;
			}
		}

		pEle = pEmissionNew->pMatElement+pEmissionNew->nWidth;
		for(ni=1; ni<pEmissionNew->nHeight; ni++)
		{
			for(nj=0; nj<pEmissionNew->nWidth; nj++)
			{
				*pEle = 0.1;
				pEle++;
			}
		}

		pEle = pStationaryNew->pMatElement;
		for(ni=0; ni<pStationaryNew->nHeight; ni++)
		{
			for(nj=0; nj<pStationaryNew->nWidth; nj++)
			{
				*pEle = 0.1;
				pEle++;
			}
		}

		pEle = pStationaryTrue->pMatElement;
		for(ni=0; ni<pStationaryTrue->nHeight; ni++)
		{
			for(nj=0; nj<pStationaryTrue->nWidth; nj++)
			{
				*pEle = 0.1;
				pEle++;
			}
		}


		/* first base */
		dSum = 0.0;
		for(ni=0; ni<nStateNum; ni++)
		{
			dTemp = vForwardSum[ni]->pMatElement[0]+vBackwardSum[ni]->pMatElement[0]-dSumF;
			dTemp = exp(dTemp);
			dSum += dTemp;
			pDP->pMatElement[ni] = dTemp;
		}

		nk = Tiling_BindingRegionSelection_HMM_GetEmissionId(pEmissionOld, pScore->pMatElement[0], nEqualLenInterval, dIntS, dIntE, dIntStep);

		for(ni=0; ni<nStateNum; ni++)
		{
			pDP->pMatElement[ni] /= dSum;
			/* count transition */
			pStationaryNew->pMatElement[ni] += pDP->pMatElement[ni];
			pStationaryTrue->pMatElement[ni] += pDP->pMatElement[ni];
			/* count emission */
			dTemp = DMGETAT(pEmissionNew, ni+1, nk)+pDP->pMatElement[ni];
			DMSETAT(pEmissionNew, ni+1, nk, dTemp);
		}

		/* following bases */
		for(nj=1; nj<nProbeNum; nj++)
		{
			dDist = pPosition->pMatElement[nj]-pPosition->pMatElement[nj-1];
			nk = Tiling_BindingRegionSelection_HMM_GetEmissionId(pEmissionOld, pScore->pMatElement[nj], nEqualLenInterval, dIntS, dIntE, dIntStep);

			dSum = 0.0;
			for(ni=0; ni<nStateNum; ni++)
			{
				dTemp = vForwardSum[ni]->pMatElement[nj]+vBackwardSum[ni]->pMatElement[nj]-dSumF;
				dTemp = exp(dTemp);
				dSum += dTemp;
				pDP->pMatElement[ni] = dTemp;
			}
			if(dDist > dGapDist)
			{
				for(ni=0; ni<nStateNum; ni++)
				{
					pDP->pMatElement[ni] /= dSum;
					/* count transition */
					pStationaryNew->pMatElement[ni] += pDP->pMatElement[ni];
					pStationaryTrue->pMatElement[ni] += pDP->pMatElement[ni];
					/* count emission */
					dTemp = DMGETAT(pEmissionNew, ni+1, nk)+pDP->pMatElement[ni];
					DMSETAT(pEmissionNew, ni+1, nk, dTemp);
				}
			}
			else
			{
				for(ni=0; ni<nStateNum; ni++)
				{
					pDP->pMatElement[ni] /= dSum;
					/* count transition */
					pStationaryTrue->pMatElement[ni] += pDP->pMatElement[ni];
					/* count emission */
					dTemp = DMGETAT(pEmissionNew, ni+1, nk)+pDP->pMatElement[ni];
					DMSETAT(pEmissionNew, ni+1, nk, dTemp);
				}

				/* count transition */
				dSum = 0.0;
				for(ni=0; ni<nStateNum; ni++)
				{
					for(nx=0; nx<nStateNum; nx++)
					{
						dTemp = vForwardSum[ni]->pMatElement[nj-1]+vBackwardSum[nx]->pMatElement[nj]
							+Tiling_BindingRegionSelection_HMM_GetTransition(pStationaryOld, pTransitionOld, ni, nx, dDist, dGapDist) 
							+Tiling_BindingRegionSelection_HMM_GetEmission(pEmissionOld, nx, pScore->pMatElement[nj], nEqualLenInterval, dIntS, dIntE, dIntStep)-dSumF;
						dTemp = exp(dTemp);
						DMSETAT(pTP, ni, nx, dTemp);
						dSum += dTemp;
					}
				}
				for(ni=0; ni<nStateNum; ni++)
				{
					for(nx=0; nx<nStateNum; nx++)
					{
						dTemp = DMGETAT(pTP, ni, nx) / dSum;
						dTemp += DMGETAT(pTransitionNew, ni, nx);
						DMSETAT(pTransitionNew, ni, nx, dTemp);
					}
				}
			}
		}

		/* normalization */
		pEle = pTransitionNew->pMatElement;
		for(ni=0; ni<pTransitionNew->nHeight; ni++)
		{
			dTemp = 0.0;
			pEle2 = pEle;
			for(nj=0; nj<pTransitionNew->nWidth; nj++)
			{
				dTemp += *pEle2;
				pEle2++;
			}
			for(nj=0; nj<pTransitionNew->nWidth; nj++)
			{
				*pEle = *pEle/dTemp;
				pEle++;
			}
		}

		pEle = pEmissionNew->pMatElement+pEmissionNew->nWidth;
		for(ni=1; ni<pEmissionNew->nHeight; ni++)
		{
			dTemp = 0.0;
			pEle2 = pEle;
			for(nj=0; nj<pEmissionNew->nWidth; nj++)
			{
				dTemp += *pEle2;
				pEle2++;
			}
			for(nj=0; nj<pEmissionNew->nWidth; nj++)
			{
				*pEle = *pEle/dTemp;
				pEle++;
			}
		}

		pEle = pStationaryNew->pMatElement;
		for(ni=0; ni<pStationaryNew->nHeight; ni++)
		{
			dTemp = 0.0;
			pEle2 = pEle;
			for(nj=0; nj<pStationaryNew->nWidth; nj++)
			{
				dTemp += *pEle2;
				pEle2++;
			}
			for(nj=0; nj<pStationaryNew->nWidth; nj++)
			{
				*pEle = *pEle/dTemp;
				pEle++;
			}
		}

		pEle = pStationaryTrue->pMatElement;
		for(ni=0; ni<pStationaryTrue->nHeight; ni++)
		{
			dTemp = 0.0;
			pEle2 = pEle;
			for(nj=0; nj<pStationaryTrue->nWidth; nj++)
			{
				dTemp += *pEle2;
				pEle2++;
			}
			for(nj=0; nj<pStationaryTrue->nWidth; nj++)
			{
				*pEle = *pEle/dTemp;
				pEle++;
			}
		}

		/* prepare for the next step */
		dError = fabs(dPrevLike - dSumF);
		if(dError < dStopCut)
			break;

		dPrevLike = dSumF;
		DestroyDoubleMatrix(pStationaryOld);
		DestroyDoubleMatrix(pEmissionOld);
		DestroyDoubleMatrix(pTransitionOld);
		pStationaryOld = NULL;
		pStationaryOld = DMCLONE(pStationaryNew);
		if(pStationaryOld == NULL)
		{
			printf("Error: Tiling_BindingRegionSelection_HMM_BaumWelch, cannot transfer parameter correctly!\n");
			exit(EXIT_FAILURE);
		}
		pEmissionOld = NULL;
		pEmissionOld = DMCLONE(pEmissionNew);
		if(pEmissionOld == NULL)
		{
			printf("Error: Tiling_BindingRegionSelection_HMM_BaumWelch, cannot transfer parameter correctly!\n");
			exit(EXIT_FAILURE);
		}
		pTransitionOld = NULL;
		pTransitionOld = DMCLONE(pTransitionNew);
		if(pTransitionOld == NULL)
		{
			printf("Error: Tiling_BindingRegionSelection_HMM_BaumWelch, cannot transfer parameter correctly!\n");
			exit(EXIT_FAILURE);
		}

		/* transform parameters */
		pEle = pTransitionOld->pMatElement;
		for(ni=0; ni<pTransitionOld->nHeight; ni++)
		{
			for(nj=0; nj<pTransitionOld->nWidth; nj++)
			{
				*pEle = log(*pEle);
				pEle++;
			}
		}

		pEle = pEmissionOld->pMatElement+pEmissionOld->nWidth;
		for(ni=1; ni<pEmissionOld->nHeight; ni++)
		{
			for(nj=0; nj<pEmissionOld->nWidth; nj++)
			{
				*pEle = log(*pEle);
				pEle++;
			}
		}

		pEle = pStationaryOld->pMatElement;
		for(ni=0; ni<pStationaryOld->nHeight; ni++)
		{
			for(nj=0; nj<pStationaryOld->nWidth; nj++)
			{
				*pEle = log(*pEle);
				pEle++;
			}
		}

	}

	/* destroy space */
	for(ni=0; ni<nStateNum; ni++)
	{
		DestroyDoubleMatrix(vForwardSum[ni]);
		vForwardSum[ni] = NULL;
		DestroyDoubleMatrix(vBackwardSum[ni]);
		vBackwardSum[ni] = NULL;
	}
	free(vForwardSum);
	free(vBackwardSum);

	DestroyDoubleMatrix(pDP);
	DestroyDoubleMatrix(pTP);
	DestroyDoubleMatrix(pStationaryOld);
	DestroyDoubleMatrix(pEmissionOld);
	DestroyDoubleMatrix(pTransitionOld);

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  Tiling_BindingRegionSelection_HMM_GetEmission()                        */
/*  Get emission probability for HMM                                       */
/* ----------------------------------------------------------------------- */ 
int Tiling_BindingRegionSelection_HMM_GetEmissionId(struct DOUBLEMATRIX *pEmission, 
					double dScore, int nEqualLenInterval, 
					double dIntS, double dIntE, double dIntStep)
{
	/* define */
	int nIndex;

	/* process */
	if(nEqualLenInterval == 1)
	{
		if( dScore <= dIntS)
		{
			nIndex = 0;
		}
		else
		{
			nIndex = (int)((dScore-dIntS)/dIntStep)+1;
		}
		if(nIndex >= pEmission->nWidth)
			nIndex--;
	}
	else
	{
		for(nIndex=0; nIndex<pEmission->nWidth; nIndex++)
		{
			if(dScore <= DMGETAT(pEmission, 0, nIndex))
				break;
		}
		if(nIndex == pEmission->nWidth)
			nIndex--;
	}

	/* return */
	return nIndex;
}

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
				char strOutPath[])
{
	/* define */
	struct DOUBLEMATRIX *pScore;
	struct DOUBLEMATRIX *pNewScore;
	struct DOUBLEMATRIX *pPosition;
	struct DOUBLEMATRIX *pTransition;
	struct DOUBLEMATRIX *pEmission;
	struct DOUBLEMATRIX *pStationary;
	int nStateNum;
	struct DOUBLEMATRIX **vPosterior;
	double dTemp;
	int ni,nj;
	double *pEle;

	/* for loading data */
	char strLine[LINE_LENGTH];
	char strChr[LINE_LENGTH];
	int nPos;
	double dScore;
	FILE *fpIn;
	FILE *fpOut;

	
	/* check */

	/* init */
	pTransition = NULL;
	pTransition = DMLOAD(strTransitionPath);
	if(pTransition == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_ExpLen_Main, cannot load transition probability!\n");
		exit(EXIT_FAILURE);
	}

	pEmission = NULL;
	pEmission = DMLOAD(strEmissionPath);
	if(pEmission == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_ExpLen_Main, cannot load emission probability!\n");
		exit(EXIT_FAILURE);
	}

	if((pTransition->nHeight+1) != pEmission->nHeight)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_ExpLen_Main, state number not match!\n");
		exit(EXIT_FAILURE);
	}

	nStateNum = pTransition->nHeight;
	if(nStateNum != 2)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_ExpLen_Main, this function only support 2 states now!\n");
		exit(EXIT_FAILURE);
	}

	pStationary = NULL;
	pStationary = CreateDoubleMatrix(1, nStateNum);
	dTemp = DMGETAT(pTransition, 0, 1);
	pStationary->pMatElement[0] = log(1.0-dTemp);
	pStationary->pMatElement[1] = log(dTemp);

	pEle = pTransition->pMatElement;
	for(ni=0; ni<pTransition->nHeight-1; ni++)
	{
		for(nj=0; nj<pTransition->nWidth; nj++)
		{
			*pEle = log(*pEle);
			pEle++;
		}
	}

	pEle = pEmission->pMatElement+pEmission->nWidth;
	for(ni=1; ni<pEmission->nHeight; ni++)
	{
		for(nj=0; nj<pEmission->nWidth; nj++)
		{
			*pEle = log(*pEle);
			pEle++;
		}
	}

	

	/* prepare space */
	pScore = NULL;
	pScore = CreateDoubleMatrix(1, nProbeNum);
	if(pScore == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_ExpLen_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}

	pPosition = NULL;
	pPosition = CreateDoubleMatrix(1, nProbeNum);
	if(pPosition == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_ExpLen_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}

	vPosterior = NULL;
	vPosterior = (struct DOUBLEMATRIX **)calloc(nStateNum, sizeof(struct DOUBLEMATRIX *));
	if(vPosterior == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_ExpLen_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nStateNum; ni++)
	{
		vPosterior[ni] = CreateDoubleMatrix(1, nProbeNum);
		if(vPosterior[ni] == NULL)
		{
			printf("Error: Tiling_BindingRegionSelection_HMM_ExpLen_Main, cannot create data space!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* load data */
	fpIn = NULL;
	fpIn = fopen(strDataPath, "rt");
	if(fpIn == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_ExpLen_Main, cannot open data file!\n");
		exit(EXIT_FAILURE);
	}

	/* fgets(strLine, LINE_LENGTH, fpIn); */

	ni = 0;
	while(fgets(strLine, LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		/* sscanf(strLine, "%s %d %lf", strChr, &nPos, &dScore); */
		sscanf(strLine, "%d %lf", &nPos, &dScore);
		pScore->pMatElement[ni] = dScore;
		pPosition->pMatElement[ni] = nPos;
		ni++;
	}

	fclose(fpIn);

	if(ni != nProbeNum)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_ExpLen_Main, probe number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* transform score */
	pNewScore = NULL;
	pNewScore = Tiling_BindingRegionSelection_ScoreTransform(pScore, dScoreResolution, strTransformType);
	DestroyDoubleMatrix(pScore);

	/* calculate posterior probability */
	Tiling_BindingRegionSelection_HMM_ExpLen(nProbeNum, nStateNum,
				pStationary, pTransition, pEmission, dGapDist,
				pNewScore, pPosition, vPosterior);

	/* save */
	sprintf(strLine, "%s_pl.txt", strOutPath);
	fpOut = NULL;
	fpOut = fopen(strLine, "wt");
	if(fpOut == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_ExpLen_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nProbeNum; ni++)
	{
		for(nj=0; nj<nStateNum; nj++)
		{
			fprintf(fpOut, "%f\t", vPosterior[nj]->pMatElement[ni]);
		}
		fprintf(fpOut, "\n");
	}

	fclose(fpOut);

	/* save to *.bed file */
	sprintf(strLine, "%s_rg.bed", strOutPath);
	Tiling_BindingRegionSelection_OutputToBed(nProbeNum, pPosition, vPosterior, 1, 
		dPosteriorCutoff, dGapDist, strLine);


	/* release memory */
	DestroyDoubleMatrix(pTransition);
	DestroyDoubleMatrix(pEmission);
	DestroyDoubleMatrix(pStationary);
	DestroyDoubleMatrix(pNewScore);
	DestroyDoubleMatrix(pPosition);
	for(ni=0; ni<nStateNum; ni++)
	{
		DestroyDoubleMatrix(vPosterior[ni]);
		vPosterior[ni] = NULL;
	}
	free(vPosterior);

	/* return */
	return PROC_SUCCESS;
}

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
				struct DOUBLEMATRIX **vPosterior)
{
	/* define */
	struct DOUBLEMATRIX **vForwardSum;
	struct DOUBLEMATRIX **vBackwardSum;
	struct DOUBLEMATRIX *pDP;
	double dSum,dSumF,dSumB;
	double dTemp,dMax;
	double dDist;
	int ni,nj,nk;
	double dInitMax;

	int nEqualLenInterval;
	double dIntS,dIntE,dIntStep;

	/* check */
	if( (nProbeNum <= 0) || (nStateNum <= 0) )
	{
		printf("Error: Tiling_BindingRegionSelection_HMM, nProbeNum/nStateNum <=0!\n");
		exit(EXIT_FAILURE);
	}
	if((pScore == NULL) || (pStationary == NULL) || (pTransition == NULL) 
		|| (pEmission == NULL) || (pPosition == NULL) || (vPosterior == NULL))
	{
		printf("Error: Tiling_BindingRegionSelection_HMM, no input data/parameters!\n");
		exit(EXIT_FAILURE);
	}
	if(pScore->nWidth != nProbeNum)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM, probe number not match!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nStateNum; ni++)
	{
		if(vPosterior[ni]->nWidth != nProbeNum)
		{
			printf("Error: Tiling_BindingRegionSelection_HMM, probe number not match!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* init */
	dInitMax = -DM_ACCESS_VIOLATION;
	
	/* by default, use intervals of equal length for likelihood calculation */
	nEqualLenInterval = 1;
	if(pEmission->nWidth <= 2)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM, emission probability need to be set for more than two intervals!\n");
		exit(EXIT_FAILURE);
	}

	dIntS = DMGETAT(pEmission, 0, 0);
	dIntE = DMGETAT(pEmission, 0, (pEmission->nWidth-2));
	
	dIntStep = (dIntE-dIntS)/(double)(pEmission->nWidth-2);
	if(dIntStep < 0.0)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM, negative interval step length!\n");
		exit(EXIT_FAILURE);
	}


	/* prepare space */
	vForwardSum = NULL;
	vForwardSum = (struct DOUBLEMATRIX **)calloc(nStateNum, sizeof(struct DOUBLEMATRIX *));
	if(vForwardSum == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM, cannot create space for HMM calculation!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nStateNum; ni++)
	{
		vForwardSum[ni] = CreateDoubleMatrix(1, nProbeNum);
		if(vForwardSum[ni] == NULL)
		{
			printf("Error: Tiling_BindingRegionSelection_HMM, cannot create space for HMM calculation!\n");
			exit(EXIT_FAILURE);
		}
	}

	vBackwardSum = NULL;
	vBackwardSum = (struct DOUBLEMATRIX **)calloc(nStateNum, sizeof(struct DOUBLEMATRIX *));
	if(vBackwardSum == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM, cannot create space for HMM calculation!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nStateNum; ni++)
	{
		vBackwardSum[ni] = CreateDoubleMatrix(1, nProbeNum);
		if(vBackwardSum[ni] == NULL)
		{
			printf("Error: Tiling_BindingRegionSelection_HMM, cannot create space for HMM calculation!\n");
			exit(EXIT_FAILURE);
		}
	}

	pDP = NULL;
	pDP = CreateDoubleMatrix(nStateNum, 1);
	if(pDP == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM, cannot create space for HMM calculation!\n");
		exit(EXIT_FAILURE);
	}

	/* forward summation */
	for(ni=0; ni<nStateNum; ni++)
	{
		dTemp = pStationary->pMatElement[ni]+Tiling_BindingRegionSelection_HMM_GetEmission(pEmission, ni, pScore->pMatElement[0], nEqualLenInterval, dIntS, dIntE, dIntStep);
		vForwardSum[ni]->pMatElement[0] = dTemp;
	}

	for(nj=1; nj<nProbeNum; nj++)
	{
		dDist = pPosition->pMatElement[nj] - pPosition->pMatElement[nj-1];
		for(ni=0; ni<nStateNum; ni++)
		{
			dMax = dInitMax;
			for(nk=0; nk<nStateNum; nk++)
			{
				dTemp = vForwardSum[nk]->pMatElement[nj-1] + Tiling_BindingRegionSelection_HMM_GetTransition_ExpLen(pStationary, pTransition, nk, ni, dDist, dGapDist); 
				if(dTemp > dMax)
				{
					dMax = dTemp;
				}
				pDP->pMatElement[nk] = dTemp;
			}

			dSum = 0.0;
			for(nk=0; nk<nStateNum; nk++)
			{
				dTemp = pDP->pMatElement[nk]-dMax;
				dSum += exp(dTemp);
			}
			dSum = log(dSum)+dMax+Tiling_BindingRegionSelection_HMM_GetEmission(pEmission, ni, pScore->pMatElement[nj], nEqualLenInterval, dIntS, dIntE, dIntStep);
			vForwardSum[ni]->pMatElement[nj] = dSum;
		}
	}

	
	dMax = dInitMax;
	for(ni=0; ni<nStateNum; ni++)
	{
		dTemp = vForwardSum[ni]->pMatElement[nj-1]; 
		if(dTemp > dMax)
		{
			dMax = dTemp;
		}
		pDP->pMatElement[ni] = dTemp;
	}

	dSum = 0.0;
	for(ni=0; ni<nStateNum; ni++)
	{
		dTemp = pDP->pMatElement[ni]-dMax;
		dSum += exp(dTemp);
	}
	dSumF = log(dSum)+dMax;


	/* backward summation */
	for(ni=0; ni<nStateNum; ni++)
	{
		vBackwardSum[ni]->pMatElement[nProbeNum-1] = 0.0;
	}
	for(nj=nProbeNum-2; nj>=0; nj--)
	{
		dDist = pPosition->pMatElement[nj+1]-pPosition->pMatElement[nj];
		for(ni=0; ni<nStateNum; ni++)
		{
			dMax = dInitMax;
			for(nk=0; nk<nStateNum; nk++)
			{
				dTemp = vBackwardSum[nk]->pMatElement[nj+1] 
					+ Tiling_BindingRegionSelection_HMM_GetTransition_ExpLen(pStationary, pTransition, ni, nk, dDist, dGapDist) 
					+ Tiling_BindingRegionSelection_HMM_GetEmission(pEmission, nk, pScore->pMatElement[nj+1], nEqualLenInterval, dIntS, dIntE, dIntStep); 
				if(dTemp > dMax)
				{
					dMax = dTemp;
				}
				pDP->pMatElement[nk] = dTemp;
			}

			dSum = 0.0;
			for(nk=0; nk<nStateNum; nk++)
			{
				dTemp = pDP->pMatElement[nk]-dMax;
				dSum += exp(dTemp);
			}
			dSum = log(dSum)+dMax;
			vBackwardSum[ni]->pMatElement[nj] = dSum;
		}
	}

	dMax = dInitMax;
	for(ni=0; ni<nStateNum; ni++)
	{
		dTemp = pStationary->pMatElement[ni] 
			+ Tiling_BindingRegionSelection_HMM_GetEmission(pEmission, ni, pScore->pMatElement[0], nEqualLenInterval, dIntS, dIntE, dIntStep)
			+ vBackwardSum[ni]->pMatElement[0]; 
		if(dTemp > dMax)
		{
			dMax = dTemp;
		}
		pDP->pMatElement[ni] = dTemp;
	}

	dSum = 0.0;
	for(ni=0; ni<nStateNum; ni++)
	{
		dTemp = pDP->pMatElement[ni]-dMax;
		dSum += exp(dTemp);
	}
	dSumB = log(dSum)+dMax;

	if(fabs(dSumF-dSumB) > 1e-3)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM, forward and backbward summation not match!\n");
		exit(EXIT_FAILURE);
	}
	dSumF = (dSumF+dSumB)/2.0;

	/* posterior calculation */
	for(nj=0; nj<nProbeNum; nj++)
	{
		dSum = 0.0;
		for(ni=0; ni<nStateNum; ni++)
		{
			dTemp = vForwardSum[ni]->pMatElement[nj]+vBackwardSum[ni]->pMatElement[nj]-dSumF;
			dTemp = exp(dTemp);
			dSum += dTemp;
			vPosterior[ni]->pMatElement[nj] = dTemp;
		}
		for(ni=0; ni<nStateNum; ni++)
		{
			vPosterior[ni]->pMatElement[nj] /= dSum;
		}
	}

	/* destroy space */
	for(ni=0; ni<nStateNum; ni++)
	{
		DestroyDoubleMatrix(vForwardSum[ni]);
		vForwardSum[ni] = NULL;
		DestroyDoubleMatrix(vBackwardSum[ni]);
		vBackwardSum[ni] = NULL;
	}
	free(vForwardSum);
	free(vBackwardSum);

	DestroyDoubleMatrix(pDP);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Tiling_BindingRegionSelection_HMM_GetTransition_ExpLen()               */
/*  Get transition probability for HMM. Transition probability is gap      */
/*  length dependent and exponentially decaying.                           */
/* ----------------------------------------------------------------------- */ 
double Tiling_BindingRegionSelection_HMM_GetTransition_ExpLen(struct DOUBLEMATRIX *pStationary,
					struct DOUBLEMATRIX *pTransition, 
					int nFromS, int nToS, double dDist, double dGapDist)
{
	/* define */
	double dTP;

	/* process */
	if(nFromS == 0)
	{
		dTP = DMGETAT(pTransition, nFromS, nToS);
	}
	else
	{
		if(nToS == 1)
		{
			dTP = dDist*DMGETAT(pTransition, nFromS, nToS);
		}
		else
		{
			dTP = log(1.0-exp(dDist*DMGETAT(pTransition, nFromS, 1)));
		}
	}
	
	
	/* return */
	return dTP;
}


/* ----------------------------------------------------------------------- */ 
/*  Tiling_BindingRegionSelection_HMM_ConstLen_Main()                      */
/*  HMM Forward-Backward algorithm for calculating the posterior prob. of  */
/*  binding for each position. The binding region has a fixed length.      */
/* ----------------------------------------------------------------------- */ 
int Tiling_BindingRegionSelection_HMM_ConstLen_Main(int nProbeNum, char strDataPath[],
				char strTransformType[], double dScoreResolution,
				char strTransitionPath[], char strEmissionPath[], 
				double dGapDist, double dPosteriorCutoff,
				char strOutPath[])
{
	/* define */
	struct DOUBLEMATRIX *pScore;
	struct DOUBLEMATRIX *pNewScore;
	struct DOUBLEMATRIX *pPosition;
	struct DOUBLEMATRIX *pTransition;
	struct DOUBLEMATRIX *pEmission;
	struct DOUBLEMATRIX *pStationary;
	int nStateNum;
	struct DOUBLEMATRIX **vPosterior;
	double dTemp;
	int ni,nj;
	double *pEle;

	/* for loading data */
	char strLine[LINE_LENGTH];
	char strChr[LINE_LENGTH];
	int nPos;
	double dScore;
	FILE *fpIn;
	FILE *fpOut;

	
	/* check */

	/* init */
	pTransition = NULL;
	pTransition = DMLOAD(strTransitionPath);
	if(pTransition == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_ExpLen_Main, cannot load transition probability!\n");
		exit(EXIT_FAILURE);
	}

	pEmission = NULL;
	pEmission = DMLOAD(strEmissionPath);
	if(pEmission == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_ExpLen_Main, cannot load emission probability!\n");
		exit(EXIT_FAILURE);
	}

	if((pTransition->nHeight+1) != pEmission->nHeight)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_ExpLen_Main, state number not match!\n");
		exit(EXIT_FAILURE);
	}

	nStateNum = pTransition->nHeight;
	if(nStateNum != 2)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_ExpLen_Main, this function only support 2 states now!\n");
		exit(EXIT_FAILURE);
	}

	pStationary = NULL;
	pStationary = CreateDoubleMatrix(1, nStateNum);
	dTemp = DMGETAT(pTransition, 0, 1);
	pStationary->pMatElement[0] = log(1.0-dTemp);
	pStationary->pMatElement[1] = log(dTemp);

	pEle = pTransition->pMatElement;
	for(ni=0; ni<pTransition->nHeight-1; ni++)
	{
		for(nj=0; nj<pTransition->nWidth; nj++)
		{
			*pEle = log(*pEle);
			pEle++;
		}
	}

	pEle = pEmission->pMatElement+pEmission->nWidth;
	for(ni=1; ni<pEmission->nHeight; ni++)
	{
		for(nj=0; nj<pEmission->nWidth; nj++)
		{
			*pEle = log(*pEle);
			pEle++;
		}
	}

	

	/* prepare space */
	pScore = NULL;
	pScore = CreateDoubleMatrix(1, nProbeNum);
	if(pScore == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_ExpLen_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}

	pPosition = NULL;
	pPosition = CreateDoubleMatrix(1, nProbeNum);
	if(pPosition == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_ExpLen_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}

	vPosterior = NULL;
	vPosterior = (struct DOUBLEMATRIX **)calloc(nStateNum, sizeof(struct DOUBLEMATRIX *));
	if(vPosterior == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_ExpLen_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nStateNum; ni++)
	{
		vPosterior[ni] = CreateDoubleMatrix(1, nProbeNum);
		if(vPosterior[ni] == NULL)
		{
			printf("Error: Tiling_BindingRegionSelection_HMM_ExpLen_Main, cannot create data space!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* load data */
	fpIn = NULL;
	fpIn = fopen(strDataPath, "rt");
	if(fpIn == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_ExpLen_Main, cannot open data file!\n");
		exit(EXIT_FAILURE);
	}

	/* fgets(strLine, LINE_LENGTH, fpIn); */

	ni = 0;
	while(fgets(strLine, LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		/* sscanf(strLine, "%s %d %lf", strChr, &nPos, &dScore); */
		sscanf(strLine, "%d %lf", &nPos, &dScore);
		pScore->pMatElement[ni] = dScore;
		pPosition->pMatElement[ni] = nPos;
		ni++;
	}

	fclose(fpIn);

	if(ni != nProbeNum)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_ExpLen_Main, probe number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* transform score */
	pNewScore = NULL;
	pNewScore = Tiling_BindingRegionSelection_ScoreTransform(pScore, dScoreResolution, strTransformType);
	DestroyDoubleMatrix(pScore);

	/* calculate posterior probability */
	Tiling_BindingRegionSelection_HMM_ConstLen(nProbeNum, nStateNum,
				pStationary, pTransition, pEmission, dGapDist,
				pNewScore, pPosition, vPosterior);

	/* save */
	sprintf(strLine, "%s_pl.txt", strOutPath);
	fpOut = NULL;
	fpOut = fopen(strLine, "wt");
	if(fpOut == NULL)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM_ExpLen_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nProbeNum; ni++)
	{
		for(nj=0; nj<nStateNum; nj++)
		{
			fprintf(fpOut, "%f\t", vPosterior[nj]->pMatElement[ni]);
		}
		fprintf(fpOut, "\n");
	}

	fclose(fpOut);

	/* save to *.bed file */
	sprintf(strLine, "%s_rg.bed", strOutPath);
	Tiling_BindingRegionSelection_OutputToBed(nProbeNum, pPosition, vPosterior, 1, 
		dPosteriorCutoff, dGapDist, strLine);


	/* release memory */
	DestroyDoubleMatrix(pTransition);
	DestroyDoubleMatrix(pEmission);
	DestroyDoubleMatrix(pStationary);
	DestroyDoubleMatrix(pNewScore);
	DestroyDoubleMatrix(pPosition);
	for(ni=0; ni<nStateNum; ni++)
	{
		DestroyDoubleMatrix(vPosterior[ni]);
		vPosterior[ni] = NULL;
	}
	free(vPosterior);

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  Tiling_BindingRegionSelection_HMM_ConstLen()                           */
/*  HMM Forward-Backward algorithm for calculating the posterior prob. of  */
/*  binding for each position. The length of the binding region is fixed.  */
/* ----------------------------------------------------------------------- */ 
int Tiling_BindingRegionSelection_HMM_ConstLen(int nProbeNum, int nStateNum,
				struct DOUBLEMATRIX *pStationary, struct DOUBLEMATRIX *pTransition, 
				struct DOUBLEMATRIX *pEmission, double dGapDist, 
				struct DOUBLEMATRIX *pScore, struct DOUBLEMATRIX *pPosition, 
				struct DOUBLEMATRIX **vPosterior)
{
	/* define */
	struct DOUBLEMATRIX **vForwardSum;
	struct DOUBLEMATRIX **vBackwardSum;
	struct DOUBLEMATRIX *pDP;
	double dSum,dSumF,dSumB;
	double dTemp,dMax;
	double dDist;
	int ni,nj,nk,nx;
	double dInitMax;
	double dBindingLen;

	int nEqualLenInterval;
	double dIntS,dIntE,dIntStep;
	int nStartId,nEndId,nIgnore;

	/* check */
	if( (nProbeNum <= 0) || (nStateNum <= 0) )
	{
		printf("Error: Tiling_BindingRegionSelection_HMM, nProbeNum/nStateNum <=0!\n");
		exit(EXIT_FAILURE);
	}
	if((pScore == NULL) || (pStationary == NULL) || (pTransition == NULL) 
		|| (pEmission == NULL) || (pPosition == NULL) || (vPosterior == NULL))
	{
		printf("Error: Tiling_BindingRegionSelection_HMM, no input data/parameters!\n");
		exit(EXIT_FAILURE);
	}
	if(pScore->nWidth != nProbeNum)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM, probe number not match!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nStateNum; ni++)
	{
		if(vPosterior[ni]->nWidth != nProbeNum)
		{
			printf("Error: Tiling_BindingRegionSelection_HMM, probe number not match!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* init */
	dBindingLen = DMGETAT(pTransition, 1, 1);
	dInitMax = -DM_ACCESS_VIOLATION;
	
	/* by default, use intervals of equal length for likelihood calculation */
	nEqualLenInterval = 1;
	if(pEmission->nWidth <= 2)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM, emission probability need to be set for more than two intervals!\n");
		exit(EXIT_FAILURE);
	}

	dIntS = DMGETAT(pEmission, 0, 0);
	dIntE = DMGETAT(pEmission, 0, (pEmission->nWidth-2));
	
	dIntStep = (dIntE-dIntS)/(double)(pEmission->nWidth-2);
	if(dIntStep < 0.0)
	{
		printf("Error: Tiling_BindingRegionSelection_HMM, negative interval step length!\n");
		exit(EXIT_FAILURE);
	}


	/* prepare space */
	for(nj=0; nj<nProbeNum; nj++)
	{
		for(nk=nj; nk>=0; nk--)
		{
			if(pPosition->pMatElement[nj]-pPosition->pMatElement[nk] <= dBindingLen)
			{
				vPosterior[0]->pMatElement[nj] += Tiling_BindingRegionSelection_HMM_GetEmission(pEmission, 0, pScore->pMatElement[nk], nEqualLenInterval, dIntS, dIntE, dIntStep);
				vPosterior[1]->pMatElement[nj] += Tiling_BindingRegionSelection_HMM_GetEmission(pEmission, 1, pScore->pMatElement[nk], nEqualLenInterval, dIntS, dIntE, dIntStep);
			}
			else
			{
				break;
			}
		}
		for(nk=nj+1; nk<nProbeNum; nk++)
		{
			if(pPosition->pMatElement[nk]-pPosition->pMatElement[nj] <= dBindingLen)
			{
				vPosterior[0]->pMatElement[nj] += Tiling_BindingRegionSelection_HMM_GetEmission(pEmission, 0, pScore->pMatElement[nk], nEqualLenInterval, dIntS, dIntE, dIntStep);
				vPosterior[1]->pMatElement[nj] += Tiling_BindingRegionSelection_HMM_GetEmission(pEmission, 1, pScore->pMatElement[nk], nEqualLenInterval, dIntS, dIntE, dIntStep);
			}
			else
			{
				break;
			}
		}
		vPosterior[0]->pMatElement[nj] += pTransition->pMatElement[0];
		vPosterior[1]->pMatElement[nj] += pTransition->pMatElement[1];

		dTemp = vPosterior[1]->pMatElement[nj]-vPosterior[0]->pMatElement[nj];
		dTemp = exp(dTemp);
		vPosterior[0]->pMatElement[nj] = 1.0/(1.0+dTemp);
		vPosterior[1]->pMatElement[nj] = 1.0-vPosterior[0]->pMatElement[nj];
	}

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  Tiling_UMS_FDR_Main()                                                  */
/*  UMS for FDR estimation.                                                */
/* ----------------------------------------------------------------------- */ 
int Tiling_UMS_FDR_Main(char strSelectPath[], char strScorePath[], 
						double dPcut, double dQcut,
						int nStepSize, int nIntervalNum,
						char strOutPath[])
{
	/* definition */
	int nProbeNum;
	FILE *fpSelect;
	FILE *fpScore;
	FILE *fpOut;
	struct DOUBLEMATRIX *pPos;
	struct DOUBLEMATRIX *pSelect;
	struct DOUBLEMATRIX *pScore;
	struct DOUBLEMATRIX *pFDR;
	double dPos,dScore;
	char strLine[LINE_LENGTH];
	int ni;

	/* init */
	nProbeNum = 0;
	fpScore = NULL;
	fpScore = fopen(strScorePath, "r");
	if(fpScore == NULL)
	{
		printf("Error: Tiling_UMS_FDR_Main, cannot open selection file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LINE_LENGTH, fpScore) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		ni++;
	}

	fclose(fpScore);
	nProbeNum = ni;

	pPos = NULL;
	pPos = CreateDoubleMatrix(1, nProbeNum);
	if(pPos == NULL)
	{
		printf("Error: Tiling_UMS_FDR_Main, cannot create position matrix!\n");
		exit(EXIT_FAILURE);
	}

	pSelect = NULL;
	pSelect = CreateDoubleMatrix(1, nProbeNum);
	if(pSelect == NULL)
	{
		printf("Error: Tiling_UMS_FDR_Main, cannot create position matrix!\n");
		exit(EXIT_FAILURE);
	}

	pScore = NULL;
	pScore = CreateDoubleMatrix(1, nProbeNum);
	if(pScore == NULL)
	{
		printf("Error: Tiling_UMS_FDR_Main, cannot create position matrix!\n");
		exit(EXIT_FAILURE);
	}

	pFDR = NULL;
	pFDR = CreateDoubleMatrix(1, nProbeNum);
	if(pFDR == NULL)
	{
		printf("Error: Tiling_UMS_FDR_Main, cannot create position matrix!\n");
		exit(EXIT_FAILURE);
	}

	/* load data */
	fpSelect = NULL;
	fpSelect = fopen(strSelectPath, "r");
	if(fpSelect == NULL)
	{
		printf("Error: Tiling_UMS_FDR_Main, cannot open selection file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LINE_LENGTH, fpSelect) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%lf %lf", &dPos, &dScore);
		
		if(ni>=nProbeNum)
		{
			printf("Error: Tiling_UMS_FDR_Main, probe number not match!\n");
			exit(EXIT_FAILURE);
		}
		
		pPos->pMatElement[ni] = dPos;
		pSelect->pMatElement[ni] = dScore;
		ni++;
	}

	fclose(fpSelect);
	
	if(ni != nProbeNum)
	{
		printf("Error: Tiling_UMS_FDR_Main, probe number not match!\n");
		exit(EXIT_FAILURE);
	}

	fpScore = NULL;
	fpScore = fopen(strScorePath, "r");
	if(fpScore == NULL)
	{
		printf("Error: Tiling_UMS_FDR_Main, cannot open selection file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LINE_LENGTH, fpScore) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%lf %lf", &dPos, &dScore);
		
		if(ni>=nProbeNum)
		{
			printf("Error: Tiling_UMS_FDR_Main, probe number not match!\n");
			exit(EXIT_FAILURE);
		}
		
		if((int)dPos != (int)(pPos->pMatElement[ni]))
		{
			printf("Error: Tiling_UMS_FDR_Main, position not match!\n");
			exit(EXIT_FAILURE);
		}
		
		pScore->pMatElement[ni] = dScore;
		ni++;
	}

	fclose(fpScore);
	
	if(ni != nProbeNum)
	{
		printf("Error: Tiling_UMS_FDR_Main, probe number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* UMS */
	Tiling_UMS_FDR(pSelect, dPcut, dQcut, nStepSize, pScore, nIntervalNum, pFDR);

	/* write */
	fpOut = NULL;
	fpOut = fopen(strOutPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: Tiling_UMS_FDR_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "position\tsignal\tums_lfdr\n");
	for(ni=0; ni<nProbeNum; ni++)
	{
		fprintf(fpOut, "%d\t%f\t%f\n", (int)(pPos->pMatElement[ni]),
			pScore->pMatElement[ni], pFDR->pMatElement[ni]);
	}

	fclose(fpOut);


	/* destroy memory */
	DestroyDoubleMatrix(pPos);
	DestroyDoubleMatrix(pSelect);
	DestroyDoubleMatrix(pScore);
	DestroyDoubleMatrix(pFDR);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Tiling_UMS_FDR()                                                       */
/*  UMS for FDR estimation.                                                */
/* ----------------------------------------------------------------------- */ 
int Tiling_UMS_FDR(struct DOUBLEMATRIX *pSelect, double dPcut, double dQcut,
				   int nStepSize, struct DOUBLEMATRIX *pScore, 
				   int nIntervalNum, struct DOUBLEMATRIX *pFDR)
{
	/* define */
	struct DOUBLEMATRIX *pH;
	struct DOUBLEMATRIX *pG0;
	struct DOUBLEMATRIX *pG1;
	struct DOUBLEMATRIX *pF1;
	struct DOUBLEMATRIX *pLikeRatio;
	struct DOUBLEMATRIX *pLfdr;

	int nPrcNum = 100;
	struct DOUBLEMATRIX *pG0Prc;
	struct DOUBLEMATRIX *pG1Prc;
	struct DOUBLEMATRIX *pR;
	double r;
	struct DOUBLEMATRIX *pScorePrctile;
	double dZeroCut = 1e-6;

	struct DOUBLEMATRIX *pSelectSort;
	struct DOUBLEMATRIX *pScoreSort;

	int ni,nj,nk;
	int nSelectWid;
	double dTotNum,dG0Num,dG1Num;
	double dTp,dTq;
	double dScore;
	double dSum;
	double theta;

	/* check parameters */
	if((pSelect == NULL) || (pScore == NULL) || (pFDR == NULL))
	{
		printf("Error: Tiling_UMS_FDR, score/selection/output matrices not specified!\n");
		exit(EXIT_FAILURE);
	}
	if(nStepSize < 0)
	{
		printf("Error: Tiling_UMS_FDR, please set a nonnegative interval number as stepsize of selection criteria!\n");
		exit(EXIT_FAILURE);
	}
	if((pScore->nWidth != pSelect->nWidth) || (pScore->nWidth <= nStepSize))
	{
		printf("Error: Tiling_UMS_FDR, matrix size not match!\n");
		exit(EXIT_FAILURE);
	}
	nSelectWid = pSelect->nWidth - nStepSize;
	if(nIntervalNum <= 0)
	{
		printf("Error: Tiling_UMS_FDR, please set a positive interval number for dividing [0,1]!\n");
		exit(EXIT_FAILURE);
	}

	/* prepare space */
	pH = NULL;
	pH = CreateDoubleMatrix(1, nIntervalNum);
	if(pH == NULL)
	{
		printf("Error: Tiling_UMS_FDR, cannot create memory for UMS!\n");
		exit(EXIT_FAILURE);
	}

	pG0 = NULL;
	pG0 = CreateDoubleMatrix(1, nIntervalNum);
	if(pG0 == NULL)
	{
		printf("Error: Tiling_UMS_FDR, cannot create memory for UMS!\n");
		exit(EXIT_FAILURE);
	}

	pG1 = NULL;
	pG1 = CreateDoubleMatrix(1, nIntervalNum);
	if(pG1 == NULL)
	{
		printf("Error: Tiling_UMS_FDR, cannot create memory for UMS!\n");
		exit(EXIT_FAILURE);
	}

	pF1 = NULL;
	pF1 = CreateDoubleMatrix(1, nIntervalNum);
	if(pF1 == NULL)
	{
		printf("Error: Tiling_UMS_FDR, cannot create memory for UMS!\n");
		exit(EXIT_FAILURE);
	}

	pLikeRatio = NULL;
	pLikeRatio = CreateDoubleMatrix(1, nIntervalNum);
	if(pLikeRatio == NULL)
	{
		printf("Error: Tiling_UMS_FDR, cannot create memory for UMS!\n");
		exit(EXIT_FAILURE);
	}

	pLfdr = NULL;
	pLfdr = CreateDoubleMatrix(1, nIntervalNum);
	if(pLfdr == NULL)
	{
		printf("Error: Tiling_UMS_FDR, cannot create memory for UMS!\n");
		exit(EXIT_FAILURE);
	}


	for(ni=0; ni<nIntervalNum; ni++)
	{
		pH->pMatElement[ni] = 1.0;
		pG0->pMatElement[ni] = 1.0;
		pG1->pMatElement[ni] = 1.0;
	}

	pScorePrctile = NULL;
	pScorePrctile = CreateDoubleMatrix(1, (nPrcNum-1));
	if(pScorePrctile == NULL)
	{
		printf("Error: Tiling_UMS_FDR, cannot create memory for UMS!\n");
		exit(EXIT_FAILURE);
	}

	pG0Prc = NULL;
	pG0Prc = CreateDoubleMatrix(1, (nPrcNum-1));
	if(pG0Prc == NULL)
	{
		printf("Error: Tiling_UMS_FDR, cannot create memory for UMS!\n");
		exit(EXIT_FAILURE);
	}

	pG1Prc = NULL;
	pG1Prc = CreateDoubleMatrix(1, (nPrcNum-1));
	if(pG1Prc == NULL)
	{
		printf("Error: Tiling_UMS_FDR, cannot create memory for UMS!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<(nPrcNum-1); ni++)
	{
		pG0Prc->pMatElement[ni] = 1.0;
		pG1Prc->pMatElement[ni] = 1.0;
	}
	
	/* get percentiles */
	pSelectSort = NULL;
	DMSORTMERGEA_0(pSelect, &pSelectSort, NULL);
	nj = (int)(pSelect->nWidth*dPcut);
	if(nj == pSelect->nWidth)
		nj--;
	dTp = pSelectSort->pMatElement[nj];

	nj = (int)(pSelect->nWidth*dQcut);
	if(nj == pSelect->nWidth)
		nj--;
	dTq = pSelectSort->pMatElement[nj];
	DestroyDoubleMatrix(pSelectSort);

	pScoreSort = NULL;
	DMSORTMERGEA_0(pScore, &pScoreSort, NULL);
	for(ni=1; ni<nPrcNum; ni++)
	{
		nj = (int)((pScore->nWidth)*(double)ni/(double)nPrcNum);
		if(nj == pScore->nWidth)
			nj--;
		pScorePrctile->pMatElement[ni-1] = pScoreSort->pMatElement[nj];
	}

	/* get h(t) */
	for(nj=0; nj<pScore->nWidth; nj++)
	{
		ni = (int)((pScore->pMatElement[nj]) * nIntervalNum);
		if(ni == nIntervalNum)
			ni--;

		pH->pMatElement[ni] = pH->pMatElement[ni]+1.0; 
	}
	dTotNum = pScore->nWidth+nIntervalNum;
	DMPDIVTS(pH, dTotNum);

	/* DMSAVE(pH, "testdistn.txt"); */

	/* get g0(t), g1(t) */
	dG0Num = nIntervalNum;
	dG1Num = nIntervalNum;
	for(nj=0; nj<nSelectWid; nj++)
	{
		/* if g0(t) */
		if(pSelect->pMatElement[nj] > dTp)
		{
			dScore = pScore->pMatElement[nj+nStepSize];
			ni = (int)(dScore * nIntervalNum);
			if(ni == nIntervalNum)
				ni--;

			pG0->pMatElement[ni] = pG0->pMatElement[ni]+1.0;
			dG0Num += 1.0;

			for(nk=(nPrcNum-2); nk>=0; nk--)
			{
				if(dScore >= pScorePrctile->pMatElement[nk])
				{
					pG0Prc->pMatElement[nk] += 1.0; 
				}
			}
		}

		/* if g1(t) */
		if(pSelect->pMatElement[nj] <= dTq)
		{
			dScore = pScore->pMatElement[nj+nStepSize];
			ni = (int)(dScore * nIntervalNum);
			if(ni == nIntervalNum)
				ni--;

			pG1->pMatElement[ni] = pG1->pMatElement[ni]+1.0;
			dG1Num += 1.0;

			for(nk=(nPrcNum-2); nk>=0; nk--)
			{
				if(dScore >= pScorePrctile->pMatElement[nk])
				{
					pG1Prc->pMatElement[nk] += 1.0; 
				}
			}
		}
	}
	DMPDIVTS(pG0, dG0Num);
	DMPDIVTS(pG1, dG1Num);
	DMPDIVTS(pG0Prc, (dG0Num-nIntervalNum+1.0));
	DMPDIVTS(pG1Prc, (dG1Num-nIntervalNum+1.0));
	
	/* DMSAVE(pG0, "testdistng0.txt");
	DMSAVE(pG1, "testdistng1.txt");
	DMSAVE(pG0Prc, "testdistnr0.txt");
	DMSAVE(pG1Prc, "testdistnr1.txt");
	*/

	/* get lim g1(t)/g0(t) */
	for(ni=0; ni<(nPrcNum-1); ni++)
	{
		pG1Prc->pMatElement[ni] /= pG0Prc->pMatElement[ni];
	}
	
	/* DMSAVE(pG1Prc, "testdistnr.txt"); */

	pR = NULL;
	DMSORTMERGEA_0(pG1Prc, &pR, NULL);
	ni = (int)(0.5*(nPrcNum-1));
	r = pR->pMatElement[ni];
	if(r >= 1.0)
	{
		printf("Warning: Tiling_UMS_FDR, r>=1.0, the estimate might be wrong!\n");
		r = 0.999;
	}
	DestroyDoubleMatrix(pR);


	/* estimate f1(t) */
	for(ni=0; ni<nIntervalNum; ni++)
	{
		pF1->pMatElement[ni] = (pG1->pMatElement[ni]-r*pG0->pMatElement[ni])/(1.0-r);
		if(pF1->pMatElement[ni] < dZeroCut)
			pF1->pMatElement[ni] = dZeroCut;
		if(pF1->pMatElement[ni] > 1.0-dZeroCut)
			pF1->pMatElement[ni] = 1.0-dZeroCut;

		pLikeRatio->pMatElement[ni] = pF1->pMatElement[ni]/pG0->pMatElement[ni];
	}
	
	for(ni=1; ni<nIntervalNum; ni++)
	{
		if(pLikeRatio->pMatElement[ni] <= 0.0)
		{
			pLikeRatio->pMatElement[ni] = pLikeRatio->pMatElement[ni-1];
		}
		else if(pLikeRatio->pMatElement[ni] > pLikeRatio->pMatElement[ni-1])
		{
			pLikeRatio->pMatElement[ni] = pLikeRatio->pMatElement[ni-1];
		}
	}

	/* DMSAVE(pLikeRatio, "testlikeratio.txt"); */

	dSum = 0.0;
	for(ni=0; ni<nIntervalNum; ni++)
	{
		pF1->pMatElement[ni] = pG0->pMatElement[ni]*pLikeRatio->pMatElement[ni];
		dSum += pF1->pMatElement[ni];
	}
	DMPDIVTS(pF1, dSum);
	/* DMSAVE(pF1, "testdistf1.txt"); */

	/* estimate theta */
	dScore = 0.0;
	dSum = 0.0;
	for(ni=0; ni<nIntervalNum; ni++)
	{
		dScore += (pH->pMatElement[ni]-pG0->pMatElement[ni])*(pF1->pMatElement[ni]-pG0->pMatElement[ni]);
		dSum += (pF1->pMatElement[ni]-pG0->pMatElement[ni])*(pF1->pMatElement[ni]-pG0->pMatElement[ni]);
	}
	theta = dScore/dSum;
	if(theta > 1.0)
		theta = 1.0;
	if(theta < 0.0)
		theta = 0.0;

	/* estimate lfdr */
	for(ni=0; ni<nIntervalNum; ni++)
	{
		dSum = theta*pF1->pMatElement[ni]+(1.0-theta)*pG0->pMatElement[ni];
		pLfdr->pMatElement[ni] = 1.0-theta*pF1->pMatElement[ni]/dSum;
	}

	/* assign fdr to the original points */
	for(nj=0; nj<pScore->nWidth; nj++)
	{
		ni = (int)(pScore->pMatElement[nj]*nIntervalNum);
		if(ni == nIntervalNum)
			ni--;
		pFDR->pMatElement[nj] = pLfdr->pMatElement[ni];
	}

	/* release memory */
	DestroyDoubleMatrix(pH);
	DestroyDoubleMatrix(pG0);
	DestroyDoubleMatrix(pG1);
	DestroyDoubleMatrix(pF1);
	DestroyDoubleMatrix(pLikeRatio);
	DestroyDoubleMatrix(pLfdr);
	DestroyDoubleMatrix(pG0Prc);
	DestroyDoubleMatrix(pG1Prc);
	DestroyDoubleMatrix(pScorePrctile);
	DestroyDoubleMatrix(pScoreSort);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMap_ImportAffy_Main()                                              */
/*  TileMap loading data from affymetrix's *.CEL files                     */
/* ----------------------------------------------------------------------- */ 
int TileMap_ImportAffy_Main(char strParamPath[])
{
	/* define */
	int nArrayNum = 0;
	int nTotalProbeNum = 0;
	int ni,nlen;
	char strLine[MED_LINE_LENGTH];
	char strWorkPath[MED_LINE_LENGTH];
	char strBpmapPath[MED_LINE_LENGTH];
	struct tagString **vCELPath;
	struct tagString **vAlias;
	int nIntensityType = 0;
	int nLogTransform = 0;
	double dLowerBound = 1.0;
	char strCELName[LINE_LENGTH];
	char strArrayName[LINE_LENGTH];
	char strOutPath[MED_LINE_LENGTH];
	char strCombinedDataPath[MED_LINE_LENGTH];
	char strPosPath[MED_LINE_LENGTH];

	FILE *fpIn;
	char *chSep;
	int nError = 0;

	/* load array list */
	fpIn = NULL;
	fpIn = fopen(strParamPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileMap_ImportAffy_Main, cannot open the parameter file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strstr(strLine, "[Array number]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			nArrayNum = atoi(chSep);
			if(nArrayNum <= 0)
			{
				printf("Error: TileMap_ImportAffy_Main, no arrays available!\n");
				nError = 1;
				break;
			}

			vCELPath = NULL;
			vCELPath = (struct tagString **)calloc(nArrayNum, sizeof(struct tagString *));
			if(vCELPath == NULL)
			{
				printf("Error: TileMap_ImportAffy_Main, cannot allocate memory for tracking *.CEL files!\n");
				exit(EXIT_FAILURE);
			}

			vAlias = NULL;
			vAlias = (struct tagString **)calloc(nArrayNum, sizeof(struct tagString *));
			if(vAlias == NULL)
			{
				printf("Error: TileMap_ImportAffy_Main, cannot allocate memory for tracking *.CEL files!\n");
				exit(EXIT_FAILURE);
			}
		}
		else if(strstr(strLine, "[BPMAP file]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strBpmapPath, chSep);
		}
		else if(strstr(strLine, "[Output file]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strCombinedDataPath, chSep);
		}
		else if(strstr(strLine, "[Working directory]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strWorkPath, chSep);
			
			if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
			{
				nlen = strlen(strWorkPath);
				if(strWorkPath[nlen-1] != '\\')
				{
					strWorkPath[nlen] = '\\';
					strWorkPath[nlen+1] = '\0';
				}
			}
			else
			{
				nlen = strlen(strWorkPath);
				if(strWorkPath[nlen-1] != '/')
				{
					strWorkPath[nlen] = '/';
					strWorkPath[nlen+1] = '\0';
				}
			}
		}

		else if(strstr(strLine, "[How to compute intensity]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			nIntensityType = atoi(chSep);
		}
		else if(strstr(strLine, "[Truncate lower bound]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			dLowerBound = atof(chSep);
		}
		else if(strstr(strLine, "[Take log2 transformation]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			nLogTransform = atoi(chSep);
		}
		else if(strstr(strLine, "[Arrays]") == strLine)
		{
			ni = 0;
			while(ni < nArrayNum)
			{
				fgets(strLine, MED_LINE_LENGTH, fpIn);
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] == '\0')
					continue;

				sscanf(strLine, "%s %s", strCELName, strArrayName);
				StringAddTail(vCELPath+ni, strCELName);
				StringAddTail(vAlias+ni, strArrayName);
				ni++;
			}
		}
		else
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	fclose(fpIn);

	if(nError > 0)
	{
		printf("Error: TileMap_ImportAffy_Main, cannot load the parameter file correctly!\n");
		exit(EXIT_FAILURE);
	}

	/* prepare local repeat filters */
	printf("Prepare Local Repeat Mask...\n");
	sprintf(strOutPath, "%s%s.refmask", strWorkPath, strBpmapPath);
	sprintf(strLine, "%s%s", strWorkPath, strBpmapPath);
	sprintf(strPosPath, "%s%s.pbpos", strWorkPath, strBpmapPath);
	Affy_LoadBPMAP_TileMap(strLine, strPosPath, strOutPath, &nTotalProbeNum);
	printf("\narray number = %d\n", nArrayNum);
	printf("probe number = %d\n", nTotalProbeNum);

	/* load *.cel file and sort the probes according to *.bpmap file */
	for(ni=0; ni<nArrayNum; ni++)
	{
		printf("Loading %s...\n", vCELPath[ni]->m_pString);
		TileMap_ImportCELIntensity(strWorkPath, strBpmapPath, 
			vCELPath[ni]->m_pString, vAlias[ni]->m_pString,
			nIntensityType, dLowerBound, nLogTransform);
	}

	/* combine separate intensity files into one file */
	printf("Merging files...\n");
	TileMap_CombineIntensity(nArrayNum, nTotalProbeNum, strWorkPath, strBpmapPath, vAlias, strCombinedDataPath);

	/* remove files */
	if(strcmp(OS_SYSTEM, "WINDOWS")==0)
	{
		for(ni=0; ni<nArrayNum; ni++)
		{
			sprintf(strLine, "del %s%s.intensity", strWorkPath, vAlias[ni]->m_pString);
			system(strLine);
		}
		sprintf(strLine, "del %s%s.pbpos", strWorkPath, strBpmapPath);
		system(strLine);
	}
	else
	{
		for(ni=0; ni<nArrayNum; ni++)
		{
			sprintf(strLine, "rm %s%s.intensity", strWorkPath, vAlias[ni]->m_pString);
			system(strLine);
		}
		sprintf(strLine, "rm %s%s.pbpos", strWorkPath, strBpmapPath);
		system(strLine);
	}

	/* destroy */
	for(ni=0; ni<nArrayNum; ni++)
	{
		DeleteString(vCELPath[ni]);
		DeleteString(vAlias[ni]);
	}
	free(vCELPath);
	free(vAlias);

	printf("Done!\n");

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMap_ImportCELIntensity()                                           */
/*  TileMap loading data from a single affymetrix's *.CEL file             */
/* ----------------------------------------------------------------------- */ 
int TileMap_ImportCELIntensity(char strWorkPath[], char strBpmapPath[],  
			char strCelFile[], char strAlias[],
			int nIntensityType, double dLowerBound, int nLogTransform)
{
	/* define */
	FILE *fpCel;
	FILE *fpBpmap;
	FILE *fpOut;
	struct DOUBLEMATRIX *pMean;
	struct DOUBLEMATRIX *pSD;
	
	char strLine[MED_LINE_LENGTH];
	double dVersion;
	int nCols, nRows;
	int nTotalX,nTotalY;
	int nOffsetX,nOffsetY;
	int nULx,nURx,nLRx,nLLx,nULy,nURy,nLRy,nLLy;
	int nInvertX, nInvertY;
	int swapXY;
	char strDatHeader[MED_LINE_LENGTH];
	char strAlgorithm[MED_LINE_LENGTH];
	char strAlgorithmParameters[MED_LINE_LENGTH];
	int nNumberCells;

	char *chSep;

	int ni,nj,nX,nY,nP;
	double dM,dS;

	/* variables for loading *.bpmap */
	char strFileType[9];
	float fVersion;
	unsigned long nSeqNum;

	/* seq info */
	struct tagString **vSeqName;
	struct INTMATRIX *vProbePairNum;
	unsigned long nSeqNameLen;
	unsigned long nProbePairNum;

	/* map info */
	struct INTMATRIX *vSeqID;
	unsigned long nSeqID;
	int nProbeNum;
	struct tagAffyBpMapUnit *pNewUnit;
	double dPM,dMM;

	/* init */
	nTotalX = 0;
	nTotalY = 0;

	/* load *.CEL */
	sprintf(strLine, "%s%s", strWorkPath, strCelFile);
	fpCel = NULL;
	fpCel = fopen(strLine, "r");
	if(fpCel == NULL)
	{
		printf("Error: TileMap_ImportCELIntensity, cannot open *.CEL file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, MED_LINE_LENGTH, fpCel)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strstr(strLine, "[CEL]") == strLine)
		{
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			chSep = strchr(strLine, '=');
			chSep++;
			dVersion = atof(chSep);
		}
		else if(strstr(strLine, "[HEADER]") == strLine)
		{
			/* cols */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "Cols") != strLine)
			{
				printf("Error: TileMap_ImportCELIntensity, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			nCols = atoi(chSep);

			/* rows */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "Rows") != strLine)
			{
				printf("Error: TileMap_ImportCELIntensity, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			nRows = atoi(chSep);

			/* total x */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "TotalX") != strLine)
			{
				printf("Error: TileMap_ImportCELIntensity, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			nTotalX = atoi(chSep);

			/* total y */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "TotalY") != strLine)
			{
				printf("Error: TileMap_ImportCELIntensity, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			nTotalY = atoi(chSep);

			/* offset x */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "OffsetX") != strLine)
			{
				printf("Error: TileMap_ImportCELIntensity, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			nOffsetX = atoi(chSep);

			/* offset y */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "OffsetY") != strLine)
			{
				printf("Error: TileMap_ImportCELIntensity, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			nOffsetY = atoi(chSep);

			/* UL */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "GridCornerUL") != strLine)
			{
				printf("Error: TileMap_ImportCELIntensity, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			sscanf(chSep, "%d %d", &nULx, &nULy);

			/* UR */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "GridCornerUR") != strLine)
			{
				printf("Error: TileMap_ImportCELIntensity, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			sscanf(chSep, "%d %d", &nURx, &nURy);

			/* LR */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "GridCornerLR") != strLine)
			{
				printf("Error: TileMap_ImportCELIntensity, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			sscanf(chSep, "%d %d", &nLRx, &nLRy);

			/* LL */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "GridCornerLL") != strLine)
			{
				printf("Error: TileMap_ImportCELIntensity, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			sscanf(chSep, "%d %d", &nLLx, &nLLy);

			/* invertx */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "Axis") != strLine)
			{
				printf("Error: TileMap_ImportCELIntensity, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			nInvertX = atoi(chSep);

			/* inverty */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "Axis") != strLine)
			{
				printf("Error: TileMap_ImportCELIntensity, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			nInvertY = atoi(chSep);

			/* swapxy */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "swapXY") != strLine)
			{
				printf("Error: TileMap_ImportCELIntensity, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			swapXY = atoi(chSep);

			/* DatHeader */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "DatHeader") != strLine)
			{
				printf("Error: TileMap_ImportCELIntensity, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			strcpy(strDatHeader, chSep);

			/* Algorithm */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "Algorithm") != strLine)
			{
				printf("Error: TileMap_ImportCELIntensity, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			strcpy(strAlgorithm, chSep);

			/* Algorithm Parameters */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "AlgorithmParameters") != strLine)
			{
				printf("Error: TileMap_ImportCELIntensity, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			strcpy(strAlgorithmParameters, chSep);
		}

		else if(strstr(strLine, "[INTENSITY]") == strLine)
		{
			/* NumberCells */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "NumberCells") != strLine)
			{
				printf("Error: TileMap_ImportCELIntensity, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			nNumberCells = atoi(chSep);

			/* cell header */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "CellHeader") != strLine)
			{
				printf("Error: TileMap_ImportCELIntensity, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}

			/* allocate memory */
			if((int)(nTotalX*nTotalY) != nNumberCells)
			{
				printf("Error: TileMap_ImportCELIntensity, loading *.CEL file incorrectly, cell number not match!\n");
				exit(EXIT_FAILURE);
			}

			pMean = NULL;
			pMean = CreateDoubleMatrix(nTotalX, nTotalY);
			if(pMean == NULL)
			{
				printf("Error: TileMap_ImportCELIntensity, cannot allocate enough memory for loading *.CEL file!\n");
				exit(EXIT_FAILURE);
			}

			pSD = NULL;
			pSD = CreateDoubleMatrix(nTotalX, nTotalY);
			if(pSD == NULL)
			{
				printf("Error: TileMap_ImportCELIntensity, cannot allocate enough memory for loading *.CEL file!\n");
				exit(EXIT_FAILURE);
			}

			/* load intensity */
			for(ni=0; ni<nNumberCells; ni++)
			{
				fgets(strLine, MED_LINE_LENGTH, fpCel);
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				sscanf(strLine, "%d %d %lf %lf %d", &nX, &nY, &dM, &dS, &nP);
				DMSETAT(pMean, nX, nY, dM);
				DMSETAT(pSD, nX, nY, dS);
			}
		}

		else if(strstr(strLine, "[MASKS]") == strLine)
		{
			/* NumberCells */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "NumberCells") != strLine)
			{
				printf("Error: TileMap_ImportCELIntensity, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			nNumberCells = atoi(chSep);

			/* cell header */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "CellHeader") != strLine)
			{
				printf("Error: TileMap_ImportCELIntensity, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}

			/* load cells */
			for(ni=0; ni<nNumberCells; ni++)
			{
				fgets(strLine, MED_LINE_LENGTH, fpCel);
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
			}
		}

		else if(strstr(strLine, "[OUTLIERS]") == strLine)
		{
			/* NumberCells */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "NumberCells") != strLine)
			{
				printf("Error: TileMap_ImportCELIntensity, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			nNumberCells = atoi(chSep);

			/* cell header */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "CellHeader") != strLine)
			{
				printf("Error: TileMap_ImportCELIntensity, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}

			/* load cells */
			for(ni=0; ni<nNumberCells; ni++)
			{
				fgets(strLine, MED_LINE_LENGTH, fpCel);
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				sscanf(strLine, "%d %d", &nX, &nY);
			}
		}

		else if(strstr(strLine, "[MODIFIED]") == strLine)
		{
			/* NumberCells */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "NumberCells") != strLine)
			{
				printf("Error: TileMap_ImportCELIntensity, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			nNumberCells = atoi(chSep);

			/* cell header */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "CellHeader") != strLine)
			{
				printf("Error: TileMap_ImportCELIntensity, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}

			/* load cells */
			for(ni=0; ni<nNumberCells; ni++)
			{
				fgets(strLine, MED_LINE_LENGTH, fpCel);
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
			}
		}

		else
		{
		}
	}

	fclose(fpCel);

	/* load .bpmap file */
	sprintf(strLine, "%s%s", strWorkPath, strBpmapPath);
	fpBpmap = NULL;
	fpBpmap = fopen(strLine, "rb");
	if(fpBpmap == NULL)
	{
		printf("Error: TileMap_ImportCELIntensity, cannot open bpmap file!\n");
		exit(EXIT_FAILURE);
	}

	sprintf(strLine, "%s%s.intensity", strWorkPath, strAlias);
	fpOut = NULL;
	fpOut = fopen(strLine, "w");
	if(fpOut == NULL)
	{
		printf("Error: TileMap_ImportCELIntensity, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	
	/* load head */
	fread( strFileType, sizeof(char), 8, fpBpmap );
	strFileType[8] = '\0';
	AFFYBAR_READ_FLOAT(fpBpmap, &fVersion);
	AFFYBAR_READ_ULONG(fpBpmap, &nSeqNum);

	/* load sequence names */
	vSeqName = NULL;
	vSeqName = (struct tagString **)calloc(nSeqNum, sizeof(struct tagString *));
	if(vSeqName == NULL)
	{
		printf("Error: TileMap_ImportCELIntensity, cannot load sequence name!\n");
		exit(EXIT_FAILURE);
	}
	vProbePairNum = NULL;
	vProbePairNum = CreateIntMatrix(nSeqNum,1);
	if(vProbePairNum == NULL)
	{
		printf("Error: TileMap_ImportCELIntensity, cannot load probe pair number!\n");
		exit(EXIT_FAILURE);
	}
	
	for(ni=0; ni<(int)nSeqNum; ni++)
	{
		AFFYBAR_READ_ULONG(fpBpmap, &nSeqNameLen);
		vSeqName[ni] = CreateString(nSeqNameLen);
		fread( vSeqName[ni]->m_pString, sizeof(char), nSeqNameLen, fpBpmap );
		vSeqName[ni]->m_pString[nSeqNameLen] = '\0';
		AFFYBAR_READ_ULONG(fpBpmap, &nProbePairNum);
		vProbePairNum->pMatElement[ni] = (int)nProbePairNum;
	}

	/* load map */
	vSeqID = NULL;
	vSeqID = CreateIntMatrix(nSeqNum,1);
	if(vSeqID == NULL)
	{
		printf("Error: Affy_LoadBPMAP, cannot load sequence id!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<(int)nSeqNum; ni++)
	{
		/* load seq id */
		AFFYBAR_READ_ULONG(fpBpmap, &nSeqID);
		nProbeNum = vProbePairNum->pMatElement[ni];
			
		/* load probes */
		for(nj=0; nj<nProbeNum; nj++)
		{
			/* load new unit */
			pNewUnit = NULL;
			pNewUnit = AffyBpMapUnitCreate();
			AffyBpMapUnitLoad(pNewUnit, fpBpmap);

			dPM = DMGETAT(pMean, pNewUnit->nPMX, pNewUnit->nPMY);
			dMM = DMGETAT(pMean, pNewUnit->nMMX, pNewUnit->nMMY);

			/* if PM-MM */
			if(nIntensityType == 1)
			{
				dPM -= dMM;
			}
			
			if(dPM < dLowerBound)
			{
				dPM = dLowerBound;
			}
			if(nLogTransform == 1)
			{
				dPM = log(dPM)/log(2.0);
			}

			fprintf(fpOut, "%f\n", dPM);

			AffyBpMapUnitDestroy(pNewUnit);			
		}
	}

	/* load tail if any */
	while(feof(fpBpmap) == 0)
	{
		fread(strFileType, sizeof(char), 1, fpBpmap);
		/* printf("%c", strFileType[0]); */
	}

	/* clear memeory */
	for(ni=0; ni<(int)nSeqNum; ni++)
	{
		DeleteString(vSeqName[ni]);
		vSeqName[ni] = NULL;
	}
	free(vSeqName);
	DestroyIntMatrix(vProbePairNum);
	DestroyIntMatrix(vSeqID);

	/* close file */
	fclose(fpBpmap);
	fclose(fpOut);


	/* destroy */
	DestroyDoubleMatrix(pMean);
	DestroyDoubleMatrix(pSD);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMap_CombineIntensity()                                             */
/*  TileMap combine intensities into one file.                             */
/* ----------------------------------------------------------------------- */ 
int TileMap_CombineIntensity(int nArrayNum, int nTotalProbeNum, char strWorkPath[], 
							 char strBpmapPath[], struct tagString **vAlias, 
							 char strOutFile[])
{
	/* define */
	FILE **vfpIn;
	FILE *fpPos;
	FILE *fpOut;
	int ni,nj;
	char strLine[MED_LINE_LENGTH];

	/* open files */
	vfpIn = NULL;
	vfpIn = (FILE **)calloc(nArrayNum, sizeof(FILE*));
	if(vfpIn == NULL)
	{
		printf("Error: TileMap_CombineIntensity, cannot create file pointers to separate intensity files!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nArrayNum; ni++)
	{
		sprintf(strLine, "%s%s.intensity", strWorkPath, vAlias[ni]->m_pString);
		vfpIn[ni] = fopen(strLine, "r");
		if(vfpIn[ni] == NULL)
		{
			printf("Error: TileMap_CombineIntensity, cannot open intensity files!\n");
			exit(EXIT_FAILURE);
		}
	}

	sprintf(strLine, "%s%s.pbpos", strWorkPath, strBpmapPath);
	fpPos = NULL;
	fpPos = fopen(strLine, "r");
	if(fpPos == NULL)
	{
		printf("Error: TileMap_CombineIntensity, cannot open the position file!\n");
		exit(EXIT_FAILURE);
	}

	sprintf(strLine, "%s%s", strWorkPath, strOutFile);
	fpOut = NULL;
	fpOut = fopen(strLine, "w");
	if(fpOut == NULL)
	{
		printf("Error: TileMap_CombineIntensity, cannot open the output file!\n");
		exit(EXIT_FAILURE);
	}

	/* write */
	fprintf(fpOut, "chromosome\tposition");
	for(ni=0; ni<nArrayNum; ni++)
	{
		fprintf(fpOut, "\t%s", vAlias[ni]->m_pString);
	}
	fprintf(fpOut, "\n");

	for(nj=0; nj<nTotalProbeNum; nj++)
	{
		fgets(strLine, MED_LINE_LENGTH, fpPos);
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		fprintf(fpOut, "%s", strLine);

		for(ni=0; ni<nArrayNum; ni++)
		{
			fgets(strLine, MED_LINE_LENGTH, vfpIn[ni]);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			fprintf(fpOut, "\t%s", strLine);
		}
		fprintf(fpOut, "\n");
	}

	/* close files */
	for(ni=0; ni<nArrayNum; ni++)
	{
		fclose(vfpIn[ni]);
	}
	free(vfpIn);

	fclose(fpPos);
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  TileMap_ImportAffy_Normalization_Main()                                */
/*  TileMap loading data from affymetrix's *.CEL files, do normalizations, */
/*  and compute intensities.                                               */
/* ----------------------------------------------------------------------- */ 
int TileMap_ImportAffy_Normalization_Main(char strParamPath[])
{
	/* define */

	/* working path */
	char strWorkPath[MED_LINE_LENGTH];
	char strBpmapPath[MED_LINE_LENGTH];
	char strCombinedDataPath[MED_LINE_LENGTH];
	
	/* arrays */
	int nArrayNum = 0;
	struct tagString **vCELPath;
	struct tagString **vAlias;
	int nTotalProbeNum = 0;
	int nRealProbeNum = 0;
	struct DOUBLEMATRIX **vArray;
	struct DOUBLEMATRIX **vSD;
	int nTotalX, nTotalY, nArrayX, nArrayY;
	
	/* normalization */
	int nIncludeNormalization = 1;
	double dNormLowerBound = 0.0;
	int nNormLogTransform = 0;
	double dLog2;
	
	/* intensity computation */
	int nIntensityType = 0;
	double dIntLowerBound = 1.0;
	int nIntLogTransform = 1;

	/* others */
	int ni,nj,nlen;
	char strLine[MED_LINE_LENGTH];
	char strCELName[LINE_LENGTH];
	char strArrayName[LINE_LENGTH];
	char strOutPath[MED_LINE_LENGTH];
	char strMaskPath[MED_LINE_LENGTH];

	FILE *fpIn;
	char *chSep;
	int nError = 0;

	/* init */
	nTotalX = 0;
	nTotalY = 0;
	nArrayX = 0;
	nArrayY = 0;

	/* load parameters */
	fpIn = NULL;
	fpIn = fopen(strParamPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileMap_ImportAffy_Normalization_Main, cannot open the parameter file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strstr(strLine, "[Working directory]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strWorkPath, chSep);
			
			if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
			{
				nlen = strlen(strWorkPath);
				if(strWorkPath[nlen-1] != '\\')
				{
					strWorkPath[nlen] = '\\';
					strWorkPath[nlen+1] = '\0';
				}
			}
			else
			{
				nlen = strlen(strWorkPath);
				if(strWorkPath[nlen-1] != '/')
				{
					strWorkPath[nlen] = '/';
					strWorkPath[nlen+1] = '\0';
				}
			}
		}

		else if(strstr(strLine, "[BPMAP file]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strBpmapPath, chSep);
		}
		
		else if(strstr(strLine, "[Export file]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strCombinedDataPath, chSep);
		}

		else if(strstr(strLine, "[Array number]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			nArrayNum = atoi(chSep);
			if(nArrayNum <= 0)
			{
				printf("Error: TileMap_ImportAffy_Main, no arrays available!\n");
				nError = 1;
				break;
			}

			vCELPath = NULL;
			vCELPath = (struct tagString **)calloc(nArrayNum, sizeof(struct tagString *));
			if(vCELPath == NULL)
			{
				printf("Error: TileMap_ImportAffy_Main, cannot allocate memory for tracking *.CEL files!\n");
				exit(EXIT_FAILURE);
			}

			vAlias = NULL;
			vAlias = (struct tagString **)calloc(nArrayNum, sizeof(struct tagString *));
			if(vAlias == NULL)
			{
				printf("Error: TileMap_ImportAffy_Main, cannot allocate memory for tracking *.CEL files!\n");
				exit(EXIT_FAILURE);
			}
		}
		
		else if(strstr(strLine, "[Arrays]") == strLine)
		{
			ni = 0;
			while(ni < nArrayNum)
			{
				fgets(strLine, MED_LINE_LENGTH, fpIn);
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] == '\0')
					continue;

				sscanf(strLine, "%s %s", strCELName, strArrayName);
				StringAddTail(vCELPath+ni, strCELName);
				StringAddTail(vAlias+ni, strArrayName);
				ni++;
			}
		}

		else if(strstr(strLine, "[Apply normalization before computing intensity]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			nIncludeNormalization = atoi(chSep);
		}

		else if(strstr(strLine, "[Truncation lower bound before normalization]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			dNormLowerBound = atof(chSep);
		}
		
		else if(strstr(strLine, "[Take log2 transformation before normalization]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			nNormLogTransform = atoi(chSep);
		}
	
		else if(strstr(strLine, "[How to compute intensity]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			nIntensityType = atoi(chSep);
		}

		else if(strstr(strLine, "[Truncation lower bound after intensity computation]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			dIntLowerBound = atof(chSep);
		}
		
		else if(strstr(strLine, "[Take log2 transformation after intensity computation]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			nIntLogTransform = atoi(chSep);
		}
		
		else if(strLine[0] == '#')
		{
		}

		else
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	fclose(fpIn);

	if(nError > 0)
	{
		printf("Error: TileMap_ImportAffy_Normalization_Main, cannot load the parameter file correctly!\n");
		exit(EXIT_FAILURE);
	}

	/* load *.cel file and sort the probes according to *.bpmap file */
	vArray = NULL;
	vArray = (struct DOUBLEMATRIX **)calloc(nArrayNum, sizeof(struct DOUBLEMATRIX *));
	if(vArray == NULL)
	{
		printf("Error: TileMap_ImportAffy_Normalization_Main, allocate memory for loading arrays!\n");
		exit(EXIT_FAILURE);
	}

	vSD = NULL;
	vSD = (struct DOUBLEMATRIX **)calloc(nArrayNum, sizeof(struct DOUBLEMATRIX *));
	if(vSD == NULL)
	{
		printf("Error: TileMap_ImportAffy_Normalization_Main, allocate memory for loading arrays!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nArrayNum; ni++)
	{
		printf("Loading %s...\n", vCELPath[ni]->m_pString);
		TileMap_LoadCEL(strWorkPath, vCELPath[ni]->m_pString, &nArrayX, &nArrayY, (vArray+ni), (vSD+ni));
		DestroyDoubleMatrix(vSD[ni]);

		if(ni==0)
		{
			nTotalX = nArrayX;
			nTotalY = nArrayY;
			nTotalProbeNum = vArray[ni]->nWidth;
		}
		else
		{
			if( (nTotalX != nArrayX) || (nTotalY != nArrayY) || (nTotalProbeNum != vArray[ni]->nWidth) )
			{
				printf("Error: TileMap_ImportAffy_Normalization_Main, array size not match!\n");
				exit(EXIT_FAILURE);
			}
		}
	}

	/* normalization */
	dLog2 = log(2.0);
	if(nIncludeNormalization == 1)
	{
		printf("Normalization...\n");
		
		/* preprocessing */
		if(nNormLogTransform == 1)
		{
			for(ni=0; ni<nArrayNum; ni++)
			{
				for(nj=0; nj<nTotalProbeNum; nj++)
				{
					/* truncate */
					if(vArray[ni]->pMatElement[nj] < dNormLowerBound)
						vArray[ni]->pMatElement[nj] = dNormLowerBound;

					/* log transformation */
					vArray[ni]->pMatElement[nj] = log(vArray[ni]->pMatElement[nj])/dLog2;
				}
			}
		}
		else
		{
			for(ni=0; ni<nArrayNum; ni++)
			{
				for(nj=0; nj<nTotalProbeNum; nj++)
				{
					/* truncate */
					if(vArray[ni]->pMatElement[nj] < dNormLowerBound)
						vArray[ni]->pMatElement[nj] = dNormLowerBound;
				}
			}
		}
		
		/* normalization */
		Expression_Normalization_Quantile(nArrayNum, nTotalProbeNum, vArray);
	}

	/* export intensities and prepare local repeat filters */
	printf("Export Intensities and Local Repeat Mask...\n");
	sprintf(strLine, "%s%s", strWorkPath, strBpmapPath);
	sprintf(strMaskPath, "%s%s.refmask", strWorkPath, strBpmapPath);
	sprintf(strOutPath, "%s%s", strWorkPath, strCombinedDataPath);
	TileMap_ExportAffyIntensity(nArrayNum, nTotalProbeNum, &nRealProbeNum,
		vArray, vAlias, nTotalX, nTotalY, 
		strLine, strMaskPath, strOutPath, nIntensityType, 
		dIntLowerBound, nIntLogTransform);
	printf("\nArray Number = %d\n", nArrayNum);
	printf("Probe number = %d\n", nRealProbeNum);
	
	/* destroy */
	for(ni=0; ni<nArrayNum; ni++)
	{
		DeleteString(vCELPath[ni]);
		DeleteString(vAlias[ni]);
		DestroyDoubleMatrix(vArray[ni]);
	}
	free(vCELPath);
	free(vAlias);
	free(vArray);
	free(vSD);

	printf("Done!\n");

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  TileMap_LoadCEL()                                                      */
/*  TileMap loading raw data from a single affymetrix's *.CEL file         */
/* ----------------------------------------------------------------------- */ 
int TileMap_LoadCEL(char strWorkPath[], char strCelFile[], 
					int *pTotalX, int *pTotalY,
					struct DOUBLEMATRIX **pMean, struct DOUBLEMATRIX **pSD)
{
	/* define */
	FILE *fpCel;
		
	char strLine[MED_LINE_LENGTH];
	double dVersion;
	int nCols, nRows;
	int nTotalX,nTotalY;
	int nOffsetX,nOffsetY;
	int nULx,nURx,nLRx,nLLx,nULy,nURy,nLRy,nLLy;
	int nInvertX, nInvertY;
	int swapXY;
	char strDatHeader[MED_LINE_LENGTH];
	char strAlgorithm[MED_LINE_LENGTH];
	char strAlgorithmParameters[MED_LINE_LENGTH];
	int nNumberCells;

	char *chSep;

	int ni,nX,nY,nP,nidx;
	double dM,dS;

	/* init */
	nTotalX = 0;
	nTotalY = 0;

	/* load *.CEL */
	sprintf(strLine, "%s%s", strWorkPath, strCelFile);
	fpCel = NULL;
	fpCel = fopen(strLine, "r");
	if(fpCel == NULL)
	{
		printf("Error: TileMap_LoadCEL, cannot open *.CEL file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, MED_LINE_LENGTH, fpCel)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strstr(strLine, "[CEL]") == strLine)
		{
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			chSep = strchr(strLine, '=');
			chSep++;
			dVersion = atof(chSep);
		}
		else if(strstr(strLine, "[HEADER]") == strLine)
		{
			/* cols */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "Cols") != strLine)
			{
				printf("Error: TileMap_LoadCEL, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			nCols = atoi(chSep);

			/* rows */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "Rows") != strLine)
			{
				printf("Error: TileMap_LoadCEL, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			nRows = atoi(chSep);

			/* total x */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "TotalX") != strLine)
			{
				printf("Error: TileMap_LoadCEL, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			nTotalX = atoi(chSep);

			/* total y */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "TotalY") != strLine)
			{
				printf("Error: TileMap_LoadCEL, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			nTotalY = atoi(chSep);

			/* offset x */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "OffsetX") != strLine)
			{
				printf("Error: TileMap_LoadCEL, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			nOffsetX = atoi(chSep);

			/* offset y */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "OffsetY") != strLine)
			{
				printf("Error: TileMap_LoadCEL, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			nOffsetY = atoi(chSep);

			/* UL */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "GridCornerUL") != strLine)
			{
				printf("Error: TileMap_LoadCEL, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			sscanf(chSep, "%d %d", &nULx, &nULy);

			/* UR */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "GridCornerUR") != strLine)
			{
				printf("Error: TileMap_LoadCEL, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			sscanf(chSep, "%d %d", &nURx, &nURy);

			/* LR */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "GridCornerLR") != strLine)
			{
				printf("Error: TileMap_LoadCEL, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			sscanf(chSep, "%d %d", &nLRx, &nLRy);

			/* LL */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "GridCornerLL") != strLine)
			{
				printf("Error: TileMap_LoadCEL, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			sscanf(chSep, "%d %d", &nLLx, &nLLy);

			/* invertx */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "Axis") != strLine)
			{
				printf("Error: TileMap_LoadCEL, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			nInvertX = atoi(chSep);

			/* inverty */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "Axis") != strLine)
			{
				printf("Error: TileMap_LoadCEL, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			nInvertY = atoi(chSep);

			/* swapxy */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "swapXY") != strLine)
			{
				printf("Error: TileMap_LoadCEL, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			swapXY = atoi(chSep);

			/* DatHeader */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "DatHeader") != strLine)
			{
				printf("Error: TileMap_LoadCEL, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			strcpy(strDatHeader, chSep);

			/* Algorithm */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "Algorithm") != strLine)
			{
				printf("Error: TileMap_LoadCEL, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			strcpy(strAlgorithm, chSep);

			/* Algorithm Parameters */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "AlgorithmParameters") != strLine)
			{
				printf("Error: TileMap_LoadCEL, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			strcpy(strAlgorithmParameters, chSep);
		}

		else if(strstr(strLine, "[INTENSITY]") == strLine)
		{
			/* NumberCells */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "NumberCells") != strLine)
			{
				printf("Error: TileMap_LoadCEL, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			nNumberCells = atoi(chSep);

			/* cell header */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "CellHeader") != strLine)
			{
				printf("Error: TileMap_LoadCEL, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}

			/* allocate memory */
			if((int)(nTotalX*nTotalY) != nNumberCells)
			{
				printf("Error: TileMap_LoadCEL, loading *.CEL file incorrectly, cell number not match!\n");
				exit(EXIT_FAILURE);
			}

			*pTotalX = nTotalX;
			*pTotalY = nTotalY;

			*pMean = NULL;
			*pMean = CreateDoubleMatrix(1, nNumberCells);
			if(*pMean == NULL)
			{
				printf("Error: TileMap_LoadCEL, cannot allocate enough memory for loading *.CEL file!\n");
				exit(EXIT_FAILURE);
			}

			*pSD = NULL;
			*pSD = CreateDoubleMatrix(1, nNumberCells);
			if(*pSD == NULL)
			{
				printf("Error: TileMap_LoadCEL, cannot allocate enough memory for loading *.CEL file!\n");
				exit(EXIT_FAILURE);
			}

			/* load intensity */
			for(ni=0; ni<nNumberCells; ni++)
			{
				fgets(strLine, MED_LINE_LENGTH, fpCel);
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				sscanf(strLine, "%d %d %lf %lf %d", &nX, &nY, &dM, &dS, &nP);
				
				nidx = nY*nTotalX+nX;
				(*pMean)->pMatElement[nidx] = dM;
				(*pSD)->pMatElement[nidx] = dS;
			}
		}

		else if(strstr(strLine, "[MASKS]") == strLine)
		{
			/* NumberCells */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "NumberCells") != strLine)
			{
				printf("Error: TileMap_LoadCEL, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			nNumberCells = atoi(chSep);

			/* cell header */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "CellHeader") != strLine)
			{
				printf("Error: TileMap_LoadCEL, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}

			/* load cells */
			for(ni=0; ni<nNumberCells; ni++)
			{
				fgets(strLine, MED_LINE_LENGTH, fpCel);
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
			}
		}

		else if(strstr(strLine, "[OUTLIERS]") == strLine)
		{
			/* NumberCells */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "NumberCells") != strLine)
			{
				printf("Error: TileMap_LoadCEL, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			nNumberCells = atoi(chSep);

			/* cell header */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "CellHeader") != strLine)
			{
				printf("Error: TileMap_LoadCEL, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}

			/* load cells */
			for(ni=0; ni<nNumberCells; ni++)
			{
				fgets(strLine, MED_LINE_LENGTH, fpCel);
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				sscanf(strLine, "%d %d", &nX, &nY);
			}
		}

		else if(strstr(strLine, "[MODIFIED]") == strLine)
		{
			/* NumberCells */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "NumberCells") != strLine)
			{
				printf("Error: TileMap_LoadCEL, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			nNumberCells = atoi(chSep);

			/* cell header */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "CellHeader") != strLine)
			{
				printf("Error: TileMap_LoadCEL, loading *.CEL file incorrectly!\n");
				exit(EXIT_FAILURE);
			}

			/* load cells */
			for(ni=0; ni<nNumberCells; ni++)
			{
				fgets(strLine, MED_LINE_LENGTH, fpCel);
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
			}
		}

		else
		{
		}
	}

	fclose(fpCel);

	/* return */
	return PROC_SUCCESS;
}

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
								double dLowerBound, int nLogTransform)
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	FILE *fpExport;

	/* variables */
	char strFileType[9];
	float fVersion;
	unsigned long nSeqNum;

	/* seq info */
	struct tagString **vSeqName;
	unsigned long nSeqNameLen;

	struct INTMATRIX *vProbeMappingType;
	unsigned long nProbeMappingType;
	
	struct INTMATRIX *vSequenceFileOffset;
	unsigned long nSequenceFileOffset;
	
	struct INTMATRIX *vProbePairNum;
	unsigned long nProbePairNum;
	
	struct tagString **vGroupName;
	struct tagString **vVersionName;

	struct INTMATRIX *vParamNum;
	unsigned long nParamNum;
	struct tagString ***vParamName;
	struct tagString ***vParamValue;
	
	
	/* map info */
	struct INTMATRIX *vSeqID;
	unsigned long nSeqID;
	int nProbeNum;
	struct tagAffyBpMapUnit *pNewUnit,*pCUnit,*pPUnit;
	struct tagAffyBpMapUnit *pUnitList;

	/* count */
	int ni,nj,nk,nidx1,nidx2;
	int nEndPos;
	int nIgnore;
	double dPM,dMM;
	double dLog2 = log(2.0);


	/* load */
	fpIn = NULL;
	fpIn = fopen(strBpmapPath, "rb");
	if(fpIn == NULL)
	{
		printf("Error: TileMap_ExportAffyIntensity, cannot open .bpmap file!\n");
		exit(EXIT_FAILURE);
	}
	fpOut = NULL;
	fpOut = fopen(strMaskPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: TileMap_ExportAffyIntensity, cannot open .refmask file!\n");
		exit(EXIT_FAILURE);
	}
	fprintf(fpOut, "chromosome\tposition\tprobe_num\trepeat_num\tprobe_seq\n");
	fpExport = NULL;
	fpExport = fopen(strOutPath, "w");
	if(fpExport == NULL)
	{
		printf("Error: TileMap_ExportAffyIntensity, cannot open export file!\n");
		exit(EXIT_FAILURE);
	}
	fprintf(fpExport, "chromosome\tposition");
	for(nk=0; nk<nArrayNum; nk++)
	{
		fprintf(fpExport, "\t%s", vAlias[nk]->m_pString);
	}
	fprintf(fpExport, "\n");

	/* load head */
	fread( strFileType, sizeof(char), 8, fpIn );
	strFileType[8] = '\0';
	AFFYBAR_READ_VERSION(fpIn, &fVersion);
	AFFYBAR_READ_ULONG(fpIn, &nSeqNum);
	printf("BPMAP Version = %f\n", fVersion);
	printf("Number of Sequences = %d\n", nSeqNum); 

	/* load sequence names */
	vSeqName = NULL;
	vSeqName = (struct tagString **)calloc(nSeqNum, sizeof(struct tagString *));
	if(vSeqName == NULL)
	{
		printf("Error: TileMap_ExportAffyIntensity, cannot load sequence name!\n");
		exit(EXIT_FAILURE);
	}

	vProbeMappingType = NULL;
	vProbeMappingType = CreateIntMatrix(nSeqNum,1);
	if(vProbeMappingType == NULL)
	{
		printf("Error: TileMap_ExportAffyIntensity, cannot load probe mapping type!\n");
		exit(EXIT_FAILURE);
	}

	if(fVersion > 2.5)
	{
		vSequenceFileOffset = NULL;
		vSequenceFileOffset = CreateIntMatrix(nSeqNum,1);
		if(vSequenceFileOffset == NULL)
		{
			printf("Error: TileMap_ExportAffyIntensity, cannot load sequence file offset!\n");
			exit(EXIT_FAILURE);
		}
	}

	vProbePairNum = NULL;
	vProbePairNum = CreateIntMatrix(nSeqNum,1);
	if(vProbePairNum == NULL)
	{
		printf("Error: TileMap_ExportAffyIntensity, cannot load probe pair number!\n");
		exit(EXIT_FAILURE);
	}
	
	if(fVersion > 1.5)
	{
		vGroupName = NULL;
		vGroupName = (struct tagString **)calloc(nSeqNum, sizeof(struct tagString *));
		if(vGroupName == NULL)
		{
			printf("Error: TileMap_ExportAffyIntensity, cannot load group name!\n");
			exit(EXIT_FAILURE);
		}
		
		vVersionName = NULL;
		vVersionName = (struct tagString **)calloc(nSeqNum, sizeof(struct tagString *));
		if(vVersionName == NULL)
		{
			printf("Error: TileMap_ExportAffyIntensity, cannot load version name!\n");
			exit(EXIT_FAILURE);
		}

		vParamNum = NULL;
		vParamNum = CreateIntMatrix(nSeqNum, 1);
		if(vParamNum == NULL)
		{
			printf("Error: TileMap_ExportAffyIntensity, cannot allocate memory for loading paramter number!\n");
			exit(EXIT_FAILURE);
		}

		vParamName = NULL;
		vParamName = (struct tagString ***)calloc(nSeqNum, sizeof(struct tagString **));
		if(vParamName == NULL)
		{
			printf("Error: TileMap_ExportAffyIntensity, cannot load sequence parameter!\n");
			exit(EXIT_FAILURE);
		}

		vParamValue = NULL;
		vParamValue = (struct tagString ***)calloc(nSeqNum, sizeof(struct tagString **));
		if(vParamValue == NULL)
		{
			printf("Error: TileMap_ExportAffyIntensity, cannot load sequence parameter!\n");
			exit(EXIT_FAILURE);
		}
	}

	*pnTotalProbeNum = 0;
	for(ni=0; ni<(int)nSeqNum; ni++)
	{
		/* for version 1.0 or later */
		AFFYBAR_READ_ULONG(fpIn, &nSeqNameLen);
		vSeqName[ni] = CreateString(nSeqNameLen);
		fread( vSeqName[ni]->m_pString, sizeof(char), nSeqNameLen, fpIn );
		vSeqName[ni]->m_pString[nSeqNameLen] = '\0';
		
		/* for version 3.0 or later */
		if(fVersion > 2.5)
		{
			AFFYBAR_READ_ULONG(fpIn,&nProbeMappingType);
			vProbeMappingType->pMatElement[ni] = (int)nProbeMappingType;
			AFFYBAR_READ_ULONG(fpIn, &nSequenceFileOffset);
			vSequenceFileOffset->pMatElement[ni] = (int)nSequenceFileOffset;
		}

		/* for version 1.0 or later */
		AFFYBAR_READ_ULONG(fpIn, &nProbePairNum);
		vProbePairNum->pMatElement[ni] = (int)nProbePairNum;
		*pnTotalProbeNum += (int)nProbePairNum;

		/* for version 2.0 or later */
		if(fVersion > 1.5)
		{
			/* read group name */
			AFFYBAR_READ_ULONG(fpIn, &nSeqNameLen);
			vGroupName[ni] = CreateString(nSeqNameLen);
			fread( vGroupName[ni]->m_pString, sizeof(char), nSeqNameLen, fpIn );
			vGroupName[ni]->m_pString[nSeqNameLen] = '\0';

			/* read version name */
			AFFYBAR_READ_ULONG(fpIn, &nSeqNameLen);
			vVersionName[ni] = CreateString(nSeqNameLen);
			fread( vVersionName[ni]->m_pString, sizeof(char), nSeqNameLen, fpIn );
			vVersionName[ni]->m_pString[nSeqNameLen] = '\0';

			/* read paramters */
			AFFYBAR_READ_ULONG(fpIn, &nParamNum);
			vParamNum->pMatElement[ni] = nParamNum;
			if(nParamNum > 0)
			{
				vParamName[ni] = (struct tagString **)calloc(nParamNum, sizeof(struct tagString *));
				if(vParamName[ni] == NULL)
				{
					printf("Error: TileMap_ExportAffyIntensity, cannot load sequence parameter!\n");
					exit(EXIT_FAILURE);
				}
				vParamValue[ni] = (struct tagString **)calloc(nParamNum, sizeof(struct tagString *));
				if(vParamValue[ni] == NULL)
				{
					printf("Error: TileMap_ExportAffyIntensity, cannot load sequence parameter!\n");
					exit(EXIT_FAILURE);
				}
			}
			else
			{
				vParamName[ni] = NULL;
				vParamValue[ni] = NULL;
			}

			for(nj=0; nj<nParamNum; nj++)
			{
				AFFYBAR_READ_ULONG(fpIn, &nSeqNameLen);
				vParamName[ni][nj] = CreateString(nSeqNameLen);
				fread( vParamName[ni][nj]->m_pString, sizeof(char), nSeqNameLen, fpIn );
				vParamName[ni][nj]->m_pString[nSeqNameLen] = '\0';

				AFFYBAR_READ_ULONG(fpIn, &nSeqNameLen);
				vParamValue[ni][nj] = CreateString(nSeqNameLen);
				fread( vParamValue[ni][nj]->m_pString, sizeof(char), nSeqNameLen, fpIn );
				vParamValue[ni][nj]->m_pString[nSeqNameLen] = '\0';
			}
		}
	}

	/* load map */
	vSeqID = NULL;
	vSeqID = CreateIntMatrix(nSeqNum,1);
	if(vSeqID == NULL)
	{
		printf("Error: TileMap_ExportAffyIntensity, cannot load sequence id!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<(int)nSeqNum; ni++)
	{
		/* load seq id */
		AFFYBAR_READ_ULONG(fpIn, &nSeqID);
		nProbeNum = vProbePairNum->pMatElement[ni];
		pUnitList = NULL;
		
		/* load probes */
		for(nj=0; nj<nProbeNum; nj++)
		{
			/* load new unit */
			pNewUnit = NULL;
			pNewUnit = AffyBpMapUnitCreate();
			if( fVersion > 2.5)
			{
				AffyBpMapUnitLoad_v3(pNewUnit, fpIn, vProbeMappingType->pMatElement[ni]);
			}
			else
			{
				AffyBpMapUnitLoad(pNewUnit, fpIn);
			}

			/* export intensities */
			nidx1 = pNewUnit->nPMY*nTotalX+pNewUnit->nPMX;
			if(nidx1 >= nCellNum)
			{
				printf("Error: TileMap_ExportAffyIntensity, index out of range\n");
				exit(EXIT_FAILURE);
			}
			nidx2 = pNewUnit->nMMY*nTotalX+pNewUnit->nMMX;
			if(nidx2 >= nCellNum)
			{
				printf("Error: TileMap_ExportAffyIntensity, index out of range\n");
				exit(EXIT_FAILURE);
			}

			fprintf(fpExport, "%s\t%d", vSeqName[ni]->m_pString, pNewUnit->nPos);
			for(nk=0; nk<nArrayNum; nk++)
			{
				dPM = vArray[nk]->pMatElement[nidx1];
				if(vProbeMappingType->pMatElement[ni] == 0)
				{
					dMM = vArray[nk]->pMatElement[nidx2];
				}
				else
				{
					dMM = 0.0;
				}

				if(nIntensityType == 1)
				{
					dPM -= dMM;
				}
				
				if(dPM < dLowerBound)
					dPM = dLowerBound;

				if(nLogTransform == 1)
				{
					dPM = log(dPM)/dLog2;
				}

				fprintf(fpExport, "\t%f", dPM);
			}
			fprintf(fpExport, "\n");

			/* add to unitlist */
			if(pUnitList == NULL)
			{
				pUnitList = pNewUnit;
			}
			else
			{
				nIgnore = 0;
				pCUnit = pUnitList;
				while(pCUnit != NULL)
				{
					/* if same position */
					if(pNewUnit->nPos == pCUnit->nPos)
					{
						pCUnit->nDepthNum += 1;
						nIgnore = 1;
					}
					pCUnit = pCUnit->pNext;
				}

				if(nIgnore == 1)
				{
					AffyBpMapUnitDestroy(pNewUnit);
				}
				else
				{
					pCUnit = pUnitList;
					pPUnit = NULL;
					while(pCUnit != NULL)
					{
						/* if local repeat */
						if(pNewUnit->nPos-pCUnit->nPos <= LOCAL_REPEAT_RANGE)
						{
							if(strcmp(pCUnit->strProbeSeq, pNewUnit->strProbeSeq) == 0)
							{
								pCUnit->nRepeatNum += 1;
								pNewUnit->nRepeatNum += 1;
							}
						}

						/* get next */
						pPUnit = pCUnit;
						pCUnit = pPUnit->pNext;
					}

					pPUnit->pNext = pNewUnit;
					nEndPos = pNewUnit->nPos;

					/* write old unit */
					while(pUnitList != NULL)
					{
						/* if local repeat */
						if(nEndPos-pUnitList->nPos > LOCAL_REPEAT_RANGE)
						{
							pCUnit = pUnitList;
							pUnitList = pCUnit->pNext;
							pCUnit->pNext = NULL;

							fprintf(fpOut, "%s\t%d\t%d\t%d\t%s\n", vSeqName[ni]->m_pString, pCUnit->nPos, pCUnit->nDepthNum, pCUnit->nRepeatNum, pCUnit->strProbeSeq);
							AffyBpMapUnitDestroy(pCUnit);
						}
						else
						{
							break;
						}
					}
				}
			}			
		}

		/* clear unit list */
		while(pUnitList != NULL)
		{
			pCUnit = pUnitList;
			pUnitList = pCUnit->pNext;
			pCUnit->pNext = NULL;

			fprintf(fpOut, "%s\t%d\t%d\t%d\t%s\n", vSeqName[ni]->m_pString, pCUnit->nPos, pCUnit->nDepthNum, pCUnit->nRepeatNum, pCUnit->strProbeSeq);
			AffyBpMapUnitDestroy(pCUnit);
		}
	}

	/* load tail if any */
	while(feof(fpIn) == 0)
	{
		fread(strFileType, sizeof(char), 1, fpIn);
		printf("%c", strFileType[0]);
	}

	/* clear memeory */
	for(ni=0; ni<(int)nSeqNum; ni++)
	{
		DeleteString(vSeqName[ni]);
		vSeqName[ni] = NULL;
	}
	free(vSeqName);
	DestroyIntMatrix(vProbeMappingType);
	if(fVersion > 2.5)
	{
		DestroyIntMatrix(vSequenceFileOffset);
	}
	DestroyIntMatrix(vProbePairNum);
	if(fVersion > 1.5)
	{
		for(ni=0; ni<(int)nSeqNum; ni++)
		{
			DeleteString(vGroupName[ni]);
			vGroupName[ni] = NULL;

			DeleteString(vVersionName[ni]);
			vVersionName[ni] = NULL;
		
			if(vParamNum->pMatElement[ni] > 0)
			{
				for(nj=0; nj<vParamNum->pMatElement[ni]; nj++)
				{
					DeleteString(vParamName[ni][nj]);
					vParamName[ni][nj] = NULL;
					DeleteString(vParamValue[ni][nj]);
					vParamValue[ni][nj] = NULL;
				}
				free(vParamName[ni]);
				vParamName[ni] = NULL;
				free(vParamValue[ni]);
				vParamValue[ni] = NULL;
			}
		}
		free(vGroupName);
		free(vVersionName);
		free(vParamName);
		free(vParamValue);
		DestroyIntMatrix(vParamNum);
	}
	DestroyIntMatrix(vSeqID);

	/* close file */
	fclose(fpIn);
	fclose(fpOut);
	fclose(fpExport);

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  TileMap_Normalization_Main()                                           */
/*  TileMap normalization module.                                          */
/* ----------------------------------------------------------------------- */ 
int TileMap_Normalization_Main(char strParamPath[])
{
	/* define */

	/* working path */
	char strWorkPath[MED_LINE_LENGTH];
	char strInFile[MED_LINE_LENGTH];
	char strOutFile[MED_LINE_LENGTH];
	char strDataPath[MED_LINE_LENGTH];
	char strExportPath[MED_LINE_LENGTH];
	
	/* arrays */
	int nArrayNum = 0;
	int nProbeNum = 0;
	
	/* normalization */
	double dNormLowerBound = 0.0;
	int nNormLogTransform = 0;
	
	/* others */
	char strLine[MED_LINE_LENGTH];
	char strDataLine[LONG_LINE_LENGTH];
	FILE *fpIn;
	
	int nlen;
	char *chSep;
	int nError = 0;

	/* load parameters */
	fpIn = NULL;
	fpIn = fopen(strParamPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileMap_Normalization_Main, cannot open the parameter file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strstr(strLine, "[Working directory]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strWorkPath, chSep);
			
			if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
			{
				nlen = strlen(strWorkPath);
				if(strWorkPath[nlen-1] != '\\')
				{
					strWorkPath[nlen] = '\\';
					strWorkPath[nlen+1] = '\0';
				}
			}
			else
			{
				nlen = strlen(strWorkPath);
				if(strWorkPath[nlen-1] != '/')
				{
					strWorkPath[nlen] = '/';
					strWorkPath[nlen+1] = '\0';
				}
			}
		}

		else if(strstr(strLine, "[Raw Data file]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strInFile, chSep);
		}
		
		else if(strstr(strLine, "[Export file]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strOutFile, chSep);
		}

		else if(strstr(strLine, "[Array number]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			nArrayNum = atoi(chSep);
			if(nArrayNum <= 0)
			{
				printf("Error: TileMap_Normalization_Main, no arrays available!\n");
				nError = 1;
				break;
			}
		}
		
		else if(strstr(strLine, "[Truncation lower bound before normalization]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			dNormLowerBound = atof(chSep);
		}
		
		else if(strstr(strLine, "[Take log2 transformation before normalization]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			nNormLogTransform = atoi(chSep);
		}
	
		else if(strLine[0] == '#')
		{
		}

		else
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	fclose(fpIn);

	if(nError > 0)
	{
		printf("Error: TileMap_Normalization_Main, cannot load the parameter file correctly!\n");
		exit(EXIT_FAILURE);
	}

	/* get probe number */
	sprintf(strDataPath, "%s%s", strWorkPath, strInFile);
	fpIn = NULL;
	fpIn = fopen(strDataPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileMap_Normalization_Main, cannot open raw data file!\n");
		exit(EXIT_FAILURE);
	}

	fgets(strDataLine, LONG_LINE_LENGTH, fpIn);
	while(fgets(strDataLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strDataLine);
		StrTrimRight(strDataLine);
		if(strDataLine[0] == '\0')
			continue;

		nProbeNum++;
	}

	fclose(fpIn);

	printf("Array Number = %d\n", nArrayNum);
	printf("Probe Number = %d\n", nProbeNum);
	printf("Normalization...\n");

	/* normalization */
	sprintf(strExportPath, "%s%s", strWorkPath, strOutFile);
	TileMap_Normalization_Quantile(nArrayNum, nProbeNum, 
								strDataPath, strExportPath,
								nNormLogTransform, dNormLowerBound);
	
	printf("Done!\n");

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  TileMap_Normalization_Quantile()                                       */
/*  Quantile normalization for tilemap.                                    */
/*  Return PROC_SUCCESS if success.                                        */
/* ----------------------------------------------------------------------- */ 
int TileMap_Normalization_Quantile(int nArrayNum, int nProbeNum, 
								char strDataFile[], char strOutFile[],
								int nTakeLog, double dTruncLow)
{
	/* define */
	struct tagString **vChr;
	struct INTMATRIX *pPosition;
	struct DOUBLEMATRIX **vArray;
	struct DOUBLEMATRIX **vSortArray;
	struct LONGMATRIX **vSortIndex;
	FILE *fpData;
	FILE *fpOut;
	int ni,nj,nk;
	char strLongLine[LONG_LINE_LENGTH];
	char *chSep,*chSep2;
	double dValue;
	double dLog2 = log(2.0);
	

	/* check */
	if( (nArrayNum <= 0) || (nProbeNum <= 0) )
	{
		printf("Warning: TileMap_Normalization_Quantile, array or probe number <= 0!\n");
		return PROC_SUCCESS;
	}

	/* init */
	pPosition = NULL;
	pPosition = CreateIntMatrix(1, nProbeNum);
	if(pPosition == NULL)
	{
		printf("Error: TileMap_Normalization_Quantile, cannot create position matrix!\n");
		exit(EXIT_FAILURE);
	}

	vChr = NULL;
	vChr = (struct tagString **)calloc(nProbeNum, sizeof(struct tagString *));
	if(vChr == NULL)
	{
		printf("Error: TileMap_Normalization_Quantile, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}

	vArray = NULL;
	vArray = (struct DOUBLEMATRIX **)calloc(nArrayNum, sizeof(struct DOUBLEMATRIX *));
	if(vArray == NULL)
	{
		printf("Error: TileMap_Normalization_Quantile, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nArrayNum; ni++)
	{
		vArray[ni] = CreateDoubleMatrix(1, nProbeNum);
		if(vArray[ni] == NULL)
		{
			printf("Error: TileMap_Normalization_Quantile, cannot create data space!\n");
			exit(EXIT_FAILURE);
		}
	}

	vSortArray = NULL;
	vSortArray = (struct DOUBLEMATRIX **)calloc(nArrayNum, sizeof(struct DOUBLEMATRIX *));
	if(vSortArray == NULL)
	{
		printf("Error: TileMap_Normalization_Quantile, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}
	vSortIndex = NULL;
	vSortIndex = (struct LONGMATRIX **)calloc(nArrayNum, sizeof(struct LONGMATRIX *));
	if(vSortIndex == NULL)
	{
		printf("Error: TileMap_Normalization_Quantile, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}

	/* data */
	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: TileMap_Normalization_Quantile, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fpData = NULL;
	fpData = fopen(strDataFile, "r");
	if(fpData == NULL)
	{
		printf("Error: TileMap_Normalization_Quantile, cannot open data file!\n");
		exit(EXIT_FAILURE);
	}
	
	fgets(strLongLine, LONG_LINE_LENGTH, fpData);
	fprintf(fpOut, "%s", strLongLine);

	nj = 0;
	while(fgets(strLongLine, LONG_LINE_LENGTH, fpData) != NULL)
	{
		StrTrimLeft(strLongLine);
		StrTrimRight(strLongLine);
		if(strLongLine[0] == '\0')
			continue;

		if(nj >= nProbeNum)
		{
			printf("Error: Normalization_Quantile_Main, probeset number not match!\n");
			exit(EXIT_FAILURE);
		}

		ni = 0;
		
		/* chromosome */
		chSep = strchr(strLongLine, '\t');
		if(chSep != NULL)
			*chSep = '\0';
		StringAddTail((vChr+nj), strLongLine);
		
		/* position */
		chSep++;
		chSep2 = strchr(chSep, '\t');
		if(chSep2 != NULL)
			*chSep2 = '\0';
		pPosition->pMatElement[nj] = (int)(atof(chSep));
		chSep = chSep2;

		while(chSep != NULL)
		{
			if( ni >= nArrayNum )
			{
				printf("Error: Normalization_Quantile_Main, array number not match!\n");
				exit(EXIT_FAILURE);
			}

			chSep++;
			chSep2 = strchr(chSep, '\t');
		
			/* middle number */
			if(chSep2 != NULL)
			{
				*chSep2 = '\0';

				if(chSep == chSep2)
				{
					dValue = 0.0;
				}
				else
				{
					dValue = atof(chSep);
				}
			}
			/* last number */
			else
			{
				if(chSep == chSep2)
				{
					dValue = 0.0;
				}
				else
				{
					dValue = atof(chSep);
				}
			}

			if(dValue < dTruncLow)
				dValue = dTruncLow;
			if(nTakeLog == 1)
				dValue = log(dValue)/dLog2;
		
			DMSETAT(vArray[ni], 0, nj, dValue);
			ni++;

			/* get next */
			chSep = chSep2;
		}

		if(ni!=nArrayNum)
		{
			printf("Error: Normalization_Quantile_Main, array number not match!\n");
			exit(EXIT_FAILURE);
		}

		/* get next line */
		nj++;
	}

	fclose(fpData);

	if(nj != nProbeNum)
	{
		printf("Error: Normalization_Quantile_Main, probe number not match!\n");
		exit(EXIT_FAILURE);
	}


	/* normalization */
	for(ni=0; ni<nArrayNum; ni++)
	{
		DMSORTMERGEA_0(vArray[ni], (vSortArray+ni), (vSortIndex+ni));
	}
	for(nj=0; nj<nProbeNum; nj++)
	{
		dValue = 0.0;
		for(ni=0; ni<nArrayNum; ni++)
		{
			dValue += vSortArray[ni]->pMatElement[nj];
		}
		dValue /= (double)nArrayNum;
		for(ni=0; ni<nArrayNum; ni++)
		{
			vSortArray[ni]->pMatElement[nj] = dValue;
		}
	}

	/* save */
	for(nj=0; nj<nProbeNum; nj++)
	{
		for(ni=0; ni<nArrayNum; ni++)
		{
			nk = (int)(vSortIndex[ni]->pMatElement[nj]);
			vArray[ni]->pMatElement[nk] = vSortArray[ni]->pMatElement[nj];
		}
	}

	for(nj=0; nj<nProbeNum; nj++)
	{
		fprintf(fpOut, "%s\t%d", vChr[nj]->m_pString, pPosition->pMatElement[nj]);
		for(ni=0; ni<nArrayNum; ni++)
		{
			fprintf(fpOut, "\t%9.7e", vArray[ni]->pMatElement[nj]);
		}
		fprintf(fpOut, "\n");
	}


	/* release memory */
	fclose(fpOut);
	DestroyIntMatrix(pPosition);
	for(ni=0; ni<nArrayNum; ni++)
	{
		DestroyDoubleMatrix(vArray[ni]);
		vArray[ni] = NULL;
		DestroyDoubleMatrix(vSortArray[ni]);
		vSortArray[ni] = NULL;
		DestroyLongMatrix(vSortIndex[ni]);
		vSortIndex[ni] = NULL;
	}
	free(vArray);
	free(vSortArray);
	free(vSortIndex);

	for(ni=0; ni<nProbeNum; ni++)
	{
		DeleteString(vChr[ni]);
		vChr[ni] = NULL;
	}
	free(vChr);

	/* return */
	return PROC_SUCCESS;
}



/* ----------------------------------------------------------------------- */ 
/*  TileMap_Main()                                                         */
/*  TileMap pipeline                                                       */
/* ----------------------------------------------------------------------- */ 
int TileMap_Main(char strParamPath[])
{
	/* ------ */
	/* define */
	/* ------ */

	/* working directory */
	char strWorkPath[MED_LINE_LENGTH];
	char strProjectTitle[LINE_LENGTH];
	int nTotalProbeNum = 0;
	int nRealProbeNum = 0;
	int nPermNum = 0;
	int nApplyPerm = 0;

	/* probe summary */
	int nApplyProbeSummary = 1;
	char strRawDataFile[LINE_LENGTH];
	int nProbeSummaryRange = 0;
	double dZeroCut = 1e-6;

	/* repeat filter */
	int nApplyFilter = 0;
	char strRefFilterFile[LINE_LENGTH];

	/* region summary */
	int nApplyRegionSummary = 1;
	int nRegionSummaryType = 0;

	/* UMS */
	char strUMSSelectPath[MED_LINE_LENGTH];
	char strUMSScorePath[MED_LINE_LENGTH];

	/* HMM */
	double dHMMPostCut = 0.5;
	int nHMMMaxGap = 1000;
	int nHMMParameterType = 0;

	int nHMMProvideSelectionV = 0;
	char strHMMSelectPath[MED_LINE_LENGTH];
	double dHMMTp = 0.01;
	double dHMMTq = 0.05;
	int nHMMoff = 1;
	int nHMMGrid = 1000;
	int nHMMFragLen = 28;

	char strHMMTransitionFile[LINE_LENGTH];
	char strHMMEmissionFile[LINE_LENGTH];

	/* MA */
	double dMAFDRCut = 0.5;
	int nMAMaxGap = 500;
	int nMAW = 5;
	int nMAParameterType = 0;

	int nMAProvideSelectionV = 0;
	char strMASelectPath[MED_LINE_LENGTH];
	double dMATp = 0.01;
	double dMATq = 0.05;
	int nMAoff = 6;
	int nMAGrid = 1000;

	/* output files */
	char strDataPath[MED_LINE_LENGTH];
	char strRefMaskPath[MED_LINE_LENGTH];
	char strCmpInfoPath[MED_LINE_LENGTH];
	char strPBsumPath[MED_LINE_LENGTH];
	char strMAsumPath[MED_LINE_LENGTH];
	char strHMMsumPath[MED_LINE_LENGTH];
	char strCommand[MED_LINE_LENGTH];

	
	/* for loading parameter settings */
	FILE *fpIn;
	char strLine[MED_LINE_LENGTH];
	char strDataLine[LONG_LINE_LENGTH];
	char strTransformType[LINE_LENGTH];
	int nlen;
	char *chSep;
	int nError = 0;
	int ni;

	/* --------------- */
	/* load parameters */
	/* --------------- */
	fpIn = NULL;
	fpIn = fopen(strParamPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileMap_Main, cannot open the parameter file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strstr(strLine, "O.1-[Working directory]") == strLine)
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
			
			if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
			{
				nlen = strlen(strWorkPath);
				if(nlen == 0)
				{
					sprintf(strWorkPath, ".\\"); 
				}
				else if(strWorkPath[nlen-1] != '\\')
				{
					strWorkPath[nlen] = '\\';
					strWorkPath[nlen+1] = '\0';
				}
			}
			else
			{
				nlen = strlen(strWorkPath);
				if(nlen == 0)
				{
					sprintf(strWorkPath, "./"); 
				}
				else if(strWorkPath[nlen-1] != '/')
				{
					strWorkPath[nlen] = '/';
					strWorkPath[nlen+1] = '\0';
				}
			}
		}

		else if(strstr(strLine, "O.2-[Project Title]") == strLine)
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
				strcpy(strProjectTitle, "unnamed_project");
			}
			else
			{
				strcpy(strProjectTitle, chSep);
			}
		}

		else if(strstr(strLine, "I.1-[Compute probe level test-statistics?]") == strLine)
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
				nApplyProbeSummary = atoi(chSep);
			}
		}


		else if(strstr(strLine, "I.2-[Raw data file]") == strLine)
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
				printf("Error: TileMap_Main, you have to specify the data file in I.2!\n");
				nError = 1;
				break;
			}
			else
			{
				strcpy(strRawDataFile, chSep);
			}
		}
		
		else if(strstr(strLine, "I.3-[Range of test-statistics]") == strLine)
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
				nProbeSummaryRange = atoi(chSep);
			}
		}

		else if(strstr(strLine, "I.4-[Zero cut]") == strLine)
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
				dZeroCut = atof(chSep);
			}
		}

		else if(strstr(strLine, "II.1-[Apply local repeat filter?]") == strLine)
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
				nApplyFilter = atoi(chSep);
			}
		}

		else if(strstr(strLine, "II.2-[*.refmask file]") == strLine)
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
				strcpy(strRefFilterFile, "NULL");
			}
			else
			{
				strcpy(strRefFilterFile, chSep);
			}
		}

		else if(strstr(strLine, "III.1-[Combine neighboring probes?]") == strLine)
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
				nApplyRegionSummary = atoi(chSep);
			}
		}

		else if(strstr(strLine, "III.2-[Method to combine neighboring probes]") == strLine)
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
				nRegionSummaryType = atoi(chSep);
			}
		}

		else if(strstr(strLine, "IV.1-[Posterior probability >]") == strLine)
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
				dHMMPostCut = atof(chSep);
			}
		}

		else if(strstr(strLine, "IV.2-[Maximal gap allowed]") == strLine)
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
				nHMMMaxGap = atoi(chSep);
			}
		}

		else if(strstr(strLine, "IV.3-[Method to set HMM parameters]") == strLine)
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
				nHMMParameterType = atoi(chSep);
			}
		}

		else if(strstr(strLine, "IV.4-[Provide your own selection statistics?]") == strLine)
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
				nHMMProvideSelectionV = atoi(chSep);
			}
		}

		else if(strstr(strLine, "IV.5-[If Yes to IV.4, selection statistics file]") == strLine)
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
				strcpy(strHMMSelectPath, chSep);
			}
			else
			{
				strcpy(strHMMSelectPath, "NULL");
			}
		}

		else if(strstr(strLine, "IV.6-[G0 Selection Criteria, p%]") == strLine)
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
				dHMMTp = atof(chSep);
			}
		}

		else if(strstr(strLine, "IV.7-[G1 Selection Criteria, q%]") == strLine)
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
				dHMMTq = atof(chSep);
			}
		}

		else if(strstr(strLine, "IV.8-[Selection Offset]") == strLine)
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
				nHMMoff = atoi(chSep);
			}
		}

		else if(strstr(strLine, "IV.9-[Grid Size]") == strLine)
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
				nHMMGrid = atoi(chSep);
			}
		}

		else if(strstr(strLine, "IV.10-[Expected hybridization length]") == strLine)
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
				nHMMFragLen = atoi(chSep);
			}
		}

		else if(strstr(strLine, "IV.11-[Path to transition probability matrix]") == strLine)
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
				strcpy(strHMMTransitionFile, "NULL");
			}
			else
			{
				strcpy(strHMMTransitionFile, chSep);
			}
		}

		else if(strstr(strLine, "IV.12-[Path to emission probability matrix]") == strLine)
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
				strcpy(strHMMEmissionFile, "NULL");
			}
			else
			{
				strcpy(strHMMEmissionFile, chSep);
			}
		}

		else if(strstr(strLine, "V.1-[Local FDR <]") == strLine)
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
				dMAFDRCut = atof(chSep);
			}
		}

		else if(strstr(strLine, "V.2-[Maximal gap allowed]") == strLine)
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
				nMAMaxGap = atoi(chSep);
			}
		}

		else if(strstr(strLine, "V.3-[W]") == strLine)
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
				nMAW = atoi(chSep);
			}
		}

		else if(strstr(strLine, "V.4-[Method to compute local FDR]") == strLine)
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
				nMAParameterType = atoi(chSep);
			}
		}

		else if(strstr(strLine, "V.5-[Provide your own selection statistics?]") == strLine)
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
				nMAProvideSelectionV = atoi(chSep);
			}
		}

		else if(strstr(strLine, "V.6-[If Yes to IV.4, selection statistics file]") == strLine)
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
				strcpy(strMASelectPath, chSep);
			}
			else
			{
				strcpy(strMASelectPath, "NULL");
			}
		}

		else if(strstr(strLine, "V.7-[G0 Selection Criteria, p%]") == strLine)
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
				dMATp = atof(chSep);
			}
		}

		else if(strstr(strLine, "V.8-[G1 Selection Criteria, q%]") == strLine)
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
				dMATq = atof(chSep);
			}
		}

		else if(strstr(strLine, "V.9-[Selection Offset]") == strLine)
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
				nMAoff = atoi(chSep);
			}
		}
		
		else if(strstr(strLine, "V.10-[Grid Size]") == strLine)
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
				nMAGrid = atoi(chSep);
			}
		}

		else if(strLine[0] == '#')
		{
		}

		else
		{
			printf("Error: TileMap_Main, unknown parameters!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	fclose(fpIn);

	if(nError > 0)
	{
		printf("Error: TileMap_Main, cannot load the parameter file correctly!\n");
		exit(EXIT_FAILURE);
	}

	
	/* -------------------- */
	/* Step I:              */
	/* probe level summary  */
	/* -------------------- */

	/* get probe number */
	printf("######################\n");
	printf("       TileMap        \n");
	printf("######################\n");
	printf("Counting Probes...\n");
	nTotalProbeNum = 0;
	sprintf(strLine, "%s%s", strWorkPath, strRawDataFile);
	fpIn = NULL;
	fpIn = fopen(strLine, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileMap_Main, cannot open raw data file!\n");
		exit(EXIT_FAILURE);
	}

	fgets(strDataLine, LONG_LINE_LENGTH, fpIn);
	while(fgets(strDataLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strDataLine);
		StrTrimRight(strDataLine);
		if(strDataLine[0] == '\0')
			continue;

		nTotalProbeNum++;
	}

	fclose(fpIn);
	printf("Number of probes in the raw data = %d\n", nTotalProbeNum);

	/* import user provided statistics */
	
	/* compute statistics */
	if((nApplyRegionSummary == 1) && (nRegionSummaryType==1) && (nMAParameterType == 1))
	{
		nApplyPerm = 1;
	}
	else
	{
		nApplyPerm = 0;
	}

	if(nApplyProbeSummary == 0)
	{
		/* cannot apply permutation test for MA */
		if(nApplyPerm == 1)
		{
			printf("Error: TileMap_Main, permutation test cannot be done if you provide your own probe level statistics. UMS will be used instead!\n");
			nMAParameterType = 0;
			nApplyPerm = 0;
		}

		/* must specify the range of statistics */
		if(nProbeSummaryRange == 0)
		{
			printf("Error: TileMap_Main, you have to specify the range of the probe level statistics!\n");
			nMAParameterType = 0;
			nApplyPerm = 0;
		}
	}

	/* or compute probe level statistics by tilemap */
	else
	{
		/* prepare path */
		sprintf(strDataPath, "%s%s", strWorkPath, strRawDataFile);
		sprintf(strCmpInfoPath, "%s%s.cmpinfo", strWorkPath, strProjectTitle);
		sprintf(strPBsumPath, "%s%s", strWorkPath, strProjectTitle);
				
		TileMap_ProbeSelection_Main(strDataPath, "NULL",
							   strCmpInfoPath, strPBsumPath,
							   nTotalProbeNum, nApplyPerm, &nPermNum,
							   &nProbeSummaryRange,	dZeroCut,
							   nApplyFilter);
	}

	/* -------------------- */
	/* Step II:             */
	/* repeat filtering     */
	/* -------------------- */
	printf("\n");
	printf("/* ---------------------- */\n");
	printf("/* TileMap Repeat Filter  */\n");
	printf("/* ---------------------- */\n");


	/* apply filtering */
	if(nApplyFilter == 1)
	{
		printf("Filtering...\n");

		if(nApplyProbeSummary == 0)
		{
			sprintf(strDataPath, "%s%s", strWorkPath, strRawDataFile);
			sprintf(strPBsumPath, "%s%s_f_pb.sum", strWorkPath, strProjectTitle);
			if(strcmp(strDataPath, strPBsumPath) == 0)
			{
				printf("Error: TileMap_Main, your raw data file name cannot be {Project Title}_f_pb.sum. Please change the name of your raw data file!\n");
				exit(EXIT_FAILURE);
			}
			sprintf(strRefMaskPath, "%s%s", strWorkPath, strRefFilterFile);
			nRealProbeNum = TileMap_BPMAPFilter_Pbsum(strDataPath, strRefMaskPath, strPBsumPath);
		}
		else
		{
			/* filter score */
			sprintf(strDataPath, "%s%s_pb.sum", strWorkPath, strProjectTitle);
			sprintf(strPBsumPath, "%s%s_f_pb.sum", strWorkPath, strProjectTitle);
			sprintf(strRefMaskPath, "%s%s", strWorkPath, strRefFilterFile);
			nRealProbeNum = TileMap_BPMAPFilter_Pbsum(strDataPath, strRefMaskPath, strPBsumPath);

			if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
			{
				sprintf(strCommand, "del %s", strDataPath);
				system(strCommand);
			}
			else
			{
				sprintf(strCommand, "rm %s", strDataPath);
				system(strCommand);
			}

			/* filter permutations */
			if(nApplyPerm == 1)
			{
				for(ni=0; ni<nPermNum; ni++)
				{
					sprintf(strDataPath, "%s%s_perm%d_pb.sum", strWorkPath, strProjectTitle, ni);
					sprintf(strPBsumPath, "%s%s_perm%d_f_pb.sum", strWorkPath, strProjectTitle, ni);
					sprintf(strRefMaskPath, "%s%s", strWorkPath, strRefFilterFile);
					TileMap_BPMAPFilter_Pbsum(strDataPath, strRefMaskPath, strPBsumPath);

					if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
					{
						sprintf(strCommand, "del %s", strDataPath);
						system(strCommand);
					}
					else
					{
						sprintf(strCommand, "rm %s", strDataPath);
						system(strCommand);
					}
				}
			}
		}
	}
	/* do not apply filtering */
	else
	{
		printf("No Filtering...\n");
		nRealProbeNum = nTotalProbeNum;
	}
	printf("Number of probes after filtering = %d\n", nRealProbeNum);
	

	/* -------------------- */
	/* Step III:            */
	/* region summary       */
	/* -------------------- */
	if(nApplyProbeSummary == 0)
	{
		if(nApplyFilter == 0)
			sprintf(strPBsumPath, "%s%s", strWorkPath, strRawDataFile);
		else
			sprintf(strPBsumPath, "%s%s_f_pb.sum", strWorkPath, strProjectTitle);
	}
	else
	{
		sprintf(strPBsumPath, "%s%s_f_pb.sum", strWorkPath, strProjectTitle);
	}


	if(nApplyRegionSummary == 0)
	{
		/* return */
		return PROC_SUCCESS;
	}

	/* -------------------- */
	/* MA                   */
	/* -------------------- */
	if(nRegionSummaryType == 1)
	{
		/* -------------------- */
		/* MA                   */
		/* -------------------- */
		if(nMAParameterType == 1)
		{
			/* -------------------- */
			/* PERMUTATION Based MA */
			/* -------------------- */
			TileMap_BindingRegionSelection_MA_PERM_Main(nRealProbeNum, 
				strPBsumPath, strWorkPath, strProjectTitle, nMAW,
				nProbeSummaryRange,	dZeroCut, nMAGrid, nPermNum, nMAMaxGap, dMAFDRCut);
		}
		else
		{
			/* -------------------- */
			/* UMS Based MA         */
			/* -------------------- */
			if(nMAProvideSelectionV == 1)
			{
				if(strcmp(strMASelectPath, "NULL") == 0)
				{
					printf("Error: TileMap_Main, you have to specify the file containing selection statistics!\n");
					exit(EXIT_FAILURE);
				}
				sprintf(strUMSSelectPath, "%s%s", strWorkPath, strMASelectPath);
			}
			else
			{
				sprintf(strUMSSelectPath, "%s", strPBsumPath);
			}

			TileMap_BindingRegionSelection_MA_UMS_Main(nRealProbeNum, 
				strPBsumPath, strWorkPath, strProjectTitle, nMAW,
				nProbeSummaryRange,	dZeroCut, strUMSSelectPath, 
				dMATp, dMATq, nMAoff, nMAGrid, nMAMaxGap, dMAFDRCut);
		}
	}

	/* -------------------- */
	/* HMM                  */
	/* -------------------- */
	else
	{
		/* if use UMS to estimate parameters */
		if(nHMMParameterType == 0)
		{
			/* -------------------- */
			/* UMS                  */
			/* -------------------- */
			printf("\n");
			printf("/* ---------------------- */\n");
			printf("/* TileMap UMS            */\n");
			printf("/* ---------------------- */\n");

			if(nHMMProvideSelectionV == 1)
			{
				if(strcmp(strHMMSelectPath, "NULL") == 0)
				{
					printf("Error: TileMap_Main, you have to specify the file containing selection statistics!\n");
					exit(EXIT_FAILURE);
				}
				sprintf(strUMSSelectPath, "%s%s", strWorkPath, strHMMSelectPath);
			}
			else
			{
				sprintf(strUMSSelectPath, "%s", strPBsumPath);
			}

			
			sprintf(strUMSScorePath, "%s", strPBsumPath);
			sprintf(strLine, "%s%s", strWorkPath, strProjectTitle);
			TileMap_UMS_HMM_Main(strUMSSelectPath, strUMSScorePath, 
							nRealProbeNum,  nProbeSummaryRange,
							dHMMTp, dHMMTq, nHMMoff, nHMMGrid,
							nHMMFragLen, strLine);
			sprintf(strHMMTransitionFile, "%s_transp.txt", strLine);
			sprintf(strHMMEmissionFile, "%s_emissp.txt", strLine);
		}
		/* if use user specified parameters */
		else
		{
			if((strcmp(strHMMTransitionFile, "NULL") == 0) || (strcmp(strHMMEmissionFile, "NULL") == 0))
			{
				printf("Error: TileMap_Main, you have to specify the transition and emission probability matrix for HMM!\n");
				exit(EXIT_FAILURE);
			}

			sprintf(strLine, "%s%s", strWorkPath, strHMMTransitionFile);
			strcpy(strHMMTransitionFile, strLine);

			sprintf(strLine, "%s%s", strWorkPath, strHMMEmissionFile);
			strcpy(strHMMEmissionFile, strLine);
		}


		/* -------------------- */
		/* HMM                  */
		/* -------------------- */
		printf("\n");
		printf("/* ---------------------- */\n");
		printf("/* TileMap HMM            */\n");
		printf("/* ---------------------- */\n");
		if(nProbeSummaryRange == 2)
		{
			strcpy(strTransformType, "invlogit");
		}
		else
		{
			strcpy(strTransformType, "identity");
		}
		sprintf(strHMMsumPath, "%s%s", strWorkPath, strProjectTitle);
		TileMap_BindingRegionSelection_HMM_Main(nRealProbeNum, strPBsumPath,
				strTransformType, dZeroCut,
				strHMMTransitionFile, strHMMEmissionFile, 
				nHMMMaxGap, dHMMPostCut,
				strHMMsumPath);
		
		/* -------------------- */
		/* RANK binding regions */
		/* -------------------- */
	}

	/* -------------------- */
	/* Done                 */
	/* -------------------- */
	printf("Done!\n");

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  TileMap_ProbeSelection_Main()                                          */
/*  Probe selection based on criteria specified in the criteria file       */
/*  Return PROC_SUCCESS if success.                                        */
/* ----------------------------------------------------------------------- */ 
int TileMap_ProbeSelection_Main(char strDataFile[], char strGeneInfoFile[],
							   char strCriteriaFile[], char strOutFile[],
							   int nProbeNum, int nApplyPerm, int *nPermNum,
							   int *nProbeSummaryRange,	double dZeroCut,
							   int nApplyFilter)
{
	/* ------- */
	/* define  */
	/* ------- */

	/* array */
	int nArrayNum;
	struct DOUBLEMATRIX **vArray;

	/* class id */
	int nClassNum;
	struct INTMATRIX *pClassSize;
	struct INTMATRIX *pClassID;
	struct INTMATRIX *pPermClassID;
	struct INTMATRIX *pDataClassID;

	/* variance group */
	int nVargroupNum;
	struct INTMATRIX *pVargroupSize;
	struct INTMATRIX **vVargroupMap;

	/* permutation group */
	int nPermgroupNum = 0;
	struct INTMATRIX *pPermgroupSize;
	struct INTMATRIX **vPermgroupMap;

	/* files */
	FILE *fpData;
	FILE *fpCriteria;
	FILE *fpProbe;

	/* strings */
	char strFileName[MED_LINE_LENGTH];
	char strPermOutFile[MED_LINE_LENGTH];
	char strLine[MED_LINE_LENGTH];
	char strLongLine[LONG_LINE_LENGTH];
	char strComparisons[MED_LINE_LENGTH];
	char strCommand[MED_LINE_LENGTH];
	int nPairwiseComp;
	
	/* scores */
	struct DOUBLEMATRIX *pScore;
	struct DOUBLEMATRIX *pPermScore;
	
	/* parameters */
	double dTruncLow = 0.0;
	int nTakeLog = 0;
	double dLog2 = log(2.0);
	int nFDRPermNum = 0;
	int nCycPermNum = 0;

	/* error */
	int nError;

	/* count */
	int ni,nj,nk;
	char *chSep,*chSep2;
	double dValue;
	int nTemp;


	/* --------- */
	/* load data */
	/* --------- */
	printf("\n");
	printf("/* ---------------------- */\n");
	printf("/* TileMap Probe Summary  */\n");
	printf("/* ---------------------- */\n");
	printf("Loading data...\n");
	
	/* header info */
	nError = 0;
	fpCriteria = NULL;
	fpCriteria = fopen(strCriteriaFile, "rt");
	if(fpCriteria == NULL)
	{
		printf("Error: TileMap_ProbeSelection_Main, cannot open comparison info file!\n");
		exit(EXIT_FAILURE);
	}
	while(fgets(strLine, MED_LINE_LENGTH, fpCriteria) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		/* load basic info */

		/* array number */
		if(strstr(strLine, "[Array number]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			if(*chSep == '\0')
			{
				printf("Error: TileMap_ProbeSelection_Main, you have to specify the number of arrays!\n");
				exit(EXIT_FAILURE);
			}
			nArrayNum = atoi(chSep);
		}
		
		/* group number */
		else if(strstr(strLine, "[Group number]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			if(*chSep == '\0')
			{
				printf("Error: TileMap_ProbeSelection_Main, you have to specify the number of groups!\n");
				exit(EXIT_FAILURE);
			}
			nClassNum = atoi(chSep);
		}
			
		/* load class info */
		else if(strcmp(strLine, "[Group ID]") == 0)
		{
			fgets(strLongLine, LONG_LINE_LENGTH, fpCriteria);
			pDataClassID = Expression_GeneSelection_LoadGroupId(strLongLine);
			if(pDataClassID == NULL)
			{
				
				printf("Error: TileMap_ProbeSelection_Main, you have to specify group id for each array!\n");
				nError = 1;
				break;
			}
		}

		/* load comparisons */
		else if(strcmp(strLine, "[Comparisons]") == 0)
		{
			fgets(strComparisons, MED_LINE_LENGTH, fpCriteria);
			StrTrimLeft(strComparisons);
			StrTrimRight(strComparisons);
			if(strComparisons[0] == '\0')
			{
				printf("Error: TileMap_Probeselection_Main, you have to specify patterns you want to select!\n");
				nError = 1;
				break;
			}

			nPairwiseComp = 0;
			chSep = strpbrk(strComparisons, "<>" );
			while(chSep != NULL)
			{
				nPairwiseComp++;
				chSep = strpbrk((chSep+1), "<>" );
			}

			if(nPairwiseComp > 1)
				nPairwiseComp = 0;
		}

		/* load preprocessing */
		else if(strstr(strLine, "[Truncation lower bound]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			if(*chSep == '\0')
			{
				printf("Error: TileMap_ProbeSelection_Main, you have to specify the truncation lower bound!\n");
				exit(EXIT_FAILURE);
			}
			dTruncLow = atof(chSep);
		}

		else if(strstr(strLine, "[Take log2 before calculation?]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			if(*chSep == '\0')
			{
				printf("Error: TileMap_ProbeSelection_Main, you have to specify whether to take log2 transformation!\n");
				exit(EXIT_FAILURE);
			}
			nTakeLog = atoi(chSep);
		}


		else if(strstr(strLine, "[Monte Carlo draws for posterior prob.]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			if(*chSep == '\0')
			{
				if(nPairwiseComp == 0)
				{
					printf("Error: TileMap_ProbeSelection_Main, you have to specify the number of MC draws for computing posterior probabilities!\n");
					exit(EXIT_FAILURE);
				}
			}
			else
			{
				nCycPermNum = atoi(chSep);
			}
		}

		/* variance group */
		else if(strstr(strLine, "[Common variance groups]") == strLine)
		{
			/* group number */
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			if(*chSep == '\0')
			{
				printf("Error: TileMap_ProbeSelection_Main, you have to specify the common variance groups!\n");
				exit(EXIT_FAILURE);
			}
			else
			{
				nVargroupNum = atoi(chSep);
			}
			

			if(nVargroupNum <= 0)
			{
				printf("Error: TileMap_ProbeSelection_Main, there must be at least one common variance group!\n");
				exit(EXIT_FAILURE);
			}

			/* groups */
			pVargroupSize = NULL;
			pVargroupSize = CreateIntMatrix(nVargroupNum, 1);
			if(pVargroupSize == NULL)
			{
				nError = 1;
				break;
			}
			vVargroupMap = NULL;
			vVargroupMap = (struct INTMATRIX **)calloc(nVargroupNum, sizeof(struct INTMATRIX *));
			if(vVargroupMap == NULL)
			{
				DestroyIntMatrix(pVargroupSize);
				nError = 1;
				break;
			}

			for(ni=0; ni<nVargroupNum; ni++)
			{
				fgets(strLine, MED_LINE_LENGTH, fpCriteria);
				vVargroupMap[ni] = Expression_GeneSelection_LoadGroupId(strLine);
				if(vVargroupMap[ni] == NULL)
				{
					nError = 1;
					break;
				}
				IMSETAT(pVargroupSize, ni, 0, vVargroupMap[ni]->nWidth);
			}
		}


		/* number of permutations */
		else if(strstr(strLine, "[Number of permutations]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			if(*chSep == '\0')
			{
				if(nApplyPerm == 1)
				{
					printf("Error: TileMap_ProbeSelection_Main, you have to specify the number of permutations for estimating FDR!\n");
					exit(EXIT_FAILURE);
				}
			}
			else
			{
				nFDRPermNum = atoi(chSep);
			}

			*nPermNum = nFDRPermNum;
		}
			

		/* permutation group */
		else if(strstr(strLine, "[Exchangeable groups]") == strLine)
		{
			/* group number */
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			if(*chSep == '\0')
			{
				if(nApplyPerm == 1)
				{
					printf("Error: TileMap_ProbeSelection_Main, you have to specify the number of exchangeble groups in permutation test!\n");
					exit(EXIT_FAILURE);
				}
			}
			else
			{
				nPermgroupNum = atoi(chSep);
			}

			/* groups */
			if(nPermgroupNum > 0)
			{
				pPermgroupSize = NULL;
				pPermgroupSize = CreateIntMatrix(nPermgroupNum, 1);
				if(pPermgroupSize == NULL)
				{
					nError = 1;
					break;
				}
				vPermgroupMap = NULL;
				vPermgroupMap = (struct INTMATRIX **)calloc(nPermgroupNum, sizeof(struct INTMATRIX *));
				if(vPermgroupMap == NULL)
				{
					DestroyIntMatrix(pPermgroupSize);
					nError = 1;
					break;
				}

				for(ni=0; ni<nPermgroupNum; ni++)
				{
					fgets(strLine, MED_LINE_LENGTH, fpCriteria);
					vPermgroupMap[ni] = Expression_GeneSelection_LoadGroupId(strLine);
					if(vPermgroupMap[ni] == NULL)
					{
						nError = 1;
						break;
					}
					IMSETAT(pPermgroupSize, ni, 0, vPermgroupMap[ni]->nWidth);
				}
			}
		}

		
		/* do nothing */
		else if(strLine[0] == '#')
		{
		}

		/* else error */
		else
		{
			printf("Error: TileMap_Probeselection_Main, unknown parameters in *.cmpinfo file!\n");
			exit(EXIT_FAILURE);
		}
	}

	fclose(fpCriteria);

	if(nError > 0)
	{
		printf("Error: TileMap_ProbeSelection_Main, comparison info file format wrong!\n");
		exit(EXIT_FAILURE);
	}

	/* prepare space */
	pClassID = NULL;
	pClassID = CreateIntMatrix(1, nArrayNum);
	pClassSize = NULL;
	pClassSize = CreateIntMatrix(1, nClassNum);
	nj = 0;
	for(ni=0; ni<pDataClassID->nWidth; ni++)
	{
		nk = IMGETAT(pDataClassID, 0, ni);
		if( nk > nClassNum )
		{
			printf("Error: TileMap_ProbeSelection_Main, group id out of range!\n");
			exit(EXIT_FAILURE);
		}
		if( nk>0 )
		{
			if(nj >= nArrayNum)
			{
				printf("Error: TileMap_ProbeSelection_Main, group id out of range!\n");
				exit(EXIT_FAILURE);
			}
			IMSETAT(pClassID, 0, nj, nk);
			nTemp = IMGETAT(pClassSize, 0, (nk-1))+1;
			IMSETAT(pClassSize, 0, (nk-1), nTemp);
			nj++;
		}
	}
	if(nj != nArrayNum)
	{
		printf("Error: TileMap_ProbeSelection_Main, array number not match!\n");
		exit(EXIT_FAILURE);
	}

	vArray = NULL;
	vArray = (struct DOUBLEMATRIX **)calloc(nArrayNum, sizeof(struct DOUBLEMATRIX *));
	if(vArray == NULL)
	{
		printf("Error: TileMap_ProbeSelection_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nArrayNum; ni++)
	{
		vArray[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vArray[ni] == NULL)
		{
			printf("Error: TileMap_ProbeSelection_Main, cannot create data space!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* data */
	fpData = NULL;
	fpData = fopen(strDataFile, "r");
	if(fpData == NULL)
	{
		printf("Error: TileMap_ProbeSelection_Main, cannot open data file!\n");
		exit(EXIT_FAILURE);
	}
	sprintf(strFileName, "%s.pb", strOutFile);
	fpProbe = NULL;
	fpProbe = fopen(strFileName, "w");
	if(fpProbe == NULL)
	{
		printf("Error: TileMap_ProbeSelection_Main, cannot open probe file!\n");
		exit(EXIT_FAILURE);
	}
	fgets(strLongLine, LONG_LINE_LENGTH, fpData);
	nj = 0;
	while(fgets(strLongLine, LONG_LINE_LENGTH, fpData) != NULL)
	{
		StrTrimLeft(strLongLine);
		StrTrimRight(strLongLine);
		if(strLongLine[0] == '\0')
			continue;

		if(nj >= nProbeNum)
		{
			printf("Error: TileMap_ProbeSelection_Main, probeset number not match!\n");
			exit(EXIT_FAILURE);
		}

		ni = 0;
		nk = 0;
		chSep = strchr(strLongLine, '\t');
		if(chSep != NULL)
			*chSep = '\0';
		fprintf(fpProbe, "%s\t", strLongLine);

		chSep++;
		chSep2 = strchr(chSep, '\t');
		if(chSep2 != NULL)
			*chSep2 = '\0';
		fprintf(fpProbe, "%s\n", chSep);
		chSep = chSep2;

		while(chSep != NULL)
		{
			if( nk >= pDataClassID->nWidth )
			{
				printf("Error: TileMap_ProbeSelection_Main, array number not match!\n");
				exit(EXIT_FAILURE);
			}

			chSep++;
			chSep2 = strchr(chSep, '\t');
		
			if( pDataClassID->pMatElement[nk] > 0 )
			{
				/* middle number */
				if(chSep2 != NULL)
				{
					*chSep2 = '\0';

					if(chSep == chSep2)
					{
						dValue = 0.0;
					}
					else
					{
						dValue = atof(chSep);
					}
				}
				/* last number */
				else
				{
					if(chSep == chSep2)
					{
						dValue = 0.0;
					}
					else
					{
						dValue = atof(chSep);
					}
				}

				if(dValue < dTruncLow)
					dValue = dTruncLow;
				if(nTakeLog == 1)
					dValue = log(dValue)/dLog2;
			
				DMSETAT(vArray[ni], nj, 0, dValue);
				ni++;
			}

			/* get next */
			nk++;
			chSep = chSep2;
		}

		/* get next line */
		nj++;
	}

	fclose(fpData);
	fclose(fpProbe);
	DestroyIntMatrix(pDataClassID);

	if(nj != nProbeNum)
	{
		printf("Error: TileMap_ProbeSelection_Main, probeset number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* -------------------- */
	/* get original ranking */
	/* -------------------- */
	printf("Computing probe level score...\n");
	pScore = NULL;
	if(nPairwiseComp == 1)
	{
		*nProbeSummaryRange = 2;
		pScore = TileMap_ProbeSelection_tTest(nProbeNum, nArrayNum, vArray,
									nClassNum, pClassSize, pClassID,
									nVargroupNum, pVargroupSize, vVargroupMap, 
									strComparisons);
	}
	else
	{
		*nProbeSummaryRange = 1;
		pScore = TileMap_ProbeSelection_MonteCarlo(nProbeNum, nArrayNum, vArray,
									nClassNum, pClassSize, pClassID,
									nVargroupNum, pVargroupSize, vVargroupMap, 
									nCycPermNum, strComparisons);
	}


	/* -------------------------------- */
	/* save scores                      */
	/* -------------------------------- */
	printf("Exporting score...\n");
	sprintf(strFileName, "%s.pb", strOutFile);
	if(nApplyFilter == 0)
	{
		sprintf(strPermOutFile, "%s_f", strOutFile);
	}
	else
	{
		sprintf(strPermOutFile, "%s", strOutFile);
	}
	TileMap_ProbeSelection_Output(nProbeNum, pScore, strFileName, strPermOutFile);

	/* ---------------------------- */
	/* permutations to estimate FDR */
	/* ---------------------------- */
	if((nApplyPerm == 1) && (nFDRPermNum > 0))
	{
		printf("Permutations...\n");
		
		/* cycles */
		for(ni=0; ni<nFDRPermNum; ni++)
		{
			printf("perm %d...\n", ni);
			/* permute the class label */
			pPermClassID = NULL;
			pPermClassID = Expression_GeneSelection_ClassPerm(nClassNum, pClassID, pClassSize,
				nPermgroupNum, pPermgroupSize, vPermgroupMap);

			/* rerun a single cycle */
			pPermScore = NULL;
			if(nPairwiseComp == 1)
			{
				pPermScore = TileMap_ProbeSelection_tTest(nProbeNum, nArrayNum, vArray,
									nClassNum, pClassSize, pPermClassID,
									nVargroupNum, pVargroupSize, vVargroupMap, 
									strComparisons);
			}
			else
			{
				pPermScore = TileMap_ProbeSelection_MonteCarlo(nProbeNum, nArrayNum, vArray,
									nClassNum, pClassSize, pPermClassID,
									nVargroupNum, pVargroupSize, vVargroupMap, 
									nCycPermNum, strComparisons);
			}

			/* save scores */
			/* printf("exporting score...\n"); */
			sprintf(strFileName, "%s.pb", strOutFile);
			if(nApplyFilter == 0)
			{
				sprintf(strPermOutFile, "%s_perm%d_f", strOutFile, ni);
			}
			else
			{
				sprintf(strPermOutFile, "%s_perm%d", strOutFile, ni);
			}
			TileMap_ProbeSelection_Output(nProbeNum, pPermScore, strFileName, strPermOutFile);
			
			/* release memory */
			DestroyIntMatrix(pPermClassID);
			DestroyDoubleMatrix(pPermScore);
		}
	}

	/* -------------- */
	/* release memory */
	/* -------------- */

	sprintf(strFileName, "%s.pb", strOutFile);
	if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
	{
		sprintf(strCommand, "del %s", strFileName);
		system(strCommand);
	}
	else
	{
		sprintf(strCommand, "rm %s", strFileName);
		system(strCommand);
	}

	/* destroy arrays */
	for(ni=0; ni<nArrayNum; ni++)
	{
		DestroyDoubleMatrix(vArray[ni]);
		vArray[ni] = NULL;
	}
	free(vArray);

	/* destroy class ids */
	DestroyIntMatrix(pClassSize);
	DestroyIntMatrix(pClassID);

	/* destroy variance group ids */
	for(ni=0; ni<nVargroupNum; ni++)
	{
		DestroyIntMatrix(vVargroupMap[ni]);
		vVargroupMap[ni] = NULL;
	}
	free(vVargroupMap);
	DestroyIntMatrix(pVargroupSize);

	/* destroy permutation group ids */
	for(ni=0; ni<nPermgroupNum; ni++)
	{
		DestroyIntMatrix(vPermgroupMap[ni]);
		vPermgroupMap[ni] = NULL;
	}
	free(vPermgroupMap);
	DestroyIntMatrix(pPermgroupSize);
	

	/* destroy results */
	DestroyDoubleMatrix(pScore);

	/* ------ */
	/* return */
	/* ------ */
	return PROC_SUCCESS;
}


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
					char strComparisons[])
{
	/* ------ */
	/* define */
	/* ------ */
	struct INTMATRIX *pClassSizeCopy;
	struct INTMATRIX *pDf;
	struct DOUBLEMATRIX *pScore;
	struct DOUBLEMATRIX **vMeans;
	struct DOUBLEMATRIX **vVars;
	struct DOUBLEMATRIX **vSDs;
	
	struct DOUBLEMATRIX *pSigCoef;
	double *vExp,*vSum,*vAve,*vSigma,*vMu,*vSigma2,*vMu2;
	double dExp,dExp2,dVarTemp;
	int nClustId;
	double dDenom,dTemp;
	double dVarmean,dVarsst;
	double dB,dV2,dK,dN;

	char vLogic[MED_LINE_LENGTH];
	char vTNumber[LINE_LENGTH];
	double vGid[MED_LINE_LENGTH];
	
	int nLogicLen;
	int nLeftNum,nRightNum,nTNumLen;

	int ni,nj,nk,nLen;

	int ng1,ng2;

	/* ------------- */
	/* expression    */
	/* ------------- */
	nLogicLen = 0;
	StrTrimLeft(strComparisons);
	StrTrimRight(strComparisons);
	if(strComparisons[0] == '\0')
	{
		printf("Error: Expression_GeneSelection_t-Test, no logic expression formula available!\n");
		exit(EXIT_FAILURE);
	}
	nLen = strlen(strComparisons);
	nj = 0;
	nLeftNum = 0;
	nRightNum = 0;
	nTNumLen = 0;
	for(ni=0; ni<nLen; ni++)
	{
		if( strComparisons[ni] == '(' )
		{
			if(nTNumLen != 0)
			{
				printf("Error: Expression_GeneSelection_t-Test, logic expression wrong!\n");
				exit(EXIT_FAILURE);
			}

			nLeftNum++;
		}
		else if(strComparisons[ni] == ')')
		{
			if(nTNumLen > 0)
			{
				vTNumber[nTNumLen] = '\0';
				vLogic[nj] = 'G';
				vGid[nj] = atof(vTNumber);
				nTNumLen = 0;
				nj++;
			}
			nRightNum++;
		}
		else if( (strComparisons[ni] == '<') || (strComparisons[ni] == '>') 
			|| (strComparisons[ni] == '&') || (strComparisons[ni] == '|') )
		{
			if(nTNumLen > 0)
			{
				vTNumber[nTNumLen] = '\0';
				vLogic[nj] = 'G';
				vGid[nj] = atof(vTNumber);
				nTNumLen = 0;
				nj++;
			}
			vLogic[nj] = strComparisons[ni];
			vGid[nj] = -1.0;
			nj++;	
		}
		else if( (strComparisons[ni] >= '0') && (strComparisons[ni] <= '9') )
		{
			vTNumber[nTNumLen] = strComparisons[ni];
			nTNumLen++;
		}
		else if( (strComparisons[ni] == ' ') || (strComparisons[ni] == '\t') )
		{
			/* do nothing */
		}
		else
		{
			printf("Error: Expression_GeneSelection_t-Test, %c not supported in logic expressions!\n", strComparisons[ni]);
			exit(EXIT_FAILURE);
		}
	}
	if(nTNumLen > 0)
	{
		vTNumber[nTNumLen] = '\0';
		vLogic[nj] = 'G';
		vGid[nj] = atof(vTNumber);
		nTNumLen = 0;
		nj++;
	}
	nLogicLen = nj;

	if(nLeftNum != nRightNum)
	{
		printf("Error: Expression_GeneSelection_t-Test, lack '(' or ')'!\n");
		exit(EXIT_FAILURE);
	}

	if(nLogicLen != 3)
	{
		printf("Error: Expression_GeneSelection_t-Test, logic expression wrong!\n");
		exit(EXIT_FAILURE);
	}
	if((vLogic[0] != 'G') || (vLogic[2] != 'G'))
	{
		printf("Error: Expression_GeneSelection_t-Test, logic expression wrong!\n");
		exit(EXIT_FAILURE);
	}
	if((vLogic[1] != '<') && (vLogic[1] != '>'))
	{
		printf("Error: Expression_GeneSelection_t-Test, logic expression wrong!\n");
		exit(EXIT_FAILURE);
	}

	if(vLogic[1] == '<')
	{
		ng1 = (int)vGid[0]-1;
		ng2 = (int)vGid[2]-1;
	}
	else if(vLogic[1] == '>')
	{
		ng1 = (int)vGid[2]-1;
		ng2 = (int)vGid[0]-1;
	}
	else
	{
		printf("Error: Expression_GeneSelection_t-Test, logic expression wrong!\n");
		exit(EXIT_FAILURE);
	}

	if((ng1<0) || (ng1>=nClassNum) || (ng2<0) || (ng2>=nClassNum))
	{
		printf("Error: Expression_GeneSelection_t-Test, group id out of range!\n");
		exit(EXIT_FAILURE);
	}

	/* for debug purpose*/
	/* for(ni=0; ni<nLogicLen; ni++)
	{
		printf("%c %d\n", vLogic[ni], (int)(vGid[ni]));
	} */
	
	/* ---- */
	/* init */
	/* ---- */
	pScore = NULL;
	pScore = CreateDoubleMatrix(1, nProbeNum);
	pClassSizeCopy = NULL;
	pClassSizeCopy = CreateIntMatrix(1, nClassNum);
	pDf = NULL;
	pDf = CreateIntMatrix(1, nVargroupNum);


	/* ------------- */
	/* create memory */
	/* ------------- */
	vMeans = NULL;
	vMeans = (struct DOUBLEMATRIX **)calloc(nClassNum, sizeof(struct DOUBLEMATRIX *));
	if(vMeans == NULL)
	{
		printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
		exit(EXIT_FAILURE);
	}
	vVars = NULL;
	vVars = (struct DOUBLEMATRIX **)calloc(nClassNum, sizeof(struct DOUBLEMATRIX *));
	if(vVars == NULL)
	{
		printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
		exit(EXIT_FAILURE);
	}
	vSDs = NULL;
	vSDs = (struct DOUBLEMATRIX **)calloc(nVargroupNum, sizeof(struct DOUBLEMATRIX *));
	if(vSDs == NULL)
	{
		printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
		exit(EXIT_FAILURE);
	}
		
	for(ni=0; ni<nClassNum; ni++)
	{
		vMeans[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vMeans[ni] == NULL)
		{
			printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
			exit(EXIT_FAILURE);
		}

		vVars[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vVars[ni] == NULL)
		{
			printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	for(ni=0; ni<nVargroupNum; ni++)
	{
		vSDs[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vSDs[ni] == NULL)
		{
			printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* ---------------- */
	/* get mean and var */
	/* ---------------- */

	/* mean */
	for(ni=0; ni<nArrayNum; ni++)
	{
		nClustId = (pClassID->pMatElement)[ni]-1;
		(pClassSizeCopy->pMatElement)[nClustId] += 1;
		vSum = vMeans[nClustId]->pMatElement;
		vExp = vArray[ni]->pMatElement;
		for(nj=0; nj<nProbeNum; nj++)
		{
			vSum[nj] += vExp[nj];
		}
	}
	
	for(ni=0; ni<nClassNum; ni++)
	{
		if((pClassSizeCopy->pMatElement)[ni] != (pClassSize->pMatElement)[ni])
		{
			printf("Error: Expression_GeneSelection_t-Test, class size not match!\n");
			exit(EXIT_FAILURE);
		}
		if((pClassSizeCopy->pMatElement)[ni] > 0)
		{
			dDenom = (double)((pClassSizeCopy->pMatElement)[ni]);
			vSum = vMeans[ni]->pMatElement;
			for(nj=0; nj<nProbeNum; nj++)
			{
				vSum[nj] /= dDenom;
			}
		}
	}

	/* variance: sum of squares */
	for(ni=0; ni<nArrayNum; ni++)
	{
		nClustId = (pClassID->pMatElement)[ni]-1;
		vSum = vVars[nClustId]->pMatElement;
		vAve = vMeans[nClustId]->pMatElement;
		vExp = vArray[ni]->pMatElement;
		for(nj=0; nj<nProbeNum; nj++)
		{
			dTemp = (vExp[nj]-vAve[nj]);		
			vSum[nj] += dTemp*dTemp;
		}
	}

	/* variance: estimates */
	for(ni=0; ni<nVargroupNum; ni++)
	{
		if(vVargroupMap[ni]->nWidth != pVargroupSize->pMatElement[ni])
		{
			printf("Error: Expression_GeneSelection_t-Test, variance group number not match!\n");
			exit(EXIT_FAILURE);
		}

		for(nj=0; nj<vVargroupMap[ni]->nWidth; nj++)
		{
			nClustId = (vVargroupMap[ni]->pMatElement)[nj]-1;
			if((pClassSize->pMatElement)[nClustId] > 0)
				(pDf->pMatElement)[ni] += (pClassSize->pMatElement)[nClustId]-1;
			vSum = vSDs[ni]->pMatElement;
			vExp = vVars[nClustId]->pMatElement;
			for(nk=0; nk<nProbeNum; nk++)
			{
				vSum[nk] += vExp[nk];
			}
		}
	}
	

	/* ------------- */
	/* shrinking var */
	/* ------------- */
	for(ni=0; ni<nVargroupNum; ni++)
	{
		dVarmean = 0.0;
		dVarsst = 0.0;
		dB = 0.0;

		if((pDf->pMatElement)[ni] > 0)
		{
			dDenom = (double)(pDf->pMatElement)[ni];
			vSum = vSDs[ni]->pMatElement;

			/* xbar */
			for(nj=0; nj<nProbeNum; nj++)
			{
				vSum[nj] /= dDenom;
				dVarmean += vSum[nj];
			}
			dVarmean /= (double)nProbeNum;

			/* sst */
			for(nj=0; nj<nProbeNum; nj++)
			{
				dTemp = vSum[nj]-dVarmean;
				dVarsst += dTemp*dTemp;
			}

			/* shrinkage factor */
			dV2 = 2.0/dDenom;
			dK = (double)nProbeNum;
			dN = 1.0;
			if((dK > 0.0) && (dVarsst > 0.0))
				dB = dV2*(dK-1)/((dN+dV2)*dK)+dN*(dK-1)*dV2*dVarmean*dVarmean/((dN+dV2)*dN*dVarsst);
			else
				dB = 0.0;

			if(dB < 0.0)
				dB = 0.0;
			else if(dB > 1.0)
				dB = 1.0;

			/* shrink variance */
			for(nj=0; nj<nProbeNum; nj++)
			{
				vSum[nj] = (1-dB)*vSum[nj]+dB*dVarmean+1e-16;
				/* vSum[nj] += 1e-16; */
			}
		}
	}

	for(ni=0; ni<nClassNum; ni++)
	{
		DestroyDoubleMatrix(vVars[ni]);
		vVars[ni] = NULL;
	}

	for(ni=0; ni<nVargroupNum; ni++)
	{
		for(nj=0; nj<vVargroupMap[ni]->nWidth; nj++)
		{
			nClustId = (vVargroupMap[ni]->pMatElement)[nj]-1;
			vVars[nClustId] = vSDs[ni];
		}
	}

	/* ----------- */
	/* simulation  */
	/* ----------- */
	pSigCoef = NULL;
	pSigCoef = CreateDoubleMatrix(1, nClassNum);
	vExp = pSigCoef->pMatElement;
	for(ni=0; ni<nClassNum; ni++)
	{
		if((pClassSize->pMatElement)[ni] > 0)
		{
			dDenom = (double)((pClassSize->pMatElement)[ni]);
			vExp[ni] = 1.0/dDenom;
		}
		else
		{
			vExp[ni] = 0.0;
		}
	}


	vMu = vMeans[ng1]->pMatElement;
	vSigma = vVars[ng1]->pMatElement;
	vMu2 = vMeans[ng2]->pMatElement;
	vSigma2 = vVars[ng2]->pMatElement;
	dExp = pSigCoef->pMatElement[ng1];
	dExp2 = pSigCoef->pMatElement[ng2];
	vSum = pScore->pMatElement;

	/* probeset by probeset */
	for(nj=0; nj<nProbeNum; nj++)
	{
		if((vSigma[nj] <= 0.0) || (vSigma2[nj] <= 0.0))
		{
			printf("Warning: Expression_GeneSelection_t-Test, variance=0, may not have enough sample to estimate variance!\n");
		}
		dVarTemp = sqrt(vSigma[nj]*dExp+vSigma2[nj]*dExp2);
		if(dVarTemp > 0.0)
			vSum[nj] = (vMu[nj]-vMu2[nj])/dVarTemp;
		else
			vSum[nj] = 0.0;
	}
	
	DestroyDoubleMatrix(pSigCoef);

	/* -------------- */
	/* release memory */
	/* -------------- */
	for(ni=0; ni<nClassNum; ni++)
	{
		DestroyDoubleMatrix(vMeans[ni]);
		vMeans[ni] = NULL;
	}

	for(ni=0; ni<nVargroupNum; ni++)
	{
		DestroyDoubleMatrix(vSDs[ni]);
		vSDs[ni] = NULL;		
	}

	free(vMeans);
	free(vVars);
	free(vSDs);

	DestroyIntMatrix(pClassSizeCopy);
	DestroyIntMatrix(pDf);

	/* ------ */
	/* return */
	/* ------ */
	return pScore;
}


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
					int nCycPermNum, char strComparisons[])
{
	/* ------ */
	/* define */
	/* ------ */
	struct INTMATRIX *pClassSizeCopy;
	struct INTMATRIX *pDf;
	struct DOUBLEMATRIX *pScore;
	struct DOUBLEMATRIX **vMeans;
	struct DOUBLEMATRIX **vVars;
	struct DOUBLEMATRIX **vSDs;
	/* struct DOUBLEMATRIX *pMeanDraws; */
	struct DOUBLEMATRIX **vMeanDraws;
	struct DOUBLEMATRIX *pSigCoef;
	double *vExp,*vSum,*vAve,*vSigma,*vMu;
	unsigned char *vEval;
	int nClustId;
	double dDenom,dTemp;
	double dVarmean,dVarsst;
	double dB,dV2,dK,dN;

	char vLogic[MED_LINE_LENGTH];
	char vTNumber[LINE_LENGTH];
	double vGid[MED_LINE_LENGTH];
	struct DOUBLEMATRIX *vVid[MED_LINE_LENGTH];
	struct BYTEMATRIX *vLid[MED_LINE_LENGTH];
	int nLogicLen;
	int nLeftNum,nRightNum,nTNumLen;
	int nTrue;
	struct BYTEMATRIX *pTrueVec;
	struct DOUBLEMATRIX *pTrueVal;

	int ni,nj,nk,nLen;

	/* ------------- */
	/* expression    */
	/* ------------- */
	nLogicLen = 0;
	StrTrimLeft(strComparisons);
	StrTrimRight(strComparisons);
	if(strComparisons[0] == '\0')
	{
		printf("Error: Expression_GeneSelection_MonteCarlo, no logic expression formula available!\n");
		exit(EXIT_FAILURE);
	}
	nLen = strlen(strComparisons);
	nj = 0;
	nLeftNum = 0;
	nRightNum = 0;
	nTNumLen = 0;
	for(ni=0; ni<nLen; ni++)
	{
		if( strComparisons[ni] == '(' )
		{
			if(nTNumLen != 0)
			{
				printf("Error: Expression_GeneSelection_MonteCarlo, logic expression wrong!\n");
				exit(EXIT_FAILURE);
			}

			vLogic[nj] = strComparisons[ni];
			vGid[nj] = -1.0;
			nj++;
			nLeftNum++;
		}
		else if(strComparisons[ni] == ')')
		{
			if(nTNumLen > 0)
			{
				vTNumber[nTNumLen] = '\0';
				vLogic[nj] = 'G';
				vGid[nj] = atof(vTNumber);
				nTNumLen = 0;
				nj++;
			}
			vLogic[nj] = strComparisons[ni];
			vGid[nj] = -1.0;
			nj++;
			nRightNum++;
		}
		else if( (strComparisons[ni] == '<') || (strComparisons[ni] == '>') 
			|| (strComparisons[ni] == '&') || (strComparisons[ni] == '|') )
		{
			if(nTNumLen > 0)
			{
				vTNumber[nTNumLen] = '\0';
				vLogic[nj] = 'G';
				vGid[nj] = atof(vTNumber);
				nTNumLen = 0;
				nj++;
			}
			vLogic[nj] = strComparisons[ni];
			vGid[nj] = -1.0;
			nj++;
			
		}
		else if( (strComparisons[ni] >= '0') && (strComparisons[ni] <= '9') )
		{
			vTNumber[nTNumLen] = strComparisons[ni];
			nTNumLen++;
		}
		else if( (strComparisons[ni] == ' ') || (strComparisons[ni] == '\t') )
		{
			/* do nothing */
		}
		else
		{
			printf("Error: Expression_GeneSelection_MonteCarlo, %c not supported in logic expressions!\n", strComparisons[ni]);
			exit(EXIT_FAILURE);
		}
	}
	if(nTNumLen > 0)
	{
		vTNumber[nTNumLen] = '\0';
		vLogic[nj] = 'G';
		vGid[nj] = atof(vTNumber);
		nTNumLen = 0;
		nj++;
	}
	nLogicLen = nj;

	if(nLeftNum != nRightNum)
	{
		printf("Error: Expression_GeneSelection_MonteCarlo, lack '(' or ')'!\n");
		exit(EXIT_FAILURE);
	}

	/* for debug purpose*/
	/* for(ni=0; ni<nLogicLen; ni++)
	{
		printf("%c %d\n", vLogic[ni], (int)(vGid[ni]));
	} */
	
	/* ---- */
	/* init */
	/* ---- */
	pScore = NULL;
	pScore = CreateDoubleMatrix(1, nProbeNum);
	pClassSizeCopy = NULL;
	pClassSizeCopy = CreateIntMatrix(1, nClassNum);
	pDf = NULL;
	pDf = CreateIntMatrix(1, nVargroupNum);


	/* ------------- */
	/* create memory */
	/* ------------- */
	vMeans = NULL;
	vMeans = (struct DOUBLEMATRIX **)calloc(nClassNum, sizeof(struct DOUBLEMATRIX *));
	if(vMeans == NULL)
	{
		printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
		exit(EXIT_FAILURE);
	}
	vVars = NULL;
	vVars = (struct DOUBLEMATRIX **)calloc(nClassNum, sizeof(struct DOUBLEMATRIX *));
	if(vVars == NULL)
	{
		printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
		exit(EXIT_FAILURE);
	}
	vSDs = NULL;
	vSDs = (struct DOUBLEMATRIX **)calloc(nVargroupNum, sizeof(struct DOUBLEMATRIX *));
	if(vSDs == NULL)
	{
		printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
		exit(EXIT_FAILURE);
	}
	vMeanDraws = NULL;
	vMeanDraws = (struct DOUBLEMATRIX **)calloc(nClassNum, sizeof(struct DOUBLEMATRIX *));
	if(vMeanDraws == NULL)
	{
		printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
		exit(EXIT_FAILURE);
	}
	
	for(ni=0; ni<nClassNum; ni++)
	{
		vMeans[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vMeans[ni] == NULL)
		{
			printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
			exit(EXIT_FAILURE);
		}

		vVars[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vVars[ni] == NULL)
		{
			printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
			exit(EXIT_FAILURE);
		}

		vMeanDraws[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vMeanDraws[ni] == NULL)
		{
			printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation.!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	for(ni=0; ni<nVargroupNum; ni++)
	{
		vSDs[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vSDs[ni] == NULL)
		{
			printf("Error: Expression_GeneSelection_MonteCarlo, cannot create space for simulation!\n");
			exit(EXIT_FAILURE);
		}
	}

	/*pMeanDraws = NULL;
	pMeanDraws = CreateDoubleMatrix(1,nClassNum);
	*/

	/* ---------------- */
	/* get mean and var */
	/* ---------------- */

	/* mean */
	for(ni=0; ni<nArrayNum; ni++)
	{
		nClustId = (pClassID->pMatElement)[ni]-1;
		(pClassSizeCopy->pMatElement)[nClustId] += 1;
		vSum = vMeans[nClustId]->pMatElement;
		vExp = vArray[ni]->pMatElement;
		for(nj=0; nj<nProbeNum; nj++)
		{
			vSum[nj] += vExp[nj];
		}
	}
	
	for(ni=0; ni<nClassNum; ni++)
	{
		if((pClassSizeCopy->pMatElement)[ni] != (pClassSize->pMatElement)[ni])
		{
			printf("Error: Expression_GeneSelection_MonteCarlo, class size not match!\n");
			exit(EXIT_FAILURE);
		}
		if((pClassSizeCopy->pMatElement)[ni] > 0)
		{
			dDenom = (double)((pClassSizeCopy->pMatElement)[ni]);
			vSum = vMeans[ni]->pMatElement;
			for(nj=0; nj<nProbeNum; nj++)
			{
				vSum[nj] /= dDenom;
			}
		}
	}

	/* variance: sum of squares */
	for(ni=0; ni<nArrayNum; ni++)
	{
		nClustId = (pClassID->pMatElement)[ni]-1;
		vSum = vVars[nClustId]->pMatElement;
		vAve = vMeans[nClustId]->pMatElement;
		vExp = vArray[ni]->pMatElement;
		for(nj=0; nj<nProbeNum; nj++)
		{
			dTemp = (vExp[nj]-vAve[nj]);		
			vSum[nj] += dTemp*dTemp;
		}
	}

	/* variance: estimates */
	for(ni=0; ni<nVargroupNum; ni++)
	{
		if(vVargroupMap[ni]->nWidth != pVargroupSize->pMatElement[ni])
		{
			printf("Error: Expression_GeneSelection_MonteCarlo, variance group number not match!\n");
			exit(EXIT_FAILURE);
		}

		for(nj=0; nj<vVargroupMap[ni]->nWidth; nj++)
		{
			nClustId = (vVargroupMap[ni]->pMatElement)[nj]-1;
			if((pClassSize->pMatElement)[nClustId] > 0)
				(pDf->pMatElement)[ni] += (pClassSize->pMatElement)[nClustId]-1;
			vSum = vSDs[ni]->pMatElement;
			vExp = vVars[nClustId]->pMatElement;
			for(nk=0; nk<nProbeNum; nk++)
			{
				vSum[nk] += vExp[nk];
			}
		}
	}
	

	/* ------------- */
	/* shrinking var */
	/* ------------- */
	for(ni=0; ni<nVargroupNum; ni++)
	{
		dVarmean = 0.0;
		dVarsst = 0.0;
		dB = 0.0;

		if((pDf->pMatElement)[ni] > 0)
		{
			dDenom = (double)(pDf->pMatElement)[ni];
			vSum = vSDs[ni]->pMatElement;

			/* xbar */
			for(nj=0; nj<nProbeNum; nj++)
			{
				vSum[nj] /= dDenom;
				dVarmean += vSum[nj];
			}
			dVarmean /= (double)nProbeNum;

			/* sst */
			for(nj=0; nj<nProbeNum; nj++)
			{
				dTemp = vSum[nj]-dVarmean;
				dVarsst += dTemp*dTemp;
			}

			/* shrinkage factor */
			dV2 = 2.0/dDenom;
			dK = (double)nProbeNum;
			dN = 1.0;
			if((dK > 0.0) && (dVarsst > 0.0))
				dB = dV2*(dK-1)/((dN+dV2)*dK)+dN*(dK-1)*dV2*dVarmean*dVarmean/((dN+dV2)*dN*dVarsst);
			else
				dB = 0.0;

			if(dB < 0.0)
				dB = 0.0;
			else if(dB > 1.0)
				dB = 1.0;

			/* shrink variance */
			for(nj=0; nj<nProbeNum; nj++)
			{
				vSum[nj] = (1-dB)*vSum[nj]+dB*dVarmean;
				/* vSum[nj] = sqrt(vSum[nj]); */
 				vSum[nj] = sqrt(vSum[nj]) + 1e-16;
			}
		}
	}

	for(ni=0; ni<nClassNum; ni++)
	{
		DestroyDoubleMatrix(vVars[ni]);
		vVars[ni] = NULL;
	}

	for(ni=0; ni<nVargroupNum; ni++)
	{
		for(nj=0; nj<vVargroupMap[ni]->nWidth; nj++)
		{
			nClustId = (vVargroupMap[ni]->pMatElement)[nj]-1;
			vVars[nClustId] = vSDs[ni];
		}
	}

	/* ----------- */
	/* simulation  */
	/* ----------- */
	pSigCoef = NULL;
	pSigCoef = CreateDoubleMatrix(1, nClassNum);
	vExp = pSigCoef->pMatElement;
	for(ni=0; ni<nClassNum; ni++)
	{
		if((pClassSize->pMatElement)[ni] > 0)
		{
			dDenom = (double)((pClassSize->pMatElement)[ni]);
			vExp[ni] = sqrt(1.0/dDenom);
		}
		else
		{
			vExp[ni] = 0.0;
		}
	}

	for(ni=0; ni<nLogicLen; ni++)
	{
		if(vGid[ni] > 0.0)
		{
			nk = (int)(vGid[ni]-1.0);
			vVid[ni] = vMeanDraws[nk];
		}
		else
		{
			vVid[ni] = NULL;
		}
		vLid[ni] = NULL;
	}

	for(ni=0; ni<nCycPermNum; ni++)
	{
		if(ni%100 == 0)
		{
			printf("iter %d...\n", ni);
		}

		vExp = pSigCoef->pMatElement;
		
		/* probeset by probeset */
		for(nk=0; nk<nClassNum; nk++)
		{
			vAve = vMeanDraws[nk]->pMatElement;
			vMu = vMeans[nk]->pMatElement;
			vSigma = vVars[nk]->pMatElement;
			dTemp = normrnd(0.0, 1.0);
			for(nj=0; nj<nProbeNum; nj++)
			{
				/* simulate */
				vAve[nj] = vMu[nj]+vExp[nk]*vSigma[nj]*dTemp;
			}
		}

		

		/* evaluate */
		pTrueVec = NULL;
		pTrueVal = NULL;
		nTrue = Expression_GeneSelection_EvaluateVec(nProbeNum, vLogic, vVid, vLid, nLogicLen, 0, &pTrueVec, &pTrueVal);
		if(nTrue != 4)
		{
			printf("Error: logic evaluation wrong!\n");
			exit(EXIT_FAILURE);
		}

		/* nTrue = Expression_GeneSelection_Evaluate(pMeanDraws, vLogic, vGid, nLogicLen, 0); */

		/* add */
		vSum = pScore->pMatElement;
		/* vSum[nj] += (double)(nTrue); */
		
		vEval = pTrueVec->pMatElement;
		for(nj=0; nj<nProbeNum; nj++)
		{
			vSum[nj] += (double)(vEval[nj]);
		}
		DestroyByteMatrix(pTrueVec);
	}
	DestroyDoubleMatrix(pSigCoef);

	/* normalize */
	dDenom = (double)nCycPermNum;
	vSum = pScore->pMatElement;
	for(nj=0; nj<nProbeNum; nj++)
	{
		vSum[nj] = 1.0-vSum[nj]/dDenom;
	}

	/* -------------- */
	/* release memory */
	/* -------------- */
	for(ni=0; ni<nClassNum; ni++)
	{
		DestroyDoubleMatrix(vMeans[ni]);
		vMeans[ni] = NULL;
		
		DestroyDoubleMatrix(vMeanDraws[ni]);
		vMeanDraws[ni] = NULL;
	}

	for(ni=0; ni<nVargroupNum; ni++)
	{
		DestroyDoubleMatrix(vSDs[ni]);
		vSDs[ni] = NULL;		
	}

	free(vMeans);
	free(vVars);
	free(vSDs);
	free(vMeanDraws);

	/* DestroyDoubleMatrix(pMeanDraws); */
	DestroyIntMatrix(pClassSizeCopy);
	DestroyIntMatrix(pDf);

	/* ------ */
	/* return */
	/* ------ */
	return pScore;
}


/* ----------------------------------------------------------------------- */ 
/*  TileMap_ProbeSelection_Output()                                        */
/*  Output results, a group of random control will be selected.            */
/* ----------------------------------------------------------------------- */ 
int TileMap_ProbeSelection_Output(int nProbeNum, struct DOUBLEMATRIX *pScore, 
				char strProbeFile[], char strOutFile[])
{
	/* define */
	struct tagString **vProbeName;
	int ni;
	char strOutPath[LINE_LENGTH];
	char strLine[LINE_LENGTH];
	FILE *fpIn;
	FILE *fpOut;
	double *vSS;

	/* init */
	vProbeName = NULL;
	vProbeName = (struct tagString **)calloc(nProbeNum, sizeof(struct tagString *));

	/* load probe name */
	fpIn = NULL;
	fpIn = fopen(strProbeFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileMap_ProbeSelection_Output, cannot open probeset name file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(ni >= nProbeNum)
		{
			printf("Error: TileMap_ProbeSelection_Output, probeset number not match!\n");
			exit(EXIT_FAILURE);
		}

		StringAddTail((vProbeName+ni), strLine);
		ni++;
	}

	fclose(fpIn);
	if(ni != nProbeNum)
	{
		printf("Error: TileMap_ProbeSelection_Output, probeset number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* output unsorted */
	sprintf(strOutPath, "%s_pb.sum", strOutFile);
	fpOut = NULL;
	fpOut = fopen(strOutPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: Tiling_ProbeSelection_Output, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	vSS = pScore->pMatElement;
	
	
	fprintf(fpOut, "chromosome\tposition\tscore\n");
	for(ni=0; ni<nProbeNum; ni++)
	{
		if(vProbeName[ni] != NULL)
		{
			fprintf(fpOut, "%s\t", vProbeName[ni]->m_pString);
		}
		else
		{
			fprintf(fpOut, "NA\t-1\t");
		}
		
		fprintf(fpOut, "% 9.7e\n", vSS[ni]);
	}

	fclose(fpOut);

	/* release memory */
	for(ni=0; ni<nProbeNum; ni++)
	{
		DeleteString(vProbeName[ni]);
		vProbeName[ni] = NULL;
	}
	free(vProbeName);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMap_BPMAPFilter_Pbsum()                                            */
/*  Filter out repeat and bad probes with bpmap informaiton.               */
/* ----------------------------------------------------------------------- */ 
int TileMap_BPMAPFilter_Pbsum(char strInFile[], char strRefFile[], char strOutFile[])
{
	/* define */
	FILE *fpIn;
	FILE *fpRef;
	FILE *fpOut;

	/* load line */
	char strLine[LINE_LENGTH];
	int nRefPos,nPos;
	char strChr[LINE_LENGTH];
	char strRefChr[LINE_LENGTH];
	double dPos;
	int nDepthNum,nRepeatNum;
	char strProbeSeq[LINE_LENGTH];
	double dValue;
	double dMean,dN;
	int nInFileEnd,nRefFileEnd;
	int nWritePos;
	int nDiffChr;
	int nRefCatch;
	int nRealProbeNum = 0;
	
	/* open file */
	fpIn = NULL;
	fpIn = fopen(strInFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileMap_BPMAPFilter_Pbsum, cannot read input file!\n");
		exit(EXIT_FAILURE);
	}
	fpRef = NULL;
	fpRef = fopen(strRefFile, "r");
	if(fpRef == NULL)
	{
		printf("Error: TileMap_BPMAPFilter_Pbsum, cannot read reference file!\n");
		exit(EXIT_FAILURE);
	}
	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: TileMap_BPMAPFilter_Pbsum, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	fprintf(fpOut, "chromosome\tposition\tscore\n");

	/* read files */
	nInFileEnd = 1;
	nRefFileEnd = 1;
	strcpy(strChr, "");
	strcpy(strRefChr, "");
	nPos = 0;
	nRefPos = 0;
	dValue = 0.0;
	nDepthNum = 0;
	nRepeatNum = 0;
	
	if(fgets(strLine, LINE_LENGTH, fpIn) == NULL)
	{
		printf("Error: TileMap_BPMAPFilter_Pbsum, empty *_pb.sum file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if((strLine[0] != '#') && (strLine[0] != '\0'))
		{
			nInFileEnd = 0;
			break;
		}
	}
	
	if(nInFileEnd == 0)
	{
		sscanf(strLine, "%s %lf %lf", strChr, &dPos, &dValue);
		nPos = (int)dPos;
	}

	if(fgets(strLine, LINE_LENGTH, fpRef) == NULL)
	{
		printf("Error: TileMap_BPMAPFilter_Pbsum, empty *.refmask file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LINE_LENGTH, fpRef) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] != '\0')
		{
			nRefFileEnd = 0;
			break;
		}
	}
	
	if(nRefFileEnd == 0)
	{
		sscanf(strLine, "%s %d %d %d %s", strRefChr, &nRefPos, &nDepthNum, &nRepeatNum, strProbeSeq);
	}
	
	dMean = 0.0;
	dN = 0.0;
	nRefCatch = 0;
	while((nInFileEnd == 0) && (nRefFileEnd == 0))
	{
		nDiffChr = strcmp(strChr, strRefChr);
		/* if equal position */
		if((nDiffChr == 0) && (nPos == nRefPos))
		{
			dMean += dValue;
			dN += 1.0;

			nInFileEnd = 1;
			nRefCatch = 0;
			while(fgets(strLine, LINE_LENGTH, fpIn) != NULL)
			{
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if((strLine[0] != '#') && (strLine[0] != '\0'))
				{
					nInFileEnd = 0;
					break;
				}
			}
			
			if(nInFileEnd == 0)
			{
				sscanf(strLine, "%s %lf %lf", strChr, &dPos, &dValue);
				nPos = (int)dPos;
			}
		}
		/* if pos > refpos */
		else if( ((nDiffChr == 0) && (nPos > nRefPos)) || (nDiffChr != 0))
		{
			if(nRefCatch == 1)
			{
				printf("Error: TileMap_BPMAPFilter_Pbsum, the *.refmask file does not contain all the unique probes!\n");
				exit(EXIT_FAILURE);
			}

			if(dN > 0.0)
			{
				if((int)dN != nDepthNum)
				{
					printf("Warning: TileMap_BPMAPFilter_Pbsum, depth number not consistent!\n");
				}

				dMean /= dN;
				nWritePos = nRefPos;
				if(nRepeatNum <= 1)
				{
					fprintf(fpOut, "%s\t%d\t% 9.7e\n", strRefChr, nRefPos, dMean);
					nRealProbeNum++;
				}
			}

			dMean = 0.0;
			dN = 0.0;

			nRefFileEnd = 1;
			while(fgets(strLine, LINE_LENGTH, fpRef) != NULL)
			{
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] != '\0')
				{
					nRefFileEnd = 0;
					break;
				}
			}
			
			if(nRefFileEnd == 0)
			{
				sscanf(strLine, "%s %d %d %d %s", strRefChr, &nRefPos, &nDepthNum, &nRepeatNum, strProbeSeq);
			}
			nRefCatch = 1;
		}
		else
		{
			printf("Error: TileMap_BPMAPFilter_Pbsum, the *.refmask file does not contain all the unique probes!\n");
			exit(EXIT_FAILURE);
		}
	}

	if((nInFileEnd == 1) && (nRefFileEnd ==0))
	{
		if(dN > 0.0)
		{
			if((int)dN != nDepthNum)
			{
				printf("Warning: TileMap_BPMAPFilter_Pbsum, depth number not consistent!\n");
			}
			dMean /= dN;
			if(nRepeatNum <= 1)
			{
				fprintf(fpOut, "%s\t%d\t% 9.7e\n", strRefChr, nRefPos, dMean);
				nRealProbeNum++;
			}
		}
	}
	else
	{
		printf("Error: TileMap_BPMAPFilter_Pbsum, the *.refmask file does not contain all the unique probes!\n");
		exit(EXIT_FAILURE);
	}

	/* close file */
	fclose(fpIn);
	fclose(fpRef);
	fclose(fpOut);


	/* return */
	return nRealProbeNum;
}


/* ----------------------------------------------------------------------- */ 
/*  TileMap_UMS_HMM_Main()                                                 */
/*  UMS for estimating HMM parameters.                                     */
/* ----------------------------------------------------------------------- */ 
int TileMap_UMS_HMM_Main(char strSelectPath[], char strScorePath[], 
						int nProbeNum,  int nScoreRange,
						double dPcut, double dQcut,
						int nStepSize, int nIntervalNum,
						int nFragLen, char strOutPath[])
{
	/* definition */
	FILE *fpSelect;
	FILE *fpScore;
	FILE *fpOut;
	struct DOUBLEMATRIX *pPos;
	struct DOUBLEMATRIX *pSelect;
	struct DOUBLEMATRIX *pScore;
	struct DOUBLEMATRIX *pIntx;
	struct DOUBLEMATRIX *pF0;
	struct DOUBLEMATRIX *pF1;
	struct DOUBLEMATRIX *pTrans;
	double dTheta;
	double a0,a1;

	double dPos,dScore;
	char strLine[LINE_LENGTH];
	char strChr[LINE_LENGTH];
	char strFileName[MED_LINE_LENGTH];
	int ni;

	/* init */
	pPos = NULL;
	pPos = CreateDoubleMatrix(1, nProbeNum);
	if(pPos == NULL)
	{
		printf("Error: TileMap_UMS_HMM_Main, cannot create position matrix!\n");
		exit(EXIT_FAILURE);
	}

	pSelect = NULL;
	pSelect = CreateDoubleMatrix(1, nProbeNum);
	if(pSelect == NULL)
	{
		printf("Error: TileMap_UMS_HMM_Main, cannot create selection matrix!\n");
		exit(EXIT_FAILURE);
	}

	pScore = NULL;
	pScore = CreateDoubleMatrix(1, nProbeNum);
	if(pScore == NULL)
	{
		printf("Error: TileMap_UMS_HMM_Main, cannot create score matrix!\n");
		exit(EXIT_FAILURE);
	}

	/* load data */
	fpSelect = NULL;
	fpSelect = fopen(strSelectPath, "r");
	if(fpSelect == NULL)
	{
		printf("Error: TileMap_UMS_HMM_Main, cannot open selection file!\n");
		exit(EXIT_FAILURE);
	}

	if(fgets(strLine, LINE_LENGTH, fpSelect) == NULL)
	{
		printf("Error: TileMap_UMS_HMM_Main, selection file is empty!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LINE_LENGTH, fpSelect) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%s %lf %lf", strChr, &dPos, &dScore);
		
		if(ni>=nProbeNum)
		{
			printf("Error: TileMap_UMS_HMM_Main, probe number not match!\n");
			exit(EXIT_FAILURE);
		}
		
		pPos->pMatElement[ni] = dPos;
		pSelect->pMatElement[ni] = dScore;
		ni++;
	}

	fclose(fpSelect);
	
	if(ni != nProbeNum)
	{
		printf("Error: TileMap_UMS_HMM_Main, probe number not match!\n");
		exit(EXIT_FAILURE);
	}

	fpScore = NULL;
	fpScore = fopen(strScorePath, "r");
	if(fpScore == NULL)
	{
		printf("Error: TileMap_UMS_HMM_Main, cannot open selection file!\n");
		exit(EXIT_FAILURE);
	}

	if(fgets(strLine, LINE_LENGTH, fpScore) == NULL)
	{
		printf("Error: TileMap_UMS_HMM_Main, score file is empty!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LINE_LENGTH, fpScore) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%s %lf %lf", strChr, &dPos, &dScore);
		
		if(ni>=nProbeNum)
		{
			printf("Error: TileMap_UMS_HMM_Main, probe number not match!\n");
			exit(EXIT_FAILURE);
		}
		
		if((int)dPos != (int)(pPos->pMatElement[ni]))
		{
			printf("Error: TileMap_UMS_HMM_Main, position not match!\n");
			exit(EXIT_FAILURE);
		}
		
		if(nScoreRange == 2)
		{
			dScore = exp(dScore);
			dScore = dScore/(1.0+dScore);
		}

		pScore->pMatElement[ni] = dScore;
		ni++;
	}

	fclose(fpScore);
	
	if(ni != nProbeNum)
	{
		printf("Error: TileMap_UMS_HMM_Main, probe number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* UMS */
	pIntx = NULL;
	pIntx = CreateDoubleMatrix(1, nIntervalNum);
	if(pIntx == NULL)
	{
		printf("Error: TileMap_UMS_HMM_Main, cannot create interval matrix!\n");
		exit(EXIT_FAILURE);
	}

	pF0 = NULL;
	pF0 = CreateDoubleMatrix(1, nIntervalNum);
	if(pF0 == NULL)
	{
		printf("Error: TileMap_UMS_HMM_Main, cannot create f0 matrix!\n");
		exit(EXIT_FAILURE);
	}

	pF1 = NULL;
	pF1 = CreateDoubleMatrix(1, nIntervalNum);
	if(pF1 == NULL)
	{
		printf("Error: TileMap_UMS_HMM_Main, cannot create f1 matrix!\n");
		exit(EXIT_FAILURE);
	}

	TileMap_UMS(pSelect, dPcut, dQcut, nStepSize, pScore, nIntervalNum, pIntx, pF0, pF1, &dTheta);

	/* estimate emission */
	pTrans = NULL;
	pTrans = CreateDoubleMatrix(2, 2);
	if(pTrans == NULL)
	{
		printf("Error: TileMap_UMS_HMM_Main, cannot create transition matrix!\n");
		exit(EXIT_FAILURE);
	}

	a1 = 1.0/(double)nFragLen;
	if((dTheta < ZERO_BOUND) || (dTheta > 1.0-ZERO_BOUND))
	{
		printf("Error: TileMap_UMS_HMM_Main, theta=0, cannot create HMM parameters!\n");
		exit(EXIT_FAILURE);
	}
	a0 = a1*dTheta/(1.0-dTheta);

	DMSETAT(pTrans, 0, 0, (1.0-a0));
	DMSETAT(pTrans, 0, 1, a0);
	DMSETAT(pTrans, 1, 0, a1);
	DMSETAT(pTrans, 1, 1, (1.0-a1));

	/* write */
	sprintf(strFileName, "%s_transp.txt", strOutPath);
	DMSAVE(pTrans, strFileName);

	sprintf(strFileName, "%s_emissp.txt", strOutPath);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: TileMap_UMS_HMM_Main, cannot export emission probability!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nIntervalNum; ni++)
	{
		fprintf(fpOut, "% 9.7e\t", pIntx->pMatElement[ni]);
	}
	fprintf(fpOut, "\n");
	for(ni=0; ni<nIntervalNum; ni++)
	{
		fprintf(fpOut, "% 9.7e\t", pF0->pMatElement[ni]);
	}
	fprintf(fpOut, "\n");
	for(ni=0; ni<nIntervalNum; ni++)
	{
		fprintf(fpOut, "% 9.7e\t", pF1->pMatElement[ni]);
	}
	fprintf(fpOut, "\n");

	fclose(fpOut);


	/* destroy memory */
	DestroyDoubleMatrix(pPos);
	DestroyDoubleMatrix(pSelect);
	DestroyDoubleMatrix(pScore);
	DestroyDoubleMatrix(pIntx);
	DestroyDoubleMatrix(pF0);
	DestroyDoubleMatrix(pF1);
	DestroyDoubleMatrix(pTrans);

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  TileMap_UMS()                                                          */
/*  Unbalanced Mixture Subtraction.                                        */
/* ----------------------------------------------------------------------- */ 
int TileMap_UMS(struct DOUBLEMATRIX *pSelect, double dPcut, double dQcut,
				   int nStepSize, struct DOUBLEMATRIX *pScore, 
				   int nIntervalNum, struct DOUBLEMATRIX *pIntx,
				   struct DOUBLEMATRIX *pFhat0, struct DOUBLEMATRIX *pFhat1,
				   double *dTheta)
{
	/* define */
	struct DOUBLEMATRIX *pH;
	struct DOUBLEMATRIX *pG0;
	struct DOUBLEMATRIX *pG1;
	struct DOUBLEMATRIX *pF1;
	struct DOUBLEMATRIX *pLikeRatio;
	
	int nPrcNum = 100;
	struct DOUBLEMATRIX *pG0Prc;
	struct DOUBLEMATRIX *pG1Prc;
	struct DOUBLEMATRIX *pR;
	double r;
	struct DOUBLEMATRIX *pScorePrctile;
	double dZeroCut = 1e-6;

	struct DOUBLEMATRIX *pSelectSort;
	struct DOUBLEMATRIX *pScoreSort;

	int ni,nj,nk;
	int nSelectWid;
	double dTotNum,dG0Num,dG1Num;
	double dTp,dTq;
	double dScore;
	double dSum;
	double theta;

	/* check parameters */
	if((pSelect == NULL) || (pScore == NULL))
	{
		printf("Error: Tiling_UMS, score/selection matrices not specified!\n");
		exit(EXIT_FAILURE);
	}
	if((pIntx == NULL) || (pFhat0 == NULL) || (pFhat1 == NULL))
	{
		printf("Error: Tiling_UMS, output matrices not specified!\n");
		exit(EXIT_FAILURE);
	}
	if( (pIntx->nWidth != nIntervalNum) || (pFhat0->nWidth != nIntervalNum) || (pFhat1->nWidth != nIntervalNum) )
	{
		printf("Error: Tiling_UMS, matrices dimensions and number of intervals specified do not match!\n");
		exit(EXIT_FAILURE);
	}
	if(nStepSize < 0)
	{
		printf("Error: Tiling_UMS, please set a nonnegative interval number as stepsize of selection criteria!\n");
		exit(EXIT_FAILURE);
	}
	if((pScore->nWidth != pSelect->nWidth) || (pScore->nWidth <= nStepSize))
	{
		printf("Error: Tiling_UMS, matrix size not match!\n");
		exit(EXIT_FAILURE);
	}
	nSelectWid = pSelect->nWidth - nStepSize;
	if(nIntervalNum <= 0)
	{
		printf("Error: Tiling_UMS, please set a positive interval number for dividing [0,1]!\n");
		exit(EXIT_FAILURE);
	}
	if( (dPcut < 0.0) || (dPcut >= 1.0) || (dQcut <= 0.0) || (dQcut > 1.0) )
	{
		printf("Error: Tiling_UMS, p% and q% in selection criteria should fall within [0,1]!\n");
		exit(EXIT_FAILURE);
	}


	/* prepare space */
	pH = NULL;
	pH = CreateDoubleMatrix(1, nIntervalNum);
	if(pH == NULL)
	{
		printf("Error: Tiling_UMS, cannot create memory for UMS!\n");
		exit(EXIT_FAILURE);
	}

	pG0 = NULL;
	pG0 = CreateDoubleMatrix(1, nIntervalNum);
	if(pG0 == NULL)
	{
		printf("Error: Tiling_UMS, cannot create memory for UMS!\n");
		exit(EXIT_FAILURE);
	}

	pG1 = NULL;
	pG1 = CreateDoubleMatrix(1, nIntervalNum);
	if(pG1 == NULL)
	{
		printf("Error: Tiling_UMS, cannot create memory for UMS!\n");
		exit(EXIT_FAILURE);
	}

	pF1 = NULL;
	pF1 = CreateDoubleMatrix(1, nIntervalNum);
	if(pF1 == NULL)
	{
		printf("Error: Tiling_UMS, cannot create memory for UMS!\n");
		exit(EXIT_FAILURE);
	}

	pLikeRatio = NULL;
	pLikeRatio = CreateDoubleMatrix(1, nIntervalNum);
	if(pLikeRatio == NULL)
	{
		printf("Error: Tiling_UMS, cannot create memory for UMS!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nIntervalNum; ni++)
	{
		pH->pMatElement[ni] = 1.0;
		pG0->pMatElement[ni] = 1.0;
		pG1->pMatElement[ni] = 1.0;
	}

	pScorePrctile = NULL;
	pScorePrctile = CreateDoubleMatrix(1, (nPrcNum-1));
	if(pScorePrctile == NULL)
	{
		printf("Error: Tiling_UMS, cannot create memory for UMS!\n");
		exit(EXIT_FAILURE);
	}

	pG0Prc = NULL;
	pG0Prc = CreateDoubleMatrix(1, (nPrcNum-1));
	if(pG0Prc == NULL)
	{
		printf("Error: Tiling_UMS, cannot create memory for UMS!\n");
		exit(EXIT_FAILURE);
	}

	pG1Prc = NULL;
	pG1Prc = CreateDoubleMatrix(1, (nPrcNum-1));
	if(pG1Prc == NULL)
	{
		printf("Error: Tiling_UMS, cannot create memory for UMS!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<(nPrcNum-1); ni++)
	{
		pG0Prc->pMatElement[ni] = 1.0;
		pG1Prc->pMatElement[ni] = 1.0;
	}
	
	/* get percentiles */
	pSelectSort = NULL;
	DMSORTMERGEA_0(pSelect, &pSelectSort, NULL);
	nj = (int)(pSelect->nWidth*dPcut);
	if(nj == pSelect->nWidth)
		nj--;
	dTp = pSelectSort->pMatElement[nj];

	nj = (int)(pSelect->nWidth*dQcut);
	if(nj == pSelect->nWidth)
		nj--;
	dTq = pSelectSort->pMatElement[nj];
	DestroyDoubleMatrix(pSelectSort);

	pScoreSort = NULL;
	DMSORTMERGEA_0(pScore, &pScoreSort, NULL);
	for(ni=1; ni<nPrcNum; ni++)
	{
		nj = (int)((pScore->nWidth)*(double)ni/(double)nPrcNum);
		if(nj == pScore->nWidth)
			nj--;
		pScorePrctile->pMatElement[ni-1] = pScoreSort->pMatElement[nj];
	}

	/* get h(t) */
	for(nj=0; nj<pScore->nWidth; nj++)
	{
		ni = (int)((pScore->pMatElement[nj]) * nIntervalNum);
		if(ni == nIntervalNum)
			ni--;

		pH->pMatElement[ni] = pH->pMatElement[ni]+1.0; 
	}
	dTotNum = pScore->nWidth+nIntervalNum;
	DMPDIVTS(pH, dTotNum);

	/* DMSAVE(pH, "testdistn.txt"); */

	/* get g0(t), g1(t) */
	dG0Num = nIntervalNum;
	dG1Num = nIntervalNum;
	for(nj=0; nj<nSelectWid; nj++)
	{
		/* if g0(t) */
		if(pSelect->pMatElement[nj] > dTp)
		{
			dScore = pScore->pMatElement[nj+nStepSize];
			ni = (int)(dScore * nIntervalNum);
			if(ni == nIntervalNum)
				ni--;

			pG0->pMatElement[ni] = pG0->pMatElement[ni]+1.0;
			dG0Num += 1.0;

			for(nk=(nPrcNum-2); nk>=0; nk--)
			{
				if(dScore >= pScorePrctile->pMatElement[nk])
				{
					pG0Prc->pMatElement[nk] += 1.0; 
				}
			}
		}

		/* if g1(t) */
		if(pSelect->pMatElement[nj] <= dTq)
		{
			dScore = pScore->pMatElement[nj+nStepSize];
			ni = (int)(dScore * nIntervalNum);
			if(ni == nIntervalNum)
				ni--;

			pG1->pMatElement[ni] = pG1->pMatElement[ni]+1.0;
			dG1Num += 1.0;

			for(nk=(nPrcNum-2); nk>=0; nk--)
			{
				if(dScore >= pScorePrctile->pMatElement[nk])
				{
					pG1Prc->pMatElement[nk] += 1.0; 
				}
			}
		}
	}
	DMPDIVTS(pG0, dG0Num);
	DMPDIVTS(pG1, dG1Num);
	DMPDIVTS(pG0Prc, (dG0Num-nIntervalNum+1.0));
	DMPDIVTS(pG1Prc, (dG1Num-nIntervalNum+1.0));
	
	/*DMSAVE(pG0, "testdistng0.txt");
	DMSAVE(pG1, "testdistng1.txt");
	DMSAVE(pG0Prc, "testdistnr0.txt");
	DMSAVE(pG1Prc, "testdistnr1.txt");
	*/

	/* get lim g1(t)/g0(t) */
	for(ni=0; ni<(nPrcNum-1); ni++)
	{
		pG1Prc->pMatElement[ni] /= pG0Prc->pMatElement[ni];
	}
	
	/* DMSAVE(pG1Prc, "testdistnr.txt"); */

	pR = NULL;
	DMSORTMERGEA_0(pG1Prc, &pR, NULL);
	ni = (int)(0.5*(nPrcNum-1));
	r = pR->pMatElement[ni];
	if(r >= 1.0)
	{
		printf("Error: Tiling_UMS, r>=1.0, the estimate might be wrong!\n");
		r = 0.999;
	}
	DestroyDoubleMatrix(pR);


	/* estimate f1(t) */
	for(ni=0; ni<nIntervalNum; ni++)
	{
		pF1->pMatElement[ni] = (pG1->pMatElement[ni]-r*pG0->pMatElement[ni])/(1.0-r);
		if(pF1->pMatElement[ni] < dZeroCut)
			pF1->pMatElement[ni] = dZeroCut;
		if(pF1->pMatElement[ni] > 1.0-dZeroCut)
			pF1->pMatElement[ni] = 1.0-dZeroCut;

		pLikeRatio->pMatElement[ni] = pF1->pMatElement[ni]/pG0->pMatElement[ni];
	}
	
	for(ni=1; ni<nIntervalNum; ni++)
	{
		if(pLikeRatio->pMatElement[ni] <= 0.0)
		{
			pLikeRatio->pMatElement[ni] = pLikeRatio->pMatElement[ni-1];
		}
		else if(pLikeRatio->pMatElement[ni] > pLikeRatio->pMatElement[ni-1])
		{
			pLikeRatio->pMatElement[ni] = pLikeRatio->pMatElement[ni-1];
		}
	}

	/* DMSAVE(pLikeRatio, "testlikeratio.txt"); */

	dSum = 0.0;
	for(ni=0; ni<nIntervalNum; ni++)
	{
		pF1->pMatElement[ni] = pG0->pMatElement[ni]*pLikeRatio->pMatElement[ni];
		dSum += pF1->pMatElement[ni];
	}
	DMPDIVTS(pF1, dSum);
	/* DMSAVE(pF1, "testdistf1.txt"); */

	/* estimate theta */
	dScore = 0.0;
	dSum = 0.0;
	for(ni=0; ni<nIntervalNum; ni++)
	{
		dScore += (pH->pMatElement[ni]-pG0->pMatElement[ni])*(pF1->pMatElement[ni]-pG0->pMatElement[ni]);
		dSum += (pF1->pMatElement[ni]-pG0->pMatElement[ni])*(pF1->pMatElement[ni]-pG0->pMatElement[ni]);
	}
	theta = dScore/dSum;
	if(theta > 1.0)
		theta = 1.0;
	if(theta < 0.0)
		theta = 0.0;

	/* return parameters */
	*dTheta = theta;
	for(ni=0; ni<nIntervalNum; ni++)
	{
		pIntx->pMatElement[ni] = (double)(ni+1)/(double)nIntervalNum;
		pFhat0->pMatElement[ni] = pG0->pMatElement[ni];
		pFhat1->pMatElement[ni] = pF1->pMatElement[ni];
	}


	/* release memory */
	DestroyDoubleMatrix(pH);
	DestroyDoubleMatrix(pG0);
	DestroyDoubleMatrix(pG1);
	DestroyDoubleMatrix(pF1);
	DestroyDoubleMatrix(pLikeRatio);
	DestroyDoubleMatrix(pG0Prc);
	DestroyDoubleMatrix(pG1Prc);
	DestroyDoubleMatrix(pScorePrctile);
	DestroyDoubleMatrix(pScoreSort);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMap_BindingRegionSelection_HMM_Main()                              */
/*  HMM Forward-Backward algorithm for calculating the posterior prob. of  */
/*  binding for each position.                                             */
/* ----------------------------------------------------------------------- */ 
int TileMap_BindingRegionSelection_HMM_Main(int nProbeNum, char strDataPath[],
				char strTransformType[], double dScoreResolution,
				char strTransitionPath[], char strEmissionPath[], 
				double dGapDist, double dPosteriorCutoff,
				char strOutPath[])
{
	/* define */
	struct DOUBLEMATRIX *pScore;
	struct DOUBLEMATRIX *pNewScore;
	struct INTMATRIX *pChr;
	int nChr;
	struct tagString **vChrName;
	struct DOUBLEMATRIX *pPosition;
	struct DOUBLEMATRIX *pTransition;
	struct DOUBLEMATRIX *pEmission;
	struct DOUBLEMATRIX *pStationary;
	int nStateNum;
	struct DOUBLEMATRIX **vPosterior;
	double dTemp;
	int ni,nj;
	double *pEle;

	/* for loading data */
	char strLine[LINE_LENGTH];
	char strChr[LINE_LENGTH];
	double dPos,dScore;
	FILE *fpIn;
	FILE *fpOut;

	
	/* check */

	/* init */
	pTransition = NULL;
	pTransition = DMLOAD(strTransitionPath);
	if(pTransition == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_HMM_Main, cannot load transition probability!\n");
		exit(EXIT_FAILURE);
	}

	pEmission = NULL;
	pEmission = DMLOAD(strEmissionPath);
	if(pEmission == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_HMM_Main, cannot load emission probability!\n");
		exit(EXIT_FAILURE);
	}

	if((pTransition->nHeight+1) != pEmission->nHeight)
	{
		printf("Error: TileMap_BindingRegionSelection_HMM_Main, state number not match!\n");
		exit(EXIT_FAILURE);
	}

	nStateNum = pTransition->nHeight;
	if(nStateNum != 2)
	{
		printf("Error: TileMap_BindingRegionSelection_HMM_Main, this function only support 2 states now!\n");
		exit(EXIT_FAILURE);
	}

	pStationary = NULL;
	pStationary = CreateDoubleMatrix(1, nStateNum);
	dTemp = DMGETAT(pTransition, 1, 0)+DMGETAT(pTransition, 0, 1);
	pStationary->pMatElement[0] = log(DMGETAT(pTransition, 1, 0)/dTemp);
	pStationary->pMatElement[1] = log(DMGETAT(pTransition, 0, 1)/dTemp);

	pEle = pTransition->pMatElement;
	for(ni=0; ni<pTransition->nHeight; ni++)
	{
		for(nj=0; nj<pTransition->nWidth; nj++)
		{
			*pEle = log(*pEle);
			pEle++;
		}
	}

	pEle = pEmission->pMatElement+pEmission->nWidth;
	for(ni=1; ni<pEmission->nHeight; ni++)
	{
		for(nj=0; nj<pEmission->nWidth; nj++)
		{
			*pEle = log(*pEle);
			pEle++;
		}
	}

	

	/* prepare space */
	pScore = NULL;
	pScore = CreateDoubleMatrix(1, nProbeNum);
	if(pScore == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_HMM_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}

	pChr = NULL;
	pChr = CreateIntMatrix(1, nProbeNum);
	if(pChr == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_HMM_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}

	vChrName = NULL;
	vChrName = (struct tagString **)calloc(nProbeNum, sizeof(struct tagString *));
	if(vChrName == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_HMM_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}

	pPosition = NULL;
	pPosition = CreateDoubleMatrix(1, nProbeNum);
	if(pPosition == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_HMM_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}

	vPosterior = NULL;
	vPosterior = (struct DOUBLEMATRIX **)calloc(nStateNum, sizeof(struct DOUBLEMATRIX *));
	if(vPosterior == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_HMM_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nStateNum; ni++)
	{
		vPosterior[ni] = CreateDoubleMatrix(1, nProbeNum);
		if(vPosterior[ni] == NULL)
		{
			printf("Error: TileMap_BindingRegionSelection_HMM_Main, cannot create data space!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* load data */
	fpIn = NULL;
	fpIn = fopen(strDataPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_HMM_Main, cannot open data file!\n");
		exit(EXIT_FAILURE);
	}

	if(fgets(strLine, LINE_LENGTH, fpIn) == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_HMM_Main, empty data file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	nChr = 0;
	while(fgets(strLine, LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%s %lf %lf", strChr, &dPos, &dScore);
		/* sscanf(strLine, "%lf %lf", &dPos, &dScore); */
		pScore->pMatElement[ni] = dScore;
		pPosition->pMatElement[ni] = dPos;
		StringAddTail((vChrName+ni), strChr);

		if(ni == 0)
		{
			pChr->pMatElement[ni] = 0;
		}
		else
		{
			if(strcmp(vChrName[ni-1]->m_pString, vChrName[ni]->m_pString) == 0)
			{
				pChr->pMatElement[ni] = nChr;
			}
			else
			{
				nChr++;
				pChr->pMatElement[ni] = nChr;
			}
		}
		ni++;
	}

	fclose(fpIn);

	if(ni != nProbeNum)
	{
		printf("Error: TileMap_BindingRegionSelection_HMM_Main, probe number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* transform score */
	pNewScore = NULL;
	pNewScore = Tiling_BindingRegionSelection_ScoreTransform(pScore, dScoreResolution, strTransformType);
	DestroyDoubleMatrix(pScore);

	/* calculate posterior probability */
	TileMap_BindingRegionSelection_HMM(nProbeNum, nStateNum,
				pStationary, pTransition, pEmission, dGapDist,
				pNewScore, pChr, pPosition, vPosterior);

	/* save */
	sprintf(strLine, "%s_hmm.sum", strOutPath);
	fpOut = NULL;
	fpOut = fopen(strLine, "w");
	if(fpOut == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_HMM_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "chromosome\tposition\tposterior_prob\n");
	for(ni=0; ni<nProbeNum; ni++)
	{
		fprintf(fpOut, "%s\t%d", vChrName[ni]->m_pString, (int)(pPosition->pMatElement[ni]));
		for(nj=1; nj<nStateNum; nj++)
		{
			fprintf(fpOut, "\t%9.7e", vPosterior[nj]->pMatElement[ni]);
		}
		fprintf(fpOut, "\n");
	}

	fclose(fpOut);

	/* save to *.bed file */
	TileMap_CallBindingRegion_HMM(nProbeNum, vPosterior[0],
				pChr, pPosition, vChrName, (1.0-dPosteriorCutoff), dGapDist,
				strOutPath);

	/* release memory */
	DestroyDoubleMatrix(pTransition);
	DestroyDoubleMatrix(pEmission);
	DestroyDoubleMatrix(pStationary);
	DestroyDoubleMatrix(pNewScore);
	DestroyIntMatrix(pChr);
	DestroyDoubleMatrix(pPosition);
	for(ni=0; ni<nStateNum; ni++)
	{
		DestroyDoubleMatrix(vPosterior[ni]);
		vPosterior[ni] = NULL;
	}
	free(vPosterior);
	for(ni=0; ni<nProbeNum; ni++)
	{
		DeleteString(vChrName[ni]);
		vChrName[ni] = NULL;
	}
	free(vChrName);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMap_BindingRegionSelection_HMM()                                   */
/*  HMM Forward-Backward algorithm for calculating the posterior prob. of  */
/*  binding for each position.                                             */
/* ----------------------------------------------------------------------- */ 
int TileMap_BindingRegionSelection_HMM(int nProbeNum, int nStateNum,
				struct DOUBLEMATRIX *pStationary, struct DOUBLEMATRIX *pTransition, 
				struct DOUBLEMATRIX *pEmission, double dGapDist, 
				struct DOUBLEMATRIX *pScore, struct INTMATRIX *pChr,
				struct DOUBLEMATRIX *pPosition, struct DOUBLEMATRIX **vPosterior)
{
	/* define */
	struct DOUBLEMATRIX **vForwardSum;
	struct DOUBLEMATRIX **vBackwardSum;
	struct DOUBLEMATRIX *pDP;
	double dSum,dSumF,dSumB;
	double dTemp,dMax;
	double dDist;
	int ni,nj,nk;
	double dInitMax;

	int nEqualLenInterval;
	double dIntS,dIntE,dIntStep;

	/* check */
	if( (nProbeNum <= 0) || (nStateNum <= 0) )
	{
		printf("Error: TileMap_BindingRegionSelection_HMM, nProbeNum/nStateNum <=0!\n");
		exit(EXIT_FAILURE);
	}
	if((pScore == NULL) || (pStationary == NULL) || (pTransition == NULL) || (pChr == NULL)
		|| (pEmission == NULL) || (pPosition == NULL) || (vPosterior == NULL))
	{
		printf("Error: TileMap_BindingRegionSelection_HMM, no input data/parameters!\n");
		exit(EXIT_FAILURE);
	}
	if(pScore->nWidth != nProbeNum)
	{
		printf("Error: TileMap_BindingRegionSelection_HMM, probe number not match!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nStateNum; ni++)
	{
		if(vPosterior[ni]->nWidth != nProbeNum)
		{
			printf("Error: TileMap_BindingRegionSelection_HMM, probe number not match!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* init */
	dInitMax = -DM_ACCESS_VIOLATION;
	
	/* by default, use intervals of equal length for likelihood calculation */
	nEqualLenInterval = 1;
	if(pEmission->nWidth <= 2)
	{
		printf("Error: TileMap_BindingRegionSelection_HMM, emission probability need to be set for more than two intervals!\n");
		exit(EXIT_FAILURE);
	}

	dIntS = DMGETAT(pEmission, 0, 0);
	dIntE = DMGETAT(pEmission, 0, (pEmission->nWidth-2));
	
	dIntStep = (dIntE-dIntS)/(double)(pEmission->nWidth-2);
	if(dIntStep < 0.0)
	{
		printf("Error: TileMap_BindingRegionSelection_HMM, negative interval step length!\n");
		exit(EXIT_FAILURE);
	}


	/* prepare space */
	vForwardSum = NULL;
	vForwardSum = (struct DOUBLEMATRIX **)calloc(nStateNum, sizeof(struct DOUBLEMATRIX *));
	if(vForwardSum == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_HMM, cannot create space for HMM calculation!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nStateNum; ni++)
	{
		vForwardSum[ni] = CreateDoubleMatrix(1, nProbeNum);
		if(vForwardSum[ni] == NULL)
		{
			printf("Error: TileMap_BindingRegionSelection_HMM, cannot create space for HMM calculation!\n");
			exit(EXIT_FAILURE);
		}
	}

	vBackwardSum = NULL;
	vBackwardSum = (struct DOUBLEMATRIX **)calloc(nStateNum, sizeof(struct DOUBLEMATRIX *));
	if(vBackwardSum == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_HMM, cannot create space for HMM calculation!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nStateNum; ni++)
	{
		vBackwardSum[ni] = CreateDoubleMatrix(1, nProbeNum);
		if(vBackwardSum[ni] == NULL)
		{
			printf("Error: TileMap_BindingRegionSelection_HMM, cannot create space for HMM calculation!\n");
			exit(EXIT_FAILURE);
		}
	}

	pDP = NULL;
	pDP = CreateDoubleMatrix(nStateNum, 1);
	if(pDP == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_HMM, cannot create space for HMM calculation!\n");
		exit(EXIT_FAILURE);
	}

	/* forward summation */
	for(ni=0; ni<nStateNum; ni++)
	{
		dTemp = pStationary->pMatElement[ni]+Tiling_BindingRegionSelection_HMM_GetEmission(pEmission, ni, pScore->pMatElement[0], nEqualLenInterval, dIntS, dIntE, dIntStep);
		vForwardSum[ni]->pMatElement[0] = dTemp;
	}

	for(nj=1; nj<nProbeNum; nj++)
	{
		dDist = pPosition->pMatElement[nj] - pPosition->pMatElement[nj-1];
		if(pChr->pMatElement[nj] != pChr->pMatElement[nj-1])
			dDist = dGapDist+1e6;

		for(ni=0; ni<nStateNum; ni++)
		{
			dMax = dInitMax;
			for(nk=0; nk<nStateNum; nk++)
			{
				dTemp = vForwardSum[nk]->pMatElement[nj-1] + Tiling_BindingRegionSelection_HMM_GetTransition(pStationary, pTransition, nk, ni, dDist, dGapDist); 
				if(dTemp > dMax)
				{
					dMax = dTemp;
				}
				pDP->pMatElement[nk] = dTemp;
			}

			dSum = 0.0;
			for(nk=0; nk<nStateNum; nk++)
			{
				dTemp = pDP->pMatElement[nk]-dMax;
				dSum += exp(dTemp);
			}
			dSum = log(dSum)+dMax+Tiling_BindingRegionSelection_HMM_GetEmission(pEmission, ni, pScore->pMatElement[nj], nEqualLenInterval, dIntS, dIntE, dIntStep);
			vForwardSum[ni]->pMatElement[nj] = dSum;
		}
	}

	
	dMax = dInitMax;
	for(ni=0; ni<nStateNum; ni++)
	{
		dTemp = vForwardSum[ni]->pMatElement[nj-1]; 
		if(dTemp > dMax)
		{
			dMax = dTemp;
		}
		pDP->pMatElement[ni] = dTemp;
	}

	dSum = 0.0;
	for(ni=0; ni<nStateNum; ni++)
	{
		dTemp = pDP->pMatElement[ni]-dMax;
		dSum += exp(dTemp);
	}
	dSumF = log(dSum)+dMax;


	/* backward summation */
	for(ni=0; ni<nStateNum; ni++)
	{
		vBackwardSum[ni]->pMatElement[nProbeNum-1] = 0.0;
	}
	for(nj=nProbeNum-2; nj>=0; nj--)
	{
		dDist = pPosition->pMatElement[nj+1]-pPosition->pMatElement[nj];
		if(pChr->pMatElement[nj+1] != pChr->pMatElement[nj])
			dDist = dGapDist+1e6;

		for(ni=0; ni<nStateNum; ni++)
		{
			dMax = dInitMax;
			for(nk=0; nk<nStateNum; nk++)
			{
				dTemp = vBackwardSum[nk]->pMatElement[nj+1] 
					+ Tiling_BindingRegionSelection_HMM_GetTransition(pStationary, pTransition, ni, nk, dDist, dGapDist) 
					+ Tiling_BindingRegionSelection_HMM_GetEmission(pEmission, nk, pScore->pMatElement[nj+1], nEqualLenInterval, dIntS, dIntE, dIntStep); 
				if(dTemp > dMax)
				{
					dMax = dTemp;
				}
				pDP->pMatElement[nk] = dTemp;
			}

			dSum = 0.0;
			for(nk=0; nk<nStateNum; nk++)
			{
				dTemp = pDP->pMatElement[nk]-dMax;
				dSum += exp(dTemp);
			}
			dSum = log(dSum)+dMax;
			vBackwardSum[ni]->pMatElement[nj] = dSum;
		}
	}

	dMax = dInitMax;
	for(ni=0; ni<nStateNum; ni++)
	{
		dTemp = pStationary->pMatElement[ni] 
			+ Tiling_BindingRegionSelection_HMM_GetEmission(pEmission, ni, pScore->pMatElement[0], nEqualLenInterval, dIntS, dIntE, dIntStep)
			+ vBackwardSum[ni]->pMatElement[0]; 
		if(dTemp > dMax)
		{
			dMax = dTemp;
		}
		pDP->pMatElement[ni] = dTemp;
	}

	dSum = 0.0;
	for(ni=0; ni<nStateNum; ni++)
	{
		dTemp = pDP->pMatElement[ni]-dMax;
		dSum += exp(dTemp);
	}
	dSumB = log(dSum)+dMax;

	if(fabs(dSumF-dSumB) > 1e-3)
	{
		printf("Error: TileMap_BindingRegionSelection_HMM, forward and backbward summation not match!\n");
		exit(EXIT_FAILURE);
	}
	dSumF = (dSumF+dSumB)/2.0;

	/* posterior calculation */
	for(nj=0; nj<nProbeNum; nj++)
	{
		dSum = 0.0;
		for(ni=0; ni<nStateNum; ni++)
		{
			dTemp = vForwardSum[ni]->pMatElement[nj]+vBackwardSum[ni]->pMatElement[nj]-dSumF;
			dTemp = exp(dTemp);
			dSum += dTemp;
			vPosterior[ni]->pMatElement[nj] = dTemp;
		}
		for(ni=0; ni<nStateNum; ni++)
		{
			vPosterior[ni]->pMatElement[nj] /= dSum;
		}
	}

	/* destroy space */
	for(ni=0; ni<nStateNum; ni++)
	{
		DestroyDoubleMatrix(vForwardSum[ni]);
		vForwardSum[ni] = NULL;
		DestroyDoubleMatrix(vBackwardSum[ni]);
		vBackwardSum[ni] = NULL;
	}
	free(vForwardSum);
	free(vBackwardSum);

	DestroyDoubleMatrix(pDP);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMap_BindingRegionSelection_MA_UMS_Main()                           */
/*  Moving Average algorithm for calling binding regions in TileMap.       */
/*  FDR will be estimated by UMS.                                          */
/* ----------------------------------------------------------------------- */ 
int TileMap_BindingRegionSelection_MA_UMS_Main(int nProbeNum, 
				char strScorePath[], char strWorkPath[], char strProjectTitle[],
				int nMAW, int nScoreRange,	double dZeroCut, char strUMSSelectPath[], 
				double dTp, double dTq, int nOffset, int nIntervalNum, 
				int nMaxGap, double dFDRCut)
{
	/* define */
	struct DOUBLEMATRIX *pScore;
	struct DOUBLEMATRIX *pSelect;
	struct DOUBLEMATRIX *pMA;
	struct DOUBLEMATRIX *pFDR;

	struct DOUBLEMATRIX *pNewScore;
	struct INTMATRIX *pChr;
	int nChr;
	struct tagString **vChrName;
	struct DOUBLEMATRIX *pPosition;

	/* UMS */
	struct DOUBLEMATRIX *pIntx;
	struct DOUBLEMATRIX *pF0;
	struct DOUBLEMATRIX *pF1;
	struct DOUBLEMATRIX *pLfdr;
	double dSum;
	double dTheta;

	/* for loading data */
	char strFileName[MED_LINE_LENGTH];
	char strLine[LINE_LENGTH];
	char strChr[LINE_LENGTH];
	double dPos,dScore;
	FILE *fpIn;
	FILE *fpOut;
	int ni,nj;

	/* init */
	pMA = NULL;
	pFDR = NULL;

	printf("\n");
	printf("/* ---------------------- */\n");
	printf("/* TileMap MA             */\n");
	printf("/* ---------------------- */\n");

	/* prepare space */
	pScore = NULL;
	pScore = CreateDoubleMatrix(1, nProbeNum);
	if(pScore == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_MA_UMS_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}

	pChr = NULL;
	pChr = CreateIntMatrix(1, nProbeNum);
	if(pChr == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_MA_UMS_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}

	vChrName = NULL;
	vChrName = (struct tagString **)calloc(nProbeNum, sizeof(struct tagString *));
	if(vChrName == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_MA_UMS_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}

	pPosition = NULL;
	pPosition = CreateDoubleMatrix(1, nProbeNum);
	if(pPosition == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_MA_UMS_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}

	/* load data */
	fpIn = NULL;
	fpIn = fopen(strScorePath, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_MA_UMS_Main, cannot open data file!\n");
		exit(EXIT_FAILURE);
	}

	if(fgets(strLine, LINE_LENGTH, fpIn) == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_MA_UMS_Main, empty data file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	nChr = 0;
	while(fgets(strLine, LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%s %lf %lf", strChr, &dPos, &dScore);
		/* sscanf(strLine, "%lf %lf", &dPos, &dScore); */
		pScore->pMatElement[ni] = dScore;
		pPosition->pMatElement[ni] = dPos;
		StringAddTail((vChrName+ni), strChr);

		if(ni == 0)
		{
			pChr->pMatElement[ni] = nChr;
		}
		else
		{
			if(strcmp(vChrName[ni-1]->m_pString, vChrName[ni]->m_pString) == 0)
			{
				pChr->pMatElement[ni] = nChr;
			}
			else
			{
				nChr++;
				pChr->pMatElement[ni] = nChr;
			}
		}
		ni++;
	}

	fclose(fpIn);

	if(ni != nProbeNum)
	{
		printf("Error: TileMap_BindingRegionSelection_MA_UMS_Main, probe number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* score transformation */
	pNewScore = NULL;
	if(nScoreRange == 1)
	{
		pNewScore = Tiling_BindingRegionSelection_ScoreTransform(pScore, dZeroCut, "logit");
		DestroyDoubleMatrix(pScore);
		pScore = pNewScore;
		pNewScore = NULL;
	}

	/* compute MA statistics */
	pMA = NULL;
	pMA = CreateDoubleMatrix(1, nProbeNum);
	if(pMA == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_MA_UMS_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}
	TileMap_MA(pScore, pChr, pPosition, nMAW, pMA);

	/* UMS */
	printf("\n");
	printf("/* ---------------------- */\n");
	printf("/* TileMap UMS            */\n");
	printf("/* ---------------------- */\n");

	pNewScore = NULL;
	pNewScore = Tiling_BindingRegionSelection_ScoreTransform(pMA, dZeroCut, "invlogit");
	DestroyDoubleMatrix(pMA);
	pMA = pNewScore;
	pNewScore = NULL;

	pSelect = NULL;
	pSelect = CreateDoubleMatrix(1, nProbeNum);
	if(pSelect == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_MA_UMS_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}

	fpIn = NULL;
	fpIn = fopen(strUMSSelectPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_MA_UMS_Main, cannot open selection file!\n");
		exit(EXIT_FAILURE);
	}

	if(fgets(strLine, LINE_LENGTH, fpIn) == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_MA_UMS_Main, empty selection file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%s %lf %lf", strChr, &dPos, &dScore);
		pSelect->pMatElement[ni] = dScore;
		if( ((int)(pPosition->pMatElement[ni]) != (int)dPos) || (strcmp(strChr, vChrName[ni]->m_pString) != 0) )
		{
			printf("Error: TileMap_BindingRegionSelection_MA_UMS_Main, probe position not match between selection file and probe summary file!\n");
			exit(EXIT_FAILURE);
		}
		
		ni++;
	}

	fclose(fpIn);

	if(ni != nProbeNum)
	{
		printf("Error: TileMap_BindingRegionSelection_MA_UMS_Main, probe number not match!\n");
		exit(EXIT_FAILURE);
	}

	pIntx = NULL;
	pIntx = CreateDoubleMatrix(1, nIntervalNum);
	if(pIntx == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_MA_UMS_Main, cannot create memory for UMS!\n");
		exit(EXIT_FAILURE);
	}
	pF0 = NULL;
	pF0 = CreateDoubleMatrix(1, nIntervalNum);
	if(pF0 == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_MA_UMS_Main, cannot create memory for UMS!\n");
		exit(EXIT_FAILURE);
	}
	pF1 = NULL;
	pF1 = CreateDoubleMatrix(1, nIntervalNum);
	if(pF1 == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_MA_UMS_Main, cannot create memory for UMS!\n");
		exit(EXIT_FAILURE);
	}
	pLfdr = NULL;
	pLfdr = CreateDoubleMatrix(1, nIntervalNum);
	if(pLfdr == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_MA_UMS_Main, cannot create memory for UMS!\n");
		exit(EXIT_FAILURE);
	}
	
	dTheta = 0.0;

	TileMap_UMS(pSelect, dTp, dTq, nOffset, pMA, nIntervalNum, pIntx, pF0, pF1, &dTheta);

	DestroyDoubleMatrix(pSelect);

	/* get FDR */
	for(ni=0; ni<nIntervalNum; ni++)
	{
		dSum = dTheta*pF1->pMatElement[ni]+(1.0-dTheta)*pF0->pMatElement[ni];
		pLfdr->pMatElement[ni] = 1.0-dTheta*pF1->pMatElement[ni]/dSum;
	}

	/* assign fdr to the original points */
	pFDR = NULL;
	pFDR = CreateDoubleMatrix(1, nProbeNum);
	if(pFDR == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_MA_UMS_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}
	for(nj=0; nj<pScore->nWidth; nj++)
	{
		ni = (int)(pMA->pMatElement[nj]*nIntervalNum);
		if(ni == nIntervalNum)
			ni--;
		pFDR->pMatElement[nj] = pLfdr->pMatElement[ni];
	}

	/* save */
	sprintf(strFileName, "%s%s_ma.sum", strWorkPath, strProjectTitle);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_MA_UMS_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "chromosome\tposition\tMA_stat\tlocal_fdr\n");
	for(ni=0; ni<nProbeNum; ni++)
	{
		/* dScore = log(pMA->pMatElement[ni]/(1.0-pMA->pMatElement[ni])); */
		fprintf(fpOut, "%s\t%d\t%9.7e\t%9.7e\n", vChrName[ni]->m_pString, 
			(int)(pPosition->pMatElement[ni]), pMA->pMatElement[ni], 
			pFDR->pMatElement[ni]);
	}

	fclose(fpOut);

	/* save to *.bed file */
	sprintf(strFileName, "%s%s", strWorkPath, strProjectTitle);
	TileMap_CallBindingRegion_MA(nProbeNum, pFDR,
				pChr, pPosition, vChrName, dFDRCut, nMaxGap,
				strFileName);


	/* release memory */
	DestroyDoubleMatrix(pLfdr);
	DestroyDoubleMatrix(pFDR);
	DestroyDoubleMatrix(pIntx);
	DestroyDoubleMatrix(pF0);
	DestroyDoubleMatrix(pF1);
	DestroyDoubleMatrix(pScore);
	DestroyIntMatrix(pChr);
	DestroyDoubleMatrix(pPosition);
	DestroyDoubleMatrix(pMA);
	for(ni=0; ni<nProbeNum; ni++)
	{
		DeleteString(vChrName[ni]);
		vChrName[ni] = NULL;
	}
	free(vChrName);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMap_MA()                                                           */
/*  Moving Average of TileMap.                                             */
/* ----------------------------------------------------------------------- */ 
int TileMap_MA(struct DOUBLEMATRIX *pScore, struct INTMATRIX * pChr, 
			   struct DOUBLEMATRIX *pPosition, int nW, 
			   struct DOUBLEMATRIX *pMA)
{
	/* define */
	int nPos,nA,nP;
	int ni;
	int nProbeNum;
	double dSum;

	/* init */
	if( (pScore == NULL) || (pChr == NULL) || (pPosition == NULL) || (pMA == NULL) )
	{
		printf("Error: TileMap_MA, incomplete probe level information!\n");
		exit(EXIT_FAILURE);
	}
	if(nW < 0)
	{
		printf("Error: TileMap_MA, W must be >=0!\n");
		exit(EXIT_FAILURE);
	}
	nProbeNum = pScore->nWidth;
	if( (pChr->nWidth != nProbeNum) || (pPosition->nWidth != nProbeNum) || (pMA->nWidth != nProbeNum) )
	{
		printf("Error: number of probes not match!\n");
		exit(EXIT_FAILURE);
	}

	/* if nW == 0 */
	if(nW == 0)
	{
		for(ni=0; ni<nProbeNum; ni++)
		{
			pMA->pMatElement[ni] = pScore->pMatElement[ni];
		}

		return PROC_SUCCESS;
	}

	/* get first position */
	nPos = 0;
	for(ni=1; ni<=nW; ni++)
	{
		if( (nPos+ni) >= nProbeNum)
			break;

		if(pChr->pMatElement[nPos+ni] != pChr->pMatElement[nPos])
		{
			break;
		}
	}
	nP = ni-1;
	nA = 0;

	dSum = 0.0;
	for(ni=0; ni<=nP; ni++)
	{
		dSum += pScore->pMatElement[nPos+ni];
	}
	pMA->pMatElement[nPos] = dSum/(double)(nA+nP+1);

	/* following positions */
	for(nPos=1; nPos<nProbeNum; nPos++)
	{
		/* if new chromosom */
		if(pChr->pMatElement[nPos] != pChr->pMatElement[nPos-1])
		{
			nA = 0;
			for(ni=1; ni<=nW; ni++)
			{
				if( (nPos+ni) >= nProbeNum)
					break;

				if(pChr->pMatElement[nPos+ni] != pChr->pMatElement[nPos])
				{
					break;
				}
			}
			nP = ni-1;

			dSum = 0.0;
			for(ni=0; ni<=nP; ni++)
			{
				dSum += pScore->pMatElement[nPos+ni];
			}
			pMA->pMatElement[nPos] = dSum/(double)(nA+nP+1);
		}

		/* if old chromosome */
		else
		{
			/* subtract the last one */
			if(nA == nW)
			{
				dSum -= pScore->pMatElement[nPos-nA-1];
			}
			else
			{
				nA++;
			}

			/* get the newst one */
			if( (nPos+nP) >= nProbeNum)
			{
				nP--;
			}
			else if(pChr->pMatElement[nPos+nP] != pChr->pMatElement[nPos+nP-1])
			{
				nP--;
			}
			else
			{
				dSum += pScore->pMatElement[nPos+nP];
			}

			if((nA > nW) || (nP > nW))
			{
				printf("Error: TileMap_MA, moving window miscount!\n");
				exit(EXIT_FAILURE);
			}

			pMA->pMatElement[nPos] = dSum/(double)(nA+nP+1);
		}
	}


	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  TileMap_BindingRegionSelection_MA_PERM_Main()                          */
/*  Moving Average algorithm for calling binding regions in TileMap.       */
/*  FDR will be estimated by permutation test.                             */
/* ----------------------------------------------------------------------- */ 
int TileMap_BindingRegionSelection_MA_PERM_Main(int nProbeNum, 
				char strScorePath[], char strWorkPath[], char strProjectTitle[],
				int nW, int nScoreRange, double dZeroCut, int nIntervalNum,
				int nPermNum, int nMaxGap, double dFDRCut)
{
	/* define */
	struct DOUBLEMATRIX *pScore;
	struct DOUBLEMATRIX *pMA;
	struct DOUBLEMATRIX *pFDR;

	struct DOUBLEMATRIX *pPermScore;
	struct DOUBLEMATRIX *pPermMA;
	
	struct DOUBLEMATRIX *pNewScore;
	struct INTMATRIX *pChr;
	int nChr;
	struct tagString **vChrName;
	struct DOUBLEMATRIX *pPosition;

	/* UMS */
	struct DOUBLEMATRIX *pF0;
	struct DOUBLEMATRIX *pF1;
	struct DOUBLEMATRIX *pLfdr;
	double dSum0,dSum1;
	
	/* for loading data */
	char strFileName[MED_LINE_LENGTH];
	char strCommand[MED_LINE_LENGTH];
	char strLine[LINE_LENGTH];
	char strChr[LINE_LENGTH];
	double dPos,dScore;
	FILE *fpIn;
	FILE *fpOut;
	int ni,nj;

	/* init */
	pMA = NULL;
	pFDR = NULL;

	printf("\n");
	printf("/* ---------------------- */\n");
	printf("/* TileMap MA             */\n");
	printf("/* ---------------------- */\n");

	/* prepare space */
	pScore = NULL;
	pScore = CreateDoubleMatrix(1, nProbeNum);
	if(pScore == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_MA_PERM_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}

	pChr = NULL;
	pChr = CreateIntMatrix(1, nProbeNum);
	if(pChr == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_MA_PERM_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}

	vChrName = NULL;
	vChrName = (struct tagString **)calloc(nProbeNum, sizeof(struct tagString *));
	if(vChrName == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_MA_PERM_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}

	pPosition = NULL;
	pPosition = CreateDoubleMatrix(1, nProbeNum);
	if(pPosition == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_MA_PERM_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}

	/* load data */
	fpIn = NULL;
	fpIn = fopen(strScorePath, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_MA_PERM_Main, cannot open data file!\n");
		exit(EXIT_FAILURE);
	}

	if(fgets(strLine, LINE_LENGTH, fpIn) == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_MA_PERM_Main, empty data file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	nChr = 0;
	while(fgets(strLine, LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%s %lf %lf", strChr, &dPos, &dScore);
		/* sscanf(strLine, "%lf %lf", &dPos, &dScore); */
		pScore->pMatElement[ni] = dScore;
		pPosition->pMatElement[ni] = dPos;
		StringAddTail((vChrName+ni), strChr);

		if(ni == 0)
		{
			pChr->pMatElement[ni] = nChr;
		}
		else
		{
			if(strcmp(vChrName[ni-1]->m_pString, vChrName[ni]->m_pString) == 0)
			{
				pChr->pMatElement[ni] = nChr;
			}
			else
			{
				nChr++;
				pChr->pMatElement[ni] = nChr;
			}
		}
		ni++;
	}

	fclose(fpIn);

	if(ni != nProbeNum)
	{
		printf("Error: TileMap_BindingRegionSelection_MA_PERM_Main, probe number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* score transformation */
	pNewScore = NULL;
	if(nScoreRange == 1)
	{
		pNewScore = Tiling_BindingRegionSelection_ScoreTransform(pScore, dZeroCut, "logit");
		DestroyDoubleMatrix(pScore);
		pScore = pNewScore;
		pNewScore = NULL;
	}

	/* compute MA statistics */
	pMA = NULL;
	pMA = CreateDoubleMatrix(1, nProbeNum);
	if(pMA == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_MA_PERM_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}
	TileMap_MA(pScore, pChr, pPosition, nW, pMA);

	pNewScore = NULL;
	pNewScore = Tiling_BindingRegionSelection_ScoreTransform(pMA, dZeroCut, "invlogit");
	DestroyDoubleMatrix(pMA);
	pMA = pNewScore;
	pNewScore = NULL;

	/* Permutation */
	printf("\n");
	printf("/* ---------------------- */\n");
	printf("/* TileMap Permutation    */\n");
	printf("/* ---------------------- */\n");

	pF0 = NULL;
	pF0 = CreateDoubleMatrix(1, nIntervalNum);
	if(pF0 == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_MA_PERM_Main, cannot create memory for UMS!\n");
		exit(EXIT_FAILURE);
	}
	pF1 = NULL;
	pF1 = CreateDoubleMatrix(1, nIntervalNum);
	if(pF1 == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_MA_PERM_Main, cannot create memory for UMS!\n");
		exit(EXIT_FAILURE);
	}
	pLfdr = NULL;
	pLfdr = CreateDoubleMatrix(1, nIntervalNum);
	if(pLfdr == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_MA_PERM_Main, cannot create memory for UMS!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nIntervalNum; ni++)
	{
		pF0->pMatElement[ni] = 1.0;
		pF1->pMatElement[ni] = 1.0;
	}
	
	/* count observed */
	TileMap_CountScore(pMA, nIntervalNum, pF0);
	
	/* count permutation distribution */
	for(nj=0; nj<nPermNum; nj++)
	{
		/* prepare space */
		pPermScore = NULL;
		pPermScore = CreateDoubleMatrix(1, nProbeNum);
		if(pPermScore == NULL)
		{
			printf("Error: TileMap_BindingRegionSelection_MA_PERM_Main, cannot create data space!\n");
			exit(EXIT_FAILURE);
		}

		/* load data */
		sprintf(strFileName, "%s%s_perm%d_f_pb.sum", strWorkPath, strProjectTitle, nj);
		fpIn = NULL;
		fpIn = fopen(strFileName, "r");
		if(fpIn == NULL)
		{
			printf("Error: TileMap_BindingRegionSelection_MA_PERM_Main, cannot open data file!\n");
			exit(EXIT_FAILURE);
		}

		if(fgets(strLine, LINE_LENGTH, fpIn) == NULL)
		{
			printf("Error: TileMap_BindingRegionSelection_MA_PERM_Main, empty data file!\n");
			exit(EXIT_FAILURE);
		}

		ni = 0;
		nChr = 0;
		while(fgets(strLine, LINE_LENGTH, fpIn)!=NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;

			sscanf(strLine, "%s %lf %lf", strChr, &dPos, &dScore);
			/* sscanf(strLine, "%lf %lf", &dPos, &dScore); */
			pPermScore->pMatElement[ni] = dScore;
			if( (int)(pPosition->pMatElement[ni]) != (int)dPos)
			{
				printf("Error: TileMap_BindingRegionSelection_MA_PERM_Main, probe position not match!\n");
				exit(EXIT_FAILURE);
			}
			if( strcmp(strChr, vChrName[ni]->m_pString) != 0)
			{
				printf("Error: TileMap_BindingRegionSelection_MA_PERM_Main, probe position not match!\n");
				exit(EXIT_FAILURE);
			}

			ni++;
		}

		fclose(fpIn);

		if(ni != nProbeNum)
		{
			printf("Error: TileMap_BindingRegionSelection_MA_PERM_Main, probe number not match!\n");
			exit(EXIT_FAILURE);
		}

		/* score transformation */
		pNewScore = NULL;
		if(nScoreRange == 1)
		{
			pNewScore = Tiling_BindingRegionSelection_ScoreTransform(pPermScore, dZeroCut, "logit");
			DestroyDoubleMatrix(pPermScore);
			pPermScore = pNewScore;
			pNewScore = NULL;
		}

		/* compute MA statistics */
		pPermMA = NULL;
		pPermMA = CreateDoubleMatrix(1, nProbeNum);
		if(pPermMA == NULL)
		{
			printf("Error: TileMap_BindingRegionSelection_MA_PERM_Main, cannot create data space!\n");
			exit(EXIT_FAILURE);
		}
		TileMap_MA(pPermScore, pChr, pPosition, nW, pPermMA);

		pNewScore = NULL;
		pNewScore = Tiling_BindingRegionSelection_ScoreTransform(pPermMA, dZeroCut, "invlogit");
		DestroyDoubleMatrix(pPermMA);
		pPermMA = pNewScore;
		pNewScore = NULL;


		/* get statistics */
		TileMap_CountScore(pPermMA, nIntervalNum, pF1);

		DestroyDoubleMatrix(pPermScore);
		DestroyDoubleMatrix(pPermMA);

		/* remove permutation files */
		if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
		{
			sprintf(strCommand, "del %s", strFileName);
			system(strCommand);
		}
		else
		{
			sprintf(strCommand, "rm %s", strFileName);
			system(strCommand);
		}
	}


	/* get fdr */
	dSum0 = 0.0;
	dSum1 = 0.0;
	for(ni=0; ni<nIntervalNum; ni++)
	{
		dSum0 += pF0->pMatElement[ni];
		dSum1 += pF1->pMatElement[ni];
	}
	for(ni=0; ni<nIntervalNum; ni++)
	{
		pF0->pMatElement[ni] = pF0->pMatElement[ni]/dSum0;
		pF1->pMatElement[ni] = pF1->pMatElement[ni]/dSum1;
		pLfdr->pMatElement[ni] = pF1->pMatElement[ni]/pF0->pMatElement[ni];
		if(pLfdr->pMatElement[ni] > 1.0)
			pLfdr->pMatElement[ni] = 1.0;
		
		if(ni>0)
		{
			if(pLfdr->pMatElement[ni] < pLfdr->pMatElement[ni-1])
			{
				pLfdr->pMatElement[ni] = pLfdr->pMatElement[ni-1];
			}
		}
	}

	/* assign fdr to the original points */
	pFDR = NULL;
	pFDR = CreateDoubleMatrix(1, nProbeNum);
	if(pFDR == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_MA_PERM_Main, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}
	for(nj=0; nj<pMA->nWidth; nj++)
	{
		ni = (int)(pMA->pMatElement[nj]*nIntervalNum);
		if(ni == nIntervalNum)
			ni--;
		pFDR->pMatElement[nj] = pLfdr->pMatElement[ni];
	}

	/* save */
	sprintf(strFileName, "%s%s_ma.sum", strWorkPath, strProjectTitle);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: TileMap_BindingRegionSelection_MA_PERM_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "chromosome\tposition\tMA_stat\tlocal_fdr\n");
	for(ni=0; ni<nProbeNum; ni++)
	{
		/* dScore = log(pMA->pMatElement[ni]/(1.0-pMA->pMatElement[ni])); */
		fprintf(fpOut, "%s\t%d\t%9.7e\t%9.7e\n", vChrName[ni]->m_pString, 
			(int)(pPosition->pMatElement[ni]), pMA->pMatElement[ni], 
			pFDR->pMatElement[ni]);
	}

	fclose(fpOut);

	/* save to *.bed file */
	sprintf(strFileName, "%s%s", strWorkPath, strProjectTitle);
	TileMap_CallBindingRegion_MA(nProbeNum, pFDR,
				pChr, pPosition, vChrName, dFDRCut, nMaxGap,
				strFileName);

	/* release memory */
	DestroyDoubleMatrix(pLfdr);
	DestroyDoubleMatrix(pFDR);
	DestroyDoubleMatrix(pF0);
	DestroyDoubleMatrix(pF1);
	DestroyDoubleMatrix(pScore);
	DestroyIntMatrix(pChr);
	DestroyDoubleMatrix(pPosition);
	DestroyDoubleMatrix(pMA);
	for(ni=0; ni<nProbeNum; ni++)
	{
		DeleteString(vChrName[ni]);
		vChrName[ni] = NULL;
	}
	free(vChrName);


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMap_CountScore()                                                   */
/*  Count empirical distributions of scores.                               */
/* ----------------------------------------------------------------------- */ 
int TileMap_CountScore(struct DOUBLEMATRIX *pPermMA, int nIntervalNum, 
					   struct DOUBLEMATRIX *pF)
{
	/* define */
	int ni,nidx;

	/* init */
	if( (pPermMA == NULL) || (pF == NULL) )
	{
		printf("Error: TileMap_CountScore, score does not exist!\n");
		exit(EXIT_FAILURE);
	}
	if(pF->nWidth != nIntervalNum)
	{
		printf("Error: TileMap_CountScore, bin number of empirical distributions not match!\n");
		exit(EXIT_FAILURE);
	}

	/* count */
	for(ni=0; ni<pPermMA->nWidth; ni++)
	{
		nidx = (int)(pPermMA->pMatElement[ni]*nIntervalNum);
		if(nidx == nIntervalNum)
			nidx--;

		pF->pMatElement[nidx] = pF->pMatElement[nidx]+1.0;
	}

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  TileMap_CallBindingRegion_HMM()                                        */
/*  Call binding regions based on HMM posterior probability.               */
/* ----------------------------------------------------------------------- */ 
int TileMap_CallBindingRegion_HMM(int nProbeNum, struct DOUBLEMATRIX *pScore,
				struct INTMATRIX *pChr, struct DOUBLEMATRIX *pPosition, 
				struct tagString **vChrName, double dCutoff, double dGapDist, 
				char strOutPath[])
{
	/* define */
	FILE *fpOut;
	FILE *fpIn;
	FILE *fpBed;
	double dDist;
	int ni,nidx;
	int nStart,nEnd,nCurrentChr;
	int nP1,nP2;
	double dMeanPost;
	double dMinPost;
	int nProbeCount;
	char strChr[LINE_LENGTH];
	char strFileName[MED_LINE_LENGTH];
	char strLine[MED_LINE_LENGTH];
	char strCommand[MED_LINE_LENGTH];
	int nRegionCount;

	/* for resorting regions */
	struct tagString **vSChr;
	struct DOUBLEMATRIX *pSInfo;
	struct DOUBLEMATRIX *pSScore;
	struct DOUBLEMATRIX *pSSortScore;
	struct LONGMATRIX *pSSortIndex;

	/* init */
	if( (pScore->nWidth != nProbeNum) || (pPosition->nWidth != nProbeNum) ||
		(pChr->nWidth != nProbeNum) )
	{
		printf("Error: TileMap_CallBindingRegion_HMM, probe number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* write initial */
	sprintf(strFileName, "%s_hmm.tmpout", strOutPath);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: TileMap_CallBindingRegion_HMM, cannot open *.tmpout file for output!\n");
	}
	
	sprintf(strFileName, "%s_hmm.bed", strOutPath);
	fpBed = NULL;
	fpBed = fopen(strFileName, "w");
	if(fpBed == NULL)
	{
		printf("Error: TileMap_CallBindingRegion_HMM, cannot open *.bed file for output!\n");
	}
	fprintf(fpBed, "browser position chr1:1-1000\n");
	fprintf(fpBed, "track name=TileMap description=\"TileMap track\"\n");


	nRegionCount = 0;
	dDist = 0.0;
	nProbeCount = 0;
	dMeanPost = 0.0;
	dMinPost = 1e20;
	nStart = -1;
	nEnd = -1;
	nP1 = -1;
	nP2 = -1;
	
	if(pScore->pMatElement[0] < dCutoff)
	{
		nCurrentChr = pChr->pMatElement[0];
		nStart = (int)(pPosition->pMatElement[0]);
		nEnd = (int)(pPosition->pMatElement[0]);
		nP1 = 0;
		nP2 = 0;
		nProbeCount = 1;
		dMeanPost = pScore->pMatElement[0];
		dMinPost = pScore->pMatElement[0];
	}

	for(ni=1; ni<nProbeNum; ni++)
	{
		dDist = pPosition->pMatElement[ni]-pPosition->pMatElement[ni-1];
		if(pChr->pMatElement[ni] != pChr->pMatElement[ni-1])
			dDist = dGapDist+1e6;

		if(dDist > dGapDist)
		{
			if(nProbeCount > 0)
			{
				dMeanPost /= (double)nProbeCount;
				fprintf(fpOut, "%s\t%d\t%d\t%d\t%d\t% 9.7e\t% 9.7e\n", 
					vChrName[nP1]->m_pString, nStart, nEnd, (nP1+1), (nP2+1), dMinPost, dMeanPost);
				fprintf(fpBed, "%s\t%d\t%d\ttarget\t%d\t+\n", 
					vChrName[nP1]->m_pString, nStart, nEnd, (int)(1000*(1.0-dMinPost)));
				nRegionCount++;
				nProbeCount = 0;
				dMeanPost = 0.0;
				dMinPost = 1e20;
				nStart = -1;
				nEnd = -1;
				nP1 = -1;
				nP2 = -1;
			}
			else
			{
			}

			if(pScore->pMatElement[ni] < dCutoff)
			{
				nCurrentChr = pChr->pMatElement[ni];
				nStart = (int)(pPosition->pMatElement[ni]);
				nEnd = (int)(pPosition->pMatElement[ni]);
				nP1 = ni;
				nP2 = ni;
				nProbeCount = 1;
				dMeanPost = pScore->pMatElement[ni];
				dMinPost = pScore->pMatElement[ni];
			}
		}
		else
		{
			if(pScore->pMatElement[ni] < dCutoff)
			{
				if(nProbeCount > 0)
				{
					nEnd = (int)(pPosition->pMatElement[ni]);
					nP2 = ni;
					nProbeCount += 1;
					dMeanPost += pScore->pMatElement[ni];
					if(pScore->pMatElement[ni] < dMinPost)
						dMinPost = pScore->pMatElement[ni];
				}
				else
				{
					nCurrentChr = pChr->pMatElement[ni];
					nStart = (int)(pPosition->pMatElement[ni]);
					nEnd = (int)(pPosition->pMatElement[ni]);
					nP1 = ni;
					nP2 = ni;
					nProbeCount = 1;
					dMeanPost = pScore->pMatElement[ni];
					dMinPost = pScore->pMatElement[ni];
				}
			}
			else
			{
				if(nProbeCount > 0)
				{
					dMeanPost /= (double)nProbeCount;
					fprintf(fpOut, "%s\t%d\t%d\t%d\t%d\t% 9.7e\t% 9.7e\n", 
						vChrName[nP1]->m_pString, nStart, nEnd, (nP1+1), (nP2+1), dMinPost, dMeanPost);
					fprintf(fpBed, "%s\t%d\t%d\ttarget\t%d\t+\n", 
						vChrName[nP1]->m_pString, nStart, nEnd, (int)(1000*(1.0-dMinPost)));
					nRegionCount++;
					nProbeCount = 0;
					dMeanPost = 0.0;
					dMinPost = 1e20;
					nStart = -1;
					nEnd = -1;
					nP1 = -1;
					nP2 = -1;
				}
			}
		}
	}

	if(nProbeCount > 0)
	{
		dMeanPost /= (double)nProbeCount;
		fprintf(fpOut, "%s\t%d\t%d\t%d\t%d\t% 9.7e\t% 9.7e\n", 
				vChrName[nP1]->m_pString, nStart, nEnd, (nP1+1), (nP2+1), dMinPost, dMeanPost);
		fprintf(fpBed, "%s\t%d\t%d\ttarget\t%d\t+\n", 
					vChrName[nP1]->m_pString, nStart, nEnd, (int)(1000*(1.0-dMinPost)));
		nRegionCount++;
	}

	/* close file */
	fclose(fpOut);
	fclose(fpBed);

	/* load and resort regions */
	if(nRegionCount <= 0)
	{
		if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
		{
			sprintf(strCommand, "del %s_hmm.tmpout", strOutPath);
			system(strCommand);
		}
		else
		{
			sprintf(strCommand, "rm %s_hmm.tmpout", strOutPath);
			system(strCommand);
		}
		printf("No regions found using the specified criteria!\n");
		return PROC_SUCCESS;
	}
	
	vSChr = NULL;
	vSChr = (struct tagString **)calloc(nRegionCount, sizeof(struct tagString *));
	if(vSChr == NULL)
	{
		printf("Error: TileMap_CallBindingRegion_HMM, cannot re-sort target region\n");
		exit(EXIT_FAILURE);
	}

	pSInfo = NULL;
	pSInfo = CreateDoubleMatrix(nRegionCount, 5);
	if(pSInfo == NULL)
	{
		printf("Error: TileMap_CallBindingRegion_HMM, cannot re-sort target region\n");
		exit(EXIT_FAILURE);
	}

	pSScore = NULL;
	pSScore = CreateDoubleMatrix(1, nRegionCount);
	if(pSScore == NULL)
	{
		printf("Error: TileMap_CallBindingRegion_HMM, cannot re-sort target region\n");
		exit(EXIT_FAILURE);
	}


	/* reload */
	sprintf(strFileName, "%s_hmm.tmpout", strOutPath);
	fpIn = NULL;
	fpIn = fopen(strFileName, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileMap_CallBindingRegion_HMM, cannot open *.tmpout file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%s %d %d %d %d %lf %lf", strChr, &nStart, &nEnd, &nP1, &nP2,
			&dMinPost, &dMeanPost);
		StringAddTail((vSChr+ni), strChr);
		pSScore->pMatElement[ni] = dMinPost;
		DMSETAT(pSInfo, ni, 0, (double)nStart);
		DMSETAT(pSInfo, ni, 1, (double)nEnd);
		DMSETAT(pSInfo, ni, 2, (double)nP1);
		DMSETAT(pSInfo, ni, 3, (double)nP2);
		DMSETAT(pSInfo, ni, 4, dMeanPost);
		ni++;
	}

	fclose(fpIn);

	/* sort */
	pSSortScore = NULL;
	pSSortIndex = NULL;
	DMSORTMERGEA_0(pSScore, &pSSortScore, &pSSortIndex);

	/* rewrite */
	sprintf(strFileName, "%s_hmm.reg", strOutPath);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: TileMap_CallBindingRegion_HMM, cannot open *.reg file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "chromosome\tstart\tend\tstart_line\tend_line\tmax_score\t%mean_score\n");

	for(ni=0; ni<nRegionCount; ni++)
	{
		nidx = pSSortIndex->pMatElement[ni];
		fprintf(fpOut, "%s\t%d\t%d\t%d\t%d\t%9.7e\t%9.7e\n", vSChr[nidx]->m_pString,
			(int)(DMGETAT(pSInfo, nidx, 0)), (int)(DMGETAT(pSInfo, nidx, 1)),
			(int)(DMGETAT(pSInfo, nidx, 2)), (int)(DMGETAT(pSInfo, nidx, 3)),
			1.0-(pSScore->pMatElement[nidx]), DMGETAT(pSInfo, nidx, 4));
	}
	fclose(fpOut);


	/* release memory */
	if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
	{
		sprintf(strCommand, "del %s_hmm.tmpout", strOutPath);
		system(strCommand);
	}
	else
	{
		sprintf(strCommand, "rm %s_hmm.tmpout", strOutPath);
		system(strCommand);
	}

	DestroyDoubleMatrix(pSInfo);
	DestroyDoubleMatrix(pSScore);
	DestroyDoubleMatrix(pSSortScore);
	DestroyLongMatrix(pSSortIndex);
	for(ni=0; ni<nRegionCount; ni++)
	{
		DeleteString(vSChr[ni]);
		vSChr[ni] = NULL;
	}
	free(vSChr);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMap_CallBindingRegion_MA()                                         */
/*  Call binding regions based on MA FDR.                                  */
/* ----------------------------------------------------------------------- */ 
int TileMap_CallBindingRegion_MA(int nProbeNum, struct DOUBLEMATRIX *pScore,
				struct INTMATRIX *pChr, struct DOUBLEMATRIX *pPosition, 
				struct tagString **vChrName, double dCutoff, double dGapDist, 
				char strOutPath[])
{
	/* define */
	FILE *fpOut;
	FILE *fpIn;
	FILE *fpBed;
	double dDist;
	int ni,nidx;
	int nStart,nEnd,nCurrentChr;
	int nLineStart,nLineEnd;
	int nP1,nP2;
	double dP1,dP2;
	double dMeanPost;
	double dMinPost;
	double dMin,dMean;
	int nProbeCount;
	char strChr[LINE_LENGTH];
	char strFileName[MED_LINE_LENGTH];
	char strLine[MED_LINE_LENGTH];
	char strCommand[MED_LINE_LENGTH];
	char strLastChr[LINE_LENGTH];
	int nRegionCount;

	/* for resorting regions */
	struct tagString **vSChr;
	struct DOUBLEMATRIX *pSInfo;
	struct DOUBLEMATRIX *pSScore;
	struct DOUBLEMATRIX *pSSortScore;
	struct LONGMATRIX *pSSortIndex;

	/* init */
	if( (pScore->nWidth != nProbeNum) || (pPosition->nWidth != nProbeNum) ||
		(pChr->nWidth != nProbeNum) )
	{
		printf("Error: TileMap_CallBindingRegion_MA, probe number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* write initial */
	sprintf(strFileName, "%s_ma.tmpout", strOutPath);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: TileMap_CallBindingRegion_MA, cannot open *.tmpout file for output!\n");
	}
	
	nRegionCount = 0;
	dDist = 0.0;
	nProbeCount = 0;
	dMeanPost = 0.0;
	dMinPost = 1e20;
	nStart = -1;
	nEnd = -1;
	nP1 = -1;
	nP2 = -1;
	
	if(pScore->pMatElement[0] < dCutoff)
	{
		nCurrentChr = pChr->pMatElement[0];
		nStart = (int)(pPosition->pMatElement[0]);
		nEnd = (int)(pPosition->pMatElement[0]);
		nP1 = 0;
		nP2 = 0;
		nProbeCount = 1;
		dMeanPost = pScore->pMatElement[0];
		dMinPost = pScore->pMatElement[0];
	}

	for(ni=1; ni<nProbeNum; ni++)
	{
		dDist = pPosition->pMatElement[ni]-pPosition->pMatElement[ni-1];
		if(pChr->pMatElement[ni] != pChr->pMatElement[ni-1])
			dDist = dGapDist+1e6;

		if(dDist > dGapDist)
		{
			if(nProbeCount > 0)
			{
				dMeanPost /= (double)nProbeCount;
				fprintf(fpOut, "%s\t%d\t%d\t%d\t%d\t% 9.7e\t% 9.7e\n", 
					vChrName[nP1]->m_pString, nStart, nEnd, (nP1+1), (nP2+1), dMinPost, dMeanPost);
				
				nRegionCount++;
				nProbeCount = 0;
				dMeanPost = 0.0;
				dMinPost = 1e20;
				nStart = -1;
				nEnd = -1;
				nP1 = -1;
				nP2 = -1;
			}
			else
			{
			}

			if(pScore->pMatElement[ni] < dCutoff)
			{
				nCurrentChr = pChr->pMatElement[ni];
				nStart = (int)(pPosition->pMatElement[ni]);
				nEnd = (int)(pPosition->pMatElement[ni]);
				nP1 = ni;
				nP2 = ni;
				nProbeCount = 1;
				dMeanPost = pScore->pMatElement[ni];
				dMinPost = pScore->pMatElement[ni];
			}
		}
		else
		{
			if(pScore->pMatElement[ni] < dCutoff)
			{
				if(nProbeCount > 0)
				{
					nEnd = (int)(pPosition->pMatElement[ni]);
					nP2 = ni;
					nProbeCount += 1;
					dMeanPost += pScore->pMatElement[ni];
					if(pScore->pMatElement[ni] < dMinPost)
						dMinPost = pScore->pMatElement[ni];
				}
				else
				{
					nCurrentChr = pChr->pMatElement[ni];
					nStart = (int)(pPosition->pMatElement[ni]);
					nEnd = (int)(pPosition->pMatElement[ni]);
					nP1 = ni;
					nP2 = ni;
					nProbeCount = 1;
					dMeanPost = pScore->pMatElement[ni];
					dMinPost = pScore->pMatElement[ni];
				}
			}
			else
			{
				if(nProbeCount > 0)
				{
					dMeanPost /= (double)nProbeCount;
					fprintf(fpOut, "%s\t%d\t%d\t%d\t%d\t% 9.7e\t% 9.7e\n", 
						vChrName[nP1]->m_pString, nStart, nEnd, (nP1+1), (nP2+1), dMinPost, dMeanPost);
					nRegionCount++;
					nProbeCount = 0;
					dMeanPost = 0.0;
					dMinPost = 1e20;
					nStart = -1;
					nEnd = -1;
					nP1 = -1;
					nP2 = -1;
				}
			}
		}
	}

	if(nProbeCount > 0)
	{
		dMeanPost /= (double)nProbeCount;
		fprintf(fpOut, "%s\t%d\t%d\t%d\t%d\t% 9.7e\t% 9.7e\n", 
				vChrName[nP1]->m_pString, nStart, nEnd, (nP1+1), (nP2+1), dMinPost, dMeanPost);
		nRegionCount++;
	}

	/* close file */
	fclose(fpOut);

	/* merge signals */
	if(nRegionCount <= 0)
	{
		/* release memory */
		if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
		{
			sprintf(strCommand, "del %s_ma.tmpout", strOutPath);
			system(strCommand);
		}
		else
		{
			sprintf(strCommand, "rm %s_ma.tmpout", strOutPath);
			system(strCommand);
		}
		printf("No regions found using the specified criteria!\n");
		return PROC_SUCCESS;
	}
	
	vSChr = NULL;
	vSChr = (struct tagString **)calloc(nRegionCount, sizeof(struct tagString *));
	if(vSChr == NULL)
	{
		printf("Error: TileMap_CallBindingRegion_MA, cannot re-sort target region\n");
		exit(EXIT_FAILURE);
	}

	pSInfo = NULL;
	pSInfo = CreateDoubleMatrix(nRegionCount, 6);
	if(pSInfo == NULL)
	{
		printf("Error: TileMap_CallBindingRegion_MA, cannot re-sort target region\n");
		exit(EXIT_FAILURE);
	}

	
	/* reload */
	sprintf(strFileName, "%s_ma.tmpout", strOutPath);
	fpIn = NULL;
	fpIn = fopen(strFileName, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileMap_CallBindingRegion_MA, cannot open *.tmpout file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%s %d %d %d %d %lf %lf", strChr, &nStart, &nEnd, &nP1, &nP2,
			&dMinPost, &dMeanPost);
		StringAddTail((vSChr+ni), strChr);
		DMSETAT(pSInfo, ni, 0, (double)nStart);
		DMSETAT(pSInfo, ni, 1, (double)nEnd);
		DMSETAT(pSInfo, ni, 2, (double)nP1);
		DMSETAT(pSInfo, ni, 3, (double)nP2);
		DMSETAT(pSInfo, ni, 4, dMinPost);
		DMSETAT(pSInfo, ni, 5, dMeanPost);
		ni++;
	}

	fclose(fpIn);

	/* write to bed */
	sprintf(strFileName, "%s_ma.tmpout2", strOutPath);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: TileMap_CallBindingRegion_MA, cannot open *.tmpout file for output!\n");
	}

	sprintf(strFileName, "%s_ma.bed", strOutPath);
	fpBed = NULL;
	fpBed = fopen(strFileName, "w");
	if(fpBed == NULL)
	{
		printf("Error: TileMap_CallBindingRegion_MA, cannot open *.bed file for output!\n");
	}
	
	fprintf(fpBed, "browser position chr1:1-1000\n");
	fprintf(fpBed, "track name=TileMap description=\"TileMap track\"\n");

	nStart = -1000000;
	nEnd = -1000000;
	nLineStart = -1000000;
	nLineEnd = -1000000;
	dMin = 1e25;
	dMean = 1e25;
	strcpy(strLastChr, "chr_NA");

	nRegionCount = 0;
	for(ni=0; ni<pSInfo->nHeight; ni++)
	{
		nP1 = (int)(DMGETAT(pSInfo, ni, 0));
		nP2 = (int)(DMGETAT(pSInfo, ni, 1));

		if( ((double)(nP1-nEnd) < dGapDist) && (strcmp(vSChr[ni]->m_pString, strLastChr) == 0) )
		{
			nEnd = nP2;
			nLineEnd = (int)(DMGETAT(pSInfo, ni, 3));
			dP1 = DMGETAT(pSInfo, ni, 4);
			if(dP1 < dMin);
				dMin = dP1;
			dP2 = DMGETAT(pSInfo, ni, 5);
			if(dP2 < dMean);
				dMean = dP2;
		}
		else
		{
			if(nStart >= 0)
			{
				fprintf(fpOut, "%s\t%d\t%d\t%d\t%d\t% 9.7e\t% 9.7e\n", 
					strLastChr, nStart, nEnd, nLineStart, nLineEnd, dMinPost, dMeanPost);
				fprintf(fpBed, "%s\t%d\t%d\ttarget\t%d\t+\n", 
					strLastChr, nStart, nEnd, (int)(1000*(1.0-dMinPost)));
				nRegionCount++;
			}

			strcpy(strLastChr, vSChr[ni]->m_pString);
			nStart = nP1;
			nEnd = nP2;
			nLineStart = (int)(DMGETAT(pSInfo, ni, 2));
			nLineEnd = (int)(DMGETAT(pSInfo, ni, 3));
			dMin = DMGETAT(pSInfo, ni, 4);
			dMean = DMGETAT(pSInfo, ni, 5);
		}
	}

	if(nStart >= 0)
	{
		fprintf(fpOut, "%s\t%d\t%d\t%d\t%d\t% 9.7e\t% 9.7e\n", 
			strLastChr, nStart, nEnd, nLineStart, nLineEnd, dMinPost, dMeanPost);
		fprintf(fpBed, "%s\t%d\t%d\ttarget\t%d\t+\n", 
			strLastChr, nStart, nEnd, (int)(1000*(1.0-dMinPost)));
		nRegionCount++;
	}

	fclose(fpOut);
	fclose(fpBed);

	/* first round of clearing memory */
	DestroyDoubleMatrix(pSInfo);
	for(ni=0; ni<nRegionCount; ni++)
	{
		DeleteString(vSChr[ni]);
		vSChr[ni] = NULL;
	}
	free(vSChr);

	/* load and resort regions */
	if(nRegionCount <= 0)
	{
		/* release memory */
		if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
		{
			sprintf(strCommand, "del %s_ma.tmpout", strOutPath);
			system(strCommand);
			sprintf(strCommand, "del %s_ma.tmpout2", strOutPath);
			system(strCommand);
		}
		else
		{
			sprintf(strCommand, "rm %s_ma.tmpout", strOutPath);
			system(strCommand);
			sprintf(strCommand, "rm %s_ma.tmpout2", strOutPath);
			system(strCommand);
		}
		printf("No regions found using the specified criteria!\n");
		return PROC_SUCCESS;
	}
	
	vSChr = NULL;
	vSChr = (struct tagString **)calloc(nRegionCount, sizeof(struct tagString *));
	if(vSChr == NULL)
	{
		printf("Error: TileMap_CallBindingRegion_MA, cannot re-sort target region\n");
		exit(EXIT_FAILURE);
	}

	pSInfo = NULL;
	pSInfo = CreateDoubleMatrix(nRegionCount, 5);
	if(pSInfo == NULL)
	{
		printf("Error: TileMap_CallBindingRegion_MA, cannot re-sort target region\n");
		exit(EXIT_FAILURE);
	}

	pSScore = NULL;
	pSScore = CreateDoubleMatrix(1, nRegionCount);
	if(pSScore == NULL)
	{
		printf("Error: TileMap_CallBindingRegion_MA, cannot re-sort target region\n");
		exit(EXIT_FAILURE);
	}


	/* reload */
	sprintf(strFileName, "%s_ma.tmpout2", strOutPath);
	fpIn = NULL;
	fpIn = fopen(strFileName, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileMap_CallBindingRegion_MA, cannot open *.tmpout file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%s %d %d %d %d %lf %lf", strChr, &nStart, &nEnd, &nP1, &nP2,
			&dMinPost, &dMeanPost);
		StringAddTail((vSChr+ni), strChr);
		pSScore->pMatElement[ni] = dMinPost;
		DMSETAT(pSInfo, ni, 0, (double)nStart);
		DMSETAT(pSInfo, ni, 1, (double)nEnd);
		DMSETAT(pSInfo, ni, 2, (double)nP1);
		DMSETAT(pSInfo, ni, 3, (double)nP2);
		DMSETAT(pSInfo, ni, 4, dMeanPost);
		ni++;
	}

	fclose(fpIn);

	/* sort */
	pSSortScore = NULL;
	pSSortIndex = NULL;
	DMSORTMERGEA_0(pSScore, &pSSortScore, &pSSortIndex);

	/* rewrite */
	sprintf(strFileName, "%s_ma.reg", strOutPath);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: TileMap_CallBindingRegion_MA, cannot open *.reg file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "chromosome\tstart\tend\tstart_line\tend_line\tmin_score\t%mean_score\n");

	for(ni=0; ni<nRegionCount; ni++)
	{
		nidx = pSSortIndex->pMatElement[ni];
		fprintf(fpOut, "%s\t%d\t%d\t%d\t%d\t%9.7e\t%9.7e\n", vSChr[nidx]->m_pString,
			(int)(DMGETAT(pSInfo, nidx, 0)), (int)(DMGETAT(pSInfo, nidx, 1)),
			(int)(DMGETAT(pSInfo, nidx, 2)), (int)(DMGETAT(pSInfo, nidx, 3)),
			pSScore->pMatElement[nidx], DMGETAT(pSInfo, nidx, 4));
	}
	fclose(fpOut);


	/* release memory */
	if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
	{
		sprintf(strCommand, "del %s_ma.tmpout", strOutPath);
		system(strCommand);
		sprintf(strCommand, "del %s_ma.tmpout2", strOutPath);
		system(strCommand);
	}
	else
	{
		sprintf(strCommand, "rm %s_ma.tmpout", strOutPath);
		system(strCommand);
		sprintf(strCommand, "rm %s_ma.tmpout2", strOutPath);
		system(strCommand);
	}

	DestroyDoubleMatrix(pSInfo);
	DestroyDoubleMatrix(pSScore);
	DestroyDoubleMatrix(pSSortScore);
	DestroyLongMatrix(pSSortIndex);
	for(ni=0; ni<nRegionCount; ni++)
	{
		DeleteString(vSChr[ni]);
		vSChr[ni] = NULL;
	}
	free(vSChr);

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  TileMap_Extract_Main()                                                 */
/*  TileMap extract data for plot.                                         */
/* ----------------------------------------------------------------------- */ 
int TileMap_Extract_Main(char strParamPath[])
{
	/* ------ */
	/* define */
	/* ------ */

	/* working directory */
	char strWorkPath[MED_LINE_LENGTH];
	char strProjectTitle[LINE_LENGTH];
	char strProbeFile[LINE_LENGTH];
	char strRawData[LINE_LENGTH];
	char strChr[LINE_LENGTH];
	int nStart,nEnd;
	double dP1,dP2;
	FILE *fpIn;
	
	char strLine[MED_LINE_LENGTH];
	int nlen;
	char *chSep;
	int nError=0;

	
	/* --------------- */
	/* load parameters */
	/* --------------- */
	fpIn = NULL;
	fpIn = fopen(strParamPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileMap_Extract_Main, cannot open the parameter file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strstr(strLine, "[Working directory]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strWorkPath, chSep);
			
			if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
			{
				nlen = strlen(strWorkPath);
				if(nlen == 0)
				{
					sprintf(strWorkPath, ".\\"); 
				}
				else if(strWorkPath[nlen-1] != '\\')
				{
					strWorkPath[nlen] = '\\';
					strWorkPath[nlen+1] = '\0';
				}
			}
			else
			{
				nlen = strlen(strWorkPath);
				if(nlen == 0)
				{
					sprintf(strWorkPath, "./"); 
				}
				else if(strWorkPath[nlen-1] != '/')
				{
					strWorkPath[nlen] = '/';
					strWorkPath[nlen+1] = '\0';
				}
			}
		}

		else if(strstr(strLine, "[Project Title]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strProjectTitle, chSep);
		}

		else if(strstr(strLine, "[Probe Level Summary]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strProbeFile, chSep);
		}


		else if(strstr(strLine, "[Raw Data]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strRawData, chSep);
		}

		else if(strstr(strLine, "[Regions]") == strLine)
		{
			break;
		}
	}

	/* get regions */
	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%s %lf %lf", strChr, &dP1, &dP2);
		nStart = (int)dP1;
		nEnd = (int)dP2;
		
		/* extract regions */
		printf("Extract %s:%d-%d\n", strChr, (int)dP1, (int)dP2);
		TileMap_Extract(strWorkPath, strProjectTitle, strProbeFile,
			strRawData, strChr, nStart, nEnd);
	}

	fclose(fpIn);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMap_Extract()                                                      */
/*  TileMap extract data for plot.                                         */
/*  position, probe-sum, hmm, ma-stat, ma-fdr, raw data.                   */
/* ----------------------------------------------------------------------- */ 
int TileMap_Extract(char strWorkPath[], char strProjectTitle[], 
					char strProbeFile[], char strRawData[],
					char strChr[], int nStart, int nEnd)
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	int nProbeNum;
	struct DOUBLEMATRIX *pRInfo;
	char strLine[MED_LINE_LENGTH];
	char strDataLine[LONG_LINE_LENGTH];
	char strFilePath[MED_LINE_LENGTH];
	char strLChr[LINE_LENGTH];
	int nCPos;
	double dP,dScore,dFDR;
	int nP,nL1,nL2,ni,nj;
	int nRegionStart;
	char *chSep,*chSep2;
	struct tagString **vData;

	/* count probes */
	sprintf(strFilePath, "%s%s", strWorkPath, strProbeFile);
	fpIn = NULL;
	fpIn = fopen(strFilePath, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileMap_Extract, cannot open probe summary file!\n");
		exit(EXIT_FAILURE);
	}

	nProbeNum = 0;

	if(fgets(strLine, MED_LINE_LENGTH, fpIn) == NULL)
	{
		printf("Error: TileMap_Extract, empty probe summary file!\n");
		exit(EXIT_FAILURE);
	}

	nRegionStart = 0;
	ni = 0;
	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%s %lf", strLChr, &dP);
		nP = (int)dP;
		
		if(strcmp(strLChr, strChr) == 0)
		{
			if( (nP >= nStart) && (nP <= nEnd) )
			{
				nProbeNum++;

				if(nRegionStart == 0)
				{
					nL1 = ni;
					nRegionStart = 1;
				}
			}
			else
			{
				if(nRegionStart == 1)
				{
					nL2 = ni-1;
					nRegionStart = 2;
					break;
				}
			}
		}

		ni++;

	}

	fclose(fpIn);


	/* create space */
	if(nProbeNum <= 0)
	{
		printf("Warning: TileMap_Extract, no match to %s:%d-%d!\n", strChr, nStart, nEnd);
		return PROC_SUCCESS;
	}

	pRInfo = NULL;
	pRInfo = CreateDoubleMatrix(nProbeNum, 5);
	if(pRInfo == NULL)
	{
		printf("Error: TileMap_Extract, cannot create memory for storing data!\n");
		exit(EXIT_FAILURE);
	}

	
	/* load pb.sum */
	sprintf(strFilePath, "%s%s", strWorkPath, strProbeFile);
	fpIn = NULL;
	fpIn = fopen(strFilePath, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileMap_Extract, cannot open probe summary file!\n");
		exit(EXIT_FAILURE);
	}

	if(fgets(strLine, MED_LINE_LENGTH, fpIn) == NULL)
	{
		printf("Error: TileMap_Extract, empty probe summary file!\n");
		exit(EXIT_FAILURE);
	}

	nRegionStart = 0;
	ni = 0;
	nj = 0;
	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(ni <= nL2)
		{
			if(ni >= nL1)
			{
				sscanf(strLine, "%s %lf %lf", strLChr, &dP, &dScore);
				DMSETAT(pRInfo, nj, 0, dP);
				DMSETAT(pRInfo, nj, 1, dScore);
				nj++;
			}
		}
		else
		{
			 break;
		}

		ni++;
	}

	fclose(fpIn);
	
	if(nj != nProbeNum)
	{
		printf("Error: TileMap_Extract, probe number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* load hmm.sum*/
	sprintf(strFilePath, "%s%s_hmm.sum", strWorkPath, strProjectTitle);
	fpIn = NULL;
	fpIn = fopen(strFilePath, "r");
	if(fpIn != NULL)
	{
		
		if(fgets(strLine, MED_LINE_LENGTH, fpIn) == NULL)
		{
			printf("Warning: TileMap_Extract, empty hmm file!\n");
		}

		else
		{
			nRegionStart = 0;
			ni = 0;
			nj = 0;
			while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
			{
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] == '\0')
					continue;

				if(ni <= nL2)
				{
					if(ni >= nL1)
					{
						sscanf(strLine, "%s %lf %lf", strLChr, &dP, &dScore);
						
						if((int)(DMGETAT(pRInfo, nj, 0)) != (int)dP)
						{
							printf("Error: TileMap_Extract, probe position not match!\n");
							exit(EXIT_FAILURE);
						}

						DMSETAT(pRInfo, nj, 2, dScore);
						nj++;
					}
				}
				else
				{
					 break;
				}

				ni++;
			}

			if(nj != nProbeNum)
			{
				printf("Error: TileMap_Extract, probe number not match!\n");
				exit(EXIT_FAILURE);
			}
		}

		fclose(fpIn);
	}

	/* load ma.sum */
	sprintf(strFilePath, "%s%s_ma.sum", strWorkPath, strProjectTitle);
	fpIn = NULL;
	fpIn = fopen(strFilePath, "r");
	if(fpIn != NULL)
	{
		
		if(fgets(strLine, MED_LINE_LENGTH, fpIn) == NULL)
		{
			printf("Warning: TileMap_Extract, empty ma file!\n");
		}

		else
		{
			nRegionStart = 0;
			ni = 0;
			nj = 0;
			while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
			{
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] == '\0')
					continue;

				if(ni <= nL2)
				{
					if(ni >= nL1)
					{
						sscanf(strLine, "%s %lf %lf %lf", strLChr, &dP, &dScore, &dFDR);
						
						if((int)(DMGETAT(pRInfo, nj, 0)) != (int)dP)
						{
							printf("Error: TileMap_Extract, probe position not match!\n");
							exit(EXIT_FAILURE);
						}

						DMSETAT(pRInfo, nj, 3, dScore);
						DMSETAT(pRInfo, nj, 4, dFDR);
						nj++;
					}
				}
				else
				{
					 break;
				}

				ni++;
			}

			if(nj != nProbeNum)
			{
				printf("Error: TileMap_Extract, probe number not match!\n");
				exit(EXIT_FAILURE);
			}
		}

		fclose(fpIn);
	}

	/* load raw data */
	vData = NULL;
	vData = (struct tagString **)calloc(nProbeNum, sizeof(struct tagString *));
	if(vData == NULL)
	{
		printf("Error: TileMap_Extract, cannot create enough memory for data!\n");
		exit(EXIT_FAILURE);
	}

	sprintf(strFilePath, "%s%s", strWorkPath, strRawData);
	fpIn = NULL;
	fpIn = fopen(strFilePath, "r");
	if(fpIn != NULL)
	{
		
		if(fgets(strDataLine, LONG_LINE_LENGTH, fpIn) == NULL)
		{
			printf("Warning: TileMap_Extract, empty raw data file!\n");
		}

		else
		{
			nRegionStart = 0;
			ni = 0;
			nj = 0;
			nCPos = (int)(DMGETAT(pRInfo, nj, 0));
			while(fgets(strDataLine, LONG_LINE_LENGTH, fpIn) != NULL)
			{
				StrTrimLeft(strDataLine);
				StrTrimRight(strDataLine);
				if(strLine[0] == '\0')
					continue;

				sscanf(strDataLine, "%s %lf", strLChr, &dP);
				nP = (int)dP;
		
				if(strcmp(strLChr, strChr) == 0)
				{
					if( (nP >= nStart) && (nP <= nEnd) )
					{
						if(nRegionStart == 0)
						{
							nRegionStart = 1;
						}

						if(nj == nProbeNum)
							break;

						if((int)(DMGETAT(pRInfo, nj, 0)) == nP)
						{
							chSep = strchr(strDataLine, '\t');
							chSep2 = chSep+1;
							chSep = strchr(chSep2, '\t');
							chSep2 = chSep+1;
							StrTrimLeft(chSep2);
							StringAddTail(vData+nj, chSep2);
							nj++;
						}
						
					}
					else
					{
						if(nRegionStart == 1)
						{
							nRegionStart = 2;
							break;
						}
					}
				}

				ni++;

			}	

			
			if(nj != nProbeNum)
			{
				printf("Error: TileMap_Extract, probe number not match!\n");
				exit(EXIT_FAILURE);
			}
		}

		fclose(fpIn);
	}

	/* save */
	sprintf(strFilePath, "%s%s_%s_%d_%d.txt", strWorkPath, strProjectTitle, strChr, nStart, nEnd);
	fpOut = NULL;
	fpOut = fopen(strFilePath, "w");
	if(fpOut == NULL)
	{
		printf("Error: TileMap_Extract, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nProbeNum; ni++)
	{
		fprintf(fpOut, "%d\t% 9.7e\t% 9.7e\t% 9.7e\t% 9.7e", (int)(DMGETAT(pRInfo, ni, 0)),
			DMGETAT(pRInfo, ni, 1), DMGETAT(pRInfo, ni, 2), DMGETAT(pRInfo, ni, 3),
			DMGETAT(pRInfo, ni, 4));
		if(vData[ni] != NULL)
		{
			fprintf(fpOut, "\t%s", vData[ni]->m_pString);
		}
		fprintf(fpOut, "\n");
	}

	fclose(fpOut);

	/* release memory */
	DestroyDoubleMatrix(pRInfo);
	for(ni=0; ni<nProbeNum; ni++)
	{
		DeleteString(vData[ni]);
		vData[ni] = NULL;
	}
	free(vData);

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  TransLoc_Main()                                                        */
/*  Transloc Main function.                                                */
/*  get clusters of expression units.                                      */
/* ----------------------------------------------------------------------- */ 
int TransLoc_Main(char strParamPath[])
{
	/* define */
	/* parameters */
	char strWorkPath[MED_LINE_LENGTH];
	char strArrayAnnotPath[MED_LINE_LENGTH];
	char strArrayDataPath[MED_LINE_LENGTH];
	char strGenomeAnnotPath[MED_LINE_LENGTH];
	char strOutputFile[MED_LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	char strLine[MED_LINE_LENGTH];

	int nArrayNum;
	double dTruncLow;
	int nTakeLog;
	double dCVCutoff;
	int nGeneNum = 0;
	int nFinalProbeNum = 0;

	/* file */
	FILE *fpIn;
	char *chSep;
	int nError,nlen;

	/* --------------- */
	/* load parameters */
	/* --------------- */
	fpIn = NULL;
	fpIn = fopen(strParamPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: TransLoc_Main, cannot open the parameter file!\n");
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
			
			if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
			{
				nlen = strlen(strWorkPath);
				if(nlen == 0)
				{
					sprintf(strWorkPath, ".\\"); 
				}
				else if(strWorkPath[nlen-1] != '\\')
				{
					strWorkPath[nlen] = '\\';
					strWorkPath[nlen+1] = '\0';
				}
			}
			else
			{
				nlen = strlen(strWorkPath);
				if(nlen == 0)
				{
					sprintf(strWorkPath, "./"); 
				}
				else if(strWorkPath[nlen-1] != '/')
				{
					strWorkPath[nlen] = '/';
					strWorkPath[nlen+1] = '\0';
				}
			}
		}

		else if(strstr(strLine, "[Species]") == strLine)
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
				printf("Error: no species specified!\n");
				nError = 1;
				break;
			}
			else
			{
				strcpy(strSpecies, chSep);
			}
		}

		else if(strstr(strLine, "[Array Number]") == strLine)
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
				nArrayNum = 0;
				printf("Error: zero arrays!\n");
				nError = 1;
				break;
			}
			else
			{
				nArrayNum = atoi(chSep);
			}
		}


		else if(strstr(strLine, "[Array Data]") == strLine)
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
				printf("Error: no input array data!\n");
				nError = 1;
				break;
			}
			else
			{
				strcpy(strArrayDataPath, chSep);
			}
		}

		else if(strstr(strLine, "[Array Annotation]") == strLine)
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
				printf("No input array annotations!\n");
				nError = 1;
				break;
			}
			else
			{
				strcpy(strArrayAnnotPath, chSep);
			}
		}

		else if(strstr(strLine, "[Genome Annotation]") == strLine)
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
				printf("No input genome annotations!\n");
				nError = 1;
				break;
			}
			else
			{
				strcpy(strGenomeAnnotPath, chSep);
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
				printf("No output files!\n");
				nError = 1;
				break;
			}
			else
			{
				strcpy(strOutputFile, chSep);
			}
		}

		else if(strstr(strLine, "[Truncate lower bound]") == strLine)
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
				dTruncLow = 1.0;
				printf("Error: truncation lower bound not set!\n");
				nError = 1;
				break;
			}
			else
			{
				dTruncLow = atof(chSep);
			}
		}

		else if(strstr(strLine, "[Take log2 before calculation (1:yes; 0:no)]") == strLine)
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
				nTakeLog = 0;
			}
			else
			{
				nTakeLog = atoi(chSep);
			}
		}

		else if(strstr(strLine, "[CV cutoff]") == strLine)
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
				dCVCutoff = 0.0;
			}
			else
			{
				dCVCutoff = atof(chSep);
			}
		}

		else if(strLine[0] == '#')
		{
		}

		else
		{
			printf("Error: TransLoc_Main, unknown parameters!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	fclose(fpIn);

	if(nError > 0)
	{
		printf("Error: TransLoc_Main, cannot load the parameter file correctly!\n");
		exit(EXIT_FAILURE);
	}

	/* ---------------------------------- */
	/* get genomic coordinates for probes */
	/* ---------------------------------- */
	/* TransLoc_GetProbeCoordinates(strArrayAnnotPath, strWorkPath, strOutputFile, strSpecies);
	*/

	/* TransLoc_GroupRefGene(strGenomeAnnotPath, strWorkPath, strOutputFile, strSpecies);
	*/
	
	/* ---------------------------------- */
	/* filter out probes                  */
	/* ---------------------------------- */
	TransLoc_FilterProbes(strWorkPath, strOutputFile, strArrayDataPath, nArrayNum,
		dTruncLow, nTakeLog, dCVCutoff, &nGeneNum, &nFinalProbeNum);
	TransLoc_NeighborDistance(strWorkPath, strOutputFile, nArrayNum,
		nGeneNum, nFinalProbeNum);
	

	/* TransLoc_ImportScore(strWorkPath, strOutputFile, strArrayDataPath);
	*/

	/* ---------------------------------- */
	/* sort probes on chromosomes         */
	/* ---------------------------------- */

	/* ---------------------------------- */
	/* compute gene correlations          */
	/* ---------------------------------- */

	/* ---------------------------------- */
	/* get clusters of co-expression and  */
	/* co-localization genes              */
	/* ---------------------------------- */

	/* ---------------------------------- */
	/* Output results                     */
	/* ---------------------------------- */


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TransLoc_GetProbeCoordinates()                                         */
/*  get probe coordinates.                                                 */
/* ----------------------------------------------------------------------- */ 
int TransLoc_GetProbeCoordinates(char strArrayAnnotPath[], char strWorkPath[], 
								 char strOutputFile[], char strSpecies[])
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;

	char strLine[LONG_LINE_LENGTH];
	char strTemp[MEDLONG_LINE_LENGTH];
	char strChrPos[MEDLONG_LINE_LENGTH];
	char strRefLine[MEDLONG_LINE_LENGTH];
	char *chi,*chs;
	char *chSep,*chp;
	char strProbe[LINE_LENGTH];
	char strGenSym[MEDLONG_LINE_LENGTH];
	char strUniGene[LINE_LENGTH];
	char strLocusID[LINE_LENGTH];
	int nLocusID;

	struct tagAffyGenomeAlign *vAlign;
	struct tagAffyGenomeAlign *vTrans;
	struct tagAffyGenomeAlign *pNewAlign,*pPrevAlign;
	struct tagString *vRefId[1024];
	int nRefCount;
	int ni,nLen,nLineId;
	
	
	/* init */
	vAlign = NULL;
	vTrans = NULL;
	fpIn = NULL;
	fpOut = NULL;

	/* open file */
	fpIn = fopen(strArrayAnnotPath,"r");
	if(fpIn == NULL)
	{
		printf("Error: TransLoc_GetProbeCoordinates, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}
	sprintf(strLine, "%s%s.pbcod", strWorkPath, strOutputFile);
	fpOut = fopen(strLine,"w");
	if(fpOut == NULL)
	{
		fclose(fpIn);
		printf("Error: TransLoc_GetProbeCoordinates, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	/* read */
	nLineId = 0;
	fgets(strLine, LONG_LINE_LENGTH, fpIn);
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		nLineId++;

		/* parse */
		vAlign = NULL;
		vTrans = NULL;

		/* 1: probe */
		if(strchr(strLine, '\"') == strLine)
		{
			chs = strLine+1;
			chi = strstr(chs, "\",");
			*chi = '\0';
			chi++;
		}
		else
		{
			chs = strLine;
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		strcpy(strProbe, chs);
		chs = chi+1;
		/* printf("%s\n", strProbe); */
		

		/* 2: GeneChip Array */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 3: Species Scientific Name */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 4: Annotation Date */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 5: Sequence Type */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 6: Sequence Source */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 7: Transcript ID */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 8: Target Description */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;
		

		/* 9: Representative Public ID */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 10-a: Archival UniGene Cluster */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		
		/* 10-b: UniGene ID */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		strcpy(strUniGene, chs);
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 11: Genome Version */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 12: Alignments */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chs++;
			chi = strstr(chs, "\",");
			*chi = '\0';
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
			*chi = '\0';
		}		
		strcpy(strChrPos, chs);
		chs = chi+1;

		if(strcmp(strChrPos, "---") == 0)
		{	
		}
		else
		{
			chp = strChrPos;
			chSep = strstr(chp, "///");
			while(chSep != NULL)
			{
				*chSep = '\0';
				strcpy(strTemp, chp);

				pNewAlign = NULL;
				pNewAlign = AffyGenomeAlignCreate();
				AffyGenomeAlignInit(pNewAlign, strTemp, strSpecies);
				if(vAlign == NULL)
				{
					vAlign = pNewAlign;
					pPrevAlign = pNewAlign;
				}
				else
				{
					pPrevAlign->pNext = pNewAlign;
					pPrevAlign = pNewAlign;
				}

				chp = chSep+3;
				chSep = strstr(chp, "///");
			}


			strcpy(strTemp, chp);
			pNewAlign = NULL;
			pNewAlign = AffyGenomeAlignCreate();
			AffyGenomeAlignInit(pNewAlign, strTemp, strSpecies);
			if(vAlign == NULL)
			{
				vAlign = pNewAlign;
				pPrevAlign = pNewAlign;
			}
			else
			{
				pPrevAlign->pNext = pNewAlign;
				pPrevAlign = pNewAlign;
			}
		}
			

		/* 13: Overlapping Transcripts */
		/* StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chs++;
			chi = strstr(chs, "\",");
			*chi = '\0';
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
			*chi = '\0';
		}
		strcpy(strChrPos, chs);

		chs = chi+1;

		if(strcmp(strChrPos, "---") == 0)
		{	
		}
		else
		{
			chp = strChrPos;
			chSep = strstr(chp, "///");
			while(chSep != NULL)
			{
				*chSep = '\0';
				strcpy(strTemp, chp);

				pNewAlign = NULL;
				pNewAlign = AffyGenomeAlignCreate();
				AffyGenomeTransInit(pNewAlign, strTemp, strSpecies);
				if(vTrans == NULL)
				{
					vTrans = pNewAlign;
					pPrevAlign = pNewAlign;
				}
				else
				{
					pPrevAlign->pNext = pNewAlign;
					pPrevAlign = pNewAlign;
				}

				chp = chSep+3;
				chSep = strstr(chp, "///");
			}


			strcpy(strTemp, chp);
			pNewAlign = NULL;
			pNewAlign = AffyGenomeAlignCreate();
			AffyGenomeTransInit(pNewAlign, strTemp, strSpecies);
			if(vTrans == NULL)
			{
				vTrans = pNewAlign;
				pPrevAlign = pNewAlign;
			}
			else
			{
				pPrevAlign->pNext = pNewAlign;
				pPrevAlign = pNewAlign;
			}
		}
		*/

		/* 14: Gene Title */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 15: Gene Symbol */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chs++;
			chi = strstr(chs, "\",");
			*chi = '\0';
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
			*chi = '\0';
		}
		strcpy(strGenSym, chs);
		chs = chi+1;

		if(strcmp(strGenSym, "---") == 0)
			strcpy(strGenSym, "NA");

		/* 16: Chromosomal Location */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 17: UniGene ClusterType */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chs++;
			chi = strstr(chs, "\",");
			*chi = '\0';
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
			*chi = '\0';
		}
		/* strcpy(strUniGene, chs); */
		chs = chi+1;

		/* 18: Ensembl */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 19: LocusLink or Entrez Gene */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chs++;
			chi = strstr(chs, "\",");
			*chi = '\0';
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
			*chi = '\0';
		}
		strcpy(strLocusID, chs);
		chs = chi+1;

		if(strcmp(strLocusID, "---") == 0)
		{
			nLocusID = -1;
		}
		else if(strstr(strLocusID, "///") != NULL)
		{
			nLocusID = -1;
		}
		else
		{
			nLocusID = atoi(strLocusID);
		}

		/* 20: SwissProt */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 21: EC */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 22: OMIM */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 23: RefSeq Protein ID */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 24: RefSeq Transcript ID */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chs++;
			chi = strstr(chs, "\",");
			*chi = '\0';
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
			*chi = '\0';
		}
		strcpy(strRefLine, chs);
		chs = chi+1;

		nRefCount = 0;
		ni = 0;
		if(strcmp(strRefLine, "---") == 0)
		{	
		}
		else
		{
			chp = strRefLine;
			chSep = strstr(chp, "///");
			while(chSep != NULL)
			{
				*chSep = '\0';
				strcpy(strTemp, chp);
				StrTrimRight(strTemp);
				StrTrimLeft(strTemp);
				nLen = strlen(strTemp)+1;
				vRefId[ni] = CreateString(nLen);
				strcpy(vRefId[ni]->m_pString, strTemp);
				ni++;

				chp = chSep+3;
				chSep = strstr(chp, "///");
			}

			strcpy(strTemp, chp);
			StrTrimRight(strTemp);
			StrTrimLeft(strTemp);
			nLen = strlen(strTemp)+1;
			vRefId[ni] = CreateString(nLen);
			strcpy(vRefId[ni]->m_pString, strTemp);
			ni++;

			nRefCount = ni;
		}

		
		/* 25: FlyBase */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 26: AGI */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 27: WormBase */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 28: MGI Name */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 29: RGD Name */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 30: SGD accession number */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 31: Gene Ontology Biological Process */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 32: Gene Ontology Cellular Component */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 33: Gene Ontology Molecular Function */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 34: Pathway */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 35: Protein Families */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 36: Protein Domains */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 37: InterPro */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 38: Trans Membrane */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 39: QTL */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 40: Annotation Description */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;
		
		/* 41: Annotation Transcript Cluster */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 42: Transcript Assignments */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 43: Annotation Notes */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
		}
		else
		{
			chi = strchr(chs, ',');
		}
		if(chi != NULL)
		{
			printf("Error: Affy_CSVANNOT_To_Reduced_200408, csv annotation loading error!\n");
			exit(EXIT_FAILURE);
		}
		/* strcpy(strTemp, chs); */


		/* write */
		if(vAlign == NULL)
		{
			
		}
		else
		{
			pNewAlign = vAlign;
			while(pNewAlign != NULL)
			{
				if( (pNewAlign->nChr >= 0) && (pNewAlign->dIdentity > 75.0) )
				{
					fprintf(fpOut, "%d\t%d\t%d\t%c\t%f\t%s\t%d\t%s\n",
						pNewAlign->nChr, pNewAlign->nStart, pNewAlign->nEnd-1,
						pNewAlign->chRc, pNewAlign->dIdentity,
						strProbe, nLineId-1, strGenSym);
				}
				pNewAlign = pNewAlign->pNext;
			}
		}

	
		/* destroy */
		while(vAlign != NULL)
		{
			pNewAlign = vAlign;
			vAlign = vAlign->pNext;
			AffyGenomeAlignDestroy(pNewAlign);
		}
		while(vTrans != NULL)
		{
			pNewAlign = vTrans;
			vTrans = vTrans->pNext;
			AffyGenomeAlignDestroy(pNewAlign);
		}
		for(ni=0; ni<nRefCount; ni++)
		{
			DeleteString(vRefId[ni]);
			vRefId[ni] = NULL;
		}
	}


	/* close file */
	fclose(fpIn);
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TransLoc_FilterProbes()                                                */
/*  filter probes and sort them along chromosomes.                         */
/* ----------------------------------------------------------------------- */ 
int TransLoc_FilterProbes(char strWorkPath[], char strOutputFile[], 
						  char strArrayDataPath[], int nArrayNum,
						  double dTruncLow, int nTakeLog, double dCVCutoff,
						  int *pGeneNum, int *pFinalProbeNum)
{
	/* define */
	int nProbeNum;
	struct DOUBLEMATRIX **vData;
	struct tagString **vProbe;
	struct DOUBLEMATRIX *pCV;
	struct DOUBLEMATRIX *pCVSort;

	/* files */
	FILE *fpData;
	FILE *fpIn;
	FILE *fpOut;

	/* strings */
	char strLongLine[LONG_LINE_LENGTH];
	char strLine[MED_LINE_LENGTH];

	/* info */
	int nChr,nStart,nEnd;
	char chStrand;
	double dIdentity;
	char strProbe[LINE_LENGTH];
	int nLineId;
	int nGeneId,nLastGeneId;

	/* others */
	int ni,nj;
	char *chSep,*chSep2;
	double dValue;
	double dCVh;

	/* -------------------------------- */
	/* get probe number                 */
	/* -------------------------------- */
	fpData = NULL;
	fpData = fopen(strArrayDataPath, "r");
	if(fpData == NULL)
	{
		printf("Error: TransLoc_FilterProbes, cannot open array data file!\n");
		exit(EXIT_FAILURE);
	}

	nProbeNum = 0;
	fgets(strLongLine, LONG_LINE_LENGTH, fpData);
	while(fgets(strLongLine, LONG_LINE_LENGTH, fpData) != NULL)
	{
		StrTrimLeft(strLongLine);
		StrTrimRight(strLongLine);
		if(strLongLine[0] == '\0')
			continue;
		nProbeNum++;
	}
	fclose(fpData);

	printf("Probe number = %d\n", nProbeNum);

	/* -------------------------------- */
	/* prepare data space               */
	/* -------------------------------- */
	vData = NULL;
	vData = (struct DOUBLEMATRIX **)calloc(nProbeNum, sizeof(struct DOUBLEMATRIX *));
	if(vData == NULL)
	{
		printf("Error: TransLoc_FilterProbes, cannot allocate memory for bring in data!\n");
		exit(EXIT_FAILURE);
	}

	vProbe = NULL;
	vProbe = (struct tagString **)calloc(nProbeNum, sizeof(struct tagString *));
	if(vProbe == NULL)
	{
		printf("Error: TransLoc_FilterProbes, cannot allocate memory for bring in data!\n");
		exit(EXIT_FAILURE);
	}

	/* -------------------------------- */
	/* load data                        */
	/* -------------------------------- */
	fpData = NULL;
	fpData = fopen(strArrayDataPath, "r");
	if(fpData == NULL)
	{
		printf("Error: TransLoc_FilterProbes, cannot open array data file!\n");
		exit(EXIT_FAILURE);
	}
	
	nj = 0;
	fgets(strLongLine, LONG_LINE_LENGTH, fpData);
	while(fgets(strLongLine, LONG_LINE_LENGTH, fpData) != NULL)
	{
		StrTrimLeft(strLongLine);
		StrTrimRight(strLongLine);
		if(strLongLine[0] == '\0')
			continue;

		if(nj >= nProbeNum)
		{
			printf("Error: TransLoc_FilterProbes, probeset number not match!\n");
			exit(EXIT_FAILURE);
		}

		chSep = strchr(strLongLine, '\t');
		if(chSep != NULL)
			*chSep = '\0';
		StringAddTail(vProbe+nj, strLongLine);

		vData[nj] = NULL;
		vData[nj] = CreateDoubleMatrix(1, nArrayNum);
		if( vData[nj] == NULL)
		{
			printf("Error: TransLoc_FilterProbes, cannot create space to bring in data!\n");
			exit(EXIT_FAILURE);
		}

		ni = 0;
		while(chSep != NULL)
		{
			if( ni >= nArrayNum)
			{
				printf("Error: TransLoc_FilterProbes, array number not match!\n");
				exit(EXIT_FAILURE);
			}

			chSep++;
			chSep2 = strchr(chSep, '\t');
		
			/* middle number */
			if(chSep2 != NULL)
			{
				*chSep2 = '\0';

				if(chSep == chSep2)
				{
					dValue = 0.0;
				}
				else
				{
					dValue = atof(chSep);
				}
			}
			/* last number */
			else
			{
				if(chSep == chSep2)
				{
					dValue = 0.0;
				}
				else
				{
					dValue = atof(chSep);
				}
			}

			if(dValue < dTruncLow)
				dValue = dTruncLow;
			if(nTakeLog == 1)
				dValue = log(dValue)/log(2.0);
		
			vData[nj]->pMatElement[ni] = dValue;
			
			/* get next */
			ni++;
			chSep = chSep2;
		}

		if(ni != nArrayNum)
		{
			printf("Error: TransLoc_FilterProbes, array number not match!\n");
			exit(EXIT_FAILURE);
		}

		/* get next line */
		nj++;
	}

	fclose(fpData);

	if(nj != nProbeNum)
	{
		printf("Error: TransLoc_FilterProbes, probeset number not match!\n");
		exit(EXIT_FAILURE);
	}


	/* -------------------------------- */
	/* compute coefficient of variation */
	/* -------------------------------- */
	pCV = NULL;
	pCV = CreateDoubleMatrix(1, nProbeNum);
	if(pCV == NULL)
	{
		printf("Error: TransLoc_FilterProbes, cannot allocate memory for cv!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nProbeNum; ni++)
	{
		pCV->pMatElement[ni] = TransLoc_ComputeCV(vData[ni]);
	}

	/* DMSAVE(pCV, "shh_cv.txt"); */
	pCVSort = NULL;
	DMSORTMERGEA_0(pCV, &pCVSort, NULL);

	ni = (int)(nProbeNum*dCVCutoff);
	dCVh = pCVSort->pMatElement[ni];
	DestroyDoubleMatrix(pCVSort);

	/* ----------------------------------- */
	/* sort filtered data along chromosome */
	/* ----------------------------------- */
	sprintf(strLine, "%s%s_f.pbcod", strWorkPath, strOutputFile);
	fpIn = NULL;
	fpIn = fopen(strLine, "r");
	if(fpIn == NULL)
	{
		printf("Error: TransLoc_FilterProbes, cannot open coordinates file!\n");
		exit(EXIT_FAILURE);
	}

	*pGeneNum = 0;
	*pFinalProbeNum = 0;

	sprintf(strLine, "%s%s.fdat", strWorkPath, strOutputFile);
	fpOut = NULL;
	fpOut = fopen(strLine, "w");
	if(fpOut == NULL)
	{
		printf("Error: TransLoc_FilterProbes, cannot open file for exporting filtered data!\n");
		exit(EXIT_FAILURE);
	}

	nLastGeneId = -1;
	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%d %d %d %d %c %lf %s %d", &nGeneId, &nChr, &nStart, &nEnd,
					&chStrand, &dIdentity, strProbe, &nLineId);
		if(strcmp(strProbe, vProbe[nLineId]->m_pString) != 0)
		{
			printf("Error: TransLoc_FilterProbes, probe not match!\n");
			exit(EXIT_FAILURE);
		}

		if(pCV->pMatElement[nLineId] >= dCVh)
		{
			fprintf(fpOut, "%d\t%d\t%d\t%d\t%c\t%s", nGeneId, nChr, nStart, nEnd, chStrand, strProbe);
			/* fprintf(fpOut, "%s(%d|chr%d:%d-%d)", strProbe, nGeneId, nChr, nStart, nEnd); */
			for(ni=0; ni<nArrayNum; ni++)
			{
				fprintf(fpOut, "\t%9.7e", vData[nLineId]->pMatElement[ni]);
			}
			fprintf(fpOut, "\n");
			*pFinalProbeNum = *pFinalProbeNum+1;
			if(nGeneId != nLastGeneId)
			{
				*pGeneNum = *pGeneNum+1;
				nLastGeneId = nGeneId;
			}
		}
	}

	fclose(fpIn);
	fclose(fpOut);

	/* -------------------------------- */
	/* release data space               */
	/* -------------------------------- */
	DestroyDoubleMatrix(pCV);

	for(ni=0; ni<nProbeNum; ni++)
	{
		DestroyDoubleMatrix(vData[ni]);
		vData[ni] = NULL;
		DeleteString(vProbe[ni]);
		vProbe[ni] = NULL;
	}
	free(vData);
	free(vProbe);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TransLoc_ComputeCV()                                                   */
/*  Compute coefficient of variation.                                      */
/* ----------------------------------------------------------------------- */ 
double TransLoc_ComputeCV(struct DOUBLEMATRIX *pData)
{
	/* define */
	double dM,dS,dCV,dTemp;
	int nN;
	int ni;

	/* init */
	dM = 0.0;
	dS = 0.0;
	dCV = 0.0;

	if(pData == NULL)
	{
		return dCV;
	}

	nN = pData->nWidth;
	if(nN <= 1)
	{
		return dCV;
	}

	for(ni=0; ni<nN; ni++)
	{
		dM += pData->pMatElement[ni];
	}
	dM /= (double)nN;

	for(ni=0; ni<nN; ni++)
	{
		dTemp = pData->pMatElement[ni]-dM;
		dS += dTemp*dTemp;
	}
	dS /= (double)(nN-1);
	dS = sqrt(dS);

	if( fabs(dM) > 1e-10 )
	{
		dCV = dS/dM;
	}
	else
	{
		dCV = 0.0;
	}

	/* return */
	return dCV;
}

/* ----------------------------------------------------------------------- */ 
/*  TransLoc_ImportScore()                                                 */
/*  import scores and sort them along chromosomes.                         */
/* ----------------------------------------------------------------------- */ 
int TransLoc_ImportScore(char strWorkPath[], char strOutputFile[], char strArrayDataPath[])
{
	/* define */
	int nProbeNum;
	struct DOUBLEMATRIX *pData;
	struct tagString **vProbe;
	int nGeneId,nLastGeneId;
	double dOptScore;
	
	/* files */
	FILE *fpData;
	FILE *fpIn;
	FILE *fpOut;

	/* strings */
	char strLongLine[LONG_LINE_LENGTH];
	char strLine[MED_LINE_LENGTH];

	/* info */
	int nChr,nStart,nEnd;
	char chStrand;
	double dIdentity;
	char strProbe[LINE_LENGTH];
	int nLineId;

	int nLastChr,nLastStart,nLastEnd;
	char chLastStrand;
	double dLastIdentity;
	char strLastProbe[LINE_LENGTH];
	int nLastLineId;

	/* others */
	int ni,nj;
	char *chSep;
	double dValue;
	

	/* -------------------------------- */
	/* get probe number                 */
	/* -------------------------------- */
	fpData = NULL;
	fpData = fopen(strArrayDataPath, "r");
	if(fpData == NULL)
	{
		printf("Error: TransLoc_ImportScore, cannot open array score file!\n");
		exit(EXIT_FAILURE);
	}

	nProbeNum = 0;
	fgets(strLongLine, LONG_LINE_LENGTH, fpData);
	while(fgets(strLongLine, LONG_LINE_LENGTH, fpData) != NULL)
	{
		StrTrimLeft(strLongLine);
		StrTrimRight(strLongLine);
		if(strLongLine[0] == '\0')
			continue;
		nProbeNum++;
	}
	fclose(fpData);

	printf("Probe number = %d\n", nProbeNum);

	/* -------------------------------- */
	/* prepare data space               */
	/* -------------------------------- */
	pData = NULL;
	pData = CreateDoubleMatrix(1, nProbeNum);
	if(pData == NULL)
	{
		printf("Error: TransLoc_ImportScore, cannot allocate memory for bring in data!\n");
		exit(EXIT_FAILURE);
	}

	vProbe = NULL;
	vProbe = (struct tagString **)calloc(nProbeNum, sizeof(struct tagString *));
	if(vProbe == NULL)
	{
		printf("Error: TransLoc_ImportScore, cannot allocate memory for bring in data!\n");
		exit(EXIT_FAILURE);
	}

	/* -------------------------------- */
	/* load data                        */
	/* -------------------------------- */
	fpData = NULL;
	fpData = fopen(strArrayDataPath, "r");
	if(fpData == NULL)
	{
		printf("Error: TransLoc_ImportScore, cannot open array data file!\n");
		exit(EXIT_FAILURE);
	}
	
	nj = 0;
	fgets(strLongLine, LONG_LINE_LENGTH, fpData);
	while(fgets(strLongLine, LONG_LINE_LENGTH, fpData) != NULL)
	{
		StrTrimLeft(strLongLine);
		StrTrimRight(strLongLine);
		if(strLongLine[0] == '\0')
			continue;

		if(nj >= nProbeNum)
		{
			printf("Error: TransLoc_ImportScore, probeset number not match!\n");
			exit(EXIT_FAILURE);
		}

		chSep = strchr(strLongLine, '\t');
		if(chSep != NULL)
			*chSep = '\0';
		StringAddTail(vProbe+nj, strLongLine);

		chSep++;
		sscanf(chSep, "%lf", &dValue);
		pData->pMatElement[nj] = dValue;
			
		/* get next line */
		nj++;
	}

	fclose(fpData);

	if(nj != nProbeNum)
	{
		printf("Error: TransLoc_ImportScore, probeset number not match!\n");
		exit(EXIT_FAILURE);
	}


	/* ----------------------------------- */
	/* sort filtered data along chromosome */
	/* ----------------------------------- */
	sprintf(strLine, "%s%s_f.pbcod", strWorkPath, strOutputFile);
	fpIn = NULL;
	fpIn = fopen(strLine, "r");
	if(fpIn == NULL)
	{
		printf("Error: TransLoc_ImportScore, cannot open coordinates file!\n");
		exit(EXIT_FAILURE);
	}

	sprintf(strLine, "%s%s.fdat", strWorkPath, strOutputFile);
	fpOut = NULL;
	fpOut = fopen(strLine, "w");
	if(fpOut == NULL)
	{
		printf("Error: TransLoc_ImportScore, cannot open file for exporting filtered data!\n");
		exit(EXIT_FAILURE);
	}

	nLastGeneId = -1;
	dOptScore = 0.0;
	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%d %d %d %d %c %lf %s %d", &nGeneId, &nChr, &nStart, &nEnd,
					&chStrand, &dIdentity, strProbe, &nLineId);
		if(strcmp(strProbe, vProbe[nLineId]->m_pString) != 0)
		{
			printf("Error: TransLoc_FilterProbes, probe not match!\n");
			exit(EXIT_FAILURE);
		}

		if(nGeneId != nLastGeneId)
		{
			if(nLastGeneId >= 0)
			{
				fprintf(fpOut, "%d\t%d\t%d\t%f\n", nLastChr, nLastStart, nLastEnd, dOptScore);
			}

			nLastGeneId = nGeneId;
			dOptScore = pData->pMatElement[nLineId];
			nLastChr = nChr;
			nLastStart = nStart;
			nLastEnd = nEnd;
			chLastStrand = chStrand;
			dLastIdentity = dIdentity;
			strcpy(strLastProbe, strProbe);
			nLastLineId = nLineId;
		}
		else
		{
			if(pData->pMatElement[nLineId] < dOptScore)
			{
				nLastGeneId = nGeneId;
				dOptScore = pData->pMatElement[nLineId];
				nLastChr = nChr;
				nLastStart = nStart;
				nLastEnd = nEnd;
				chLastStrand = chStrand;
				dLastIdentity = dIdentity;
				strcpy(strLastProbe, strProbe);
				nLastLineId = nLineId;
			}
		}

		
		/* fprintf(fpOut, "%s(chr%d:%d-%d)", strProbe, nChr, nStart, nEnd); */
	}

	if(nLastGeneId >= 0)
	{
		fprintf(fpOut, "%d\t%d\t%d\t%f\n", nLastChr, nLastStart, nLastEnd, dOptScore);
	}

	fclose(fpIn);
	fclose(fpOut);

	/* -------------------------------- */
	/* release data space               */
	/* -------------------------------- */
	for(ni=0; ni<nProbeNum; ni++)
	{
		DeleteString(vProbe[ni]);
		vProbe[ni] = NULL;
	}
	free(vProbe);
	DestroyDoubleMatrix(pData);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TransLoc_GroupRefGene()                                                */
/*  group refgene annotations.                                             */
/* ----------------------------------------------------------------------- */ 
int TransLoc_GroupRefGene(char strGenomeAnnotPath[], char strWorkPath[], 
						  char strOutputFile[], char strSpecies[])
{
	/* define */

	/* for database */
	FILE *fpRefGene;
	struct tagRefGene *pRefGeneList;
	struct tagRefGene **vRefGene;
	int nRefNum;
	struct INTMATRIX *pGeneId;

	/* for output */
	FILE *fpIn;
	FILE *fpOut;

	/* for affy probes */
	int nProbeNum;
	struct BYTEMATRIX *pMask;
	struct INTMATRIX *pChr;
	struct INTMATRIX *pStart;
	struct INTMATRIX *pEnd;
	struct BYTEMATRIX *pStrand;
	struct DOUBLEMATRIX *pIdentity;
	struct tagString **vProbe;
	struct INTMATRIX *pLineId;
	struct tagString **vGeneName;
	int nChr,nStart,nEnd,nLineId;
	char chStrand;
	double dIdentity;
	char strProbe[LINE_LENGTH];
	char strGenSym[LONG_LINE_LENGTH];
	char *chSep1,*chSep2;
	int nRepeatMask;
	int nDistDiff;

	/* for access */
	char strLine[MED_LINE_LENGTH];
	struct tagRefGene *pRefGene, *pCurrentRefGene;
	char strRefLine[LONG_LINE_LENGTH];
	int ni,nj,nk,nLastId;
	int nCurrentGeneId;
	int nPreGeneNum;
	int vPreGeneId[TRANSLOC_PREGENENUM];
	int nMatched;
	double dMatchScore,dMaxScore;
	int nMaxId;
	int nMinEnd,nMaxStart;

	/* init */
	pRefGeneList = NULL;
	vRefGene = NULL;
	nRefNum = 0;
	
	/* load source refgene */
	fpRefGene = NULL;
	fpRefGene = fopen(strGenomeAnnotPath, "r");
	if(fpRefGene == NULL)
	{
		printf("Error: TransLoc_GroupRefGene, cannot open source refgene file!\n");
		exit(EXIT_FAILURE);
	}

	pCurrentRefGene = NULL;
	while(fgets(strRefLine, LONG_LINE_LENGTH, fpRefGene) != NULL)
	{
		StrTrimLeft(strRefLine);
		StrTrimRight(strRefLine);
		if(strRefLine[0] == '\0')
			continue;

		pRefGene = NULL;
		pRefGene = RefGeneCreate();
		RefGeneInit_FromGenomeLabFormat(pRefGene, strRefLine, strSpecies);
		if(pRefGeneList == NULL)
		{
			pRefGeneList = pRefGene;
			pCurrentRefGene = pRefGene;
		}
		else
		{
			pCurrentRefGene->pNext = pRefGene;
			pCurrentRefGene = pRefGene;
		}
		nRefNum++;
	}

	fclose(fpRefGene);

	vRefGene = (struct tagRefGene **)calloc(nRefNum, sizeof(struct tagRefGene*));
	if(vRefGene == NULL)
	{
		printf("Error: TransLoc_GroupRefGene, organize refgene into array!\n");
		exit(EXIT_FAILURE);
	}

	ni=0;
	while(pRefGeneList != NULL)
	{
		pRefGene = pRefGeneList;
		pRefGeneList = pRefGene->pNext;
		pRefGene->pNext = NULL;
		vRefGene[ni] = pRefGene;
		ni++;
	}

	if(ni != nRefNum)
	{
		printf("Error: TransLoc_GroupRefGene, refgene number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* init gene id */
	pGeneId = NULL;
	pGeneId = CreateIntMatrix(1, nRefNum);
	if(pGeneId == NULL)
	{
		printf("Error: TransLoc_GroupRefGene, cannot create space for gene id!\n");
		exit(EXIT_FAILURE);
	}

	/* assign gene id */
	pGeneId->pMatElement[0] = 0;
	nCurrentGeneId = 0;
	for(ni=1; ni<nRefNum; ni++)
	{
		nPreGeneNum = 0;
		dMaxScore = 0.0;
		nMaxId = -1;
		for(nj=ni-1; nj>=0; nj--)
		{
			nMatched = 0;
			for(nk=nPreGeneNum-1; nk>=0; nk--)
			{
				if(pGeneId->pMatElement[nj] == vPreGeneId[nk])
				{
					nMatched = 1;
					break;
				}
			}
			if(nMatched == 0)
			{
				if(nPreGeneNum >= TRANSLOC_PREGENENUM)
				{
					break;
				}
				vPreGeneId[nPreGeneNum] = pGeneId->pMatElement[nj];
				nPreGeneNum++;
			}

			dMatchScore = RefGene_Match(vRefGene[ni], vRefGene[nj]);
			if(dMatchScore > dMaxScore)
			{
				dMaxScore = dMatchScore;
				nMaxId = pGeneId->pMatElement[nj];
			}
		}

		if(dMaxScore > TRANSLOC_REFGENEMATCH_TH)
		{
			pGeneId->pMatElement[ni] = nMaxId;
		}
		else
		{
			nCurrentGeneId++;
			pGeneId->pMatElement[ni] = nCurrentGeneId;
		}
	}

	/* save gene id */
	sprintf(strLine, "%s%s.gid", strWorkPath, strOutputFile);
	fpOut = NULL;
	fpOut = fopen(strLine, "w");
	if(fpOut == NULL)
	{
		printf("Error: TransLoc_GroupRefGene, cannot open file to save the gene-id map!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nRefNum; ni++)
	{
		fprintf(fpOut, "%d\t", pGeneId->pMatElement[ni]);
		RefGeneWrite(vRefGene[ni], fpOut);
	}

	fclose(fpOut);

	/* --------------------- */
	/* load affy coordinates */
	/* --------------------- */

	/* get probe number */
	sprintf(strLine, "%s%s_sorted.pbcod", strWorkPath, strOutputFile);
	fpIn = NULL;
	fpIn = fopen(strLine, "r");
	if(fpIn == NULL)
	{
		printf("Error: TransLoc_GroupRefGene, cannot open file to load affy probes!\n");
		exit(EXIT_FAILURE);
	}

	nProbeNum = 0;
	while(fgets(strRefLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strRefLine);
		StrTrimRight(strRefLine);
		if(strRefLine[0] == '\0')
			continue;

		nProbeNum++;
	}

	fclose(fpIn);

	/* init space */
	pMask = NULL;
	pMask = CreateByteMatrix(1, nProbeNum);
	if(pMask == NULL)
	{
		printf("Error: TransLoc_GroupRefGene, cannot create memory for loading affy data!\n");
		exit(EXIT_FAILURE);
	}

	pChr = NULL;
	pChr = CreateIntMatrix(1, nProbeNum);
	if(pChr == NULL)
	{
		printf("Error: TransLoc_GroupRefGene, cannot create memory for loading affy data!\n");
		exit(EXIT_FAILURE);
	}

	pStart = NULL;
	pStart = CreateIntMatrix(1, nProbeNum);
	if(pStart == NULL)
	{
		printf("Error: TransLoc_GroupRefGene, cannot create memory for loading affy data!\n");
		exit(EXIT_FAILURE);
	}

	pEnd = NULL;
	pEnd = CreateIntMatrix(1, nProbeNum);
	if(pEnd == NULL)
	{
		printf("Error: TransLoc_GroupRefGene, cannot create memory for loading affy data!\n");
		exit(EXIT_FAILURE);
	}

	pStrand = NULL;
	pStrand = CreateByteMatrix(1, nProbeNum);
	if(pStrand == NULL)
	{
		printf("Error: TransLoc_GroupRefGene, cannot create memory for loading affy data!\n");
		exit(EXIT_FAILURE);
	}

	pIdentity = NULL;
	pIdentity = CreateDoubleMatrix(1, nProbeNum);
	if(pIdentity == NULL)
	{
		printf("Error: TransLoc_GroupRefGene, cannot create memory for loading affy data!\n");
		exit(EXIT_FAILURE);
	}

	pLineId = NULL;
	pLineId = CreateIntMatrix(1, nProbeNum);
	if(pLineId == NULL)
	{
		printf("Error: TransLoc_GroupRefGene, cannot create memory for loading affy data!\n");
		exit(EXIT_FAILURE);
	}

	vProbe = NULL;
	vProbe = (struct tagString **)calloc(nProbeNum, sizeof(struct tagString *));
	if(vProbe == NULL)
	{
		printf("Error: TransLoc_GroupRefGene, cannot create memory for loading affy data!\n");
		exit(EXIT_FAILURE);
	}

	vGeneName = NULL;
	vGeneName = (struct tagString **)calloc(nProbeNum, sizeof(struct tagString *));
	if(vGeneName == NULL)
	{
		printf("Error: TransLoc_GroupRefGene, cannot create memory for loading affy data!\n");
		exit(EXIT_FAILURE);
	}


	/* load */
	sprintf(strLine, "%s%s_sorted.pbcod", strWorkPath, strOutputFile);
	fpIn = NULL;
	fpIn = fopen(strLine, "r");
	if(fpIn == NULL)
	{
		printf("Error: TransLoc_GroupRefGene, cannot open file to load affy probes!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strRefLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strRefLine);
		StrTrimRight(strRefLine);
		if(strRefLine[0] == '\0')
			continue;

		chSep1 = strchr(strRefLine, '\t');
		*chSep1 = '\0';
		nChr = atoi(strRefLine);

		chSep1++;
		chSep2 = strchr(chSep1, '\t');
		*chSep2 = '\0';
		nStart = atoi(chSep1);

		chSep1 = chSep2+1;
		chSep2 = strchr(chSep1, '\t');
		*chSep2 = '\0';
		nEnd = atoi(chSep1);

		chSep1 = chSep2+1;
		chSep2 = strchr(chSep1, '\t');
		*chSep2 = '\0';
		StrTrimLeft(chSep1);
		chStrand = *chSep1;

		chSep1 = chSep2+1;
		chSep2 = strchr(chSep1, '\t');
		*chSep2 = '\0';
		dIdentity = atof(chSep1);

		chSep1 = chSep2+1;
		chSep2 = strchr(chSep1, '\t');
		*chSep2 = '\0';
		strcpy(strProbe, chSep1);
		StrTrimLeft(strProbe);
		StrTrimRight(strProbe);
		
		
		chSep1 = chSep2+1;
		chSep2 = strchr(chSep1, '\t');
		*chSep2 = '\0';
		nLineId = atoi(chSep1);

		
		chSep1 = chSep2+1;
		strcpy(strGenSym, chSep1);
		StrTrimLeft(strGenSym);
		StrTrimRight(strGenSym);

		pChr->pMatElement[ni] = nChr;
		pStart->pMatElement[ni] = nStart;
		pEnd->pMatElement[ni] = nEnd;
		if(chStrand == '-')
			pStrand->pMatElement[ni] = 1;
		else
			pStrand->pMatElement[ni] = 0;

		pIdentity->pMatElement[ni] = dIdentity;
		pLineId->pMatElement[ni] = nLineId;

		StringAddTail(vProbe+ni, strProbe);
		StringAddTail(vGeneName+ni, strGenSym);

		ni++;
	}

	fclose(fpIn);
	if(ni != nProbeNum)
	{
		printf("Error: TransLoc_GroupRefGene, probe number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* mask repeated probes */
	for(ni=0; ni<nProbeNum; ni++)
	{
		nRepeatMask = 0;
		for(nj=ni-1; nj>=0; nj--)
		{
			if(pChr->pMatElement[nj] != pChr->pMatElement[ni])
			{
				break;
			}
			nDistDiff = (pStart->pMatElement[ni]+pEnd->pMatElement[ni])/2-(pStart->pMatElement[nj]+pEnd->pMatElement[nj])/2;
			if(nDistDiff > TRANSLOC_REPEATMASK_WIN)
			{
				break;
			}
			if(strcmp(vProbe[ni]->m_pString, vProbe[nj]->m_pString) == 0)
			{
				nRepeatMask = 1;
				break;
			}
		}

		if(nRepeatMask == 1)
		{
			pMask->pMatElement[ni] = 1;
			continue;
		}

		for(nj=ni+1; nj<nProbeNum; nj++)
		{
			if(pChr->pMatElement[nj] != pChr->pMatElement[ni])
			{
				break;
			}
			nDistDiff = (pStart->pMatElement[nj]+pEnd->pMatElement[nj])/2-(pStart->pMatElement[ni]+pEnd->pMatElement[ni])/2;
			if(nDistDiff > TRANSLOC_REPEATMASK_WIN)
			{
				break;
			}
			if(strcmp(vProbe[ni]->m_pString, vProbe[nj]->m_pString) == 0)
			{
				nRepeatMask = 1;
				break;
			}
		}

		if(nRepeatMask == 1)
		{
			pMask->pMatElement[ni] = 1;
		}
	}
	

	/* map affy probes to gene id */
	sprintf(strLine, "%s%s_f.pbcod", strWorkPath, strOutputFile);
	fpOut = NULL;
	fpOut = fopen(strLine, "w");
	if(fpOut == NULL)
	{
		printf("Error: TransLoc_GroupRefGene, cannot open file to write affy probes!\n");
		exit(EXIT_FAILURE);
	}

	nLastId = 0;
	for(ni=0; ni<nProbeNum; ni++)
	{
		if(pMask->pMatElement[ni] == 1)
			continue;

		dMaxScore = 0.0;
		nMaxId = -1;
		
		if(pStrand->pMatElement[ni] == 1)
		{
			chStrand = '-';
		}
		else
		{
			chStrand = '+';
		}

		for(nj=nLastId; nj>=0; nj--)
		{
			if(vRefGene[nj]->nChrom > pChr->pMatElement[ni])
				continue;
			else if(vRefGene[nj]->nChrom < pChr->pMatElement[ni])
				break;
			
			if(vRefGene[nj]->nTxEnd < pStart->pMatElement[ni])
				break;
			if(vRefGene[nj]->nTxStart > pEnd->pMatElement[ni])
				continue;

			nMaxStart = pStart->pMatElement[ni];
			if(vRefGene[nj]->nTxStart > nMaxStart)
				nMaxStart = vRefGene[nj]->nTxStart;

			nMinEnd = pEnd->pMatElement[ni];
			if(vRefGene[nj]->nTxEnd < nMinEnd)
				nMinEnd = vRefGene[nj]->nTxEnd;

			if(nMinEnd < nMaxStart)
				break;

			if(chStrand == vRefGene[nj]->chStrand)	
				dMatchScore = (double)(nMinEnd-nMaxStart+1)/(double)(vRefGene[nj]->nTxEnd-vRefGene[nj]->nTxStart+1);
			else
				dMatchScore = 0.0;

			if(dMatchScore > dMaxScore)
			{
				dMaxScore = dMatchScore;
				nMaxId = nj;
			}
		}

		for(nj=nLastId+1; nj<nRefNum; nj++)
		{
			if(vRefGene[nj]->nChrom < pChr->pMatElement[ni])
				continue;
			else if(vRefGene[nj]->nChrom > pChr->pMatElement[ni])
				break;
			if(vRefGene[nj]->nTxStart > pEnd->pMatElement[ni])
				break;
			if(vRefGene[nj]->nTxEnd < pStart->pMatElement[ni])
				continue;

			nMaxStart = pStart->pMatElement[ni];
			if(vRefGene[nj]->nTxStart > nMaxStart)
				nMaxStart = vRefGene[nj]->nTxStart;

			nMinEnd = pEnd->pMatElement[ni];
			if(vRefGene[nj]->nTxEnd < nMinEnd)
				nMinEnd = vRefGene[nj]->nTxEnd;

			if(nMinEnd < nMaxStart)
				continue;

			if(chStrand == vRefGene[nj]->chStrand)	
				dMatchScore = (double)(nMinEnd-nMaxStart+1)/(double)(vRefGene[nj]->nTxEnd-vRefGene[nj]->nTxStart+1);
			else
				dMatchScore = 0.0;

			if(dMatchScore > dMaxScore)
			{
				dMaxScore = dMatchScore;
				nMaxId = nj;
			}
		}

		if(dMaxScore > TRANSLOC_REFGENEMATCH_TH)
		{
			nLastId = nMaxId;
			fprintf(fpOut, "%d\t%d\t%d\t%d\t%c\t%f\t%s\t%d\t%s\n",
				pGeneId->pMatElement[nMaxId], pChr->pMatElement[ni],
				pStart->pMatElement[ni], pEnd->pMatElement[ni],
				chStrand, pIdentity->pMatElement[ni], 
				vProbe[ni]->m_pString, pLineId->pMatElement[ni],
				vGeneName[ni]->m_pString);
		}
	}

	fclose(fpOut);

	/* release memory */
	DestroyByteMatrix(pMask);
	DestroyIntMatrix(pChr);
	DestroyIntMatrix(pStart);
	DestroyIntMatrix(pEnd);
	DestroyByteMatrix(pStrand);
	DestroyDoubleMatrix(pIdentity);
	DestroyIntMatrix(pLineId);
	
	for(ni=0; ni<nProbeNum; ni++)
	{
		DeleteString(vProbe[ni]);
		DeleteString(vGeneName[ni]);
	}
	free(vProbe);
	free(vGeneName);

	for(ni=0; ni<nRefNum; ni++)
	{
		RefGeneDestroy(vRefGene[ni]);
		vRefGene[ni] = NULL;
	}
	free(vRefGene);
	DestroyIntMatrix(pGeneId);

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  TransLoc_NeighborDistance()                                            */
/*  compute distance of neighboring genes.                                 */
/* ----------------------------------------------------------------------- */ 
int TransLoc_NeighborDistance(char strWorkPath[], char strOutputFile[], 
		int nArrayNum, int nGeneNum, int nProbeNum)
{
	/* define */
	struct INTMATRIX *pGeneId;
	struct INTMATRIX *pGeneP0;
	struct INTMATRIX *pGeneP1;
	struct DOUBLEMATRIX **vData;
	struct INTMATRIX *pChr;
	struct INTMATRIX *pStart;
	struct INTMATRIX *pEnd;
	struct BYTEMATRIX *pStrand;
	struct tagString **vProbe;
	int ni,nj,nk,nx,ny,nz;
	int nGeneCorrHalfWin = TRANSLOC_GENECORR_WIN/2;

	/* data file */
	FILE *fpData;
	char strLongLine[LONG_LINE_LENGTH];
	char strLine[MED_LINE_LENGTH];
	char *chSep,*chSep2;
	int nLastGeneId,nGeneId;
	int nGP0,nGP1;
	char strProbe[LINE_LENGTH];
	double dValue;

	/* correlation distance */
	struct DOUBLEMATRIX **vCorr;
	/* genomic distance */
	struct DOUBLEMATRIX **vDist;
	/* temp corr */
	struct DOUBLEMATRIX *pTempCorr;
	struct DOUBLEMATRIX *pTempDist;

	struct DOUBLEMATRIX *pFCorr;
	struct DOUBLEMATRIX *pFDist;
	struct DOUBLEMATRIX *pSortCorr;
	struct LONGMATRIX *pSortIndex;

	int nS0,nS1,nST;
	double dCorr,dDist,dTemp;
	int nEffectCorrNum;
	

	/* create space */
	pGeneId = NULL;
	pGeneId = CreateIntMatrix(1, nGeneNum);
	if(pGeneId == NULL)
	{
		printf("Error: TransLoc_NeighborDistance, cannot create space for loading data");
		exit(EXIT_FAILURE);
	}
	pGeneP0 = NULL;
	pGeneP0 = CreateIntMatrix(1, nGeneNum);
	if(pGeneP0 == NULL)
	{
		printf("Error: TransLoc_NeighborDistance, cannot create space for loading data");
		exit(EXIT_FAILURE);
	}
	pGeneP1 = NULL;
	pGeneP1 = CreateIntMatrix(1, nGeneNum);
	if(pGeneP1 == NULL)
	{
		printf("Error: TransLoc_NeighborDistance, cannot create space for loading data");
		exit(EXIT_FAILURE);
	}

	pChr = NULL;
	pChr = CreateIntMatrix(1, nProbeNum);
	if(pChr == NULL)
	{
		printf("Error: TransLoc_NeighborDistance, cannot create space for loading data");
		exit(EXIT_FAILURE);
	}

	pStart = NULL;
	pStart = CreateIntMatrix(1, nProbeNum);
	if(pStart == NULL)
	{
		printf("Error: TransLoc_NeighborDistance, cannot create space for loading data");
		exit(EXIT_FAILURE);
	}

	pEnd = NULL;
	pEnd = CreateIntMatrix(1, nProbeNum);
	if(pEnd == NULL)
	{
		printf("Error: TransLoc_NeighborDistance, cannot create space for loading data");
		exit(EXIT_FAILURE);
	}

	pStrand = NULL;
	pStrand = CreateByteMatrix(1, nProbeNum);
	if(pStrand == NULL)
	{
		printf("Error: TransLoc_NeighborDistance, cannot create space for loading data");
		exit(EXIT_FAILURE);
	}

	vProbe = NULL;
	vProbe = (struct tagString **)calloc(nProbeNum, sizeof(struct tagString*));
	if(vProbe == NULL)
	{
		printf("Error: TransLoc_NeighborDistance, cannot create space for loading data");
		exit(EXIT_FAILURE);
	}

	vData = NULL;
	vData = (struct DOUBLEMATRIX **)calloc(nProbeNum, sizeof(struct DOUBLEMATRIX*));
	if(vData == NULL)
	{
		printf("Error: TransLoc_NeighborDistance, cannot create space for loading data");
		exit(EXIT_FAILURE);
	}

	/* ------------------- */
	/* load data from file */
	/* ------------------- */
	sprintf(strLine, "%s%s.fdat", strWorkPath, strOutputFile);
	fpData = NULL;
	fpData = fopen(strLine, "r");
	if(fpData == NULL)
	{
		printf("Error: TransLoc_NeighborDistance, cannot open array data file!\n");
		exit(EXIT_FAILURE);
	}
	
	nj = 0;
	nk = 0;
	nLastGeneId = -1;
	while(fgets(strLongLine, LONG_LINE_LENGTH, fpData) != NULL)
	{
		StrTrimLeft(strLongLine);
		StrTrimRight(strLongLine);
		if(strLongLine[0] == '\0')
			continue;

		if(nj >= nProbeNum)
		{
			printf("Error: TransLoc_NeighborDistance, probeset number not match!\n");
			exit(EXIT_FAILURE);
		}

		chSep = strchr(strLongLine, '\t');
		*chSep = '\0';
		nGeneId = atoi(strLongLine);
		if(nGeneId != nLastGeneId)
		{
			if(nLastGeneId >= 0)
			{
				pGeneId->pMatElement[nk] = nLastGeneId;
				pGeneP0->pMatElement[nk] = nGP0;
				pGeneP1->pMatElement[nk] = nGP1;
				nk++;
			}
			nLastGeneId = nGeneId;
			nGP0 = nj;
			nGP1 = nj;
		}
		else
		{
			nGP1 = nj;
		}

		chSep++;
		chSep2 = strchr(chSep, '\t');
		*chSep2 = '\0';
		pChr->pMatElement[nj] = atoi(chSep);

		chSep = chSep2+1;
		chSep2 = strchr(chSep, '\t');
		*chSep2 = '\0';
		pStart->pMatElement[nj] = atoi(chSep);

		chSep = chSep2+1;
		chSep2 = strchr(chSep, '\t');
		*chSep2 = '\0';
		pEnd->pMatElement[nj] = atoi(chSep);

		chSep = chSep2+1;
		chSep2 = strchr(chSep, '\t');
		*chSep2 = '\0';
		StrTrimLeft(chSep);
		if(*chSep == '-')
			pStrand->pMatElement[nj] = 1;
		else
			pStrand->pMatElement[nj] = 0;

		chSep = chSep2+1;
		chSep2 = strchr(chSep, '\t');
		*chSep2 = '\0';
		strcpy(strProbe, chSep);
		StrTrimLeft(strProbe);
		StrTrimRight(strProbe);
		StringAddTail(vProbe+nj, strProbe);

		vData[nj] = NULL;
		vData[nj] = CreateDoubleMatrix(1, nArrayNum);
		if( vData[nj] == NULL)
		{
			printf("Error: TransLoc_NeighborDistance, cannot create space to bring in data!\n");
			exit(EXIT_FAILURE);
		}

		ni = 0;
		while(chSep != NULL)
		{
			if( ni >= nArrayNum)
			{
				printf("Error: TransLoc_NeighborDistance, array number not match!\n");
				exit(EXIT_FAILURE);
			}

			chSep = chSep2+1;
			chSep2 = strchr(chSep, '\t');
		
			/* middle number */
			if(chSep2 != NULL)
			{
				*chSep2 = '\0';

				if(chSep == chSep2)
				{
					dValue = 0.0;
				}
				else
				{
					dValue = atof(chSep);
				}
			}
			/* last number */
			else
			{
				if(chSep == chSep2)
				{
					dValue = 0.0;
				}
				else
				{
					dValue = atof(chSep);
				}
			}

			vData[nj]->pMatElement[ni] = dValue;
			
			/* get next */
			ni++;
			chSep = chSep2;
		}

		if(ni != nArrayNum)
		{
			printf("Error: TransLoc_NeighborDistance, array number not match!\n");
			exit(EXIT_FAILURE);
		}

		/* get next line */
		nj++;
	}

	if(nLastGeneId >= 0)
	{
		pGeneId->pMatElement[nk] = nLastGeneId;
		pGeneP0->pMatElement[nk] = nGP0;
		pGeneP1->pMatElement[nk] = nGP1;
		nk++;
	}

	fclose(fpData);

	if(nj != nProbeNum)
	{
		printf("Error: TransLoc_NeighborDistance, probeset number not match!\n");
		exit(EXIT_FAILURE);
	}

	if(nk != nGeneNum)
	{
		printf("Error: TransLoc_NeighborDistance, probeset number not match!\n");
		exit(EXIT_FAILURE);
	}


	/* compute correlation and distance */
	vCorr = NULL;
	vCorr = (struct DOUBLEMATRIX **)calloc(nGeneNum, sizeof(struct DOUBLEMATRIX*));
	if(vCorr == NULL)
	{
		printf("Error: TransLoc_NeighborDistance, cannot create space for computing!\n");
		exit(EXIT_FAILURE);
	}

	vDist = NULL;
	vDist = (struct DOUBLEMATRIX **)calloc(nGeneNum, sizeof(struct DOUBLEMATRIX*));
	if(vDist == NULL)
	{
		printf("Error: TransLoc_NeighborDistance, cannot create space for computing!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nGeneNum; ni++)
	{
		vCorr[ni] = NULL;
		vCorr[ni] = CreateDoubleMatrix(1,TRANSLOC_GENECORR_WIN);
		if(vCorr[ni] == NULL)
		{
			printf("Error: TransLoc_NeighborDistance, cannot create space for computing!\n");
			exit(EXIT_FAILURE);
		}

		vDist[ni] = NULL;
		vDist[ni] = CreateDoubleMatrix(1,TRANSLOC_GENECORR_WIN);
		if(vDist[ni] == NULL)
		{
			printf("Error: TransLoc_NeighborDistance, cannot create space for computing!\n");
			exit(EXIT_FAILURE);
		}
	}

	for(ni=0; ni<(nGeneNum-nGeneCorrHalfWin); ni++)
	{
		for(nj=1; nj<=nGeneCorrHalfWin; nj++)
		{
			/* compute correlation */
			nk = ni+nj;
			nS0 = pGeneP1->pMatElement[ni]-pGeneP0->pMatElement[ni]+1;
			nS1 = pGeneP1->pMatElement[nk]-pGeneP0->pMatElement[nk]+1;
			nST = (int)(nS0*nS1);
			nEffectCorrNum = 0;

			pTempCorr = NULL;
			pTempCorr = CreateDoubleMatrix(1, nST);
			if(pTempCorr == NULL)
			{
				printf("Error: TransLoc_NeighborDistance, cannot create space for computing!\n");
				exit(EXIT_FAILURE);
			}

			pTempDist = NULL;
			pTempDist = CreateDoubleMatrix(1, nST);
			if(pTempDist == NULL)
			{
				printf("Error: TransLoc_NeighborDistance, cannot create space for computing!\n");
				exit(EXIT_FAILURE);
			}

			/* individual pair */
			nz = 0;
			for(nx=pGeneP0->pMatElement[ni]; nx<=pGeneP1->pMatElement[ni]; nx++)
			{
				for(ny=pGeneP0->pMatElement[nk]; ny<=pGeneP1->pMatElement[nk]; ny++)
				{
					dTemp = TransLoc_ComputeCorrelation(vData[nx],vData[ny]);
					pTempCorr->pMatElement[nz] = 1.0-fabs(dTemp);

					dTemp = (pStart->pMatElement[ny]+pEnd->pMatElement[ny])/2-(pStart->pMatElement[nx]+pEnd->pMatElement[nx])/2;
					if(pChr->pMatElement[ny] != pChr->pMatElement[nx])
						dTemp = TRANSLOC_DIST_UPPERBOUND+1.0;

					pTempDist->pMatElement[nz] = dTemp;
					if(pTempDist->pMatElement[nz] <= TRANSLOC_DIST_UPPERBOUND)
						nEffectCorrNum++;

					nz++;
				}
			}

			/* keep correlations for probes within the distance limit */
			if(nEffectCorrNum > 0)
			{
				pFCorr = NULL;
				pFCorr = CreateDoubleMatrix(1, nEffectCorrNum);
				if(pFCorr == NULL)
				{
					printf("Error: TransLoc_NeighborDistance, cannot create space for computing!\n");
					exit(EXIT_FAILURE);
				}

				pFDist = NULL;
				pFDist = CreateDoubleMatrix(1, nEffectCorrNum);
				if(pFDist == NULL)
				{
					printf("Error: TransLoc_NeighborDistance, cannot create space for computing!\n");
					exit(EXIT_FAILURE);
				}

				ny = 0;
				for(nx=0; nx<nST; nx++)
				{
					if(pTempDist->pMatElement[nx] <= TRANSLOC_DIST_UPPERBOUND)
					{
						pFCorr->pMatElement[ny] = pTempCorr->pMatElement[nx];
						pFDist->pMatElement[ny] = pTempDist->pMatElement[nx];
						ny++;
					}
				}

				/* sort */
				pSortCorr = NULL;
				pSortIndex = NULL;
				DMSORTMERGEA_0(pFCorr, &pSortCorr, &pSortIndex);
				nz = nEffectCorrNum/2;
				dCorr = pSortCorr->pMatElement[nz];
				ny = pSortIndex->pMatElement[nz];
				dDist = pFDist->pMatElement[ny];
			}
			else
			{
				/* todo */
				dCorr = 1.0;
				dDist = TRANSLOC_DIST_UPPERBOUND+1.0;
			}

			/* release memory */
			DestroyDoubleMatrix(pTempCorr);
			DestroyDoubleMatrix(pTempDist);
			DestroyDoubleMatrix(pFCorr);
			DestroyDoubleMatrix(pFDist);
			DestroyDoubleMatrix(pSortCorr);
			DestroyLongMatrix(pSortIndex);
		}
	}

	/* save */

	/* release memory */
	DestroyIntMatrix(pGeneId);
	DestroyIntMatrix(pGeneP0);
	DestroyIntMatrix(pGeneP1);
	DestroyIntMatrix(pChr);
	DestroyIntMatrix(pStart);
	DestroyIntMatrix(pEnd);
	DestroyByteMatrix(pStrand);
	for(ni=0; ni<nProbeNum; ni++)
	{
		DeleteString(vProbe[ni]);
		vProbe[ni] = NULL;

		DestroyDoubleMatrix(vData[ni]);
		vData[ni] = NULL;
	}
	free(vProbe);
	free(vData);

	for(ni=0; ni<nGeneNum; ni++)
	{
		DestroyDoubleMatrix(vCorr[ni]);
		vCorr[ni] = NULL;
		DestroyDoubleMatrix(vDist[ni]);
		vDist[ni] = NULL;
	}
	free(vCorr);
	free(vDist);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TransLoc_ComputeCorrelation()                                          */
/*  compute correlation of two expression vectors.                         */
/* ----------------------------------------------------------------------- */ 
double TransLoc_ComputeCorrelation(struct DOUBLEMATRIX *pVec1, struct DOUBLEMATRIX *pVec2)
{
	/* define */
	double dDist = 0.0;
	int ni,nj;
	double dM1,dM2,dS1,dS2,dCov;
	double dTemp1,dTemp2;

	/* init check */
	if(pVec1->nWidth != pVec2->nWidth)
	{
		printf("Error: TransLoc_ComputeCorrelation, vector dimension not match!\n");
		exit(EXIT_FAILURE);
	}

	/* compute */
	dM1 = 0.0;
	for(ni=0; ni<pVec1->nWidth; ni++)
	{
		dM1 += pVec1->pMatElement[ni];
	}
	dM1 /= (double)(pVec1->nWidth);

	dM2 = 0.0;
	for(ni=0; ni<pVec2->nWidth; ni++)
	{
		dM2 += pVec2->pMatElement[ni];
	}
	dM2 /= (double)(pVec2->nWidth);

	dS1 = 0.0;
	dS2 = 0.0;
	dCov = 0.0;
	for(ni=0; ni<pVec1->nWidth; ni++)
	{
		dTemp1 = pVec1->pMatElement[ni]-dM1;
		dS1 += dTemp1*dTemp1;
		dTemp2 = pVec2->pMatElement[ni]-dM2;
		dS2 += dTemp2*dTemp2;
		dCov += dTemp1*dTemp2;
	}

	if( (fabs(dS1) < ZERO_BOUND) || (fabs(dS2) < ZERO_BOUND) )
	{
		dDist = 0.0;
	}
	else
	{
		dDist = dCov/sqrt(dS1*dS2);
	}

	/* return */
	return dDist;
}