/* ----------------------------------------------------------------------- */
/*  TilingArrayLib.c : implementation of the tiling array library          */
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
#include "MotifLib.h"
#include "AffyLib.h"
#include "SequenceLib.h"
#include "TilingArrayLib.h"

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ProbeReMapping_Main()                                        */
/*  Remap probes to a new genome version.                                  */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_ProbeReMapping_Main(char strBpmapFile[], char strGenomePath[],
						char strOutputFile[], char strSpecies[],
						int nKeyLen)
{
	/* define */
	char strGenomePathC[MED_LINE_LENGTH];
	char strFileName[MED_LINE_LENGTH];
	struct INTMATRIX *pChrLen = NULL;
	int nChrNum;
	struct tagString **vChrName = NULL; 
	int ni,nj;
	FILE *fpIn;
	char strLine[MED_LINE_LENGTH];

	/* mapping list */
	struct tagAffyBpMapUnit ***vvProbeMap = NULL;
	struct INTMATRIX *pProbeNum = NULL;
	struct tagAffyBpMapUnit *pProbeList = NULL;
	
	/* hash table for genome sequences */
	int *vGenomeHashIndex = NULL;
	int *vGenomeHashTable = NULL;
	int nBaseTypeNum = 4;
	int nDim,nNumRead;
	int nTotalWord = 0;

	/* step0: initialize */
	nDim = (int)pow((double)nBaseTypeNum, (double)nKeyLen);
	strcpy(strGenomePathC, strGenomePath);
	if(nKeyLen <= 0)
	{
		printf("Error: TileMapv2_ProbeReMapping_Main, key length should be >0!\n");
		exit(EXIT_FAILURE);
	}

	AdjustDirectoryPath(strGenomePathC);
	sprintf(strFileName, "%schrlen.txt", strGenomePathC);
	pChrLen = IMLOAD(strFileName);
	if(pChrLen == NULL)
	{
		printf("Error: TileMapv2_ProbeReMapping_Main, cannot load chromosome length!\n");
		exit(EXIT_FAILURE);
	}
	nChrNum = pChrLen->nHeight;
	if(nChrNum <= 0)
	{
		printf("Error: TileMapv2_ProbeReMapping_Main, cannot load chromosome length!\n");
		exit(EXIT_FAILURE);
	}

	vChrName = NULL;
	vChrName = (struct tagString **)calloc(nChrNum, sizeof(struct tagString *));
	if(vChrName == NULL)
	{
		printf("Error: TileMapv2_ProbeReMapping_Main, cannot load chromosome name!\n");
		exit(EXIT_FAILURE);
	}

	sprintf(strFileName, "%schrlist.txt", strGenomePathC);
	fpIn = NULL;
	fpIn = fopen(strFileName, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileMapv2_ProbeReMapping_Main, cannot load chromosome name!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, MED_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(ni >= nChrNum)
		{
			printf("Error: TileMapv2_ProbeReMapping_Main, chromosome number not match!\n");
			exit(EXIT_FAILURE);
		}

		StringAddTail(vChrName+ni, strLine);
		ni++;
	}

	fclose(fpIn);

	if(ni != nChrNum)
	{
		printf("Error: TileMapv2_ProbeReMapping_Main, chromosome number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* step1: create a mapping list */
	vvProbeMap = (struct tagAffyBpMapUnit ***)calloc(nChrNum, sizeof(struct tagAffyBpMapUnit **));
	if(vvProbeMap == NULL)
	{
		printf("Error: TileMapv2_ProbeReMapping_Main, cannot create map!\n");
		exit(EXIT_FAILURE);
	}
	pProbeNum = CreateIntMatrix(nChrNum, 1);
	if(pProbeNum == NULL)
	{
		printf("Error: TileMapv2_ProbeReMapping_Main, cannot create vector for recording probe number!\n");
		exit(EXIT_FAILURE);
	}

	/* step2: process chromomsome by chromosome */
	for(ni=0; ni<nChrNum; ni++)
	{
		printf("Processing %s ...\n", vChrName[ni]->m_pString);

		/* step3: for each chromosome, create a hash table */
		vGenomeHashIndex = NULL;
		vGenomeHashIndex = (int *)calloc((nDim+1), sizeof(int));
		if(vGenomeHashIndex == NULL)
		{
			printf("Error: TileMapv2_ProbeReMapping_Main, cannot create hash index!\n");
			exit(EXIT_FAILURE);
		}

		sprintf(strFileName, "%s%s.hashidx", strGenomePathC, vChrName[ni]->m_pString);
		fpIn = NULL;
		fpIn = fopen(strFileName, "rb");
		if(fpIn == NULL)
		{
			printf("Error: TileMapv2_ProbeReMapping_Main, cannot load hash index!\n");
			exit(EXIT_FAILURE);
		}
		nNumRead = fread(vGenomeHashIndex, sizeof(int), nDim+1, fpIn);
		fclose(fpIn);
		if( nNumRead != (nDim+1) )
		{
			printf("Error: TileMapv2_ProbeReMapping_Main, cannot load hash index correctly!\n");
			exit(EXIT_FAILURE);
		}


		nTotalWord = vGenomeHashIndex[nDim];
		vGenomeHashTable = NULL;
		vGenomeHashTable = (int *)calloc(nTotalWord, sizeof(int));
		if(vGenomeHashTable == NULL)
		{
			printf("Error: TileMapv2_ProbeReMapping_Main, cannot create hash map!\n");
			exit(EXIT_FAILURE);
		}

		sprintf(strFileName, "%s%s.hashmap", strGenomePathC, vChrName[ni]->m_pString);
		fpIn = NULL;
		fpIn = fopen(strFileName, "rb");
		if(fpIn == NULL)
		{
			printf("Error: TileMapv2_ProbeReMapping_Main, cannot load hash table!\n");
			exit(EXIT_FAILURE);
		}
		nNumRead = fread(vGenomeHashTable, sizeof(int), nTotalWord, fpIn);
		fclose(fpIn);
		if( nNumRead != nTotalWord )
		{
			printf("Error: TileMapv2_ProbeReMapping_Main, cannot load hash map correctly!\n");
			exit(EXIT_FAILURE);
		}


		/* step4: load the bpmap file; for each probe, search the genomic match */
		sprintf(strFileName, "%s%s.sq", strGenomePathC, vChrName[ni]->m_pString);
		pProbeList = TileMapv2_ProbeReMapping_ProbeAlign(strBpmapFile, strSpecies, 
			strFileName, pChrLen->pMatElement[ni], nKeyLen, nBaseTypeNum,
			vGenomeHashIndex, nDim,
			vGenomeHashTable, nTotalWord, pProbeNum->pMatElement+ni);

		/* step5: resort the mapping list */
		vvProbeMap[ni] = TileMapv2_ProbeReMapping_SortProbeList(pProbeNum->pMatElement[ni], &pProbeList);
		
		/* step 6: destroy the hash table */
		free(vGenomeHashIndex);
		free(vGenomeHashTable);
	}

	/* step7: rewrite the bpmap file */
	TileMapv2_ProbeReMapping_WriteBpmapFile(strOutputFile, 
			strSpecies, nChrNum, vChrName, vvProbeMap, pProbeNum);

	/* release memory */
	for(ni=0; ni<nChrNum; ni++)
	{
		for(nj=0; nj<pProbeNum->pMatElement[ni]; nj++)
		{
			AffyBpMapUnitDestroy(*(vvProbeMap[ni]+nj));
		}
		free(vvProbeMap[ni]);
		vvProbeMap[ni] = NULL;

		DeleteString(vChrName[ni]);
		vChrName[ni] = NULL;
	}
	free(vvProbeMap);
	free(vChrName);

	DestroyIntMatrix(pChrLen);
	DestroyIntMatrix(pProbeNum);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ProbeReMapping_ProbeAlign()                                  */
/*  Align probes to genome.                                                */
/* ----------------------------------------------------------------------- */ 
struct tagAffyBpMapUnit *TileMapv2_ProbeReMapping_ProbeAlign(char strBpmapFile[], 
				char strSpecies[], char strSeqFile[], int nChrLen, 
				int nKeyLen, int nBaseTypeNum,
				int *vGenomeHashIndex, int nIndexNum,
				int *vGenomeHashTable, int nTableNum,
				int *pProbeNum)
{
	/* define */
	int pnTotalProbeNum;
	FILE *fpIn;

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
	struct tagAffyBpMapUnit *pUnitList,*pMList;

	/* count */
	int ni,nj;

	/* load */
	pUnitList = NULL;
	pCUnit = NULL;

	fpIn = NULL;
	fpIn = fopen(strBpmapFile, "rb");
	if(fpIn == NULL)
	{
		printf("Error: TileMapv2_ProbeReMapping_ProbeAlign, cannot open .bpmap file!\n");
		exit(EXIT_FAILURE);
	}
	
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

	pnTotalProbeNum = 0;
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

			/* map probe sequence */
			pMList = TileMapv2_ProbeReMapping_HashAlign(pNewUnit,
				strSeqFile, nChrLen, nKeyLen, nBaseTypeNum, 
				vGenomeHashIndex, nIndexNum,
				vGenomeHashTable, nTableNum, 1);
			
			while(pMList != NULL)
			{
				pPUnit = pMList;
				pMList = pPUnit->pNext;
				pPUnit->pNext = NULL;

				if(pUnitList == NULL)
				{
					pUnitList = pPUnit;
					pCUnit = pPUnit;
				}
				else
				{
					pCUnit->pNext = pPUnit;
					pCUnit = pPUnit;
				}

				pnTotalProbeNum += 1;
			}

			pMList = TileMapv2_ProbeReMapping_HashAlign(pNewUnit,
				strSeqFile, nChrLen, nKeyLen, nBaseTypeNum, 
				vGenomeHashIndex, nIndexNum,
				vGenomeHashTable, nTableNum, 0);
			while(pMList != NULL)
			{
				pPUnit = pMList;
				pMList = pPUnit->pNext;
				pPUnit->pNext = NULL;

				if(pUnitList == NULL)
				{
					pUnitList = pPUnit;
					pCUnit = pPUnit;
				}
				else
				{
					pCUnit->pNext = pPUnit;
					pCUnit = pPUnit;
				}

				pnTotalProbeNum += 1;
			}

			AffyBpMapUnitDestroy(pNewUnit);
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

	*pProbeNum = pnTotalProbeNum;

	/* return */
	return pUnitList;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ProbeReMapping_HashAlign()                                   */
/*  Align a probe to a chromosome.                                         */
/* ----------------------------------------------------------------------- */ 
struct tagAffyBpMapUnit *TileMapv2_ProbeReMapping_HashAlign(struct tagAffyBpMapUnit *pUnit,
			char strSeqFile[], int nChrLen, int nKeyLen, int nBaseTypeNum,
			int *vGenomeHashIndex, int nIndexNum,
			int *vGenomeHashTable, int nTableNum,
			int nForward)
{
	/* define */
	struct tagAffyBpMapUnit *pUnitMap = NULL;
	struct tagAffyBpMapUnit *pCUnit;
	struct tagAffyBpMapUnit *pNewUnit;
	char strSeq[TILE_PROBE_LEN];
	int ni,nz;
	int nWordID;
	int nBadNum;
	char chBase;
	int nLine1,nLine2;
	int nPos1,nPos2;
	struct tagSequence *pSeq;

	/* code sequence */
	if(2*nKeyLen-1 > pUnit->bProbeLen)
	{
		printf("Error: TileMapv2_ProbeReMapping_HashAlign, probe seq too short for hash search!\n");
		exit(EXIT_FAILURE);
	}

	if(nForward == 0)
	{
		for(ni=0; ni<pUnit->bProbeLen; ni++)
		{
			chBase = pUnit->strProbeSeq[pUnit->bProbeLen-1-ni];
			switch(chBase)
			{
				case 'A': strSeq[ni] = 'T';
					break;
				case 'C': strSeq[ni] = 'G';
					break;
				case 'G': strSeq[ni] = 'C';
					break;
				case 'T': strSeq[ni] = 'A';
					break;
				default: strSeq[ni] = 'N';
					break;
			}
		}
		strSeq[pUnit->bProbeLen]= '\0';
	}
	else
	{
		strcpy(strSeq, pUnit->strProbeSeq);
	}

	nWordID = 0;
	nBadNum = 0;
	nz = 0;
	for(ni=0; ni<nKeyLen; ni++)
	{
		switch(strSeq[nz])
		{
			case 'A': nWordID = nWordID*nBaseTypeNum;
					break;
			case 'C': nWordID = nWordID*nBaseTypeNum+1;
					break;
			case 'G': nWordID = nWordID*nBaseTypeNum+2;
					break;
			case 'T': nWordID = nWordID*nBaseTypeNum+3;
					break;
			default: nWordID = nWordID*nBaseTypeNum;
					nBadNum += 1;
					break;
		}
		nz = nz+2;
	}

	/* probe alignment */
	if(nBadNum > 0)
	{
		return NULL;
	}
	if(nWordID >= nIndexNum)
	{
		printf("Error: TileMapv2_ProbeReMapping_HashAlign, word index out of range!\n");
		exit(EXIT_FAILURE);
	}

	nLine1 = vGenomeHashIndex[nWordID];
	nLine2 = vGenomeHashIndex[nWordID+1];

	for(ni=nLine1; ni<nLine2; ni++)
	{
		if(ni >= nTableNum)
		{
			printf("Error: TileMapv2_ProbeReMapping_HashAlign, hash index/map inconsistent!\n");
			exit(EXIT_FAILURE);
		}

		nPos1 = vGenomeHashTable[ni];
		nPos2 = nPos1+pUnit->bProbeLen-1;
		if( (nPos1 >= nChrLen) || (nPos2 >= nChrLen) )
		{
			continue;
		}

		pSeq = NULL;
		pSeq = Genome_Code_4bit_GetSeq(strSeqFile, nPos1, nPos2);
		if(pSeq != NULL)
		{
			if(strcmp(strSeq, pSeq->m_pSequence->m_pString) == 0)
			{
				pNewUnit = NULL;
				pNewUnit = AffyBpMapUnitCreate();
				if(pNewUnit == NULL)
				{
					printf("Error: TileMapv2_ProbeReMapping_HashAlign, cannot create new probe/bpmap unit!\n");
					exit(EXIT_FAILURE);
				}

				pNewUnit->bProbeLen = pUnit->bProbeLen;
				pNewUnit->bStrand = (unsigned char)nForward;
				pNewUnit->fMatchScore = pUnit->fMatchScore;
				pNewUnit->nDepthNum = pUnit->nDepthNum;
				pNewUnit->nMMX = pUnit->nMMX;
				pNewUnit->nMMY = pUnit->nMMY;
				pNewUnit->nPMX = pUnit->nPMX;
				pNewUnit->nPMY = pUnit->nPMY;
				pNewUnit->nPos = nPos1;
				pNewUnit->nRepeatNum = pUnit->nRepeatNum;
				pNewUnit->pNext = NULL;
				strcpy(pNewUnit->strProbeSeq, pUnit->strProbeSeq);

				if(pUnitMap == NULL)
				{
					pUnitMap = pNewUnit;
					pCUnit = pNewUnit;
				}
				else
				{
					pCUnit->pNext = pNewUnit;
					pCUnit = pNewUnit;
				}
			}
			SequenceDelete(pSeq);
		}
	}

	/* return */
	return pUnitMap;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ProbeReMapping_SortProbeList()                               */
/*  Resorting probes according to their position.                          */
/* ----------------------------------------------------------------------- */ 
struct tagAffyBpMapUnit **TileMapv2_ProbeReMapping_SortProbeList(int nProbeNum, 
		struct tagAffyBpMapUnit **ppProbeList)
{
	/* define */
	int ni,nj;
	struct tagAffyBpMapUnit *pProbe,*pPrevProbe;
	struct tagAffyBpMapUnit **vProbeList;
	struct tagAffyBpMapUnit **vSortProbeList;
	struct INTMATRIX *pPos;
	struct INTMATRIX *pSortedPos = NULL;
	struct LONGMATRIX *pSortIndex = NULL;
	
	/* init */
	if(ppProbeList == NULL)
		return NULL;

	/* preparation for sorting positions */
	vProbeList = NULL;
	vProbeList = (struct tagAffyBpMapUnit **)calloc(nProbeNum, sizeof(struct tagAffyBpMapUnit *));
	if(vProbeList == NULL)
	{
		printf("Error: TileMapv2_ProbeReMapping_SortProbeList, cannot create memory for resorting probes!\n");
		exit(EXIT_FAILURE);
	}

	pPos = NULL;
	pPos = CreateIntMatrix(1, nProbeNum);
	if(pPos == NULL)
	{
		printf("Error: TileMapv2_ProbeReMapping_SortProbeList, cannot create memory for resorting probes!\n");
		exit(EXIT_FAILURE);
	}

	/* convert positions to vector */
	ni = 0;
	while(*ppProbeList != NULL)
	{
		if(ni >= nProbeNum)
		{
			printf("Error: TileMapv2_ProbeReMapping_SortProbeList, probe number not match!\n");
			exit(EXIT_FAILURE);
		}

		pProbe = *ppProbeList;
		*ppProbeList = pProbe->pNext;
		pProbe->pNext = NULL;
		vProbeList[ni] = pProbe;
		pPos->pMatElement[ni] = pProbe->nPos;
		ni++;
	}

	if(ni != nProbeNum)
	{
		printf("Error: TileMapv2_ProbeReMapping_SortProbeList, probe number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* sort */
	IMSORTMERGEA_0(pPos, &pSortedPos, &pSortIndex);

	vSortProbeList = NULL;
	vSortProbeList = (struct tagAffyBpMapUnit **)calloc(nProbeNum, sizeof(struct tagAffyBpMapUnit *));
	if(vSortProbeList == NULL)
	{
		printf("Error: TileMapv2_ProbeReMapping_SortProbeList, cannot create memory for resorting probes!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nProbeNum; ni++)
	{
		nj = pSortIndex->pMatElement[ni];
		vSortProbeList[ni] = vProbeList[nj];
		vProbeList[nj] = NULL;
	}

	/* release memory */
	free(vProbeList);
	DestroyIntMatrix(pPos);
	DestroyIntMatrix(pSortedPos);
	DestroyLongMatrix(pSortIndex);

	/* return */
	return vSortProbeList;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ProbeReMapping_WriteBpmapFile()                              */
/*  Write remapped probe to a new bpmap file.                              */
/* ----------------------------------------------------------------------- */
int TileMapv2_ProbeReMapping_WriteBpmapFile(char strOutputFile[], 
			char strSpecies[], int nChrNum, struct tagString **vChrName, 
			struct tagAffyBpMapUnit ***vvProbeMap, struct INTMATRIX *pProbeNum)
{
	/* define */
	FILE *fpOut;
	int ni,nj;
	struct tagAffyBpMapUnit *pUnit;

	/* open file */
	fpOut = NULL;
	fpOut = fopen(strOutputFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: TileMapv2_ProbeReMapping_WriteBpmapFile, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	/* write file */
	for(ni=0; ni<nChrNum; ni++)
	{
		for(nj=0; nj<pProbeNum->pMatElement[ni]; nj++)
		{
			pUnit = *(vvProbeMap[ni]+nj);
			fprintf(fpOut, "%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\n",
				vChrName[ni]->m_pString, pUnit->nPos, pUnit->bStrand,
				pUnit->strProbeSeq, pUnit->nPMX, pUnit->nPMY, pUnit->nMMX, pUnit->nMMY);
		}
	}

	/* close file */
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  ProbeGenomeInfoCreate()                                                */
/*  Create Probe Genome Info object.                                       */
/* ----------------------------------------------------------------------- */ 
struct tagProbeGenomeInfo *ProbeGenomeInfoCreate()
{
	/* define */
	struct tagProbeGenomeInfo *pInfo = NULL;

	/* create */
	pInfo = (struct tagProbeGenomeInfo *)calloc(1, sizeof(struct tagProbeGenomeInfo));
	if(pInfo == NULL)
	{
		printf("Error: ProbeGenomeInfoCreate, cannot create info object!\n");
		exit(EXIT_FAILURE);
	}

	pInfo->nChr = -1;
	pInfo->nPos = -1;
	pInfo->nX = -1;
	pInfo->nY = -1;
	pInfo->pProbe = NULL;
	pInfo->pNext = NULL;

	/* return */
	return pInfo;
}

/* ----------------------------------------------------------------------- */ 
/*  ProbeGenomeInfoDestroy()                                               */
/*  Destroy Probe Genome Info object.                                      */
/* ----------------------------------------------------------------------- */ 
void ProbeGenomeInfoDestroy(struct tagProbeGenomeInfo **pInfo)
{
	if(pInfo != NULL)
	{
		if(*pInfo != NULL)
		{
			(*pInfo)->pNext = NULL;
			DeleteString((*pInfo)->pProbe);
			(*pInfo)->pProbe = NULL;

			free(*pInfo);
			*pInfo = NULL;
		}
	}
}

/* ----------------------------------------------------------------------- */ 
/*  ChrPosCreate()                                                         */
/*  Create ChrPos object.                                                  */
/* ----------------------------------------------------------------------- */ 
struct tagChrPos *ChrPosCreate()
{
	/* define */
	struct tagChrPos *pChrPos = NULL;

	/* create */
	pChrPos = (struct tagChrPos *)calloc(1, sizeof(struct tagChrPos));
	if(pChrPos == NULL)
	{
		printf("Error: ChrPosCreate, cannot create ChrPos object!\n");
		exit(EXIT_FAILURE);
	}

	pChrPos->nPos = -1;
	pChrPos->pNext = NULL;

	/* return */
	return pChrPos;
}

/* ----------------------------------------------------------------------- */ 
/*  ChrPosDestroy()                                                        */
/*  Destroy ChrPos object.                                                 */
/* ----------------------------------------------------------------------- */ 
void ChrPosDestroy(struct tagChrPos **pChrPos)
{
	if(pChrPos != NULL)
	{
		if(*pChrPos != NULL)
		{
			free(*pChrPos);
			*pChrPos = NULL;
		}
	}
}

/* ----------------------------------------------------------------------- */ 
/*  ChrPosListCreate()                                                     */
/*  Create ChrPos List object.                                             */
/* ----------------------------------------------------------------------- */ 
struct tagChrPosList *ChrPosListCreate()
{
	/* define */
	struct tagChrPosList *pChrPosList = NULL;

	/* create */
	pChrPosList = (struct tagChrPosList *)calloc(1, sizeof(struct tagChrPosList));
	if(pChrPosList == NULL)
	{
		printf("Error: ChrPosListCreate, cannot create ChrPosList object!\n");
		exit(EXIT_FAILURE);
	}

	pChrPosList->nNodeNum = 0;
	pChrPosList->pPosList = NULL;

	/* return */
	return pChrPosList;
}

/* ----------------------------------------------------------------------- */ 
/*  ChrPosListDestroy()                                                    */
/*  Destroy ChrPosList object.                                             */
/* ----------------------------------------------------------------------- */ 
void ChrPosListDestroy(struct tagChrPosList **pChrPosList)
{
	/* define */
	struct tagChrPos *pChrPos;

	if(pChrPosList != NULL)
	{
		if(*pChrPosList != NULL)
		{
			while( (*pChrPosList)->pPosList != NULL)
			{
				pChrPos = (*pChrPosList)->pPosList;
				(*pChrPosList)->pPosList = pChrPos->pNext;
				pChrPos->pNext = NULL;
				(*pChrPosList)->nNodeNum -= 1;

				ChrPosDestroy(&pChrPos);
			}

			free(*pChrPosList);
			*pChrPosList = NULL;
		}
	}
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ImportAffy_Main()                                            */
/*  Create a tiling array analysis project by loading affymetrix CEL data  */
/*  The data will be normalized, mapped to chromosomes, and saved to *.bar */
/*  files. If specified by users, intensities will also be exported to a   */
/*  combined .txt file.                                                    */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_ImportAffy_Main(char strParamPath[])
{
	/* define */

	/* working path */
	char strProjectTitle[MED_LINE_LENGTH];
	char strCELPath[MED_LINE_LENGTH];
	char strBpmapPath[MED_LINE_LENGTH];
	char strWorkPath[MED_LINE_LENGTH];
	
	/* arrays */
	int nLibNum = 0;
	struct tagString **vLibName = NULL;
	int nSampleNum = 0;
	int nArrayNum = 0;
	struct tagString **vCELPath;
	struct tagString **vArrayPath;
	struct tagString **vAlias;
	
	/* CEL files */
	int nTotalProbeNum = 0;
	int nRealProbeNum = 0;
	int nCELNumberCells = 0;
	int nCELTotalX = 0;
	int nCELTotalY = 0;
	struct tagCELData *pCELData;
	struct DOUBLEMATRIX *pArray;
	struct LONGMATRIX *pSortId;
	struct BYTEMATRIX *pMask;
	struct INTMATRIX *pNumMaskCells;
	struct DOUBLEMATRIX *pSortMean;
	struct DOUBLEMATRIX *pPercentile;
	int nTotalX, nTotalY, nArrayX, nArrayY;
	
	/* mask */
	int nRemoveMaskedCells = 0;
	int nRemoveOutlierCells = 0;

	/* normalization */
	int nIncludeNormalization = 1;
	double dNormLowerBound = 0.0;
	int nNormLogTransform = 0;
	double dLog2 = log(2.0);
	int nMissingNum;
	double dMissingValue = -1.0;
	
	/* intensity computation */
	struct tagBARData *pBARPos;
	int nIntensityType = 0;
	double dIntLowerBound = 1.0;
	int nIntLogTransform = 1;
	int nExportMode = 0;

	/* others */
	int ni,nj;
	char strLine[LONG_LINE_LENGTH];
	char strMaskPath[MED_LINE_LENGTH];
	char strPrcPath[MED_LINE_LENGTH];

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
		printf("Error: TileMapv2_ImportAffy_Main, cannot open the parameter file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strstr(strLine, "[Project title]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_ImportAffy_Main, cannot load project title!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strProjectTitle, chSep);
		}
		else if(strstr(strLine, "[CEL directory]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_ImportAffy_Main, cannot load CEL directory!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strCELPath, chSep);
			AdjustDirectoryPath(strCELPath);
		}
		else if(strstr(strLine, "[BPMAP directory]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_ImportAffy_Main, cannot load BPMAP directory!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strBpmapPath, chSep);
			AdjustDirectoryPath(strBpmapPath);
		}
		else if(strstr(strLine, "[Working directory]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_ImportAffy_Main, cannot load working directory!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strWorkPath, chSep);
			AdjustDirectoryPath(strWorkPath);
		}
		
		else if(strstr(strLine, "[No. of Libraries]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_ImportAffy_Main, cannot load no. of libraries!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			nLibNum = atoi(chSep);
			if(nLibNum <= 0)
			{
				printf("Warning: No BPMAP libraries provided!");
				return PROC_SUCCESS;
			}

			vLibName = NULL;
			vLibName = (struct tagString **)calloc(nLibNum, sizeof(struct tagString *));
			if(vLibName == NULL)
			{
				printf("Error: TileMapv2_ImportAffy_Main, cannot allocate memory for loading BPMAP lists!\n");
				exit(EXIT_FAILURE);
			}
		}
		else if(strstr(strLine, "[Libraries]") == strLine)
		{
			ni = 0;
			while(ni < nLibNum)
			{
				fgets(strLine, LONG_LINE_LENGTH, fpIn);
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] == '\0')
					continue;

				StringAddTail(vLibName+ni, strLine);
				ni++;
			}
		}
		else if(strstr(strLine, "[No. of Samples]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_ImportAffy_Main, cannot load no. of samples!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			nSampleNum = atoi(chSep);
            if(nSampleNum <= 0)
			{
				printf("Error: TileMapv2_ImportAffy_Main, no arrays available!\n");
				exit(EXIT_FAILURE);
			}

			nArrayNum = (int)(nSampleNum*nLibNum);

			vCELPath = NULL;
			vCELPath = (struct tagString **)calloc(nArrayNum, sizeof(struct tagString *));
			if(vCELPath == NULL)
			{
				printf("Error: TileMapv2_ImportAffy_Main, cannot allocate memory for tracking *.CEL files!\n");
				exit(EXIT_FAILURE);
			}

			vArrayPath = NULL;
			vArrayPath = (struct tagString **)calloc(nArrayNum, sizeof(struct tagString *));
			if(vArrayPath == NULL)
			{
				printf("Error: TileMapv2_ImportAffy_Main, cannot allocate memory for tracking array files!\n");
				exit(EXIT_FAILURE);
			}

			vAlias = NULL;
			vAlias = (struct tagString **)calloc(nSampleNum, sizeof(struct tagString *));
			if(vAlias == NULL)
			{
				printf("Error: TileMapv2_ImportAffy_Main, cannot allocate memory for tracking *.CEL files!\n");
				exit(EXIT_FAILURE);
			}

			pNumMaskCells = NULL;
			pNumMaskCells = CreateIntMatrix(1,nArrayNum);
			if(pNumMaskCells == NULL)
			{
				printf("Error: TileMapv2_ImportAffy_Main, cannot allocate memory for tracking number of masked cells!\n");
				exit(EXIT_FAILURE);
			}
		}
		
		else if(strstr(strLine, "[Arrays]") == strLine)
		{
			ni = 0;
			while(ni < nSampleNum)
			{
				fgets(strLine, LONG_LINE_LENGTH, fpIn);
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] == '\0')
					continue;
				
				if(strLine[0] != '>')
				{
					printf("Error: TileMapv2_ImportAffy_Main, error when loading samples!\n");
					exit(EXIT_FAILURE);
				}

				chSep = strLine+1;
				StrTrimLeft(chSep);
				StringAddTail(vAlias+ni, chSep);
				
				nj = 0;
				while(nj < nLibNum)
				{
					fgets(strLine, LONG_LINE_LENGTH, fpIn);
					StrTrimLeft(strLine);
					StrTrimRight(strLine);
					if(strLine[0] == '\0')
						continue;

					StringAddTail(vCELPath+(ni*nLibNum)+nj, strLine);
					chSep = NULL;
					chSep = strrchr(strLine, '.');
					if(chSep != NULL)
						*chSep = '\0';
					StringAddTail(vArrayPath+(ni*nLibNum)+nj, strLine);
					nj++;
				}
				ni++;
			}
		}

		else if(strstr(strLine, "[Remove masked cells in the CEL files]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_ImportAffy_Main, cannot load masked cell option!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			nRemoveMaskedCells = atoi(chSep);
		}

		else if(strstr(strLine, "[Remove outlier cells in the CEL files]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_ImportAffy_Main, cannot load masked cell option!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			nRemoveOutlierCells = atoi(chSep);
		}

		else if(strstr(strLine, "[Apply normalization before computing intensity]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_ImportAffy_Main, cannot load normalization option!\n");
				exit(EXIT_FAILURE);
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
				printf("Error: TileMapv2_ImportAffy_Main, cannot load normalization truncation lower bound!\n");
				exit(EXIT_FAILURE);
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
				printf("Error: TileMapv2_ImportAffy_Main, cannot load normalization log-transformation option!\n");
				exit(EXIT_FAILURE);
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
				printf("Error: TileMapv2_ImportAffy_Main, cannot load intensity option!\n");
				exit(EXIT_FAILURE);
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
				printf("Error: TileMapv2_ImportAffy_Main, cannot load intensity truncation option!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			dIntLowerBound = atof(chSep);
			dMissingValue = dIntLowerBound-1e6;
		}
		
		else if(strstr(strLine, "[Take log2 transformation after intensity computation]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_ImportAffy_Main, cannot load intensity log-transformation option!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			nIntLogTransform = atoi(chSep);
		}

		else if(strstr(strLine, "[Output mode]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_ImportAffy_Main, cannot load export option!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			nExportMode = atoi(chSep);
		}

		

		else if(strLine[0] == '#')
		{
		}

		else
		{
			printf("Error: TileMapv2_ImportAffy_Main, unknown parameters!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	fclose(fpIn);

	printf("########################################\n");
	printf("# Loading Raw Data                     #\n");
	printf("########################################\n");

	for(ni=0; ni<nArrayNum; ni++)
	{
		/* load values */
		pArray = NULL;
		nMissingNum = 0;
	
		printf("Loading %s ", vCELPath[ni]->m_pString);
		sprintf(strLine, "%s%s", strCELPath, vCELPath[ni]->m_pString);
		pCELData = TileMapv2_LoadCEL(strLine);
		if(pCELData == NULL)
		{
			printf("Error: TileMapv2_ImportAffy_Main, cannot load *.CEL file!\n");
			exit(EXIT_FAILURE);
		}

		/* check if the number of cells can be matched */
		if(ni == 0)
		{
			nCELNumberCells = pCELData->nNumberCells;
			nCELTotalX = pCELData->nCols;
			nCELTotalY = pCELData->nRows;
			if(nCELTotalX*nCELTotalY != nCELNumberCells)
			{
				printf("Error: TileMapv2_ImportAffy_Main, CEL file Cols*Rows != NumberCells!\n");
				exit(EXIT_FAILURE);
			}
		}
		else
		{
			if( (nCELNumberCells != pCELData->nNumberCells) || (nCELTotalX != pCELData->nCols)
				|| (nCELTotalY != pCELData->nRows) )
			{
				printf("Error: TileMapv2_ImportAffy_Main, CEL file dimensions do not match!\n");
				exit(EXIT_FAILURE);
			}
		}
		
		/* create and export mask */
		pMask = TileMapv2_GetMask(pCELData, nRemoveMaskedCells, nRemoveOutlierCells, 0);
		pNumMaskCells->pMatElement[ni] = (int)(BMSUM(pMask));
		
		if(nIncludeNormalization == 1)
		{
			/* load *.cel file and preprocess values */
			if(ni == 0)
			{
				pSortMean = NULL;
				pSortMean = CreateDoubleMatrix(1, pCELData->nNumberCells);
				if(pSortMean == NULL)
				{
					printf("Error: cannot allocate memory for quantile normalization!\n");
					exit(EXIT_FAILURE);
				}
			}

			pArray = NULL;
			pSortId = NULL;

			/* preprocessing values */
			if(nNormLogTransform == 1)
			{
				for(nj=0; nj<pCELData->nNumberCells; nj++)
				{
					/* truncate */
					if(pCELData->pIntensity->pMatElement[nj] < dNormLowerBound)
						pCELData->pIntensity->pMatElement[nj] = dNormLowerBound;

					/* log transformation */
					pCELData->pIntensity->pMatElement[nj] = log(pCELData->pIntensity->pMatElement[nj])/dLog2;
				}
			}
			else
			{
				for(nj=0; nj<pCELData->nNumberCells; nj++)
				{
					/* truncate */
					if(pCELData->pIntensity->pMatElement[nj] < dNormLowerBound)
						pCELData->pIntensity->pMatElement[nj] = dNormLowerBound;
				}
			}

			/* export sorted values */
			printf("  Sorting %s...\n", vCELPath[ni]->m_pString);
			DMSORTMERGEA_0(pCELData->pIntensity, &pArray, &pSortId);

			/* compute percentiles */
			pPercentile = TileMapv2_QuantileNorm_AddQuantile(pSortMean, pArray, pSortId, pMask, pNumMaskCells->pMatElement[ni]);
            			
			/* export raw values */
			sprintf(strLine, "%s%s.prc", strWorkPath, vCELPath[ni]->m_pString);
			TileMapv2_SaveToBinaryFile(pPercentile->pMatElement, sizeof(double), pCELData->nNumberCells, strLine);

			DestroyDoubleMatrix(pPercentile);
			DestroyDoubleMatrix(pArray);
			DestroyLongMatrix(pSortId);
		}
		else
		{
			/* export raw values */
			sprintf(strLine, "%s%s.int", strWorkPath, vCELPath[ni]->m_pString);
			TileMapv2_SaveToBinaryFile(pCELData->pIntensity->pMatElement, sizeof(double), pCELData->nNumberCells, strLine);
		}

		/* save masks */
		sprintf(strLine, "%s%s.mask", strWorkPath, vCELPath[ni]->m_pString);
		TileMapv2_SaveToBinaryFile(pMask->pMatElement, sizeof(unsigned char), pCELData->nNumberCells, strLine);
		DestroyByteMatrix(pMask);

		/* clear memory */
		Affy_CELData_Destroy(&pCELData);
	}

	printf("\n");

	/* quantile normalization */
	if(nIncludeNormalization == 1)
	{
		printf("########################################\n");
		printf("# Quantile Normalization               #\n");
		printf("########################################\n");

		if(nArrayNum > 0)
			DMPDIVTS(pSortMean, (double)nArrayNum);

		/* quantile normalization */
		for(ni=0; ni<nArrayNum; ni++)
		{
			printf("Scaling %s...\n", vCELPath[ni]->m_pString);
			sprintf(strPrcPath, "%s%s.prc", strWorkPath, vCELPath[ni]->m_pString);
			sprintf(strLine, "%s%s.int", strWorkPath, vCELPath[ni]->m_pString);
			TileMapv2_QuantileNorm_Rescale(strPrcPath, pSortMean, strLine);
			RemoveFiles(strPrcPath);
		}

		/* release memory */
		DestroyDoubleMatrix(pSortMean);
		printf("\n");
	}

	/* export intensities and prepare local repeat filters */
	printf("########################################\n");
	printf("# Exporting Intensities                #\n");
	printf("########################################\n");
	for(ni=0; ni<nLibNum; ni++)
	{
		printf("Loading %s...\n", vLibName[ni]->m_pString);
		sprintf(strLine, "%s%s", strBpmapPath, vLibName[ni]->m_pString);
		sprintf(strMaskPath, "%s%s.refmask", strWorkPath, vLibName[ni]->m_pString);
		pBARPos = NULL;
		pBARPos = TileMapv2_BpmapToBAR(strLine, strMaskPath);
		if(pBARPos == NULL)
		{
			printf("Error: TileMapv2_ImportAffy_Main, empty bpmap file!\n");
			exit(EXIT_FAILURE);
		}

		TileMapv2_ExportAffyIntensity(strWorkPath, strProjectTitle, vCELPath, 
			nCELNumberCells, nCELTotalX, nCELTotalY,
			nSampleNum, vAlias, nLibNum, ni, 
			pBARPos, nIntensityType, dIntLowerBound, nIntLogTransform, 1, nExportMode);
		
		Affy_BARData_Destroy(&pBARPos);
	}

	/* destroy */
	for(ni=0; ni<nArrayNum; ni++)
	{
		DeleteString(vCELPath[ni]);
		DeleteString(vArrayPath[ni]);
	}
	free(vCELPath);
	free(vArrayPath);

	for(ni=0; ni<nSampleNum; ni++)
	{
		DeleteString(vAlias[ni]);
	}
	free(vAlias);

	for(ni=0; ni<nLibNum; ni++)
	{
		DeleteString(vLibName[ni]);
	}
	free(vLibName);

	DestroyIntMatrix(pNumMaskCells);

	printf("Done!\n");

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_LoadCEL()                                                    */
/*  Load CEL files.                                                        */
/* ----------------------------------------------------------------------- */ 
struct tagCELData *TileMapv2_LoadCEL(char strFileName[])
{
	struct tagCELData *pCELData = NULL;

	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));
	
	FILE *fpCel;
	unsigned char bMagicnumber;
	int nMagicnumber = 0;
	int nVersionnumber = 0;
	int nIsV = 0;
	char *chSep;
	char strLine[MED_LINE_LENGTH];

	/* try generic file */
	fpCel = NULL;
	fpCel = fopen(strFileName, "rb");
	if(fpCel == NULL)
	{
		printf("Error: TileMapv2_LoadCEL, %s does not exist!\n", strFileName);
		exit(EXIT_FAILURE);
	}

	if(fread(&bMagicnumber, 1, 1, fpCel) != 1)
	{
		nIsV = 0;
	}
	else if(bMagicnumber != 59)
	{
		nIsV = 0;
	}
	else
	{
		nIsV = 5;
	}

	fclose(fpCel);

	if(nIsV == 5)
	{
		printf("(command console format)...\n");
        pCELData = Affy_LoadCEL_CmdCslv1(strFileName);
		return pCELData;
	}


	/* try version 4 */
	fpCel = NULL;
	fpCel = fopen(strFileName, "rb");
	if(fpCel == NULL)
	{
		printf("Error: TileMapv2_LoadCEL, %s does not exist!\n", strFileName);
		exit(EXIT_FAILURE);
	}

	if(little_endian_fread(&nMagicnumber, INT_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		nIsV = 0;
	}
	else if(nMagicnumber != 64)
	{
		nIsV = 0;
	}
	else
	{
		/* load version number */
		if(little_endian_fread(&nVersionnumber, INT_SIZE, 1, fpCel, little_endian_machine) != 1)
		{
			nIsV = 0;
		}
		else if(nVersionnumber != 4)
		{
			nIsV = 0;
		}
		else
		{
			nIsV = 4;
		}
	}

	fclose(fpCel);

	/* if CEL version = 4, load using Affy_LoadCELv4 */
	if(nIsV == 4)
	{
		printf("(v4)...\n");
        pCELData = Affy_LoadCELv4_Fast(strFileName);
		return pCELData;
	}

	/* try version 3 */
	fpCel = NULL;
	fpCel = fopen(strFileName, "r");
	if(fpCel == NULL)
	{
		printf("Error: TileMapv2_LoadCEL, %s does not exist!\n", strFileName);
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, MED_LINE_LENGTH, fpCel)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strstr(strLine, "[CEL]") != strLine)
		{
			nIsV = 0;
			break;
		}
		else
		{
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			chSep = strchr(strLine, '=');
			chSep++;
			if( atoi(chSep) != 3)
			{
				nIsV = 0;
			}
			else
			{
				nIsV = 3;
			}
			break;
		}
	}

	fclose(fpCel);

	/* if CEL version = 3, load using Affy_LoadCELv3 */
	if(nIsV == 3)
	{
		printf("(v3)...\n");
		pCELData = Affy_LoadCELv3(strFileName);
		return pCELData;
	}

	/* if not recognized, print error message */
	printf("Error: TileMapv2_LoadCEL, cannot recognize %s as a CEL file!\n", strFileName);
	exit(EXIT_FAILURE);

	return NULL;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_GetMask()                                                    */
/*  Get mask indicators.                                                   */
/* ----------------------------------------------------------------------- */ 
struct BYTEMATRIX *TileMapv2_GetMask(struct tagCELData *pCELData, 
	int nIncludeMasked, int nIncludeOutlier,  int nIncludeModified)
{
	/* define */
	struct BYTEMATRIX *pMask = NULL;
	int ni;
	int nidx;
	int nTotalX;

	/* create mask */
	if(pCELData == NULL)
		return NULL;

	pMask = CreateByteMatrix(1, pCELData->nNumberCells);
	if(pMask == NULL)
	{
		printf("Error: TileMapv2_GetMask, cannot allocate memory for mask indicators!\n");
		exit(EXIT_FAILURE);
	}

	nTotalX = pCELData->nTotalX;

	if(nIncludeMasked == 1)
	{
		for(ni=0; ni<(int)(pCELData->nMaskedCells); ni++)
		{
			nidx = pCELData->pMaskedY->pMatElement[ni]*nTotalX+pCELData->pMaskedX->pMatElement[ni];
			pMask->pMatElement[nidx] = 1;
		}
	}

	if(nIncludeOutlier == 1)
	{
		for(ni=0; ni<(int)(pCELData->nOutlierCells); ni++)
		{
			nidx = pCELData->pOutlierY->pMatElement[ni]*nTotalX+pCELData->pOutlierX->pMatElement[ni];
			pMask->pMatElement[nidx] = 1;
		}
	}

	if(nIncludeModified == 1)
	{
		for(ni=0; ni<(int)(pCELData->nModifiedCells); ni++)
		{
			nidx = pCELData->pModifiedY->pMatElement[ni]*nTotalX+pCELData->pModifiedX->pMatElement[ni];
			pMask->pMatElement[nidx] = 1;
		}
	}

	/* return */
	return pMask;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_SaveToBinaryFile()                                           */
/*  Write data to binary files.                                            */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_SaveToBinaryFile(const void *buffer, size_t size, size_t count, char strFileName[])
{
	/* define */
	FILE *fpOut;
	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));

	/* write */
	fpOut = NULL;
	fpOut = fopen(strFileName, "wb");
	if(fpOut == NULL)
	{
		printf("Error: TileMapv2_SaveToBinaryFile, cannot save data to binary files!\n");
		exit(EXIT_FAILURE);
	}

	if( little_endian_fwrite(buffer, size, count, fpOut, little_endian_machine) != count)
	{
		printf("Error: TileMapv2_SaveToBinaryFile, cannot save data to binary files correctly!\n");
		exit(EXIT_FAILURE);
	}

	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_LoadFromBinaryFile()                                         */
/*  Load data from binary files.                                           */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_LoadFromBinaryFile(void *buffer, size_t size, size_t count, char strFileName[])
{
	/* define */
	FILE *fpIn;
	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));

	/* load */
	fpIn = NULL;
	fpIn = fopen(strFileName, "rb");
	if(fpIn == NULL)
	{
		printf("Error: TileMapv2_LoadFromBinaryFile, cannot load data from binary files!\n");
		exit(EXIT_FAILURE);
	}

	if( little_endian_fread(buffer, size, count, fpIn, little_endian_machine) != count)
	{
		printf("Error: TileMapv2_LoadFromBinaryFile, cannot load data from binary files correctly!\n");
		exit(EXIT_FAILURE);
	}

	fclose(fpIn);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_QuantileNorm_AddQuantile()                                   */
/*  Add quantiles for quantile normalization.                              */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *TileMapv2_QuantileNorm_AddQuantile(struct DOUBLEMATRIX *pSortMean, 
			struct DOUBLEMATRIX *pArray, struct LONGMATRIX *pSortId, 
			struct BYTEMATRIX *pMask, int nMaskedCells)
{
	/* define */
	struct DOUBLEMATRIX *pPrc = NULL;
	int ni,nj,nk,nP1,nP2;
	long nId;
	int nPrcDenom;
	int nArrayDenom;
	double dR;
	double dR1,dR2,dR3;
	double dTemp,dLambda;

	/* init */
	if( (pSortMean == NULL) || (pArray == NULL) || (pSortId == NULL) )
	{
		return NULL;
	}

	/* compute percentiles */
	pPrc = CreateDoubleMatrix(pArray->nHeight, pArray->nWidth);
	if(pPrc == NULL)
	{
		printf("Error: TileMapv2_QuantileNorm_AddQuantile, cannot create percentile matrix!\n");
		exit(EXIT_FAILURE);
	}

	nPrcDenom = pPrc->nWidth;
	nArrayDenom = pArray->nWidth-nMaskedCells;
	if(nArrayDenom <= 0)
	{
		printf("Error: TileMapv2_QuantileNorm_AddQuantile, there is an array with all cells masked!\n");
		exit(EXIT_FAILURE);
	}

	nk = -1;
	nP1 = 0;
	nP2 = 0;
	dR2 = (double)nP2/(double)nArrayDenom;
	dR1 = dR2;
	nj = 1;

	for(ni=0; ni<pArray->nWidth; ni++)
	{
		nId = pSortId->pMatElement[ni];
		if( pMask->pMatElement[nId] == 1)
		{
			pPrc->pMatElement[nId] = dR2;
		}
		else
		{
			nP1 = nP2;
			nP2 += 1;
			dR1 = dR2;
			dR2 = (double)nP2/(double)nArrayDenom;
			pPrc->pMatElement[nId] = dR2;

			dR = (double)nP2*(double)nPrcDenom/(double)nArrayDenom;
			
			for(; nj<=dR; nj++)
			{
				if(nP1 == 0)
				{
					dTemp = pArray->pMatElement[ni];
				}
				else
				{
					dR3 = (double)nj/(double)nPrcDenom;
					dLambda = (dR3-dR1)/(dR2-dR1);
					dTemp = (1.0-dLambda)*pArray->pMatElement[nk]+dLambda*pArray->pMatElement[ni];
				}
				pSortMean->pMatElement[nj-1] += dTemp;
			}

			nk = ni;
		}
	}

	for(; nj<=nPrcDenom; nj++)
	{
		dTemp = pArray->pMatElement[nk];
		pSortMean->pMatElement[nj-1] += dTemp;
	}

	/* return */
	return pPrc;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_QuantileNorm_Rescale()                                       */
/*  Rescale arrays according to the quantile normalized values.            */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_QuantileNorm_Rescale(char strPrcPath[], struct DOUBLEMATRIX *pSortMean, 
								   char strOutPath[])
{
	/* define */
	struct DOUBLEMATRIX *pV = NULL;
	int ni;
	double dR;
	int nj,nk;

	/* init */
	if(pSortMean == NULL)
	{
		printf("Warning: TileMapv2_QuantileNorm_Rescale, empty quantiles!\n");
		return PROC_FAILURE;
	}

	/* load data */
	pV = CreateDoubleMatrix(pSortMean->nHeight, pSortMean->nWidth);
	if(pV == NULL)
	{
		printf("Error: TileMapv2_QuantileNorm_Rescale, cannot allocate memory for rescaling array intensities!\n");
		exit(EXIT_FAILURE);
	}

	TileMapv2_LoadFromBinaryFile(pV->pMatElement, sizeof(double), pSortMean->nWidth, strPrcPath);

	for(ni=0; ni<pV->nWidth; ni++)
	{
		dR = (pV->pMatElement[ni]*pSortMean->nWidth)-1.0;
		if(dR < 0.0)
			dR = 0.0;

		nj = (int)dR;
		nk = nj+1;
		if(nk >= pV->nWidth)
			nk = pV->nWidth-1;

		dR = dR-nj;
		pV->pMatElement[ni] = (1.0-dR)*pSortMean->pMatElement[nj]+dR*pSortMean->pMatElement[nk];
	}
	
	
	TileMapv2_SaveToBinaryFile(pV->pMatElement, sizeof(double), pV->nWidth, strOutPath);

	DestroyDoubleMatrix(pV);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_BpmapToBAR()                                                 */
/*  Create BAR template from bpmap file. The function will create local    */
/*  repeat mask at the same time.                                          */
/* ----------------------------------------------------------------------- */
struct tagBARData *TileMapv2_BpmapToBAR(char strBpmapFile[], char strMaskPath[])
{
	/* define */
	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));

	struct tagBARData *pBARPos = NULL;

	FILE *fpIn;
	FILE *fpOut;

	/* variables */
	char strFileType[9];
	float fVersion;
	unsigned int nSeqNum;

	/* seq info */
	unsigned int nSeqNameLen;

	struct INTMATRIX *vProbeMappingType;
	unsigned int nProbeMappingType;
	
	struct INTMATRIX *vSequenceFileOffset;
	unsigned int nSequenceFileOffset;
	
	struct INTMATRIX *vProbePairNum;
	unsigned int nProbePairNum;
	
	unsigned int nParamNum;
		
	/* map info */
	struct INTMATRIX *vSeqID;
	unsigned int nSeqID;
	int nProbeNum;
	int nTotalProbeNum = 0;
	int nMaskedProbeNum = 0;
	struct tagAffyBpMapUnit *pNewUnit,*pCUnit,*pPUnit;
	struct tagAffyBpMapUnit *pUnitList;
	struct DOUBLEMATRIX *pResizeMat;
	int nIgnore = 0;

	/* count */
	int ni,nj,nk;
	int nEndPos;

	/* load */
	fpIn = NULL;
	fpIn = fopen(strBpmapFile, "rb");
	if(fpIn == NULL)
	{
		printf("Error: TileMapv2_BpmapToBAR, cannot open .bpmap file!\n");
		exit(EXIT_FAILURE);
	}
	fpOut = NULL;
	fpOut = fopen(strMaskPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: TileMapv2_BpmapToBAR, cannot open .refmask file!\n");
		exit(EXIT_FAILURE);
	}
    fprintf(fpOut, "chromosome\tposition\tprobe_num\trepeat_num\tprobe_seq\tPMx\tPMy\tMMx\tMMy\n");
	
	/* load head */
	fread( strFileType, sizeof(char), 8, fpIn );
	strFileType[8] = '\0';
	AFFYBAR_READ_VERSION(fpIn, &fVersion);
	if(big_endian_fread(&nSeqNum, 4, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: TileMapv2_BpmapToBAR, cannot load sequence number.\n");
		exit(EXIT_FAILURE);
	}
	printf("  BPMAP Version = %f\n", fVersion);
	printf("  Number of Sequences = %d\n", nSeqNum);
	if(nSeqNum <= 0)
	{
		printf("Warning: empty sequences!\n");
		return NULL;
	}

	/* create BAR object */
	pBARPos = Affy_BARData_Create();
	if(pBARPos == NULL)
	{
		printf("Error: TileMapv2_BpmapToBAR, BARData object was not created!\n");
		exit(EXIT_FAILURE);
	}
	strcpy(pBARPos->strMagicnumber, "barr\r\n\032\n");
	pBARPos->fVersionnumber = 2.0;
    pBARPos->nSeqNum = nSeqNum;
	pBARPos->vSeqData = (struct tagBARSeq **)calloc(pBARPos->nSeqNum, sizeof(struct tagBARSeq *));
	if(pBARPos->vSeqData == NULL)
	{
		printf("Error: TileMapv2_BpmapToBAR, cannot allocate memory for loading sequence data.\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<pBARPos->nSeqNum; ni++)
	{
		/* create BARSeq object */
		pBARPos->vSeqData[ni] = Affy_BARSeq_Create();
		if(pBARPos->vSeqData[ni] == NULL)
		{
			printf("Error: TileMapv2_BpmapToBAR, cannot create BARSeq object.\n");
			exit(EXIT_FAILURE);
		}
	}

	vProbeMappingType = NULL;
	vProbeMappingType = CreateIntMatrix(nSeqNum,1);
	if(vProbeMappingType == NULL)
	{
		printf("Error: TileMapv2_BpmapToBAR, cannot load probe mapping type!\n");
		exit(EXIT_FAILURE);
	}

	if(fVersion > 2.5)
	{
		vSequenceFileOffset = NULL;
		vSequenceFileOffset = CreateIntMatrix(nSeqNum,1);
		if(vSequenceFileOffset == NULL)
		{
			printf("Error: TileMapv2_BpmapToBAR, cannot load sequence file offset!\n");
			exit(EXIT_FAILURE);
		}
	}

	vProbePairNum = NULL;
	vProbePairNum = CreateIntMatrix(nSeqNum,1);
	if(vProbePairNum == NULL)
	{
		printf("Error: TileMapv2_BpmapToBAR, cannot load probe pair number!\n");
		exit(EXIT_FAILURE);
	}
	
	for(ni=0; ni<(int)nSeqNum; ni++)
	{
		/* for version 1.0 or later */
		/* load sequence name */
		if(big_endian_fread(&nSeqNameLen, 4, 1, fpIn, little_endian_machine) != 1)
		{
			printf("Error: TileMapv2_BpmapToBAR, cannot load sequence name length.\n");
			exit(EXIT_FAILURE);
		}
		if(nSeqNameLen > 0)
		{
			pBARPos->vSeqData[ni]->pSeqName = CreateString(nSeqNameLen);
			if(pBARPos->vSeqData[ni]->pSeqName == NULL)
			{
				printf("Error: TileMapv2_BpmapToBAR, cannot load sequence name.\n");
				exit(EXIT_FAILURE);
			}
			if(big_endian_fread(pBARPos->vSeqData[ni]->pSeqName->m_pString, 1, nSeqNameLen, fpIn, little_endian_machine) != nSeqNameLen)
			{
				printf("Error: TileMapv2_BpmapToBAR, cannot load sequence name.\n");
				exit(EXIT_FAILURE);
			}
			pBARPos->vSeqData[ni]->pSeqName->m_pString[nSeqNameLen] = '\0';
		}

		/* for version 3.0 or later */
		/* load probe mapping type and sequence file offset */
		if(fVersion > 2.5)
		{
			if(big_endian_fread(&nProbeMappingType, 4, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: TileMapv2_BpmapToBAR, cannot load probe mapping type.\n");
				exit(EXIT_FAILURE);
			}
			vProbeMappingType->pMatElement[ni] = (int)nProbeMappingType;
			if(big_endian_fread(&nSequenceFileOffset, 4, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: TileMapv2_BpmapToBAR, cannot load sequence file offset.\n");
				exit(EXIT_FAILURE);
			}
			vSequenceFileOffset->pMatElement[ni] = (int)nSequenceFileOffset;
		}

		/* for version 1.0 or later */
		if(big_endian_fread(&nProbePairNum, 4, 1, fpIn, little_endian_machine) != 1)
		{
			printf("Error: TileMapv2_BpmapToBAR, cannot load probe pair number.\n");
			exit(EXIT_FAILURE);
		}
		vProbePairNum->pMatElement[ni] = (int)nProbePairNum;
		nTotalProbeNum += (int)nProbePairNum;

		/* for version 2.0 or later */
		if(fVersion > 1.5)
		{
			/* read group name */
			if(big_endian_fread(&nSeqNameLen, 4, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: TileMapv2_BpmapToBAR, cannot load group name length.\n");
				exit(EXIT_FAILURE);
			}
			if(nSeqNameLen > 0)
			{
				pBARPos->vSeqData[ni]->pSeqGroupName = CreateString(nSeqNameLen);
				if(pBARPos->vSeqData[ni]->pSeqGroupName == NULL)
				{
					printf("Error: TileMapv2_BpmapToBAR, cannot load group name.\n");
					exit(EXIT_FAILURE);
				}
				if(big_endian_fread(pBARPos->vSeqData[ni]->pSeqGroupName->m_pString, 1, nSeqNameLen, fpIn, little_endian_machine) != nSeqNameLen)
				{
					printf("Error: TileMapv2_BpmapToBAR, cannot load group name.\n");
					exit(EXIT_FAILURE);
				}
				pBARPos->vSeqData[ni]->pSeqGroupName->m_pString[nSeqNameLen] = '\0';
			}

			/* read version name */
			if(big_endian_fread(&nSeqNameLen, 4, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: TileMapv2_BpmapToBAR, cannot load version number length.\n");
				exit(EXIT_FAILURE);
			}
			if(nSeqNameLen > 0)
			{
				pBARPos->vSeqData[ni]->pSeqVersion = CreateString(nSeqNameLen);
				if(pBARPos->vSeqData[ni]->pSeqVersion == NULL)
				{
					printf("Error: TileMapv2_BpmapToBAR, cannot load version number.\n");
					exit(EXIT_FAILURE);
				}
				if(big_endian_fread(pBARPos->vSeqData[ni]->pSeqVersion->m_pString, 1, nSeqNameLen, fpIn, little_endian_machine) != nSeqNameLen)
				{
					printf("Error: TileMapv2_BpmapToBAR, cannot load version number.\n");
					exit(EXIT_FAILURE);
				}
				pBARPos->vSeqData[ni]->pSeqVersion->m_pString[nSeqNameLen] = '\0';
			}

			/* read paramters */
			if(big_endian_fread(&nParamNum, 4, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: TileMapv2_BpmapToBAR, cannot load number of parameters.\n");
				exit(EXIT_FAILURE);
			}
			pBARPos->vSeqData[ni]->nParamNum = (int)nParamNum;
			if(nParamNum > 0)
			{
				pBARPos->vSeqData[ni]->vParamName = (struct tagString **)calloc(pBARPos->vSeqData[ni]->nParamNum, sizeof(struct tagString *));
				pBARPos->vSeqData[ni]->vParamValue = (struct tagString **)calloc(pBARPos->vSeqData[ni]->nParamNum, sizeof(struct tagString *));

				if( (pBARPos->vSeqData[ni]->vParamName == NULL) || (pBARPos->vSeqData[ni]->vParamValue == NULL) )
				{
					printf("Error: TileMapv2_BpmapToBAR, cannot allocate memory for loading parameter name/value pairs\n");
					exit(EXIT_FAILURE);
				}
			}
			
			for(nj=0; nj<(int)nParamNum; nj++)
			{
				if(big_endian_fread(&nSeqNameLen, 4, 1, fpIn, little_endian_machine) != 1)
				{
					printf("Error: TileMapv2_BpmapToBAR, cannot load parameter name length.\n");
					exit(EXIT_FAILURE);
				}
				if(nSeqNameLen > 0)
				{
					pBARPos->vSeqData[ni]->vParamName[nj] = CreateString(nSeqNameLen);
					if(pBARPos->vSeqData[ni]->vParamName[nj] == NULL)
					{
						printf("Error: TileMapv2_BpmapToBAR, cannot load parameter name.\n");
						exit(EXIT_FAILURE);
					}
					if(big_endian_fread(pBARPos->vSeqData[ni]->vParamName[nj]->m_pString, 1, nSeqNameLen, fpIn, little_endian_machine) != nSeqNameLen)
					{
						printf("Error: TileMapv2_BpmapToBAR, cannot load parameter name.\n");
						exit(EXIT_FAILURE);
					}
					pBARPos->vSeqData[ni]->vParamName[nj]->m_pString[nSeqNameLen] = '\0';
				}

				if(big_endian_fread(&nSeqNameLen, 4, 1, fpIn, little_endian_machine) != 1)
				{
					printf("Error: TileMapv2_BpmapToBAR, cannot load parameter value length.\n");
					exit(EXIT_FAILURE);
				}
				if(nSeqNameLen > 0)
				{
					pBARPos->vSeqData[ni]->vParamValue[nj] = CreateString(nSeqNameLen);
					if(pBARPos->vSeqData[ni]->vParamValue[nj] == NULL)
					{
						printf("Error: TileMapv2_BpmapToBAR, cannot load parameter value.\n");
						exit(EXIT_FAILURE);
					}
					if(big_endian_fread(pBARPos->vSeqData[ni]->vParamValue[nj]->m_pString, 1, nSeqNameLen, fpIn, little_endian_machine) != nSeqNameLen)
					{
						printf("Error: TileMapv2_BpmapToBAR, cannot load parameter value.\n");
						exit(EXIT_FAILURE);
					}
					pBARPos->vSeqData[ni]->vParamValue[nj]->m_pString[nSeqNameLen] = '\0';
				}
			}
		}
	}

	/* load map */
	vSeqID = NULL;
	vSeqID = CreateIntMatrix(nSeqNum,1);
	if(vSeqID == NULL)
	{
		printf("Error: TileMapv2_BpmapToBAR, cannot load sequence id!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<(int)nSeqNum; ni++)
	{
		/* load seq id */
		if(big_endian_fread(&nSeqID, 4, 1, fpIn, little_endian_machine) != 1)
		{
			printf("Error: TileMapv2_BpmapToBAR, cannot load sequence ID.\n");
			exit(EXIT_FAILURE);
		}
		vSeqID->pMatElement[ni] = nSeqID;

		/* create initial memory */
		pBARPos->vSeqData[ni]->nDataNum = 0;
		if(vProbeMappingType->pMatElement[ni] == 1)
		{
			pBARPos->vSeqData[ni]->nColNum = 3;
		}
		else
		{
			pBARPos->vSeqData[ni]->nColNum = 5;
		}

		pBARPos->vSeqData[ni]->vData = (struct DOUBLEMATRIX **)calloc(pBARPos->vSeqData[ni]->nColNum, sizeof(struct DOUBLEMATRIX *));
		if(pBARPos->vSeqData[ni]->vData == NULL)
		{
			printf("Error: TileMapv2_BpmapToBAR, cannot allocate memory for loading bpmap coordinate data!\n");
			exit(EXIT_FAILURE);
		}
		nProbeNum = vProbePairNum->pMatElement[ni];
		if(nProbeNum > 0)
		{
			for(nj=0; nj<pBARPos->vSeqData[ni]->nColNum; nj++)
			{
				pBARPos->vSeqData[ni]->vData[nj] = CreateDoubleMatrix(1, nProbeNum);
				if(pBARPos->vSeqData[ni]->vData[nj] == NULL)
				{
					printf("Error: TileMapv2_BpmapToBAR, cannot allocate memory for loading bpmap coordinate data!\n");
					exit(EXIT_FAILURE);
				}
			}
		}

		pUnitList = NULL;
		
		/* load probes */
		for(nj=0; nj<nProbeNum; nj++)
		{
			/* load new unit */
			pNewUnit = NULL;
			pNewUnit = AffyBpMapUnitCreate();
			AffyBpMapUnitLoad_v3m2(pNewUnit, fpIn, vProbeMappingType->pMatElement[ni], little_endian_machine);
			
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

							fprintf(fpOut, "%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%f\t%d\n", 
								pBARPos->vSeqData[ni]->pSeqName->m_pString, pCUnit->nPos, pCUnit->nDepthNum, 
								pCUnit->nRepeatNum, pCUnit->strProbeSeq, 
								pCUnit->nPMX, pCUnit->nPMY, pCUnit->nMMX, pCUnit->nMMY,
								pCUnit->fMatchScore, pCUnit->bStrand);

							if(pCUnit->nRepeatNum < 2)
							{
								nk = pBARPos->vSeqData[ni]->nDataNum;
								pBARPos->vSeqData[ni]->vData[0]->pMatElement[nk] = pCUnit->nPos;
								pBARPos->vSeqData[ni]->vData[1]->pMatElement[nk] = pCUnit->nPMX;
								pBARPos->vSeqData[ni]->vData[2]->pMatElement[nk] = pCUnit->nPMY;
								if(vProbeMappingType->pMatElement[ni] == 0)
								{
									pBARPos->vSeqData[ni]->vData[3]->pMatElement[nk] = pCUnit->nMMX;
									pBARPos->vSeqData[ni]->vData[4]->pMatElement[nk] = pCUnit->nMMY;
								}
								pBARPos->vSeqData[ni]->nDataNum += 1;
							}

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

			fprintf(fpOut, "%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%f\t%d\n", 
						pBARPos->vSeqData[ni]->pSeqName->m_pString, pCUnit->nPos, pCUnit->nDepthNum, 
						pCUnit->nRepeatNum, pCUnit->strProbeSeq, 
						pCUnit->nPMX, pCUnit->nPMY, pCUnit->nMMX, pCUnit->nMMY,
						pCUnit->fMatchScore, pCUnit->bStrand);

			if(pCUnit->nRepeatNum < 2)
			{
				nk = pBARPos->vSeqData[ni]->nDataNum;
				pBARPos->vSeqData[ni]->vData[0]->pMatElement[nk] = pCUnit->nPos;
				pBARPos->vSeqData[ni]->vData[1]->pMatElement[nk] = pCUnit->nPMX;
				pBARPos->vSeqData[ni]->vData[2]->pMatElement[nk] = pCUnit->nPMY;
				if(vProbeMappingType->pMatElement[ni] == 0)
				{
					pBARPos->vSeqData[ni]->vData[3]->pMatElement[nk] = pCUnit->nMMX;
					pBARPos->vSeqData[ni]->vData[4]->pMatElement[nk] = pCUnit->nMMY;
				}
				pBARPos->vSeqData[ni]->nDataNum += 1;
			}

			AffyBpMapUnitDestroy(pCUnit);
		}


		for(nj=0; nj<pBARPos->vSeqData[ni]->nColNum; nj++)
		{
			if(pBARPos->vSeqData[ni]->nDataNum > 0)
			{
				pResizeMat = NULL;
				pResizeMat = CreateDoubleMatrix(1, pBARPos->vSeqData[ni]->nDataNum);
				if(pResizeMat == NULL)
				{
					printf("Error: TileMapv2_BpmapToBAR, cannot create an intermediate matrix to move data!\n");
					exit(EXIT_FAILURE);
				}
				memcpy(pResizeMat->pMatElement, pBARPos->vSeqData[ni]->vData[nj]->pMatElement, pBARPos->vSeqData[ni]->nDataNum*sizeof(double));
				DestroyDoubleMatrix(pBARPos->vSeqData[ni]->vData[nj]);
				pBARPos->vSeqData[ni]->vData[nj] = pResizeMat;
			}
			else
			{
				DestroyDoubleMatrix(pBARPos->vSeqData[ni]->vData[nj]);
				pBARPos->vSeqData[ni]->vData[nj] = NULL;
			}
		}

		/* update masked probe number */
		nMaskedProbeNum += pBARPos->vSeqData[ni]->nDataNum;
	}

	/* load tail if any */
	printf("  ");
	while(feof(fpIn) == 0)
	{
		fread(strFileType, sizeof(char), 1, fpIn);
		printf("%c", strFileType[0]);
	}
	printf("\n");

	/* clear memeory */
	DestroyIntMatrix(vProbeMappingType);
	if(fVersion > 2.5)
	{
		DestroyIntMatrix(vSequenceFileOffset);
	}
    DestroyIntMatrix(vProbePairNum);
	DestroyIntMatrix(vSeqID);

	/* close file */
	fclose(fpIn);
	fclose(fpOut);

	printf("  Probe number before masking = %d\n", nTotalProbeNum);
	printf("  Probe number after masking = %d\n", nMaskedProbeNum);

	/* return */
	return pBARPos;
}

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
			int nApplyMask, int nExportMode)
{
	/* define */
	FILE **vfpIn;
	FILE **vfpMask;
	char strLine[MED_LINE_LENGTH];
	int ni;
	int nidx;
	struct tagString **vExportPath;
	char strExportPath[MED_LINE_LENGTH];
	
	/* init */
	if( (nSampleNum <= 0) || (nLibNum <= 0) )
	{
		return PROC_SUCCESS;
	}
	if( (nLibId < 0) || (nLibId >= nLibNum) )
	{
		printf("Error: TileMapv2_ExportAffyIntensity, library id out of range!\n");
		exit(EXIT_FAILURE);
	}
	
	/* open cel files */
	vfpIn = NULL;
	vfpIn = (FILE **)calloc(nSampleNum, sizeof(FILE *));
	if(vfpIn == NULL)
	{
		printf("Error: TileMapv2_ExportAffyIntensity, cannot allocate memory for tracking file pointers!\n");
		exit(EXIT_FAILURE);
	}

	vfpMask = NULL;
	vfpMask = (FILE **)calloc(nSampleNum, sizeof(FILE *));
	if(vfpMask == NULL)
	{
		printf("Error: TileMapv2_ExportAffyIntensity, cannot allocate memory for tracking file pointers!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nSampleNum; ni++)
	{
		nidx = ni*nLibNum+nLibId;
		sprintf(strLine, "%s%s.int", strWorkPath, vCELPath[nidx]->m_pString);
		vfpIn[ni] = fopen(strLine, "rb");
		if(vfpIn[ni] == NULL)
		{
			printf("Error: TileMapv2_ExportAffyIntensity, cannot open %s\n", strLine);
			exit(EXIT_FAILURE);
		}

		if(nApplyMask == 1)
		{
			sprintf(strLine, "%s%s.mask", strWorkPath, vCELPath[nidx]->m_pString);
			vfpMask[ni] = fopen(strLine, "rb");
			if(vfpMask[ni] == NULL)
			{
				printf("Error: TileMapv2_ExportAffyIntensity, cannot open %s\n", strLine);
				exit(EXIT_FAILURE);
			}
		}
	}

	/* export */
	if( nExportMode == 0 )
	{
		/* separate bar files */
		vExportPath = NULL;
		vExportPath = (struct tagString **)calloc(nSampleNum, sizeof(struct tagString *));
		if(vExportPath == NULL)
		{
			printf("Error: TileMapv2_ExportAffyIntensity, cannot allocate memory for naming output files!\n");
			exit(EXIT_FAILURE);
		}
		
		for(ni=0; ni<nSampleNum; ni++)
		{
			nidx = ni*nLibNum+nLibId;
			sprintf(strLine, "%s%s.bar", strWorkPath, vCELPath[nidx]->m_pString);
			StringAddTail(vExportPath+ni, strLine);
		}

		TileMapv2_ExportAffyIntensity_MultiBAR(vExportPath, 
			nSampleNum, vSampleAlias, vfpIn, 
			nCELNumberCells, nCELTotalX, nCELTotalY,
			nApplyMask, vfpMask, pBARPos, 
			nIntensityType,	dIntLowerBound, nIntLogTransform);

		for(ni=0; ni<nSampleNum; ni++)
		{
			DeleteString(vExportPath[ni]);
		}
		free(vExportPath);
	}
	else if( nExportMode == 1)
	{
		/* a combined bar file */
		sprintf(strExportPath, "%s%s_%d.bar", strWorkPath, strProjectTitle, (nLibId+1));
		printf("Exporting %s...\n", strExportPath);
		TileMapv2_ExportAffyIntensity_SingleBAR(strExportPath, 
			nSampleNum, vSampleAlias, vfpIn, 
			nCELNumberCells, nCELTotalX, nCELTotalY,
			nApplyMask, vfpMask, pBARPos, 
			nIntensityType,	dIntLowerBound, nIntLogTransform);
	}
	else
	{
		/* a combined txt file */
		sprintf(strExportPath, "%s%s_%d.txt", strWorkPath, strProjectTitle, (nLibId+1));
		printf("Exporting %s...\n", strExportPath);
		TileMapv2_ExportAffyIntensity_SingleTXT(strExportPath, 
			nSampleNum, vSampleAlias, vfpIn, 
			nCELNumberCells, nCELTotalX, nCELTotalY,
			nApplyMask, vfpMask, pBARPos, 
			nIntensityType,	dIntLowerBound, nIntLogTransform);
	}

	/* close files */
	for(ni=0; ni<nSampleNum; ni++)
	{
		fclose(vfpIn[ni]);
		if(nApplyMask == 1)
		{
			fclose(vfpMask[ni]);
		}
	}

	/* release memory */
	free(vfpIn);
	free(vfpMask);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ExportAffyIntensity_MultiBAR()                               */
/*  Export intensity data to separate bar files.                           */
/* ----------------------------------------------------------------------- */
int TileMapv2_ExportAffyIntensity_MultiBAR(struct tagString **vExportPath, 
			int nSampleNum, struct tagString **vSampleAlias, FILE **vfpIn, 
			int nCELNumberCells, int nCELTotalX, int nCELTotalY,
			int nApplyMask, FILE **vfpMask,	struct tagBARData *pBARPos, 
			int nIntensityType,	double dIntLowerBound, int nIntLogTransform)
{
	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));
	FILE **vfpOut;
	FILE **vfpOutMask;
	char strLine[MED_LINE_LENGTH];
	int vFieldType[2];
	int ni,nj,nk;
	float fV;
	int nV;
	struct DOUBLEMATRIX *pCEL;
	struct BYTEMATRIX *pMask;
	int nidx1,nidx2;
	double dLog2 = log(2.0);
	int nMask;
	unsigned char bMask;
	int nLen;
	double dPM,dMM;
	
	/* init */
	if( (nSampleNum <= 0) || (nCELNumberCells <= 0))
	{
		return PROC_SUCCESS;
	}

	pCEL = NULL;
	pCEL = CreateDoubleMatrix(1, nCELNumberCells);
	if(pCEL == NULL)
	{
		printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot allocate memory for loading data!\n");
		exit(EXIT_FAILURE);
	}

	pMask = NULL;
	pMask = CreateByteMatrix(1, nCELNumberCells);
	if(pMask == NULL)
	{
		printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot allocate memory for loading data!\n");
		exit(EXIT_FAILURE);
	}

	
	/* create cel file memory */
	vfpOut = NULL;
	vfpOut = (FILE **)calloc(nSampleNum, sizeof(FILE *));
	if(vfpOut == NULL)
	{
		printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot allocate memory for tracking file pointers!\n");
		exit(EXIT_FAILURE);
	}

	vfpOutMask = NULL;
	vfpOutMask = (FILE **)calloc(nSampleNum, sizeof(FILE *));
	if(vfpOutMask == NULL)
	{
		printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot allocate memory for tracking file pointers!\n");
		exit(EXIT_FAILURE);
	}

	/* write files */
	for(ni=0; ni<nSampleNum; ni++)
	{
		/* load CEL data */
		if( little_endian_fread(pCEL->pMatElement, sizeof(double), nCELNumberCells, vfpIn[ni], little_endian_machine) != nCELNumberCells)
		{
			printf("Error: TileMapv2_LoadFromBinaryFile, cannot load int data from binary files correctly!\n");
			exit(EXIT_FAILURE);
		}

		if(nApplyMask == 1)
		{
			if( little_endian_fread(pMask->pMatElement, sizeof(unsigned char), nCELNumberCells, vfpMask[ni], little_endian_machine) != nCELNumberCells)
			{
				printf("Error: TileMapv2_LoadFromBinaryFile, cannot load mask data from binary files correctly!\n");
				exit(EXIT_FAILURE);
			}
		}

		/* open files */
		printf("Exporting %s...\n", vExportPath[ni]->m_pString);
		vfpOut[ni] = fopen(vExportPath[ni]->m_pString, "wb");
		if(vfpOut[ni] == NULL)
		{
			printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot open output file!\n");
			exit(EXIT_FAILURE);
		}

		if(nApplyMask == 1)
		{
			sprintf(strLine, "%s.mask", vExportPath[ni]->m_pString);
			vfpOutMask[ni] = fopen(strLine, "wb");
			if(vfpOutMask[ni] == NULL)
			{
				printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot open output file!\n");
				exit(EXIT_FAILURE);
			}
		}

		/* export */
		strcpy(strLine, "barr\r\n\032\n");
		if(big_endian_fwrite(strLine, 1, 8, vfpOut[ni], little_endian_machine) != 8)
		{
			printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write magic number.\n");
			exit(EXIT_FAILURE);
		}
		if(nApplyMask == 1)
		{
			if(big_endian_fwrite(strLine, 1, 8, vfpOutMask[ni], little_endian_machine) != 8)
			{
				printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write magic number.\n");
				exit(EXIT_FAILURE);
			}
		}

		/* write version number */
		if(big_endian_fwrite(&(pBARPos->fVersionnumber), FLOAT_SIZE, 1, vfpOut[ni], little_endian_machine) != 1)
		{
			printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write version number.\n");
			exit(EXIT_FAILURE);
		}
		if(nApplyMask == 1)
		{
			if(big_endian_fwrite(&(pBARPos->fVersionnumber), FLOAT_SIZE, 1, vfpOutMask[ni], little_endian_machine) != 1)
			{
				printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write version number.\n");
				exit(EXIT_FAILURE);
			}
		}

		/* write sequence number */
		if(big_endian_fwrite(&(pBARPos->nSeqNum), INT_SIZE, 1, vfpOut[ni], little_endian_machine) != 1)
		{
			printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write sequence number.\n");
			exit(EXIT_FAILURE);
		}
		if(nApplyMask == 1)
		{
			if(big_endian_fwrite(&(pBARPos->nSeqNum), INT_SIZE, 1, vfpOutMask[ni], little_endian_machine) != 1)
			{
				printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write sequence number.\n");
				exit(EXIT_FAILURE);
			}
		}

		/* write column number */
		nV = 2;
		if(big_endian_fwrite(&nV, INT_SIZE, 1, vfpOut[ni], little_endian_machine) != 1)
		{
			printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write column number.\n");
			exit(EXIT_FAILURE);
		}
		if(nApplyMask == 1)
		{
			if(big_endian_fwrite(&nV, INT_SIZE, 1, vfpOutMask[ni], little_endian_machine) != 1)
			{
				printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write column number.\n");
				exit(EXIT_FAILURE);
			}
		}

		/* write column type */
		vFieldType[0] = 2;
		vFieldType[1] = 1;
		if((int)big_endian_fwrite(&vFieldType, INT_SIZE, 2, vfpOut[ni], little_endian_machine) != 2)
		{
			printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write field type correctly.\n");
			exit(EXIT_FAILURE);
		}
		if(nApplyMask == 1)
		{
			vFieldType[0] = 2;
			vFieldType[1] = 7;
			if((int)big_endian_fwrite(&vFieldType, INT_SIZE, 2, vfpOutMask[ni], little_endian_machine) != 2)
			{
				printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write field type correctly.\n");
				exit(EXIT_FAILURE);
			}
		}
		
		/* write parameter name/value pairs */
		if(big_endian_fwrite(&(pBARPos->nParamNum), INT_SIZE, 1, vfpOut[ni], little_endian_machine) != 1)
		{
			printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write number of name/value pairs.\n");
			exit(EXIT_FAILURE);
		}
		if(nApplyMask == 1)
		{
			if(big_endian_fwrite(&(pBARPos->nParamNum), INT_SIZE, 1, vfpOutMask[ni], little_endian_machine) != 1)
			{
				printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write number of name/value pairs.\n");
				exit(EXIT_FAILURE);
			}
		}

		if(pBARPos->nParamNum > 0)
		{
			if( (pBARPos->vParamName == NULL) || (pBARPos->vParamValue == NULL) )
			{
				printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot find parameter name/value pairs\n");
				exit(EXIT_FAILURE);
			}

			for(nj=0; nj<pBARPos->nParamNum; nj++)
			{
				if(pBARPos->vParamName[nj] == NULL)
				{
					nLen = 0;
					if(big_endian_fwrite(&nLen, INT_SIZE, 1, vfpOut[ni], little_endian_machine) != 1)
					{
						printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write parameter name length.\n");
						exit(EXIT_FAILURE);
					}
					if(nApplyMask == 1)
					{
						if(big_endian_fwrite(&nLen, INT_SIZE, 1, vfpOutMask[ni], little_endian_machine) != 1)
						{
							printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write parameter name length.\n");
							exit(EXIT_FAILURE);
						}
					}
				}
				else
				{
					if(big_endian_fwrite(&(pBARPos->vParamName[nj]->m_nLength), INT_SIZE, 1, vfpOut[ni], little_endian_machine) != 1)
					{
						printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write parameter name length.\n");
						exit(EXIT_FAILURE);
					}
					
					if(big_endian_fwrite(pBARPos->vParamName[nj]->m_pString, 1, pBARPos->vParamName[nj]->m_nLength, vfpOut[ni], little_endian_machine) != pBARPos->vParamName[nj]->m_nLength)
					{
						printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write parameter name.\n");
						exit(EXIT_FAILURE);
					}

					if(nApplyMask == 1)
					{
						if(big_endian_fwrite(&(pBARPos->vParamName[nj]->m_nLength), INT_SIZE, 1, vfpOutMask[ni], little_endian_machine) != 1)
						{
							printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write parameter name length.\n");
							exit(EXIT_FAILURE);
						}
						
						if(big_endian_fwrite(pBARPos->vParamName[nj]->m_pString, 1, pBARPos->vParamName[nj]->m_nLength, vfpOutMask[ni], little_endian_machine) != pBARPos->vParamName[nj]->m_nLength)
						{
							printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write parameter name.\n");
							exit(EXIT_FAILURE);
						}
					}
				}

				if(pBARPos->vParamValue[nj] == NULL)
				{
					nLen = 0;
					if(big_endian_fwrite(&nLen, INT_SIZE, 1, vfpOut[ni], little_endian_machine) != 1)
					{
						printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write parameter name length.\n");
						exit(EXIT_FAILURE);
					}
					if(nApplyMask == 1)
					{
						if(big_endian_fwrite(&nLen, INT_SIZE, 1, vfpOutMask[ni], little_endian_machine) != 1)
						{
							printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write parameter name length.\n");
							exit(EXIT_FAILURE);
						}
					}
				}
				else
				{
					if(big_endian_fwrite(&(pBARPos->vParamValue[nj]->m_nLength), INT_SIZE, 1, vfpOut[ni], little_endian_machine) != 1)
					{
						printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write parameter value length.\n");
						exit(EXIT_FAILURE);
					}

					if(big_endian_fwrite(pBARPos->vParamValue[nj]->m_pString, 1, pBARPos->vParamValue[nj]->m_nLength, vfpOut[ni], little_endian_machine) != pBARPos->vParamValue[nj]->m_nLength)
					{
						printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write parameter value.\n");
						exit(EXIT_FAILURE);
					}

					if(nApplyMask == 1)
					{
						if(big_endian_fwrite(&(pBARPos->vParamValue[nj]->m_nLength), INT_SIZE, 1, vfpOutMask[ni], little_endian_machine) != 1)
						{
							printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write parameter value length.\n");
							exit(EXIT_FAILURE);
						}

						if(big_endian_fwrite(pBARPos->vParamValue[nj]->m_pString, 1, pBARPos->vParamValue[nj]->m_nLength, vfpOutMask[ni], little_endian_machine) != pBARPos->vParamValue[nj]->m_nLength)
						{
							printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write parameter value.\n");
							exit(EXIT_FAILURE);
						}
					}
				}
			}
		}

		/* write sequence data */
		for(nj=0; nj<pBARPos->nSeqNum; nj++)
		{
			/* create BARSeq object */
			if(pBARPos->vSeqData[nj] == NULL)
			{
				printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot find BARSeq object.\n");
				exit(EXIT_FAILURE);
			}

			/* write sequence name */
			if(pBARPos->vSeqData[nj]->pSeqName == NULL)
			{
				nLen = 0;
				if(big_endian_fwrite(&nLen, INT_SIZE, 1, vfpOut[ni], little_endian_machine) != 1)
				{
					printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write sequence name length.\n");
					exit(EXIT_FAILURE);
				}
				if(nApplyMask == 1)
				{
					if(big_endian_fwrite(&nLen, INT_SIZE, 1, vfpOutMask[ni], little_endian_machine) != 1)
					{
						printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write sequence name length.\n");
						exit(EXIT_FAILURE);
					}
				}
			}
			else
			{
				if(big_endian_fwrite(&(pBARPos->vSeqData[nj]->pSeqName->m_nLength), INT_SIZE, 1, vfpOut[ni], little_endian_machine) != 1)
				{
					printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write sequence name length.\n");
					exit(EXIT_FAILURE);
				}
				
				if(big_endian_fwrite(pBARPos->vSeqData[nj]->pSeqName->m_pString, 1, pBARPos->vSeqData[nj]->pSeqName->m_nLength, vfpOut[ni], little_endian_machine) != pBARPos->vSeqData[nj]->pSeqName->m_nLength)
				{
					printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write sequence name.\n");
					exit(EXIT_FAILURE);
				}

				if(nApplyMask == 1)
				{
					if(big_endian_fwrite(&(pBARPos->vSeqData[nj]->pSeqName->m_nLength), INT_SIZE, 1, vfpOutMask[ni], little_endian_machine) != 1)
					{
						printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write sequence name length.\n");
						exit(EXIT_FAILURE);
					}
					
					if(big_endian_fwrite(pBARPos->vSeqData[nj]->pSeqName->m_pString, 1, pBARPos->vSeqData[nj]->pSeqName->m_nLength, vfpOutMask[ni], little_endian_machine) != pBARPos->vSeqData[nj]->pSeqName->m_nLength)
					{
						printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write sequence name.\n");
						exit(EXIT_FAILURE);
					}
				}
			}
			
			/* write group name */
			if(pBARPos->fVersionnumber > 1.5)
			{
				if(pBARPos->vSeqData[nj]->pSeqGroupName == NULL)
				{
					nLen = 0;
					if(big_endian_fwrite(&nLen, INT_SIZE, 1, vfpOut[ni], little_endian_machine) != 1)
					{
						printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write sequence group name length.\n");
						exit(EXIT_FAILURE);
					}

					if(nApplyMask == 1)
					{
						if(big_endian_fwrite(&nLen, INT_SIZE, 1, vfpOutMask[ni], little_endian_machine) != 1)
						{
							printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write sequence group name length.\n");
							exit(EXIT_FAILURE);
						}
					}
				}
				else
				{
					if(big_endian_fwrite(&(pBARPos->vSeqData[nj]->pSeqGroupName->m_nLength), INT_SIZE, 1, vfpOut[ni], little_endian_machine) != 1)
					{
						printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write sequence group name length.\n");
						exit(EXIT_FAILURE);
					}

					if(big_endian_fwrite(pBARPos->vSeqData[nj]->pSeqGroupName->m_pString, 1, pBARPos->vSeqData[nj]->pSeqGroupName->m_nLength, vfpOut[ni], little_endian_machine) != pBARPos->vSeqData[nj]->pSeqGroupName->m_nLength)
					{
						printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write sequence group name.\n");
						exit(EXIT_FAILURE);
					}

					if(nApplyMask == 1)
					{
						if(big_endian_fwrite(&(pBARPos->vSeqData[nj]->pSeqGroupName->m_nLength), INT_SIZE, 1, vfpOutMask[ni], little_endian_machine) != 1)
						{
							printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write sequence group name length.\n");
							exit(EXIT_FAILURE);
						}

						if(big_endian_fwrite(pBARPos->vSeqData[nj]->pSeqGroupName->m_pString, 1, pBARPos->vSeqData[nj]->pSeqGroupName->m_nLength, vfpOutMask[ni], little_endian_machine) != pBARPos->vSeqData[nj]->pSeqGroupName->m_nLength)
						{
							printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write sequence group name.\n");
							exit(EXIT_FAILURE);
						}
					}
				}
			}

			/* write sequence version */
			if(pBARPos->vSeqData[nj]->pSeqVersion == NULL)
			{
				nLen = 0;
				if(big_endian_fwrite(&nLen, INT_SIZE, 1, vfpOut[ni], little_endian_machine) != 1)
				{
					printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write sequence version length.\n");
					exit(EXIT_FAILURE);
				}
				if(nApplyMask == 1)
				{
					if(big_endian_fwrite(&nLen, INT_SIZE, 1, vfpOutMask[ni], little_endian_machine) != 1)
					{
						printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write sequence version length.\n");
						exit(EXIT_FAILURE);
					}
				}
			}
			else
			{
				if(big_endian_fwrite(&(pBARPos->vSeqData[nj]->pSeqVersion->m_nLength), INT_SIZE, 1, vfpOut[ni], little_endian_machine) != 1)
				{
					printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write sequence version length.\n");
					exit(EXIT_FAILURE);
				}
				
				if(big_endian_fwrite(pBARPos->vSeqData[nj]->pSeqVersion->m_pString, 1, pBARPos->vSeqData[nj]->pSeqVersion->m_nLength, vfpOut[ni], little_endian_machine) != pBARPos->vSeqData[nj]->pSeqVersion->m_nLength)
				{
					printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write sequence version.\n");
					exit(EXIT_FAILURE);
				}

				if(nApplyMask == 1)
				{
					if(big_endian_fwrite(&(pBARPos->vSeqData[nj]->pSeqVersion->m_nLength), INT_SIZE, 1, vfpOutMask[ni], little_endian_machine) != 1)
					{
						printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write sequence version length.\n");
						exit(EXIT_FAILURE);
					}
					
					if(big_endian_fwrite(pBARPos->vSeqData[nj]->pSeqVersion->m_pString, 1, pBARPos->vSeqData[nj]->pSeqVersion->m_nLength, vfpOutMask[ni], little_endian_machine) != pBARPos->vSeqData[nj]->pSeqVersion->m_nLength)
					{
						printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write sequence version.\n");
						exit(EXIT_FAILURE);
					}
				}
			}
			
			/* write parameter name/value pairs */
			if(pBARPos->fVersionnumber > 1.5)
			{
				if(big_endian_fwrite(&(pBARPos->vSeqData[nj]->nParamNum), INT_SIZE, 1, vfpOut[ni], little_endian_machine) != 1)
				{
					printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write number of name/value pairs for BARseq object.\n");
					exit(EXIT_FAILURE);
				}
				if(nApplyMask == 1)
				{
					if(big_endian_fwrite(&(pBARPos->vSeqData[nj]->nParamNum), INT_SIZE, 1, vfpOutMask[ni], little_endian_machine) != 1)
					{
						printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write number of name/value pairs for BARseq object.\n");
						exit(EXIT_FAILURE);
					}
				}

				if(pBARPos->vSeqData[nj]->nParamNum > 0)
				{
					if( (pBARPos->vSeqData[nj]->vParamName == NULL) || (pBARPos->vSeqData[nj]->vParamValue == NULL) )
					{
						printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot find sequence parameter name/value pairs\n");
						exit(EXIT_FAILURE);
					}

					for(nk=0; nk<pBARPos->vSeqData[nj]->nParamNum; nk++)
					{
						if(pBARPos->vSeqData[nj]->vParamName[nk] == NULL)
						{
							nLen = 0;
							if(big_endian_fwrite(&nLen, INT_SIZE, 1, vfpOut[ni], little_endian_machine) != 1)
							{
								printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write sequence parameter name length.\n");
								exit(EXIT_FAILURE);
							}

							if(nApplyMask == 1)
							{
								if(big_endian_fwrite(&nLen, INT_SIZE, 1, vfpOutMask[ni], little_endian_machine) != 1)
								{
									printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write sequence parameter name length.\n");
									exit(EXIT_FAILURE);
								}
							}
						}
						else
						{
							if(big_endian_fwrite(&(pBARPos->vSeqData[nj]->vParamName[nk]->m_nLength), INT_SIZE, 1, vfpOut[ni], little_endian_machine) != 1)
							{
								printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write sequence parameter name length.\n");
								exit(EXIT_FAILURE);
							}

							if(big_endian_fwrite(pBARPos->vSeqData[nj]->vParamName[nk]->m_pString, 1, pBARPos->vSeqData[nj]->vParamName[nk]->m_nLength, vfpOut[ni], little_endian_machine) != pBARPos->vSeqData[nj]->vParamName[nk]->m_nLength)
							{
								printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write sequence parameter name.\n");
								exit(EXIT_FAILURE);
							}

							if(nApplyMask == 1)
							{
								if(big_endian_fwrite(&(pBARPos->vSeqData[nj]->vParamName[nk]->m_nLength), INT_SIZE, 1, vfpOutMask[ni], little_endian_machine) != 1)
								{
									printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write sequence parameter name length.\n");
									exit(EXIT_FAILURE);
								}

								if(big_endian_fwrite(pBARPos->vSeqData[nj]->vParamName[nk]->m_pString, 1, pBARPos->vSeqData[nj]->vParamName[nk]->m_nLength, vfpOutMask[ni], little_endian_machine) != pBARPos->vSeqData[nj]->vParamName[nk]->m_nLength)
								{
									printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write sequence parameter name.\n");
									exit(EXIT_FAILURE);
								}
							}
						}
						
						if(pBARPos->vSeqData[nj]->vParamValue[nk] == NULL)
						{
							nLen = 0;
							if(big_endian_fwrite(&nLen, INT_SIZE, 1, vfpOut[ni], little_endian_machine) != 1)
							{
								printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write sequence parameter value length.\n");
								exit(EXIT_FAILURE);
							}

							if(nApplyMask == 1)
							{
								if(big_endian_fwrite(&nLen, INT_SIZE, 1, vfpOutMask[ni], little_endian_machine) != 1)
								{
									printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write sequence parameter value length.\n");
									exit(EXIT_FAILURE);
								}
							}
						}
						else
						{
							if(big_endian_fwrite(&(pBARPos->vSeqData[nj]->vParamValue[nk]->m_nLength), INT_SIZE, 1, vfpOut[ni], little_endian_machine) != 1)
							{
								printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write sequence parameter value length.\n");
								exit(EXIT_FAILURE);
							}

							if(big_endian_fwrite(pBARPos->vSeqData[nj]->vParamValue[nk]->m_pString, 1, pBARPos->vSeqData[nj]->vParamValue[nk]->m_nLength, vfpOut[ni], little_endian_machine) != pBARPos->vSeqData[nj]->vParamValue[nk]->m_nLength)
							{
								printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write parameter value.\n");
								exit(EXIT_FAILURE);
							}

							if(nApplyMask == 1)
							{
								if(big_endian_fwrite(&(pBARPos->vSeqData[nj]->vParamValue[nk]->m_nLength), INT_SIZE, 1, vfpOutMask[ni], little_endian_machine) != 1)
								{
									printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write sequence parameter value length.\n");
									exit(EXIT_FAILURE);
								}

								if(big_endian_fwrite(pBARPos->vSeqData[nj]->vParamValue[nk]->m_pString, 1, pBARPos->vSeqData[nj]->vParamValue[nk]->m_nLength, vfpOutMask[ni], little_endian_machine) != pBARPos->vSeqData[nj]->vParamValue[nk]->m_nLength)
								{
									printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write parameter value.\n");
									exit(EXIT_FAILURE);
								}
							}
						}
					}
				}
			}

			/* write data points */
			if(big_endian_fwrite(&(pBARPos->vSeqData[nj]->nDataNum), INT_SIZE, 1, vfpOut[ni], little_endian_machine) != 1)
			{
				printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write number of data points for a sequence.\n");
				exit(EXIT_FAILURE);
			}

			if(nApplyMask == 1)
			{
				if(big_endian_fwrite(&(pBARPos->vSeqData[nj]->nDataNum), INT_SIZE, 1, vfpOutMask[ni], little_endian_machine) != 1)
				{
					printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot write number of data points for a sequence.\n");
					exit(EXIT_FAILURE);
				}
			}

			if(pBARPos->vSeqData[nj]->nColNum > 0)
			{
				if(pBARPos->vSeqData[nj]->vData == NULL)
				{
					printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot find sequence data.\n");
					exit(EXIT_FAILURE);
				}
			}

			if(pBARPos->vSeqData[nj]->nDataNum > 0)
			{
				for(nk=0; nk<pBARPos->vSeqData[nj]->nColNum; nk++)
				{
					if(pBARPos->vSeqData[nj]->vData[nk] == NULL)
					{
						printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, cannot find sequence data.\n");
						exit(EXIT_FAILURE);
					}
				}

				for(nk=0; nk<pBARPos->vSeqData[nj]->nDataNum; nk++)
				{
					nidx1 = (int)(pBARPos->vSeqData[nj]->vData[2]->pMatElement[nk])*nCELTotalX+(int)(pBARPos->vSeqData[nj]->vData[1]->pMatElement[nk]);
					if(nidx1 >= nCELNumberCells)
					{
						printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, index out of range\n");
						exit(EXIT_FAILURE);
					}

					dPM = pCEL->pMatElement[nidx1];

					if(pBARPos->vSeqData[nj]->nColNum >= 5)
					{
						nidx2 = (int)(pBARPos->vSeqData[nj]->vData[4]->pMatElement[nk])*nCELTotalX+(int)(pBARPos->vSeqData[nj]->vData[3]->pMatElement[nk]);
						if(nidx2 >= nCELNumberCells)
						{
							printf("Error: TileMapv2_ExportAffyIntensity_MultiBAR, index out of range\n");
							exit(EXIT_FAILURE);
						}

						dMM = pCEL->pMatElement[nidx2];
					}
					else
					{
						dMM = 0.0;
					}

					if(nIntensityType == 1)
					{
						dPM -= dMM;
					}
							
					if(dPM < dIntLowerBound)
						dPM = dIntLowerBound;

					if(nIntLogTransform == 1)
					{
						dPM = log(dPM)/dLog2;
					}

					nV = (int)(pBARPos->vSeqData[nj]->vData[0]->pMatElement[nk]);
					big_endian_fwrite(&nV, INT_SIZE, 1, vfpOut[ni], little_endian_machine);
					fV = (float)dPM;
					big_endian_fwrite(&fV, FLOAT_SIZE, 1, vfpOut[ni], little_endian_machine);

					/* mask */
					if(nApplyMask == 1)
					{
						if(pMask->pMatElement[nidx1] == 1)
							nMask = 1;
						else
							nMask = 0;

						if(pBARPos->vSeqData[nj]->nColNum >= 5)
						{
							if( (pMask->pMatElement[nidx2] == 1) && (nIntensityType == 1) )
								nMask = 1;
						}

						if(nMask == 1)
							bMask = 1;
						else
							bMask = 0;

						big_endian_fwrite(&nV, INT_SIZE, 1, vfpOutMask[ni], little_endian_machine);
						big_endian_fwrite(&bMask, 1, 1, vfpOutMask[ni], little_endian_machine);
					}
				}
			}
		}

		/* close files */
		fclose(vfpOut[ni]);
		if(nApplyMask == 1)
		{
			fclose(vfpOutMask[ni]);
		}
	}

	/* release memory */
	DestroyDoubleMatrix(pCEL);
	DestroyByteMatrix(pMask);
	free(vfpOut);
	free(vfpOutMask);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ExportAffyIntensity_SingleBAR()                              */
/*  Export intensity data to a single bar files.                           */
/* ----------------------------------------------------------------------- */
int TileMapv2_ExportAffyIntensity_SingleBAR(char strExportPath[], 
			int nSampleNum, struct tagString **vSampleAlias, FILE **vfpIn,
			int nCELNumberCells, int nCELTotalX, int nCELTotalY,
			int nApplyMask, FILE **vfpMask,	struct tagBARData *pBARPos, 
			int nIntensityType,	double dIntLowerBound, int nIntLogTransform)
{
	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));
	FILE *fpOut;
	FILE *fpOutMask;
	char strLine[MED_LINE_LENGTH];
	int nFieldType;
	int ni,nj,nk;
	float fV;
	int nV;
	int nidx1,nidx2;
	double dLog2 = log(2.0);
	int nMask;
	int nLen;
	double dPM,dMM;
	unsigned char bMask;
	
	/* init */
	if( (nSampleNum <= 0) || (nCELNumberCells <= 0))
	{
		return PROC_SUCCESS;
	}

	/* open files */
	fpOut = NULL;
	fpOutMask = NULL;

	fpOut = fopen(strExportPath, "wb");
	if(fpOut == NULL)
	{
		printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	
	if(nApplyMask == 1)
	{
		sprintf(strLine, "%s.mask", strExportPath);
		fpOutMask = fopen(strLine, "wb");
		if(fpOutMask == NULL)
		{
			printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot open output file!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	
	/* write files */
	strcpy(strLine, "barr\r\n\032\n");
	if(big_endian_fwrite(strLine, 1, 8, fpOut, little_endian_machine) != 8)
	{
		printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write magic number.\n");
		exit(EXIT_FAILURE);
	}
	if(nApplyMask == 1)
	{
		if(big_endian_fwrite(strLine, 1, 8, fpOutMask, little_endian_machine) != 8)
		{
			printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write magic number.\n");
			exit(EXIT_FAILURE);
		}
	}

	/* write version number */
	if(big_endian_fwrite(&(pBARPos->fVersionnumber), FLOAT_SIZE, 1, fpOut, little_endian_machine) != 1)
	{
		printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write version number.\n");
		exit(EXIT_FAILURE);
	}
	if(nApplyMask == 1)
	{
		if(big_endian_fwrite(&(pBARPos->fVersionnumber), FLOAT_SIZE, 1, fpOutMask, little_endian_machine) != 1)
		{
			printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write version number.\n");
			exit(EXIT_FAILURE);
		}
	}

	/* write sequence number */
	if(big_endian_fwrite(&(pBARPos->nSeqNum), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
	{
		printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write sequence number.\n");
		exit(EXIT_FAILURE);
	}
	if(nApplyMask == 1)
	{
		if(big_endian_fwrite(&(pBARPos->nSeqNum), INT_SIZE, 1, fpOutMask, little_endian_machine) != 1)
		{
			printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write sequence number.\n");
			exit(EXIT_FAILURE);
		}
	}

	/* write column number */
	nV = 1+nSampleNum;
	if(big_endian_fwrite(&nV, INT_SIZE, 1, fpOut, little_endian_machine) != 1)
	{
		printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write column number.\n");
		exit(EXIT_FAILURE);
	}
	if(nApplyMask == 1)
	{
		if(big_endian_fwrite(&nV, INT_SIZE, 1, fpOutMask, little_endian_machine) != 1)
		{
			printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write column number.\n");
			exit(EXIT_FAILURE);
		}
	}

	/* write column type */
	nFieldType = 2;
	if((int)big_endian_fwrite(&nFieldType, INT_SIZE, 1, fpOut, little_endian_machine) != 1)
	{
		printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write field type correctly.\n");
		exit(EXIT_FAILURE);
	}
	nFieldType = 1;
	for(ni=0; ni<nSampleNum; ni++)
	{
		if((int)big_endian_fwrite(&nFieldType, INT_SIZE, 1, fpOut, little_endian_machine) != 1)
		{
			printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write field type correctly.\n");
			exit(EXIT_FAILURE);
		}
	}

	if(nApplyMask == 1)
	{
		nFieldType = 2;
		if((int)big_endian_fwrite(&nFieldType, INT_SIZE, 1, fpOutMask, little_endian_machine) != 1)
		{
			printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write field type correctly.\n");
			exit(EXIT_FAILURE);
		}
		nFieldType = 7;
		for(ni=0; ni<nSampleNum; ni++)
		{
			if((int)big_endian_fwrite(&nFieldType, INT_SIZE, 1, fpOutMask, little_endian_machine) != 1)
			{
				printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write field type correctly.\n");
				exit(EXIT_FAILURE);
			}
		}
	}
	
	/* write parameter name/value pairs */
	if(big_endian_fwrite(&(pBARPos->nParamNum), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
	{
		printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write number of name/value pairs.\n");
		exit(EXIT_FAILURE);
	}
	if(nApplyMask == 1)
	{
		if(big_endian_fwrite(&(pBARPos->nParamNum), INT_SIZE, 1, fpOutMask, little_endian_machine) != 1)
		{
			printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write number of name/value pairs.\n");
			exit(EXIT_FAILURE);
		}
	}

	if(pBARPos->nParamNum > 0)
	{
		if( (pBARPos->vParamName == NULL) || (pBARPos->vParamValue == NULL) )
		{
			printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot find parameter name/value pairs\n");
			exit(EXIT_FAILURE);
		}

		for(nj=0; nj<pBARPos->nParamNum; nj++)
		{
			if(pBARPos->vParamName[nj] == NULL)
			{
				nLen = 0;
				if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOut, little_endian_machine) != 1)
				{
					printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write parameter name length.\n");
					exit(EXIT_FAILURE);
				}
				if(nApplyMask == 1)
				{
					if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOutMask, little_endian_machine) != 1)
					{
						printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write parameter name length.\n");
						exit(EXIT_FAILURE);
					}
				}
			}
			else
			{
				if(big_endian_fwrite(&(pBARPos->vParamName[nj]->m_nLength), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
				{
					printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write parameter name length.\n");
					exit(EXIT_FAILURE);
				}
				
				if(big_endian_fwrite(pBARPos->vParamName[nj]->m_pString, 1, pBARPos->vParamName[nj]->m_nLength, fpOut, little_endian_machine) != pBARPos->vParamName[nj]->m_nLength)
				{
					printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write parameter name.\n");
					exit(EXIT_FAILURE);
				}

				if(nApplyMask == 1)
				{
					if(big_endian_fwrite(&(pBARPos->vParamName[nj]->m_nLength), INT_SIZE, 1, fpOutMask, little_endian_machine) != 1)
					{
						printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write parameter name length.\n");
						exit(EXIT_FAILURE);
					}
					
					if(big_endian_fwrite(pBARPos->vParamName[nj]->m_pString, 1, pBARPos->vParamName[nj]->m_nLength, fpOutMask, little_endian_machine) != pBARPos->vParamName[nj]->m_nLength)
					{
						printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write parameter name.\n");
						exit(EXIT_FAILURE);
					}
				}
			}

			if(pBARPos->vParamValue[nj] == NULL)
			{
				nLen = 0;
				if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOut, little_endian_machine) != 1)
				{
					printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write parameter name length.\n");
					exit(EXIT_FAILURE);
				}
				if(nApplyMask == 1)
				{
					if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOutMask, little_endian_machine) != 1)
					{
						printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write parameter name length.\n");
						exit(EXIT_FAILURE);
					}
				}
			}
			else
			{
				if(big_endian_fwrite(&(pBARPos->vParamValue[nj]->m_nLength), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
				{
					printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write parameter value length.\n");
					exit(EXIT_FAILURE);
				}

				if(big_endian_fwrite(pBARPos->vParamValue[nj]->m_pString, 1, pBARPos->vParamValue[nj]->m_nLength, fpOut, little_endian_machine) != pBARPos->vParamValue[nj]->m_nLength)
				{
					printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write parameter value.\n");
					exit(EXIT_FAILURE);
				}

				if(nApplyMask == 1)
				{
					if(big_endian_fwrite(&(pBARPos->vParamValue[nj]->m_nLength), INT_SIZE, 1, fpOutMask, little_endian_machine) != 1)
					{
						printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write parameter value length.\n");
						exit(EXIT_FAILURE);
					}

					if(big_endian_fwrite(pBARPos->vParamValue[nj]->m_pString, 1, pBARPos->vParamValue[nj]->m_nLength, fpOutMask, little_endian_machine) != pBARPos->vParamValue[nj]->m_nLength)
					{
						printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write parameter value.\n");
						exit(EXIT_FAILURE);
					}
				}
			}
		}
	}

	/* write sequence data */
	for(nj=0; nj<pBARPos->nSeqNum; nj++)
	{
		/* create BARSeq object */
		if(pBARPos->vSeqData[nj] == NULL)
		{
			printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot find BARSeq object.\n");
			exit(EXIT_FAILURE);
		}

		/* write sequence name */
		if(pBARPos->vSeqData[nj]->pSeqName == NULL)
		{
			nLen = 0;
			if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOut, little_endian_machine) != 1)
			{
				printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write sequence name length.\n");
				exit(EXIT_FAILURE);
			}
			if(nApplyMask == 1)
			{
				if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOutMask, little_endian_machine) != 1)
				{
					printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write sequence name length.\n");
					exit(EXIT_FAILURE);
				}
			}
		}
		else
		{
			if(big_endian_fwrite(&(pBARPos->vSeqData[nj]->pSeqName->m_nLength), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
			{
				printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write sequence name length.\n");
				exit(EXIT_FAILURE);
			}
			
			if(big_endian_fwrite(pBARPos->vSeqData[nj]->pSeqName->m_pString, 1, pBARPos->vSeqData[nj]->pSeqName->m_nLength, fpOut, little_endian_machine) != pBARPos->vSeqData[nj]->pSeqName->m_nLength)
			{
				printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write sequence name.\n");
				exit(EXIT_FAILURE);
			}

			if(nApplyMask == 1)
			{
				if(big_endian_fwrite(&(pBARPos->vSeqData[nj]->pSeqName->m_nLength), INT_SIZE, 1, fpOutMask, little_endian_machine) != 1)
				{
					printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write sequence name length.\n");
					exit(EXIT_FAILURE);
				}
				
				if(big_endian_fwrite(pBARPos->vSeqData[nj]->pSeqName->m_pString, 1, pBARPos->vSeqData[nj]->pSeqName->m_nLength, fpOutMask, little_endian_machine) != pBARPos->vSeqData[nj]->pSeqName->m_nLength)
				{
					printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write sequence name.\n");
					exit(EXIT_FAILURE);
				}
			}
		}
		
		/* write group name */
		if(pBARPos->fVersionnumber > 1.5)
		{
			if(pBARPos->vSeqData[nj]->pSeqGroupName == NULL)
			{
				nLen = 0;
				if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOut, little_endian_machine) != 1)
				{
					printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write sequence group name length.\n");
					exit(EXIT_FAILURE);
				}

				if(nApplyMask == 1)
				{
					if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOutMask, little_endian_machine) != 1)
					{
						printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write sequence group name length.\n");
						exit(EXIT_FAILURE);
					}
				}
			}
			else
			{
				if(big_endian_fwrite(&(pBARPos->vSeqData[nj]->pSeqGroupName->m_nLength), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
				{
					printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write sequence group name length.\n");
					exit(EXIT_FAILURE);
				}

				if(big_endian_fwrite(pBARPos->vSeqData[nj]->pSeqGroupName->m_pString, 1, pBARPos->vSeqData[nj]->pSeqGroupName->m_nLength, fpOut, little_endian_machine) != pBARPos->vSeqData[nj]->pSeqGroupName->m_nLength)
				{
					printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write sequence group name.\n");
					exit(EXIT_FAILURE);
				}

				if(nApplyMask == 1)
				{
					if(big_endian_fwrite(&(pBARPos->vSeqData[nj]->pSeqGroupName->m_nLength), INT_SIZE, 1, fpOutMask, little_endian_machine) != 1)
					{
						printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write sequence group name length.\n");
						exit(EXIT_FAILURE);
					}

					if(big_endian_fwrite(pBARPos->vSeqData[nj]->pSeqGroupName->m_pString, 1, pBARPos->vSeqData[nj]->pSeqGroupName->m_nLength, fpOutMask, little_endian_machine) != pBARPos->vSeqData[nj]->pSeqGroupName->m_nLength)
					{
						printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write sequence group name.\n");
						exit(EXIT_FAILURE);
					}
				}
			}
		}

		/* write sequence version */
		if(pBARPos->vSeqData[nj]->pSeqVersion == NULL)
		{
			nLen = 0;
			if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOut, little_endian_machine) != 1)
			{
				printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write sequence version length.\n");
				exit(EXIT_FAILURE);
			}
			if(nApplyMask == 1)
			{
				if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOutMask, little_endian_machine) != 1)
				{
					printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write sequence version length.\n");
					exit(EXIT_FAILURE);
				}
			}
		}
		else
		{
			if(big_endian_fwrite(&(pBARPos->vSeqData[nj]->pSeqVersion->m_nLength), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
			{
				printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write sequence version length.\n");
				exit(EXIT_FAILURE);
			}
			
			if(big_endian_fwrite(pBARPos->vSeqData[nj]->pSeqVersion->m_pString, 1, pBARPos->vSeqData[nj]->pSeqVersion->m_nLength, fpOut, little_endian_machine) != pBARPos->vSeqData[nj]->pSeqVersion->m_nLength)
			{
				printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write sequence version.\n");
				exit(EXIT_FAILURE);
			}

			if(nApplyMask == 1)
			{
				if(big_endian_fwrite(&(pBARPos->vSeqData[nj]->pSeqVersion->m_nLength), INT_SIZE, 1, fpOutMask, little_endian_machine) != 1)
				{
					printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write sequence version length.\n");
					exit(EXIT_FAILURE);
				}
				
				if(big_endian_fwrite(pBARPos->vSeqData[nj]->pSeqVersion->m_pString, 1, pBARPos->vSeqData[nj]->pSeqVersion->m_nLength, fpOutMask, little_endian_machine) != pBARPos->vSeqData[nj]->pSeqVersion->m_nLength)
				{
					printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write sequence version.\n");
					exit(EXIT_FAILURE);
				}
			}
		}
		
		/* write parameter name/value pairs */
		if(pBARPos->fVersionnumber > 1.5)
		{
			if(big_endian_fwrite(&(pBARPos->vSeqData[nj]->nParamNum), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
			{
				printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write number of name/value pairs for BARseq object.\n");
				exit(EXIT_FAILURE);
			}
			if(nApplyMask == 1)
			{
				if(big_endian_fwrite(&(pBARPos->vSeqData[nj]->nParamNum), INT_SIZE, 1, fpOutMask, little_endian_machine) != 1)
				{
					printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write number of name/value pairs for BARseq object.\n");
					exit(EXIT_FAILURE);
				}
			}

			if(pBARPos->vSeqData[nj]->nParamNum > 0)
			{
				if( (pBARPos->vSeqData[nj]->vParamName == NULL) || (pBARPos->vSeqData[nj]->vParamValue == NULL) )
				{
					printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot find sequence parameter name/value pairs\n");
					exit(EXIT_FAILURE);
				}

				for(nk=0; nk<pBARPos->vSeqData[nj]->nParamNum; nk++)
				{
					if(pBARPos->vSeqData[nj]->vParamName[nk] == NULL)
					{
						nLen = 0;
						if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOut, little_endian_machine) != 1)
						{
							printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write sequence parameter name length.\n");
							exit(EXIT_FAILURE);
						}

						if(nApplyMask == 1)
						{
							if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOutMask, little_endian_machine) != 1)
							{
								printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write sequence parameter name length.\n");
								exit(EXIT_FAILURE);
							}
						}
					}
					else
					{
						if(big_endian_fwrite(&(pBARPos->vSeqData[nj]->vParamName[nk]->m_nLength), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
						{
							printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write sequence parameter name length.\n");
							exit(EXIT_FAILURE);
						}

						if(big_endian_fwrite(pBARPos->vSeqData[nj]->vParamName[nk]->m_pString, 1, pBARPos->vSeqData[nj]->vParamName[nk]->m_nLength, fpOut, little_endian_machine) != pBARPos->vSeqData[nj]->vParamName[nk]->m_nLength)
						{
							printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write sequence parameter name.\n");
							exit(EXIT_FAILURE);
						}

						if(nApplyMask == 1)
						{
							if(big_endian_fwrite(&(pBARPos->vSeqData[nj]->vParamName[nk]->m_nLength), INT_SIZE, 1, fpOutMask, little_endian_machine) != 1)
							{
								printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write sequence parameter name length.\n");
								exit(EXIT_FAILURE);
							}

							if(big_endian_fwrite(pBARPos->vSeqData[nj]->vParamName[nk]->m_pString, 1, pBARPos->vSeqData[nj]->vParamName[nk]->m_nLength, fpOutMask, little_endian_machine) != pBARPos->vSeqData[nj]->vParamName[nk]->m_nLength)
							{
								printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write sequence parameter name.\n");
								exit(EXIT_FAILURE);
							}
						}
					}
					
					if(pBARPos->vSeqData[nj]->vParamValue[nk] == NULL)
					{
						nLen = 0;
						if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOut, little_endian_machine) != 1)
						{
							printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write sequence parameter value length.\n");
							exit(EXIT_FAILURE);
						}

						if(nApplyMask == 1)
						{
							if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOutMask, little_endian_machine) != 1)
							{
								printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write sequence parameter value length.\n");
								exit(EXIT_FAILURE);
							}
						}
					}
					else
					{
						if(big_endian_fwrite(&(pBARPos->vSeqData[nj]->vParamValue[nk]->m_nLength), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
						{
							printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write sequence parameter value length.\n");
							exit(EXIT_FAILURE);
						}

						if(big_endian_fwrite(pBARPos->vSeqData[nj]->vParamValue[nk]->m_pString, 1, pBARPos->vSeqData[nj]->vParamValue[nk]->m_nLength, fpOut, little_endian_machine) != pBARPos->vSeqData[nj]->vParamValue[nk]->m_nLength)
						{
							printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write parameter value.\n");
							exit(EXIT_FAILURE);
						}

						if(nApplyMask == 1)
						{
							if(big_endian_fwrite(&(pBARPos->vSeqData[nj]->vParamValue[nk]->m_nLength), INT_SIZE, 1, fpOutMask, little_endian_machine) != 1)
							{
								printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write sequence parameter value length.\n");
								exit(EXIT_FAILURE);
							}

							if(big_endian_fwrite(pBARPos->vSeqData[nj]->vParamValue[nk]->m_pString, 1, pBARPos->vSeqData[nj]->vParamValue[nk]->m_nLength, fpOutMask, little_endian_machine) != pBARPos->vSeqData[nj]->vParamValue[nk]->m_nLength)
							{
								printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write parameter value.\n");
								exit(EXIT_FAILURE);
							}
						}
					}
				}
			}
		}

		/* write data points */
		if(big_endian_fwrite(&(pBARPos->vSeqData[nj]->nDataNum), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
		{
			printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write number of data points for a sequence.\n");
			exit(EXIT_FAILURE);
		}

		if(nApplyMask == 1)
		{
			if(big_endian_fwrite(&(pBARPos->vSeqData[nj]->nDataNum), INT_SIZE, 1, fpOutMask, little_endian_machine) != 1)
			{
				printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot write number of data points for a sequence.\n");
				exit(EXIT_FAILURE);
			}
		}

		if(pBARPos->vSeqData[nj]->nColNum > 0)
		{
			if(pBARPos->vSeqData[nj]->vData == NULL)
			{
				printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot find sequence data.\n");
				exit(EXIT_FAILURE);
			}
		}

		if(pBARPos->vSeqData[nj]->nDataNum > 0)
		{
			for(nk=0; nk<pBARPos->vSeqData[nj]->nColNum; nk++)
			{
				if(pBARPos->vSeqData[nj]->vData[nk] == NULL)
				{
					printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot find sequence data.\n");
					exit(EXIT_FAILURE);
				}
			}

			for(nk=0; nk<pBARPos->vSeqData[nj]->nDataNum; nk++)
			{
				nidx1 = (int)(pBARPos->vSeqData[nj]->vData[2]->pMatElement[nk])*nCELTotalX+(int)(pBARPos->vSeqData[nj]->vData[1]->pMatElement[nk]);
				if(nidx1 >= nCELNumberCells)
				{
					printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, index out of range\n");
					exit(EXIT_FAILURE);
				}

				if(pBARPos->vSeqData[nj]->nColNum >= 5)
				{
					nidx2 = (int)(pBARPos->vSeqData[nj]->vData[4]->pMatElement[nk])*nCELTotalX+(int)(pBARPos->vSeqData[nj]->vData[3]->pMatElement[nk]);
					if(nidx2 >= nCELNumberCells)
					{
						printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, index out of range\n");
						exit(EXIT_FAILURE);
					}
				}

				nV = (int)(pBARPos->vSeqData[nj]->vData[0]->pMatElement[nk]);
				big_endian_fwrite(&nV, INT_SIZE, 1, fpOut, little_endian_machine);
				if(nApplyMask == 1)
				{
					big_endian_fwrite(&nV, INT_SIZE, 1, fpOutMask, little_endian_machine);
				}
								
				for(ni=0; ni<nSampleNum; ni++)
				{
					if( fseek( vfpIn[ni], nidx1*sizeof(double), SEEK_SET ) != 0 )
					{
						printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot locate the required data point!\n");
						exit(EXIT_FAILURE);
					}

					if( little_endian_fread(&dPM, sizeof(double), 1, vfpIn[ni], little_endian_machine) != 1)
					{
						printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot load data from binary files correctly!\n");
						exit(EXIT_FAILURE);
					}

					if(pBARPos->vSeqData[nj]->nColNum >= 5)
					{
						if( fseek( vfpIn[ni], nidx2*sizeof(double), SEEK_SET ) != 0 )
						{
							printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot locate the required data point!\n");
							exit(EXIT_FAILURE);
						}

						if( little_endian_fread(&dMM, sizeof(double), 1, vfpIn[ni], little_endian_machine) != 1)
						{
							printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot load data from binary files correctly!\n");
							exit(EXIT_FAILURE);
						}
					}
					else
					{
						dMM = 0.0;
					}

					if(nIntensityType == 1)
					{
						dPM -= dMM;
					}
						
					if(dPM < dIntLowerBound)
						dPM = dIntLowerBound;

					if(nIntLogTransform == 1)
					{
						dPM = log(dPM)/dLog2;
					}

					fV = (float)dPM;
					big_endian_fwrite(&fV, FLOAT_SIZE, 1, fpOut, little_endian_machine);

					/* mask */
					if(nApplyMask == 1)
					{
						if( fseek( vfpMask[ni], nidx1, SEEK_SET ) != 0 )
						{
							printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot locate the required data point!\n");
							exit(EXIT_FAILURE);
						}

						if( little_endian_fread(&bMask, sizeof(unsigned char), 1, vfpMask[ni], little_endian_machine) != 1)
						{
							printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot load data from binary files correctly!\n");
							exit(EXIT_FAILURE);
						}

						if(bMask == 1)
							nMask = 1;
						else
							nMask = 0;

						if(pBARPos->vSeqData[nj]->nColNum >= 5)
						{
							if( fseek( vfpMask[ni], nidx2, SEEK_SET ) != 0 )
							{
								printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot locate the required data point!\n");
								exit(EXIT_FAILURE);
							}

							if( little_endian_fread(&bMask, sizeof(unsigned char), 1, vfpMask[ni], little_endian_machine) != 1)
							{
								printf("Error: TileMapv2_ExportAffyIntensity_SingleBAR, cannot load data from binary files correctly!\n");
								exit(EXIT_FAILURE);
							}

							if( (bMask == 1) && (nIntensityType == 1) )
								nMask = 1;
						}

						if(nMask == 1)
							bMask = 1;
						else
							bMask = 0;

						big_endian_fwrite(&bMask, 1, 1, fpOutMask, little_endian_machine);
					}

				}
			}
		}
	}

	/* close files */
	fclose(fpOut);
	if(nApplyMask == 1)
	{
		fclose(fpOutMask);
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ExportAffyIntensity_SingleTXT()                              */
/*  Export intensity data to a single txt files.                           */
/* ----------------------------------------------------------------------- */
int TileMapv2_ExportAffyIntensity_SingleTXT(char strExportPath[], 
			int nSampleNum, struct tagString **vSampleAlias, FILE **vfpIn,
			int nCELNumberCells, int nCELTotalX, int nCELTotalY,
			int nApplyMask, FILE **vfpMask,	struct tagBARData *pBARPos, 
			int nIntensityType,	double dIntLowerBound, int nIntLogTransform)
{
	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));
	FILE *fpOut;
	FILE *fpOutMask;
	char strLine[MED_LINE_LENGTH];
	char strSeqAlias[MED_LINE_LENGTH];
	int ni,nj,nk;
	float fV;
	int nV;
	int nidx1,nidx2;
	double dLog2 = log(2.0);
	int nMask;
	double dPM,dMM;
	unsigned char bMask;
	
	/* init */
	if( (nSampleNum <= 0) || (nCELNumberCells <= 0))
	{
		return PROC_SUCCESS;
	}

	/* open files */
	fpOut = NULL;
	fpOutMask = NULL;

	fpOut = fopen(strExportPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: TileMapv2_ExportAffyIntensity_SingleTXT, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	
	if(nApplyMask == 1)
	{
		sprintf(strLine, "%s.mask", strExportPath);
		fpOutMask = fopen(strLine, "w");
		if(fpOutMask == NULL)
		{
			printf("Error: TileMapv2_ExportAffyIntensity_SingleTXT, cannot open output file!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	/* write sequence data */
	for(nj=0; nj<pBARPos->nSeqNum; nj++)
	{
		/* create BARSeq object */
		if(pBARPos->vSeqData[nj] == NULL)
		{
			printf("Error: TileMapv2_ExportAffyIntensity_SingleTXT, cannot find BARSeq object.\n");
			exit(EXIT_FAILURE);
		}

		/* write sequence name */
		if(pBARPos->vSeqData[nj]->pSeqName == NULL)
		{
			printf("Error: TileMapv2_ExportAffyIntensity_SingleTXT, cannot find sequence name.\n");
			exit(EXIT_FAILURE);
		}
		if(pBARPos->vSeqData[nj]->pSeqVersion != NULL)
		{
			sprintf(strSeqAlias, "%s:%s", pBARPos->vSeqData[nj]->pSeqVersion->m_pString,
				pBARPos->vSeqData[nj]->pSeqName->m_pString);
		}
		else
		{
			sprintf(strSeqAlias, "%s", pBARPos->vSeqData[nj]->pSeqName->m_pString);
		}

		if(pBARPos->fVersionnumber > 1.5)
		{
			if(pBARPos->vSeqData[nj]->pSeqGroupName != NULL)
			{
				if(pBARPos->vSeqData[nj]->pSeqVersion != NULL)
				{
					sprintf(strSeqAlias, "%s:%s:%s", pBARPos->vSeqData[nj]->pSeqGroupName->m_pString,
						pBARPos->vSeqData[nj]->pSeqVersion->m_pString,
						pBARPos->vSeqData[nj]->pSeqName->m_pString);
				}
				else
				{
					sprintf(strSeqAlias, "%s:%s", pBARPos->vSeqData[nj]->pSeqGroupName->m_pString,
						pBARPos->vSeqData[nj]->pSeqName->m_pString);
				}
			}	
		}

		if(pBARPos->vSeqData[nj]->nDataNum > 0)
		{
			for(nk=0; nk<pBARPos->vSeqData[nj]->nColNum; nk++)
			{
				if(pBARPos->vSeqData[nj]->vData[nk] == NULL)
				{
					printf("Error: TileMapv2_ExportAffyIntensity_SingleTXT, cannot find sequence data.\n");
					exit(EXIT_FAILURE);
				}
			}

			for(nk=0; nk<pBARPos->vSeqData[nj]->nDataNum; nk++)
			{
				nidx1 = (int)(pBARPos->vSeqData[nj]->vData[2]->pMatElement[nk])*nCELTotalX+(int)(pBARPos->vSeqData[nj]->vData[1]->pMatElement[nk]);
				if(nidx1 >= nCELNumberCells)
				{
					printf("Error: TileMapv2_ExportAffyIntensity_SingleTXT, index out of range\n");
					exit(EXIT_FAILURE);
				}

				if(pBARPos->vSeqData[nj]->nColNum >= 5)
				{
					nidx2 = (int)(pBARPos->vSeqData[nj]->vData[4]->pMatElement[nk])*nCELTotalX+(int)(pBARPos->vSeqData[nj]->vData[3]->pMatElement[nk]);
					if(nidx2 >= nCELNumberCells)
					{
						printf("Error: TileMapv2_ExportAffyIntensity_SingleTXT, index out of range\n");
						exit(EXIT_FAILURE);
					}
				}

				nV = (int)(pBARPos->vSeqData[nj]->vData[0]->pMatElement[nk]);
				fprintf(fpOut, "%s\t%d", strSeqAlias, nV);
				if(nApplyMask == 1)
				{
					fprintf(fpOutMask, "%s\t%d", strSeqAlias, nV);
				}
								
				for(ni=0; ni<nSampleNum; ni++)
				{
					if( fseek( vfpIn[ni], nidx1*sizeof(double), SEEK_SET ) != 0 )
					{
						printf("Error: TileMapv2_ExportAffyIntensity_SingleTXT, cannot locate the required data point!\n");
						exit(EXIT_FAILURE);
					}

					if( little_endian_fread(&dPM, sizeof(double), 1, vfpIn[ni], little_endian_machine) != 1)
					{
						printf("Error: TileMapv2_ExportAffyIntensity_SingleTXT, cannot load data from binary files correctly!\n");
						exit(EXIT_FAILURE);
					}

					if(pBARPos->vSeqData[nj]->nColNum >= 5)
					{
						if( fseek( vfpIn[ni], nidx2*sizeof(double), SEEK_SET ) != 0 )
						{
							printf("Error: TileMapv2_ExportAffyIntensity_SingleTXT, cannot locate the required data point!\n");
							exit(EXIT_FAILURE);
						}

						if( little_endian_fread(&dMM, sizeof(double), 1, vfpIn[ni], little_endian_machine) != 1)
						{
							printf("Error: TileMapv2_ExportAffyIntensity_SingleTXT, cannot load data from binary files correctly!\n");
							exit(EXIT_FAILURE);
						}
					}
					else
					{
						dMM = 0.0;
					}

					if(nIntensityType == 1)
					{
						dPM -= dMM;
					}
						
					if(dPM < dIntLowerBound)
						dPM = dIntLowerBound;

					if(nIntLogTransform == 1)
					{
						dPM = log(dPM)/dLog2;
					}

					fV = (float)dPM;
					fprintf(fpOut, "\t%f", fV);
					
					/* mask */
					if(nApplyMask == 1)
					{
						if( fseek( vfpMask[ni], nidx1, SEEK_SET ) != 0 )
						{
							printf("Error: TileMapv2_ExportAffyIntensity_SingleTXT, cannot locate the required data point!\n");
							exit(EXIT_FAILURE);
						}

						if( little_endian_fread(&bMask, sizeof(unsigned char), 1, vfpMask[ni], little_endian_machine) != 1)
						{
							printf("Error: TileMapv2_ExportAffyIntensity_SingleTXT, cannot load data from binary files correctly!\n");
							exit(EXIT_FAILURE);
						}

						if(bMask == 1)
							nMask = 1;
						else
							nMask = 0;

						if(pBARPos->vSeqData[nj]->nColNum >= 5)
						{
							if( fseek( vfpMask[ni], nidx2, SEEK_SET ) != 0 )
							{
								printf("Error: TileMapv2_ExportAffyIntensity_SingleTXT, cannot locate the required data point!\n");
								exit(EXIT_FAILURE);
							}

							if( little_endian_fread(&bMask, sizeof(unsigned char), 1, vfpMask[ni], little_endian_machine) != 1)
							{
								printf("Error: TileMapv2_ExportAffyIntensity_SingleTXT, cannot load data from binary files correctly!\n");
								exit(EXIT_FAILURE);
							}

							if( (bMask == 1) && (nIntensityType == 1) )
								nMask = 1;
						}

						fprintf(fpOutMask, "\t%d", nMask);
					}
				}

				fprintf(fpOut, "\n");
				if(nApplyMask == 1)
				{
					fprintf(fpOutMask, "\n");
				}
			}
		}
	}

	/* close files */
	fclose(fpOut);
	if(nApplyMask == 1)
	{
		fclose(fpOutMask);
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2Param_Create()                                                */
/*  Create TileMapv2Param object.                                          */
/* ----------------------------------------------------------------------- */ 
struct tagTileMapv2Param *TileMapv2Param_Create()
{
	/* define */
	struct tagTileMapv2Param *pParam = NULL;

	/* create */
	pParam = (struct tagTileMapv2Param *)calloc(1, sizeof(struct tagTileMapv2Param));
	if(pParam == NULL)
	{
		printf("Error: TileMapv2Param_Create, cannot create tilemap parameters!\n");
		exit(EXIT_FAILURE);
	}

	pParam->nComparisonType = 2;

	pParam->nNoiseMask = 0;
	pParam->dLower = -1e6;
	pParam->dUpper = 1e6;
	pParam->nTransformType = 0;

	pParam->nMonteCarloNum = 1000;
	
	pParam->nW = 5;
	pParam->nWindowBoundary = 200;
	pParam->nMAStandardize = 1;
	pParam->dBaseStd = 3.0;
	
	pParam->nExpLen = 28;
	pParam->dPostCut = 0.5;
	pParam->nHMMParamUserSpecified = 0;
	strcpy(pParam->strTransitionPath, "");
	strcpy(pParam->strEmissionPath, "");
	
	pParam->dTp = 0.01;
	pParam->dTq = 0.05;
	pParam->nOffset = 1;
	pParam->nGridSize = 1000;

	pParam->nGap = 400;
	pParam->nGapW = 5;
	pParam->nMinRegLen = 1;
	pParam->nMinRegProbeNum = 1;

	/* return */
	return pParam;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2Param_Destroy()                                               */
/*  Create TileMapv2Param object.                                          */
/* ----------------------------------------------------------------------- */ 
int TileMapv2Param_Destroy(struct tagTileMapv2Param **pParam)
{
	/* define */
	int ni;

	/* init */
	if(pParam == NULL)
		return PROC_SUCCESS;

	if(*pParam == NULL)
		return PROC_SUCCESS;

	/* Sample name */
	DestroyIntMatrix((*pParam)->pGroupLabel);
	DestroyIntMatrix((*pParam)->pGroupSize);
	for(ni=0; ni<(*pParam)->nSampleNum; ni++)
	{
		DeleteString((*pParam)->vSampleAlias[ni]);
	}
	free((*pParam)->vSampleAlias);
	for(ni=0; ni<(*pParam)->nArrayNum; ni++)
	{
		DeleteString((*pParam)->vArrayFile[ni]);
	}
	free((*pParam)->vArrayFile);
	
	/* Probe summary */
	for(ni=0; ni<(*pParam)->nVargroupNum; ni++)
	{
		DestroyIntMatrix((*pParam)->vVargroupMap[ni]);
	}
	free((*pParam)->vVargroupMap);
	DestroyIntMatrix((*pParam)->pVargroupSize);

	/* permutation parameters */
	for(ni=0; ni<(*pParam)->nPermgroupNum; ni++)
	{
		DestroyIntMatrix((*pParam)->vPermgroupMap[ni]);
	}
	free((*pParam)->vPermgroupMap);
	DestroyIntMatrix((*pParam)->pPermgroupSize);

	/* free the whole structure */
	free(*pParam);
	*pParam = NULL;
	
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_Main()                                                       */
/*  TileMapv2 pipeline                                                     */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_Main(char strParamPath[])
{
	/* ------ */
	/* define */
	/* ------ */
	struct tagTileMapv2Param *pParam = NULL;

	/* --------------- */
	/* load parameters */
	/* --------------- */
	pParam = TileMapv2_LoadParameters(strParamPath);
	if(pParam == NULL)
	{
		printf("Error: TileMapv2_Main, cannot load the parameter file correctly!\n");
		exit(EXIT_FAILURE);
	}
	
	if(pParam->nComparisonType == 1)
	{
		if(pParam->nPatternType != 2)
		{
			printf("Error: TileMapv2_Main, parameter file wrongly specified!\n");
			exit(EXIT_FAILURE);
		}
		pParam->nPatternType = 1;
	}
	else if(pParam->nComparisonType == 2)
	{
		if(pParam->nPatternType != 2)
		{
			printf("Error: TileMapv2_Main, parameter file wrongly specified!\n");
			exit(EXIT_FAILURE);
		}	
	}
	else
	{
		
	}
	
	TileMapv2(pParam);

	/* -------------------- */
	/* Done                 */
	/* -------------------- */
	TileMapv2Param_Destroy(&pParam);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_LoadParameters()                                             */
/*  Load tilemap parameters.                                               */
/* ----------------------------------------------------------------------- */ 
struct tagTileMapv2Param *TileMapv2_LoadParameters(char strParamFile[])
{
	/* define */
	struct tagTileMapv2Param *pParam = NULL;
	FILE *fpIn = NULL;
	char strLine[MED_LINE_LENGTH];
	int ni,nj,nk,nz;
	char *chSep;

	/* init */
	fpIn = fopen(strParamFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileMapv2_LoadParameters, cannot open the parameter file!\n");
		exit(EXIT_FAILURE);
	}

	/* create */
	pParam = TileMapv2Param_Create();
	if(pParam == NULL)
	{
		printf("Error: TileMapv2_LoadParameters, cannot load tilemap parameters!\n");
		exit(EXIT_FAILURE);
	}

	/* load */
	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strstr(strLine, "[Comparison Type]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_LoadParameters, cannot load tilemap parameters!\n");
				printf("%s\n", strLine);
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			if(*chSep != '\0')
				pParam->nComparisonType = atoi(chSep);
			else
			{
				printf("Error: TileMapv2_LoadParameters, unrecognized comparison type!\n");
				exit(EXIT_FAILURE);
			}

			if( (pParam->nComparisonType <= 0) || (pParam->nComparisonType > 3) )
			{
				printf("Error: TileMapv2_LoadParameters, unrecognized comparison type!\n");
				exit(EXIT_FAILURE);
			}

		}

		else if(strstr(strLine, "[Working Directory]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_LoadParameters, cannot load tilemap parameters!\n");
				printf("%s\n", strLine);
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			if(*chSep != '\0')
				strcpy(pParam->strWorkPath, chSep);
			else
				strcpy(pParam->strWorkPath, ".");
			AdjustDirectoryPath(pParam->strWorkPath);
		}

		else if(strstr(strLine, "[Project Title]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_LoadParameters, cannot load tilemap parameters!\n");
				printf("%s\n", strLine);
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			if(*chSep == '\0')
			{
				strcpy(pParam->strProjectTitle, "unnamed_project");
			}
			else
			{
				strcpy(pParam->strProjectTitle, chSep);
			}
		}

		else if(strstr(strLine, "[No. of Libraries]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_LoadParameters, cannot load tilemap parameters!\n");
				printf("%s\n", strLine);
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			if(*chSep != '\0')
				pParam->nLibNum = atoi(chSep);
			pParam->nArrayNum = (pParam->nLibNum)*(pParam->nSampleNum);
		}

		else if(strstr(strLine, "[No. of Samples]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_LoadParameters, cannot load tilemap parameters!\n");
				printf("%s\n", strLine);
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
				pParam->nSampleNum = atoi(chSep);
			pParam->nArrayNum = (pParam->nLibNum)*(pParam->nSampleNum);
		}

		else if(strstr(strLine, "[No. of Groups]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_LoadParameters, cannot load tilemap parameters!\n");
				printf("%s\n", strLine);
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);

			if(*chSep != '\0')
				pParam->nGroupNum = atoi(chSep);
		}

		else if(strstr(strLine, "[Data]") == strLine)
		{
			if( (pParam->nSampleNum>0) && (pParam->nLibNum>0) && (pParam->nGroupNum>0))
			{
				pParam->pGroupLabel = NULL;
				pParam->pGroupLabel = CreateIntMatrix(1,pParam->nSampleNum);
				if(pParam->pGroupLabel == NULL)
				{
					printf("Error: TileMapv2_LoadParameters, cannot allocate memory for group labels!\n");
					exit(EXIT_FAILURE);
				}
				pParam->pGroupSize = NULL;
				pParam->pGroupSize = CreateIntMatrix(1,pParam->nGroupNum);
				if(pParam->pGroupSize == NULL)
				{
					printf("Error: TileMapv2_LoadParameters, cannot allocate memory for group labels!\n");
					exit(EXIT_FAILURE);
				}
				pParam->vSampleAlias = NULL;
				pParam->vSampleAlias = (struct tagString **)calloc(pParam->nSampleNum, sizeof(struct tagString *));
				if(pParam->vSampleAlias == NULL)
				{
					printf("Error: TileMapv2_LoadParameters, cannot allocate memory for sample alias!\n");
					exit(EXIT_FAILURE);
				}
				pParam->vArrayFile = NULL;
				pParam->vArrayFile = (struct tagString **)calloc(pParam->nArrayNum, sizeof(struct tagString *));
				if(pParam->vArrayFile == NULL)
				{
					printf("Error: TileMapv2_LoadParameters, cannot allocate memory for array file names!\n");
					exit(EXIT_FAILURE);
				}

				ni = 0;
				nj = 0;
				nk = 0;

				while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
				{
					StrTrimLeft(strLine);
					StrTrimRight(strLine);
					if(strLine[0] == '\0')
						continue;

					if(nj == 0)
					{
						chSep = strstr(strLine, "->");
						if(chSep == NULL)
						{
							printf("Error: TileMapv2_LoadParameters, sample name format error!\n");
							exit(EXIT_FAILURE);
						}

						*chSep = '\0';
						chSep = chSep+2;

						StrTrimLeft(strLine);
						StrTrimRight(strLine);
						pParam->pGroupLabel->pMatElement[ni] = atoi(strLine);
						if( pParam->pGroupLabel->pMatElement[ni] > pParam->nGroupNum )
						{
							printf("Error: TileMapv2_LoadParameters, Group ID out of range!\n");
							exit(EXIT_FAILURE);
						}
						if( pParam->pGroupLabel->pMatElement[ni] > 0)
						{
							nz = pParam->pGroupLabel->pMatElement[ni]-1;
							pParam->pGroupSize->pMatElement[nz] +=1; 
						}
					
						StrTrimLeft(chSep);
						StrTrimRight(chSep);
						StringAddTail(pParam->vSampleAlias+ni, chSep);
					}
					else
					{
						StringAddTail((pParam->vArrayFile+(ni*pParam->nLibNum+nj-1)), strLine);
					}

					nj++;
					if(nj > pParam->nLibNum)
					{
						ni++;
						nj = 0;
					}
					nk++;
					if(nk == (pParam->nArrayNum+pParam->nSampleNum))
						break;
				}

				if( nk != (pParam->nArrayNum+pParam->nSampleNum) )
				{
					printf("Error: TileMapv2_LoadParameters, array number is not correct!\n");
					printf("Check if LibNum = %d; SampleNum = %d? \n", pParam->nLibNum, pParam->nSampleNum);
					exit(EXIT_FAILURE);
				}
			}

		}

		else if(strstr(strLine, "[Patterns of Interest]") == strLine)
		{
			if(fgets(strLine, MED_LINE_LENGTH, fpIn) == NULL)
			{
				printf("Error: TileMapv2_LoadParameters, cannot load tilemap parameters!\n");
				printf("%s\n", strLine);
				exit(EXIT_FAILURE);
			}
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			
			if(strLine[0] == '\0')
			{
				printf("Error: TileMapv2_LoadParameters, cannot load tilemap parameters!\n");
				printf("%s\n", strLine);
				exit(EXIT_FAILURE);
			}
			
			strcpy(pParam->strPattern, strLine);

			pParam->nPatternType = 1;
			chSep = strpbrk(pParam->strPattern, "<>" );
			while(chSep != NULL)
			{
				pParam->nPatternType++;
				chSep = strpbrk((chSep+1), "<>" );
			}
		}

		else if(strstr(strLine, "[Masking Bad Data Points]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_LoadParameters, cannot load tilemap parameters!\n");
				printf("%s\n", strLine);
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);

			if(*chSep != '\0')
				pParam->nNoiseMask = atoi(chSep);
		}

		else if(strstr(strLine, "[Truncation Lower Bound]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_LoadParameters, cannot load tilemap parameters!\n");
				printf("%s\n", strLine);
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);

			if(*chSep != '\0')
				pParam->dLower = atof(chSep);
		}

		else if(strstr(strLine, "[Truncation Upper Bound]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_LoadParameters, cannot load tilemap parameters!\n");
				printf("%s\n", strLine);
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);

			if(*chSep != '\0')
				pParam->dUpper = atof(chSep);
		}

		else if(strstr(strLine, "[Transformation]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_LoadParameters, cannot load tilemap parameters!\n");
				printf("%s\n", strLine);
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);

			if(*chSep != '\0')
				pParam->nTransformType = atoi(chSep);
		}

		else if(strstr(strLine, "[Number of Monte Carlo Draws]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_LoadParameters, cannot load tilemap parameters!\n");
				printf("%s\n", strLine);
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);

			if(*chSep != '\0')
				pParam->nMonteCarloNum = atoi(chSep);
		}

		else if(strstr(strLine, "[Common Variance Groups]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_LoadParameters, cannot load tilemap parameters!\n");
				printf("%s\n", strLine);
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
				pParam->nVargroupNum = atoi(chSep);

			if(pParam->nVargroupNum <= 0)
			{
				printf("Error: TileMapv2_LoadParameters, there must be at least one common variance group!\n");
				exit(EXIT_FAILURE);
			}

			/* groups */
			pParam->pVargroupSize = NULL;
			pParam->pVargroupSize = CreateIntMatrix(pParam->nVargroupNum, 1);
			if(pParam->pVargroupSize == NULL)
			{
				printf("Error: TileMapv2_LoadParameters, cannot track common variance group!\n");
				exit(EXIT_FAILURE);
			}
			pParam->vVargroupMap = NULL;
			pParam->vVargroupMap = (struct INTMATRIX **)calloc(pParam->nVargroupNum, sizeof(struct INTMATRIX *));
			if(pParam->vVargroupMap == NULL)
			{
				printf("Error: TileMapv2_LoadParameters, cannot track common variance group!\n");
				exit(EXIT_FAILURE);
			}

			for(ni=0; ni<pParam->nVargroupNum; ni++)
			{
				fgets(strLine, MED_LINE_LENGTH, fpIn);
				pParam->vVargroupMap[ni] = Expression_GeneSelection_LoadGroupId(strLine);
				if(pParam->vVargroupMap[ni] == NULL)
				{
					printf("Error: TileMapv2_LoadParameters, cannot track common variance group!\n");
					exit(EXIT_FAILURE);
				}
				IMSETAT(pParam->pVargroupSize, ni, 0, pParam->vVargroupMap[ni]->nWidth);
			}
		}

		else if(strstr(strLine, "[Method to Combine Neighboring Probes]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_LoadParameters, cannot load tilemap parameters!\n");
				printf("%s\n", strLine);
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
				pParam->nRegionSummaryType = atoi(chSep);
		}

		else if(strstr(strLine, "[W]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_LoadParameters, cannot load tilemap parameters!\n");
				printf("%s\n", strLine);
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
				pParam->nW = atoi(chSep);
		}

		
		else if(strstr(strLine, "[Window Boundary]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_LoadParameters, cannot load tilemap parameters!\n");
				printf("%s\n", strLine);
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
				pParam->nWindowBoundary = atoi(chSep);
		}

		else if(strstr(strLine, "[Standardize MA Statistics]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_LoadParameters, cannot load tilemap parameters!\n");
				printf("%s\n", strLine);
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
				pParam->nMAStandardize = atoi(chSep);
		}

		else if(strstr(strLine, "[Region Boundary Cutoff, MA>]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_LoadParameters, cannot load tilemap parameters!\n");
				printf("%s\n", strLine);
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
				pParam->dBaseStd = atof(chSep);
		}
		
		else if(strstr(strLine, "[Method to Compute FDR]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_LoadParameters, cannot load tilemap parameters!\n");
				printf("%s\n", strLine);
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
				pParam->nFDRType = atoi(chSep);
		}

		else if(strstr(strLine, "[Expected Hybridization Length]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_LoadParameters, cannot load tilemap parameters!\n");
				printf("%s\n", strLine);
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
				pParam->nExpLen = atoi(chSep);
		}

		else if(strstr(strLine, "[Posterior Probability Cutoff, P>]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_LoadParameters, cannot load tilemap parameters!\n");
				printf("%s\n", strLine);
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
				pParam->dPostCut = atof(chSep);
		}

		else if(strstr(strLine, "[G0 Selection Criteria, p%]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_LoadParameters, cannot load tilemap parameters!\n");
				printf("%s\n", strLine);
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
				pParam->dTp = atof(chSep);
		}

		else if(strstr(strLine, "[G1 Selection Criteria, q%]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_LoadParameters, cannot load tilemap parameters!\n");
				printf("%s\n", strLine);
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
				pParam->dTq = atof(chSep);
		}

		else if(strstr(strLine, "[Selection Offset]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_LoadParameters, cannot load tilemap parameters!\n");
				printf("%s\n", strLine);
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
				pParam->nOffset = atoi(chSep);
		}

		else if(strstr(strLine, "[Grid Size]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_LoadParameters, cannot load tilemap parameters!\n");
				printf("%s\n", strLine);
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
				pParam->nGridSize = atoi(chSep);
		}

		else if(strstr(strLine, "[Number of Permutations]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_LoadParameters, cannot load tilemap parameters!\n");
				printf("%s\n", strLine);
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
				pParam->nPermutationNum = atoi(chSep);
		}

		else if(strstr(strLine, "[Exchangeable Groups]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_LoadParameters, cannot load tilemap parameters!\n");
				printf("%s\n", strLine);
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
				pParam->nPermgroupNum = atoi(chSep);

			/* groups */
			if(pParam->nPermgroupNum > 0)
			{
				pParam->pPermgroupSize = NULL;
				pParam->pPermgroupSize = CreateIntMatrix(pParam->nPermgroupNum, 1);
				if(pParam->pPermgroupSize == NULL)
				{
					printf("Error: TileMapv2_LoadParameters, cannot track permutation group!\n");
					exit(EXIT_FAILURE);
				}
				pParam->vPermgroupMap = NULL;
				pParam->vPermgroupMap = (struct INTMATRIX **)calloc(pParam->nPermgroupNum, sizeof(struct INTMATRIX *));
				if(pParam->vPermgroupMap == NULL)
				{
					printf("Error: TileMapv2_LoadParameters, cannot track permutation group!\n");
					exit(EXIT_FAILURE);
				}

				for(ni=0; ni<pParam->nPermgroupNum; ni++)
				{
					fgets(strLine, MED_LINE_LENGTH, fpIn);
					pParam->vPermgroupMap[ni] = Expression_GeneSelection_LoadGroupId(strLine);
					if(pParam->vPermgroupMap[ni] == NULL)
					{
						printf("Error: TileMapv2_LoadParameters, cannot track permutation group!\n");
						exit(EXIT_FAILURE);
					}
					IMSETAT(pParam->pPermgroupSize, ni, 0, pParam->vPermgroupMap[ni]->nWidth);
				}
			}
		}

		else if(strstr(strLine, "[Max Gap within a Region]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_LoadParameters, cannot load tilemap parameters!\n");
				printf("%s\n", strLine);
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
				pParam->nGap = atoi(chSep);
		}

		else if(strstr(strLine, "[Max Run of Insignificant Probes within a Region]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_LoadParameters, cannot load tilemap parameters!\n");
				printf("%s\n", strLine);
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
				pParam->nGapW = atoi(chSep);
		}

		else if(strstr(strLine, "[Min Region Length]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_LoadParameters, cannot load tilemap parameters!\n");
				printf("%s\n", strLine);
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
				pParam->nMinRegLen = atoi(chSep);
		}

		else if(strstr(strLine, "[Min No. of Significant Probes within a Region]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileMapv2_LoadParameters, cannot load tilemap parameters!\n");
				printf("%s\n", strLine);
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
				pParam->nMinRegProbeNum = atoi(chSep);
		}

		else if(strLine[0] == '#')
		{
		}

		else
		{
			printf("Error: TileMapv2_LoadParameters, unknown parameters!\n");
			exit(EXIT_FAILURE);
		}
	}


	/* close files */
	fclose(fpIn);

	/* return */
	return pParam;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2()                                                            */
/*  TileMap for two sample comparisons.                                    */
/* ----------------------------------------------------------------------- */ 
int TileMapv2(struct tagTileMapv2Param *pParam)
{
	/* define */
	int nLibId;
	int nRegionNum = 0;

	/* init */
	if(pParam == NULL)
	{
		return PROC_SUCCESS;
	}

	/* get probe number */
	printf("############################\n");
	printf("       TileMap_v2.0         \n");
	printf("############################\n");
	
	for(nLibId=0; nLibId<pParam->nLibNum; nLibId++)
	{
		printf("\n");
		printf("############################\n");
		printf(" Processing Arrayset %d ... \n", (nLibId+1));
		printf("############################\n");

		printf("\n");
		printf("/* ---------------------- */\n");
		printf("/* Probe Level Summary    */\n");
		printf("/* ---------------------- */\n");
		TileMapv2_ProbeSelection_Main(pParam, nLibId);

		if(pParam->nRegionSummaryType == 1)
		{
			printf("\n");
			printf("/* ---------------------- */\n");
			printf("/* TileMap MA             */\n");
			printf("/* ---------------------- */\n");
			nRegionNum += TileMapv2_RegionDetection_MA_Main(pParam, nLibId);
		}
		else
		{
			printf("\n");
			printf("/* ---------------------- */\n");
			printf("/* TileMap HMM            */\n");
			printf("/* ---------------------- */\n");
			nRegionNum += TileMapv2_RegionDetection_HMM_Main(pParam, nLibId);
		}
	}

	/* -------------------- */
	/* combine results      */
	/* -------------------- */
	printf("\n");
	printf("############################\n");
	printf(" Summarizing Results ...    \n");
	printf("############################\n");
	TileMapv2_SummarizeResults_Main(pParam, nRegionNum);
	
	/* -------------------- */
	/* Done                 */
	/* -------------------- */
	printf("Done!\n");

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ProbeSelection_Main()                                        */
/*  TileMapv2 probe level summary.                                         */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_ProbeSelection_Main(struct tagTileMapv2Param *pParam, 
						int nLibId)
{
	/* define */
	char strFileName[MED_LINE_LENGTH];
	char strFCFileName[MED_LINE_LENGTH];
	struct tagBARData *pData = NULL;
	struct tagBARData *pDataNew = NULL;
	struct tagBARData *pMask = NULL;
	struct tagBARData *pMaskNew = NULL;
	struct INTMATRIX *pFieldType = NULL;
	struct DOUBLEMATRIX **vDataVec = NULL;
	int nTotalProbeNum = 0;
	int ni,nj;
	struct INTMATRIX *pPermGroupLabel = NULL;
	struct INTMATRIX *pOriGroupLabel = NULL;

	/* init */
	if(pParam == NULL)
	{
		return PROC_SUCCESS;
	}

	/* ------------- */
	/* load raw data */
	/* ------------- */
	for(ni=0; ni<pParam->nSampleNum; ni++)
	{
		/* load individual array */
		nj = ni*pParam->nLibNum+nLibId;
		printf("Loading %s ...\n", pParam->vArrayFile[nj]->m_pString);
		sprintf(strFileName, "%s%s", pParam->strWorkPath, pParam->vArrayFile[nj]->m_pString);
		pDataNew = NULL;
		pDataNew = Affy_LoadBAR_Fast(strFileName);
		if(pDataNew == NULL)
		{
			printf("Error: TileMapv2_ProbeSelection_Main, cannot load raw data!\n");
			exit(EXIT_FAILURE);
		}

		/* if the first array, prepare enough space for nSampleNum arrays + 
			1 positions + 1 test statistics
		   + nGroupNum means + nGroupNum sample variances + nGroupNum sample sizes 
		   + nVargroupNum pooled variances + nVargroupNum pooled df + 1 FC */ 
		if(ni == 0)
		{
			pData = pDataNew;
			pData->nColNum = pParam->nSampleNum+2+4*pParam->nGroupNum+2*pParam->nVargroupNum+1;
			pFieldType = pData->pFieldType;
			pData->pFieldType = NULL;
			pData->pFieldType = CreateIntMatrix(1,pData->nColNum);
			if(pData->pFieldType == NULL)
			{
				printf("Error: TileMapv2_ProbeSelection_Main, cannot create memory for field types!\n");
				exit(EXIT_FAILURE);
			}
			pData->pFieldType->pMatElement[0] = pFieldType->pMatElement[0];
			pData->pFieldType->pMatElement[1] = pFieldType->pMatElement[1];
			DestroyIntMatrix(pFieldType);
			
			for(nj=0; nj<pData->nSeqNum; nj++)
			{
				nTotalProbeNum += pData->vSeqData[nj]->nDataNum;
				pData->vSeqData[nj]->nColNum = pData->nColNum;
				vDataVec = pData->vSeqData[nj]->vData;
				pData->vSeqData[nj]->vData = NULL;
				pData->vSeqData[nj]->vData = (struct DOUBLEMATRIX **)calloc(pData->vSeqData[nj]->nColNum, sizeof(struct DOUBLEMATRIX *));
				if(pData->vSeqData[nj]->vData == NULL)
				{
					printf("Error: TileMapv2_ProbeSelection_Main, cannot create memory for tracking intensity data!\n");
					exit(EXIT_FAILURE);
				}
				pData->vSeqData[nj]->vData[0] = vDataVec[0];
				pData->vSeqData[nj]->vData[1] = vDataVec[1];
				free(vDataVec);
			}
		}
		/* if not the first array, transfer to the first array */
		else
		{
			if(pDataNew->nSeqNum != pData->nSeqNum)
			{
				printf("Error: TileMapv2_ProbeSelection_Main, array types do not match!\n");
				exit(EXIT_FAILURE);
			}

			for(nj=0; nj<pData->nSeqNum; nj++)
			{
				if(pData->vSeqData[nj]->nDataNum != pDataNew->vSeqData[nj]->nDataNum)
				{
					printf("Error: TileMapv2_ProbeSelection_Main, array types do not match!\n");
					exit(EXIT_FAILURE);
				}
				pData->vSeqData[nj]->vData[ni+1] = pDataNew->vSeqData[nj]->vData[1];
				pDataNew->vSeqData[nj]->vData[1] = NULL;
			}

			pData->pFieldType->pMatElement[ni+1] = pDataNew->pFieldType->pMatElement[1];

			Affy_BARData_Destroy(&pDataNew);
		}
	}

	/* -------------------- */
	/* load masks if needed */
	/* -------------------- */
	if(pParam->nNoiseMask == 1)
	{
		for(ni=0; ni<pParam->nSampleNum; ni++)
		{
			/* load individual array */
			nj = ni*pParam->nLibNum+nLibId;
			printf("Loading %s.mask ...\n", pParam->vArrayFile[nj]->m_pString);
			sprintf(strFileName, "%s%s.mask", pParam->strWorkPath, pParam->vArrayFile[nj]->m_pString);
			pMaskNew = NULL;
			pMaskNew = Affy_LoadBAR_Fast(strFileName);
			if(pMaskNew == NULL)
			{
				printf("Error: TileMapv2_ProbeSelection_Main, cannot load mask data!\n");
				exit(EXIT_FAILURE);
			}

			/* if the first array, prepare enough space for nSampleNum arrays + 
				1 positions + 1 final mask */ 
			if(ni == 0)
			{
				pMask = pMaskNew;
				pMask->nColNum = pParam->nSampleNum+2;
				pFieldType = pMask->pFieldType;
				pMask->pFieldType = NULL;
				pMask->pFieldType = CreateIntMatrix(1,pMask->nColNum);
				if(pMask->pFieldType == NULL)
				{
					printf("Error: TileMapv2_ProbeSelection_Main, cannot create memory for field types!\n");
					exit(EXIT_FAILURE);
				}
				pMask->pFieldType->pMatElement[0] = pFieldType->pMatElement[0];
				pMask->pFieldType->pMatElement[1] = pFieldType->pMatElement[1];
				DestroyIntMatrix(pFieldType);
				
				for(nj=0; nj<pMask->nSeqNum; nj++)
				{
					pMask->vSeqData[nj]->nColNum = pMask->nColNum;
					vDataVec = pMask->vSeqData[nj]->vData;
					pMask->vSeqData[nj]->vData = NULL;
					pMask->vSeqData[nj]->vData = (struct DOUBLEMATRIX **)calloc(pMask->vSeqData[nj]->nColNum, sizeof(struct DOUBLEMATRIX *));
					if(pMask->vSeqData[nj]->vData == NULL)
					{
						printf("Error: TileMapv2_ProbeSelection_Main, cannot create memory for tracking intensity data!\n");
						exit(EXIT_FAILURE);
					}
					pMask->vSeqData[nj]->vData[0] = vDataVec[0];
					pMask->vSeqData[nj]->vData[1] = vDataVec[1];
					free(vDataVec);
				}
			}
			/* if not the first array, transfer to the first array */
			else
			{
				if(pMaskNew->nSeqNum != pMask->nSeqNum)
				{
					printf("Error: TileMapv2_ProbeSelection_Main, array types do not match!\n");
					exit(EXIT_FAILURE);
				}

				for(nj=0; nj<pMask->nSeqNum; nj++)
				{
					if(pMask->vSeqData[nj]->nDataNum != pMaskNew->vSeqData[nj]->nDataNum)
					{
						printf("Error: TileMapv2_ProbeSelection_Main, array types do not match!\n");
						exit(EXIT_FAILURE);
					}
					pMask->vSeqData[nj]->vData[ni+1] = pMaskNew->vSeqData[nj]->vData[1];
					pMaskNew->vSeqData[nj]->vData[1] = NULL;
				}

				pMask->pFieldType->pMatElement[ni+1] = pMaskNew->pFieldType->pMatElement[1];

				Affy_BARData_Destroy(&pMaskNew);
			}
		}
	}

	printf("Total Number of Probes = %d\n", nTotalProbeNum);

	/* -------------------- */
	/* transformation       */
	/* -------------------- */
	TileMapv2_ProbeSelection_DataTransformation(pData, pParam);

	/* -------------------- */
	/* probe level summary  */
	/* -------------------- */
	if(pParam->nPatternType == 1)
	{
		/* compute test statistics */
		sprintf(strFileName, "%s%s_%d.pb.bar", pParam->strWorkPath, pParam->strProjectTitle, nLibId+1);
		sprintf(strFCFileName, "%s%s_%d.fc.bar", pParam->strWorkPath, pParam->strProjectTitle, nLibId+1);
		TileMapv2_ProbeSelection_OneSample_WithMask(pData, pMask, pParam, strFileName, strFCFileName);
	}
	else if( (pParam->nPatternType == 2) && (pParam->nNoiseMask == 1) )
	{
		/* compute test statistics */
		sprintf(strFileName, "%s%s_%d.pb.bar", pParam->strWorkPath, pParam->strProjectTitle, nLibId+1);
		sprintf(strFCFileName, "%s%s_%d.fc.bar", pParam->strWorkPath, pParam->strProjectTitle, nLibId+1);
		TileMapv2_ProbeSelection_TwoSample_WithMask(pData, pMask, pParam, strFileName, strFCFileName);
	}
	else if( (pParam->nPatternType == 2) && (pParam->nNoiseMask == 0) )
	{
		/* compute test statistics */
		sprintf(strFileName, "%s%s_%d.pb.bar", pParam->strWorkPath, pParam->strProjectTitle, nLibId+1);
		sprintf(strFCFileName, "%s%s_%d.fc.bar", pParam->strWorkPath, pParam->strProjectTitle, nLibId+1);
		TileMapv2_ProbeSelection_TwoSample(pData, pParam, strFileName, strFCFileName);
	}
	else if(pParam->nNoiseMask == 1)
	{
		/* compute test statistics */
		sprintf(strFileName, "%s%s_%d.pb.bar", pParam->strWorkPath, pParam->strProjectTitle, nLibId+1);
		strcpy(strFCFileName, "");
		TileMapv2_ProbeSelection_MultiSample_WithMask_Fast(pData, pMask, pParam, strFileName, strFCFileName);
	}
	else
	{
		/* compute test statistics */
		sprintf(strFileName, "%s%s_%d.pb.bar", pParam->strWorkPath, pParam->strProjectTitle, nLibId+1);
		strcpy(strFCFileName, "");
		TileMapv2_ProbeSelection_MultiSample_Fast(pData, pParam, strFileName, strFCFileName);
	}

	/* ---------------------------- */
	/* permutations to estimate FDR */
	/* ---------------------------- */
	if( (pParam->nFDRType == 1) && (pParam->nPermutationNum > 0) )
	{
		printf("Permutations...\n");
		
		/* remove mean intensities from each group */
		TileMapv2_ProbeSelection_RemoveMean(pData, pParam, pParam->nNoiseMask, pMask);

		/* cycles */
		for(ni=0; ni<pParam->nPermutationNum; ni++)
		{
			printf("perm %d...\n", ni);

			/* permute the class label */
			pPermGroupLabel = NULL;
			pPermGroupLabel = Expression_GeneSelection_ClassPerm(pParam->nGroupNum, 
				pParam->pGroupLabel, pParam->pGroupSize,
				pParam->nPermgroupNum, pParam->pPermgroupSize, pParam->vPermgroupMap);
			pOriGroupLabel = pParam->pGroupLabel;
			pParam->pGroupLabel = pPermGroupLabel;

			/* compute score */
			if(pParam->nPatternType == 1)
			{
				/* compute test statistics */
				sprintf(strFileName, "%s%s_%d.perm%d.pb.bar", pParam->strWorkPath, pParam->strProjectTitle, nLibId+1, ni);
				sprintf(strFCFileName, "%s%s_%d.perm%d.fc.bar", pParam->strWorkPath, pParam->strProjectTitle, nLibId+1, ni);
				TileMapv2_ProbeSelection_OneSample_WithMask(pData, pMask, pParam, strFileName, strFCFileName);
			}
			else if( (pParam->nPatternType == 2) && (pParam->nNoiseMask == 1) )
			{
				/* compute test statistics */
				sprintf(strFileName, "%s%s_%d.perm%d.pb.bar", pParam->strWorkPath, pParam->strProjectTitle, nLibId+1, ni);
				sprintf(strFCFileName, "%s%s_%d.perm%d.fc.bar", pParam->strWorkPath, pParam->strProjectTitle, nLibId+1, ni);
				TileMapv2_ProbeSelection_TwoSample_WithMask(pData, pMask, pParam, strFileName, strFCFileName);
			}
			else if( (pParam->nPatternType == 2) && (pParam->nNoiseMask == 0) )
			{
				/* compute test statistics */
				sprintf(strFileName, "%s%s_%d.perm%d.pb.bar", pParam->strWorkPath, pParam->strProjectTitle, nLibId+1, ni);
				sprintf(strFCFileName, "%s%s_%d.perm%d.fc.bar", pParam->strWorkPath, pParam->strProjectTitle, nLibId+1, ni);
				TileMapv2_ProbeSelection_TwoSample(pData, pParam, strFileName, strFCFileName);
			}
			else if(pParam->nNoiseMask == 1)
			{
				/* compute test statistics */
				sprintf(strFileName, "%s%s_%d.perm%d.pb.bar", pParam->strWorkPath, pParam->strProjectTitle, nLibId+1, ni);
				strcpy(strFCFileName, "");
				TileMapv2_ProbeSelection_MultiSample_WithMask_Fast(pData, pMask, pParam, strFileName, strFCFileName);
			}
			else
			{
				/* compute test statistics */
				sprintf(strFileName, "%s%s_%d.perm%d.pb.bar", pParam->strWorkPath, pParam->strProjectTitle, nLibId+1, ni);
				strcpy(strFCFileName, "");
				TileMapv2_ProbeSelection_MultiSample_Fast(pData, pParam, strFileName, strFCFileName);
			}
			
			/* release memory */
			pPermGroupLabel = pParam->pGroupLabel;
			pParam->pGroupLabel = pOriGroupLabel;
			DestroyIntMatrix(pPermGroupLabel);
		}
	}

	/* -------------- */
	/* release memory */
	/* -------------- */
	Affy_BARData_Destroy(&pData);
	if(pParam->nNoiseMask == 1)
	{
		Affy_BARData_Destroy(&pMask);
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ProbeSelection_DataTransformation()                          */
/*  TileMapv2 transform raw data.                                          */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_ProbeSelection_DataTransformation(struct tagBARData *pData,
						struct tagTileMapv2Param *pParam)
{
	/* define */
	int ni,nj,nk;
	double dTemp;
	double dLog2 = log(2.0);

	/* init */
	if( (pData == NULL) || (pParam == NULL) )
	{
		return PROC_SUCCESS;
	}

	/* transformation */
	if(pParam->nTransformType == 1)
	{
		/* log2 */
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			for(nj=1; nj<=pParam->nSampleNum; nj++)
			{
				for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
				{
					dTemp = pData->vSeqData[ni]->vData[nj]->pMatElement[nk];
					if(dTemp < pParam->dLower)
						dTemp = pParam->dLower;
					if(dTemp > pParam->dUpper)
						dTemp = pParam->dUpper;
					if(dTemp <= 0.0)
					{
						printf("Error: TileMapv2_ProbeSelection_DataTransformation, log(%f)!\n", dTemp);
						exit(EXIT_FAILURE);
					}
					pData->vSeqData[ni]->vData[nj]->pMatElement[nk] = log(dTemp)/dLog2;
				}
			}
		}
	}
	else if(pParam->nTransformType == 2)
	{
		/* logit */
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			for(nj=1; nj<=pParam->nSampleNum; nj++)
			{
				for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
				{
					dTemp = pData->vSeqData[ni]->vData[nj]->pMatElement[nk];
					if(dTemp < pParam->dLower)
						dTemp = pParam->dLower;
					if(dTemp > pParam->dUpper)
						dTemp = pParam->dUpper;
					if(dTemp <= 0.0)
					{
						printf("Error: TileMapv2_ProbeSelection_DataTransformation, log(%f)!\n", dTemp);
						exit(EXIT_FAILURE);
					}
					if(dTemp >= 1.0)
					{
						printf("Error: TileMapv2_ProbeSelection_DataTransformation, log(%f)!\n", 1.0-dTemp);
						exit(EXIT_FAILURE);
					}
					pData->vSeqData[ni]->vData[nj]->pMatElement[nk] = log(dTemp/(1.0-dTemp));
				}
			}
		}
	}
	else if(pParam->nTransformType == 3)
	{
		/* exp(t)/[1+exp(t)] */
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			for(nj=1; nj<=pParam->nSampleNum; nj++)
			{
				for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
				{
					dTemp = pData->vSeqData[ni]->vData[nj]->pMatElement[nk];
					if(dTemp < pParam->dLower)
						dTemp = pParam->dLower;
					if(dTemp > pParam->dUpper)
						dTemp = pParam->dUpper;
					dTemp = exp(dTemp);
					pData->vSeqData[ni]->vData[nj]->pMatElement[nk] = dTemp/(1.0+dTemp);
				}
			}
		}
	}
	else
	{
		/* identity */
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			for(nj=1; nj<=pParam->nSampleNum; nj++)
			{
				for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
				{
					dTemp = pData->vSeqData[ni]->vData[nj]->pMatElement[nk];
					if(dTemp < pParam->dLower)
						dTemp = pParam->dLower;
					if(dTemp > pParam->dUpper)
						dTemp = pParam->dUpper;
					pData->vSeqData[ni]->vData[nj]->pMatElement[nk] = dTemp;
				}
			}
		}
	}
		
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ProbeSelection_OneSample_WithMask()                          */
/*  TileMapv2 probe level summary: two sample comparisons.                 */
/*  this function can handle masked cells.                                 */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_ProbeSelection_OneSample_WithMask(struct tagBARData *pData, 
						struct tagBARData *pMask, struct tagTileMapv2Param *pParam, 
						char strOutFile[], char strFCFile[])
{
	/* ------ */
	/* define */
	/* ------ */
	char *chp;
	int ni,nj,nk,nl;
	int ng;
	double dRefValue = 0.0;
	int nB = 1;

	int nGroupNum;
	int nSampleNum;
	int nVarGroupNum;
	struct INTMATRIX *pGroupSizeCopy;
	struct INTMATRIX *pDf;

	int nClustId,nDfId;
	double *vExp,*vSum,*vAve,*vSigma,*vMu,*vSize,*vSize2;
	double *vMask;
	double dTemp;
	double dVarTemp;
	struct INTMATRIX *pCol;

	char strLine[MED_LINE_LENGTH];
	
	/* ------------- */
	/* expression    */
	/* ------------- */
	if(pParam->strPattern[0] == '\0')
	{
		printf("Error: TileMapv2_ProbeSelection_OneSample_WithMask, no logic expression formula available!\n");
		exit(EXIT_FAILURE);
	}
	strcpy(strLine, pParam->strPattern);
	chp = strpbrk(strLine, "<>");
	if(chp == NULL)
	{
		printf("Error: TileMapv2_ProbeSelection_OneSample_WithMask, unrecognized pattern!\n");
		exit(EXIT_FAILURE);
	}

	if(*chp == '<')
		nB = 0;
	else
		nB = 1;

	*chp = '\0';
	ng = atoi(strLine)-1;
	chp++;
	dRefValue = atof(chp);

	if( (ng<0) || (ng>=pParam->nGroupNum) )
	{
		printf("Error: TileMapv2_ProbeSelection_OneSample_WithMask, group id out of range!\n");
		exit(EXIT_FAILURE);
	}
	if( pParam->pGroupSize->pMatElement[ng] <= 0 )
	{
		printf("Error: TileMapv2_ProbeSelection_OneSample_WithMask, the group %d is empty!\n", ng+1);
		exit(EXIT_FAILURE);
	}
	
	/* ---- */
	/* init */
	/* ---- */
	nSampleNum = pParam->nSampleNum;
	nGroupNum = pParam->nGroupNum;
	nVarGroupNum = pParam->nVargroupNum;
	pGroupSizeCopy = NULL;
	pGroupSizeCopy = CreateIntMatrix(1, pParam->nGroupNum);
	pDf = NULL;
	pDf = CreateIntMatrix(1, pParam->nVargroupNum);
	
	/* ------------- */
	/* create memory */
	/* ------------- */
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		/* chromosome position, ID = 0 */
		/* raw data, ID = (1:nSampleNum) */
		/* test statistics, ID = nSampleNum+1 */
		nClustId = nSampleNum+1;
		pData->vSeqData[ni]->vData[nClustId] = NULL;
		pData->vSeqData[ni]->vData[nClustId] = CreateDoubleMatrix(1, pData->vSeqData[ni]->nDataNum);
		if( (pData->vSeqData[ni]->nDataNum > 0) && (pData->vSeqData[ni]->vData[nClustId] == NULL))
		{
			printf("Error: TileMapv2_ProbeSelection_OneSample_WithMask, cannot create space for storing test statistics!\n");
			exit(EXIT_FAILURE);
		}

		if(pMask != NULL)
		{
			/* mask, ID = nSampleNum+1 */
			nClustId = nSampleNum+1;
			pMask->vSeqData[ni]->vData[nClustId] = NULL;
			pMask->vSeqData[ni]->vData[nClustId] = CreateDoubleMatrix(1, pMask->vSeqData[ni]->nDataNum);
			if( (pMask->vSeqData[ni]->nDataNum > 0) && (pMask->vSeqData[ni]->vData[nClustId] == NULL))
			{
				printf("Error: TileMapv2_ProbeSelection_OneSample_WithMask, cannot create space for storing masks!\n");
				exit(EXIT_FAILURE);
			}
		}

		/* within sample mean, ID = nSampleNum+2+(0:nGroupNum-1) */
		for(nj=0; nj<nGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nj;
			pData->vSeqData[ni]->vData[nClustId] = NULL;
			pData->vSeqData[ni]->vData[nClustId] = CreateDoubleMatrix(1, pData->vSeqData[ni]->nDataNum);
			if( (pData->vSeqData[ni]->nDataNum > 0) && (pData->vSeqData[ni]->vData[nClustId] == NULL))
			{
				printf("Error: TileMapv2_ProbeSelection_OneSample_WithMask, cannot create space for storing sample mean!\n");
				exit(EXIT_FAILURE);
			}
		}

		/* within sample variance, ID = nSampleNum+2+nGroupNum+(0:nGroupNum-1) */
		for(nj=0; nj<nGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nGroupNum+nj;
			pData->vSeqData[ni]->vData[nClustId] = NULL;
			pData->vSeqData[ni]->vData[nClustId] = CreateDoubleMatrix(1, pData->vSeqData[ni]->nDataNum);
			if( (pData->vSeqData[ni]->nDataNum > 0) && (pData->vSeqData[ni]->vData[nClustId] == NULL))
			{
				printf("Error: TileMapv2_ProbeSelection_OneSample_WithMask, cannot create space for storing within sample variance!\n");
				exit(EXIT_FAILURE);
			}
		}

		/* sample size, ID = nSampleNum+2+2*nGroupNum+(0:nGroupNum-1) */
		for(nj=0; nj<nGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nGroupNum+nGroupNum+nj;
			pData->vSeqData[ni]->vData[nClustId] = NULL;
			pData->vSeqData[ni]->vData[nClustId] = CreateDoubleMatrix(1, pData->vSeqData[ni]->nDataNum);
			if( (pData->vSeqData[ni]->nDataNum > 0) && (pData->vSeqData[ni]->vData[nClustId] == NULL))
			{
				printf("Error: TileMapv2_ProbeSelection_OneSample_WithMask, cannot create space for storing sample size!\n");
				exit(EXIT_FAILURE);
			}
		}

		/* pooled variance, ID = nSampleNum+2+3*nGroupNum+(0:nVarGroupNum-1) */
		for(nj=0; nj<nVarGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nGroupNum+nGroupNum+nGroupNum+nj;
			pData->vSeqData[ni]->vData[nClustId] = NULL;
			pData->vSeqData[ni]->vData[nClustId] = CreateDoubleMatrix(1, pData->vSeqData[ni]->nDataNum);
			if( (pData->vSeqData[ni]->nDataNum > 0) && (pData->vSeqData[ni]->vData[nClustId] == NULL))
			{
				printf("Error: TileMapv2_ProbeSelection_OneSample_WithMask, cannot create space for storing pooled variance!\n");
				exit(EXIT_FAILURE);
			}
		}

		/* pooled degree of freedom, ID = nSampleNum+2+3*nGroupNum+nVarGroupNum+(0:nVarGroupNum-1) */
		for(nj=0; nj<nVarGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nGroupNum+nGroupNum+nGroupNum+nVarGroupNum+nj;
			pData->vSeqData[ni]->vData[nClustId] = NULL;
			pData->vSeqData[ni]->vData[nClustId] = CreateDoubleMatrix(1, pData->vSeqData[ni]->nDataNum);
			if( (pData->vSeqData[ni]->nDataNum > 0) && (pData->vSeqData[ni]->vData[nClustId] == NULL))
			{
				printf("Error: TileMapv2_ProbeSelection_OneSample_WithMask, cannot create space for storing pooled d.f.!\n");
				exit(EXIT_FAILURE);
			}
		}

		/* fold change, ID = nSampleNum+2+3*nGroupNum+2*nVarGroupNum */
		nClustId = nSampleNum+2+nGroupNum+nGroupNum+nGroupNum+nVarGroupNum+nVarGroupNum;
		pData->vSeqData[ni]->vData[nClustId] = NULL;
		pData->vSeqData[ni]->vData[nClustId] = CreateDoubleMatrix(1, pData->vSeqData[ni]->nDataNum);
		if( (pData->vSeqData[ni]->nDataNum > 0) && (pData->vSeqData[ni]->vData[nClustId] == NULL))
		{
			printf("Error: TileMapv2_ProbeSelection_OneSample_WithMask, cannot create space for storing fold change!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	/* ---------------- */
	/* get mean         */
	/* ---------------- */

	/* mean */
	for(nj=0; nj<nSampleNum; nj++)
	{
		nClustId = pParam->pGroupLabel->pMatElement[nj]-1;
		if( (nClustId<0) || (nClustId>=nGroupNum) )
			continue;

		pGroupSizeCopy->pMatElement[nClustId] += 1;
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			if(pData->vSeqData[ni]->nDataNum <= 0)
				continue;

			vSum = pData->vSeqData[ni]->vData[nSampleNum+2+nClustId]->pMatElement;
			vSize = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+nClustId]->pMatElement;
			vExp = pData->vSeqData[ni]->vData[1+nj]->pMatElement;
			if(pMask != NULL)
			{
				vMask = pMask->vSeqData[ni]->vData[1+nj]->pMatElement;
				for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
				{
					if(vMask[nk] < 0.5)
					{
						vSum[nk] += vExp[nk];
						vSize[nk] += 1.0;
					}
				}
			}
			else
			{
				for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
				{
					vSum[nk] += vExp[nk];
					vSize[nk] += 1.0;
				}
			}
		}
	}

	for(nj=0; nj<nGroupNum; nj++)
	{
		if(pGroupSizeCopy->pMatElement[nj] != pParam->pGroupSize->pMatElement[nj])
		{
			printf("Error: TileMapv2_ProbeSelection_OneSample_WithMask, class size not match!\n");
			exit(EXIT_FAILURE);
		}

		if(pGroupSizeCopy->pMatElement[nj] > 0)
		{
            for(ni=0; ni<pData->nSeqNum; ni++)
			{
				if(pData->vSeqData[ni]->nDataNum <= 0)
					continue;
				
				vSum = pData->vSeqData[ni]->vData[nSampleNum+2+nj]->pMatElement;
				vSize = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+nj]->pMatElement;
				for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
				{
					if( vSize[nk] > 0.5)
						vSum[nk] /= vSize[nk];
				}
			}
		}
	}

	/* ------------------------------------------------------------------------ */
	/* if each group has only one replicate, compute the difference of the mean */
	/* as the test statistics.                                                  */
	/* ------------------------------------------------------------------------ */
	if( pParam->pGroupSize->pMatElement[ng] == 1 )
	{
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			if(pData->vSeqData[ni]->nDataNum <= 0)
				continue;
			
			vSum = pData->vSeqData[ni]->vData[nSampleNum+1]->pMatElement;
			vExp = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+nGroupNum+nVarGroupNum+nVarGroupNum]->pMatElement;
			vMu = pData->vSeqData[ni]->vData[nSampleNum+2+ng]->pMatElement;
			vSize = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+ng]->pMatElement;
			
			/* with mask */
			if(pMask != NULL)
			{
				vMask = pMask->vSeqData[ni]->vData[nSampleNum+1]->pMatElement;
				
				for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
				{
					if(vSize[nk] > 0.5)
					{
						if(nB == 0)
							vSum[nk] = dRefValue-vMu[nk];
						else
							vSum[nk] = vMu[nk]-dRefValue;
							
						vExp[nk] = vSum[nk];
						vMask[nk] = 0.0;
					}
					else
					{
						vMask[nk] = 1.0;
					}
				}				
			}
			/* without mask */
			else
			{
				for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
				{
					if(nB == 0)
						vSum[nk] = dRefValue-vMu[nk];
					else
						vSum[nk] = vMu[nk]-dRefValue;
							
					vExp[nk] = vSum[nk];
				}				
			}
		}

		/* destroy unnecessary memory and return */
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			/* chromosome position, ID = 0 */
			/* raw data, ID = (1:nSampleNum) */
			/* test statistics, ID = nSampleNum+1 */
		
			/* within sample mean, ID = nSampleNum+2+(0:nGroupNum-1) */
			for(nj=0; nj<nGroupNum; nj++)
			{
				nClustId = nSampleNum+2+nj;
				DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
				pData->vSeqData[ni]->vData[nClustId] = NULL;
			}
			/* within sample variance, ID = nSampleNum+2+nGroupNum+(0:nGroupNum-1) */
			for(nj=0; nj<nGroupNum; nj++)
			{
				nClustId = nSampleNum+2+nGroupNum+nj;
				DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
				pData->vSeqData[ni]->vData[nClustId] = NULL;
			}
			/* sample size, ID = nSampleNum+2+nGroupNum+nGroupNum+(0:nGroupNum-1) */
			for(nj=0; nj<nGroupNum; nj++)
			{
				nClustId = nSampleNum+2+nGroupNum+nGroupNum+nj;
				DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
				pData->vSeqData[ni]->vData[nClustId] = NULL;
			}
			/* pooled variance, ID = nSampleNum+2+3*nGroupNum+(0:nVarGroupNum-1) */
			for(nj=0; nj<nVarGroupNum; nj++)
			{
				nClustId = nSampleNum+2+nGroupNum+nGroupNum+nGroupNum+nj;
				DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
				pData->vSeqData[ni]->vData[nClustId] = NULL;
			}
			/* pooled d.f., ID = nSampleNum+2+3*nGroupNum+nVarGroupNum+(0:nVarGroupNum-1) */
			for(nj=0; nj<nVarGroupNum; nj++)
			{
				nClustId = nSampleNum+2+nGroupNum+nGroupNum+nGroupNum+nVarGroupNum+nj;
				DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
				pData->vSeqData[ni]->vData[nClustId] = NULL;
			}

			nClustId = nSampleNum+2+nGroupNum+nGroupNum+nGroupNum+nVarGroupNum+nVarGroupNum;
			pData->vSeqData[ni]->vData[nSampleNum+2] = pData->vSeqData[ni]->vData[nClustId];
			pData->vSeqData[ni]->vData[nClustId] = NULL;
		}

		pData->pFieldType->pMatElement[nSampleNum+1] = 1;
        pData->pFieldType->pMatElement[nSampleNum+2] = 1;
		if(pMask != NULL)
			pMask->pFieldType->pMatElement[nSampleNum+1] = 7;
			
		DestroyIntMatrix(pGroupSizeCopy);
		DestroyIntMatrix(pDf);

		/* -------------- */
		/* save data      */
		/* -------------- */
		pCol = NULL;
		pCol = CreateIntMatrix(1,pData->nColNum);
		if(pCol == NULL)
		{
			printf("Error: TileMapv2_ProbeSelection_OneSample_WithMask, cannot create output column information!\n");
			exit(EXIT_FAILURE);
		}
		
		printf("Exporting probe level test statistics...\n");
		pCol->pMatElement[0] = 1;
		pCol->pMatElement[nSampleNum+1] = 1;
		Affy_SaveBAR_Columns_Fast(strOutFile, pData, pCol);

		printf("Exporting fold changes...\n");
		pCol->pMatElement[nSampleNum+1] = 0;
		pCol->pMatElement[nSampleNum+2] = 1;
		Affy_SaveBAR_Columns_Fast(strFCFile, pData, pCol);

		DestroyIntMatrix(pCol);

		if(pMask != NULL)
		{
			pCol = NULL;
			pCol = CreateIntMatrix(1,pMask->nColNum);
			if(pCol == NULL)
			{
				printf("Error: TileMapv2_ProbeSelection_OneSample_WithMask, cannot create output column information!\n");
				exit(EXIT_FAILURE);
			}
			
			printf("Exporting masks...\n");
			pCol->pMatElement[0] = 1;
			pCol->pMatElement[nSampleNum+1] = 1;
			sprintf(strLine, "%s.mask", strOutFile);
			Affy_SaveBAR_Columns_Fast(strLine, pMask, pCol);

			DestroyIntMatrix(pCol);
		}

		/* -------------- */
		/* release memory */
		/* -------------- */
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			/* chromosome position, ID = 0 */
			/* raw data, ID = (1:nSampleNum) */
			/* test statistics, ID = nSampleNum+1 */
			nClustId = nSampleNum+1;
			DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
			pData->vSeqData[ni]->vData[nClustId] = NULL;

			nClustId = nSampleNum+2;
			DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
			pData->vSeqData[ni]->vData[nClustId] = NULL;
		}
		pData->pFieldType->pMatElement[nSampleNum+1] = 0;
		pData->pFieldType->pMatElement[nSampleNum+2] = 0;

		if(pMask != NULL)
		{
			for(ni=0; ni<pMask->nSeqNum; ni++)
			{
				nClustId = nSampleNum+1;
				DestroyDoubleMatrix(pMask->vSeqData[ni]->vData[nClustId]);
				pMask->vSeqData[ni]->vData[nClustId] = NULL;
			}
			pMask->pFieldType->pMatElement[nSampleNum+1] = 0;
		}

		/* return */
		return PROC_SUCCESS;
	}
		
	/* --------------------------------------------- */
	/* if one has extra d.f., compute the variance.  */
	/* --------------------------------------------- */
	for(nj=0; nj<nSampleNum; nj++)
	{
		nClustId = pParam->pGroupLabel->pMatElement[nj]-1;
		if( (nClustId<0) || (nClustId>=nGroupNum) )
			continue;
		
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			if(pData->vSeqData[ni]->nDataNum <= 0)
				continue;
			
			vSize = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+nClustId]->pMatElement;
			vSum = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nClustId]->pMatElement;
			vAve = pData->vSeqData[ni]->vData[nSampleNum+2+nClustId]->pMatElement;
			vExp = pData->vSeqData[ni]->vData[1+nj]->pMatElement;

			if(pMask != NULL)
			{
				vMask = pMask->vSeqData[ni]->vData[1+nj]->pMatElement;
				for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
				{
					if(vMask[nk] < 0.5)
					{
						dTemp = (vExp[nk]-vAve[nk]);		
						vSum[nk] += dTemp*dTemp;
					}
				}
			}
			else
			{
				for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
				{
					dTemp = (vExp[nk]-vAve[nk]);		
					vSum[nk] += dTemp*dTemp;
				}
			}
		}
	}

	/* variance: estimates */
	for(nj=0; nj<nVarGroupNum; nj++)
	{
		if(pParam->vVargroupMap[nj]->nWidth != pParam->pVargroupSize->pMatElement[nj])
		{
			printf("Error: TileMapv2_ProbeSelection_OneSample_WithMask, variance group number not match!\n");
			exit(EXIT_FAILURE);
		}
	}

	for(nj=0; nj<nVarGroupNum; nj++)
	{
		for(nl=0; nl<pParam->vVargroupMap[nj]->nWidth; nl++)
		{
			nClustId = pParam->vVargroupMap[nj]->pMatElement[nl]-1;
			if(pParam->pGroupSize->pMatElement[nClustId] > 0)
				pDf->pMatElement[nj] += pParam->pGroupSize->pMatElement[nClustId]-1;
			
			for(ni=0; ni<pData->nSeqNum; ni++)
			{
				if(pData->vSeqData[ni]->nDataNum <= 0)
					continue;

				vSum = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+nGroupNum+nj]->pMatElement;
				vSize2 = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+nGroupNum+nVarGroupNum+nj]->pMatElement;
				vExp = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nClustId]->pMatElement;
				vSize = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+nClustId]->pMatElement;
				for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
				{
					if(vSize[nk] > 0.5)
					{
						vSum[nk] += vExp[nk];
						vSize2[nk] += vSize[nk]-1.0;
					}
				}
			}
		}
	}

	for(nj=0; nj<nVarGroupNum; nj++)
	{
		for(nl=0; nl<pParam->vVargroupMap[nj]->nWidth; nl++)
		{
			nClustId = pParam->vVargroupMap[nj]->pMatElement[nl]-1;
			if( nClustId == ng )
			{
				if(pDf->pMatElement[nj] <= 0)
				{
					printf("Error: TileMapv2_ProbeSelection_OneSample_WithMask, group %d does not have enough d.f. to estimate variance!\n", (nClustId+1));
					exit(EXIT_FAILURE);
				}
			}
		}
	}

	/* ------------- */
	/* shrinking var */
	/* ------------- */
	for(nj=0; nj<nVarGroupNum; nj++)
	{
		nClustId = nSampleNum+2+nGroupNum+nGroupNum+nGroupNum+nj;
		nDfId = nSampleNum+2+nGroupNum+nGroupNum+nGroupNum+nVarGroupNum+nj;
		TileMapv2_GetShrinkingVar_FromSSR(pData, pDf->pMatElement[nj], nClustId, nDfId);
	}

	/* destroy old variances */
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		for(nj=0; nj<pParam->nGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nGroupNum+nj;
			DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
			pData->vSeqData[ni]->vData[nClustId] = NULL;
		}
	}


	/* assign new variances */
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		for(nj=0; nj<nVarGroupNum; nj++)
		{
			for(nl=0; nl<pParam->vVargroupMap[nj]->nWidth; nl++)
			{
				nClustId = pParam->vVargroupMap[nj]->pMatElement[nl]-1;
				pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nClustId] = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+nGroupNum+nj];
			}
		}
	}

	/* ----------- */
	/* simulation  */
	/* ----------- */
	/* compute probeset level summary */
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		if(pData->vSeqData[ni]->nDataNum <= 0)
			continue;

		vMu = pData->vSeqData[ni]->vData[nSampleNum+2+ng]->pMatElement;
		vSigma = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+ng]->pMatElement;
		vSize = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+ng]->pMatElement;
		
		vSum = pData->vSeqData[ni]->vData[nSampleNum+1]->pMatElement;
		vExp = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+nGroupNum+nVarGroupNum+nVarGroupNum]->pMatElement;
		
		if(pMask != NULL)
		{
			vMask = pMask->vSeqData[ni]->vData[nSampleNum+1]->pMatElement;
			for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
			{
				if( (vSize[nk]<0.5) || (vSigma[nk]<=0.0) )
				{
					vMask[nk] = 1.0;
				}
				else
				{
					vMask[nk] = 0.0;
					if(nB == 0)
                        vExp[nk] = dRefValue-vMu[nk];
					else
						vExp[nk] = vMu[nk]-dRefValue;

					dVarTemp = sqrt(vSigma[nk]/vSize[nk]);
					vSum[nk] = vExp[nk]/dVarTemp;
				}
			}
		}
		else
		{
			for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
			{
				if(nB == 0)
					vExp[nk] = dRefValue-vMu[nk];
				else
					vExp[nk] = vMu[nk]-dRefValue;

				dVarTemp = sqrt(vSigma[nk]/vSize[nk]);
				vSum[nk] = vExp[nk]/dVarTemp;
			}
		}
	}
	
	/* -------------- */
	/* release memory */
	/* -------------- */
	/* destroy unnecessary memory and return */
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		/* chromosome position, ID = 0 */
		/* raw data, ID = (1:nSampleNum) */
		/* test statistics, ID = nSampleNum+1 */
	
		/* within sample mean, ID = nSampleNum+2+(0:nGroupNum-1) */
		for(nj=0; nj<nGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nj;
			DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
			pData->vSeqData[ni]->vData[nClustId] = NULL;
		}
		/* within sample variance, ID = nSampleNum+2+nGroupNum+(0:nGroupNum-1) */
		for(nj=0; nj<nGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nGroupNum+nj;
			pData->vSeqData[ni]->vData[nClustId] = NULL;
		}
		/* sample size, ID = nSampleNum+2+nGroupNum+nGroupNum+(0:nGroupNum-1) */
		for(nj=0; nj<nGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nGroupNum+nGroupNum+nj;
			DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
			pData->vSeqData[ni]->vData[nClustId] = NULL;
		}
		/* pooled variance, ID = nSampleNum+2+3*nGroupNum+(0:nVarGroupNum-1) */
		for(nj=0; nj<nVarGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nGroupNum+nGroupNum+nGroupNum+nj;
			DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
			pData->vSeqData[ni]->vData[nClustId] = NULL;
		}
		/* pooled d.f., ID = nSampleNum+2+3*nGroupNum+nVarGroupNum+(0:nVarGroupNum-1) */
		for(nj=0; nj<nVarGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nGroupNum+nGroupNum+nGroupNum+nVarGroupNum+nj;
			DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
			pData->vSeqData[ni]->vData[nClustId] = NULL;
		}

		nClustId = nSampleNum+2+nGroupNum+nGroupNum+nGroupNum+nVarGroupNum+nVarGroupNum;
		pData->vSeqData[ni]->vData[nSampleNum+2] = pData->vSeqData[ni]->vData[nClustId];
		pData->vSeqData[ni]->vData[nClustId] = NULL;
	}

	pData->pFieldType->pMatElement[nSampleNum+1] = 1;
	pData->pFieldType->pMatElement[nSampleNum+2] = 1;
	if(pMask != NULL)
		pMask->pFieldType->pMatElement[nSampleNum+1] = 7;

	DestroyIntMatrix(pGroupSizeCopy);
	DestroyIntMatrix(pDf);

	/* -------------- */
	/* save data      */
	/* -------------- */
	pCol = NULL;
	pCol = CreateIntMatrix(1,pData->nColNum);
	if(pCol == NULL)
	{
		printf("Error: TileMapv2_ProbeSelection_OneSample_WithMask, cannot create output column information!\n");
		exit(EXIT_FAILURE);
	}
	
	printf("Exporting probe level test statistics...\n");
	pCol->pMatElement[0] = 1;
	pCol->pMatElement[nSampleNum+1] = 1;
	Affy_SaveBAR_Columns_Fast(strOutFile, pData, pCol);

    printf("Exporting fold changes...\n");
	pCol->pMatElement[nSampleNum+1] = 0;
	pCol->pMatElement[nSampleNum+2] = 1;
	Affy_SaveBAR_Columns_Fast(strFCFile, pData, pCol);

	DestroyIntMatrix(pCol);

	if(pMask != NULL)
	{
		pCol = NULL;
		pCol = CreateIntMatrix(1,pMask->nColNum);
		if(pCol == NULL)
		{
			printf("Error: TileMapv2_ProbeSelection_OneSample_WithMask, cannot create output column information!\n");
			exit(EXIT_FAILURE);
		}
		
		printf("Exporting masks...\n");
		pCol->pMatElement[0] = 1;
		pCol->pMatElement[nSampleNum+1] = 1;
		sprintf(strLine, "%s.mask", strOutFile);
		Affy_SaveBAR_Columns_Fast(strLine, pMask, pCol);

		DestroyIntMatrix(pCol);
	}

	/* -------------- */
	/* release memory */
	/* -------------- */
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		/* chromosome position, ID = 0 */
		/* raw data, ID = (1:nSampleNum) */
		/* test statistics, ID = nSampleNum+1 */
		nClustId = nSampleNum+1;
		DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
		pData->vSeqData[ni]->vData[nClustId] = NULL;

		nClustId = nSampleNum+2;
		DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
		pData->vSeqData[ni]->vData[nClustId] = NULL;
	}
	pData->pFieldType->pMatElement[nSampleNum+1] = 0;
	pData->pFieldType->pMatElement[nSampleNum+2] = 0;

	if(pMask != NULL)
	{
		for(ni=0; ni<pMask->nSeqNum; ni++)
		{
			nClustId = nSampleNum+1;
			DestroyDoubleMatrix(pMask->vSeqData[ni]->vData[nClustId]);
			pMask->vSeqData[ni]->vData[nClustId] = NULL;
		}
		pMask->pFieldType->pMatElement[nSampleNum+1] = 0;
	}

	/* ------ */
	/* return */
	/* ------ */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ProbeSelection_TwoSample()                                   */
/*  TileMapv2 probe level summary: two sample comparisons.                 */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_ProbeSelection_TwoSample(struct tagBARData *pData,
						struct tagTileMapv2Param *pParam, char strOutFile[],
						char strFCFile[])
{
	/* ------ */
	/* define */
	/* ------ */
	int nLogicLen;
	int ni,nj,nk,nl,nLen;
	int nLeftNum,nRightNum,nTNumLen;
	char vLogic[MED_LINE_LENGTH];
	char vTNumber[LINE_LENGTH];
	double vGid[MED_LINE_LENGTH];
	int ng1,ng2;

	int nGroupNum;
	int nSampleNum;
	int nVarGroupNum;
	struct INTMATRIX *pGroupSizeCopy;
	struct INTMATRIX *pDf;

	int nClustId;
	double *vExp,*vSum,*vAve,*vSigma,*vMu,*vSigma2,*vMu2;
	double dDenom;
	double dTemp;
	double dVarmean,dVarsst;
	double dB,dV2,dK,dN;
	int nTotalProbeNum;
	struct DOUBLEMATRIX *pSigCoef;
	double dExp,dExp2,dVarTemp;
	struct INTMATRIX *pCol;
	
	/* ------------- */
	/* expression    */
	/* ------------- */
	nLogicLen = 0;
	if(pParam->strPattern[0] == '\0')
	{
		printf("Error: TileMapv2_ProbeSelection_TwoSample, no logic expression formula available!\n");
		exit(EXIT_FAILURE);
	}
	nLen = (int)strlen(pParam->strPattern);
	nj = 0;
	nLeftNum = 0;
	nRightNum = 0;
	nTNumLen = 0;

	for(ni=0; ni<nLen; ni++)
	{
		if( pParam->strPattern[ni] == '(' )
		{
			if(nTNumLen != 0)
			{
				printf("Error: TileMapv2_ProbeSelection_TwoSample, logic expression wrong!\n");
				exit(EXIT_FAILURE);
			}

			nLeftNum++;
		}
		else if(pParam->strPattern[ni] == ')')
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
		else if( (pParam->strPattern[ni] == '<') || (pParam->strPattern[ni] == '>') 
			|| (pParam->strPattern[ni] == '&') || (pParam->strPattern[ni] == '|') )
		{
			if(nTNumLen > 0)
			{
				vTNumber[nTNumLen] = '\0';
				vLogic[nj] = 'G';
				vGid[nj] = atof(vTNumber);
				nTNumLen = 0;
				nj++;
			}
			vLogic[nj] = pParam->strPattern[ni];
			vGid[nj] = -1.0;
			nj++;	
		}
		else if( (pParam->strPattern[ni] >= '0') && (pParam->strPattern[ni] <= '9') )
		{
			vTNumber[nTNumLen] = pParam->strPattern[ni];
			nTNumLen++;
		}
		else if( (pParam->strPattern[ni] == ' ') || (pParam->strPattern[ni] == '\t') )
		{
			/* do nothing */
		}
		else
		{
			printf("Error: TileMapv2_ProbeSelection_TwoSample, %c not supported in logic expressions!\n", pParam->strPattern[ni]);
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
		printf("Error: TileMapv2_ProbeSelection_TwoSample, lack '(' or ')'!\n");
		exit(EXIT_FAILURE);
	}

	if(nLogicLen != 3)
	{
		printf("Error: TileMapv2_ProbeSelection_TwoSample, logic expression wrong!\n");
		printf("       Not a two sample comparison!\n");
		exit(EXIT_FAILURE);
	}
	if((vLogic[0] != 'G') || (vLogic[2] != 'G'))
	{
		printf("Error: TileMapv2_ProbeSelection_TwoSample, logic expression wrong!\n");
		exit(EXIT_FAILURE);
	}
	if((vLogic[1] != '<') && (vLogic[1] != '>'))
	{
		printf("Error: TileMapv2_ProbeSelection_TwoSample, logic expression wrong!\n");
		exit(EXIT_FAILURE);
	}
	if(vLogic[1] == '<')
	{
		ng1 = (int)vGid[2]-1;
		ng2 = (int)vGid[0]-1;
	}
	else if(vLogic[1] == '>')
	{
		ng1 = (int)vGid[0]-1;
		ng2 = (int)vGid[2]-1;
	}
	else
	{
		printf("Error: TileMapv2_ProbeSelection_TwoSample, logic expression wrong!\n");
		exit(EXIT_FAILURE);
	}

	if((ng1<0) || (ng1>=pParam->nGroupNum) || (ng2<0) || (ng2>=pParam->nGroupNum))
	{
		printf("Error: TileMapv2_ProbeSelection_TwoSample, group id out of range!\n");
		exit(EXIT_FAILURE);
	}
	if( (pParam->pGroupSize->pMatElement[ng1] <= 0) || (pParam->pGroupSize->pMatElement[ng2] <= 0) )
	{
		printf("Error: TileMapv2_ProbeSelection_TwoSample, at least one group is empty!\n");
		exit(EXIT_FAILURE);
	}
	
	/* ---- */
	/* init */
	/* ---- */
	nSampleNum = pParam->nSampleNum;
	nGroupNum = pParam->nGroupNum;
	nVarGroupNum = pParam->nVargroupNum;
	pGroupSizeCopy = NULL;
	pGroupSizeCopy = CreateIntMatrix(1, pParam->nGroupNum);
	pDf = NULL;
	pDf = CreateIntMatrix(1, pParam->nVargroupNum);

	/* ------------- */
	/* create memory */
	/* ------------- */
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		/* chromosome position, ID = 0 */
		/* raw data, ID = (1:nSampleNum) */
		/* test statistics, ID = nSampleNum+1 */
		nClustId = nSampleNum+1;
		pData->vSeqData[ni]->vData[nClustId] = NULL;
		pData->vSeqData[ni]->vData[nClustId] = CreateDoubleMatrix(1, pData->vSeqData[ni]->nDataNum);
		if( (pData->vSeqData[ni]->nDataNum > 0) && (pData->vSeqData[ni]->vData[nClustId] == NULL))
		{
			printf("Error: TileMapv2_ProbeSelection_TwoSample, cannot create space for storing test statistics!\n");
			exit(EXIT_FAILURE);
		}

		/* within sample mean, ID = nSampleNum+2+(0:nGroupNum-1) */
		for(nj=0; nj<nGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nj;
			pData->vSeqData[ni]->vData[nClustId] = NULL;
			pData->vSeqData[ni]->vData[nClustId] = CreateDoubleMatrix(1, pData->vSeqData[ni]->nDataNum);
			if( (pData->vSeqData[ni]->nDataNum > 0) && (pData->vSeqData[ni]->vData[nClustId] == NULL))
			{
				printf("Error: TileMapv2_ProbeSelection_TwoSample, cannot create space for storing sample mean!\n");
				exit(EXIT_FAILURE);
			}
		}

		/* within sample variance, ID = nSampleNum+2+nGroupNum+(0:nGroupNum-1) */
		for(nj=0; nj<nGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nGroupNum+nj;
			pData->vSeqData[ni]->vData[nClustId] = NULL;
			pData->vSeqData[ni]->vData[nClustId] = CreateDoubleMatrix(1, pData->vSeqData[ni]->nDataNum);
			if( (pData->vSeqData[ni]->nDataNum > 0) && (pData->vSeqData[ni]->vData[nClustId] == NULL))
			{
				printf("Error: TileMapv2_ProbeSelection_TwoSample, cannot create space for storing within sample variance!\n");
				exit(EXIT_FAILURE);
			}
		}

		/* pooled variance, ID = nSampleNum+2+2*nGroupNum+(0:nVarGroupNum-1) */
		for(nj=0; nj<nVarGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nGroupNum+nGroupNum+nj;
			pData->vSeqData[ni]->vData[nClustId] = NULL;
			pData->vSeqData[ni]->vData[nClustId] = CreateDoubleMatrix(1, pData->vSeqData[ni]->nDataNum);
			if( (pData->vSeqData[ni]->nDataNum > 0) && (pData->vSeqData[ni]->vData[nClustId] == NULL))
			{
				printf("Error: TileMapv2_ProbeSelection_TwoSample, cannot create space for storing pooled variance!\n");
				exit(EXIT_FAILURE);
			}
		}

		/* fold change, ID = nSampleNum+2+2*nGroupNum+nVarGroupNum */
		nClustId = nSampleNum+2+nGroupNum+nGroupNum+nVarGroupNum;
		pData->vSeqData[ni]->vData[nClustId] = NULL;
		pData->vSeqData[ni]->vData[nClustId] = CreateDoubleMatrix(1, pData->vSeqData[ni]->nDataNum);
		if( (pData->vSeqData[ni]->nDataNum > 0) && (pData->vSeqData[ni]->vData[nClustId] == NULL))
		{
			printf("Error: TileMapv2_ProbeSelection_TwoSample, cannot create space for storing fold change!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	/* ---------------- */
	/* get mean         */
	/* ---------------- */

	/* mean */
	for(nj=0; nj<nSampleNum; nj++)
	{
		nClustId = pParam->pGroupLabel->pMatElement[nj]-1;
		if( (nClustId<0) || (nClustId>=nGroupNum) )
			continue;

		pGroupSizeCopy->pMatElement[nClustId] += 1;
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			if(pData->vSeqData[ni]->nDataNum <= 0)
				continue;

			vSum = pData->vSeqData[ni]->vData[nSampleNum+2+nClustId]->pMatElement;
			vExp = pData->vSeqData[ni]->vData[1+nj]->pMatElement;
			for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
			{
				vSum[nk] += vExp[nk];
			}
		}
	}

	for(nj=0; nj<nGroupNum; nj++)
	{
		if(pGroupSizeCopy->pMatElement[nj] != pParam->pGroupSize->pMatElement[nj])
		{
			printf("Error: TileMapv2_ProbeSelection_TwoSample, class size not match!\n");
			exit(EXIT_FAILURE);
		}

		if(pGroupSizeCopy->pMatElement[nj] > 0)
		{
            dDenom = (double)(pGroupSizeCopy->pMatElement[nj]);
			for(ni=0; ni<pData->nSeqNum; ni++)
			{
				if(pData->vSeqData[ni]->nDataNum <= 0)
					continue;
				
				vSum = pData->vSeqData[ni]->vData[nSampleNum+2+nj]->pMatElement;
				for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
				{
					vSum[nk] /= dDenom;
				}
			}
		}
	}

	/* ------------------------------------------------------------------------ */
	/* if each group has only one replicate, compute the difference of the mean */
	/* as the test statistics.                                                  */
	/* ------------------------------------------------------------------------ */
	if( (pParam->pGroupSize->pMatElement[ng1] == 1) && (pParam->pGroupSize->pMatElement[ng2] == 1) )
	{
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			if(pData->vSeqData[ni]->nDataNum <= 0)
				continue;

			vSum = pData->vSeqData[ni]->vData[nSampleNum+1]->pMatElement;
			vExp = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+nVarGroupNum]->pMatElement;
			vMu = pData->vSeqData[ni]->vData[nSampleNum+2+ng1]->pMatElement;
			vMu2 = pData->vSeqData[ni]->vData[nSampleNum+2+ng2]->pMatElement;
			for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
			{
				vSum[nk] = vMu[nk]-vMu2[nk];
				vExp[nk] = vSum[nk];
			}
		}

		/* destroy unnecessary memory and return */
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			/* chromosome position, ID = 0 */
			/* raw data, ID = (1:nSampleNum) */
			/* test statistics, ID = nSampleNum+1 */
		
			/* within sample mean, ID = nSampleNum+2+(0:nGroupNum-1) */
			for(nj=0; nj<nGroupNum; nj++)
			{
				nClustId = nSampleNum+2+nj;
				DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
				pData->vSeqData[ni]->vData[nClustId] = NULL;
			}
			/* within sample variance, ID = nSampleNum+2+nGroupNum+(0:nGroupNum-1) */
			for(nj=0; nj<nGroupNum; nj++)
			{
				nClustId = nSampleNum+2+nGroupNum+nj;
				DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
				pData->vSeqData[ni]->vData[nClustId] = NULL;
			}
			/* pooled variance, ID = nSampleNum+2+2*nGroupNum+(0:nVarGroupNum-1) */
			for(nj=0; nj<nVarGroupNum; nj++)
			{
				nClustId = nSampleNum+2+nGroupNum+nGroupNum+nj;
				DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
				pData->vSeqData[ni]->vData[nClustId] = NULL;
			}

			nClustId = nSampleNum+2+nGroupNum+nGroupNum+nVarGroupNum;
			pData->vSeqData[ni]->vData[nSampleNum+2] = pData->vSeqData[ni]->vData[nClustId];
			pData->vSeqData[ni]->vData[nClustId] = NULL;
		}

		pData->pFieldType->pMatElement[nSampleNum+1] = 1;
		pData->pFieldType->pMatElement[nSampleNum+2] = 1;
			
		DestroyIntMatrix(pGroupSizeCopy);
		DestroyIntMatrix(pDf);

		/* -------------- */
		/* save data      */
		/* -------------- */
		pCol = NULL;
		pCol = CreateIntMatrix(1,pData->nColNum);
		if(pCol == NULL)
		{
			printf("Error: TileMapv2_ProbeSelection_TwoSample, cannot create output column information!\n");
			exit(EXIT_FAILURE);
		}
		
		printf("Exporting probe level test statistics...\n");
		pCol->pMatElement[0] = 1;
		pCol->pMatElement[nSampleNum+1] = 1;
		Affy_SaveBAR_Columns_Fast(strOutFile, pData, pCol);

		printf("Exporting fold changes...\n");
		pCol->pMatElement[nSampleNum+1] = 0;
		pCol->pMatElement[nSampleNum+2] = 1;
		Affy_SaveBAR_Columns_Fast(strFCFile, pData, pCol);

		DestroyIntMatrix(pCol);

		/* -------------- */
		/* release memory */
		/* -------------- */
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			/* chromosome position, ID = 0 */
			/* raw data, ID = (1:nSampleNum) */
			/* test statistics, ID = nSampleNum+1 */
			nClustId = nSampleNum+1;
			DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
			pData->vSeqData[ni]->vData[nClustId] = NULL;

			nClustId = nSampleNum+2;
			DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
			pData->vSeqData[ni]->vData[nClustId] = NULL;
		}
		pData->pFieldType->pMatElement[nSampleNum+1] = 0;
		pData->pFieldType->pMatElement[nSampleNum+2] = 0;


		/* return */
		return PROC_SUCCESS;
	}
		
	/* --------------------------------------------- */
	/* if one has extra d.f., compute the variance.  */
	/* --------------------------------------------- */
	for(nj=0; nj<nSampleNum; nj++)
	{
		nClustId = pParam->pGroupLabel->pMatElement[nj]-1;
		if( (nClustId<0) || (nClustId>=nGroupNum) )
			continue;
		
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			if(pData->vSeqData[ni]->nDataNum <= 0)
				continue;

			vSum = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nClustId]->pMatElement;
			vAve = pData->vSeqData[ni]->vData[nSampleNum+2+nClustId]->pMatElement;
			vExp = pData->vSeqData[ni]->vData[1+nj]->pMatElement;
						
			for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
			{
				dTemp = (vExp[nk]-vAve[nk]);		
				vSum[nk] += dTemp*dTemp;
			}
		}
	}

	/* variance: estimates */
	for(nj=0; nj<nVarGroupNum; nj++)
	{
		if(pParam->vVargroupMap[nj]->nWidth != pParam->pVargroupSize->pMatElement[nj])
		{
			printf("Error: TileMapv2_ProbeSelection_TwoSample, variance group number not match!\n");
			exit(EXIT_FAILURE);
		}
	}

	for(nj=0; nj<nVarGroupNum; nj++)
	{
		for(nl=0; nl<pParam->vVargroupMap[nj]->nWidth; nl++)
		{
			nClustId = pParam->vVargroupMap[nj]->pMatElement[nl]-1;
			if(pParam->pGroupSize->pMatElement[nClustId] > 0)
				pDf->pMatElement[nj] += pParam->pGroupSize->pMatElement[nClustId]-1;
			
			for(ni=0; ni<pData->nSeqNum; ni++)
			{
				if(pData->vSeqData[ni]->nDataNum <= 0)
					continue;
				
				vSum = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+nj]->pMatElement;
				vExp = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nClustId]->pMatElement;
				for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
				{
					vSum[nk] += vExp[nk];
				}
			}
		}
	}

	for(nj=0; nj<nVarGroupNum; nj++)
	{
		for(nl=0; nl<pParam->vVargroupMap[nj]->nWidth; nl++)
		{
			nClustId = pParam->vVargroupMap[nj]->pMatElement[nl]-1;
			if( (nClustId == ng1) || (nClustId == ng2) )
			{
				if(pDf->pMatElement[nj] <= 0)
				{
					printf("Error: TileMapv2_ProbeSelection_TwoSample, group %d does not have enough d.f. to estimate variance!\n", (nClustId+1));
					exit(EXIT_FAILURE);
				}
			}
		}
	}

	/* ------------- */
	/* shrinking var */
	/* ------------- */
	for(nj=0; nj<nVarGroupNum; nj++)
	{
		dVarmean = 0.0;
		dVarsst = 0.0;
		dB = 0.0;

		if(pDf->pMatElement[nj] > 0)
		{
			dDenom = (double)(pDf->pMatElement[nj]);
			nTotalProbeNum = 0;
			for(ni=0; ni<pData->nSeqNum; ni++)
			{
				if(pData->vSeqData[ni]->nDataNum <= 0)
					continue;
				
				nTotalProbeNum += pData->vSeqData[ni]->nDataNum;
                vSum = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+nj]->pMatElement;
			
				/* xbar */
				for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
				{
					vSum[nk] /= dDenom;
					dVarmean += vSum[nk];
				}
			}
			
			dVarmean /= (double)nTotalProbeNum;

			/* sst */
			for(ni=0; ni<pData->nSeqNum; ni++)
			{
				if(pData->vSeqData[ni]->nDataNum <= 0)
					continue;
				
				vSum = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+nj]->pMatElement;
			
				/* xbar */
				for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
				{
					dTemp = vSum[nk]-dVarmean;
					dVarsst += dTemp*dTemp;
				}
			}

			/* shrinkage factor */
			dV2 = 2.0/dDenom;
			dK = (double)nTotalProbeNum;
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
			for(ni=0; ni<pData->nSeqNum; ni++)
			{
				if(pData->vSeqData[ni]->nDataNum <= 0)
					continue;
				
				vSum = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+nj]->pMatElement;
			
				/* xbar */
				for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
				{
					vSum[nk] = (1-dB)*vSum[nk]+dB*dVarmean+1e-16;
				}
			}
		}
	}

	/* destroy old variances */
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		for(nj=0; nj<pParam->nGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nGroupNum+nj;
			DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
			pData->vSeqData[ni]->vData[nClustId] = NULL;
		}
	}


	/* assign new variances */
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		for(nj=0; nj<nVarGroupNum; nj++)
		{
			for(nl=0; nl<pParam->vVargroupMap[nj]->nWidth; nl++)
			{
				nClustId = pParam->vVargroupMap[nj]->pMatElement[nl]-1;
				pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nClustId] = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+nj];
			}
		}
	}

	/* ----------- */
	/* simulation  */
	/* ----------- */
	pSigCoef = NULL;
	pSigCoef = CreateDoubleMatrix(1, pParam->nGroupNum);
	vExp = pSigCoef->pMatElement;
	for(ni=0; ni<pParam->nGroupNum; ni++)
	{
		if(pParam->pGroupSize->pMatElement[ni] > 0)
		{
			dDenom = (double)(pParam->pGroupSize->pMatElement[ni]);
			vExp[ni] = 1.0/dDenom;
		}
		else
		{
			vExp[ni] = 0.0;
		}
	}

	dExp = pSigCoef->pMatElement[ng1];
	dExp2 = pSigCoef->pMatElement[ng2];
	
	/* compute probeset level summary */
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		if(pData->vSeqData[ni]->nDataNum <= 0)
			continue;
				
		vMu = pData->vSeqData[ni]->vData[nSampleNum+2+ng1]->pMatElement;
		vSigma = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+ng1]->pMatElement;
		vMu2 = pData->vSeqData[ni]->vData[nSampleNum+2+ng2]->pMatElement;
		vSigma2 = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+ng2]->pMatElement;
		vSum = pData->vSeqData[ni]->vData[nSampleNum+1]->pMatElement;
		vExp = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+nVarGroupNum]->pMatElement;
		
		for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
		{
			vExp[nk] = vMu[nk]-vMu2[nk];
			if((vSigma[nk] <= 0.0) || (vSigma2[nk] <= 0.0))
			{
				printf("Warning: TileMapv2_ProbeSelection_TwoSample, variance=0, may not have enough sample to estimate variance!\n");
			}
			dVarTemp = sqrt(vSigma[nk]*dExp+vSigma2[nk]*dExp2);
			if(dVarTemp > 0.0)
				vSum[nk] = (vMu[nk]-vMu2[nk])/dVarTemp;
			else
				vSum[nk] = 0.0;
		}
	}
	
	DestroyDoubleMatrix(pSigCoef);

	/* -------------- */
	/* release memory */
	/* -------------- */
	/* destroy unnecessary memory and return */
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		/* chromosome position, ID = 0 */
		/* raw data, ID = (1:nSampleNum) */
		/* test statistics, ID = nSampleNum+1 */
	
		/* within sample mean, ID = nSampleNum+2+(0:nGroupNum-1) */
		for(nj=0; nj<nGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nj;
			DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
			pData->vSeqData[ni]->vData[nClustId] = NULL;
		}
		/* within sample variance, ID = nSampleNum+2+nGroupNum+(0:nGroupNum-1) */
		for(nj=0; nj<nGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nGroupNum+nj;
			pData->vSeqData[ni]->vData[nClustId] = NULL;
		}
		/* pooled variance, ID = nSampleNum+2+2*nGroupNum+(0:nVarGroupNum-1) */
		for(nj=0; nj<nVarGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nGroupNum+nGroupNum+nj;
			DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
			pData->vSeqData[ni]->vData[nClustId] = NULL;
		}

		nClustId = nSampleNum+2+nGroupNum+nGroupNum+nVarGroupNum;
		pData->vSeqData[ni]->vData[nSampleNum+2] = pData->vSeqData[ni]->vData[nClustId];
		pData->vSeqData[ni]->vData[nClustId] = NULL;
	}

	pData->pFieldType->pMatElement[nSampleNum+1] = 1;
	pData->pFieldType->pMatElement[nSampleNum+2] = 1;

	DestroyIntMatrix(pGroupSizeCopy);
	DestroyIntMatrix(pDf);

	/* -------------- */
	/* save data      */
	/* -------------- */
	pCol = NULL;
	pCol = CreateIntMatrix(1,pData->nColNum);
	if(pCol == NULL)
	{
		printf("Error: TileMapv2_ProbeSelection_TwoSample, cannot create output column information!\n");
		exit(EXIT_FAILURE);
	}
	
	printf("Exporting probe level test statistics...\n");
	pCol->pMatElement[0] = 1;
	pCol->pMatElement[nSampleNum+1] = 1;
	Affy_SaveBAR_Columns_Fast(strOutFile, pData, pCol);

    printf("Exporting fold changes...\n");
	pCol->pMatElement[nSampleNum+1] = 0;
	pCol->pMatElement[nSampleNum+2] = 1;
	Affy_SaveBAR_Columns_Fast(strFCFile, pData, pCol);

	DestroyIntMatrix(pCol);

	/* -------------- */
	/* release memory */
	/* -------------- */
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		/* chromosome position, ID = 0 */
		/* raw data, ID = (1:nSampleNum) */
		/* test statistics, ID = nSampleNum+1 */
		nClustId = nSampleNum+1;
		DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
		pData->vSeqData[ni]->vData[nClustId] = NULL;

		nClustId = nSampleNum+2;
		DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
		pData->vSeqData[ni]->vData[nClustId] = NULL;
	}
	pData->pFieldType->pMatElement[nSampleNum+1] = 0;
	pData->pFieldType->pMatElement[nSampleNum+2] = 0;

	/* ------ */
	/* return */
	/* ------ */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ProbeSelection_TwoSample_WithMask()                          */
/*  TileMapv2 probe level summary: two sample comparisons.                 */
/*  this function can handle masked cells.                                 */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_ProbeSelection_TwoSample_WithMask(struct tagBARData *pData, 
						struct tagBARData *pMask, struct tagTileMapv2Param *pParam, 
						char strOutFile[], char strFCFile[])
{
	/* ------ */
	/* define */
	/* ------ */
	int nLogicLen;
	int ni,nj,nk,nl,nLen;
	int nLeftNum,nRightNum,nTNumLen;
	char vLogic[MED_LINE_LENGTH];
	char vTNumber[LINE_LENGTH];
	double vGid[MED_LINE_LENGTH];
	int ng1,ng2;

	int nGroupNum;
	int nSampleNum;
	int nVarGroupNum;
	struct INTMATRIX *pGroupSizeCopy;
	struct INTMATRIX *pDf;

	int nClustId,nDfId;
	double *vExp,*vSum,*vAve,*vSigma,*vMu,*vSigma2,*vMu2,*vSize,*vSize2;
	double *vMask;
	double dTemp;
	double dVarTemp;
	struct INTMATRIX *pCol;

	char strLine[MED_LINE_LENGTH];
	
	/* ------------- */
	/* expression    */
	/* ------------- */
	nLogicLen = 0;
	if(pParam->strPattern[0] == '\0')
	{
		printf("Error: TileMapv2_ProbeSelection_TwoSample_WithMask, no logic expression formula available!\n");
		exit(EXIT_FAILURE);
	}
	nLen = (int)strlen(pParam->strPattern);
	nj = 0;
	nLeftNum = 0;
	nRightNum = 0;
	nTNumLen = 0;

	for(ni=0; ni<nLen; ni++)
	{
		if( pParam->strPattern[ni] == '(' )
		{
			if(nTNumLen != 0)
			{
				printf("Error: TileMapv2_ProbeSelection_TwoSample_WithMask, logic expression wrong!\n");
				exit(EXIT_FAILURE);
			}

			nLeftNum++;
		}
		else if(pParam->strPattern[ni] == ')')
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
		else if( (pParam->strPattern[ni] == '<') || (pParam->strPattern[ni] == '>') 
			|| (pParam->strPattern[ni] == '&') || (pParam->strPattern[ni] == '|') )
		{
			if(nTNumLen > 0)
			{
				vTNumber[nTNumLen] = '\0';
				vLogic[nj] = 'G';
				vGid[nj] = atof(vTNumber);
				nTNumLen = 0;
				nj++;
			}
			vLogic[nj] = pParam->strPattern[ni];
			vGid[nj] = -1.0;
			nj++;	
		}
		else if( (pParam->strPattern[ni] >= '0') && (pParam->strPattern[ni] <= '9') )
		{
			vTNumber[nTNumLen] = pParam->strPattern[ni];
			nTNumLen++;
		}
		else if( (pParam->strPattern[ni] == ' ') || (pParam->strPattern[ni] == '\t') )
		{
			/* do nothing */
		}
		else
		{
			printf("Error: TileMapv2_ProbeSelection_TwoSample_WithMask, %c not supported in logic expressions!\n", pParam->strPattern[ni]);
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
		printf("Error: TileMapv2_ProbeSelection_TwoSample_WithMask, lack '(' or ')'!\n");
		exit(EXIT_FAILURE);
	}

	if(nLogicLen != 3)
	{
		printf("Error: TileMapv2_ProbeSelection_TwoSample_WithMask, logic expression wrong!\n");
		printf("       Not a two sample comparison!\n");
		exit(EXIT_FAILURE);
	}
	if((vLogic[0] != 'G') || (vLogic[2] != 'G'))
	{
		printf("Error: TileMapv2_ProbeSelection_TwoSample_WithMask, logic expression wrong!\n");
		exit(EXIT_FAILURE);
	}
	if((vLogic[1] != '<') && (vLogic[1] != '>'))
	{
		printf("Error: TileMapv2_ProbeSelection_TwoSample_WithMask, logic expression wrong!\n");
		exit(EXIT_FAILURE);
	}
	if(vLogic[1] == '<')
	{
		ng1 = (int)vGid[2]-1;
		ng2 = (int)vGid[0]-1;
	}
	else if(vLogic[1] == '>')
	{
		ng1 = (int)vGid[0]-1;
		ng2 = (int)vGid[2]-1;
	}
	else
	{
		printf("Error: TileMapv2_ProbeSelection_TwoSample_WithMask, logic expression wrong!\n");
		exit(EXIT_FAILURE);
	}

	if((ng1<0) || (ng1>=pParam->nGroupNum) || (ng2<0) || (ng2>=pParam->nGroupNum))
	{
		printf("Error: TileMapv2_ProbeSelection_TwoSample_WithMask, group id out of range!\n");
		exit(EXIT_FAILURE);
	}
	if( (pParam->pGroupSize->pMatElement[ng1] <= 0) || (pParam->pGroupSize->pMatElement[ng2] <= 0) )
	{
		printf("Error: TileMapv2_ProbeSelection_TwoSample_WithMask, at least one group is empty!\n");
		exit(EXIT_FAILURE);
	}
	
	/* ---- */
	/* init */
	/* ---- */
	nSampleNum = pParam->nSampleNum;
	nGroupNum = pParam->nGroupNum;
	nVarGroupNum = pParam->nVargroupNum;
	pGroupSizeCopy = NULL;
	pGroupSizeCopy = CreateIntMatrix(1, pParam->nGroupNum);
	pDf = NULL;
	pDf = CreateIntMatrix(1, pParam->nVargroupNum);
	
	/* ------------- */
	/* create memory */
	/* ------------- */
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		/* chromosome position, ID = 0 */
		/* raw data, ID = (1:nSampleNum) */
		/* test statistics, ID = nSampleNum+1 */
		nClustId = nSampleNum+1;
		pData->vSeqData[ni]->vData[nClustId] = NULL;
		pData->vSeqData[ni]->vData[nClustId] = CreateDoubleMatrix(1, pData->vSeqData[ni]->nDataNum);
		if( (pData->vSeqData[ni]->nDataNum > 0) && (pData->vSeqData[ni]->vData[nClustId] == NULL))
		{
			printf("Error: TileMapv2_ProbeSelection_TwoSample_WithMask, cannot create space for storing test statistics!\n");
			exit(EXIT_FAILURE);
		}

		/* mask, ID = nSampleNum+1 */
		nClustId = nSampleNum+1;
		pMask->vSeqData[ni]->vData[nClustId] = NULL;
		pMask->vSeqData[ni]->vData[nClustId] = CreateDoubleMatrix(1, pMask->vSeqData[ni]->nDataNum);
		if( (pMask->vSeqData[ni]->nDataNum > 0) && (pMask->vSeqData[ni]->vData[nClustId] == NULL))
		{
			printf("Error: TileMapv2_ProbeSelection_TwoSample_WithMask, cannot create space for storing masks!\n");
			exit(EXIT_FAILURE);
		}

		/* within sample mean, ID = nSampleNum+2+(0:nGroupNum-1) */
		for(nj=0; nj<nGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nj;
			pData->vSeqData[ni]->vData[nClustId] = NULL;
			pData->vSeqData[ni]->vData[nClustId] = CreateDoubleMatrix(1, pData->vSeqData[ni]->nDataNum);
			if( (pData->vSeqData[ni]->nDataNum > 0) && (pData->vSeqData[ni]->vData[nClustId] == NULL))
			{
				printf("Error: TileMapv2_ProbeSelection_TwoSample_WithMask, cannot create space for storing sample mean!\n");
				exit(EXIT_FAILURE);
			}
		}

		/* within sample variance, ID = nSampleNum+2+nGroupNum+(0:nGroupNum-1) */
		for(nj=0; nj<nGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nGroupNum+nj;
			pData->vSeqData[ni]->vData[nClustId] = NULL;
			pData->vSeqData[ni]->vData[nClustId] = CreateDoubleMatrix(1, pData->vSeqData[ni]->nDataNum);
			if( (pData->vSeqData[ni]->nDataNum > 0) && (pData->vSeqData[ni]->vData[nClustId] == NULL))
			{
				printf("Error: TileMapv2_ProbeSelection_TwoSample_WithMask, cannot create space for storing within sample variance!\n");
				exit(EXIT_FAILURE);
			}
		}

		/* sample size, ID = nSampleNum+2+2*nGroupNum+(0:nGroupNum-1) */
		for(nj=0; nj<nGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nGroupNum+nGroupNum+nj;
			pData->vSeqData[ni]->vData[nClustId] = NULL;
			pData->vSeqData[ni]->vData[nClustId] = CreateDoubleMatrix(1, pData->vSeqData[ni]->nDataNum);
			if( (pData->vSeqData[ni]->nDataNum > 0) && (pData->vSeqData[ni]->vData[nClustId] == NULL))
			{
				printf("Error: TileMapv2_ProbeSelection_TwoSample_WithMask, cannot create space for storing sample size!\n");
				exit(EXIT_FAILURE);
			}
		}

		/* pooled variance, ID = nSampleNum+2+3*nGroupNum+(0:nVarGroupNum-1) */
		for(nj=0; nj<nVarGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nGroupNum+nGroupNum+nGroupNum+nj;
			pData->vSeqData[ni]->vData[nClustId] = NULL;
			pData->vSeqData[ni]->vData[nClustId] = CreateDoubleMatrix(1, pData->vSeqData[ni]->nDataNum);
			if( (pData->vSeqData[ni]->nDataNum > 0) && (pData->vSeqData[ni]->vData[nClustId] == NULL))
			{
				printf("Error: TileMapv2_ProbeSelection_TwoSample_WithMask, cannot create space for storing pooled variance!\n");
				exit(EXIT_FAILURE);
			}
		}

		/* pooled degree of freedom, ID = nSampleNum+2+3*nGroupNum+nVarGroupNum+(0:nVarGroupNum-1) */
		for(nj=0; nj<nVarGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nGroupNum+nGroupNum+nGroupNum+nVarGroupNum+nj;
			pData->vSeqData[ni]->vData[nClustId] = NULL;
			pData->vSeqData[ni]->vData[nClustId] = CreateDoubleMatrix(1, pData->vSeqData[ni]->nDataNum);
			if( (pData->vSeqData[ni]->nDataNum > 0) && (pData->vSeqData[ni]->vData[nClustId] == NULL))
			{
				printf("Error: TileMapv2_ProbeSelection_TwoSample_WithMask, cannot create space for storing pooled d.f.!\n");
				exit(EXIT_FAILURE);
			}
		}

		/* fold change, ID = nSampleNum+2+3*nGroupNum+2*nVarGroupNum */
		nClustId = nSampleNum+2+nGroupNum+nGroupNum+nGroupNum+nVarGroupNum+nVarGroupNum;
		pData->vSeqData[ni]->vData[nClustId] = NULL;
		pData->vSeqData[ni]->vData[nClustId] = CreateDoubleMatrix(1, pData->vSeqData[ni]->nDataNum);
		if( (pData->vSeqData[ni]->nDataNum > 0) && (pData->vSeqData[ni]->vData[nClustId] == NULL))
		{
			printf("Error: TileMapv2_ProbeSelection_TwoSample_WithMask, cannot create space for storing fold change!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	/* ---------------- */
	/* get mean         */
	/* ---------------- */

	/* mean */
	for(nj=0; nj<nSampleNum; nj++)
	{
		nClustId = pParam->pGroupLabel->pMatElement[nj]-1;
		if( (nClustId<0) || (nClustId>=nGroupNum) )
			continue;

		pGroupSizeCopy->pMatElement[nClustId] += 1;
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			if(pData->vSeqData[ni]->nDataNum <= 0)
				continue;
				
			vSum = pData->vSeqData[ni]->vData[nSampleNum+2+nClustId]->pMatElement;
			vSize = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+nClustId]->pMatElement;
			vExp = pData->vSeqData[ni]->vData[1+nj]->pMatElement;
			vMask = pMask->vSeqData[ni]->vData[1+nj]->pMatElement;
			for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
			{
				if(vMask[nk] < 0.5)
				{
					vSum[nk] += vExp[nk];
					vSize[nk] += 1.0;
				}
			}
		}
	}

	for(nj=0; nj<nGroupNum; nj++)
	{
		if(pGroupSizeCopy->pMatElement[nj] != pParam->pGroupSize->pMatElement[nj])
		{
			printf("Error: TileMapv2_ProbeSelection_TwoSample_WithMask, class size not match!\n");
			exit(EXIT_FAILURE);
		}

		if(pGroupSizeCopy->pMatElement[nj] > 0)
		{
            for(ni=0; ni<pData->nSeqNum; ni++)
			{
				if(pData->vSeqData[ni]->nDataNum <= 0)
					continue;
				
				vSum = pData->vSeqData[ni]->vData[nSampleNum+2+nj]->pMatElement;
				vSize = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+nj]->pMatElement;
				for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
				{
					if( vSize[nk] > 0.5)
						vSum[nk] /= vSize[nk];
				}
			}
		}
	}

	/* ------------------------------------------------------------------------ */
	/* if each group has only one replicate, compute the difference of the mean */
	/* as the test statistics.                                                  */
	/* ------------------------------------------------------------------------ */
	if( (pParam->pGroupSize->pMatElement[ng1] == 1) && (pParam->pGroupSize->pMatElement[ng2] == 1) )
	{
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			if(pData->vSeqData[ni]->nDataNum <= 0)
				continue;
				
			vSum = pData->vSeqData[ni]->vData[nSampleNum+1]->pMatElement;
			vMask = pMask->vSeqData[ni]->vData[nSampleNum+1]->pMatElement;
			vExp = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+nGroupNum+nVarGroupNum+nVarGroupNum]->pMatElement;
			vMu = pData->vSeqData[ni]->vData[nSampleNum+2+ng1]->pMatElement;
			vMu2 = pData->vSeqData[ni]->vData[nSampleNum+2+ng2]->pMatElement;
			vSize = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+ng1]->pMatElement;
			vSize2 = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+ng2]->pMatElement;
			for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
			{
				if( (vSize[nk] > 0.5) && (vSize2[nk] > 0.5) )
				{
					vSum[nk] = vMu[nk]-vMu2[nk];
					vExp[nk] = vSum[nk];
					vMask[nk] = 0.0;
				}
				else
				{
					vMask[nk] = 1.0;
				}
			}
		}

		/* destroy unnecessary memory and return */
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			/* chromosome position, ID = 0 */
			/* raw data, ID = (1:nSampleNum) */
			/* test statistics, ID = nSampleNum+1 */
		
			/* within sample mean, ID = nSampleNum+2+(0:nGroupNum-1) */
			for(nj=0; nj<nGroupNum; nj++)
			{
				nClustId = nSampleNum+2+nj;
				DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
				pData->vSeqData[ni]->vData[nClustId] = NULL;
			}
			/* within sample variance, ID = nSampleNum+2+nGroupNum+(0:nGroupNum-1) */
			for(nj=0; nj<nGroupNum; nj++)
			{
				nClustId = nSampleNum+2+nGroupNum+nj;
				DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
				pData->vSeqData[ni]->vData[nClustId] = NULL;
			}
			/* sample size, ID = nSampleNum+2+nGroupNum+nGroupNum+(0:nGroupNum-1) */
			for(nj=0; nj<nGroupNum; nj++)
			{
				nClustId = nSampleNum+2+nGroupNum+nGroupNum+nj;
				DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
				pData->vSeqData[ni]->vData[nClustId] = NULL;
			}
			/* pooled variance, ID = nSampleNum+2+3*nGroupNum+(0:nVarGroupNum-1) */
			for(nj=0; nj<nVarGroupNum; nj++)
			{
				nClustId = nSampleNum+2+nGroupNum+nGroupNum+nGroupNum+nj;
				DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
				pData->vSeqData[ni]->vData[nClustId] = NULL;
			}
			/* pooled d.f., ID = nSampleNum+2+3*nGroupNum+nVarGroupNum+(0:nVarGroupNum-1) */
			for(nj=0; nj<nVarGroupNum; nj++)
			{
				nClustId = nSampleNum+2+nGroupNum+nGroupNum+nGroupNum+nVarGroupNum+nj;
				DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
				pData->vSeqData[ni]->vData[nClustId] = NULL;
			}

			nClustId = nSampleNum+2+nGroupNum+nGroupNum+nGroupNum+nVarGroupNum+nVarGroupNum;
			pData->vSeqData[ni]->vData[nSampleNum+2] = pData->vSeqData[ni]->vData[nClustId];
			pData->vSeqData[ni]->vData[nClustId] = NULL;
		}

		pData->pFieldType->pMatElement[nSampleNum+1] = 1;
        pData->pFieldType->pMatElement[nSampleNum+2] = 1;
		pMask->pFieldType->pMatElement[nSampleNum+1] = 7;
			
		DestroyIntMatrix(pGroupSizeCopy);
		DestroyIntMatrix(pDf);

		/* -------------- */
		/* save data      */
		/* -------------- */
		pCol = NULL;
		pCol = CreateIntMatrix(1,pData->nColNum);
		if(pCol == NULL)
		{
			printf("Error: TileMapv2_ProbeSelection_TwoSample_WithMask, cannot create output column information!\n");
			exit(EXIT_FAILURE);
		}
		
		printf("Exporting probe level test statistics...\n");
		pCol->pMatElement[0] = 1;
		pCol->pMatElement[nSampleNum+1] = 1;
		Affy_SaveBAR_Columns_Fast(strOutFile, pData, pCol);

		printf("Exporting fold changes...\n");
		pCol->pMatElement[nSampleNum+1] = 0;
		pCol->pMatElement[nSampleNum+2] = 1;
		Affy_SaveBAR_Columns_Fast(strFCFile, pData, pCol);

		DestroyIntMatrix(pCol);

		pCol = NULL;
		pCol = CreateIntMatrix(1,pMask->nColNum);
		if(pCol == NULL)
		{
			printf("Error: TileMapv2_ProbeSelection_TwoSample_WithMask, cannot create output column information!\n");
			exit(EXIT_FAILURE);
		}
		
		printf("Exporting masks...\n");
		pCol->pMatElement[0] = 1;
		pCol->pMatElement[nSampleNum+1] = 1;
		sprintf(strLine, "%s.mask", strOutFile);
		Affy_SaveBAR_Columns_Fast(strLine, pMask, pCol);

		DestroyIntMatrix(pCol);

		/* -------------- */
		/* release memory */
		/* -------------- */
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			/* chromosome position, ID = 0 */
			/* raw data, ID = (1:nSampleNum) */
			/* test statistics, ID = nSampleNum+1 */
			nClustId = nSampleNum+1;
			DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
			pData->vSeqData[ni]->vData[nClustId] = NULL;

			nClustId = nSampleNum+2;
			DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
			pData->vSeqData[ni]->vData[nClustId] = NULL;
		}
		pData->pFieldType->pMatElement[nSampleNum+1] = 0;
		pData->pFieldType->pMatElement[nSampleNum+2] = 0;

		for(ni=0; ni<pMask->nSeqNum; ni++)
		{
			nClustId = nSampleNum+1;
			DestroyDoubleMatrix(pMask->vSeqData[ni]->vData[nClustId]);
			pMask->vSeqData[ni]->vData[nClustId] = NULL;
		}
		pMask->pFieldType->pMatElement[nSampleNum+1] = 0;

		/* return */
		return PROC_SUCCESS;
	}
		
	/* --------------------------------------------- */
	/* if one has extra d.f., compute the variance.  */
	/* --------------------------------------------- */
	for(nj=0; nj<nSampleNum; nj++)
	{
		nClustId = pParam->pGroupLabel->pMatElement[nj]-1;
		if( (nClustId<0) || (nClustId>=nGroupNum) )
			continue;
		
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			if(pData->vSeqData[ni]->nDataNum <= 0)
				continue;
				
			vSize = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+nClustId]->pMatElement;
			vSum = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nClustId]->pMatElement;
			vAve = pData->vSeqData[ni]->vData[nSampleNum+2+nClustId]->pMatElement;
			vExp = pData->vSeqData[ni]->vData[1+nj]->pMatElement;
			vMask = pMask->vSeqData[ni]->vData[1+nj]->pMatElement;

			for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
			{
				if(vMask[nk] < 0.5)
				{
					dTemp = (vExp[nk]-vAve[nk]);		
					vSum[nk] += dTemp*dTemp;
				}
			}
		}
	}

	/* variance: estimates */
	for(nj=0; nj<nVarGroupNum; nj++)
	{
		if(pParam->vVargroupMap[nj]->nWidth != pParam->pVargroupSize->pMatElement[nj])
		{
			printf("Error: TileMapv2_ProbeSelection_TwoSample_WithMask, variance group number not match!\n");
			exit(EXIT_FAILURE);
		}
	}

	for(nj=0; nj<nVarGroupNum; nj++)
	{
		for(nl=0; nl<pParam->vVargroupMap[nj]->nWidth; nl++)
		{
			nClustId = pParam->vVargroupMap[nj]->pMatElement[nl]-1;
			if(pParam->pGroupSize->pMatElement[nClustId] > 0)
				pDf->pMatElement[nj] += pParam->pGroupSize->pMatElement[nClustId]-1;
			
			for(ni=0; ni<pData->nSeqNum; ni++)
			{
				if(pData->vSeqData[ni]->nDataNum <= 0)
					continue;
				
				vSum = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+nGroupNum+nj]->pMatElement;
				vSize2 = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+nGroupNum+nVarGroupNum+nj]->pMatElement;
				vExp = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nClustId]->pMatElement;
				vSize = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+nClustId]->pMatElement;
				for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
				{
					if(vSize[nk] > 0.5)
					{
						vSum[nk] += vExp[nk];
						vSize2[nk] += vSize[nk]-1.0;
					}
				}
			}
		}
	}

	for(nj=0; nj<nVarGroupNum; nj++)
	{
		for(nl=0; nl<pParam->vVargroupMap[nj]->nWidth; nl++)
		{
			nClustId = pParam->vVargroupMap[nj]->pMatElement[nl]-1;
			if( (nClustId == ng1) || (nClustId == ng2) )
			{
				if(pDf->pMatElement[nj] <= 0)
				{
					printf("Error: TileMapv2_ProbeSelection_TwoSample_WithMask, group %d does not have enough d.f. to estimate variance!\n", (nClustId+1));
					exit(EXIT_FAILURE);
				}
			}
		}
	}

	/* ------------- */
	/* shrinking var */
	/* ------------- */
	for(nj=0; nj<nVarGroupNum; nj++)
	{
		nClustId = nSampleNum+2+nGroupNum+nGroupNum+nGroupNum+nj;
		nDfId = nSampleNum+2+nGroupNum+nGroupNum+nGroupNum+nVarGroupNum+nj;
		TileMapv2_GetShrinkingVar_FromSSR(pData, pDf->pMatElement[nj], nClustId, nDfId);
	}

	/* destroy old variances */
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		for(nj=0; nj<pParam->nGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nGroupNum+nj;
			DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
			pData->vSeqData[ni]->vData[nClustId] = NULL;
		}
	}


	/* assign new variances */
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		for(nj=0; nj<nVarGroupNum; nj++)
		{
			for(nl=0; nl<pParam->vVargroupMap[nj]->nWidth; nl++)
			{
				nClustId = pParam->vVargroupMap[nj]->pMatElement[nl]-1;
				pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nClustId] = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+nGroupNum+nj];
			}
		}
	}

	/* ----------- */
	/* simulation  */
	/* ----------- */
	/* compute probeset level summary */
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		if(pData->vSeqData[ni]->nDataNum <= 0)
			continue;
				
		vMu = pData->vSeqData[ni]->vData[nSampleNum+2+ng1]->pMatElement;
		vSigma = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+ng1]->pMatElement;
		vSize = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+ng1]->pMatElement;
		vMu2 = pData->vSeqData[ni]->vData[nSampleNum+2+ng2]->pMatElement;
		vSigma2 = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+ng2]->pMatElement;
		vSize2 = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+ng2]->pMatElement;
		vSum = pData->vSeqData[ni]->vData[nSampleNum+1]->pMatElement;
		vMask = pMask->vSeqData[ni]->vData[nSampleNum+1]->pMatElement;
		vExp = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+nGroupNum+nVarGroupNum+nVarGroupNum]->pMatElement;
		
		for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
		{
			if( (vSize[nk]<0.5) || (vSize2[nk]<0.5) || (vSigma[nk]<=0.0) || (vSigma2[nk]<=0.0) )
			{
				vMask[nk] = 1.0;
			}
			else
			{
				vMask[nk] = 0.0;
				vExp[nk] = vMu[nk]-vMu2[nk];
				dVarTemp = sqrt(vSigma[nk]/vSize[nk]+vSigma2[nk]/vSize2[nk]);
				vSum[nk] = (vMu[nk]-vMu2[nk])/dVarTemp;
			}
		}
	}
	
	/* -------------- */
	/* release memory */
	/* -------------- */
	/* destroy unnecessary memory and return */
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		/* chromosome position, ID = 0 */
		/* raw data, ID = (1:nSampleNum) */
		/* test statistics, ID = nSampleNum+1 */
	
		/* within sample mean, ID = nSampleNum+2+(0:nGroupNum-1) */
		for(nj=0; nj<nGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nj;
			DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
			pData->vSeqData[ni]->vData[nClustId] = NULL;
		}
		/* within sample variance, ID = nSampleNum+2+nGroupNum+(0:nGroupNum-1) */
		for(nj=0; nj<nGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nGroupNum+nj;
			pData->vSeqData[ni]->vData[nClustId] = NULL;
		}
		/* sample size, ID = nSampleNum+2+nGroupNum+nGroupNum+(0:nGroupNum-1) */
		for(nj=0; nj<nGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nGroupNum+nGroupNum+nj;
			DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
			pData->vSeqData[ni]->vData[nClustId] = NULL;
		}
		/* pooled variance, ID = nSampleNum+2+3*nGroupNum+(0:nVarGroupNum-1) */
		for(nj=0; nj<nVarGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nGroupNum+nGroupNum+nGroupNum+nj;
			DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
			pData->vSeqData[ni]->vData[nClustId] = NULL;
		}
		/* pooled d.f., ID = nSampleNum+2+3*nGroupNum+nVarGroupNum+(0:nVarGroupNum-1) */
		for(nj=0; nj<nVarGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nGroupNum+nGroupNum+nGroupNum+nVarGroupNum+nj;
			DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
			pData->vSeqData[ni]->vData[nClustId] = NULL;
		}

		nClustId = nSampleNum+2+nGroupNum+nGroupNum+nGroupNum+nVarGroupNum+nVarGroupNum;
		pData->vSeqData[ni]->vData[nSampleNum+2] = pData->vSeqData[ni]->vData[nClustId];
		pData->vSeqData[ni]->vData[nClustId] = NULL;
	}

	pData->pFieldType->pMatElement[nSampleNum+1] = 1;
	pData->pFieldType->pMatElement[nSampleNum+2] = 1;
	pMask->pFieldType->pMatElement[nSampleNum+1] = 7;

	DestroyIntMatrix(pGroupSizeCopy);
	DestroyIntMatrix(pDf);

	/* -------------- */
	/* save data      */
	/* -------------- */
	pCol = NULL;
	pCol = CreateIntMatrix(1,pData->nColNum);
	if(pCol == NULL)
	{
		printf("Error: TileMapv2_ProbeSelection_TwoSample_WithMask, cannot create output column information!\n");
		exit(EXIT_FAILURE);
	}
	
	printf("Exporting probe level test statistics...\n");
	pCol->pMatElement[0] = 1;
	pCol->pMatElement[nSampleNum+1] = 1;
	Affy_SaveBAR_Columns_Fast(strOutFile, pData, pCol);

    printf("Exporting fold changes...\n");
	pCol->pMatElement[nSampleNum+1] = 0;
	pCol->pMatElement[nSampleNum+2] = 1;
	Affy_SaveBAR_Columns_Fast(strFCFile, pData, pCol);

	DestroyIntMatrix(pCol);

	pCol = NULL;
	pCol = CreateIntMatrix(1,pMask->nColNum);
	if(pCol == NULL)
	{
		printf("Error: TileMapv2_ProbeSelection_TwoSample_WithMask, cannot create output column information!\n");
		exit(EXIT_FAILURE);
	}
	
	printf("Exporting masks...\n");
	pCol->pMatElement[0] = 1;
	pCol->pMatElement[nSampleNum+1] = 1;
	sprintf(strLine, "%s.mask", strOutFile);
	Affy_SaveBAR_Columns_Fast(strLine, pMask, pCol);

	DestroyIntMatrix(pCol);

	/* -------------- */
	/* release memory */
	/* -------------- */
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		/* chromosome position, ID = 0 */
		/* raw data, ID = (1:nSampleNum) */
		/* test statistics, ID = nSampleNum+1 */
		nClustId = nSampleNum+1;
		DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
		pData->vSeqData[ni]->vData[nClustId] = NULL;

		nClustId = nSampleNum+2;
		DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
		pData->vSeqData[ni]->vData[nClustId] = NULL;
	}
	pData->pFieldType->pMatElement[nSampleNum+1] = 0;
	pData->pFieldType->pMatElement[nSampleNum+2] = 0;

	for(ni=0; ni<pMask->nSeqNum; ni++)
	{
		nClustId = nSampleNum+1;
		DestroyDoubleMatrix(pMask->vSeqData[ni]->vData[nClustId]);
		pMask->vSeqData[ni]->vData[nClustId] = NULL;
	}
	pMask->pFieldType->pMatElement[nSampleNum+1] = 0;

	/* ------ */
	/* return */
	/* ------ */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ProbeSelection_MultiSample()                                 */
/*  TileMapv2 probe level summary: multiple sample comparisons.            */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_ProbeSelection_MultiSample(struct tagBARData *pData,
						struct tagTileMapv2Param *pParam, char strOutFile[],
						char strFCFile[])
{
	/* ------ */
	/* define */
	/* ------ */
	int nLogicLen;
	int ni,nj,nk,nl,nLen,nIter;
	int nLeftNum,nRightNum,nTNumLen;
	char vLogic[MED_LINE_LENGTH];
	char vTNumber[LINE_LENGTH];
	double vGid[MED_LINE_LENGTH];
	int ng,nAvailableDF,nHit;
	struct DOUBLEMATRIX *vVid[MED_LINE_LENGTH];
	struct BYTEMATRIX *vLid[MED_LINE_LENGTH];
	int nTrue;
	struct BYTEMATRIX *pTrueVec;
	struct DOUBLEMATRIX *pTrueVal;

	int nGroupNum;
	int nSampleNum;
	int nVarGroupNum;
	struct INTMATRIX *pGroupSizeCopy;
	struct INTMATRIX *pDf;

	int nClustId;
	double *vExp,*vSum,*vAve,*vSigma,*vMu;
	unsigned char *vEval;
	double dDenom;
	double dTemp;
	double dVarmean,dVarsst;
	double dB,dV2,dK,dN;
	int nTotalProbeNum;
	struct DOUBLEMATRIX *pSigCoef;
	struct INTMATRIX *pCol;
	double dLowLimit,dHighLimit;
	
	
	/* ------------- */
	/* expression    */
	/* ------------- */
	nLogicLen = 0;
	if(pParam->strPattern[0] == '\0')
	{
		printf("Error: TileMapv2_ProbeSelection_MultiSample, no logic expression formula available!\n");
		exit(EXIT_FAILURE);
	}
	nLen = (int)strlen(pParam->strPattern);
	nj = 0;
	nLeftNum = 0;
	nRightNum = 0;
	nTNumLen = 0;

	for(ni=0; ni<nLen; ni++)
	{
		if( pParam->strPattern[ni] == '(' )
		{
			if(nTNumLen != 0)
			{
				printf("Error: TileMapv2_ProbeSelection_MultiSample, logic expression wrong!\n");
				exit(EXIT_FAILURE);
			}

			vLogic[nj] = pParam->strPattern[ni];
			vGid[nj] = -1.0;
			nj++;
			nLeftNum++;
		}

		else if(pParam->strPattern[ni] == ')')
		{
			if(nTNumLen > 0)
			{
				vTNumber[nTNumLen] = '\0';
				vLogic[nj] = 'G';
				vGid[nj] = atof(vTNumber);
				nTNumLen = 0;
				nj++;
			}
			
			vLogic[nj] = pParam->strPattern[ni];
			vGid[nj] = -1.0;
			nj++;
			nRightNum++;
		}
		else if( (pParam->strPattern[ni] == '<') || (pParam->strPattern[ni] == '>') 
			|| (pParam->strPattern[ni] == '&') || (pParam->strPattern[ni] == '|') )
		{
			if(nTNumLen > 0)
			{
				vTNumber[nTNumLen] = '\0';
				vLogic[nj] = 'G';
				vGid[nj] = atof(vTNumber);
				nTNumLen = 0;
				nj++;
			}
			vLogic[nj] = pParam->strPattern[ni];
			vGid[nj] = -1.0;
			nj++;	
		}

		else if( (pParam->strPattern[ni] >= '0') && (pParam->strPattern[ni] <= '9') )
		{
			vTNumber[nTNumLen] = pParam->strPattern[ni];
			nTNumLen++;
		}

		else if( (pParam->strPattern[ni] == ' ') || (pParam->strPattern[ni] == '\t') )
		{
			/* do nothing */
		}
		else
		{
			printf("Error: TileMapv2_ProbeSelection_MultiSample, %c not supported in logic expressions!\n", pParam->strPattern[ni]);
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
		printf("Error: TileMapv2_ProbeSelection_MultiSample, lack '(' or ')'!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nLogicLen; ni++)
	{
		if(vGid[ni] > 0.5)
		{
			ng = (int)(vGid[ni])-1;
	
			for(nj=0; nj<pParam->nVargroupNum; nj++)
			{
				nAvailableDF = 0;
				nHit = 0;
				for(nl=0; nl<pParam->vVargroupMap[nj]->nWidth; nl++)
				{
					nClustId = pParam->vVargroupMap[nj]->pMatElement[nl]-1;
					
					if(pParam->pGroupSize->pMatElement[nClustId]>0)
						nAvailableDF += (pParam->pGroupSize->pMatElement[nClustId]-1);

					if(nClustId == ng)
					{
						nHit = 1;
					}
				}

				if(nHit == 1)
				{
					if(nAvailableDF <= 0)
					{
						printf("Error: TileMapv2_ProbeSelection_MultiSample, group %d does not have enough d.f. to estimate variance!\n", (nClustId+1));
						exit(EXIT_FAILURE);
					}
					break;
				}
			}
		}
	}

	/* ---- */
	/* init */
	/* ---- */
	nSampleNum = pParam->nSampleNum;
	nGroupNum = pParam->nGroupNum;
	nVarGroupNum = pParam->nVargroupNum;
	pGroupSizeCopy = NULL;
	pGroupSizeCopy = CreateIntMatrix(1, pParam->nGroupNum);
	pDf = NULL;
	pDf = CreateIntMatrix(1, pParam->nVargroupNum);

	/* ------------- */
	/* create memory */
	/* ------------- */
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		/* chromosome position, ID = 0 */
		/* raw data, ID = (1:nSampleNum) */
		/* test statistics, ID = nSampleNum+1 */
		nClustId = nSampleNum+1;
		pData->vSeqData[ni]->vData[nClustId] = NULL;
		pData->vSeqData[ni]->vData[nClustId] = CreateDoubleMatrix(1, pData->vSeqData[ni]->nDataNum);
		if( (pData->vSeqData[ni]->nDataNum > 0) && (pData->vSeqData[ni]->vData[nClustId] == NULL))
		{
			printf("Error: TileMapv2_ProbeSelection_MultiSample, cannot create space for storing test statistics!\n");
			exit(EXIT_FAILURE);
		}

		/* within sample mean, ID = nSampleNum+2+(0:nGroupNum-1) */
		for(nj=0; nj<nGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nj;
			pData->vSeqData[ni]->vData[nClustId] = NULL;
			pData->vSeqData[ni]->vData[nClustId] = CreateDoubleMatrix(1, pData->vSeqData[ni]->nDataNum);
			if( (pData->vSeqData[ni]->nDataNum > 0) && (pData->vSeqData[ni]->vData[nClustId] == NULL))
			{
				printf("Error: TileMapv2_ProbeSelection_MultiSample, cannot create space for storing sample mean!\n");
				exit(EXIT_FAILURE);
			}
		}

		/* within sample variance, ID = nSampleNum+2+nGroupNum+(0:nGroupNum-1) */
		for(nj=0; nj<nGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nGroupNum+nj;
			pData->vSeqData[ni]->vData[nClustId] = NULL;
			pData->vSeqData[ni]->vData[nClustId] = CreateDoubleMatrix(1, pData->vSeqData[ni]->nDataNum);
			if( (pData->vSeqData[ni]->nDataNum > 0) && (pData->vSeqData[ni]->vData[nClustId] == NULL))
			{
				printf("Error: TileMapv2_ProbeSelection_MultiSample, cannot create space for storing within sample variance!\n");
				exit(EXIT_FAILURE);
			}
		}

		/* pooled variance, ID = nSampleNum+2+2*nGroupNum+(0:nVarGroupNum-1) */
		for(nj=0; nj<nVarGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nGroupNum+nGroupNum+nj;
			pData->vSeqData[ni]->vData[nClustId] = NULL;
			pData->vSeqData[ni]->vData[nClustId] = CreateDoubleMatrix(1, pData->vSeqData[ni]->nDataNum);
			if( (pData->vSeqData[ni]->nDataNum > 0) && (pData->vSeqData[ni]->vData[nClustId] == NULL))
			{
				printf("Error: TileMapv2_ProbeSelection_MultiSample, cannot create space for storing pooled variance!\n");
				exit(EXIT_FAILURE);
			}
		}

		/* sampled means, ID = nSampleNum+2+2*nGroupNum+nVarGroupNum+(0:nGroupNum-1) */
		for(nj=0; nj<nGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nGroupNum+nGroupNum+nVarGroupNum+nj;
			pData->vSeqData[ni]->vData[nClustId] = NULL;
			pData->vSeqData[ni]->vData[nClustId] = CreateDoubleMatrix(1, pData->vSeqData[ni]->nDataNum);
			if( (pData->vSeqData[ni]->nDataNum > 0) && (pData->vSeqData[ni]->vData[nClustId] == NULL))
			{
				printf("Error: TileMapv2_ProbeSelection_MultiSample, cannot create space for simulation!\n");
				exit(EXIT_FAILURE);
			}
		}
	}
	
	/* ---------------- */
	/* get mean         */
	/* ---------------- */

	/* mean */
	for(nj=0; nj<nSampleNum; nj++)
	{
		nClustId = pParam->pGroupLabel->pMatElement[nj]-1;
		if( (nClustId<0) || (nClustId>=nGroupNum) )
			continue;

		pGroupSizeCopy->pMatElement[nClustId] += 1;
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			if(pData->vSeqData[ni]->nDataNum <= 0)
				continue;
			
			vSum = pData->vSeqData[ni]->vData[nSampleNum+2+nClustId]->pMatElement;
			vExp = pData->vSeqData[ni]->vData[1+nj]->pMatElement;
			for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
			{
				vSum[nk] += vExp[nk];
			}
		}
	}

	for(nj=0; nj<nGroupNum; nj++)
	{
		if(pGroupSizeCopy->pMatElement[nj] != pParam->pGroupSize->pMatElement[nj])
		{
			printf("Error: TileMapv2_ProbeSelection_MultiSample, class size not match!\n");
			exit(EXIT_FAILURE);
		}

		if(pGroupSizeCopy->pMatElement[nj] > 0)
		{
            dDenom = (double)(pGroupSizeCopy->pMatElement[nj]);
			for(ni=0; ni<pData->nSeqNum; ni++)
			{
				if(pData->vSeqData[ni]->nDataNum <= 0)
					continue;
				
				vSum = pData->vSeqData[ni]->vData[nSampleNum+2+nj]->pMatElement;
				for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
				{
					vSum[nk] /= dDenom;
				}
			}
		}
	}

	/* --------------------------------------------- */
	/* if one has extra d.f., compute the variance.  */
	/* --------------------------------------------- */
	for(nj=0; nj<nSampleNum; nj++)
	{
		nClustId = pParam->pGroupLabel->pMatElement[nj]-1;
		if( (nClustId<0) || (nClustId>=nGroupNum) )
			continue;
		
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			if(pData->vSeqData[ni]->nDataNum <= 0)
				continue;
			
			vSum = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nClustId]->pMatElement;
			vAve = pData->vSeqData[ni]->vData[nSampleNum+2+nClustId]->pMatElement;
			vExp = pData->vSeqData[ni]->vData[1+nj]->pMatElement;
						
			for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
			{
				dTemp = (vExp[nk]-vAve[nk]);		
				vSum[nk] += dTemp*dTemp;
			}
		}
	}

	/* variance: estimates */
	for(nj=0; nj<nVarGroupNum; nj++)
	{
		if(pParam->vVargroupMap[nj]->nWidth != pParam->pVargroupSize->pMatElement[nj])
		{
			printf("Error: TileMapv2_ProbeSelection_MultiSample, variance group number not match!\n");
			exit(EXIT_FAILURE);
		}
	}

	for(nj=0; nj<nVarGroupNum; nj++)
	{
		for(nl=0; nl<pParam->vVargroupMap[nj]->nWidth; nl++)
		{
			nClustId = pParam->vVargroupMap[nj]->pMatElement[nl]-1;
			if(pParam->pGroupSize->pMatElement[nClustId] > 0)
				pDf->pMatElement[nj] += pParam->pGroupSize->pMatElement[nClustId]-1;
			
			for(ni=0; ni<pData->nSeqNum; ni++)
			{
				if(pData->vSeqData[ni]->nDataNum <= 0)
					continue;
				
				vSum = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+nj]->pMatElement;
				vExp = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nClustId]->pMatElement;
				for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
				{
					vSum[nk] += vExp[nk];
				}
			}
		}
	}

	/* ------------- */
	/* shrinking var */
	/* ------------- */
	for(nj=0; nj<nVarGroupNum; nj++)
	{
		dVarmean = 0.0;
		dVarsst = 0.0;
		dB = 0.0;

		if(pDf->pMatElement[nj] > 0)
		{
			dDenom = (double)(pDf->pMatElement[nj]);
			nTotalProbeNum = 0;
			for(ni=0; ni<pData->nSeqNum; ni++)
			{
				if(pData->vSeqData[ni]->nDataNum <= 0)
					continue;
				
				nTotalProbeNum += pData->vSeqData[ni]->nDataNum;
                vSum = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+nj]->pMatElement;
			
				/* xbar */
				for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
				{
					vSum[nk] /= dDenom;
					dVarmean += vSum[nk];
				}
			}
			
			dVarmean /= (double)nTotalProbeNum;

			/* sst */
			for(ni=0; ni<pData->nSeqNum; ni++)
			{
				if(pData->vSeqData[ni]->nDataNum <= 0)
					continue;
				
				vSum = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+nj]->pMatElement;
			
				/* xbar */
				for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
				{
					dTemp = vSum[nk]-dVarmean;
					dVarsst += dTemp*dTemp;
				}
			}

			/* shrinkage factor */
			dV2 = 2.0/dDenom;
			dK = (double)nTotalProbeNum;
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
			for(ni=0; ni<pData->nSeqNum; ni++)
			{
				if(pData->vSeqData[ni]->nDataNum <= 0)
					continue;
				
				vSum = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+nj]->pMatElement;
			
				/* xbar */
				for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
				{
					vSum[nk] = (1-dB)*vSum[nk]+dB*dVarmean+1e-16;
					vSum[nk] = sqrt(vSum[nk]) + 1e-16;
				}
			}
		}
	}

	/* destroy old variances */
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		for(nj=0; nj<pParam->nGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nGroupNum+nj;
			DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
			pData->vSeqData[ni]->vData[nClustId] = NULL;
		}
	}


	/* assign new variances */
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		for(nj=0; nj<nVarGroupNum; nj++)
		{
			for(nl=0; nl<pParam->vVargroupMap[nj]->nWidth; nl++)
			{
				nClustId = pParam->vVargroupMap[nj]->pMatElement[nl]-1;
				pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nClustId] = pData->vSeqData[ni]->vData[nSampleNum+2+nGroupNum+nGroupNum+nj];
			}
		}
	}

	/* ----------- */
	/* simulation  */
	/* ----------- */
	pSigCoef = NULL;
	pSigCoef = CreateDoubleMatrix(1, pParam->nGroupNum);
	vExp = pSigCoef->pMatElement;
	for(ni=0; ni<nGroupNum; ni++)
	{
		if(pParam->pGroupSize->pMatElement[ni] > 0)
		{
			dDenom = (double)(pParam->pGroupSize->pMatElement[ni]);
			vExp[ni] = sqrt(1.0/dDenom);
		}
		else
		{
			vExp[ni] = 0.0;
		}
	}

	/* compute probeset level summary */
	for(nIter=0; nIter<pParam->nMonteCarloNum; nIter++)
	{
		if(nIter%100 == 0)
		{
			printf("iter %d...\n", nIter);
		}

		/* sample probeset by probeset */
		for(nj=0; nj<nGroupNum; nj++)
		{
			dTemp = normrnd(0.0, 1.0);
			for(ni=0; ni<pData->nSeqNum; ni++)
			{
				if(pData->vSeqData[ni]->nDataNum <= 0)
					continue;
				
				nClustId = nSampleNum+2+nGroupNum+nGroupNum+nVarGroupNum+nj;
				vAve = pData->vSeqData[ni]->vData[nClustId]->pMatElement;
				nClustId = nSampleNum+2+nj;
				vMu = pData->vSeqData[ni]->vData[nClustId]->pMatElement;
				nClustId = nSampleNum+2+nGroupNum+nj;
				vSigma = pData->vSeqData[ni]->vData[nClustId]->pMatElement;
				for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
				{
					/* simulate */
					vAve[nk] = vMu[nk]+vExp[nj]*vSigma[nk]*dTemp;
				}
			}
		}

		/* evaluate probeset by probeset */
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			/* attach labels */
			for(nj=0; nj<nLogicLen; nj++)
			{
				if(vGid[nj] > 0.5)
				{
					nk = (int)(vGid[nj]-1.0);
					nClustId = nSampleNum+2+nGroupNum+nGroupNum+nVarGroupNum+nk;
					vVid[nj] = pData->vSeqData[ni]->vData[nClustId];
				}
				else
				{
					vVid[nj] = NULL;
				}
				vLid[nj] = NULL;
			}

			/* logic evaluation */
			pTrueVec = NULL;
			pTrueVal = NULL;
			nTrue = Expression_GeneSelection_EvaluateVec(pData->vSeqData[ni]->nDataNum, vLogic, vVid, vLid, nLogicLen, 0, &pTrueVec, &pTrueVal);
			if(nTrue != 4)
			{
				printf("Error: logic evaluation wrong!\n");
				exit(EXIT_FAILURE);
			}

			/* add */
			vSum = pData->vSeqData[ni]->vData[nSampleNum+1]->pMatElement;
			vEval = pTrueVec->pMatElement;
			for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
			{
				vSum[nk] += (double)(vEval[nk]);
			}
			DestroyByteMatrix(pTrueVec);
		}
	}
	
	DestroyDoubleMatrix(pSigCoef);

	/* normalize */
	dDenom = (double)pParam->nMonteCarloNum;
	dLowLimit = 0.5/dDenom;
	dHighLimit = 1.0-dLowLimit;
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		if(pData->vSeqData[ni]->nDataNum <= 0)
			continue;

		nClustId = nSampleNum+1;
		vSum = pData->vSeqData[ni]->vData[nClustId]->pMatElement;
		for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
		{
			/* simulate */
			vSum[nk] = vSum[nk]/dDenom;

			if(vSum[nk] < dLowLimit)
			{
				vSum[nk] = dLowLimit;
			}
			else if(vSum[nk] > dHighLimit)
			{
				vSum[nk] = dHighLimit;
			}
			vSum[nk] = log(vSum[nk]/(1.0-vSum[nk]));
		}
	}

	/* -------------- */
	/* release memory */
	/* -------------- */
	/* destroy unnecessary memory and return */
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		/* chromosome position, ID = 0 */
		/* raw data, ID = (1:nSampleNum) */
		/* test statistics, ID = nSampleNum+1 */
	
		/* within sample mean, ID = nSampleNum+2+(0:nGroupNum-1) */
		for(nj=0; nj<nGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nj;
			DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
			pData->vSeqData[ni]->vData[nClustId] = NULL;
		}
		/* within sample variance, ID = nSampleNum+2+nGroupNum+(0:nGroupNum-1) */
		for(nj=0; nj<nGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nGroupNum+nj;
			pData->vSeqData[ni]->vData[nClustId] = NULL;
		}
		/* pooled variance, ID = nSampleNum+2+2*nGroupNum+(0:nVarGroupNum-1) */
		for(nj=0; nj<nVarGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nGroupNum+nGroupNum+nj;
			DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
			pData->vSeqData[ni]->vData[nClustId] = NULL;
		}

		/* sampled means, ID = nSampleNum+2+2*nGroupNum+nVarGroupNum+(0:nGroupNum-1) */
		for(nj=0; nj<nGroupNum; nj++)
		{
			nClustId = nSampleNum+2+nGroupNum+nGroupNum+nVarGroupNum+nj;
			DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
			pData->vSeqData[ni]->vData[nClustId] = NULL;
		}
	}

	pData->pFieldType->pMatElement[nSampleNum+1] = 1;
	DestroyIntMatrix(pGroupSizeCopy);
	DestroyIntMatrix(pDf);

	/* -------------- */
	/* save data      */
	/* -------------- */
	pCol = NULL;
	pCol = CreateIntMatrix(1,pData->nColNum);
	if(pCol == NULL)
	{
		printf("Error: TileMapv2_ProbeSelection_MultiSample, cannot create output column information!\n");
		exit(EXIT_FAILURE);
	}
	
	printf("Exporting probe level test statistics...\n");
	pCol->pMatElement[0] = 1;
	pCol->pMatElement[nSampleNum+1] = 1;
	Affy_SaveBAR_Columns_Fast(strOutFile, pData, pCol);

	DestroyIntMatrix(pCol);

	/* -------------- */
	/* release memory */
	/* -------------- */
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		/* chromosome position, ID = 0 */
		/* raw data, ID = (1:nSampleNum) */
		/* test statistics, ID = nSampleNum+1 */
		nClustId = nSampleNum+1;
		DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
		pData->vSeqData[ni]->vData[nClustId] = NULL;
	}
	pData->pFieldType->pMatElement[nSampleNum+1] = 0;

	/* ------ */
	/* return */
	/* ------ */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ProbeSelection_MultiSample_Fast()                            */
/*  TileMapv2 probe level summary: multiple sample comparisons.            */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_ProbeSelection_MultiSample_Fast(struct tagBARData *pData,
						struct tagTileMapv2Param *pParam, char strOutFile[],
						char strFCFile[])
{
	/* ------ */
	/* define */
	/* ------ */
	int nLogicLen;
	int ni,nj,nk,nl,nLen;
	int nLeftNum,nRightNum,nTNumLen;
	char vLogic[MED_LINE_LENGTH];
	char vTNumber[LINE_LENGTH];
	double vGid[MED_LINE_LENGTH];
	int ng,nAvailableDF,nHit;
	struct DOUBLEMATRIX *vVid[MED_LINE_LENGTH];
	struct BYTEMATRIX *vLid[MED_LINE_LENGTH];
	int nTrue;
	struct BYTEMATRIX *pTrueVec;
	struct DOUBLEMATRIX *pTrueVal;

	int nProbeNum;
	int nGroupNum;
	int nSampleNum;
	int nVarGroupNum;
	struct INTMATRIX *pGroupSizeCopy;
	struct INTMATRIX *pDf;
	struct DOUBLEMATRIX *pScore;
	struct DOUBLEMATRIX **vMeans;
	struct DOUBLEMATRIX **vVars;
	struct DOUBLEMATRIX **vSDs;
	struct DOUBLEMATRIX **vMeanDraws;
	
	int nClustId;
	double *vExp,*vSum,*vAve,*vSigma,*vMu;
	unsigned char *vEval;
	double dDenom;
	double dTemp;
	double dVarmean,dVarsst;
	double dB,dV2,dK,dN;
	int nTotalProbeNum;
	struct DOUBLEMATRIX *pSigCoef;
	struct INTMATRIX *pCol;
	double dLowLimit, dHighLimit;
	
	/* ------------- */
	/* expression    */
	/* ------------- */
	nLogicLen = 0;
	if(pParam->strPattern[0] == '\0')
	{
		printf("Error: TileMapv2_ProbeSelection_MultiSample, no logic expression formula available!\n");
		exit(EXIT_FAILURE);
	}
	nLen = (int)strlen(pParam->strPattern);
	nj = 0;
	nLeftNum = 0;
	nRightNum = 0;
	nTNumLen = 0;

	for(ni=0; ni<nLen; ni++)
	{
		if( pParam->strPattern[ni] == '(' )
		{
			if(nTNumLen != 0)
			{
				printf("Error: TileMapv2_ProbeSelection_MultiSample, logic expression wrong!\n");
				exit(EXIT_FAILURE);
			}

			vLogic[nj] = pParam->strPattern[ni];
			vGid[nj] = -1.0;
			nj++;
			nLeftNum++;
		}

		else if(pParam->strPattern[ni] == ')')
		{
			if(nTNumLen > 0)
			{
				vTNumber[nTNumLen] = '\0';
				vLogic[nj] = 'G';
				vGid[nj] = atof(vTNumber);
				nTNumLen = 0;
				nj++;
			}
			
			vLogic[nj] = pParam->strPattern[ni];
			vGid[nj] = -1.0;
			nj++;
			nRightNum++;
		}
		else if( (pParam->strPattern[ni] == '<') || (pParam->strPattern[ni] == '>') 
			|| (pParam->strPattern[ni] == '&') || (pParam->strPattern[ni] == '|') )
		{
			if(nTNumLen > 0)
			{
				vTNumber[nTNumLen] = '\0';
				vLogic[nj] = 'G';
				vGid[nj] = atof(vTNumber);
				nTNumLen = 0;
				nj++;
			}
			vLogic[nj] = pParam->strPattern[ni];
			vGid[nj] = -1.0;
			nj++;	
		}

		else if( (pParam->strPattern[ni] >= '0') && (pParam->strPattern[ni] <= '9') )
		{
			vTNumber[nTNumLen] = pParam->strPattern[ni];
			nTNumLen++;
		}

		else if( (pParam->strPattern[ni] == ' ') || (pParam->strPattern[ni] == '\t') )
		{
			/* do nothing */
		}
		else
		{
			printf("Error: TileMapv2_ProbeSelection_MultiSample, %c not supported in logic expressions!\n", pParam->strPattern[ni]);
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
		printf("Error: TileMapv2_ProbeSelection_MultiSample, lack '(' or ')'!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nLogicLen; ni++)
	{
		if(vGid[ni] > 0.5)
		{
			ng = (int)(vGid[ni])-1;
	
			for(nj=0; nj<pParam->nVargroupNum; nj++)
			{
				nAvailableDF = 0;
				nHit = 0;
				for(nl=0; nl<pParam->vVargroupMap[nj]->nWidth; nl++)
				{
					nClustId = pParam->vVargroupMap[nj]->pMatElement[nl]-1;
					
					if(pParam->pGroupSize->pMatElement[nClustId]>0)
						nAvailableDF += (pParam->pGroupSize->pMatElement[nClustId]-1);

					if(nClustId == ng)
					{
						nHit = 1;
					}
				}

				if(nHit == 1)
				{
					if(nAvailableDF <= 0)
					{
						printf("Error: TileMapv2_ProbeSelection_MultiSample, group %d does not have enough d.f. to estimate variance!\n", (nClustId+1));
						exit(EXIT_FAILURE);
					}
					break;
				}
			}
		}
	}

	/* ---- */
	/* init */
	/* ---- */
	nProbeNum = 0;
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		nProbeNum += pData->vSeqData[ni]->nDataNum;
	}
	nSampleNum = pParam->nSampleNum;
	nGroupNum = pParam->nGroupNum;
	nVarGroupNum = pParam->nVargroupNum;
	
	/* ------------- */
	/* create memory */
	/* ------------- */
	pScore = NULL;
	pScore = CreateDoubleMatrix(1, nProbeNum);
	pGroupSizeCopy = NULL;
	pGroupSizeCopy = CreateIntMatrix(1, nGroupNum);
	pDf = NULL;
	pDf = CreateIntMatrix(1, nVarGroupNum);

	vMeans = NULL;
	vMeans = (struct DOUBLEMATRIX **)calloc(nGroupNum, sizeof(struct DOUBLEMATRIX *));
	if(vMeans == NULL)
	{
		printf("Error: TileMapv2_ProbeSelection_MultiSample, cannot create space for simulation.!\n");
		exit(EXIT_FAILURE);
	}
	vVars = NULL;
	vVars = (struct DOUBLEMATRIX **)calloc(nGroupNum, sizeof(struct DOUBLEMATRIX *));
	if(vVars == NULL)
	{
		printf("Error: TileMapv2_ProbeSelection_MultiSample, cannot create space for simulation.!\n");
		exit(EXIT_FAILURE);
	}
	vSDs = NULL;
	vSDs = (struct DOUBLEMATRIX **)calloc(nVarGroupNum, sizeof(struct DOUBLEMATRIX *));
	if(vSDs == NULL)
	{
		printf("Error: TileMapv2_ProbeSelection_MultiSample, cannot create space for simulation.!\n");
		exit(EXIT_FAILURE);
	}
	vMeanDraws = NULL;
	vMeanDraws = (struct DOUBLEMATRIX **)calloc(nGroupNum, sizeof(struct DOUBLEMATRIX *));
	if(vMeanDraws == NULL)
	{
		printf("Error: TileMapv2_ProbeSelection_MultiSample, cannot create space for simulation.!\n");
		exit(EXIT_FAILURE);
	}
	
	for(ni=0; ni<nGroupNum; ni++)
	{
		vMeans[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vMeans[ni] == NULL)
		{
			printf("Error: TileMapv2_ProbeSelection_MultiSample, cannot create space for simulation.!\n");
			exit(EXIT_FAILURE);
		}

		vVars[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vVars[ni] == NULL)
		{
			printf("Error: TileMapv2_ProbeSelection_MultiSample, cannot create space for simulation.!\n");
			exit(EXIT_FAILURE);
		}

		vMeanDraws[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vMeanDraws[ni] == NULL)
		{
			printf("Error: TileMapv2_ProbeSelection_MultiSample, cannot create space for simulation.!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	for(ni=0; ni<nVarGroupNum; ni++)
	{
		vSDs[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vSDs[ni] == NULL)
		{
			printf("Error: TileMapv2_ProbeSelection_MultiSample, cannot create space for simulation!\n");
			exit(EXIT_FAILURE);
		}
	}

	
	/* ---------------- */
	/* get mean         */
	/* ---------------- */

	/* mean */
	for(nj=0; nj<nSampleNum; nj++)
	{
		nClustId = pParam->pGroupLabel->pMatElement[nj]-1;
		if( (nClustId<0) || (nClustId>=nGroupNum) )
			continue;

		pGroupSizeCopy->pMatElement[nClustId] += 1;
		vSum = vMeans[nClustId]->pMatElement;
		nTotalProbeNum = 0;
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			if(pData->vSeqData[ni]->nDataNum <= 0)
				continue;

			vExp = pData->vSeqData[ni]->vData[1+nj]->pMatElement;
			for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
			{
				vSum[nTotalProbeNum+nk] += vExp[nk];
			}
			nTotalProbeNum += pData->vSeqData[ni]->nDataNum;
		}
		if(nTotalProbeNum != nProbeNum)
		{
			printf("Error: TileMapv2_ProbeSelection_MultiSample, probe number do not match!\n");
			exit(EXIT_FAILURE);
		}
	}

	for(nj=0; nj<nGroupNum; nj++)
	{
		if(pGroupSizeCopy->pMatElement[nj] != pParam->pGroupSize->pMatElement[nj])
		{
			printf("Error: TileMapv2_ProbeSelection_MultiSample, class size not match!\n");
			exit(EXIT_FAILURE);
		}

		if(pGroupSizeCopy->pMatElement[nj] > 0)
		{
            dDenom = (double)(pGroupSizeCopy->pMatElement[nj]);
			vSum = vMeans[nj]->pMatElement;
			for(nk=0; nk<nProbeNum; nk++)
			{
				vSum[nk] /= dDenom;
			}
		}
	}

	/* --------------------------------------------- */
	/* if one has extra d.f., compute the variance.  */
	/* --------------------------------------------- */
	/* variance: sum of squares */
	for(nj=0; nj<nSampleNum; nj++)
	{
		nClustId = pParam->pGroupLabel->pMatElement[nj]-1;
		if( (nClustId<0) || (nClustId>=nGroupNum) )
			continue;
		
		vSum = vVars[nClustId]->pMatElement;
		vAve = vMeans[nClustId]->pMatElement;
		nTotalProbeNum = 0;
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			if(pData->vSeqData[ni]->nDataNum <= 0)
				continue;
			
			vExp = pData->vSeqData[ni]->vData[1+nj]->pMatElement;
			
			for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
			{
				dTemp = (vExp[nk]-vAve[nTotalProbeNum+nk]);		
				vSum[nTotalProbeNum+nk] += dTemp*dTemp;
			}
			nTotalProbeNum += pData->vSeqData[ni]->nDataNum;
		}
		if(nTotalProbeNum != nProbeNum)
		{
			printf("Error: TileMapv2_ProbeSelection_MultiSample, probe number do not match!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* variance: estimates */
	for(nj=0; nj<nVarGroupNum; nj++)
	{
		if(pParam->vVargroupMap[nj]->nWidth != pParam->pVargroupSize->pMatElement[nj])
		{
			printf("Error: TileMapv2_ProbeSelection_MultiSample, variance group number not match!\n");
			exit(EXIT_FAILURE);
		}
	
		for(nl=0; nl<pParam->vVargroupMap[nj]->nWidth; nl++)
		{
			nClustId = pParam->vVargroupMap[nj]->pMatElement[nl]-1;
			if(pParam->pGroupSize->pMatElement[nClustId] > 0)
				pDf->pMatElement[nj] += pParam->pGroupSize->pMatElement[nClustId]-1;
			
			vSum = vSDs[nj]->pMatElement;
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
	for(ni=0; ni<nVarGroupNum; ni++)
	{
		dVarmean = 0.0;
		dVarsst = 0.0;
		dB = 0.0;

		if(pDf->pMatElement[ni] > 0)
		{
			dDenom = (double)(pDf->pMatElement[ni]);
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
				vSum[nj] = sqrt(vSum[nj]) + 1e-16;
			}
		}
	}

	
	/* destroy old variances */
	for(ni=0; ni<nGroupNum; ni++)
	{
		DestroyDoubleMatrix(vVars[ni]);
		vVars[ni] = NULL;
	}

	for(ni=0; ni<nVarGroupNum; ni++)
	{
		for(nj=0; nj<pParam->vVargroupMap[ni]->nWidth; nj++)
		{
			nClustId = (pParam->vVargroupMap[ni]->pMatElement)[nj]-1;
			vVars[nClustId] = vSDs[ni];
		}
	}

	/* ----------- */
	/* simulation  */
	/* ----------- */
	pSigCoef = NULL;
	pSigCoef = CreateDoubleMatrix(1, nGroupNum);
	vExp = pSigCoef->pMatElement;
	for(ni=0; ni<nGroupNum; ni++)
	{
		if(pParam->pGroupSize->pMatElement[ni] > 0)
		{
			dDenom = (double)(pParam->pGroupSize->pMatElement[ni]);
			vExp[ni] = sqrt(1.0/dDenom);
		}
		else
		{
			vExp[ni] = 0.0;
		}
	}

	for(ni=0; ni<nLogicLen; ni++)
	{
		if(vGid[ni] > 0.5)
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

	for(ni=0; ni<pParam->nMonteCarloNum; ni++)
	{
		if(ni%100 == 0)
		{
			printf("iter %d...\n", ni);
		}

		vExp = pSigCoef->pMatElement;
		
		/* probeset by probeset */
		for(nk=0; nk<nGroupNum; nk++)
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

		/* add */
		vSum = pScore->pMatElement;
		vEval = pTrueVec->pMatElement;
		for(nj=0; nj<nProbeNum; nj++)
		{
			vSum[nj] += (double)(vEval[nj]);
		}
		DestroyByteMatrix(pTrueVec);
	}
	DestroyDoubleMatrix(pSigCoef);

	/* normalize */
	dDenom = (double)pParam->nMonteCarloNum;
	dLowLimit = 0.5/dDenom;
	dHighLimit = 1.0-dLowLimit;
	vSum = pScore->pMatElement;
	for(nj=0; nj<nProbeNum; nj++)
	{
		vSum[nj] = vSum[nj]/dDenom;

		if(vSum[nj] < dLowLimit)
		{
			vSum[nj] = dLowLimit;
		}
		else if(vSum[nj] > dHighLimit)
		{
			vSum[nj] = dHighLimit;
		}
		vSum[nj] = log(vSum[nj]/(1.0-vSum[nj]));
	}

	/* -------------- */
	/* release memory */
	/* -------------- */
	for(ni=0; ni<nGroupNum; ni++)
	{
		DestroyDoubleMatrix(vMeans[ni]);
		vMeans[ni] = NULL;
		
		DestroyDoubleMatrix(vMeanDraws[ni]);
		vMeanDraws[ni] = NULL;
	}

	for(ni=0; ni<nVarGroupNum; ni++)
	{
		DestroyDoubleMatrix(vSDs[ni]);
		vSDs[ni] = NULL;		
	}

	free(vMeans);
	free(vVars);
	free(vSDs);
	free(vMeanDraws);

	DestroyIntMatrix(pGroupSizeCopy);
	DestroyIntMatrix(pDf);


	/* destroy unnecessary memory and return */
	nClustId = nSampleNum+1;
	nTotalProbeNum = 0;
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		/* chromosome position, ID = 0 */
		/* raw data, ID = (1:nSampleNum) */
		/* test statistics, ID = nSampleNum+1 */
		pData->vSeqData[ni]->vData[nClustId] = NULL;
		pData->vSeqData[ni]->vData[nClustId] = CreateDoubleMatrix(1, pData->vSeqData[ni]->nDataNum);
		if( (pData->vSeqData[ni]->nDataNum > 0) && (pData->vSeqData[ni]->vData[nClustId] == NULL))
		{
			printf("Error: TileMapv2_ProbeSelection_MultiSample, cannot create space for storing test statistics!\n");
			exit(EXIT_FAILURE);
		}

		if(pData->vSeqData[ni]->nDataNum <= 0)
			continue;

		memcpy(pData->vSeqData[ni]->vData[nClustId]->pMatElement, pScore->pMatElement+nTotalProbeNum, sizeof(double)*pData->vSeqData[ni]->nDataNum);
		nTotalProbeNum += pData->vSeqData[ni]->nDataNum;
	}

	pData->pFieldType->pMatElement[nSampleNum+1] = 1;
	DestroyDoubleMatrix(pScore);

	/* -------------- */
	/* save data      */
	/* -------------- */
	pCol = NULL;
	pCol = CreateIntMatrix(1,pData->nColNum);
	if(pCol == NULL)
	{
		printf("Error: TileMapv2_ProbeSelection_MultiSample, cannot create output column information!\n");
		exit(EXIT_FAILURE);
	}
	
	printf("Exporting probe level test statistics...\n");
	pCol->pMatElement[0] = 1;
	pCol->pMatElement[nSampleNum+1] = 1;
	Affy_SaveBAR_Columns_Fast(strOutFile, pData, pCol);

	DestroyIntMatrix(pCol);

	/* -------------- */
	/* release memory */
	/* -------------- */
	nClustId = nSampleNum+1;
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		/* chromosome position, ID = 0 */
		/* raw data, ID = (1:nSampleNum) */
		/* test statistics, ID = nSampleNum+1 */
		DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
		pData->vSeqData[ni]->vData[nClustId] = NULL;
	}
	pData->pFieldType->pMatElement[nSampleNum+1] = 0;

	/* ------ */
	/* return */
	/* ------ */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ProbeSelection_MultiSample_WithMask_Fast()                   */
/*  TileMapv2 probe level summary: multiple sample comparisons.            */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_ProbeSelection_MultiSample_WithMask_Fast(struct tagBARData *pData,
						struct tagBARData *pMask, struct tagTileMapv2Param *pParam, 
						char strOutFile[], char strFCFile[])
{
	/* ------ */
	/* define */
	/* ------ */
	int nLogicLen;
	int ni,nj,nk,nl,nLen;
	int nLeftNum,nRightNum,nTNumLen;
	char vLogic[MED_LINE_LENGTH];
	char strLine[MED_LINE_LENGTH];
	char vTNumber[LINE_LENGTH];
	double vGid[MED_LINE_LENGTH];
	int ng,nAvailableDF,nHit;
	struct DOUBLEMATRIX *vVid[MED_LINE_LENGTH];
	struct BYTEMATRIX *vLid[MED_LINE_LENGTH];
	int nTrue;
	struct BYTEMATRIX *pTrueVec;
	struct DOUBLEMATRIX *pTrueVal;

	int nProbeNum;
	int nGroupNum;
	int nSampleNum;
	int nVarGroupNum;
	struct INTMATRIX *pGroupSizeCopy;
	struct INTMATRIX *pDf;
	struct DOUBLEMATRIX *pScore;
	struct DOUBLEMATRIX *pFinalMask;

	struct DOUBLEMATRIX **vMeans;
	struct DOUBLEMATRIX **vSampleSizes;
	struct DOUBLEMATRIX **vVars;
	struct DOUBLEMATRIX **vSDs;
	struct DOUBLEMATRIX **vPoolDfs;
	struct DOUBLEMATRIX **vMeanDraws;
	
	int nClustId;
	double *vExp,*vSum,*vAve,*vSigma,*vMu,*vSize,*vMask,*vSize2;
	unsigned char *vEval;
	double dDenom;
	double dTemp;
	int nTotalProbeNum;
	struct INTMATRIX *pCol;
	double dLowLimit,dHighLimit;
	
	/* ------------- */
	/* expression    */
	/* ------------- */
	nLogicLen = 0;
	if(pParam->strPattern[0] == '\0')
	{
		printf("Error: TileMapv2_ProbeSelection_MultiSample_WithMask, no logic expression formula available!\n");
		exit(EXIT_FAILURE);
	}
	nLen = (int)strlen(pParam->strPattern);
	nj = 0;
	nLeftNum = 0;
	nRightNum = 0;
	nTNumLen = 0;

	for(ni=0; ni<nLen; ni++)
	{
		if( pParam->strPattern[ni] == '(' )
		{
			if(nTNumLen != 0)
			{
				printf("Error: TileMapv2_ProbeSelection_MultiSample_WithMask, logic expression wrong!\n");
				exit(EXIT_FAILURE);
			}

			vLogic[nj] = pParam->strPattern[ni];
			vGid[nj] = -1.0;
			nj++;
			nLeftNum++;
		}

		else if(pParam->strPattern[ni] == ')')
		{
			if(nTNumLen > 0)
			{
				vTNumber[nTNumLen] = '\0';
				vLogic[nj] = 'G';
				vGid[nj] = atof(vTNumber);
				nTNumLen = 0;
				nj++;
			}
			
			vLogic[nj] = pParam->strPattern[ni];
			vGid[nj] = -1.0;
			nj++;
			nRightNum++;
		}
		else if( (pParam->strPattern[ni] == '<') || (pParam->strPattern[ni] == '>') 
			|| (pParam->strPattern[ni] == '&') || (pParam->strPattern[ni] == '|') )
		{
			if(nTNumLen > 0)
			{
				vTNumber[nTNumLen] = '\0';
				vLogic[nj] = 'G';
				vGid[nj] = atof(vTNumber);
				nTNumLen = 0;
				nj++;
			}
			vLogic[nj] = pParam->strPattern[ni];
			vGid[nj] = -1.0;
			nj++;	
		}

		else if( (pParam->strPattern[ni] >= '0') && (pParam->strPattern[ni] <= '9') )
		{
			vTNumber[nTNumLen] = pParam->strPattern[ni];
			nTNumLen++;
		}

		else if( (pParam->strPattern[ni] == ' ') || (pParam->strPattern[ni] == '\t') )
		{
			/* do nothing */
		}
		else
		{
			printf("Error: TileMapv2_ProbeSelection_MultiSample_WithMask, %c not supported in logic expressions!\n", pParam->strPattern[ni]);
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
		printf("Error: TileMapv2_ProbeSelection_MultiSample_WithMask, lack '(' or ')'!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nLogicLen; ni++)
	{
		if(vGid[ni] > 0.5)
		{
			ng = (int)(vGid[ni])-1;
	
			for(nj=0; nj<pParam->nVargroupNum; nj++)
			{
				nAvailableDF = 0;
				nHit = 0;
				for(nl=0; nl<pParam->vVargroupMap[nj]->nWidth; nl++)
				{
					nClustId = pParam->vVargroupMap[nj]->pMatElement[nl]-1;
					
					if(pParam->pGroupSize->pMatElement[nClustId]>0)
						nAvailableDF += (pParam->pGroupSize->pMatElement[nClustId]-1);

					if(nClustId == ng)
					{
						nHit = 1;
					}
				}

				if(nHit == 1)
				{
					if(nAvailableDF <= 0)
					{
						printf("Error: TileMapv2_ProbeSelection_MultiSample_WithMask, group %d does not have enough d.f. to estimate variance!\n", (nClustId+1));
						exit(EXIT_FAILURE);
					}
					break;
				}
			}
		}
	}

	/* ---- */
	/* init */
	/* ---- */
	nProbeNum = 0;
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		nProbeNum += pData->vSeqData[ni]->nDataNum;
	}
	nSampleNum = pParam->nSampleNum;
	nGroupNum = pParam->nGroupNum;
	nVarGroupNum = pParam->nVargroupNum;
	
	/* ------------- */
	/* create memory */
	/* ------------- */
	pScore = NULL;
	pScore = CreateDoubleMatrix(1, nProbeNum);
	pFinalMask = NULL;
	pFinalMask = CreateDoubleMatrix(1, nProbeNum);
	pGroupSizeCopy = NULL;
	pGroupSizeCopy = CreateIntMatrix(1, nGroupNum);
	pDf = NULL;
	pDf = CreateIntMatrix(1, nVarGroupNum);

	vMeans = NULL;
	vMeans = (struct DOUBLEMATRIX **)calloc(nGroupNum, sizeof(struct DOUBLEMATRIX *));
	if(vMeans == NULL)
	{
		printf("Error: TileMapv2_ProbeSelection_MultiSample_WithMask, cannot create space for simulation.!\n");
		exit(EXIT_FAILURE);
	}
	vVars = NULL;
	vVars = (struct DOUBLEMATRIX **)calloc(nGroupNum, sizeof(struct DOUBLEMATRIX *));
	if(vVars == NULL)
	{
		printf("Error: TileMapv2_ProbeSelection_MultiSample_WithMask, cannot create space for simulation.!\n");
		exit(EXIT_FAILURE);
	}
	vSampleSizes = NULL;
	vSampleSizes = (struct DOUBLEMATRIX **)calloc(nGroupNum, sizeof(struct DOUBLEMATRIX *));
	if(vSampleSizes == NULL)
	{
		printf("Error: TileMapv2_ProbeSelection_MultiSample_WithMask, cannot create space for simulation.!\n");
		exit(EXIT_FAILURE);
	}
	vSDs = NULL;
	vSDs = (struct DOUBLEMATRIX **)calloc(nVarGroupNum, sizeof(struct DOUBLEMATRIX *));
	if(vSDs == NULL)
	{
		printf("Error: TileMapv2_ProbeSelection_MultiSample_WithMask, cannot create space for simulation.!\n");
		exit(EXIT_FAILURE);
	}
	vPoolDfs = NULL;
	vPoolDfs = (struct DOUBLEMATRIX **)calloc(nVarGroupNum, sizeof(struct DOUBLEMATRIX *));
	if(vPoolDfs == NULL)
	{
		printf("Error: TileMapv2_ProbeSelection_MultiSample_WithMask, cannot create space for simulation.!\n");
		exit(EXIT_FAILURE);
	}
	vMeanDraws = NULL;
	vMeanDraws = (struct DOUBLEMATRIX **)calloc(nGroupNum, sizeof(struct DOUBLEMATRIX *));
	if(vMeanDraws == NULL)
	{
		printf("Error: TileMapv2_ProbeSelection_MultiSample_WithMask, cannot create space for simulation.!\n");
		exit(EXIT_FAILURE);
	}
	
	for(ni=0; ni<nGroupNum; ni++)
	{
		vMeans[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vMeans[ni] == NULL)
		{
			printf("Error: TileMapv2_ProbeSelection_MultiSample_WithMask, cannot create space for simulation.!\n");
			exit(EXIT_FAILURE);
		}

		vVars[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vVars[ni] == NULL)
		{
			printf("Error: TileMapv2_ProbeSelection_MultiSample_WithMask, cannot create space for simulation.!\n");
			exit(EXIT_FAILURE);
		}

		vSampleSizes[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vSampleSizes[ni] == NULL)
		{
			printf("Error: TileMapv2_ProbeSelection_MultiSample_WithMask, cannot create space for simulation.!\n");
			exit(EXIT_FAILURE);
		}

		vMeanDraws[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vMeanDraws[ni] == NULL)
		{
			printf("Error: TileMapv2_ProbeSelection_MultiSample_WithMask, cannot create space for simulation.!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	for(ni=0; ni<nVarGroupNum; ni++)
	{
		vSDs[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vSDs[ni] == NULL)
		{
			printf("Error: TileMapv2_ProbeSelection_MultiSample_WithMask, cannot create space for simulation!\n");
			exit(EXIT_FAILURE);
		}

		vPoolDfs[ni] = CreateDoubleMatrix(nProbeNum,1);
		if(vPoolDfs[ni] == NULL)
		{
			printf("Error: TileMapv2_ProbeSelection_MultiSample_WithMask, cannot create space for simulation!\n");
			exit(EXIT_FAILURE);
		}
	}

	
	/* ---------------- */
	/* get mean         */
	/* ---------------- */

	/* mean */
	for(nj=0; nj<nSampleNum; nj++)
	{
		nClustId = pParam->pGroupLabel->pMatElement[nj]-1;
		if( (nClustId<0) || (nClustId>=nGroupNum) )
			continue;

		pGroupSizeCopy->pMatElement[nClustId] += 1;
		vSum = vMeans[nClustId]->pMatElement;
		vSize = vSampleSizes[nClustId]->pMatElement;
		nTotalProbeNum = 0;
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			if(pData->vSeqData[ni]->nDataNum <= 0)
				continue;
			
			vExp = pData->vSeqData[ni]->vData[1+nj]->pMatElement;
			vMask = pMask->vSeqData[ni]->vData[1+nj]->pMatElement;
			for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
			{
				if(vMask[nk] < 0.5)
				{
					vSum[nk] += vExp[nk];
					vSize[nk] += 1.0;
				}
			}
			nTotalProbeNum += pData->vSeqData[ni]->nDataNum;
			vSum += pData->vSeqData[ni]->nDataNum;
			vSize += pData->vSeqData[ni]->nDataNum;
		}
		if(nTotalProbeNum != nProbeNum)
		{
			printf("Error: TileMapv2_ProbeSelection_MultiSample_WithMask, probe number do not match!\n");
			exit(EXIT_FAILURE);
		}
	}

	for(nj=0; nj<nGroupNum; nj++)
	{
		if(pGroupSizeCopy->pMatElement[nj] != pParam->pGroupSize->pMatElement[nj])
		{
			printf("Error: TileMapv2_ProbeSelection_MultiSample_WithMask, class size not match!\n");
			exit(EXIT_FAILURE);
		}

		if(pGroupSizeCopy->pMatElement[nj] > 0)
		{
            dDenom = (double)(pGroupSizeCopy->pMatElement[nj]);
			vSum = vMeans[nj]->pMatElement;
			vSize = vSampleSizes[nj]->pMatElement;
			for(nk=0; nk<nProbeNum; nk++)
			{
				if(vSize[nk] > 0.5)
					vSum[nk] /= vSize[nk];
			}
		}
	}

	/* --------------------------------------------- */
	/* if one has extra d.f., compute the variance.  */
	/* --------------------------------------------- */
	/* variance: sum of squares */
	for(nj=0; nj<nSampleNum; nj++)
	{
		nClustId = pParam->pGroupLabel->pMatElement[nj]-1;
		if( (nClustId<0) || (nClustId>=nGroupNum) )
			continue;
		
		vSum = vVars[nClustId]->pMatElement;
		vAve = vMeans[nClustId]->pMatElement;
		nTotalProbeNum = 0;
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			if(pData->vSeqData[ni]->nDataNum <= 0)
				continue;
			
			vExp = pData->vSeqData[ni]->vData[1+nj]->pMatElement;
			vMask = pMask->vSeqData[ni]->vData[1+nj]->pMatElement;
			
			for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
			{
				if(vMask[nk] < 0.5)
				{
					dTemp = (vExp[nk]-vAve[nk]);		
					vSum[nk] += dTemp*dTemp;
				}
			}
			nTotalProbeNum += pData->vSeqData[ni]->nDataNum;
			vAve += pData->vSeqData[ni]->nDataNum;
			vSum += pData->vSeqData[ni]->nDataNum;
		}
		if(nTotalProbeNum != nProbeNum)
		{
			printf("Error: TileMapv2_ProbeSelection_MultiSample_WithMask, probe number do not match!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* variance: estimates */
	for(nj=0; nj<nVarGroupNum; nj++)
	{
		if(pParam->vVargroupMap[nj]->nWidth != pParam->pVargroupSize->pMatElement[nj])
		{
			printf("Error: TileMapv2_ProbeSelection_MultiSample_WithMask, variance group number not match!\n");
			exit(EXIT_FAILURE);
		}
	
		for(nl=0; nl<pParam->vVargroupMap[nj]->nWidth; nl++)
		{
			nClustId = pParam->vVargroupMap[nj]->pMatElement[nl]-1;
			if(pParam->pGroupSize->pMatElement[nClustId] > 0)
				pDf->pMatElement[nj] += pParam->pGroupSize->pMatElement[nClustId]-1;
			
			vSum = vSDs[nj]->pMatElement;
			vSize2 = vPoolDfs[nj]->pMatElement;
			vExp = vVars[nClustId]->pMatElement;
			vSize = vSampleSizes[nClustId]->pMatElement;
			for(nk=0; nk<nProbeNum; nk++)
			{
				if(vSize[nk] > 0.5)
				{
					vSum[nk] += vExp[nk];
					vSize2[nk] += vSize[nk]-1.0;
				}
			}
		}
	}

	/* ------------- */
	/* shrinking var */
	/* ------------- */
	for(ni=0; ni<nVarGroupNum; ni++)
	{
		TileMapv2_GetShrinkingVar_UnequalDF(vSDs[ni], vPoolDfs[ni], pDf->pMatElement[ni], nProbeNum);
	}

	/* destroy old variances */
	for(ni=0; ni<nGroupNum; ni++)
	{
		DestroyDoubleMatrix(vVars[ni]);
		vVars[ni] = NULL;
	}

	for(ni=0; ni<nVarGroupNum; ni++)
	{
		for(nj=0; nj<pParam->vVargroupMap[ni]->nWidth; nj++)
		{
			nClustId = (pParam->vVargroupMap[ni]->pMatElement)[nj]-1;
			vVars[nClustId] = vSDs[ni];
		}
	}

	/* ------------- */
	/* prepare mask  */
	/* ------------- */
	for(ni=0; ni<nLogicLen; ni++)
	{
		if(vGid[ni] > 0.5)
		{
			nk = (int)(vGid[ni]-1.0);
			vSize = vSampleSizes[nk]->pMatElement;
			vSigma = vVars[nk]->pMatElement;
			vMask = pFinalMask->pMatElement;
			for(nj=0; nj<nProbeNum; nj++)
			{
				if( (vSize[nj] < 0.5) || (vSigma[nj] < 0.0) )
					vMask[nj] = 1.0;
			}
		}
	}

	/* -------------------- */
	/* prepare coefficient  */
	/* -------------------- */
	for(ni=0; ni<nGroupNum; ni++)
	{
		vExp = vSampleSizes[ni]->pMatElement;

		for(nj=0; nj<nProbeNum; nj++)
		{
			if(vExp[nj] > 0.5)
				vExp[nj] = sqrt(1.0/vExp[nj]);
			else
				vExp[nj] = 0.0;
		}
	}

	/* ----------- */
	/* simulation  */
	/* ----------- */
	for(ni=0; ni<nLogicLen; ni++)
	{
		if(vGid[ni] > 0.5)
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

	for(ni=0; ni<pParam->nMonteCarloNum; ni++)
	{
		if(ni%100 == 0)
		{
			printf("iter %d...\n", ni);
		}
		
		/* probeset by probeset */
		for(nk=0; nk<nGroupNum; nk++)
		{
			vAve = vMeanDraws[nk]->pMatElement;
			vMu = vMeans[nk]->pMatElement;
			vSigma = vVars[nk]->pMatElement;
			vExp = vSampleSizes[nk]->pMatElement;
			dTemp = normrnd(0.0, 1.0);
			for(nj=0; nj<nProbeNum; nj++)
			{
				/* simulate */
				vAve[nj] = vMu[nj]+vExp[nj]*vSigma[nj]*dTemp;
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

		/* add */
		vSum = pScore->pMatElement;
		vEval = pTrueVec->pMatElement;
		for(nj=0; nj<nProbeNum; nj++)
		{
			vSum[nj] += (double)(vEval[nj]);
		}
		DestroyByteMatrix(pTrueVec);
	}

	/* normalize */
	dDenom = (double)pParam->nMonteCarloNum;
	dLowLimit = 0.5/dDenom;
	dHighLimit = 1.0-dLowLimit;
	vSum = pScore->pMatElement;
	vMask = pFinalMask->pMatElement;
	for(nj=0; nj<nProbeNum; nj++)
	{
		if(vMask[nj] < 0.5)
			vSum[nj] = vSum[nj]/dDenom;
		else
			vSum[nj] = 0.0;

		if(vSum[nj] < dLowLimit)
		{
			vSum[nj] = dLowLimit;
		}
		else if(vSum[nj] > dHighLimit)
		{
			vSum[nj] = dHighLimit;
		}
		vSum[nj] = log(vSum[nj]/(1.0-vSum[nj]));
	}

	/* -------------- */
	/* release memory */
	/* -------------- */
	for(ni=0; ni<nGroupNum; ni++)
	{
		DestroyDoubleMatrix(vMeans[ni]);
		vMeans[ni] = NULL;

		DestroyDoubleMatrix(vSampleSizes[ni]);
		vSampleSizes[ni] = NULL;
		
		DestroyDoubleMatrix(vMeanDraws[ni]);
		vMeanDraws[ni] = NULL;
	}

	for(ni=0; ni<nVarGroupNum; ni++)
	{
		DestroyDoubleMatrix(vSDs[ni]);
		vSDs[ni] = NULL;	

		DestroyDoubleMatrix(vPoolDfs[ni]);
		vPoolDfs[ni] = NULL;	
	}

	free(vMeans);
	free(vVars);
	free(vSampleSizes);
	free(vSDs);
	free(vPoolDfs);
	free(vMeanDraws);
	
	DestroyIntMatrix(pGroupSizeCopy);
	DestroyIntMatrix(pDf);

	/* destroy unnecessary memory and return */
	nClustId = nSampleNum+1;
	nTotalProbeNum = 0;
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		/* chromosome position, ID = 0 */
		/* raw data, ID = (1:nSampleNum) */
		/* test statistics, ID = nSampleNum+1 */
		pData->vSeqData[ni]->vData[nClustId] = NULL;
		pData->vSeqData[ni]->vData[nClustId] = CreateDoubleMatrix(1, pData->vSeqData[ni]->nDataNum);
		if( (pData->vSeqData[ni]->nDataNum > 0) && (pData->vSeqData[ni]->vData[nClustId] == NULL))
		{
			printf("Error: TileMapv2_ProbeSelection_MultiSample_WithMask, cannot create space for storing test statistics!\n");
			exit(EXIT_FAILURE);
		}

		if(pData->vSeqData[ni]->nDataNum <= 0)
			continue;

		memcpy(pData->vSeqData[ni]->vData[nClustId]->pMatElement, pScore->pMatElement+nTotalProbeNum, sizeof(double)*pData->vSeqData[ni]->nDataNum);
		nTotalProbeNum += pData->vSeqData[ni]->nDataNum;
	}

	pData->pFieldType->pMatElement[nSampleNum+1] = 1;
	DestroyDoubleMatrix(pScore);

	nClustId = nSampleNum+1;
	nTotalProbeNum = 0;
	for(ni=0; ni<pMask->nSeqNum; ni++)
	{
		/* chromosome position, ID = 0 */
		/* raw data, ID = (1:nSampleNum) */
		/* test statistics, ID = nSampleNum+1 */
		pMask->vSeqData[ni]->vData[nClustId] = NULL;
		pMask->vSeqData[ni]->vData[nClustId] = CreateDoubleMatrix(1, pMask->vSeqData[ni]->nDataNum);
		if( (pMask->vSeqData[ni]->nDataNum > 0) && (pMask->vSeqData[ni]->vData[nClustId] == NULL))
		{
			printf("Error: TileMapv2_ProbeSelection_MultiSample_WithMask, cannot create space for storing test statistics!\n");
			exit(EXIT_FAILURE);
		}

		if(pMask->vSeqData[ni]->nDataNum <= 0)
			continue;

		memcpy(pMask->vSeqData[ni]->vData[nClustId]->pMatElement, pFinalMask->pMatElement+nTotalProbeNum, sizeof(double)*pMask->vSeqData[ni]->nDataNum);
		nTotalProbeNum += pMask->vSeqData[ni]->nDataNum;
	}

	pMask->pFieldType->pMatElement[nSampleNum+1] = 7;
	DestroyDoubleMatrix(pFinalMask);

	/* -------------- */
	/* save data      */
	/* -------------- */
	pCol = NULL;
	pCol = CreateIntMatrix(1,pData->nColNum);
	if(pCol == NULL)
	{
		printf("Error: TileMapv2_ProbeSelection_MultiSample_WithMask, cannot create output column information!\n");
		exit(EXIT_FAILURE);
	}
	
	printf("Exporting probe level test statistics...\n");
	pCol->pMatElement[0] = 1;
	pCol->pMatElement[nSampleNum+1] = 1;
	Affy_SaveBAR_Columns_Fast(strOutFile, pData, pCol);

	DestroyIntMatrix(pCol);

	pCol = NULL;
	pCol = CreateIntMatrix(1,pMask->nColNum);
	if(pCol == NULL)
	{
		printf("Error: TileMapv2_ProbeSelection_TwoSample_WithMask, cannot create output column information!\n");
		exit(EXIT_FAILURE);
	}
	
	printf("Exporting masks...\n");
	pCol->pMatElement[0] = 1;
	pCol->pMatElement[nSampleNum+1] = 1;
	sprintf(strLine, "%s.mask", strOutFile);
	Affy_SaveBAR_Columns_Fast(strLine, pMask, pCol);

	DestroyIntMatrix(pCol);

	/* -------------- */
	/* release memory */
	/* -------------- */
	nClustId = nSampleNum+1;
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		/* chromosome position, ID = 0 */
		/* raw data, ID = (1:nSampleNum) */
		/* test statistics, ID = nSampleNum+1 */
		DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
		pData->vSeqData[ni]->vData[nClustId] = NULL;
	}
	pData->pFieldType->pMatElement[nSampleNum+1] = 0;

	nClustId = nSampleNum+1;
	for(ni=0; ni<pMask->nSeqNum; ni++)
	{
		DestroyDoubleMatrix(pMask->vSeqData[ni]->vData[nClustId]);
		pMask->vSeqData[ni]->vData[nClustId] = NULL;
	}
	pMask->pFieldType->pMatElement[nSampleNum+1] = 0;

	/* ------ */
	/* return */
	/* ------ */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_ProbeSelection_RemoveMean()                                  */
/*  Remove mean from each group, this is a preparation for permutations.   */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_ProbeSelection_RemoveMean(struct tagBARData *pData, 
				struct tagTileMapv2Param *pParam, int nNoiseMask, 
				struct tagBARData *pMask)
{
	/* define */
	int ni,nj,nk;
	int nSampleNum;
	int nGroupNum;
	struct INTMATRIX *pGroupSizeCopy;
	int nClustId;
	double *vSum,*vExp,*vSize,*vMask;
    
	/* ---- */
	/* init */
	/* ---- */
	nSampleNum = pParam->nSampleNum;
	nGroupNum = pParam->nGroupNum;
	pGroupSizeCopy = NULL;
	pGroupSizeCopy = CreateIntMatrix(1, pParam->nGroupNum);

	/* ------------- */
	/* create memory */
	/* ------------- */
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		/* chromosome position, ID = 0 */
		/* raw data, ID = (1:nSampleNum) */
		
		/* within sample mean, ID = nSampleNum+1+(0:nGroupNum-1) */
		for(nj=0; nj<nGroupNum; nj++)
		{
			nClustId = nSampleNum+1+nj;
			pData->vSeqData[ni]->vData[nClustId] = NULL;
			pData->vSeqData[ni]->vData[nClustId] = CreateDoubleMatrix(1, pData->vSeqData[ni]->nDataNum);
			if( (pData->vSeqData[ni]->nDataNum > 0) && (pData->vSeqData[ni]->vData[nClustId] == NULL))
			{
				printf("Error: TileMapv2_ProbeSelection_RemoveMean, cannot create space for storing sample mean!\n");
				exit(EXIT_FAILURE);
			}
		}

		/* sample size, ID = nSampleNum+1+nGroupNum+(0:nGroupNum-1) */
		for(nj=0; nj<nGroupNum; nj++)
		{
			nClustId = nSampleNum+1+nGroupNum+nj;
			pData->vSeqData[ni]->vData[nClustId] = NULL;
			pData->vSeqData[ni]->vData[nClustId] = CreateDoubleMatrix(1, pData->vSeqData[ni]->nDataNum);
			if( (pData->vSeqData[ni]->nDataNum > 0) && (pData->vSeqData[ni]->vData[nClustId] == NULL))
			{
				printf("Error: TileMapv2_ProbeSelection_RemoveMean, cannot create space for storing sample size!\n");
				exit(EXIT_FAILURE);
			}
		}
	}
	
	/* ---------------- */
	/* get mean         */
	/* ---------------- */
	/* mean */
	for(nj=0; nj<nSampleNum; nj++)
	{
		nClustId = pParam->pGroupLabel->pMatElement[nj]-1;
		if( (nClustId<0) || (nClustId>=nGroupNum) )
			continue;

		pGroupSizeCopy->pMatElement[nClustId] += 1;
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			if(pData->vSeqData[ni]->nDataNum <= 0)
				continue;
			
			vSum = pData->vSeqData[ni]->vData[nSampleNum+1+nClustId]->pMatElement;
			vSize = pData->vSeqData[ni]->vData[nSampleNum+1+nGroupNum+nClustId]->pMatElement;
			vExp = pData->vSeqData[ni]->vData[1+nj]->pMatElement;
			if(nNoiseMask == 1)
			{
				vMask = pMask->vSeqData[ni]->vData[1+nj]->pMatElement;
			}
			else
			{
				vMask = NULL;
			}

			for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
			{
				if(vMask == NULL)
				{
					vSum[nk] += vExp[nk];
					vSize[nk] += 1.0;
				}
				else
				{
					if(vMask[nk] < 0.5)
					{
						vSum[nk] += vExp[nk];
						vSize[nk] += 1.0;
					}
				}
			}
		}
	}

	for(nj=0; nj<nGroupNum; nj++)
	{
		if(pGroupSizeCopy->pMatElement[nj] != pParam->pGroupSize->pMatElement[nj])
		{
			printf("Error: TileMapv2_ProbeSelection_RemoveMean, class size not match!\n");
			exit(EXIT_FAILURE);
		}

		if(pGroupSizeCopy->pMatElement[nj] > 0)
		{
            for(ni=0; ni<pData->nSeqNum; ni++)
			{
				if(pData->vSeqData[ni]->nDataNum <= 0)
					continue;
				
				vSum = pData->vSeqData[ni]->vData[nSampleNum+1+nj]->pMatElement;
				vSize = pData->vSeqData[ni]->vData[nSampleNum+1+nGroupNum+nj]->pMatElement;
				for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
				{
					if( vSize[nk] > 0.5)
						vSum[nk] /= vSize[nk];
				}
			}
		}
	}

	/* remove mean */
	for(nj=0; nj<nSampleNum; nj++)
	{
		nClustId = pParam->pGroupLabel->pMatElement[nj]-1;
		if( (nClustId<0) || (nClustId>=nGroupNum) )
			continue;

		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			if(pData->vSeqData[ni]->nDataNum <= 0)
				continue;
			
			vSum = pData->vSeqData[ni]->vData[nSampleNum+1+nClustId]->pMatElement;
			vSize = pData->vSeqData[ni]->vData[nSampleNum+1+nGroupNum+nClustId]->pMatElement;
			vExp = pData->vSeqData[ni]->vData[1+nj]->pMatElement;
			if(nNoiseMask == 1)
			{
				vMask = pMask->vSeqData[ni]->vData[1+nj]->pMatElement;
			}
			else
			{
				vMask = NULL;
			}

			for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
			{
				if(vSize[nk] > 0.5)
					vExp[nk] -= vSum[nk];
				else
					vExp[nk] = 0.0;
			}
		}
	}

	/* ------------- */
	/* destroy memory */
	/* ------------- */
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		/* chromosome position, ID = 0 */
		/* raw data, ID = (1:nSampleNum) */
		
		/* within sample mean, ID = nSampleNum+1+(0:nGroupNum-1) */
		for(nj=0; nj<nGroupNum; nj++)
		{
			nClustId = nSampleNum+1+nj;
			DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
			pData->vSeqData[ni]->vData[nClustId] = NULL;
		}

		/* sample size, ID = nSampleNum+1+nGroupNum+(0:nGroupNum-1) */
		for(nj=0; nj<nGroupNum; nj++)
		{
			nClustId = nSampleNum+1+nGroupNum+nj;
			DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nClustId]);
			pData->vSeqData[ni]->vData[nClustId] = NULL;
		}
	}

	DestroyIntMatrix(pGroupSizeCopy);

	/* ------ */
	/* return */
	/* ------ */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_GetShrinkingVar_UnequalDF()                                  */
/*  Variance shrinking: unequal d.f. case.                                 */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_GetShrinkingVar_UnequalDF(struct DOUBLEMATRIX *pSSR, 
					struct DOUBLEMATRIX *pDf, int nMaxDf, int nProbeNum)
{
	/* ------ */
	/* define */
	/* ------ */
	struct DOUBLEMATRIX *pVarmean = NULL;
	struct DOUBLEMATRIX *pVarsst = NULL;
	struct DOUBLEMATRIX *pB = NULL;
	struct DOUBLEMATRIX *pV = NULL;
	struct DOUBLEMATRIX *pK = NULL;
	double *vS,*vDf;
	int ni,nk,ndf;
	double dVarmean = 0.0;
	double dV2,dK,dN;
	double dV0;
	double dTemp;
	double dB;
	int nTotalProbeNum;

	/* init */
	if(nMaxDf <= 0)
	{
		printf("Warning: TileMapv2_GetShrinkingVar_UnequalDF, maximum d.f. = 0!\n");
		return PROC_SUCCESS;
	}

	pVarmean = CreateDoubleMatrix(1,nMaxDf);
	pVarsst = CreateDoubleMatrix(1,nMaxDf);
	pB = CreateDoubleMatrix(1,nMaxDf);
	pV = CreateDoubleMatrix(1,nMaxDf);
	pK = CreateDoubleMatrix(1,nMaxDf);
	if( (pVarmean == NULL) || (pVarsst == NULL) || (pB == NULL) || (pV == NULL) ||
		 (pK == NULL) )
	{
		printf("Error: TileMapv2_GetShrinkingVar_UnequalDF, cannot create memory for variance computation!\n");
		exit(EXIT_FAILURE);
	}

	/* variance */
	vS = pSSR->pMatElement;
	vDf = pDf->pMatElement;
	for(nk=0; nk<nProbeNum; nk++)
	{
		if(vDf[nk] > 0.5)
		{
			vS[nk] /= vDf[nk];
			ndf = (int)(vDf[nk])-1;
			dVarmean += vS[nk];
			pVarmean->pMatElement[ndf] += vS[nk];
			pK->pMatElement[ndf] += 1.0;
		}
	}
	
	nTotalProbeNum = 0;
	for(ni=0; ni<nMaxDf; ni++)
	{
		if(pK->pMatElement[ni] > 0.5)
		{
			pVarmean->pMatElement[ni] /= pK->pMatElement[ni];
		}
		nTotalProbeNum += (int)(pK->pMatElement[ni]);
	}
	if(nTotalProbeNum == 0)
	{
		printf("Warning: TileMapv2_GetShrinkingVar_UnequalDF, d.f. are all zero!\n");
		return PROC_SUCCESS;
	}
	dVarmean /= (double)nTotalProbeNum;

	/* sst */
	vS = pSSR->pMatElement;
	vDf = pDf->pMatElement;
	for(nk=0; nk<nProbeNum; nk++)
	{
		if(vDf[nk] > 0.5)
		{
			ndf = (int)(vDf[nk])-1;
			dTemp = vS[nk]-pVarmean->pMatElement[ndf];
			pVarsst->pMatElement[ndf] += dTemp*dTemp;
		}
	}

	for(ni=0; ni<nMaxDf; ni++)
	{
		/* shrinkage factor */
		ndf = ni+1;		
		dV2 = 2.0/(double)ndf;
		dK = pK->pMatElement[ni];
		dN = 1.0;
		
		if((dK > 0.0) && (pVarsst->pMatElement[ni] > 0.0))
			dB = dV2*(dK-1)/((dN+dV2)*dK)+dN*(dK-1)*dV2*(pVarmean->pMatElement[ni])*(pVarmean->pMatElement[ni])/((dN+dV2)*dN*(pVarsst->pMatElement[ni]));
		else
			dB = 0.0;

		if(dB < 0.0)
			dB = 0.0;
		else if(dB > 1.0)
			dB = 1.0;

		pB->pMatElement[ni] = dB;
		if(dB > (1.0-ZERO_BOUND))
			pV->pMatElement[ni] = 2.0+ndf*dB/ZERO_BOUND;
		else
			pV->pMatElement[ni] = 2.0+ndf*dB/(1.0-dB);
	}

	dV0 = 0.0;
	dN = 0.0;
	dK = 0.0;
	for(ni=0; ni<nMaxDf; ni++)
	{
		if( dK/(double)nTotalProbeNum < 0.6 )
		{
			dV0 += pV->pMatElement[nMaxDf-1-ni];
			dN += 1.0;
			dK += pK->pMatElement[nMaxDf-1-ni];
		}
	}
	dV0 /= dN;
	if(dV0 < 2.0)
		dV0 = 2.0;

	for(ni=0; ni<nMaxDf; ni++)
	{
		ndf = ni+1;
		pB->pMatElement[ni] = (dV0-2.0)/(ndf+dV0-2.0);
		if(pB->pMatElement[ni] < 0.0)
			pB->pMatElement[ni] = 0.0;
	}

	/* shrink variance */
	vS = pSSR->pMatElement;
	vDf = pDf->pMatElement;
		        
	/* xbar */
	for(nk=0; nk<nProbeNum; nk++)
	{
		if(vDf[nk] > 0.5)
		{
			ndf = (int)(vDf[nk])-1;
			vS[nk] = (1.0-pB->pMatElement[ndf])*vS[nk]+pB->pMatElement[ndf]*dVarmean;
			vS[nk] = sqrt(vS[nk])+1e-16;
		}
		else
		{
			vS[nk] = -1e6;
		}
	}

	/* release memory */
	DestroyDoubleMatrix(pVarmean);
	DestroyDoubleMatrix(pVarsst);
	DestroyDoubleMatrix(pB);
	DestroyDoubleMatrix(pV);
	DestroyDoubleMatrix(pK);

	/* ------ */
	/* return */
	/* ------ */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_GetShrinkingVar_FromSSR()                                    */
/*  Variance shrinking: unequal d.f. case.                                 */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_GetShrinkingVar_FromSSR(struct tagBARData *pData, int nMaxDf,
						int nVarId, int nDfId)
{
	/* ------ */
	/* define */
	/* ------ */
	struct DOUBLEMATRIX *pVarmean = NULL;
	struct DOUBLEMATRIX *pVarsst = NULL;
	struct DOUBLEMATRIX *pB = NULL;
	struct DOUBLEMATRIX *pV = NULL;
	struct DOUBLEMATRIX *pK = NULL;
	int nTotalProbeNum;
	double *vS,*vDf;
	int ni,nk,ndf;
	double dVarmean = 0.0;
	double dV2,dK,dN;
	double dV0;
	double dTemp;
	double dB;

	/* init */
	if(nMaxDf <= 0)
	{
		printf("Warning: TileMapv2_GetShrinkingVar_FromSSR, maximum d.f. = 0!\n");
		return PROC_SUCCESS;
	}

	pVarmean = CreateDoubleMatrix(1,nMaxDf);
	pVarsst = CreateDoubleMatrix(1,nMaxDf);
	pB = CreateDoubleMatrix(1,nMaxDf);
	pV = CreateDoubleMatrix(1,nMaxDf);
	pK = CreateDoubleMatrix(1,nMaxDf);
	if( (pVarmean == NULL) || (pVarsst == NULL) || (pB == NULL) || (pV == NULL) ||
		 (pK == NULL) )
	{
		printf("Error: TileMapv2_GetShrinkingVar_FromSSR, cannot create memory for variance computation!\n");
		exit(EXIT_FAILURE);
	}

	/* variance */
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		if(pData->vSeqData[ni]->nDataNum <= 0)
			continue;
				
		vS = pData->vSeqData[ni]->vData[nVarId]->pMatElement;
		vDf = pData->vSeqData[ni]->vData[nDfId]->pMatElement;
		
        
		/* xbar */
		for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
		{
			if(vDf[nk] > 0.5)
			{
				vS[nk] /= vDf[nk];
				ndf = (int)(vDf[nk])-1;
				dVarmean += vS[nk];
				pVarmean->pMatElement[ndf] += vS[nk];
				pK->pMatElement[ndf] += 1.0;
			}
		}
	}

	nTotalProbeNum = 0;
	for(ni=0; ni<nMaxDf; ni++)
	{
		if(pK->pMatElement[ni] > 0.5)
		{
			pVarmean->pMatElement[ni] /= pK->pMatElement[ni];
		}
		nTotalProbeNum += (int)(pK->pMatElement[ni]);
	}
	if(nTotalProbeNum == 0)
	{
		printf("Warning: TileMapv2_GetShrinkingVar_FromSSR, d.f. are all zero!\n");
		return PROC_SUCCESS;
	}
	dVarmean /= (double)nTotalProbeNum;

	/* sst */
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		if(pData->vSeqData[ni]->nDataNum <= 0)
			continue;
				
		vS = pData->vSeqData[ni]->vData[nVarId]->pMatElement;
		vDf = pData->vSeqData[ni]->vData[nDfId]->pMatElement;
			
		/* xbar */
		for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
		{
			if(vDf[nk] > 0.5)
			{
				ndf = (int)(vDf[nk])-1;
				dTemp = vS[nk]-pVarmean->pMatElement[ndf];
				pVarsst->pMatElement[ndf] += dTemp*dTemp;
			}
		}
	}

	for(ni=0; ni<nMaxDf; ni++)
	{
		/* shrinkage factor */
		ndf = ni+1;		
		dV2 = 2.0/(double)ndf;
		dK = pK->pMatElement[ni];
		dN = 1.0;
		
		if((dK > 0.0) && (pVarsst->pMatElement[ni] > 0.0))
			dB = dV2*(dK-1)/((dN+dV2)*dK)+dN*(dK-1)*dV2*(pVarmean->pMatElement[ni])*(pVarmean->pMatElement[ni])/((dN+dV2)*dN*(pVarsst->pMatElement[ni]));
		else
			dB = 0.0;

		if(dB < 0.0)
			dB = 0.0;
		else if(dB > 1.0)
			dB = 1.0;

		pB->pMatElement[ni] = dB;
		if(dB > (1.0-ZERO_BOUND))
			pV->pMatElement[ni] = 2.0+ndf*dB/ZERO_BOUND;
		else
			pV->pMatElement[ni] = 2.0+ndf*dB/(1.0-dB);
	}

	dV0 = 0.0;
	dN = 0.0;
	dK = 0.0;
	for(ni=0; ni<nMaxDf; ni++)
	{
		if( dK/(double)nTotalProbeNum < 0.6 )
		{
			dV0 += pV->pMatElement[nMaxDf-1-ni];
			dN += 1.0;
			dK += pK->pMatElement[nMaxDf-1-ni];
		}
	}
	dV0 /= dN;
	if(dV0 < 2.0)
		dV0 = 2.0;

	for(ni=0; ni<nMaxDf; ni++)
	{
		ndf = ni+1;
		pB->pMatElement[ni] = (dV0-2.0)/(ndf+dV0-2.0);
		if(pB->pMatElement[ni] < 0.0)
			pB->pMatElement[ni] = 0.0;
	}

	/* shrink variance */
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		if(pData->vSeqData[ni]->nDataNum <= 0)
			continue;
				
		vS = pData->vSeqData[ni]->vData[nVarId]->pMatElement;
		vDf = pData->vSeqData[ni]->vData[nDfId]->pMatElement;
		        
		/* xbar */
		for(nk=0; nk<pData->vSeqData[ni]->nDataNum; nk++)
		{
			if(vDf[nk] > 0.5)
			{
				ndf = (int)(vDf[nk])-1;
				vS[nk] = (1.0-pB->pMatElement[ndf])*vS[nk]+pB->pMatElement[ndf]*dVarmean+1e-16;
			}
			else
			{
				vS[nk] = -1e6;
			}
		}
	}

	/* release memory */
	DestroyDoubleMatrix(pVarmean);
	DestroyDoubleMatrix(pVarsst);
	DestroyDoubleMatrix(pB);
	DestroyDoubleMatrix(pV);
	DestroyDoubleMatrix(pK);

	/* ------ */
	/* return */
	/* ------ */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_RegionDetection_MA_Main()                                    */
/*  TileMapv2 region detection, MA.                                        */
/*  return number of regions detected.                                     */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_RegionDetection_MA_Main(struct tagTileMapv2Param *pParam, int nLibId)
{
	/* define */
	int nRegionNum = 0;
	char strProbeFile[MED_LINE_LENGTH];
	char strMaskFile[MED_LINE_LENGTH];
	char strFCFile[MED_LINE_LENGTH];
	char strFileName[MED_LINE_LENGTH];
	struct INTMATRIX *pCol = NULL;
	struct INTMATRIX *pPermGroupLabel = NULL;
	struct INTMATRIX *pOriGroupLabel = NULL;
	struct tagBARData *pData = NULL;
	struct DOUBLEMATRIX *pRegion = NULL;
	struct DOUBLEMATRIX *pNewRegion = NULL;
	struct DOUBLEMATRIX *pRegionMaxSort = NULL;
	struct LONGMATRIX *pRegionMaxSid = NULL;
	struct DOUBLEMATRIX *pRegionMaxFDR = NULL;
	struct DOUBLEMATRIX *pRegionSumSort = NULL;
	struct LONGMATRIX *pRegionSumSid = NULL;
	struct DOUBLEMATRIX *pRegionSumFDR = NULL;
	struct INTMATRIX *pType = NULL;
	struct INTMATRIX *pPriority = NULL;

	struct DOUBLEMATRIX *pRegionCT = NULL;
	struct DOUBLEMATRIX *pRegionCTMaxSort = NULL;
	struct LONGMATRIX *pRegionCTMaxSid = NULL;
	struct DOUBLEMATRIX *pRegionCTSumSort = NULL;
	struct LONGMATRIX *pRegionCTSumSid = NULL;
	FILE *fpOut;
	char strTemp[MED_LINE_LENGTH];

	/* init */
	if(pParam == NULL)
	{
		return 0;
	}

	/* ------------- */
	/* load raw data */
	/* ------------- */
	strcpy(strMaskFile, "");
	strcpy(strFCFile, "");
	sprintf(strProbeFile, "%s%s_%d.pb.bar", pParam->strWorkPath, pParam->strProjectTitle, nLibId+1);
	if(pParam->nNoiseMask == 1)
	{
		sprintf(strMaskFile, "%s%s_%d.pb.bar.mask", pParam->strWorkPath, pParam->strProjectTitle, nLibId+1);
	}
	if(pParam->nPatternType <= 2)
	{
		sprintf(strFCFile, "%s%s_%d.fc.bar", pParam->strWorkPath, pParam->strProjectTitle, nLibId+1);
	}
	
	pData = NULL;
	pData = TileMapv2_RegionDetection_MA_PrepareData(strProbeFile, strMaskFile, strFCFile);
	if(pData == NULL)
	{
		printf("Error: TileMapv2_RegionDetection_MA_Main, cannot load raw data!\n");
		exit(EXIT_FAILURE);
	}

	/* -------------------- */
	/* MA                   */
	/* -------------------- */
	/* MA */
	TileMapv2_MA(pData, pParam);

	/* Save MA statistics */
	pCol = NULL;
	pCol = CreateIntMatrix(1,pData->nColNum);
	if(pCol == NULL)
	{
		printf("Error: TileMapv2_RegionDetection_MA_Main, cannot create output column information!\n");
		exit(EXIT_FAILURE);
	}
	printf("Exporting MA statistics...\n");
	pCol->pMatElement[0] = 1;
	pCol->pMatElement[2] = 1;
	sprintf(strFileName, "%s%s_%d.ma.bar", pParam->strWorkPath, pParam->strProjectTitle, nLibId+1);
	Affy_SaveBAR_Columns_Fast(strFileName, pData, pCol);

	pCol->pMatElement[2] = 0;
	pCol->pMatElement[3] = 1;
	sprintf(strFileName, "%s%s_%d.maN.bar", pParam->strWorkPath, pParam->strProjectTitle, nLibId+1);
	Affy_SaveBAR_Columns_Fast(strFileName, pData, pCol);
	DestroyIntMatrix(pCol);

	/* Call binding regions */
	pRegion = NULL;
	pRegionMaxFDR = NULL;
	pRegionSumFDR = NULL;
	pRegion = TileMapv2_MA_CallRegion(pData, pParam, nLibId);
	
	/* ---------------------------- */
	/* permutations to estimate FDR */
	/* ---------------------------- */
	if(pRegion != NULL)
	{
		pRegionMaxFDR = CreateDoubleMatrix(1, pRegion->nHeight);
		pRegionSumFDR = CreateDoubleMatrix(1, pRegion->nHeight);
		if( (pRegionMaxFDR == NULL) || (pRegionSumFDR == NULL) )
		{
			printf("Error: TileMapv2_RegionDetection_MA_Main, cannot allocate memory for FDR computation!\n");
			exit(EXIT_FAILURE);
		}

		pType = NULL;
		pType = CreateIntMatrix(1, pRegion->nWidth);
		if(pType == NULL)
		{
			printf("Error: TileMapv2_RegionDetection_MA_Main, cannot create column data type!\n");
			exit(EXIT_FAILURE);
		}
		
		/* "%d\t%d\t%d\t% 9.7e\t%d\t% 9.7e\t%d\t% 9.7e\t%d\n", 
						nSeqId, nStart, nEnd, dMaxT, nMaxTPos,
						dMaxFC, nMaxFCPos, dSumT, nProbeCount); */

		pType->pMatElement[0] = 2;
		pType->pMatElement[1] = 2;
		pType->pMatElement[2] = 2;
		pType->pMatElement[3] = 1;
		pType->pMatElement[4] = 2;
		pType->pMatElement[5] = 1;
		pType->pMatElement[6] = 2;
		pType->pMatElement[7] = 1;
		pType->pMatElement[8] = 2;

		pPriority = NULL;
		pPriority = CreateIntMatrix(1, 1);
		if(pPriority == NULL)
		{
			printf("Error: TileMapv2_RegionDetection_MA_Main, cannot create sorting key!\n");
			exit(EXIT_FAILURE);
		}

		pPriority->pMatElement[0] = 3;
		pRegionMaxSort = NULL;
		pRegionMaxSid = NULL;
		DMSORTROWS(pRegion, pType, pPriority, &pRegionMaxSort, &pRegionMaxSid);
		
		pPriority->pMatElement[0] = 7;
		pRegionSumSort = NULL;
		pRegionSumSid = NULL;
		DMSORTROWS(pRegion, pType, pPriority, &pRegionSumSort, &pRegionSumSid);
		
		DestroyIntMatrix(pPriority);
	
		if( (pParam->nFDRType == 1) && (pParam->nPermutationNum > 0) )
		{
			/* permutations */
			TileMapv2_MA_RegionFDR_Permutation(pType, 
				pRegionMaxSort, pRegionMaxFDR,
				pRegionSumSort, pRegionSumFDR, 
				pParam, nLibId);

			pNewRegion = NULL;
			pNewRegion = TileMapv2_RegionFDR_Reorganize(pRegion, pRegionMaxFDR, pRegionMaxSid, 
				pRegionSumFDR, pRegionSumSid);
		}
		else if(pParam->nFDRType == 0)
		{
			/* convert test statistics */
			TileMapv2_MA_RegionFDR_LeftTail(pType,
				pRegionMaxSort, pRegionMaxFDR,
				pRegionSumSort, pRegionSumFDR, 
				pData, pParam, nLibId);

			pNewRegion = NULL;
			pNewRegion = TileMapv2_RegionFDR_Reorganize(pRegion, pRegionMaxFDR, pRegionMaxSid, 
				pRegionSumFDR, pRegionSumSid); 
		}
		else if (pParam->nFDRType == 2)
		{
			pNewRegion = NULL;
			pNewRegion = TileMapv2_MA_RegionFDR_UMS(pRegion, pData, pParam, nLibId);
		}
		else
		{
			pNewRegion = NULL;
			pNewRegion = TileMapv2_RegionFDR_Reorganize(pRegion, pRegionMaxFDR, pRegionMaxSid, 
				pRegionSumFDR, pRegionSumSid);
		}

		/* ---------------------------- */
		/* exporting results            */
		/* ---------------------------- */
		/* seqid, chr, start, end, strand, length, maxMA, maxMA Pos, maxMA FDR, maxFC, maxFC Pos,
		sumMA, sumMA good probe count, sumMA FDR, libid */
		nRegionNum = pNewRegion->nHeight;
		TileMapv2_MA_ExportRegion(pNewRegion, pData, pParam, nLibId);
	
		/* -------------- */
		/* release memory */
		/* -------------- */
		DestroyDoubleMatrix(pRegion);
		DestroyIntMatrix(pType);
		DestroyDoubleMatrix(pRegionMaxFDR);
		DestroyDoubleMatrix(pRegionSumFDR);
		DestroyDoubleMatrix(pRegionMaxSort);
		DestroyDoubleMatrix(pRegionSumSort);
		DestroyLongMatrix(pRegionMaxSid);
		DestroyLongMatrix(pRegionSumSid);
		DestroyDoubleMatrix(pNewRegion);
	}
	else
	{
		sprintf(strTemp, "%s%s_%d.cod", pParam->strWorkPath, pParam->strProjectTitle, (nLibId+1));
		fpOut = NULL;
		fpOut = fopen(strTemp, "w");
		if(fpOut == NULL)
		{
			printf("Error: TileMapv2_RegionDetection_MA_Main, cannot open the output file!\n");
			exit(EXIT_FAILURE);
		}

		fclose(fpOut);
	}
	

	/* -------------- */
	/* release memory */
	/* -------------- */
	Affy_BARData_Destroy(&pData);
	
	/* return */
	return nRegionNum;
}


/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_RegionDetection_MA_PrepareData()                             */
/*  TileMapv2 MA preparation: loading probe level statistics.              */
/* ----------------------------------------------------------------------- */ 
struct tagBARData *TileMapv2_RegionDetection_MA_PrepareData(char strProbeFile[],
							char strMaskFile[], char strFCFile[])
{
	struct tagBARData *pData = NULL;
	struct tagBARData *pMask = NULL;
	struct INTMATRIX *pFieldType = NULL;
	struct DOUBLEMATRIX **vDataVec = NULL;
	int nTotalProbeNum = 0;
	int nj;
	
	/* -------------------- */
	/* load probe level t   */
	/* -------------------- */
    pData = NULL;
	pData = Affy_LoadBAR_Fast(strProbeFile);
	if(pData == NULL)
	{
		printf("Error: TileMapv2_RegionDetection_MA_PrepareData, cannot load raw data!\n");
		exit(EXIT_FAILURE);
	}

	/* 1 position + 1 probe level statistics + 1 MA + 1 N + 1 mask + 1 Fold change */ 
	pData->nColNum = 2+4;
	pFieldType = pData->pFieldType;
	pData->pFieldType = NULL;
	pData->pFieldType = CreateIntMatrix(1,pData->nColNum);
	if(pData->pFieldType == NULL)
	{
		printf("Error: TileMapv2_RegionDetection_MA_PrepareData, cannot create memory for field types!\n");
		exit(EXIT_FAILURE);
	}
	pData->pFieldType->pMatElement[0] = pFieldType->pMatElement[0];
	pData->pFieldType->pMatElement[1] = pFieldType->pMatElement[1];
	DestroyIntMatrix(pFieldType);
	pData->pFieldType->pMatElement[2] = 1;
	pData->pFieldType->pMatElement[3] = 2;
	
	for(nj=0; nj<pData->nSeqNum; nj++)
	{
		nTotalProbeNum += pData->vSeqData[nj]->nDataNum;
		pData->vSeqData[nj]->nColNum = pData->nColNum;
		vDataVec = pData->vSeqData[nj]->vData;
		pData->vSeqData[nj]->vData = NULL;
		pData->vSeqData[nj]->vData = (struct DOUBLEMATRIX **)calloc(pData->vSeqData[nj]->nColNum, sizeof(struct DOUBLEMATRIX *));
		if(pData->vSeqData[nj]->vData == NULL)
		{
			printf("Error: TileMapv2_RegionDetection_MA_PrepareData, cannot create memory for tracking intensity data!\n");
			exit(EXIT_FAILURE);
		}
		pData->vSeqData[nj]->vData[0] = vDataVec[0];
		pData->vSeqData[nj]->vData[1] = vDataVec[1];
		free(vDataVec);

		if(pData->vSeqData[nj]->nDataNum <= 0)
			continue;

		pData->vSeqData[nj]->vData[2] = CreateDoubleMatrix(1, pData->vSeqData[nj]->nDataNum);
		if(pData->vSeqData[nj]->vData[2] == NULL)
		{
			printf("Error: TileMapv2_RegionDetection_MA_PrepareData, cannot create memory for storing MA statistics!\n");
			exit(EXIT_FAILURE);
		}

		pData->vSeqData[nj]->vData[3] = CreateDoubleMatrix(1, pData->vSeqData[nj]->nDataNum);
		if(pData->vSeqData[nj]->vData[3] == NULL)
		{
			printf("Error: TileMapv2_RegionDetection_MA_PrepareData, cannot create memory for storing MA number!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* -------------------- */
	/* load masks if needed */
	/* -------------------- */
	if(strcmp(strMaskFile, "") != 0)
	{
		pMask = NULL;
		pMask = Affy_LoadBAR_Fast(strMaskFile);
		if(pMask == NULL)
		{
			printf("Error: TileMapv2_RegionDetection_MA_PrepareData, cannot load mask data!\n");
			exit(EXIT_FAILURE);
		}

		if(pMask->nSeqNum != pData->nSeqNum)
		{
			printf("Error: TileMapv2_RegionDetection_MA_PrepareData, array types do not match!\n");
			exit(EXIT_FAILURE);
		}

		for(nj=0; nj<pMask->nSeqNum; nj++)
		{
			if(pMask->vSeqData[nj]->nDataNum != pData->vSeqData[nj]->nDataNum)
			{
				printf("Error: TileMapv2_RegionDetection_MA_PrepareData, array types do not match!\n");
				exit(EXIT_FAILURE);
			}
			pData->vSeqData[nj]->vData[4] = pMask->vSeqData[nj]->vData[1];
			pMask->vSeqData[nj]->vData[1] = NULL;
		}
        
		pData->pFieldType->pMatElement[4] = pMask->pFieldType->pMatElement[1];

		Affy_BARData_Destroy(&pMask);
	}


	/* -------------------- */
	/* load FC if needed    */
	/* -------------------- */
	if(strcmp(strFCFile, "") != 0)
	{
		pMask = NULL;
		pMask = Affy_LoadBAR_Fast(strFCFile);
		if(pMask == NULL)
		{
			/* return */
			return pData;
		}

		if(pMask->nSeqNum != pData->nSeqNum)
		{
			printf("Error: TileMapv2_RegionDetection_MA_PrepareData, array types do not match!\n");
			exit(EXIT_FAILURE);
		}

		for(nj=0; nj<pMask->nSeqNum; nj++)
		{
			if(pMask->vSeqData[nj]->nDataNum != pData->vSeqData[nj]->nDataNum)
			{
				printf("Error: TileMapv2_RegionDetection_MA_PrepareData, array types do not match!\n");
				exit(EXIT_FAILURE);
			}
			pData->vSeqData[nj]->vData[5] = pMask->vSeqData[nj]->vData[1];
			pMask->vSeqData[nj]->vData[1] = NULL;
		}
        
		pData->pFieldType->pMatElement[5] = pMask->pFieldType->pMatElement[1];

		Affy_BARData_Destroy(&pMask);
	}

	/* return */
	return pData;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_MA()                                                         */
/*  TileMapv2 MA.                                                          */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_MA(struct tagBARData *pData, struct tagTileMapv2Param *pParam)
{
	/* define */
	int ni,nj,nk;
	double *vT,*vM,*vPos,*vMask,*vN,*vF;
	int nPos,nA,nP;
	double dSum;
	int nN;
	int nContribute;
	int nTotalProbeNum;
	double dMean,dSD;

	/* init */
	if( (pParam->nW < 0) || (pParam->nWindowBoundary < 0) ) 
	{
		printf("Error: TileMapv2_MA, W must be >=0!\n");
		exit(EXIT_FAILURE);
	}
	
	/* if nW == 0 */
	if( (pParam->nW == 0) || (pParam->nWindowBoundary == 0) )
	{
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			if(pData->vSeqData[ni]->nDataNum <= 0)
				continue;
		
			vT = pData->vSeqData[ni]->vData[1]->pMatElement;
			vM = pData->vSeqData[ni]->vData[2]->pMatElement;
			vN = pData->vSeqData[ni]->vData[3]->pMatElement;
			for(nj=0; nj<pData->vSeqData[ni]->nDataNum; nj++)
			{
				vM[nj] = vT[nj];
				vN[nj] = 1.0;
			}
		}
	}

	/* if nW > 0 */
	else
	{

		/* FOR DEBUG */
		/* for(ni=0; ni<pData->nSeqNum; ni++)
		{
			vPos = pData->vSeqData[ni]->vData[0]->pMatElement;
			vT = pData->vSeqData[ni]->vData[1]->pMatElement;
			vM = pData->vSeqData[ni]->vData[5]->pMatElement;
			vN = pData->vSeqData[ni]->vData[3]->pMatElement;
			if(pData->vSeqData[ni]->vData[4] != NULL)
			{
				vMask = pData->vSeqData[ni]->vData[4]->pMatElement;
			}
			else
			{
				vMask = NULL;
			}

			for(nPos=0; nPos<pData->vSeqData[ni]->nDataNum; nPos++)
			{
				dSum = 0.0;
				nN = 0;
				nA = 0;
				nP = 0;
				if(vMask == NULL)
				{
					dSum = vT[nPos];
					nN = 1;
				}
				else
				{
					if(vMask[nPos] < 0.5)
					{
						dSum = vT[nPos];
						nN = 1;
					}
				}

				nk = nPos;
				while(nA<pParam->nW)
				{
					nk--;
					if(nk<0)
						break;
					if( (vPos[nPos]-vPos[nk]) > pParam->nWindowBoundary)
						break;

					if(vMask == NULL)
					{
						dSum += vT[nk];
						nA++;
						nN++;
					}
					else
					{
						if(vMask[nk] < 0.5)
						{
							dSum += vT[nk];
							nA++;
							nN++;
						}
					}
				}

				nk = nPos;
				while(nP<pParam->nW)
				{
					nk++;
					if(nk>=pData->vSeqData[ni]->nDataNum)
						break;
					if( (vPos[nk]-vPos[nPos]) > pParam->nWindowBoundary)
						break;

					if(vMask == NULL)
					{
						dSum += vT[nk];
						nP++;
						nN++;
					}
					else
					{
						if(vMask[nk] < 0.5)
						{
							dSum += vT[nk];
							nP++;
							nN++;
						}
					}
				}

				if(vMask == NULL)
				{
					if(nN != (nA+nP+1))
					{
						printf("Error!\n");
						exit(EXIT_FAILURE);
					}
				}
				else
				{
					if(vMask[nPos] < 0.5)
					{
						if(nN != (nA+nP+1))
						{
							printf("Error!\n");
							exit(EXIT_FAILURE);
						}
					}
					else
					{
						if(nN != (nA+nP))
						{
							printf("Error!\n");
							exit(EXIT_FAILURE);
						}
					}
				}

				if(nN > 0)
				{
					vM[nPos] = dSum/(double)nN;
					vN[nPos] = nN;
				}
				else
				{
					vM[nPos] = 0.0;
					vN[nPos] = 0;
				}
			}
		} */

		/* if nW > 0 */
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			if(pData->vSeqData[ni]->nDataNum <= 0)
				continue;
			
			vPos = pData->vSeqData[ni]->vData[0]->pMatElement;
			vT = pData->vSeqData[ni]->vData[1]->pMatElement;
			vM = pData->vSeqData[ni]->vData[2]->pMatElement;
			vN = pData->vSeqData[ni]->vData[3]->pMatElement;
			if(pData->vSeqData[ni]->vData[4] != NULL)
			{
				vMask = pData->vSeqData[ni]->vData[4]->pMatElement;
			}
			else
			{
				vMask = NULL;
			}
			if(pData->vSeqData[ni]->vData[5] != NULL)
			{
				vF = pData->vSeqData[ni]->vData[5]->pMatElement;
			}
			else
			{
				vF = NULL;
			}

			/* the first position */
			nPos = 0;
			nA = 0;
			nP = 0;
			nN = 0;
			dSum = 0.0;
			if(vMask == NULL)
			{
				dSum = vT[nPos];
				nN = 1;
			}
			else
			{
				if(vMask[nPos] < 0.5)
				{
					dSum = vT[nPos];
					nN = 1;
				}
			}
			
			nj = 1;
			nk = 0;
			while(nP<pParam->nW) 
			{
				if( nj >= pData->vSeqData[ni]->nDataNum )
					break;
				if( (vPos[nj]-vPos[nPos]) > pParam->nWindowBoundary )
					break;

				if(vMask == NULL)
				{
					dSum += vT[nj];
					nP++;
					nN++;
				}
				else
				{
					if(vMask[nj] < 0.5)
					{
						dSum += vT[nj];
						nP++;
						nN++;
					}
				}

				nj++;
			}

			if(vMask == NULL)
			{
				if(nN != (nA+nP+1))
				{
					printf("Error: TileMapv2_MA, MA sum was wrong!\n");
					exit(EXIT_FAILURE);
				}
				vM[nPos] = dSum/(double)nN;
			}
			else
			{
				if(vMask[nPos] < 0.5)
				{
					if(nN != (nA+nP+1))
					{
						printf("Error: TileMapv2_MA, MA sum was wrong!\n");
						exit(EXIT_FAILURE);
					}
					vM[nPos] = dSum/(double)nN;
				}
				else
				{
					if(nN != (nA+nP))
					{
						printf("Error: TileMapv2_MA, MA sum was wrong!\n");
						exit(EXIT_FAILURE);
					}
					if(nN > 0)
					{
						vM[nPos] = dSum/(double)nN;
					}
					else
					{
						nN = 0;
						vM[nPos] = 0.0;
					}
				}
			}

			/* FOR DEBUG */
			/* 
			if(nN!=(int)(vN[nPos]))
			{
				printf("Error!\n");
				exit(EXIT_FAILURE);
			}
			if(fabs(vM[nPos]-vF[nPos]) > ZERO_BOUND )
			{
				printf("Error!\n");
				exit(EXIT_FAILURE);
			} */

			vN[nPos] = nN;

			/* following positions */
			for(nPos=1; nPos<pData->vSeqData[ni]->nDataNum; nPos++)
			{
				nContribute = 0;
				if(vMask == NULL)
				{
					nContribute = 1;
				}
				else
				{
					if(vMask[nPos-1] < 0.5)
						nContribute = 1;
				}

				if(nContribute == 1)
					nA++;

				nContribute = 0;
				if(nPos < nj)
				{
					if(vMask == NULL)
					{
						nContribute = 1;
						nP--;
					}
					else
					{
						if(vMask[nPos] < 0.5)
						{
							nContribute = 1;
							nP--;
						}
						else
						{
							nContribute = 0;
						}
					}
				}
				else
				{
					nContribute = 0;
					if(vMask == NULL)
					{
						dSum += vT[nPos];
						nN++;
					}
					else
					{
						if(vMask[nPos] < 0.5)
						{
							dSum += vT[nPos];
							nN++;
						}
					}

					if(nj<=nPos)
						nj = nPos+1;
				}

				/* remove previous */
				while( (nA>pParam->nW) || ((vPos[nPos]-vPos[nk]) > pParam->nWindowBoundary) )
				{
					if(vMask == NULL)
					{
						dSum -= vT[nk];
						nA--;
						nN--;
					}
					else
					{
						if(vMask[nk] < 0.5)
						{
							dSum -= vT[nk];
							nA--;
							nN--;
						}
					}

					nk++;
					if(nk == nPos)
						break;
				}

				/* add after */
				while(nP<pParam->nW)
				{
					if(nj >= pData->vSeqData[ni]->nDataNum)
						break;
					if( (vPos[nj]-vPos[nPos]) > pParam->nWindowBoundary)
						break;

					if(vMask == NULL)
					{
						dSum += vT[nj];
						nP++;
						nN++;
					}
					else
					{
						if(vMask[nj] < 0.5)
						{
							dSum += vT[nj];
							nP++;
							nN++;
						}
					}

					nj++;
				}

				/* average */
				if(vMask == NULL)
				{
					if(nN != (nA+nP+1))
					{
						printf("Error: TileMapv2_MA, MA sum was wrong!\n");
						exit(EXIT_FAILURE);
					}
					vM[nPos] = dSum/(double)nN;
				}
				else
				{
					if(vMask[nPos] < 0.5)
					{
						if(nN != (nA+nP+1))
						{
							printf("Error: TileMapv2_MA, MA sum was wrong!\n");
							exit(EXIT_FAILURE);
						}
						vM[nPos] = dSum/(double)nN;
					}
					else
					{
						if(nN != (nA+nP))
						{
							printf("Error: TileMapv2_MA, MA sum was wrong!\n");
							exit(EXIT_FAILURE);
						}
						if(nN > 0)
						{
							vM[nPos] = dSum/(double)nN;
						}
						else
						{
							nN = 0;
							vM[nPos] = 0.0;
						}
					}
				}

				/* FOR DEBUG */
				/*
				if(nN!=(int)(vN[nPos]))
				{
					printf("Error!\n");
					exit(EXIT_FAILURE);
				}
				if(fabs(vM[nPos]-vF[nPos]) > ZERO_BOUND)
				{
					printf("Error!\n");
					exit(EXIT_FAILURE);
				} */

				vN[nPos] = nN;
			}
		}
	}

	/* compute mean */
	dMean = 0.0;
	nTotalProbeNum = 0;
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		if(pData->vSeqData[ni]->nDataNum <= 0)
			continue;

		vM = pData->vSeqData[ni]->vData[2]->pMatElement;
		vN = pData->vSeqData[ni]->vData[3]->pMatElement;
		if(pData->vSeqData[ni]->vData[4] != NULL)
		{
			vMask = pData->vSeqData[ni]->vData[4]->pMatElement;
		}
		else
		{
			vMask = NULL;
		}

		if(vMask == NULL)
		{
			for(nj=0; nj<pData->vSeqData[ni]->nDataNum; nj++)
			{
				if(vN[nj] > pParam->nW)
				{
					dMean += vM[nj];
					nTotalProbeNum += 1;
				}
			}
		}
		else
		{
			for(nj=0; nj<pData->vSeqData[ni]->nDataNum; nj++)
			{
				if( (vMask[nj] < 0.5) && (vN[nj] > pParam->nW) )
				{
					dMean += vM[nj];
					nTotalProbeNum += 1;
				}
			}
		}
	}

	printf("Total Number of Good Probes = %d\n", nTotalProbeNum);

	if(nTotalProbeNum > 0)
	{
		dMean /= (double)nTotalProbeNum;
	}
	else
	{
		printf("Warning: TileMapv2_MA, no usable probes!\n");
	}

	/* compute variance */
	dSD = 0.0;
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		if(pData->vSeqData[ni]->nDataNum <= 0)
			continue;

		vM = pData->vSeqData[ni]->vData[2]->pMatElement;
		vN = pData->vSeqData[ni]->vData[3]->pMatElement;
		if(pData->vSeqData[ni]->vData[4] != NULL)
		{
			vMask = pData->vSeqData[ni]->vData[4]->pMatElement;
		}
		else
		{
			vMask = NULL;
		}

		if(vMask == NULL)
		{
			for(nj=0; nj<pData->vSeqData[ni]->nDataNum; nj++)
			{
				if(vN[nj] > pParam->nW)
				{
					dSD += (vM[nj]-dMean)*(vM[nj]-dMean);
				}
			}
		}
		else
		{
			for(nj=0; nj<pData->vSeqData[ni]->nDataNum; nj++)
			{
				if( (vMask[nj] < 0.5) && (vN[nj] > pParam->nW) )
				{
					dSD += (vM[nj]-dMean)*(vM[nj]-dMean);
				}
			}
		}
	}

	if(nTotalProbeNum > 1)
	{
		dSD /= (double)(nTotalProbeNum-1);
		dSD = sqrt(dSD);
	}
	else
	{
		printf("Warning: TileMapv2_MA, no usable probes!\n");
	}

	/* standardize */
	if(pParam->nMAStandardize == 1)
	{
		/* compute std */
		if( (nTotalProbeNum > 1) && (dSD > 0.0) )
		{
			for(ni=0; ni<pData->nSeqNum; ni++)
			{
				if(pData->vSeqData[ni]->nDataNum <= 0)
					continue;
				
				vM = pData->vSeqData[ni]->vData[2]->pMatElement;
				vN = pData->vSeqData[ni]->vData[3]->pMatElement;
				if(pData->vSeqData[ni]->vData[4] != NULL)
				{
					vMask = pData->vSeqData[ni]->vData[4]->pMatElement;
				}
				else
				{
					vMask = NULL;
				}

				if(vMask == NULL)
				{
					for(nj=0; nj<pData->vSeqData[ni]->nDataNum; nj++)
					{
						if(vN[nj] > pParam->nW)
							vM[nj] = (vM[nj]-dMean)/dSD;
						else
							vM[nj] = 0.0;
					}
				}
				else
				{
					for(nj=0; nj<pData->vSeqData[ni]->nDataNum; nj++)
					{
						if( (vMask[nj]<0.5) && (vN[nj] > pParam->nW) )
							vM[nj] = (vM[nj]-dMean)/dSD;
						else
							vM[nj] = 0.0;
					}
				}
			}
		}
		else
		{
			for(ni=0; ni<pData->nSeqNum; ni++)
			{
				if(pData->vSeqData[ni]->nDataNum <= 0)
					continue;

				vM = pData->vSeqData[ni]->vData[2]->pMatElement;
				for(nj=0; nj<pData->vSeqData[ni]->nDataNum; nj++)
				{
					vM[nj] = 0.0;
				}
			}
		}
	}

	else
	{
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			if(pData->vSeqData[ni]->nDataNum <= 0)
				continue;
			
			vM = pData->vSeqData[ni]->vData[2]->pMatElement;
			vN = pData->vSeqData[ni]->vData[3]->pMatElement;
			if(pData->vSeqData[ni]->vData[4] != NULL)
			{
				vMask = pData->vSeqData[ni]->vData[4]->pMatElement;
			}
			else
			{
				vMask = NULL;
			}

			if(vMask == NULL)
			{
				for(nj=0; nj<pData->vSeqData[ni]->nDataNum; nj++)
				{
					if(vN[nj] <= pParam->nW)
						vM[nj] = 0.0;
				}
			}
			else
			{
				for(nj=0; nj<pData->vSeqData[ni]->nDataNum; nj++)
				{
					if( (vMask[nj]>0.5) || (vN[nj] <= pParam->nW) )
						vM[nj] = 0.0;
				}
			}
		}
	}

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_MA_CallRegion()                                              */
/*  TileMapv2 call binding regions based on MA statistics.                 */
/* ----------------------------------------------------------------------- */
struct DOUBLEMATRIX *TileMapv2_MA_CallRegion(struct tagBARData *pData, 
					struct tagTileMapv2Param *pParam, int nLibId)
{
	/* define */
	struct DOUBLEMATRIX *pRegion = NULL;
	int nRegionCount;
	double dDist;
	double dSumT;
	double dMaxT,dMaxFC;
	int nMaxTPos,nMaxFCPos;
	int nProbeCount;
	int nStart,nEnd;
	int nLastRegionProbe;
	int nSeqId;
	double *vM,*vN,*vMask,*vFC,*vPos;
	int nGoodProbe;
	int nInterGoodProbe;
	int ni;
	
	/* define */
	FILE *fpOut;
	char strFileName[MED_LINE_LENGTH];
	
	/* init */
	sprintf(strFileName, "%s%s_%d.mareg.tmpout", pParam->strWorkPath, 
		pParam->strProjectTitle, nLibId+1);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: TileMapv2_MA_CallRegion, cannot open *.tmpout file for output!\n");
		exit(EXIT_FAILURE);
	}
	
	nRegionCount = 0;
	for(nSeqId=0; nSeqId<pData->nSeqNum; nSeqId++)
	{
		if(pData->vSeqData[nSeqId]->nDataNum <= 0)
			continue;

		dDist = 0.0;
		nProbeCount = 0;
		dSumT = 0.0;
		dMaxT = -1e20;
		dMaxFC = -1e20;
		nStart = -1;
		nEnd = -1;
		nMaxTPos = -1;
		nMaxFCPos = -1;
		nLastRegionProbe = -1;
		nInterGoodProbe = 0;

		vPos = pData->vSeqData[nSeqId]->vData[0]->pMatElement;
		vM = pData->vSeqData[nSeqId]->vData[2]->pMatElement;
		vN = pData->vSeqData[nSeqId]->vData[3]->pMatElement;
		if(pData->vSeqData[nSeqId]->vData[4] != NULL)
		{
			vMask = pData->vSeqData[nSeqId]->vData[4]->pMatElement;
		}
		else
		{
			vMask = NULL;
		}
		if(pData->vSeqData[nSeqId]->vData[5] != NULL)
		{
			vFC = pData->vSeqData[nSeqId]->vData[5]->pMatElement;
		}
		else
		{
			vFC = NULL;
		}
	
		for(ni=0; ni<pData->vSeqData[nSeqId]->nDataNum; ni++)
		{
			/* judge whether the probe is good or not */
			nGoodProbe = 0;
			if(vMask == NULL)
			{
				if(vN[ni] > pParam->nW)
					nGoodProbe = 1;
			}
			else
			{
				if( (vMask[ni] < 0.5) && (vN[ni] > pParam->nW) )
					nGoodProbe = 1;
			}

			if(nGoodProbe == 0)
				continue;

			/* if the probe does not pass the cutoff value */
			if(vM[ni] < pParam->dBaseStd)
			{
				nInterGoodProbe += 1;
			}
			/* if the probe passes the cutoff value */
			else
			{
				/* if not the first probe, does it belong to the previous region? */
				if(nLastRegionProbe >= 0)
				{
					dDist = vPos[ni]-vPos[nLastRegionProbe];
				}
				else
				{
					dDist = 1e20;
				}

				/* if it belongs to the previous region, update the previous region */
				if( (dDist <= pParam->nGap) && (nInterGoodProbe <= pParam->nGapW) )
				{
					nLastRegionProbe = ni;
					nEnd = (int)(vPos[ni]);
					nProbeCount += 1;
					dSumT += vM[ni];
					if(dMaxT < vM[ni])
					{
						dMaxT = vM[ni];
						nMaxTPos = (int)(vPos[ni]);
					}
					
					if(vFC == NULL)
					{
						dMaxFC = dMaxT;
						nMaxFCPos = nMaxTPos;
					}
					else
					{
						if(dMaxFC < vFC[ni])
						{
							dMaxFC = vFC[ni];
							nMaxFCPos = (int)(vPos[ni]);
						}
					}
				}

				/* if it is the first probe, or 
				   if it does not belong to the previous region, save the previous region
				   and initiate a new region */
				else
				{
					if( (nProbeCount >= pParam->nMinRegProbeNum) &&
						((nEnd-nStart+1) >= pParam->nMinRegLen) )
					{
						fprintf(fpOut, "%d\t%d\t%d\t% 9.7e\t%d\t% 9.7e\t%d\t% 9.7e\t%d\n", 
							nSeqId, nStart, nEnd, dMaxT, nMaxTPos,
							dMaxFC, nMaxFCPos, dSumT, nProbeCount);
						
						nRegionCount++;
						
						nProbeCount = 0;
						dSumT = 0.0;
						dMaxT = -1e20;
						dMaxFC = -1e20;
						nStart = -1;
						nEnd = -1;
						nMaxTPos = -1;
						nMaxFCPos = -1;
						nLastRegionProbe = -1;
					}

					nStart = (int)(vPos[ni]);
					nEnd = (int)(vPos[ni]);
					nProbeCount = 1;
					dSumT = vM[ni];
					dMaxT = vM[ni];
					nMaxTPos = (int)(vPos[ni]);
					if(vFC == NULL)
					{
						dMaxFC = dMaxT;
						nMaxFCPos = nMaxTPos;
					}
					else
					{
						dMaxFC = vFC[ni];
						nMaxFCPos = (int)(vPos[ni]);
					}
					nLastRegionProbe = ni;
				}

				/* reset intergoodprobe number */
				nInterGoodProbe = 0;
			}
		}
		
		if( (nProbeCount >= pParam->nMinRegProbeNum) &&
						((nEnd-nStart+1) >= pParam->nMinRegLen) )
		{
			fprintf(fpOut, "%d\t%d\t%d\t% 9.7e\t%d\t% 9.7e\t%d\t% 9.7e\t%d\n", 
						nSeqId, nStart, nEnd, dMaxT, nMaxTPos,
						dMaxFC, nMaxFCPos, dSumT, nProbeCount);

			nRegionCount++;
			nProbeCount = 0;
			dSumT = 0.0;
			dMaxT = -1e20;
			dMaxFC = -1e20;
			nStart = -1;
			nEnd = -1;
			nMaxTPos = -1;
			nMaxFCPos = -1;
			nLastRegionProbe = -1;
		}
	}

	/* close file */
	fclose(fpOut);

	/* merge signals */
	if(nRegionCount <= 0)
	{
		/* release memory */
		RemoveFiles(strFileName);
		printf("No regions found using the specified criteria!\n");
		return NULL;
	}

	pRegion = NULL;
	pRegion = DMLOAD(strFileName);
	if(pRegion == NULL)
	{
		printf("Error: TileMapv2_MA_CallRegion, cannot reload binding regions!\n");
		exit(EXIT_FAILURE);
	}
	if(pRegion->nHeight != nRegionCount)
	{
		printf("Error: TileMapv2_MA_CallRegion, cannot reload binding regions correctly!\n");
		exit(EXIT_FAILURE);
	}
	RemoveFiles(strFileName);

	/* return */
	return pRegion;
} 

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_RegionDetection_MA_ResetData()                               */
/*  TileMapv2 MA reset data.                                               */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_RegionDetection_MA_ResetData(struct tagBARData *pData, 
							struct INTMATRIX *pCol)
{
	/* define */
	int ni,nj;
	
	/* init */
	if( (pData == NULL) || (pCol == NULL) )
		return PROC_SUCCESS;

	if(pCol->nWidth != pData->nColNum)
	{
		printf("Error: TileMapv2_RegionDetection_MA_ResetData, column number not match!\n");
		exit(EXIT_FAILURE);
	}
	
	/* 1 position + 1 probe level statistics + 1 MA + 1 N + 1 mask + 1 Fold change */ 
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		if(pData->vSeqData[ni]->vData != NULL)
		{
			for(nj=0; nj<pCol->nWidth; nj++)
			{
				if( (pCol->pMatElement[nj] == 1) && (pData->vSeqData[ni]->vData[nj] != NULL) )
				{
					DestroyDoubleMatrix(pData->vSeqData[ni]->vData[nj]);
					pData->vSeqData[ni]->vData[nj] = NULL;
					pData->vSeqData[ni]->vData[nj] = CreateDoubleMatrix(1, pData->vSeqData[ni]->nDataNum);
					if(pData->vSeqData[ni]->vData[nj] == NULL)
					{
						printf("Error: TileMapv2_RegionDetection_MA_ResetData, cannot reset bar data!\n");
						exit(EXIT_FAILURE);
					}
				}
			}
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_RegionFDR_Count()                                            */
/*  Count random observations for computing FDR based on column nCol.      */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_RegionFDR_Count(struct DOUBLEMATRIX *pRegionSort, 
							  struct DOUBLEMATRIX *pRegionCTSort, 
							  struct DOUBLEMATRIX *pRegionFDR, int nCol)
{
	/* define */
	int ni,nj;
	double dObs,dExp;

	/* init check */
	if( (pRegionSort == NULL) || (pRegionCTSort == NULL) || (pRegionFDR == NULL) )
	{
		printf("Warning: TileMapv2_RegionFDR_Count, empty regions!\n");
		return PROC_SUCCESS;
	}

	if( pRegionSort->nHeight != pRegionFDR->nWidth )
	{
		printf("Error: TileMapv2_RegionFDR_Count, region number not match!\n");
		exit(EXIT_FAILURE);
	}

	if( (nCol<0) || (nCol>=pRegionSort->nWidth) )
	{
		printf("Error: TileMapv2_RegionFDR_Count, column index out of range!\n");
		exit(EXIT_FAILURE);
	}

	nj = 0;
	dExp = DMGETAT(pRegionCTSort, (pRegionCTSort->nHeight-1-nj), nCol);
	for(ni=0; ni<pRegionSort->nHeight; ni++)
	{
		dObs = DMGETAT(pRegionSort, (pRegionSort->nHeight-1-ni), nCol);
		while( (dExp >= dObs) && (nj<pRegionCTSort->nHeight) )
		{
			nj++;
			if(nj == pRegionCTSort->nHeight)
				break;
			dExp = DMGETAT(pRegionCTSort, (pRegionCTSort->nHeight-1-nj), nCol);
		}

		pRegionFDR->pMatElement[pRegionFDR->nWidth-1-ni] += nj;
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_RegionFDR_Compute()                                          */
/*  Compute FDR from the count.                                            */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_RegionFDR_Compute(struct DOUBLEMATRIX *pFDR)
{
	/* define */
	int ni;

	/* compute */
	if(pFDR == NULL)
		return PROC_SUCCESS;

	for(ni=0; ni<pFDR->nWidth; ni++)
	{
		pFDR->pMatElement[ni] = pFDR->pMatElement[ni]/(double)(pFDR->nWidth-ni);
		if(pFDR->pMatElement[ni] > 1.0)
			pFDR->pMatElement[ni] = 1.0;

		if(ni>0)
		{
			if(pFDR->pMatElement[ni] > pFDR->pMatElement[ni-1])
			{
				pFDR->pMatElement[ni] = pFDR->pMatElement[ni-1];
			}
		}
	}

	/* return */
	return PROC_SUCCESS;
}
			
/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_RegionFDR_Reorganize()                                       */
/*  Reorganize regions by adding FDR.                                      */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *TileMapv2_RegionFDR_Reorganize(struct DOUBLEMATRIX *pRegion, 
				struct DOUBLEMATRIX *pRegionMaxFDR, struct LONGMATRIX *pRegionMaxSid, 
				struct DOUBLEMATRIX *pRegionSumFDR, struct LONGMATRIX *pRegionSumSid)
{
	/* define */
	struct DOUBLEMATRIX *pNewRegion = NULL;
	int ni,nj;
	double dTemp;
	int nRow;
	struct DOUBLEMATRIX *pRegionMaxLfdr = NULL;
	struct DOUBLEMATRIX *pRegionSumLfdr = NULL;

	/* init */
	if(pRegion == NULL)
		return NULL;

	if( (pRegionMaxFDR==NULL) || (pRegionMaxSid==NULL) 
		|| (pRegionSumFDR == NULL) || (pRegionSumSid==NULL) )
	{
		printf("Error: TileMapv2_RegionFDR_Reorganize, no fdr available!\n");
		exit(EXIT_FAILURE);
	}

	pRegionMaxLfdr = TileMapv2_RegionFDR_LocalFDR(pRegionMaxFDR);
	pRegionSumLfdr = TileMapv2_RegionFDR_LocalFDR(pRegionSumFDR);

	if( (pRegionMaxLfdr==NULL) || (pRegionSumLfdr==NULL) )
	{
		printf("Error: TileMapv2_RegionFDR_Reorganize, no lfdr available!\n");
		exit(EXIT_FAILURE);
	}

	pNewRegion = CreateDoubleMatrix(pRegion->nHeight, pRegion->nWidth+4);
	if(pNewRegion == NULL)
	{
		printf("Error: TileMapv2_RegionFDR_Reorganize, cannot allocate memory for new regions!\n");
		exit(EXIT_FAILURE);
	}

	/* remap */
	for(ni=0; ni<pRegion->nHeight; ni++)
	{
		/* id, start, end, maxT, maxTPos */
		for(nj=0; nj<5; nj++)
		{
			dTemp = DMGETAT(pRegion, ni, nj);
			DMSETAT(pNewRegion, ni, nj, dTemp);
		}

		/* maxFC, maxFC Pos, sumMA, sumMA probe count */
		for(; nj<pRegion->nWidth; nj++)
		{
			dTemp = DMGETAT(pRegion, ni, nj);
			DMSETAT(pNewRegion, ni, (nj+2), dTemp);
		}
	}

	for(ni=0; ni<pRegionMaxSid->nHeight; ni++)
	{
		nRow = pRegionMaxSid->pMatElement[ni];
		dTemp = pRegionMaxFDR->pMatElement[ni];
		DMSETAT(pNewRegion, nRow, 5, dTemp);
		dTemp = pRegionMaxLfdr->pMatElement[ni];
		DMSETAT(pNewRegion, nRow, 6, dTemp);
	}

	for(ni=0; ni<pRegionSumSid->nHeight; ni++)
	{
		nRow = pRegionSumSid->pMatElement[ni];
		dTemp = pRegionSumFDR->pMatElement[ni];
		DMSETAT(pNewRegion, nRow, (pRegion->nWidth+2), dTemp);
		dTemp = pRegionSumLfdr->pMatElement[ni];
		DMSETAT(pNewRegion, nRow, (pRegion->nWidth+3), dTemp);
	}

	/* release memory */
	DestroyDoubleMatrix(pRegionMaxLfdr);
	DestroyDoubleMatrix(pRegionSumLfdr);

	/* return */
	return pNewRegion;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_RegionFDR_LocalFDR()                                         */
/*  Compute local fdr.                                                     */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX* TileMapv2_RegionFDR_LocalFDR(struct DOUBLEMATRIX* pFDR)
{
	/* define */
	struct DOUBLEMATRIX *pLfdr = NULL;
	double dStep = 0.01;
	double dCurrent;
	int nStep = 100;
	int ni,nj;
	int nP1;
	double dlfdr;
	double dLastlfdr;

	/* init */
	if(pFDR == NULL)
		return NULL;

	/* compute lfdr */
	pLfdr = CreateDoubleMatrix(pFDR->nHeight, pFDR->nWidth);
	if(pLfdr == NULL)
	{
		printf("Error: TileMapv2_RegionFDR_LocalFDR, cannot compute lfdr!\n");
		exit(EXIT_FAILURE);
	}

	nP1 = 0;
	dCurrent = 1.0;
	dLastlfdr = 1.0;
	for(ni=0; ni<pFDR->nWidth; ni++)
	{
		if(pFDR->pMatElement[ni]<dCurrent)
		{
			if( (ni-nP1) >= nStep)
			{
				dlfdr = ((pFDR->nWidth-nP1)*pFDR->pMatElement[nP1]-(pFDR->nWidth-ni)*pFDR->pMatElement[ni])/(double)(ni-nP1);
				if(dlfdr > dLastlfdr)
					dlfdr = dLastlfdr;

				for(nj=nP1; nj<ni; nj++)
				{
					pLfdr->pMatElement[nj] = dlfdr;
				}

				dLastlfdr = dlfdr;
				nP1 = ni;
			}

			while(pFDR->pMatElement[ni]<dCurrent)
			{
				dCurrent -= dStep;
			}
		}
	}

	dlfdr = pFDR->pMatElement[nP1];
	if(dlfdr > dLastlfdr)
		dlfdr = dLastlfdr;
	for(nj=nP1; nj<pFDR->nWidth; nj++)
	{
		pLfdr->pMatElement[nj] = dlfdr;
	}
	
	/* return */
	return pLfdr;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_MA_ExportRegion()                                            */
/*  Export regions.                                                        */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_MA_ExportRegion(struct DOUBLEMATRIX *pRegion, struct tagBARData *pData, 
							  struct tagTileMapv2Param *pParam, int nLibId)
{
	/* define */
	FILE *fpOut;
	char strFileName[MED_LINE_LENGTH];
	int ni,nSeqId;
	char strSeqName[MED_LINE_LENGTH];
	char strGroupName[MED_LINE_LENGTH];
	int nStart,nEnd,nLen;
	double dMaxT,dMaxTFDR,dMaxFC,dMaxSum,dMaxSumFDR;
	double dMaxTLfdr,dMaxSumLfdr;
	int nMaxTPos,nMaxFCPos,nSumProbe;

	/* init */
	sprintf(strFileName, "%s%s_%d.cod", pParam->strWorkPath, pParam->strProjectTitle, (nLibId+1));
	if(pRegion == NULL)
	{
		fpOut = NULL;
		fpOut = fopen(strFileName, "w");
		if(fpOut == NULL)
		{
			printf("Error: TileMapv2_MA_ExportRegion, cannot open the output file!\n");
			exit(EXIT_FAILURE);
		}

		fclose(fpOut);
		return PROC_SUCCESS;
	}

	if(pData == NULL)
	{
		printf("Error: TileMapv2_MA_ExportRegion, empty data!\n");
		exit(EXIT_FAILURE);
	}

	/* open file */
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: TileMapv2_MA_ExportRegion, cannot open the output file!\n");
		exit(EXIT_FAILURE);
	}
	
	/* write */
	for(ni=0; ni<pRegion->nHeight; ni++)
	{
		nSeqId = (int)(DMGETAT(pRegion, ni, 0));

		if(pData->vSeqData[nSeqId]->pSeqVersion != NULL)
		{
			sprintf(strSeqName, "%s", pData->vSeqData[nSeqId]->pSeqName->m_pString);
			sprintf(strGroupName, "%s", pData->vSeqData[nSeqId]->pSeqVersion->m_pString);
		}
		else
		{
			sprintf(strSeqName, "%s", pData->vSeqData[nSeqId]->pSeqName->m_pString);
			sprintf(strGroupName, "NA");
		}
		if(pData->fVersionnumber > 1.5)
		{
			if(pData->vSeqData[nSeqId]->pSeqGroupName != NULL)
			{
				if(pData->vSeqData[nSeqId]->pSeqVersion != NULL)
				{
					sprintf(strSeqName, "%s", pData->vSeqData[nSeqId]->pSeqName->m_pString);
					sprintf(strGroupName, "%s:%s", pData->vSeqData[nSeqId]->pSeqGroupName->m_pString,
						pData->vSeqData[nSeqId]->pSeqVersion->m_pString);
				}
				else
				{
					sprintf(strSeqName, "%s", pData->vSeqData[nSeqId]->pSeqName->m_pString);
					sprintf(strGroupName, "%s", pData->vSeqData[nSeqId]->pSeqGroupName->m_pString);
				}
			}	
		}

		nStart = (int)(DMGETAT(pRegion, ni, 1));
		nEnd = (int)(DMGETAT(pRegion, ni, 2));
		nLen = nEnd-nStart+1;

		
		dMaxT = DMGETAT(pRegion, ni, 3);
		nMaxTPos = (int)(DMGETAT(pRegion, ni, 4));
		dMaxTFDR = DMGETAT(pRegion, ni, 5);
		dMaxTLfdr = DMGETAT(pRegion, ni, 6);
		dMaxFC = DMGETAT(pRegion, ni, 7);
		nMaxFCPos = (int)(DMGETAT(pRegion, ni, 8));
		dMaxSum = DMGETAT(pRegion, ni, 9);
		nSumProbe = (int)(DMGETAT(pRegion, ni, 10));
		dMaxSumFDR = DMGETAT(pRegion, ni, 11);
		dMaxSumLfdr = DMGETAT(pRegion, ni, 12);

		/* seqid, chr, start, end, strand, length, maxMA, maxMA Pos, maxMA FDR, maxFC, maxFC Pos,
			sumMA, sumMA good probe count, sumMA FDR, libid */
		fprintf(fpOut, "%d\t%s\t%d\t%d\t+\t%d\t%f\t%d\t%f\t%f\t%f\t%d\t%f\t%d\t%f\t%f\t%d\t%s\n",
			ni, strSeqName, nStart, nEnd, nLen, dMaxT, nMaxTPos, dMaxTFDR, dMaxTLfdr,
			dMaxFC, nMaxFCPos, dMaxSum, nSumProbe, dMaxSumFDR, dMaxSumLfdr, (nLibId+1), strGroupName);
	}

	/* close file */
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_MA_RegionFDR_Permutation()                                   */
/*  TileMap Region FDR.                                                    */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_MA_RegionFDR_Permutation(struct INTMATRIX *pType,
				struct DOUBLEMATRIX *pRegionMaxSort, struct DOUBLEMATRIX *pRegionMaxFDR,
				struct DOUBLEMATRIX *pRegionSumSort, struct DOUBLEMATRIX *pRegionSumFDR, 
				struct tagTileMapv2Param *pParam, int nLibId)
{
	/* define */
	char strProbeFile[MED_LINE_LENGTH];
	char strMaskFile[MED_LINE_LENGTH];
	char strFCFile[MED_LINE_LENGTH];
	struct tagBARData *pData = NULL;

	struct INTMATRIX *pPriority = NULL;
	struct DOUBLEMATRIX *pRegionCT = NULL;
	struct DOUBLEMATRIX *pRegionCTMaxSort = NULL;
	struct LONGMATRIX *pRegionCTMaxSid = NULL;
	struct DOUBLEMATRIX *pRegionCTSumSort = NULL;
	struct LONGMATRIX *pRegionCTSumSid = NULL;
	int nCol;
	int nPermId;
	struct INTMATRIX *pCol = NULL;
	

	/* permutation */
	for(nPermId=0; nPermId<pParam->nPermutationNum; nPermId++)
	{
		printf("perm %d...\n", nPermId);

		/* ------------- */
		/* load raw data */
		/* ------------- */
		strcpy(strMaskFile, "");
		strcpy(strFCFile, "");
		sprintf(strProbeFile, "%s%s_%d.perm%d.pb.bar", 
			pParam->strWorkPath, pParam->strProjectTitle, nLibId+1, nPermId);
		if(pParam->nNoiseMask == 1)
		{
			sprintf(strMaskFile, "%s%s_%d.perm%d.pb.bar.mask", 
				pParam->strWorkPath, pParam->strProjectTitle, nLibId+1, nPermId);
		}
		if(pParam->nPatternType == 2)
		{
			sprintf(strFCFile, "%s%s_%d.perm%d.fc.bar", 
				pParam->strWorkPath, pParam->strProjectTitle, nLibId+1, nPermId);
		}
		
		pData = NULL;
		pData = TileMapv2_RegionDetection_MA_PrepareData(strProbeFile, strMaskFile, strFCFile);
		if(pData == NULL)
		{
			printf("Error: TileMapv2_MA_RegionFDR_Permutation, cannot load raw data!\n");
			exit(EXIT_FAILURE);
		}

		/* -------------------- */
		/* MA                   */
		/* -------------------- */
		/* MA */
		TileMapv2_MA(pData, pParam);

		/* Call binding regions */
		pRegionCT = NULL;
		pRegionCT = TileMapv2_MA_CallRegion(pData, pParam, nLibId);

		if(pRegionCT != NULL)
		{
			/* sort region, regionCT */
			
			/* maxT */
			pPriority = NULL;
			pPriority = CreateIntMatrix(1, 1);
			if(pPriority == NULL)
			{
				printf("Error: TileMapv2_MA_RegionFDR_Permutation, cannot create sorting key!\n");
				exit(EXIT_FAILURE);
			}
			pPriority->pMatElement[0] = 3;
			pRegionCTMaxSort = NULL;
			pRegionCTMaxSid = NULL;
			DMSORTROWS(pRegionCT, pType, pPriority, &pRegionCTMaxSort, &pRegionCTMaxSid);

			/* sumT */
			pPriority->pMatElement[0] = 7;
			pRegionCTSumSort = NULL;
			pRegionCTSumSid = NULL;
			DMSORTROWS(pRegionCT, pType, pPriority, &pRegionCTSumSort, &pRegionCTSumSid);
			
			DestroyIntMatrix(pPriority);

			/* count */
			nCol = 3;
			TileMapv2_RegionFDR_Count(pRegionMaxSort, pRegionCTMaxSort, pRegionMaxFDR, nCol);

			nCol = 7;
			TileMapv2_RegionFDR_Count(pRegionSumSort, pRegionCTSumSort, pRegionSumFDR, nCol);

			/* -------------- */
			/* release memory */
			/* -------------- */
			DestroyDoubleMatrix(pRegionCTMaxSort);
			DestroyDoubleMatrix(pRegionCTSumSort);
			DestroyLongMatrix(pRegionCTMaxSid);
			DestroyLongMatrix(pRegionCTSumSid);
			DestroyDoubleMatrix(pRegionCT);
		}
		else
		{
			/* FDR == 0*/
		}

		/* -------------------- */
		/* Destroy              */
		/* -------------------- */
		Affy_BARData_Destroy(&pData);
		RemoveFiles(strProbeFile);
		if(strcmp(strMaskFile, "") != 0)
			RemoveFiles(strMaskFile);
		if(strcmp(strFCFile, "") != 0)
			RemoveFiles(strFCFile);
	}

	/* FDR */
	DMPDIVTS(pRegionMaxFDR, (double)pParam->nPermutationNum);
	DMPDIVTS(pRegionSumFDR, (double)pParam->nPermutationNum);
	TileMapv2_RegionFDR_Compute(pRegionMaxFDR);
	TileMapv2_RegionFDR_Compute(pRegionSumFDR);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_MA_RegionFDR_LeftTail()                                      */
/*  TileMap Region FDR.                                                    */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_MA_RegionFDR_LeftTail(struct INTMATRIX *pType,
				struct DOUBLEMATRIX *pRegionMaxSort, struct DOUBLEMATRIX *pRegionMaxFDR,
				struct DOUBLEMATRIX *pRegionSumSort, struct DOUBLEMATRIX *pRegionSumFDR, 
				struct tagBARData *pData, struct tagTileMapv2Param *pParam, int nLibId)
{
	/* define */
	struct INTMATRIX *pCol = NULL;
	struct INTMATRIX *pPriority = NULL;

	struct DOUBLEMATRIX *pRegionCT = NULL;
	struct DOUBLEMATRIX *pRegionCTMaxSort = NULL;
	struct LONGMATRIX *pRegionCTMaxSid = NULL;
	struct DOUBLEMATRIX *pRegionCTSumSort = NULL;
	struct LONGMATRIX *pRegionCTSumSid = NULL;
	int nCol;

	printf("FDR from left tail ...\n");

	pCol = NULL;
	pCol = CreateIntMatrix(1, pData->nColNum);
	pCol->pMatElement[2] = 1;
	pCol->pMatElement[5] = 1;
	if(pCol == NULL)
	{
		printf("Error: TileMapv2_RegionFDR_CountFromLeftTail, cannot create column information!\n");
		exit(EXIT_FAILURE);
	}

	Affy_BARData_SubTS(pData, 0.0, pCol);
	DestroyIntMatrix(pCol);
	
	/* call regions */
	pRegionCT = NULL;
	pRegionCT = TileMapv2_MA_CallRegion(pData, pParam, nLibId);

    if(pRegionCT != NULL)
	{
		/* sort region, regionCT */
		
		/* maxT */
		pPriority = NULL;
		pPriority = CreateIntMatrix(1, 1);
		if(pPriority == NULL)
		{
			printf("Error: TileMapv2_RegionFDR_CountFromLeftTail, cannot create sorting key!\n");
			exit(EXIT_FAILURE);
		}
		pPriority->pMatElement[0] = 3;
		pRegionCTMaxSort = NULL;
		pRegionCTMaxSid = NULL;
		DMSORTROWS(pRegionCT, pType, pPriority, &pRegionCTMaxSort, &pRegionCTMaxSid);

		/* sumT */
		pPriority->pMatElement[0] = 7;
		pRegionCTSumSort = NULL;
		pRegionCTSumSid = NULL;
		DMSORTROWS(pRegionCT, pType, pPriority, &pRegionCTSumSort, &pRegionCTSumSid);
		
		DestroyIntMatrix(pPriority);

		/* count */
		nCol = 3;
		TileMapv2_RegionFDR_Count(pRegionMaxSort, pRegionCTMaxSort, pRegionMaxFDR, nCol);

		nCol = 7;
		TileMapv2_RegionFDR_Count(pRegionSumSort, pRegionCTSumSort, pRegionSumFDR, nCol);

		/* -------------- */
		/* release memory */
		/* -------------- */
		DestroyDoubleMatrix(pRegionCTMaxSort);
		DestroyDoubleMatrix(pRegionCTSumSort);
		DestroyLongMatrix(pRegionCTMaxSid);
		DestroyLongMatrix(pRegionCTSumSid);
		DestroyDoubleMatrix(pRegionCT);
	}
	else
	{
		/* FDR == 0*/
	}

	/* FDR */
	TileMapv2_RegionFDR_Compute(pRegionMaxFDR);
	TileMapv2_RegionFDR_Compute(pRegionSumFDR);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_MA_RegionFDR_UMS()                                           */
/*  TileMap Region FDR: UMS                                                */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *TileMapv2_MA_RegionFDR_UMS(struct DOUBLEMATRIX *pRegion,
				struct tagBARData *pData, struct tagTileMapv2Param *pParam, 
				int nLibId)
{
	/* define */
	struct DOUBLEMATRIX *pNewRegion;
	struct DOUBLEMATRIX *pSelect;
	struct DOUBLEMATRIX *pScore;
	struct DOUBLEMATRIX *pTemp;
	struct DOUBLEMATRIX *pIntx;
	struct DOUBLEMATRIX *pF0;
	struct DOUBLEMATRIX *pF1;
	struct DOUBLEMATRIX *pLfdr;
	double dTheta = 0.0;
	double dTemp;
	double dSum;

	int nProbeNum,nTotalProbeNum;
	int ni,nj,nk;
	double *vT,*vM,*vN,*vMask;

	/* init */
	printf("UMS ...\n");

	/* get probe number */
	nProbeNum = 0;
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		nProbeNum += pData->vSeqData[ni]->nDataNum;
	}

	if(nProbeNum <= 0)
	{
		printf("Error: TileMapv2_MA_RegionFDR_UMS, empty MA input!\n");
		exit(EXIT_FAILURE);
	}
	
	pSelect = NULL;
	pSelect = CreateDoubleMatrix(1, nProbeNum);
	if(pSelect == NULL)
	{
		printf("Error: TileMapv2_MA_RegionFDR_UMS, cannot create selection matrix!\n");
		exit(EXIT_FAILURE);
	}

	pScore = NULL;
	pScore = CreateDoubleMatrix(1, nProbeNum);
	if(pScore == NULL)
	{
		printf("Error: TileMapv2_MA_RegionFDR_UMS, cannot create score matrix!\n");
		exit(EXIT_FAILURE);
	}

	/* load data */
	if(pParam->nNoiseMask == 1)
	{
		nTotalProbeNum = 0;
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			if(pData->vSeqData[ni]->nDataNum <= 0)
				continue;
			
			vT = pData->vSeqData[ni]->vData[1]->pMatElement;
			vM = pData->vSeqData[ni]->vData[2]->pMatElement;
			vN = pData->vSeqData[ni]->vData[3]->pMatElement;
			vMask = pData->vSeqData[ni]->vData[4]->pMatElement;
			for(nj=0; nj<pData->vSeqData[ni]->nDataNum; nj++)
			{
				if( (vMask[nj] < 0.5) && (vN[nj] > pParam->nW) )
				{
					pSelect->pMatElement[nTotalProbeNum] = vT[nj];
					pScore->pMatElement[nTotalProbeNum] = vM[nj];
					nTotalProbeNum++;
				}
			}
		}

		pTemp = NULL;
		pTemp = CreateDoubleMatrix(1, nTotalProbeNum);
		if(pTemp == NULL)
		{
			printf("Error: TileMapv2_MA_RegionFDR_UMS, cannot create selection matrix!\n");
			exit(EXIT_FAILURE);
		}

		memcpy(pTemp->pMatElement, pSelect->pMatElement, sizeof(double)*nTotalProbeNum);
		DestroyDoubleMatrix(pSelect);
		pSelect = pTemp;

		pTemp = NULL;
		pTemp = CreateDoubleMatrix(1, nTotalProbeNum);
		if(pTemp == NULL)
		{
			printf("Error: TileMapv2_MA_RegionFDR_UMS, cannot create score matrix!\n");
			exit(EXIT_FAILURE);
		}

		memcpy(pTemp->pMatElement, pScore->pMatElement, sizeof(double)*nTotalProbeNum);
		DestroyDoubleMatrix(pScore);
		pScore = pTemp;
	}
	else
	{
		nTotalProbeNum = 0;
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			if(pData->vSeqData[ni]->nDataNum <= 0)
				continue;
			
			memcpy(pSelect->pMatElement+nTotalProbeNum, pData->vSeqData[ni]->vData[1]->pMatElement, sizeof(double)*pData->vSeqData[ni]->nDataNum);
			memcpy(pScore->pMatElement+nTotalProbeNum, pData->vSeqData[ni]->vData[2]->pMatElement, sizeof(double)*pData->vSeqData[ni]->nDataNum);
			nTotalProbeNum += pData->vSeqData[ni]->nDataNum;
		}

		if(nTotalProbeNum != nProbeNum)
		{
			printf("Error: TileMapv2_MA_RegionFDR_UMS, probe number not match!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* transformation */
	for(ni=0; ni<nTotalProbeNum; ni++)
	{
		dTemp = exp(pScore->pMatElement[ni]);
		pScore->pMatElement[ni] = 1.0/(1.0+dTemp);
		pSelect->pMatElement[ni] = -pSelect->pMatElement[ni];
	}


	/* UMS */
	pIntx = NULL;
	pIntx = CreateDoubleMatrix(1, pParam->nGridSize);
	if(pIntx == NULL)
	{
		printf("Error: TileMapv2_MA_RegionFDR_UMS, cannot create interval matrix!\n");
		exit(EXIT_FAILURE);
	}

	pF0 = NULL;
	pF0 = CreateDoubleMatrix(1, pParam->nGridSize);
	if(pF0 == NULL)
	{
		printf("Error: TileMapv2_MA_RegionFDR_UMS, cannot create f0 matrix!\n");
		exit(EXIT_FAILURE);
	}

	pF1 = NULL;
	pF1 = CreateDoubleMatrix(1, pParam->nGridSize);
	if(pF1 == NULL)
	{
		printf("Error: TileMapv2_MA_RegionFDR_UMS, cannot create f1 matrix!\n");
		exit(EXIT_FAILURE);
	}

	pLfdr = NULL;
	pLfdr = CreateDoubleMatrix(1, pParam->nGridSize);
	if(pLfdr == NULL)
	{
		printf("Error: TileMapv2_MA_RegionFDR_UMS, cannot create memory for UMS!\n");
		exit(EXIT_FAILURE);
	}
	
	dTheta = 0.0;
	TileMap_UMS(pSelect, pParam->dTp, pParam->dTq, pParam->nOffset, pScore, pParam->nGridSize, pIntx, pF0, pF1, &dTheta);

	/* get FDR */
	for(ni=0; ni<pParam->nGridSize; ni++)
	{
		dSum = dTheta*pF1->pMatElement[ni]+(1.0-dTheta)*pF0->pMatElement[ni];
		if(dSum <= 0.0)
		{
			printf("Error: TileMapv2_MA_RegionFDR_UMS, divide by zero when compute lfdr!\n");
			exit(EXIT_FAILURE);
		}
		pLfdr->pMatElement[ni] = 1.0-dTheta*pF1->pMatElement[ni]/dSum;
	}

	DestroyDoubleMatrix(pSelect);
	DestroyDoubleMatrix(pScore);

	/* assign fdr to the original regions */
	pNewRegion = NULL;
	pNewRegion = CreateDoubleMatrix(pRegion->nHeight, pRegion->nWidth+4);
	if(pNewRegion == NULL)
	{
		printf("Error: TileMapv2_MA_RegionFDR_UMS, cannot allocate memory for new regions!\n");
		exit(EXIT_FAILURE);
	}

	/* remap */
	for(ni=0; ni<pRegion->nHeight; ni++)
	{
		/* id, start, end, maxT, maxTPos */
		for(nj=0; nj<5; nj++)
		{
			dTemp = DMGETAT(pRegion, ni, nj);
			DMSETAT(pNewRegion, ni, nj, dTemp);
		}

		dTemp = DMGETAT(pRegion, ni, 3);
		dTemp = exp(dTemp);
		dTemp = 1.0/(1.0+dTemp);
		nk = (int)(dTemp*pParam->nGridSize);
		if(nk == pParam->nGridSize)
			nk--;
		dTemp = pLfdr->pMatElement[nk];
		DMSETAT(pNewRegion, ni, 5, dTemp);
		DMSETAT(pNewRegion, ni, 6, dTemp);

		/* maxFC, maxFC Pos, sumMA, sumMA probe count */
		for(; nj<pRegion->nWidth; nj++)
		{
			dTemp = DMGETAT(pRegion, ni, nj);
			DMSETAT(pNewRegion, ni, (nj+2), dTemp);
		}
	}

	/* destroy memory */
	DestroyDoubleMatrix(pIntx);
	DestroyDoubleMatrix(pF0);
	DestroyDoubleMatrix(pF1);
	DestroyDoubleMatrix(pLfdr);

	/* return */
	return pNewRegion;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_SummarizeResults_Main()                                      */
/*  Summarize region detection results.                                    */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_SummarizeResults_Main(struct tagTileMapv2Param *pParam, int nRegionNum)
{
	/* define */
	struct tagString **vReg = NULL;
	struct tagString **vGroup = NULL;
	int ni,nRowId;
	int nLibId;
	char strFileName[MED_LINE_LENGTH];
	char strLine[MED_LINE_LENGTH];
	FILE *fpIn,*fpOut;
	struct DOUBLEMATRIX *pRegion = NULL;
	struct DOUBLEMATRIX *pRegionSort = NULL;
	struct LONGMATRIX *pRegionSid = NULL;
	struct INTMATRIX *pType = NULL;
	struct INTMATRIX *pPriority = NULL;
	char strId[LINE_LENGTH];
	char strChr[MED_LINE_LENGTH];
	char strGroupName[MED_LINE_LENGTH];
	int nStart, nEnd, nMaxTPos, nMaxFCPos;
	char chStrand;
	int nLen,nSumProbe,nArrayId;
	double dMaxT,dMaxTFDR,dMaxFC,dSumT,dSumTFDR,dMaxTLfdr,dSumTLfdr;
	
	/* init */
	printf("Total Number of Regions Detected = %d\n", nRegionNum);
	if(nRegionNum <= 0)
	{
		sprintf(strFileName, "%s%s_all.cod", pParam->strWorkPath, pParam->strProjectTitle);
		fpOut = NULL;
		fpOut = fopen(strFileName, "w");
		if(fpOut == NULL)
		{
			printf("Error: TileMapv2_SummarizeResults_Main, cannot open the output file!\n");
			exit(EXIT_FAILURE);
		}

		fprintf(fpOut, "#rank\tchromosome\tstart\tend\tstrand\tregion_length\tmaxM/P\tposition_of_maxM/P\tFDR\tlocal_FDR\tmaxFC(log2)\tposition_of_maxFC\tsumM/P\tno_of_goodprobes\tsumM/P_FDR\tsumM/P_local_FDR\tlibrary_id\tgroup_name\n");
		fclose(fpOut);
		return PROC_SUCCESS;
	}

	vReg = NULL;
	vReg = (struct tagString **)calloc(nRegionNum, sizeof(struct tagString *));
	if(vReg == NULL)
	{
		printf("Error: TileMapv2_SummarizeResults_Main, cannot allocate memory for summarizing results!\n");
		exit(EXIT_FAILURE);
	}

	vGroup = NULL;
	vGroup = (struct tagString **)calloc(nRegionNum, sizeof(struct tagString *));
	if(vGroup == NULL)
	{
		printf("Error: TileMapv2_SummarizeResults_Main, cannot allocate memory for summarizing results!\n");
		exit(EXIT_FAILURE);
	}

	pRegion = NULL;
	pRegion = CreateDoubleMatrix(nRegionNum, 16);
	if(pRegion == NULL)
	{
		printf("Error: TileMapv2_SummarizeResults_Main, cannot allocate memory for summarizing results!\n");
		exit(EXIT_FAILURE);
	}

	/* "%s %d %d %c %d %lf %d %lf %lf %lf %d %lf %d %lf %lf %d",
				strChr, &nStart, &nEnd, &chStrand, &nLen, &dMaxT, &nMaxTPos,
				&dMaxTFDR, &dMaxTLfdr, &dMaxFC, &nMaxFCPos, &dSumT, &nSumProbe, &dSumTFDR, &dSumTLfdr, &nArrayId);
	*/
	pType = NULL;
	pType = CreateIntMatrix(1, 16);
	if(pType == NULL)
	{
		printf("Error: TileMapv2_SummarizeResults_Main, cannot allocate memory for summarizing results!\n");
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
	pType->pMatElement[8] = 1;
	pType->pMatElement[9] = 1;
	pType->pMatElement[10] = 2;
	pType->pMatElement[11] = 1;
	pType->pMatElement[12] = 2;
	pType->pMatElement[13] = 1;
	pType->pMatElement[14] = 1;
	pType->pMatElement[15] = 2;

	pPriority = NULL;
	pPriority = CreateIntMatrix(1, 4);
	if(pPriority == NULL)
	{
		printf("Error: TileMapv2_SummarizeResults_Main, cannot allocate memory for summarizing results!\n");
		exit(EXIT_FAILURE);
	}
	pPriority->pMatElement[0] = 7;
	pPriority->pMatElement[1] = 5;
	pPriority->pMatElement[2] = 9;
	pPriority->pMatElement[3] = 15;

	/* load regions */
	ni = 0;
	for(nLibId=0; nLibId<pParam->nLibNum; nLibId++)
	{
		sprintf(strFileName, "%s%s_%d.cod", pParam->strWorkPath, pParam->strProjectTitle, nLibId+1);
		fpIn = NULL;
		fpIn = fopen(strFileName, "r");
		if(fpIn == NULL)
			continue;

		while(fgets(strLine, MED_LINE_LENGTH, fpIn)!=NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;

			if(ni >= nRegionNum)
			{
				printf("Error: TileMapv2_SummarizeResults_Main, region number inconsistent!\n");
				exit(EXIT_FAILURE);
			}

			sscanf(strLine, "%s %s %d %d %c %d %lf %d %lf %lf %lf %d %lf %d %lf %lf %d %s",
				strId, strChr, &nStart, &nEnd, &chStrand, &nLen, &dMaxT, &nMaxTPos,
				&dMaxTFDR, &dMaxTLfdr, &dMaxFC, &nMaxFCPos, &dSumT, &nSumProbe, 
				&dSumTFDR, &dSumTLfdr, &nArrayId, strGroupName);

			DMSETAT(pRegion, ni, 1, nStart);
			DMSETAT(pRegion, ni, 2, nEnd);
			DMSETAT(pRegion, ni, 3, -nLen);
			if(chStrand == '+')
				DMSETAT(pRegion, ni, 4, 0.0);
			else
				DMSETAT(pRegion, ni, 4, 1.0);
			DMSETAT(pRegion, ni, 5, -dMaxT);
			DMSETAT(pRegion, ni, 6, nMaxTPos);
			DMSETAT(pRegion, ni, 7, dMaxTFDR);
			DMSETAT(pRegion, ni, 8, dMaxTLfdr);
			DMSETAT(pRegion, ni, 9, -dMaxFC);
			DMSETAT(pRegion, ni, 10, nMaxFCPos);
			DMSETAT(pRegion, ni, 11, -dSumT);
			DMSETAT(pRegion, ni, 12, -nSumProbe);
			DMSETAT(pRegion, ni, 13, dSumTFDR);
			DMSETAT(pRegion, ni, 14, dSumTLfdr);
			DMSETAT(pRegion, ni, 15, nArrayId);

			StringAddTail(vReg+ni, strChr);
			StringAddTail(vGroup+ni, strGroupName);
			ni++;
		}

		fclose(fpIn);
	}

	if(ni != nRegionNum)
	{
		printf("Error: TileMapv2_SummarizeResults_Main, region number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* sorting according to maxT */
	pRegionSort = NULL;
	pRegionSid = NULL;
	DMSORTROWS(pRegion, pType, pPriority, &pRegionSort, &pRegionSid);


	/* output */
	sprintf(strFileName, "%s%s_all.cod", pParam->strWorkPath, pParam->strProjectTitle);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: TileMapv2_SummarizeResults_Main, cannot open the output file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "#rank\tchromosome\tstart\tend\tstrand\tregion_length\tmaxM/P\tposition_of_maxM/P\tFDR\tlocal_FDR\tmaxFC(log2)\tposition_of_maxFC\tsumM/P\tno_of_goodprobes\tsumM/P_FDR\tsumM/P_local_FDR\tlibrary_id\tgroup_name\n");

	for(ni=0; ni<nRegionNum; ni++)
	{
		nRowId = pRegionSid->pMatElement[ni];
		strcpy(strChr, vReg[nRowId]->m_pString);
		strcpy(strGroupName, vGroup[nRowId]->m_pString);

		nStart = (int)(DMGETAT(pRegionSort, ni, 1));
		nEnd = (int)(DMGETAT(pRegionSort, ni, 2));
		nLen = -(int)(DMGETAT(pRegionSort, ni, 3));
		if( DMGETAT(pRegionSort, ni, 4) > 0.5)
			chStrand = '-';
		else
			chStrand = '+';
		dMaxT = -DMGETAT(pRegionSort, ni, 5);
		nMaxTPos = (int)(DMGETAT(pRegionSort, ni, 6));
		dMaxTFDR = DMGETAT(pRegionSort, ni, 7);
		dMaxTLfdr = DMGETAT(pRegionSort, ni, 8);
		dMaxFC = -DMGETAT(pRegionSort, ni, 9);
		nMaxFCPos = (int)(DMGETAT(pRegionSort, ni, 10));
		dSumT = -DMGETAT(pRegionSort, ni, 11);
		nSumProbe = -(int)(DMGETAT(pRegionSort, ni, 12));
		dSumTFDR = DMGETAT(pRegionSort, ni, 13);
		dSumTLfdr = DMGETAT(pRegionSort, ni, 14);
		nArrayId = (int)(DMGETAT(pRegionSort, ni, 15));

		fprintf(fpOut, "%d\t%s\t%d\t%d\t%c\t%d\t%f\t%d\t%f\t%f\t%f\t%d\t%f\t%d\t%f\t%f\t%d\t%s\n",
				(ni+1), strChr, nStart, nEnd, chStrand, nLen, dMaxT, nMaxTPos,
				dMaxTFDR, dMaxTLfdr, dMaxFC, nMaxFCPos, dSumT, nSumProbe, 
				dSumTFDR, dSumTLfdr, nArrayId, strGroupName);

	}

	fclose(fpOut);
	
	/* release memory */
	DestroyDoubleMatrix(pRegion);
	DestroyDoubleMatrix(pRegionSort);
	DestroyLongMatrix(pRegionSid);
	DestroyIntMatrix(pType);
	DestroyIntMatrix(pPriority);

	for(ni=0; ni<nRegionNum; ni++)
	{
		DeleteString(vReg[ni]);
		vReg[ni] = NULL;
		DeleteString(vGroup[ni]);
		vGroup[ni] = NULL;
	}
	free(vReg);
	free(vGroup);

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_RegionDetection_HMM_Main()                                   */
/*  TileMapv2 region detection, HMM.                                       */
/*  return number of regions detected.                                     */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_RegionDetection_HMM_Main(struct tagTileMapv2Param *pParam, int nLibId)
{
	/* define */
	int nRegionNum = 0;
	char strProbeFile[MED_LINE_LENGTH];
	char strMaskFile[MED_LINE_LENGTH];
	char strFCFile[MED_LINE_LENGTH];
	char strFileName[MED_LINE_LENGTH];
	struct INTMATRIX *pCol = NULL;
	struct INTMATRIX *pPermGroupLabel = NULL;
	struct INTMATRIX *pOriGroupLabel = NULL;
	struct tagBARData *pData = NULL;
	struct DOUBLEMATRIX *pRegion = NULL;
	struct DOUBLEMATRIX *pNewRegion = NULL;
	struct DOUBLEMATRIX *pRegionMaxSort = NULL;
	struct LONGMATRIX *pRegionMaxSid = NULL;
	struct DOUBLEMATRIX *pRegionMaxFDR = NULL;
	struct DOUBLEMATRIX *pRegionSumSort = NULL;
	struct LONGMATRIX *pRegionSumSid = NULL;
	struct DOUBLEMATRIX *pRegionSumFDR = NULL;
	struct INTMATRIX *pType = NULL;
	struct INTMATRIX *pPriority = NULL;

	struct DOUBLEMATRIX *pRegionCT = NULL;
	struct DOUBLEMATRIX *pRegionCTMaxSort = NULL;
	struct LONGMATRIX *pRegionCTMaxSid = NULL;
	struct DOUBLEMATRIX *pRegionCTSumSort = NULL;
	struct LONGMATRIX *pRegionCTSumSid = NULL;
	FILE *fpOut;
	char strTemp[MED_LINE_LENGTH];

	/* init */
	if(pParam == NULL)
	{
		return 0;
	}

	/* ------------- */
	/* load raw data */
	/* ------------- */
	strcpy(strMaskFile, "");
	strcpy(strFCFile, "");
	sprintf(strProbeFile, "%s%s_%d.pb.bar", pParam->strWorkPath, pParam->strProjectTitle, nLibId+1);
	if(pParam->nNoiseMask == 1)
	{
		sprintf(strMaskFile, "%s%s_%d.pb.bar.mask", pParam->strWorkPath, pParam->strProjectTitle, nLibId+1);
	}
	if(pParam->nPatternType <= 2)
	{
		sprintf(strFCFile, "%s%s_%d.fc.bar", pParam->strWorkPath, pParam->strProjectTitle, nLibId+1);
	}
	
	pData = NULL;
	pData = TileMapv2_RegionDetection_HMM_PrepareData(strProbeFile, strMaskFile, strFCFile);
	if(pData == NULL)
	{
		printf("Error: TileMapv2_RegionDetection_HMM_Main, cannot load raw data!\n");
		exit(EXIT_FAILURE);
	}

	/* -------------------- */
	/* Data Transformation  */
	/* -------------------- */
	TileMapv2_DataTransform(pData, 1, "oneminusinvlogit");
	/* if(pParam->nPatternType <= 2)
	{
		TileMapv2_DataTransform(pData, 1, "oneminusinvlogit");
	}
	else
	{
		TileMapv2_DataTransform(pData, 1, "oneminus");
	} */

	/* -------------------- */
	/* UMS                  */
	/* -------------------- */
	if(pParam->nHMMParamUserSpecified == 0)
	{
		printf("UMS ...\n");
		sprintf(pParam->strTransitionPath, "%s%s.transp", 
			pParam->strWorkPath, pParam->strProjectTitle);
		sprintf(pParam->strEmissionPath, "%s%s.emissp", 
			pParam->strWorkPath, pParam->strProjectTitle);
		TileMapv2_UMS_HMM_Main(pData, pParam, nLibId);
	}

	/* -------------------- */
	/* HMM                  */
	/* -------------------- */
	printf("HMM Decoding ...\n");
	TileMapv2_HMM_Decoding(pData, pParam, nLibId);

	/* Save HMM statistics */
	pCol = NULL;
	pCol = CreateIntMatrix(1,pData->nColNum);
	if(pCol == NULL)
	{
		printf("Error: TileMapv2_RegionDetection_HMM_Main, cannot create output column information!\n");
		exit(EXIT_FAILURE);
	}
	printf("Exporting HMM statistics...\n");
	pCol->pMatElement[0] = 1;
	pCol->pMatElement[2] = 1;
	sprintf(strFileName, "%s%s_%d.hmm.bar", pParam->strWorkPath, pParam->strProjectTitle, nLibId+1);
	Affy_SaveBAR_Columns_Fast(strFileName, pData, pCol);

	DestroyIntMatrix(pCol);

	/* Call binding regions */
	pRegion = NULL;
	pRegionMaxFDR = NULL;
	pRegionSumFDR = NULL;
	pRegion = TileMapv2_HMM_CallRegion(pData, pParam, nLibId);
	
	/* ---------------------------- */
	/* permutations to estimate FDR */
	/* ---------------------------- */
	if(pRegion != NULL)
	{
		pRegionMaxFDR = CreateDoubleMatrix(1, pRegion->nHeight);
		pRegionSumFDR = CreateDoubleMatrix(1, pRegion->nHeight);
		if( (pRegionMaxFDR == NULL) || (pRegionSumFDR == NULL) )
		{
			printf("Error: TileMapv2_RegionDetection_HMM_Main, cannot allocate memory for FDR computation!\n");
			exit(EXIT_FAILURE);
		}

		pType = NULL;
		pType = CreateIntMatrix(1, pRegion->nWidth);
		if(pType == NULL)
		{
			printf("Error: TileMapv2_RegionDetection_HMM_Main, cannot create column data type!\n");
			exit(EXIT_FAILURE);
		}
		
		/* "%d\t%d\t%d\t% 9.7e\t%d\t% 9.7e\t%d\t% 9.7e\t%d\n", 
						nSeqId, nStart, nEnd, dMaxT, nMaxTPos,
						dMaxFC, nMaxFCPos, dSumT, nProbeCount); */

		pType->pMatElement[0] = 2;
		pType->pMatElement[1] = 2;
		pType->pMatElement[2] = 2;
		pType->pMatElement[3] = 1;
		pType->pMatElement[4] = 2;
		pType->pMatElement[5] = 1;
		pType->pMatElement[6] = 2;
		pType->pMatElement[7] = 1;
		pType->pMatElement[8] = 2;

		pPriority = NULL;
		pPriority = CreateIntMatrix(1, 1);
		if(pPriority == NULL)
		{
			printf("Error: TileMapv2_RegionDetection_HMM_Main, cannot create sorting key!\n");
			exit(EXIT_FAILURE);
		}

		pPriority->pMatElement[0] = 3;
		pRegionMaxSort = NULL;
		pRegionMaxSid = NULL;
		DMSORTROWS(pRegion, pType, pPriority, &pRegionMaxSort, &pRegionMaxSid);
		
		pPriority->pMatElement[0] = 7;
		pRegionSumSort = NULL;
		pRegionSumSid = NULL;
		DMSORTROWS(pRegion, pType, pPriority, &pRegionSumSort, &pRegionSumSid);
		
		DestroyIntMatrix(pPriority);
	
		if( (pParam->nFDRType == 1) && (pParam->nPermutationNum > 0) )
		{
			/* permutations */
			TileMapv2_HMM_RegionFDR_Permutation(pType, 
				pRegionMaxSort, pRegionMaxFDR,
				pRegionSumSort, pRegionSumFDR, 
				pParam, nLibId);

			pNewRegion = NULL;
			pNewRegion = TileMapv2_RegionFDR_Reorganize(pRegion, pRegionMaxFDR, pRegionMaxSid, 
				pRegionSumFDR, pRegionSumSid);
		}
		else if(pParam->nFDRType == 0)
		{
			/* convert test statistics */
			TileMapv2_HMM_RegionFDR_LeftTail(pType,
				pRegionMaxSort, pRegionMaxFDR,
				pRegionSumSort, pRegionSumFDR, 
				pData, pParam, nLibId);

			pNewRegion = NULL;
			pNewRegion = TileMapv2_RegionFDR_Reorganize(pRegion, pRegionMaxFDR, pRegionMaxSid, 
				pRegionSumFDR, pRegionSumSid); 
		}
		else
		{
			pNewRegion = NULL;
			pNewRegion = TileMapv2_RegionFDR_Reorganize(pRegion, pRegionMaxFDR, pRegionMaxSid, 
				pRegionSumFDR, pRegionSumSid);
		}

		/* ---------------------------- */
		/* exporting results            */
		/* ---------------------------- */
		/* seqid, chr, start, end, strand, length, maxP, maxP Pos, maxP FDR, maxFC, maxFC Pos,
		sumP, sumP good probe count, sumP FDR, libid */
		nRegionNum = pNewRegion->nHeight;
		TileMapv2_MA_ExportRegion(pNewRegion, pData, pParam, nLibId);
	
		/* -------------- */
		/* release memory */
		/* -------------- */
		DestroyDoubleMatrix(pRegion);
		DestroyIntMatrix(pType);
		DestroyDoubleMatrix(pRegionMaxFDR);
		DestroyDoubleMatrix(pRegionSumFDR);
		DestroyDoubleMatrix(pRegionMaxSort);
		DestroyDoubleMatrix(pRegionSumSort);
		DestroyLongMatrix(pRegionMaxSid);
		DestroyLongMatrix(pRegionSumSid);
		DestroyDoubleMatrix(pNewRegion);
	}
	else
	{
		sprintf(strTemp, "%s%s_%d.cod", pParam->strWorkPath, pParam->strProjectTitle, (nLibId+1));
		fpOut = NULL;
		fpOut = fopen(strTemp, "w");
		if(fpOut == NULL)
		{
			printf("Error: TileMapv2_RegionDetection_HMM_Main, cannot open the output file!\n");
			exit(EXIT_FAILURE);
		}

		fclose(fpOut);
	}

	/* -------------- */
	/* release memory */
	/* -------------- */
	Affy_BARData_Destroy(&pData);
	
	/* return */
	return nRegionNum;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_RegionDetection_HMM_PrepareData()                            */
/*  TileMapv2 HMM preparation: loading probe level statistics.             */
/* ----------------------------------------------------------------------- */ 
struct tagBARData *TileMapv2_RegionDetection_HMM_PrepareData(char strProbeFile[],
							char strMaskFile[], char strFCFile[])
{
	struct tagBARData *pData = NULL;
	struct tagBARData *pMask = NULL;
	struct INTMATRIX *pFieldType = NULL;
	struct DOUBLEMATRIX **vDataVec = NULL;
	int nTotalProbeNum = 0;
	int nj;

	
	/* -------------------- */
	/* load probe level t   */
	/* -------------------- */
    pData = NULL;
	pData = Affy_LoadBAR_Fast(strProbeFile);
	if(pData == NULL)
	{
		printf("Error: TileMapv2_RegionDetection_HMM_PrepareData, cannot load raw data!\n");
		exit(EXIT_FAILURE);
	}

	/* 1 position + 1 probe level statistics + 1 HMM + 1 mask + 1 Fold change */ 
	pData->nColNum = 2+3;
	pFieldType = pData->pFieldType;
	pData->pFieldType = NULL;
	pData->pFieldType = CreateIntMatrix(1,pData->nColNum);
	if(pData->pFieldType == NULL)
	{
		printf("Error: TileMapv2_RegionDetection_HMM_PrepareData, cannot create memory for field types!\n");
		exit(EXIT_FAILURE);
	}
	pData->pFieldType->pMatElement[0] = pFieldType->pMatElement[0];
	pData->pFieldType->pMatElement[1] = pFieldType->pMatElement[1];
	DestroyIntMatrix(pFieldType);
	pData->pFieldType->pMatElement[2] = 1;
	
	for(nj=0; nj<pData->nSeqNum; nj++)
	{
		nTotalProbeNum += pData->vSeqData[nj]->nDataNum;
		pData->vSeqData[nj]->nColNum = pData->nColNum;
		vDataVec = pData->vSeqData[nj]->vData;
		pData->vSeqData[nj]->vData = NULL;
		pData->vSeqData[nj]->vData = (struct DOUBLEMATRIX **)calloc(pData->vSeqData[nj]->nColNum, sizeof(struct DOUBLEMATRIX *));
		if(pData->vSeqData[nj]->vData == NULL)
		{
			printf("Error: TileMapv2_RegionDetection_HMM_PrepareData, cannot create memory for tracking intensity data!\n");
			exit(EXIT_FAILURE);
		}
		pData->vSeqData[nj]->vData[0] = vDataVec[0];
		pData->vSeqData[nj]->vData[1] = vDataVec[1];
		free(vDataVec);

		if(pData->vSeqData[nj]->nDataNum <= 0)
			continue;

		pData->vSeqData[nj]->vData[2] = CreateDoubleMatrix(1, pData->vSeqData[nj]->nDataNum);
		if(pData->vSeqData[nj]->vData[2] == NULL)
		{
			printf("Error: TileMapv2_RegionDetection_HMM_PrepareData, cannot create memory for storing MA statistics!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* -------------------- */
	/* load masks if needed */
	/* -------------------- */
	if(strcmp(strMaskFile, "") != 0)
	{
		pMask = NULL;
		pMask = Affy_LoadBAR_Fast(strMaskFile);
		if(pMask == NULL)
		{
			printf("Error: TileMapv2_RegionDetection_HMM_PrepareData, cannot load mask data!\n");
			exit(EXIT_FAILURE);
		}

		if(pMask->nSeqNum != pData->nSeqNum)
		{
			printf("Error: TileMapv2_RegionDetection_HMM_PrepareData, array types do not match!\n");
			exit(EXIT_FAILURE);
		}

		for(nj=0; nj<pMask->nSeqNum; nj++)
		{
			if(pMask->vSeqData[nj]->nDataNum != pData->vSeqData[nj]->nDataNum)
			{
				printf("Error: TileMapv2_RegionDetection_HMM_PrepareData, array types do not match!\n");
				exit(EXIT_FAILURE);
			}
			pData->vSeqData[nj]->vData[3] = pMask->vSeqData[nj]->vData[1];
			pMask->vSeqData[nj]->vData[1] = NULL;
		}
        
		pData->pFieldType->pMatElement[3] = pMask->pFieldType->pMatElement[1];

		Affy_BARData_Destroy(&pMask);
	}


	/* -------------------- */
	/* load FC if needed    */
	/* -------------------- */
	if(strcmp(strFCFile, "") != 0)
	{
		pMask = NULL;
		pMask = Affy_LoadBAR_Fast(strFCFile);
		if(pMask == NULL)
		{
			/* return */
			return pData;
		}

		if(pMask->nSeqNum != pData->nSeqNum)
		{
			printf("Error: TileMapv2_RegionDetection_HMM_PrepareData, array types do not match!\n");
			exit(EXIT_FAILURE);
		}

		for(nj=0; nj<pMask->nSeqNum; nj++)
		{
			if(pMask->vSeqData[nj]->nDataNum != pData->vSeqData[nj]->nDataNum)
			{
				printf("Error: TileMapv2_RegionDetection_HMM_PrepareData, array types do not match!\n");
				exit(EXIT_FAILURE);
			}
			pData->vSeqData[nj]->vData[4] = pMask->vSeqData[nj]->vData[1];
			pMask->vSeqData[nj]->vData[1] = NULL;
		}
        
		pData->pFieldType->pMatElement[4] = pMask->pFieldType->pMatElement[1];

		Affy_BARData_Destroy(&pMask);
	}

	/* return */
	return pData;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_HMM()                                                        */
/*  TileMapv2 HMM.                                                         */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_HMM(struct tagBARData *pData, struct tagTileMapv2Param *pParam, int nLibId)
{
	/* -------------------- */
	/* UMS                  */
	/* -------------------- */
	if(pParam->nHMMParamUserSpecified == 0)
	{
		printf("UMS ...\n");
		TileMapv2_UMS_HMM_Main(pData, pParam, nLibId);
	}

	/* -------------------- */
	/* HMM                  */
	/* -------------------- */
	printf("HMM Decoding ...\n");
	TileMapv2_HMM_Decoding(pData, pParam, nLibId);
	
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_UMS_HMM_Main()                                               */
/*  TileMapv2 UMS HMM.                                                     */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_UMS_HMM_Main(struct tagBARData *pData, struct tagTileMapv2Param *pParam, int nLibId)
{
	/* define */
	FILE *fpOut;
	struct DOUBLEMATRIX *pSelect;
	struct DOUBLEMATRIX *pScore;
	struct DOUBLEMATRIX *pTemp;
	struct DOUBLEMATRIX *pIntx;
	struct DOUBLEMATRIX *pF0;
	struct DOUBLEMATRIX *pF1;
	struct DOUBLEMATRIX *pTrans;
	double dTheta = 0.0;
	double a0,a1;

	char strFileName[MED_LINE_LENGTH];
	int nProbeNum,nTotalProbeNum;
	int ni,nj;
	double *vT,*vMask;

	/* init */

	/* get probe number */
	nProbeNum = 0;
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		nProbeNum += pData->vSeqData[ni]->nDataNum;
	}

	if(nProbeNum <= 0)
	{
		printf("Error: TileMapv2_UMS_HMM_Main, empty HMM input!\n");
		exit(EXIT_FAILURE);
	}
	
	pScore = NULL;
	pScore = CreateDoubleMatrix(1, nProbeNum);
	if(pScore == NULL)
	{
		printf("Error: TileMapv2_UMS_HMM_Main, cannot create score matrix!\n");
		exit(EXIT_FAILURE);
	}

	/* load data */
	if(pParam->nNoiseMask == 1)
	{
		nTotalProbeNum = 0;
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			if(pData->vSeqData[ni]->nDataNum <= 0)
				continue;
			
			vT = pData->vSeqData[ni]->vData[1]->pMatElement;
			vMask = pData->vSeqData[ni]->vData[3]->pMatElement;
			for(nj=0; nj<pData->vSeqData[ni]->nDataNum; nj++)
			{
				if(vMask[nj] < 0.5)
				{
					pScore->pMatElement[nTotalProbeNum] = vT[nj];
					nTotalProbeNum++;
				}
			}
		}

		pTemp = NULL;
		pTemp = CreateDoubleMatrix(1, nTotalProbeNum);
		if(pTemp == NULL)
		{
			printf("Error: TileMapv2_UMS_HMM_Main, cannot create selection matrix!\n");
			exit(EXIT_FAILURE);
		}

		memcpy(pTemp->pMatElement, pScore->pMatElement, sizeof(double)*nTotalProbeNum);
		DestroyDoubleMatrix(pScore);
		pScore = pTemp;
	}
	else
	{
		nTotalProbeNum = 0;
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			if(pData->vSeqData[ni]->nDataNum <= 0)
				continue;
			
			memcpy(pScore->pMatElement+nTotalProbeNum, pData->vSeqData[ni]->vData[1]->pMatElement, sizeof(double)*pData->vSeqData[ni]->nDataNum);
			nTotalProbeNum += pData->vSeqData[ni]->nDataNum;
		}

		if(nTotalProbeNum != nProbeNum)
		{
			printf("Error: TileMapv2_UMS_HMM_Main, probe number not match!\n");
			exit(EXIT_FAILURE);
		}
	}

	pSelect = NULL;
	pSelect = DMCLONE(pScore);
	if(pSelect == NULL)
	{
		printf("Error: TileMapv2_UMS_HMM_Main, cannot create selection matrix!\n");
		exit(EXIT_FAILURE);
	}

	/* UMS */
	pIntx = NULL;
	pIntx = CreateDoubleMatrix(1, pParam->nGridSize);
	if(pIntx == NULL)
	{
		printf("Error: TileMapv2_UMS_HMM_Main, cannot create interval matrix!\n");
		exit(EXIT_FAILURE);
	}

	pF0 = NULL;
	pF0 = CreateDoubleMatrix(1, pParam->nGridSize);
	if(pF0 == NULL)
	{
		printf("Error: TileMapv2_UMS_HMM_Main, cannot create f0 matrix!\n");
		exit(EXIT_FAILURE);
	}

	pF1 = NULL;
	pF1 = CreateDoubleMatrix(1, pParam->nGridSize);
	if(pF1 == NULL)
	{
		printf("Error: TileMapv2_UMS_HMM_Main, cannot create f1 matrix!\n");
		exit(EXIT_FAILURE);
	}

	TileMap_UMS(pSelect, pParam->dTp, pParam->dTq, pParam->nOffset, pScore, pParam->nGridSize, pIntx, pF0, pF1, &dTheta);

	/* estimate emission */
	pTrans = NULL;
	pTrans = CreateDoubleMatrix(2, 2);
	if(pTrans == NULL)
	{
		printf("Error: TileMapv2_UMS_HMM_Main, cannot create transition matrix!\n");
		exit(EXIT_FAILURE);
	}

	a1 = 1.0/(double)pParam->nExpLen;
	if((dTheta < ZERO_BOUND) || (dTheta > 1.0-ZERO_BOUND))
	{
		printf("Error: TileMapv2_UMS_HMM_Main, theta=0, cannot create HMM parameters!\n");
		exit(EXIT_FAILURE);
	}
	a0 = a1*dTheta/(1.0-dTheta);

	DMSETAT(pTrans, 0, 0, (1.0-a0));
	DMSETAT(pTrans, 0, 1, a0);
	DMSETAT(pTrans, 1, 0, a1);
	DMSETAT(pTrans, 1, 1, (1.0-a1));

	/* write */
	sprintf(strFileName, "%s.%d.txt", pParam->strTransitionPath, (nLibId+1));
	DMSAVE(pTrans, strFileName);

	sprintf(strFileName, "%s.%d.txt", pParam->strEmissionPath, (nLibId+1));
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: TileMapv2_UMS_HMM_Main, cannot export emission probability!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<pParam->nGridSize; ni++)
	{
		fprintf(fpOut, "% 9.7e\t", pIntx->pMatElement[ni]);
	}
	fprintf(fpOut, "\n");
	for(ni=0; ni<pParam->nGridSize; ni++)
	{
		fprintf(fpOut, "% 9.7e\t", pF0->pMatElement[ni]);
	}
	fprintf(fpOut, "\n");
	for(ni=0; ni<pParam->nGridSize; ni++)
	{
		fprintf(fpOut, "% 9.7e\t", pF1->pMatElement[ni]);
	}
	fprintf(fpOut, "\n");

	fclose(fpOut);


	/* destroy memory */
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
/*  TileMapv2_HMM_Decoding()                                               */
/*  TileMapv2 HMM decoding.                                                */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_HMM_Decoding(struct tagBARData *pData, struct tagTileMapv2Param *pParam, 
						   int nLibId)
{
	/* define */
	char strFileName[MED_LINE_LENGTH];
	struct DOUBLEMATRIX *pTransition;
	struct DOUBLEMATRIX *pEmission;
	struct DOUBLEMATRIX *pStationary;
	int nStateNum;
	double *pEle;
	int ni,nj;
	double dTemp;

	/* init */
	if(pParam->nHMMParamUserSpecified == 1)
	{
		sprintf(strFileName, "%s", pParam->strTransitionPath);
	}
	else
	{
		sprintf(strFileName, "%s.%d.txt", pParam->strTransitionPath, (nLibId+1));
	}
	pTransition = NULL;
	pTransition = DMLOAD(strFileName);
	if(pTransition == NULL)
	{
		printf("Error: TileMapv2_HMM_Decoding, cannot load transition probability!\n");
		exit(EXIT_FAILURE);
	}

	if(pParam->nHMMParamUserSpecified == 1)
	{
		sprintf(strFileName, "%s", pParam->strEmissionPath);
	}
	else
	{
		sprintf(strFileName, "%s.%d.txt", pParam->strEmissionPath, (nLibId+1));
	}
	pEmission = NULL;
	pEmission = DMLOAD(strFileName);
	if(pEmission == NULL)
	{
		printf("Error: TileMapv2_HMM_Decoding, cannot load emission probability!\n");
		exit(EXIT_FAILURE);
	}

	if((pTransition->nHeight+1) != pEmission->nHeight)
	{
		printf("Error: TileMapv2_HMM_Decoding, state number not match!\n");
		exit(EXIT_FAILURE);
	}

	nStateNum = pTransition->nHeight;
	if(nStateNum != 2)
	{
		printf("Error: TileMapv2_HMM_Decoding, this function only support 2 states now!\n");
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

	/* calculate posterior probability */
	for(ni=0; ni<pData->nSeqNum; ni++)
	{
		if(pData->vSeqData[ni]->nDataNum <= 0)
			continue;

		TileMapv2_HMM_Decoding_Chr(pData->vSeqData[ni]->nDataNum,
			nStateNum, pStationary, pTransition, pEmission, pParam->nGap, 
			pData->vSeqData[ni]->vData[0], pData->vSeqData[ni]->vData[1],
			pData->vSeqData[ni]->vData[2], pData->vSeqData[ni]->vData[3]);
	}

	/* release memory */
	DestroyDoubleMatrix(pTransition);
	DestroyDoubleMatrix(pEmission);
	DestroyDoubleMatrix(pStationary);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_HMM_Decoding_Chr()                                           */
/*  TileMapv2 HMM decoding.                                                */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_HMM_Decoding_Chr(int nProbeNum, int nStateNum, 
			struct DOUBLEMATRIX *pStationary, struct DOUBLEMATRIX *pTransition, 
			struct DOUBLEMATRIX *pEmission, double dGapDist, 
			struct DOUBLEMATRIX *pPosition, struct DOUBLEMATRIX *pScore,
			struct DOUBLEMATRIX *pPosterior, struct DOUBLEMATRIX *pMask)
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
	int nSkipF,nSkipB;
	int nLastId;
	int nEqualLenInterval;
	double dIntS,dIntE,dIntStep;

	/* check */
	if( (nProbeNum <= 0) || (nStateNum <= 0) )
	{
		printf("Error: TileMapv2_HMM_Decoding_Chr, nProbeNum/nStateNum <=0!\n");
		exit(EXIT_FAILURE);
	}
	if((pScore == NULL) || (pStationary == NULL) || (pTransition == NULL) 
		|| (pEmission == NULL) || (pPosition == NULL) || (pPosterior == NULL))
	{
		printf("Error: TileMapv2_HMM_Decoding_Chr, no input data/parameters!\n");
		exit(EXIT_FAILURE);
	}
	if(pScore->nWidth != nProbeNum)
	{
		printf("Error: TileMapv2_HMM_Decoding_Chr, probe number not match!\n");
		exit(EXIT_FAILURE);
	}
	if(pPosterior->nWidth != nProbeNum)
	{
		printf("Error: TileMapv2_HMM_Decoding_Chr, probe number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* init */
	dInitMax = -DM_ACCESS_VIOLATION;
	
	/* by default, use intervals of equal length for likelihood calculation */
	nEqualLenInterval = 1;
	if(pEmission->nWidth <= 2)
	{
		printf("Error: TileMapv2_HMM_Decoding_Chr, emission probability need to be set for more than two intervals!\n");
		exit(EXIT_FAILURE);
	}

	dIntS = DMGETAT(pEmission, 0, 0);
	dIntE = DMGETAT(pEmission, 0, (pEmission->nWidth-2));
	
	dIntStep = (dIntE-dIntS)/(double)(pEmission->nWidth-2);
	if(dIntStep < 0.0)
	{
		printf("Error: TileMapv2_HMM_Decoding_Chr, negative interval step length!\n");
		exit(EXIT_FAILURE);
	}

	/* prepare space */
	vForwardSum = NULL;
	vForwardSum = (struct DOUBLEMATRIX **)calloc(nStateNum, sizeof(struct DOUBLEMATRIX *));
	if(vForwardSum == NULL)
	{
		printf("Error: TileMapv2_HMM_Decoding_Chr, cannot create space for HMM calculation!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nStateNum; ni++)
	{
		vForwardSum[ni] = CreateDoubleMatrix(1, nProbeNum);
		if(vForwardSum[ni] == NULL)
		{
			printf("Error: TileMapv2_HMM_Decoding_Chr, cannot create space for HMM calculation!\n");
			exit(EXIT_FAILURE);
		}
	}

	vBackwardSum = NULL;
	vBackwardSum = (struct DOUBLEMATRIX **)calloc(nStateNum, sizeof(struct DOUBLEMATRIX *));
	if(vBackwardSum == NULL)
	{
		printf("Error: TileMapv2_HMM_Decoding_Chr, cannot create space for HMM calculation!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nStateNum; ni++)
	{
		vBackwardSum[ni] = CreateDoubleMatrix(1, nProbeNum);
		if(vBackwardSum[ni] == NULL)
		{
			printf("Error: TileMapv2_HMM_Decoding_Chr, cannot create space for HMM calculation!\n");
			exit(EXIT_FAILURE);
		}
	}

	pDP = NULL;
	pDP = CreateDoubleMatrix(nStateNum, 1);
	if(pDP == NULL)
	{
		printf("Error: TileMapv2_HMM_Decoding_Chr, cannot create space for HMM calculation!\n");
		exit(EXIT_FAILURE);
	}

	/* forward summation */
	nLastId = -1;

	for(nj=0; nj<nProbeNum; nj++)
	{
		if(pMask != NULL)
		{
			if(pMask->pMatElement[nj] > 0.5)
			{
				continue;
			}
		}

		/* if the first position */
		if(nLastId < 0)
		{
			for(ni=0; ni<nStateNum; ni++)
			{
				dTemp = pStationary->pMatElement[ni]+Tiling_BindingRegionSelection_HMM_GetEmission(pEmission, ni, pScore->pMatElement[nj], nEqualLenInterval, dIntS, dIntE, dIntStep);
				vForwardSum[ni]->pMatElement[nj] = dTemp;
			}
		}
		/* if not the first position */
		else
		{
			dDist = pPosition->pMatElement[nj] - pPosition->pMatElement[nLastId];

			for(ni=0; ni<nStateNum; ni++)
			{
				dMax = dInitMax;
				for(nk=0; nk<nStateNum; nk++)
				{
					dTemp = vForwardSum[nk]->pMatElement[nLastId] + Tiling_BindingRegionSelection_HMM_GetTransition(pStationary, pTransition, nk, ni, dDist, dGapDist); 
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

		/* update last id */
		nLastId = nj;
	}

	if(nLastId >= 0)
	{
		dMax = dInitMax;
		for(ni=0; ni<nStateNum; ni++)
		{
			dTemp = vForwardSum[ni]->pMatElement[nLastId]; 
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
		nSkipF = 0;
	}
	else
	{
		dSumF = 0.0;
		nSkipF = 1;
	}

	/* backward summation */
	nLastId = nProbeNum;
	for(nj=nProbeNum-1; nj>=0; nj--)
	{
		if(pMask != NULL)
		{
			if(pMask->pMatElement[nj] > 0.5)
			{
				continue;
			}
		}

		/* if the last position */
		if(nLastId >= nProbeNum)
		{
			for(ni=0; ni<nStateNum; ni++)
			{
				vBackwardSum[ni]->pMatElement[nj] = 0.0;
			}
		}
		/* if not the last position */
		else
		{
			dDist = pPosition->pMatElement[nLastId]-pPosition->pMatElement[nj];
		
			for(ni=0; ni<nStateNum; ni++)
			{
				dMax = dInitMax;
				for(nk=0; nk<nStateNum; nk++)
				{
					dTemp = vBackwardSum[nk]->pMatElement[nLastId] 
						+ Tiling_BindingRegionSelection_HMM_GetTransition(pStationary, pTransition, ni, nk, dDist, dGapDist) 
						+ Tiling_BindingRegionSelection_HMM_GetEmission(pEmission, nk, pScore->pMatElement[nLastId], nEqualLenInterval, dIntS, dIntE, dIntStep); 
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

		/* update last id */
		nLastId = nj;
	}

	if(nLastId < nProbeNum)
	{
		dMax = dInitMax;
		for(ni=0; ni<nStateNum; ni++)
		{
			dTemp = pStationary->pMatElement[ni] 
				+ Tiling_BindingRegionSelection_HMM_GetEmission(pEmission, ni, pScore->pMatElement[nLastId], nEqualLenInterval, dIntS, dIntE, dIntStep)
				+ vBackwardSum[ni]->pMatElement[nLastId]; 
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
		nSkipB = 0;
	}
	else
	{
		dSumB = 0.0;
		nSkipB = 1;
	}

	if( (nSkipF == 1) || (nSkipB == 1))
	{
		/* printf("Warning: TileMapv2_HMM_Decoding_Chr, empty chain!\n"); */
		return PROC_SUCCESS;
	}

	if(fabs(dSumF-dSumB) > 1e-3)
	{
		printf("Error: TileMapv2_HMM_Decoding_Chr, forward and backbward summation not match!\n");
		exit(EXIT_FAILURE);
	}
	dSumF = (dSumF+dSumB)/2.0;

	/* posterior calculation */
	for(nj=0; nj<nProbeNum; nj++)
	{
		if(pMask != NULL)
		{
			if(pMask->pMatElement[nj] > 0.5)
			{
				pPosterior->pMatElement[nj] = 0.0;
				continue;
			}
		}

		dSum = 0.0;
		for(ni=0; ni<nStateNum; ni++)
		{
			dTemp = vForwardSum[ni]->pMatElement[nj]+vBackwardSum[ni]->pMatElement[nj]-dSumF;
			dTemp = exp(dTemp);
			dSum += dTemp;
			if(ni == 1)
				pPosterior->pMatElement[nj] = dTemp;
		}
		pPosterior->pMatElement[nj] /= dSum;
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
/*  TileMapv2_HMM_CallRegion()                                             */
/*  TileMapv2 call binding regions based on HMM statistics.                */
/* ----------------------------------------------------------------------- */
struct DOUBLEMATRIX *TileMapv2_HMM_CallRegion(struct tagBARData *pData, 
					struct tagTileMapv2Param *pParam, int nLibId)
{
	/* define */
	struct DOUBLEMATRIX *pRegion = NULL;
	int nRegionCount;
	double dDist;
	double dSumT;
	double dMaxT,dMaxFC;
	int nMaxTPos,nMaxFCPos;
	int nProbeCount;
	int nStart,nEnd;
	int nLastRegionProbe;
	int nSeqId;
	double *vM,*vMask,*vFC,*vPos;
	int nGoodProbe;
	int nInterGoodProbe;
	int ni;
	
	/* define */
	FILE *fpOut;
	char strFileName[MED_LINE_LENGTH];
	
	/* init */
	sprintf(strFileName, "%s%s_%d.hmmreg.tmpout", pParam->strWorkPath, 
		pParam->strProjectTitle, nLibId+1);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: TileMapv2_HMM_CallRegion, cannot open *.tmpout file for output!\n");
		exit(EXIT_FAILURE);
	}
	
	nRegionCount = 0;
	for(nSeqId=0; nSeqId<pData->nSeqNum; nSeqId++)
	{
		if(pData->vSeqData[nSeqId]->nDataNum <= 0)
			continue;

		dDist = 0.0;
		nProbeCount = 0;
		dSumT = 0.0;
		dMaxT = -1e20;
		dMaxFC = -1e20;
		nStart = -1;
		nEnd = -1;
		nMaxTPos = -1;
		nMaxFCPos = -1;
		nLastRegionProbe = -1;
		nInterGoodProbe = 0;

		vPos = pData->vSeqData[nSeqId]->vData[0]->pMatElement;
		vM = pData->vSeqData[nSeqId]->vData[2]->pMatElement;
		if(pData->vSeqData[nSeqId]->vData[3] != NULL)
		{
			vMask = pData->vSeqData[nSeqId]->vData[3]->pMatElement;
		}
		else
		{
			vMask = NULL;
		}
		if(pData->vSeqData[nSeqId]->vData[4] != NULL)
		{
			vFC = pData->vSeqData[nSeqId]->vData[4]->pMatElement;
		}
		else
		{
			vFC = NULL;
		}
	
		for(ni=0; ni<pData->vSeqData[nSeqId]->nDataNum; ni++)
		{
			/* judge whether the probe is good or not */
			nGoodProbe = 1;
			if(vMask != NULL)
			{
				if(vMask[ni] > 0.5)
					nGoodProbe = 0;
			}

			if(nGoodProbe == 0)
				continue;

			/* if the probe does not pass the cutoff value */
			if(vM[ni] < pParam->dPostCut)
			{
				nInterGoodProbe += 1;
			}
			/* if the probe passes the cutoff value */
			else
			{
				/* if not the first probe, does it belong to the previous region? */
				if(nLastRegionProbe >= 0)
				{
					dDist = vPos[ni]-vPos[nLastRegionProbe];
				}
				else
				{
					dDist = 1e20;
				}

				/* if it belongs to the previous region, update the previous region */
				if( (dDist <= pParam->nGap) && (nInterGoodProbe <= pParam->nGapW) )
				{
					nLastRegionProbe = ni;
					nEnd = (int)(vPos[ni]);
					nProbeCount += 1;
					dSumT += vM[ni];
					if(dMaxT < vM[ni])
					{
						dMaxT = vM[ni];
						nMaxTPos = (int)(vPos[ni]);
					}
					
					if(vFC == NULL)
					{
						dMaxFC = dMaxT;
						nMaxFCPos = nMaxTPos;
					}
					else
					{
						if(dMaxFC < vFC[ni])
						{
							dMaxFC = vFC[ni];
							nMaxFCPos = (int)(vPos[ni]);
						}
					}
				}

				/* if it is the first probe, or 
				   if it does not belong to the previous region, save the previous region
				   and initiate a new region */
				else
				{
					if( (nProbeCount >= pParam->nMinRegProbeNum) &&
						((nEnd-nStart+1) >= pParam->nMinRegLen) )
					{
						fprintf(fpOut, "%d\t%d\t%d\t% 9.7e\t%d\t% 9.7e\t%d\t% 9.7e\t%d\n", 
							nSeqId, nStart, nEnd, dMaxT, nMaxTPos,
							dMaxFC, nMaxFCPos, dSumT, nProbeCount);
						
						nRegionCount++;
						
						nProbeCount = 0;
						dSumT = 0.0;
						dMaxT = -1e20;
						dMaxFC = -1e20;
						nStart = -1;
						nEnd = -1;
						nMaxTPos = -1;
						nMaxFCPos = -1;
						nLastRegionProbe = -1;
					}

					nStart = (int)(vPos[ni]);
					nEnd = (int)(vPos[ni]);
					nProbeCount = 1;
					dSumT = vM[ni];
					dMaxT = vM[ni];
					nMaxTPos = (int)(vPos[ni]);
					if(vFC == NULL)
					{
						dMaxFC = dMaxT;
						nMaxFCPos = nMaxTPos;
					}
					else
					{
						dMaxFC = vFC[ni];
						nMaxFCPos = (int)(vPos[ni]);
					}
					nLastRegionProbe = ni;
				}

				/* reset intergoodprobe number */
				nInterGoodProbe = 0;
			}
		}
		
		if( (nProbeCount >= pParam->nMinRegProbeNum) &&
						((nEnd-nStart+1) >= pParam->nMinRegLen) )
		{
			fprintf(fpOut, "%d\t%d\t%d\t% 9.7e\t%d\t% 9.7e\t%d\t% 9.7e\t%d\n", 
						nSeqId, nStart, nEnd, dMaxT, nMaxTPos,
						dMaxFC, nMaxFCPos, dSumT, nProbeCount);

			nRegionCount++;
			nProbeCount = 0;
			dSumT = 0.0;
			dMaxT = -1e20;
			dMaxFC = -1e20;
			nStart = -1;
			nEnd = -1;
			nMaxTPos = -1;
			nMaxFCPos = -1;
			nLastRegionProbe = -1;
		}
	}

	/* close file */
	fclose(fpOut);

	/* merge signals */
	if(nRegionCount <= 0)
	{
		/* release memory */
		RemoveFiles(strFileName);
		printf("No regions found using the specified criteria!\n");
		return NULL;
	}

	pRegion = NULL;
	pRegion = DMLOAD(strFileName);
	if(pRegion == NULL)
	{
		printf("Error: TileMapv2_HMM_CallRegion, cannot reload binding regions!\n");
		exit(EXIT_FAILURE);
	}
	if(pRegion->nHeight != nRegionCount)
	{
		printf("Error: TileMapv2_HMM_CallRegion, cannot reload binding regions correctly!\n");
		exit(EXIT_FAILURE);
	}
	RemoveFiles(strFileName);

	/* return */
	return pRegion;
} 

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_HMM_RegionFDR_Permutation()                                  */
/*  TileMap Region FDR.                                                    */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_HMM_RegionFDR_Permutation(struct INTMATRIX *pType,
				struct DOUBLEMATRIX *pRegionMaxSort, struct DOUBLEMATRIX *pRegionMaxFDR,
				struct DOUBLEMATRIX *pRegionSumSort, struct DOUBLEMATRIX *pRegionSumFDR, 
				struct tagTileMapv2Param *pParam, int nLibId)
{
	/* define */
	char strProbeFile[MED_LINE_LENGTH];
	char strMaskFile[MED_LINE_LENGTH];
	char strFCFile[MED_LINE_LENGTH];
	char strFileName[MED_LINE_LENGTH];
	struct tagBARData *pData = NULL;

	struct INTMATRIX *pPriority = NULL;
	struct DOUBLEMATRIX *pRegionCT = NULL;
	struct DOUBLEMATRIX *pRegionCTMaxSort = NULL;
	struct LONGMATRIX *pRegionCTMaxSid = NULL;
	struct DOUBLEMATRIX *pRegionCTSumSort = NULL;
	struct LONGMATRIX *pRegionCTSumSid = NULL;
	int nCol;
	int nPermId;
	struct INTMATRIX *pCol = NULL;
	

	/* permutation */
	for(nPermId=0; nPermId<pParam->nPermutationNum; nPermId++)
	{
		printf("perm %d...\n", nPermId);

		/* ------------- */
		/* load raw data */
		/* ------------- */
		strcpy(strMaskFile, "");
		strcpy(strFCFile, "");
		sprintf(strProbeFile, "%s%s_%d.perm%d.pb.bar", 
			pParam->strWorkPath, pParam->strProjectTitle, nLibId+1, nPermId);
		if(pParam->nNoiseMask == 1)
		{
			sprintf(strMaskFile, "%s%s_%d.perm%d.pb.bar.mask", 
				pParam->strWorkPath, pParam->strProjectTitle, nLibId+1, nPermId);
		}
		if(pParam->nPatternType <= 2)
		{
			sprintf(strFCFile, "%s%s_%d.perm%d.fc.bar", 
				pParam->strWorkPath, pParam->strProjectTitle, nLibId+1, nPermId);
		}
		
		pData = NULL;
		pData = TileMapv2_RegionDetection_HMM_PrepareData(strProbeFile, strMaskFile, strFCFile);
		if(pData == NULL)
		{
			printf("Error: TileMapv2_HMM_RegionFDR_Permutation, cannot load raw data!\n");
			exit(EXIT_FAILURE);
		}

		/* -------------------- */
		/* HMM                  */
		/* -------------------- */
		TileMapv2_DataTransform(pData, 1, "oneminusinvlogit");
		/* if(pParam->nPatternType <= 2)
		{
			TileMapv2_DataTransform(pData, 1, "oneminusinvlogit");
		}
		else
		{
			TileMapv2_DataTransform(pData, 1, "oneminus");
		} */

		if(pParam->nHMMParamUserSpecified == 0)
		{
			printf("UMS ...\n");
			sprintf(pParam->strTransitionPath, "%s%s.transp.perm%d", 
				pParam->strWorkPath, pParam->strProjectTitle, nPermId);
			sprintf(pParam->strEmissionPath, "%s%s.emissp.perm%d", 
				pParam->strWorkPath, pParam->strProjectTitle, nPermId);
			TileMapv2_UMS_HMM_Main(pData, pParam, nLibId);
		}

		/* -------------------- */
		/* HMM                  */
		/* -------------------- */
		printf("HMM Decoding ...\n");
		TileMapv2_HMM_Decoding(pData, pParam, nLibId);
		sprintf(strFileName, "%s*", pParam->strTransitionPath);
		RemoveFiles(strFileName);
		sprintf(strFileName, "%s*", pParam->strEmissionPath);
		RemoveFiles(strFileName);

		/* Call binding regions */
		pRegionCT = NULL;
		pRegionCT = TileMapv2_HMM_CallRegion(pData, pParam, nLibId);
	
		if(pRegionCT != NULL)
		{
			/* sort region, regionCT */
			
			/* maxT */
			pPriority = NULL;
			pPriority = CreateIntMatrix(1, 1);
			if(pPriority == NULL)
			{
				printf("Error: TileMapv2_HMM_RegionFDR_Permutation, cannot create sorting key!\n");
				exit(EXIT_FAILURE);
			}
			pPriority->pMatElement[0] = 3;
			pRegionCTMaxSort = NULL;
			pRegionCTMaxSid = NULL;
			DMSORTROWS(pRegionCT, pType, pPriority, &pRegionCTMaxSort, &pRegionCTMaxSid);

			/* sumT */
			pPriority->pMatElement[0] = 7;
			pRegionCTSumSort = NULL;
			pRegionCTSumSid = NULL;
			DMSORTROWS(pRegionCT, pType, pPriority, &pRegionCTSumSort, &pRegionCTSumSid);
			
			DestroyIntMatrix(pPriority);

			/* count */
			nCol = 3;
			TileMapv2_RegionFDR_Count(pRegionMaxSort, pRegionCTMaxSort, pRegionMaxFDR, nCol);

			nCol = 7;
			TileMapv2_RegionFDR_Count(pRegionSumSort, pRegionCTSumSort, pRegionSumFDR, nCol);

			/* -------------- */
			/* release memory */
			/* -------------- */
			DestroyDoubleMatrix(pRegionCTMaxSort);
			DestroyDoubleMatrix(pRegionCTSumSort);
			DestroyLongMatrix(pRegionCTMaxSid);
			DestroyLongMatrix(pRegionCTSumSid);
			DestroyDoubleMatrix(pRegionCT);
		}
		else
		{
			/* FDR == 0*/
		}

		/* -------------------- */
		/* Destroy              */
		/* -------------------- */
		Affy_BARData_Destroy(&pData);
		RemoveFiles(strProbeFile);
		if(strcmp(strMaskFile, "") != 0)
			RemoveFiles(strMaskFile);
		if(strcmp(strFCFile, "") != 0)
			RemoveFiles(strFCFile);
	}

	/* FDR */
	DMPDIVTS(pRegionMaxFDR, (double)pParam->nPermutationNum);
	DMPDIVTS(pRegionSumFDR, (double)pParam->nPermutationNum);
	TileMapv2_RegionFDR_Compute(pRegionMaxFDR);
	TileMapv2_RegionFDR_Compute(pRegionSumFDR);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_HMM_RegionFDR_LeftTail()                                     */
/*  TileMap Region FDR.                                                    */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_HMM_RegionFDR_LeftTail(struct INTMATRIX *pType,
				struct DOUBLEMATRIX *pRegionMaxSort, struct DOUBLEMATRIX *pRegionMaxFDR,
				struct DOUBLEMATRIX *pRegionSumSort, struct DOUBLEMATRIX *pRegionSumFDR, 
				struct tagBARData *pData, struct tagTileMapv2Param *pParam, int nLibId)
{
	/* define */
	struct INTMATRIX *pCol = NULL;
	struct INTMATRIX *pPriority = NULL;

	struct DOUBLEMATRIX *pRegionCT = NULL;
	struct DOUBLEMATRIX *pRegionCTMaxSort = NULL;
	struct LONGMATRIX *pRegionCTMaxSid = NULL;
	struct DOUBLEMATRIX *pRegionCTSumSort = NULL;
	struct LONGMATRIX *pRegionCTSumSid = NULL;
	int nCol;

	printf("FDR from left tail ...\n");
	TileMapv2_DataTransform(pData, 1, "oneminus");

	pCol = NULL;
	pCol = CreateIntMatrix(1, pData->nColNum);
	pCol->pMatElement[4] = 1;
	if(pCol == NULL)
	{
		printf("Error: TileMapv2_HMM_RegionFDR_LeftTail, cannot create column information!\n");
		exit(EXIT_FAILURE);
	}

	Affy_BARData_SubTS(pData, 0.0, pCol);
	DestroyIntMatrix(pCol);
	
	/* -------------------- */
	/* HMM                  */
	/* -------------------- */
	printf("HMM Decoding ...\n");
	TileMapv2_HMM_Decoding(pData, pParam, nLibId);

	/* Call binding regions */
	pRegionCT = NULL;
	pRegionCT = TileMapv2_HMM_CallRegion(pData, pParam, nLibId);

    if(pRegionCT != NULL)
	{
		/* sort region, regionCT */
		
		/* maxT */
		pPriority = NULL;
		pPriority = CreateIntMatrix(1, 1);
		if(pPriority == NULL)
		{
			printf("Error: TileMapv2_HMM_RegionFDR_LeftTail, cannot create sorting key!\n");
			exit(EXIT_FAILURE);
		}
		pPriority->pMatElement[0] = 3;
		pRegionCTMaxSort = NULL;
		pRegionCTMaxSid = NULL;
		DMSORTROWS(pRegionCT, pType, pPriority, &pRegionCTMaxSort, &pRegionCTMaxSid);

		/* sumT */
		pPriority->pMatElement[0] = 7;
		pRegionCTSumSort = NULL;
		pRegionCTSumSid = NULL;
		DMSORTROWS(pRegionCT, pType, pPriority, &pRegionCTSumSort, &pRegionCTSumSid);
		
		DestroyIntMatrix(pPriority);

		/* count */
		nCol = 3;
		TileMapv2_RegionFDR_Count(pRegionMaxSort, pRegionCTMaxSort, pRegionMaxFDR, nCol);

		nCol = 7;
		TileMapv2_RegionFDR_Count(pRegionSumSort, pRegionCTSumSort, pRegionSumFDR, nCol);

		/* -------------- */
		/* release memory */
		/* -------------- */
		DestroyDoubleMatrix(pRegionCTMaxSort);
		DestroyDoubleMatrix(pRegionCTSumSort);
		DestroyLongMatrix(pRegionCTMaxSid);
		DestroyLongMatrix(pRegionCTSumSid);
		DestroyDoubleMatrix(pRegionCT);
	}
	else
	{
		/* FDR == 0*/
	}

	/* FDR */
	TileMapv2_RegionFDR_Compute(pRegionMaxFDR);
	TileMapv2_RegionFDR_Compute(pRegionSumFDR);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_DataTransform()                                              */
/*  Transform column nCol according to strTransform.                       */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_DataTransform(struct tagBARData *pData, int nCol, 
							char strTransform[])
{
	/* define */
	int ni,nj;
	double *vT;

	/* init */
	if( (nCol<0) || (nCol>=pData->nColNum) )
	{
		printf("Error: TileMapv2_DataTransform, column index out of range!\n");
		exit(EXIT_FAILURE);
	}

	/* operation */
	if(strcmp(strTransform, "invlogit") == 0)
	{
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			if(pData->vSeqData[ni]->vData[nCol] != NULL)
			{
				vT = pData->vSeqData[ni]->vData[nCol]->pMatElement;
				for(nj=0; nj<pData->vSeqData[ni]->nDataNum; nj++)
				{
					vT[nj] = exp(vT[nj]);
					vT[nj] = vT[nj]/(1.0+vT[nj]);
				}
			}
		}
	}
	else if(strcmp(strTransform, "oneminusinvlogit") == 0)
	{
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			if(pData->vSeqData[ni]->vData[nCol] != NULL)
			{
				vT = pData->vSeqData[ni]->vData[nCol]->pMatElement;
				for(nj=0; nj<pData->vSeqData[ni]->nDataNum; nj++)
				{
					vT[nj] = exp(vT[nj]);
					vT[nj] = 1.0/(1.0+vT[nj]);
				}
			}
		}
	}
	else if(strcmp(strTransform, "oneminus") == 0)
	{
		for(ni=0; ni<pData->nSeqNum; ni++)
		{
			if(pData->vSeqData[ni]->vData[nCol] != NULL)
			{
				vT = pData->vSeqData[ni]->vData[nCol]->pMatElement;
				for(nj=0; nj<pData->vSeqData[ni]->nDataNum; nj++)
				{
					vT[nj] = 1.0-vT[nj];
				}
			}
		}
	}
	else
	{
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_RegionInfo_Main()                                            */
/*  TileMapv2 get region information                                       */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_RegionInfo_Main(char strRegionPath[], char strInfoPath[], char strOutPath[])
{
	/* ------ */
	/* define */
	/* ------ */
	int nTrackNum;
	int nRegionNum;
	int nRegionId;
	int ni;
	FILE *fpReg;
	FILE *fpOut;
	struct tagTileMapv2InfoTrack *pTrackList = NULL;
	struct tagTileMapv2InfoTrack *pTrack = NULL;
	struct DOUBLEMATRIX *pInfo;
	char strLine[MED_LINE_LENGTH];
	struct tagString *pHeader = NULL;
	char strTemp[MED_LINE_LENGTH];

	/* get region number */
	nRegionNum = 0;
	fpReg = NULL;
	fpReg = fopen(strRegionPath, "r");
	if(fpReg == NULL)
	{
		printf("Error: TileMapv2_RegionInfo_Main, cannot open region file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, MED_LINE_LENGTH, fpReg) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		nRegionNum++;
	}

	fclose(fpReg);

	/* get track info */
	nTrackNum = TileMapv2_RegionInfo_LoadTrackInfo(strInfoPath, &pTrackList);

	if( (nRegionNum<=0) || (nTrackNum<=0) )
	{
		printf("Warning: TileMapv2_RegionInfo_Main, empty regions or tracks!\n");
		return PROC_SUCCESS;
	}

	/* create matrix */
	pInfo = NULL;
	pInfo = CreateDoubleMatrix(nRegionNum, nTrackNum);
	if(pInfo == NULL)
	{
		printf("Error: TileMapv2_RegionInfo_Main, cannot create memory for loading region info!\n");
		exit(EXIT_FAILURE);
	}

	/* load info */
	ni = 0;
	while(pTrackList != NULL)
	{
		printf("Processing track %d ...\n", ni+1);
		pTrack = pTrackList;
		pTrackList = pTrack->pNext;
		pTrack->pNext = NULL;

		if(ni >= nTrackNum)
		{
			printf("Error: TileMapv2_RegionInfo_Main, inconsistent track number!\n");
			exit(EXIT_FAILURE);
		}

		sprintf(strTemp, "\t%s", pTrack->strTrackName);
		StringAddTail(&pHeader, strTemp);

		TileMapv2_RegionInfo_CollectTrackInfo(strRegionPath, pTrack, pInfo, ni);

		TileMapv2_InfoTrack_Destroy(&pTrack);
		ni++;
	}

	if(ni != nTrackNum)
	{
		printf("Error: TileMapv2_RegionInfo_Main, inconsistent track number!\n");
		exit(EXIT_FAILURE);
	}

	/* export results */
	nRegionId = 0;
	fpReg = NULL;
	fpReg = fopen(strRegionPath, "r");
	if(fpReg == NULL)
	{
		printf("Error: TileMapv2_RegionInfo_Main, cannot open region file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: TileMapv2_RegionInfo_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, MED_LINE_LENGTH, fpReg) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
		{
			if(pHeader != NULL)
				fprintf(fpOut, "%s%s\n", strLine, pHeader->m_pString);
			else
				fprintf(fpOut, "%s\n", strLine);
			continue;
		}

		fprintf(fpOut, "%s", strLine);
		for(ni=0; ni<nTrackNum; ni++)
		{
			fprintf(fpOut, "\t%f", DMGETAT(pInfo, nRegionId, ni));
		}
		fprintf(fpOut, "\n");

		nRegionId++;
	}

	fclose(fpReg);
	fclose(fpOut);
	
	if(nRegionId != nRegionNum)
	{
		printf("Error: TileMapv2_RegionInfo_Main, inconsistent region number!\n");
		exit(EXIT_FAILURE);
	}

	/* release memory */
	DestroyDoubleMatrix(pInfo);
	DeleteString(pHeader);
	pHeader = NULL;

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_InfoTrack_Create()                                           */
/*  TileMapv2 create information track.                                    */
/* ----------------------------------------------------------------------- */ 
struct tagTileMapv2InfoTrack *TileMapv2_InfoTrack_Create(int nFilterNum)
{
	/* define */
	struct tagTileMapv2InfoTrack *pTrack = NULL;

	/* create */
	pTrack = (struct tagTileMapv2InfoTrack *)calloc(1, sizeof(struct tagTileMapv2InfoTrack));
	if(pTrack == NULL)
	{
		printf("Error: TileMapv2_InfoTrack_Create, cannot create info track!\n");
		exit(EXIT_FAILURE);
	}

	pTrack->nFilterNum = nFilterNum;
	
	if(nFilterNum > 0)
	{
		pTrack->pFilterType = CreateByteMatrix(1, nFilterNum);
		if(pTrack->pFilterType == NULL)
		{
			printf("Error: TileMapv2_InfoTrack_Create, cannot create info track memory!\n");
			exit(EXIT_FAILURE);
		}

		pTrack->pFilterValue = CreateDoubleMatrix(1, nFilterNum);
		if(pTrack->pFilterValue == NULL)
		{
			printf("Error: TileMapv2_InfoTrack_Create, cannot create info track memory!\n");
			exit(EXIT_FAILURE);
		}

		pTrack->vFilterFile = (struct tagString **)calloc(nFilterNum, sizeof(struct tagString*));
		if(pTrack->vFilterFile == NULL)
		{
			printf("Error: TileMapv2_InfoTrack_Create, cannot create info track memory!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* return */
	return pTrack;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_InfoTrack_Destroy()                                          */
/*  TileMapv2 destroy information track.                                   */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_InfoTrack_Destroy(struct tagTileMapv2InfoTrack **pTrack)
{
	/* define */
	int ni;

	/* release memory */
	if(pTrack == NULL)
		return PROC_SUCCESS;

	if(*pTrack == NULL)
		return PROC_SUCCESS;

	DestroyByteMatrix((*pTrack)->pFilterType);
	DestroyDoubleMatrix((*pTrack)->pFilterValue);
	for(ni=0; ni<(*pTrack)->nFilterNum; ni++)
	{
		DeleteString((*pTrack)->vFilterFile[ni]);
		(*pTrack)->vFilterFile[ni] = NULL;
	}
	free((*pTrack)->vFilterFile);
	free(*pTrack);
	*pTrack = NULL;

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_RegionInfo_LoadTrackInfo()                                   */
/*  TileMapv2 load info track information                                  */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_RegionInfo_LoadTrackInfo(char strInfoPath[], 
				struct tagTileMapv2InfoTrack **pTrackList)
{
	/* define */
	int nTrackNum = 0;
	FILE *fpIn = NULL;
	char strLine[MED_LINE_LENGTH];
	char strTrackName[MED_LINE_LENGTH];
	char *chp;
	int nFilterNum = 0;
	int nFilterId = 0;
	struct tagTileMapv2InfoTrack *pTrack = NULL;
	struct tagTileMapv2InfoTrack *pLastTrack = NULL;
	int nF1 = 0;
	int nF2 = 0;
	int nF3 = 0;

	/* init */
	fpIn = fopen(strInfoPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileMapv2_RegionInfo_LoadTrackInfo, cannot find the information file!\n");
		exit(EXIT_FAILURE);
	}

	if(pTrackList == NULL)
	{
		printf("Error: TileMapv2_RegionInfo_LoadTrackInfo, invalid info pointer!\n");
		exit(EXIT_FAILURE);
	}
	if(*pTrackList != NULL)
	{
		printf("Error: TileMapv2_RegionInfo_LoadTrackInfo, invalid info pointer!\n");
		exit(EXIT_FAILURE);
	}

	/* read file */
	while(fgets(strLine, MED_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strLine[0] == '>')
		{
			/* save old track */
			if(pTrack != NULL)
			{
				if(nFilterId != nFilterNum)
				{
					printf("Error: TileMapv2_RegionInfo_LoadTrackInfo, unrecognized file format!\n");
					exit(EXIT_FAILURE);
				}

				if(*pTrackList == NULL)
				{
					*pTrackList = pTrack;
				}
				else
				{
					pLastTrack->pNext = pTrack;
				}
				pLastTrack = pTrack;
				nTrackNum++;
				pTrack = NULL;
				nFilterId = 0;
				nFilterNum = 0;
				nF1 = 0;
				nF2 = 0;
				nF3 = 0;
			}

			strcpy(strTrackName, strLine+1);
		}
		else if(strstr(strLine, "[Filters]") == strLine)
		{
			chp = strstr(strLine, "=");
			chp++;
			StrTrimLeft(chp);
			nFilterNum = atoi(chp);

			if(pTrack != NULL)
			{
				printf("Error: TileMapv2_RegionInfo_LoadTrackInfo, unrecognized file format!\n");
				exit(EXIT_FAILURE);
			}

			pTrack = TileMapv2_InfoTrack_Create(nFilterNum);
			if(pTrack == NULL)
			{
				printf("Error: TileMapv2_RegionInfo_LoadTrackInfo, cannot create track info object!\n");
				exit(EXIT_FAILURE);
			}
			strcpy(pTrack->strTrackName, strTrackName);
			nFilterId = 0;
		}
		else if(strstr(strLine, "[Data Path]") == strLine)
		{
			chp = strstr(strLine, "=");
			chp++;
			StrTrimLeft(chp);

			if(pTrack == NULL)
			{
				printf("Error: TileMapv2_RegionInfo_LoadTrackInfo, cannot create track info object!\n");
				exit(EXIT_FAILURE);
			}
			strcpy(pTrack->strDataFile, chp);
		}
		else if(strstr(strLine, "[Filter Path]") == strLine)
		{
			chp = strstr(strLine, "=");
			chp++;
			StrTrimLeft(chp);

			if( (pTrack == NULL) || (nFilterId >= nFilterNum) )
			{
				printf("Error: TileMapv2_RegionInfo_LoadTrackInfo, cannot create track info object!\n");
				exit(EXIT_FAILURE);
			}

			StringAddTail(pTrack->vFilterFile+nFilterId, chp);
			nF1 = 1;
			if( (nF1 == 1) && (nF2 == 1) && (nF3 ==1) )
			{
				nFilterId++;
				nF1 = 0;
				nF2 = 0;
				nF3 = 0;
			}
		}
		else if(strstr(strLine, "[Filter Type]") == strLine)
		{
			chp = strstr(strLine, "=");
			chp++;
			StrTrimLeft(chp);

			if( (pTrack == NULL) || (nFilterId >= nFilterNum) )
			{
				printf("Error: TileMapv2_RegionInfo_LoadTrackInfo, cannot create track info object!\n");
				exit(EXIT_FAILURE);
			}

			if(*chp == '<')
				pTrack->pFilterType->pMatElement[nFilterId] = 1;
			else
				pTrack->pFilterType->pMatElement[nFilterId] = 0;

			nF2 = 1;
			if( (nF1 == 1) && (nF2 == 1) && (nF3 ==1) )
			{
				nFilterId++;
				nF1 = 0;
				nF2 = 0;
				nF3 = 0;
			}
		}
		else if(strstr(strLine, "[Filter Value]") == strLine)
		{
			chp = strstr(strLine, "=");
			chp++;
			StrTrimLeft(chp);

			if( (pTrack == NULL) || (nFilterId >= nFilterNum) )
			{
				printf("Error: TileMapv2_RegionInfo_LoadTrackInfo, cannot create track info object!\n");
				exit(EXIT_FAILURE);
			}

			pTrack->pFilterValue->pMatElement[nFilterId] = atof(chp);
			nF3 = 1;
			if( (nF1 == 1) && (nF2 == 1) && (nF3 ==1) )
			{
				nFilterId++;
				nF1 = 0;
				nF2 = 0;
				nF3 = 0;
			}
		}
		else
		{
			printf("Error: TileMapv2_RegionInfo_LoadTrackInfo, unrecognized file format!\n");
			exit(EXIT_FAILURE);
		}
	}

	if(pTrack != NULL)
	{
		if(nFilterId != nFilterNum)
		{
			printf("Error: TileMapv2_RegionInfo_LoadTrackInfo, unrecognized file format!\n");
			exit(EXIT_FAILURE);
		}

		if(*pTrackList == NULL)
		{
			*pTrackList = pTrack;
		}
		else
		{
			pLastTrack->pNext = pTrack;
		}
		pLastTrack = pTrack;
		nTrackNum++;
		pTrack = NULL;
		nFilterId = 0;
		nFilterNum = 0;
		nF1 = 0;
		nF2 = 0;
		nF3 = 0;
	}

	/* close file */
	fclose(fpIn);

	/* return */
	return nTrackNum;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_RegionInfo_CollectTrackInfo()                                */
/*  TileMapv2 collect track information                                    */
/* ----------------------------------------------------------------------- */
int TileMapv2_RegionInfo_CollectTrackInfo(char strRegionPath[], 
			struct tagTileMapv2InfoTrack *pTrack, 
			struct DOUBLEMATRIX *pInfo, int nCol)
{
	/* define */
	struct tagBARData *pData = NULL;
	struct tagBARData *pMask = NULL;
	struct INTMATRIX *pFieldType = NULL;
	int nDataCol;
	int ni,nj;
	FILE *fpReg;
	int nRegionId;
	char strAlias[LINE_LENGTH];
	char strChr[LINE_LENGTH];
	int nStart,nEnd;
	int nSeqId,nSeqP1,nSeqP2;
	double dValue;
	int nPassFilter;
	char strSeqName[LINE_LENGTH];
	char strLine[MED_LINE_LENGTH];
	int nP0,nP1;
	struct DOUBLEMATRIX ** vDataVec;

	/* init */
	if(pTrack == NULL)
	{
		printf("Warning: TileMapv2_RegionInfo_CollectTrackInfo, empty track!\n");
		return PROC_SUCCESS;
	}

	/* -------------------- */
	/* load probe level t   */
	/* -------------------- */
    pData = NULL;
	pData = Affy_LoadBAR_Fast(pTrack->strDataFile);
	if(pData == NULL)
	{
		printf("Error: TileMapv2_RegionInfo_CollectTrackInfo, cannot load track data!\n");
		exit(EXIT_FAILURE);
	}

	/* nColNum + nTrackFilters */ 
	nDataCol = pData->nColNum;
	pData->nColNum =pData->nColNum+pTrack->nFilterNum;
	pFieldType = pData->pFieldType;
	pData->pFieldType = NULL;
	pData->pFieldType = CreateIntMatrix(1,pData->nColNum);
	if(pData->pFieldType == NULL)
	{
		printf("Error: TileMapv2_RegionInfo_CollectTrackInfo, cannot create memory for field types!\n");
		exit(EXIT_FAILURE);
	}
	memcpy(pData->pFieldType->pMatElement, pFieldType->pMatElement, sizeof(int)*pFieldType->nWidth);
	DestroyIntMatrix(pFieldType);
	
	for(nj=0; nj<pData->nSeqNum; nj++)
	{
		pData->vSeqData[nj]->nColNum = pData->nColNum;
		vDataVec = pData->vSeqData[nj]->vData;
		pData->vSeqData[nj]->vData = NULL;
		pData->vSeqData[nj]->vData = (struct DOUBLEMATRIX **)calloc(pData->vSeqData[nj]->nColNum, sizeof(struct DOUBLEMATRIX *));
		if(pData->vSeqData[nj]->vData == NULL)
		{
			printf("Error: TileMapv2_RegionInfo_CollectTrackInfo, cannot create memory for tracking intensity data!\n");
			exit(EXIT_FAILURE);
		}
		for(ni=0; ni<nDataCol; ni++)
			pData->vSeqData[nj]->vData[ni] = vDataVec[ni];
		free(vDataVec);
	}

	for(ni=0; ni<pTrack->nFilterNum; ni++)
	{
		if(pTrack->vFilterFile[ni] == NULL)
			continue;

		pMask = NULL;
		pMask = Affy_LoadBAR_Fast(pTrack->vFilterFile[ni]->m_pString);
		if(pMask == NULL)
		{
			printf("Error: TileMapv2_RegionInfo_CollectTrackInfo, cannot load mask data!\n");
			exit(EXIT_FAILURE);
		}

		if(pMask->nSeqNum != pData->nSeqNum)
		{
			printf("Error: TileMapv2_RegionInfo_CollectTrackInfo, array types do not match!\n");
			exit(EXIT_FAILURE);
		}

		for(nj=0; nj<pMask->nSeqNum; nj++)
		{
			if(pMask->vSeqData[nj]->nDataNum != pData->vSeqData[nj]->nDataNum)
			{
				printf("Error: TileMapv2_RegionInfo_CollectTrackInfo, array types do not match!\n");
				exit(EXIT_FAILURE);
			}
			pData->vSeqData[nj]->vData[nDataCol+ni] = pMask->vSeqData[nj]->vData[1];
			pMask->vSeqData[nj]->vData[1] = NULL;
		}
        
		pData->pFieldType->pMatElement[nDataCol+ni] = pMask->pFieldType->pMatElement[1];

		Affy_BARData_Destroy(&pMask);
	}


	/* process regions */
	nRegionId = 0;
	fpReg = NULL;
	fpReg = fopen(strRegionPath, "r");
	if(fpReg == NULL)
	{
		printf("Error: TileMapv2_RegionInfo_CollectTrackInfo, cannot open region file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, MED_LINE_LENGTH, fpReg) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		sscanf(strLine, "%s %s %d %d", strAlias, strChr, &nStart, &nEnd);

		dValue = -1e6;

		for(nSeqId=0; nSeqId<pData->nSeqNum; nSeqId++)
		{
			/* find seq */
			sprintf(strSeqName, "%s", pData->vSeqData[nSeqId]->pSeqName->m_pString);
			if(pData->fVersionnumber > 1.5)
			{
				if(pData->vSeqData[nSeqId]->pSeqGroupName != NULL)
				{
					if(pData->vSeqData[nSeqId]->pSeqVersion != NULL)
					{
						if(strcmp(pData->vSeqData[nSeqId]->pSeqVersion->m_pString, "") == 0)
						{
							sprintf(strSeqName, "%s", pData->vSeqData[nSeqId]->pSeqName->m_pString);
						}
						else if(strstr(pData->vSeqData[nSeqId]->pSeqVersion->m_pString, "NCBI") == pData->vSeqData[nSeqId]->pSeqVersion->m_pString)
						{
							sprintf(strSeqName, "%s", pData->vSeqData[nSeqId]->pSeqName->m_pString);
						}
						else
						{
							sprintf(strSeqName, "");
						}
					}
				}
			}

			/* if(pData->vSeqData[nSeqId]->pSeqVersion != NULL)
			{
				sprintf(strSeqName, "%s:%s", pData->vSeqData[nSeqId]->pSeqVersion->m_pString,
					pData->vSeqData[nSeqId]->pSeqName->m_pString);
			}
			else
			{
				sprintf(strSeqName, "%s", pData->vSeqData[nSeqId]->pSeqName->m_pString);
			}
			if(pData->fVersionnumber > 1.5)
			{
				if(pData->vSeqData[nSeqId]->pSeqGroupName != NULL)
				{
					if(pData->vSeqData[nSeqId]->pSeqVersion != NULL)
					{
						sprintf(strSeqName, "%s:%s:%s", pData->vSeqData[nSeqId]->pSeqGroupName->m_pString,
							pData->vSeqData[nSeqId]->pSeqVersion->m_pString,
							pData->vSeqData[nSeqId]->pSeqName->m_pString);
					}
					else
					{
						sprintf(strSeqName, "%s:%s", pData->vSeqData[nSeqId]->pSeqGroupName->m_pString,
							pData->vSeqData[nSeqId]->pSeqName->m_pString);
					}
				}	
			} */

			if(strcmp(strSeqName, strChr) != 0)
			{
				continue;
			}
			if(pData->vSeqData[nSeqId]->nDataNum <= 0)
			{
				continue;
			}

			/* find start */
			nP0 = 0;
			nP1 = pData->vSeqData[nSeqId]->nDataNum-1;
			
			if(nStart <= (int)(pData->vSeqData[nSeqId]->vData[0]->pMatElement[nP0]))
			{
				nSeqP1 = nP0;
			}
			else if(nStart > (int)(pData->vSeqData[nSeqId]->vData[0]->pMatElement[nP1]))
			{
				nSeqP1 = pData->vSeqData[nSeqId]->nDataNum;
			}
			else
			{
				while((nP1-nP0) > 1)
				{
					nSeqP1 = (nP0+nP1)/2;
					if((int)(pData->vSeqData[nSeqId]->vData[0]->pMatElement[nSeqP1]) >= nStart)
					{
						nP1 = nSeqP1;
					}
					else
					{
						nP0 = nSeqP1;
					}
				}
				
				nSeqP1 = nP1;
			}

			/* find end */
			nP0 = 0;
			nP1 = pData->vSeqData[nSeqId]->nDataNum-1;
			
			if(nEnd < (int)(pData->vSeqData[nSeqId]->vData[0]->pMatElement[nP0]))
			{
				nSeqP2 = -1;
			}
			else if(nEnd >= (int)(pData->vSeqData[nSeqId]->vData[0]->pMatElement[nP1]))
			{
				nSeqP2 = nP1;
			}
			else
			{
				while((nP1-nP0) > 1)
				{
					nSeqP2 = (nP0+nP1)/2;
					if((int)(pData->vSeqData[nSeqId]->vData[0]->pMatElement[nSeqP2]) > nEnd)
					{
						nP1 = nSeqP2;
					}
					else
					{
						nP0 = nSeqP2;
					}
				}
				
				nSeqP2 = nP0;
			}

			/* get info */
			for(ni=nSeqP1; ni<=nSeqP2; ni++)
			{
				nPassFilter = 1;
				for(nj=0; nj<pTrack->nFilterNum; nj++)
				{
					if(pData->vSeqData[nSeqId]->vData[nDataCol+nj] == NULL)
						continue;

					if(pTrack->pFilterType->pMatElement[nj] == 1)
					{
						if(pData->vSeqData[nSeqId]->vData[nDataCol+nj]->pMatElement[ni] >= pTrack->pFilterValue->pMatElement[nj])
						{
							nPassFilter = 0;
							break;
						}
					}
					else
					{
						if(pData->vSeqData[nSeqId]->vData[nDataCol+nj]->pMatElement[ni] <= pTrack->pFilterValue->pMatElement[nj])
						{
							nPassFilter = 0;
							break;
						}
					}
				}

				if(nPassFilter == 0)
					continue;

				if(pData->vSeqData[nSeqId]->vData[1]->pMatElement[ni] > dValue)
					dValue = pData->vSeqData[nSeqId]->vData[1]->pMatElement[ni];
			}

		}

		DMSETAT(pInfo, nRegionId, nCol, dValue);
		nRegionId++;
	}

	fclose(fpReg);

	/* release memory */
	Affy_BARData_Destroy(&pData);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_RegionInfo_Integral_Main()                                   */
/*  TileMapv2 get region information                                       */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_RegionInfo_Integral_Main(char strRegionPath[], char strInfoPath[], char strOutPath[])
{
	/* ------ */
	/* define */
	/* ------ */
	int nTrackNum;
	int nRegionNum;
	int nRegionId;
	int ni;
	FILE *fpReg;
	FILE *fpOut;
	struct tagTileMapv2InfoTrack *pTrackList = NULL;
	struct tagTileMapv2InfoTrack *pTrack = NULL;
	struct DOUBLEMATRIX *pInfo;
	char strLine[MED_LINE_LENGTH];
	struct tagString *pHeader = NULL;
	char strTemp[MED_LINE_LENGTH];

	/* get region number */
	nRegionNum = 0;
	fpReg = NULL;
	fpReg = fopen(strRegionPath, "r");
	if(fpReg == NULL)
	{
		printf("Error: TileMapv2_RegionInfo_Integral_Main, cannot open region file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, MED_LINE_LENGTH, fpReg) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		nRegionNum++;
	}

	fclose(fpReg);

	/* get track info */
	nTrackNum = TileMapv2_RegionInfo_LoadTrackInfo(strInfoPath, &pTrackList);

	if( (nRegionNum<=0) || (nTrackNum<=0) )
	{
		printf("Warning: TileMapv2_RegionInfo_Integral_Main, empty regions or tracks!\n");
		return PROC_SUCCESS;
	}

	/* create matrix */
	pInfo = NULL;
	pInfo = CreateDoubleMatrix(nRegionNum, 3*nTrackNum);
	if(pInfo == NULL)
	{
		printf("Error: TileMapv2_RegionInfo_Integral_Main, cannot create memory for loading region info!\n");
		exit(EXIT_FAILURE);
	}

	/* load info */
	ni = 0;
	while(pTrackList != NULL)
	{
		printf("Processing track %d ...\n", ni+1);
		pTrack = pTrackList;
		pTrackList = pTrack->pNext;
		pTrack->pNext = NULL;

		if(ni >= nTrackNum)
		{
			printf("Error: TileMapv2_RegionInfo_Integral_Main, inconsistent track number!\n");
			exit(EXIT_FAILURE);
		}

		sprintf(strTemp, "\t%s_sum\t%s_probenum\t%s_mean", pTrack->strTrackName, pTrack->strTrackName, pTrack->strTrackName);
		StringAddTail(&pHeader, strTemp);

		TileMapv2_RegionInfo_Integral_CollectTrackInfo(strRegionPath, pTrack, pInfo, ni);

		TileMapv2_InfoTrack_Destroy(&pTrack);
		ni++;
	}

	if(ni != nTrackNum)
	{
		printf("Error: TileMapv2_RegionInfo_Integral_Main, inconsistent track number!\n");
		exit(EXIT_FAILURE);
	}

	/* export results */
	nRegionId = 0;
	fpReg = NULL;
	fpReg = fopen(strRegionPath, "r");
	if(fpReg == NULL)
	{
		printf("Error: TileMapv2_RegionInfo_Integral_Main, cannot open region file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: TileMapv2_RegionInfo_Integral_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, MED_LINE_LENGTH, fpReg) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
		{
			if(pHeader != NULL)
				fprintf(fpOut, "%s%s\n", strLine, pHeader->m_pString);
			else
				fprintf(fpOut, "%s\n", strLine);
			continue;
		}

		fprintf(fpOut, "%s", strLine);
		for(ni=0; ni<nTrackNum; ni++)
		{
			fprintf(fpOut, "\t%f\t%d\t%f", DMGETAT(pInfo, nRegionId, 3*ni), (int)(DMGETAT(pInfo, nRegionId, 3*ni+1)), DMGETAT(pInfo, nRegionId, 3*ni+2));
		}
		fprintf(fpOut, "\n");

		nRegionId++;
	}

	fclose(fpReg);
	fclose(fpOut);
	
	if(nRegionId != nRegionNum)
	{
		printf("Error: TileMapv2_RegionInfo_Integral_Main, inconsistent region number!\n");
		exit(EXIT_FAILURE);
	}

	/* release memory */
	DestroyDoubleMatrix(pInfo);
	DeleteString(pHeader);
	pHeader = NULL;

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_RegionInfo_Integral_CollectTrackInfo()                       */
/*  TileMapv2 collect track information                                    */
/* ----------------------------------------------------------------------- */
int TileMapv2_RegionInfo_Integral_CollectTrackInfo(char strRegionPath[], 
			struct tagTileMapv2InfoTrack *pTrack, 
			struct DOUBLEMATRIX *pInfo, int nCol)
{
	/* define */
	struct tagBARData *pData = NULL;
	struct tagBARData *pMask = NULL;
	struct INTMATRIX *pFieldType = NULL;
	int nDataCol;
	int ni,nj;
	FILE *fpReg;
	int nRegionId;
	char strAlias[LINE_LENGTH];
	char strChr[LINE_LENGTH];
	int nStart,nEnd;
	int nSeqId,nSeqP1,nSeqP2;
	double dValue,dMean;
	int nNum;
	int nPassFilter;
	char strSeqName[LINE_LENGTH];
	char strLine[MED_LINE_LENGTH];
	int nP0,nP1;
	struct DOUBLEMATRIX ** vDataVec;

	/* init */
	if(pTrack == NULL)
	{
		printf("Warning: TileMapv2_RegionInfo_Integral_CollectTrackInfo, empty track!\n");
		return PROC_SUCCESS;
	}

	/* -------------------- */
	/* load probe level t   */
	/* -------------------- */
    pData = NULL;
	pData = Affy_LoadBAR_Fast(pTrack->strDataFile);
	if(pData == NULL)
	{
		printf("Error: TileMapv2_RegionInfo_Integral_CollectTrackInfo, cannot load track data!\n");
		exit(EXIT_FAILURE);
	}

	/* nColNum + nTrackFilters */ 
	nDataCol = pData->nColNum;
	pData->nColNum =pData->nColNum+pTrack->nFilterNum;
	pFieldType = pData->pFieldType;
	pData->pFieldType = NULL;
	pData->pFieldType = CreateIntMatrix(1,pData->nColNum);
	if(pData->pFieldType == NULL)
	{
		printf("Error: TileMapv2_RegionInfo_CollectTrackInfo, cannot create memory for field types!\n");
		exit(EXIT_FAILURE);
	}
	memcpy(pData->pFieldType->pMatElement, pFieldType->pMatElement, sizeof(int)*pFieldType->nWidth);
	DestroyIntMatrix(pFieldType);
	
	for(nj=0; nj<pData->nSeqNum; nj++)
	{
		pData->vSeqData[nj]->nColNum = pData->nColNum;
		vDataVec = pData->vSeqData[nj]->vData;
		pData->vSeqData[nj]->vData = NULL;
		pData->vSeqData[nj]->vData = (struct DOUBLEMATRIX **)calloc(pData->vSeqData[nj]->nColNum, sizeof(struct DOUBLEMATRIX *));
		if(pData->vSeqData[nj]->vData == NULL)
		{
			printf("Error: TileMapv2_RegionInfo_CollectTrackInfo, cannot create memory for tracking intensity data!\n");
			exit(EXIT_FAILURE);
		}
		for(ni=0; ni<nDataCol; ni++)
			pData->vSeqData[nj]->vData[ni] = vDataVec[ni];
		free(vDataVec);
	}

	for(ni=0; ni<pTrack->nFilterNum; ni++)
	{
		if(pTrack->vFilterFile[ni] == NULL)
			continue;

		pMask = NULL;
		pMask = Affy_LoadBAR_Fast(pTrack->vFilterFile[ni]->m_pString);
		if(pMask == NULL)
		{
			printf("Error: TileMapv2_RegionInfo_CollectTrackInfo, cannot load mask data!\n");
			exit(EXIT_FAILURE);
		}

		if(pMask->nSeqNum != pData->nSeqNum)
		{
			printf("Error: TileMapv2_RegionInfo_CollectTrackInfo, array types do not match!\n");
			exit(EXIT_FAILURE);
		}

		for(nj=0; nj<pMask->nSeqNum; nj++)
		{
			if(pMask->vSeqData[nj]->nDataNum != pData->vSeqData[nj]->nDataNum)
			{
				printf("Error: TileMapv2_RegionInfo_CollectTrackInfo, array types do not match!\n");
				exit(EXIT_FAILURE);
			}
			pData->vSeqData[nj]->vData[nDataCol+ni] = pMask->vSeqData[nj]->vData[1];
			pMask->vSeqData[nj]->vData[1] = NULL;
		}
        
		pData->pFieldType->pMatElement[nDataCol+ni] = pMask->pFieldType->pMatElement[1];

		Affy_BARData_Destroy(&pMask);
	}


	/* process regions */
	nRegionId = 0;
	fpReg = NULL;
	fpReg = fopen(strRegionPath, "r");
	if(fpReg == NULL)
	{
		printf("Error: TileMapv2_RegionInfo_CollectTrackInfo, cannot open region file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, MED_LINE_LENGTH, fpReg) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		sscanf(strLine, "%s %s %d %d", strAlias, strChr, &nStart, &nEnd);

		dValue = 0.0;
		nNum = 0;
		dMean = 0.0;

		for(nSeqId=0; nSeqId<pData->nSeqNum; nSeqId++)
		{
			/* find seq */
			sprintf(strSeqName, "%s", pData->vSeqData[nSeqId]->pSeqName->m_pString);
			if(pData->fVersionnumber > 1.5)
			{
				if(pData->vSeqData[nSeqId]->pSeqGroupName != NULL)
				{
					if(pData->vSeqData[nSeqId]->pSeqVersion != NULL)
					{
						if(strcmp(pData->vSeqData[nSeqId]->pSeqVersion->m_pString, "") == 0)
						{
							sprintf(strSeqName, "%s", pData->vSeqData[nSeqId]->pSeqName->m_pString);
						}
						else if(strstr(pData->vSeqData[nSeqId]->pSeqVersion->m_pString, "NCBI") == pData->vSeqData[nSeqId]->pSeqVersion->m_pString)
						{
							sprintf(strSeqName, "%s", pData->vSeqData[nSeqId]->pSeqName->m_pString);
						}
						else
						{
							sprintf(strSeqName, "");
						}
					}
				}
			}

			/* if(pData->vSeqData[nSeqId]->pSeqVersion != NULL)
			{
				sprintf(strSeqName, "%s:%s", pData->vSeqData[nSeqId]->pSeqVersion->m_pString,
					pData->vSeqData[nSeqId]->pSeqName->m_pString);
			}
			else
			{
				sprintf(strSeqName, "%s", pData->vSeqData[nSeqId]->pSeqName->m_pString);
			}
			if(pData->fVersionnumber > 1.5)
			{
				if(pData->vSeqData[nSeqId]->pSeqGroupName != NULL)
				{
					if(pData->vSeqData[nSeqId]->pSeqVersion != NULL)
					{
						sprintf(strSeqName, "%s:%s:%s", pData->vSeqData[nSeqId]->pSeqGroupName->m_pString,
							pData->vSeqData[nSeqId]->pSeqVersion->m_pString,
							pData->vSeqData[nSeqId]->pSeqName->m_pString);
					}
					else
					{
						sprintf(strSeqName, "%s:%s", pData->vSeqData[nSeqId]->pSeqGroupName->m_pString,
							pData->vSeqData[nSeqId]->pSeqName->m_pString);
					}
				}	
			} */

			if(strcmp(strSeqName, strChr) != 0)
			{
				continue;
			}
			if(pData->vSeqData[nSeqId]->nDataNum <= 0)
			{
				continue;
			}

			/* find start */
			nP0 = 0;
			nP1 = pData->vSeqData[nSeqId]->nDataNum-1;
			
			if(nStart <= (int)(pData->vSeqData[nSeqId]->vData[0]->pMatElement[nP0]))
			{
				nSeqP1 = nP0;
			}
			else if(nStart > (int)(pData->vSeqData[nSeqId]->vData[0]->pMatElement[nP1]))
			{
				nSeqP1 = pData->vSeqData[nSeqId]->nDataNum;
			}
			else
			{
				while((nP1-nP0) > 1)
				{
					nSeqP1 = (nP0+nP1)/2;
					if((int)(pData->vSeqData[nSeqId]->vData[0]->pMatElement[nSeqP1]) >= nStart)
					{
						nP1 = nSeqP1;
					}
					else
					{
						nP0 = nSeqP1;
					}
				}
				
				nSeqP1 = nP1;
			}

			/* find end */
			nP0 = 0;
			nP1 = pData->vSeqData[nSeqId]->nDataNum-1;
			
			if(nEnd < (int)(pData->vSeqData[nSeqId]->vData[0]->pMatElement[nP0]))
			{
				nSeqP2 = -1;
			}
			else if(nEnd >= (int)(pData->vSeqData[nSeqId]->vData[0]->pMatElement[nP1]))
			{
				nSeqP2 = nP1;
			}
			else
			{
				while((nP1-nP0) > 1)
				{
					nSeqP2 = (nP0+nP1)/2;
					if((int)(pData->vSeqData[nSeqId]->vData[0]->pMatElement[nSeqP2]) > nEnd)
					{
						nP1 = nSeqP2;
					}
					else
					{
						nP0 = nSeqP2;
					}
				}
				
				nSeqP2 = nP0;
			}

			/* get info */
			for(ni=nSeqP1; ni<=nSeqP2; ni++)
			{
				nPassFilter = 1;
				for(nj=0; nj<pTrack->nFilterNum; nj++)
				{
					if(pData->vSeqData[nSeqId]->vData[nDataCol+nj] == NULL)
						continue;

					if(pTrack->pFilterType->pMatElement[nj] == 1)
					{
						if(pData->vSeqData[nSeqId]->vData[nDataCol+nj]->pMatElement[ni] >= pTrack->pFilterValue->pMatElement[nj])
						{
							nPassFilter = 0;
							break;
						}
					}
					else
					{
						if(pData->vSeqData[nSeqId]->vData[nDataCol+nj]->pMatElement[ni] <= pTrack->pFilterValue->pMatElement[nj])
						{
							nPassFilter = 0;
							break;
						}
					}
				}

				if(nPassFilter == 0)
					continue;

				dValue += pData->vSeqData[nSeqId]->vData[1]->pMatElement[ni];
				nNum += 1;
			}

		}

		if(nNum > 0)
			dMean = dValue/(double)(nNum);
		else
			dMean = 0.0;
		DMSETAT(pInfo, nRegionId, 3*nCol, dValue);
		DMSETAT(pInfo, nRegionId, 3*nCol+1, nNum);
		DMSETAT(pInfo, nRegionId, 3*nCol+2, dMean);
		nRegionId++;
	}

	fclose(fpReg);

	/* release memory */
	Affy_BARData_Destroy(&pData);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_CollectProbes_Main()                                         */
/*  TileMapv2 get probe information                                        */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_CollectProbes_Main(char strRegionPath[], char strInfoPath[], char strOutPath[])
{
	/* ------ */
	/* define */
	/* ------ */
	int nTrackNum;
	int nRegionId;
	int ni;
	FILE *fpReg;
	FILE *fpOut;
	struct tagTileMapv2InfoTrack *pTrackList = NULL;
	struct tagTileMapv2InfoTrack *pTrack = NULL;
	char strLine[LONG_LINE_LENGTH];
	/* struct tagString *pHeader = NULL; */
	char strTemp[MED_LINE_LENGTH];
	char strChr[MED_LINE_LENGTH];
	int nStart,nEnd;
	struct tagBARData **vData;

	
	/* get track info */
	nTrackNum = TileMapv2_RegionInfo_LoadTrackInfo(strInfoPath, &pTrackList);

	if( nTrackNum<=0 )
	{
		printf("Warning: TileMapv2_CollectProbes_Main, empty tracks!\n");
		return PROC_SUCCESS;
	}

	vData = NULL;
	vData = (struct tagBARData **)calloc(nTrackNum, sizeof(struct tagBARData *));
	if(vData == NULL)
	{
		printf("Error: TileMapv2_CollectProbes_Main, cannot allocate memory!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	pTrack = pTrackList;
	while(pTrack != NULL)
	{
		printf("Loading %s...\n", pTrack->strDataFile);
		vData[ni] = Affy_LoadBAR_Fast(pTrack->strDataFile);
		if(vData[ni] == NULL)
		{
			printf("Error: TileMapv2_CollectProbes_Main, cannot load data!\n");
			exit(EXIT_FAILURE);
		}
		pTrack = pTrack->pNext;
		ni++;
	}
	if(ni != nTrackNum)
	{
		printf("Error: TileMapv2_CollectProbes_Main, inconsistent track number!\n");
		exit(EXIT_FAILURE);
	}

	/* load region */
	fpReg = NULL;
	fpReg = fopen(strRegionPath, "r");
	if(fpReg == NULL)
	{
		printf("Error: TileMapv2_CollectProbes_Main, cannot open region file!\n");
		exit(EXIT_FAILURE);
	}

	nRegionId = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpReg) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		sscanf(strLine, "%s %s %d %d", strTemp, strChr, &nStart, &nEnd);
		sprintf(strTemp, "%s_%d.txt", strOutPath, nRegionId);
		fpOut = NULL;
		fpOut = fopen(strTemp, "w");
		if(fpOut == NULL)
		{
			printf("Error: TileMapv2_CollectProbes_Main, cannot open output file!\n");
			exit(EXIT_FAILURE);
		}
		
		TileMapv2_CollectProbes_ProcessRegion(strChr, nStart, nEnd, vData, nTrackNum, fpOut);
		
		fclose(fpOut);

		nRegionId++;
	}

	fclose(fpReg);


	/* release memory */
	while(pTrackList != NULL)
	{
		pTrack = pTrackList;
		pTrackList = pTrack->pNext;
		pTrack->pNext = NULL;
		TileMapv2_InfoTrack_Destroy(&pTrack);
	}

	for(ni=0; ni<nTrackNum; ni++)
	{
		Affy_BARData_Destroy(vData+ni);
	}
	free(vData);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_CollectProbes_ProcessRegion()                                */
/*  TileMapv2 get probe information                                        */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_CollectProbes_ProcessRegion(char strChr[], int nStart, int nEnd, 
							struct tagBARData **vData, int nTrackNum, FILE *fpOut)
{
	/* define */
	int ni,nk;
	int nSID;
	int nP0,nP1,nSeqP1,nSeqP2;

	/* init check */
	if(nTrackNum <= 0)
		return PROC_SUCCESS;


	/* search */
	nSID = -1;
	for(ni=0; ni<vData[0]->nSeqNum; ni++)
	{
		if(strcmp(strChr, vData[0]->vSeqData[ni]->pSeqName->m_pString) == 0)
		{
			nSID = ni;
			break;
		}
	}

	if(nSID < 0)
	{
		printf("Warning: no probes have been found for %s:%d-%d!\n", strChr, nStart, nEnd);
		return PROC_SUCCESS;
	}

	for(ni=0; ni<nTrackNum; ni++)
	{
		if(strcmp(strChr, vData[ni]->vSeqData[nSID]->pSeqName->m_pString) != 0)
		{
			printf("Error: bar files are not consistent.\n");
			exit(EXIT_FAILURE);
		}
	}

	if(vData[0]->vSeqData[nSID]->nDataNum <= 0)
		return PROC_SUCCESS;
	

	/* find start */
	nP0 = 0;
	nP1 = vData[0]->vSeqData[nSID]->nDataNum-1;
	
	if(nStart <= (int)(vData[0]->vSeqData[nSID]->vData[0]->pMatElement[nP0]))
	{
		nSeqP1 = nP0;
	}
	else if(nStart > (int)(vData[0]->vSeqData[nSID]->vData[0]->pMatElement[nP1]))
	{
		nSeqP1 = vData[0]->vSeqData[nSID]->nDataNum;
	}
	else
	{
		while((nP1-nP0) > 1)
		{
			nSeqP1 = (nP0+nP1)/2;
			if((int)(vData[0]->vSeqData[nSID]->vData[0]->pMatElement[nSeqP1]) >= nStart)
			{
				nP1 = nSeqP1;
			}
			else
			{
				nP0 = nSeqP1;
			}
		}
		
		nSeqP1 = nP1;
	}

	/* find end */
	nP0 = 0;
	nP1 = vData[0]->vSeqData[nSID]->nDataNum-1;
	
	if(nEnd < (int)(vData[0]->vSeqData[nSID]->vData[0]->pMatElement[nP0]))
	{
		nSeqP2 = -1;
	}
	else if(nEnd >= (int)(vData[0]->vSeqData[nSID]->vData[0]->pMatElement[nP1]))
	{
		nSeqP2 = nP1;
	}
	else
	{
		while((nP1-nP0) > 1)
		{
			nSeqP2 = (nP0+nP1)/2;
			if((int)(vData[0]->vSeqData[nSID]->vData[0]->pMatElement[nSeqP2]) > nEnd)
			{
				nP1 = nSeqP2;
			}
			else
			{
				nP0 = nSeqP2;
			}
		}
		
		nSeqP2 = nP0;
	}

	if(nSeqP1 > nSeqP2)
		return PROC_SUCCESS;

	for(nk=nSeqP1; nk<=nSeqP2; nk++)
	{
		fprintf(fpOut, "%d", (int)(vData[0]->vSeqData[nSID]->vData[0]->pMatElement[nk]));
		for(ni=0; ni<nTrackNum; ni++)
		{
			if((int)(vData[ni]->vSeqData[nSID]->vData[0]->pMatElement[nk]) != (int)(vData[0]->vSeqData[nSID]->vData[0]->pMatElement[nk]))
			{
				printf("Error: bar files are not consistent.\n");
				exit(EXIT_FAILURE);
			}
			fprintf(fpOut, "\t%f", vData[ni]->vSeqData[nSID]->vData[1]->pMatElement[nk]);
		}
		fprintf(fpOut, "\n");
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_TXT2BAR()                                                    */
/*  Convert a text file to BAR tiling array project.                       */
/* ----------------------------------------------------------------------- */
int TileMapv2_TXT2BAR(char strTXTFile[], char strExportFolder[], char strProjectTitle[])
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
	struct tagString **vSampleName = NULL;
	struct INTMATRIX *vGroupId = NULL;
	int nTotalProbeNum = 0;
	char strLine[LONG_LINE_LENGTH];
	char strTemp[LONG_LINE_LENGTH];
	int nAnnotLine = 0;
	char *chp1,*chp2;
	char strChr[MED_LINE_LENGTH];
	char strLastChr[MED_LINE_LENGTH];
	char strFileName[LONG_LINE_LENGTH];
	struct INTMATRIX *pCol = NULL;
	
	/* count */
	int ni,nj,nk;
	int nGroupNum = 0;
	int nGroupIdWrong = 0;
	struct INTMATRIX *pGroupOK;
	
	/* load */
	AdjustDirectoryPath(strExportFolder);

	fpIn = NULL;
	fpIn = fopen(strTXTFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileMapv2_TXT2BAR, cannot open source file!\n");
		exit(EXIT_FAILURE);
	}

	/* load input */
	pSeqInfo = NULL;
	strcpy(strLastChr, "");
	nProbeNum = 0;
	nColNum = 0;
	nSeqNum = 0;
	nFieldNum = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '0')
			continue;
		if(strLine[0] == '#')
		{
			if(nAnnotLine == 0)
			{
				/* count number of columns */
				chp1 = strLine;
				chp2 = strchr(chp1, '\t');
				while(chp2 != NULL)
				{
					nFieldNum++;
					chp1 = chp2+1;
					chp2 = strchr(chp1, '\t');
				}
				nFieldNum++;
				nColNum = nFieldNum-2;
				if(nColNum < 0)
				{
					printf("Error: TileMapv2_TXT2BAR, wrong input file format, too few columns!\n");
					exit(EXIT_FAILURE);
				}

				if(nColNum == 0)
					nColNum++;

				/* read sample name */
				vSampleName = NULL;
				vSampleName = (struct tagString **)calloc(nColNum, sizeof(struct tagString *));
				if(vSampleName == NULL)
				{
					printf("Error: TileMapv2_TXT2BAR, cannot create memory for sample name!\n");
					exit(EXIT_FAILURE);
				}

				if(nFieldNum == 2)
				{
					sprintf(strTemp, "%s", strProjectTitle);
					StringAddTail(vSampleName+0, strTemp);
				}
				else
				{
					chp1 = strLine;
					for(ni=0; ni<2; ni++)
					{
						chp2 = strchr(chp1, '\t');
						chp1 = chp2+1;
					}
					
					ni  = 0;
					chp2 = strchr(chp1, '\t');
					while(chp2 != NULL)
					{
						*chp2 = '\0';
						StrTrimRight(chp1);
						strcpy(strTemp, chp1);
						StringAddTail(vSampleName+ni, strTemp);

						chp1 = chp2+1;
						ni++;
						chp2 = strchr(chp1, '\t');
					}
					strcpy(strTemp, chp1);
					StringAddTail(vSampleName+ni, strTemp);
				}

				vGroupId = NULL;
				vGroupId = CreateIntMatrix(1, nColNum);
				if(vGroupId == NULL)
				{
					printf("Error: TileMapv2_TXT2BAR, cannot create memory for group id!\n");
					exit(EXIT_FAILURE);
				}

				for(ni=0; ni<nColNum; ni++)
					vGroupId->pMatElement[ni] = 1;
			}
			else if(nAnnotLine == 1)
			{
				if(nFieldNum > 2)
				{
					/* read group id */
					chp1 = strLine;
					for(ni=0; ni<2; ni++)
					{
						chp2 = strchr(chp1, '\t');
						chp1 = chp2+1;
					}
					ni  = 0;
					chp2 = strchr(chp1, '\t');
					while(chp2 != NULL)
					{
						*chp2 = '\0';
						StrTrimRight(chp1);
						vGroupId->pMatElement[ni] = atoi(chp1);
						
						chp1 = chp2+1;
						ni++;
						chp2 = strchr(chp1, '\t');
					}
					vGroupId->pMatElement[ni] = atoi(chp1);
				}
			}
			
			nAnnotLine++;
			continue;
		}

		chp1 = strLine;
		chp2 = strchr(chp1, '\t');
		if(chp2 == NULL)
		{
			printf("Error: TileMapv2_TXT2BAR, wrong input file format!\n");
			exit(EXIT_FAILURE);
		}
		*chp2 = '\0';
		StrTrimRight(chp1);
		strcpy(strChr, chp1);

		if(strcmp(strChr, strLastChr) != 0)
		{
			if(strcmp(strLastChr, "") != 0)
			{
				sprintf(strTemp, "%s\t%d\t", strLastChr, nProbeNum);
				StringAddTail(&pSeqInfo, strTemp);
			}

			nSeqNum++;
			strcpy(strLastChr, strChr);
			nProbeNum = 1;
		}
		else
		{
			nProbeNum++;
		}
	}

	if(strcmp(strLastChr, "") != 0)
	{
		sprintf(strTemp, "%s\t%d\t", strLastChr, nProbeNum);
		StringAddTail(&pSeqInfo, strTemp);
	}

	fclose(fpIn);

	/* load head */
	if(nAnnotLine <= 0)
	{
		printf("Error: TileMapv2_TXT2BAR, sample names must be specified in the first line!\n");
		exit(EXIT_FAILURE);
	}

	if(nSeqNum <= 0)
	{
		printf("Warning: empty sequences!\n");
		return PROC_FAILURE;
	}

	/* create BAR object */
	pBARData = Affy_BARData_Create();
	if(pBARData == NULL)
	{
		printf("Error: TileMapv2_TXT2BAR, BARData object was not created!\n");
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
		printf("Error: TileMapv2_TXT2BAR, cannot allocate memory for field type.\n");
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
			printf("Error: TileMapv2_TXT2BAR, cannot allocate memory for loading bpmap coordinate data!\n");
			exit(EXIT_FAILURE);
		}
		if(pBARData->vSeqData[ni]->nDataNum > 0)
		{
			for(nj=0; nj<pBARData->vSeqData[ni]->nColNum; nj++)
			{
				pBARData->vSeqData[ni]->vData[nj] = CreateDoubleMatrix(1, pBARData->vSeqData[ni]->nDataNum );
				if(pBARData->vSeqData[ni]->vData[nj] == NULL)
				{
					printf("Error: TileMapv2_TXT2BAR, cannot allocate memory for loading bpmap coordinate data!\n");
					exit(EXIT_FAILURE);
				}
			}
		}
	}
	
	/* load data */
	fpIn = NULL;
	fpIn = fopen(strTXTFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileMapv2_TXT2BAR, cannot open source file!\n");
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
			printf("Error: TileMapv2_TXT2BAR, wrong input file format!\n");
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
					printf("Error: TileMapv2_TXT2BAR, cannot read input file correctly!");
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
				printf("Error: TileMapv2_TXT2BAR, input file format error, column number inconsistent!");
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
			printf("Error: TileMapv2_TXT2BAR, input file format error, column number inconsistent!");
			exit(EXIT_FAILURE);
		}
		
		StrTrimLeft(chp1);
		StrTrimRight(chp1);
		pBARData->vSeqData[nSeqNum]->vData[nk]->pMatElement[nProbeNum] = atof(chp1);
		nk++;

		if(nFieldNum == 2)
		{
			if(nk != (pBARData->vSeqData[nSeqNum]->nColNum-1))
			{
				printf("Error: TileMapv2_TXT2BAR, input file format error, column number inconsistent!");
				exit(EXIT_FAILURE);
			}
			else
			{
				pBARData->vSeqData[nSeqNum]->vData[nk]->pMatElement[nProbeNum] = 1.0;
			}
		}
		else
		{
			if(nk != pBARData->vSeqData[nSeqNum]->nColNum)
			{
				printf("Error: TileMapv2_TXT2BAR, input file format error, column number inconsistent!");
				exit(EXIT_FAILURE);
			}
		}
	}

	
	/* close file */
	fclose(fpIn);
	if(nProbeNum != (pBARData->vSeqData[nSeqNum]->nDataNum-1))
	{
		printf("Error: TileMapv2_TXT2BAR, cannot read input file correctly!");
		exit(EXIT_FAILURE);
	}

	if(nSeqNum != pBARData->nSeqNum-1)
	{
		printf("Error: TileMapv2_TXT2BAR, cannot read input file correctly!");
		exit(EXIT_FAILURE);
	}
	
	/* write to BAR file */
	pCol = NULL;
	pCol = CreateIntMatrix(1,pBARData->nColNum);
	if(pCol == NULL)
	{
		printf("Error: TileMapv2_TXT2BAR, cannot create output column information!\n");
		exit(EXIT_FAILURE);
	}
	
	pCol->pMatElement[0] = 1;
	for(ni=0; ni<nColNum; ni++)
	{
		pCol->pMatElement[1+ni] = 1;
		sprintf(strFileName, "%s%s.bar", strExportFolder, vSampleName[ni]->m_pString);
		Affy_SaveBAR_Columns_Fast(strFileName, pBARData, pCol);
		pCol->pMatElement[1+ni] = 0;
	}
	DestroyIntMatrix(pCol);

	/* check group id */
	nGroupNum = 0;
	for(ni=0; ni<nColNum; ni++)
	{
		if(vGroupId->pMatElement[ni] > nGroupNum)
			nGroupNum = vGroupId->pMatElement[ni];
	}

	if(nGroupNum > 0)
	{
		nGroupIdWrong = 0;
		pGroupOK = NULL;
		pGroupOK = CreateIntMatrix(1, nGroupNum);
		if(pGroupOK == NULL)
		{
			printf("Error: TileMapv2_TXT2BAR, cannot create memory for group id check!\n");
			exit(EXIT_FAILURE);
		}

		for(ni=0; ni<nColNum; ni++)
		{
			if(vGroupId->pMatElement[ni] > 0)
			{
				nk = vGroupId->pMatElement[ni]-1;
				pGroupOK->pMatElement[nk] = 1;
			}
		}

		for(ni=0; ni<nGroupNum; ni++)
		{
			if(pGroupOK->pMatElement[ni] == 0)
			{
				nGroupIdWrong = 1;
				break;
			}
		}

		if(nGroupIdWrong == 1)
		{
			printf("Warning: TileMapv2_TXT2BAR, group id setting was wrong, all group ids are now set to 1!\n");
			for(ni=0; ni<nColNum; ni++)
			{
				vGroupId->pMatElement[ni] = 1;
			}
			nGroupNum = 1;
		}

		DestroyIntMatrix(pGroupOK);
	}

	/* write to CGW file */
	sprintf(strFileName, "%s%s", strExportFolder, strProjectTitle);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: TileMapv2_TXT2BAR, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "[item1]\n");
	fprintf(fpOut, "type=bartilingexp\n");
	fprintf(fpOut, "item_name=%s\n", strProjectTitle);
	fprintf(fpOut, "bar_folder=%s\n", strExportFolder);
	fprintf(fpOut, "lib_num=1\n");
	fprintf(fpOut, "sample_num=%d\n", nColNum);
	fprintf(fpOut, "group_num=%d\n", nGroupNum);
	for(ni=0; ni<nColNum; ni++)
	{
		fprintf(fpOut, "sample_alias_%d=%s\n", ni+1, vSampleName[ni]->m_pString);
		fprintf(fpOut, "sample_group_%d=%d\n", ni+1, vGroupId->pMatElement[ni]);
		fprintf(fpOut, "bar_file_%d_1=%s.bar\n", ni+1, vSampleName[ni]->m_pString);
	}
	fclose(fpOut);
	
	/* clear memeory */
	for(ni=0; ni<nColNum; ni++)
	{
		DeleteString(vSampleName[ni]);
		vSampleName[ni] = NULL;
	}
	free(vSampleName);
	
	DestroyIntMatrix(vGroupId);
	vGroupId = NULL;

	DeleteString(pSeqInfo);
	pSeqInfo = NULL;

	Affy_BARData_Destroy(&pBARData);
		
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileMapv2_TXT_QuantileNormalization()                                  */
/*  Quantile normalization for a tab-delimited txt file.                   */
/*  Return PROC_SUCCESS if success.                                        */
/*  nTransform: 0=identity (default); 1=log2.                              */
/* ----------------------------------------------------------------------- */ 
int TileMapv2_TXT_QuantileNormalization(char strDataFile[], char strOutFile[],
			int nSkipColNum, int nTransform, double dTruncLow, double dTruncHigh)
{
	/* define */
	struct tagString **vChr;
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
	int nArrayNum = 0;
	int nProbeNum = 0;
	
	/* check */
	fpData = NULL;
	fpData = fopen(strDataFile, "r");
	if(fpData == NULL)
	{
		printf("Error: TileMapv2_TXT_QuantileNormalization, cannot find input file!\n");
		exit(EXIT_FAILURE);
	}

	nArrayNum = 0;
	nProbeNum = 0;
	while(fgets(strLongLine, LONG_LINE_LENGTH, fpData) != NULL)
	{
		StrTrimLeft(strLongLine);
		StrTrimRight(strLongLine);
		if(strLongLine[0] == '\0')
			continue;
		if(strLongLine[0] == '#')
			continue;

		if(nProbeNum == 0)
		{
			
			chSep = strLongLine;
			chSep2 = strchr(chSep, '\t');
			while(chSep2 != NULL)
			{
				nArrayNum++;
				chSep = chSep2+1;
				chSep2 = strchr(chSep, '\t');
			}
			nArrayNum++;
		}

		nProbeNum++;
	}

	fclose(fpData);

	nArrayNum -= nSkipColNum;

	if( (nArrayNum <= 0) || (nProbeNum <= 0) )
	{
		printf("Warning: TileMapv2_TXT_QuantileNormalization, array or probe number <= 0!\n");
		return PROC_SUCCESS;
	}

	/* init */
	vChr = NULL;
	vChr = (struct tagString **)calloc(nProbeNum, sizeof(struct tagString *));
	if(vChr == NULL)
	{
		printf("Error: TileMapv2_TXT_QuantileNormalization, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}

	vArray = NULL;
	vArray = (struct DOUBLEMATRIX **)calloc(nArrayNum, sizeof(struct DOUBLEMATRIX *));
	if(vArray == NULL)
	{
		printf("Error: TileMapv2_TXT_QuantileNormalization, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nArrayNum; ni++)
	{
		vArray[ni] = CreateDoubleMatrix(1, nProbeNum);
		if(vArray[ni] == NULL)
		{
			printf("Error: TileMapv2_TXT_QuantileNormalization, cannot create data space!\n");
			exit(EXIT_FAILURE);
		}
	}

	vSortArray = NULL;
	vSortArray = (struct DOUBLEMATRIX **)calloc(nArrayNum, sizeof(struct DOUBLEMATRIX *));
	if(vSortArray == NULL)
	{
		printf("Error: TileMapv2_TXT_QuantileNormalization, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}
	vSortIndex = NULL;
	vSortIndex = (struct LONGMATRIX **)calloc(nArrayNum, sizeof(struct LONGMATRIX *));
	if(vSortIndex == NULL)
	{
		printf("Error: TileMapv2_TXT_QuantileNormalization, cannot create data space!\n");
		exit(EXIT_FAILURE);
	}

	/* data */
	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: TileMapv2_TXT_QuantileNormalization, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fpData = NULL;
	fpData = fopen(strDataFile, "r");
	if(fpData == NULL)
	{
		printf("Error: TileMapv2_TXT_QuantileNormalization, cannot open data file!\n");
		exit(EXIT_FAILURE);
	}
	
	nj = 0;
	while(fgets(strLongLine, LONG_LINE_LENGTH, fpData) != NULL)
	{
		StrTrimLeft(strLongLine);
		StrTrimRight(strLongLine);
		if(strLongLine[0] == '\0')
			continue;
		if(strLongLine[0] == '#')
		{
			fprintf(fpOut, "%s\n", strLongLine);
			continue;
		}

		if(nj >= nProbeNum)
		{
			printf("Error: TileMapv2_TXT_QuantileNormalization, probe number not match!\n");
			exit(EXIT_FAILURE);
		}

		chSep = strLongLine;
		if(nSkipColNum > 0)
		{
			for(ni=0; ni<nSkipColNum; ni++)
			{
				chSep2 = strchr(chSep, '\t');
				chSep = chSep2+1;
			}
			*chSep2 = '\0';
			StringAddTail((vChr+nj), strLongLine);
		}
		
		ni = 0;
		chSep2 = strchr(chSep, '\t');
		while(chSep2 != NULL)
		{
			if( ni >= nArrayNum )
			{
				printf("Error: TileMapv2_TXT_QuantileNormalization, array number not match!\n");
				exit(EXIT_FAILURE);
			}

			/* middle number */
			*chSep2 = '\0';
			if(chSep == chSep2)
			{
				dValue = 0.0;
			}
			else
			{
				dValue = atof(chSep);
			}

			if(dValue < dTruncLow)
				dValue = dTruncLow;
			if(dValue > dTruncHigh)
				dValue = dTruncHigh;
			if(nTransform == 1)
				dValue = log(dValue)/dLog2;
		
			vArray[ni]->pMatElement[nj] = dValue;
			ni++;

			/* get next */
			chSep = chSep2+1;
			chSep2 = strchr(chSep, '\t');
		}

		if( ni >= nArrayNum )
		{
			printf("Error: TileMapv2_TXT_QuantileNormalization, array number not match!\n");
			exit(EXIT_FAILURE);
		}

		/* middle number */
		dValue = atof(chSep);

		if(dValue < dTruncLow)
			dValue = dTruncLow;
		if(dValue > dTruncHigh)
			dValue = dTruncHigh;
		if(nTransform == 1)
			dValue = log(dValue)/dLog2;
	
		vArray[ni]->pMatElement[nj] = dValue;
		ni++;

		if(ni!=nArrayNum)
		{
			printf("Error: TileMapv2_TXT_QuantileNormalization, array number not match!\n");
			exit(EXIT_FAILURE);
		}

		/* get next line */
		nj++;
	}

	fclose(fpData);

	if(nj != nProbeNum)
	{
		printf("Error: TileMapv2_TXT_QuantileNormalization, probe number not match!\n");
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
		fprintf(fpOut, "%s", vChr[nj]->m_pString);
		for(ni=0; ni<nArrayNum; ni++)
		{
			fprintf(fpOut, "\t%f", vArray[ni]->pMatElement[nj]);
		}
		fprintf(fpOut, "\n");
	}


	/* release memory */
	fclose(fpOut);
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
/*  TileProbe_Build_Main()                                                 */
/*  Build a probe background model for microarray CEL files.               */
/* ----------------------------------------------------------------------- */ 
int TileProbe_Build_Main(char strInputPath[], char strOutputPath[], 
						 int nTakeLog, double dNormLowerBound, int nShrink, int nTest)
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	FILE *fpOut2;

	/* arrays */
	int nArrayNum = 0;
	struct tagCELData **vCELData = NULL;
	int nCELNumberCells = 0;
	int nCELTotalX = 0;
	int nCELTotalY = 0;
	struct DOUBLEMATRIX *pSortMean;
	struct DOUBLEMATRIX *pArray;
	struct LONGMATRIX *pSortId;
	struct DOUBLEMATRIX *pProbeMedian;
	struct DOUBLEMATRIX *pProbeInterR;
	struct DOUBLEMATRIX *pProbeData;
	struct DOUBLEMATRIX *pProbeDataSort;
	struct LONGMATRIX *pProbeDataSortId;
	struct DOUBLEMATRIX *pOldIntensity;
	struct DOUBLEMATRIX *pOldSD;
	
	/* others */
	int ni,nj;
	double dLog2 = log(2.0);
	int nId;
	int nQ50,nQ25,nQ05;
	char strOutFile[MED_LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];

	/* for IQR shrinking */
	double n1,n2,n3,nt,d1,d2,varv,varu,v0,u0,alpha,beta,lambda,dB;
	double dpi = 4.0*atan(1.0);
	double dQm,dQv,dQm2,dQmin,dQmax;

	/* get array number */
	nArrayNum = 0;
	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileProbe_Build_Main, cannot open the input array list!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		nArrayNum++;
	}
	
	fclose(fpIn);
	
	printf("########################################\n");
	printf("# Loading Raw Data                     #\n");
	printf("########################################\n");
	printf(" Array number = %d\n", nArrayNum);

	if(nArrayNum == 0)
	{
		return PROC_SUCCESS;
	}

	/* load arrays */
	vCELData = NULL;
	vCELData = (struct tagCELData **)calloc(nArrayNum, sizeof(struct tagCELData *));
	if(vCELData == NULL)
	{
		printf("Error: TileProbe_Build_Main, cannot allocate memory for the array list!\n");
		exit(EXIT_FAILURE);
	}

	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileProbe_Build_Main, cannot open the input array list!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		vCELData[ni] = TileMapv2_LoadCEL(strLine);
		if(vCELData[ni] == NULL)
		{
			printf("Error: TileProbe_Build_Main, cannot load *.CEL file!\n");
			exit(EXIT_FAILURE);
		}

		if(ni == 0)
		{
			nCELNumberCells = vCELData[ni]->nNumberCells;
			nCELTotalX = vCELData[ni]->nCols;
			nCELTotalY = vCELData[ni]->nRows;
			if( nCELNumberCells != nCELTotalX*nCELTotalY)
			{
				printf("Error: TileProbe_Build_Main, CEL file Cols*Rows != NumberCells!\n");
				exit(EXIT_FAILURE);
			}
		}
		else
		{
			if( (nCELNumberCells != vCELData[ni]->nNumberCells) || (nCELTotalX != vCELData[ni]->nCols)
				|| (nCELTotalY != vCELData[ni]->nRows) )
			{
				printf("Error: TileProbe_Build_Main, CEL file dimensions do not match!\n");
				exit(EXIT_FAILURE);
			}
		}

		if(nTakeLog == 1)
		{
			for(nj=0; nj<vCELData[ni]->nNumberCells; nj++)
			{
				/* truncate */
				if(vCELData[ni]->pIntensity->pMatElement[nj] < dNormLowerBound)
					vCELData[ni]->pIntensity->pMatElement[nj] = dNormLowerBound;

				/* log transformation */
				vCELData[ni]->pIntensity->pMatElement[nj] = log(vCELData[ni]->pIntensity->pMatElement[nj])/dLog2;
			}
		}
		else
		{
			for(nj=0; nj<vCELData[ni]->nNumberCells; nj++)
			{
				/* truncate */
				if(vCELData[ni]->pIntensity->pMatElement[nj] < dNormLowerBound)
					vCELData[ni]->pIntensity->pMatElement[nj] = dNormLowerBound;
			}
		}

		ni++;
	}
	
	fclose(fpIn);

	if(ni != nArrayNum)
	{
		printf("Error: TileProbe_Build_Main, array number inconsistent!\n");
		exit(EXIT_FAILURE);
	}

	printf("########################################\n");
	printf("# Quantile Normalization               #\n");
	printf("########################################\n");
	
	pSortMean = NULL;
	pSortMean = CreateDoubleMatrix(1, nCELNumberCells);
	if(pSortMean == NULL)
	{
		printf("Error: TileProbe_Build_Main, cannot allocate memory for quantile normalization!\n");
		exit(EXIT_FAILURE);
	}
				
	for(ni=0; ni<nArrayNum; ni++)
	{
		pArray = NULL;
		pSortId = NULL;

		/* export sorted values */
		printf("  Sorting array %d...\n", ni);
		DMSORTMERGEA_0(vCELData[ni]->pIntensity, &pArray, &pSortId);

		/* compute percentiles */
		for(nj=0; nj<nCELNumberCells; nj++)
		{
			pSortMean->pMatElement[nj] += pArray->pMatElement[nj];
			nId = pSortId->pMatElement[nj];
			vCELData[ni]->pIntensity->pMatElement[nId] = nj;
		}

		DestroyDoubleMatrix(pArray);
		DestroyLongMatrix(pSortId);
	}

	for(nj=0; nj<nCELNumberCells; nj++)
	{
		pSortMean->pMatElement[nj] /= nArrayNum;
	}

	for(ni=0; ni<nArrayNum; ni++)
	{
		for(nj=0; nj<nCELNumberCells; nj++)
		{
			nId = (long)(vCELData[ni]->pIntensity->pMatElement[nj]);
			vCELData[ni]->pIntensity->pMatElement[nj] = pSortMean->pMatElement[nId];
		}
	}

	/* build probe model */
	printf("########################################\n");
	printf("# Build Probe Model                    #\n");
	printf("########################################\n");

	if(nTest == 1)
	{
		sprintf(strOutFile, "%s.testq", strOutputPath);
		fpOut = NULL;
		fpOut = fopen(strOutFile, "w");
		if(fpOut == NULL)
		{
			printf("Error: TileProbe_Build_Main, cannot open the output argument file!\n");
			exit(EXIT_FAILURE);
		}
		sprintf(strOutFile, "%s.testprobe", strOutputPath);
		fpOut2 = NULL;
		fpOut2 = fopen(strOutFile, "w");
		if(fpOut2 == NULL)
		{
			printf("Error: TileProbe_Build_Main, cannot open the output argument file!\n");
			exit(EXIT_FAILURE);
		}
	}

	pProbeMedian = NULL;
	pProbeMedian = CreateDoubleMatrix(1,nCELNumberCells);
	if(pProbeMedian == NULL)
	{
		printf("Error: TileProbe_Build_Main, cannot allocate memory for probe median!\n");
		exit(EXIT_FAILURE);
	}

	pProbeInterR = NULL;
	pProbeInterR = CreateDoubleMatrix(1,nCELNumberCells);
	if(pProbeInterR == NULL)
	{
		printf("Error: TileProbe_Build_Main, cannot allocate memory for probe interquatile range!\n");
		exit(EXIT_FAILURE);
	}

	pProbeData = NULL;
	pProbeData = CreateDoubleMatrix(1,nArrayNum);
	if(pProbeData == NULL)
	{
		printf("Error: TileProbe_Build_Main, cannot allocate memory for modeling probe effect!\n");
		exit(EXIT_FAILURE);
	}

	nQ50 = (int)(nArrayNum/2.0);
	nQ25 = (int)(nArrayNum/4.0);
	nQ05 = (int)(nArrayNum/20.0);

	for(nj=0; nj<nCELNumberCells; nj++)
	{
		if(nj%100000 == 0)
		{
			printf("%d ...\n", nj);
		}

		pProbeDataSort = NULL;
		pProbeDataSortId = NULL;

		for(ni=0; ni<nArrayNum; ni++)
		{
			pProbeData->pMatElement[ni] = vCELData[ni]->pIntensity->pMatElement[nj];
		}

		DMSORTMERGEA_0(pProbeData, &pProbeDataSort, &pProbeDataSortId);
		pProbeMedian->pMatElement[nj] = pProbeDataSort->pMatElement[nQ50];
		pProbeInterR->pMatElement[nj] = pProbeDataSort->pMatElement[nQ50]-pProbeDataSort->pMatElement[nQ25];
		if(pProbeInterR->pMatElement[nj] < 0.0)
		{
			pProbeInterR->pMatElement[nj] = pProbeInterR->pMatElement[nj];
		}

		if(nTest == 1)
		{
			fprintf(fpOut, "%e\t%e\n", pProbeInterR->pMatElement[nj], pProbeDataSort->pMatElement[nQ50]-pProbeDataSort->pMatElement[nQ05]); 
			if(rand_u() < 0.01)
			{
				fprintf(fpOut2, "%f", pProbeData->pMatElement[0]);
				for(ni=1; ni<nArrayNum; ni++)
				{
					fprintf(fpOut2, "\t%f", pProbeData->pMatElement[ni]);
				}
				fprintf(fpOut2, "\n");
			}
		}

		DestroyDoubleMatrix(pProbeDataSort);
		DestroyLongMatrix(pProbeDataSortId);
	}

	/* IQR shrinkage */
	if(nShrink != 0)
	{
		n1 = nQ25;
		n2 = nQ50-nQ25;
		n3 = nArrayNum-nQ50;
		nt = nArrayNum;
		u0 = n1/(double)(nQ50);
		v0 = (double)(nQ50)/(nt+1);
		d1 = norminv(0.0,1.0,v0);
		d2 = norminv(0.0,1.0,u0*v0);
		varu = n1*n2/pow((n1+n2),2.0)/(n1+n2+1);
		varv = (n1+n2)*(n3+1)/pow((nt+1.0),2.0)/(nt+2.0);
		alpha = 1.0+(dpi*d1*exp(d1*d1)*varv-dpi*d2*exp(d2*d2)*(v0*v0*varu+u0*u0*varv))/(d1-d2);
		beta = 2*dpi*(v0*v0*exp(d2*d2)*varu+(exp(d1*d1)+u0*u0*exp(d2*d2)-2*u0*exp(d1*d1/2.0+d2*d2/2.0))*varv)/pow((d1-d2),2.0);
		lambda = beta/(alpha*alpha);
		
		dQmin = 1e10;
		dQmax = -1e10;
		dQm = 0.0;
		dQm2 = 0.0;
		dQv = 0.0;
		for(nj=0; nj<nCELNumberCells; nj++)
		{
			pProbeInterR->pMatElement[nj] /= alpha;

			if(pProbeInterR->pMatElement[nj] < dQmin)
				dQmin = pProbeInterR->pMatElement[nj];
			if(pProbeInterR->pMatElement[nj] > dQmax)
				dQmax = pProbeInterR->pMatElement[nj];

			dQm += pProbeInterR->pMatElement[nj];
			dQm2 += pProbeInterR->pMatElement[nj]*pProbeInterR->pMatElement[nj];
		}
		dQm /= nCELNumberCells;
		dQm2 /= nCELNumberCells;
		for(nj=0; nj<nCELNumberCells; nj++)
		{
			dQv += (pProbeInterR->pMatElement[nj]-dQm)*(pProbeInterR->pMatElement[nj]-dQm);
		}
		dQv /= (nCELNumberCells-1);

		
		dB = (lambda/(1.0+lambda))*(dQm2/dQv);

		if(dB>1.0)
			dB=1.0;
		else if(dB<0.0)
			dB=0.0;

		if(nShrink == -1)
			dB = 0.0;

		for(nj=0; nj<nCELNumberCells; nj++)
		{
			pProbeInterR->pMatElement[nj] = pProbeInterR->pMatElement[nj]*(1.0-dB)+dQm*dB;
		}
		
		printf("min IQR = %f; max IQR = %f\n", dQmin, dQmax);
		printf("B = %f\n", dB);
	}

	if(nTest == 1)
	{
		fclose(fpOut);
		fclose(fpOut2);
	}

	/* save file */
	printf("########################################\n");
	printf("# Export Probe Model                   #\n");
	printf("########################################\n");
	pOldIntensity = NULL;
	pOldSD = NULL;
	pOldIntensity = vCELData[0]->pIntensity;
	pOldSD = vCELData[0]->pSD;
	vCELData[0]->pIntensity = pProbeMedian;
	vCELData[0]->pSD = pProbeInterR;

	sprintf(strOutFile, "%s.prbm", strOutputPath);
	Affy_SaveCELv4(strOutFile, vCELData[0]);

	vCELData[0]->pIntensity = pSortMean;
	sprintf(strOutFile, "%s.quan", strOutputPath);
	Affy_SaveCELv4(strOutFile, vCELData[0]);

	vCELData[0]->pIntensity = pOldIntensity;
	vCELData[0]->pSD = pOldSD;


	/* release memory */
	for(ni=0; ni<nArrayNum; ni++)
	{
		Affy_CELData_Destroy(vCELData+ni);
	}
	free(vCELData);

	DestroyDoubleMatrix(pSortMean);
	DestroyDoubleMatrix(pProbeMedian);
	DestroyDoubleMatrix(pProbeInterR);
	DestroyDoubleMatrix(pProbeData);

	printf("Done!\n");

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_Norm_Main()                                                  */
/*  Adjust probe intensities based on the probe background model to remove */
/*  probe effects.                                                         */
/*  nInputType: 0=single input cel file; 1=a file of array list            */
/*  strOutputPath: folder of output                                        */
/*  dB shrinkage factor for Q50-Q25.                                       */
/* ----------------------------------------------------------------------- */ 
int TileProbe_Norm_Main(char strInputPath[], int nInputType,
						char strOutputPath[], char strModelPath[],
						int nTakeLog, double dNormLowerBound,
						double dB, int nLogAfterNorm)
{
	/* define */
	FILE *fpIn;
	struct tagCELData *pProbeModel;
	struct tagCELData *pProbeQuan;
	struct tagString **vCELFile;

	/* arrays */
	int nArrayNum = 0;
	struct tagCELData *pCELData = NULL;
	struct DOUBLEMATRIX *pArray;
	struct LONGMATRIX *pSortId;

	/* others */
	int ni,nj;
	double dLog2 = log(2.0);
	int nId;
	char strOutFile[MED_LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	char strTemp[MED_LINE_LENGTH];
	char strOutFolder[MED_LINE_LENGTH];
	double dSDM;

	/* init */
	strcpy(strOutFolder, strOutputPath);
	AdjustDirectoryPath(strOutFolder);


	/* get array number */
	if(nInputType == 1)
	{
		nArrayNum = 0;
		fpIn = NULL;
		fpIn = fopen(strInputPath, "r");
		if(fpIn == NULL)
		{
			printf("Error: TileProbe_Norm_Main, cannot open the input array list!\n");
			exit(EXIT_FAILURE);
		}

		while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;

			nArrayNum++;
		}
		
		fclose(fpIn);
	}
	else
	{
		nArrayNum = 1;
	}

	vCELFile = NULL;
	vCELFile = (struct tagString **)calloc(nArrayNum, sizeof(struct tagString *));
	if(vCELFile == NULL)
	{
		printf("Error: TileProbe_Norm_Main, cannot allocate memory for the array list!\n");
		exit(EXIT_FAILURE);
	}

	if(nInputType == 1)
	{
		fpIn = NULL;
		fpIn = fopen(strInputPath, "r");
		if(fpIn == NULL)
		{
			printf("Error: TileProbe_Norm_Main, cannot open the input array list!\n");
			exit(EXIT_FAILURE);
		}

		ni = 0;
		while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;

			StringAddTail(vCELFile+ni, strLine);
			ni++;
		}
		
		fclose(fpIn);

		if(ni != nArrayNum)
		{
			printf("Error: TileProbe_Norm_Main, array number inconsistent!\n");
			exit(EXIT_FAILURE);
		}
	}
	else
	{
		StringAddTail(vCELFile+0, strInputPath);
	}

	
	/* load model */
	sprintf(strLine, "%s.prbm", strModelPath);
	pProbeModel = NULL;
	pProbeModel = TileMapv2_LoadCEL(strLine);
	if(pProbeModel == NULL)
	{
		printf("Error: TileProbe_Build_Main, cannot load probe model!\n");
		exit(EXIT_FAILURE);
	}

	sprintf(strLine, "%s.quan", strModelPath);
	pProbeQuan = NULL;
	pProbeQuan = TileMapv2_LoadCEL(strLine);
	if(pProbeQuan == NULL)
	{
		printf("Error: TileProbe_Build_Main, cannot load probe model!\n");
		exit(EXIT_FAILURE);
	}

	/* shrink variance */
	if(dB > 0.0)
	{
		pArray = NULL;
		pSortId = NULL;

		DMSORTMERGEA_0(pProbeModel->pSD, &pArray, &pSortId);
		ni = (int)(pProbeModel->nNumberCells*0.5);
		dSDM = pArray->pMatElement[ni];
		for(ni=0; ni<pProbeModel->nNumberCells; ni++)
		{
			pProbeModel->pSD->pMatElement[ni] = (1-dB)*pProbeModel->pSD->pMatElement[ni]+dB*dSDM;
		}

		DestroyDoubleMatrix(pArray);
		DestroyLongMatrix(pSortId);
	}

	printf("########################################\n");
	printf("# Processing                           #\n");
	printf("########################################\n");
	printf(" Array number = %d\n", nArrayNum);

	if(nArrayNum == 0)
	{
		return PROC_SUCCESS;
	}


	for(ni=0; ni<nArrayNum; ni++)
	{
		/* load data */
		printf("Processing %s ...\n", vCELFile[ni]->m_pString);
		pCELData = NULL;
		pCELData = TileMapv2_LoadCEL(vCELFile[ni]->m_pString);
		if(pCELData == NULL)
		{
			printf("Error: TileProbe_Norm_Main, cannot load *.CEL file!\n");
			exit(EXIT_FAILURE);
		}

		if( (pProbeModel->nNumberCells != pCELData->nNumberCells) || (pProbeModel->nCols != pCELData->nCols)
			|| (pProbeModel->nRows != pCELData->nRows) )
		{
			printf("Error: TileProbe_Norm_Main, CEL file dimension (%s) does not match with the model!\n", vCELFile[ni]->m_pString);
			exit(EXIT_FAILURE);
		}

		/* pre-processing */
		if(nLogAfterNorm == 0)
		{
			if(nTakeLog == 1)
			{
				for(nj=0; nj<pCELData->nNumberCells; nj++)
				{
					/* truncate */
					if(pCELData->pIntensity->pMatElement[nj] < dNormLowerBound)
						pCELData->pIntensity->pMatElement[nj] = dNormLowerBound;

					/* log transformation */
					pCELData->pIntensity->pMatElement[nj] = log(pCELData->pIntensity->pMatElement[nj])/dLog2;
				}
			}
			else
			{
				for(nj=0; nj<pCELData->nNumberCells; nj++)
				{
					/* truncate */
					if(pCELData->pIntensity->pMatElement[nj] < dNormLowerBound)
						pCELData->pIntensity->pMatElement[nj] = dNormLowerBound;
				}
			}
		}

		/* quantile normalization */
		pArray = NULL;
		pSortId = NULL;

		/* export sorted values */
		DMSORTMERGEA_0(pCELData->pIntensity, &pArray, &pSortId);

		/* reset cel values */
		for(nj=0; nj<pCELData->nNumberCells; nj++)
		{
			nId = pSortId->pMatElement[nj];
			pCELData->pIntensity->pMatElement[nId] = pProbeQuan->pIntensity->pMatElement[nj];
			if(nLogAfterNorm != 0)
			{
				if(pCELData->pIntensity->pMatElement[nId] < dNormLowerBound)
					pCELData->pIntensity->pMatElement[nId] = dNormLowerBound;

				if(nTakeLog == 1)
				{
					pCELData->pIntensity->pMatElement[nId] = log(pCELData->pIntensity->pMatElement[nId])/dLog2;
				}
			}
			pCELData->pIntensity->pMatElement[nId] = (pCELData->pIntensity->pMatElement[nId]-pProbeModel->pIntensity->pMatElement[nId])/(pProbeModel->pSD->pMatElement[nId]+1e-6);
		}

		/* export results */
		GetFileName(vCELFile[ni]->m_pString, strTemp);
		sprintf(strOutFile, "%s%s.tpnorm.CEL", strOutFolder, strTemp);
		Affy_SaveCELv4(strOutFile, pCELData);

		DestroyDoubleMatrix(pArray);
		DestroyLongMatrix(pSortId);
		Affy_CELData_Destroy(&pCELData);
	}
		
	/* release memory */
	for(ni=0; ni<nArrayNum; ni++)
	{
		DeleteString(vCELFile[ni]);
		vCELFile[ni] = NULL;
	}
	free(vCELFile);

	Affy_CELData_Destroy(&pProbeModel);
	Affy_CELData_Destroy(&pProbeQuan);
	
	printf("Done!\n");

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_Peak_Main()                                                  */
/*  Find peaks.                                                            */
/* ----------------------------------------------------------------------- */ 
int TileProbe_Peak_Main(char strInputPath[])
{
	/* define */
	FILE *fpIn;
	char strGrpName[255];
	int nW = 300;
	int nMaxGap = 300;
	int nMinProbe = 10;
	int nUseVar = 0;
	int nIPNum = 0;
	int nCTNum = 0;
	double dCut = 3.0;
	char strOutputPath[MED_LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	char *chp;
	int ni,nj,nk;
	int nTotalProbeNum = 0;

	struct tagBARData *vIP = NULL;
	struct tagBARData *vCT = NULL;
	struct tagBARData *pArray = NULL;
	struct INTMATRIX *pFieldType = NULL;
	struct DOUBLEMATRIX **vDataVec = NULL;

	/* load parameter */
	strcpy(strGrpName, "");
	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileProbe_Peak_Main, cannot open parameter file!\n");
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

		if(strstr(strLine, "Output_Path") == strLine)
		{
			chp = strchr(strLine, '=');
			if(chp != NULL)
			{
				chp++;
				strcpy(strOutputPath, chp);
			}
			else
			{
				strcpy(strOutputPath, "tileprobepeak1");
			}
		}
		else if(strstr(strLine, "GenomeGrp") == strLine)
		{
			chp = strchr(strLine, '=');
			if(chp != NULL)
			{
				chp++;
				strcpy(strGrpName, chp);
				StrTrimLeft(strGrpName);
				StrTrimRight(strGrpName);
			}
		}
		else if(strstr(strLine, "BandWidth") == strLine)
		{
			chp = strchr(strLine, '=');
			if(chp != NULL)
			{
				chp++;
				nW = atoi(chp);
			}
		}
		else if(strstr(strLine, "MaxGap") == strLine)
		{
			chp = strchr(strLine, '=');
			if(chp != NULL)
			{
				chp++;
				nMaxGap = atoi(chp);
			}
		}
		else if(strstr(strLine, "MinProbe") == strLine)
		{
			chp = strchr(strLine, '=');
			if(chp != NULL)
			{
				chp++;
				nMinProbe = atoi(chp);
			}
		}
		else if(strstr(strLine, "Var") == strLine)
		{
			chp = strchr(strLine, '=');
			if(chp != NULL)
			{
				chp++;
				nUseVar = atoi(chp);
			}
		}
		else if(strstr(strLine, "Cutoff") == strLine)
		{
			chp = strchr(strLine, '=');
			if(chp != NULL)
			{
				chp++;
				dCut = atof(chp);
			}
		}
		else if(strstr(strLine, "IP") == strLine)
		{
			chp = strchr(strLine, '=');
			if(chp != NULL)
			{
				chp++;
				nIPNum = atoi(chp);
			}

			if(nIPNum <= 0)
			{
				printf("Error: TileProbe_Peak_Main, at least one IP sample need to be provided!\n");
				exit(EXIT_FAILURE);
			}

			ni = 0;
			while((ni<nIPNum) && (fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL))
			{
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] == '\0')
					continue;
				if(strLine[0] == '#')
					continue;

				pArray = NULL;
				pArray = Affy_LoadBAR_Fast(strLine);
				if(pArray == NULL)
				{
					printf("Error: TileProbe_Peak_Main, cannot load raw data!\n");
					exit(EXIT_FAILURE);
				}

				/* if the first array, prepare enough space for nIPNum arrays + 
					1 positions + 1 group trimed mean + 1 group variance + 1 test statistics */ 
				if(ni == 0)
				{
					vIP = pArray;
					vIP->nColNum = nIPNum+4;
					pFieldType = vIP->pFieldType;
					vIP->pFieldType = NULL;
					vIP->pFieldType = CreateIntMatrix(1,vIP->nColNum);
					if(vIP->pFieldType == NULL)
					{
						printf("Error: TileProbe_Peak_Main, cannot create memory for field types!\n");
						exit(EXIT_FAILURE);
					}
					for(nj=0; nj<vIP->nColNum; nj++)
						vIP->pFieldType->pMatElement[nj] = 1;
					vIP->pFieldType->pMatElement[0] = pFieldType->pMatElement[0];
					vIP->pFieldType->pMatElement[1] = pFieldType->pMatElement[1];
					DestroyIntMatrix(pFieldType);
					
					for(nj=0; nj<vIP->nSeqNum; nj++)
					{
						nTotalProbeNum += vIP->vSeqData[nj]->nDataNum;
						vIP->vSeqData[nj]->nColNum = vIP->nColNum;
						vDataVec = vIP->vSeqData[nj]->vData;
						vIP->vSeqData[nj]->vData = NULL;
						vIP->vSeqData[nj]->vData = (struct DOUBLEMATRIX **)calloc(vIP->vSeqData[nj]->nColNum, sizeof(struct DOUBLEMATRIX *));
						if(vIP->vSeqData[nj]->vData == NULL)
						{
							printf("Error: TileProbe_Peak_Main, cannot create memory for tracking intensity data!\n");
							exit(EXIT_FAILURE);
						}
						vIP->vSeqData[nj]->vData[0] = vDataVec[0];
						vIP->vSeqData[nj]->vData[1] = vDataVec[1];
						free(vDataVec);

						for(nk=nIPNum+1; nk<vIP->nColNum; nk++)
						{
							vIP->vSeqData[nj]->vData[nk] = CreateDoubleMatrix(1, vIP->vSeqData[nj]->nDataNum);
							if((vIP->vSeqData[nj]->vData[nk] == NULL) && (vIP->vSeqData[nj]->nDataNum > 0))
							{
								printf("Error: TileProbe_Peak_Main, cannot create memory for data processing!\n");
								exit(EXIT_FAILURE);
							}
						}
					}
				}
				/* if not the first array, transfer to the first array */
				else
				{
					if(pArray->nSeqNum != vIP->nSeqNum)
					{
						printf("Error: TileProbe_Peak_Main, array types do not match!\n");
						exit(EXIT_FAILURE);
					}

					for(nj=0; nj<vIP->nSeqNum; nj++)
					{
						if(vIP->vSeqData[nj]->nDataNum != pArray->vSeqData[nj]->nDataNum)
						{
							printf("Error: TileProbe_Peak_Main, array types do not match!\n");
							exit(EXIT_FAILURE);
						}
						vIP->vSeqData[nj]->vData[ni+1] = pArray->vSeqData[nj]->vData[1];
						pArray->vSeqData[nj]->vData[1] = NULL;
					}

					vIP->pFieldType->pMatElement[ni+1] = pArray->pFieldType->pMatElement[1];

					Affy_BARData_Destroy(&pArray);
				}


				ni++;
			}

			
			if(ni!=nIPNum)
			{
				printf("Error: TileProbe_Peak_Main, IP sample number inconsistent!\n");
				exit(EXIT_FAILURE);
			}
		}

		else if(strstr(strLine, "CT") == strLine)
		{
			chp = strchr(strLine, '=');
			if(chp != NULL)
			{
				chp++;
				nCTNum = atoi(chp);
			}

			ni = 0;
			while((ni<nCTNum) && (fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL))
			{
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] == '\0')
					continue;
				if(strLine[0] == '#')
					continue;

				pArray = NULL;
				pArray = Affy_LoadBAR_Fast(strLine);
				if(pArray == NULL)
				{
					printf("Error: TileProbe_Peak_Main, cannot load raw data!\n");
					exit(EXIT_FAILURE);
				}

				/* if the first array, prepare enough space for nIPNum arrays + 
					1 positions + 1 group trimed mean + 1 group variance + 1 test statistics */ 
				if(ni == 0)
				{
					vCT = pArray;
					vCT->nColNum = nCTNum+4;
					pFieldType = vCT->pFieldType;
					vCT->pFieldType = NULL;
					vCT->pFieldType = CreateIntMatrix(1,vCT->nColNum);
					if(vCT->pFieldType == NULL)
					{
						printf("Error: TileProbe_Peak_Main, cannot create memory for field types!\n");
						exit(EXIT_FAILURE);
					}
					for(nj=0; nj<vCT->nColNum; nj++)
						vCT->pFieldType->pMatElement[nj] = 1;
					vCT->pFieldType->pMatElement[0] = pFieldType->pMatElement[0];
					vCT->pFieldType->pMatElement[1] = pFieldType->pMatElement[1];
					DestroyIntMatrix(pFieldType);
					
					for(nj=0; nj<vCT->nSeqNum; nj++)
					{
						vCT->vSeqData[nj]->nColNum = vCT->nColNum;
						vDataVec = vCT->vSeqData[nj]->vData;
						vCT->vSeqData[nj]->vData = NULL;
						vCT->vSeqData[nj]->vData = (struct DOUBLEMATRIX **)calloc(vCT->vSeqData[nj]->nColNum, sizeof(struct DOUBLEMATRIX *));
						if(vCT->vSeqData[nj]->vData == NULL)
						{
							printf("Error: TileProbe_Peak_Main, cannot create memory for tracking intensity data!\n");
							exit(EXIT_FAILURE);
						}
						vCT->vSeqData[nj]->vData[0] = vDataVec[0];
						vCT->vSeqData[nj]->vData[1] = vDataVec[1];
						free(vDataVec);

						for(nk=nCTNum+1; nk<vCT->nColNum; nk++)
						{
							vCT->vSeqData[nj]->vData[nk] = CreateDoubleMatrix(1, vCT->vSeqData[nj]->nDataNum);
							if( (vCT->vSeqData[nj]->vData[nk] == NULL) && (vCT->vSeqData[nj]->nDataNum > 0) )
							{
								printf("Error: TileProbe_Peak_Main, cannot create memory for data processing!\n");
								exit(EXIT_FAILURE);
							}
						}
					}
				}
				/* if not the first array, transfer to the first array */
				else
				{
					if(pArray->nSeqNum != vCT->nSeqNum)
					{
						printf("Error: TileProbe_Peak_Main, array types do not match!\n");
						exit(EXIT_FAILURE);
					}

					for(nj=0; nj<vCT->nSeqNum; nj++)
					{
						if(vCT->vSeqData[nj]->nDataNum != pArray->vSeqData[nj]->nDataNum)
						{
							printf("Error: TileProbe_Peak_Main, array types do not match!\n");
							exit(EXIT_FAILURE);
						}
						vCT->vSeqData[nj]->vData[ni+1] = pArray->vSeqData[nj]->vData[1];
						pArray->vSeqData[nj]->vData[1] = NULL;
					}

					vCT->pFieldType->pMatElement[ni+1] = pArray->pFieldType->pMatElement[1];

					Affy_BARData_Destroy(&pArray);
				}


				ni++;
			}

			
			if(ni!=nCTNum)
			{
				printf("Error: TileProbe_Peak_Main, CT sample number inconsistent!\n");
				exit(EXIT_FAILURE);
			}
		}

	}

	fclose(fpIn);

	TileProbe_Peak(nIPNum, nCTNum, vIP, vCT,
				   nW, nMaxGap, nMinProbe, nUseVar, dCut, strOutputPath, strGrpName);


	/* release memory */
	Affy_BARData_Destroy(&vIP);
	Affy_BARData_Destroy(&vCT);


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_Peak()                                                       */
/*  Find peaks.                                                            */
/* ----------------------------------------------------------------------- */ 
int TileProbe_Peak(int nIPNum, int nCTNum, struct tagBARData *vIP, struct tagBARData *vCT,
				   int nBandWidth, int nMaxGap, int nMinProbe, int nUseVar, double dCut,
				   char strOutputPath[], char strGrpName[])
{
	/* define */
	int ni,nj,nk,ny,nz,nrank;
	int nCurrentPos;
	int nP1,nP2;
	int nQ1,nQ2;
	double dPrc = 0.1;
	int nCellNum;
	struct DOUBLEMATRIX *pWD;
	struct DOUBLEMATRIX *pWDSort;
	struct LONGMATRIX *pWDId;
	double dTemp,dTemp2;
	int nIPMATCol,nCTMATCol;
	struct DOUBLEMATRIX *pPosRegion;
	struct DOUBLEMATRIX *pNegRegion;
	struct DOUBLEMATRIX *pFDR;
	struct INTMATRIX *pCol;
	char strOutFile[MED_LINE_LENGTH];
	FILE *fpOut;

	/* init */
	nIPMATCol = 1+nIPNum;
	nCTMATCol = 1+nCTNum;

	/* scan the genome to compute statistics */
	for(ni=0; ni<vIP->nSeqNum; ni++)
	{
		nj = 0;
		
		nP1 = 0;
		nP2 = 0;

		for(nj=0; nj<vIP->vSeqData[ni]->nDataNum; nj++)
		{
			nCurrentPos = (int)(vIP->vSeqData[ni]->vData[0]->pMatElement[nj]);
			for(; nP1<nj; nP1++)
			{
				if( (int)(vIP->vSeqData[ni]->vData[0]->pMatElement[nP1]) >= (nCurrentPos-nBandWidth) )
					break;
			}

			for(; nP2<vIP->vSeqData[ni]->nDataNum; nP2++)
			{
				if( (int)(vIP->vSeqData[ni]->vData[0]->pMatElement[nP2]) > (nCurrentPos+nBandWidth) )
				{
					nP2--;
					break;
				}
			}
			if(nP2 >= vIP->vSeqData[ni]->nDataNum)
				nP2 = vIP->vSeqData[ni]->nDataNum-1;

			nCellNum = (nP2-nP1+1);
			if(nCellNum < nMinProbe)
			{
				continue;
			}
			
			nCellNum = (nP2-nP1+1)*nIPNum;
			nQ1 = (int)(nCellNum*dPrc);
			nQ2 = (int)(nCellNum*(1-dPrc))-1;

			
			pWD = NULL;
			pWDSort = NULL;
			pWDId = NULL;
			pWD = CreateDoubleMatrix(1, nCellNum);
			if(pWD == NULL)
			{
				printf("Error: TileProbe_Peak, cannot create matrix for sliding window!\n");
				exit(EXIT_FAILURE);
			}
			ny = 0;
			for(nz=0; nz<nIPNum; nz++)
			{
				for(nk=nP1; nk<=nP2; nk++)
				{
					pWD->pMatElement[ny] = vIP->vSeqData[ni]->vData[1+nz]->pMatElement[nk];
					ny++;
				}
			}
			DMSORTMERGEA_0(pWD, &pWDSort, &pWDId);

			dTemp = 0.0;
			for(nz=nQ1; nz<=nQ2; nz++)
			{
				dTemp += pWDSort->pMatElement[nz];
			}
			dTemp /= (nQ2-nQ1+1);
			vIP->vSeqData[ni]->vData[nIPMATCol]->pMatElement[nj] = dTemp;
			vIP->vSeqData[ni]->vData[nIPMATCol+2]->pMatElement[nj] = dTemp*sqrt((double)(nQ2-nQ1+1));

			DestroyDoubleMatrix(pWD);
			DestroyDoubleMatrix(pWDSort);
			DestroyLongMatrix(pWDId);

			if(vCT != NULL)
			{
				if( ((int)(vIP->vSeqData[ni]->vData[0]->pMatElement[nP1]) != (int)(vCT->vSeqData[ni]->vData[0]->pMatElement[nP1])) || 
					((int)(vIP->vSeqData[ni]->vData[0]->pMatElement[nP2]) != (int)(vCT->vSeqData[ni]->vData[0]->pMatElement[nP2])) )
				{
					printf("Error: TileProbe_Peak, coordinates in input bar files do not match!\n");
					exit(EXIT_FAILURE);
				}

				nCellNum = (nP2-nP1+1)*nCTNum;
				nQ1 = (int)(nCellNum*dPrc);
				nQ2 = (int)(nCellNum*(1-dPrc))-1;

				pWD = NULL;
				pWDSort = NULL;
				pWDId = NULL;
				pWD = CreateDoubleMatrix(1, nCellNum);
				if(pWD == NULL)
				{
					printf("Error: TileProbe_Peak, cannot create matrix for sliding window!\n");
					exit(EXIT_FAILURE);
				}
				ny = 0;
				for(nz=0; nz<nCTNum; nz++)
				{
					for(nk=nP1; nk<=nP2; nk++)
					{
						pWD->pMatElement[ny] = vCT->vSeqData[ni]->vData[1+nz]->pMatElement[nk];
						ny++;
					}
				}
				DMSORTMERGEA_0(pWD, &pWDSort, &pWDId);

				dTemp = 0.0;
				for(nz=nQ1; nz<=nQ2; nz++)
				{
					dTemp += pWDSort->pMatElement[nz];
				}
				dTemp /= (nQ2-nQ1+1);
				vCT->vSeqData[ni]->vData[nCTMATCol]->pMatElement[nj] = dTemp;
				dTemp = 0.0;
				for(nz=nQ1; nz<=nQ2; nz++)
				{
					dTemp2 = (pWDSort->pMatElement[nz]-vCT->vSeqData[ni]->vData[nCTMATCol]->pMatElement[nj]);
					dTemp += dTemp2*dTemp2;
				}
				dTemp /= (nQ2-nQ1);

				vCT->vSeqData[ni]->vData[nCTMATCol+1]->pMatElement[nj] = sqrt(dTemp);
				vIP->vSeqData[ni]->vData[nIPMATCol+2]->pMatElement[nj] = (vIP->vSeqData[ni]->vData[nIPMATCol]->pMatElement[nj]-vCT->vSeqData[ni]->vData[nCTMATCol]->pMatElement[nj])*sqrt((double)(nQ2-nQ1+1));
				if(nUseVar == 1)
				{
					vIP->vSeqData[ni]->vData[nIPMATCol+2]->pMatElement[nj] /= (vCT->vSeqData[ni]->vData[nCTMATCol+1]->pMatElement[nj]+1e-3);
				}

				DestroyDoubleMatrix(pWD);
				DestroyDoubleMatrix(pWDSort);
				DestroyLongMatrix(pWDId);
			}
		}
	}

	/* find positive peaks */
	pPosRegion = NULL;
	pPosRegion = TileProbe_Peak_Call(nIPMATCol+2, vIP, dCut, nMaxGap, 1, strOutputPath);
	if(pPosRegion == NULL)
	{
		printf("No peaks found!\n");
		return PROC_SUCCESS;
	}

	/* find negative peaks */
	pNegRegion = NULL;
	pNegRegion = TileProbe_Peak_Call(nIPMATCol+2, vIP, dCut, nMaxGap, 0, strOutputPath);

	/* compute FDR */
	pFDR = NULL;
	pFDR = CreateDoubleMatrix(1, pPosRegion->nHeight);
	if(pFDR == NULL)
	{
		printf("Error: TileProbe_Peak, cannot create matrix for computing FDR!\n");
		exit(EXIT_FAILURE);
	}

	TileMapv2_RegionFDR_Count(pPosRegion, pNegRegion, pFDR, 3); 
	TileMapv2_RegionFDR_Compute(pFDR);

	/* export results */
	pCol = NULL;
	pCol = CreateIntMatrix(1,vIP->nColNum);
	if(pCol == NULL)
	{
		printf("Error: TileProbe_Peak, cannot create output column information!\n");
		exit(EXIT_FAILURE);
	}
	
	printf("Exporting results...\n");
	sprintf(strOutFile, "%s_peak.bar", strOutputPath);
	pCol->pMatElement[0] = 1;
	pCol->pMatElement[nIPMATCol+2] = 1;
	Affy_SaveBAR_Columns_Fast(strOutFile, vIP, pCol);

	DestroyIntMatrix(pCol);

	/* export regions */
	sprintf(strOutFile, "%s_peak.cod", strOutputPath);
	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: TileProbe_Peak, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "#rank\tchr\tstart\tend\tstrand\tpeak_length\tmax_score\tmax_score_pos\tFDR\n");
	nrank = 0;
	for(ni=0; ni<pPosRegion->nHeight; ni++)
	{
		nj = pPosRegion->nHeight-1-ni;
		nk = (int)(DMGETAT(pPosRegion, nj, 0));
		ny = (int)(DMGETAT(pPosRegion, nj, 1))-12;
		nz = (int)(DMGETAT(pPosRegion, nj, 2))+12;
		
		if(strcmp(strGrpName, "")==0)
		{
			fprintf(fpOut, "%d\t%s\t%d\t%d\t+\t%d\t%f\t%d\t%f\n", (nrank+1), vIP->vSeqData[nk]->pSeqName->m_pString, 
				ny, nz, (nz-ny+1), DMGETAT(pPosRegion, nj, 3), (int)(DMGETAT(pPosRegion, nj, 4)), pFDR->pMatElement[nj]);
			nrank++;
		}
		else if(strcmp(strGrpName, vIP->vSeqData[nk]->pSeqGroupName->m_pString) == 0)
		{
			fprintf(fpOut, "%d\t%s\t%d\t%d\t+\t%d\t%f\t%d\t%f\n", (nrank+1), vIP->vSeqData[nk]->pSeqName->m_pString, 
				ny, nz, (nz-ny+1), DMGETAT(pPosRegion, nj, 3), (int)(DMGETAT(pPosRegion, nj, 4)), pFDR->pMatElement[nj]);
			nrank++;
		}
	}

	fclose(fpOut);
	
	/* release memory */
	DestroyDoubleMatrix(pFDR);
	DestroyDoubleMatrix(pPosRegion);
	DestroyDoubleMatrix(pNegRegion);


	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  TileProbe_Peak_Call()                                                  */
/*  Call peaks.                                                            */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *TileProbe_Peak_Call(int nTCol, struct tagBARData *vIP, 
	double dCut, int nMaxGap, int nPosVal, char strOutputPath[])
{
	/* define */
	FILE *fpOut;
	char strFileName[MED_LINE_LENGTH];
	int nStart,nEnd,nMaxPos;
	int ni,nj;
	int nPass;
	double dT,dMaxT;
	struct DOUBLEMATRIX *pReg;
	struct DOUBLEMATRIX *pRegSort;
	struct INTMATRIX *pType;
	struct INTMATRIX *pPriority;
	struct LONGMATRIX *pRegSortId;

	/* open temp file */
	sprintf(strFileName, "%s_peak.tmp", strOutputPath);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: TileProbe_Peak_Call, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	/* find peaks */
	for(ni=0; ni<vIP->nSeqNum; ni++)
	{
		nStart = -1;
		nEnd = -1;
		nPass = 0;
		dT = 0.0;
		dMaxT = 0.0;
		nMaxPos = 0;
		for(nj=0; nj<vIP->vSeqData[ni]->nDataNum; nj++)
		{
			if( (nPosVal == 1) & (vIP->vSeqData[ni]->vData[nTCol]->pMatElement[nj] >= dCut) )
			{
				nPass = 1;
				dT = vIP->vSeqData[ni]->vData[nTCol]->pMatElement[nj];
			}
			else if( (nPosVal == 0) & (vIP->vSeqData[ni]->vData[nTCol]->pMatElement[nj] <= -dCut) )
			{
				nPass = 1;
				dT = -vIP->vSeqData[ni]->vData[nTCol]->pMatElement[nj];
			}
			else
			{
				nPass = 0;
				dT = 0.0;
			}

			if(nPass == 1)
			{
				if(nStart < 0)
				{
					nStart = (int)(vIP->vSeqData[ni]->vData[0]->pMatElement[nj]);
					nEnd = nStart;
					dMaxT = dT;
					nMaxPos = nStart;
				}
				else
				{
					if( ((int)(vIP->vSeqData[ni]->vData[0]->pMatElement[nj])-nEnd) <= nMaxGap)
					{
						nEnd = (int)(vIP->vSeqData[ni]->vData[0]->pMatElement[nj]);
						if(dT > dMaxT)
						{
							dMaxT = dT;
							nMaxPos = nEnd;
						}
					}
					else
					{
						fprintf(fpOut, "%d\t%d\t%d\t%e\t%d\n", ni, nStart, nEnd, dMaxT, nMaxPos);

						nStart = (int)(vIP->vSeqData[ni]->vData[0]->pMatElement[nj]);
						nEnd = nStart;
						dMaxT = dT;
						nMaxPos = nStart;
					}
				}
			}
		}

		if(nStart > 0)
		{
			fprintf(fpOut, "%d\t%d\t%d\t%e\t%d\n", ni, nStart, nEnd, dMaxT, nMaxPos);
		}
	}


	/* close temp file */
	fclose(fpOut);

	/* load and sort the regions */
	pReg = NULL;
	pReg = DMLOAD(strFileName);
	if(pReg == NULL)
	{
		return NULL;
	}

	pType = NULL;
	pType = CreateIntMatrix(1, pReg->nWidth);
	if(pType == NULL)
	{
		printf("Error: TileProbe_Peak_Call, cannot create column data type!\n");
		exit(EXIT_FAILURE);
	}
		
	pType->pMatElement[0] = 2;
	pType->pMatElement[1] = 2;
	pType->pMatElement[2] = 2;
	pType->pMatElement[3] = 1;
	pType->pMatElement[4] = 2;

	pPriority = NULL;
	pPriority = CreateIntMatrix(1, 1);
	if(pPriority == NULL)
	{
		printf("Error: TileProbe_Peak_Call, cannot create column data type!\n");
		exit(EXIT_FAILURE);
	}
		
	pPriority->pMatElement[0] = 3;

	DMSORTROWS(pReg, pType, pPriority, &pRegSort, &pRegSortId);

	/* release memory */
	DestroyDoubleMatrix(pReg);
	DestroyIntMatrix(pType);
	DestroyIntMatrix(pPriority);
	DestroyLongMatrix(pRegSortId);
	RemoveFiles(strFileName);


	/* return */
	return pRegSort;
}


/* ----------------------------------------------------------------------- */ 
/*  TileProbe_MAT_Main()                                                   */
/*  MAT background correction for tileprobe.                               */
/* ----------------------------------------------------------------------- */ 
int TileProbe_MAT_Main(char strParamPath[])
{
	/* define */

	/* working path */
	char strProjectTitle[MED_LINE_LENGTH];
	char strCELPath[MED_LINE_LENGTH];
	char strBpmapPath[MED_LINE_LENGTH];
	char strWorkPath[MED_LINE_LENGTH];
	char strMaskPath[MED_LINE_LENGTH];
	
	/* arrays */
	int nLibNum = 0;
	struct tagString **vLibName = NULL;
	int nSampleNum = 0;
	int nArrayNum = 0;
	struct tagString **vCELPath;
	struct tagString **vArrayPath;
	struct tagString **vAlias;
	
	/* CEL files */
	int nTotalProbeNum = 0;
	int nRealProbeNum = 0;
	int nCELNumberCells = 0;
	int nCELTotalX = 0;
	int nCELTotalY = 0;
	struct tagCELData *pCELData;
	
	/* mask */
	int nRemoveMaskedCells = 0;
	int nRemoveOutlierCells = 0;
	struct INTMATRIX *pNumMaskCells;

	/* normalization */
	double dNormLowerBound = 1.0;
	int nNormLogTransform = 1;
	double dLog2 = log(2.0);
	
	/* intensity computation */
	struct tagBARData *pBARPos;
	int nExportMode = 0;

	/* LS matrix */
	struct DOUBLEMATRIX *pXX,*pXY,*pBeta;

	/* others */
	int ni,nj,nk;
	char strLine[LONG_LINE_LENGTH];
	
	FILE *fpIn;
	char *chSep;
	int nError = 0;

	/* init */

	/* load parameters */
	fpIn = NULL;
	fpIn = fopen(strParamPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileProbe_MAT_Main, cannot open the parameter file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strstr(strLine, "[Project title]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot load project title!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strProjectTitle, chSep);
		}
		else if(strstr(strLine, "[CEL directory]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot load CEL directory!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strCELPath, chSep);
			AdjustDirectoryPath(strCELPath);
		}
		else if(strstr(strLine, "[BPMAP directory]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot load BPMAP directory!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strBpmapPath, chSep);
			AdjustDirectoryPath(strBpmapPath);
		}
		else if(strstr(strLine, "[Working directory]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot load working directory!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strWorkPath, chSep);
			AdjustDirectoryPath(strWorkPath);
		}
		
		else if(strstr(strLine, "[No. of Libraries]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot load no. of libraries!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			nLibNum = atoi(chSep);
			if(nLibNum <= 0)
			{
				printf("Warning: No BPMAP libraries provided!");
				return PROC_SUCCESS;
			}

			vLibName = NULL;
			vLibName = (struct tagString **)calloc(nLibNum, sizeof(struct tagString *));
			if(vLibName == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot allocate memory for loading BPMAP lists!\n");
				exit(EXIT_FAILURE);
			}
		}
		else if(strstr(strLine, "[Libraries]") == strLine)
		{
			ni = 0;
			while(ni < nLibNum)
			{
				fgets(strLine, LONG_LINE_LENGTH, fpIn);
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] == '\0')
					continue;

				StringAddTail(vLibName+ni, strLine);
				ni++;
			}
		}
		else if(strstr(strLine, "[No. of Samples]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot load no. of samples!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			nSampleNum = atoi(chSep);
            if(nSampleNum <= 0)
			{
				printf("Error: TileProbe_MAT_Main, no arrays available!\n");
				exit(EXIT_FAILURE);
			}

			nArrayNum = (int)(nSampleNum*nLibNum);

			vCELPath = NULL;
			vCELPath = (struct tagString **)calloc(nArrayNum, sizeof(struct tagString *));
			if(vCELPath == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot allocate memory for tracking *.CEL files!\n");
				exit(EXIT_FAILURE);
			}

			vArrayPath = NULL;
			vArrayPath = (struct tagString **)calloc(nArrayNum, sizeof(struct tagString *));
			if(vArrayPath == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot allocate memory for tracking array files!\n");
				exit(EXIT_FAILURE);
			}

			vAlias = NULL;
			vAlias = (struct tagString **)calloc(nSampleNum, sizeof(struct tagString *));
			if(vAlias == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot allocate memory for tracking *.CEL files!\n");
				exit(EXIT_FAILURE);
			}

			pNumMaskCells = NULL;
			pNumMaskCells = CreateIntMatrix(1,nArrayNum);
			if(pNumMaskCells == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot allocate memory for tracking number of masked cells!\n");
				exit(EXIT_FAILURE);
			}
		}
		
		else if(strstr(strLine, "[Arrays]") == strLine)
		{
			ni = 0;
			while(ni < nSampleNum)
			{
				fgets(strLine, LONG_LINE_LENGTH, fpIn);
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] == '\0')
					continue;
				
				if(strLine[0] != '>')
				{
					printf("Error: TileProbe_MAT_Main, error when loading samples!\n");
					exit(EXIT_FAILURE);
				}

				chSep = strLine+1;
				StrTrimLeft(chSep);
				StringAddTail(vAlias+ni, chSep);
				
				nj = 0;
				while(nj < nLibNum)
				{
					fgets(strLine, LONG_LINE_LENGTH, fpIn);
					StrTrimLeft(strLine);
					StrTrimRight(strLine);
					if(strLine[0] == '\0')
						continue;

					StringAddTail(vCELPath+(ni*nLibNum)+nj, strLine);
					chSep = NULL;
					chSep = strrchr(strLine, '.');
					if(chSep != NULL)
						*chSep = '\0';
					StringAddTail(vArrayPath+(ni*nLibNum)+nj, strLine);
					nj++;
				}
				ni++;
			}
		}

		else if(strstr(strLine, "[Remove masked cells in the CEL files]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot load masked cell option!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			nRemoveMaskedCells = atoi(chSep);
		}

		else if(strstr(strLine, "[Remove outlier cells in the CEL files]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot load masked cell option!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			nRemoveOutlierCells = atoi(chSep);
		}
		else if(strLine[0] == '#')
		{
		}

		else
		{
			printf("Error: TileProbe_MAT_Main, unknown parameters!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	fclose(fpIn);

	printf("########################################\n");
	printf("# MAT background correction            #\n");
	printf("########################################\n");

	for(ni=0; ni<nLibNum; ni++)
	{
		printf("Processing %s...\n", vLibName[ni]->m_pString);
		sprintf(strLine, "%s%s", strBpmapPath, vLibName[ni]->m_pString);
		sprintf(strMaskPath, "%s%s.refmask", strWorkPath, vLibName[ni]->m_pString);
		
		/* load bpmap */
		pBARPos = NULL;
		pBARPos = TileProbe_BpmapToBAR(strLine, strMaskPath);
		if(pBARPos == NULL)
		{
			printf("Error: TileProbe_MAT_Main, empty bpmap file!\n");
			exit(EXIT_FAILURE);
		}

		/* prepare X'X */
		pXX = NULL;
		pXX = TileProbe_MAT_XX(pBARPos);
		if(pXX == NULL)
		{
			printf("Error: TileProbe_MAT_Main, cannot obtain X'X!\n");
			exit(EXIT_FAILURE);
		}

		/* mat */
		for(nj=0; nj<nSampleNum; nj++)
		{
			nk = nj*nLibNum+ni;
			printf("Loading %s \n", vCELPath[nk]->m_pString);
			sprintf(strLine, "%s%s", strCELPath, vCELPath[nk]->m_pString);

			/* load CEL */
			pCELData = TileMapv2_LoadCEL(strLine);
			if(pCELData == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot load *.CEL file!\n");
				exit(EXIT_FAILURE);
			}

			/* map intensity to probes, and obtain pXY */
			pXY = NULL;
			pXY = TileProbe_MAT_XY(pBARPos, pCELData);
			if(pXY == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot obtain X'Y!\n");
				exit(EXIT_FAILURE);
			}

			/* solve (X'X)Beta=X'Y */
			pBeta = NULL;
			/* TODO: solve AX=Y */
			pBeta = DMSOLVEEQU(pXX, pXY);

			if(pBeta == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot solve beta!\n");
				exit(EXIT_FAILURE);
			}

			/* clear memory */
			Affy_CELData_Destroy(&pCELData);

			/* MAT adjustment */
			sprintf(strLine, "%s%s.bar", strWorkPath, vArrayPath[nk]->m_pString);
			TileProbe_MAT_Correction(pBARPos, pBeta, strLine);

			/* clear memory */
			DestroyDoubleMatrix(pXY);
			DestroyDoubleMatrix(pBeta);
		}

		/* write cgw file */
		
		/* clear memory */
		Affy_BARData_Destroy(&pBARPos);
		DestroyDoubleMatrix(pXX);
	}


	/* destroy */
	for(ni=0; ni<nArrayNum; ni++)
	{
		DeleteString(vCELPath[ni]);
		DeleteString(vArrayPath[ni]);
	}
	free(vCELPath);
	free(vArrayPath);

	for(ni=0; ni<nSampleNum; ni++)
	{
		DeleteString(vAlias[ni]);
	}
	free(vAlias);

	for(ni=0; ni<nLibNum; ni++)
	{
		DeleteString(vLibName[ni]);
	}
	free(vLibName);

	DestroyIntMatrix(pNumMaskCells);

	printf("Done!\n");

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_BpmapToBAR()                                                 */
/*  Create BAR template from bpmap file. The function will create local    */
/*  repeat mask at the same time.                                          */
/* ----------------------------------------------------------------------- */
struct tagBARData *TileProbe_BpmapToBAR(char strBpmapFile[], char strMaskPath[])
{
	/* define */
	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));

	struct tagBARData *pBARPos = NULL;

	FILE *fpIn;
	FILE *fpOut;

	/* variables */
	char strFileType[9];
	float fVersion;
	unsigned int nSeqNum;

	/* seq info */
	unsigned int nSeqNameLen;

	struct INTMATRIX *vProbeMappingType;
	unsigned int nProbeMappingType;
	
	struct INTMATRIX *vSequenceFileOffset;
	unsigned int nSequenceFileOffset;
	
	struct INTMATRIX *vProbePairNum;
	unsigned int nProbePairNum;
	
	unsigned int nParamNum;
		
	/* map info */
	struct INTMATRIX *vSeqID;
	unsigned int nSeqID;
	int nProbeNum;
	int nInfoCol;
	int nTotalProbeNum = 0;
	int nMaskedProbeNum = 0;
	struct tagAffyBpMapUnit *pNewUnit,*pCUnit,*pPUnit;
	struct tagAffyBpMapUnit *pUnitList;
	struct DOUBLEMATRIX *pResizeMat;
	int nIgnore = 0;

	/* count */
	unsigned char cA,cC,cG,cT;
	unsigned char bChar;
	int ni,nj,nk,nx,ny,nz,nu;
	int nEndPos;
	struct BYTEMATRIX *pOB,*pResizeB;
	struct INTMATRIX *pOI,*pResizeI;

	/* load */
	fpIn = NULL;
	fpIn = fopen(strBpmapFile, "rb");
	if(fpIn == NULL)
	{
		printf("Error: TileProbe_BpmapToBAR, cannot open .bpmap file!\n");
		exit(EXIT_FAILURE);
	}
	fpOut = NULL;
	fpOut = fopen(strMaskPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: TileProbe_BpmapToBAR, cannot open .refmask file!\n");
		exit(EXIT_FAILURE);
	}
    fprintf(fpOut, "chromosome\tposition\tprobe_num\trepeat_num\tprobe_seq\tPMx\tPMy\tMMx\tMMy\n");
	
	/* load head */
	fread( strFileType, sizeof(char), 8, fpIn );
	strFileType[8] = '\0';
	AFFYBAR_READ_VERSION(fpIn, &fVersion);
	if(big_endian_fread(&nSeqNum, 4, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: TileProbe_BpmapToBAR, cannot load sequence number.\n");
		exit(EXIT_FAILURE);
	}
	printf("  BPMAP Version = %f\n", fVersion);
	printf("  Number of Sequences = %d\n", nSeqNum);
	if(nSeqNum <= 0)
	{
		printf("Warning: empty sequences!\n");
		return NULL;
	}

	/* create BAR object */
	pBARPos = Affy_BARData_Create();
	if(pBARPos == NULL)
	{
		printf("Error: TileProbe_BpmapToBAR, BARData object was not created!\n");
		exit(EXIT_FAILURE);
	}
	strcpy(pBARPos->strMagicnumber, "barr\r\n\032\n");
	pBARPos->fVersionnumber = 2.0;
    pBARPos->nSeqNum = nSeqNum;
	/* columns: 1 for coord, 1 for PMX, 1 for PMY, 1 for logPM value, 75/8 for A,C,G, 1 for nT, 4 for (A,C,G,T)^2 */
	pBARPos->nColNum = 4+10+1+4;
	nInfoCol = 4;
	pBARPos->pFieldType = CreateIntMatrix(1, pBARPos->nColNum);
	if(pBARPos->pFieldType == NULL)
	{
		printf("Error: TileProbe_BpmapToBAR, cannot allocate memory for saving field type.\n");
		exit(EXIT_FAILURE);
	}
	pBARPos->pFieldType->pMatElement[0] = 2;
	pBARPos->pFieldType->pMatElement[1] = 2;
	pBARPos->pFieldType->pMatElement[2] = 2;
	pBARPos->pFieldType->pMatElement[3] = 1;
	for(ni=4; ni<15; ni++)
		pBARPos->pFieldType->pMatElement[ni] = 7;
	for(; ni<19; ni++)
		pBARPos->pFieldType->pMatElement[ni] = 2;

	pBARPos->vSeqData = (struct tagBARSeq **)calloc(pBARPos->nSeqNum, sizeof(struct tagBARSeq *));
	if(pBARPos->vSeqData == NULL)
	{
		printf("Error: TileProbe_BpmapToBAR, cannot allocate memory for loading sequence data.\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<pBARPos->nSeqNum; ni++)
	{
		/* create BARSeq object */
		pBARPos->vSeqData[ni] = Affy_BARSeq_Create();
		if(pBARPos->vSeqData[ni] == NULL)
		{
			printf("Error: TileProbe_BpmapToBAR, cannot create BARSeq object.\n");
			exit(EXIT_FAILURE);
		}
	}

	vProbeMappingType = NULL;
	vProbeMappingType = CreateIntMatrix(nSeqNum,1);
	if(vProbeMappingType == NULL)
	{
		printf("Error: TileProbe_BpmapToBAR, cannot load probe mapping type!\n");
		exit(EXIT_FAILURE);
	}

	if(fVersion > 2.5)
	{
		vSequenceFileOffset = NULL;
		vSequenceFileOffset = CreateIntMatrix(nSeqNum,1);
		if(vSequenceFileOffset == NULL)
		{
			printf("Error: TileProbe_BpmapToBAR, cannot load sequence file offset!\n");
			exit(EXIT_FAILURE);
		}
	}

	vProbePairNum = NULL;
	vProbePairNum = CreateIntMatrix(nSeqNum,1);
	if(vProbePairNum == NULL)
	{
		printf("Error: TileProbe_BpmapToBAR, cannot load probe pair number!\n");
		exit(EXIT_FAILURE);
	}
	
	for(ni=0; ni<(int)nSeqNum; ni++)
	{
		/* for version 1.0 or later */
		/* load sequence name */
		if(big_endian_fread(&nSeqNameLen, 4, 1, fpIn, little_endian_machine) != 1)
		{
			printf("Error: TileProbe_BpmapToBAR, cannot load sequence name length.\n");
			exit(EXIT_FAILURE);
		}
		if(nSeqNameLen > 0)
		{
			pBARPos->vSeqData[ni]->pSeqName = CreateString(nSeqNameLen);
			if(pBARPos->vSeqData[ni]->pSeqName == NULL)
			{
				printf("Error: TileProbe_BpmapToBAR, cannot load sequence name.\n");
				exit(EXIT_FAILURE);
			}
			if(big_endian_fread(pBARPos->vSeqData[ni]->pSeqName->m_pString, 1, nSeqNameLen, fpIn, little_endian_machine) != nSeqNameLen)
			{
				printf("Error: TileProbe_BpmapToBAR, cannot load sequence name.\n");
				exit(EXIT_FAILURE);
			}
			pBARPos->vSeqData[ni]->pSeqName->m_pString[nSeqNameLen] = '\0';
		}

		/* for version 3.0 or later */
		/* load probe mapping type and sequence file offset */
		if(fVersion > 2.5)
		{
			if(big_endian_fread(&nProbeMappingType, 4, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: TileProbe_BpmapToBAR, cannot load probe mapping type.\n");
				exit(EXIT_FAILURE);
			}
			vProbeMappingType->pMatElement[ni] = (int)nProbeMappingType;
			if(big_endian_fread(&nSequenceFileOffset, 4, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: TileProbe_BpmapToBAR, cannot load sequence file offset.\n");
				exit(EXIT_FAILURE);
			}
			vSequenceFileOffset->pMatElement[ni] = (int)nSequenceFileOffset;
		}

		/* for version 1.0 or later */
		if(big_endian_fread(&nProbePairNum, 4, 1, fpIn, little_endian_machine) != 1)
		{
			printf("Error: TileProbe_BpmapToBAR, cannot load probe pair number.\n");
			exit(EXIT_FAILURE);
		}
		vProbePairNum->pMatElement[ni] = (int)nProbePairNum;
		nTotalProbeNum += (int)nProbePairNum;

		/* for version 2.0 or later */
		if(fVersion > 1.5)
		{
			/* read group name */
			if(big_endian_fread(&nSeqNameLen, 4, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: TileProbe_BpmapToBAR, cannot load group name length.\n");
				exit(EXIT_FAILURE);
			}
			if(nSeqNameLen > 0)
			{
				pBARPos->vSeqData[ni]->pSeqGroupName = CreateString(nSeqNameLen);
				if(pBARPos->vSeqData[ni]->pSeqGroupName == NULL)
				{
					printf("Error: TileProbe_BpmapToBAR, cannot load group name.\n");
					exit(EXIT_FAILURE);
				}
				if(big_endian_fread(pBARPos->vSeqData[ni]->pSeqGroupName->m_pString, 1, nSeqNameLen, fpIn, little_endian_machine) != nSeqNameLen)
				{
					printf("Error: TileProbe_BpmapToBAR, cannot load group name.\n");
					exit(EXIT_FAILURE);
				}
				pBARPos->vSeqData[ni]->pSeqGroupName->m_pString[nSeqNameLen] = '\0';
			}

			/* read version name */
			if(big_endian_fread(&nSeqNameLen, 4, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: TileProbe_BpmapToBAR, cannot load version number length.\n");
				exit(EXIT_FAILURE);
			}
			if(nSeqNameLen > 0)
			{
				pBARPos->vSeqData[ni]->pSeqVersion = CreateString(nSeqNameLen);
				if(pBARPos->vSeqData[ni]->pSeqVersion == NULL)
				{
					printf("Error: TileProbe_BpmapToBAR, cannot load version number.\n");
					exit(EXIT_FAILURE);
				}
				if(big_endian_fread(pBARPos->vSeqData[ni]->pSeqVersion->m_pString, 1, nSeqNameLen, fpIn, little_endian_machine) != nSeqNameLen)
				{
					printf("Error: TileProbe_BpmapToBAR, cannot load version number.\n");
					exit(EXIT_FAILURE);
				}
				pBARPos->vSeqData[ni]->pSeqVersion->m_pString[nSeqNameLen] = '\0';
			}

			/* read paramters */
			if(big_endian_fread(&nParamNum, 4, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: TileProbe_BpmapToBAR, cannot load number of parameters.\n");
				exit(EXIT_FAILURE);
			}
			pBARPos->vSeqData[ni]->nParamNum = (int)nParamNum;
			if(nParamNum > 0)
			{
				pBARPos->vSeqData[ni]->vParamName = (struct tagString **)calloc(pBARPos->vSeqData[ni]->nParamNum, sizeof(struct tagString *));
				pBARPos->vSeqData[ni]->vParamValue = (struct tagString **)calloc(pBARPos->vSeqData[ni]->nParamNum, sizeof(struct tagString *));

				if( (pBARPos->vSeqData[ni]->vParamName == NULL) || (pBARPos->vSeqData[ni]->vParamValue == NULL) )
				{
					printf("Error: TileProbe_BpmapToBAR, cannot allocate memory for loading parameter name/value pairs\n");
					exit(EXIT_FAILURE);
				}
			}
			
			for(nj=0; nj<(int)nParamNum; nj++)
			{
				if(big_endian_fread(&nSeqNameLen, 4, 1, fpIn, little_endian_machine) != 1)
				{
					printf("Error: TileProbe_BpmapToBAR, cannot load parameter name length.\n");
					exit(EXIT_FAILURE);
				}
				if(nSeqNameLen > 0)
				{
					pBARPos->vSeqData[ni]->vParamName[nj] = CreateString(nSeqNameLen);
					if(pBARPos->vSeqData[ni]->vParamName[nj] == NULL)
					{
						printf("Error: TileProbe_BpmapToBAR, cannot load parameter name.\n");
						exit(EXIT_FAILURE);
					}
					if(big_endian_fread(pBARPos->vSeqData[ni]->vParamName[nj]->m_pString, 1, nSeqNameLen, fpIn, little_endian_machine) != nSeqNameLen)
					{
						printf("Error: TileProbe_BpmapToBAR, cannot load parameter name.\n");
						exit(EXIT_FAILURE);
					}
					pBARPos->vSeqData[ni]->vParamName[nj]->m_pString[nSeqNameLen] = '\0';
				}

				if(big_endian_fread(&nSeqNameLen, 4, 1, fpIn, little_endian_machine) != 1)
				{
					printf("Error: TileProbe_BpmapToBAR, cannot load parameter value length.\n");
					exit(EXIT_FAILURE);
				}
				if(nSeqNameLen > 0)
				{
					pBARPos->vSeqData[ni]->vParamValue[nj] = CreateString(nSeqNameLen);
					if(pBARPos->vSeqData[ni]->vParamValue[nj] == NULL)
					{
						printf("Error: TileProbe_BpmapToBAR, cannot load parameter value.\n");
						exit(EXIT_FAILURE);
					}
					if(big_endian_fread(pBARPos->vSeqData[ni]->vParamValue[nj]->m_pString, 1, nSeqNameLen, fpIn, little_endian_machine) != nSeqNameLen)
					{
						printf("Error: TileProbe_BpmapToBAR, cannot load parameter value.\n");
						exit(EXIT_FAILURE);
					}
					pBARPos->vSeqData[ni]->vParamValue[nj]->m_pString[nSeqNameLen] = '\0';
				}
			}
		}
	}

	/* load map */
	vSeqID = NULL;
	vSeqID = CreateIntMatrix(nSeqNum,1);
	if(vSeqID == NULL)
	{
		printf("Error: TileProbe_BpmapToBAR, cannot load sequence id!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<(int)nSeqNum; ni++)
	{
		/* load seq id */
		if(big_endian_fread(&nSeqID, 4, 1, fpIn, little_endian_machine) != 1)
		{
			printf("Error: TileProbe_BpmapToBAR, cannot load sequence ID.\n");
			exit(EXIT_FAILURE);
		}
		vSeqID->pMatElement[ni] = nSeqID;

		/* create initial memory, 1 for logPM value, 75/8 for A,C,G, 1 for nT, 4 for (A,C,G,T)^2 */
		pBARPos->vSeqData[ni]->nDataNum = 0;
		pBARPos->vSeqData[ni]->nColNum = pBARPos->nColNum;
		/* 
		pBARPos->vSeqData[ni]->nColNum = 4;
		if(vProbeMappingType->pMatElement[ni] == 1)
		{
			pBARPos->vSeqData[ni]->nColNum = 4;
		}
		else
		{
			pBARPos->vSeqData[ni]->nColNum = 6;
		}
		nInfoCol = pBARPos->vSeqData[ni]->nColNum;		
		pBARPos->vSeqData[ni]->nColNum = pBARPos->vSeqData[ni]->nColNum+10+1+4;
		*/

		pBARPos->vSeqData[ni]->vData = (struct DOUBLEMATRIX **)calloc(pBARPos->vSeqData[ni]->nColNum, sizeof(struct DOUBLEMATRIX *));
		if(pBARPos->vSeqData[ni]->vData == NULL)
		{
			printf("Error: TileProbe_BpmapToBAR, cannot allocate memory for loading bpmap coordinate data!\n");
			exit(EXIT_FAILURE);
		}
		nProbeNum = vProbePairNum->pMatElement[ni];
		if(nProbeNum > 0)
		{
			for(nj=0; nj<nInfoCol; nj++)
			{
				pBARPos->vSeqData[ni]->vData[nj] = CreateDoubleMatrix(1, nProbeNum);
				if(pBARPos->vSeqData[ni]->vData[nj] == NULL)
				{
					printf("Error: TileProbe_BpmapToBAR, cannot allocate memory for loading bpmap coordinate data!\n");
					exit(EXIT_FAILURE);
				}
			}

			for(; nj<nInfoCol+11; nj++)
			{
				pBARPos->vSeqData[ni]->vData[nj] = (struct DOUBLEMATRIX *)(CreateByteMatrix(1, nProbeNum));
				if(pBARPos->vSeqData[ni]->vData[nj] == NULL)
				{
					printf("Error: TileProbe_BpmapToBAR, cannot allocate memory for loading bpmap sequence data!\n");
					exit(EXIT_FAILURE);
				}
			}

			for(; nj<pBARPos->vSeqData[ni]->nColNum; nj++)
			{
				pBARPos->vSeqData[ni]->vData[nj] = (struct DOUBLEMATRIX *)(CreateIntMatrix(1, nProbeNum));
				if(pBARPos->vSeqData[ni]->vData[nj] == NULL)
				{
					printf("Error: TileProbe_BpmapToBAR, cannot allocate memory for loading bpmap sequence-square data!\n");
					exit(EXIT_FAILURE);
				}
			}
		}


		pUnitList = NULL;
		
		/* load probes */
		for(nj=0; nj<nProbeNum; nj++)
		{
			/* load new unit */
			pNewUnit = NULL;
			pNewUnit = AffyBpMapUnitCreate();
			AffyBpMapUnitLoad_v3m2(pNewUnit, fpIn, vProbeMappingType->pMatElement[ni], little_endian_machine);
			
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

							fprintf(fpOut, "%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%f\t%d\n", 
								pBARPos->vSeqData[ni]->pSeqName->m_pString, pCUnit->nPos, pCUnit->nDepthNum, 
								pCUnit->nRepeatNum, pCUnit->strProbeSeq, 
								pCUnit->nPMX, pCUnit->nPMY, pCUnit->nMMX, pCUnit->nMMY,
								pCUnit->fMatchScore, pCUnit->bStrand);

							if(pCUnit->nRepeatNum < 2)
							{
								nk = pBARPos->vSeqData[ni]->nDataNum;
								pBARPos->vSeqData[ni]->vData[0]->pMatElement[nk] = pCUnit->nPos;
								pBARPos->vSeqData[ni]->vData[1]->pMatElement[nk] = pCUnit->nPMX;
								pBARPos->vSeqData[ni]->vData[2]->pMatElement[nk] = pCUnit->nPMY;
								/* 
								nInfoCol = 4;
								if(vProbeMappingType->pMatElement[ni] == 0)
								{
									pBARPos->vSeqData[ni]->vData[3]->pMatElement[nk] = pCUnit->nMMX;
									pBARPos->vSeqData[ni]->vData[4]->pMatElement[nk] = pCUnit->nMMY;
									nInfoCol = 6;
								}
								*/

								cA = 0;
								cC = 0;
								cG = 0;
								cT = 0;

								ny = 0;
								for(nx=0; nx<25; nx++)
								{
									if(pCUnit->strProbeSeq[nx] == 'A')
									{
										nz = ny/8;
										nu = 7-ny%8;
										bChar = 0x01;
										bChar <<=nu;
										pOB = (struct BYTEMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+nz]);
										pOB->pMatElement[nk] |= bChar;
										cA++;
									}
									else if(pCUnit->strProbeSeq[nx] == 'C')
									{
										nz = (ny+1)/8;
										nu = 7-(ny+1)%8;
										bChar = 0x01;
										bChar <<=nu;
										pOB = (struct BYTEMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+nz]);
										pOB->pMatElement[nk] |= bChar;
										cC++;
									}
									else if(pCUnit->strProbeSeq[nx] == 'G')
									{
										nz = (ny+2)/8;
										nu = 7-(ny+2)%8;
										bChar = 0x01;
										bChar <<=nu;
										pOB = (struct BYTEMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+nz]);
										pOB->pMatElement[nk] |= bChar;
										cG++;
									}
									else if(pCUnit->strProbeSeq[nx] == 'T')
									{
										cT++;
									}
									else
									{
										printf("Error: TileProbe_BpmapToBAR, cannot recognize the base in the probe sequence. \n");
										exit(EXIT_FAILURE);
									}

									ny += 3;
								}

								pOB = (struct BYTEMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+10]);
								pOB->pMatElement[nk] = cT;
								pOI = (struct INTMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+11]);
								pOI->pMatElement[nk] = (int)(cA*cA);
								pOI = (struct INTMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+12]);
								pOI->pMatElement[nk] = (int)(cC*cC);
								pOI = (struct INTMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+13]);
								pOI->pMatElement[nk] = (int)(cG*cG);
								pOI = (struct INTMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+14]);
								pOI->pMatElement[nk] = (int)(cT*cT);

								pBARPos->vSeqData[ni]->nDataNum += 1;
							}

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

			fprintf(fpOut, "%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%f\t%d\n", 
						pBARPos->vSeqData[ni]->pSeqName->m_pString, pCUnit->nPos, pCUnit->nDepthNum, 
						pCUnit->nRepeatNum, pCUnit->strProbeSeq, 
						pCUnit->nPMX, pCUnit->nPMY, pCUnit->nMMX, pCUnit->nMMY,
						pCUnit->fMatchScore, pCUnit->bStrand);

			if(pCUnit->nRepeatNum < 2)
			{
				nk = pBARPos->vSeqData[ni]->nDataNum;
				pBARPos->vSeqData[ni]->vData[0]->pMatElement[nk] = pCUnit->nPos;
				pBARPos->vSeqData[ni]->vData[1]->pMatElement[nk] = pCUnit->nPMX;
				pBARPos->vSeqData[ni]->vData[2]->pMatElement[nk] = pCUnit->nPMY;
				
				/* 
				nInfoCol = 4;
				if(vProbeMappingType->pMatElement[ni] == 0)
				{
					pBARPos->vSeqData[ni]->vData[3]->pMatElement[nk] = pCUnit->nMMX;
					pBARPos->vSeqData[ni]->vData[4]->pMatElement[nk] = pCUnit->nMMY;
					nInfoCol = 6;
				}
				*/

				cA = 0;
				cC = 0;
				cG = 0;
				cT = 0;

				ny = 0;
				for(nx=0; nx<25; nx++)
				{
					if(pCUnit->strProbeSeq[nx] == 'A')
					{
						nz = ny/8;
						nu = 7-ny%8;
						bChar = 0x01;
						bChar <<=nu;
						pOB = (struct BYTEMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+nz]);
						pOB->pMatElement[nk] |= bChar;
						cA++;
					}
					else if(pCUnit->strProbeSeq[nx] == 'C')
					{
						nz = (ny+1)/8;
						nu = 7-(ny+1)%8;
						bChar = 0x01;
						bChar <<=nu;
						pOB = (struct BYTEMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+nz]);
						pOB->pMatElement[nk] |= bChar;
						cC++;
					}
					else if(pCUnit->strProbeSeq[nx] == 'G')
					{
						nz = (ny+2)/8;
						nu = 7-(ny+2)%8;
						bChar = 0x01;
						bChar <<=nu;
						pOB = (struct BYTEMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+nz]);
						pOB->pMatElement[nk] |= bChar;
						cG++;
					}
					else if(pCUnit->strProbeSeq[nx] == 'T')
					{
						cT++;
					}
					else
					{
						printf("Error: TileProbe_BpmapToBAR, cannot recognize the base in the probe sequence. \n");
						exit(EXIT_FAILURE);
					}

					ny += 3;
				}

				pOB = (struct BYTEMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+10]);
				pOB->pMatElement[nk] = cT;
				pOI = (struct INTMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+11]);
				pOI->pMatElement[nk] = (int)(cA*cA);
				pOI = (struct INTMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+12]);
				pOI->pMatElement[nk] = (int)(cC*cC);
				pOI = (struct INTMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+13]);
				pOI->pMatElement[nk] = (int)(cG*cG);
				pOI = (struct INTMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+14]);
				pOI->pMatElement[nk] = (int)(cT*cT);

				pBARPos->vSeqData[ni]->nDataNum += 1;
			}

			AffyBpMapUnitDestroy(pCUnit);
		}

		/* 
		nInfoCol = 4;
		if(vProbeMappingType->pMatElement[ni] == 0)
		{
			pBARPos->vSeqData[ni]->vData[3]->pMatElement[nk] = pCUnit->nMMX;
			pBARPos->vSeqData[ni]->vData[4]->pMatElement[nk] = pCUnit->nMMY;
			nInfoCol = 6;
		}
		*/

		for(nj=0; nj<nInfoCol; nj++)
		{
			if(pBARPos->vSeqData[ni]->nDataNum > 0)
			{
				pResizeMat = NULL;
				pResizeMat = CreateDoubleMatrix(1, pBARPos->vSeqData[ni]->nDataNum);
				if(pResizeMat == NULL)
				{
					printf("Error: TileProbe_BpmapToBAR, cannot create an intermediate matrix to move data!\n");
					exit(EXIT_FAILURE);
				}
				memcpy(pResizeMat->pMatElement, pBARPos->vSeqData[ni]->vData[nj]->pMatElement, pBARPos->vSeqData[ni]->nDataNum*sizeof(double));
				DestroyDoubleMatrix(pBARPos->vSeqData[ni]->vData[nj]);
				pBARPos->vSeqData[ni]->vData[nj] = pResizeMat;
			}
			else
			{
				DestroyDoubleMatrix(pBARPos->vSeqData[ni]->vData[nj]);
				pBARPos->vSeqData[ni]->vData[nj] = NULL;
			}
		}

		for(; nj<nInfoCol+11; nj++)
		{
			if(pBARPos->vSeqData[ni]->nDataNum > 0)
			{
				pResizeB = NULL;
				pResizeB = CreateByteMatrix(1, pBARPos->vSeqData[ni]->nDataNum);
				if(pResizeB == NULL)
				{
					printf("Error: TileProbe_BpmapToBAR, cannot create an intermediate matrix to move data!\n");
					exit(EXIT_FAILURE);
				}
				pOB = (struct BYTEMATRIX *)(pBARPos->vSeqData[ni]->vData[nj]);
				memcpy(pResizeB->pMatElement, pOB->pMatElement, pBARPos->vSeqData[ni]->nDataNum*sizeof(unsigned char));
				DestroyByteMatrix(pOB);
				pBARPos->vSeqData[ni]->vData[nj] = (struct DOUBLEMATRIX *)pResizeB;
			}
			else
			{
				pOB = (struct BYTEMATRIX *)(pBARPos->vSeqData[ni]->vData[nj]);
				DestroyByteMatrix(pOB);
				pBARPos->vSeqData[ni]->vData[nj] = NULL;
			}
		}

		for(; nj<pBARPos->vSeqData[ni]->nColNum; nj++)
		{
			if(pBARPos->vSeqData[ni]->nDataNum > 0)
			{
				pResizeI = NULL;
				pResizeI = CreateIntMatrix(1, pBARPos->vSeqData[ni]->nDataNum);
				if(pResizeI == NULL)
				{
					printf("Error: TileProbe_BpmapToBAR, cannot create an intermediate matrix to move data!\n");
					exit(EXIT_FAILURE);
				}
				pOI = (struct INTMATRIX *)(pBARPos->vSeqData[ni]->vData[nj]);
				memcpy(pResizeI->pMatElement, pOI->pMatElement, pBARPos->vSeqData[ni]->nDataNum*sizeof(int));
				DestroyIntMatrix(pOI);
				pBARPos->vSeqData[ni]->vData[nj] = (struct DOUBLEMATRIX *)pResizeI;
			}
			else
			{
				pOI = (struct INTMATRIX *)(pBARPos->vSeqData[ni]->vData[nj]);
				DestroyIntMatrix(pOI);
				pBARPos->vSeqData[ni]->vData[nj] = NULL;
			}
		}

		/* update masked probe number */
		nMaskedProbeNum += pBARPos->vSeqData[ni]->nDataNum;
	}

	/* load tail if any */
	printf("  ");
	while(feof(fpIn) == 0)
	{
		fread(strFileType, sizeof(char), 1, fpIn);
		printf("%c", strFileType[0]);
	}
	printf("\n");

	/* clear memeory */
	DestroyIntMatrix(vProbeMappingType);
	if(fVersion > 2.5)
	{
		DestroyIntMatrix(vSequenceFileOffset);
	}
    DestroyIntMatrix(vProbePairNum);
	DestroyIntMatrix(vSeqID);

	/* close file */
	fclose(fpIn);
	fclose(fpOut);

	printf("  Probe number before masking = %d\n", nTotalProbeNum);
	printf("  Probe number after masking = %d\n", nMaskedProbeNum);

	/* return */
	return pBARPos;
}

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_MAT_XX()                                                     */
/*  Obtain X'X.                                                            */
/* ----------------------------------------------------------------------- */
struct DOUBLEMATRIX *TileProbe_MAT_XX(struct tagBARData *pBARPos)
{
	/* define */
	struct DOUBLEMATRIX *pXX = NULL;
	int nS = 75+1+4;
	int nBN = 75;
	int nBWN = (nBN/8)+1;
	int nCN = nS-nBN;
	int nInfoCol = 4;
	int ni,nj,nk,nx,ny,nz,nu;
	int nFinish;
	/* unsigned char *vB; */
	int vBI[25];
	int *vC;
	struct BYTEMATRIX *pBM;
	struct INTMATRIX *pIM;
	double dTemp;
	/* unsigned char cTemp; */

	/* init */
	if(pBARPos == NULL)
		return NULL;

	pXX = CreateDoubleMatrix(nS,nS);
	if(pXX == NULL)
	{
		printf("Error: TileProbe_XX, cannot allocate memory for X'X\n");
		exit(EXIT_FAILURE);
	}

	/* vB = NULL;
	vB = (unsigned char *)calloc(nBN, sizeof(unsigned char));
	if(vB == NULL)
	{
		printf("Error: TileProbe_XX, cannot allocate memory for vB\n");
		exit(EXIT_FAILURE);
	} */

	vC = NULL;
	vC = (int *)calloc(nCN, sizeof(int));
	if(vC == NULL)
	{
		printf("Error: TileProbe_XX, cannot allocate memory for vC\n");
		exit(EXIT_FAILURE);
	}

	/* compute */
	for(ni=0; ni<pBARPos->nSeqNum; ni++)
	{
		if(pBARPos->vSeqData[ni]->nDataNum <= 0)
			continue;

		for(nj=0; nj<pBARPos->vSeqData[ni]->nDataNum; nj++)
		{
			/* get vB */
			nx = 0;
			nu = 0;
			nFinish = 0;
			for(nk=0; nk<nBWN; nk++)
			{
				pBM = (struct BYTEMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+nk]);
				for(ny=0; ny<8; ny++)
				{
					nz = 7-ny;
					
					/* vB[nx] = (pBM->pMatElement[nj] >> nz) & 0x01; */
					if( ((pBM->pMatElement[nj] >> nz) & 0x01) != 0 )
					{
						vBI[nu] = nx;
						nu++;
					}

					nx++;
					if(nx >= nBN)
					{
						nFinish = 1;
						break;
					}
				}

				if(nFinish == 1)
					break;
			}

			/* get vC */
			pBM = (struct BYTEMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+nBWN]);
			vC[0] = pBM->pMatElement[nj];
			for(nk=1; nk<nCN; nk++)
			{
				pIM = (struct INTMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+nBWN+nk]);
				vC[nk] = pIM->pMatElement[nj];
			}

			if((nu+vC[0]) != 25)
			{
				printf("Error: TileProbe_XX, probe sequence matching wrong!\n");
				exit(EXIT_FAILURE);
			}

			/* get X'X */
			/* for(nx=0; nx<nBN; nx++)
			{
				for(ny=nx; ny<nBN; ny++)
				{
					cTemp = vB[nx] & vB[ny];
					dTemp = cTemp+DMGETAT(pXX, nx, ny);
					DMSETAT(pXX, nx, ny, dTemp);
				}

				for(ny=0; ny<nCN; ny++)
				{
					dTemp = vB[nx]*vC[ny]+DMGETAT(pXX, nx, nBN+ny);
					DMSETAT(pXX, nx, nBN+ny, dTemp);
				}

			} */

			for(nx=0; nx<nu; nx++)
			{
				for(ny=nx; ny<nu; ny++)
				{
					dTemp = 1+DMGETAT(pXX, vBI[nx], vBI[ny]);
					DMSETAT(pXX, vBI[nx], vBI[ny], dTemp);
				}

				for(ny=0; ny<nCN; ny++)
				{
					dTemp = vC[ny]+DMGETAT(pXX, vBI[nx], nBN+ny);
					DMSETAT(pXX, vBI[nx], nBN+ny, dTemp);
				}

			}

			for(nx=0; nx<nCN; nx++)
			{
				for(ny=nx; ny<nCN; ny++)
				{
					dTemp = vC[nx]*vC[ny]+DMGETAT(pXX, nBN+nx, nBN+ny);
					DMSETAT(pXX, nBN+nx, nBN+ny, dTemp);
				}
			}
		}
	}

	/* make X'X symmetrix */
	/* DMSAVE(pXX, "G:\\Projects\\TileProbe_Project\\mattest\\testXX.txt"); */
	for(nx=0; nx<nS; nx++)
	{
		for(ny=0; ny<nx; ny++)
		{
			dTemp = DMGETAT(pXX, ny, nx);
			DMSETAT(pXX, nx, ny, dTemp);
		}
	}
	/* DMSAVE(pXX, "G:\\Projects\\TileProbe_Project\\mattest\\testXXsym.txt"); */


	/* clear memory */
	/* free(vB); */
	free(vC);

	/* return */
	return pXX;
}

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_MAT_XY()                                                     */
/*  Obtain X'Y.                                                            */
/* ----------------------------------------------------------------------- */
struct DOUBLEMATRIX *TileProbe_MAT_XY(struct tagBARData *pBARPos, struct tagCELData *pCELData)
{
	/* define */
	double dNormLowerBound = 1.0;
	double dLog2 = log(2.0);
	struct DOUBLEMATRIX *pXY = NULL;
	int nS = 75+1+4;
	int nBN = 75;
	int nBWN = (nBN/8)+1;
	int nCN = nS-nBN;
	int nInfoCol = 4;
	int ni,nj,nk,nx,ny,nz,nu,nidx1;
	int nFinish;
	/* unsigned char *vB; */ 
	int vBI[25];
	int *vC;
	struct BYTEMATRIX *pBM;
	struct INTMATRIX *pIM;
	int nCELNumberCells;
	int nCELTotalX,nCELTotalY;
	double dPM;

	/* init */
	if(pBARPos == NULL)
		return NULL;

	nCELNumberCells = pCELData->nNumberCells;
	nCELTotalX = pCELData->nCols;
	nCELTotalY = pCELData->nRows;
	if(nCELTotalX*nCELTotalY != nCELNumberCells)
	{
		printf("Error: TileProbe_MAT_XY, CEL file Cols*Rows != NumberCells!\n");
		exit(EXIT_FAILURE);
	}

	pXY = CreateDoubleMatrix(nS,1);
	if(pXY == NULL)
	{
		printf("Error: TileProbe_XY, cannot allocate memory for X'Y\n");
		exit(EXIT_FAILURE);
	}

	/* vB = NULL;
	vB = (unsigned char *)calloc(nBN, sizeof(unsigned char));
	if(vB == NULL)
	{
		printf("Error: TileProbe_XY, cannot allocate memory for vB\n");
		exit(EXIT_FAILURE);
	} */

	vC = NULL;
	vC = (int *)calloc(nCN, sizeof(int));
	if(vC == NULL)
	{
		printf("Error: TileProbe_XY, cannot allocate memory for vC\n");
		exit(EXIT_FAILURE);
	}

	/* compute */
	for(ni=0; ni<pBARPos->nSeqNum; ni++)
	{
		if(pBARPos->vSeqData[ni]->nDataNum <= 0)
			continue;

		for(nj=0; nj<pBARPos->vSeqData[ni]->nDataNum; nj++)
		{
			/* get intensity */
			nidx1 = (int)(pBARPos->vSeqData[ni]->vData[2]->pMatElement[nj])*nCELTotalX+(int)(pBARPos->vSeqData[ni]->vData[1]->pMatElement[nj]);
			if(nidx1 >= nCELNumberCells)
			{
				printf("Error: TileProbe_MAT_XY, index out of range\n");
				exit(EXIT_FAILURE);
			}

			dPM = pCELData->pIntensity->pMatElement[nidx1];

			if(dPM < dNormLowerBound)
				dPM = dNormLowerBound;

			dPM = log(dPM)/dLog2;
			pBARPos->vSeqData[ni]->vData[3]->pMatElement[nj] = dPM;

			/* get vB */
			nx = 0;
			nu = 0;
			nFinish = 0;
			for(nk=0; nk<nBWN; nk++)
			{
				pBM = (struct BYTEMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+nk]);
				for(ny=0; ny<8; ny++)
				{
					nz = 7-ny;
					/* vB[nx] = (pBM->pMatElement[nj] >> nz) & 0x01; */
					if( ((pBM->pMatElement[nj] >> nz) & 0x01) != 0)
					{
						vBI[nu] = nx;
						nu++;
					}

					nx++;
					if(nx >= nBN)
					{
						nFinish = 1;
						break;
					}
				}

				if(nFinish == 1)
					break;
			}

			/* get vC */
			pBM = (struct BYTEMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+nBWN]);
			vC[0] = pBM->pMatElement[nj];
			for(nk=1; nk<nCN; nk++)
			{
				pIM = (struct INTMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+nBWN+nk]);
				vC[nk] = pIM->pMatElement[nj];
			}

			if( (nu+vC[0]) != 25)
			{
				printf("Error: TileProbe_MAT_XY, probe sequence matching problem!\n");
				exit(EXIT_FAILURE);
			}

			/* get X'Y */
			for(nx=0; nx<nu; nx++)
			{
				ny = vBI[nx];
				pXY->pMatElement[ny] += dPM; 
			}
			/* for(nx=0; nx<nBN; nx++)
			{
				pXY->pMatElement[nx] += vB[nx]*dPM; 
			} */
			for(nx=0; nx<nCN; nx++)
			{
				pXY->pMatElement[nBN+nx] += vC[nx]*dPM;
			}
		}
	}

	
	/* clear memory */
	/* free(vB); */
	free(vC);

	/* return */
	return pXY;
}

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_MAT_Correction()                                             */
/*  Background correction.                                                 */
/* ----------------------------------------------------------------------- */
int TileProbe_MAT_Correction(struct tagBARData *pBARPos, struct DOUBLEMATRIX *pBeta,
							 char strOutFile[])
{
	/* define */
	int nDataTotNum = 0;
	struct DOUBLEMATRIX *pE;
	struct DOUBLEMATRIX *pV;
	struct DOUBLEMATRIX *pESort;
	struct LONGMATRIX *pESortIndex;
	int nBinSize = 3000;
	int nS = 75+1+4;
	int nBN = 75;
	int nBWN = (nBN/8)+1;
	int nCN = nS-nBN;
	int nInfoCol = 4;
	int ni,nj,nk,nx,ny,nz,nh,nu;
	int nFinish;
	/* unsigned char *vB; */
	int vBI[25];
	int *vC;
	struct BYTEMATRIX *pBM;
	struct INTMATRIX *pIM;
	double dTemp,dSD;
	struct INTMATRIX *pCol;

	/* init */
	if(pBARPos == NULL || pBeta == NULL)
	{
		printf("Error: TileProbe_MAT_Correction, cannot correct background for %s!\n", strOutFile);
		return PROC_FAILURE;
	}

	/* 1: Create a single matrix for MAT correction */
	for(ni=0; ni<pBARPos->nSeqNum; ni++)
	{
		nDataTotNum += pBARPos->vSeqData[ni]->nDataNum;
	}
	pE = NULL;
	pE = CreateDoubleMatrix(1, nDataTotNum);
	if(pE == NULL)
	{
		printf("Error: TileProbe_MAT_Correction, cannot allocate memory for MAT adjustment!\n");
		exit(EXIT_FAILURE);
	}
	pV = NULL;
	pV = CreateDoubleMatrix(1, nDataTotNum);
	if(pV == NULL)
	{
		printf("Error: TileProbe_MAT_Correction, cannot allocate memory for MAT adjustment!\n");
		exit(EXIT_FAILURE);
	}

	/* 2: Sort fitted background */
	/* vB = NULL;
	vB = (unsigned char *)calloc(nBN, sizeof(unsigned char));
	if(vB == NULL)
	{
		printf("Error: TileProbe_XX, cannot allocate memory for vB\n");
		exit(EXIT_FAILURE);
	} */

	vC = NULL;
	vC = (int *)calloc(nCN, sizeof(int));
	if(vC == NULL)
	{
		printf("Error: TileProbe_MAT_Correction, cannot allocate memory for vC\n");
		exit(EXIT_FAILURE);
	}

	/* compute */
	nh = 0;
	for(ni=0; ni<pBARPos->nSeqNum; ni++)
	{
		if(pBARPos->vSeqData[ni]->nDataNum <= 0)
			continue;

		for(nj=0; nj<pBARPos->vSeqData[ni]->nDataNum; nj++)
		{
			/* get vB */
			nx = 0;
			nu = 0;
			nFinish = 0;
			for(nk=0; nk<nBWN; nk++)
			{
				pBM = (struct BYTEMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+nk]);
				for(ny=0; ny<8; ny++)
				{
					nz = 7-ny;
					/* vB[nx] = (pBM->pMatElement[nj] >> nz) & 0x01; */
					if( ((pBM->pMatElement[nj] >> nz) & 0x01) != 0 )
					{
						vBI[nu] = nx;
						nu++;
					}

					nx++;
					if(nx >= nBN)
					{
						nFinish = 1;
						break;
					}
				}

				if(nFinish == 1)
					break;
			}

			/* get vC */
			pBM = (struct BYTEMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+nBWN]);
			vC[0] = pBM->pMatElement[nj];
			for(nk=1; nk<nCN; nk++)
			{
				pIM = (struct INTMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+nBWN+nk]);
				vC[nk] = pIM->pMatElement[nj];
			}

			if( (nu+vC[0]) != 25 )
			{
				printf("Error: TileProbe_MAT_Correction, probe sequence matching problem!\n");
				exit(EXIT_FAILURE);
			}

			/* get expected value */
			dTemp = 0.0;
			for(nx=0; nx<nu; nx++)
			{
				ny = vBI[nx];
				dTemp += pBeta->pMatElement[ny];
			}
			for(nx=0; nx<nCN; nx++)
			{
				dTemp += vC[nx]*pBeta->pMatElement[nBN+nx];
			}
			pE->pMatElement[nh] = dTemp;
			pV->pMatElement[nh] = pBARPos->vSeqData[ni]->vData[3]->pMatElement[nj];
			pBARPos->vSeqData[ni]->vData[3]->pMatElement[nj] -= dTemp;
			
			/* nh ++ */
			nh++;
		}
	}

	/* clear memory */
	/* free(vB); */
	free(vC);

	/* 3: Stratified variance correction */
	DMSORTMERGEA_0(pE, &pESort, &pESortIndex);
	
	ni = 0;
	
	while(ni<nDataTotNum)
	{
		nj = ni+nBinSize-1;
		if(nj >= nDataTotNum)
			nj = nDataTotNum-1;
		else if( (nDataTotNum-nj)< (nBinSize/2) )
			nj = nDataTotNum-1;

		dTemp = 0.0;
		for(nk=ni; nk<=nj; nk++)
		{
			nx = pESortIndex->pMatElement[nk];
			dTemp += pV->pMatElement[nx];
		}
		dTemp /= (nj-ni+1);
		dSD = 0.0;
		for(nk=ni; nk<=nj; nk++)
		{
			nx = pESortIndex->pMatElement[nk];
			dSD += (pV->pMatElement[nx]-dTemp)*(pV->pMatElement[nx]-dTemp);
		}
		dSD /= (nj-ni);
		dSD = sqrt(dSD);

		for(nk=ni; nk<=nj; nk++)
		{
			nx = pESortIndex->pMatElement[nk];
			pV->pMatElement[nx] = dSD;
		}

		ni = nj+1;
	}

	DestroyDoubleMatrix(pESort);
	DestroyLongMatrix(pESortIndex);

	nh = 0;
	for(ni=0; ni<pBARPos->nSeqNum; ni++)
	{
		if(pBARPos->vSeqData[ni]->nDataNum <= 0)
			continue;

		for(nj=0; nj<pBARPos->vSeqData[ni]->nDataNum; nj++)
		{
			pBARPos->vSeqData[ni]->vData[3]->pMatElement[nj] /= pV->pMatElement[nh];
			nh++;
		}
	}

	/* 4: Export */
	pCol = NULL;
	pCol = CreateIntMatrix(1, pBARPos->nColNum);
	if(pCol == NULL)
	{
		printf("Error: TileProbe_MAT_Correction, cannot write background corrected intensities to file!\n");
		exit(EXIT_FAILURE);
	}
	pCol->pMatElement[0] = 1;
	pCol->pMatElement[3] = 1;
	Affy_SaveBAR_Columns_Fast(strOutFile, pBARPos, pCol);
	DestroyIntMatrix(pCol);

	/* release memory */
	DestroyDoubleMatrix(pE);
	DestroyDoubleMatrix(pV);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_BARBuild_Main()                                              */
/*  Build a probe background model for microarray BAR files.               */
/* ----------------------------------------------------------------------- */ 
int TileProbe_BARBuild_Main(char strInputPath[], char strOutputPath[], 
						 int nTakeLog, double dNormLowerBound, int nShrink, int nTest)
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	FILE *fpOut2;

	/* arrays */
	int nArrayNum = 0;
	int nSeqNum = 0;
	int nProbeNum = 0;
	int nProbeNumC;

	struct tagBARData **vBARData = NULL;
	struct DOUBLEMATRIX *pProbeMedian;
	struct DOUBLEMATRIX *pProbeInterR;
	struct DOUBLEMATRIX *pProbeData;
	struct DOUBLEMATRIX *pProbeDataSort;
	struct LONGMATRIX *pProbeDataSortId;
	
	/* others */
	int ni,nj,nk,nu;
	double dLog2 = log(2.0);
	int nQ50,nQ25,nQ05;
	char strOutFile[MED_LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	struct INTMATRIX *pCol;

	/* for IQR shrinking */
	double nt,d1,d2,v0,u0,beta,dB;
	double dpi = 4.0*atan(1.0);
	double dQm,dQv,dQm2,dQmin,dQmax;

	/* get array number */
	nArrayNum = 0;
	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileProbe_BARBuild_Main, cannot open the input array list!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		nArrayNum++;
	}
	
	fclose(fpIn);
	
	printf("########################################\n");
	printf("# Loading Raw Data                     #\n");
	printf("########################################\n");
	printf(" Array number = %d\n", nArrayNum);

	if(nArrayNum == 0)
	{
		return PROC_SUCCESS;
	}

	/* load arrays */
	vBARData = NULL;
	vBARData = (struct tagBARData **)calloc(nArrayNum, sizeof(struct tagBARData *));
	if(vBARData == NULL)
	{
		printf("Error: TileProbe_BARBuild_Main, cannot allocate memory for the array list!\n");
		exit(EXIT_FAILURE);
	}

	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileProbe_BARBuild_Main, cannot open the input array list!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		printf("Loading %s ...\n", strLine);
		vBARData[ni] = Affy_LoadBAR_Fast(strLine);
		if(vBARData[ni] == NULL)
		{
			printf("Error: TileProbe_BARBuild_Main, cannot load *.bar file!\n");
			exit(EXIT_FAILURE);
		}

		nProbeNumC = 0;
		for(nj=0; nj<vBARData[ni]->nSeqNum; nj++)
		{
			nProbeNumC += vBARData[ni]->vSeqData[nj]->nDataNum;
		}
		
		if(ni == 0)
		{
			nSeqNum = vBARData[ni]->nSeqNum;
			nProbeNum = nProbeNumC;
			
		}
		else
		{
			if( (nSeqNum != vBARData[ni]->nSeqNum) || (nProbeNum != nProbeNumC) )
			{
				printf("Error: TileProbe_BARBuild_Main, file dimensions do not match!\n");
				exit(EXIT_FAILURE);
			}
		}

		if(nTakeLog == 1)
		{
			for(nj=0; nj<vBARData[ni]->nSeqNum; nj++)
			{
				for(nk=0; nk<vBARData[ni]->vSeqData[nj]->nDataNum; nk++)
				{
					/* truncate */
					if(vBARData[ni]->vSeqData[nj]->vData[1]->pMatElement[nk] < dNormLowerBound)
						vBARData[ni]->vSeqData[nj]->vData[1]->pMatElement[nk] = dNormLowerBound;

					/* log transformation */
					vBARData[ni]->vSeqData[nj]->vData[1]->pMatElement[nk] = log(vBARData[ni]->vSeqData[nj]->vData[1]->pMatElement[nk])/dLog2;
				}
			}
		}
		else
		{
			for(nj=0; nj<vBARData[ni]->nSeqNum; nj++)
			{
				for(nk=0; nk<vBARData[ni]->vSeqData[nj]->nDataNum; nk++)
				{
					/* truncate */
					if(vBARData[ni]->vSeqData[nj]->vData[1]->pMatElement[nk] < dNormLowerBound)
						vBARData[ni]->vSeqData[nj]->vData[1]->pMatElement[nk] = dNormLowerBound;
				}
			}
		}

		ni++;
	}
	
	fclose(fpIn);

	if(ni != nArrayNum)
	{
		printf("Error: TileProbe_BARBuild_Main, array number inconsistent!\n");
		exit(EXIT_FAILURE);
	}

	/* build probe model */
	printf("########################################\n");
	printf("# Build Probe Model                    #\n");
	printf("########################################\n");

	if(nTest == 1)
	{
		sprintf(strOutFile, "%s.testq", strOutputPath);
		fpOut = NULL;
		fpOut = fopen(strOutFile, "w");
		if(fpOut == NULL)
		{
			printf("Error: TileProbe_BARBuild_Main, cannot open the output argument file!\n");
			exit(EXIT_FAILURE);
		}
		sprintf(strOutFile, "%s.testprobe", strOutputPath);
		fpOut2 = NULL;
		fpOut2 = fopen(strOutFile, "w");
		if(fpOut2 == NULL)
		{
			printf("Error: TileProbe_BARBuild_Main, cannot open the output argument file!\n");
			exit(EXIT_FAILURE);
		}
	}

	pProbeMedian = NULL;
	pProbeMedian = CreateDoubleMatrix(1,nProbeNum);
	if(pProbeMedian == NULL)
	{
		printf("Error: TileProbe_BARBuild_Main, cannot allocate memory for probe median!\n");
		exit(EXIT_FAILURE);
	}

	pProbeInterR = NULL;
	pProbeInterR = CreateDoubleMatrix(1,nProbeNum);
	if(pProbeInterR == NULL)
	{
		printf("Error: TileProbe_BARBuild_Main, cannot allocate memory for probe interquatile range!\n");
		exit(EXIT_FAILURE);
	}

	pProbeData = NULL;
	pProbeData = CreateDoubleMatrix(1,nArrayNum);
	if(pProbeData == NULL)
	{
		printf("Error: TileProbe_BARBuild_Main, cannot allocate memory for modeling probe effect!\n");
		exit(EXIT_FAILURE);
	}

	nQ50 = (int)(nArrayNum*0.5);
	if( (nArrayNum*0.5)-nQ50 > 1e-6 )
		nQ50++;
	nQ50--;

	nQ25 = (int)(nArrayNum*0.25);
	if( (nArrayNum*0.25)-nQ25 > 1e-6 )
		nQ25++;
	nQ25--;

	nQ05 = (int)(nArrayNum*0.05);
	if( (nArrayNum*0.05)-nQ05 > 1e-6 )
		nQ05++;
	nQ05--;

	if(nQ50 < 0)
		nQ50 = 0;
	if(nQ25 < 0)
		nQ25 = 0;
	if(nQ05 < 0)
		nQ05 = 0;
	if(nQ50 >= nArrayNum)
		nQ50 = nArrayNum-1;
	if(nQ25 >= nArrayNum)
		nQ25 = nArrayNum-1;
	if(nQ05 >= nArrayNum)
		nQ05 = nArrayNum-1;

	nk = 0;
	for(ni=0; ni<vBARData[0]->nSeqNum; ni++)
	{
		for(nj=0; nj<vBARData[0]->vSeqData[ni]->nDataNum; nj++)
		{
			if(nk%100000 == 0)
			{
				printf("%d ...\n", nk);
			}

			pProbeDataSort = NULL;
			pProbeDataSortId = NULL;

			for(nu=0; nu<nArrayNum; nu++)
			{
				pProbeData->pMatElement[nu] = vBARData[nu]->vSeqData[ni]->vData[1]->pMatElement[nj];
			}

			DMSORTMERGEA_0(pProbeData, &pProbeDataSort, &pProbeDataSortId);
			pProbeMedian->pMatElement[nk] = pProbeDataSort->pMatElement[nQ50];
			pProbeInterR->pMatElement[nk] = pProbeDataSort->pMatElement[nQ50]-pProbeDataSort->pMatElement[nQ25];
			if(pProbeInterR->pMatElement[nk] < 0.0)
			{
				/* pProbeInterR->pMatElement[nk] = pProbeInterR->pMatElement[nk]; */
				printf("Warning: IQR<0!\n");
			}

			if(nTest == 1)
			{
				fprintf(fpOut, "%e\t%e\n", pProbeInterR->pMatElement[nk], pProbeDataSort->pMatElement[nQ50]-pProbeDataSort->pMatElement[nQ05]); 
				if(rand_u() < 0.01)
				{
					fprintf(fpOut2, "%f", pProbeData->pMatElement[0]);
					for(nu=1; nu<nArrayNum; nu++)
					{
						fprintf(fpOut2, "\t%f", pProbeData->pMatElement[nu]);
					}
					fprintf(fpOut2, "\n");
				}
			}
			
			DestroyDoubleMatrix(pProbeDataSort);
			DestroyLongMatrix(pProbeDataSortId);

			nk++;
		}
	}

	if(nk != nProbeNum)
	{
		printf("Error: TileProbe_BARBuild_Main, probe number do not match!\n");
		exit(EXIT_FAILURE);
	}

	/* IQR shrinkage */
	if(nShrink != 0)
	{
		nt = nArrayNum;
		u0 = 0.25;
		v0 = 0.5;
		d1 = norminv(0.0,1.0,u0);
		d2 = norminv(0.0,1.0,v0);
		beta = 2*dpi*(u0*(1-u0)*exp(d1*d1)+v0*(1-v0)*exp(d2*d2)-2*u0*(1-v0)*exp((d1*d1+d2*d2)/2.0))/nt/pow((d2-d1),2.0);
		
		dQmin = 1e10;
		dQmax = -1e10;
		dQm = 0.0;
		dQm2 = 0.0;
		dQv = 0.0;
		for(nj=0; nj<nProbeNum; nj++)
		{
			if(pProbeInterR->pMatElement[nj] < dQmin)
				dQmin = pProbeInterR->pMatElement[nj];
			if(pProbeInterR->pMatElement[nj] > dQmax)
				dQmax = pProbeInterR->pMatElement[nj];

			dQm += pProbeInterR->pMatElement[nj];
			dQm2 += pProbeInterR->pMatElement[nj]*pProbeInterR->pMatElement[nj];
		}
		dQm /= nProbeNum;
		dQm2 /= nProbeNum;
		for(nj=0; nj<nProbeNum; nj++)
		{
			dQv += (pProbeInterR->pMatElement[nj]-dQm)*(pProbeInterR->pMatElement[nj]-dQm);
		}
		dQv /= (nProbeNum-1);

		
		dB = (beta/(1.0+beta))*(dQm2/dQv);

		if(dB>1.0)
			dB=1.0;
		else if(dB<0.0)
			dB=0.0;

		if(nShrink == -1)
			dB = 0.0;

		for(nj=0; nj<nProbeNum; nj++)
		{
			pProbeInterR->pMatElement[nj] = pProbeInterR->pMatElement[nj]*(1.0-dB)+dQm*dB;
		}
		
		printf("min IQR = %f; max IQR = %f\n", dQmin, dQmax);
		printf("B = %f\n", dB);
	}

	if(nTest == 1)
	{
		fclose(fpOut);
		fclose(fpOut2);
	}

	/* save file */
	printf("########################################\n");
	printf("# Export Probe Model                   #\n");
	printf("########################################\n");

	pCol = NULL;
	pCol = CreateIntMatrix(1, vBARData[0]->nColNum);
	if(pCol == NULL)
	{
		printf("Error: TileProbe_BARBuild_Main, cannot create column indicator for exporting fitted model!\n");
		exit(EXIT_FAILURE);
	}
	pCol->pMatElement[0] = 1;
	pCol->pMatElement[1] = 1;

	sprintf(strOutFile, "%s.prm", strOutputPath);
	nk = 0;
	for(ni=0; ni<vBARData[0]->nSeqNum; ni++)
	{
		for(nj=0; nj<vBARData[0]->vSeqData[ni]->nDataNum; nj++)
		{
			vBARData[0]->vSeqData[ni]->vData[1]->pMatElement[nj] = pProbeMedian->pMatElement[nk];
			nk++;
		}
	}
	if(nk != nProbeNum)
	{
		printf("Error: TileProbe_BARBuild_Main, probe number do not match!\n");
		exit(EXIT_FAILURE);
	}
	Affy_SaveBAR_Columns_Fast(strOutFile, vBARData[0], pCol);

	sprintf(strOutFile, "%s.prv", strOutputPath);
	nk = 0;
	for(ni=0; ni<vBARData[0]->nSeqNum; ni++)
	{
		for(nj=0; nj<vBARData[0]->vSeqData[ni]->nDataNum; nj++)
		{
			vBARData[0]->vSeqData[ni]->vData[1]->pMatElement[nj] = pProbeInterR->pMatElement[nk];
			nk++;
		}
	}
	if(nk != nProbeNum)
	{
		printf("Error: TileProbe_BARBuild_Main, probe number do not match!\n");
		exit(EXIT_FAILURE);
	}
	Affy_SaveBAR_Columns_Fast(strOutFile, vBARData[0], pCol);

	DestroyIntMatrix(pCol);

	/* release memory */
	for(ni=0; ni<nArrayNum; ni++)
	{
		Affy_BARData_Destroy(vBARData+ni);
	}
	free(vBARData);

	DestroyDoubleMatrix(pProbeMedian);
	DestroyDoubleMatrix(pProbeInterR);
	DestroyDoubleMatrix(pProbeData);

	printf("Done!\n");

	/* return */
	return PROC_SUCCESS;
}

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
						double dB)
{
	/* define */
	FILE *fpIn;
	struct tagBARData *pProbeModel;
	struct tagBARData *pProbeQuan;
	struct tagString **vBARFile;

	struct DOUBLEMATRIX *pArray = NULL;
	struct DOUBLEMATRIX *pArrayS = NULL;
	struct LONGMATRIX *pSortId = NULL;
	struct INTMATRIX *pCol = NULL;

	/* arrays */
	int nArrayNum = 0;
	int nProbeNum = 0;
	int nProbeNumC;
	struct tagBARData *pBARData = NULL;
	
	/* others */
	int ni,nj,nk;
	double dLog2 = log(2.0);
	char strOutFile[MED_LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	char strTemp[MED_LINE_LENGTH];
	char strOutFolder[MED_LINE_LENGTH];
	double dSDM;

	/* init */
	strcpy(strOutFolder, strOutputPath);
	AdjustDirectoryPath(strOutFolder);


	/* get array number */
	if(nInputType == 1)
	{
		nArrayNum = 0;
		fpIn = NULL;
		fpIn = fopen(strInputPath, "r");
		if(fpIn == NULL)
		{
			printf("Error: TileProbe_BARNorm_Main, cannot open the input array list!\n");
			exit(EXIT_FAILURE);
		}

		while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;

			nArrayNum++;
		}
		
		fclose(fpIn);
	}
	else
	{
		nArrayNum = 1;
	}

	vBARFile = NULL;
	vBARFile = (struct tagString **)calloc(nArrayNum, sizeof(struct tagString *));
	if(vBARFile == NULL)
	{
		printf("Error: TileProbe_BARNorm_Main, cannot allocate memory for the array list!\n");
		exit(EXIT_FAILURE);
	}

	if(nInputType == 1)
	{
		fpIn = NULL;
		fpIn = fopen(strInputPath, "r");
		if(fpIn == NULL)
		{
			printf("Error: TileProbe_BARNorm_Main, cannot open the input array list!\n");
			exit(EXIT_FAILURE);
		}

		ni = 0;
		while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;

			StringAddTail(vBARFile+ni, strLine);
			ni++;
		}
		
		fclose(fpIn);

		if(ni != nArrayNum)
		{
			printf("Error: TileProbe_BARNorm_Main, array number inconsistent!\n");
			exit(EXIT_FAILURE);
		}
	}
	else
	{
		StringAddTail(vBARFile+0, strInputPath);
	}

	
	/* load model */
	sprintf(strLine, "%s.prm", strModelPath);
	pProbeModel = NULL;
	pProbeModel = Affy_LoadBAR_Fast(strLine);
	if(pProbeModel == NULL)
	{
		printf("Error: TileProbe_BARNorm_Main, cannot load probe model!\n");
		exit(EXIT_FAILURE);
	}

	sprintf(strLine, "%s.prv", strModelPath);
	pProbeQuan = NULL;
	pProbeQuan = Affy_LoadBAR_Fast(strLine);
	if(pProbeQuan == NULL)
	{
		printf("Error: TileProbe_BARNorm_Main, cannot load probe model!\n");
		exit(EXIT_FAILURE);
	}

	nProbeNum = 0;
	for(ni=0; ni<pProbeQuan->nSeqNum; ni++)
	{
		nProbeNum += pProbeQuan->vSeqData[ni]->nDataNum;
	}

	/* shrink variance */
	if(dB > 0.0)
	{
		pArray = NULL;
		pArray = CreateDoubleMatrix(1, nProbeNum);
		if(pArray == NULL)
		{
			printf("Error: TileProbe_BARNorm_Main, cannot create memory for shrinking IQR!\n");
			exit(EXIT_FAILURE);
		}
		
		nk = 0;
		for(ni=0; ni<pProbeQuan->nSeqNum; ni++)
		{
			for(nj=0; nj<pProbeQuan->vSeqData[ni]->nDataNum; nj++)
			{
				pArray->pMatElement[nk] = pProbeQuan->vSeqData[ni]->vData[1]->pMatElement[nj];
				nk++;
			}
		}
		if(nk != nProbeNum)
		{
			printf("Error: TileProbe_BARNorm_Main, inconsistent probe number!\n");
			exit(EXIT_FAILURE);
		}

		pArrayS = NULL;
		pSortId = NULL;
		DMSORTMERGEA_0(pArray, &pArrayS, &pSortId);

		ni = (int)(nProbeNum*0.5);
		dSDM = pArrayS->pMatElement[ni];
		
		for(ni=0; ni<pProbeQuan->nSeqNum; ni++)
		{
			for(nj=0; nj<pProbeQuan->vSeqData[ni]->nDataNum; nj++)
			{
				pProbeQuan->vSeqData[ni]->vData[1]->pMatElement[nj] = (1.0-dB)*pProbeQuan->vSeqData[ni]->vData[1]->pMatElement[nj]+dB*dSDM;
			}
		}

		DestroyDoubleMatrix(pArray);
		DestroyDoubleMatrix(pArrayS);
		DestroyLongMatrix(pSortId);
	}

	printf("########################################\n");
	printf("# Processing                           #\n");
	printf("########################################\n");
	printf(" Array number = %d\n", nArrayNum);

	if(nArrayNum == 0)
	{
		return PROC_SUCCESS;
	}

	pCol = NULL;
	pCol = CreateIntMatrix(1,2);
	if(pCol == NULL)
	{
		printf("Error: TileProbe_BARNorm_Main, cannot create column index for exporting bar files!\n");
		exit(EXIT_FAILURE);
	}
	pCol->pMatElement[0] = 1;
	pCol->pMatElement[1] = 1;

	for(ni=0; ni<nArrayNum; ni++)
	{
		/* load data */
		printf("Processing %s ...\n", vBARFile[ni]->m_pString);
		pBARData = NULL;
		pBARData = Affy_LoadBAR_Fast(vBARFile[ni]->m_pString);
		if(pBARData == NULL)
		{
			printf("Error: TileProbe_BARNorm_Main, cannot load *.bar file!\n");
			exit(EXIT_FAILURE);
		}

		nProbeNumC = 0;
		for(nj=0; nj<pBARData->nSeqNum; nj++)
			nProbeNumC += pBARData->vSeqData[nj]->nDataNum;

		if( (pProbeModel->nSeqNum != pBARData->nSeqNum) || (nProbeNum != nProbeNumC) )
		{
			printf("Error: TileProbe_BARNorm_Main, BAR file dimension (%s) does not match with the model!\n", vBARFile[ni]->m_pString);
			exit(EXIT_FAILURE);
		}

		/* pre-processing */
		if(nTakeLog == 1)
		{
			for(nj=0; nj<pBARData->nSeqNum; nj++)
			{
				for(nk=0; nk<pBARData->vSeqData[nj]->nDataNum; nk++)
				{
					/* truncate */
					if(pBARData->vSeqData[nj]->vData[1]->pMatElement[nk] < dNormLowerBound)
						pBARData->vSeqData[nj]->vData[1]->pMatElement[nk] = dNormLowerBound;

					/* log transformation */
					pBARData->vSeqData[nj]->vData[1]->pMatElement[nk] = log(pBARData->vSeqData[nj]->vData[1]->pMatElement[nk])/dLog2;
				}
			}
		}
		else
		{
			for(nj=0; nj<pBARData->nSeqNum; nj++)
			{
				for(nk=0; nk<pBARData->vSeqData[nj]->nDataNum; nk++)
				{
					/* truncate */
					if(pBARData->vSeqData[nj]->vData[1]->pMatElement[nk] < dNormLowerBound)
						pBARData->vSeqData[nj]->vData[1]->pMatElement[nk] = dNormLowerBound;
				}
			}
		}

		/* reset bar values */
		for(nj=0; nj<pBARData->nSeqNum; nj++)
		{
			for(nk=0; nk<pBARData->vSeqData[nj]->nDataNum; nk++)
			{
				/* truncate */
				pBARData->vSeqData[nj]->vData[1]->pMatElement[nk] = (pBARData->vSeqData[nj]->vData[1]->pMatElement[nk]-pProbeModel->vSeqData[nj]->vData[1]->pMatElement[nk])/(pProbeQuan->vSeqData[nj]->vData[1]->pMatElement[nk]+1e-6);
			}
		}

		/* export results */
		GetFileName(vBARFile[ni]->m_pString, strTemp);
		sprintf(strOutFile, "%s%s.tpnorm.bar", strOutFolder, strTemp);
		Affy_SaveBAR_Columns_Fast(strOutFile, pBARData, pCol);

		Affy_BARData_Destroy(&pBARData);
	}
		
	/* release memory */
	DestroyIntMatrix(pCol);
	for(ni=0; ni<nArrayNum; ni++)
	{
		DeleteString(vBARFile[ni]);
		vBARFile[ni] = NULL;
	}
	free(vBARFile);

	Affy_BARData_Destroy(&pProbeModel);
	Affy_BARData_Destroy(&pProbeQuan);
	
	printf("Done!\n");

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  TileProbe_MATv2_Main()                                                 */
/*  MAT background correction for tileprobe.                               */
/* ----------------------------------------------------------------------- */ 
int TileProbe_MATv2_Main(char strParamPath[])
{
	/* define */

	/* working path */
	char strProjectTitle[MED_LINE_LENGTH];
	char strCELPath[MED_LINE_LENGTH];
	char strBpmapPath[MED_LINE_LENGTH];
	char strWorkPath[MED_LINE_LENGTH];
	char strMaskPath[MED_LINE_LENGTH];
	char strGenomeGrp[255];
	int nIncludeNonGrp = 0;
	
	/* arrays */
	int nLibNum = 0;
	struct tagString **vLibName = NULL;
	int nSampleNum = 0;
	int nArrayNum = 0;
	struct tagString **vCELPath;
	struct tagString **vArrayPath;
	struct tagString **vAlias;
	
	/* CEL files */
	int nTotalProbeNum = 0;
	int nRealProbeNum = 0;
	int nCELNumberCells = 0;
	int nCELTotalX = 0;
	int nCELTotalY = 0;
	struct tagCELData *pCELData;
	
	/* mask */
	int nRemoveMaskedCells = 0;
	int nRemoveOutlierCells = 0;
	struct INTMATRIX *pNumMaskCells;

	/* normalization */
	double dNormLowerBound = 1.0;
	int nNormLogTransform = 1;
	double dLog2 = log(2.0);
	
	/* intensity computation */
	struct tagBARData *pBARPos;
	int nExportMode = 0;

	/* LS matrix */
	struct DOUBLEMATRIX *pXX,*pXY,*pBeta;

	/* others */
	int ni,nj,nk;
	char strLine[LONG_LINE_LENGTH];
	
	FILE *fpIn;
	char *chSep;
	int nError = 0;

	/* init */
	strcpy(strGenomeGrp, "");

	/* load parameters */
	fpIn = NULL;
	fpIn = fopen(strParamPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileProbe_MAT_Main, cannot open the parameter file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strstr(strLine, "[Project title]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot load project title!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strProjectTitle, chSep);
		}
		else if(strstr(strLine, "[CEL directory]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot load CEL directory!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strCELPath, chSep);
			AdjustDirectoryPath(strCELPath);
		}
		else if(strstr(strLine, "[BPMAP directory]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot load BPMAP directory!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strBpmapPath, chSep);
			AdjustDirectoryPath(strBpmapPath);
		}
		else if(strstr(strLine, "[GenomeGrp]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot load Genome Group!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strGenomeGrp, chSep);
		}
		else if(strstr(strLine, "[Include Probes not in GenomeGrp]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot load Genome Group!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			nIncludeNonGrp = atoi(chSep);
		}
		else if(strstr(strLine, "[Working directory]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot load working directory!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strWorkPath, chSep);
			AdjustDirectoryPath(strWorkPath);
		}
		
		else if(strstr(strLine, "[No. of Libraries]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot load no. of libraries!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			nLibNum = atoi(chSep);
			if(nLibNum <= 0)
			{
				printf("Warning: No BPMAP libraries provided!");
				return PROC_SUCCESS;
			}

			vLibName = NULL;
			vLibName = (struct tagString **)calloc(nLibNum, sizeof(struct tagString *));
			if(vLibName == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot allocate memory for loading BPMAP lists!\n");
				exit(EXIT_FAILURE);
			}
		}
		else if(strstr(strLine, "[Libraries]") == strLine)
		{
			ni = 0;
			while(ni < nLibNum)
			{
				fgets(strLine, LONG_LINE_LENGTH, fpIn);
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] == '\0')
					continue;

				StringAddTail(vLibName+ni, strLine);
				ni++;
			}
		}
		else if(strstr(strLine, "[No. of Samples]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot load no. of samples!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			nSampleNum = atoi(chSep);
            if(nSampleNum <= 0)
			{
				printf("Error: TileProbe_MAT_Main, no arrays available!\n");
				exit(EXIT_FAILURE);
			}

			nArrayNum = (int)(nSampleNum*nLibNum);

			vCELPath = NULL;
			vCELPath = (struct tagString **)calloc(nArrayNum, sizeof(struct tagString *));
			if(vCELPath == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot allocate memory for tracking *.CEL files!\n");
				exit(EXIT_FAILURE);
			}

			vArrayPath = NULL;
			vArrayPath = (struct tagString **)calloc(nArrayNum, sizeof(struct tagString *));
			if(vArrayPath == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot allocate memory for tracking array files!\n");
				exit(EXIT_FAILURE);
			}

			vAlias = NULL;
			vAlias = (struct tagString **)calloc(nSampleNum, sizeof(struct tagString *));
			if(vAlias == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot allocate memory for tracking *.CEL files!\n");
				exit(EXIT_FAILURE);
			}

			pNumMaskCells = NULL;
			pNumMaskCells = CreateIntMatrix(1,nArrayNum);
			if(pNumMaskCells == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot allocate memory for tracking number of masked cells!\n");
				exit(EXIT_FAILURE);
			}
		}
		
		else if(strstr(strLine, "[Arrays]") == strLine)
		{
			ni = 0;
			while(ni < nSampleNum)
			{
				fgets(strLine, LONG_LINE_LENGTH, fpIn);
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] == '\0')
					continue;
				
				if(strLine[0] != '>')
				{
					printf("Error: TileProbe_MAT_Main, error when loading samples!\n");
					exit(EXIT_FAILURE);
				}

				chSep = strLine+1;
				StrTrimLeft(chSep);
				StringAddTail(vAlias+ni, chSep);
				
				nj = 0;
				while(nj < nLibNum)
				{
					fgets(strLine, LONG_LINE_LENGTH, fpIn);
					StrTrimLeft(strLine);
					StrTrimRight(strLine);
					if(strLine[0] == '\0')
						continue;

					StringAddTail(vCELPath+(ni*nLibNum)+nj, strLine);
					chSep = NULL;
					chSep = strrchr(strLine, '.');
					if(chSep != NULL)
						*chSep = '\0';
					StringAddTail(vArrayPath+(ni*nLibNum)+nj, strLine);
					nj++;
				}
				ni++;
			}
		}

		else if(strstr(strLine, "[Remove masked cells in the CEL files]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot load masked cell option!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			nRemoveMaskedCells = atoi(chSep);
		}

		else if(strstr(strLine, "[Remove outlier cells in the CEL files]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot load masked cell option!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			nRemoveOutlierCells = atoi(chSep);
		}
		else if(strLine[0] == '#')
		{
		}

		else
		{
			printf("Error: TileProbe_MAT_Main, unknown parameters!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	fclose(fpIn);

	printf("########################################\n");
	printf("# MAT background correction            #\n");
	printf("########################################\n");

	for(ni=0; ni<nLibNum; ni++)
	{
		printf("Processing %s...\n", vLibName[ni]->m_pString);
		sprintf(strLine, "%s%s", strBpmapPath, vLibName[ni]->m_pString);
		sprintf(strMaskPath, "%s%s.refmask", strWorkPath, vLibName[ni]->m_pString);
		
		/* load bpmap */
		pBARPos = NULL;
		pBARPos = TileProbe_BpmapToBARv2(strLine, strMaskPath, strGenomeGrp, nIncludeNonGrp);
		if(pBARPos == NULL)
		{
			printf("Error: TileProbe_MAT_Main, empty bpmap file!\n");
			exit(EXIT_FAILURE);
		}

		/* prepare X'X */
		pXX = NULL;
		pXX = TileProbe_MATv2_XX(pBARPos, strGenomeGrp, nIncludeNonGrp);
		if(pXX == NULL)
		{
			printf("Error: TileProbe_MAT_Main, cannot obtain X'X!\n");
			exit(EXIT_FAILURE);
		}

		/* mat */
		for(nj=0; nj<nSampleNum; nj++)
		{
			nk = nj*nLibNum+ni;
			printf("Loading %s \n", vCELPath[nk]->m_pString);
			sprintf(strLine, "%s%s", strCELPath, vCELPath[nk]->m_pString);

			/* load CEL */
			pCELData = TileMapv2_LoadCEL(strLine);
			if(pCELData == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot load *.CEL file!\n");
				exit(EXIT_FAILURE);
			}

			/* map intensity to probes, and obtain pXY */
			pXY = NULL;
			pXY = TileProbe_MATv2_XY(pBARPos, pCELData, strGenomeGrp, nIncludeNonGrp);
			if(pXY == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot obtain X'Y!\n");
				exit(EXIT_FAILURE);
			}

			/* solve (X'X)Beta=X'Y */
			pBeta = NULL;
			/* TODO: solve AX=Y */
			pBeta = DMSOLVEEQU(pXX, pXY);

			if(pBeta == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot solve beta!\n");
				exit(EXIT_FAILURE);
			}

			/* clear memory */
			Affy_CELData_Destroy(&pCELData);

			/* MAT adjustment */
			sprintf(strLine, "%s%s.bar", strWorkPath, vArrayPath[nk]->m_pString);
			TileProbe_MATv2_Correction(pBARPos, pBeta, strLine, strGenomeGrp, nIncludeNonGrp);

			/* clear memory */
			DestroyDoubleMatrix(pXY);
			DestroyDoubleMatrix(pBeta);
		}

		/* write cgw file */
		
		/* clear memory */
		Affy_BARData_Destroy(&pBARPos);
		DestroyDoubleMatrix(pXX);
	}


	/* destroy */
	for(ni=0; ni<nArrayNum; ni++)
	{
		DeleteString(vCELPath[ni]);
		DeleteString(vArrayPath[ni]);
	}
	free(vCELPath);
	free(vArrayPath);

	for(ni=0; ni<nSampleNum; ni++)
	{
		DeleteString(vAlias[ni]);
	}
	free(vAlias);

	for(ni=0; ni<nLibNum; ni++)
	{
		DeleteString(vLibName[ni]);
	}
	free(vLibName);

	DestroyIntMatrix(pNumMaskCells);

	printf("Done!\n");

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_BpmapToBARv2()                                               */
/*  Create BAR template from bpmap file. The function will create local    */
/*  repeat mask at the same time.                                          */
/* ----------------------------------------------------------------------- */
struct tagBARData *TileProbe_BpmapToBARv2(char strBpmapFile[], char strMaskPath[],
			char strGenomeGrp[], int nIncludeNonGrp)
{
	/* define */
	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));
	double dLog2 = log(2.0);

	struct tagBARData *pBARPos = NULL;

	FILE *fpIn;
	FILE *fpOut;

	/* variables */
	char strFileType[9];
	float fVersion;
	unsigned int nSeqNum;

	/* seq info */
	unsigned int nSeqNameLen;

	struct INTMATRIX *vProbeMappingType;
	unsigned int nProbeMappingType;
	
	struct INTMATRIX *vSequenceFileOffset;
	unsigned int nSequenceFileOffset;
	
	struct INTMATRIX *vProbePairNum;
	unsigned int nProbePairNum;
	
	unsigned int nParamNum;
		
	/* map info */
	struct INTMATRIX *vSeqID;
	unsigned int nSeqID;
	int nProbeNum;
	int nInfoCol;
	int nTotalProbeNum = 0;
	int nMaskedProbeNum = 0;
	struct tagAffyBpMapUnit *pNewUnit,*pCUnit,*pPUnit;
	struct tagAffyBpMapUnit *pUnitList;
	struct DOUBLEMATRIX *pResizeMat;
	int nIgnore = 0;
	int nDivide1M = 0;

	/* count */
	unsigned char cA,cC,cG,cT;
	unsigned char bChar;
	int ni,nj,nk,nx,ny,nz,nu;
	int nEndPos;
	struct BYTEMATRIX *pOB,*pResizeB;
	struct INTMATRIX *pOI,*pResizeI;

	/* load */
	fpIn = NULL;
	fpIn = fopen(strBpmapFile, "rb");
	if(fpIn == NULL)
	{
		printf("Error: TileProbe_BpmapToBAR, cannot open .bpmap file!\n");
		exit(EXIT_FAILURE);
	}
	fpOut = NULL;
	fpOut = fopen(strMaskPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: TileProbe_BpmapToBAR, cannot open .refmask file!\n");
		exit(EXIT_FAILURE);
	}
    fprintf(fpOut, "chromosome\tposition\tprobe_num\trepeat_num\tprobe_seq\tPMx\tPMy\tMMx\tMMy\n");
	
	/* load head */
	fread( strFileType, sizeof(char), 8, fpIn );
	strFileType[8] = '\0';
	AFFYBAR_READ_VERSION(fpIn, &fVersion);
	if(big_endian_fread(&nSeqNum, 4, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: TileProbe_BpmapToBAR, cannot load sequence number.\n");
		exit(EXIT_FAILURE);
	}
	printf("  BPMAP Version = %f\n", fVersion);
	printf("  Number of Sequences = %d\n", nSeqNum);
	if(nSeqNum <= 0)
	{
		printf("Warning: empty sequences!\n");
		return NULL;
	}

	/* create BAR object */
	pBARPos = Affy_BARData_Create();
	if(pBARPos == NULL)
	{
		printf("Error: TileProbe_BpmapToBAR, BARData object was not created!\n");
		exit(EXIT_FAILURE);
	}
	strcpy(pBARPos->strMagicnumber, "barr\r\n\032\n");
	pBARPos->fVersionnumber = 2.0;
    pBARPos->nSeqNum = nSeqNum;
	/* columns: 1 for coord, 1 for PMX, 1 for PMY, 1 for logPM value, 75/8 for A,C,G, 1 for nT, 4 for (A,C,G,T)^2, 1 for copy number */
	pBARPos->nColNum = 4+10+1+4+1;
	nInfoCol = 4;
	pBARPos->pFieldType = CreateIntMatrix(1, pBARPos->nColNum);
	if(pBARPos->pFieldType == NULL)
	{
		printf("Error: TileProbe_BpmapToBAR, cannot allocate memory for saving field type.\n");
		exit(EXIT_FAILURE);
	}
	pBARPos->pFieldType->pMatElement[0] = 2;
	pBARPos->pFieldType->pMatElement[1] = 2;
	pBARPos->pFieldType->pMatElement[2] = 2;
	pBARPos->pFieldType->pMatElement[3] = 1;
	for(ni=4; ni<15; ni++)
		pBARPos->pFieldType->pMatElement[ni] = 7;
	for(; ni<19; ni++)
		pBARPos->pFieldType->pMatElement[ni] = 2;
	pBARPos->pFieldType->pMatElement[ni] = 1;

	pBARPos->vSeqData = (struct tagBARSeq **)calloc(pBARPos->nSeqNum, sizeof(struct tagBARSeq *));
	if(pBARPos->vSeqData == NULL)
	{
		printf("Error: TileProbe_BpmapToBAR, cannot allocate memory for loading sequence data.\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<pBARPos->nSeqNum; ni++)
	{
		/* create BARSeq object */
		pBARPos->vSeqData[ni] = Affy_BARSeq_Create();
		if(pBARPos->vSeqData[ni] == NULL)
		{
			printf("Error: TileProbe_BpmapToBAR, cannot create BARSeq object.\n");
			exit(EXIT_FAILURE);
		}
	}

	vProbeMappingType = NULL;
	vProbeMappingType = CreateIntMatrix(nSeqNum,1);
	if(vProbeMappingType == NULL)
	{
		printf("Error: TileProbe_BpmapToBAR, cannot load probe mapping type!\n");
		exit(EXIT_FAILURE);
	}

	if(fVersion > 2.5)
	{
		vSequenceFileOffset = NULL;
		vSequenceFileOffset = CreateIntMatrix(nSeqNum,1);
		if(vSequenceFileOffset == NULL)
		{
			printf("Error: TileProbe_BpmapToBAR, cannot load sequence file offset!\n");
			exit(EXIT_FAILURE);
		}
	}

	vProbePairNum = NULL;
	vProbePairNum = CreateIntMatrix(nSeqNum,1);
	if(vProbePairNum == NULL)
	{
		printf("Error: TileProbe_BpmapToBAR, cannot load probe pair number!\n");
		exit(EXIT_FAILURE);
	}
	
	for(ni=0; ni<(int)nSeqNum; ni++)
	{
		/* for version 1.0 or later */
		/* load sequence name */
		if(big_endian_fread(&nSeqNameLen, 4, 1, fpIn, little_endian_machine) != 1)
		{
			printf("Error: TileProbe_BpmapToBAR, cannot load sequence name length.\n");
			exit(EXIT_FAILURE);
		}
		if(nSeqNameLen > 0)
		{
			pBARPos->vSeqData[ni]->pSeqName = CreateString(nSeqNameLen);
			if(pBARPos->vSeqData[ni]->pSeqName == NULL)
			{
				printf("Error: TileProbe_BpmapToBAR, cannot load sequence name.\n");
				exit(EXIT_FAILURE);
			}
			if(big_endian_fread(pBARPos->vSeqData[ni]->pSeqName->m_pString, 1, nSeqNameLen, fpIn, little_endian_machine) != nSeqNameLen)
			{
				printf("Error: TileProbe_BpmapToBAR, cannot load sequence name.\n");
				exit(EXIT_FAILURE);
			}
			pBARPos->vSeqData[ni]->pSeqName->m_pString[nSeqNameLen] = '\0';
		}

		/* for version 3.0 or later */
		/* load probe mapping type and sequence file offset */
		if(fVersion > 2.5)
		{
			if(big_endian_fread(&nProbeMappingType, 4, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: TileProbe_BpmapToBAR, cannot load probe mapping type.\n");
				exit(EXIT_FAILURE);
			}
			vProbeMappingType->pMatElement[ni] = (int)nProbeMappingType;
			if(big_endian_fread(&nSequenceFileOffset, 4, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: TileProbe_BpmapToBAR, cannot load sequence file offset.\n");
				exit(EXIT_FAILURE);
			}
			vSequenceFileOffset->pMatElement[ni] = (int)nSequenceFileOffset;
		}

		/* for version 1.0 or later */
		if(big_endian_fread(&nProbePairNum, 4, 1, fpIn, little_endian_machine) != 1)
		{
			printf("Error: TileProbe_BpmapToBAR, cannot load probe pair number.\n");
			exit(EXIT_FAILURE);
		}
		vProbePairNum->pMatElement[ni] = (int)nProbePairNum;
		nTotalProbeNum += (int)nProbePairNum;

		/* for version 2.0 or later */
		if(fVersion > 1.5)
		{
			/* read group name */
			if(big_endian_fread(&nSeqNameLen, 4, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: TileProbe_BpmapToBAR, cannot load group name length.\n");
				exit(EXIT_FAILURE);
			}
			if(nSeqNameLen > 0)
			{
				pBARPos->vSeqData[ni]->pSeqGroupName = CreateString(nSeqNameLen);
				if(pBARPos->vSeqData[ni]->pSeqGroupName == NULL)
				{
					printf("Error: TileProbe_BpmapToBAR, cannot load group name.\n");
					exit(EXIT_FAILURE);
				}
				if(big_endian_fread(pBARPos->vSeqData[ni]->pSeqGroupName->m_pString, 1, nSeqNameLen, fpIn, little_endian_machine) != nSeqNameLen)
				{
					printf("Error: TileProbe_BpmapToBAR, cannot load group name.\n");
					exit(EXIT_FAILURE);
				}
				pBARPos->vSeqData[ni]->pSeqGroupName->m_pString[nSeqNameLen] = '\0';
			}

			/* read version name */
			if(big_endian_fread(&nSeqNameLen, 4, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: TileProbe_BpmapToBAR, cannot load version number length.\n");
				exit(EXIT_FAILURE);
			}
			if(nSeqNameLen > 0)
			{
				pBARPos->vSeqData[ni]->pSeqVersion = CreateString(nSeqNameLen);
				if(pBARPos->vSeqData[ni]->pSeqVersion == NULL)
				{
					printf("Error: TileProbe_BpmapToBAR, cannot load version number.\n");
					exit(EXIT_FAILURE);
				}
				if(big_endian_fread(pBARPos->vSeqData[ni]->pSeqVersion->m_pString, 1, nSeqNameLen, fpIn, little_endian_machine) != nSeqNameLen)
				{
					printf("Error: TileProbe_BpmapToBAR, cannot load version number.\n");
					exit(EXIT_FAILURE);
				}
				pBARPos->vSeqData[ni]->pSeqVersion->m_pString[nSeqNameLen] = '\0';
			}

			/* read paramters */
			if(big_endian_fread(&nParamNum, 4, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: TileProbe_BpmapToBAR, cannot load number of parameters.\n");
				exit(EXIT_FAILURE);
			}
			pBARPos->vSeqData[ni]->nParamNum = (int)nParamNum;
			if(nParamNum > 0)
			{
				pBARPos->vSeqData[ni]->vParamName = (struct tagString **)calloc(pBARPos->vSeqData[ni]->nParamNum, sizeof(struct tagString *));
				pBARPos->vSeqData[ni]->vParamValue = (struct tagString **)calloc(pBARPos->vSeqData[ni]->nParamNum, sizeof(struct tagString *));

				if( (pBARPos->vSeqData[ni]->vParamName == NULL) || (pBARPos->vSeqData[ni]->vParamValue == NULL) )
				{
					printf("Error: TileProbe_BpmapToBAR, cannot allocate memory for loading parameter name/value pairs\n");
					exit(EXIT_FAILURE);
				}
			}
			
			for(nj=0; nj<(int)nParamNum; nj++)
			{
				if(big_endian_fread(&nSeqNameLen, 4, 1, fpIn, little_endian_machine) != 1)
				{
					printf("Error: TileProbe_BpmapToBAR, cannot load parameter name length.\n");
					exit(EXIT_FAILURE);
				}
				if(nSeqNameLen > 0)
				{
					pBARPos->vSeqData[ni]->vParamName[nj] = CreateString(nSeqNameLen);
					if(pBARPos->vSeqData[ni]->vParamName[nj] == NULL)
					{
						printf("Error: TileProbe_BpmapToBAR, cannot load parameter name.\n");
						exit(EXIT_FAILURE);
					}
					if(big_endian_fread(pBARPos->vSeqData[ni]->vParamName[nj]->m_pString, 1, nSeqNameLen, fpIn, little_endian_machine) != nSeqNameLen)
					{
						printf("Error: TileProbe_BpmapToBAR, cannot load parameter name.\n");
						exit(EXIT_FAILURE);
					}
					pBARPos->vSeqData[ni]->vParamName[nj]->m_pString[nSeqNameLen] = '\0';
				}

				if(big_endian_fread(&nSeqNameLen, 4, 1, fpIn, little_endian_machine) != 1)
				{
					printf("Error: TileProbe_BpmapToBAR, cannot load parameter value length.\n");
					exit(EXIT_FAILURE);
				}
				if(nSeqNameLen > 0)
				{
					pBARPos->vSeqData[ni]->vParamValue[nj] = CreateString(nSeqNameLen);
					if(pBARPos->vSeqData[ni]->vParamValue[nj] == NULL)
					{
						printf("Error: TileProbe_BpmapToBAR, cannot load parameter value.\n");
						exit(EXIT_FAILURE);
					}
					if(big_endian_fread(pBARPos->vSeqData[ni]->vParamValue[nj]->m_pString, 1, nSeqNameLen, fpIn, little_endian_machine) != nSeqNameLen)
					{
						printf("Error: TileProbe_BpmapToBAR, cannot load parameter value.\n");
						exit(EXIT_FAILURE);
					}
					pBARPos->vSeqData[ni]->vParamValue[nj]->m_pString[nSeqNameLen] = '\0';
				}
			}
		}
	}

	/* load map */
	vSeqID = NULL;
	vSeqID = CreateIntMatrix(nSeqNum,1);
	if(vSeqID == NULL)
	{
		printf("Error: TileProbe_BpmapToBAR, cannot load sequence id!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<(int)nSeqNum; ni++)
	{
		nDivide1M = 0;
		if(strcmp(strGenomeGrp, "") != 0)
		{
			if(pBARPos->vSeqData[ni]->pSeqGroupName == NULL)
			{
				nDivide1M = 1;
			}
			else
			{
				if(strcmp(strGenomeGrp, pBARPos->vSeqData[ni]->pSeqGroupName->m_pString) != 0)
					nDivide1M = 1;
			}
		}

		/* load seq id */
		if(big_endian_fread(&nSeqID, 4, 1, fpIn, little_endian_machine) != 1)
		{
			printf("Error: TileProbe_BpmapToBAR, cannot load sequence ID.\n");
			exit(EXIT_FAILURE);
		}
		vSeqID->pMatElement[ni] = nSeqID;

		/* create initial memory, 1 for logPM value, 75/8 for A,C,G, 1 for nT, 4 for (A,C,G,T)^2, 1 for copy number */
		pBARPos->vSeqData[ni]->nDataNum = 0;
		pBARPos->vSeqData[ni]->nColNum = pBARPos->nColNum;
		/* 
		pBARPos->vSeqData[ni]->nColNum = 4;
		if(vProbeMappingType->pMatElement[ni] == 1)
		{
			pBARPos->vSeqData[ni]->nColNum = 4;
		}
		else
		{
			pBARPos->vSeqData[ni]->nColNum = 6;
		}
		nInfoCol = pBARPos->vSeqData[ni]->nColNum;		
		pBARPos->vSeqData[ni]->nColNum = pBARPos->vSeqData[ni]->nColNum+10+1+4;
		*/

		pBARPos->vSeqData[ni]->vData = (struct DOUBLEMATRIX **)calloc(pBARPos->vSeqData[ni]->nColNum, sizeof(struct DOUBLEMATRIX *));
		if(pBARPos->vSeqData[ni]->vData == NULL)
		{
			printf("Error: TileProbe_BpmapToBAR, cannot allocate memory for loading bpmap coordinate data!\n");
			exit(EXIT_FAILURE);
		}
		nProbeNum = vProbePairNum->pMatElement[ni];
		if(nProbeNum > 0)
		{
			for(nj=0; nj<nInfoCol; nj++)
			{
				pBARPos->vSeqData[ni]->vData[nj] = CreateDoubleMatrix(1, nProbeNum);
				if(pBARPos->vSeqData[ni]->vData[nj] == NULL)
				{
					printf("Error: TileProbe_BpmapToBAR, cannot allocate memory for loading bpmap coordinate data!\n");
					exit(EXIT_FAILURE);
				}
			}

			for(; nj<nInfoCol+11; nj++)
			{
				pBARPos->vSeqData[ni]->vData[nj] = (struct DOUBLEMATRIX *)(CreateByteMatrix(1, nProbeNum));
				if(pBARPos->vSeqData[ni]->vData[nj] == NULL)
				{
					printf("Error: TileProbe_BpmapToBAR, cannot allocate memory for loading bpmap sequence data!\n");
					exit(EXIT_FAILURE);
				}
			}

			for(; nj<nInfoCol+15; nj++)
			{
				pBARPos->vSeqData[ni]->vData[nj] = (struct DOUBLEMATRIX *)(CreateIntMatrix(1, nProbeNum));
				if(pBARPos->vSeqData[ni]->vData[nj] == NULL)
				{
					printf("Error: TileProbe_BpmapToBAR, cannot allocate memory for loading bpmap sequence-square data!\n");
					exit(EXIT_FAILURE);
				}
			}

			pBARPos->vSeqData[ni]->vData[nj] = CreateDoubleMatrix(1, nProbeNum);
			if(pBARPos->vSeqData[ni]->vData[nj] == NULL)
			{
				printf("Error: TileProbe_BpmapToBAR, cannot allocate memory for loading copy number data!\n");
				exit(EXIT_FAILURE);
			}
		}


		pUnitList = NULL;
		
		/* load probes */
		for(nj=0; nj<nProbeNum; nj++)
		{
			/* load new unit */
			pNewUnit = NULL;
			pNewUnit = AffyBpMapUnitCreate();
			AffyBpMapUnitLoad_v3m2(pNewUnit, fpIn, vProbeMappingType->pMatElement[ni], little_endian_machine);
			if(nDivide1M == 1)
			{
				pNewUnit->fMatchScore /= 1000000;
			}
			if(pNewUnit->fMatchScore < 1e-7)
				pNewUnit->fMatchScore = 1e-6;

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

							fprintf(fpOut, "%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%f\t%d\t%f\n", 
								pBARPos->vSeqData[ni]->pSeqName->m_pString, pCUnit->nPos, pCUnit->nDepthNum, 
								pCUnit->nRepeatNum, pCUnit->strProbeSeq, 
								pCUnit->nPMX, pCUnit->nPMY, pCUnit->nMMX, pCUnit->nMMY,
								pCUnit->fMatchScore, pCUnit->bStrand);

							if(pCUnit->nRepeatNum < 2)
							{
								nk = pBARPos->vSeqData[ni]->nDataNum;
								pBARPos->vSeqData[ni]->vData[0]->pMatElement[nk] = pCUnit->nPos;
								pBARPos->vSeqData[ni]->vData[1]->pMatElement[nk] = pCUnit->nPMX;
								pBARPos->vSeqData[ni]->vData[2]->pMatElement[nk] = pCUnit->nPMY;
								/* 
								nInfoCol = 4;
								if(vProbeMappingType->pMatElement[ni] == 0)
								{
									pBARPos->vSeqData[ni]->vData[3]->pMatElement[nk] = pCUnit->nMMX;
									pBARPos->vSeqData[ni]->vData[4]->pMatElement[nk] = pCUnit->nMMY;
									nInfoCol = 6;
								}
								*/

								cA = 0;
								cC = 0;
								cG = 0;
								cT = 0;

								ny = 0;
								for(nx=0; nx<25; nx++)
								{
									if(pCUnit->strProbeSeq[nx] == 'A')
									{
										nz = ny/8;
										nu = 7-ny%8;
										bChar = 0x01;
										bChar <<=nu;
										pOB = (struct BYTEMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+nz]);
										pOB->pMatElement[nk] |= bChar;
										cA++;
									}
									else if(pCUnit->strProbeSeq[nx] == 'C')
									{
										nz = (ny+1)/8;
										nu = 7-(ny+1)%8;
										bChar = 0x01;
										bChar <<=nu;
										pOB = (struct BYTEMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+nz]);
										pOB->pMatElement[nk] |= bChar;
										cC++;
									}
									else if(pCUnit->strProbeSeq[nx] == 'G')
									{
										nz = (ny+2)/8;
										nu = 7-(ny+2)%8;
										bChar = 0x01;
										bChar <<=nu;
										pOB = (struct BYTEMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+nz]);
										pOB->pMatElement[nk] |= bChar;
										cG++;
									}
									else if(pCUnit->strProbeSeq[nx] == 'T')
									{
										cT++;
									}
									else
									{
										printf("Error: TileProbe_BpmapToBAR, cannot recognize the base in the probe sequence. \n");
										exit(EXIT_FAILURE);
									}

									ny += 3;
								}

								pOB = (struct BYTEMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+10]);
								pOB->pMatElement[nk] = cT;
								pOI = (struct INTMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+11]);
								pOI->pMatElement[nk] = (int)(cA*cA);
								pOI = (struct INTMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+12]);
								pOI->pMatElement[nk] = (int)(cC*cC);
								pOI = (struct INTMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+13]);
								pOI->pMatElement[nk] = (int)(cG*cG);
								pOI = (struct INTMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+14]);
								pOI->pMatElement[nk] = (int)(cT*cT);
								pBARPos->vSeqData[ni]->vData[nInfoCol+15]->pMatElement[nk] = log(pCUnit->fMatchScore)/dLog2;

								pBARPos->vSeqData[ni]->nDataNum += 1;
							}

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

			fprintf(fpOut, "%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%f\t%d\n", 
						pBARPos->vSeqData[ni]->pSeqName->m_pString, pCUnit->nPos, pCUnit->nDepthNum, 
						pCUnit->nRepeatNum, pCUnit->strProbeSeq, 
						pCUnit->nPMX, pCUnit->nPMY, pCUnit->nMMX, pCUnit->nMMY,
						pCUnit->fMatchScore, pCUnit->bStrand);

			if(pCUnit->nRepeatNum < 2)
			{
				nk = pBARPos->vSeqData[ni]->nDataNum;
				pBARPos->vSeqData[ni]->vData[0]->pMatElement[nk] = pCUnit->nPos;
				pBARPos->vSeqData[ni]->vData[1]->pMatElement[nk] = pCUnit->nPMX;
				pBARPos->vSeqData[ni]->vData[2]->pMatElement[nk] = pCUnit->nPMY;
				
				/* 
				nInfoCol = 4;
				if(vProbeMappingType->pMatElement[ni] == 0)
				{
					pBARPos->vSeqData[ni]->vData[3]->pMatElement[nk] = pCUnit->nMMX;
					pBARPos->vSeqData[ni]->vData[4]->pMatElement[nk] = pCUnit->nMMY;
					nInfoCol = 6;
				}
				*/

				cA = 0;
				cC = 0;
				cG = 0;
				cT = 0;

				ny = 0;
				for(nx=0; nx<25; nx++)
				{
					if(pCUnit->strProbeSeq[nx] == 'A')
					{
						nz = ny/8;
						nu = 7-ny%8;
						bChar = 0x01;
						bChar <<=nu;
						pOB = (struct BYTEMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+nz]);
						pOB->pMatElement[nk] |= bChar;
						cA++;
					}
					else if(pCUnit->strProbeSeq[nx] == 'C')
					{
						nz = (ny+1)/8;
						nu = 7-(ny+1)%8;
						bChar = 0x01;
						bChar <<=nu;
						pOB = (struct BYTEMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+nz]);
						pOB->pMatElement[nk] |= bChar;
						cC++;
					}
					else if(pCUnit->strProbeSeq[nx] == 'G')
					{
						nz = (ny+2)/8;
						nu = 7-(ny+2)%8;
						bChar = 0x01;
						bChar <<=nu;
						pOB = (struct BYTEMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+nz]);
						pOB->pMatElement[nk] |= bChar;
						cG++;
					}
					else if(pCUnit->strProbeSeq[nx] == 'T')
					{
						cT++;
					}
					else
					{
						printf("Error: TileProbe_BpmapToBAR, cannot recognize the base in the probe sequence. \n");
						exit(EXIT_FAILURE);
					}

					ny += 3;
				}

				pOB = (struct BYTEMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+10]);
				pOB->pMatElement[nk] = cT;
				pOI = (struct INTMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+11]);
				pOI->pMatElement[nk] = (int)(cA*cA);
				pOI = (struct INTMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+12]);
				pOI->pMatElement[nk] = (int)(cC*cC);
				pOI = (struct INTMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+13]);
				pOI->pMatElement[nk] = (int)(cG*cG);
				pOI = (struct INTMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+14]);
				pOI->pMatElement[nk] = (int)(cT*cT);
				pBARPos->vSeqData[ni]->vData[nInfoCol+15]->pMatElement[nk] = log(pCUnit->fMatchScore)/dLog2;


				pBARPos->vSeqData[ni]->nDataNum += 1;
			}

			AffyBpMapUnitDestroy(pCUnit);
		}

		/* 
		nInfoCol = 4;
		if(vProbeMappingType->pMatElement[ni] == 0)
		{
			pBARPos->vSeqData[ni]->vData[3]->pMatElement[nk] = pCUnit->nMMX;
			pBARPos->vSeqData[ni]->vData[4]->pMatElement[nk] = pCUnit->nMMY;
			nInfoCol = 6;
		}
		*/

		for(nj=0; nj<nInfoCol; nj++)
		{
			if(pBARPos->vSeqData[ni]->nDataNum > 0)
			{
				pResizeMat = NULL;
				pResizeMat = CreateDoubleMatrix(1, pBARPos->vSeqData[ni]->nDataNum);
				if(pResizeMat == NULL)
				{
					printf("Error: TileProbe_BpmapToBAR, cannot create an intermediate matrix to move data!\n");
					exit(EXIT_FAILURE);
				}
				memcpy(pResizeMat->pMatElement, pBARPos->vSeqData[ni]->vData[nj]->pMatElement, pBARPos->vSeqData[ni]->nDataNum*sizeof(double));
				DestroyDoubleMatrix(pBARPos->vSeqData[ni]->vData[nj]);
				pBARPos->vSeqData[ni]->vData[nj] = pResizeMat;
			}
			else
			{
				DestroyDoubleMatrix(pBARPos->vSeqData[ni]->vData[nj]);
				pBARPos->vSeqData[ni]->vData[nj] = NULL;
			}
		}

		for(; nj<nInfoCol+11; nj++)
		{
			if(pBARPos->vSeqData[ni]->nDataNum > 0)
			{
				pResizeB = NULL;
				pResizeB = CreateByteMatrix(1, pBARPos->vSeqData[ni]->nDataNum);
				if(pResizeB == NULL)
				{
					printf("Error: TileProbe_BpmapToBAR, cannot create an intermediate matrix to move data!\n");
					exit(EXIT_FAILURE);
				}
				pOB = (struct BYTEMATRIX *)(pBARPos->vSeqData[ni]->vData[nj]);
				memcpy(pResizeB->pMatElement, pOB->pMatElement, pBARPos->vSeqData[ni]->nDataNum*sizeof(unsigned char));
				DestroyByteMatrix(pOB);
				pBARPos->vSeqData[ni]->vData[nj] = (struct DOUBLEMATRIX *)pResizeB;
			}
			else
			{
				pOB = (struct BYTEMATRIX *)(pBARPos->vSeqData[ni]->vData[nj]);
				DestroyByteMatrix(pOB);
				pBARPos->vSeqData[ni]->vData[nj] = NULL;
			}
		}

		for(; nj<nInfoCol+15; nj++)
		{
			if(pBARPos->vSeqData[ni]->nDataNum > 0)
			{
				pResizeI = NULL;
				pResizeI = CreateIntMatrix(1, pBARPos->vSeqData[ni]->nDataNum);
				if(pResizeI == NULL)
				{
					printf("Error: TileProbe_BpmapToBAR, cannot create an intermediate matrix to move data!\n");
					exit(EXIT_FAILURE);
				}
				pOI = (struct INTMATRIX *)(pBARPos->vSeqData[ni]->vData[nj]);
				memcpy(pResizeI->pMatElement, pOI->pMatElement, pBARPos->vSeqData[ni]->nDataNum*sizeof(int));
				DestroyIntMatrix(pOI);
				pBARPos->vSeqData[ni]->vData[nj] = (struct DOUBLEMATRIX *)pResizeI;
			}
			else
			{
				pOI = (struct INTMATRIX *)(pBARPos->vSeqData[ni]->vData[nj]);
				DestroyIntMatrix(pOI);
				pBARPos->vSeqData[ni]->vData[nj] = NULL;
			}
		}

		if(pBARPos->vSeqData[ni]->nDataNum > 0)
		{
			pResizeMat = NULL;
			pResizeMat = CreateDoubleMatrix(1, pBARPos->vSeqData[ni]->nDataNum);
			if(pResizeMat == NULL)
			{
				printf("Error: TileProbe_BpmapToBAR, cannot create an intermediate matrix to move data!\n");
				exit(EXIT_FAILURE);
			}
			memcpy(pResizeMat->pMatElement, pBARPos->vSeqData[ni]->vData[nj]->pMatElement, pBARPos->vSeqData[ni]->nDataNum*sizeof(double));
			DestroyDoubleMatrix(pBARPos->vSeqData[ni]->vData[nj]);
			pBARPos->vSeqData[ni]->vData[nj] = pResizeMat;
		}
		else
		{
			DestroyDoubleMatrix(pBARPos->vSeqData[ni]->vData[nj]);
			pBARPos->vSeqData[ni]->vData[nj] = NULL;
		}

		/* update masked probe number */
		nMaskedProbeNum += pBARPos->vSeqData[ni]->nDataNum;
	}

	/* load tail if any */
	printf("  ");
	while(feof(fpIn) == 0)
	{
		fread(strFileType, sizeof(char), 1, fpIn);
		printf("%c", strFileType[0]);
	}
	printf("\n");

	/* clear memeory */
	DestroyIntMatrix(vProbeMappingType);
	if(fVersion > 2.5)
	{
		DestroyIntMatrix(vSequenceFileOffset);
	}
    DestroyIntMatrix(vProbePairNum);
	DestroyIntMatrix(vSeqID);

	/* close file */
	fclose(fpIn);
	fclose(fpOut);

	printf("  Probe number before masking = %d\n", nTotalProbeNum);
	printf("  Probe number after masking = %d\n", nMaskedProbeNum);

	/* return */
	return pBARPos;
}

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_MATv2_XX()                                                   */
/*  Obtain X'X.                                                            */
/* ----------------------------------------------------------------------- */
struct DOUBLEMATRIX *TileProbe_MATv2_XX(struct tagBARData *pBARPos, 
			char strGenomeGrp[], int nIncludeNonGrp)
{
	/* define */
	struct DOUBLEMATRIX *pXX = NULL;
	int nS = 75+1+4+1;
	int nBN = 75;
	int nBWN = (nBN/8)+1;
	int nCN = nS-nBN-1;
	int nInfoCol = 4;
	int ni,nj,nk,nx,ny,nz,nu;
	int nFinish;
	/* unsigned char *vB; */
	int vBI[25];
	int *vC;
	double dMatchScore;
	struct BYTEMATRIX *pBM;
	struct INTMATRIX *pIM;
	double dTemp;
	/* unsigned char cTemp; */

	/* init */
	if(pBARPos == NULL)
		return NULL;

	pXX = CreateDoubleMatrix(nS,nS);
	if(pXX == NULL)
	{
		printf("Error: TileProbe_XX, cannot allocate memory for X'X\n");
		exit(EXIT_FAILURE);
	}

	/* vB = NULL;
	vB = (unsigned char *)calloc(nBN, sizeof(unsigned char));
	if(vB == NULL)
	{
		printf("Error: TileProbe_XX, cannot allocate memory for vB\n");
		exit(EXIT_FAILURE);
	} */

	vC = NULL;
	vC = (int *)calloc(nCN, sizeof(int));
	if(vC == NULL)
	{
		printf("Error: TileProbe_XX, cannot allocate memory for vC\n");
		exit(EXIT_FAILURE);
	}
	dMatchScore = 0.0;

	/* compute */
	for(ni=0; ni<pBARPos->nSeqNum; ni++)
	{
		if(pBARPos->vSeqData[ni]->nDataNum <= 0)
			continue;

		if(nIncludeNonGrp == 0)
		{
			if(strcmp(strGenomeGrp, "") != 0)
			{
				if(pBARPos->vSeqData[ni]->pSeqGroupName == NULL)
					continue;

				if(strcmp(strGenomeGrp, pBARPos->vSeqData[ni]->pSeqGroupName->m_pString) != 0)
					continue;
			}
		}

		for(nj=0; nj<pBARPos->vSeqData[ni]->nDataNum; nj++)
		{
			/* get vB */
			nx = 0;
			nu = 0;
			nFinish = 0;
			for(nk=0; nk<nBWN; nk++)
			{
				pBM = (struct BYTEMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+nk]);
				for(ny=0; ny<8; ny++)
				{
					nz = 7-ny;
					
					/* vB[nx] = (pBM->pMatElement[nj] >> nz) & 0x01; */
					if( ((pBM->pMatElement[nj] >> nz) & 0x01) != 0 )
					{
						vBI[nu] = nx;
						nu++;
					}

					nx++;
					if(nx >= nBN)
					{
						nFinish = 1;
						break;
					}
				}

				if(nFinish == 1)
					break;
			}

			/* get vC */
			pBM = (struct BYTEMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+nBWN]);
			vC[0] = pBM->pMatElement[nj];
			for(nk=1; nk<nCN; nk++)
			{
				pIM = (struct INTMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+nBWN+nk]);
				vC[nk] = pIM->pMatElement[nj];
			}
			dMatchScore = pBARPos->vSeqData[ni]->vData[nInfoCol+nBWN+nk]->pMatElement[nj];

			if((nu+vC[0]) != 25)
			{
				printf("Error: TileProbe_XX, probe sequence matching wrong!\n");
				exit(EXIT_FAILURE);
			}

			/* get X'X */
			/* for(nx=0; nx<nBN; nx++)
			{
				for(ny=nx; ny<nBN; ny++)
				{
					cTemp = vB[nx] & vB[ny];
					dTemp = cTemp+DMGETAT(pXX, nx, ny);
					DMSETAT(pXX, nx, ny, dTemp);
				}

				for(ny=0; ny<nCN; ny++)
				{
					dTemp = vB[nx]*vC[ny]+DMGETAT(pXX, nx, nBN+ny);
					DMSETAT(pXX, nx, nBN+ny, dTemp);
				}

			} */

			for(nx=0; nx<nu; nx++)
			{
				for(ny=nx; ny<nu; ny++)
				{
					dTemp = 1+DMGETAT(pXX, vBI[nx], vBI[ny]);
					DMSETAT(pXX, vBI[nx], vBI[ny], dTemp);
				}

				for(ny=0; ny<nCN; ny++)
				{
					dTemp = vC[ny]+DMGETAT(pXX, vBI[nx], nBN+ny);
					DMSETAT(pXX, vBI[nx], nBN+ny, dTemp);
				}

				dTemp = dMatchScore+DMGETAT(pXX, vBI[nx], nBN+nCN);
				DMSETAT(pXX, vBI[nx], nBN+nCN, dTemp);
			}

			for(nx=0; nx<nCN; nx++)
			{
				for(ny=nx; ny<nCN; ny++)
				{
					dTemp = vC[nx]*vC[ny]+DMGETAT(pXX, nBN+nx, nBN+ny);
					DMSETAT(pXX, nBN+nx, nBN+ny, dTemp);
				}

				dTemp = vC[nx]*dMatchScore+DMGETAT(pXX, nBN+nx, nBN+nCN);
				DMSETAT(pXX, nBN+nx, nBN+nCN, dTemp);
			}

			dTemp = dMatchScore*dMatchScore+DMGETAT(pXX, nBN+nCN, nBN+nCN);
			DMSETAT(pXX, nBN+nCN, nBN+nCN, dTemp);
		}
	}

	/* make X'X symmetrix */
	/* DMSAVE(pXX, "G:\\Projects\\TileProbe_Project\\mattest\\testXX.txt"); */
	for(nx=0; nx<nS; nx++)
	{
		for(ny=0; ny<nx; ny++)
		{
			dTemp = DMGETAT(pXX, ny, nx);
			DMSETAT(pXX, nx, ny, dTemp);
		}
	}
	/* DMSAVE(pXX, "G:\\Projects\\TileProbe_Project\\mattest\\testXXsym.txt"); */


	/* clear memory */
	/* free(vB); */
	free(vC);

	/* return */
	return pXX;
}

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_MATv2_XY()                                                   */
/*  Obtain X'Y.                                                            */
/* ----------------------------------------------------------------------- */
struct DOUBLEMATRIX *TileProbe_MATv2_XY(struct tagBARData *pBARPos, struct tagCELData *pCELData, 
				char strGenomeGrp[], int nIncludeNonGrp)
{
	/* define */
	double dNormLowerBound = 1.0;
	double dLog2 = log(2.0);
	struct DOUBLEMATRIX *pXY = NULL;
	int nS = 75+1+4+1;
	int nBN = 75;
	int nBWN = (nBN/8)+1;
	int nCN = nS-nBN-1;
	int nInfoCol = 4;
	int ni,nj,nk,nx,ny,nz,nu,nidx1;
	int nFinish;
	/* unsigned char *vB; */ 
	int vBI[25];
	int *vC;
	double dMatchScore;
	struct BYTEMATRIX *pBM;
	struct INTMATRIX *pIM;
	int nCELNumberCells;
	int nCELTotalX,nCELTotalY;
	double dPM;

	/* init */
	if(pBARPos == NULL)
		return NULL;

	nCELNumberCells = pCELData->nNumberCells;
	nCELTotalX = pCELData->nCols;
	nCELTotalY = pCELData->nRows;
	if(nCELTotalX*nCELTotalY != nCELNumberCells)
	{
		printf("Error: TileProbe_MAT_XY, CEL file Cols*Rows != NumberCells!\n");
		exit(EXIT_FAILURE);
	}

	pXY = CreateDoubleMatrix(nS,1);
	if(pXY == NULL)
	{
		printf("Error: TileProbe_XY, cannot allocate memory for X'Y\n");
		exit(EXIT_FAILURE);
	}

	/* vB = NULL;
	vB = (unsigned char *)calloc(nBN, sizeof(unsigned char));
	if(vB == NULL)
	{
		printf("Error: TileProbe_XY, cannot allocate memory for vB\n");
		exit(EXIT_FAILURE);
	} */

	vC = NULL;
	vC = (int *)calloc(nCN, sizeof(int));
	if(vC == NULL)
	{
		printf("Error: TileProbe_XY, cannot allocate memory for vC\n");
		exit(EXIT_FAILURE);
	}
	dMatchScore = 0.0;

	/* compute */
	for(ni=0; ni<pBARPos->nSeqNum; ni++)
	{
		if(pBARPos->vSeqData[ni]->nDataNum <= 0)
			continue;

		if(nIncludeNonGrp == 0)
		{
			if(strcmp(strGenomeGrp, "") != 0)
			{
				if(pBARPos->vSeqData[ni]->pSeqGroupName == NULL)
					continue;

				if(strcmp(strGenomeGrp, pBARPos->vSeqData[ni]->pSeqGroupName->m_pString) != 0)
					continue;
			}
		}

		for(nj=0; nj<pBARPos->vSeqData[ni]->nDataNum; nj++)
		{
			/* get intensity */
			nidx1 = (int)(pBARPos->vSeqData[ni]->vData[2]->pMatElement[nj])*nCELTotalX+(int)(pBARPos->vSeqData[ni]->vData[1]->pMatElement[nj]);
			if(nidx1 >= nCELNumberCells)
			{
				printf("Error: TileProbe_MAT_XY, index out of range\n");
				exit(EXIT_FAILURE);
			}

			dPM = pCELData->pIntensity->pMatElement[nidx1];

			if(dPM < dNormLowerBound)
				dPM = dNormLowerBound;

			dPM = log(dPM)/dLog2;
			pBARPos->vSeqData[ni]->vData[3]->pMatElement[nj] = dPM;

			/* get vB */
			nx = 0;
			nu = 0;
			nFinish = 0;
			for(nk=0; nk<nBWN; nk++)
			{
				pBM = (struct BYTEMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+nk]);
				for(ny=0; ny<8; ny++)
				{
					nz = 7-ny;
					/* vB[nx] = (pBM->pMatElement[nj] >> nz) & 0x01; */
					if( ((pBM->pMatElement[nj] >> nz) & 0x01) != 0)
					{
						vBI[nu] = nx;
						nu++;
					}

					nx++;
					if(nx >= nBN)
					{
						nFinish = 1;
						break;
					}
				}

				if(nFinish == 1)
					break;
			}

			/* get vC */
			pBM = (struct BYTEMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+nBWN]);
			vC[0] = pBM->pMatElement[nj];
			for(nk=1; nk<nCN; nk++)
			{
				pIM = (struct INTMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+nBWN+nk]);
				vC[nk] = pIM->pMatElement[nj];
			}
			dMatchScore = pBARPos->vSeqData[ni]->vData[nInfoCol+nBWN+nCN]->pMatElement[nj];

			if( (nu+vC[0]) != 25)
			{
				printf("Error: TileProbe_MAT_XY, probe sequence matching problem!\n");
				exit(EXIT_FAILURE);
			}

			/* get X'Y */
			for(nx=0; nx<nu; nx++)
			{
				ny = vBI[nx];
				pXY->pMatElement[ny] += dPM; 
			}
			/* for(nx=0; nx<nBN; nx++)
			{
				pXY->pMatElement[nx] += vB[nx]*dPM; 
			} */
			for(nx=0; nx<nCN; nx++)
			{
				pXY->pMatElement[nBN+nx] += vC[nx]*dPM;
			}

			pXY->pMatElement[nBN+nCN] += dMatchScore*dPM;
		}
	}

	
	/* clear memory */
	/* free(vB); */
	free(vC);

	/* return */
	return pXY;
}

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_MATv2_Correction()                                           */
/*  Background correction.                                                 */
/* ----------------------------------------------------------------------- */
int TileProbe_MATv2_Correction(struct tagBARData *pBARPos, struct DOUBLEMATRIX *pBeta,
							 char strOutFile[], char strGenomeGrp[], int nIncludeNonGrp)
{
	/* define */
	int nDataTotNum = 0;
	struct DOUBLEMATRIX *pE;
	struct DOUBLEMATRIX *pV;
	struct DOUBLEMATRIX *pESort;
	struct LONGMATRIX *pESortIndex;
	int nBinSize = 3000;
	int nS = 75+1+4+1;
	int nBN = 75;
	int nBWN = (nBN/8)+1;
	int nCN = nS-nBN-1;
	int nInfoCol = 4;
	int ni,nj,nk,nx,ny,nz,nh,nu;
	int nFinish;
	/* unsigned char *vB; */
	int vBI[25];
	int *vC;
	double dMatchScore;
	struct BYTEMATRIX *pBM;
	struct INTMATRIX *pIM;
	double dTemp,dSD;
	struct INTMATRIX *pCol;

	/* init */
	if(pBARPos == NULL || pBeta == NULL)
	{
		printf("Error: TileProbe_MAT_Correction, cannot correct background for %s!\n", strOutFile);
		return PROC_FAILURE;
	}

	/* 1: Create a single matrix for MAT correction */
	for(ni=0; ni<pBARPos->nSeqNum; ni++)
	{
		if(nIncludeNonGrp == 0)
		{
			if(strcmp(strGenomeGrp, "") != 0)
			{
				if(pBARPos->vSeqData[ni]->pSeqGroupName == NULL)
					continue;

				if(strcmp(strGenomeGrp, pBARPos->vSeqData[ni]->pSeqGroupName->m_pString) != 0)
					continue;
			}
		}
		nDataTotNum += pBARPos->vSeqData[ni]->nDataNum;
	}
	pE = NULL;
	pE = CreateDoubleMatrix(1, nDataTotNum);
	if(pE == NULL)
	{
		printf("Error: TileProbe_MAT_Correction, cannot allocate memory for MAT adjustment!\n");
		exit(EXIT_FAILURE);
	}
	pV = NULL;
	pV = CreateDoubleMatrix(1, nDataTotNum);
	if(pV == NULL)
	{
		printf("Error: TileProbe_MAT_Correction, cannot allocate memory for MAT adjustment!\n");
		exit(EXIT_FAILURE);
	}

	/* 2: Sort fitted background */
	/* vB = NULL;
	vB = (unsigned char *)calloc(nBN, sizeof(unsigned char));
	if(vB == NULL)
	{
		printf("Error: TileProbe_XX, cannot allocate memory for vB\n");
		exit(EXIT_FAILURE);
	} */

	vC = NULL;
	vC = (int *)calloc(nCN, sizeof(int));
	if(vC == NULL)
	{
		printf("Error: TileProbe_MAT_Correction, cannot allocate memory for vC\n");
		exit(EXIT_FAILURE);
	}
	dMatchScore = 0.0;

	/* compute */
	nh = 0;
	for(ni=0; ni<pBARPos->nSeqNum; ni++)
	{
		if(pBARPos->vSeqData[ni]->nDataNum <= 0)
			continue;

		if(nIncludeNonGrp == 0)
		{
			if(strcmp(strGenomeGrp, "") != 0)
			{
				if(pBARPos->vSeqData[ni]->pSeqGroupName == NULL)
					continue;

				if(strcmp(strGenomeGrp, pBARPos->vSeqData[ni]->pSeqGroupName->m_pString) != 0)
					continue;
			}
		}

		for(nj=0; nj<pBARPos->vSeqData[ni]->nDataNum; nj++)
		{
			/* get vB */
			nx = 0;
			nu = 0;
			nFinish = 0;
			for(nk=0; nk<nBWN; nk++)
			{
				pBM = (struct BYTEMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+nk]);
				for(ny=0; ny<8; ny++)
				{
					nz = 7-ny;
					/* vB[nx] = (pBM->pMatElement[nj] >> nz) & 0x01; */
					if( ((pBM->pMatElement[nj] >> nz) & 0x01) != 0 )
					{
						vBI[nu] = nx;
						nu++;
					}

					nx++;
					if(nx >= nBN)
					{
						nFinish = 1;
						break;
					}
				}

				if(nFinish == 1)
					break;
			}

			/* get vC */
			pBM = (struct BYTEMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+nBWN]);
			vC[0] = pBM->pMatElement[nj];
			for(nk=1; nk<nCN; nk++)
			{
				pIM = (struct INTMATRIX *)(pBARPos->vSeqData[ni]->vData[nInfoCol+nBWN+nk]);
				vC[nk] = pIM->pMatElement[nj];
			}
			dMatchScore = pBARPos->vSeqData[ni]->vData[nInfoCol+nBWN+nCN]->pMatElement[nj];

			if( (nu+vC[0]) != 25 )
			{
				printf("Error: TileProbe_MAT_Correction, probe sequence matching problem!\n");
				exit(EXIT_FAILURE);
			}

			/* get expected value */
			dTemp = 0.0;
			for(nx=0; nx<nu; nx++)
			{
				ny = vBI[nx];
				dTemp += pBeta->pMatElement[ny];
			}
			for(nx=0; nx<nCN; nx++)
			{
				dTemp += vC[nx]*pBeta->pMatElement[nBN+nx];
			}
			dTemp += dMatchScore*pBeta->pMatElement[nBN+nCN];

			pE->pMatElement[nh] = dTemp;
			pV->pMatElement[nh] = pBARPos->vSeqData[ni]->vData[3]->pMatElement[nj];
			pBARPos->vSeqData[ni]->vData[3]->pMatElement[nj] -= dTemp;
			
			/* nh ++ */
			nh++;
		}
	}

	/* clear memory */
	/* free(vB); */
	free(vC);

	if(nh != nDataTotNum)
	{
		printf("Error: TileProbe_MAT_Correction, probe number not match!\n");
		exit(EXIT_FAILURE);
	}
	
	/* 3: Stratified variance correction */
	DMSORTMERGEA_0(pE, &pESort, &pESortIndex);
	
	ni = 0;
	
	while(ni<nDataTotNum)
	{
		nj = ni+nBinSize-1;
		if(nj >= nDataTotNum)
			nj = nDataTotNum-1;
		else if( (nDataTotNum-nj)< (nBinSize/2) )
			nj = nDataTotNum-1;

		dTemp = 0.0;
		for(nk=ni; nk<=nj; nk++)
		{
			nx = pESortIndex->pMatElement[nk];
			dTemp += pV->pMatElement[nx];
		}
		dTemp /= (nj-ni+1);
		dSD = 0.0;
		for(nk=ni; nk<=nj; nk++)
		{
			nx = pESortIndex->pMatElement[nk];
			dSD += (pV->pMatElement[nx]-dTemp)*(pV->pMatElement[nx]-dTemp);
		}
		dSD /= (nj-ni);
		dSD = sqrt(dSD);

		for(nk=ni; nk<=nj; nk++)
		{
			nx = pESortIndex->pMatElement[nk];
			pV->pMatElement[nx] = dSD;
		}

		ni = nj+1;
	}

	DestroyDoubleMatrix(pESort);
	DestroyLongMatrix(pESortIndex);

	nh = 0;
	for(ni=0; ni<pBARPos->nSeqNum; ni++)
	{
		if(pBARPos->vSeqData[ni]->nDataNum <= 0)
			continue;

		if(nIncludeNonGrp == 0)
		{
			if(strcmp(strGenomeGrp, "") != 0)
			{
				if(pBARPos->vSeqData[ni]->pSeqGroupName == NULL)
					continue;

				if(strcmp(strGenomeGrp, pBARPos->vSeqData[ni]->pSeqGroupName->m_pString) != 0)
					continue;
			}
		}

		for(nj=0; nj<pBARPos->vSeqData[ni]->nDataNum; nj++)
		{
			pBARPos->vSeqData[ni]->vData[3]->pMatElement[nj] /= pV->pMatElement[nh];
			nh++;
		}
	}

	if(nh != nDataTotNum)
	{
		printf("Error: TileProbe_MAT_Correction, probe number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* 4: Export */
	pCol = NULL;
	pCol = CreateIntMatrix(1, pBARPos->nColNum);
	if(pCol == NULL)
	{
		printf("Error: TileProbe_MAT_Correction, cannot write background corrected intensities to file!\n");
		exit(EXIT_FAILURE);
	}
	pCol->pMatElement[0] = 1;
	pCol->pMatElement[3] = 1;
	Affy_SaveFilteredBAR_Columns_Fast(strOutFile, pBARPos, pCol, strGenomeGrp);
	DestroyIntMatrix(pCol);

	/* release memory */
	DestroyDoubleMatrix(pE);
	DestroyDoubleMatrix(pV);


	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  TileProbe_BARBuildv2_Main()                                            */
/*  Build a probe background model for microarray BAR files.               */
/*  The variance of each probe is the mean within group variance           */
/* ----------------------------------------------------------------------- */ 
int TileProbe_BARBuildv2_Main(char strInputPath[], char strOutputPath[], 
						 int nTakeLog, double dNormLowerBound, int nShrink, int nTest)
{
	/* define */
	FILE *fpIn;

	/* arrays */
	int nArrayNum = 0;
	int nSeqNum = 0;
	int nProbeNum = 0;
	int nProbeNumC;

	struct tagBARData **vBARData = NULL;
	struct DOUBLEMATRIX *pProbeMedian;
	struct DOUBLEMATRIX *pProbeInterR;
	struct DOUBLEMATRIX *pProbeData;
	struct DOUBLEMATRIX *pProbeDataSort;
	struct LONGMATRIX *pProbeDataSortId;
	
	/* others */
	int ni,nj,nk,nu,nx,ny,nz;
	double dLog2 = log(2.0);
	int nQ50;
	char strOutFile[MED_LINE_LENGTH];
	char strFileName[MED_LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	struct INTMATRIX *pCol;
	int *vGroupID;

	/* for IQR shrinking */
	double dB;
	double dM,dV;
	int nGN,nDf;
	double dQmin,dQmax;

	/* get array number */
	nArrayNum = 0;
	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileProbe_BARBuild_Main, cannot open the input array list!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		nArrayNum++;
	}
	
	fclose(fpIn);
	
	printf("########################################\n");
	printf("# Loading Raw Data                     #\n");
	printf("########################################\n");
	printf(" Array number = %d\n", nArrayNum);

	if(nArrayNum == 0)
	{
		return PROC_SUCCESS;
	}

	/* memory for group id */
	vGroupID = NULL;
	vGroupID = (int *)calloc(nArrayNum, sizeof(int));
	if(vGroupID == NULL)
	{
		printf("Error: TileProbe_BARBuild_Main, cannot allocate memory for group id!\n");
		exit(EXIT_FAILURE);
	}

	/* load arrays */
	vBARData = NULL;
	vBARData = (struct tagBARData **)calloc(nArrayNum, sizeof(struct tagBARData *));
	if(vBARData == NULL)
	{
		printf("Error: TileProbe_BARBuild_Main, cannot allocate memory for the array list!\n");
		exit(EXIT_FAILURE);
	}

	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileProbe_BARBuild_Main, cannot open the input array list!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%d %s", &nu, strFileName);
		vGroupID[ni] = nu;

		printf("Loading %s ...\n", strFileName);
		vBARData[ni] = Affy_LoadBAR_Fast(strFileName);
		if(vBARData[ni] == NULL)
		{
			printf("Error: TileProbe_BARBuild_Main, cannot load *.bar file!\n");
			exit(EXIT_FAILURE);
		}

		nProbeNumC = 0;
		for(nj=0; nj<vBARData[ni]->nSeqNum; nj++)
		{
			nProbeNumC += vBARData[ni]->vSeqData[nj]->nDataNum;
		}
		
		if(ni == 0)
		{
			nSeqNum = vBARData[ni]->nSeqNum;
			nProbeNum = nProbeNumC;
			
		}
		else
		{
			if( (nSeqNum != vBARData[ni]->nSeqNum) || (nProbeNum != nProbeNumC) )
			{
				printf("Error: TileProbe_BARBuild_Main, file dimensions do not match!\n");
				exit(EXIT_FAILURE);
			}
		}

		if(nTakeLog == 1)
		{
			for(nj=0; nj<vBARData[ni]->nSeqNum; nj++)
			{
				for(nk=0; nk<vBARData[ni]->vSeqData[nj]->nDataNum; nk++)
				{
					/* truncate */
					if(vBARData[ni]->vSeqData[nj]->vData[1]->pMatElement[nk] < dNormLowerBound)
						vBARData[ni]->vSeqData[nj]->vData[1]->pMatElement[nk] = dNormLowerBound;

					/* log transformation */
					vBARData[ni]->vSeqData[nj]->vData[1]->pMatElement[nk] = log(vBARData[ni]->vSeqData[nj]->vData[1]->pMatElement[nk])/dLog2;
				}
			}
		}
		else
		{
			for(nj=0; nj<vBARData[ni]->nSeqNum; nj++)
			{
				for(nk=0; nk<vBARData[ni]->vSeqData[nj]->nDataNum; nk++)
				{
					/* truncate */
					if(vBARData[ni]->vSeqData[nj]->vData[1]->pMatElement[nk] < dNormLowerBound)
						vBARData[ni]->vSeqData[nj]->vData[1]->pMatElement[nk] = dNormLowerBound;
				}
			}
		}

		ni++;
	}
	
	fclose(fpIn);

	if(ni != nArrayNum)
	{
		printf("Error: TileProbe_BARBuild_Main, array number inconsistent!\n");
		exit(EXIT_FAILURE);
	}

	/* build probe model */
	printf("########################################\n");
	printf("# Build Probe Model                    #\n");
	printf("########################################\n");

	pProbeMedian = NULL;
	pProbeMedian = CreateDoubleMatrix(1,nProbeNum);
	if(pProbeMedian == NULL)
	{
		printf("Error: TileProbe_BARBuild_Main, cannot allocate memory for probe median!\n");
		exit(EXIT_FAILURE);
	}

	pProbeInterR = NULL;
	pProbeInterR = CreateDoubleMatrix(1,nProbeNum);
	if(pProbeInterR == NULL)
	{
		printf("Error: TileProbe_BARBuild_Main, cannot allocate memory for probe interquatile range!\n");
		exit(EXIT_FAILURE);
	}

	pProbeData = NULL;
	pProbeData = CreateDoubleMatrix(1,nArrayNum);
	if(pProbeData == NULL)
	{
		printf("Error: TileProbe_BARBuild_Main, cannot allocate memory for modeling probe effect!\n");
		exit(EXIT_FAILURE);
	}

	nQ50 = (int)(nArrayNum*0.5);
	if( (nArrayNum*0.5)-nQ50 > 1e-6 )
		nQ50++;
	nQ50--;

	if(nQ50 < 0)
		nQ50 = 0;
	if(nQ50 >= nArrayNum)
		nQ50 = nArrayNum-1;

	nk = 0;
	nDf = 0;
	for(ni=0; ni<vBARData[0]->nSeqNum; ni++)
	{
		for(nj=0; nj<vBARData[0]->vSeqData[ni]->nDataNum; nj++)
		{
			if(nk%100000 == 0)
			{
				printf("%d ...\n", nk);
			}

			pProbeDataSort = NULL;
			pProbeDataSortId = NULL;

			for(nu=0; nu<nArrayNum; nu++)
			{
				pProbeData->pMatElement[nu] = vBARData[nu]->vSeqData[ni]->vData[1]->pMatElement[nj];
			}

			nx = 0;
			dV = 0.0;
			nDf = 0;
			while(nx < nArrayNum)
			{
				dM = 0.0;
				
				nGN = 0;
				for(ny=nx; ny<nArrayNum; ny++)
				{
					if(vGroupID[ny] != vGroupID[nx])
						break;
					dM += pProbeData->pMatElement[ny];
					nGN++;
				}
				dM = dM/nGN;
				for(nz=nx; nz<ny; nz++)
				{
					dV += (pProbeData->pMatElement[nz]-dM)*(pProbeData->pMatElement[nz]-dM);
				}
				nDf = nDf+nGN-1;
				nx = ny;
			}
			dV /= nDf;
			pProbeInterR->pMatElement[nk] = dV;

			DMSORTMERGEA_0(pProbeData, &pProbeDataSort, &pProbeDataSortId);
			pProbeMedian->pMatElement[nk] = pProbeDataSort->pMatElement[nQ50];
			

			DestroyDoubleMatrix(pProbeDataSort);
			DestroyLongMatrix(pProbeDataSortId);

			nk++;
		}
	}

	if(nk != nProbeNum)
	{
		printf("Error: TileProbe_BARBuild_Main, probe number do not match!\n");
		exit(EXIT_FAILURE);
	}
	if(nDf == 0)
	{
		printf("Error: TileProbe_BARBuild_Main, no degree of freedom to compute variance!\n");
		exit(EXIT_FAILURE);
	}

	/* Variance shrinkage */
	if(nShrink != 0)
	{
		dM = 0.0;
		dV = 0.0;
		dQmin = 1e10;
		dQmax = -1e10;

		for(nj=0; nj<nProbeNum; nj++)
		{
			if(pProbeInterR->pMatElement[nj] < dQmin)
				dQmin = pProbeInterR->pMatElement[nj];
			if(pProbeInterR->pMatElement[nj] > dQmax)
				dQmax = pProbeInterR->pMatElement[nj];

			dM += pProbeInterR->pMatElement[nj];
			dV += pProbeInterR->pMatElement[nj]*pProbeInterR->pMatElement[nj];
		}

		dM /= nProbeNum;
		dV = dV-nProbeNum*dM*dM;
		dB = 2.0*(nProbeNum-1)/(2.0+nDf)/nProbeNum+2.0*dM*dM*(nProbeNum-1)/(2.0+nDf)/dV;
		
		if(dB>1.0)
			dB=1.0;
		else if(dB<0.0)
			dB=0.0;

		if(nShrink == -1)
			dB = 0.0;

		for(nj=0; nj<nProbeNum; nj++)
		{
			pProbeInterR->pMatElement[nj] = sqrt(pProbeInterR->pMatElement[nj]*(1.0-dB)+dM*dB);
		}
		
		printf("min Var = %f; max Var = %f\n", dQmin, dQmax);
		printf("B = %f\n", dB);
	}

	/* save file */
	printf("########################################\n");
	printf("# Export Probe Model                   #\n");
	printf("########################################\n");

	pCol = NULL;
	pCol = CreateIntMatrix(1, vBARData[0]->nColNum);
	if(pCol == NULL)
	{
		printf("Error: TileProbe_BARBuild_Main, cannot create column indicator for exporting fitted model!\n");
		exit(EXIT_FAILURE);
	}
	pCol->pMatElement[0] = 1;
	pCol->pMatElement[1] = 1;

	sprintf(strOutFile, "%s.prm", strOutputPath);
	nk = 0;
	for(ni=0; ni<vBARData[0]->nSeqNum; ni++)
	{
		for(nj=0; nj<vBARData[0]->vSeqData[ni]->nDataNum; nj++)
		{
			vBARData[0]->vSeqData[ni]->vData[1]->pMatElement[nj] = pProbeMedian->pMatElement[nk];
			nk++;
		}
	}
	if(nk != nProbeNum)
	{
		printf("Error: TileProbe_BARBuild_Main, probe number do not match!\n");
		exit(EXIT_FAILURE);
	}
	Affy_SaveBAR_Columns_Fast(strOutFile, vBARData[0], pCol);

	sprintf(strOutFile, "%s.prv", strOutputPath);
	nk = 0;
	for(ni=0; ni<vBARData[0]->nSeqNum; ni++)
	{
		for(nj=0; nj<vBARData[0]->vSeqData[ni]->nDataNum; nj++)
		{
			vBARData[0]->vSeqData[ni]->vData[1]->pMatElement[nj] = pProbeInterR->pMatElement[nk];
			nk++;
		}
	}
	if(nk != nProbeNum)
	{
		printf("Error: TileProbe_BARBuild_Main, probe number do not match!\n");
		exit(EXIT_FAILURE);
	}
	Affy_SaveBAR_Columns_Fast(strOutFile, vBARData[0], pCol);

	DestroyIntMatrix(pCol);

	/* release memory */
	free(vGroupID);

	for(ni=0; ni<nArrayNum; ni++)
	{
		Affy_BARData_Destroy(vBARData+ni);
	}
	free(vBARData);

	DestroyDoubleMatrix(pProbeMedian);
	DestroyDoubleMatrix(pProbeInterR);
	DestroyDoubleMatrix(pProbeData);

	printf("Done!\n");

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_Buildv2_Main()                                               */
/*  Build a probe background model for microarray CEL files.               */
/* ----------------------------------------------------------------------- */ 
int TileProbe_Buildv2_Main(char strInputPath[], char strOutputPath[], 
						 int nTakeLog, double dNormLowerBound, int nShrink, 
						 int nLogAfterNorm, int nTest)
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	FILE *fpOut2;

	/* arrays */
	int nArrayNum = 0;
	int *vGroupID;
	struct tagCELData **vCELData = NULL;
	char strFileName[MED_LINE_LENGTH];
	int ni,nj,nu,nx,ny,nz;
	int nCELNumberCells = 0;
	int nCELTotalX = 0;
	int nCELTotalY = 0;

	struct DOUBLEMATRIX *pSortMean;
	struct DOUBLEMATRIX *pArray;
	struct LONGMATRIX *pSortId;

	struct DOUBLEMATRIX *pProbeMedian;
	struct DOUBLEMATRIX *pProbeInterR;
	struct DOUBLEMATRIX *pProbeData;
	struct DOUBLEMATRIX *pProbeDataSort;
	struct LONGMATRIX *pProbeDataSortId;

	double dM,dV,dB;
	int	nDf,nGN;

	struct DOUBLEMATRIX *pOldIntensity;
	struct DOUBLEMATRIX *pOldSD;
	
	/* others */
	double dLog2 = log(2.0);
	int nId;
	int nQ50;
	char strOutFile[MED_LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	double dQmin,dQmax;

	/* get array number */
	nArrayNum = 0;
	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileProbe_Build_Main, cannot open the input array list!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		nArrayNum++;
	}
	
	fclose(fpIn);
	
	printf("########################################\n");
	printf("# Loading Raw Data                     #\n");
	printf("########################################\n");
	printf(" Array number = %d\n", nArrayNum);

	if(nArrayNum == 0)
	{
		return PROC_SUCCESS;
	}

	/* memory for group id */
	vGroupID = NULL;
	vGroupID = (int *)calloc(nArrayNum, sizeof(int));
	if(vGroupID == NULL)
	{
		printf("Error: TileProbe_Build_Main, cannot allocate memory for group id!\n");
		exit(EXIT_FAILURE);
	}

	/* load arrays */
	vCELData = NULL;
	vCELData = (struct tagCELData **)calloc(nArrayNum, sizeof(struct tagCELData *));
	if(vCELData == NULL)
	{
		printf("Error: TileProbe_Build_Main, cannot allocate memory for the array list!\n");
		exit(EXIT_FAILURE);
	}

	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileProbe_Build_Main, cannot open the input array list!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%d %s", &nu, strFileName);
		vGroupID[ni] = nu;
		vCELData[ni] = TileMapv2_LoadCEL(strFileName);
		if(vCELData[ni] == NULL)
		{
			printf("Error: TileProbe_Build_Main, cannot load *.CEL file!\n");
			exit(EXIT_FAILURE);
		}

		if(ni == 0)
		{
			nCELNumberCells = vCELData[ni]->nNumberCells;
			nCELTotalX = vCELData[ni]->nCols;
			nCELTotalY = vCELData[ni]->nRows;
			if( nCELNumberCells != nCELTotalX*nCELTotalY)
			{
				printf("Error: TileProbe_Build_Main, CEL file Cols*Rows != NumberCells!\n");
				exit(EXIT_FAILURE);
			}
		}
		else
		{
			if( (nCELNumberCells != vCELData[ni]->nNumberCells) || (nCELTotalX != vCELData[ni]->nCols)
				|| (nCELTotalY != vCELData[ni]->nRows) )
			{
				printf("Error: TileProbe_Build_Main, CEL file dimensions do not match!\n");
				exit(EXIT_FAILURE);
			}
		}

		if(nLogAfterNorm == 0)
		{
			if(nTakeLog == 1)
			{
				for(nj=0; nj<vCELData[ni]->nNumberCells; nj++)
				{
					/* truncate */
					if(vCELData[ni]->pIntensity->pMatElement[nj] < dNormLowerBound)
						vCELData[ni]->pIntensity->pMatElement[nj] = dNormLowerBound;

					/* log transformation */
					vCELData[ni]->pIntensity->pMatElement[nj] = log(vCELData[ni]->pIntensity->pMatElement[nj])/dLog2;
				}
			}
			else
			{
				for(nj=0; nj<vCELData[ni]->nNumberCells; nj++)
				{
					/* truncate */
					if(vCELData[ni]->pIntensity->pMatElement[nj] < dNormLowerBound)
						vCELData[ni]->pIntensity->pMatElement[nj] = dNormLowerBound;
				}
			}
		}

		ni++;
	}
	
	fclose(fpIn);

	if(ni != nArrayNum)
	{
		printf("Error: TileProbe_Build_Main, array number inconsistent!\n");
		exit(EXIT_FAILURE);
	}

	printf("########################################\n");
	printf("# Quantile Normalization               #\n");
	printf("########################################\n");
	
	pSortMean = NULL;
	pSortMean = CreateDoubleMatrix(1, nCELNumberCells);
	if(pSortMean == NULL)
	{
		printf("Error: TileProbe_Build_Main, cannot allocate memory for quantile normalization!\n");
		exit(EXIT_FAILURE);
	}
				
	for(ni=0; ni<nArrayNum; ni++)
	{
		pArray = NULL;
		pSortId = NULL;

		/* export sorted values */
		printf("  Sorting array %d...\n", ni);
		DMSORTMERGEA_0(vCELData[ni]->pIntensity, &pArray, &pSortId);

		/* compute percentiles */
		for(nj=0; nj<nCELNumberCells; nj++)
		{
			pSortMean->pMatElement[nj] += pArray->pMatElement[nj];
			nId = pSortId->pMatElement[nj];
			vCELData[ni]->pIntensity->pMatElement[nId] = nj;
		}

		DestroyDoubleMatrix(pArray);
		DestroyLongMatrix(pSortId);
	}

	for(nj=0; nj<nCELNumberCells; nj++)
	{
		pSortMean->pMatElement[nj] /= nArrayNum;
	}

	for(ni=0; ni<nArrayNum; ni++)
	{
		for(nj=0; nj<nCELNumberCells; nj++)
		{
			nId = (long)(vCELData[ni]->pIntensity->pMatElement[nj]);
			vCELData[ni]->pIntensity->pMatElement[nj] = pSortMean->pMatElement[nId];
		}

		if(nLogAfterNorm != 0)
		{
			if(nTakeLog == 1)
			{
				for(nj=0; nj<vCELData[ni]->nNumberCells; nj++)
				{
					/* truncate */
					if(vCELData[ni]->pIntensity->pMatElement[nj] < dNormLowerBound)
						vCELData[ni]->pIntensity->pMatElement[nj] = dNormLowerBound;

					/* log transformation */
					vCELData[ni]->pIntensity->pMatElement[nj] = log(vCELData[ni]->pIntensity->pMatElement[nj])/dLog2;
				}
			}
			else
			{
				for(nj=0; nj<vCELData[ni]->nNumberCells; nj++)
				{
					/* truncate */
					if(vCELData[ni]->pIntensity->pMatElement[nj] < dNormLowerBound)
						vCELData[ni]->pIntensity->pMatElement[nj] = dNormLowerBound;
				}
			}
		}
	}

	/* build probe model */
	printf("########################################\n");
	printf("# Build Probe Model                    #\n");
	printf("########################################\n");

	if(nTest == 1)
	{
		sprintf(strOutFile, "%s.testq", strOutputPath);
		fpOut = NULL;
		fpOut = fopen(strOutFile, "w");
		if(fpOut == NULL)
		{
			printf("Error: TileProbe_Build_Main, cannot open the output argument file!\n");
			exit(EXIT_FAILURE);
		}
		sprintf(strOutFile, "%s.testprobe", strOutputPath);
		fpOut2 = NULL;
		fpOut2 = fopen(strOutFile, "w");
		if(fpOut2 == NULL)
		{
			printf("Error: TileProbe_Build_Main, cannot open the output argument file!\n");
			exit(EXIT_FAILURE);
		}
	}

	pProbeMedian = NULL;
	pProbeMedian = CreateDoubleMatrix(1,nCELNumberCells);
	if(pProbeMedian == NULL)
	{
		printf("Error: TileProbe_Build_Main, cannot allocate memory for probe median!\n");
		exit(EXIT_FAILURE);
	}

	pProbeInterR = NULL;
	pProbeInterR = CreateDoubleMatrix(1,nCELNumberCells);
	if(pProbeInterR == NULL)
	{
		printf("Error: TileProbe_Build_Main, cannot allocate memory for probe interquatile range!\n");
		exit(EXIT_FAILURE);
	}

	pProbeData = NULL;
	pProbeData = CreateDoubleMatrix(1,nArrayNum);
	if(pProbeData == NULL)
	{
		printf("Error: TileProbe_Build_Main, cannot allocate memory for modeling probe effect!\n");
		exit(EXIT_FAILURE);
	}

	nQ50 = (int)(nArrayNum*0.5);
	if( (nArrayNum*0.5)-nQ50 > 1e-6 )
		nQ50++;
	nQ50--;

	if(nQ50 < 0)
		nQ50 = 0;
	if(nQ50 >= nArrayNum)
		nQ50 = nArrayNum-1;

	nDf = 0;
	for(nj=0; nj<nCELNumberCells; nj++)
	{
		if(nj%100000 == 0)
		{
			printf("%d ...\n", nj);
		}

		pProbeDataSort = NULL;
		pProbeDataSortId = NULL;

		for(ni=0; ni<nArrayNum; ni++)
		{
			pProbeData->pMatElement[ni] = vCELData[ni]->pIntensity->pMatElement[nj];
		}

		nx = 0;
		dV = 0.0;
		nDf = 0;
		while(nx < nArrayNum)
		{
			dM = 0.0;
			
			nGN = 0;
			for(ny=nx; ny<nArrayNum; ny++)
			{
				if(vGroupID[ny] != vGroupID[nx])
					break;
				dM += pProbeData->pMatElement[ny];
				nGN++;
			}
			dM = dM/nGN;
			for(nz=nx; nz<ny; nz++)
			{
				dV += (pProbeData->pMatElement[nz]-dM)*(pProbeData->pMatElement[nz]-dM);
			}
			nDf = nDf+nGN-1;
			nx = ny;
		}
		dV /= nDf;
		pProbeInterR->pMatElement[nj] = dV;

		DMSORTMERGEA_0(pProbeData, &pProbeDataSort, &pProbeDataSortId);
		pProbeMedian->pMatElement[nj] = pProbeDataSort->pMatElement[nQ50];

		if(nTest == 1)
		{
			fprintf(fpOut, "%e\n", pProbeInterR->pMatElement[nj]); 
			if(rand_u() < 0.01)
			{
				fprintf(fpOut2, "%f", pProbeData->pMatElement[0]);
				for(ni=1; ni<nArrayNum; ni++)
				{
					fprintf(fpOut2, "\t%f", pProbeData->pMatElement[ni]);
				}
				fprintf(fpOut2, "\n");
			}
		}

		DestroyDoubleMatrix(pProbeDataSort);
		DestroyLongMatrix(pProbeDataSortId);
	}

	/* IQR shrinkage */
		if(nDf == 0)
	{
		printf("Error: TileProbe_BARBuild_Main, no degree of freedom to compute variance!\n");
		exit(EXIT_FAILURE);
	}

	/* Variance shrinkage */
	if(nShrink != 0)
	{
		dM = 0.0;
		dV = 0.0;
		dQmin = 1e10;
		dQmax = -1e10;

		for(nj=0; nj<nCELNumberCells; nj++)
		{
			if(pProbeInterR->pMatElement[nj] < dQmin)
				dQmin = pProbeInterR->pMatElement[nj];
			if(pProbeInterR->pMatElement[nj] > dQmax)
				dQmax = pProbeInterR->pMatElement[nj];

			dM += pProbeInterR->pMatElement[nj];
			dV += pProbeInterR->pMatElement[nj]*pProbeInterR->pMatElement[nj];
		}

		dM /= nCELNumberCells;
		dV = dV-nCELNumberCells*dM*dM;
		dB = 2.0*(nCELNumberCells-1)/(2.0+nDf)/nCELNumberCells+2.0*dM*dM*(nCELNumberCells-1)/(2.0+nDf)/dV;
		
		if(dB>1.0)
			dB=1.0;
		else if(dB<0.0)
			dB=0.0;

		if(nShrink == -1)
			dB = 0.0;

		for(nj=0; nj<nCELNumberCells; nj++)
		{
			pProbeInterR->pMatElement[nj] = sqrt(pProbeInterR->pMatElement[nj]*(1.0-dB)+dM*dB);
		}
		
		printf("min Var = %f; max Var = %f\n", dQmin, dQmax);
		printf("B = %f\n", dB);
	}

	if(nTest == 1)
	{
		fclose(fpOut);
		fclose(fpOut2);
	}

	/* save file */
	printf("########################################\n");
	printf("# Export Probe Model                   #\n");
	printf("########################################\n");
	pOldIntensity = NULL;
	pOldSD = NULL;
	pOldIntensity = vCELData[0]->pIntensity;
	pOldSD = vCELData[0]->pSD;
	vCELData[0]->pIntensity = pProbeMedian;
	vCELData[0]->pSD = pProbeInterR;

	sprintf(strOutFile, "%s.prbm", strOutputPath);
	Affy_SaveCELv4(strOutFile, vCELData[0]);

	vCELData[0]->pIntensity = pSortMean;
	sprintf(strOutFile, "%s.quan", strOutputPath);
	Affy_SaveCELv4(strOutFile, vCELData[0]);

	vCELData[0]->pIntensity = pOldIntensity;
	vCELData[0]->pSD = pOldSD;


	/* release memory */
	free(vGroupID);
	for(ni=0; ni<nArrayNum; ni++)
	{
		Affy_CELData_Destroy(vCELData+ni);
	}
	free(vCELData);

	DestroyDoubleMatrix(pSortMean);
	DestroyDoubleMatrix(pProbeMedian);
	DestroyDoubleMatrix(pProbeInterR);
	DestroyDoubleMatrix(pProbeData);

	printf("Done!\n");

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_BuildHMMT_Main()                                             */
/*  Build a probe background model for microarray CEL files using HMM      */
/*  tiling algorithm.                                                      */
/* ----------------------------------------------------------------------- */ 
int TileProbe_BuildHMMT_Main(char strInputPath[], char strOutputPath[], 
						 int nTakeLog, double dNormLowerBound, 
						 double dFilterTopPrc, int nTest)
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	FILE *fpOut2;
	double dNAValue = -1e10;

	/* arrays */
	int nArrayNum = 0;
	int *vGroupID;
	struct tagCELData **vCELData = NULL;
	char strFileName[MED_LINE_LENGTH];
	int ni,nj,nu;
	int nCELNumberCells = 0;
	int nCELTotalX = 0;
	int nCELTotalY = 0;
	int nFilterCutId = 0;

	struct DOUBLEMATRIX *pSortMean;
	struct DOUBLEMATRIX *pArray;
	struct LONGMATRIX *pSortId;

	struct DOUBLEMATRIX *pProbeMedian;
	struct DOUBLEMATRIX *pProbeInterR;
	struct DOUBLEMATRIX *pProbeData;
	struct DOUBLEMATRIX *pProbeDataSort;
	struct LONGMATRIX *pProbeDataSortId;

	double dM,dV;
	int	nDf;

	struct DOUBLEMATRIX *pOldIntensity;
	struct DOUBLEMATRIX *pOldSD;
	
	/* others */
	double dLog2 = log(2.0);
	int nId;
	int nQ50;
	char strOutFile[MED_LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	double dQmin,dQmax;

	/* get array number */
	nArrayNum = 0;
	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileProbe_BuildHMMT_Main, cannot open the input array list!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		nArrayNum++;
	}
	
	fclose(fpIn);
	
	printf("########################################\n");
	printf("# Loading Raw Data                     #\n");
	printf("########################################\n");
	printf(" Array number = %d\n", nArrayNum);

	if(nArrayNum == 0)
	{
		return PROC_SUCCESS;
	}

	/* memory for group id */
	vGroupID = NULL;
	vGroupID = (int *)calloc(nArrayNum, sizeof(int));
	if(vGroupID == NULL)
	{
		printf("Error: TileProbe_BuildHMMT_Main, cannot allocate memory for group id!\n");
		exit(EXIT_FAILURE);
	}

	/* load arrays */
	vCELData = NULL;
	vCELData = (struct tagCELData **)calloc(nArrayNum, sizeof(struct tagCELData *));
	if(vCELData == NULL)
	{
		printf("Error: TileProbe_BuildHMMT_Main, cannot allocate memory for the array list!\n");
		exit(EXIT_FAILURE);
	}

	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileProbe_BuildHMMT_Main, cannot open the input array list!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%d %s", &nu, strFileName);
		vGroupID[ni] = nu;
		vCELData[ni] = TileMapv2_LoadCEL(strFileName);
		if(vCELData[ni] == NULL)
		{
			printf("Error: TileProbe_BuildHMMT_Main, cannot load *.CEL file!\n");
			exit(EXIT_FAILURE);
		}

		if(ni == 0)
		{
			nCELNumberCells = vCELData[ni]->nNumberCells;
			nCELTotalX = vCELData[ni]->nCols;
			nCELTotalY = vCELData[ni]->nRows;
			if( nCELNumberCells != nCELTotalX*nCELTotalY)
			{
				printf("Error: TileProbe_BuildHMMT_Main, CEL file Cols*Rows != NumberCells!\n");
				exit(EXIT_FAILURE);
			}
		}
		else
		{
			if( (nCELNumberCells != vCELData[ni]->nNumberCells) || (nCELTotalX != vCELData[ni]->nCols)
				|| (nCELTotalY != vCELData[ni]->nRows) )
			{
				printf("Error: TileProbe_BuildHMMT_Main, CEL file dimensions do not match!\n");
				exit(EXIT_FAILURE);
			}
		}

		if(nTakeLog == 1)
		{
			for(nj=0; nj<vCELData[ni]->nNumberCells; nj++)
			{
				/* truncate */
				if(vCELData[ni]->pIntensity->pMatElement[nj] < dNormLowerBound)
					vCELData[ni]->pIntensity->pMatElement[nj] = dNormLowerBound;

				/* log transformation */
				vCELData[ni]->pIntensity->pMatElement[nj] = log(vCELData[ni]->pIntensity->pMatElement[nj])/dLog2;
			}
		}
		else
		{
			for(nj=0; nj<vCELData[ni]->nNumberCells; nj++)
			{
				/* truncate */
				if(vCELData[ni]->pIntensity->pMatElement[nj] < dNormLowerBound)
					vCELData[ni]->pIntensity->pMatElement[nj] = dNormLowerBound;
			}
		}

		ni++;
	}
	
	fclose(fpIn);

	if(ni != nArrayNum)
	{
		printf("Error: TileProbe_BuildHMMT_Main, array number inconsistent!\n");
		exit(EXIT_FAILURE);
	}

	printf("########################################\n");
	printf("# Quantile Normalization               #\n");
	printf("########################################\n");
	
	pSortMean = NULL;
	pSortMean = CreateDoubleMatrix(1, nCELNumberCells);
	if(pSortMean == NULL)
	{
		printf("Error: TileProbe_BuildHMMT_Main, cannot allocate memory for quantile normalization!\n");
		exit(EXIT_FAILURE);
	}
				
	for(ni=0; ni<nArrayNum; ni++)
	{
		pArray = NULL;
		pSortId = NULL;

		/* export sorted values */
		printf("  Sorting array %d...\n", ni);
		DMSORTMERGEA_0(vCELData[ni]->pIntensity, &pArray, &pSortId);

		/* compute percentiles */
		for(nj=0; nj<nCELNumberCells; nj++)
		{
			pSortMean->pMatElement[nj] += pArray->pMatElement[nj];
			nId = pSortId->pMatElement[nj];
			vCELData[ni]->pIntensity->pMatElement[nId] = nj;
		}

		DestroyDoubleMatrix(pArray);
		DestroyLongMatrix(pSortId);
	}

	for(nj=0; nj<nCELNumberCells; nj++)
	{
		pSortMean->pMatElement[nj] /= nArrayNum;
	}

	nFilterCutId = nCELNumberCells-(int)(nCELNumberCells*dFilterTopPrc);
	if(nFilterCutId < 0)
		nFilterCutId = 0;

	for(ni=0; ni<nArrayNum; ni++)
	{
		if(vGroupID[ni] > 0)
		{
			for(nj=0; nj<nCELNumberCells; nj++)
			{
				nId = (long)(vCELData[ni]->pIntensity->pMatElement[nj]);
				if(nId < nFilterCutId)
					vCELData[ni]->pIntensity->pMatElement[nj] = pSortMean->pMatElement[nId];
				else
					vCELData[ni]->pIntensity->pMatElement[nj] = dNAValue;

			}
		}
		else
		{
			for(nj=0; nj<nCELNumberCells; nj++)
			{
				nId = (long)(vCELData[ni]->pIntensity->pMatElement[nj]);
				vCELData[ni]->pIntensity->pMatElement[nj] = pSortMean->pMatElement[nId];
			}
		}
	}
	dNAValue = dNAValue+1.0;

	/* build probe model */
	printf("########################################\n");
	printf("# Build Probe Model                    #\n");
	printf("########################################\n");

	if(nTest == 1)
	{
		sprintf(strOutFile, "%s.testq", strOutputPath);
		fpOut = NULL;
		fpOut = fopen(strOutFile, "w");
		if(fpOut == NULL)
		{
			printf("Error: TileProbe_BuildHMMT_Main, cannot open the output argument file!\n");
			exit(EXIT_FAILURE);
		}
		sprintf(strOutFile, "%s.testprobe", strOutputPath);
		fpOut2 = NULL;
		fpOut2 = fopen(strOutFile, "w");
		if(fpOut2 == NULL)
		{
			printf("Error: TileProbe_BuildHMMT_Main, cannot open the output argument file!\n");
			exit(EXIT_FAILURE);
		}
	}

	pProbeMedian = NULL;
	pProbeMedian = CreateDoubleMatrix(1,nCELNumberCells);
	if(pProbeMedian == NULL)
	{
		printf("Error: TileProbe_BuildHMMT_Main, cannot allocate memory for probe median!\n");
		exit(EXIT_FAILURE);
	}

	pProbeInterR = NULL;
	pProbeInterR = CreateDoubleMatrix(1,nCELNumberCells);
	if(pProbeInterR == NULL)
	{
		printf("Error: TileProbe_BuildHMMT_Main, cannot allocate memory for probe interquatile range!\n");
		exit(EXIT_FAILURE);
	}

	pProbeData = NULL;
	pProbeData = CreateDoubleMatrix(1,nArrayNum);
	if(pProbeData == NULL)
	{
		printf("Error: TileProbe_BuildHMMT_Main, cannot allocate memory for modeling probe effect!\n");
		exit(EXIT_FAILURE);
	}

	nQ50 = (int)(nArrayNum*0.5);
	if( (nArrayNum*0.5)-nQ50 > 1e-6 )
		nQ50++;
	nQ50--;

	if(nQ50 < 0)
		nQ50 = 0;
	if(nQ50 >= nArrayNum)
		nQ50 = nArrayNum-1;

	dQmin = 1e10;
	dQmax = -1e10;
	for(nj=0; nj<nCELNumberCells; nj++)
	{
		if(nj%100000 == 0)
		{
			printf("%d ...\n", nj);
		}

		pProbeDataSort = NULL;
		pProbeDataSortId = NULL;

		dM = 0.0;
		dV = 0.0;
		nDf = 0;
		for(ni=0; ni<nArrayNum; ni++)
		{
			pProbeData->pMatElement[ni] = vCELData[ni]->pIntensity->pMatElement[nj];
			if(pProbeData->pMatElement[ni] > dNAValue)
			{
				nDf += 1;
				dM += pProbeData->pMatElement[ni];
			}
		}

		if(nDf <= 1)
		{
			pProbeMedian->pMatElement[nj] = 0.0;
			pProbeInterR->pMatElement[nj] = 1e6;
			continue;
		}

		dM = dM/nDf;
		pProbeMedian->pMatElement[nj] = dM;
		nDf = 0;
		for(ni=0; ni<nArrayNum; ni++)
		{
			if(pProbeData->pMatElement[ni] > dNAValue)
			{
				dV += (pProbeData->pMatElement[ni]-dM)*(pProbeData->pMatElement[ni]-dM);
				nDf += 1;
			}
		}
		nDf -= 1;
		dV /= nDf;
		pProbeInterR->pMatElement[nj] = sqrt(dV);
		if(pProbeInterR->pMatElement[nj] < dQmin)
			dQmin = pProbeInterR->pMatElement[nj];
		if(pProbeInterR->pMatElement[nj] > dQmax)
			dQmax = pProbeInterR->pMatElement[nj];

		if(nTest == 1)
		{
			fprintf(fpOut, "%e\n", pProbeInterR->pMatElement[nj]); 
			if(rand_u() < 0.01)
			{
				fprintf(fpOut2, "%f", pProbeData->pMatElement[0]);
				for(ni=1; ni<nArrayNum; ni++)
				{
					fprintf(fpOut2, "\t%f", pProbeData->pMatElement[ni]);
				}
				fprintf(fpOut2, "\n");
			}
		}
	}

	if(nTest == 1)
	{
		fclose(fpOut);
		fclose(fpOut2);
	}

	printf("min Var = %f; max Var = %f\n", dQmin, dQmax);

	/* save file */
	printf("########################################\n");
	printf("# Export Probe Model                   #\n");
	printf("########################################\n");
	pOldIntensity = NULL;
	pOldSD = NULL;
	pOldIntensity = vCELData[0]->pIntensity;
	pOldSD = vCELData[0]->pSD;
	vCELData[0]->pIntensity = pProbeMedian;
	vCELData[0]->pSD = pProbeInterR;

	sprintf(strOutFile, "%s.prbm", strOutputPath);
	Affy_SaveCELv4(strOutFile, vCELData[0]);

	vCELData[0]->pIntensity = pSortMean;
	sprintf(strOutFile, "%s.quan", strOutputPath);
	Affy_SaveCELv4(strOutFile, vCELData[0]);

	vCELData[0]->pIntensity = pOldIntensity;
	vCELData[0]->pSD = pOldSD;


	/* release memory */
	free(vGroupID);
	for(ni=0; ni<nArrayNum; ni++)
	{
		Affy_CELData_Destroy(vCELData+ni);
	}
	free(vCELData);

	DestroyDoubleMatrix(pSortMean);
	DestroyDoubleMatrix(pProbeMedian);
	DestroyDoubleMatrix(pProbeInterR);
	DestroyDoubleMatrix(pProbeData);

	printf("Done!\n");

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  TileProbe_BuildHMMB_Main()                                             */
/*  Build a probe background model for microarray CEL files using HMM      */
/*  tiling algorithm with variance shrinking.                              */
/* ----------------------------------------------------------------------- */ 
int TileProbe_BuildHMMB_Main(char strInputPath[], char strOutputPath[], 
						 int nTakeLog, double dNormLowerBound, 
						 double dFilterTopPrc, int nTest)
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	FILE *fpOut2;
	double dNAValue = -1e10;

	/* arrays */
	int nArrayNum = 0;
	int *vGroupID;
	struct tagCELData **vCELData = NULL;
	char strFileName[MED_LINE_LENGTH];
	int ni,nj,nu;
	int nCELNumberCells = 0;
	int nCELTotalX = 0;
	int nCELTotalY = 0;
	int nFilterCutId = 0;

	struct DOUBLEMATRIX *pSortMean;
	struct DOUBLEMATRIX *pArray;
	struct LONGMATRIX *pSortId;

	struct DOUBLEMATRIX *pProbeMedian;
	struct DOUBLEMATRIX *pProbeInterR;
	struct INTMATRIX *pDf;
	struct DOUBLEMATRIX *pProbeData;
	struct DOUBLEMATRIX *pProbeDataSort;
	struct LONGMATRIX *pProbeDataSortId;

	struct INTMATRIX *pN;
	struct DOUBLEMATRIX *pS;
	struct DOUBLEMATRIX *pX;
	struct DOUBLEMATRIX *pV;

	double dM,dV,dV0,dGM,dB;
	int	nDf,nMaxDf,nD0,nGN;

	struct DOUBLEMATRIX *pOldIntensity;
	struct DOUBLEMATRIX *pOldSD;
	
	/* others */
	double dLog2 = log(2.0);
	int nId;
	int nQ50;
	char strOutFile[MED_LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	double dQmin,dQmax;

	/* get array number */
	nArrayNum = 0;
	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileProbe_BuildHMMB_Main, cannot open the input array list!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		nArrayNum++;
	}
	
	fclose(fpIn);
	
	printf("########################################\n");
	printf("# Loading Raw Data                     #\n");
	printf("########################################\n");
	printf(" Array number = %d\n", nArrayNum);

	if(nArrayNum == 0)
	{
		return PROC_SUCCESS;
	}

	/* memory for group id */
	vGroupID = NULL;
	vGroupID = (int *)calloc(nArrayNum, sizeof(int));
	if(vGroupID == NULL)
	{
		printf("Error: TileProbe_BuildHMMB_Main, cannot allocate memory for group id!\n");
		exit(EXIT_FAILURE);
	}

	/* load arrays */
	vCELData = NULL;
	vCELData = (struct tagCELData **)calloc(nArrayNum, sizeof(struct tagCELData *));
	if(vCELData == NULL)
	{
		printf("Error: TileProbe_BuildHMMB_Main, cannot allocate memory for the array list!\n");
		exit(EXIT_FAILURE);
	}

	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileProbe_BuildHMMB_Main, cannot open the input array list!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%d %s", &nu, strFileName);
		vGroupID[ni] = nu;
		vCELData[ni] = TileMapv2_LoadCEL(strFileName);
		if(vCELData[ni] == NULL)
		{
			printf("Error: TileProbe_BuildHMMB_Main, cannot load *.CEL file!\n");
			exit(EXIT_FAILURE);
		}

		if(ni == 0)
		{
			nCELNumberCells = vCELData[ni]->nNumberCells;
			nCELTotalX = vCELData[ni]->nCols;
			nCELTotalY = vCELData[ni]->nRows;
			if( nCELNumberCells != nCELTotalX*nCELTotalY)
			{
				printf("Error: TileProbe_BuildHMMB_Main, CEL file Cols*Rows != NumberCells!\n");
				exit(EXIT_FAILURE);
			}
		}
		else
		{
			if( (nCELNumberCells != vCELData[ni]->nNumberCells) || (nCELTotalX != vCELData[ni]->nCols)
				|| (nCELTotalY != vCELData[ni]->nRows) )
			{
				printf("Error: TileProbe_BuildHMMT_Main, CEL file dimensions do not match!\n");
				exit(EXIT_FAILURE);
			}
		}

		if(nTakeLog == 1)
		{
			for(nj=0; nj<vCELData[ni]->nNumberCells; nj++)
			{
				/* truncate */
				if(vCELData[ni]->pIntensity->pMatElement[nj] < dNormLowerBound)
					vCELData[ni]->pIntensity->pMatElement[nj] = dNormLowerBound;

				/* log transformation */
				vCELData[ni]->pIntensity->pMatElement[nj] = log(vCELData[ni]->pIntensity->pMatElement[nj])/dLog2;
			}
		}
		else
		{
			for(nj=0; nj<vCELData[ni]->nNumberCells; nj++)
			{
				/* truncate */
				if(vCELData[ni]->pIntensity->pMatElement[nj] < dNormLowerBound)
					vCELData[ni]->pIntensity->pMatElement[nj] = dNormLowerBound;
			}
		}

		ni++;
	}
	
	fclose(fpIn);

	if(ni != nArrayNum)
	{
		printf("Error: TileProbe_BuildHMMB_Main, array number inconsistent!\n");
		exit(EXIT_FAILURE);
	}

	printf("########################################\n");
	printf("# Quantile Normalization               #\n");
	printf("########################################\n");
	
	pSortMean = NULL;
	pSortMean = CreateDoubleMatrix(1, nCELNumberCells);
	if(pSortMean == NULL)
	{
		printf("Error: TileProbe_BuildHMMB_Main, cannot allocate memory for quantile normalization!\n");
		exit(EXIT_FAILURE);
	}
				
	for(ni=0; ni<nArrayNum; ni++)
	{
		pArray = NULL;
		pSortId = NULL;

		/* export sorted values */
		printf("  Sorting array %d...\n", ni);
		DMSORTMERGEA_0(vCELData[ni]->pIntensity, &pArray, &pSortId);

		/* compute percentiles */
		for(nj=0; nj<nCELNumberCells; nj++)
		{
			pSortMean->pMatElement[nj] += pArray->pMatElement[nj];
			nId = pSortId->pMatElement[nj];
			vCELData[ni]->pIntensity->pMatElement[nId] = nj;
		}

		DestroyDoubleMatrix(pArray);
		DestroyLongMatrix(pSortId);
	}

	for(nj=0; nj<nCELNumberCells; nj++)
	{
		pSortMean->pMatElement[nj] /= nArrayNum;
	}

	nFilterCutId = nCELNumberCells-(int)(nCELNumberCells*dFilterTopPrc);
	if(nFilterCutId < 0)
		nFilterCutId = 0;

	for(ni=0; ni<nArrayNum; ni++)
	{
		if(vGroupID[ni] > 0)
		{
			for(nj=0; nj<nCELNumberCells; nj++)
			{
				nId = (long)(vCELData[ni]->pIntensity->pMatElement[nj]);
				if(nId < nFilterCutId)
					vCELData[ni]->pIntensity->pMatElement[nj] = pSortMean->pMatElement[nId];
				else
					vCELData[ni]->pIntensity->pMatElement[nj] = dNAValue;

			}
		}
		else
		{
			for(nj=0; nj<nCELNumberCells; nj++)
			{
				nId = (long)(vCELData[ni]->pIntensity->pMatElement[nj]);
				vCELData[ni]->pIntensity->pMatElement[nj] = pSortMean->pMatElement[nId];
			}
		}
	}
	dNAValue = dNAValue+1.0;

	/* build probe model */
	printf("########################################\n");
	printf("# Build Probe Model                    #\n");
	printf("########################################\n");

	if(nTest == 1)
	{
		sprintf(strOutFile, "%s.testq", strOutputPath);
		fpOut = NULL;
		fpOut = fopen(strOutFile, "w");
		if(fpOut == NULL)
		{
			printf("Error: TileProbe_BuildHMMB_Main, cannot open the output argument file!\n");
			exit(EXIT_FAILURE);
		}
		sprintf(strOutFile, "%s.testprobe", strOutputPath);
		fpOut2 = NULL;
		fpOut2 = fopen(strOutFile, "w");
		if(fpOut2 == NULL)
		{
			printf("Error: TileProbe_BuildHMMB_Main, cannot open the output argument file!\n");
			exit(EXIT_FAILURE);
		}
	}

	pProbeMedian = NULL;
	pProbeMedian = CreateDoubleMatrix(1,nCELNumberCells);
	if(pProbeMedian == NULL)
	{
		printf("Error: TileProbe_BuildHMMB_Main, cannot allocate memory for probe median!\n");
		exit(EXIT_FAILURE);
	}

	pProbeInterR = NULL;
	pProbeInterR = CreateDoubleMatrix(1,nCELNumberCells);
	if(pProbeInterR == NULL)
	{
		printf("Error: TileProbe_BuildHMMB_Main, cannot allocate memory for probe interquatile range!\n");
		exit(EXIT_FAILURE);
	}

	pDf = NULL;
	pDf = CreateIntMatrix(1,nCELNumberCells);
	if(pDf == NULL)
	{
		printf("Error: TileProbe_BuildHMMB_Main, cannot allocate memory for storing degrees of freedom!\n");
		exit(EXIT_FAILURE);
	}

	pProbeData = NULL;
	pProbeData = CreateDoubleMatrix(1,nArrayNum);
	if(pProbeData == NULL)
	{
		printf("Error: TileProbe_BuildHMMB_Main, cannot allocate memory for modeling probe effect!\n");
		exit(EXIT_FAILURE);
	}

	nQ50 = (int)(nArrayNum*0.5);
	if( (nArrayNum*0.5)-nQ50 > 1e-6 )
		nQ50++;
	nQ50--;

	if(nQ50 < 0)
		nQ50 = 0;
	if(nQ50 >= nArrayNum)
		nQ50 = nArrayNum-1;

	dQmin = 1e10;
	dQmax = -1e10;
	nMaxDf = 0;
	for(nj=0; nj<nCELNumberCells; nj++)
	{
		if(nj%100000 == 0)
		{
			printf("%d ...\n", nj);
		}

		pProbeDataSort = NULL;
		pProbeDataSortId = NULL;

		dM = 0.0;
		dV = 0.0;
		nDf = 0;
		for(ni=0; ni<nArrayNum; ni++)
		{
			pProbeData->pMatElement[ni] = vCELData[ni]->pIntensity->pMatElement[nj];
			if(pProbeData->pMatElement[ni] > dNAValue)
			{
				nDf += 1;
				dM += pProbeData->pMatElement[ni];
			}
		}

		if(nDf <= 1)
		{
			pDf->pMatElement[nj] = 0;
			pProbeMedian->pMatElement[nj] = 0.0;
			pProbeInterR->pMatElement[nj] = 1e6;
			continue;
		}

		dM = dM/nDf;
		pProbeMedian->pMatElement[nj] = dM;
		nDf = 0;
		for(ni=0; ni<nArrayNum; ni++)
		{
			if(pProbeData->pMatElement[ni] > dNAValue)
			{
				dV += (pProbeData->pMatElement[ni]-dM)*(pProbeData->pMatElement[ni]-dM);
				nDf += 1;
			}
		}
		nDf -= 1;
		dV /= nDf;
		pProbeInterR->pMatElement[nj] = dV;
		if(pProbeInterR->pMatElement[nj] < dQmin)
			dQmin = pProbeInterR->pMatElement[nj];
		if(pProbeInterR->pMatElement[nj] > dQmax)
			dQmax = pProbeInterR->pMatElement[nj];
		
		pDf->pMatElement[nj] = nDf;
		if(nDf > nMaxDf)
			nMaxDf = nDf;
	}

	printf("min Var = %f; max Var = %f\n", dQmin, dQmax);

	if(nMaxDf <= 0)
	{
		printf("Error: TileProbe_BuildHMMB_Main, don't have enough degrees of freedom to estimate shrinkage factor!\n");
		exit(EXIT_FAILURE);
	}

	pN = NULL;
	pN = CreateIntMatrix(1, nMaxDf+1);

	pS = NULL;
	pS = CreateDoubleMatrix(1, nMaxDf+1);

	pX = NULL;
	pX = CreateDoubleMatrix(1, nMaxDf+1);

	pV = NULL;
	pV = CreateDoubleMatrix(1, nMaxDf+1);

	if( (pN == NULL) || (pS == NULL) || (pX == NULL) || (pV == NULL) )
	{
		printf("Error: TileProbe_BuildHMMB_Main, cannot allocate memory for computing shrinkage factor!\n");
		exit(EXIT_FAILURE);
	}

	dGM = 0.0;
	nGN = 0;
	for(nj=0; nj<nCELNumberCells; nj++)
	{
		nDf = pDf->pMatElement[nj];
		if(nDf <= 0)
			continue;

		pN->pMatElement[nDf] += 1;
		pX->pMatElement[nDf] += pProbeInterR->pMatElement[nj];
		pS->pMatElement[nDf] += pProbeInterR->pMatElement[nj]*pProbeInterR->pMatElement[nj];
		dGM += pProbeInterR->pMatElement[nj];
		nGN += 1;
	}
	dGM /= nGN;

	nD0 = 0;
	dV0 = 0.0;
	for(ni=1; ni<=nMaxDf; ni++)
	{
		if(pN->pMatElement[ni] < 10)
			continue;

		dM = pX->pMatElement[ni]/pN->pMatElement[ni];
		dV = pS->pMatElement[ni]-pN->pMatElement[ni]*dM*dM;
		dB = 2.0*(pN->pMatElement[ni]-1)/(2.0+ni)/pN->pMatElement[ni]+2.0*dM*dM*(pN->pMatElement[ni]-1)/(2.0+ni)/dV;
		if(dB>1.0)
			dB=1.0;
		else if(dB<0.0)
			dB=0.0;

		pV->pMatElement[ni] = (1.0/dB-1.0)/ni;
		printf("%d\t%f\n", pN->pMatElement[ni], pV->pMatElement[ni]); 

		dV0 += pV->pMatElement[ni]*pN->pMatElement[ni];
		nD0 += pN->pMatElement[ni];
	}

	if(nD0 <= 0)
	{
		printf("Error: TileProbe_BuildHMMB_Main, don't have enough degrees of freedom to estimate shrinkage factor!\n");
		exit(EXIT_FAILURE);
	}

	dV0 /= nD0;
	printf("Prior C = %f\n", dV0); 
	for(ni=0; ni<=nMaxDf; ni++)
	{
		pV->pMatElement[ni] = 1.0/(1.0+ni*dV0);
	}
		
	for(nj=0; nj<nCELNumberCells; nj++)
	{
		nDf = pDf->pMatElement[nj];
		dB = pV->pMatElement[nDf];

		pProbeInterR->pMatElement[nj] = sqrt(pProbeInterR->pMatElement[nj]*(1.0-dB)+dGM*dB);

		if(nTest == 1)
		{
			fprintf(fpOut, "%e\n", pProbeInterR->pMatElement[nj]); 
			if(rand_u() < 0.01)
			{
				fprintf(fpOut2, "%f", pProbeData->pMatElement[0]);
				for(ni=1; ni<nArrayNum; ni++)
				{
					fprintf(fpOut2, "\t%f", pProbeData->pMatElement[ni]);
				}
				fprintf(fpOut2, "\n");
			}
		}
	}

	if(nTest == 1)
	{
		fclose(fpOut);
		fclose(fpOut2);
	}

	DestroyIntMatrix(pN);
	DestroyDoubleMatrix(pS);
	DestroyDoubleMatrix(pX);
	DestroyDoubleMatrix(pV);

	/* save file */
	printf("########################################\n");
	printf("# Export Probe Model                   #\n");
	printf("########################################\n");
	pOldIntensity = NULL;
	pOldSD = NULL;
	pOldIntensity = vCELData[0]->pIntensity;
	pOldSD = vCELData[0]->pSD;
	vCELData[0]->pIntensity = pProbeMedian;
	vCELData[0]->pSD = pProbeInterR;

	sprintf(strOutFile, "%s.prbm", strOutputPath);
	Affy_SaveCELv4(strOutFile, vCELData[0]);

	vCELData[0]->pIntensity = pSortMean;
	sprintf(strOutFile, "%s.quan", strOutputPath);
	Affy_SaveCELv4(strOutFile, vCELData[0]);

	vCELData[0]->pIntensity = pOldIntensity;
	vCELData[0]->pSD = pOldSD;


	/* release memory */
	free(vGroupID);
	for(ni=0; ni<nArrayNum; ni++)
	{
		Affy_CELData_Destroy(vCELData+ni);
	}
	free(vCELData);

	DestroyDoubleMatrix(pSortMean);
	DestroyDoubleMatrix(pProbeMedian);
	DestroyDoubleMatrix(pProbeInterR);
	DestroyIntMatrix(pDf);
	DestroyDoubleMatrix(pProbeData);

	printf("Done!\n");

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_BuildHMMM_Main()                                             */
/*  Build a probe background model for microarray CEL files using HMM      */
/*  tiling algorithm with variance shrinking.                              */
/* ----------------------------------------------------------------------- */ 
int TileProbe_BuildHMMM_Main(char strInputPath[], char strOutputPath[], 
						 int nTakeLog, double dNormLowerBound, 
						 double dFilterTopPrc, int nTest)
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	FILE *fpOut2;
	double dNAValue = -1e10;

	/* arrays */
	int nArrayNum = 0;
	int *vGroupID;
	struct tagCELData **vCELData = NULL;
	char strFileName[MED_LINE_LENGTH];
	int ni,nj,nu;
	int nCELNumberCells = 0;
	int nCELTotalX = 0;
	int nCELTotalY = 0;
	int nFilterCutId = 0;

	struct DOUBLEMATRIX *pSortMean;
	struct DOUBLEMATRIX *pArray;
	struct LONGMATRIX *pSortId;

	double dQCut;
	struct DOUBLEMATRIX *pProbeMedian;
	struct DOUBLEMATRIX *pProbeInterR;
	struct INTMATRIX *pDf;
	struct DOUBLEMATRIX *pProbeData;
	struct DOUBLEMATRIX *pProbeDataSort;
	struct LONGMATRIX *pProbeDataSortId;

	struct INTMATRIX *pN;
	struct DOUBLEMATRIX *pS;
	struct DOUBLEMATRIX *pX;
	struct DOUBLEMATRIX *pV;

	double dM,dV,dV0,dGM,dB;
	int	nDf,nMaxDf,nD0,nGN;

	struct DOUBLEMATRIX *pOldIntensity;
	struct DOUBLEMATRIX *pOldSD;
	
	/* others */
	double dLog2 = log(2.0);
	int nId;
	int nQ50;
	char strOutFile[MED_LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	double dQmin,dQmax;

	/* get array number */
	nArrayNum = 0;
	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileProbe_BuildHMMB_Main, cannot open the input array list!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		nArrayNum++;
	}
	
	fclose(fpIn);
	
	printf("########################################\n");
	printf("# Loading Raw Data                     #\n");
	printf("########################################\n");
	printf(" Array number = %d\n", nArrayNum);

	if(nArrayNum == 0)
	{
		return PROC_SUCCESS;
	}

	/* memory for group id */
	vGroupID = NULL;
	vGroupID = (int *)calloc(nArrayNum, sizeof(int));
	if(vGroupID == NULL)
	{
		printf("Error: TileProbe_BuildHMMB_Main, cannot allocate memory for group id!\n");
		exit(EXIT_FAILURE);
	}

	/* load arrays */
	vCELData = NULL;
	vCELData = (struct tagCELData **)calloc(nArrayNum, sizeof(struct tagCELData *));
	if(vCELData == NULL)
	{
		printf("Error: TileProbe_BuildHMMB_Main, cannot allocate memory for the array list!\n");
		exit(EXIT_FAILURE);
	}

	/* store cutoffs */
	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileProbe_BuildHMMB_Main, cannot open the input array list!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%d %s", &nu, strFileName);
		vGroupID[ni] = nu;
		vCELData[ni] = TileMapv2_LoadCEL(strFileName);
		if(vCELData[ni] == NULL)
		{
			printf("Error: TileProbe_BuildHMMB_Main, cannot load *.CEL file!\n");
			exit(EXIT_FAILURE);
		}

		if(ni == 0)
		{
			nCELNumberCells = vCELData[ni]->nNumberCells;
			nCELTotalX = vCELData[ni]->nCols;
			nCELTotalY = vCELData[ni]->nRows;
			if( nCELNumberCells != nCELTotalX*nCELTotalY)
			{
				printf("Error: TileProbe_BuildHMMB_Main, CEL file Cols*Rows != NumberCells!\n");
				exit(EXIT_FAILURE);
			}
		}
		else
		{
			if( (nCELNumberCells != vCELData[ni]->nNumberCells) || (nCELTotalX != vCELData[ni]->nCols)
				|| (nCELTotalY != vCELData[ni]->nRows) )
			{
				printf("Error: TileProbe_BuildHMMT_Main, CEL file dimensions do not match!\n");
				exit(EXIT_FAILURE);
			}
		}

		if(nTakeLog == 1)
		{
			for(nj=0; nj<vCELData[ni]->nNumberCells; nj++)
			{
				/* truncate */
				if(vCELData[ni]->pIntensity->pMatElement[nj] < dNormLowerBound)
					vCELData[ni]->pIntensity->pMatElement[nj] = dNormLowerBound;

				/* log transformation */
				vCELData[ni]->pIntensity->pMatElement[nj] = log(vCELData[ni]->pIntensity->pMatElement[nj])/dLog2;
			}
		}
		else
		{
			for(nj=0; nj<vCELData[ni]->nNumberCells; nj++)
			{
				/* truncate */
				if(vCELData[ni]->pIntensity->pMatElement[nj] < dNormLowerBound)
					vCELData[ni]->pIntensity->pMatElement[nj] = dNormLowerBound;
			}
		}

		ni++;
	}
	
	fclose(fpIn);

	if(ni != nArrayNum)
	{
		printf("Error: TileProbe_BuildHMMB_Main, array number inconsistent!\n");
		exit(EXIT_FAILURE);
	}

	printf("########################################\n");
	printf("# Quantile Normalization               #\n");
	printf("########################################\n");
	
	pSortMean = NULL;
	pSortMean = CreateDoubleMatrix(1, nCELNumberCells);
	if(pSortMean == NULL)
	{
		printf("Error: TileProbe_BuildHMMB_Main, cannot allocate memory for quantile normalization!\n");
		exit(EXIT_FAILURE);
	}
				
	for(ni=0; ni<nArrayNum; ni++)
	{
		pArray = NULL;
		pSortId = NULL;

		/* export sorted values */
		printf("  Sorting array %d...\n", ni);
		DMSORTMERGEA_0(vCELData[ni]->pIntensity, &pArray, &pSortId);

		/* compute percentiles */
		for(nj=0; nj<nCELNumberCells; nj++)
		{
			pSortMean->pMatElement[nj] += pArray->pMatElement[nj];
			nId = pSortId->pMatElement[nj];
			vCELData[ni]->pIntensity->pMatElement[nId] = nj;
		}

		DestroyDoubleMatrix(pArray);
		DestroyLongMatrix(pSortId);
	}

	for(nj=0; nj<nCELNumberCells; nj++)
	{
		pSortMean->pMatElement[nj] /= nArrayNum;
	}

	nFilterCutId = nCELNumberCells-(int)(nCELNumberCells*dFilterTopPrc);
	if(nFilterCutId < 0)
		nFilterCutId = 0;
	if(nFilterCutId < nCELNumberCells)
		dQCut = pSortMean->pMatElement[nFilterCutId];
	else
		dQCut = 1e8;


	for(ni=0; ni<nArrayNum; ni++)
	{
		for(nj=0; nj<nCELNumberCells; nj++)
		{
			nId = (long)(vCELData[ni]->pIntensity->pMatElement[nj]);
			vCELData[ni]->pIntensity->pMatElement[nj] = pSortMean->pMatElement[nId];
		}
	}
	dNAValue = dNAValue+1.0;

	/* build probe model */
	printf("########################################\n");
	printf("# Build Probe Model                    #\n");
	printf("########################################\n");

	if(nTest == 1)
	{
		sprintf(strOutFile, "%s.testq", strOutputPath);
		fpOut = NULL;
		fpOut = fopen(strOutFile, "w");
		if(fpOut == NULL)
		{
			printf("Error: TileProbe_BuildHMMB_Main, cannot open the output argument file!\n");
			exit(EXIT_FAILURE);
		}
		sprintf(strOutFile, "%s.testprobe", strOutputPath);
		fpOut2 = NULL;
		fpOut2 = fopen(strOutFile, "w");
		if(fpOut2 == NULL)
		{
			printf("Error: TileProbe_BuildHMMB_Main, cannot open the output argument file!\n");
			exit(EXIT_FAILURE);
		}
	}

	pProbeMedian = NULL;
	pProbeMedian = CreateDoubleMatrix(1,nCELNumberCells);
	if(pProbeMedian == NULL)
	{
		printf("Error: TileProbe_BuildHMMB_Main, cannot allocate memory for probe median!\n");
		exit(EXIT_FAILURE);
	}

	pProbeInterR = NULL;
	pProbeInterR = CreateDoubleMatrix(1,nCELNumberCells);
	if(pProbeInterR == NULL)
	{
		printf("Error: TileProbe_BuildHMMB_Main, cannot allocate memory for probe interquatile range!\n");
		exit(EXIT_FAILURE);
	}

	pDf = NULL;
	pDf = CreateIntMatrix(1,nCELNumberCells);
	if(pDf == NULL)
	{
		printf("Error: TileProbe_BuildHMMB_Main, cannot allocate memory for storing degrees of freedom!\n");
		exit(EXIT_FAILURE);
	}

	pProbeData = NULL;
	pProbeData = CreateDoubleMatrix(1,nArrayNum);
	if(pProbeData == NULL)
	{
		printf("Error: TileProbe_BuildHMMB_Main, cannot allocate memory for modeling probe effect!\n");
		exit(EXIT_FAILURE);
	}

	nQ50 = (int)(nArrayNum*0.5);
	if( (nArrayNum*0.5)-nQ50 > 1e-6 )
		nQ50++;
	nQ50--;

	if(nQ50 < 0)
		nQ50 = 0;
	if(nQ50 >= nArrayNum)
		nQ50 = nArrayNum-1;

	dQmin = 1e10;
	dQmax = -1e10;
	nMaxDf = 0;
	for(nj=0; nj<nCELNumberCells; nj++)
	{
		if(nj%100000 == 0)
		{
			printf("%d ...\n", nj);
		}

		pProbeDataSort = NULL;
		pProbeDataSortId = NULL;

		dM = 0.0;
		dV = 0.0;
		nDf = 0;
		for(ni=0; ni<nArrayNum; ni++)
		{
			pProbeData->pMatElement[ni] = vCELData[ni]->pIntensity->pMatElement[nj];

			if(vGroupID[ni] <= 0)
			{
				nDf += 1;
				dM += pProbeData->pMatElement[ni];
			}
			else
			{
				if(pProbeData->pMatElement[ni] <= dQCut)
				{
					nDf += 1;
					dM += pProbeData->pMatElement[ni];
				}
			}
		}

		if(nDf <= 1)
		{
			pDf->pMatElement[nj] = 0;
			pProbeMedian->pMatElement[nj] = 0.0;
			pProbeInterR->pMatElement[nj] = 1e6;
			continue;
		}

		dM = dM/nDf;
		
		DMSORTMERGEA_0(pProbeData, &pProbeDataSort, &pProbeDataSortId);
		pProbeMedian->pMatElement[nj] = pProbeDataSort->pMatElement[nQ50];

		nDf = 0;
		for(ni=0; ni<nArrayNum; ni++)
		{
			if(vGroupID[ni] <= 0)
			{
				dV += (pProbeData->pMatElement[ni]-dM)*(pProbeData->pMatElement[ni]-dM);
				nDf += 1;
			}
			else
			{
				if(pProbeData->pMatElement[ni] <= dQCut)
				{
					dV += (pProbeData->pMatElement[ni]-dM)*(pProbeData->pMatElement[ni]-dM);
					nDf += 1;
				}
			}
		}
		nDf -= 1;
		dV /= nDf;
		pProbeInterR->pMatElement[nj] = dV;
		if(pProbeInterR->pMatElement[nj] < dQmin)
			dQmin = pProbeInterR->pMatElement[nj];
		if(pProbeInterR->pMatElement[nj] > dQmax)
			dQmax = pProbeInterR->pMatElement[nj];
		
		pDf->pMatElement[nj] = nDf;
		if(nDf > nMaxDf)
			nMaxDf = nDf;

		DestroyDoubleMatrix(pProbeDataSort);
		DestroyLongMatrix(pProbeDataSortId);
	}

	printf("min Var = %f; max Var = %f\n", dQmin, dQmax);

	if(nMaxDf <= 0)
	{
		printf("Error: TileProbe_BuildHMMB_Main, don't have enough degrees of freedom to estimate shrinkage factor!\n");
		exit(EXIT_FAILURE);
	}

	pN = NULL;
	pN = CreateIntMatrix(1, nMaxDf+1);

	pS = NULL;
	pS = CreateDoubleMatrix(1, nMaxDf+1);

	pX = NULL;
	pX = CreateDoubleMatrix(1, nMaxDf+1);

	pV = NULL;
	pV = CreateDoubleMatrix(1, nMaxDf+1);

	if( (pN == NULL) || (pS == NULL) || (pX == NULL) || (pV == NULL) )
	{
		printf("Error: TileProbe_BuildHMMB_Main, cannot allocate memory for computing shrinkage factor!\n");
		exit(EXIT_FAILURE);
	}

	dGM = 0.0;
	nGN = 0;
	for(nj=0; nj<nCELNumberCells; nj++)
	{
		nDf = pDf->pMatElement[nj];
		if(nDf <= 0)
			continue;

		pN->pMatElement[nDf] += 1;
		pX->pMatElement[nDf] += pProbeInterR->pMatElement[nj];
		pS->pMatElement[nDf] += pProbeInterR->pMatElement[nj]*pProbeInterR->pMatElement[nj];
		dGM += pProbeInterR->pMatElement[nj];
		nGN += 1;
	}
	dGM /= nGN;

	nD0 = 0;
	dV0 = 0.0;
	for(ni=1; ni<=nMaxDf; ni++)
	{
		if(pN->pMatElement[ni] < 10)
			continue;

		dM = pX->pMatElement[ni]/pN->pMatElement[ni];
		dV = pS->pMatElement[ni]-pN->pMatElement[ni]*dM*dM;
		dB = 2.0*(pN->pMatElement[ni]-1)/(2.0+ni)/pN->pMatElement[ni]+2.0*dM*dM*(pN->pMatElement[ni]-1)/(2.0+ni)/dV;
		if(dB>1.0)
			dB=1.0;
		else if(dB<0.0)
			dB=0.0;

		pV->pMatElement[ni] = (1.0/dB-1.0)/ni;
		printf("%d\t%f\n", pN->pMatElement[ni], pV->pMatElement[ni]); 

		dV0 += pV->pMatElement[ni]*pN->pMatElement[ni];
		nD0 += pN->pMatElement[ni];
	}

	if(nD0 <= 0)
	{
		printf("Error: TileProbe_BuildHMMB_Main, don't have enough degrees of freedom to estimate shrinkage factor!\n");
		exit(EXIT_FAILURE);
	}

	dV0 /= nD0;
	printf("Prior C = %f\n", dV0); 
	for(ni=0; ni<=nMaxDf; ni++)
	{
		pV->pMatElement[ni] = 1.0/(1.0+ni*dV0);
	}
		
	for(nj=0; nj<nCELNumberCells; nj++)
	{
		nDf = pDf->pMatElement[nj];
		dB = pV->pMatElement[nDf];

		pProbeInterR->pMatElement[nj] = sqrt(pProbeInterR->pMatElement[nj]*(1.0-dB)+dGM*dB);

		if(nTest == 1)
		{
			fprintf(fpOut, "%e\n", pProbeInterR->pMatElement[nj]); 
			if(rand_u() < 0.01)
			{
				fprintf(fpOut2, "%f", pProbeData->pMatElement[0]);
				for(ni=1; ni<nArrayNum; ni++)
				{
					fprintf(fpOut2, "\t%f", pProbeData->pMatElement[ni]);
				}
				fprintf(fpOut2, "\n");
			}
		}
	}

	if(nTest == 1)
	{
		fclose(fpOut);
		fclose(fpOut2);
	}

	DestroyIntMatrix(pN);
	DestroyDoubleMatrix(pS);
	DestroyDoubleMatrix(pX);
	DestroyDoubleMatrix(pV);

	/* save file */
	printf("########################################\n");
	printf("# Export Probe Model                   #\n");
	printf("########################################\n");
	pOldIntensity = NULL;
	pOldSD = NULL;
	pOldIntensity = vCELData[0]->pIntensity;
	pOldSD = vCELData[0]->pSD;
	vCELData[0]->pIntensity = pProbeMedian;
	vCELData[0]->pSD = pProbeInterR;

	sprintf(strOutFile, "%s.prbm", strOutputPath);
	Affy_SaveCELv4(strOutFile, vCELData[0]);

	vCELData[0]->pIntensity = pSortMean;
	sprintf(strOutFile, "%s.quan", strOutputPath);
	Affy_SaveCELv4(strOutFile, vCELData[0]);

	vCELData[0]->pIntensity = pOldIntensity;
	vCELData[0]->pSD = pOldSD;


	/* release memory */
	free(vGroupID);
	for(ni=0; ni<nArrayNum; ni++)
	{
		Affy_CELData_Destroy(vCELData+ni);
	}
	free(vCELData);

	DestroyDoubleMatrix(pSortMean);
	DestroyDoubleMatrix(pProbeMedian);
	DestroyDoubleMatrix(pProbeInterR);
	DestroyIntMatrix(pDf);
	DestroyDoubleMatrix(pProbeData);

	printf("Done!\n");

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_BuildHMMW_Main()                                             */
/*  Build a probe background model for microarray CEL files using HMM      */
/*  tiling algorithm with variance shrinking.                              */
/* ----------------------------------------------------------------------- */ 
int TileProbe_BuildHMMW_Main(char strInputPath[], char strOutputPath[], 
						 int nTakeLog, double dNormLowerBound, 
						 double dFilterTopPrc, int nShrink, int nTest)
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	FILE *fpOut2;
	double dNAValue = -1e10;

	/* arrays */
	int nArrayNum = 0;
	int *vGroupID;
	struct tagCELData **vCELData = NULL;
	char strFileName[MED_LINE_LENGTH];
	int ni,nj,nu,nx,ny,nz;
	int nCELNumberCells = 0;
	int nCELTotalX = 0;
	int nCELTotalY = 0;
	int nFilterCutId = 0;

	struct DOUBLEMATRIX *pSortMean;
	struct DOUBLEMATRIX *pArray;
	struct LONGMATRIX *pSortId;

	double dQCut;
	struct DOUBLEMATRIX *pProbeMedian;
	struct DOUBLEMATRIX *pProbeInterR;
	struct DOUBLEMATRIX *pProbeData;
	struct DOUBLEMATRIX *pProbeDataSort;
	struct LONGMATRIX *pProbeDataSortId;

	double dM,dV,dB;
	int	nDf,nMaxDf,nGN;

	struct DOUBLEMATRIX *pOldIntensity;
	struct DOUBLEMATRIX *pOldSD;
	
	/* others */
	double dLog2 = log(2.0);
	int nId;
	int nQ50;
	char strOutFile[MED_LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	double dQmin,dQmax;

	/* get array number */
	nArrayNum = 0;
	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileProbe_BuildHMMB_Main, cannot open the input array list!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		nArrayNum++;
	}
	
	fclose(fpIn);
	
	printf("########################################\n");
	printf("# Loading Raw Data                     #\n");
	printf("########################################\n");
	printf(" Array number = %d\n", nArrayNum);

	if(nArrayNum == 0)
	{
		return PROC_SUCCESS;
	}

	/* memory for group id */
	vGroupID = NULL;
	vGroupID = (int *)calloc(nArrayNum, sizeof(int));
	if(vGroupID == NULL)
	{
		printf("Error: TileProbe_BuildHMMB_Main, cannot allocate memory for group id!\n");
		exit(EXIT_FAILURE);
	}

	/* load arrays */
	vCELData = NULL;
	vCELData = (struct tagCELData **)calloc(nArrayNum, sizeof(struct tagCELData *));
	if(vCELData == NULL)
	{
		printf("Error: TileProbe_BuildHMMB_Main, cannot allocate memory for the array list!\n");
		exit(EXIT_FAILURE);
	}

	/* store cutoffs */
	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileProbe_BuildHMMB_Main, cannot open the input array list!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%d %s", &nu, strFileName);
		vGroupID[ni] = nu;
		vCELData[ni] = TileMapv2_LoadCEL(strFileName);
		if(vCELData[ni] == NULL)
		{
			printf("Error: TileProbe_BuildHMMB_Main, cannot load *.CEL file!\n");
			exit(EXIT_FAILURE);
		}

		if(ni == 0)
		{
			nCELNumberCells = vCELData[ni]->nNumberCells;
			nCELTotalX = vCELData[ni]->nCols;
			nCELTotalY = vCELData[ni]->nRows;
			if( nCELNumberCells != nCELTotalX*nCELTotalY)
			{
				printf("Error: TileProbe_BuildHMMB_Main, CEL file Cols*Rows != NumberCells!\n");
				exit(EXIT_FAILURE);
			}
		}
		else
		{
			if( (nCELNumberCells != vCELData[ni]->nNumberCells) || (nCELTotalX != vCELData[ni]->nCols)
				|| (nCELTotalY != vCELData[ni]->nRows) )
			{
				printf("Error: TileProbe_BuildHMMT_Main, CEL file dimensions do not match!\n");
				exit(EXIT_FAILURE);
			}
		}

		if(nTakeLog == 1)
		{
			for(nj=0; nj<vCELData[ni]->nNumberCells; nj++)
			{
				/* truncate */
				if(vCELData[ni]->pIntensity->pMatElement[nj] < dNormLowerBound)
					vCELData[ni]->pIntensity->pMatElement[nj] = dNormLowerBound;

				/* log transformation */
				vCELData[ni]->pIntensity->pMatElement[nj] = log(vCELData[ni]->pIntensity->pMatElement[nj])/dLog2;
			}
		}
		else
		{
			for(nj=0; nj<vCELData[ni]->nNumberCells; nj++)
			{
				/* truncate */
				if(vCELData[ni]->pIntensity->pMatElement[nj] < dNormLowerBound)
					vCELData[ni]->pIntensity->pMatElement[nj] = dNormLowerBound;
			}
		}

		ni++;
	}
	
	fclose(fpIn);

	if(ni != nArrayNum)
	{
		printf("Error: TileProbe_BuildHMMB_Main, array number inconsistent!\n");
		exit(EXIT_FAILURE);
	}

	printf("########################################\n");
	printf("# Quantile Normalization               #\n");
	printf("########################################\n");
	
	pSortMean = NULL;
	pSortMean = CreateDoubleMatrix(1, nCELNumberCells);
	if(pSortMean == NULL)
	{
		printf("Error: TileProbe_BuildHMMB_Main, cannot allocate memory for quantile normalization!\n");
		exit(EXIT_FAILURE);
	}
				
	for(ni=0; ni<nArrayNum; ni++)
	{
		pArray = NULL;
		pSortId = NULL;

		/* export sorted values */
		printf("  Sorting array %d...\n", ni);
		DMSORTMERGEA_0(vCELData[ni]->pIntensity, &pArray, &pSortId);

		/* compute percentiles */
		for(nj=0; nj<nCELNumberCells; nj++)
		{
			pSortMean->pMatElement[nj] += pArray->pMatElement[nj];
			nId = pSortId->pMatElement[nj];
			vCELData[ni]->pIntensity->pMatElement[nId] = nj;
		}

		DestroyDoubleMatrix(pArray);
		DestroyLongMatrix(pSortId);
	}

	for(nj=0; nj<nCELNumberCells; nj++)
	{
		pSortMean->pMatElement[nj] /= nArrayNum;
	}

	nFilterCutId = nCELNumberCells-(int)(nCELNumberCells*dFilterTopPrc);
	if(nFilterCutId < 0)
		nFilterCutId = 0;
	if(nFilterCutId < nCELNumberCells)
		dQCut = pSortMean->pMatElement[nFilterCutId];
	else
		dQCut = 1e8;


	for(ni=0; ni<nArrayNum; ni++)
	{
		for(nj=0; nj<nCELNumberCells; nj++)
		{
			nId = (long)(vCELData[ni]->pIntensity->pMatElement[nj]);
			vCELData[ni]->pIntensity->pMatElement[nj] = pSortMean->pMatElement[nId];
		}
	}
	dNAValue = dNAValue+1.0;

	/* build probe model */
	printf("########################################\n");
	printf("# Build Probe Model                    #\n");
	printf("########################################\n");

	if(nTest == 1)
	{
		sprintf(strOutFile, "%s.testq", strOutputPath);
		fpOut = NULL;
		fpOut = fopen(strOutFile, "w");
		if(fpOut == NULL)
		{
			printf("Error: TileProbe_BuildHMMB_Main, cannot open the output argument file!\n");
			exit(EXIT_FAILURE);
		}
		sprintf(strOutFile, "%s.testprobe", strOutputPath);
		fpOut2 = NULL;
		fpOut2 = fopen(strOutFile, "w");
		if(fpOut2 == NULL)
		{
			printf("Error: TileProbe_BuildHMMB_Main, cannot open the output argument file!\n");
			exit(EXIT_FAILURE);
		}
	}

	pProbeMedian = NULL;
	pProbeMedian = CreateDoubleMatrix(1,nCELNumberCells);
	if(pProbeMedian == NULL)
	{
		printf("Error: TileProbe_BuildHMMB_Main, cannot allocate memory for probe median!\n");
		exit(EXIT_FAILURE);
	}

	pProbeInterR = NULL;
	pProbeInterR = CreateDoubleMatrix(1,nCELNumberCells);
	if(pProbeInterR == NULL)
	{
		printf("Error: TileProbe_BuildHMMB_Main, cannot allocate memory for probe interquatile range!\n");
		exit(EXIT_FAILURE);
	}

	pProbeData = NULL;
	pProbeData = CreateDoubleMatrix(1,nArrayNum);
	if(pProbeData == NULL)
	{
		printf("Error: TileProbe_BuildHMMB_Main, cannot allocate memory for modeling probe effect!\n");
		exit(EXIT_FAILURE);
	}

	nQ50 = (int)(nArrayNum*0.5);
	if( (nArrayNum*0.5)-nQ50 > 1e-6 )
		nQ50++;
	nQ50--;

	if(nQ50 < 0)
		nQ50 = 0;
	if(nQ50 >= nArrayNum)
		nQ50 = nArrayNum-1;

	dQmin = 1e10;
	dQmax = -1e10;
	nMaxDf = 0;
	for(nj=0; nj<nCELNumberCells; nj++)
	{
		if(nj%100000 == 0)
		{
			printf("%d ...\n", nj);
		}

		pProbeDataSort = NULL;
		pProbeDataSortId = NULL;

		for(ni=0; ni<nArrayNum; ni++)
		{
			pProbeData->pMatElement[ni] = vCELData[ni]->pIntensity->pMatElement[nj];
		}
		
		nx = 0;
		dV = 0.0;
		nDf = 0;
		while(nx < nArrayNum)
		{
			dM = 0.0;
			
			nGN = 0;
			for(ny=nx; ny<nArrayNum; ny++)
			{
				if(vGroupID[ny] != vGroupID[nx])
					break;
				dM += pProbeData->pMatElement[ny];
				nGN++;
			}
			dM = dM/nGN;
			for(nz=nx; nz<ny; nz++)
			{
				dV += (pProbeData->pMatElement[nz]-dM)*(pProbeData->pMatElement[nz]-dM);
			}
			nDf = nDf+nGN-1;
			nx = ny;
		}
		dV /= nDf;
		pProbeInterR->pMatElement[nj] = dV;
		nMaxDf = nDf;


		dM = 0.0;
		nDf = 0;
		for(ni=0; ni<nArrayNum; ni++)
		{
			if(vGroupID[ni] <= 0)
			{
				nDf += 1;
				dM += pProbeData->pMatElement[ni];
			}
			else
			{
				if(pProbeData->pMatElement[ni] <= dQCut)
				{
					nDf += 1;
					dM += pProbeData->pMatElement[ni];
				}
			}
		}

		if(nDf <= 0)
		{
			pProbeMedian->pMatElement[nj] = 0.0;
			pProbeInterR->pMatElement[nj] = 1e6;
			continue;
		}

		dM = dM/nDf;
		pProbeMedian->pMatElement[nj] = dM;
	}

	nDf = nMaxDf;

	if(nShrink != 0)
	{
		dM = 0.0;
		dV = 0.0;
		dQmin = 1e10;
		dQmax = -1e10;

		for(nj=0; nj<nCELNumberCells; nj++)
		{
			if(pProbeInterR->pMatElement[nj] < dQmin)
				dQmin = pProbeInterR->pMatElement[nj];
			if(pProbeInterR->pMatElement[nj] > dQmax)
				dQmax = pProbeInterR->pMatElement[nj];

			dM += pProbeInterR->pMatElement[nj];
			dV += pProbeInterR->pMatElement[nj]*pProbeInterR->pMatElement[nj];
		}

		dM /= nCELNumberCells;
		dV = dV-nCELNumberCells*dM*dM;
		dB = 2.0*(nCELNumberCells-1)/(2.0+nDf)/nCELNumberCells+2.0*dM*dM*(nCELNumberCells-1)/(2.0+nDf)/dV;
		
		if(dB>1.0)
			dB=1.0;
		else if(dB<0.0)
			dB=0.0;

		if(nShrink == -1)
			dB = 0.0;

		for(nj=0; nj<nCELNumberCells; nj++)
		{
			pProbeInterR->pMatElement[nj] = sqrt(pProbeInterR->pMatElement[nj]*(1.0-dB)+dM*dB);
		}
		
		printf("min Var = %f; max Var = %f\n", dQmin, dQmax);
		printf("B = %f\n", dB);
	}

	if(nTest == 1)
	{
		fclose(fpOut);
		fclose(fpOut2);
	}

	/* save file */
	printf("########################################\n");
	printf("# Export Probe Model                   #\n");
	printf("########################################\n");
	pOldIntensity = NULL;
	pOldSD = NULL;
	pOldIntensity = vCELData[0]->pIntensity;
	pOldSD = vCELData[0]->pSD;
	vCELData[0]->pIntensity = pProbeMedian;
	vCELData[0]->pSD = pProbeInterR;

	sprintf(strOutFile, "%s.prbm", strOutputPath);
	Affy_SaveCELv4(strOutFile, vCELData[0]);

	vCELData[0]->pIntensity = pSortMean;
	sprintf(strOutFile, "%s.quan", strOutputPath);
	Affy_SaveCELv4(strOutFile, vCELData[0]);

	vCELData[0]->pIntensity = pOldIntensity;
	vCELData[0]->pSD = pOldSD;


	/* release memory */
	free(vGroupID);
	for(ni=0; ni<nArrayNum; ni++)
	{
		Affy_CELData_Destroy(vCELData+ni);
	}
	free(vCELData);

	DestroyDoubleMatrix(pSortMean);
	DestroyDoubleMatrix(pProbeMedian);
	DestroyDoubleMatrix(pProbeInterR);
	DestroyDoubleMatrix(pProbeData);

	printf("Done!\n");

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_HuberData_Main()                                             */
/*  Get data for testing Huber et al.                                      */
/* ----------------------------------------------------------------------- */ 
int TileProbe_HuberData_Main(char strParamPath[])
{
	/* define */

	/* working path */
	char strProjectTitle[MED_LINE_LENGTH];
	char strCELPath[MED_LINE_LENGTH];
	char strBpmapPath[MED_LINE_LENGTH];
	char strWorkPath[MED_LINE_LENGTH];
	char strMaskPath[MED_LINE_LENGTH];
	char strGenomeGrp[255];
	int nIncludeNonGrp = 0;
	
	/* arrays */
	int nLibNum = 0;
	struct tagString **vLibName = NULL;
	int nSampleNum = 0;
	int nArrayNum = 0;
	struct tagString **vCELPath;
	struct tagString **vArrayPath;
	struct tagString **vAlias;
	
	/* CEL files */
	int nTotalProbeNum = 0;
	int nRealProbeNum = 0;
	int nCELNumberCells = 0;
	int nCELTotalX = 0;
	int nCELTotalY = 0;
	struct tagCELData **vCELData;
	
	/* mask */
	int nRemoveMaskedCells = 0;
	int nRemoveOutlierCells = 0;
	struct INTMATRIX *pNumMaskCells;

	/* normalization */
	double dNormLowerBound = 1.0;
	int nNormLogTransform = 1;
	double dLog2 = log(2.0);
	
	/* intensity computation */
	struct tagBARData *pBARPos;
	int nExportMode = 0;

	/* others */
	int ni,nj,nk,nx,ny;
	char strLine[LONG_LINE_LENGTH];
	
	FILE *fpIn;
	char *chSep;
	int nError = 0;
	FILE *fpOut;
	char strOutPath[MED_LINE_LENGTH];
	int nidx1;

	/* init */
	strcpy(strGenomeGrp, "");

	/* load parameters */
	fpIn = NULL;
	fpIn = fopen(strParamPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileProbe_HuberData_Main, cannot open the parameter file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strstr(strLine, "[Project title]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileProbe_HuberData_Main, cannot load project title!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strProjectTitle, chSep);
		}
		else if(strstr(strLine, "[CEL directory]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileProbe_HuberData_Main, cannot load CEL directory!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strCELPath, chSep);
			AdjustDirectoryPath(strCELPath);
		}
		else if(strstr(strLine, "[BPMAP directory]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileProbe_HuberData_Main, cannot load BPMAP directory!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strBpmapPath, chSep);
			AdjustDirectoryPath(strBpmapPath);
		}
		else if(strstr(strLine, "[GenomeGrp]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileProbe_HuberData_Main, cannot load Genome Group!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strGenomeGrp, chSep);
		}
		else if(strstr(strLine, "[Include Probes not in GenomeGrp]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileProbe_MAT_Main, cannot load Genome Group!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			nIncludeNonGrp = atoi(chSep);
		}
		else if(strstr(strLine, "[Working directory]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileProbe_HuberData_Main, cannot load working directory!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strWorkPath, chSep);
			AdjustDirectoryPath(strWorkPath);
		}
		
		else if(strstr(strLine, "[No. of Libraries]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileProbe_HuberData_Main, cannot load no. of libraries!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			nLibNum = atoi(chSep);
			if(nLibNum <= 0)
			{
				printf("Warning: No BPMAP libraries provided!");
				return PROC_SUCCESS;
			}

			vLibName = NULL;
			vLibName = (struct tagString **)calloc(nLibNum, sizeof(struct tagString *));
			if(vLibName == NULL)
			{
				printf("Error: TileProbe_HuberData_Main, cannot allocate memory for loading BPMAP lists!\n");
				exit(EXIT_FAILURE);
			}
		}
		else if(strstr(strLine, "[Libraries]") == strLine)
		{
			ni = 0;
			while(ni < nLibNum)
			{
				fgets(strLine, LONG_LINE_LENGTH, fpIn);
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] == '\0')
					continue;

				StringAddTail(vLibName+ni, strLine);
				ni++;
			}
		}
		else if(strstr(strLine, "[No. of Samples]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileProbe_HuberData_Main, cannot load no. of samples!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			nSampleNum = atoi(chSep);
            if(nSampleNum <= 0)
			{
				printf("Error: TileProbe_HuberData_Main, no arrays available!\n");
				exit(EXIT_FAILURE);
			}

			nArrayNum = (int)(nSampleNum*nLibNum);

			vCELPath = NULL;
			vCELPath = (struct tagString **)calloc(nArrayNum, sizeof(struct tagString *));
			if(vCELPath == NULL)
			{
				printf("Error: TileProbe_HuberData_Main, cannot allocate memory for tracking *.CEL files!\n");
				exit(EXIT_FAILURE);
			}

			vArrayPath = NULL;
			vArrayPath = (struct tagString **)calloc(nArrayNum, sizeof(struct tagString *));
			if(vArrayPath == NULL)
			{
				printf("Error: TileProbe_HuberData_Main, cannot allocate memory for tracking array files!\n");
				exit(EXIT_FAILURE);
			}

			vAlias = NULL;
			vAlias = (struct tagString **)calloc(nSampleNum, sizeof(struct tagString *));
			if(vAlias == NULL)
			{
				printf("Error: TileProbe_HuberData_Main, cannot allocate memory for tracking *.CEL files!\n");
				exit(EXIT_FAILURE);
			}

			pNumMaskCells = NULL;
			pNumMaskCells = CreateIntMatrix(1,nArrayNum);
			if(pNumMaskCells == NULL)
			{
				printf("Error: TileProbe_HuberData_Main, cannot allocate memory for tracking number of masked cells!\n");
				exit(EXIT_FAILURE);
			}
		}
		
		else if(strstr(strLine, "[Arrays]") == strLine)
		{
			ni = 0;
			while(ni < nSampleNum)
			{
				fgets(strLine, LONG_LINE_LENGTH, fpIn);
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] == '\0')
					continue;
				
				if(strLine[0] != '>')
				{
					printf("Error: TileProbe_MAT_Main, error when loading samples!\n");
					exit(EXIT_FAILURE);
				}

				chSep = strLine+1;
				StrTrimLeft(chSep);
				StringAddTail(vAlias+ni, chSep);
				
				nj = 0;
				while(nj < nLibNum)
				{
					fgets(strLine, LONG_LINE_LENGTH, fpIn);
					StrTrimLeft(strLine);
					StrTrimRight(strLine);
					if(strLine[0] == '\0')
						continue;

					StringAddTail(vCELPath+(ni*nLibNum)+nj, strLine);
					chSep = NULL;
					chSep = strrchr(strLine, '.');
					if(chSep != NULL)
						*chSep = '\0';
					StringAddTail(vArrayPath+(ni*nLibNum)+nj, strLine);
					nj++;
				}
				ni++;
			}
		}

		else if(strstr(strLine, "[Remove masked cells in the CEL files]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileProbe_HuberData_Main, cannot load masked cell option!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			nRemoveMaskedCells = atoi(chSep);
		}

		else if(strstr(strLine, "[Remove outlier cells in the CEL files]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("Error: TileProbe_HuberData_Main, cannot load masked cell option!\n");
				exit(EXIT_FAILURE);
			}
			chSep++;
			StrTrimLeft(chSep);
			nRemoveOutlierCells = atoi(chSep);
		}
		else if(strLine[0] == '#')
		{
		}

		else
		{
			printf("Error: TileProbe_HuberData_Main, unknown parameters!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	fclose(fpIn);

	printf("########################################\n");
	printf("# MAT background correction            #\n");
	printf("########################################\n");

	for(ni=0; ni<nLibNum; ni++)
	{
		printf("Processing %s...\n", vLibName[ni]->m_pString);
		sprintf(strLine, "%s%s", strBpmapPath, vLibName[ni]->m_pString);
		sprintf(strMaskPath, "%s%s.refmask", strWorkPath, vLibName[ni]->m_pString);
		
		/* load bpmap */
		pBARPos = NULL;
		pBARPos = TileProbe_BpmapToBARHuber(strLine, strMaskPath, strGenomeGrp, nIncludeNonGrp);
		if(pBARPos == NULL)
		{
			printf("Error: TileProbe_HuberData_Main, empty bpmap file!\n");
			exit(EXIT_FAILURE);
		}

		vCELData = NULL;
		vCELData = (struct tagCELData **)calloc(nSampleNum, sizeof(struct tagCELData*));
		if(vCELData == NULL)
		{
			printf("Error: TileProbe_HuberData_Main, cannot create memory for cel data!\n");
			exit(EXIT_FAILURE);
		}
		

		/* load cel files */
		for(nj=0; nj<nSampleNum; nj++)
		{
			nk = nj*nLibNum+ni;
			printf("Loading %s \n", vCELPath[nk]->m_pString);
			sprintf(strLine, "%s%s", strCELPath, vCELPath[nk]->m_pString);

			/* load CEL */
			vCELData[nj] = TileMapv2_LoadCEL(strLine);
			if(vCELData[nj] == NULL)
			{
				printf("Error: TileProbe_HuberData_Main, cannot load *.CEL file!\n");
				exit(EXIT_FAILURE);
			}

			if(nj == 0)
			{
				nCELTotalX = vCELData[nj]->nCols;
				nCELTotalY = vCELData[nj]->nRows;
				nCELNumberCells = vCELData[nj]->nNumberCells;
				if( (int)(nCELTotalX*nCELTotalY) != nCELNumberCells)
				{
					printf("Error: TileProbe_HuberData_Main, CEL dimension not match!\n");
					exit(EXIT_FAILURE);
				}
			}
			else
			{
				if( (nCELTotalX != vCELData[nj]->nCols) || (nCELTotalY != vCELData[nj]->nRows)) 
				{
					printf("Error: TileProbe_HuberData_Main, CEL dimension not match!\n");
					exit(EXIT_FAILURE);
				}
			}
		}

		/* write data */
		sprintf(strOutPath, "%s%s_%d.txt", strWorkPath, strProjectTitle, ni+1);
		fpOut = NULL;
		fpOut = fopen(strOutPath, "w");
		if(fpOut == NULL)
		{
			printf("Error: TileProbe_HuberData_Main, cannot open output file!\n");
			exit(EXIT_FAILURE);
		}

		/* fprintf(fpOut, "#ProbeID"); */
		nk = 0*nLibNum+ni;
		fprintf(fpOut, "%s", vCELPath[nk]->m_pString);
		for(nj=1; nj<nSampleNum; nj++)
		{
			nk = nj*nLibNum+ni;
			fprintf(fpOut, "\t%s", vCELPath[nk]->m_pString);
		}
		fprintf(fpOut, "\n");

		for(nx=0; nx<pBARPos->nSeqNum; nx++)
		{
			if(nIncludeNonGrp == 0)
			{
				if(pBARPos->vSeqData[nx]->pSeqGroupName == NULL)
				{
					if(strcmp(strGenomeGrp, "") != 0)
					{
						continue;
					}
				}
				else if(strcmp(pBARPos->vSeqData[nx]->pSeqGroupName->m_pString, strGenomeGrp) != 0)
				{
					continue;
				}
			}

			for(ny=0; ny<pBARPos->vSeqData[nx]->nDataNum; ny++)
			{
				if( (pBARPos->vSeqData[nx]->vData[3]->pMatElement[ny] > 1.5) || (pBARPos->vSeqData[nx]->vData[3]->pMatElement[ny] < 0.5))
					continue;

				fprintf(fpOut, "%s:%d", pBARPos->vSeqData[nx]->pSeqName->m_pString, (int)(pBARPos->vSeqData[nx]->vData[0]->pMatElement[ny]));

				/* get coordinate */
				nidx1 = (int)(pBARPos->vSeqData[nx]->vData[2]->pMatElement[ny])*nCELTotalX+(int)(pBARPos->vSeqData[nx]->vData[1]->pMatElement[ny]);
				if(nidx1 >= nCELNumberCells)
				{
					printf("Error: index out of range\n");
					exit(EXIT_FAILURE);
				}

				for(nj=0; nj<nSampleNum; nj++)
				{
					fprintf(fpOut, "\t%f", vCELData[nj]->pIntensity->pMatElement[nidx1]);
				}
				fprintf(fpOut, "\n");
			}
		}

		fclose(fpOut);

		/* release memory */
		for(nj=0; nj<nSampleNum; nj++)
		{
			Affy_CELData_Destroy(vCELData+nj);
		}
		free(vCELData);

		Affy_BARData_Destroy(&pBARPos);
	}


	/* destroy */
	for(ni=0; ni<nArrayNum; ni++)
	{
		DeleteString(vCELPath[ni]);
		DeleteString(vArrayPath[ni]);
	}
	free(vCELPath);
	free(vArrayPath);

	for(ni=0; ni<nSampleNum; ni++)
	{
		DeleteString(vAlias[ni]);
	}
	free(vAlias);

	for(ni=0; ni<nLibNum; ni++)
	{
		DeleteString(vLibName[ni]);
	}
	free(vLibName);

	DestroyIntMatrix(pNumMaskCells);

	printf("Done!\n");

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_BpmapToBARv2()                                               */
/*  Create BAR template from bpmap file. The function will create local    */
/*  repeat mask at the same time.                                          */
/* ----------------------------------------------------------------------- */
struct tagBARData *TileProbe_BpmapToBARHuber(char strBpmapFile[], char strMaskPath[],
			char strGenomeGrp[], int nIncludeNonGrp)
{
	/* define */
	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));
	double dLog2 = log(2.0);

	struct tagBARData *pBARPos = NULL;

	FILE *fpIn;
	FILE *fpOut;

	/* variables */
	char strFileType[9];
	float fVersion;
	unsigned int nSeqNum;

	/* seq info */
	unsigned int nSeqNameLen;

	struct INTMATRIX *vProbeMappingType;
	unsigned int nProbeMappingType;
	
	struct INTMATRIX *vSequenceFileOffset;
	unsigned int nSequenceFileOffset;
	
	struct INTMATRIX *vProbePairNum;
	unsigned int nProbePairNum;
	
	unsigned int nParamNum;
		
	/* map info */
	struct INTMATRIX *vSeqID;
	unsigned int nSeqID;
	int nProbeNum;
	int nInfoCol;
	int nTotalProbeNum = 0;
	int nMaskedProbeNum = 0;
	struct tagAffyBpMapUnit *pNewUnit,*pCUnit,*pPUnit;
	struct tagAffyBpMapUnit *pUnitList;
	struct DOUBLEMATRIX *pResizeMat;
	int nIgnore = 0;
	int nDivide1M = 0;

	/* count */
	int ni,nj,nk;
	int nEndPos;
	
	/* load */
	fpIn = NULL;
	fpIn = fopen(strBpmapFile, "rb");
	if(fpIn == NULL)
	{
		printf("Error: TileProbe_BpmapToBAR, cannot open .bpmap file!\n");
		exit(EXIT_FAILURE);
	}
	fpOut = NULL;
	fpOut = fopen(strMaskPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: TileProbe_BpmapToBAR, cannot open .refmask file!\n");
		exit(EXIT_FAILURE);
	}
    fprintf(fpOut, "chromosome\tposition\tprobe_num\trepeat_num\tprobe_seq\tPMx\tPMy\tMMx\tMMy\n");
	
	/* load head */
	fread( strFileType, sizeof(char), 8, fpIn );
	strFileType[8] = '\0';
	AFFYBAR_READ_VERSION(fpIn, &fVersion);
	if(big_endian_fread(&nSeqNum, 4, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: TileProbe_BpmapToBAR, cannot load sequence number.\n");
		exit(EXIT_FAILURE);
	}
	printf("  BPMAP Version = %f\n", fVersion);
	printf("  Number of Sequences = %d\n", nSeqNum);
	if(nSeqNum <= 0)
	{
		printf("Warning: empty sequences!\n");
		return NULL;
	}

	/* create BAR object */
	pBARPos = Affy_BARData_Create();
	if(pBARPos == NULL)
	{
		printf("Error: TileProbe_BpmapToBAR, BARData object was not created!\n");
		exit(EXIT_FAILURE);
	}
	strcpy(pBARPos->strMagicnumber, "barr\r\n\032\n");
	pBARPos->fVersionnumber = 2.0;
    pBARPos->nSeqNum = nSeqNum;
	/* columns: 1 for coord, 1 for PMX, 1 for PMY, 1 for copy number */
	pBARPos->nColNum = 4;
	nInfoCol = 4;
	pBARPos->pFieldType = CreateIntMatrix(1, pBARPos->nColNum);
	if(pBARPos->pFieldType == NULL)
	{
		printf("Error: TileProbe_BpmapToBAR, cannot allocate memory for saving field type.\n");
		exit(EXIT_FAILURE);
	}
	pBARPos->pFieldType->pMatElement[0] = 2;
	pBARPos->pFieldType->pMatElement[1] = 2;
	pBARPos->pFieldType->pMatElement[2] = 2;
	pBARPos->pFieldType->pMatElement[3] = 1;
	
	pBARPos->vSeqData = (struct tagBARSeq **)calloc(pBARPos->nSeqNum, sizeof(struct tagBARSeq *));
	if(pBARPos->vSeqData == NULL)
	{
		printf("Error: TileProbe_BpmapToBAR, cannot allocate memory for loading sequence data.\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<pBARPos->nSeqNum; ni++)
	{
		/* create BARSeq object */
		pBARPos->vSeqData[ni] = Affy_BARSeq_Create();
		if(pBARPos->vSeqData[ni] == NULL)
		{
			printf("Error: TileProbe_BpmapToBAR, cannot create BARSeq object.\n");
			exit(EXIT_FAILURE);
		}
	}

	vProbeMappingType = NULL;
	vProbeMappingType = CreateIntMatrix(nSeqNum,1);
	if(vProbeMappingType == NULL)
	{
		printf("Error: TileProbe_BpmapToBAR, cannot load probe mapping type!\n");
		exit(EXIT_FAILURE);
	}

	if(fVersion > 2.5)
	{
		vSequenceFileOffset = NULL;
		vSequenceFileOffset = CreateIntMatrix(nSeqNum,1);
		if(vSequenceFileOffset == NULL)
		{
			printf("Error: TileProbe_BpmapToBAR, cannot load sequence file offset!\n");
			exit(EXIT_FAILURE);
		}
	}

	vProbePairNum = NULL;
	vProbePairNum = CreateIntMatrix(nSeqNum,1);
	if(vProbePairNum == NULL)
	{
		printf("Error: TileProbe_BpmapToBAR, cannot load probe pair number!\n");
		exit(EXIT_FAILURE);
	}
	
	for(ni=0; ni<(int)nSeqNum; ni++)
	{
		/* for version 1.0 or later */
		/* load sequence name */
		if(big_endian_fread(&nSeqNameLen, 4, 1, fpIn, little_endian_machine) != 1)
		{
			printf("Error: TileProbe_BpmapToBAR, cannot load sequence name length.\n");
			exit(EXIT_FAILURE);
		}
		if(nSeqNameLen > 0)
		{
			pBARPos->vSeqData[ni]->pSeqName = CreateString(nSeqNameLen);
			if(pBARPos->vSeqData[ni]->pSeqName == NULL)
			{
				printf("Error: TileProbe_BpmapToBAR, cannot load sequence name.\n");
				exit(EXIT_FAILURE);
			}
			if(big_endian_fread(pBARPos->vSeqData[ni]->pSeqName->m_pString, 1, nSeqNameLen, fpIn, little_endian_machine) != nSeqNameLen)
			{
				printf("Error: TileProbe_BpmapToBAR, cannot load sequence name.\n");
				exit(EXIT_FAILURE);
			}
			pBARPos->vSeqData[ni]->pSeqName->m_pString[nSeqNameLen] = '\0';
		}

		/* for version 3.0 or later */
		/* load probe mapping type and sequence file offset */
		if(fVersion > 2.5)
		{
			if(big_endian_fread(&nProbeMappingType, 4, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: TileProbe_BpmapToBAR, cannot load probe mapping type.\n");
				exit(EXIT_FAILURE);
			}
			vProbeMappingType->pMatElement[ni] = (int)nProbeMappingType;
			if(big_endian_fread(&nSequenceFileOffset, 4, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: TileProbe_BpmapToBAR, cannot load sequence file offset.\n");
				exit(EXIT_FAILURE);
			}
			vSequenceFileOffset->pMatElement[ni] = (int)nSequenceFileOffset;
		}

		/* for version 1.0 or later */
		if(big_endian_fread(&nProbePairNum, 4, 1, fpIn, little_endian_machine) != 1)
		{
			printf("Error: TileProbe_BpmapToBAR, cannot load probe pair number.\n");
			exit(EXIT_FAILURE);
		}
		vProbePairNum->pMatElement[ni] = (int)nProbePairNum;
		nTotalProbeNum += (int)nProbePairNum;

		/* for version 2.0 or later */
		if(fVersion > 1.5)
		{
			/* read group name */
			if(big_endian_fread(&nSeqNameLen, 4, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: TileProbe_BpmapToBAR, cannot load group name length.\n");
				exit(EXIT_FAILURE);
			}
			if(nSeqNameLen > 0)
			{
				pBARPos->vSeqData[ni]->pSeqGroupName = CreateString(nSeqNameLen);
				if(pBARPos->vSeqData[ni]->pSeqGroupName == NULL)
				{
					printf("Error: TileProbe_BpmapToBAR, cannot load group name.\n");
					exit(EXIT_FAILURE);
				}
				if(big_endian_fread(pBARPos->vSeqData[ni]->pSeqGroupName->m_pString, 1, nSeqNameLen, fpIn, little_endian_machine) != nSeqNameLen)
				{
					printf("Error: TileProbe_BpmapToBAR, cannot load group name.\n");
					exit(EXIT_FAILURE);
				}
				pBARPos->vSeqData[ni]->pSeqGroupName->m_pString[nSeqNameLen] = '\0';
			}

			/* read version name */
			if(big_endian_fread(&nSeqNameLen, 4, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: TileProbe_BpmapToBAR, cannot load version number length.\n");
				exit(EXIT_FAILURE);
			}
			if(nSeqNameLen > 0)
			{
				pBARPos->vSeqData[ni]->pSeqVersion = CreateString(nSeqNameLen);
				if(pBARPos->vSeqData[ni]->pSeqVersion == NULL)
				{
					printf("Error: TileProbe_BpmapToBAR, cannot load version number.\n");
					exit(EXIT_FAILURE);
				}
				if(big_endian_fread(pBARPos->vSeqData[ni]->pSeqVersion->m_pString, 1, nSeqNameLen, fpIn, little_endian_machine) != nSeqNameLen)
				{
					printf("Error: TileProbe_BpmapToBAR, cannot load version number.\n");
					exit(EXIT_FAILURE);
				}
				pBARPos->vSeqData[ni]->pSeqVersion->m_pString[nSeqNameLen] = '\0';
			}

			/* read paramters */
			if(big_endian_fread(&nParamNum, 4, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: TileProbe_BpmapToBAR, cannot load number of parameters.\n");
				exit(EXIT_FAILURE);
			}
			pBARPos->vSeqData[ni]->nParamNum = (int)nParamNum;
			if(nParamNum > 0)
			{
				pBARPos->vSeqData[ni]->vParamName = (struct tagString **)calloc(pBARPos->vSeqData[ni]->nParamNum, sizeof(struct tagString *));
				pBARPos->vSeqData[ni]->vParamValue = (struct tagString **)calloc(pBARPos->vSeqData[ni]->nParamNum, sizeof(struct tagString *));

				if( (pBARPos->vSeqData[ni]->vParamName == NULL) || (pBARPos->vSeqData[ni]->vParamValue == NULL) )
				{
					printf("Error: TileProbe_BpmapToBAR, cannot allocate memory for loading parameter name/value pairs\n");
					exit(EXIT_FAILURE);
				}
			}
			
			for(nj=0; nj<(int)nParamNum; nj++)
			{
				if(big_endian_fread(&nSeqNameLen, 4, 1, fpIn, little_endian_machine) != 1)
				{
					printf("Error: TileProbe_BpmapToBAR, cannot load parameter name length.\n");
					exit(EXIT_FAILURE);
				}
				if(nSeqNameLen > 0)
				{
					pBARPos->vSeqData[ni]->vParamName[nj] = CreateString(nSeqNameLen);
					if(pBARPos->vSeqData[ni]->vParamName[nj] == NULL)
					{
						printf("Error: TileProbe_BpmapToBAR, cannot load parameter name.\n");
						exit(EXIT_FAILURE);
					}
					if(big_endian_fread(pBARPos->vSeqData[ni]->vParamName[nj]->m_pString, 1, nSeqNameLen, fpIn, little_endian_machine) != nSeqNameLen)
					{
						printf("Error: TileProbe_BpmapToBAR, cannot load parameter name.\n");
						exit(EXIT_FAILURE);
					}
					pBARPos->vSeqData[ni]->vParamName[nj]->m_pString[nSeqNameLen] = '\0';
				}

				if(big_endian_fread(&nSeqNameLen, 4, 1, fpIn, little_endian_machine) != 1)
				{
					printf("Error: TileProbe_BpmapToBAR, cannot load parameter value length.\n");
					exit(EXIT_FAILURE);
				}
				if(nSeqNameLen > 0)
				{
					pBARPos->vSeqData[ni]->vParamValue[nj] = CreateString(nSeqNameLen);
					if(pBARPos->vSeqData[ni]->vParamValue[nj] == NULL)
					{
						printf("Error: TileProbe_BpmapToBAR, cannot load parameter value.\n");
						exit(EXIT_FAILURE);
					}
					if(big_endian_fread(pBARPos->vSeqData[ni]->vParamValue[nj]->m_pString, 1, nSeqNameLen, fpIn, little_endian_machine) != nSeqNameLen)
					{
						printf("Error: TileProbe_BpmapToBAR, cannot load parameter value.\n");
						exit(EXIT_FAILURE);
					}
					pBARPos->vSeqData[ni]->vParamValue[nj]->m_pString[nSeqNameLen] = '\0';
				}
			}
		}
	}

	/* load map */
	vSeqID = NULL;
	vSeqID = CreateIntMatrix(nSeqNum,1);
	if(vSeqID == NULL)
	{
		printf("Error: TileProbe_BpmapToBAR, cannot load sequence id!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<(int)nSeqNum; ni++)
	{
		nDivide1M = 0;
		if(strcmp(strGenomeGrp, "") != 0)
		{
			if(pBARPos->vSeqData[ni]->pSeqGroupName == NULL)
			{
				nDivide1M = 1;
			}
			else
			{
				if(strcmp(strGenomeGrp, pBARPos->vSeqData[ni]->pSeqGroupName->m_pString) != 0)
					nDivide1M = 1;
			}
		}

		/* load seq id */
		if(big_endian_fread(&nSeqID, 4, 1, fpIn, little_endian_machine) != 1)
		{
			printf("Error: TileProbe_BpmapToBAR, cannot load sequence ID.\n");
			exit(EXIT_FAILURE);
		}
		vSeqID->pMatElement[ni] = nSeqID;

		/* create initial memory, 1 for logPM value, 75/8 for A,C,G, 1 for nT, 4 for (A,C,G,T)^2, 1 for copy number */
		pBARPos->vSeqData[ni]->nDataNum = 0;
		pBARPos->vSeqData[ni]->nColNum = pBARPos->nColNum;
	
		pBARPos->vSeqData[ni]->vData = (struct DOUBLEMATRIX **)calloc(pBARPos->vSeqData[ni]->nColNum, sizeof(struct DOUBLEMATRIX *));
		if(pBARPos->vSeqData[ni]->vData == NULL)
		{
			printf("Error: TileProbe_BpmapToBAR, cannot allocate memory for loading bpmap coordinate data!\n");
			exit(EXIT_FAILURE);
		}
		nProbeNum = vProbePairNum->pMatElement[ni];
		if(nProbeNum > 0)
		{
			for(nj=0; nj<nInfoCol; nj++)
			{
				pBARPos->vSeqData[ni]->vData[nj] = CreateDoubleMatrix(1, nProbeNum);
				if(pBARPos->vSeqData[ni]->vData[nj] == NULL)
				{
					printf("Error: TileProbe_BpmapToBAR, cannot allocate memory for loading bpmap coordinate data!\n");
					exit(EXIT_FAILURE);
				}
			}
		}

		pUnitList = NULL;
		
		/* load probes */
		for(nj=0; nj<nProbeNum; nj++)
		{
			/* load new unit */
			pNewUnit = NULL;
			pNewUnit = AffyBpMapUnitCreate();
			AffyBpMapUnitLoad_v3m2(pNewUnit, fpIn, vProbeMappingType->pMatElement[ni], little_endian_machine);
			if(nDivide1M != 1)
			{
				pNewUnit->fMatchScore *= 1000000;
			}
			if( (pNewUnit->fMatchScore < 1.5) && (pNewUnit->fMatchScore > 0.5))
				pNewUnit->fMatchScore = 1.0;

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

							fprintf(fpOut, "%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%f\t%d\t%f\n", 
								pBARPos->vSeqData[ni]->pSeqName->m_pString, pCUnit->nPos, pCUnit->nDepthNum, 
								pCUnit->nRepeatNum, pCUnit->strProbeSeq, 
								pCUnit->nPMX, pCUnit->nPMY, pCUnit->nMMX, pCUnit->nMMY,
								pCUnit->fMatchScore, pCUnit->bStrand);

							if(pCUnit->nRepeatNum < 2)
							{
								nk = pBARPos->vSeqData[ni]->nDataNum;
								pBARPos->vSeqData[ni]->vData[0]->pMatElement[nk] = pCUnit->nPos;
								pBARPos->vSeqData[ni]->vData[1]->pMatElement[nk] = pCUnit->nPMX;
								pBARPos->vSeqData[ni]->vData[2]->pMatElement[nk] = pCUnit->nPMY;
								pBARPos->vSeqData[ni]->vData[3]->pMatElement[nk] = pCUnit->fMatchScore;
								pBARPos->vSeqData[ni]->nDataNum += 1;
							}

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

			fprintf(fpOut, "%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%f\t%d\n", 
						pBARPos->vSeqData[ni]->pSeqName->m_pString, pCUnit->nPos, pCUnit->nDepthNum, 
						pCUnit->nRepeatNum, pCUnit->strProbeSeq, 
						pCUnit->nPMX, pCUnit->nPMY, pCUnit->nMMX, pCUnit->nMMY,
						pCUnit->fMatchScore, pCUnit->bStrand);

			if(pCUnit->nRepeatNum < 2)
			{
				nk = pBARPos->vSeqData[ni]->nDataNum;
				pBARPos->vSeqData[ni]->vData[0]->pMatElement[nk] = pCUnit->nPos;
				pBARPos->vSeqData[ni]->vData[1]->pMatElement[nk] = pCUnit->nPMX;
				pBARPos->vSeqData[ni]->vData[2]->pMatElement[nk] = pCUnit->nPMY;
				pBARPos->vSeqData[ni]->vData[3]->pMatElement[nk] = pCUnit->fMatchScore;

				pBARPos->vSeqData[ni]->nDataNum += 1;
			}

			AffyBpMapUnitDestroy(pCUnit);
		}

		for(nj=0; nj<nInfoCol; nj++)
		{
			if(pBARPos->vSeqData[ni]->nDataNum > 0)
			{
				pResizeMat = NULL;
				pResizeMat = CreateDoubleMatrix(1, pBARPos->vSeqData[ni]->nDataNum);
				if(pResizeMat == NULL)
				{
					printf("Error: TileProbe_BpmapToBAR, cannot create an intermediate matrix to move data!\n");
					exit(EXIT_FAILURE);
				}
				memcpy(pResizeMat->pMatElement, pBARPos->vSeqData[ni]->vData[nj]->pMatElement, pBARPos->vSeqData[ni]->nDataNum*sizeof(double));
				DestroyDoubleMatrix(pBARPos->vSeqData[ni]->vData[nj]);
				pBARPos->vSeqData[ni]->vData[nj] = pResizeMat;
			}
			else
			{
				DestroyDoubleMatrix(pBARPos->vSeqData[ni]->vData[nj]);
				pBARPos->vSeqData[ni]->vData[nj] = NULL;
			}
		}

		/* update masked probe number */
		nMaskedProbeNum += pBARPos->vSeqData[ni]->nDataNum;
	}

	/* load tail if any */
	printf("  ");
	while(feof(fpIn) == 0)
	{
		fread(strFileType, sizeof(char), 1, fpIn);
		printf("%c", strFileType[0]);
	}
	printf("\n");

	/* clear memeory */
	DestroyIntMatrix(vProbeMappingType);
	if(fVersion > 2.5)
	{
		DestroyIntMatrix(vSequenceFileOffset);
	}
    DestroyIntMatrix(vProbePairNum);
	DestroyIntMatrix(vSeqID);

	/* close file */
	fclose(fpIn);
	fclose(fpOut);

	printf("  Probe number before masking = %d\n", nTotalProbeNum);
	printf("  Probe number after masking = %d\n", nMaskedProbeNum);

	/* return */
	return pBARPos;
}


/* ----------------------------------------------------------------------- */ 
/*  TileProbe_Resample_Main()                                              */
/*  Test robustness by resampling                                          */
/* ----------------------------------------------------------------------- */
int TileProbe_Resample_Main(char strParamFile[])
{
	/* define */
	FILE *fpIn;
	char strOutPath[MED_LINE_LENGTH];
	char *chp;
	char strFileName[MED_LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];

	/* train */
	int nStudyNum = 0;
	int nArrayNum = 0;
	double dNormLowerBound = -100000.0;
	int nTakeLog = 0;
	int nShrinkVar = 1;
	int nResampleStudy = 1;
	int nResampleNum = 10;
	int nW = 300;
	int nMaxGap = 300;
	int nMinProbe = 10;
	int nUseVar = 0;
	int nIPNum = 0;
	int nCTNum = 0;

	/* data */
	int *vStudyID = NULL;
	int *vGroupID = NULL;
	struct tagBARData **vTrain = NULL;
	struct tagBARData *vIP = NULL;
	struct tagBARData *vCT = NULL;
	struct tagBARData *pArray = NULL;
	int nRi,nGnum;
	int ni,nj,nk,nu,nw,nx,ny;
	double dRandNum;

	struct tagBARData **vTrainResamp = NULL;
	int *vGroupIDResamp = NULL;
	int nArrayNumResamp = 0;

	/* results */
	struct DOUBLEMATRIX *pM = NULL;
	struct DOUBLEMATRIX *pV = NULL;
	struct INTMATRIX *pFieldType = NULL;
	struct DOUBLEMATRIX **vDataVec = NULL;
	int nTotalProbeNum = 0;
	int nSeqNum = 0;
	int nProbeNum = 0;
	int nProbeNumC = 0;
	double dLog2 = log(2.0);

	struct DOUBLEMATRIX *pTPM,*pTPV,*pTP;
	struct DOUBLEMATRIX *pTPS;
	struct LONGMATRIX *pTPSI,*pTPMI,*pTPVI;

	/* counts */
	int nStep = 3000;
	int nBin = 50;
	struct INTMATRIX *pCountM = NULL;
	struct INTMATRIX *pCountV = NULL;
	struct INTMATRIX *pCountO = NULL;
	int nrmax;

	/* load parameters */
	fpIn = NULL;
	fpIn = fopen(strParamFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: TileProbe_Resample_Main, cannot open the input file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strstr(strLine, "Outpath") == strLine)
		{
			chp = strchr(strLine, '=');
			strcpy(strOutPath, chp+1);
			StrTrimLeft(strOutPath);
		}
		else if(strstr(strLine, "NormTrunc") == strLine)
		{
			chp = strchr(strLine, '=');
			dNormLowerBound = atof(chp+1);
		}
		else if(strstr(strLine, "NormLog") == strLine)
		{
			chp = strchr(strLine, '=');
			nTakeLog = atoi(chp+1);
		}
		else if(strstr(strLine, "ShrinkVar") == strLine)
		{
			chp = strchr(strLine, '=');
			nShrinkVar = atoi(chp+1);
		}
		else if(strstr(strLine, "StudyNum") == strLine)
		{
			chp = strchr(strLine, '=');
			nStudyNum = atoi(chp+1);
		}
		else if(strstr(strLine, "ArrayNum") == strLine)
		{
			chp = strchr(strLine, '=');
			nArrayNum = atoi(chp+1);

			if(nArrayNum == 0)
			{
				fclose(fpIn);
				return PROC_SUCCESS;
			}

			/* memory for group id */
			vGroupID = NULL;
			vGroupID = (int *)calloc(nArrayNum, sizeof(int));
			if(vGroupID == NULL)
			{
				printf("Error: TileProbe_Resample_Main, cannot allocate memory for group id!\n");
				exit(EXIT_FAILURE);
			}

			vStudyID = NULL;
			vStudyID = (int *)calloc(nArrayNum, sizeof(int));
			if(vStudyID == NULL)
			{
				printf("Error: TileProbe_Resample_Main, cannot allocate memory for study id!\n");
				exit(EXIT_FAILURE);
			}

			vTrain = NULL;
			vTrain = (struct tagBARData **)calloc(nArrayNum, sizeof(struct tagBARData *));
			if(vTrain == NULL)
			{
				printf("Error: TileProbe_Resample_Main, cannot allocate memory for training list!\n");
				exit(EXIT_FAILURE);
			}

			for(ni=0; ni<nArrayNum; ni++)
			{
				fgets(strLine, LONG_LINE_LENGTH, fpIn);
				sscanf(strLine, "%d %d %s", &nw, &nu, strFileName);
				vStudyID[ni] = nw;
				vGroupID[ni] = nu;

				printf("Loading %s ...\n", strFileName);
				vTrain[ni] = Affy_LoadBAR_Fast(strFileName);
				if(vTrain[ni] == NULL)
				{
					printf("Error: TileProbe_Resample_Main, cannot load *.bar file!\n");
					exit(EXIT_FAILURE);
				}

				nProbeNumC = 0;
				for(nj=0; nj<vTrain[ni]->nSeqNum; nj++)
				{
					nProbeNumC += vTrain[ni]->vSeqData[nj]->nDataNum;
				}
				
				if(ni == 0)
				{
					nSeqNum = vTrain[ni]->nSeqNum;
					nProbeNum = nProbeNumC;
					
				}
				else
				{
					if( (nSeqNum != vTrain[ni]->nSeqNum) || (nProbeNum != nProbeNumC) )
					{
						printf("Error: TileProbe_Resample_Main, file dimensions do not match!\n");
						exit(EXIT_FAILURE);
					}
				}

				if(nTakeLog == 1)
				{
					for(nj=0; nj<vTrain[ni]->nSeqNum; nj++)
					{
						for(nk=0; nk<vTrain[ni]->vSeqData[nj]->nDataNum; nk++)
						{
							/* truncate */
							if(vTrain[ni]->vSeqData[nj]->vData[1]->pMatElement[nk] < dNormLowerBound)
								vTrain[ni]->vSeqData[nj]->vData[1]->pMatElement[nk] = dNormLowerBound;

							/* log transformation */
							vTrain[ni]->vSeqData[nj]->vData[1]->pMatElement[nk] = log(vTrain[ni]->vSeqData[nj]->vData[1]->pMatElement[nk])/dLog2;
						}
					}
				}
				else
				{
					for(nj=0; nj<vTrain[ni]->nSeqNum; nj++)
					{
						for(nk=0; nk<vTrain[ni]->vSeqData[nj]->nDataNum; nk++)
						{
							/* truncate */
							if(vTrain[ni]->vSeqData[nj]->vData[1]->pMatElement[nk] < dNormLowerBound)
								vTrain[ni]->vSeqData[nj]->vData[1]->pMatElement[nk] = dNormLowerBound;
						}
					}
				}
			}
		}
		else if(strstr(strLine, "ResampleStudy") == strLine)
		{
			chp = strchr(strLine, '=');
			nResampleStudy = atoi(chp+1);
		}
		else if(strstr(strLine, "ResampleNum") == strLine)
		{
			chp = strchr(strLine, '=');
			nResampleNum = atoi(chp+1);
		}
		else if(strstr(strLine, "BandWidth") == strLine)
		{
			chp = strchr(strLine, '=');
			nW = atoi(chp+1);
		}
		else if(strstr(strLine, "MaxGap") == strLine)
		{
			chp = strchr(strLine, '=');
			nMaxGap = atoi(chp+1);
		}
		else if(strstr(strLine, "MinProbe") == strLine)
		{
			chp = strchr(strLine, '=');
			nMinProbe = atoi(chp+1);
		}
		else if(strstr(strLine, "Var") == strLine)
		{
			chp = strchr(strLine, '=');
			nUseVar = atoi(chp+1);
		}
		else if(strstr(strLine, "IP") == strLine)
		{
			chp = strchr(strLine, '=');
			nIPNum = atoi(chp+1);
			if(nIPNum <= 0)
			{
				printf("Error: TileProbe_Resample_Main, at least one IP sample need to be provided!\n");
				fclose(fpIn);
				exit(EXIT_FAILURE);
			}

			ni = 0;
			while((ni<nIPNum) && (fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL))
			{
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] == '\0')
					continue;
				if(strLine[0] == '#')
					continue;

				pArray = NULL;
				pArray = Affy_LoadBAR_Fast(strLine);
				if(pArray == NULL)
				{
					printf("Error: TileProbe_Resample_Main, cannot load raw data!\n");
					exit(EXIT_FAILURE);
				}
				if(nTakeLog == 1)
				{
					for(nj=0; nj<pArray->nSeqNum; nj++)
					{
						for(nk=0; nk<pArray->vSeqData[nj]->nDataNum; nk++)
						{
							/* truncate */
							if(pArray->vSeqData[nj]->vData[1]->pMatElement[nk] < dNormLowerBound)
								pArray->vSeqData[nj]->vData[1]->pMatElement[nk] = dNormLowerBound;

							/* log transformation */
							pArray->vSeqData[nj]->vData[1]->pMatElement[nk] = log(pArray->vSeqData[nj]->vData[1]->pMatElement[nk])/dLog2;
						}
					}
				}
				else
				{
					for(nj=0; nj<pArray->nSeqNum; nj++)
					{
						for(nk=0; nk<pArray->vSeqData[nj]->nDataNum; nk++)
						{
							/* truncate */
							if(pArray->vSeqData[nj]->vData[1]->pMatElement[nk] < dNormLowerBound)
								pArray->vSeqData[nj]->vData[1]->pMatElement[nk] = dNormLowerBound;
						}
					}
				}

				/* if the first array, prepare enough space for nIPNum arrays + 
					1 positions + 1 group trimed mean + 1 group variance + 1 test statistics */ 
				if(ni == 0)
				{
					vIP = pArray;
					vIP->nColNum = nIPNum+4;
					pFieldType = vIP->pFieldType;
					vIP->pFieldType = NULL;
					vIP->pFieldType = CreateIntMatrix(1,vIP->nColNum);
					if(vIP->pFieldType == NULL)
					{
						printf("Error: TileProbe_Resample_Main, cannot create memory for field types!\n");
						exit(EXIT_FAILURE);
					}
					for(nj=0; nj<vIP->nColNum; nj++)
						vIP->pFieldType->pMatElement[nj] = 1;
					vIP->pFieldType->pMatElement[0] = pFieldType->pMatElement[0];
					vIP->pFieldType->pMatElement[1] = pFieldType->pMatElement[1];
					DestroyIntMatrix(pFieldType);
					
					nTotalProbeNum = 0;
					for(nj=0; nj<vIP->nSeqNum; nj++)
					{
						nTotalProbeNum += vIP->vSeqData[nj]->nDataNum;
						vIP->vSeqData[nj]->nColNum = vIP->nColNum;
						vDataVec = vIP->vSeqData[nj]->vData;
						vIP->vSeqData[nj]->vData = NULL;
						vIP->vSeqData[nj]->vData = (struct DOUBLEMATRIX **)calloc(vIP->vSeqData[nj]->nColNum, sizeof(struct DOUBLEMATRIX *));
						if(vIP->vSeqData[nj]->vData == NULL)
						{
							printf("Error: TileProbe_Resample_Main, cannot create memory for tracking intensity data!\n");
							exit(EXIT_FAILURE);
						}
						vIP->vSeqData[nj]->vData[0] = vDataVec[0];
						vIP->vSeqData[nj]->vData[1] = vDataVec[1];
						free(vDataVec);

						for(nk=nIPNum+1; nk<vIP->nColNum; nk++)
						{
							vIP->vSeqData[nj]->vData[nk] = CreateDoubleMatrix(1, vIP->vSeqData[nj]->nDataNum);
							if((vIP->vSeqData[nj]->vData[nk] == NULL) && (vIP->vSeqData[nj]->nDataNum > 0))
							{
								printf("Error: TileProbe_Resample_Main, cannot create memory for data processing!\n");
								exit(EXIT_FAILURE);
							}
						}
					}

					if(nTotalProbeNum != nProbeNum)
					{
						printf("Error: TileProbe_Resample_Main, probe number do not match!\n");
						exit(EXIT_FAILURE);
					}
				}
				/* if not the first array, transfer to the first array */
				else
				{
					if(pArray->nSeqNum != vIP->nSeqNum)
					{
						printf("Error: TileProbe_Resample_Main, array types do not match!\n");
						exit(EXIT_FAILURE);
					}

					for(nj=0; nj<vIP->nSeqNum; nj++)
					{
						if(vIP->vSeqData[nj]->nDataNum != pArray->vSeqData[nj]->nDataNum)
						{
							printf("Error: TileProbe_Resample_Main, array types do not match!\n");
							exit(EXIT_FAILURE);
						}
						vIP->vSeqData[nj]->vData[ni+1] = pArray->vSeqData[nj]->vData[1];
						pArray->vSeqData[nj]->vData[1] = NULL;
					}

					vIP->pFieldType->pMatElement[ni+1] = pArray->pFieldType->pMatElement[1];

					Affy_BARData_Destroy(&pArray);
				}


				ni++;
			}

			
			if(ni!=nIPNum)
			{
				printf("Error: TileProbe_Resample_Main, IP sample number inconsistent!\n");
				exit(EXIT_FAILURE);
			}
		}

		else if(strstr(strLine, "CT") == strLine)
		{
			chp = strchr(strLine, '=');
			nCTNum = atoi(chp+1);
		
			ni = 0;
			while((ni<nCTNum) && (fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL))
			{
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] == '\0')
					continue;
				if(strLine[0] == '#')
					continue;

				pArray = NULL;
				pArray = Affy_LoadBAR_Fast(strLine);
				if(pArray == NULL)
				{
					printf("Error: TileProbe_Resample_Main, cannot load raw data!\n");
					exit(EXIT_FAILURE);
				}
				if(nTakeLog == 1)
				{
					for(nj=0; nj<pArray->nSeqNum; nj++)
					{
						for(nk=0; nk<pArray->vSeqData[nj]->nDataNum; nk++)
						{
							/* truncate */
							if(pArray->vSeqData[nj]->vData[1]->pMatElement[nk] < dNormLowerBound)
								pArray->vSeqData[nj]->vData[1]->pMatElement[nk] = dNormLowerBound;

							/* log transformation */
							pArray->vSeqData[nj]->vData[1]->pMatElement[nk] = log(pArray->vSeqData[nj]->vData[1]->pMatElement[nk])/dLog2;
						}
					}
				}
				else
				{
					for(nj=0; nj<pArray->nSeqNum; nj++)
					{
						for(nk=0; nk<pArray->vSeqData[nj]->nDataNum; nk++)
						{
							/* truncate */
							if(pArray->vSeqData[nj]->vData[1]->pMatElement[nk] < dNormLowerBound)
								pArray->vSeqData[nj]->vData[1]->pMatElement[nk] = dNormLowerBound;
						}
					}
				}

				/* if the first array, prepare enough space for nIPNum arrays + 
					1 positions + 1 group trimed mean + 1 group variance + 1 test statistics */ 
				if(ni == 0)
				{
					vCT = pArray;
					vCT->nColNum = nCTNum+4;
					pFieldType = vCT->pFieldType;
					vCT->pFieldType = NULL;
					vCT->pFieldType = CreateIntMatrix(1,vCT->nColNum);
					if(vCT->pFieldType == NULL)
					{
						printf("Error: TileProbe_Peak_Main, cannot create memory for field types!\n");
						exit(EXIT_FAILURE);
					}
					for(nj=0; nj<vCT->nColNum; nj++)
						vCT->pFieldType->pMatElement[nj] = 1;
					vCT->pFieldType->pMatElement[0] = pFieldType->pMatElement[0];
					vCT->pFieldType->pMatElement[1] = pFieldType->pMatElement[1];
					DestroyIntMatrix(pFieldType);
					
					for(nj=0; nj<vCT->nSeqNum; nj++)
					{
						vCT->vSeqData[nj]->nColNum = vCT->nColNum;
						vDataVec = vCT->vSeqData[nj]->vData;
						vCT->vSeqData[nj]->vData = NULL;
						vCT->vSeqData[nj]->vData = (struct DOUBLEMATRIX **)calloc(vCT->vSeqData[nj]->nColNum, sizeof(struct DOUBLEMATRIX *));
						if(vCT->vSeqData[nj]->vData == NULL)
						{
							printf("Error: TileProbe_Resample_Main, cannot create memory for tracking intensity data!\n");
							exit(EXIT_FAILURE);
						}
						vCT->vSeqData[nj]->vData[0] = vDataVec[0];
						vCT->vSeqData[nj]->vData[1] = vDataVec[1];
						free(vDataVec);

						for(nk=nCTNum+1; nk<vCT->nColNum; nk++)
						{
							vCT->vSeqData[nj]->vData[nk] = CreateDoubleMatrix(1, vCT->vSeqData[nj]->nDataNum);
							if( (vCT->vSeqData[nj]->vData[nk] == NULL) && (vCT->vSeqData[nj]->nDataNum > 0) )
							{
								printf("Error: TileProbe_Peak_Main, cannot create memory for data processing!\n");
								exit(EXIT_FAILURE);
							}
						}
					}
				}
				/* if not the first array, transfer to the first array */
				else
				{
					if(pArray->nSeqNum != vCT->nSeqNum)
					{
						printf("Error: TileProbe_Resample_Main, array types do not match!\n");
						exit(EXIT_FAILURE);
					}

					for(nj=0; nj<vCT->nSeqNum; nj++)
					{
						if(vCT->vSeqData[nj]->nDataNum != pArray->vSeqData[nj]->nDataNum)
						{
							printf("Error: TileProbe_Resample_Main, array types do not match!\n");
							exit(EXIT_FAILURE);
						}
						vCT->vSeqData[nj]->vData[ni+1] = pArray->vSeqData[nj]->vData[1];
						pArray->vSeqData[nj]->vData[1] = NULL;
					}

					vCT->pFieldType->pMatElement[ni+1] = pArray->pFieldType->pMatElement[1];

					Affy_BARData_Destroy(&pArray);
				}


				ni++;
			}

			
			if(ni!=nCTNum)
			{
				printf("Error: TileProbe_Resample_Main, CT sample number inconsistent!\n");
				exit(EXIT_FAILURE);
			}
		}

	}
	
	fclose(fpIn);

	/* compute original model */
	/* need a function, return bar file */
	pM = NULL;
	pM = CreateDoubleMatrix(1, nProbeNum);
	if(pM == NULL)
	{
		printf("Error: TileProbe_Resample_Main, cannot allocate memory for median!\n");
		exit(EXIT_FAILURE);
	}
	pV = NULL;
	pV = CreateDoubleMatrix(1, nProbeNum);
	if(pV == NULL)
	{
		printf("Error: TileProbe_Resample_Main, cannot allocate memory for variance!\n");
		exit(EXIT_FAILURE);
	}

	TileProbe_Resample_Model(nArrayNum, nProbeNum, nShrinkVar, vGroupID, vTrain, pM, pV);

	/* compute TPM, apply MAT */
	/* need a function for MAT */
	pTPM = NULL;
	pTPM = TileProbe_Resample_MATPeak(nIPNum, nCTNum, nProbeNum, vIP, vCT, nW, nMaxGap, nMinProbe, nUseVar, pM, pV, 0);
	if(pTPM == NULL)
	{
		printf("Error: TileProbe_Resample_Main, cannot compute TPM!\n");
		exit(EXIT_FAILURE);
	}
	
	/* compute TPV, apply MAT */
	pTPV = NULL;
	pTPV = TileProbe_Resample_MATPeak(nIPNum, nCTNum, nProbeNum, vIP, vCT, nW, nMaxGap, nMinProbe, nUseVar, pM, pV, 1);
	if(pTPV == NULL)
	{
		printf("Error: TileProbe_Resample_Main, cannot compute TPV!\n");
		exit(EXIT_FAILURE);
	}

	/* sort and get rank */
	pCountO = NULL;
	pCountO = CreateIntMatrix(nBin, 1);
	if(pCountO == NULL)
	{
		printf("Error: TileProbe_Resample_Main, cannot create count!\n");
		exit(EXIT_FAILURE);
	}

	pTPMI = NULL;
	pTPMI = CreateLongMatrix(1, nProbeNum);
	if(pTPMI == NULL)
	{
		printf("Error: TileProbe_Resample_Main, cannot create count!\n");
		exit(EXIT_FAILURE);
	}

	pTPVI = NULL;
	pTPVI = CreateLongMatrix(1, nProbeNum);
	if(pTPVI == NULL)
	{
		printf("Error: TileProbe_Resample_Main, cannot create count!\n");
		exit(EXIT_FAILURE);
	}

	pTPS = NULL;
	pTPSI = NULL;
	DMSORTMERGEA_0(pTPM, &pTPS, &pTPSI);
	for(ni=0; ni<nProbeNum; ni++)
	{
		nk = pTPSI->pMatElement[ni];
		pTPMI->pMatElement[nk] = nProbeNum-ni-1;
	}
	DestroyDoubleMatrix(pTPS);
	DestroyLongMatrix(pTPSI);

	pTPS = NULL;
	pTPSI = NULL;
	DMSORTMERGEA_0(pTPV, &pTPS, &pTPSI);
	for(ni=0; ni<nProbeNum; ni++)
	{
		nk = pTPSI->pMatElement[ni];
		pTPVI->pMatElement[nk] = nProbeNum-ni-1;

		nrmax = pTPVI->pMatElement[nk];
		if(pTPMI->pMatElement[nk] > nrmax)
			nrmax = pTPMI->pMatElement[nk];
		nk = nrmax/nStep;
		if(nk < nBin)
			pCountO->pMatElement[nk] += 1;
	}
	DestroyDoubleMatrix(pTPS);
	DestroyLongMatrix(pTPSI);

	/* resample */
	if(nResampleStudy == 1)
	{
		pCountM = NULL;
		pCountM = CreateIntMatrix(nBin, nStudyNum);
		if(pCountM == NULL)
		{
			printf("Error: TileProbe_Resample_Main, cannot create count!\n");
			exit(EXIT_FAILURE);
		}
		pCountV = NULL;
		pCountV = CreateIntMatrix(nBin, nStudyNum);
		if(pCountV == NULL)
		{
			printf("Error: TileProbe_Resample_Main, cannot create count!\n");
			exit(EXIT_FAILURE);
		}


		/* resample study */
		for(nRi=1; nRi<=nStudyNum; nRi++)
		{
			/* compute model after leaving one study out */
			nArrayNumResamp = 0;
			for(ni=0; ni<nArrayNum; ni++)
			{
				if(vStudyID[ni] != nRi)
					nArrayNumResamp += 1;
			}
			
			if(nArrayNumResamp == 0)
			{
				printf("Warning: resample contains 0 arrays!\n");
				continue;
			}

			vGroupIDResamp = NULL;
			vGroupIDResamp = (int *)calloc(nArrayNumResamp, sizeof(int));
			if(vGroupIDResamp == NULL)
			{
				printf("Error: TileProbe_Resample_Main, cannot allocate memory for group id!\n");
				exit(EXIT_FAILURE);
			}

			vTrainResamp = NULL;
			vTrainResamp = (struct tagBARData **)calloc(nArrayNumResamp, sizeof(struct tagBARData *));
			if(vTrainResamp == NULL)
			{
				printf("Error: TileProbe_Resample_Main, cannot allocate memory for training list!\n");
				exit(EXIT_FAILURE);
			}

			nj = 0;
			for(ni=0; ni<nArrayNum; ni++)
			{
				if(vStudyID[ni] != nRi)
				{
					vGroupIDResamp[nj] = vGroupID[ni];
					vTrainResamp[nj] = vTrain[ni];
					nj++;
				}
			}
			if(nj != nArrayNumResamp)
			{
				printf("Error: TileProbe_Resample_Main, resample incorrect!\n");
				exit(EXIT_FAILURE);
			}

			DestroyDoubleMatrix(pM);
			DestroyDoubleMatrix(pV);
			pM = NULL;
			pV = NULL;
			pM = CreateDoubleMatrix(1, nProbeNum);
			pV = CreateDoubleMatrix(1, nProbeNum);
			if((pM == NULL) || (pV == NULL))
			{
				printf("Error: TileProbe_Resample_Main, cannot allocate memory for probe model!\n");
				exit(EXIT_FAILURE);
			}
			TileProbe_Resample_Model(nArrayNumResamp, nProbeNum, nShrinkVar, vGroupIDResamp, vTrainResamp, pM, pV);

			free(vGroupIDResamp);
			free(vTrainResamp);

			/* compute TPV, apply MAT */
			pTP = NULL;
			pTP = TileProbe_Resample_MATPeak(nIPNum, nCTNum, nProbeNum, vIP, vCT, nW, nMaxGap, nMinProbe, nUseVar, pM, pV, 1);
			if(pTP == NULL)
			{
				printf("Error: TileProbe_Resample_Main, cannot compute TPV!\n");
				exit(EXIT_FAILURE);
			}

			/* count overlap */
			pTPS = NULL;
			pTPSI = NULL;
			DMSORTMERGEA_0(pTP, &pTPS, &pTPSI);
			for(ni=0; ni<nProbeNum; ni++)
			{
				nk = pTPSI->pMatElement[ni];
				nu = nProbeNum-ni-1;

				nrmax = nu;
				if(pTPMI->pMatElement[nk] > nrmax)
					nrmax = pTPMI->pMatElement[nk];
				nx = nrmax/nStep;
				if(nx < nBin)
				{
					ny = IMGETAT(pCountM, nx, nRi-1)+1;
					IMSETAT(pCountM, nx, nRi-1, ny);
				}

				nrmax = nu;
				if(pTPVI->pMatElement[nk] > nrmax)
					nrmax = pTPVI->pMatElement[nk];
				nx = nrmax/nStep;
				if(nx < nBin)
				{
					ny = IMGETAT(pCountV, nx, nRi-1)+1;
					IMSETAT(pCountV, nx, nRi-1, ny);
				}
			}

			DestroyDoubleMatrix(pTPS);
			DestroyLongMatrix(pTPSI);

			/* destroy */
			DestroyDoubleMatrix(pTP);
		}
	}
	else
	{
		pCountM = NULL;
		pCountM = CreateIntMatrix(nBin, nResampleNum);
		if(pCountM == NULL)
		{
			printf("Error: TileProbe_Resample_Main, cannot create count!\n");
			exit(EXIT_FAILURE);
		}
		pCountV = NULL;
		pCountV = CreateIntMatrix(nBin, nResampleNum);
		if(pCountV == NULL)
		{
			printf("Error: TileProbe_Resample_Main, cannot create count!\n");
			exit(EXIT_FAILURE);
		}

		/* keep median, compute residual */
		TileProbe_Resample_Residual(nArrayNum, nProbeNum, vGroupID, vTrain);

		for(nRi=0; nRi<nResampleNum; nRi++)
		{
			/* resample array id */
			vTrainResamp = NULL;
			vTrainResamp = (struct tagBARData **)calloc(nArrayNum, sizeof(struct tagBARData *));
			if(vTrainResamp == NULL)
			{
				printf("Error: TileProbe_Resample_Main, cannot allocate memory for training list!\n");
				exit(EXIT_FAILURE);
			}

			ni = 0;
			nj = 0;
			while(ni < nArrayNum)
			{
				for(nj=ni; nj<nArrayNum; nj++)
				{
					if(vGroupID[nj] != vGroupID[ni])
						break;
				}
				nGnum = nj-ni;

				for(nk=ni; nk<nj; nk++)
				{
					dRandNum = rand_u();
					nx = (int)(dRandNum*nGnum);
					if(nx == nGnum)
						nx = nGnum-1;
					vTrainResamp[nk] = vTrain[ni+nx];
				}

				ni = nj;
			}
			if(nj != nArrayNum)
			{
				printf("Error: TileProbe_Resample_Main, resample incorrect!\n");
				exit(EXIT_FAILURE);
			}

			/* using residual to compute variance */
			DestroyDoubleMatrix(pV);
			pV = NULL;
			pV = CreateDoubleMatrix(1, nProbeNum);
			TileProbe_Resample_VarModel(nArrayNum, nProbeNum, nShrinkVar, vGroupID, vTrainResamp, pV);
			free(vTrainResamp);

			/* apply MAT */
			pTP = NULL;
			pTP = TileProbe_Resample_MATPeak(nIPNum, nCTNum, nProbeNum, vIP, vCT, nW, nMaxGap, nMinProbe, nUseVar, pM, pV, 1);
			if(pTP == NULL)
			{
				printf("Error: TileProbe_Resample_Main, cannot compute TPV!\n");
				exit(EXIT_FAILURE);
			}

			/* count overlap */
			pTPS = NULL;
			pTPSI = NULL;
			DMSORTMERGEA_0(pTP, &pTPS, &pTPSI);
			for(ni=0; ni<nProbeNum; ni++)
			{
				nk = pTPSI->pMatElement[ni];
				nu = nProbeNum-ni-1;

				nrmax = nu;
				if(pTPMI->pMatElement[nk] > nrmax)
					nrmax = pTPMI->pMatElement[nk];
				nx = nrmax/nStep;
				if(nx < nBin)
				{
					ny = IMGETAT(pCountM, nx, nRi)+1;
					IMSETAT(pCountM, nx, nRi, ny);
				}

				nrmax = nu;
				if(pTPVI->pMatElement[nk] > nrmax)
					nrmax = pTPVI->pMatElement[nk];
				nx = nrmax/nStep;
				if(nx < nBin)
				{
					ny = IMGETAT(pCountV, nx, nRi)+1;
					IMSETAT(pCountV, nx, nRi, ny);
				}
			}

			DestroyDoubleMatrix(pTPS);
			DestroyLongMatrix(pTPSI);

			/* destroy */
			DestroyDoubleMatrix(pTP);
		}
	}

	/* output */
	for(ni=1; ni<pCountM->nHeight; ni++)
	{
		nk = pCountO->pMatElement[ni]+pCountO->pMatElement[ni-1];
		pCountO->pMatElement[ni] = nk;

		for(nj=0; nj<pCountM->nWidth; nj++)
		{
			nk = IMGETAT(pCountM, ni-1, nj)+IMGETAT(pCountM, ni, nj);
			IMSETAT(pCountM, ni, nj, nk);
			nk = IMGETAT(pCountV, ni-1, nj)+IMGETAT(pCountV, ni, nj);
			IMSETAT(pCountV, ni, nj, nk);
		}
	}

	sprintf(strFileName, "%s_M.txt", strOutPath);
	IMSAVE(pCountM, strFileName);
	sprintf(strFileName, "%s_V.txt", strOutPath);
	IMSAVE(pCountV, strFileName);
	sprintf(strFileName, "%s_O.txt", strOutPath);
	IMSAVE(pCountO, strFileName);


	/* release memory */
	free(vGroupID);
	free(vStudyID);
	for(ni=0; ni<nArrayNum; ni++)
	{
		Affy_BARData_Destroy(vTrain+ni);
	}
	free(vTrain);
	Affy_BARData_Destroy(&vIP);
	Affy_BARData_Destroy(&vCT);
	DestroyDoubleMatrix(pM);
	DestroyDoubleMatrix(pV);
	DestroyDoubleMatrix(pTPM);
	DestroyDoubleMatrix(pTPV);
	DestroyLongMatrix(pTPMI);
	DestroyLongMatrix(pTPVI);
	DestroyIntMatrix(pCountM);
	DestroyIntMatrix(pCountV);
	DestroyIntMatrix(pCountO);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_Resample_Model()                                             */
/*  Build probe model                                                      */
/* ----------------------------------------------------------------------- */
int TileProbe_Resample_Model(int nArrayNum, int nProbeNum, int nShrink,
							 int *vGroupID, struct tagBARData **vBARData, 
							 struct DOUBLEMATRIX *pM, struct DOUBLEMATRIX *pV)
{
	/* define */
	struct DOUBLEMATRIX *pProbeData;
	struct DOUBLEMATRIX *pProbeDataSort;
	struct LONGMATRIX *pProbeDataSortId;
	
	/* others */
	int ni,nj,nk,nu,nx,ny,nz;
	int nQ50;
	
	/* for IQR shrinking */
	double dB;
	double dM,dV;
	int nGN,nDf;
	double dQmin,dQmax;

	
	pProbeData = NULL;
	pProbeData = CreateDoubleMatrix(1,nArrayNum);
	if(pProbeData == NULL)
	{
		printf("Error: TileProbe_Resample_Model, cannot allocate memory for modeling probe effect!\n");
		exit(EXIT_FAILURE);
	}

	nQ50 = (int)(nArrayNum*0.5);
	if( (nArrayNum*0.5)-nQ50 > 1e-6 )
		nQ50++;
	nQ50--;

	if(nQ50 < 0)
		nQ50 = 0;
	if(nQ50 >= nArrayNum)
		nQ50 = nArrayNum-1;

	nk = 0;
	nDf = 0;
	for(ni=0; ni<vBARData[0]->nSeqNum; ni++)
	{
		for(nj=0; nj<vBARData[0]->vSeqData[ni]->nDataNum; nj++)
		{
			if(nk%100000 == 0)
			{
				printf("%d ...\n", nk);
			}

			pProbeDataSort = NULL;
			pProbeDataSortId = NULL;

			for(nu=0; nu<nArrayNum; nu++)
			{
				pProbeData->pMatElement[nu] = vBARData[nu]->vSeqData[ni]->vData[1]->pMatElement[nj];
			}

			nx = 0;
			dV = 0.0;
			nDf = 0;
			while(nx < nArrayNum)
			{
				dM = 0.0;
				
				nGN = 0;
				for(ny=nx; ny<nArrayNum; ny++)
				{
					if(vGroupID[ny] != vGroupID[nx])
						break;
					dM += pProbeData->pMatElement[ny];
					nGN++;
				}
				dM = dM/nGN;
				for(nz=nx; nz<ny; nz++)
				{
					dV += (pProbeData->pMatElement[nz]-dM)*(pProbeData->pMatElement[nz]-dM);
				}
				nDf = nDf+nGN-1;
				nx = ny;
			}
			dV /= nDf;
			pV->pMatElement[nk] = dV;

			DMSORTMERGEA_0(pProbeData, &pProbeDataSort, &pProbeDataSortId);
			pM->pMatElement[nk] = pProbeDataSort->pMatElement[nQ50];
			

			DestroyDoubleMatrix(pProbeDataSort);
			DestroyLongMatrix(pProbeDataSortId);

			nk++;
		}
	}

	if(nk != nProbeNum)
	{
		printf("Error: TileProbe_Resample_Model, probe number do not match!\n");
		exit(EXIT_FAILURE);
	}
	if(nDf == 0)
	{
		printf("Error: TileProbe_Resample_Model, no degree of freedom to compute variance!\n");
		exit(EXIT_FAILURE);
	}

	/* Variance shrinkage */
	if(nShrink != 0)
	{
		dM = 0.0;
		dV = 0.0;
		dQmin = 1e10;
		dQmax = -1e10;

		for(nj=0; nj<nProbeNum; nj++)
		{
			if(pV->pMatElement[nj] < dQmin)
				dQmin = pV->pMatElement[nj];
			if(pV->pMatElement[nj] > dQmax)
				dQmax = pV->pMatElement[nj];

			dM += pV->pMatElement[nj];
			dV += pV->pMatElement[nj]*pV->pMatElement[nj];
		}

		dM /= nProbeNum;
		dV = dV-nProbeNum*dM*dM;
		dB = 2.0*(nProbeNum-1)/(2.0+nDf)/nProbeNum+2.0*dM*dM*(nProbeNum-1)/(2.0+nDf)/dV;
		
		if(dB>1.0)
			dB=1.0;
		else if(dB<0.0)
			dB=0.0;

		if(nShrink == -1)
			dB = 0.0;

		for(nj=0; nj<nProbeNum; nj++)
		{
			pV->pMatElement[nj] = sqrt(pV->pMatElement[nj]*(1.0-dB)+dM*dB);
		}
		
		printf("min Var = %f; max Var = %f\n", dQmin, dQmax);
		printf("B = %f\n", dB);
	}
	else
	{
		for(nj=0; nj<nProbeNum; nj++)
		{
			pV->pMatElement[nj] = sqrt(pV->pMatElement[nj]);
		}
	}

	/* release memory */
	DestroyDoubleMatrix(pProbeData);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_Resample_MATPeak()                                           */
/*  Find peaks.                                                            */
/* ----------------------------------------------------------------------- */ 
struct DOUBLEMATRIX *TileProbe_Resample_MATPeak(int nIPNum, int nCTNum, int nProbeNum,
				   struct tagBARData *vIP, struct tagBARData *vCT,
				   int nBandWidth, int nMaxGap, int nMinProbe, int nUseVar,
				   struct DOUBLEMATRIX *pM, struct DOUBLEMATRIX *pV, int nDivideV)
{
	/* define */
	int ni,nj,nk,ny,nz,nu;
	int nCurrentPos;
	int nP1,nP2;
	int nQ1,nQ2;
	double dPrc = 0.1;
	int nCellNum;
	struct DOUBLEMATRIX *pWD;
	struct DOUBLEMATRIX *pWDSort;
	struct LONGMATRIX *pWDId;
	double dTemp,dTemp2;
	int nIPMATCol,nCTMATCol;
	/* 
	int nrank;
	struct DOUBLEMATRIX *pPosRegion;
	struct DOUBLEMATRIX *pNegRegion;
	struct DOUBLEMATRIX *pFDR;
	struct INTMATRIX *pCol;
	char strOutFile[MED_LINE_LENGTH];
	FILE *fpOut; 
	*/
	struct DOUBLEMATRIX *pMATScore;
	int nTCol;

	/* init */
	nIPMATCol = 1+nIPNum;
	nCTMATCol = 1+nCTNum;

	for(ni=0; ni<vIP->nSeqNum; ni++)
	{
		for(nj=0; nj<vIP->vSeqData[ni]->nDataNum; nj++)
		{
			vIP->vSeqData[ni]->vData[nIPMATCol]->pMatElement[nj] = 0.0;
			vIP->vSeqData[ni]->vData[nIPMATCol+2]->pMatElement[nj] = 0.0;
		}
	}
	
	if(vCT != NULL)
	{
		for(ni=0; ni<vCT->nSeqNum; ni++)
		{
			for(nj=0; nj<vCT->vSeqData[ni]->nDataNum; nj++)
			{
				vCT->vSeqData[ni]->vData[nCTMATCol]->pMatElement[nj] = 0.0;
				vCT->vSeqData[ni]->vData[nCTMATCol+1]->pMatElement[nj] = 0.0;
			}
		}
	}

	/* scan the genome to compute statistics */
	nu = 0;
	for(ni=0; ni<vIP->nSeqNum; ni++)
	{
		nj = 0;
		
		nP1 = 0;
		nP2 = 0;

		for(nj=0; nj<vIP->vSeqData[ni]->nDataNum; nj++)
		{
			nCurrentPos = (int)(vIP->vSeqData[ni]->vData[0]->pMatElement[nj]);
			for(; nP1<nj; nP1++)
			{
				if( (int)(vIP->vSeqData[ni]->vData[0]->pMatElement[nP1]) >= (nCurrentPos-nBandWidth) )
					break;
			}

			for(; nP2<vIP->vSeqData[ni]->nDataNum; nP2++)
			{
				if( (int)(vIP->vSeqData[ni]->vData[0]->pMatElement[nP2]) > (nCurrentPos+nBandWidth) )
				{
					nP2--;
					break;
				}
			}
			if(nP2 >= vIP->vSeqData[ni]->nDataNum)
				nP2 = vIP->vSeqData[ni]->nDataNum-1;

			nCellNum = (nP2-nP1+1);
			if(nCellNum < nMinProbe)
			{
				continue;
			}
			
			nCellNum = (nP2-nP1+1)*nIPNum;
			nQ1 = (int)(nCellNum*dPrc);
			nQ2 = (int)(nCellNum*(1-dPrc))-1;

			
			pWD = NULL;
			pWDSort = NULL;
			pWDId = NULL;
			pWD = CreateDoubleMatrix(1, nCellNum);
			if(pWD == NULL)
			{
				printf("Error: TileProbe_Peak, cannot create matrix for sliding window!\n");
				exit(EXIT_FAILURE);
			}
			ny = 0;
			for(nz=0; nz<nIPNum; nz++)
			{
				for(nk=nP1; nk<=nP2; nk++)
				{
					pWD->pMatElement[ny] = vIP->vSeqData[ni]->vData[1+nz]->pMatElement[nk]-pM->pMatElement[nu+nk];
					if(nDivideV == 1)
					{
						pWD->pMatElement[ny] /= (pV->pMatElement[nu+nk]+1e-6);
					}
					ny++;
				}
			}
			DMSORTMERGEA_0(pWD, &pWDSort, &pWDId);

			dTemp = 0.0;
			for(nz=nQ1; nz<=nQ2; nz++)
			{
				dTemp += pWDSort->pMatElement[nz];
			}
			dTemp /= (nQ2-nQ1+1);
			vIP->vSeqData[ni]->vData[nIPMATCol]->pMatElement[nj] = dTemp;
			vIP->vSeqData[ni]->vData[nIPMATCol+2]->pMatElement[nj] = dTemp*sqrt((double)(nQ2-nQ1+1));

			DestroyDoubleMatrix(pWD);
			DestroyDoubleMatrix(pWDSort);
			DestroyLongMatrix(pWDId);

			if(vCT != NULL)
			{
				if( ((int)(vIP->vSeqData[ni]->vData[0]->pMatElement[nP1]) != (int)(vCT->vSeqData[ni]->vData[0]->pMatElement[nP1])) || 
					((int)(vIP->vSeqData[ni]->vData[0]->pMatElement[nP2]) != (int)(vCT->vSeqData[ni]->vData[0]->pMatElement[nP2])) )
				{
					printf("Error: TileProbe_Peak, coordinates in input bar files do not match!\n");
					exit(EXIT_FAILURE);
				}

				nCellNum = (nP2-nP1+1)*nCTNum;
				nQ1 = (int)(nCellNum*dPrc);
				nQ2 = (int)(nCellNum*(1-dPrc))-1;

				pWD = NULL;
				pWDSort = NULL;
				pWDId = NULL;
				pWD = CreateDoubleMatrix(1, nCellNum);
				if(pWD == NULL)
				{
					printf("Error: TileProbe_Peak, cannot create matrix for sliding window!\n");
					exit(EXIT_FAILURE);
				}
				ny = 0;
				for(nz=0; nz<nCTNum; nz++)
				{
					for(nk=nP1; nk<=nP2; nk++)
					{
						pWD->pMatElement[ny] = vCT->vSeqData[ni]->vData[1+nz]->pMatElement[nk]-pM->pMatElement[nu+nk];
						if(nDivideV == 1)
						{
							pWD->pMatElement[ny] /= (pV->pMatElement[nu+nk]+1e-6);
						}
						ny++;
					}
				}
				DMSORTMERGEA_0(pWD, &pWDSort, &pWDId);

				dTemp = 0.0;
				for(nz=nQ1; nz<=nQ2; nz++)
				{
					dTemp += pWDSort->pMatElement[nz];
				}
				dTemp /= (nQ2-nQ1+1);
				vCT->vSeqData[ni]->vData[nCTMATCol]->pMatElement[nj] = dTemp;
				dTemp = 0.0;
				for(nz=nQ1; nz<=nQ2; nz++)
				{
					dTemp2 = (pWDSort->pMatElement[nz]-vCT->vSeqData[ni]->vData[nCTMATCol]->pMatElement[nj]);
					dTemp += dTemp2*dTemp2;
				}
				dTemp /= (nQ2-nQ1);

				vCT->vSeqData[ni]->vData[nCTMATCol+1]->pMatElement[nj] = sqrt(dTemp);
				vIP->vSeqData[ni]->vData[nIPMATCol+2]->pMatElement[nj] = (vIP->vSeqData[ni]->vData[nIPMATCol]->pMatElement[nj]-vCT->vSeqData[ni]->vData[nCTMATCol]->pMatElement[nj])*sqrt((double)(nQ2-nQ1+1));
				if(nUseVar == 1)
				{
					vIP->vSeqData[ni]->vData[nIPMATCol+2]->pMatElement[nj] /= (vCT->vSeqData[ni]->vData[nCTMATCol+1]->pMatElement[nj]+1e-3);
				}

				DestroyDoubleMatrix(pWD);
				DestroyDoubleMatrix(pWDSort);
				DestroyLongMatrix(pWDId);
			}
		}

		nu += vIP->vSeqData[ni]->nDataNum;
	}

	/* find positive peaks */
	pMATScore = NULL;
	pMATScore = CreateDoubleMatrix(1, nProbeNum);
	if(pMATScore == NULL)
	{
		printf("Error: cannot create memory for MAT score!\n");
		exit(EXIT_FAILURE);
	}

	nTCol = nIPMATCol+2;
	nk = 0;
	for(ni=0; ni<vIP->nSeqNum; ni++)
	{
		for(nj=0; nj<vIP->vSeqData[ni]->nDataNum; nj++)
		{
			pMATScore->pMatElement[nk] = vIP->vSeqData[ni]->vData[nTCol]->pMatElement[nj];
			nk++;
		}
	}

	if(nk != nProbeNum)
	{
		printf("Error: probe number do not match!\n");
		exit(EXIT_FAILURE);
	}

	/* return */
	return pMATScore;
}

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_Resample_Residual()                                          */
/*  Build probe model                                                      */
/* ----------------------------------------------------------------------- */
int TileProbe_Resample_Residual(int nArrayNum, int nProbeNum,
							 int *vGroupID, struct tagBARData **vBARData)
{
	/* define */
	int ni,nj,nk,nx,ny,nz;
	
	/* for IQR shrinking */
	double dM;
	int nGN,nDf;
	
	nk = 0;
	nDf = 0;
	for(ni=0; ni<vBARData[0]->nSeqNum; ni++)
	{
		for(nj=0; nj<vBARData[0]->vSeqData[ni]->nDataNum; nj++)
		{
			nx = 0;
			while(nx < nArrayNum)
			{
				dM = 0.0;
				nGN = 0;
				for(ny=nx; ny<nArrayNum; ny++)
				{
					if(vGroupID[ny] != vGroupID[nx])
						break;
					dM += vBARData[ny]->vSeqData[ni]->vData[1]->pMatElement[nj];
					nGN++;
				}
				dM = dM/nGN;
				for(nz=nx; nz<ny; nz++)
				{
					vBARData[nz]->vSeqData[ni]->vData[1]->pMatElement[nj] = vBARData[nz]->vSeqData[ni]->vData[1]->pMatElement[nj]-dM;
				}
				nx = ny;
			}
			nk++;
		}
	}

	if(nk != nProbeNum)
	{
		printf("Error: TileProbe_Resample_Residual, probe number do not match!\n");
		exit(EXIT_FAILURE);
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  TileProbe_Resample_VarModel()                                          */
/*  Build probe model                                                      */
/* ----------------------------------------------------------------------- */
int TileProbe_Resample_VarModel(int nArrayNum, int nProbeNum, int nShrink,
							 int *vGroupID, struct tagBARData **vBARData, 
							 struct DOUBLEMATRIX *pV)
{
	/* define */
	struct DOUBLEMATRIX *pProbeData;
	
	/* others */
	int ni,nj,nk,nu,nx,ny,nz;
	
	/* for IQR shrinking */
	double dB;
	double dM,dV;
	int nGN,nDf;
	double dQmin,dQmax;

	
	pProbeData = NULL;
	pProbeData = CreateDoubleMatrix(1,nArrayNum);
	if(pProbeData == NULL)
	{
		printf("Error: TileProbe_Resample_VarModel, cannot allocate memory for modeling probe effect!\n");
		exit(EXIT_FAILURE);
	}

	nk = 0;
	nDf = 0;
	for(ni=0; ni<vBARData[0]->nSeqNum; ni++)
	{
		for(nj=0; nj<vBARData[0]->vSeqData[ni]->nDataNum; nj++)
		{
			if(nk%100000 == 0)
			{
				printf("%d ...\n", nk);
			}

			for(nu=0; nu<nArrayNum; nu++)
			{
				pProbeData->pMatElement[nu] = vBARData[nu]->vSeqData[ni]->vData[1]->pMatElement[nj];
			}

			nx = 0;
			dV = 0.0;
			nDf = 0;
			while(nx < nArrayNum)
			{
				dM = 0.0;
				
				nGN = 0;
				for(ny=nx; ny<nArrayNum; ny++)
				{
					if(vGroupID[ny] != vGroupID[nx])
						break;

					dM += pProbeData->pMatElement[ny];
					nGN++;
				}
				dM = dM/nGN;

				for(nz=nx; nz<ny; nz++)
				{
					/* dV += (pProbeData->pMatElement[nz]-dM)*(pProbeData->pMatElement[nz]-dM); */
					dV += pProbeData->pMatElement[nz]*pProbeData->pMatElement[nz];
				}
				nDf = nDf+nGN-1;
				nx = ny;
			}
			dV /= nDf;
			pV->pMatElement[nk] = dV;

			nk++;
		}
	}

	if(nk != nProbeNum)
	{
		printf("Error: TileProbe_Resample_VarModel, probe number do not match!\n");
		exit(EXIT_FAILURE);
	}
	if(nDf == 0)
	{
		printf("Error: TileProbe_Resample_VarModel, no degree of freedom to compute variance!\n");
		exit(EXIT_FAILURE);
	}

	/* Variance shrinkage */
	if(nShrink != 0)
	{
		dM = 0.0;
		dV = 0.0;
		dQmin = 1e10;
		dQmax = -1e10;

		for(nj=0; nj<nProbeNum; nj++)
		{
			if(pV->pMatElement[nj] < dQmin)
				dQmin = pV->pMatElement[nj];
			if(pV->pMatElement[nj] > dQmax)
				dQmax = pV->pMatElement[nj];

			dM += pV->pMatElement[nj];
			dV += pV->pMatElement[nj]*pV->pMatElement[nj];
		}

		dM /= nProbeNum;
		dV = dV-nProbeNum*dM*dM;
		dB = 2.0*(nProbeNum-1)/(2.0+nDf)/nProbeNum+2.0*dM*dM*(nProbeNum-1)/(2.0+nDf)/dV;
		
		if(dB>1.0)
			dB=1.0;
		else if(dB<0.0)
			dB=0.0;

		if(nShrink == -1)
			dB = 0.0;

		for(nj=0; nj<nProbeNum; nj++)
		{
			pV->pMatElement[nj] = sqrt(pV->pMatElement[nj]*(1.0-dB)+dM*dB);
		}
		
		printf("min Var = %f; max Var = %f\n", dQmin, dQmax);
		printf("B = %f\n", dB);
	}
	else
	{
		for(nj=0; nj<nProbeNum; nj++)
		{
			pV->pMatElement[nj] = sqrt(pV->pMatElement[nj]);
		}
	}

	/* release memory */
	DestroyDoubleMatrix(pProbeData);

	/* return */
	return PROC_SUCCESS;
}