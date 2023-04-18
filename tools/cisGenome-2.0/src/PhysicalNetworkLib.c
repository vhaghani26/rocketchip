/* ----------------------------------------------------------------------- */
/*  PhysicalNetworkLib.c : implementation of the physical network library  */
/*  Author : Ji HongKai ; Time: 2006.10                                    */
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
#include "PhysicalNetworkLib.h"

/* ----------------------------------------------------------------------- */ 
/* Network_LinkBINDEntrez_Main:                                            */
/* Link entrez gene id.                                                    */
/* ----------------------------------------------------------------------- */ 
int Network_LinkBINDEntrez_Main(char strBindFile[], char strEntrezFile[], 
						   char strOutFile[])
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	int ni,nj;
	char *chp1,*chp2;
	int nEntrezNum;
	struct INTMATRIX *pGeneId = NULL;
	struct DOUBLEMATRIX *pRNAId = NULL;
	struct DOUBLEMATRIX *pProteinId = NULL;
	struct DOUBLEMATRIX *pGenomicId = NULL;
	struct tagString **vRNAAcc = NULL;
	struct tagString **vProteinAcc = NULL;
	struct tagString **vGenomicAcc = NULL;
	struct DOUBLEMATRIX *pRNASort = NULL;
	struct DOUBLEMATRIX *pProteinSort = NULL;
	struct DOUBLEMATRIX *pGenomicSort = NULL;
	struct LONGMATRIX *pRNASortId = NULL;
	struct LONGMATRIX *pProteinSortId = NULL;
	struct LONGMATRIX *pGenomicSortId = NULL;

	char strLine[LONG_LINE_LENGTH];
	int nTaxId;
	char strStatus[LINE_LENGTH];
	char strRNAAcc[LINE_LENGTH];
	char strRNAGI[LINE_LENGTH];
	char strProteinAcc[LINE_LENGTH];
	char strProteinGI[LINE_LENGTH];
	char strGenomicAcc[LINE_LENGTH];
	char strGenomicGI[LINE_LENGTH];

	int nBINDAcc;
	char strMT1[LINE_LENGTH];
	char strMD1[LINE_LENGTH];
	char strMA1[LINE_LENGTH];
	int nMGI1;
	int nMS1;
	char strMT2[LINE_LENGTH];
	char strMD2[LINE_LENGTH];
	char strMA2[LINE_LENGTH];
	int nMGI2;
	int nMS2;
	int nGID1,nGID2;
	char strIntType[LINE_LENGTH];
	char strTemp1[LINE_LENGTH];
	char strTemp2[LINE_LENGTH];


	/* get entrez database size */
	fpIn = NULL;
	fpIn = fopen(strEntrezFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: Network_LinkBINDEntrez_Main, cannot open entrez gene file!\n");
		exit(EXIT_FAILURE);
	}
	
	nEntrezNum = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		nEntrezNum++;
	}

	fclose(fpIn);

	/* allocate memory */
	pGeneId = CreateIntMatrix(1, nEntrezNum);
	pRNAId = CreateDoubleMatrix(1, nEntrezNum);
	pProteinId = CreateDoubleMatrix(1, nEntrezNum);
	pGenomicId = CreateDoubleMatrix(1, nEntrezNum);
	if( (pGeneId == NULL) || (pRNAId == NULL) || (pProteinId == NULL)
		|| (pGenomicId == NULL) )
	{
		printf("Error: Network_LinkBINDEntrez_Main, cannot allocate memory!\n");
		exit(EXIT_FAILURE);
	}

	vRNAAcc = (struct tagString **)calloc(nEntrezNum, sizeof(struct tagString *));
	vProteinAcc = (struct tagString **)calloc(nEntrezNum, sizeof(struct tagString *));
	vGenomicAcc = (struct tagString **)calloc(nEntrezNum, sizeof(struct tagString *));
	if( (vRNAAcc == NULL) || (vProteinAcc == NULL) || (vGenomicAcc == NULL) )
	{
		printf("Error: Network_LinkBINDEntrez_Main, cannot allocate memory!\n");
		exit(EXIT_FAILURE);
	}

	/* load entrez database */
	fpIn = NULL;
	fpIn = fopen(strEntrezFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: Network_LinkBINDEntrez_Main, cannot open entrez gene file!\n");
		exit(EXIT_FAILURE);
	}
	
	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		chp1 = strLine;
		chp2 = strchr(chp1, '\t');
		*chp2 = '\0';
		nTaxId = atoi(chp1);

		chp1 = chp2+1;
		chp2 = strchr(chp1, '\t');
		*chp2 = '\0';
		pGeneId->pMatElement[ni] = atoi(chp1);

		chp1 = chp2+1;
		chp2 = strchr(chp1, '\t');
		*chp2 = '\0';
		strcpy(strStatus, chp1);

		chp1 = chp2+1;
		chp2 = strchr(chp1, '\t');
		*chp2 = '\0';
		strcpy(strRNAAcc, chp1);
		StrTrimLeft(strRNAAcc);
		chp1 = strchr(strRNAAcc, '.');
		if(chp1 != NULL)
			*chp1 = '\0';
		StringAddTail(vRNAAcc+ni, strRNAAcc);

		chp1 = chp2+1;
		chp2 = strchr(chp1, '\t');
		*chp2 = '\0';
		strcpy(strRNAGI, chp1);
		StrTrimLeft(strRNAGI);
		StrTrimRight(strRNAGI);
		if(strcmp(strRNAGI, "-") == 0)
		{
			pRNAId->pMatElement[ni] = -1;
		}
		else if(strRNAGI[0] == '\0')
		{
			pRNAId->pMatElement[ni] = -1;
		}
		else
		{
			pRNAId->pMatElement[ni] = atoi(strRNAGI);
		}

		chp1 = chp2+1;
		chp2 = strchr(chp1, '\t');
		*chp2 = '\0';
		strcpy(strProteinAcc, chp1);
		StrTrimLeft(strProteinAcc);
		chp1 = strchr(strProteinAcc, '.');
		if(chp1 != NULL)
			*chp1 = '\0';
		StringAddTail(vProteinAcc+ni, strProteinAcc);

		chp1 = chp2+1;
		chp2 = strchr(chp1, '\t');
		*chp2 = '\0';
		strcpy(strProteinGI, chp1);
		StrTrimLeft(strProteinGI);
		StrTrimRight(strProteinGI);
		if(strcmp(strProteinGI, "-") == 0)
		{
			pProteinId->pMatElement[ni] = -1;
		}
		else if(strProteinGI[0] == '\0')
		{
			pProteinId->pMatElement[ni] = -1;
		}
		else
		{
			pProteinId->pMatElement[ni] = atoi(strProteinGI);
		}

		chp1 = chp2+1;
		chp2 = strchr(chp1, '\t');
		*chp2 = '\0';
		strcpy(strGenomicAcc, chp1);
		StrTrimLeft(strGenomicAcc);
		chp1 = strchr(strGenomicAcc, '.');
		if(chp1 != NULL)
			*chp1 = '\0';
		StringAddTail(vGenomicAcc+ni, strGenomicAcc);

		chp1 = chp2+1;
		chp2 = strchr(chp1, '\t');
		*chp2 = '\0';
		strcpy(strGenomicGI, chp1);
		StrTrimLeft(strGenomicGI);
		StrTrimRight(strGenomicGI);
		if(strcmp(strGenomicGI, "-") == 0)
		{
			pGenomicId->pMatElement[ni] = -1;
		}
		else if(strGenomicGI[0] == '\0')
		{
			pGenomicId->pMatElement[ni] = -1;
		}
		else
		{
			pGenomicId->pMatElement[ni] = atoi(strGenomicGI);
		}

		ni++;
	}

	fclose(fpIn);
	
	if(ni != nEntrezNum)
	{
		printf("Error: Network_LinkBINDEntrez_Main, entrez gene number do not match!\n");
		exit(EXIT_FAILURE);
	}

	/* sort */
	DMSORTMERGEA_0(pRNAId, &pRNASort, &pRNASortId);
	DMSORTMERGEA_0(pProteinId, &pProteinSort, &pProteinSortId);
	DMSORTMERGEA_0(pGenomicId, &pGenomicSort, &pGenomicSortId);

	/* convert BIND data */
	fpIn = NULL;
	fpIn = fopen(strBindFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: Network_LinkBINDEntrez_Main, cannot open BIND file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: Network_LinkBINDEntrez_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	nj = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		if(nj%1000 == 0)
		{
			printf("interaction %d...\n", nj);
		}

		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		chp1 = strLine;
		chp2 = strchr(chp1, '\t');
		*chp2 = '\0';
		nBINDAcc = atoi(chp1);

		chp1 = chp2+1;
		chp2 = strchr(chp1, '\t');
		*chp2 = '\0';
		strcpy(strMT1, chp1);
		StrTrimLeft(strMT1);
		StrTrimRight(strMT1);

		chp1 = chp2+1;
		chp2 = strchr(chp1, '\t');
		*chp2 = '\0';
		strcpy(strMD1, chp1);
		StrTrimLeft(strMD1);
		StrTrimRight(strMD1);

		chp1 = chp2+1;
		chp2 = strchr(chp1, '\t');
		*chp2 = '\0';
		strcpy(strMA1, chp1);
		StrTrimLeft(strMA1);
		StrTrimRight(strMA1);

		chp1 = chp2+1;
		chp2 = strchr(chp1, '\t');
		*chp2 = '\0';
		nMGI1 = atoi(chp1);

		chp1 = chp2+1;
		chp2 = strchr(chp1, '\t');
		*chp2 = '\0';
		nMS1 = atoi(chp1);
	
		chp1 = chp2+1;
		chp2 = strchr(chp1, '\t');
		*chp2 = '\0';
		strcpy(strMT2, chp1);
		StrTrimLeft(strMT2);
		StrTrimRight(strMT2);

		chp1 = chp2+1;
		chp2 = strchr(chp1, '\t');
		*chp2 = '\0';
		strcpy(strMD2, chp1);
		StrTrimLeft(strMD2);
		StrTrimRight(strMD2);

		chp1 = chp2+1;
		chp2 = strchr(chp1, '\t');
		*chp2 = '\0';
		strcpy(strMA2, chp1);
		StrTrimLeft(strMA2);
		StrTrimRight(strMA2);

		chp1 = chp2+1;
		chp2 = strchr(chp1, '\t');
		*chp2 = '\0';
		nMGI2 = atoi(chp1);

		chp1 = chp2+1;
		chp2 = strchr(chp1, '\t');
		if(chp2 != NULL)
			*chp2 = '\0';
		nMS2 = atoi(chp1);

		/* search molecule one */
		nGID1 = Network_LinkBINDEntrez_GetGeneID(nMGI1, strMT1, strMA1, 
			nEntrezNum, pGeneId, 
			pRNASortId, pRNASort, pRNAId, vRNAAcc, 
			pProteinSortId, pProteinSort, pProteinId, vProteinAcc,
			pGenomicSortId, pGenomicSort, pGenomicId, vGenomicAcc);

		/* search molecule two */
		nGID2 = Network_LinkBINDEntrez_GetGeneID(nMGI2, strMT2, strMA2, 
			nEntrezNum, pGeneId, 
			pRNASortId, pRNASort, pRNAId, vRNAAcc, 
			pProteinSortId, pProteinSort, pProteinId, vProteinAcc,
			pGenomicSortId, pGenomicSort, pGenomicId, vGenomicAcc);

		/* write */
		if(strcmp(strMT1, "protein") == 0)
		{
			strcpy(strTemp1, "p");
		}
		else if(strcmp(strMT1, "DNA") == 0)
		{
			strcpy(strTemp1, "d");
		}
		else if(strcmp(strMT1, "RNA") == 0)
		{
			strcpy(strTemp1, "r");
		}
		else
		{
			strcpy(strTemp1, "-");
		}

		if(strcmp(strMT2, "protein") == 0)
		{
			strcpy(strTemp2, "p");
		}
		else if(strcmp(strMT2, "DNA") == 0)
		{
			strcpy(strTemp2, "d");
		}
		else if(strcmp(strMT2, "RNA") == 0)
		{
			strcpy(strTemp2, "r");
		}
		else
		{
			strcpy(strTemp2, "-");
		}

		sprintf(strIntType, "%s%s", strTemp1, strTemp2);
		if( (strcmp(strTemp1, "-") != 0) && (strcmp(strTemp2, "-") != 0) &&
			(nGID1 >= 0) && (nGID2 >=0) )
		{
			fprintf(fpOut, "%d\t%s\t%d\n", nGID1, strIntType, nGID2);
		}

		nj++;
	}

	fclose(fpIn);
	fclose(fpOut);

	/* release memory */
	DestroyDoubleMatrix(pRNASort);
	DestroyDoubleMatrix(pProteinSort);
	DestroyDoubleMatrix(pGenomicSort);
	DestroyLongMatrix(pRNASortId);
	DestroyLongMatrix(pProteinSortId);
	DestroyLongMatrix(pGenomicSortId);

	DestroyIntMatrix(pGeneId);
	DestroyDoubleMatrix(pRNAId);
	DestroyDoubleMatrix(pProteinId);
	DestroyDoubleMatrix(pGenomicId);
	for(ni=0; ni<nEntrezNum; ni++)
	{
		DeleteString(vRNAAcc[ni]);
		DeleteString(vProteinAcc[ni]);
		DeleteString(vGenomicAcc[ni]);
		vRNAAcc[ni] = NULL;
		vProteinAcc[ni] = NULL;
		vGenomicAcc[ni] = NULL;
	}
	free(vRNAAcc);
	free(vProteinAcc);
	free(vGenomicAcc);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/* Network_LinkBINDEntrez_GetGeneID:                                       */
/* Get entrez gene id.                                                     */
/* ----------------------------------------------------------------------- */ 
int Network_LinkBINDEntrez_GetGeneID(int nMGI, char strMT[], char strMA[], 
			int nGeneNum, struct INTMATRIX *pGeneId,
			struct LONGMATRIX *pRNASortId, struct DOUBLEMATRIX *pRNASort,
			struct DOUBLEMATRIX *pRNAId, struct tagString **vRNAAcc, 
			struct LONGMATRIX *pProteinSortId, struct DOUBLEMATRIX *pProteinSort,
			struct DOUBLEMATRIX *pProteinId, struct tagString **vProteinAcc,
			struct LONGMATRIX *pGenomicSortId, struct DOUBLEMATRIX *pGenomicSort,
			struct DOUBLEMATRIX *pGenomicId, struct tagString **vGenomicAcc)
{
	/* define */
	int nGeneID = -1;
	int ni,nj,nk;
	int nRID,nDID,nPID;

	/* search RNA */
	/* if(strcmp(strMT, "RNA") == 0)
	{ */
		ni = 0;
		nj = nGeneNum-1;
		nk = (ni+nj)/2;

		while(nj-ni>1)
		{
			if( (int)(pRNASort->pMatElement[nk]) < nMGI)
			{
				ni = nk;
			}
			else
			{
				nj = nk;
			}
			nk = (ni+nj)/2;
		}

		if( (int)(pRNASort->pMatElement[ni]) == nMGI)
			nk = ni;
		else if((int)(pRNASort->pMatElement[nj]) == nMGI)
			nk = nj;
		else
			nk = -1;

		if(nk > 0)
		{
			ni = pRNASortId->pMatElement[nk];
			if(nMGI != (int)(pRNAId->pMatElement[ni]) )
			{
				printf("Error: matching error!\n");
				exit(EXIT_FAILURE);
			}
			if(strcmp(strMA, vRNAAcc[ni]->m_pString) != 0)
			{
				printf("Warning %d: %s vs %s not match!\n", nMGI, strMA, vRNAAcc[ni]->m_pString);
				/* nRID = -1; */
				nRID = pGeneId->pMatElement[ni];
			}
			else
			{
				nRID = pGeneId->pMatElement[ni];
			}
		}
		else
		{
			nRID = -1;
		}
	/* } */

	/* search Protein */
	/* else if(strcmp(strMT, "protein") == 0)
	{ */
		ni = 0;
		nj = nGeneNum-1;
		nk = (ni+nj)/2;

		while(nj-ni>1)
		{
			if( (int)(pProteinSort->pMatElement[nk]) < nMGI)
			{
				ni = nk;
			}
			else
			{
				nj = nk;
			}
			nk = (ni+nj)/2;
		}

		if( (int)(pProteinSort->pMatElement[ni]) == nMGI)
			nk = ni;
		else if((int)(pProteinSort->pMatElement[nj]) == nMGI)
			nk = nj;
		else
			nk = -1;

		if(nk > 0)
		{
			ni = pProteinSortId->pMatElement[nk];
			if(nMGI != (int)(pProteinId->pMatElement[ni]) )
			{
				printf("Error: matching error!\n");
				exit(EXIT_FAILURE);
			}
			if(strcmp(strMA, vProteinAcc[ni]->m_pString) != 0)
			{
				printf("Warning %d: %s vs %s not match!\n", nMGI, strMA, vProteinAcc[ni]->m_pString);
				/* nPID = -1; */
				nPID = pGeneId->pMatElement[ni];
			}
			else
			{
				nPID = pGeneId->pMatElement[ni];
			}
		}
		else
		{
			nPID = -1;
		}
	/* } */

	/* search DNA */
	/* else if(strcmp(strMT, "DNA") == 0)
	{ */
		ni = 0;
		nj = nGeneNum-1;
		nk = (ni+nj)/2;

		while(nj-ni>1)
		{
			if( (int)(pGenomicSort->pMatElement[nk]) < nMGI)
			{
				ni = nk;
			}
			else
			{
				nj = nk;
			}
			nk = (ni+nj)/2;
		}

		if( (int)(pGenomicSort->pMatElement[ni]) == nMGI)
			nk = ni;
		else if((int)(pGenomicSort->pMatElement[nj]) == nMGI)
			nk = nj;
		else
			nk = -1;

		if(nk > 0)
		{
			ni = pGenomicSortId->pMatElement[nk];
			if(nMGI != (int)(pGenomicId->pMatElement[ni]) )
			{
				printf("Error: matching error!\n");
				exit(EXIT_FAILURE);
			}
			if(strcmp(strMA, vGenomicAcc[ni]->m_pString) != 0)
			{
				printf("Warning %d: %s vs %s not match!\n", nMGI, strMA, vGenomicAcc[ni]->m_pString);
				/* nDID = -1; */
				nDID = pGeneId->pMatElement[ni];
			}
			else
			{
				nDID = pGeneId->pMatElement[ni];
			}
		}
		else
		{
			nDID = -1;
		}
	/* } */

	/* else
	{
		nGeneID = -1;
		return nGeneID;
	} */

	if((nRID == -1) && (nPID == -1) && (nDID == -1))
	{
		nGeneID = -1;
	}
	else
	{
		if(nRID > 0)
			nGeneID = nRID;
		if(nPID > 0)
			nGeneID = nPID;
		if(nDID > 0)
			nGeneID = nDID;
	}

	/* return */
	return nGeneID;
}

/* ----------------------------------------------------------------------- */ 
/* Network_GetSubNet_Main:                                                 */
/* Get subnetwork that contains given nodes.                               */
/* ----------------------------------------------------------------------- */ 
int Network_GetSubNet_Main(char strNetWorkFile[], char strNodeFile[], 
						   int nExportAnnotation, char strAnnotationFile[], 
						   char strOutFile[])
{
	/* define */
	struct NETWORKPHYSICAL *pNet = NULL;
	struct INTMATRIX *pNodeId = NULL;
	FILE *fpOut;
	int nTemp;
	int ni,nId;
	int nNodeId1,nNodeId2;

	/* load data */
	pNet = Network_LoadFromSIF(strNetWorkFile);
	if(pNet == NULL)
	{
		printf("Warning: Network_GetSubNet_Main, Empty network!\n");
		return PROC_SUCCESS;
	}
	if(nExportAnnotation == 1)
		Network_LoadAnnotation(pNet, strAnnotationFile);
	
	pNodeId = IMLOAD(strNodeFile);
	if(pNodeId == NULL)
	{
		printf("Warning: Network_GetSubNet_Main, Empty target nodes!\n");
		return PROC_SUCCESS;
	}

	if( (pNodeId->nWidth == 1) && (pNodeId->nHeight > 1) )
	{
		nTemp = pNodeId->nHeight;
		pNodeId->nHeight = pNodeId->nWidth;
		pNodeId->nWidth = nTemp;
	}

	/* label target nodes */
	for(ni=0; ni<pNodeId->nWidth; ni++)
	{
		nId = pNodeId->pMatElement[ni];
		if( (nId<0) || (nId>=pNet->nMaxNodeNum) )
		{
			continue;
		}
		if(pNet->vNodes[nId] == NULL)
		{
			continue;
		}

		pNet->vNodes[nId]->nTag = ni;
	}

	/* export */
	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: Network_GetSubNet_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<pNet->nEdgeNum; ni++)
	{
		nNodeId1 = pNet->vEdges[ni]->nNodeID1;
		nNodeId2 = pNet->vEdges[ni]->nNodeID2;
		if( (pNet->vNodes[nNodeId1] != NULL) && (pNet->vNodes[nNodeId2] != NULL) )
		{
			if( (pNet->vNodes[nNodeId1]->nTag >=0) && (pNet->vNodes[nNodeId2]->nTag >=0) )
			{
				if(nExportAnnotation == 1)
				{
					fprintf(fpOut, "%s\tpp\t%s\n", pNet->vNodes[nNodeId1]->strName->m_pString, 
						pNet->vNodes[nNodeId2]->strName->m_pString);
				}
				else
				{
					fprintf(fpOut, "%d\tpp\t%d\n", pNet->vEdges[ni]->nNodeID1, pNet->vEdges[ni]->nNodeID2); 
				}
			}
		}
	}

	fclose(fpOut);



	/* release memory */
	NETWORKPHYSICALDESTROY(&pNet);
	DestroyIntMatrix(pNodeId);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/* Network_CreateOrthogloNet_Main:                                         */
/* Create an ortholog network.                                             */
/* ----------------------------------------------------------------------- */ 
int Network_CreateOrthogloNet_Main(char strNetWorkFile[], char strHomoloFile[], 
								   int nSrcSpecies, int nDestSpecies, 
								   char strOutFile[])
{
	/* define */
	struct NETWORKEDGE *pEdge;
	struct NETWORKNODE *pNode1,*pNode2;
	struct NETWORKPHYSICAL *pNet = NULL;
	FILE *fpIn;
	FILE *fpOut;
	char strLine[LONG_LINE_LENGTH];
	int nHID,nSID,nGID;
	char strGeneName[LINE_LENGTH];
	int nLastHID = -1;
	int nSrcOK,nDstOK;
	int nSrcID,nDstID;
	int ni;

	/* load data */
	pNet = Network_LoadFromSIF(strNetWorkFile);
	if(pNet == NULL)
	{
		printf("Warning: Empty network!\n");
		return PROC_SUCCESS;
	}

	/* open file */
	fpIn = NULL;
	fpIn = fopen(strHomoloFile, "r");
	if(fpIn == NULL)
	{
		printf("Warning: Empty ortholog mapping!\n");
		NETWORKPHYSICALDESTROY(&pNet);
		return PROC_SUCCESS;
	}

	/* map nodes to orthologs */
	nSrcOK = 0;
	nDstOK = 0;
	nSrcID = -1;
	nDstID = -1;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%d %d %d %s", &nHID, &nSID, &nGID, strGeneName);
		if(nHID != nLastHID)
		{
			if( (nSrcOK==1) && (nDstOK == 1) )
			{
				if(nSrcID < pNet->nMaxNodeNum)
				{
					if(pNet->vNodes[nSrcID] != NULL)
					{
						pNet->vNodes[nSrcID]->nTag = nDstID;
					}
				}
			}
			nLastHID = nHID;
			nSrcOK = 0;
			nDstOK = 0;
			nSrcID = -1;
			nDstID = -1;
		}
		
		if(nSID == nSrcSpecies)
		{
			nSrcID = nGID;
			nSrcOK = 1;
		}
		else if(nSID == nDestSpecies)
		{
			nDstID = nGID;
			nDstOK = 1;
		}
	}

	if( (nSrcOK==1) && (nDstOK == 1) )
	{
		if(nSrcID < pNet->nMaxNodeNum)
		{
			if(pNet->vNodes[nSrcID] != NULL)
			{
				pNet->vNodes[nSrcID]->nTag = nDstID;
			}
		}
	}
	nLastHID = nHID;
	nSrcOK = 0;
	nDstOK = 0;
	nSrcID = -1;
	nDstID = -1;

	/* close file */
	fclose(fpIn);

	/* map edges to orthologs */
	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: Empty ortholog mapping!\n");
		NETWORKPHYSICALDESTROY(&pNet);
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<pNet->nEdgeNum; ni++)
	{
		pEdge = pNet->vEdges[ni];
		pNode1 = pNet->vNodes[pEdge->nNodeID1];
		pNode2 = pNet->vNodes[pEdge->nNodeID2];
		if( (pNode1->nTag<0) || (pNode2->nTag<0) )
			continue;

		fprintf(fpOut, "%d\tpp\t%d\n", pNode1->nTag, pNode2->nTag);
	}

	/* close file */
	fclose(fpOut);

	/* release memory */
	NETWORKPHYSICALDESTROY(&pNet);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/* Network_FindShortestPath_Main:                                          */
/* Find shortest path of a network.                                        */
/* ----------------------------------------------------------------------- */ 
int Network_FindShortestPath_Main(char strNetWorkFile[], char strAnnotationFile[],
								  char strSourceNodeFile[], char strDestNodeFile[], 
								  char strOutFile[], int nMaxIterNum)
{
	/* define */
	struct NETWORKSHORTESTPATH *pPath = NULL;
	struct NETWORKPHYSICAL *pNet = NULL;
	struct INTMATRIX *pSourceNodes = NULL;
	struct INTMATRIX *pDestNodes = NULL;
	long nTemp;

	/* load data */
	pNet = Network_LoadFromSIF(strNetWorkFile);
	Network_LoadAnnotation(pNet, strAnnotationFile);
	pSourceNodes = IMLOAD(strSourceNodeFile);
	pDestNodes = IMLOAD(strDestNodeFile);
	if( (pNet == NULL) || (pSourceNodes == NULL) || (pDestNodes == NULL) )
	{
		printf("Error: Network_FindShortestPath_Main, empty network or nodes!\n");
		exit(EXIT_FAILURE);
	}
	if( (pSourceNodes->nWidth == 1) && (pSourceNodes->nHeight > 1) )
	{
		nTemp = pSourceNodes->nHeight;
		pSourceNodes->nHeight = pSourceNodes->nWidth;
		pSourceNodes->nWidth = nTemp;
	}
	if( (pDestNodes->nWidth == 1) && (pDestNodes->nHeight > 1) )
	{
		nTemp = pDestNodes->nHeight;
		pDestNodes->nHeight = pDestNodes->nWidth;
		pDestNodes->nWidth = nTemp;
	}

	/* find shortest path */
	pPath = Network_FindShortestPath(pNet, pDestNodes, nMaxIterNum);
	if(pPath == NULL)
	{
		printf("Warning: Network_FindShortestPath_Main, empty path objects!\n");
		
		/* release memory */
		DestroyIntMatrix(pSourceNodes);
		DestroyIntMatrix(pDestNodes);
		NETWORKPHYSICALDESTROY(&pNet);

		return PROC_SUCCESS;
	}

	/* export results */
	Network_FindShortestPath_ExportResults(strOutFile, pNet, pPath, 
				pSourceNodes, pDestNodes);


	/* release memory */
	DestroyIntMatrix(pSourceNodes);
	DestroyIntMatrix(pDestNodes);
	NETWORKPHYSICALDESTROY(&pNet);
	NETWORKSHORTESTPATHDESTROY(&pPath);


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/* Network_FindShortestPath:                                               */
/* Find shortest path.                                                     */
/* ----------------------------------------------------------------------- */ 
struct NETWORKSHORTESTPATH *Network_FindShortestPath(struct NETWORKPHYSICAL *pNet, 
					struct INTMATRIX *pDestNode, int nMaxIterNum)
{
	/* define */
	struct NETWORKSHORTESTPATH *pPath = NULL;
	int nMaxIter = 0;
	int nDestNum = 0;
	int ni,nj,nk,nID,nNeighborID,nDistID;
	struct NETWORKNODE *pNode;
	struct NETWORKEDGE *pEdge;
	struct NODEEDGEINTERFACE *pInterface;


	/* init check */
	if( (pNet == NULL) || (pDestNode == NULL) || (nMaxIterNum <= 0) )
	{
		printf("Warning: Network_FindShortestPath, empty network, destination nodes or maximum iterations!\n");
		return NULL;
	}

	if( (pNet->nRealNodeNum <= 0) || (pDestNode->nWidth <= 0) )
	{
		printf("Warning: Network_FindShortestPath, empty network or destination nodes!\n");
		return NULL;
	}

	nMaxIter = pNet->nRealNodeNum;
	if(nMaxIter > nMaxIterNum)
		nMaxIter = nMaxIterNum;
	nDestNum = pDestNode->nWidth;

	/* create path */
	pPath = NETWORKSHORTESTPATHCREATE(pNet->nRealNodeNum, nDestNum);
	
	/* init target node */
	for(ni=0; ni<pDestNode->nWidth; ni++)
	{
		nID = pDestNode->pMatElement[ni];
		if( (nID<0) || (nID>=pNet->nMaxNodeNum) )
		{
			printf("Warning: Network_FindShortestPath, node %d not covered by the network!\n", nID);
			continue;
		}
		if(pNet->vNodes[nID] != NULL)
			pNet->vNodes[nID]->nTag = ni;
	}

	/* init node id */
	Network_FindShortestPath_InitDist(pNet, pPath);

	/* find shorted distance */
	for(ni=0; ni<nMaxIter; ni++)
	{
		for(nj=0; nj<pPath->nNodeNum; nj++)
		{
			nID = pPath->pNodeID->pMatElement[nj];
			if(pNet->vNodes[nID] == NULL)
			{
				printf("Error: Network_FindShortestPath, node ids do not match!\n");
				exit(EXIT_FAILURE);
			}
			if(nID != pNet->vNodes[nID]->nID)
			{
				printf("Error: Network_FindShortestPath, node ids do not match!\n");
				exit(EXIT_FAILURE);
			}
			pNode = pNet->vNodes[nID];

			pInterface = pNode->pEdgeList;
			while(pInterface != NULL)
			{
				pEdge = pNet->vEdges[pInterface->nEdgeID];
				if(nID == pEdge->nNodeID1)
					nNeighborID = pEdge->nNodeID2;
				else
					nNeighborID = pEdge->nNodeID1;

				nDistID = pNet->vNodes[nNeighborID]->nNetIndex;
				if(pPath->pNodeID->pMatElement[nDistID] != nNeighborID)
				{
					printf("Error: Network_FindShortestPath, node ids do not match!\n");
					exit(EXIT_FAILURE);
				}

				for(nk=0; nk<pPath->vShortDist[nj]->nWidth; nk++)
				{
					if(pPath->vShortDist[nj]->pMatElement[nk] > (1.0+pPath->vShortDist[nDistID]->pMatElement[nk]))
					{
						pPath->vShortDist[nj]->pMatElement[nk] = (1.0+pPath->vShortDist[nDistID]->pMatElement[nk]);
						pPath->vNeighborNode[nj]->pMatElement[nk] = nNeighborID;
					}
				}

				pInterface = pInterface->pNext;
			}
		}

	}

	/* return */
	return pPath;
}

/* ----------------------------------------------------------------------- */ 
/* Network_FindShortestPath_InitDist:                                      */
/* Init distance for shortest path finding.                                */
/* ----------------------------------------------------------------------- */ 
int Network_FindShortestPath_InitDist(struct NETWORKPHYSICAL *pNet, struct NETWORKSHORTESTPATH *pPath)
{
	/* define */
	struct NETWORKNODE *pNode;
	int ni,nj;
	double dMaxDist;

	/* init */
	if( (pNet == NULL) || (pPath == NULL) )
	{
		printf("Warning: empty network or path element!\n");
		return PROC_SUCCESS;
	}
	dMaxDist = (double)(pNet->nRealNodeNum+1);

	/* set values */
	ni = 0;
	pNode = pNet->pHeadNode;
	while(pNode != NULL)
	{
		/* set node ID */
		if(ni >= pPath->nNodeNum)
		{
			printf("Error: net/path node number do not match!\n");
			exit(EXIT_FAILURE);
		}

		pPath->pNodeID->pMatElement[ni] = pNode->nID;

		/* set init path */
		if(pPath->vNeighborNode[ni]->nWidth != pPath->vShortDist[ni]->nWidth)
		{
			printf("Error: path register lengths do not match!\n");
			exit(EXIT_FAILURE);
		}

		for(nj=0; nj<pPath->vNeighborNode[ni]->nWidth; nj++)
		{
			pPath->vNeighborNode[ni]->pMatElement[nj] = -1;
			pPath->vShortDist[ni]->pMatElement[nj] = dMaxDist;
		}

		/* set target id */
		if(pNode->nTag >= 0)
		{
			pPath->vShortDist[ni]->pMatElement[pNode->nTag] = 0;
		}

		ni++;
		pNode = pNode->pNext;
	}

	if(ni != pPath->nNodeNum)
	{
		printf("Error: node numbers do not match!\n");
		exit(EXIT_FAILURE);
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/* Network_FindShortestPath_ExportResults:                                 */
/* Export shortest path.                                                   */
/* ----------------------------------------------------------------------- */ 
int Network_FindShortestPath_ExportResults(char strFileName[], 
				struct NETWORKPHYSICAL *pNet, struct NETWORKSHORTESTPATH *pPath, 
				struct INTMATRIX *pSourceNode, struct INTMATRIX *pDestNode)
{
	/* define */
	FILE *fpOut;
	int ni,nj;
	int nPathID;
	int nNodeID1,nNodeID2;
	int nLastNodeID;
	int nTempPathID;
	int nPathScore;

	/* init check */
	if( (pNet == NULL) || (pPath == NULL) || (pSourceNode == NULL) || (pDestNode == NULL))
	{
		printf("Warning: Network_FindShortestPath_ExportResults, empty network components!\n");
		return PROC_SUCCESS;
	}

	/* open file */
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: Network_FindShortestPath_ExportResults, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	/* write */
	for(ni=0; ni<pSourceNode->nWidth; ni++)
	{
		nNodeID1 = pSourceNode->pMatElement[ni];
		if( (nNodeID1 < 0) || (nNodeID1 >= pNet->nMaxNodeNum) )
		{
			printf("Warning: source node %d is not covered by the network!\n", nNodeID1);
			continue;
		}
		if(pNet->vNodes[nNodeID1] == NULL)
		{
			printf("Warning: source node %d is not covered by the network!\n", nNodeID1);
			continue;
		}
		nPathID = pNet->vNodes[nNodeID1]->nNetIndex;
		
		for(nj=0; nj<pDestNode->nWidth; nj++)
		{
			if(pPath->vShortDist[nPathID]->pMatElement[nj] >= (double)(pNet->nRealNodeNum+0.5))
				continue;

			nPathScore = 0;
			nNodeID2 = pDestNode->pMatElement[nj];
			if( (nNodeID2 < 0) || (nNodeID2 >= pNet->nMaxNodeNum) )
			{
				printf("Warning: source node %d is not covered by the network!\n", nNodeID2);
				continue;
			}
			if(pNet->vNodes[nNodeID2] == NULL)
			{
				printf("Warning: source node %d is not covered by the network!\n", nNodeID2);
				continue;
			}

			fprintf(fpOut, "%s\t%s\t%f\t", pNet->vNodes[nNodeID1]->strName->m_pString, 
				pNet->vNodes[nNodeID2]->strName->m_pString, pPath->vShortDist[nPathID]->pMatElement[nj]); 

			if(pNet->vNodes[nNodeID1]->nTag >= 0)
			{
				fprintf(fpOut, "(%s)", pNet->vNodes[nNodeID1]->strName->m_pString);
				nPathScore++;
			}
			else
			{
				fprintf(fpOut, "%s", pNet->vNodes[nNodeID1]->strName->m_pString);
			}
			
			nLastNodeID = pPath->vNeighborNode[nPathID]->pMatElement[nj];
			while(nLastNodeID != -1)
			{
				if(pNet->vNodes[nLastNodeID]->nTag >= 0)
				{
					fprintf(fpOut, "-->(%s)", pNet->vNodes[nLastNodeID]->strName->m_pString);
					nPathScore++;
				}
				else
				{
					fprintf(fpOut, "-->%s", pNet->vNodes[nLastNodeID]->strName->m_pString);
				}
				nTempPathID = pNet->vNodes[nLastNodeID]->nNetIndex;
				nLastNodeID = pPath->vNeighborNode[nTempPathID]->pMatElement[nj];
			}

			fprintf(fpOut, "\t%d\n", nPathScore);
		}

	}

	/* close file */
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/* Network_FindShortestPath_V:                                             */
/* Find shortest path.                                                     */
/* ----------------------------------------------------------------------- */ 
struct NETWORKSHORTESTPATH *Network_FindShortestPath_V(struct NETWORKPHYSICAL *pNet, 
					struct INTMATRIX *pDestNode, int nRow, int nCol, double dMaxNodeValue,
					int nMaxIterNum)
{
	/* define */
	struct NETWORKSHORTESTPATH *pPath = NULL;
	int nMaxIter = 0;
	int nDestNum = 0;
	int ni,nj,nk,nID,nNeighborID,nDistID;
	struct NETWORKNODE *pNode;
	struct NETWORKEDGE *pEdge;
	struct NODEEDGEINTERFACE *pInterface;
	double dTemp;

	/* init check */
	if( (pNet == NULL) || (pDestNode == NULL) || (nMaxIterNum <= 0) )
	{
		printf("Warning: Network_FindShortestPath_V, empty network, destination nodes or maximum iterations!\n");
		return NULL;
	}

	if( (pNet->nRealNodeNum <= 0) || (pDestNode->nWidth <= 0) )
	{
		printf("Warning: Network_FindShortestPath_V, empty network or destination nodes!\n");
		return NULL;
	}

	nMaxIter = pNet->nRealNodeNum;
	if(nMaxIter > nMaxIterNum)
		nMaxIter = nMaxIterNum;
	nDestNum = pDestNode->nWidth;

	/* create path */
	pPath = NETWORKSHORTESTPATHCREATE(pNet->nRealNodeNum, nDestNum);
	
	/* init target node */
	for(ni=0; ni<pDestNode->nWidth; ni++)
	{
		nID = pDestNode->pMatElement[ni];
		if( (nID<0) || (nID>=pNet->nMaxNodeNum) )
		{
			printf("Warning: Network_FindShortestPath, node %d not covered by the network!\n", nID);
			continue;
		}
		if(pNet->vNodes[nID] != NULL)
		{
			pNet->vNodes[nID]->nTag = ni;
		}
	}

	/* init node id */
	Network_FindShortestPath_V_InitDist(pNet, nRow, nCol, pPath, dMaxNodeValue);

	/* find shorted distance */
	for(ni=0; ni<nMaxIter; ni++)
	{
		for(nj=0; nj<pPath->nNodeNum; nj++)
		{
			nID = pPath->pNodeID->pMatElement[nj];
			if(pNet->vNodes[nID] == NULL)
			{
				printf("Error: Network_FindShortestPath, node ids do not match!\n");
				exit(EXIT_FAILURE);
			}
			if(nID != pNet->vNodes[nID]->nID)
			{
				printf("Error: Network_FindShortestPath, node ids do not match!\n");
				exit(EXIT_FAILURE);
			}
			pNode = pNet->vNodes[nID];
			dTemp = DMGETAT(pNode->pDV, nRow, nCol);

			pInterface = pNode->pEdgeList;
			while(pInterface != NULL)
			{
				pEdge = pNet->vEdges[pInterface->nEdgeID];
				if(nID == pEdge->nNodeID1)
					nNeighborID = pEdge->nNodeID2;
				else
					nNeighborID = pEdge->nNodeID1;

				nDistID = pNet->vNodes[nNeighborID]->nNetIndex;
				if(pPath->pNodeID->pMatElement[nDistID] != nNeighborID)
				{
					printf("Error: Network_FindShortestPath, node ids do not match!\n");
					exit(EXIT_FAILURE);
				}

				for(nk=0; nk<pPath->vShortDist[nj]->nWidth; nk++)
				{
					if(pPath->vShortDist[nj]->pMatElement[nk] > (dTemp+pPath->vShortDist[nDistID]->pMatElement[nk]))
					{
						pPath->vShortDist[nj]->pMatElement[nk] = (dTemp+pPath->vShortDist[nDistID]->pMatElement[nk]);
						pPath->vNeighborNode[nj]->pMatElement[nk] = nNeighborID;
					}
				}

				pInterface = pInterface->pNext;
			}
		}

	}

	/* return */
	return pPath;
}

/* ----------------------------------------------------------------------- */ 
/* Network_FindShortestPath_V_ExportResults:                               */
/* Export shortest path.                                                   */
/* ----------------------------------------------------------------------- */ 
int Network_FindShortestPath_V_ExportResults(char strFileName[], 
				struct NETWORKPHYSICAL *pNet, int nRow, int nCol, double dMaxNodeValue,
				struct NETWORKSHORTESTPATH *pPath, 
				struct INTMATRIX *pSourceNode, struct INTMATRIX *pDestNode)
{
	/* define */
	FILE *fpOut;
	int ni,nj;
	int nPathID;
	int nNodeID1,nNodeID2;
	int nLastNodeID;
	int nTempPathID;
	double dPathScore,dTemp;

	/* init check */
	if( (pNet == NULL) || (pPath == NULL) || (pSourceNode == NULL) || (pDestNode == NULL))
	{
		printf("Warning: Network_FindShortestPath_ExportResults, empty network components!\n");
		return PROC_SUCCESS;
	}

	/* open file */
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: Network_FindShortestPath_ExportResults, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	/* write */
	for(ni=0; ni<pSourceNode->nWidth; ni++)
	{
		nNodeID1 = pSourceNode->pMatElement[ni];
		if( (nNodeID1 < 0) || (nNodeID1 >= pNet->nMaxNodeNum) )
		{
			printf("Warning: source node %d is not covered by the network!\n", nNodeID1);
			continue;
		}
		if(pNet->vNodes[nNodeID1] == NULL)
		{
			printf("Warning: source node %d is not covered by the network!\n", nNodeID1);
			continue;
		}
		nPathID = pNet->vNodes[nNodeID1]->nNetIndex;
		
		for(nj=0; nj<pDestNode->nWidth; nj++)
		{
			if(pPath->vShortDist[nPathID]->pMatElement[nj] >= ((double)(pNet->nRealNodeNum)*dMaxNodeValue+0.5))
				continue;

			dPathScore = 0.0;
			nNodeID2 = pDestNode->pMatElement[nj];
			if( (nNodeID2 < 0) || (nNodeID2 >= pNet->nMaxNodeNum) )
			{
				printf("Warning: source node %d is not covered by the network!\n", nNodeID2);
				continue;
			}
			if(pNet->vNodes[nNodeID2] == NULL)
			{
				printf("Warning: source node %d is not covered by the network!\n", nNodeID2);
				continue;
			}

			fprintf(fpOut, "%s\t%s\t%f\t", pNet->vNodes[nNodeID1]->strName->m_pString, 
				pNet->vNodes[nNodeID2]->strName->m_pString, pPath->vShortDist[nPathID]->pMatElement[nj]); 

			dTemp = DMGETAT(pNet->vNodes[nNodeID1]->pDV, nRow, nCol);
			if(pNet->vNodes[nNodeID1]->nTag >= 0)
			{
				fprintf(fpOut, "(%s,%f)", pNet->vNodes[nNodeID1]->strName->m_pString, dTemp);
				dPathScore += dTemp;
			}
			else
			{
				fprintf(fpOut, "%s,%f", pNet->vNodes[nNodeID1]->strName->m_pString, dTemp);
			}
			
			nLastNodeID = pPath->vNeighborNode[nPathID]->pMatElement[nj];
			while(nLastNodeID != -1)
			{
				dTemp = DMGETAT(pNet->vNodes[nNodeID1]->pDV, nRow, nCol);
				if(pNet->vNodes[nLastNodeID]->nTag >= 0)
				{
					fprintf(fpOut, "-->(%s,%f)", pNet->vNodes[nLastNodeID]->strName->m_pString, dTemp);
					dPathScore += dTemp;
				}
				else
				{
					fprintf(fpOut, "-->%s,%f", pNet->vNodes[nLastNodeID]->strName->m_pString, dTemp);
				}
				nTempPathID = pNet->vNodes[nLastNodeID]->nNetIndex;
				nLastNodeID = pPath->vNeighborNode[nTempPathID]->pMatElement[nj];
			}

			fprintf(fpOut, "\t%f\n", dPathScore);
		}

	}

	/* close file */
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/* Network_FindShortestPath_V_InitDist:                                    */
/* Init distance for shortest path finding.                                */
/* ----------------------------------------------------------------------- */ 
int Network_FindShortestPath_V_InitDist(struct NETWORKPHYSICAL *pNet, int nRow, int nCol,
				struct NETWORKSHORTESTPATH *pPath, double dMaxNodeValue)
{
	/* define */
	struct NETWORKNODE *pNode;
	int ni,nj;
	double dMaxDist;

	/* init */
	if( (pNet == NULL) || (pPath == NULL) )
	{
		printf("Warning: empty network or path element!\n");
		return PROC_SUCCESS;
	}
	dMaxDist = (double)(pNet->nRealNodeNum)*dMaxNodeValue+1.0;

	/* set values */
	ni = 0;
	pNode = pNet->pHeadNode;
	while(pNode != NULL)
	{
		/* set node ID */
		if(ni >= pPath->nNodeNum)
		{
			printf("Error: net/path node number do not match!\n");
			exit(EXIT_FAILURE);
		}

		pPath->pNodeID->pMatElement[ni] = pNode->nID;

		/* set init path */
		if(pPath->vNeighborNode[ni]->nWidth != pPath->vShortDist[ni]->nWidth)
		{
			printf("Error: path register lengths do not match!\n");
			exit(EXIT_FAILURE);
		}

		for(nj=0; nj<pPath->vNeighborNode[ni]->nWidth; nj++)
		{
			pPath->vNeighborNode[ni]->pMatElement[nj] = -1;
			pPath->vShortDist[ni]->pMatElement[nj] = dMaxDist;
		}

		/* set target id */
		if(pNode->nTag >= 0)
		{
			pPath->vShortDist[ni]->pMatElement[pNode->nTag] = DMGETAT(pNode->pDV, nRow, nCol);
		}

		ni++;
		pNode = pNode->pNext;
	}

	if(ni != pPath->nNodeNum)
	{
		printf("Error: node numbers do not match!\n");
		exit(EXIT_FAILURE);
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/* Network_LoadFromSIF:                                                    */
/* Load a physical network from a SIF file.                                */
/* ----------------------------------------------------------------------- */ 
struct NETWORKPHYSICAL *Network_LoadFromSIF(char strFileName[])
{
	/* define */
	struct NETWORKPHYSICAL *pNet = NULL;
	int nMaxNodeNum = 0;
	int nEdgeNum = 0;
	struct NETWORKEDGE *pEdge = NULL;
	struct NETWORKNODE *pNode = NULL;
	struct NETWORKNODE *pPrevNode = NULL;

	/* file manipulation */
	FILE *fpIn;
	char strLine[LONG_LINE_LENGTH];
	int nNodeID1,nNodeID2;
	char strLinkType[LINE_LENGTH];
	int ni;
	char strTemp[LINE_LENGTH];

	/* determine edge number and max node number */
	fpIn = NULL;
	fpIn = fopen(strFileName, "r");
	if(fpIn == NULL)
	{
		printf("Warning: cannot open the file to load network!\n");
		return NULL;
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%d %s %d", &nNodeID1, strLinkType, &nNodeID2);
		if(nNodeID1 >= nMaxNodeNum)
			nMaxNodeNum = nNodeID1+1;
		if(nNodeID2 >= nMaxNodeNum)
			nMaxNodeNum = nNodeID2+1;

		nEdgeNum++;
	}

	fclose(fpIn);

	if(nEdgeNum <= 0)
	{
		printf("Warning: the network is empty!\n");
		return NULL;
	}

	/* create network */
	pNet = NULL;
	pNet = NETWORKPHYSICALCREATE();
	if(pNet == NULL)
	{
		printf("Error: cannot create network!\n");
		exit(EXIT_FAILURE);
	}

	pNet->nMaxNodeNum = nMaxNodeNum;
	pNet->nEdgeNum = nEdgeNum;
	pNet->vNodes = (struct NETWORKNODE **)calloc(nMaxNodeNum, sizeof(struct NETWORKNODE *));
	if(pNet->vNodes == NULL)
	{
		printf("Error: cannot allocate memory for netwok nodes!\n");
		exit(EXIT_FAILURE);
	}
	pNet->vEdges = (struct NETWORKEDGE **)calloc(nEdgeNum, sizeof(struct NETWORKEDGE *));
	if(pNet->vEdges == NULL)
	{
		printf("Error: cannot allocate memory for netwok edges!\n");
		exit(EXIT_FAILURE);
	}

	/* initialize network */
	fpIn = NULL;
	fpIn = fopen(strFileName, "r");
	if(fpIn == NULL)
	{
		printf("Error: cannot open the file to load network!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%d %s %d", &nNodeID1, strLinkType, &nNodeID2);

		/* create edge */
		pEdge = NULL;
		pEdge = NETWORKEDGECREATE();
		if(pEdge == NULL)
		{
			printf("Error: cannot create network edge!\n");
			exit(EXIT_FAILURE);
		}
		pEdge->nEdgeID = ni;
		pEdge->nNodeID1 = nNodeID1;
		pEdge->nNodeID2 = nNodeID2;
		if(strcmp(strLinkType, "pp") == 0)
		{
			pEdge->nInteractionType = 1;
		}
		else if(strcmp(strLinkType, "pd") == 0)
		{
			pEdge->nInteractionType = 2;
		}
		else if(strcmp(strLinkType, "dp") == 0)
		{
			pEdge->nInteractionType = 3;
		}
		else if(strcmp(strLinkType, "pr") == 0)
		{
			pEdge->nInteractionType = 4;
		}
		else if(strcmp(strLinkType, "rp") == 0)
		{
			pEdge->nInteractionType = 5;
		}
		else
		{
			pEdge->nInteractionType = 0;
		}
		pNet->vEdges[ni] = pEdge;

		/* update node1 */
		NETWORKNODE_ADDEDGE(pNet->vNodes+nNodeID1, nNodeID1, ni);

		/* update node2 */
		NETWORKNODE_ADDEDGE(pNet->vNodes+nNodeID2, nNodeID2, ni);

		ni++;
	}

	fclose(fpIn);
	if(ni != nEdgeNum)
	{
		printf("Error: network edge number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* link nodes */
	pPrevNode = NULL;
	for(ni=0; ni<nMaxNodeNum; ni++)
	{
		if(pNet->vNodes[ni] != NULL)
		{
			pNet->vNodes[ni]->nNetIndex = pNet->nRealNodeNum;
			pNet->nRealNodeNum += 1;
			sprintf(strTemp, "%d", pNet->vNodes[ni]->nID);
			StringAddTail(&(pNet->vNodes[ni]->strName), strTemp);
			if(pNet->pHeadNode == NULL)
			{
				pNet->pHeadNode = pNet->vNodes[ni];
			}
			if(pPrevNode != NULL)
			{
				pPrevNode->pNext = pNet->vNodes[ni];
			}
			pPrevNode = pNet->vNodes[ni];
		}
	}

	/* return */
	return pNet;
}

/* ----------------------------------------------------------------------- */ 
/* Network_InitNodeDV:                                                     */
/* Create DV vector for network nodes.                                     */
/* ----------------------------------------------------------------------- */
int Network_InitNodeDV(struct NETWORKPHYSICAL *pNet, int nHeight, int nWidth)
{
	/* define */
	struct NETWORKNODE *pNode;

	/* init */
	if(pNet == NULL)
	{
		return PROC_SUCCESS;
	}

	pNode = pNet->pHeadNode;
	while(pNode != NULL)
	{
		if(pNode->pDV != NULL)
		{
			DestroyDoubleMatrix(pNode->pDV);
			pNode->pDV = NULL;
		}

		pNode->pDV = CreateDoubleMatrix(nHeight, nWidth);
		if(pNode->pDV == NULL)
		{
			printf("Error: Network_InitNodeDV, cannot create value matrix for network!\n");
			exit(EXIT_FAILURE);
		}

		pNode = pNode->pNext;
	}


	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/* Network_InitNodeBV:                                                     */
/* Create BV vector for network nodes.                                     */
/* ----------------------------------------------------------------------- */
int Network_InitNodeBV(struct NETWORKPHYSICAL *pNet, int nHeight, int nWidth)
{
	/* define */
	struct NETWORKNODE *pNode;

	/* init */
	if(pNet == NULL)
	{
		return PROC_SUCCESS;
	}

	pNode = pNet->pHeadNode;
	while(pNode != NULL)
	{
		if(pNode->pBV != NULL)
		{
			DestroyByteMatrix(pNode->pBV);
			pNode->pBV = NULL;
		}

		pNode->pBV = CreateByteMatrix(nHeight, nWidth);
		if(pNode->pBV == NULL)
		{
			printf("Error: Network_InitNodeBV, cannot create value matrix for network!\n");
			exit(EXIT_FAILURE);
		}

		pNode = pNode->pNext;
	}


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/* Network_InitNodeIV:                                                     */
/* Create IV vector for network nodes.                                     */
/* ----------------------------------------------------------------------- */
int Network_InitNodeIV(struct NETWORKPHYSICAL *pNet, int nHeight, int nWidth)
{
	/* define */
	struct NETWORKNODE *pNode;

	/* init */
	if(pNet == NULL)
	{
		return PROC_SUCCESS;
	}

	pNode = pNet->pHeadNode;
	while(pNode != NULL)
	{
		if(pNode->pIV != NULL)
		{
			DestroyIntMatrix(pNode->pIV);
			pNode->pIV = NULL;
		}

		pNode->pIV = CreateIntMatrix(nHeight, nWidth);
		if(pNode->pIV == NULL)
		{
			printf("Error: Network_InitNodeIV, cannot create value matrix for network!\n");
			exit(EXIT_FAILURE);
		}

		pNode = pNode->pNext;
	}


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/* Network_InitNodeSV:                                                     */
/* Create SV vector for network nodes.                                     */
/* ----------------------------------------------------------------------- */
int Network_InitNodeSV(struct NETWORKPHYSICAL *pNet, int nSVnum)
{
	/* define */
	struct NETWORKNODE *pNode;
	int ni;

	/* init */
	if(pNet == NULL)
	{
		return PROC_SUCCESS;
	}

	pNode = pNet->pHeadNode;
	while(pNode != NULL)
	{
		if(pNode->vSV != NULL)
		{
			for(ni=0; ni<pNode->nSVnum; ni++)
			{
				DeleteString(pNode->vSV[ni]);
				pNode->vSV[ni] = NULL;
			}
			free(pNode->vSV);
			pNode->nSVnum = 0;	
			pNode->vSV = NULL;
		}

		pNode->vSV = (struct tagString**)calloc(nSVnum, sizeof(struct tagString*));
		if(pNode->vSV == NULL)
		{
			printf("Error: Network_InitNodeSV, cannot create node strings for network!\n");
			exit(EXIT_FAILURE);
		}
		pNode->nSVnum = nSVnum;

		pNode = pNode->pNext;
	}


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/* Network_LoadAnnotation:                                                 */
/* Annotate network nodes from a file.                                     */
/* ----------------------------------------------------------------------- */
int Network_LoadAnnotation(struct NETWORKPHYSICAL *pNet, char strFileName[])
{
	/* define */
	FILE *fpIn;
	char strLine[LONG_LINE_LENGTH];
	int nNodeId;
	char strNodeName[LONG_LINE_LENGTH];

	/* init check */
	if(pNet == NULL)
	{
		printf("Warning: Network_LoadAnnotation, empty network!\n");
		return PROC_SUCCESS;
	}

	/* open */
	fpIn = NULL;
	fpIn = fopen(strFileName, "r");
	if(fpIn == NULL)
	{
		printf("Warning: Network_LoadAnnotation, no annotation file!\n");
		return PROC_SUCCESS;
	}

	/* read file */
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%d %s", &nNodeId, strNodeName);
		if( (nNodeId>=0) && (nNodeId<pNet->nMaxNodeNum) )
		{
			if(pNet->vNodes[nNodeId] != NULL)
			{
				DeleteString(pNet->vNodes[nNodeId]->strName);
				pNet->vNodes[nNodeId]->strName = NULL;
				StringAddTail(&(pNet->vNodes[nNodeId]->strName), strNodeName);
			}
		}
	}

	/* close */
	fclose(fpIn);


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/* Network_WriteToSIF:                                                     */
/* Write a physical network to a SIF file.                                 */
/* ----------------------------------------------------------------------- */ 
int Network_WriteToSIF(struct NETWORKPHYSICAL *pNet, char strFileName[])
{
	/* file manipulation */
	FILE *fpOut;
	int ni;
	/* struct NETWORKNODE *pNode;
	struct NETWORKEDGE *pEdge;
	struct NODEEDGEINTERFACE *pNodeEdge; */
	
	/* init check */
	if(pNet == NULL)
	{
		printf("Warning: Network_WriteToSIF, empty network!\n");
		return PROC_SUCCESS;
	}

	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: Network_WriteToSIF, cannot open the file to write network!\n");
		exit(EXIT_FAILURE);
	}

	/* fprintf(fpOut, "NodeNum = %d\n", pNet->nRealNodeNum);
	pNode = pNet->pHeadNode;
	while(pNode != NULL)
	{
		fprintf(fpOut, "%d\t", pNode->nID);
		pNodeEdge = pNode->pEdgeList;
		while(pNodeEdge != NULL)
		{
			pEdge = pNet->vEdges[pNodeEdge->nEdgeID];

			if(pEdge->nNodeID1 == pNode->nID)
				fprintf(fpOut, "%d,", pEdge->nNodeID2);
			else
				fprintf(fpOut, "%d,", pEdge->nNodeID1);

			pNodeEdge = pNodeEdge->pNext;
		}

		fprintf(fpOut, "\n");

		pNode = pNode->pNext;
	}

	fprintf(fpOut, "EdgeNum = %d\n", pNet->nEdgeNum);
	*/

	for(ni=0; ni<pNet->nEdgeNum; ni++)
	{
		fprintf(fpOut, "%d\tpp\t%d\n", pNet->vEdges[ni]->nNodeID1, pNet->vEdges[ni]->nNodeID2); 
	}

	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/* Network_LoadValue_FromRefGene:                                          */
/* Load double values from a refgene database.                             */
/* ----------------------------------------------------------------------- */ 
int Network_LoadValue_FromRefGene(struct NETWORKPHYSICAL *pNet, int nRow, int nCol,
				struct tagRefGene **vRefGeneDatabase, int nRefGeneNum, 
				struct BYTEMATRIX *pHasValue, struct DOUBLEMATRIX *pSigValue,
				int nTakeAbsoluteValue, int nChooseMax)
{
	/* define */
	int ni;
	int nLocusID;
	double dValue,dTemp;

	/* init check */
	if( (pNet == NULL) || (vRefGeneDatabase == NULL) || (nRefGeneNum <= 0) )
	{
		printf("Warning: Network_LoadValue_FromRefGene, empty network or refgene database!\n");
		return PROC_SUCCESS;
	}

	/* assign value */
	for(ni=0; ni<nRefGeneNum; ni++)
	{
		if(vRefGeneDatabase[ni] == NULL)
			continue;

		nLocusID = vRefGeneDatabase[ni]->nGeneID;
		if(nLocusID < 0)
			continue;

		if(nLocusID >= pNet->nMaxNodeNum)
			continue;

		if(pNet->vNodes[nLocusID] == NULL)
			continue;

		if(pNet->vNodes[nLocusID]->nID != nLocusID)
		{
			printf("Error: Network_LoadValue_FromRefGene, gene id not match!\n");
			exit(EXIT_FAILURE);
		}


		if(pHasValue->pMatElement[ni] == 1)
		{
			dValue = pSigValue->pMatElement[ni];
			if(nTakeAbsoluteValue == 1)
				dValue = fabs(dValue);

			if(BMGETAT(pNet->vNodes[nLocusID]->pBV, nRow, nCol) == 0)
			{
				DMSETAT(pNet->vNodes[nLocusID]->pDV, nRow, nCol, dValue);
			}
			else
			{
				dTemp = DMGETAT(pNet->vNodes[nLocusID]->pDV, nRow, nCol);
				if(nChooseMax == 1)
				{
					if(dTemp < dValue)
						dTemp = dValue;
				}
				else if(nChooseMax == 2)
				{
					if(fabs(dTemp) < fabs(dValue))
						dTemp = dValue;
				}
				else if(nChooseMax == 3)
				{
					if(fabs(dTemp) > fabs(dValue))
						dTemp = dValue;
				}
				else
				{
					if(dTemp > dValue)
						dTemp = dValue;
				}
				DMSETAT(pNet->vNodes[nLocusID]->pDV, nRow, nCol, dTemp);
			}
			BMSETAT(pNet->vNodes[nLocusID]->pBV, nRow, nCol, 1);
		}
	}


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/* Network_LoadValue_FromRefGene_ForShortestPath:                          */
/* Load double values from a refgene database.                             */
/* ----------------------------------------------------------------------- */ 
int Network_LoadValue_FromRefGene_ForShortestPath(struct NETWORKPHYSICAL *pNet, 
				int nRow, int nCol,
				struct tagRefGene **vRefGeneDatabase, int nRefGeneNum, 
				struct BYTEMATRIX *pHasValue, struct DOUBLEMATRIX *pSigValue,
				int nTakeAbsoluteValue, int nChooseMax, double dMaxNodeValue)
{
	/* define */
	int ni;
	int nLocusID;
	double dValue,dTemp;
	double dLog10 = log(10.0);
	struct NETWORKNODE *pNode;

	/* init check */
	if( (pNet == NULL) || (vRefGeneDatabase == NULL) || (nRefGeneNum <= 0) )
	{
		printf("Warning: Network_LoadValue_FromRefGene, empty network or refgene database!\n");
		return PROC_SUCCESS;
	}

	/* assign value */
	for(ni=0; ni<nRefGeneNum; ni++)
	{
		if(vRefGeneDatabase[ni] == NULL)
			continue;

		nLocusID = vRefGeneDatabase[ni]->nGeneID;
		if(nLocusID < 0)
			continue;

		if(nLocusID >= pNet->nMaxNodeNum)
			continue;

		if(pNet->vNodes[nLocusID] == NULL)
			continue;

		if(pNet->vNodes[nLocusID]->nID != nLocusID)
		{
			printf("Error: Network_LoadValue_FromRefGene, gene id not match!\n");
			exit(EXIT_FAILURE);
		}


		if(pHasValue->pMatElement[ni] == 1)
		{
			dValue = pSigValue->pMatElement[ni];
			if(nTakeAbsoluteValue == 1)
				dValue = fabs(dValue);
			dValue = -log(2.0*normcdf(0.0, 1.0, dValue)-1.0)/dLog10;
			if(dValue > dMaxNodeValue)
				dValue = dMaxNodeValue;

			if(BMGETAT(pNet->vNodes[nLocusID]->pBV, nRow, nCol) == 0)
			{
				DMSETAT(pNet->vNodes[nLocusID]->pDV, nRow, nCol, dValue);
			}
			else
			{
				dTemp = DMGETAT(pNet->vNodes[nLocusID]->pDV, nRow, nCol);
				if(nChooseMax == 1)
				{
					if(dTemp < dValue)
						dTemp = dValue;
				}
				else
				{
					if(dTemp > dValue)
						dTemp = dValue;
				}
				DMSETAT(pNet->vNodes[nLocusID]->pDV, nRow, nCol, dTemp);
			}
			BMSETAT(pNet->vNodes[nLocusID]->pBV, nRow, nCol, 1);
		}
	}

	/* other nodes */
	pNode = pNet->pHeadNode;
	while(pNode != NULL)
	{
		if(BMGETAT(pNode->pBV, nRow, nCol) == 0)
		{
			DMSETAT(pNode->pDV, nRow, nCol, dMaxNodeValue);
		}
		pNode = pNode->pNext;
	}


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/* Network_DepthSearch_GetMaxScore:                                        */
/* Get maximun cumulative scores in the depth search.                      */
/* ----------------------------------------------------------------------- */ 
int Network_DepthSearch_GetMaxScore(struct NETWORKPHYSICAL *pNet, int nRow,
									int nCol, int nDepth)
{
	/* define */
	int ni;
	struct NETWORKNODE *pNode;
		
	struct DOUBLEMATRIX *pMaxScore;
	struct BYTEMATRIX *pMaxHas;
	struct tagString **vMaxPath;
	struct tagString *pPath;
	
	/* init check */
	if(pNet == NULL)
	{
		printf("Warning: Network_DepthSearch_GetMaxScore, empty network!\n");
		return PROC_SUCCESS;
	}

	if(nDepth <= 0)
	{
		printf("Warning: Network_DepthSearch_GetMaxScore, search depth <=0!\n");
		return PROC_SUCCESS;
	}

	/* process node by node */
	pNode = pNet->pHeadNode;
	while(pNode != NULL)
	{
		if(pNode->nNetIndex % 100 == 0)
		{
			printf("network depth search %d ...\n", pNode->nNetIndex);
		}

		/* prepare */
		pPath = NULL;
		pMaxScore = NULL;
		pMaxScore = CreateDoubleMatrix(1, nDepth);
		pMaxHas = NULL;
		pMaxHas = CreateByteMatrix(1, nDepth);
		vMaxPath = NULL;
		vMaxPath = (struct tagString **)calloc(nDepth, sizeof(struct tagString *));
		if( (pMaxScore == NULL) || (pMaxHas == NULL) || (vMaxPath == NULL) )
		{
			printf("Error: Network_DepthSearch_GetMaxScore, cannot create memory for depth search!\n");
			return PROC_SUCCESS;
		}

		/* search one node */
		Network_DepthSearch_GetMaxScore_Recursive(pNet, pNode, nRow, nCol, nDepth,
			0, pPath, 0.0, pMaxScore, pMaxHas, vMaxPath);

		/* assign maxscore */
		for(ni=0; ni<nDepth; ni++)
		{
			BMSETAT(pNode->pBV, nRow, (nCol+ni+1), pMaxHas->pMatElement[ni]);
			if(pMaxHas->pMatElement[ni] == 1)
			{
				DMSETAT(pNode->pDV, nRow, (nCol+ni+1), pMaxScore->pMatElement[ni]);
				StringAddTail(pNode->vSV+ni+1, vMaxPath[ni]->m_pString);
			}
		}

		/* release memory */
		DestroyDoubleMatrix(pMaxScore);
		DestroyByteMatrix(pMaxHas);
		for(ni=0; ni<nDepth; ni++)
		{
			DeleteString(vMaxPath[ni]);
			vMaxPath[ni] = NULL;
		}
		free(vMaxPath);

		/* get next node */
		pNode = pNode->pNext;
	}
	
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/* Network_DepthSearch_GetMaxScore_Recursive:                              */
/* Perform recursive depth search.                                         */
/* ----------------------------------------------------------------------- */ 
int Network_DepthSearch_GetMaxScore_Recursive(struct NETWORKPHYSICAL *pNet,
		struct NETWORKNODE *pNode, int nRow, int nCol, int nDepth,
		int nLayer, struct tagString *pPastPath, double dPastScore, 
		struct DOUBLEMATRIX *pMaxScore, struct BYTEMATRIX *pMaxHas,
		struct tagString **vMaxPath)
{
	/* define */
	double dScore;
	struct tagString *pPath = NULL;
	struct NODEEDGEINTERFACE *pInterface;
	int nEdgeID;
	int nNodeID;
	char strTemp[LINE_LENGTH];
	
	/* init check */
	if( (pNet == NULL) || (pNode == NULL) )
	{
		return PROC_SUCCESS;
	}

	if(nLayer > nDepth)
	{
		printf("Error: Network_DepthSearch_GetMaxScore, recursive overflow!\n");
		exit(EXIT_FAILURE);
	}

	if(BMGETAT(pNode->pBV, nRow, nCol) == 0)
	{
		return PROC_SUCCESS;
	}

	/* perform path check */
	dScore = DMGETAT(pNode->pDV, nRow, nCol);
	sprintf(strTemp, "-->%s(%d,%f)", pNode->strName->m_pString, pNode->nID, dScore);
	if(pPastPath == NULL)
	{
		dScore += dPastScore;
		StringAddTail(&pPath, strTemp);
	}
	else
	{
		if(strstr(pPastPath->m_pString, strTemp) != NULL)
		{
			return PROC_SUCCESS;
		}
		else
		{
			dScore += dPastScore;
			StringAddTail(&pPath, pPastPath->m_pString);
			StringAddTail(&pPath, strTemp);
		}
	}

	if(nLayer > 0)
	{
		if(pMaxHas->pMatElement[nLayer-1] == 0)
		{
			pMaxScore->pMatElement[nLayer-1] = dScore;
			DeleteString(vMaxPath[nLayer-1]);
			vMaxPath[nLayer-1] = NULL;
			StringAddTail(vMaxPath+(nLayer-1), pPath->m_pString);
			pMaxHas->pMatElement[nLayer-1] = 1;
		}
		else
		{
			if(dScore > pMaxScore->pMatElement[nLayer-1])
			{
				pMaxScore->pMatElement[nLayer-1] = dScore;
				DeleteString(vMaxPath[nLayer-1]);
				vMaxPath[nLayer-1] = NULL;
				StringAddTail(vMaxPath+(nLayer-1), pPath->m_pString);
			}
		}
	}

	/* if last node, only perform path check. */
	if(nLayer == nDepth)
	{
		DeleteString(pPath);
		return PROC_SUCCESS;
	}

	/* if intermediate node, perform recursive search */
	pInterface = pNode->pEdgeList;
	while(pInterface != NULL)
	{
		nEdgeID = pInterface->nEdgeID;
		if(pNet->vEdges[nEdgeID] != NULL)
		{
			if(pNet->vEdges[nEdgeID]->nNodeID1 != pNet->vEdges[nEdgeID]->nNodeID2)
			{
				if(pNode->nID == pNet->vEdges[nEdgeID]->nNodeID1)
				{
					nNodeID = pNet->vEdges[nEdgeID]->nNodeID2;
				}
				else
				{
					nNodeID = pNet->vEdges[nEdgeID]->nNodeID1;
				}

				if( (nNodeID>=0) && (nNodeID<pNet->nMaxNodeNum) )
				{
					Network_DepthSearch_GetMaxScore_Recursive(pNet, pNet->vNodes[nNodeID], 
						nRow, nCol, nDepth,	(nLayer+1), pPath, dScore, 
						pMaxScore, pMaxHas, vMaxPath);
				}
			}
		}

		pInterface = pInterface->pNext;
	}
	
	/* release memory */
	DeleteString(pPath);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/* NETWORKNODECREATE:                                                      */
/* Create an instance of NETWORKNODE.                                      */
/* ----------------------------------------------------------------------- */ 
struct NETWORKNODE *NETWORKNODECREATE()
{
	struct NETWORKNODE *pNode = NULL;

	pNode = (struct NETWORKNODE *)calloc(1, sizeof(struct NETWORKNODE));
	if(pNode == NULL)
	{
		printf("Warning: Cannot create network node!\n");
		return NULL;
	}

	pNode->nID = -1;
	pNode->nNetIndex = -1;
	pNode->strName = NULL;
	pNode->nTag = -1;
	pNode->pEdgeList = NULL;
	pNode->pNext = NULL;
	pNode->pBV = NULL;
	pNode->pDV = NULL;
	pNode->pIV = NULL;
	pNode->vSV = NULL;
	pNode->nSVnum = 0;

	return pNode;
}

/* ----------------------------------------------------------------------- */ 
/* NETWORKNODEDESTROY:                                                     */
/* Destroy an instance of NETWORKNODE.                                     */
/* ----------------------------------------------------------------------- */ 
void NETWORKNODEDESTROY(struct NETWORKNODE **ppNode)
{
	/* define */
	int ni;

	if(ppNode != NULL)
	{
		if(*ppNode != NULL)
		{
			/* clear name */
			DeleteString((*ppNode)->strName);
			(*ppNode)->strName = NULL;
			
			/* clear edge list */
			NODEEDGEINTERFACELIST_CLEARALL(&((*ppNode)->pEdgeList));

			/* clear values */
			if((*ppNode)->pBV != NULL)
			{
				DestroyByteMatrix((*ppNode)->pBV);
			}
			if((*ppNode)->pDV != NULL)
			{
				DestroyDoubleMatrix((*ppNode)->pDV);
			}
			if((*ppNode)->pIV != NULL)
			{
				DestroyIntMatrix((*ppNode)->pIV);
			}
			if((*ppNode)->vSV != NULL)
			{
				for(ni=0; ni<(*ppNode)->nSVnum; ni++)
				{
					DeleteString((*ppNode)->vSV[ni]);
					(*ppNode)->vSV[ni] = NULL;
				}
				free((*ppNode)->vSV);
			}
			
			/* free node */
			free(*ppNode);
			*ppNode = NULL;
		}
	}
}

/* ----------------------------------------------------------------------- */ 
/* NETWORKEDGECREATE:                                                      */
/* Create an instance of NETWORKEDGE.                                      */
/* ----------------------------------------------------------------------- */ 
struct NETWORKEDGE *NETWORKEDGECREATE()
{
	struct NETWORKEDGE *pEdge = NULL;

	pEdge = (struct NETWORKEDGE *)calloc(1, sizeof(struct NETWORKEDGE));
	if(pEdge == NULL)
	{
		printf("Warning: Cannot create network edge!\n");
		return NULL;
	}

	pEdge->nEdgeID = -1;
	pEdge->nNodeID1 = -1;
	pEdge->nNodeID2 = -1;
	pEdge->nDirection = 0;
	pEdge->nInteractionType = 0;

	return pEdge;
}

/* ----------------------------------------------------------------------- */ 
/* NETWORKEDGEDESTROY:                                                     */
/* Destroy an instance of NETWORKEDGE.                                     */
/* ----------------------------------------------------------------------- */ 
void NETWORKEDGEDESTROY(struct NETWORKEDGE **ppEdge)
{
	if(ppEdge != NULL)
	{
		if(*ppEdge != NULL)
		{
			free(*ppEdge);
			*ppEdge = NULL;
		}
	}
}

/* ----------------------------------------------------------------------- */ 
/* NODEEDGEINTERFACECREATE:                                                */
/* Create an instance of NODEEDGEINTERFACE.                                */
/* ----------------------------------------------------------------------- */ 
struct NODEEDGEINTERFACE *NODEEDGEINTERFACECREATE()
{
	struct NODEEDGEINTERFACE *pInterface = NULL;

	pInterface = (struct NODEEDGEINTERFACE *)calloc(1, sizeof(struct NODEEDGEINTERFACE));
	if(pInterface == NULL)
	{
		printf("Warning: Cannot create network interface!\n");
		return NULL;
	}

	pInterface->nEdgeID = -1;
	pInterface->pNext = NULL;

	return pInterface;
}

/* ----------------------------------------------------------------------- */ 
/* NODEEDGEINTERFACEDESTROY:                                               */
/* Destroy an instance of NODEEDGEINTERFACE.                               */
/* ----------------------------------------------------------------------- */ 
void NODEEDGEINTERFACEDESTROY(struct NODEEDGEINTERFACE **ppInterface)
{
	if(ppInterface != NULL)
	{
		if(*ppInterface != NULL)
		{
			(*ppInterface)->pNext = NULL;
			free(*ppInterface);
			*ppInterface = NULL;
		}
	}
}

/* ----------------------------------------------------------------------- */ 
/* NODEEDGEINTERFACELIST_CLEARALL:                                         */
/* Clear a list of NODEEDGEINTERFACE.                                      */
/* ----------------------------------------------------------------------- */ 
void NODEEDGEINTERFACELIST_CLEARALL(struct NODEEDGEINTERFACE **ppInterfaceList)
{
	/* define */
	struct NODEEDGEINTERFACE *pX = NULL;

	if(ppInterfaceList == NULL)
		return;

    while(*ppInterfaceList != NULL)
	{
		pX = *ppInterfaceList;
		*ppInterfaceList = pX->pNext;
		NODEEDGEINTERFACEDESTROY(&pX);
	}
}

/* ----------------------------------------------------------------------- */ 
/* NETWORKPHYSICALCREATE:                                                  */
/* Create an instance of NETWORKPHYSICAL.                                  */
/* ----------------------------------------------------------------------- */ 
struct NETWORKPHYSICAL *NETWORKPHYSICALCREATE()
{
	struct NETWORKPHYSICAL *pNet = NULL;

	pNet = (struct NETWORKPHYSICAL *)calloc(1, sizeof(struct NETWORKPHYSICAL));
	if(pNet == NULL)
	{
		printf("Warning: Cannot create a physical network!\n");
		return NULL;
	}

	pNet->nMaxNodeNum = 0;
	pNet->nRealNodeNum = 0;
	pNet->nEdgeNum = 0;
	pNet->pHeadNode = NULL;
	pNet->vNodes = NULL;
	pNet->vEdges = NULL;

	return pNet;
}

/* ----------------------------------------------------------------------- */ 
/* NETWORKPHYSICALDESTROY:                                                 */
/* Destroy an instance of NETWORKPHYSICAL.                                 */
/* ----------------------------------------------------------------------- */ 
void NETWORKPHYSICALDESTROY(struct NETWORKPHYSICAL **ppNetwork)
{
	/* define */
	struct NETWORKNODE *pNode = NULL;
	struct NETWORKEDGE *pEdge = NULL;
	int ni;

	if(ppNetwork == NULL)
		return;

    if(*ppNetwork == NULL)
		return;

	for(ni=0; ni<(*ppNetwork)->nEdgeNum; ni++)
	{
		NETWORKEDGEDESTROY((*ppNetwork)->vEdges+ni);
	}
	free((*ppNetwork)->vEdges);
	(*ppNetwork)->vEdges = NULL;

	for(ni=0; ni<(*ppNetwork)->nMaxNodeNum; ni++)
	{
		NETWORKNODEDESTROY((*ppNetwork)->vNodes+ni);
	}
	free((*ppNetwork)->vNodes);
	(*ppNetwork)->vNodes = NULL;

	free(*ppNetwork);
	*ppNetwork = NULL;
}

/* ----------------------------------------------------------------------- */ 
/* NETWORKSHORTESTPATHCREATE:                                              */
/* Create an object for performing shortest path search.                   */
/* ----------------------------------------------------------------------- */ 
struct NETWORKSHORTESTPATH *NETWORKSHORTESTPATHCREATE(int nNodeNum, int nTargetNum)
{
	/* define */
	struct NETWORKSHORTESTPATH *pNetPath = NULL;
	int ni;

	/* init check */
	if( (nNodeNum <= 0) || (nTargetNum <= 0) )
	{
		printf("Warning: NETWORKSHORTESTPATHCREATE, empty node list or target list!\n");
		return NULL;
	}

	/* create */
	pNetPath = (struct NETWORKSHORTESTPATH *)calloc(1, sizeof(struct NETWORKSHORTESTPATH));
	if(pNetPath == NULL)
	{
		printf("Error: NETWORKSHORTESTPATHCREATE, cannot create network shortest path object!\n");
		exit(EXIT_FAILURE);
	}

	/* init */
	pNetPath->nNodeNum = nNodeNum;
	pNetPath->pNodeID = CreateIntMatrix(1, nNodeNum);
	if(pNetPath->pNodeID == NULL)
	{
		printf("Error: NETWORKSHORTESTPATHCREATE, cannot create node IDs in network shortest path object!\n");
		exit(EXIT_FAILURE);
	}
	pNetPath->vNeighborNode = (struct INTMATRIX **)calloc(nNodeNum, sizeof(struct INTMATRIX *));
	if(pNetPath->vNeighborNode == NULL)
	{
		printf("Error: NETWORKSHORTESTPATHCREATE, cannot create neighbor node list in network shortest path object!\n");
		exit(EXIT_FAILURE);
	}
	pNetPath->vShortDist = (struct DOUBLEMATRIX **)calloc(nNodeNum, sizeof(struct DOUBLEMATRIX *));
	if(pNetPath->vShortDist == NULL)
	{
		printf("Error: NETWORKSHORTESTPATHCREATE, cannot create short dist list in network shortest path object!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nNodeNum; ni++)
	{
		pNetPath->vNeighborNode[ni] = CreateIntMatrix(1, nTargetNum);
		pNetPath->vShortDist[ni] = CreateDoubleMatrix(1, nTargetNum);
		if((pNetPath->vNeighborNode[ni] == NULL) || (pNetPath->vShortDist[ni] == NULL))
		{
			printf("Error: NETWORKSHORTESTPATHCREATE, cannot create shortest path elements in network shortest path object!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* return */
	return pNetPath;
}

/* ----------------------------------------------------------------------- */ 
/* NETWORKSHORTESTPATHDESTROY:                                             */
/* Destroy a shortest path object.                                         */
/* ----------------------------------------------------------------------- */ 
void NETWORKSHORTESTPATHDESTROY(struct NETWORKSHORTESTPATH **ppNetPath)
{
	/* define */
	int ni;

	if(ppNetPath == NULL)
		return;

    if(*ppNetPath == NULL)
		return;

	for(ni=0; ni<(*ppNetPath)->nNodeNum; ni++)
	{
		DestroyIntMatrix((*ppNetPath)->vNeighborNode[ni]);
		(*ppNetPath)->vNeighborNode[ni] = NULL;
		DestroyDoubleMatrix((*ppNetPath)->vShortDist[ni]);
		(*ppNetPath)->vShortDist[ni] = NULL;
	}
	free((*ppNetPath)->vNeighborNode);
	free((*ppNetPath)->vShortDist);

	DestroyIntMatrix((*ppNetPath)->pNodeID);

	free(*ppNetPath);
	*ppNetPath = NULL;
}

/* ----------------------------------------------------------------------- */ 
/* NETWORKNODE_ADDEDGE:                                                    */
/* Add an edge to a network node.                                          */
/* ----------------------------------------------------------------------- */ 
int NETWORKNODE_ADDEDGE(struct NETWORKNODE **ppNode, int nNodeID, int nEdgeID)
{
	/* define */
	struct NODEEDGEINTERFACE *pNewInterface;

	/* process */
	if(ppNode == NULL)
	{
		printf("Warning: there is no valid network handle!\n");
		return PROC_FAILURE;
	}

	if(*ppNode == NULL)
	{
		*ppNode = NETWORKNODECREATE();
		if(*ppNode == NULL)
		{
			printf("Error: there is no valid network handle!\n");
			exit(EXIT_FAILURE);
		}
		(*ppNode)->nID = nNodeID;
	}

	if((*ppNode)->nID != nNodeID)
	{
		printf("Error: network node IDs do not match!\n");
		exit(EXIT_FAILURE);
	}

	pNewInterface = NULL;
	pNewInterface = NODEEDGEINTERFACECREATE();
	if(pNewInterface == NULL)
	{
		printf("Error: cannot create network node edge interface!\n");
		exit(EXIT_FAILURE);
	}
	pNewInterface->nEdgeID = nEdgeID;
	if(NODEEDGEINTERFACELIST_ADDHEAD(&((*ppNode)->pEdgeList), pNewInterface) == PROC_FAILURE)
	{
		NODEEDGEINTERFACEDESTROY(&pNewInterface);
	}

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/* NODEEDGEINTERFACELIST_ADDHEAD:                                          */
/* Add a node-edge interface object to a node-edge interface list.         */
/* ----------------------------------------------------------------------- */ 
int NODEEDGEINTERFACELIST_ADDHEAD(struct NODEEDGEINTERFACE **ppInterfaceList, 
		struct NODEEDGEINTERFACE *pInterface)
{
	/* add */
	if(ppInterfaceList == NULL)
	{
		printf("Warning: there is no valid interface list handle!\n");
		return PROC_FAILURE;
	}

	if(pInterface == NULL)
	{
		printf("Warning: the node-edge interface object is empty!\n");
		return PROC_SUCCESS;
	}

	pInterface->pNext = *ppInterfaceList;
	*ppInterfaceList = pInterface;

	/* return */
	return PROC_SUCCESS;
}
