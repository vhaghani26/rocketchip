#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
#include "limits.h"
#include "time.h"

#include "StringLib.h"
#include "MatrixLib.h"
#include "RandomLib.h"
#include "MathLib.h"
#include "MotifLib.h"
#include "SequenceLib.h"
#include "GenomeLib.h"
#include "MicroarrayLib.h"
#include "AffyLib.h"
#include "Alignment.h"
#include "WorkLib.h"

int GenomeGetConservedSeg_Main(int argv, char **argc);
int GenomeGetConservedSeg_Write(struct SEQLINKMAP **pSegList, int nMinSegLen, 
							 int nMaxGap, int nMinTotLen, FILE *fpOut, int *pId);


int main(int argv, char **argc)
{
	int nLen;
	int nseed;

	/* init rand */
	srand( (unsigned)time( NULL ) );
	if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
	{
		nseed = (int)(rand()*1000/RAND_MAX);
	}
	else
	{
		nseed = rand()%1000;
	}
	rand_u_init(nseed);


	/* ---- */
	/* menu */
	/* ---- */
	GenomeGetConservedSeg_Main(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int GenomeGetConservedSeg_Main(int argv, char **argc)
{
	/* define */
	char strOutFile[LINE_LENGTH];
	char strFileName[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	char strGenomePath[LINE_LENGTH];
	char strCSPath[LINE_LENGTH];
	char strCSFile[LINE_LENGTH];
	char strCDSPath[LINE_LENGTH];
	char strCDSFile[LINE_LENGTH];
	int nId = 0;
	int nC = 132;
	int nInitSegLen = 50;
	int nMinSegLen = 100;
	int nMaxGap = 50;
	int nMinTotLen = 200;
	struct INTMATRIX *pChrLen;
	int nChrNum = 21;
	char strChrName[LINE_LENGTH];
	int nChr;
	FILE *fpIn;
	FILE *fpCDS;
	FILE *fpOut;
	unsigned char bScore,bCDS;

	int nNewSeg = 0;
	struct SEQLINKMAP *pSeg = NULL;
	struct SEQLINKMAP *pSegList = NULL;
	struct SEQLINKMAP *pLinkSegList = NULL;
	struct SEQLINKMAP *pCombinedSeg = NULL;
	struct SEQLINKMAP *pPrev = NULL;
	int ni,numread;


	/* init */
	if(argv != 2)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	strcpy(strOutFile, argc[1]);

	strcpy(strSpecies, "mouse");
	strcpy(strGenomePath, "/Volumes/Server HD 2/genomes/mouse/mm7/");
	strcpy(strCSPath, "/Volumes/Server HD 2/genomes/mouse/mm7/conservation/phastcons/");
	strcpy(strCDSPath, "/Volumes/Server HD 2/genomes/mouse/mm7/cds/");
	
	/* load chromosome length */
	sprintf(strFileName, "%schrlen.txt", strGenomePath);
	pChrLen = NULL;
	pChrLen = IMLOAD(strFileName);
	if(pChrLen == NULL)
	{
		printf("Error: cannot load chromosome length!\n");
		exit(EXIT_FAILURE);
	}

	/* open export file */
	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	/* process one by one */
	for(nChr=1; nChr<=nChrNum; nChr++)
	{
		/* get chr name */
		Genome_Index_To_ChromosomeName(strChrName, strSpecies, nChr);
			
		/* get conserved score */
		sprintf(strCSFile, "%s%s.cs", strCSPath, strChrName);
		fpIn = NULL;
		fpIn = fopen(strCSFile, "rb");
		if(fpIn == NULL)
		{
			printf("Error: cannot open conservation score file!\n");
			exit(EXIT_FAILURE);
		}

		sprintf(strCDSFile, "%s%s.cds", strCDSPath, strChrName);
		fpCDS = NULL;
		fpCDS = fopen(strCDSFile, "rb");
		if(fpCDS == NULL)
		{
			printf("Error: cannot open cds file!\n");
			exit(EXIT_FAILURE);
		}

		/* read */
		for(ni=0; ni<pChrLen->pMatElement[nChr-1]; ni++)
		{
			numread = fread(&bScore, sizeof(unsigned char), 1, fpIn);
			if(numread != 1)
			{
				printf("Error: loading error!\n");
				exit(EXIT_FAILURE);
			}

			numread = fread(&bCDS, sizeof(unsigned char), 1, fpCDS);
			if(numread != 1)
			{
				printf("Error: loading error!\n");
				exit(EXIT_FAILURE);
			}
	
			if( ((int)bScore > nC) && (bCDS != 1) )
			{
				if(nNewSeg == 0)
				{
					pSeg = NULL;
					pSeg = SEQLINKMAPCREATE();
					pSeg->nGStart = ni;
					pSeg->nGEnd = ni;
					strcpy(pSeg->strGenome, strChrName);
					pSeg->chGStrand = '+';
					nNewSeg = 1;
				}
				else
				{
					pSeg->nGEnd = ni;
				}
			}
			else
			{
				if(nNewSeg == 1)
				{
					if((pSeg->nGEnd-pSeg->nGStart+1) >= nInitSegLen)
					{
						if(pSegList == NULL)
						{
							pSegList = pSeg;
							pPrev = pSeg;
						}
						else
						{
							pPrev->pNext = pSeg;
							pPrev = pSeg;
						}
					}
					else
					{
						SEQLINKMAPDESTROY(pSeg);
					}
					nNewSeg = 0;
				}
			}
		}
		
		if(nNewSeg == 1)
		{
			if((pSeg->nGEnd-pSeg->nGStart+1) >= nInitSegLen)
			{
				if(pSegList == NULL)
				{
					pSegList = pSeg;
					pPrev = pSeg;
				}
				else
				{
					pPrev->pNext = pSeg;
					pPrev = pSeg;
				}
			}
			else
			{
				SEQLINKMAPDESTROY(pSeg);
			}
			nNewSeg = 0;
		}

		fclose(fpIn);
		fclose(fpCDS);

		/* write */
		GenomeGetConservedSeg_Write(&pSegList, nMinSegLen, nMaxGap, nMinTotLen,
							fpOut, &nId);
	}

	/* close files */		
	fclose(fpOut);

	/* release memory */
	DestroyIntMatrix(pChrLen);

	/* return */
	return PROC_SUCCESS;
}

int GenomeGetConservedSeg_Write(struct SEQLINKMAP **pSegList, int nMinSegLen, 
							 int nMaxGap, int nMinTotLen, FILE *fpOut, int *pId)
{
	/* define */
	struct SEQLINKMAP *pLinkSegList = NULL;
	struct SEQLINKMAP *pSeg = NULL;
	struct SEQLINKMAP *pCombinedSeg = NULL;
	struct SEQLINKMAP *pPrev = NULL;
	int nNewSeg = 0;
	int nAddCombined;
	int nMaxSegLen = 0;
	int nSegLen = 0;

	/* link segments */
	pCombinedSeg = NULL;
	while(*pSegList != NULL)
	{
		pSeg = *pSegList;
		*pSegList = pSeg->pNext;
		pSeg->pNext = NULL;
		
		if(pCombinedSeg == NULL)
		{
			pCombinedSeg = pSeg;
			nMaxSegLen = pSeg->nGEnd-pSeg->nGStart+1;
		}
		else
		{
			if( (pSeg->nGStart-pCombinedSeg->nGEnd) <= nMaxGap )
			{
					pCombinedSeg->nGEnd = pSeg->nGEnd;
					nSegLen = pSeg->nGEnd-pSeg->nGStart+1;
					if(nSegLen > nMaxSegLen)
						nMaxSegLen = nSegLen;

					nAddCombined = 0;
					SEQLINKMAPDESTROY(pSeg);
			}
			else
			{
				nAddCombined = 1;
			}
				
			if(nAddCombined == 1)
			{
				if( ((pCombinedSeg->nGEnd-pCombinedSeg->nGStart+1) < nMinTotLen)
					|| (nMaxSegLen < nMinSegLen) )
				{
					SEQLINKMAPDESTROY(pCombinedSeg);
					pCombinedSeg = NULL;
					nMaxSegLen = 0;
				}
				else
				{
					if(pLinkSegList == NULL)
					{
						pLinkSegList = pCombinedSeg;
						pPrev = pCombinedSeg;
					}
					else
					{
						pPrev->pNext = pCombinedSeg;
						pPrev = pCombinedSeg;
					}
				}

				pCombinedSeg = pSeg;
				nMaxSegLen = pSeg->nGEnd-pSeg->nGStart+1;
			}
		}
	}

	if(pCombinedSeg != NULL)
	{
		if( ((pCombinedSeg->nGEnd-pCombinedSeg->nGStart+1) < nMinTotLen)
			|| (nMaxSegLen < nMinSegLen) )
		{
			SEQLINKMAPDESTROY(pCombinedSeg);
			pCombinedSeg = NULL;
			nMaxSegLen = 0;
		}
		else
		{
			if(pLinkSegList == NULL)
			{
				pLinkSegList = pCombinedSeg;
				pPrev = pCombinedSeg;
			}
			else
			{
				pPrev->pNext = pCombinedSeg;
				pPrev = pCombinedSeg;
			}
		}
	}


	/* write */
	while(pLinkSegList != NULL)
	{
		pSeg = pLinkSegList;
		pLinkSegList = pSeg->pNext;
		pSeg->pNext = NULL;

		fprintf(fpOut, "%d\t%s\t%d\t%d\t+\n", *pId, pSeg->strGenome, pSeg->nGStart, pSeg->nGEnd);
		(*pId) += 1;

		SEQLINKMAPDESTROY(pSeg);
	}

	/* return */
	return PROC_SUCCESS;
}

