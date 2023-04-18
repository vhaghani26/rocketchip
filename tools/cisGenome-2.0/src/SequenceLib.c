/* ----------------------------------------------------------------------- */
/*  SequenceLib.c : implementation of the sequence library                 */
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

/* ----------------------------------------------------------------------------- */ 
/*                       struct tagSequence *SequenceCreate()                    */
/* This function is used for creating a sequence structure.                      */
/* ----------------------------------------------------------------------------- */ 
struct tagSequence *SequenceCreate()
{
	struct tagSequence *pNewSequence;

	/* Init */
	pNewSequence = NULL;

	/* new */
	pNewSequence = (struct tagSequence *)malloc(sizeof(struct tagSequence));
	if( pNewSequence == NULL )
		return NULL;
	pNewSequence->m_nIndex = -1;
	pNewSequence->m_nLength = 0;
	pNewSequence->m_pSequence = NULL;
	strcpy(pNewSequence->m_strAlias, "");
	pNewSequence->m_nStart = 0;
	pNewSequence->m_nEnd = 0;
	pNewSequence->m_pNext = NULL;

	/* return */
	return pNewSequence;
}

/* ----------------------------------------------------------------------------- */ 
/*                           void SequenceDelete()                               */
/* This function is used for deleting a sequence.                                */
/* ----------------------------------------------------------------------------- */ 
void SequenceDelete(struct tagSequence *pDelSequence)
{
	if( pDelSequence != NULL )
	{
		pDelSequence->m_pNext = NULL;
		DeleteString(pDelSequence->m_pSequence);
		free(pDelSequence);
	}	
}

/* ----------------------------------------------------------------------------- */ 
/*                            int  SequenceAddTail()                             */
/* This function is used for adding a segment to the tail of a sequence.         */
/* If the process is successful, it will return 1;                               */
/* else it will return 0.                                                        */            
/* ----------------------------------------------------------------------------- */ 
int SequenceAddTail(struct tagSequence *pSequence, char strLine[])
{
	/* nLen is used for recording the length of a string */
	int nLen;
	/* pString is used for manipulating the string */
	char *pString,*pTempString;

	/* check parameter */
	if(pSequence == NULL)
		return PROC_FAILURE;
	nLen = (int)strlen(strLine);
	if(nLen == 0)
		return PROC_SUCCESS;

	/* add */
	if(pSequence->m_pSequence == NULL)
	{
		pSequence->m_pSequence = CreateString(nLen);
		if(pSequence->m_pSequence == NULL)
		{
			printf("Error in SequenceAddTail: cannot create string file...\n");
			return PROC_FAILURE;
		}
		pString = pSequence->m_pSequence->m_pString;
		strcat(pString, strLine);
		pSequence->m_nLength = nLen;
		pSequence->m_nStart = 1;
		pSequence->m_nEnd = nLen;
	}
	else
	{
		pSequence->m_pSequence->m_nLength += nLen;
		pString = pSequence->m_pSequence->m_pString;
		pTempString = (char *)malloc(sizeof(char)*(pSequence->m_pSequence->m_nLength+1));
		if(pTempString == NULL)
			return PROC_FAILURE;
		strcpy(pTempString, pString);
		free(pString);
		strcat(pTempString, strLine);
		pSequence->m_pSequence->m_pString = pTempString;
		pSequence->m_nLength += nLen;
		pSequence->m_nEnd += nLen;
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------------- */ 
/*                           void SequenceListClear()                            */
/* This function is used for clearing sequence list.                             */
/* ----------------------------------------------------------------------------- */ 
void SequenceListClear(struct tagSequence **pSequenceList)
{
	struct tagSequence *pHeadNode;

	/* clear sequence list. */
	pHeadNode = *pSequenceList;
	while( pHeadNode != NULL )
	{
		*pSequenceList = pHeadNode->m_pNext;
		SequenceDelete(pHeadNode);
		pHeadNode = *pSequenceList;
	}
}


/* ----------------------------------------------------------------------------- */ 
/*                          int SequenceListGetSize()                            */
/* This function is used for getting the size of a sequence list.                */
/* ----------------------------------------------------------------------------- */ 
int SequenceListGetSize(struct tagSequence *pSequenceList)
{
	int nCount;
	struct tagSequence *pSeq;

	/* init */
	nCount = 0;

	/* counting */
	pSeq = pSequenceList;
	while( pSeq != NULL )
	{
		nCount++;
		pSeq = pSeq->m_pNext;
	}

	/* return */
	return nCount;
}

/* ----------------------------------------------------------------------------- */ 
/*                     struct tagSequence *SequenceListGetAt()                   */
/* This function is used for getting a sequence from a sequence list.            */
/* ----------------------------------------------------------------------------- */ 
struct tagSequence *SequenceListGetAt(struct tagSequence *pSequenceList, int nIndex)
{
	int nCount;
	struct tagSequence *pSeq,*pTarget;

	/* init */
	nCount = SequenceListGetSize(pSequenceList);
	if((nIndex < 0) || (nIndex >= nCount))
	{
		printf("Error in struct tagSequence *SequenceListGetAt(): array index out of bound...\n");
		exit(EXIT_FAILURE);
	}

	pTarget = NULL;
	nCount = 0;

	/* searching */
	pSeq = pSequenceList;
	while( pSeq != NULL )
	{
		if(nCount == nIndex)
		{
			pTarget = pSeq;
			break;
		}
		nCount++;
		pSeq = pSeq->m_pNext;
	}

	/* return */
	return pTarget;
}

/* ----------------------------------------------------------------------------- */ 
/*                   struct tagSequence *SequenceListGetTail()                   */
/* This function is used for getting the tail sequence from a sequence list.     */
/* ----------------------------------------------------------------------------- */ 
struct tagSequence *SequenceListGetTail(struct tagSequence *pSequenceList)
{
	struct tagSequence *pSeq,*pTarget;

	/* searching */
	pTarget = NULL;
	pSeq = pSequenceList;
	while( pSeq != NULL )
	{
		pTarget = pSeq;		
		pSeq = pSeq->m_pNext;
	}

	/* return */
	return pTarget;
}


/* ----------------------------------------------------------------------------- */ 
/*                      int LoadFullSequenceList()                               */
/* This function loads all sequences from the input FASTA file.                  */
/* If the process is successful, it will return the count of sequences;          */
/* else it will return 0;                                                        */            
/* ----------------------------------------------------------------------------- */
int LoadFullSequenceList(char strInFilePath[], struct tagSequence **pSeqList)
{
	/* define */
	FILE *fpInput;
	char strLine[LONG_LINE_LENGTH];
	struct tagSequence *pNewSeq;
	struct tagSequence *pTailSeq;
	int nCount;
	
	/* get tail seq */
	pNewSeq = NULL;
	nCount = 0;
	if(*pSeqList == NULL)
		pTailSeq = NULL;
	else
	{
		pTailSeq = SequenceListGetTail(*pSeqList);
	}

	/* load file */
	fpInput = NULL;
	fpInput = fopen(strInFilePath, "rt");
	if(fpInput == NULL)
	{
        printf("Error in LoadConserveSequenceList: Cannot open input file...\n");
		return 0;
	}
	
	/* read */
	while(fgets(strLine, LONG_LINE_LENGTH, fpInput) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strLine[0] == '>')
		{
			if(pNewSeq != NULL)
			{
				if(pTailSeq == NULL)
				{
					*pSeqList = pNewSeq;
					pTailSeq = pNewSeq;
				}
				else
				{
					pTailSeq->m_pNext = pNewSeq;
					pTailSeq = pNewSeq;
				}

				nCount++;
			}

			pNewSeq = NULL;
			pNewSeq = SequenceCreate();
			if(pNewSeq == NULL)
			{
				fclose(fpInput);
				printf("Error in LoadConserveSequenceList: Cannot create a new sequence node...\n");
				return 0;
			}

			pNewSeq->m_nIndex = nCount;
			strcpy(pNewSeq->m_strAlias, (strLine+1));
		}
		else
		{
			if(SequenceAddTail(pNewSeq, strLine) == PROC_FAILURE)
			{
				fclose(fpInput);
				printf("Error in LoadConserveSequenceList: cannot add string to sequence...\n");
				return 0;
			}			
		}
	}

	if(pNewSeq != NULL)
	{
		if(pTailSeq == NULL)
		{
			*pSeqList = pNewSeq;
			pTailSeq = pNewSeq;
		}
		else
		{
			pTailSeq->m_pNext = pNewSeq;
			pTailSeq = pNewSeq;
		}
		nCount++;
	}

	/* close files */
	fclose(fpInput);
	
	/* return */
	return nCount;
}

/* ----------------------------------------------------------------------------- */ 
/*                      FILE *SequenceWriteToFasta()                             */
/* This function is used for writing a sequences into a fasta file.              */
/* If the process is successful, it will return the end pointer of the file;     */
/* else it will return NULL;                                                     */
/* ----------------------------------------------------------------------------- */
FILE *SequenceWriteToFasta(struct tagSequence *pSequence, FILE *fpOut)
{
	/* nLinePos and pSeqPos are used for accessing the sequence. */
	int nLinePos;
	char *pSeqPos;

	/* check parameter */
	if(pSequence == NULL)
		return fpOut;
	if(fpOut == NULL)
		return NULL;

	/* write */
	if(pSequence->m_pSequence != NULL)
	{
		/* write head */
		fprintf(fpOut, ">%d\n", pSequence->m_nIndex);
		
		/* write sequence */
		nLinePos = 0;
		pSeqPos = pSequence->m_pSequence->m_pString;
		while(*pSeqPos != '\0')
		{
			fprintf(fpOut, "%c", *pSeqPos);
			nLinePos++;
			pSeqPos++;

			if(nLinePos == FASTA_LINE_LEN)
			{
				fprintf(fpOut, "\n");
				nLinePos = 0;
			}
		}
		if(nLinePos != 0)
		{
			fprintf(fpOut, "\n");
		}
	}

	/* return */
	return fpOut;
}

/* ----------------------------------------------------------------------------- */ 
/*                    SequenceWriteToFasta_ByStrand()                            */
/* This function is used for writing a sequences into a fasta file, according    */
/* to the strand specified.                                                      */
/* If the process is successful, it will return the end pointer of the file;     */
/* else it will return NULL;                                                     */
/* if nAliasType = 0: output the m_nIndex of the sequence as the header after >  */
/*                 1: output the strAlias as the header after >                  */
/* ----------------------------------------------------------------------------- */
FILE *SequenceWriteToFasta_ByStrand(struct tagSequence *pSequence, FILE *fpOut, 
									char chStrand, int nAliasType)
{
	/* nLinePos and pSeqPos are used for accessing the sequence. */
	int nLinePos;
	char *pSeqPos;
	int ni;
	char chBase;

	/* check parameter */
	if(pSequence == NULL)
		return fpOut;
	if(fpOut == NULL)
		return NULL;

	/* write */
	if(pSequence->m_pSequence != NULL)
	{
		/* write head */
		if(nAliasType == 1)
		{
			fprintf(fpOut, ">%s\n", pSequence->m_strAlias);
		}
		else
		{
			fprintf(fpOut, ">%d\n", pSequence->m_nIndex);
		}

		/* write sequence */
		if(chStrand == '-')
		{
			nLinePos = 0;
			pSeqPos = pSequence->m_pSequence->m_pString;
			for(ni=(pSequence->m_nLength-1); ni>=0; ni--)
			{
				switch(pSeqPos[ni])
				{
					case 'a': chBase = 't';
						break;
					case 'A': chBase = 'T';
						break;
					case 'c': chBase = 'g';
						break;
					case 'C': chBase = 'G';
						break;
					case 'g': chBase = 'c';
						break;
					case 'G': chBase = 'C';
						break;
					case 't': chBase = 'a';
						break;
					case 'T': chBase = 'A';
						break;
					case 'n': chBase = 'n';
						break;
					default: chBase = 'N';
				}

				fprintf(fpOut, "%c", chBase);
				nLinePos++;
				
				if(nLinePos == FASTA_LINE_LEN)
				{
					fprintf(fpOut, "\n");
					nLinePos = 0;
				}
			}
			if(nLinePos != 0)
			{
				fprintf(fpOut, "\n");
			}
		}
		else
		{
			nLinePos = 0;
			pSeqPos = pSequence->m_pSequence->m_pString;
			while(*pSeqPos != '\0')
			{
				fprintf(fpOut, "%c", *pSeqPos);
				nLinePos++;
				pSeqPos++;

				if(nLinePos == FASTA_LINE_LEN)
				{
					fprintf(fpOut, "\n");
					nLinePos = 0;
				}
			}
			if(nLinePos != 0)
			{
				fprintf(fpOut, "\n");
			}
		}
		
	}

	/* return */
	return fpOut;
}

/* ----------------------------------------------------------------------------- */ 
/*                      int FastaSequenceMask_Main()                             */
/* This function is used for masking specified regions from a FASTA file.        */
/* The Maskfile should have the following format:                                */
/*   Col1: numerical seqid                                                       */
/*   Col2: start                                                                 */
/*   Col3: end                                                                   */
/*  If nMaskType == 0, masked sequence will in little letters a,c,g,t.           */
/*  If nMaskType == 1, masked sequence will be converted to N.                   */
/* ----------------------------------------------------------------------------- */
int FastaSequenceMask_Main(char strSeqFile[], char strMaskFile[],
			int nMaskType, char strOutFile[])
{
	/* define */
	FILE *fpMask;
	FILE *fpOut;
	struct tagSequence *pSeqList = NULL;
	int nSeqCount = 0;
	struct tagSequence **vSeqList = NULL;
	struct tagSequence *pSeq;
	char strLine[LONG_LINE_LENGTH];
	int nSeqId,nStart,nEnd,nTemp;
	int ni;
	char *vB;

	/* load sequence */
	nSeqCount = LoadFullSequenceList(strSeqFile, &pSeqList);
	if(nSeqCount <= 0)
	{
		printf("Warning: No sequences loaded!\n");
		SequenceListClear(&pSeqList);
		return PROC_SUCCESS;
	}

	/* reorganize sequence */
	vSeqList = NULL;
	vSeqList = (struct tagSequence **)calloc(nSeqCount, sizeof(struct tagSequence *));
	if(vSeqList == NULL)
	{
		printf("Error: cannot allocate memory for sequence list!\n");
		SequenceListClear(&pSeqList);
		exit(EXIT_FAILURE);
	}
	
	ni = 0;
	pSeq = pSeqList;
	while(pSeq != NULL)
	{
		vSeqList[ni] = pSeq;
		ni++;
		pSeq = pSeq->m_pNext;
	}
	if(ni != nSeqCount)
	{
		printf("Error: sequence number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* masking */
	fpMask = NULL;
	fpMask = fopen(strMaskFile, "r");
	if(fpMask == NULL)
	{
		printf("Error: cannot open mask file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpMask) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		sscanf(strLine, "%d %d %d", &nSeqId, &nStart, &nEnd);
		if( (nSeqId < 0) || (nSeqId >= nSeqCount) )
		{
			printf("Warning: mask out of range: %s!\n", strLine);
			continue;
		}
		vB = vSeqList[nSeqId]->m_pSequence->m_pString;
		if( (nStart < 0) || (nEnd < 0) || (nStart >= vSeqList[nSeqId]->m_nLength) ||
			(nEnd >= vSeqList[nSeqId]->m_nLength) )
		{
			printf("Warning: mask out of range: %s!\n", strLine);
			continue;
		}
		if(nStart > nEnd)
		{
			nTemp = nStart;
			nStart = nEnd;
			nEnd = nTemp;
		}

		if(nMaskType == 0)
		{
			for(ni=nStart; ni<=nEnd; ni++)
			{
				if(vB[ni] == 'A')
				{
					vB[ni] = 'a';
				}
				else if(vB[ni] == 'C')
				{
					vB[ni] = 'c';
				}
				else if(vB[ni] == 'G')
				{
					vB[ni] = 'g';
				}
				else if(vB[ni] == 'T')
				{
					vB[ni] = 't';
				}
				else if(vB[ni] == 'a')
				{
				}
				else if(vB[ni] == 'c')
				{	
				}
				else if(vB[ni] == 'g')
				{	
				}
				else if(vB[ni] == 't')
				{	
				}
				else
				{
					vB[ni] = 'N';
				}
			}
		}
		else
		{
			for(ni=nStart; ni<=nEnd; ni++)
			{
				vB[ni] = 'N';
			}
		}
	}

	fclose(fpMask);

	/* export */
	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nSeqCount; ni++)
	{
		SequenceWriteToFasta_ByStrand(vSeqList[ni], fpOut, '+', 1);
	}
	
	fclose(fpOut);

	/* clear memory */
	SequenceListClear(&pSeqList);
	free(vSeqList);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------------- */ 
/*                     int FastaSequenceExtract_Main()                           */
/* This function is used for extracting sequences for specified regions from a   */
/* FASTA file. The Maskfile should have the following format:                    */
/*   Col1: numerical seqid                                                       */
/*   Col2: start                                                                 */
/*   Col3: end                                                                   */
/*   Col4: strand                                                                */
/* ----------------------------------------------------------------------------- */
int FastaSequenceExtract_Main(char strFASTAPath[], char strTargetFile[], char strOutputFile[], 
					  int nUp, int nDown, int nStrandType)
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	struct tagSequence *pSeqList = NULL;
	int nSeqCount = 0;
	struct tagSequence **vSeqList = NULL;
	struct tagSequence *pSeq;
	char strLine[LONG_LINE_LENGTH];
	int nSeqId,nStart,nEnd;
	char chStrand;
	int ni,nj;
	char *vB;

	/* load sequence */
	nSeqCount = LoadFullSequenceList(strFASTAPath, &pSeqList);
	if(nSeqCount <= 0)
	{
		printf("Warning: No sequences loaded!\n");
		SequenceListClear(&pSeqList);
		return PROC_SUCCESS;
	}

	/* reorganize sequence */
	vSeqList = NULL;
	vSeqList = (struct tagSequence **)calloc(nSeqCount, sizeof(struct tagSequence *));
	if(vSeqList == NULL)
	{
		printf("Error: cannot allocate memory for sequence list!\n");
		SequenceListClear(&pSeqList);
		exit(EXIT_FAILURE);
	}
	
	ni = 0;
	pSeq = pSeqList;
	while(pSeq != NULL)
	{
		vSeqList[ni] = pSeq;
		ni++;
		pSeq = pSeq->m_pNext;
	}
	if(ni != nSeqCount)
	{
		printf("Error: sequence number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* extracting sequences */
	fpIn = NULL;
	fpIn = fopen(strTargetFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: cannot open input file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutputFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%d %d %d %c", &nSeqId, &nStart, &nEnd, &chStrand);
		if( (nSeqId < 0) || (nSeqId >= nSeqCount) )
		{
			printf("Warning: target out of range: %s!\n", strLine);
			continue;
		}

		if(nStrandType == 0)
		{
			nStart -= nUp;
			nEnd += nDown;
		}
		else
		{
			if(chStrand == '+')
			{
				nStart -= nUp;
				nEnd += nDown;
			}
			else
			{
				nEnd += nUp;
				nStart -= nDown;
			}
		}

		if( (nStart > nEnd) || (nStart < 0) || (nEnd >= vSeqList[nSeqId]->m_nLength) )
		{
			printf("Warning: target out of range: %s!\n", strLine);
			continue;
		}

		fprintf(fpOut, "%d\t%d\t%d\t%c\t", nSeqId, nStart, nEnd, chStrand, vSeqList[nSeqId]->m_strAlias);
		
		vB = vSeqList[nSeqId]->m_pSequence->m_pString;
		if( (nStrandType == 0) || (chStrand == '+') )
		{
			nj = 0;
			for(ni=nStart; ni<=nEnd; ni++)
			{
				fprintf(fpOut, "%c", vB[ni]);
				nj++;
				if( nj == FASTA_LINE_LEN )
				{
					fprintf(fpOut, "\n");
					nj = 0;
				}
			}

			if(nj != 0)
			{
				fprintf(fpOut, "\n");
			}
		}
		else
		{
			nj = 0;
			for(ni=nEnd; ni>=nStart; ni--)
			{
				if(vB[ni] == 'A')
				{
					fprintf(fpOut, "T");
				}
				else if(vB[ni] == 'C')
				{
					fprintf(fpOut, "G");
				}
				else if(vB[ni] == 'G')
				{
					fprintf(fpOut, "C");
				}
				else if(vB[ni] == 'T')
				{
					fprintf(fpOut, "A");
				}
				else if(vB[ni] == 'a')
				{
					fprintf(fpOut, "t");
				}
				else if(vB[ni] == 'c')
				{	
					fprintf(fpOut, "g");
				}
				else if(vB[ni] == 'g')
				{	
					fprintf(fpOut, "c");
				}
				else if(vB[ni] == 't')
				{	
					fprintf(fpOut, "a");
				}
				else
				{
					vB[ni] = 'N';
				}

				nj++;
				if( nj == FASTA_LINE_LEN )
				{
					fprintf(fpOut, "\n");
					nj = 0;
				}
			}

			if(nj != 0)
			{
				fprintf(fpOut, "\n");
			}
		}
	}

	fclose(fpIn);
	fclose(fpOut);

	/* clear memory */
	SequenceListClear(&pSeqList);
	free(vSeqList);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------------- */ 
/*                int FastaSequenceSoft2HardMask_Main()                          */
/* This function is used for hard masking a FASTA file.                          */
/* The Maskfile should have the following format:                                */
/* All a,c,g,t in little case will be converted to N.                            */
/* ----------------------------------------------------------------------------- */
int FastaSequenceSoft2HardMask_Main(char strSeqFile[], char strOutFile[])
{
	FILE *fpOut;
	struct tagSequence *pSeqList = NULL;
	int nSeqCount = 0;
	struct tagSequence **vSeqList = NULL;
	struct tagSequence *pSeq;
	int ni;
	char *vB;

	/* load sequence */
	nSeqCount = LoadFullSequenceList(strSeqFile, &pSeqList);
	if(nSeqCount <= 0)
	{
		printf("Warning: No sequences loaded!\n");
		SequenceListClear(&pSeqList);
		return PROC_SUCCESS;
	}

	/* reorganize sequence */
	vSeqList = NULL;
	vSeqList = (struct tagSequence **)calloc(nSeqCount, sizeof(struct tagSequence *));
	if(vSeqList == NULL)
	{
		printf("Error: cannot allocate memory for sequence list!\n");
		SequenceListClear(&pSeqList);
		exit(EXIT_FAILURE);
	}
	
	ni = 0;
	pSeq = pSeqList;
	while(pSeq != NULL)
	{
		vSeqList[ni] = pSeq;
		ni++;
		pSeq = pSeq->m_pNext;
	}
	if(ni != nSeqCount)
	{
		printf("Error: sequence number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* mask */
	for(ni=0; ni<nSeqCount; ni++)
	{
		vB = vSeqList[ni]->m_pSequence->m_pString;
		while(*vB != '\0')
		{
			if( (*vB != 'A') && (*vB != 'C') && (*vB != 'G') && (*vB != 'T') )
				*vB = 'N';
			vB++;
		}
	}

	/* export */
	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nSeqCount; ni++)
	{
		SequenceWriteToFasta_ByStrand(vSeqList[ni], fpOut, '+', 1);
	}
	
	fclose(fpOut);

	/* clear memory */
	SequenceListClear(&pSeqList);
	free(vSeqList);

	/* return */
	return PROC_SUCCESS;
}