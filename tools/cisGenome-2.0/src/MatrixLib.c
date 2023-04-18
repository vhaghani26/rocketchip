/* ----------------------------------------------------------------------- */
/*  MatrixLib.c : implementation of the matrix process library             */
/*  Author : Ji HongKai ; Time: 1999.11                                    */
/* ----------------------------------------------------------------------- */

#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
#include "limits.h"
#include "time.h"

#include "StringLib.h"
#include "MatrixLib.h"
#include "RandomLib.h"

/* ----------------------------------------------------------------------- */
/* CreateByteMatrix Function                                               */
/* This function is used for creating a new BYTEMATRIX .                   */
/* The size of the matrix is given by two parameters : nHeight , nWidth    */
/* In other words , a nHeight*nWidth matrix will be created.               */
/* If it fails to create new matrix , the return value will be NULL .      */
/* Else the function will return the pointer to the new BYTEMATRIX .       */
/* ----------------------------------------------------------------------- */
struct BYTEMATRIX* CreateByteMatrix(long nHeight, long nWidth)
{
	/* pNewMatrix is the pointer to the new matrix. */
	struct BYTEMATRIX* pNewMatrix;
	/* pNewElement is the pointer to the block of matrix elements. */
	unsigned char * pNewElement;
	/* nSize is the size of matrix elements. */
	long lnSize;

	/* if nHeight<=0 or nWidth<=0 , return NULL */
	if((nHeight<=0)||(nWidth<=0))
		return NULL;
	/* allocates memory to struct BYTEMATRIX. */
	pNewMatrix = NULL;
	pNewMatrix = (struct BYTEMATRIX*)malloc(sizeof(struct BYTEMATRIX));
	if (pNewMatrix == NULL)
	{
		printf("Insufficient memory available for creating bytematrix.\n");
	}
	else
	{
		/* initialize new BYTEMATRIX. */
		pNewMatrix->nHeight = nHeight;
		pNewMatrix->nWidth = nWidth;

		/* calculate size of matrix elements. */
		lnSize = (long)(nHeight*nWidth);
		/* allocate memory for matrix elements. */
		pNewElement = NULL;
		pNewElement = (unsigned char *)calloc(lnSize, sizeof(unsigned char));
		if( pNewElement != NULL )
		{
			pNewMatrix->pMatElement = pNewElement;
		}
		else
		{
			printf( "Can't allocate memory for bytematrix elements.\n" );
			free(pNewElement);
			free(pNewMatrix);
			pNewMatrix = NULL;
		}
	}

	/* return pNewMatrix. */
	return pNewMatrix;
}

/* ----------------------------------------------------------------------- */
/* DestroyByteMatrix Function                                              */
/* This function is used for destroying a BYTEMATRIX .                     */
/* First, the memory block of elements will be freed .                     */
/* Then the memory block of the BYTEMATRIX struct will be freed .          */
/* ----------------------------------------------------------------------- */
void DestroyByteMatrix(struct BYTEMATRIX* pDelMatrix)
{
	/* pDelElement is the pointer to the block of matrix elements. */
	unsigned char * pDelElement;

	/* delete pDelMatrix if it exists. */
	if (pDelMatrix != NULL)
	{
		/* free memory of matrix elements. */
		pDelElement = pDelMatrix->pMatElement;
		free (pDelElement);
		/* free memory of matrix strunct. */
		free (pDelMatrix);
		pDelMatrix = NULL;
	}
}

/* ----------------------------------------------------------------------- */
/* BMGETAT Function                                                        */
/* This function is used for getting element indicated by row index and    */
/* column index .                                                          */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateByteMatrix function to create it          */
/* before you call this function !                                         */
/* If the row(col) indexes is larger than or equal to the row(col) number  */
/* of the matrix or smaller than 0 or some other exceptions occur, the     */
/* function will return BM_ACCESS_VIOLATION;                               */  
/* Else the function will return matrix[row][col] .                        */
/* ----------------------------------------------------------------------- */
unsigned char BMGETAT(struct BYTEMATRIX* pSourceMatrix,long nRow,long nCol)
{
	/* nHeight and nWidth are dimension informations of source matrix. */
	long nHeight,nWidth;
	/* pSourceElement is used for accessing the element of matrixes . */
	unsigned char *pSourceElement;
	/* nResult is the return value */
	unsigned char bResult;

	/* if source matrix doesn't exist , return BM_ACCESS_VIOLATION */
	if(pSourceMatrix==NULL)
	{
		printf("Access Violation!\n"); 
		return BM_ACCESS_VIOLATION;
	}
	/* if source matrix exists, do next */
	else
	{
		/* get dimension information */
		nHeight = pSourceMatrix->nHeight;
		nWidth = pSourceMatrix->nWidth;
		/* if nRow(nCol) < 0 or >=nHeight(nWidth),*/
		/* return BM_ACCESS_VIOLATION */
		if((nRow<0)||(nRow>=nHeight)||(nCol<0)||(nCol>=nWidth))
		{
			printf("Access Violation!\n");
			return BM_ACCESS_VIOLATION;
		}
		/* if 0<= nRow(nCol) <nHeight(nWidth), */
		/* return matrix[nRow][nCol] .  */
		else
		{
			pSourceElement = pSourceMatrix->pMatElement;
			bResult = *(pSourceElement+nRow*nWidth+nCol);
			/* return */
			return bResult;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* BMSETAT Function                                                        */
/* This function is used for setting element indicated by row index and    */
/* column index .                                                          */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateByteMatrix function to create it          */
/* before you call this function !                                         */
/* If the row(col) indexes is larger than or equal to the row(col) number  */
/* of the matrix or smaller than 0 or some other exceptions occur, the     */
/* function will return PROC_FAILURE;                                      */  
/* Else the function will let matrix[row][col] = nValue and                */
/* return PROC_SUCCESS.                                                    */
/* ----------------------------------------------------------------------- */
int BMSETAT(struct BYTEMATRIX* pSourceMatrix,long nRow,long nCol,unsigned char bValue)
{
	/* nHeight and nWidth are dimension informations of source matrix. */
	long nHeight,nWidth;
	/* pSourceElement is used for accessing the element of matrixes . */
	unsigned char *pSourceElement;
	
	/* if source matrix doesn't exist , return PROC_FAILURE */
	if(pSourceMatrix==NULL)
	{
		printf("Access Violation!\n");
		return PROC_FAILURE;
	}
	/* if source matrix exists, do next */
	else
	{
		/* get dimension information */
		nHeight = pSourceMatrix->nHeight;
		nWidth = pSourceMatrix->nWidth;
		/* if nRow(nCol) < 0 or >=nHeight(nWidth),*/
		/* return DM_ACCESS_VIOLATION */
		if((nRow<0)||(nRow>=nHeight)||(nCol<0)||(nCol>=nWidth))
		{
			printf("Access Violation!\n");
			return PROC_FAILURE;
		}
		/* if 0<= nRow(nCol) <nHeight(nWidth), */
		/* let matrix[nRow][nCol]=nValue and .  */
		/* return PROC_SUCCESS */
		else
		{
			pSourceElement = pSourceMatrix->pMatElement;
			*(pSourceElement+nRow*nWidth+nCol) = bValue;
			/* return */
			return PROC_SUCCESS;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* BMSAVE Function                                                         */
/* This function is used for saving a BYTEMATRIX.                          */
/* If successful, return PROC_SUCCESS;                                     */
/* else return PROC_FAILURE.                                               */
/* ----------------------------------------------------------------------- */
int BMSAVE(struct BYTEMATRIX* pByteMatrix, char strFilePath[])
{
	/* fpSaveFile is used for saving. */
	FILE *fpSaveFile;
	/* nRow, nCol is used for access the matrix */
	long nRow, nCol;
	/* pElement is used for access the elements of the matrix. */
	unsigned char *pElement;
	
	/* open file */
	fpSaveFile = NULL;
	fpSaveFile = fopen(strFilePath, "wt");
	if(fpSaveFile == NULL)
	{
		printf("Can't open output file!\n");
		return PROC_FAILURE;
	}

	/* check parameter */
	if(pByteMatrix == NULL)
		return PROC_SUCCESS;

	/* save */
	pElement = pByteMatrix->pMatElement;
	for(nRow=0; nRow<pByteMatrix->nHeight; nRow++)
	{
		for(nCol=0; nCol<pByteMatrix->nWidth; nCol++)
		{
			fprintf(fpSaveFile, "%d ", *pElement);
			pElement++;
		}
		fprintf(fpSaveFile, "\n");
	}

	/* close file */
	fclose(fpSaveFile);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */
/* BMLOAD Function                                                         */
/* This function is used for loading a BYTEMATRIX.                         */
/* If successful, return the pointer of the matrix;                        */
/* else return NULL.                                                       */
/* ----------------------------------------------------------------------- */
struct BYTEMATRIX* BMLOAD(char strFilePath[])
{
	/* pMatrix is used for loading matrix */
	struct BYTEMATRIX* pMatrix;
	/* fpLoadFile is used for loading. */
	FILE *fpLoadFile;
	/* pElement and nElement is used for access the elements of the matrix. */
	unsigned char bElement[BM_LOAD_MAX_WIDTH];
	long nElementCount;
	/* strInLine is used for loading file. */
	char strInLine[LONG_LINE_LENGTH];

	/* init */
	pMatrix = NULL;

	/* open file */
	fpLoadFile = NULL;
	fpLoadFile = fopen(strFilePath, "rt");
	if(fpLoadFile == NULL)
	{
		printf("Can't open input file!\n");
		return NULL;
	}

	/* load the first line */
	while(fgets(strInLine, LONG_LINE_LENGTH, fpLoadFile) != NULL)
	{
		/* trim right and left */
		nElementCount = ByteLoadRowVector(strInLine, bElement);
		if(nElementCount == 0)
			continue;
		if(BMADDROW(&pMatrix, bElement, nElementCount) == PROC_FAILURE)
		{
			DestroyByteMatrix(pMatrix);
			pMatrix = NULL;
			printf("Can't load int matrix!\n");
			break;
		}
	}

	/* close file */
	fclose(fpLoadFile);

	/* return */
	return pMatrix;
}

/* ----------------------------------------------------------------------- */
/* ByteLoadRowVector                                                       */
/* This function is used for loading a byte row vector from a string.      */
/* If successful, return the count of loaded elements;                     */
/* else return 0.                                                          */
/* ----------------------------------------------------------------------- */
long ByteLoadRowVector(char strInLine[], unsigned char bElement[])
{
	/* nCount is count of loaded elements. */
	long nCount;
	/* startp and endp are pointers to the string */
	char *startp;
	char *endp;
	double dTemp;
 
	/* init */
	nCount = 0;
	StrTrimRight(strInLine);
	StrTrimLeft(strInLine);
	if(strInLine[0] == '\0')
		return nCount;

	/* load */
	startp = strInLine;
	endp = strpbrk(startp, WORD_SEPARATORS);

	while(endp != NULL)
	{
		*endp = '\0';
		dTemp = atof(startp);
		bElement[nCount] = (unsigned char)dTemp;
		nCount++;
		startp = endp+1;
		StrTrimLeft(startp);
		endp = strpbrk(startp, WORD_SEPARATORS);
	}

	dTemp = atof(startp);
	bElement[nCount] = (unsigned char)dTemp;
	nCount++;

	/* return */
	return nCount;
}

/* ----------------------------------------------------------------------- */
/* BMADDROW Function                                                       */
/* This function is used for adding a row to BYTEMATRIX.                   */
/* If successful, return PROC_SUCCESS;                                     */
/* else return PROC_FAILURE.                                               */
/* ----------------------------------------------------------------------- */
int BMADDROW(struct BYTEMATRIX **pByteMatrix, unsigned char bElement[], long nElementCount)
{
	/* define */
	unsigned char *pElement;
	long lnSize;

	/* check parameter */
	if(pByteMatrix == NULL)
		return PROC_FAILURE;

	/* if pByteMatrix doesn't exist */
	if(*pByteMatrix == NULL)
	{
		*pByteMatrix = CreateByteMatrix(1, nElementCount);
		if(*pByteMatrix == NULL)
			return PROC_FAILURE;
		
		pElement = (*pByteMatrix)->pMatElement;
		memcpy(pElement, bElement, nElementCount*sizeof(unsigned char));

		return PROC_SUCCESS; 
	}

	/* if pByteMatrix does exist */
	if((*pByteMatrix)->nWidth != nElementCount)
		return PROC_FAILURE;

	lnSize = (long)(((*pByteMatrix)->nHeight)*((*pByteMatrix)->nWidth));
	pElement = (*pByteMatrix)->pMatElement;
	(*pByteMatrix)->pMatElement = (unsigned char*)realloc(pElement, (lnSize+nElementCount)*sizeof(unsigned char));
	if((*pByteMatrix)->pMatElement == NULL)
		exit(EXIT_FAILURE);

	pElement = (*pByteMatrix)->pMatElement+lnSize;
	memcpy(pElement, bElement, nElementCount*sizeof(unsigned char));
	(*pByteMatrix)->nHeight += 1;

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */
/* BM_T Function                                                           */
/* This function creates a new BYTEMATRIX pY, and lets it equal to pX'.    */
/* !note: pX must exist .                                                  */
/* Make sure you have used CreateByteMatrix function to create it          */
/* before you call this function !                                         */
/* If the new matrix can't be created or other exceptions occur , the      */
/* function will return NULL ;                                             */
/* Else the function will translocate the given matrix and return the      */
/* pointer of the new matrix .                                             */
/* ----------------------------------------------------------------------- */
struct BYTEMATRIX* BM_T(struct BYTEMATRIX* pX)
{
	/* pY is the pointer to the new matrix. */
	struct BYTEMATRIX* pY;
	/* nHeight and nWidth are dimension informations of matrix. */
	long nHeight,nWidth;
	long nRow,nCol;
	unsigned char *pXelement,*pYelement;

	/* if source matrix doesn't exist , return NULL */
	if(pX==NULL)
	{
		return NULL;
	}
	/* if sorce matrix exists , create it's clone */
	/* get dimension information */
	nHeight = pX->nWidth;
	nWidth = pX->nHeight;

	/* create a new matrix */
	pY = NULL;
	pY = CreateByteMatrix(nHeight,nWidth);
	if(pY == NULL)
	{
		return NULL;
	}

	/* Y = X' */
	pYelement = pY->pMatElement;
	for(nRow=0; nRow<nHeight; nRow++)
	{
		pXelement = pX->pMatElement + nRow;
		for(nCol=0; nCol<nWidth; nCol++)
		{
			*pYelement = *pXelement;
			pYelement++;
			pXelement = pXelement + nHeight;
		}
	}

	/* return */
	return pY;
}

/* ----------------------------------------------------------------------- */
/* BMCOPY Function                                                         */
/* This function lets pDestByteMatrix = pSourceByteMatrix .                */
/* !note: pDestByteMatrix and pSourceByteMatrix must exist .               */
/* Make sure you have used CreateByteMatrix function to create them        */
/* before you call this function !                                         */
/* If the dimensions of the two matrix are not equal , the function will   */
/* return PROC_FAILURE;                                                    */
/* Else the function will copy pSourceByteMatrix to pDestByteMatrix        */
/* and return PROC_SUCCESS .                                               */
/* ----------------------------------------------------------------------- */
int BMCOPY(struct BYTEMATRIX* pDestByteMatrix,struct BYTEMATRIX* pSourceByteMatrix)
{
	long lnSize;
	unsigned char *pDestElement,*pSourceElement;

	/* if one of the two matrixes doesn't exist , return PROC_FAILURE */
	if((pDestByteMatrix==NULL)||(pSourceByteMatrix==NULL))
	{
		return PROC_FAILURE;
	}
	/* if both matrixes exist , do next . */
	else 
	{
		/* if the dimensions of two matrixes are equal , copy source to dest*/
		/* and return PROC_SUCCESS */
		if((pDestByteMatrix->nHeight==pSourceByteMatrix->nHeight)&&
			(pDestByteMatrix->nWidth==pSourceByteMatrix->nWidth))
		{
			pDestElement = pDestByteMatrix->pMatElement;
			pSourceElement = pSourceByteMatrix->pMatElement;
			lnSize = (long)((pSourceByteMatrix->nHeight)*(pSourceByteMatrix->nWidth));
			memcpy(pDestElement, pSourceElement, lnSize*sizeof(unsigned char));

			return PROC_SUCCESS;
		}
		/* if the dimensions of two matrixes are not equal ,*/
		/* return PROC_FAILURE */
		else
		{
			return PROC_FAILURE;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* BMCLONE Function                                                        */
/* This function creates a new BYTEMATRIX , and lets it equal to the       */
/* given BYTEMATRIX .                                                      */
/* !note: pSourceByteMatrix must exist .                                   */
/* Make sure you have used CreateByteMatrix function to create it          */
/* before you call this function !                                         */
/* If the new matrix can't be created or other exceptions occur , the      */
/* function will return NULL ;                                             */
/* Else the function will clone the given matrix and return the pointer of */
/* the new matrix .                                                        */
/* ----------------------------------------------------------------------- */
struct BYTEMATRIX* BMCLONE(struct BYTEMATRIX* pSourceByteMatrix)
{
	/* pNewMatrix is the pointer to the new matrix. */
	struct BYTEMATRIX* pNewMatrix;
	/* nHeight and nWidth are dimension informations of matrix. */
	long nHeight,nWidth;

	/* if source matrix doesn't exist , return NULL */
	if(pSourceByteMatrix==NULL)
	{
		return NULL;
	}
	/* if sorce matrix exists , create it's clone */
	else
	{
		/* get dimension information */
		nHeight = pSourceByteMatrix->nHeight;
		nWidth = pSourceByteMatrix->nWidth;
		/* create a new matrix */
		pNewMatrix = CreateByteMatrix(nHeight,nWidth);
		/* copy the source to dest */
		if(BMCOPY(pNewMatrix,pSourceByteMatrix)==PROC_SUCCESS)
		{
			return pNewMatrix;
		}
		else
		{
			DestroyByteMatrix(pNewMatrix);
			return NULL;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* BMSUM Function                                                          */
/* This function computes the sum of all elements of a BM matrix.          */
/* ----------------------------------------------------------------------- */
double BMSUM(struct BYTEMATRIX* pByteMatrix)
{
	long ni,nj;
	unsigned char *pElement;
	double dSum = 0.0;

	/* if the matrix doesn't exist , return 0 */
	if(pByteMatrix==NULL)
	{
		printf("Warning: BMSUM, matrix does not exist!\n");
		return dSum;
	}
	/* if both matrixes exist , do next . */
	else 
	{
		/* if the dimensions of two matrixes are equal , copy source to dest*/
		/* and return PROC_SUCCESS */
		pElement = pByteMatrix->pMatElement;
		for(ni=0; ni<pByteMatrix->nHeight; ni++)
		{
			for(nj=0; nj<pByteMatrix->nWidth; nj++)
			{
				dSum += (*pElement);
				pElement++;
			}
		}

		return dSum;
	}
}

/* ----------------------------------------------------------------------- */
/* BMPERMUTEROWS Function                                                  */
/* This function permutes the rows of the source matrix, and retrun a new  */
/* matrix which is the permutated one.                                     */
/* If the new matrix can't be created or other exceptions occur , the      */
/* function will return NULL ;                                             */
/* Else the function will permute the source matrix and return the pointer */
/* of the new matrix.                                                      */
/* ----------------------------------------------------------------------- */
struct BYTEMATRIX* BMPERMUTEROWS(struct BYTEMATRIX* pSourceByteMatrix)
{
	/* pNewMatrix is the pointer to the new matrix. */
	struct BYTEMATRIX* pNewMatrix;
	/* nHeight and nWidth are dimension informations of matrix. */
	long nHeight,nWidth,ni;
	/* pRandMatrix, pSortMatrix and pSortIndex are used to create a permutation. */
	struct DOUBLEMATRIX* pRandMatrix;
	struct DOUBLEMATRIX* pSortMatrix;
	struct LONGMATRIX* pSortIndex;
	/* pSourceStart, pDestStart and pIdx are used for copying. */
	unsigned char *pSourceStart,*pDestStart;
	long *pIdx;

	/* if source matrix doesn't exist , return NULL */
	if(pSourceByteMatrix==NULL)
	{
		return NULL;
	}
	
	/* if sorce matrix exists , permute its rows */
	pNewMatrix = NULL;

	/* get dimension information */
	nHeight = pSourceByteMatrix->nHeight;
	nWidth = pSourceByteMatrix->nWidth;
	
	/* create a new matrix */
	pNewMatrix = CreateByteMatrix(nHeight,nWidth);
	if (pNewMatrix == NULL)
	{
		printf("Insufficient memory available for creating bytematrix while doing permutation.\n");
		return NULL;
	}

	/* create a permutation */
	pRandMatrix = NULL;
	pSortMatrix = NULL;
	pSortIndex = NULL;
	pRandMatrix = DMRANDU(1,nHeight);
	if(pRandMatrix == NULL)
	{
		DestroyByteMatrix(pNewMatrix);
		printf("Can't create permutation.\n");
		return NULL;
	}
	if(DMSORTMERGEA_0(pRandMatrix, &pSortMatrix, &pSortIndex) == PROC_FAILURE)
	{
		DestroyDoubleMatrix(pRandMatrix);
		DestroyByteMatrix(pNewMatrix);
		DestroyDoubleMatrix(pSortMatrix);
		DestroyLongMatrix(pSortIndex);
		printf("Can't create permutation.\n");
		return NULL;
	}

	/* release memory */
	DestroyDoubleMatrix(pRandMatrix);
	DestroyDoubleMatrix(pSortMatrix);
	

	/* create new matrix */
	pDestStart = pNewMatrix->pMatElement;
	pIdx = pSortIndex->pMatElement;
	for(ni=0; ni<nHeight; ni++)
	{
		pSourceStart = pSourceByteMatrix->pMatElement+(*pIdx)*pSourceByteMatrix->nWidth;
		memcpy(pDestStart, pSourceStart, nWidth*sizeof(unsigned char));

		/* get next */
		pDestStart = pDestStart+pNewMatrix->nWidth;
		pIdx++;
	}

	/* release memory */
	DestroyLongMatrix(pSortIndex);

	/* return */
	return pNewMatrix;
}


/* ----------------------------------------------------------------------- */
/* BMPERMUTEELEMENTINONECOLUMN Function                                    */
/* This function permutes the elements of the specific column of the       */
/* source matrix, and retrun a new matrix which is the permutated one.     */
/* If the new matrix can't be created or other exceptions occur , the      */
/* function will return NULL ;                                             */
/* Else the function will permute the source matrix and return the pointer */
/* of the new matrix.                                                      */
/* ----------------------------------------------------------------------- */
struct BYTEMATRIX* BMPERMUTEELEMENTSINONECOLUMN(struct BYTEMATRIX* pSourceByteMatrix, long nColumnIndex)
{
	/* pNewMatrix is the pointer to the new matrix. */
	struct BYTEMATRIX* pNewMatrix;
	/* nHeight and nWidth are dimension informations of matrix. */
	long nHeight,nWidth,ni;
	/* pRandMatrix, pSortMatrix and pSortIndex are used to create a permutation. */
	struct DOUBLEMATRIX* pRandMatrix;
	struct DOUBLEMATRIX* pSortMatrix;
	struct LONGMATRIX* pSortIndex;
	/* pSourceStart, pDestStart and pIdx are used for copying. */
	unsigned char *pSourceStart,*pDestStart;
	long *pIdx;

	/* if source matrix doesn't exist , return NULL */
	if(pSourceByteMatrix==NULL)
	{
		return NULL;
	}
	
	/* get dimension information */
	nHeight = pSourceByteMatrix->nHeight;
	nWidth = pSourceByteMatrix->nWidth;

	/* if sorce matrix exists , permute its rows */
	pNewMatrix = NULL;

	/* create a new matrix */
	pNewMatrix = BMCLONE(pSourceByteMatrix);
	if (pNewMatrix == NULL)
	{
		printf("Insufficient memory available for creating bytematrix while doing permutation.\n");
		return NULL;
	}

	if((nColumnIndex >= nWidth) || (nColumnIndex < 0))
	{
		printf("Warning: Permutation column index overflows!\n"); 
		return pNewMatrix;
	}

	/* create a permutation */
	pRandMatrix = NULL;
	pSortMatrix = NULL;
	pSortIndex = NULL;
	pRandMatrix = DMRANDU(1,nHeight);
	if(pRandMatrix == NULL)
	{
		DestroyByteMatrix(pNewMatrix);
		printf("Can't create permutation.\n");
		return NULL;
	}
	if(DMSORTMERGEA_0(pRandMatrix, &pSortMatrix, &pSortIndex) == PROC_FAILURE)
	{
		DestroyDoubleMatrix(pRandMatrix);
		DestroyByteMatrix(pNewMatrix);
		DestroyDoubleMatrix(pSortMatrix);
		DestroyLongMatrix(pSortIndex);
		printf("Can't create permutation.\n");
		return NULL;
	}

	/* release memory */
	DestroyDoubleMatrix(pRandMatrix);
	DestroyDoubleMatrix(pSortMatrix);
	

	/* create new matrix */
	pDestStart = pNewMatrix->pMatElement+nColumnIndex;
	pIdx = pSortIndex->pMatElement;
	for(ni=0; ni<nHeight; ni++)
	{
		pSourceStart = pSourceByteMatrix->pMatElement+(*pIdx)*pSourceByteMatrix->nWidth+nColumnIndex;
		*pDestStart = *pSourceStart;

		/* get next */
		pDestStart = pDestStart+pNewMatrix->nWidth;
		pIdx++;
	}

	/* release memory */
	DestroyLongMatrix(pSortIndex);

	/* return */
	return pNewMatrix;
}


/* ----------------------------------------------------------------------- */
/* BMBOOTSTRAP Function                                                    */
/* This function lets pDestByteMatrix = bootstrap(pSourceByteMatrix).      */
/* !note: pDestByteMatrix and pSourceByteMatrix must exist .               */
/* Make sure you have used CreateByteMatrix function to create them        */
/* before you call this function !                                         */
/* If the dimensions of the two matrix are not equal , the function will   */
/* return PROC_FAILURE;                                                    */
/* Else return PROC_SUCCESS .                                              */
/* ----------------------------------------------------------------------- */
int BMBOOTSTRAP(struct BYTEMATRIX* pDestByteMatrix,struct BYTEMATRIX* pSourceByteMatrix)
{
	/* define */
	unsigned char *pDestElement,*pSourceElement;
	long nRow1,nRow2;
	double d_Rand;

	/* if one of the two matrixes doesn't exist , return PROC_FAILURE */
	if((pDestByteMatrix==NULL)||(pSourceByteMatrix==NULL))
	{
		return PROC_FAILURE;
	}
	/* if both matrixes exist , do next . */
	else 
	{
		/* if the dimensions of two matrixes are equal , create bootstrap data */
		/* and return PROC_SUCCESS */
		if(pDestByteMatrix->nWidth==pSourceByteMatrix->nWidth)
		{
			pDestElement = pDestByteMatrix->pMatElement;
			for(nRow1=0; nRow1<pDestByteMatrix->nHeight; nRow1++)
			{
				d_Rand = rand_u();
				nRow2 = (long)(d_Rand*pSourceByteMatrix->nHeight);
				if(nRow2 == pSourceByteMatrix->nHeight)
					nRow2--;

				pSourceElement = pSourceByteMatrix->pMatElement+nRow2*pSourceByteMatrix->nWidth;
				memcpy(pDestElement, pSourceElement, pDestByteMatrix->nWidth*sizeof(unsigned char));

				/* get next */
				pDestElement = pDestElement+pDestByteMatrix->nWidth;
			}
			
			return PROC_SUCCESS;
		}
		/* if the dimensions of two matrixes are not equal ,*/
		/* return PROC_FAILURE */
		else
		{
			return PROC_FAILURE;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* BMBOOTSTRAPC Function                                                   */
/* This function lets pDestByteMatrix = bootstrap(pSourceByteMatrix) along */
/* colum.                                                                  */
/* !note: pDestByteMatrix and pSourceByteMatrix must exist .               */
/* Make sure you have used CreateByteMatrix function to create them        */
/* before you call this function !                                         */
/* If the dimensions of the two matrix are not equal , the function will   */
/* return PROC_FAILURE;                                                    */
/* Else return PROC_SUCCESS .                                              */
/* ----------------------------------------------------------------------- */
int BMBOOTSTRAPC(struct BYTEMATRIX* pDestByteMatrix,struct BYTEMATRIX* pSourceByteMatrix)
{
	/* define */
	unsigned char *pDestElement,*ch,*pSourceElement;
	long nCol1,nCol2;
	double d_Rand;
	long nj;

	/* if one of the two matrixes doesn't exist , return PROC_FAILURE */
	if((pDestByteMatrix==NULL)||(pSourceByteMatrix==NULL))
	{
		return PROC_FAILURE;
	}
	/* if both matrixes exist , do next . */
	else 
	{
		/* if the dimensions of two matrixes are equal , create bootstrap data */
		/* and return PROC_SUCCESS */
		if(pDestByteMatrix->nHeight==pSourceByteMatrix->nHeight)
		{
			pDestElement = pDestByteMatrix->pMatElement;
			for(nCol1=0; nCol1<pDestByteMatrix->nWidth; nCol1++)
			{
				d_Rand = rand_u();
				nCol2 = (long)(d_Rand*pSourceByteMatrix->nWidth);
				if(nCol2 == pSourceByteMatrix->nWidth)
					nCol2--;

				pSourceElement = pSourceByteMatrix->pMatElement+nCol2;
				
				ch = pDestElement;
				for(nj=0; nj<pDestByteMatrix->nHeight; nj++)
				{
					*ch = *pSourceElement;
					ch = ch+pDestByteMatrix->nWidth;
					pSourceElement = pSourceElement+pSourceByteMatrix->nWidth;
				}

				/* get next */
				pDestElement++;
			}
			
			return PROC_SUCCESS;
		}
		/* if the dimensions of two matrixes are not equal ,*/
		/* return PROC_FAILURE */
		else
		{
			return PROC_FAILURE;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* CreateIntMatrix Function                                                */
/* This function is used for creating a new INTMATRIX .                    */
/* The size of the matrix is given by two parameters : nHeight , nWidth    */
/* In other words , a nHeight*nWidth matrix will be created.               */
/* If it fails to create new matrix , the return value will be NULL .      */
/* Else the function will return the pointer to the new INTMATRIX .        */
/* ----------------------------------------------------------------------- */
struct INTMATRIX* CreateIntMatrix(long nHeight, long nWidth)
{
	/* pNewMatrix is the pointer to the new matrix. */
	struct INTMATRIX* pNewMatrix;
	/* pNewElement is the pointer to the block of matrix elements. */
	int *pNewElement;
	/* nSize is the size of matrix elements. */
	long lnSize;

	/* if nHeight=0 or nWidth=0 , return NULL */
	if((nHeight<=0)||(nWidth<=0))
		return NULL;
	/* allocates memory to struct INTMATRIX. */
	pNewMatrix = NULL;
	pNewMatrix = (struct INTMATRIX*)malloc(sizeof(struct INTMATRIX));
	if (pNewMatrix == NULL)
	{
		printf("Insufficient memory available for creating intmatrix.\n"); 
	}
	else
	{
		/* initialize new INTMATRIX. */
		pNewMatrix->nHeight = nHeight;
		pNewMatrix->nWidth = nWidth;

		/* calculate size of matrix elements. */
		lnSize = (long)(nHeight*nWidth);
		/* allocate memory for matrix elements. */
		pNewElement = NULL;
		pNewElement = (int *)calloc(lnSize, sizeof(int));
		if( pNewElement != NULL )
		{
			pNewMatrix->pMatElement = pNewElement;
		}
		else
		{
			printf( "Can't allocate memory for intmatrix elements.\n" );
			free(pNewElement);
			free(pNewMatrix);
			pNewMatrix = NULL;
		}
	}

	/* return pNewMatrix. */
	return pNewMatrix;
}

/* ----------------------------------------------------------------------- */
/* DestroyIntMatrix Function                                               */
/* This function is used for destroying a INTMATRIX .                      */
/* First, the memory block of elements will be freed .                     */
/* Then the memory block of the INTMATRIX struct will be freed .           */
/* ----------------------------------------------------------------------- */
void DestroyIntMatrix(struct INTMATRIX* pDelMatrix)
{
	/* pDelElement is the pointer to the block of matrix elements. */
	int *pDelElement;

	/* delete pDelMatrix if it exists. */
	if (pDelMatrix != NULL)
	{
		/* free memory of matrix elements. */
		pDelElement = pDelMatrix->pMatElement;
		free (pDelElement);
		/* free memory of matrix strunct. */
		free (pDelMatrix);
		pDelMatrix = NULL;
	}
} 

/* ----------------------------------------------------------------------- */
/* IMGETAT Function                                                        */
/* This function is used for getting element indicated by row index and    */
/* column index .                                                          */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateINTMatrix function to create it           */
/* before you call this function !                                         */
/* If the row(col) indexes is larger than or equal to the row(col) number  */
/* of the matrix or smaller than 0 or some other exceptions occur, the     */
/* function will return IM_ACCESS_VIOLATION;                               */  
/* Else the function will return matrix[row][col] .                        */
/* ----------------------------------------------------------------------- */
int IMGETAT(struct INTMATRIX* pSourceMatrix,long nRow,long nCol)
{
	/* nHeight and nWidth are dimension informations of source matrix. */
	long nHeight,nWidth;
	/* pSourceElement is used for accessing the element of matrixes . */
	int *pSourceElement;
	/* nResult is the return value */
	int nResult;

	/* if source matrix doesn't exist , return IM_ACCESS_VIOLATION */
	if(pSourceMatrix==NULL)
	{
		printf("Access Violation!\n"); 
		return IM_ACCESS_VIOLATION;
	}
	/* if source matrix exists, do next */
	else
	{
		/* get dimension information */
		nHeight = pSourceMatrix->nHeight;
		nWidth = pSourceMatrix->nWidth;
		/* if nRow(nCol) < 0 or >=nHeight(nWidth),*/
		/* return DM_ACCESS_VIOLATION */
		if((nRow<0)||(nRow>=nHeight)||(nCol<0)||(nCol>=nWidth))
		{
			printf("Access Violation!\n");
			return IM_ACCESS_VIOLATION;
		}
		/* if 0<= nRow(nCol) <nHeight(nWidth), */
		/* return matrix[nRow][nCol] .  */
		else
		{
			pSourceElement = pSourceMatrix->pMatElement;
			nResult = *(pSourceElement+nRow*nWidth+nCol);
			/* return */
			return nResult;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* IMSETAT Function                                                        */
/* This function is used for setting element indicated by row index and    */
/* column index .                                                          */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateIntMatrix function to create it           */
/* before you call this function !                                         */
/* If the row(col) indexes is larger than or equal to the row(col) number  */
/* of the matrix or smaller than 0 or some other exceptions occur, the     */
/* function will return PROC_FAILURE;                                      */  
/* Else the function will let matrix[row][col] = nValue and                */
/* return PROC_SUCCESS.                                                    */
/* ----------------------------------------------------------------------- */
int IMSETAT(struct INTMATRIX* pSourceMatrix,long nRow,long nCol,int nValue)
{
	/* nHeight and nWidth are dimension informations of source matrix. */
	long nHeight,nWidth;
	/* pSourceElement is used for accessing the element of matrixes . */
	int *pSourceElement;
	
	/* if source matrix doesn't exist , return PROC_FAILURE */
	if(pSourceMatrix==NULL)
	{
		printf("Access Violation!\n");
		return PROC_FAILURE;
	}
	/* if source matrix exists, do next */
	else
	{
		/* get dimension information */
		nHeight = pSourceMatrix->nHeight;
		nWidth = pSourceMatrix->nWidth;
		/* if nRow(nCol) < 0 or >=nHeight(nWidth),*/
		/* return DM_ACCESS_VIOLATION */
		if((nRow<0)||(nRow>=nHeight)||(nCol<0)||(nCol>=nWidth))
		{
			printf("Access Violation!\n");
			return PROC_FAILURE;
		}
		/* if 0<= nRow(nCol) <nHeight(nWidth), */
		/* let matrix[nRow][nCol]=nValue and .  */
		/* return PROC_SUCCESS */
		else
		{
			pSourceElement = pSourceMatrix->pMatElement;
			*(pSourceElement+nRow*nWidth+nCol) = nValue;
			/* return */
			return PROC_SUCCESS;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* IMSAVE Function                                                         */
/* This function is used for saving a INTMATRIX.                           */
/* If successful, return PROC_SUCCESS;                                     */
/* else return PROC_FAILURE.                                               */
/* ----------------------------------------------------------------------- */
int IMSAVE(struct INTMATRIX* pIntMatrix, char strFilePath[])
{
	/* fpSaveFile is used for saving. */
	FILE *fpSaveFile;
	/* nRow, nCol is used for access the matrix */
	long nRow, nCol;
	/* pElement is used for access the elements of the matrix. */
	int *pElement;
	
	/* open file */
	fpSaveFile = NULL;
	fpSaveFile = fopen(strFilePath, "wt");
	if(fpSaveFile == NULL)
	{
		printf("Can't open output file!\n");
		return PROC_FAILURE;
	}

	/* check parameter */
	if(pIntMatrix == NULL)
		return PROC_SUCCESS;

	/* save */
	pElement = pIntMatrix->pMatElement;
	for(nRow=0; nRow<pIntMatrix->nHeight; nRow++)
	{
		for(nCol=0; nCol<pIntMatrix->nWidth; nCol++)
		{
			fprintf(fpSaveFile, "%d ", *pElement);
			pElement++;
		}
		fprintf(fpSaveFile, "\n");
	}

	/* close file */
	fclose(fpSaveFile);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */
/* IMLOAD Function                                                         */
/* This function is used for loading a INTMATRIX.                          */
/* If successful, return the pointer of the matrix;                        */
/* else return NULL.                                                       */
/* ----------------------------------------------------------------------- */
struct INTMATRIX* IMLOAD(char strFilePath[])
{
	/* pMatrix is used for loading matrix */
	struct INTMATRIX* pMatrix;
	/* fpLoadFile is used for loading. */
	FILE *fpLoadFile;
	/* pElement and nElement is used for access the elements of the matrix. */
	int nElement[IM_LOAD_MAX_WIDTH];
	long nElementCount;
	/* strInLine is used for loading file. */
	char strInLine[LONG_LINE_LENGTH];

	/* init */
	pMatrix = NULL;

	/* open file */
	fpLoadFile = NULL;
	fpLoadFile = fopen(strFilePath, "rt");
	if(fpLoadFile == NULL)
	{
		printf("Can't open input file!\n");
		return NULL;
	}

	/* load the first line */
	while(fgets(strInLine, LONG_LINE_LENGTH, fpLoadFile) != NULL)
	{
		/* trim right and left */
		nElementCount = IntLoadRowVector(strInLine, nElement);
		if(nElementCount == 0)
			continue;
		if(IMADDROW(&pMatrix, nElement, nElementCount) == PROC_FAILURE)
		{
			DestroyIntMatrix(pMatrix);
			pMatrix = NULL;
			printf("Can't load int matrix!\n");
			break;
		}
	}

	/* close file */
	fclose(fpLoadFile);

	/* return */
	return pMatrix;
}

/* ----------------------------------------------------------------------- */
/* IntLoadRowVector                                                        */
/* This function is used for loading a int row vector from a string.       */
/* If successful, return the count of loaded elements;                     */
/* else return 0.                                                          */
/* ----------------------------------------------------------------------- */
long IntLoadRowVector(char strInLine[], int nElement[])
{
	/* nCount is count of loaded elements. */
	long nCount;
	/* startp and endp are pointers to the string */
	char *startp;
	char *endp;
	double dTemp;
 
	/* init */
	nCount = 0;
	StrTrimRight(strInLine);
	StrTrimLeft(strInLine);
	if(strInLine[0] == '\0')
		return nCount;

	/* load */
	startp = strInLine;
	endp = strpbrk(startp, WORD_SEPARATORS);

	while(endp != NULL)
	{
		*endp = '\0';
		dTemp = atof(startp);
		nElement[nCount] = (int)dTemp;
		nCount++;
		startp = endp+1;
		StrTrimLeft(startp);
		endp = strpbrk(startp, WORD_SEPARATORS);
	}

	dTemp = atof(startp);
	nElement[nCount] = (int)dTemp;
	nCount++;

	/* return */
	return nCount;
}

/* ----------------------------------------------------------------------- */
/* IMADDROW Function                                                       */
/* This function is used for adding a row to INTMATRIX.                    */
/* If successful, return PROC_SUCCESS;                                     */
/* else return PROC_FAILURE.                                               */
/* ----------------------------------------------------------------------- */
int IMADDROW(struct INTMATRIX **pIntMatrix, int nElement[], long nElementCount)
{
	/* define */
	int *pElement;
	long lnSize;

	/* check parameter */
	if(pIntMatrix == NULL)
		return PROC_FAILURE;

	/* if pIntMatrix doesn't exist */
	if(*pIntMatrix == NULL)
	{
		*pIntMatrix = CreateIntMatrix(1, nElementCount);
		if(*pIntMatrix == NULL)
			return PROC_FAILURE;
		
		pElement = (*pIntMatrix)->pMatElement;
		memcpy(pElement, nElement, nElementCount*sizeof(int));

		return PROC_SUCCESS; 
	}

	/* if pIntMatrix does exist */
	if((*pIntMatrix)->nWidth != nElementCount)
		return PROC_FAILURE;

	lnSize = (long)(((*pIntMatrix)->nHeight)*((*pIntMatrix)->nWidth));
	pElement = (*pIntMatrix)->pMatElement;
	(*pIntMatrix)->pMatElement = (int *)realloc(pElement, (lnSize+nElementCount)*sizeof(int));
	if((*pIntMatrix)->pMatElement == NULL)
		exit(EXIT_FAILURE);

	pElement = (*pIntMatrix)->pMatElement+lnSize;
	memcpy(pElement, nElement, nElementCount*sizeof(int));
	(*pIntMatrix)->nHeight += 1;

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */
/* IM_T Function                                                           */
/* This function creates a new INTMATRIX pY, and lets it equal to pX'.     */
/* !note: pX must exist .                                                  */
/* Make sure you have used CreateIntMatrix function to create it           */
/* before you call this function !                                         */
/* If the new matrix can't be created or other exceptions occur , the      */
/* function will return NULL ;                                             */
/* Else the function will translocate the given matrix and return the      */
/* pointer of the new matrix .                                             */
/* ----------------------------------------------------------------------- */
struct INTMATRIX* IM_T(struct INTMATRIX* pX)
{
	/* pY is the pointer to the new matrix. */
	struct INTMATRIX* pY;
	/* nHeight and nWidth are dimension informations of matrix. */
	long nHeight,nWidth;
	long nRow,nCol;
	int *pXelement,*pYelement;

	/* if source matrix doesn't exist , return NULL */
	if(pX==NULL)
	{
		return NULL;
	}
	/* if sorce matrix exists , create it's clone */
	/* get dimension information */
	nHeight = pX->nWidth;
	nWidth = pX->nHeight;

	/* create a new matrix */
	pY = NULL;
	pY = CreateIntMatrix(nHeight,nWidth);
	if(pY == NULL)
	{
		return NULL;
	}

	/* Y = X' */
	pYelement = pY->pMatElement;
	for(nRow=0; nRow<nHeight; nRow++)
	{
		pXelement = pX->pMatElement + nRow;
		for(nCol=0; nCol<nWidth; nCol++)
		{
			*pYelement = *pXelement;
			pYelement++;
			pXelement = pXelement + nHeight;
		}
	}

	/* return */
	return pY;
}

/* ----------------------------------------------------------------------- */
/* IMCOPY Function                                                         */
/* This function lets pDestIntMatrix = pSourceIntMatrix .                  */
/* !note: pDestIntMatrix and pSourceIntMatrix must exist .                 */
/* Make sure you have used CreateIntMatrix function to create them         */
/* before you call this function !                                         */
/* If the dimensions of the two matrix are not equal , the function will   */
/* return PROC_FAILURE;                                                    */
/* Else the function will copy pSourceIntMatrix to pDestIntMatrix          */
/* and return PROC_SUCCESS .                                               */
/* ----------------------------------------------------------------------- */
int IMCOPY(struct INTMATRIX* pDestIntMatrix,struct INTMATRIX* pSourceIntMatrix)
{
	long lnSize;
	int *pDestElement,*pSourceElement;

	/* if one of the two matrixes doesn't exist , return PROC_FAILURE */
	if((pDestIntMatrix==NULL)||(pSourceIntMatrix==NULL))
	{
		return PROC_FAILURE;
	}
	/* if both matrixes exist , do next . */
	else 
	{
		/* if the dimensions of two matrixes are equal , copy source to dest*/
		/* and return PROC_SUCCESS */
		if((pDestIntMatrix->nHeight==pSourceIntMatrix->nHeight)&&
			(pDestIntMatrix->nWidth==pSourceIntMatrix->nWidth))
		{
			pDestElement = pDestIntMatrix->pMatElement;
			pSourceElement = pSourceIntMatrix->pMatElement;
			lnSize = (long)((pSourceIntMatrix->nHeight)*(pSourceIntMatrix->nWidth));
			memcpy(pDestElement, pSourceElement, lnSize*sizeof(int));

			/*for(nRow=0;nRow<pSourceDoubleMatrix->nHeight;nRow++)
			{
				for(nCol=0;nCol<pSourceDoubleMatrix->nWidth;nCol++)
				{
					*pDestElement = *pSourceElement;
					pDestElement++;
					pSourceElement++;
				}
			}*/
			
			return PROC_SUCCESS;
		}
		/* if the dimensions of two matrixes are not equal ,*/
		/* return PROC_FAILURE */
		else
		{
			return PROC_FAILURE;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* IMCLONE Function                                                        */
/* This function creates a new INTMATRIX , and lets it equal to the        */
/* given INTMATRIX .                                                       */
/* !note: pSourceIntMatrix must exist .                                    */
/* Make sure you have used CreateIntMatrix function to create it           */
/* before you call this function !                                         */
/* If the new matrix can't be created or other exceptions occur , the      */
/* function will return NULL ;                                             */
/* Else the function will clone the given matrix and return the pointer of */
/* the new matrix .                                                        */
/* ----------------------------------------------------------------------- */
struct INTMATRIX* IMCLONE(struct INTMATRIX* pSourceIntMatrix)
{
	/* pNewMatrix is the pointer to the new matrix. */
	struct INTMATRIX* pNewMatrix;
	/* nHeight and nWidth are dimension informations of matrix. */
	long nHeight,nWidth;

	/* if source matrix doesn't exist , return NULL */
	if(pSourceIntMatrix==NULL)
	{
		return NULL;
	}
	/* if sorce matrix exists , create it's clone */
	else
	{
		/* get dimension information */
		nHeight = pSourceIntMatrix->nHeight;
		nWidth = pSourceIntMatrix->nWidth;
		/* create a new matrix */
		pNewMatrix = CreateIntMatrix(nHeight,nWidth);
		/* copy the source to dest */
		if(IMCOPY(pNewMatrix,pSourceIntMatrix)==PROC_SUCCESS)
		{
			return pNewMatrix;
		}
		else
		{
			DestroyIntMatrix(pNewMatrix);
			return NULL;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* IMSORTMERGEA_0 Function                                                 */
/* This function is used for sorting rows of a int matrix using            */
/* ameliorated two-way merge algorithm.                                    */
/* The function will return PROC_SUCCESS;                                  */
/* if any error, return PROC_FAILURE.                                      */
/* ----------------------------------------------------------------------- */
int IMSORTMERGEA_0(struct INTMATRIX* pSourceMatrix, struct INTMATRIX** pSortMatrix, struct LONGMATRIX** pSortIndexMatrix)
{
	/* dStart is used for recording start point */
	int *arraystart;
	int *sortarraystart;
	long *indexarray;
	long *nidx;
	long *sortindexstart;
	/* nRow is used for recording current row */
	long nRow,nCol;

	/* init */
	indexarray = NULL;

	/* check parameters */
	if(pSortMatrix != NULL)
		*pSortMatrix = NULL;
	if(pSortIndexMatrix != NULL)
		*pSortIndexMatrix = NULL;

	if(pSourceMatrix == NULL)
		return PROC_SUCCESS;

	if(pSortMatrix == NULL)
		return PROC_FAILURE;

	/* create dest matrix */
	*pSortMatrix = CreateIntMatrix(pSourceMatrix->nHeight, pSourceMatrix->nWidth);
	if(*pSortMatrix == NULL)
	{
		printf("Can't allocate memory for sorting!\n");
		return PROC_FAILURE;
	}
	if(pSortIndexMatrix != NULL)
	{
		*pSortIndexMatrix = CreateLongMatrix(pSourceMatrix->nHeight, pSourceMatrix->nWidth);
		if(*pSortIndexMatrix == NULL)
		{
			DestroyIntMatrix(*pSortMatrix);
			*pSortMatrix = NULL;
			printf("Can't allocate memory for sorting!\n");
			return PROC_FAILURE;
		}
	}

	/* sort */
	arraystart = pSourceMatrix->pMatElement;
	sortarraystart = (*pSortMatrix)->pMatElement;
	if(pSortIndexMatrix == NULL)
	{
		for(nRow=0; nRow<pSourceMatrix->nHeight; nRow++)
		{
			if(SORT_MERGE_INT_A(arraystart, NULL, sortarraystart, NULL, 0, pSourceMatrix->nWidth-1) == -1)
			{
				DestroyIntMatrix(*pSortMatrix);
				*pSortMatrix = NULL;
				return PROC_FAILURE;
			}

			arraystart = arraystart + pSourceMatrix->nWidth;
			sortarraystart = sortarraystart + (*pSortMatrix)->nWidth;
		}
	}
	else
	{
		/* create index array */
		indexarray = (long *)calloc(pSourceMatrix->nWidth, sizeof(long));
		if(indexarray == NULL)
		{
			DestroyIntMatrix(*pSortMatrix);
			*pSortMatrix = NULL;
			DestroyLongMatrix(*pSortIndexMatrix);
			*pSortIndexMatrix = NULL;
			
			return PROC_FAILURE;
		}
		nidx = indexarray;
		for(nCol=0; nCol<pSourceMatrix->nWidth; nCol++)
		{
			*nidx = nCol;
			nidx++;
		}

		/* sort */
		sortindexstart = (*pSortIndexMatrix)->pMatElement;
		for(nRow=0; nRow<pSourceMatrix->nHeight; nRow++)
		{
			if(SORT_MERGE_INT_A(arraystart, indexarray, sortarraystart, sortindexstart, 0, pSourceMatrix->nWidth-1) == -1)
			{
				DestroyIntMatrix(*pSortMatrix);
				DestroyLongMatrix(*pSortIndexMatrix);
				*pSortMatrix = NULL;
				*pSortIndexMatrix = NULL;
				free(indexarray);
				indexarray = NULL;

				return PROC_FAILURE;
			}

			arraystart = arraystart + pSourceMatrix->nWidth;
			sortarraystart = sortarraystart + (*pSortMatrix)->nWidth;
			sortindexstart = sortindexstart + (*pSortIndexMatrix)->nWidth;
		}

		/* free index array */
		free(indexarray);
		indexarray = NULL;
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */
/* CreateLongMatrix Function                                               */
/* This function is used for creating a new LONGMATRIX .                   */
/* The size of the matrix is given by two parameters : nHeight , nWidth    */
/* In other words , a nHeight*nWidth matrix will be created.               */
/* If it fails to create new matrix , the return value will be NULL .      */
/* Else the function will return the pointer to the new LONGMATRIX .       */
/* ----------------------------------------------------------------------- */
struct LONGMATRIX* CreateLongMatrix(long nHeight, long nWidth)
{
	/* pNewMatrix is the pointer to the new matrix. */
	struct LONGMATRIX* pNewMatrix;
	/* pNewElement is the pointer to the block of matrix elements. */
	long *pNewElement;
	/* nSize is the size of matrix elements. */
	long lnSize;

	/* if nHeight<=0 or nWidth<=0 , return NULL */
	if((nHeight<=0)||(nWidth<=0))
		return NULL;
	/* allocates memory to struct LONGMATRIX. */
	pNewMatrix = NULL;
	pNewMatrix = (struct LONGMATRIX*)malloc(sizeof(struct LONGMATRIX));
	if (pNewMatrix == NULL)
	{
		printf( "Can't allocate memory for creating long matrix.\n" ); 
	}
	else
	{
		/* initialize new INTMATRIX. */
		pNewMatrix->nHeight = nHeight;
		pNewMatrix->nWidth = nWidth;

		/* calculate size of matrix elements. */
		lnSize = (long)(nHeight*nWidth);
		/* allocate memory for matrix elements. */
		pNewElement = NULL;
		pNewElement = (long *)calloc(lnSize, sizeof(long));
		if( pNewElement != NULL )
		{
			pNewMatrix->pMatElement = pNewElement;
		}
		else
		{
			printf( "Can't allocate memory for long matrix elements.\n" );
			free(pNewElement);
			free(pNewMatrix);
			pNewMatrix = NULL;
		}
	}

	/* return pNewMatrix. */
	return pNewMatrix;
}

/* ----------------------------------------------------------------------- */
/* DestroyLongMatrix Function                                              */
/* This function is used for destroying a LONGMATRIX .                     */
/* First, the memory block of elements will be freed .                     */
/* Then the memory block of the LONGMATRIX struct will be freed.           */
/* ----------------------------------------------------------------------- */
void DestroyLongMatrix(struct LONGMATRIX* pDelMatrix)
{
	/* pDelElement is the pointer to the block of matrix elements. */
	long *pDelElement;

	/* delete pDelMatrix if it exists. */
	if (pDelMatrix != NULL)
	{
		/* free memory of matrix elements. */
		pDelElement = pDelMatrix->pMatElement;
		free (pDelElement);
		/* free memory of matrix strunct. */
		free (pDelMatrix);
		pDelMatrix = NULL;
	}
} 

/* ----------------------------------------------------------------------- */
/* LMSAVE Function                                                         */
/* This function is used for saving a LONGMATRIX.                          */
/* If successful, return PROC_SUCCESS;                                     */
/* else return PROC_FAILURE.                                               */
/* ----------------------------------------------------------------------- */
int LMSAVE(struct LONGMATRIX* pLongMatrix, char strFilePath[])
{
	/* fpSaveFile is used for saving. */
	FILE *fpSaveFile;
	/* nRow, nCol is used for access the matrix */
	long nRow, nCol;
	/* pElement is used for access the elements of the matrix. */
	long *pElement;
	
	/* open file */
	fpSaveFile = NULL;
	fpSaveFile = fopen(strFilePath, "wt");
	if(fpSaveFile == NULL)
	{
		printf("Can't open output file!\n");
		return PROC_FAILURE;
	}

	/* check parameter */
	if(pLongMatrix == NULL)
		return PROC_SUCCESS;

	/* save */
	pElement = pLongMatrix->pMatElement;
	for(nRow=0; nRow<pLongMatrix->nHeight; nRow++)
	{
		for(nCol=0; nCol<pLongMatrix->nWidth; nCol++)
		{
			fprintf(fpSaveFile, "%ld ", *pElement);
			pElement++;
		}
		fprintf(fpSaveFile, "\n");
	}

	/* close file */
	fclose(fpSaveFile);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */
/* LMGETAT Function                                                        */
/* This function is used for getting element indicated by row index and    */
/* column index .                                                          */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateLongMatrix function to create it          */
/* before you call this function !                                         */
/* If the row(col) indexes is larger than or equal to the row(col) number  */
/* of the matrix or smaller than 0 or some other exceptions occur, the     */
/* function will return LM_ACCESS_VIOLATION;                               */  
/* Else the function will return matrix[row][col] .                        */
/* ----------------------------------------------------------------------- */
long LMGETAT(struct LONGMATRIX* pSourceMatrix,long nRow,long nCol)
{
	/* nHeight and nWidth are dimension informations of source matrix. */
	long nHeight,nWidth;
	/* pSourceElement is used for accessing the element of matrixes . */
	long *pSourceElement;
	/* nResult is the return value */
	long nResult;

	/* check parameter */
	if((nRow<0) || (nCol<0))
	{
		printf("Access Violation!\n");
		return LM_ACCESS_VIOLATION;
	}

	/* if source matrix doesn't exist , return LM_ACCESS_VIOLATION */
	if(pSourceMatrix==NULL)
	{
		printf("Access Violation!\n");
		return LM_ACCESS_VIOLATION;
	}
	
	/* if source matrix exists, do next */
	/* get dimension information */
	nHeight = pSourceMatrix->nHeight;
	nWidth = pSourceMatrix->nWidth;

	/* if nRow(nCol) >= nHeight(nWidth),*/
	/* return LM_ACCESS_VIOLATION */
	if((nRow>=nHeight)|| (nCol>=nWidth))
	{
		printf("Access Violation!\n");
		return LM_ACCESS_VIOLATION;
	}
	/* if 0<= nRow(nCol) <nHeight(nWidth), */
	/* return matrix[nRow][nCol] .  */
	else
	{
		pSourceElement = pSourceMatrix->pMatElement;
		nResult = *(pSourceElement+nRow*nWidth+nCol);
		/* return */
		return nResult;
	}
}

/* ----------------------------------------------------------------------- */
/* LMSETAT Function                                                        */
/* This function is used for setting element indicated by row index and    */
/* column index .                                                          */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateLongMatrix function to create it          */
/* before you call this function !                                         */
/* If the row(col) indexes is larger than or equal to the row(col) number  */
/* of the matrix or smaller than 0 or some other exceptions occur, the     */
/* function will return PROC_FAILURE;                                      */  
/* Else the function will let matrix[row][col] = nValue and                */
/* return PROC_SUCCESS.                                                    */
/* ----------------------------------------------------------------------- */
int LMSETAT(struct LONGMATRIX* pSourceMatrix,long nRow,long nCol,long nValue)
{
	/* nHeight and nWidth are dimension informations of source matrix. */
	long nHeight,nWidth;
	/* pSourceElement is used for accessing the element of matrixes . */
	long *pSourceElement;
	
	/* check parameter */
	if((nRow<0) || (nCol<0))
	{
		printf("Access Violation!\n");
		return PROC_FAILURE;
	}

	/* if source matrix doesn't exist , return PROC_FAILURE */
	if(pSourceMatrix==NULL)
	{
		printf("Access Violation!\n");
		return PROC_FAILURE;
	}

	/* if source matrix exists, do next */
	else
	{
		/* get dimension information */
		nHeight = pSourceMatrix->nHeight;
		nWidth = pSourceMatrix->nWidth;
		/* if nRow(nCol) >= nHeight(nWidth),*/
		/* return PROC_FAILURE */
		if((nRow>=nHeight) || (nCol>=nWidth))
		{
			printf("Access Violation!\n");
			return PROC_FAILURE;
		}
		/* if 0<= nRow(nCol) <nHeight(nWidth), */
		/* let matrix[nRow][nCol]=nValue and .  */
		/* return PROC_SUCCESS */
		else
		{
			pSourceElement = pSourceMatrix->pMatElement;
			*(pSourceElement+nRow*nWidth+nCol) = nValue;
			/* return */
			return PROC_SUCCESS;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* LMCOPY Function                                                         */
/* This function lets pDestLongMatrix = pSourceLongMatrix .                */
/* !note: pDestLongMatrix and pSourceLongMatrix must exist .               */
/* Make sure you have used CreateLongMatrix function to create them        */
/* before you call this function !                                         */
/* If the dimensions of the two matrix are not equal , the function will   */
/* return PROC_FAILURE;                                                    */
/* Else the function will copy pSourceLongMatrix to pDestLongMatrix        */
/* and return PROC_SUCCESS .                                               */
/* ----------------------------------------------------------------------- */
int LMCOPY(struct LONGMATRIX* pDestLongMatrix,struct LONGMATRIX* pSourceLongMatrix)
{
	long lnSize;
	long *pDestElement,*pSourceElement;

	/* if one of the two matrixes doesn't exist , return PROC_FAILURE */
	if((pDestLongMatrix==NULL)||(pSourceLongMatrix==NULL))
	{
		return PROC_FAILURE;
	}
	/* if both matrixes exist , do next . */
	else 
	{
		/* if the dimensions of two matrixes are equal , copy source to dest*/
		/* and return PROC_SUCCESS */
		if((pDestLongMatrix->nHeight==pSourceLongMatrix->nHeight)&&
			(pDestLongMatrix->nWidth==pSourceLongMatrix->nWidth))
		{
			pDestElement = pDestLongMatrix->pMatElement;
			pSourceElement = pSourceLongMatrix->pMatElement;
			lnSize = (long)((pSourceLongMatrix->nHeight)*(pSourceLongMatrix->nWidth));
			memcpy(pDestElement, pSourceElement, lnSize*sizeof(long));

			/*for(nRow=0;nRow<pSourceDoubleMatrix->nHeight;nRow++)
			{
				for(nCol=0;nCol<pSourceDoubleMatrix->nWidth;nCol++)
				{
					*pDestElement = *pSourceElement;
					pDestElement++;
					pSourceElement++;
				}
			}*/
			
			return PROC_SUCCESS;
		}
		/* if the dimensions of two matrixes are not equal ,*/
		/* return PROC_FAILURE */
		else
		{
			return PROC_FAILURE;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* LMLOAD Function                                                         */
/* This function is used for loading a LONGMATRIX.                         */
/* If successful, return the pointer of the matrix;                        */
/* else return NULL.                                                       */
/* ----------------------------------------------------------------------- */
struct LONGMATRIX* LMLOAD(char strFilePath[])
{
	/* pMatrix is used for loading matrix */
	struct LONGMATRIX* pMatrix;
	/* fpLoadFile is used for loading. */
	FILE *fpLoadFile;
	/* pElement and nElement is used for access the elements of the matrix. */
	long nElement[LM_LOAD_MAX_WIDTH];
	long nElementCount;
	/* strInLine is used for loading file. */
	char strInLine[LONG_LINE_LENGTH];

	/* init */
	pMatrix = NULL;

	/* open file */
	fpLoadFile = NULL;
	fpLoadFile = fopen(strFilePath, "rt");
	if(fpLoadFile == NULL)
	{
		printf("Can't open input file!\n");
		return NULL;
	}

	/* load the first line */
	while(fgets(strInLine, LONG_LINE_LENGTH, fpLoadFile) != NULL)
	{
		/* trim right and left */
		nElementCount = LongLoadRowVector(strInLine, nElement);
		if(nElementCount == 0)
			continue;
		if(LMADDROW(&pMatrix, nElement, nElementCount) == PROC_FAILURE)
		{
			DestroyLongMatrix(pMatrix);
			pMatrix = NULL;
			printf("Can't load long matrix!\n");
			break;
		}
	}

	/* close file */
	fclose(fpLoadFile);

	/* return */
	return pMatrix;
}

/* ----------------------------------------------------------------------- */
/* LongLoadRowVector                                                       */
/* This function is used for loading a long row vector from a string.      */
/* If successful, return the count of loaded elements;                     */
/* else return 0.                                                          */
/* ----------------------------------------------------------------------- */
long LongLoadRowVector(char strInLine[], long nElement[])
{
	/* nCount is count of loaded elements. */
	long nCount;
	/* startp and endp are pointers to the string */
	char *startp;
	char *endp;
	double dTemp;
 
	/* init */
	nCount = 0;
	StrTrimRight(strInLine);
	StrTrimLeft(strInLine);
	if(strInLine[0] == '\0')
		return nCount;

	/* load */
	startp = strInLine;
	endp = strpbrk(startp, WORD_SEPARATORS);

	while(endp != NULL)
	{
		*endp = '\0';
		dTemp = atof(startp);
		nElement[nCount] = (long)dTemp;
		nCount++;
		startp = endp+1;
		StrTrimLeft(startp);
		endp = strpbrk(startp, WORD_SEPARATORS);
	}

	dTemp = atof(startp);
	nElement[nCount] = (long)dTemp;
	nCount++;

	/* return */
	return nCount;
}

/* ----------------------------------------------------------------------- */
/* LMADDROW Function                                                       */
/* This function is used for adding a row to LONGMATRIX.                   */
/* If successful, return PROC_SUCCESS;                                     */
/* else return PROC_FAILURE.                                               */
/* ----------------------------------------------------------------------- */
int LMADDROW(struct LONGMATRIX **pLongMatrix, long nElement[], long nElementCount)
{
	/* define */
	long *pElement;
	long lnSize;

	/* check parameter */
	if(pLongMatrix == NULL)
		return PROC_FAILURE;

	/* if pLongMatrix doesn't exist */
	if(*pLongMatrix == NULL)
	{
		*pLongMatrix = CreateLongMatrix(1, nElementCount);
		if(*pLongMatrix == NULL)
			return PROC_FAILURE;
		
		pElement = (*pLongMatrix)->pMatElement;
		memcpy(pElement, nElement, nElementCount*sizeof(long));

		return PROC_SUCCESS; 
	}

	/* if pLongMatrix does exist */
	if((*pLongMatrix)->nWidth != nElementCount)
		return PROC_FAILURE;

	lnSize = (long)(((*pLongMatrix)->nHeight)*((*pLongMatrix)->nWidth));
	pElement = (*pLongMatrix)->pMatElement;
	(*pLongMatrix)->pMatElement = (long *)realloc(pElement, (lnSize+nElementCount)*sizeof(long));
	if((*pLongMatrix)->pMatElement == NULL)
		exit(EXIT_FAILURE);

	pElement = (*pLongMatrix)->pMatElement+lnSize;
	memcpy(pElement, nElement, nElementCount*sizeof(long));
	(*pLongMatrix)->nHeight += 1;

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */
/* LM_T Function                                                           */
/* This function creates a new LONGMATRIX pY, and lets it equal to pX'.    */
/* !note: pX must exist .                                                  */
/* Make sure you have used CreateLongMatrix function to create it          */
/* before you call this function !                                         */
/* If the new matrix can't be created or other exceptions occur , the      */
/* function will return NULL ;                                             */
/* Else the function will translocate the given matrix and return the      */
/* pointer of the new matrix .                                             */
/* ----------------------------------------------------------------------- */
struct LONGMATRIX* LM_T(struct LONGMATRIX* pX)
{
	/* pNewMatrix is the pointer to the new matrix. */
	struct LONGMATRIX* pY;
	/* nHeight and nWidth are dimension informations of matrix. */
	long nHeight,nWidth;
	long nRow,nCol;
	long *pXelement,*pYelement;

	/* if source matrix doesn't exist , return NULL */
	if(pX==NULL)
	{
		return NULL;
	}
	/* if sorce matrix exists , create it's clone */
	/* get dimension information */
	nHeight = pX->nWidth;
	nWidth = pX->nHeight;

	/* create a new matrix */
	pY = NULL;
	pY = CreateLongMatrix(nHeight,nWidth);
	if(pY == NULL)
	{
		return NULL;
	}

	/* Y = X' */
	pYelement = pY->pMatElement;
	for(nRow=0; nRow<nHeight; nRow++)
	{
		pXelement = pX->pMatElement + nRow;
		for(nCol=0; nCol<nWidth; nCol++)
		{
			*pYelement = *pXelement;
			pYelement++;
			pXelement = pXelement + nHeight;
		}
	}

	/* return */
	return pY;
}

/* ----------------------------------------------------------------------- */
/* CreateDoubleMatrix Function                                             */
/* This function is used for creating a new DOUBLEMATRIX .                 */
/* The size of the matrix is given by two parameters : nHeight , nWidth    */
/* In other words , a nHeight*nWidth matrix will be created.               */
/* If it fails to create new matrix , the return value will be NULL .      */
/* Else the function will return the pointer to the new DOUBLEMATRIX .     */
/* ----------------------------------------------------------------------- */
struct DOUBLEMATRIX* CreateDoubleMatrix(long nHeight, long nWidth)
{
	/* pNewMatrix is the pointer to the new matrix. */
	struct DOUBLEMATRIX* pNewMatrix;
	/* pNewElement is the pointer to the block of matrix elements. */
	double *pNewElement;
	/* nSize is the size of matrix elements. */
	long lnSize;

	/* if nHeight<=0 or nWidth<=0 , return NULL */
	if((nHeight<=0)||(nWidth<=0))
		return NULL;
	/* allocates memory to struct DOUBLEMATRIX. */
	pNewMatrix = NULL;
	pNewMatrix = (struct DOUBLEMATRIX*)malloc(sizeof(struct DOUBLEMATRIX));
	if (pNewMatrix == NULL)
	{
		printf("Insufficient memory available for creating double matrix.\n");
	}
	else
	{
		/* initialize new DOUBLEMATRIX. */
		pNewMatrix->nHeight = nHeight;
		pNewMatrix->nWidth = nWidth;

		/* calculate size of matrix elements. */
		lnSize = (long)(nHeight*nWidth);
		/* allocate memory for matrix elements. */
		pNewElement = NULL;
		pNewElement = (double *)calloc(lnSize, sizeof(double));
		if( pNewElement != NULL )
		{
			pNewMatrix->pMatElement = pNewElement;
		}
		else
		{
			printf( "Can't allocate memory for double matrix elements.\n" );
			free(pNewElement);
			free(pNewMatrix);
			pNewMatrix = NULL;
		}
	}

	/* return pNewMatrix. */
	return pNewMatrix;
}

/* ----------------------------------------------------------------------- */
/* DestroyDoubleMatrix Function                                            */
/* This function is used for destroying a DOUBLEMATRIX .                   */
/* First, the memory block of elements will be freed .                     */
/* Then the memory block of the DOUBLEMATRIX struct will be freed .        */
/* ----------------------------------------------------------------------- */
void DestroyDoubleMatrix(struct DOUBLEMATRIX* pDelMatrix)
{
	/* pDelElement is the pointer to the block of matrix elements. */
	double *pDelElement;

	/* delete pDelMatrix if it exists. */
	if (pDelMatrix != NULL)
	{
		/* free memory of matrix elements. */
		pDelElement = pDelMatrix->pMatElement;
		free(pDelElement);
		/* free memory of matrix strunct. */
		free (pDelMatrix);
		pDelMatrix = NULL;
	}
}

/* ----------------------------------------------------------------------- */
/* DMRANDU Function                                                        */
/* This function is used for creating a new rand DOUBLEMATRIX .            */
/* The size of the matrix is given by two parameters : nHeight , nWidth    */
/* In other words , a nHeight*nWidth matrix will be created, and the       */
/* elements of the matrix is uniformly distributed on (0,1)                */ 
/* If it fails to create new matrix , the return value will be NULL .      */
/* Else the function will return the pointer to the new DOUBLEMATRIX .     */
/* ----------------------------------------------------------------------- */
struct DOUBLEMATRIX* DMRANDU(long nHeight, long nWidth)
{
	/* pNewMatrix is the pointer to the new matrix. */
	struct DOUBLEMATRIX* pNewMatrix;
	/* pNewElement is the pointer to the block of matrix elements. */
	double *pNewElement;
	/* nSize is the size of matrix elements. */
	long lnSize;
	long ni,nj;

	/* if nHeight<=0 or nWidth<=0 , return NULL */
	if((nHeight<=0)||(nWidth<=0))
		return NULL;
	/* allocates memory to struct DOUBLEMATRIX. */
	pNewMatrix = NULL;
	pNewMatrix = (struct DOUBLEMATRIX*)malloc(sizeof(struct DOUBLEMATRIX));
	if (pNewMatrix == NULL)
	{
		printf("Insufficient memory available for creating rand matrix.\n");
	}
	else
	{
		/* initialize new DOUBLEMATRIX. */
		pNewMatrix->nHeight = nHeight;
		pNewMatrix->nWidth = nWidth;

		/* calculate size of matrix elements. */
		lnSize = (long)(nHeight*nWidth);
		/* allocate memory for matrix elements. */
		pNewElement = NULL;
		pNewElement = (double *)calloc(lnSize, sizeof(double));
		if( pNewElement != NULL )
		{
			pNewMatrix->pMatElement = pNewElement;
			for(ni=0; ni<nHeight; ni++)
			{
				for(nj=0; nj<nWidth; nj++)
				{
					*pNewElement = rand_u();
					pNewElement++;
				}
			}
		}
		else
		{
			printf( "Can't allocate memory for double matrix elements.\n" );
			free(pNewElement);
			free(pNewMatrix);
			pNewMatrix = NULL;
		}
	}

	/* return pNewMatrix. */
	return pNewMatrix;
}

/* ----------------------------------------------------------------------- */
/* DMRANDPRODDIRICHLET Function                                            */
/* This function is used for creating a new product dirichlet rand         */
/* DOUBLEMATRIX of parameter pParamMatrix.                                 */
/* Each row of the pParamMatrix is a parameter vector for one component of */
/* product dirichlet.                                                      */
/* If it fails to create new matrix , the return value will be NULL .      */
/* Else the function will return the pointer to the new DOUBLEMATRIX .     */
/* ----------------------------------------------------------------------- */
struct DOUBLEMATRIX* DMRANDPRODDIRICHLET(struct DOUBLEMATRIX* pParamMatrix)
{
	/* pNewMatrix is the pointer to the new matrix. */
	struct DOUBLEMATRIX* pNewMatrix;
	/* pNewElement is the pointer to the block of matrix elements. */
	double *pNewElement;
	/* nSize is the size of matrix elements. */
	long lnSize;
	double *pElement,*pParam;
	long ni;
	int nLen;


	/* if no parameters, return NULL */
	if(pParamMatrix == NULL)
		return NULL;
	/* if Height<=0 or Width<=0 , return NULL */
	if((pParamMatrix->nHeight<=0)||(pParamMatrix->nWidth<=0))
		return NULL;
	/* allocates memory to struct DOUBLEMATRIX. */
	pNewMatrix = NULL;
	pNewMatrix = (struct DOUBLEMATRIX*)malloc(sizeof(struct DOUBLEMATRIX));
	if (pNewMatrix == NULL)
	{
		printf("Insufficient memory available for creating rand matrix.\n");
	}
	else
	{
		/* initialize new DOUBLEMATRIX. */
		pNewMatrix->nHeight = pParamMatrix->nHeight;
		pNewMatrix->nWidth = pParamMatrix->nWidth;

		/* calculate size of matrix elements. */
		lnSize = (long)(pNewMatrix->nHeight*pNewMatrix->nWidth);
		/* allocate memory for matrix elements. */
		pNewElement = NULL;
		pNewElement = (double *)calloc(lnSize, sizeof(double));
		if( pNewElement != NULL )
		{
			pNewMatrix->pMatElement = pNewElement;
			
			/* create dirichlet random variables */
			nLen = pParamMatrix->nWidth;
			pParam = pParamMatrix->pMatElement;
			pElement = pNewMatrix->pMatElement;
			for(ni=0; ni<pParamMatrix->nHeight; ni++)
			{
				dirichletrnd(pParam, nLen, pElement);
				pParam += nLen;
				pElement += nLen;
			}
		}
		else
		{
			printf( "Can't allocate memory for double matrix elements.\n" );
			free(pNewElement);
			free(pNewMatrix);
			pNewMatrix = NULL;
		}
	}

	/* return pNewMatrix. */
	return pNewMatrix;
}

/* ----------------------------------------------------------------------- */
/* DMPRODDIRICHLETRND Function                                             */
/* This function is used for creating product dirichlet random numbers of  */
/* of parameter pParamMatrix, and assigning the value to pOutMatrix.       */
/* Each row of the pParamMatrix is a parameter vector for one component of */
/* product dirichlet.                                                      */
/* If it fails to create new matrix, the return value will be PROC_FAILURE.*/
/* Else the function will return PROC_SUCCESS.                             */
/* ----------------------------------------------------------------------- */
int DMPRODDIRICHLETRND(struct DOUBLEMATRIX* pOutMatrix, struct DOUBLEMATRIX* pParamMatrix)
{
	double *pElement,*pParam;
	long ni;
	int nLen;

	/* if one of the matricies does not exist, return PROC_FAILURE */
	if( (pOutMatrix==NULL) || (pParamMatrix==NULL) )
	{
		return PROC_FAILURE;
	}
	/* if dimension does not match, return PROC_FAILURE */
	if( (pOutMatrix->nHeight != pParamMatrix->nHeight) || (pOutMatrix->nWidth != pParamMatrix->nWidth) )
	{
		return PROC_FAILURE;
	}

	/* create dirichlet random variables */
	nLen = pParamMatrix->nWidth;
	pParam = pParamMatrix->pMatElement;
	pElement = pOutMatrix->pMatElement;
	for(ni=0; ni<pParamMatrix->nHeight; ni++)
	{
		dirichletrnd(pParam, nLen, pElement);
		pParam += nLen;
		pElement += nLen;
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */
/* DMRANDNORM Function                                                     */
/* This function is used for creating a new normal rand DOUBLEMATRIX .     */
/* The size of the matrix is given by two parameters : nHeight , nWidth    */
/* In other words , a nHeight*nWidth matrix will be created, and the       */
/* elements of the matrix is uniformly distributed on (0,1)                */ 
/* If it fails to create new matrix , the return value will be NULL .      */
/* Else the function will return the pointer to the new DOUBLEMATRIX .     */
/* ----------------------------------------------------------------------- */
struct DOUBLEMATRIX* DMRANDNORM(long nHeight, long nWidth, double mu, double sigma)
{
	/* pNewMatrix is the pointer to the new matrix. */
	struct DOUBLEMATRIX* pNewMatrix;
	/* pNewElement is the pointer to the block of matrix elements. */
	double *pNewElement;
	/* nSize is the size of matrix elements. */
	long lnSize;
	long ni,nj;

	/* if nHeight<=0 or nWidth<=0 , return NULL */
	if((nHeight<=0)||(nWidth<=0))
		return NULL;
	/* allocates memory to struct DOUBLEMATRIX. */
	pNewMatrix = NULL;
	pNewMatrix = (struct DOUBLEMATRIX*)malloc(sizeof(struct DOUBLEMATRIX));
	if (pNewMatrix == NULL)
	{
		printf("Insufficient memory available for creating rand matrix.\n");
	}
	else
	{
		/* initialize new DOUBLEMATRIX. */
		pNewMatrix->nHeight = nHeight;
		pNewMatrix->nWidth = nWidth;

		/* calculate size of matrix elements. */
		lnSize = (long)(nHeight*nWidth);
		/* allocate memory for matrix elements. */
		pNewElement = NULL;
		pNewElement = (double *)calloc(lnSize, sizeof(double));
		if( pNewElement != NULL )
		{
			pNewMatrix->pMatElement = pNewElement;
			for(ni=0; ni<nHeight; ni++)
			{
				for(nj=0; nj<nWidth; nj++)
				{
					*pNewElement = normrnd(mu, sigma);
					pNewElement++;
				}
			}
		}
		else
		{
			printf( "Can't allocate memory for double matrix elements.\n" );
			free(pNewElement);
			free(pNewMatrix);
			pNewMatrix = NULL;
		}
	}

	/* return pNewMatrix. */
	return pNewMatrix;
}

/* ----------------------------------------------------------------------- */
/* DMSAVE Function                                                         */
/* This function is used for saving a DOUBLEMATRIX.                        */
/* If successful, return PROC_SUCCESS;                                     */
/* else return PROC_FAILURE.                                               */
/* ----------------------------------------------------------------------- */
int DMSAVE(struct DOUBLEMATRIX* pDoubleMatrix, char strFilePath[])
{
	/* fpSaveFile is used for saving. */
	FILE *fpSaveFile;
	/* nRow, nCol is used for access the matrix */
	long nRow, nCol;
	/* pElement is used for access the elements of the matrix. */
	double *pElement;
	
	/* open file */
	fpSaveFile = NULL;
	fpSaveFile = fopen(strFilePath, "wt");
	if(fpSaveFile == NULL)
	{
		printf("Can't open output file!\n");
		return PROC_FAILURE;
	}

	/* check parameter */
	if(pDoubleMatrix == NULL)
		return PROC_SUCCESS;

	/* save */
	pElement = pDoubleMatrix->pMatElement;
	for(nRow=0; nRow<pDoubleMatrix->nHeight; nRow++)
	{
		for(nCol=0; nCol<pDoubleMatrix->nWidth; nCol++)
		{
			fprintf(fpSaveFile, "% 9.7e ", *pElement);
			pElement++;
		}
		fprintf(fpSaveFile, "\n");
	}

	/* close file */
	fclose(fpSaveFile);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */
/* DMLOAD Function                                                         */
/* This function is used for loading a DOUBLEMATRIX.                       */
/* If successful, return the pointer of the matrix;                        */
/* else return NULL.                                                       */
/* ----------------------------------------------------------------------- */
struct DOUBLEMATRIX* DMLOAD(char strFilePath[])
{
	/* pDoubleMatrix is used for loading matrix */
	struct DOUBLEMATRIX* pDoubleMatrix;
	/* fpLoadFile is used for loading. */
	FILE *fpLoadFile;
	/* pElement and dElement is used for access the elements of the matrix. */
	double dElement[DM_LOAD_MAX_WIDTH];
	long nElementCount;
	/* strInLine is used for loading file. */
	char strInLine[LONG_LINE_LENGTH];

	/* init */
	pDoubleMatrix = NULL;

	/* open file */
	fpLoadFile = NULL;
	fpLoadFile = fopen(strFilePath, "rt");
	if(fpLoadFile == NULL)
	{
		printf("Can't open input file!\n");
		return NULL;
	}

	/* load the first line */
	while(fgets(strInLine, LONG_LINE_LENGTH, fpLoadFile) != NULL)
	{
		/* trim right and left */
		nElementCount = DoubleLoadRowVector(strInLine, dElement);
		if(nElementCount == 0)
			continue;
		if(DMADDROW(&pDoubleMatrix, dElement, nElementCount) == PROC_FAILURE)
		{
			DestroyDoubleMatrix(pDoubleMatrix);
			pDoubleMatrix = NULL;
			printf("Can't load double matrix!\n");
			break;
		}
	}

	/* close file */
	fclose(fpLoadFile);

	/* return */
	return pDoubleMatrix;
}

/* ----------------------------------------------------------------------- */
/* DoubleLoadRowVector                                                     */
/* This function is used for loading a double row vector from a string.    */
/* If successful, return the count of loaded elements;                     */
/* else return 0.                                                          */
/* ----------------------------------------------------------------------- */
long DoubleLoadRowVector(char strInLine[], double dElement[])
{
	/* nCount is count of loaded elements. */
	long nCount;
	/* startp and endp are pointers to the string */
	char *startp;
	char *endp;
 
	/* init */
	nCount = 0;
	StrTrimRight(strInLine);
	StrTrimLeft(strInLine);
	if(strInLine[0] == '\0')
		return nCount;

	/* load */
	startp = strInLine;
	endp = strpbrk(startp, WORD_SEPARATORS);

	while(endp != NULL)
	{
		*endp = '\0';
		dElement[nCount] = atof(startp);
		nCount++;
		startp = endp+1;
		StrTrimLeft(startp);
		endp = strpbrk(startp, WORD_SEPARATORS);
	}

	dElement[nCount] = atof(startp);
	nCount++;

	/* return */
	return nCount;
}

/* ----------------------------------------------------------------------- */
/* DMADDROW Function                                                       */
/* This function is used for adding a row to DOUBLEMATRIX.                 */
/* If successful, return PROC_SUCCESS;                                     */
/* else return PROC_FAILURE.                                               */
/* ----------------------------------------------------------------------- */
int DMADDROW(struct DOUBLEMATRIX **pDoubleMatrix, double dElement[], long nElementCount)
{
	/* define */
	double *pElement;
	long lnSize;

	/* check parameter */
	if(pDoubleMatrix == NULL)
		return PROC_FAILURE;

	/* if pDoubleMatrix doesn't exist */
	if(*pDoubleMatrix == NULL)
	{
		*pDoubleMatrix = CreateDoubleMatrix(1, nElementCount);
		if(*pDoubleMatrix == NULL)
			return PROC_FAILURE;
		
		pElement = (*pDoubleMatrix)->pMatElement;
		memcpy(pElement, dElement, nElementCount*sizeof(double));

		return PROC_SUCCESS; 
	}

	/* if pDoubleMatrix does exist */
	if((*pDoubleMatrix)->nWidth != nElementCount)
		return PROC_FAILURE;

	lnSize = (long)(((*pDoubleMatrix)->nHeight)*((*pDoubleMatrix)->nWidth));
	pElement = (*pDoubleMatrix)->pMatElement;
	(*pDoubleMatrix)->pMatElement = (double *)realloc(pElement, (lnSize+nElementCount)*sizeof(double));
	if((*pDoubleMatrix)->pMatElement == NULL)
		exit(EXIT_FAILURE);

	pElement = (*pDoubleMatrix)->pMatElement+lnSize;
	memcpy(pElement, dElement, nElementCount*sizeof(double));
	(*pDoubleMatrix)->nHeight += 1;

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */
/* DMCOPY Function                                                         */
/* This function lets pDestDoubleMatrix = pSourceDoubleMatrix .            */
/* !note: pDestDoubleMatrix and pSourceDoubleMatrix must exist .           */
/* Make sure you have used CreateDoubleMatrix function to create them      */
/* before you call this function !                                         */
/* If the dimensions of the two matrix are not equal , the function will   */
/* return PROC_FAILURE;                                                    */
/* Else the function will copy pSourceDoubleMatrix to pDestDoubleMatrix    */
/* and return PROC_SUCCESS .                                               */
/* ----------------------------------------------------------------------- */
int DMCOPY(struct DOUBLEMATRIX* pDestDoubleMatrix,struct DOUBLEMATRIX* pSourceDoubleMatrix)
{
	long lnSize;
	double *pDestElement,*pSourceElement;

	/* if one of the two matrixes doesn't exist , return PROC_FAILURE */
	if((pDestDoubleMatrix==NULL)||(pSourceDoubleMatrix==NULL))
	{
		return PROC_FAILURE;
	}
	/* if both matrixes exist , do next . */
	else 
	{
		/* if the dimensions of two matrixes are equal , copy source to dest*/
		/* and return PROC_SUCCESS */
		if((pDestDoubleMatrix->nHeight==pSourceDoubleMatrix->nHeight)&&
			(pDestDoubleMatrix->nWidth==pSourceDoubleMatrix->nWidth))
		{
			pDestElement = pDestDoubleMatrix->pMatElement;
			pSourceElement = pSourceDoubleMatrix->pMatElement;
			lnSize = (long)((pSourceDoubleMatrix->nHeight)*(pSourceDoubleMatrix->nWidth));
			memcpy(pDestElement, pSourceElement, lnSize*sizeof(double));

			/*for(nRow=0;nRow<pSourceDoubleMatrix->nHeight;nRow++)
			{
				for(nCol=0;nCol<pSourceDoubleMatrix->nWidth;nCol++)
				{
					*pDestElement = *pSourceElement;
					pDestElement++;
					pSourceElement++;
				}
			}*/
			
			return PROC_SUCCESS;
		}
		/* if the dimensions of two matrixes are not equal ,*/
		/* return PROC_FAILURE */
		else
		{
			return PROC_FAILURE;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* DMCLONE Function                                                        */
/* This function creates a new DOUBLEMATRIX , and lets it equal to the     */
/* given DOUBLEMATRIX .                                                    */
/* !note: pSourceDoubleMatrix must exist .                                 */
/* Make sure you have used CreateDoubleMatrix function to create it        */
/* before you call this function !                                         */
/* If the new matrix can't be created or other exceptions occur , the      */
/* function will return NULL ;                                             */
/* Else the function will clone the given matrix and return the pointer of */
/* the new matrix .                                                        */
/* ----------------------------------------------------------------------- */
struct DOUBLEMATRIX* DMCLONE(struct DOUBLEMATRIX* pSourceDoubleMatrix)
{
	/* pNewMatrix is the pointer to the new matrix. */
	struct DOUBLEMATRIX* pNewMatrix;
	/* nHeight and nWidth are dimension informations of matrix. */
	long nHeight,nWidth;

	/* if source matrix doesn't exist , return NULL */
	if(pSourceDoubleMatrix==NULL)
	{
		return NULL;
	}
	/* if sorce matrix exists , create it's clone */
	else
	{
		/* get dimension information */
		nHeight = pSourceDoubleMatrix->nHeight;
		nWidth = pSourceDoubleMatrix->nWidth;
		/* create a new matrix */
		pNewMatrix = CreateDoubleMatrix(nHeight,nWidth);
		/* copy the source to dest */
		if(DMCOPY(pNewMatrix,pSourceDoubleMatrix)==PROC_SUCCESS)
		{
			return pNewMatrix;
		}
		else
		{
			DestroyDoubleMatrix(pNewMatrix);
			return NULL;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* DM_T Function                                                           */
/* This function creates a new DOUBLEMATRIX pY, and lets it equal to pX'.  */
/* !note: pX must exist .                                                  */
/* Make sure you have used CreateDoubleMatrix function to create it        */
/* before you call this function !                                         */
/* If the new matrix can't be created or other exceptions occur , the      */
/* function will return NULL ;                                             */
/* Else the function will translocate the given matrix and return the      */
/* pointer of the new matrix .                                             */
/* ----------------------------------------------------------------------- */
struct DOUBLEMATRIX* DM_T(struct DOUBLEMATRIX* pX)
{
	/* pNewMatrix is the pointer to the new matrix. */
	struct DOUBLEMATRIX* pY;
	/* nHeight and nWidth are dimension informations of matrix. */
	long nHeight,nWidth;
	long nRow,nCol;
	double *pXelement,*pYelement;

	/* if source matrix doesn't exist , return NULL */
	if(pX==NULL)
	{
		return NULL;
	}
	/* if sorce matrix exists , create it's clone */
	/* get dimension information */
	nHeight = pX->nWidth;
	nWidth = pX->nHeight;

	/* create a new matrix */
	pY = NULL;
	pY = CreateDoubleMatrix(nHeight,nWidth);
	if(pY == NULL)
	{
		return NULL;
	}

	/* Y = X' */
	pYelement = pY->pMatElement;
	for(nRow=0; nRow<nHeight; nRow++)
	{
		pXelement = pX->pMatElement + nRow;
		for(nCol=0; nCol<nWidth; nCol++)
		{
			*pYelement = *pXelement;
			pYelement++;
			pXelement = pXelement + nHeight;
		}
	}

	/* return */
	return pY;
}

/* ----------------------------------------------------------------------- */
/* DMADD Function                                                          */
/* This function adds two DOUBLEMATRIXs and records the result to a new    */
/* matrix .                                                                */
/* !note: pDoubleMatrix1 and pDoubleMatrix2 must exist .                   */
/* Make sure you have used CreateDoubleMatrix function to create them      */
/* before you call this function !                                         */
/* If the result matrix can't be created or the dimensions of given        */
/* matrixes are not equal , the function will return NULL ;                */
/* Else the function will let ResultMatrix=DoubleMatrix1+DoubleMatrix2 and */
/* return the pointer of ResultMatrix .                                    */
/* ----------------------------------------------------------------------- */
struct DOUBLEMATRIX* DMADD(struct DOUBLEMATRIX* pDoubleMatrix1,struct DOUBLEMATRIX* pDoubleMatrix2)
{
	/* pResultMatrix is the pointer to the result matrix. */
	struct DOUBLEMATRIX* pResultMatrix;
	/* nHeight and nWidth are dimension informations of matrix. */
	long nHeight,nWidth;
	/* pResultElement,pSourceElement1 and pSourceElement2 are used for */
	/* accessing the element of matrixes . */
	double *pResultElement,*pSourceElement1,*pSourceElement2;
	long nRow,nCol;

	/* if source matrixes don't exist , return NULL */
	if((pDoubleMatrix1==NULL)||(pDoubleMatrix2==NULL))
	{
		return NULL;
	}
	/* if sorce matrix exist , add them */
	else
	{
		/* if the dimensions of given matrixes are equal , add them */
		if((pDoubleMatrix1->nHeight==pDoubleMatrix2->nHeight)&&
			(pDoubleMatrix1->nWidth==pDoubleMatrix2->nWidth))
		{
			/* get dimension information */
			nHeight = pDoubleMatrix1->nHeight;
			nWidth = pDoubleMatrix1->nWidth;
			/* create result matrix */
			pResultMatrix = CreateDoubleMatrix(nHeight,nWidth);
			/* if result matrix can't be created , return NULL */
			if(pResultMatrix==NULL)
			{
				return NULL;
			}
			/* if result matrix has been created , let it equal to */
			/* DoubleMatrix1+DoubleMatrix2 and return it's pointer */
			else
			{
				pResultElement = pResultMatrix->pMatElement;
				pSourceElement1 = pDoubleMatrix1->pMatElement;
				pSourceElement2 = pDoubleMatrix2->pMatElement;
				for(nRow=0;nRow<nHeight;nRow++)
				{
					for(nCol=0;nCol<nWidth;nCol++)
					{
						*pResultElement = (*pSourceElement1)+
							(*pSourceElement2);
						pResultElement++;
						pSourceElement1++;
						pSourceElement2++;
					}
				}

				return pResultMatrix;
			}
		}
		/* if the dimensions of given matrixes aren't equal , return NULL */
		else
		{
			printf("Dimensions not match, while add two matrixes\n!");
			return NULL;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* DMSUB Function                                                          */
/* This function subtracts two DOUBLEMATRIXs and records the result to a   */
/* new matrix .                                                            */
/* !note: pDoubleMatrix1 and pDoubleMatrix2 must exist .                   */
/* Make sure you have used CreateDoubleMatrix function to create them      */
/* before you call this function !                                         */
/* If the result matrix can't be created or the dimensions of given        */
/* matrixes are not equal , the function will return NULL ;                */
/* Else the function will let ResultMatrix=DoubleMatrix1-DoubleMatrix2 and */
/* return the pointer of ResultMatrix .                                    */
/* ----------------------------------------------------------------------- */
struct DOUBLEMATRIX* DMSUB(struct DOUBLEMATRIX* pDoubleMatrix1,struct DOUBLEMATRIX* pDoubleMatrix2)
{
	/* pResultMatrix is the pointer to the result matrix. */
	struct DOUBLEMATRIX* pResultMatrix;
	/* nHeight and nWidth are dimension informations of matrix. */
	long nHeight,nWidth;
	/* pResultElement,pSourceElement1 and pSourceElement2 are used for */
	/* accessing the element of matrixes . */
	double *pResultElement,*pSourceElement1,*pSourceElement2;
	long nRow,nCol;

	/* if source matrixes don't exist , return NULL */
	if((pDoubleMatrix1==NULL)||(pDoubleMatrix2==NULL))
	{
		return NULL;
	}
	/* if sorce matrix exist , subtract them */
	else
	{
		/* if the dimensions of given matrixes are equal , subtract them */
		if((pDoubleMatrix1->nHeight==pDoubleMatrix2->nHeight)&&
			(pDoubleMatrix1->nWidth==pDoubleMatrix2->nWidth))
		{
			/* get dimension information */
			nHeight = pDoubleMatrix1->nHeight;
			nWidth = pDoubleMatrix1->nWidth;
			/* create result matrix */
			pResultMatrix = CreateDoubleMatrix(nHeight,nWidth);
			/* if result matrix can't be created , return NULL */
			if(pResultMatrix==NULL)
			{
				return NULL;
			}
			/* if result matrix has been created , let it equal to */
			/* DoubleMatrix1-DoubleMatrix2 and return it's pointer */
			else
			{
				pResultElement = pResultMatrix->pMatElement;
				pSourceElement1 = pDoubleMatrix1->pMatElement;
				pSourceElement2 = pDoubleMatrix2->pMatElement;
				for(nRow=0;nRow<nHeight;nRow++)
				{
					for(nCol=0;nCol<nWidth;nCol++)
					{
						*pResultElement = (*pSourceElement1)-
							(*pSourceElement2);
						pResultElement++;
						pSourceElement1++;
						pSourceElement2++;
					}
				}

				return pResultMatrix;
			}
		}
		/* if the dimensions of given matrixes aren't equal , return NULL */
		else
		{
			printf("Dimensions not match, while subtract two matrixes\n!");
			return NULL;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* DMADDTF Function                                                        */
/* This function adds two DOUBLEMATRIXs and records the result to the      */
/* first one .                                                             */
/* !note: pDoubleMatrix1 and pDoubleMatrix2 must exist .                   */
/* Make sure you have used CreateDoubleMatrix function to create them      */
/* before you call this function !                                         */
/* If the dimensions of given matrixes are not equal , the function will   */
/* return PROC_FAILURE ;                                                   */
/* Else the function will let DoubleMatrix1=DoubleMatrix1+DoubleMatrix2    */
/* and return PROC_SUCCESS .                                               */
/* ----------------------------------------------------------------------- */
int DMADDTF(struct DOUBLEMATRIX* pDoubleMatrix1,struct DOUBLEMATRIX* pDoubleMatrix2)
{
	/* nHeight and nWidth are dimension informations of matrix. */
	long nHeight,nWidth;
	/* pSourceElement1 and pSourceElement2 are used for */
	/* accessing the element of matrixes . */
	double *pSourceElement1,*pSourceElement2;
	long nRow,nCol;

	/* if source matrixes don't exist , return PROC_FAILURE */
	if((pDoubleMatrix1==NULL)||(pDoubleMatrix2==NULL))
	{
		return PROC_FAILURE;
	}
	/* if sorce matrix exist , add them */
	else
	{
		/* if the dimensions of given matrixes are equal , add them */
		if((pDoubleMatrix1->nHeight==pDoubleMatrix2->nHeight)&&
			(pDoubleMatrix1->nWidth==pDoubleMatrix2->nWidth))
		{
			/* get dimension information */
			nHeight = pDoubleMatrix1->nHeight;
			nWidth = pDoubleMatrix1->nWidth;
			
			/* let DoubleMatrix1 equal to DoubleMatrix1+DoubleMatrix2 */
			/* and return PROC_SUCCESS */
			pSourceElement1 = pDoubleMatrix1->pMatElement;
			pSourceElement2 = pDoubleMatrix2->pMatElement;
			for(nRow=0;nRow<nHeight;nRow++)
			{
				for(nCol=0;nCol<nWidth;nCol++)
				{
					(*pSourceElement1)+=(*pSourceElement2);
					pSourceElement1++;
					pSourceElement2++;
				}
			}

			return PROC_SUCCESS;
		}
		/* if the dimensions of given matrixes aren't equal , */
		/* return PROC_FAILURE */
		else
		{
			printf("Dimensions not match\n!");
			return PROC_FAILURE;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* DMPADDTS Function                                                       */ 
/* This function add a number to a DOUBLEMATRIX                            */
/* records the result to the original matrix .                             */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateDoubleMatrix function to create it        */
/* before you call this function !                                         */
/* The function will let SourceMatrix = SourceMatrix + dLamda and          */
/* return PROC_SUCCESS .                                                   */
/* ----------------------------------------------------------------------- */
int DMPADDTS(struct DOUBLEMATRIX* pSourceMatrix, double dLamda)
{
	/* nHeight and nWidth are dimension informations of matrix. */
	long nHeight,nWidth;
	/* pSourceElement is used for accessing the element of matrixes . */
	double *pSourceElement;
	long nRow,nCol;

	/* if source matrixes doesn't exist , return PROC_FAILURE */
	if(pSourceMatrix==NULL)
	{
		return PROC_FAILURE;
	}
	
	/* if source matrixes exists , let it equal to pSourceMatrix/dLamda . */
	else
	{
		/* get dimension information */
		nHeight = pSourceMatrix->nHeight;
		nWidth = pSourceMatrix->nWidth;

		/* let it equal to pSourceMatrix/dLamda and return PROC_SUCCESS . */
		pSourceElement = pSourceMatrix->pMatElement;
		for(nRow=0;nRow<nHeight;nRow++)
		{
			for(nCol=0;nCol<nWidth;nCol++)
			{
				(*pSourceElement) += dLamda;
				pSourceElement++;
			}
		}

		/* return */
		return PROC_SUCCESS;
	}
}

/* ----------------------------------------------------------------------- */
/* DMSUBTF Function                                                        */
/* This function subtracts two DOUBLEMATRIXs and records the result to the */
/* first one .                                                             */
/* !note: pDoubleMatrix1 and pDoubleMatrix2 must exist .                   */
/* Make sure you have used CreateDoubleMatrix function to create them      */
/* before you call this function !                                         */
/* If the dimensions of given matrixes are not equal , the function will   */
/* return PROC_FAILURE ;                                                   */
/* Else the function will let DoubleMatrix1=DoubleMatrix1-DoubleMatrix2    */
/* and return PROC_SUCCESS .                                               */
/* ----------------------------------------------------------------------- */
int DMSUBTF(struct DOUBLEMATRIX* pDoubleMatrix1,struct DOUBLEMATRIX* pDoubleMatrix2)
{
	/* nHeight and nWidth are dimension informations of matrix. */
	long nHeight,nWidth;
	/* pSourceElement1 and pSourceElement2 are used for */
	/* accessing the element of matrixes . */
	double *pSourceElement1,*pSourceElement2;
	long nRow,nCol;

	/* if source matrixes don't exist , return PROC_FAILURE */
	if((pDoubleMatrix1==NULL)||(pDoubleMatrix2==NULL))
	{
		return PROC_FAILURE;
	}
	/* if sorce matrix exist , subtract them */
	else
	{
		/* if the dimensions of given matrixes are equal , subtract them */
		if((pDoubleMatrix1->nHeight==pDoubleMatrix2->nHeight)&&
			(pDoubleMatrix1->nWidth==pDoubleMatrix2->nWidth))
		{
			/* get dimension information */
			nHeight = pDoubleMatrix1->nHeight;
			nWidth = pDoubleMatrix1->nWidth;
			
			/* let DoubleMatrix1 equal to DoubleMatrix1-DoubleMatrix2 */
			/* and return PROC_SUCCESS */
			pSourceElement1 = pDoubleMatrix1->pMatElement;
			pSourceElement2 = pDoubleMatrix2->pMatElement;
			for(nRow=0;nRow<nHeight;nRow++)
			{
				for(nCol=0;nCol<nWidth;nCol++)
				{
					(*pSourceElement1)-=(*pSourceElement2);
					pSourceElement1++;
					pSourceElement2++;
				}
			}

			return PROC_SUCCESS;
		}
		/* if the dimensions of given matrixes aren't equal , */ 
		/* return PROC_FAILURE */
		else
		{
			printf("Dimensions not match\n!");
			return PROC_FAILURE;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* DMSUBTL Function                                                        */
/* This function subtracts two DOUBLEMATRIXs and records the result to the */
/* second one .                                                            */
/* !note: pDoubleMatrix1 and pDoubleMatrix2 must exist .                   */
/* Make sure you have used CreateDoubleMatrix function to create them      */
/* before you call this function !                                         */
/* If the dimensions of given matrixes are not equal , the function will   */
/* return PROC_FAILURE ;                                                   */
/* Else the function will let DoubleMatrix2=DoubleMatrix1-DoubleMatrix2    */
/* and return PROC_SUCCESS .                                               */
/* ----------------------------------------------------------------------- */
int DMSUBTL(struct DOUBLEMATRIX* pDoubleMatrix1,struct DOUBLEMATRIX* pDoubleMatrix2)
{
	/* nHeight and nWidth are dimension informations of matrix. */
	long nHeight,nWidth;
	/* pSourceElement1 and pSourceElement2 are used for */
	/* accessing the element of matrixes . */
	double *pSourceElement1,*pSourceElement2;
	long nRow,nCol;

	/* if source matrixes don't exist , return PROC_FAILURE */
	if((pDoubleMatrix1==NULL)||(pDoubleMatrix2==NULL))
	{
		return PROC_FAILURE;
	}
	/* if sorce matrix exist , subtract them */
	else
	{
		/* if the dimensions of given matrixes are equal , subtract them */
		if((pDoubleMatrix1->nHeight==pDoubleMatrix2->nHeight)&&
			(pDoubleMatrix1->nWidth==pDoubleMatrix2->nWidth))
		{
			/* get dimension information */
			nHeight = pDoubleMatrix1->nHeight;
			nWidth = pDoubleMatrix1->nWidth;
			
			/* let DoubleMatrix1 equal to DoubleMatrix1-DoubleMatrix2 */
			/* and return PROC_SUCCESS */
			pSourceElement1 = pDoubleMatrix1->pMatElement;
			pSourceElement2 = pDoubleMatrix2->pMatElement;
			for(nRow=0;nRow<nHeight;nRow++)
			{
				for(nCol=0;nCol<nWidth;nCol++)
				{
					(*pSourceElement2)= (*pSourceElement1)-(*pSourceElement2);
					pSourceElement1++;
					pSourceElement2++;
				}
			}

			return PROC_SUCCESS;
		}
		/* if the dimensions of given matrixes aren't equal , */ 
		/* return PROC_FAILURE */
		else
		{
			printf("Dimensions not match\n!");
			return PROC_FAILURE;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* DMPSUBTS Function                                                       */ 
/* This function subract a DOUBLEMATRIX from a number.                     */
/* records the result to the original matrix .                             */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateDoubleMatrix function to create it        */
/* before you call this function !                                         */
/* The function will let SourceMatrix = dLamda - SourceMatrix and          */
/* return PROC_SUCCESS .                                                   */
/* ----------------------------------------------------------------------- */
int DMPSUBTS(double dLamda, struct DOUBLEMATRIX* pSourceMatrix)
{
	/* nHeight and nWidth are dimension informations of matrix. */
	long nHeight,nWidth;
	/* pSourceElement is used for accessing the element of matrixes . */
	double *pSourceElement;
	long nRow,nCol;

	/* if source matrixes doesn't exist , return PROC_FAILURE */
	if(pSourceMatrix==NULL)
	{
		return PROC_FAILURE;
	}
	
	/* if source matrixes exists , let it equal to pSourceMatrix/dLamda . */
	else
	{
		/* get dimension information */
		nHeight = pSourceMatrix->nHeight;
		nWidth = pSourceMatrix->nWidth;

		/* let it equal to pSourceMatrix/dLamda and return PROC_SUCCESS . */
		pSourceElement = pSourceMatrix->pMatElement;
		for(nRow=0;nRow<nHeight;nRow++)
		{
			for(nCol=0;nCol<nWidth;nCol++)
			{
				(*pSourceElement) = dLamda-(*pSourceElement);
				pSourceElement++;
			}
		}

		/* return */
		return PROC_SUCCESS;
	}
}

/* ----------------------------------------------------------------------- */
/* DMPMUL Function                                                         */ 
/* This function multiplies a scalar quantity with a DOUBLEMATRIX and      */
/* records the result to a new DOUBLEMATRIX .                              */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateDoubleMatrix function to create it        */
/* before you call this function !                                         */
/* If the new matrix can't be created or other exceptions occur , the      */
/* function will return NULL ;                                             */
/* Else the function will let DestMatrix = dLamda * SourceMatrix and       */
/* return the pointer to DestMatrix .                                      */
/* ----------------------------------------------------------------------- */
struct DOUBLEMATRIX* DMPMUL(double dLamda, struct DOUBLEMATRIX* pSourceMatrix)
{
	/* pDestMatrix is the pointer to the new matrix . */
	struct DOUBLEMATRIX* pDestMatrix;
	/* nHeight and nWidth are dimension informations of matrix. */
	long nHeight,nWidth;
	/* pSourceElement and pDestElement are used for */
	/* accessing the element of matrixes . */
	double *pSourceElement,*pDestElement;
	long nRow,nCol;

	/* if source matrixes doesn't exist , return NULL */
	if(pSourceMatrix==NULL)
	{
		return NULL;
	}
	/* if source matrixes exists , create a new matrix and let it */
	/* equal to dLamda*pSourceMatrix . */
	else
	{
		/* get dimension information */
		nHeight = pSourceMatrix->nHeight;
		nWidth = pSourceMatrix->nWidth;

		/* create a new matrix */
		pDestMatrix = CreateDoubleMatrix(nHeight,nWidth);
		/* if dest matrix can't be created , return NULL */
		if(pDestMatrix==NULL)
		{
			return NULL;
		}
		/* if dest matrix has been created , let it equal to */
		/* dLamda*pSourceMatrix and return it's pointer . */
		else
		{
			pSourceElement = pSourceMatrix->pMatElement;
			pDestElement = pDestMatrix->pMatElement;
			for(nRow=0;nRow<nHeight;nRow++)
			{
				for(nCol=0;nCol<nWidth;nCol++)
				{
					*pDestElement = (*pSourceElement)*dLamda;
					pDestElement++;
					pSourceElement++;
				}
			}

			/* return */
			return pDestMatrix;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* DMPMULTS Function                                                       */ 
/* This function multiplies a scalar quantity with a DOUBLEMATRIX and      */
/* records the result to the original matrix .                             */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateDoubleMatrix function to create it        */
/* before you call this function !                                         */
/* If some exceptions occur , the function will return PROC_FAILURE ;      */
/* Else the function will let SourceMatrix = dLamda * SourceMatrix and     */
/* return PROC_SUCCESS .                                                   */
/* ----------------------------------------------------------------------- */
int DMPMULTS(double dLamda, struct DOUBLEMATRIX* pSourceMatrix)
{
	/* nHeight and nWidth are dimension informations of matrix. */
	long nHeight,nWidth;
	/* pSourceElement is used for accessing the element of matrixes . */
	double *pSourceElement;
	long nRow,nCol;

	/* if source matrixes doesn't exist , return PROC_FAILURE */
	if(pSourceMatrix==NULL)
	{
		return PROC_FAILURE;
	}
	/* if source matrixes exists , let it equal to dLamda*pSourceMatrix . */
	else
	{
		/* get dimension information */
		nHeight = pSourceMatrix->nHeight;
		nWidth = pSourceMatrix->nWidth;

		/* let it equal to dLamda*pSourceMatrix and return PROC_SUCCESS . */
		pSourceElement = pSourceMatrix->pMatElement;
		for(nRow=0;nRow<nHeight;nRow++)
		{
			for(nCol=0;nCol<nWidth;nCol++)
			{
				(*pSourceElement)*=dLamda;
				pSourceElement++;
			}
		}

		/* return */
		return PROC_SUCCESS;
	}
}

/* ----------------------------------------------------------------------- */
/* DMPMULTF Function                                                       */
/* This function multiply corresponding elements of two DOUBLEMATRIXs and  */
/* records the result to the first one.                                    */
/* !note: pDoubleMatrix1 and pDoubleMatrix2 must exist .                   */
/* Make sure you have used CreateDoubleMatrix function to create them      */
/* before you call this function !                                         */
/* If the dimensions of given matrixes are not equal , the function will   */
/* return PROC_FAILURE ;                                                   */
/* Else the function will let DoubleMatrix1=DoubleMatrix1.*DoubleMatrix2   */
/* and return PROC_SUCCESS .                                               */
/* ----------------------------------------------------------------------- */
int DMPMULTF(struct DOUBLEMATRIX* pDoubleMatrix1,struct DOUBLEMATRIX* pDoubleMatrix2)
{
	/* nHeight and nWidth are dimension informations of matrix. */
	long nHeight,nWidth;
	/* pSourceElement1 and pSourceElement2 are used for */
	/* accessing the element of matrixes . */
	double *pSourceElement1,*pSourceElement2;
	long nRow,nCol;

	/* if source matrixes don't exist , return PROC_FAILURE */
	if((pDoubleMatrix1==NULL)||(pDoubleMatrix2==NULL))
	{
		return PROC_FAILURE;
	}
	/* if sorce matrix exist , add them */
	else
	{
		/* if the dimensions of given matrixes are equal , add them */
		if((pDoubleMatrix1->nHeight==pDoubleMatrix2->nHeight)&&
			(pDoubleMatrix1->nWidth==pDoubleMatrix2->nWidth))
		{
			/* get dimension information */
			nHeight = pDoubleMatrix1->nHeight;
			nWidth = pDoubleMatrix1->nWidth;
			
			/* let DoubleMatrix1 equal to DoubleMatrix1+DoubleMatrix2 */
			/* and return PROC_SUCCESS */
			pSourceElement1 = pDoubleMatrix1->pMatElement;
			pSourceElement2 = pDoubleMatrix2->pMatElement;
			for(nRow=0;nRow<nHeight;nRow++)
			{
				for(nCol=0;nCol<nWidth;nCol++)
				{
					(*pSourceElement1) = (*pSourceElement1)*(*pSourceElement2);
					pSourceElement1++;
					pSourceElement2++;
				}
			}

			return PROC_SUCCESS;
		}
		/* if the dimensions of given matrixes aren't equal , */
		/* return PROC_FAILURE */
		else
		{
			printf("Dimensions not match\n!");
			return PROC_FAILURE;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* DMPDIV Function                                                         */ 
/* This function divides a DOUBLEMATRIX by a scalar quantity and           */
/* records the result to a new DOUBLEMATRIX .                              */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateDoubleMatrix function to create it        */
/* before you call this function !                                         */
/* If the new matrix can't be created or the absolute value of the scalar  */
/* quantity is smaller than ZERO_BOUND or other exceptions occur , the     */
/* function will return NULL ;                                             */
/* Else the function will let DestMatrix = SourceMatrix /dLamda and return */
/* the pointer to DestMatrix .                                             */
/* ----------------------------------------------------------------------- */
struct DOUBLEMATRIX* DMPDIV(struct DOUBLEMATRIX* pSourceMatrix, double dLamda)
{
	/* pDestMatrix is the pointer to the new matrix . */
	struct DOUBLEMATRIX* pDestMatrix;
	/* nHeight and nWidth are dimension informations of matrix. */
	long nHeight,nWidth;
	/* pSourceElement and pDestElement are used for */
	/* accessing the element of matrixes . */
	double *pSourceElement,*pDestElement;
	long nRow,nCol;

	/* if source matrixes doesn't exist , return NULL */
	if(pSourceMatrix==NULL)
	{
		return NULL;
	}
	/* if |dLamda| < ZERO_BOUND , return NULL */
	else if(fabs(dLamda) < ZERO_BOUND)
	{
		return NULL;
	}
	/* if source matrixes exists , create a new matrix and let it */
	/* equal to pSourceMatrix/dLamda . */
	else
	{
		/* get dimension information */
		nHeight = pSourceMatrix->nHeight;
		nWidth = pSourceMatrix->nWidth;

		/* create a new matrix */
		pDestMatrix = CreateDoubleMatrix(nHeight,nWidth);
		/* if dest matrix can't be created , return NULL */
		if(pDestMatrix==NULL)
		{
			return NULL;
		}
		/* if dest matrix has been created , let it equal to */
		/* pSourceMatrix/dLamda and return it's pointer . */
		else
		{
			pSourceElement = pSourceMatrix->pMatElement;
			pDestElement = pDestMatrix->pMatElement;
			for(nRow=0;nRow<nHeight;nRow++)
			{
				for(nCol=0;nCol<nWidth;nCol++)
				{
					*pDestElement = (*pSourceElement)/dLamda;
					pDestElement++;
					pSourceElement++;
				}
			}

			/* return */
			return pDestMatrix;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* DMPDIVTS Function                                                       */ 
/* This function divides a DOUBLEMATRIX by a scalar quantity and           */
/* records the result to the original matrix .                             */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateDoubleMatrix function to create it        */
/* before you call this function !                                         */
/* If the absolute value of the scalar quantity is smaller than ZERO_BOUND */
/* or other exceptions occur , the function will return PROC_FAILURE ;     */
/* Else the function will let SourceMatrix = SourceMatrix / dLamda and     */
/* return PROC_SUCCESS .                                                   */
/* ----------------------------------------------------------------------- */
int DMPDIVTS(struct DOUBLEMATRIX* pSourceMatrix, double dLamda)
{
	/* nHeight and nWidth are dimension informations of matrix. */
	long nHeight,nWidth;
	/* pSourceElement is used for accessing the element of matrixes . */
	double *pSourceElement;
	long nRow,nCol;

	/* if source matrixes doesn't exist , return PROC_FAILURE */
	if(pSourceMatrix==NULL)
	{
		return PROC_FAILURE;
	}
	/* if |dLamda| < ZERO_BOUND , return NULL */
	else if(fabs(dLamda) < ZERO_BOUND)
	{
		return PROC_FAILURE;
	}
	/* if source matrixes exists , let it equal to pSourceMatrix/dLamda . */
	else
	{
		/* get dimension information */
		nHeight = pSourceMatrix->nHeight;
		nWidth = pSourceMatrix->nWidth;

		/* let it equal to pSourceMatrix/dLamda and return PROC_SUCCESS . */
		pSourceElement = pSourceMatrix->pMatElement;
		for(nRow=0;nRow<nHeight;nRow++)
		{
			for(nCol=0;nCol<nWidth;nCol++)
			{
				(*pSourceElement)/=dLamda;
				pSourceElement++;
			}
		}

		/* return */
		return PROC_SUCCESS;
	}
}

/* ----------------------------------------------------------------------- */
/* DMPRECTS Function                                                       */ 
/* This function divides a number by a DOUBLEMATRIX and                    */
/* records the result to the original matrix .                             */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateDoubleMatrix function to create it        */
/* before you call this function !                                         */
/* If the absolute value of the matrix is smaller than ZERO_BOUND          */
/* or other exceptions occur , the function will return PROC_FAILURE ;     */
/* Else the function will let SourceMatrix = dLamda / SourceMatrix and     */
/* return PROC_SUCCESS .                                                   */
/* ----------------------------------------------------------------------- */
int DMPRECTS(double dLamda, struct DOUBLEMATRIX* pSourceMatrix)
{
	/* nHeight and nWidth are dimension informations of matrix. */
	long nHeight,nWidth;
	/* pSourceElement is used for accessing the element of matrixes . */
	double *pSourceElement;
	long nRow,nCol;

	/* if source matrixes doesn't exist , return PROC_FAILURE */
	if(pSourceMatrix==NULL)
	{
		return PROC_FAILURE;
	}
	/* if source matrixes exists , let it equal to pSourceMatrix/dLamda . */
	else
	{
		/* get dimension information */
		nHeight = pSourceMatrix->nHeight;
		nWidth = pSourceMatrix->nWidth;

		/* let it equal to pSourceMatrix/dLamda and return PROC_SUCCESS . */
		pSourceElement = pSourceMatrix->pMatElement;
		for(nRow=0;nRow<nHeight;nRow++)
		{
			for(nCol=0;nCol<nWidth;nCol++)
			{
				/* if |denominator| < ZERO_BOUND , return NULL */
				if(fabs(*pSourceElement) < ZERO_BOUND)
				{
					printf("Error: DMPRECTS, divide by zero!\n");
					exit(EXIT_FAILURE);
				}

				(*pSourceElement) = dLamda/(*pSourceElement);
				pSourceElement++;
			}
		}

		/* return */
		return PROC_SUCCESS;
	}
}

/* ----------------------------------------------------------------------- */
/* DMPDIVTF Function                                                       */
/* This function divides two DOUBLEMATRIXs and records the result to the   */
/* first one .                                                             */
/* !note: pDoubleMatrix1 and pDoubleMatrix2 must exist .                   */
/* Make sure you have used CreateDoubleMatrix function to create them      */
/* before you call this function !                                         */
/* If the dimensions of given matrixes are not equal , the function will   */
/* return PROC_FAILURE ;                                                   */
/* Else the function will let DoubleMatrix1=DoubleMatrix1./DoubleMatrix2   */
/* and return PROC_SUCCESS .                                               */
/* ----------------------------------------------------------------------- */
int DMPDIVTF(struct DOUBLEMATRIX* pDoubleMatrix1,struct DOUBLEMATRIX* pDoubleMatrix2)
{
	/* nHeight and nWidth are dimension informations of matrix. */
	long nHeight,nWidth;
	/* pSourceElement1 and pSourceElement2 are used for */
	/* accessing the element of matrixes . */
	double *pSourceElement1,*pSourceElement2;
	long nRow,nCol;

	/* if source matrixes don't exist , return PROC_FAILURE */
	if((pDoubleMatrix1==NULL)||(pDoubleMatrix2==NULL))
	{
		return PROC_FAILURE;
	}
	/* if sorce matrix exist , subtract them */
	else
	{
		/* if the dimensions of given matrixes are equal , subtract them */
		if((pDoubleMatrix1->nHeight==pDoubleMatrix2->nHeight)&&
			(pDoubleMatrix1->nWidth==pDoubleMatrix2->nWidth))
		{
			/* get dimension information */
			nHeight = pDoubleMatrix1->nHeight;
			nWidth = pDoubleMatrix1->nWidth;
			
			/* let DoubleMatrix1 equal to DoubleMatrix1-DoubleMatrix2 */
			/* and return PROC_SUCCESS */
			pSourceElement1 = pDoubleMatrix1->pMatElement;
			pSourceElement2 = pDoubleMatrix2->pMatElement;
			for(nRow=0;nRow<nHeight;nRow++)
			{
				for(nCol=0;nCol<nWidth;nCol++)
				{
					/* if |denominator| < ZERO_BOUND , return NULL */
					if(fabs(*pSourceElement2) < ZERO_BOUND)
					{
						printf("Error: divide by zero!\n");
						exit(EXIT_FAILURE);
					}

					(*pSourceElement1) /= (*pSourceElement2);
					pSourceElement1++;
					pSourceElement2++;
				}
			}

			return PROC_SUCCESS;
		}
		/* if the dimensions of given matrixes aren't equal , */ 
		/* return PROC_FAILURE */
		else
		{
			printf("Dimensions not match\n!");
			return PROC_FAILURE;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* DMPDIVTL Function                                                       */
/* This function divides two DOUBLEMATRIXs and records the result to the   */
/* second one .                                                            */
/* !note: pDoubleMatrix1 and pDoubleMatrix2 must exist .                   */
/* Make sure you have used CreateDoubleMatrix function to create them      */
/* before you call this function !                                         */
/* If the dimensions of given matrixes are not equal , the function will   */
/* return PROC_FAILURE ;                                                   */
/* Else the function will let DoubleMatrix2=DoubleMatrix1./DoubleMatrix2   */
/* and return PROC_SUCCESS .                                               */
/* ----------------------------------------------------------------------- */
int DMPDIVTL(struct DOUBLEMATRIX* pDoubleMatrix1,struct DOUBLEMATRIX* pDoubleMatrix2)
{
	/* nHeight and nWidth are dimension informations of matrix. */
	long nHeight,nWidth;
	/* pSourceElement1 and pSourceElement2 are used for */
	/* accessing the element of matrixes . */
	double *pSourceElement1,*pSourceElement2;
	long nRow,nCol;

	/* if source matrixes don't exist , return PROC_FAILURE */
	if((pDoubleMatrix1==NULL)||(pDoubleMatrix2==NULL))
	{
		return PROC_FAILURE;
	}
	/* if sorce matrix exist , subtract them */
	else
	{
		/* if the dimensions of given matrixes are equal , subtract them */
		if((pDoubleMatrix1->nHeight==pDoubleMatrix2->nHeight)&&
			(pDoubleMatrix1->nWidth==pDoubleMatrix2->nWidth))
		{
			/* get dimension information */
			nHeight = pDoubleMatrix1->nHeight;
			nWidth = pDoubleMatrix1->nWidth;
			
			/* let DoubleMatrix1 equal to DoubleMatrix1-DoubleMatrix2 */
			/* and return PROC_SUCCESS */
			pSourceElement1 = pDoubleMatrix1->pMatElement;
			pSourceElement2 = pDoubleMatrix2->pMatElement;
			for(nRow=0;nRow<nHeight;nRow++)
			{
				for(nCol=0;nCol<nWidth;nCol++)
				{
					/* if |denominator| < ZERO_BOUND , return NULL */
					if(fabs(*pSourceElement2) < ZERO_BOUND)
					{
						printf("Error: divide by zero!\n");
						exit(EXIT_FAILURE);
					}

					(*pSourceElement2) = (*pSourceElement1)/(*pSourceElement2);
					pSourceElement1++;
					pSourceElement2++;
				}
			}

			return PROC_SUCCESS;
		}
		/* if the dimensions of given matrixes aren't equal , */ 
		/* return PROC_FAILURE */
		else
		{
			printf("Dimensions not match\n!");
			return PROC_FAILURE;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* DMMUL Function                                                          */ 
/* This function multiplies two DOUBLEMATRIXs and records the result to a  */
/* new DOUBLEMATRIX .                                                      */
/* !note: pSourceMatrix1 and pSourceMatrix2 must exist .                   */
/* Make sure you have used CreateDoubleMatrix function to create them      */
/* before you call this function !                                         */
/* If the column number of the source matrix1 does not equal to the row    */
/* number of the source matrix2 or the new matrix can't be created or      */
/* other exceptions occur , the function will return NULL ;                */
/* Else the function will let DestMatrix = SourceMatrix1 * SourceMatrix2   */
/* and return the pointer to DestMatrix .                                  */
/* ----------------------------------------------------------------------- */
struct DOUBLEMATRIX* DMMUL(struct DOUBLEMATRIX* pSourceMatrix1,struct DOUBLEMATRIX* pSourceMatrix2)
{
	/* pDestMatrix is the pointer to the new matrix. */
	struct DOUBLEMATRIX* pDestMatrix;
	/* nHeight1 and nWidth1 are dimension informations of source matrix1. */
	long nHeight1,nWidth1;
	/* nHeight2 and nWidth2 are dimension informations of source matrix2. */
	long nHeight2,nWidth2;
	/* pSourceElement1, pSourceElement2 and pDestElement are used for */
	/* accessing the element of matrixes . */
	double *pSourceElement1,*pSourceElement2,*pDestElement;
	/* nRow1,nCol1 are counters for source matrix1. */
	long nRow1,nCol1;
	/* nCol2 is counter for source matrix2. */
	long nCol2;

	/* if source matrixes don't exist , return NULL */
	if((pSourceMatrix1==NULL)||(pSourceMatrix2==NULL))
	{
		return NULL;
	}
	/* if source matrixes exist, do next */
	else
	{
		/* get dimension information */
		nHeight1 = pSourceMatrix1->nHeight;
		nWidth1 = pSourceMatrix1->nWidth;
		nHeight2 = pSourceMatrix2->nHeight;
		nWidth2 = pSourceMatrix2->nWidth;

		/* if the column number of the source matrix1 does not */
		/* equal to the row number of the source matrix2 ,     */
		/* return NULL . */
		if(nWidth1 != nHeight2)
		{
			printf("Dimensions not match while multiplying two matrixes!\n");
			return NULL;
		}
		/* if the column number of the source matrix1 equals to */
		/* the row number of the source matrix2 , multiply them */
		/* and record the result to a new matrix . */
		else
		{
			/* create a new matrix */
			pDestMatrix = CreateDoubleMatrix(nHeight1,nWidth2);
			/* if the dest matrix can't be created , return NULL */
			if(pDestMatrix==NULL)
			{
				return NULL;
			}
			/* if the dest matrix has been created , let it equal to */
			/* sourcematrix1*sourcematrix2 and return the pointer to */
			/* the dest matrix . */
			else
			{
				pSourceElement1 = pSourceMatrix1->pMatElement;
				for(nRow1=0;nRow1<nHeight1;nRow1++)
				{
					pSourceElement2 = pSourceMatrix2->pMatElement;
					pDestElement = (pDestMatrix->pMatElement)
						+nRow1*nWidth2;
					for(nCol1=0;nCol1<nWidth1;nCol1++)
					{
						for(nCol2=0;nCol2<nWidth2;nCol2++)
						{
							(*pDestElement)+=(*pSourceElement1)*
								(*pSourceElement2);
							pDestElement++;
							pSourceElement2++;
						}
						pDestElement-=nWidth2;
						pSourceElement1++;
					}
				}

				/* return */
				return pDestMatrix;
			}
		}
	}
}

/* ----------------------------------------------------------------------- */
/* DMROWEXCH Function                                                      */
/* This function exchanges two rows of the given DOUBLEMATRIX .            */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateDoubleMatrix function to create it        */
/* before you call this function !                                         */
/* If one of the row indexes is larger than or equal to the row number of  */
/* the matrix or smaller than 0 or some other exceptions occur, the        */
/* function will return PROC_FAILURE ;                                     */  
/* Else the function will exchange the indicated rows and return           */
/* PROC_SUCCESS .                                                          */
/* ----------------------------------------------------------------------- */
int DMROWEXCH(struct DOUBLEMATRIX* pSourceMatrix,long nRow1,long nRow2)
{
	/* nHeight and nWidth are dimension informations of source matrix. */
	long nHeight,nWidth;
	/* pSourceElement1 and pSourceElement2 are used for */
	/* accessing the element of matrixes . */
	double *pSourceElement1,*pSourceElement2;
	/* dTempReg is the temp register for exchanged element */
	double dTempReg;
	long nCol;

	/* if source matrixes don't exist , return PROC_FAILURE */
	if(pSourceMatrix==NULL)
	{
		return PROC_FAILURE;
	}
	/* if source matrixes exist, do next */
	else
	{
		/* get dimension information */
		nHeight = pSourceMatrix->nHeight;
		nWidth = pSourceMatrix->nWidth;
		/* if (nRow1 or nRow2) < 0 or >=nHeight , return PROC_FAILURE */
		if((nRow1<0)||(nRow1>=nHeight)||(nRow2<0)||(nRow2>=nHeight))
		{
			printf("Invalid row index!\n");
			return PROC_FAILURE;
		}
		/* if 0<= nRow1 and nRow2 <nHeight ,exchange these two rows */
		/* and return PROC_SUCCESS . */
		else
		{
			pSourceElement1 = (pSourceMatrix->pMatElement)+
				nRow1*nWidth;
			pSourceElement2 = (pSourceMatrix->pMatElement)+
				nRow2*nWidth;
			for(nCol=0;nCol<nWidth;nCol++)
			{
				dTempReg = *pSourceElement1;
				*pSourceElement1 = *pSourceElement2;
				*pSourceElement2 = dTempReg;
				pSourceElement1++;
				pSourceElement2++;
			}

			/* return PROC_SUCCESS */
			return PROC_SUCCESS;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* DMROWPMUL Function                                                      */
/* This function multiplies the given row of a DOUBLEMATRIX with a scalar  */
/* quantity .                                                              */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateDoubleMatrix function to create it        */
/* before you call this function !                                         */
/* If the row index is larger than or equal to the row number of           */
/* the matrix or smaller than 0 or some other exceptions occur, the        */
/* function will return PROC_FAILURE ;                                     */  
/* Else the function will multiply the indicated row with dLamda and       */
/* return PROC_SUCCESS .                                                   */
/* ----------------------------------------------------------------------- */
int DMROWPMUL(struct DOUBLEMATRIX* pSourceMatrix,long nRow,double dLamda)
{
	/* nHeight and nWidth are dimension informations of source matrix. */
	long nHeight,nWidth;
	/* pSourceElement is used for accessing the element of matrixes */
	double *pSourceElement;
	long nCol;

	/* if source matrix doesn't exist , return PROC_FAILURE */
	if(pSourceMatrix==NULL)
	{
		return PROC_FAILURE;
	}
	/* if source matrix exists, do next */
	else
	{
		/* get dimension information */
		nHeight = pSourceMatrix->nHeight;
		nWidth = pSourceMatrix->nWidth;
		/* if nRow1 < 0 or >=nHeight , return PROC_FAILURE */
		if((nRow<0)||(nRow>=nHeight))
		{
			printf("Invalid row index!\n");
			return PROC_FAILURE;
		}
		/* if 0<= nRow <nHeight ,multiply the row with dLamda */
		/* and return PROC_SUCCESS . */
		else
		{
			pSourceElement = (pSourceMatrix->pMatElement)+
				nRow*nWidth;
			for(nCol=0;nCol<nWidth;nCol++)
			{
				(*pSourceElement)*=dLamda;
				pSourceElement++;
			}

			/* return PROC_SUCCESS */
			return PROC_SUCCESS;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* DMROWADDTF Function                                                     */
/* This function adds two rows of a DOUBLEMATRIX and records the result to */
/* the first row in parameter list .                                       */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateDoubleMatrix function to create it        */
/* before you call this function !                                         */
/* If one of the row indexes is larger than or equal to the row number of  */
/* the matrix or smaller than 0 or some other exceptions occur, the        */
/* function will return PROC_FAILURE ;                                     */  
/* Else the function will let row1 = row1 + row2 , and                     */
/* return PROC_SUCCESS .                                                   */
/* ----------------------------------------------------------------------- */
int DMROWADDTF(struct DOUBLEMATRIX* pSourceMatrix,long nRow1,long nRow2)
{
	/* nHeight and nWidth are dimension informations of source matrix. */
	long nHeight,nWidth;
	/* pSourceElement1 and pSourceElement2 are used for */
	/* accessing the element of matrixes . */
	double *pSourceElement1,*pSourceElement2;
	long nCol;

	/* if source matrixes don't exist , return PROC_FAILURE */
	if(pSourceMatrix==NULL)
	{
		return PROC_FAILURE;
	}
	/* if source matrixes exist, do next */
	else
	{
		/* get dimension information */
		nHeight = pSourceMatrix->nHeight;
		nWidth = pSourceMatrix->nWidth;
		/* if (nRow1 or nRow2) < 0 or >=nHeight , return PROC_FAILURE */
		if((nRow1<0)||(nRow1>=nHeight)||(nRow2<0)||(nRow2>=nHeight))
		{
			printf("Invalid row index!\n");
			return PROC_FAILURE;
		}
		/* if 0<= nRow1 and nRow2 <nHeight ,add these two rows */
		/* and return PROC_SUCCESS . */
		else
		{
			pSourceElement1 = (pSourceMatrix->pMatElement)+
				nRow1*nWidth;
			pSourceElement2 = (pSourceMatrix->pMatElement)+
				nRow2*nWidth;
			for(nCol=0;nCol<nWidth;nCol++)
			{
				(*pSourceElement1) += *pSourceElement2;
				pSourceElement1++;
				pSourceElement2++;
			}

			/* return PROC_SUCCESS */
			return PROC_SUCCESS;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* DMROWSUBTF Function                                                     */
/* This function subtracts two rows of a DOUBLEMATRIX and records the      */
/* result to the first row in parameter list .                             */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateDoubleMatrix function to create it        */
/* before you call this function !                                         */
/* If one of the row indexes is larger than or equal to the row number of  */
/* the matrix or smaller than 0 or some other exceptions occur, the        */
/* function will return PROC_FAILURE ;                                     */  
/* Else the function will let row1 = row1 - row2 , and                     */
/* return PROC_SUCCESS .                                                   */
/* ----------------------------------------------------------------------- */
int DMROWSUBTF(struct DOUBLEMATRIX* pSourceMatrix,long nRow1,long nRow2)
{
	/* nHeight and nWidth are dimension informations of source matrix. */
	long nHeight,nWidth;
	/* pSourceElement1 and pSourceElement2 are used for */
	/* accessing the element of matrixes . */
	double *pSourceElement1,*pSourceElement2;
	long nCol;

	/* if source matrixes don't exist , return PROC_FAILURE */
	if(pSourceMatrix==NULL)
	{
		return PROC_FAILURE;
	}
	/* if source matrixes exist, do next */
	else
	{
		/* get dimension information */
		nHeight = pSourceMatrix->nHeight;
		nWidth = pSourceMatrix->nWidth;
		/* if (nRow1 or nRow2) < 0 or >=nHeight , return PROC_FAILURE */
		if((nRow1<0)||(nRow1>=nHeight)||(nRow2<0)||(nRow2>=nHeight))
		{
			printf("Invalid row index!\n");
			return PROC_FAILURE;
		}
		/* if 0<= nRow1 and nRow2 <nHeight ,subtract these two rows */
		/* and return PROC_SUCCESS . */
		else
		{
			pSourceElement1 = (pSourceMatrix->pMatElement)+
				nRow1*nWidth;
			pSourceElement2 = (pSourceMatrix->pMatElement)+
				nRow2*nWidth;
			for(nCol=0;nCol<nWidth;nCol++)
			{
				(*pSourceElement1) -= *pSourceElement2;
				pSourceElement1++;
				pSourceElement2++;
			}

			/* return PROC_SUCCESS */
			return PROC_SUCCESS;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* DMROWTRNTF Function                                                     */
/* This function combines two rows of a DOUBLEMATRIX and records the       */
/* result to the first row in parameter list .                             */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateDoubleMatrix function to create it        */
/* before you call this function !                                         */
/* If one of the row indexes is larger than or equal to the row number of  */
/* the matrix or smaller than 0 or some other exceptions occur, the        */
/* function will return PROC_FAILURE ;                                     */  
/* Else the function will let row1 = row1 + dLamda2*row2 , and             */
/* return PROC_SUCCESS .                                                   */
/* ----------------------------------------------------------------------- */
int DMROWTRNTF(struct DOUBLEMATRIX* pSourceMatrix,long nRow1,double dLamda2,
			   long nRow2)
{
	/* nHeight and nWidth are dimension informations of source matrix. */
	long nHeight,nWidth;
	/* pSourceElement1 and pSourceElement2 are used for */
	/* accessing the element of matrixes . */
	double *pSourceElement1,*pSourceElement2;
	long nCol;

	/* if source matrix doesn't exist , return PROC_FAILURE */
	if(pSourceMatrix==NULL)
	{
		return PROC_FAILURE;
	}
	/* if source matrix exists, do next */
	else
	{
		/* get dimension information */
		nHeight = pSourceMatrix->nHeight;
		nWidth = pSourceMatrix->nWidth;
		/* if (nRow1 or nRow2) < 0 or >=nHeight , return PROC_FAILURE */
		if((nRow1<0)||(nRow1>=nHeight)||(nRow2<0)||(nRow2>=nHeight))
		{
			/* printf("Invalid row index!\n"); */
			return PROC_FAILURE;
		}
		/* if 0<= nRow1 and nRow2 <nHeight ,let row1=row1+dLamda2*row2 */
		/* and return PROC_SUCCESS . */
		else
		{
			pSourceElement1 = (pSourceMatrix->pMatElement)+
				nRow1*nWidth;
			pSourceElement2 = (pSourceMatrix->pMatElement)+
				nRow2*nWidth;
			for(nCol=0;nCol<nWidth;nCol++)
			{
				(*pSourceElement1) += (*pSourceElement2)*dLamda2;
				pSourceElement1++;
				pSourceElement2++;
			}

			/* return PROC_SUCCESS */
			return PROC_SUCCESS;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* DMROWNORM Function                                                      */
/* This function normalizes each row of a matrix, so that the sum of each  */
/* row is 1. The function can be used to convert a count matrix to a       */
/* Markov chain transition matrix.                                         */
/* !note: pSourceMatrix must exist .                                       */
/* If the sum of a row is 0, or some other exceptions occur, the           */
/* function will return PROC_FAILURE;                                      */
/* Else the function will return PROC_SUCCESS.                             */
/* ----------------------------------------------------------------------- */
int DMROWNORM(struct DOUBLEMATRIX* pSourceMatrix)
{
	/* pEle1 and pEle2 are used for accessing the element of matrixes. */
	double *pEle1,*pEle2;
	double dSum;
	long ni,nj;

	/* if source matrix doesn't exist , return PROC_FAILURE */
	if(pSourceMatrix==NULL)
	{
		return PROC_FAILURE;
	}
	/* if source matrix exists, do next */
	else
	{
		pEle1 = pSourceMatrix->pMatElement;
		for(ni=0; ni<pSourceMatrix->nHeight; ni++)
		{
			dSum = 0.0;
			pEle2 = pEle1;
			for(nj=0; nj<pSourceMatrix->nWidth; nj++)
			{
				dSum += *pEle2;
				pEle2++;
			}

			if( fabs(dSum) < ZERO_BOUND )
			{
				printf("Error: DMROWNORM, divide by zero since sum of row elements is equal to 0!\n"); 
				return PROC_FAILURE;
			}
			
			for(nj=0; nj<pSourceMatrix->nWidth; nj++)
			{
				*pEle1 = *pEle1/dSum;
				pEle1++;
			}
		}
	}

	/* return PROC_SUCCESS */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */
/* DMLOGTS Function                                                        */ 
/* This function takes natural logarithm of a DOUBLEMATRIX and records     */
/* the result to itself.                                                   */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateDoubleMatrix function to create it        */
/* before you call this function !                                         */
/* If log(<=0) or other exceptions occur , the function will return        */
/* PROC_FAILURE;                                                           */
/* Else the function will return PROC_SUCCESS.                             */
/* ----------------------------------------------------------------------- */
int DMLOGTS(struct DOUBLEMATRIX* pSourceMatrix)
{
	/* pSourceElement and pDestElement are used for */
	/* accessing the element of matrixes . */
	double *pSourceElement;
	long nRow,nCol;

	/* if source matrixes doesn't exist , return NULL */
	if(pSourceMatrix==NULL)
	{
		return PROC_FAILURE;
	}
	/* if source matrixes exists , create a new matrix and let it */
	/* equal to dLamda*pSourceMatrix . */
	else
	{
		pSourceElement = pSourceMatrix->pMatElement;
		for(nRow=0; nRow<pSourceMatrix->nHeight; nRow++)
		{
			for(nCol=0; nCol<pSourceMatrix->nWidth; nCol++)
			{
				if( (*pSourceElement) <= 0.0 )
				{
					printf("Error: trying to take logarithm of a non-positive number!\n");
					return PROC_FAILURE;
				}
				(*pSourceElement) = log(*pSourceElement);
				pSourceElement++;
			}
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */
/* DMPLOGTS Function                                                       */ 
/* This function takes logarithm of a DOUBLEMATRIX and records             */
/* the result to itself.                                                   */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateDoubleMatrix function to create it        */
/* before you call this function !                                         */
/* If log(<=0) or other exceptions occur , the function will return        */
/* PROC_FAILURE;                                                           */
/* Else the function will return PROC_SUCCESS.                             */
/* ----------------------------------------------------------------------- */
int DMPLOGTS(struct DOUBLEMATRIX* pSourceMatrix, double dBase)
{
	/* pSourceElement and pDestElement are used for */
	/* accessing the element of matrixes . */
	double *pSourceElement;
	long nRow,nCol;
	double dTemp;

	/* if source matrixes doesn't exist , return NULL */
	if(pSourceMatrix==NULL)
	{
		return PROC_FAILURE;
	}
	else if(dBase < 0.0)
	{
		printf("Error: DMPLOGTS, the base number in log_base() should be greater than zero!\n");
		exit(EXIT_FAILURE);
	}
	/* if source matrixes exists , create a new matrix and let it */
	/* equal to dLamda*pSourceMatrix . */
	else
	{
		dTemp = log(dBase);
		pSourceElement = pSourceMatrix->pMatElement;
		for(nRow=0; nRow<pSourceMatrix->nHeight; nRow++)
		{
			for(nCol=0; nCol<pSourceMatrix->nWidth; nCol++)
			{
				if( (*pSourceElement) <= 0.0 )
				{
					printf("Error: DMPLOGTS, trying to take logarithm of a non-positive number!\n");
					exit(EXIT_FAILURE);
				}
				(*pSourceElement) = log(*pSourceElement)/dTemp;
				pSourceElement++;
			}
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */
/* DMPPOWTS Function                                                       */ 
/* This function computes power of a DOUBLEMATRIX and records              */
/* the result to itself.                                                   */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateDoubleMatrix function to create it        */
/* before you call this function !                                         */
/* The function will compute matrix^dPow and return PROC_SUCCESS.          */
/* ----------------------------------------------------------------------- */
int DMPPOWTS(struct DOUBLEMATRIX* pSourceMatrix, double dPow)
{
	/* pSourceElement and pDestElement are used for */
	/* accessing the element of matrixes . */
	double *pSourceElement;
	long nRow,nCol;

	/* if source matrixes doesn't exist , return NULL */
	if(pSourceMatrix==NULL)
	{
		return PROC_FAILURE;
	}
	/* if source matrixes exists , create a new matrix and let it */
	/* equal to dLamda*pSourceMatrix . */
	else
	{
		pSourceElement = pSourceMatrix->pMatElement;
		for(nRow=0; nRow<pSourceMatrix->nHeight; nRow++)
		{
			for(nCol=0; nCol<pSourceMatrix->nWidth; nCol++)
			{
				(*pSourceElement) = pow((*pSourceElement),dPow);
				pSourceElement++;
			}
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */
/* DMPPOWTF Function                                                       */
/* This function compute power using two DOUBLEMATRIXs and records the     */
/* result to the first one .                                               */
/* !note: pDoubleMatrix1 and pDoubleMatrix2 must exist .                   */
/* Make sure you have used CreateDoubleMatrix function to create them      */
/* before you call this function !                                         */
/* If the dimensions of given matrixes are not equal , the function will   */
/* return PROC_FAILURE ;                                                   */
/* Else the function will let DoubleMatrix1=DoubleMatrix1.^DoubleMatrix2   */
/* and return PROC_SUCCESS .                                               */
/* ----------------------------------------------------------------------- */
int DMPPOWTF(struct DOUBLEMATRIX* pDoubleMatrix1,struct DOUBLEMATRIX* pDoubleMatrix2)
{
	/* nHeight and nWidth are dimension informations of matrix. */
	long nHeight,nWidth;
	/* pSourceElement1 and pSourceElement2 are used for */
	/* accessing the element of matrixes . */
	double *pSourceElement1,*pSourceElement2;
	long nRow,nCol;

	/* if source matrixes don't exist , return PROC_FAILURE */
	if((pDoubleMatrix1==NULL)||(pDoubleMatrix2==NULL))
	{
		return PROC_FAILURE;
	}
	/* if sorce matrix exist , subtract them */
	else
	{
		/* if the dimensions of given matrixes are equal , subtract them */
		if((pDoubleMatrix1->nHeight==pDoubleMatrix2->nHeight)&&
			(pDoubleMatrix1->nWidth==pDoubleMatrix2->nWidth))
		{
			/* get dimension information */
			nHeight = pDoubleMatrix1->nHeight;
			nWidth = pDoubleMatrix1->nWidth;
			
			/* let DoubleMatrix1 equal to DoubleMatrix1-DoubleMatrix2 */
			/* and return PROC_SUCCESS */
			pSourceElement1 = pDoubleMatrix1->pMatElement;
			pSourceElement2 = pDoubleMatrix2->pMatElement;
			for(nRow=0;nRow<nHeight;nRow++)
			{
				for(nCol=0;nCol<nWidth;nCol++)
				{
					(*pSourceElement1) = pow((*pSourceElement1),(*pSourceElement2));
					pSourceElement1++;
					pSourceElement2++;
				}
			}

			return PROC_SUCCESS;
		}
		/* if the dimensions of given matrixes aren't equal , */ 
		/* return PROC_FAILURE */
		else
		{
			printf("Dimensions not match\n!");
			return PROC_FAILURE;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* DMPPOWTL Function                                                       */
/* This function compute power using two DOUBLEMATRIXs and records the     */
/* result to the second one .                                              */
/* !note: pDoubleMatrix1 and pDoubleMatrix2 must exist .                   */
/* Make sure you have used CreateDoubleMatrix function to create them      */
/* before you call this function !                                         */
/* If the dimensions of given matrixes are not equal , the function will   */
/* return PROC_FAILURE ;                                                   */
/* Else the function will let DoubleMatrix2=DoubleMatrix1.^DoubleMatrix2   */
/* and return PROC_SUCCESS .                                               */
/* ----------------------------------------------------------------------- */
int DMPPOWTL(struct DOUBLEMATRIX* pDoubleMatrix1,struct DOUBLEMATRIX* pDoubleMatrix2)
{
	/* nHeight and nWidth are dimension informations of matrix. */
	long nHeight,nWidth;
	/* pSourceElement1 and pSourceElement2 are used for */
	/* accessing the element of matrixes . */
	double *pSourceElement1,*pSourceElement2;
	long nRow,nCol;

	/* if source matrixes don't exist , return PROC_FAILURE */
	if((pDoubleMatrix1==NULL)||(pDoubleMatrix2==NULL))
	{
		return PROC_FAILURE;
	}
	/* if sorce matrix exist , subtract them */
	else
	{
		/* if the dimensions of given matrixes are equal , subtract them */
		if((pDoubleMatrix1->nHeight==pDoubleMatrix2->nHeight)&&
			(pDoubleMatrix1->nWidth==pDoubleMatrix2->nWidth))
		{
			/* get dimension information */
			nHeight = pDoubleMatrix1->nHeight;
			nWidth = pDoubleMatrix1->nWidth;
			
			/* let DoubleMatrix1 equal to DoubleMatrix1-DoubleMatrix2 */
			/* and return PROC_SUCCESS */
			pSourceElement1 = pDoubleMatrix1->pMatElement;
			pSourceElement2 = pDoubleMatrix2->pMatElement;
			for(nRow=0;nRow<nHeight;nRow++)
			{
				for(nCol=0;nCol<nWidth;nCol++)
				{
					(*pSourceElement2) = pow((*pSourceElement1),(*pSourceElement2));
					pSourceElement1++;
					pSourceElement2++;
				}
			}

			return PROC_SUCCESS;
		}
		/* if the dimensions of given matrixes aren't equal , */ 
		/* return PROC_FAILURE */
		else
		{
			printf("Dimensions not match\n!");
			return PROC_FAILURE;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* DMPUPOWTS Function                                                      */ 
/* This function uses a DOUBLEMATRIX to compute power of a number and      */
/* records the result to the matrix itself.                                */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateDoubleMatrix function to create it        */
/* before you call this function !                                         */
/* The function will compute dBase^matrix and return PROC_SUCCESS.         */
/* ----------------------------------------------------------------------- */
int DMPUPOWTS(double dBase, struct DOUBLEMATRIX* pSourceMatrix)
{
	/* pSourceElement and pDestElement are used for */
	/* accessing the element of matrixes . */
	double *pSourceElement;
	long nRow,nCol;

	/* if source matrixes doesn't exist , return NULL */
	if(pSourceMatrix==NULL)
	{
		return PROC_FAILURE;
	}
	/* if source matrixes exists , create a new matrix and let it */
	/* equal to dLamda*pSourceMatrix . */
	else
	{
		pSourceElement = pSourceMatrix->pMatElement;
		for(nRow=0; nRow<pSourceMatrix->nHeight; nRow++)
		{
			for(nCol=0; nCol<pSourceMatrix->nWidth; nCol++)
			{
				(*pSourceElement) = pow(dBase, (*pSourceElement));
				pSourceElement++;
			}
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */
/* DMPEXPTS Function                                                       */ 
/* This function uses a DOUBLEMATRIX to compute exp(matrix) and            */
/* records the result to the matrix itself.                                */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateDoubleMatrix function to create it        */
/* before you call this function !                                         */
/* The function will compute exp(matrix) and return PROC_SUCCESS.          */
/* ----------------------------------------------------------------------- */
int DMPEXPTS(struct DOUBLEMATRIX* pSourceMatrix)
{
	/* pSourceElement and pDestElement are used for */
	/* accessing the element of matrixes . */
	double *pSourceElement;
	long nRow,nCol;

	/* if source matrixes doesn't exist , return NULL */
	if(pSourceMatrix==NULL)
	{
		return PROC_FAILURE;
	}
	/* if source matrixes exists , create a new matrix and let it */
	/* equal to dLamda*pSourceMatrix . */
	else
	{
		pSourceElement = pSourceMatrix->pMatElement;
		for(nRow=0; nRow<pSourceMatrix->nHeight; nRow++)
		{
			for(nCol=0; nCol<pSourceMatrix->nWidth; nCol++)
			{
				(*pSourceElement) = exp(*pSourceElement);
				pSourceElement++;
			}
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */
/* DMPABSTS Function                                                       */ 
/* This function computes absolute value of a matrix and                   */
/* records the result to the matrix itself.                                */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateDoubleMatrix function to create it        */
/* before you call this function !                                         */
/* The function will compute |matrix| and return PROC_SUCCESS.             */
/* ----------------------------------------------------------------------- */
int DMPABSTS(struct DOUBLEMATRIX* pSourceMatrix)
{
	/* pSourceElement and pDestElement are used for */
	/* accessing the element of matrixes . */
	double *pSourceElement;
	long nRow,nCol;

	/* if source matrixes doesn't exist , return NULL */
	if(pSourceMatrix==NULL)
	{
		return PROC_FAILURE;
	}
	/* if source matrixes exists , create a new matrix and let it */
	/* equal to dLamda*pSourceMatrix . */
	else
	{
		pSourceElement = pSourceMatrix->pMatElement;
		for(nRow=0; nRow<pSourceMatrix->nHeight; nRow++)
		{
			for(nCol=0; nCol<pSourceMatrix->nWidth; nCol++)
			{
				(*pSourceElement) = fabs(*pSourceElement);
				pSourceElement++;
			}
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */
/* DMGETAT Function                                                        */
/* This function is used for getting element indicated by row index and    */
/* column index .                                                          */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateDoubleMatrix function to create it        */
/* before you call this function !                                         */
/* If the row(col) indexes is larger than or equal to the row(col) number  */
/* of the matrix or smaller than 0 or some other exceptions occur, the     */
/* function will return DM_ACCESS_VIOLATION;                               */  
/* Else the function will return matrix[row][col] .                        */
/* ----------------------------------------------------------------------- */
double DMGETAT(struct DOUBLEMATRIX* pSourceMatrix,long nRow,long nCol)
{
	/* nHeight and nWidth are dimension informations of source matrix. */
	long nHeight,nWidth;
	/* pSourceElement is used for accessing the element of matrixes . */
	double *pSourceElement;
	/* dResult is the return value */
	double dResult;

	/* if source matrix doesn't exist , return DM_ACCESS_VIOLATION */
	if(pSourceMatrix==NULL)
	{
		printf("Access Violation!\n");
		return DM_ACCESS_VIOLATION;
	}
	/* if source matrix exists, do next */
	else
	{
		/* get dimension information */
		nHeight = pSourceMatrix->nHeight;
		nWidth = pSourceMatrix->nWidth;
		/* if nRow(nCol) < 0 or >=nHeight(nWidth),*/
		/* return DM_ACCESS_VIOLATION */
		if((nRow<0)||(nRow>=nHeight)||(nCol<0)||(nCol>=nWidth))
		{
			printf("Access Violation!\n");
			return DM_ACCESS_VIOLATION;
		}
		/* if 0<= nRow(nCol) <nHeight(nWidth), */
		/* return matrix[nRow][nCol] .  */
		else
		{
			pSourceElement = pSourceMatrix->pMatElement;
			dResult = *(pSourceElement+nRow*nWidth+nCol);
			/* return */
			return dResult;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* DMSETAT Function                                                        */
/* This function is used for setting element indicated by row index and    */
/* column index .                                                          */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateDoubleMatrix function to create it        */
/* before you call this function !                                         */
/* If the row(col) indexes is larger than or equal to the row(col) number  */
/* of the matrix or smaller than 0 or some other exceptions occur, the     */
/* function will return PROC_FAILURE;                                      */  
/* Else the function will let matrix[row][col] = dValue and                */
/* return PROC_SUCCESS.                                                    */
/* ----------------------------------------------------------------------- */
int DMSETAT(struct DOUBLEMATRIX* pSourceMatrix,long nRow,long nCol,double dValue)
{
	/* nHeight and nWidth are dimension informations of source matrix. */
	long nHeight,nWidth;
	/* pSourceElement is used for accessing the element of matrixes . */
	double *pSourceElement;
	
	/* if source matrix doesn't exist , return PROC_FAILURE */
	if(pSourceMatrix==NULL)
	{
		printf("Access Violation!\n"); 
		return PROC_FAILURE;
	}
	/* if source matrix exists, do next */
	else
	{
		/* get dimension information */
		nHeight = pSourceMatrix->nHeight;
		nWidth = pSourceMatrix->nWidth;
		/* if nRow(nCol) < 0 or >=nHeight(nWidth),*/
		/* return DM_ACCESS_VIOLATION */
		if((nRow<0)||(nRow>=nHeight)||(nCol<0)||(nCol>=nWidth))
		{
			printf("Access Violation!\n");
			return PROC_FAILURE;
		}
		/* if 0<= nRow(nCol) <nHeight(nWidth), */
		/* let matrix[nRow][nCol]=dValue and .  */
		/* return PROC_SUCCESS */
		else
		{
			pSourceElement = pSourceMatrix->pMatElement;
			*(pSourceElement+nRow*nWidth+nCol) = dValue;
			/* return */
			return PROC_SUCCESS;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* DMSOLVEEQU Funtion                                                      */
/* This function is used for solving the matrix equation A*X=B .           */
/* !note: pMatrixA and pMatrixB must exist .                               */
/* Make sure you have used CreateDoubleMatrix function to create them      */
/* before you call this function !                                         */
/* If the row number of A and B are not equal or the equation can't be     */
/* solved , the function will return NULL ;                                */
/* Else the function will return the pointer of the MatrixX which          */
/* satisfies A*X =B ;                                                      */
/* ----------------------------------------------------------------------- */
struct DOUBLEMATRIX* DMSOLVEEQU(struct DOUBLEMATRIX* pMatrixA,struct DOUBLEMATRIX* pMatrixB)
{
	/* pMA and pMX are pointers to matrix A and X. */
	struct DOUBLEMATRIX *pMA,*pMX;
	/* nHeight1 and nWidth1 are dimension informations of matrix A. */
	long nHeight1,nWidth1;
	/* nHeight2 and nWidth2 are dimension informations of matrix B. */
	long nHeight2,nWidth2;
	/* dElement1,dElement2,dLamda are used for manipulating */
	/* the matrixes */
	double dElement1,dElement2,dLamda;
	/* nRow1 and nRow2 are counters for matrixes. */
	long nRow1;
	long nRow2;
	/* nChanged serves as a state register */
	int nChanged;

	/* if source matrixes don't exist , return NULL */
	if((pMatrixA==NULL)||(pMatrixB==NULL))
	{
		return NULL;
	}
	/* if source matrixes exist, do next */
	else
	{
		/* get dimension information */
		nHeight1 = pMatrixA->nHeight;
		nWidth1 = pMatrixA->nWidth;
		nHeight2 = pMatrixB->nHeight;
		nWidth2 = pMatrixB->nWidth;

		/* if the row number of matrix A does not */
		/* equal to the row number of matrix B ,  */
		/* or the row number and column number of */
		/* matrix A are not equal ,the function will */
		/* return NULL . */
		if((nHeight1!=nHeight2)||(nHeight1!=nWidth1))
		{
			printf("Dimensions not match while multiplying two matrixes!\n");
			return NULL;
		}
		/* if the row number of matrix A equal to */
		/* the row number of matrix B ,and the row */
		/* number and column number of matrix A are */
		/* equal ,solve the equation . */
		else
		{
			/* create transition matrix pMA and result matrix */
			/* pMX . */
			pMA = DMCLONE(pMatrixA);
			pMX = DMCLONE(pMatrixB);
			/* if pMA or pMX can't be created , destroy them */
			/* and return NULL . */
			if((pMA==NULL)||(pMX==NULL))
			{
				DestroyDoubleMatrix(pMA);
				DestroyDoubleMatrix(pMX);
				return NULL;
			}
			/* if pMA or pMX have been created , solve the */
			/* equation */
			else
			{
				/* convert MA into upper triangle matrix */
				for(nRow1=0;nRow1<nHeight1;nRow1++)
				{
					dElement1 = DMGETAT(pMA,nRow1,nRow1);
					/* if dElement1=0 ,search another */
					/* row whose first Element!=0, */
					/* and exchange these two rows */
					if(fabs(dElement1)<ZERO_BOUND)
					{
						nChanged = 0;
						for(nRow2=nRow1+1;nRow2<nHeight1;nRow2++)
						{
							dElement1 = DMGETAT(pMA,nRow2,nRow1);
							if(fabs(dElement1)>=ZERO_BOUND)
							{
								nChanged = 1;
								break;
							}
						}
						if(nChanged==1)
						{
							DMROWEXCH(pMA,nRow1,nRow2);
							DMROWEXCH(pMX,nRow1,nRow2);
						}
						/* if A is a singular matrix , return NULL */
						else
						{
							/* printf("Singular Matrix!\n"); */
							DestroyDoubleMatrix(pMA);
							DestroyDoubleMatrix(pMX);
							return NULL;
						}
					}
					/* now dElement1!=0 ,let A[nRow1][nRow1]=1,*/
					/* and let A[row][nRow1]=0, where row>nRow1 */
					dLamda = 1/dElement1;
					DMROWPMUL(pMA,nRow1,dLamda);
					DMROWPMUL(pMX,nRow1,dLamda);
					for(nRow2=nRow1+1;nRow2<nHeight1;nRow2++)
					{
						dElement2 = DMGETAT(pMA,nRow2,nRow1);
						dLamda = -dElement2;
						DMROWTRNTF(pMA,nRow2,dLamda,nRow1);
						DMROWTRNTF(pMX,nRow2,dLamda,nRow1);
					}	
				}

				/* convert MA into diagonal triangle matrix */
				for(nRow1=nHeight1-1;nRow1>0;nRow1--)
				{
					for(nRow2=nRow1-1;nRow2>=0;nRow2--)
					{
						dElement2 = DMGETAT(pMA,nRow2,nRow1);
						dLamda = -dElement2;
						DMROWTRNTF(pMA,nRow2,dLamda,nRow1);
						DMROWTRNTF(pMX,nRow2,dLamda,nRow1);
					}
				}

				/* return */
				DestroyDoubleMatrix(pMA);
				return pMX;
			}
		}
	}
}

/* ----------------------------------------------------------------------- */
/* DMNORM1 Function                                                        */
/* This function calculates the 1_Norm of the DOUBLEMATRIX , in other      */
/* words , the sum of absolute values of all elements .                    */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateDoubleMatrix function to create it        */
/* before you call this function !                                         */
/* If some exceptions occur or pSourceMatrix==NULL , the function will     */
/* return 0 ;                                                              */
/* Else the function will return the sum of absolute values of all         */
/* elements .                                                              */
/* ----------------------------------------------------------------------- */
double DMNORM1(struct DOUBLEMATRIX* pSourceMatrix)
{
	/* nHeight and nWidth are dimension informations of source matrix. */
	long nHeight,nWidth;
	/* pSourceElement is used for accessing the element of matrixes . */
	double *pSourceElement;
	/* dNorm is the 1_Norm of the matrix */
	double dNorm;
	long nRow,nCol;
	
	/* Initialize */
	dNorm = 0.0;
	/* if source matrix doesn't exist , return 0 */
	if(pSourceMatrix==NULL)
	{
		printf("Access Violation!\n");
		return dNorm;
	}
	/* if source matrix exists, do next */
	else
	{
		/* get dimension information */
		nHeight = pSourceMatrix->nHeight;
		nWidth = pSourceMatrix->nWidth;
		/* calculate 1_Norm */
		pSourceElement = pSourceMatrix->pMatElement;
		for(nRow=0;nRow<nHeight;nRow++)
		{
			for(nCol=0;nCol<nWidth;nCol++)
			{
				dNorm+=fabs(*pSourceElement);
				pSourceElement++;
			}
		}

		/* return */
		return dNorm;
	}
}

/* ----------------------------------------------------------------------- */
/* DMGETMAX Function                                                       */
/* This function is used for getting the largest element of the            */
/* DOUBLEMATRIX .                                                          */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateDoubleMatrix function to create it        */
/* before you call this function !                                         */
/* If some exceptions occur or pSourceMatrix==NULL , the function will     */
/* return 0 ;                                                              */
/* Else the function will return the largest element of the matrix .       */
/* ----------------------------------------------------------------------- */
double DMGETMAX(struct DOUBLEMATRIX* pSourceMatrix)
{
	/* nHeight and nWidth are dimension informations of source matrix. */
	long nHeight,nWidth;
	/* pSourceElement is used for accessing the element of matrixes . */
	double *pSourceElement;
	/* dMax is the largest element of the matrix */
	double dMax;
	long nRow,nCol;
	
	/* Initialize */
	dMax = 0.0;
	/* if source matrix doesn't exist , return 0 */
	if(pSourceMatrix==NULL)
	{
		printf("Access Violation!\n");
		return dMax;
	}
	/* if source matrix exists, do next */
	else
	{
		/* get dimension information */
		nHeight = pSourceMatrix->nHeight;
		nWidth = pSourceMatrix->nWidth;
		/* get the largest element */
		pSourceElement = pSourceMatrix->pMatElement;
		dMax = *pSourceElement;
		for(nRow=0;nRow<nHeight;nRow++)
		{
			for(nCol=0;nCol<nWidth;nCol++)
			{
				if((*pSourceElement)>dMax)
					dMax = *pSourceElement;
				pSourceElement++;
			}
		}

		/* return */
		return dMax;
	}
}

/* ----------------------------------------------------------------------- */
/* DMGETMIN Function                                                       */
/* This function is used for getting the smallest element of the           */
/* DOUBLEMATRIX .                                                          */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateDoubleMatrix function to create it        */
/* before you call this function !                                         */
/* If some exceptions occur or pSourceMatrix==NULL , the function will     */
/* return 0 ;                                                              */
/* Else the function will return the smallest element of the matrix .      */
/* ----------------------------------------------------------------------- */
double DMGETMIN(struct DOUBLEMATRIX* pSourceMatrix)
{
	/* nHeight and nWidth are dimension informations of source matrix. */
	long nHeight,nWidth;
	/* pSourceElement is used for accessing the element of matrixes . */
	double *pSourceElement;
	/* dMin is the smallest element of the matrix */
	double dMin;
	long nRow,nCol;
	
	/* Initialize */
	dMin = 0.0;
	/* if source matrix doesn't exist , return 0 */
	if(pSourceMatrix==NULL)
	{
		printf("Access Violation!\n"); 
		return dMin;
	}
	/* if source matrix exists, do next */
	else
	{
		/* get dimension information */
		nHeight = pSourceMatrix->nHeight;
		nWidth = pSourceMatrix->nWidth;
		/* get the smallest element */
		pSourceElement = pSourceMatrix->pMatElement;
		dMin = *pSourceElement;
		for(nRow=0;nRow<nHeight;nRow++)
		{
			for(nCol=0;nCol<nWidth;nCol++)
			{
				if((*pSourceElement)<dMin)
					dMin = *pSourceElement;
				pSourceElement++;
			}
		}

		/* return */
		return dMin;
	}
}

/* ----------------------------------------------------------------------- */
/* DMCONVTOBM Function                                                     */
/* This function is used for converting DOUBLEMATRIX to BYTEMATRIX .       */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateDoubleMatrix function to create it        */
/* before you call this function !                                         */
/* If pSourceMatrix==NULL or the BYTEMATRIX can't be created or other      */
/* exceptions occur , the function will return NULL ;                      */
/* Else the function will convert source matrix to BYTEMATRIX , and        */
/* return the pointer of the BYTEMATRIX .                                  */
/* ----------------------------------------------------------------------- */
struct BYTEMATRIX* DMCONVTOBM(struct DOUBLEMATRIX* pSourceMatrix)
{
	/* nHeight and nWidth are dimension informations of source matrix. */
	long nHeight,nWidth;
	/* pSourceElement is used for accessing the element of source matrixes . */
	double *pSourceElement;
	/* pDestMatrix is used for manipulating the dest matrix */
	struct BYTEMATRIX* pDestMatrix;
	/* pDestElement is used for accessing the element of BYTEMATRIX . */
	unsigned char * pDestElement;
	/* dMax and dMin are the largest and the smallest element */
	/* of source matrix respectively */
	double dMax,dMin,dGap;
	long nRow,nCol;
	double dTemp;
	int nTemp;
	
	/* if source matrix doesn't exist , return NULL */
	if(pSourceMatrix==NULL)
	{
		printf("Access Violation!\n");
		return NULL;
	}
	/* if source matrix exists, do next */
	else
	{
		/* get dimension information */
		nHeight = pSourceMatrix->nHeight;
		nWidth = pSourceMatrix->nWidth;

		/* create the dest BYTEMATRIX */
		pDestMatrix = CreateByteMatrix(nHeight,nWidth);
		/* if dest matrix can't be created , return NULL */
		if(pDestMatrix==NULL)
		{
			return NULL;
		}
		/* if dest matrix has been created , convert source matrix */
		/* to BYTEMATRIX */
		else
		{
			/* get dMax and dMin */
			dMax = DMGETMAX(pSourceMatrix);
			dMin = DMGETMIN(pSourceMatrix);
			
			/* convert */
			dGap = dMax-dMin;
			if(dGap>0)
			{
				pSourceElement = pSourceMatrix->pMatElement;
				pDestElement = pDestMatrix->pMatElement;
				for(nRow=0;nRow<nHeight;nRow++)
				{
					for(nCol=0;nCol<nWidth;nCol++)
					{
						dTemp = ((*pSourceElement)-dMin)*255/dGap;
						nTemp = (int)(floor(dTemp));
						*pDestElement = (unsigned char)nTemp;

						pSourceElement++;
						pDestElement++;
					}
				}
			}

			/* return */
			return pDestMatrix;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* DMSORTMERGE_0 Function                                                  */
/* This function is used for sorting rows of a double matrix using two-way */
/* merge algorithm.                                                        */
/* The function will return PROC_SUCCESS;                                  */
/* if any error, return PROC_FAILURE.                                      */
/* ----------------------------------------------------------------------- */
int DMSORTMERGE_0(struct DOUBLEMATRIX* pSourceMatrix, struct DOUBLEMATRIX** pSortMatrix, struct LONGMATRIX** pSortIndexMatrix)
{
	/* sortarray and sortindexarray are used for recording sorted elements. */
	double *sortarray;
	long *sortindexarray;
	/* dStart is used for recording start point */
	double *dStart;
	long *indexarray;
	long *nidx;
	/* nRow is used for recording current row */
	long nRow,nCol;

	/* init */
	sortarray = NULL;
	sortindexarray = NULL;
	indexarray = NULL;

	/* check parameters */
	if(pSortMatrix != NULL)
		*pSortMatrix = NULL;
	if(pSortIndexMatrix != NULL)
		*pSortIndexMatrix = NULL;

	if(pSourceMatrix == NULL)
		return PROC_SUCCESS;

	if(pSortMatrix == NULL)
		return PROC_FAILURE;

	/* sort */
	dStart = pSourceMatrix->pMatElement;
	if(pSortIndexMatrix == NULL)
	{
		for(nRow=0; nRow<pSourceMatrix->nHeight; nRow++)
		{
			if(SORT_MERGE_DOUBLE(dStart, NULL, &sortarray, NULL, 0, pSourceMatrix->nWidth-1) == -1)
			{
				free(sortarray);
				DestroyDoubleMatrix(*pSortMatrix);
				*pSortMatrix = NULL;
				return PROC_FAILURE;
			}

			DMADDROW(pSortMatrix, sortarray, pSourceMatrix->nWidth);
			free(sortarray);
			sortarray = NULL;
			dStart = dStart + pSourceMatrix->nWidth;
		}
	}
	else
	{
		/* create index array */
		indexarray = (long *)calloc(pSourceMatrix->nWidth, sizeof(long));
		if(indexarray == NULL)
			return PROC_FAILURE;
		nidx = indexarray;
		for(nCol=0; nCol<pSourceMatrix->nWidth; nCol++)
		{
			*nidx = nCol;
			nidx++;
		}

		for(nRow=0; nRow<pSourceMatrix->nHeight; nRow++)
		{
			if(SORT_MERGE_DOUBLE(dStart, indexarray, &sortarray, &sortindexarray, 0, pSourceMatrix->nWidth-1) == -1)
			{
				free(sortarray);
				free(sortindexarray);
				free(indexarray);
				DestroyDoubleMatrix(*pSortMatrix);
				DestroyLongMatrix(*pSortIndexMatrix);
				*pSortMatrix = NULL;
				*pSortIndexMatrix = NULL;

				return PROC_FAILURE;
			}

			DMADDROW(pSortMatrix, sortarray, pSourceMatrix->nWidth);
			LMADDROW(pSortIndexMatrix, sortindexarray, pSourceMatrix->nWidth);
			free(sortarray);
			free(sortindexarray);
			sortarray = NULL;
			sortindexarray = NULL;
			dStart = dStart + pSourceMatrix->nWidth;
		}

		/* free index array */
		free(indexarray);
		indexarray = NULL;
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */
/* SORT_MERGE_DOUBLE Function                                              */
/* This function is used for sorting double array using two-way merge sort */
/* algorithm.                                                              */
/* Time complexity: n*log2(n)                                              */
/* The function will return the length of sortarray;                       */
/* if any error, return -1.                                                */
/* ----------------------------------------------------------------------- */
long SORT_MERGE_DOUBLE(double *array, long *indexarray, double **sortarray, long **sortindexarray, long nstart, long nend)
{
	long nlen,nlen1,nlen2;
	double *array1;
	double *array2;
	long *indexarray1;
	long *indexarray2;

	/* check parameter */
	if(sortarray != NULL)
		*sortarray = NULL;
	if(sortindexarray != NULL)
		*sortindexarray = NULL;
	
	if(array == NULL)
	{
		printf("Wrong parameters for sorting!\n");
		return -1;
	}
	
	if(sortarray == NULL)
	{
		printf("Access Violation!\n");
		return -1;
	}
	
	if((nstart < 0) || (nend < nstart))
	{
		printf("Access Violation!\n");
		return -1;
	}
	
	if((indexarray != NULL) && (sortindexarray == NULL))
	{
		printf("Access Violation!\n");
		return -1;
	}

	nlen = nend-nstart+1;

	/* allocate memory for sortarray */
	*sortarray = (double *)calloc(nlen, sizeof(double));
	if(*sortarray == NULL)
	{
		printf("Can't allocate memory for sorting!\n");
		return -1;
	}

	if(indexarray != NULL)
	{
		*sortindexarray = (long *)calloc(nlen, sizeof(long));
		if(*sortindexarray == NULL)
		{
			free(*sortarray);
			*sortarray = NULL;
			printf("Can't allocate memory for sorting!\n");
			return -1;
		}
	}

	/* sort if only one element */
	if(nstart == nend)
	{
		(*sortarray)[0] = array[nstart];
		if(indexarray != NULL)
			(*sortindexarray)[0] = indexarray[nstart];

		return nlen;
	}

	/* if there are more than one elements */
	/* create additional arrays */
	nlen1 = (long)(nlen/2);
	nlen2 = nlen - nlen1;

	/* sort */
	if(SORT_MERGE_DOUBLE(array, indexarray, &array1, &indexarray1, nstart, nstart+nlen1-1) != nlen1)
	{
		free(*sortarray);
		*sortarray = NULL;
		if(indexarray != NULL)
		{
			free(*sortindexarray);
			*sortindexarray = NULL;
		}

		return -1;
	}
	
	if(SORT_MERGE_DOUBLE(array, indexarray, &array2, &indexarray2, nstart+nlen1, nend) != nlen2)
	{
		free(*sortarray);
		*sortarray = NULL;
		if(indexarray != NULL)
		{
			free(*sortindexarray);
			*sortindexarray = NULL;
		}

		return -1;
	}

	if(indexarray == NULL)
	{
		if(SORT_MERGE_DOUBLE_MERGE_1(array1, nlen1, array2, nlen2, sortarray, nlen) == PROC_FAILURE)
		{	
			free(array1);
			free(array2);

			free(*sortarray);
			*sortarray = NULL;
			
			return -1;
		}
	}
	else
	{
		if(SORT_MERGE_DOUBLE_MERGE_2(array1, indexarray1, nlen1, array2, indexarray2, nlen2, sortarray, sortindexarray, nlen) == PROC_FAILURE)
		{	
			free(array1);
			free(array2);
			free(indexarray1);
			free(indexarray2);

			free(*sortarray);
			*sortarray = NULL;
			free(*sortindexarray);
			*sortindexarray = NULL;
			
			return -1;
		}
	}

	/* release memory */
	free(array1);
	free(array2);
	free(indexarray1);
	free(indexarray2);

	/* return */
	return nlen;
}

/* ----------------------------------------------------------------------- */
/* SORT_MERGE_DOUBLE_MERGE_1 Function                                      */
/* This function is used for merging and sorting two double arrays.        */
/* The function will return PROC_SUCCESS;                                  */
/* if any error, return PROC_FAILURE.                                      */
/* ----------------------------------------------------------------------- */
int SORT_MERGE_DOUBLE_MERGE_1(double *array1, long nlen1, double *array2, long nlen2, double **sortarray, long nlen)
{
	long i,j,k;
	
	/* check parameter */
	i = 0;
	j = 0;
	k = 0;

	while((i < nlen1) && (j < nlen2))
	{
		if(array1[i] <= array2[j])
		{
			(*sortarray)[k] = array1[i];
			i++;
		}
		else
		{
			(*sortarray)[k] = array2[j];
			j++;
		}
		k++;
	}

	if(i == nlen1)
	{
		while(j < nlen2)
		{
			(*sortarray)[k] = array2[j];
			j++;
			k++;
		}
	}
	else
	{
		while(i < nlen1)
		{
			(*sortarray)[k] = array1[i];
			i++;
			k++;
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */
/* SORT_MERGE_DOUBLE_MERGE_2 Function                                      */
/* This function is used for merging and sorting two double arrays.        */
/* The function will return PROC_SUCCESS;                                  */
/* if any error, return PROC_FAILURE.                                      */
/* ----------------------------------------------------------------------- */
int SORT_MERGE_DOUBLE_MERGE_2(double *array1, long *indexarray1, long nlen1, double *array2, long *indexarray2, long nlen2, double **sortarray, long **sortindexarray, long nlen)
{
	long i,j,k;
	
	/* check parameter */
	i = 0;
	j = 0;
	k = 0;

	while((i < nlen1) && (j < nlen2))
	{
		if(array1[i] <= array2[j])
		{
			(*sortarray)[k] = array1[i];
			(*sortindexarray)[k] = indexarray1[i];
			i++;
		}
		else
		{
			(*sortarray)[k] = array2[j];
			(*sortindexarray)[k] = indexarray2[j];
			j++;
		}
		k++;
	}

	if(i == nlen1)
	{
		while(j < nlen2)
		{
			(*sortarray)[k] = array2[j];
			(*sortindexarray)[k] = indexarray2[j];
			j++;
			k++;
		}
	}
	else
	{
		while(i < nlen1)
		{
			(*sortarray)[k] = array1[i];
			(*sortindexarray)[k] = indexarray1[i];
			i++;
			k++;
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */
/* DMSORTMERGEA_0 Function                                                 */
/* This function is used for sorting rows of a double matrix using         */
/* ameliorated two-way merge algorithm.                                    */
/* The function will return PROC_SUCCESS;                                  */
/* if any error, return PROC_FAILURE.                                      */
/* ----------------------------------------------------------------------- */
int DMSORTMERGEA_0(struct DOUBLEMATRIX* pSourceMatrix, struct DOUBLEMATRIX** pSortMatrix, struct LONGMATRIX** pSortIndexMatrix)
{
	/* dStart is used for recording start point */
	double *arraystart;
	double *sortarraystart;
	long *indexarray;
	long *nidx;
	long *sortindexstart;
	/* nRow is used for recording current row */
	long nRow,nCol;

	/* init */
	indexarray = NULL;

	/* check parameters */
	if(pSortMatrix != NULL)
		*pSortMatrix = NULL;
	if(pSortIndexMatrix != NULL)
		*pSortIndexMatrix = NULL;

	if(pSourceMatrix == NULL)
		return PROC_SUCCESS;

	if(pSortMatrix == NULL)
		return PROC_FAILURE;

	/* create dest matrix */
	*pSortMatrix = CreateDoubleMatrix(pSourceMatrix->nHeight, pSourceMatrix->nWidth);
	if(*pSortMatrix == NULL)
	{
		printf("Can't allocate memory for sorting!\n");
		return PROC_FAILURE;
	}
	if(pSortIndexMatrix != NULL)
	{
		*pSortIndexMatrix = CreateLongMatrix(pSourceMatrix->nHeight, pSourceMatrix->nWidth);
		if(*pSortIndexMatrix == NULL)
		{
			DestroyDoubleMatrix(*pSortMatrix);
			*pSortMatrix = NULL;
			printf("Can't allocate memory for sorting!\n");
			return PROC_FAILURE;
		}
	}

	/* sort */
	arraystart = pSourceMatrix->pMatElement;
	sortarraystart = (*pSortMatrix)->pMatElement;
	if(pSortIndexMatrix == NULL)
	{
		for(nRow=0; nRow<pSourceMatrix->nHeight; nRow++)
		{
			if(SORT_MERGE_DOUBLE_A(arraystart, NULL, sortarraystart, NULL, 0, pSourceMatrix->nWidth-1) == -1)
			{
				DestroyDoubleMatrix(*pSortMatrix);
				*pSortMatrix = NULL;
				return PROC_FAILURE;
			}

			arraystart = arraystart + pSourceMatrix->nWidth;
			sortarraystart = sortarraystart + (*pSortMatrix)->nWidth;
		}
	}
	else
	{
		/* create index array */
		indexarray = (long *)calloc(pSourceMatrix->nWidth, sizeof(long));
		if(indexarray == NULL)
		{
			DestroyDoubleMatrix(*pSortMatrix);
			*pSortMatrix = NULL;
			DestroyLongMatrix(*pSortIndexMatrix);
			*pSortIndexMatrix = NULL;
			
			return PROC_FAILURE;
		}
		nidx = indexarray;
		for(nCol=0; nCol<pSourceMatrix->nWidth; nCol++)
		{
			*nidx = nCol;
			nidx++;
		}

		/* sort */
		sortindexstart = (*pSortIndexMatrix)->pMatElement;
		for(nRow=0; nRow<pSourceMatrix->nHeight; nRow++)
		{
			if(SORT_MERGE_DOUBLE_A(arraystart, indexarray, sortarraystart, sortindexstart, 0, pSourceMatrix->nWidth-1) == -1)
			{
				DestroyDoubleMatrix(*pSortMatrix);
				DestroyLongMatrix(*pSortIndexMatrix);
				*pSortMatrix = NULL;
				*pSortIndexMatrix = NULL;
				free(indexarray);
				indexarray = NULL;

				return PROC_FAILURE;
			}

			arraystart = arraystart + pSourceMatrix->nWidth;
			sortarraystart = sortarraystart + (*pSortMatrix)->nWidth;
			sortindexstart = sortindexstart + (*pSortIndexMatrix)->nWidth;
		}

		/* free index array */
		free(indexarray);
		indexarray = NULL;
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */
/* SORT_MERGE_DOUBLE_A Function                                            */
/* This function is used for sorting double array using ameliorated        */
/* two-way merge sort algorithm.                                           */
/* Time complexity: n*log2(n)                                              */
/* The function will return the length of sortarray;                       */
/* if any error, return -1.                                                */
/* ----------------------------------------------------------------------- */
long SORT_MERGE_DOUBLE_A(double *array, long *indexarray, double *sortarray, long *sortindexarray, long nstart, long nend)
{
	long nlen,nlen1,nlen2;
	double *array1;
	double *array2;
	long *indexarray1;
	long *indexarray2;

	/* check parameter */
	if(array == NULL)
	{
		printf("Wrong parameters for sorting!\n");
		return -1;
	}
	
	if(sortarray == NULL)
	{
		printf("Access Violation!\n");
		return -1;
	}
	
	if((nstart < 0) || (nend < nstart))
	{
		printf("Access Violation!\n");
		return -1;
	}
	
	if((indexarray != NULL) && (sortindexarray == NULL))
	{
		printf("Access Violation!\n");
		return -1;
	}

	nlen = nend-nstart+1;

	/* sort if only one element */
	if(nstart == nend)
	{
		sortarray[0] = array[nstart];
		if(indexarray != NULL)
			sortindexarray[0] = indexarray[nstart];

		return nlen;
	}

	/* if there are more than one elements */
	/* create additional arrays */
	nlen1 = (long)(nlen/2);
	nlen2 = nlen - nlen1;

	/* allocate memory */
	array1 = NULL;
	array2 = NULL;
	indexarray1 = NULL;
	indexarray2 = NULL;

	array1 = (double *)calloc(nlen1, sizeof(double));
	array2 = (double *)calloc(nlen2, sizeof(double));
	if((array1 == NULL) || (array2 == NULL))
	{
		free(array1);
		free(array2);

		printf("Can't allocate memory for sorting!\n");
		return -1;
	}
	if(indexarray != NULL)
	{
		indexarray1 = (long *)calloc(nlen1, sizeof(long));
		indexarray2 = (long *)calloc(nlen2, sizeof(long));
		if((indexarray1 == NULL) || (indexarray2 == NULL))
		{
			free(array1);
			free(array2);
			free(indexarray1);
			free(indexarray2);

			printf("Can't allocate memory for sorting!\n");
			return -1;
		}
	}

	/* sort */
	if(SORT_MERGE_DOUBLE_A(array, indexarray, array1, indexarray1, nstart, nstart+nlen1-1) != nlen1)
	{
		free(array1);
		free(array2);
		free(indexarray1);
		free(indexarray2);

		return -1;
	}
	
	if(SORT_MERGE_DOUBLE_A(array, indexarray, array2, indexarray2, nstart+nlen1, nend) != nlen2)
	{
		free(array1);
		free(array2);
		free(indexarray1);
		free(indexarray2);

		return -1;
	}

	if(indexarray == NULL)
	{
		if(SORT_MERGE_DOUBLE_MERGE_A_1(array1, nlen1, array2, nlen2, sortarray, nlen) == PROC_FAILURE)
		{	
			free(array1);
			free(array2);

			return -1;
		}
	}
	else
	{
		if(SORT_MERGE_DOUBLE_MERGE_A_2(array1, indexarray1, nlen1, array2, indexarray2, nlen2, sortarray, sortindexarray, nlen) == PROC_FAILURE)
		{	
			free(array1);
			free(array2);
			free(indexarray1);
			free(indexarray2);

			return -1;
		}
	}

	/* release memory */
	free(array1);
	free(array2);
	free(indexarray1);
	free(indexarray2);

	/* return */
	return nlen;
}

/* ----------------------------------------------------------------------- */
/* SORT_MERGE_DOUBLE_MERGE_A_1 Function                                    */
/* This function is used for merging and sorting two double arrays.        */
/* The function will return PROC_SUCCESS;                                  */
/* if any error, return PROC_FAILURE.                                      */
/* ----------------------------------------------------------------------- */
int SORT_MERGE_DOUBLE_MERGE_A_1(double *array1, long nlen1, double *array2, long nlen2, double *sortarray, long nlen)
{
	long i,j,k;
	
	/* check parameter */
	i = 0;
	j = 0;
	k = 0;

	while((i < nlen1) && (j < nlen2))
	{
		if(array1[i] <= array2[j])
		{
			sortarray[k] = array1[i];
			i++;
		}
		else
		{
			sortarray[k] = array2[j];
			j++;
		}
		k++;
	}

	if(i == nlen1)
	{
		while(j < nlen2)
		{
			sortarray[k] = array2[j];
			j++;
			k++;
		}
	}
	else
	{
		while(i < nlen1)
		{
			sortarray[k] = array1[i];
			i++;
			k++;
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */
/* SORT_MERGE_DOUBLE_MERGE_A_2 Function                                    */
/* This function is used for merging and sorting two double arrays.        */
/* The function will return PROC_SUCCESS;                                  */
/* if any error, return PROC_FAILURE.                                      */
/* ----------------------------------------------------------------------- */
int SORT_MERGE_DOUBLE_MERGE_A_2(double *array1, long *indexarray1, long nlen1, double *array2, long *indexarray2, long nlen2, double *sortarray, long *sortindexarray, long nlen)
{
	long i,j,k;
	
	/* check parameter */
	i = 0;
	j = 0;
	k = 0;

	while((i < nlen1) && (j < nlen2))
	{
		if(array1[i] <= array2[j])
		{
			sortarray[k] = array1[i];
			sortindexarray[k] = indexarray1[i];
			i++;
		}
		else
		{
			sortarray[k] = array2[j];
			sortindexarray[k] = indexarray2[j];
			j++;
		}
		k++;
	}

	if(i == nlen1)
	{
		while(j < nlen2)
		{
			sortarray[k] = array2[j];
			sortindexarray[k] = indexarray2[j];
			j++;
			k++;
		}
	}
	else
	{
		while(i < nlen1)
		{
			sortarray[k] = array1[i];
			sortindexarray[k] = indexarray1[i];
			i++;
			k++;
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */
/* SORT_MERGE_INT_A Function                                               */
/* This function is used for sorting int array using ameliorated           */
/* two-way merge sort algorithm.                                           */
/* Time complexity: n*log2(n)                                              */
/* The function will return the length of sortarray;                       */
/* if any error, return -1.                                                */
/* ----------------------------------------------------------------------- */
long SORT_MERGE_INT_A(int *array, long *indexarray, int *sortarray, long *sortindexarray, long nstart, long nend)
{
	long nlen,nlen1,nlen2;
	int *array1;
	int *array2;
	long *indexarray1;
	long *indexarray2;

	/* check parameter */
	if(array == NULL)
	{
		printf("Wrong parameters for sorting!\n");
		return -1;
	}
	
	if(sortarray == NULL)
	{
		printf("Access Violation!\n");
		return -1;
	}
	
	if((nstart < 0) || (nend < nstart))
	{
		printf("Access Violation!\n");
		return -1;
	}
	
	if((indexarray != NULL) && (sortindexarray == NULL))
	{
		printf("Access Violation!\n");
		return -1;
	}

	nlen = nend-nstart+1;

	/* sort if only one element */
	if(nstart == nend)
	{
		sortarray[0] = array[nstart];
		if(indexarray != NULL)
			sortindexarray[0] = indexarray[nstart];

		return nlen;
	}

	/* if there are more than one elements */
	/* create additional arrays */
	nlen1 = (long)(nlen/2);
	nlen2 = nlen - nlen1;

	/* allocate memory */
	array1 = NULL;
	array2 = NULL;
	indexarray1 = NULL;
	indexarray2 = NULL;

	array1 = (int *)calloc(nlen1, sizeof(int));
	array2 = (int *)calloc(nlen2, sizeof(int));
	if((array1 == NULL) || (array2 == NULL))
	{
		free(array1);
		free(array2);

		printf("Can't allocate memory for sorting!\n");
		return -1;
	}
	if(indexarray != NULL)
	{
		indexarray1 = (long *)calloc(nlen1, sizeof(long));
		indexarray2 = (long *)calloc(nlen2, sizeof(long));
		if((indexarray1 == NULL) || (indexarray2 == NULL))
		{
			free(array1);
			free(array2);
			free(indexarray1);
			free(indexarray2);

			printf("Can't allocate memory for sorting!\n");
			return -1;
		}
	}

	/* sort */
	if(SORT_MERGE_INT_A(array, indexarray, array1, indexarray1, nstart, nstart+nlen1-1) != nlen1)
	{
		free(array1);
		free(array2);
		free(indexarray1);
		free(indexarray2);

		return -1;
	}
	
	if(SORT_MERGE_INT_A(array, indexarray, array2, indexarray2, nstart+nlen1, nend) != nlen2)
	{
		free(array1);
		free(array2);
		free(indexarray1);
		free(indexarray2);

		return -1;
	}

	if(indexarray == NULL)
	{
		if(SORT_MERGE_INT_MERGE_A_1(array1, nlen1, array2, nlen2, sortarray, nlen) == PROC_FAILURE)
		{	
			free(array1);
			free(array2);

			return -1;
		}
	}
	else
	{
		if(SORT_MERGE_INT_MERGE_A_2(array1, indexarray1, nlen1, array2, indexarray2, nlen2, sortarray, sortindexarray, nlen) == PROC_FAILURE)
		{	
			free(array1);
			free(array2);
			free(indexarray1);
			free(indexarray2);

			return -1;
		}
	}

	/* release memory */
	free(array1);
	free(array2);
	free(indexarray1);
	free(indexarray2);

	/* return */
	return nlen;
}

/* ----------------------------------------------------------------------- */
/* SORT_MERGE_INT_MERGE_A_1 Function                                       */
/* This function is used for merging and sorting two int arrays.           */
/* The function will return PROC_SUCCESS;                                  */
/* if any error, return PROC_FAILURE.                                      */
/* ----------------------------------------------------------------------- */
int SORT_MERGE_INT_MERGE_A_1(int *array1, long nlen1, int *array2, long nlen2, int *sortarray, long nlen)
{
	long i,j,k;
	
	/* check parameter */
	i = 0;
	j = 0;
	k = 0;

	while((i < nlen1) && (j < nlen2))
	{
		if(array1[i] <= array2[j])
		{
			sortarray[k] = array1[i];
			i++;
		}
		else
		{
			sortarray[k] = array2[j];
			j++;
		}
		k++;
	}

	if(i == nlen1)
	{
		while(j < nlen2)
		{
			sortarray[k] = array2[j];
			j++;
			k++;
		}
	}
	else
	{
		while(i < nlen1)
		{
			sortarray[k] = array1[i];
			i++;
			k++;
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */
/* SORT_MERGE_INT_MERGE_A_2 Function                                       */
/* This function is used for merging and sorting two int arrays.           */
/* The function will return PROC_SUCCESS;                                  */
/* if any error, return PROC_FAILURE.                                      */
/* ----------------------------------------------------------------------- */
int SORT_MERGE_INT_MERGE_A_2(int *array1, long *indexarray1, long nlen1, int *array2, long *indexarray2, long nlen2, int *sortarray, long *sortindexarray, long nlen)
{
	long i,j,k;
	
	/* check parameter */
	i = 0;
	j = 0;
	k = 0;

	while((i < nlen1) && (j < nlen2))
	{
		if(array1[i] <= array2[j])
		{
			sortarray[k] = array1[i];
			sortindexarray[k] = indexarray1[i];
			i++;
		}
		else
		{
			sortarray[k] = array2[j];
			sortindexarray[k] = indexarray2[j];
			j++;
		}
		k++;
	}

	if(i == nlen1)
	{
		while(j < nlen2)
		{
			sortarray[k] = array2[j];
			sortindexarray[k] = indexarray2[j];
			j++;
			k++;
		}
	}
	else
	{
		while(i < nlen1)
		{
			sortarray[k] = array1[i];
			sortindexarray[k] = indexarray1[i];
			i++;
			k++;
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */
/* DMSORTROWS Function                                                     */
/* This function sorts rows of a matrix pData.                             */
/* pPriority specifies which columns are used as sorting key, e.g. [2,5,0] */
/* pType specifies data types of each column: 0-double, 1-float, 2-int,    */
/*  3-2byte int, 4-1byte int, 5-uint, 6-2byte uint, 7-1byte uint.          */
/* pDataSort will store the sorted matrix, and pDataSid will store the     */
/* original index, i.e., pDataSort = pData[pDataSid.                       */
/* ----------------------------------------------------------------------- */
int DMSORTROWS(struct DOUBLEMATRIX *pData, struct INTMATRIX *pType, 
			   struct INTMATRIX *pPriority, struct DOUBLEMATRIX **pDataSort, 
			   struct LONGMATRIX **pDataSid)
{
	/* define */
	int ni,nj,nk,nId;
	int nTied;
	int nCol,nRow;
	struct DOUBLEMATRIX *pTempData;
	struct DOUBLEMATRIX *pTempSort;
	struct LONGMATRIX *pTempSid;
	double *vSrc,*vDst;
	long *vISrc,*vIDst;
	int nStart,nEnd;

	struct DOUBLEMATRIX *pX;
	struct LONGMATRIX *pXId;

	
	/* init check */
	if( (pDataSort == NULL) || (pDataSid == NULL) )
	{
		printf("Error: DMSORTROWS, invalid pointer to save sorting results!\n");
		exit(EXIT_FAILURE);
	}

	*pDataSort = NULL;
	*pDataSid = NULL;

	if(pData == NULL)
		return PROC_SUCCESS;

	if(pType == NULL)
	{
		printf("Error: DMSORTROWS, no data type specified!\n");
		return PROC_FAILURE;
	}

	if(pType->nWidth != pData->nWidth)
	{
		printf("Error: DMSORTROWS, data type specified incorrectly!\n");
		exit(EXIT_FAILURE);
	}

	if(pPriority == NULL)
	{
		printf("Warning: DMSORTROWS, sorting key not specified!\n");
		return PROC_SUCCESS;
	}

	if(pPriority->nWidth < 1)
	{
		printf("Warning: DMSORTROWS, sorting key not specified!\n");
		return PROC_SUCCESS;
	}

	/* create dest matrix */
	*pDataSort = CreateDoubleMatrix(pData->nHeight, pData->nWidth);
	if(*pDataSort == NULL)
	{
		printf("Error: DMSORTROWS, can't allocate memory for sorting!\n");
		exit(EXIT_FAILURE);
	}
	
	*pDataSid = CreateLongMatrix(pData->nHeight, 1);
	if(*pDataSid == NULL)
	{
		DestroyDoubleMatrix(*pDataSort);
		*pDataSort = NULL;
		printf("Error: DMSORTROWS, can't allocate memory for sorting!\n");
		exit(EXIT_FAILURE);
	}

	/* sort the first key */
	nCol = pPriority->pMatElement[0];
	if( (nCol<0) || (nCol>=pData->nWidth) )
	{
		printf("Error: DMSORTROWS, sorting key out of range!\n");
		exit(EXIT_FAILURE);
	}

	pTempData = NULL;
	pTempSort = NULL;
	pTempSid = NULL;
	
	pTempData = CreateDoubleMatrix(1, pData->nHeight);
	if(pTempData == NULL)
	{
		printf("Error: DMSORTROWS, can't allocate memory for storing sorting key values!\n");
		exit(EXIT_FAILURE);
	}

	for(nj=0; nj<pData->nHeight; nj++)
	{
		pTempData->pMatElement[nj] = DMGETAT(pData, nj, nCol);
	}

	DMSORTMERGEA_0(pTempData, &pTempSort, &pTempSid);

	vDst = (*pDataSort)->pMatElement;
	for(nj=0; nj<pTempSid->nWidth; nj++)
	{
		nRow = pTempSid->pMatElement[nj];
		vSrc = pData->pMatElement+nRow*pData->nWidth;
		memcpy(vDst, vSrc, sizeof(double)*pData->nWidth);
		vDst += pData->nWidth;
	}
	vIDst = (*pDataSid)->pMatElement;
	vISrc = pTempSid->pMatElement;
	memcpy(vIDst, vISrc, sizeof(long)*pTempSid->nWidth);

	DestroyDoubleMatrix(pTempData);
	DestroyDoubleMatrix(pTempSort);
	DestroyLongMatrix(pTempSid);


	/* sort the remaining keys */
	for(ni=1; ni<pPriority->nWidth; ni++)
	{
		nCol = pPriority->pMatElement[ni];
		if( (nCol<0) || (nCol>=pData->nWidth) )
		{
			printf("Error: DMSORTROWS, sorting key out of range!\n");
			exit(EXIT_FAILURE);
		}

		nStart = 0;
		nEnd = 0;
		for(nj=1; nj<(*pDataSort)->nHeight; nj++)
		{
			/* test if tied */
			nTied = 1;
			for(nk=0; nk<ni; nk++)
			{
				nId = pPriority->pMatElement[nk];
				switch(pType->pMatElement[nId])
				{
					case 0: 
						if( DMGETAT(*pDataSort, nj, nId) > DMGETAT(*pDataSort, nj-1, nId) )
							nTied = 0;
						break;
					case 1:
						if( fabs(DMGETAT(*pDataSort, nj, nId)-DMGETAT(*pDataSort, nj-1, nId)) > 1.175494351e-38 )
							nTied = 0;
						break;
					case 2:
						if( (int)(DMGETAT(*pDataSort, nj, nId)) > (int)(DMGETAT(*pDataSort, nj-1, nId)) )
							nTied = 0;
						break;
					case 3:
						if( (short)(DMGETAT(*pDataSort, nj, nId)) > (short)(DMGETAT(*pDataSort, nj-1, nId)) )
							nTied = 0;
						break;
					case 4:
						if( (char)(DMGETAT(*pDataSort, nj, nId)) > (char)(DMGETAT(*pDataSort, nj-1, nId)) )
							nTied = 0;
						break;
					case 5:
						if( (unsigned int)(DMGETAT(*pDataSort, nj, nId)) > (unsigned int)(DMGETAT(*pDataSort, nj-1, nId)) )
							nTied = 0;
						break;
					case 6:
						if( (unsigned short)(DMGETAT(*pDataSort, nj, nId)) > (unsigned short)(DMGETAT(*pDataSort, nj-1, nId)) )
							nTied = 0;
						break;
					case 7:
						if( (unsigned char)(DMGETAT(*pDataSort, nj, nId)) > (unsigned char)(DMGETAT(*pDataSort, nj-1, nId)) )
							nTied = 0;
						break;
					default:
						if( DMGETAT(*pDataSort, nj, nId) > DMGETAT(*pDataSort, nj-1, nId) )
							nTied = 0;
				}

				if(nTied == 0)
					break;
			}

			/* if tied, update nEnd */
			if(nTied == 1)
			{
				nEnd = nj;
			}

			/* if not tied, if there is an unsorted segment, sort the segment,
			   then update nStart and nEnd */
			else
			{
				if(nStart < nEnd)
				{
					pTempData = NULL;
					pTempSort = NULL;
					pTempSid = NULL;
					
					pTempData = CreateDoubleMatrix(1, nEnd-nStart+1);
					if(pTempData == NULL)
					{
						printf("Error: DMSORTROWS, can't allocate memory for storing sorting key values!\n");
						exit(EXIT_FAILURE);
					}

					for(nk=0; nk<pTempData->nWidth; nk++)
					{
						pTempData->pMatElement[nk] = DMGETAT(*pDataSort, (nStart+nk), nCol);
					}

					DMSORTMERGEA_0(pTempData, &pTempSort, &pTempSid);

					pX = NULL;
					pX = CreateDoubleMatrix(pTempData->nWidth, (*pDataSort)->nWidth);
					if(pX == NULL)
					{
						printf("Error: DMSORTROWS, can't allocate memory for storing sorting key values!\n");
						exit(EXIT_FAILURE);
					}
					pXId = NULL;
					pXId = CreateLongMatrix(1, pTempData->nWidth);
					if(pXId == NULL)
					{
						printf("Error: DMSORTROWS, can't allocate memory for storing sorting key values!\n");
						exit(EXIT_FAILURE);
					}

					vDst = pX->pMatElement;
					for(nk=0; nk<pTempSid->nWidth; nk++)
					{
						nRow = nStart+pTempSid->pMatElement[nk];
						vSrc = (*pDataSort)->pMatElement+nRow*(*pDataSort)->nWidth;
						memcpy(vDst, vSrc, sizeof(double)*(*pDataSort)->nWidth);
						vDst += pX->nWidth;
						pXId->pMatElement[nk] = (*pDataSid)->pMatElement[nRow];
					}

					vDst = (*pDataSort)->pMatElement+nStart*(*pDataSort)->nWidth;
					vSrc = pX->pMatElement;
					memcpy(vDst, vSrc, sizeof(double)*pX->nHeight*pX->nWidth);

					vIDst = (*pDataSid)->pMatElement+nStart;
					vISrc = pXId->pMatElement;
					memcpy(vIDst, vISrc, sizeof(long)*pXId->nWidth);

					DestroyDoubleMatrix(pX);
					DestroyLongMatrix(pXId);

					DestroyDoubleMatrix(pTempData);
					DestroyDoubleMatrix(pTempSort);
					DestroyLongMatrix(pTempSid);
				}

				nStart = nj;
				nEnd = nj;
			}
		}

		if(nStart < nEnd)
		{
			pTempData = NULL;
			pTempSort = NULL;
			pTempSid = NULL;
			
			pTempData = CreateDoubleMatrix(1, nEnd-nStart+1);
			if(pTempData == NULL)
			{
				printf("Error: DMSORTROWS, can't allocate memory for storing sorting key values!\n");
				exit(EXIT_FAILURE);
			}

			for(nk=0; nk<pTempData->nWidth; nk++)
			{
				pTempData->pMatElement[nk] = DMGETAT(*pDataSort, (nStart+nk), nCol);
			}

			DMSORTMERGEA_0(pTempData, &pTempSort, &pTempSid);

			pX = NULL;
			pX = CreateDoubleMatrix(pTempData->nWidth, (*pDataSort)->nWidth);
			if(pX == NULL)
			{
				printf("Error: DMSORTROWS, can't allocate memory for storing sorting key values!\n");
				exit(EXIT_FAILURE);
			}
			pXId = NULL;
			pXId = CreateLongMatrix(1, pTempData->nWidth);
			if(pXId == NULL)
			{
				printf("Error: DMSORTROWS, can't allocate memory for storing sorting key values!\n");
				exit(EXIT_FAILURE);
			}

			vDst = pX->pMatElement;
			for(nk=0; nk<pTempSid->nWidth; nk++)
			{
				nRow = nStart+pTempSid->pMatElement[nk];
				vSrc = (*pDataSort)->pMatElement+nRow*(*pDataSort)->nWidth;
				memcpy(vDst, vSrc, sizeof(double)*(*pDataSort)->nWidth);
				vDst += pX->nWidth;
				pXId->pMatElement[nk] = (*pDataSid)->pMatElement[nRow];
			}

			vDst = (*pDataSort)->pMatElement+nStart*(*pDataSort)->nWidth;
			vSrc = pX->pMatElement;
			memcpy(vDst, vSrc, sizeof(double)*pX->nHeight*pX->nWidth);

			vIDst = (*pDataSid)->pMatElement+nStart;
			vISrc = pXId->pMatElement;
			memcpy(vIDst, vISrc, sizeof(long)*pXId->nWidth);

			DestroyDoubleMatrix(pX);
			DestroyLongMatrix(pXId);

			DestroyDoubleMatrix(pTempData);
			DestroyDoubleMatrix(pTempSort);
			DestroyLongMatrix(pTempSid);
		}
	}

		
	
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */
/* DVNORM2 Function                                                        */
/* This function calculate the 2-norm of a double vector.                  */
/* Make sure you have allocated the required memory to the vetor pointer   */
/* before you call this function !                                         */
/* The function will return the 2-norm value;                              */
/* ----------------------------------------------------------------------- */
double DVNORM2(double *Vectorarray, long nLen)
{
	/* define */
	double dResult;
	double *pElement;
	long ni;

	/* check parameter */
	if(Vectorarray == NULL)
	{
		printf("Vector does not exist!\n");
		return 0.0;
	}

	/* init */
	dResult = 0.0;

	/* calculate */
	pElement = Vectorarray;
	for(ni=0; ni<nLen; ni++)
	{
		dResult += (*pElement)*(*pElement);
		pElement++;
	}

	dResult = sqrt(dResult);

	/* return */
	return dResult;
}

/* ----------------------------------------------------------------------- */
/* DVNORMPRODUCT Function                                                  */
/* This function calculate the normalized product of two double vectors.   */
/* Make sure you have allocated the required and identical memories to     */
/* vetors before you call this function !                                  */
/* The function will return the product value;                             */
/* ----------------------------------------------------------------------- */
double DVNORMPRODUCT(double *Vectorarray1, double *Vectorarray2, long nLen)
{
	/* define */
	double dResult;
	double *pElement1,*pElement2;
	double dNorm1,dNorm2;
	long ni;

	/* check parameter */
	if((Vectorarray1 == NULL) || (Vectorarray2 == NULL))
	{
		printf("Vector does not exist!\n");
		return 0.0;
	}

	/* init */
	dResult = 0.0;
	dNorm1 = 0.0;
	dNorm2 = 0.0;

	/* calculate */
	pElement1 = Vectorarray1;
	pElement2 = Vectorarray2;
	for(ni=0; ni<nLen; ni++)
	{
		dResult += (*pElement1)*(*pElement2);
		dNorm1 += (*pElement1)*(*pElement1);
		dNorm2 += (*pElement2)*(*pElement2);
		pElement1++;
		pElement2++;
	}

	dNorm1 = sqrt(dNorm1);
	dNorm2 = sqrt(dNorm2);

	if((dNorm1 < ZERO_BOUND) || (dNorm2 < ZERO_BOUND))
	{
		return 0.0;
	}

	dResult /= (dNorm1*dNorm2);
	if(dResult > 1.0)
		dResult = 1.0;
	if(dResult < -1.0)
		dResult = -1.0;

	/* return */
	return dResult;
}

/* ----------------------------------------------------------------------- */
/* BVNORM2 Function                                                        */
/* This function calculate the 2-norm of a byte vector.                    */
/* Make sure you have allocated the required memory to the vetor pointer   */
/* before you call this function !                                         */
/* The function will return the 2-norm value;                              */
/* ----------------------------------------------------------------------- */
double BVNORM2(unsigned char *Vectorarray, long nLen)
{
	/* define */
	double dResult;
	unsigned char *pElement;
	long ni;

	/* check parameter */
	if(Vectorarray == NULL)
	{
		printf("Vector does not exist!\n");
		return 0.0;
	}

	/* init */
	dResult = 0.0;

	/* calculate */
	pElement = Vectorarray;
	for(ni=0; ni<nLen; ni++)
	{
		dResult += (*pElement)*(*pElement);
		pElement++;
	}

	dResult = sqrt(dResult);

	/* return */
	return dResult;
}

/* ----------------------------------------------------------------------- */
/* BVNORMPRODUCT Function                                                  */
/* This function calculate the normalized product of two byte vectors.     */
/* Make sure you have allocated the required and identical memories to     */
/* vetors before you call this function !                                  */
/* The function will return the product value;                             */
/* ----------------------------------------------------------------------- */
double BVNORMPRODUCT(unsigned char *Vectorarray1, unsigned char *Vectorarray2, long nLen)
{
	/* define */
	double dResult;
	unsigned char *pElement1,*pElement2;
	double dNorm1,dNorm2;
	long ni;

	/* check parameter */
	if((Vectorarray1 == NULL) || (Vectorarray2 == NULL))
	{
		printf("Vector does not exist!\n");
		return 0.0;
	}

	/* init */
	dResult = 0.0;
	dNorm1 = 0.0;
	dNorm2 = 0.0;

	/* calculate */
	pElement1 = Vectorarray1;
	pElement2 = Vectorarray2;
	for(ni=0; ni<nLen; ni++)
	{
		dResult += (*pElement1)*(*pElement2);
		dNorm1 += (*pElement1)*(*pElement1);
		dNorm2 += (*pElement2)*(*pElement2);
		pElement1++;
		pElement2++;
	}

	dNorm1 = sqrt(dNorm1);
	dNorm2 = sqrt(dNorm2);

	if((dNorm1 < ZERO_BOUND) || (dNorm2 < ZERO_BOUND))
	{
		return 0.0;
	}

	dResult /= (dNorm1*dNorm2);

	/* return */
	return dResult;
}
/* ----------------------------------------------------------------------- */
/* pythag Function                                                         */
/* Compute (a^2+b^2)^1/2 without destructive underflow or overflow         */
/* ----------------------------------------------------------------------- */
double pythag(double a, double b)
{
	double absa,absb,temp;
	absa = fabs(a);
	absb = fabs(b);
	if(absa > absb)
	{
		temp = absb/absa;
		return absa*sqrt(1.0+temp*temp);
	}
	else
	{
		temp = absa/absb;
		return (absb == 0.0 ? 0.0: absb*sqrt(1.0+temp*temp));
	}
}



/* ----------------------------------------------------------------------- */
/* tred2 Function                                                          */
/* Householder reduction of a real, symmetric matrix a[1:n][1:n].          */
/* On output, a is replaced by the orthogonal matrix Q effecting the       */
/* transformation. d[1:n] returns the diagonal elements of the             */
/* resulting tridiagonal matrix, and e[1:n] the off-diagonal elements,     */
/* with e[1] = 0.                                                          */
/* ----------------------------------------------------------------------- */
void tred2(double **a, int n, double d[], double e[])
{
	int l,k,j,i;
	double scale,hh,h,g,f;

	for(i=(n-1); i>=1; i--)
	{
		l = i-1;
		scale = 0.0;
		h = scale;
		if( l > 0 )
		{
			for(k=0; k<=l; k++)
				scale += fabs(a[i][k]);

			if((float)scale == (float)0.0) /* skip transformation */
				e[i] = a[i][l];
			else
			{
				for(k=0; k<=l; k++)
				{
					a[i][k] /= scale; /* use scaled a's for transformation */
					h += a[i][k]*a[i][k]; /* form sigma in h */
				}

				f = a[i][l];
				g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i] = scale*g;
				h -= f*g;
				a[i][l] = f-g;
				f = 0.0;
				for(j=0; j<=l; j++)
				{
					a[j][i] = a[i][j]/h;
					g = 0.0;
					for(k=0; k<=j; k++)
						g += a[j][k]*a[i][k];
					for(k=j+1; k<=l; k++)
						g += a[k][j]*a[i][k];
					e[j] = g/h; /* form element of p in temporarily unused element of e */
					f += e[j]*a[i][j];
				}
				hh = f/(h+h); /* form K */
				for(j=0; j<=l; j++)
				{
					f = a[i][j];
					g = e[j]-hh*f;
					e[j] = g;
					for(k=0; k<=j; k++)
						a[j][k] -= (f*e[k]+g*a[i][k]);
				}
			}
		}
		else
		{
			e[i] = a[i][l];
		}
		d[i] = h;
	}

	d[0] = 0.0;
	e[0] = 0.0;
	for(i=0; i<n; i++)
	{
		l = i-1;
		if(d[i]) 
		{
			for(j=0; j<=l; j++)
			{
				g = 0.0;
				for(k=0; k<=l; k++)
					g += a[i][k]*a[k][j];
				for(k=0; k<=l; k++)
					a[k][j] -= g*a[k][i];
			}
		}
		d[i] = a[i][i];
		a[i][i] = 1.0;
		for(j=0; j<=l; j++)
		{
			a[j][i] = 0.0; /* reset row and column of a to identity matrix for next iteration*/
			a[i][j] = 0.0;
		}
	}
}

/* ----------------------------------------------------------------------- */
/* tqli Function                                                           */
/* QL algorithm with implicit shifts, to determine the eigenvalues and     */
/* eigenvectors of a real, symmetric tridiagonal matrix previously         */
/* reduced by tred2.                                                       */
/* ----------------------------------------------------------------------- */
void tqli(double d[], double e[], int n, double **z)
{
	int m,l,iter,i,k;
	double s,r,p,g,f,dd,c,b;

	for(i=1; i<n; i++)
		e[i-1] = e[i];
	e[n-1] = 0.0;

	for(l=0; l<n; l++)
	{
		iter = 0;
		do {
			for(m=l; m<=n-2; m++)
			{
				dd = fabs(d[m])+fabs(d[m+1]);
				if( (float)(fabs(e[m])+dd) == (float)dd )
					break;
			}
			if(m != l)
			{
				if(iter++ == 30)
				{
					printf("Error: too many iterations in tqli");
					exit(EXIT_FAILURE);
				}
				g = (d[l+1]-d[l])/(2*e[l]);
				r = pythag(g,1.0);
				g = d[m]-d[l]+e[l]/(g+SIGN(r,g));
				s = 1.0;
				c = 1.0;
				p = 0.0;

				for(i=m-1; i>=l; i--)
				{
					f = s*e[i];
					b = c*e[i];
					r = pythag(f,g);
					e[i+1] = r;
					if(r == 0.0)
					{
						d[i+1] -= p;
						e[m] = 0.0;
						break;
					}

					s = f/r;
					c = g/r;
					g = d[i+1]-p;
					r = (d[i]-g)*s+2.0*c*b;
					p = s*r;
					d[i+1] = g+p;
					g = c*r-b;

					for(k=0; k<n; k++)
					{
						f = z[k][i+1];
						z[k][i+1] = s*z[k][i]+c*f;
						z[k][i] = c*z[k][i]-s*f;
					}
				}

				if( ((float)r == 0.0) && (i>=l) )
					continue;

				d[l] -= p;
				e[l] = g;
				e[m] = 0.0;
	 		} 
		} while(m != l);
	}
}


/* ----------------------------------------------------------------------- */
/* eigsrt Function                                                         */
/* Given the eigenvalues d[1:n] and eigenvectors v[1:n][1:n], this         */
/* function sorts the eigenvalues into descending order, and rearranges    */
/* the columns of v correspondingly.                                       */
/* ----------------------------------------------------------------------- */
void eigsrt(double d[], double **v, int n)
{
	int k,j,i;
	double p;

	for(i=0; i<(n-1); i++)
	{
		k = i;
		p = d[k];

		for(j=i+1; j<n; j++)
		{
			if(d[j] >= p)
			{
				k = j;
				p = d[k];
			}
		}

		if(k != i)
		{
			d[k] = d[i];
			d[i] = p;
			for(j=0; j<n; j++)
			{
				p = v[j][i];
				v[j][i] = v[j][k];
				v[j][k] = p;
			}
		}
	}
}

/* ----------------------------------------------------------------------- */
/* jacobi Function                                                         */
/* Computes eigenvalues and eigenvectors of a real symmetrix matrix        */
/* a[1:n][1:n]. On output, elements of a above the diagonal are destroyed. */
/* d[1:n] returns the eigenvalues; v[1:n][1:n] is a matrix whose columns   */
/* contain the normalized eigenvectors. nrot returns the number of jacobi  */
/* rotations that were required.                                           */
/* ----------------------------------------------------------------------- */
void jacobi(double **a, int n, double d[], double **v, int *nrot)
{
	int j,iq,ip,i;
	double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

	b = NULL;
	z = NULL;
	b = (double *)calloc(n, sizeof(double));
	z = (double *)calloc(n, sizeof(double));
	if((b == NULL) || (z == NULL))
	{
		printf("Error: jacobi, cannot create memory for b and z!\n");
		exit(EXIT_FAILURE);
	}

	/* initialize to the identity matrix */
	for(ip=0; ip<n; ip++)
	{
		for(iq=0; iq<n; iq++)
			v[ip][iq] = 0.0;

		v[ip][ip] = 1.0;
	}

	/* initialize b and d to the diagonal of a */
	for(ip=0; ip<n; ip++)
	{
		b[ip] = a[ip][ip];
		d[ip] = a[ip][ip];
		z[ip] = 0.0;
	}

	*nrot = 0;
	for(i=1; i<=50; i++)
	{
		/* sum of off-diagonal elements */
		sm = 0.0;
		for(ip=0; ip<(n-1); ip++)
		{
			for(iq=ip+1; iq<n; iq++)
				sm += fabs(a[ip][iq]);
		}

		/* normal return */
		if( (float)sm == (float)0.0 )
		{
			free(b);
			free(z);
			return;
		}

		/* on the first three sweeps */
		if(i<4)
			tresh = 0.2*sm/(n*n);
		else
			tresh = 0.0;

		for(ip=0; ip<=(n-2); ip++)
		{
			for(iq=ip+1; iq<n; iq++)
			{
				g = 100.0*fabs(a[ip][iq]);
				if( (i>4) && ((float)(fabs(d[ip]+g)) == (float)(fabs(d[ip]))) && ((float)(fabs(d[iq]+g)) == (float)(fabs(d[iq]))) )
					a[ip][iq] = 0.0;
				else if(fabs(a[ip][iq]) > tresh)
				{
					h = d[iq]-d[ip];
					if( (float)(fabs(h)+g) == (float)(fabs(h)) )
						t = a[ip][iq]/h;
					else
					{
						theta = 0.5*h/a[ip][iq];
						t = 1.0/(fabs(theta)+sqrt(1.0+theta*theta));
						if(theta < 0.0)
							t = -t;
					}

					c = 1.0/sqrt(1+t*t);
					s = t*c;
					tau = s/(1.0+c);
					h = t*a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[ip][iq] = 0.0;
					for(j=0; j<=(ip-1); j++)
					{
						ROTATE(a,j,ip,j,iq);
					}
					for(j=ip+1; j<=iq-1; j++)
					{
						ROTATE(a,ip,j,j,iq);
					}
					for(j=iq+1; j<n; j++)
					{
						ROTATE(a,ip,j,iq,j);
					}
					for(j=0; j<n; j++)
					{
						ROTATE(v,j,ip,j,iq);
					}
					++(*nrot);
				}
			}
		}

		for(ip=0; ip<n; ip++)
		{
			b[ip] += z[ip];
			d[ip] = b[ip];
			z[ip] = 0.0;
		}
	}

	printf("Error: jacobi, too many iterations in routine jacobi!\n");

}

/* ----------------------------------------------------------------------- */
/* DMSYMEIGEN Function                                                     */
/* This function computes eigen values and eigen vectors of a symmetrix    */
/* matrix.                                                                 */
/* pX is the input matrix, pD is a vector that contains eigen values       */
/* pV is a matrix whose columns are eigen vectors with norm = 1.           */
/* ----------------------------------------------------------------------- */
int DMSYMEIGEN(struct DOUBLEMATRIX *pX,struct DOUBLEMATRIX **pD,struct DOUBLEMATRIX **pV)
{
	/* define */
	struct DOUBLEMATRIX *pY;
	double **a;
	double **v;
	int ni;
	int nrot = 0;

	/* init */
	pY = NULL;
	*pD = NULL;
	*pV = NULL;

	/* if source matrixes don't exist , return PROC_FAILURE */
	if(pX==NULL)
	{
		return PROC_FAILURE;
	}

	pY = DMCLONE(pX);
	*pV = CreateDoubleMatrix(pX->nHeight, pX->nWidth);
	*pD = CreateDoubleMatrix(1, pX->nHeight);

	if((pY == NULL) || (*pD == NULL) || (*pV == NULL) )
	{
		printf("Error: DMSYMEIGEN, cannot create memory for finding eigenvalues and eigenvectors!\n");
		exit(EXIT_FAILURE);
	}

	a = NULL;
	a = (double **)calloc(pX->nHeight, sizeof(double *));
	v = NULL;
	v = (double **)calloc(pX->nHeight, sizeof(double *));
	if( (a == NULL) || (v == NULL) )
	{
		printf("Error: DMSYMEIGEN, cannot create memory for organizing input data!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<pX->nHeight; ni++)
	{
		a[ni] = pY->pMatElement+ni*(pX->nWidth);
		v[ni] = (*pV)->pMatElement+ni*(pX->nWidth);
	}

	jacobi(a, pX->nHeight, (*pD)->pMatElement, v, &nrot);
	eigsrt((*pD)->pMatElement, v, pX->nHeight);

	free(a);
	free(v);

	DestroyDoubleMatrix(pY);
	pY = NULL;

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */
/* DMSYMINVEIGEN Function                                                  */
/* This function computes inverse of a symmetric matrix based on its eigen */
/* decomposition. pD and pV are output of DMSYMEIGEN                       */
/* ----------------------------------------------------------------------- */
int DMSYMINVEIGEN(struct DOUBLEMATRIX *pD,struct DOUBLEMATRIX *pV, struct DOUBLEMATRIX **pI)
{
	/* define */
	struct DOUBLEMATRIX *pVC,*pVT;
	double *pElement;
	int ni,nj;

	/* init */
	*pI = NULL;
	pVC = NULL;
	pVT = NULL;

	if(pD == NULL || pV == NULL)
	{
		return PROC_FAILURE;
	}

	for(ni=0; ni<pV->nHeight; ni++)
	{
		if( (float)(pD->pMatElement[ni]) == (float)0.0 )
		{
			printf("Error: DMSYMINVEIGEN, singular matrix, cannot be inverted\n");
			exit(EXIT_FAILURE);
		}

		if( fabs(pD->pMatElement[ni]) < 1e-10 )
		{
			printf("Warning: DMSYMINVEIGEN, fabs(eigenvalue_%d) < 1e-10, close to a singular matrix\n", (ni+1));
		}
	}

	pVC = DMCLONE(pV);
	pVT = DM_T(pV);

	pElement = pVC->pMatElement;
	for(ni=0; ni<pVC->nHeight; ni++)
	{
		for(nj=0; nj<pVC->nWidth; nj++)
		{
			(*pElement) = (*pElement)/pD->pMatElement[nj];
			pElement++;
		}
	}

	*pI = DMMUL(pVC, pVT);

	DestroyDoubleMatrix(pVC);
	DestroyDoubleMatrix(pVT);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */
/* DMSYMEIGENQR Function                                                   */
/* This function computes eigen values and eigen vectors of a symmetrix    */
/* matrix using QR decomposition                                           */
/* pX is the input matrix, pD is a vector that contains eigen values       */
/* pV is a matrix whose columns are eigen vectors with norm = 1.           */
/* ----------------------------------------------------------------------- */
int DMSYMEIGENQR(struct DOUBLEMATRIX *pX,struct DOUBLEMATRIX **pD,struct DOUBLEMATRIX **pV)
{
	/* define */
	double **a;
	double *e;
	int ni;
	int nrot = 0;

	/* init */
	*pD = NULL;
	*pV = NULL;

	/* if source matrixes don't exist , return PROC_FAILURE */
	if(pX==NULL)
	{
		return PROC_FAILURE;
	}

	*pV = DMCLONE(pX);
	*pD = CreateDoubleMatrix(1, pX->nHeight);

	if((*pD == NULL) || (*pV == NULL))
	{
		printf("Error: DMSYMEIGEN, cannot create memory for finding eigenvalues and eigenvectors!\n");
		exit(EXIT_FAILURE);
	}

	a = NULL;
	a = (double **)calloc(pX->nHeight, sizeof(double *));
	e = NULL;
	e = (double *)calloc(pX->nHeight, sizeof(double));
	if( (a == NULL) || (e == NULL))
	{
		printf("Error: DMSYMEIGEN, cannot create memory for organizing input data!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<pX->nHeight; ni++)
	{
		a[ni] = (*pV)->pMatElement+ni*(pX->nWidth);
	}

	tred2(a, pX->nHeight, (*pD)->pMatElement, e);
	tqli((*pD)->pMatElement, e, pX->nHeight, a);
	eigsrt((*pD)->pMatElement, a, pX->nHeight);

	free(a);
	free(e);

	/* return */
	return PROC_SUCCESS;
}
