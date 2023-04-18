/* ----------------------------------------------------------------------- */
/*  MatrixLib.h : interface of the matrix process library                  */
/*  Author : Ji HongKai ; Time: 1999.11                                    */
/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------- */
/* Declaration                                                             */
/* ----------------------------------------------------------------------- */
#define ZERO_BOUND 1e-250
#define DM_ACCESS_VIOLATION 1e250
#define IM_ACCESS_VIOLATION INT_MAX
#define LM_ACCESS_VIOLATION LONG_MAX
#define BM_ACCESS_VIOLATION 255
#define DM_LOAD_MAX_WIDTH 4096
#define LM_LOAD_MAX_WIDTH 4096
#define IM_LOAD_MAX_WIDTH 4096
#define BM_LOAD_MAX_WIDTH 4096

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);


/* ----------------------------------------------------------------------- */
/*  Definition of BYTEMATRIX struct.                                       */
/*  BYTEMATRIX is used for manipulating matrix , the elements of which are */
/*  of byte(unsigned char) type.                                           */ 
/* ----------------------------------------------------------------------- */
struct BYTEMATRIX
{
	/* nWidth is the width of the matrix associated with pMatElement. */
	long nWidth;
	/* nHeight is the height of the matrix associated with pMatElement. */
	long nHeight;
	/* pMatElement is the pointer to the element array. */
	unsigned char * pMatElement;
};

/* ----------------------------------------------------------------------- */
/*  Definition of INTMATRIX struct.                                        */
/*  INTMATRIX is used for manipulating matrix , the elements of which are  */
/*  of integer(int) type.                                                  */ 
/* ----------------------------------------------------------------------- */
struct INTMATRIX
{
	/* nWidth is the width of the matrix associated with pMatElement. */
	long nWidth;
	/* nHeight is the height of the matrix associated with pMatElement. */
	long nHeight;
	/* pMatElement is the pointer to the element array. */
	int * pMatElement;
};

/* ----------------------------------------------------------------------- */
/*  Definition of LONGMATRIX struct.                                       */
/*  LONGMATRIX is used for manipulating matrix , the elements of which are */
/*  of long type.                                                          */ 
/* ----------------------------------------------------------------------- */
struct LONGMATRIX
{
	/* nWidth is the width of the matrix associated with pMatElement. */
	long nWidth;
	/* nHeight is the height of the matrix associated with pMatElement. */
	long nHeight;
	/* pMatElement is the pointer to the element array. */
	long * pMatElement;
};

/* ----------------------------------------------------------------------- */
/*  Definition of DOUBLEMATRIX struct.                                     */
/*  DOUBLEMATRIX is used for manipulating matrix , the elements of which   */
/*  are of double type.											           */ 
/* ----------------------------------------------------------------------- */
struct  DOUBLEMATRIX
{
	/* nWidth is the width of the matrix associated with pMatElement. */
	long nWidth;
	/* nHeight is the height of the matrix associated with pMatElement. */
	long nHeight;
	/* pMatElement is the pointer to the element array. */
	double * pMatElement;
};

/* ----------------------------------------------------------------------- */
/* CreateByteMatrix Function                                               */
/* This function is used for creating a new BYTEMATRIX .                   */
/* The size of the matrix is given by two parameters : nHeight , nWidth    */
/* In other words , a nHeight*nWidth matrix will be created.               */
/* If it fails to create new matrix , the return value will be NULL .      */
/* Else the function will return the pointer to the new BYTEMATRIX .       */
/* ----------------------------------------------------------------------- */
struct BYTEMATRIX* CreateByteMatrix(long nHeight, long nWidth);

/* ----------------------------------------------------------------------- */
/* DestroyByteMatrix Function                                              */
/* This function is used for destroying a BYTEMATRIX .                     */
/* First, the memory block of elements will be freed .                     */
/* Then the memory block of the BYTEMATRIX struct will be freed .          */
/* ----------------------------------------------------------------------- */
void DestroyByteMatrix(struct BYTEMATRIX* pDelMatrix); 

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
unsigned char BMGETAT(struct BYTEMATRIX* pSourceMatrix,long nRow,long nCol);

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
int BMSETAT(struct BYTEMATRIX* pSourceMatrix,long nRow,long nCol,unsigned char bValue);

/* ----------------------------------------------------------------------- */
/* BMSAVE Function                                                         */
/* This function is used for saving a BYTEMATRIX.                          */
/* If successful, return PROC_SUCCESS;                                     */
/* else return PROC_FAILURE.                                               */
/* ----------------------------------------------------------------------- */
int BMSAVE(struct BYTEMATRIX* pByteMatrix, char strFilePath[]);

/* ----------------------------------------------------------------------- */
/* BMLOAD Function                                                         */
/* This function is used for loading a BYTEMATRIX.                         */
/* If successful, return the pointer of the matrix;                        */
/* else return NULL.                                                       */
/* ----------------------------------------------------------------------- */
struct BYTEMATRIX* BMLOAD(char strFilePath[]);

/* ----------------------------------------------------------------------- */
/* ByteLoadRowVector                                                       */
/* This function is used for loading a byte row vector from a string.      */
/* If successful, return the count of loaded elements;                     */
/* else return 0.                                                          */
/* ----------------------------------------------------------------------- */
long ByteLoadRowVector(char strInLine[], unsigned char bElement[]);

/* ----------------------------------------------------------------------- */
/* BMADDROW Function                                                       */
/* This function is used for adding a row to BYTEMATRIX.                   */
/* If successful, return PROC_SUCCESS;                                     */
/* else return PROC_FAILURE.                                               */
/* ----------------------------------------------------------------------- */
int BMADDROW(struct BYTEMATRIX **pByteMatrix, unsigned char bElement[], long nElementCount);

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
struct BYTEMATRIX* BM_T(struct BYTEMATRIX* pX);

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
int BMCOPY(struct BYTEMATRIX* pDestByteMatrix,struct BYTEMATRIX* pSourceByteMatrix);

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
struct BYTEMATRIX* BMCLONE(struct BYTEMATRIX* pSourceByteMatrix);

/* ----------------------------------------------------------------------- */
/* BMSUM Function                                                          */
/* This function computes the sum of all elements of a BM matrix.          */
/* ----------------------------------------------------------------------- */
double BMSUM(struct BYTEMATRIX* pByteMatrix);

/* ----------------------------------------------------------------------- */
/* BMPERMUTEROWS Function                                                  */
/* This function permutes the rows of the source matrix, and retrun a new  */
/* matrix which is the permutated one.                                     */
/* If the new matrix can't be created or other exceptions occur , the      */
/* function will return NULL ;                                             */
/* Else the function will permute the source matrix and return the pointer */
/* of the new matrix.                                                      */
/* ----------------------------------------------------------------------- */
struct BYTEMATRIX* BMPERMUTEROWS(struct BYTEMATRIX* pSourceByteMatrix);

/* ----------------------------------------------------------------------- */
/* BMPERMUTEELEMENTINONECOLUMN Function                                    */
/* This function permutes the elements of the specific column of the       */
/* source matrix, and retrun a new matrix which is the permutated one.     */
/* If the new matrix can't be created or other exceptions occur , the      */
/* function will return NULL ;                                             */
/* Else the function will permute the source matrix and return the pointer */
/* of the new matrix.                                                      */
/* ----------------------------------------------------------------------- */
struct BYTEMATRIX* BMPERMUTEELEMENTSINONECOLUMN(struct BYTEMATRIX* pSourceByteMatrix, long nColumnIndex);

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
int BMBOOTSTRAP(struct BYTEMATRIX* pDestByteMatrix,struct BYTEMATRIX* pSourceByteMatrix);

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
int BMBOOTSTRAPC(struct BYTEMATRIX* pDestByteMatrix,struct BYTEMATRIX* pSourceByteMatrix);

/* ----------------------------------------------------------------------- */
/* CreateIntMatrix Function                                                */
/* This function is used for creating a new INTMATRIX .                    */
/* The size of the matrix is given by two parameters : nHeight , nWidth    */
/* In other words , a nHeight*nWidth matrix will be created.               */
/* If it fails to create new matrix , the return value will be NULL .      */
/* Else the function will return the pointer to the new INTMATRIX .        */
/* ----------------------------------------------------------------------- */
struct INTMATRIX* CreateIntMatrix(long nHeight, long nWidth);

/* ----------------------------------------------------------------------- */
/* DestroyIntMatrix Function                                               */
/* This function is used for destroying a INTMATRIX .                      */
/* First, the memory block of elements will be freed .                     */
/* Then the memory block of the INTMATRIX struct will be freed .           */
/* ----------------------------------------------------------------------- */
void DestroyIntMatrix(struct INTMATRIX* pDelMatrix); 

/* ----------------------------------------------------------------------- */
/* IMGETAT Function                                                        */
/* This function is used for getting element indicated by row index and    */
/* column index .                                                          */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateINTMatrix function to create it           */
/* before you call this function !                                         */
/* If the row(col) indexes is larger than or equal to the row(col) number  */
/* of the matrix or smaller than 0 or some other exceptions occur, the     */
/* function will return INT_MAX;                                           */  
/* Else the function will return matrix[row][col] .                        */
/* ----------------------------------------------------------------------- */
int IMGETAT(struct INTMATRIX* pSourceMatrix,long nRow,long nCol); 

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
int IMSETAT(struct INTMATRIX* pSourceMatrix,long nRow,long nCol,int nValue);

/* ----------------------------------------------------------------------- */
/* IMSAVE Function                                                         */
/* This function is used for saving a INTMATRIX.                           */
/* If successful, return PROC_SUCCESS;                                     */
/* else return PROC_FAILURE.                                               */
/* ----------------------------------------------------------------------- */
int IMSAVE(struct INTMATRIX* pIntMatrix, char strFilePath[]);

/* ----------------------------------------------------------------------- */
/* IMLOAD Function                                                         */
/* This function is used for loading a INTMATRIX.                          */
/* If successful, return the pointer of the matrix;                        */
/* else return NULL.                                                       */
/* ----------------------------------------------------------------------- */
struct INTMATRIX* IMLOAD(char strFilePath[]);

/* ----------------------------------------------------------------------- */
/* IntLoadRowVector                                                        */
/* This function is used for loading a int row vector from a string.       */
/* If successful, return the count of loaded elements;                     */
/* else return 0.                                                          */
/* ----------------------------------------------------------------------- */
long IntLoadRowVector(char strInLine[], int nElement[]);

/* ----------------------------------------------------------------------- */
/* IMADDROW Function                                                       */
/* This function is used for adding a row to INTMATRIX.                    */
/* If successful, return PROC_SUCCESS;                                     */
/* else return PROC_FAILURE.                                               */
/* ----------------------------------------------------------------------- */
int IMADDROW(struct INTMATRIX **pIntMatrix, int nElement[], long nElementCount);

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
struct INTMATRIX* IM_T(struct INTMATRIX* pX);

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
int IMCOPY(struct INTMATRIX* pDestIntMatrix,struct INTMATRIX* pSourceIntMatrix);

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
struct INTMATRIX* IMCLONE(struct INTMATRIX* pSourceIntMatrix);

/* ----------------------------------------------------------------------- */
/* IMSORTMERGEA_0 Function                                                 */
/* This function is used for sorting rows of a int matrix using            */
/* ameliorated two-way merge algorithm.                                    */
/* The function will return PROC_SUCCESS;                                  */
/* if any error, return PROC_FAILURE.                                      */
/* ----------------------------------------------------------------------- */
int IMSORTMERGEA_0(struct INTMATRIX* pSourceMatrix, struct INTMATRIX** pSortMatrix, struct LONGMATRIX** pSortIndexMatrix);

/* ----------------------------------------------------------------------- */
/* CreateLongMatrix Function                                               */
/* This function is used for creating a new LONGMATRIX .                   */
/* The size of the matrix is given by two parameters : nHeight , nWidth    */
/* In other words , a nHeight*nWidth matrix will be created.               */
/* If it fails to create new matrix , the return value will be NULL .      */
/* Else the function will return the pointer to the new LONGMATRIX .       */
/* ----------------------------------------------------------------------- */
struct LONGMATRIX* CreateLongMatrix(long nHeight, long nWidth);

/* ----------------------------------------------------------------------- */
/* DestroyLongMatrix Function                                              */
/* This function is used for destroying a LONGMATRIX .                     */
/* First, the memory block of elements will be freed .                     */
/* Then the memory block of the LONGMATRIX struct will be freed.           */
/* ----------------------------------------------------------------------- */
void DestroyLongMatrix(struct LONGMATRIX* pDelMatrix);

/* ----------------------------------------------------------------------- */
/* LMSAVE Function                                                         */
/* This function is used for saving a LONGMATRIX.                          */
/* If successful, return PROC_SUCCESS;                                     */
/* else return PROC_FAILURE.                                               */
/* ----------------------------------------------------------------------- */
int LMSAVE(struct LONGMATRIX* pLongMatrix, char strFilePath[]);

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
long LMGETAT(struct LONGMATRIX* pSourceMatrix,long nRow,long nCol); 

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
int LMSETAT(struct LONGMATRIX* pSourceMatrix,long nRow,long nCol,long nValue);

/* ----------------------------------------------------------------------- */
/* LMLOAD Function                                                         */
/* This function is used for loading a LONGMATRIX.                         */
/* If successful, return the pointer of the matrix;                        */
/* else return NULL.                                                       */
/* ----------------------------------------------------------------------- */
struct LONGMATRIX* LMLOAD(char strFilePath[]);

/* ----------------------------------------------------------------------- */
/* LongLoadRowVector                                                       */
/* This function is used for loading a long row vector from a string.      */
/* If successful, return the count of loaded elements;                     */
/* else return 0.                                                          */
/* ----------------------------------------------------------------------- */
long LongLoadRowVector(char strInLine[], long nElement[]);

/* ----------------------------------------------------------------------- */
/* LMADDROW Function                                                       */
/* This function is used for adding a row to LONGMATRIX.                   */
/* If successful, return PROC_SUCCESS;                                     */
/* else return PROC_FAILURE.                                               */
/* ----------------------------------------------------------------------- */
int LMADDROW(struct LONGMATRIX **pLongMatrix, long nElement[], long nElementCount);

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
struct LONGMATRIX* LM_T(struct LONGMATRIX* pX);

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
int LMCOPY(struct LONGMATRIX* pDestLongMatrix,struct LONGMATRIX* pSourceLongMatrix);

/* ----------------------------------------------------------------------- */
/* CreateDoubleMatrix Function                                             */
/* This function is used for creating a new DOUBLEMATRIX .                 */
/* The size of the matrix is given by two parameters : nHeight , nWidth    */
/* In other words , a nHeight*nWidth matrix will be created.               */
/* If it fails to create new matrix , the return value will be NULL .      */
/* Else the function will return the pointer to the new DOUBLEMATRIX .     */
/* ----------------------------------------------------------------------- */
struct DOUBLEMATRIX* CreateDoubleMatrix(long nHeight, long nWidth);

/* ----------------------------------------------------------------------- */
/* DestroyDoubleMatrix Function                                            */
/* This function is used for destroying a DOUBLEMATRIX .                   */
/* First, the memory block of elements will be freed .                     */
/* Then the memory block of the DOUBLEMATRIX struct will be freed .        */
/* ----------------------------------------------------------------------- */
void DestroyDoubleMatrix(struct DOUBLEMATRIX* pDelMatrix);

/* ----------------------------------------------------------------------- */
/* DMRANDU Function                                                        */
/* This function is used for creating a new rand DOUBLEMATRIX .            */
/* The size of the matrix is given by two parameters : nHeight , nWidth    */
/* In other words , a nHeight*nWidth matrix will be created, and the       */
/* elements of the matrix is uniformly distributed on (0,1)                */ 
/* If it fails to create new matrix , the return value will be NULL .      */
/* Else the function will return the pointer to the new DOUBLEMATRIX .     */
/* ----------------------------------------------------------------------- */
struct DOUBLEMATRIX* DMRANDU(long nHeight, long nWidth);

/* ----------------------------------------------------------------------- */
/* DMRANDNORM Function                                                     */
/* This function is used for creating a new normal rand DOUBLEMATRIX .     */
/* The size of the matrix is given by two parameters : nHeight , nWidth    */
/* In other words , a nHeight*nWidth matrix will be created, and the       */
/* elements of the matrix is uniformly distributed on (0,1)                */ 
/* If it fails to create new matrix , the return value will be NULL .      */
/* Else the function will return the pointer to the new DOUBLEMATRIX .     */
/* ----------------------------------------------------------------------- */
struct DOUBLEMATRIX* DMRANDNORM(long nHeight, long nWidth, double mu, double sigma);

/* ----------------------------------------------------------------------- */
/* DMRANDPRODDIRICHLET Function                                            */
/* This function is used for creating a new product dirichlet rand         */
/* DOUBLEMATRIX of parameter pParamMatrix.                                 */
/* Each row of the pParamMatrix is a parameter vector for one component of */
/* product dirichlet.                                                      */
/* If it fails to create new matrix , the return value will be NULL .      */
/* Else the function will return the pointer to the new DOUBLEMATRIX .     */
/* ----------------------------------------------------------------------- */
struct DOUBLEMATRIX* DMRANDPRODDIRICHLET(struct DOUBLEMATRIX* pParamMatrix);

/* ----------------------------------------------------------------------- */
/* DMPRODDIRICHLETRND Function                                             */
/* This function is used for creating product dirichlet random numbers of  */
/* of parameter pParamMatrix, and assigning the value to pOutMatrix.       */
/* Each row of the pParamMatrix is a parameter vector for one component of */
/* product dirichlet.                                                      */
/* If it fails to create new matrix, the return value will be PROC_FAILURE.*/
/* Else the function will return PROC_SUCCESS.                             */
/* ----------------------------------------------------------------------- */
int DMPRODDIRICHLETRND(struct DOUBLEMATRIX* pOutMatrix, struct DOUBLEMATRIX* pParamMatrix);

/* ----------------------------------------------------------------------- */
/* DMSAVE Function                                                         */
/* This function is used for saving a DOUBLEMATRIX.                        */
/* If successful, return PROC_SUCCESS;                                     */
/* else return PROC_FAILURE.                                               */
/* ----------------------------------------------------------------------- */
int DMSAVE(struct DOUBLEMATRIX* pDoubleMatrix, char strFilePath[]);

/* ----------------------------------------------------------------------- */
/* DMLOAD Function                                                         */
/* This function is used for loading a DOUBLEMATRIX.                       */
/* If successful, return the pointer of the matrix;                        */
/* else return NULL.                                                       */
/* ----------------------------------------------------------------------- */
struct DOUBLEMATRIX* DMLOAD(char strFilePath[]);

/* ----------------------------------------------------------------------- */
/* DoubleLoadRowVector                                                     */
/* This function is used for loading a double row vector from a string.    */
/* If successful, return the count of loaded elements;                     */
/* else return 0.                                                          */
/* ----------------------------------------------------------------------- */
long DoubleLoadRowVector(char strInLine[], double dElement[]);

/* ----------------------------------------------------------------------- */
/* DMADDROW Function                                                       */
/* This function is used for adding a row to DOUBLEMATRIX.                 */
/* If successful, return PROC_SUCCESS;                                     */
/* else return PROC_FAILURE.                                               */
/* ----------------------------------------------------------------------- */
int DMADDROW(struct DOUBLEMATRIX **pDoubleMatrix, double dElement[], long nElementCount);

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
int DMCOPY(struct DOUBLEMATRIX* pDestDoubleMatrix,struct DOUBLEMATRIX* pSourceDoubleMatrix);

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
struct DOUBLEMATRIX* DMCLONE(struct DOUBLEMATRIX* pSourceDoubleMatrix);

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
struct DOUBLEMATRIX* DM_T(struct DOUBLEMATRIX* pX);

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
struct DOUBLEMATRIX* DMADD(struct DOUBLEMATRIX* pDoubleMatrix1,struct DOUBLEMATRIX* pDoubleMatrix2); 

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
struct DOUBLEMATRIX* DMSUB(struct DOUBLEMATRIX* pDoubleMatrix1,struct DOUBLEMATRIX* pDoubleMatrix2); 
    
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
int DMADDTF(struct DOUBLEMATRIX* pDoubleMatrix1,struct DOUBLEMATRIX* pDoubleMatrix2); 

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
int DMSUBTF(struct DOUBLEMATRIX* pDoubleMatrix1,struct DOUBLEMATRIX* pDoubleMatrix2); 

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
int DMSUBTL(struct DOUBLEMATRIX* pDoubleMatrix1,struct DOUBLEMATRIX* pDoubleMatrix2);

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
int DMPADDTS(struct DOUBLEMATRIX* pSourceMatrix, double dLamda);

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
int DMPSUBTS(double dLamda, struct DOUBLEMATRIX* pSourceMatrix);

/* ----------------------------------------------------------------------- */
/* DMPMUL Function                                                         */ 
/* This function multiplies a scalar quantity with a DOUBLEMATRIX and      */
/* records the result to a new DOUBLEMATRIX .                              */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateDoubleMatrix function to create it        */
/* before you call this function !                                         */
/* If the new matrix can't be created or other exceptions occur , the      */
/* function will return NULL ;                                             */
/* Else the function will let DestMatrix = dLamda *SourceMatrix and return */
/* the pointer to DestMatrix .                                             */
/* ----------------------------------------------------------------------- */
struct DOUBLEMATRIX* DMPMUL(double dLamda, struct DOUBLEMATRIX* pSourceMatrix);

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
int DMPMULTS(double dLamda, struct DOUBLEMATRIX* pSourceMatrix);

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
int DMPMULTF(struct DOUBLEMATRIX* pDoubleMatrix1,struct DOUBLEMATRIX* pDoubleMatrix2);

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
struct DOUBLEMATRIX* DMPDIV(struct DOUBLEMATRIX* pSourceMatrix, double dLamda);

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
int DMPDIVTS(struct DOUBLEMATRIX* pSourceMatrix, double dLamda);

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
int DMPRECTS(double dLamda, struct DOUBLEMATRIX* pSourceMatrix);

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
int DMPDIVTF(struct DOUBLEMATRIX* pDoubleMatrix1,struct DOUBLEMATRIX* pDoubleMatrix2);

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
int DMPDIVTL(struct DOUBLEMATRIX* pDoubleMatrix1,struct DOUBLEMATRIX* pDoubleMatrix2);

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
struct DOUBLEMATRIX* DMMUL(struct DOUBLEMATRIX* pSourceMatrix1,struct DOUBLEMATRIX* pSourceMatrix2);

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
int DMROWEXCH(struct DOUBLEMATRIX* pSourceMatrix,long nRow1,long nRow2);  

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
int DMROWPMUL(struct DOUBLEMATRIX* pSourceMatrix,long nRow,double dLamda);

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
int DMROWADDTF(struct DOUBLEMATRIX* pSourceMatrix,long nRow1,long nRow2);

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
int DMROWSUBTF(struct DOUBLEMATRIX* pSourceMatrix,long nRow1,long nRow2);

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
			   long nRow2);

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
int DMROWNORM(struct DOUBLEMATRIX* pSourceMatrix);

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
int DMLOGTS(struct DOUBLEMATRIX* pSourceMatrix);

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
int DMPLOGTS(struct DOUBLEMATRIX* pSourceMatrix, double dBase);

/* ----------------------------------------------------------------------- */
/* DMPPOWTS Function                                                       */ 
/* This function computes power of a DOUBLEMATRIX and records              */
/* the result to itself.                                                   */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateDoubleMatrix function to create it        */
/* before you call this function !                                         */
/* The function will compute matrix^dPow and return PROC_SUCCESS.          */
/* ----------------------------------------------------------------------- */
int DMPPOWTS(struct DOUBLEMATRIX* pSourceMatrix, double dPow);

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
int DMPPOWTF(struct DOUBLEMATRIX* pDoubleMatrix1,struct DOUBLEMATRIX* pDoubleMatrix2);

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
int DMPPOWTL(struct DOUBLEMATRIX* pDoubleMatrix1,struct DOUBLEMATRIX* pDoubleMatrix2);

/* ----------------------------------------------------------------------- */
/* DMPUPOWTS Function                                                      */ 
/* This function uses a DOUBLEMATRIX to compute power of a number and      */
/* records the result to the matrix itself.                                */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateDoubleMatrix function to create it        */
/* before you call this function !                                         */
/* The function will compute dBase^matrix and return PROC_SUCCESS.         */
/* ----------------------------------------------------------------------- */
int DMPUPOWTS(double dBase, struct DOUBLEMATRIX* pSourceMatrix);

/* ----------------------------------------------------------------------- */
/* DMPEXPTS Function                                                       */ 
/* This function uses a DOUBLEMATRIX to compute exp(matrix) and            */
/* records the result to the matrix itself.                                */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateDoubleMatrix function to create it        */
/* before you call this function !                                         */
/* The function will compute exp(matrix) and return PROC_SUCCESS.          */
/* ----------------------------------------------------------------------- */
int DMPEXPTS(struct DOUBLEMATRIX* pSourceMatrix);

/* ----------------------------------------------------------------------- */
/* DMPABSTS Function                                                       */ 
/* This function computes absolute value of a matrix and                   */
/* records the result to the matrix itself.                                */
/* !note: pSourceMatrix must exist .                                       */
/* Make sure you have used CreateDoubleMatrix function to create it        */
/* before you call this function !                                         */
/* The function will compute |matrix| and return PROC_SUCCESS.             */
/* ----------------------------------------------------------------------- */
int DMPABSTS(struct DOUBLEMATRIX* pSourceMatrix);

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
double DMGETAT(struct DOUBLEMATRIX* pSourceMatrix,long nRow,long nCol); 

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
int DMSETAT(struct DOUBLEMATRIX* pSourceMatrix,long nRow,long nCol,double dValue); 

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
struct DOUBLEMATRIX* DMSOLVEEQU(struct DOUBLEMATRIX* pMatrixA,struct DOUBLEMATRIX* pMatrixB);

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
double DMNORM1(struct DOUBLEMATRIX* pSourceMatrix);

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
double DMGETMAX(struct DOUBLEMATRIX* pSourceMatrix); 

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
double DMGETMIN(struct DOUBLEMATRIX* pSourceMatrix); 

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
struct BYTEMATRIX* DMCONVTOBM(struct DOUBLEMATRIX* pSourceMatrix);

/* ----------------------------------------------------------------------- */
/* DMSORTMERGE_0 Function                                                  */
/* This function is used for sorting rows of a double matrix using two-way */
/* merge algorithm.                                                        */
/* The function will return PROC_SUCCESS;                                  */
/* if any error, return PROC_FAILURE.                                      */
/* ----------------------------------------------------------------------- */
int DMSORTMERGE_0(struct DOUBLEMATRIX* pSourceMatrix, struct DOUBLEMATRIX** pSortMatrix, struct LONGMATRIX** pSortIndexMatrix);

/* ----------------------------------------------------------------------- */
/* DMSORTMERGEA_0 Function                                                 */
/* This function is used for sorting rows of a double matrix using         */
/* ameliorated two-way merge algorithm.                                    */
/* The function will return PROC_SUCCESS;                                  */
/* if any error, return PROC_FAILURE.                                      */
/* ----------------------------------------------------------------------- */
int DMSORTMERGEA_0(struct DOUBLEMATRIX* pSourceMatrix, struct DOUBLEMATRIX** pSortMatrix, struct LONGMATRIX** pSortIndexMatrix);

/* ----------------------------------------------------------------------- */
/* SORT_MERGE_DOUBLE Function                                              */
/* This function is used for sorting double array using two-way merge sort */
/* algorithm.                                                              */
/* Time complexity: n*log2(n)                                              */
/* The function will return the length of sortarray;                       */
/* if any error, return -1.                                                */
/* ----------------------------------------------------------------------- */
long SORT_MERGE_DOUBLE(double *array, long *indexarray, double **sortarray, long **sortindexarray, long nstart, long nend);

/* ----------------------------------------------------------------------- */
/* SORT_MERGE_DOUBLE_A Function                                            */
/* This function is used for sorting double array using ameliorated        */
/* two-way merge sort algorithm.                                           */
/* Time complexity: n*log2(n)                                              */
/* The function will return the length of sortarray;                       */
/* if any error, return -1.                                                */
/* ----------------------------------------------------------------------- */
long SORT_MERGE_DOUBLE_A(double *array, long *indexarray, double *sortarray, long *sortindexarray, long nstart, long nend);

/* ----------------------------------------------------------------------- */
/* SORT_MERGE_DOUBLE_MERGE_1 Function                                      */
/* This function is used for merging and sorting two double arrays.        */
/* The function will return PROC_SUCCESS;                                  */
/* if any error, return PROC_FAILURE.                                      */
/* ----------------------------------------------------------------------- */
int SORT_MERGE_DOUBLE_MERGE_1(double *array1, long nlen1, double *array2, long nlen2, double **sortarray, long nlen);

/* ----------------------------------------------------------------------- */
/* SORT_MERGE_DOUBLE_MERGE_2 Function                                      */
/* This function is used for merging and sorting two double arrays.        */
/* The function will return PROC_SUCCESS;                                  */
/* if any error, return PROC_FAILURE.                                      */
/* ----------------------------------------------------------------------- */
int SORT_MERGE_DOUBLE_MERGE_2(double *array1, long *indexarray1, long nlen1, double *array2, long *indexarray2, long nlen2, double **sortarray, long **sortindexarray, long nlen);

/* ----------------------------------------------------------------------- */
/* SORT_MERGE_DOUBLE_MERGE_A_1 Function                                    */
/* This function is used for merging and sorting two double arrays.        */
/* The function will return PROC_SUCCESS;                                  */
/* if any error, return PROC_FAILURE.                                      */
/* ----------------------------------------------------------------------- */
int SORT_MERGE_DOUBLE_MERGE_A_1(double *array1, long nlen1, double *array2, long nlen2, double *sortarray, long nlen);

/* ----------------------------------------------------------------------- */
/* SORT_MERGE_DOUBLE_MERGE_A_2 Function                                    */
/* This function is used for merging and sorting two double arrays.        */
/* The function will return PROC_SUCCESS;                                  */
/* if any error, return PROC_FAILURE.                                      */
/* ----------------------------------------------------------------------- */
int SORT_MERGE_DOUBLE_MERGE_A_2(double *array1, long *indexarray1, long nlen1, double *array2, long *indexarray2, long nlen2, double *sortarray, long *sortindexarray, long nlen);

/* ----------------------------------------------------------------------- */
/* SORT_MERGE_INT_A Function                                               */
/* This function is used for sorting int array using ameliorated           */
/* two-way merge sort algorithm.                                           */
/* Time complexity: n*log2(n)                                              */
/* The function will return the length of sortarray;                       */
/* if any error, return -1.                                                */
/* ----------------------------------------------------------------------- */
long SORT_MERGE_INT_A(int *array, long *indexarray, int *sortarray, long *sortindexarray, long nstart, long nend);

/* ----------------------------------------------------------------------- */
/* SORT_MERGE_INT_MERGE_A_1 Function                                       */
/* This function is used for merging and sorting two int arrays.           */
/* The function will return PROC_SUCCESS;                                  */
/* if any error, return PROC_FAILURE.                                      */
/* ----------------------------------------------------------------------- */
int SORT_MERGE_INT_MERGE_A_1(int *array1, long nlen1, int *array2, long nlen2, int *sortarray, long nlen);

/* ----------------------------------------------------------------------- */
/* SORT_MERGE_INT_MERGE_A_2 Function                                       */
/* This function is used for merging and sorting two int arrays.           */
/* The function will return PROC_SUCCESS;                                  */
/* if any error, return PROC_FAILURE.                                      */
/* ----------------------------------------------------------------------- */
int SORT_MERGE_INT_MERGE_A_2(int *array1, long *indexarray1, long nlen1, int *array2, long *indexarray2, long nlen2, int *sortarray, long *sortindexarray, long nlen);

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
			   struct LONGMATRIX **pDataSid);

/* ----------------------------------------------------------------------- */
/* DVNORM2 Function                                                        */
/* This function calculate the 2-norm of a double vector.                  */
/* Make sure you have allocated the required memory to the vetor pointer   */
/* before you call this function !                                         */
/* The function will return the 2-norm value;                              */
/* ----------------------------------------------------------------------- */
double DVNORM2(double *Vectorarray, long nLen);

/* ----------------------------------------------------------------------- */
/* DVNORMPRODUCT Function                                                  */
/* This function calculate the normalized product of two double vectors.   */
/* Make sure you have allocated the required and identical memories to     */
/* vetors before you call this function !                                  */
/* The function will return the product value;                             */
/* ----------------------------------------------------------------------- */
double DVNORMPRODUCT(double *Vectorarray1, double *Vectorarray2, long nLen);

/* ----------------------------------------------------------------------- */
/* BVNORM2 Function                                                        */
/* This function calculate the 2-norm of a byte vector.                    */
/* Make sure you have allocated the required memory to the vetor pointer   */
/* before you call this function !                                         */
/* The function will return the 2-norm value;                              */
/* ----------------------------------------------------------------------- */
double BVNORM2(unsigned char *Vectorarray, long nLen);

/* ----------------------------------------------------------------------- */
/* BVNORMPRODUCT Function                                                  */
/* This function calculate the normalized product of two byte vectors.     */
/* Make sure you have allocated the required and identical memories to     */
/* vetors before you call this function !                                  */
/* The function will return the product value;                             */
/* ----------------------------------------------------------------------- */
double BVNORMPRODUCT(unsigned char *Vectorarray1, unsigned char *Vectorarray2, long nLen);

/* ----------------------------------------------------------------------- */
/* pythag Function                                                         */
/* Compute (a^2+b^2)^1/2 without destructive underflow or overflow         */
/* ----------------------------------------------------------------------- */
double pythag(double a, double b);

/* ----------------------------------------------------------------------- */
/* tred2 Function                                                          */
/* Householder reduction of a real, symmetric matrix a[1:n][1:n].          */
/* On output, a is replaced by the orthogonal matrix Q effecting the       */
/* transformation. d[1:n] returns the diagonal elements of the             */
/* resulting tridiagonal matrix, and e[1:n] the off-diagonal elements,     */
/* with e[1] = 0.                                                          */
/* ----------------------------------------------------------------------- */
void tred2(double **a, int n, double d[], double e[]);

/* ----------------------------------------------------------------------- */
/* tqli Function                                                           */
/* QL algorithm with implicit shifts, to determine the eigenvalues and     */
/* eigenvectors of a real, symmetric tridiagonal matrix previously         */
/* reduced by tred2.                                                       */
/* ----------------------------------------------------------------------- */
void tqli(double d[], double e[], int n, double **z);

/* ----------------------------------------------------------------------- */
/* eigsrt Function                                                         */
/* Given the eigenvalues d[1:n] and eigenvectors v[1:n][1:n], this         */
/* function sorts the eigenvalues into descending order, and rearranges    */
/* the columns of v correspondingly.                                       */
/* ----------------------------------------------------------------------- */
void eigsrt(double d[], double **v, int n);

/* ----------------------------------------------------------------------- */
/* jacobi Function                                                         */
/* Computes eigenvalues and eigenvectors of a real symmetrix matrix        */
/* a[1:n][1:n]. On output, elements of a above the diagonal are destroyed. */
/* d[1:n] returns the eigenvalues; v[1:n][1:n] is a matrix whose columns   */
/* contain the normalized eigenvectors. nrot returns the number of jacobi  */
/* rotations that were required.                                           */
/* ----------------------------------------------------------------------- */
void jacobi(double **a, int n, double d[], double **v, int *nrot);

/* ----------------------------------------------------------------------- */
/* DMSYMEIGEN Function                                                     */
/* This function computes eigen values and eigen vectors of a symmetrix    */
/* matrix.                                                                 */
/* pX is the input matrix, pD is a vector that contains eigen values       */
/* pV is a matrix whose columns are eigen vectors with norm = 1.           */
/* ----------------------------------------------------------------------- */
int DMSYMEIGEN(struct DOUBLEMATRIX *pX,struct DOUBLEMATRIX **pD,struct DOUBLEMATRIX **pV);

/* ----------------------------------------------------------------------- */
/* DMSYMINVEIGEN Function                                                  */
/* This function computes inverse of a symmetric matrix based on its eigen */
/* decomposition. pD and pV are output of DMSYMEIGEN                       */
/* ----------------------------------------------------------------------- */
int DMSYMINVEIGEN(struct DOUBLEMATRIX *pD,struct DOUBLEMATRIX *pV, struct DOUBLEMATRIX **pI);

/* ----------------------------------------------------------------------- */
/* DMSYMEIGENQR Function                                                   */
/* This function computes eigen values and eigen vectors of a symmetrix    */
/* matrix using QR decomposition                                           */
/* pX is the input matrix, pD is a vector that contains eigen values       */
/* pV is a matrix whose columns are eigen vectors with norm = 1.           */
/* ----------------------------------------------------------------------- */
int DMSYMEIGENQR(struct DOUBLEMATRIX *pX,struct DOUBLEMATRIX **pD,struct DOUBLEMATRIX **pV);