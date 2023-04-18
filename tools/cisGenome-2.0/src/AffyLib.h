/* ----------------------------------------------------------------------- */
/*  AffyLib.h : interface of the affymetrix library                        */
/*  Author : Ji HongKai ; Time: 2004.08                                    */
/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------- */
/*                              Declaration                                */
/* ----------------------------------------------------------------------- */

#define MAX_MATCH_GAP 10000
#define TILE_PROBE_LEN 26
#define PROBE_BPMAP_LEN 7
#define DNA_BASE_A     0
#define DNA_BASE_C     1
#define DNA_BASE_G     2
#define DNA_BASE_T     3
#define DNA_BASE_OTHER 4
#define LOCAL_REPEAT_RANGE 1000

#define FLOAT_SIZE sizeof(float)
#define DOUBLE_SIZE sizeof(double)
#define INT_SIZE sizeof(int)
#define DWORD_SIZE sizeof(unsigned int)
#define ULONG_SIZE sizeof(unsigned long)
#define LONG_SIZE sizeof(long)
#define SHORT_SIZE sizeof(short)
#define USHORT_SIZE sizeof(unsigned short)
#define CHAR_SIZE sizeof(char)
#define UCHAR_SIZE sizeof(unsigned char)
#define WCHAR_SIZE sizeof(wchar_t)
#define WCHAR16_SIZE sizeof(unsigned short)

typedef unsigned short wchar_t16;

/* ----------------------------------------------------------------------- */
/*                              Structures                                 */
/* ----------------------------------------------------------------------- */
struct tagAffyGenomeAlign
{
	/* chromosome */
	char strChr[30];
	int nChr;
	/* start */
	int nStart;
	/* end */
	int nEnd;
	/* rc */
	char chRc;
	/* identity */
	double dIdentity;
	/* physical map */
	char strPhyMap[30];
	/* refseq id */
	char strRefId[30];

	/* next */
	struct tagAffyGenomeAlign *pNext;
};


struct tagAffyBpMapUnit
{
	/* X coordinate of PM */
	unsigned long nPMX;
	/* Y coordinate of PM */
	unsigned long nPMY;
	/* X coordinate of MM */
	unsigned long nMMX;
	/* Y coordinate of MM */
	unsigned long nMMY;
	/* length of PM or MM probe */
	unsigned char bProbeLen;
	/* probe sequence */
	char strProbeSeq[TILE_PROBE_LEN];
	/* match score */
	float fMatchScore;
	/* position in the sequence */
	unsigned long nPos;
	/* strand */
	unsigned char bStrand;

	/* local repeat times */
	int nRepeatNum;
	/* same postion map times */
	int nDepthNum;

	/* next */
	struct tagAffyBpMapUnit *pNext;
};

struct tagCELSubGrid
{
	/* row and column number */
	int nRows;
	int nCols;

	/* coordinate in pixels */
	float fUpperLeft[2];
	float fUpperRight[2];
	float fLowerLeft[2];
	float fLowerRight[2];

	/* cell position */
	int nCellLeft;
	int nCellTop;
	int nCellRight;
	int nCellBottom;
};

struct tagCELData
{
	/* magic number; version number; number of columns, rows and cells */
    int nMagicnumber;
	int nVersionnumber;
	int nCols;
	int nRows;
	int nNumberCells;

	/* header information */
	int nHeaderLength;
	struct tagString *vHeader;

	/* algorithm information */
	int nAlgorithmLength;
	struct tagString *vAlgorithm;

	/* algorithm parameter information */
	int nAlgorithmParamLength;
	struct tagString *vAlgorithmParam;

	/* data */
	struct DOUBLEMATRIX *pIntensity;
	struct DOUBLEMATRIX *pSD;
	struct INTMATRIX *pPixelNum;

	/* masked cells */
	unsigned int nMaskedCells;
	struct INTMATRIX *pMaskedX;
	struct INTMATRIX *pMaskedY;

	/* outlier cells */
	unsigned int nOutlierCells;
	struct INTMATRIX *pOutlierX;
	struct INTMATRIX *pOutlierY;

	/* v3 info */
	int nTotalX;
	int nTotalY;
	int nOffsetX;
	int nOffsetY;
	int nUL[2];
	int nUR[2];
	int nLR[2];
	int nLL[2];
	int nInvertX, nInvertY;
	int swapXY;
	char strDatHeader[MED_LINE_LENGTH];

	unsigned int nModifiedCells;
	struct INTMATRIX *pModifiedX;
	struct INTMATRIX *pModifiedY;
	struct DOUBLEMATRIX *pModifiedOrig;

	/* v4 info */
    int nCellMargin;
	int nSubGrids;
	struct tagCELSubGrid **vSubGrids;
};

struct tagBARSeq
{
	/* sequence name */
	struct tagString *pSeqName;

	/* sequence group name */
	struct tagString *pSeqGroupName;
	
	/* sequence version */
	struct tagString *pSeqVersion;

	/* paramenter name/value pairs */
	int nParamNum;
	struct tagString **vParamName;
	struct tagString **vParamValue;

	/* data points */
	int nDataNum;
	int nColNum;
	struct DOUBLEMATRIX **vData;
};

struct tagBARData
{
	/* magic number; version number */
    char strMagicnumber[9];
	float fVersionnumber;
	
	/* number of sequences and columns, */
	int nSeqNum;
	int nColNum;

	/* column type */
	struct INTMATRIX *pFieldType;

	/* paramenter name/value pairs */
	int nParamNum;
	struct tagString **vParamName;
	struct tagString **vParamValue;

	/* sequence information */
	struct tagBARSeq **vSeqData;
};

struct tagGenericDataSet
{
	/* The file position of the first data element in the data set. 
	This is the first byte after the data set header. */
	unsigned int nFirstDataElementPos;
	/* The file position of the next data set within the data group. 
	When this is the last data set in the data group the value shall 
	be 1 byte past the end of the data set. This way the size of the 
	data set may be determined.*/
	unsigned int nNextDataSetPos;
	/* The data set name. */
	struct tagWString *pDataSetName;
	/* The number of name/value/type parameters. */
	int nParamNum;
	/* Array of name/value/type parameters. (WSTRING / VALUE / TYPE) [ ] */
	struct tagWString **vParamName;
	struct tagString **vParamValue;
	struct tagWString **vParamType;
	/* Number of columns in the data set.Example: For expression arrays, 
	columns may include signal, p-value, detection call and for genotyping 
	arrays columns may include allele call, and confidence value. For universal 
	arrays, columns may include probe set intensities and background. UINT */
	unsigned int nColNum;
	/* An array of column names, column value types and column type sizes (one per column).
	The value type shall be represented by the value from the value type table. 
	The size shall be the size of the type in bytes. For strings, this value shall be 
	the size of the string in bytes plus 4 bytes for the string length written before 
	the string in the file. (WSTRING / BYTE / INT) [ ] */
	struct tagWString **vColName;
	unsigned char *vColType;
	int *vColSize;
	/* The number of rows in the data set. UINT */
	unsigned int nRowNum;
	/* The data set table, consisting of rows of columns (data values). 
	The specific type and size of each column is described by the data and size types above. ROW [ ] */
	char *vData;
};

struct tagGenericDataGroup
{
	/* File position of the next data group. When this is the last data 
	group in the file, the value should be 0. UINT */
	unsigned int nNextDataGroupPos;
	/* File position of the first data set within the data group. UINT */
	unsigned int nFirstDataSetPos;
	/* The number of data sets within the data group. INT */
	int nDataSetNum;
	/* The data group name. WSTRING */
	struct tagWString *pDataGroupName;
	/* Data Sets */
	struct tagGenericDataSet **vDataSets;
};

struct tagGenericDataHeader
{
	 /* The data type identifier. This is used to identify the type of data stored in the file. For example:
		acquisition data (affymetrix-calvin-scan-acquisition) 
		intensity data (tbd) 
		expression results (tbd) 
		genotyping results (tbd) GUID
	*/
	struct tagString *pDataTypeID;
	
	/* Unique file identifier. This is the identifier to use to link the file with parent files. 
	This identifier will be updated whenever the contents of the file change.
	Example: When a  user manually aligns the grid in a DAT file the grid coordinates are updated in the DAT 
	file and the file is given a new file identifier. GUID */
	struct tagString *pUniqueFileID;

	/* Date and time of file creation. DATETIME */
	struct tagWString *pDateTime;

	/* The locale of the operating system that the file was created on. LOCALE */
	struct tagWString *pLocale;

	/* The number of name/type/value parameters. INT */
	int nParamNum;

	/* Array of parameters stored as name/value/type triplets. (WSTRING / VALUE / TYPE ) [ ] */
	struct tagWString **vParamName;
	struct tagString **vParamValue;
	struct tagWString **vParamType;
	 
	/* Number of parent file headers. INT */
	int nParentNum;

	/* Array of parent file headers. Generic Data Header [ ] */
	struct tagGenericDataHeader **vParent;
};

struct tagGenericFileHeader
{
	/* Magic number. A value to identify that this is a Command Console data file. 
	The value will be fixed to 59. UBYTE */
	unsigned char nMagicNumber;
	/* The version number of the file. This is the version of the file format. 
	It is currently fixed to 1. */
	unsigned char nVersionNumber;
	/* The number of data groups. INT */
	int nDataGroupNum;
	/* File position of the first data group. UINT */
	unsigned int nFirstDataGroupPos;
};

struct tagGenericFileObject
{
	char strFilePath[MED_LINE_LENGTH];
	struct tagGenericFileHeader *pFileHeader;
	struct tagGenericDataHeader *pDataHeader;
	struct tagGenericDataGroup **vDataGroup;
};

/* ----------------------------------------------------------------------- */ 
/*                                Functions                                */
/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------- */ 
/* Load BIG_ENDIAN and LITTLE_ENDIAN files.                                */
/* ----------------------------------------------------------------------- */ 
void reverse_buf(char *buf, int size);

size_t big_endian_fread(void *ptr, size_t size, size_t count, FILE *stream, 
						int little_endian_machine);

size_t big_endian_fwrite(const void *ptr, size_t size, size_t count, FILE *stream,
						 int little_endian_machine);

size_t little_endian_fread(void *ptr, size_t size, size_t count, FILE *stream, 
						int little_endian_machine);

size_t little_endian_fwrite(const void *ptr, size_t size, size_t count, FILE *stream,
						 int little_endian_machine);


/* ----------------------------------------------------------------------- */ 
/*  Affy_CELData_Create()                                                  */
/*  create a CELData object.                                               */
/* ----------------------------------------------------------------------- */ 
struct tagCELData *Affy_CELData_Create();

/* ----------------------------------------------------------------------- */ 
/*  Affy_CELData_Destroy()                                                 */
/*  delete a CELData object.                                               */
/* ----------------------------------------------------------------------- */ 
void Affy_CELData_Destroy(struct tagCELData **pCELData);

/* ----------------------------------------------------------------------- */ 
/*  Affy_CELSubGrid_Create()                                               */
/*  create a CELSubGrid object.                                            */
/* ----------------------------------------------------------------------- */ 
struct tagCELSubGrid *Affy_CELSubGrid_Create();

/* ----------------------------------------------------------------------- */ 
/*  Affy_CELSubGrid_Destroy()                                              */
/*  delete a CELSubGrid object.                                            */
/* ----------------------------------------------------------------------- */ 
void Affy_CELSubGrid_Destroy(struct tagCELSubGrid **pCELSubGrid);

/* ----------------------------------------------------------------------- */ 
/*  Affy_LoadCELv4()                                                       */
/*  Loading raw data from a single affymetrix's *.CEL (v4) file            */
/* ----------------------------------------------------------------------- */ 
struct tagCELData *Affy_LoadCELv4(char strFileName[]);

/* ----------------------------------------------------------------------- */ 
/*  Affy_LoadCELv4_Fast()                                                  */
/*  Loading raw data from a single affymetrix's *.CEL (v4) file            */
/* ----------------------------------------------------------------------- */ 
struct tagCELData *Affy_LoadCELv4_Fast(char strFileName[]);

/* ----------------------------------------------------------------------- */ 
/*  Affy_SaveCELv4()                                                       */
/*  Saving raw data to a affymetrix's *.CEL (v4) file                      */
/* ----------------------------------------------------------------------- */ 
int Affy_SaveCELv4(char strFileName[], struct tagCELData *pCELData);

/* ----------------------------------------------------------------------- */ 
/*  Affy_LoadCEL_CmdCslv1()                                                */
/*  Loading raw data from a single affymetrix's command console *.CEL file */
/* ----------------------------------------------------------------------- */ 
struct tagCELData *Affy_LoadCEL_CmdCslv1(char strFileName[]);

/* ----------------------------------------------------------------------- */ 
/*  Affy_LoadCELv3()                                                       */
/*  Loading raw data from a single affymetrix's *.CEL (v3) file            */
/* ----------------------------------------------------------------------- */ 
struct tagCELData *Affy_LoadCELv3(char strFileName[]);

/* ----------------------------------------------------------------------- */ 
/*  Affy_SaveCELv3()                                                       */
/*  Exporting raw data from a single affymetrix's *.CEL (v3) file          */
/* ----------------------------------------------------------------------- */ 
int Affy_SaveCELv3(char strFileName[], struct tagCELData *pCELData);

/* ----------------------------------------------------------------------- */ 
/*  Affy_BARSeq_Create()                                                   */
/*  create a BARSeq object.                                                */
/* ----------------------------------------------------------------------- */ 
struct tagBARSeq *Affy_BARSeq_Create();

/* ----------------------------------------------------------------------- */ 
/*  Affy_BARSeq_Destroy()                                                  */
/*  delete a BARSeq object.                                                */
/* ----------------------------------------------------------------------- */ 
void Affy_BARSeq_Destroy(struct tagBARSeq **pBARSeq);

/* ----------------------------------------------------------------------- */ 
/*  Affy_BARData_Create()                                                  */
/*  create a BARData object.                                               */
/* ----------------------------------------------------------------------- */ 
struct tagBARData *Affy_BARData_Create();

/* ----------------------------------------------------------------------- */ 
/*  Affy_BARData_Destroy()                                                 */
/*  delete a BARData object.                                               */
/* ----------------------------------------------------------------------- */ 
void Affy_BARData_Destroy(struct tagBARData **pBARData);

/* ----------------------------------------------------------------------- */ 
/*  Affy_LoadBAR()                                                         */
/*  Loading raw data from a single affymetrix's *.bar file                 */
/* ----------------------------------------------------------------------- */ 
struct tagBARData *Affy_LoadBAR(char strFileName[]);

/* ----------------------------------------------------------------------- */ 
/*  Affy_LoadBAR_Fast()                                                    */
/*  Loading raw data from a single affymetrix's *.bar file                 */
/* ----------------------------------------------------------------------- */ 
struct tagBARData *Affy_LoadBAR_Fast(char strFileName[]);

/* ----------------------------------------------------------------------- */ 
/*  Affy_SaveBAR()                                                         */
/*  Saving raw data to a affymetrix's *.bar file                           */
/* ----------------------------------------------------------------------- */ 
int Affy_SaveBAR(char strFileName[], struct tagBARData *pBARData);

/* ----------------------------------------------------------------------- */ 
/*  Affy_SaveBAR_Columns_Fast()                                            */
/*  Saving specified columns of a bar data object to a *.bar file          */
/* ----------------------------------------------------------------------- */ 
int Affy_SaveBAR_Columns_Fast(char strFileName[], struct tagBARData *pBARData,
							  struct INTMATRIX *pCol);

/* ----------------------------------------------------------------------- */ 
/*  Affy_SaveFilteredBAR_Columns_Fast()                                    */
/*  Saving specified columns of a bar data object to a *.bar file          */
/*  Exporting data will be filtered according to GenomeGrp name.           */
/* ----------------------------------------------------------------------- */ 
int Affy_SaveFilteredBAR_Columns_Fast(char strFileName[], struct tagBARData *pBARData,
							  struct INTMATRIX *pCol, char strGenomeGrp[]);

/* ----------------------------------------------------------------------- */ 
/*  Affy_BAR2TXT()                                                         */
/*  Convert *.bar file to *.txt file.                                      */
/* ----------------------------------------------------------------------- */ 
int Affy_BAR2TXT(char strBARFile[], char strTXTFile[]);

/* ----------------------------------------------------------------------- */ 
/*  Affy_BAR2TXT_Fast()                                                    */
/*  Convert *.bar file to *.txt file.                                      */
/* ----------------------------------------------------------------------- */ 
int Affy_BAR2TXT_Fast(char strBARFile[], char strTXTFile[]);

/* ----------------------------------------------------------------------- */ 
/*  Affy_BAR2WIG_Main()                                                    */
/*  Convert BAR file to WIG file.                                          */
/* ----------------------------------------------------------------------- */ 
int Affy_BAR2WIG_Main(char strInputFile[], char strOutputFile[]);

/* ----------------------------------------------------------------------- */ 
/*  Affy_BARSeq_Clone()                                                    */
/*  clone a BARSeq object to a new BARData object.                         */
/* ----------------------------------------------------------------------- */ 
struct tagBARSeq *Affy_BARSeq_Clone(struct tagBARSeq *pBARSeq);

/* ----------------------------------------------------------------------- */ 
/*  Affy_BARData_Clone()                                                   */
/*  clone a BARData object to a new BARData object.                        */
/* ----------------------------------------------------------------------- */ 
struct tagBARData *Affy_BARData_Clone(struct tagBARData *pBARData);

/* ----------------------------------------------------------------------- */ 
/*  Affy_BARData_AddTS()                                                   */
/*  Add a number with given columns of a BARData object.                   */
/* ----------------------------------------------------------------------- */ 
int Affy_BARData_AddTS(struct tagBARData *pBARData, double dLambda,
		struct INTMATRIX *pColIndicator);

/* ----------------------------------------------------------------------- */ 
/*  Affy_BARData_SubTS()                                                   */
/*  subtract a number by given columns of a BARData object, i.e.           */
/*  dLambda-pBARData[nCol]                                                 */
/* ----------------------------------------------------------------------- */ 
int Affy_BARData_SubTS(struct tagBARData *pBARData, double dLambda,
		struct INTMATRIX *pColIndicator);

/* ----------------------------------------------------------------------- */ 
/*  Affy_BARData_MultTS()                                                  */
/*  multiply a number with given columns of a BARData object.              */
/* ----------------------------------------------------------------------- */ 
int Affy_BARData_MultTS(struct tagBARData *pBARData, double dLambda,
		struct INTMATRIX *pColIndicator);

/* ----------------------------------------------------------------------- */ 
/*  Affy_BARData_RecTS()                                                   */
/*  divide a number by given columns of a BARData object.                  */
/* ----------------------------------------------------------------------- */ 
int Affy_BARData_RecTS(struct tagBARData *pBARData, double dLambda,
		struct INTMATRIX *pColIndicator);

/* ----------------------------------------------------------------------- */ 
/*  Affy_BARData_LogTS()                                                   */
/*  take logarithm of given columns of a BARData object.                   */
/* ----------------------------------------------------------------------- */ 
int Affy_BARData_LogTS(struct tagBARData *pBARData, double dBase, 
		struct INTMATRIX *pColIndicator);

/* ----------------------------------------------------------------------- */ 
/*  Affy_BARData_PowTS()                                                   */
/*  take power of given columns of a BARData object.					   */
/* ----------------------------------------------------------------------- */ 
int Affy_BARData_PowTS(struct tagBARData *pBARData, double dPow, 
		struct INTMATRIX *pColIndicator);

/* ----------------------------------------------------------------------- */ 
/*  Affy_BARData_UPowTS()                                                  */
/*  take dBase^power where power = given columns of a BARData object.      */
/* ----------------------------------------------------------------------- */ 
int Affy_BARData_UPowTS(struct tagBARData *pBARData, double dBase, 
		struct INTMATRIX *pColIndicator);

/* ----------------------------------------------------------------------- */ 
/*  Affy_BARData_ExpTS()                                                   */
/*  take exponential of given columns of a BARData object.                 */
/* ----------------------------------------------------------------------- */ 
int Affy_BARData_ExpTS(struct tagBARData *pBARData, struct INTMATRIX *pColIndicator);

/* ----------------------------------------------------------------------- */ 
/*  Affy_BARData_AbsTS()                                                   */
/*  take absolute value of given columns of a BARData object.              */
/* ----------------------------------------------------------------------- */ 
int Affy_BARData_AbsTS(struct tagBARData *pBARData, struct INTMATRIX *pColIndicator);

/* ----------------------------------------------------------------------- */ 
/*  Affy_BARData_Add()                                                     */
/*  Add specified columns of two BARData objects and save results to one   */
/*  of the two BARData objects as specified.                               */
/*  If nSaveTo == 0, save results to pD1;                                  */
/*  If nSaveTo == 1, save results to pD2.                                  */
/* ----------------------------------------------------------------------- */ 
int Affy_BARData_Add(struct tagBARData *pD1, struct tagBARData *pD2,
		struct INTMATRIX *pColIndicator, int nSaveTo);

/* ----------------------------------------------------------------------- */ 
/*  Affy_BARData_Sub()                                                     */
/*  Subtract specified columns of two BARData objects and save results to  */
/*  one of the two BARData objects as specified.                           */
/*  return pD1-pD2                                                         */
/*  pColIndicator == 1 specifies which column will be involved in the      */
/*   operation;                                                            */
/*  pColIndicator == 2 specifies which column will be used as an internal  */
/*   control to match the two BARData objects so that |pD1-pD2|<1e-6 .     */
/*  If nSaveTo == 0, save results to pD1;                                  */
/*  If nSaveTo == 1, save results to pD2.                                  */
/* ----------------------------------------------------------------------- */ 
int Affy_BARData_Sub(struct tagBARData *pD1, struct tagBARData *pD2,
		struct INTMATRIX *pColIndicator, int nSaveTo);

/* ----------------------------------------------------------------------- */ 
/*  Affy_BARData_Mult()                                                    */
/*  Multiply specified columns of two BARData objects and save results to  */
/*  one of the two BARData objects as specified.                           */
/*  pColIndicator == 1 specifies which column will be involved in the      */
/*   operation;                                                            */
/*  pColIndicator == 2 specifies which column will be used as an internal  */
/*   control to match the two BARData objects so that |pD1-pD2|<1e-6 .     */
/*  If nSaveTo == 0, save results to pD1;                                  */
/*  If nSaveTo == 1, save results to pD2.                                  */
/* ----------------------------------------------------------------------- */ 
int Affy_BARData_Mult(struct tagBARData *pD1, struct tagBARData *pD2,
		struct INTMATRIX *pColIndicator, int nSaveTo);
	
/* ----------------------------------------------------------------------- */ 
/*  Affy_BARData_Div()                                                     */
/*  Divide specified columns of two BARData objects and save results to    */
/*  one of the two BARData objects as specified.                           */
/*  return pD1/pD2                                                         */
/*  pColIndicator == 1 specifies which column will be involved in the      */
/*   operation;                                                            */
/*  pColIndicator == 2 specifies which column will be used as an internal  */
/*   control to match the two BARData objects so that |pD1-pD2|<1e-6 .     */
/*  If nSaveTo == 0, save results to pD1;                                  */
/*  If nSaveTo == 1, save results to pD2.                                  */
/* ----------------------------------------------------------------------- */ 
int Affy_BARData_Div(struct tagBARData *pD1, struct tagBARData *pD2,
		struct INTMATRIX *pColIndicator, int nSaveTo);

/* ----------------------------------------------------------------------- */ 
/*  Affy_BARData_Pow()                                                     */
/*  Compute pD1.^pD2 for specified columns of two BARData objects and save */
/*  results to one of the two BARData objects as specified.                */
/*  pColIndicator == 1 specifies which column will be involved in the      */
/*   operation;                                                            */
/*  pColIndicator == 2 specifies which column will be used as an internal  */
/*   control to match the two BARData objects so that |pD1-pD2|<1e-6 .     */
/*  If nSaveTo == 0, save results to pD1;                                  */
/*  If nSaveTo == 1, save results to pD2.                                  */
/* ----------------------------------------------------------------------- */ 
int Affy_BARData_Pow(struct tagBARData *pD1, struct tagBARData *pD2,
		struct INTMATRIX *pColIndicator, int nSaveTo);

/* ----------------------------------------------------------------------- */ 
/*  Affy_CSVANNOT_To_GenomeAlign_200506()                                  */
/*  convert affy csv annotation file to genome mapping coordinates.        */
/* ----------------------------------------------------------------------- */ 
int Affy_CSVANNOT_To_GenomeAlign_200506(char strInFile[], char strOutFile[], char strSpecies[]);

/* ----------------------------------------------------------------------- */ 
/*  Affy_CSVANNOT_To_Reduced_200408()                                      */
/*  convert affy csv annotation file to reduced annotation                 */
/* ----------------------------------------------------------------------- */ 
int Affy_CSVANNOT_To_Reduced_200408(char strInFile[], char strOutFile[], char strSpecies[]);

/* ----------------------------------------------------------------------- */ 
/*  AffyGenomeAlignCreate()                                                */
/*  create genome alignment element                                        */
/* ----------------------------------------------------------------------- */ 
struct tagAffyGenomeAlign *AffyGenomeAlignCreate();

/* ----------------------------------------------------------------------- */ 
/*  AffyGenomeAlignDestroy()                                               */
/*  create genome alignment element                                        */
/* ----------------------------------------------------------------------- */ 
int AffyGenomeAlignDestroy(struct tagAffyGenomeAlign *pGA);

/* ----------------------------------------------------------------------- */ 
/*  AffyGenomeAlignInit()                                                  */
/*  init alignment element                                                 */
/* ----------------------------------------------------------------------- */ 
int AffyGenomeAlignInit(struct tagAffyGenomeAlign *pGA, char strParam[], char strSpecies[]);

/* ----------------------------------------------------------------------- */ 
/*  AffyGenomeTransInit()                                                  */
/*  init alignment element                                                 */
/* ----------------------------------------------------------------------- */ 
int AffyGenomeTransInit(struct tagAffyGenomeAlign *pGA, char strParam[], char strSpecies[]);

/* ----------------------------------------------------------------------- */ 
/*  Affy_PickRandomControlFromCSV()                                        */
/*  pick random control probesets from affymetrix's CSV annotation file.   */
/* ----------------------------------------------------------------------- */ 
int Affy_PickRandomControlFromCSV(char strInFile[], char strOutFile[], int nControlNum);

/* ----------------------------------------------------------------------- */ 
/*  Affy_MapProbesetToRefSeq_1()                                           */
/*  map probeset to refseq.                                                */
/* ----------------------------------------------------------------------- */ 
int Affy_LinkProbesetToAnnotation(char strProbesetFile[], char strAnnotationFile[], 
								  char strOutputFile[], int nClassId);

/* ----------------------------------------------------------------------- */ 
/*  Affy_CreateProbeRefGeneMap()                                           */
/*  map probeset to refseq.                                                */
/* ----------------------------------------------------------------------- */ 
int Affy_CreateProbeRefGeneMap(char strProbeLists[], char strRefGenePath[], 
						char strMapFile[], char strOutFile[], char strSpecies[]);

/* ----------------------------------------------------------------------- */ 
/*  Affy_CreateProbeRefGeneMap_BothDirection()                             */
/*  map probeset to refseq.                                                */
/* ----------------------------------------------------------------------- */ 
int Affy_CreateProbeRefGeneMap_BothDirection(char strProbeLists[], char strRefGenePath[], 
						char strOutFile[], char strSpecies[]);

/* ----------------------------------------------------------------------- */ 
/*  Affy_LinkRefGeneScoreBackToProbeset()                                  */
/*  map refseq back to probeset.                                           */
/* ----------------------------------------------------------------------- */ 
int Affy_LinkRefGeneScoreBackToProbeset(char strProbeLists[], char strProbeRefGeneMapFile[],
										char strInFile[], char strOutFile[]);


/* ----------------------------------------------------------------------- */ 
/*  Affy_LoadBar_Intensity()                                               */
/*  load intensity data from *.bar file                                    */
/* ----------------------------------------------------------------------- */ 
int Affy_LoadBar_Intensity(char strBarFile[], char strOutPath[]);

/* ----------------------------------------------------------------------- */ 
/*  Affy_LoadBar_Intensity_Group()                                         */
/*  load a group of intensity data from *.bar file                         */
/* ----------------------------------------------------------------------- */ 
int Affy_LoadBar_Intensity_Group(char strDataPath[], char strBarFileList[], char strOutPath[]);

/* ----------------------------------------------------------------------- */ 
/*  Affy_BPMAPFilter_GTrans()                                              */
/*  Filter out repeat and bad probes with bpmap informaiton.               */
/* ----------------------------------------------------------------------- */ 
int Affy_BPMAPFilter_GTrans(char strInFile[], char strRefFile[], char strOutFile[]);

/* ----------------------------------------------------------------------- */ 
/*  Affy_BPMAPFilter_ORI()                                                 */
/*  Filter out repeat and bad probes with bpmap informaiton.               */
/* ----------------------------------------------------------------------- */ 
int Affy_BPMAPFilter_ORI(char strInFile[], char strRefFile[], char strOutFile[]);

/* ----------------------------------------------------------------------- */ 
/*  Affy_BPMAPFilter_DATA()                                                 */
/*  Filter out repeat and bad probes with bpmap informaiton.               */
/* ----------------------------------------------------------------------- */ 
int Affy_BPMAPFilter_DATA(char strInFile[], char strRefFile[], char strOutFile[]);


/* ----------------------------------------------------------------------- */ 
/*  Affy_LoadBPMAP()                                                       */
/*  load bpmap data from *.bar file.                                       */
/* ----------------------------------------------------------------------- */ 
int Affy_LoadBPMAP(char strInFile[], char strOutFile[]);

/* ----------------------------------------------------------------------- */ 
/*  Affy_LoadBPMAP_TileMap()                                               */
/*  load bpmap data from *.bar file, tilemap format                        */
/* ----------------------------------------------------------------------- */ 
int Affy_LoadBPMAP_TileMap(char strInFile[], char strPosFile[], 
						   char strMaskFile[], int *pnTotalProbeNum);


/* ----------------------------------------------------------------------- */ 
/*  AffyBpMapUnitCreate()                                                  */
/*  create bpmap unit.                                                     */
/* ----------------------------------------------------------------------- */ 
struct tagAffyBpMapUnit *AffyBpMapUnitCreate();

/* ----------------------------------------------------------------------- */ 
/*  AffyBpMapUnitDestroy()                                                 */
/*  delete bpmap unit.                                                     */
/* ----------------------------------------------------------------------- */ 
int AffyBpMapUnitDestroy(struct tagAffyBpMapUnit *pUnit);

/* ----------------------------------------------------------------------- */ 
/*  AffyBpMapUnitLoad()                                                    */
/*  load bpmap unit.                                                       */
/* ----------------------------------------------------------------------- */ 
int AffyBpMapUnitLoad(struct tagAffyBpMapUnit *pUnit, FILE *fpIn);

/* ----------------------------------------------------------------------- */ 
/*  AffyBpMapUnitLoad_v3()                                                 */
/*  load bpmap unit.                                                       */
/* ----------------------------------------------------------------------- */ 
int AffyBpMapUnitLoad_v3(struct tagAffyBpMapUnit *pUnit, FILE *fpIn, int nPMType);

/* ----------------------------------------------------------------------- */ 
/*  AffyBpMapUnitLoad_v3m2()                                               */
/*  load bpmap unit.                                                       */
/* ----------------------------------------------------------------------- */ 
int AffyBpMapUnitLoad_v3m2(struct tagAffyBpMapUnit *pUnit, FILE *fpIn, 
						   int nPMType, int little_endian_machine);

/* ----------------------------------------------------------------------- */ 
/*  Functions for loading affy binary files.                               */
/* ----------------------------------------------------------------------- */ 
void AFFYBAR_READ_INT(FILE *fpIn, int *value);
void AFFYBAR_READ_ULONG(FILE *fpIn, unsigned long *value);
void AFFYBAR_READ_FLOAT(FILE *fpIn, float *value);
void AFFYBAR_READ_VERSION(FILE *fpIn, float *value);

/* ----------------------------------------------------------------------- */ 
/*  Affy_LoadDatabase_Affy2Refid()                                         */
/*  Load affy probe id to refid map.                                       */
/* ----------------------------------------------------------------------- */ 
struct tagStringPair **Affy_LoadDatabase_Affy2Refid(char strDatabasePath[], 
				int *pPairNum);

/* ----------------------------------------------------------------------- */ 
/*  Affy_LoadDatabase_ExpressionData()                                     */
/*  Load expression data for probe selection                               */
/* ----------------------------------------------------------------------- */ 
struct tagStringPair **Affy_LoadDatabase_ExpressionData(char strDatabasePath[], 
				int *pPairNum, char strHeader[]);

/* ----------------------------------------------------------------------- */ 
/*  Affy_GenericFileObject_Create()                                        */
/*  Create a object for loading generic file object.                       */
/* ----------------------------------------------------------------------- */ 
struct tagGenericFileObject *Affy_GenericFileObject_Create();

/* ----------------------------------------------------------------------- */ 
/*  Affy_GenericFileObject_Delete()                                        */
/*  Delete an object for loading generic file object.                      */
/* ----------------------------------------------------------------------- */ 
void Affy_GenericFileObject_Delete(struct tagGenericFileObject **pGenericObj);

/* ----------------------------------------------------------------------- */ 
/*  Affy_GenericFileHeader_Create()                                        */
/*  Create a generic file header object.                                   */
/* ----------------------------------------------------------------------- */ 
struct tagGenericFileHeader *Affy_GenericFileHeader_Create();

/* ----------------------------------------------------------------------- */ 
/*  Affy_GenericFileHeader_Delete()                                        */
/*  Delete a generic file header object.                                   */
/* ----------------------------------------------------------------------- */ 
void Affy_GenericFileHeader_Delete(struct tagGenericFileHeader **pHeader);

/* ----------------------------------------------------------------------- */ 
/*  Affy_GenericDataHeader_Create()                                        */
/*  Create a generic data header object.                                   */
/* ----------------------------------------------------------------------- */ 
struct tagGenericDataHeader *Affy_GenericDataHeader_Create();

/* ----------------------------------------------------------------------- */ 
/*  Affy_GenericDataHeader_Delete()                                        */
/*  Delete a generic  data header object.                                  */
/* ----------------------------------------------------------------------- */ 
void Affy_GenericDataHeader_Delete(struct tagGenericDataHeader **pHeader);

/* ----------------------------------------------------------------------- */ 
/*  Affy_GenericDataGroup_Create()                                         */
/*  Create a generic file data group object.                               */
/* ----------------------------------------------------------------------- */ 
struct tagGenericDataGroup *Affy_GenericDataGroup_Create();

/* ----------------------------------------------------------------------- */ 
/*  Affy_GenericDataGroup_Delete()                                         */
/*  Delete a generic file data group object.                               */
/* ----------------------------------------------------------------------- */ 
void Affy_GenericDataGroup_Delete(struct tagGenericDataGroup **pDataGroup);

/* ----------------------------------------------------------------------- */ 
/*  Affy_GenericDataSet_Create()                                           */
/*  Create a generic file data set object.                                 */
/* ----------------------------------------------------------------------- */ 
struct tagGenericDataSet *Affy_GenericDataSet_Create();

/* ----------------------------------------------------------------------- */ 
/*  Affy_GenericDataSet_Delete()                                           */
/*  Delete a generic file data set object.                                 */
/* ----------------------------------------------------------------------- */ 
void Affy_GenericDataSet_Delete(struct tagGenericDataSet **pDataSet);

/* ----------------------------------------------------------------------- */ 
/*  Affy_GenericFileObject_Load()                                          */
/*  Load a generic file object from a file.                                */
/* ----------------------------------------------------------------------- */ 
struct tagGenericFileObject *Affy_GenericFileObject_Load(char strFilePath[]);

/* ----------------------------------------------------------------------- */ 
/*  Affy_GenericFile_LoadFileHeader()                                      */
/*  Load generic file header.                                              */
/* ----------------------------------------------------------------------- */ 
struct tagGenericFileHeader *Affy_GenericFile_LoadFileHeader(FILE *fpIn);

/* ----------------------------------------------------------------------- */ 
/*  Affy_GenericFile_LoadDataGroup()                                       */
/*  Load generic data group.                                               */
/* ----------------------------------------------------------------------- */ 
struct tagGenericDataHeader *Affy_GenericFile_LoadDataHeader(FILE *fpIn);

/* ----------------------------------------------------------------------- */ 
/*  Affy_GenericFile_LoadDataGroup()                                       */
/*  Load generic data group.                                               */
/* ----------------------------------------------------------------------- */ 
struct tagGenericDataGroup *Affy_GenericFile_LoadDataGroup(FILE *fpIn);

/* ----------------------------------------------------------------------- */ 
/*  Affy_GenericFile_LoadDataSet()                                         */
/*  Load generic data set.                                                 */
/* ----------------------------------------------------------------------- */ 
struct tagGenericDataSet *Affy_GenericFile_LoadDataSet(FILE *fpIn);

/* ----------------------------------------------------------------------- */ 
/*  Affy_GenericFile_LoadWString()                                         */
/*  Load wchar_t string from generic file.                                 */
/* ----------------------------------------------------------------------- */ 
struct tagWString *Affy_GenericFile_LoadWString(FILE *fpIn, int little_endian_machine);

/* ----------------------------------------------------------------------- */ 
/*  Affy_GenericFile_LoadString()                                          */
/*  Load 1-byte character string from generic file.                        */
/* ----------------------------------------------------------------------- */ 
struct tagString *Affy_GenericFile_LoadString(FILE *fpIn, int little_endian_machine);

/* ----------------------------------------------------------------------- */ 
/*  Affy_MIME_Value2Int()                                                  */
/*  Convert MIME Value to Int.                                             */
/* ----------------------------------------------------------------------- */ 
int Affy_MIME_Value2Int(char *strValue, int nBufferLen);

/* ----------------------------------------------------------------------- */ 
/*  Affy_MIME_Value2Float()                                                */
/*  Convert MIME Value to Float.                                           */
/* ----------------------------------------------------------------------- */ 
float Affy_MIME_Value2Float(char *strValue, int nBufferLen);

/* ----------------------------------------------------------------------- */ 
/*  Affy_ParamValueType2String()                                           */
/*  Convert affymetrix Param/Value/Type to a 1 byte character string.      */
/* ----------------------------------------------------------------------- */ 
struct tagString *Affy_ParamValueType2String(wchar_t strParamName[], char strParamValue[], 
	int nBufferLen, wchar_t strParamType[]);