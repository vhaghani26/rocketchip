/* ----------------------------------------------------------------------- */
/*  SequenceLib.h : interface of the sequence library                      */
/*  Author : Ji HongKai ; Time: 2004.07                                    */
/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------- */
/*                              Declaration                                */
/* ----------------------------------------------------------------------- */
#define FASTA_LINE_LEN 60
/* the length of sequence names */
#define ALIAS_LENGTH 255

/* ----------------------------------------------------------------------- */ 
/*                                 Struct                                  */
/* ----------------------------------------------------------------------- */ 

/* tagSequence is used for manipulating (nucleotide) sequences. For example,    */
/* extract a segment from a given sequence, get complementary sequence etc.     */
struct tagSequence
{
	/* m_strAlias is a alias of the sequence. */
	char m_strAlias[ALIAS_LENGTH];
	/* m_nIndex is a index for this sequence. */
	int m_nIndex;
	/* m_nLength is the length of the sequence. */
	int m_nLength;
	/* m_pSequence is the real sequence. */
	struct tagString *m_pSequence;
	/* m_nStart and m_nEnd is the position in the segment named by alias */
	int m_nStart;
	int m_nEnd;
	/* m_pNext is the pointer to the next sequence */
	struct tagSequence *m_pNext;
};

/* tagSegment is used for manipulating HSP -- result of BLAST.                  */
struct tagSegment
{
	/* m_strAlias is a alias of the source sequence. */
	char m_strAlias[LINE_LENGTH];
	/* m_nDirection is direction of the segment, when compared with another one.*/
	int m_nDirection;
	/* m_nStart and m_nEnd are used for indicating the range of the segment.    */
	int m_nStart;
	int m_nEnd;
	/* m_pSequence is the sequence in the alignment. */
	struct tagString *m_pSequence;
};

/* tagBlastResult is used for recording the align result of BLAST.              */
struct tagBlastResult
{
	/* m_QuerySeg and m_SbjctSeg are query segment and subject segment          */
	/* respectively.                                                            */
	struct tagSegment m_QuerySeg;
	struct tagSegment m_SbjctSeg;
	/* m_dScore is the alignment score. */
	double m_dScore;
	/* m_dExpect is the E-Value of alignment. */
	double m_dExpect;
	/* m_nIdentitiesLength is the length of identities. */
	int m_nIdentitiesLength;
	/* m_dIdentitiesRatio is the ratio of identities. */
	double m_dIdentitiesRatio;
	/* m_pNextResult is the pointer to next result node. */
	struct tagBlastResult *m_pNextResult;
	/* m_nCommonCount is the count of common sequences */
	int m_nCommonCount;
};


/* ----------------------------------------------------------------------- */ 
/*                                Functions                                */
/* ----------------------------------------------------------------------- */ 

/* new and delete actions */
struct tagSequence *SequenceCreate();
void SequenceDelete(struct tagSequence *pDelSequence);

/* sequence */
int SequenceAddTail(struct tagSequence *pSequence, char strLine[]);
FILE *SequenceWriteToFasta(struct tagSequence *pSequence, FILE *fpOut);
int SequenceAddStringToTail(struct tagSequence *pSequence, struct tagString *pString);

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
									char chStrand, int nAliasType);

/* sequence list */
void SequenceListClear(struct tagSequence **pSequenceList);
int SequenceListGetSize(struct tagSequence *pSequenceList);
struct tagSequence *SequenceListGetAt(struct tagSequence *pSequenceList, int nIndex);
struct tagSequence *SequenceListGetTail(struct tagSequence *pSequenceList);
int SequenceListWriteToFastaSeparate(struct tagSequence *pSequenceList);

/* function */
void PrintHelp();
int LoadFullSequenceList(char strInFilePath[], struct tagSequence **pSeqList);

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
			int nMaskType, char strOutFile[]);

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
					  int nUp, int nDown, int nStrandType);

/* ----------------------------------------------------------------------------- */ 
/*                int FastaSequenceSoft2HardMask_Main()                          */
/* This function is used for hard masking a FASTA file.                          */
/* The Maskfile should have the following format:                                */
/* All a,c,g,t in little case will be converted to N.                            */
/* ----------------------------------------------------------------------------- */
int FastaSequenceSoft2HardMask_Main(char strSeqFile[], char strOutFile[]);