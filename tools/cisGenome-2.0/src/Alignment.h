/* ----------------------------------------------------------------------- */
/*  Alignment.h : interface of the alignment library                       */
/*  Author : Ji HongKai ; Time: 2005.10                                    */
/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------- */
/*                              Declaration                                */
/* ----------------------------------------------------------------------- */
#define BLAST_SEQNAME_LEN 255
#define SEQLINKMAP_FILLGAP_LEN 20
#define MOTIFSYMBOL_LEN 10
#define ALIGN_LINE_LEN 50
#define ALIGN_NAME_LEN 10
#define ALIGN_START_LEN 10
#define ALIGN_END_LEN 10

/* ----------------------------------------------------------------------- */ 
/*                                 Struct                                  */
/* ----------------------------------------------------------------------- */ 

/* ----------------------------------------------------------------------- */ 
/*              BLASTHIT: for recording blast results                      */
/* ----------------------------------------------------------------------- */ 
struct BLASTHIT
{
	/* query sequence info */
	/* name of query sequence */
	char strQuery[BLAST_SEQNAME_LEN];
	/* start position in query sequence */
	int nQStart;
	/* end position in query sequence */
	int nQEnd;
	/* strand of query */
	char chQStrand;

	/* hit sequence info */
	/* name of hit sequence */
	char strHit[BLAST_SEQNAME_LEN];
	/* start position in hit sequence */
	int nHStart;
	/* end position in hit sequence */
	int nHEnd;
	/* strand of hit */
	char chHStrand;

	/* dScore is the alignment score. */
	double dScore;
	double dScore2;
	/* dValue is the E-Value of alignment. */
	double dEValue;
	/* nAlnLen is the length of the alignment. */
	int nAlnLen;
	/* dIdentities is number of matches/alignment length */
	double dIdentities;

	/* query in alignment form */
	struct tagString *pQAln;
	/* hit in alignment form */
	struct tagString *pHAln;

	/* pNext is the pointer to the next blast hit. */
	struct BLASTHIT *pNext;
};

/* ----------------------------------------------------------------------- */ 
/*  SEQLINKMAP: for tracking sequences in different coordination systems.  */
/* ----------------------------------------------------------------------- */ 
struct SEQLINKMAP
{
	/* dScore is the alignment score. */
	double dScore;
	/* nAlnLen is the length of the original alignment. */
	int nAlnLen;
	/* dIdentities is number of matches/alignment length */
	double dIdentities;

	/* query sequence info */
	/* name of query sequence */
	char strQuery[BLAST_SEQNAME_LEN];
	/* start position in query sequence */
	int nQStart;
	/* end position in query sequence */
	int nQEnd;
	/* strand of query */
	char chQStrand;

	/* hit sequence info */
	/* name of hit sequence */
	char strHit[BLAST_SEQNAME_LEN];
	/* start position in hit sequence */
	int nHStart;
	/* end position in hit sequence */
	int nHEnd;
	/* strand of hit */
	char chHStrand;

	/* genome sequence info */
	/* name of genome sequence */
	char strGenome[BLAST_SEQNAME_LEN];
	/* start position in genome sequence */
	int nGStart;
	/* end position in genome sequence */
	int nGEnd;
	/* strand of genome */
	char chGStrand;

	/* pNext is the pointer to the next seqlinkmap. */
	struct SEQLINKMAP *pNext;
};

/* ----------------------------------------------------------------------- */ 
/*  MOTIFSITELINKMAP: for tracking motif sites in different coordination   */
/*  systems.                                                               */
/* ----------------------------------------------------------------------- */ 
struct MOTIFSITELINKMAP
{
	/* dScore is the alignment score. */
	double dScore;
	/* site sequence */
	char strSiteSeq[LINE_LENGTH];
	/* motif name */
	char strMotifName[LINE_LENGTH];
	/* motif symbol */
	char strMotifSymbol[MOTIFSYMBOL_LEN];
	/* map sequence id */
	int nSeqId;

	/* query sequence info */
	/* name of query sequence */
	char strQuery[BLAST_SEQNAME_LEN];
	/* start position in query sequence */
	int nQStart;
	/* end position in query sequence */
	int nQEnd;
	/* strand of query */
	char chQStrand;

	/* hit sequence info */
	/* name of hit sequence */
	char strHit[BLAST_SEQNAME_LEN];
	/* start position in hit sequence */
	int nHStart;
	/* end position in hit sequence */
	int nHEnd;
	/* strand of hit */
	char chHStrand;

	/* genome sequence info */
	/* name of genome sequence */
	char strGenome[BLAST_SEQNAME_LEN];
	/* start position in genome sequence */
	int nGStart;
	/* end position in genome sequence */
	int nGEnd;
	/* strand of genome */
	char chGStrand;

	/* pNext is the pointer to the next seqlinkmap. */
	struct MOTIFSITELINKMAP *pNext;
};

/* ----------------------------------------------------------------------- */ 
/*                                Functions                                */
/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------- */ 
/*  MAlign_BlastHit_Main: MAlign_BlastHit main function                    */
/*  Using blast to find orthologous segments of a set of sequences.        */
/* ----------------------------------------------------------------------- */ 
int MAlign_BlastHit_Main(char strParamPath[]);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_BlastHit_PrepareIndividualFasta: prepare a fasta file for each  */
/*  single sequence, so that MAlign_BlastHit can search for orthologous    */
/*  segments.                                                              */
/* ----------------------------------------------------------------------- */ 
int MAlign_BlastHit_PrepareIndividualFasta(char strWorkPath[], int nDatabaseNum,
		struct tagString **vDataAlias, struct tagString **vDataFasta, 
		char strOutFile[]);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_BlastHit_PrepareIndividualCod: prepare a genomic coordinates    */
/*  file for each cluster.                                                 */
/* ----------------------------------------------------------------------- */ 
int MAlign_BlastHit_PrepareIndividualCod(char strWorkPath[], int nDatabaseNum,
		struct tagString **vDataAlias, struct tagString **vDataCod, 
		struct tagString **vDataDirec, 
		char strOutFile[], int nClusterNum);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_BlastHit_PrepareOrthologSegments: construct ortholog segments   */
/*  for each query sequence.                                               */
/* ----------------------------------------------------------------------- */ 
int MAlign_BlastHit_PrepareOrthologSegments(char strWorkPath[], int nDatabaseNum,
			struct tagString **vDataAlias, 
			char strOutFile[], int nClusterId,
			char strBlastPath[], double dEValueCut, int nSmallMask,
			int nMinLenCut, double dMinIdnCut,
			int nColinear, double dExtendPercent, int nMaxExtension,
			int nLinkGap, int nKeepBest);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_BlastHit_GetHitSegments: get a set of consistent blast results  */
/* ----------------------------------------------------------------------- */ 
struct BLASTHIT *MAlign_BlastHit_GetHitSegments(char strBlsOutFile[], 
			struct tagSequence *pSeedSeq, struct tagSequence *pDataSeq,
			int nMinLenCut, double dMinIdnCut, int nColinear);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_BlastHit_CheckColinearity: check colinearity with a list of hit */
/* ----------------------------------------------------------------------- */ 
int MAlign_BlastHit_CheckColinearity(struct BLASTHIT *pBlastHit, struct BLASTHIT *pHitList);


/* ----------------------------------------------------------------------- */ 
/*  MAlign_BlastHit_TrimOverlap: trim segments of alignment that are       */
/*  already covered by pBlastHit.                                          */
/* ----------------------------------------------------------------------- */ 
int MAlign_BlastHit_TrimOverlap(struct BLASTHIT *pBlastHit, struct BLASTHIT **pInitList,
								int nMinLenCut, double dMinIdnCut);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_BlastHit_LinkHitSegments: link all hits into a sequence.        */
/* ----------------------------------------------------------------------- */ 
struct SEQLINKMAP *MAlign_BlastHit_LinkHitSegments(struct BLASTHIT **pHitList, 
				double dExtendPercent, int nMaxExtension, int nSeqLen);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_BlastHit_SortLinkMapByQuery: sort all linkmap elements by their */
/*  locations on query sequence.                                           */
/* ----------------------------------------------------------------------- */ 
int MAlign_BlastHit_SortLinkMapByQuery(struct SEQLINKMAP **pSeqMapList, 
				int nLinkGap, int nKeepBest);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_BlastHit_FillGapsInLinkMap: fill gaps in linkmap.               */
/*  Gaps will be filled using a series of N's, the length is specified by  */
/*  SEQLINKMAP_FILLGAP_LEN.                                                */
/* ----------------------------------------------------------------------- */ 
int MAlign_BlastHit_FillGapsInLinkMap(struct SEQLINKMAP **pSeqMapList);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_BlastHit_ComputeGenomicCod: compute genomic coordinates.        */
/* ----------------------------------------------------------------------- */ 
int MAlign_BlastHit_ComputeGenomicCod(struct SEQLINKMAP *pSeqMapList, 
			char strChr[], int nStart, int nEnd, char chStrand);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_BlastHit_ExportLinkMapSeq: export sequence and link maps        */
/* ----------------------------------------------------------------------- */ 
int MAlign_BlastHit_ExportLinkMapSeq(struct tagSequence *pDataSeq, char strSeqAlias[],
			struct SEQLINKMAP *pSeqMapList, char strSeqOutFile[], char strCodOutFile[]);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_MLAGAN_Align: create mlagan alignment                           */
/* ----------------------------------------------------------------------- */ 
int MAlign_MLAGAN_Align(char strWorkPath[], int nDatabaseNum,
			struct tagString **vDataAlias, 
			char strOutFile[], int nClusterId,
			char strMlaganPath[], char strMlaganParam[]);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_MotifMap_Main: MAlign_MotifMap main function                    */
/*  Map motifs to aligned orthologous seqeunces.                           */
/* ----------------------------------------------------------------------- */ 
int MAlign_MotifMap_Main(char strParamPath[]);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_MotifMap_PrepareFasta: prepare a single fasta file for each     */
/*  species used for MAlign_MotifMap.                                      */
/* ----------------------------------------------------------------------- */ 
int MAlign_MotifMap_PrepareFasta(char strWorkPath[], char strFileHeader[], 
			int nClusterNum, char strSpeciesName[]);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_MotifMap_GroupMotifSiteByClusters: group motif sites according  */
/*  to sequence clusters.                                                  */
/* ----------------------------------------------------------------------- */ 
int MAlign_MotifMap_GroupMotifSiteByClusters(char strWorkPath[], char strFileHeader[], 
			char strSpeciesName[], int nClusterNum, int nMotifNum, 
			struct tagString **vMotifSymbol, struct tagString **vMotifName);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_MotifMap_ComputeGenomicCod: compute genomic coordinates of      */
/*  mapped sites.                                                          */
/* ----------------------------------------------------------------------- */ 
int MAlign_MotifMap_ComputeGenomicCod(char strWorkPath[], char strOriginMap[], 
				char strFileHeader[], int nClusterNum, char strSpeciesName[],
				char strMotifName[]);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_MotifMap_GenerateAlignment: generate alignment including        */
/*  mapped sites.                                                          */
/* ----------------------------------------------------------------------- */ 
int MAlign_MotifMap_GenerateAlignment(char strWorkPath[], char strFileHeader[], 
			int nClusterId, int nSpeciesNum, struct tagString **vSpeciesName,  
			int nMotifNum, struct tagString **vMotifSymbol, 
			struct tagString **vMotifName);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_MotifMap_LoadSeqLinkMap: load sequence link map.                */
/* ----------------------------------------------------------------------- */ 
int MAlign_MotifMap_LoadSeqLinkMap(struct SEQLINKMAP **pLinkMapList, char strFileName[]);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_MotifMap_UpdateAlignAndMotif: update alignment to reflect motif */
/*  sites and repeats.                                                     */
/* ----------------------------------------------------------------------- */ 
int MAlign_MotifMap_UpdateAlignAndMotif(struct tagSequence *pSeq, 
			struct tagString *pMS, char strWorkPath[], char strFileHeader[], 
			int nClusterId, char strSpeciesName[],  
			int nMotifNum, struct tagString **vMotifSymbol, 
			struct tagString **vMotifName);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_MotifMap_ExportAlignment: export alignment.                     */
/* ----------------------------------------------------------------------- */ 
int MAlign_MotifMap_ExportAlignment(int nSpeciesNum, 
		struct tagString **vSpeciesName, struct tagSequence *pAln, 
		struct tagString **vMS, struct SEQLINKMAP **vSeqLinkMap,
		char strOutFile[], int nMotifNum, struct tagString **vMotifSymbol, 
		struct tagString **vMotifName);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_MotifMap_TRUNCSTRING: truncate string to a given length.        */
/* ----------------------------------------------------------------------- */ 
int MAlign_MotifMap_TRUNCSTRING(char strLine[], int nTruncLen);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_MotifMap_GetGenomicCod: get genomic coordinate of a given pos.  */
/* ----------------------------------------------------------------------- */ 
int MAlign_MotifMap_GetGenomicCod(int nPos, struct SEQLINKMAP *pSeqLinkMap);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_GenMotifAln_Main: generate motif alignment from existing   */
/*   motif mapping results.                                                */
/* ----------------------------------------------------------------------- */ 
int MAlign_GenMotifAln_Main(char strParamPath[]);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_GenerateMotifAln_GenerateAlignment: generate alignment          */
/*   including mapped sites.                                               */
/* ----------------------------------------------------------------------- */ 
int MAlign_GenerateMotifAln_GenerateAlignment(char strWorkPath[], char strFileHeader[], 
			int nClusterId, int nSpeciesNum, struct tagString **vSpeciesName);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_GenerateMotifAln_UpdateAlignAndMotif: update alignment to       */
/*   reflect motif sites and repeats.                                      */
/* ----------------------------------------------------------------------- */ 
int MAlign_GenerateMotifAln_UpdateAlignAndMotif(struct tagSequence *pSeq, 
			struct tagString *pMS, char strWorkPath[], char strFileHeader[], 
			int nClusterId, char strSpeciesName[]);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_GenerateMotifAln_ExportAlignment: export alignment.             */
/* ----------------------------------------------------------------------- */ 
int MAlign_GenerateMotifAln_ExportAlignment(int nSpeciesNum, 
		struct tagString **vSpeciesName, struct tagSequence *pAln, 
		struct tagString **vMS, struct SEQLINKMAP **vSeqLinkMap,
		char strOutFile[]);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_Genome_PrepareOrtholog: prepare ortholog region for malign      */
/* ----------------------------------------------------------------------- */
int MAlign_Genome_PrepareOrtholog(char strInPath[], char strOutPath[],
			int nSpeciesNum, int nSkip, int nRefType, int nUp, int nDown);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_Genome_BlastHit_Main: MAlign_BlastHit main function             */
/*  Using blast to find orthologous segments of a set of sequences.        */
/* ----------------------------------------------------------------------- */ 
int MAlign_Genome_BlastHit_Main(char strParamPath[]);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_Genome_PrepareIndividualFasta: Prepare individual fasta file.   */
/* ----------------------------------------------------------------------- */ 
int MAlign_Genome_PrepareIndividualFasta(char strRefLine[], char strWorkPath[], char strOutFile[],
			int nDataId, int nClusterId, struct tagString **vDataAlias, struct tagString **vDataDirec, 
			int nSpeciesNum, struct tagString **vSpeciesName, 
			struct tagString **vSpeciesSeq, struct tagString **vSpeciesCons, 
			struct tagString **vSpeciesAnnot, struct INTMATRIX **vChrLen, 
			FILE *fpCod);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_ModuleMap_Main: MAlign_ModuleMap main function                  */
/*  Search for modules to aligned orthologous seqeunces.                   */
/* ----------------------------------------------------------------------- */ 
int MAlign_ModuleMap_Main(char strParamPath[]);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_ModuleMap_SearchModule:                                         */
/*  Search for modules to aligned orthologous seqeunces.                   */
/* ----------------------------------------------------------------------- */ 
int MAlign_ModuleMap_SearchModule(char strWorkPath[], char strFileHeader[], 
			int nClusterId, int nSpeciesNum, struct tagString **vSpeciesName,  
			int nMotifNum, struct tagString **vMotifSymbol, struct tagString **vMotifName,
			struct INTMATRIX *pMotifC, int nModuleLen, 
			struct INTMATRIX *pModuleC, struct BYTEMATRIX *pModuleM, 
			FILE *fpOut);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_ModuleMap_LoadRawMotifMapping:                                  */
/*  Load raw motif mapping.                                                */
/* ----------------------------------------------------------------------- */ 
struct MOTIFSITELINKMAP *MAlign_ModuleMap_LoadRawMotifMapping(char strInFile[], char strSpecies[]);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_ModuleMap_LoadRawMotifMapping:                                  */
/*  Load raw motif mapping.                                                */
/* ----------------------------------------------------------------------- */ 
struct MOTIFSITELINKMAP *MAlign_ModuleMap_LoadRawMotifMapping(char strInFile[], char strSpecies[]);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_ModuleMap_GetConservedMotif:                                    */
/*  Get conserved motif.                                                   */
/* ----------------------------------------------------------------------- */ 
struct MOTIFSITELINKMAP *MAlign_ModuleMap_GetConservedMotif(int nSpeciesNum, 
			struct tagString **vSpeciesName, struct tagSequence **vAln, 
			int nAlnLen, int nMotifNum, struct tagString **vMotifSymbol, 
			struct tagString **vMotifName, struct INTMATRIX *pMotifC, 
			struct MOTIFSITELINKMAP **vMotifList);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_ModuleMap_ExportModules:                                        */
/*  Find and export modules.                                               */
/* ----------------------------------------------------------------------- */ 
int MAlign_ModuleMap_ExportModules(struct MOTIFSITELINKMAP *pConsMotifList, 
			int nMotifNum, int nModuleLen, struct INTMATRIX *pModuleC, 
			struct BYTEMATRIX *pModuleM, int nClusterId, FILE *fpOut);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_ModuleMap_JudgeModule:                                          */
/*  judge if a segment satisfy module condition.                           */
/* ----------------------------------------------------------------------- */ 
int MAlign_ModuleMap_JudgeModule(int nMotifNum, struct INTMATRIX *pModuleStat, 
			struct INTMATRIX *pModuleC, struct BYTEMATRIX *pModuleM);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_KMerStat_Main: MAlign_KMerStat main function                    */
/*  Count frequency of kmers.                                              */
/* ----------------------------------------------------------------------- */ 
int MAlign_KMerStat_Main(char strParamPath[]);

/* ----------------------------------------------------------------------- */ 
/*  MAlign_KMerStat_Count:                                                 */
/*  Count kmers in aligned orthologous seqeunces.                          */
/* ----------------------------------------------------------------------- */ 
int MAlign_KMerStat_Count(char strWorkPath[], char strFileHeader[], 
			int nClusterId, int nSpeciesNum, struct tagString **vSpeciesName,  
			struct INTMATRIX *pMotifC, int nKmerLen, struct DOUBLEMATRIX *pKmerCount);

/* ----------------------------------------------------------------------- */ 
/*  BLASTHITCREATE: create blast hit structure.                            */
/* ----------------------------------------------------------------------- */ 
struct BLASTHIT *BLASTHITCREATE();

/* ----------------------------------------------------------------------- */ 
/*  BLASTHITDESTROY: destroy an instance of blast hit structure.           */
/* ----------------------------------------------------------------------- */ 
void BLASTHITDESTROY(struct BLASTHIT *pBlastHit);

/* ----------------------------------------------------------------------- */ 
/*  BLASTHITCLEARLIST: clear all elements in a linear list of blast hits.  */
/* ----------------------------------------------------------------------- */ 
void BLASTHITCLEARLIST(struct BLASTHIT **pHitList);

/* ----------------------------------------------------------------------- */ 
/*  BLASTHITLOADFROMBLS: load blast results from a normal blast output.    */
/* ----------------------------------------------------------------------- */ 
struct BLASTHIT *BLASTHITLOADFROMBLS(char strFileName[]);

/* ----------------------------------------------------------------------- */ 
/*  BLASTHIT_COLINEAR: check colinearity of two blast hits.                */
/*  return 1 if they are colinear, otherwise return 0.                     */
/* ----------------------------------------------------------------------- */ 
int BLASTHIT_COLINEAR(struct BLASTHIT *pHit1, struct BLASTHIT *pHit2);

/* ----------------------------------------------------------------------- */ 
/*  BLASTHIT_INSERTBYQUERYLOC: insert a blast hit into hitlist. The list   */
/*  of hits should be sorted according to their location in query sequence */
/* ----------------------------------------------------------------------- */ 
int BLASTHIT_INSERTBYQUERYLOC(struct BLASTHIT *pBlastHit, struct BLASTHIT **pHitList);

/* ----------------------------------------------------------------------- */ 
/*  BLASTHIT_RELATIVELOCBYQUERY: compare relative position of two hits     */
/*  according to their positions in query sequences.                       */
/*  return -1 if pHit1 < pHit2;                                            */
/*          0 if pHit1 == pHit2;                                           */
/*          1 if pHit1 > pHit2                                             */
/* ----------------------------------------------------------------------- */ 
int BLASTHIT_RELATIVELOCBYQUERY(struct BLASTHIT *pHit1, struct BLASTHIT *pHit2);

/* ----------------------------------------------------------------------- */ 
/*  BLASTHIT_TRIMOVERLAP: trim overlap between two blast results. The Gold */
/*  will be kept in its original form, while pBlastHit will be trimed.     */ 
/* ----------------------------------------------------------------------- */ 
int BLASTHIT_TRIMOVERLAP(struct BLASTHIT *pGoldHit, struct BLASTHIT **pBlastHit);

/* ----------------------------------------------------------------------- */ 
/*  BLASTHITUPDATEALNLENANDIDENTITY: update alignment length and identity. */
/*  of a blast hit.                                                        */
/* ----------------------------------------------------------------------- */ 
int BLASTHITUPDATEALNLENANDIDENTITY(struct BLASTHIT *pBlastHit);

/* ----------------------------------------------------------------------- */ 
/*  BLASTHITTRIMEND: trim non-sense aligment (i.e. gap to gap alignment)   */
/*  at both ends of a blast hit.                                           */
/* ----------------------------------------------------------------------- */ 
int BLASTHITTRIMEND(struct BLASTHIT *pBlastHit);

/* ----------------------------------------------------------------------- */ 
/*  BLASTHIT_INSERTBYALNLENANDIDENTITY: insert a blast hit into hitlist.   */
/*  The list of hits should be sorted according to their length and        */
/*  identity.                                                              */
/* ----------------------------------------------------------------------- */ 
int BLASTHIT_INSERTBYALNLENANDIDENTITY(struct BLASTHIT *pBlastHit, struct BLASTHIT **pHitList);

/* ----------------------------------------------------------------------- */ 
/*  BLASTHIT_RELATIVELOCBYALNLENANDIDENTITY: compare relative position of  */
/*  two hits according to their length and identity.                       */
/*  return -1 if similarity level of pHit1 > pHit2;                        */
/*          0 if pHit1 == pHit2;                                           */
/*          1 if pHit1 < pHit2                                             */
/* ----------------------------------------------------------------------- */ 
int BLASTHIT_RELATIVELOCBYALNLENANDIDENTITY(struct BLASTHIT *pHit1, struct BLASTHIT *pHit2);

/* ----------------------------------------------------------------------- */ 
/*  BLASTHIT_INSERTBYHITLOC: insert a blast hit into hitlist. The list     */
/*  of hits should be sorted according to their location in hit sequence   */
/* ----------------------------------------------------------------------- */ 
int BLASTHIT_INSERTBYHITLOC(struct BLASTHIT *pBlastHit, struct BLASTHIT **pHitList);

/* ----------------------------------------------------------------------- */ 
/*  BLASTHIT_RELATIVELOCBYHIT: compare relative position of two hits       */
/*  according to their positions in hit sequences.                         */
/*  return -1 if pHit1 < pHit2;                                            */
/*          0 if pHit1 == pHit2;                                           */
/*          1 if pHit1 > pHit2                                             */
/* ----------------------------------------------------------------------- */ 
int BLASTHIT_RELATIVELOCBYHIT(struct BLASTHIT *pHit1, struct BLASTHIT *pHit2);

/* ----------------------------------------------------------------------- */ 
/*  SEQLINKMAPCREATE: create sequence link structure.                      */
/* ----------------------------------------------------------------------- */ 
struct SEQLINKMAP *SEQLINKMAPCREATE();

/* ----------------------------------------------------------------------- */ 
/*  SEQLINKMAPDESTROY: destroy sequence link structure.                    */
/* ----------------------------------------------------------------------- */ 
void SEQLINKMAPDESTROY(struct SEQLINKMAP *pLinkMap);

/* ----------------------------------------------------------------------- */ 
/*  SEQLINKMAPCLEARLIST: clear all elements in a linear list of linkmaps.  */
/* ----------------------------------------------------------------------- */ 
void SEQLINKMAPCLEARLIST(struct SEQLINKMAP **pLinkMapList);

/* ----------------------------------------------------------------------- */ 
/*  SEQLINKMAP_INSERTBYQUERYLOC: insert a seq link map into a list. The    */
/*  list should be sorted according to their location in query sequence    */
/* ----------------------------------------------------------------------- */ 
int SEQLINKMAP_INSERTBYQUERYLOC(struct SEQLINKMAP *pLinkMap, struct SEQLINKMAP **pLinkMapList);

/* ----------------------------------------------------------------------- */ 
/*  SEQLINKMAP_RELATIVELOCBYQUERY: compare relative position of two links  */
/*  according to their positions in query sequences.                       */
/*  return -1 if pLink1 < pLink2;                                          */
/*          0 if pLink1 == pLink2;                                         */
/*          1 if pLink1 > pLink2                                           */
/* ----------------------------------------------------------------------- */ 
int SEQLINKMAP_RELATIVELOCBYQUERY(struct SEQLINKMAP *pLink1, struct SEQLINKMAP *pLink2);

/* ----------------------------------------------------------------------- */ 
/*  MOTIFSITELINKMAPCREATE: create motif site link structure.              */
/* ----------------------------------------------------------------------- */ 
struct MOTIFSITELINKMAP *MOTIFSITELINKMAPCREATE();

/* ----------------------------------------------------------------------- */ 
/*  MOTIFSITELINKMAPDESTROY: destroy motif site link structure.            */
/* ----------------------------------------------------------------------- */ 
void MOTIFSITELINKMAPDESTROY(struct MOTIFSITELINKMAP *pLinkMap);

/* ----------------------------------------------------------------------- */ 
/*  MOTIFSITELINKMAPCLEARLIST: clear all elements in a linear list of      */
/*  motif site linkmaps.                                                   */
/* ----------------------------------------------------------------------- */ 
void MOTIFSITELINKMAPCLEARLIST(struct MOTIFSITELINKMAP **pLinkMapList);

/* ----------------------------------------------------------------------- */ 
/*  MOTIFSITELINKMAP_INSERTBYQUERYLOC: insert a motif site link map into a */
/*  list. The list should be sorted according to their location in query   */
/*  sequence.                                                              */
/* ----------------------------------------------------------------------- */ 
int MOTIFSITELINKMAP_INSERTBYQUERYLOC(struct MOTIFSITELINKMAP *pLinkMap, 
				struct MOTIFSITELINKMAP **pLinkMapList);

/* ----------------------------------------------------------------------- */ 
/*  MOTIFSITELINKMAP_RELATIVELOCBYQUERY: compare relative position of two  */
/*  motif site links according to their positions in query sequences.      */
/*  return -1 if pLink1 < pLink2;                                          */
/*          0 if pLink1 == pLink2;                                         */
/*          1 if pLink1 > pLink2                                           */
/* ----------------------------------------------------------------------- */ 
int MOTIFSITELINKMAP_RELATIVELOCBYQUERY(struct MOTIFSITELINKMAP *pLink1, 
				struct MOTIFSITELINKMAP *pLink2);

