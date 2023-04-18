/* ----------------------------------------------------------------------- */
/*  GenomeLib.h : interface of the genome library                          */
/*  Author : Ji HongKai ; Time: 2004.08                                    */
/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------- */
/*                              Declaration                                */
/* ----------------------------------------------------------------------- */
#define GENOME_CONTIG_LEN 500000
#define GENOME_BUFFER_LEN 100000

#define GENOME_MED_CONTIG_LEN 100000
#define GENOME_MED_BUFFER_LEN 20000

#define CODE_BATCH_LEN 60
#define ALN_LENGTH 1000000
#define NAME_LENGTH 20
#define GENE_NAME_LENGTH 50

#define MAX_CHROMOSOME_NUM 1024

/* ORTHOLOG_COLINEAR_NUM can only be an odd number */
#define ORTHOLOG_COLINEAR_NUM 5
#define ORTHOLOG_MAP_COLINEAR_NUM 10
#define ORTHOLOG_COLINEAR_WIN 3000000
#define ORTHOLOG_MATCHCOUNT_LIM 2
#define ORTHOLOG_ORIENCOUNT_LIM 2
#define ORTHOLOG_COMBCOUNT_COEF 0.6

/* REFGENE REDUNDANCY CRITERIA */
#define REFGENE_OVERLAP_TH 0.3
#define REFGENE_MATCH_TH 0.5


#define OS_SYSTEM "UNIX"

/* ----------------------------------------------------------------------- */
/*                              Structures                                 */
/* ----------------------------------------------------------------------- */
/* for organizing RefGene annotation */
struct tagRefGene
{
	/* "Name of Gene" */
	char strGene[GENE_NAME_LENGTH];
	/* "Name of refid" */
    char  strName[NAME_LENGTH];
	/* Gene ID */
	int nGeneID;
	/* "Chromosome name" */
    char  strChrom[NAME_LENGTH];
	/* "Chromosome index" */
	int   nChrom;
	/* "+ or - for strand" */
    char  chStrand;
	/* "Transcription start position" */
    int   nTxStart;     
	/* "Transcription end position" */
    int   nTxEnd;
	/* "Coding region start" */
    int   nCdsStart;      
	/* "Coding region end" */
    int   nCdsEnd;            
	/* "Number of exons" */
    int   nExonCount;        
	/* "Exon start positions": col1: starts, col2: ends. */
    struct INTMATRIX *pmatExonStartsEnds;   
	/* next */
	struct tagRefGene *pNext;
};

/* for organizing genomic region */
struct tagGenomicRegion
{
	/* "Chromosome name" */
    char  strChrom[NAME_LENGTH];
	/* "Chromosome index" */
	int   nChrom;
	/* "+ or - for strand" */
    char  chStrand;
	/* "start position" */
    int   nStart;     
	/* "end position" */
    int   nEnd;
	/* next */
	struct tagGenomicRegion *pNext;
};

/* ----------------------------------------------------------------------- */ 
/*                                Functions                                */
/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------- */ 
/*  RemoveFiles:                                                           */
/*  remove files.                                                          */
/* ----------------------------------------------------------------------- */ 
void RemoveFiles(char strPath[]);

/* ----------------------------------------------------------------------- */ 
/*  AdjustDirectoryPath:                                                   */
/*  add appropriate '/' or '\' to directory path based on OS system.       */
/* ----------------------------------------------------------------------- */ 
void AdjustDirectoryPath(char strPath[]);

/* ----------------------------------------------------------------------- */ 
/*  Retrive file name:                                                     */
/*  Get file name from a file path.                                        */
/* ----------------------------------------------------------------------- */ 
int GetFileName(char strFilePath[], char strFileName[]);

/* ----------------------------------------------------------------------- */ 
/*  File_Bed2Cod_Main:                                                     */
/*  Convert a BED file to a COD file.                                      */
/* ----------------------------------------------------------------------- */ 
int File_Bed2Cod_Main(char strInputPath[], char strOutputPath[]);

/* ----------------------------------------------------------------------- */ 
/*  File_Cod2Bed_Main:                                                     */
/*  Convert a COD file to a BED file.                                      */
/* ----------------------------------------------------------------------- */ 
int File_Cod2Bed_Main(char strInputPath[], char strOutputPath[]);

/* ----------------------------------------------------------------------- */ 
/*  File_Bed2Aln_Main:                                                     */
/*  Convert a BED file to a ALN file.                                      */
/* ----------------------------------------------------------------------- */ 
int File_Bed2Aln_Main(char strInputPath[], char strOutputPath[]);

/* ----------------------------------------------------------------------- */ 
/*  Genome_Fasta_To_Code_4bit_Main:                                        */
/*  Convert Genome Sequence from Fasta files to coding files.              */
/*  1 byte for 2 bases.                                                    */
/*  0000: A; 0001: C; 0010: G; 0011: T.                                    */
/*  0100: a; 0101: c; 0110: g; 0111: t.                                    */
/*  1000: N;                                                               */ 
/*  return number of fasta files coded.                                    */
/* ----------------------------------------------------------------------- */ 
int Genome_Fasta_To_Code_4bit_Main(char strInPath[], char strOutPath[]);

/* ----------------------------------------------------------------------- */ 
/*  Genome_Fasta_To_Code_4bit:                                             */
/*  Convert Genome Sequence from Fasta files to coding files.              */
/*  1 byte for 2 bases.                                                    */
/*  0000: A; 0001: C; 0010: G; 0011: T.                                    */
/*  0100: a; 0101: c; 0110: g; 0111: t.                                    */
/*  1000: N;                                                               */ 
/* ----------------------------------------------------------------------- */ 
int Genome_Fasta_To_Code_4bit(char strInFile[], char strOutFile[]);

/* ----------------------------------------------------------------------- */ 
/*  Genome_PhastCons_To_Code_8bit_Main:                                    */
/*  Convert PhastCons score to genomelab coded score.                      */
/*  1 byte per base.                                                       */
/*  255 highest score.                                                     */
/*  return length of chromosome coded.                                     */
/* ----------------------------------------------------------------------- */ 
int Genome_PhastCons_To_Code_8bit_Main(char strGenomePath[], char strPhastPath[], 
									   char strOutPath[], char strExt[]);

/* ----------------------------------------------------------------------- */ 
/*  Genome_PhastCons_To_Code_8bit:                                         */
/*  Convert PhastCons score to genomelab coded score.                      */
/*  1 byte per base.                                                       */
/*  255 highest score.                                                     */
/*  return length of chromosome coded.                                     */
/* ----------------------------------------------------------------------- */ 
int Genome_PhastCons_To_Code_8bit(char strPhastFile[], char strOutFile[], int nChrLen);

/* ----------------------------------------------------------------------- */ 
/*  Genome_PhastCons_To_Code_8bit_Main_v2:                                 */
/*  To deal with phastcons Format released Nov. 2004 and later.            */ 
/*  Convert PhastCons score to genomelab coded score.                      */
/*  1 byte per base.                                                       */
/*  255 highest score.                                                     */
/*  return length of chromosome coded.                                     */
/* ----------------------------------------------------------------------- */ 
int Genome_PhastCons_To_Code_8bit_Main_v2(char strGenomePath[], char strPhastPath[], 
									   char strOutPath[], char strExt[]);

/* ----------------------------------------------------------------------- */ 
/*  Genome_PhastCons_To_Code_8bit_v2:                                      */
/*  To deal with phastcons Format released Nov. 2004 and later.            */ 
/*  Convert PhastCons score to genomelab coded score.                      */
/*  1 byte per base.                                                       */
/*  255 highest score.                                                     */
/*  return length of chromosome coded.                                     */
/* ----------------------------------------------------------------------- */ 
int Genome_PhastCons_To_Code_8bit_v2(char strPhastFile[], char strOutFile[], int nChrLen);

/* ----------------------------------------------------------------------- */ 
/*  Genome_CodeCDS_FromRefGene_Main:                                       */
/*  Code CDS regions.                                                      */
/*  1 byte for 1 base.                                                     */
/*  0: intergenic.                                                         */
/*  1: CDS                                                                 */
/* ----------------------------------------------------------------------- */ 
int Genome_CodeCDS_FromRefGene_Main(char strGenomePath[], char strRefGenePath[], 
					int nGType, char strSpecies[], char strOutputPath[]);

/* ----------------------------------------------------------------------- */ 
/*  Genome_Code_4bit_GetSeq_Main:                                          */
/*  Get sequences from Genome Sequence coding files.                       */
/*  1 byte for 2 bases.                                                    */
/*  0000: A; 0001: C; 0010: G; 0011: T.                                    */
/*  0100: a; 0101: c; 0110: g; 0111: t.                                    */
/*  1000: N;                                                               */
/*  nStrandType = 0: get + strand in relative to the .sq file.             */
/*                1: get strand according to Infile, if '+', + strand of   */
/*                   .sq file; if'-', - strand of .sq file.                */ 
/* ----------------------------------------------------------------------- */ 
int Genome_Code_4bit_GetSeq_Main(char strGenomePath[], char strSpecies[],
								 char strInFile[], char strOutFile[], 
								 int nStrandType);

/* ----------------------------------------------------------------------- */ 
/*  Genome_Code_4bit_GetSeq_C_Main:                                        */
/*  Get sequences from Genome Sequence coding files.                       */
/*  1 byte for 2 bases.                                                    */
/*  0000: A; 0001: C; 0010: G; 0011: T.                                    */
/*  0100: a; 0101: c; 0110: g; 0111: t.                                    */
/*  1000: N;                                                               */
/*  nStrandType = 0: get + strand in relative to the .sq file.             */
/*                1: get strand according to Infile, if '+', + strand of   */
/*                   .sq file; if'-', - strand of .sq file.                */ 
/* ----------------------------------------------------------------------- */ 
int Genome_Code_4bit_GetSeq_C_Main(char strGenomePath[], char strSpecies[],
								 char strInFile[], char strOutFile[], 
								 int nStrandType);

/* ----------------------------------------------------------------------- */ 
/*  Genome_Code_4bit_GetSeq:                                               */
/*  Get sequences from Genome Sequence coding files.                       */
/*  1 byte for 2 bases.                                                    */
/*  0000: A; 0001: C; 0010: G; 0011: T.                                    */
/*  0100: a; 0101: c; 0110: g; 0111: t.                                    */
/*  1000: N;                                                               */ 
/* ----------------------------------------------------------------------- */ 
struct tagSequence *Genome_Code_4bit_GetSeq(char strInFile[], int nStart, int nEnd);

/* ----------------------------------------------------------------------- */ 
/*  Genome_Code_4bit_GetUncompressed:                                      */
/*  Get uncompressed code from Genome Sequence coding files.               */
/*  1 byte for 2 bases.                                                    */
/*  0000: A; 0001: C; 0010: G; 0011: T.                                    */
/*  0100: a; 0101: c; 0110: g; 0111: t.                                    */
/*  1000: N;                                                               */ 
/*  Return number of bases read.                                           */
/* ----------------------------------------------------------------------- */ 
int Genome_Code_4bit_GetUncompressed(struct BYTEMATRIX *pSeq, char strInFile[], int nStart, int nEnd);

/* ----------------------------------------------------------------------- */ 
/*  Genome_Code_4bit_GetSeqCS_Main:                                        */
/*  Get sequences and conservation score from Genome Sequence coding files */
/*  1 byte for 2 bases.                                                    */
/*  0000: A; 0001: C; 0010: G; 0011: T.                                    */
/*  0100: a; 0101: c; 0110: g; 0111: t.                                    */
/*  1000: N;                                                               */
/*  nStrandType = 0: get + strand in relative to the .sq file.             */
/*                1: get strand according to Infile, if '+', + strand of   */
/*                   .sq file; if'-', - strand of .sq file.                */ 
/* ----------------------------------------------------------------------- */ 
int Genome_Code_4bit_GetSeqCS_Main(char strGenomePath[], char strConsPath[],
								 char strSpecies[],
								 char strInFile[], char strOutPath[],
								 char strOutTitle[],
								 int nStrandType, char strCSFileType[]);

/* ----------------------------------------------------------------------- */ 
/*  Genome_Code_4bit_GetSeqCS_C_Main:                                      */
/*  Get sequences and conservation score from Genome Sequence coding files */
/*  1 byte for 2 bases.                                                    */
/*  0000: A; 0001: C; 0010: G; 0011: T.                                    */
/*  0100: a; 0101: c; 0110: g; 0111: t.                                    */
/*  1000: N;                                                               */
/*  nStrandType = 0: get + strand in relative to the .sq file.             */
/*                1: get strand according to Infile, if '+', + strand of   */
/*                   .sq file; if'-', - strand of .sq file.                */ 
/* ----------------------------------------------------------------------- */ 
int Genome_Code_4bit_GetSeqCS_C_Main(char strGenomePath[], char strConsPath[],
								 char strSpecies[],
								 char strInFile[], char strOutPath[],
								 char strOutTitle[],
								 int nStrandType, char strCSFileType[]);

/* ----------------------------------------------------------------------- */ 
/*  Genome_MafAlign_To_FootPrint_Chr:                                      */
/*  Convert maf alignments to footprint score.                             */
/*  return PROC_SUCCESS if successful.                                     */
/* ----------------------------------------------------------------------- */ 
int Genome_MafAlign_To_FootPrint_Chr(char strAlnFile[], 
				char strChrName[], FILE **vfpScore, 
				int nSpeciesNum, struct tagString **vTag, 
				int nRefid, int nChrSize, 
				int nTarid, struct INTMATRIX *pTargetChrSize,
				int nConditionNum, struct INTMATRIX **vCondition,
				int nWindowSize);

/* ----------------------------------------------------------------------- */ 
/*  Genome_Code_4bit_GetMaskedSeq_Main:                                    */
/*  Get sequences from Genome Sequence coding files.                       */
/*  All repeats will be masked with N.                                     */
/*  If nUseCS==1, base pairs with conservation score < dC will be masked   */
/*  with N.                                                                */
/*  If nUseCDS==1, base pairs in coding regions will be masked with N.     */
/*  1 byte for 2 bases.                                                    */
/*  0000: A; 0001: C; 0010: G; 0011: T.                                    */
/*  0100: a; 0101: c; 0110: g; 0111: t.                                    */
/*  1000: N;                                                               */
/*  nStrandType = 0: get + strand in relative to the .sq file.             */
/*                1: get strand according to Infile, if '+', + strand of   */
/*                   .sq file; if'-', - strand of .sq file.                */ 
/* ----------------------------------------------------------------------- */ 
int Genome_Code_4bit_GetMaskedSeq_Main(char strGenomePath[], char strSpecies[],
			int nUseCS, double dC, char strConsPath[], 
			int nUseCDS, char strCDSPath[],
			char strInFile[], char strOutFile[], int nStrandType);

/* ----------------------------------------------------------------------- */ 
/*  Genome_Code_4bit_GetMaskedSeq_C_Main:                                  */
/*  Get sequences from Genome Sequence coding files.                       */
/*  All repeats will be masked with N.                                     */
/*  If nUseCS==1, base pairs with conservation score < dC will be masked   */
/*  with N.                                                                */
/*  If nUseCDS==1, base pairs in coding regions will be masked with N.     */
/*  1 byte for 2 bases.                                                    */
/*  0000: A; 0001: C; 0010: G; 0011: T.                                    */
/*  0100: a; 0101: c; 0110: g; 0111: t.                                    */
/*  1000: N;                                                               */
/*  nStrandType = 0: get + strand in relative to the .sq file.             */
/*                1: get strand according to Infile, if '+', + strand of   */
/*                   .sq file; if'-', - strand of .sq file.                */ 
/* ----------------------------------------------------------------------- */ 
int Genome_Code_4bit_GetMaskedSeq_C_Main(char strGenomePath[], char strSpecies[],
			int nUseCS, double dC, char strConsPath[], 
			int nUseCDS, char strCDSPath[],
			char strInFile[], char strOutFile[], int nStrandType);

/* ----------------------------------------------------------------------- */ 
/*  Genome_Code_4bit_GetMaskedReg_Main:                                    */
/*  Get masked regions.                                                    */
/*  If nUseRepMask == 1, regions will be kept if more than dR*100%         */
/*     base pairs are not repeats.                                         */
/*  If nUseCS==1, base pairs with conservation score >= dC are defined as  */
/*     conserved bp. Regions will be kept if more than dCR*100% bps        */
/*     are conserved.                                                      */
/*  If nUseCDS==1, regions will be kept if more than dCDS*100% bps are not */
/*     in coding regions.                                                  */
/* ----------------------------------------------------------------------- */ 
int Genome_Code_4bit_GetMaskedReg_Main(char strGenomePath[], char strSpecies[],
			int nUseRepMask, double dR, 
			int nUseCS, double dC, double dCR, char strConsPath[], 
			int nUseCDS, double dCDS, char strCdsPath[],
			char strInputFile[], char strOutputFile[]);

/* ----------------------------------------------------------------------- */ 
/*  Genome_Code_4bit_GetMaskedReg_C_Main:                                  */
/*  Get masked regions.                                                    */
/*  If nUseRepMask == 1, regions will be kept if more than dR*100%         */
/*     base pairs are not repeats.                                         */
/*  If nUseCS==1, base pairs with conservation score >= dC are defined as  */
/*     conserved bp. Regions will be kept if more than dCR*100% bps        */
/*     are conserved.                                                      */
/*  If nUseCDS==1, regions will be kept if more than dCDS*100% bps are not */
/*     in coding regions.                                                  */
/* ----------------------------------------------------------------------- */ 
int Genome_Code_4bit_GetMaskedReg_C_Main(char strGenomePath[], char strSpecies[],
			int nUseRepMask, double dR, 
			int nUseCS, double dC, double dCR, char strConsPath[], 
			int nUseCDS, double dCDS, char strCdsPath[],
			char strInputFile[], char strOutputFile[]);

/* ----------------------------------------------------------------------- */ 
/*  Genome_GetConsScore:                                                   */
/*  Get conservation scores from Genome Sequence coding files.             */
/*  1 byte for 2 bases.                                                    */
/*  0000: A; 0001: C; 0010: G; 0011: T.                                    */
/*  0100: a; 0101: c; 0110: g; 0111: t.                                    */
/*  1000: N;                                                               */ 
/* ----------------------------------------------------------------------- */ 
struct BYTEMATRIX *Genome_GetConsScore(char strConsFile[], int nStart, int nEnd);

/* ----------------------------------------------------------------------- */ 
/*  Genome_GetCDS:                                                         */
/*  Get CDS status from Genome CDS coding files.                           */
/*  1 byte for each base. 0=non cds, 1=cds.                                */
/* ----------------------------------------------------------------------- */ 
struct BYTEMATRIX *Genome_GetCDS(char strCDSFile[], int nStart, int nEnd);

/* ----------------------------------------------------------------------- */ 
/*  ConsScoreWriteToFile_ByStrand:                                         */
/*  write conservation score to files                                      */
/* ----------------------------------------------------------------------- */ 
int ConsScoreWriteToBinaryFile_ByStrand(struct BYTEMATRIX *pCS, char strConsFile[], 
								  char chStrand);

/* ----------------------------------------------------------------------- */ 
/*  ConsScoreWriteToTextFile_ByStrand:                                   */
/*  write conservation score to files                                      */
/* ----------------------------------------------------------------------- */ 
int ConsScoreWriteToTextFile_ByStrand(struct BYTEMATRIX *pCS, char strConsFile[], 
								  char chStrand);

/* ----------------------------------------------------------------------- */ 
/*  ConsScoreWriteToBedFile_ByStrand:                                      */
/*  write conservation score to files                                      */
/* ----------------------------------------------------------------------- */ 
int ConsScoreWriteToBedFile_ByStrand(char strChrName[], int nStart, int nEnd, 
									 struct BYTEMATRIX *pCS, char strConsFile[], 
									 char chStrand);

/* ----------------------------------------------------------------------- */ 
/*  Genome_Index_To_ChromosomeName:                                        */
/*  convert numeric chromosome id to chromsome name                        */
/* ----------------------------------------------------------------------- */ 
char *Genome_Index_To_ChromosomeName(char *strChrName, char strSpecies[], int nChr);

/* ----------------------------------------------------------------------- */ 
/*  Genome_Index_To_ChromosomeName:                                        */
/*  convert chromsome name to numeric chromosome id.                       */
/* ----------------------------------------------------------------------- */ 
int Genome_ChromosomeName_To_Index(char strChrName[], char strSpecies[]);

/* ----------------------------------------------------------------------- */ 
/*  Genome_SpeciesAbbr_To_Species:                                         */
/*  convert species abbreviation name to species name.                     */
/* ----------------------------------------------------------------------- */ 
int Genome_SpeciesAbbr_To_Species(char strAbbr[], char strSpecies[]);

/* ----------------------------------------------------------------------- */ 
/*  Genome_Species_To_SpeciesAbbr:                                         */
/*  convert species name to species abbreviation name.                     */
/* ----------------------------------------------------------------------- */ 
int Genome_Species_To_SpeciesAbbr(char strSpecies[], char strAbbr[]);

/* ----------------------------------------------------------------------- */ 
/*  Genome_MafAlign_To_BG_Main:                                            */
/*  Convert maf alignments to background variation.                        */
/*  return number of pairs of species compared.                            */
/* ----------------------------------------------------------------------- */ 
int Genome_MafAlign_To_BG_Main(char strInfoPath[]);

/* ----------------------------------------------------------------------- */ 
/*  Genome_MafAlign_To_BG_Chr:                                             */
/*  Convert maf alignments to background variation.                        */
/*  return PROC_SUCCESS if successful.                                     */
/* ----------------------------------------------------------------------- */ 
int Genome_MafAlign_To_BG_Chr(char strAlnFile[], char strOutFile[], 
				int nSpeciesNum, struct tagString **vTag, 
				int nRefid, int nChrSize, 
				int nCompNum, struct INTMATRIX *pCompPair,
				int nWindowSize, int nStepSize);

/* ----------------------------------------------------------------------- */ 
/*  Genome_MafAlign_To_CS_Main:                                            */
/*  Convert maf alignments to conservation score.                          */
/*  return number of pairs of species compared.                            */
/* ----------------------------------------------------------------------- */ 
int Genome_MafAlign_To_CS_Main(char strInfoPath[]);

/* ----------------------------------------------------------------------- */ 
/*  Genome_MafAlign_To_CS_Chr:                                             */
/*  Convert maf alignments to conservation score.                          */
/*  return PROC_SUCCESS if successful.                                     */
/* ----------------------------------------------------------------------- */ 
int Genome_MafAlign_To_CS_Chr(char strAlnFile[], char strBGPath[], 
				char strChrName[], FILE **vfpScore, 
				int nSpeciesNum, struct tagString **vTag, 
				int nRefid, int nChrSize, 
				int nTarid, struct INTMATRIX *pTargetChrSize,
				int nCompNum, struct INTMATRIX *pCompPair,
				int nWindowSize, int nStepSize, int nConserveWinSize,
				double dMinPcutoff);

/* ----------------------------------------------------------------------- */ 
/*  Genome_MafAlign_To_FootPrint_Main:                                     */
/*  Convert maf alignments to phylogenetic footprints.                     */
/*  Footprints are specified by conditions such as human+mouse+dog>=2 and  */
/*  window size = 3, etc.                                                  */
/*  Allowed conditions:                                                    */
/*  0: ==; 1: >=; 2: <=; 3: !=; 4: >; 5: <.                                */
/*  return number of alignment files processed.                            */
/* ----------------------------------------------------------------------- */ 
int Genome_MafAlign_To_FootPrint_Main(char strParamPath[]);

/* ----------------------------------------------------------------------- */ 
/*  Genome_CS_Get_Distribution:                                            */
/*  get distribution of conservation score.                                */
/*  return PROC_SUCCESS if successful.                                     */
/* ----------------------------------------------------------------------- */ 
int Genome_CS_Get_Distribution(char strCSPath[], char strFileListName[],
							   char strChrLenPath[], char strOutPath[],
							   char strSpecies[]);

/* ----------------------------------------------------------------------- */ 
/*  Genome_GetGCCS_Summary_Main()                                          */
/*  Get GC content and conservation score distributions for specified      */
/*  regions.                                                               */
/* ----------------------------------------------------------------------- */ 
int Genome_GetGCCS_Summary_Main(char strGenomePath[], char strCodPath[],
					char strOutputPath[], int nUseCS, double dC, 
					char strCSPath[]);

/* ----------------------------------------------------------------------- */ 
/*  Genome_GetRegionCS_Summary_Main()                                      */
/*  Get conservation score summaries for each region specified in the      */
/*  input file.                                                            */
/* ----------------------------------------------------------------------- */ 
int Genome_GetRegionCS_Summary_Main(char strGenomePath[], char strCodPath[],
					char strOutputPath[], int nUseCS, double dC, 
					char strCSPath[]);

/* ----------------------------------------------------------------------- */ 
/*  Genome_GetRegionCS_Summary()                                           */
/*  Get conservation score summaries for a genomic region.                 */
/* ----------------------------------------------------------------------- */ 
int Genome_GetRegionCS_Summary(char strGenomePath[], char strCSPath[], 
					char strChr[], int nStart, int nEnd, double dC, 
					int *pnTotPos, double *pdTotScore, double *pdMeanScore, 
					double *pdSD);

/* ----------------------------------------------------------------------- */ 
/*  Genome_GetMarkovBGForRegion:                                           */
/*  Get markovian background for specific genomic regions.                 */
/* ----------------------------------------------------------------------- */ 
int Genome_GetMarkovBGForRegion(char strGenomePath[], char strGenomicTargetFile[],
								char strSpecies[], int nMCOrder, char strOutFile[]);

/* ----------------------------------------------------------------------- */ 
/*  Genome_MapMotifToRegion()                                              */
/*  Map motif to specific genomic regions.                                 */
/* ----------------------------------------------------------------------- */ 
int Genome_MapMotifToRegion(char strGenomePath[], char strGenomicTargetFile[], 
							char strSpecies[], char strVersion[],
							char strMotifPath[], char strMotifList[], 
							int nBGOrder, char strMarkovBGFile[],
							char strLikehoodRatioCutoffList[],
							char strOutPath[]);

/* ----------------------------------------------------------------------- */ 
/*  Genome_FilterMotifSites()                                              */
/*  Filter motif sites by conservation score and repeatmasker.             */
/*  Sites in initial map (specified by strInitMapPath) will be filtered,   */
/*  and only those with average conservation score > dCSCutoff will be     */
/*  preserved. If nRepeatMask == 1, sites in repeats will be discarded.    */
/*  The remaining sites will be written to strFinalMapPath.                */
/* ----------------------------------------------------------------------- */ 
int Genome_FilterMotifSites(char strGenomePath[], char strCScorePath[], 
					   char strSpecies[], char strVersion[],
					   char strMotifPath[], char strMotifList[], 
					   char strInitMapPath[], char strFinalMapPath[], 
					   double dCSCutoff, int nRepeatMask,
					   int nTypeIIMask, char strTypeIICutoffList[]);

/* ----------------------------------------------------------------------- */ 
/*  RefGeneCreate():                                                       */
/*  Create a new object of tagRefGene.                                     */
/* ----------------------------------------------------------------------- */ 
struct tagRefGene * RefGeneCreate();

/* ----------------------------------------------------------------------- */ 
/*  RefGeneDestroy():                                                      */
/*  Destroy a new object of tagRefGene.                                    */
/* ----------------------------------------------------------------------- */ 
int RefGeneDestroy(struct tagRefGene *pRefGene);

/* ----------------------------------------------------------------------- */ 
/*  RefGeneClone():                                                        */
/*  Clone an object of tagRefGene.                                         */
/* ----------------------------------------------------------------------- */ 
struct tagRefGene *RefGeneClone(struct tagRefGene *pRefGene);

/* ----------------------------------------------------------------------- */ 
/*  RefGeneInit():                                                         */
/*  Initialize reference gene.                                             */
/* ----------------------------------------------------------------------- */ 
int RefGeneInit(struct tagRefGene * pRefGene, char strParamLine[], char strSpecies[]);

/* ----------------------------------------------------------------------- */ 
/*  RefFlatInit():                                                         */
/*  Initialize reference gene.                                             */
/* ----------------------------------------------------------------------- */ 
int RefFlatInit(struct tagRefGene * pRefGene, char strParamLine[], char strSpecies[]);

/* ----------------------------------------------------------------------- */ 
/*  RefGeneInit_FromGenomeLabFormat():                                     */
/*  Initialize reference gene from a genome lab refgene format.            */
/* ----------------------------------------------------------------------- */ 
int RefGeneInit_FromGenomeLabFormat(struct tagRefGene * pRefGene, char strParamLine[], char strSpecies[]);

/* ----------------------------------------------------------------------- */ 
/*  RefFlatInit_FromGenomeLabFormat():                                     */
/*  Initialize reference gene from a genome lab refgene format.            */
/* ----------------------------------------------------------------------- */ 
int RefFlatInit_FromGenomeLabFormat(struct tagRefGene * pRefGene, char strParamLine[], char strSpecies[]);

/* ----------------------------------------------------------------------- */ 
/*  RefLocusInit_FromGenomeLabFormat():                                    */
/*  Initialize reference gene from a genome lab reflocus format.           */
/* ----------------------------------------------------------------------- */ 
int RefLocusInit_FromGenomeLabFormat(struct tagRefGene * pRefGene, char strParamLine[], char strSpecies[]);

/* ----------------------------------------------------------------------- */ 
/*  RefGeneWrite():                                                        */
/*  Write refgene to file.                                                 */
/* ----------------------------------------------------------------------- */ 
int RefGeneWrite(struct tagRefGene * pRefGene, FILE *fpOut);

/* ----------------------------------------------------------------------- */ 
/*  RefFlatWrite():                                                        */
/*  Write refgene to file.                                                 */
/* ----------------------------------------------------------------------- */ 
int RefFlatWrite(struct tagRefGene *pRefGene, FILE *fpOut);

/* ----------------------------------------------------------------------- */ 
/*  RefLocusWrite():                                                       */
/*  Write refgene to file.                                                 */
/* ----------------------------------------------------------------------- */ 
int RefLocusWrite(struct tagRefGene *pRefGene, FILE *fpOut);

/* ----------------------------------------------------------------------- */ 
/*  RefGeneInsert():                                                       */
/*  Insert refgene.                                                        */
/* ----------------------------------------------------------------------- */ 
int RefGeneInsert(struct tagRefGene **vList, struct tagRefGene *pRefGene);

/* ----------------------------------------------------------------------- */ 
/*  RefGene_Convert_UCSC_To_Lab():                                         */
/*  Convert ucsc refgene to lab format.                                    */
/* ----------------------------------------------------------------------- */ 
int RefGene_Convert_UCSC_To_Lab(char strInPath[], char strOutPath[], char strSpecies[], int nChrNum);

/* ----------------------------------------------------------------------- */ 
/*  RefFlat_Convert_UCSC_To_Lab():                                         */
/*  Convert Ucsc refgene to lab format.                                    */
/* ----------------------------------------------------------------------- */ 
int RefFlat_Convert_UCSC_To_Lab(char strInPath[], char strOutPath[], char strSpecies[], int nChrNum);

/* ----------------------------------------------------------------------- */ 
/*  RefGene_AnnotateWithLocusID():                                         */
/*  Annotate refgene with gene ID. The input file is a sorted file that    */
/*  contains the refflat annotations and refgene2geneid annotations.       */
/* ----------------------------------------------------------------------- */ 
int RefGene_AnnotateWithLocusID(char strInPath[], char strOutPath[], char strSpecies[], int nChrNum, int nChangeGeneName);

/* ----------------------------------------------------------------------- */ 
/*  RefGene_AnnotateWithExonArrayID():                                     */
/*  Annotate refgene with exon array ID. The input file is a sorted file   */
/*  that contains the refflat annotations and refgene2exonid annotations.  */
/* ----------------------------------------------------------------------- */ 
int RefGene_AnnotateWithExonArrayID(char strInPath[], char strOutPath[], char strSpecies[], int nChrNum);

/* ----------------------------------------------------------------------- */ 
/*  RefGene_GetNearestGene_Main                                            */
/*  get nearest gene from a refgene database.                              */
/* ----------------------------------------------------------------------- */ 
int RefGene_GetNearestGene_Main(char strDatabasePath[], int nDatabaseType,
			char strSpecies[], char strInputPath[], char strOutputPath[],
			int nRefType, int nUP, int nDOWN);

/* ----------------------------------------------------------------------- */ 
/*  RefGene_GetNearestGene                                                 */
/*  get nearest gene from a refgene database.                              */
/* ----------------------------------------------------------------------- */ 
int RefGene_GetNearestGene(int nChr, int nStart, int nEnd, char chStrand, 
			struct tagRefGene **vSourceRefGene, int nSourceRefNum, 
			int nRefType, int nUP, int nDOWN);

/* ----------------------------------------------------------------------- */ 
/*  RefGene_GetLocationSummary_Main                                        */
/*  get nearest gene from a refgene database.                              */
/* ----------------------------------------------------------------------- */ 
int RefGene_GetLocationSummary_Main(char strDatabasePath[], int nDatabaseType,
			char strSpecies[], char strInputPath[], int nInputType, char strOutputPath[]);

/* ----------------------------------------------------------------------- */ 
/*  RefGene_GetLocationCode                                                */
/*  get location code.                                                     */
/*  exon(5'UTR), exon(3'UTR), exon(other), intron, TSS-up1k, TES-down1K,   */
/*  intergenic(other).                                                     */
/* ----------------------------------------------------------------------- */ 
int RefGene_GetLocationCode(int nChr, int nPos, struct tagRefGene **vSourceRefGene, 
				int nSourceRefNum, struct BYTEMATRIX *pLocInfo);

/* ----------------------------------------------------------------------- */ 
/*  RefGene_GetMatchedControl_Main                                         */
/*  get matched genomic controls.                                          */
/* ----------------------------------------------------------------------- */ 
int RefGene_GetMatchedControl_Main(char strDatabasePath[], int nDatabaseType,
			char strSpecies[], char strChrLenPath[], 
			char strInputPath[], char strOutputPath[],
			int nRepNum, int nRegionLen, int nRemoveRedundancy);

/* ----------------------------------------------------------------------- */ 
/*  RefGene_LoadDatabase:                                                  */
/*  init refgene database from refGene(0)/refFlat(1) file.                 */
/* ----------------------------------------------------------------------- */ 
struct tagRefGene **RefGene_LoadDatabase(char strDatabasePath[], 
					int nDatabaseType, char strSpecies[], int *pSourceRefNum);

/* ----------------------------------------------------------------------- */ 
/*  RefGene_ClearDatabase:                                                 */
/*  clear refgene database.                                                */
/* ----------------------------------------------------------------------- */ 
int RefGene_ClearDatabase(struct tagRefGene ***vSourceRefGene, int nSourceRefNum);

/* ----------------------------------------------------------------------- */ 
/*  RefGene_GetTargetGenomicRegion_TSSTES                                  */
/*  get non-redundant genomic region from refgene list in genomelab        */
/*  refgene format.                                                        */
/* ----------------------------------------------------------------------- */ 
int RefGene_GetTargetGenomicRegion_TSSTES(char strRefGenePath[], 
			char strOutPath[], int nTSSUP, int nTESDOWN, 
			char strSpecies[], int nChrNum, char strChrLenFile[]);

/* ----------------------------------------------------------------------- */ 
/*  RefGene_GetTargetTSSAround():                                          */
/*  get coordinates according to refgene structure.                        */
/* ----------------------------------------------------------------------- */ 
int RefGene_GetTargetTSSAround(char strRefGenePath[], char strTargetPath[],
			int nTSSUP, int nTSSDOWN, char strSpecies[], char strChrLen[],
			char strOutPath[]);

/* ----------------------------------------------------------------------- */ 
/*  RefGene_GetOrtholog():                                                 */
/*  get ortholog genes according to colinear refgene structure & alignmetn */
/* ----------------------------------------------------------------------- */ 
int RefGene_GetOrtholog(char strRefGenePath[], 
			char strSourceSpecies[], char strSourceRefGenePath[],
			char strDestSpecies[], char strMapRefGenePath[],
			char strDestRefGenePath[], 
			char strOutPath[]);

/* ----------------------------------------------------------------------- */ 
/*  RefGene_MergeMultiOrtholog():                                          */
/*  Merge ortholog genes from multiple pairwise ortholog mapping.          */
/* ----------------------------------------------------------------------- */ 
int RefGene_MergeMultiOrtholog(int nSpeciesNum, struct tagString **vSpeciesName,
							   struct tagString **vOrthologPath, char strOutPath[]);

/* ----------------------------------------------------------------------- */ 
/*  RefGene_GetMultiOrtholog1way_Main()                                    */
/*  Get ortholog genes from 1way refgene mapping.                          */
/* ----------------------------------------------------------------------- */ 
int RefGene_GetMultiOrtholog1way_Main(char strTargetPath[], char strParamPath[], 
								  char strOutPath[]);

/* ----------------------------------------------------------------------- */ 
/*  RefGene_GetBestOverlap()                                               */
/*  Get best overlap from a refgene database.                              */
/* ----------------------------------------------------------------------- */ 
int RefGene_GetBestOverlap(struct tagRefGene *pSeedRefGene, 
				struct tagRefGene **vRefGeneDatabase, int nRefGeneNum);

/* ----------------------------------------------------------------------- */ 
/*  RefGene_GetSeededRelaxedOrtholog()                                     */
/*  Get ortholog genes for a seed refgene. No colinearity is required.     */
/* ----------------------------------------------------------------------- */ 
int RefGene_GetSeededRelaxedOrtholog(int nSeedMatchId, 
					struct tagRefGene **vSourceRefGene, int nSourceRefNum, 
					struct tagRefGene **vMapRefGene, int nMapRefNum,
					double *pOptOrienCount, double *pOptMatchCount);

/* ----------------------------------------------------------------------- */ 
/*  RefGene_GetMultiOrtholog1way()                                         */
/*  Get ortholog genes from transcript mapping. This mapping does not      */
/*  require colinearity as getmultiortholog and is therefore much less     */
/*  stringent.                                                             */
/* ----------------------------------------------------------------------- */ 
int RefGene_GetMultiOrtholog1way(char strInPath[], char strOutPath[],
			char strRefSpeciesName[], char strMapSpeciesName[],
			char strRefSpeciesMapbase[], char strMapSpeciesMapbase[],
			char strMapSpeciesDatabase[]);

/* ----------------------------------------------------------------------- */ 
/*  RefGene_GetMultiOrtholog_Main()                                        */
/*  Get ortholog genes from multiple pairwise ortholog mapping.            */
/* ----------------------------------------------------------------------- */ 
int RefGene_GetMultiOrtholog_Main(char strTargetPath[], char strParamPath[], 
								  char strOutPath[]);

/* ----------------------------------------------------------------------- */ 
/*  RefGene_MapBetweenDatabase()                                           */
/*  map a list of refid from one database to another database.             */
/* ----------------------------------------------------------------------- */ 
int RefGene_MapBetweenDatabase(char strInPath[], char strOutPath[],
			char strSpecies[], char strFromDatabase[], char strToDatabase[]);

/* ----------------------------------------------------------------------- */ 
/*  RefGene_RemoveRedundancy()                                             */
/*  Remove redundancy in the target list.                                  */
/* ----------------------------------------------------------------------- */ 
int RefGene_RemoveRedundancy(char strTargetPath[], char strSpecies[], 
							 char strRefGenePath[], char strOutPath[]);

/* ----------------------------------------------------------------------- */ 
/*  RefFlat_RemoveRedundancy()                                             */
/*  Remove redundancy in the target list.                                  */
/* ----------------------------------------------------------------------- */ 
int RefFlat_RemoveRedundancy(char strTargetPath[], char strSpecies[], 
							 char strRefGenePath[], char strOutPath[]);

/* ----------------------------------------------------------------------- */ 
/*  RefGene_Overlap()                                                      */
/*  Judge if two refgene entry are redundant. If yes, return 1; else 0.    */
/* ----------------------------------------------------------------------- */ 
int RefGene_Overlap(struct tagRefGene *pRefGene1, struct tagRefGene *pRefGene2);

/* ----------------------------------------------------------------------- */ 
/*  RefGene_Match()                                                        */
/*  Evaluate how two refgene entries are similar to each other.            */
/* ----------------------------------------------------------------------- */ 
double RefGene_Match(struct tagRefGene *pRefGene1, struct tagRefGene *pRefGene2);

/* ----------------------------------------------------------------------- */ 
/*  RefGene_MapMultiOrthologToOriginalSpecies()                            */
/*  Map refgene orthologs back to the original species.                    */
/* ----------------------------------------------------------------------- */ 
int RefGene_MapMultiOrthologToOriginalSpecies(int nSpeciesNum, struct tagString **vSpeciesName, 
							struct tagString **vSpeciesMapbase, struct tagString **vSpeciesDatabase, 
							int nOrthoGroupNum, char strInPath[], char strOutPath[]);

/* ----------------------------------------------------------------------- */ 
/*  RefGene_MultiOrthologGetTSSAround()                                    */
/*  Get TSS UP and Down sequence for multiortholog.                        */
/* ----------------------------------------------------------------------- */ 
int RefGene_MultiOrthologGetTSSAround(char strInPath[], int nTSSUP, int nTSSDOWN,
									  char strParamPath[], char strOutPath[],
									  char strSeqFile[]);

/* ----------------------------------------------------------------------- */ 
/*  RefGene_MultiOrthologGetTSSAroundExonMasked()                          */
/*  Get TSS UP and Down sequence for multiortholog, exons will be masked   */
/*  in little cases.                                                       */
/* ----------------------------------------------------------------------- */ 
int RefGene_MultiOrthologGetTSSAroundExonMasked(char strInPath[], int nTSSUP, int nTSSDOWN,
									  char strParamPath[], char strOutPath[],
									  char strSeqFile[]);

/* ----------------------------------------------------------------------- */ 
/*  RefGene_PickSpeciesSpecific()                                          */
/*  Pick up specified refgenes from a database.                            */
/*  This can be used to create species-species map using xenoRefGene map   */
/* ----------------------------------------------------------------------- */ 
int RefGene_PickSpeciesSpecific(char strTargetPath[], char strDatabasePath[], 
								char strOutPath[]);

/* ----------------------------------------------------------------------- */ 
/*  RefGene_GetAffy_Main()                                                 */
/*  Link refgene id to affy probeset id.                                   */
/* ----------------------------------------------------------------------- */ 
int RefGene_GetAffy_Main(char strDatabasePath[], char strInputPath[], 
			int nColumn, char strOutputPath[]);

/* ----------------------------------------------------------------------- */ 
/*  RefFlex_GetMultiOrtholog_Main()                                        */
/*  Get ortholog genes from multiple pairwise ortholog mapping.            */
/*  Difference from refgene_getmultiortholog:                              */
/*    (1) refflex_get* uses known coordinates of refgene, whereas          */
/*    refgene_get* only uses refid and has to do the map first.            */
/*    (2) refflex can specify from which column refgene annotation starts  */
/*    whereas refgene cannot.                                              */
/*    (3) reffelx does not ignore anything, nor does it trim redundancy    */
/*    whereas refgene does both.                                           */
/* ----------------------------------------------------------------------- */ 
int RefFlex_GetMultiOrtholog_Main(char strTargetPath[], int nColumn, 
					char strParamPath[], char strOutPath[]);

/* ----------------------------------------------------------------------- */ 
/*  RefFlex_LoadOriginRefGene()                                            */
/*  load refgene from the nColumn-th column of a file.                     */
/* ----------------------------------------------------------------------- */ 
struct tagRefGene **RefFlex_LoadOriginRefGene(char strInPath[], int nColumn, 
						char strSpecies[], int *pRefNum);

/* ----------------------------------------------------------------------- */ 
/*  RefFlex_GetOrtholog():                                                 */
/*  get ortholog genes according to colinear refgene structure & alignmetn */
/* ----------------------------------------------------------------------- */ 
int RefFlex_GetOrtholog(struct tagRefGene **vOriginRefGene, int nOriginRefNum, 
			char strSourceSpecies[], char strSourceRefGenePath[],
			char strDestSpecies[], char strMapRefGenePath[],
			char strDestRefGenePath[], 
			char strOutPath[]);

/* ----------------------------------------------------------------------- */ 
/*  RefFlex_MergeMultiOrtholog():                                          */
/*  Merge ortholog genes from multiple pairwise ortholog mapping.          */
/*  Return the number of ortholog groups.                                  */
/* ----------------------------------------------------------------------- */ 
int RefFlex_MergeMultiOrtholog(int nSpeciesNum, struct tagString **vSpeciesName,
			char strOriginPath[], struct tagRefGene **vOriginRefGene, int nOriginRefNum,
			struct tagString **vOrthologPath, char strOutPath[]);

/* ----------------------------------------------------------------------- */ 
/*  RefFlex_MapMultiOrthologToOriginalSpecies()                            */
/*  Map refgene orthologs back to the original species.                    */
/* ----------------------------------------------------------------------- */ 
int RefFlex_MapMultiOrthologToOriginalSpecies(int nSpeciesNum, struct tagString **vSpeciesName, 
			struct tagString **vSpeciesMapbase, struct tagString **vSpeciesDatabase,
			struct tagRefGene **vOriginRefGene, int nOriginRefNum,
			char strInPath[], char strOutPath[]);

/* ----------------------------------------------------------------------- */ 
/*  RefFlex_GetMultiOrtholog1way_Main()                                    */
/*  Get ortholog genes from 1way refgene mapping.                          */
/* ----------------------------------------------------------------------- */ 
int RefFlex_GetMultiOrtholog1way_Main(char strTargetPath[], char strParamPath[], 
								  char strOutPath[]);

/* ----------------------------------------------------------------------- */ 
/*  RefFlex_GetMultiOrtholog1way()                                         */
/*  Get ortholog genes from transcript mapping. This mapping does not      */
/*  require colinearity as getmultiortholog and is therefore much less     */
/*  stringent.                                                             */
/* ----------------------------------------------------------------------- */ 
int RefFlex_GetMultiOrtholog1way(char strInPath[], char strOutPath[],
			char strRefSpeciesName[], char strMapSpeciesName[],
			char strRefSpeciesMapbase[], char strMapSpeciesMapbase[],
			char strMapSpeciesDatabase[]);

/* ----------------------------------------------------------------------- */ 
/*  GenomicRegionCreate():                                                 */
/*  Create a new object of tagGenomicRegion.                               */
/* ----------------------------------------------------------------------- */ 
struct tagGenomicRegion * GenomicRegionCreate();

/* ----------------------------------------------------------------------- */ 
/*  GenomicRegionDestroy():                                                */
/*  Destroy a new object of tagGenomicRegion.                              */
/* ----------------------------------------------------------------------- */ 
int GenomicRegionDestroy(struct tagGenomicRegion *pRegion);

/* ----------------------------------------------------------------------- */ 
/*  GenomicRegionWrite():                                                  */
/*  Write a genomic region to file.                                        */
/*  If nChromAsNumber == 1, write chromosome as a numeric index.           */
/*  Else write chromosome as a string index.                               */
/* ----------------------------------------------------------------------- */ 
int GenomicRegionWrite(struct tagGenomicRegion *pRegion, FILE *fpOut, int nChromAsNumber);

/* ----------------------------------------------------------------------- */ 
/*  GenomicRegionInsert():                                                 */
/*  Insert a genomic region to a ordered region list.                      */
/* ----------------------------------------------------------------------- */ 
int GenomicRegionInsert(struct tagGenomicRegion **vList, struct tagGenomicRegion *pRegion);

/* ----------------------------------------------------------------------- */ 
/*  GenomicRegionGetUnique():                                              */
/*  get a uniquee genomic region list, and discard any redundant part.     */
/* ----------------------------------------------------------------------- */ 
int GenomicRegionGetUnique(struct tagGenomicRegion **vList);

/* ----------------------------------------------------------------------- */ 
/*  RefGeneSelectRegion():                                                 */
/*  select regions according to refgene structure.                         */
/* ----------------------------------------------------------------------- */ 
struct tagGenomicRegion *RefGeneSelectRegion(struct tagRefGene *pRefGene,
			int nTSSUP, int nTSSDOWN, int nTESUP, int nTESDOWN, 
			int nIncludeIntron, int nIncludeExon, struct INTMATRIX *pChrLen);

/* ----------------------------------------------------------------------- */ 
/*  RefGene_Database_AssignValues_Main()                                   */
/*  Assign significance values for individual genes in a refgene database  */
/* ----------------------------------------------------------------------- */ 
int RefGene_Database_AssignValues_Main(char strDatabasePath[], char strSpecies[], 
			char strOutputPath[], char strValuePath[], int nNormalize,
			double dTruncateLowerBound, double dTruncateUpperBound, char strTransform[], 
			double dPostLowerBound, double dPostUpperBound, char strPostTransform[], 
			int nTakeAbsoluteValue,
			char strMapPath[], char strMapValidationPath[],
			char strNetworkPath[], char strNetworkAnnotationPath[], int nNetDepth);

/* ----------------------------------------------------------------------- */ 
/*  RefGene_GetNeighborGenes_Main                                          */
/*  get neighboring genes from a refgene database.                         */
/* ----------------------------------------------------------------------- */ 
int RefGene_GetNeighborGenes_Main(char strDatabasePath[], int nDatabaseType,
			char strSpecies[], char strInputPath[], char strOutputPath[],
			char strAnnotationPath[], int nUPNum, int nDOWNNum, 
			int nDistanceUpperLimit);

/* ----------------------------------------------------------------------- */ 
/*  RefGene_GetNeighborGenes                                               */
/*  get neighbor genes from a refgene database.                            */
/* ----------------------------------------------------------------------- */ 
int RefGene_GetNeighborGenes(int nChr, int nStart, int nEnd, char chStrand, 
			struct tagRefGene **vSourceRefGene, int nSourceRefNum,
			int nUPNum, int nDOWNNum, int nDistanceUpperLimit,
			int *pUpK, int *pDownK);

/* ----------------------------------------------------------------------- */ 
/*  RefGene_GetLocationInfo                                                */
/*  get relative location info of a region.                                */
/* ----------------------------------------------------------------------- */ 
int RefGene_GetLocationInfo(int nChr, int nStart, int nEnd, char chStrand, 
			struct tagRefGene *pRefGene, int *pTSSDist, int *pTESDist,
			char strLocType[]);

/* ----------------------------------------------------------------------- */ 
/*  RefGene_AnnotateWithMicroArrayID():                                    */
/*  Annotate refgene with exon array ID. The input file is a sorted file   */
/*  that contains the refflat annotations and refgene2exonid annotations.  */
/* ----------------------------------------------------------------------- */ 
int RefGene_AnnotateWithMicroArrayID(char strInPath[], char strOutPath[], char strSpecies[], int nChrNum);

/* ----------------------------------------------------------------------- */ 
/*  RefGene_Database_AssignEdgeValue_Main()                                */
/*  Assign significance values for edges in a protein-protein network.     */
/* ----------------------------------------------------------------------- */ 
int RefGene_Database_AssignEdgeValue_Main(char strDatabasePath[], char strSpecies[], 
			char strOutputPath[], char strValuePath[], int nNormalize, 
			double dTruncateLowerBound, double dTruncateUpperBound, char strTransform[],
			double dPostLowerBound, double dPostUpperBound, char strPostTransform[], 
			int nTakeAbsoluteValue,
			char strMapPath[], char strMapValidationPath[],
			char strNetworkPath[], char strNetworkAnnotationPath[], int nOperationType);

/* ----------------------------------------------------------------------- */ 
/*  Genome_CreateHash_Main                                                 */
/*  Create Hash Table for genome.                                          */
/* ----------------------------------------------------------------------- */ 
int Genome_CreateHash_Main(char strGenomePath[], char strOutPath[], 
						   int nKeyLen);

/* ----------------------------------------------------------------------- */ 
/*  Genome_CreateHash_Chr_v0                                                  */
/*  Create Hash Table for a chromosome.                                    */
/* ----------------------------------------------------------------------- */ 
int Genome_CreateHash_Chr_v0(char strGenomePath[],	char strChr[], char strOutPath[],
				int nChrLen, int nBaseTypeNum, int nKey1Len, int nKey2Len);

/* ----------------------------------------------------------------------- */ 
/*  Genome_CreateHash_Chr                                                  */
/*  Create Hash Table for a chromosome.                                    */
/* ----------------------------------------------------------------------- */ 
int Genome_CreateHash_Chr(char strGenomePath[],	char strChr[], char strOutPath[],
				int nChrLen, int nBaseTypeNum, int nKeyLen);

/* ----------------------------------------------------------------------- */ 
/*  Genome_RegionExtend_Main: extend genomic regions                       */
/* ----------------------------------------------------------------------- */ 
int Genome_RegionExtend_Main(char strInputPath[], char strGenomePath[], 
			char strOutputPath[], char strSpecies[], int nL, int nR, 
			int nCN, int nA, int nUseStrand);
