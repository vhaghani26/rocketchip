/* ----------------------------------------------------------------------- */
/*  WorkLib.h : interface of the research projects                         */
/*  Author : Ji HongKai ; Time: 2004.10                                    */
/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------- */
/*                              Declaration                                */
/* ----------------------------------------------------------------------- */
#define CMYC_CONS_CUTOFF 18.0
/* #define CMYC_CONS_CUTOFF 18.0 */
#define CMYC_WIN_SIZE 200

struct tagSequence;
struct SEQLINKMAP;

/* ----------------------------------------------------------------------- */ 
/*                                Functions                                */
/* ----------------------------------------------------------------------- */

/* mapping RPL RNA to genome */
int menu_MapRPL(char strAlnPath[], char strRefGenePath[], char strOutPath[], 
				char strSpecies[]);


/* sonic hedgehog insitu */
int menu_shhinsitu();

/* map cMyc binding site to chromosome 21 and 22 */
int Map_cMyc_Chr2122();
int MapcMyc_LOADSEQFROMGENOME(struct SEQMOTIF *pSeqMtf, int nSeqId, 
							 FILE *fpSeq, int nStart, int nEnd, char chStrand);

/* count cMyc in target region */
int Count_cMyc_In_TargetRegion_Main();
int Count_cMyc_In_TargetRegion(char strGenomePath[], int nChr, int nMotifLen, 
							   char strSiteFile[], 
							   char strRegionFile[], char strOutFile[]);


/* map cMyc binding site to chromosome 21 and 22 */
int Map_cMyc_Chr2122_Conserved();
int MapcMyc_LOADSEQFROMGENOME_Conserved(struct SEQMOTIF *pSeqMtf, int nSeqId, 
							 FILE *fpSeq, FILE *fpCS, int nStart, int nEnd, char chStrand);

/* count cMyc in target region */
int Count_cMyc_In_ConservedTargetRegion_Main();
int Count_cMyc_In_ConservedTargetRegion(char strGenomePath[], char strCSPath[],
							   int nChr, int nMotifLen, 
							   char strSiteFile[], char strRegionFile[], 
							   char strOutFile[], double dCutoff);

int Oct4GetPromoter(int argv, char **argc);
int ESClustGetPromoter_Write(struct SEQLINKMAP **pSegList, struct BYTEMATRIX *pCS,
							 int nC, int nMinSegLen, int nMaxGap, int nMinTotLen,
							 struct tagRefGene *pRefGene, FILE *fpOut, int *pId);

/* get genomic conserved segments */
int GenomeGetConservedSeg_Main(int argv, char **argc);
int GenomeGetConservedSeg_Write(struct SEQLINKMAP **pSegList, int nMinSegLen, 
							 int nMaxGap, int nMinTotLen, FILE *fpOut, int *pId);

int HG17GetTSS_Main(int argv, char **argc);

int MM6GetGeneCover_Main(int argv, char **argc);

/* ----------------------------------------------------------------------- */ 
/*  EEL_PrepareOrtholog_Main.                                              */
/* ----------------------------------------------------------------------- */ 
int EEL_PrepareOrtholog_Main(int argv, char **argc);
/* ----------------------------------------------------------------------- */ 
/*  EEL_GetReverseComplement_Main.                                         */
/* ----------------------------------------------------------------------- */ 
int EEL_GetReverseComplement_Main();

int GliArrayGenerateRand();

/* Gene set analysis */
int menu_geneset(int argv, char **argc);

/* RNAi gene map */
int MapHomoloGene();
int MapHomoloGene_ByHID();
int MapGeneID2RefGene();

int Feinberg_CollectGenoType();
int Feinberg_CollectGenoType_mouse();
int Feinberg_CollectGenoType_ASM();
int Feinberg_CollectGenoType_Tumor();

int HTS_selectread(char strInputFile[], char strOutputFile[], double dR);

int IQR_ProbeMatch();

int cMyc_AffyMatch();

int cMyc_AgilentMatch();

int cMyc_ExonMatch();

int book_norm();

int shhmbtest(int argv, char **argc);