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

int EEL_PrepareOrtholog_Main(int argv, char **argc);

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
	EEL_PrepareOrtholog_Main(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

/* ----------------------------------------------------------------------- */ 
/*  EEL_PrepareOrtholog_Main.                                              */
/* ----------------------------------------------------------------------- */ 
int EEL_PrepareOrtholog_Main(int argv, char **argc)
{
	/* define */
	char strInFile[MED_LINE_LENGTH];
	char strOutFile[MED_LINE_LENGTH];
	char strGenomeFile[MED_LINE_LENGTH];
	char strWorkPath[MED_LINE_LENGTH];
		
	FILE *fpIn;
	FILE *fpOut;
	int nUp = 0;
	int nDown = 0;

	char strLine[MED_LINE_LENGTH];
	char strLine0[MED_LINE_LENGTH];
	char strRefLine[LONG_LINE_LENGTH];
	char *chSep,*chp;
	int nError = 0;
	int ni,nj,nk;

	/* for genome databases */
	char strConservationType[LINE_LENGTH];
	int nSpeciesNum = 0;
	struct tagString **vSpeciesName = NULL;
	struct tagString **vSpeciesSeq = NULL;
	struct tagString **vSpeciesCons = NULL;
	struct tagString **vSpeciesAnnot = NULL;
	struct INTMATRIX **vChrLen = NULL;

	/* statistics obtained during processing */
	int nSeqNum = 0;
	int nDatabaseNum = 0;
	char strSpecies[LINE_LENGTH];
	struct tagRefGene *pRefGene = NULL;
	int nStart;
	int nEnd;
	int nOK;
	int nChr;
	struct tagSequence *pSeq;
	char strSeqFile[MED_LINE_LENGTH];
		
	/* ---------------------------------- */
	/* Step I: Load Parameters            */
	/* ---------------------------------- */
	if(argv != 7)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	strcpy(strInFile, argc[1]);
	nDatabaseNum = atoi(argc[2]);
	nUp = atoi(argc[3]);
	nDown = atoi(argc[4]);
	strcpy(strGenomeFile, argc[5]);
	strcpy(strWorkPath, argc[6]);
	AdjustDirectoryPath(strWorkPath);
	
	
	/* ---------------------------------- */
	/* Step II: Load Genome Settings      */
	/* ---------------------------------- */
	fpIn = NULL;
	fpIn = fopen(strGenomeFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: cannot load environment settings!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strstr(strLine, "[Strand Type]") == strLine)
		{
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);	
		}
		else if(strstr(strLine, "[Conservation Type]") == strLine)
		{
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);
			sprintf(strConservationType, "%s", chp);
		}
		else if(strstr(strLine, "[Species Number]") == strLine)
		{
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);
			nSpeciesNum = atoi(chp);

			if(nSpeciesNum <= 0)
			{
				printf("Error: MAlign_Genome_BlastHit_Main, no species!\n");
				exit(EXIT_FAILURE);
			}

			vSpeciesName = NULL;
			vSpeciesName = (struct tagString **)calloc(nSpeciesNum, sizeof(struct tagString *));
			if(vSpeciesName == NULL)
			{
				printf("Error: MAlign_Genome_BlastHit_Main, cannot load environment settings!\n");
				exit(EXIT_FAILURE);
			}

			vSpeciesSeq = NULL;
			vSpeciesSeq = (struct tagString **)calloc(nSpeciesNum, sizeof(struct tagString *));
			if(vSpeciesSeq == NULL)
			{
				printf("Error: MAlign_Genome_BlastHit_Main, cannot load environment settings!\n");
				exit(EXIT_FAILURE);
			}

			vSpeciesCons = NULL;
			vSpeciesCons = (struct tagString **)calloc(nSpeciesNum, sizeof(struct tagString *));
			if(vSpeciesCons == NULL)
			{
				printf("Error: MAlign_Genome_BlastHit_Main, cannot load environment settings!\n");
				exit(EXIT_FAILURE);
			}

			vSpeciesAnnot = NULL;
			vSpeciesAnnot = (struct tagString **)calloc(nSpeciesNum, sizeof(struct tagString *));
			if(vSpeciesAnnot == NULL)
			{
				printf("Error: MAlign_Genome_BlastHit_Main, cannot load environment settings!\n");
				exit(EXIT_FAILURE);
			}

			vChrLen = NULL;
			vChrLen = (struct INTMATRIX **)calloc(nSpeciesNum, sizeof(struct INTMATRIX *));
			if(vChrLen == NULL)
			{
				printf("Error: MAlign_Genome_BlastHit_Main, cannot load environment settings!\n");
				exit(EXIT_FAILURE);
			}
		}
		else if(strstr(strLine, "[Species Name]") == strLine)
		{
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);
			StringAddTail(vSpeciesName+ni, chp);
		}
		else if(strstr(strLine, "[Species Genome]") == strLine)
		{
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);
			StringAddTail(vSpeciesSeq+ni, chp);

			sprintf(strLine0, "%schrlen.txt", chp);
			vChrLen[ni] = IMLOAD(strLine0);
			if(vChrLen[ni] == NULL)
			{
				printf("Error: MAlign_Genome_BlastHit_Main, cannot load environment settings!\n");
				exit(EXIT_FAILURE);
			}
		}
		else if(strstr(strLine, "[Species Conservation]") == strLine)
		{
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);
			StringAddTail(vSpeciesCons+ni, chp);
		}
		else if(strstr(strLine, "[Species Annotation]") == strLine)
		{
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);
			StringAddTail(vSpeciesAnnot+ni, chp);
			ni++;
		}
		else
		{
			printf("Error: MAlign_Genome_BlastHit_Main, unknown environment settings!\n");
			exit(EXIT_FAILURE);
		}

	}

	fclose(fpIn);

	if(ni != nSpeciesNum)
	{
		printf("Error: MAlign_Genome_BlastHit_Main, species number not match!\n");
		exit(EXIT_FAILURE);
	}
	
	/* ---------------------------------- */
	/* Step III: Get Ortholog Segments    */
	/* ---------------------------------- */
	fpIn = NULL;
	fpIn = fopen(strInFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: cannot open input file!\n");
	}
	
	nSeqNum = 0;
	while(fgets(strRefLine, LONG_LINE_LENGTH, fpIn) !=NULL )
	{
		StrTrimLeft(strRefLine);
		StrTrimRight(strRefLine);
		if(strRefLine[0] == '\0')
			continue;

		if(strRefLine[0] == '>')
		{
			/* ---------------------------------- */
			/* Step III.1: Prepare Seq&Cod        */
			/* ---------------------------------- */
			fgets(strRefLine, LONG_LINE_LENGTH, fpIn);

			for(ni=0; ni<nDatabaseNum; ni++)
			{
				fgets(strRefLine, LONG_LINE_LENGTH, fpIn);
				StrTrimLeft(strRefLine);
				StrTrimRight(strRefLine);
				
				chSep = strchr(strRefLine, '\t');
				*chSep = '\0';
				chp = chSep+1;
				strcpy(strSpecies, strRefLine);
				
				chSep = strchr(chp, '\t');
				*chSep = '\0';
				chp = chSep+1;

				chSep = strchr(chp, '\t');
				*chSep = '\0';
				chp = chSep+1;

				if(strstr(chp, "---") == chp)
				{
					sprintf(strOutFile, "%s%d_%s.fa", strWorkPath, nSeqNum, strSpecies);
					fpOut = NULL;
					fpOut = fopen(strOutFile, "w");
					if(fpOut == NULL)
					{
						printf("Error: cannot load output file!\n");
						exit(EXIT_FAILURE);
					}
					fclose(fpOut);
					continue;
				}

				pRefGene = NULL;
				pRefGene = RefGeneCreate();
				RefGeneInit_FromGenomeLabFormat(pRefGene, chp, strSpecies);
				if(pRefGene->chStrand == '-')
				{
					nStart = pRefGene->nTxStart-nDown;
					nEnd = pRefGene->nTxEnd+nUp;
				}
				else
				{
					nStart = pRefGene->nTxStart-nUp;
					nEnd = pRefGene->nTxEnd+nDown;
				}
				if(nStart < 0)
					nStart = 0;

				/* adjust end */
				nOK = 0;
				for(nj=0; nj<nSpeciesNum; nj++)
				{
					if(strcmp(strSpecies, vSpeciesName[nj]->m_pString) == 0)
					{
						nOK = 1;
						break;
					}
				}
				if(nOK == 0)
				{
					printf("Error: cannot find database information for %s\n", strSpecies);
					exit(EXIT_FAILURE);
				}
				else
				{
					nChr = pRefGene->nChrom-1;
					if(nEnd >= vChrLen[nj]->pMatElement[nChr])
						nEnd = vChrLen[nj]->pMatElement[nChr]-1;
				}

				/* load seq */
				sprintf(strSeqFile, "%s%s.sq", vSpeciesSeq[nj]->m_pString, pRefGene->strChrom);
				pSeq = NULL;
				pSeq = Genome_Code_4bit_GetSeq(strSeqFile, nStart, nEnd);
				if(pSeq == NULL)
				{
					printf("Error: cannot load sequence!\n");
					exit(EXIT_FAILURE);
				}

				/* mask repeats by N */
				for(nk=0; nk<pSeq->m_nLength; nk++)
				{
					if( (pSeq->m_pSequence->m_pString[nk] == 'a') | (pSeq->m_pSequence->m_pString[nk] == 'c') |
						(pSeq->m_pSequence->m_pString[nk] == 'g') | (pSeq->m_pSequence->m_pString[nk] == 't') )
					{
						pSeq->m_pSequence->m_pString[nk] = 'N';
					}
				}

				/* export */
				sprintf(strOutFile, "%s%d_%s.fa", strWorkPath, nSeqNum, strSpecies);
				fpOut = NULL;
				fpOut = fopen(strOutFile, "w");
				if(fpOut == NULL)
				{
					printf("Error: cannot load output file!\n");
					exit(EXIT_FAILURE);
				}
				sprintf(pSeq->m_strAlias, "%d|%s|%s|%d|%d|%d", nSeqNum, strSpecies, pRefGene->strName,
					pRefGene->nChrom, nStart, nEnd);
				SequenceWriteToFasta_ByStrand(pSeq, fpOut, '+', 1);
				fclose(fpOut);

				SequenceDelete(pSeq);
				RefGeneDestroy(pRefGene);
			}
			
			/* ---------------------------------- */
			/* Step III.2: Blast hit              */
			/* ---------------------------------- */
			
			/* ---------------------------------- */
			/* Step III.3: Remove Temporary file  */
			/* ---------------------------------- */
			nSeqNum++;
		}
	}

	fclose(fpIn);

	/* ---------------------------------- */
	/* Step VI: Release memory            */
	/* ---------------------------------- */
	for(ni=0; ni<nSpeciesNum; ni++)
	{
		DeleteString(vSpeciesName[ni]);
		DeleteString(vSpeciesSeq[ni]);
		DeleteString(vSpeciesCons[ni]);
		DeleteString(vSpeciesAnnot[ni]);
		DestroyIntMatrix(vChrLen[ni]);
	}
	free(vSpeciesName);
	free(vSpeciesSeq);
	free(vSpeciesCons);
	free(vSpeciesAnnot);
	free(vChrLen);

	/* return */
	return PROC_SUCCESS;
}