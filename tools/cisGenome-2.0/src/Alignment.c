/* ----------------------------------------------------------------------- */
/*  Alignment.c : implementation of the alignment library                  */
/*  Author : Ji HongKai ; Time: 2005.10                                    */
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
#include "MotifLib.h"
#include "GenomeLib.h"
#include "AffyLib.h"
#include "Alignment.h"
#include "FlexModuleInterface.h"


/* ----------------------------------------------------------------------- */ 
/*  MAlign_BlastHit_Main: MAlign_BlastHit main function                    */
/*  Using blast to find orthologous segments of a set of sequences.        */
/* ----------------------------------------------------------------------- */ 
int MAlign_BlastHit_Main(char strParamPath[])
{
	/* define */

	/* for loading parameters */
	FILE *fpIn;
	char strLine[MED_LINE_LENGTH];
	char *chSep;
	int nError = 0;
	int ni;

	/* parameters */
	char strBlastPath[MED_LINE_LENGTH];
	double dEValueCut = 1e-5;
	int nSmallMask = 0;
	int nMinLenCut = 30;
	double dMinIdnCut = 0.6;
	int nColinear = 1;
	double dExtendPercent = 0.25;
	int nMaxExtension = 100;
	int nLinkGap = 100;
	int nKeepBest = 1;

	char strMlaganPath[MED_LINE_LENGTH];
	char strMlaganParam[MED_LINE_LENGTH];

	char strWorkPath[MED_LINE_LENGTH];
	char strOutFile[MED_LINE_LENGTH];
	int nDatabaseNum = 0;
	struct tagString **vDataAlias = NULL;
	struct tagString **vDataFasta = NULL;
	struct tagString **vDataCod = NULL;
	struct tagString **vDataDirec = NULL;
	char strDataAlias[MED_LINE_LENGTH];
	char strDataFasta[MED_LINE_LENGTH];
	char strDataCod[MED_LINE_LENGTH];
	char strDataDirec[MED_LINE_LENGTH];

	/* statistics obtained during processing */
	int nSeqNum = 0;
	
	/* ---------------------------------- */
	/* Step I: Load Parameters            */
	/* ---------------------------------- */
	fpIn = NULL;
	fpIn = fopen(strParamPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: MAlign_BlastHit_Main, cannot open the parameter file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strstr(strLine, "[Path of Blast]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strBlastPath, chSep);
			AdjustDirectoryPath(strBlastPath);
		}

		else if(strstr(strLine, "[E-value Cutoff]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				dEValueCut = atof(chSep);
			}
		}

		else if(strstr(strLine, "[Mask Small Case]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nSmallMask = atoi(chSep);
			}
		}

		else if(strstr(strLine, "[Minimal Length]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nMinLenCut = atoi(chSep);
			}
		}

		else if(strstr(strLine, "[Minimal Identity]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				dMinIdnCut = atof(chSep);
			}
		}

		else if(strstr(strLine, "[Colinearity]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nColinear = atoi(chSep);
			}
		}

		else if(strstr(strLine, "[Extension Percentage]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				dExtendPercent = atof(chSep);
			}
		}

		else if(strstr(strLine, "[Maximal Extension]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nMaxExtension = atoi(chSep);
			}
		}

		else if(strstr(strLine, "[Maximal Link Gap]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nLinkGap = atoi(chSep);
			}
		}

		else if(strstr(strLine, "[Keep Best vs. All]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nKeepBest = atoi(chSep);
			}
		}

		else if(strstr(strLine, "[Path of MLAGAN]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strMlaganPath, chSep);
			AdjustDirectoryPath(strMlaganPath);
		}

		else if(strstr(strLine, "[MLAGAN Parameter]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strMlaganParam, chSep);
		}

		else if(strstr(strLine, "[Working Directory]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strWorkPath, chSep);
			AdjustDirectoryPath(strWorkPath);
		}

		else if(strstr(strLine, "[Export File]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			strcpy(strOutFile, chSep);
			
			if(strcmp(strOutFile, "") == 0)
			{
				printf("Error: output file not speciefied!\n");
				nError = 1;
				break;
			}
		}

		else if(strstr(strLine, "[Number of Databases]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nDatabaseNum = atoi(chSep);
			}

			if(nDatabaseNum <= 0)
			{
				printf("Error: you have to speciefy at least 1 database!\n");
				nError = 1;
				break;
			}

			nDatabaseNum += 1;
			
			vDataAlias = (struct tagString **)calloc(nDatabaseNum, sizeof(struct tagString *));
			if(vDataAlias == NULL)
			{
				printf("Error: MAlign_BlastHit_Main, cannot allocate memory for loading sequence alias!\n");
				nError = 1;
				break;
			}

			vDataFasta = (struct tagString **)calloc(nDatabaseNum, sizeof(struct tagString *));
			if(vDataFasta == NULL)
			{
				printf("Error: MAlign_BlastHit_Main, cannot allocate memory for loading Fasta file path!\n");
				nError = 1;
				break;
			}

			vDataCod = (struct tagString **)calloc(nDatabaseNum, sizeof(struct tagString *));
			if(vDataCod == NULL)
			{
				printf("Error: MAlign_BlastHit_Main, cannot allocate memory for loading Cod file path!\n");
				nError = 1;
				break;
			}

			vDataDirec = (struct tagString **)calloc(nDatabaseNum, sizeof(struct tagString *));
			if(vDataDirec == NULL)
			{
				printf("Error: MAlign_BlastHit_Main, cannot allocate memory for loading direction file path!\n");
				nError = 1;
				break;
			}
		}

		else if(strstr(strLine, "[Target & Databases]") == strLine)
		{
			ni = 0;
			while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
			{
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] == '\0')
					continue;

				if(ni >= nDatabaseNum)
				{
					printf("Error: MAlign_BlastHit_Main, database number does not match the database provided!\n");
					nError = 1;
					break;
				}

				sscanf(strLine, "%s %s %s %s", strDataAlias, strDataFasta, strDataCod, strDataDirec);
				StringAddTail(vDataAlias+ni, strDataAlias);
				StringAddTail(vDataFasta+ni, strDataFasta);
				StringAddTail(vDataCod+ni, strDataCod);
				StringAddTail(vDataDirec+ni, strDataDirec);

				ni++;
			}

			if( (ni != nDatabaseNum) || (nError == 1) ) 
			{
				printf("Error: MAlign_BlastHit_Main, database number does not match the database provided!\n");
				nError = 1;
				break;
			}
		}
		
		else if(strLine[0] == '#')
		{
		}

		else
		{
			printf("Error: MAlign_BlastHit_Main, unknown parameters!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	fclose(fpIn);

	if(nError > 0)
	{
		printf("Error: MAlign_BlastHit_Main, cannot load the parameter file correctly!\n");
		exit(EXIT_FAILURE);
	}

	/* ---------------------------------- */
	/* Step II: Prepare Fasta file        */
	/* ---------------------------------- */
	nSeqNum = MAlign_BlastHit_PrepareIndividualFasta(strWorkPath, nDatabaseNum,
		vDataAlias, vDataFasta, strOutFile);

	if(nSeqNum > 0)
	{
		MAlign_BlastHit_PrepareIndividualCod(strWorkPath, nDatabaseNum,
			vDataAlias, vDataCod, vDataDirec, strOutFile, nSeqNum);
	}

	/* ---------------------------------- */
	/* Step III: Blast hit                */
	/* ---------------------------------- */
	for(ni=0; ni<nSeqNum; ni++)
	{
		printf("Preprocessing cluster %d...\n", ni);
		MAlign_BlastHit_PrepareOrthologSegments(strWorkPath, nDatabaseNum,
			vDataAlias, strOutFile, ni,
			strBlastPath, dEValueCut, nSmallMask, nMinLenCut, dMinIdnCut,
			nColinear, dExtendPercent, nMaxExtension, nLinkGap, nKeepBest);
	}

	/* ---------------------------------- */
	/* Step IV: MLAGAN                    */
	/* ---------------------------------- */
	for(ni=0; ni<nSeqNum; ni++)
	{
		printf("Align cluster %d...\n", ni);
		MAlign_MLAGAN_Align(strWorkPath, nDatabaseNum,
			vDataAlias, strOutFile, ni,
			strMlaganPath, strMlaganParam);
	}

	/* ---------------------------------- */
	/* Step V: Remove Files               */
	/* ---------------------------------- */
	sprintf(strLine, "%s%s*.tmp*", strWorkPath, strOutFile);
	RemoveFiles(strLine);
	sprintf(strLine, "%s*.anchors", strWorkPath);
	RemoveFiles(strLine);

	/* ---------------------------------- */
	/* Step VI: Release memory            */
	/* ---------------------------------- */
	for(ni=0; ni<nDatabaseNum; ni++)
	{
		DeleteString(vDataAlias[ni]);
		DeleteString(vDataFasta[ni]);
		DeleteString(vDataCod[ni]);
		DeleteString(vDataDirec[ni]);
	}
	free(vDataAlias);
	free(vDataFasta);
	free(vDataCod);
	free(vDataDirec);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_BlastHit_PrepareIndividualFasta: prepare a fasta file for each  */
/*  single sequence, so that MAlign_BlastHit can search for orthologous    */
/*  segments.                                                              */
/* ----------------------------------------------------------------------- */ 
int MAlign_BlastHit_PrepareIndividualFasta(char strWorkPath[], int nDatabaseNum,
		struct tagString **vDataAlias, struct tagString **vDataFasta, 
		char strOutFile[])
{
	/* define */
	int nSeqNum, nDataSeqNum;
	FILE *fpIn;
	FILE *fpOut;
	char strInFile[MED_LINE_LENGTH];
	char strOutPath[MED_LINE_LENGTH];
	int ni;
	char strLine[LONG_LINE_LENGTH];

	/* divide file one by one */
	for(ni=0; ni<nDatabaseNum; ni++)
	{
		/* open file */
		sprintf(strInFile, "%s%s", strWorkPath, vDataFasta[ni]->m_pString);
		fpIn = NULL;
		fpIn = fopen(strInFile, "r");
		if(fpIn == NULL)
		{
			printf("Error: MAlign_BlastHit_PrepareIndividualFasta, cannot open %s\n", strInFile);
			exit(EXIT_FAILURE);
		}

		/* divide file */
		nDataSeqNum = 0;
		fpOut = NULL;
		while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;

			if(strLine[0] == '>')
			{
				if(fpOut != NULL)
				{
					fclose(fpOut);
					nDataSeqNum++;
				}
				
				sprintf(strOutPath, "%s%s_%s_%d_fa.tmp", strWorkPath, strOutFile, vDataAlias[ni]->m_pString, nDataSeqNum);
				fpOut = NULL;
				fpOut = fopen(strOutPath, "w");
				if(fpOut == NULL)
				{
					printf("Error: MAlign_BlastHit_PrepareIndividualFasta, cannot open %s\n", strOutFile);
					exit(EXIT_FAILURE);
				}
				fprintf(fpOut, "%s\n", strLine);
			}
			else
			{
				StrTrimSpace(strLine);
				fprintf(fpOut, "%s\n", strLine);
			}
		}

		if(fpOut != NULL)
		{
			fclose(fpOut);
			nDataSeqNum++;
		}

		/* close file */
		fclose(fpIn);

		/* check sequence number */
		if(ni == 0)
		{
			nSeqNum = nDataSeqNum;
		}
		else
		{
			if(nDataSeqNum != nSeqNum)
			{
				printf("Error: MAlign_BlastHit_PrepareIndividualFasta, sequence numbers across databases do not match!\n");
				exit(EXIT_FAILURE);
			}
		}
	}
	
	/* return */
	return nSeqNum;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_BlastHit_PrepareIndividualCod: prepare a genomic coordinates    */
/*  file for each cluster.                                                 */
/* ----------------------------------------------------------------------- */ 
int MAlign_BlastHit_PrepareIndividualCod(char strWorkPath[], int nDatabaseNum,
		struct tagString **vDataAlias, struct tagString **vDataCod, 
		struct tagString **vDataDirec, 
		char strOutFile[], int nClusterNum)
{
	/* define */
	FILE **vfpIn;
	FILE *fpOut;
	char strInFile[MED_LINE_LENGTH];
	char strOutPath[MED_LINE_LENGTH];
	int ni,nId;
	char strLine[MED_LINE_LENGTH];
	int nEnd = 0;
	int nNullNum = 0; 
	int nFileEnd;
	struct INTMATRIX *pCodNum;
	struct tagString **vCodString;
	char strSeqName[LINE_LENGTH];
	char strChr[LINE_LENGTH];
	int nGStart,nGEnd;
	char chStrand;

	/* init */
	vfpIn = NULL;
	vfpIn = (FILE **)calloc(nDatabaseNum, sizeof(FILE *));
	if(vfpIn == NULL)
	{
		printf("Error: MAlign_BlastHit_PrepareIndividualCod, cannot open file for exporting coordinates!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nDatabaseNum; ni++)
	{
		vfpIn[ni] = NULL;
		if(strcmp(vDataCod[ni]->m_pString, "NULL") != 0)
		{
			sprintf(strInFile, "%s%s", strWorkPath, vDataCod[ni]->m_pString);
			vfpIn[ni] = fopen(strInFile, "r");
			if(vfpIn[ni] == NULL)
			{
				printf("Warning: MAlign_BlastHit_PrepareIndividualCod, cannot open %s!\n", strInFile);
				nNullNum++;
			}
		}
		else
		{
			nNullNum++;
		}
	}

	vCodString = NULL;
	vCodString = (struct tagString **)calloc(nDatabaseNum, sizeof(struct tagString *));
	if(vCodString == NULL)
	{
		printf("Error: MAlign_BlastHit_PrepareIndividualCod, cannot track cluster coordinates!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nDatabaseNum; ni++)
	{
		vCodString[ni] = CreateString(MED_LINE_LENGTH);
		if(vCodString[ni] == NULL)
		{
			printf("Error: MAlign_BlastHit_PrepareIndividualCod, cannot track cluster coordinates!\n");
			exit(EXIT_FAILURE);
		}
	}


	pCodNum = NULL;
	pCodNum = CreateIntMatrix(1, nDatabaseNum);
	if(pCodNum == NULL)
	{
		printf("Error: MAlign_BlastHit_PrepareIndividualCod, cannot track cluster number!\n");
		exit(EXIT_FAILURE);
	}

	/* write cods files */
	if(nNullNum == nDatabaseNum)
	{
		for(nId=0; nId<nClusterNum; nId++)
		{
			sprintf(strOutPath, "%s%s_%d_cod.tmp", strWorkPath, strOutFile, nId);
			fpOut = NULL;
			fpOut = fopen(strOutPath, "w");
			if(fpOut == NULL)
			{
				printf("Error: MAlign_BlastHit_PrepareIndividualCod, cannot open file for exporting coordinates!\n");
				exit(EXIT_FAILURE);
			}

			for(ni=0; ni<nDatabaseNum; ni++)
			{
				sprintf(strLine, "NA\t-1\t-1\t-1\t?");  
				fprintf(fpOut, "%s\t%s\n", vDataAlias[ni]->m_pString, strLine);
				pCodNum->pMatElement[ni] += 1;
			}

			fclose(fpOut);
		}
	}
	else
	{
		nId = 0;
		while(nEnd == 0)
		{
			for(ni=0; ni<nDatabaseNum; ni++)
			{
				if(vfpIn[ni] != NULL)
				{
					nFileEnd = 1;
					while(fgets(strLine, MED_LINE_LENGTH, vfpIn[ni]) != NULL)
					{
						StrTrimLeft(strLine);
						StrTrimRight(strLine);
						if(strLine[0] != '\0')
						{
							nFileEnd = 0;
							break;
						}
					}

					if(nFileEnd == 1)
					{
						nEnd = 1;
					}
					else
					{
						sscanf(strLine, "%s %s %d %d %c", strSeqName, strChr, &nGStart, &nGEnd, &chStrand);
						if(strcmp(vDataDirec[ni]->m_pString, "genewise") == 0)
						{
							sprintf(vCodString[ni]->m_pString, "%s\t%s\t%s\t%d\t%d\t%c", 
								vDataAlias[ni]->m_pString, strSeqName, strChr, nGStart, nGEnd, chStrand);
						}
						else
						{
							sprintf(vCodString[ni]->m_pString, "%s\t%s\t%s\t%d\t%d\t+", 
								vDataAlias[ni]->m_pString, strSeqName, strChr, nGStart, nGEnd);
						}
					}
				}
				else
				{
					sprintf(vCodString[ni]->m_pString, "%s\tNA\t-1\t-1\t-1\t?", vDataAlias[ni]->m_pString);  
				}
			}

			if(nEnd == 0)
			{
				sprintf(strOutPath, "%s%s_%d_cod.tmp", strWorkPath, strOutFile, nId);
				fpOut = NULL;
				fpOut = fopen(strOutPath, "w");
				if(fpOut == NULL)
				{
					printf("Error: MAlign_BlastHit_PrepareIndividualCod, cannot open file for exporting coordinates!\n");
					exit(EXIT_FAILURE);
				}

				for(ni=0; ni<nDatabaseNum; ni++)
				{
					fprintf(fpOut, "%s\n", vCodString[ni]->m_pString);
					pCodNum->pMatElement[ni] += 1;
				}

				fclose(fpOut);
				nId++;
			}
		}
	}
	
	for(ni=0; ni<nDatabaseNum; ni++)
	{
		if(vfpIn[ni] != NULL)
		{
			fclose(vfpIn[ni]);
		}

		DeleteString(vCodString[ni]);
	}
	free(vfpIn);
	free(vCodString);

	/* check seq */
	for(ni=0; ni<nDatabaseNum; ni++)
	{
		if(pCodNum->pMatElement[ni] != pCodNum->pMatElement[0])
		{
			printf("Error: MAlign_BlastHit_PrepareIndividualCod, cluster number do not match!\n");
			exit(EXIT_FAILURE);
		}
	}
	DestroyIntMatrix(pCodNum);

	if(nId != nClusterNum)
	{
		printf("Error: MAlign_BlastHit_PrepareIndividualCod, cluster number do not match!\n");
		exit(EXIT_FAILURE);
	}

	/* return */
	return PROC_SUCCESS;
}

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
			int nLinkGap, int nKeepBest)
{
	/* define */
	int ni;
	char strSeedFile[MED_LINE_LENGTH];
	char strSeqFile[MED_LINE_LENGTH];
	char strTmpOutFile[MED_LINE_LENGTH];
	char strFileName[MED_LINE_LENGTH];
	char strExportSeqFile[MED_LINE_LENGTH];
	char strExportCodFile[MED_LINE_LENGTH];
	FILE *fpSeqOut = NULL;
	FILE *fpCodOut = NULL;
	struct tagSequence *pSeedSeq = NULL;
	struct tagSequence *pDataSeq = NULL;
	int nP1,nP2,nL1,nL2,nG1,nG2;
	char chP,chL,chG;
	char strCommand[MED_LINE_LENGTH];
	
	/* a set of consistent blast hit */
	struct BLASTHIT *pHitList;
	struct SEQLINKMAP *pSeqMapList;

	/* for recording genomic coordinates */
	struct tagString **vGSeqAlias;
	struct tagString **vGSeqName;
	struct tagString **vGChr;
	struct INTMATRIX *pGStart;
	struct INTMATRIX *pGEnd;
	struct BYTEMATRIX *pGStrand;
	char strGSeqAlias[LINE_LENGTH];
	char strGSeqName[LINE_LENGTH];
	char strGChr[LINE_LENGTH];
	int nGStart,nGEnd;
	char chGStrand;
	FILE *fpIn;
	char strLine[MED_LINE_LENGTH];

	/* init */
	vGSeqAlias = NULL;
	vGSeqAlias = (struct tagString **)calloc(nDatabaseNum, sizeof(struct tagString *));
	if(vGSeqAlias == NULL)
	{
		printf("Error: MAlign_BlastHit_PrepareOrthologSegments, cannot allocate memory for processing genomic coordinates!\n");
		exit(EXIT_FAILURE);
	}

	vGSeqName = NULL;
	vGSeqName = (struct tagString **)calloc(nDatabaseNum, sizeof(struct tagString *));
	if(vGSeqName == NULL)
	{
		printf("Error: MAlign_BlastHit_PrepareOrthologSegments, cannot allocate memory for processing genomic coordinates!\n");
		exit(EXIT_FAILURE);
	}

	vGChr = NULL;
	vGChr = (struct tagString **)calloc(nDatabaseNum, sizeof(struct tagString *));
	if(vGChr == NULL)
	{
		printf("Error: MAlign_BlastHit_PrepareOrthologSegments, cannot allocate memory for processing genomic coordinates!\n");
		exit(EXIT_FAILURE);
	}

	pGStart = NULL;
	pGStart = CreateIntMatrix(1, nDatabaseNum);
	if(pGStart == NULL)
	{
		printf("Error: MAlign_BlastHit_PrepareOrthologSegments, cannot allocate memory for processing genomic coordinates!\n");
		exit(EXIT_FAILURE);
	}

	pGEnd = NULL;
	pGEnd = CreateIntMatrix(1, nDatabaseNum);
	if(pGEnd == NULL)
	{
		printf("Error: MAlign_BlastHit_PrepareOrthologSegments, cannot allocate memory for processing genomic coordinates!\n");
		exit(EXIT_FAILURE);
	}

	pGStrand = NULL;
	pGStrand = CreateByteMatrix(1, nDatabaseNum);
	if(pGStrand == NULL)
	{
		printf("Error: MAlign_BlastHit_PrepareOrthologSegments, cannot allocate memory for processing genomic coordinates!\n");
		exit(EXIT_FAILURE);
	}

	/* load genomic coordinates */
	sprintf(strFileName, "%s%s_%d_cod.tmp", strWorkPath, strOutFile, nClusterId);
	fpIn = NULL;
	fpIn = fopen(strFileName, "r");
	if(fpIn == NULL)
	{
		printf("Error: MAlign_BlastHit_PrepareOrthologSegments, cannot open genomic coordinate file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%s %s %s %d %d %c", strGSeqAlias, strGSeqName, strGChr,
			&nGStart, &nGEnd, &chGStrand);
		StringAddTail(vGSeqAlias+ni, strGSeqAlias);
		StringAddTail(vGSeqName+ni, strGSeqName);
		StringAddTail(vGChr+ni, strGChr);
		pGStart->pMatElement[ni] = nGStart;
		pGEnd->pMatElement[ni] = nGEnd;
		if(chGStrand == '+')
		{
			pGStrand->pMatElement[ni] = 0;
		}
		else if(chGStrand == '-')
		{
			pGStrand->pMatElement[ni] = 1;
		}
		else
		{
			pGStrand->pMatElement[ni] = 2;
		}

		ni++;
	}

	fclose(fpIn);
	if(ni != nDatabaseNum)
	{
		printf("Error: MAlign_BlastHit_PrepareOrthologSegments, number of genomic coordinates != database number!\n");
		exit(EXIT_FAILURE);
	}

	
	/* load the first sequence */
	sprintf(strTmpOutFile, "%s%s_blsout.tmp", strWorkPath, strOutFile);
	sprintf(strSeedFile, "%s%s_%s_%d_fa.tmp", strWorkPath, strOutFile, vDataAlias[0]->m_pString, nClusterId);
	LoadFullSequenceList(strSeedFile, &pSeedSeq);

	/* init file */
	sprintf(strFileName, "%s%s_%d_%s.fa", strWorkPath, strOutFile, nClusterId, vDataAlias[0]->m_pString);
	fpSeqOut = NULL;
	fpSeqOut = fopen(strFileName, "w");
	if(fpSeqOut == NULL)
	{
		printf("Error: MAlign_BlastHit_PrepareOrthologSegments, cannot open sequence output file!\n");
		exit(EXIT_FAILURE);
	}

	sprintf(strFileName, "%s%s_%d_%s.gln", strWorkPath, strOutFile, nClusterId, vDataAlias[0]->m_pString);
	fpCodOut = NULL;
	fpCodOut = fopen(strFileName, "w");
	if(fpCodOut == NULL)
	{
		printf("Error: MAlign_BlastHit_PrepareOrthologSegments, cannot open coordinates output file!\n");
		exit(EXIT_FAILURE);
	}

	/* export first sequence */
	nP1 = 0;
	nP2 = pSeedSeq->m_nLength-1;
	chP = '+';
	nL1 = 0;
	nL2 = nP2;
	chL = '+';
	nG1 = pGStart->pMatElement[0];
	nG2 = pGEnd->pMatElement[0];
	if(pGStrand->pMatElement[0] == 0)
	{
		chG = '+';
	}
	else if(pGStrand->pMatElement[0] == 1)
	{
		chG = '-';
	}
	else
	{
		chG = '?';
	}
	if(strcmp(vDataAlias[0]->m_pString, vGSeqAlias[0]->m_pString) != 0)
	{
		printf("Error: MAlign_BlastHit_PrepareOrthologSegments, the order of genomic coordinates is wrong!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpCodOut, "%s\t%d\t%d\t%c\t%s\t%d\t%d\t%c\t%s\t%d\t%d\t%c\n", 
		vDataAlias[0]->m_pString, nP1, nP2, chP,
		pSeedSeq->m_strAlias, nL1, nL2, chL,
		vGChr[0]->m_pString, nG1, nG2, chG);

	strcpy(pSeedSeq->m_strAlias, vDataAlias[0]->m_pString);
	SequenceWriteToFasta_ByStrand(pSeedSeq, fpSeqOut, '+', 1);
	
	/* close files */
	fclose(fpSeqOut);
	fclose(fpCodOut);

	/* process sequences one by one */
	for(ni=1; ni<nDatabaseNum; ni++)
	{
		/* ------------------------------- */
		/* Step O: load sequence           */
		/* ------------------------------- */
		pDataSeq = NULL;
		sprintf(strSeqFile, "%s%s_%s_%d_fa.tmp", strWorkPath, strOutFile, vDataAlias[ni]->m_pString, nClusterId);
		LoadFullSequenceList(strSeqFile, &pDataSeq);

		/* ------------------------------- */
		/* Step I: blast                   */
		/* ------------------------------- */
		sprintf(strCommand, "%sformatdb -i %s -p F -o T", strBlastPath, strSeqFile);
		system(strCommand);
		if(nSmallMask == 1)
		{
			sprintf(strCommand, "%sblastall -p blastn -d %s -i %s -e %e -o %s -U T", strBlastPath,
				strSeqFile, strSeedFile, dEValueCut, strTmpOutFile);
		}
		else
		{
			sprintf(strCommand, "%sblastall -p blastn -d %s -i %s -e %e -o %s", strBlastPath,
				strSeqFile, strSeedFile, dEValueCut, strTmpOutFile);
		}
		system(strCommand);

		/* ------------------------------- */
		/* Step II: Parse blast result     */
		/* ------------------------------- */
		pHitList = NULL;
		pHitList = MAlign_BlastHit_GetHitSegments(strTmpOutFile, pSeedSeq, pDataSeq,
			nMinLenCut, dMinIdnCut, nColinear);

		/* ------------------------------- */
		/* Step III: Get Link Map          */
		/* ------------------------------- */
		pSeqMapList = NULL;
		pSeqMapList = MAlign_BlastHit_LinkHitSegments(&pHitList, dExtendPercent, 
			nMaxExtension, pDataSeq->m_nLength);
		
		/* clear hitlist */
		BLASTHITCLEARLIST(&pHitList);

		/* sort */
		MAlign_BlastHit_SortLinkMapByQuery(&pSeqMapList, nLinkGap, nKeepBest);
		MAlign_BlastHit_FillGapsInLinkMap(&pSeqMapList);

		if(pGStrand->pMatElement[ni] == 0)
		{
			chG = '+';
		}
		else if(pGStrand->pMatElement[ni] == 1)
		{
			chG = '-';
		}
		else
		{
			chG = '?';
		}
		MAlign_BlastHit_ComputeGenomicCod(pSeqMapList, vGChr[ni]->m_pString, 
			pGStart->pMatElement[ni], pGEnd->pMatElement[ni], chG);

		/* ------------------------------- */
		/* Step IV: Export Seq&LinkMap     */
		/* ------------------------------- */
		sprintf(strExportSeqFile, "%s%s_%d_%s.fa", strWorkPath, strOutFile, nClusterId, vDataAlias[ni]->m_pString);
		sprintf(strExportCodFile, "%s%s_%d_%s.gln", strWorkPath, strOutFile, nClusterId, vDataAlias[ni]->m_pString);
		MAlign_BlastHit_ExportLinkMapSeq(pDataSeq, vDataAlias[ni]->m_pString, 
			pSeqMapList, strExportSeqFile, strExportCodFile);

		/* ------------------------------- */
		/* Step V: Destroy Sequences       */
		/* ------------------------------- */
		SEQLINKMAPCLEARLIST(&pSeqMapList);
		SequenceDelete(pDataSeq);
	}

	/* release memory */
	SequenceDelete(pSeedSeq);
	for(ni=0; ni<nDatabaseNum; ni++)
	{
		DeleteString(vGSeqAlias[ni]);
		DeleteString(vGSeqName[ni]);
		DeleteString(vGChr[ni]);
	}
	free(vGSeqAlias);
	free(vGSeqName);
	free(vGChr);
	DestroyIntMatrix(pGStart);
	DestroyIntMatrix(pGEnd);
	DestroyByteMatrix(pGStrand);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_BlastHit_GetHitSegments: get a set of consistent blast results  */
/* ----------------------------------------------------------------------- */ 
struct BLASTHIT *MAlign_BlastHit_GetHitSegments(char strBlsOutFile[], 
			struct tagSequence *pSeedSeq, struct tagSequence *pDataSeq,
			int nMinLenCut, double dMinIdnCut, int nColinear)
{
	/* define */
	struct BLASTHIT *pHitList = NULL;
	struct BLASTHIT *pInitList = NULL;
	struct BLASTHIT *pBlastHit;
	int nConsistency = 1;

	/* load initial blast results */
	pInitList = BLASTHITLOADFROMBLS(strBlsOutFile);

	/* iteratively add consistent segment to the filtered list */
	while(pInitList != NULL)
	{
		/* remove the first hit from init list */
		pBlastHit = pInitList;
		pInitList = pBlastHit->pNext;
		pBlastHit->pNext = NULL;

		/* check basic length and identify */
		if( (pBlastHit->nAlnLen < nMinLenCut) || (pBlastHit->dIdentities < dMinIdnCut) )
		{
			BLASTHITDESTROY(pBlastHit);
			continue;
		}

		/* check its consistency with current list */
		nConsistency = 1;
		if(nColinear == 1)
		{
			nConsistency = MAlign_BlastHit_CheckColinearity(pBlastHit, pHitList);
		}

		/* if consistent, add it to current list */
		if(nConsistency == 1)
		{
			BLASTHIT_INSERTBYQUERYLOC(pBlastHit, &pHitList);
		
			/* adjust all hits in init hit list and resort the init list according to length and identity */
			MAlign_BlastHit_TrimOverlap(pBlastHit, &pInitList, nMinLenCut, dMinIdnCut);
		}
		/* else delete the hit */
		else
		{
			BLASTHITDESTROY(pBlastHit);
		}
	}

	/* return */
	return pHitList;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_BlastHit_CheckColinearity: check colinearity with a list of hit */
/* ----------------------------------------------------------------------- */ 
int MAlign_BlastHit_CheckColinearity(struct BLASTHIT *pBlastHit, struct BLASTHIT *pHitList)
{
	/* define */
	int nColinear = 1;
	struct BLASTHIT *pHit;

	/* judge */
	if( (pBlastHit == NULL) || (pHitList == NULL) )
	{
		return nColinear;
	}

	pHit = pHitList;
	while(pHit != NULL)
	{
		if( BLASTHIT_COLINEAR(pBlastHit, pHit) != 1)
		{
			nColinear = 0;
			break;
		}
		pHit = pHit->pNext;
	}


	/* return */
	return nColinear;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_BlastHit_TrimOverlap: trim segments of alignment that are       */
/*  already covered by pBlastHit.                                          */
/* ----------------------------------------------------------------------- */ 
int MAlign_BlastHit_TrimOverlap(struct BLASTHIT *pBlastHit, struct BLASTHIT **pInitList,
								int nMinLenCut, double dMinIdnCut)
{
	/* define */
	struct BLASTHIT *pHit,*pHit2;
	struct BLASTHIT *pNewList = NULL;

	/* process one by one */
	if(pInitList == NULL)
	{
		return PROC_SUCCESS;
	}

	if(*pInitList == NULL)
	{
		return PROC_SUCCESS;
	}

	
	/* resort with length & identity */
	while((*pInitList) != NULL)
	{
		pHit = *pInitList;
		*pInitList = pHit->pNext;
		pHit->pNext = NULL;

		/* trim overlap */
		BLASTHIT_TRIMOVERLAP(pBlastHit, &pHit);

		while(pHit != NULL)
		{
			pHit2 = pHit;
			pHit = pHit2->pNext;
			pHit2->pNext = NULL;

			/* insert to new list */
			if( (pHit2->nAlnLen < nMinLenCut) || (pHit2->dIdentities < dMinIdnCut) )
			{
				BLASTHITDESTROY(pHit2);
			}
			else
			{
				BLASTHIT_INSERTBYALNLENANDIDENTITY(pHit2, &pNewList);
			}
		}
	}

	/* update new list */
	*pInitList = pNewList;

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_BlastHit_LinkHitSegments: link all hits into a sequence.        */
/* ----------------------------------------------------------------------- */ 
struct SEQLINKMAP *MAlign_BlastHit_LinkHitSegments(struct BLASTHIT **pHitList, 
				double dExtendPercent, int nMaxExtension, int nSeqLen)
{
	/* define */
	struct SEQLINKMAP *pLinkMapList = NULL;
	struct SEQLINKMAP *pLinkMap,*pPrev;
	struct BLASTHIT *pBlastHit,*pPrevHit;
	struct BLASTHIT *pNewHitList = NULL;
	int nAlnLen,nExtLen;
	int nLen1,nLen2,nDist;
	double dLambda;

	/* sort hitlist by hit sequence */
	if(pHitList == NULL)
	{
		return NULL;
	}
	if(*pHitList == NULL)
	{
		return NULL;
	}

	while(*pHitList != NULL)
	{
		pBlastHit = *pHitList;
		*pHitList = pBlastHit->pNext;
		pBlastHit->pNext = NULL;

		BLASTHIT_INSERTBYHITLOC(pBlastHit, &pNewHitList);
	}
	*pHitList = pNewHitList;

	/* create linkmap */
	pBlastHit = pNewHitList;
	pPrevHit = NULL;
	while(pBlastHit != NULL)
	{
		pLinkMap = NULL;
		pLinkMap = SEQLINKMAPCREATE();
		if(pLinkMap == NULL)
		{
			printf("Error: MAlign_BlastHit_LinkHitSegments, cannot create link map!\n");
			exit(EXIT_FAILURE);
		}

		pLinkMap->dScore = pBlastHit->dEValue;
		pLinkMap->nAlnLen = pBlastHit->nAlnLen;
		pLinkMap->dIdentities = pBlastHit->dIdentities;
		pLinkMap->nQStart = pBlastHit->nQStart;
		pLinkMap->nQEnd = pBlastHit->nQEnd;
		pLinkMap->chQStrand = pBlastHit->chQStrand;
		strcpy(pLinkMap->strQuery, pBlastHit->strQuery);
		pLinkMap->nHStart = pBlastHit->nHStart;
		pLinkMap->nHEnd = pBlastHit->nHEnd;
		strcpy(pLinkMap->strHit, pBlastHit->strHit);
		pLinkMap->chHStrand = pBlastHit->chHStrand;

		if(pLinkMap->chQStrand == '-')
		{
			pLinkMap->chQStrand = '+';
			if(pLinkMap->chHStrand == '-')
				pLinkMap->chHStrand = '+';
			else if(pLinkMap->chHStrand == '+')
				pLinkMap->chHStrand = '-';
			else
				pLinkMap->chHStrand = '?';
		}
		
		nAlnLen = pLinkMap->nHEnd-pLinkMap->nHStart+1;
		nExtLen = (int)(dExtendPercent*nAlnLen);
		if(nExtLen > nMaxExtension)
			nExtLen = nMaxExtension;

		pLinkMap->nHStart -= nExtLen;
		pLinkMap->nHEnd += nExtLen;
		if(pLinkMap->nHStart < 1)
			pLinkMap->nHStart = 1;
		if(pLinkMap->nHEnd > nSeqLen)
			pLinkMap->nHEnd = nSeqLen;

		/* add link map to list */
		if(pLinkMapList == NULL)
		{
			pLinkMapList = pLinkMap;
			pPrev = pLinkMap;
		}
		else
		{
			if( (strcmp(pPrev->strHit, pLinkMap->strHit) == 0) && 
				(pPrev->nHEnd >= pLinkMap->nHStart) )
			{
				nDist = pBlastHit->nHStart-pPrevHit->nHEnd-1;
				if(nDist < 0)
				{
					printf("Error: MAlign_BlastHit_LinkHitSegments, extension conflict!\n");
					exit(EXIT_FAILURE);
				}
				nLen1 = (int)(pPrevHit->nAlnLen*pPrevHit->dIdentities);
				nLen2 = (int)(pBlastHit->nAlnLen*pBlastHit->dIdentities);
				dLambda = (double)nLen1/(double)(nLen1+nLen2);

				nExtLen = (int)(dLambda*nDist);
				if(nExtLen > nMaxExtension)
					nExtLen = nMaxExtension;
				pPrev->nHEnd = pPrevHit->nHEnd+nExtLen;

				nExtLen = (int)((1.0-dLambda)*nDist);
				if(nExtLen > nMaxExtension)
					nExtLen = nMaxExtension;
				pLinkMap->nHStart = pBlastHit->nHStart-nExtLen;
			}

			pPrev->pNext = pLinkMap;
			pPrev = pLinkMap;
		}

		/* get next */
		pPrevHit = pBlastHit;
		pBlastHit = pBlastHit->pNext;
	}

	/* return */
	return pLinkMapList;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_BlastHit_SortLinkMapByQuery: sort all linkmap elements by their */
/*  locations on query sequence.                                           */
/* ----------------------------------------------------------------------- */ 
int MAlign_BlastHit_SortLinkMapByQuery(struct SEQLINKMAP **pSeqMapList, 
				int nLinkGap, int nKeepBest)
{
	/* define */
	struct SEQLINKMAP *pNewMapList = NULL;
	struct SEQLINKMAP *pLinkMap = NULL;
	struct SEQLINKMAP *pPrev = NULL;
	struct SEQLINKMAP *pBest = NULL;

	/* sort by query */
	if(pSeqMapList == NULL)
	{
		return PROC_SUCCESS;
	}
	if(*pSeqMapList == NULL)
	{
		return PROC_SUCCESS;
	}

	/* merge gaps */
	while(*pSeqMapList != NULL)
	{
		pLinkMap = *pSeqMapList;
		*pSeqMapList = pLinkMap->pNext;
		pLinkMap->pNext = NULL;

		SEQLINKMAP_INSERTBYQUERYLOC(pLinkMap, &pNewMapList);
	}

	*pSeqMapList = pNewMapList;

	/* merge gaps */
	pNewMapList = NULL;
	while(*pSeqMapList != NULL)
	{
		pLinkMap = *pSeqMapList;
		*pSeqMapList = pLinkMap->pNext;
		pLinkMap->pNext = NULL;

		if(pNewMapList == NULL)
		{
			pNewMapList = pLinkMap;
			pPrev = pLinkMap;
		}
		else
		{
			if( (strcmp(pPrev->strHit, pLinkMap->strHit) == 0) &&
				(pLinkMap->chHStrand == pPrev->chHStrand) )
			{
				if( (pPrev->chHStrand == '-') && ((pPrev->nHStart-pLinkMap->nHEnd)<nLinkGap) )
				{
					if(pLinkMap->nHStart < pPrev->nHStart)
						pPrev->nHStart = pLinkMap->nHStart;
					if(pLinkMap->dScore < pPrev->dScore)
						pPrev->dScore = pLinkMap->dScore;
					if( (pLinkMap->nAlnLen*pLinkMap->dIdentities) > (pPrev->nAlnLen*pPrev->dIdentities) )
					{
						pPrev->nAlnLen = pLinkMap->nAlnLen;
						pPrev->dIdentities = pLinkMap->dIdentities;
					}

					SEQLINKMAPDESTROY(pLinkMap);
				}
				else if( (pPrev->chHStrand == '+') && ((pLinkMap->nHStart-pPrev->nHEnd)<nLinkGap) )
				{
					if(pLinkMap->nHEnd > pPrev->nHEnd)
						pPrev->nHEnd = pLinkMap->nHEnd;
					if(pLinkMap->dScore < pPrev->dScore)
						pPrev->dScore = pLinkMap->dScore;
					if( (pLinkMap->nAlnLen*pLinkMap->dIdentities) > (pPrev->nAlnLen*pPrev->dIdentities) )
					{
						pPrev->nAlnLen = pLinkMap->nAlnLen;
						pPrev->dIdentities = pLinkMap->dIdentities;
					}

					SEQLINKMAPDESTROY(pLinkMap);
				}
				else
				{
					pPrev->pNext = pLinkMap;
					pPrev = pLinkMap;
				}
			}
			else
			{
				pPrev->pNext = pLinkMap;
				pPrev = pLinkMap;
			}
		}
	}

	*pSeqMapList = pNewMapList;

	/* if keep best == 0, return the best one and clear others */
	if(nKeepBest == 0)
	{
		/* find the best */
		pLinkMap = pNewMapList;
		while(pLinkMap != NULL)
		{
			if(pBest == NULL)
			{
				pBest = pLinkMap;
			}
			else
			{
				if(pLinkMap->dScore < pBest->dScore)
				{
					pBest = pLinkMap;
				}
			}

			pLinkMap = pLinkMap->pNext;
		}

		/* remove everything except for the best one */
		pNewMapList = NULL;
		while(*pSeqMapList != NULL)
		{
			pLinkMap = *pSeqMapList;
			*pSeqMapList = pLinkMap->pNext;
			pLinkMap->pNext = NULL;

			if(pLinkMap != pBest)
			{
				SEQLINKMAPDESTROY(pLinkMap);
			}
		}

		*pSeqMapList = pBest;
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_BlastHit_FillGapsInLinkMap: fill gaps in linkmap.               */
/*  Gaps will be filled using a series of N's, the length is specified by  */
/*  SEQLINKMAP_FILLGAP_LEN.                                                */
/* ----------------------------------------------------------------------- */ 
int MAlign_BlastHit_FillGapsInLinkMap(struct SEQLINKMAP **pSeqMapList)
{
	/* define */
	struct SEQLINKMAP *pGap;
	struct SEQLINKMAP *pLinkMap;
	struct SEQLINKMAP *pPrev;
	int nTotLen,nAlnLen;

	/* init */
	if(pSeqMapList == NULL)
	{
		return PROC_SUCCESS;
	}
	if(*pSeqMapList == NULL)
	{
		pGap = NULL;
		pGap = SEQLINKMAPCREATE();
		strcpy(pGap->strQuery, "GAP");
		strcpy(pGap->strHit, "GAP");
		strcpy(pGap->strGenome, "GAP");
		pGap->nHStart = -1;
		pGap->nHEnd = -1;
		pGap->chHStrand = '?';
		pGap->nQStart = 0;
		pGap->chQStrand = '+';
		pGap->nQEnd = pGap->nQStart+SEQLINKMAP_FILLGAP_LEN-1;
		*pSeqMapList = pGap;

		return PROC_SUCCESS;
	}

	/* fill gaps */
	nTotLen = 0;
	pPrev = *pSeqMapList;
	pPrev->nQStart = nTotLen;
	nAlnLen = pPrev->nHEnd-pPrev->nHStart;
	pPrev->nQEnd = pPrev->nQStart+nAlnLen;
	pPrev->chQStrand = '+';
	pPrev->nHStart -= 1;
	pPrev->nHEnd -= 1;
	nTotLen += nAlnLen+1;
	
	pLinkMap = pPrev->pNext;
	while(pLinkMap != NULL)
	{
		/* fill gap */
		pGap = NULL;
		pGap = SEQLINKMAPCREATE();
		strcpy(pGap->strQuery, "GAP");
		strcpy(pGap->strHit, "GAP");
		strcpy(pGap->strGenome, "GAP");
		pGap->nHStart = -1;
		pGap->nHEnd = -1;
		pGap->chHStrand = '?';
		pGap->nQStart = nTotLen;
		pGap->chQStrand = '+';
		pGap->nQEnd = pGap->nQStart+SEQLINKMAP_FILLGAP_LEN-1;
		nTotLen += SEQLINKMAP_FILLGAP_LEN;
		
		pGap->pNext = pLinkMap;
		pPrev->pNext = pGap;
		pPrev = pGap;

		/* adjust current node */
		pLinkMap->nQStart = nTotLen;
		nAlnLen = pLinkMap->nHEnd-pLinkMap->nHStart;
		pLinkMap->nQEnd = pLinkMap->nQStart+nAlnLen;
		pLinkMap->chQStrand = '+';
		pLinkMap->nHStart -= 1;
		pLinkMap->nHEnd -= 1;
		nTotLen += nAlnLen+1;

		/* get next */
		pPrev = pLinkMap;
		pLinkMap = pPrev->pNext;
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_BlastHit_ComputeGenomicCod: compute genomic coordinates.        */
/* ----------------------------------------------------------------------- */ 
int MAlign_BlastHit_ComputeGenomicCod(struct SEQLINKMAP *pSeqMapList, 
			char strChr[], int nStart, int nEnd, char chStrand)
{
	/* define */
	struct SEQLINKMAP *pLinkMap;

	/* process */
	pLinkMap = pSeqMapList;
	while(pLinkMap != NULL)
	{
		if(strcmp(pLinkMap->strHit, "GAP") == 0)
		{
			strcpy(pLinkMap->strGenome, "GAP");
			pLinkMap->nGStart = -1;
			pLinkMap->nGEnd = -1;
			pLinkMap->chGStrand = '?';
		}
		else
		{
			if(chStrand == '+')
			{
				strcpy(pLinkMap->strGenome, strChr);
				pLinkMap->nGStart = nStart+pLinkMap->nHStart;
				pLinkMap->nGEnd = nStart+pLinkMap->nHEnd;
				if(pLinkMap->chHStrand == '+')
				{
					pLinkMap->chGStrand = '+';
				}
				else if(pLinkMap->chHStrand == '-')
				{
					pLinkMap->chGStrand = '-';
				}
				else
				{
					pLinkMap->chGStrand = '?';
				}
			}
			else if(chStrand == '-')
			{
				strcpy(pLinkMap->strGenome, strChr);
				pLinkMap->nGStart = nEnd-pLinkMap->nHEnd;
				pLinkMap->nGEnd = nEnd-pLinkMap->nHStart;
				if(pLinkMap->chHStrand == '+')
				{
					pLinkMap->chGStrand = '-';
				}
				else if(pLinkMap->chHStrand == '-')
				{
					pLinkMap->chGStrand = '+';
				}
				else
				{
					pLinkMap->chGStrand = '?';
				}
			}
			else
			{
				strcpy(pLinkMap->strGenome, "GAP");
				pLinkMap->nGStart = -1;
				pLinkMap->nGEnd = -1;
				pLinkMap->chGStrand = '?';
			}
		}

		/* get next */
		pLinkMap = pLinkMap->pNext;
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_BlastHit_ExportLinkMapSeq: export sequence and link maps        */
/* ----------------------------------------------------------------------- */ 
int MAlign_BlastHit_ExportLinkMapSeq(struct tagSequence *pDataSeq, char strSeqAlias[],
			struct SEQLINKMAP *pSeqMapList, char strSeqOutFile[], char strCodOutFile[])
{
	/* define */
	FILE *fpSeq = NULL;
	FILE *fpCod = NULL;
	struct SEQLINKMAP *pLinkMap;
	struct tagString *pSeq = NULL;
	struct tagString *pTail;
	int ni;
	int nLinePos;
	char *pSeqPos,*chBase;
	
	/* init */
	if(pDataSeq == NULL)
	{
		return PROC_SUCCESS;	
	}

	fpSeq = fopen(strSeqOutFile, "w");
	if(fpSeq == NULL)
	{
		printf("Error: MAlign_BlastHit_ExportLinkMapSeq, cannot open sequence output file!\n");
		exit(EXIT_FAILURE);
	}
	fpCod = fopen(strCodOutFile, "w");
	if(fpCod == NULL)
	{
		printf("Error: MAlign_BlastHit_ExportLinkMapSeq, cannot open linkmap output file!\n");
		exit(EXIT_FAILURE);
	}

	
	/* export */
	fprintf(fpSeq, ">%s\n", strSeqAlias);
	pLinkMap = pSeqMapList;
	while(pLinkMap != NULL)
	{
		fprintf(fpCod, "%s\t%d\t%d\t%c\t%s\t%d\t%d\t%c\t%s\t%d\t%d\t%c\n", strSeqAlias,
			pLinkMap->nQStart, pLinkMap->nQEnd, pLinkMap->chQStrand,
			pLinkMap->strHit, pLinkMap->nHStart, pLinkMap->nHEnd, pLinkMap->chHStrand,
			pLinkMap->strGenome, pLinkMap->nGStart, pLinkMap->nGEnd, pLinkMap->chGStrand);
		
		if(strcmp(pLinkMap->strHit, "GAP") == 0)
		{
			pTail = NULL;
			pTail = CreateString(pLinkMap->nQEnd-pLinkMap->nQStart+1);
			for(ni=0; ni<pTail->m_nLength; ni++)
			{
				pTail->m_pString[ni] = 'N';
			}
			pTail->m_pString[ni] = '\0';

			StringAddTail(&pSeq, pTail->m_pString);
			DeleteString(pTail);
		}
		else
		{
			pTail = NULL;
			pTail = CreateString(pLinkMap->nQEnd-pLinkMap->nQStart+1);
			if(pLinkMap->chHStrand == '-')
			{
				chBase = pDataSeq->m_pSequence->m_pString+pLinkMap->nHEnd;
				pSeqPos = pTail->m_pString;
				for(ni=0; ni<pTail->m_nLength; ni++)
				{
					switch(*chBase)
					{
						case 'a': pSeqPos[ni] = 't';
							break;
						case 'A': pSeqPos[ni] = 'T';
							break;
						case 'c': pSeqPos[ni] = 'g';
							break;
						case 'C': pSeqPos[ni] = 'G';
							break;
						case 'g': pSeqPos[ni] = 'c';
							break;
						case 'G': pSeqPos[ni] = 'C';
							break;
						case 't': pSeqPos[ni] = 'a';
							break;
						case 'T': pSeqPos[ni] = 'A';
							break;
						case 'n': pSeqPos[ni] = 'n';
							break;
						default: pSeqPos[ni] = 'N';
					}
					chBase--;
				}

				pSeqPos[pTail->m_nLength] = '\0';
			}
			else
			{
				memcpy(pTail->m_pString, pDataSeq->m_pSequence->m_pString+pLinkMap->nHStart, pTail->m_nLength);
				pTail->m_pString[pTail->m_nLength] = '\0';
			}

			StringAddTail(&pSeq, pTail->m_pString);
			DeleteString(pTail);
		}

		/* get next */
		pLinkMap = pLinkMap->pNext;
	}

	/* write sequence */
	if(pSeq != NULL)
	{
		nLinePos = 0;
		pSeqPos = pSeq->m_pString;
		while(*pSeqPos != '\0')
		{
			fprintf(fpSeq, "%c", *pSeqPos);
			nLinePos++;
			pSeqPos++;

			if(nLinePos == FASTA_LINE_LEN)
			{
				fprintf(fpSeq, "\n");
				nLinePos = 0;
			}
		}
		if(nLinePos != 0)
		{
			fprintf(fpSeq, "\n");
		}

		DeleteString(pSeq);
	}
	

	/* close files */
	fclose(fpSeq);
	fclose(fpCod);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_MLAGAN_Align: create mlagan alignment                           */
/* ----------------------------------------------------------------------- */ 
int MAlign_MLAGAN_Align(char strWorkPath[], int nDatabaseNum,
			struct tagString **vDataAlias, 
			char strOutFile[], int nClusterId,
			char strMlaganPath[], char strMlaganParam[])
{
	/* define */
	char strCommand[LONG_LINE_LENGTH];
	char strInputFile[LONG_LINE_LENGTH];
	char strOutputFile[MED_LINE_LENGTH];
	char strFileName[LONG_LINE_LENGTH];
	int ni;

	/* align */
	strcpy(strInputFile, "");
	for(ni=0; ni<nDatabaseNum; ni++)
	{
		sprintf(strFileName, "%s%s_%d_%s.fa", strWorkPath, strOutFile, nClusterId, vDataAlias[ni]->m_pString);
		sprintf(strCommand, "%s %s", strInputFile, strFileName);
		strcpy(strInputFile, strCommand);
	}

	sprintf(strOutputFile, "%s%s_%d_mlagan", strWorkPath, strOutFile, nClusterId);
	sprintf(strCommand, "%smlagan %s %s >%s.out", strMlaganPath, strInputFile, strMlaganParam,
		strOutputFile);
	system(strCommand);

	sprintf(strFileName, "%sutils", strMlaganPath);
	AdjustDirectoryPath(strFileName);
	sprintf(strCommand, "%smpretty.pl %s.out >%s.aln", strFileName, strOutputFile, strOutputFile);
	system(strCommand);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_MotifMap_Main: MAlign_MotifMap main function                    */
/*  Map motifs to aligned orthologous seqeunces.                           */
/* ----------------------------------------------------------------------- */ 
int MAlign_MotifMap_Main(char strParamPath[])
{
	/* define */

	/* for loading parameters */
	FILE *fpIn;
	char strLine[MED_LINE_LENGTH];
	char strFileName[MED_LINE_LENGTH];
	char strMapOutName[MED_LINE_LENGTH];
	char strCommand[LONG_LINE_LENGTH];
	char *chSep,*chSep2;
	int nError = 0;
	int ni,nj;

	/* parameters */
	char strCisGenomePath[MED_LINE_LENGTH];
	char strWorkPath[MED_LINE_LENGTH];
	char strFileHeader[MED_LINE_LENGTH];
	int nClusterNum = 0;
	int nSpeciesNum = 0;
	int nMotifNum = 0;

	struct tagString **vSpeciesName = NULL;
	struct tagString **vMotifSymbol = NULL;
	struct tagString **vMotifName = NULL;
	struct tagString **vMotifFile = NULL;
	struct tagString **vMapProgram = NULL;

	/* ---------------------------------- */
	/* Step I: Load Parameters            */
	/* ---------------------------------- */
	fpIn = NULL;
	fpIn = fopen(strParamPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: MAlign_MotifMap_Main, cannot open the parameter file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strstr(strLine, "[Path of CisGenome]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strCisGenomePath, chSep);
			AdjustDirectoryPath(strCisGenomePath);
		}

		else if(strstr(strLine, "[Working Directory]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strWorkPath, chSep);
			AdjustDirectoryPath(strWorkPath);
		}

		else if(strstr(strLine, "[File Header]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			strcpy(strFileHeader, chSep);
			
			if(strcmp(strFileHeader, "") == 0)
			{
				printf("Error: file header not speciefied!\n");
				nError = 1;
				break;
			}
		}

		else if(strstr(strLine, "[Number of Clusters]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nClusterNum = atoi(chSep);
			}

			if(nClusterNum <= 0)
			{
				printf("Number of clusters = %d\n", nClusterNum);
				nError = 1;
				break;
			}
		}

		else if(strstr(strLine, "[Number of Species]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nSpeciesNum = atoi(chSep);
			}

			if(nSpeciesNum <= 0)
			{
				printf("Number of species = %d\n", nSpeciesNum);
				nError = 1;
				break;
			}

			vSpeciesName = (struct tagString **)calloc(nSpeciesNum, sizeof(struct tagString *));
			if(vSpeciesName == NULL)
			{
				printf("Error: MAlign_MotifMap_Main, cannot allocate memory for loading sequence alias!\n");
				nError = 1;
				break;
			}
		}

		else if(strstr(strLine, "[Species]") == strLine)
		{
			ni = 0;
			while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
			{
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] == '\0')
					continue;

				StringAddTail(vSpeciesName+ni, strLine);
			
				ni++;

				if(ni == nSpeciesNum)
					break;
			}
		}

		else if(strstr(strLine, "[Number of Motifs]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nMotifNum = atoi(chSep);
			}

			if(nMotifNum <= 0)
			{
				printf("Number of motifs = %d\n", nMotifNum);
				nError = 1;
				break;
			}

			vMotifSymbol = (struct tagString **)calloc(nMotifNum, sizeof(struct tagString *));
			if(vMotifSymbol == NULL)
			{
				printf("Error: MAlign_MotifMap_Main, cannot allocate memory for loading motif mapping parameters!\n");
				nError = 1;
				break;
			}

			vMotifName = (struct tagString **)calloc(nMotifNum, sizeof(struct tagString *));
			if(vMotifName == NULL)
			{
				printf("Error: MAlign_MotifMap_Main, cannot allocate memory for loading motif mapping parameters!\n");
				nError = 1;
				break;
			}
			
			vMotifFile = (struct tagString **)calloc(nMotifNum, sizeof(struct tagString *));
			if(vMotifFile == NULL)
			{
				printf("Error: MAlign_MotifMap_Main, cannot allocate memory for loading motif mapping parameters!\n");
				nError = 1;
				break;
			}

			vMapProgram = (struct tagString **)calloc(nMotifNum, sizeof(struct tagString *));
			if(vMapProgram == NULL)
			{
				printf("Error: MAlign_MotifMap_Main, cannot allocate memory for loading motif mapping parameters!\n");
				nError = 1;
				break;
			}

		}

		else if(strstr(strLine, "[Motif Mapping]") == strLine)
		{
			ni = 0;
			while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
			{
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] == '\0')
					continue;

				if(ni >= nMotifNum)
				{
					printf("Error: MAlign_MotifMap_Main, database number does not match the database provided!\n");
					nError = 1;
					break;
				}

				chSep2 = strchr(strLine, '\t');
				*chSep2 = '\0';
				StringAddTail(vMotifSymbol+ni, strLine);

				chSep = chSep2+1;
				chSep2 = strchr(chSep, '\t');
				*chSep2 = '\0';
				StringAddTail(vMotifName+ni, chSep);

				chSep = chSep2+1;
				chSep2 = strchr(chSep, '\t');
				*chSep2 = '\0';
				StringAddTail(vMotifFile+ni, chSep);

				chSep = chSep2+1;
				StringAddTail(vMapProgram+ni, chSep);
			
				ni++;
			}

			if( (ni != nMotifNum) || (nError == 1) ) 
			{
				printf("Error: MAlign_MotifMap_Main, motif number does not match the motif mapping parameters provided!\n");
				nError = 1;
				break;
			}
		}
		
		else if(strLine[0] == '#')
		{
		}

		else
		{
			printf("Error: MAlign_MotifMap_Main, unknown parameters!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	fclose(fpIn);

	if(nError > 0)
	{
		printf("Error: MAlign_MotifMap_Main, cannot load the parameter file correctly!\n");
		exit(EXIT_FAILURE);
	}

	/* ------------------------------------ */
	/* Step II: Prepare Fasta File          */
	/* ------------------------------------ */
	for(ni=0; ni<nSpeciesNum; ni++)
	{
		MAlign_MotifMap_PrepareFasta(strWorkPath, strFileHeader, 
			nClusterNum, vSpeciesName[ni]->m_pString);
	}

	/* ------------------------------------ */
	/* Step III: Map motifs to each species */
	/* ------------------------------------ */
	for(ni=0; ni<nSpeciesNum; ni++)
	{
		for(nj=0; nj<nMotifNum; nj++)
		{
			/* ------------------------------------ */
			/* Step III.1: Map motifs               */
			/* ------------------------------------ */
			sprintf(strFileName, "%s_%s_combined.fa", strFileHeader, 
				vSpeciesName[ni]->m_pString);
			sprintf(strMapOutName, "%s_%s_%s.map.tmp", strFileHeader, 
				vSpeciesName[ni]->m_pString, vMotifName[nj]->m_pString);
			sprintf(strCommand, "%s%s -m %s%s -d %s -i %s -o %s%s", strCisGenomePath,
				vMapProgram[nj]->m_pString, strWorkPath, vMotifFile[nj]->m_pString, 
				strWorkPath, strFileName, strWorkPath, strMapOutName);
			system(strCommand);

			/* ------------------------------------ */
			/* Step III.2: compute genomic cods.    */
			/* ------------------------------------ */
			MAlign_MotifMap_ComputeGenomicCod(strWorkPath, strMapOutName, 
				strFileHeader, nClusterNum, vSpeciesName[ni]->m_pString,
				vMotifName[nj]->m_pString);
		}

		/* ------------------------------------ */
		/* Step III.3: divide into clusters.    */
		/* ------------------------------------ */
		MAlign_MotifMap_GroupMotifSiteByClusters(strWorkPath, strFileHeader, 
			vSpeciesName[ni]->m_pString, nClusterNum, 
			nMotifNum, vMotifSymbol, vMotifName);
	}

	/* ------------------------------------ */
	/* Step IV: Map motifs to alignment     */
	/* ------------------------------------ */
	for(ni=0; ni<nClusterNum; ni++)
	{
		MAlign_MotifMap_GenerateAlignment(strWorkPath, strFileHeader, 
			ni, nSpeciesNum, vSpeciesName,  
			nMotifNum, vMotifSymbol, vMotifName);
	}

	/* ------------------------------------ */
	/* Step V: Release memory               */
	/* ------------------------------------ */
	for(ni=0; ni<nSpeciesNum; ni++)
	{
		DeleteString(vSpeciesName[ni]);
	}
	free(vSpeciesName);

	for(ni=0; ni<nMotifNum; ni++)
	{
		DeleteString(vMotifSymbol[ni]);
		DeleteString(vMotifName[ni]);
		DeleteString(vMotifFile[ni]);
		DeleteString(vMapProgram[ni]);
	}
	free(vMotifSymbol);
	free(vMotifName);
	free(vMotifFile);
	free(vMapProgram);

	/* ---------------------------------- */
	/* Step VI: Remove Files               */
	/* ---------------------------------- */
	sprintf(strFileName, "%s%s*.map.tmp*", strWorkPath, strFileHeader);
	RemoveFiles(strFileName);
	
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_MotifMap_PrepareFasta: prepare a single fasta file for each     */
/*  species used for MAlign_MotifMap.                                      */
/* ----------------------------------------------------------------------- */ 
int MAlign_MotifMap_PrepareFasta(char strWorkPath[], char strFileHeader[], 
			int nClusterNum, char strSpeciesName[])
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	int ni;
	char strLine[LONG_LINE_LENGTH];
	char strOutFile[MED_LINE_LENGTH];
	char strFileName[MED_LINE_LENGTH];
	
	/* init */
	sprintf(strOutFile, "%s%s_%s_combined.fa", strWorkPath, 
		strFileHeader, strSpeciesName);
	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: MAlign_MotifMap_PrepareFasta, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	/* combine */
	for(ni=0; ni<nClusterNum; ni++)
	{
		sprintf(strFileName, "%s%s_%d_%s.fa", strWorkPath, 
			strFileHeader, ni, strSpeciesName);

		fpIn = NULL;
		fpIn = fopen(strFileName, "r");
		if(fpIn == NULL)
		{
			printf("Error: MAlign_MotifMap_PrepareFasta, cannot open input sequence file!\n");
			exit(EXIT_FAILURE);
		}

		while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;

			if(strLine[0] == '>')
			{
				fprintf(fpOut, ">%d\n", ni);
			}
			else
			{
				fprintf(fpOut, "%s\n", strLine);
			}
		}

		fclose(fpIn);
	}

	/* close files */
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_MotifMap_GroupMotifSiteByClusters: group motif sites according  */
/*  to sequence clusters.                                                  */
/* ----------------------------------------------------------------------- */ 
int MAlign_MotifMap_GroupMotifSiteByClusters(char strWorkPath[], char strFileHeader[], 
			char strSpeciesName[], int nClusterNum, int nMotifNum, 
			struct tagString **vMotifSymbol, struct tagString **vMotifName)
{
	/* define */
	FILE *fpIn;
	FILE **vfpOut;
	int ni,nId;
	char strFileName[MED_LINE_LENGTH];
	char strLine[MED_LINE_LENGTH];
	char *chp,*chp2;

	/* init */
	vfpOut = NULL;
	vfpOut = (FILE **)calloc(nClusterNum, sizeof(FILE *));
	if(vfpOut == NULL)
	{
		printf("Error: MAlign_MotifMap_GroupMotifSiteByClusters, cannot open output file, memory not enough!\n");
		printf("%d\n", nClusterNum);
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nClusterNum; ni++)
	{
		sprintf(strFileName, "%s%s_%d_%s.map", strWorkPath, strFileHeader, ni, strSpeciesName);
		vfpOut[ni] = fopen(strFileName, "w");
		if(vfpOut[ni] == NULL)
		{
			printf("Error: MAlign_MotifMap_GroupMotifSiteByClusters, cannot open output file!\n");
			printf("%d\t%d\n", ni, nClusterNum);
			exit(EXIT_FAILURE);
		}
	}

	/* divide */
	for(ni=0; ni<nMotifNum; ni++)
	{
		sprintf(strFileName, "%s%s_%s_%s.map", strWorkPath, strFileHeader, 
			strSpeciesName, vMotifName[ni]->m_pString);
		fpIn = NULL;
		fpIn = fopen(strFileName, "r");
		if(fpIn == NULL)
		{
			printf("Error: MAlign_MotifMap_GroupMotifSiteByClusters, cannot open input file!\n");
			exit(EXIT_FAILURE);
		}

		while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;

			sscanf(strLine, "%d", &nId);
			
			if(nId < nClusterNum)
			{
				chp = strchr(strLine, '\t');
				chp2 = chp+1;
				chp = strchr(chp2, '\t');
				chp2 = chp+1;
				fprintf(vfpOut[nId], "%s\t%s\t%s\n", 
					vMotifSymbol[ni]->m_pString, vMotifName[ni]->m_pString, chp2);
			}
		}

		fclose(fpIn);
	}


	/* close files */
	for(ni=0; ni<nClusterNum; ni++)
	{
		fclose(vfpOut[ni]);
	}
	free(vfpOut);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_MotifMap_ComputeGenomicCod: compute genomic coordinates of      */
/*  mapped sites.                                                          */
/* ----------------------------------------------------------------------- */ 
int MAlign_MotifMap_ComputeGenomicCod(char strWorkPath[], char strOriginMap[], 
				char strFileHeader[], int nClusterNum, char strSpeciesName[],
				char strMotifName[])
{
	/* define */
	struct MOTIFSITELINKMAP *pSiteList = NULL;
	struct MOTIFSITELINKMAP *pSite,*pPrev;
	FILE *fpIn;
	FILE *fpOut;
	char strLine[MED_LINE_LENGTH];
	char strFileName[MED_LINE_LENGTH];
	int nClusterId = -1;
	int nFindRef;

	/* coordinate map */
	struct SEQLINKMAP *pRefCodList = NULL;
	struct SEQLINKMAP *pRefCod,*pRefPrev;

	/* load all sites */
	sprintf(strFileName, "%s%s", strWorkPath, strOriginMap);
	fpIn = NULL;
	fpIn = fopen(strFileName, "r");
	if(fpIn == NULL)
	{
		printf("Error: MAlign_MotifMap_ComputeGenomicCod, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		pSite = NULL;
		pSite = MOTIFSITELINKMAPCREATE();
		if(pSite == NULL)
		{
			printf("Error: MAlign_MotifMap_ComputeGenomicCod, cannot create motif site link map!\n");
			exit(EXIT_FAILURE);
		}

		sscanf(strLine, "%d %s %d %d %c %lf %s", &(pSite->nSeqId),
			pSite->strQuery, &(pSite->nQStart), &(pSite->nQEnd),
			&(pSite->chQStrand), &(pSite->dScore), pSite->strSiteSeq);

		if(pSiteList == NULL)
		{
			pSiteList = pSite;
			pPrev = pSite;
		}
		else
		{
			pPrev->pNext = pSite;
			pPrev = pSite;
		}
	}

	fclose(fpIn);

	/* compute relative coordination */
	sprintf(strFileName, "%s%s_%s_%s.map", strWorkPath, strFileHeader, strSpeciesName, strMotifName);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: MAlign_MotifMap_ComputeGenomicCod, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	nClusterId = -1;
	pSite = pSiteList;
	while(pSite != NULL)
	{
		/* update reference coordinate */
		if(pSite->nSeqId != nClusterId)
		{
			if(pRefCodList != NULL)
			{
				SEQLINKMAPCLEARLIST(&pRefCodList);
			}

			nClusterId = pSite->nSeqId;
			pRefCodList = NULL;

			sprintf(strFileName, "%s%s_%d_%s.gln", strWorkPath, strFileHeader, nClusterId, strSpeciesName);
			fpIn = NULL;
			fpIn = fopen(strFileName, "r");
			if(fpIn == NULL)
			{
				printf("Error: MAlign_MotifMap_ComputeGenomicCod, cannot open output file!\n");
				exit(EXIT_FAILURE);
			}

			while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
			{
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] == '\0')
					continue;

				pRefCod = NULL;
				pRefCod = SEQLINKMAPCREATE();
				if(pRefCod == NULL)
				{
					printf("Error: MAlign_MotifMap_ComputeGenomicCod, cannot create reference coordinates!\n");
					exit(EXIT_FAILURE);
				}

				sscanf(strLine, "%s %d %d %c %s %d %d %c %s %d %d %c", 
					pRefCod->strQuery, &(pRefCod->nQStart), &(pRefCod->nQEnd), &(pRefCod->chQStrand),
					pRefCod->strHit, &(pRefCod->nHStart), &(pRefCod->nHEnd), &(pRefCod->chHStrand),
					pRefCod->strGenome, &(pRefCod->nGStart), &(pRefCod->nGEnd), &(pRefCod->chGStrand));

				if(pRefCodList == NULL)
				{
					pRefCodList = pRefCod;
					pRefPrev = pRefCod;
				}
				else
				{
					pRefPrev->pNext = pRefCod;
					pRefPrev = pRefCod;
				}
			}

			fclose(fpIn);	
		}

		/* compute coordinates */
		if(pRefCodList == NULL)
		{
			strcpy(pSite->strHit, "GAP");
			pSite->nHStart = -1;
			pSite->nHEnd = -1;
			pSite->chHStrand = '?';
			
			strcpy(pSite->strGenome, "GAP");
			pSite->nGStart = -1;
			pSite->nGEnd = -1;
			pSite->chGStrand = '?';
		}
		else
		{
			nFindRef = 0;
			pRefCod = pRefCodList;
			while(pRefCod != NULL)
			{
				if( (pSite->nQStart>=pRefCod->nQStart) && (pSite->nQEnd <= pRefCod->nQEnd) )
				{
					nFindRef = 1;
					break;
				}

				/* get next */
				pRefCod = pRefCod->pNext;
			}


			if(nFindRef == 1)
			{
				/* check query */
				if(pRefCod->chQStrand != '+')
				{
					printf("Error: MAlign_MotifMap_ComputeGenomicCod, original strand must be +!\n");
					exit(EXIT_FAILURE);
				}

				/* adjust hit */
				if(strcmp(pRefCod->strHit, "GAP") == 0)
				{
					strcpy(pSite->strHit, "GAP");
					pSite->nHStart = -1;
					pSite->nHEnd = -1;
					pSite->chHStrand = '?';
				}
				else
				{
					if(pRefCod->chHStrand == '+')
					{
						strcpy(pSite->strHit, pRefCod->strHit);
						pSite->nHStart = pSite->nQStart-pRefCod->nQStart+pRefCod->nHStart;
						pSite->nHEnd = pSite->nQEnd-pRefCod->nQStart+pRefCod->nHStart;
						if(pSite->chQStrand == '+')
						{
							pSite->chHStrand = '+';
						}
						else if(pSite->chQStrand == '-')
						{
							pSite->chHStrand = '-';
						}
						else
						{
							pSite->chHStrand = '?';
						}
					}
					else if(pRefCod->chHStrand == '-')
					{
						strcpy(pSite->strHit, pRefCod->strHit);
						pSite->nHStart = pRefCod->nHEnd-(pSite->nQEnd-pRefCod->nQStart);
						pSite->nHEnd = pRefCod->nHEnd-(pSite->nQStart-pRefCod->nQStart);
						if(pSite->chQStrand == '+')
						{
							pSite->chHStrand = '-';
						}
						else if(pSite->chQStrand == '-')
						{
							pSite->chHStrand = '+';
						}
						else
						{
							pSite->chHStrand = '?';
						}
					}
					else
					{
						strcpy(pSite->strHit, "GAP");
						pSite->nHStart = -1;
						pSite->nHEnd = -1;
						pSite->chHStrand = '?';
					}
				}

				/* adjust genome */
				if(strcmp(pRefCod->strGenome, "GAP") == 0)
				{
					strcpy(pSite->strGenome, "GAP");
					pSite->nGStart = -1;
					pSite->nGEnd = -1;
					pSite->chGStrand = '?';
				}
				else
				{
					if(pRefCod->chGStrand == '+')
					{
						strcpy(pSite->strGenome, pRefCod->strGenome);
						pSite->nGStart = pSite->nQStart-pRefCod->nQStart+pRefCod->nGStart;
						pSite->nGEnd = pSite->nQEnd-pRefCod->nQStart+pRefCod->nGStart;
						if(pSite->chQStrand == '+')
						{
							pSite->chGStrand = '+';
						}
						else if(pSite->chQStrand == '-')
						{
							pSite->chGStrand = '-';
						}
						else
						{
							pSite->chGStrand = '?';
						}
					}
					else if(pRefCod->chGStrand == '-')
					{
						strcpy(pSite->strGenome, pRefCod->strGenome);
						pSite->nGStart = pRefCod->nGEnd-(pSite->nQEnd-pRefCod->nQStart);
						pSite->nGEnd = pRefCod->nGEnd-(pSite->nQStart-pRefCod->nQStart);
						if(pSite->chQStrand == '+')
						{
							pSite->chGStrand = '-';
						}
						else if(pSite->chQStrand == '-')
						{
							pSite->chGStrand = '+';
						}
						else
						{
							pSite->chGStrand = '?';
						}
					}
					else
					{
						strcpy(pSite->strGenome, "GAP");
						pSite->nGStart = -1;
						pSite->nGEnd = -1;
						pSite->chGStrand = '?';
					}
				}
			}
			else
			{
				strcpy(pSite->strHit, "GAP");
				pSite->nHStart = -1;
				pSite->nHEnd = -1;
				pSite->chHStrand = '?';
				
				strcpy(pSite->strGenome, "GAP");
				pSite->nGStart = -1;
				pSite->nGEnd = -1;
				pSite->chGStrand = '?';
			}
		}

		/* write */
		fprintf(fpOut, "%d\t%s\t%d\t%d\t%c\t%f\t%s\t%s\t%d\t%d\t%c\t%s\t%d\t%d\t%c\n", 
			pSite->nSeqId, pSite->strQuery, pSite->nQStart, pSite->nQEnd,
			pSite->chQStrand, pSite->dScore, pSite->strSiteSeq,
			pSite->strHit, pSite->nHStart, pSite->nHEnd, pSite->chHStrand,
			pSite->strGenome, pSite->nGStart, pSite->nGEnd, pSite->chGStrand);

		/* get next */
		pSite = pSite->pNext;
	}

	fclose(fpOut);

	/* clear reference coordinate */
	if(pRefCodList != NULL)
	{
		SEQLINKMAPCLEARLIST(&pRefCodList);
	}

	/* clear list */
	MOTIFSITELINKMAPCLEARLIST(&pSiteList);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_MotifMap_GenerateAlignment: generate alignment including        */
/*  mapped sites.                                                          */
/* ----------------------------------------------------------------------- */ 
int MAlign_MotifMap_GenerateAlignment(char strWorkPath[], char strFileHeader[], 
			int nClusterId, int nSpeciesNum, struct tagString **vSpeciesName,  
			int nMotifNum, struct tagString **vMotifSymbol, 
			struct tagString **vMotifName)
{
	/* define */
	int nSeqNum;
	struct tagSequence *pAln;
	struct tagString **vMS;
	struct tagSequence *pSeq;
	char chBase;
	struct SEQLINKMAP **vSeqLinkMap;
	char strFileName[MED_LINE_LENGTH];
	int ni,nj,nLen;
	char *vS;
	int nAlnLen;
	int nIdn;

	/* load alignment */
	sprintf(strFileName, "%s%s_%d_mlagan.out", strWorkPath, strFileHeader, 
		nClusterId);
	pAln = NULL;
	nSeqNum = 0;
	nSeqNum = LoadFullSequenceList(strFileName, &pAln);
	if(nSeqNum != nSpeciesNum)
	{
		printf("Error: MAlign_MotifMap_GenerateAlignment, sequence number does not match!\n");
		exit(EXIT_FAILURE);
	}
	if(nSeqNum <= 1)
	{
		printf("Error: MAlign_MotifMap_GenerateAlignment, <=1 sequence available, require at least two sequences!\n");
		exit(EXIT_FAILURE);
	}

	/* init motif status */
	vMS = NULL;
	vMS = (struct tagString **)calloc(nSeqNum, sizeof(struct tagString*));
	if(vMS == NULL)
	{
		printf("Error: MAlign_MotifMap_GenerateAlignment, cannot allocate space for recording motif status!\n");
		exit(EXIT_FAILURE);
	}

	pSeq = pAln;
	for(ni=0; ni<nSeqNum; ni++)
	{
		nLen = pSeq->m_nLength;
		vMS[ni] = NULL;
		vMS[ni] = CreateString(nLen);
		if(vMS[ni] == NULL)
		{
			printf("Error: MAlign_MotifMap_GenerateAlignment, cannot allocate space for recording motif status!\n");
			exit(EXIT_FAILURE);
		}

		vS = vMS[ni]->m_pString;
		for(nj=0; nj<nLen; nj++)
		{
			vS[nj] = ' ';
		}
		vS[nj] = '\0';

		pSeq = pSeq->m_pNext;
	}

	for(ni=0; ni<nSpeciesNum; ni++)
	{
		if(vMS[ni]->m_nLength != vMS[0]->m_nLength)
		{
			printf("Error: MAlign_MotifMap_GenerateAlignment, sequence length not the same!\n");
			exit(EXIT_FAILURE);
		}
	}
	nAlnLen = pAln->m_nLength;


	/* load seqlink map */
	vSeqLinkMap = NULL;
	vSeqLinkMap = (struct SEQLINKMAP **)calloc(nSeqNum, sizeof(struct SEQLINKMAP *));
	if(vSeqLinkMap == NULL)
	{
		printf("Error: MAlign_MotifMap_GenerateAlignment, cannot allocate space for recording sequence link status!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nSeqNum; ni++)
	{
		sprintf(strFileName, "%s%s_%d_%s.gln", strWorkPath, strFileHeader, 
				nClusterId, vSpeciesName[ni]->m_pString);
		MAlign_MotifMap_LoadSeqLinkMap(vSeqLinkMap+ni, strFileName);
	}

	
	/* find identity */
	for(ni=0; ni<nAlnLen; ni++)
	{
		nIdn = 1;
		pSeq = pAln;
		chBase = pSeq->m_pSequence->m_pString[ni];
		if(chBase == '-')
		{
			nIdn = 0;
		}
		else
		{
			for(nj=1; nj<nSpeciesNum; nj++)
			{
				pSeq = pSeq->m_pNext;
				if(pSeq->m_pSequence->m_pString[ni] == chBase)
				{
					vMS[nj]->m_pString[ni] = '.';
				}
				else
				{
					nIdn = 0;
				}
			}
		}

		if(nIdn == 1)
		{
			for(nj=1; nj<nSpeciesNum; nj++)
			{
				vMS[nj]->m_pString[ni] = ':';
			}
		}
	}


	/* update motif status */
	pSeq = pAln;
	for(ni=0; ni<nSpeciesNum; ni++)
	{
		MAlign_MotifMap_UpdateAlignAndMotif(pSeq, vMS[ni], strWorkPath, strFileHeader, 
			nClusterId, vSpeciesName[ni]->m_pString,  
			nMotifNum, vMotifSymbol, vMotifName);
		pSeq = pSeq->m_pNext;
	}

	/* export */
	sprintf(strFileName, "%s%s_%d_motif.aln", strWorkPath, strFileHeader, nClusterId);
	MAlign_MotifMap_ExportAlignment(nSpeciesNum, vSpeciesName, pAln, vMS, vSeqLinkMap,
		strFileName, nMotifNum, vMotifSymbol, vMotifName);

	/* clear memory */
	for(ni=0; ni<nSeqNum; ni++)
	{
		DeleteString(vMS[ni]);
		vMS[ni] = NULL;
		SEQLINKMAPCLEARLIST(vSeqLinkMap+ni);
	}
	free(vMS);
	free(vSeqLinkMap);
	SequenceListClear(&pAln);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_MotifMap_LoadSeqLinkMap: load sequence link map.                */
/* ----------------------------------------------------------------------- */ 
int MAlign_MotifMap_LoadSeqLinkMap(struct SEQLINKMAP **pLinkMapList, char strFileName[])
{
	/* define */
	FILE *fpIn;
	char strLine[MED_LINE_LENGTH]; 
	struct SEQLINKMAP *pRefCod,*pRefPrev;

	/* init */
	if(pLinkMapList == NULL)
	{
		printf("Warning: MAlign_MotifMap_LoadSeqLinkMap, empty link map list!\n");
		return PROC_SUCCESS;
	}
	
	/* load */
	fpIn = NULL;
	fpIn = fopen(strFileName, "r");
	if(fpIn == NULL)
	{
		printf("Error: MAlign_MotifMap_LoadSeqLinkMap, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		pRefCod = NULL;
		pRefCod = SEQLINKMAPCREATE();
		if(pRefCod == NULL)
		{
			printf("Error: MAlign_MotifMap_LoadSeqLinkMap, cannot create reference coordinates!\n");
			exit(EXIT_FAILURE);
		}

		sscanf(strLine, "%s %d %d %c %s %d %d %c %s %d %d %c", 
			pRefCod->strQuery, &(pRefCod->nQStart), &(pRefCod->nQEnd), &(pRefCod->chQStrand),
			pRefCod->strHit, &(pRefCod->nHStart), &(pRefCod->nHEnd), &(pRefCod->chHStrand),
			pRefCod->strGenome, &(pRefCod->nGStart), &(pRefCod->nGEnd), &(pRefCod->chGStrand));

		if(*pLinkMapList == NULL)
		{
			*pLinkMapList = pRefCod;
			pRefPrev = pRefCod;
		}
		else
		{
			pRefPrev->pNext = pRefCod;
			pRefPrev = pRefCod;
		}
	}

	fclose(fpIn);	

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_MotifMap_UpdateAlignAndMotif: update alignment to reflect motif */
/*  sites and repeats.                                                     */
/* ----------------------------------------------------------------------- */ 
int MAlign_MotifMap_UpdateAlignAndMotif(struct tagSequence *pSeq, 
			struct tagString *pMS, char strWorkPath[], char strFileHeader[], 
			int nClusterId, char strSpeciesName[],  
			int nMotifNum, struct tagString **vMotifSymbol, 
			struct tagString **vMotifName)
{
	/* define */
	FILE *fpIn;
	struct tagSequence *pOriginSeq;
	struct tagString *pOriginStatus;
	struct DOUBLEMATRIX *pOriginScore;
	struct MOTIFSITELINKMAP *pSiteList = NULL;
	struct MOTIFSITELINKMAP *pSite;
	char strLine[MED_LINE_LENGTH];
	char strFileName[MED_LINE_LENGTH];
	int nCount;
	char *chA,*chS;
	char *chOA,*chOS;
	double *pScore;
	int ni,nAlnLen,nPos;

	/* load sequences */
	sprintf(strFileName, "%s%s_%d_%s.fa", strWorkPath, strFileHeader, 
		nClusterId, strSpeciesName);
	pOriginSeq = NULL;
	nCount = LoadFullSequenceList(strFileName, &pOriginSeq);
	if( (nCount != 1) || (pOriginSeq == NULL))
	{
		printf("Error: MAlign_MotifMap_UpdateAlignAndMotif, wrong number of sequence loaded!\n");
		exit(EXIT_FAILURE);
	}
	pOriginStatus = NULL;
	pOriginStatus = CreateString(pOriginSeq->m_nLength);
	pOriginScore = NULL;
	pOriginScore = CreateDoubleMatrix(1, pOriginSeq->m_nLength);
	if( (pOriginStatus == NULL) || (pOriginScore == NULL) )
	{
		printf("Error: MAlign_MotifMap_UpdateAlignAndMotif, cannot create memory for origin status&score!\n");
		exit(EXIT_FAILURE);
	}
	chS = pOriginStatus->m_pString;
	pScore = pOriginScore->pMatElement;
	for(ni=0; ni<pOriginSeq->m_nLength; ni++)
	{
		chS[ni] = ' ';
		pScore[ni] = -DM_ACCESS_VIOLATION;
	}

	/* load sites */
	sprintf(strFileName, "%s%s_%d_%s.map", strWorkPath, strFileHeader, 
		nClusterId, strSpeciesName);
	fpIn = NULL;
	fpIn = fopen(strFileName, "r");
	if(fpIn == NULL)
	{
		printf("Error: MAlign_MotifMap_UpdateAlignAndMotif, cannot open site coordinates file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		pSite = NULL;
		pSite = MOTIFSITELINKMAPCREATE();
		if(pSite == NULL)
		{
			printf("Error: MAlign_MotifMap_UpdateAlignAndMotif, cannot create site link map!\n");
			exit(EXIT_FAILURE);
		}

		sscanf(strLine, "%s %s %d %d %c %lf %s %s %d %d %c %s %d %d %c",
			pSite->strMotifSymbol, pSite->strMotifName, 
			&(pSite->nQStart), &(pSite->nQEnd), &(pSite->chQStrand),
			&(pSite->dScore), &(pSite->strSiteSeq),
			pSite->strHit, &(pSite->nHStart), &(pSite->nHEnd), &(pSite->chHStrand),
			pSite->strGenome, &(pSite->nGStart), &(pSite->nGEnd), &(pSite->chGStrand));

		strcpy(pSite->strQuery, strSpeciesName);

		MOTIFSITELINKMAP_INSERTBYQUERYLOC(pSite, &pSiteList);
	}

	fclose(fpIn);

	/* update status */
	while(pSiteList != NULL)
	{
		pSite = pSiteList;
		pSiteList = pSite->pNext;
		pSite->pNext = NULL;
		
		for(ni=pSite->nQStart; ni<=pSite->nQEnd; ni++)
		{
			if(pSite->dScore > pOriginScore->pMatElement[ni])
			{
				pOriginStatus->m_pString[ni] = pSite->strMotifSymbol[0];
				pOriginScore->pMatElement[ni] = pSite->dScore;
			}
		}
		
		MOTIFSITELINKMAPDESTROY(pSite);
	}


	/* generate aln status */
	nAlnLen = pSeq->m_nLength;
	chA = pSeq->m_pSequence->m_pString;
	chS = pMS->m_pString;
	chOA = pOriginSeq->m_pSequence->m_pString;
	chOS = pOriginStatus->m_pString;
	nPos = 0;
	for(ni=0; ni<nAlnLen; ni++)
	{
		if(chA[ni] != '-')
		{
			chA[ni] = chOA[nPos];
			if(chOS[nPos] != ' ')
				chS[ni] = chOS[nPos];
			nPos++;
		}
	}

	if(pSiteList != NULL)
	{
		printf("Error: MAlign_MotifMap_UpdateAlignAndMotif, not all sites are added!\n");
		exit(EXIT_FAILURE);
	}

	/* release memory */
	SequenceDelete(pOriginSeq);
	DeleteString(pOriginStatus);
	DestroyDoubleMatrix(pOriginScore);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_MotifMap_ExportAlignment: export alignment.                     */
/* ----------------------------------------------------------------------- */ 
int MAlign_MotifMap_ExportAlignment(int nSpeciesNum, 
		struct tagString **vSpeciesName, struct tagSequence *pAln, 
		struct tagString **vMS, struct SEQLINKMAP **vSeqLinkMap,
		char strOutFile[], int nMotifNum, struct tagString **vMotifSymbol, 
		struct tagString **vMotifName)
{
	/* define */
	FILE *fpOut;
	struct INTMATRIX *pPos;
	struct tagSequence *pCurrentAln;
	char *vSeq,*vSeq2;
	int ni,nj;
	int nLineLen;
	int nAlnLen;
	int nRemainLen;
	char strLine[LINE_LENGTH];
	int nRefStart;
	int nPos;

	/* init */
	if(nSpeciesNum <= 0)
	{
		printf("Warning: MAlign_MotifMap_ExportAlignment, no alignment sequence!\n");
		return PROC_SUCCESS;
	}

	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: MAlign_MotifMap_ExportAlignment, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	pPos = NULL;
	pPos = CreateIntMatrix(1, nSpeciesNum);
	if(pPos == NULL)
	{
		printf("Error: MAlign_MotifMap_ExportAlignment, cannot track alignment position!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nSpeciesNum; ni++)
	{
		if(vMS[ni]->m_nLength != vMS[0]->m_nLength)
		{
			printf("Error: MAlign_MotifMap_ExportAlignment, sequence length not the same!\n");
			exit(EXIT_FAILURE);
		}
	}
	nAlnLen = pAln->m_nLength;

	/* export */
	nRefStart = 0;
	nRemainLen = nAlnLen;
	while(nRemainLen > 0)
	{
		nLineLen = ALIGN_LINE_LEN;
		if(nLineLen > nRemainLen)
			nLineLen = nRemainLen;

		pCurrentAln = pAln;
		for(ni=0; ni<nSpeciesNum; ni++)
		{
			/* export motif site */
			if( (ni == 0) && (vSeqLinkMap[0] != NULL) )
			{
				strcpy(strLine, vSeqLinkMap[0]->strGenome);
				MAlign_MotifMap_TRUNCSTRING(strLine, ALIGN_NAME_LEN);
				fprintf(fpOut, "%s ", strLine);

				nPos = pPos->pMatElement[ni];
				nPos = MAlign_MotifMap_GetGenomicCod(nPos, vSeqLinkMap[0]); 
				sprintf(strLine, "%d", nPos);
				MAlign_MotifMap_TRUNCSTRING(strLine, ALIGN_START_LEN);
				fprintf(fpOut, "%s ", strLine);
			}
			else
			{
				strcpy(strLine, "");
				MAlign_MotifMap_TRUNCSTRING(strLine, ALIGN_NAME_LEN);
				fprintf(fpOut, "%s ", strLine);

				sprintf(strLine, "");
				MAlign_MotifMap_TRUNCSTRING(strLine, ALIGN_START_LEN);
				fprintf(fpOut, "%s ", strLine);
			}

			vSeq = vMS[ni]->m_pString+nRefStart;
			for(nj=0; nj<nLineLen; nj++)
			{
				fprintf(fpOut, "%c", vSeq[nj]);
			}

			fprintf(fpOut, " ");

			if( (ni == 0) && (vSeqLinkMap[0] != NULL) )
			{
				nPos = pPos->pMatElement[ni];
				vSeq2 = pCurrentAln->m_pSequence->m_pString+nRefStart;
				for(nj=0; nj<nLineLen; nj++)
				{
					if(vSeq2[nj] != '-')
						nPos += 1;
				}
				nPos -= 1;
				if(nPos < pPos->pMatElement[ni])
					nPos = pPos->pMatElement[ni];

				nPos = MAlign_MotifMap_GetGenomicCod(nPos, vSeqLinkMap[0]); 
				sprintf(strLine, "%d", nPos);
				MAlign_MotifMap_TRUNCSTRING(strLine, ALIGN_START_LEN);
				fprintf(fpOut, "%s\n", strLine);
			}
			else
			{
				sprintf(strLine, "");
				MAlign_MotifMap_TRUNCSTRING(strLine, ALIGN_START_LEN);
				fprintf(fpOut, "%s\n", strLine);
			}

			/* export sequence name */
			strcpy(strLine, vSpeciesName[ni]->m_pString);
			MAlign_MotifMap_TRUNCSTRING(strLine, ALIGN_NAME_LEN);
			fprintf(fpOut, "%s ", strLine);

			/* export start site */
			sprintf(strLine, "%d", pPos->pMatElement[ni]);
			MAlign_MotifMap_TRUNCSTRING(strLine, ALIGN_START_LEN);
			fprintf(fpOut, "%s ", strLine);

			/* export alignment */
			vSeq = pCurrentAln->m_pSequence->m_pString+nRefStart;
			for(nj=0; nj<nLineLen; nj++)
			{
				fprintf(fpOut, "%c", vSeq[nj]);
				if(vSeq[nj] != '-')
					pPos->pMatElement[ni] += 1;
			}

			/* export end site */
			fprintf(fpOut, " ");

			sprintf(strLine, "%d", pPos->pMatElement[ni]);
			MAlign_MotifMap_TRUNCSTRING(strLine, ALIGN_END_LEN);
			fprintf(fpOut, "%s\n", strLine);

			/* get next */
			pCurrentAln = pCurrentAln->m_pNext;
		}

		fprintf(fpOut, "\n\n");
		nRemainLen -= nLineLen;
		nRefStart += nLineLen;
	}

	/* close file */
	fclose(fpOut);
	DestroyIntMatrix(pPos);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_MotifMap_TRUNCSTRING: truncate string to a given length.        */
/* ----------------------------------------------------------------------- */ 
int MAlign_MotifMap_TRUNCSTRING(char strLine[], int nTruncLen)
{
	/* define */
	int nLen;
	int ni;

	/* truncation */
	nLen = strlen(strLine);
	if(nLen > nTruncLen)
	{
		strLine[nTruncLen] = '\0';
	}
	else
	{
		for(ni=nLen; ni<nTruncLen; ni++)
		{
			strLine[ni] = ' ';
		}
		strLine[nTruncLen] = '\0';
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_MotifMap_GetGenomicCod: get genomic coordinate of a given pos.  */
/* ----------------------------------------------------------------------- */ 
int MAlign_MotifMap_GetGenomicCod(int nPos, struct SEQLINKMAP *pSeqLinkMap)
{
	/* define */
	int nGPos = -1;
	int nFindRef;
	struct SEQLINKMAP *pRefCod;

	/* compute */
	nFindRef = 0;
	pRefCod = pSeqLinkMap;
	while(pRefCod != NULL)
	{
		if( (nPos>=pRefCod->nQStart) && (nPos <= pRefCod->nQEnd) )
		{
			nFindRef = 1;
			break;
		}

		/* get next */
		pRefCod = pRefCod->pNext;
	}

	if(nFindRef == 1)
	{
		/* adjust genome */
		if(strcmp(pRefCod->strGenome, "GAP") == 0)
		{
			return nGPos;
		}
		else
		{
			if(pRefCod->chGStrand == '+')
			{
				nGPos = nPos-pRefCod->nQStart+pRefCod->nGStart;
			}
			else if(pRefCod->chGStrand == '-')
			{
				nGPos = pRefCod->nGEnd-(nPos-pRefCod->nQStart);
			}
			else
			{
				nGPos = -1;
			}
		}
	}
	else
	{
		nGPos = -1;
	}

	/* return */
	return nGPos;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_GenMotifAln_Main: generate motif alignment from existing   */
/*   motif mapping results.                                                */
/* ----------------------------------------------------------------------- */ 
int MAlign_GenMotifAln_Main(char strParamPath[])
{
	/* define */

	/* for loading parameters */
	FILE *fpIn;
	char strLine[MED_LINE_LENGTH];
	char *chSep;
	int nError = 0;
	int ni;

	/* parameters */
	char strWorkPath[MED_LINE_LENGTH];
	char strFileHeader[MED_LINE_LENGTH];
	int nClusterNum = 0;
	int nSpeciesNum = 0;
	
	struct tagString **vSpeciesName = NULL;
	struct tagString **vMotifFile = NULL;
	
	/* ---------------------------------- */
	/* Step I: Load Parameters            */
	/* ---------------------------------- */
	fpIn = NULL;
	fpIn = fopen(strParamPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: MAlign_GenMotifAln_Main, cannot open the parameter file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strstr(strLine, "[Working Directory]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strWorkPath, chSep);
			AdjustDirectoryPath(strWorkPath);
		}

		else if(strstr(strLine, "[File Header]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			strcpy(strFileHeader, chSep);
			
			if(strcmp(strFileHeader, "") == 0)
			{
				printf("Error: file header not speciefied!\n");
				nError = 1;
				break;
			}
		}

		else if(strstr(strLine, "[Number of Clusters]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nClusterNum = atoi(chSep);
			}

			if(nClusterNum <= 0)
			{
				printf("Number of clusters = %d\n", nClusterNum);
				nError = 1;
				break;
			}
		}

		else if(strstr(strLine, "[Number of Species]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nSpeciesNum = atoi(chSep);
			}

			if(nSpeciesNum <= 0)
			{
				printf("Number of species = %d\n", nSpeciesNum);
				nError = 1;
				break;
			}

			vSpeciesName = (struct tagString **)calloc(nSpeciesNum, sizeof(struct tagString *));
			if(vSpeciesName == NULL)
			{
				printf("Error: MAlign_MotifMap_Main, cannot allocate memory for loading sequence alias!\n");
				nError = 1;
				break;
			}
		}

		else if(strstr(strLine, "[Species]") == strLine)
		{
			ni = 0;
			while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
			{
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] == '\0')
					continue;

				StringAddTail(vSpeciesName+ni, strLine);
			
				ni++;

				if(ni == nSpeciesNum)
					break;
			}
		}

		else if(strLine[0] == '#')
		{
		}

		else
		{
			printf("Error: MAlign_GenMotifAln_Main, unknown parameters!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	fclose(fpIn);

	if(nError > 0)
	{
		printf("Error: MAlign_GenMotifAln_Main, cannot load the parameter file correctly!\n");
		exit(EXIT_FAILURE);
	}

	/* ------------------------------------ */
	/* Step II: Map motifs to alignment     */
	/* ------------------------------------ */
	for(ni=0; ni<nClusterNum; ni++)
	{
		MAlign_GenerateMotifAln_GenerateAlignment(strWorkPath, strFileHeader, 
			ni, nSpeciesNum, vSpeciesName);
	}

	/* ------------------------------------ */
	/* Step III: Release memory             */
	/* ------------------------------------ */
	for(ni=0; ni<nSpeciesNum; ni++)
	{
		DeleteString(vSpeciesName[ni]);
	}
	free(vSpeciesName);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_GenerateMotifAln_GenerateAlignment: generate alignment          */
/*   including mapped sites.                                               */
/* ----------------------------------------------------------------------- */ 
int MAlign_GenerateMotifAln_GenerateAlignment(char strWorkPath[], char strFileHeader[], 
			int nClusterId, int nSpeciesNum, struct tagString **vSpeciesName)
{
	/* define */
	int nSeqNum;
	struct tagSequence *pAln;
	struct tagString **vMS;
	struct tagSequence *pSeq;
	char chBase;
	struct SEQLINKMAP **vSeqLinkMap;
	char strFileName[MED_LINE_LENGTH];
	int ni,nj,nLen;
	char *vS;
	int nAlnLen;
	int nIdn;

	/* load alignment */
	sprintf(strFileName, "%s%s_%d_mlagan.out", strWorkPath, strFileHeader, 
		nClusterId);
	pAln = NULL;
	nSeqNum = 0;
	nSeqNum = LoadFullSequenceList(strFileName, &pAln);
	if(nSeqNum != nSpeciesNum)
	{
		printf("Error: MAlign_GenerateMotifAln_GenerateAlignment, sequence number does not match!\n");
		exit(EXIT_FAILURE);
	}
	if(nSeqNum <= 1)
	{
		printf("Error: MAlign_GenerateMotifAln_GenerateAlignment, <=1 sequence available, require at least two sequences!\n");
		exit(EXIT_FAILURE);
	}

	/* init motif status */
	vMS = NULL;
	vMS = (struct tagString **)calloc(nSeqNum, sizeof(struct tagString*));
	if(vMS == NULL)
	{
		printf("Error: MAlign_GenerateMotifAln_GenerateAlignment, cannot allocate space for recording motif status!\n");
		exit(EXIT_FAILURE);
	}

	pSeq = pAln;
	for(ni=0; ni<nSeqNum; ni++)
	{
		nLen = pSeq->m_nLength;
		vMS[ni] = NULL;
		vMS[ni] = CreateString(nLen);
		if(vMS[ni] == NULL)
		{
			printf("Error: MAlign_GenerateMotifAln_GenerateAlignment, cannot allocate space for recording motif status!\n");
			exit(EXIT_FAILURE);
		}

		vS = vMS[ni]->m_pString;
		for(nj=0; nj<nLen; nj++)
		{
			vS[nj] = ' ';
		}
		vS[nj] = '\0';

		pSeq = pSeq->m_pNext;
	}

	for(ni=0; ni<nSpeciesNum; ni++)
	{
		if(vMS[ni]->m_nLength != vMS[0]->m_nLength)
		{
			printf("Error: MAlign_GenerateMotifAln_GenerateAlignment, sequence length not the same!\n");
			exit(EXIT_FAILURE);
		}
	}
	nAlnLen = pAln->m_nLength;


	/* load seqlink map */
	vSeqLinkMap = NULL;
	vSeqLinkMap = (struct SEQLINKMAP **)calloc(nSeqNum, sizeof(struct SEQLINKMAP *));
	if(vSeqLinkMap == NULL)
	{
		printf("Error: MAlign_GenerateMotifAln_GenerateAlignment, cannot allocate space for recording sequence link status!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nSeqNum; ni++)
	{
		sprintf(strFileName, "%s%s_%d_%s.gln", strWorkPath, strFileHeader, 
				nClusterId, vSpeciesName[ni]->m_pString);
		MAlign_MotifMap_LoadSeqLinkMap(vSeqLinkMap+ni, strFileName);
	}

	
	/* find identity */
	for(ni=0; ni<nAlnLen; ni++)
	{
		nIdn = 1;
		pSeq = pAln;
		chBase = pSeq->m_pSequence->m_pString[ni];
		if(chBase == '-')
		{
			nIdn = 0;
		}
		else
		{
			for(nj=1; nj<nSpeciesNum; nj++)
			{
				pSeq = pSeq->m_pNext;
				if(pSeq->m_pSequence->m_pString[ni] == chBase)
				{
					vMS[nj]->m_pString[ni] = '.';
				}
				else
				{
					nIdn = 0;
				}
			}
		}

		if(nIdn == 1)
		{
			for(nj=1; nj<nSpeciesNum; nj++)
			{
				vMS[nj]->m_pString[ni] = ':';
			}
		}
	}


	/* update motif status */
	pSeq = pAln;
	for(ni=0; ni<nSpeciesNum; ni++)
	{
		MAlign_GenerateMotifAln_UpdateAlignAndMotif(pSeq, vMS[ni], strWorkPath, strFileHeader, 
			nClusterId, vSpeciesName[ni]->m_pString);
		pSeq = pSeq->m_pNext;
	}

	/* export */
	sprintf(strFileName, "%s%s_%d_motif.aln", strWorkPath, strFileHeader, nClusterId);
	MAlign_GenerateMotifAln_ExportAlignment(nSpeciesNum, vSpeciesName, pAln, vMS, vSeqLinkMap,
		strFileName);

	/* clear memory */
	for(ni=0; ni<nSeqNum; ni++)
	{
		DeleteString(vMS[ni]);
		vMS[ni] = NULL;
		SEQLINKMAPCLEARLIST(vSeqLinkMap+ni);
	}
	free(vMS);
	free(vSeqLinkMap);
	SequenceListClear(&pAln);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_GenerateMotifAln_UpdateAlignAndMotif: update alignment to       */
/*   reflect motif sites and repeats.                                      */
/* ----------------------------------------------------------------------- */ 
int MAlign_GenerateMotifAln_UpdateAlignAndMotif(struct tagSequence *pSeq, 
			struct tagString *pMS, char strWorkPath[], char strFileHeader[], 
			int nClusterId, char strSpeciesName[])
{
	/* define */
	FILE *fpIn;
	struct tagSequence *pOriginSeq;
	struct tagString *pOriginStatus;
	struct DOUBLEMATRIX *pOriginScore;
	struct MOTIFSITELINKMAP *pSiteList = NULL;
	struct MOTIFSITELINKMAP *pSite;
	char strLine[MED_LINE_LENGTH];
	char strFileName[MED_LINE_LENGTH];
	int nCount;
	char *chA,*chS;
	char *chOA,*chOS;
	double *pScore;
	int ni,nAlnLen,nPos;

	/* load sequences */
	sprintf(strFileName, "%s%s_%d_%s.fa", strWorkPath, strFileHeader, 
		nClusterId, strSpeciesName);
	pOriginSeq = NULL;
	nCount = LoadFullSequenceList(strFileName, &pOriginSeq);
	if( (nCount != 1) || (pOriginSeq == NULL))
	{
		printf("Error: MAlign_GenerateMotifAln_UpdateAlignAndMotif, wrong number of sequence loaded!\n");
		exit(EXIT_FAILURE);
	}
	pOriginStatus = NULL;
	pOriginStatus = CreateString(pOriginSeq->m_nLength);
	pOriginScore = NULL;
	pOriginScore = CreateDoubleMatrix(1, pOriginSeq->m_nLength);
	if( (pOriginStatus == NULL) || (pOriginScore == NULL) )
	{
		printf("Error: MAlign_GenerateMotifAln_UpdateAlignAndMotif, cannot create memory for origin status&score!\n");
		exit(EXIT_FAILURE);
	}
	chS = pOriginStatus->m_pString;
	pScore = pOriginScore->pMatElement;
	for(ni=0; ni<pOriginSeq->m_nLength; ni++)
	{
		chS[ni] = ' ';
		pScore[ni] = -DM_ACCESS_VIOLATION;
	}

	/* load sites */
	sprintf(strFileName, "%s%s_%d_%s.map", strWorkPath, strFileHeader, 
		nClusterId, strSpeciesName);
	fpIn = NULL;
	fpIn = fopen(strFileName, "r");
	if(fpIn == NULL)
	{
		printf("Error: MAlign_GenerateMotifAln_UpdateAlignAndMotif, cannot open site coordinates file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		pSite = NULL;
		pSite = MOTIFSITELINKMAPCREATE();
		if(pSite == NULL)
		{
			printf("Error: MAlign_GenerateMotifAln_UpdateAlignAndMotif, cannot create site link map!\n");
			exit(EXIT_FAILURE);
		}

		sscanf(strLine, "%s %s %d %d %c %lf %s %s %d %d %c %s %d %d %c",
			pSite->strMotifSymbol, pSite->strMotifName, 
			&(pSite->nQStart), &(pSite->nQEnd), &(pSite->chQStrand),
			&(pSite->dScore), &(pSite->strSiteSeq),
			pSite->strHit, &(pSite->nHStart), &(pSite->nHEnd), &(pSite->chHStrand),
			pSite->strGenome, &(pSite->nGStart), &(pSite->nGEnd), &(pSite->chGStrand));

		strcpy(pSite->strQuery, strSpeciesName);

		MOTIFSITELINKMAP_INSERTBYQUERYLOC(pSite, &pSiteList);
	}

	fclose(fpIn);

	/* update status */
	while(pSiteList != NULL)
	{
		pSite = pSiteList;
		pSiteList = pSite->pNext;
		pSite->pNext = NULL;
		
		for(ni=pSite->nQStart; ni<=pSite->nQEnd; ni++)
		{
			if(pSite->dScore > pOriginScore->pMatElement[ni])
			{
				pOriginStatus->m_pString[ni] = pSite->strMotifSymbol[0];
				pOriginScore->pMatElement[ni] = pSite->dScore;
			}
		}
		
		MOTIFSITELINKMAPDESTROY(pSite);
	}


	/* generate aln status */
	nAlnLen = pSeq->m_nLength;
	chA = pSeq->m_pSequence->m_pString;
	chS = pMS->m_pString;
	chOA = pOriginSeq->m_pSequence->m_pString;
	chOS = pOriginStatus->m_pString;
	nPos = 0;
	for(ni=0; ni<nAlnLen; ni++)
	{
		if(chA[ni] != '-')
		{
			chA[ni] = chOA[nPos];
			if(chOS[nPos] != ' ')
				chS[ni] = chOS[nPos];
			nPos++;
		}
	}

	if(pSiteList != NULL)
	{
		printf("Error: MAlign_GenerateMotifAln_UpdateAlignAndMotif, not all sites are added!\n");
		exit(EXIT_FAILURE);
	}

	/* release memory */
	SequenceDelete(pOriginSeq);
	DeleteString(pOriginStatus);
	DestroyDoubleMatrix(pOriginScore);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_GenerateMotifAln_ExportAlignment: export alignment.             */
/* ----------------------------------------------------------------------- */ 
int MAlign_GenerateMotifAln_ExportAlignment(int nSpeciesNum, 
		struct tagString **vSpeciesName, struct tagSequence *pAln, 
		struct tagString **vMS, struct SEQLINKMAP **vSeqLinkMap,
		char strOutFile[])
{
	/* define */
	FILE *fpOut;
	struct INTMATRIX *pPos;
	struct tagSequence *pCurrentAln;
	char *vSeq,*vSeq2;
	int ni,nj;
	int nLineLen;
	int nAlnLen;
	int nRemainLen;
	char strLine[LINE_LENGTH];
	int nRefStart;
	int nPos;

	/* init */
	if(nSpeciesNum <= 0)
	{
		printf("Warning: MAlign_GenerateMotifAln_ExportAlignment, no alignment sequence!\n");
		return PROC_SUCCESS;
	}

	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: MAlign_GenerateMotifAln_ExportAlignment, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	pPos = NULL;
	pPos = CreateIntMatrix(1, nSpeciesNum);
	if(pPos == NULL)
	{
		printf("Error: MAlign_GenerateMotifAln_ExportAlignment, cannot track alignment position!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nSpeciesNum; ni++)
	{
		if(vMS[ni]->m_nLength != vMS[0]->m_nLength)
		{
			printf("Error: MAlign_GenerateMotifAln_ExportAlignment, sequence length not the same!\n");
			exit(EXIT_FAILURE);
		}
	}
	nAlnLen = pAln->m_nLength;

	/* export */
	nRefStart = 0;
	nRemainLen = nAlnLen;
	while(nRemainLen > 0)
	{
		nLineLen = ALIGN_LINE_LEN;
		if(nLineLen > nRemainLen)
			nLineLen = nRemainLen;

		pCurrentAln = pAln;
		for(ni=0; ni<nSpeciesNum; ni++)
		{
			/* export motif site */
			if( (ni == 0) && (vSeqLinkMap[0] != NULL) )
			{
				strcpy(strLine, vSeqLinkMap[0]->strGenome);
				MAlign_MotifMap_TRUNCSTRING(strLine, ALIGN_NAME_LEN);
				fprintf(fpOut, "%s ", strLine);

				nPos = pPos->pMatElement[ni];
				nPos = MAlign_MotifMap_GetGenomicCod(nPos, vSeqLinkMap[0]); 
				sprintf(strLine, "%d", nPos);
				MAlign_MotifMap_TRUNCSTRING(strLine, ALIGN_START_LEN);
				fprintf(fpOut, "%s ", strLine);
			}
			else
			{
				strcpy(strLine, "");
				MAlign_MotifMap_TRUNCSTRING(strLine, ALIGN_NAME_LEN);
				fprintf(fpOut, "%s ", strLine);

				sprintf(strLine, "");
				MAlign_MotifMap_TRUNCSTRING(strLine, ALIGN_START_LEN);
				fprintf(fpOut, "%s ", strLine);
			}

			vSeq = vMS[ni]->m_pString+nRefStart;
			for(nj=0; nj<nLineLen; nj++)
			{
				fprintf(fpOut, "%c", vSeq[nj]);
			}

			fprintf(fpOut, " ");

			if( (ni == 0) && (vSeqLinkMap[0] != NULL) )
			{
				nPos = pPos->pMatElement[ni];
				vSeq2 = pCurrentAln->m_pSequence->m_pString+nRefStart;
				for(nj=0; nj<nLineLen; nj++)
				{
					if(vSeq2[nj] != '-')
						nPos += 1;
				}
				nPos -= 1;
				if(nPos < pPos->pMatElement[ni])
					nPos = pPos->pMatElement[ni];

				nPos = MAlign_MotifMap_GetGenomicCod(nPos, vSeqLinkMap[0]); 
				sprintf(strLine, "%d", nPos);
				MAlign_MotifMap_TRUNCSTRING(strLine, ALIGN_START_LEN);
				fprintf(fpOut, "%s\n", strLine);
			}
			else
			{
				sprintf(strLine, "");
				MAlign_MotifMap_TRUNCSTRING(strLine, ALIGN_START_LEN);
				fprintf(fpOut, "%s\n", strLine);
			}

			/* export sequence name */
			strcpy(strLine, vSpeciesName[ni]->m_pString);
			MAlign_MotifMap_TRUNCSTRING(strLine, ALIGN_NAME_LEN);
			fprintf(fpOut, "%s ", strLine);

			/* export start site */
			sprintf(strLine, "%d", pPos->pMatElement[ni]);
			MAlign_MotifMap_TRUNCSTRING(strLine, ALIGN_START_LEN);
			fprintf(fpOut, "%s ", strLine);

			/* export alignment */
			vSeq = pCurrentAln->m_pSequence->m_pString+nRefStart;
			for(nj=0; nj<nLineLen; nj++)
			{
				fprintf(fpOut, "%c", vSeq[nj]);
				if(vSeq[nj] != '-')
					pPos->pMatElement[ni] += 1;
			}

			/* export end site */
			fprintf(fpOut, " ");

			sprintf(strLine, "%d", pPos->pMatElement[ni]);
			MAlign_MotifMap_TRUNCSTRING(strLine, ALIGN_END_LEN);
			fprintf(fpOut, "%s\n", strLine);

			/* get next */
			pCurrentAln = pCurrentAln->m_pNext;
		}

		fprintf(fpOut, "\n\n");
		nRemainLen -= nLineLen;
		nRefStart += nLineLen;
	}

	/* close file */
	fclose(fpOut);
	DestroyIntMatrix(pPos);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_Genome_PrepareOrtholog: prepare ortholog region for malign      */
/* ----------------------------------------------------------------------- */
int MAlign_Genome_PrepareOrtholog(char strInPath[], char strOutPath[],
			int nSpeciesNum, int nSkip, int nRefType, int nUp, int nDown)
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	char strLine[LONG_LINE_LENGTH];
	char strRefLine[LONG_LINE_LENGTH];
	int ni;
	char strSpecies[LINE_LENGTH];
	int nMatchCount;
	int nOrienCount;
	struct tagRefGene *pRefGene;
	char *chp1,*chp2;
	int nStart,nEnd;

	/* open */
	fpIn = NULL;
	fpIn = fopen(strInPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: MAlign_Genome_PrepareOrtholog, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}
	fpOut = NULL;
	fpOut = fopen(strOutPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: MAlign_Genome_PrepareOrtholog, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	/* load */
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strLine[0] == '>')
		{
			fprintf(fpOut, "%s\n", strLine);

			/* load first line */
			fgets(strLine, LONG_LINE_LENGTH, fpIn);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "---") == strLine)
			{
				sprintf(strRefLine, "NA\t-1\t-1\t-1\t?");
			}
			else
			{
				chp1 = strchr(strLine, '\t');
				chp1++;
				sprintf(strRefLine, "%s", chp1);
			}

			/* load other lines */
			for(ni=0; ni<nSpeciesNum; ni++)
			{
				fgets(strLine, LONG_LINE_LENGTH, fpIn);
				StrTrimLeft(strLine);
				StrTrimRight(strLine);

				chp1 = strLine;
				chp2 = strchr(chp1, '\t');
				*chp2 = '\0';
				strcpy(strSpecies, chp1);
				chp1 = chp2+1;
				chp2 = strchr(chp1, '\t');
				*chp2 = '\0';
				nMatchCount = atoi(chp1);
				chp1 = chp2+1;
				chp2 = strchr(chp1, '\t');
				*chp2 = '\0';
				nOrienCount = atoi(chp1);

				if(ni == 0)
				{
					fprintf(fpOut, "%s\t%s\n", strSpecies, strRefLine);
				}
				
				if( (nSkip == 1) && (ni == 0) )
				{
					continue;
				}
			
				chp1 = chp2+1;
				if(strstr(chp1, "---") == chp1)
				{
					fprintf(fpOut, "%s\tNA\t-1\t-1\t-1\t?\n", strSpecies);  
				}
				else
				{
					pRefGene = NULL;
					pRefGene = RefGeneCreate();
					if(pRefGene == NULL)
					{
						printf("Error: MAlign_Genome_PrepareOrtholog, cannot read refgene correctly!\n");
						exit(EXIT_FAILURE);
					}

					RefGeneInit_FromGenomeLabFormat(pRefGene, chp1, strSpecies);
					if(nRefType == 0)
					{
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
					}
					else if(nRefType == 1)
					{
						if(pRefGene->chStrand == '-')
						{
							nStart = pRefGene->nTxEnd-nDown;
							nEnd = pRefGene->nTxEnd+nUp;
						}
						else
						{
							nStart = pRefGene->nTxStart-nUp;
							nEnd = pRefGene->nTxStart+nDown;
						}
					}
					else if(nRefType == 2)
					{
						if(pRefGene->chStrand == '-')
						{
							nStart = pRefGene->nTxStart-nDown;
							nEnd = pRefGene->nTxStart+nUp;
						}
						else
						{
							nStart = pRefGene->nTxEnd-nUp;
							nEnd = pRefGene->nTxEnd+nDown;
						}
					}
					else if(nRefType == 3)
					{
						if(pRefGene->chStrand == '-')
						{
							nStart = pRefGene->nCdsStart-nDown;
							nEnd = pRefGene->nCdsEnd+nUp;
						}
						else
						{
							nStart = pRefGene->nCdsStart-nUp;
							nEnd = pRefGene->nCdsEnd+nDown;
						}
					}
					else if(nRefType == 4)
					{
						if(pRefGene->chStrand == '-')
						{
							nStart = pRefGene->nCdsEnd-nDown;
							nEnd = pRefGene->nCdsEnd+nUp;
						}
						else
						{
							nStart = pRefGene->nCdsStart-nUp;
							nEnd = pRefGene->nCdsStart+nDown;
						}
					}
					else if(nRefType == 5)
					{
						if(pRefGene->chStrand == '-')
						{
							nStart = pRefGene->nCdsStart-nDown;
							nEnd = pRefGene->nCdsStart+nUp;
						}
						else
						{
							nStart = pRefGene->nCdsEnd-nUp;
							nEnd = pRefGene->nCdsEnd+nDown;
						}
					}
					else
					{
						printf("Error: MAlign_Genome_PrepareOrtholog, cannot incorrect -r parameter!\n");
						exit(EXIT_FAILURE);
					}

					fprintf(fpOut, "%s\t%s\t%d\t%d\t+\n", strSpecies, pRefGene->strChrom, nStart, nEnd);
					RefGeneDestroy(pRefGene);
				}
			}



			/* write last line */
			fprintf(fpOut, "\n");
		}
	}


	/* close files */
	fclose(fpIn);
	fclose(fpOut);


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_Genome_BlastHit_Main: MAlign_BlastHit main function             */
/*  Using blast to find orthologous segments of a set of sequences.        */
/* ----------------------------------------------------------------------- */ 
int MAlign_Genome_BlastHit_Main(char strParamPath[])
{
	/* define */

	/* for loading parameters */
	FILE *fpIn;
	FILE *fpCod;
	char strLine[MED_LINE_LENGTH];
	char strLine0[MED_LINE_LENGTH];
	char strRefLine[LONG_LINE_LENGTH];
	char *chSep,*chp;
	int nError = 0;
	int ni,nj;

	/* parameters */
	char strBlastPath[MED_LINE_LENGTH];
	double dEValueCut = 1e-5;
	int nSmallMask = 0;
	int nMinLenCut = 30;
	double dMinIdnCut = 0.6;
	int nColinear = 1;
	double dExtendPercent = 0.25;
	int nMaxExtension = 100;
	int nLinkGap = 100;
	int nKeepBest = 1;

	char strMlaganPath[MED_LINE_LENGTH];
	char strMlaganParam[MED_LINE_LENGTH];

	char strWorkPath[MED_LINE_LENGTH];
	char strOutFile[MED_LINE_LENGTH];
	char strInFile[MED_LINE_LENGTH];
	char strGenomeFile[MED_LINE_LENGTH];

	int nDatabaseNum = 0;
	struct tagString **vDataAlias = NULL;
	struct tagString **vDataDirec = NULL;
	char strDataAlias[MED_LINE_LENGTH];

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
	
	/* ---------------------------------- */
	/* Step I: Load Parameters            */
	/* ---------------------------------- */
	fpIn = NULL;
	fpIn = fopen(strParamPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: MAlign_Genome_BlastHit_Main, cannot open the parameter file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strstr(strLine, "[Path of Blast]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strBlastPath, chSep);
			AdjustDirectoryPath(strBlastPath);
		}

		else if(strstr(strLine, "[E-value Cutoff]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				dEValueCut = atof(chSep);
			}
		}

		else if(strstr(strLine, "[Mask Small Case]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nSmallMask = atoi(chSep);
			}
		}

		else if(strstr(strLine, "[Minimal Length]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nMinLenCut = atoi(chSep);
			}
		}

		else if(strstr(strLine, "[Minimal Identity]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				dMinIdnCut = atof(chSep);
			}
		}

		else if(strstr(strLine, "[Colinearity]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nColinear = atoi(chSep);
			}
		}

		else if(strstr(strLine, "[Extension Percentage]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				dExtendPercent = atof(chSep);
			}
		}

		else if(strstr(strLine, "[Maximal Extension]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nMaxExtension = atoi(chSep);
			}
		}

		else if(strstr(strLine, "[Maximal Link Gap]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nLinkGap = atoi(chSep);
			}
		}

		else if(strstr(strLine, "[Keep Best vs. All]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nKeepBest = atoi(chSep);
			}
		}

		else if(strstr(strLine, "[Path of MLAGAN]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strMlaganPath, chSep);
			AdjustDirectoryPath(strMlaganPath);
		}

		else if(strstr(strLine, "[MLAGAN Parameter]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strMlaganParam, chSep);
		}

		else if(strstr(strLine, "[Working Directory]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strWorkPath, chSep);
			AdjustDirectoryPath(strWorkPath);
		}

		else if(strstr(strLine, "[Export File]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			strcpy(strOutFile, chSep);
			
			if(strcmp(strOutFile, "") == 0)
			{
				printf("Error: output file not speciefied!\n");
				nError = 1;
				break;
			}
		}

		else if(strstr(strLine, "[Input File]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			strcpy(strInFile, chSep);
			
			if(strcmp(strInFile, "") == 0)
			{
				printf("Error: input file not speciefied!\n");
				nError = 1;
				break;
			}
		}

		else if(strstr(strLine, "[Genome Setting]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			strcpy(strGenomeFile, chSep);
			
			if(strcmp(strGenomeFile, "") == 0)
			{
				printf("Error: genome settings not speciefied!\n");
				nError = 1;
				break;
			}
		}

		else if(strstr(strLine, "[Number of Databases]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nDatabaseNum = atoi(chSep);
			}

			if(nDatabaseNum <= 0)
			{
				printf("Error: you have to speciefy at least 1 database!\n");
				nError = 1;
				break;
			}

			nDatabaseNum += 1;
			
			vDataAlias = (struct tagString **)calloc(nDatabaseNum, sizeof(struct tagString *));
			if(vDataAlias == NULL)
			{
				printf("Error: MAlign_Genome_BlastHit_Main, cannot allocate memory for loading sequence alias!\n");
				nError = 1;
				break;
			}

			vDataDirec = (struct tagString **)calloc(nDatabaseNum, sizeof(struct tagString *));
			if(vDataDirec == NULL)
			{
				printf("Error: MAlign_Genome_BlastHit_Main, cannot allocate memory for loading direction file path!\n");
				nError = 1;
				break;
			}
		}

		else if(strstr(strLine, "[Target & Databases]") == strLine)
		{
			ni = 0;
			while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
			{
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] == '\0')
					continue;

				if(ni >= nDatabaseNum)
				{
					printf("Error: MAlign_Genome_BlastHit_Main, database number does not match the database provided!\n");
					nError = 1;
					break;
				}

				sscanf(strLine, "%s", strDataAlias);
				StringAddTail(vDataAlias+ni, strDataAlias);
				
				ni++;
			}

			if( (ni != nDatabaseNum) || (nError == 1) ) 
			{
				printf("Error: MAlign_Genome_BlastHit_Main, database number does not match the database provided!\n");
				nError = 1;
				break;
			}
		}

		else if(strLine[0] == '#')
		{
		}

		else
		{
			printf("Error: MAlign_Genome_BlastHit_Main, unknown parameters!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	fclose(fpIn);

	if(nError > 0)
	{
		printf("Error: MAlign_Genome_BlastHit_Main, cannot load the parameter file correctly!\n");
		exit(EXIT_FAILURE);
	}

	/* ---------------------------------- */
	/* Step II: Load Genome Settings      */
	/* ---------------------------------- */
	sprintf(strLine, "%s%s", strWorkPath, strGenomeFile);
	fpIn = NULL;
	fpIn = fopen(strLine, "r");
	if(fpIn == NULL)
	{
		printf("Error: MAlign_Genome_BlastHit_Main, cannot load environment settings!\n");
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
			for(nj=0; nj<nDatabaseNum; nj++)
			{
				StringAddTail(vDataDirec+nj, chp);
			}		
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
	sprintf(strLine, "%s%s", strWorkPath, strInFile);
	fpIn = NULL;
	fpIn = fopen(strLine, "r");
	if(fpIn == NULL)
	{
		printf("Error: MAlign_Genome_BlastHit_Main, cannot open input file!\n");
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
			sprintf(strLine, "%s%s_%d_cod.tmp", strWorkPath, strOutFile, nSeqNum);
			fpCod = NULL;
			fpCod = fopen(strLine, "w");
			if(fpCod == NULL)
			{
				printf("Error: cannot open temporary cod file!\n");
				exit(EXIT_FAILURE);
			}
			for(ni=0; ni<nDatabaseNum; ni++)
			{
				fgets(strRefLine, LONG_LINE_LENGTH, fpIn);
				StrTrimLeft(strRefLine);
				StrTrimRight(strRefLine);
				MAlign_Genome_PrepareIndividualFasta(strRefLine, strWorkPath, strOutFile, ni, 
					nSeqNum, vDataAlias, vDataDirec, nSpeciesNum, vSpeciesName, vSpeciesSeq, 
					vSpeciesCons, vSpeciesAnnot, vChrLen, fpCod);
			}
			fclose(fpCod);

			/* ---------------------------------- */
			/* Step III.2: Blast hit              */
			/* ---------------------------------- */
			printf("Preprocessing cluster %d...\n", nSeqNum);
			MAlign_BlastHit_PrepareOrthologSegments(strWorkPath, nDatabaseNum,
				vDataAlias, strOutFile, nSeqNum,
				strBlastPath, dEValueCut, nSmallMask, nMinLenCut, dMinIdnCut,
				nColinear, dExtendPercent, nMaxExtension, nLinkGap, nKeepBest);
		
			/* ---------------------------------- */
			/* Step III.3: Remove Temporary file  */
			/* ---------------------------------- */
			sprintf(strLine, "%s%s*.tmp*", strWorkPath, strOutFile);
			RemoveFiles(strLine);
			nSeqNum++;
		}
	}

	fclose(fpIn);

	/* ---------------------------------- */
	/* Step IV: MLAGAN                    */
	/* ---------------------------------- */
	for(ni=0; ni<nSeqNum; ni++)
	{
		printf("Align cluster %d...\n", ni);
		MAlign_MLAGAN_Align(strWorkPath, nDatabaseNum,
			vDataAlias, strOutFile, ni,
			strMlaganPath, strMlaganParam);
	}

	/* ---------------------------------- */
	/* Step V: Remove Files               */
	/* ---------------------------------- */
	sprintf(strLine, "%s%s*.tmp*", strWorkPath, strOutFile);
	RemoveFiles(strLine);
	sprintf(strLine, "%s*.anchors", strWorkPath);
	RemoveFiles(strLine);

	/* ---------------------------------- */
	/* Step VI: Release memory            */
	/* ---------------------------------- */
	for(ni=0; ni<nDatabaseNum; ni++)
	{
		DeleteString(vDataAlias[ni]);
		DeleteString(vDataDirec[ni]);
	}
	free(vDataAlias);
	free(vDataDirec);

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

/* ----------------------------------------------------------------------- */ 
/*  MAlign_Genome_PrepareIndividualFasta: Prepare individual fasta file.   */
/* ----------------------------------------------------------------------- */ 
int MAlign_Genome_PrepareIndividualFasta(char strRefLine[], char strWorkPath[], char strOutFile[],
			int nDataId, int nClusterId, struct tagString **vDataAlias, struct tagString **vDataDirec, 
			int nSpeciesNum, struct tagString **vSpeciesName, 
			struct tagString **vSpeciesSeq, struct tagString **vSpeciesCons, 
			struct tagString **vSpeciesAnnot, struct INTMATRIX **vChrLen, 
			FILE *fpCod)
{
	/* define */
	FILE *fpSeq;
	char strFilePath[MED_LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	char strChr[LINE_LENGTH];
	char strTemp[LINE_LENGTH];
	int nChr,nStart,nEnd;
	char chStrand;
	int nSpeciesId;
	int nFind;
	struct tagSequence *pSeq;

	/* init */
	sprintf(strFilePath, "%s%s_%s_%d_fa.tmp", strWorkPath, strOutFile, vDataAlias[nDataId]->m_pString, nClusterId);
	fpSeq = NULL;
	fpSeq = fopen(strFilePath, "w");
	if(fpSeq == NULL)
	{
		printf("Error: MAlign_Genome_PrepareIndividualFasta, cannot open fasta output file!\n");
		exit(EXIT_FAILURE);
	}
	
	/* write */
	if(strstr(strRefLine, "---") == strRefLine)
	{
		fprintf(fpSeq, ">%s\n", vDataAlias[nDataId]->m_pString);
		fprintf(fpSeq, "N\n");
		fprintf(fpCod, "%s\tNA\t-1\t-1\t-1\t?\n", vDataAlias[nDataId]->m_pString);
		fclose(fpSeq);
		return PROC_SUCCESS;
	}
	
	sscanf(strRefLine, "%s %s %d %d %c", strSpecies, strChr, &nStart, &nEnd, &chStrand);
	nChr = Genome_ChromosomeName_To_Index(strChr, strSpecies);
	
	if( (strcmp(strChr, "NA") == 0) || (strcmp(strChr, "---") == 0) || (nChr<0) || (chStrand == '?') )
	{
		fprintf(fpSeq, ">%s\n", vDataAlias[nDataId]->m_pString);
		fprintf(fpSeq, "N\n");
		fprintf(fpCod, "%s\tNA\t-1\t-1\t-1\t?\n", vDataAlias[nDataId]->m_pString);
		fclose(fpSeq);
		return PROC_SUCCESS;
	}

	nFind = 0;
	for(nSpeciesId=0; nSpeciesId<nSpeciesNum; nSpeciesId++)
	{
		if(strcmp(strSpecies, vSpeciesName[nSpeciesId]->m_pString) == 0)
		{
			nFind = 1;
			break;
		}
	}

	if(nFind == 0)
	{
		fprintf(fpSeq, ">%s\n", vDataAlias[nDataId]->m_pString);
		fprintf(fpSeq, "N\n");
		fprintf(fpCod, "%s\tNA\t-1\t-1\t-1\t?\n", vDataAlias[nDataId]->m_pString);
		fclose(fpSeq);
		return PROC_SUCCESS;
	}

	if(nStart < 0)
		nStart = 0;
	if(nEnd >= vChrLen[nSpeciesId]->pMatElement[nChr-1])
		nEnd = vChrLen[nSpeciesId]->pMatElement[nChr-1]-1;

	sprintf(strFilePath, "%s%s.sq", vSpeciesSeq[nSpeciesId]->m_pString, strChr);
	pSeq = NULL;
	pSeq = Genome_Code_4bit_GetSeq(strFilePath, nStart, nEnd);
	if(pSeq == NULL)
	{
		printf("Error: cannot load sequence!\n");
		exit(EXIT_FAILURE);
	}
	
	pSeq->m_nIndex = nClusterId;
	strcpy(pSeq->m_strAlias, vDataAlias[nDataId]->m_pString);
	

	strcpy(strTemp, vDataDirec[nDataId]->m_pString);
	StrMakeUpper(strTemp);
	if(strcmp(strTemp, "GENEWISE") == 0)
	{
		SequenceWriteToFasta_ByStrand(pSeq, fpSeq, chStrand, 1);
		fprintf(fpCod, "%s\t%s\t%s\t%d\t%d\t%c\n", vDataAlias[nDataId]->m_pString, strSpecies, strChr, nStart, nEnd, chStrand);
	}
	else
	{
		SequenceWriteToFasta_ByStrand(pSeq, fpSeq, '+', 1);
		fprintf(fpCod, "%s\t%s\t%s\t%d\t%d\t+\n", vDataAlias[nDataId]->m_pString, strSpecies, strChr, nStart, nEnd);
	}

	/* close files */
	fclose(fpSeq);
	SequenceDelete(pSeq);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_ModuleMap_Main: MAlign_ModuleMap main function                  */
/*  Search for modules to aligned orthologous seqeunces.                   */
/* ----------------------------------------------------------------------- */ 
int MAlign_ModuleMap_Main(char strParamPath[])
{
	/* define */

	/* for loading parameters */
	FILE *fpIn;
	FILE *fpOut;
	char strLine[MED_LINE_LENGTH];
	char strFileName[MED_LINE_LENGTH];
	char *chSep,*chSep2;
	int nError = 0;
	int ni,nj;

	/* parameters */
	char strWorkPath[MED_LINE_LENGTH];
	char strFileHeader[MED_LINE_LENGTH];
	char strOutputFile[LINE_LENGTH];
	int nClusterNum = 0;
	int nSpeciesNum = 0;
	int nMotifNum = 0;
	int nModuleLen = 100;
	struct tagString **vSpeciesName = NULL;
	struct tagString **vMotifSymbol = NULL;
	struct tagString **vMotifName = NULL;
	struct INTMATRIX *pMotifC = NULL;
	struct INTMATRIX *pModuleC = NULL;
	struct BYTEMATRIX *pModuleM = NULL;
	char strSpecies[LINE_LENGTH];
	char strTemp[LINE_LENGTH];
	int nTemp;

	/* ---------------------------------- */
	/* Step I: Load Parameters            */
	/* ---------------------------------- */
	fpIn = NULL;
	fpIn = fopen(strParamPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: MAlign_ModuleMap_Main, cannot open the parameter file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strstr(strLine, "[Working Directory]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strWorkPath, chSep);
			AdjustDirectoryPath(strWorkPath);
		}

		else if(strstr(strLine, "[File Header]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			strcpy(strFileHeader, chSep);
			
			if(strcmp(strFileHeader, "") == 0)
			{
				printf("Error: file header not speciefied!\n");
				nError = 1;
				break;
			}
		}

		else if(strstr(strLine, "[Export File]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			strcpy(strOutputFile, chSep);
			
			if(strcmp(strOutputFile, "") == 0)
			{
				printf("Error: output file not speciefied!\n");
				nError = 1;
				break;
			}
		}

		else if(strstr(strLine, "[Number of Clusters]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nClusterNum = atoi(chSep);
			}

			if(nClusterNum <= 0)
			{
				printf("Number of clusters = %d\n", nClusterNum);
				nError = 1;
				break;
			}
		}

		else if(strstr(strLine, "[Number of Species]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nSpeciesNum = atoi(chSep);
			}

			if(nSpeciesNum <= 0)
			{
				printf("Number of species = %d\n", nSpeciesNum);
				nError = 1;
				break;
			}

			vSpeciesName = (struct tagString **)calloc(nSpeciesNum, sizeof(struct tagString *));
			if(vSpeciesName == NULL)
			{
				printf("Error: MAlign_ModuleMap_Main, cannot allocate memory for loading sequence alias!\n");
				nError = 1;
				break;
			}
		}

		else if(strstr(strLine, "[Number of Motifs]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nMotifNum = atoi(chSep);
			}

			if(nMotifNum <= 0)
			{
				printf("Number of motifs = %d\n", nMotifNum);
				nError = 1;
				break;
			}

			vMotifSymbol = (struct tagString **)calloc(nMotifNum, sizeof(struct tagString *));
			if(vMotifSymbol == NULL)
			{
				printf("Error: MAlign_ModuleMap_Main, cannot allocate memory for loading motif mapping parameters!\n");
				nError = 1;
				break;
			}

			vMotifName = (struct tagString **)calloc(nMotifNum, sizeof(struct tagString *));
			if(vMotifName == NULL)
			{
				printf("Error: MAlign_ModuleMap_Main, cannot allocate memory for loading motif mapping parameters!\n");
				nError = 1;
				break;
			}
		}

		else if(strstr(strLine, "[Motif Criteria]") == strLine)
		{
			/* load species */
			while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
			{
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] != '\0')
					break;
			}

			chSep2 = strchr(strLine, '\t');
			chSep = chSep2+1;
			chSep2 = strchr(chSep, '\t');

			chSep = chSep2+1;
			for(ni=0; ni<nSpeciesNum; ni++)
			{
				chSep2 = strchr(chSep, '\t');
				if(chSep2 == NULL)
				{
					printf("Error: MAlign_ModuleMap_Main, species number not correctly specified!\n");
					exit(EXIT_FAILURE);
				}
				*chSep2 = '\0';
				strcpy(strSpecies, chSep);
				StrTrimLeft(strSpecies);
				StrTrimRight(strSpecies);
				StringAddTail(vSpeciesName+ni, strSpecies);
				chSep = chSep2+1;
			}

			/* load criteria */
			pMotifC = NULL;
			pMotifC = CreateIntMatrix(nMotifNum, (nSpeciesNum+1));
			if(pMotifC == NULL)
			{
				printf("Error: MAlign_ModuleMap_Main, cannot allocate memory for loading motif criteria!\n");
				nError = 1;
				break;
			}

			for(nj=0; nj<nMotifNum; nj++)
			{
				while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
				{
					StrTrimLeft(strLine);
					StrTrimRight(strLine);
					if(strLine[0] != '\0')
						break;
				}

				chSep2 = strchr(strLine, '\t');
				*chSep2 = '\0';
				strcpy(strTemp, strLine);
				StrTrimLeft(strTemp);
				StrTrimRight(strTemp);
				StringAddTail(vMotifSymbol+nj, strTemp);

				chSep = chSep2+1;
				chSep2 = strchr(chSep, '\t');
				*chSep2 = '\0';
				strcpy(strTemp, chSep);
				StrTrimLeft(strTemp);
				StrTrimRight(strTemp);
				StringAddTail(vMotifName+nj, strTemp);

				chSep = chSep2+1;
				for(ni=0; ni<nSpeciesNum; ni++)
				{
					chSep2 = strchr(chSep, '\t');
					if(chSep2 == NULL)
					{
						printf("Error: MAlign_ModuleMap_Main, species number not correctly specified!\n");
						exit(EXIT_FAILURE);
					}
					*chSep2 = '\0';
					strcpy(strTemp, chSep);
					StrTrimLeft(strTemp);
					StrTrimRight(strTemp);
					nTemp = atoi(strTemp);
					IMSETAT(pMotifC, nj, ni, nTemp);
					chSep = chSep2+1;
				}
				
				strcpy(strTemp, chSep);
				StrTrimLeft(strTemp);
				StrTrimRight(strTemp);
				nTemp = atoi(strTemp);
				IMSETAT(pMotifC, nj, ni, nTemp);
			}
		}
		
		else if(strstr(strLine, "[Module Length]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nModuleLen = atoi(chSep);
			}
		}

		else if(strstr(strLine, "[Module Criteria]") == strLine)
		{
			/* load conditions */
			while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
			{
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] != '\0')
					break;
			}

			/* load criteria */
			pModuleC = NULL;
			pModuleC = CreateIntMatrix(nMotifNum, 2);
			if(pModuleC == NULL)
			{
				printf("Error: MAlign_ModuleMap_Main, cannot allocate memory for loading module criteria!\n");
				nError = 1;
				break;
			}

			pModuleM = NULL;
			pModuleM = CreateByteMatrix(nMotifNum, 2);
			if(pModuleM == NULL)
			{
				printf("Error: MAlign_ModuleMap_Main, cannot allocate memory for loading module criteria!\n");
				nError = 1;
				break;
			}


			for(nj=0; nj<nMotifNum; nj++)
			{
				while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
				{
					StrTrimLeft(strLine);
					StrTrimRight(strLine);
					if(strLine[0] != '\0')
						break;
				}

				chSep2 = strchr(strLine, '\t');
				*chSep2 = '\0';
				strcpy(strTemp, strLine);
				StrTrimLeft(strTemp);
				StrTrimRight(strTemp);
				if(strcmp(strTemp, vMotifSymbol[nj]->m_pString) != 0)
				{
					printf("Error: MAlign_ModuleMap_Main, motif symbols in motif & module sections are not in the same order!\n");
					exit(EXIT_FAILURE);
				}
				
				chSep = chSep2+1;
				for(ni=0; ni<1; ni++)
				{
					chSep2 = strchr(chSep, '\t');
					if(chSep2 == NULL)
					{
						printf("Error: MAlign_ModuleMap_Main, conditions not correctly specified!\n");
						exit(EXIT_FAILURE);
					}
					*chSep2 = '\0';
					strcpy(strTemp, chSep);
					StrTrimLeft(strTemp);
					StrTrimRight(strTemp);
					StrMakeUpper(strTemp);
					if(strcmp(strTemp, "NA") == 0)
					{
						BMSETAT(pModuleM, nj, ni, 0);
					}
					else
					{
						nTemp = atoi(strTemp);
						IMSETAT(pModuleC, nj, ni, nTemp);
						BMSETAT(pModuleM, nj, ni, 1);
					}
					
					chSep = chSep2+1;
				}
				
				strcpy(strTemp, chSep);
				StrTrimLeft(strTemp);
				StrTrimRight(strTemp);
				StrMakeUpper(strTemp);
				if(strcmp(strTemp, "NA") == 0)
				{
					BMSETAT(pModuleM, nj, ni, 0);
				}
				else
				{
					nTemp = atoi(strTemp);
					IMSETAT(pModuleC, nj, ni, nTemp);
					BMSETAT(pModuleM, nj, ni, 1);
				}
			}
		}

		else if(strLine[0] == '#')
		{
		}

		else
		{
			printf("Error: MAlign_ModuleMap_Main, unknown parameters!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	fclose(fpIn);

	if(nError > 0)
	{
		printf("Error: MAlign_ModuleMap_Main, cannot load the parameter file correctly!\n");
		exit(EXIT_FAILURE);
	}

	/* ------------------------------------ */
	/* Step II: search for modules          */
	/* ------------------------------------ */
	sprintf(strFileName, "%s%s", strWorkPath, strOutputFile);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: MAlign_ModuleMap_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nClusterNum; ni++)
	{
		printf("Processing cluster %d...\n", ni);
		MAlign_ModuleMap_SearchModule(strWorkPath, strFileHeader, 
			ni, nSpeciesNum, vSpeciesName,  
			nMotifNum, vMotifSymbol, vMotifName,
			pMotifC, nModuleLen, pModuleC, pModuleM, fpOut);
	}

	fclose(fpOut);

	/* ------------------------------------ */
	/* Step III: Release memory             */
	/* ------------------------------------ */
	for(ni=0; ni<nSpeciesNum; ni++)
	{
		DeleteString(vSpeciesName[ni]);
	}
	free(vSpeciesName);

	for(ni=0; ni<nMotifNum; ni++)
	{
		DeleteString(vMotifSymbol[ni]);
		DeleteString(vMotifName[ni]);
	}
	free(vMotifSymbol);
	free(vMotifName);
	DestroyIntMatrix(pMotifC);
	DestroyIntMatrix(pModuleC);
	DestroyByteMatrix(pModuleM);
	
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_ModuleMap_SearchModule:                                         */
/*  Search for modules to aligned orthologous seqeunces.                   */
/* ----------------------------------------------------------------------- */ 
int MAlign_ModuleMap_SearchModule(char strWorkPath[], char strFileHeader[], 
			int nClusterId, int nSpeciesNum, struct tagString **vSpeciesName,  
			int nMotifNum, struct tagString **vMotifSymbol, struct tagString **vMotifName,
			struct INTMATRIX *pMotifC, int nModuleLen, 
			struct INTMATRIX *pModuleC, struct BYTEMATRIX *pModuleM, 
			FILE *fpOut)
{
	/* define */
	char strFileName[MED_LINE_LENGTH];
	int nSeqNum;
	struct tagSequence *pAln;
	struct tagSequence *pSeq;
	struct tagSequence **vAln;
	struct MOTIFSITELINKMAP **vMotifList;
	struct MOTIFSITELINKMAP *pConsMotifList;

	int ni;
	int nAlnLen;
	
	/* load alignment */
	sprintf(strFileName, "%s%s_%d_mlagan.out", strWorkPath, strFileHeader, 
		nClusterId);
	pAln = NULL;
	nSeqNum = 0;
	nSeqNum = LoadFullSequenceList(strFileName, &pAln);
	if(nSeqNum != nSpeciesNum)
	{
		printf("Error: MAlign_ModuleMap_SearchModule, sequence number does not match!\n");
		exit(EXIT_FAILURE);
	}
	if(nSeqNum <= 0)
	{
		printf("Error: MAlign_ModuleMap_SearchModule, no sequence available!\n");
		exit(EXIT_FAILURE);
	}

	vAln = NULL;
	vAln = (struct tagSequence **)calloc(nSpeciesNum, sizeof(struct tagSequence *));
	if(vAln == NULL)
	{
		printf("Error: MAlign_ModuleMap_SearchModule, cannot allocate memory for organizing alignment!\n");
		exit(EXIT_FAILURE);
	}

	nAlnLen = pAln->m_nLength;
	while(pAln != NULL)
	{
		pSeq = pAln;
		pAln = pSeq->m_pNext;
		pSeq->m_pNext = NULL;

		if(pSeq->m_nLength != nAlnLen)
		{
			printf("Error: MAlign_ModuleMap_SearchModule, alignment length do not match!\n");
			exit(EXIT_FAILURE);
		}

		for(ni=0; ni<nSpeciesNum; ni++)
		{
			if(strstr(pSeq->m_strAlias, vSpeciesName[ni]->m_pString) == pSeq->m_strAlias)
				break;
		}
		
		if(ni >= nSpeciesNum)
		{
			printf("Error: MAlign_ModuleMap_SearchModule, cannot find %s in alignment!\n", vSpeciesName[ni]->m_pString);
			exit(EXIT_FAILURE);
		}
		else
		{
			vAln[ni] = pSeq;
		}
	}

	/* load motif sites */
	vMotifList = NULL;
	vMotifList = (struct MOTIFSITELINKMAP **)calloc(nSpeciesNum, sizeof(struct MOTIFSITELINKMAP *));
	if(vMotifList == NULL)
	{
		printf("Error: MAlign_ModuleMap_SearchModule, cannot allocate memory for loading motifs!\n");
		exit(EXIT_FAILURE);
	}

	/* update motif status */
	for(ni=0; ni<nSpeciesNum; ni++)
	{
		/* load sites */
		sprintf(strFileName, "%s%s_%d_%s.map", strWorkPath, strFileHeader, 
			nClusterId, vSpeciesName[ni]->m_pString);
		vMotifList[ni] = MAlign_ModuleMap_LoadRawMotifMapping(strFileName, vSpeciesName[ni]->m_pString);
	
	}

	/* find a set of conserved motif sites */
	pConsMotifList = MAlign_ModuleMap_GetConservedMotif(nSpeciesNum, vSpeciesName, vAln, nAlnLen, 
			nMotifNum, vMotifSymbol, vMotifName, pMotifC, vMotifList);

	/* search for modules */
	MAlign_ModuleMap_ExportModules(pConsMotifList, nMotifNum, nModuleLen, pModuleC, pModuleM, 
		nClusterId, fpOut);
	
	/* clear memory */
	for(ni=0; ni<nSpeciesNum; ni++)
	{
		SequenceDelete(vAln[ni]);
		MOTIFSITELINKMAPCLEARLIST(vMotifList+ni);
	}
	free(vAln);
	free(vMotifList);
	MOTIFSITELINKMAPCLEARLIST(&pConsMotifList);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_ModuleMap_LoadRawMotifMapping:                                  */
/*  Load raw motif mapping.                                                */
/* ----------------------------------------------------------------------- */ 
struct MOTIFSITELINKMAP *MAlign_ModuleMap_LoadRawMotifMapping(char strInFile[], char strSpecies[])
{
	/* define */
	FILE *fpIn;
	char strLine[MED_LINE_LENGTH];
	struct MOTIFSITELINKMAP *pSiteList = NULL;
	struct MOTIFSITELINKMAP *pSite;

	/* open file */
	fpIn = NULL;
	fpIn = fopen(strInFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: MAlign_ModuleMap_LoadRawMotifMapping, cannot open site coordinates file!\n");
		exit(EXIT_FAILURE);
	}

	/* load */
	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		pSite = NULL;
		pSite = MOTIFSITELINKMAPCREATE();
		if(pSite == NULL)
		{
			printf("Error: MAlign_ModuleMap_LoadRawMotifMapping, cannot create site link map!\n");
			exit(EXIT_FAILURE);
		}

		sscanf(strLine, "%s %s %d %d %c %lf %s %s %d %d %c %s %d %d %c",
			pSite->strMotifSymbol, pSite->strMotifName, 
			&(pSite->nQStart), &(pSite->nQEnd), &(pSite->chQStrand),
			&(pSite->dScore), &(pSite->strSiteSeq),
			pSite->strHit, &(pSite->nHStart), &(pSite->nHEnd), &(pSite->chHStrand),
			pSite->strGenome, &(pSite->nGStart), &(pSite->nGEnd), &(pSite->chGStrand));

		strcpy(pSite->strQuery, strSpecies);

		MOTIFSITELINKMAP_INSERTBYQUERYLOC(pSite, &pSiteList);
	}

	/* close file */
	fclose(fpIn);

	/* return */
	return pSiteList;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_ModuleMap_GetConservedMotif:                                    */
/*  Get conserved motif.                                                   */
/* ----------------------------------------------------------------------- */ 
struct MOTIFSITELINKMAP *MAlign_ModuleMap_GetConservedMotif(int nSpeciesNum, 
			struct tagString **vSpeciesName, struct tagSequence **vAln, 
			int nAlnLen, int nMotifNum, struct tagString **vMotifSymbol, 
			struct tagString **vMotifName, struct INTMATRIX *pMotifC, 
			struct MOTIFSITELINKMAP **vMotifList)
{
	/* define */
	struct MOTIFSITELINKMAP *pConsSiteList = NULL;
	struct MOTIFSITELINKMAP *pSite;
	struct MOTIFSITELINKMAP *pStorage = NULL;
	int ni,nj,nk;
	struct INTMATRIX *pPos;
	struct INTMATRIX *pSiteStat;
	int nTemp;
	int nOK;
	int nTargetSite;

	/* init */
	pPos = NULL;
	pPos = CreateIntMatrix(1, nSpeciesNum);
	if(pPos == NULL)
	{
		printf("Error: MAlign_ModuleMap_GetConservedMotif, cannot track positions in the alignment!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nSpeciesNum; ni++)
	{
		pPos->pMatElement[ni] = -1;
	}

	/* scan */
	for(ni=0; ni<nAlnLen; ni++)
	{
		/* prepare space for checking motif */
		pSiteStat = NULL;
		pSiteStat = CreateIntMatrix(nMotifNum, (nSpeciesNum+1));
		if(pSiteStat == NULL)
		{
			printf("Error: MAlign_ModuleMap_GetConservedMotif, cannot track motif sites in the alignment!\n");
			exit(EXIT_FAILURE);
		}

		/* compute position */
		pStorage = NULL;
		for(nj=0; nj<nSpeciesNum; nj++)
		{
			if(vAln[nj]->m_pSequence->m_pString[ni] != '-')
				pPos->pMatElement[nj] += 1;

			while(vMotifList[nj] != NULL)
			{
				if(vMotifList[nj]->nQStart > pPos->pMatElement[nj])
					break;

				pSite = vMotifList[nj];
				vMotifList[nj] = pSite->pNext;
				pSite->pNext = NULL;

				if(pSite->nQStart == pPos->pMatElement[nj])
				{
					for(nk=0; nk<nMotifNum; nk++)
					{
						if(strcmp(pSite->strMotifSymbol, vMotifSymbol[nk]->m_pString) == 0)
							break;
					}
					if(nk < nMotifNum)
					{
						nTargetSite = 1;
						IMSETAT(pSiteStat, nk, nj, 1);
						nTemp = IMGETAT(pSiteStat, nk, nSpeciesNum)+1;
						IMSETAT(pSiteStat, nk, nSpeciesNum, nTemp);
					}
					else
					{
						nTargetSite = 0;
						MOTIFSITELINKMAPDESTROY(pSite);
					}

					if(nTargetSite == 1)
					{
						if(nj == 0)
						{
							MOTIFSITELINKMAP_INSERTBYQUERYLOC(pSite, &pStorage);
						}
						else
						{
							MOTIFSITELINKMAPDESTROY(pSite);
						}
					}
				}
				else
				{
					MOTIFSITELINKMAPDESTROY(pSite);
				}
			}
		}

		/* call conserved site */
		while(pStorage != NULL)
		{
			pSite = pStorage;
			pStorage = pSite->pNext;
			pSite->pNext = NULL;

			for(nk=0; nk<nMotifNum; nk++)
			{
				if(strcmp(pSite->strMotifSymbol, vMotifSymbol[nk]->m_pString) == 0)
					break;
			}

			if(nk < nMotifNum)
			{
				nTargetSite = 1;
				pSite->nSeqId = nk;
			}
			else
			{
				nTargetSite = 0;
				MOTIFSITELINKMAPDESTROY(pSite);
			}

			if(nTargetSite == 1)
			{
				nOK = 1;
				for(nj=0; nj<=nSpeciesNum; nj++)
				{
					if(IMGETAT(pSiteStat, nk, nj) < IMGETAT(pMotifC, nk, nj))
					{
						nOK = 0;
						break;
					}
				}

				if(nOK == 1)
				{
					MOTIFSITELINKMAP_INSERTBYQUERYLOC(pSite, &pConsSiteList);
				}
				else
				{
					MOTIFSITELINKMAPDESTROY(pSite);
				}
			}
		}

		/* release memory */
		DestroyIntMatrix(pSiteStat);
	}


	/* release memory */
	DestroyIntMatrix(pPos);

	/* return */
	return pConsSiteList;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_ModuleMap_ExportModules:                                        */
/*  Find and export modules.                                               */
/* ----------------------------------------------------------------------- */ 
int MAlign_ModuleMap_ExportModules(struct MOTIFSITELINKMAP *pConsMotifList, 
			int nMotifNum, int nModuleLen, struct INTMATRIX *pModuleC, 
			struct BYTEMATRIX *pModuleM, int nClusterId, FILE *fpOut)
{
	/* define */
	struct MOTIFSITELINKMAP *pS1,*pS2,*pSC;
	int nModuleStart = -1;
	int nModuleEnd = -1;
	struct INTMATRIX *pModuleStat;
	int ni,nMidPos;
	int nTrueModule;
	struct MOTIFSITELINKMAP *pModuleList = NULL;
	struct MOTIFSITELINKMAP *pNRList = NULL;
	struct MOTIFSITELINKMAP *pModule = NULL;
	struct MOTIFSITELINKMAP *pModulePrev = NULL;
	struct MOTIFSITELINKMAP *pNRPrev = NULL;

	/* check */
	if(pConsMotifList == NULL)
	{
		return PROC_SUCCESS;
	}

	/* init */
	pModuleStat = NULL;
	pModuleStat = CreateIntMatrix(1, nMotifNum);

	/* search for modules */
	pS1 = pConsMotifList;
	pS2 = pS1;
	nModuleStart = (pS1->nQStart+pS1->nQEnd)/2;
	nModuleEnd = nModuleStart;
	ni = pS1->nSeqId;
	pModuleStat->pMatElement[ni] += 1;

	/* construct the first module */
	pSC = pS1->pNext;
	while(pSC != NULL)
	{
		/* form module */
		nMidPos = (pSC->nQStart+pSC->nQEnd)/2;
		if( (nMidPos-nModuleStart) <= nModuleLen)
		{
			ni = pSC->nSeqId;
			pModuleStat->pMatElement[ni] += 1;
			nModuleEnd = nMidPos;
			pS2 = pSC;
			pSC = pSC->pNext;
		}
		else
		{
			break;
		}
	}

	nTrueModule = MAlign_ModuleMap_JudgeModule(nMotifNum, pModuleStat, pModuleC, pModuleM);
	if(nTrueModule == 1)
	{
		pModule = MOTIFSITELINKMAPCREATE();
		pModule->nQStart = pS1->nQStart;
		pModule->nQEnd = pS2->nQEnd;
		pModule->chQStrand = pS1->chQStrand;
		strcpy(pModule->strQuery, pS1->strQuery);
		pModule->nHStart = pS1->nHStart;
		pModule->nHEnd = pS2->nHEnd;
		pModule->chHStrand = pS1->chHStrand;
		strcpy(pModule->strHit, pS1->strHit);
		pModule->nGStart = pS1->nGStart;
		pModule->nGEnd = pS2->nGEnd;
		pModule->chGStrand = pS1->chGStrand;
		strcpy(pModule->strGenome, pS1->strGenome);

		if(pModuleList == NULL)
		{
			pModuleList = pModule;
			pModulePrev = pModule;
		}
		else
		{
			pModulePrev->pNext = pModule;
			pModulePrev = pModule;
		}
	}

	/* construct othre modules */
	while(pSC != NULL)
	{
		/* move old sites */
		nMidPos = (pSC->nQStart+pSC->nQEnd)/2;
		while( (nMidPos-nModuleStart) > nModuleLen )
		{
			ni = pS1->nSeqId;
			pModuleStat->pMatElement[ni] -= 1;

			pS1 = pS1->pNext;
			nModuleStart = (pS1->nQStart+pS1->nQEnd)/2;
		}

		/* update new sites */
		while(pSC != NULL)
		{
			/* form module */
			nMidPos = (pSC->nQStart+pSC->nQEnd)/2;
			if( (nMidPos-nModuleStart) <= nModuleLen )
			{
				ni = pSC->nSeqId;
				pModuleStat->pMatElement[ni] += 1;
				nModuleEnd = nMidPos;
				pS2 = pSC;
				pSC = pSC->pNext;
			}
			else
			{
				break;
			}
		}

		/* create new module */
		nTrueModule = MAlign_ModuleMap_JudgeModule(nMotifNum, pModuleStat, pModuleC, pModuleM);
		if(nTrueModule == 1)
		{
			pModule = MOTIFSITELINKMAPCREATE();
			pModule->nQStart = pS1->nQStart;
			pModule->nQEnd = pS2->nQEnd;
			pModule->chQStrand = pS1->chQStrand;
			strcpy(pModule->strQuery, pS1->strQuery);
			pModule->nHStart = pS1->nHStart;
			pModule->nHEnd = pS2->nHEnd;
			pModule->chHStrand = pS1->chHStrand;
			strcpy(pModule->strHit, pS1->strHit);
			pModule->nGStart = pS1->nGStart;
			pModule->nGEnd = pS2->nGEnd;
			pModule->chGStrand = pS1->chGStrand;
			strcpy(pModule->strGenome, pS1->strGenome);

			if(pModuleList == NULL)
			{
				pModuleList = pModule;
				pModulePrev = pModule;
			}
			else
			{
				pModulePrev->pNext = pModule;
				pModulePrev = pModule;
			}
		}
	}

	/* get nr modules */
	pModulePrev = NULL;
	while(pModuleList != NULL)
	{
		pModule = pModuleList;
		pModuleList = pModule->pNext;
		pModule->pNext = NULL;

		if(pNRList == NULL)
		{
			pNRList = pModule;
			pNRPrev = pModule;
		}
		else
		{
			if(pModule->nQStart <= pNRPrev->nQEnd)
			{
				pNRPrev->nQEnd = pModule->nQEnd;
				pNRPrev->nHEnd = pModule->nHEnd;
				pNRPrev->nGEnd = pModule->nGEnd;
				MOTIFSITELINKMAPDESTROY(pModule);
			}
			else
			{
				pNRPrev->pNext = pModule;
				pNRPrev = pModule;
			}
		}
	}

	/* export */
	while(pNRList != NULL)
	{
		pModule = pNRList;
		pNRList = pModule->pNext;
		pModule->pNext = NULL;

		fprintf(fpOut, "%d\t%s\t%d\t%d\t+\t%s\t%d\t%d\t+\n", nClusterId, 
			pModule->strQuery, pModule->nQStart, pModule->nQEnd,
			pModule->strGenome, pModule->nGStart, pModule->nGEnd);
		MOTIFSITELINKMAPDESTROY(pModule);
	}

	/* release memory */
	DestroyIntMatrix(pModuleStat);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_ModuleMap_JudgeModule:                                          */
/*  judge if a segment satisfy module condition.                           */
/* ----------------------------------------------------------------------- */ 
int MAlign_ModuleMap_JudgeModule(int nMotifNum, struct INTMATRIX *pModuleStat, 
			struct INTMATRIX *pModuleC, struct BYTEMATRIX *pModuleM)
{
	/* define */
	int nIsM = 1;
	int ni;

	/* judge */
	for(ni=0; ni<nMotifNum; ni++)
	{
		if(BMGETAT(pModuleM, ni, 0) == 1)
		{
			if(pModuleStat->pMatElement[ni] >= IMGETAT(pModuleC, ni, 0))
			{
			}
			else
			{
				nIsM = 0;
				break;
			}
		}

		if(BMGETAT(pModuleM, ni, 1) == 1)
		{
			if(pModuleStat->pMatElement[ni] <= IMGETAT(pModuleC, ni, 1))
			{
			}
			else
			{
				nIsM = 0;
				break;
			}
		}
	}

	/* return */
	return nIsM;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_KMerStat_Main: MAlign_KMerStat main function                    */
/*  Count frequency of kmers.                                              */
/* ----------------------------------------------------------------------- */ 
int MAlign_KMerStat_Main(char strParamPath[])
{
	/* define */

	/* for loading parameters */
	FILE *fpIn;
	FILE *fpOut;
	char strLine[MED_LINE_LENGTH];
	char strFileName[MED_LINE_LENGTH];
	char *chSep,*chSep2;
	int nError = 0;
	int ni,nj;

	/* parameters */
	char strWorkPath[MED_LINE_LENGTH];
	char strFileHeader[MED_LINE_LENGTH];
	char strOutputFile[LINE_LENGTH];
	int nClusterNum = 0;
	int nSpeciesNum = 0;
	int nKmerLen = 4;
	struct tagString **vSpeciesName = NULL;
	struct INTMATRIX *pMotifC = NULL;
	struct DOUBLEMATRIX *pKmerCount = NULL;
	char strSpecies[LINE_LENGTH];
	char strTemp[LINE_LENGTH];
	int nTemp;
	int nKmerNum;
	double dTotal = 0.0;
	int nWordId;
	char strWord[LINE_LENGTH];
	int nBaseTypeNum = 4;

	/* ---------------------------------- */
	/* Step I: Load Parameters            */
	/* ---------------------------------- */
	fpIn = NULL;
	fpIn = fopen(strParamPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: MAlign_KMerStat_Main, cannot open the parameter file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strstr(strLine, "[Working Directory]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			strcpy(strWorkPath, chSep);
			AdjustDirectoryPath(strWorkPath);
		}

		else if(strstr(strLine, "[File Header]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			strcpy(strFileHeader, chSep);
			
			if(strcmp(strFileHeader, "") == 0)
			{
				printf("Error: file header not speciefied!\n");
				nError = 1;
				break;
			}
		}

		else if(strstr(strLine, "[Export File]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			strcpy(strOutputFile, chSep);
			
			if(strcmp(strOutputFile, "") == 0)
			{
				printf("Error: output file not speciefied!\n");
				nError = 1;
				break;
			}
		}

		else if(strstr(strLine, "[Number of Clusters]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nClusterNum = atoi(chSep);
			}

			if(nClusterNum <= 0)
			{
				printf("Number of clusters = %d\n", nClusterNum);
				nError = 1;
				break;
			}
		}

		else if(strstr(strLine, "[Number of Species]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nSpeciesNum = atoi(chSep);
			}

			if(nSpeciesNum <= 0)
			{
				printf("Number of species = %d\n", nSpeciesNum);
				nError = 1;
				break;
			}

			vSpeciesName = (struct tagString **)calloc(nSpeciesNum, sizeof(struct tagString *));
			if(vSpeciesName == NULL)
			{
				printf("Error: MAlign_ModuleMap_Main, cannot allocate memory for loading sequence alias!\n");
				nError = 1;
				break;
			}
		}

		else if(strstr(strLine, "[k-mer Length]") == strLine)
		{
			chSep = strstr(strLine, "=");
			if(chSep == NULL)
			{
				printf("%s\n", strLine);
				nError = 1;
				break;
			}
			chSep++;
			StrTrimLeft(chSep);
			
			if(*chSep != '\0')
			{
				nKmerLen = atoi(chSep);
			}

			if(nKmerLen <= 0)
			{
				printf("kmer length = %d\n", nKmerLen);
				nError = 1;
				break;
			}
		}

		else if(strstr(strLine, "[Motif Criteria]") == strLine)
		{
			/* load species */
			while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
			{
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] != '\0')
					break;
			}

			chSep = strLine;
			for(ni=0; ni<nSpeciesNum; ni++)
			{
				chSep2 = strchr(chSep, '\t');
				if(chSep2 == NULL)
				{
					printf("Error: MAlign_KMerStat_Main, species number not correctly specified!\n");
					exit(EXIT_FAILURE);
				}
				*chSep2 = '\0';
				strcpy(strSpecies, chSep);
				StrTrimLeft(strSpecies);
				StrTrimRight(strSpecies);
				StringAddTail(vSpeciesName+ni, strSpecies);
				chSep = chSep2+1;
			}

			/* load criteria */
			pMotifC = NULL;
			pMotifC = CreateIntMatrix(1, (nSpeciesNum+1));
			if(pMotifC == NULL)
			{
				printf("Error: MAlign_KMerStat_Main, cannot allocate memory for loading motif criteria!\n");
				nError = 1;
				break;
			}

			while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
			{
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] != '\0')
					break;
			}

			chSep = strLine;
			for(ni=0; ni<nSpeciesNum; ni++)
			{
				chSep2 = strchr(chSep, '\t');
				if(chSep2 == NULL)
				{
					printf("Error: MAlign_KMerStat_Main, species number not correctly specified!\n");
					exit(EXIT_FAILURE);
				}
				*chSep2 = '\0';
				strcpy(strTemp, chSep);
				StrTrimLeft(strTemp);
				StrTrimRight(strTemp);
				nTemp = atoi(strTemp);
				IMSETAT(pMotifC, 0, ni, nTemp);
				chSep = chSep2+1;
			}
			
			strcpy(strTemp, chSep);
			StrTrimLeft(strTemp);
			StrTrimRight(strTemp);
			nTemp = atoi(strTemp);
			IMSETAT(pMotifC, 0, ni, nTemp);
		}
		
		else if(strLine[0] == '#')
		{
		}

		else
		{
			printf("Error: MAlign_KMerStat_Main, unknown parameters!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	fclose(fpIn);

	if(nError > 0)
	{
		printf("Error: MAlign_KMerStat_Main, cannot load the parameter file correctly!\n");
		exit(EXIT_FAILURE);
	}

	/* ------------------------------------ */
	/* Step II: count kmers                 */
	/* ------------------------------------ */
	nKmerNum = (int)(pow((double)nBaseTypeNum, (double)nKmerLen));
	pKmerCount = NULL;
	pKmerCount = CreateDoubleMatrix(1, nKmerNum);
	if(pKmerCount == NULL)
	{
		printf("Error: MAlign_KMerStat_Main, cannot create matrix for counting!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nClusterNum; ni++)
	{
		printf("Processing cluster %d...\n", ni);
		MAlign_KMerStat_Count(strWorkPath, strFileHeader, 
			ni, nSpeciesNum, vSpeciesName,  
			pMotifC, nKmerLen, pKmerCount);
	}

	/* ------------------------------------ */
	/* Step III: Export                     */
	/* ------------------------------------ */
	sprintf(strFileName, "%s%s", strWorkPath, strOutputFile);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: MAlign_KMerStat_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	dTotal = 0.0;
	for(ni=0; ni<nKmerNum; ni++)
	{
		dTotal += pKmerCount->pMatElement[ni];

		nWordId = ni;
		for(nj=1; nj<nKmerLen; nj++)
		{
			nTemp = nWordId%nBaseTypeNum;
			switch(nTemp)
			{
				case 0: strWord[nKmerLen-nj] = 'A';
					break;
				case 1: strWord[nKmerLen-nj] = 'C';
					break;
				case 2: strWord[nKmerLen-nj] = 'G';
					break;
				case 3: strWord[nKmerLen-nj] = 'T';
					break;
			}
			nWordId -= nTemp;
			nWordId /= nBaseTypeNum;
		}

		nTemp = nWordId;
		switch(nTemp)
		{
			case 0: strWord[nKmerLen-nj] = 'A';
				break;
			case 1: strWord[nKmerLen-nj] = 'C';
				break;
			case 2: strWord[nKmerLen-nj] = 'G';
				break;
			case 3: strWord[nKmerLen-nj] = 'T';
				break;
		}

		strWord[nKmerLen] = '\0';

		fprintf(fpOut, "%s\t%d\n", strWord, (int)(pKmerCount->pMatElement[ni]));
	}

	fprintf(fpOut, "Total\t%d\n", (int)dTotal);
	fclose(fpOut);

	/* ------------------------------------ */
	/* Step IV: Release memory              */
	/* ------------------------------------ */
	DestroyDoubleMatrix(pKmerCount);

	for(ni=0; ni<nSpeciesNum; ni++)
	{
		DeleteString(vSpeciesName[ni]);
	}
	free(vSpeciesName);
	DestroyIntMatrix(pMotifC);
	
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MAlign_KMerStat_Count:                                                 */
/*  Count kmers in aligned orthologous seqeunces.                          */
/* ----------------------------------------------------------------------- */ 
int MAlign_KMerStat_Count(char strWorkPath[], char strFileHeader[], 
			int nClusterId, int nSpeciesNum, struct tagString **vSpeciesName,  
			struct INTMATRIX *pMotifC, int nKmerLen, struct DOUBLEMATRIX *pKmerCount)
{
	/* define */
	char strFileName[MED_LINE_LENGTH];
	int nSeqNum;
	struct tagSequence *pAln;
	struct tagSequence *pSeq;
	struct tagSequence **vAln;
	struct INTMATRIX *pStatus;
	struct BYTEMATRIX *pCons;

	int ni,nj;
	int nAlnLen;
	int nWordId = 0;
	int nBadNum = 0;
	int nPass = 1;
	int nBaseTypeNum = 4;
	char chBase;
	int nScale = (int)pow((double)nBaseTypeNum, (double)(nKmerLen-1));
	
	/* load alignment */
	sprintf(strFileName, "%s%s_%d_mlagan.out", strWorkPath, strFileHeader, 
		nClusterId);
	pAln = NULL;
	nSeqNum = 0;
	nSeqNum = LoadFullSequenceList(strFileName, &pAln);
	if(nSeqNum != nSpeciesNum)
	{
		printf("Error: MAlign_KMerStat_Count, sequence number does not match!\n");
		exit(EXIT_FAILURE);
	}
	if(nSeqNum <= 0)
	{
		printf("Error: MAlign_KMerStat_Count, no sequence available!\n");
		exit(EXIT_FAILURE);
	}

	vAln = NULL;
	vAln = (struct tagSequence **)calloc(nSpeciesNum, sizeof(struct tagSequence *));
	if(vAln == NULL)
	{
		printf("Error: MAlign_KMerStat_Count, cannot allocate memory for organizing alignment!\n");
		exit(EXIT_FAILURE);
	}

	nAlnLen = pAln->m_nLength;
	while(pAln != NULL)
	{
		pSeq = pAln;
		pAln = pSeq->m_pNext;
		pSeq->m_pNext = NULL;

		if(pSeq->m_nLength != nAlnLen)
		{
			printf("Error: MAlign_KMerStat_Count, alignment length do not match!\n");
			exit(EXIT_FAILURE);
		}

		for(ni=0; ni<nSpeciesNum; ni++)
		{
			if(strstr(pSeq->m_strAlias, vSpeciesName[ni]->m_pString) == pSeq->m_strAlias)
				break;
		}
		
		if(ni >= nSpeciesNum)
		{
			printf("Error: MAlign_KMerStat_Count, cannot find %s in alignment!\n", vSpeciesName[ni]->m_pString);
			exit(EXIT_FAILURE);
		}
		else
		{
			vAln[ni] = pSeq;
			StrMakeUpper(pSeq->m_pSequence->m_pString);
		}
	}

	/* counting */
	if(nAlnLen < nKmerLen)
	{
		for(ni=0; ni<nSpeciesNum; ni++)
		{
			SequenceDelete(vAln[ni]);
		}
		free(vAln);

		/* return */
		return PROC_SUCCESS;
	}

	/* init conservation status */
	pStatus = NULL;
	pStatus = CreateIntMatrix(1, (nSpeciesNum+1));
	if(pStatus == NULL)
	{
		printf("Error: MAlign_KMerStat_Count, cannot trace status!\n");
		exit(EXIT_FAILURE);
	}
	pCons = NULL;
	pCons = CreateByteMatrix(1, nAlnLen);
	if(pCons == NULL)
	{
		printf("Error: MAlign_KMerStat_Count, cannot trace status!\n");
		exit(EXIT_FAILURE);
	}

	/* prepare the first word */
	nWordId = 0;
	for(ni=0; ni<(nKmerLen-1); ni++)
	{
		chBase = vAln[0]->m_pSequence->m_pString[ni];
		nPass = 1;
		switch(chBase)
		{
			case 'A': nWordId = nWordId*nBaseTypeNum;
				break;
			case 'a': nWordId = nWordId*nBaseTypeNum;
				break;
			case 'C': nWordId = nWordId*nBaseTypeNum+1;
				break;
			case 'c': nWordId = nWordId*nBaseTypeNum+1;
				break;
			case 'G': nWordId = nWordId*nBaseTypeNum+2;
				break;
			case 'g': nWordId = nWordId*nBaseTypeNum+2;
				break;
			case 'T': nWordId = nWordId*nBaseTypeNum+3;
				break;
			case 't': nWordId = nWordId*nBaseTypeNum+3;
				break;
			default: nWordId = nWordId*nBaseTypeNum;
				nPass = 0;
		}

		if(nPass == 1)
		{
			pStatus->pMatElement[0] = 1;
			pStatus->pMatElement[nSpeciesNum] = 1;
			for(nj=1; nj<nSpeciesNum; nj++)
			{
				if(vAln[nj]->m_pSequence->m_pString[ni] != chBase)
				{
					pStatus->pMatElement[nj] = 0;
					if(pStatus->pMatElement[nj] < pMotifC->pMatElement[nj])
					{
						nPass = 0;
						break;
					}
				}
				else
				{
					pStatus->pMatElement[nj] = 1;
					pStatus->pMatElement[nSpeciesNum] += 1;
				}
			}

			if(nPass == 1)
			{
				if(pStatus->pMatElement[nSpeciesNum] < pMotifC->pMatElement[nSpeciesNum])
				{
					nPass = 0;
				}
			}
		}
		
		if(nPass == 0)
		{
			pCons->pMatElement[ni] = 1;
			nBadNum += 1;
		}
	}

	/* scan */
	for(; ni<nAlnLen; ni++)
	{
		/* add new */
		chBase = vAln[0]->m_pSequence->m_pString[ni];
		nPass = 1;
		switch(chBase)
		{
			case 'A': nWordId = nWordId*nBaseTypeNum;
				break;
			case 'a': nWordId = nWordId*nBaseTypeNum;
				break;
			case 'C': nWordId = nWordId*nBaseTypeNum+1;
				break;
			case 'c': nWordId = nWordId*nBaseTypeNum+1;
				break;
			case 'G': nWordId = nWordId*nBaseTypeNum+2;
				break;
			case 'g': nWordId = nWordId*nBaseTypeNum+2;
				break;
			case 'T': nWordId = nWordId*nBaseTypeNum+3;
				break;
			case 't': nWordId = nWordId*nBaseTypeNum+3;
				break;
			default: nWordId = nWordId*nBaseTypeNum;
				nPass = 0;
		}

		if(nPass == 1)
		{
			pStatus->pMatElement[0] = 1;
			pStatus->pMatElement[nSpeciesNum] = 1;
			for(nj=1; nj<nSpeciesNum; nj++)
			{
				if(vAln[nj]->m_pSequence->m_pString[ni] != chBase)
				{
					pStatus->pMatElement[nj] = 0;
					if(pStatus->pMatElement[nj] < pMotifC->pMatElement[nj])
					{
						nPass = 0;
						break;
					}
				}
				else
				{
					pStatus->pMatElement[nj] = 1;
					pStatus->pMatElement[nSpeciesNum] += 1;
				}
			}

			if(nPass == 1)
			{
				if(pStatus->pMatElement[nSpeciesNum] < pMotifC->pMatElement[nSpeciesNum])
				{
					nPass = 0;
				}
			}
		}
		
		if(nPass == 0)
		{
			pCons->pMatElement[ni] = 1;
			nBadNum += 1;
		}

		/* count */
		if(nBadNum == 0)
		{
			pKmerCount->pMatElement[nWordId] += 1.0;
		}

		/* remove old */
		chBase = vAln[0]->m_pSequence->m_pString[ni-nKmerLen+1];
		switch(chBase)
		{
			case 'A': nWordId = nWordId;
				break;
			case 'a': nWordId = nWordId;
				break;
			case 'C': nWordId = nWordId-nScale;
				break;
			case 'c': nWordId = nWordId-nScale;
				break;
			case 'G': nWordId = nWordId-2*nScale;
				break;
			case 'g': nWordId = nWordId-2*nScale;
				break;
			case 'T': nWordId = nWordId-3*nScale;
				break;
			case 't': nWordId = nWordId-3*nScale;
				break;
			default: nWordId = nWordId;
		}
		if(pCons->pMatElement[ni-nKmerLen+1] == 1)
		{
			nBadNum -= 1;
		}
	}
	
	/* clear memory */
	DestroyIntMatrix(pStatus);
	DestroyByteMatrix(pCons);
	for(ni=0; ni<nSpeciesNum; ni++)
	{
		SequenceDelete(vAln[ni]);
	}
	free(vAln);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  BLASTHITCREATE: create blast hit structure.                            */
/* ----------------------------------------------------------------------- */ 
struct BLASTHIT *BLASTHITCREATE()
{
	/* define */
	struct BLASTHIT *pBlastHit = NULL;

	/* create */
	pBlastHit = (struct BLASTHIT *)calloc(1, sizeof(struct BLASTHIT));
	pBlastHit->pQAln = NULL;
	pBlastHit->pHAln = NULL;
	pBlastHit->pNext = NULL;
	pBlastHit->nQStart = -1;
	pBlastHit->nQEnd = -1;
	pBlastHit->chQStrand = '?';
	pBlastHit->nHStart = -1;
	pBlastHit->nHEnd = -1;
	pBlastHit->chHStrand = '?';

	/* return */
	return pBlastHit;
}

/* ----------------------------------------------------------------------- */ 
/*  BLASTHITDESTROY: destroy an instance of blast hit structure.           */
/* ----------------------------------------------------------------------- */ 
void BLASTHITDESTROY(struct BLASTHIT *pBlastHit)
{
	if(pBlastHit != NULL)
	{
		DeleteString(pBlastHit->pQAln);
		DeleteString(pBlastHit->pHAln);
		pBlastHit->pNext = NULL;
		free(pBlastHit);
	}
}

/* ----------------------------------------------------------------------- */ 
/*  BLASTHITCLEARLIST: clear all elements in a linear list of blast hits.  */
/* ----------------------------------------------------------------------- */ 
void BLASTHITCLEARLIST(struct BLASTHIT **pHitList)
{
	/* define */
	struct BLASTHIT *pBlastHit;

	if(pHitList != NULL)
	{
		while(*pHitList != NULL)
		{
			pBlastHit = *pHitList;
			*pHitList = pBlastHit->pNext;
			pBlastHit->pNext = NULL;
			BLASTHITDESTROY(pBlastHit);
		}
	}
}

/* ----------------------------------------------------------------------- */ 
/*  BLASTHITLOADFROMBLS: load blast results from a normal blast output.    */
/* ----------------------------------------------------------------------- */ 
struct BLASTHIT *BLASTHITLOADFROMBLS(char strFileName[])
{
	/* define */
	struct BLASTHIT *pHitList = NULL;
	struct BLASTHIT *pNewHit = NULL;
	struct BLASTHIT *pPrev = NULL;
	FILE *fpIn;
	char strLine[MED_LINE_LENGTH];
	char strQuery[LINE_LENGTH];
	char strHit[LINE_LENGTH];
	char strTemp1[LINE_LENGTH];
	char strTemp2[LINE_LENGTH];
	char strTemp3[LINE_LENGTH];
	char strTemp4[LINE_LENGTH];
	char strTemp5[LINE_LENGTH];
	char strTemp6[LINE_LENGTH];
	char strTemp7[LINE_LENGTH];
	char strTemp8[LINE_LENGTH];
	char *chp1,*chp2;

	/* open file */
	fpIn = NULL;
	fpIn = fopen(strFileName, "r");
	if(fpIn == NULL)
	{
		printf("Error: BLASTHITLOADFROMBLS, cannot open blast result file!\n");
		exit(EXIT_FAILURE);
	}

	/* load blast result */
	while( fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL )
	{
		/* load query */
		if( strstr(strLine, "Query=") != NULL)
		{
			sscanf(strLine, "%s %s", strTemp1, strQuery);
		}

		/* load subject */
		else if( strLine[0] == '>' )
		{
			StrTrimRight(strLine);
			strcpy(strHit, strLine+1);
		}

		/* load result */
		else if( strstr(strLine, "Score =") != NULL )
		{
			/* add old hit to list */
			if(pNewHit != NULL)
			{
				if( (pNewHit->pQAln->m_nLength != pNewHit->nAlnLen) ||
					(pNewHit->pHAln->m_nLength != pNewHit->nAlnLen) )
				{
					printf("Error: BLASTHITLOADFROMBLS, blast result loaded incorrectly!\n");
					exit(EXIT_FAILURE);
				}

				if(pHitList == NULL)
				{
					pHitList = pNewHit;
					pPrev = pNewHit;
				}
				else
				{
					pPrev->pNext = pNewHit;
					pPrev = pNewHit;
				}
			}

			/* create a new hit */
			pNewHit = NULL;
			pNewHit = BLASTHITCREATE();
			if(pNewHit == NULL)
			{
				printf("Error: BLASTHITLOADFROMBLS, cannot create blast result!\n");
				exit(EXIT_FAILURE);
			}
			strcpy(pNewHit->strQuery, strQuery);
			strcpy(pNewHit->strHit, strHit);

			/* read score */
			sscanf(strLine, "%s %s %s %s %s %s %s %s", strTemp1, strTemp2, strTemp3, strTemp4,
				strTemp5, strTemp6, strTemp7, strTemp8);
			pNewHit->dScore = atof(strTemp3);
			pNewHit->dEValue = atof(strTemp8);
			
			chp1 = strstr(strTemp5, "(");
			chp1 += 1;
			chp2 = strstr(strTemp5, "),");
			*chp2 = '\0';
			pNewHit->dScore2 = atof(chp1);

			/* read "Identities =" */
			fgets(strLine, MED_LINE_LENGTH, fpIn);
			sscanf(strLine, "%s %s %s %s", strTemp1, strTemp2, strTemp3, strTemp4);
			chp1 = strstr(strTemp3, "/");
			chp2 = chp1+1;
			*chp1 = '\0'; 
			pNewHit->nAlnLen = atoi(chp2);
			pNewHit->dIdentities = atof(strTemp3);
			pNewHit->dIdentities /= (double)(pNewHit->nAlnLen);

			chp1 = strstr(strTemp4, "(");
			chp1 += 1;
			chp2 = strstr(strTemp4, "%)");
			*chp2 = '\0';
			if( fabs(pNewHit->dIdentities*100.0-atof(chp1)) > 2.0)
			{
				printf("Error: BLASTHITLOADFROMBLS, blast result loaded incorrectly!\n");
				exit(EXIT_FAILURE);
			}

			/* read "Strand =" */
			fgets(strLine, MED_LINE_LENGTH, fpIn);
			sscanf(strLine, "%s %s %s %s %s", strTemp1, strTemp2, strTemp3, strTemp4, strTemp5);
			if( strcmp(strTemp3, "Plus") == 0 )
			{
				pNewHit->chQStrand = '+';
			}
			else
			{
				pNewHit->chQStrand = '-';
			}

			if( strcmp(strTemp5, "Plus") == 0 )
			{
				pNewHit->chHStrand = '+';
			}
			else
			{
				pNewHit->chHStrand = '-';
			}
		}

		else if(strstr( strLine, "Query:") == strLine)
		{
			sscanf(strLine, "%s %s %s %s", strTemp1, strTemp2, strTemp3, strTemp4);
			if(pNewHit->chQStrand == '+')
			{
				if(pNewHit->nQStart < 0)
					pNewHit->nQStart = atoi(strTemp2);
				pNewHit->nQEnd = atoi(strTemp4);
			}
			else
			{
				if(pNewHit->nQEnd < 0)
					pNewHit->nQEnd = atoi(strTemp2);
				pNewHit->nQStart = atoi(strTemp4);
			}

			StringAddTail(&(pNewHit->pQAln), strTemp3);
		}

		else if(strstr( strLine, "Sbjct:") == strLine)
		{
			sscanf(strLine, "%s %s %s %s", strTemp1, strTemp2, strTemp3, strTemp4);
			if(pNewHit->chHStrand == '+')
			{
				if(pNewHit->nHStart < 0)
					pNewHit->nHStart = atoi(strTemp2);
				pNewHit->nHEnd = atoi(strTemp4);
			}
			else
			{
				if(pNewHit->nHEnd < 0)
					pNewHit->nHEnd = atoi(strTemp2);
				pNewHit->nHStart = atoi(strTemp4);
			}

			StringAddTail(&(pNewHit->pHAln), strTemp3);
		}

		else if(strstr( strLine, "Database:") != NULL)
		{
		}

		else
		{
		}
	}		

	/* add old hit to list */
	if(pNewHit != NULL)
	{
		if( (pNewHit->pQAln->m_nLength != pNewHit->nAlnLen) ||
			(pNewHit->pHAln->m_nLength != pNewHit->nAlnLen) )
		{
			printf("Error: BLASTHITLOADFROMBLS, blast result loaded incorrectly!\n");
			exit(EXIT_FAILURE);
		}

		if(pHitList == NULL)
		{
			pHitList = pNewHit;
			pPrev = pNewHit;
		}
		else
		{
			pPrev->pNext = pNewHit;
			pPrev = pNewHit;
		}
	}

	/* close file */
	fclose(fpIn);

	/* return */
	return pHitList;
}

/* ----------------------------------------------------------------------- */ 
/*  BLASTHIT_COLINEAR: check colinearity of two blast hits.                */
/*  return 1 if they are colinear, otherwise return 0.                     */
/* ----------------------------------------------------------------------- */ 
int BLASTHIT_COLINEAR(struct BLASTHIT *pHit1, struct BLASTHIT *pHit2)
{
	/* define */
	int nColinear = 1;
	int nQUp = 0;
	int nHUp = 0;
	int nDirec1,nDirec2;

	/* check colinearity */
	if( (pHit1 == NULL) || (pHit2 == NULL) )
	{
		return nColinear;
	}

	/* seq name */
	if( (strcmp(pHit1->strQuery, pHit2->strQuery) != 0) || 
		(strcmp(pHit1->strHit, pHit2->strHit) != 0) )
	{
		nColinear = 0;
		return nColinear;
	}

	/* seq direction */
	if(pHit1->chQStrand == pHit1->chHStrand)
	{
		nDirec1 = 0;
	}
	else
	{
		nDirec1 = 1;
	}

	if(pHit2->chQStrand == pHit2->chHStrand)
	{
		nDirec2 = 0;
	}
	else
	{
		nDirec2 = 1;
	}

	if(nDirec1 != nDirec2)
	{
		nColinear = 0;
		return nColinear;
	}

	/* overlap */
	if(pHit2->nQEnd < pHit1->nQStart)
	{
		nQUp = 1;
	}
	else if(pHit1->nQEnd < pHit2->nQStart)
	{
		nQUp = 0;
	}
	else
	{
		nColinear = 0;
		return nColinear;
	}

	if(pHit2->nHEnd < pHit1->nHStart)
	{
		nHUp = 1;
	}
	else if(pHit1->nHEnd < pHit2->nHStart)
	{
		nHUp = 0;
	}
	else
	{
		nColinear = 0;
		return nColinear;
	}

	/* cross */
	if(nDirec1 == 0)
	{
		if(nQUp != nHUp)
		{
			nColinear = 0;
			return nColinear;
		}
	}
	else
	{
		if(nQUp == nHUp)
		{
			nColinear = 0;
			return nColinear;
		}
	}

	/* return */
	return nColinear;
}

/* ----------------------------------------------------------------------- */ 
/*  BLASTHIT_INSERTBYQUERYLOC: insert a blast hit into hitlist. The list   */
/*  of hits should be sorted according to their location in query sequence */
/* ----------------------------------------------------------------------- */ 
int BLASTHIT_INSERTBYQUERYLOC(struct BLASTHIT *pBlastHit, struct BLASTHIT **pHitList)
{
	/* define */
	struct BLASTHIT *pHit,*pPrev;

	/* insert */
	if(pHitList == NULL)
	{
		printf("Error: BLASTHIT_INSERTBYQUERYLOC, null hitlist!\n");
		return PROC_FAILURE;
	}

	if(*pHitList == NULL)
	{
		*pHitList = pBlastHit;
	}
	else
	{
		pPrev = NULL;
		pHit = *pHitList;

		while(pHit != NULL)
		{
			if(BLASTHIT_RELATIVELOCBYQUERY(pHit, pBlastHit) > 0)
			{
				break;
			}

			pPrev = pHit;
			pHit = pHit->pNext;
		}

		/* insert */
		if(pPrev == NULL)
		{
			pBlastHit->pNext = *pHitList;
			*pHitList = pBlastHit;
		}
		else
		{
			pBlastHit->pNext = pPrev->pNext;
			pPrev->pNext = pBlastHit;
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  BLASTHIT_RELATIVELOCBYQUERY: compare relative position of two hits     */
/*  according to their positions in query sequences.                       */
/*  return -1 if pHit1 < pHit2;                                            */
/*          0 if pHit1 == pHit2;                                           */
/*          1 if pHit1 > pHit2                                             */
/* ----------------------------------------------------------------------- */ 
int BLASTHIT_RELATIVELOCBYQUERY(struct BLASTHIT *pHit1, struct BLASTHIT *pHit2)
{
	/* define */
	int nRelativePos = 0;
	int nSeqName;

	if( (pHit1 == NULL) || (pHit2 == NULL) )
	{
		printf("Error: BLASTHIT_RELATIVELOCBYQUERY, null hits!\n");
		exit(EXIT_FAILURE);
	}

	/* seq name */
	nSeqName = strcmp(pHit1->strQuery, pHit2->strQuery);
	if(nSeqName < 0)
	{
		nRelativePos = -1;
		return nRelativePos;
	}
	else if(nSeqName > 0)
	{
		nRelativePos = 1;
		return nRelativePos;
	}

	/* seq direction */
	if(pHit1->nQStart < pHit2->nQStart)
	{
		nRelativePos = -1;
	}
	else if(pHit1->nQStart > pHit2->nQStart)
	{
		nRelativePos = 1;
	}
	else
	{
		if(pHit1->nQEnd < pHit2->nQEnd)
		{
			nRelativePos = -1;
		}
		else if(pHit1->nQEnd > pHit2->nQEnd)
		{
			nRelativePos = 1;
		}
		else
		{
			nRelativePos = 0;
		}
	}

	/* return */
	return nRelativePos;
}

/* ----------------------------------------------------------------------- */ 
/*  BLASTHIT_TRIMOVERLAP: trim overlap between two blast results. The Gold */
/*  will be kept in its original form, while pBlastHit will be trimed.     */ 
/* ----------------------------------------------------------------------- */ 
int BLASTHIT_TRIMOVERLAP(struct BLASTHIT *pGoldHit, struct BLASTHIT **pBlastHit)
{
	/* define */
	struct BLASTHIT *pHitList = NULL;
	struct BLASTHIT *pPrev = NULL;
	struct BLASTHIT *pNewHit;
	int nQName,nHName;
	int nQOverlap,nHOverlap;
	int nPQOver,nPHOver;
	int nQStart,nQEnd,nHStart,nHEnd,nQP,nHP;
	int ni,nA1,nA2;
	char *vQAln,*vHAln;
	int nNewSegment = 0;
	int nNewAlnLen;


	/* init check */
	if(pGoldHit == NULL)
	{
		return PROC_SUCCESS;
	}

	if(pBlastHit == NULL)
	{
		return PROC_SUCCESS;
	}

	if(*pBlastHit == NULL)
	{
		return PROC_SUCCESS;
	}

	/* if both sequences are different, then return  */
	nQName = strcmp(pGoldHit->strQuery, (*pBlastHit)->strQuery);
	nHName = strcmp(pGoldHit->strHit, (*pBlastHit)->strHit);
	if( (nQName != 0) && (nHName !=0) )
	{
		return PROC_SUCCESS;
	}

	/* if not overlap, then return */
	if(nQName == 0)
	{
		if((*pBlastHit)->nQEnd < pGoldHit->nQStart)
		{
			nQOverlap = 0;
		}
		else if(pGoldHit->nQEnd < (*pBlastHit)->nQStart)
		{
			nQOverlap = 0;
		}
		else if( ((*pBlastHit)->nQStart>=pGoldHit->nQStart) && ((*pBlastHit)->nQEnd<=pGoldHit->nQEnd) )
		{
			nQOverlap = 2;
		}
		else
		{
			nQOverlap = 1;
		}
	}
	else
	{
		nQOverlap = 0;
	}

	if(nHName == 0)
	{
		if((*pBlastHit)->nHEnd < pGoldHit->nHStart)
		{
			nHOverlap = 0;
		}
		else if(pGoldHit->nHEnd < (*pBlastHit)->nHStart)
		{
			nHOverlap = 0;
		}
		else if( ((*pBlastHit)->nHStart>=pGoldHit->nHStart) && ((*pBlastHit)->nHEnd<=pGoldHit->nHEnd) )
		{
			nHOverlap = 2;
		}
		else
		{
			nHOverlap = 1;
		}
	}
	else
	{
		nHOverlap = 0;
	}

	/* if blasthit is included in gold hit, then trim all the hit */
	if( (nQOverlap == 2) || (nHOverlap == 2) )
	{
		BLASTHITDESTROY(*pBlastHit);
		*pBlastHit = NULL;
		return PROC_SUCCESS;
	}

	/* if the two hits are not overlapping */
	if( (nQOverlap == 0) && (nHOverlap == 0) )
	{
		return PROC_SUCCESS;
	}

	/* if overlap, then check base by base */
	if( ((*pBlastHit)->nAlnLen != (*pBlastHit)->pQAln->m_nLength) ||
		((*pBlastHit)->nAlnLen != (*pBlastHit)->pHAln->m_nLength) )
	{
		printf("Error: BLASTHIT_TRIMOVERLAP, alignment length not match!\n");
		exit(EXIT_FAILURE);
	}


	if((*pBlastHit)->chQStrand == '+')
	{
		nQP = (*pBlastHit)->nQStart-1;
	}
	else
	{
		nQP = (*pBlastHit)->nQEnd+1;
	}

	if((*pBlastHit)->chHStrand == '+')
	{
		nHP = (*pBlastHit)->nHStart-1;
	}
	else
	{
		nHP = (*pBlastHit)->nHEnd+1;
	}

	vQAln = (*pBlastHit)->pQAln->m_pString;
	vHAln = (*pBlastHit)->pHAln->m_pString;

	for(ni=0; ni<(*pBlastHit)->nAlnLen; ni++)
	{
		/* update current position */
		if(vQAln[ni] != '-')
		{
			if((*pBlastHit)->chQStrand == '+')
			{
				nQP += 1;
			}
			else
			{
				nQP -= 1;
			}
		}

		if(vHAln[ni] != '-')
		{
			if((*pBlastHit)->chHStrand == '+')
			{
				nHP += 1;
			}
			else
			{
				nHP -= 1;
			}
		}

		/* check if position is overlapping with gold hit */
		nPQOver = 0;
		nPHOver = 0;
		if(nQOverlap == 1)
		{
			if( (nQP>=pGoldHit->nQStart) && (nQP<=pGoldHit->nQEnd) && (vQAln[ni] != '-'))
			{
				nPQOver = 1;
			}
		}

		if(nHOverlap == 1)
		{
			if( (nHP>=pGoldHit->nHStart) && (nHP<=pGoldHit->nHEnd) && (vHAln[ni] != '-'))
			{
				nPHOver = 1;
			}
		}

		/* if point overlap */
		if((nPQOver == 1) || (nPHOver == 1))
		{
			/* create a new blast hit to record unconflicted alignment */
			if(nNewSegment == 1)
			{
				pNewHit = NULL;
				pNewHit = BLASTHITCREATE();
				strcpy(pNewHit->strQuery, (*pBlastHit)->strQuery);
				strcpy(pNewHit->strHit, (*pBlastHit)->strHit);
				pNewHit->chHStrand = (*pBlastHit)->chHStrand;
				pNewHit->chQStrand = (*pBlastHit)->chQStrand;
				pNewHit->dEValue = (*pBlastHit)->dEValue;
				pNewHit->dScore = (*pBlastHit)->dScore;
				pNewHit->dScore2 = (*pBlastHit)->dScore2;
				pNewHit->pNext = NULL;
				
				if((*pBlastHit)->chQStrand == '-')
				{
					pNewHit->nQStart = nQEnd;
					pNewHit->nQEnd = nQStart;
				}
				else
				{
					pNewHit->nQStart = nQStart;
					pNewHit->nQEnd = nQEnd;
				}

				if((*pBlastHit)->chHStrand == '-')
				{
					pNewHit->nHStart = nHEnd;
					pNewHit->nHEnd = nHStart;
				}
				else
				{
					pNewHit->nHStart = nHStart;
					pNewHit->nHEnd = nHEnd;	
				}

				nNewAlnLen = nA2-nA1+1;
				pNewHit->nAlnLen = nNewAlnLen;
				pNewHit->pQAln = CreateString(nNewAlnLen);
				pNewHit->pHAln = CreateString(nNewAlnLen);
				if( (pNewHit->pHAln == NULL) || (pNewHit->pQAln == NULL) )
				{
					printf("Error: BLASTHIT_TRIMOVERLAP, cannot create new blast hit!\n");
					exit(EXIT_FAILURE);
				}
				memcpy(pNewHit->pQAln->m_pString, vQAln+nA1, nNewAlnLen);
				memcpy(pNewHit->pHAln->m_pString, vHAln+nA1, nNewAlnLen);
				(pNewHit->pQAln->m_pString)[nNewAlnLen] = '\0';
				(pNewHit->pHAln->m_pString)[nNewAlnLen] = '\0';
				
				if(pNewHit->pQAln->m_pString[0] == '-')
				{
					if(pNewHit->chQStrand == '+')
					{
						pNewHit->nQStart += 1;
					}
					else
					{
						pNewHit->nQEnd -= 1;
					}
				}

				if(pNewHit->pHAln->m_pString[0] == '-')
				{
					if(pNewHit->chHStrand == '+')
					{
						pNewHit->nHStart += 1;
					}
					else
					{
						pNewHit->nHEnd -= 1;
					}
				}
	
				BLASTHITTRIMEND(pNewHit);
				BLASTHITUPDATEALNLENANDIDENTITY(pNewHit);
				
				if(pHitList == NULL)
				{
					pHitList = pNewHit;
					pPrev = pNewHit;
				}
				else
				{
					pPrev->pNext = pNewHit;
					pPrev = pNewHit;
				}
			}

			/* reinitialize new segment indicator */
			nNewSegment = 0;
		}
		else
		{
			if(nNewSegment == 0)
			{
				nQStart = nQP; 
				nHStart = nHP;
				nQEnd = nQP; 
				nHEnd = nHP;
				nNewSegment = 1;
				nA1 = ni;
				nA2 = ni;
			}
			else
			{
				nQEnd = nQP;
				nHEnd = nHP;
				nA2 = ni;
			}
		}
	}

	/* add last */
	/* create a new blast hit to record unconflicted alignment */
	if(nNewSegment == 1)
	{
		/* add the last one */
		pNewHit = NULL;
		pNewHit = BLASTHITCREATE();
		strcpy(pNewHit->strQuery, (*pBlastHit)->strQuery);
		strcpy(pNewHit->strHit, (*pBlastHit)->strHit);
		pNewHit->chHStrand = (*pBlastHit)->chHStrand;
		pNewHit->chQStrand = (*pBlastHit)->chQStrand;
		pNewHit->dEValue = (*pBlastHit)->dEValue;
		pNewHit->dScore = (*pBlastHit)->dScore;
		pNewHit->dScore2 = (*pBlastHit)->dScore2;
		pNewHit->pNext = NULL;
		
		if((*pBlastHit)->chQStrand == '-')
		{
			pNewHit->nQStart = nQEnd;
			pNewHit->nQEnd = nQStart;
		}
		else
		{
			pNewHit->nQStart = nQStart;
			pNewHit->nQEnd = nQEnd;
		}

		if((*pBlastHit)->chHStrand == '-')
		{
			pNewHit->nHStart = nHEnd;
			pNewHit->nHEnd = nHStart;
		}
		else
		{
			pNewHit->nHStart = nHStart;
			pNewHit->nHEnd = nHEnd;	
		}

		nNewAlnLen = nA2-nA1+1;
		pNewHit->nAlnLen = nNewAlnLen;
		pNewHit->pQAln = CreateString(nNewAlnLen);
		pNewHit->pHAln = CreateString(nNewAlnLen);
		if( (pNewHit->pHAln == NULL) || (pNewHit->pQAln == NULL) )
		{
			printf("Error: BLASTHIT_TRIMOVERLAP, cannot create new blast hit!\n");
			exit(EXIT_FAILURE);
		}
		memcpy(pNewHit->pQAln->m_pString, vQAln+nA1, nNewAlnLen);
		memcpy(pNewHit->pHAln->m_pString, vHAln+nA1, nNewAlnLen);
		(pNewHit->pQAln->m_pString)[nNewAlnLen] = '\0';
		(pNewHit->pHAln->m_pString)[nNewAlnLen] = '\0';
		
		if(pNewHit->pQAln->m_pString[0] == '-')
		{
			if(pNewHit->chQStrand == '+')
			{
				pNewHit->nQStart += 1;
			}
			else
			{
				pNewHit->nQEnd -= 1;
			}
		}

		if(pNewHit->pHAln->m_pString[0] == '-')
		{
			if(pNewHit->chHStrand == '+')
			{
				pNewHit->nHStart += 1;
			}
			else
			{
				pNewHit->nHEnd -= 1;
			}
		}

		BLASTHITTRIMEND(pNewHit);
		BLASTHITUPDATEALNLENANDIDENTITY(pNewHit);

		if(pHitList == NULL)
		{
			pHitList = pNewHit;
			pPrev = pNewHit;
		}
		else
		{
			pPrev->pNext = pNewHit;
			pPrev = pNewHit;
		}
	}

	nNewSegment = 0;

	/* update the blast hit list */
	BLASTHITDESTROY(*pBlastHit);
	*pBlastHit = pHitList;

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  BLASTHITUPDATEALNLENANDIDENTITY: update alignment length and identity. */
/*  of a blast hit.                                                        */
/* ----------------------------------------------------------------------- */ 
int BLASTHITUPDATEALNLENANDIDENTITY(struct BLASTHIT *pBlastHit)
{
	/* define */
	int nAlnLen;
	int nIdnNum;
	char *vQAln,*vHAln;
	int ni;

	/* check */
	if(pBlastHit == NULL)
	{
		return PROC_SUCCESS;
	}
	if(pBlastHit->pQAln == NULL)
	{
		if(pBlastHit->pHAln == NULL)
		{
			pBlastHit->nAlnLen = 0;
			pBlastHit->dIdentities = 0.0;
			return PROC_SUCCESS;
		}
		else
		{
			printf("Error: BLASTHITUPDATEALNLENANDIDENTITY, alignment lengths do not match!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* update */
	nAlnLen = strlen(pBlastHit->pQAln->m_pString);
	if( (nAlnLen != (int)strlen(pBlastHit->pHAln->m_pString)) ||
		(nAlnLen != pBlastHit->pQAln->m_nLength) || 
		(nAlnLen != pBlastHit->pHAln->m_nLength) )
	{
		printf("Error: BLASTHITUPDATEALNLENANDIDENTITY, alignment lengths do not match!\n");
		exit(EXIT_FAILURE);
	}

	pBlastHit->nAlnLen = nAlnLen;
	if(nAlnLen == 0)
	{
		pBlastHit->dIdentities = 0.0;
		return PROC_SUCCESS;
	}
	
	vQAln = pBlastHit->pQAln->m_pString;
	vHAln = pBlastHit->pHAln->m_pString;
	nIdnNum = 0;
	for(ni=0; ni<nAlnLen; ni++)
	{
		if( (vQAln[ni] == '-') || (vHAln[ni] == '-') )
		{
		}
		else if( (vQAln[ni] == vHAln[ni]) && (vQAln[ni] != 'n') && (vQAln[ni] != 'N'))
		{
			nIdnNum++;
		}
	}

	pBlastHit->dIdentities = (double)nIdnNum/(double)nAlnLen;

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  BLASTHITTRIMEND: trim non-sense aligment (i.e. gap to gap alignment)   */
/*  at both ends of a blast hit.                                           */
/* ----------------------------------------------------------------------- */ 
int BLASTHITTRIMEND(struct BLASTHIT *pBlastHit)
{
	/* define */
	int ni,nj,nA1,nA2;
	struct tagString *pNewQAln,*pNewHAln;
	char *vQAln,*vHAln;
	int nAlnLen;

	/* check */
	if(pBlastHit == NULL)
	{
		return PROC_SUCCESS;
	}

	if(pBlastHit->nAlnLen == 0)
	{
		return PROC_SUCCESS;
	}

	if( (pBlastHit->nAlnLen != pBlastHit->pQAln->m_nLength) ||
		(pBlastHit->nAlnLen != pBlastHit->pHAln->m_nLength) )
	{
		printf("Error: BLASTHITTRIMEND, blast hit alignment info incorrect!\n");
		exit(EXIT_FAILURE);
	}

	/* adjust both ends */
	vQAln = pBlastHit->pQAln->m_pString;
	vHAln = pBlastHit->pHAln->m_pString;
	nA1 = 0;
	nA2 = pBlastHit->nAlnLen-1;

	for(ni=0; ni<pBlastHit->nAlnLen; ni++)
	{
		nA1 = ni;
		if( (vQAln[ni] != '-') && (vHAln[ni] != '-') )
		{
			break;
		}
		else
		{
			if(vQAln[ni] != '-')
			{
				if(pBlastHit->chQStrand == '-')
				{
					pBlastHit->nQEnd -= 1;
				}
				else
				{
					pBlastHit->nQStart += 1;
				}
			}

			if(vHAln[ni] != '-')
			{
				if(pBlastHit->chHStrand == '-')
				{
					pBlastHit->nHEnd -= 1;
				}
				else
				{
					pBlastHit->nHStart += 1;
				}
			}
		}
	}

	for(nj=(pBlastHit->nAlnLen-1); nj>=0; nj--)
	{
		nA2 = nj;
		if( (vQAln[nj] != '-') && (vHAln[nj] != '-') )
		{
			break;
		}
		else
		{
			if(vQAln[nj] != '-')
			{
				if(pBlastHit->chQStrand == '-')
				{
					pBlastHit->nQStart += 1;
				}
				else
				{
					pBlastHit->nQEnd -= 1;
				}
			}

			if(vHAln[nj] != '-')
			{
				if(pBlastHit->chHStrand == '-')
				{
					pBlastHit->nHStart += 1;
				}
				else
				{
					pBlastHit->nHEnd -= 1;
				}
			}
		}
	}

	if(nj<ni)
	{
		pBlastHit->nAlnLen = 0;
		DeleteString(pBlastHit->pQAln);
		pBlastHit->pQAln = NULL;
		DeleteString(pBlastHit->pHAln);
		pBlastHit->pHAln = NULL;
		
		/* return */
		return PROC_SUCCESS;
	}


	nAlnLen = nA2-nA1+1;
	pNewQAln = NULL;
	pNewQAln = CreateString(nAlnLen);
	pNewHAln = NULL;
	pNewHAln = CreateString(nAlnLen);
	if( (pNewQAln == NULL) || (pNewHAln == NULL) )
	{
		printf("Error: BLASTHITTRIMEND, cannot create new blast hits, nAlnLen = %d (%d-%d)!\n", nAlnLen, nA1, nA2);
		exit(EXIT_FAILURE);
	}


	memcpy(pNewQAln->m_pString, vQAln+nA1, nAlnLen);
	pNewQAln->m_pString[nAlnLen] = '\0';
	memcpy(pNewHAln->m_pString, vHAln+nA1, nAlnLen);
	pNewHAln->m_pString[nAlnLen] = '\0';

	pBlastHit->nAlnLen = nAlnLen;
	DeleteString(pBlastHit->pQAln);
	pBlastHit->pQAln = pNewQAln;
	DeleteString(pBlastHit->pHAln);
	pBlastHit->pHAln = pNewHAln;

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  BLASTHIT_INSERTBYALNLENANDIDENTITY: insert a blast hit into hitlist.   */
/*  The list of hits should be sorted according to their length and        */
/*  identity.                                                              */
/* ----------------------------------------------------------------------- */ 
int BLASTHIT_INSERTBYALNLENANDIDENTITY(struct BLASTHIT *pBlastHit, struct BLASTHIT **pHitList)
{
	/* define */
	struct BLASTHIT *pHit,*pPrev;

	/* insert */
	if(pHitList == NULL)
	{
		printf("Error: BLASTHIT_INSERTBYALNLENANDIDENTITY, null hitlist!\n");
		return PROC_FAILURE;
	}

	if(*pHitList == NULL)
	{
		*pHitList = pBlastHit;
	}
	else
	{
		pPrev = NULL;
		pHit = *pHitList;

		while(pHit != NULL)
		{
			if(BLASTHIT_RELATIVELOCBYALNLENANDIDENTITY(pHit, pBlastHit) > 0)
			{
				break;
			}

			pPrev = pHit;
			pHit = pHit->pNext;
		}

		/* insert */
		if(pPrev == NULL)
		{
			pBlastHit->pNext = *pHitList;
			*pHitList = pBlastHit;
		}
		else
		{
			pBlastHit->pNext = pPrev->pNext;
			pPrev->pNext = pBlastHit;
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  BLASTHIT_RELATIVELOCBYALNLENANDIDENTITY: compare relative position of  */
/*  two hits according to their length and identity.                       */
/*  return -1 if similarity level of pHit1 > pHit2;                        */
/*          0 if pHit1 == pHit2;                                           */
/*          1 if pHit1 < pHit2                                             */
/* ----------------------------------------------------------------------- */ 
int BLASTHIT_RELATIVELOCBYALNLENANDIDENTITY(struct BLASTHIT *pHit1, struct BLASTHIT *pHit2)
{
	/* define */
	int nRelativePos = 0;
	double dI1,dI2;
	
	/* init check */
	if( (pHit1 == NULL) || (pHit2 == NULL) )
	{
		printf("Error: BLASTHIT_RELATIVELOCBYALNLENANDIDENTITY, null hits!\n");
		exit(EXIT_FAILURE);
	}

	/* compare */
	dI1 = (double)(pHit1->nAlnLen)*(pHit1->dIdentities);
	dI2 = (double)(pHit2->nAlnLen)*(pHit2->dIdentities);
	
	if(dI1 > dI2)
	{
		nRelativePos = -1;
	}
	else if(dI1 < dI2)
	{
		nRelativePos = 1;
	}
	else
	{
		if(pHit1->nAlnLen < pHit2->nAlnLen)
		{
			nRelativePos = -1;
		}
		else if(pHit1->nAlnLen > pHit2->nAlnLen)
		{
			nRelativePos = 1;
		}
		else
		{
			nRelativePos = 0;
		}
	}
	

	/* return */
	return nRelativePos;
}

/* ----------------------------------------------------------------------- */ 
/*  BLASTHIT_INSERTBYHITLOC: insert a blast hit into hitlist. The list     */
/*  of hits should be sorted according to their location in hit sequence   */
/* ----------------------------------------------------------------------- */ 
int BLASTHIT_INSERTBYHITLOC(struct BLASTHIT *pBlastHit, struct BLASTHIT **pHitList)
{
	/* define */
	struct BLASTHIT *pHit,*pPrev;

	/* insert */
	if(pHitList == NULL)
	{
		printf("Error: BLASTHIT_INSERTBYHITLOC, null hitlist!\n");
		return PROC_FAILURE;
	}

	if(*pHitList == NULL)
	{
		*pHitList = pBlastHit;
	}
	else
	{
		pPrev = NULL;
		pHit = *pHitList;

		while(pHit != NULL)
		{
			if(BLASTHIT_RELATIVELOCBYHIT(pHit, pBlastHit) > 0)
			{
				break;
			}

			pPrev = pHit;
			pHit = pHit->pNext;
		}

		/* insert */
		if(pPrev == NULL)
		{
			pBlastHit->pNext = *pHitList;
			*pHitList = pBlastHit;
		}
		else
		{
			pBlastHit->pNext = pPrev->pNext;
			pPrev->pNext = pBlastHit;
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  BLASTHIT_RELATIVELOCBYHIT: compare relative position of two hits       */
/*  according to their positions in hit sequences.                         */
/*  return -1 if pHit1 < pHit2;                                            */
/*          0 if pHit1 == pHit2;                                           */
/*          1 if pHit1 > pHit2                                             */
/* ----------------------------------------------------------------------- */ 
int BLASTHIT_RELATIVELOCBYHIT(struct BLASTHIT *pHit1, struct BLASTHIT *pHit2)
{
	/* define */
	int nRelativePos = 0;
	int nSeqName;

	if( (pHit1 == NULL) || (pHit2 == NULL) )
	{
		printf("Error: BLASTHIT_RELATIVELOCBYHIT, null hits!\n");
		exit(EXIT_FAILURE);
	}

	/* seq name */
	nSeqName = strcmp(pHit1->strHit, pHit2->strHit);
	if(nSeqName < 0)
	{
		nRelativePos = -1;
		return nRelativePos;
	}
	else if(nSeqName > 0)
	{
		nRelativePos = 1;
		return nRelativePos;
	}

	/* seq direction */
	if(pHit1->nHStart < pHit2->nHStart)
	{
		nRelativePos = -1;
	}
	else if(pHit1->nHStart > pHit2->nHStart)
	{
		nRelativePos = 1;
	}
	else
	{
		if(pHit1->nHEnd < pHit2->nHEnd)
		{
			nRelativePos = -1;
		}
		else if(pHit1->nHEnd > pHit2->nHEnd)
		{
			nRelativePos = 1;
		}
		else
		{
			nRelativePos = 0;
		}
	}

	/* return */
	return nRelativePos;
}

/* ----------------------------------------------------------------------- */ 
/*  SEQLINKMAPCREATE: create sequence link structure.                      */
/* ----------------------------------------------------------------------- */ 
struct SEQLINKMAP *SEQLINKMAPCREATE()
{
	/* define */
	struct SEQLINKMAP *pLinkMap = NULL;

	/* create */
	pLinkMap = (struct SEQLINKMAP *)calloc(1, sizeof(struct SEQLINKMAP));
	if(pLinkMap == NULL)
	{
		printf("Error: SEQLINKMAPCREATE, cannot create linkmap structure!\n");
		exit(EXIT_FAILURE);
	}
	pLinkMap->pNext = NULL;
	pLinkMap->nQStart = -1;
	pLinkMap->nQEnd = -1;
	pLinkMap->chQStrand = '?';
	pLinkMap->nHStart = -1;
	pLinkMap->nHEnd = -1;
	pLinkMap->chHStrand = '?';
	pLinkMap->nGStart = -1;
	pLinkMap->nGEnd = -1;
	pLinkMap->chGStrand = '?';

	/* return */
	return pLinkMap;
}

/* ----------------------------------------------------------------------- */ 
/*  SEQLINKMAPDESTROY: destroy sequence link structure.                    */
/* ----------------------------------------------------------------------- */ 
void SEQLINKMAPDESTROY(struct SEQLINKMAP *pLinkMap)
{
	if(pLinkMap != NULL)
	{
		pLinkMap->pNext = NULL;
		free(pLinkMap);
	}
}

/* ----------------------------------------------------------------------- */ 
/*  SEQLINKMAPCLEARLIST: clear all elements in a linear list of linkmaps.  */
/* ----------------------------------------------------------------------- */ 
void SEQLINKMAPCLEARLIST(struct SEQLINKMAP **pLinkMapList)
{
	/* define */
	struct SEQLINKMAP *pLinkMap;

	if(pLinkMapList != NULL)
	{
		while(*pLinkMapList != NULL)
		{
			pLinkMap = *pLinkMapList;
			*pLinkMapList = pLinkMap->pNext;
			pLinkMap->pNext = NULL;
			SEQLINKMAPDESTROY(pLinkMap);
		}
	}
}


/* ----------------------------------------------------------------------- */ 
/*  SEQLINKMAP_INSERTBYQUERYLOC: insert a seq link map into a list. The    */
/*  list should be sorted according to their location in query sequence    */
/* ----------------------------------------------------------------------- */ 
int SEQLINKMAP_INSERTBYQUERYLOC(struct SEQLINKMAP *pLinkMap, struct SEQLINKMAP **pLinkMapList)
{
	/* define */
	struct SEQLINKMAP *pLink,*pPrev;

	/* insert */
	if(pLinkMapList == NULL)
	{
		printf("Error: SEQLINKMAP_INSERTBYQUERYLOC, null hitlist!\n");
		return PROC_FAILURE;
	}

	if(*pLinkMapList == NULL)
	{
		*pLinkMapList = pLinkMap;
	}
	else
	{
		pPrev = NULL;
		pLink = *pLinkMapList;

		while(pLink != NULL)
		{
			if(SEQLINKMAP_RELATIVELOCBYQUERY(pLink, pLinkMap) > 0)
			{
				break;
			}

			pPrev = pLink;
			pLink = pLink->pNext;
		}

		/* insert */
		if(pPrev == NULL)
		{
			pLinkMap->pNext = *pLinkMapList;
			*pLinkMapList = pLinkMap;
		}
		else
		{
			pLinkMap->pNext = pPrev->pNext;
			pPrev->pNext = pLinkMap;
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  SEQLINKMAP_RELATIVELOCBYQUERY: compare relative position of two links  */
/*  according to their positions in query sequences.                       */
/*  return -1 if pLink1 < pLink2;                                          */
/*          0 if pLink1 == pLink2;                                         */
/*          1 if pLink1 > pLink2                                           */
/* ----------------------------------------------------------------------- */ 
int SEQLINKMAP_RELATIVELOCBYQUERY(struct SEQLINKMAP *pLink1, struct SEQLINKMAP *pLink2)
{
	/* define */
	int nRelativePos = 0;
	int nSeqName;

	if( (pLink1 == NULL) || (pLink2 == NULL) )
	{
		printf("Error: SEQLINKMAP_RELATIVELOCBYQUERY, null hits!\n");
		exit(EXIT_FAILURE);
	}

	/* seq name */
	nSeqName = strcmp(pLink1->strQuery, pLink2->strQuery);
	if(nSeqName < 0)
	{
		nRelativePos = -1;
		return nRelativePos;
	}
	else if(nSeqName > 0)
	{
		nRelativePos = 1;
		return nRelativePos;
	}

	/* seq direction */
	if(pLink1->nQStart < pLink2->nQStart)
	{
		nRelativePos = -1;
	}
	else if(pLink1->nQStart > pLink2->nQStart)
	{
		nRelativePos = 1;
	}
	else
	{
		if(pLink1->nQEnd < pLink2->nQEnd)
		{
			nRelativePos = -1;
		}
		else if(pLink1->nQEnd > pLink2->nQEnd)
		{
			nRelativePos = 1;
		}
		else
		{
			nRelativePos = 0;
		}
	}

	/* return */
	return nRelativePos;
}

/* ----------------------------------------------------------------------- */ 
/*  MOTIFSITELINKMAPCREATE: create motif site link structure.              */
/* ----------------------------------------------------------------------- */ 
struct MOTIFSITELINKMAP *MOTIFSITELINKMAPCREATE()
{
	/* define */
	struct MOTIFSITELINKMAP *pLinkMap = NULL;

	/* create */
	pLinkMap = (struct MOTIFSITELINKMAP *)calloc(1, sizeof(struct MOTIFSITELINKMAP));
	if(pLinkMap == NULL)
	{
		printf("Error: MOTIFSITELINKMAPCREATE, cannot create linkmap structure!\n");
		exit(EXIT_FAILURE);
	}
	pLinkMap->pNext = NULL;
	pLinkMap->nSeqId = -1;
	strcpy(pLinkMap->strSiteSeq, "");
	strcpy(pLinkMap->strMotifSymbol, "");
	strcpy(pLinkMap->strMotifName, "");
	strcpy(pLinkMap->strQuery, "");
	strcpy(pLinkMap->strHit, "");
	strcpy(pLinkMap->strGenome, "");
	pLinkMap->nQStart = -1;
	pLinkMap->nQEnd = -1;
	pLinkMap->chQStrand = '?';
	pLinkMap->nHStart = -1;
	pLinkMap->nHEnd = -1;
	pLinkMap->chHStrand = '?';
	pLinkMap->nGStart = -1;
	pLinkMap->nGEnd = -1;
	pLinkMap->chGStrand = '?';

	/* return */
	return pLinkMap;
}

/* ----------------------------------------------------------------------- */ 
/*  MOTIFSITELINKMAPDESTROY: destroy motif site link structure.            */
/* ----------------------------------------------------------------------- */ 
void MOTIFSITELINKMAPDESTROY(struct MOTIFSITELINKMAP *pLinkMap)
{
	if(pLinkMap != NULL)
	{
		pLinkMap->pNext = NULL;
		free(pLinkMap);
	}
}

/* ----------------------------------------------------------------------- */ 
/*  MOTIFSITELINKMAPCLEARLIST: clear all elements in a linear list of      */
/*  motif site linkmaps.                                                   */
/* ----------------------------------------------------------------------- */ 
void MOTIFSITELINKMAPCLEARLIST(struct MOTIFSITELINKMAP **pLinkMapList)
{
	/* define */
	struct MOTIFSITELINKMAP *pLinkMap;

	if(pLinkMapList != NULL)
	{
		while(*pLinkMapList != NULL)
		{
			pLinkMap = *pLinkMapList;
			*pLinkMapList = pLinkMap->pNext;
			pLinkMap->pNext = NULL;
			MOTIFSITELINKMAPDESTROY(pLinkMap);
		}
	}
}

/* ----------------------------------------------------------------------- */ 
/*  MOTIFSITELINKMAP_INSERTBYQUERYLOC: insert a motif site link map into a */
/*  list. The list should be sorted according to their location in query   */
/*  sequence.                                                              */
/* ----------------------------------------------------------------------- */ 
int MOTIFSITELINKMAP_INSERTBYQUERYLOC(struct MOTIFSITELINKMAP *pLinkMap, 
				struct MOTIFSITELINKMAP **pLinkMapList)
{
	/* define */
	struct MOTIFSITELINKMAP *pLink,*pPrev;

	/* insert */
	if(pLinkMapList == NULL)
	{
		printf("Error: MOTIFSITELINKMAP_INSERTBYQUERYLOC, null list!\n");
		return PROC_FAILURE;
	}

	if(*pLinkMapList == NULL)
	{
		*pLinkMapList = pLinkMap;
	}
	else
	{
		pPrev = NULL;
		pLink = *pLinkMapList;

		while(pLink != NULL)
		{
			if(MOTIFSITELINKMAP_RELATIVELOCBYQUERY(pLink, pLinkMap) > 0)
			{
				break;
			}

			pPrev = pLink;
			pLink = pLink->pNext;
		}

		/* insert */
		if(pPrev == NULL)
		{
			pLinkMap->pNext = *pLinkMapList;
			*pLinkMapList = pLinkMap;
		}
		else
		{
			pLinkMap->pNext = pPrev->pNext;
			pPrev->pNext = pLinkMap;
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  MOTIFSITELINKMAP_RELATIVELOCBYQUERY: compare relative position of two  */
/*  motif site links according to their positions in query sequences.      */
/*  return -1 if pLink1 < pLink2;                                          */
/*          0 if pLink1 == pLink2;                                         */
/*          1 if pLink1 > pLink2                                           */
/* ----------------------------------------------------------------------- */ 
int MOTIFSITELINKMAP_RELATIVELOCBYQUERY(struct MOTIFSITELINKMAP *pLink1, 
				struct MOTIFSITELINKMAP *pLink2)
{
	/* define */
	int nRelativePos = 0;
	int nSeqName;

	if( (pLink1 == NULL) || (pLink2 == NULL) )
	{
		printf("Error: MOTIFSITELINKMAP_RELATIVELOCBYQUERY, null hits!\n");
		exit(EXIT_FAILURE);
	}

	/* seq name */
	nSeqName = strcmp(pLink1->strQuery, pLink2->strQuery);
	if(nSeqName < 0)
	{
		nRelativePos = -1;
		return nRelativePos;
	}
	else if(nSeqName > 0)
	{
		nRelativePos = 1;
		return nRelativePos;
	}

	/* seq direction */
	if(pLink1->nQStart < pLink2->nQStart)
	{
		nRelativePos = -1;
	}
	else if(pLink1->nQStart > pLink2->nQStart)
	{
		nRelativePos = 1;
	}
	else
	{
		if(pLink1->nQEnd < pLink2->nQEnd)
		{
			nRelativePos = -1;
		}
		else if(pLink1->nQEnd > pLink2->nQEnd)
		{
			nRelativePos = 1;
		}
		else
		{
			nRelativePos = 0;
		}
	}

	/* return */
	return nRelativePos;
}


