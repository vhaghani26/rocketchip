/* ----------------------------------------------------------------------- */
/*  AffyLib.c : implementation of the affymetrix library                   */
/*  Author : Ji HongKai ; Time: 2004.08                                    */
/* ----------------------------------------------------------------------- */

#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
#include "limits.h"

#include "MathLib.h"
#include "MatrixLib.h"
#include "StringLib.h"
#include "SequenceLib.h"
#include "GenomeLib.h"
#include "AffyLib.h"

/* ----------------------------------------------------------------------- */ 
/*                            Global Variables                             */
/* ----------------------------------------------------------------------- */ 
int AFFY_IS_BIG_ENDIAN = 1;

/* ----------------------------------------------------------------------- */ 
/* Load BIG_ENDIAN and LITTLE_ENDIAN files.                                */
/* ----------------------------------------------------------------------- */ 
void reverse_buf(char *buf, int size)
{
	int i;
	char temp;

	for (i=0; i <size/2; i++)
	{
		temp = buf[i];
		buf[i] = buf[size-1-i];
		buf[size-1-i] = temp;
	}
}

size_t big_endian_fread(void *ptr, size_t size, size_t count, FILE *stream, 
						int little_endian_machine)
{
	size_t result;
	char *temp;
	int i;

	temp = (char *)calloc(count, size);
	result = fread(temp, size, count, stream);

	if(little_endian_machine == 1)
	{
		for(i=0; i<(int)count; i++)
		{
			reverse_buf(temp+(int)size*i, (int)size);
		}
	}

	memcpy(ptr, temp, size*count);
	free(temp);
	
	return result;
}

size_t big_endian_fwrite(const void *ptr, size_t size, size_t count, FILE *stream,
						 int little_endian_machine)
{
	size_t result;
	char *temp;
	int i;

	temp = (char *)calloc(count, size);
	memcpy(temp, ptr, size*count);

	if (little_endian_machine == 1)
	{
		for(i = 0; i<(int)count; i++)
		{
			reverse_buf(temp+(int)size*i, (int)size);
		}
	}
	result = fwrite(temp, size, count, stream);
	free(temp);

	return result;
}

size_t little_endian_fread(void *ptr, size_t size, size_t count, FILE *stream, 
						int little_endian_machine)
{
	size_t result;
	char *temp;
	int i;

	temp = (char *)calloc(count, size);
	result = fread(temp, size, count, stream);

	if(little_endian_machine == 0)
	{
		for(i=0; i<(int)count; i++)
		{
			reverse_buf(temp+(int)size*i, (int)size);
		}
	}

	memcpy(ptr, temp, size*count);
	free(temp);
	
	return result;
}

size_t little_endian_fwrite(const void *ptr, size_t size, size_t count, FILE *stream,
						 int little_endian_machine)
{
	size_t result;
	char *temp;
	int i;

	temp = (char *)calloc(count, size);
	memcpy(temp, ptr, size*count);

	if (little_endian_machine == 0)
	{
		for(i = 0; i<(int)count; i++)
		{
			reverse_buf(temp+(int)size*i, (int)size);
		}
	}
	result = fwrite(temp, size, count, stream);
	free(temp);

	return result;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_CELData_Create()                                                  */
/*  create a CELData object.                                               */
/* ----------------------------------------------------------------------- */ 
struct tagCELData *Affy_CELData_Create()
{
	/* define */
	struct tagCELData *pCELData = NULL;

	/* create */
	pCELData = (struct tagCELData *)calloc(1, sizeof(struct tagCELData));
	if(pCELData == NULL)
	{
		printf("Error: Affy_CELData_Create, cannot create CELData object!\n");
		return NULL;
	}

	/* return */
	return pCELData;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_CELData_Destroy()                                                 */
/*  delete a CELData object.                                               */
/* ----------------------------------------------------------------------- */ 
void Affy_CELData_Destroy(struct tagCELData **pCELData)
{
	int ni;

	if(pCELData != NULL)
	{
		if( (*pCELData) != NULL)
		{
			DestroyDoubleMatrix((*pCELData)->pIntensity);
			DestroyDoubleMatrix((*pCELData)->pSD);
			DestroyIntMatrix((*pCELData)->pPixelNum);

			DeleteString((*pCELData)->vHeader);
			DeleteString((*pCELData)->vAlgorithm);
			DeleteString((*pCELData)->vAlgorithmParam);

			DestroyIntMatrix((*pCELData)->pMaskedX);
			DestroyIntMatrix((*pCELData)->pMaskedY);

			DestroyIntMatrix((*pCELData)->pOutlierX);
			DestroyIntMatrix((*pCELData)->pOutlierY);

			DestroyIntMatrix((*pCELData)->pModifiedX);
			DestroyIntMatrix((*pCELData)->pModifiedY);
			DestroyDoubleMatrix((*pCELData)->pModifiedOrig);

			for(ni=0; ni<(*pCELData)->nSubGrids; ni++)
			{
				Affy_CELSubGrid_Destroy((*pCELData)->vSubGrids+ni);
			}
			free((*pCELData)->vSubGrids);

			free(*pCELData);
			*pCELData = NULL;
		}
	}
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_CELSubGrid_Create()                                               */
/*  create a CELSubGrid object.                                            */
/* ----------------------------------------------------------------------- */ 
struct tagCELSubGrid *Affy_CELSubGrid_Create()
{
	/* define */
	struct tagCELSubGrid *pCELSubGrid = NULL;

	/* create */
	pCELSubGrid = (struct tagCELSubGrid *)calloc(1, sizeof(struct tagCELSubGrid));
	if(pCELSubGrid == NULL)
	{
		printf("Error: Affy_CELSubGrid_Create, cannot create CELSubGrid object!\n");
		return NULL;
	}

	/* return */
	return pCELSubGrid;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_CELSubGrid_Destroy()                                              */
/*  delete a CELSubGrid object.                                            */
/* ----------------------------------------------------------------------- */ 
void Affy_CELSubGrid_Destroy(struct tagCELSubGrid **pCELSubGrid)
{
	if(pCELSubGrid != NULL)
	{
		if( (*pCELSubGrid) != NULL)
		{
			free(*pCELSubGrid);
			*pCELSubGrid = NULL;
		}
	}
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_CSVANNOT_To_GenomeAlign_200506()                                  */
/*  convert affy csv annotation file to genome mapping coordinates.        */
/* ----------------------------------------------------------------------- */ 
int Affy_CSVANNOT_To_GenomeAlign_200506(char strInFile[], char strOutFile[], char strSpecies[])
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;

	char strLine[LONG_LINE_LENGTH];
	char strTemp[MEDLONG_LINE_LENGTH];
	char strChrPos[MEDLONG_LINE_LENGTH];
	char strRefLine[MEDLONG_LINE_LENGTH];
	char *chi,*chs;
	char *chSep,*chp;
	char strProbe[LINE_LENGTH];
	char strGenSym[LINE_LENGTH];
	char strUniGene[LINE_LENGTH];
	char strLocusID[LINE_LENGTH];
	int nLocusID;

	struct tagAffyGenomeAlign *vAlign;
	struct tagAffyGenomeAlign *vTrans;
	struct tagAffyGenomeAlign *pNewAlign,*pPrevAlign;
	struct tagString *vRefId[1024];
	int nRefCount;
	int ni,nLen,nLineId;
	
	
	/* init */
	vAlign = NULL;
	vTrans = NULL;
	fpIn = NULL;
	fpOut = NULL;

	/* open file */
	fpIn = fopen(strInFile,"r");
	if(fpIn == NULL)
	{
		printf("Error: Affy_CSVANNOT_To_GenomeAlign_200506, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}
	fpOut = fopen(strOutFile,"w");
	if(fpOut == NULL)
	{
		fclose(fpIn);
		printf("Error: Affy_CSVANNOT_To_GenomeAlign_200506, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	/* read */
	nLineId = 0;
	fgets(strLine, LONG_LINE_LENGTH, fpIn);
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		nLineId++;

		/* parse */
		vAlign = NULL;
		vTrans = NULL;

		/* 1: probe */
		if(strchr(strLine, '\"') == strLine)
		{
			chs = strLine+1;
			chi = strstr(chs, "\",");
			*chi = '\0';
			chi++;
		}
		else
		{
			chs = strLine;
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		strcpy(strProbe, chs);
		chs = chi+1;
		/* printf("%s\n", strProbe); */
		

		/* 2: GeneChip Array */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 3: Species Scientific Name */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 4: Annotation Date */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 5: Sequence Type */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 6: Sequence Source */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 7: Transcript ID */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 8: Target Description */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;
		

		/* 9: Representative Public ID */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 10-a: Archival UniGene Cluster */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		
		/* 10-b: UniGene ID */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		strcpy(strUniGene, chs);
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 11: Genome Version */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 12: Alignments */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chs++;
			chi = strstr(chs, "\",");
			*chi = '\0';
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
			*chi = '\0';
		}		
		strcpy(strChrPos, chs);
		chs = chi+1;

		if(strcmp(strChrPos, "---") == 0)
		{	
		}
		else
		{
			chp = strChrPos;
			chSep = strstr(chp, "///");
			while(chSep != NULL)
			{
				*chSep = '\0';
				strcpy(strTemp, chp);

				pNewAlign = NULL;
				pNewAlign = AffyGenomeAlignCreate();
				AffyGenomeAlignInit(pNewAlign, strTemp, strSpecies);
				if(vAlign == NULL)
				{
					vAlign = pNewAlign;
					pPrevAlign = pNewAlign;
				}
				else
				{
					pPrevAlign->pNext = pNewAlign;
					pPrevAlign = pNewAlign;
				}

				chp = chSep+3;
				chSep = strstr(chp, "///");
			}


			strcpy(strTemp, chp);
			pNewAlign = NULL;
			pNewAlign = AffyGenomeAlignCreate();
			AffyGenomeAlignInit(pNewAlign, strTemp, strSpecies);
			if(vAlign == NULL)
			{
				vAlign = pNewAlign;
				pPrevAlign = pNewAlign;
			}
			else
			{
				pPrevAlign->pNext = pNewAlign;
				pPrevAlign = pNewAlign;
			}
		}
			

		/* 13: Overlapping Transcripts */
		/* StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chs++;
			chi = strstr(chs, "\",");
			*chi = '\0';
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
			*chi = '\0';
		}
		strcpy(strChrPos, chs);

		chs = chi+1;

		if(strcmp(strChrPos, "---") == 0)
		{	
		}
		else
		{
			chp = strChrPos;
			chSep = strstr(chp, "///");
			while(chSep != NULL)
			{
				*chSep = '\0';
				strcpy(strTemp, chp);

				pNewAlign = NULL;
				pNewAlign = AffyGenomeAlignCreate();
				AffyGenomeTransInit(pNewAlign, strTemp, strSpecies);
				if(vTrans == NULL)
				{
					vTrans = pNewAlign;
					pPrevAlign = pNewAlign;
				}
				else
				{
					pPrevAlign->pNext = pNewAlign;
					pPrevAlign = pNewAlign;
				}

				chp = chSep+3;
				chSep = strstr(chp, "///");
			}


			strcpy(strTemp, chp);
			pNewAlign = NULL;
			pNewAlign = AffyGenomeAlignCreate();
			AffyGenomeTransInit(pNewAlign, strTemp, strSpecies);
			if(vTrans == NULL)
			{
				vTrans = pNewAlign;
				pPrevAlign = pNewAlign;
			}
			else
			{
				pPrevAlign->pNext = pNewAlign;
				pPrevAlign = pNewAlign;
			}
		}
		*/

		/* 14: Gene Title */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 15: Gene Symbol */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chs++;
			chi = strstr(chs, "\",");
			*chi = '\0';
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
			*chi = '\0';
		}
		strcpy(strGenSym, chs);
		chs = chi+1;

		if(strcmp(strGenSym, "---") == 0)
			strcpy(strGenSym, "NA");

		/* 16: Chromosomal Location */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 17: UniGene ClusterType */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chs++;
			chi = strstr(chs, "\",");
			*chi = '\0';
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
			*chi = '\0';
		}
		/* strcpy(strUniGene, chs); */
		chs = chi+1;

		/* 18: Ensembl */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 19: LocusLink or Entrez Gene */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chs++;
			chi = strstr(chs, "\",");
			*chi = '\0';
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
			*chi = '\0';
		}
		strcpy(strLocusID, chs);
		chs = chi+1;

		if(strcmp(strLocusID, "---") == 0)
		{
			nLocusID = -1;
		}
		else if(strstr(strLocusID, "///") != NULL)
		{
			nLocusID = -1;
		}
		else
		{
			nLocusID = atoi(strLocusID);
		}

		/* 20: SwissProt */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 21: EC */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 22: OMIM */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 23: RefSeq Protein ID */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 24: RefSeq Transcript ID */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chs++;
			chi = strstr(chs, "\",");
			*chi = '\0';
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
			*chi = '\0';
		}
		strcpy(strRefLine, chs);
		chs = chi+1;

		nRefCount = 0;
		ni = 0;
		if(strcmp(strRefLine, "---") == 0)
		{	
		}
		else
		{
			chp = strRefLine;
			chSep = strstr(chp, "///");
			while(chSep != NULL)
			{
				*chSep = '\0';
				strcpy(strTemp, chp);
				StrTrimRight(strTemp);
				StrTrimLeft(strTemp);
				nLen = (int)strlen(strTemp)+1;
				vRefId[ni] = CreateString(nLen);
				strcpy(vRefId[ni]->m_pString, strTemp);
				ni++;

				chp = chSep+3;
				chSep = strstr(chp, "///");
			}

			strcpy(strTemp, chp);
			StrTrimRight(strTemp);
			StrTrimLeft(strTemp);
			nLen = (int)strlen(strTemp)+1;
			vRefId[ni] = CreateString(nLen);
			strcpy(vRefId[ni]->m_pString, strTemp);
			ni++;

			nRefCount = ni;
		}

		
		/* 25: FlyBase */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 26: AGI */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 27: WormBase */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 28: MGI Name */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 29: RGD Name */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 30: SGD accession number */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 31: Gene Ontology Biological Process */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 32: Gene Ontology Cellular Component */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 33: Gene Ontology Molecular Function */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 34: Pathway */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 35: Protein Families */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 36: Protein Domains */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 37: InterPro */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 38: Trans Membrane */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 39: QTL */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 40: Annotation Description */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;
		
		/* 41: Annotation Transcript Cluster */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 42: Transcript Assignments */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 43: Annotation Notes */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
		}
		else
		{
			chi = strchr(chs, ',');
		}
		if(chi != NULL)
		{
			printf("Error: Affy_CSVANNOT_To_Reduced_200408, csv annotation loading error!\n");
			exit(EXIT_FAILURE);
		}
		/* strcpy(strTemp, chs); */


		/* write */
		if(vAlign == NULL)
		{
			fprintf(fpOut, "%d\t%s\t-1\t-1\t-1\t?\t-1.0\t%d\t%s\n", 
						nLineId, strProbe, nLocusID, strGenSym);
		}
		else
		{
			pNewAlign = vAlign;
			while(pNewAlign != NULL)
			{
				fprintf(fpOut, "%d\t%s\t%d\t%d\t%d\t%c\t%f\t%d\t%s\n", 
						nLineId, strProbe, pNewAlign->nChr, pNewAlign->nStart, pNewAlign->nEnd-1,
						pNewAlign->chRc, pNewAlign->dIdentity, 
						nLocusID, strGenSym);
				pNewAlign = pNewAlign->pNext;
			}
		}

	
		/* destroy */
		while(vAlign != NULL)
		{
			pNewAlign = vAlign;
			vAlign = vAlign->pNext;
			AffyGenomeAlignDestroy(pNewAlign);
		}
		while(vTrans != NULL)
		{
			pNewAlign = vTrans;
			vTrans = vTrans->pNext;
			AffyGenomeAlignDestroy(pNewAlign);
		}
		for(ni=0; ni<nRefCount; ni++)
		{
			DeleteString(vRefId[ni]);
			vRefId[ni] = NULL;
		}
	}


	/* close file */
	fclose(fpIn);
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_CSVANNOT_To_Reduced_200408()                                      */
/*  convert affy csv annotation file to reduced annotation                 */
/* ----------------------------------------------------------------------- */ 
int Affy_CSVANNOT_To_Reduced_200408(char strInFile[], char strOutFile[], char strSpecies[])
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;

	char strLine[LONG_LINE_LENGTH];
	char strTemp[MEDLONG_LINE_LENGTH];
	char strChrPos[MEDLONG_LINE_LENGTH];
	char strRefLine[MEDLONG_LINE_LENGTH];
	char *chi,*chs;
	char *chSep,*chp;
	char strProbe[LINE_LENGTH];
	char strGenSym[MEDLONG_LINE_LENGTH];
	char strUniGene[LINE_LENGTH];
	char strLocusID[LINE_LENGTH];
	int nLocusID;

	struct tagAffyGenomeAlign *vAlign;
	struct tagAffyGenomeAlign *vTrans;
	struct tagAffyGenomeAlign *pNewAlign,*pPrevAlign,*pTarAlign;
	struct tagString *vRefId[1024];
	int nRefCount;
	int ni,nLen,nLineId;
	int nAlignWritten;
	
	/* init */
	vAlign = NULL;
	vTrans = NULL;
	fpIn = NULL;
	fpOut = NULL;

	/* open file */
	fpIn = fopen(strInFile,"rt");
	if(fpIn == NULL)
	{
		printf("Error: Affy_CSVANNOT_To_Reduced_200408, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}
	fpOut = fopen(strOutFile,"wt");
	if(fpOut == NULL)
	{
		fclose(fpIn);
		printf("Error: Affy_CSVANNOT_To_Reduced_200408, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	/* read */
	nLineId = 0;
	fgets(strLine, LONG_LINE_LENGTH, fpIn);
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		nLineId++;

		/* parse */
		vAlign = NULL;
		vTrans = NULL;

		/* 1: probe */
		if(strchr(strLine, '\"') == strLine)
		{
			chs = strLine+1;
			chi = strstr(chs, "\",");
			*chi = '\0';
			chi++;
		}
		else
		{
			chs = strLine;
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		strcpy(strProbe, chs);
		chs = chi+1;
		printf("%s\n", strProbe); 
		if(strcmp(strProbe, "AFFX-MURINE_B2_at") == 0)
		{
			chs=chs;
		}

		/* debug */
		/* if(strcmp(strProbe, "31807_at") == 0)
		{
			chs = chs;
		} */

		/* 2: GeneChip Array */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 3: Species Scientific Name */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 4: Annotation Date */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 5: Sequence Type */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 6: Sequence Source */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 7: Transcript ID */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 8: Target Description */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;
		

		/* 9: Representative Public ID */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 10-a: Archival UniGene Cluster */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		
		/* 10-b: UniGene ID */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		strcpy(strUniGene, chs);
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 11: Genome Version */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 12: Alignments */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chs++;
			chi = strstr(chs, "\",");
			*chi = '\0';
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
			*chi = '\0';
		}		
		strcpy(strChrPos, chs);
		chs = chi+1;

		if(strcmp(strChrPos, "---") == 0)
		{	
		}
		else
		{
			chp = strChrPos;
			chSep = strstr(chp, "///");
			while(chSep != NULL)
			{
				*chSep = '\0';
				strcpy(strTemp, chp);

				pNewAlign = NULL;
				pNewAlign = AffyGenomeAlignCreate();
				AffyGenomeAlignInit(pNewAlign, strTemp, strSpecies);
				if(vAlign == NULL)
				{
					vAlign = pNewAlign;
					pPrevAlign = pNewAlign;
				}
				else
				{
					pPrevAlign->pNext = pNewAlign;
					pPrevAlign = pNewAlign;
				}

				chp = chSep+3;
				chSep = strstr(chp, "///");
			}


			strcpy(strTemp, chp);
			pNewAlign = NULL;
			pNewAlign = AffyGenomeAlignCreate();
			AffyGenomeAlignInit(pNewAlign, strTemp, strSpecies);
			if(vAlign == NULL)
			{
				vAlign = pNewAlign;
				pPrevAlign = pNewAlign;
			}
			else
			{
				pPrevAlign->pNext = pNewAlign;
				pPrevAlign = pNewAlign;
			}
		}
			

		/* 13: Overlapping Transcripts */
		/* StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chs++;
			chi = strstr(chs, "\",");
			*chi = '\0';
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
			*chi = '\0';
		}
		strcpy(strChrPos, chs);

		chs = chi+1;

		if(strcmp(strChrPos, "---") == 0)
		{	
		}
		else
		{
			chp = strChrPos;
			chSep = strstr(chp, "///");
			while(chSep != NULL)
			{
				*chSep = '\0';
				strcpy(strTemp, chp);

				pNewAlign = NULL;
				pNewAlign = AffyGenomeAlignCreate();
				AffyGenomeTransInit(pNewAlign, strTemp, strSpecies);
				if(vTrans == NULL)
				{
					vTrans = pNewAlign;
					pPrevAlign = pNewAlign;
				}
				else
				{
					pPrevAlign->pNext = pNewAlign;
					pPrevAlign = pNewAlign;
				}

				chp = chSep+3;
				chSep = strstr(chp, "///");
			}


			strcpy(strTemp, chp);
			pNewAlign = NULL;
			pNewAlign = AffyGenomeAlignCreate();
			AffyGenomeTransInit(pNewAlign, strTemp, strSpecies);
			if(vTrans == NULL)
			{
				vTrans = pNewAlign;
				pPrevAlign = pNewAlign;
			}
			else
			{
				pPrevAlign->pNext = pNewAlign;
				pPrevAlign = pNewAlign;
			}
		}
		*/

		/* 14: Gene Title */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 15: Gene Symbol */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chs++;
			chi = strstr(chs, "\",");
			*chi = '\0';
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
			*chi = '\0';
		}
		strcpy(strGenSym, chs);
		chs = chi+1;

		if(strcmp(strGenSym, "---") == 0)
			strcpy(strGenSym, "NA");

		/* 16: Chromosomal Location */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 17: UniGene ClusterType */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chs++;
			chi = strstr(chs, "\",");
			*chi = '\0';
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
			*chi = '\0';
		}
		/* strcpy(strUniGene, chs); */
		chs = chi+1;

		/* 18: Ensembl */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 19: LocusLink or Entrez Gene */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chs++;
			chi = strstr(chs, "\",");
			*chi = '\0';
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
			*chi = '\0';
		}
		strcpy(strLocusID, chs);
		chs = chi+1;

		if(strcmp(strLocusID, "---") == 0)
		{
			nLocusID = -1;
		}
		else if(strstr(strLocusID, "///") != NULL)
		{
			nLocusID = -1;
		}
		else
		{
			nLocusID = atoi(strLocusID);
		}

		/* 20: SwissProt */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 21: EC */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 22: OMIM */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 23: RefSeq Protein ID */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 24: RefSeq Transcript ID */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chs++;
			chi = strstr(chs, "\",");
			*chi = '\0';
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
			*chi = '\0';
		}
		strcpy(strRefLine, chs);
		chs = chi+1;

		nRefCount = 0;
		ni = 0;
		if(strcmp(strRefLine, "---") == 0)
		{	
		}
		else
		{
			chp = strRefLine;
			chSep = strstr(chp, "///");
			while(chSep != NULL)
			{
				*chSep = '\0';
				strcpy(strTemp, chp);
				StrTrimRight(strTemp);
				StrTrimLeft(strTemp);
				nLen = (int)strlen(strTemp)+1;
				vRefId[ni] = CreateString(nLen);
				strcpy(vRefId[ni]->m_pString, strTemp);
				ni++;

				chp = chSep+3;
				chSep = strstr(chp, "///");
			}

			strcpy(strTemp, chp);
			StrTrimRight(strTemp);
			StrTrimLeft(strTemp);
			nLen = (int)strlen(strTemp)+1;
			vRefId[ni] = CreateString(nLen);
			strcpy(vRefId[ni]->m_pString, strTemp);
			ni++;

			nRefCount = ni;
		}

		
		/* 25: FlyBase */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 26: AGI */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 27: WormBase */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 28: MGI Name */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 29: RGD Name */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 30: SGD accession number */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 31: Gene Ontology Biological Process */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 32: Gene Ontology Cellular Component */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 33: Gene Ontology Molecular Function */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 34: Pathway */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 35: Protein Families */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 36: Protein Domains */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 37: InterPro */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 38: Trans Membrane */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 39: QTL */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 40: Annotation Description */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;
		
		/* 41: Annotation Transcript Cluster */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 42: Transcript Assignments */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
			chi++;
		}
		else
		{
			chi = strchr(chs, ',');
		}
		*chi = '\0';
		/* strcpy(strTemp, chs); */
		chs = chi+1;

		/* 43: Annotation Notes */
		StrTrimLeft(chs);
		if(strchr(chs, '\"') == chs)
		{
			chi = strstr((chs+1), "\",");
		}
		else
		{
			chi = strchr(chs, ',');
		}
		if(chi != NULL)
		{
			printf("Error: Affy_CSVANNOT_To_Reduced_200408, csv annotation loading error!\n");
			exit(EXIT_FAILURE);
		}
		/* strcpy(strTemp, chs); */


		/* write */
		if(nRefCount == 0)
		{
			fprintf(fpOut, "%s\t---\t%d\t%s\n", 
					strProbe, nLocusID, strGenSym);
		}
		else
		{
			for(ni=0; ni<nRefCount; ni++)
			{
				fprintf(fpOut, "%s\t%s\t%d\t%s\n", 
					strProbe, vRefId[ni]->m_pString, nLocusID, strGenSym);
			}
		}


		/* if(vAlign == NULL)
		{
			if(nRefCount == 0)
			{
				fprintf(fpOut, "%s\t-1\t-1\t-1\t?\t-1.0\tNA\t%d\t%s\n", 
					strProbe, nLocusID, strGenSym);
			}
			else
			{
				for(ni=0; ni<nRefCount; ni++)
				{
					fprintf(fpOut, "%s\t-1\t-1\t-1\t?\t-1.0\t%s\t%d\t%s\n", 
						strProbe, vRefId[ni]->m_pString, nLocusID, strGenSym);
				}
			}
		}
		else
		{
			pNewAlign = vAlign;
			while(pNewAlign != NULL)
			{
				if(nRefCount == 0)
				{
					fprintf(fpOut, "%s\t%d\t%d\t%d\t%c\t%f\tNA\t%d\t%s\n", 
						strProbe, pNewAlign->nChr, pNewAlign->nStart, pNewAlign->nEnd-1,
						pNewAlign->chRc, pNewAlign->dIdentity,
						nLocusID, strGenSym);
				}
				else
				{
					nAlignWritten = 0;
					for(ni=0; ni<nRefCount; ni++)
					{
						pTarAlign = NULL;
						pPrevAlign = vTrans;
						while(pPrevAlign != NULL)
						{
							if(strcmp(vRefId[ni]->m_pString, pPrevAlign->strRefId) == 0)
							{
								pTarAlign = pPrevAlign;
								break;
							}
							pPrevAlign = pPrevAlign->pNext;
						}
						
						if(pTarAlign != NULL)
						{
							if( (pTarAlign->nChr == pNewAlign->nChr) &&
								(pTarAlign->chRc == pNewAlign->chRc) )
							{
								if((pTarAlign->nStart >= pNewAlign->nStart) &&
									(pTarAlign->nStart <= pNewAlign->nEnd))
								{
									fprintf(fpOut, "%s\t%d\t%d\t%d\t%c\t%f\t%s\t%d\t%s\n", 
										strProbe, pNewAlign->nChr, pNewAlign->nStart, pNewAlign->nEnd-1,
										pNewAlign->chRc, pNewAlign->dIdentity, vRefId[ni]->m_pString,
										nLocusID, strGenSym);
									nAlignWritten = 1;
								}
								else if((pNewAlign->nStart >= pTarAlign->nStart) &&
									(pNewAlign->nStart <= pTarAlign->nEnd))
								{
									fprintf(fpOut, "%s\t%d\t%d\t%d\t%c\t%f\t%s\t%d\t%s\n", 
										strProbe, pNewAlign->nChr, pNewAlign->nStart, pNewAlign->nEnd-1,
										pNewAlign->chRc, pNewAlign->dIdentity, vRefId[ni]->m_pString,
										nLocusID, strGenSym);
									nAlignWritten = 1;
								}
							}
						}
						else
						{
							fprintf(fpOut, "%s\t%d\t%d\t%d\t%c\t%f\t%s\t%d\t%s\n", 
								strProbe, pNewAlign->nChr, pNewAlign->nStart, pNewAlign->nEnd-1,
								pNewAlign->chRc, pNewAlign->dIdentity, vRefId[ni]->m_pString,
								nLocusID, strGenSym);
							nAlignWritten = 1;
						}						
					}

					if(nAlignWritten == 0)
					{
						fprintf(fpOut, "%s\t%d\t%d\t%d\t%c\t%f\tNA\t%d\t%s\n", 
							strProbe, pNewAlign->nChr, pNewAlign->nStart, pNewAlign->nEnd-1,
							pNewAlign->chRc, pNewAlign->dIdentity,
							nLocusID, strGenSym);
					}
				}

				

				pNewAlign = pNewAlign->pNext;
			}
		} 
		*/
			
		/* destroy */
		while(vAlign != NULL)
		{
			pNewAlign = vAlign;
			vAlign = vAlign->pNext;
			AffyGenomeAlignDestroy(pNewAlign);
		}
		while(vTrans != NULL)
		{
			pNewAlign = vTrans;
			vTrans = vTrans->pNext;
			AffyGenomeAlignDestroy(pNewAlign);
		}
		for(ni=0; ni<nRefCount; ni++)
		{
			DeleteString(vRefId[ni]);
			vRefId[ni] = NULL;
		}
	}


	/* close file */
	fclose(fpIn);
	fclose(fpOut);

	/* return */
	return PROC_FAILURE;
}

/* ----------------------------------------------------------------------- */ 
/*  AffyGenomeAlignCreate()                                                */
/*  create genome alignment element                                        */
/* ----------------------------------------------------------------------- */ 
struct tagAffyGenomeAlign *AffyGenomeAlignCreate()
{
	/* define */
	struct tagAffyGenomeAlign *pGA;

	/* create */
	pGA = NULL;
	pGA = (struct tagAffyGenomeAlign *)calloc(1, sizeof(struct tagAffyGenomeAlign));
	if(pGA == NULL)
	{
		printf("Error: cannot create affygenomealign element!\n");
		return NULL;
	}

	/* return */
	return pGA;
}

/* ----------------------------------------------------------------------- */ 
/*  AffyGenomeAlignDestroy()                                               */
/*  create genome alignment element                                        */
/* ----------------------------------------------------------------------- */ 
int AffyGenomeAlignDestroy(struct tagAffyGenomeAlign *pGA)
{
	if(pGA == NULL)
		return PROC_SUCCESS;

	pGA->pNext = NULL;
	free(pGA);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  AffyGenomeAlignInit()                                                  */
/*  init alignment element                                                 */
/* ----------------------------------------------------------------------- */ 
int AffyGenomeAlignInit(struct tagAffyGenomeAlign *pGA, char strParam[], char strSpecies[])
{
	/* define */
	char *p1,*p2;
	char strLine[MED_LINE_LENGTH];
	char strTemp[LINE_LENGTH];

	/* init */
	if(pGA == NULL)
		return PROC_FAILURE;

	strcpy(strLine, strParam);

	p1 = strstr(strLine, "//");
	*p1 = '\0';
	p1 = p1+2;
	p2 = strstr(p1, "//");
	*p2 = '\0';
	p2 = p2+2;

	strcpy(strTemp, p1);
	StrTrimLeft(strTemp);
	StrTrimRight(strTemp);
	pGA->dIdentity = atof(strTemp);
	
	strcpy(strTemp, p2);
	StrTrimLeft(strTemp);
	StrTrimRight(strTemp);
	strcpy(pGA->strPhyMap, strTemp);

	StrTrimLeft(strLine);
	StrTrimRight(strLine);

	p1 = strstr(strLine, "(");
	*p1 = '\0';
	p1++;
	p2 = strstr(p1, ")");
	*p2 = '\0';
	strcpy(strTemp, p1);
	StrTrimLeft(strTemp);
	StrTrimRight(strTemp);
	if(strcmp(strTemp, "+") == 0)
	{
		pGA->chRc = '+';
	}
	else if(strcmp(strTemp, "-") == 0)
	{
		pGA->chRc = '-';
	}
	else
	{
		pGA->chRc = '?';
	}


	StrTrimRight(strLine);
	p1 = strstr(strLine, ":");
	*p1 = '\0';
	p1++;
	strcpy(strTemp, strLine);
	StrTrimRight(strTemp);
	strcpy(pGA->strChr, strTemp);
	pGA->nChr = Genome_ChromosomeName_To_Index(strTemp, strSpecies);
	p2 = strstr(p1, "-");
	*p2 = '\0';
	p2++;
	strcpy(strTemp, p1);
	StrTrimLeft(strTemp);
	StrTrimRight(strTemp);
	pGA->nStart = atoi(strTemp);
	strcpy(strTemp, p2);
	StrTrimLeft(strTemp);
	StrTrimRight(strTemp);
	pGA->nEnd = atoi(strTemp);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  AffyGenomeTransInit()                                                  */
/*  init alignment element                                                 */
/* ----------------------------------------------------------------------- */ 
int AffyGenomeTransInit(struct tagAffyGenomeAlign *pGA, char strParam[], char strSpecies[])
{
	/* define */
	char *p1,*p2;
	char strLine[MED_LINE_LENGTH];
	char strTemp[LINE_LENGTH];

	/* init */
	if(pGA == NULL)
		return PROC_FAILURE;

	strcpy(strLine, strParam);

	p1 = strstr(strLine, "//");
	*p1 = '\0';
	p1 = p1+2;
	strcpy(strTemp, strLine);
	StrTrimLeft(strTemp);
	StrTrimRight(strTemp);
	strcpy(pGA->strRefId, strTemp);


	p2 = strstr(p1, "//");
	*p2 = '\0';
	p2 = p2+2;
	strcpy(strLine, p2);

	StrTrimLeft(strLine);
	StrTrimRight(strLine);

	p1 = strstr(strLine, "(");
	*p1 = '\0';
	p1++;
	p2 = strstr(p1, ")");
	*p2 = '\0';
	strcpy(strTemp, p1);
	StrTrimLeft(strTemp);
	StrTrimRight(strTemp);
	if(strcmp(strTemp, "+") == 0)
	{
		pGA->chRc = '+';
	}
	else if(strcmp(strTemp, "-") == 0)
	{
		pGA->chRc = '-';
	}
	else
	{
		pGA->chRc = '?';
	}


	StrTrimRight(strLine);
	p1 = strstr(strLine, ":");
	*p1 = '\0';
	p1++;
	strcpy(strTemp, strLine);
	StrTrimRight(strTemp);
	strcpy(pGA->strChr, strTemp);
	pGA->nChr = Genome_ChromosomeName_To_Index(strTemp, strSpecies);
	p2 = strstr(p1, "-");
	*p2 = '\0';
	p2++;
	strcpy(strTemp, p1);
	StrTrimLeft(strTemp);
	StrTrimRight(strTemp);
	pGA->nStart = atoi(strTemp);
	strcpy(strTemp, p2);
	StrTrimLeft(strTemp);
	StrTrimRight(strTemp);
	pGA->nEnd = atoi(strTemp);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_PickRandomControlFromCSV()                                        */
/*  pick random control probesets from affymetrix's CSV annotation file.   */
/* ----------------------------------------------------------------------- */ 
int Affy_PickRandomControlFromCSV(char strInFile[], char strOutFile[], int nControlNum)
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;

	char strLine[LONG_LINE_LENGTH];
	int nProbeNum;
	int ni,nj;
	struct tagString **vProbeName;
	char strProbe[LINE_LENGTH];
	char *chi;

	struct DOUBLEMATRIX *pRandMat;
	struct DOUBLEMATRIX *pSortMat;
	struct LONGMATRIX *pSortId;
	
	/* init */
	fpIn = NULL;
	fpOut = NULL;
	vProbeName = NULL;
	nProbeNum = 0;

	/* open file */
	fpIn = fopen(strInFile,"rt");
	if(fpIn == NULL)
	{
		printf("Error: Affy_PickRandomControlFromCSV, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}
	/* get probenum */
	fgets(strLine, LONG_LINE_LENGTH, fpIn);
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		/* parse */
		nProbeNum++;
	}

	fclose(fpIn);

	/* init probe name vector */
	if(nProbeNum <= 0)
		return PROC_SUCCESS;

	if(nControlNum > nProbeNum)
	{
		printf("Error: Affy_PickRandomControlFromCSV, control probeset number should be <= total probeset number!\n");
		exit(EXIT_FAILURE);
	}

	vProbeName = (struct tagString **)calloc(nProbeNum, sizeof(struct tagString *));
	if(vProbeName == NULL)
	{
		printf("Error: Affy_PickRandomControlFromCSV, cannot allocate memory for loading probeset id's!\n");
		exit(EXIT_FAILURE);
	}

	/* open file again and load probeset id */
	fpIn = NULL;
	fpIn = fopen(strInFile,"rt");
	if(fpIn == NULL)
	{
		printf("Error: Affy_PickRandomControlFromCSV, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}

	/* read probeset id */
	ni = 0;
	fgets(strLine, LONG_LINE_LENGTH, fpIn);
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		/* parse */
		chi = strchr(strLine, ',');
		*chi = '\0';
		strcpy(strProbe, strLine);
		StrTrimRight(strProbe);

		StringAddTail((vProbeName+ni), strProbe);

		/* get next index */
		ni++;
	}

	/* close */
	fclose(fpIn);
	
	if(ni != nProbeNum)
	{
		printf("Error: Affy_PickRandomControlFromCSV, probeset number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* select probesets randomly */
	pRandMat = NULL;
	pSortMat = NULL;
	pSortId = NULL;
	pRandMat = DMRANDU(1, nProbeNum);
	DMSORTMERGEA_0(pRandMat, &pSortMat, &pSortId);

	/* write */
	fpOut = fopen(strOutFile,"wt");
	if(fpOut == NULL)
	{
		printf("Error: Affy_PickRandomControlFromCSV, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "probeset_id\n");
	for(ni=0; ni<nControlNum; ni++)
	{
		nj = (int)(LMGETAT(pSortId, 0, ni));
		fprintf(fpOut, "%s\n", vProbeName[nj]->m_pString);
	}

	fclose(fpOut);


	/* release probeset ids */
	DestroyDoubleMatrix(pRandMat);
	DestroyDoubleMatrix(pSortMat);
	DestroyLongMatrix(pSortId);

	for(ni=0; ni<nProbeNum; ni++)
	{
		DeleteString(vProbeName[ni]);
		vProbeName[ni] = NULL;
	}
	free(vProbeName);

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  Affy_LinkProbesetToAnnotation()                                        */
/*  Link probeset with annotation.                                         */
/* ----------------------------------------------------------------------- */ 
int Affy_LinkProbesetToAnnotation(char strProbesetFile[], char strAnnotationFile[], 
								  char strOutputFile[], int nClassId)
{
	/* define */
	FILE *fpIn = NULL;
	FILE *fpMap = NULL;
	FILE *fpOut = NULL;

	char strLine[MED_LINE_LENGTH];
	char strLine2[MED_LINE_LENGTH];
	char strProbe[LINE_LENGTH];
	char strProbe2[LINE_LENGTH];
	int nLen;
	char *chStart;

	int nLineId, nRank;
	double dScore, dFDR;

	/* open input list */
	fpIn = fopen(strProbesetFile, "rt");
	if(fpIn == NULL)
	{
		printf("Error: Affy_LinkProbesetToAnnotation, cannot open probeset file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = fopen(strOutputFile, "wt");
	if(fpOut == NULL)
	{
		fclose(fpIn);
		printf("Error: Affy_LinkProbesetToAnnotation, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	/* load probeset one by one */
	fgets(strLine, MED_LINE_LENGTH, fpIn);
	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%s %d %d %lf %lf ", strProbe, &nLineId, &nRank, &dScore, &dFDR);
		
		fpMap = NULL;
		fpMap = fopen(strAnnotationFile, "rt");
		if(fpMap == NULL)
		{
			printf("Warning: Affy_LinkProbesetToAnnotation, cannot open annotation file!\n");
			continue;
		}

		while(fgets(strLine2, MED_LINE_LENGTH, fpMap) != NULL)
		{
			StrTrimLeft(strLine2);
			StrTrimRight(strLine2);
			if(strLine2[0] == '\0')
				continue;

			sscanf(strLine2, "%s ", strProbe2);
			if(strcmp(strProbe, strProbe2) == 0)
			{
				nLen = (int)strlen(strProbe2);
				chStart = strLine2+nLen;
				StrTrimLeft(chStart);

				fprintf(fpOut, "%s\t%d\t%9.7e\t%9.7e\t%s\n", 
					strProbe, nClassId, dScore, dFDR, chStart);
			}
		}
		fclose(fpMap);
	}

	/* close output file */
	fclose(fpIn);
	fclose(fpOut);

	/* exit */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_CreateProbeRefGeneMap()                                           */
/*  map probeset to refseq.                                                */
/* ----------------------------------------------------------------------- */ 
int Affy_CreateProbeRefGeneMap(char strProbeLists[], char strRefGenePath[], 
						char strMapFile[], char strOutFile[], char strSpecies[])
{
	/* define */
	int nRefGeneNum;
	struct tagRefGene **vRefGene;
	struct BYTEMATRIX *pActive;
	
	int nProbeRefPairNum;
	struct tagString **vProbeRefPair;
	struct INTMATRIX *pProbeRefMap;
	struct INTMATRIX *pClassLabel;
	
	int nUniqRefNum;
	struct INTMATRIX *pRefGeneUniqMap;

	int nLastChrom;
	int nChromNum;
	int vChromDelimits[MAX_CHROMOSOME_NUM][2];
	int nD1;

	FILE *fpList;
	FILE *fpIn;
	FILE *fpOut;

	char strFileName[MED_LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	int ni,nj,nk;

	char strProbe[LINE_LENGTH];
	char strRefId[LINE_LENGTH];
	int nChr;
	int nStart,nEnd;
	char chStrand;
	int nP1,nP2;
	int nClassId;
	double dIdentity;

	/* ------------------------ */
	/* load refgene annotations */
	/* ------------------------ */

	/* get refgene number */
	fpIn = NULL;
	fpIn = fopen(strRefGenePath, "rt");
	if(fpIn == NULL)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, cannot open the refgene file!\n");
		exit(EXIT_FAILURE);
	}

	nRefGeneNum = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		nRefGeneNum++;
	}

	fclose(fpIn);

	if(nRefGeneNum <= 0)
	{
		printf("Warning: no refgenes available!\n");
		return PROC_SUCCESS;
	}

	/* init memory for storing refgene */
	vRefGene = NULL;
	vRefGene = (struct tagRefGene **)calloc(nRefGeneNum, sizeof(struct tagRefGene *));
	if(vRefGene == NULL)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, cannot create the memory for refgenes!\n");
		exit(EXIT_FAILURE);
	}

	/* load all refgene annotations */
	fpIn = NULL;
	fpIn = fopen(strRefGenePath, "rt");
	if(fpIn == NULL)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, cannot open the refgene file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	nChromNum = 0;
	nLastChrom = -MAX_CHROMOSOME_NUM;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		vRefGene[ni] = RefGeneCreate();
		RefGeneInit_FromGenomeLabFormat(vRefGene[ni], strLine, strSpecies);

		if(vRefGene[ni]->nChrom != nLastChrom)
		{
			if(nLastChrom != -MAX_CHROMOSOME_NUM)
			{
				vChromDelimits[nChromNum][0] = nD1;
				vChromDelimits[nChromNum][1] = ni-1;
				nChromNum++;
			}

			nLastChrom = vRefGene[ni]->nChrom;
			nD1 = ni;
		}

		ni++;
	}

	if(nLastChrom != -MAX_CHROMOSOME_NUM)
	{
		vChromDelimits[nChromNum][0] = nD1;
		vChromDelimits[nChromNum][1] = ni-1;
		nChromNum++;
	}

	fclose(fpIn);
	if(ni != nRefGeneNum)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, refgene number not match!\n");
		exit(EXIT_FAILURE);
	}

	pActive = NULL;
	pActive = CreateByteMatrix(1, nRefGeneNum);
	if(pActive == NULL)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, cannot create lable matrix for refgene!\n");
		exit(EXIT_FAILURE);
	}

	/* ---------------------------------------- */
	/* get active refgene and create intial map */
	/* ---------------------------------------- */
	
	/* get probe-refgene pair number */
	nProbeRefPairNum = 0;
	fpList = NULL;
	fpList = fopen(strProbeLists, "rt");
	if(fpList == NULL)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, cannot open probelists!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpList) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%d %s", &nClassId, strFileName);
		StrTrimLeft(strFileName);
		
		fpIn = NULL;
		fpIn = fopen(strFileName, "rt");
		if(fpIn == NULL)
		{
			printf("Error: Affy_CreateProbeRefGeneMap, cannot open probelists!\n");
			exit(EXIT_FAILURE);
		}

		while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;

			nProbeRefPairNum++;
		}

		fclose(fpIn);
	}

	fclose(fpList);

	/* init storage space */
	vProbeRefPair = NULL;
	vProbeRefPair = (struct tagString **)calloc(nProbeRefPairNum, sizeof(struct tagString *));
	if(vProbeRefPair == NULL)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, cannot create enough storage space for probe-refgene pair!\n");
		exit(EXIT_FAILURE);
	}

	pClassLabel = NULL;
	pClassLabel = CreateIntMatrix(1, nProbeRefPairNum);
	if(pClassLabel == NULL)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, cannot create enough storage space for class label!\n");
		exit(EXIT_FAILURE);
	}

	pProbeRefMap = NULL;
	pProbeRefMap = CreateIntMatrix(1, nProbeRefPairNum);
	if(pProbeRefMap == NULL)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, cannot create enough storage space for probe-refgene pair!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nProbeRefPairNum; ni++)
	{
		pProbeRefMap->pMatElement[ni] = -1;
	}


	/* create initial map */
	ni = 0;
	fpList = NULL;
	fpList = fopen(strProbeLists, "rt");
	if(fpList == NULL)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, cannot open probelists!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpList) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%d %s", &nClassId, strFileName);
		StrTrimLeft(strFileName);
		
		fpIn = NULL;
		fpIn = fopen(strFileName, "rt");
		if(fpIn == NULL)
		{
			printf("Error: Affy_CreateProbeRefGeneMap, cannot open probelists!\n");
			exit(EXIT_FAILURE);
		}

		while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;

			sscanf(strLine, "%s %d %d %d %c %lf %s ", strProbe, &nChr, &nStart, 
				&nEnd, &chStrand, &dIdentity, strRefId);

			pClassLabel->pMatElement[ni] = nClassId;
			StringAddTail((vProbeRefPair+ni) , strProbe);

			if( (nChr <= 0) || (nChr > nChromNum) || 
				(dIdentity < 80.0) || (strcmp(strRefId, "NA") == 0))
			{
				ni++;
				continue;
			}

			nChr = nChr-1;
			nP1 = vChromDelimits[nChr][0];
			nP2 = vChromDelimits[nChr][1];

			for(nk=nP1; nk<=nP2; nk++)
			{
				if(strcmp(strRefId, vRefGene[nk]->strName) == 0)
				{
					if( (chStrand == vRefGene[nk]->chStrand) && 
						((nEnd-vRefGene[nk]->nTxStart) > -MAX_MATCH_GAP) &&
						((vRefGene[nk]->nTxEnd-nStart) > -MAX_MATCH_GAP) )
					{
						pActive->pMatElement[nk] = 1;
						pProbeRefMap->pMatElement[ni] = nk;
						break;
					}
				}
			}

			ni++;
		}

		fclose(fpIn);
	}

	fclose(fpList);

	if(ni != nProbeRefPairNum)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, probe-refgene pair number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* ---------------------------------------- */
	/* get unique refgene list and save map     */
	/* ---------------------------------------- */
	/* get the number of unique refgene */
	nUniqRefNum = 0;
	for(ni=0; ni<nRefGeneNum; ni++)
	{
		if(pActive->pMatElement[ni] == 1)
			nUniqRefNum++;
	}

	if(nUniqRefNum <= 0)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, unique refgene number = 0!\n");
		exit(EXIT_FAILURE);
	}

	/* get uniq gene map */
	pRefGeneUniqMap = NULL;
	pRefGeneUniqMap = CreateIntMatrix(1, nRefGeneNum);
	if(pRefGeneUniqMap == NULL)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, cannot create uniq refgene map!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutFile, "wt");
	if(fpOut == NULL)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	nj = 0;
	for(ni=0; ni<nRefGeneNum; ni++)
	{
		if(pActive->pMatElement[ni] == 1)
		{
			pRefGeneUniqMap->pMatElement[ni] = nj;
			RefGeneWrite(vRefGene[ni], fpOut);
			nj++;
		}
		else
		{
			pRefGeneUniqMap->pMatElement[ni] = -1;
		}
	}

	if(nj != nUniqRefNum)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, unique refgene number not match!\n");
		exit(EXIT_FAILURE);
	}

	fclose(fpOut);

	/* output map */
	fpOut = NULL;
	fpOut = fopen(strMapFile, "wt");
	if(fpOut == NULL)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nProbeRefPairNum; ni++)
	{
		nj = pProbeRefMap->pMatElement[ni];
		if(nj >= 0)
		{
			nk = pRefGeneUniqMap->pMatElement[nj];
			if(nk < 0)
			{
				printf("Error: Affy_CreateProbeRefGeneMap, logic error in mapping!\n");
				exit(EXIT_FAILURE);
			}
			fprintf(fpOut, "%s\t%d\t%d\n", vProbeRefPair[ni]->m_pString,
				pClassLabel->pMatElement[ni], nk);
		}
	}

	fclose(fpOut);
		
	/* ---------------------------------------- */
	/* release memory                           */
	/* ---------------------------------------- */
	for(ni=0; ni<nRefGeneNum; ni++)
	{
		RefGeneDestroy(vRefGene[ni]);
		vRefGene[ni] = NULL;
	}
	free(vRefGene);
	
	DestroyByteMatrix(pActive);

	for(ni=0; ni<nProbeRefPairNum; ni++)
	{
		DeleteString(vProbeRefPair[ni]);
		vProbeRefPair[ni] = NULL;
	}
	free(vProbeRefPair);

	DestroyIntMatrix(pClassLabel);
	DestroyIntMatrix(pProbeRefMap);
	DestroyIntMatrix(pRefGeneUniqMap);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_CreateProbeRefGeneMap_BothDirection()                             */
/*  map probeset to refseq.                                                */
/* ----------------------------------------------------------------------- */ 
int Affy_CreateProbeRefGeneMap_BothDirection(char strProbeLists[], char strRefGenePath[], 
						char strOutFile[], char strSpecies[])
{
	/* define */
	int nRefGeneNum;
	int nClassNum;
	struct tagRefGene **vRefGene;
	
	int nProbeRefPairNum;
	struct tagString **vProbeRefPair;
	struct DOUBLEMATRIX *pProbeVec;
	struct DOUBLEMATRIX *pRefGeneVec;
	struct INTMATRIX *pClassVec;
	struct DOUBLEMATRIX *pScoreVec;
	struct DOUBLEMATRIX *pIdentityVec;
	struct BYTEMATRIX *pFilterVec;
	struct BYTEMATRIX *pFinalVec;

	struct INTMATRIX *pClassIndicator;

	struct DOUBLEMATRIX *pSortScore;
	struct LONGMATRIX *pSortIndex;

	int nLastChrom;
	int nChromNum;
	int vChromDelimits[MAX_CHROMOSOME_NUM][2];
	int nD1;


	FILE *fpList;
	FILE *fpIn;
	FILE *fpOut;
	FILE *fpRefGeneOut;

	char strFileName[MED_LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	int ni,nj,nk,nx,ny,nz;
	int nLastStart;

	char strProbe[LINE_LENGTH];
	char strRefId[LINE_LENGTH];
	int nChr;
	int nStart,nEnd;
	char chStrand;
	int nP1,nP2;
	int nClassId;
	double dScore,dFDR;
	double dIdentity;
	double dMinScore;
	int nMinId;
	int nOverlap;
	int nMaxStart,nMinLen,nTempLen;

	/* ------------------------ */
	/* load refgene annotations */
	/* ------------------------ */

	/* get refgene number */
	fpIn = NULL;
	fpIn = fopen(strRefGenePath, "rt");
	if(fpIn == NULL)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, cannot open the refgene file!\n");
		exit(EXIT_FAILURE);
	}

	nRefGeneNum = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		nRefGeneNum++;
	}

	fclose(fpIn);

	if(nRefGeneNum <= 0)
	{
		printf("Warning: no refgenes available!\n");
		return PROC_SUCCESS;
	}

	/* init memory for storing refgene */
	vRefGene = NULL;
	vRefGene = (struct tagRefGene **)calloc(nRefGeneNum, sizeof(struct tagRefGene *));
	if(vRefGene == NULL)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, cannot create the memory for refgenes!\n");
		exit(EXIT_FAILURE);
	}

	/* load all refgene annotations */
	fpIn = NULL;
	fpIn = fopen(strRefGenePath, "rt");
	if(fpIn == NULL)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, cannot open the refgene file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	nChromNum = 0;
	nLastChrom = -MAX_CHROMOSOME_NUM;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		vRefGene[ni] = RefGeneCreate();
		RefGeneInit_FromGenomeLabFormat(vRefGene[ni], strLine, strSpecies);

		if(vRefGene[ni]->nChrom != nLastChrom)
		{
			if(nLastChrom != -MAX_CHROMOSOME_NUM)
			{
				vChromDelimits[nChromNum][0] = nD1;
				vChromDelimits[nChromNum][1] = ni-1;
				nChromNum++;
			}

			nLastChrom = vRefGene[ni]->nChrom;
			nD1 = ni;
		}

		ni++;
	}

	if(nLastChrom != -MAX_CHROMOSOME_NUM)
	{
		vChromDelimits[nChromNum][0] = nD1;
		vChromDelimits[nChromNum][1] = ni-1;
		nChromNum++;
	}

	fclose(fpIn);
	if(ni != nRefGeneNum)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, refgene number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* ---------------------------------------- */
	/* get active refgene and create intial map */
	/* ---------------------------------------- */
	
	/* get probe-refgene pair number */
	nProbeRefPairNum = 0;
	fpList = NULL;
	fpList = fopen(strProbeLists, "rt");
	if(fpList == NULL)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, cannot open probelists!\n");
		exit(EXIT_FAILURE);
	}

	nClassNum = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpList) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%d %s", &nClassId, strFileName);
		StrTrimLeft(strFileName);
		
		nClassNum++;

		fpIn = NULL;
		fpIn = fopen(strFileName, "rt");
		if(fpIn == NULL)
		{
			printf("Error: Affy_CreateProbeRefGeneMap, cannot open probelists!\n");
			exit(EXIT_FAILURE);
		}

		while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;

			nProbeRefPairNum++;
		}

		fclose(fpIn);
	}

	fclose(fpList);

	/* init storage space */
	vProbeRefPair = NULL;
	vProbeRefPair = (struct tagString **)calloc(nProbeRefPairNum, sizeof(struct tagString *));
	if(vProbeRefPair == NULL)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, cannot create enough storage space for probe-refgene pair!\n");
		exit(EXIT_FAILURE);
	}

	pClassVec = NULL;
	pClassVec = CreateIntMatrix(1, nProbeRefPairNum);
	if(pClassVec == NULL)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, cannot create enough storage space for class label!\n");
		exit(EXIT_FAILURE);
	}

	pProbeVec = NULL;
	pProbeVec = CreateDoubleMatrix(1, nProbeRefPairNum);
	if(pProbeVec == NULL)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, cannot create enough storage space for probe-refgene pair!\n");
		exit(EXIT_FAILURE);
	}
	
	pRefGeneVec = NULL;
	pRefGeneVec = CreateDoubleMatrix(1, nProbeRefPairNum);
	if(pRefGeneVec == NULL)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, cannot create enough storage space for probe-refgene pair!\n");
		exit(EXIT_FAILURE);
	}

	pScoreVec = NULL;
	pScoreVec = CreateDoubleMatrix(1, nProbeRefPairNum);
	if(pScoreVec == NULL)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, cannot create enough storage space for probe-refgene pair!\n");
		exit(EXIT_FAILURE);
	}


	pIdentityVec = NULL;
	pIdentityVec = CreateDoubleMatrix(1, nProbeRefPairNum);
	if(pIdentityVec == NULL)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, cannot create enough storage space for probe-refgene pair!\n");
		exit(EXIT_FAILURE);
	}

	pFilterVec = NULL;
	pFilterVec = CreateByteMatrix(1, nProbeRefPairNum);
	if(pFilterVec == NULL)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, cannot create enough storage space for probe-refgene pair!\n");
		exit(EXIT_FAILURE);
	}

	pClassIndicator = NULL;
	pClassIndicator = CreateIntMatrix(1, nClassNum);
	if(pClassIndicator == NULL)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, cannot create enough storage space for probe-refgene pair!\n");
		exit(EXIT_FAILURE);
	}

	pFinalVec = NULL;
	pFinalVec = CreateByteMatrix(1, nProbeRefPairNum);
	if(pFinalVec == NULL)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, cannot create enough storage space for probe-refgene pair!\n");
		exit(EXIT_FAILURE);
	}

	/* create initial map */
	ni = 0;
	fpList = NULL;
	fpList = fopen(strProbeLists, "rt");
	if(fpList == NULL)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, cannot open probelists!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpList) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%d %s", &nClassId, strFileName);
		StrTrimLeft(strFileName);
		
		fpIn = NULL;
		fpIn = fopen(strFileName, "rt");
		if(fpIn == NULL)
		{
			printf("Error: Affy_CreateProbeRefGeneMap, cannot open probelists!\n");
			exit(EXIT_FAILURE);
		}

		while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;

			sscanf(strLine, "%s %d %lf %lf %d %d %d %c %lf %s ", 
				strProbe, &nClassId,
				&dScore, &dFDR, 
				&nChr, &nStart, &nEnd, 
				&chStrand, &dIdentity, strRefId);

			pClassVec->pMatElement[ni] = nClassId;
			pProbeVec->pMatElement[ni] = (double)ni;
			pRefGeneVec->pMatElement[ni] = -1.0;
			pScoreVec->pMatElement[ni] = dScore;
			StringAddTail((vProbeRefPair+ni), strProbe);
			pFilterVec->pMatElement[ni] = 0;

			if( (nChr <= 0) || (nChr > nChromNum) || 
				(dIdentity < 80.0) || (strcmp(strRefId, "NA") == 0))
			{
				ni++;
				continue;
			}

			nChr = nChr-1;
			nP1 = vChromDelimits[nChr][0];
			nP2 = vChromDelimits[nChr][1];

			for(nk=nP1; nk<=nP2; nk++)
			{
				if(strcmp(strRefId, vRefGene[nk]->strName) == 0)
				{
					if( (chStrand == vRefGene[nk]->chStrand) && 
						((nEnd-vRefGene[nk]->nTxStart) > -MAX_MATCH_GAP) &&
						((vRefGene[nk]->nTxEnd-nStart) > -MAX_MATCH_GAP) )
					{
						pRefGeneVec->pMatElement[ni] = nk;
						pFilterVec->pMatElement[ni] = 1;
						break;
					}
				}
			}

			ni++;
		}

		fclose(fpIn);
	}

	fclose(fpList);

	if(ni != nProbeRefPairNum)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, probe-refgene pair number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* ---------------------------------------- */
	/* get unique refgene list and save map     */
	/* ---------------------------------------- */
	
	/* sort refgenes */
	pSortScore = NULL;
	pSortIndex = NULL;
	DMSORTMERGEA_0(pRefGeneVec, &pSortScore, &pSortIndex);
	
	/* determine classes of refgenes */
	nj = -1;
	for(ni=0; ni<nProbeRefPairNum; ni++)
	{
		if((int)(pSortScore->pMatElement[ni]) == -1)
			continue;

		if((int)(pSortScore->pMatElement[ni]) != nj)
		{
			if(nj != -1)
			{
				/* judge a refgene */
				ny = 0;
				
				for(nx=0; nx<nClassNum; nx++)
				{
					if(pClassIndicator->pMatElement[nx] > 0)
					{
						ny++;
					}
				}

				/* if the refgene cannot be consistently resolved, discard it */
				if(ny > 1)
				{
					for(nz=nLastStart; nz<ni; nz++)
					{
						nk = pSortIndex->pMatElement[nz];
						
						/* for debug */
						if((int)(pRefGeneVec->pMatElement[nk]) != nj)
						{
							printf("Error: logic error!\n");
							exit(EXIT_FAILURE);
						}

						pFilterVec->pMatElement[nk] = 0;
					}
				}
				else if(ny == 1)
				{
					dMinScore = 1e6;
					nMinId = pSortIndex->pMatElement[ni-1];
					for(nz=nLastStart; nz<ni; nz++)
					{
						nk = pSortIndex->pMatElement[nz];
						if(pFilterVec->pMatElement[nk] == 1)
						{
							if(pScoreVec->pMatElement[nk] < dMinScore)
							{
								dMinScore = pScoreVec->pMatElement[nk];
								nMinId = nk;
							}
						}
					}

					pFinalVec->pMatElement[nMinId] = 1;
				}
				else
				{
					
				}

				for(nx=0; nx<nClassNum; nx++)
				{
					pClassIndicator->pMatElement[nx] = 0;
				}
			}

			nLastStart = ni;
			nj = (int)(pSortScore->pMatElement[ni]);
			nk = pSortIndex->pMatElement[ni];
			if(pFilterVec->pMatElement[nk] == 1)
			{
				nClassId = pClassVec->pMatElement[nk];
				pClassIndicator->pMatElement[nClassId] = pClassIndicator->pMatElement[nClassId]+1;
			}
		}
		else
		{
			nk = pSortIndex->pMatElement[ni];
			nClassId = pClassVec->pMatElement[nk];
			if(pFilterVec->pMatElement[nk] == 1)
			{
				nClassId = pClassVec->pMatElement[nk];
				pClassIndicator->pMatElement[nClassId] = pClassIndicator->pMatElement[nClassId]+1;
			}
		}
	}

	if(nj != -1)
	{
		/* judge a refgene */
		ny = 0;
		
		for(nx=0; nx<nClassNum; nx++)
		{
			if(pClassIndicator->pMatElement[nx] > 0)
			{
				ny++;
			}
		}

		/* if the refgene cannot be consistently resolved, discard it */
		if(ny > 1)
		{
			for(nz=nLastStart; nz<ni; nz++)
			{
				nk = pSortIndex->pMatElement[nz];
				
				/* for debug */
				if((int)(pRefGeneVec->pMatElement[nk]) != nj)
				{
					printf("Error: logic error!\n");
					exit(EXIT_FAILURE);
				}

				pFilterVec->pMatElement[nk] = 0;
			}
		}
		else if(ny == 1)
		{
			dMinScore = 1e6;
			nMinId = pSortIndex->pMatElement[ni-1];
			for(nz=nLastStart; nz<ni; nz++)
			{
				nk = pSortIndex->pMatElement[nz];
				if(pFilterVec->pMatElement[nk] == 1)
				{
					if(pScoreVec->pMatElement[nk] < dMinScore)
					{
						dMinScore = pScoreVec->pMatElement[nk];
						nMinId = nk;
					}
				}
			}

			pFinalVec->pMatElement[nMinId] = 1;
		}
		else
		{
			
		}

		for(nx=0; nx<nClassNum; nx++)
		{
			pClassIndicator->pMatElement[nx] = 0;
		}
	}

	/* output refgene->probeset map */
	sprintf(strFileName, "%s_ref2probe.txt", strOutFile);
	fpOut = NULL;
	fpOut = fopen(strFileName, "wt");
	if(fpOut == NULL)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	sprintf(strFileName, "%s_refgene.txt", strOutFile);
	fpRefGeneOut = NULL;
	fpRefGeneOut = fopen(strFileName, "wt");
	if(fpRefGeneOut == NULL)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nProbeRefPairNum; ni++)
	{
		nk = pSortIndex->pMatElement[ni];
		if(pFinalVec->pMatElement[nk] == 1)
		{
			/* output */
			nj = (int)(pRefGeneVec->pMatElement[nk]);
			fprintf(fpOut, "%s\t%d\t%s\t%d\t%d\t%9.7e\n",
				vRefGene[nj]->strName, nj, vProbeRefPair[nk]->m_pString, nk, 
				pClassVec->pMatElement[nk], pScoreVec->pMatElement[nk]);
			RefGeneWrite(vRefGene[nj], fpRefGeneOut);
		}
	}
	
	fclose(fpOut);
	fclose(fpRefGeneOut);

	/* output probeset->refgene map */
	sprintf(strFileName, "%s_probe2ref.txt", strOutFile);
	fpOut = NULL;
	fpOut = fopen(strFileName, "wt");
	if(fpOut == NULL)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nProbeRefPairNum; ni++)
	{
		if(pFilterVec->pMatElement[ni] == 1)
		{
			/* output */
			nj = (int)(pRefGeneVec->pMatElement[ni]);
			fprintf(fpOut, "%s\t%d\t%s\t%d\t%d\t%9.7e\n",
				vProbeRefPair[ni]->m_pString, ni, vRefGene[nj]->strName, nj, 
				pClassVec->pMatElement[ni], pScoreVec->pMatElement[ni]);
		}
	}
	
	fclose(fpOut);

	/* output nonredundant refgene map */
	sprintf(strFileName, "%s_uniref.txt", strOutFile);
	fpRefGeneOut = NULL;
	fpRefGeneOut = fopen(strFileName, "wt");
	if(fpRefGeneOut == NULL)
	{
		printf("Error: Affy_CreateProbeRefGeneMap, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	nLastStart = -1;
	for(ni=0; ni<nProbeRefPairNum; ni++)
	{
		nk = pSortIndex->pMatElement[ni];
		if(pFinalVec->pMatElement[nk] == 1)
		{
			/* output */
			nj = (int)(pRefGeneVec->pMatElement[nk]);
			if(nLastStart == -1)
			{
				fprintf(fpRefGeneOut, "%s\t%d\t%s\t%d\t%d\t%9.7e\t",
					vRefGene[nj]->strName, nj, vProbeRefPair[nk]->m_pString, nk, 
					pClassVec->pMatElement[nk], pScoreVec->pMatElement[nk]);
				RefGeneWrite(vRefGene[nj], fpRefGeneOut);
				nLastStart = nj;
			}
			else
			{
				if(vRefGene[nLastStart]->nChrom == vRefGene[nj]->nChrom)
				{
					nOverlap = 0;
					nMaxStart = vRefGene[nj]->nTxStart;
					if(vRefGene[nLastStart]->nTxStart > nMaxStart)
						nMaxStart = vRefGene[nLastStart]->nTxStart;
					if(vRefGene[nj]->nTxEnd < vRefGene[nLastStart]->nTxEnd)
					{
						nOverlap = vRefGene[nj]->nTxEnd-nMaxStart+1;
					}
					else
					{
						nOverlap = vRefGene[nLastStart]->nTxEnd-nMaxStart+1;
					}

					nMinLen = vRefGene[nj]->nTxEnd-vRefGene[nj]->nTxStart+1;
					nTempLen = vRefGene[nLastStart]->nTxEnd-vRefGene[nLastStart]->nTxStart+1;
					if(nTempLen < nMinLen)
						nMinLen = nTempLen;

					if( (nOverlap <=0) )
					{
						nLastStart = nj;
						fprintf(fpRefGeneOut, "%s\t%d\t%s\t%d\t%d\t%9.7e\t",
							vRefGene[nj]->strName, nj, vProbeRefPair[nk]->m_pString, nk, 
							pClassVec->pMatElement[nk], pScoreVec->pMatElement[nk]);
						RefGeneWrite(vRefGene[nj], fpRefGeneOut);
					}
				}
				else
				{
					nLastStart = nj;
					fprintf(fpRefGeneOut, "%s\t%d\t%s\t%d\t%d\t%9.7e\t",
						vRefGene[nj]->strName, nj, vProbeRefPair[nk]->m_pString, nk, 
						pClassVec->pMatElement[nk], pScoreVec->pMatElement[nk]);
					RefGeneWrite(vRefGene[nj], fpRefGeneOut);
				}
			}
		}
	}
	
	fclose(fpRefGeneOut);
		
	/* ---------------------------------------- */
	/* release memory                           */
	/* ---------------------------------------- */
	for(ni=0; ni<nRefGeneNum; ni++)
	{
		RefGeneDestroy(vRefGene[ni]);
		vRefGene[ni] = NULL;
	}
	free(vRefGene);
	
	for(ni=0; ni<nProbeRefPairNum; ni++)
	{
		DeleteString(vProbeRefPair[ni]);
		vProbeRefPair[ni] = NULL;
	}
	free(vProbeRefPair);

	DestroyIntMatrix(pClassVec);
	DestroyDoubleMatrix(pProbeVec);
	DestroyDoubleMatrix(pRefGeneVec);
	DestroyDoubleMatrix(pScoreVec);
	DestroyDoubleMatrix(pIdentityVec);
	DestroyByteMatrix(pFilterVec);
	DestroyByteMatrix(pFinalVec);
	DestroyIntMatrix(pClassIndicator);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_LinkRefGeneScoreBackToProbeset()                                  */
/*  map refseq back to probeset.                                           */
/* ----------------------------------------------------------------------- */ 
int Affy_LinkRefGeneScoreBackToProbeset(char strProbeLists[], char strProbeRefGeneMapFile[],
										char strInFile[], char strOutFile[])
{
	/* define */
	FILE *fpList;
	FILE *fpIn;
	FILE *fpOut;
	struct tagString **vProbeName;
	struct tagString **vRefName;
	struct tagString **vAnnotation;
	int nProbeRefPairNum;
	char strFileName[LINE_LENGTH];
	int ni,nj,nk,nClassNum;
	char strLine[LONG_LINE_LENGTH];
	char strProbe[LINE_LENGTH];
	int nClassId;
	double dScore, dFDR;
	int nChr, nStart, nEnd, nLocusLink;
	char chStrand;
	double dIdentity;
	char strRefId[LINE_LENGTH];
	char strAnnotation[LONG_LINE_LENGTH];


	/* get proberef pair number */
	nProbeRefPairNum = 0;
	fpList = NULL;
	fpList = fopen(strProbeLists, "rt");
	if(fpList == NULL)
	{
		printf("Error: Affy_LinkRefGeneScoreBackToProbeset, cannot open probelists!\n");
		exit(EXIT_FAILURE);
	}

	nClassNum = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpList) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%d %s", &nClassId, strFileName);
		StrTrimLeft(strFileName);
		
		nClassNum++;

		fpIn = NULL;
		fpIn = fopen(strFileName, "rt");
		if(fpIn == NULL)
		{
			printf("Error: Affy_LinkRefGeneScoreBackToProbeset, cannot open probelists!\n");
			exit(EXIT_FAILURE);
		}

		while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;

			nProbeRefPairNum++;
		}

		fclose(fpIn);
	}

	fclose(fpList);

	/* prepare space */
	vProbeName = NULL;
	vProbeName = (struct tagString **)calloc(nProbeRefPairNum, sizeof(struct tagString *));
	if(vProbeName == NULL)
	{
		printf("Error: Affy_LinkRefGeneScoreBackToProbeset, cannot create space for loading string!\n");
		exit(EXIT_FAILURE);
	}

	vRefName = NULL;
	vRefName = (struct tagString **)calloc(nProbeRefPairNum, sizeof(struct tagString *));
	if(vRefName == NULL)
	{
		printf("Error: Affy_LinkRefGeneScoreBackToProbeset, cannot create space for loading string!\n");
		exit(EXIT_FAILURE);
	}

	vAnnotation = NULL;
	vAnnotation = (struct tagString **)calloc(nProbeRefPairNum, sizeof(struct tagString *));
	if(vAnnotation == NULL)
	{
		printf("Error: Affy_LinkRefGeneScoreBackToProbeset, cannot create space for loading string!\n");
		exit(EXIT_FAILURE);
	}


	/* load annotation */
	ni = 0;
	fpList = NULL;
	fpList = fopen(strProbeLists, "rt");
	if(fpList == NULL)
	{
		printf("Error: Affy_LinkRefGeneScoreBackToProbeset, cannot open probelists!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpList) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%d %s", &nClassId, strFileName);
		StrTrimLeft(strFileName);
		
		fpIn = NULL;
		fpIn = fopen(strFileName, "rt");
		if(fpIn == NULL)
		{
			printf("Error: Affy_LinkRefGeneScoreBackToProbeset, cannot open probelists!\n");
			exit(EXIT_FAILURE);
		}

		while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;

			sscanf(strLine, "%s %d %lf %lf %d %d %d %c %lf %s %d %s", 
				strProbe, &nClassId,
				&dScore, &dFDR, 
				&nChr, &nStart, &nEnd, 
				&chStrand, &dIdentity, strRefId, 
				&nLocusLink, strAnnotation);

			StringAddTail((vProbeName+ni), strProbe);
			StringAddTail((vRefName+ni), strRefId);
			StringAddTail((vAnnotation+ni), strAnnotation);

			ni++;
		}

		fclose(fpIn);
	}

	fclose(fpList);

	if(ni != nProbeRefPairNum)
	{
		printf("Error: Affy_LinkRefGeneScoreBackToProbeset, probe number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* write the mapping result */
	fpIn = NULL;
	fpIn = fopen(strInFile, "rt");
	if(fpIn == NULL)
	{
		printf("Error: Affy_LinkRefGeneScoreBackToProbeset, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutFile, "wt");
	if(fpOut == NULL)
	{
		printf("Error: Affy_LinkRefGeneScoreBackToProbeset, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fgets(strLine, LONG_LINE_LENGTH, fpIn);
	/* fprintf(fpOut, "Probeset\tRefId\tGene\t%s", strLine);
	*/
	fprintf(fpOut, "Gene\t%s", strLine);

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%s %d %s %d ", strRefId, &nj, strProbe, &nk);
		/*fprintf(fpOut, "%s\t%s\t%s\t%s\n", vProbeName[nk]->m_pString, 
			vRefName[nk]->m_pString, vAnnotation[nk]->m_pString, strLine);
		*/
		fprintf(fpOut, "%s\t%s\n", vAnnotation[nk]->m_pString, strLine);
	}


	fclose(fpIn);
	fclose(fpOut);


	/* release memory */
	for(ni=0; ni<nProbeRefPairNum; ni++)
	{
		DeleteString(vProbeName[ni]);
		vProbeName[ni] = NULL;
		DeleteString(vRefName[ni]);
		vRefName[ni] = NULL;
		DeleteString(vAnnotation[ni]);
		vAnnotation[ni] = NULL;
	}
	free(vProbeName);
	free(vRefName);
	free(vAnnotation);

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  Affy_LoadBar_Intensity()                                               */
/*  load intensity data from *.bar file                                    */
/* ----------------------------------------------------------------------- */ 
int Affy_LoadBar_Intensity(char strBarFile[], char strOutPath[])
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	struct INTMATRIX *pPos;
	struct DOUBLEMATRIX *pPM;
	struct DOUBLEMATRIX *pMM;
	char strOutFile[LINE_LENGTH];


	char strFiletype[9];
	float fVersion;
	int nSeq;
	int nCol;
	int *vFieldTypes;
	int nTagValPairNum;
	int nTagLen;
	char *vTagName;
	int nTagNameLen;
	char *vTagValPairName;
	int nSeqLen;
	char *vSeqName;
	int nSeqVerLen;
	char *vSeqVersion;
	int nDataPoint;

	int ni,nj,nk;

	float fTemp;
	int nTemp;

	/* load */
	fpIn = NULL;
	fpIn = fopen(strBarFile, "rb");
	if(fpIn == NULL)
	{
		printf("Error: Affy_LoadBar_Intensity, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}

	/* load head */
	fread( strFiletype, sizeof(char), 8, fpIn );
	strFiletype[8] = '\0';
	AFFYBAR_READ_FLOAT(fpIn, &fVersion);
	AFFYBAR_READ_INT(fpIn, &nSeq);
	AFFYBAR_READ_INT(fpIn, &nCol);

		
	/* load field types */
	vFieldTypes = NULL;
	vFieldTypes = (int *)calloc(nCol, sizeof(int));
	if(vFieldTypes == NULL)
	{
		printf("Error: Affy_LoadBar_Intensity, cannot create memory for loading!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nCol; ni++)
	{
		AFFYBAR_READ_INT(fpIn, (vFieldTypes+ni));
	}

	/* load tag-value pair */
	AFFYBAR_READ_INT(fpIn, &nTagValPairNum);
	/* for(ni = 0; ni<nTagValPairNum; ni++) */
	for(ni = 0; ni<1; ni++)
	{
		AFFYBAR_READ_INT(fpIn, &nTagLen);
		vTagName = NULL;
		vTagName = (char *)calloc((nTagLen+1), sizeof(char));
		fread( vTagName, sizeof(char), nTagLen, fpIn);
		vTagName[nTagLen] = '\0';
		free(vTagName);
		AFFYBAR_READ_INT(fpIn, &nTagNameLen);
		vTagValPairName = NULL;
		vTagValPairName = (char *)calloc((nTagNameLen+1), sizeof(char));
		fread( vTagValPairName, sizeof(char), nTagNameLen, fpIn);
		vTagValPairName[nTagNameLen] = '\0';
		free(vTagValPairName);
	}

	for(ni=0; ni<nSeq; ni++)
	{
		AFFYBAR_READ_INT(fpIn, &nSeqLen);
		vSeqName = NULL;
		vSeqName = (char *)calloc((nSeqLen+1), sizeof(char));
		fread( vSeqName, sizeof(char), nSeqLen, fpIn);
		vSeqName[nSeqLen] = '\0'; 
		

		AFFYBAR_READ_INT(fpIn, &nSeqVerLen);
		vSeqVersion = NULL;
		vSeqVersion = (char *)calloc((nSeqVerLen+1), sizeof(char));
		fread( vSeqVersion, sizeof(char), nSeqVerLen, fpIn);
		vSeqVersion[nSeqVerLen] ='\0';
		
		AFFYBAR_READ_INT(fpIn, &nDataPoint);

		/* prepare storage space */
		pPos = NULL;
		pPos = CreateIntMatrix(1, nDataPoint);
		if(pPos == NULL)
		{
			printf("Error: Affy_LoadBar_Intensity, cannot create memory for storing position!\n");
			exit(EXIT_FAILURE);
		}

		pPM = NULL;
		pPM = CreateDoubleMatrix(1, nDataPoint);
		if(pPM == NULL)
		{
			printf("Error: Affy_LoadBar_Intensity, cannot create memory for storing intensity values!\n");
			exit(EXIT_FAILURE);
		}

		pMM = NULL;
		pMM = CreateDoubleMatrix(1, nDataPoint);
		if(pMM == NULL)
		{
			printf("Error: Affy_LoadBar_Intensity, cannot create memory for storing intensity values!\n");
			exit(EXIT_FAILURE);
		}

		/* read data point */
		for(nj=0; nj<nDataPoint; nj++)
		{
			for(nk=0; nk<nCol; nk++)
			{
				switch(vFieldTypes[nk])
				{
					case 1: AFFYBAR_READ_FLOAT(fpIn, &fTemp);
						break;
					case 2: AFFYBAR_READ_INT(fpIn, &nTemp);
						break;
					default: printf("Error: Affy_LoadBar_Intensity, loading error!\n");
						exit(EXIT_FAILURE);
				}

				switch(nk)
				{
					case 0: pPos->pMatElement[nj] = nTemp;
						break;
					case 1: pPM->pMatElement[nj] = fTemp;
						break;
					case 2: pMM->pMatElement[nj] = fTemp;
						break;
					default: printf("Error: Affy_LoadBar_Intensity, loading error!\n");
						exit(EXIT_FAILURE);
				}
			}
		}

		/* save */
		sprintf(strOutFile, "%s%s_intensity.txt", strOutPath, vSeqName);
		fpOut = NULL;
		fpOut = fopen(strOutFile, "wt");
		if(fpOut == NULL)
		{
			printf("Error: Affy_LoadBar_Intensity, cannot open output files!\n");
			exit(EXIT_FAILURE);
		}
		for(nj=0; nj<nDataPoint; nj++)
		{
			fprintf(fpOut, "%d\t%f\t%f\n", pPos->pMatElement[nj],
				pPM->pMatElement[nj], pMM->pMatElement[nj]);
		}

		fclose(fpOut);

		/* destroy */
		free(vSeqName);
		free(vSeqVersion);

		DestroyIntMatrix(pPos);
		DestroyDoubleMatrix(pPM);
		DestroyDoubleMatrix(pMM);
	}

	while(feof(fpIn) == 0)
	{
		fread(strFiletype, sizeof(char), 1, fpIn);
		printf("%c", strFiletype[0]);
	}

	/* close file */
	fclose(fpIn);

	/* destroy */
	free(vFieldTypes);

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  Affy_LoadBar_Intensity_Group()                                         */
/*  load a group of intensity data from *.bar file                         */
/* ----------------------------------------------------------------------- */ 
int Affy_LoadBar_Intensity_Group(char strDataPath[], char strBarFileList[], char strOutPath[])
{
	/* define */
	FILE *fpList;
	FILE *fpOut;

	FILE **vfpIn;
	struct tagString **vCondition;
	int nArrayNum;

	char strLine[LINE_LENGTH];
	char strAlias[LINE_LENGTH];
	char strFilePath[LINE_LENGTH];
	
	struct INTMATRIX *pPos;
	struct DOUBLEMATRIX *pPM;
	struct DOUBLEMATRIX *pMM;

	char strFiletype[9];
	float fVersion;
	int nSeq;
	int nCol;
	char strFiletype0[9];
	float fVersion0;
	int nSeq0;
	int nCol0;


	int *vFieldTypes;
	int nTagValPairNum;
	int nTagLen;
	char *vTagName;
	int nTagNameLen;
	char *vTagValPairName;
	
	int nSeqLen;
	char *vSeqName;
	int nSeqGroupNameLen;
	char *vSeqGroupName;
	
	int nSeqVerLen;
	char *vSeqVersion;

	struct tagString *pSeqName;
	struct tagString *pSeqGroupName;
	struct tagString *pSeqVersion;

	int nNameValPairNum;
	int nParamNameLen;
	char *vParamName;
	int nParamValueLen;
	char *vParamValue;

	int nDataPoint, nDataPoint0;

	int ni,nj,nk,na;

	float fTemp;
	double fdiff;
	int nTemp;

	/* load file list */
	sprintf(strLine, "%s%s", strDataPath, strBarFileList);
	fpList = NULL;
	fpList = fopen(strLine, "r");
	if(fpList == NULL)
	{
		printf("Error: Affy_LoadBar_Intensity_Group, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}

	fgets(strLine, LINE_LENGTH, fpList);
	StrTrimLeft(strLine);
	StrTrimRight(strLine);
	nArrayNum = atoi(strLine);
	if(nArrayNum <= 0)
	{
		fclose(fpList);
		printf("Affy_LoadBar_Intensity_Group, no arrays were loaded!\n");
		return PROC_SUCCESS;
	}

	vCondition = NULL;
	vCondition = (struct tagString **)calloc(nArrayNum, sizeof(struct tagString*));
	if(vCondition == NULL)
	{
		printf("Error: Affy_LoadBar_Intensity_Group, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}

	vfpIn = NULL;
	vfpIn = (FILE **)calloc(nArrayNum, sizeof(FILE *));
	if(vfpIn == NULL)
	{
		printf("Error: Affy_LoadBar_Intensity_Group, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LINE_LENGTH, fpList) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(ni >= nArrayNum)
		{
			printf("Error: Affy_LoadBar_Intensity_Group, array number not match!\n");
			exit(EXIT_FAILURE);
		}

		sscanf(strLine, "%s %s", strFilePath, strAlias);
		StringAddTail((vCondition+ni), strAlias);
		sprintf(strLine, "%s%s", strDataPath, strFilePath);
		vfpIn[ni] = fopen(strLine, "rb");
		if(vfpIn[ni] == NULL)
		{
			printf("Error: Affy_LoadBar_Intensity_Group, cannot open input file!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if(ni != nArrayNum)
	{
		printf("Error: Affy_LoadBar_Intensity_Group, array number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* close list file */
	fclose(fpList);

	/* load head */
	fread( strFiletype, sizeof(char), 8, vfpIn[0] );
	strFiletype[8] = '\0';
	AFFYBAR_READ_FLOAT(vfpIn[0], &fVersion);
	AFFYBAR_READ_INT(vfpIn[0], &nSeq);
	AFFYBAR_READ_INT(vfpIn[0], &nCol);

	for(ni=1; ni<nArrayNum; ni++)
	{
		fread( strFiletype0, sizeof(char), 8, vfpIn[ni] );
		strFiletype0[8] = '\0';
		AFFYBAR_READ_FLOAT(vfpIn[ni], &fVersion0);
		AFFYBAR_READ_INT(vfpIn[ni], &nSeq0);
		AFFYBAR_READ_INT(vfpIn[ni], &nCol0);

		if( (strcmp(strFiletype, strFiletype0) != 0) || (fabs(fVersion0-fVersion)>1e-6)
			|| (nSeq != nSeq0) || (nCol != nCol0) )
		{
			printf("Error: Affy_LoadBar_Intensity_Group, *.bar file not match!\n");
			exit(EXIT_FAILURE);
		}
	}
		
	/* load field types */
	vFieldTypes = NULL;
	vFieldTypes = (int *)calloc(nCol, sizeof(int));
	if(vFieldTypes == NULL)
	{
		printf("Error: Affy_LoadBar_Intensity_Group, cannot create memory for loading!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nArrayNum; ni++)
	{
		for(nj=0; nj<nCol; nj++)
		{
			AFFYBAR_READ_INT(vfpIn[ni], (vFieldTypes+nj));
		}
	}

	/* load tag-value pair */
	for(ni=0; ni<nArrayNum; ni++)
	{
		AFFYBAR_READ_INT(vfpIn[ni], &nTagValPairNum);
		for(nj = 0; nj<nTagValPairNum; nj++)
		/* for(nj=0; nj<1; nj++) */
		{
			AFFYBAR_READ_INT(vfpIn[ni], &nTagLen);
			vTagName = NULL;
			vTagName = (char *)calloc((nTagLen+1), sizeof(char));
			fread( vTagName, sizeof(char), nTagLen, vfpIn[ni]);
			vTagName[nTagLen] = '\0';
			free(vTagName);
			AFFYBAR_READ_INT(vfpIn[ni], &nTagNameLen);
			vTagValPairName = NULL;
			vTagValPairName = (char *)calloc((nTagNameLen+1), sizeof(char));
			fread( vTagValPairName, sizeof(char), nTagNameLen, vfpIn[ni]);
			vTagValPairName[nTagNameLen] = '\0';
			free(vTagValPairName);
		}
	}


	/* prepare storage space */
	pPos = NULL;
	pPos = CreateIntMatrix(1, nArrayNum);
	if(pPos == NULL)
	{
		printf("Error: Affy_LoadBar_Intensity_Group, cannot create memory for storing position!\n");
		exit(EXIT_FAILURE);
	}

	pPM = NULL;
	pPM = CreateDoubleMatrix(1, nArrayNum);
	if(pPM == NULL)
	{
		printf("Error: Affy_LoadBar_Intensity_Group, cannot create memory for storing intensity values!\n");
		exit(EXIT_FAILURE);
	}

	pMM = NULL;
	pMM = CreateDoubleMatrix(1, nArrayNum);
	if(pMM == NULL)
	{
		printf("Error: Affy_LoadBar_Intensity_Group, cannot create memory for storing intensity values!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutPath, "wt");
	if(fpOut == NULL)
	{
		printf("Error: Affy_LoadBar_Intensity_Group, cannot open output files!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "chromosome\tposition");
	for(ni=0; ni<nArrayNum; ni++)
	{
		fprintf(fpOut, "\t%s", vCondition[ni]->m_pString);
	}
	fprintf(fpOut, "\n");
	
	
	/* load value */
	for(ni=0; ni<nSeq; ni++)
	{
		nDataPoint = -1;
		pSeqName = NULL;
		pSeqGroupName = NULL;
		pSeqVersion = NULL;

		for(na=0; na<nArrayNum; na++)
		{
			AFFYBAR_READ_INT(vfpIn[na], &nSeqLen);
			vSeqName = NULL;
			vSeqName = (char *)calloc((nSeqLen+1), sizeof(char));
			fread( vSeqName, sizeof(char), nSeqLen, vfpIn[na]);
			vSeqName[nSeqLen] = '\0'; 

			if(pSeqName == NULL)
				StringAddTail(&pSeqName, vSeqName);
			else
			{
				if(strcmp(pSeqName->m_pString, vSeqName) != 0)
				{
					printf("Error: Affy_LoadBar_Intensity_Group, *.bar file sequence name not match!\n");
					exit(EXIT_FAILURE);
				}
			}
			free(vSeqName);


			if(fVersion>1.99)
			{
				AFFYBAR_READ_INT(vfpIn[na], &nSeqGroupNameLen);
				vSeqGroupName = NULL;
				vSeqGroupName = (char *)calloc((nSeqGroupNameLen+1), sizeof(char));
				fread( vSeqGroupName, sizeof(char), nSeqGroupNameLen, vfpIn[na]);
				vSeqGroupName[nSeqGroupNameLen] = '\0'; 

				if(pSeqGroupName == NULL)
					StringAddTail(&pSeqGroupName, vSeqGroupName);
				else
				{
					if(strcmp(pSeqGroupName->m_pString, vSeqGroupName) != 0)
					{
						printf("Error: Affy_LoadBar_Intensity_Group, *.bar file sequence name not match!\n");
						exit(EXIT_FAILURE);
					}
				}
				free(vSeqGroupName);
			}
		
			AFFYBAR_READ_INT(vfpIn[na], &nSeqVerLen);
			vSeqVersion = NULL;
			vSeqVersion = (char *)calloc((nSeqVerLen+1), sizeof(char));
			fread( vSeqVersion, sizeof(char), nSeqVerLen, vfpIn[na]);
			vSeqVersion[nSeqVerLen] ='\0';

			if(pSeqVersion == NULL)
				StringAddTail(&pSeqVersion, vSeqVersion);
			else
			{
				if(strcmp(pSeqVersion->m_pString, vSeqVersion) != 0)
				{
					printf("Error: Affy_LoadBar_Intensity_Group, *.bar file sequence version not match!\n");
					exit(EXIT_FAILURE);
				}
			}
			free(vSeqVersion);

			
			/* load tag-value pair */
			if(fVersion>1.99)
			{
				AFFYBAR_READ_INT(vfpIn[na], &nNameValPairNum);
				for(nj = 0; nj<nNameValPairNum; nj++)
				{
					AFFYBAR_READ_INT(vfpIn[na], &nParamNameLen);
					vParamName = NULL;
					vParamName = (char *)calloc((nParamNameLen+1), sizeof(char));
					fread( vParamName, sizeof(char), nParamNameLen, vfpIn[na]);
					vParamName[nParamNameLen] = '\0';
					free(vParamName);
					AFFYBAR_READ_INT(vfpIn[na], &nParamValueLen);
					vParamValue = NULL;
					vParamValue = (char *)calloc((nParamValueLen+1), sizeof(char));
					fread( vParamValue, sizeof(char), nParamValueLen, vfpIn[na]);
					vParamValue[nParamValueLen] = '\0';
					free(vParamValue);
				}
			}

			AFFYBAR_READ_INT(vfpIn[na], &nDataPoint0);
			if(nDataPoint < 0)
			{
				nDataPoint = nDataPoint0;
			}
			else
			{
				if(nDataPoint != nDataPoint0)
				{
					printf("Error: Affy_LoadBar_Intensity_Group, *.bar file number of data point not match!\n");
					exit(EXIT_FAILURE);
				}
			}
		}

		/* read data point */
		for(nj=0; nj<nDataPoint; nj++)
		{
			for(na=0; na<nArrayNum; na++)
			{
				for(nk=0; nk<nCol; nk++)
				{
					switch(vFieldTypes[nk])
					{
						case 1: AFFYBAR_READ_FLOAT(vfpIn[na], &fTemp);
							break;
						case 2: AFFYBAR_READ_INT(vfpIn[na], &nTemp);
							break;
						default: printf("Error: Affy_LoadBar_Intensity_Group, loading error!\n");
							exit(EXIT_FAILURE);
					}

					switch(nk)
					{
						case 0: pPos->pMatElement[na] = nTemp;
							break;
						case 1: pPM->pMatElement[na] = fTemp;
							break;
						case 2: pMM->pMatElement[na] = fTemp;
							break;
						default: printf("Error: Affy_LoadBar_Intensity_Group, loading error!\n");
							exit(EXIT_FAILURE);
					}
				}

				if(pPos->pMatElement[na] != pPos->pMatElement[0])
				{
					printf("Error: Affy_LoadBar_Intensity_Group, chromosomal location not match!\n");
					exit(EXIT_FAILURE);
				}
			}

			/* save */
			fprintf(fpOut, "%s\t%d", pSeqName->m_pString, pPos->pMatElement[0]);
			/* fprintf(fpOut, "%d", pPos->pMatElement[0]); */
			for(na=0; na<nArrayNum; na++)
			{
				/* fdiff = pPM->pMatElement[na]-pMM->pMatElement[na];
				if(fdiff < 1.0)
					fdiff = 1.0;
				fprintf(fpOut, "\t%9.7e", fdiff); */
				fprintf(fpOut, "\t%9.7e", pPM->pMatElement[na]);
			}
			fprintf(fpOut, "\n");

		}
		
		/* destroy */
		DeleteString(pSeqName);
		DeleteString(pSeqVersion);
		if(fVersion > 1.99)
		{
			DeleteString(pSeqGroupName);
		}
	}

	fclose(fpOut);

	for(ni=0; ni<nArrayNum; ni++)
	{
		printf("file %d (%s) ", ni+1, vCondition[ni]->m_pString);
		while(feof(vfpIn[ni]) == 0)
		{
			fread(strFiletype, sizeof(char), 1, vfpIn[ni]);
			printf("%c", strFiletype[0]);
		}
	}
	
	DestroyIntMatrix(pPos);
	DestroyDoubleMatrix(pPM);
	DestroyDoubleMatrix(pMM);


	/* close file */
	for(ni=0; ni<nArrayNum; ni++)
	{
		fclose(vfpIn[ni]);
		DeleteString(vCondition[ni]);
		vCondition[ni] = NULL;
	}
	free(vfpIn);
	free(vCondition);

	/* destroy */
	free(vFieldTypes);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_LoadBPMAP()                                                       */
/*  load bpmap data from *.bar file.                                       */
/* ----------------------------------------------------------------------- */ 
int Affy_LoadBPMAP(char strInFile[], char strOutFile[])
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;

	/* variables */
	char strFileType[9];
	float fVersion;
	unsigned long nSeqNum;

	/* seq info */
	struct tagString **vSeqName;
	struct INTMATRIX *vProbePairNum;
	unsigned long nSeqNameLen;
	unsigned long nProbePairNum;

	/* map info */
	struct INTMATRIX *vSeqID;
	unsigned long nSeqID;
	int nProbeNum;
	struct tagAffyBpMapUnit *pNewUnit,*pCUnit,*pPUnit;
	struct tagAffyBpMapUnit *pUnitList;

	/* count */
	int ni,nj;
	int nEndPos;
	int nIgnore;

	/* load */
	fpIn = NULL;
	fpIn = fopen(strInFile, "rb");
	if(fpIn == NULL)
	{
		printf("Error: Affy_LoadBPMAP, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}
	fpOut = NULL;
	fpOut = fopen(strOutFile, "wt");
	if(fpOut == NULL)
	{
		printf("Error: Affy_LoadBPMAP, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	/* load head */
	fread( strFileType, sizeof(char), 8, fpIn );
	strFileType[8] = '\0';
	AFFYBAR_READ_FLOAT(fpIn, &fVersion);
	AFFYBAR_READ_ULONG(fpIn, &nSeqNum);

	/* load sequence names */
	vSeqName = NULL;
	vSeqName = (struct tagString **)calloc(nSeqNum, sizeof(struct tagString *));
	if(vSeqName == NULL)
	{
		printf("Error: Affy_LoadBPMAP, cannot load sequence name!\n");
		exit(EXIT_FAILURE);
	}
	vProbePairNum = NULL;
	vProbePairNum = CreateIntMatrix(nSeqNum,1);
	if(vProbePairNum == NULL)
	{
		printf("Error: Affy_LoadBPMAP, cannot load probe pair number!\n");
		exit(EXIT_FAILURE);
	}
	
	for(ni=0; ni<(int)nSeqNum; ni++)
	{
		AFFYBAR_READ_ULONG(fpIn, &nSeqNameLen);
		vSeqName[ni] = CreateString(nSeqNameLen);
		fread( vSeqName[ni]->m_pString, sizeof(char), nSeqNameLen, fpIn );
		vSeqName[ni]->m_pString[nSeqNameLen] = '\0';
		AFFYBAR_READ_ULONG(fpIn, &nProbePairNum);
		vProbePairNum->pMatElement[ni] = (int)nProbePairNum;
	}

	/* load map */
	vSeqID = NULL;
	vSeqID = CreateIntMatrix(nSeqNum,1);
	if(vSeqID == NULL)
	{
		printf("Error: Affy_LoadBPMAP, cannot load sequence id!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<(int)nSeqNum; ni++)
	{
		/* load seq id */
		AFFYBAR_READ_ULONG(fpIn, &nSeqID);
		nProbeNum = vProbePairNum->pMatElement[ni];
		pUnitList = NULL;
		
		/* load probes */
		for(nj=0; nj<nProbeNum; nj++)
		{
			/* load new unit */
			pNewUnit = NULL;
			pNewUnit = AffyBpMapUnitCreate();
			AffyBpMapUnitLoad(pNewUnit, fpIn);
			
			/* add to unitlist */
			if(pUnitList == NULL)
			{
				pUnitList = pNewUnit;
			}
			else
			{
				nIgnore = 0;
				pCUnit = pUnitList;
				while(pCUnit != NULL)
				{
					/* if same position */
					if(pNewUnit->nPos == pCUnit->nPos)
					{
						pCUnit->nDepthNum += 1;
						nIgnore = 1;
					}
					pCUnit = pCUnit->pNext;
				}

				if(nIgnore == 1)
				{
					AffyBpMapUnitDestroy(pNewUnit);
				}
				else
				{
					pCUnit = pUnitList;
					pPUnit = NULL;
					while(pCUnit != NULL)
					{
						/* if local repeat */
						if(pNewUnit->nPos-pCUnit->nPos <= LOCAL_REPEAT_RANGE)
						{
							if(strcmp(pCUnit->strProbeSeq, pNewUnit->strProbeSeq) == 0)
							{
								pCUnit->nRepeatNum += 1;
								pNewUnit->nRepeatNum += 1;
							}
						}

						/* get next */
						pPUnit = pCUnit;
						pCUnit = pPUnit->pNext;
					}

					pPUnit->pNext = pNewUnit;
					nEndPos = pNewUnit->nPos;

					/* write old unit */
					while(pUnitList != NULL)
					{
						/* if local repeat */
						if(nEndPos-pUnitList->nPos > LOCAL_REPEAT_RANGE)
						{
							pCUnit = pUnitList;
							pUnitList = pCUnit->pNext;
							pCUnit->pNext = NULL;

							fprintf(fpOut, "%d\t%d\t%d\t%s\n", pCUnit->nPos, pCUnit->nDepthNum, pCUnit->nRepeatNum, pCUnit->strProbeSeq);
							AffyBpMapUnitDestroy(pCUnit);
						}
						else
						{
							break;
						}
					}
				}
			}			
		}

		/* clear unit list */
		while(pUnitList != NULL)
		{
			pCUnit = pUnitList;
			pUnitList = pCUnit->pNext;
			pCUnit->pNext = NULL;

			fprintf(fpOut, "%d\t%d\t%d\t%s\n", pCUnit->nPos, pCUnit->nDepthNum, pCUnit->nRepeatNum, pCUnit->strProbeSeq);
			AffyBpMapUnitDestroy(pCUnit);
		}
	}

	/* load tail if any */
	while(feof(fpIn) == 0)
	{
		fread(strFileType, sizeof(char), 1, fpIn);
		printf("%c", strFileType[0]);
	}

	/* clear memeory */
	for(ni=0; ni<(int)nSeqNum; ni++)
	{
		DeleteString(vSeqName[ni]);
		vSeqName[ni] = NULL;
	}
	free(vSeqName);
	DestroyIntMatrix(vProbePairNum);
	DestroyIntMatrix(vSeqID);

	/* close file */
	fclose(fpIn);
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_LoadBPMAP_TileMap()                                               */
/*  load bpmap data from *.bar file, tilemap format                        */
/* ----------------------------------------------------------------------- */ 
int Affy_LoadBPMAP_TileMap(char strInFile[], char strPosFile[], 
						   char strMaskFile[], int *pnTotalProbeNum)
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	FILE *fpPos;

	/* variables */
	char strFileType[9];
	float fVersion;
	unsigned long nSeqNum;

	/* seq info */
	struct tagString **vSeqName;
	struct INTMATRIX *vProbePairNum;
	unsigned long nSeqNameLen;
	unsigned long nProbePairNum;

	/* map info */
	struct INTMATRIX *vSeqID;
	unsigned long nSeqID;
	int nProbeNum;
	struct tagAffyBpMapUnit *pNewUnit,*pCUnit,*pPUnit;
	struct tagAffyBpMapUnit *pUnitList;

	/* count */
	int ni,nj;
	int nEndPos;
	int nIgnore;

	/* load */
	fpIn = NULL;
	fpIn = fopen(strInFile, "rb");
	if(fpIn == NULL)
	{
		printf("Error: Affy_LoadBPMAP, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}
	fpOut = NULL;
	fpOut = fopen(strMaskFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: Affy_LoadBPMAP, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	fpPos = NULL;
	fpPos = fopen(strPosFile, "w");
	if(fpPos == NULL)
	{
		printf("Error: Affy_LoadBPMAP, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	/* load head */
	fread( strFileType, sizeof(char), 8, fpIn );
	strFileType[8] = '\0';
	AFFYBAR_READ_FLOAT(fpIn, &fVersion);
	AFFYBAR_READ_ULONG(fpIn, &nSeqNum);

	/* load sequence names */
	vSeqName = NULL;
	vSeqName = (struct tagString **)calloc(nSeqNum, sizeof(struct tagString *));
	if(vSeqName == NULL)
	{
		printf("Error: Affy_LoadBPMAP, cannot load sequence name!\n");
		exit(EXIT_FAILURE);
	}
	vProbePairNum = NULL;
	vProbePairNum = CreateIntMatrix(nSeqNum,1);
	if(vProbePairNum == NULL)
	{
		printf("Error: Affy_LoadBPMAP, cannot load probe pair number!\n");
		exit(EXIT_FAILURE);
	}
	
	*pnTotalProbeNum = 0;
	for(ni=0; ni<(int)nSeqNum; ni++)
	{
		AFFYBAR_READ_ULONG(fpIn, &nSeqNameLen);
		vSeqName[ni] = CreateString(nSeqNameLen);
		fread( vSeqName[ni]->m_pString, sizeof(char), nSeqNameLen, fpIn );
		vSeqName[ni]->m_pString[nSeqNameLen] = '\0';
		AFFYBAR_READ_ULONG(fpIn, &nProbePairNum);
		vProbePairNum->pMatElement[ni] = (int)nProbePairNum;
		*pnTotalProbeNum += (int)nProbePairNum;
	}

	/* load map */
	vSeqID = NULL;
	vSeqID = CreateIntMatrix(nSeqNum,1);
	if(vSeqID == NULL)
	{
		printf("Error: Affy_LoadBPMAP, cannot load sequence id!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<(int)nSeqNum; ni++)
	{
		/* load seq id */
		AFFYBAR_READ_ULONG(fpIn, &nSeqID);
		nProbeNum = vProbePairNum->pMatElement[ni];
		pUnitList = NULL;
		
		/* load probes */
		for(nj=0; nj<nProbeNum; nj++)
		{
			/* load new unit */
			pNewUnit = NULL;
			pNewUnit = AffyBpMapUnitCreate();
			AffyBpMapUnitLoad(pNewUnit, fpIn);
			fprintf(fpPos, "%s\t%d\n", vSeqName[ni]->m_pString, pNewUnit->nPos);
			
			/* add to unitlist */
			if(pUnitList == NULL)
			{
				pUnitList = pNewUnit;
			}
			else
			{
				nIgnore = 0;
				pCUnit = pUnitList;
				while(pCUnit != NULL)
				{
					/* if same position */
					if(pNewUnit->nPos == pCUnit->nPos)
					{
						pCUnit->nDepthNum += 1;
						nIgnore = 1;
					}
					pCUnit = pCUnit->pNext;
				}

				if(nIgnore == 1)
				{
					AffyBpMapUnitDestroy(pNewUnit);
				}
				else
				{
					pCUnit = pUnitList;
					pPUnit = NULL;
					while(pCUnit != NULL)
					{
						/* if local repeat */
						if(pNewUnit->nPos-pCUnit->nPos <= LOCAL_REPEAT_RANGE)
						{
							if(strcmp(pCUnit->strProbeSeq, pNewUnit->strProbeSeq) == 0)
							{
								pCUnit->nRepeatNum += 1;
								pNewUnit->nRepeatNum += 1;
							}
						}

						/* get next */
						pPUnit = pCUnit;
						pCUnit = pPUnit->pNext;
					}

					pPUnit->pNext = pNewUnit;
					nEndPos = pNewUnit->nPos;

					/* write old unit */
					while(pUnitList != NULL)
					{
						/* if local repeat */
						if(nEndPos-pUnitList->nPos > LOCAL_REPEAT_RANGE)
						{
							pCUnit = pUnitList;
							pUnitList = pCUnit->pNext;
							pCUnit->pNext = NULL;

							fprintf(fpOut, "%s\t%d\t%d\t%d\t%s\n", vSeqName[ni]->m_pString, pCUnit->nPos, pCUnit->nDepthNum, pCUnit->nRepeatNum, pCUnit->strProbeSeq);
							AffyBpMapUnitDestroy(pCUnit);
						}
						else
						{
							break;
						}
					}
				}
			}			
		}

		/* clear unit list */
		while(pUnitList != NULL)
		{
			pCUnit = pUnitList;
			pUnitList = pCUnit->pNext;
			pCUnit->pNext = NULL;

			fprintf(fpOut, "%s\t%d\t%d\t%d\t%s\n", vSeqName[ni]->m_pString, pCUnit->nPos, pCUnit->nDepthNum, pCUnit->nRepeatNum, pCUnit->strProbeSeq);
			AffyBpMapUnitDestroy(pCUnit);
		}
	}

	/* load tail if any */
	while(feof(fpIn) == 0)
	{
		fread(strFileType, sizeof(char), 1, fpIn);
		printf("%c", strFileType[0]);
	}

	/* clear memeory */
	for(ni=0; ni<(int)nSeqNum; ni++)
	{
		DeleteString(vSeqName[ni]);
		vSeqName[ni] = NULL;
	}
	free(vSeqName);
	DestroyIntMatrix(vProbePairNum);
	DestroyIntMatrix(vSeqID);

	/* close file */
	fclose(fpIn);
	fclose(fpOut);
	fclose(fpPos);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  AffyBpMapUnitCreate()                                                  */
/*  create bpmap unit.                                                     */
/* ----------------------------------------------------------------------- */ 
struct tagAffyBpMapUnit *AffyBpMapUnitCreate()
{
	/* define */
	struct tagAffyBpMapUnit *pUnit;
	
	/* create */
	pUnit = NULL;
	pUnit = (struct tagAffyBpMapUnit *)calloc(1, sizeof(struct tagAffyBpMapUnit));
	if(pUnit == NULL)
	{
		printf("Error: AffyBpMapUnitCreate, cannot create bpmap structure!\n");
		exit(EXIT_FAILURE);
	}

	/* return */
	return pUnit;
}

/* ----------------------------------------------------------------------- */ 
/*  AffyBpMapUnitDestroy()                                                 */
/*  delete bpmap unit.                                                     */
/* ----------------------------------------------------------------------- */ 
int AffyBpMapUnitDestroy(struct tagAffyBpMapUnit *pUnit)
{
	/* destroy */
	if(pUnit != NULL)
	{
		pUnit->pNext = NULL;
		free(pUnit);
		pUnit = NULL;
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  AffyBpMapUnitLoad()                                                    */
/*  load bpmap unit.                                                       */
/* ----------------------------------------------------------------------- */ 
int AffyBpMapUnitLoad(struct tagAffyBpMapUnit *pUnit, FILE *fpIn)
{
	/* define */
	char strTemp[PROBE_BPMAP_LEN];
	int ni,nj,nk;
	int nProbeLen;
	char bChar;

	/* load */
	if((pUnit == NULL) || (fpIn == NULL))
	{
		printf("Error: cannot read bpmap from file!\n");
		exit(EXIT_FAILURE);
	}

	nProbeLen = TILE_PROBE_LEN-1;
	AFFYBAR_READ_ULONG(fpIn, &(pUnit->nPMX));
	AFFYBAR_READ_ULONG(fpIn, &(pUnit->nPMY));
	AFFYBAR_READ_ULONG(fpIn, &(pUnit->nMMX));
	AFFYBAR_READ_ULONG(fpIn, &(pUnit->nMMY));
	fread( &(pUnit->bProbeLen), sizeof(unsigned char), 1, fpIn );
	fread( strTemp, sizeof(char), PROBE_BPMAP_LEN, fpIn );
	nk = 0;
	for(ni=0; ni<PROBE_BPMAP_LEN; ni++)
	{
		for(nj=0; (nj<4 && nk<nProbeLen); nj++)
		{
			bChar = (strTemp[ni] >> (2*(3-nj))) & 0x03;
			switch(bChar)
			{
				case DNA_BASE_A:
				  pUnit->strProbeSeq[nk] = 'A';
				  break;
				case DNA_BASE_C:
				  pUnit->strProbeSeq[nk] = 'C';
				  break;
				case DNA_BASE_G:
				  pUnit->strProbeSeq[nk] = 'G';
				  break;
				case DNA_BASE_T:
				  pUnit->strProbeSeq[nk] = 'T';
				  break;
				default:
				  pUnit->strProbeSeq[nk] = 'N';
				  break;
			}
			nk++;
		}
	}
	pUnit->strProbeSeq[nk] = '\0';
	AFFYBAR_READ_FLOAT(fpIn, &(pUnit->fMatchScore));
	AFFYBAR_READ_ULONG(fpIn, &(pUnit->nPos));
	fread( &(pUnit->bStrand), sizeof(unsigned char), 1, fpIn );
	pUnit->nRepeatNum = 1;
	pUnit->nDepthNum = 1;


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  AffyBpMapUnitLoad_v3()                                                 */
/*  load bpmap unit.                                                       */
/* ----------------------------------------------------------------------- */ 
int AffyBpMapUnitLoad_v3(struct tagAffyBpMapUnit *pUnit, FILE *fpIn, int nPMType)
{
	/* define */
	char strTemp[PROBE_BPMAP_LEN];
	int ni,nj,nk;
	int nProbeLen;
	char bChar;

	/* load */
	if((pUnit == NULL) || (fpIn == NULL))
	{
		printf("Error: cannot read bpmap from file!\n");
		exit(EXIT_FAILURE);
	}

	nProbeLen = TILE_PROBE_LEN-1;
	AFFYBAR_READ_ULONG(fpIn, &(pUnit->nPMX));
	AFFYBAR_READ_ULONG(fpIn, &(pUnit->nPMY));
	if(nPMType == 0)
	{
		AFFYBAR_READ_ULONG(fpIn, &(pUnit->nMMX));
		AFFYBAR_READ_ULONG(fpIn, &(pUnit->nMMY));
	}
	fread( &(pUnit->bProbeLen), sizeof(unsigned char), 1, fpIn );
	fread( strTemp, sizeof(char), PROBE_BPMAP_LEN, fpIn );
	nk = 0;
	for(ni=0; ni<PROBE_BPMAP_LEN; ni++)
	{
		for(nj=0; (nj<4 && nk<nProbeLen); nj++)
		{
			bChar = (strTemp[ni] >> (2*(3-nj))) & 0x03;
			switch(bChar)
			{
				case DNA_BASE_A:
				  pUnit->strProbeSeq[nk] = 'A';
				  break;
				case DNA_BASE_C:
				  pUnit->strProbeSeq[nk] = 'C';
				  break;
				case DNA_BASE_G:
				  pUnit->strProbeSeq[nk] = 'G';
				  break;
				case DNA_BASE_T:
				  pUnit->strProbeSeq[nk] = 'T';
				  break;
				default:
				  pUnit->strProbeSeq[nk] = 'N';
				  break;
			}
			nk++;
		}
	}
	pUnit->strProbeSeq[nk] = '\0';
	AFFYBAR_READ_FLOAT(fpIn, &(pUnit->fMatchScore));
	AFFYBAR_READ_ULONG(fpIn, &(pUnit->nPos));
	fread( &(pUnit->bStrand), sizeof(unsigned char), 1, fpIn );
	pUnit->nRepeatNum = 1;
	pUnit->nDepthNum = 1;


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  AffyBpMapUnitLoad_v3m2()                                               */
/*  load bpmap unit.                                                       */
/* ----------------------------------------------------------------------- */ 
int AffyBpMapUnitLoad_v3m2(struct tagAffyBpMapUnit *pUnit, FILE *fpIn, 
						   int nPMType, int little_endian_machine)
{
	/* define */
	char strTemp[PROBE_BPMAP_LEN];
	int ni,nj,nk;
	int nProbeLen;
	char bChar;
	unsigned int nX;

	/* load */
	if((pUnit == NULL) || (fpIn == NULL))
	{
		printf("Error: cannot read bpmap from file!\n");
		exit(EXIT_FAILURE);
	}

	nProbeLen = TILE_PROBE_LEN-1;

	if(big_endian_fread(&nX, 4, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: AffyBpMapUnitLoad_v3m2, cannot load PMx.\n");
		exit(EXIT_FAILURE);
	}
	pUnit->nPMX = nX;
	if(big_endian_fread(&nX, 4, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: AffyBpMapUnitLoad_v3m2, cannot load PMy.\n");
		exit(EXIT_FAILURE);
	}
	pUnit->nPMY = nX;

	if(nPMType == 0)
	{
		if(big_endian_fread(&nX, 4, 1, fpIn, little_endian_machine) != 1)
		{
			printf("Error: AffyBpMapUnitLoad_v3m2, cannot load MMx.\n");
			exit(EXIT_FAILURE);
		}
		pUnit->nMMX = nX;
		if(big_endian_fread(&nX, 4, 1, fpIn, little_endian_machine) != 1)
		{
			printf("Error: AffyBpMapUnitLoad_v3m2, cannot load MMy.\n");
			exit(EXIT_FAILURE);
		}
		pUnit->nMMY = nX;
	}
	fread( &(pUnit->bProbeLen), sizeof(unsigned char), 1, fpIn );
	fread( strTemp, sizeof(char), PROBE_BPMAP_LEN, fpIn );
	nk = 0;
	for(ni=0; ni<PROBE_BPMAP_LEN; ni++)
	{
		for(nj=0; (nj<4 && nk<nProbeLen); nj++)
		{
			bChar = (strTemp[ni] >> (2*(3-nj))) & 0x03;
			switch(bChar)
			{
				case DNA_BASE_A:
				  pUnit->strProbeSeq[nk] = 'A';
				  break;
				case DNA_BASE_C:
				  pUnit->strProbeSeq[nk] = 'C';
				  break;
				case DNA_BASE_G:
				  pUnit->strProbeSeq[nk] = 'G';
				  break;
				case DNA_BASE_T:
				  pUnit->strProbeSeq[nk] = 'T';
				  break;
				default:
				  pUnit->strProbeSeq[nk] = 'N';
				  break;
			}
			nk++;
		}
	}
	pUnit->strProbeSeq[nk] = '\0';
	if(big_endian_fread(&(pUnit->fMatchScore), 4, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: AffyBpMapUnitLoad_v3m2, cannot load match score.\n");
		exit(EXIT_FAILURE);
	}

	if(big_endian_fread(&nX, 4, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: AffyBpMapUnitLoad_v3m2, cannot load position of the probe.\n");
		exit(EXIT_FAILURE);
	}
	pUnit->nPos = nX;
	fread( &(pUnit->bStrand), sizeof(unsigned char), 1, fpIn );
	pUnit->nRepeatNum = 1;
	pUnit->nDepthNum = 1;

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_BPMAPFilter_GTrans()                                              */
/*  Filter out repeat and bad probes with bpmap informaiton.               */
/* ----------------------------------------------------------------------- */ 
int Affy_BPMAPFilter_GTrans(char strInFile[], char strRefFile[], char strOutFile[])
{
	/* define */
	FILE *fpIn;
	FILE *fpRef;
	FILE *fpOut;

	/* load line */
	char strLine[LINE_LENGTH];
	int nRefPos,nPos;
	int nDepthNum,nRepeatNum;
	char strProbeSeq[LINE_LENGTH];
	double dValue;
	double dMean,dN;
	int nInFileEnd,nRefFileEnd;
	int nWritePos;
	
	/* open file */
	fpIn = NULL;
	fpIn = fopen(strInFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: Affy_BPMAPFilter_GTrans, cannot read input file!\n");
		exit(EXIT_FAILURE);
	}
	fpRef = NULL;
	fpRef = fopen(strRefFile, "r");
	if(fpRef == NULL)
	{
		printf("Error: Affy_BPMAPFilter_GTrans, cannot read reference file!\n");
		exit(EXIT_FAILURE);
	}
	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: Affy_BPMAPFilter_GTrans, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	/* read files */
	nInFileEnd = 1;
	nRefFileEnd = 1;
	nPos = 0;
	nRefPos = 0;
	dValue = 0.0;
	nDepthNum = 0;
	nRepeatNum = 0;
	
	while(fgets(strLine, LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if((strLine[0] != '#') && (strLine[0] != '\0'))
		{
			nInFileEnd = 0;
			break;
		}
	}
	
	if(nInFileEnd == 0)
	{
		sscanf(strLine, "%d %lf", &nPos, &dValue);
	}

	while(fgets(strLine, LINE_LENGTH, fpRef) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] != '\0')
		{
			nRefFileEnd = 0;
			break;
		}
	}
	
	if(nRefFileEnd == 0)
	{
		sscanf(strLine, "%d %d %d %s", &nRefPos, &nDepthNum, &nRepeatNum, strProbeSeq);
	}
	
	dMean = 0.0;
	dN = 0.0;
	while((nInFileEnd == 0) && (nRefFileEnd == 0))
	{
		/* if equal position */
		if(nPos == nRefPos)
		{
			dMean += dValue;
			dN += 1.0;

			nInFileEnd = 1;
			while(fgets(strLine, LINE_LENGTH, fpIn) != NULL)
			{
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if((strLine[0] != '#') && (strLine[0] != '\0'))
				{
					nInFileEnd = 0;
					break;
				}
			}
			
			if(nInFileEnd == 0)
			{
				sscanf(strLine, "%d %lf", &nPos, &dValue);
			}
		}
		/* if pos > refpos */
		else if(nPos > nRefPos)
		{
			if(dN > 0.0)
			{
				dMean /= dN;
				nWritePos = nRefPos;
				if(nRepeatNum <= 1)
				{
					fprintf(fpOut, "%d\t%lf\n", nRefPos, dMean);
				}
			}

			dMean = 0.0;
			dN = 0.0;

			nRefFileEnd = 1;
			while(fgets(strLine, LINE_LENGTH, fpRef) != NULL)
			{
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] != '\0')
				{
					nRefFileEnd = 0;
					break;
				}
			}
			
			if(nRefFileEnd == 0)
			{
				sscanf(strLine, "%d %d %d %s", &nRefPos, &nDepthNum, &nRepeatNum, strProbeSeq);
			}
		}
		else
		{
			nInFileEnd = 1;
			while(fgets(strLine, LINE_LENGTH, fpIn) != NULL)
			{
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if((strLine[0] != '#') && (strLine[0] != '\0'))
				{
					nInFileEnd = 0;
					break;
				}
			}
			
			if(nInFileEnd == 0)
			{
				sscanf(strLine, "%d %lf", &nPos, &dValue);
			}
		}
	}

	if( (nWritePos < nRefPos) && (dN > 0.0))
	{
		dMean /= dN;
		if(nRepeatNum <= 1)
		{
			fprintf(fpOut, "%d\t%lf\n", nRefPos, dMean);
		}
	}

	/* close file */
	fclose(fpIn);
	fclose(fpRef);
	fclose(fpOut);


	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  Affy_BPMAPFilter_ORI()                                                 */
/*  Filter out repeat and bad probes with bpmap informaiton.               */
/* ----------------------------------------------------------------------- */ 
int Affy_BPMAPFilter_ORI(char strInFile[], char strRefFile[], char strOutFile[])
{
	/* define */
	FILE *fpIn;
	FILE *fpRef;
	FILE *fpOut;

	/* load line */
	char strLine[LINE_LENGTH];
	int nRefPos,nPos;
	double dPos;
	int nDepthNum,nRepeatNum;
	char strProbeSeq[LINE_LENGTH];
	double dValue;
	double dMean,dN;
	int nInFileEnd,nRefFileEnd;
	int nWritePos;
	
	/* open file */
	fpIn = NULL;
	fpIn = fopen(strInFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: Affy_BPMAPFilter_GTrans, cannot read input file!\n");
		exit(EXIT_FAILURE);
	}
	fpRef = NULL;
	fpRef = fopen(strRefFile, "r");
	if(fpRef == NULL)
	{
		printf("Error: Affy_BPMAPFilter_GTrans, cannot read reference file!\n");
		exit(EXIT_FAILURE);
	}
	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: Affy_BPMAPFilter_GTrans, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	/* read files */
	nInFileEnd = 1;
	nRefFileEnd = 1;
	nPos = 0;
	nRefPos = 0;
	dValue = 0.0;
	nDepthNum = 0;
	nRepeatNum = 0;
	
	while(fgets(strLine, LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if((strLine[0] != '#') && (strLine[0] != '\0'))
		{
			nInFileEnd = 0;
			break;
		}
	}
	
	if(nInFileEnd == 0)
	{
		sscanf(strLine, "%lf %lf", &dPos, &dValue);
		nPos = (int)dPos;
	}

	while(fgets(strLine, LINE_LENGTH, fpRef) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] != '\0')
		{
			nRefFileEnd = 0;
			break;
		}
	}
	
	if(nRefFileEnd == 0)
	{
		sscanf(strLine, "%d %d %d %s", &nRefPos, &nDepthNum, &nRepeatNum, strProbeSeq);
	}
	
	dMean = 0.0;
	dN = 0.0;
	while((nInFileEnd == 0) && (nRefFileEnd == 0))
	{
		/* if equal position */
		if(nPos == nRefPos)
		{
			dMean += dValue;
			dN += 1.0;

			nInFileEnd = 1;
			while(fgets(strLine, LINE_LENGTH, fpIn) != NULL)
			{
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if((strLine[0] != '#') && (strLine[0] != '\0'))
				{
					nInFileEnd = 0;
					break;
				}
			}
			
			if(nInFileEnd == 0)
			{
				sscanf(strLine, "%lf %lf", &dPos, &dValue);
				nPos = (int)dPos;
			}
		}
		/* if pos > refpos */
		else if(nPos > nRefPos)
		{
			if(dN > 0.0)
			{
				if((int)dN != nDepthNum)
				{
					printf("Warning: Affy_BPMAPFilter_ORI, depth number not consistent!\n");
				}

				dMean /= dN;
				nWritePos = nRefPos;
				if(nRepeatNum <= 1)
				{
					fprintf(fpOut, "%d\t%lf\n", nRefPos, dMean);
				}
			}

			dMean = 0.0;
			dN = 0.0;

			nRefFileEnd = 1;
			while(fgets(strLine, LINE_LENGTH, fpRef) != NULL)
			{
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] != '\0')
				{
					nRefFileEnd = 0;
					break;
				}
			}
			
			if(nRefFileEnd == 0)
			{
				sscanf(strLine, "%d %d %d %s", &nRefPos, &nDepthNum, &nRepeatNum, strProbeSeq);
			}
		}
		else
		{
			nInFileEnd = 1;
			while(fgets(strLine, LINE_LENGTH, fpIn) != NULL)
			{
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if((strLine[0] != '#') && (strLine[0] != '\0'))
				{
					nInFileEnd = 0;
					break;
				}
			}
			
			if(nInFileEnd == 0)
			{
				sscanf(strLine, "%lf %lf", &dPos, &dValue);
				nPos = (int)dPos;
			}
		}
	}

	if( (nWritePos < nRefPos) && (dN > 0.0))
	{
		if((int)dN != nDepthNum)
		{
			printf("Warning: Affy_BPMAPFilter_ORI, depth number not consistent!\n");
		}
		dMean /= dN;
		if(nRepeatNum <= 1)
		{
			fprintf(fpOut, "%d\t%lf\n", nRefPos, dMean);
		}
	}

	/* close file */
	fclose(fpIn);
	fclose(fpRef);
	fclose(fpOut);


	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  Affy_BPMAPFilter_DATA()                                                */
/*  Filter out repeat and bad probes with bpmap informaiton.               */
/* ----------------------------------------------------------------------- */ 
int Affy_BPMAPFilter_DATA(char strInFile[], char strRefFile[], char strOutFile[])
{
	/* define */
	FILE *fpIn;
	FILE *fpRef;
	FILE *fpOut;

	/* load line */
	char strLine[LONG_LINE_LENGTH];
	char strDataLine[LONG_LINE_LENGTH];
	int nRefPos,nPos;
	double dPos;
	int nDepthNum,nRepeatNum;
	char strProbeSeq[LINE_LENGTH];
	double dValue;
	int nInFileEnd,nRefFileEnd;
	
	/* open file */
	fpIn = NULL;
	fpIn = fopen(strInFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: Affy_BPMAPFilter_DATA, cannot read input file!\n");
		exit(EXIT_FAILURE);
	}
	fpRef = NULL;
	fpRef = fopen(strRefFile, "r");
	if(fpRef == NULL)
	{
		printf("Error: Affy_BPMAPFilter_DATA, cannot read reference file!\n");
		exit(EXIT_FAILURE);
	}
	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: Affy_BPMAPFilter_DATA, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	/* read files */
	nInFileEnd = 1;
	nRefFileEnd = 1;
	nPos = 0;
	nRefPos = 0;
	dValue = 0.0;
	nDepthNum = 0;
	nRepeatNum = 0;
	
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if((strLine[0] != '#') && (strLine[0] != '\0'))
		{
			nInFileEnd = 0;
			break;
		}
	}
	
	if(nInFileEnd == 0)
	{
		sscanf(strLine, "%lf", &dPos);
		strcpy(strDataLine, strLine);
		nPos = (int)dPos;
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpRef) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] != '\0')
		{
			nRefFileEnd = 0;
			break;
		}
	}
	
	if(nRefFileEnd == 0)
	{
		sscanf(strLine, "%d %d %d %s", &nRefPos, &nDepthNum, &nRepeatNum, strProbeSeq);
	}
	
	while((nInFileEnd == 0) && (nRefFileEnd == 0))
	{
		if(nRefPos == nPos)
		{
			if(nRepeatNum <= 1)
			{
				fprintf(fpOut, "%s\n", strDataLine);
			}
			if(fgets(strLine, LONG_LINE_LENGTH, fpIn) == NULL)
			{
				nInFileEnd = 1;
				break;
			}
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			sscanf(strLine, "%lf", &dPos);
			strcpy(strDataLine, strLine);
			nPos = (int)dPos;

			if(fgets(strLine, LONG_LINE_LENGTH, fpRef) == NULL)
			{
				nRefFileEnd = 1;
				break;
			}
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			sscanf(strLine, "%d %d %d %s", &nRefPos, &nDepthNum, &nRepeatNum, strProbeSeq);
		}
		else if(nPos < nRefPos)
		{
			while(nPos < nRefPos)
			{
				if(fgets(strLine, LONG_LINE_LENGTH, fpIn) == NULL)
				{
					nInFileEnd = 1;
					break;
				}
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				sscanf(strLine, "%lf", &dPos);
				strcpy(strDataLine, strLine);
				nPos = (int)dPos;
			}
		}
		else
		{
			while(nRefPos < nPos)
			{
				if(fgets(strLine, LONG_LINE_LENGTH, fpRef) == NULL)
				{
					nRefFileEnd = 1;
					break;
				}
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				sscanf(strLine, "%d %d %d %s", &nRefPos, &nDepthNum, &nRepeatNum, strProbeSeq);
			}
		}
	}
		

	/* close file */
	fclose(fpIn);
	fclose(fpRef);
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}



/* ----------------------------------------------------------------------- */ 
/*  Functions for loading affy binary files.                               */
/* ----------------------------------------------------------------------- */ 
void AFFYBAR_READ_INT(FILE *fpIn, int *value)
{
	char str[INT_SIZE];
	char str2[INT_SIZE];
	int ni;
	
	fread(str, sizeof(char), INT_SIZE, fpIn);
	if(AFFY_IS_BIG_ENDIAN == 1)
	{
		for(ni=0; ni<INT_SIZE; ni++)
			str2[ni] = str[INT_SIZE-1-ni];
		memcpy(value, str2, INT_SIZE);
	}
	else
	{
		memcpy(value, str, INT_SIZE);
	}
}

void AFFYBAR_READ_ULONG(FILE *fpIn, unsigned long *value)
{
	char str[ULONG_SIZE];
	char str2[ULONG_SIZE];
	int ni;
	
	fread(str, sizeof(char), ULONG_SIZE, fpIn);
	if(AFFY_IS_BIG_ENDIAN == 1)
	{
		for(ni=0; ni<ULONG_SIZE; ni++)
			str2[ni] = str[ULONG_SIZE-1-ni];
		memcpy(value, str2, ULONG_SIZE);
	}
	else
	{
		memcpy(value, str, ULONG_SIZE);
	}
}

void AFFYBAR_READ_FLOAT(FILE *fpIn, float *value)
{
	char str[FLOAT_SIZE];
	char str2[FLOAT_SIZE];
	int ni;
	
	fread(str, sizeof(char), FLOAT_SIZE, fpIn);
	if(AFFY_IS_BIG_ENDIAN == 1)
	{
		for(ni=0; ni<FLOAT_SIZE; ni++)
			str2[ni] = str[FLOAT_SIZE-1-ni];
		memcpy(value, str2, FLOAT_SIZE);
	}
	else
	{
		memcpy(value, str, FLOAT_SIZE);
	}
}

void AFFYBAR_READ_VERSION(FILE *fpIn, float *value)
{
	size_t num_read;
	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));
	float version;
	int affy_bug = 0;
	int temp;

	num_read = big_endian_fread(&version, 4, 1, fpIn, little_endian_machine);

	if (1 != num_read) 
	{
		printf("Error reading bpmap file: no version!\n");
		exit(EXIT_FAILURE);
	 }

	 if(version != 3.0f && version != 1.0f && version != 2.0f) 
	 {
		reverse_buf((char*)&version, 4);
		temp = (int)version;
		reverse_buf((char*)&temp, 4);
		version = (float)temp;

		affy_bug = 1;
		if(version != 1.0f && version != 2.0f)
		{
			printf("Error reading bpmap file: bad version!\n");
			exit(EXIT_FAILURE);
		}
	 }

	 *value = version;


	/* old */
	/* char str[FLOAT_SIZE];
	char str2[FLOAT_SIZE];
	int x,y,ni;
	float z;
	
	fread(str, sizeof(char), FLOAT_SIZE, fpIn);
	
	x = *((int *)(str));
	z = *((float*)&x);
	x = (int)z;
	x = ((((x) & 0xff000000) >> 24) | (((x) & 0x00ff0000) >>  8) | (((x) & 0x0000ff00) <<  8) | (((x) & 0x000000ff) << 24));
	if( (x==1) || (x==2) || (x==3) )
	{
		AFFY_IS_BIG_ENDIAN = 1;
		*value = (float)x;	
	}
	else
	{
		x = *((int *)(str));
		x = ((((x) & 0xff000000) >> 24) | (((x) & 0x00ff0000) >>  8) | (((x) & 0x0000ff00) <<  8) | (((x) & 0x000000ff) << 24));
		y = (int)(*(float*)(&x));
		y = ((((y) & 0xff000000) >> 24) | (((y) & 0x00ff0000) >>  8) | (((y) & 0x0000ff00) <<  8) | (((y) & 0x000000ff) << 24));
		
		if( (y==1) || (y==2) || (y==3) )
		{
			AFFY_IS_BIG_ENDIAN = 0;
			*value = (float)y;
		}

		else
		{
			for(ni=0; ni<FLOAT_SIZE; ni++)
				str2[ni] = str[FLOAT_SIZE-1-ni];
			z = *((float *)str2);
			y = (int)z;
			
			if( (y==1) || (y==2) || (y==3) )
			{
				AFFY_IS_BIG_ENDIAN = 1;
				*value = (float)y;
			}

			else
			{
				z = *((float *)str);
				y = (int)z;
				
				if( (y==1) || (y==2) || (y==3) )
				{
					AFFY_IS_BIG_ENDIAN = 0;
					*value = (float)y;
				}
				else
				{
					printf("Error: AFFYBAR_READ_VERSION, cannot read verion information of *.bpmap file!\n");
					exit(EXIT_FAILURE);
				}
			}
		}
	}

	printf("AFFY_IS_BIG_ENDIAN = %d\n", AFFY_IS_BIG_ENDIAN); */
}

/* void affy_bug_fix(float &data)
{
	reverse_buf((char*)&data, 4);
	int temp = (int)data;
	reverse_buf((char*)&temp, 4);
	data = (float)temp;
}
*/

/* ----------------------------------------------------------------------- */ 
/*  Affy_LoadDatabase_Affy2Refid()                                         */
/*  Load affy probe id to refid map.                                       */
/* ----------------------------------------------------------------------- */ 
struct tagStringPair **Affy_LoadDatabase_Affy2Refid(char strDatabasePath[], int *pPairNum)
{
	/* define */
	FILE *fpIn;
	struct tagStringPair *pPair,*pCurrentPair;
	struct tagStringPair *pPairList;
	struct tagStringPair **vPair = NULL;
	int nPairNum = 0;
	char strLine[LONG_LINE_LENGTH];
	char strProbe[LINE_LENGTH];
	char strRefid[LINE_LENGTH];
	int ni;
	char *chp1,*chp2;

	/* load */
	fpIn = NULL;
	fpIn = fopen(strDatabasePath, "r");
	if(fpIn == NULL)
	{
		printf("Error: Affy_LoadDatabase_Affy2Refid, cannot open mapping database file!\n");
		exit(EXIT_FAILURE);
	}

	nPairNum = 0;
	pCurrentPair = NULL;
	pPairList = NULL;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		/* loading pairs */
		chp1 = strchr(strLine, '\t');
		if(chp1 == NULL)
			continue;

		*chp1 = '\0';
		strcpy(strProbe, strLine);
		StrTrimRight(strProbe);
		if( (strcmp(strProbe, "") == 0) || (strcmp(strProbe, "---") == 0) || (strcmp(strProbe, "NA") == 0) )
			continue;
		
		chp2 = chp1+1;
		chp1 = strstr(chp2, "///");
		while(chp1 != NULL)
		{
			*chp1 = '\0';
			strcpy(strRefid, chp2);
			StrTrimLeft(strRefid);
			StrTrimRight(strRefid);

			if( (strcmp(strRefid, "") == 0) || (strcmp(strRefid, "---") == 0)
				|| (strcmp(strRefid, "NA") == 0) )
			{
			}
			else
			{
				pPair = NULL;
				pPair = StringPairCreate();
				if(pPair == NULL)
				{
					printf("Error: Affy_LoadDatabase_Affy2Refid, cannot create string pair!\n");
					exit(EXIT_FAILURE);
				}
				StringAddTail( &(pPair->m_pStr1), strProbe);
				StringAddTail( &(pPair->m_pStr2), strRefid);
				if(pPairList == NULL)
				{
					pPairList = pPair;
					pCurrentPair = pPair;
				}
				else
				{
					pCurrentPair->m_pNext = pPair;
					pCurrentPair = pPair;
				}
				nPairNum++;
			}

			chp2 = chp1+3;
			chp1 = strstr(chp2, "///");
		}

		/* process the last one */
		strcpy(strRefid, chp2);
		StrTrimLeft(strRefid);
		StrTrimRight(strRefid);

		if( (strcmp(strRefid, "") == 0) || (strcmp(strRefid, "---") == 0)
			|| (strcmp(strRefid, "NA") == 0) )
		{
		}
		else
		{
			pPair = NULL;
			pPair = StringPairCreate();
			if(pPair == NULL)
			{
				printf("Error: Affy_LoadDatabase_Affy2Refid, cannot create string pair!\n");
				exit(EXIT_FAILURE);
			}
			StringAddTail( &(pPair->m_pStr1), strProbe);
			StringAddTail( &(pPair->m_pStr2), strRefid);
			if(pPairList == NULL)
			{
				pPairList = pPair;
				pCurrentPair = pPair;
			}
			else
			{
				pCurrentPair->m_pNext = pPair;
				pCurrentPair = pPair;
			}
			nPairNum++;
		}
	}
	fclose(fpIn);

	if(nPairNum <= 0)
	{
		printf("Warning: Affy_LoadDatabase_Affy2Refid, null database!\n");
		*pPairNum = 0;
		return NULL;
	}

	vPair = NULL;
	vPair = (struct tagStringPair **)calloc(nPairNum, sizeof(struct tagStringPair *));
	if(vPair == NULL)
	{
		printf("Error: Affy_LoadDatabase_Affy2Refid, cannot organize pairs into array!\n");
		exit(EXIT_FAILURE);
	}

	ni=0;
	while(pPairList != NULL)
	{
		pPair = pPairList;
		pPairList = pPair->m_pNext;
		pPair->m_pNext = NULL;
		vPair[ni] = pPair;
		ni++;
	}

	if(ni != nPairNum)
	{
		printf("Error: Affy_LoadDatabase_Affy2Refid, pair number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* return */
	*pPairNum = nPairNum;
	return vPair;
}


/* ----------------------------------------------------------------------- */ 
/*  Affy_LoadDatabase_ExpressionData()                                     */
/*  Load expression data for probe selection                               */
/* ----------------------------------------------------------------------- */ 
struct tagStringPair **Affy_LoadDatabase_ExpressionData(char strDatabasePath[], 
					int *pPairNum, char strHeader[])
{
	/* define */
	FILE *fpIn;
	struct tagStringPair *pPair,*pCurrentPair;
	struct tagStringPair *pPairList;
	struct tagStringPair **vPair = NULL;
	int nPairNum = 0;
	char strLine[LONG_LINE_LENGTH];
	char strProbe[LINE_LENGTH];
	int ni;

	/* load */
	fpIn = NULL;
	fpIn = fopen(strDatabasePath, "r");
	if(fpIn == NULL)
	{
		printf("Error: Affy_LoadDatabase_ExpressionData, cannot open mapping database file!\n");
		exit(EXIT_FAILURE);
	}

	if(fgets(strLine, LONG_LINE_LENGTH, fpIn) == NULL)
	{
		*pPairNum = 0;
		return NULL;
	}

	StrTrimLeft(strLine);
	StrTrimRight(strLine);
	strcpy(strHeader, strLine);

	nPairNum = 0;
	pCurrentPair = NULL;
	pPairList = NULL;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		/* loading pairs */
		sscanf(strLine, "%s", strProbe);
		if( (strcmp(strProbe, "") == 0) || (strcmp(strProbe, "---") == 0) || (strcmp(strProbe, "NA") == 0) )
			continue;
		
		pPair = NULL;
		pPair = StringPairCreate();
		if(pPair == NULL)
		{
			printf("Error: Affy_LoadDatabase_ExpressionData, cannot create string pair!\n");
			exit(EXIT_FAILURE);
		}
		StringAddTail( &(pPair->m_pStr1), strProbe);
		StringAddTail( &(pPair->m_pStr2), strLine);
		if(pPairList == NULL)
		{
			pPairList = pPair;
			pCurrentPair = pPair;
		}
		else
		{
			pCurrentPair->m_pNext = pPair;
			pCurrentPair = pPair;
		}
		nPairNum++;
	}

	
	fclose(fpIn);

	if(nPairNum <= 0)
	{
		printf("Warning: Affy_LoadDatabase_ExpressionData, null database!\n");
		*pPairNum = 0;
		return NULL;
	}

	vPair = NULL;
	vPair = (struct tagStringPair **)calloc(nPairNum, sizeof(struct tagStringPair *));
	if(vPair == NULL)
	{
		printf("Error: Affy_LoadDatabase_ExpressionData, cannot organize pairs into array!\n");
		exit(EXIT_FAILURE);
	}

	ni=0;
	while(pPairList != NULL)
	{
		pPair = pPairList;
		pPairList = pPair->m_pNext;
		pPair->m_pNext = NULL;
		vPair[ni] = pPair;
		ni++;
	}

	if(ni != nPairNum)
	{
		printf("Error: Affy_LoadDatabase_ExpressionData, pair number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* return */
	*pPairNum = nPairNum;
	return vPair;
}


/* ----------------------------------------------------------------------- */ 
/*  Affy_LoadCELv4()                                                       */
/*  Loading raw data from a single affymetrix's *.CEL (v4) file            */
/* ----------------------------------------------------------------------- */ 
struct tagCELData *Affy_LoadCELv4(char strFileName[])
{
	/* define */
	struct tagCELData *pCELData = NULL;

	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));
	
	FILE *fpCel;
	int nMagicnumber;
	int nVersionnumber;
	int nCols, nRows;
	int nNumberCells;

	int nHeaderLength;
	int nAlgorithmLength;
	int nAlgorithmParamLength;

	int nCellMargin;
	unsigned int nOutlierCells;
	unsigned int nMaskedCells;
	int nSubGrids;

	float fI[2];
	short int nCount;
	short int nEntry[2];
	int nPos[2];
	int ni;

	int nTotalX,nTotalY;
	
	/* init */
	nTotalX = 0;
	nTotalY = 0;

	/* load *.CEL */
	fpCel = NULL;
	fpCel = fopen(strFileName, "rb");
	if(fpCel == NULL)
	{
		printf("Error: Affy_LoadCELv4, cannot open *.CEL file!\n");
		exit(EXIT_FAILURE);
	}

	/* load magic number */
	if(little_endian_fread(&nMagicnumber, INT_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_LoadCELv4, cannot load magic number.\n");
		return NULL;
	}
	
	if(nMagicnumber != 64)
	{
		printf("Error: Affy_LoadCELv4, cannot load magic number correctly.\n");
		return NULL;
	}

	/* load version number */
	if(little_endian_fread(&nVersionnumber, INT_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_LoadCELv4, cannot load version number.\n");
		return NULL;
	}
	
	if(nVersionnumber != 4)
	{
		printf("Error: Affy_LoadCELv4, cannot load version number correctly.\n");
		return NULL;
	}

	/* load number of columns, rows and cells */
	if(little_endian_fread(&nCols, INT_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_LoadCELv4, cannot load number of columns.\n");
		return NULL;
	}
	if(little_endian_fread(&nRows, INT_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_LoadCELv4, cannot load number of rows.\n");
		return NULL;
	}
	if(little_endian_fread(&nNumberCells, INT_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_LoadCELv4, cannot load number of cells.\n");
		return NULL;
	}
	if(nNumberCells != nCols*nRows)
	{
		printf("Error: Affy_LoadCELv4, cannot load number of columns, rows and cells correctly.\n");
		return NULL;
	}

	/* create CELData object */
	pCELData = Affy_CELData_Create();
	if(pCELData == NULL)
	{
		printf("Error: Affy_LoadCELv4, cannot create CELData structure.\n");
		return NULL;
	}

	pCELData->nMagicnumber = nMagicnumber;
	pCELData->nVersionnumber = nVersionnumber;
	pCELData->nCols = nCols;
	pCELData->nRows = nRows;
	pCELData->nNumberCells = nNumberCells;
	pCELData->nTotalX = nCols;
	pCELData->nTotalY = nRows;

	/* load header info */
	if(little_endian_fread(&nHeaderLength, INT_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_LoadCELv4, cannot load header length.\n");
		exit(EXIT_FAILURE);
	}
	pCELData->nHeaderLength = nHeaderLength;
	pCELData->vHeader = NULL;
	pCELData->vHeader = CreateString(nHeaderLength);
	if( (int)little_endian_fread(pCELData->vHeader->m_pString, sizeof(char), nHeaderLength, fpCel,  little_endian_machine) != nHeaderLength)
	{
		printf("Error: Affy_LoadCELv4, cannot load header.\n");
		exit(EXIT_FAILURE);
	}
	pCELData->vHeader->m_pString[nHeaderLength] = '\0';
	
	/* load algorithm */
	if(little_endian_fread(&nAlgorithmLength, INT_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_LoadCELv4, cannot load algorithm length.\n");
		exit(EXIT_FAILURE);
	}
	pCELData->nAlgorithmLength = nAlgorithmLength;
	pCELData->vAlgorithm = NULL;
	pCELData->vAlgorithm = CreateString(nAlgorithmLength);
	if( (int)little_endian_fread(pCELData->vAlgorithm->m_pString, sizeof(char), nAlgorithmLength, fpCel, little_endian_machine) != nAlgorithmLength)
	{
		printf("Error: Affy_LoadCELv4, cannot load algorithm.\n");
		exit(EXIT_FAILURE);
	}
	pCELData->vAlgorithm->m_pString[nAlgorithmLength] = '\0';

	/* load algorithm parameters */
	if(little_endian_fread(&nAlgorithmParamLength, INT_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_LoadCELv4, cannot load algorithm parameter length.\n");
		exit(EXIT_FAILURE);
	}
	pCELData->nAlgorithmParamLength = nAlgorithmParamLength;
	pCELData->vAlgorithmParam = NULL;
	pCELData->vAlgorithmParam = CreateString(nAlgorithmParamLength);
	if( (int)little_endian_fread(pCELData->vAlgorithmParam->m_pString, sizeof(char), nAlgorithmParamLength, fpCel, little_endian_machine) != nAlgorithmParamLength)
	{
		printf("Error: Affy_LoadCELv4, cannot load algorithm parameters.\n");
		exit(EXIT_FAILURE);
	}
	pCELData->vAlgorithmParam->m_pString[nAlgorithmParamLength] = '\0';

	/* load cell margin */
	if(little_endian_fread(&nCellMargin, INT_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_LoadCELv4, cannot load cell margin.\n");
		exit(EXIT_FAILURE);
	}
	pCELData->nCellMargin = nCellMargin;

	/* load number of outlier cells */
	if(little_endian_fread(&nOutlierCells, DWORD_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_LoadCELv4, cannot load number of outlier cells.\n");
		exit(EXIT_FAILURE);
	}
	pCELData->nOutlierCells = nOutlierCells;

	/* load number of masked cells */
	if(little_endian_fread(&nMaskedCells, DWORD_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_LoadCELv4, cannot load number of masked cells.\n");
		exit(EXIT_FAILURE);
	}
	pCELData->nMaskedCells = nMaskedCells;

	/* load number of sub-grids */
	if(little_endian_fread(&nSubGrids, INT_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_LoadCELv4, cannot load number of sub-grids.\n");
		exit(EXIT_FAILURE);
	}
	pCELData->nSubGrids = nSubGrids;

	/* allocate memory */
	pCELData->pIntensity = NULL;
	pCELData->pIntensity = CreateDoubleMatrix(1, nNumberCells);
	if((pCELData->pIntensity == NULL) && (nNumberCells>0))
	{
		printf("Error: Affy_LoadCELv4, cannot allocate enough memory for loading intensities!\n");
		exit(EXIT_FAILURE);
	}

	pCELData->pSD = NULL;
	pCELData->pSD = CreateDoubleMatrix(1, nNumberCells);
	if((pCELData->pSD == NULL) && (nNumberCells>0))
	{
		printf("Error: Affy_LoadCELv4, cannot allocate enough memory for loading SD!\n");
		exit(EXIT_FAILURE);
	}

	pCELData->pPixelNum = NULL;
	pCELData->pPixelNum = CreateIntMatrix(1, nNumberCells);
	if((pCELData->pPixelNum == NULL) && (nNumberCells>0))
	{
		printf("Error: Affy_LoadCELv4, cannot allocate enough memory for loading pixel numbers!\n");
		exit(EXIT_FAILURE);
	}

	pCELData->pMaskedX = NULL;
	pCELData->pMaskedX = CreateIntMatrix(1, nMaskedCells);
	if((pCELData->pMaskedX == NULL) && (nMaskedCells>0))
	{
		printf("Error: Affy_LoadCELv4, cannot allocate enough memory for loading masked cells!\n");
		exit(EXIT_FAILURE);
	}

	pCELData->pMaskedY = NULL;
	pCELData->pMaskedY = CreateIntMatrix(1, nMaskedCells);
	if((pCELData->pMaskedY == NULL) && (nMaskedCells>0))
	{
		printf("Error: Affy_LoadCELv4, cannot allocate enough memory for loading masked cells!\n");
		exit(EXIT_FAILURE);
	}

	pCELData->pOutlierX = NULL;
	pCELData->pOutlierX = CreateIntMatrix(1, nOutlierCells);
	if((pCELData->pOutlierX == NULL) && (nOutlierCells>0))
	{
		printf("Error: Affy_LoadCELv4, cannot allocate enough memory for loading outlier cells!\n");
		exit(EXIT_FAILURE);
	}

	pCELData->pOutlierY = NULL;
	pCELData->pOutlierY = CreateIntMatrix(1, nOutlierCells);
	if((pCELData->pOutlierY == NULL) && (nOutlierCells>0))
	{
		printf("Error: Affy_LoadCELv4, cannot allocate enough memory for loading outlier cells!\n");
		exit(EXIT_FAILURE);
	}

	pCELData->nModifiedCells = 0;
	pCELData->pModifiedX = NULL;
	pCELData->pModifiedY = NULL;
	pCELData->pModifiedOrig = NULL;

	/* load intensity */
	for(ni=0; ni<nNumberCells; ni++)
	{
		if(little_endian_fread(fI, FLOAT_SIZE, 2, fpCel, little_endian_machine) != 2)
		{
			printf("Error: Affy_LoadCELv4, cannot load intensity values!\n");
			exit(EXIT_FAILURE);
		}
		if(little_endian_fread(&nCount, SHORT_SIZE, 1, fpCel, little_endian_machine) != 1)
		{
			printf("Error: Affy_LoadCELv4, cannot load pixel count!\n");
			exit(EXIT_FAILURE);
		}

		pCELData->pIntensity->pMatElement[ni] = fI[0];
		pCELData->pSD->pMatElement[ni] = fI[1];
		pCELData->pPixelNum->pMatElement[ni] = nCount;
	}

	/* load masked entries */
	for(ni=0; ni<(int)nMaskedCells; ni++)
	{
		if(little_endian_fread(nEntry, SHORT_SIZE, 2, fpCel, little_endian_machine) != 2)
		{
			printf("Error: Affy_LoadCELv4, cannot load masked entries!\n");
			exit(EXIT_FAILURE);
		}
		pCELData->pMaskedX->pMatElement[ni] = nEntry[0];
		pCELData->pMaskedY->pMatElement[ni] = nEntry[1];
	}

	/* load outlier entries */
	for(ni=0; ni<(int)nOutlierCells; ni++)
	{
		if(little_endian_fread(nEntry, SHORT_SIZE, 2, fpCel, little_endian_machine) != 2)
		{
			printf("Error: Affy_LoadCELv4, cannot load outlier entries!\n");
			exit(EXIT_FAILURE);
		}
		pCELData->pOutlierX->pMatElement[ni] = nEntry[0];
		pCELData->pOutlierY->pMatElement[ni] = nEntry[1];
	}

	/* load sub-grids */
	pCELData->vSubGrids = NULL;
	if(nSubGrids > 0)
	{
		pCELData->vSubGrids = (struct tagCELSubGrid **)calloc(nSubGrids, sizeof(struct tagCELSubGrid *));
		if(pCELData->vSubGrids == NULL)
		{
			printf("Error: Affy_LoadCELv4, cannot create memory for loading subgrids!\n");
			exit(EXIT_FAILURE);
		}

		for(ni=0; ni<nSubGrids; ni++)
		{
			pCELData->vSubGrids[ni] = NULL;
			pCELData->vSubGrids[ni] = (struct tagCELSubGrid *)calloc(1, sizeof(struct tagCELSubGrid));
			if(pCELData->vSubGrids[ni] == NULL)
			{
				printf("Error: Affy_LoadCELv4, cannot create memory for loading subgrids!\n");
				exit(EXIT_FAILURE);
			}

			if(little_endian_fread(nPos, INT_SIZE, 2, fpCel, little_endian_machine) != 2)
			{
				printf("Error: Affy_LoadCELv4, cannot load row and column number!\n");
				exit(EXIT_FAILURE);
			}
			pCELData->vSubGrids[ni]->nRows = nPos[0];
			pCELData->vSubGrids[ni]->nCols = nPos[1];
			
			if(little_endian_fread(fI, FLOAT_SIZE, 2, fpCel, little_endian_machine) != 2)
			{
				printf("Error: Affy_LoadCELv4, cannot load upper left coordinates!\n");
				exit(EXIT_FAILURE);
			}
			pCELData->vSubGrids[ni]->fUpperLeft[0] = fI[0];
			pCELData->vSubGrids[ni]->fUpperLeft[1] = fI[1];
			
			if(little_endian_fread(fI, FLOAT_SIZE, 2, fpCel, little_endian_machine) != 2)
			{
				printf("Error: Affy_LoadCELv4, cannot load upper right coordinates!\n");
				exit(EXIT_FAILURE);
			}
			pCELData->vSubGrids[ni]->fUpperRight[0] = fI[0];
			pCELData->vSubGrids[ni]->fUpperRight[1] = fI[1];
						
			if(little_endian_fread(fI, FLOAT_SIZE, 2, fpCel, little_endian_machine) != 2)
			{
				printf("Error: Affy_LoadCELv4, cannot load lower left coordinates!\n");
				exit(EXIT_FAILURE);
			}
			pCELData->vSubGrids[ni]->fLowerLeft[0] = fI[0];
			pCELData->vSubGrids[ni]->fLowerLeft[1] = fI[1];
						
			if(little_endian_fread(fI, FLOAT_SIZE, 2, fpCel, little_endian_machine) != 2)
			{
				printf("Error: Affy_LoadCELv4, cannot load lower right coordinates!\n");
				exit(EXIT_FAILURE);
			}
			pCELData->vSubGrids[ni]->fLowerRight[0] = fI[0];
			pCELData->vSubGrids[ni]->fLowerRight[1] = fI[1];
			
			if(little_endian_fread(nPos, INT_SIZE, 2, fpCel, little_endian_machine) != 2)
			{
				printf("Error: Affy_LoadCELv4, cannot load left and top position!\n");
				exit(EXIT_FAILURE);
			}
			pCELData->vSubGrids[ni]->nCellLeft = nPos[0];
			pCELData->vSubGrids[ni]->nCellTop = nPos[1];
						
			if(little_endian_fread(nPos, INT_SIZE, 2, fpCel, little_endian_machine) != 2)
			{
				printf("Error: Affy_LoadCELv4, cannot load right and bottom position!\n");
				exit(EXIT_FAILURE);
			}
			pCELData->vSubGrids[ni]->nCellRight = nPos[0];
			pCELData->vSubGrids[ni]->nCellBottom = nPos[1];
		}
	}

	/* close file */
	fclose(fpCel);

	/* return */
	return pCELData;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_LoadCELv4_Fast()                                                  */
/*  Loading raw data from a single affymetrix's *.CEL (v4) file            */
/* ----------------------------------------------------------------------- */ 
struct tagCELData *Affy_LoadCELv4_Fast(char strFileName[])
{
	/* define */
	struct tagCELData *pCELData = NULL;

	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));
	size_t unitsize;
	char *tempmem;
	char *tempptr;
	
	FILE *fpCel;
	int nMagicnumber;
	int nVersionnumber;
	int nCols, nRows;
	int nNumberCells;

	int nHeaderLength;
	int nAlgorithmLength;
	int nAlgorithmParamLength;

	int nCellMargin;
	unsigned int nOutlierCells;
	unsigned int nMaskedCells;
	int nSubGrids;

	float fI[2];
	short int nCount;
	short int nEntry[2];
	int nPos[2];
	int ni;

	int nTotalX,nTotalY;
	
	/* init */
	nTotalX = 0;
	nTotalY = 0;

	/* load *.CEL */
	fpCel = NULL;
	fpCel = fopen(strFileName, "rb");
	if(fpCel == NULL)
	{
		printf("Error: Affy_LoadCELv4, cannot open *.CEL file!\n");
		exit(EXIT_FAILURE);
	}

	/* load magic number */
	if(little_endian_fread(&nMagicnumber, INT_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_LoadCELv4, cannot load magic number.\n");
		return NULL;
	}
	
	if(nMagicnumber != 64)
	{
		printf("Error: Affy_LoadCELv4, cannot load magic number correctly.\n");
		return NULL;
	}

	/* load version number */
	if(little_endian_fread(&nVersionnumber, INT_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_LoadCELv4, cannot load version number.\n");
		return NULL;
	}
	
	if(nVersionnumber != 4)
	{
		printf("Error: Affy_LoadCELv4, cannot load version number correctly.\n");
		return NULL;
	}

	/* load number of columns, rows and cells */
	if(little_endian_fread(&nCols, INT_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_LoadCELv4, cannot load number of columns.\n");
		return NULL;
	}
	if(little_endian_fread(&nRows, INT_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_LoadCELv4, cannot load number of rows.\n");
		return NULL;
	}
	if(little_endian_fread(&nNumberCells, INT_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_LoadCELv4, cannot load number of cells.\n");
		return NULL;
	}
	if(nNumberCells != nCols*nRows)
	{
		printf("Error: Affy_LoadCELv4, cannot load number of columns, rows and cells correctly.\n");
		return NULL;
	}

	/* create CELData object */
	pCELData = Affy_CELData_Create();
	if(pCELData == NULL)
	{
		printf("Error: Affy_LoadCELv4, cannot create CELData structure.\n");
		return NULL;
	}

	pCELData->nMagicnumber = nMagicnumber;
	pCELData->nVersionnumber = nVersionnumber;
	pCELData->nCols = nCols;
	pCELData->nRows = nRows;
	pCELData->nNumberCells = nNumberCells;
	pCELData->nTotalX = nCols;
	pCELData->nTotalY = nRows;

	/* load header info */
	if(little_endian_fread(&nHeaderLength, INT_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_LoadCELv4, cannot load header length.\n");
		exit(EXIT_FAILURE);
	}
	pCELData->nHeaderLength = nHeaderLength;
	pCELData->vHeader = NULL;
	pCELData->vHeader = CreateString(nHeaderLength);
	if( (int)little_endian_fread(pCELData->vHeader->m_pString, sizeof(char), nHeaderLength, fpCel,  little_endian_machine) != nHeaderLength)
	{
		printf("Error: Affy_LoadCELv4, cannot load header.\n");
		exit(EXIT_FAILURE);
	}
	pCELData->vHeader->m_pString[nHeaderLength] = '\0';
	
	/* load algorithm */
	if(little_endian_fread(&nAlgorithmLength, INT_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_LoadCELv4, cannot load algorithm length.\n");
		exit(EXIT_FAILURE);
	}
	pCELData->nAlgorithmLength = nAlgorithmLength;
	pCELData->vAlgorithm = NULL;
	pCELData->vAlgorithm = CreateString(nAlgorithmLength);
	if( (int)little_endian_fread(pCELData->vAlgorithm->m_pString, sizeof(char), nAlgorithmLength, fpCel, little_endian_machine) != nAlgorithmLength)
	{
		printf("Error: Affy_LoadCELv4, cannot load algorithm.\n");
		exit(EXIT_FAILURE);
	}
	pCELData->vAlgorithm->m_pString[nAlgorithmLength] = '\0';

	/* load algorithm parameters */
	if(little_endian_fread(&nAlgorithmParamLength, INT_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_LoadCELv4, cannot load algorithm parameter length.\n");
		exit(EXIT_FAILURE);
	}
	pCELData->nAlgorithmParamLength = nAlgorithmParamLength;
	pCELData->vAlgorithmParam = NULL;
	pCELData->vAlgorithmParam = CreateString(nAlgorithmParamLength);
	if( (int)little_endian_fread(pCELData->vAlgorithmParam->m_pString, sizeof(char), nAlgorithmParamLength, fpCel, little_endian_machine) != nAlgorithmParamLength)
	{
		printf("Error: Affy_LoadCELv4, cannot load algorithm parameters.\n");
		exit(EXIT_FAILURE);
	}
	pCELData->vAlgorithmParam->m_pString[nAlgorithmParamLength] = '\0';

	/* load cell margin */
	if(little_endian_fread(&nCellMargin, INT_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_LoadCELv4, cannot load cell margin.\n");
		exit(EXIT_FAILURE);
	}
	pCELData->nCellMargin = nCellMargin;

	/* load number of outlier cells */
	if(little_endian_fread(&nOutlierCells, DWORD_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_LoadCELv4, cannot load number of outlier cells.\n");
		exit(EXIT_FAILURE);
	}
	pCELData->nOutlierCells = nOutlierCells;

	/* load number of masked cells */
	if(little_endian_fread(&nMaskedCells, DWORD_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_LoadCELv4, cannot load number of masked cells.\n");
		exit(EXIT_FAILURE);
	}
	pCELData->nMaskedCells = nMaskedCells;

	/* load number of sub-grids */
	if(little_endian_fread(&nSubGrids, INT_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_LoadCELv4, cannot load number of sub-grids.\n");
		exit(EXIT_FAILURE);
	}
	pCELData->nSubGrids = nSubGrids;

	/* allocate memory */
	pCELData->pIntensity = NULL;
	pCELData->pIntensity = CreateDoubleMatrix(1, nNumberCells);
	if((pCELData->pIntensity == NULL) && (nNumberCells>0))
	{
		printf("Error: Affy_LoadCELv4, cannot allocate enough memory for loading intensities!\n");
		exit(EXIT_FAILURE);
	}

	pCELData->pSD = NULL;
	pCELData->pSD = CreateDoubleMatrix(1, nNumberCells);
	if((pCELData->pSD == NULL) && (nNumberCells>0))
	{
		printf("Error: Affy_LoadCELv4, cannot allocate enough memory for loading SD!\n");
		exit(EXIT_FAILURE);
	}

	pCELData->pPixelNum = NULL;
	pCELData->pPixelNum = CreateIntMatrix(1, nNumberCells);
	if((pCELData->pPixelNum == NULL) && (nNumberCells>0))
	{
		printf("Error: Affy_LoadCELv4, cannot allocate enough memory for loading pixel numbers!\n");
		exit(EXIT_FAILURE);
	}

	pCELData->pMaskedX = NULL;
	pCELData->pMaskedX = CreateIntMatrix(1, nMaskedCells);
	if((pCELData->pMaskedX == NULL) && (nMaskedCells>0))
	{
		printf("Error: Affy_LoadCELv4, cannot allocate enough memory for loading masked cells!\n");
		exit(EXIT_FAILURE);
	}

	pCELData->pMaskedY = NULL;
	pCELData->pMaskedY = CreateIntMatrix(1, nMaskedCells);
	if((pCELData->pMaskedY == NULL) && (nMaskedCells>0))
	{
		printf("Error: Affy_LoadCELv4, cannot allocate enough memory for loading masked cells!\n");
		exit(EXIT_FAILURE);
	}

	pCELData->pOutlierX = NULL;
	pCELData->pOutlierX = CreateIntMatrix(1, nOutlierCells);
	if((pCELData->pOutlierX == NULL) && (nOutlierCells>0))
	{
		printf("Error: Affy_LoadCELv4, cannot allocate enough memory for loading outlier cells!\n");
		exit(EXIT_FAILURE);
	}

	pCELData->pOutlierY = NULL;
	pCELData->pOutlierY = CreateIntMatrix(1, nOutlierCells);
	if((pCELData->pOutlierY == NULL) && (nOutlierCells>0))
	{
		printf("Error: Affy_LoadCELv4, cannot allocate enough memory for loading outlier cells!\n");
		exit(EXIT_FAILURE);
	}

	pCELData->nModifiedCells = 0;
	pCELData->pModifiedX = NULL;
	pCELData->pModifiedY = NULL;
	pCELData->pModifiedOrig = NULL;

	/* load intensity */
	unitsize = FLOAT_SIZE+FLOAT_SIZE+SHORT_SIZE;
	tempmem = NULL;
	tempmem = (char *)calloc(nNumberCells, unitsize);
	if(tempmem == NULL)
	{
		printf("Error: Affy_LoadCELv4, cannot allocate memory for loading intensity data.\n");
		exit(EXIT_FAILURE);
	}
	if( (int)fread(tempmem, unitsize, nNumberCells, fpCel) != nNumberCells)
	{
		printf("Error: Affy_LoadCELv4, cannot load intensity data correctly.\n");
		exit(EXIT_FAILURE);
	}
	tempptr = tempmem;
	for(ni=0; ni<nNumberCells; ni++)
	{
		if(little_endian_machine == 0)
		{
			reverse_buf(tempptr, FLOAT_SIZE);
		}
        memcpy(fI, tempptr, FLOAT_SIZE);
		tempptr += FLOAT_SIZE;

		if(little_endian_machine == 0)
		{
			reverse_buf(tempptr, FLOAT_SIZE);
		}
        memcpy(fI+1, tempptr, FLOAT_SIZE);
		tempptr += FLOAT_SIZE;

		if(little_endian_machine == 0)
		{
			reverse_buf(tempptr, SHORT_SIZE);
		}
        memcpy(&nCount, tempptr, SHORT_SIZE);
		tempptr += SHORT_SIZE;

		pCELData->pIntensity->pMatElement[ni] = fI[0];
		pCELData->pSD->pMatElement[ni] = fI[1];
		pCELData->pPixelNum->pMatElement[ni] = nCount;
	}
	free(tempmem);

	/* load masked entries */
	if(nMaskedCells > 0)
	{
		unitsize = SHORT_SIZE+SHORT_SIZE;
		tempmem = NULL;
		tempmem = (char *)calloc(nMaskedCells, unitsize);
		if(tempmem == NULL)
		{
			printf("Error: Affy_LoadCELv4, cannot allocate memory for loading masked cells.\n");
			exit(EXIT_FAILURE);
		}
		if( (int)fread(tempmem, unitsize, nMaskedCells, fpCel) != nMaskedCells)
		{
			printf("Error: Affy_LoadCELv4, cannot load masked cells correctly.\n");
			exit(EXIT_FAILURE);
		}
		tempptr = tempmem;
		for(ni=0; ni<(int)nMaskedCells; ni++)
		{
			if(little_endian_machine == 0)
			{
				reverse_buf(tempptr, SHORT_SIZE);
			}
			memcpy(nEntry, tempptr, SHORT_SIZE);
			tempptr += SHORT_SIZE;

			if(little_endian_machine == 0)
			{
				reverse_buf(tempptr, SHORT_SIZE);
			}
			memcpy(nEntry+1, tempptr, SHORT_SIZE);
			tempptr += SHORT_SIZE;

			pCELData->pMaskedX->pMatElement[ni] = nEntry[0];
			pCELData->pMaskedY->pMatElement[ni] = nEntry[1];
		}
		free(tempmem);
	}

	/* load outlier entries */
	if(nOutlierCells > 0)
	{
		unitsize = SHORT_SIZE+SHORT_SIZE;
		tempmem = NULL;
		tempmem = (char *)calloc(nOutlierCells, unitsize);
		if(tempmem == NULL)
		{
			printf("Error: Affy_LoadCELv4, cannot allocate memory for loading outlier cells.\n");
			exit(EXIT_FAILURE);
		}
		if( (int)fread(tempmem, unitsize, nOutlierCells, fpCel) != nOutlierCells)
		{
			printf("Error: Affy_LoadCELv4, cannot load outlier cells correctly.\n");
			exit(EXIT_FAILURE);
		}
		tempptr = tempmem;
		for(ni=0; ni<(int)nOutlierCells; ni++)
		{
			if(little_endian_machine == 0)
			{
				reverse_buf(tempptr, SHORT_SIZE);
			}
			memcpy(nEntry, tempptr, SHORT_SIZE);
			tempptr += SHORT_SIZE;

			if(little_endian_machine == 0)
			{
				reverse_buf(tempptr, SHORT_SIZE);
			}
			memcpy(nEntry+1, tempptr, SHORT_SIZE);
			tempptr += SHORT_SIZE;

			pCELData->pOutlierX->pMatElement[ni] = nEntry[0];
			pCELData->pOutlierY->pMatElement[ni] = nEntry[1];
		}
		free(tempmem);
	}

	/* load sub-grids */
	pCELData->vSubGrids = NULL;
	if(nSubGrids > 0)
	{
		pCELData->vSubGrids = (struct tagCELSubGrid **)calloc(nSubGrids, sizeof(struct tagCELSubGrid *));
		if(pCELData->vSubGrids == NULL)
		{
			printf("Error: Affy_LoadCELv4, cannot create memory for loading subgrids!\n");
			exit(EXIT_FAILURE);
		}

		for(ni=0; ni<nSubGrids; ni++)
		{
			pCELData->vSubGrids[ni] = NULL;
			pCELData->vSubGrids[ni] = (struct tagCELSubGrid *)calloc(1, sizeof(struct tagCELSubGrid));
			if(pCELData->vSubGrids[ni] == NULL)
			{
				printf("Error: Affy_LoadCELv4, cannot create memory for loading subgrids!\n");
				exit(EXIT_FAILURE);
			}

			if(little_endian_fread(nPos, INT_SIZE, 2, fpCel, little_endian_machine) != 2)
			{
				printf("Error: Affy_LoadCELv4, cannot load row and column number!\n");
				exit(EXIT_FAILURE);
			}
			pCELData->vSubGrids[ni]->nRows = nPos[0];
			pCELData->vSubGrids[ni]->nCols = nPos[1];
			
			if(little_endian_fread(fI, FLOAT_SIZE, 2, fpCel, little_endian_machine) != 2)
			{
				printf("Error: Affy_LoadCELv4, cannot load upper left coordinates!\n");
				exit(EXIT_FAILURE);
			}
			pCELData->vSubGrids[ni]->fUpperLeft[0] = fI[0];
			pCELData->vSubGrids[ni]->fUpperLeft[1] = fI[1];
			
			if(little_endian_fread(fI, FLOAT_SIZE, 2, fpCel, little_endian_machine) != 2)
			{
				printf("Error: Affy_LoadCELv4, cannot load upper right coordinates!\n");
				exit(EXIT_FAILURE);
			}
			pCELData->vSubGrids[ni]->fUpperRight[0] = fI[0];
			pCELData->vSubGrids[ni]->fUpperRight[1] = fI[1];
						
			if(little_endian_fread(fI, FLOAT_SIZE, 2, fpCel, little_endian_machine) != 2)
			{
				printf("Error: Affy_LoadCELv4, cannot load lower left coordinates!\n");
				exit(EXIT_FAILURE);
			}
			pCELData->vSubGrids[ni]->fLowerLeft[0] = fI[0];
			pCELData->vSubGrids[ni]->fLowerLeft[1] = fI[1];
						
			if(little_endian_fread(fI, FLOAT_SIZE, 2, fpCel, little_endian_machine) != 2)
			{
				printf("Error: Affy_LoadCELv4, cannot load lower right coordinates!\n");
				exit(EXIT_FAILURE);
			}
			pCELData->vSubGrids[ni]->fLowerRight[0] = fI[0];
			pCELData->vSubGrids[ni]->fLowerRight[1] = fI[1];
			
			if(little_endian_fread(nPos, INT_SIZE, 2, fpCel, little_endian_machine) != 2)
			{
				printf("Error: Affy_LoadCELv4, cannot load left and top position!\n");
				exit(EXIT_FAILURE);
			}
			pCELData->vSubGrids[ni]->nCellLeft = nPos[0];
			pCELData->vSubGrids[ni]->nCellTop = nPos[1];
						
			if(little_endian_fread(nPos, INT_SIZE, 2, fpCel, little_endian_machine) != 2)
			{
				printf("Error: Affy_LoadCELv4, cannot load right and bottom position!\n");
				exit(EXIT_FAILURE);
			}
			pCELData->vSubGrids[ni]->nCellRight = nPos[0];
			pCELData->vSubGrids[ni]->nCellBottom = nPos[1];
		}
	}

	/* close file */
	fclose(fpCel);

	/* return */
	return pCELData;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_SaveCELv4()                                                       */
/*  Saving raw data to a affymetrix's *.CEL (v4) file                      */
/* ----------------------------------------------------------------------- */ 
int Affy_SaveCELv4(char strFileName[], struct tagCELData *pCELData)
{
	/* define */
	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));
	
	FILE *fpCel;

	int nMagicnumber;
	int nVersionnumber;
	int nCols, nRows;
	int nNumberCells;

	int nHeaderLength;
	int nAlgorithmLength;
	int nAlgorithmParamLength;

	unsigned int nOutlierCells;
	unsigned int nMaskedCells;
	int nSubGrids;

	float fI[2];
	short int nCount;
	short int nEntry[2];
	int nPos[2];
	int ni;

	/* initial check */
	if(pCELData == NULL)
	{
		printf("Warning: empty CEL data!\n");
		return PROC_SUCCESS;
	}
	
	/* load *.CEL */
	fpCel = NULL;
	fpCel = fopen(strFileName, "wb");
	if(fpCel == NULL)
	{
		printf("Error: Affy_SaveCELv4, cannot open *.CEL file!\n");
		return PROC_FAILURE;
	}

	/* write magic number */
	nMagicnumber = 64;
	if(little_endian_fwrite(&nMagicnumber, INT_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_SaveCELv4, cannot write magic number.\n");
		return PROC_FAILURE;
	}
	
	/* write version number */
	nVersionnumber = 4;
	if(little_endian_fwrite(&nVersionnumber, INT_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_SaveCELv4, cannot write version number.\n");
		return PROC_FAILURE;
	}
	
	/* write number of columns, rows and cells */
	nCols = pCELData->nCols;
	nRows = pCELData->nRows;
	nNumberCells = pCELData->nNumberCells;
	if(nNumberCells != nCols*nRows)
	{
		printf("Error: Affy_SaveCELv4, cannot write number of columns, rows and cells correctly.\n");
		return PROC_FAILURE;
	}

	if(little_endian_fwrite(&nCols, INT_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_SaveCELv4, cannot write number of columns.\n");
		return PROC_FAILURE;
	}
	if(little_endian_fwrite(&nRows, INT_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_SaveCELv4, cannot write number of rows.\n");
		return PROC_FAILURE;
	}
	if(little_endian_fwrite(&nNumberCells, INT_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_SaveCELv4, cannot write number of cells.\n");
		return PROC_FAILURE;
	}

	/* write header info */
	nHeaderLength = pCELData->nHeaderLength;
	if(little_endian_fwrite(&nHeaderLength, INT_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_SaveCELv4, cannot write header length.\n");
		return PROC_FAILURE;
	}
	if( (int)little_endian_fwrite(pCELData->vHeader->m_pString, sizeof(char), nHeaderLength, fpCel,  little_endian_machine) != nHeaderLength)
	{
		printf("Error: Affy_SaveCELv4, cannot write header.\n");
		return PROC_FAILURE;
	}
	
	/* write algorithm */
	nAlgorithmLength = pCELData->nAlgorithmLength;
	if(little_endian_fwrite(&nAlgorithmLength, INT_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_SaveCELv4, cannot write algorithm length.\n");
		return PROC_FAILURE;
	}
	if( (int)little_endian_fwrite(pCELData->vAlgorithm->m_pString, sizeof(char), nAlgorithmLength, fpCel, little_endian_machine) != nAlgorithmLength)
	{
		printf("Error: Affy_SaveCELv4, cannot write algorithm.\n");
		return PROC_FAILURE;
	}
	
	/* write algorithm parameters */
	nAlgorithmParamLength = pCELData->nAlgorithmParamLength;
	if(little_endian_fwrite(&nAlgorithmParamLength, INT_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_SaveCELv4, cannot write algorithm parameter length.\n");
		return PROC_FAILURE;
	}
	if( (int)little_endian_fwrite(pCELData->vAlgorithmParam->m_pString, sizeof(char), nAlgorithmParamLength, fpCel, little_endian_machine) != nAlgorithmParamLength)
	{
		printf("Error: Affy_SaveCELv4, cannot write algorithm parameters.\n");
		return PROC_FAILURE;
	}
	
	/* write cell margin */
	if(little_endian_fwrite(&(pCELData->nCellMargin), INT_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_SaveCELv4, cannot write cell margin.\n");
		return PROC_FAILURE;
	}

	/* write number of outlier cells */
	nOutlierCells = pCELData->nOutlierCells;
	if(little_endian_fwrite(&nOutlierCells, DWORD_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_SaveCELv4, cannot write number of outlier cells.\n");
		return PROC_FAILURE;
	}
	
	/* write number of masked cells */
	nMaskedCells = pCELData->nMaskedCells;
	if(little_endian_fwrite(&nMaskedCells, DWORD_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_SaveCELv4, cannot write number of masked cells.\n");
		return PROC_FAILURE;
	}

	/* write number of sub-grids */
	nSubGrids = pCELData->nSubGrids;
	if(little_endian_fwrite(&nSubGrids, INT_SIZE, 1, fpCel, little_endian_machine) != 1)
	{
		printf("Error: Affy_SaveCELv4, cannot write number of sub-grids.\n");
		return PROC_FAILURE;
	}
	
	/* write intensity */
	for(ni=0; ni<nNumberCells; ni++)
	{
		fI[0] = (float)(pCELData->pIntensity->pMatElement[ni]);
		fI[1] = (float)(pCELData->pSD->pMatElement[ni]);
		nCount = (short int)(pCELData->pPixelNum->pMatElement[ni]);
		if(little_endian_fwrite(fI, FLOAT_SIZE, 2, fpCel, little_endian_machine) != 2)
		{
			printf("Error: Affy_SaveCELv4, cannot write intensity values!\n");
			exit(EXIT_FAILURE);
		}
		if(little_endian_fwrite(&nCount, SHORT_SIZE, 1, fpCel, little_endian_machine) != 1)
		{
			printf("Error: Affy_SaveCELv4, cannot write pixel count!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* write masked entries */
	for(ni=0; ni<(int)nMaskedCells; ni++)
	{
		nEntry[0] = (short int)(pCELData->pMaskedX->pMatElement[ni]);
		nEntry[1] = (short int)(pCELData->pMaskedY->pMatElement[ni]);
		if(little_endian_fwrite(nEntry, SHORT_SIZE, 2, fpCel, little_endian_machine) != 2)
		{
			printf("Error: Affy_SaveCELv4, cannot write masked entries!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* write outlier entries */
	for(ni=0; ni<(int)nOutlierCells; ni++)
	{
		nEntry[0] = (short int)(pCELData->pOutlierX->pMatElement[ni]);
		nEntry[1] = (short int)(pCELData->pOutlierY->pMatElement[ni]);
		if(little_endian_fwrite(nEntry, SHORT_SIZE, 2, fpCel, little_endian_machine) != 2)
		{
			printf("Error: Affy_SaveCELv4, cannot write outlier entries!\n");
			exit(EXIT_FAILURE);
		}
		
	}

	/* write sub-grids */
	for(ni=0; ni<nSubGrids; ni++)
	{
		nPos[0] = pCELData->vSubGrids[ni]->nRows;
		nPos[1] = pCELData->vSubGrids[ni]->nCols;

		if(little_endian_fwrite(nPos, INT_SIZE, 2, fpCel, little_endian_machine) != 2)
		{
			printf("Error: Affy_SaveCELv4, cannot write row and column number!\n");
			exit(EXIT_FAILURE);
		}
		
		fI[0] = pCELData->vSubGrids[ni]->fUpperLeft[0];
		fI[1] = pCELData->vSubGrids[ni]->fUpperLeft[1];
		if(little_endian_fwrite(fI, FLOAT_SIZE, 2, fpCel, little_endian_machine) != 2)
		{
			printf("Error: Affy_SaveCELv4, cannot write upper left coordinates!\n");
			exit(EXIT_FAILURE);
		}
		
		fI[0] = pCELData->vSubGrids[ni]->fUpperRight[0];
		fI[1] = pCELData->vSubGrids[ni]->fUpperRight[1];
		if(little_endian_fwrite(fI, FLOAT_SIZE, 2, fpCel, little_endian_machine) != 2)
		{
			printf("Error: Affy_SaveCELv4, cannot write upper right coordinates!\n");
			exit(EXIT_FAILURE);
		}
				
		fI[0] = pCELData->vSubGrids[ni]->fLowerLeft[0];
		fI[1] = pCELData->vSubGrids[ni]->fLowerLeft[1];
		if(little_endian_fwrite(fI, FLOAT_SIZE, 2, fpCel, little_endian_machine) != 2)
		{
			printf("Error: Affy_SaveCELv4, cannot write lower left coordinates!\n");
			exit(EXIT_FAILURE);
		}
		
		fI[0] = pCELData->vSubGrids[ni]->fLowerRight[0];
		fI[1] = pCELData->vSubGrids[ni]->fLowerRight[1];
		if(little_endian_fwrite(fI, FLOAT_SIZE, 2, fpCel, little_endian_machine) != 2)
		{
			printf("Error: Affy_SaveCELv4, cannot write lower right coordinates!\n");
			exit(EXIT_FAILURE);
		}
		
		nPos[0] = pCELData->vSubGrids[ni]->nCellLeft;
		nPos[1] = pCELData->vSubGrids[ni]->nCellTop;
		if(little_endian_fwrite(nPos, INT_SIZE, 2, fpCel, little_endian_machine) != 2)
		{
			printf("Error: Affy_SaveCELv4, cannot write left and top position!\n");
			exit(EXIT_FAILURE);
		}
		
		nPos[0] = pCELData->vSubGrids[ni]->nCellRight;
		nPos[1] = pCELData->vSubGrids[ni]->nCellBottom;
		if(little_endian_fwrite(nPos, INT_SIZE, 2, fpCel, little_endian_machine) != 2)
		{
			printf("Error: Affy_SaveCELv4, cannot write right and bottom position!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* close file */
	fclose(fpCel);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_LoadCEL_CmdCslv1()                                                */
/*  Loading raw data from a single affymetrix's command console *.CEL file */
/* ----------------------------------------------------------------------- */ 
struct tagCELData *Affy_LoadCEL_CmdCslv1(char strFileName[])
{
	/* define */
	struct tagString *pParamString = NULL;
	struct tagString *pHeaderString = NULL;
	char strNewLine[LINE_LENGTH];
	char *chp;

	struct tagCELData *pCELData = NULL;
	struct tagGenericFileObject *pGenericData = NULL;

	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));
	char *tempptr1,*tempptr2,*tempptr3;
	
	int nCols, nRows;
	int nNumberCells = 0;

	int nHeaderLength = 0;
	char strAlgorithm[MED_LINE_LENGTH];
	int nAlgorithmLength = 0;
	int nAlgorithmParamLength = 0;

	int nCellMargin = 0;
	unsigned int nOutlierCells = 0;
	unsigned int nMaskedCells = 0;
	int nSubGrids = 0;

	float fI[2];
	short int nCount;
	short int nEntry[2];
	int ni;

	int nTotalX,nTotalY;
	
	/* init */
	nTotalX = 0;
	nTotalY = 0;

	/* load generic data */
	/* printf("char=%d\n", sizeof(char));
	printf("wchar=%d\n", sizeof(wchar_t));
	printf("unsigned char=%d\n", sizeof(unsigned char));
	printf("short=%d\n", sizeof(short));
	printf("unsigned short=%d\n", sizeof(unsigned short));
	printf("int=%d\n", sizeof(int));
	printf("unsigned int=%d\n", sizeof(unsigned int));
	printf("long=%d\n", sizeof(long));
	printf("unsigned long=%d\n", sizeof(unsigned long));
	printf("float=%d\n", sizeof(float));
	printf("double=%d\n", sizeof(double)); */
	/* printf("Load Generic Data!\n"); */ 

	pGenericData = Affy_GenericFileObject_Load(strFileName);
	
	/* printf("Generic Data Loaded!\n"); */
	
	if(pGenericData == NULL)
	{
		printf("Error: Affy_LoadCEL_CmdCslv1, cannot load %s! The file is either empty or does not exist.\n", strFileName);
		return NULL;
	}

	if( (pGenericData->pFileHeader == NULL) || (pGenericData->pDataHeader == NULL) || (pGenericData->vDataGroup == NULL))
	{
		printf("Error: Affy_LoadCEL_CmdCslv1, cannot load file correctly.\n");
		Affy_GenericFileObject_Delete(&pGenericData);
		return NULL;
	}
	else
	{
		if( (pGenericData->pFileHeader->nDataGroupNum != 1) ||
			(strcmp(pGenericData->pDataHeader->pDataTypeID->m_pString, "affymetrix-calvin-intensity") != 0) )
		{
			printf("Error: Affy_LoadCEL_CmdCslv1, cannot load file correctly.\n");
			Affy_GenericFileObject_Delete(&pGenericData);
			return NULL;
		}
		else if(pGenericData->vDataGroup[0]->nDataSetNum < 5)
		{
			printf("Error: Affy_LoadCEL_CmdCslv1, cannot load file correctly.\n");
			Affy_GenericFileObject_Delete(&pGenericData);
			return NULL;
		}
	}

	for(ni=0; ni<5; ni++)
	{
		if(pGenericData->vDataGroup[0]->vDataSets[ni] == NULL)
		{
			printf("Error: Affy_LoadCEL_CmdCslv1, cannot load file correctly.\n");
			Affy_GenericFileObject_Delete(&pGenericData);
			return NULL;
		}
	}

	nNumberCells = pGenericData->vDataGroup[0]->vDataSets[0]->nRowNum;
	if( (nNumberCells != pGenericData->vDataGroup[0]->vDataSets[1]->nRowNum) ||
		(nNumberCells != pGenericData->vDataGroup[0]->vDataSets[2]->nRowNum) )
	{
		printf("Error: Affy_LoadCEL_CmdCslv1, cannot load file correctly.\n");
		Affy_GenericFileObject_Delete(&pGenericData);
		return NULL;
	}
	nOutlierCells = pGenericData->vDataGroup[0]->vDataSets[3]->nRowNum;
	nMaskedCells = pGenericData->vDataGroup[0]->vDataSets[4]->nRowNum;


	/* printf("Load Parameter!\n"); */
	strcpy(strNewLine, "\n");
	for(ni=0; ni<pGenericData->pDataHeader->nParamNum; ni++)
	{
		pParamString = NULL;
		pParamString = Affy_ParamValueType2String(pGenericData->pDataHeader->vParamName[ni]->m_pWString, 
			pGenericData->pDataHeader->vParamValue[ni]->m_pString, pGenericData->pDataHeader->vParamValue[ni]->m_nLength,
			pGenericData->pDataHeader->vParamType[ni]->m_pWString);
		StringAddTail(&pHeaderString, pParamString->m_pString);
		StringAddTail(&pHeaderString, strNewLine);

		/* printf("%s\n", pParamString->m_pString); */
		chp = strchr(pParamString->m_pString, '=');
		if(chp != NULL)
		{
			chp++;
			if( wcscmp(pGenericData->pDataHeader->vParamName[ni]->m_pWString, L"affymetrix-cel-cols") == 0)
			{
				nCols = atoi(chp);
			}
			else if( wcscmp(pGenericData->pDataHeader->vParamName[ni]->m_pWString, L"affymetrix-cel-rows") == 0)
			{
				nRows = atoi(chp);
			}
			else if( wcscmp(pGenericData->pDataHeader->vParamName[ni]->m_pWString, L"affymetrix-algorithm-name") == 0 )
			{
				strcpy(strAlgorithm, chp);
			}
			else if( wcscmp(pGenericData->pDataHeader->vParamName[ni]->m_pWString, L"affymetrix-algorithm-param-CellMargin") == 0 )
			{
				nCellMargin = atoi(chp);
			}
		}

		DeleteString(pParamString);
		pParamString = NULL;
	}
	StringAddTail(&pHeaderString, strNewLine);

	if(nNumberCells != nCols*nRows)
	{
		printf("Error: Affy_LoadCEL_CmdCslv1, cannot write number of columns, rows and cells correctly.\n");
		return PROC_FAILURE;
	}

	/* create CELData object */
	/* printf("Create CEL!\n"); */
	pCELData = Affy_CELData_Create();
	if(pCELData == NULL)
	{
		printf("Error: Affy_LoadCEL_CmdCslv1, cannot create CELData structure.\n");
		return NULL;
	}

	pCELData->nMagicnumber = 64;
	pCELData->nVersionnumber = 4;
	pCELData->nCols = nCols;
	pCELData->nRows = nRows;
	pCELData->nNumberCells = nNumberCells;
	pCELData->nTotalX = nCols;
	pCELData->nTotalY = nRows;

	/* load header info */
	nHeaderLength = pHeaderString->m_nLength;
	pCELData->nHeaderLength = nHeaderLength;
	pCELData->vHeader = NULL;
	pCELData->vHeader = pHeaderString;
	
	/* load algorithm */
	pCELData->nAlgorithmLength = strlen(strAlgorithm);
	pCELData->vAlgorithm = NULL;
	StringAddTail(&(pCELData->vAlgorithm), strAlgorithm);

	/* load algorithm parameters */
	pCELData->nAlgorithmParamLength = 0;
	pCELData->vAlgorithmParam = NULL;
	pCELData->vAlgorithmParam = CreateString(nAlgorithmParamLength);
	pCELData->vAlgorithmParam->m_pString[nAlgorithmParamLength] = '\0';

	/* load cell margin */
	pCELData->nCellMargin = nCellMargin;

	/* load number of outlier cells */
	pCELData->nOutlierCells = nOutlierCells;

	/* load number of masked cells */
	pCELData->nMaskedCells = nMaskedCells;

	/* load number of sub-grids */
	nSubGrids = 0;
	pCELData->nSubGrids = nSubGrids;

	/* allocate memory */
	/* printf("Load Intensity!\n"); */
	pCELData->pIntensity = NULL;
	pCELData->pIntensity = CreateDoubleMatrix(1, nNumberCells);
	if((pCELData->pIntensity == NULL) && (nNumberCells>0))
	{
		printf("Error: Affy_LoadCEL_CmdCslv1, cannot allocate enough memory for loading intensities!\n");
		exit(EXIT_FAILURE);
	}

	pCELData->pSD = NULL;
	pCELData->pSD = CreateDoubleMatrix(1, nNumberCells);
	if((pCELData->pSD == NULL) && (nNumberCells>0))
	{
		printf("Error: Affy_LoadCEL_CmdCslv1, cannot allocate enough memory for loading SD!\n");
		exit(EXIT_FAILURE);
	}

	pCELData->pPixelNum = NULL;
	pCELData->pPixelNum = CreateIntMatrix(1, nNumberCells);
	if((pCELData->pPixelNum == NULL) && (nNumberCells>0))
	{
		printf("Error: Affy_LoadCEL_CmdCslv1, cannot allocate enough memory for loading pixel numbers!\n");
		exit(EXIT_FAILURE);
	}

	pCELData->pMaskedX = NULL;
	pCELData->pMaskedX = CreateIntMatrix(1, nMaskedCells);
	if((pCELData->pMaskedX == NULL) && (nMaskedCells>0))
	{
		printf("Error: Affy_LoadCEL_CmdCslv1, cannot allocate enough memory for loading masked cells!\n");
		exit(EXIT_FAILURE);
	}

	pCELData->pMaskedY = NULL;
	pCELData->pMaskedY = CreateIntMatrix(1, nMaskedCells);
	if((pCELData->pMaskedY == NULL) && (nMaskedCells>0))
	{
		printf("Error: Affy_LoadCEL_CmdCslv1, cannot allocate enough memory for loading masked cells!\n");
		exit(EXIT_FAILURE);
	}

	pCELData->pOutlierX = NULL;
	pCELData->pOutlierX = CreateIntMatrix(1, nOutlierCells);
	if((pCELData->pOutlierX == NULL) && (nOutlierCells>0))
	{
		printf("Error: Affy_LoadCEL_CmdCslv1, cannot allocate enough memory for loading outlier cells!\n");
		exit(EXIT_FAILURE);
	}

	pCELData->pOutlierY = NULL;
	pCELData->pOutlierY = CreateIntMatrix(1, nOutlierCells);
	if((pCELData->pOutlierY == NULL) && (nOutlierCells>0))
	{
		printf("Error: Affy_LoadCEL_CmdCslv1, cannot allocate enough memory for loading outlier cells!\n");
		exit(EXIT_FAILURE);
	}

	pCELData->nModifiedCells = 0;
	pCELData->pModifiedX = NULL;
	pCELData->pModifiedY = NULL;
	pCELData->pModifiedOrig = NULL;

	/* load intensity */
	tempptr1 = pGenericData->vDataGroup[0]->vDataSets[0]->vData;
	tempptr2 = pGenericData->vDataGroup[0]->vDataSets[1]->vData;
	tempptr3 = pGenericData->vDataGroup[0]->vDataSets[2]->vData;
	for(ni=0; ni<nNumberCells; ni++)
	{
		memcpy(fI, tempptr1, FLOAT_SIZE);
		tempptr1 += FLOAT_SIZE;

		memcpy(fI+1, tempptr2, FLOAT_SIZE);
		tempptr2 += FLOAT_SIZE;

	    memcpy(&nCount, tempptr3, SHORT_SIZE);
		tempptr3 += SHORT_SIZE;

		pCELData->pIntensity->pMatElement[ni] = fI[0];
		pCELData->pSD->pMatElement[ni] = fI[1];
		pCELData->pPixelNum->pMatElement[ni] = nCount;
	}

	/* load masked entries */
	if(nMaskedCells > 0)
	{
		tempptr1 = pGenericData->vDataGroup[0]->vDataSets[4]->vData;
		for(ni=0; ni<(int)nMaskedCells; ni++)
		{
			memcpy(nEntry, tempptr1, SHORT_SIZE);
			tempptr1 += SHORT_SIZE;

			memcpy(nEntry+1, tempptr1, SHORT_SIZE);
			tempptr1 += SHORT_SIZE;

			pCELData->pMaskedX->pMatElement[ni] = nEntry[0];
			pCELData->pMaskedY->pMatElement[ni] = nEntry[1];
		}
	}

	/* load outlier entries */
	if(nOutlierCells > 0)
	{
		tempptr1 = pGenericData->vDataGroup[0]->vDataSets[3]->vData;
		for(ni=0; ni<(int)nOutlierCells; ni++)
		{
			memcpy(nEntry, tempptr1, SHORT_SIZE);
			tempptr1 += SHORT_SIZE;

			memcpy(nEntry+1, tempptr1, SHORT_SIZE);
			tempptr1 += SHORT_SIZE;

			pCELData->pOutlierX->pMatElement[ni] = nEntry[0];
			pCELData->pOutlierY->pMatElement[ni] = nEntry[1];
		}
	}

	/* load sub-grids */
	pCELData->vSubGrids = NULL;
	
	/* delete generic data object */
	Affy_GenericFileObject_Delete(&pGenericData);

	/* return */
	return pCELData;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_LoadCELv3()                                                       */
/*  Loading raw data from a single affymetrix's *.CEL (v3) file            */
/* ----------------------------------------------------------------------- */ 
struct tagCELData *Affy_LoadCELv3(char strFileName[])
{
	/* define */
	struct tagCELData *pCELData = NULL;

	FILE *fpCel;
	
	int nCols, nRows;
	int nNumberCells;

	unsigned int nOutlierCells;
	unsigned int nMaskedCells;
	unsigned int nModifiedCells;

	int ni;

	int nULx,nURx,nLRx,nLLx,nULy,nURy,nLRy,nLLy;
	
	char strLine[MED_LINE_LENGTH];
	char strAlgorithm[MED_LINE_LENGTH];
	char strAlgorithmParameters[MED_LINE_LENGTH];

	int nX,nY,nP,nidx;
	double dM,dS;
	char *chSep;

	/* init */

	/* create CELData object */
	pCELData = Affy_CELData_Create();
	if(pCELData == NULL)
	{
		printf("Error: Affy_LoadCELv3, cannot create CELData structure.\n");
		return NULL;
	}

	pCELData->nMagicnumber = 64;
	
	/* load *.CEL */
	fpCel = NULL;
	fpCel = fopen(strFileName, "r");
	if(fpCel == NULL)
	{
		printf("Error: Affy_LoadCELv3, cannot open *.CEL file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, MED_LINE_LENGTH, fpCel)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strstr(strLine, "[CEL]") == strLine)
		{
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			chSep = strchr(strLine, '=');
			chSep++;
			pCELData->nVersionnumber = atoi(chSep);
		}
		else if(strstr(strLine, "[HEADER]") == strLine)
		{
			/* cols */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StringAddTail(&(pCELData->vHeader),strLine);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "Cols") != strLine)
			{
				printf("Error: Affy_LoadCELv3, cannot load Cols!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			nCols = atoi(chSep);
			pCELData->nCols = nCols;

			/* rows */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StringAddTail(&(pCELData->vHeader),strLine);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "Rows") != strLine)
			{
				printf("Error: Affy_LoadCELv3, cannot load Rows!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			nRows = atoi(chSep);
			pCELData->nRows = nRows;

			/* total x */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StringAddTail(&(pCELData->vHeader),strLine);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "TotalX") != strLine)
			{
				printf("Error: Affy_LoadCELv3, cannot load TotalX!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			pCELData->nTotalX = atoi(chSep);

			/* total y */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StringAddTail(&(pCELData->vHeader),strLine);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "TotalY") != strLine)
			{
				printf("Error: Affy_LoadCELv3, cannot load TotalY!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			pCELData->nTotalY = atoi(chSep);

			/* offset x */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StringAddTail(&(pCELData->vHeader),strLine);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "OffsetX") != strLine)
			{
				printf("Error: Affy_LoadCELv3, cannot load OffsetX!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			pCELData->nOffsetX = atoi(chSep);

			/* offset y */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StringAddTail(&(pCELData->vHeader),strLine);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "OffsetY") != strLine)
			{
				printf("Error: Affy_LoadCELv3, cannot load OffsetY!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			pCELData->nOffsetY = atoi(chSep);

			/* UL */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StringAddTail(&(pCELData->vHeader),strLine);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "GridCornerUL") != strLine)
			{
				printf("Error: Affy_LoadCELv3, cannot load GridCornerUL!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			sscanf(chSep, "%d %d", &nULx, &nULy);
			pCELData->nUL[0] = nULx;
			pCELData->nUL[1] = nULy;

			/* UR */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StringAddTail(&(pCELData->vHeader),strLine);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "GridCornerUR") != strLine)
			{
				printf("Error: Affy_LoadCELv3, cannot load GridCornerUR!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			sscanf(chSep, "%d %d", &nURx, &nURy);
			pCELData->nUR[0] = nURx;
			pCELData->nUR[1] = nURy;

			/* LR */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StringAddTail(&(pCELData->vHeader),strLine);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "GridCornerLR") != strLine)
			{
				printf("Error: Affy_LoadCELv3, cannot load GridCornerLR!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			sscanf(chSep, "%d %d", &nLRx, &nLRy);
			pCELData->nLR[0] = nLRx;
			pCELData->nLR[1] = nLRy;

			/* LL */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StringAddTail(&(pCELData->vHeader),strLine);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "GridCornerLL") != strLine)
			{
				printf("Error: Affy_LoadCELv3, cannot load GridCornerLL!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			sscanf(chSep, "%d %d", &nLLx, &nLLy);
			pCELData->nLL[0] = nLLx;
			pCELData->nLL[1] = nLLy;

			/* invertx */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StringAddTail(&(pCELData->vHeader),strLine);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "Axis") != strLine)
			{
				printf("Error: Affy_LoadCELv3, cannot load InvertX!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			pCELData->nInvertX = atoi(chSep);

			/* inverty */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StringAddTail(&(pCELData->vHeader),strLine);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "Axis") != strLine)
			{
				printf("Error: Affy_LoadCELv3, cannot load InvertY!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			pCELData->nInvertY = atoi(chSep);

			/* swapxy */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StringAddTail(&(pCELData->vHeader),strLine);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "swapXY") != strLine)
			{
				printf("Error: Affy_LoadCELv3, cannot load swapXY!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			pCELData->swapXY = atoi(chSep);

			/* DatHeader */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StringAddTail(&(pCELData->vHeader),strLine);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "DatHeader") != strLine)
			{
				printf("Error: Affy_LoadCELv3, cannot load DatHeader!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			strcpy(pCELData->strDatHeader, chSep);

			/* Algorithm */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StringAddTail(&(pCELData->vHeader),strLine);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "Algorithm") != strLine)
			{
				printf("Error: Affy_LoadCELv3, cannot load Algorithm!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			strcpy(strAlgorithm, chSep);
			pCELData->nAlgorithmLength = (int)strlen(strAlgorithm);
			StringAddTail(&(pCELData->vAlgorithm), strAlgorithm);

			/* Algorithm Parameters */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StringAddTail(&(pCELData->vHeader),strLine);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "AlgorithmParameters") != strLine)
			{
				printf("Error: Affy_LoadCELv3, cannot load Algorithm Parameters!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			strcpy(strAlgorithmParameters, chSep);
			pCELData->nAlgorithmParamLength = (int)strlen(strAlgorithmParameters);
			StringAddTail(&(pCELData->vAlgorithmParam), strAlgorithmParameters);

			pCELData->nHeaderLength = pCELData->vHeader->m_nLength;
		}

		else if(strstr(strLine, "[INTENSITY]") == strLine)
		{
			/* NumberCells */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "NumberCells") != strLine)
			{
				printf("Error: Affy_LoadCELv3, cannot load NumberCells!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			nNumberCells = atoi(chSep);
			pCELData->nNumberCells = nNumberCells;
			if(pCELData->nNumberCells != pCELData->nCols*pCELData->nRows)
			{
				printf("Error: Affy_LoadCELv3, NumberCells != Cols*Rows!\n");
				exit(EXIT_FAILURE);
			}

			/* cell header */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "CellHeader") != strLine)
			{
				printf("Error: Affy_LoadCELv3, cannot load CellHeader!\n");
				exit(EXIT_FAILURE);
			}

			/* allocate memory */
			pCELData->pIntensity = NULL;
			pCELData->pIntensity = CreateDoubleMatrix(1, nNumberCells);
			if((pCELData->pIntensity == NULL) && (nNumberCells>0))
			{
				printf("Error: Affy_LoadCELv3, cannot allocate enough memory for loading intensities!\n");
				exit(EXIT_FAILURE);
			}

			pCELData->pSD = NULL;
			pCELData->pSD = CreateDoubleMatrix(1, nNumberCells);
			if((pCELData->pSD == NULL) && (nNumberCells>0))
			{
				printf("Error: Affy_LoadCELv3, cannot allocate enough memory for loading SD!\n");
				exit(EXIT_FAILURE);
			}

			pCELData->pPixelNum = NULL;
			pCELData->pPixelNum = CreateIntMatrix(1, nNumberCells);
			if((pCELData->pPixelNum == NULL) && (nNumberCells>0))
			{
				printf("Error: Affy_LoadCELv3, cannot allocate enough memory for loading pixel numbers!\n");
				exit(EXIT_FAILURE);
			}

			/* load intensity */
			for(ni=0; ni<nNumberCells; ni++)
			{
				fgets(strLine, MED_LINE_LENGTH, fpCel);
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				sscanf(strLine, "%d %d %lf %lf %d", &nX, &nY, &dM, &dS, &nP);
				
				nidx = nY*(pCELData->nTotalX)+nX;
				pCELData->pIntensity->pMatElement[nidx] = dM;
				pCELData->pSD->pMatElement[nidx] = dS;
				pCELData->pPixelNum->pMatElement[nidx] = nP;
			}
		}

		else if(strstr(strLine, "[MASKS]") == strLine)
		{
			/* NumberCells */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "NumberCells") != strLine)
			{
				printf("Error: Affy_LoadCELv3, cannot load MaskedCells!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			nMaskedCells = atoi(chSep);
			pCELData->nMaskedCells = nMaskedCells;

			/* cell header */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "CellHeader") != strLine)
			{
				printf("Error: Affy_LoadCELv3, cannot load MASK CellHeader!\n");
				exit(EXIT_FAILURE);
			}

			pCELData->pMaskedX = NULL;
			pCELData->pMaskedX = CreateIntMatrix(1, nMaskedCells);
			if((pCELData->pMaskedX == NULL) && (nMaskedCells>0))
			{
				printf("Error: Affy_LoadCELv3, cannot allocate enough memory for loading masked cells!\n");
				exit(EXIT_FAILURE);
			}

			pCELData->pMaskedY = NULL;
			pCELData->pMaskedY = CreateIntMatrix(1, nMaskedCells);
			if((pCELData->pMaskedY == NULL) && (nMaskedCells>0))
			{
				printf("Error: Affy_LoadCELv3, cannot allocate enough memory for loading masked cells!\n");
				exit(EXIT_FAILURE);
			}

			/* load cells */
			for(ni=0; ni<(int)nMaskedCells; ni++)
			{
				fgets(strLine, MED_LINE_LENGTH, fpCel);
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				sscanf(strLine, "%d %d", &nX, &nY);
				
				pCELData->pMaskedX->pMatElement[ni] = nX;
				pCELData->pMaskedY->pMatElement[ni] = nY;
			}
		}

		else if(strstr(strLine, "[OUTLIERS]") == strLine)
		{
			/* NumberCells */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "NumberCells") != strLine)
			{
				printf("Error: Affy_LoadCELv3, cannot load OutlierCells!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			nOutlierCells = atoi(chSep);
			pCELData->nOutlierCells = nOutlierCells;


			/* cell header */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "CellHeader") != strLine)
			{
				printf("Error: Affy_LoadCELv3, cannot load OUTLIER CellHeader!\n");
				exit(EXIT_FAILURE);
			}

			pCELData->pOutlierX = NULL;
			pCELData->pOutlierX = CreateIntMatrix(1, nOutlierCells);
			if((pCELData->pOutlierX == NULL) && (nOutlierCells>0))
			{
				printf("Error: Affy_LoadCELv3, cannot allocate enough memory for loading outlier cells!\n");
				exit(EXIT_FAILURE);
			}

			pCELData->pOutlierY = NULL;
			pCELData->pOutlierY = CreateIntMatrix(1, nOutlierCells);
			if((pCELData->pOutlierY == NULL) && (nOutlierCells>0))
			{
				printf("Error: Affy_LoadCELv3, cannot allocate enough memory for loading outlier cells!\n");
				exit(EXIT_FAILURE);
			}

			/* load cells */
			for(ni=0; ni<(int)nOutlierCells; ni++)
			{
				fgets(strLine, MED_LINE_LENGTH, fpCel);
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				sscanf(strLine, "%d %d", &nX, &nY);
				
				pCELData->pOutlierX->pMatElement[ni] = nX;
				pCELData->pOutlierY->pMatElement[ni] = nY;
			}
		}

		else if(strstr(strLine, "[MODIFIED]") == strLine)
		{
			/* NumberCells */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "NumberCells") != strLine)
			{
				printf("Error: Affy_LoadCELv3, cannot load ModifiedCells!\n");
				exit(EXIT_FAILURE);
			}
			chSep = strchr(strLine, '=');
			chSep++;
			nModifiedCells = atoi(chSep);
			pCELData->nModifiedCells = nModifiedCells;

			/* cell header */
			fgets(strLine, MED_LINE_LENGTH, fpCel);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "CellHeader") != strLine)
			{
				printf("Error: Affy_LoadCELv3, cannot load MODIFIED CellHeader!\n");
				exit(EXIT_FAILURE);
			}

			pCELData->pModifiedX = NULL;
			pCELData->pModifiedX = CreateIntMatrix(1, nModifiedCells);
			if((pCELData->pModifiedX == NULL) && (nModifiedCells>0))
			{
				printf("Error: Affy_LoadCELv3, cannot allocate enough memory for loading modified cells!\n");
				exit(EXIT_FAILURE);
			}

			pCELData->pModifiedY = NULL;
			pCELData->pModifiedY = CreateIntMatrix(1, nModifiedCells);
			if((pCELData->pModifiedY == NULL) && (nModifiedCells>0))
			{
				printf("Error: Affy_LoadCELv3, cannot allocate enough memory for loading modified cells!\n");
				exit(EXIT_FAILURE);
			}

			pCELData->pModifiedOrig = NULL;
			pCELData->pModifiedOrig = CreateDoubleMatrix(1, nModifiedCells);
			if((pCELData->pModifiedOrig == NULL) && (nModifiedCells>0))
			{
				printf("Error: Affy_LoadCELv3, cannot allocate enough memory for loading modified cells!\n");
				exit(EXIT_FAILURE);
			}

			/* load cells */
			for(ni=0; ni<(int)nModifiedCells; ni++)
			{
				fgets(strLine, MED_LINE_LENGTH, fpCel);
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				sscanf(strLine, "%d %d %lf", &nX, &nY, &dM);
				pCELData->pModifiedX->pMatElement[ni] = nX;
				pCELData->pModifiedY->pMatElement[ni] = nY;
				pCELData->pModifiedOrig->pMatElement[ni] = dM;
			}
		}

		else
		{
		}
	}

	fclose(fpCel);

	/* return */
	return pCELData;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_SaveCELv3()                                                       */
/*  Exporting raw data from a single affymetrix's *.CEL (v3) file          */
/* ----------------------------------------------------------------------- */ 
int Affy_SaveCELv3(char strFileName[], struct tagCELData *pCELData)
{
	/* define */
	FILE *fpCel;
	
	int ni;
	int nX,nY;
	double *pIntensity,*pSD;
	int *pPixel;

	/* init */

	/* initial check */
	if(pCELData == NULL)
	{
		printf("Warning: empty CEL data!\n");
		return PROC_SUCCESS;
	}
	
	/* open *.CEL */
	fpCel = NULL;
	fpCel = fopen(strFileName, "w");
	if(fpCel == NULL)
	{
		printf("Error: Affy_SaveCELv3, cannot open *.CEL file!\n");
		exit(EXIT_FAILURE);
	}

	/* write [CEL] */
	fprintf(fpCel, "[CEL]\n");
	fprintf(fpCel, "Version=3\n");
	fprintf(fpCel, "\n");

	/* write [HEADER] */
	fprintf(fpCel, "[HEADER]\n");
	fprintf(fpCel, "%s\n", pCELData->vHeader->m_pString);
	/*if(pCELData->nVersionnumber == 3)
	{
		fprintf(fpCel, "[HEADER]\n");
		fprintf(fpCel, "Cols=%d\n", pCELData->nCols);
		fprintf(fpCel, "Rows=%d\n", pCELData->nRows);
		fprintf(fpCel, "TotalX=%d\n", pCELData->nTotalX);
		fprintf(fpCel, "TotalY=%d\n", pCELData->nTotalY);
		fprintf(fpCel, "OffsetX=%d\n", pCELData->nOffsetX);
		fprintf(fpCel, "OffsetY=%d\n", pCELData->nOffsetY);
		fprintf(fpCel, "GridCornerUL=%d %d\n", pCELData->nUL[0], pCELData->nUL[1]);
		fprintf(fpCel, "GridCornerUR=%d %d\n", pCELData->nUR[0], pCELData->nUR[1]);
		fprintf(fpCel, "GridCornerLR=%d %d\n", pCELData->nLR[0], pCELData->nLR[1]);
		fprintf(fpCel, "GridCornerLL=%d %d\n", pCELData->nLL[0], pCELData->nLL[1]);
		fprintf(fpCel, "Axis-invertX=%d\n", pCELData->nInvertX);
		fprintf(fpCel, "AxisInvertY=%d\n", pCELData->nInvertY);
		fprintf(fpCel, "swapXY=%d\n", pCELData->swapXY);
		fprintf(fpCel, "DatHeader=%s\n", pCELData->strDatHeader);
		fprintf(fpCel, "Algorithm=%s\n", pCELData->vAlgorithm->m_pString);
		fprintf(fpCel, "AlgorithmParameters=%s\n", pCELData->vAlgorithmParam->m_pString);
		fprintf(fpCel, "\n");
	}
	else if(pCELData->nVersionnumber == 4)
	{
		fprintf(fpCel, "[HEADER]\n");
		fprintf(fpCel, "%s\n", pCELData->vHeader->m_pString);
	}*/

	/* write [INTENSITY]*/
	fprintf(fpCel, "[INTENSITY]\n");
	fprintf(fpCel, "NumberCells=%d\n", pCELData->nNumberCells);
	fprintf(fpCel, "CellHeader=X\tY\tMEAN\tSTDV\tNPIXELS\n");
	pIntensity = pCELData->pIntensity->pMatElement;
	pSD = pCELData->pSD->pMatElement;
	pPixel = pCELData->pPixelNum->pMatElement;
	for(nY=0; nY<pCELData->nRows; nY++)
	{
		for(nX=0; nX<pCELData->nCols; nX++)
		{
			fprintf(fpCel, "%3d\t%3d\t%f\t%f\t%3d\n", nX, nY,
				*pIntensity, *pSD, *pPixel);
			pIntensity++;
			pSD++;
			pPixel++;
		}
	}
	fprintf(fpCel, "\n");
	
	/* write [MASKS]*/
	fprintf(fpCel, "[MASKS]\n");
	fprintf(fpCel, "NumberCells=%d\n", pCELData->nMaskedCells);
	fprintf(fpCel, "CellHeader=X\tY\n");
	for(ni=0; ni<(int)(pCELData->nMaskedCells); ni++)
	{
		fprintf(fpCel, "%d\t%d\n", pCELData->pMaskedX->pMatElement[ni],
			pCELData->pMaskedY->pMatElement[ni]);
	}
	fprintf(fpCel, "\n");

	/* write [OUTLIERS]*/
	fprintf(fpCel, "[OUTLIERS]\n");
	fprintf(fpCel, "NumberCells=%d\n", pCELData->nOutlierCells);
	fprintf(fpCel, "CellHeader=X\tY\n");
	for(ni=0; ni<(int)(pCELData->nOutlierCells); ni++)
	{
		fprintf(fpCel, "%d\t%d\n", pCELData->pOutlierX->pMatElement[ni],
			pCELData->pOutlierY->pMatElement[ni]);
	}
	fprintf(fpCel, "\n");

	/* write [MODIFIED]*/
	fprintf(fpCel, "[MODIFIED]\n");
	fprintf(fpCel, "NumberCells=%d\n", pCELData->nModifiedCells);
	fprintf(fpCel, "CellHeader=X\tY\tORIGMEAN\n");
	for(ni=0; ni<(int)(pCELData->nModifiedCells); ni++)
	{
		fprintf(fpCel, "%d\t%d\t%f\n", pCELData->pModifiedX->pMatElement[ni],
			pCELData->pModifiedY->pMatElement[ni], 
			pCELData->pModifiedOrig->pMatElement[ni]);
	}
	
	fclose(fpCel);

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  Affy_BARSeq_Create()                                                   */
/*  create a BARSeq object.                                                */
/* ----------------------------------------------------------------------- */ 
struct tagBARSeq *Affy_BARSeq_Create()
{
	/* define */
	struct tagBARSeq *pBARSeq = NULL;

	/* create */
	pBARSeq = (struct tagBARSeq *)calloc(1, sizeof(struct tagBARSeq));
	if(pBARSeq == NULL)
	{
		printf("Error: Affy_BARSeq_Create, cannot create BARSeq object!\n");
		return NULL;
	}

	/* return */
	return pBARSeq;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_BARSeq_Destroy()                                                  */
/*  delete a BARSeq object.                                                */
/* ----------------------------------------------------------------------- */ 
void Affy_BARSeq_Destroy(struct tagBARSeq **pBARSeq)
{
	int ni;

	if(pBARSeq != NULL)
	{
		if( (*pBARSeq) != NULL)
		{
			DeleteString((*pBARSeq)->pSeqName);
			DeleteString((*pBARSeq)->pSeqGroupName);
			DeleteString((*pBARSeq)->pSeqVersion);

			for(ni=0; ni<(*pBARSeq)->nParamNum; ni++)
			{
				DeleteString((*pBARSeq)->vParamName[ni]);
				DeleteString((*pBARSeq)->vParamValue[ni]);
			}
			free((*pBARSeq)->vParamName);
			free((*pBARSeq)->vParamValue);

			for(ni=0; ni<(*pBARSeq)->nColNum; ni++)
			{
				DestroyDoubleMatrix((*pBARSeq)->vData[ni]);
			}
			free((*pBARSeq)->vData);

			free(*pBARSeq);
			*pBARSeq = NULL;
		}
	}
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_BARData_Create()                                                  */
/*  create a BARData object.                                               */
/* ----------------------------------------------------------------------- */ 
struct tagBARData *Affy_BARData_Create()
{
	/* define */
	struct tagBARData *pBARData = NULL;

	/* create */
	pBARData = (struct tagBARData *)calloc(1, sizeof(struct tagBARData));
	if(pBARData == NULL)
	{
		printf("Error: Affy_BARData_Create, cannot create BARData object!\n");
		return NULL;
	}

	/* return */
	return pBARData;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_BARData_Destroy()                                                 */
/*  delete a BARData object.                                               */
/* ----------------------------------------------------------------------- */ 
void Affy_BARData_Destroy(struct tagBARData **pBARData)
{
	int ni;

	if(pBARData != NULL)
	{
		if( (*pBARData) != NULL)
		{
			DestroyIntMatrix((*pBARData)->pFieldType);
			for(ni=0; ni<(*pBARData)->nParamNum; ni++)
			{
				DeleteString((*pBARData)->vParamName[ni]);
				DeleteString((*pBARData)->vParamValue[ni]);
			}
			free((*pBARData)->vParamName);
			free((*pBARData)->vParamValue);

			for(ni=0; ni<(*pBARData)->nSeqNum; ni++)
			{
				Affy_BARSeq_Destroy((*pBARData)->vSeqData+ni);
			}
			free((*pBARData)->vSeqData);

			free(*pBARData);
			*pBARData = NULL;
		}
	}
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_LoadBAR()                                                         */
/*  Loading raw data from a single affymetrix's *.bar file                 */
/* ----------------------------------------------------------------------- */ 
struct tagBARData *Affy_LoadBAR(char strFileName[])
{
	/* define */
	struct tagBARData *pBARData = NULL;

	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));
	
	FILE *fpIn;
	int ni,nj,nk;
	int nLen;

	double dV;
	float fV;
	int nV;
	short sV;
	char cV;
	unsigned int unV;
	unsigned short usV;
	unsigned char ucV;

	/* init */
	fpIn = NULL;
	fpIn = fopen(strFileName, "rb");
	if(fpIn == NULL)
	{
		printf("Error: Affy_LoadBAR, cannot open *.bar file!\n");
		exit(EXIT_FAILURE);
	}

	pBARData = Affy_BARData_Create();
	if(pBARData == NULL)
	{
		printf("Warning: Affy_LoadBAR, BARData object was not created!\n");
		exit(EXIT_FAILURE);
	}

	/* load magic number */
	if(big_endian_fread(pBARData->strMagicnumber, 1, 8, fpIn, little_endian_machine) != 8)
	{
		printf("Error: Affy_LoadBAR, cannot load magic number.\n");
		exit(EXIT_FAILURE);
	}
	pBARData->strMagicnumber[8] = '\0';
	
	if(strcmp(pBARData->strMagicnumber, "barr\r\n\032\n") != 0)
	{
		printf("Error: Affy_LoadBAR, cannot load magic number correctly.\n");
		exit(EXIT_FAILURE);
	}

	/* load version number */
	if(big_endian_fread(&(pBARData->fVersionnumber), FLOAT_SIZE, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: Affy_LoadBAR, cannot load version number.\n");
		exit(EXIT_FAILURE);
	}

	/* load sequence number */
	if(big_endian_fread(&(pBARData->nSeqNum), INT_SIZE, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: Affy_LoadBAR, cannot load sequence number.\n");
		exit(EXIT_FAILURE);
	}

	/* load column number */
	if(big_endian_fread(&(pBARData->nColNum), INT_SIZE, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: Affy_LoadBAR, cannot load column number.\n");
		exit(EXIT_FAILURE);
	}

	/* load column type */
	if(pBARData->nColNum>0)
	{
		pBARData->pFieldType = CreateIntMatrix(1, pBARData->nColNum);
		if(pBARData->pFieldType == NULL)
		{
			printf("Error: Affy_LoadBAR, cannot allocate memory for field type.\n");
			exit(EXIT_FAILURE);
		}
	
		if((int)big_endian_fread(pBARData->pFieldType->pMatElement, INT_SIZE, pBARData->nColNum, fpIn, little_endian_machine) != pBARData->nColNum)
		{
			printf("Error: Affy_LoadBAR, cannot load field type correctly.\n");
			exit(EXIT_FAILURE);
		}
	}

	/* load parameter name/value pairs */
	if(big_endian_fread(&(pBARData->nParamNum), INT_SIZE, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: Affy_LoadBAR, cannot load number of name/value pairs.\n");
		exit(EXIT_FAILURE);
	}

	if(pBARData->nParamNum > 0)
	{
		pBARData->vParamName = (struct tagString **)calloc(pBARData->nParamNum, sizeof(struct tagString *));
		pBARData->vParamValue = (struct tagString **)calloc(pBARData->nParamNum, sizeof(struct tagString *));

		if( (pBARData->vParamName == NULL) || (pBARData->vParamValue == NULL) )
		{
			printf("Error: Affy_LoadBAR, cannot allocate memory for loading parameter name/value pairs\n");
			exit(EXIT_FAILURE);
		}

		for(ni=0; ni<pBARData->nParamNum; ni++)
		{
			if(big_endian_fread(&nLen, INT_SIZE, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: Affy_LoadBAR, cannot load parameter name length.\n");
				exit(EXIT_FAILURE);
			}

			if(nLen > 0)
			{
				pBARData->vParamName[ni] = CreateString(nLen);
				if(pBARData->vParamName[ni] == NULL)
				{
					printf("Error: Affy_LoadBAR, cannot load parameter name.\n");
					exit(EXIT_FAILURE);
				}
			
				if(big_endian_fread(pBARData->vParamName[ni]->m_pString, 1, nLen, fpIn, little_endian_machine) != nLen)
				{
					printf("Error: Affy_LoadBAR, cannot load parameter name.\n");
					exit(EXIT_FAILURE);
				}
				pBARData->vParamName[ni]->m_pString[nLen] = '\0';
			}

			if(big_endian_fread(&nLen, INT_SIZE, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: Affy_LoadBAR, cannot load parameter value length.\n");
				exit(EXIT_FAILURE);
			}

			if(nLen > 0)
			{
				pBARData->vParamValue[ni] = CreateString(nLen);
				if(pBARData->vParamValue[ni] == NULL)
				{
					printf("Error: Affy_LoadBAR, cannot load parameter value.\n");
					exit(EXIT_FAILURE);
				}

				if(big_endian_fread(pBARData->vParamValue[ni]->m_pString, 1, nLen, fpIn, little_endian_machine) != nLen)
				{
					printf("Error: Affy_LoadBAR, cannot load parameter value.\n");
					exit(EXIT_FAILURE);
				}
				pBARData->vParamValue[ni]->m_pString[nLen] = '\0';
			}
		}
	}

	/* load sequence data */
	if(pBARData->nSeqNum <= 0)
		return pBARData;

	pBARData->vSeqData = (struct tagBARSeq **)calloc(pBARData->nSeqNum, sizeof(struct tagBARSeq *));
	if(pBARData->vSeqData == NULL)
	{
		printf("Error: Affy_LoadBAR, cannot allocate memory for loading sequence data.\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<pBARData->nSeqNum; ni++)
	{
		/* create BARSeq object */
		pBARData->vSeqData[ni] = Affy_BARSeq_Create();
		if(pBARData->vSeqData[ni] == NULL)
		{
			printf("Error: Affy_LoadBAR, cannot create BARSeq object.\n");
			exit(EXIT_FAILURE);
		}
		pBARData->vSeqData[ni]->nColNum = pBARData->nColNum;


		/* load sequence name */
		if(big_endian_fread(&nLen, INT_SIZE, 1, fpIn, little_endian_machine) != 1)
		{
			printf("Error: Affy_LoadBAR, cannot load sequence name length.\n");
			exit(EXIT_FAILURE);
		}

		if(nLen > 0)
		{
			pBARData->vSeqData[ni]->pSeqName = CreateString(nLen);
			if(pBARData->vSeqData[ni]->pSeqName == NULL)
			{
				printf("Error: Affy_LoadBAR, cannot load sequence name.\n");
				exit(EXIT_FAILURE);
			}
		
			if(big_endian_fread(pBARData->vSeqData[ni]->pSeqName->m_pString, 1, nLen, fpIn, little_endian_machine) != nLen)
			{
				printf("Error: Affy_LoadBAR, cannot load sequence name.\n");
				exit(EXIT_FAILURE);
			}
			pBARData->vSeqData[ni]->pSeqName->m_pString[nLen] = '\0';
		}

		/* load group name */
		if(pBARData->fVersionnumber > 1.5)
		{
			if(big_endian_fread(&nLen, INT_SIZE, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: Affy_LoadBAR, cannot load sequence group name length.\n");
				exit(EXIT_FAILURE);
			}

			if(nLen > 0)
			{
				pBARData->vSeqData[ni]->pSeqGroupName = CreateString(nLen);
				if(pBARData->vSeqData[ni]->pSeqGroupName == NULL)
				{
					printf("Error: Affy_LoadBAR, cannot load sequence group name.\n");
					exit(EXIT_FAILURE);
				}
			
				if(big_endian_fread(pBARData->vSeqData[ni]->pSeqGroupName->m_pString, 1, nLen, fpIn, little_endian_machine) != nLen)
				{
					printf("Error: Affy_LoadBAR, cannot load sequence group name.\n");
					exit(EXIT_FAILURE);
				}
				pBARData->vSeqData[ni]->pSeqGroupName->m_pString[nLen] = '\0';
			}
		}

		/* load sequence version */
		if(big_endian_fread(&nLen, INT_SIZE, 1, fpIn, little_endian_machine) != 1)
		{
			printf("Error: Affy_LoadBAR, cannot load sequence version.\n");
			exit(EXIT_FAILURE);
		}

		if(nLen > 0)
		{
			pBARData->vSeqData[ni]->pSeqVersion = CreateString(nLen);
			if(pBARData->vSeqData[ni]->pSeqVersion == NULL)
			{
				printf("Error: Affy_LoadBAR, cannot load sequence version.\n");
				exit(EXIT_FAILURE);
			}
		
			if(big_endian_fread(pBARData->vSeqData[ni]->pSeqVersion->m_pString, 1, nLen, fpIn, little_endian_machine) != nLen)
			{
				printf("Error: Affy_LoadBAR, cannot load sequence version.\n");
				exit(EXIT_FAILURE);
			}
			pBARData->vSeqData[ni]->pSeqVersion->m_pString[nLen] = '\0';
		}

		/* load parameter name/value pairs */
		if(pBARData->fVersionnumber > 1.5)
		{
			if(big_endian_fread(&(pBARData->vSeqData[ni]->nParamNum), INT_SIZE, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: Affy_LoadBAR, cannot load number of name/value pairs for BARseq object.\n");
				exit(EXIT_FAILURE);
			}

			if(pBARData->vSeqData[ni]->nParamNum > 0)
			{
				pBARData->vSeqData[ni]->vParamName = (struct tagString **)calloc(pBARData->vSeqData[ni]->nParamNum, sizeof(struct tagString *));
				pBARData->vSeqData[ni]->vParamValue = (struct tagString **)calloc(pBARData->vSeqData[ni]->nParamNum, sizeof(struct tagString *));

				if( (pBARData->vSeqData[ni]->vParamName == NULL) || (pBARData->vSeqData[ni]->vParamValue == NULL) )
				{
					printf("Error: Affy_LoadBAR, cannot allocate memory for loading sequence parameter name/value pairs\n");
					exit(EXIT_FAILURE);
				}

				for(nj=0; nj<pBARData->vSeqData[ni]->nParamNum; nj++)
				{
					if(big_endian_fread(&nLen, INT_SIZE, 1, fpIn, little_endian_machine) != 1)
					{
						printf("Error: Affy_LoadBAR, cannot load sequence parameter name length.\n");
						exit(EXIT_FAILURE);
					}

					if(nLen > 0)
					{
						pBARData->vSeqData[ni]->vParamName[nj] = CreateString(nLen);
						if(pBARData->vSeqData[ni]->vParamName[nj] == NULL)
						{
							printf("Error: Affy_LoadBAR, cannot load sequence parameter name.\n");
							exit(EXIT_FAILURE);
						}

						if(big_endian_fread(pBARData->vSeqData[ni]->vParamName[nj]->m_pString, 1, nLen, fpIn, little_endian_machine) != nLen)
						{
							printf("Error: Affy_LoadBAR, cannot load sequence parameter name.\n");
							exit(EXIT_FAILURE);
						}
						pBARData->vSeqData[ni]->vParamName[nj]->m_pString[nLen] = '\0';
					}

					if(big_endian_fread(&nLen, INT_SIZE, 1, fpIn, little_endian_machine) != 1)
					{
						printf("Error: Affy_LoadBAR, cannot load sequence parameter value length.\n");
						exit(EXIT_FAILURE);
					}

					if(nLen > 0)
					{					
						pBARData->vSeqData[ni]->vParamValue[nj] = CreateString(nLen);
						if(pBARData->vSeqData[ni]->vParamValue[nj] == NULL)
						{
							printf("Error: Affy_LoadBAR, cannot load sequence parameter value.\n");
							exit(EXIT_FAILURE);
						}

						if(big_endian_fread(pBARData->vSeqData[ni]->vParamValue[nj]->m_pString, 1, nLen, fpIn, little_endian_machine) != nLen)
						{
							printf("Error: Affy_LoadBAR, cannot load parameter value.\n");
							exit(EXIT_FAILURE);
						}
						pBARData->vSeqData[ni]->vParamValue[nj]->m_pString[nLen] = '\0';
					}
				}
			}
		}

		/* load data points */
		if(big_endian_fread(&(pBARData->vSeqData[ni]->nDataNum), INT_SIZE, 1, fpIn, little_endian_machine) != 1)
		{
			printf("Error: Affy_LoadBAR, cannot load number of data points for a sequence.\n");
			exit(EXIT_FAILURE);
		}

		if(pBARData->vSeqData[ni]->nColNum > 0)
		{
			pBARData->vSeqData[ni]->vData = (struct DOUBLEMATRIX **)calloc(pBARData->vSeqData[ni]->nColNum, sizeof(struct DOUBLEMATRIX *));
			if(pBARData->vSeqData[ni]->vData == NULL)
			{
				printf("Error: Affy_LoadBAR, cannot allocate memory for loading sequence data.\n");
				exit(EXIT_FAILURE);
			}

			if(pBARData->vSeqData[ni]->nDataNum > 0)
			{
				for(nj=0; nj<pBARData->nColNum; nj++)
				{
					pBARData->vSeqData[ni]->vData[nj] = CreateDoubleMatrix(1, pBARData->vSeqData[ni]->nDataNum);
					if(pBARData->vSeqData[ni]->vData[nj] == NULL)
					{
						printf("Error: Affy_LoadBAR, cannot allocate memory for loading sequence data.\n");
						exit(EXIT_FAILURE);
					}
				}

				for(nk=0; nk<pBARData->vSeqData[ni]->nDataNum; nk++)
				{
					for(nj=0; nj<pBARData->nColNum; nj++)
					{
						switch(pBARData->pFieldType->pMatElement[nj])
						{
							case 0: big_endian_fread(&dV, DOUBLE_SIZE, 1, fpIn, little_endian_machine);
								pBARData->vSeqData[ni]->vData[nj]->pMatElement[nk] = dV;
								break;
							case 1: big_endian_fread(&fV, FLOAT_SIZE, 1, fpIn, little_endian_machine);
								pBARData->vSeqData[ni]->vData[nj]->pMatElement[nk] = fV;
								break;
							case 2: big_endian_fread(&nV, INT_SIZE, 1, fpIn, little_endian_machine);
								pBARData->vSeqData[ni]->vData[nj]->pMatElement[nk] = nV;
								break;
							case 3: big_endian_fread(&sV, SHORT_SIZE, 1, fpIn, little_endian_machine);
								pBARData->vSeqData[ni]->vData[nj]->pMatElement[nk] = sV;
								break;
							case 4: big_endian_fread(&cV, CHAR_SIZE, 1, fpIn, little_endian_machine);
								pBARData->vSeqData[ni]->vData[nj]->pMatElement[nk] = cV;
								break;
							case 5: big_endian_fread(&unV, DWORD_SIZE, 1, fpIn, little_endian_machine);
								pBARData->vSeqData[ni]->vData[nj]->pMatElement[nk] = unV;
								break;
							case 6: big_endian_fread(&usV, USHORT_SIZE, 1, fpIn, little_endian_machine);
								pBARData->vSeqData[ni]->vData[nj]->pMatElement[nk] = usV;
								break;
							case 7: big_endian_fread(&ucV, UCHAR_SIZE, 1, fpIn, little_endian_machine);
								pBARData->vSeqData[ni]->vData[nj]->pMatElement[nk] = ucV;
								break;
						}
					}
				}
			}
		}
	}

	/* close file */
	fclose(fpIn);

	/* return */
	return pBARData;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_LoadBAR_Fast()                                                    */
/*  Loading raw data from a single affymetrix's *.bar file                 */
/* ----------------------------------------------------------------------- */ 
struct tagBARData *Affy_LoadBAR_Fast(char strFileName[])
{
	/* define */
	struct tagBARData *pBARData = NULL;

	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));
	size_t unitsize;
	char *tempmem;
	char *tempptr;
	
	FILE *fpIn;
	int ni,nj,nk;
	int nLen;

	double dV;
	float fV;
	int nV;
	short sV;
	char cV;
	unsigned int unV;
	unsigned short usV;
	unsigned char ucV;

	/* init */
	fpIn = NULL;
	fpIn = fopen(strFileName, "rb");
	if(fpIn == NULL)
	{
		printf("Error: Affy_LoadBAR, cannot open *.bar file!\n");
		exit(EXIT_FAILURE);
	}

	pBARData = Affy_BARData_Create();
	if(pBARData == NULL)
	{
		printf("Warning: Affy_LoadBAR, BARData object was not created!\n");
		exit(EXIT_FAILURE);
	}

	/* load magic number */
	if(big_endian_fread(pBARData->strMagicnumber, 1, 8, fpIn, little_endian_machine) != 8)
	{
		printf("Error: Affy_LoadBAR, cannot load magic number.\n");
		exit(EXIT_FAILURE);
	}
	pBARData->strMagicnumber[8] = '\0';
	
	if(strcmp(pBARData->strMagicnumber, "barr\r\n\032\n") != 0)
	{
		printf("Error: Affy_LoadBAR, cannot load magic number correctly.\n");
		exit(EXIT_FAILURE);
	}

	/* load version number */
	if(big_endian_fread(&(pBARData->fVersionnumber), FLOAT_SIZE, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: Affy_LoadBAR, cannot load version number.\n");
		exit(EXIT_FAILURE);
	}

	/* load sequence number */
	if(big_endian_fread(&(pBARData->nSeqNum), INT_SIZE, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: Affy_LoadBAR, cannot load sequence number.\n");
		exit(EXIT_FAILURE);
	}

	/* load column number */
	if(big_endian_fread(&(pBARData->nColNum), INT_SIZE, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: Affy_LoadBAR, cannot load column number.\n");
		exit(EXIT_FAILURE);
	}

	/* load column type */
	if(pBARData->nColNum>0)
	{
		pBARData->pFieldType = CreateIntMatrix(1, pBARData->nColNum);
		if(pBARData->pFieldType == NULL)
		{
			printf("Error: Affy_LoadBAR, cannot allocate memory for field type.\n");
			exit(EXIT_FAILURE);
		}
	
		if((int)big_endian_fread(pBARData->pFieldType->pMatElement, INT_SIZE, pBARData->nColNum, fpIn, little_endian_machine) != pBARData->nColNum)
		{
			printf("Error: Affy_LoadBAR, cannot load field type correctly.\n");
			exit(EXIT_FAILURE);
		}
	}

	/* load parameter name/value pairs */
	if(big_endian_fread(&(pBARData->nParamNum), INT_SIZE, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: Affy_LoadBAR, cannot load number of name/value pairs.\n");
		exit(EXIT_FAILURE);
	}

	if(pBARData->nParamNum > 0)
	{
		pBARData->vParamName = (struct tagString **)calloc(pBARData->nParamNum, sizeof(struct tagString *));
		pBARData->vParamValue = (struct tagString **)calloc(pBARData->nParamNum, sizeof(struct tagString *));

		if( (pBARData->vParamName == NULL) || (pBARData->vParamValue == NULL) )
		{
			printf("Error: Affy_LoadBAR, cannot allocate memory for loading parameter name/value pairs\n");
			exit(EXIT_FAILURE);
		}

		for(ni=0; ni<pBARData->nParamNum; ni++)
		{
			if(big_endian_fread(&nLen, INT_SIZE, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: Affy_LoadBAR, cannot load parameter name length.\n");
				exit(EXIT_FAILURE);
			}

			if(nLen > 0)
			{
				pBARData->vParamName[ni] = CreateString(nLen);
				if(pBARData->vParamName[ni] == NULL)
				{
					printf("Error: Affy_LoadBAR, cannot load parameter name.\n");
					exit(EXIT_FAILURE);
				}
			
				if(big_endian_fread(pBARData->vParamName[ni]->m_pString, 1, nLen, fpIn, little_endian_machine) != nLen)
				{
					printf("Error: Affy_LoadBAR, cannot load parameter name.\n");
					exit(EXIT_FAILURE);
				}
				pBARData->vParamName[ni]->m_pString[nLen] = '\0';
			}

			if(big_endian_fread(&nLen, INT_SIZE, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: Affy_LoadBAR, cannot load parameter value length.\n");
				exit(EXIT_FAILURE);
			}

			if(nLen > 0)
			{
				pBARData->vParamValue[ni] = CreateString(nLen);
				if(pBARData->vParamValue[ni] == NULL)
				{
					printf("Error: Affy_LoadBAR, cannot load parameter value.\n");
					exit(EXIT_FAILURE);
				}

				if(big_endian_fread(pBARData->vParamValue[ni]->m_pString, 1, nLen, fpIn, little_endian_machine) != nLen)
				{
					printf("Error: Affy_LoadBAR, cannot load parameter value.\n");
					exit(EXIT_FAILURE);
				}
				pBARData->vParamValue[ni]->m_pString[nLen] = '\0';
			}
		}
	}

	/* load sequence data */
	if(pBARData->nSeqNum <= 0)
		return pBARData;

	pBARData->vSeqData = (struct tagBARSeq **)calloc(pBARData->nSeqNum, sizeof(struct tagBARSeq *));
	if(pBARData->vSeqData == NULL)
	{
		printf("Error: Affy_LoadBAR, cannot allocate memory for loading sequence data.\n");
		exit(EXIT_FAILURE);
	}

	unitsize = 0;
	for(nj=0; nj<pBARData->nColNum; nj++)
	{
		switch(pBARData->pFieldType->pMatElement[nj])
		{
			case 0: unitsize += sizeof(double);
				break;
			case 1: unitsize += sizeof(float);
				break;
			case 2: unitsize += 4;
				break;
			case 3: unitsize += 2;
				break;
			case 4: unitsize += 1;
				break;
			case 5: unitsize += 4;
				break;
			case 6: unitsize += 2;
				break;
			case 7: unitsize += 1;
				break;
		}
	}

	for(ni=0; ni<pBARData->nSeqNum; ni++)
	{
		/* create BARSeq object */
		pBARData->vSeqData[ni] = Affy_BARSeq_Create();
		if(pBARData->vSeqData[ni] == NULL)
		{
			printf("Error: Affy_LoadBAR, cannot create BARSeq object.\n");
			exit(EXIT_FAILURE);
		}
		pBARData->vSeqData[ni]->nColNum = pBARData->nColNum;


		/* load sequence name */
		if(big_endian_fread(&nLen, INT_SIZE, 1, fpIn, little_endian_machine) != 1)
		{
			printf("Error: Affy_LoadBAR, cannot load sequence name length.\n");
			exit(EXIT_FAILURE);
		}

		if(nLen > 0)
		{
			pBARData->vSeqData[ni]->pSeqName = CreateString(nLen);
			if(pBARData->vSeqData[ni]->pSeqName == NULL)
			{
				printf("Error: Affy_LoadBAR, cannot load sequence name.\n");
				exit(EXIT_FAILURE);
			}
		
			if(big_endian_fread(pBARData->vSeqData[ni]->pSeqName->m_pString, 1, nLen, fpIn, little_endian_machine) != nLen)
			{
				printf("Error: Affy_LoadBAR, cannot load sequence name.\n");
				exit(EXIT_FAILURE);
			}
			pBARData->vSeqData[ni]->pSeqName->m_pString[nLen] = '\0';
		}

		/* load group name */
		if(pBARData->fVersionnumber > 1.5)
		{
			if(big_endian_fread(&nLen, INT_SIZE, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: Affy_LoadBAR, cannot load sequence group name length.\n");
				exit(EXIT_FAILURE);
			}

			if(nLen > 0)
			{
				pBARData->vSeqData[ni]->pSeqGroupName = CreateString(nLen);
				if(pBARData->vSeqData[ni]->pSeqGroupName == NULL)
				{
					printf("Error: Affy_LoadBAR, cannot load sequence group name.\n");
					exit(EXIT_FAILURE);
				}
			
				if(big_endian_fread(pBARData->vSeqData[ni]->pSeqGroupName->m_pString, 1, nLen, fpIn, little_endian_machine) != nLen)
				{
					printf("Error: Affy_LoadBAR, cannot load sequence group name.\n");
					exit(EXIT_FAILURE);
				}
				pBARData->vSeqData[ni]->pSeqGroupName->m_pString[nLen] = '\0';
			}
		}

		/* load sequence version */
		if(big_endian_fread(&nLen, INT_SIZE, 1, fpIn, little_endian_machine) != 1)
		{
			printf("Error: Affy_LoadBAR, cannot load sequence version.\n");
			exit(EXIT_FAILURE);
		}

		if(nLen > 0)
		{
			pBARData->vSeqData[ni]->pSeqVersion = CreateString(nLen);
			if(pBARData->vSeqData[ni]->pSeqVersion == NULL)
			{
				printf("Error: Affy_LoadBAR, cannot load sequence version.\n");
				exit(EXIT_FAILURE);
			}
		
			if(big_endian_fread(pBARData->vSeqData[ni]->pSeqVersion->m_pString, 1, nLen, fpIn, little_endian_machine) != nLen)
			{
				printf("Error: Affy_LoadBAR, cannot load sequence version.\n");
				exit(EXIT_FAILURE);
			}
			pBARData->vSeqData[ni]->pSeqVersion->m_pString[nLen] = '\0';
		}

		/* load parameter name/value pairs */
		if(pBARData->fVersionnumber > 1.5)
		{
			if(big_endian_fread(&(pBARData->vSeqData[ni]->nParamNum), INT_SIZE, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: Affy_LoadBAR, cannot load number of name/value pairs for BARseq object.\n");
				exit(EXIT_FAILURE);
			}

			if(pBARData->vSeqData[ni]->nParamNum > 0)
			{
				pBARData->vSeqData[ni]->vParamName = (struct tagString **)calloc(pBARData->vSeqData[ni]->nParamNum, sizeof(struct tagString *));
				pBARData->vSeqData[ni]->vParamValue = (struct tagString **)calloc(pBARData->vSeqData[ni]->nParamNum, sizeof(struct tagString *));

				if( (pBARData->vSeqData[ni]->vParamName == NULL) || (pBARData->vSeqData[ni]->vParamValue == NULL) )
				{
					printf("Error: Affy_LoadBAR, cannot allocate memory for loading sequence parameter name/value pairs\n");
					exit(EXIT_FAILURE);
				}

				for(nj=0; nj<pBARData->vSeqData[ni]->nParamNum; nj++)
				{
					if(big_endian_fread(&nLen, INT_SIZE, 1, fpIn, little_endian_machine) != 1)
					{
						printf("Error: Affy_LoadBAR, cannot load sequence parameter name length.\n");
						exit(EXIT_FAILURE);
					}

					if(nLen > 0)
					{
						pBARData->vSeqData[ni]->vParamName[nj] = CreateString(nLen);
						if(pBARData->vSeqData[ni]->vParamName[nj] == NULL)
						{
							printf("Error: Affy_LoadBAR, cannot load sequence parameter name.\n");
							exit(EXIT_FAILURE);
						}

						if(big_endian_fread(pBARData->vSeqData[ni]->vParamName[nj]->m_pString, 1, nLen, fpIn, little_endian_machine) != nLen)
						{
							printf("Error: Affy_LoadBAR, cannot load sequence parameter name.\n");
							exit(EXIT_FAILURE);
						}
						pBARData->vSeqData[ni]->vParamName[nj]->m_pString[nLen] = '\0';
					}

					if(big_endian_fread(&nLen, INT_SIZE, 1, fpIn, little_endian_machine) != 1)
					{
						printf("Error: Affy_LoadBAR, cannot load sequence parameter value length.\n");
						exit(EXIT_FAILURE);
					}

					if(nLen > 0)
					{					
						pBARData->vSeqData[ni]->vParamValue[nj] = CreateString(nLen);
						if(pBARData->vSeqData[ni]->vParamValue[nj] == NULL)
						{
							printf("Error: Affy_LoadBAR, cannot load sequence parameter value.\n");
							exit(EXIT_FAILURE);
						}

						if(big_endian_fread(pBARData->vSeqData[ni]->vParamValue[nj]->m_pString, 1, nLen, fpIn, little_endian_machine) != nLen)
						{
							printf("Error: Affy_LoadBAR, cannot load parameter value.\n");
							exit(EXIT_FAILURE);
						}
						pBARData->vSeqData[ni]->vParamValue[nj]->m_pString[nLen] = '\0';
					}
				}
			}
		}

		/* load data points */
		if(big_endian_fread(&(pBARData->vSeqData[ni]->nDataNum), INT_SIZE, 1, fpIn, little_endian_machine) != 1)
		{
			printf("Error: Affy_LoadBAR, cannot load number of data points for a sequence.\n");
			exit(EXIT_FAILURE);
		}

		if(pBARData->vSeqData[ni]->nColNum > 0)
		{
			pBARData->vSeqData[ni]->vData = (struct DOUBLEMATRIX **)calloc(pBARData->vSeqData[ni]->nColNum, sizeof(struct DOUBLEMATRIX *));
			if(pBARData->vSeqData[ni]->vData == NULL)
			{
				printf("Error: Affy_LoadBAR, cannot allocate memory for loading sequence data.\n");
				exit(EXIT_FAILURE);
			}

			if(pBARData->vSeqData[ni]->nDataNum > 0)
			{
				for(nj=0; nj<pBARData->nColNum; nj++)
				{
					pBARData->vSeqData[ni]->vData[nj] = CreateDoubleMatrix(1, pBARData->vSeqData[ni]->nDataNum);
					if(pBARData->vSeqData[ni]->vData[nj] == NULL)
					{
						printf("Error: Affy_LoadBAR, cannot allocate memory for loading sequence data.\n");
						exit(EXIT_FAILURE);
					}
				}

				tempmem = NULL;
				tempmem = (char *)calloc(pBARData->vSeqData[ni]->nDataNum, unitsize);
				if(tempmem == NULL)
				{
					printf("Error: Affy_LoadBAR, cannot allocate memory for loading sequence data.\n");
					exit(EXIT_FAILURE);
				}

				if( (int)fread(tempmem, unitsize, pBARData->vSeqData[ni]->nDataNum, fpIn) != pBARData->vSeqData[ni]->nDataNum)
				{
					printf("Error: Affy_LoadBAR, cannot load sequence data correctly.\n");
					exit(EXIT_FAILURE);
				}
				tempptr = tempmem;

				for(nk=0; nk<pBARData->vSeqData[ni]->nDataNum; nk++)
				{
					for(nj=0; nj<pBARData->nColNum; nj++)
					{
						switch(pBARData->pFieldType->pMatElement[nj])
						{
							case 0: 
								if(little_endian_machine == 1)
								{
									reverse_buf(tempptr, sizeof(double));
								}
								memcpy(&dV, tempptr, sizeof(double));
								tempptr += sizeof(double);
								pBARData->vSeqData[ni]->vData[nj]->pMatElement[nk] = dV;
								break;
							case 1: 
								if(little_endian_machine == 1)
								{
									reverse_buf(tempptr, sizeof(float));
								}
								memcpy(&fV, tempptr, sizeof(float));
								tempptr += sizeof(float);
								pBARData->vSeqData[ni]->vData[nj]->pMatElement[nk] = fV;
								break;
							case 2: 
								if(little_endian_machine == 1)
								{
									reverse_buf(tempptr, 4);
								}
								memcpy(&nV, tempptr, 4);
								tempptr += 4;
								pBARData->vSeqData[ni]->vData[nj]->pMatElement[nk] = nV;
								break;
							case 3: 
								if(little_endian_machine == 1)
								{
									reverse_buf(tempptr, 2);
								}
								memcpy(&sV, tempptr, 2);
								tempptr += 2;
								pBARData->vSeqData[ni]->vData[nj]->pMatElement[nk] = sV;
								break;
							case 4: 
								if(little_endian_machine == 1)
								{
									reverse_buf(tempptr, 1);
								}
								memcpy(&cV, tempptr, 1);
								tempptr += 1;
								pBARData->vSeqData[ni]->vData[nj]->pMatElement[nk] = cV;
								break;
							case 5: 
								if(little_endian_machine == 1)
								{
									reverse_buf(tempptr, 4);
								}
								memcpy(&unV, tempptr, 4);
								tempptr += 4;
								pBARData->vSeqData[ni]->vData[nj]->pMatElement[nk] = unV;
								break;
							case 6: 
								if(little_endian_machine == 1)
								{
									reverse_buf(tempptr, 2);
								}
								memcpy(&usV, tempptr, 2);
								tempptr += 2;
								pBARData->vSeqData[ni]->vData[nj]->pMatElement[nk] = usV;
								break;
							case 7: 
								if(little_endian_machine == 1)
								{
									reverse_buf(tempptr, 1);
								}
								memcpy(&ucV, tempptr, 1);
								tempptr += 1;
								pBARData->vSeqData[ni]->vData[nj]->pMatElement[nk] = ucV;
								break;
						}
					}
				}

				free(tempmem);
			}
		}
	}

	/* close file */
	fclose(fpIn);

	/* return */
	return pBARData;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_SaveBAR()                                                         */
/*  Saving raw data to a affymetrix's *.bar file                           */
/* ----------------------------------------------------------------------- */ 
int Affy_SaveBAR(char strFileName[], struct tagBARData *pBARData)
{
	/* define */
	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));
	
	FILE *fpOut;
	int ni,nj,nk;
	int nLen;

	double dV;
	float fV;
	int nV;
	short sV;
	char cV;
	unsigned int unV;
	unsigned short usV;
	unsigned char ucV;

	/* initial check */
	if(pBARData == NULL)
	{
		printf("Warning: empty BAR data!\n");
		return PROC_SUCCESS;
	}
	
	/* load *.CEL */
	fpOut = NULL;
	fpOut = fopen(strFileName, "wb");
	if(fpOut == NULL)
	{
		printf("Error: Affy_SaveBAR, cannot open output file!\n");
		return PROC_FAILURE;
	}

	/* write magic number */
	if(strcmp(pBARData->strMagicnumber, "barr\r\n\032\n") != 0)
	{
		printf("Error: Affy_SaveBAR, magic number wrong.\n");
		exit(EXIT_FAILURE);
	}

	if(big_endian_fwrite(pBARData->strMagicnumber, 1, 8, fpOut, little_endian_machine) != 8)
	{
		printf("Error: Affy_SaveBAR, cannot write magic number.\n");
		exit(EXIT_FAILURE);
	}
	
	/* write version number */
	if(big_endian_fwrite(&(pBARData->fVersionnumber), FLOAT_SIZE, 1, fpOut, little_endian_machine) != 1)
	{
		printf("Error: Affy_SaveBAR, cannot write version number.\n");
		exit(EXIT_FAILURE);
	}

	/* write sequence number */
	if(big_endian_fwrite(&(pBARData->nSeqNum), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
	{
		printf("Error: Affy_SaveBAR, cannot write sequence number.\n");
		exit(EXIT_FAILURE);
	}

	/* write column number */
	if(big_endian_fwrite(&(pBARData->nColNum), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
	{
		printf("Error: Affy_SaveBAR, cannot write column number.\n");
		exit(EXIT_FAILURE);
	}

	/* write column type */
	if(pBARData->nColNum>0)
	{
		if(pBARData->pFieldType == NULL)
		{
			printf("Error: Affy_SaveBAR, cannot find field type.\n");
			exit(EXIT_FAILURE);
		}

		if((int)big_endian_fwrite(pBARData->pFieldType->pMatElement, INT_SIZE, pBARData->nColNum, fpOut, little_endian_machine) != pBARData->nColNum)
		{
			printf("Error: Affy_SaveBAR, cannot write field type correctly.\n");
			exit(EXIT_FAILURE);
		}
	}

	/* write parameter name/value pairs */
	if(big_endian_fwrite(&(pBARData->nParamNum), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
	{
		printf("Error: Affy_SaveBAR, cannot write number of name/value pairs.\n");
		exit(EXIT_FAILURE);
	}

	if(pBARData->nParamNum > 0)
	{
		if( (pBARData->vParamName == NULL) || (pBARData->vParamValue == NULL) )
		{
			printf("Error: Affy_SaveBAR, cannot find parameter name/value pairs\n");
			exit(EXIT_FAILURE);
		}

		for(ni=0; ni<pBARData->nParamNum; ni++)
		{
			if(pBARData->vParamName[ni] == NULL)
			{
				nLen = 0;
				if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOut, little_endian_machine) != 1)
				{
					printf("Error: Affy_SaveBAR, cannot write parameter name length.\n");
					exit(EXIT_FAILURE);
				}
			}
			else
			{
				if(big_endian_fwrite(&(pBARData->vParamName[ni]->m_nLength), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
				{
					printf("Error: Affy_SaveBAR, cannot write parameter name length.\n");
					exit(EXIT_FAILURE);
				}

				if(big_endian_fwrite(pBARData->vParamName[ni]->m_pString, 1, pBARData->vParamName[ni]->m_nLength, fpOut, little_endian_machine) != pBARData->vParamName[ni]->m_nLength)
				{
					printf("Error: Affy_SaveBAR, cannot write parameter name.\n");
					exit(EXIT_FAILURE);
				}
			}

			if(pBARData->vParamValue[ni] == NULL)
			{
				nLen = 0;
				if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOut, little_endian_machine) != 1)
				{
					printf("Error: Affy_SaveBAR, cannot write parameter name length.\n");
					exit(EXIT_FAILURE);
				}
			}
			else
			{
				if(big_endian_fwrite(&(pBARData->vParamValue[ni]->m_nLength), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
				{
					printf("Error: Affy_SaveBAR, cannot write parameter value length.\n");
					exit(EXIT_FAILURE);
				}

				if(big_endian_fwrite(pBARData->vParamValue[ni]->m_pString, 1, pBARData->vParamValue[ni]->m_nLength, fpOut, little_endian_machine) != pBARData->vParamValue[ni]->m_nLength)
				{
					printf("Error: Affy_SaveBAR, cannot write parameter value.\n");
					exit(EXIT_FAILURE);
				}
			}
		}
	}

	/* write sequence data */
	for(ni=0; ni<pBARData->nSeqNum; ni++)
	{
		/* create BARSeq object */
		if(pBARData->vSeqData[ni] == NULL)
		{
			printf("Error: Affy_SaveBAR, cannot find BARSeq object.\n");
			exit(EXIT_FAILURE);
		}

		/* write sequence name */
		if(pBARData->vSeqData[ni]->pSeqName == NULL)
		{
			nLen = 0;
			if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOut, little_endian_machine) != 1)
			{
				printf("Error: Affy_SaveBAR, cannot write sequence name length.\n");
				exit(EXIT_FAILURE);
			}
		}
		else
		{
			if(big_endian_fwrite(&(pBARData->vSeqData[ni]->pSeqName->m_nLength), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
			{
				printf("Error: Affy_SaveBAR, cannot write sequence name length.\n");
				exit(EXIT_FAILURE);
			}
			
			if(big_endian_fwrite(pBARData->vSeqData[ni]->pSeqName->m_pString, 1, pBARData->vSeqData[ni]->pSeqName->m_nLength, fpOut, little_endian_machine) != pBARData->vSeqData[ni]->pSeqName->m_nLength)
			{
				printf("Error: Affy_SaveBAR, cannot write sequence name.\n");
				exit(EXIT_FAILURE);
			}
		}
		
		/* write group name */
		if(pBARData->fVersionnumber > 1.5)
		{
			if(pBARData->vSeqData[ni]->pSeqGroupName == NULL)
			{
				nLen = 0;
				if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOut, little_endian_machine) != 1)
				{
					printf("Error: Affy_SaveBAR, cannot write sequence group name length.\n");
					exit(EXIT_FAILURE);
				}
			}
			else
			{
				if(big_endian_fwrite(&(pBARData->vSeqData[ni]->pSeqGroupName->m_nLength), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
				{
					printf("Error: Affy_SaveBAR, cannot write sequence group name length.\n");
					exit(EXIT_FAILURE);
				}

				if(big_endian_fwrite(pBARData->vSeqData[ni]->pSeqGroupName->m_pString, 1, pBARData->vSeqData[ni]->pSeqGroupName->m_nLength, fpOut, little_endian_machine) != pBARData->vSeqData[ni]->pSeqGroupName->m_nLength)
				{
					printf("Error: Affy_SaveBAR, cannot write sequence group name.\n");
					exit(EXIT_FAILURE);
				}
			}
		}

		/* write sequence version */
		if(pBARData->vSeqData[ni]->pSeqVersion == NULL)
		{
			nLen = 0;
			if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOut, little_endian_machine) != 1)
			{
				printf("Error: Affy_SaveBAR, cannot write sequence version length.\n");
				exit(EXIT_FAILURE);
			}
		}
		else
		{
			if(big_endian_fwrite(&(pBARData->vSeqData[ni]->pSeqVersion->m_nLength), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
			{
				printf("Error: Affy_SaveBAR, cannot write sequence version length.\n");
				exit(EXIT_FAILURE);
			}
			
			if(big_endian_fwrite(pBARData->vSeqData[ni]->pSeqVersion->m_pString, 1, pBARData->vSeqData[ni]->pSeqVersion->m_nLength, fpOut, little_endian_machine) != pBARData->vSeqData[ni]->pSeqVersion->m_nLength)
			{
				printf("Error: Affy_SaveBAR, cannot write sequence version.\n");
				exit(EXIT_FAILURE);
			}
		}
		
		/* write parameter name/value pairs */
		if(pBARData->fVersionnumber > 1.5)
		{
			if(big_endian_fwrite(&(pBARData->vSeqData[ni]->nParamNum), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
			{
				printf("Error: Affy_SaveBAR, cannot write number of name/value pairs for BARseq object.\n");
				exit(EXIT_FAILURE);
			}

			if(pBARData->vSeqData[ni]->nParamNum > 0)
			{
				if( (pBARData->vSeqData[ni]->vParamName == NULL) || (pBARData->vSeqData[ni]->vParamValue == NULL) )
				{
					printf("Error: Affy_SaveBAR, cannot find sequence parameter name/value pairs\n");
					exit(EXIT_FAILURE);
				}

				for(nj=0; nj<pBARData->vSeqData[ni]->nParamNum; nj++)
				{
					if(pBARData->vSeqData[ni]->vParamName[nj] == NULL)
					{
						nLen = 0;
						if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOut, little_endian_machine) != 1)
						{
							printf("Error: Affy_SaveBAR, cannot write sequence parameter name length.\n");
							exit(EXIT_FAILURE);
						}
					}
					else
					{
						if(big_endian_fwrite(&(pBARData->vSeqData[ni]->vParamName[nj]->m_nLength), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
						{
							printf("Error: Affy_SaveBAR, cannot write sequence parameter name length.\n");
							exit(EXIT_FAILURE);
						}

						if(big_endian_fwrite(pBARData->vSeqData[ni]->vParamName[nj]->m_pString, 1, pBARData->vSeqData[ni]->vParamName[nj]->m_nLength, fpOut, little_endian_machine) != pBARData->vSeqData[ni]->vParamName[nj]->m_nLength)
						{
							printf("Error: Affy_SaveBAR, cannot write sequence parameter name.\n");
							exit(EXIT_FAILURE);
						}
					}
					
					if(pBARData->vSeqData[ni]->vParamValue[nj] == NULL)
					{
						nLen = 0;
						if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOut, little_endian_machine) != 1)
						{
							printf("Error: Affy_SaveBAR, cannot write sequence parameter value length.\n");
							exit(EXIT_FAILURE);
						}
					}
					else
					{
						if(big_endian_fwrite(&(pBARData->vSeqData[ni]->vParamValue[nj]->m_nLength), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
						{
							printf("Error: Affy_SaveBAR, cannot write sequence parameter value length.\n");
							exit(EXIT_FAILURE);
						}

						if(big_endian_fwrite(pBARData->vSeqData[ni]->vParamValue[nj]->m_pString, 1, pBARData->vSeqData[ni]->vParamValue[nj]->m_nLength, fpOut, little_endian_machine) != pBARData->vSeqData[ni]->vParamValue[nj]->m_nLength)
						{
							printf("Error: Affy_SaveBAR, cannot write parameter value.\n");
							exit(EXIT_FAILURE);
						}
					}
				}
			}
		}

		/* write data points */
		if(big_endian_fwrite(&(pBARData->vSeqData[ni]->nDataNum), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
		{
			printf("Error: Affy_SaveBAR, cannot write number of data points for a sequence.\n");
			exit(EXIT_FAILURE);
		}

		if(pBARData->vSeqData[ni]->nColNum > 0)
		{
			if(pBARData->vSeqData[ni]->vData == NULL)
			{
				printf("Error: Affy_SaveBAR, cannot find sequence data.\n");
				exit(EXIT_FAILURE);
			}
		}

		if(pBARData->vSeqData[ni]->nDataNum > 0)
		{
			for(nj=0; nj<pBARData->nColNum; nj++)
			{
				if(pBARData->vSeqData[ni]->vData[nj] == NULL)
				{
					printf("Error: Affy_SaveBAR, cannot find sequence data.\n");
					exit(EXIT_FAILURE);
				}
			}

			for(nk=0; nk<pBARData->vSeqData[ni]->nDataNum; nk++)
			{
				for(nj=0; nj<pBARData->nColNum; nj++)
				{
					switch(pBARData->pFieldType->pMatElement[nj])
					{
						case 0: dV = (double)(pBARData->vSeqData[ni]->vData[nj]->pMatElement[nk]);
							big_endian_fwrite(&dV, DOUBLE_SIZE, 1, fpOut, little_endian_machine);
							break;
						case 1: fV = (float)(pBARData->vSeqData[ni]->vData[nj]->pMatElement[nk]);
							big_endian_fwrite(&fV, FLOAT_SIZE, 1, fpOut, little_endian_machine);
							break;
						case 2: nV = (int)(pBARData->vSeqData[ni]->vData[nj]->pMatElement[nk]);
							big_endian_fwrite(&nV, INT_SIZE, 1, fpOut, little_endian_machine);
							break;
						case 3: sV = (short)(pBARData->vSeqData[ni]->vData[nj]->pMatElement[nk]);
							big_endian_fwrite(&sV, SHORT_SIZE, 1, fpOut, little_endian_machine);
							break;
						case 4: cV = (char)(pBARData->vSeqData[ni]->vData[nj]->pMatElement[nk]);
							big_endian_fwrite(&cV, CHAR_SIZE, 1, fpOut, little_endian_machine);
							break;
						case 5: unV = (unsigned int)(pBARData->vSeqData[ni]->vData[nj]->pMatElement[nk]);
							big_endian_fwrite(&unV, DWORD_SIZE, 1, fpOut, little_endian_machine);
							break;
						case 6: usV = (unsigned short)(pBARData->vSeqData[ni]->vData[nj]->pMatElement[nk]);
							big_endian_fwrite(&usV, USHORT_SIZE, 1, fpOut, little_endian_machine);
							break;
						case 7: ucV = (unsigned char)(pBARData->vSeqData[ni]->vData[nj]->pMatElement[nk]);
							big_endian_fwrite(&ucV, UCHAR_SIZE, 1, fpOut, little_endian_machine);
							break;
					}
				}
			}
		}
	}

	/* close file */
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_SaveBAR_Columns_Fast()                                            */
/*  Saving specified columns of a bar data object to a *.bar file          */
/* ----------------------------------------------------------------------- */ 
int Affy_SaveBAR_Columns_Fast(char strFileName[], struct tagBARData *pBARData,
							  struct INTMATRIX *pCol)
{
	/* define */
	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));
	size_t unitsize;
	char *tempmem;
	char *tempptr;
	
	FILE *fpOut;
	int ni,nj,nk;
	int nLen;

	double dV;
	float fV;
	int nV;
	short sV;
	char cV;
	unsigned int unV;
	unsigned short usV;
	unsigned char ucV;

	int nNewColNum;
	int *vNewFieldType;
	int *vNewColId;
	int nColId;

	/* initial check */
	if( (pBARData == NULL) || (pCol == NULL) )
	{
		printf("Warning: empty BAR/Column data!\n");
		return PROC_SUCCESS;
	}

	/* get column number */
	nNewColNum = 0;
	for(ni=0; ni<pCol->nWidth; ni++)
	{
		if(pCol->pMatElement[ni] == 1)
		{
			nNewColNum += 1;
		}
	}
	if(nNewColNum == 0)
	{
		printf("Warning: Affy_SaveBAR_Columns_Fast, no columns need to be exported!\n");
		return PROC_SUCCESS;
	}
	
	vNewFieldType = NULL;
	vNewFieldType = (int *)calloc(nNewColNum, sizeof(int));
	vNewColId = NULL;
	vNewColId = (int *)calloc(nNewColNum, sizeof(int));
	if( (vNewFieldType == NULL) || (vNewColId == NULL) )
	{
		printf("Error: Affy_SaveBAR_Columns_Fast, cannot allocate memory for saving columns!\n");
		return PROC_SUCCESS;
	}
	nj = 0;
	unitsize = 0;
	for(ni=0; ni<pCol->nWidth; ni++)
	{
		if(pCol->pMatElement[ni] == 1)
		{
			vNewColId[nj] = ni;
			vNewFieldType[nj] = pBARData->pFieldType->pMatElement[ni];
			switch(vNewFieldType[nj])
			{
				case 0: unitsize += sizeof(double);
					break;
				case 1: unitsize += sizeof(float);
					break;
				case 2: unitsize += 4;
					break;
				case 3: unitsize += 2;
					break;
				case 4: unitsize += 1;
					break;
				case 5: unitsize += 4;
					break;
				case 6: unitsize += 2;
					break;
				case 7: unitsize += 1;
					break;
			}
			nj++;
		}
	}

	/* open file */
	fpOut = NULL;
	fpOut = fopen(strFileName, "wb");
	if(fpOut == NULL)
	{
		printf("Error: Affy_SaveBAR_Columns_Fast, cannot open output file!\n");
		return PROC_FAILURE;
	}

	/* write magic number */
	if(strcmp(pBARData->strMagicnumber, "barr\r\n\032\n") != 0)
	{
		printf("Error: Affy_SaveBAR_Columns_Fast, magic number wrong.\n");
		exit(EXIT_FAILURE);
	}

	if(big_endian_fwrite(pBARData->strMagicnumber, 1, 8, fpOut, little_endian_machine) != 8)
	{
		printf("Error: Affy_SaveBAR, cannot write magic number.\n");
		exit(EXIT_FAILURE);
	}
	
	/* write version number */
	if(big_endian_fwrite(&(pBARData->fVersionnumber), FLOAT_SIZE, 1, fpOut, little_endian_machine) != 1)
	{
		printf("Error: Affy_SaveBAR_Columns_Fast, cannot write version number.\n");
		exit(EXIT_FAILURE);
	}

	/* write sequence number */
	if(big_endian_fwrite(&(pBARData->nSeqNum), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
	{
		printf("Error: Affy_SaveBAR_Columns_Fast, cannot write sequence number.\n");
		exit(EXIT_FAILURE);
	}

	/* write column number */
	if(big_endian_fwrite(&(nNewColNum), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
	{
		printf("Error: Affy_SaveBAR_Columns_Fast, cannot write column number.\n");
		exit(EXIT_FAILURE);
	}

	/* write column type */
	if(nNewColNum>0)
	{
		if((int)big_endian_fwrite(vNewFieldType, INT_SIZE, nNewColNum, fpOut, little_endian_machine) != nNewColNum)
		{
			printf("Error: Affy_SaveBAR_Columns_Fast, cannot write field type correctly.\n");
			exit(EXIT_FAILURE);
		}
	}

	/* write parameter name/value pairs */
	if(big_endian_fwrite(&(pBARData->nParamNum), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
	{
		printf("Error: Affy_SaveBAR_Columns_Fast, cannot write number of name/value pairs.\n");
		exit(EXIT_FAILURE);
	}

	if(pBARData->nParamNum > 0)
	{
		if( (pBARData->vParamName == NULL) || (pBARData->vParamValue == NULL) )
		{
			printf("Error: Affy_SaveBAR_Columns_Fast, cannot find parameter name/value pairs\n");
			exit(EXIT_FAILURE);
		}

		for(ni=0; ni<pBARData->nParamNum; ni++)
		{
			if(pBARData->vParamName[ni] == NULL)
			{
				nLen = 0;
				if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOut, little_endian_machine) != 1)
				{
					printf("Error: Affy_SaveBAR_Columns_Fast, cannot write parameter name length.\n");
					exit(EXIT_FAILURE);
				}
			}
			else
			{
				if(big_endian_fwrite(&(pBARData->vParamName[ni]->m_nLength), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
				{
					printf("Error: Affy_SaveBAR_Columns_Fast, cannot write parameter name length.\n");
					exit(EXIT_FAILURE);
				}

				if(big_endian_fwrite(pBARData->vParamName[ni]->m_pString, 1, pBARData->vParamName[ni]->m_nLength, fpOut, little_endian_machine) != pBARData->vParamName[ni]->m_nLength)
				{
					printf("Error: Affy_SaveBAR_Columns_Fast, cannot write parameter name.\n");
					exit(EXIT_FAILURE);
				}
			}

			if(pBARData->vParamValue[ni] == NULL)
			{
				nLen = 0;
				if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOut, little_endian_machine) != 1)
				{
					printf("Error: Affy_SaveBAR_Columns_Fast, cannot write parameter name length.\n");
					exit(EXIT_FAILURE);
				}
			}
			else
			{
				if(big_endian_fwrite(&(pBARData->vParamValue[ni]->m_nLength), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
				{
					printf("Error: Affy_SaveBAR_Columns_Fast, cannot write parameter value length.\n");
					exit(EXIT_FAILURE);
				}

				if(big_endian_fwrite(pBARData->vParamValue[ni]->m_pString, 1, pBARData->vParamValue[ni]->m_nLength, fpOut, little_endian_machine) != pBARData->vParamValue[ni]->m_nLength)
				{
					printf("Error: Affy_SaveBAR_Columns_Fast, cannot write parameter value.\n");
					exit(EXIT_FAILURE);
				}
			}
		}
	}

	/* write sequence data */
	for(ni=0; ni<pBARData->nSeqNum; ni++)
	{
		/* create BARSeq object */
		if(pBARData->vSeqData[ni] == NULL)
		{
			printf("Error: Affy_SaveBAR_Columns_Fast, cannot find BARSeq object.\n");
			exit(EXIT_FAILURE);
		}

		/* write sequence name */
		if(pBARData->vSeqData[ni]->pSeqName == NULL)
		{
			nLen = 0;
			if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOut, little_endian_machine) != 1)
			{
				printf("Error: Affy_SaveBAR_Columns_Fast, cannot write sequence name length.\n");
				exit(EXIT_FAILURE);
			}
		}
		else
		{
			if(big_endian_fwrite(&(pBARData->vSeqData[ni]->pSeqName->m_nLength), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
			{
				printf("Error: Affy_SaveBAR_Columns_Fast, cannot write sequence name length.\n");
				exit(EXIT_FAILURE);
			}
			
			if(big_endian_fwrite(pBARData->vSeqData[ni]->pSeqName->m_pString, 1, pBARData->vSeqData[ni]->pSeqName->m_nLength, fpOut, little_endian_machine) != pBARData->vSeqData[ni]->pSeqName->m_nLength)
			{
				printf("Error: Affy_SaveBAR_Columns_Fast, cannot write sequence name.\n");
				exit(EXIT_FAILURE);
			}
		}
		
		/* write group name */
		if(pBARData->fVersionnumber > 1.5)
		{
			if(pBARData->vSeqData[ni]->pSeqGroupName == NULL)
			{
				nLen = 0;
				if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOut, little_endian_machine) != 1)
				{
					printf("Error: Affy_SaveBAR_Columns_Fast, cannot write sequence group name length.\n");
					exit(EXIT_FAILURE);
				}
			}
			else
			{
				if(big_endian_fwrite(&(pBARData->vSeqData[ni]->pSeqGroupName->m_nLength), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
				{
					printf("Error: Affy_SaveBAR_Columns_Fast, cannot write sequence group name length.\n");
					exit(EXIT_FAILURE);
				}

				if(big_endian_fwrite(pBARData->vSeqData[ni]->pSeqGroupName->m_pString, 1, pBARData->vSeqData[ni]->pSeqGroupName->m_nLength, fpOut, little_endian_machine) != pBARData->vSeqData[ni]->pSeqGroupName->m_nLength)
				{
					printf("Error: Affy_SaveBAR_Columns_Fast, cannot write sequence group name.\n");
					exit(EXIT_FAILURE);
				}
			}
		}

		/* write sequence version */
		if(pBARData->vSeqData[ni]->pSeqVersion == NULL)
		{
			nLen = 0;
			if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOut, little_endian_machine) != 1)
			{
				printf("Error: Affy_SaveBAR_Columns_Fast, cannot write sequence version length.\n");
				exit(EXIT_FAILURE);
			}
		}
		else
		{
			if(big_endian_fwrite(&(pBARData->vSeqData[ni]->pSeqVersion->m_nLength), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
			{
				printf("Error: Affy_SaveBAR_Columns_Fast, cannot write sequence version length.\n");
				exit(EXIT_FAILURE);
			}
			
			if(big_endian_fwrite(pBARData->vSeqData[ni]->pSeqVersion->m_pString, 1, pBARData->vSeqData[ni]->pSeqVersion->m_nLength, fpOut, little_endian_machine) != pBARData->vSeqData[ni]->pSeqVersion->m_nLength)
			{
				printf("Error: Affy_SaveBAR_Columns_Fast, cannot write sequence version.\n");
				exit(EXIT_FAILURE);
			}
		}
		
		/* write parameter name/value pairs */
		if(pBARData->fVersionnumber > 1.5)
		{
			if(big_endian_fwrite(&(pBARData->vSeqData[ni]->nParamNum), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
			{
				printf("Error: Affy_SaveBAR_Columns_Fast, cannot write number of name/value pairs for BARseq object.\n");
				exit(EXIT_FAILURE);
			}

			if(pBARData->vSeqData[ni]->nParamNum > 0)
			{
				if( (pBARData->vSeqData[ni]->vParamName == NULL) || (pBARData->vSeqData[ni]->vParamValue == NULL) )
				{
					printf("Error: Affy_SaveBAR_Columns_Fast, cannot find sequence parameter name/value pairs\n");
					exit(EXIT_FAILURE);
				}

				for(nj=0; nj<pBARData->vSeqData[ni]->nParamNum; nj++)
				{
					if(pBARData->vSeqData[ni]->vParamName[nj] == NULL)
					{
						nLen = 0;
						if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOut, little_endian_machine) != 1)
						{
							printf("Error: Affy_SaveBAR_Columns_Fast, cannot write sequence parameter name length.\n");
							exit(EXIT_FAILURE);
						}
					}
					else
					{
						if(big_endian_fwrite(&(pBARData->vSeqData[ni]->vParamName[nj]->m_nLength), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
						{
							printf("Error: Affy_SaveBAR_Columns_Fast, cannot write sequence parameter name length.\n");
							exit(EXIT_FAILURE);
						}

						if(big_endian_fwrite(pBARData->vSeqData[ni]->vParamName[nj]->m_pString, 1, pBARData->vSeqData[ni]->vParamName[nj]->m_nLength, fpOut, little_endian_machine) != pBARData->vSeqData[ni]->vParamName[nj]->m_nLength)
						{
							printf("Error: Affy_SaveBAR_Columns_Fast, cannot write sequence parameter name.\n");
							exit(EXIT_FAILURE);
						}
					}
					
					if(pBARData->vSeqData[ni]->vParamValue[nj] == NULL)
					{
						nLen = 0;
						if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOut, little_endian_machine) != 1)
						{
							printf("Error: Affy_SaveBAR_Columns_Fast, cannot write sequence parameter value length.\n");
							exit(EXIT_FAILURE);
						}
					}
					else
					{
						if(big_endian_fwrite(&(pBARData->vSeqData[ni]->vParamValue[nj]->m_nLength), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
						{
							printf("Error: Affy_SaveBAR_Columns_Fast, cannot write sequence parameter value length.\n");
							exit(EXIT_FAILURE);
						}

						if(big_endian_fwrite(pBARData->vSeqData[ni]->vParamValue[nj]->m_pString, 1, pBARData->vSeqData[ni]->vParamValue[nj]->m_nLength, fpOut, little_endian_machine) != pBARData->vSeqData[ni]->vParamValue[nj]->m_nLength)
						{
							printf("Error: Affy_SaveBAR_Columns_Fast, cannot write parameter value.\n");
							exit(EXIT_FAILURE);
						}
					}
				}
			}
		}

		/* write data points */
		if(big_endian_fwrite(&(pBARData->vSeqData[ni]->nDataNum), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
		{
			printf("Error: Affy_SaveBAR_Columns_Fast, cannot write number of data points for a sequence.\n");
			exit(EXIT_FAILURE);
		}

		if(pBARData->vSeqData[ni]->nColNum > 0)
		{
			if(pBARData->vSeqData[ni]->vData == NULL)
			{
				printf("Error: Affy_SaveBAR_Columns_Fast, cannot find sequence data.\n");
				exit(EXIT_FAILURE);
			}
		}

		if(pBARData->vSeqData[ni]->nDataNum > 0)
		{
			for(nj=0; nj<nNewColNum; nj++)
			{
				nColId = vNewColId[nj];
				if(pBARData->vSeqData[ni]->vData[nColId] == NULL)
				{
					printf("Error: Affy_SaveBAR_Columns_Fast, cannot find sequence data.\n");
					exit(EXIT_FAILURE);
				}
			}

			tempmem = NULL;
			tempmem = (char *)calloc(pBARData->vSeqData[ni]->nDataNum, unitsize);
			if(tempmem == NULL)
			{
				printf("Error: Affy_SaveBAR_Columns_Fast, cannot allocate memory for saving sequence data.\n");
				exit(EXIT_FAILURE);
			}

			tempptr = tempmem;
			for(nk=0; nk<pBARData->vSeqData[ni]->nDataNum; nk++)
			{
				for(nj=0; nj<nNewColNum; nj++)
				{
					nColId = vNewColId[nj];
					switch(vNewFieldType[nj])
					{
						case 0:
							dV = (double)(pBARData->vSeqData[ni]->vData[nColId]->pMatElement[nk]);
							memcpy(tempptr, &dV, sizeof(double));
							if(little_endian_machine == 1)
							{
								reverse_buf(tempptr, sizeof(double));
							}
							tempptr += sizeof(double);
							break;
						case 1: 
							fV = (float)(pBARData->vSeqData[ni]->vData[nColId]->pMatElement[nk]);
							memcpy(tempptr, &fV, sizeof(float));
							if(little_endian_machine == 1)
							{
								reverse_buf(tempptr, sizeof(float));
							}
							tempptr += sizeof(float);
							break;
						case 2: 
							nV = (int)(pBARData->vSeqData[ni]->vData[nColId]->pMatElement[nk]);
							memcpy(tempptr, &nV, 4);
							if(little_endian_machine == 1)
							{
								reverse_buf(tempptr, 4);
							}
							tempptr += 4;			
							break;
						case 3: 
							sV = (short)(pBARData->vSeqData[ni]->vData[nColId]->pMatElement[nk]);
							memcpy(tempptr, &sV, 2);
							if(little_endian_machine == 1)
							{
								reverse_buf(tempptr, 2);
							}
							tempptr += 2;	
							break;
						case 4: 
							cV = (char)(pBARData->vSeqData[ni]->vData[nColId]->pMatElement[nk]);
							memcpy(tempptr, &cV, 1);
							if(little_endian_machine == 1)
							{
								reverse_buf(tempptr, 1);
							}
							tempptr += 1;	
							break;
						case 5: 
							unV = (unsigned int)(pBARData->vSeqData[ni]->vData[nColId]->pMatElement[nk]);
							memcpy(tempptr, &unV, 4);
							if(little_endian_machine == 1)
							{
								reverse_buf(tempptr, 4);
							}
							tempptr += 4;
							break;
						case 6: 
							usV = (unsigned short)(pBARData->vSeqData[ni]->vData[nColId]->pMatElement[nk]);
							memcpy(tempptr, &usV, 2);
							if(little_endian_machine == 1)
							{
								reverse_buf(tempptr, 2);
							}
							tempptr += 2;	
							break;
						case 7: 
							ucV = (unsigned char)(pBARData->vSeqData[ni]->vData[nColId]->pMatElement[nk]);
							memcpy(tempptr, &ucV, 1);
							if(little_endian_machine == 1)
							{
								reverse_buf(tempptr, 1);
							}
							tempptr += 1;
							break;
					}
				}
			}

			if( (int)fwrite(tempmem, 1, unitsize*pBARData->vSeqData[ni]->nDataNum, fpOut) != unitsize*pBARData->vSeqData[ni]->nDataNum)
			{
				printf("Error: Affy_SaveBAR_Columns_Fast, cannot save sequence data correctly.\n");
				exit(EXIT_FAILURE);
			}
			free(tempmem);
		}
	}

	/* close file */
	fclose(fpOut);
	free(vNewFieldType);
	free(vNewColId);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_SaveFilteredBAR_Columns_Fast()                                    */
/*  Saving specified columns of a bar data object to a *.bar file          */
/*  Exporting data will be filtered according to GenomeGrp name.           */
/* ----------------------------------------------------------------------- */ 
int Affy_SaveFilteredBAR_Columns_Fast(char strFileName[], struct tagBARData *pBARData,
							  struct INTMATRIX *pCol, char strGenomeGrp[])
{
	/* define */
	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));
	size_t unitsize;
	char *tempmem;
	char *tempptr;
	
	FILE *fpOut;
	int ni,nj,nk;
	int nLen;

	double dV;
	float fV;
	int nV;
	short sV;
	char cV;
	unsigned int unV;
	unsigned short usV;
	unsigned char ucV;

	int nNewColNum;
	int *vNewFieldType;
	int *vNewColId;
	int nColId;

	int nNewSeqNum;

	/* initial check */
	if( (pBARData == NULL) || (pCol == NULL) )
	{
		printf("Warning: empty BAR/Column data!\n");
		return PROC_SUCCESS;
	}

	/* get column number */
	nNewColNum = 0;
	for(ni=0; ni<pCol->nWidth; ni++)
	{
		if(pCol->pMatElement[ni] == 1)
		{
			nNewColNum += 1;
		}
	}
	if(nNewColNum == 0)
	{
		printf("Warning: Affy_SaveBAR_Columns_Fast, no columns need to be exported!\n");
		return PROC_SUCCESS;
	}

	nNewSeqNum = 0;
	for(ni=0; ni<pBARData->nSeqNum; ni++)
	{
		if(strcmp(strGenomeGrp, "") != 0)
		{
			if(pBARData->vSeqData[ni]->pSeqGroupName == NULL)
				continue;

			if(strcmp(strGenomeGrp, pBARData->vSeqData[ni]->pSeqGroupName->m_pString) != 0)
				continue;
		}
		
		nNewSeqNum++;
	}
	
	vNewFieldType = NULL;
	vNewFieldType = (int *)calloc(nNewColNum, sizeof(int));
	vNewColId = NULL;
	vNewColId = (int *)calloc(nNewColNum, sizeof(int));
	if( (vNewFieldType == NULL) || (vNewColId == NULL) )
	{
		printf("Error: Affy_SaveBAR_Columns_Fast, cannot allocate memory for saving columns!\n");
		return PROC_SUCCESS;
	}
	nj = 0;
	unitsize = 0;
	for(ni=0; ni<pCol->nWidth; ni++)
	{
		if(pCol->pMatElement[ni] == 1)
		{
			vNewColId[nj] = ni;
			vNewFieldType[nj] = pBARData->pFieldType->pMatElement[ni];
			switch(vNewFieldType[nj])
			{
				case 0: unitsize += sizeof(double);
					break;
				case 1: unitsize += sizeof(float);
					break;
				case 2: unitsize += 4;
					break;
				case 3: unitsize += 2;
					break;
				case 4: unitsize += 1;
					break;
				case 5: unitsize += 4;
					break;
				case 6: unitsize += 2;
					break;
				case 7: unitsize += 1;
					break;
			}
			nj++;
		}
	}

	/* open file */
	fpOut = NULL;
	fpOut = fopen(strFileName, "wb");
	if(fpOut == NULL)
	{
		printf("Error: Affy_SaveBAR_Columns_Fast, cannot open output file!\n");
		return PROC_FAILURE;
	}

	/* write magic number */
	if(strcmp(pBARData->strMagicnumber, "barr\r\n\032\n") != 0)
	{
		printf("Error: Affy_SaveBAR_Columns_Fast, magic number wrong.\n");
		exit(EXIT_FAILURE);
	}

	if(big_endian_fwrite(pBARData->strMagicnumber, 1, 8, fpOut, little_endian_machine) != 8)
	{
		printf("Error: Affy_SaveBAR, cannot write magic number.\n");
		exit(EXIT_FAILURE);
	}
	
	/* write version number */
	if(big_endian_fwrite(&(pBARData->fVersionnumber), FLOAT_SIZE, 1, fpOut, little_endian_machine) != 1)
	{
		printf("Error: Affy_SaveBAR_Columns_Fast, cannot write version number.\n");
		exit(EXIT_FAILURE);
	}

	/* write sequence number */
	if(big_endian_fwrite(&nNewSeqNum, INT_SIZE, 1, fpOut, little_endian_machine) != 1)
	{
		printf("Error: Affy_SaveBAR_Columns_Fast, cannot write sequence number.\n");
		exit(EXIT_FAILURE);
	}

	/* write column number */
	if(big_endian_fwrite(&(nNewColNum), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
	{
		printf("Error: Affy_SaveBAR_Columns_Fast, cannot write column number.\n");
		exit(EXIT_FAILURE);
	}

	/* write column type */
	if(nNewColNum>0)
	{
		if((int)big_endian_fwrite(vNewFieldType, INT_SIZE, nNewColNum, fpOut, little_endian_machine) != nNewColNum)
		{
			printf("Error: Affy_SaveBAR_Columns_Fast, cannot write field type correctly.\n");
			exit(EXIT_FAILURE);
		}
	}

	/* write parameter name/value pairs */
	if(big_endian_fwrite(&(pBARData->nParamNum), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
	{
		printf("Error: Affy_SaveBAR_Columns_Fast, cannot write number of name/value pairs.\n");
		exit(EXIT_FAILURE);
	}

	if(pBARData->nParamNum > 0)
	{
		if( (pBARData->vParamName == NULL) || (pBARData->vParamValue == NULL) )
		{
			printf("Error: Affy_SaveBAR_Columns_Fast, cannot find parameter name/value pairs\n");
			exit(EXIT_FAILURE);
		}

		for(ni=0; ni<pBARData->nParamNum; ni++)
		{
			if(pBARData->vParamName[ni] == NULL)
			{
				nLen = 0;
				if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOut, little_endian_machine) != 1)
				{
					printf("Error: Affy_SaveBAR_Columns_Fast, cannot write parameter name length.\n");
					exit(EXIT_FAILURE);
				}
			}
			else
			{
				if(big_endian_fwrite(&(pBARData->vParamName[ni]->m_nLength), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
				{
					printf("Error: Affy_SaveBAR_Columns_Fast, cannot write parameter name length.\n");
					exit(EXIT_FAILURE);
				}

				if(big_endian_fwrite(pBARData->vParamName[ni]->m_pString, 1, pBARData->vParamName[ni]->m_nLength, fpOut, little_endian_machine) != pBARData->vParamName[ni]->m_nLength)
				{
					printf("Error: Affy_SaveBAR_Columns_Fast, cannot write parameter name.\n");
					exit(EXIT_FAILURE);
				}
			}

			if(pBARData->vParamValue[ni] == NULL)
			{
				nLen = 0;
				if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOut, little_endian_machine) != 1)
				{
					printf("Error: Affy_SaveBAR_Columns_Fast, cannot write parameter name length.\n");
					exit(EXIT_FAILURE);
				}
			}
			else
			{
				if(big_endian_fwrite(&(pBARData->vParamValue[ni]->m_nLength), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
				{
					printf("Error: Affy_SaveBAR_Columns_Fast, cannot write parameter value length.\n");
					exit(EXIT_FAILURE);
				}

				if(big_endian_fwrite(pBARData->vParamValue[ni]->m_pString, 1, pBARData->vParamValue[ni]->m_nLength, fpOut, little_endian_machine) != pBARData->vParamValue[ni]->m_nLength)
				{
					printf("Error: Affy_SaveBAR_Columns_Fast, cannot write parameter value.\n");
					exit(EXIT_FAILURE);
				}
			}
		}
	}

	/* write sequence data */
	for(ni=0; ni<pBARData->nSeqNum; ni++)
	{
		/* create BARSeq object */
		if(pBARData->vSeqData[ni] == NULL)
		{
			printf("Error: Affy_SaveBAR_Columns_Fast, cannot find BARSeq object.\n");
			exit(EXIT_FAILURE);
		}

		if(strcmp(strGenomeGrp, "") != 0)
		{
			if(pBARData->vSeqData[ni]->pSeqGroupName == NULL)
				continue;

			if(strcmp(strGenomeGrp, pBARData->vSeqData[ni]->pSeqGroupName->m_pString) != 0)
				continue;
		}

		/* write sequence name */
		if(pBARData->vSeqData[ni]->pSeqName == NULL)
		{
			nLen = 0;
			if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOut, little_endian_machine) != 1)
			{
				printf("Error: Affy_SaveBAR_Columns_Fast, cannot write sequence name length.\n");
				exit(EXIT_FAILURE);
			}
		}
		else
		{
			if(big_endian_fwrite(&(pBARData->vSeqData[ni]->pSeqName->m_nLength), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
			{
				printf("Error: Affy_SaveBAR_Columns_Fast, cannot write sequence name length.\n");
				exit(EXIT_FAILURE);
			}
			
			if(big_endian_fwrite(pBARData->vSeqData[ni]->pSeqName->m_pString, 1, pBARData->vSeqData[ni]->pSeqName->m_nLength, fpOut, little_endian_machine) != pBARData->vSeqData[ni]->pSeqName->m_nLength)
			{
				printf("Error: Affy_SaveBAR_Columns_Fast, cannot write sequence name.\n");
				exit(EXIT_FAILURE);
			}
		}
		
		/* write group name */
		if(pBARData->fVersionnumber > 1.5)
		{
			nLen = 0;
			if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOut, little_endian_machine) != 1)
			{
				printf("Error: Affy_SaveBAR_Columns_Fast, cannot write sequence group name length.\n");
				exit(EXIT_FAILURE);
			}

			/*if(pBARData->vSeqData[ni]->pSeqGroupName == NULL)
			{
				nLen = 0;
				if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOut, little_endian_machine) != 1)
				{
					printf("Error: Affy_SaveBAR_Columns_Fast, cannot write sequence group name length.\n");
					exit(EXIT_FAILURE);
				}
			}
			else
			{
				if(big_endian_fwrite(&(pBARData->vSeqData[ni]->pSeqGroupName->m_nLength), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
				{
					printf("Error: Affy_SaveBAR_Columns_Fast, cannot write sequence group name length.\n");
					exit(EXIT_FAILURE);
				}

				if(big_endian_fwrite(pBARData->vSeqData[ni]->pSeqGroupName->m_pString, 1, pBARData->vSeqData[ni]->pSeqGroupName->m_nLength, fpOut, little_endian_machine) != pBARData->vSeqData[ni]->pSeqGroupName->m_nLength)
				{
					printf("Error: Affy_SaveBAR_Columns_Fast, cannot write sequence group name.\n");
					exit(EXIT_FAILURE);
				}
			} */
		}

		/* write sequence version */
		nLen = 0;
		if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOut, little_endian_machine) != 1)
		{
			printf("Error: Affy_SaveBAR_Columns_Fast, cannot write sequence version length.\n");
			exit(EXIT_FAILURE);
		}
		
		/* if(pBARData->vSeqData[ni]->pSeqVersion == NULL)
		{
			nLen = 0;
			if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOut, little_endian_machine) != 1)
			{
				printf("Error: Affy_SaveBAR_Columns_Fast, cannot write sequence version length.\n");
				exit(EXIT_FAILURE);
			}
		}
		else
		{
			if(big_endian_fwrite(&(pBARData->vSeqData[ni]->pSeqVersion->m_nLength), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
			{
				printf("Error: Affy_SaveBAR_Columns_Fast, cannot write sequence version length.\n");
				exit(EXIT_FAILURE);
			}
			
			if(big_endian_fwrite(pBARData->vSeqData[ni]->pSeqVersion->m_pString, 1, pBARData->vSeqData[ni]->pSeqVersion->m_nLength, fpOut, little_endian_machine) != pBARData->vSeqData[ni]->pSeqVersion->m_nLength)
			{
				printf("Error: Affy_SaveBAR_Columns_Fast, cannot write sequence version.\n");
				exit(EXIT_FAILURE);
			}
		} */
		
		/* write parameter name/value pairs */
		if(pBARData->fVersionnumber > 1.5)
		{
			if(big_endian_fwrite(&(pBARData->vSeqData[ni]->nParamNum), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
			{
				printf("Error: Affy_SaveBAR_Columns_Fast, cannot write number of name/value pairs for BARseq object.\n");
				exit(EXIT_FAILURE);
			}

			if(pBARData->vSeqData[ni]->nParamNum > 0)
			{
				if( (pBARData->vSeqData[ni]->vParamName == NULL) || (pBARData->vSeqData[ni]->vParamValue == NULL) )
				{
					printf("Error: Affy_SaveBAR_Columns_Fast, cannot find sequence parameter name/value pairs\n");
					exit(EXIT_FAILURE);
				}

				for(nj=0; nj<pBARData->vSeqData[ni]->nParamNum; nj++)
				{
					if(pBARData->vSeqData[ni]->vParamName[nj] == NULL)
					{
						nLen = 0;
						if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOut, little_endian_machine) != 1)
						{
							printf("Error: Affy_SaveBAR_Columns_Fast, cannot write sequence parameter name length.\n");
							exit(EXIT_FAILURE);
						}
					}
					else
					{
						if(big_endian_fwrite(&(pBARData->vSeqData[ni]->vParamName[nj]->m_nLength), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
						{
							printf("Error: Affy_SaveBAR_Columns_Fast, cannot write sequence parameter name length.\n");
							exit(EXIT_FAILURE);
						}

						if(big_endian_fwrite(pBARData->vSeqData[ni]->vParamName[nj]->m_pString, 1, pBARData->vSeqData[ni]->vParamName[nj]->m_nLength, fpOut, little_endian_machine) != pBARData->vSeqData[ni]->vParamName[nj]->m_nLength)
						{
							printf("Error: Affy_SaveBAR_Columns_Fast, cannot write sequence parameter name.\n");
							exit(EXIT_FAILURE);
						}
					}
					
					if(pBARData->vSeqData[ni]->vParamValue[nj] == NULL)
					{
						nLen = 0;
						if(big_endian_fwrite(&nLen, INT_SIZE, 1, fpOut, little_endian_machine) != 1)
						{
							printf("Error: Affy_SaveBAR_Columns_Fast, cannot write sequence parameter value length.\n");
							exit(EXIT_FAILURE);
						}
					}
					else
					{
						if(big_endian_fwrite(&(pBARData->vSeqData[ni]->vParamValue[nj]->m_nLength), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
						{
							printf("Error: Affy_SaveBAR_Columns_Fast, cannot write sequence parameter value length.\n");
							exit(EXIT_FAILURE);
						}

						if(big_endian_fwrite(pBARData->vSeqData[ni]->vParamValue[nj]->m_pString, 1, pBARData->vSeqData[ni]->vParamValue[nj]->m_nLength, fpOut, little_endian_machine) != pBARData->vSeqData[ni]->vParamValue[nj]->m_nLength)
						{
							printf("Error: Affy_SaveBAR_Columns_Fast, cannot write parameter value.\n");
							exit(EXIT_FAILURE);
						}
					}
				}
			}
		}

		/* write data points */
		if(big_endian_fwrite(&(pBARData->vSeqData[ni]->nDataNum), INT_SIZE, 1, fpOut, little_endian_machine) != 1)
		{
			printf("Error: Affy_SaveBAR_Columns_Fast, cannot write number of data points for a sequence.\n");
			exit(EXIT_FAILURE);
		}

		if(pBARData->vSeqData[ni]->nColNum > 0)
		{
			if(pBARData->vSeqData[ni]->vData == NULL)
			{
				printf("Error: Affy_SaveBAR_Columns_Fast, cannot find sequence data.\n");
				exit(EXIT_FAILURE);
			}
		}

		if(pBARData->vSeqData[ni]->nDataNum > 0)
		{
			for(nj=0; nj<nNewColNum; nj++)
			{
				nColId = vNewColId[nj];
				if(pBARData->vSeqData[ni]->vData[nColId] == NULL)
				{
					printf("Error: Affy_SaveBAR_Columns_Fast, cannot find sequence data.\n");
					exit(EXIT_FAILURE);
				}
			}

			tempmem = NULL;
			tempmem = (char *)calloc(pBARData->vSeqData[ni]->nDataNum, unitsize);
			if(tempmem == NULL)
			{
				printf("Error: Affy_SaveBAR_Columns_Fast, cannot allocate memory for saving sequence data.\n");
				exit(EXIT_FAILURE);
			}

			tempptr = tempmem;
			for(nk=0; nk<pBARData->vSeqData[ni]->nDataNum; nk++)
			{
				for(nj=0; nj<nNewColNum; nj++)
				{
					nColId = vNewColId[nj];
					switch(vNewFieldType[nj])
					{
						case 0:
							dV = (double)(pBARData->vSeqData[ni]->vData[nColId]->pMatElement[nk]);
							memcpy(tempptr, &dV, sizeof(double));
							if(little_endian_machine == 1)
							{
								reverse_buf(tempptr, sizeof(double));
							}
							tempptr += sizeof(double);
							break;
						case 1: 
							fV = (float)(pBARData->vSeqData[ni]->vData[nColId]->pMatElement[nk]);
							memcpy(tempptr, &fV, sizeof(float));
							if(little_endian_machine == 1)
							{
								reverse_buf(tempptr, sizeof(float));
							}
							tempptr += sizeof(float);
							break;
						case 2: 
							nV = (int)(pBARData->vSeqData[ni]->vData[nColId]->pMatElement[nk]);
							memcpy(tempptr, &nV, 4);
							if(little_endian_machine == 1)
							{
								reverse_buf(tempptr, 4);
							}
							tempptr += 4;			
							break;
						case 3: 
							sV = (short)(pBARData->vSeqData[ni]->vData[nColId]->pMatElement[nk]);
							memcpy(tempptr, &sV, 2);
							if(little_endian_machine == 1)
							{
								reverse_buf(tempptr, 2);
							}
							tempptr += 2;	
							break;
						case 4: 
							cV = (char)(pBARData->vSeqData[ni]->vData[nColId]->pMatElement[nk]);
							memcpy(tempptr, &cV, 1);
							if(little_endian_machine == 1)
							{
								reverse_buf(tempptr, 1);
							}
							tempptr += 1;	
							break;
						case 5: 
							unV = (unsigned int)(pBARData->vSeqData[ni]->vData[nColId]->pMatElement[nk]);
							memcpy(tempptr, &unV, 4);
							if(little_endian_machine == 1)
							{
								reverse_buf(tempptr, 4);
							}
							tempptr += 4;
							break;
						case 6: 
							usV = (unsigned short)(pBARData->vSeqData[ni]->vData[nColId]->pMatElement[nk]);
							memcpy(tempptr, &usV, 2);
							if(little_endian_machine == 1)
							{
								reverse_buf(tempptr, 2);
							}
							tempptr += 2;	
							break;
						case 7: 
							ucV = (unsigned char)(pBARData->vSeqData[ni]->vData[nColId]->pMatElement[nk]);
							memcpy(tempptr, &ucV, 1);
							if(little_endian_machine == 1)
							{
								reverse_buf(tempptr, 1);
							}
							tempptr += 1;
							break;
					}
				}
			}

			if( (int)fwrite(tempmem, 1, unitsize*pBARData->vSeqData[ni]->nDataNum, fpOut) != unitsize*pBARData->vSeqData[ni]->nDataNum)
			{
				printf("Error: Affy_SaveBAR_Columns_Fast, cannot save sequence data correctly.\n");
				exit(EXIT_FAILURE);
			}
			free(tempmem);
		}
	}

	/* close file */
	fclose(fpOut);
	free(vNewFieldType);
	free(vNewColId);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_BAR2TXT()                                                         */
/*  Convert *.bar file to *.txt file.                                      */
/* ----------------------------------------------------------------------- */ 
int Affy_BAR2TXT(char strBARFile[], char strTXTFile[])
{
	/* define */
	struct tagBARData *pBARData = NULL;

	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));
	
	FILE *fpIn,*fpOut;
	int ni,nj,nk;
	int nLen;

	double dV;
	float fV;
	int nV;
	short sV;
	char cV;
	unsigned int unV;
	unsigned short usV;
	unsigned char ucV;
	char strSeqAlias[MED_LINE_LENGTH];

	/* init */
	fpIn = NULL;
	fpIn = fopen(strBARFile, "rb");
	if(fpIn == NULL)
	{
		printf("Error: Affy_BAR2TXT, cannot open input *.bar file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strTXTFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: Affy_BAR2TXT, cannot open output *.txt file!\n");
		exit(EXIT_FAILURE);
	}

	pBARData = Affy_BARData_Create();
	if(pBARData == NULL)
	{
		printf("Warning: Affy_BAR2TXT, BARData object was not created!\n");
		exit(EXIT_FAILURE);
	}

	/* load magic number */
	if(big_endian_fread(pBARData->strMagicnumber, 1, 8, fpIn, little_endian_machine) != 8)
	{
		printf("Error: Affy_BAR2TXT, cannot load magic number.\n");
		exit(EXIT_FAILURE);
	}
	pBARData->strMagicnumber[8] = '\0';
	
	if(strcmp(pBARData->strMagicnumber, "barr\r\n\032\n") != 0)
	{
		printf("Error: Affy_BAR2TXT, cannot load magic number correctly.\n");
		exit(EXIT_FAILURE);
	}

	/* load version number */
	if(big_endian_fread(&(pBARData->fVersionnumber), FLOAT_SIZE, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: Affy_BAR2TXT, cannot load version number.\n");
		exit(EXIT_FAILURE);
	}

	/* load sequence number */
	if(big_endian_fread(&(pBARData->nSeqNum), INT_SIZE, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: Affy_BAR2TXT, cannot load sequence number.\n");
		exit(EXIT_FAILURE);
	}

	/* load column number */
	if(big_endian_fread(&(pBARData->nColNum), INT_SIZE, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: Affy_BAR2TXT, cannot load column number.\n");
		exit(EXIT_FAILURE);
	}

	/* load column type */
	if(pBARData->nColNum>0)
	{
		pBARData->pFieldType = CreateIntMatrix(1, pBARData->nColNum);
		if(pBARData->pFieldType == NULL)
		{
			printf("Error: Affy_BAR2TXT, cannot allocate memory for field type.\n");
			exit(EXIT_FAILURE);
		}
	
		if((int)big_endian_fread(pBARData->pFieldType->pMatElement, INT_SIZE, pBARData->nColNum, fpIn, little_endian_machine) != pBARData->nColNum)
		{
			printf("Error: Affy_BAR2TXT, cannot load field type correctly.\n");
			exit(EXIT_FAILURE);
		}
	}

	/* load parameter name/value pairs */
	if(big_endian_fread(&(pBARData->nParamNum), INT_SIZE, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: Affy_BAR2TXT, cannot load number of name/value pairs.\n");
		exit(EXIT_FAILURE);
	}

	if(pBARData->nParamNum > 0)
	{
		pBARData->vParamName = (struct tagString **)calloc(pBARData->nParamNum, sizeof(struct tagString *));
		pBARData->vParamValue = (struct tagString **)calloc(pBARData->nParamNum, sizeof(struct tagString *));

		if( (pBARData->vParamName == NULL) || (pBARData->vParamValue == NULL) )
		{
			printf("Error: Affy_BAR2TXT, cannot allocate memory for loading parameter name/value pairs\n");
			exit(EXIT_FAILURE);
		}

		for(ni=0; ni<pBARData->nParamNum; ni++)
		{
			if(big_endian_fread(&nLen, INT_SIZE, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: Affy_BAR2TXT, cannot load parameter name length.\n");
				exit(EXIT_FAILURE);
			}

			if(nLen > 0)
			{
				pBARData->vParamName[ni] = CreateString(nLen);
				if(pBARData->vParamName[ni] == NULL)
				{
					printf("Error: Affy_BAR2TXT, cannot load parameter name.\n");
					exit(EXIT_FAILURE);
				}
			
				if(big_endian_fread(pBARData->vParamName[ni]->m_pString, 1, nLen, fpIn, little_endian_machine) != nLen)
				{
					printf("Error: Affy_BAR2TXT, cannot load parameter name.\n");
					exit(EXIT_FAILURE);
				}
				pBARData->vParamName[ni]->m_pString[nLen] = '\0';
			}

			if(big_endian_fread(&nLen, INT_SIZE, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: Affy_BAR2TXT, cannot load parameter value length.\n");
				exit(EXIT_FAILURE);
			}

			if(nLen > 0)
			{
				pBARData->vParamValue[ni] = CreateString(nLen);
				if(pBARData->vParamValue[ni] == NULL)
				{
					printf("Error: Affy_BAR2TXT, cannot load parameter value.\n");
					exit(EXIT_FAILURE);
				}

				if(big_endian_fread(pBARData->vParamValue[ni]->m_pString, 1, nLen, fpIn, little_endian_machine) != nLen)
				{
					printf("Error: Affy_BAR2TXT, cannot load parameter value.\n");
					exit(EXIT_FAILURE);
				}
				pBARData->vParamValue[ni]->m_pString[nLen] = '\0';
			}
		}
	}

	/* load sequence data */
	if(pBARData->nSeqNum <= 0)
	{
		/* close file */
		Affy_BARData_Destroy(&pBARData);
		fclose(fpIn);
		fclose(fpOut);
		return PROC_SUCCESS;
	}

	pBARData->vSeqData = (struct tagBARSeq **)calloc(pBARData->nSeqNum, sizeof(struct tagBARSeq *));
	if(pBARData->vSeqData == NULL)
	{
		printf("Error: Affy_BAR2TXT, cannot allocate memory for loading sequence data.\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<pBARData->nSeqNum; ni++)
	{
		/* create BARSeq object */
		pBARData->vSeqData[ni] = Affy_BARSeq_Create();
		if(pBARData->vSeqData[ni] == NULL)
		{
			printf("Error: Affy_BAR2TXT, cannot create BARSeq object.\n");
			exit(EXIT_FAILURE);
		}
		pBARData->vSeqData[ni]->nColNum = pBARData->nColNum;


		/* load sequence name */
		if(big_endian_fread(&nLen, INT_SIZE, 1, fpIn, little_endian_machine) != 1)
		{
			printf("Error: Affy_BAR2TXT, cannot load sequence name length.\n");
			exit(EXIT_FAILURE);
		}

		if(nLen > 0)
		{
			pBARData->vSeqData[ni]->pSeqName = CreateString(nLen);
			if(pBARData->vSeqData[ni]->pSeqName == NULL)
			{
				printf("Error: Affy_BAR2TXT, cannot load sequence name.\n");
				exit(EXIT_FAILURE);
			}
		
			if(big_endian_fread(pBARData->vSeqData[ni]->pSeqName->m_pString, 1, nLen, fpIn, little_endian_machine) != nLen)
			{
				printf("Error: Affy_BAR2TXT, cannot load sequence name.\n");
				exit(EXIT_FAILURE);
			}
			pBARData->vSeqData[ni]->pSeqName->m_pString[nLen] = '\0';
		}

		/* load group name */
		if(pBARData->fVersionnumber > 1.5)
		{
			if(big_endian_fread(&nLen, INT_SIZE, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: Affy_BAR2TXT, cannot load sequence group name length.\n");
				exit(EXIT_FAILURE);
			}

			if(nLen > 0)
			{
				pBARData->vSeqData[ni]->pSeqGroupName = CreateString(nLen);
				if(pBARData->vSeqData[ni]->pSeqGroupName == NULL)
				{
					printf("Error: Affy_BAR2TXT, cannot load sequence group name.\n");
					exit(EXIT_FAILURE);
				}
			
				if(big_endian_fread(pBARData->vSeqData[ni]->pSeqGroupName->m_pString, 1, nLen, fpIn, little_endian_machine) != nLen)
				{
					printf("Error: Affy_BAR2TXT, cannot load sequence group name.\n");
					exit(EXIT_FAILURE);
				}
				pBARData->vSeqData[ni]->pSeqGroupName->m_pString[nLen] = '\0';
			}
		}

		/* load sequence version */
		if(big_endian_fread(&nLen, INT_SIZE, 1, fpIn, little_endian_machine) != 1)
		{
			printf("Error: Affy_BAR2TXT, cannot load sequence version.\n");
			exit(EXIT_FAILURE);
		}

		if(nLen > 0)
		{
			pBARData->vSeqData[ni]->pSeqVersion = CreateString(nLen);
			if(pBARData->vSeqData[ni]->pSeqVersion == NULL)
			{
				printf("Error: Affy_BAR2TXT, cannot load sequence version.\n");
				exit(EXIT_FAILURE);
			}
		
			if(big_endian_fread(pBARData->vSeqData[ni]->pSeqVersion->m_pString, 1, nLen, fpIn, little_endian_machine) != nLen)
			{
				printf("Error: Affy_BAR2TXT, cannot load sequence version.\n");
				exit(EXIT_FAILURE);
			}
			pBARData->vSeqData[ni]->pSeqVersion->m_pString[nLen] = '\0';
		}

		/* load parameter name/value pairs */
		if(pBARData->fVersionnumber > 1.5)
		{
			if(big_endian_fread(&(pBARData->vSeqData[ni]->nParamNum), INT_SIZE, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: Affy_BAR2TXT, cannot load number of name/value pairs for BARseq object.\n");
				exit(EXIT_FAILURE);
			}

			if(pBARData->vSeqData[ni]->nParamNum > 0)
			{
				pBARData->vSeqData[ni]->vParamName = (struct tagString **)calloc(pBARData->vSeqData[ni]->nParamNum, sizeof(struct tagString *));
				pBARData->vSeqData[ni]->vParamValue = (struct tagString **)calloc(pBARData->vSeqData[ni]->nParamNum, sizeof(struct tagString *));

				if( (pBARData->vSeqData[ni]->vParamName == NULL) || (pBARData->vSeqData[ni]->vParamValue == NULL) )
				{
					printf("Error: Affy_BAR2TXT, cannot allocate memory for loading sequence parameter name/value pairs\n");
					exit(EXIT_FAILURE);
				}

				for(nj=0; nj<pBARData->vSeqData[ni]->nParamNum; nj++)
				{
					if(big_endian_fread(&nLen, INT_SIZE, 1, fpIn, little_endian_machine) != 1)
					{
						printf("Error: Affy_BAR2TXT, cannot load sequence parameter name length.\n");
						exit(EXIT_FAILURE);
					}

					if(nLen > 0)
					{
						pBARData->vSeqData[ni]->vParamName[nj] = CreateString(nLen);
						if(pBARData->vSeqData[ni]->vParamName[nj] == NULL)
						{
							printf("Error: Affy_BAR2TXT, cannot load sequence parameter name.\n");
							exit(EXIT_FAILURE);
						}

						if(big_endian_fread(pBARData->vSeqData[ni]->vParamName[nj]->m_pString, 1, nLen, fpIn, little_endian_machine) != nLen)
						{
							printf("Error: Affy_BAR2TXT, cannot load sequence parameter name.\n");
							exit(EXIT_FAILURE);
						}
						pBARData->vSeqData[ni]->vParamName[nj]->m_pString[nLen] = '\0';
					}

					if(big_endian_fread(&nLen, INT_SIZE, 1, fpIn, little_endian_machine) != 1)
					{
						printf("Error: Affy_BAR2TXT, cannot load sequence parameter value length.\n");
						exit(EXIT_FAILURE);
					}

					if(nLen > 0)
					{					
						pBARData->vSeqData[ni]->vParamValue[nj] = CreateString(nLen);
						if(pBARData->vSeqData[ni]->vParamValue[nj] == NULL)
						{
							printf("Error: Affy_BAR2TXT, cannot load sequence parameter value.\n");
							exit(EXIT_FAILURE);
						}

						if(big_endian_fread(pBARData->vSeqData[ni]->vParamValue[nj]->m_pString, 1, nLen, fpIn, little_endian_machine) != nLen)
						{
							printf("Error: Affy_BAR2TXT, cannot load parameter value.\n");
							exit(EXIT_FAILURE);
						}
						pBARData->vSeqData[ni]->vParamValue[nj]->m_pString[nLen] = '\0';
					}
				}
			}
		}

		/* load data points */
		if(big_endian_fread(&(pBARData->vSeqData[ni]->nDataNum), INT_SIZE, 1, fpIn, little_endian_machine) != 1)
		{
			printf("Error: Affy_BAR2TXT, cannot load number of data points for a sequence.\n");
			exit(EXIT_FAILURE);
		}

		if(pBARData->vSeqData[ni]->pSeqName == NULL)
		{
			printf("Error: TileMapv2_ExportAffyIntensity_SingleTXT, cannot find sequence name.\n");
			exit(EXIT_FAILURE);
		}
		if(pBARData->vSeqData[ni]->pSeqVersion != NULL)
		{
			sprintf(strSeqAlias, "%s:%s", pBARData->vSeqData[ni]->pSeqVersion->m_pString,
				pBARData->vSeqData[ni]->pSeqName->m_pString);
		}
		else
		{
			sprintf(strSeqAlias, "%s", pBARData->vSeqData[ni]->pSeqName->m_pString);
		}

		if(pBARData->fVersionnumber > 1.5)
		{
			if(pBARData->vSeqData[ni]->pSeqGroupName != NULL)
			{
				if(pBARData->vSeqData[ni]->pSeqVersion != NULL)
				{
					sprintf(strSeqAlias, "%s:%s:%s", pBARData->vSeqData[ni]->pSeqGroupName->m_pString,
						pBARData->vSeqData[ni]->pSeqVersion->m_pString,
						pBARData->vSeqData[ni]->pSeqName->m_pString);
				}
				else
				{
					sprintf(strSeqAlias, "%s:%s", pBARData->vSeqData[ni]->pSeqGroupName->m_pString,
						pBARData->vSeqData[ni]->pSeqName->m_pString);
				}
			}	
		}

		if(pBARData->vSeqData[ni]->nColNum > 0)
		{
			pBARData->vSeqData[ni]->vData = (struct DOUBLEMATRIX **)calloc(pBARData->vSeqData[ni]->nColNum, sizeof(struct DOUBLEMATRIX *));
			if(pBARData->vSeqData[ni]->vData == NULL)
			{
				printf("Error: Affy_BAR2TXT, cannot allocate memory for loading sequence data.\n");
				exit(EXIT_FAILURE);
			}

			if(pBARData->vSeqData[ni]->nDataNum > 0)
			{
				for(nk=0; nk<pBARData->vSeqData[ni]->nDataNum; nk++)
				{
					fprintf(fpOut, "%s", strSeqAlias);
					for(nj=0; nj<pBARData->nColNum; nj++)
					{
						switch(pBARData->pFieldType->pMatElement[nj])
						{
							case 0: big_endian_fread(&dV, DOUBLE_SIZE, 1, fpIn, little_endian_machine);
								fprintf(fpOut, "\t%lf", dV);
								break;
							case 1: big_endian_fread(&fV, FLOAT_SIZE, 1, fpIn, little_endian_machine);
								fprintf(fpOut, "\t%f", fV);
								break;
							case 2: big_endian_fread(&nV, INT_SIZE, 1, fpIn, little_endian_machine);
								fprintf(fpOut, "\t%d", nV);
								break;
							case 3: big_endian_fread(&sV, SHORT_SIZE, 1, fpIn, little_endian_machine);
								fprintf(fpOut, "\t%d", sV);
								break;
							case 4: big_endian_fread(&cV, CHAR_SIZE, 1, fpIn, little_endian_machine);
								fprintf(fpOut, "\t%d", cV);
								break;
							case 5: big_endian_fread(&unV, DWORD_SIZE, 1, fpIn, little_endian_machine);
								fprintf(fpOut, "\t%d", unV);
								break;
							case 6: big_endian_fread(&usV, USHORT_SIZE, 1, fpIn, little_endian_machine);
								fprintf(fpOut, "\t%lf", usV);
								break;
							case 7: big_endian_fread(&ucV, UCHAR_SIZE, 1, fpIn, little_endian_machine);
								fprintf(fpOut, "\t%d", ucV);
								break;
						}
					}
					fprintf(fpOut, "\n");
				}
			}
		}
	}

	/* close file */
	Affy_BARData_Destroy(&pBARData);
	fclose(fpIn);
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_BAR2TXT_Fast()                                                    */
/*  Convert *.bar file to *.txt file.                                      */
/* ----------------------------------------------------------------------- */ 
int Affy_BAR2TXT_Fast(char strBARFile[], char strTXTFile[])
{
	/* define */
	struct tagBARData *pBARData = NULL;

	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));
	size_t unitsize;
	char *tempmem;
	char *tempptr;
	
	FILE *fpIn,*fpOut;
	int ni,nj,nk;
	int nLen;

	double dV;
	float fV;
	int nV;
	short sV;
	char cV;
	unsigned int unV;
	unsigned short usV;
	unsigned char ucV;
	char strSeqAlias[MED_LINE_LENGTH];

	/* init */
	fpIn = NULL;
	fpIn = fopen(strBARFile, "rb");
	if(fpIn == NULL)
	{
		printf("Error: Affy_BAR2TXT, cannot open input *.bar file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strTXTFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: Affy_BAR2TXT, cannot open output *.txt file!\n");
		exit(EXIT_FAILURE);
	}

	pBARData = Affy_BARData_Create();
	if(pBARData == NULL)
	{
		printf("Warning: Affy_BAR2TXT, BARData object was not created!\n");
		exit(EXIT_FAILURE);
	}

	/* load magic number */
	if(big_endian_fread(pBARData->strMagicnumber, 1, 8, fpIn, little_endian_machine) != 8)
	{
		printf("Error: Affy_BAR2TXT, cannot load magic number.\n");
		exit(EXIT_FAILURE);
	}
	pBARData->strMagicnumber[8] = '\0';
	
	if(strcmp(pBARData->strMagicnumber, "barr\r\n\032\n") != 0)
	{
		printf("Error: Affy_BAR2TXT, cannot load magic number correctly.\n");
		exit(EXIT_FAILURE);
	}

	/* load version number */
	if(big_endian_fread(&(pBARData->fVersionnumber), FLOAT_SIZE, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: Affy_BAR2TXT, cannot load version number.\n");
		exit(EXIT_FAILURE);
	}

	/* load sequence number */
	if(big_endian_fread(&(pBARData->nSeqNum), INT_SIZE, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: Affy_BAR2TXT, cannot load sequence number.\n");
		exit(EXIT_FAILURE);
	}

	/* load column number */
	if(big_endian_fread(&(pBARData->nColNum), INT_SIZE, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: Affy_BAR2TXT, cannot load column number.\n");
		exit(EXIT_FAILURE);
	}

	/* load column type */
	if(pBARData->nColNum>0)
	{
		pBARData->pFieldType = CreateIntMatrix(1, pBARData->nColNum);
		if(pBARData->pFieldType == NULL)
		{
			printf("Error: Affy_BAR2TXT, cannot allocate memory for field type.\n");
			exit(EXIT_FAILURE);
		}
	
		if((int)big_endian_fread(pBARData->pFieldType->pMatElement, INT_SIZE, pBARData->nColNum, fpIn, little_endian_machine) != pBARData->nColNum)
		{
			printf("Error: Affy_BAR2TXT, cannot load field type correctly.\n");
			exit(EXIT_FAILURE);
		}
	}

	/* load parameter name/value pairs */
	if(big_endian_fread(&(pBARData->nParamNum), INT_SIZE, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: Affy_BAR2TXT, cannot load number of name/value pairs.\n");
		exit(EXIT_FAILURE);
	}

	if(pBARData->nParamNum > 0)
	{
		pBARData->vParamName = (struct tagString **)calloc(pBARData->nParamNum, sizeof(struct tagString *));
		pBARData->vParamValue = (struct tagString **)calloc(pBARData->nParamNum, sizeof(struct tagString *));

		if( (pBARData->vParamName == NULL) || (pBARData->vParamValue == NULL) )
		{
			printf("Error: Affy_BAR2TXT, cannot allocate memory for loading parameter name/value pairs\n");
			exit(EXIT_FAILURE);
		}

		for(ni=0; ni<pBARData->nParamNum; ni++)
		{
			if(big_endian_fread(&nLen, INT_SIZE, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: Affy_BAR2TXT, cannot load parameter name length.\n");
				exit(EXIT_FAILURE);
			}

			if(nLen > 0)
			{
				pBARData->vParamName[ni] = CreateString(nLen);
				if(pBARData->vParamName[ni] == NULL)
				{
					printf("Error: Affy_BAR2TXT, cannot load parameter name.\n");
					exit(EXIT_FAILURE);
				}
			
				if(big_endian_fread(pBARData->vParamName[ni]->m_pString, 1, nLen, fpIn, little_endian_machine) != nLen)
				{
					printf("Error: Affy_BAR2TXT, cannot load parameter name.\n");
					exit(EXIT_FAILURE);
				}
				pBARData->vParamName[ni]->m_pString[nLen] = '\0';
			}

			if(big_endian_fread(&nLen, INT_SIZE, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: Affy_BAR2TXT, cannot load parameter value length.\n");
				exit(EXIT_FAILURE);
			}

			if(nLen > 0)
			{
				pBARData->vParamValue[ni] = CreateString(nLen);
				if(pBARData->vParamValue[ni] == NULL)
				{
					printf("Error: Affy_BAR2TXT, cannot load parameter value.\n");
					exit(EXIT_FAILURE);
				}

				if(big_endian_fread(pBARData->vParamValue[ni]->m_pString, 1, nLen, fpIn, little_endian_machine) != nLen)
				{
					printf("Error: Affy_BAR2TXT, cannot load parameter value.\n");
					exit(EXIT_FAILURE);
				}
				pBARData->vParamValue[ni]->m_pString[nLen] = '\0';
			}
		}
	}

	/* load sequence data */
	if(pBARData->nSeqNum <= 0)
	{
		/* close file */
		Affy_BARData_Destroy(&pBARData);
		fclose(fpIn);
		fclose(fpOut);
		return PROC_SUCCESS;
	}

	unitsize = 0;
	for(nj=0; nj<pBARData->nColNum; nj++)
	{
		switch(pBARData->pFieldType->pMatElement[nj])
		{
			case 0: unitsize += sizeof(double);
				break;
			case 1: unitsize += sizeof(float);
				break;
			case 2: unitsize += 4;
				break;
			case 3: unitsize += 2;
				break;
			case 4: unitsize += 1;
				break;
			case 5: unitsize += 4;
				break;
			case 6: unitsize += 2;
				break;
			case 7: unitsize += 1;
				break;
		}
	}

	pBARData->vSeqData = (struct tagBARSeq **)calloc(pBARData->nSeqNum, sizeof(struct tagBARSeq *));
	if(pBARData->vSeqData == NULL)
	{
		printf("Error: Affy_BAR2TXT, cannot allocate memory for loading sequence data.\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<pBARData->nSeqNum; ni++)
	{
		/* create BARSeq object */
		pBARData->vSeqData[ni] = Affy_BARSeq_Create();
		if(pBARData->vSeqData[ni] == NULL)
		{
			printf("Error: Affy_BAR2TXT, cannot create BARSeq object.\n");
			exit(EXIT_FAILURE);
		}
		pBARData->vSeqData[ni]->nColNum = pBARData->nColNum;


		/* load sequence name */
		if(big_endian_fread(&nLen, INT_SIZE, 1, fpIn, little_endian_machine) != 1)
		{
			printf("Error: Affy_BAR2TXT, cannot load sequence name length.\n");
			exit(EXIT_FAILURE);
		}

		if(nLen > 0)
		{
			pBARData->vSeqData[ni]->pSeqName = CreateString(nLen);
			if(pBARData->vSeqData[ni]->pSeqName == NULL)
			{
				printf("Error: Affy_BAR2TXT, cannot load sequence name.\n");
				exit(EXIT_FAILURE);
			}
		
			if(big_endian_fread(pBARData->vSeqData[ni]->pSeqName->m_pString, 1, nLen, fpIn, little_endian_machine) != nLen)
			{
				printf("Error: Affy_BAR2TXT, cannot load sequence name.\n");
				exit(EXIT_FAILURE);
			}
			pBARData->vSeqData[ni]->pSeqName->m_pString[nLen] = '\0';
		}

		/* load group name */
		if(pBARData->fVersionnumber > 1.5)
		{
			if(big_endian_fread(&nLen, INT_SIZE, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: Affy_BAR2TXT, cannot load sequence group name length.\n");
				exit(EXIT_FAILURE);
			}

			if(nLen > 0)
			{
				pBARData->vSeqData[ni]->pSeqGroupName = CreateString(nLen);
				if(pBARData->vSeqData[ni]->pSeqGroupName == NULL)
				{
					printf("Error: Affy_BAR2TXT, cannot load sequence group name.\n");
					exit(EXIT_FAILURE);
				}
			
				if(big_endian_fread(pBARData->vSeqData[ni]->pSeqGroupName->m_pString, 1, nLen, fpIn, little_endian_machine) != nLen)
				{
					printf("Error: Affy_BAR2TXT, cannot load sequence group name.\n");
					exit(EXIT_FAILURE);
				}
				pBARData->vSeqData[ni]->pSeqGroupName->m_pString[nLen] = '\0';
			}
		}

		/* load sequence version */
		if(big_endian_fread(&nLen, INT_SIZE, 1, fpIn, little_endian_machine) != 1)
		{
			printf("Error: Affy_BAR2TXT, cannot load sequence version.\n");
			exit(EXIT_FAILURE);
		}

		if(nLen > 0)
		{
			pBARData->vSeqData[ni]->pSeqVersion = CreateString(nLen);
			if(pBARData->vSeqData[ni]->pSeqVersion == NULL)
			{
				printf("Error: Affy_BAR2TXT, cannot load sequence version.\n");
				exit(EXIT_FAILURE);
			}
		
			if(big_endian_fread(pBARData->vSeqData[ni]->pSeqVersion->m_pString, 1, nLen, fpIn, little_endian_machine) != nLen)
			{
				printf("Error: Affy_BAR2TXT, cannot load sequence version.\n");
				exit(EXIT_FAILURE);
			}
			pBARData->vSeqData[ni]->pSeqVersion->m_pString[nLen] = '\0';
		}

		/* load parameter name/value pairs */
		if(pBARData->fVersionnumber > 1.5)
		{
			if(big_endian_fread(&(pBARData->vSeqData[ni]->nParamNum), INT_SIZE, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: Affy_BAR2TXT, cannot load number of name/value pairs for BARseq object.\n");
				exit(EXIT_FAILURE);
			}

			if(pBARData->vSeqData[ni]->nParamNum > 0)
			{
				pBARData->vSeqData[ni]->vParamName = (struct tagString **)calloc(pBARData->vSeqData[ni]->nParamNum, sizeof(struct tagString *));
				pBARData->vSeqData[ni]->vParamValue = (struct tagString **)calloc(pBARData->vSeqData[ni]->nParamNum, sizeof(struct tagString *));

				if( (pBARData->vSeqData[ni]->vParamName == NULL) || (pBARData->vSeqData[ni]->vParamValue == NULL) )
				{
					printf("Error: Affy_BAR2TXT, cannot allocate memory for loading sequence parameter name/value pairs\n");
					exit(EXIT_FAILURE);
				}

				for(nj=0; nj<pBARData->vSeqData[ni]->nParamNum; nj++)
				{
					if(big_endian_fread(&nLen, INT_SIZE, 1, fpIn, little_endian_machine) != 1)
					{
						printf("Error: Affy_BAR2TXT, cannot load sequence parameter name length.\n");
						exit(EXIT_FAILURE);
					}

					if(nLen > 0)
					{
						pBARData->vSeqData[ni]->vParamName[nj] = CreateString(nLen);
						if(pBARData->vSeqData[ni]->vParamName[nj] == NULL)
						{
							printf("Error: Affy_BAR2TXT, cannot load sequence parameter name.\n");
							exit(EXIT_FAILURE);
						}

						if(big_endian_fread(pBARData->vSeqData[ni]->vParamName[nj]->m_pString, 1, nLen, fpIn, little_endian_machine) != nLen)
						{
							printf("Error: Affy_BAR2TXT, cannot load sequence parameter name.\n");
							exit(EXIT_FAILURE);
						}
						pBARData->vSeqData[ni]->vParamName[nj]->m_pString[nLen] = '\0';
					}

					if(big_endian_fread(&nLen, INT_SIZE, 1, fpIn, little_endian_machine) != 1)
					{
						printf("Error: Affy_BAR2TXT, cannot load sequence parameter value length.\n");
						exit(EXIT_FAILURE);
					}

					if(nLen > 0)
					{					
						pBARData->vSeqData[ni]->vParamValue[nj] = CreateString(nLen);
						if(pBARData->vSeqData[ni]->vParamValue[nj] == NULL)
						{
							printf("Error: Affy_BAR2TXT, cannot load sequence parameter value.\n");
							exit(EXIT_FAILURE);
						}

						if(big_endian_fread(pBARData->vSeqData[ni]->vParamValue[nj]->m_pString, 1, nLen, fpIn, little_endian_machine) != nLen)
						{
							printf("Error: Affy_BAR2TXT, cannot load parameter value.\n");
							exit(EXIT_FAILURE);
						}
						pBARData->vSeqData[ni]->vParamValue[nj]->m_pString[nLen] = '\0';
					}
				}
			}
		}

		/* load data points */
		if(big_endian_fread(&(pBARData->vSeqData[ni]->nDataNum), INT_SIZE, 1, fpIn, little_endian_machine) != 1)
		{
			printf("Error: Affy_BAR2TXT, cannot load number of data points for a sequence.\n");
			exit(EXIT_FAILURE);
		}

		if(pBARData->vSeqData[ni]->pSeqName == NULL)
		{
			printf("Error: TileMapv2_ExportAffyIntensity_SingleTXT, cannot find sequence name.\n");
			exit(EXIT_FAILURE);
		}
		if(pBARData->vSeqData[ni]->pSeqVersion != NULL)
		{
			sprintf(strSeqAlias, "%s:%s", pBARData->vSeqData[ni]->pSeqVersion->m_pString,
				pBARData->vSeqData[ni]->pSeqName->m_pString);
		}
		else
		{
			sprintf(strSeqAlias, "%s", pBARData->vSeqData[ni]->pSeqName->m_pString);
		}

		if(pBARData->fVersionnumber > 1.5)
		{
			if(pBARData->vSeqData[ni]->pSeqGroupName != NULL)
			{
				if(pBARData->vSeqData[ni]->pSeqVersion != NULL)
				{
					sprintf(strSeqAlias, "%s:%s:%s", pBARData->vSeqData[ni]->pSeqGroupName->m_pString,
						pBARData->vSeqData[ni]->pSeqVersion->m_pString,
						pBARData->vSeqData[ni]->pSeqName->m_pString);
				}
				else
				{
					sprintf(strSeqAlias, "%s:%s", pBARData->vSeqData[ni]->pSeqGroupName->m_pString,
						pBARData->vSeqData[ni]->pSeqName->m_pString);
				}
			}	
		}

		if(pBARData->vSeqData[ni]->nColNum > 0)
		{
			pBARData->vSeqData[ni]->vData = (struct DOUBLEMATRIX **)calloc(pBARData->vSeqData[ni]->nColNum, sizeof(struct DOUBLEMATRIX *));
			if(pBARData->vSeqData[ni]->vData == NULL)
			{
				printf("Error: Affy_BAR2TXT, cannot allocate memory for loading sequence data.\n");
				exit(EXIT_FAILURE);
			}

			if(pBARData->vSeqData[ni]->nDataNum > 0)
			{
				tempmem = NULL;
				tempmem = (char *)calloc(pBARData->vSeqData[ni]->nDataNum, unitsize);
				if(tempmem == NULL)
				{
					printf("Error: Affy_LoadBAR, cannot allocate memory for loading sequence data.\n");
					exit(EXIT_FAILURE);
				}

				if( (int)fread(tempmem, unitsize, pBARData->vSeqData[ni]->nDataNum, fpIn) != pBARData->vSeqData[ni]->nDataNum)
				{
					printf("Error: Affy_LoadBAR, cannot load sequence data correctly.\n");
					exit(EXIT_FAILURE);
				}

				tempptr = tempmem;
				for(nk=0; nk<pBARData->vSeqData[ni]->nDataNum; nk++)
				{
					fprintf(fpOut, "%s", strSeqAlias);
					for(nj=0; nj<pBARData->nColNum; nj++)
					{
						switch(pBARData->pFieldType->pMatElement[nj])
						{
							case 0: 
								if(little_endian_machine == 1)
								{
									reverse_buf(tempptr, sizeof(double));
								}
								memcpy(&dV, tempptr, sizeof(double));
								fprintf(fpOut, "\t%lf", dV);
								tempptr += sizeof(double);
								break;
							case 1: 
								if(little_endian_machine == 1)
								{
									reverse_buf(tempptr, sizeof(float));
								}
								memcpy(&fV, tempptr, sizeof(float));
								fprintf(fpOut, "\t%f", fV);
								tempptr += sizeof(float);
								break;
							case 2: 
								if(little_endian_machine == 1)
								{
									reverse_buf(tempptr, 4);
								}
								memcpy(&nV, tempptr, 4);
								fprintf(fpOut, "\t%d", nV);
								tempptr += 4;
								break;
							case 3: 
								if(little_endian_machine == 1)
								{
									reverse_buf(tempptr, 2);
								}
								memcpy(&sV, tempptr, 2);
								fprintf(fpOut, "\t%d", sV);
								tempptr += 2;
								break;
							case 4: 
								if(little_endian_machine == 1)
								{
									reverse_buf(tempptr, 1);
								}
								memcpy(&cV, tempptr, 1);
								fprintf(fpOut, "\t%d", cV);
								tempptr += 1;
								break;
							case 5: 
								if(little_endian_machine == 1)
								{
									reverse_buf(tempptr, 4);
								}
								memcpy(&unV, tempptr, 4);
								fprintf(fpOut, "\t%d", unV);
								tempptr += 4;
								break;
							case 6: 
								if(little_endian_machine == 1)
								{
									reverse_buf(tempptr, 2);
								}
								memcpy(&usV, tempptr, 2);
								fprintf(fpOut, "\t%lf", usV);
								tempptr += 2;
								break;
							case 7: 
								if(little_endian_machine == 1)
								{
									reverse_buf(tempptr, 1);
								}
								memcpy(&ucV, tempptr, 1);
								fprintf(fpOut, "\t%d", ucV);
								tempptr += 1;
								break;
						}
					}
					fprintf(fpOut, "\n");
				}

				free(tempmem);
			}
		}
	}

	/* close file */
	Affy_BARData_Destroy(&pBARData);
	fclose(fpIn);
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_BAR2WIG_Main()                                                    */
/*  Convert BAR file to WIG file.                                          */
/* ----------------------------------------------------------------------- */ 
int Affy_BAR2WIG_Main(char strInputFile[], char strOutputFile[])
{
	/* define */
	struct tagBARData *pData = NULL;
	int nSeqId = 0;
	int nStartPos = 0;
	int nEndPos = 0;
	int ni;
	FILE *fpOut = NULL;
	char strFileName[MED_LINE_LENGTH];


	/* load bar file */
	GetFileName(strOutputFile, strFileName);

	pData = Affy_LoadBAR_Fast(strInputFile);
	if(pData == NULL)
	{
		printf("Error: Affy_BAR2WIG_Main, cannot open the BAR file!\n");
		exit(EXIT_FAILURE);
	}

	/* open output file */
	fpOut = fopen(strOutputFile, "wt");
	if(fpOut == NULL)
	{
		printf("Error: Affy_BAR2WIG_Main, cannot open the output WIG file!\n");
		exit(EXIT_FAILURE);
	}

	/* write */
	fprintf(fpOut, "track type=wiggle_0 name=%s description=%s graphType=bar\n", strFileName, strFileName);
	for(nSeqId = 0; nSeqId<pData->nSeqNum; nSeqId++) 
	{
		fprintf(fpOut, "\nvariableStep chrom=%s\n", pData->vSeqData[nSeqId]->pSeqName->m_pString);  
		nStartPos = 0;
		nEndPos = pData->vSeqData[nSeqId]->nDataNum - 1;
		for(ni=nStartPos; ni<=nEndPos; ni++)
		{
			fprintf(fpOut, "%d\t%f\n", (int)(pData->vSeqData[nSeqId]->vData[0]->pMatElement[ni]),
				pData->vSeqData[nSeqId]->vData[1]->pMatElement[ni]);
		}
	}

	/* close file */
	fclose(fpOut);
	Affy_BARData_Destroy(&pData);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_BARSeq_Clone()                                                    */
/*  clone a BARSeq object to a new BARData object.                         */
/* ----------------------------------------------------------------------- */ 
struct tagBARSeq *Affy_BARSeq_Clone(struct tagBARSeq *pBARSeq)
{
	/* define */
	struct tagBARSeq *pNewBARSeq = NULL;
	int ni, nj;

	/* init */
	if(pBARSeq == NULL)
		return NULL;

	pNewBARSeq = Affy_BARSeq_Create();
	if(pNewBARSeq == NULL)
	{
		printf("Error: Affy_BARSeq_Clone, cannot create a new BARSeq object.\n");
		exit(EXIT_FAILURE);
	}
	pNewBARSeq->nColNum = pBARSeq->nColNum;

	if(pBARSeq->pSeqName != NULL)
	{
		StringAddTail(&(pNewBARSeq->pSeqName), pBARSeq->pSeqName->m_pString);
	}

	if(pBARSeq->pSeqGroupName != NULL)
	{
		StringAddTail(&(pNewBARSeq->pSeqGroupName), pBARSeq->pSeqGroupName->m_pString);
	}

	if(pBARSeq->pSeqVersion != NULL)
	{
		StringAddTail(&(pNewBARSeq->pSeqVersion), pBARSeq->pSeqVersion->m_pString);
	}
	
	pNewBARSeq->nParamNum = pBARSeq->nParamNum;
	if(pNewBARSeq->nParamNum > 0)
	{
		pNewBARSeq->vParamName = (struct tagString **)calloc(pNewBARSeq->nParamNum, sizeof(struct tagString *));
		pNewBARSeq->vParamValue = (struct tagString **)calloc(pNewBARSeq->nParamNum, sizeof(struct tagString *));

		if( (pNewBARSeq->vParamName == NULL) || (pNewBARSeq->vParamValue == NULL) )
		{
			printf("Error: Affy_BARSeq_Clone, cannot allocate memory for loading sequence parameter name/value pairs\n");
			exit(EXIT_FAILURE);
		}

		for(ni=0; ni<pNewBARSeq->nParamNum; ni++)
		{
			if(pBARSeq->vParamName[ni] != NULL)
			{
				StringAddTail(pNewBARSeq->vParamName+ni, pBARSeq->vParamName[ni]->m_pString);
			}

			if(pBARSeq->vParamValue[ni] != NULL)
			{
				StringAddTail(pNewBARSeq->vParamValue+ni, pBARSeq->vParamValue[ni]->m_pString);
			}
		}
	}

	pNewBARSeq->nDataNum = pBARSeq->nDataNum;

	if(pNewBARSeq->nColNum > 0)
	{
		pNewBARSeq->vData = (struct DOUBLEMATRIX **)calloc(pNewBARSeq->nColNum, sizeof(struct DOUBLEMATRIX *));
		if(pNewBARSeq->vData == NULL)
		{
			printf("Error: Affy_BARSeq_Clone, cannot allocate memory for loading sequence data.\n");
			exit(EXIT_FAILURE);
		}

		if(pNewBARSeq->nDataNum > 0)
		{
			for(nj=0; nj<pNewBARSeq->nColNum; nj++)
			{
				pNewBARSeq->vData[nj] = DMCLONE(pBARSeq->vData[nj]);
			}
		}
	}

	/* return */
	return pNewBARSeq;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_BARData_Clone()                                                   */
/*  clone a BARData object to a new BARData object.                        */
/* ----------------------------------------------------------------------- */ 
struct tagBARData *Affy_BARData_Clone(struct tagBARData *pBARData)
{
	/* define */
	struct tagBARData *pNewBARData = NULL;
	int ni;

	/* init */
	if(pBARData == NULL)
		return NULL;

	/* copy */
	pNewBARData = Affy_BARData_Create();
	if(pNewBARData == NULL)
	{
		printf("Error: Affy_BARData_Clone, cannot create a new BARData object!\n");
		exit(EXIT_FAILURE);
	}

	/* header info */
	strcpy(pNewBARData->strMagicnumber, pBARData->strMagicnumber);
	pNewBARData->fVersionnumber = pBARData->fVersionnumber;
	pNewBARData->nSeqNum = pBARData->nSeqNum;
	pNewBARData->nColNum = pBARData->nColNum;
	pNewBARData->pFieldType = IMCLONE(pBARData->pFieldType);
	pNewBARData->nParamNum = pBARData->nParamNum;

	/* load parameter name/value pairs */
	if(pNewBARData->nParamNum > 0)
	{
		pNewBARData->vParamName = (struct tagString **)calloc(pNewBARData->nParamNum, sizeof(struct tagString *));
		pNewBARData->vParamValue = (struct tagString **)calloc(pNewBARData->nParamNum, sizeof(struct tagString *));

		if( (pNewBARData->vParamName == NULL) || (pNewBARData->vParamValue == NULL) )
		{
			printf("Error: Affy_BARData_Clone, cannot allocate memory for loading parameter name/value pairs\n");
			exit(EXIT_FAILURE);
		}

		for(ni=0; ni<pNewBARData->nParamNum; ni++)
		{
			if(pBARData->vParamName[ni] != NULL)
			{
				StringAddTail(pNewBARData->vParamName+ni, pBARData->vParamName[ni]->m_pString);
			}
			if(pBARData->vParamValue[ni] != NULL)
			{
				StringAddTail(pNewBARData->vParamValue+ni, pBARData->vParamValue[ni]->m_pString);
			}
		}
	}

	/* load sequence data */
	if(pNewBARData->nSeqNum <= 0)
		return pNewBARData;

	pNewBARData->vSeqData = (struct tagBARSeq **)calloc(pNewBARData->nSeqNum, sizeof(struct tagBARSeq *));
	if(pNewBARData->vSeqData == NULL)
	{
		printf("Error: Affy_BARData_Clone, cannot allocate memory for loading sequence data.\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<pNewBARData->nSeqNum; ni++)
	{
		/* create BARSeq object */
		pNewBARData->vSeqData[ni] = Affy_BARSeq_Clone(pBARData->vSeqData[ni]);		
	}

	/* return */
	return pNewBARData;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_BARData_AddTS()                                                   */
/*  Add a number with given columns of a BARData object.                   */
/* ----------------------------------------------------------------------- */ 
int Affy_BARData_AddTS(struct tagBARData *pBARData, double dLambda,
		struct INTMATRIX *pColIndicator)
{
	/* define */
	int ni,nj;

	/* init */
	if(pBARData == NULL)
		return PROC_SUCCESS;

	if(pColIndicator == NULL)
	{
		printf("Warning: no column indicators for carrying out BARData operations!\n");
		return PROC_SUCCESS;
	}

	if(pColIndicator->nWidth != pBARData->nColNum)
	{
		printf("Error: column indicator's dimension does not match BARData columns!\n");
		exit(EXIT_FAILURE);
	}

	/* operation */
	if(pBARData->nSeqNum > 0)
	{
		if(pBARData->vSeqData == NULL)
		{
			printf("Error: empty BARSeq data!\n");
			exit(EXIT_FAILURE);
		}

		for(ni=0; ni<pBARData->nSeqNum; ni++)
		{
			if(pBARData->vSeqData[ni] != NULL)
			{
				if(pBARData->vSeqData[ni]->nColNum > 0)
				{
					if(pBARData->vSeqData[ni]->vData == NULL)
					{
						printf("Error: empty BARSeq data!\n");
						exit(EXIT_FAILURE);
					}

					for(nj=0; nj<pBARData->vSeqData[ni]->nColNum; nj++)
					{
						if(pColIndicator->pMatElement[nj] == 1)
						{
							DMPADDTS(pBARData->vSeqData[ni]->vData[nj], dLambda);
						}
					}
				}
			}
		}
	}


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_BARData_SubTS()                                                   */
/*  subtract a number by given columns of a BARData object, i.e.           */
/*  dLambda-pBARData[nCol]                                                 */
/* ----------------------------------------------------------------------- */ 
int Affy_BARData_SubTS(struct tagBARData *pBARData, double dLambda,
		struct INTMATRIX *pColIndicator)
{
	/* define */
	int ni,nj;

	/* init */
	if(pBARData == NULL)
		return PROC_SUCCESS;

	if(pColIndicator == NULL)
	{
		printf("Warning: no column indicators for carrying out BARData operations!\n");
		return PROC_SUCCESS;
	}

	if(pColIndicator->nWidth != pBARData->nColNum)
	{
		printf("Error: column indicator's dimension does not match BARData columns!\n");
		exit(EXIT_FAILURE);
	}

	/* operation */
	if(pBARData->nSeqNum > 0)
	{
		if(pBARData->vSeqData == NULL)
		{
			printf("Error: empty BARSeq data!\n");
			exit(EXIT_FAILURE);
		}

		for(ni=0; ni<pBARData->nSeqNum; ni++)
		{
			if(pBARData->vSeqData[ni] != NULL)
			{
				if(pBARData->vSeqData[ni]->nColNum > 0)
				{
					if(pBARData->vSeqData[ni]->vData == NULL)
					{
						printf("Error: empty BARSeq data!\n");
						exit(EXIT_FAILURE);
					}

					for(nj=0; nj<pBARData->vSeqData[ni]->nColNum; nj++)
					{
						if(pColIndicator->pMatElement[nj] == 1)
						{
							DMPSUBTS(dLambda, pBARData->vSeqData[ni]->vData[nj]);
						}
					}
				}
			}
		}
	}


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_BARData_MultTS()                                                  */
/*  multiply a number with given columns of a BARData object.              */
/* ----------------------------------------------------------------------- */ 
int Affy_BARData_MultTS(struct tagBARData *pBARData, double dLambda,
		struct INTMATRIX *pColIndicator)
{
	/* define */
	int ni,nj;

	/* init */
	if(pBARData == NULL)
		return PROC_SUCCESS;

	if(pColIndicator == NULL)
	{
		printf("Warning: no column indicators for carrying out BARData operations!\n");
		return PROC_SUCCESS;
	}

	if(pColIndicator->nWidth != pBARData->nColNum)
	{
		printf("Error: column indicator's dimension does not match BARData columns!\n");
		exit(EXIT_FAILURE);
	}

	/* operation */
	if(pBARData->nSeqNum > 0)
	{
		if(pBARData->vSeqData == NULL)
		{
			printf("Error: empty BARSeq data!\n");
			exit(EXIT_FAILURE);
		}

		for(ni=0; ni<pBARData->nSeqNum; ni++)
		{
			if(pBARData->vSeqData[ni] != NULL)
			{
				if(pBARData->vSeqData[ni]->nColNum > 0)
				{
					if(pBARData->vSeqData[ni]->vData == NULL)
					{
						printf("Error: empty BARSeq data!\n");
						exit(EXIT_FAILURE);
					}

					for(nj=0; nj<pBARData->vSeqData[ni]->nColNum; nj++)
					{
						if(pColIndicator->pMatElement[nj] == 1)
						{
							DMPMULTS(dLambda, pBARData->vSeqData[ni]->vData[nj]);
						}
					}
				}
			}
		}
	}


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_BARData_RecTS()                                                   */
/*  divide a number by given columns of a BARData object.                  */
/* ----------------------------------------------------------------------- */ 
int Affy_BARData_RecTS(struct tagBARData *pBARData, double dLambda,
		struct INTMATRIX *pColIndicator)
{
	/* define */
	int ni,nj;

	/* init */
	if(pBARData == NULL)
		return PROC_SUCCESS;

	if(pColIndicator == NULL)
	{
		printf("Warning: no column indicators for carrying out BARData operations!\n");
		return PROC_SUCCESS;
	}

	if(pColIndicator->nWidth != pBARData->nColNum)
	{
		printf("Error: column indicator's dimension does not match BARData columns!\n");
		exit(EXIT_FAILURE);
	}

	/* operation */
	if(pBARData->nSeqNum > 0)
	{
		if(pBARData->vSeqData == NULL)
		{
			printf("Error: empty BARSeq data!\n");
			exit(EXIT_FAILURE);
		}

		for(ni=0; ni<pBARData->nSeqNum; ni++)
		{
			if(pBARData->vSeqData[ni] != NULL)
			{
				if(pBARData->vSeqData[ni]->nColNum > 0)
				{
					if(pBARData->vSeqData[ni]->vData == NULL)
					{
						printf("Error: empty BARSeq data!\n");
						exit(EXIT_FAILURE);
					}

					for(nj=0; nj<pBARData->vSeqData[ni]->nColNum; nj++)
					{
						if(pColIndicator->pMatElement[nj] == 1)
						{
							DMPRECTS(dLambda, pBARData->vSeqData[ni]->vData[nj]);
						}
					}
				}
			}
		}
	}


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_BARData_LogTS()                                                   */
/*  take logarithm of given columns of a BARData object.                   */
/* ----------------------------------------------------------------------- */ 
int Affy_BARData_LogTS(struct tagBARData *pBARData, double dBase, 
		struct INTMATRIX *pColIndicator)
{
	/* define */
	int ni,nj;

	/* init */
	if(pBARData == NULL)
		return PROC_SUCCESS;

	if(pColIndicator == NULL)
	{
		printf("Warning: no column indicators for carrying out BARData operations!\n");
		return PROC_SUCCESS;
	}

	if(pColIndicator->nWidth != pBARData->nColNum)
	{
		printf("Error: column indicator's dimension does not match BARData columns!\n");
		exit(EXIT_FAILURE);
	}

	/* operation */
	if(pBARData->nSeqNum > 0)
	{
		if(pBARData->vSeqData == NULL)
		{
			printf("Error: empty BARSeq data!\n");
			exit(EXIT_FAILURE);
		}

		for(ni=0; ni<pBARData->nSeqNum; ni++)
		{
			if(pBARData->vSeqData[ni] != NULL)
			{
				if(pBARData->vSeqData[ni]->nColNum > 0)
				{
					if(pBARData->vSeqData[ni]->vData == NULL)
					{
						printf("Error: empty BARSeq data!\n");
						exit(EXIT_FAILURE);
					}

					for(nj=0; nj<pBARData->vSeqData[ni]->nColNum; nj++)
					{
						if(pColIndicator->pMatElement[nj] == 1)
						{
							DMPLOGTS(pBARData->vSeqData[ni]->vData[nj], dBase);
						}
					}
				}
			}
		}
	}


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_BARData_PowTS()                                                   */
/*  take power of given columns of a BARData object.					   */
/* ----------------------------------------------------------------------- */ 
int Affy_BARData_PowTS(struct tagBARData *pBARData, double dPow, 
		struct INTMATRIX *pColIndicator)
{
	/* define */
	int ni,nj;

	/* init */
	if(pBARData == NULL)
		return PROC_SUCCESS;

	if(pColIndicator == NULL)
	{
		printf("Warning: no column indicators for carrying out BARData operations!\n");
		return PROC_SUCCESS;
	}

	if(pColIndicator->nWidth != pBARData->nColNum)
	{
		printf("Error: column indicator's dimension does not match BARData columns!\n");
		exit(EXIT_FAILURE);
	}

	/* operation */
	if(pBARData->nSeqNum > 0)
	{
		if(pBARData->vSeqData == NULL)
		{
			printf("Error: empty BARSeq data!\n");
			exit(EXIT_FAILURE);
		}

		for(ni=0; ni<pBARData->nSeqNum; ni++)
		{
			if(pBARData->vSeqData[ni] != NULL)
			{
				if(pBARData->vSeqData[ni]->nColNum > 0)
				{
					if(pBARData->vSeqData[ni]->vData == NULL)
					{
						printf("Error: empty BARSeq data!\n");
						exit(EXIT_FAILURE);
					}

					for(nj=0; nj<pBARData->vSeqData[ni]->nColNum; nj++)
					{
						if(pColIndicator->pMatElement[nj] == 1)
						{
							DMPPOWTS(pBARData->vSeqData[ni]->vData[nj], dPow);
						}
					}
				}
			}
		}
	}


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_BARData_UPowTS()                                                  */
/*  take dBase^power where power = given columns of a BARData object.      */
/* ----------------------------------------------------------------------- */ 
int Affy_BARData_UPowTS(struct tagBARData *pBARData, double dBase, 
		struct INTMATRIX *pColIndicator)
{
	/* define */
	int ni,nj;

	/* init */
	if(pBARData == NULL)
		return PROC_SUCCESS;

	if(pColIndicator == NULL)
	{
		printf("Warning: no column indicators for carrying out BARData operations!\n");
		return PROC_SUCCESS;
	}

	if(pColIndicator->nWidth != pBARData->nColNum)
	{
		printf("Error: column indicator's dimension does not match BARData columns!\n");
		exit(EXIT_FAILURE);
	}

	/* operation */
	if(pBARData->nSeqNum > 0)
	{
		if(pBARData->vSeqData == NULL)
		{
			printf("Error: empty BARSeq data!\n");
			exit(EXIT_FAILURE);
		}

		for(ni=0; ni<pBARData->nSeqNum; ni++)
		{
			if(pBARData->vSeqData[ni] != NULL)
			{
				if(pBARData->vSeqData[ni]->nColNum > 0)
				{
					if(pBARData->vSeqData[ni]->vData == NULL)
					{
						printf("Error: empty BARSeq data!\n");
						exit(EXIT_FAILURE);
					}

					for(nj=0; nj<pBARData->vSeqData[ni]->nColNum; nj++)
					{
						if(pColIndicator->pMatElement[nj] == 1)
						{
							DMPUPOWTS(dBase, pBARData->vSeqData[ni]->vData[nj]);
						}
					}
				}
			}
		}
	}


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_BARData_ExpTS()                                                   */
/*  take exponential of given columns of a BARData object.                 */
/* ----------------------------------------------------------------------- */ 
int Affy_BARData_ExpTS(struct tagBARData *pBARData, struct INTMATRIX *pColIndicator)
{
	/* define */
	int ni,nj;

	/* init */
	if(pBARData == NULL)
		return PROC_SUCCESS;

	if(pColIndicator == NULL)
	{
		printf("Warning: no column indicators for carrying out BARData operations!\n");
		return PROC_SUCCESS;
	}

	if(pColIndicator->nWidth != pBARData->nColNum)
	{
		printf("Error: column indicator's dimension does not match BARData columns!\n");
		exit(EXIT_FAILURE);
	}

	/* operation */
	if(pBARData->nSeqNum > 0)
	{
		if(pBARData->vSeqData == NULL)
		{
			printf("Error: empty BARSeq data!\n");
			exit(EXIT_FAILURE);
		}

		for(ni=0; ni<pBARData->nSeqNum; ni++)
		{
			if(pBARData->vSeqData[ni] != NULL)
			{
				if(pBARData->vSeqData[ni]->nColNum > 0)
				{
					if(pBARData->vSeqData[ni]->vData == NULL)
					{
						printf("Error: empty BARSeq data!\n");
						exit(EXIT_FAILURE);
					}

					for(nj=0; nj<pBARData->vSeqData[ni]->nColNum; nj++)
					{
						if(pColIndicator->pMatElement[nj] == 1)
						{
							DMPEXPTS(pBARData->vSeqData[ni]->vData[nj]);
						}
					}
				}
			}
		}
	}


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_BARData_AbsTS()                                                   */
/*  take absolute value of given columns of a BARData object.              */
/* ----------------------------------------------------------------------- */ 
int Affy_BARData_AbsTS(struct tagBARData *pBARData, struct INTMATRIX *pColIndicator)
{
	/* define */
	int ni,nj;

	/* init */
	if(pBARData == NULL)
		return PROC_SUCCESS;

	if(pColIndicator == NULL)
	{
		printf("Warning: no column indicators for carrying out BARData operations!\n");
		return PROC_SUCCESS;
	}

	if(pColIndicator->nWidth != pBARData->nColNum)
	{
		printf("Error: column indicator's dimension does not match BARData columns!\n");
		exit(EXIT_FAILURE);
	}

	/* operation */
	if(pBARData->nSeqNum > 0)
	{
		if(pBARData->vSeqData == NULL)
		{
			printf("Error: empty BARSeq data!\n");
			exit(EXIT_FAILURE);
		}

		for(ni=0; ni<pBARData->nSeqNum; ni++)
		{
			if(pBARData->vSeqData[ni] != NULL)
			{
				if(pBARData->vSeqData[ni]->nColNum > 0)
				{
					if(pBARData->vSeqData[ni]->vData == NULL)
					{
						printf("Error: empty BARSeq data!\n");
						exit(EXIT_FAILURE);
					}

					for(nj=0; nj<pBARData->vSeqData[ni]->nColNum; nj++)
					{
						if(pColIndicator->pMatElement[nj] == 1)
						{
							DMPABSTS(pBARData->vSeqData[ni]->vData[nj]);
						}
					}
				}
			}
		}
	}


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_BARData_Add()                                                     */
/*  Add specified columns of two BARData objects and save results to one   */
/*  of the two BARData objects as specified.                               */
/*  pColIndicator == 1 specifies which column will be involved in the      */
/*   operation;                                                            */
/*  pColIndicator == 2 specifies which column will be used as an internal  */
/*   control to match the two BARData objects so that |pD1-pD2|<1e-6 .     */
/*  If nSaveTo == 0, save results to pD1;                                  */
/*  If nSaveTo == 1, save results to pD2.                                  */
/* ----------------------------------------------------------------------- */ 
int Affy_BARData_Add(struct tagBARData *pD1, struct tagBARData *pD2,
		struct INTMATRIX *pColIndicator, int nSaveTo)
{
	/* define */
	int ni,nj,nk;

	/* init */
	if( (pD1 == NULL) || (pD2 == NULL) )
		return PROC_SUCCESS;

	if(pColIndicator == NULL)
	{
		printf("Warning: no column indicators for carrying out BARData operations!\n");
		return PROC_FAILURE;
	}

	if( (pColIndicator->nWidth != pD1->nColNum) || (pColIndicator->nWidth != pD2->nColNum) )
	{
		printf("Error: column indicator's dimension does not match BARData columns!\n");
		exit(EXIT_FAILURE);
	}

	if( pD1->nSeqNum != pD2->nSeqNum )
	{
		printf("Error: BARData dimension mismatch!\n");
		exit(EXIT_FAILURE);
	}

	if( (pD1->vSeqData == NULL) || (pD2->vSeqData == NULL) )
	{
		if(pD1->nSeqNum > 0)
		{
			printf("Error: empty BARSeq data!\n");
			exit(EXIT_FAILURE);
		}
	}

	for(ni=0; ni<pD2->nSeqNum; ni++)
	{
		if( (pD1->vSeqData[ni] == NULL) || (pD2->vSeqData[ni] == NULL) )
		{
			printf("Error: empty BARSeq data!\n");
			exit(EXIT_FAILURE);
		}

		if(pD1->vSeqData[ni]->nDataNum != pD2->vSeqData[ni]->nDataNum)
		{
			printf("Error: BARData dimension mismatch!\n");
			exit(EXIT_FAILURE);
		}

		for(nj=0; nj<pD2->nColNum; nj++)
		{
			if(pColIndicator->pMatElement[nj] == 1)
			{
				if(nSaveTo == 0)
				{
					DMADDTF(pD1->vSeqData[ni]->vData[nj], pD2->vSeqData[ni]->vData[nj]);
				}
				else if(nSaveTo == 1)
				{
					DMADDTF(pD2->vSeqData[ni]->vData[nj], pD1->vSeqData[ni]->vData[nj]);
				}
				else
				{
					printf("Error: nSaveTo index out of range!\n");
					exit(EXIT_FAILURE);
				}
			}
			else if(pColIndicator->pMatElement[nj] == 2)
			{
				for(nk=0; nk<pD1->vSeqData[ni]->nDataNum; nk++)
				{
					if(fabs(pD1->vSeqData[ni]->vData[nj]->pMatElement[nk]-pD2->vSeqData[ni]->vData[nj]->pMatElement[nk]) > 1e-6)
					{
						printf("Error: internal matching check failed!\n");
						exit(EXIT_FAILURE);
					}
				}
			}
		}
	}

	/* return */
	return PROC_SUCCESS;
}

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
		struct INTMATRIX *pColIndicator, int nSaveTo)
{
	/* define */
	int ni,nj,nk;

	/* init */
	if( (pD1 == NULL) || (pD2 == NULL) )
		return PROC_SUCCESS;

	if(pColIndicator == NULL)
	{
		printf("Warning: no column indicators for carrying out BARData operations!\n");
		return PROC_FAILURE;
	}

	if( (pColIndicator->nWidth != pD1->nColNum) || (pColIndicator->nWidth != pD2->nColNum) )
	{
		printf("Error: column indicator's dimension does not match BARData columns!\n");
		exit(EXIT_FAILURE);
	}

	if( pD1->nSeqNum != pD2->nSeqNum )
	{
		printf("Error: BARData dimension mismatch!\n");
		exit(EXIT_FAILURE);
	}

	if( (pD1->vSeqData == NULL) || (pD2->vSeqData == NULL) )
	{
		if(pD1->nSeqNum > 0)
		{
			printf("Error: empty BARSeq data!\n");
			exit(EXIT_FAILURE);
		}
	}

	for(ni=0; ni<pD2->nSeqNum; ni++)
	{
		if( (pD1->vSeqData[ni] == NULL) || (pD2->vSeqData[ni] == NULL) )
		{
			printf("Error: empty BARSeq data!\n");
			exit(EXIT_FAILURE);
		}

		if(pD1->vSeqData[ni]->nDataNum != pD2->vSeqData[ni]->nDataNum)
		{
			printf("Error: BARData dimension mismatch!\n");
			exit(EXIT_FAILURE);
		}

		for(nj=0; nj<pD2->nColNum; nj++)
		{
			if(pColIndicator->pMatElement[nj] == 1)
			{
				if(nSaveTo == 0)
				{
					DMSUBTF(pD1->vSeqData[ni]->vData[nj], pD2->vSeqData[ni]->vData[nj]);
				}
				else if(nSaveTo == 1)
				{
					DMSUBTL(pD1->vSeqData[ni]->vData[nj], pD2->vSeqData[ni]->vData[nj]);
				}
				else
				{
					printf("Error: nSaveTo index out of range!\n");
					exit(EXIT_FAILURE);
				}
			}
			else if(pColIndicator->pMatElement[nj] == 2)
			{
				for(nk=0; nk<pD1->vSeqData[ni]->nDataNum; nk++)
				{
					if(fabs(pD1->vSeqData[ni]->vData[nj]->pMatElement[nk]-pD2->vSeqData[ni]->vData[nj]->pMatElement[nk]) > 1e-6)
					{
						printf("Error: internal matching check failed!\n");
						exit(EXIT_FAILURE);
					}
				}
			}
		}
	}

	/* return */
	return PROC_SUCCESS;
}

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
		struct INTMATRIX *pColIndicator, int nSaveTo)
{
	/* define */
	int ni,nj,nk;

	/* init */
	if( (pD1 == NULL) || (pD2 == NULL) )
		return PROC_SUCCESS;

	if(pColIndicator == NULL)
	{
		printf("Warning: no column indicators for carrying out BARData operations!\n");
		return PROC_FAILURE;
	}

	if( (pColIndicator->nWidth != pD1->nColNum) || (pColIndicator->nWidth != pD2->nColNum) )
	{
		printf("Error: column indicator's dimension does not match BARData columns!\n");
		exit(EXIT_FAILURE);
	}

	if( pD1->nSeqNum != pD2->nSeqNum )
	{
		printf("Error: BARData dimension mismatch!\n");
		exit(EXIT_FAILURE);
	}

	if( (pD1->vSeqData == NULL) || (pD2->vSeqData == NULL) )
	{
		if(pD1->nSeqNum > 0)
		{
			printf("Error: empty BARSeq data!\n");
			exit(EXIT_FAILURE);
		}
	}

	for(ni=0; ni<pD2->nSeqNum; ni++)
	{
		if( (pD1->vSeqData[ni] == NULL) || (pD2->vSeqData[ni] == NULL) )
		{
			printf("Error: empty BARSeq data!\n");
			exit(EXIT_FAILURE);
		}

		if(pD1->vSeqData[ni]->nDataNum != pD2->vSeqData[ni]->nDataNum)
		{
			printf("Error: BARData dimension mismatch!\n");
			exit(EXIT_FAILURE);
		}

		for(nj=0; nj<pD2->nColNum; nj++)
		{
			if(pColIndicator->pMatElement[nj] == 1)
			{
				if(nSaveTo == 0)
				{
					DMPMULTF(pD1->vSeqData[ni]->vData[nj], pD2->vSeqData[ni]->vData[nj]);
				}
				else if(nSaveTo == 1)
				{
					DMPMULTF(pD2->vSeqData[ni]->vData[nj], pD1->vSeqData[ni]->vData[nj]);
				}
				else
				{
					printf("Error: nSaveTo index out of range!\n");
					exit(EXIT_FAILURE);
				}
			}
			else if(pColIndicator->pMatElement[nj] == 2)
			{
				for(nk=0; nk<pD1->vSeqData[ni]->nDataNum; nk++)
				{
					if(fabs(pD1->vSeqData[ni]->vData[nj]->pMatElement[nk]-pD2->vSeqData[ni]->vData[nj]->pMatElement[nk]) > 1e-6)
					{
						printf("Error: internal matching check failed!\n");
						exit(EXIT_FAILURE);
					}
				}
			}
		}
	}

	/* return */
	return PROC_SUCCESS;
}

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
		struct INTMATRIX *pColIndicator, int nSaveTo)
{
	/* define */
	int ni,nj,nk;

	/* init */
	if( (pD1 == NULL) || (pD2 == NULL) )
		return PROC_SUCCESS;

	if(pColIndicator == NULL)
	{
		printf("Warning: no column indicators for carrying out BARData operations!\n");
		return PROC_FAILURE;
	}

	if( (pColIndicator->nWidth != pD1->nColNum) || (pColIndicator->nWidth != pD2->nColNum) )
	{
		printf("Error: column indicator's dimension does not match BARData columns!\n");
		exit(EXIT_FAILURE);
	}

	if( pD1->nSeqNum != pD2->nSeqNum )
	{
		printf("Error: BARData dimension mismatch!\n");
		exit(EXIT_FAILURE);
	}

	if( (pD1->vSeqData == NULL) || (pD2->vSeqData == NULL) )
	{
		if(pD1->nSeqNum > 0)
		{
			printf("Error: empty BARSeq data!\n");
			exit(EXIT_FAILURE);
		}
	}

	for(ni=0; ni<pD2->nSeqNum; ni++)
	{
		if( (pD1->vSeqData[ni] == NULL) || (pD2->vSeqData[ni] == NULL) )
		{
			printf("Error: empty BARSeq data!\n");
			exit(EXIT_FAILURE);
		}

		if(pD1->vSeqData[ni]->nDataNum != pD2->vSeqData[ni]->nDataNum)
		{
			printf("Error: BARData dimension mismatch!\n");
			exit(EXIT_FAILURE);
		}

		for(nj=0; nj<pD2->nColNum; nj++)
		{
			if(pColIndicator->pMatElement[nj] == 1)
			{
				if(nSaveTo == 0)
				{
					DMPDIVTF(pD1->vSeqData[ni]->vData[nj], pD2->vSeqData[ni]->vData[nj]);
				}
				else if(nSaveTo == 1)
				{
					DMPDIVTL(pD1->vSeqData[ni]->vData[nj], pD2->vSeqData[ni]->vData[nj]);
				}
				else
				{
					printf("Error: nSaveTo index out of range!\n");
					exit(EXIT_FAILURE);
				}
			}
			else if(pColIndicator->pMatElement[nj] == 2)
			{
				for(nk=0; nk<pD1->vSeqData[ni]->nDataNum; nk++)
				{
					if(fabs(pD1->vSeqData[ni]->vData[nj]->pMatElement[nk]-pD2->vSeqData[ni]->vData[nj]->pMatElement[nk]) > 1e-6)
					{
						printf("Error: internal matching check failed!\n");
						exit(EXIT_FAILURE);
					}
				}
			}
		}
	}

	/* return */
	return PROC_SUCCESS;
}

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
		struct INTMATRIX *pColIndicator, int nSaveTo)
{
	/* define */
	int ni,nj,nk;

	/* init */
	if( (pD1 == NULL) || (pD2 == NULL) )
		return PROC_SUCCESS;

	if(pColIndicator == NULL)
	{
		printf("Warning: no column indicators for carrying out BARData operations!\n");
		return PROC_FAILURE;
	}

	if( (pColIndicator->nWidth != pD1->nColNum) || (pColIndicator->nWidth != pD2->nColNum) )
	{
		printf("Error: column indicator's dimension does not match BARData columns!\n");
		exit(EXIT_FAILURE);
	}

	if( pD1->nSeqNum != pD2->nSeqNum )
	{
		printf("Error: BARData dimension mismatch!\n");
		exit(EXIT_FAILURE);
	}

	if( (pD1->vSeqData == NULL) || (pD2->vSeqData == NULL) )
	{
		if(pD1->nSeqNum > 0)
		{
			printf("Error: empty BARSeq data!\n");
			exit(EXIT_FAILURE);
		}
	}

	for(ni=0; ni<pD2->nSeqNum; ni++)
	{
		if( (pD1->vSeqData[ni] == NULL) || (pD2->vSeqData[ni] == NULL) )
		{
			printf("Error: empty BARSeq data!\n");
			exit(EXIT_FAILURE);
		}

		if(pD1->vSeqData[ni]->nDataNum != pD2->vSeqData[ni]->nDataNum)
		{
			printf("Error: BARData dimension mismatch!\n");
			exit(EXIT_FAILURE);
		}

		for(nj=0; nj<pD2->nColNum; nj++)
		{
			if(pColIndicator->pMatElement[nj] == 1)
			{
				if(nSaveTo == 0)
				{
					DMPPOWTF(pD1->vSeqData[ni]->vData[nj], pD2->vSeqData[ni]->vData[nj]);
				}
				else if(nSaveTo == 1)
				{
					DMPPOWTL(pD1->vSeqData[ni]->vData[nj], pD2->vSeqData[ni]->vData[nj]);
				}
				else
				{
					printf("Error: nSaveTo index out of range!\n");
					exit(EXIT_FAILURE);
				}
			}
			else if(pColIndicator->pMatElement[nj] == 2)
			{
				for(nk=0; nk<pD1->vSeqData[ni]->nDataNum; nk++)
				{
					if(fabs(pD1->vSeqData[ni]->vData[nj]->pMatElement[nk]-pD2->vSeqData[ni]->vData[nj]->pMatElement[nk]) > 1e-6)
					{
						printf("Error: internal matching check failed!\n");
						exit(EXIT_FAILURE);
					}
				}
			}
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_GenericFileObject_Create()                                        */
/*  Create a object for loading generic file object.                       */
/* ----------------------------------------------------------------------- */ 
struct tagGenericFileObject *Affy_GenericFileObject_Create()
{
	struct tagGenericFileObject *pNewData = NULL;

	pNewData = (struct tagGenericFileObject *)calloc(1, sizeof(struct tagGenericFileObject));
	if(pNewData == NULL)
	{
		printf("Error: Affy_GenericFileData_Create, cannot create generic file data object!\n");
		exit(EXIT_FAILURE);
	}

	pNewData->pFileHeader = NULL;
	pNewData->pDataHeader = NULL;
	pNewData->vDataGroup = NULL;

	return pNewData;
};


/* ----------------------------------------------------------------------- */ 
/*  Affy_GenericFileObject_Delete()                                        */
/*  Delete an object for loading generic file object.                      */
/* ----------------------------------------------------------------------- */ 
void Affy_GenericFileObject_Delete(struct tagGenericFileObject **pGenericObj)
{
	/* define */
	int ni;

	/* free memory */
	if(pGenericObj == NULL)
		return;

	if(*pGenericObj == NULL)
		return;

	if((*pGenericObj)->pFileHeader != NULL)
	{
		for(ni=0; ni<(*pGenericObj)->pFileHeader->nDataGroupNum; ni++)
		{
			Affy_GenericDataGroup_Delete((*pGenericObj)->vDataGroup+ni);
		}
		free((*pGenericObj)->vDataGroup);
	}

	Affy_GenericFileHeader_Delete(&((*pGenericObj)->pFileHeader));
	Affy_GenericDataHeader_Delete(&((*pGenericObj)->pDataHeader));

	free(*pGenericObj);

	*pGenericObj = NULL;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_GenericFileHeader_Create()                                        */
/*  Create a generic file header object.                                   */
/* ----------------------------------------------------------------------- */ 
struct tagGenericFileHeader *Affy_GenericFileHeader_Create()
{
	struct tagGenericFileHeader *pNewHeader = NULL;

	pNewHeader = (struct tagGenericFileHeader *)calloc(1, sizeof(struct tagGenericFileHeader));
	if(pNewHeader == NULL)
	{
		printf("Error: Affy_GenericFileFileHeader_Create, cannot create generic file file header object!\n");
		exit(EXIT_FAILURE);
	}

	return pNewHeader;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_GenericFileHeader_Delete()                                        */
/*  Delete a generic file header object.                                   */
/* ----------------------------------------------------------------------- */ 
void Affy_GenericFileHeader_Delete(struct tagGenericFileHeader **pHeader)
{
	/* free memory */
	if(pHeader == NULL)
		return;

	if(*pHeader == NULL)
		return;

	free(*pHeader);

	*pHeader = NULL;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_GenericDataHeader_Create()                                        */
/*  Create a generic data header object.                                   */
/* ----------------------------------------------------------------------- */ 
struct tagGenericDataHeader *Affy_GenericDataHeader_Create()
{
	struct tagGenericDataHeader *pNewHeader = NULL;

	pNewHeader = (struct tagGenericDataHeader *)calloc(1, sizeof(struct tagGenericDataHeader));
	if(pNewHeader == NULL)
	{
		printf("Error: Affy_GenericFileDataHeader_Create, cannot create generic file data header object!\n");
		exit(EXIT_FAILURE);
	}

	return pNewHeader;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_GenericDataHeader_Delete()                                        */
/*  Delete a generic  data header object.                                  */
/* ----------------------------------------------------------------------- */ 
void Affy_GenericDataHeader_Delete(struct tagGenericDataHeader **pHeader)
{
	/* define */
	int ni;

	/* free memory */
	if(pHeader == NULL)
		return;

	if(*pHeader == NULL)
		return;

	for(ni=0; ni<(*pHeader)->nParentNum; ni++)
	{
		Affy_GenericDataHeader_Delete((*pHeader)->vParent+ni);
	}
	free((*pHeader)->vParent);

	DeleteString((*pHeader)->pDataTypeID);
	(*pHeader)->pDataTypeID = NULL;

	DeleteString((*pHeader)->pUniqueFileID);
	(*pHeader)->pUniqueFileID = NULL;

	DeleteWString(&((*pHeader)->pDateTime));
	DeleteWString(&((*pHeader)->pLocale));

	for(ni=0; ni<(*pHeader)->nParamNum; ni++)
	{
		DeleteWString((*pHeader)->vParamName+ni);
		DeleteString((*pHeader)->vParamValue[ni]);
		(*pHeader)->vParamValue[ni] = NULL;
		DeleteWString((*pHeader)->vParamType+ni);
	}
	free((*pHeader)->vParamName);
	free((*pHeader)->vParamValue);
	free((*pHeader)->vParamType);

	free(*pHeader);

	*pHeader = NULL;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_GenericDataGroup_Create()                                         */
/*  Create a generic file data group object.                               */
/* ----------------------------------------------------------------------- */ 
struct tagGenericDataGroup *Affy_GenericDataGroup_Create()
{
	struct tagGenericDataGroup *pNewDataGroup = NULL;

	pNewDataGroup = (struct tagGenericDataGroup *)calloc(1, sizeof(struct tagGenericDataGroup));
	if(pNewDataGroup == NULL)
	{
		printf("Error: Affy_GenericFileDataGroup_Create, cannot create generic file data group object!\n");
		exit(EXIT_FAILURE);
	}

	return pNewDataGroup;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_GenericDataGroup_Delete()                                         */
/*  Delete a generic file data group object.                               */
/* ----------------------------------------------------------------------- */ 
void Affy_GenericDataGroup_Delete(struct tagGenericDataGroup **pDataGroup)
{
	/* define */
	int ni;

	/* free memory */
	if(pDataGroup == NULL)
		return;

	if(*pDataGroup == NULL)
		return;

	DeleteWString(&((*pDataGroup)->pDataGroupName));
	for(ni=0; ni<(*pDataGroup)->nDataSetNum; ni++)
	{
		Affy_GenericDataSet_Delete((*pDataGroup)->vDataSets+ni);
	}
	free((*pDataGroup)->vDataSets);

	free(*pDataGroup);

	*pDataGroup = NULL;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_GenericDataSet_Create()                                           */
/*  Create a generic file data set object.                                 */
/* ----------------------------------------------------------------------- */ 
struct tagGenericDataSet *Affy_GenericDataSet_Create()
{
	struct tagGenericDataSet *pNewDataSet = NULL;

	pNewDataSet = (struct tagGenericDataSet *)calloc(1, sizeof(struct tagGenericDataSet));
	if(pNewDataSet == NULL)
	{
		printf("Error: Affy_GenericFileDataSet_Create, cannot create generic file data set object!\n");
		exit(EXIT_FAILURE);
	}

	return pNewDataSet;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_GenericDataSet_Delete()                                           */
/*  Delete a generic file data set object.                                 */
/* ----------------------------------------------------------------------- */ 
void Affy_GenericDataSet_Delete(struct tagGenericDataSet **pDataSet)
{
	/* define */
	int ni;
	unsigned int nj;

	/* free memory */
	if(pDataSet == NULL)
		return;

	if(*pDataSet == NULL)
		return;

	DeleteWString(&((*pDataSet)->pDataSetName));
	for(ni=0; ni<(*pDataSet)->nParamNum; ni++)
	{
		DeleteWString((*pDataSet)->vParamName+ni);
		DeleteString((*pDataSet)->vParamValue[ni]);
		(*pDataSet)->vParamValue[ni] = NULL;
		DeleteWString((*pDataSet)->vParamType+ni);
	}
	free((*pDataSet)->vParamName);
	free((*pDataSet)->vParamValue);
	free((*pDataSet)->vParamType);

	for(nj=0; nj<(*pDataSet)->nColNum; nj++)
	{
		DeleteWString((*pDataSet)->vColName+nj);
	}
	free((*pDataSet)->vColName);
	free((*pDataSet)->vColType);
	free((*pDataSet)->vColSize);
	free((*pDataSet)->vData);

	free(*pDataSet);

	*pDataSet = NULL;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_GenericFileObject_Load()                                          */
/*  Load a generic file object from a file.                                */
/* ----------------------------------------------------------------------- */ 
struct tagGenericFileObject *Affy_GenericFileObject_Load(char strFilePath[])
{
	/* define */
	struct tagGenericFileObject *pGenericObj = NULL;

	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));

	int ni;	
	FILE *fpIn;
	/* char chbuf[1]; */

	/* load file */
	fpIn = NULL;
	fpIn = fopen(strFilePath, "rb");
	if(fpIn == NULL)
	{
		printf("Error: Affy_GenericFileObject_Load, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}

	/* prepare object */
	pGenericObj = Affy_GenericFileObject_Create();
	if(pGenericObj == NULL)
	{
		printf("Error: Affy_GenericFileObject_Load, cannot create buffer for loading generic file data!\n");
		exit(EXIT_FAILURE);
	}

	/* load file header */
	/* printf("Load Generic File Header!\n"); */
	pGenericObj->pFileHeader = Affy_GenericFile_LoadFileHeader(fpIn);
	if(pGenericObj->pFileHeader == NULL)
	{
		printf("Error: Affy_GenericFileObject_Load, cannot load generic file header!\n");
		Affy_GenericFileObject_Delete(&pGenericObj);
		exit(EXIT_FAILURE);
	}
	
	pGenericObj->vDataGroup = (struct tagGenericDataGroup **)calloc(pGenericObj->pFileHeader->nDataGroupNum, sizeof(struct tagGenericDataGroup *));
	if(pGenericObj->vDataGroup == NULL)
	{
		printf("Error: Affy_GenericFileObject_Load, cannot allocate memory for generic data groups!\n");
		pGenericObj->pFileHeader->nDataGroupNum = 0;
		Affy_GenericFileObject_Delete(&pGenericObj);
		exit(EXIT_FAILURE);
	}

	/* load data header */
	/* printf("Load Generic Data Header!\n"); */
	pGenericObj->pDataHeader = Affy_GenericFile_LoadDataHeader(fpIn);
	if(pGenericObj->pDataHeader == NULL)
	{
		printf("Error: Affy_GenericFileObject_Load, cannot load generic data header!\n");
		Affy_GenericFileObject_Delete(&pGenericObj);
		exit(EXIT_FAILURE);
	}
	
	/* load data group */
	/* printf("Load Generic Data Group!\n"); */
	if(ftell(fpIn) != pGenericObj->pFileHeader->nFirstDataGroupPos)
	{
		printf("Warning: file pointer not match!\n");
	}

	if(fseek(fpIn, pGenericObj->pFileHeader->nFirstDataGroupPos, SEEK_SET) != 0)
	{
		printf("Error: Affy_GenericFileObject_Load, cannot load generic file data group!\n");
		Affy_GenericFileObject_Delete(&pGenericObj);
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<pGenericObj->pFileHeader->nDataGroupNum; ni++)
	{
		pGenericObj->vDataGroup[ni] = Affy_GenericFile_LoadDataGroup(fpIn);
		
		if(pGenericObj->vDataGroup[ni] == NULL)
		{
			printf("Error: Affy_GenericFileObject_Load, cannot load generic file data group!\n");
			Affy_GenericFileObject_Delete(&pGenericObj);
			exit(EXIT_FAILURE);
		}

		if(ni<pGenericObj->pFileHeader->nDataGroupNum-1)
		{
			if(fseek(fpIn, pGenericObj->vDataGroup[ni]->nNextDataGroupPos, SEEK_SET) != 0)
			{
				printf("Error: Affy_GenericFileObject_Load, cannot load generic file data group correctly!\n");
				Affy_GenericFileObject_Delete(&pGenericObj);
				exit(EXIT_FAILURE);
			}
		}
		else if(ni == pGenericObj->pFileHeader->nDataGroupNum-1)
		{
			if( (pGenericObj->vDataGroup[ni]->nNextDataGroupPos != 0) &&
				(pGenericObj->vDataGroup[ni]->nNextDataGroupPos != ftell(fpIn)) )
			{
				printf("Error: Affy_GenericFileObject_Load, cannot load generic file data group correctly!\n");
				Affy_GenericFileObject_Delete(&pGenericObj);
				exit(EXIT_FAILURE);
			}
		}

		/* if(ni<pGenericObj->pFileHeader->nDataGroupNum-1)
		{
			if(ftell(fpIn) != pGenericObj->vDataGroup[ni]->nNextDataGroupPos)
			{
				printf("Error: Affy_GenericFileObject_Load, cannot load generic file data group correctly!\n");
				Affy_GenericFileObject_Delete(&pGenericObj);
				exit(EXIT_FAILURE);
			}
		}
		else if(ni == pGenericObj->pFileHeader->nDataGroupNum-1)
		{
			if(pGenericObj->vDataGroup[ni]->nNextDataGroupPos != 0)
			{
				printf("Error: Affy_GenericFileObject_Load, cannot load generic file data group correctly!\n");
				Affy_GenericFileObject_Delete(&pGenericObj);
				exit(EXIT_FAILURE);
			}
		}*/
	}

	/* close file */
	fclose(fpIn);

	/* return */
	return pGenericObj;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_GenericFile_LoadFileHeader()                                      */
/*  Load generic file header.                                              */
/* ----------------------------------------------------------------------- */ 
struct tagGenericFileHeader *Affy_GenericFile_LoadFileHeader(FILE *fpIn)
{
	/* define */
	struct tagGenericFileHeader *pHeader = NULL;
	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));

	
	/* create header */
	pHeader = Affy_GenericFileHeader_Create();
	if(pHeader == NULL)
	{
		printf("Warning: Affy_GenericFile_LoadFileHeader, cannot create buffer for file header!\n");
		return pHeader;
	}

	/* load magic number */
	if(fread(&(pHeader->nMagicNumber), 1, 1, fpIn) != 1)
	{
		printf("Error: Affy_GenericFile_LoadFileHeader, cannot load magic number!\n");
		Affy_GenericFileHeader_Delete(&pHeader);
		return NULL;
	}

	if(pHeader->nMagicNumber != 59)
	{
		printf("Error: Affy_GenericFile_LoadFileHeader, cannot load magic number correctly!\n");
		Affy_GenericFileHeader_Delete(&pHeader);
		return NULL;
	}

	/* load version number */
	if(fread(&(pHeader->nVersionNumber), 1, 1, fpIn) != 1)
	{
		printf("Error: Affy_GenericFile_LoadFileHeader, cannot load version number!\n");
		Affy_GenericFileHeader_Delete(&pHeader);
		return NULL;
	}

	if(pHeader->nVersionNumber != 1)
	{
		printf("Error: Affy_GenericFile_LoadFileHeader, sorry currently we do not support verion %d!\n", pHeader->nVersionNumber);
		Affy_GenericFileHeader_Delete(&pHeader);
		return NULL;
	}

	/* load data group number */
	if(big_endian_fread(&(pHeader->nDataGroupNum), INT_SIZE, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: Affy_GenericFile_LoadFileHeader, cannot load data group number!\n");
		Affy_GenericFileHeader_Delete(&pHeader);
		return NULL;
	}
	
	/* load first data group pos */
	if(big_endian_fread(&(pHeader->nFirstDataGroupPos), DWORD_SIZE, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: Affy_GenericFile_LoadFileHeader, cannot load file position of the first data group!\n");
		Affy_GenericFileHeader_Delete(&pHeader);
		return NULL;
	}


	/* return */
	return pHeader;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_GenericFile_LoadDataGroup()                                       */
/*  Load generic data group.                                               */
/* ----------------------------------------------------------------------- */ 
struct tagGenericDataHeader *Affy_GenericFile_LoadDataHeader(FILE *fpIn)
{
	/* define */
	struct tagGenericDataHeader *pHeader = NULL;
	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));
	int ni;
	
	/* create header */
	/* printf("Generic Header!\n"); */
	pHeader = Affy_GenericDataHeader_Create();
	if(pHeader == NULL)
	{
		printf("Warning: Affy_GenericFile_LoadDataHeader, cannot create buffer for data header!\n");
		return pHeader;
	}

	/* load data type identifier */
	pHeader->pDataTypeID = Affy_GenericFile_LoadString(fpIn, little_endian_machine);
	/* printf("%s\n", pHeader->pDataTypeID->m_pString); */

	/* load unique file identifier */
	pHeader->pUniqueFileID = Affy_GenericFile_LoadString(fpIn, little_endian_machine);
	/* printf("%s\n", pHeader->pUniqueFileID->m_pString); */

	/* load date time */
	pHeader->pDateTime = Affy_GenericFile_LoadWString(fpIn, little_endian_machine);
	/* printf("%S\n", pHeader->pDateTime->m_pWString); */

	/* load locale */
	pHeader->pLocale = Affy_GenericFile_LoadWString(fpIn, little_endian_machine);
	/* printf("%S\n", pHeader->pLocale->m_pWString); */

	/* load number of parameters */
	if(big_endian_fread(&(pHeader->nParamNum), INT_SIZE, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: Affy_GenericFile_LoadDataHeader, cannot load number of parameters!\n");
		Affy_GenericDataHeader_Delete(&pHeader);
		return NULL;
	}

	/* printf("Param Num = %d!\n", pHeader->nParamNum); */

	if(pHeader->nParamNum > 0)
	{
		pHeader->vParamName = (struct tagWString **)calloc(pHeader->nParamNum, sizeof(struct tagWString *));
		pHeader->vParamValue = (struct tagString **)calloc(pHeader->nParamNum, sizeof(struct tagString *));
		pHeader->vParamType = (struct tagWString **)calloc(pHeader->nParamNum, sizeof(struct tagWString *));

		if((pHeader->vParamName == NULL) || (pHeader->vParamType == NULL) || (pHeader->vParamValue == NULL))
		{
			printf("Error: Affy_GenericFile_LoadDataHeader, cannot allocate memory for loading parameters!\n");
			pHeader->nParamNum = 0;
			Affy_GenericDataHeader_Delete(&pHeader);
			return NULL;
		}

		/* load parameters */
		for(ni=0; ni<pHeader->nParamNum; ni++)
		{
			pHeader->vParamName[ni] = Affy_GenericFile_LoadWString(fpIn, little_endian_machine);
			pHeader->vParamValue[ni] = Affy_GenericFile_LoadString(fpIn, little_endian_machine);
			pHeader->vParamType[ni] = Affy_GenericFile_LoadWString(fpIn, little_endian_machine);
		}
	}

	/* load number of parent headers */
	if(big_endian_fread(&(pHeader->nParentNum), INT_SIZE, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: Affy_GenericFile_LoadDataHeader, cannot load number of parents!\n");
		Affy_GenericDataHeader_Delete(&pHeader);
		return NULL;
	}
	/* printf("Parent Num = %d!\n", pHeader->nParentNum); */

	/* load parent header */
	if(pHeader->nParentNum > 0)
	{
		pHeader->vParent = (struct tagGenericDataHeader **)calloc(pHeader->nParentNum, sizeof(struct tagGenericDataHeader *));
		if(pHeader->vParent == NULL)
		{
			printf("Error: Affy_GenericFile_LoadDataHeader, cannot allocate memory for parent headers!\n");
			pHeader->nParentNum = 0;
			Affy_GenericDataHeader_Delete(&pHeader);
			return NULL;
		}

		for(ni=0; ni<pHeader->nParentNum; ni++)
		{
			pHeader->vParent[ni] = Affy_GenericFile_LoadDataHeader(fpIn);
			if(pHeader->vParent[ni] == NULL)
			{
				printf("Error: Affy_GenericFile_LoadDataHeader, cannot load parent data header!\n");
				Affy_GenericDataHeader_Delete(&pHeader);
				return NULL;
			}
		}
	}
	
	/* return */
	return pHeader;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_GenericFile_LoadDataGroup()                                       */
/*  Load generic data group.                                               */
/* ----------------------------------------------------------------------- */ 
struct tagGenericDataGroup *Affy_GenericFile_LoadDataGroup(FILE *fpIn)
{
	/* define */
	struct tagGenericDataGroup *pDataGroup = NULL;

	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));
	int ni;
	/* char chbuf[1]; */
	int nSeekRet;
	
	/* create data group */
	pDataGroup = Affy_GenericDataGroup_Create();
	if(pDataGroup == NULL)
	{
		printf("Warning: Affy_GenericFile_LoadDataGroup, cannot create buffer for data group!\n");
		return pDataGroup;
	}

	/* load next data group pos */
	if(big_endian_fread(&(pDataGroup->nNextDataGroupPos), DWORD_SIZE, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: Affy_GenericFile_LoadDataGroup, cannot load position of the next data group!\n");
		Affy_GenericDataGroup_Delete(&pDataGroup);
		return NULL;
	}
	
	/* load first data set pos */
	if(big_endian_fread(&(pDataGroup->nFirstDataSetPos), DWORD_SIZE, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: Affy_GenericFile_LoadDataGroup, cannot load position of the first dataset!\n");
		Affy_GenericDataGroup_Delete(&pDataGroup);
		return NULL;
	}

	/* load data set number */
	if(big_endian_fread(&(pDataGroup->nDataSetNum), INT_SIZE, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: Affy_GenericFile_LoadDataGroup, cannot load dataset number!\n");
		Affy_GenericDataGroup_Delete(&pDataGroup);
		return NULL;
	}

	if(pDataGroup->nDataSetNum < 0)
	{
		printf("Error: Affy_GenericFile_LoadDataGroup, cannot load dataset number!\n");
		pDataGroup->nDataSetNum = 0;
		Affy_GenericDataGroup_Delete(&pDataGroup);
		return NULL;
	}

	pDataGroup->vDataSets = (struct tagGenericDataSet **)calloc(pDataGroup->nDataSetNum, sizeof(struct tagGenericDataSet *));
	if(pDataGroup->vDataSets == NULL)
	{
		printf("Error: Affy_GenericFile_LoadDataGroup, cannot allocate memory for data sets!\n");
		pDataGroup->nDataSetNum = 0;
		Affy_GenericDataGroup_Delete(&pDataGroup);
		return NULL;
	}

	/* load data group name */
	pDataGroup->pDataGroupName = Affy_GenericFile_LoadWString(fpIn, little_endian_machine);

	if(pDataGroup->nDataSetNum > 0)
	{
		if(ftell(fpIn) != pDataGroup->nFirstDataSetPos)
		{
			printf("Error: Affy_GenericFile_LoadDataGroup, cannot match position of the first data set correctly!\n");
			Affy_GenericDataGroup_Delete(&pDataGroup);
			return NULL;
		}
	}

	/* load data sets */
	for(ni=0; ni<pDataGroup->nDataSetNum; ni++)
	{
		pDataGroup->vDataSets[ni] = Affy_GenericFile_LoadDataSet(fpIn);
		if(pDataGroup->vDataSets[ni] == NULL)
		{
			printf("Error: Affy_GenericFile_LoadDataGroup, cannot load data set!\n");
			Affy_GenericDataGroup_Delete(&pDataGroup);
			return NULL;
		}

		nSeekRet = fseek(fpIn, pDataGroup->vDataSets[ni]->nNextDataSetPos, SEEK_SET);
		if(ni<pDataGroup->nDataSetNum-1)
		{
			if(nSeekRet != 0)
			{
				printf("Error: Affy_GenericFile_LoadDataGroup, cannot load data set correctly!\n");
				Affy_GenericDataGroup_Delete(&pDataGroup);
				return NULL;
			}
		}
		
		/* if(ni<pDataGroup->nDataSetNum-1)
		{
			if(pDataGroup->vDataSets[ni]->nNextDataSetPos != ftell(fpIn))
			{
				printf("Error: Affy_GenericFile_LoadDataGroup, cannot load data set correctly, file pointer offset wrong!\n");
				Affy_GenericDataGroup_Delete(&pDataGroup);
				return NULL;
			}
		}*/
	}

	/* if(fread(chbuf, 1, 1, fpIn) != 1)
	{
		printf("Error: Affy_GenericFile_LoadDataGroup, cannot load data correctly!\n");
		Affy_GenericDataGroup_Delete(&pDataGroup);
		return NULL;
	}*/

	/* return */
	return pDataGroup;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_GenericFile_LoadDataSet()                                         */
/*  Load generic data set.                                                 */
/* ----------------------------------------------------------------------- */ 
struct tagGenericDataSet *Affy_GenericFile_LoadDataSet(FILE *fpIn)
{
	/* define */
	struct tagGenericDataSet *pDataSet = NULL;

	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));
	int ni,nj,nk;
	int unitsize;
	int memsize;
	char *buf;
	int nLen;
	/* char chbuf[1]; */

	
	/* create data set */
	pDataSet = Affy_GenericDataSet_Create();
	if(pDataSet == NULL)
	{
		printf("Warning: Affy_GenericFile_LoadDataSet, cannot create buffer for data set!\n");
		return pDataSet;
	}

	/* load position of the first data element */
	if(big_endian_fread(&(pDataSet->nFirstDataElementPos), DWORD_SIZE, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: Affy_GenericFile_LoadDataSet, cannot load position of the first data element!\n");
		Affy_GenericDataSet_Delete(&pDataSet);
		return NULL;
	}

	/* load position of the next data set */
	if(big_endian_fread(&(pDataSet->nNextDataSetPos), DWORD_SIZE, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: Affy_GenericFile_LoadDataSet, cannot load position of the next dataset!\n");
		Affy_GenericDataSet_Delete(&pDataSet);
		return NULL;
	}
	
	/* load data set name */
	pDataSet->pDataSetName = Affy_GenericFile_LoadWString(fpIn, little_endian_machine);

	/* load number of parameters */
	if(big_endian_fread(&(pDataSet->nParamNum), INT_SIZE, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: Affy_GenericFile_LoadDataSet, cannot load parameter number!\n");
		Affy_GenericDataSet_Delete(&pDataSet);
		return NULL;
	}

	if(pDataSet->nParamNum < 0)
	{
		printf("Error: Affy_GenericFile_LoadDataSet, parameter number < 0!\n");
		pDataSet->nParamNum = 0;
		Affy_GenericDataSet_Delete(&pDataSet);
		return NULL;
	}

	if(pDataSet->nParamNum > 0)
	{
		pDataSet->vParamName = (struct tagWString **)calloc(pDataSet->nParamNum, sizeof(struct tagWString *));
		pDataSet->vParamValue = (struct tagString **)calloc(pDataSet->nParamNum, sizeof(struct tagString *));
		pDataSet->vParamType = (struct tagWString **)calloc(pDataSet->nParamNum, sizeof(struct tagWString *));

		if((pDataSet->vParamName == NULL) || (pDataSet->vParamType == NULL) || (pDataSet->vParamValue == NULL))
		{
			printf("Error: Affy_GenericFile_LoadDataSet, cannot allocate memory for loading parameters!\n");
			pDataSet->nParamNum = 0;
			Affy_GenericDataSet_Delete(&pDataSet);
			return NULL;
		}

		/* load parameters */
		for(ni=0; ni<pDataSet->nParamNum; ni++)
		{
			pDataSet->vParamName[ni] = Affy_GenericFile_LoadWString(fpIn, little_endian_machine);
			pDataSet->vParamValue[ni] = Affy_GenericFile_LoadString(fpIn, little_endian_machine);
			pDataSet->vParamType[ni] = Affy_GenericFile_LoadWString(fpIn, little_endian_machine);
		}
	}


	/* load column number */
	if(big_endian_fread(&(pDataSet->nColNum), DWORD_SIZE, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: Affy_GenericFile_LoadDataSet, cannot load column number!\n");
		Affy_GenericDataSet_Delete(&pDataSet);
		return NULL;
	}

	if(pDataSet->nColNum > 0)
	{
		pDataSet->vColName = (struct tagWString **)calloc(pDataSet->nColNum, sizeof(struct tagWString *));
		pDataSet->vColType = (unsigned char *)calloc(pDataSet->nColNum, sizeof(unsigned char));
		pDataSet->vColSize = (int *)calloc(pDataSet->nColNum, sizeof(int));

		if((pDataSet->vColName == NULL) || (pDataSet->vColType == NULL) || (pDataSet->vColSize == NULL))
		{
			printf("Error: Affy_GenericFile_LoadDataSet, cannot allocate memory for loading column name/type/size!\n");
			pDataSet->nColNum = 0;
			Affy_GenericDataSet_Delete(&pDataSet);
			return NULL;
		}

		/* load column parameters */
		for(ni=0; ni<(int)(pDataSet->nColNum); ni++)
		{
			pDataSet->vColName[ni] = Affy_GenericFile_LoadWString(fpIn, little_endian_machine);
			if(fread( pDataSet->vColType+ni, 1, 1, fpIn) != 1)
			{
				printf("Error: Affy_GenericFile_LoadDataSet, cannot allocate memory for loading column name/type/size!\n");
				Affy_GenericDataSet_Delete(&pDataSet);
				return NULL;
			}

			if(big_endian_fread(pDataSet->vColSize+ni, INT_SIZE, 1, fpIn, little_endian_machine) != 1)
			{
				printf("Error: Affy_GenericFile_LoadDataSet, cannot allocate memory for loading column name/type/size!\n");
				Affy_GenericDataSet_Delete(&pDataSet);
				return NULL;
			}
		}
	}

	/* load row number */
	if(big_endian_fread(&(pDataSet->nRowNum), DWORD_SIZE, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: Affy_GenericFile_LoadDataSet, cannot load row number!\n");
		Affy_GenericDataSet_Delete(&pDataSet);
		return NULL;
	}

	unitsize = 0;
	for(ni=0; ni<(int)(pDataSet->nColNum); ni++)
		unitsize += pDataSet->vColSize[ni];

	memsize = unitsize*(pDataSet->nRowNum);

	/* load data */
	if(memsize > 0)
	{
		pDataSet->vData = (char *)calloc(memsize, sizeof(char));
		if(pDataSet->vData == NULL)
		{
			printf("Error: Affy_GenericFile_LoadDataSet, cannot allocate memory for loading data!\n");
			Affy_GenericDataSet_Delete(&pDataSet);
			return NULL;
		}

		if(ftell(fpIn) != pDataSet->nFirstDataElementPos)
		{
			printf("Error: Affy_GenericFile_LoadDataSet, cannot match first data element position correctly!\n");
			Affy_GenericDataSet_Delete(&pDataSet);
			return NULL;
		}

		if( fread(pDataSet->vData, sizeof(char), memsize, fpIn) != memsize)
		{
			printf("Error: Affy_GenericFile_LoadDataSet, cannot load data correctly!\n");
			Affy_GenericDataSet_Delete(&pDataSet);
			return NULL;
		}

		/* transform data from big_endian to host byte order */
		if(little_endian_machine == 1)
		{
			buf = pDataSet->vData;
			for(ni=0; ni<(int)(pDataSet->nRowNum); ni++)
			{
				for(nj=0; nj<(int)(pDataSet->nColNum); nj++)
				{
					switch(pDataSet->vColType[nj])
					{
						/* short */
						case 2:
							reverse_buf(buf, pDataSet->vColSize[nj]);
							buf += pDataSet->vColSize[nj];
							break;
						/* ushort */
						case 3:
							reverse_buf(buf, pDataSet->vColSize[nj]);
							buf += pDataSet->vColSize[nj];
							break;
						/* int */
						case 4:
							reverse_buf(buf, pDataSet->vColSize[nj]);
							buf += pDataSet->vColSize[nj];
							break;
						/* uint */
						case 5:
							reverse_buf(buf, pDataSet->vColSize[nj]);
							buf += pDataSet->vColSize[nj];
							break;
						/* float */
						case 6:
							reverse_buf(buf, pDataSet->vColSize[nj]);
							buf += pDataSet->vColSize[nj];
							break;
						/* double */
						case 7:
							reverse_buf(buf, pDataSet->vColSize[nj]);
							buf += pDataSet->vColSize[nj];
							break;
						/* string */
						case 8:
							reverse_buf(buf, INT_SIZE);
							buf += pDataSet->vColSize[nj];
							break;
						/* wstring */
						case 9:
							reverse_buf(buf, INT_SIZE);
							memcpy(&nLen, buf, INT_SIZE);
							buf += INT_SIZE;
							for(nk=0; nk<nLen; nk++)
							{
								reverse_buf(buf, WCHAR16_SIZE);
								buf += WCHAR16_SIZE;
							}
							break;
						/* byte, ubyte */
						default:
							buf += pDataSet->vColSize[nj];
							break;
					}
				}
			}
		}
	}


	/* if(fread(chbuf, 1, 1, fpIn) != 1)
	{
		printf("Error: Affy_GenericFile_LoadDataSet, cannot load data correctly!\n");
		Affy_GenericDataSet_Delete(&pDataSet);
		return NULL;
	}*/

	/* return */
	return pDataSet;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_GenericFile_LoadWString()                                         */
/*  Load wchar_t string from generic file.                                 */
/* ----------------------------------------------------------------------- */ 
struct tagWString *Affy_GenericFile_LoadWString(FILE *fpIn, int little_endian_machine)
{
	/* define */
	int nLen = 0;
	struct tagWString *pString = NULL;
	char *buf = NULL;
	int ni;

	/* load */
	if(big_endian_fread(&nLen, INT_SIZE, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: Affy_GenericFile_LoadWString, cannot load wchar_t string length!\n");
		return NULL;
	}

	if(nLen < 0)
	{
		printf("Error: Affy_GenericFile_LoadWString, string length < 0!\n");
		return NULL;
	}

	pString = CreateWString(nLen);
	if(pString == NULL)
	{
		printf("Error: Affy_GenericFile_LoadWString, cannot allocate memory for loading wchar_t string!\n");
		return NULL;
	}

	buf = (char *)calloc(nLen+1, WCHAR16_SIZE);
	if(buf == NULL)
	{
		printf("Error: Affy_GenericFile_LoadWString, cannot create buffer for loading wchar_t string!\n");
		return NULL;
	}
	if(big_endian_fread(buf, WCHAR16_SIZE, nLen, fpIn, little_endian_machine) != nLen)
	{
		printf("Error: Affy_GenericFile_LoadWString, cannot load wchar_t string correctly!\n");
		DeleteWString(&pString);
		return NULL;
	}

	for(ni=0; ni<nLen; ni++)
	{
		pString->m_pWString[ni] = (wchar_t)(*((wchar_t16 *)(buf+ni*WCHAR16_SIZE)));
	}

	free(buf);

	/* for debug */
	/* printf("%S\n", pString->m_pWString); */

	/* return */
	return pString;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_GenericFile_LoadString()                                          */
/*  Load 1-byte character string from generic file.                        */
/* ----------------------------------------------------------------------- */ 
struct tagString *Affy_GenericFile_LoadString(FILE *fpIn, int little_endian_machine)
{
	/* define */
	int nLen = 0;
	struct tagString *pString = NULL;

	/* load */
	if(big_endian_fread(&nLen, INT_SIZE, 1, fpIn, little_endian_machine) != 1)
	{
		printf("Error: Affy_GenericFile_LoadString, cannot load 1 byte character string length!\n");
		return NULL;
	}

	if(nLen < 0)
	{
		printf("Error: Affy_GenericFile_LoadString, string length < 0!\n");
		return NULL;
	}

	pString = CreateString(nLen);
	if(pString == NULL)
	{
		printf("Error: Affy_GenericFile_LoadString, cannot allocate memory for loading 1 byte character string!\n");
		return NULL;
	}

	if(fread(pString->m_pString, 1, nLen, fpIn) != nLen)
	{
		printf("Error: Affy_GenericFile_LoadString, cannot load 1 byte character string correctly!\n");
		DeleteString(pString);
		pString = NULL;
		return NULL;
	}

	pString->m_pString[nLen] = '\0';

	/* for debug */
	/* printf("%s\n", pString->m_pString); */

	/* return */
	return pString;
}


/* ----------------------------------------------------------------------- */ 
/*  Affy_MIME_Value2Int()                                                  */
/*  Convert MIME Value to Int.                                             */
/* ----------------------------------------------------------------------- */ 
int Affy_MIME_Value2Int(char *strValue, int nBufferLen)
{
	/* define */
	int nVal = 0;
	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));
	char *buf = NULL;

	if(nBufferLen < DWORD_SIZE)
	{
		printf("Error: Affy_MIME_Value2Int, buffer length < DWORD_SIZE!\n");
		exit(EXIT_FAILURE);
	}

	buf = (char *)calloc(nBufferLen, sizeof(char));
	if(buf == NULL)
	{
		printf("Error: Affy_MIME_Value2Int, cannot create buffer!\n");
		exit(EXIT_FAILURE);
	}
	memcpy(buf, strValue, nBufferLen);
	if(little_endian_machine ==1)
		reverse_buf(buf, DWORD_SIZE);
	memcpy(&nVal, buf, DWORD_SIZE);
	
	free(buf);

	/* return */
	return nVal;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_MIME_Value2Float()                                                */
/*  Convert MIME Value to Float.                                           */
/* ----------------------------------------------------------------------- */ 
float Affy_MIME_Value2Float(char *strValue, int nBufferLen)
{
	/* define */
	float fVal = 0.0;	
	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));
	char *buf = NULL;

	if(nBufferLen < FLOAT_SIZE)
	{
		printf("Error: Affy_MIME_Value2Float, buffer length < FLOAT_SIZE!\n");
		exit(EXIT_FAILURE);
	}

	buf = (char *)calloc(nBufferLen, sizeof(char));
	if(buf == NULL)
	{
		printf("Error: Affy_MIME_Value2Float, cannot create buffer!\n");
		exit(EXIT_FAILURE);
	}
	memcpy(buf, strValue, nBufferLen);
	if(little_endian_machine ==1)
		reverse_buf(buf, FLOAT_SIZE);
	memcpy(&fVal, buf,FLOAT_SIZE);
	
	free(buf);
	
	/* return */
	return fVal;
}

/* ----------------------------------------------------------------------- */ 
/*  Affy_ParamValueType2String()                                           */
/*  Convert affymetrix Param/Value/Type to a 1 byte character string.      */
/* ----------------------------------------------------------------------- */ 
struct tagString *Affy_ParamValueType2String(wchar_t strParamName[], char strParamValue[], 
	int nBufferLen, wchar_t strParamType[])
{
	/* define */
	unsigned int x = 1;
	int little_endian_machine = (1 == *((char*)(&x)));
	char strTemp[LONG_LINE_LENGTH];
	struct tagString *pNewString = NULL;
	char *buf = NULL;
	int nTemp;
	float fTemp;
	wchar_t *wvbuf = NULL;
	char *buf2;
	char *vbuf;
	int nj;

	buf = (char *)calloc(nBufferLen, sizeof(char));
	if(buf == NULL)
	{
		printf("Warning: Affy_ParamValueType2String, empty buffer!\n");
		return NULL;
	}
	memcpy(buf, strParamValue, nBufferLen);

	sprintf(strTemp, "%S=", strParamName);
	StringAddTail(&pNewString, strTemp);

	if( wcscmp(strParamType, L"text/x-calvin-integer-8") ==0 )
	{
		nTemp = Affy_MIME_Value2Int(buf, nBufferLen);
		sprintf(strTemp, "%d", (char)nTemp);
	}
	else if( wcscmp(strParamType, L"text/x-calvin-unsigned-integer-8") ==0 )
	{
		nTemp = Affy_MIME_Value2Int(buf, nBufferLen);
		sprintf(strTemp, "%d", (unsigned char)nTemp);
	}
	else if( wcscmp(strParamType, L"text/x-calvin-integer-16") ==0 )
	{
		nTemp = Affy_MIME_Value2Int(buf, nBufferLen);
		sprintf(strTemp, "%d", (short)nTemp);
	}
	else if( wcscmp(strParamType, L"text/x-calvin-unsigned-integer-16") ==0 )
	{
		nTemp = Affy_MIME_Value2Int(buf, nBufferLen);
		sprintf(strTemp, "%d", (unsigned short)nTemp);
	}
	else if( wcscmp(strParamType, L"text/x-calvin-integer-32") ==0 )
	{
		nTemp = Affy_MIME_Value2Int(buf, nBufferLen);
		sprintf(strTemp, "%d", (int)nTemp);
	}
	else if( wcscmp(strParamType, L"text/x-calvin-unsigned-integer-32") ==0 )
	{
		nTemp = Affy_MIME_Value2Int(buf, nBufferLen);
		sprintf(strTemp, "%d", (unsigned int)nTemp);
	}
	else if( wcscmp(strParamType, L"text/x-calvin-float") ==0 )
	{
		fTemp = Affy_MIME_Value2Float(buf, nBufferLen);
		sprintf(strTemp, "%f", fTemp);
	}
	else if( wcscmp(strParamType, L"text/plain") ==0 )
	{
		wvbuf = (wchar_t *)calloc(nBufferLen/WCHAR16_SIZE+1, sizeof(wchar_t));
		buf2 = (char *)buf;
		for(nj=0; nj<nBufferLen/WCHAR16_SIZE; nj++)
		{
			if(little_endian_machine ==1)
				reverse_buf(buf2, WCHAR16_SIZE);
			wvbuf[nj] = (wchar_t)(*((wchar_t16 *)buf2));
			buf2 += WCHAR16_SIZE;
		}
		sprintf(strTemp, "%S", wvbuf);
		free(wvbuf);
	}
	else
	{
		vbuf = (char *)calloc(nBufferLen+1, sizeof(char));
		memcpy(vbuf, buf, nBufferLen);
		sprintf(strTemp, "%s", vbuf);
		free(vbuf);
	}

	free(buf);

	StringAddTail(&pNewString, strTemp);

	/* return */
	return pNewString;
}

