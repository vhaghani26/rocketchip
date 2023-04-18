/* ----------------------------------------------------------------------- */
/*  GenomeLib.c : implementation of the genome library                     */
/*  Author : Ji HongKai ; Time: 2004.08                                    */
/* ----------------------------------------------------------------------- */

#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
#include "limits.h"

#include "MathLib.h"
#include "MatrixLib.h"
#include "RandomLib.h"
#include "StringLib.h"
#include "SequenceLib.h"
#include "GenomeLib.h"
#include "MotifLib.h"
#include "AffyLib.h"
#include "PhysicalNetworkLib.h"

/* ----------------------------------------------------------------------- */ 
/*  RemoveFiles:                                                           */
/*  remove files.                                                          */
/* ----------------------------------------------------------------------- */ 
void RemoveFiles(char strPath[])
{
	/* define */
	char strCommand[LONG_LINE_LENGTH];

	/* remove */
	if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
	{
		sprintf(strCommand, "del %s", strPath);
		system(strCommand);
	}
	else
	{
		sprintf(strCommand, "rm %s", strPath);
		system(strCommand);
	}
}

/* ----------------------------------------------------------------------- */ 
/*  AdjustDirectoryPath:                                                   */
/*  add appropriate '/' or '\' to directory path based on OS system.       */
/* ----------------------------------------------------------------------- */ 
void AdjustDirectoryPath(char strPath[])
{
	/* define */
	int nlen;

	if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
	{
		nlen = (int)strlen(strPath);
		if(nlen == 0)
		{
			sprintf(strPath, ".\\"); 
		}
		else if(strPath[nlen-1] != '\\')
		{
			strPath[nlen] = '\\';
			strPath[nlen+1] = '\0';
		}
	}
	else
	{
		nlen = (int)strlen(strPath);
		if(nlen == 0)
		{
			sprintf(strPath, "./"); 
		}
		else if(strPath[nlen-1] != '/')
		{
			strPath[nlen] = '/';
			strPath[nlen+1] = '\0';
		}
	}
}

/* ----------------------------------------------------------------------- */ 
/*  Retrive file name:                                                     */
/*  Get file name from a file path.                                        */
/* ----------------------------------------------------------------------- */ 
int GetFileName(char strFilePath[], char strFileName[])
{
	/* define */
	int nLen;
	char *chp;

	/* operation */
	nLen = (int)(strlen(strFilePath));
	if(nLen == 0)
	{
		strcpy(strFileName, "");
		return PROC_SUCCESS;
	}

	if(strcmp(OS_SYSTEM, "WINDOWS") == 0)
		chp = strrchr(strFilePath, '\\');
	else
		chp = strrchr(strFilePath, '/');

	if(chp == NULL)
	{
		strcpy(strFileName, strFilePath);
		return PROC_SUCCESS;
	}

	if(chp == strFilePath+nLen-1)
	{
		strcpy(strFileName, "");
	}
	else
	{
		strcpy(strFileName, chp+1);
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  File_Bed2Cod_Main:                                                     */
/*  Convert a BED file to a COD file.                                      */
/* ----------------------------------------------------------------------- */ 
int File_Bed2Cod_Main(char strInputPath[], char strOutputPath[])
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	char strLine[LONG_LINE_LENGTH];
	char strAlias[LINE_LENGTH];
	char strChr[LINE_LENGTH];
	int nStart,nEnd;
	char chStrand;
	double dScore;
	int ni;
	char *chp1,*chp2;
	int nNoMoreField;

	/* open files */
	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: File_Bed2Cod_Main, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: File_Bed2Cod_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	fprintf(fpOut, "#seq_id\tchromosome\tstart\tend\tstrand\n");

	/* convert */
	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;
		if(strstr(strLine, "browser") == strLine)
			continue;
		if(strstr(strLine, "track") == strLine)
			continue;

		strcpy(strChr, "-1");
		nStart = -1;
		nEnd = -1;
		nNoMoreField = 0;
		sprintf(strAlias, "%d", ni+1);
		chStrand = '+';
		dScore = 0.0;

		chp1 = strLine;
		chp2 = strpbrk(chp1, WORD_SEPARATORS);
		if(chp2 == NULL)
			continue;
		*chp2 = '\0';
		strcpy(strChr, chp1);
		chp1 = chp2+1;

		chp2 = strpbrk(chp1, WORD_SEPARATORS);
		if(chp2 == NULL)
			continue;
		*chp2 = '\0';
		nStart = atoi(chp1);
		chp1 = chp2+1;

		chp2 = strpbrk(chp1, WORD_SEPARATORS);
		if(chp2 == NULL)
		{
			nEnd = atoi(chp1);
			nNoMoreField = 1;
		}
		else
		{
			*chp2 = '\0';
			nEnd = atoi(chp1);
			chp1 = chp2+1;
		}
		
		if(nNoMoreField == 0)
		{
			chp2 = strpbrk(chp1, WORD_SEPARATORS);
			if(chp2 == NULL)
			{
				strcpy(strAlias, chp1);
				nNoMoreField = 1;
			}
			else
			{
				*chp2 = '\0';
				strcpy(strAlias, chp1);
				chp1 = chp2+1;
			}
		}

		if(nNoMoreField == 0)
		{
			chp2 = strpbrk(chp1, WORD_SEPARATORS);
			if(chp2 == NULL)
			{
				dScore = atof(chp1);
				nNoMoreField = 1;
			}
			else
			{
				*chp2 = '\0';
				dScore = atof(chp1);
				chp1 = chp2+1;
			}
		}

		if(nNoMoreField == 0)
		{
			StrTrimLeft(chp1);
			chp2 = strpbrk(chp1, WORD_SEPARATORS);
			if(chp2 == NULL)
			{
				chStrand = *chp1;
				nNoMoreField = 1;
			}
			else
			{
				*chp2 = '\0';
				chStrand = *chp1;
				chp1 = chp2+1;
			}
		}

		fprintf(fpOut, "%s\t%s\t%d\t%d\t%c\n", strAlias, strChr, nStart, nEnd, chStrand);
		ni++;
	}


	/* close files */
	fclose(fpIn);
	fclose(fpOut);
	
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  File_Cod2Bed_Main:                                                     */
/*  Convert a COD file to a BED file.                                      */
/* ----------------------------------------------------------------------- */ 
int File_Cod2Bed_Main(char strInputPath[], char strOutputPath[])
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	char strFileName[MED_LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	char strAlias[LINE_LENGTH];
	char strChr[LINE_LENGTH];
	int nStart,nEnd;
	char chStrand;
	int ni;

	/* open files */
	GetFileName(strInputPath, strFileName);

	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: File_Cod2Bed_Main, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: File_Cod2Bed_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	/* convert */
	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		sscanf(strLine, "%s %s %d %d %c", strAlias, strChr, &nStart, &nEnd, &chStrand);

		if(ni == 0)
		{
			fprintf(fpOut, "browser position %s:%d-%d\n", strChr, nStart, nEnd);
			fprintf(fpOut, "track name=\"%s\" description=\"%s\"\n", strFileName, strFileName);
		}

		fprintf(fpOut, "%s\t%d\t%d\t%s\t1000\t%c\n", strChr, nStart, nEnd, strAlias, chStrand);
		ni++;
	}


	/* close files */
	fclose(fpIn);
	fclose(fpOut);
	
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  File_Bed2Aln_Main:                                                     */
/*  Convert a BED file to a ALN file.                                      */
/* ----------------------------------------------------------------------- */ 
int File_Bed2Aln_Main(char strInputPath[], char strOutputPath[])
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	char strLine[LONG_LINE_LENGTH];
	char strAlias[LINE_LENGTH];
	char strChr[LINE_LENGTH];
	int nStart,nEnd;
	char chStrand;
	double dScore;
	int ni;
	char *chp1,*chp2;
	int nNoMoreField;

	/* open files */
	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: File_Bed2Aln_Main, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: File_Bed2Aln_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	/* convert */
	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;
		if(strstr(strLine, "browser") == strLine)
			continue;
		if(strstr(strLine, "track") == strLine)
			continue;

		strcpy(strChr, "-1");
		nStart = -1;
		nEnd = -1;
		nNoMoreField = 0;
		sprintf(strAlias, "%d", ni+1);
		chStrand = '+';
		dScore = 0.0;

		chp1 = strLine;
		chp2 = strpbrk(chp1, WORD_SEPARATORS);
		if(chp2 == NULL)
			continue;
		*chp2 = '\0';
		strcpy(strChr, chp1);
		chp1 = chp2+1;

		chp2 = strpbrk(chp1, WORD_SEPARATORS);
		if(chp2 == NULL)
			continue;
		*chp2 = '\0';
		nStart = atoi(chp1);
		chp1 = chp2+1;

		chp2 = strpbrk(chp1, WORD_SEPARATORS);
		if(chp2 == NULL)
		{
			nEnd = atoi(chp1);
			nNoMoreField = 1;
		}
		else
		{
			*chp2 = '\0';
			nEnd = atoi(chp1);
			chp1 = chp2+1;
		}
		
		if(nNoMoreField == 0)
		{
			chp2 = strpbrk(chp1, WORD_SEPARATORS);
			if(chp2 == NULL)
			{
				strcpy(strAlias, chp1);
				nNoMoreField = 1;
			}
			else
			{
				*chp2 = '\0';
				strcpy(strAlias, chp1);
				chp1 = chp2+1;
			}
		}

		if(nNoMoreField == 0)
		{
			chp2 = strpbrk(chp1, WORD_SEPARATORS);
			if(chp2 == NULL)
			{
				dScore = atof(chp1);
				nNoMoreField = 1;
			}
			else
			{
				*chp2 = '\0';
				dScore = atof(chp1);
				chp1 = chp2+1;
			}
		}

		if(nNoMoreField == 0)
		{
			StrTrimLeft(chp1);
			chp2 = strpbrk(chp1, WORD_SEPARATORS);
			if(chp2 == NULL)
			{
				chStrand = *chp1;
				nNoMoreField = 1;
			}
			else
			{
				*chp2 = '\0';
				chStrand = *chp1;
				chp1 = chp2+1;
			}
		}

		fprintf(fpOut, "%s\t%d\t%c\n", strChr, (nStart+nEnd)/2, chStrand);
		ni++;
	}


	/* close files */
	fclose(fpIn);
	fclose(fpOut);
	
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Genome_Fasta_To_Code_4bit_Main:                                        */
/*  Convert Genome Sequence from Fasta files to coding files.              */
/*  1 byte for 2 bases.                                                    */
/*  0000: A; 0001: C; 0010: G; 0011: T.                                    */
/*  0100: a; 0101: c; 0110: g; 0111: t.                                    */
/*  1000: N;                                                               */ 
/*  return number of fasta files coded.                                    */
/* ----------------------------------------------------------------------- */ 
int Genome_Fasta_To_Code_4bit_Main(char strInPath[], char strOutPath[])
{
	/* define */
	int nCount;
	FILE *fpIn;
	FILE *fpChrLen;
	char strLine[LINE_LENGTH];
	char strInFile[LINE_LENGTH];
	char strOutFile[LINE_LENGTH];
	int nLen;

	/* init */
	nCount = 0;
	AdjustDirectoryPath(strInPath);
	AdjustDirectoryPath(strOutPath);

	sprintf(strInFile, "%schrlist.txt", strInPath);
	fpIn = NULL;
	fpIn = fopen(strInFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: cannot open chromosome list file!\n");
		return 0;
	}

	sprintf(strOutFile, "%schrlen.txt", strOutPath);
	fpChrLen = fopen(strOutFile, "w");
	if(fpChrLen == NULL)
	{
		printf("Error: cannot open chromosome length file!\n");
		return 0;
	}

	/* convert */
	while(fgets(strLine, LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sprintf(strInFile, "%s%s.fa", strInPath, strLine);
		sprintf(strOutFile, "%s%s.sq", strOutPath, strLine);
		nLen = Genome_Fasta_To_Code_4bit(strInFile, strOutFile);
		fprintf(fpChrLen, "%d\n", nLen);
	}
	
	/* close file */
	fclose(fpIn);
	fclose(fpChrLen);

	/* return */
	return nCount;
}

/* ----------------------------------------------------------------------- */ 
/*  Genome_Fasta_To_Code_4bit:                                             */
/*  Convert Genome Sequence from Fasta files to coding files.              */
/*  1 byte for 2 bases.                                                    */
/*  0000: A; 0001: C; 0010: G; 0011: T.                                    */
/*  0100: a; 0101: c; 0110: g; 0111: t.                                    */
/*  1000: N;                                                               */ 
/* ----------------------------------------------------------------------- */ 
int Genome_Fasta_To_Code_4bit(char strInFile[], char strOutFile[])
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	char strLine[LINE_LENGTH];
	unsigned char vBase[LINE_LENGTH];
	unsigned char bCode;
	int nLen,ni,nj,nk;
	int numwritten;
	int nTotalLen;
	int nEmptyFile;

	/* open file */
	nEmptyFile = 0;
	fpIn = NULL;
	fpIn = fopen(strInFile, "r");
	if(fpIn == NULL)
	{
		printf("Warning: Genome_Fasta_To_Code_4bit, cannot open Fasta file!\n");
		nEmptyFile = 1;
	}

	fpOut = NULL;
	fpOut = fopen(strOutFile, "wb");
	if(fpOut == NULL)
	{
		if(nEmptyFile == 0)
		{
			fclose(fpIn);
		}
		printf("Error: Genome_Fasta_To_Code_4bit, cannot open Code file!\n");
		exit(EXIT_FAILURE);
	}

	/* if the file is empty */
	if(nEmptyFile == 1)
	{
		bCode = 8;
		vBase[0] = bCode << 4;
		vBase[0] = (vBase[0] & 0xF0);
		vBase[0] = ( vBase[0] | bCode );
		numwritten = fwrite(vBase , sizeof(unsigned char), 1, fpOut);
		nTotalLen = 2*numwritten;
		fclose(fpOut);
		return nTotalLen;
	}

	/* read and convert */
	nj = 0;
	nk = 0;
	nTotalLen = 0;
	while(fgets(strLine, LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '>')
			continue;

		/* get length */
		nLen = strlen(strLine);
		for(ni=0; ni<nLen; ni++)
		{
			if(nj == 0)
			{
				switch(strLine[ni])
				{
					case 'A': bCode = 0;
						break;
					case 'a': bCode = 4;
						break;
					case 'C': bCode = 1;
						break;
					case 'c': bCode = 5;
						break;
					case 'G': bCode = 2;
						break;
					case 'g': bCode = 6;
						break;
					case 'T': bCode = 3;
						break;
					case 't': bCode = 7;
						break;
					default: bCode = 8;
				}
				vBase[nk] = bCode << 4;
				vBase[nk] = (vBase[nk] & 0xF0);
				nj = 1;	
			}
			else
			{
				switch(strLine[ni])
				{
					case 'A': bCode = 0;
						break;
					case 'a': bCode = 4;
						break;
					case 'C': bCode = 1;
						break;
					case 'c': bCode = 5;
						break;
					case 'G': bCode = 2;
						break;
					case 'g': bCode = 6;
						break;
					case 'T': bCode = 3;
						break;
					case 't': bCode = 7;
						break;
					default: bCode = 8;
				}
				vBase[nk] = ( vBase[nk] | bCode );
				nj = 0;
				nk++;
				if(nk == CODE_BATCH_LEN)
				{
					numwritten = fwrite(vBase , sizeof(unsigned char), CODE_BATCH_LEN, fpOut);
					nTotalLen += 2*numwritten;
					nk = 0;
				}
			}
		}
	}

	if(nj == 0)
	{
		if(nk > 0)
		{
			numwritten = fwrite(vBase , sizeof(unsigned char), nk, fpOut);
			nTotalLen += 2*numwritten;
		}
	}
	else
	{
		numwritten = fwrite(vBase , sizeof(unsigned char), (nk+1), fpOut);
		nTotalLen += 2*numwritten-1;
	}


	/* close file */
	fclose(fpIn);
	fclose(fpOut);

	/* return */
	return nTotalLen;
}

/* ----------------------------------------------------------------------- */ 
/*  Genome_PhastCons_To_Code_8bit_Main:                                    */
/*  Convert PhastCons score to genomelab coded score.                      */
/*  1 byte per base.                                                       */
/*  255 highest score.                                                     */
/*  return length of chromosome coded.                                     */
/* ----------------------------------------------------------------------- */ 
int Genome_PhastCons_To_Code_8bit_Main(char strGenomePath[], char strPhastPath[], 
									   char strOutPath[], char strExt[])
{
	/* define */
	int nCount;
	FILE *fpIn;
	struct INTMATRIX *pChrLen;

	char strLine[LINE_LENGTH];
	char strInFile[LINE_LENGTH];
	char strPhastFile[LINE_LENGTH];
	char strOutFile[LINE_LENGTH];
	int nChrLen,nLen;

	AdjustDirectoryPath(strGenomePath);
	AdjustDirectoryPath(strPhastPath);
	AdjustDirectoryPath(strOutPath);

	/* init */
	sprintf(strInFile, "%schrlen.txt", strGenomePath);
	pChrLen = NULL;
	pChrLen = IMLOAD(strInFile);
	if(pChrLen == NULL)
	{
		printf("Error: cannot open chromosome length file!\n");
		return 0;
	}

	sprintf(strInFile, "%schrlist.txt", strGenomePath);
	fpIn = NULL;
	fpIn = fopen(strInFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: cannot open chromosome list file!\n");
		return 0;
	}

	/* convert */
	nCount = 0;
	while(fgets(strLine, LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		nChrLen = pChrLen->pMatElement[nCount];
		sprintf(strPhastFile, "%s%s%s", strPhastPath, strLine, strExt);
		sprintf(strOutFile, "%s%s.cs", strOutPath, strLine);
		nLen = Genome_PhastCons_To_Code_8bit(strPhastFile, strOutFile, nChrLen);
		printf("%s %d %d\n", strLine, nLen, nChrLen);
		nCount++;
	}
	
	/* close file */
	fclose(fpIn);
	DestroyIntMatrix(pChrLen);

	/* return */
	return nCount;
}

/* ----------------------------------------------------------------------- */ 
/*  Genome_PhastCons_To_Code_8bit:                                         */
/*  Convert PhastCons score to genomelab coded score.                      */
/*  1 byte per base.                                                       */
/*  255 highest score.                                                     */
/*  return length of chromosome coded.                                     */
/* ----------------------------------------------------------------------- */ 
int Genome_PhastCons_To_Code_8bit(char strPhastFile[], char strOutFile[], int nChrLen)
{
	/* define */
	int ni;
	FILE *fpPhast;
	FILE *fpCS;
	int numwritten;
	char strLine[LINE_LENGTH];
	int nPos;
	double dScore;
	unsigned char bScore;
	unsigned char zScore;
	int nEmptyFile;

	/* init */
	nEmptyFile = 0;
	fpPhast = NULL;
	fpPhast = fopen(strPhastFile, "r");
	if(fpPhast == NULL)
	{
		printf("Warning: cannot open phastCons file!\n");
		nEmptyFile = 1;
	}

	fpCS = NULL;
	fpCS = fopen(strOutFile, "w+b");
	if(fpCS == NULL)
	{
		printf("Error: cannot open .cs file!\n");
		exit(EXIT_FAILURE);
	}

	if(nEmptyFile == 1)
	{
		zScore = 0;
		for(ni=0; ni<nChrLen; ni++)
		{
			numwritten = fwrite(&zScore, sizeof(unsigned char), 1, fpCS);
		}

		fclose(fpCS);
	}

	/* load file */
	zScore = 0;
	ni = 0;
	while(fgets(strLine, LINE_LENGTH, fpPhast) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%d %lf", &nPos, &dScore);
		nPos--;
		bScore = (unsigned char)(dScore*255);
		while(ni<nPos)
		{
			numwritten = fwrite(&zScore, sizeof(unsigned char), 1, fpCS);
			ni++;
		}

		numwritten = fwrite(&bScore, sizeof(unsigned char), 1, fpCS);
		ni++;
	}

	while(ni<nChrLen)
	{
		numwritten = fwrite(&zScore, sizeof(unsigned char), 1, fpCS);
		ni++;
	}

	if(ni != nChrLen)
	{
		printf("Error: coordinate wrong!\n");
		exit(EXIT_FAILURE);
	}

	/* close files */
	fclose(fpPhast);
	fclose(fpCS);

	/* return */
	return ni;
}

/* ----------------------------------------------------------------------- */ 
/*  Genome_PhastCons_To_Code_8bit_Main_v2:                                 */
/*  To deal with phastcons Format released Nov. 2004 and later.            */ 
/*  Convert PhastCons score to genomelab coded score.                      */
/*  1 byte per base.                                                       */
/*  255 highest score.                                                     */
/*  return length of chromosome coded.                                     */
/* ----------------------------------------------------------------------- */ 
int Genome_PhastCons_To_Code_8bit_Main_v2(char strGenomePath[], char strPhastPath[], 
									   char strOutPath[], char strExt[])
{
	/* define */
	int nCount;
	FILE *fpIn;
	struct INTMATRIX *pChrLen;

	char strLine[LINE_LENGTH];
	char strInFile[LINE_LENGTH];
	char strPhastFile[LINE_LENGTH];
	char strOutFile[LINE_LENGTH];
	int nChrLen,nLen;

	/* init */
	AdjustDirectoryPath(strGenomePath);
	AdjustDirectoryPath(strPhastPath);
	AdjustDirectoryPath(strOutPath);

	sprintf(strInFile, "%schrlen.txt", strGenomePath);
	pChrLen = NULL;
	pChrLen = IMLOAD(strInFile);
	if(pChrLen == NULL)
	{
		printf("Error: cannot open chromosome length file!\n");
		return 0;
	}

	sprintf(strInFile, "%schrlist.txt", strGenomePath);
	fpIn = NULL;
	fpIn = fopen(strInFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: cannot open chromosome list file!\n");
		return 0;
	}

	/* convert */
	nCount = 0;
	while(fgets(strLine, LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		nChrLen = pChrLen->pMatElement[nCount];
		sprintf(strPhastFile, "%s%s%s", strPhastPath, strLine, strExt);
		sprintf(strOutFile, "%s%s.cs", strOutPath, strLine);
		nLen = Genome_PhastCons_To_Code_8bit_v2(strPhastFile, strOutFile, nChrLen);
		printf("%s %d %d\n", strLine, nLen, nChrLen);
		nCount++;
	}
	
	/* close file */
	fclose(fpIn);
	DestroyIntMatrix(pChrLen);

	/* return */
	return nCount;
}


/* ----------------------------------------------------------------------- */ 
/*  Genome_PhastCons_To_Code_8bit_v2:                                      */
/*  Convert PhastCons score to genomelab coded score.                      */
/*  1 byte per base.                                                       */
/*  255 highest score.                                                     */
/*  return length of chromosome coded.                                     */
/* ----------------------------------------------------------------------- */ 
int Genome_PhastCons_To_Code_8bit_v2(char strPhastFile[], char strOutFile[], int nChrLen)
{
	/* define */
	int ni,nj;
	FILE *fpPhast;
	FILE *fpCS;
	int numwritten;
	char strLine[LINE_LENGTH];
	int nPos;
	double dScore;
	unsigned char bScore;
	unsigned char zScore;
	char strHead[LINE_LENGTH];
	char strChrom[LINE_LENGTH];
	char strStart[LINE_LENGTH];
	char strStep[LINE_LENGTH];
	int nStart,nStep;
	int nEmptyFile;

	/* init */
	nEmptyFile = 0;
	fpPhast = NULL;
	fpPhast = fopen(strPhastFile, "r");
	if(fpPhast == NULL)
	{
		printf("Warning: cannot open phastCons file!\n");
		nEmptyFile = 1;
	}

	fpCS = NULL;
	fpCS = fopen(strOutFile, "w+b");
	if(fpCS == NULL)
	{
		printf("Error: cannot open .cs file!\n");
		exit(EXIT_FAILURE);
	}

	if(nEmptyFile == 1)
	{
		zScore = 0;
		for(ni=0; ni<nChrLen; ni++)
		{
			numwritten = fwrite(&zScore, sizeof(unsigned char), 1, fpCS);
		}

		fclose(fpCS);
		return ni;
	}

	/* load file */
	zScore = 0;
	ni = 0;
	while(fgets(strLine, LINE_LENGTH, fpPhast) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		/* deal with head info */
		if(strstr(strLine, "fixedStep") == strLine)
		{
			sscanf(strLine, "%s %s %s %s", strHead, strChrom, strStart, strStep);

			if(strstr(strStart, "start=") == strStart)
			{
				nStart = atoi(strStart+6);
			}
			else
			{
				printf("error: Genome_PhastCons_To_Code_8bit_v2, cannot read head info correctly!\n");
			}

			if(strstr(strStep, "step=") == strStep)
			{
				nStep = atoi(strStep+5);
			}
			else
			{
				printf("error: Genome_PhastCons_To_Code_8bit_v2, cannot read head info correctly!\n");
			}

			nPos = nStart-1;
			while(ni<nPos)
			{
				numwritten = fwrite(&zScore, sizeof(unsigned char), 1, fpCS);
				ni++;
			}
		}

		else
		{
			sscanf(strLine, "%lf", &dScore);
			bScore = (unsigned char)(dScore*255);
			for(nj=0; nj<nStep; nj++)
			{
				numwritten = fwrite(&bScore, sizeof(unsigned char), 1, fpCS);
				ni++;
			}
		}
	}

	while(ni<nChrLen)
	{
		numwritten = fwrite(&zScore, sizeof(unsigned char), 1, fpCS);
		ni++;
	}

	if(ni != nChrLen)
	{
		printf("Error: coordinate wrong!\n");
		exit(EXIT_FAILURE);
	}

	/* close files */
	fclose(fpPhast);
	fclose(fpCS);

	/* return */
	return ni;
}

/* ----------------------------------------------------------------------- */ 
/*  Genome_CodeCDS_FromRefGene_Main:                                       */
/*  Code CDS regions.                                                      */
/*  1 byte for 1 base.                                                     */
/*  0: intergenic.                                                         */
/*  1: CDS                                                                 */
/* ----------------------------------------------------------------------- */ 
int Genome_CodeCDS_FromRefGene_Main(char strGenomePath[], char strRefGenePath[], 
					int nGType, char strSpecies[], char strOutputPath[])
{
	/* define */
	int nCount;
	FILE *fpIn;
	FILE *fpCDS;
	struct INTMATRIX *pChrLen;
	struct tagRefGene **vRefGeneDatabase;
	int nRefGeneNum;

	char strLine[LINE_LENGTH];
	char strInFile[LINE_LENGTH];
	char strOutFile[LINE_LENGTH];
	int nChrLen;
	int ni,nj,nk,nStart,nEnd,numwritten;
	unsigned char zScore;

	/* init */
	/* AdjustDirectoryPath(strGenomePath);
	AdjustDirectoryPath(strOutputPath); */
	sprintf(strInFile, "%schrlen.txt", strGenomePath);
	pChrLen = NULL;
	pChrLen = IMLOAD(strInFile);
	if(pChrLen == NULL)
	{
		printf("Error: Genome_CodeCDS_FromRefGene_Main, cannot open chromosome length file!\n");
		return 0;
	}

	sprintf(strInFile, "%schrlist.txt", strGenomePath);
	fpIn = NULL;
	fpIn = fopen(strInFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: Genome_CodeCDS_FromRefGene_Main, cannot open chromosome list file!\n");
		return 0;
	}

	/* convert */
	nCount = 0;
	while(fgets(strLine, LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		nChrLen = pChrLen->pMatElement[nCount];
		
		sprintf(strOutFile, "%s%s.cds", strOutputPath, strLine);
		fpCDS = NULL;
		fpCDS = fopen(strOutFile, "wb");
		if(fpCDS == NULL)
		{
			printf("Error: Genome_CodeCDS_FromRefGene_Main, cannot open .cs file!\n");
			exit(EXIT_FAILURE);
		}

		zScore = 0;
		for(ni=0; ni<nChrLen; ni++)
		{
			numwritten = fwrite(&zScore, sizeof(unsigned char), 1, fpCDS);
		}

		fclose(fpCDS);
		
		nCount++;
	}
	
	/* close file */
	fclose(fpIn);

	/* code */
	vRefGeneDatabase = RefGene_LoadDatabase(strRefGenePath, 
					nGType, strSpecies, &nRefGeneNum);

	/* process one by one */
	for(ni=0; ni<nRefGeneNum; ni++)
	{
		sprintf(strOutFile, "%s%s.cds", strOutputPath, vRefGeneDatabase[ni]->strChrom);
		fpCDS = NULL;
		fpCDS = fopen(strOutFile, "r+b");
		if(fpCDS == NULL)
		{
			printf("Error: Genome_CodeCDS_FromRefGene_Main, cannot open .cds file!\n");
			exit(EXIT_FAILURE);
		}

		/* process exon by exon */
		for(nj=0; nj<vRefGeneDatabase[ni]->nExonCount; nj++)
		{
			nStart = IMGETAT(vRefGeneDatabase[ni]->pmatExonStartsEnds, nj, 0);
			nEnd = IMGETAT(vRefGeneDatabase[ni]->pmatExonStartsEnds, nj, 1);
			if(nStart < vRefGeneDatabase[ni]->nCdsStart)
				nStart = vRefGeneDatabase[ni]->nCdsStart;
			if(nEnd > vRefGeneDatabase[ni]->nCdsEnd)
				nEnd = vRefGeneDatabase[ni]->nCdsEnd;

			if(nStart <= nEnd)
			{
				if( fseek( fpCDS, nStart, SEEK_SET ) != 0 )
				{
					printf("Error: Genome_CodeCDS_FromRefGene_Main, cannot locate the required sequence!\n");
					exit(EXIT_FAILURE);
				}

				zScore = 1;
				for(nk=nStart; nk<=nEnd; nk++)
				{
					numwritten = fwrite(&zScore, sizeof(unsigned char), 1, fpCDS);
				}
			}
		}

		fclose(fpCDS);
	}

	/* release memory */
	RefGene_ClearDatabase(&vRefGeneDatabase, nRefGeneNum);
	DestroyIntMatrix(pChrLen);

	/* return */
	return nCount;
}

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
								 int nStrandType)
{
	/* define */
	int nCount;
	char strSeqFile[LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	char strChrName[LINE_LENGTH];
	char strAlias[LINE_LENGTH];
	FILE *fpIn;
	FILE *fpOut;

	int nChr,nStart,nEnd;
	char chStrand;
	struct tagSequence *pSeq;
	struct INTMATRIX *pChrLen;

	/* init */
	nCount = 0;
	AdjustDirectoryPath(strGenomePath);

	/* load chrmosome length */
	sprintf(strLine, "%schrlen.txt", strGenomePath);
	pChrLen = NULL;
	pChrLen = IMLOAD(strLine);
	if(pChrLen == NULL)
	{
		printf("Error: cannot load chromosome length!\n");
		exit(EXIT_FAILURE);
	}

	/* open files */
	fpIn = NULL;
	fpIn = fopen(strInFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: Genome_Code_4bit_GetSeq_Main, cannot open input file.\n");
		return 0;
	}

	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: Genome_Code_4bit_GetSeq_Main, cannot open output file.\n");
		return 0;
	}

	/* get seq */
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;
		
		sscanf(strLine, "%s %d %d %d %c", strAlias, &nChr, &nStart, &nEnd, &chStrand);
		Genome_Index_To_ChromosomeName(strChrName, strSpecies, nChr);

		if( (nChr <= 0) || (nChr > pChrLen->nHeight))
		{
			printf("Warning: (%d:%d-%d) is out of normal genome range, sequence not obtained!\n", nChr, nStart, nEnd);
			continue;
		}

		if( (nStart < 0) || (nEnd >= pChrLen->pMatElement[nChr-1]) )
		{
			printf("Warning: (%d:%d-%d) is out of normal genome range, sequence not obtained!\n", nChr, nStart, nEnd);
			continue;
		}

		sprintf(strSeqFile, "%s%s.sq", strGenomePath, strChrName);
		pSeq = NULL;
		pSeq = Genome_Code_4bit_GetSeq(strSeqFile, nStart, nEnd);
		if(pSeq == NULL)
			continue;

		pSeq->m_nIndex = nCount;
		strcpy(pSeq->m_strAlias, strAlias);
		if(nStrandType == 1)
			SequenceWriteToFasta_ByStrand(pSeq, fpOut, chStrand, 1);
		else
			SequenceWriteToFasta_ByStrand(pSeq, fpOut, '+', 1);
		SequenceDelete(pSeq);
		nCount++;
	}

	/* close file */
	fclose(fpIn);
	fclose(fpOut);
	DestroyIntMatrix(pChrLen);

	/* return */
	return nCount;
}


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
								 int nStrandType)
{
	/* define */
	int nCount;
	char strSeqFile[LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	char strChrName[LINE_LENGTH];
	char strAlias[LINE_LENGTH];
	FILE *fpIn;
	FILE *fpOut;

	int nChr,nStart,nEnd;
	char chStrand;
	struct tagSequence *pSeq;
	struct INTMATRIX *pChrLen;

	/* init */
	nCount = 0;
	AdjustDirectoryPath(strGenomePath);

	/* load chrmosome length */
	sprintf(strLine, "%schrlen.txt", strGenomePath);
	pChrLen = NULL;
	pChrLen = IMLOAD(strLine);
	if(pChrLen == NULL)
	{
		printf("Error: cannot load chromosome length!\n");
		exit(EXIT_FAILURE);
	}

	/* open files */
	fpIn = NULL;
	fpIn = fopen(strInFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: Genome_Code_4bit_GetSeq_Main, cannot open input file.\n");
		return 0;
	}

	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: Genome_Code_4bit_GetSeq_Main, cannot open output file.\n");
		return 0;
	}

	/* get seq */
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;
		
		sscanf(strLine, "%s %s %d %d %c", strAlias, strChrName, &nStart, &nEnd, &chStrand);
		nChr = Genome_ChromosomeName_To_Index(strChrName, strSpecies);
		
		if( (nChr <= 0) || (nChr > pChrLen->nHeight))
		{
			printf("Warning: (%s:%d-%d) is out of normal genome range, sequence not obtained!\n", strChrName, nStart, nEnd);
			continue;
		}

		if( (nStart < 0) || (nEnd >= pChrLen->pMatElement[nChr-1]) )
		{
			printf("Warning: (%s:%d-%d) is out of normal genome range, sequence not obtained!\n", strChrName, nStart, nEnd);
			continue;
		}

		sprintf(strSeqFile, "%s%s.sq", strGenomePath, strChrName);
		pSeq = NULL;
		pSeq = Genome_Code_4bit_GetSeq(strSeqFile, nStart, nEnd);
		if(pSeq == NULL)
			continue;

		pSeq->m_nIndex = nCount;
		strcpy(pSeq->m_strAlias, strAlias);
		if(nStrandType == 1)
			SequenceWriteToFasta_ByStrand(pSeq, fpOut, chStrand, 1);
		else
			SequenceWriteToFasta_ByStrand(pSeq, fpOut, '+', 1);
		SequenceDelete(pSeq);
		nCount++;
	}

	/* close file */
	fclose(fpIn);
	fclose(fpOut);
	DestroyIntMatrix(pChrLen);

	/* return */
	return nCount;
}


/* ----------------------------------------------------------------------- */ 
/*  Genome_Code_4bit_GetSeq:                                               */
/*  Get sequences from Genome Sequence coding files.                       */
/*  1 byte for 2 bases.                                                    */
/*  0000: A; 0001: C; 0010: G; 0011: T.                                    */
/*  0100: a; 0101: c; 0110: g; 0111: t.                                    */
/*  1000: N;                                                               */ 
/* ----------------------------------------------------------------------- */ 
struct tagSequence *Genome_Code_4bit_GetSeq(char strInFile[], int nStart, int nEnd)
{
	/* define */
	struct tagSequence *pSeq;
	FILE *fpIn;
	int nP1,nP2,nR1,nR2;
	int numread;

	char strLine[LINE_LENGTH];
	int ni,nk,nl,nLen;
	unsigned char bBase,bChar;

	/* init */
	fpIn = NULL;
	fpIn = fopen(strInFile, "rb");
	if(fpIn == NULL)
	{
		printf("Error: Genome_Code_4bit_GetSeq, cannot open sequence file!\n");
		exit(EXIT_FAILURE);
	}
	if(nStart > nEnd)
	{
		printf("Error: Genome_Code_4bit_GetSeq, start > end!\n");
		exit(EXIT_FAILURE);
	}
	nLen = nEnd-nStart+1;

	pSeq = NULL;
	pSeq = SequenceCreate();
	if(pSeq == NULL)
	{
		printf("Error: Genome_Code_4bit_GetSeq, cannot create sequences.\n");
		exit(EXIT_FAILURE);
	}

	/* load seq */
	nP1 = (int)(nStart/2);
	nR1 = nStart%2;
	nP2 = (int)(nEnd/2);
	nR2 = nEnd%2;

	if( fseek( fpIn, nP1, SEEK_SET ) != 0 )
	{
		printf("Error: Genome_Code_4bit_GetSeq, cannot locate the required sequence!\n");
		exit(EXIT_FAILURE);
	}

	nk = 0;
	nl = 0;

	/* first base */
	numread = fread( &bChar, sizeof( unsigned char ), 1, fpIn );
	if(nR1 == 0)
	{
		bBase = bChar >> 4;
		switch(bBase)
		{
			case 0: strLine[nk] = 'A';
				break;
			case 1: strLine[nk] = 'C';
				break;
			case 2: strLine[nk] = 'G';
				break;
			case 3: strLine[nk] = 'T';
				break;
			case 4: strLine[nk] = 'a';
				break;
			case 5: strLine[nk] = 'c';
				break;
			case 6: strLine[nk] = 'g';
				break;
			case 7: strLine[nk] = 't';
				break;
			default: strLine[nk] = 'N';
		}
		nk++;
		nl++;
	}

	if(nl < nLen)
	{
		bBase = bChar & 0x0F;
		switch(bBase)
		{
			case 0: strLine[nk] = 'A';
				break;
			case 1: strLine[nk] = 'C';
				break;
			case 2: strLine[nk] = 'G';
				break;
			case 3: strLine[nk] = 'T';
				break;
			case 4: strLine[nk] = 'a';
				break;
			case 5: strLine[nk] = 'c';
				break;
			case 6: strLine[nk] = 'g';
				break;
			case 7: strLine[nk] = 't';
				break;
			default: strLine[nk] = 'N';
		}
		nk++;
		nl++;
	}

	/* middle bases */
	if(nl < nLen)
	{
		for(ni=(nP1+1); ni<nP2; ni++)
		{
			numread = fread( &bChar, sizeof( unsigned char ), 1, fpIn );
			bBase = bChar >> 4;
			switch(bBase)
			{
				case 0: strLine[nk] = 'A';
					break;
				case 1: strLine[nk] = 'C';
					break;
				case 2: strLine[nk] = 'G';
					break;
				case 3: strLine[nk] = 'T';
					break;
				case 4: strLine[nk] = 'a';
					break;
				case 5: strLine[nk] = 'c';
					break;
				case 6: strLine[nk] = 'g';
					break;
				case 7: strLine[nk] = 't';
					break;
				default: strLine[nk] = 'N';
			}
			nk++;
			nl++;

			if(nk == CODE_BATCH_LEN)
			{
				strLine[nk] = '\0';
				SequenceAddTail(pSeq, strLine);
				nk = 0;
			}

			bBase = bChar & 0x0F;
			switch(bBase)
			{
				case 0: strLine[nk] = 'A';
					break;
				case 1: strLine[nk] = 'C';
					break;
				case 2: strLine[nk] = 'G';
					break;
				case 3: strLine[nk] = 'T';
					break;
				case 4: strLine[nk] = 'a';
					break;
				case 5: strLine[nk] = 'c';
					break;
				case 6: strLine[nk] = 'g';
					break;
				case 7: strLine[nk] = 't';
					break;
				default: strLine[nk] = 'N';
			}
			nk++;
			nl++;

			if(nk == CODE_BATCH_LEN)
			{
				strLine[nk] = '\0';
				SequenceAddTail(pSeq, strLine);
				nk = 0;
			}
		}
	}
	
	/* last base */
	if(nl < nLen)
	{
		numread = fread( &bChar, sizeof( unsigned char ), 1, fpIn );
		bBase = bChar >> 4;
		switch(bBase)
		{
			case 0: strLine[nk] = 'A';
				break;
			case 1: strLine[nk] = 'C';
				break;
			case 2: strLine[nk] = 'G';
				break;
			case 3: strLine[nk] = 'T';
				break;
			case 4: strLine[nk] = 'a';
				break;
			case 5: strLine[nk] = 'c';
				break;
			case 6: strLine[nk] = 'g';
				break;
			case 7: strLine[nk] = 't';
				break;
			default: strLine[nk] = 'N';
		}
		nk++;
		nl++;
	
		if(nk == CODE_BATCH_LEN)
		{
			strLine[nk] = '\0';
			SequenceAddTail(pSeq, strLine);
			nk = 0;
		}
	}

	if(nl < nLen)
	{
		if(nR2 == 1)
		{
			bBase = bChar & 0x0F;
			switch(bBase)
			{
				case 0: strLine[nk] = 'A';
					break;
				case 1: strLine[nk] = 'C';
					break;
				case 2: strLine[nk] = 'G';
					break;
				case 3: strLine[nk] = 'T';
					break;
				case 4: strLine[nk] = 'a';
					break;
				case 5: strLine[nk] = 'c';
					break;
				case 6: strLine[nk] = 'g';
					break;
				case 7: strLine[nk] = 't';
					break;
				default: strLine[nk] = 'N';
			}
			nk++;
			nl++;

			if(nk == CODE_BATCH_LEN)
			{
				strLine[nk] = '\0';
				SequenceAddTail(pSeq, strLine);
				nk = 0;
			}
		}
	}

	if(nk > 0)
	{
		strLine[nk] = '\0';
		SequenceAddTail(pSeq, strLine);
		nk = 0;
	}

	/* close file */
	fclose(fpIn);

	if(nl != nLen)
	{
		printf("Error: Genome_Code_4bit_GetSeq, incorrect sequence length!\n");
		exit(EXIT_FAILURE);
	}

	/* return */
	return pSeq;
}

/* ----------------------------------------------------------------------- */ 
/*  Genome_Code_4bit_GetUncompressed:                                      */
/*  Get uncompressed code from Genome Sequence coding files.               */
/*  1 byte for 2 bases.                                                    */
/*  0000: A; 0001: C; 0010: G; 0011: T.                                    */
/*  0100: a; 0101: c; 0110: g; 0111: t.                                    */
/*  1000: N;                                                               */ 
/*  Return number of bases read.                                           */
/* ----------------------------------------------------------------------- */ 
int Genome_Code_4bit_GetUncompressed(struct BYTEMATRIX *pSeq, char strInFile[], int nStart, int nEnd)
{
	/* define */
	FILE *fpIn;
	int nP1,nP2,nR1,nR2;
	int numread;
	unsigned char *pVec;

	int ni,nk,nLen;
	unsigned char bChar;

	/* check */
	if(pSeq == NULL)
	{
		printf("Warning: Genome_Code_4bit_GetUncompressed, null sequence vector!\n");
		return 0;
	}
	nLen = nEnd-nStart+1;
	if(nLen <= 0)
	{
		printf("Error: Genome_Code_4bit_GetUncompressed, sequence length <= 0!\n");
		exit(EXIT_FAILURE);
	}
	if(pSeq->nWidth < nLen )
	{
		printf("Error: Genome_Code_4bit_GetUncompressed, sequence length not enough!\n");
		exit(EXIT_FAILURE);
	}
	
	/* init */
	pVec = pSeq->pMatElement;
	fpIn = NULL;
	fpIn = fopen(strInFile, "rb");
	if(fpIn == NULL)
	{
		printf("Error: Genome_Code_4bit_GetUncompressed, cannot open sequence file!\n");
		exit(EXIT_FAILURE);
	}

	/* load seq */
	nP1 = (int)(nStart/2);
	nR1 = nStart%2;
	nP2 = (int)(nEnd/2);
	nR2 = nEnd%2;

	if( fseek( fpIn, nP1, SEEK_SET ) != 0 )
	{
		printf("Error: Genome_Code_4bit_GetUncompressed, cannot locate the required sequence!\n");
		exit(EXIT_FAILURE);
	}

	nk = 0;

	/* first base */
	numread = fread( &bChar, sizeof( unsigned char ), 1, fpIn );
	if(nR1 == 0)
	{
		pVec[nk] = bChar >> 4;
		nk++;
	}

	if(nk < nLen)
	{
		pVec[nk] = bChar & 0x0F;
		nk++;
	}

	/* middle bases */
	if(nk < nLen)
	{
		for(ni=(nP1+1); ni<nP2; ni++)
		{
			numread = fread( &bChar, sizeof( unsigned char ), 1, fpIn );
			pVec[nk] = bChar >> 4;
			nk++;

			pVec[nk] = bChar & 0x0F;
			nk++;
		}
	}
	
	/* last base */
	if(nk < nLen)
	{
		numread = fread( &bChar, sizeof( unsigned char ), 1, fpIn );
		pVec[nk] = bChar >> 4;
		nk++;
	}

	if(nk < nLen)
	{
		if(nR2 == 1)
		{
			pVec[nk] = bChar & 0x0F;
			nk++;
		}
	}

	/* close file */
	fclose(fpIn);

	if(nk != nLen)
	{
		printf("Error: Genome_Code_4bit_GetUncompressed, incorrect sequence length!\n");
		exit(EXIT_FAILURE);
	}

	/* return */
	return nk;
}


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
								 int nStrandType, char strCSFileType[])
{
	/* define */
	int nCount;
	char strSeqFile[LINE_LENGTH];
	char strOutFile[LINE_LENGTH];
	char strConsFile[LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	char strChrName[LINE_LENGTH];
	char strAlias[LINE_LENGTH];
	FILE *fpIn;
	FILE *fpOut;

	int nChr,nStart,nEnd;
	char chStrand;
	struct tagSequence *pSeq;
	struct BYTEMATRIX *pCS;
	struct INTMATRIX *pChrLen;

	/* init */
	nCount = 0;
	AdjustDirectoryPath(strGenomePath);
	AdjustDirectoryPath(strConsPath);
	AdjustDirectoryPath(strOutPath);
	
	/* load chrmosome length */
	sprintf(strLine, "%schrlen.txt", strGenomePath);
	pChrLen = NULL;
	pChrLen = IMLOAD(strLine);
	if(pChrLen == NULL)
	{
		printf("Error: cannot load chromosome length!\n");
		exit(EXIT_FAILURE);
	}

	/* open files */
	fpIn = NULL;
	fpIn = fopen(strInFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: Genome_Code_4bit_GetSeq_Main, cannot open input file.\n");
		return 0;
	}

	sprintf(strOutFile, "%s%s.fa", strOutPath, strOutTitle);
	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: Genome_Code_4bit_GetSeq_Main, cannot open output file.\n");
		return 0;
	}

	/* get seq */
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;
		
		sscanf(strLine, "%s %d %d %d %c", strAlias, &nChr, &nStart, &nEnd, &chStrand);
		Genome_Index_To_ChromosomeName(strChrName, strSpecies, nChr);

		if( (nChr <= 0) || (nChr > pChrLen->nHeight))
		{
			printf("Warning: (%d:%d-%d) is out of normal genome range, sequence not obtained!\n", nChr, nStart, nEnd);
			continue;
		}

		if( (nStart < 0) || (nEnd >= pChrLen->pMatElement[nChr-1]) )
		{
			printf("Warning: (%d:%d-%d) is out of normal genome range, sequence not obtained!\n", nChr, nStart, nEnd);
			continue;
		}

		sprintf(strSeqFile, "%s%s.sq", strGenomePath, strChrName);
		pSeq = NULL;
		pSeq = Genome_Code_4bit_GetSeq(strSeqFile, nStart, nEnd);
		if(pSeq == NULL)
			continue;

		pSeq->m_nIndex = nCount;
		strcpy(pSeq->m_strAlias, strAlias);
		if(nStrandType == 1)
			SequenceWriteToFasta_ByStrand(pSeq, fpOut, chStrand, 1);
		else
			SequenceWriteToFasta_ByStrand(pSeq, fpOut, '+', 1);
		SequenceDelete(pSeq);

		sprintf(strConsFile, "%s%s.cs", strConsPath, strChrName);
		pCS = NULL;
		pCS = Genome_GetConsScore(strConsFile, nStart, nEnd);
		if(pCS == NULL)
			continue;

		
		if(nStrandType == 1)
		{
			if(strcmp(strCSFileType, "cs") == 0)
			{
				sprintf(strConsFile, "%s%s_%d_%s.cs", strOutPath, strOutTitle, nCount, strAlias);
				ConsScoreWriteToBinaryFile_ByStrand(pCS, strConsFile, chStrand);
			}
			else if(strcmp(strCSFileType, "txt") == 0)
			{
				sprintf(strConsFile, "%s%s_%d_%s.txt", strOutPath, strOutTitle, nCount, strAlias);
				ConsScoreWriteToTextFile_ByStrand(pCS, strConsFile, chStrand);
			}
			else if(strcmp(strCSFileType, "bed") == 0)
			{
				sprintf(strConsFile, "%s%s_%d_%s.bed", strOutPath, strOutTitle, nCount, strAlias);
				ConsScoreWriteToBedFile_ByStrand(strChrName, nStart, nEnd, pCS, strConsFile, chStrand);
			}
		}
		else
		{
			if(strcmp(strCSFileType, "cs") == 0)
			{
				sprintf(strConsFile, "%s%s_%d_%s.cs", strOutPath, strOutTitle, nCount, strAlias);
				ConsScoreWriteToBinaryFile_ByStrand(pCS, strConsFile, '+');
			}
			else if(strcmp(strCSFileType, "txt") == 0)
			{
				sprintf(strConsFile, "%s%s_%d_%s.txt", strOutPath, strOutTitle, nCount, strAlias);
				ConsScoreWriteToTextFile_ByStrand(pCS, strConsFile, '+');
			}
			else if(strcmp(strCSFileType, "bed") == 0)
			{
				sprintf(strConsFile, "%s%s_%d_%s.bed", strOutPath, strOutTitle, nCount, strAlias);
				ConsScoreWriteToBedFile_ByStrand(strChrName, nStart, nEnd, pCS, strConsFile, '+');
			}
		}

		DestroyByteMatrix(pCS);

		nCount++;
	}

	/* close file */
	fclose(fpIn);
	fclose(fpOut);
	DestroyIntMatrix(pChrLen);

	/* return */
	return nCount;
}

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
								 int nStrandType, char strCSFileType[])
{
	/* define */
	int nCount;
	char strSeqFile[LINE_LENGTH];
	char strOutFile[LINE_LENGTH];
	char strConsFile[LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	char strChrName[LINE_LENGTH];
	char strAlias[LINE_LENGTH];
	FILE *fpIn;
	FILE *fpOut;

	int nChr,nStart,nEnd;
	char chStrand;
	struct tagSequence *pSeq;
	struct BYTEMATRIX *pCS;
	struct INTMATRIX *pChrLen;

	/* init */
	nCount = 0;
	AdjustDirectoryPath(strGenomePath);
	AdjustDirectoryPath(strConsPath);
	AdjustDirectoryPath(strOutPath);
	
	/* load chrmosome length */
	sprintf(strLine, "%schrlen.txt", strGenomePath);
	pChrLen = NULL;
	pChrLen = IMLOAD(strLine);
	if(pChrLen == NULL)
	{
		printf("Error: cannot load chromosome length!\n");
		exit(EXIT_FAILURE);
	}

	/* open files */
	fpIn = NULL;
	fpIn = fopen(strInFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: Genome_Code_4bit_GetSeq_Main, cannot open input file.\n");
		return 0;
	}

	sprintf(strOutFile, "%s%s.fa", strOutPath, strOutTitle);
	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: Genome_Code_4bit_GetSeq_Main, cannot open output file.\n");
		return 0;
	}

	/* get seq */
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;
		
		sscanf(strLine, "%s %s %d %d %c", strAlias, strChrName, &nStart, &nEnd, &chStrand);
		nChr = Genome_ChromosomeName_To_Index(strChrName, strSpecies);

		if( (nChr <= 0) || (nChr > pChrLen->nHeight))
		{
			printf("Warning: (%s:%d-%d) is out of normal genome range, sequence not obtained!\n", strChrName, nStart, nEnd);
			continue;
		}

		if( (nStart < 0) || (nEnd >= pChrLen->pMatElement[nChr-1]) )
		{
			printf("Warning: (%s:%d-%d) is out of normal genome range, sequence not obtained!\n", strChrName, nStart, nEnd);
			continue;
		}

		sprintf(strSeqFile, "%s%s.sq", strGenomePath, strChrName);
		pSeq = NULL;
		pSeq = Genome_Code_4bit_GetSeq(strSeqFile, nStart, nEnd);
		if(pSeq == NULL)
			continue;

		pSeq->m_nIndex = nCount;
		strcpy(pSeq->m_strAlias, strAlias);
		if(nStrandType == 1)
			SequenceWriteToFasta_ByStrand(pSeq, fpOut, chStrand, 1);
		else
			SequenceWriteToFasta_ByStrand(pSeq, fpOut, '+', 1);
		SequenceDelete(pSeq);

		sprintf(strConsFile, "%s%s.cs", strConsPath, strChrName);
		pCS = NULL;
		pCS = Genome_GetConsScore(strConsFile, nStart, nEnd);
		if(pCS == NULL)
			continue;

		
		if(nStrandType == 1)
		{
			if(strcmp(strCSFileType, "cs") == 0)
			{
				sprintf(strConsFile, "%s%s_%d_%s.cs", strOutPath, strOutTitle, nCount, strAlias);
				ConsScoreWriteToBinaryFile_ByStrand(pCS, strConsFile, chStrand);
			}
			else if(strcmp(strCSFileType, "txt") == 0)
			{
				sprintf(strConsFile, "%s%s_%d_%s.txt", strOutPath, strOutTitle, nCount, strAlias);
				ConsScoreWriteToTextFile_ByStrand(pCS, strConsFile, chStrand);
			}
			else if(strcmp(strCSFileType, "bed") == 0)
			{
				sprintf(strConsFile, "%s%s_%d_%s.bed", strOutPath, strOutTitle, nCount, strAlias);
				ConsScoreWriteToBedFile_ByStrand(strChrName, nStart, nEnd, pCS, strConsFile, chStrand);
			}
		}
		else
		{
			if(strcmp(strCSFileType, "cs") == 0)
			{
				sprintf(strConsFile, "%s%s_%d_%s.cs", strOutPath, strOutTitle, nCount, strAlias);
				ConsScoreWriteToBinaryFile_ByStrand(pCS, strConsFile, '+');
			}
			else if(strcmp(strCSFileType, "txt") == 0)
			{
				sprintf(strConsFile, "%s%s_%d_%s.txt", strOutPath, strOutTitle, nCount, strAlias);
				ConsScoreWriteToTextFile_ByStrand(pCS, strConsFile, '+');
			}
			else if(strcmp(strCSFileType, "bed") == 0)
			{
				sprintf(strConsFile, "%s%s_%d_%s.bed", strOutPath, strOutTitle, nCount, strAlias);
				ConsScoreWriteToBedFile_ByStrand(strChrName, nStart, nEnd, pCS, strConsFile, '+');
			}
		}

		DestroyByteMatrix(pCS);

		nCount++;
	}

	/* close file */
	fclose(fpIn);
	fclose(fpOut);
	DestroyIntMatrix(pChrLen);

	/* return */
	return nCount;
}

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
			char strInFile[], char strOutFile[], int nStrandType)
{
	/* define */
	int nCount;
	char strSeqFile[LINE_LENGTH];
	char strConsFile[LINE_LENGTH];
	char strCDSFile[LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	char strChrName[LINE_LENGTH];
	char strAlias[LINE_LENGTH];
	FILE *fpIn;
	FILE *fpOut;
	int ni;
	char *vSeq;

	int nChr,nStart,nEnd;
	char chStrand;
	struct tagSequence *pSeq;
	struct BYTEMATRIX *pCS;
	struct BYTEMATRIX *pCDS;
	struct INTMATRIX *pChrLen;

	/* init */
	nCount = 0;
	AdjustDirectoryPath(strGenomePath);
	if(nUseCS == 1)
		AdjustDirectoryPath(strConsPath);
	if(nUseCDS == 1)
		AdjustDirectoryPath(strCDSPath);

	/* load chrmosome length */
	sprintf(strLine, "%schrlen.txt", strGenomePath);
	pChrLen = NULL;
	pChrLen = IMLOAD(strLine);
	if(pChrLen == NULL)
	{
		printf("Error: cannot load chromosome length!\n");
		exit(EXIT_FAILURE);
	}

	/* open files */
	fpIn = NULL;
	fpIn = fopen(strInFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: Genome_Code_4bit_GetMaskedSeq_Main, cannot open input file.\n");
		return 0;
	}

	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: Genome_Code_4bit_GetMaskedSeq_Main, cannot open output file.\n");
		return 0;
	}

	/* get seq */
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;
		
		sscanf(strLine, "%s %d %d %d %c", strAlias, &nChr, &nStart, &nEnd, &chStrand);
		Genome_Index_To_ChromosomeName(strChrName, strSpecies, nChr);
		
		if( (nChr <= 0) || (nChr > pChrLen->nHeight))
		{
			printf("Warning: (%d:%d-%d) is out of normal genome range, sequence not obtained!\n", nChr, nStart, nEnd);
			continue;
		}

		if( (nStart < 0) || (nEnd >= pChrLen->pMatElement[nChr-1]) )
		{
			printf("Warning: (%d:%d-%d) is out of normal genome range, sequence not obtained!\n", nChr, nStart, nEnd);
			continue;
		}
		
		sprintf(strSeqFile, "%s%s.sq", strGenomePath, strChrName);
		pSeq = NULL;
		pSeq = Genome_Code_4bit_GetSeq(strSeqFile, nStart, nEnd);
		if(pSeq == NULL)
			continue;

		pSeq->m_nIndex = nCount;
		strcpy(pSeq->m_strAlias, strAlias);

		vSeq = pSeq->m_pSequence->m_pString;
		
		if(nUseCS == 1)
		{
			sprintf(strConsFile, "%s%s.cs", strConsPath, strChrName);
			pCS = NULL;
			pCS = Genome_GetConsScore(strConsFile, nStart, nEnd);
			if(pCS == NULL)
				continue;

			for(ni=0; ni<pSeq->m_nLength; ni++)
			{
				if( (double)(pCS->pMatElement[ni]) < dC )
					vSeq[ni] = 'N';
			}

			DestroyByteMatrix(pCS);
		}

		if(nUseCDS == 1)
		{
			sprintf(strCDSFile, "%s%s.cds", strCDSPath, strChrName);
			pCDS = NULL;
			pCDS = Genome_GetCDS(strCDSFile, nStart, nEnd);
			if(pCDS == NULL)
				continue;
			
			for(ni=0; ni<pSeq->m_nLength; ni++)
			{
				if( pCDS->pMatElement[ni] != 0 )
					vSeq[ni] = 'N';
			}

			DestroyByteMatrix(pCDS);
		}

		/* hard mask repeats */
		for(ni=0; ni<pSeq->m_nLength; ni++)
		{
			if( (vSeq[ni] != 'A') && (vSeq[ni] != 'C') && (vSeq[ni] != 'G') && (vSeq[ni] != 'T') )
				vSeq[ni] = 'N';
		}

		if(nStrandType == 1)
			SequenceWriteToFasta_ByStrand(pSeq, fpOut, chStrand, 1);
		else
			SequenceWriteToFasta_ByStrand(pSeq, fpOut, '+', 1);
		SequenceDelete(pSeq);

		nCount++;
	}

	/* close file */
	fclose(fpIn);
	fclose(fpOut);
	DestroyIntMatrix(pChrLen);

	/* return */
	return nCount;
}

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
			char strInFile[], char strOutFile[], int nStrandType)
{
	/* define */
	int nCount;
	char strSeqFile[LINE_LENGTH];
	char strConsFile[LINE_LENGTH];
	char strCDSFile[LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	char strChrName[LINE_LENGTH];
	char strAlias[LINE_LENGTH];
	FILE *fpIn;
	FILE *fpOut;
	int ni;
	char *vSeq;

	int nChr,nStart,nEnd;
	char chStrand;
	struct tagSequence *pSeq;
	struct BYTEMATRIX *pCS;
	struct BYTEMATRIX *pCDS;
	struct INTMATRIX *pChrLen;

	/* init */
	nCount = 0;
	AdjustDirectoryPath(strGenomePath);
	if(nUseCS == 1)
		AdjustDirectoryPath(strConsPath);
	if(nUseCDS == 1)
		AdjustDirectoryPath(strCDSPath);

	/* load chrmosome length */
	sprintf(strLine, "%schrlen.txt", strGenomePath);
	pChrLen = NULL;
	pChrLen = IMLOAD(strLine);
	if(pChrLen == NULL)
	{
		printf("Error: cannot load chromosome length!\n");
		exit(EXIT_FAILURE);
	}

	/* open files */
	fpIn = NULL;
	fpIn = fopen(strInFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: Genome_Code_4bit_GetMaskedSeq_Main, cannot open input file.\n");
		return 0;
	}

	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: Genome_Code_4bit_GetMaskedSeq_Main, cannot open output file.\n");
		return 0;
	}

	/* get seq */
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;
		
		sscanf(strLine, "%s %s %d %d %c", strAlias, strChrName, &nStart, &nEnd, &chStrand);
		nChr = Genome_ChromosomeName_To_Index(strChrName, strSpecies);
		
		if( (nChr <= 0) || (nChr > pChrLen->nHeight))
		{
			printf("Warning: (%s:%d-%d) is out of normal genome range, sequence not obtained!\n", strChrName, nStart, nEnd);
			continue;
		}

		if( (nStart < 0) || (nEnd >= pChrLen->pMatElement[nChr-1]) )
		{
			printf("Warning: (%s:%d-%d) is out of normal genome range, sequence not obtained!\n", strChrName, nStart, nEnd);
			continue;
		}
		
		sprintf(strSeqFile, "%s%s.sq", strGenomePath, strChrName);
		pSeq = NULL;
		pSeq = Genome_Code_4bit_GetSeq(strSeqFile, nStart, nEnd);
		if(pSeq == NULL)
			continue;

		pSeq->m_nIndex = nCount;
		strcpy(pSeq->m_strAlias, strAlias);

		vSeq = pSeq->m_pSequence->m_pString;
		
		if(nUseCS == 1)
		{
			sprintf(strConsFile, "%s%s.cs", strConsPath, strChrName);
			pCS = NULL;
			pCS = Genome_GetConsScore(strConsFile, nStart, nEnd);
			if(pCS == NULL)
				continue;

			for(ni=0; ni<pSeq->m_nLength; ni++)
			{
				if( (double)(pCS->pMatElement[ni]) < dC )
					vSeq[ni] = 'N';
			}

			DestroyByteMatrix(pCS);
		}

		if(nUseCDS == 1)
		{
			sprintf(strCDSFile, "%s%s.cds", strCDSPath, strChrName);
			pCDS = NULL;
			pCDS = Genome_GetCDS(strCDSFile, nStart, nEnd);
			if(pCDS == NULL)
				continue;
			
			for(ni=0; ni<pSeq->m_nLength; ni++)
			{
				if( pCDS->pMatElement[ni] != 0 )
					vSeq[ni] = 'N';
			}

			DestroyByteMatrix(pCDS);
		}

		/* hard mask repeats */
		for(ni=0; ni<pSeq->m_nLength; ni++)
		{
			if( (vSeq[ni] != 'A') && (vSeq[ni] != 'C') && (vSeq[ni] != 'G') && (vSeq[ni] != 'T') )
				vSeq[ni] = 'N';
		}

		if(nStrandType == 1)
			SequenceWriteToFasta_ByStrand(pSeq, fpOut, chStrand, 1);
		else
			SequenceWriteToFasta_ByStrand(pSeq, fpOut, '+', 1);
		SequenceDelete(pSeq);

		nCount++;
	}

	/* close file */
	fclose(fpIn);
	fclose(fpOut);
	DestroyIntMatrix(pChrLen);

	/* return */
	return nCount;
}

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
			char strInputFile[], char strOutputFile[])
{
	/* define */
	int nCount;
	char strSeqFile[LINE_LENGTH];
	char strConsFile[LINE_LENGTH];
	char strCDSFile[LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	char strChrName[LINE_LENGTH];
	char strAlias[LINE_LENGTH];
	FILE *fpIn;
	FILE *fpOut;
	int ni;
	char *vSeq;
	int nNonRep;
	double dNonRep;
	int nNonCS;
	double dNonCS;
	int nNonCDS;
	double dNonCDS;
	int nSeqLen;
	int nOK;

	int nChr,nStart,nEnd;
	char chStrand;
	struct tagSequence *pSeq;
	struct BYTEMATRIX *pCS;
	struct BYTEMATRIX *pCDS;
	struct INTMATRIX *pChrLen;

	/* init */
	nCount = 0;
	AdjustDirectoryPath(strGenomePath);
	if(nUseCS == 1)
		AdjustDirectoryPath(strConsPath);
	if(nUseCDS == 1)
		AdjustDirectoryPath(strCdsPath);

	/* load chrmosome length */
	sprintf(strLine, "%schrlen.txt", strGenomePath);
	pChrLen = NULL;
	pChrLen = IMLOAD(strLine);
	if(pChrLen == NULL)
	{
		printf("Error: cannot load chromosome length!\n");
		exit(EXIT_FAILURE);
	}

	/* open files */
	fpIn = NULL;
	fpIn = fopen(strInputFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: Genome_Code_4bit_GetMaskedReg_Main, cannot open input file.\n");
		return 0;
	}

	fpOut = NULL;
	fpOut = fopen(strOutputFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: Genome_Code_4bit_GetMaskedReg_Main, cannot open output file.\n");
		return 0;
	}

	/* get seq */
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;
		
		sscanf(strLine, "%s %d %d %d %c", strAlias, &nChr, &nStart, &nEnd, &chStrand);
		Genome_Index_To_ChromosomeName(strChrName, strSpecies, nChr);

		if( (nChr <= 0) || (nChr > pChrLen->nHeight))
		{
			printf("Warning: (%s:%d-%d) is out of normal genome range, sequence not obtained!\n", strChrName, nStart, nEnd);
			continue;
		}

		if( (nStart < 0) || (nEnd >= pChrLen->pMatElement[nChr-1]) || (nStart > nEnd) )
		{
			printf("Warning: (%s:%d-%d) is out of normal genome range, sequence not obtained!\n", strChrName, nStart, nEnd);
			continue;
		}
		
		nOK = 1;
		nSeqLen = nEnd-nStart+1;


		if(nUseRepMask == 1)
		{
			sprintf(strSeqFile, "%s%s.sq", strGenomePath, strChrName);
			pSeq = NULL;
			pSeq = Genome_Code_4bit_GetSeq(strSeqFile, nStart, nEnd);
			if(pSeq == NULL)
				continue;

			pSeq->m_nIndex = nCount;
			strcpy(pSeq->m_strAlias, strAlias);

			vSeq = pSeq->m_pSequence->m_pString;
			
			nNonRep = 0;
			dNonRep = 0.0;
			for(ni=0; ni<pSeq->m_nLength; ni++)
			{
				if( (vSeq[ni] == 'A') || (vSeq[ni] == 'C') || (vSeq[ni] == 'G') || (vSeq[ni] == 'T') )
					nNonRep++;
			}
			if(pSeq->m_nLength > 0)
			{
				dNonRep = ((double)nNonRep)/((double)(pSeq->m_nLength));
			}
			if(dNonRep < dR)
			{
				nOK = 0;
			}

			SequenceDelete(pSeq);
		}

		if( (nOK == 1) && (nUseCS == 1) )
		{
			sprintf(strConsFile, "%s%s.cs", strConsPath, strChrName);
			pCS = NULL;
			pCS = Genome_GetConsScore(strConsFile, nStart, nEnd);
			if(pCS == NULL)
				continue;

			nNonCS = 0;
			dNonCS = 1.0;
			for(ni=0; ni<nSeqLen; ni++)
			{
				if( (double)(pCS->pMatElement[ni]) < dC )
					nNonCS++;
			}

			if(nSeqLen > 0)
			{
				dNonCS = ((double)nNonCS)/((double)(nSeqLen));
			}

			if((1.0-dNonCS) < dCR)
			{
				nOK = 0;
			}

			DestroyByteMatrix(pCS);
		}

		if( (nOK == 1) && (nUseCDS == 1) )
		{
			sprintf(strCDSFile, "%s%s.cds", strCdsPath, strChrName);
			pCDS = NULL;
			pCDS = Genome_GetCDS(strCDSFile, nStart, nEnd);
			if(pCDS == NULL)
				continue;
			
			nNonCDS = 0;
			dNonCDS = 0.0;

			for(ni=0; ni<nSeqLen; ni++)
			{
				if( pCDS->pMatElement[ni] == 0 )
					nNonCDS++;
			}

			if(nSeqLen > 0)
			{
				dNonCDS = ((double)nNonCDS)/((double)(nSeqLen));
			}

			if( dNonCDS < dCDS)
			{
				nOK = 0;
			}

			DestroyByteMatrix(pCDS);
		}

		/* write to files */
		if(nOK == 1)
		{
			fprintf(fpOut, "%s\n", strLine);
		}
		
		nCount++;
	}

	/* close file */
	fclose(fpIn);
	fclose(fpOut);
	DestroyIntMatrix(pChrLen);

	/* return */
	return nCount;
}

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
			char strInputFile[], char strOutputFile[])
{
	/* define */
	int nCount;
	char strSeqFile[LINE_LENGTH];
	char strConsFile[LINE_LENGTH];
	char strCDSFile[LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	char strChrName[LINE_LENGTH];
	char strAlias[LINE_LENGTH];
	FILE *fpIn;
	FILE *fpOut;
	int ni;
	char *vSeq;
	int nNonRep;
	double dNonRep;
	int nNonCS;
	double dNonCS;
	int nNonCDS;
	double dNonCDS;
	int nSeqLen;
	int nOK;

	int nChr,nStart,nEnd;
	char chStrand;
	struct tagSequence *pSeq;
	struct BYTEMATRIX *pCS;
	struct BYTEMATRIX *pCDS;
	struct INTMATRIX *pChrLen;

	/* init */
	nCount = 0;
	AdjustDirectoryPath(strGenomePath);
	if(nUseCS == 1)
		AdjustDirectoryPath(strConsPath);
	if(nUseCDS == 1)
		AdjustDirectoryPath(strCdsPath);

	/* load chrmosome length */
	sprintf(strLine, "%schrlen.txt", strGenomePath);
	pChrLen = NULL;
	pChrLen = IMLOAD(strLine);
	if(pChrLen == NULL)
	{
		printf("Error: cannot load chromosome length!\n");
		exit(EXIT_FAILURE);
	}

	/* open files */
	fpIn = NULL;
	fpIn = fopen(strInputFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: Genome_Code_4bit_GetMaskedReg_Main, cannot open input file.\n");
		return 0;
	}

	fpOut = NULL;
	fpOut = fopen(strOutputFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: Genome_Code_4bit_GetMaskedReg_Main, cannot open output file.\n");
		return 0;
	}

	/* get seq */
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		sscanf(strLine, "%s %s %d %d %c", strAlias, strChrName, &nStart, &nEnd, &chStrand);
		nChr = Genome_ChromosomeName_To_Index(strChrName, strSpecies);
		
		if( (nChr <= 0) || (nChr > pChrLen->nHeight))
		{
			printf("Warning: (%s:%d-%d) is out of normal genome range, sequence not obtained!\n", strChrName, nStart, nEnd);
			continue;
		}

		if( (nStart < 0) || (nEnd >= pChrLen->pMatElement[nChr-1]) || (nStart > nEnd) )
		{
			printf("Warning: (%s:%d-%d) is out of normal genome range, sequence not obtained!\n", strChrName, nStart, nEnd);
			continue;
		}
		
		nOK = 1;
		nSeqLen = nEnd-nStart+1;


		if(nUseRepMask == 1)
		{
			sprintf(strSeqFile, "%s%s.sq", strGenomePath, strChrName);
			pSeq = NULL;
			pSeq = Genome_Code_4bit_GetSeq(strSeqFile, nStart, nEnd);
			if(pSeq == NULL)
				continue;

			pSeq->m_nIndex = nCount;
			strcpy(pSeq->m_strAlias, strAlias);

			vSeq = pSeq->m_pSequence->m_pString;
			
			nNonRep = 0;
			dNonRep = 0.0;
			for(ni=0; ni<pSeq->m_nLength; ni++)
			{
				if( (vSeq[ni] == 'A') || (vSeq[ni] == 'C') || (vSeq[ni] == 'G') || (vSeq[ni] == 'T') )
					nNonRep++;
			}
			if(pSeq->m_nLength > 0)
			{
				dNonRep = ((double)nNonRep)/((double)(pSeq->m_nLength));
			}
			if(dNonRep < dR)
			{
				nOK = 0;
			}

			SequenceDelete(pSeq);
		}

		if( (nOK == 1) && (nUseCS == 1) )
		{
			sprintf(strConsFile, "%s%s.cs", strConsPath, strChrName);
			pCS = NULL;
			pCS = Genome_GetConsScore(strConsFile, nStart, nEnd);
			if(pCS == NULL)
				continue;

			nNonCS = 0;
			dNonCS = 1.0;
			for(ni=0; ni<nSeqLen; ni++)
			{
				if( (double)(pCS->pMatElement[ni]) < dC )
					nNonCS++;
			}

			if(nSeqLen > 0)
			{
				dNonCS = ((double)nNonCS)/((double)(nSeqLen));
			}

			if((1.0-dNonCS) < dCR)
			{
				nOK = 0;
			}

			DestroyByteMatrix(pCS);
		}

		if( (nOK == 1) && (nUseCDS == 1) )
		{
			sprintf(strCDSFile, "%s%s.cds", strCdsPath, strChrName);
			pCDS = NULL;
			pCDS = Genome_GetCDS(strCDSFile, nStart, nEnd);
			if(pCDS == NULL)
				continue;
			
			nNonCDS = 0;
			dNonCDS = 0.0;

			for(ni=0; ni<nSeqLen; ni++)
			{
				if( pCDS->pMatElement[ni] == 0 )
					nNonCDS++;
			}

			if(nSeqLen > 0)
			{
				dNonCDS = ((double)nNonCDS)/((double)(nSeqLen));
			}

			if( dNonCDS < dCDS)
			{
				nOK = 0;
			}

			DestroyByteMatrix(pCDS);
		}

		/* write to files */
		if(nOK == 1)
		{
			fprintf(fpOut, "%s\n", strLine);
		}
		
		nCount++;
	}

	/* close file */
	fclose(fpIn);
	fclose(fpOut);
	DestroyIntMatrix(pChrLen);

	/* return */
	return nCount;
}

/* ----------------------------------------------------------------------- */ 
/*  Genome_GetConsScore:                                                   */
/*  Get conservation scores from Genome Sequence coding files.             */
/*  1 byte for 2 bases.                                                    */
/*  0000: A; 0001: C; 0010: G; 0011: T.                                    */
/*  0100: a; 0101: c; 0110: g; 0111: t.                                    */
/*  1000: N;                                                               */ 
/* ----------------------------------------------------------------------- */ 
struct BYTEMATRIX *Genome_GetConsScore(char strConsFile[], int nStart, int nEnd)
{
	/* define */
	struct BYTEMATRIX *pCS;
	FILE *fpIn;
	int nLen;
	int numread;

	/* init */
	nLen = nEnd-nStart+1;
	if(nLen <= 0)
		return NULL;

	fpIn = NULL;
	fpIn = fopen(strConsFile, "rb");
	if(fpIn == NULL)
	{
		printf("Error: Genome_GetConsScore, cannot open conservation score file!\n");
		exit(EXIT_FAILURE);
	}

	pCS = NULL;
	pCS = CreateByteMatrix(1, nLen);
	if(pCS == NULL)
	{
		printf("Error: Genome_GetConsScore, cannot create scores.\n");
		exit(EXIT_FAILURE);
	}

	/* load score */
	if( fseek( fpIn, nStart, SEEK_SET ) != 0 )
	{
		printf("Error: Genome_GetConsScore, cannot locate the required sequence!\n");
		exit(EXIT_FAILURE);
	}

	/* read */
	numread = fread(pCS->pMatElement, sizeof(unsigned char), nLen, fpIn);
	if(numread != nLen)
	{
		printf("Error: Genome_GetConsScore, loading error!\n");
		exit(EXIT_FAILURE);
	}

	fclose(fpIn);

	/* return */
	return pCS;
}

/* ----------------------------------------------------------------------- */ 
/*  Genome_GetCDS:                                                         */
/*  Get CDS status from Genome CDS coding files.                           */
/*  1 byte for each base. 0=non cds, 1=cds.                                */
/* ----------------------------------------------------------------------- */ 
struct BYTEMATRIX *Genome_GetCDS(char strCDSFile[], int nStart, int nEnd)
{
	/* define */
	struct BYTEMATRIX *pCDS;
	FILE *fpIn;
	int nLen;
	int numread;

	/* init */
	nLen = nEnd-nStart+1;
	if(nLen <= 0)
		return NULL;

	fpIn = NULL;
	fpIn = fopen(strCDSFile, "rb");
	if(fpIn == NULL)
	{
		printf("Error: Genome_GetCDS, cannot open cds file!\n");
		exit(EXIT_FAILURE);
	}

	pCDS = NULL;
	pCDS = CreateByteMatrix(1, nLen);
	if(pCDS == NULL)
	{
		printf("Error: Genome_GetCDS, cannot create scores.\n");
		exit(EXIT_FAILURE);
	}

	/* load score */
	if( fseek( fpIn, nStart, SEEK_SET ) != 0 )
	{
		printf("Error: Genome_GetCDS, cannot locate the required sequence!\n");
		exit(EXIT_FAILURE);
	}

	/* read */
	numread = fread(pCDS->pMatElement, sizeof(unsigned char), nLen, fpIn);
	if(numread != nLen)
	{
		printf("Error: Genome_GetCDS, loading error!\n");
		exit(EXIT_FAILURE);
	}

	fclose(fpIn);

	/* return */
	return pCDS;
}


/* ----------------------------------------------------------------------- */ 
/*  ConsScoreWriteToBinaryFile_ByStrand:                                   */
/*  write conservation score to files                                      */
/* ----------------------------------------------------------------------- */ 
int ConsScoreWriteToBinaryFile_ByStrand(struct BYTEMATRIX *pCS, char strConsFile[], 
								  char chStrand)
{
	/* nLinePos and pSeqPos are used for accessing the sequence. */
	FILE *fpOut;
	int numwritten;
	struct BYTEMATRIX *pRCS;
	int ni;
	unsigned char *pS1,*pS2;

	/* check parameter */
	if(pCS == NULL)
		return PROC_FAILURE;

	fpOut = NULL;
	fpOut = fopen(strConsFile, "wb");
	if(fpOut == NULL)
	{
		printf("Error: cannot open output file for conservation score!\n");
		return PROC_FAILURE;
	}

	/* write */
	/* - strand */
	if(chStrand == '-')
	{
		pRCS = NULL;
		pRCS = CreateByteMatrix(1, pCS->nWidth);

		pS1 = pCS->pMatElement;
		pS2 = pRCS->pMatElement+pCS->nWidth-1;
		for(ni=0; ni<pCS->nWidth; ni++)
		{
			*pS2 = *pS1;
			pS1++;
			pS2--;
		}

		numwritten = fwrite(pRCS->pMatElement, sizeof(unsigned char), pRCS->nWidth, fpOut);
		DestroyByteMatrix(pRCS);
	}
	/* + strand */
	else
	{
		numwritten = fwrite(pCS->pMatElement, sizeof(unsigned char), pCS->nWidth, fpOut);
	}

	if(numwritten != pCS->nWidth)
	{
		printf("Error: ConsScoreWriteToFile_ByStrand, writing error!\n");
		exit(EXIT_FAILURE);
	}
	

	/* close */
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  ConsScoreWriteToTextFile_ByStrand:                                   */
/*  write conservation score to files                                      */
/* ----------------------------------------------------------------------- */ 
int ConsScoreWriteToTextFile_ByStrand(struct BYTEMATRIX *pCS, char strConsFile[], 
								  char chStrand)
{
	/* nLinePos and pSeqPos are used for accessing the sequence. */
	struct BYTEMATRIX *pRCS;
	int ni;
	unsigned char *pS1,*pS2;

	/* check parameter */
	if(pCS == NULL)
		return PROC_FAILURE;

	/* write */
	/* - strand */
	if(chStrand == '-')
	{
		pRCS = NULL;
		pRCS = CreateByteMatrix(pCS->nWidth, 1);

		pS1 = pCS->pMatElement;
		pS2 = pRCS->pMatElement+pCS->nWidth-1;
		for(ni=0; ni<pCS->nWidth; ni++)
		{
			*pS2 = *pS1;
			pS1++;
			pS2--;
		}

		BMSAVE(pRCS, strConsFile);

		DestroyByteMatrix(pRCS);
	}
	/* + strand */
	else
	{
		pRCS = NULL;
		pRCS = CreateByteMatrix(pCS->nWidth, 1);

		pS1 = pCS->pMatElement;
		pS2 = pRCS->pMatElement;
		for(ni=0; ni<pCS->nWidth; ni++)
		{
			*pS2 = *pS1;
			pS1++;
			pS2++;
		}

		BMSAVE(pRCS, strConsFile);

		DestroyByteMatrix(pRCS);
	}

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  ConsScoreWriteToBedFile_ByStrand:                                      */
/*  write conservation score to files                                      */
/* ----------------------------------------------------------------------- */ 
int ConsScoreWriteToBedFile_ByStrand(char strChrName[], int nStart, int nEnd, 
									 struct BYTEMATRIX *pCS, char strConsFile[], 
									 char chStrand)
{
	/* nLinePos and pSeqPos are used for accessing the sequence. */
	int ni,nScore;
	FILE *fpOut;

	/* check parameter */
	if(pCS == NULL)
		return PROC_FAILURE;

	/* write */
	fpOut = NULL;
	fpOut = fopen(strConsFile, "wt");
	if(fpOut == NULL)
	{
		printf("Error: cannot open conservation output file!\n");
		return PROC_FAILURE;
	}

	fprintf(fpOut, "track name=score description=\"conservation score\" useScore=1\n");
	for(ni=nStart; ni<=nEnd; ni++)
	{
		nScore = (int)(pCS->pMatElement[ni-nStart]*1000/255);
		if(nScore > 1000)
			nScore = 1000;
		fprintf(fpOut, "%s\t%d\t%d\ts\t%d\n", strChrName, ni, (ni+1), nScore);
	}

	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Genome_Index_To_ChromosomeName:                                        */
/*  convert numeric chromosome id to chromsome name                        */
/* ----------------------------------------------------------------------- */ 
char *Genome_Index_To_ChromosomeName(char *strChrName, char strSpecies[], int nChr)
{
	/* define */
	char strLine[LINE_LENGTH];
	char strChr[LINE_LENGTH];

	/* convert */
	strcpy(strLine, strSpecies);
	StrMakeUpper(strLine);

	/* human */
	if(strcmp(strLine, "HUMAN") == 0)
	{
		if( (nChr >= 1) && (nChr <= 22) )
		{
			sprintf(strChr, "chr%d", nChr);
		}
		else if(nChr == 23)
		{ 
			strcpy(strChr, "chrX");
		}
		else if(nChr == 24)
		{
			strcpy(strChr, "chrY");
		}
		else if(nChr == 25)
		{
			strcpy(strChr, "chrM");
		}
		else
		{
			strcpy(strChr, "unplaced");
		}
	}

	/* rhesus monkey */
	else if(strcmp(strLine, "RHESUS") == 0)
	{
		if( (nChr >= 1) && (nChr <= 20) )
		{
			sprintf(strChr, "chr%d", nChr);
		}
		else if(nChr == 21)
		{ 
			strcpy(strChr, "chrX");
		}
		else if(nChr == 22)
		{
			strcpy(strChr, "chrUr");
		}
		else
		{
			strcpy(strChr, "unplaced");
		}
	}

	/* mouse */
	else if(strcmp(strLine, "MOUSE") == 0)
	{
		if( (nChr >= 1) && (nChr <= 19) )
		{
			sprintf(strChr, "chr%d", nChr);
		}
		else if(nChr == 20)
		{ 
			strcpy(strChr, "chrX");
		}
		else if(nChr == 21)
		{
			strcpy(strChr, "chrY");
		}
		else
		{
			strcpy(strChr, "unplaced");
		}
	}

	/* rat */
	else if(strcmp(strLine, "RAT") == 0)
	{
		if( (nChr >= 1) && (nChr <= 20) )
		{
			sprintf(strChr, "chr%d", nChr);
		}
		else if(nChr == 21)
		{ 
			strcpy(strChr, "chrX");
		}
		else if(nChr == 22)
		{
			strcpy(strChr, "chrM");
		}
		else if(nChr == 23)
		{
			strcpy(strChr, "chrUn");
		}
		else
		{
			strcpy(strChr, "unplaced");
		}
	}

	/* dog */
	else if(strcmp(strLine, "DOG") == 0)
	{
		if( (nChr >= 1) && (nChr <= 38) )
		{
			sprintf(strChr, "chr%d", nChr);
		}
		else if(nChr == 39)
		{ 
			strcpy(strChr, "chrX");
		}
		else if(nChr == 40)
		{
			strcpy(strChr, "chrM");
		}
		else if(nChr == 41)
		{
			strcpy(strChr, "chrUn");
		}
		else
		{
			strcpy(strChr, "unplaced");
		}
	}

	/* cow */
	else if(strcmp(strLine, "COW") == 0)
	{
		if( (nChr >= 1) && (nChr <= 29) )
		{
			sprintf(strChr, "chr%d", nChr);
		}
		else if(nChr == 30)
		{ 
			strcpy(strChr, "chrX");
		}
		else if(nChr == 31)
		{
			strcpy(strChr, "chrM");
		}
		else
		{
			strcpy(strChr, "unplaced");
		}
	}


	/* chicken */
	else if(strcmp(strLine, "CHICKEN") == 0)
	{
		if( (nChr >= 1) && (nChr <= 28) )
		{
			sprintf(strChr, "chr%d", nChr);
		}
		else if(nChr == 32)
		{
			sprintf(strChr, "chr%d", nChr);
		}
		else if(nChr == 33)
		{ 
			strcpy(strChr, "chrW");
		}
		else if(nChr == 34)
		{
			strcpy(strChr, "chrZ");
		}
		else if(nChr == 29)
		{
			strcpy(strChr, "chrE22C19W28_E50C23");
		}
		else if(nChr == 30)
		{
			strcpy(strChr, "chrE64");
		}
		else if(nChr == 31)
		{
			strcpy(strChr, "chrM");
		}
		else
		{
			strcpy(strChr, "unplaced");
		}
	}

	/* zebrafish */
	else if(strcmp(strLine, "ZEBRAFISH") == 0)
	{
		if( (nChr >= 1) && (nChr <= 25) )
		{
			sprintf(strChr, "chr%d", nChr);
		}
		else if(nChr == 26)
		{
			strcpy(strChr, "chrM");
		}
		else
		{
			strcpy(strChr, "unplaced");
		}
	}

	/* fly */
	else if(strcmp(strLine, "D_MELANOGASTER") == 0)
	{
		if( nChr == 1 )
		{
			sprintf(strChr, "chr2L");
		}
		else if( nChr == 2 )
		{
			sprintf(strChr, "chr2LHet");
		}
		else if( nChr == 3 )
		{
			sprintf(strChr, "chr2R");
		}
		else if( nChr == 4 )
		{
			sprintf(strChr, "chr2RHet");
		}
		else if( nChr == 5 )
		{
			sprintf(strChr, "chr3L");
		}
		else if( nChr == 6 )
		{
			sprintf(strChr, "chr3LHet");
		}
		else if( nChr == 7 )
		{
			sprintf(strChr, "chr3R");
		}
		else if( nChr == 8 )
		{
			sprintf(strChr, "chr3RHet");
		}
		else if( nChr == 9 )
		{
			sprintf(strChr, "chr4");
		}
		else if( nChr == 10 )
		{
			sprintf(strChr, "chrX");
		}
		else if( nChr == 11 )
		{
			sprintf(strChr, "chrXHet");
		}
		else if( nChr == 12 )
		{
			sprintf(strChr, "chrYHet");
		}
		else if( nChr == 13 )
		{
			sprintf(strChr, "chrU");
		}
		else if( nChr == 14 )
		{
			sprintf(strChr, "chrUextra");
		}
		else if( nChr == 15 )
		{
			sprintf(strChr, "chrM");
		}
		else
		{
			strcpy(strChr, "unplaced");
		}
	}

	/* C elegans */
	else if(strcmp(strLine, "C_ELEGANS") ==0)
	{
		if( nChr == 1 )
		{
			sprintf(strChr, "chrI");
		}
		else if( nChr == 2 )
		{
			sprintf(strChr, "chrII");
		}
		else if( nChr == 3 )
		{
			sprintf(strChr, "chrIII");
		}
		else if( nChr == 4 )
		{
			sprintf(strChr, "chrIV");
		}
		else if( nChr == 5 )
		{
			sprintf(strChr, "chrV");
		}
		else if( nChr == 6 )
		{
			sprintf(strChr, "chrX");
		}
		else if( nChr == 7 )
		{
			sprintf(strChr, "chrM");
		}
		else
		{
			strcpy(strChr, "unplaced");
		}
	}

	/* Yeast */
	else if(strcmp(strLine, "YEAST") == 0)
	{
		if( (nChr >= 1) && (nChr <= 16) )
		{
			sprintf(strChr, "chr%d", nChr);
		}
		else if(nChr == 17)
		{
			strcpy(strChr, "chrM");
		}
		else
		{
			strcpy(strChr, "unplaced");
		}
	}

	/* Arabidopsis */
	else if(strcmp(strLine, "ARABIDOPSIS") == 0)
	{
		if( (nChr >= 1) && (nChr <= 5) )
		{
			sprintf(strChr, "chr%d", nChr);
		}
		else if(nChr == 6)
		{
			strcpy(strChr, "chrC");
		}
		else if(nChr == 7)
		{
			strcpy(strChr, "chrM");
		}
		else
		{
			strcpy(strChr, "unplaced");
		}
	}

	/* bee */
	else if(strcmp(strLine, "BEE") == 0)
	{
		if( (nChr >= 1) && (nChr <= 9260) )
		{
			sprintf(strChr, "LG%d", nChr);
		}
		/* else if(nChr == 17)
		{
			strcpy(strChr, "LGUn");
		}
		else if(nChr == 18)
		{
			strcpy(strChr, "LGMT");
		} */
		else
		{
			strcpy(strChr, "unplaced");
		}
	}

	/* others */
	else
	{
		strcpy(strChr, "unplaced");
	}
	
	strcpy(strChrName, strChr);

	/* return */
	return strChrName;
}

/* ----------------------------------------------------------------------- */ 
/*  Genome_ChromosomeName_To_Index:                                        */
/*  convert chromsome name to numeric chromosome id.                       */
/* ----------------------------------------------------------------------- */ 
int Genome_ChromosomeName_To_Index(char strChrName[], char strSpecies[])
{
	/* define */
	int nChrid;
	char strLine[LINE_LENGTH];
	char strChr[LINE_LENGTH];
	char *chp;

	/* convert */
	strcpy(strLine, strSpecies);
	StrMakeUpper(strLine);
	strcpy(strChr, strChrName);
	StrMakeUpper(strChr);
	StrTrimLeft(strChr);
	StrTrimRight(strChr);

	/* human */
	if(strcmp(strLine, "HUMAN") == 0)
	{
		if(strstr(strChr, "CHR") != strChr)
			return -1;
	
		if( (strstr(strChr, "_") != NULL) || (strstr(strChr, "RANDOM") != NULL) 
			|| (strstr(strChr, "UNPLACED") != NULL) )
			return -1;

		chp = strChr+3;
		if(*chp == 'X')
		{
			nChrid = 23;
		}
		else if(*chp == 'Y')
		{
			nChrid = 24;
		}
		else if(*chp == 'M')
		{
			nChrid = 25;
		}
		else if((*chp >= '1') && (*chp <= '9'))
		{
			nChrid = atoi(chp);
		}
		else
		{
			nChrid = -1;
		}
	}

	/* rhesus monkey */
	else if(strcmp(strLine, "RHESUS") == 0)
	{
		if(strstr(strChr, "CHR") != strChr)
			return -1;
	
		if( (strstr(strChr, "_") != NULL) || (strstr(strChr, "RANDOM") != NULL) 
			|| (strstr(strChr, "UNPLACED") != NULL) )
			return -1;

		chp = strChr+3;
		if(*chp == 'X')
		{
			nChrid = 21;
		}
		else if(strcmp(chp, "UR") == 0)
		{
			nChrid = 22;
		}
		else if((*chp >= '1') && (*chp <= '9'))
		{
			nChrid = atoi(chp);
		}
		else
		{
			nChrid = -1;
		}
	}

	/* mouse */
	else if(strcmp(strLine, "MOUSE") == 0)
	{
		if(strstr(strChr, "CHR") != strChr)
			return -1;
	
		if( (strstr(strChr, "_") != NULL) || (strstr(strChr, "RANDOM") != NULL) 
			|| (strstr(strChr, "UNPLACED") != NULL) )
			return -1;

		chp = strChr+3;
		if(*chp == 'X')
		{
			nChrid = 20;
		}
		else if(*chp == 'Y')
		{
			nChrid = 21;
		}
		else if((*chp >= '1') && (*chp <= '9'))
		{
			nChrid = atoi(chp);
		}
		else
		{
			nChrid = -1;
		}
	}

	/* rat */
	else if(strcmp(strLine, "RAT") == 0)
	{
		if(strstr(strChr, "CHR") != strChr)
			return -1;
	
		if( (strstr(strChr, "_") != NULL) || (strstr(strChr, "RANDOM") != NULL) 
			|| (strstr(strChr, "UNPLACED") != NULL) )
			return -1;

		chp = strChr+3;
		if(*chp == 'X')
		{
			nChrid = 21;
		}
		else if(*chp == 'M')
		{
			nChrid = 22;
		}
		else if(strcmp(chp, "UN") == 0)
		{
			nChrid = 23;
		}
		else if((*chp >= '1') && (*chp <= '9'))
		{
			nChrid = atoi(chp);
		}
		else
		{
			nChrid = -1;
		}
	}

	/* dog */
	else if(strcmp(strLine, "DOG") == 0)
	{
		if(strstr(strChr, "CHR") != strChr)
			return -1;
	
		if( (strstr(strChr, "_") != NULL) || (strstr(strChr, "RANDOM") != NULL) 
			|| (strstr(strChr, "UNPLACED") != NULL) )
			return -1;

		chp = strChr+3;
		if(*chp == 'X')
		{
			nChrid = 39;
		}
		else if(*chp == 'M')
		{
			nChrid = 40;
		}
		else if(strcmp(chp, "Un") == 0)
		{
			nChrid = 41;
		}
		else if((*chp >= '1') && (*chp <= '9'))
		{
			nChrid = atoi(chp);
		}
		else
		{
			nChrid = -1;
		}
	}


	/* cow */
	else if(strcmp(strLine, "COW") == 0)
	{
		if(strstr(strChr, "CHR") != strChr)
			return -1;
	
		if( (strstr(strChr, "_") != NULL) || (strstr(strChr, "RANDOM") != NULL) 
			|| (strstr(strChr, "UNPLACED") != NULL) )
			return -1;

		chp = strChr+3;
		if(*chp == 'X')
		{
			nChrid = 30;
		}
		else if(*chp == 'M')
		{
			nChrid = 31;
		}
		else if((*chp >= '1') && (*chp <= '9'))
		{
			nChrid = atoi(chp);
		}
		else
		{
			nChrid = -1;
		}
	}

	/* chicken */
	else if(strcmp(strLine, "CHICKEN") == 0)
	{
		if(strstr(strChr, "CHR") != strChr)
			return -1;
	
		if( (strstr(strChr, "_") != NULL) || (strstr(strChr, "RANDOM") != NULL) 
			|| (strstr(strChr, "UNPLACED") != NULL) )
			return -1;

		chp = strChr+3;
		if(*chp == 'W')
		{
			nChrid = 33;
		}
		else if(*chp == 'Z')
		{
			nChrid = 34;
		}
		else if(*chp == 'E')
		{
			if(strcmp(chp, "E22C19W28_E50C23") == 0)
			{
				nChrid = 29;
			}
			else if(strcmp(chp, "E64") == 0)
			{
				nChrid = 30;
			}
		}
		else if(*chp == 'M')
		{
			nChrid = 31;
		}
		else if((*chp >= '1') && (*chp <= '9'))
		{
			nChrid = atoi(chp);
		}
		else
		{
			nChrid = -1;
		}

		/*
		if(*chp == 'W')
		{
			nChrid = 39;
		}
		else if(*chp == 'Z')
		{
			nChrid = 40;
		}
		else if(*chp == 'E')
		{
			if(strcmp(chp, "E22C19W28") == 0)
			{
				nChrid = 41;
			}
			else if(strcmp(chp, "E64") == 0)
			{
				nChrid = 42;
			}
			else if(strcmp(chp, "E26C13") == 0)
			{
				nChrid = 43;
			}
			else if(strcmp(chp, "E50C23") == 0)
			{
				nChrid = 44;
			}

		}
		else if(*chp == 'M')
		{
			nChrid = 45;
		}
		else if((*chp >= '1') && (*chp <= '9'))
		{
			nChrid = atoi(chp);
		}
		else
		{
			nChrid = -1;
		}
		*/
	}

	/* zebrafish */
	else if(strcmp(strLine, "ZEBRAFISH") == 0)
	{
		if(strstr(strChr, "CHR") != strChr)
			return -1;
	
		if( (strstr(strChr, "_") != NULL) || (strstr(strChr, "RANDOM") != NULL) 
			|| (strstr(strChr, "UNPLACED") != NULL) )
			return -1;

		chp = strChr+3;
		if(*chp == 'M')
		{
			nChrid = 26;
		}
		else if((*chp >= '1') && (*chp <= '9'))
		{
			nChrid = atoi(chp);
		}
		else
		{
			nChrid = -1;
		}
	}

	/* fly */
	else if(strcmp(strLine, "D_MELANOGASTER") == 0)
	{
		if(strstr(strChr, "CHR") != strChr)
			return -1;
	
		if( (strstr(strChr, "_") != NULL) || (strstr(strChr, "RANDOM") != NULL) 
			|| (strstr(strChr, "UNPLACED") != NULL) )
			return -1;

		if(strcmp(strChr, "CHR2L") == 0)
		{
			nChrid = 1;
		}
		else if(strcmp(strChr, "CHR2LHET") == 0)
		{
			nChrid = 2;
		}
		else if(strcmp(strChr, "CHR2R") == 0)
		{
			nChrid = 3;
		}
		else if(strcmp(strChr, "CHR2RHET") == 0)
		{
			nChrid = 4;
		}
		else if(strcmp(strChr, "CHR3L") == 0)
		{
			nChrid = 5;
		}
		else if(strcmp(strChr, "CHR3LHET") == 0)
		{
			nChrid = 6;
		}
		else if(strcmp(strChr, "CHR3R") == 0)
		{
			nChrid = 7;
		}
		else if(strcmp(strChr, "CHR3RHET") == 0)
		{
			nChrid = 8;
		}
		else if(strcmp(strChr, "CHR4") == 0)
		{
			nChrid = 9;
		}
		else if(strcmp(strChr, "CHRX") == 0)
		{
			nChrid = 10;
		}
		else if(strcmp(strChr, "CHRXHET") == 0)
		{
			nChrid = 11;
		}
		else if(strcmp(strChr, "CHRYHET") == 0)
		{
			nChrid = 12;
		}
		else if(strcmp(strChr, "CHRU") == 0)
		{
			nChrid = 13;
		}
		else if(strcmp(strChr, "CHRUEXTRA") == 0)
		{
			nChrid = 14;
		}
		else if(strcmp(strChr, "CHRM") == 0)
		{
			nChrid = 15;
		}
		else
		{
			nChrid = -1;
		}
	}

	/* C elegans */
	else if(strcmp(strLine, "C_ELEGANS") ==0)
	{
		if(strstr(strChr, "CHR") != strChr)
			return -1;
	
		if( (strstr(strChr, "_") != NULL) || (strstr(strChr, "RANDOM") != NULL) 
			|| (strstr(strChr, "UNPLACED") != NULL) )
			return -1;

		chp = strChr+3;

		if(strcmp(strChr, "CHRI") == 0)
		{
			nChrid = 1;
		}
		else if(strcmp(strChr, "CHRII") == 0)
		{
			nChrid = 2;
		}
		else if(strcmp(strChr, "CHRIII") == 0)
		{
			nChrid = 3;
		}
		else if(strcmp(strChr, "CHRIV") == 0)
		{
			nChrid = 4;
		}
		else if(strcmp(strChr, "CHRV") == 0)
		{
			nChrid = 5;
		}
		else if(strcmp(strChr, "CHRX") == 0)
		{
			nChrid = 6;
		}
		else if(strcmp(strChr, "CHRM") == 0)
		{
			nChrid = 7;
		}
		else
		{
			nChrid = -1;
		}
	}

	/* yeast */
	else if(strcmp(strLine, "YEAST") == 0)
	{
		if(strstr(strChr, "CHR") != strChr)
			return -1;
	
		if( (strstr(strChr, "_") != NULL) || (strstr(strChr, "RANDOM") != NULL) 
			|| (strstr(strChr, "UNPLACED") != NULL) )
			return -1;

		chp = strChr+3;
		if(*chp == 'M')
		{
			nChrid = 17;
		}
		else if((*chp >= '1') && (*chp <= '9'))
		{
			nChrid = atoi(chp);
		}
		else
		{
			nChrid = -1;
		}
	}

	/*  Arabidopsis */
	else if(strcmp(strLine, "ARABIDOPSIS") == 0)
	{
		if(strstr(strChr, "CHR") != strChr)
			return -1;
	
		if( (strstr(strChr, "_") != NULL) || (strstr(strChr, "RANDOM") != NULL) 
			|| (strstr(strChr, "UNPLACED") != NULL) )
			return -1;

		chp = strChr+3;
		if(*chp == 'M')
		{
			nChrid = 7;
		}
		else if(*chp == 'C')
		{
			nChrid = 6;
		}
		else if((*chp >= '1') && (*chp <= '9'))
		{
			nChrid = atoi(chp);
		}
		else
		{
			nChrid = -1;
		}
	}

	/* bee */
	else if(strcmp(strLine, "BEE") == 0)
	{
		if(strstr(strChr, "LG") != strChr)
			return -1;
	
		if( (strstr(strChr, "_") != NULL) || (strstr(strChr, "RANDOM") != NULL) 
			|| (strstr(strChr, "UNPLACED") != NULL) || (strstr(strChr, "|") != NULL) )
			return -1;

		chp = strChr+2;
		/* if(strstr(chp, "MT") != NULL)
		{
			nChrid = 18;
		}
		else if(strstr(chp, "UN") == chp)
		{
			nChrid = 17;
		} 
		else if((*chp >= '1') && (*chp <= '9'))
		{
			nChrid = atoi(chp);
		} */
		if((*chp >= '1') && (*chp <= '9'))
		{
			nChrid = atoi(chp);
		}
		else
		{
			nChrid = -1;
		}
	}

	/* others */
	else
	{
		nChrid = -1;
	}
	
	
	/* return */
	return nChrid;
}

/* ----------------------------------------------------------------------- */ 
/*  Genome_SpeciesAbbr_To_Species:                                         */
/*  convert species abbreviation name to species name.                     */
/* ----------------------------------------------------------------------- */ 
int Genome_SpeciesAbbr_To_Species(char strAbbr[], char strSpecies[])
{
	/* define */
	char strAbbrevation[LINE_LENGTH];

	strcpy(strAbbrevation, strAbbr);
	StrMakeLower(strAbbrevation);

	if(strcmp(strAbbrevation, "hg") == 0)
	{
		strcpy(strSpecies, "human");
	}
	else if(strcmp(strAbbrevation, "mm") == 0)
	{
		strcpy(strSpecies, "mouse");
	}
	else if(strcmp(strAbbrevation, "rn") == 0)
	{
		strcpy(strSpecies, "rat");
	}
	else if(strcmp(strAbbrevation, "canfam") == 0)
	{
		strcpy(strSpecies, "dog");
	}
	else if(strcmp(strAbbrevation, "bostau") == 0)
	{
		strcpy(strSpecies, "cow");
	}
	else if(strcmp(strAbbrevation, "mondom") == 0)
	{
		strcpy(strSpecies, "opossum");
	}
	else if(strcmp(strAbbrevation, "galgal") == 0)
	{
		strcpy(strSpecies, "chicken");
	}
	else if(strcmp(strAbbrevation, "xentro") == 0)
	{
		strcpy(strSpecies, "xenopus");
	}
	else if(strcmp(strAbbrevation, "danrer") == 0)
	{
		strcpy(strSpecies, "zebrafish");
	}
	else if(strcmp(strAbbrevation, "tetnig") == 0)
	{
		strcpy(strSpecies, "tetraodon");
	}
	else if(strcmp(strAbbrevation, "ce") == 0)
	{
		strcpy(strSpecies, "celegans");
	}
	else if(strcmp(strAbbrevation, "saccer") == 0)
	{
		strcpy(strSpecies, "yeast");
	}
	else
	{
		strcpy(strSpecies, "NA");
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Genome_Species_To_SpeciesAbbr:                                         */
/*  convert species name to species abbreviation name.                     */
/* ----------------------------------------------------------------------- */ 
int Genome_Species_To_SpeciesAbbr(char strSpecies[], char strAbbr[])
{
	/* define */
	char strSpeciesName[LINE_LENGTH];
	char strAbbrevation[LINE_LENGTH];

	strcpy(strSpeciesName, strSpecies);
	StrMakeLower(strSpeciesName);

    if(strcmp(strSpeciesName, "human") == 0)
	{
		strcpy(strAbbrevation, "hg");
	}
	else if(strcmp(strSpeciesName, "mouse") == 0)
	{
		strcpy(strAbbrevation, "mm");
	}
	else if(strcmp(strSpeciesName, "rat") == 0)
	{
		strcpy(strAbbrevation, "rn");
	}
	else if(strcmp(strSpeciesName, "dog") == 0)
	{
		strcpy(strAbbrevation, "canFam");
	}
	else if(strcmp(strSpeciesName, "cow") == 0)
	{
		strcpy(strAbbrevation, "bosTau");
	}
	else if(strcmp(strSpeciesName, "opossum") == 0)
	{
		strcpy(strAbbrevation, "monDom");
	}
	else if(strcmp(strSpeciesName, "chicken") == 0)
	{
		strcpy(strAbbrevation, "galGal");
	}
	else if(strcmp(strSpeciesName, "xenopus") == 0)
	{
		strcpy(strAbbrevation, "xenTro");
	}
	else if(strcmp(strSpeciesName, "zebrafish") == 0)
	{
		strcpy(strAbbrevation, "danRer");
	}
	else if(strcmp(strSpeciesName, "tetraodon") == 0)
	{
		strcpy(strAbbrevation, "tetNig");
	}
	else if(strcmp(strSpeciesName, "celegans") == 0)
	{
		strcpy(strAbbrevation, "ce");
	}
	else if(strcmp(strSpeciesName, "yeast") == 0)
	{
		strcpy(strAbbrevation, "sacCer");
	}
	else
	{
		strcpy(strAbbrevation, "NA");
	}

	strcpy(strAbbr, strAbbrevation);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Genome_MafAlign_To_BG_Main:                                            */
/*  Convert maf alignments to background variation.                        */
/*  return number of pairs of species compared.                            */
/* ----------------------------------------------------------------------- */ 
int Genome_MafAlign_To_BG_Main(char strInfoPath[])
{
	/* define */
	int nCount;
	int ni,nj,nk;
	char *chp;

	/* files */
	FILE *fpInfo;
	char strLine[MED_LINE_LENGTH];
	char strAlnPath[LINE_LENGTH];
	char strAlnFile[LINE_LENGTH];
	char strOutPath[LINE_LENGTH];
	char strInFile[LINE_LENGTH];
	char strOutFile[LINE_LENGTH];
	char strRefSpecies[LINE_LENGTH];
	char strAlnSuffix[LINE_LENGTH];
	char strChrName[LINE_LENGTH];
	char strTemp[LINE_LENGTH];

	char strLabel[LINE_LENGTH];

	/* number of species */
	int nSpeciesNum;

	/* species id */
	struct tagString **vTag;

	/* reference species */
	int nRefid;
	int nChrid;
	int nChrSize;

	/* chromosome size */
	struct INTMATRIX *pChrSize;

	/* comparisons */
	int nCompNum;
	struct INTMATRIX *pCompPair;

	/* bg param */
	int nWindowSize;
	int nStepSize;


	/* init */
	fpInfo = NULL;
	fpInfo = fopen(strInfoPath, "rt");
	if(fpInfo == NULL)
	{
		printf("Error: Genome_MafAlign_To_BG_Main, cannot open parameter info file!\n");
		return 0;
	}

	nCount = 0;

	/* load parameters */
	while(fgets(strLine, MED_LINE_LENGTH, fpInfo) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		
		sscanf(strLine, "%s ", strLabel);

		/* aln path */
		if(strcmp(strLabel, "Alignment_Dir") == 0)
		{
			chp = strLine+13;
			StrTrimLeft(chp);
			strcpy(strAlnPath, chp);
			AdjustDirectoryPath(strAlnPath);
		}
		/* species number */
		else if(strcmp(strLabel, "Species_Num") == 0)
		{
			sscanf(strLine, "%s %d", strLabel, &nSpeciesNum);
			if(nSpeciesNum <= 1)
			{
				printf("Error: the species number should be >= 2!\n");
				return 0;
			}

			vTag = (struct tagString **)calloc(nSpeciesNum, sizeof(struct tagString *));
		}
		/* species tag */
		else if(strcmp(strLabel, "ID") == 0)
		{
			sscanf(strLine, "%s %d %s", strLabel, &ni, strTemp);
			if(ni >= nSpeciesNum)
			{
				printf("Error: the species id should be < species number!\n");
				return 0;
			}

			StringAddTail((vTag+ni), strTemp);
		}
		/* reference */
		else if(strcmp(strLabel, "Reference") == 0)
		{
			sscanf(strLine, "%s %d", strLabel, &nRefid);
			if((nRefid < 0) || (nRefid >=nSpeciesNum))
			{
				printf("Error: reference id out of range");
				return 0;
			}
		}
		/* reference species */
		else if(strcmp(strLabel, "Ref_Species") == 0)
		{
			sscanf(strLine, "%s %s", strLabel, strRefSpecies);
		}
		/* reference chrlen */
		else if(strcmp(strLabel, "Ref_Chr_Size") == 0)
		{
			chp = strLine+12;
			StrTrimLeft(chp);
			strcpy(strInFile, chp);
			
			pChrSize = NULL;
			pChrSize = IMLOAD(strInFile);
			if(pChrSize == NULL)
			{
				printf("Error: cannot load chromosome length!\n");
				return 0;
			}
		}
		/* comparison number */
		else if(strcmp(strLabel, "Comparisons_Num") == 0)
		{
			sscanf(strLine, "%s %d", strLabel, &nCompNum);
			if(nCompNum < 1)
			{
				printf("Error: the comparison number should be >= 1!\n");
				return 0;
			}
			pCompPair = NULL;
			pCompPair = CreateIntMatrix(nCompNum, 2);
			if(pCompPair == NULL)
			{
				printf("Error: cannot create comparison storage space!\n");
				return 0;
			}
			nk = 0;
		}
		/* comparison id */
		else if(strcmp(strLabel, "COMP") == 0)
		{
			sscanf(strLine, "%s %d %d", strLabel, &ni, &nj);
			if(nk >= nCompNum)
			{
				printf("Error: the comparison number out of range!\n");
				return 0;
			}
			IMSETAT(pCompPair, nk, 0, ni);
			IMSETAT(pCompPair, nk, 1, nj);
			nk++;
		}
		/* window size */
		else if(strcmp(strLabel, "BG_Window_Size") == 0)
		{
			sscanf(strLine, "%s %d", strLabel, &nWindowSize);
			if(nWindowSize <=0 )
			{
				printf("Error: the window size should be >= 1!\n");
				return 0;
			}
		}
		/* step size */
		else if(strcmp(strLabel, "BG_Sliding_Step") == 0)
		{
			sscanf(strLine, "%s %d", strLabel, &nStepSize);
			if((nStepSize <=0) || (nStepSize > nWindowSize))
			{
				printf("Error: the step size should be >= 1 and <= window size!\n");
				return 0;
			}
		}
		/* output path */
		else if(strcmp(strLabel, "Output_Path") == 0)
		{
			chp = strLine+11;
			StrTrimLeft(chp);
			strcpy(strOutPath, chp);
			AdjustDirectoryPath(strOutPath);
		}
		/* alignment suffix */
		else if(strcmp(strLabel, "Alignment_Suf") == 0)
		{
			sscanf(strLine, "%s %s", strLabel, strAlnSuffix);
		}
		/* file by file */
		else if(strcmp(strLabel, "Alignment_File") == 0)
		{
			if(nk != nCompNum)
			{
				printf("Error: the comparison number not match!\n");
				return 0;
			}

			sscanf(strLine, "%s %s", strLabel, strChrName);

			nChrid = Genome_ChromosomeName_To_Index(strChrName, strRefSpecies);

			if( (nChrid<=0) || (nChrid>pChrSize->nHeight) )
			{
				printf("Error: %s is not supported in reference genome!\n", strChrName);
				return 0;
			}
			nChrid = nChrid-1;
			nChrSize = IMGETAT(pChrSize, nChrid, 0);

			sprintf(strAlnFile, "%s%s%s", strAlnPath, strChrName, strAlnSuffix);
			sprintf(strOutFile, "%s%s", strOutPath, strChrName);

			/* calculate */
			Genome_MafAlign_To_BG_Chr(strAlnFile, strOutFile, 
				nSpeciesNum, vTag, 
				nRefid, nChrSize, 
				nCompNum, pCompPair,
				nWindowSize, nStepSize);
		}
		/* do nothing */
		else
		{
		}
	}

	/* clear memeory */
	fclose(fpInfo);

	for(ni=0; ni<nSpeciesNum; ni++)
	{
		if(vTag[ni] != NULL)
		{
			DeleteString(vTag[ni]);
			vTag[ni] = NULL;
		}
	}
	free(vTag);

	DestroyIntMatrix(pCompPair);
	DestroyIntMatrix(pChrSize);

	/* return */
	return nCount;
}


/* ----------------------------------------------------------------------- */ 
/*  Genome_MafAlign_To_BG_Chr:                                             */
/*  Convert maf alignments to background variation.                        */
/*  return PROC_SUCCESS if successful.                                     */
/* ----------------------------------------------------------------------- */ 
int Genome_MafAlign_To_BG_Chr(char strAlnFile[], char strOutFile[], 
				int nSpeciesNum, struct tagString **vTag, 
				int nRefid, int nChrSize, 
				int nCompNum, struct INTMATRIX *pCompPair,
				int nWindowSize, int nStepSize)
{
	/* define */

	/* files */
	char strOutName[LINE_LENGTH];
	FILE *fpIn;
	FILE *fpOut;
	int numwritten;

	/* window size */
	int nBinNum;
	int nResLen;
	int ni,nj,nk,nl,nScanLen;
	int s1,s2;
	

	/* string */
	char strAlnLine[ALN_LENGTH];
	struct tagString **vAln;
	char strTag[LINE_LENGTH];
	char strLabel[LINE_LENGTH];

	/* background */
	struct INTMATRIX *pPosInfo;
	struct BYTEMATRIX *pGapInfo;

	struct INTMATRIX **vCompTot;
	struct INTMATRIX **vCompIdn;
	struct INTMATRIX **vLenTot;
	struct INTMATRIX **vLenAln;
	struct INTMATRIX *pWinCompTot;
	struct INTMATRIX *pWinCompIdn;
	struct INTMATRIX *pWinLenTot;
	struct INTMATRIX *pWinLenAln;
	struct DOUBLEMATRIX *pIdentity;
	struct DOUBLEMATRIX *pAlignibility;

	/* loading variables */
	int nWinBinRatio;
	int nHalfWid;
	int npos,nPos,nChrLen,nRc;
	char chRc;
	int nBinindex;

	/* init */
	nBinNum = (int)(nChrSize/nStepSize);
	if(nChrSize%nStepSize != 0)
		nBinNum++;
	pIdentity = NULL;
	pAlignibility = NULL;
	pWinCompTot = NULL;
	pWinCompIdn = NULL;
	pWinLenTot = NULL;
	pWinLenAln = NULL;
	
	/* create storage space */
	pPosInfo = NULL;
	pPosInfo = CreateIntMatrix(nSpeciesNum, 5); /* col1: aln exist?; col2: pos; col3: Pos; col4: rc; col5: length */
	if(pPosInfo == NULL)
	{
		printf("Error: cannot create memory for background calculation!\n");
		exit(EXIT_FAILURE);
	}

	pGapInfo = NULL;
	pGapInfo = CreateByteMatrix(nCompNum, 2);
	if(pGapInfo == NULL)
	{
		printf("Error: cannot create memory for background calculation!\n");
		exit(EXIT_FAILURE);
	}

	vAln = NULL;
	vAln = (struct tagString**)calloc(nSpeciesNum, sizeof(struct tagString *));
	if(vAln == NULL)
	{
		printf("Error: cannot create memory for background calculation!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nSpeciesNum; ni++)
	{
		vAln[ni] = CreateString(ALN_LENGTH);
		if(vAln[ni] == NULL)
		{
			printf("Error: cannot create memory for background calculation!\n");
			exit(EXIT_FAILURE);
		}
	}

	vCompTot = NULL;
	vCompTot = (struct INTMATRIX **)calloc(nCompNum, sizeof(struct INTMATRIX *));
	if(vCompTot == NULL)
	{
		printf("Error: cannot create memory for background calculation!\n");
		exit(EXIT_FAILURE);
	}
	vCompIdn = NULL;
	vCompIdn = (struct INTMATRIX **)calloc(nCompNum, sizeof(struct INTMATRIX *));
	if(vCompIdn == NULL)
	{
		printf("Error: cannot create memory for background calculation!\n");
		exit(EXIT_FAILURE);
	}
	vLenTot = NULL;
	vLenTot = (struct INTMATRIX **)calloc(nCompNum, sizeof(struct INTMATRIX *));
	if(vLenTot == NULL)
	{
		printf("Error: cannot create memory for background calculation!\n");
		exit(EXIT_FAILURE);
	}
	vLenAln = NULL;
	vLenAln = (struct INTMATRIX **)calloc(nCompNum, sizeof(struct INTMATRIX *));
	if(vLenAln == NULL)
	{
		printf("Error: cannot create memory for background calculation!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nCompNum; ni++)
	{
		vCompTot[ni] = CreateIntMatrix(nBinNum, 1);
		if(vCompTot[ni] == NULL)
		{
			printf("Error: cannot create memory for background calculation!\n");
			exit(EXIT_FAILURE);
		}
		vCompIdn[ni] = CreateIntMatrix(nBinNum, 1);
		if(vCompIdn[ni] == NULL)
		{
			printf("Error: cannot create memory for background calculation!\n");
			exit(EXIT_FAILURE);
		}
		vLenTot[ni] = CreateIntMatrix(nBinNum, 1);
		if(vLenTot[ni] == NULL)
		{
			printf("Error: cannot create memory for background calculation!\n");
			exit(EXIT_FAILURE);
		}
		vLenAln[ni] = CreateIntMatrix(nBinNum, 1);
		if(vLenAln[ni] == NULL)
		{
			printf("Error: cannot create memory for background calculation!\n");
			exit(EXIT_FAILURE);
		}
	}


	/* read in alignment and processing */
	fpIn = NULL;
	fpIn = fopen(strAlnFile, "rt");
	if(fpIn == NULL)
	{
		printf("Error: cannot open alignment file!\n");
		return PROC_FAILURE;
	}

	while(fgets(strAlnLine, ALN_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strAlnLine);
		StrTrimRight(strAlnLine);
				
		/* ignore annotations */
		if(strAlnLine[0] == '#')
		{
			continue;
		}

		/* init */
		else if(strAlnLine[0] == 'a')
		{
			for(ni=0; ni<nSpeciesNum; ni++)
			{
				for(nj=0; nj<5; nj++)
				{
					IMSETAT(pPosInfo, ni, nj, 0);
				}
			}
			for(ni=0; ni<nCompNum; ni++)
			{
				for(nj=0; nj<2; nj++)
				{
					BMSETAT(pGapInfo, ni, nj, 0);
				}
			}
		}

		/* load alignment */
		else if(strAlnLine[0] == 's')
		{
			sscanf(strAlnLine, "%s %s", strLabel, strTag);
			for(ni=0; ni<nSpeciesNum; ni++)
			{
				if(strstr(strTag, vTag[ni]->m_pString) == strTag)
				{
					break;
				}
			}
			if(ni >= nSpeciesNum)
			{
				printf("Error: species not found!\n");
				exit(EXIT_FAILURE);
			}

			sscanf(strAlnLine, "%s %s %d %d %c %d %s",
				strLabel, strTag, &npos, &nPos, &chRc,
					&nChrLen, vAln[ni]->m_pString);
			
			StrMakeUpper(vAln[ni]->m_pString);

			nPos += npos;
			if(chRc == '-') 
				nRc = 1;
			else
				nRc = 0;

			IMSETAT(pPosInfo, ni, 0, 1);
			IMSETAT(pPosInfo, ni, 1, npos);
			IMSETAT(pPosInfo, ni, 2, nPos);
			IMSETAT(pPosInfo, ni, 3, nRc);
			IMSETAT(pPosInfo, ni, 4, nChrLen);
		}

		/* process */
		else if(strAlnLine[0] == '\0')
		{
			if(IMGETAT(pPosInfo, nRefid, 0) == 0)
				continue;

			nScanLen = strlen(vAln[nRefid]->m_pString);

			/* if + strand */
			if(IMGETAT(pPosInfo, nRefid, 3) == 0)
			{
				nl = IMGETAT(pPosInfo, nRefid, 1);
				
				/* process base by base */
				for(nj=0; nj<nScanLen; nj++)
				{
					nBinindex = nl/nStepSize;
					if(nBinindex == nBinNum)
						nBinindex--;
						
					/* pair by pair */
					for(ni=0; ni<nCompNum; ni++)
					{
						s1 = IMGETAT(pCompPair, ni, 0);
						s2 = IMGETAT(pCompPair, ni, 1);

						/* if we have alignment */
						if((IMGETAT(pPosInfo, s1, 0) == 1) && (IMGETAT(pPosInfo, s2, 0) == 1))
						{
							if(vAln[s2]->m_pString[nj] != '-')
								BMSETAT(pGapInfo, ni, 1, 0);

							/* s1 gap */
							if(vAln[s1]->m_pString[nj] == '-')
							{
								/* s1 gap vs s2 base */
								if(vAln[s2]->m_pString[nj] != '-')
								{
									/* s1 gap open */
									if(BMGETAT(pGapInfo, ni, 0) == 0)
									{
										BMSETAT(pGapInfo, ni, 0, 1);
									}
									/* s1 gap continue */
									else
									{
										
									}

									/* count */
									vCompTot[ni]->pMatElement[nBinindex] = vCompTot[ni]->pMatElement[nBinindex]+1;
								}

								/* s1 gap vs s2 gap */
								else
								{
								}
							}
							/* s1 base */
							else
							{
								BMSETAT(pGapInfo, ni, 0, 0);

								/* s1 base vs s2 gap */
								if(vAln[s2]->m_pString[nj] == '-')
								{
									/* s2 gap open */
									if(BMGETAT(pGapInfo, ni, 1) == 0)
									{
										BMSETAT(pGapInfo, ni, 1, 1);
									}
									/* s2 gap continue */
									else
									{
										
									}
								}

								/* s1 base vs s2 base */
								else
								{
									/* effective base */
									if((vAln[s1]->m_pString[nj] != 'N') && (vAln[s2]->m_pString[nj] != 'N'))
									{
										/* substitution */
										if(vAln[s1]->m_pString[nj] != vAln[s2]->m_pString[nj])
										{
										}
										/* identical */
										else
										{
											vCompIdn[ni]->pMatElement[nBinindex] = vCompIdn[ni]->pMatElement[nBinindex]+1;
										}
									}
								}

								/* count */
								vCompTot[ni]->pMatElement[nBinindex] = vCompTot[ni]->pMatElement[nBinindex]+1;
							}

							/* count */
							if(vAln[nRefid]->m_pString[nj] != '-')
							{
								vLenTot[ni]->pMatElement[nBinindex] = vLenTot[ni]->pMatElement[nBinindex]+1;
								vLenAln[ni]->pMatElement[nBinindex] = vLenAln[ni]->pMatElement[nBinindex]+1;
							}
						}

						/* if we don't have alignment */
						else
						{
							if(vAln[nRefid]->m_pString[nj] != '-')
							{
								vLenTot[ni]->pMatElement[nBinindex] = vLenTot[ni]->pMatElement[nBinindex]+1;
							}
						}
					}

					if(vAln[nRefid]->m_pString[nj] != '-')
					{
						nl++;
					}	
				}

				if(nl != IMGETAT(pPosInfo, nRefid, 2))
				{
					printf("Error: coordinate wrong!\n");
					exit(EXIT_FAILURE);
				}
			}

			/* if - strand */
			else
			{
				nl = IMGETAT(pPosInfo, nRefid, 4)-1-IMGETAT(pPosInfo, nRefid, 1);
				
				/* process base by base */
				for(nj=0; nj<nScanLen; nj++)
				{
					nBinindex = nl/nStepSize;
					if(nBinindex == nBinNum)
						nBinindex--;
						
					/* pair by pair */
					for(ni=0; ni<nCompNum; ni++)
					{
						s1 = IMGETAT(pCompPair, ni, 0);
						s2 = IMGETAT(pCompPair, ni, 1);

						/* if we have alignment */
						if((IMGETAT(pPosInfo, s1, 0) == 1) && (IMGETAT(pPosInfo, s2, 0) == 1))
						{
							if(vAln[s2]->m_pString[nj] != '-')
								BMSETAT(pGapInfo, ni, 1, 0);

							/* s1 gap */
							if(vAln[s1]->m_pString[nj] == '-')
							{
								/* s1 gap vs s2 base */
								if(vAln[s2]->m_pString[nj] != '-')
								{
									/* s1 gap open */
									if(BMGETAT(pGapInfo, ni, 0) == 0)
									{
										BMSETAT(pGapInfo, ni, 0, 1);
									}
									/* s1 gap continue */
									else
									{
										
									}

									/* count */
									vCompTot[ni]->pMatElement[nBinindex] = vCompTot[ni]->pMatElement[nBinindex]+1;
								}

								/* s1 gap vs s2 gap */
								else
								{
								}
							}
							/* s1 base */
							else
							{
								BMSETAT(pGapInfo, ni, 0, 0);

								/* s1 base vs s2 gap */
								if(vAln[s2]->m_pString[nj] == '-')
								{
									/* s2 gap open */
									if(BMGETAT(pGapInfo, ni, 1) == 0)
									{
										BMSETAT(pGapInfo, ni, 1, 1);
									}
									/* s2 gap continue */
									else
									{
										
									}
								}

								/* s1 base vs s2 base */
								else
								{
									/* effective base */
									if((vAln[s1]->m_pString[nj] != 'N') && (vAln[s2]->m_pString[nj] != 'N'))
									{
										/* substitution */
										if(vAln[s1]->m_pString[nj] != vAln[s2]->m_pString[nj])
										{
										}
										/* identical */
										else
										{
											vCompIdn[ni]->pMatElement[nBinindex] = vCompIdn[ni]->pMatElement[nBinindex]+1;
										}
									}
								}

								/* count */
								vCompTot[ni]->pMatElement[nBinindex] = vCompTot[ni]->pMatElement[nBinindex]+1;
							}

							/* count */
							if(vAln[nRefid]->m_pString[nj] != '-')
							{
								vLenTot[ni]->pMatElement[nBinindex] = vLenTot[ni]->pMatElement[nBinindex]+1;
								vLenAln[ni]->pMatElement[nBinindex] = vLenAln[ni]->pMatElement[nBinindex]+1;
							}
						}

						/* if we don't have alignment */
						else
						{
							if(vAln[nRefid]->m_pString[nj] != '-')
							{
								vLenTot[ni]->pMatElement[nBinindex] = vLenTot[ni]->pMatElement[nBinindex]+1;
							}
						}
					}

					if(vAln[nRefid]->m_pString[nj] != '-')
					{
						nl--;
					}	
				}

				if(nl != (IMGETAT(pPosInfo, nRefid, 4)-1-IMGETAT(pPosInfo, nRefid, 2)))
				{
					printf("Error: coordinate wrong!\n");
					exit(EXIT_FAILURE);
				}
			}
		}
	}

	fclose(fpIn);

	/* bin to win and output */
	nWinBinRatio = (int)(nWindowSize/nStepSize);
	nHalfWid = nWinBinRatio/2;
	nWinBinRatio = 2*nHalfWid;
	for(ni=0; ni<nCompNum; ni++)
	{
		/* prepare */
		pIdentity = NULL;
		pIdentity = CreateDoubleMatrix(nBinNum, 1);
		pAlignibility = NULL;
		pAlignibility = CreateDoubleMatrix(nBinNum, 1);
		pWinCompTot = NULL;
		pWinCompTot = CreateIntMatrix(nBinNum, 1);
		pWinCompIdn = NULL;
		pWinCompIdn = CreateIntMatrix(nBinNum, 1);
		pWinLenTot = NULL;
		pWinLenTot = CreateIntMatrix(nBinNum, 1);
		pWinLenAln = NULL;
		pWinLenAln = CreateIntMatrix(nBinNum, 1);

		/* bin to win */
		for(nj=0; nj<=nHalfWid; nj++)
		{
			pWinCompTot->pMatElement[0] = pWinCompTot->pMatElement[0]+vCompTot[ni]->pMatElement[nj];
			pWinCompIdn->pMatElement[0] = pWinCompIdn->pMatElement[0]+vCompIdn[ni]->pMatElement[nj];
			/* pWinLenTot->pMatElement[0] = pWinLenTot->pMatElement[0]+vLenTot[ni]->pMatElement[nj]; */
			pWinLenTot->pMatElement[0] = pWinLenTot->pMatElement[0]+nStepSize;
			pWinLenAln->pMatElement[0] = pWinLenAln->pMatElement[0]+vLenAln[ni]->pMatElement[nj];		
		}

		nk = 0;
		if(pWinCompTot->pMatElement[nk] != 0)
			pIdentity->pMatElement[nk] = (double)(pWinCompIdn->pMatElement[nk])/(double)(pWinCompTot->pMatElement[nk]);
		if(pWinLenTot->pMatElement[nk] != 0)
			pAlignibility->pMatElement[nk] = (double)(pWinLenAln->pMatElement[nk])/(double)(pWinLenTot->pMatElement[nk]);

		nk = 1;
		for(; nj<=nWinBinRatio; nj++)
		{
			pWinCompTot->pMatElement[nk] = pWinCompTot->pMatElement[nk-1]+vCompTot[ni]->pMatElement[nj];
			pWinCompIdn->pMatElement[nk] = pWinCompIdn->pMatElement[nk-1]+vCompIdn[ni]->pMatElement[nj];
			/* pWinLenTot->pMatElement[nk] = pWinLenTot->pMatElement[nk-1]+vLenTot[ni]->pMatElement[nj]; */
			pWinLenTot->pMatElement[nk] = pWinLenTot->pMatElement[nk-1]+nStepSize;
			pWinLenAln->pMatElement[nk] = pWinLenAln->pMatElement[nk-1]+vLenAln[ni]->pMatElement[nj];

			if(pWinCompTot->pMatElement[nk] != 0)
				pIdentity->pMatElement[nk] = (double)(pWinCompIdn->pMatElement[nk])/(double)(pWinCompTot->pMatElement[nk]);
			if(pWinLenTot->pMatElement[nk] != 0)
				pAlignibility->pMatElement[nk] = (double)(pWinLenAln->pMatElement[nk])/(double)(pWinLenTot->pMatElement[nk]);
			nk++;
		}

		nResLen = nChrSize%nStepSize;
		for(; nj<nBinNum; nj++)
		{
			pWinCompTot->pMatElement[nk] = pWinCompTot->pMatElement[nk-1]+vCompTot[ni]->pMatElement[nj]-vCompTot[ni]->pMatElement[nj-nWinBinRatio-1];
			pWinCompIdn->pMatElement[nk] = pWinCompIdn->pMatElement[nk-1]+vCompIdn[ni]->pMatElement[nj]-vCompIdn[ni]->pMatElement[nj-nWinBinRatio-1];
			/* pWinLenTot->pMatElement[nk] = pWinLenTot->pMatElement[nk-1]+vLenTot[ni]->pMatElement[nj]-vLenTot[ni]->pMatElement[nj-nWinBinRatio-1]; */
			if(nj == (nBinNum-1))
			{
				if(nResLen != 0)
					pWinLenTot->pMatElement[nk] = pWinLenTot->pMatElement[nk-1]+nResLen-nStepSize;
				else
					pWinLenTot->pMatElement[nk] = pWinLenTot->pMatElement[nk-1];
			}
			else
			{
				pWinLenTot->pMatElement[nk] = pWinLenTot->pMatElement[nk-1];
			}
			pWinLenAln->pMatElement[nk] = pWinLenAln->pMatElement[nk-1]+vLenAln[ni]->pMatElement[nj]-vLenAln[ni]->pMatElement[nj-nWinBinRatio-1];
			
			if(pWinCompTot->pMatElement[nk] != 0)
				pIdentity->pMatElement[nk] = (double)(pWinCompIdn->pMatElement[nk])/(double)(pWinCompTot->pMatElement[nk]);
			if(pWinLenTot->pMatElement[nk] != 0)
				pAlignibility->pMatElement[nk] = (double)(pWinLenAln->pMatElement[nk])/(double)(pWinLenTot->pMatElement[nk]);

			nk++;
		}

		for(; nj<(nBinNum+nHalfWid); nj++)
		{
			pWinCompTot->pMatElement[nk] = pWinCompTot->pMatElement[nk-1]-vCompTot[ni]->pMatElement[nj-nWinBinRatio-1];
			pWinCompIdn->pMatElement[nk] = pWinCompIdn->pMatElement[nk-1]-vCompIdn[ni]->pMatElement[nj-nWinBinRatio-1];
			/* pWinLenTot->pMatElement[nk] = pWinLenTot->pMatElement[nk-1]-vLenTot[ni]->pMatElement[nj-nWinBinRatio-1]; */
			pWinLenTot->pMatElement[nk] = pWinLenTot->pMatElement[nk-1]-nStepSize;
			pWinLenAln->pMatElement[nk] = pWinLenAln->pMatElement[nk-1]-vLenAln[ni]->pMatElement[nj-nWinBinRatio-1];
			
			if(pWinCompTot->pMatElement[nk] != 0)
				pIdentity->pMatElement[nk] = (double)(pWinCompIdn->pMatElement[nk])/(double)(pWinCompTot->pMatElement[nk]);
			if(pWinLenTot->pMatElement[nk] != 0)
				pAlignibility->pMatElement[nk] = (double)(pWinLenAln->pMatElement[nk])/(double)(pWinLenTot->pMatElement[nk]);

			nk++;
		}

		/* write */
		s1 = IMGETAT(pCompPair, ni, 0);
		s2 = IMGETAT(pCompPair, ni, 1);
		sprintf(strOutName, "%s_%s-%s.bg", strOutFile, vTag[s1]->m_pString, vTag[s2]->m_pString);
		/* DMSAVE(pIdentity, strOutName); */
		fpOut = NULL;
		fpOut = fopen(strOutName, "wb");
		if(fpOut == NULL)
		{
			printf("Error: Genome_MafAlign_To_BG_Chr, cannot open output file!\n");
			exit(EXIT_FAILURE);
		}
		numwritten = fwrite(pIdentity->pMatElement , sizeof(double), nBinNum, fpOut);
		if(numwritten != nBinNum)
		{
			printf("Error: write background wrong!\n");
			exit(EXIT_FAILURE);
		}
		fclose(fpOut);
		
		sprintf(strOutName, "%s_%s-%s.cov", strOutFile, vTag[s1]->m_pString, vTag[s2]->m_pString);
		/* DMSAVE(pAlignibility, strOutName); */
		fpOut = NULL;
		fpOut = fopen(strOutName, "wb");
		if(fpOut == NULL)
		{
			printf("Error: Genome_MafAlign_To_BG_Chr, cannot open output file!\n");
			exit(EXIT_FAILURE);
		}
		numwritten = fwrite(pAlignibility->pMatElement , sizeof(double), nBinNum, fpOut);
		if(numwritten != nBinNum)
		{
			printf("Error: write background wrong!\n");
			exit(EXIT_FAILURE);
		}
		fclose(fpOut);


		/* destroy */
		DestroyIntMatrix(pWinCompTot);
		DestroyIntMatrix(pWinCompIdn);
		DestroyIntMatrix(pWinLenTot);
		DestroyIntMatrix(pWinLenAln);
		DestroyDoubleMatrix(pIdentity);
		DestroyDoubleMatrix(pAlignibility);
	}
	

	/* release memory */
	DestroyIntMatrix(pPosInfo);
	DestroyByteMatrix(pGapInfo);

	for(ni=0; ni<nSpeciesNum; ni++)
	{
		DeleteString(vAln[ni]);
		vAln[ni] = NULL;
	}
	free(vAln);

	for(ni=0; ni<nCompNum; ni++)
	{
		DestroyIntMatrix(vCompTot[ni]);
		vCompTot[ni] = NULL;
		DestroyIntMatrix(vCompIdn[ni]);
		vCompIdn[ni] = NULL;
		DestroyIntMatrix(vLenTot[ni]);
		vLenTot[ni] = NULL;
		DestroyIntMatrix(vLenAln[ni]);
		vLenAln[ni] = NULL;
	}
	free(vCompTot);
	free(vCompIdn);
	free(vLenTot);
	free(vLenAln);

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  Genome_MafAlign_To_CS_Main:                                            */
/*  Convert maf alignments to conservation score.                          */
/*  return number of alignment files processed.                            */
/* ----------------------------------------------------------------------- */ 
int Genome_MafAlign_To_CS_Main(char strInfoPath[])
{
	/* define */
	int nCount;
	int ni,nj,nk;
	int nrlen;
	int numwritten;
	char *chp;

	/* files */
	FILE *fpInfo;
	char strLine[MED_LINE_LENGTH];
	char strAlnPath[LINE_LENGTH];
	char strAlnFile[LINE_LENGTH];
	char strBGPath[LINE_LENGTH];
	char strOutPath[LINE_LENGTH];
	char strInFile[LINE_LENGTH];
	char strOutFile[LINE_LENGTH];
	char strRefSpecies[LINE_LENGTH];
	char strTarSpecies[LINE_LENGTH];
	char strAlnSuffix[LINE_LENGTH];
	char strChrName[LINE_LENGTH];
	char strTemp[LINE_LENGTH];

	FILE **vfpScore;
	struct BYTEMATRIX *pInitScore;


	char strLabel[LINE_LENGTH];

	/* number of species */
	int nSpeciesNum;
	
	/* species id */
	struct tagString **vTag;

	/* reference species */
	int nRefid;
	int nChrid;
	int nChrSize;
	int nTarid;

	/* chromosome size */
	struct INTMATRIX *pChrSize;
	struct INTMATRIX *pTargetChrSize;

	/* comparisons */
	int nCompNum;
	struct INTMATRIX *pCompPair;

	/* bg param */
	int nWindowSize;
	int nStepSize;
	int nConserveWinSize;
	double dMinPcutoff;

	/* init */
	fpInfo = NULL;
	fpInfo = fopen(strInfoPath, "rt");
	if(fpInfo == NULL)
	{
		printf("Error: Genome_MafAlign_To_BG_Main, cannot open parameter info file!\n");
		return 0;
	}

	nCount = 0;

	/* load parameters */
	while(fgets(strLine, MED_LINE_LENGTH, fpInfo) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		
		sscanf(strLine, "%s ", strLabel);

		/* aln path */
		if(strcmp(strLabel, "Alignment_Dir") == 0)
		{
			chp = strLine+13;
			StrTrimLeft(chp);
			strcpy(strAlnPath, chp);
			AdjustDirectoryPath(strAlnPath);
		}
		/* species number */
		else if(strcmp(strLabel, "Species_Num") == 0)
		{
			sscanf(strLine, "%s %d", strLabel, &nSpeciesNum);
			if(nSpeciesNum <= 1)
			{
				printf("Error: the species number should be >= 2!\n");
				return 0;
			}

			vTag = (struct tagString **)calloc(nSpeciesNum, sizeof(struct tagString *));
		}
		/* species tag */
		else if(strcmp(strLabel, "ID") == 0)
		{
			sscanf(strLine, "%s %d %s", strLabel, &ni, strTemp);
			if(ni >= nSpeciesNum)
			{
				printf("Error: the species id should be < species number!\n");
				return 0;
			}

			StringAddTail((vTag+ni), strTemp);
		}
		/* reference */
		else if(strcmp(strLabel, "Reference") == 0)
		{
			sscanf(strLine, "%s %d", strLabel, &nRefid);
			if((nRefid < 0) || (nRefid >=nSpeciesNum))
			{
				printf("Error: reference id out of range");
				return 0;
			}
		}
		/* reference species */
		else if(strcmp(strLabel, "Ref_Species") == 0)
		{
			sscanf(strLine, "%s %s", strLabel, strRefSpecies);
		}
		/* reference chrlen */
		else if(strcmp(strLabel, "Ref_Chr_Size") == 0)
		{
			chp = strLine+12;
			StrTrimLeft(chp);
			strcpy(strInFile, chp);

			pChrSize = NULL;
			pChrSize = IMLOAD(strInFile);
			if(pChrSize == NULL)
			{
				printf("Error: cannot load chromosome length!\n");
				return 0;
			}
		}
		/* target */
		else if(strcmp(strLabel, "Target") == 0)
		{
			sscanf(strLine, "%s %d", strLabel, &nTarid);
			if((nTarid < 0) || (nTarid >=nSpeciesNum))
			{
				printf("Error: target id out of range");
				return 0;
			}
		}
		/* target species */
		else if(strcmp(strLabel, "Tar_Species") == 0)
		{
			sscanf(strLine, "%s %s", strLabel, strTarSpecies);
		}
		/* target chrlen */
		else if(strcmp(strLabel, "Tar_Chr_Size") == 0)
		{
			chp = strLine+12;
			StrTrimLeft(chp);
			strcpy(strInFile, chp);

			pTargetChrSize = NULL;
			pTargetChrSize = IMLOAD(strInFile);
			if(pTargetChrSize == NULL)
			{
				printf("Error: cannot load chromosome length!\n");
				return 0;
			}
		}
		/* comparison number */
		else if(strcmp(strLabel, "Comparisons_Num") == 0)
		{
			sscanf(strLine, "%s %d", strLabel, &nCompNum);
			if(nCompNum < 1)
			{
				printf("Error: the comparison number should be >= 1!\n");
				return 0;
			}
			pCompPair = NULL;
			pCompPair = CreateIntMatrix(nCompNum, 2);
			if(pCompPair == NULL)
			{
				printf("Error: cannot create comparison storage space!\n");
				return 0;
			}
			nk = 0;
		}
		/* comparison id */
		else if(strcmp(strLabel, "COMP") == 0)
		{
			sscanf(strLine, "%s %d %d", strLabel, &ni, &nj);
			if(nk >= nCompNum)
			{
				printf("Error: the comparison number out of range!\n");
				return 0;
			}
			IMSETAT(pCompPair, nk, 0, ni);
			IMSETAT(pCompPair, nk, 1, nj);
			nk++;
		}
		/* window size */
		else if(strcmp(strLabel, "BG_Window_Size") == 0)
		{
			sscanf(strLine, "%s %d", strLabel, &nWindowSize);
			if(nWindowSize <=0 )
			{
				printf("Error: the window size should be >= 1!\n");
				return 0;
			}
		}
		/* step size */
		else if(strcmp(strLabel, "BG_Sliding_Step") == 0)
		{
			sscanf(strLine, "%s %d", strLabel, &nStepSize);
			if((nStepSize <=0) || (nStepSize > nWindowSize))
			{
				printf("Error: the step size should be >= 1 and <= window size!\n");
				return 0;
			}
		}
		/* background path */
		else if(strcmp(strLabel, "BG_Path") == 0)
		{
			chp = strLine+7;
			StrTrimLeft(chp);
			strcpy(strBGPath, chp);
			AdjustDirectoryPath(strBGPath);
		}
		/* output path */
		else if(strcmp(strLabel, "Output_Path") == 0)
		{
			chp = strLine+11;
			StrTrimLeft(chp);
			strcpy(strOutPath, chp);
			AdjustDirectoryPath(strOutPath);

			vfpScore = NULL;
			vfpScore = (FILE **)calloc(pTargetChrSize->nHeight, sizeof(FILE *));
			pInitScore = NULL;
			pInitScore = CreateByteMatrix(1, nStepSize);
			for(ni=0; ni<pTargetChrSize->nHeight; ni++)
			/* for(ni=20; ni<21; ni++) */
			{
				Genome_Index_To_ChromosomeName(strChrName, strTarSpecies, (ni+1));

				if(strcmp(strChrName, "unplaced") == 0)
				{
					printf("Error: chr%d is not supported in the target genome!\n", (ni+1));
					exit(EXIT_FAILURE);
				}
				
				sprintf(strOutFile, "%s%s.cs", strOutPath, strChrName);
				
				vfpScore[ni] = fopen(strOutFile, "w+b");
				if(vfpScore[ni] == NULL)
				{
					printf("Error: Genome_MafAlign_To_CS_Main, cannot open output file!\n");
					exit(EXIT_FAILURE);
				}
				nrlen = pTargetChrSize->pMatElement[ni];
				numwritten = 0;
				while(nrlen >= nStepSize)
				{
					numwritten += fwrite(pInitScore->pMatElement, sizeof(unsigned char), nStepSize, vfpScore[ni]);
					nrlen -= nStepSize;
				}
				if(nrlen > 0)
				{
					numwritten += fwrite(pInitScore->pMatElement, sizeof(unsigned char), nrlen, vfpScore[ni]);
				}

				if(numwritten != pTargetChrSize->pMatElement[ni])
				{
					printf("Error: Genome_MafAlign_To_CS_Main, error in initializing the cs file!\n");
					exit(EXIT_FAILURE);
				}
			}

			DestroyByteMatrix(pInitScore);
		}
		/* conservation window size */
		else if(strcmp(strLabel, "Conserve_Window") == 0)
		{
			sscanf(strLine, "%s %d", strLabel, &nConserveWinSize);
		}
		/* min p-value cutoff */
		else if(strcmp(strLabel, "Min_P_Cutoff") == 0)
		{
			sscanf(strLine, "%s %lf", strLabel, &dMinPcutoff);
		}
		/* alignment suffix */
		else if(strcmp(strLabel, "Alignment_Suf") == 0)
		{
			sscanf(strLine, "%s %s", strLabel, strAlnSuffix);
		}
		/* file by file */
		else if(strcmp(strLabel, "Alignment_File") == 0)
		{
			if(nk != nCompNum)
			{
				printf("Error: the comparison number not match!\n");
				return 0;
			}

			sscanf(strLine, "%s %s", strLabel, strChrName);

			nChrid = Genome_ChromosomeName_To_Index(strChrName, strRefSpecies);
			
			if( (nChrid<=0) || (nChrid>pChrSize->nHeight) )
			{
				printf("Error: %s is not supported as a reference genome!\n", strChrName);
				return 0;
			}
			nChrid = nChrid-1;
			nChrSize = IMGETAT(pChrSize, nChrid, 0);

			sprintf(strAlnFile, "%s%s%s", strAlnPath, strChrName, strAlnSuffix);

			/* calculate */
			Genome_MafAlign_To_CS_Chr(strAlnFile, strBGPath, 
				strChrName, vfpScore, 
				nSpeciesNum, vTag, 
				nRefid, nChrSize, 
				nTarid, pTargetChrSize,
				nCompNum, pCompPair,
				nWindowSize, nStepSize, nConserveWinSize,
				dMinPcutoff);
			nCount++;
		}
		/* do nothing */
		else
		{
		}
	}

	/* clear memeory */
	fclose(fpInfo);

	for(ni=0; ni<pTargetChrSize->nHeight; ni++)
	/* for(ni=20; ni<21; ni++) */
	{
		fclose(vfpScore[ni]);
		vfpScore[ni] = NULL;
	}
	free(vfpScore);

	for(ni=0; ni<nSpeciesNum; ni++)
	{
		if(vTag[ni] != NULL)
		{
			DeleteString(vTag[ni]);
			vTag[ni] = NULL;
		}
	}
	free(vTag);

	DestroyIntMatrix(pCompPair);
	DestroyIntMatrix(pChrSize);
	DestroyIntMatrix(pTargetChrSize);

	/* return */
	return nCount;
}

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
				double dMinPcutoff)
{
	/* define */

	/* window size */
	int nBinNum;

	/* background */
	struct DOUBLEMATRIX **vBGIdn;
	struct DOUBLEMATRIX **vBGAln;

	/* background */
	struct INTMATRIX *pPosInfo;
	struct BYTEMATRIX *pGapInfo;

	/* string */
	char strAlnLine[ALN_LENGTH];
	char strAlnSpecies[LINE_LENGTH];
	struct tagString **vAln;
	struct BYTEMATRIX **vCode;
	char strTag[LINE_LENGTH];
	char strLabel[LINE_LENGTH];
	char strBGName[LINE_LENGTH];
	char strTempChr[LINE_LENGTH];
	int nTempChr;
	char *pp1;

	/* count */
	int ni,nj,nk,nl,nx,nmargin;
	int numread;
	int numwritten;
	
	/* files */
	FILE *fpIn;
	int s1,s2;
	
	/* score */
	struct BYTEMATRIX *pCS;
	struct BYTEMATRIX *pCSF;
	struct DOUBLEMATRIX *pIdnCount;
	struct DOUBLEMATRIX *pTotCount;
	double dPhat,dP,dZ;
	double dMaxScore;
	double dScore;
	unsigned char nScore;
	double dMaxCov;

	/* loading variables */
	int nScanLen;
	int npos,nPos,nChrLen,nRc;
	char chRc;
	int nBinindex;

	/* strand 0: +; 1: - */
	int nRefStrand;
	int nTarStrand;
	int nRefOffset;
	int nTarOffset;

	/* init */
	if(nConserveWinSize%2 == 0)
	{
		printf("Error: Conserve_Window should be a odd number!\n");
		exit(EXIT_FAILURE);
	}
	nmargin = nConserveWinSize/2;

	if(dMinPcutoff <= 0.0)
	{
		printf("Error: min p-value cutoff should be greater than 0!\n");
		exit(EXIT_FAILURE);
	}
	dMaxScore = -log10(dMinPcutoff);

	nBinNum = (int)(nChrSize/nStepSize);
	if(nChrSize%nStepSize != 0)
		nBinNum++;
	
	/* create storage space */
	pPosInfo = NULL;
	pPosInfo = CreateIntMatrix(nSpeciesNum, 6); /* col1: aln exist?; col2: pos; col3: Pos; col4: rc; col5: length; col6: chr */
	if(pPosInfo == NULL)
	{
		printf("Error: cannot create memory for background calculation!\n");
		exit(EXIT_FAILURE);
	}

	pGapInfo = NULL;
	pGapInfo = CreateByteMatrix(nCompNum, 2);
	if(pGapInfo == NULL)
	{
		printf("Error: cannot create memory for background calculation!\n");
		exit(EXIT_FAILURE);
	}

	vAln = NULL;
	vAln = (struct tagString**)calloc(nSpeciesNum, sizeof(struct tagString *));
	if(vAln == NULL)
	{
		printf("Error: cannot create memory for background calculation!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nSpeciesNum; ni++)
	{
		vAln[ni] = CreateString(ALN_LENGTH);
		if(vAln[ni] == NULL)
		{
			printf("Error: cannot create memory for background calculation!\n");
			exit(EXIT_FAILURE);
		}
	}

	vCode = NULL;
	vCode = (struct BYTEMATRIX**)calloc(nCompNum, sizeof(struct BYTEMATRIX*));
	if(vCode == NULL)
	{
		printf("Error: cannot create memory for alignment calculation!\n");
		exit(EXIT_FAILURE);
	}

	

	vBGIdn = NULL;
	vBGIdn = (struct DOUBLEMATRIX **)calloc(nCompNum, sizeof(struct DOUBLEMATRIX *));
	if(vBGIdn == NULL)
	{
		printf("Error: cannot create memory for background calculation!\n");
		exit(EXIT_FAILURE);
	}

	vBGAln = NULL;
	vBGAln = (struct DOUBLEMATRIX **)calloc(nCompNum, sizeof(struct DOUBLEMATRIX *));
	if(vBGAln == NULL)
	{
		printf("Error: cannot create memory for background calculation!\n");
		exit(EXIT_FAILURE);
	}


	for(ni=0; ni<nCompNum; ni++)
	{
		vBGIdn[ni] = CreateDoubleMatrix(nBinNum, 1);
		if(vBGIdn[ni] == NULL)
		{
			printf("Error: cannot create memory for background calculation!\n");
			exit(EXIT_FAILURE);
		}
		vBGAln[ni] = CreateDoubleMatrix(nBinNum, 1);
		if(vBGAln[ni] == NULL)
		{
			printf("Error: cannot create memory for background calculation!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* init */

	/* load background */
	for(ni=0; ni<nCompNum; ni++)
	{
		s1 = IMGETAT(pCompPair, ni, 0);
		s2 = IMGETAT(pCompPair, ni, 1);
		
		sprintf(strBGName, "%s%s_%s-%s.bg", strBGPath, strChrName, vTag[s1]->m_pString, vTag[s2]->m_pString);
	
		fpIn = NULL;
		fpIn = fopen(strBGName, "rb");
		if(fpIn == NULL)
		{
			printf("Error: Genome_MafAlign_To_CS_Chr, cannot open background file!\n");
			exit(EXIT_FAILURE);
		}
		numread = fread(vBGIdn[ni]->pMatElement, sizeof(double), nBinNum, fpIn);

		if(numread != nBinNum)
		{
			printf("Error: load background wrong!\n");
			exit(EXIT_FAILURE);
		}
		fclose(fpIn);
		
		sprintf(strBGName, "%s%s_%s-%s.cov", strBGPath, strChrName, vTag[s1]->m_pString, vTag[s2]->m_pString);

		fpIn = NULL;
		fpIn = fopen(strBGName, "rb");
		if(fpIn == NULL)
		{
			printf("Error: Genome_MafAlign_To_CS_Chr, cannot open background file!\n");
			exit(EXIT_FAILURE);
		}
		numread = fread(vBGAln[ni]->pMatElement, sizeof(double), nBinNum, fpIn);

		if(numread != nBinNum)
		{
			printf("Error: load background wrong!\n");
			exit(EXIT_FAILURE);
		}
		fclose(fpIn);
	}


	/* load in alignment and calculate conservation score */
	fpIn = NULL;
	fpIn = fopen(strAlnFile, "rt");
	if(fpIn == NULL)
	{
		printf("Error: cannot open alignment file!\n");
		return PROC_FAILURE;
	}

	while(fgets(strAlnLine, ALN_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strAlnLine);
		StrTrimRight(strAlnLine);
				
		/* ignore annotations */
		if(strAlnLine[0] == '#')
		{
			continue;
		}

		/* init */
		else if(strAlnLine[0] == 'a')
		{
			for(ni=0; ni<nSpeciesNum; ni++)
			{
				for(nj=0; nj<6; nj++)
				{
					IMSETAT(pPosInfo, ni, nj, 0);
				}
			}
			for(ni=0; ni<nCompNum; ni++)
			{
				for(nj=0; nj<2; nj++)
				{
					BMSETAT(pGapInfo, ni, nj, 0);
				}
			}
		}

		/* load alignment */
		else if(strAlnLine[0] == 's')
		{
			sscanf(strAlnLine, "%s %s", strLabel, strTag);
			for(ni=0; ni<nSpeciesNum; ni++)
			{
				if(strstr(strTag, vTag[ni]->m_pString) == strTag)
				{
					break;
				}
			}
			if(ni >= nSpeciesNum)
			{
				printf("Error: species not found!\n");
				exit(EXIT_FAILURE);
			}

			sscanf(strAlnLine, "%s %s %d %d %c %d %s",
				strLabel, strTag, &npos, &nPos, &chRc,
					&nChrLen, vAln[ni]->m_pString);
			
			StrMakeUpper(vAln[ni]->m_pString);

			pp1 = strchr(strTag, '.');
			pp1++;
			strcpy(strTempChr, pp1);
			Genome_SpeciesAbbr_To_Species(vTag[ni]->m_pString, strAlnSpecies);
			nTempChr = Genome_ChromosomeName_To_Index(strTempChr, strAlnSpecies)-1;
			nPos += npos;
			if(chRc == '-') 
				nRc = 1;
			else
				nRc = 0;


			IMSETAT(pPosInfo, ni, 0, 1);
			IMSETAT(pPosInfo, ni, 1, npos);
			IMSETAT(pPosInfo, ni, 2, nPos);
			IMSETAT(pPosInfo, ni, 3, nRc);
			IMSETAT(pPosInfo, ni, 4, nChrLen);
			IMSETAT(pPosInfo, ni, 5, nTempChr);
		}

		/* process */
		else if(strAlnLine[0] == '\0')
		{
			if(IMGETAT(pPosInfo, nRefid, 0) == 0)
				continue;

			if(IMGETAT(pPosInfo, nTarid, 0) == 0)
				continue;

			nScanLen = strlen(vAln[nRefid]->m_pString);

			/* create memory */
			for(ni=0; ni<nCompNum; ni++)
			{
				vCode[ni] = CreateByteMatrix(1, nScanLen);
				if(vCode[ni] == NULL)
				{
					printf("Error: cannot create memory for alignment calculation!\n");
					exit(EXIT_FAILURE);
				}
			}
			pCS = NULL;
			pCS = CreateByteMatrix(1, nScanLen);
			if(pCS == NULL)
			{
				printf("Error: cannot create memory for conservation calculation!\n");
				exit(EXIT_FAILURE);
			}
			
			/* encode alignment: 0-general, 1-match, 2-mismatch, 3-N
					10-gap open in seq1, 11-gap extension in seq1
					20-gap open in seq2, 21-gap extension in seq2 */

			/* code alignment */
			nl = IMGETAT(pPosInfo, nRefid, 1);
			for(nj=0; nj<nScanLen; nj++)
			{
				/* pair by pair */
				for(ni=0; ni<nCompNum; ni++)
				{
					s1 = IMGETAT(pCompPair, ni, 0);
					s2 = IMGETAT(pCompPair, ni, 1);

					/* if we have alignment */
					if((IMGETAT(pPosInfo, s1, 0) == 1) && (IMGETAT(pPosInfo, s2, 0) == 1))
					{
						if(vAln[s2]->m_pString[nj] != '-')
							BMSETAT(pGapInfo, ni, 1, 0);

						/* s1 gap */
						if(vAln[s1]->m_pString[nj] == '-')
						{
							/* s1 gap vs s2 base */
							if(vAln[s2]->m_pString[nj] != '-')
							{
								/* s1 gap open */
								if(BMGETAT(pGapInfo, ni, 0) == 0)
								{
									vCode[ni]->pMatElement[nj] = 10;
									BMSETAT(pGapInfo, ni, 0, 1);
								}
								/* s1 gap continue */
								else
								{
									vCode[ni]->pMatElement[nj] = 11;
								}
							}

							/* s1 gap vs s2 gap */
							else
							{
								vCode[ni]->pMatElement[nj] = 0;
							}
						}
						/* s1 base */
						else
						{
							BMSETAT(pGapInfo, ni, 0, 0);

							/* s1 base vs s2 gap */
							if(vAln[s2]->m_pString[nj] == '-')
							{
								/* s2 gap open */
								if(BMGETAT(pGapInfo, ni, 1) == 0)
								{
									vCode[ni]->pMatElement[nj] = 20;
									BMSETAT(pGapInfo, ni, 1, 1);
								}
								/* s2 gap continue */
								else
								{
									vCode[ni]->pMatElement[nj] = 21;
								}
							}

							/* s1 base vs s2 base */
							else
							{
								/* effective base */
								if((vAln[s1]->m_pString[nj] != 'N') && (vAln[s2]->m_pString[nj] != 'N'))
								{
									/* substitution */
									if(vAln[s1]->m_pString[nj] != vAln[s2]->m_pString[nj])
									{
										vCode[ni]->pMatElement[nj] = 2;
									}
									/* identical */
									else
									{
										vCode[ni]->pMatElement[nj] = 1;
									}
								}
								else
								{
									vCode[ni]->pMatElement[nj] = 3;
								}
							}
						}
					}

					/* if we don't have alignment */
					else
					{
						vCode[ni]->pMatElement[nj] = 0;
					}
				}

				if(vAln[nRefid]->m_pString[nj] != '-')
				{
					nl++;
				}	
			}

			if(nl != IMGETAT(pPosInfo, nRefid, 2))
			{
				printf("Error: coordinate wrong!\n");
				exit(EXIT_FAILURE);
			}


			/* calculate conservation score */
			pIdnCount = NULL;
			pIdnCount = CreateDoubleMatrix(nCompNum, 2);
			if(pIdnCount == NULL)
			{
				printf("Error: cannot create counting matrix!\n");
				exit(EXIT_FAILURE);
			}
			pTotCount = NULL;
			pTotCount = CreateDoubleMatrix(nCompNum, 2);
			if(pTotCount == NULL)
			{
				printf("Error: cannot create counting matrix!\n");
				exit(EXIT_FAILURE);
			}

			nl = IMGETAT(pPosInfo, nRefid, 1);
			nRefStrand = IMGETAT(pPosInfo, nRefid, 3);
			nTarStrand = IMGETAT(pPosInfo, nTarid, 3);
			nRefOffset = IMGETAT(pPosInfo, nRefid, 4)-1;
			nTarOffset = IMGETAT(pPosInfo, nTarid, 4)-1;
			for(nj=0; nj<nConserveWinSize; nj++)
			{
				if(nj == nScanLen)
					break;

				/* check every pair */
				for(ni=0; ni<nCompNum; ni++)
				{
					switch(vCode[ni]->pMatElement[nj])
					{
						case 1: pIdnCount->pMatElement[ni] = pIdnCount->pMatElement[ni]+1;
							pTotCount->pMatElement[ni] = pTotCount->pMatElement[ni]+1;
							break;
						case 2: pTotCount->pMatElement[ni] = pTotCount->pMatElement[ni]+1;
							break;
						case 10: pTotCount->pMatElement[ni] = pTotCount->pMatElement[ni]+1;
							break;
						case 11: pTotCount->pMatElement[ni] = pTotCount->pMatElement[ni]+1;
							break;
						case 20: pTotCount->pMatElement[ni] = pTotCount->pMatElement[ni]+1;
							break;
						case 21: pTotCount->pMatElement[ni] = pTotCount->pMatElement[ni]+1;
							break;
					}
				}
			}

			if(nRefStrand == 0) 
				nBinindex = nl/nStepSize;
			else
				nBinindex = (nRefOffset-nl)/nStepSize;

			
			dMaxCov = 0.0;
			for(ni=0; ni<nCompNum; ni++)
			{
				if(vBGAln[ni]->pMatElement[nBinindex] > dMaxCov)
					dMaxCov = vBGAln[ni]->pMatElement[nBinindex];
			}
			
			dScore = 0.0;
			for(ni=0; ni<nCompNum; ni++)
			{
				if(pTotCount->pMatElement[ni] >= 10)
				{
					dPhat = (double)(pIdnCount->pMatElement[ni])/(double)pTotCount->pMatElement[ni];
					dP = vBGIdn[ni]->pMatElement[nBinindex] + 1e-6;
					dZ = (dPhat-dP)/sqrt(dP*(1-dP)/(double)pTotCount->pMatElement[ni]);
					if(vBGAln[ni]->pMatElement[nBinindex] > 0.0)
						dPhat = vBGAln[ni]->pMatElement[nBinindex]/(dMaxCov*(1.0+exp(dZ)))+dMinPcutoff;
					else
						dPhat = 1.0;
					if(dPhat > 1.0)
						dPhat = 1.0;
					dScore += -log10(dPhat);
				}
			}
			dScore /= (double)nCompNum;
			nScore = (unsigned char)(dScore*255/dMaxScore);

			for(nk=0; nk<=nmargin; nk++)
			{
				if(nk == nScanLen)
					break;
				pCS->pMatElement[nk] = nScore;
			}

			/* processing continues */
			for(;nj<nScanLen; nj++)
			{
				if(vAln[nRefid]->m_pString[nj-nConserveWinSize] != '-')
				{
					nl++;
				}	
				
				/* check every pair */
				for(ni=0; ni<nCompNum; ni++)
				{
					/* subtract previous */
					switch(vCode[ni]->pMatElement[nj-nConserveWinSize])
					{
						case 1: pIdnCount->pMatElement[ni] = pIdnCount->pMatElement[ni]-1;
							pTotCount->pMatElement[ni] = pTotCount->pMatElement[ni]-1;
							break;
						case 2: pTotCount->pMatElement[ni] = pTotCount->pMatElement[ni]-1;
							break;
						case 10: pTotCount->pMatElement[ni] = pTotCount->pMatElement[ni]-1;
							break;
						case 11: pTotCount->pMatElement[ni] = pTotCount->pMatElement[ni]-1;
							break;
						case 20: pTotCount->pMatElement[ni] = pTotCount->pMatElement[ni]-1;
							break;
						case 21: pTotCount->pMatElement[ni] = pTotCount->pMatElement[ni]-1;
							break;
					}
					/* add next */
					switch(vCode[ni]->pMatElement[nj])
					{
						case 1: pIdnCount->pMatElement[ni] = pIdnCount->pMatElement[ni]+1;
							pTotCount->pMatElement[ni] = pTotCount->pMatElement[ni]+1;
							break;
						case 2: pTotCount->pMatElement[ni] = pTotCount->pMatElement[ni]+1;
							break;
						case 10: pTotCount->pMatElement[ni] = pTotCount->pMatElement[ni]+1;
							break;
						case 11: pTotCount->pMatElement[ni] = pTotCount->pMatElement[ni]+1;
							break;
						case 20: pTotCount->pMatElement[ni] = pTotCount->pMatElement[ni]+1;
							break;
						case 21: pTotCount->pMatElement[ni] = pTotCount->pMatElement[ni]+1;
							break;
					}
				}

				/* get score */
				if(nRefStrand == 0) 
					nBinindex = nl/nStepSize;
				else
					nBinindex = (nRefOffset-nl)/nStepSize;

				dMaxCov = 0.0;
				for(ni=0; ni<nCompNum; ni++)
				{
					if(vBGAln[ni]->pMatElement[nBinindex] > dMaxCov)
						dMaxCov = vBGAln[ni]->pMatElement[nBinindex];
				}

				dScore = 0.0;
				for(ni=0; ni<nCompNum; ni++)
				{
					if(pTotCount->pMatElement[ni] >= 10)
					{
						dPhat = (double)(pIdnCount->pMatElement[ni])/(double)pTotCount->pMatElement[ni];
						dP = vBGIdn[ni]->pMatElement[nBinindex] + 1e-6;
						dZ = (dPhat-dP)/sqrt(dP*(1-dP)/(double)pTotCount->pMatElement[ni]);
						if(vBGAln[ni]->pMatElement[nBinindex] > 0.0)
							dPhat = vBGAln[ni]->pMatElement[nBinindex]/(dMaxCov*(1.0+exp(dZ)))+dMinPcutoff;
						else
							dPhat = 1.0;

						if(dPhat > 1.0)
							dPhat = 1.0;
						dScore += -log10(dPhat);
					}
				}
				dScore /= (double)nCompNum;
				nScore = (unsigned char)(dScore*255/dMaxScore);

				pCS->pMatElement[nk] = nScore;
				nk++;
			}

			for(; nk<nScanLen; nk++)
			{
				pCS->pMatElement[nk] = nScore;
			}

			/* write to files */
			pCSF = NULL;
			pCSF = CreateByteMatrix(1, nScanLen);
			if(pCSF == NULL)
			{
				printf("Error: cannot create buffer for writing conservation score!\n");
			}
			/* + strand */
			if(nTarStrand == 0)
			{
				nl = IMGETAT(pPosInfo, nTarid, 1);
				nTempChr = IMGETAT(pPosInfo, nTarid, 5);

				if(fseek(vfpScore[nTempChr], nl, SEEK_SET) != 0)
				{
					printf("Error: fseek error!\n");
					exit(EXIT_FAILURE);
				}

				nx = 0;
				for(nj=0; nj<nScanLen; nj++)
				{
					if(vAln[nTarid]->m_pString[nj] != '-')
					{
						pCSF->pMatElement[nx] = pCS->pMatElement[nj];
						nx++;
						nl++;
					}	
				}
				if(nl != IMGETAT(pPosInfo, nTarid, 2))
				{
					printf("Error: coordinate wrong!\n");
					exit(EXIT_FAILURE);
				}

				

				numwritten = fwrite(pCSF->pMatElement, sizeof(unsigned char), nx, vfpScore[nTempChr]);
				if(numwritten != nx)
				{
					printf("Error: coordinate wrong!\n");
					exit(EXIT_FAILURE);
				}
			}
			/* - strand */
			else
			{
				nl = nTarOffset - IMGETAT(pPosInfo, nTarid, 2) + 1;
				nTempChr = IMGETAT(pPosInfo, nTarid, 5);

				if(fseek(vfpScore[nTempChr], nl, SEEK_SET) != 0)
				{
					printf("Error: fseek error!\n");
					exit(EXIT_FAILURE);
				}
				
				nx = 0;
				for(nj=(nScanLen-1); nj>=0; nj--)
				{
					if(vAln[nTarid]->m_pString[nj] != '-')
					{
						pCSF->pMatElement[nx] = pCS->pMatElement[nj];
						nx++;
						nl++;
					}	
				}
				if(nl != (nTarOffset - IMGETAT(pPosInfo, nTarid, 1) + 1) )
				{
					printf("Error: coordinate wrong!\n");
					exit(EXIT_FAILURE);
				}

				numwritten = fwrite(pCSF->pMatElement, sizeof(unsigned char), nx, vfpScore[nTempChr]);
				if(numwritten != nx)
				{
					printf("Error: coordinate wrong!\n");
					exit(EXIT_FAILURE);
				}
			}
			
			DestroyByteMatrix(pCSF);
	
			/* release memory */
			DestroyDoubleMatrix(pIdnCount);
			DestroyDoubleMatrix(pTotCount);
			for(ni=0; ni<nCompNum; ni++)
			{
				DestroyByteMatrix(vCode[ni]);
				vCode[ni] = NULL;
			}
			DestroyByteMatrix(pCS);
		}

		/* error */
		else
		{
			printf("fatal error: loading error!\n");
			exit(EXIT_FAILURE);
			/* continue; */
		}
	}

	fclose(fpIn);


	/* clear memory */
	DestroyIntMatrix(pPosInfo);
	DestroyByteMatrix(pGapInfo);

	for(ni=0; ni<nSpeciesNum; ni++)
	{
		DeleteString(vAln[ni]);
		vAln[ni] = NULL;
	}
	free(vAln);

	free(vCode);

	for(ni=0; ni<nCompNum; ni++)
	{
		DestroyDoubleMatrix(vBGIdn[ni]);
		vBGIdn[ni] = NULL;
		DestroyDoubleMatrix(vBGAln[ni]);
		vBGAln[ni] = NULL;
	}
	free(vBGIdn);
	free(vBGAln);
	
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Genome_MafAlign_To_FootPrint_Main:                                     */
/*  Convert maf alignments to phylogenetic footprints.                     */
/*  Footprints are specified by conditions such as human+mouse+dog>=2 and  */
/*  window size = 3, etc.                                                  */
/*  Allowed conditions:                                                    */
/*  0: ==; 1: >=; 2: <=; 3: !=; 4: >; 5: <.                                */
/*  return number of alignment files processed.                            */
/* ----------------------------------------------------------------------- */ 
int Genome_MafAlign_To_FootPrint_Main(char strParamPath[])
{
	/* define */
	int nCount;
	int ni,nk,nx;
	int nrlen;
	int numwritten;
	char *chp,*chp2,*chpn;
	int nCondLen;
	int nCompType = 0;
	
	/* files */
	FILE *fpInfo;
	char strLine[MED_LINE_LENGTH];
	char strAlnPath[LINE_LENGTH];
	char strAlnFile[LINE_LENGTH];
	char strOutPath[LINE_LENGTH];
	char strInFile[LINE_LENGTH];
	char strRefSpecies[LINE_LENGTH];
	char strTarSpecies[LINE_LENGTH];
	char strOutFile[LINE_LENGTH];
	char strAlnSuffix[LINE_LENGTH];
	char strChrName[LINE_LENGTH];
	char strTemp[LINE_LENGTH];

	FILE **vfpScore;
	struct BYTEMATRIX *pInitScore;


	char strLabel[LINE_LENGTH];

	/* number of species */
	int nSpeciesNum;
	/* species id */
	struct tagString **vTag;

	/* reference species */
	int nRefid;
	int nChrid;
	int nChrSize;
	int nTarid;

	/* chromosome size */
	struct INTMATRIX *pChrSize;
	struct INTMATRIX *pTargetChrSize;

	/* comparisons */
	int nConditionNum = 0;
	struct INTMATRIX **vCondition;

	/* bg param */
	int nStepSize = 100000;
	int nWindowSize = 3;
	
	/* init */
	fpInfo = NULL;
	fpInfo = fopen(strParamPath, "rt");
	if(fpInfo == NULL)
	{
		printf("Error: Genome_MafAlign_To_FootPrint_Main, cannot open parameter info file!\n");
		return 0;
	}

	nCount = 0;

	/* load parameters */
	while(fgets(strLine, MED_LINE_LENGTH, fpInfo) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		
		sscanf(strLine, "%s ", strLabel);

		/* aln path */
		if(strcmp(strLabel, "Alignment_Dir") == 0)
		{
			chp = strLine+13;
			StrTrimLeft(chp);
			strcpy(strAlnPath, chp);

			AdjustDirectoryPath(strAlnPath);
		}
		/* species number */
		else if(strcmp(strLabel, "Species_Num") == 0)
		{
			sscanf(strLine, "%s %d", strLabel, &nSpeciesNum);
			if(nSpeciesNum <= 1)
			{
				printf("Error: the species number should be >= 2!\n");
				return 0;
			}

			vTag = (struct tagString **)calloc(nSpeciesNum, sizeof(struct tagString *));
		}
		/* species tag */
		else if(strcmp(strLabel, "ID") == 0)
		{
			sscanf(strLine, "%s %d %s", strLabel, &ni, strTemp);
			if( (ni >= nSpeciesNum) || (ni < 0) )
			{
				printf("Error: the species id should be < species number and >=0!\n");
				return 0;
			}

			StringAddTail((vTag+ni), strTemp);
		}
		/* reference */
		else if(strcmp(strLabel, "Reference") == 0)
		{
			sscanf(strLine, "%s %d", strLabel, &nRefid);
			if((nRefid < 0) || (nRefid >=nSpeciesNum))
			{
				printf("Error: reference id out of range");
				return 0;
			}
		}
		/* reference species */
		else if(strcmp(strLabel, "Ref_Species") == 0)
		{
			sscanf(strLine, "%s %s", strLabel, strRefSpecies);
		}
		/* reference chrlen */
		else if(strcmp(strLabel, "Ref_Chr_Size") == 0)
		{
			chp = strLine+12;
			StrTrimLeft(chp);
			strcpy(strInFile, chp);

			pChrSize = NULL;
			pChrSize = IMLOAD(strInFile);
			if(pChrSize == NULL)
			{
				printf("Error: cannot load chromosome length!\n");
				return 0;
			}
		}
		/* target */
		else if(strcmp(strLabel, "Target") == 0)
		{
			sscanf(strLine, "%s %d", strLabel, &nTarid);
			if((nTarid < 0) || (nTarid >=nSpeciesNum))
			{
				printf("Error: target id out of range");
				return 0;
			}
		}
		/* target species */
		else if(strcmp(strLabel, "Tar_Species") == 0)
		{
			sscanf(strLine, "%s %s", strLabel, strTarSpecies);
		}
		/* target chrlen */
		else if(strcmp(strLabel, "Tar_Chr_Size") == 0)
		{
			chp = strLine+12;
			StrTrimLeft(chp);
			strcpy(strInFile, chp);

			pTargetChrSize = NULL;
			pTargetChrSize = IMLOAD(strInFile);
			if(pTargetChrSize == NULL)
			{
				printf("Error: cannot load chromosome length!\n");
				return 0;
			}
		}
		/* comparison number */
		else if(strcmp(strLabel, "Comparisons_Num") == 0)
		{
			sscanf(strLine, "%s %d", strLabel, &nConditionNum);
			if(nConditionNum < 1)
			{
				printf("Error: the comparison number should be >= 1!\n");
				return 0;
			}
			vCondition = NULL;
			vCondition = (struct INTMATRIX **)calloc(nConditionNum, sizeof(struct INTMATRIX*));
			if(vCondition == NULL)
			{
				printf("Error: cannot create comparison storage space!\n");
				return 0;
			}
			nk = 0;
		}
		/* comparison id */
		else if(strcmp(strLabel, "COMP") == 0)
		{
			if(nk >= nConditionNum)
			{
				printf("Error: the comparison number out of range!\n");
				return 0;
			}

			chp = strLine+4;
			StrTrimLeft(chp);
			nCondLen = 1;

			chp2 = strchr(chp, '+');
			while(chp2 != NULL)
			{
				nCondLen++;
				chp = chp2+1;
				chp2 = strchr(chp, '+');
			}

			chp2 = strpbrk(chp, "><=!");
			if(chp2 == NULL)
			{
				printf("Error: the condition is not specified correctly!\n");
				exit(EXIT_FAILURE);
			}

			nCondLen += 2;

			vCondition[nk] = NULL;
			vCondition[nk] = CreateIntMatrix(1, nCondLen);
			if(vCondition[nk] == NULL)
			{
				printf("Error: cannot load condition correctly, %s!\n", strLine);
				exit(EXIT_FAILURE);
			}

			chp = strLine+4;
			StrTrimLeft(chp);
			nx = 0;

			chp2 = strchr(chp, '+');
			while(chp2 != NULL)
			{
				*chp2 = '\0';
				StrTrimLeft(chp);
				StrTrimRight(chp);
				vCondition[nk]->pMatElement[nx] = atoi(chp);
				nx++;
				
				chp = chp2+1;
				chp2 = strchr(chp, '+');
			}

			chp2 = strpbrk(chp, "><=!");
			if(strstr(chp2, "==") == chp2)
			{
				nCompType = 0;
				chpn = chp2+2;
			}
			else if(strstr(chp2, ">=") == chp2)
			{
				nCompType = 1;
				chpn = chp2+2;
			}
			else if(strstr(chp2, "<=") == chp2)
			{
				nCompType = 2;
				chpn = chp2+2;
			}
			else if(strstr(chp2, "!=") == chp2)
			{
				nCompType = 3;
				chpn = chp2+2;
			}
			else if(strstr(chp2, ">") == chp2)
			{
				nCompType = 4;
				chpn = chp2+1;
			}
			else if(strstr(chp2, "<") == chp2)
			{
				nCompType = 5;
				chpn = chp2+1;
			}

			*chp2 = '\0';
			StrTrimLeft(chp);
			StrTrimRight(chp);
			vCondition[nk]->pMatElement[nx] = atoi(chp);
			nx++;

			vCondition[nk]->pMatElement[nx] = nCompType;
			nx++;
			StrTrimLeft(chpn);
			if(*chpn == '\0')
			{
				printf("Error: the condition is not specified correctly!\n");
				exit(EXIT_FAILURE);
			}
			vCondition[nk]->pMatElement[nx] = atoi(chpn);

			nk++;
		}
		/* output path */
		else if(strcmp(strLabel, "Output_Path") == 0)
		{
			chp = strLine+11;
			StrTrimLeft(chp);
			strcpy(strOutPath, chp);

			AdjustDirectoryPath(strOutPath);
			vfpScore = NULL;
			vfpScore = (FILE **)calloc(pTargetChrSize->nHeight, sizeof(FILE *));
			pInitScore = NULL;
			pInitScore = CreateByteMatrix(1, nStepSize);
			for(ni=0; ni<pTargetChrSize->nHeight; ni++)
			{
				Genome_Index_To_ChromosomeName(strChrName, strTarSpecies, (ni+1));
				
				if(strcmp(strChrName, "unplaced") == 0)
				{
					printf("Error: %s is not supported as a target genome!\n");
					exit(EXIT_FAILURE);
				}
				
				sprintf(strOutFile, "%s%s.cs", strOutPath, strChrName);
				
				vfpScore[ni] = fopen(strOutFile, "w+b");
				if(vfpScore[ni] == NULL)
				{
					printf("Error: Genome_MafAlign_To_FootPrint_Main, cannot open output file!\n");
					exit(EXIT_FAILURE);
				}
				nrlen = pTargetChrSize->pMatElement[ni];
				numwritten = 0;
				while(nrlen >= nStepSize)
				{
					numwritten += fwrite(pInitScore->pMatElement, sizeof(unsigned char), nStepSize, vfpScore[ni]);
					nrlen -= nStepSize;
				}
				if(nrlen > 0)
				{
					numwritten += fwrite(pInitScore->pMatElement, sizeof(unsigned char), nrlen, vfpScore[ni]);
				}

				if(numwritten != pTargetChrSize->pMatElement[ni])
				{
					printf("Error: Genome_MafAlign_To_FootPrint_Main, error in initializing the cs file!\n");
					exit(EXIT_FAILURE);
				}
			}

			DestroyByteMatrix(pInitScore);
		}
		/* conservation window size */
		else if(strcmp(strLabel, "Window_Size") == 0)
		{
			sscanf(strLine, "%s %d", strLabel, &nWindowSize);
		}
		/* alignment suffix */
		else if(strcmp(strLabel, "Alignment_Suf") == 0)
		{
			chp = strLine+13;
			StrTrimLeft(chp);
			strcpy(strAlnSuffix, chp);
		}
		/* file by file */
		else if(strcmp(strLabel, "Alignment_File") == 0)
		{
			if(nk != nConditionNum)
			{
				printf("Error: the comparison number not match!\n");
				return 0;
			}

			sscanf(strLine, "%s %s", strLabel, strChrName);

			nChrid = Genome_ChromosomeName_To_Index(strChrName, strRefSpecies);
			
			if( (nChrid<=0) || (nChrid>pChrSize->nHeight) )
			{
				printf("Error: %s is not supported as a reference genome!\n", strChrName);
				return 0;
			}
			nChrid = nChrid-1;
			nChrSize = IMGETAT(pChrSize, nChrid, 0);

			sprintf(strAlnFile, "%s%s%s", strAlnPath, strChrName, strAlnSuffix);

			/* calculate */
			Genome_MafAlign_To_FootPrint_Chr(strAlnFile, 
				strChrName, vfpScore, 
				nSpeciesNum, vTag, 
				nRefid, nChrSize, 
				nTarid, pTargetChrSize,
				nConditionNum, vCondition,
				nWindowSize);
			nCount++;
		}
		/* do nothing */
		else
		{
		}
	}

	/* clear memeory */
	fclose(fpInfo);

	for(ni=0; ni<pTargetChrSize->nHeight; ni++)
	{
		fclose(vfpScore[ni]);
		vfpScore[ni] = NULL;
	}
	free(vfpScore);

	for(ni=0; ni<nSpeciesNum; ni++)
	{
		if(vTag[ni] != NULL)
		{
			DeleteString(vTag[ni]);
			vTag[ni] = NULL;
		}
	}
	free(vTag);

	for(ni=0; ni<nConditionNum; ni++)
	{
		DestroyIntMatrix(vCondition[ni]);
	}
	free(vCondition);

	DestroyIntMatrix(pChrSize);
	DestroyIntMatrix(pTargetChrSize);

	/* return */
	return nCount;
}

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
				int nWindowSize)
{
	/* define */
	struct INTMATRIX *pPosInfo;
	struct tagString **vAln;
	
	/* string */
	char strAlnLine[ALN_LENGTH];
	char strTag[LINE_LENGTH];
	char strLabel[LINE_LENGTH];
	int npos,nPos,nChrLen,nRc;
	char chRc;
	char strTempChr[LINE_LENGTH];
	int nTempChr;
	char strAlnSpecies[LINE_LENGTH];
	int nConditionOK,nOK;
	int nMatchNum;
	int nCompNum;
	int nCompType;
	int nCompCut;
	int s1,s2;
	char *vSeq1,*vSeq2;

	char *pp1;

	/* count */
	int ni,nj,nk,nl,nx;
	int numwritten;
	
	/* files */
	FILE *fpIn;
	
	/* score */
	struct BYTEMATRIX *pCS;
	struct BYTEMATRIX *pCSF;
	
	/* loading variables */
	int nScanLen;
	int nEffecLen;
	
	/* strand 0: +; 1: - */
	int nTarStrand;
	int nTarOffset;

	/* init */
	/* create storage space */
	pPosInfo = NULL;
	pPosInfo = CreateIntMatrix(nSpeciesNum, 6); /* col1: aln exist?; col2: pos; col3: Pos; col4: rc; col5: length; col6: chr */
	if(pPosInfo == NULL)
	{
		printf("Error: cannot create memory for background calculation!\n");
		exit(EXIT_FAILURE);
	}

	vAln = NULL;
	vAln = (struct tagString**)calloc(nSpeciesNum, sizeof(struct tagString *));
	if(vAln == NULL)
	{
		printf("Error: cannot create memory for background calculation!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nSpeciesNum; ni++)
	{
		vAln[ni] = CreateString(ALN_LENGTH);
		if(vAln[ni] == NULL)
		{
			printf("Error: cannot create memory for background calculation!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* load in alignment and calculate conservation score */
	fpIn = NULL;
	fpIn = fopen(strAlnFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: cannot open alignment file!\n");
		return PROC_FAILURE;
	}

	while(fgets(strAlnLine, ALN_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strAlnLine);
		StrTrimRight(strAlnLine);
				
		/* ignore annotations */
		if(strAlnLine[0] == '#')
		{
			continue;
		}

		/* init */
		else if(strAlnLine[0] == 'a')
		{
			for(ni=0; ni<nSpeciesNum; ni++)
			{
				for(nj=0; nj<6; nj++)
				{
					IMSETAT(pPosInfo, ni, nj, 0);
				}
			}
		}

		/* load alignment */
		else if(strAlnLine[0] == 's')
		{
			sscanf(strAlnLine, "%s %s", strLabel, strTag);
			for(ni=0; ni<nSpeciesNum; ni++)
			{
				if(strstr(strTag, vTag[ni]->m_pString) == strTag)
				{
					break;
				}
			}
			if(ni >= nSpeciesNum)
			{
				printf("Error: species not found!\n");
				exit(EXIT_FAILURE);
			}

			sscanf(strAlnLine, "%s %s %d %d %c %d %s",
				strLabel, strTag, &npos, &nPos, &chRc,
					&nChrLen, vAln[ni]->m_pString);
			
			StrMakeUpper(vAln[ni]->m_pString);

			pp1 = strchr(strTag, '.');
			pp1++;
			strcpy(strTempChr, pp1);
			Genome_SpeciesAbbr_To_Species(vTag[ni]->m_pString, strAlnSpecies);
			nTempChr = Genome_ChromosomeName_To_Index(strTempChr, strAlnSpecies)-1;
			nPos += npos;
			if(chRc == '-') 
				nRc = 1;
			else
				nRc = 0;


			IMSETAT(pPosInfo, ni, 0, 1);
			IMSETAT(pPosInfo, ni, 1, npos);
			IMSETAT(pPosInfo, ni, 2, nPos);
			IMSETAT(pPosInfo, ni, 3, nRc);
			IMSETAT(pPosInfo, ni, 4, nChrLen);
			IMSETAT(pPosInfo, ni, 5, nTempChr);
		}

		/* process */
		else if(strAlnLine[0] == '\0')
		{
			if(IMGETAT(pPosInfo, nRefid, 0) == 0)
				continue;

			if(IMGETAT(pPosInfo, nTarid, 0) == 0)
				continue;

			nScanLen = strlen(vAln[nRefid]->m_pString);

			/* create memory */
			pCS = NULL;
			pCS = CreateByteMatrix(1, nScanLen);
			if(pCS == NULL)
			{
				printf("Error: cannot create memory for conservation calculation!\n");
				exit(EXIT_FAILURE);
			}
			
			/* encode alignment: 0-general, 1-match, 2-mismatch, 3-N
					10-gap open in seq1, 11-gap extension in seq1
					20-gap open in seq2, 21-gap extension in seq2 */

			/* code alignment */
			s1 = nTarid;
			if(IMGETAT(pPosInfo, s1, 0) == 1)
			{
				nEffecLen = nScanLen-nWindowSize+1;
				for(ni=0; ni<nEffecLen; ni++)
				{
					nConditionOK = 1;
					for(nj=0; nj<nConditionNum; nj++)
					{
						nMatchNum = 0;
						nCompNum = vCondition[nj]->nWidth-2;
						nCompType = vCondition[nj]->pMatElement[nCompNum];
						nCompCut = vCondition[nj]->pMatElement[nCompNum+1];
						for(nk=0; nk<nCompNum; nk++)
						{
							s2 = vCondition[nj]->pMatElement[nk];
							if(IMGETAT(pPosInfo, s2, 0) == 1)
							{
								vSeq1 = vAln[s1]->m_pString+ni;
								vSeq2 = vAln[s2]->m_pString+ni;

								nOK = 1;
								for(nl=0; nl<nWindowSize; nl++)
								{
									if(vSeq1[nl] != vSeq2[nl])
									{
										nOK = 0;
										break;
									}
									else if( (vSeq1[nl] == '-') || (vSeq1[nl] == 'N') || (vSeq1[nl] == 'n'))
									{
										nOK = 0;
										break;
									}
								}

								if(nOK == 1)
									nMatchNum++;
							}
						}

						if(nCompType == 0)
						{
							if(nMatchNum != nCompCut)
								nConditionOK = 0;
						}
						else if(nCompType == 1)
						{
							if(nMatchNum < nCompCut)
								nConditionOK = 0;
						}
						else if(nCompType == 2)
						{
							if(nMatchNum > nCompCut)
								nConditionOK = 0;
						}
						else if(nCompType == 3)
						{
							if(nMatchNum == nCompCut)
								nConditionOK = 0;
						}
						else if(nCompType == 4)
						{
							if(nMatchNum <= nCompCut)
								nConditionOK = 0;
						}
						else if(nCompType == 5)
						{
							if(nMatchNum >= nCompCut)
								nConditionOK = 0;
						}

						if(nConditionOK == 0)
							break;
					}

					/* write cs */
					if(nConditionOK == 1)
					{
						for(nl=0; nl<nWindowSize; nl++)
						{
							pCS->pMatElement[ni+nl] = 255;
						}
					}
				}
			}
					
			/* write to files */
			pCSF = NULL;
			pCSF = CreateByteMatrix(1, nScanLen);
			if(pCSF == NULL)
			{
				printf("Error: cannot create buffer for writing conservation score!\n");
			}

			nTarStrand = IMGETAT(pPosInfo, nTarid, 3);
			nTarOffset = IMGETAT(pPosInfo, nTarid, 4)-1;
			nTempChr = IMGETAT(pPosInfo, nTarid, 5);


			if( (nTempChr >= 0) && (nTempChr<pTargetChrSize->nHeight) )
			{
				/* + strand */
				if(nTarStrand == 0)
				{
					nl = IMGETAT(pPosInfo, nTarid, 1);
					

					if(fseek(vfpScore[nTempChr], nl, SEEK_SET) != 0)
					{
						printf("Error: fseek error!\n");
						exit(EXIT_FAILURE);
					}

					nx = 0;
					for(nj=0; nj<nScanLen; nj++)
					{
						if(vAln[nTarid]->m_pString[nj] != '-')
						{
							pCSF->pMatElement[nx] = pCS->pMatElement[nj];
							nx++;
							nl++;
						}	
					}
					if(nl != IMGETAT(pPosInfo, nTarid, 2))
					{
						printf("Error: coordinate wrong!\n");
						exit(EXIT_FAILURE);
					}

					

					numwritten = fwrite(pCSF->pMatElement, sizeof(unsigned char), nx, vfpScore[nTempChr]);
					if(numwritten != nx)
					{
						printf("Error: coordinate wrong!\n");
						exit(EXIT_FAILURE);
					}
				}
				/* - strand */
				else
				{
					nl = nTarOffset - IMGETAT(pPosInfo, nTarid, 2) + 1;
					nTempChr = IMGETAT(pPosInfo, nTarid, 5);

					if(fseek(vfpScore[nTempChr], nl, SEEK_SET) != 0)
					{
						printf("Error: fseek error!\n");
						exit(EXIT_FAILURE);
					}
					
					nx = 0;
					for(nj=(nScanLen-1); nj>=0; nj--)
					{
						if(vAln[nTarid]->m_pString[nj] != '-')
						{
							pCSF->pMatElement[nx] = pCS->pMatElement[nj];
							nx++;
							nl++;
						}	
					}
					if(nl != (nTarOffset - IMGETAT(pPosInfo, nTarid, 1) + 1) )
					{
						printf("Error: coordinate wrong!\n");
						exit(EXIT_FAILURE);
					}

					numwritten = fwrite(pCSF->pMatElement, sizeof(unsigned char), nx, vfpScore[nTempChr]);
					if(numwritten != nx)
					{
						printf("Error: coordinate wrong!\n");
						exit(EXIT_FAILURE);
					}
				}
			}

			DestroyByteMatrix(pCSF);
	
			/* release memory */
			DestroyByteMatrix(pCS);
		}

		/* error */
		else
		{
			printf("fatal error: loading error!\n");
			exit(EXIT_FAILURE);
			/* continue; */
		}
	}

	fclose(fpIn);


	/* clear memory */
	DestroyIntMatrix(pPosInfo);
	
	for(ni=0; ni<nSpeciesNum; ni++)
	{
		DeleteString(vAln[ni]);
		vAln[ni] = NULL;
	}
	free(vAln);

	
	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Genome_CS_Get_Distribution:                                            */
/*  get distribution of conservation score.                                */
/*  return PROC_SUCCESS if successful.                                     */
/* ----------------------------------------------------------------------- */ 
int Genome_CS_Get_Distribution(char strCSPath[], char strFileListName[],
							   char strChrLenPath[], char strOutPath[],
							   char strSpecies[])
{
	/* define */
	struct DOUBLEMATRIX *pStat;
	struct INTMATRIX *pChrSize;
	struct BYTEMATRIX *pBuffer;
	char strFileName[LINE_LENGTH];
	char strLine[LINE_LENGTH];
	FILE *fpIn;
	FILE *fpCS;
	int nChrId,nChrLen;
	int nrlen;
	int nbuffersize;
	int numread;
	int ni,nj;

	/* init */
	nbuffersize = 100000;

	pStat = NULL;
	pStat = CreateDoubleMatrix(256, 1);
	if(pStat == NULL)
	{
		printf("Error: cannot allocate enough memory for statistics matrix!\n");
		return PROC_FAILURE;
	}
	pChrSize = NULL;
	pChrSize = IMLOAD(strChrLenPath);
	if(pChrSize == NULL)
	{
		printf("Error: cannot load chromosome size!\n");
		return PROC_FAILURE;
	}

	/* load file */
	sprintf(strFileName, "%s", strFileListName);
	fpIn = NULL;
	fpIn = fopen(strFileName, "rt");
	if(fpIn == NULL)
	{
		DestroyDoubleMatrix(pStat);
		DestroyIntMatrix(pChrSize);
		printf("Error: cannot open file!\n");
		return PROC_FAILURE;
	}

	/* create buffer */
	pBuffer = NULL;
	pBuffer = CreateByteMatrix(1, nbuffersize);
	if(pBuffer == NULL)
	{
		printf("Error: cannot create buffer!\n");
		return PROC_FAILURE;
	}

	/* process chr by chr */
	while(fgets(strLine, LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		nChrId = Genome_ChromosomeName_To_Index(strLine, strSpecies)-1;
		nChrLen = pChrSize->pMatElement[nChrId];

		sprintf(strFileName, "%s%s.cs", strCSPath, strLine);
		fpCS = NULL;
		fpCS = fopen(strFileName, "rb");
		if(fpCS == NULL)
		{
			printf("Error: cannot open *.cs file!\n");
			continue;
		}

		/* load */
		nrlen = nChrLen;
		while(nrlen >= nbuffersize)
		{
			numread = fread( pBuffer->pMatElement, sizeof( unsigned char ), nbuffersize, fpCS );
			if(numread != nbuffersize)
			{
				printf("Error: load conservation score error!\n");
				exit(EXIT_FAILURE);
			}

			for(ni=0; ni<numread; ni++)
			{
				nj = (int)(pBuffer->pMatElement[ni]);
				pStat->pMatElement[nj] = pStat->pMatElement[nj]+1.0;
			}
			nrlen -= nbuffersize;
		}

		if(nrlen > 0)
		{
			numread = fread( pBuffer->pMatElement, sizeof( unsigned char ), nrlen, fpCS );
			if(numread != nrlen)
			{
				printf("Error: load conservation score error!\n");
				exit(EXIT_FAILURE);
			}
			for(ni=0; ni<numread; ni++)
			{
				nj = (int)(pBuffer->pMatElement[ni]);
				pStat->pMatElement[nj] = pStat->pMatElement[nj]+1.0;
			}
		}

		fclose(fpCS);

	}

	/* close file */
	fclose(fpIn);

	/* save */
	DMSAVE(pStat, strOutPath);
	DestroyIntMatrix(pChrSize);
	DestroyDoubleMatrix(pStat);
	DestroyByteMatrix(pBuffer);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Genome_GetGCCS_Summary_Main()                                          */
/*  Get GC content and conservation score distributions for specified      */
/*  regions.                                                               */
/* ----------------------------------------------------------------------- */ 
int Genome_GetGCCS_Summary_Main(char strGenomePath[], char strCodPath[],
					char strOutputPath[], int nUseCS, double dC, 
					char strCSPath[])
{
	/* sequences and conservation */
	/* all sequences */
	struct FLEXSEQMOTIF *pSeqMtf = NULL;
	unsigned char *pBase;
	unsigned char *pCS;
	/* sequence number */
	int nSeqCount = 0;
	/* A,C,G,T count */
	struct DOUBLEMATRIX *pBaseCount;
	/* conservation score matrix */
	struct DOUBLEMATRIX *pCSCount;
	
	/* number of base types */
	int nBaseTypeNum = 4;
	double dTotalBase;
	

	/* motif consensus */
	FILE *fpIn;
	FILE *fpOut;
	char strLine[LONG_LINE_LENGTH];
	char strSeqAlias[LINE_LENGTH];

	/* coordinates */
	char strChr[LINE_LENGTH];
	int nStart,nEnd;
	int nM1,nM2;
	int nActualStart,nActualEnd;
	
	/* other variables */
	int ni,nj;
	
	/* #################################### */
	/* initialize                           */
	/* #################################### */
	pBaseCount = NULL;
	pBaseCount = CreateDoubleMatrix(1, nBaseTypeNum);
	if(pBaseCount == NULL)
	{
		printf("Error: Genome_GetGCCS_Summary_Main, cannot create matrix for counting\n");
		exit(EXIT_FAILURE);
	}
	if(nUseCS == 1)
	{
		pCSCount = NULL;
		pCSCount = CreateDoubleMatrix(256,1);
		if(pCSCount == NULL)
		{
			printf("Error: Genome_GetGCCS_Summary_Main, cannot create matrix for counting\n");
			exit(EXIT_FAILURE);
		}
	}
	
	/* #################################### */
	/* load initial parameters              */
	/* #################################### */
	AdjustDirectoryPath(strGenomePath);
	if(nUseCS == 1)
	{
		AdjustDirectoryPath(strCSPath);
	}

	/* #################################### */
	/* load sequences and conservation      */
	/* #################################### */
	fpIn = NULL;
	fpIn = fopen(strCodPath, "r");
	if(fpIn == NULL)
	{
        printf("Error: MotifMap_ScanMatrix_Genome_Main, cannot open coordinates file!\n");
		exit(EXIT_FAILURE);
	}

	/* read */
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%s %s %d %d", strSeqAlias, strChr,
			&nStart, &nEnd);
		
		/* process segment by segment */
		nM1 = nStart;
		while(nM1 <= nEnd)
		{
			nM2 = nM1+GENOME_CONTIG_LEN-1;
			if(nM2 > nEnd)
				nM2 = nEnd;

			nActualStart = nM1;
			nActualEnd = nM2;

			/* #################################### */ 
			/* load sequence and score              */
			/* #################################### */
			pSeqMtf = NULL;
			pSeqMtf = FLEXSEQMTFLOADSEQ_FOR_COUNTMCBACKGROUND_GENOME(nSeqCount,
				strGenomePath, strChr, nActualStart, nActualEnd, 
				nUseCS, strCSPath);
			if(pSeqMtf == NULL)
			{
				break;
			}

			/* #################################### */
			/* scan gc                              */
			/* #################################### */
			pBase = pSeqMtf->vSeq[0]->pMatElement;
			if(nUseCS == 1)
			{
				pCS = pSeqMtf->vScore[0]->pMatElement;
			}
			for(ni=0; ni<pSeqMtf->nSeqLen; ni++)
			{
				if(nUseCS == 1)
				{
					if((double)(pCS[ni]) < dC)
						continue;
					nj = pCS[ni];
					pCSCount->pMatElement[nj] = pCSCount->pMatElement[nj]+1.0;
				}

				switch(pBase[ni])
				{
					case 0: pBaseCount->pMatElement[0] += 1.0;
						break;
					case 1: pBaseCount->pMatElement[1] += 1.0;
						break;
					case 2: pBaseCount->pMatElement[2] += 1.0;
						break;
					case 3: pBaseCount->pMatElement[3] += 1.0;
						break;
				}
			}
			

			/* #################################### */
			/* release memory                       */
			/* #################################### */
			FLEXSEQMTFDESTROY(pSeqMtf);
			nM1 = nM2+1;
		}

		nSeqCount++;
	}

	/* close files */
	fclose(fpIn);
	

	/* save gc */
	dTotalBase = 0.0;
	for(ni=0; ni<nBaseTypeNum; ni++)
	{
		dTotalBase += pBaseCount->pMatElement[ni];
	}

	if(dTotalBase < 1e-6)
	{
		dTotalBase = 1e-6;
	}
	
	sprintf(strLine, "%s_gc.stat", strOutputPath);
	fpOut = NULL;
	fpOut = fopen(strLine, "w");
	if(fpOut == NULL)
	{
        printf("Error: MotifMap_ScanMatrix_Genome_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nBaseTypeNum; ni++)
	{
		pBaseCount->pMatElement[ni] = pBaseCount->pMatElement[ni]/dTotalBase;
		fprintf(fpOut, "%f\t", pBaseCount->pMatElement[ni]);
	}

	fprintf(fpOut, "\n");
	fprintf(fpOut, "GC= %f", (pBaseCount->pMatElement[1]+pBaseCount->pMatElement[2]));

	fclose(fpOut);

	/* #################################### */
	/* release memory                       */
	/* #################################### */
	DestroyDoubleMatrix(pBaseCount);


	/* save conservation */
	if(nUseCS == 1)
	{
		sprintf(strLine, "%s_cs.stat", strOutputPath);
		DMSAVE(pCSCount, strLine);
		DestroyDoubleMatrix(pCSCount);
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Genome_GetRegionCS_Summary_Main()                                      */
/*  Get conservation score summaries for each region specified in the      */
/*  input file.                                                            */
/* ----------------------------------------------------------------------- */ 
int Genome_GetRegionCS_Summary_Main(char strGenomePath[], char strCodPath[],
					char strOutputPath[], int nUseCS, double dC, 
					char strCSPath[])
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	char strLine[LONG_LINE_LENGTH];
	char strSeqAlias[LINE_LENGTH];
	char strChr[LINE_LENGTH];
	int nStart,nEnd;
	int nTotPos;
	double dTotScore,dMeanScore,dSD;
	
	/* #################################### */
	/* initialize                           */
	/* #################################### */
	
	/* #################################### */
	/* load initial parameters              */
	/* #################################### */
	AdjustDirectoryPath(strGenomePath);
	AdjustDirectoryPath(strCSPath);
	
	/* #################################### */
	/* load regions                         */
	/* #################################### */
	fpIn = NULL;
	fpIn = fopen(strCodPath, "r");
	if(fpIn == NULL)
	{
        printf("Error: Genome_GetRegionCS_Summary_Main, cannot open coordinates file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: Genome_GetRegionCS_Summary_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "#seq_id\tchromosome\tstart\tend\tstrand\tregion_length\t#_conserved_non-repeat_positions\tsum_of_csscore_for_conserved_non-repeat_positions\tmean_csscore_for_conserved_non-repeat_positions\tSD_of_mean_csscore_for_conserved_non-repeat_positions\n");

	/* read */
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		sscanf(strLine, "%s %s %d %d", strSeqAlias, strChr,
			&nStart, &nEnd);

		nTotPos = 0;
		dTotScore = 0.0;
		dMeanScore = 0.0;
		dSD = 0.0;

		/* compute region level summary */
		Genome_GetRegionCS_Summary(strGenomePath, strCSPath, strChr, nStart, nEnd, 
					dC, &nTotPos, &dTotScore, &dMeanScore, &dSD);

		fprintf(fpOut, "%s\t%s\t%d\t%d\t+\t%d\t%d\t%9.7e\t%9.7e\t%9.7e\n", 
			strSeqAlias, strChr, nStart, nEnd, nEnd-nStart+1, 
			nTotPos, dTotScore, dMeanScore, dSD);
	}

	/* close files */
	fclose(fpIn);
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Genome_GetRegionCS_Summary()                                           */
/*  Get conservation score summaries for a genomic region.                 */
/* ----------------------------------------------------------------------- */ 
int Genome_GetRegionCS_Summary(char strGenomePath[], char strCSPath[], 
					char strChr[], int nStart, int nEnd, double dC, 
					int *pnTotPos, double *pdTotScore, double *pdMeanScore, 
					double *pdSD)
{
	/* define */
	/* sequences and conservation */
	struct FLEXSEQMOTIF *pSeqMtf = NULL;
	unsigned char *pBase;
	unsigned char *pCS;
	/* sequence number */
	int nSeqCount = 0;
	/* conservation score matrix */
	struct DOUBLEMATRIX *pCSCount;
	int nUseCS = 1;
	
	/* number of base types */
	int nBaseTypeNum = 4;
	
	/* coordinates */
	int nM1,nM2;
	int nActualStart,nActualEnd;
	
	/* other variables */
	int ni,nj;
	
	/* #################################### */
	/* initialize                           */
	/* #################################### */
	pCSCount = NULL;
	pCSCount = CreateDoubleMatrix(256,1);
	if(pCSCount == NULL)
	{
		printf("Error: Genome_GetGCCS_Summary_Main, cannot create matrix for counting\n");
		exit(EXIT_FAILURE);
	}
	
	/* #################################### */
	/* load initial parameters              */
	/* #################################### */
	AdjustDirectoryPath(strGenomePath);
	AdjustDirectoryPath(strCSPath);

	/* #################################### */
	/* load sequences and conservation      */
	/* #################################### */
	nM1 = nStart;
	while(nM1 <= nEnd)
	{
		nM2 = nM1+GENOME_CONTIG_LEN-1;
		if(nM2 > nEnd)
			nM2 = nEnd;

		nActualStart = nM1;
		nActualEnd = nM2;

		/* #################################### */ 
		/* load sequence and score              */
		/* #################################### */
		pSeqMtf = NULL;
		pSeqMtf = FLEXSEQMTFLOADSEQ_FOR_COUNTMCBACKGROUND_GENOME(nSeqCount,
				strGenomePath, strChr, nActualStart, nActualEnd, 
				nUseCS, strCSPath);
		if(pSeqMtf == NULL)
		{
			break;
		}

		/* #################################### */
		/* scan gc                              */
		/* #################################### */
		pBase = pSeqMtf->vSeq[0]->pMatElement;
		pCS = pSeqMtf->vScore[0]->pMatElement;
		for(ni=0; ni<pSeqMtf->nSeqLen; ni++)
		{
			if((double)(pCS[ni]) < dC)
				continue;
			if(pBase[ni] <= 3)
			{				
				nj = pCS[ni];
                pCSCount->pMatElement[nj] = pCSCount->pMatElement[nj]+1.0;
			}
		}

		/* #################################### */
		/* release memory                       */
		/* #################################### */
		FLEXSEQMTFDESTROY(pSeqMtf);
		nM1 = nM2+1;

		nSeqCount++;
	}

	/* #################################### */
	/* count cs                             */
	/* #################################### */
	*pnTotPos = 0;
	*pdTotScore = 0.0;
	*pdMeanScore = 0.0;
    *pdSD = 0.0;
	for(ni=0; ni<256; ni++)
	{
		if((double)ni >= dC)
		{
			*pnTotPos = *pnTotPos+(int)(pCSCount->pMatElement[ni]);
			*pdTotScore = *pdTotScore+ni*pCSCount->pMatElement[ni];
		}
	}
	if(*pnTotPos > 0)
	{
		*pdMeanScore = (*pdTotScore)/(*pnTotPos);
	}
	if(*pnTotPos > 1)
	{
		for(ni=0; ni<256; ni++)
		{
			if((double)ni >= dC)
			{
				*pdSD = *pdSD+(ni-(*pdMeanScore))*(ni-(*pdMeanScore))*(pCSCount->pMatElement[ni]);
			}
		}
		*pdSD = *pdSD/(*pnTotPos-1);
		*pdSD = sqrt(*pdSD);
	}
	
	/* #################################### */
	/* release memory                       */
	/* #################################### */
	DestroyDoubleMatrix(pCSCount);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Genome_GetMarkovBGForRegion()                                          */
/*  Get markovian background for specific genomic regions.                 */
/* ----------------------------------------------------------------------- */ 
int Genome_GetMarkovBGForRegion(char strGenomePath[], char strGenomicTargetFile[],
								char strSpecies[], int nMCOrder, char strOutFile[])
{
	/* define */
	FILE *fpSeq;
	FILE *fpTarget;
	char strLine[LINE_LENGTH];
	char strChrName[LINE_LENGTH];
	int nChrId;
	char strSeqFile[LINE_LENGTH];
	char strBGFFile[LINE_LENGTH];
	char strBGRFile[LINE_LENGTH];
	int nStart;
	int nEnd;
	int nP1,nP2;
	int npp1,npp2;
	int nRemainLen;
	int nBatchLen;
	int ni, nj, nstartoff, nendoff;
	int nHeight,nWidth;

	struct INTMATRIX *pChrLen;
	struct SEQMOTIF *pSeqMtf;
	struct DOUBLEMATRIX *pBGF;
	struct DOUBLEMATRIX *pBGR;
	double *pEle;
	double dTemp,dTotal;


	/* for testing */
	/* int nFastaCount;
	FILE *fpOut;
	fpOut = NULL;
	fpOut = fopen("testout.txt", "wt");
	if(fpOut == NULL)
	{
		exit(EXIT_FAILURE);
	} */

	/* init space */
	if(nMCOrder < 0)
	{
		printf("Error: Genome_GetMarkovBGForRegion, the background order should be >= 0!\n");
		exit(EXIT_FAILURE);
	}

	/* init backgroud matrix */
	nWidth = 4;
	nHeight = (int)(pow((double)nWidth, (double)nMCOrder));
	sprintf(strBGFFile, "%s_%d_f.txt", strOutFile, nMCOrder);
	sprintf(strBGRFile, "%s_%d_r.txt", strOutFile, nMCOrder);

	pBGF = NULL;
	pBGF = CreateDoubleMatrix(nHeight, nWidth);
	if(pBGF == NULL)
	{
		printf("Error: Genome_GetMarkovBGForRegion, cannot create forward background matrix!\n");
		exit(EXIT_FAILURE);
	}
	pEle = pBGF->pMatElement;
	for(ni=0; ni<nHeight; ni++)
	{
		for(nj=0; nj<nWidth; nj++)
		{
			*pEle = 1.0;
			pEle++;
		}
	}

	pBGR = NULL;
	pBGR = CreateDoubleMatrix(nHeight, nWidth);
	if(pBGR == NULL)
	{
		printf("Error: Genome_GetMarkovBGForRegion, cannot create reverse background matrix!\n");
		exit(EXIT_FAILURE);
	}
	pEle = pBGR->pMatElement;
	for(ni=0; ni<nHeight; ni++)
	{
		for(nj=0; nj<nWidth; nj++)
		{
			*pEle = 1.0;
			pEle++;
		}
	}

	/* init genome size */
	sprintf(strLine, "%schrlen.txt", strGenomePath);
	pChrLen = NULL;
	pChrLen = IMLOAD(strLine);
	if(pChrLen == NULL)
	{
		printf("Error: Genome_GetMarkovBGForRegion, cannot load chromosome length!\n");
		exit(EXIT_FAILURE);
	}

	pSeqMtf = NULL;
	pSeqMtf = SEQMTFCREATE(0, 1, 0, 0);
	if(pSeqMtf == NULL)
	{
		printf("Error: Genome_GetMarkovBGForRegion, cannot create sequence/motif complex!\n");
		exit(EXIT_FAILURE);
	}
	SEQMTFCREATESEQ(pSeqMtf, 0, (GENOME_CONTIG_LEN+2*nMCOrder));
	
	/* load candidate regions one by one and process them */
	fpTarget = NULL;
	fpTarget = fopen(strGenomicTargetFile, "rt");
	if(fpTarget == NULL)
	{
		printf("Error: Genome_GetMarkovBGForRegion, cannot open target region file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LINE_LENGTH, fpTarget) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		/* load the region */
		sscanf(strLine, "%s %d %d", strChrName, &nStart, &nEnd);
		nChrId = Genome_ChromosomeName_To_Index(strChrName, strSpecies);
		sprintf(strSeqFile, "%s%s.sq", strGenomePath, strChrName);
	
		/* open sequence file and get the sequence */
		fpSeq = NULL;
		fpSeq = fopen(strSeqFile, "rb");
		if(fpSeq == NULL)
		{
			printf("Warning: Genome_GetMarkovBGForRegion, cannot open sequence file!\n");
			continue;
		}

		/* process contig by contig */
		nP1 = nStart;
		nP2 = nP1+GENOME_CONTIG_LEN-1;
		nRemainLen = nEnd-nStart+1;
		if(nRemainLen < GENOME_CONTIG_LEN)
		{
			nP2 = nP1+nRemainLen-1;
		}

		while(nRemainLen > 0)
		{

			/* load seq */
			npp1 = nP1-nMCOrder;
			npp2 = nP2+nMCOrder;
			
			/* for test */
			/* fprintf(fpOut, "\n>%s:%d-%d\n", strChrName, npp1, npp2);
			nFastaCount = 0;
			*/

			nstartoff = 0;
			nendoff = 0;
			nBatchLen = npp2-npp1+1;
			if(npp1 < 0)
			{
				nstartoff = -npp1;
				npp1 = 0;
			}
			if(npp2 >= pChrLen->pMatElement[nChrId-1])
			{
				nendoff = npp2-pChrLen->pMatElement[nChrId-1]+1;
				npp2 = pChrLen->pMatElement[nChrId-1]-1;
			}

			for(ni=0; ni<nstartoff; ni++)
			{
				(pSeqMtf->ppSeq)[0]->pMatElement[ni] = 4;
			}
			for(ni=0; ni<nendoff; ni++)
			{
				(pSeqMtf->ppSeq)[0]->pMatElement[nBatchLen-1-ni] = 4;
			}

			SEQMTF_LOADSEQFROMGENOME(pSeqMtf, 0, nstartoff, fpSeq, npp1, npp2, '+');
			SEQMTF_ADDCOUNTTONUCLEICBGMC_BOTHSTRAND(pBGF, pBGR, 
											pSeqMtf, 0, nBatchLen,
											nMCOrder, (nBatchLen-1-nMCOrder), nMCOrder);

			/* for test */
			/* for(ni=0; ni<nBatchLen; ni++)
			{
				switch(pSeqMtf->ppSeq[0]->pMatElement[ni])
				{
					case 0: fprintf(fpOut, "A");
						break;
					case 1: fprintf(fpOut, "C");
						break;
					case 2: fprintf(fpOut, "G");
						break;
					case 3: fprintf(fpOut, "T");
						break;
					default: fprintf(fpOut, "N");
				}
				nFastaCount++;
				if(nFastaCount == FASTA_LINE_LEN)
				{
					fprintf(fpOut, "\n");
					nFastaCount = 0;
				}
			} */

			/* count processed bases */
			nRemainLen -= (nP2-nP1+1);

			/* get next */
			nP1 = nP2+1;
			if(nRemainLen >= GENOME_CONTIG_LEN)
			{
				nP2 = nP1+GENOME_CONTIG_LEN-1;
			}
			else
			{
				nP2 = nP1+nRemainLen-1;
			}
		}
	
		
		/* */

		/* close sequence file */
		fclose(fpSeq);
	}


	/* close target file */
	fclose(fpTarget);

	/* for testing */
	/* fclose(fpOut); */

	/* normalize */
	for(ni=0; ni<pBGF->nHeight; ni++)
	{
		dTotal = 0.0;
		for(nj=0; nj<pBGF->nWidth; nj++)
		{
			dTotal += DMGETAT(pBGF, ni, nj);
		}
		for(nj=0; nj<pBGF->nWidth; nj++)
		{
			dTemp = DMGETAT(pBGF, ni, nj)/dTotal;
			DMSETAT(pBGF, ni, nj, dTemp);
		}
	}

	for(ni=0; ni<pBGR->nHeight; ni++)
	{
		dTotal = 0.0;
		for(nj=0; nj<pBGR->nWidth; nj++)
		{
			dTotal += DMGETAT(pBGR, ni, nj);
		}
		for(nj=0; nj<pBGR->nWidth; nj++)
		{
			dTemp = DMGETAT(pBGR, ni, nj)/dTotal;
			DMSETAT(pBGR, ni, nj, dTemp);
		}
	}


	/* output results */
	DMSAVE(pBGF, strBGFFile);
	DMSAVE(pBGR, strBGRFile);

	/* release memory */
	SEQMTFDESTROY(pSeqMtf);
	DestroyIntMatrix(pChrLen);
	DestroyDoubleMatrix(pBGF);
	DestroyDoubleMatrix(pBGR);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Genome_MapMotifToRegion()                                              */
/*  Map motif to specific genomic regions.                                 */
/* ----------------------------------------------------------------------- */ 
int Genome_MapMotifToRegion(char strGenomePath[], char strGenomicTargetFile[], 
							char strSpecies[], char strVersion[],
							char strMotifPath[], char strMotifList[], 
							int nBGOrder, char strMarkovBGFile[], 
							char strLikehoodRatioCutoffList[],
							char strOutPath[])
{
	/* define */

	/* genome */
	FILE *fpSeq;
	FILE *fpTarget;
	char strChrName[LINE_LENGTH];
	int nChrId;
	char strSeqFile[LINE_LENGTH];
	struct INTMATRIX *pChrLen;

	/* background */
	char strBGFFile[LINE_LENGTH];
	char strBGRFile[LINE_LENGTH];
	struct DOUBLEMATRIX *pBGF;
	struct DOUBLEMATRIX *pBGR;
	struct DOUBLEMATRIX *pCutoff;
	double dCutoff;

	/* motif */
	char strMotifFile[LINE_LENGTH];
	struct DOUBLEMATRIX **vMotif;
	struct tagString **vMotifName;
	int nMotifNum,nMotifId;
	FILE *fpMotifList;
	int nMaxMotifLen;

	/* general */
	char strLine[LINE_LENGTH];
	int nStart;
	int nEnd;
	int nP1,nP2;
	int npp1,npp2;
	int nWriteOffset;
	int nInitLen;
	int nRemainLen;
	int nBatchLen;
	int ni, nj, nstartoff, nendoff;
	int nHeight,nWidth;

	/* contig */
	struct SEQMOTIF *pSeqMtf;
	double *pEle,*pEle2;

	/* output */
	FILE **vfpMapOut;
	double dTotal;

	/* for testing */
	/* int nFastaCount;
	FILE *fpOut;
	fpOut = NULL;
	fpOut = fopen("testout.txt", "wt");
	if(fpOut == NULL)
	{
		exit(EXIT_FAILURE);
	} */

	/* ---------------------------------- */
	/* init backgroud matrix and take log */
	/* ---------------------------------- */
	nWidth = 4;
	nHeight = (int)(pow((double)nWidth, (double)nBGOrder));
	
	sprintf(strBGFFile, "%s_%d_f.txt", strMarkovBGFile, nBGOrder);
	sprintf(strBGRFile, "%s_%d_r.txt", strMarkovBGFile, nBGOrder);
	pBGF = NULL;
	pBGF = DMLOAD(strBGFFile);
	if(pBGF == NULL)
	{
		printf("Error: Genome_MapMotifToRegion, cannot load background markov matrix!\n");
		exit(EXIT_FAILURE);
	}
	if((pBGF->nHeight != nHeight) || (pBGF->nWidth != nWidth))
	{
		printf("Error: Genome_MapMotifToRegion, matrix dimension not match!\n");
		exit(EXIT_FAILURE);
	}
	pEle = pBGF->pMatElement;
	for(ni=0; ni<nHeight; ni++)
	{
		for(nj=0; nj<nWidth; nj++)
		{
			if(*pEle <= 0.0)
			{
				printf("Error: Genome_MapMotifToRegion, trying to log(negative number)!\n");
				exit(EXIT_FAILURE);
			}
			*pEle = log(*pEle);
			pEle++;
		}
	}


	pBGR = NULL;
	pBGR = DMLOAD(strBGRFile);
	if(pBGR == NULL)
	{
		printf("Error: Genome_MapMotifToRegion, cannot load background markov matrix!\n");
		exit(EXIT_FAILURE);
	}
	if((pBGR->nHeight != nHeight) || (pBGR->nWidth != nWidth))
	{
		printf("Error: Genome_MapMotifToRegion, matrix dimension not match!\n");
		exit(EXIT_FAILURE);
	}
	pEle = pBGR->pMatElement;
	for(ni=0; ni<nHeight; ni++)
	{
		for(nj=0; nj<nWidth; nj++)
		{
			if(*pEle <= 0.0)
			{
				printf("Error: Genome_MapMotifToRegion, trying to log(negative number)!\n");
				exit(EXIT_FAILURE);
			}
			*pEle = log(*pEle);
			pEle++;
		}
	}


	/* ---------------------------------- */
	/* init motif matrix and take log     */
	/* ---------------------------------- */
	/* get motif number */
	nMotifNum = 0;
	sprintf(strMotifFile, "%s%s", strMotifPath, strMotifList);
	fpMotifList = NULL;
	fpMotifList = fopen(strMotifFile, "rt");
	if(fpMotifList == NULL)
	{
		printf("Error: Genome_MapMotifToRegion, cannot open motif list!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LINE_LENGTH, fpMotifList) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		nMotifNum++;
	}

	fclose(fpMotifList);

	/* init memory for motifs */
	if(nMotifNum <= 0)
	{
		printf("Warning: Genome_MapMotifToRegion, no motifs!\n");
		return PROC_FAILURE;
	}

	sprintf(strMotifFile, "%s%s", strMotifPath, strLikehoodRatioCutoffList);
	pCutoff = NULL; 
	pCutoff = DMLOAD(strMotifFile);
	if(pCutoff == NULL)
	{
		printf("Error: Genome_MapMotifToRegion, cannot load cutoff for motifs!\n");
		exit(EXIT_FAILURE);
	}
	if(pCutoff->nHeight != nMotifNum)
	{
		printf("Error: Genome_MapMotifToRegion, cutoff number and motif number not match!\n");
		exit(EXIT_FAILURE);
	}

	vMotif = NULL;
	vMotif = (struct DOUBLEMATRIX **)calloc(nMotifNum, sizeof(struct DOUBLEMATRIX *));
	if(vMotif == NULL)
	{
		printf("Error: Genome_MapMotifToRegion, cannot create memory for motifs!\n");
		exit(EXIT_FAILURE);
	}

	vMotifName = NULL;
	vMotifName = (struct tagString **)calloc(nMotifNum, sizeof(struct tagString*));
	if(vMotifName == NULL)
	{
		printf("Error: Genome_MapMotifToRegion, cannot create memory for motif names!\n");
		exit(EXIT_FAILURE);
	}

	vfpMapOut = NULL;
	vfpMapOut = (FILE **)calloc(nMotifNum, sizeof(FILE *));
	if(vfpMapOut == NULL)
	{
		printf("Error: Genome_MapMotifToRegion, cannot create memory for tracking map files!\n");
		exit(EXIT_FAILURE);
	}

	/* load motif matrices */
	nMaxMotifLen = 0;
	sprintf(strMotifFile, "%s%s", strMotifPath, strMotifList);
	fpMotifList = NULL;
	fpMotifList = fopen(strMotifFile, "rt");
	if(fpMotifList == NULL)
	{
		printf("Error: Genome_MapMotifToRegion, cannot open motif list!\n");
		exit(EXIT_FAILURE);
	}

	nMotifId = 0;
	while(fgets(strLine, LINE_LENGTH, fpMotifList) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		
		/* load motif matrix */
		sprintf(strMotifFile, "%s%s.txt", strMotifPath, strLine);
		vMotif[nMotifId] = DMLOAD(strMotifFile);
		if(vMotif[nMotifId] == NULL)
		{
			printf("Error: Genome_MapMotifToRegion, error in loading motif %s!\n", strLine);
			exit(EXIT_FAILURE);
		}

		if(vMotif[nMotifId]->nHeight > nMaxMotifLen)
		{
			nMaxMotifLen = vMotif[nMotifId]->nHeight;
		}

		pEle = vMotif[nMotifId]->pMatElement;
		for(ni=0; ni<vMotif[nMotifId]->nHeight; ni++)
		{
			pEle2 = pEle;
			dTotal = 0.0;
			for(nj=0; nj<vMotif[nMotifId]->nWidth; nj++)
			{
				*pEle2 += 0.5;
				dTotal += (*pEle2);
				pEle2++;
			}
			for(nj=0; nj<vMotif[nMotifId]->nWidth; nj++)
			{
				*pEle = log(*pEle/dTotal);
				pEle++;
			}
		}

		/* save motif name */
		StringAddTail((vMotifName+nMotifId), strLine);

		/* prepare motif file */
		sprintf(strMotifFile, "%s%s_%s_%s_map.bed", strOutPath, strLine, strSpecies, strVersion);
		vfpMapOut[nMotifId] = fopen(strMotifFile, "wt");
		if(vfpMapOut[nMotifId] == NULL)
		{
			printf("Error: Genome_MapMotifToRegion, cannot open motif mapping file to archive mapping results!\n");
			exit(EXIT_FAILURE);
		}
		fprintf(vfpMapOut[nMotifId], "browser position chr10:127416063-127516062\n");
		fprintf(vfpMapOut[nMotifId], "track name=TFBS_%s description=\"%s\"\n", strLine, strLine);

		/* get next index */
		nMotifId++;
	}
	fclose(fpMotifList);
	
	if(nMotifId != nMotifNum)
	{
		printf("Error: Genome_MapMotifToRegion, motif number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* ---------------------------------- */
	/* init genome size                   */
	/* ---------------------------------- */
	sprintf(strLine, "%schrlen.txt", strGenomePath);
	pChrLen = NULL;
	pChrLen = IMLOAD(strLine);
	if(pChrLen == NULL)
	{
		printf("Error: Genome_MapMotifToRegion, cannot load chromosome length!\n");
		exit(EXIT_FAILURE);
	}

	/* ---------------------------------- */
	/* init contig space                  */
	/* ---------------------------------- */
	pSeqMtf = NULL;
	nInitLen = GENOME_MED_CONTIG_LEN+2*nBGOrder+nMaxMotifLen-1;

	/* create 1 sequence, 2 background scores, 1 filter information */
	pSeqMtf = SEQMTFCREATE(0, 1, 2, 0);
	if(pSeqMtf == NULL)
	{
		printf("Error: Genome_MapMotifToRegion, cannot create sequence/motif complex!\n");
		exit(EXIT_FAILURE);
	}
	/* sequence */
	SEQMTFCREATESEQ(pSeqMtf, 0, nInitLen);
	/* forward mc background */
	SEQMTFCREATESCORE(pSeqMtf, 0, nInitLen);
	/* backward mc background */
	SEQMTFCREATESCORE(pSeqMtf, 1, nInitLen);
	/* filter score */
	/* SEQMTFCREATEPATH(pSeqMtf, 0, nInitLen); */
	
	/* -------------------------------------------------- */
	/* load candidate regions one by one and process them */
	/* -------------------------------------------------- */
	fpTarget = NULL;
	fpTarget = fopen(strGenomicTargetFile, "rt");
	if(fpTarget == NULL)
	{
		printf("Error: Genome_MapMotifToRegion, cannot open target region file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LINE_LENGTH, fpTarget) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		/* load the region */
		sscanf(strLine, "%s %d %d", strChrName, &nStart, &nEnd);
		printf("Processing %s:%d-%d...\n", strChrName, nStart, nEnd);
			
		nChrId = Genome_ChromosomeName_To_Index(strChrName, strSpecies);
		sprintf(strSeqFile, "%s%s.sq", strGenomePath, strChrName);
	
		/* open sequence file and get the sequence */
		fpSeq = NULL;
		fpSeq = fopen(strSeqFile, "rb");
		if(fpSeq == NULL)
		{
			printf("Warning: Genome_MapMotifToRegion, cannot open sequence file!\n");
			continue;
		}

		/* process contig by contig */
		nP1 = nStart;
		nP2 = nP1+GENOME_MED_CONTIG_LEN-1;
		nRemainLen = nEnd-nStart+1;
		if(nRemainLen < GENOME_MED_CONTIG_LEN)
		{
			nP2 = nP1+nRemainLen-1;
		}

		while(nRemainLen > 0)
		{

			/* load seq */
			npp1 = nP1-nBGOrder;
			npp2 = nP2+nBGOrder+nMaxMotifLen-1;
			nWriteOffset = npp1;
			
			/* for test */
			/* fprintf(fpOut, "\n>%s:%d-%d\n", strChrName, npp1, npp2);
			nFastaCount = 0;
			*/

			nstartoff = 0;
			nendoff = 0;
			nBatchLen = npp2-npp1+1;
			if(npp1 < 0)
			{
				nstartoff = -npp1;
				npp1 = 0;
			}
			if(npp2 >= pChrLen->pMatElement[nChrId-1])
			{
				nendoff = npp2-pChrLen->pMatElement[nChrId-1]+1;
				npp2 = pChrLen->pMatElement[nChrId-1]-1;
			}

			for(ni=0; ni<nstartoff; ni++)
			{
				(pSeqMtf->ppSeq)[0]->pMatElement[ni] = 4;
			}
			for(ni=0; ni<nendoff; ni++)
			{
				(pSeqMtf->ppSeq)[0]->pMatElement[nBatchLen-1-ni] = 4;
			}
			
			/* load sequence */
			SEQMTF_LOADSEQFROMGENOME(pSeqMtf, 0, nstartoff, fpSeq, npp1, npp2, '+');
			/* get markov background for every base */
			SEQMTF_SETBGMCSCORE_BOTHSTRAND(pSeqMtf, 0, nBatchLen, nBGOrder, (nBatchLen-1-nBGOrder), pBGF, pBGR, nBGOrder);
			/* call motifs */
			for(ni=0; ni<nMotifNum; ni++)
			{
				/* map */
				dCutoff = DMGETAT(pCutoff, ni, 0);
				SEQMTF_MAPMOTIF_BOTHSTRANDASBG(vMotif[ni], ni, dCutoff, pSeqMtf, 0, nBatchLen, nBGOrder, (nBatchLen-nBGOrder-nMaxMotifLen));

				/* write motifs */
				SEQMTF_WRITEMOTIFSITETOBED_THENCLEARSITES(pSeqMtf, nMotifNum, vMotifName, strChrName, nWriteOffset, vfpMapOut);
			}

			/* for test */
			/* for(ni=0; ni<nBatchLen; ni++)
			{
				switch(pSeqMtf->ppSeq[0]->pMatElement[ni])
				{
					case 0: fprintf(fpOut, "A");
						break;
					case 1: fprintf(fpOut, "C");
						break;
					case 2: fprintf(fpOut, "G");
						break;
					case 3: fprintf(fpOut, "T");
						break;
					default: fprintf(fpOut, "N");
				}
				nFastaCount++;
				if(nFastaCount == FASTA_LINE_LEN)
				{
					fprintf(fpOut, "\n");
					nFastaCount = 0;
				}
			} */

			/* count processed bases */
			nRemainLen -= (nP2-nP1+1);

			/* get next */
			nP1 = nP2+1;
			if(nRemainLen >= GENOME_MED_CONTIG_LEN)
			{
				nP2 = nP1+GENOME_MED_CONTIG_LEN-1;
			}
			else
			{
				nP2 = nP1+nRemainLen-1;
			}
		}
	
		
		/* */

		/* close sequence file */
		fclose(fpSeq);
	}


	/* close target file */
	fclose(fpTarget);

	/* for testing */
	/* fclose(fpOut); */

	/* -------------------------- */
	/* release memory             */
	/* -------------------------- */
	SEQMTFDESTROY(pSeqMtf);
	
	DestroyDoubleMatrix(pBGF);
	DestroyDoubleMatrix(pBGR);

	DestroyDoubleMatrix(pCutoff);

	for(ni=0; ni<nMotifNum; ni++)
	{
		DestroyDoubleMatrix(vMotif[ni]);
		vMotif[ni] = NULL;

		DeleteString(vMotifName[ni]);
		vMotifName[ni] = NULL;

		if(vfpMapOut[ni] != NULL)
		{
			fclose(vfpMapOut[ni]);
		}
	}
	free(vMotif);
	free(vMotifName);
	free(vfpMapOut);

	DestroyIntMatrix(pChrLen);

	/* return */
	return PROC_SUCCESS;
}

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
					   int nTypeIIMask, char strTypeIICutoffList[])
{
	/* define */
	FILE *fpMotifList;
	FILE *fpInitMap;
	FILE *fpFinalMap;
	char strLine[LINE_LENGTH];
	char strInitMapFile[LINE_LENGTH];
	char strFinalMapFile[LINE_LENGTH];
	char strChrName[LINE_LENGTH];
	int nStart,nEnd;
	char strAlias[LINE_LENGTH];
	double dScore;
	char chStrand;
	struct SEQMOTIF *pSeqMtf;
	int nInitLen;
	int nFiltered;
	int nMotifLen;
	int ni,nj,nk;
	double dCSTot;
	int nRepeatTot;
	double dRepeatRatio;
	unsigned char *pCSVec;
	unsigned char *pRepeatVec;
	unsigned char *pBaseVec;

	char strMotifFile[LINE_LENGTH];
	char strMotifListPath[LINE_LENGTH];
	char strTypeIICutoffPath[LINE_LENGTH];
	struct DOUBLEMATRIX *pTypeIICutoff;
	struct DOUBLEMATRIX *pLogPWM;
	double *pEle,*pEle2;
	double dTotal;
	int nMotifId;
	double dMotifScore;
	int nBadScore;


	/* init */
	sprintf(strMotifListPath, "%s%s", strMotifPath, strMotifList);
	sprintf(strTypeIICutoffPath, "%s%s", strMotifPath, strTypeIICutoffList);
	pTypeIICutoff = NULL;
	if(nTypeIIMask == 1)
	{
		pTypeIICutoff = DMLOAD(strTypeIICutoffPath);
		if(pTypeIICutoff == NULL)
		{
			printf("Error: Genome_FilterMotifSites, cannot load type II cutoff!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* create 1 sequence, 1 conservation score and 1 repeat filter information */
	pSeqMtf = NULL;
	nInitLen = MAX_MOTIF_LEN;
	pSeqMtf = SEQMTFCREATE(0, 1, 0, 2);
	if(pSeqMtf == NULL)
	{
		printf("Error: Genome_FilterMotifSites, cannot create sequence/motif complex!\n");
		exit(EXIT_FAILURE);
	}
	/* sequence */
	SEQMTFCREATESEQ(pSeqMtf, 0, nInitLen);
	/* conservation score */
	SEQMTFCREATEPATH(pSeqMtf, 0, nInitLen);
	/* repeat mask */
	SEQMTFCREATEPATH(pSeqMtf, 1, nInitLen);
	

	fpMotifList = NULL;
	fpMotifList = fopen(strMotifListPath, "rt");
	if(fpMotifList == NULL)
	{
		printf("Error: Genome_FilterMotifSites, cannot open target motif list!\n");
		exit(EXIT_FAILURE);
	}

	/* load motifs to be processed */
	nMotifId = 0;
	while(fgets(strLine, LINE_LENGTH, fpMotifList) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sprintf(strMotifFile, "%s%s.txt", strMotifPath, strLine);
		pLogPWM = NULL;
		pLogPWM = DMLOAD(strMotifFile);
		if(pLogPWM == NULL)
		{
			printf("Error: Genome_FilterMotifSites, cannot open motif PWM!\n");
			exit(EXIT_FAILURE);
		}
		pEle = pLogPWM->pMatElement;
		for(ni=0; ni<pLogPWM->nHeight; ni++)
		{
			pEle2 = pEle;
			dTotal = 0.0;
			for(nj=0; nj<pLogPWM->nWidth; nj++)
			{
				*pEle2 += 0.5;
				dTotal += (*pEle2);
				pEle2++;
			}
			for(nj=0; nj<pLogPWM->nWidth; nj++)
			{
				*pEle = log(*pEle/dTotal);
				pEle++;
			}
		}


		sprintf(strInitMapFile, "%s%s_%s_%s_map.bed", strInitMapPath, strLine, strSpecies, strVersion);
		sprintf(strFinalMapFile, "%s%s_%s_%s_map.bed", strFinalMapPath, strLine, strSpecies, strVersion);

		/* open init and final map file */
		fpInitMap = NULL;
		fpInitMap = fopen(strInitMapFile, "rt");
		if(fpInitMap == NULL)
		{
			printf("Error: Genome_FilterMotifSites, cannot open initial motif map!\n");
			exit(EXIT_FAILURE);
		}

		fpFinalMap = NULL;
		fpFinalMap = fopen(strFinalMapFile, "wt");
		if(fpFinalMap == NULL)
		{
			printf("Error: Genome_FilterMotifSites, cannot open filtered motif map!\n");
			exit(EXIT_FAILURE);
		}

		/* load head info */
		if(fgets(strLine, LINE_LENGTH, fpInitMap) == NULL)
		{
			fclose(fpInitMap);
			fclose(fpFinalMap);
			continue;
		}
		else
		{
			fprintf(fpFinalMap, "%s", strLine);
		}

		if(fgets(strLine, LINE_LENGTH, fpInitMap) == NULL)
		{
			fclose(fpInitMap);
			fclose(fpFinalMap);
			continue;
		}
		else
		{
			fprintf(fpFinalMap, "%s", strLine);
		}

		/* process motif site one by one */
		while(fgets(strLine, LINE_LENGTH, fpInitMap) != NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;

			sscanf(strLine, "%s %d %d %s %lf %c", strChrName, &nStart, &nEnd, strAlias, &dScore, &chStrand);
			
			/* load sequence and conservation score */
			SEQMTF_LOADSEQCSFROMGENOME_FORMOTIFFILTER(pSeqMtf, 0, strChrName, nStart, nEnd, 
				strGenomePath, strCScorePath);
			
			/* filter */
			nFiltered = 0;
			dCSTot = 0.0;
			dMotifScore = 0.0;
			nBadScore = 0;
			nRepeatTot = 0;
			pBaseVec = pSeqMtf->ppSeq[0]->pMatElement;
			pCSVec = pSeqMtf->ppSamplePath[0]->pMatElement;
			pRepeatVec = pSeqMtf->ppSamplePath[1]->pMatElement;
			nMotifLen = nEnd-nStart+1;
			if((nMotifLen <= 0) || (nMotifLen != pLogPWM->nHeight))
			{
				printf("Error: Genome_FilterMotifSites, motif site length <= 0 or not match!\n");
				exit(EXIT_FAILURE);
			}

			/* get statistics */
			if(nTypeIIMask == 1)
			{
				if(chStrand == '-')
				{
					for(ni=0; ni<nMotifLen; ni++)
					{
						dCSTot += (double)pCSVec[ni];
						nRepeatTot += (int)pRepeatVec[ni];
						nk = (int)pBaseVec[ni];
						if(nk >= 4)
						{
							nBadScore = 1;
						}
						else
						{
							dMotifScore += DMGETAT(pLogPWM, (nMotifLen-1-ni), (3-nk));
						}
					}
				}
				else
				{
					for(ni=0; ni<nMotifLen; ni++)
					{
						dCSTot += (double)pCSVec[ni];
						nRepeatTot += (int)pRepeatVec[ni];
						nk = (int)pBaseVec[ni];
						if(nk >= 4)
						{
							nBadScore = 1;
						}
						else
						{
							dMotifScore += DMGETAT(pLogPWM, ni, nk);
						}
					}
				}
			}
			else
			{
				for(ni=0; ni<nMotifLen; ni++)
				{
					dCSTot += (double)pCSVec[ni];
					nRepeatTot += (int)pRepeatVec[ni];
				}
			}
			
			/* conservation filter */
			dCSTot /= nMotifLen;
			if(dCSTot < dCSCutoff)
			{
				nFiltered = 1;
			}

			/* repeat filter */
			if(nRepeatMask == 1)
			{
				dRepeatRatio = (double)nRepeatTot/(double)nMotifLen;
				if(dRepeatRatio > 0.5)
				{
					nFiltered = 1;
				}
			}

			if(nTypeIIMask == 1)
			{
				if( (nBadScore == 1) || (dMotifScore < pTypeIICutoff->pMatElement[nMotifId]) )
				{
					nFiltered = 1;
				}
			}

			/* write */
			if(nFiltered == 0)
			{
				fprintf(fpFinalMap, "%s\t%d\t%d\t%s\t%f\t%c\n", 
					strChrName, nStart, nEnd, strAlias, dScore, chStrand);
			}

		}

		/* close files */
		fclose(fpInitMap);
		fclose(fpFinalMap);

		DestroyDoubleMatrix(pLogPWM);
		nMotifId++;
	}

	/* close file */
	fclose(fpMotifList);

	/* release memory */
	SEQMTFDESTROY(pSeqMtf);
	DestroyDoubleMatrix(pTypeIICutoff);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGeneCreate():                                                       */
/*  Create a new object of tagRefGene.                                     */
/* ----------------------------------------------------------------------- */ 
struct tagRefGene * RefGeneCreate()
{
	struct tagRefGene *pRefGene;

	pRefGene = NULL;
	pRefGene = (struct tagRefGene *)calloc(1, sizeof(struct tagRefGene)); 
	if(pRefGene == NULL)
	{
		printf("Error: Can't create tagRefGene structure.\n");
	}
	strcpy(pRefGene->strGene, "");
	strcpy(pRefGene->strName, "");
	strcpy(pRefGene->strChrom, "");
	pRefGene->nGeneID = -1;

	return pRefGene;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGeneDestroy():                                                      */
/*  Destroy a new object of tagRefGene.                                    */
/* ----------------------------------------------------------------------- */ 
int RefGeneDestroy(struct tagRefGene *pRefGene)
{
	if(pRefGene == NULL)
		return PROC_SUCCESS;

	if(pRefGene->pmatExonStartsEnds != NULL)
	{
		DestroyIntMatrix(pRefGene->pmatExonStartsEnds);
		pRefGene->pmatExonStartsEnds = NULL;
		pRefGene->pNext = NULL;
	}

	free(pRefGene);

	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGeneClone():                                                        */
/*  Clone an object of tagRefGene.                                         */
/* ----------------------------------------------------------------------- */ 
struct tagRefGene *RefGeneClone(struct tagRefGene *pRefGene)
{
	/* define */
	struct tagRefGene *pNewRefGene;

	/* init */
	if(pRefGene == NULL)
		return NULL;

	/* clone */
	pNewRefGene = NULL;
	pNewRefGene = RefGeneCreate();
	if(pNewRefGene == NULL)
	{
		printf("Error: cannot allocate memory for cloning refgene!\n");
		exit(EXIT_FAILURE);
	}

	pNewRefGene->chStrand = pRefGene->chStrand;
	pNewRefGene->nCdsEnd = pRefGene->nCdsEnd;
	pNewRefGene->nCdsStart = pRefGene->nCdsStart;
	pNewRefGene->nChrom = pRefGene->nChrom;
	pNewRefGene->nExonCount = pRefGene->nExonCount;
	pNewRefGene->nTxEnd = pRefGene->nTxEnd;
	pNewRefGene->nTxStart = pRefGene->nTxStart;
	pNewRefGene->pmatExonStartsEnds = CreateIntMatrix(pRefGene->pmatExonStartsEnds->nHeight,pRefGene->pmatExonStartsEnds->nWidth);
	if(pNewRefGene->pmatExonStartsEnds == NULL)
	{
		printf("Error: cannot allocate memory for cloning refgene!\n");
		exit(EXIT_FAILURE);
	}
	IMCOPY(pNewRefGene->pmatExonStartsEnds, pRefGene->pmatExonStartsEnds);
	pNewRefGene->pNext = NULL;
	strcpy(pNewRefGene->strChrom, pRefGene->strChrom);
	strcpy(pNewRefGene->strName, pRefGene->strName);
	strcpy(pNewRefGene->strGene, pRefGene->strGene);
	pNewRefGene->nGeneID = pRefGene->nGeneID;

	/* return */
	return pNewRefGene;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGeneInit():                                                         */
/*  Initialize reference gene.                                             */
/* ----------------------------------------------------------------------- */ 
int RefGeneInit(struct tagRefGene * pRefGene, char strParamLine[], char strSpecies[])
{
	/* "Name of gene" */
    char  strName[NAME_LENGTH];  
	/* "Chromosome name" */
    char  strChrom[NAME_LENGTH];
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
	/* "Name of gene" */
    char  strExonStarts[LONG_LINE_LENGTH];  
	/* "Chromosome name" */
    char  strExonEnds[LONG_LINE_LENGTH];
	/* for loading positions */
	char  *pSep1,*pSep2;
	int   ni;
	int	  nPos;

	if(pRefGene == NULL)
		return PROC_FAILURE;

	StrTrimLeft(strParamLine);
	StrTrimRight(strParamLine);

	sscanf(strParamLine, "%s %s %c %d %d %d %d %d %s %s",
		strName, strChrom, &chStrand, &nTxStart, &nTxEnd,
		&nCdsStart, &nCdsEnd, &nExonCount,
		strExonStarts, strExonEnds);

	if(nExonCount > 500)
	{
		printf("Warning: %s has more than 500 exons, loading might be wrong.\n", 
			strName);
	}

	pRefGene->pNext = NULL;
	strcpy(pRefGene->strName, strName);
	strcpy(pRefGene->strChrom, strChrom);
	pRefGene->nChrom = Genome_ChromosomeName_To_Index(strChrom, strSpecies);
	pRefGene->chStrand = chStrand;
	pRefGene->nTxStart = nTxStart;
	pRefGene->nTxEnd = nTxEnd-1;
	pRefGene->nCdsStart = nCdsStart;
	pRefGene->nCdsEnd = nCdsEnd-1;
	pRefGene->nExonCount = nExonCount;
	pRefGene->pmatExonStartsEnds = NULL;
	if(nExonCount > 0)
	{
		pRefGene->pmatExonStartsEnds = CreateIntMatrix(nExonCount,2);
		if(pRefGene->pmatExonStartsEnds == NULL)
		{
			printf("Error: Can't create exon start-end matrix for %s.\n", strName);
			return PROC_FAILURE;
		}
		/* load starts */
		ni = 0;
		pSep1 = strExonStarts;
		pSep2 = strchr(pSep1, ',');
		while(pSep2 != NULL)
		{
			*pSep2 = '\0';
			nPos = atoi(pSep1);
			IMSETAT(pRefGene->pmatExonStartsEnds, ni, 0, nPos);
			pSep1 = pSep2+1;
			pSep2 = strchr(pSep1, ',');
			ni++;
		}
		if(ni != nExonCount)
		{
			printf("Error: %s exon start-end loading error.\n", strName);
			return PROC_FAILURE;
		}

		/* load ends */
		ni = 0;
		pSep1 = strExonEnds;
		pSep2 = strchr(pSep1, ',');
		while(pSep2 != NULL)
		{
			*pSep2 = '\0';
			nPos = atoi(pSep1)-1;
			IMSETAT(pRefGene->pmatExonStartsEnds, ni, 1, nPos);
			pSep1 = pSep2+1;
			pSep2 = strchr(pSep1, ',');
			ni++;
		}
		if(ni != nExonCount)
		{
			printf("Error: %s exon start-end loading error.\n", strName);
			return PROC_FAILURE;
		}
	}

	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefFlatInit():                                                         */
/*  Initialize reference gene.                                             */
/* ----------------------------------------------------------------------- */ 
int RefFlatInit(struct tagRefGene * pRefGene, char strParamLine[], char strSpecies[])
{
	/* "Gene Name" */
	char strGene[GENE_NAME_LENGTH];
	/* "RefSeq Id" */
    char  strName[NAME_LENGTH];  
	/* "Chromosome name" */
    char  strChrom[NAME_LENGTH];
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
	/* "Name of gene" */
    char  strExonStarts[LONG_LINE_LENGTH];  
	/* "Chromosome name" */
    char  strExonEnds[LONG_LINE_LENGTH];
	/* for loading positions */
	char  *pSep1,*pSep2;
	int   ni;
	int	  nPos;

	if(pRefGene == NULL)
		return PROC_FAILURE;

	/* StrTrimLeft(strParamLine); */
	StrTrimRight(strParamLine);

	pSep1 = strchr(strParamLine, '\t');
	*pSep1 = '\0';
	strcpy(strGene, strParamLine);
	StrTrimSpace(strGene);
	if(strcmp(strGene, "") == 0)
	{
		strcpy(strGene, "---");
	}

	pSep1 = pSep1+1;
	sscanf(pSep1, "%s %s %c %d %d %d %d %d %s %s",
		strName, strChrom, &chStrand, &nTxStart, &nTxEnd,
		&nCdsStart, &nCdsEnd, &nExonCount,
		strExonStarts, strExonEnds);

	if(nExonCount > 500)
	{
		printf("Warning: %s has more than 500 exons, loading might be wrong.\n", 
			strName);
	}

	pRefGene->pNext = NULL;
	strcpy(pRefGene->strGene, strGene);
	strcpy(pRefGene->strName, strName);
	strcpy(pRefGene->strChrom, strChrom);
	pRefGene->nChrom = Genome_ChromosomeName_To_Index(strChrom, strSpecies);
	pRefGene->chStrand = chStrand;
	pRefGene->nTxStart = nTxStart;
	pRefGene->nTxEnd = nTxEnd-1;
	pRefGene->nCdsStart = nCdsStart;
	pRefGene->nCdsEnd = nCdsEnd-1;
	pRefGene->nExonCount = nExonCount;
	pRefGene->pmatExonStartsEnds = NULL;
	if(nExonCount > 0)
	{
		pRefGene->pmatExonStartsEnds = CreateIntMatrix(nExonCount,2);
		if(pRefGene->pmatExonStartsEnds == NULL)
		{
			printf("Error: Can't create exon start-end matrix for %s.\n", strName);
			return PROC_FAILURE;
		}
		/* load starts */
		ni = 0;
		pSep1 = strExonStarts;
		pSep2 = strchr(pSep1, ',');
		while(pSep2 != NULL)
		{
			*pSep2 = '\0';
			nPos = atoi(pSep1);
			IMSETAT(pRefGene->pmatExonStartsEnds, ni, 0, nPos);
			pSep1 = pSep2+1;
			pSep2 = strchr(pSep1, ',');
			ni++;
		}
		if(ni != nExonCount)
		{
			printf("Error: %s exon start-end loading error.\n", strName);
			return PROC_FAILURE;
		}

		/* load ends */
		ni = 0;
		pSep1 = strExonEnds;
		pSep2 = strchr(pSep1, ',');
		while(pSep2 != NULL)
		{
			*pSep2 = '\0';
			nPos = atoi(pSep1)-1;
			IMSETAT(pRefGene->pmatExonStartsEnds, ni, 1, nPos);
			pSep1 = pSep2+1;
			pSep2 = strchr(pSep1, ',');
			ni++;
		}
		if(ni != nExonCount)
		{
			printf("Error: %s exon start-end loading error.\n", strName);
			return PROC_FAILURE;
		}
	}

	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGeneInit_FromGenomeLabFormat():                                     */
/*  Initialize reference gene from a genome lab refgene format.            */
/* ----------------------------------------------------------------------- */ 
int RefGeneInit_FromGenomeLabFormat(struct tagRefGene * pRefGene, char strParamLine[], char strSpecies[])
{
	/* "Name of gene" */
    char  strName[NAME_LENGTH];  
	/* "Chromosome name" */
    char  strChrom[NAME_LENGTH];
	int nChrom;
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
	/* "Name of gene" */
    char  strExonStarts[LONG_LINE_LENGTH];  
	/* "Chromosome name" */
    char  strExonEnds[LONG_LINE_LENGTH];
	/* for loading positions */
	char  *pSep1,*pSep2;
	int   ni;
	int	  nPos;

	if(pRefGene == NULL)
		return PROC_FAILURE;

	StrTrimLeft(strParamLine);
	StrTrimRight(strParamLine);

	sscanf(strParamLine, "%s %d %c %d %d %d %d %d %s %s",
		strName, &nChrom, &chStrand, &nTxStart, &nTxEnd,
		&nCdsStart, &nCdsEnd, &nExonCount,
		strExonStarts, strExonEnds);

	if(nExonCount > 500)
	{
		printf("Warning: %s has more than 500 exons, loading might be wrong.\n", 
			strName);
	}

	pRefGene->pNext = NULL;
	strcpy(pRefGene->strName, strName);
	pRefGene->nChrom = nChrom;
	Genome_Index_To_ChromosomeName(strChrom, strSpecies, nChrom);
	strcpy(pRefGene->strChrom, strChrom);
	pRefGene->chStrand = chStrand;
	pRefGene->nTxStart = nTxStart;
	pRefGene->nTxEnd = nTxEnd;
	pRefGene->nCdsStart = nCdsStart;
	pRefGene->nCdsEnd = nCdsEnd;
	pRefGene->nExonCount = nExonCount;
	pRefGene->pmatExonStartsEnds = NULL;
	if(nExonCount > 0)
	{
		pRefGene->pmatExonStartsEnds = CreateIntMatrix(nExonCount,2);
		if(pRefGene->pmatExonStartsEnds == NULL)
		{
			printf("Error: Can't create exon start-end matrix for %s.\n", strName);
			return PROC_FAILURE;
		}
		/* load starts */
		ni = 0;
		pSep1 = strExonStarts;
		pSep2 = strchr(pSep1, ',');
		while(pSep2 != NULL)
		{
			*pSep2 = '\0';
			nPos = atoi(pSep1);
			IMSETAT(pRefGene->pmatExonStartsEnds, ni, 0, nPos);
			pSep1 = pSep2+1;
			pSep2 = strchr(pSep1, ',');
			ni++;
		}
		if(ni != nExonCount)
		{
			printf("Error: %s exon start-end loading error.\n", strName);
			return PROC_FAILURE;
		}

		/* load ends */
		ni = 0;
		pSep1 = strExonEnds;
		pSep2 = strchr(pSep1, ',');
		while(pSep2 != NULL)
		{
			*pSep2 = '\0';
			nPos = atoi(pSep1);
			IMSETAT(pRefGene->pmatExonStartsEnds, ni, 1, nPos);
			pSep1 = pSep2+1;
			pSep2 = strchr(pSep1, ',');
			ni++;
		}
		if(ni != nExonCount)
		{
			printf("Error: %s exon start-end loading error.\n", strName);
			return PROC_FAILURE;
		}
	}

	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefFlatInit_FromGenomeLabFormat():                                     */
/*  Initialize reference gene from a genome lab refgene format.            */
/* ----------------------------------------------------------------------- */ 
int RefFlatInit_FromGenomeLabFormat(struct tagRefGene * pRefGene, char strParamLine[], char strSpecies[])
{
	/* "Name of gene" */
    char  strGene[GENE_NAME_LENGTH];
	/* "RefSeq Id" */
	char  strName[NAME_LENGTH];  
	/* "Chromosome name" */
    char  strChrom[NAME_LENGTH];
	int nChrom;
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
	/* "Name of gene" */
    char  strExonStarts[LONG_LINE_LENGTH];  
	/* "Chromosome name" */
    char  strExonEnds[LONG_LINE_LENGTH];
	/* for loading positions */
	char  *pSep1,*pSep2;
	int   ni;
	int	  nPos;

	if(pRefGene == NULL)
		return PROC_FAILURE;

	StrTrimLeft(strParamLine);
	StrTrimRight(strParamLine);

	sscanf(strParamLine, "%s %s %d %c %d %d %d %d %d %s %s",
		strGene, strName, &nChrom, &chStrand, &nTxStart, &nTxEnd,
		&nCdsStart, &nCdsEnd, &nExonCount,
		strExonStarts, strExonEnds);

	if(nExonCount > 500)
	{
		printf("Warning: %s has more than 500 exons, loading might be wrong.\n", 
			strName);
	}

	pRefGene->pNext = NULL;
	strcpy(pRefGene->strGene, strGene);
	strcpy(pRefGene->strName, strName);
	pRefGene->nChrom = nChrom;
	Genome_Index_To_ChromosomeName(strChrom, strSpecies, nChrom);
	strcpy(pRefGene->strChrom, strChrom);
	pRefGene->chStrand = chStrand;
	pRefGene->nTxStart = nTxStart;
	pRefGene->nTxEnd = nTxEnd;
	pRefGene->nCdsStart = nCdsStart;
	pRefGene->nCdsEnd = nCdsEnd;
	pRefGene->nExonCount = nExonCount;
	pRefGene->pmatExonStartsEnds = NULL;
	if(nExonCount > 0)
	{
		pRefGene->pmatExonStartsEnds = CreateIntMatrix(nExonCount,2);
		if(pRefGene->pmatExonStartsEnds == NULL)
		{
			printf("Error: Can't create exon start-end matrix for %s.\n", strName);
			return PROC_FAILURE;
		}
		/* load starts */
		ni = 0;
		pSep1 = strExonStarts;
		pSep2 = strchr(pSep1, ',');
		while(pSep2 != NULL)
		{
			*pSep2 = '\0';
			nPos = atoi(pSep1);
			IMSETAT(pRefGene->pmatExonStartsEnds, ni, 0, nPos);
			pSep1 = pSep2+1;
			pSep2 = strchr(pSep1, ',');
			ni++;
		}
		if(ni != nExonCount)
		{
			printf("Error: %s exon start-end loading error.\n", strName);
			return PROC_FAILURE;
		}

		/* load ends */
		ni = 0;
		pSep1 = strExonEnds;
		pSep2 = strchr(pSep1, ',');
		while(pSep2 != NULL)
		{
			*pSep2 = '\0';
			nPos = atoi(pSep1);
			IMSETAT(pRefGene->pmatExonStartsEnds, ni, 1, nPos);
			pSep1 = pSep2+1;
			pSep2 = strchr(pSep1, ',');
			ni++;
		}
		if(ni != nExonCount)
		{
			printf("Error: %s exon start-end loading error.\n", strName);
			return PROC_FAILURE;
		}
	}

	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefLocusInit_FromGenomeLabFormat():                                    */
/*  Initialize reference gene from a genome lab reflocus format.           */
/* ----------------------------------------------------------------------- */ 
int RefLocusInit_FromGenomeLabFormat(struct tagRefGene * pRefGene, char strParamLine[], char strSpecies[])
{
	/* "Name of gene" */
    char  strGene[GENE_NAME_LENGTH];
	/* "RefSeq Id" */
	char  strName[NAME_LENGTH];
	/* "GeneID" */
	int nGeneID;
	/* "Chromosome name" */
    char  strChrom[NAME_LENGTH];
	int nChrom;
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
	/* "Name of gene" */
    char  strExonStarts[LONG_LINE_LENGTH];  
	/* "Chromosome name" */
    char  strExonEnds[LONG_LINE_LENGTH];
	/* for loading positions */
	char  *pSep1,*pSep2;
	int   ni;
	int	  nPos;

	if(pRefGene == NULL)
		return PROC_FAILURE;

	StrTrimLeft(strParamLine);
	StrTrimRight(strParamLine);

	sscanf(strParamLine, "%d %s %s %d %c %d %d %d %d %d %s %s",
		&nGeneID, strGene, strName, &nChrom, &chStrand, &nTxStart, &nTxEnd,
		&nCdsStart, &nCdsEnd, &nExonCount,
		strExonStarts, strExonEnds);

	if(nExonCount > 500)
	{
		printf("Warning: %s has more than 500 exons, loading might be wrong.\n", 
			strName);
	}

	pRefGene->pNext = NULL;
	pRefGene->nGeneID = nGeneID;
	strcpy(pRefGene->strGene, strGene);
	strcpy(pRefGene->strName, strName);
	pRefGene->nChrom = nChrom;
	Genome_Index_To_ChromosomeName(strChrom, strSpecies, nChrom);
	strcpy(pRefGene->strChrom, strChrom);
	pRefGene->chStrand = chStrand;
	pRefGene->nTxStart = nTxStart;
	pRefGene->nTxEnd = nTxEnd;
	pRefGene->nCdsStart = nCdsStart;
	pRefGene->nCdsEnd = nCdsEnd;
	pRefGene->nExonCount = nExonCount;
	pRefGene->pmatExonStartsEnds = NULL;
	if(nExonCount > 0)
	{
		pRefGene->pmatExonStartsEnds = CreateIntMatrix(nExonCount,2);
		if(pRefGene->pmatExonStartsEnds == NULL)
		{
			printf("Error: Can't create exon start-end matrix for %s.\n", strName);
			return PROC_FAILURE;
		}
		/* load starts */
		ni = 0;
		pSep1 = strExonStarts;
		pSep2 = strchr(pSep1, ',');
		while(pSep2 != NULL)
		{
			*pSep2 = '\0';
			nPos = atoi(pSep1);
			IMSETAT(pRefGene->pmatExonStartsEnds, ni, 0, nPos);
			pSep1 = pSep2+1;
			pSep2 = strchr(pSep1, ',');
			ni++;
		}
		if(ni != nExonCount)
		{
			printf("Error: %s exon start-end loading error.\n", strName);
			return PROC_FAILURE;
		}

		/* load ends */
		ni = 0;
		pSep1 = strExonEnds;
		pSep2 = strchr(pSep1, ',');
		while(pSep2 != NULL)
		{
			*pSep2 = '\0';
			nPos = atoi(pSep1);
			IMSETAT(pRefGene->pmatExonStartsEnds, ni, 1, nPos);
			pSep1 = pSep2+1;
			pSep2 = strchr(pSep1, ',');
			ni++;
		}
		if(ni != nExonCount)
		{
			printf("Error: %s exon start-end loading error.\n", strName);
			return PROC_FAILURE;
		}
	}

	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGeneWrite():                                                        */
/*  Write refgene to file.                                                 */
/* ----------------------------------------------------------------------- */ 
int RefGeneWrite(struct tagRefGene * pRefGene, FILE *fpOut)
{
	int ni;

	if(pRefGene == NULL)
		return PROC_FAILURE;
	if(fpOut == NULL)
		return PROC_FAILURE;

	fprintf(fpOut, "%s\t%d\t%c\t%d\t%d\t%d\t%d\t%d\t",
		pRefGene->strName, pRefGene->nChrom, pRefGene->chStrand, 
		pRefGene->nTxStart, pRefGene->nTxEnd,
		pRefGene->nCdsStart, pRefGene->nCdsEnd, pRefGene->nExonCount);
	
	for(ni=0; ni<pRefGene->nExonCount; ni++)
	{
		fprintf(fpOut, "%d,", IMGETAT(pRefGene->pmatExonStartsEnds, ni, 0));
	}
	fprintf(fpOut, "\t");
	for(ni=0; ni<pRefGene->nExonCount; ni++)
	{
		fprintf(fpOut, "%d,", IMGETAT(pRefGene->pmatExonStartsEnds, ni, 1));
	}

	fprintf(fpOut, "\n");
	

	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefFlatWrite():                                                        */
/*  Write refgene to file.                                                 */
/* ----------------------------------------------------------------------- */ 
int RefFlatWrite(struct tagRefGene *pRefGene, FILE *fpOut)
{
	int ni;

	if(pRefGene == NULL)
		return PROC_FAILURE;
	if(fpOut == NULL)
		return PROC_FAILURE;

	fprintf(fpOut, "%s\t%s\t%d\t%c\t%d\t%d\t%d\t%d\t%d\t",
		pRefGene->strGene, pRefGene->strName, pRefGene->nChrom, pRefGene->chStrand, 
		pRefGene->nTxStart, pRefGene->nTxEnd,
		pRefGene->nCdsStart, pRefGene->nCdsEnd, pRefGene->nExonCount);
	
	for(ni=0; ni<pRefGene->nExonCount; ni++)
	{
		fprintf(fpOut, "%d,", IMGETAT(pRefGene->pmatExonStartsEnds, ni, 0));
	}
	fprintf(fpOut, "\t");
	for(ni=0; ni<pRefGene->nExonCount; ni++)
	{
		fprintf(fpOut, "%d,", IMGETAT(pRefGene->pmatExonStartsEnds, ni, 1));
	}

	fprintf(fpOut, "\n");
	

	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefLocusWrite():                                                       */
/*  Write refgene to file.                                                 */
/* ----------------------------------------------------------------------- */ 
int RefLocusWrite(struct tagRefGene *pRefGene, FILE *fpOut)
{
	int ni;

	if(pRefGene == NULL)
		return PROC_FAILURE;
	if(fpOut == NULL)
		return PROC_FAILURE;

	fprintf(fpOut, "%d\t%s\t%s\t%d\t%c\t%d\t%d\t%d\t%d\t%d\t", pRefGene->nGeneID,
		pRefGene->strGene, pRefGene->strName, pRefGene->nChrom, pRefGene->chStrand, 
		pRefGene->nTxStart, pRefGene->nTxEnd,
		pRefGene->nCdsStart, pRefGene->nCdsEnd, pRefGene->nExonCount);
	
	for(ni=0; ni<pRefGene->nExonCount; ni++)
	{
		fprintf(fpOut, "%d,", IMGETAT(pRefGene->pmatExonStartsEnds, ni, 0));
	}
	fprintf(fpOut, "\t");
	for(ni=0; ni<pRefGene->nExonCount; ni++)
	{
		fprintf(fpOut, "%d,", IMGETAT(pRefGene->pmatExonStartsEnds, ni, 1));
	}

	fprintf(fpOut, "\n");
	

	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGeneInsert():                                                       */
/*  Insert refgene.                                                        */
/* ----------------------------------------------------------------------- */ 
int RefGeneInsert(struct tagRefGene **vList, struct tagRefGene *pRefGene)
{
	/* define */
	struct tagRefGene *pPrev,*pNext;
	int nOK;

	/* insert */
	if(pRefGene == NULL)
		return PROC_SUCCESS;

	/* if list is empty */
	if(*vList == NULL)
	{
		*vList = pRefGene;
	}
	/* if list is not empty */
	else
	{
		pPrev = NULL;
		pNext = *vList;
		nOK = 0;
		while(pNext != NULL)
		{
			if(pRefGene->nChrom < pNext->nChrom)
			{
				nOK = 1;
				break;
			}
			else if(pRefGene->nChrom == pNext->nChrom)
			{
				if(pRefGene->nTxStart < pNext->nTxStart)
				{
					nOK = 1;
					break;
				}
				else if(pRefGene->nTxStart == pNext->nTxStart)
				{
					if(pRefGene->nTxEnd < pNext->nTxEnd)
					{
						nOK = 1;
						break;
					}
					else if(pRefGene->nTxEnd == pNext->nTxEnd)
					{
						if(pRefGene->nCdsStart < pNext->nCdsStart)
						{
							nOK = 1;
							break;
						}
						else if(pRefGene->nCdsStart == pNext->nCdsStart)
						{
							if(pRefGene->nCdsEnd < pNext->nCdsEnd)
							{
								nOK = 1;
								break;
							}
						}
					}
				}
				
			}

			/* get next */
			pPrev = pNext;
			pNext = pNext->pNext;
		}

		if(pPrev == NULL)
		{
			pRefGene->pNext = pNext;
			*vList = pRefGene;
		}
		else
		{
			pRefGene->pNext = pNext;
			pPrev->pNext = pRefGene;
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  GenomicRegionCreate():                                                 */
/*  Create a new object of tagGenomicRegion.                               */
/* ----------------------------------------------------------------------- */ 
struct tagGenomicRegion * GenomicRegionCreate()
{
	struct tagGenomicRegion *pRegion;

	pRegion = NULL;
	pRegion = (struct tagGenomicRegion *)calloc(1, sizeof(struct tagGenomicRegion)); 
	if(pRegion == NULL)
	{
		printf("Error: Can't create tagGenomicRegion structure.\n");
	}

	return pRegion;
}

/* ----------------------------------------------------------------------- */ 
/*  GenomicRegionDestroy():                                                */
/*  Destroy a new object of tagGenomicRegion.                              */
/* ----------------------------------------------------------------------- */ 
int GenomicRegionDestroy(struct tagGenomicRegion *pRegion)
{
	if(pRegion == NULL)
		return PROC_SUCCESS;

	pRegion->pNext = NULL;
	free(pRegion);

	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  GenomicRegionWrite():                                                  */
/*  Write a genomic region to file.                                        */
/*  If nChromAsNumber == 1, write chromosome as a numeric index.           */
/*  Else write chromosome as a string index.                               */
/* ----------------------------------------------------------------------- */ 
int GenomicRegionWrite(struct tagGenomicRegion *pRegion, FILE *fpOut, int nChromAsNumber)
{
	if(pRegion == NULL)
	{
		printf("Warning: GenomicRegionWrite, null genomic region!\n");
		return PROC_FAILURE;
	}
	if(fpOut == NULL)
	{
		printf("Warning: GenomicRegionWrite, null output file!\n");
		return PROC_FAILURE;
	}

	if(nChromAsNumber == 1)
	{
		fprintf(fpOut, "%d\t", pRegion->nChrom);
	}
	else
	{
		fprintf(fpOut, "%s\t", pRegion->strChrom);
	}

	fprintf(fpOut, "%d\t%d\t%c\n", pRegion->nStart, pRegion->nEnd, pRegion->chStrand);
	
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  GenomicRegionInsert():                                                 */
/*  Insert a genomic region to a ordered region list.                      */
/* ----------------------------------------------------------------------- */ 
int GenomicRegionInsert(struct tagGenomicRegion **vList, struct tagGenomicRegion *pRegion)
{
	/* define */
	struct tagGenomicRegion *pPrev,*pNext;
	int nOK;

	/* insert */
	if(pRegion == NULL)
		return PROC_SUCCESS;

	/* if list is empty */
	if(*vList == NULL)
	{
		*vList = pRegion;
	}
	/* if list is not empty */
	else
	{
		pPrev = NULL;
		pNext = *vList;
		nOK = 0;
		while(pNext != NULL)
		{
			if(pRegion->nChrom < pNext->nChrom)
			{
				nOK = 1;
				break;
			}
			else if(pRegion->nChrom == pNext->nChrom)
			{
				if(pRegion->nStart < pNext->nStart)
				{
					nOK = 1;
					break;
				}
				else if(pRegion->nStart == pNext->nStart)
				{
					if(pRegion->nEnd < pNext->nEnd)
					{
						nOK = 1;
						break;
					}
					else if(pRegion->nEnd == pNext->nEnd)
					{
					}
				}
				
			}

			/* get next */
			pPrev = pNext;
			pNext = pNext->pNext;
		}

		if(pPrev == NULL)
		{
			pRegion->pNext = pNext;
			*vList = pRegion;
		}
		else
		{
			pRegion->pNext = pNext;
			pPrev->pNext = pRegion;
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  GenomicRegionGetUnique():                                              */
/*  get a uniquee genomic region list, and discard any redundant part.     */
/* ----------------------------------------------------------------------- */ 
int GenomicRegionGetUnique(struct tagGenomicRegion **vList)
{
	/* define */
	struct tagGenomicRegion *pPrev,*pNext;
	
	/* if list is empty */
	if(*vList == NULL)
	{
		return PROC_SUCCESS;
	}
	/* if list is not empty */
	else
	{
		pPrev = *vList;
		pNext = pPrev->pNext;
		while(pNext != NULL)
		{
			if(pNext->nChrom != pPrev->nChrom)
			{
				/* get next */
				pPrev = pNext;
				pNext = pPrev->pNext;
			}
			else
			{
				/* if overlap */
				if(pNext->nStart <= pPrev->nEnd)
				{
					if(pNext->nEnd > pPrev->nEnd)
					{
						pPrev->nEnd = pNext->nEnd;
					}

					pPrev->pNext = pNext->pNext;
					GenomicRegionDestroy(pNext);
					pNext = pPrev->pNext;
				}
				else
				{
					/* get next */
					pPrev = pNext;
					pNext = pPrev->pNext;
				}
			}
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGene_Convert_UCSC_To_Lab():                                         */
/*  Convert ucsc refgene to lab format.                                    */
/* ----------------------------------------------------------------------- */ 
int RefGene_Convert_UCSC_To_Lab(char strInPath[], char strOutPath[], char strSpecies[], int nChrNum)
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	char strLine[LONG_LINE_LENGTH];
	struct tagRefGene **vAllGene;
	struct tagRefGene *pRefGene;
	int ni;

	/* init */
	vAllGene = NULL;
	vAllGene = (struct tagRefGene **)calloc(nChrNum, sizeof(struct tagRefGene *));
	if(vAllGene == NULL)
	{
		printf("Error: can't create enough space for refgene annotation.\n");
		exit(EXIT_FAILURE);
	}

	fpIn = NULL;
	fpIn = fopen(strInPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: can't open input file.\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: can't open output file.\n");
		exit(EXIT_FAILURE);
	}

	/* load and sort */
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		pRefGene = NULL;
		pRefGene = RefGeneCreate();
		if(pRefGene != NULL)
		{
			if(RefGeneInit(pRefGene, strLine, strSpecies) == PROC_SUCCESS)
			{
				if( (pRefGene->nChrom > 0) && (pRefGene->nChrom <= nChrNum) )
				{
					ni = pRefGene->nChrom-1;
					RefGeneInsert((vAllGene+ni), pRefGene);
				}
				else
				{
					RefGeneDestroy(pRefGene);
				}
			}
			else
			{
				RefGeneDestroy(pRefGene);
			}
		}
	}

	/* write and destroy */
	for(ni=0; ni<nChrNum; ni++)
	{
		while(vAllGene[ni] != NULL)
		{
			pRefGene = vAllGene[ni];
			vAllGene[ni] = vAllGene[ni]->pNext;
			RefGeneWrite(pRefGene, fpOut);
			RefGeneDestroy(pRefGene);
		}
	}

	/* release memory */
	fclose(fpIn);
	fclose(fpOut);
	free(vAllGene);

	/* exit */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefFlat_Convert_UCSC_To_Lab():                                         */
/*  Convert Ucsc refgene to lab format.                                    */
/* ----------------------------------------------------------------------- */ 
int RefFlat_Convert_UCSC_To_Lab(char strInPath[], char strOutPath[], char strSpecies[], int nChrNum)
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	char strLine[LONG_LINE_LENGTH];
	struct tagRefGene **vAllGene;
	struct tagRefGene *pRefGene;
	int ni;

	/* init */
	vAllGene = NULL;
	vAllGene = (struct tagRefGene **)calloc(nChrNum, sizeof(struct tagRefGene *));
	if(vAllGene == NULL)
	{
		printf("Error: can't create enough space for refgene annotation.\n");
		exit(EXIT_FAILURE);
	}

	fpIn = NULL;
	fpIn = fopen(strInPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: can't open input file.\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: can't open output file.\n");
		exit(EXIT_FAILURE);
	}

	/* load and sort */
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		pRefGene = NULL;
		pRefGene = RefGeneCreate();
		if(pRefGene != NULL)
		{
			if(RefFlatInit(pRefGene, strLine, strSpecies) == PROC_SUCCESS)
			{
				if( (pRefGene->nChrom > 0) && (pRefGene->nChrom <= nChrNum) )
				{
					ni = pRefGene->nChrom-1;
					RefGeneInsert((vAllGene+ni), pRefGene);
				}
				else
				{
					RefGeneDestroy(pRefGene);
				}
			}
			else
			{
				RefGeneDestroy(pRefGene);
			}
		}
	}

	/* write and destroy */
	for(ni=0; ni<nChrNum; ni++)
	{
		while(vAllGene[ni] != NULL)
		{
			pRefGene = vAllGene[ni];
			vAllGene[ni] = vAllGene[ni]->pNext;
			RefFlatWrite(pRefGene, fpOut);
			RefGeneDestroy(pRefGene);
		}
	}

	/* release memory */
	fclose(fpIn);
	fclose(fpOut);
	free(vAllGene);

	/* exit */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGene_AnnotateWithLocusID():                                         */
/*  Annotate refgene with gene ID. The input file is a sorted file that    */
/*  contains the refflat annotations and refgene2geneid annotations.       */
/* ----------------------------------------------------------------------- */ 
int RefGene_AnnotateWithLocusID(char strInPath[], char strOutPath[], char strSpecies[], int nChrNum, int nChangeGeneName)
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	char strLine[LONG_LINE_LENGTH];
	struct tagRefGene *pRawGeneList = NULL;

	struct tagRefGene **vAllGene;
	struct tagRefGene *pRefGene;
	struct tagRefGene *pToDoRefGene = NULL;
	struct tagRefGene *pTempRefGene = NULL;
	struct tagRefGene *pLastRefGene = NULL;
	int ni;

	char strHead[LINE_LENGTH];
	char strGeneRefName[LINE_LENGTH];
	int nGeneID = -1;
	int nTaxID = -1;
	char strNewGeneName[LINE_LENGTH];

	/* init */
	fpIn = NULL;
	fpIn = fopen(strInPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: can't open input file.\n");
		exit(EXIT_FAILURE);
	}

	strcpy(strGeneRefName, "NA");
	strcpy(strNewGeneName, "-");
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%s", strHead);
		if(strcmp(strHead, "GENEID") == 0)
		{
			if(nChangeGeneName == 1)
			{
				sscanf(strLine, "%s %s %d %d %s", strHead, strGeneRefName, &nGeneID, &nTaxID, strNewGeneName);
			}
			else
			{
				sscanf(strLine, "%s %s %d", strHead, strGeneRefName, &nGeneID);
			}

			while(pToDoRefGene != NULL)
			{
				if(strcmp(pToDoRefGene->strName, strGeneRefName) == 0)
				{
					pToDoRefGene->nGeneID = nGeneID;
					if(nChangeGeneName == 1)
						strcpy(pToDoRefGene->strGene, strNewGeneName);
				}
				pToDoRefGene = pToDoRefGene->pNext;
			}
		}
		else
		{
			pRefGene = NULL;
			pRefGene = RefGeneCreate();
			if(pRefGene == NULL)
			{
				printf("Error: RefGene_AnnotateWithLocusID, cannot create refgene!\n");
				exit(EXIT_FAILURE);
			}
			RefFlatInit_FromGenomeLabFormat(pRefGene, strLine, strSpecies);
			if(pRawGeneList == NULL)
			{
				pRawGeneList = pRefGene;
				pLastRefGene = pRefGene;
			}
			else
			{
				pLastRefGene->pNext = pRefGene;
				pLastRefGene = pRefGene;
			}
			
			if(strcmp(pRefGene->strName, strGeneRefName) == 0)
			{
				pRefGene->nGeneID = nGeneID;
				if(nChangeGeneName == 1)
					strcpy(pRefGene->strGene, strNewGeneName);
				pToDoRefGene = NULL;
			}
			else
			{
				if(pToDoRefGene == NULL)
				{
					pToDoRefGene = pRefGene;
				}
			}
		}
	}

	fclose(fpIn);
	

	/* sort */
	vAllGene = NULL;
	vAllGene = (struct tagRefGene **)calloc(nChrNum, sizeof(struct tagRefGene *));
	if(vAllGene == NULL)
	{
		printf("Error: can't create enough space for refgene annotation.\n");
		exit(EXIT_FAILURE);
	}

	while(pRawGeneList != NULL)
	{
		pRefGene = pRawGeneList;
		pRawGeneList = pRefGene->pNext;
		pRefGene->pNext = NULL;

		if( (pRefGene->nChrom > 0) && (pRefGene->nChrom <= nChrNum) )
		{
			ni = pRefGene->nChrom-1;
			RefGeneInsert((vAllGene+ni), pRefGene);
		}
		else
		{
			RefGeneDestroy(pRefGene);
		}
	}

	/* write and destroy */
	fpOut = NULL;
	fpOut = fopen(strOutPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: can't open output file.\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nChrNum; ni++)
	{
		while(vAllGene[ni] != NULL)
		{
			pRefGene = vAllGene[ni];
			vAllGene[ni] = vAllGene[ni]->pNext;
			RefLocusWrite(pRefGene, fpOut);
			RefGeneDestroy(pRefGene);
		}
	}

	/* release memory */
	fclose(fpOut);
	free(vAllGene);

	/* exit */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGene_AnnotateWithExonArrayID():                                     */
/*  Annotate refgene with exon array ID. The input file is a sorted file   */
/*  that contains the refflat annotations and refgene2exonid annotations.  */
/* ----------------------------------------------------------------------- */ 
int RefGene_AnnotateWithExonArrayID(char strInPath[], char strOutPath[], char strSpecies[], int nChrNum)
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	char strLine[LONG_LINE_LENGTH];
	struct tagRefGene *pRawGeneList = NULL;

	struct tagRefGene **vAllGene;
	struct tagRefGene *pRefGene;
	struct tagRefGene *pToDoRefGene = NULL;
	struct tagRefGene *pTempRefGene = NULL;
	struct tagRefGene *pLastRefGene = NULL;
	int ni;
	int nScaleFactor = 1000000000;

	char strHead[LINE_LENGTH];
	char strGeneRefName[LINE_LENGTH];
	int nGeneID = -1;
	int nLineID = -1;

	/* init */
	fpIn = NULL;
	fpIn = fopen(strInPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: can't open input file.\n");
		exit(EXIT_FAILURE);
	}

	strcpy(strGeneRefName, "NA");
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%s", strHead);
		if(strcmp(strHead, "GENEID") == 0)
		{
			sscanf(strLine, "%s %s %d %d", strHead, strGeneRefName, &nLineID, &nGeneID);
			while(pToDoRefGene != NULL)
			{
				if(strcmp(pToDoRefGene->strName, strGeneRefName) == 0)
				{
					pToDoRefGene->nGeneID = nLineID;
				}
				pToDoRefGene = pToDoRefGene->pNext;
			}
		}
		else
		{
			pRefGene = NULL;
			pRefGene = RefGeneCreate();
			if(pRefGene == NULL)
			{
				printf("Error: RefGene_AnnotateWithLocusID, cannot create refgene!\n");
				exit(EXIT_FAILURE);
			}
			RefFlatInit_FromGenomeLabFormat(pRefGene, strLine, strSpecies);
			if(pRawGeneList == NULL)
			{
				pRawGeneList = pRefGene;
				pLastRefGene = pRefGene;
			}
			else
			{
				pLastRefGene->pNext = pRefGene;
				pLastRefGene = pRefGene;
			}
			
			if(strcmp(pRefGene->strName, strGeneRefName) == 0)
			{
				pRefGene->nGeneID = nLineID;
				pToDoRefGene = NULL;
			}
			else
			{
				if(pToDoRefGene == NULL)
				{
					pToDoRefGene = pRefGene;
				}
			}
		}
	}

	fclose(fpIn);
	

	/* sort */
	vAllGene = NULL;
	vAllGene = (struct tagRefGene **)calloc(nChrNum, sizeof(struct tagRefGene *));
	if(vAllGene == NULL)
	{
		printf("Error: can't create enough space for refgene annotation.\n");
		exit(EXIT_FAILURE);
	}

	while(pRawGeneList != NULL)
	{
		pRefGene = pRawGeneList;
		pRawGeneList = pRefGene->pNext;
		pRefGene->pNext = NULL;

		if( (pRefGene->nChrom > 0) && (pRefGene->nChrom <= nChrNum) )
		{
			ni = pRefGene->nChrom-1;
			RefGeneInsert((vAllGene+ni), pRefGene);
		}
		else
		{
			RefGeneDestroy(pRefGene);
		}
	}

	/* write and destroy */
	fpOut = NULL;
	fpOut = fopen(strOutPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: can't open output file.\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nChrNum; ni++)
	{
		while(vAllGene[ni] != NULL)
		{
			pRefGene = vAllGene[ni];
			vAllGene[ni] = vAllGene[ni]->pNext;
			if(pRefGene->nGeneID<0)
			{
				fprintf(fpOut, "-1\n");
			}
			else
			{
				fprintf(fpOut, "%d\n", pRefGene->nGeneID);
			}
			RefGeneDestroy(pRefGene);
		}
	}

	/* release memory */
	fclose(fpOut);
	free(vAllGene);

	/* exit */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGene_AnnotateWithMicroArrayID():                                    */
/*  Annotate refgene with exon array ID. The input file is a sorted file   */
/*  that contains the refflat annotations and refgene2exonid annotations.  */
/* ----------------------------------------------------------------------- */ 
int RefGene_AnnotateWithMicroArrayID(char strInPath[], char strOutPath[], char strSpecies[], int nChrNum)
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	char strLine[LONG_LINE_LENGTH];
	struct tagRefGene *pRawGeneList = NULL;

	struct tagRefGene **vAllGene;
	struct tagRefGene *pRefGene;
	struct tagRefGene *pToDoRefGene = NULL;
	struct tagRefGene *pTempRefGene = NULL;
	struct tagRefGene *pLastRefGene = NULL;
	int ni,nk;
	int nScaleFactor = 1000000000;

	char strHead[LINE_LENGTH];
	char strGeneRefName[LINE_LENGTH];
	char strAnnot[LINE_LENGTH];

	/* init */
	fpIn = NULL;
	fpIn = fopen(strInPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: can't open input file.\n");
		exit(EXIT_FAILURE);
	}

	strcpy(strGeneRefName, "NA");
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%s", strHead);
		if(strcmp(strHead, "GENEID") == 0)
		{
			sscanf(strLine, "%s %s %s", strHead, strGeneRefName, strAnnot);
			while(pToDoRefGene != NULL)
			{
				if(strcmp(pToDoRefGene->strName, strGeneRefName) == 0)
				{
					if(strlen(strAnnot) < NAME_LENGTH)
					{	
						strcpy(pToDoRefGene->strChrom, strAnnot);
						pToDoRefGene->nGeneID = 1;
					}
					else
					{
						pToDoRefGene->nGeneID = 0;
					}
				}
				else
				{
					pToDoRefGene->nGeneID = 0;
				}
				pToDoRefGene = pToDoRefGene->pNext;
			}
		}
		else
		{
			pRefGene = NULL;
			pRefGene = RefGeneCreate();
			if(pRefGene == NULL)
			{
				printf("Error: RefGene_AnnotateWithLocusID, cannot create refgene!\n");
				exit(EXIT_FAILURE);
			}
			
			RefFlatInit_FromGenomeLabFormat(pRefGene, strLine, strSpecies);
			if(pRawGeneList == NULL)
			{
				pRawGeneList = pRefGene;
				pLastRefGene = pRefGene;
			}
			else
			{
				pLastRefGene->pNext = pRefGene;
				pLastRefGene = pRefGene;
			}
			
			if(strcmp(pRefGene->strName, strGeneRefName) == 0)
			{
				if(strlen(strAnnot) < NAME_LENGTH)
				{	
					strcpy(pRefGene->strChrom, strAnnot);
					pRefGene->nGeneID = 1;
					pToDoRefGene = NULL;
				}
				else
				{
					pRefGene->nGeneID = 0;
					if(pToDoRefGene == NULL)
					{
						pToDoRefGene = pRefGene;
					}
				}
			}
			else
			{
				pRefGene->nGeneID = 0;
				if(pToDoRefGene == NULL)
				{
					pToDoRefGene = pRefGene;
				}
			}
		}
	}

	fclose(fpIn);
	

	/* sort */
	vAllGene = NULL;
	vAllGene = (struct tagRefGene **)calloc(nChrNum, sizeof(struct tagRefGene *));
	if(vAllGene == NULL)
	{
		printf("Error: can't create enough space for refgene annotation.\n");
		exit(EXIT_FAILURE);
	}

	nk = 0;
	while(pRawGeneList != NULL)
	{
		pRefGene = pRawGeneList;
		pRawGeneList = pRefGene->pNext;
		pRefGene->pNext = NULL;

		if( (pRefGene->nChrom > 0) && (pRefGene->nChrom <= nChrNum) )
		{
			ni = pRefGene->nChrom-1;
			RefGeneInsert((vAllGene+ni), pRefGene);
		}
		else
		{
			RefGeneDestroy(pRefGene);
		}

		nk++;
	}

	/* write and destroy */
	fpOut = NULL;
	fpOut = fopen(strOutPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: can't open output file.\n");
		exit(EXIT_FAILURE);
	}

	nk = 0;

	for(ni=0; ni<nChrNum; ni++)
	{
		while(vAllGene[ni] != NULL)
		{
			pRefGene = vAllGene[ni];
			vAllGene[ni] = vAllGene[ni]->pNext;
			if(pRefGene->nGeneID == 0)
			{
				fprintf(fpOut, "---\n");
			}
			else
			{
				fprintf(fpOut, "%s\n", pRefGene->strChrom);
			}
			RefGeneDestroy(pRefGene);
			nk++;
		}
	}

	/* release memory */
	fclose(fpOut);
	free(vAllGene);

	/* exit */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGene_GetNearestGene_Main                                            */
/*  get nearest gene from a refgene database.                              */
/* ----------------------------------------------------------------------- */ 
int RefGene_GetNearestGene_Main(char strDatabasePath[], int nDatabaseType,
			char strSpecies[], char strInputPath[], char strOutputPath[],
			int nRefType, int nUP, int nDOWN)
{
	/* out */
	FILE *fpIn;
	FILE *fpOut;
	int ni;
	struct tagRefGene **vSourceRefGene = NULL;
	int nSourceRefNum = 0;
	char strLine[LONG_LINE_LENGTH];
	char strSeqName[MED_LINE_LENGTH];
	char strChr[MED_LINE_LENGTH];
	int nChr,nStart,nEnd;
	char chStrand;
	struct tagRefGene *pRefGene;

	/* init */
	/* load source refgene */
	vSourceRefGene = NULL;
	vSourceRefGene = RefGene_LoadDatabase(strDatabasePath, nDatabaseType, 
		strSpecies, &nSourceRefNum);
	if(vSourceRefGene == NULL)
	{
		printf("Warning: RefGene_GetNearestGene_Main, null database!\n");
		return PROC_SUCCESS;
	}

	/* load coordinates */
	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: RefGene_GetNearestGene_Main, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: RefGene_GetNearestGene_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "#seq_id\tchromosome\tstart\tend\tstrand\tannotation\n");

	/* process one by one */
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		sscanf(strLine, "%s %s %d %d %c", strSeqName, strChr, &nStart, &nEnd, &chStrand);
		nChr = Genome_ChromosomeName_To_Index(strChr, strSpecies);
		ni = RefGene_GetNearestGene(nChr, nStart, nEnd, chStrand, 
			vSourceRefGene, nSourceRefNum, nRefType, nUP, nDOWN);
		if( (ni>=0) && (ni<nSourceRefNum) )
		{
			pRefGene = vSourceRefGene[ni];
		}
		else
		{
			pRefGene = NULL;
		}

		fprintf(fpOut, "%s\t%s\t%d\t%d\t%c\t", strSeqName, strChr, nStart, nEnd, chStrand);

		if(nDatabaseType == 2)
		{
			if(pRefGene == NULL)
			{
				fprintf(fpOut, "-1\t");
			}
			else
			{
				fprintf(fpOut, "%d\t", pRefGene->nGeneID);
			}
		}

		if(pRefGene == NULL)
		{
			fprintf(fpOut, "---\n"); 
		}
		else if(strcmp(pRefGene->strGene, "") == 0)
		{
			fprintf(fpOut, "%s\t", pRefGene->strName);
			RefGeneWrite(pRefGene, fpOut);
		}
		else
		{
			RefFlatWrite(pRefGene, fpOut);
		}
	}

	/* close files */
	fclose(fpIn);
	fclose(fpOut);
	
	/* clear database */
	RefGene_ClearDatabase(&vSourceRefGene, nSourceRefNum);

	/* exit */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGene_GetNearestGene                                                 */
/*  get nearest gene from a refgene database.                              */
/* ----------------------------------------------------------------------- */ 
int RefGene_GetNearestGene(int nChr, int nStart, int nEnd, char chStrand, 
			struct tagRefGene **vSourceRefGene, int nSourceRefNum, 
			int nRefType, int nUP, int nDOWN)
{
	/* define */
	int nOptId = -1;
	int nOptDist = IM_ACCESS_VIOLATION-1;
	int nOptLen;
	
	int ni,nMidPos,nTempLen;
	int nDistRef,nDistMin;

	/* init */
	nMidPos = (nStart+nEnd)/2;

	/* search */
	for(ni=0; ni<nSourceRefNum; ni++)
	{
		if(vSourceRefGene[ni] == NULL)
			continue;

		if(vSourceRefGene[ni]->nChrom != nChr)
			continue;

		if(nRefType == 0)
		{
			if( (nMidPos >= vSourceRefGene[ni]->nTxStart) &&
				(nMidPos <= vSourceRefGene[ni]->nTxEnd) )
			{
				nDistMin = 0;
			}
			else if(nMidPos <  vSourceRefGene[ni]->nTxStart)
			{
				nDistMin = vSourceRefGene[ni]->nTxStart-nMidPos;
				if(vSourceRefGene[ni]->chStrand == '-')
					nDistRef = nDOWN;
				else
					nDistRef = nUP;
				
				if(nDistMin > nDistRef)
					continue;
			}
			else
			{
				nDistMin = nMidPos-vSourceRefGene[ni]->nTxEnd;
				if(vSourceRefGene[ni]->chStrand == '-')
					nDistRef = nUP;
				else
					nDistRef = nDOWN;
				
				if(nDistMin > nDistRef)
					continue;
			}
		}
		else if(nRefType == 1)
		{
			if(vSourceRefGene[ni]->chStrand == '-')
			{
				if( (nMidPos <= (vSourceRefGene[ni]->nTxEnd+nUP)) &&
					(nMidPos >= (vSourceRefGene[ni]->nTxEnd-nDOWN)) )
				{
					nDistMin = (int)fabs(nMidPos-vSourceRefGene[ni]->nTxEnd);
				}
				else
				{
					continue;
				}
			}
			else
			{
				if( (nMidPos >= (vSourceRefGene[ni]->nTxStart-nUP)) &&
					(nMidPos <= (vSourceRefGene[ni]->nTxStart+nDOWN)) )
				{
					nDistMin = (int)fabs(nMidPos-vSourceRefGene[ni]->nTxStart);
				}
				else
				{
					continue;
				}
			}
		}
		else if(nRefType == 2)
		{
			if(vSourceRefGene[ni]->chStrand == '-')
			{
				if( (nMidPos <= (vSourceRefGene[ni]->nTxStart+nUP)) &&
					(nMidPos >= (vSourceRefGene[ni]->nTxStart-nDOWN)) )
				{
					nDistMin = (int)fabs(nMidPos-vSourceRefGene[ni]->nTxStart);
				}
				else
				{
					continue;
				}
			}
			else
			{
				if( (nMidPos >= (vSourceRefGene[ni]->nTxEnd-nUP)) &&
					(nMidPos <= (vSourceRefGene[ni]->nTxEnd+nDOWN)) )
				{
					nDistMin = (int)fabs(nMidPos-vSourceRefGene[ni]->nTxEnd);
				}
				else
				{
					continue;
				}
			}
		}
		else if(nRefType == 3)
		{
			if( (nMidPos >= vSourceRefGene[ni]->nCdsStart) &&
				(nMidPos <= vSourceRefGene[ni]->nCdsEnd) )
			{
				nDistMin = 0;
			}
			else if(nMidPos <  vSourceRefGene[ni]->nCdsStart)
			{
				nDistMin = vSourceRefGene[ni]->nCdsStart-nMidPos;
				if(vSourceRefGene[ni]->chStrand == '-')
					nDistRef = nDOWN;
				else
					nDistRef = nUP;
				
				if(nDistMin > nDistRef)
					continue;
			}
			else
			{
				nDistMin = nMidPos-vSourceRefGene[ni]->nCdsEnd;
				if(vSourceRefGene[ni]->chStrand == '-')
					nDistRef = nUP;
				else
					nDistRef = nDOWN;
				
				if(nDistMin > nDistRef)
					continue;
			}
		}
		else if(nRefType == 4)
		{
			if(vSourceRefGene[ni]->chStrand == '-')
			{
				if( (nMidPos <= (vSourceRefGene[ni]->nCdsEnd+nUP)) &&
					(nMidPos >= (vSourceRefGene[ni]->nCdsEnd-nDOWN)) )
				{
					nDistMin = (int)fabs(nMidPos-vSourceRefGene[ni]->nCdsEnd);
				}
				else
				{
					continue;
				}
			}
			else
			{
				if( (nMidPos >= (vSourceRefGene[ni]->nCdsStart-nUP)) &&
					(nMidPos <= (vSourceRefGene[ni]->nCdsStart+nDOWN)) )
				{
					nDistMin = (int)fabs(nMidPos-vSourceRefGene[ni]->nCdsStart);
				}
				else
				{
					continue;
				}
			}
		}
		else if(nRefType == 5)
		{
			if(vSourceRefGene[ni]->chStrand == '-')
			{
				if( (nMidPos <= (vSourceRefGene[ni]->nCdsStart+nUP)) &&
					(nMidPos >= (vSourceRefGene[ni]->nCdsStart-nDOWN)) )
				{
					nDistMin = (int)fabs(nMidPos-vSourceRefGene[ni]->nCdsStart);
				}
				else
				{
					continue;
				}
			}
			else
			{
				if( (nMidPos >= (vSourceRefGene[ni]->nCdsEnd-nUP)) &&
					(nMidPos <= (vSourceRefGene[ni]->nCdsEnd+nDOWN)) )
				{
					nDistMin = (int)fabs(nMidPos-vSourceRefGene[ni]->nCdsEnd);
				}
				else
				{
					continue;
				}
			}
		}
		else
		{
			return -1;
		}

		if(nOptId == -1)
		{
			nOptId = ni;
			nOptDist = nDistMin;
			nOptLen = vSourceRefGene[ni]->nTxEnd-vSourceRefGene[ni]->nTxStart+1;
		}
		else
		{
			if(nDistMin < nOptDist)
			{
				nOptId = ni;
				nOptDist = nDistMin;
				nOptLen = vSourceRefGene[ni]->nTxEnd-vSourceRefGene[ni]->nTxStart+1;
			}
			else if(nDistMin == nOptDist)
			{
				nTempLen = vSourceRefGene[ni]->nTxEnd-vSourceRefGene[ni]->nTxStart+1;
				if(nTempLen < nOptLen)
				{
					nOptId = ni;
					nOptDist = nDistMin;
					nOptLen = nTempLen;
				}
			}
		}
	}

	/* return */
	return nOptId;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGene_GetLocationSummary_Main                                        */
/*  get nearest gene from a refgene database.                              */
/* ----------------------------------------------------------------------- */ 
int RefGene_GetLocationSummary_Main(char strDatabasePath[], int nDatabaseType,
			char strSpecies[], char strInputPath[], int nInputType, char strOutputPath[])
{
	/* out */
	FILE *fpIn;
	FILE *fpOut;
	FILE *fpOut2;
	int ni;
	struct tagRefGene **vSourceRefGene = NULL;
	int nSourceRefNum = 0;
	char strLine[LONG_LINE_LENGTH];
	char strSeqName[MED_LINE_LENGTH];
	char strChr[MED_LINE_LENGTH];
	int nChr = -1;
	int nStart = 0;
	int nEnd = 0;
	int nPos = 0;
	char chStrand = '+';
	struct BYTEMATRIX *pLocInfo = NULL;
	struct DOUBLEMATRIX *pLocTotal = NULL;
	int nLocTypeNum = 17;
	int nRegNum = 0;

	/* init */
	/* load source refgene */
	vSourceRefGene = NULL;
	vSourceRefGene = RefGene_LoadDatabase(strDatabasePath, nDatabaseType, 
		strSpecies, &nSourceRefNum);
	if(vSourceRefGene == NULL)
	{
		printf("Warning: RefGene_GetLocationSummary_Main, null database!\n");
		return PROC_SUCCESS;
	}

	/* load coordinates */
	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: RefGene_GetLocationSummary_Main, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}

	if(nInputType == 1)
	{
		sprintf(strLine, "%s.bed", strOutputPath);
	}
	else if(nInputType == 2)
	{
		sprintf(strLine, "%s.codp", strOutputPath);
	}
	else if(nInputType == 3)
	{
		sprintf(strLine, "%s.bedp", strOutputPath);
	}
	else
	{
		sprintf(strLine, "%s.cod", strOutputPath);
	}

	
	fpOut2 = NULL;
	fpOut2 = fopen(strLine, "w");
	if(fpOut2 == NULL)
	{
		printf("Error: RefGene_GetLocationSummary_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: RefGene_GetLocationSummary_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	pLocInfo = NULL;
	pLocInfo = CreateByteMatrix(1, nLocTypeNum);

	pLocTotal = NULL;
	pLocTotal = CreateDoubleMatrix(1, nLocTypeNum);

	if( (pLocInfo == NULL) || (pLocTotal == NULL) )
	{
		printf("Error: RefGene_GetLocationSummary_Main, cannot allocate memory for counting!\n");
		exit(EXIT_FAILURE);
	}


	if(nInputType == 1)
	{
		fprintf(fpOut2, "#chromosome\tstart\tend\tintergenic\tintergenic\tintragenic\texon\tintron\tCDS\tUTR\t5'UTR\t3'UTR\tintergenic(<=1kb-from-gene)\tTSSup1k\tTESdown1k\tintergenic(<=10kb-from-gene)\tTSSup10k\tTESdown10k\tintergenic(<=100kb-from-gene)\tTSSup100k\tTESdown100k\n");
	}
	else if(nInputType == 2)
	{
		fprintf(fpOut2, "#seq_id\tchromosome\tpos\tintergenic\tintragenic\texon\tintron\tCDS\tUTR\t5'UTR\t3'UTR\tintergenic(<=1kb-from-gene)\tTSSup1k\tTESdown1k\tintergenic(<=10kb-from-gene)\tTSSup10k\tTESdown10k\tintergenic(<=100kb-from-gene)\tTSSup100k\tTESdown100k\n");
	}
	else if(nInputType == 3)
	{
		fprintf(fpOut2, "#chromosome\tpos\tintergenic\tintragenic\texon\tintron\tCDS\tUTR\t5'UTR\t3'UTR\tintergenic(<=1kb-from-gene)\tTSSup1k\tTESdown1k\tintergenic(<=10kb-from-gene)\tTSSup10k\tTESdown10k\tintergenic(<=100kb-from-gene)\tTSSup100k\tTESdown100k\n");
	}
	else
	{
		fprintf(fpOut2, "#seq_id\tchromosome\tstart\tend\tstrand\tintergenic\tintragenic\texon\tintron\tCDS\tUTR\t5'UTR\t3'UTR\tintergenic(<=1kb-from-gene)\tTSSup1k\tTESdown1k\tintergenic(<=10kb-from-gene)\tTSSup10k\tTESdown10k\tintergenic(<=100kb-from-gene)\tTSSup100k\tTESdown100k\n");
	}


	/* process one by one */
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		nChr = -1;
		nStart = 0;
		nEnd = 0;
		nPos = 0;
		chStrand = '+';
	
		if(nInputType == 1)
		{
			sscanf(strLine, "%s %d %d", strChr, &nStart, &nEnd);
			nPos = (nStart+nEnd)/2;
		}
		else if(nInputType == 2)
		{
			sscanf(strLine, "%s %s %d", strSeqName, strChr, &nPos);
		}
		else if(nInputType == 3)
		{
			sscanf(strLine, "%s %d", strChr, &nPos);
		}
		else
		{
			sscanf(strLine, "%s %s %d %d %c", strSeqName, strChr, &nStart, &nEnd, &chStrand);
			nPos = (nStart+nEnd)/2;
		}

		nChr = Genome_ChromosomeName_To_Index(strChr, strSpecies);

		memset(pLocInfo->pMatElement, 0, pLocInfo->nWidth);
		RefGene_GetLocationCode(nChr, nPos, vSourceRefGene, nSourceRefNum, pLocInfo);
		for(ni=0; ni<nLocTypeNum; ni++)
		{
			pLocTotal->pMatElement[ni] += pLocInfo->pMatElement[ni];
		}
		
		if(nInputType == 1)
		{
			fprintf(fpOut2, "%s\t%d\t%d", strChr, nStart, nEnd);
		}
		else if(nInputType == 2)
		{
			fprintf(fpOut2, "%s\t%s\t%d", strSeqName, strChr, nPos);
		}
		else if(nInputType == 3)
		{
			fprintf(fpOut2, "%s\t%d", strChr, nPos);
		}
		else
		{
			fprintf(fpOut2, "%s\t%s\t%d\t%d\t%c", strSeqName, strChr, nStart, nEnd, chStrand);
		}

		for(ni=0; ni<nLocTypeNum; ni++)
		{
			fprintf(fpOut2, "\t%d", (int)(pLocInfo->pMatElement[ni]));
		}
		fprintf(fpOut2, "\n");

		nRegNum++;
	}

	/* fprintf(fpOut, "Summary\n"); */
	fprintf(fpOut, "#total\tintergenic\tintragenic\texon\tintron\tCDS\tUTR\t5'UTR\t3'UTR\tintergenic(<=1kb-from-gene)\tTSSup1k\tTESdown1k\tintergenic(<=10kb-from-gene)\tTSSup10k\tTESdown10k\tintergenic(<=100kb-from-gene)\tTSSup100k\tTESdown100k\n");
	fprintf(fpOut, "%d", nRegNum);
	for(ni=0; ni<nLocTypeNum; ni++)
	{
		fprintf(fpOut, "\t%d", (int)(pLocTotal->pMatElement[ni]));
	}
	fprintf(fpOut, "\n");

	fprintf(fpOut, "1.0000");
	if(nRegNum > 0)
	{
		for(ni=0; ni<nLocTypeNum; ni++)
		{
			fprintf(fpOut, "\t%f", pLocTotal->pMatElement[ni]/(double)nRegNum);
		}
	}
	else
	{
		for(ni=0; ni<nLocTypeNum; ni++)
		{
			fprintf(fpOut, "\t0.0000");
		}
	}
	fprintf(fpOut, "\n");

	/* close files */
	fclose(fpIn);
	fclose(fpOut);
	fclose(fpOut2);
	
	/* clear database */
	RefGene_ClearDatabase(&vSourceRefGene, nSourceRefNum);
	DestroyByteMatrix(pLocInfo);
	DestroyDoubleMatrix(pLocTotal);
	pLocInfo = NULL;
	pLocTotal = NULL;

	/* exit */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  RefGene_GetLocationCode                                                */
/*  get location code.                                                     */
/*  intergenic, exon, intron,                                              */
/*  intergenic (1k from TSS or TES), intergenic (other)                    */
/*  exon(5'or 3' UTR), exon(other)                                         */
/*  TSS-up1k, TES-down1K.                                                  */
/*  5'UTR, 3'UTR.                                                          */
/* ----------------------------------------------------------------------- */
int RefGene_GetLocationCode(int nChr, int nPos, struct tagRefGene **vSourceRefGene, 
				int nSourceRefNum, struct BYTEMATRIX *pLocInfo)
{
	/* define */
	int ni,nj,nk,nDist,nx;
	int nP1,nP2,nIsExon;

	/* intergenic\tintragenic\texon\tintron\tCDS\tUTR\t5'UTR\t3'UTR\tintergenic(dist_to_gene<=1kb)\tTSSup1k\tTESdown1k\tintergenic(dist_to_gene<=10kb)\tTSSup10k\tTESdown10k\tintergenic(dist_to_gene<=100kb)\tTSSup100k\tTESdown100k */

	/* search */
	ni = 0;
	nj = nSourceRefNum-1;
	
	while((nj-ni)>1)
	{
		nk = (ni+nj)/2;
		if(nChr < vSourceRefGene[nk]->nChrom)
		{
			nj = nk;
		}
		else if(nChr > vSourceRefGene[nk]->nChrom)
		{
			ni = nk;
		}
		else
		{
			nDist = vSourceRefGene[nk]->nTxStart-nPos;
			if(nDist > 0)
			{
				nj = nk;
			}
			else
			{
				ni = nk;
			}
		}
	}

	if( (vSourceRefGene[ni]->nChrom == nChr) && (vSourceRefGene[nj]->nChrom == nChr) )
	{
		nk = ni;
	}
	else if(vSourceRefGene[ni]->nChrom == nChr)
	{
		nk = ni;
	}
	else if(vSourceRefGene[nj]->nChrom == nChr)
	{
		nk = nj;
	}
	else
	{
		/* if no chromosome match, set intergenic = 1, intergenic(other) = 1 */
		pLocInfo->pMatElement[0] = 1;
		return PROC_SUCCESS;
	}

	/* intergenic\tintragenic\texon\tintron\tCDS\tUTR\t5'UTR\t3'UTR\tintergenic(dist_to_gene<=1kb)\tTSSup1k\tTESdown1k\tintergenic(dist_to_gene<=10kb)\tTSSup10k\tTESdown10k\tintergenic(dist_to_gene<=100kb)\tTSSup100k\tTESdown100k */
	
	while(nk<nSourceRefNum)
	{
		if( (vSourceRefGene[nk]->nChrom != nChr) || ( (vSourceRefGene[nk]->nTxStart-100000) > nPos) )
		{
			break;
		}
		else
		{
			nk++;
		}
	} 
	nk--;

	/* set default as intergenic(other) */
	pLocInfo->pMatElement[0] = 1;

	while(nk>=0)
	{
		if( (vSourceRefGene[nk]->nChrom != nChr) )
		{
			break;
		}
		if( nPos-vSourceRefGene[nk]->nTxEnd > 100000)
		{
			nk--;
			continue;
		}
		
		/* outside gene */
		if( nPos < vSourceRefGene[nk]->nTxStart )
		{
			/* intergenic(dist_to_gene<=100kb) */
			if( (vSourceRefGene[nk]->nTxStart-nPos) <= 100000 )
			{
				pLocInfo->pMatElement[14] = 1;
				/* TSSup100k */
				if(vSourceRefGene[nk]->chStrand == '+')
				{
					pLocInfo->pMatElement[15] = 1;
				}
				else if(vSourceRefGene[nk]->chStrand == '-')
				{
					pLocInfo->pMatElement[16] = 1;
				}

				/* intergenic(dist_to_gene<=10kb) */
				if( (vSourceRefGene[nk]->nTxStart-nPos) <= 10000 )
				{
					pLocInfo->pMatElement[11] = 1;
					/* TSSup1k */
					if(vSourceRefGene[nk]->chStrand == '+')
					{
						pLocInfo->pMatElement[12] = 1;
					}
					else if(vSourceRefGene[nk]->chStrand == '-')
					{
						pLocInfo->pMatElement[13] = 1;
					}
			
					/* intergenic(dist_to_gene<=1kb) */
					if( (vSourceRefGene[nk]->nTxStart-nPos) <= 1000 )
					{
						pLocInfo->pMatElement[8] = 1;
						/* TSSup1k */
						if(vSourceRefGene[nk]->chStrand == '+')
						{
							pLocInfo->pMatElement[9] = 1;
						}
						else if(vSourceRefGene[nk]->chStrand == '-')
						{
							pLocInfo->pMatElement[10] = 1;
						}
					}
				}
			}			
		}
		/* within gene */
		else if(nPos <= vSourceRefGene[nk]->nTxEnd)
		{
			/* intergenic\tintragenic\texon\tintron\tCDS\tUTR\t5'UTR\t3'UTR\tintergenic(dist_to_gene<=1kb)\tTSSup1k\tTESdown1k\tintergenic(dist_to_gene<=10kb)\tTSSup10k\tTESdown10k\tintergenic(dist_to_gene<=100kb)\tTSSup100k\tTESdown100k */
			pLocInfo->pMatElement[1] = 1;

			nIsExon = 0;
			for(nx=0; nx<vSourceRefGene[nk]->nExonCount; nx++)
			{
				nP1 = IMGETAT(vSourceRefGene[nk]->pmatExonStartsEnds, nx, 0);
				nP2 = IMGETAT(vSourceRefGene[nk]->pmatExonStartsEnds, nx, 1);
				if(nPos < nP1)
				{
					break;
				}
				else if(nPos<=nP2)
				{
					/* exon */
					nIsExon = 1;
					break;
				}
			}

			/* exon = 1 */
			if(nIsExon == 1)
			{
				pLocInfo->pMatElement[2] = 1;

				/* CDS */
				if( (nPos >= vSourceRefGene[nk]->nCdsStart) && (nPos <= vSourceRefGene[nk]->nCdsEnd) )
				{
					pLocInfo->pMatElement[4] = 1;
				}
				/* UTR */
				else
				{
					pLocInfo->pMatElement[5] = 1;
					if(nPos < vSourceRefGene[nk]->nCdsStart)
					{
						/* 5'UTR */
						if(vSourceRefGene[nk]->chStrand == '+')
						{
							pLocInfo->pMatElement[6] = 1;
						}
						/* 3'UTR */
						else if(vSourceRefGene[nk]->chStrand == '-')
						{
							pLocInfo->pMatElement[7] = 1;
						}
					}
					else
					{
						/* 3'UTR */
						if(vSourceRefGene[nk]->chStrand == '+')
						{
							pLocInfo->pMatElement[7] = 1;
						}
						/* 5'UTR */
						else if(vSourceRefGene[nk]->chStrand == '-')
						{
							pLocInfo->pMatElement[6] = 1;
						}
					}
				}
			}
			/* intron = 1 */
			else
			{
				pLocInfo->pMatElement[3] = 1;
			}
		}
		/* outside gene */
		else
		{
			/* intergenic(dist_to_gene<=100kb) */
			if( (nPos-vSourceRefGene[nk]->nTxEnd) <= 100000 )
			{
				pLocInfo->pMatElement[14] = 1;
				/* TSSup1k */
				if(vSourceRefGene[nk]->chStrand == '-')
				{
					pLocInfo->pMatElement[15] = 1;
				}
				else if(vSourceRefGene[nk]->chStrand == '+')
				{
					pLocInfo->pMatElement[16] = 1;
				}

				/* intergenic(dist_to_gene<=10kb) */
				if( (nPos-vSourceRefGene[nk]->nTxEnd) <= 10000 )
				{
					pLocInfo->pMatElement[11] = 1;
					/* TSSup1k */
					if(vSourceRefGene[nk]->chStrand == '-')
					{
						pLocInfo->pMatElement[12] = 1;
					}
					else if(vSourceRefGene[nk]->chStrand == '+')
					{
						pLocInfo->pMatElement[13] = 1;
					}

					/* intergenic(dist_to_gene<=1kb) */
					if( (nPos-vSourceRefGene[nk]->nTxEnd) <= 1000 )
					{
						pLocInfo->pMatElement[8] = 1;
						/* TSSup1k */
						if(vSourceRefGene[nk]->chStrand == '-')
						{
							pLocInfo->pMatElement[9] = 1;
						}
						else if(vSourceRefGene[nk]->chStrand == '+')
						{
							pLocInfo->pMatElement[10] = 1;
						}
					}
				}
			}
		}

		nk--;
	}

	/* intergenic\tintragenic\texon\tintron\tCDS\tUTR\t5'UTR\t3'UTR\tintergenic(dist_to_gene<=1kb)\tTSSup1k\tTESdown1k\tintergenic(dist_to_gene<=10kb)\tTSSup10k\tTESdown10k\tintergenic(dist_to_gene<=100kb)\tTSSup100k\tTESdown100k */

	/* exon */
	if(pLocInfo->pMatElement[1] == 1)
	{
		pLocInfo->pMatElement[0] = 0;
	}
		
	/* exit */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGene_GetMatchedControl_Main                                         */
/*  get matched genomic controls.                                          */
/* ----------------------------------------------------------------------- */ 
int RefGene_GetMatchedControl_Main(char strDatabasePath[], int nDatabaseType,
			char strSpecies[], char strChrLenPath[], 
			char strInputPath[], char strOutputPath[],
			int nRepNum, int nRegionLen, int nRemoveRedundancy)
{
	/* out */
	char strFileName[MED_LINE_LENGTH];
	struct tagRefGene **vSourceRefGene = NULL;
	int nSourceRefNum = 0;
	int nRegionNum = 0;
	int nRegionCTNum = 0;
	int nRegionNRNum = 0;
	struct INTMATRIX *pInitDist = NULL;
	struct INTMATRIX *pDist = NULL;

	struct INTMATRIX *pChrLen = NULL;
	struct DOUBLEMATRIX *pRegionCT = NULL;
	struct INTMATRIX *pType = NULL;
	struct INTMATRIX *pPriority = NULL;
	struct DOUBLEMATRIX *pRegionSort = NULL;
	struct LONGMATRIX *pRegionSid = NULL;
	int ni,nj,nk;

	FILE *fpIn;
	FILE *fpOut;
	
	char strLine[LONG_LINE_LENGTH];
	char strSeqName[MED_LINE_LENGTH];
	char strChr[MED_LINE_LENGTH];
	int nChr,nStart,nEnd;
	char chStrand;
	int nLocusID = -1;
	char strGeneName[MED_LINE_LENGTH];
	char strRefID[MED_LINE_LENGTH];
	int nGChr,nGStart,nGEnd;
	char chGStrand;
	char *chp1,*chp2;

	double dRand;
	int nRand;
	int nHalfRegionLen;
	int nTemp,nTempChr;

	/* init */
	nHalfRegionLen = nRegionLen/2;

	/* get region number */
	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: RefGene_GetMatchedControl_Main, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}

	nRegionNum = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		nRegionNum++;
	}
	fclose(fpIn);

	if(nRegionNum <= 0)
	{
		printf("Warning: RefGene_GetMatchedControl_Main, no input regions!\n");
		return PROC_SUCCESS;
	}

	/* get annotation */
	sprintf(strFileName, "%s.gene", strOutputPath);
	RefGene_GetNearestGene_Main(strDatabasePath, nDatabaseType, strSpecies, 
		strInputPath, strFileName, 0, 100000000, 100000000);

	/* get reference distance */
	pInitDist = NULL;
	pInitDist = CreateIntMatrix(1, nRegionNum);
	if(pInitDist == NULL)
	{
		printf("Error: RefGene_GetMatchedControl_Main, cannot create memory for storing peak-gene distances!\n");
		return PROC_SUCCESS;
	}

	fpIn = NULL;
	fpIn = fopen(strFileName, "r");
	if(fpIn == NULL)
	{
		printf("Error: RefGene_GetMatchedControl_Main, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		chp1 = strchr(strLine, '\t');
		*chp1 = '\0';
		strcpy(strSeqName, strLine);
		chp2 = chp1+1;

		chp1 = strchr(chp2, '\t');
		*chp1 = '\0';
		strcpy(strChr, chp2);
		chp2 = chp1+1;
		nChr = Genome_ChromosomeName_To_Index(strChr, strSpecies);
		
		chp1 = strchr(chp2, '\t');
		*chp1 = '\0';
		nStart = atoi(chp2);
		chp2 = chp1+1;

		chp1 = strchr(chp2, '\t');
		*chp1 = '\0';
		nEnd = atoi(chp2);
		chp2 = chp1+1;
		
		chp1 = strchr(chp2, '\t');
		*chp1 = '\0';
		chStrand = *chp2;
		chp2 = chp1+1;

		if(strstr(chp2, "---") != NULL)
		{
			continue;
		}

		if(nDatabaseType >= 2)
		{
			chp1 = strchr(chp2, '\t');
			*chp1 = '\0';
			nLocusID = atoi(chp2);
			chp2 = chp1+1;
		}

		chp1 = strchr(chp2, '\t');
		*chp1 = '\0';
		strcpy(strGeneName, chp2);
		chp2 = chp1+1;

		sscanf(chp2, "%s %d %c %d %d", strRefID, &nGChr, &chGStrand,
			&nGStart, &nGEnd);

        if(nGChr != nChr)
		{
			printf("Error: RefGene_GetMatchedControl_Main, chromosome not match!\n");
			exit(EXIT_FAILURE);
		}

		if(chGStrand == '-')
		{
			pInitDist->pMatElement[ni] = nGEnd-(nStart+nEnd)/2;
		}
		else
		{
			pInitDist->pMatElement[ni] = (nStart+nEnd)/2-nGStart;
		}

		ni++;
	}
	fclose(fpIn);

	if( (ni <= 0) || (ni > nRegionNum) )
	{
		printf("Error: RefGene_GetMatchedControl_Main, peak-gene annotation error!\n");
		exit(EXIT_FAILURE);
	}

	pDist = NULL;
	pDist = CreateIntMatrix(1, ni);
	if(pDist == NULL)
	{
		printf("Error: RefGene_GetMatchedControl_Main, cannot create memory for storing peak-gene distances!\n");
		return PROC_SUCCESS;
	}

	memcpy(pDist->pMatElement, pInitDist->pMatElement, sizeof(int)*pDist->nWidth);
	nRegionNum = ni;
	nRegionCTNum = nRegionNum*nRepNum;
	DestroyIntMatrix(pInitDist);

	/* load source refgene */
	vSourceRefGene = NULL;
	vSourceRefGene = RefGene_LoadDatabase(strDatabasePath, nDatabaseType, 
		strSpecies, &nSourceRefNum);
	if(vSourceRefGene == NULL)
	{
		printf("Warning: RefGene_GetMatchedControl_Main, null database!\n");
		return PROC_SUCCESS;
	}

	/* load chromose length */
	pChrLen = NULL;
	pChrLen = IMLOAD(strChrLenPath);
	if(pChrLen == NULL)
	{
		printf("Error: RefGene_GetMatchedControl_Main, cannot find chromosome length!\n");
		return PROC_SUCCESS;
	}

	/* create region */
	pRegionCT = NULL;
	pRegionCT = CreateDoubleMatrix(nRegionCTNum, 3);
	if(pRegionCT == NULL)
	{
		printf("Error: RefGene_GetMatchedControl_Main, cannot create memory for regions!\n");
		return PROC_SUCCESS;
	}

	/* choose regions */
	nk = 0;
	for(ni=0; ni<nRepNum; ni++)
	{
		nj = 0;
		while(nj<nRegionNum)
		{
			dRand = rand_u();
			nRand = (int)(dRand*nSourceRefNum);
			if(nRand >= nSourceRefNum)
				nRand = nSourceRefNum-1;

			nTempChr = vSourceRefGene[nRand]->nChrom;
			if( (nTempChr <= 0) || (nTempChr > pChrLen->nHeight) )
			{
				continue;
			}

			DMSETAT(pRegionCT, nk, 0, nTempChr);
			if(vSourceRefGene[nRand]->chStrand == '-')
			{
				nTemp = vSourceRefGene[nRand]->nTxEnd-pDist->pMatElement[nj]-nHalfRegionLen;
				if(nTemp < 0)
					nTemp = 0;
				if(nTemp >= pChrLen->pMatElement[nTempChr-1])
					nTemp = pChrLen->pMatElement[nTempChr-1]-1;
				DMSETAT(pRegionCT, nk, 1, nTemp);

				nTemp = vSourceRefGene[nRand]->nTxEnd-pDist->pMatElement[nj]+nHalfRegionLen;
				if(nTemp < 0)
					nTemp = 0;
				if(nTemp >= pChrLen->pMatElement[nTempChr-1])
					nTemp = pChrLen->pMatElement[nTempChr-1]-1;
				DMSETAT(pRegionCT, nk, 2, nTemp);
			}
			else
			{
				nTemp = vSourceRefGene[nRand]->nTxStart+pDist->pMatElement[nj]-nHalfRegionLen;
				if(nTemp < 0)
					nTemp = 0;
				if(nTemp >= pChrLen->pMatElement[nTempChr-1])
					nTemp = pChrLen->pMatElement[nTempChr-1]-1;				
				DMSETAT(pRegionCT, nk, 1, nTemp);

				nTemp = vSourceRefGene[nRand]->nTxStart+pDist->pMatElement[nj]+nHalfRegionLen;
				if(nTemp < 0)
					nTemp = 0;
				if(nTemp >= pChrLen->pMatElement[nTempChr-1])
					nTemp = pChrLen->pMatElement[nTempChr-1]-1;
				DMSETAT(pRegionCT, nk, 2, nTemp);
			}

			nj++;
			nk++;
		}
	}

	if(nk != nRegionCTNum)
	{
		printf("Error: RefGene_GetMatchedControl_Main, region number do not match!\n");
		exit(EXIT_FAILURE);
	}

	/* clear database */
	DestroyIntMatrix(pDist);
	RefGene_ClearDatabase(&vSourceRefGene, nSourceRefNum);
	DestroyIntMatrix(pChrLen);

	/* sort regions and remove redundancy */
	pType = NULL;
	pType = CreateIntMatrix(1, 3);
	if(pType == NULL)
	{
		printf("Error: RefGene_GetMatchedControl_Main, cannot create memory for regions!\n");
		return PROC_SUCCESS;
	}
	for(ni=0; ni<pType->nWidth; ni++)
		pType->pMatElement[ni] = 2;

	pPriority = NULL;
	pPriority = CreateIntMatrix(1, 3);
	if(pPriority == NULL)
	{
		printf("Error: RefGene_GetMatchedControl_Main, cannot create memory for regions!\n");
		return PROC_SUCCESS;
	}
	pPriority->pMatElement[0] = 0;
	pPriority->pMatElement[1] = 1;
	pPriority->pMatElement[2] = 2;

	DMSORTROWS(pRegionCT, pType, pPriority, &pRegionSort, &pRegionSid);

	/* remove redundancy */
	if(nRemoveRedundancy == 0)
	{
		nRegionNRNum = nRegionCTNum;
	}
	else
	{
		nGChr = (int)DMGETAT(pRegionSort, 0, 0);
		nGStart = (int)DMGETAT(pRegionSort, 0, 1);
		nGEnd = (int)DMGETAT(pRegionSort, 0, 2);
		nRegionNRNum = 1;

		for(ni=1; ni<nRegionCTNum; ni++)
		{
			nChr = (int)DMGETAT(pRegionSort, ni, 0);
			nStart = (int)DMGETAT(pRegionSort, ni, 1);
			nEnd = (int)DMGETAT(pRegionSort, ni, 2);

			if(nChr != nGChr)
			{
				DMSETAT(pRegionSort, nRegionNRNum, 0, nChr);
				DMSETAT(pRegionSort, nRegionNRNum, 1, nStart);
				DMSETAT(pRegionSort, nRegionNRNum, 2, nEnd);
				nGChr = nChr;
				nGStart = nStart;
				nGEnd = nEnd;
				nRegionNRNum++;
			}
			else
			{
				if(nStart > nGEnd)
				{
					DMSETAT(pRegionSort, nRegionNRNum, 0, nChr);
					DMSETAT(pRegionSort, nRegionNRNum, 1, nStart);
					DMSETAT(pRegionSort, nRegionNRNum, 2, nEnd);
					nGChr = nChr;
					nGStart = nStart;
					nGEnd = nEnd;
					nRegionNRNum++;
				}
				else
				{
					if(nEnd > nGEnd)
					{
						nStart = nGEnd+1;
						DMSETAT(pRegionSort, nRegionNRNum, 0, nChr);
						DMSETAT(pRegionSort, nRegionNRNum, 1, nStart);
						DMSETAT(pRegionSort, nRegionNRNum, 2, nEnd);
						nGChr = nChr;
						nGStart = nStart;
						nGEnd = nEnd;
						nRegionNRNum++;
					}
					else
					{
						nRegionNRNum = nRegionNRNum;
					}
				}
			}
		}
	}

	/* output */
	fpOut = NULL;
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: RefGene_GetMatchedControl_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	fprintf(fpOut, "#seq_id\tchromosome\tstart\tend\tstrand\n");

	for(ni=0; ni<nRegionNRNum; ni++)
	{
		nTempChr = (int)DMGETAT(pRegionSort, ni, 0);
		nStart = (int)DMGETAT(pRegionSort, ni, 1);
		nEnd = (int)DMGETAT(pRegionSort, ni, 2);
		Genome_Index_To_ChromosomeName(strChr, strSpecies, nTempChr);
		fprintf(fpOut, "%d\t%s\t%d\t%d\t+\n", ni, strChr, nStart, nEnd);
	}
	
	fclose(fpOut);

	/* release memory */
	DestroyDoubleMatrix(pRegionCT);
	DestroyIntMatrix(pType);
	DestroyIntMatrix(pPriority);
	DestroyDoubleMatrix(pRegionSort);
	DestroyLongMatrix(pRegionSid);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGene_GetNeighborGenes_Main                                          */
/*  get neighboring genes from a refgene database.                         */
/* ----------------------------------------------------------------------- */ 
int RefGene_GetNeighborGenes_Main(char strDatabasePath[], int nDatabaseType,
			char strSpecies[], char strInputPath[], char strOutputPath[],
			char strAnnotationPath[], int nUPNum, int nDOWNNum, 
			int nDistanceUpperLimit)
{
	/* out */
	FILE *fpIn;
	FILE *fpOut;
	int ni,nUpK,nDownK;
	struct tagRefGene **vSourceRefGene = NULL;
	struct tagString **vAnnotation = NULL;
	int nSourceRefNum = 0;
	char strLine[LONG_LINE_LENGTH];
	char strSeqName[MED_LINE_LENGTH];
	char strChr[MED_LINE_LENGTH];
	int nChr,nStart,nEnd;
	char chStrand;
	struct tagRefGene *pRefGene;
	int nFindNum;
	int nTSSDist,nTESDist;
	char strLocType[LINE_LENGTH];

	/* init */
	/* load source refgene */
	vSourceRefGene = NULL;
	vSourceRefGene = RefGene_LoadDatabase(strDatabasePath, nDatabaseType, 
		strSpecies, &nSourceRefNum);
	if(vSourceRefGene == NULL)
	{
		printf("Warning: RefGene_GetNeighborGenes_Main, null database!\n");
		return PROC_SUCCESS;
	}

	/* load annotation */
	fpIn = NULL;
	fpIn = fopen(strAnnotationPath, "r");
	if(fpIn != NULL)
	{
		vAnnotation = (struct tagString **)calloc(nSourceRefNum, sizeof(struct tagString *));
		if(vAnnotation == NULL)
		{
			printf("Error: RefGene_GetNeighborGenes_Main, cannot create memory for loading annotations!\n");
			exit(EXIT_FAILURE);
		}

		ni = 0;
		while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;

			if(ni >= nSourceRefNum)
			{
				printf("Error: RefGene_GetNeighborGenes_Main, annotation/database dimensions do not match!\n");
				exit(EXIT_FAILURE);
			}

			StringAddTail(vAnnotation+ni, strLine);

			ni++;
		}

		if(ni != nSourceRefNum)
		{
			printf("Error: RefGene_GetNeighborGenes_Main, annotation/database dimensions do not match!\n");
			exit(EXIT_FAILURE);
		}

		fclose(fpIn);
	}
	else
	{
		printf("Proceed with no additional annotations!\n"); 
	}

	/* load coordinates */
	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: RefGene_GetNeighborGenes_Main, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: RefGene_GetNeighborGenes_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}
	fprintf(fpOut, "#seq_id\tchromosome\tstart\tend\tstrand\tgene_id\tgene_name\trefseq_id\tdistance_to_TSS\tdistance_to_TES\tlocation\tgene_chr\tgene_strand\tgene_TSS\tgene_TES\tgene_CDSS\tgene_CDSE\tannotation\n");

	/* process one by one */
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		sscanf(strLine, "%s %s %d %d %c", strSeqName, strChr, &nStart, &nEnd, &chStrand);
		nChr = Genome_ChromosomeName_To_Index(strChr, strSpecies);
		RefGene_GetNeighborGenes(nChr, nStart, nEnd, chStrand, 
			vSourceRefGene, nSourceRefNum, nUPNum, nDOWNNum, nDistanceUpperLimit,
			&nUpK, &nDownK);

		nFindNum = 0;
		for(ni=nUpK; ni<=nDownK; ni++)
		{
			if( (ni<0) || (ni>=nSourceRefNum) )
				continue;
			
			pRefGene = vSourceRefGene[ni];
			if(pRefGene == NULL)
				continue;
			
			RefGene_GetLocationInfo(nChr, nStart, nEnd, chStrand, 
			pRefGene, &nTSSDist, &nTESDist, strLocType);

			fprintf(fpOut, "%s\t%s\t%d\t%d\t%c\t%d\t%s\t%s\t%d\t%d\t%s\t%d\t%c\t%d\t%d\t%d\t%d", 
				strSeqName, strChr, nStart, nEnd, chStrand,
				pRefGene->nGeneID, pRefGene->strGene, pRefGene->strName,
				nTSSDist, nTESDist, strLocType, pRefGene->nChrom, pRefGene->chStrand,
				pRefGene->nTxStart, pRefGene->nTxEnd, 
				pRefGene->nCdsStart, pRefGene->nCdsEnd);

			if(vAnnotation != NULL)
			{
				if(vAnnotation[ni] != NULL)
				{
					fprintf(fpOut, "\t%s\n", vAnnotation[ni]->m_pString);
				}
				else
				{
					fprintf(fpOut, "\t---\n");
				}
			}
			else
			{
				fprintf(fpOut, "\n");
			}

			nFindNum++;
		}

		if(nFindNum == 0)
		{
			fprintf(fpOut, "%s\t%s\t%d\t%d\t%c\t-1\t---\t---\t---\t---\t---\t-1\t?\t---\t---\t---\t---", 
				strSeqName, strChr, nStart, nEnd, chStrand);
			if(vAnnotation != NULL)
			{
				fprintf(fpOut, "\t---\n");
			}
			else
			{
				fprintf(fpOut, "\n");
			}
		}

		fprintf(fpOut, "\n");
	}

	/* close files */
	fclose(fpIn);
	fclose(fpOut);
	
	/* clear database */
	if(vAnnotation != NULL)
	{
		for(ni=0; ni<nSourceRefNum; ni++)
		{
			DeleteString(vAnnotation[ni]);
			vAnnotation[ni] = NULL;
		}
		free(vAnnotation);
	}
	RefGene_ClearDatabase(&vSourceRefGene, nSourceRefNum);
	
	/* exit */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGene_GetNeighborGenes                                               */
/*  get neighbor genes from a refgene database.                            */
/* ----------------------------------------------------------------------- */ 
int RefGene_GetNeighborGenes(int nChr, int nStart, int nEnd, char chStrand, 
			struct tagRefGene **vSourceRefGene, int nSourceRefNum,
			int nUPNum, int nDOWNNum, int nDistanceUpperLimit,
			int *pUpK, int *pDownK)
{
	/* define */
	int nOptId = -1;
	int nOptDist = IM_ACCESS_VIOLATION-1;
	
	int ni,nj,nk,nMidPos,nDist;
	int nLastId;

	/* init */
	nMidPos = (nStart+nEnd)/2;

	/* search */
	ni = 0;
	nj = nSourceRefNum-1;
	
	while((nj-ni)>1)
	{
		nk = (ni+nj)/2;
		if(vSourceRefGene[nk]->nChrom > nChr)
		{
			nj = nk;
		}
		else if(vSourceRefGene[nk]->nChrom < nChr)
		{
			ni = nk;
		}
		else
		{
			nDist = vSourceRefGene[nk]->nTxStart-nMidPos;
			if(nDist > 0)
			{
				nj = nk;
			}
			else
			{
				ni = nk;
			}
		}
	}

	if( (vSourceRefGene[ni]->nChrom == nChr) && (vSourceRefGene[nj]->nChrom == nChr) )
	{
		nk = ni;
	}
	else if(vSourceRefGene[ni]->nChrom == nChr)
	{
		nk = ni;
	}
	else if(vSourceRefGene[nj]->nChrom == nChr)
	{
		nk = nj;
	}
	else
	{
		*pUpK = nj;
		*pDownK = nj-1;
		return PROC_SUCCESS;
	}

	while(nk<nSourceRefNum)
	{
		if( (vSourceRefGene[nk]->nChrom != nChr) || ((vSourceRefGene[nk]->nTxStart-nMidPos)>0) )
		{
			break;
		}
		else
		{
			nk++;
		}
	}
	    
	*pDownK = nk;
	nLastId = -2;
	nj = 0;
	while(nj<nDOWNNum)
	{
		if((*pDownK) >= nSourceRefNum)
		{
			break;
		}
		if(vSourceRefGene[*pDownK]->nChrom != nChr)
		{
			break;
		}
		if( (vSourceRefGene[*pDownK]->nTxStart - nMidPos) > nDistanceUpperLimit )
		{
			break;
		}
		
		if(vSourceRefGene[*pDownK]->nGeneID != nLastId)
		{
			nLastId = vSourceRefGene[*pDownK]->nGeneID;
			nj++;
		}
        
		*pDownK += 1;
	}
	*pDownK -= 1;
	
	nk--;
	
	while(nk>=0)
	{
		if( (vSourceRefGene[nk]->nChrom != nChr) || ((vSourceRefGene[nk]->nTxEnd-nMidPos)<0) )
		{
			break;
		}
		else
		{
			nk--;
		}
	}
	
	*pUpK = nk;
	nLastId = -2;
	nj = 0;
	while(nj<nUPNum)
	{
		if((*pUpK) < 0)
		{
			break;
		}
		if(vSourceRefGene[*pUpK]->nChrom != nChr)
		{
			break;
		}
		if( (nMidPos-vSourceRefGene[*pUpK]->nTxEnd) > nDistanceUpperLimit )
		{
			break;
		}

		if(vSourceRefGene[*pUpK]->nGeneID != nLastId)
		{
			nLastId = vSourceRefGene[*pUpK]->nGeneID;
			nj++;
		}
        
		*pUpK -= 1;
	}
	*pUpK += 1;


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGene_GetLocationInfo                                                */
/*  get relative location info of a region.                                */
/* ----------------------------------------------------------------------- */ 
int RefGene_GetLocationInfo(int nChr, int nStart, int nEnd, char chStrand, 
			struct tagRefGene *pRefGene, int *pTSSDist, int *pTESDist,
			char strLocType[])
{
	/* define */
	int nMidPos;
	int ni;
	int nP1,nP2;

	/* init check */
	*pTSSDist = 0;
	*pTESDist = 0;
	strcpy(strLocType, "NA");
	if(pRefGene == NULL)
	{
		return PROC_FAILURE;
	}
	if(pRefGene->nChrom != nChr)
	{
		strcpy(strLocType, "different_chr");
		return PROC_FAILURE;
	}

	nMidPos = (nStart+nEnd)/2;
	if(pRefGene->chStrand == '-')
	{
		*pTSSDist = pRefGene->nTxEnd-nMidPos;
		*pTESDist = pRefGene->nTxStart-nMidPos;
	}
	else
	{
		*pTSSDist = nMidPos-pRefGene->nTxStart;
		*pTESDist = nMidPos-pRefGene->nTxEnd;
	}

	if(nMidPos < pRefGene->nTxStart)
	{
		if(pRefGene->chStrand == '-')
		{
			strcpy(strLocType, "TES_downstream");
		}
		else
		{
			strcpy(strLocType, "TSS_upstream");
		}
	}
	else if(nMidPos > pRefGene->nTxEnd)
	{
		if(pRefGene->chStrand == '-')
		{
			strcpy(strLocType, "TSS_upstream");
		}
		else
		{
			strcpy(strLocType, "TES_downstream");
		}
	}
	else if(nMidPos < pRefGene->nCdsStart)
	{
		if(pRefGene->chStrand == '-')
		{
			strcpy(strLocType, "3'UTR");
		}
		else
		{
			strcpy(strLocType, "5'UTR");
		}
	}
	else if(nMidPos > pRefGene->nCdsEnd)
	{
		if(pRefGene->chStrand == '-')
		{
			strcpy(strLocType, "5'UTR");
		}
		else
		{
			strcpy(strLocType, "3'UTR");
		}
	}
	else
	{
		for(ni=0; ni<pRefGene->nExonCount; ni++)
		{
			nP1 = IMGETAT(pRefGene->pmatExonStartsEnds, ni, 0);
			nP2 = IMGETAT(pRefGene->pmatExonStartsEnds, ni, 1);
			if( (nMidPos>=nP1) && (nMidPos<=nP2) )
			{
				strcpy(strLocType, "exon");
				break;
			}
			else if(nMidPos < nP1)
			{
				strcpy(strLocType, "intron");
				break;
			}
		}
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGene_LoadDatabase:                                                  */
/*  init refgene database from refGene(0)/refFlat(1)/refLocus(2) file.     */
/* ----------------------------------------------------------------------- */ 
struct tagRefGene **RefGene_LoadDatabase(char strDatabasePath[], 
					int nDatabaseType, char strSpecies[], int *pSourceRefNum)
{
	/* define */
	FILE *fpRefGene;
	struct tagRefGene *pRefGene,*pCurrentRefGene;
	struct tagRefGene *pSourceRefGeneList;
	struct tagRefGene **vSourceRefGene = NULL;
	int nSourceRefNum = 0;
	char strRefLine[LONG_LINE_LENGTH];
	int ni;

	/* load */
	fpRefGene = NULL;
	fpRefGene = fopen(strDatabasePath, "r");
	if(fpRefGene == NULL)
	{
		printf("Error: RefGene_LoadDatabase, cannot open source refgene file!\n");
		exit(EXIT_FAILURE);
	}

	nSourceRefNum = 0;
	pCurrentRefGene = NULL;
	pSourceRefGeneList = NULL;
	while(fgets(strRefLine, LONG_LINE_LENGTH, fpRefGene) != NULL)
	{
		StrTrimLeft(strRefLine);
		StrTrimRight(strRefLine);
		if(strRefLine[0] == '\0')
			continue;

		pRefGene = NULL;
		pRefGene = RefGeneCreate();
		if(nDatabaseType == 1)
		{
			RefFlatInit_FromGenomeLabFormat(pRefGene, strRefLine, strSpecies);
		}
		else if(nDatabaseType == 2)
		{
			RefLocusInit_FromGenomeLabFormat(pRefGene, strRefLine, strSpecies);
		}
		else
		{
			RefGeneInit_FromGenomeLabFormat(pRefGene, strRefLine, strSpecies);
		}
		if(pSourceRefGeneList == NULL)
		{
			pSourceRefGeneList = pRefGene;
			pCurrentRefGene = pRefGene;
		}
		else
		{
			pCurrentRefGene->pNext = pRefGene;
			pCurrentRefGene = pRefGene;
		}
		nSourceRefNum++;
	}
	fclose(fpRefGene);

	if(nSourceRefNum <= 0)
	{
		printf("Warning: RefGene_LoadDatabase, null database!\n");
		return NULL;
	}

	vSourceRefGene = NULL;
	vSourceRefGene = (struct tagRefGene **)calloc(nSourceRefNum, sizeof(struct tagRefGene*));
	if(vSourceRefGene == NULL)
	{
		printf("Error: RefGene_LoadDatabase, organize refgene into array!\n");
		exit(EXIT_FAILURE);
	}

	ni=0;
	while(pSourceRefGeneList != NULL)
	{
		pRefGene = pSourceRefGeneList;
		pSourceRefGeneList = pRefGene->pNext;
		pRefGene->pNext = NULL;
		vSourceRefGene[ni] = pRefGene;
		ni++;
	}

	if(ni != nSourceRefNum)
	{
		printf("Error: RefGene_LoadDatabase, refgene number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* return */
	*pSourceRefNum = nSourceRefNum;
	return vSourceRefGene;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGene_ClearDatabase:                                                 */
/*  clear refgene database.                                                */
/* ----------------------------------------------------------------------- */ 
int RefGene_ClearDatabase(struct tagRefGene ***vSourceRefGene, int nSourceRefNum)
{
	/* define */
	int ni;
	struct tagRefGene **vDatabase;

	/* clear */
	vDatabase = *vSourceRefGene;
	for(ni=0; ni<nSourceRefNum; ni++)
	{
		RefGeneDestroy(vDatabase[ni]);
		vDatabase[ni] = NULL;
	}
	free(*vSourceRefGene);
	*vSourceRefGene = NULL;

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGene_GetTargetGenomicRegion_TSSTES                                  */
/*  get non-redundant genomic region from refgene list in genomelab        */
/*  refgene format.                                                        */
/* ----------------------------------------------------------------------- */ 
int RefGene_GetTargetGenomicRegion_TSSTES(char strRefGenePath[], 
			char strOutPath[], int nTSSUP, int nTESDOWN, 
			char strSpecies[], int nChrNum, char strChrLenFile[])
{

	/* define */
	FILE *fpIn;
	FILE *fpOut;

	struct INTMATRIX *pChrLen;
	struct tagGenomicRegion **vGenomeTarget;
	struct tagGenomicRegion *pRegion;
	int ni,nChrId;
	int nDiscard;
	int nLastStart, nLastEnd;

	struct tagRefGene *pRefGene;
	char strLine[LONG_LINE_LENGTH];
	char strChrName[LINE_LENGTH];

	/* init */
	if(nChrNum <= 0)
	{
		printf("Warning: RefGene_GetTargetGenomicRegion_TSSTES, chromosome number <= 0!\n");
		return PROC_FAILURE;
	}
	vGenomeTarget = NULL;
	vGenomeTarget = (struct tagGenomicRegion **)calloc(nChrNum, sizeof(struct tagGenomicRegion *));
	if(vGenomeTarget == NULL)
	{
		printf("Error: RefGene_GetTargetGenomicRegion_TSSTES, cannot create enough memory for genomic target!\n");
		exit(EXIT_FAILURE);
	}

	pChrLen = NULL;
	pChrLen = IMLOAD(strChrLenFile);
	if(pChrLen == NULL)
	{
		printf("Error: RefGene_GetTargetGenomicRegion_TSSTES, cannot load chromosome length!\n");
		exit(EXIT_FAILURE);
	}

	/* load refgene and get target genomic region */
	fpIn = NULL;
	fpIn = fopen(strRefGenePath, "rt");
	if(fpIn == NULL)
	{
		printf("Error: RefGene_GetTargetGenomicRegion_TSSTES, cannot open refgene file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		pRefGene = NULL;
		pRefGene = RefGeneCreate();
		RefGeneInit_FromGenomeLabFormat(pRefGene, strLine, strSpecies);

		ni = pRefGene->nChrom-1;
		if((ni<0) || (ni>=nChrNum))
		{
			printf("Warning: RefGene_GetTargetGenomicRegion_TSSTES, chromosome index out of range!\n");
			continue;
		}

		pRegion = NULL;
		pRegion = GenomicRegionCreate();
		if(pRegion == NULL)
		{
			printf("Error: RefGene_GetTargetGenomicRegion_TSSTES, cannot create new genomic region!\n");
			exit(EXIT_FAILURE);;
		}

		nDiscard = 0;

		strcpy(pRegion->strChrom, pRefGene->strChrom);
		pRegion->nChrom = pRefGene->nChrom;
		pRegion->chStrand = pRefGene->chStrand;
		if(pRefGene->chStrand == '+')
		{
			pRegion->nStart = pRefGene->nTxStart-nTSSUP;
			pRegion->nEnd = pRefGene->nTxEnd+nTESDOWN;
		}
		else if(pRefGene->chStrand == '-')
		{
			pRegion->nStart = pRefGene->nTxStart-nTESDOWN;
			pRegion->nEnd = pRefGene->nTxEnd+nTSSUP;
		}
		else
		{
			nDiscard = 1;
		}

		if(nDiscard == 1)
		{
			GenomicRegionDestroy(pRegion);
		}
		else
		{
			nChrId = pRefGene->nChrom-1;
			if(pRegion->nStart < 0)
				pRegion->nStart = 0;
			if(pRegion->nEnd >= pChrLen->pMatElement[nChrId])
				pRegion->nEnd = pChrLen->pMatElement[nChrId]-1;
			GenomicRegionInsert((vGenomeTarget+ni), pRegion);
		}

		RefGeneDestroy(pRefGene);
	}

	fclose(fpIn);

	/* write */
	fpOut = NULL;
	fpOut = fopen(strOutPath, "wt");
	if(fpOut == NULL)
	{
		printf("Error: RefGene_GetTargetGenomicRegion_TSSTES, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nChrNum; ni++)
	{
		Genome_Index_To_ChromosomeName(strChrName, strSpecies, (ni+1));
		nLastStart = -1;
		nLastEnd = -1;
		pRegion = vGenomeTarget[ni];
		while(pRegion != NULL)
		{
			if(pRegion->nStart > nLastEnd)
			{
				if(nLastStart >= 0)
				{
					fprintf(fpOut, "%s\t%d\t%d\n", strChrName, nLastStart, nLastEnd);
				}

				nLastStart = pRegion->nStart;
				nLastEnd = pRegion->nEnd;
			}
			else
			{
				if(pRegion->nEnd > nLastEnd)
					nLastEnd = pRegion->nEnd;
			}


			/* get next */
			vGenomeTarget[ni] = pRegion->pNext;
			GenomicRegionDestroy(pRegion);
			pRegion = vGenomeTarget[ni];
		}

		if(nLastStart >= 0)
		{
			fprintf(fpOut, "%s\t%d\t%d\n", strChrName, nLastStart, nLastEnd);
		}
	}

	fclose(fpOut);

	/* release memory */
	free(vGenomeTarget);
	DestroyIntMatrix(pChrLen);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGeneSelectRegion():                                                 */
/*  select regions according to refgene structure.                         */
/* ----------------------------------------------------------------------- */ 
struct tagGenomicRegion *RefGeneSelectRegion(struct tagRefGene *pRefGene,
			int nTSSUP, int nTSSDOWN, int nTESUP, int nTESDOWN, 
			int nIncludeIntron, int nIncludeExon, struct INTMATRIX *pChrLen)
{
	/* define */
	struct tagGenomicRegion *pSelRegion;
	struct tagGenomicRegion *pRegion;
	int nRegionStart,nRegionEnd;
	int nChrId;
	int npp1,npp2;
	int ni;

	/* init */
	if(pRefGene == NULL)
		return NULL;
	if((nTSSUP < 0) || (nTSSDOWN < 0) || (nTESUP < 0) || (nTESDOWN < 0) )
	{
		printf("Error: RefGeneSelectRegion, coordinates wrong!\n");
		exit(EXIT_FAILURE);
	}

	pSelRegion = NULL;

	/* get region */
	if(pRefGene->chStrand == '-')
	{
		nRegionStart = pRefGene->nTxEnd-nTSSDOWN;
		nRegionEnd = pRefGene->nTxEnd+nTSSUP;
		npp1 = pRefGene->nTxStart-nTESDOWN;
		if(npp1 < nRegionStart)
			nRegionStart = npp1;
	}
	else
	{
		nRegionStart = pRefGene->nTxStart-nTSSUP;
		nRegionEnd = pRefGene->nTxStart+nTSSDOWN;
		npp1 = pRefGene->nTxEnd+nTESDOWN;
		if(npp1 > nRegionEnd)
			nRegionEnd = npp1;
	}

	nChrId = pRefGene->nChrom-1;
	if(nRegionStart < 0)
		nRegionStart = 0;
	if(nRegionEnd >= pChrLen->pMatElement[nChrId])
		nRegionEnd = pChrLen->pMatElement[nChrId]-1;

	/* if whole gene */
	if((nIncludeIntron == 1) && (nIncludeExon == 1))
	{
		pRegion = NULL;
		pRegion = GenomicRegionCreate();
		if(pRegion == NULL)
		{
			printf("Error: RefGeneSelectRegion, cannot create the selected region!\n");
			exit(EXIT_FAILURE);
		}
		pRegion->nChrom = pRefGene->nChrom;
		pRegion->chStrand = pRefGene->chStrand;
		pRegion->nStart = nRegionStart;
		pRegion->nEnd = nRegionEnd;
		strcpy(pRegion->strChrom, pRefGene->strChrom);
		
		pSelRegion = pRegion;
	}

	/* if part of the gene */
	else
	{
		/* if + strand */
		if(pRefGene->chStrand == '+')
		{
			/* add TSS */
			if( (nTSSUP != 0) || (nTSSDOWN != 0) )
			{
				npp1 = pRefGene->nTxStart-nTSSUP;
				npp2 = pRefGene->nTxStart+nTSSDOWN;
				if(npp1 < nRegionStart)
					npp1 = nRegionStart;
				if(npp2 > nRegionEnd)
					npp2 = nRegionEnd;

				pRegion = NULL;
				pRegion = GenomicRegionCreate();
				if(pRegion == NULL)
				{
					printf("Error: RefGeneSelectRegion, cannot create the selected region!\n");
					exit(EXIT_FAILURE);
				}
				pRegion->nChrom = pRefGene->nChrom;
				pRegion->chStrand = pRefGene->chStrand;
				pRegion->nStart = npp1;
				pRegion->nEnd = npp2;
				strcpy(pRegion->strChrom, pRefGene->strChrom);

				GenomicRegionInsert(&pSelRegion, pRegion);
			}

			/* add TES */
			if( (nTESUP != 0) || (nTESDOWN != 0) )
			{
				npp1 = pRefGene->nTxEnd-nTESUP;
				npp2 = pRefGene->nTxEnd+nTESDOWN;
				if(npp1 < nRegionStart)
					npp1 = nRegionStart;
				if(npp2 > nRegionEnd)
					npp2 = nRegionEnd;

				pRegion = NULL;
				pRegion = GenomicRegionCreate();
				if(pRegion == NULL)
				{
					printf("Error: RefGeneSelectRegion, cannot create the selected region!\n");
					exit(EXIT_FAILURE);
				}
				pRegion->nChrom = pRefGene->nChrom;
				pRegion->chStrand = pRefGene->chStrand;
				pRegion->nStart = npp1;
				pRegion->nEnd = npp2;
				strcpy(pRegion->strChrom, pRefGene->strChrom);

				GenomicRegionInsert(&pSelRegion, pRegion);
			}

			/* add Exon & Intron */
			for(ni=0; ni<pRefGene->nExonCount; ni++)
			{
				/* exon */
				if(nIncludeExon == 1)
				{
					npp1 = IMGETAT(pRefGene->pmatExonStartsEnds, ni, 0);
					npp2 = IMGETAT(pRefGene->pmatExonStartsEnds, ni, 1);
					if(npp1 < nRegionStart)
						npp1 = nRegionStart;
					if(npp2 > nRegionEnd)
						npp2 = nRegionEnd;

					pRegion = NULL;
					pRegion = GenomicRegionCreate();
					if(pRegion == NULL)
					{
						printf("Error: RefGeneSelectRegion, cannot create the selected region!\n");
						exit(EXIT_FAILURE);
					}
					pRegion->nChrom = pRefGene->nChrom;
					pRegion->chStrand = pRefGene->chStrand;
					pRegion->nStart = npp1;
					pRegion->nEnd = npp2;
					strcpy(pRegion->strChrom, pRefGene->strChrom);

					GenomicRegionInsert(&pSelRegion, pRegion);
				}

				/* intron */
				if(nIncludeIntron == 1)
				{
					if( ni != (pRefGene->nExonCount-1) )
					{
						npp1 = IMGETAT(pRefGene->pmatExonStartsEnds, ni, 1)+1;
						npp2 = IMGETAT(pRefGene->pmatExonStartsEnds, (ni+1), 0)-1;
						if(npp1 < nRegionStart)
							npp1 = nRegionStart;
						if(npp2 > nRegionEnd)
							npp2 = nRegionEnd;

						pRegion = NULL;
						pRegion = GenomicRegionCreate();
						if(pRegion == NULL)
						{
							printf("Error: RefGeneSelectRegion, cannot create the selected region!\n");
							exit(EXIT_FAILURE);
						}
						pRegion->nChrom = pRefGene->nChrom;
						pRegion->chStrand = pRefGene->chStrand;
						pRegion->nStart = npp1;
						pRegion->nEnd = npp2;
						strcpy(pRegion->strChrom, pRefGene->strChrom);

						GenomicRegionInsert(&pSelRegion, pRegion);
					}
				}
			}
		}
		else if(pRefGene->chStrand == '-')
		{
			/* add TSS */
			if( (nTSSUP != 0) || (nTSSDOWN != 0) )
			{
				npp1 = pRefGene->nTxEnd-nTSSDOWN;
				npp2 = pRefGene->nTxEnd+nTSSUP;
				if(npp1 < nRegionStart)
					npp1 = nRegionStart;
				if(npp2 > nRegionEnd)
					npp2 = nRegionEnd;

				pRegion = NULL;
				pRegion = GenomicRegionCreate();
				if(pRegion == NULL)
				{
					printf("Error: RefGeneSelectRegion, cannot create the selected region!\n");
					exit(EXIT_FAILURE);
				}
				pRegion->nChrom = pRefGene->nChrom;
				pRegion->chStrand = pRefGene->chStrand;
				pRegion->nStart = npp1;
				pRegion->nEnd = npp2;
				strcpy(pRegion->strChrom, pRefGene->strChrom);

				GenomicRegionInsert(&pSelRegion, pRegion);
			}

			/* add TES */
			if( (nTESUP != 0) || (nTESDOWN != 0) )
			{
				npp1 = pRefGene->nTxStart-nTESDOWN;
				npp2 = pRefGene->nTxStart+nTESUP;
				if(npp1 < nRegionStart)
					npp1 = nRegionStart;
				if(npp2 > nRegionEnd)
					npp2 = nRegionEnd;

				pRegion = NULL;
				pRegion = GenomicRegionCreate();
				if(pRegion == NULL)
				{
					printf("Error: RefGeneSelectRegion, cannot create the selected region!\n");
					exit(EXIT_FAILURE);
				}
				pRegion->nChrom = pRefGene->nChrom;
				pRegion->chStrand = pRefGene->chStrand;
				pRegion->nStart = npp1;
				pRegion->nEnd = npp2;
				strcpy(pRegion->strChrom, pRefGene->strChrom);

				GenomicRegionInsert(&pSelRegion, pRegion);
			}

			/* add Exon & Intron */
			for(ni=(pRefGene->nExonCount-1); ni>=0; ni--)
			{
				/* exon */
				if(nIncludeExon == 1)
				{
					npp1 = IMGETAT(pRefGene->pmatExonStartsEnds, ni, 0);
					npp2 = IMGETAT(pRefGene->pmatExonStartsEnds, ni, 1);
					if(npp1 < nRegionStart)
						npp1 = nRegionStart;
					if(npp2 > nRegionEnd)
						npp2 = nRegionEnd;

					pRegion = NULL;
					pRegion = GenomicRegionCreate();
					if(pRegion == NULL)
					{
						printf("Error: RefGeneSelectRegion, cannot create the selected region!\n");
						exit(EXIT_FAILURE);
					}
					pRegion->nChrom = pRefGene->nChrom;
					pRegion->chStrand = pRefGene->chStrand;
					pRegion->nStart = npp1;
					pRegion->nEnd = npp2;
					strcpy(pRegion->strChrom, pRefGene->strChrom);

					GenomicRegionInsert(&pSelRegion, pRegion);
				}

				/* intron */
				if(nIncludeIntron == 1)
				{
					if( ni != 0 )
					{
						npp1 = IMGETAT(pRefGene->pmatExonStartsEnds, (ni-1), 1)+1;
						npp2 = IMGETAT(pRefGene->pmatExonStartsEnds, ni, 0)-1;
						if(npp1 < nRegionStart)
							npp1 = nRegionStart;
						if(npp2 > nRegionEnd)
							npp2 = nRegionEnd;

						pRegion = NULL;
						pRegion = GenomicRegionCreate();
						if(pRegion == NULL)
						{
							printf("Error: RefGeneSelectRegion, cannot create the selected region!\n");
							exit(EXIT_FAILURE);
						}
						pRegion->nChrom = pRefGene->nChrom;
						pRegion->chStrand = pRefGene->chStrand;
						pRegion->nStart = npp1;
						pRegion->nEnd = npp2;
						strcpy(pRegion->strChrom, pRefGene->strChrom);

						GenomicRegionInsert(&pSelRegion, pRegion);
					}
				}
			}
		}
		else
		{
			printf("Error: RefGeneSelectRegion, lack strand information!\n");
			exit(EXIT_FAILURE);
		}
	}

	GenomicRegionGetUnique(&pSelRegion);

	/* return */
	return pSelRegion;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGene_GetTargetTSSAround():                                          */
/*  get coordinates according to refgene structure.                        */
/* ----------------------------------------------------------------------- */ 
int RefGene_GetTargetTSSAround(char strRefGenePath[], char strTargetPath[],
			int nTSSUP, int nTSSDOWN, char strSpecies[], char strChrLen[],
			char strOutPath[])
{
	/* define */
	char strRefName[LINE_LENGTH];
	char strRefId[LINE_LENGTH];

	char strLine[LONG_LINE_LENGTH];
	char strRefLine[LONG_LINE_LENGTH];

	FILE *fpRefGene;
	FILE *fpOut;
	FILE *fpTarget;

	int nOK;
	int ni;
	int nChrLen;
	int nChrid;

	struct tagRefGene *pRefGene;
	int nPos1,nPos2;
	int nLastChr, nLastPos1, nLastPos2;
	struct INTMATRIX *pChrLen;


	/* open files */
	pChrLen = NULL;
	pChrLen = IMLOAD(strChrLen);
	if(pChrLen == NULL)
	{
		printf("Error: RefGene_GetTargetTSSAround, cannot load chromosome length!\n");
		exit(EXIT_FAILURE);
	}

	fpRefGene = NULL;
	fpRefGene = fopen(strRefGenePath, "r");
	fpOut = NULL;
	fpOut = fopen(strOutPath, "w");

	if( (fpRefGene == NULL) || (fpOut == NULL) )
	{
		printf("Error: RefGene_GetTargetTSSAround, cannot open files!\n");
		exit(EXIT_FAILURE);
	}

	/* get refgene */
	ni = 0;
	nLastChr = -1;
	nLastPos1 = -1;
	nLastPos2 = -1;
	while(fgets(strRefLine, LONG_LINE_LENGTH, fpRefGene) != NULL)
	{
		StrTrimLeft(strRefLine);
		StrTrimRight(strRefLine);

		if(strRefLine[0] == '\0')
			continue;

		sscanf(strRefLine, "%s ", strRefId);
		
		nOK = 0;
		
		fpTarget = NULL;
		fpTarget = fopen(strTargetPath, "r");
		if(fpTarget == NULL)
		{
			printf("Error: RefGene_GetTargetTSSAround, cannot open target file!\n");
			exit(EXIT_FAILURE);
		}
		while(fgets(strLine, LONG_LINE_LENGTH, fpTarget) != NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			
			if(strLine[0] == '\0')
				continue;

			sscanf(strLine, "%s ", strRefName);

			/* if found match */
			if(strcmp(strRefName, strRefId) == 0)
			{
				nOK = 1;
							
				pRefGene = NULL;
				pRefGene = RefGeneCreate();
				RefGeneInit_FromGenomeLabFormat(pRefGene, strRefLine, strSpecies);

				if(pRefGene->chStrand == '-')
				{
					nPos1 = pRefGene->nTxEnd-nTSSDOWN;
					nPos2 = pRefGene->nTxEnd+nTSSUP;
				}
				else
				{
					nPos1 = pRefGene->nTxStart-nTSSUP;
					nPos2 = pRefGene->nTxStart+nTSSDOWN;
				}
				if(nPos1 < 0)
					nPos1 = 0;
				if(nPos2 < 0)
					nPos2 = 0;
				
				nChrid = pRefGene->nChrom-1;
				nChrLen = pChrLen->pMatElement[nChrid];
				if(nPos1 >= nChrLen)
					nPos1 = nChrLen-1;
				if(nPos2 >= nChrLen)
					nPos2 = nChrLen-1;

				if(nPos1 > nPos2)
				{
					printf("Error: RefGene_GetTargetTSSAround, direction not correct (start>end?)!\n");
				}

				if( (nLastChr == pRefGene->nChrom) && (nLastPos1 == nPos1) && (nLastPos2 == nPos2))
				{
				}
				else
				{
					ni++;
					fprintf(fpOut, "%s_%d\t%s\t%d\t%d\t%c\n", 
						pRefGene->strName, ni, pRefGene->strChrom, nPos1, nPos2, pRefGene->chStrand);
					nLastPos1 = nPos1;
					nLastPos2 = nPos2;
					nLastChr = pRefGene->nChrom;
				}

				RefGeneDestroy(pRefGene);
			}			
		}
		fclose(fpTarget);
	}


	/* close files */
	fclose(fpRefGene);
	fclose(fpOut);

	DestroyIntMatrix(pChrLen);

	/* return  */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGene_GetOrtholog():                                                 */
/*  get ortholog genes according to colinear refgene structure & alignmetn */
/* ----------------------------------------------------------------------- */ 
int RefGene_GetOrtholog(char strRefGenePath[], 
			char strSourceSpecies[], char strSourceRefGenePath[],
			char strDestSpecies[], char strMapRefGenePath[],
			char strDestRefGenePath[], 
			char strOutPath[])
{
	/* define */

	/* for database */
	FILE *fpSourceRefGene;
	struct tagRefGene *pSourceRefGeneList;
	struct tagRefGene **vSourceRefGene;
	int nSourceRefNum;

	FILE *fpMapRefGene;
	struct tagRefGene *pMapRefGeneList;
	struct tagRefGene **vMapRefGene;
	int nMapRefNum;

	int nDestMap;
	FILE *fpDestRefGene;
	struct tagRefGene *pDestRefGeneList;
	struct tagRefGene **vDestRefGene;
	int nDestRefNum;

	struct tagRefGene *pRefGene, *pCurrentRefGene;
	char strRefLine[LONG_LINE_LENGTH];
	int ni,nj,nk,nflank;
	int nx,ny,nz,nlastid;

	/* for target search */
	FILE *fpIn;
	FILE *fpMapOut;
	FILE *fpDestOut;
	
	int vRefSourceIdx[ORTHOLOG_COLINEAR_NUM];
	char vRefColinear[ORTHOLOG_COLINEAR_NUM][LINE_LENGTH];

	char strLine[LINE_LENGTH];
	char strRefName[LINE_LENGTH];

	FILE *fpTemp;
	FILE *fpTemp1,*fpTemp2;
	int nSourceHitNum,nMapHitNum;
	struct INTMATRIX *pHitLimit;
	char *chp,*chp2;

	int nPos1,nPos2,nPos3;
	int nMatchNum;
	int nOrientationNum;
	int nMatched;


	/* init */
	if(strcmp(strDestRefGenePath, "NULL") == 0)
	{
		nDestMap = 0;
	}
	else
	{
		nDestMap = 1;
	}

	pSourceRefGeneList = NULL;
	pMapRefGeneList = NULL;
	pDestRefGeneList = NULL;
	nSourceRefNum = 0;
	nMapRefNum = 0;
	nDestRefNum = 0;
	vSourceRefGene = NULL;
	vMapRefGene = NULL;
	vDestRefGene = NULL;

	/* load source refgene */
	fpSourceRefGene = NULL;
	fpSourceRefGene = fopen(strSourceRefGenePath, "r");
	if(fpSourceRefGene == NULL)
	{
		printf("Error: RefGene_GetOrtholog, cannot open source refgene file!\n");
		exit(EXIT_FAILURE);
	}

	pCurrentRefGene = NULL;
	while(fgets(strRefLine, LONG_LINE_LENGTH, fpSourceRefGene) != NULL)
	{
		StrTrimLeft(strRefLine);
		StrTrimRight(strRefLine);
		if(strRefLine[0] == '\0')
			continue;

		pRefGene = NULL;
		pRefGene = RefGeneCreate();
		RefGeneInit_FromGenomeLabFormat(pRefGene, strRefLine, strSourceSpecies);
		if(pSourceRefGeneList == NULL)
		{
			pSourceRefGeneList = pRefGene;
			pCurrentRefGene = pRefGene;
		}
		else
		{
			pCurrentRefGene->pNext = pRefGene;
			pCurrentRefGene = pRefGene;
		}
		nSourceRefNum++;
	}

	fclose(fpSourceRefGene);

	vSourceRefGene = (struct tagRefGene **)calloc(nSourceRefNum, sizeof(struct tagRefGene*));
	if(vSourceRefGene == NULL)
	{
		printf("Error: RefGene_GetOrtholog, organize refgene into array!\n");
		exit(EXIT_FAILURE);
	}

	ni=0;
	while(pSourceRefGeneList != NULL)
	{
		pRefGene = pSourceRefGeneList;
		pSourceRefGeneList = pRefGene->pNext;
		pRefGene->pNext = NULL;
		vSourceRefGene[ni] = pRefGene;
		ni++;
	}

	if(ni != nSourceRefNum)
	{
		printf("Error: RefGene_GetOrtholog, refgene number not match!\n");
		exit(EXIT_FAILURE);
	}


	/* load map refgene */
	fpMapRefGene = NULL;
	fpMapRefGene = fopen(strMapRefGenePath, "r");
	if(fpMapRefGene == NULL)
	{
		printf("Error: RefGene_GetOrtholog, cannot open map refgene file!\n");
		exit(EXIT_FAILURE);
	}

	pCurrentRefGene = NULL;
	while(fgets(strRefLine, LONG_LINE_LENGTH, fpMapRefGene) != NULL)
	{
		StrTrimLeft(strRefLine);
		StrTrimRight(strRefLine);
		if(strRefLine[0] == '\0')
			continue;

		pRefGene = NULL;
		pRefGene = RefGeneCreate();
		RefGeneInit_FromGenomeLabFormat(pRefGene, strRefLine, strDestSpecies);
		if(pMapRefGeneList == NULL)
		{
			pMapRefGeneList = pRefGene;
			pCurrentRefGene = pRefGene;
		}
		else
		{
			pCurrentRefGene->pNext = pRefGene;
			pCurrentRefGene = pRefGene;
		}
		nMapRefNum++;
	}

	fclose(fpMapRefGene);

	vMapRefGene = (struct tagRefGene **)calloc(nMapRefNum, sizeof(struct tagRefGene*));
	if(vMapRefGene == NULL)
	{
		printf("Error: RefGene_GetOrtholog, organize refgene into array!\n");
		exit(EXIT_FAILURE);
	}

	ni=0;
	while(pMapRefGeneList != NULL)
	{
		pRefGene = pMapRefGeneList;
		pMapRefGeneList = pRefGene->pNext;
		pRefGene->pNext = NULL;
		vMapRefGene[ni] = pRefGene;
		ni++;
	}

	if(ni != nMapRefNum)
	{
		printf("Error: RefGene_GetOrtholog, refgene number not match!\n");
		exit(EXIT_FAILURE);
	}


	/* load dest refgene */
	if(nDestMap == 1)
	{
		fpDestRefGene = NULL;
		fpDestRefGene = fopen(strDestRefGenePath, "r");
		if(fpDestRefGene == NULL)
		{
			printf("Error: RefGene_GetOrtholog, cannot open dest refgene file!\n");
			exit(EXIT_FAILURE);
		}

		pCurrentRefGene = NULL;
		while(fgets(strRefLine, LONG_LINE_LENGTH, fpDestRefGene) != NULL)
		{
			StrTrimLeft(strRefLine);
			StrTrimRight(strRefLine);
			if(strRefLine[0] == '\0')
				continue;

			pRefGene = NULL;
			pRefGene = RefGeneCreate();
			RefGeneInit_FromGenomeLabFormat(pRefGene, strRefLine, strDestSpecies);
			if(pDestRefGeneList == NULL)
			{
				pDestRefGeneList = pRefGene;
				pCurrentRefGene = pRefGene;
			}
			else
			{
				pCurrentRefGene->pNext = pRefGene;
				pCurrentRefGene = pRefGene;
			}
			nDestRefNum++;
		}

		fclose(fpDestRefGene);

		vDestRefGene = (struct tagRefGene **)calloc(nDestRefNum, sizeof(struct tagRefGene*));
		if(vDestRefGene == NULL)
		{
			printf("Error: RefGene_GetOrtholog, organize refgene into array!\n");
			exit(EXIT_FAILURE);
		}

		ni=0;
		while(pDestRefGeneList != NULL)
		{
			pRefGene = pDestRefGeneList;
			pDestRefGeneList = pRefGene->pNext;
			pRefGene->pNext = NULL;
			vDestRefGene[ni] = pRefGene;
			ni++;
		}

		if(ni != nDestRefNum)
		{
			printf("Error: RefGene_GetOrtholog, refgene number not match!\n");
			exit(EXIT_FAILURE);
		}
	}


	/* locate ortholog */
	fpIn = NULL;
	fpIn = fopen(strRefGenePath, "r");
	if(fpIn == NULL)
	{
		printf("Error: RefGene_GetOrtholog, cannot open target refgene file!\n");
		exit(EXIT_FAILURE);
	}

	sprintf(strLine, "%s.pair", strOutPath);
	fpMapOut = NULL;
	fpMapOut = fopen(strLine, "w");
	if(fpMapOut == NULL)
	{
		printf("Error: RefGene_GetOrtholog, cannot open map output file!\n");
		exit(EXIT_FAILURE);
	}

	if(nDestMap == 1)
	{
		sprintf(strLine, "%s.map", strOutPath);
		fpDestOut = NULL;
		fpDestOut = fopen(strLine, "w");
		if(fpDestOut == NULL)
		{
			printf("Error: RefGene_GetOrtholog, cannot open dest output file!\n");
			exit(EXIT_FAILURE);
		}
	}


	/* process target refgene one by one */
	fpTemp1 = NULL;
	fpTemp1 = fopen("sourcehits2.tmp", "w");
	if(fpTemp1 == NULL)
	{
		printf("Error: cannot create temp file!\n");
		exit(EXIT_FAILURE);
	}
	fpTemp2 = NULL;
	fpTemp2 = fopen("maphits2.tmp", "w");
	if(fpTemp2 == NULL)
	{
		printf("Error: cannot create temp file!\n");
		exit(EXIT_FAILURE);
	}
	while(fgets(strRefLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strRefLine);
		StrTrimRight(strRefLine);
		if(strRefLine[0] == '\0')
			continue;

		sscanf(strRefLine, "%s", strRefName);

		/* look up source database */
		nSourceHitNum = 0;
		nflank = (ORTHOLOG_COLINEAR_NUM-1)/2;
		fpTemp = NULL;
		fpTemp = fopen("sourcehits.tmp", "w");
		if(fpTemp == NULL)
		{
			printf("Error: cannot create temp file!\n");
			exit(EXIT_FAILURE);
		}

		for(ni=0; ni<nSourceRefNum; ni++)
		{
			/* if found the location in the source organism */
			if(strcmp(strRefName, vSourceRefGene[ni]->strName) == 0)
			{
				nSourceHitNum++;
				/* init flank */
				for(nk=0; nk<ORTHOLOG_COLINEAR_NUM; nk++)
				{
					vRefSourceIdx[nk] = -1;
					strcpy(vRefColinear[nk], "---");
				}

				vRefSourceIdx[nflank] = ni;
				strcpy(vRefColinear[nflank], strRefName);

				/* get upstream flank */
				nk = 0;
				nlastid = ni;
				for(nj=ni-1; nj>=0; nj--)
				{
					if(vSourceRefGene[ni]->nChrom == vSourceRefGene[nj]->nChrom)
					{
						if((vSourceRefGene[ni]->nTxStart-vSourceRefGene[nj]->nTxEnd) <= ORTHOLOG_COLINEAR_WIN)
						{
							if(vSourceRefGene[nj]->nTxEnd < vSourceRefGene[nlastid]->nTxStart)
							{
								nk++;
								vRefSourceIdx[nflank-nk] = nj;
								nlastid = nj;
								strcpy(vRefColinear[nflank-nk], vSourceRefGene[nj]->strName);
							}

							if(nk == nflank)
								break;
						}
						else
						{
							break;
						}

					}
					else
					{
						break;
					}
				}

				/* get downstream flank */
				nk = 0;
				nlastid = ni;
				for(nj=ni+1; nj<nSourceRefNum; nj++)
				{
					if(vSourceRefGene[ni]->nChrom == vSourceRefGene[nj]->nChrom)
					{
						if((vSourceRefGene[nj]->nTxStart-vSourceRefGene[ni]->nTxEnd) <= ORTHOLOG_COLINEAR_WIN)
						{
							if( vSourceRefGene[nlastid]->nTxEnd < vSourceRefGene[nj]->nTxStart)
							{
								nk++;
								vRefSourceIdx[nflank+nk] = nj;
								nlastid = nj;
								strcpy(vRefColinear[nflank+nk], vSourceRefGene[nj]->strName);
							}

							if(nk == nflank)
								break;
						}
						else
						{
							break;
						}
					}
					else
					{
						break;
					}
				}
				
				fprintf(fpTemp, "%d\t", ni);
				fprintf(fpTemp1, "%d\t", ni);
				
				/* remove redundancy */
				for(nk=0; nk<ORTHOLOG_COLINEAR_NUM; nk++)
				{
					if(vRefSourceIdx[nk] != -1)
					{
						for(nx=0; nx<nk; nx++)
						{
							if(vRefSourceIdx[nx] != -1)
							{
								if(strcmp(vRefColinear[nk], vRefColinear[nx]) == 0)
								{
									strcpy(vRefColinear[nk], "---");
									vRefSourceIdx[nk] = -1;
									break;
								}
							}
						}
					}
					fprintf(fpTemp, "%s\t", vRefColinear[nk]);
					fprintf(fpTemp1, "%s\t", vRefColinear[nk]);
				}
				fprintf(fpTemp, "\n");
				fprintf(fpTemp1, "\n");
			}
		}

		fclose(fpTemp);


		/* look up map database */
		nMapHitNum = 0;
		fpTemp = NULL;
		fpTemp = fopen("maphits.tmp", "w");
		if(fpTemp == NULL)
		{
			printf("Error: cannot create temp file!\n");
			exit(EXIT_FAILURE);
		}

		for(ni=0; ni<nMapRefNum; ni++)
		{
			/* if found the location in the map organism */
			if(strcmp(strRefName, vMapRefGene[ni]->strName) == 0)
			{
				nMapHitNum++;

				/* get upstream flank */
				nk = 0;
				nlastid = ni;
				for(nj=ni-1; nj>=0; nj--)
				{
					if(vMapRefGene[ni]->nChrom == vMapRefGene[nj]->nChrom)
					{
						if( (vMapRefGene[ni]->nTxStart-vMapRefGene[nj]->nTxEnd) <= ORTHOLOG_COLINEAR_WIN )
						{
							if(vMapRefGene[nj]->nTxEnd < vMapRefGene[nlastid]->nTxStart)
							{
								nk++;
								if(nk > ORTHOLOG_MAP_COLINEAR_NUM)
								{
									nlastid = nj+1;
									break;
								}
								nlastid = nj;
							}
						}
						else
						{
							nlastid = nj+1;
							break;
						}
					}
					else
					{
						nlastid = nj+1;
						break;
					}
				}

				fprintf(fpTemp, "%d\t", nlastid);
				fprintf(fpTemp, "%d\t", ni);

				fprintf(fpTemp2, "%d\t", nlastid);
				fprintf(fpTemp2, "%d\t", ni);

				/* get downstream flank */
				nk = 0;
				nlastid = ni;
				for(nj=ni+1; nj<nMapRefNum; nj++)
				{
					if(vMapRefGene[ni]->nChrom == vMapRefGene[nj]->nChrom)
					{
						if( (vMapRefGene[nj]->nTxStart-vMapRefGene[ni]->nTxEnd) <= ORTHOLOG_COLINEAR_WIN )
						{
							if(vMapRefGene[nlastid]->nTxEnd < vMapRefGene[nj]->nTxStart)
							{
								nk++;
								if(nk > ORTHOLOG_MAP_COLINEAR_NUM)
								{
									nlastid = nj-1;
									break;
								}
								nlastid = nj;
							}
						}
						else
						{
							nlastid = nj-1;
							break;
						}
					}
					else
					{
						nlastid = nj-1;
						break;
					}
				}

				fprintf(fpTemp, "%d\n", nlastid);
				fprintf(fpTemp2, "%d\n", nlastid);
			}
		}

		fclose(fpTemp);
		

		/* verify collinearity */
		if(nSourceHitNum == 0)
		{
		}
		else if(nMapHitNum == 0)
		{
			fprintf(fpMapOut, ">%s\n", strRefName);
			fpTemp = NULL;
			fpTemp = fopen("sourcehits.tmp", "r");
			if(fpTemp == NULL)
			{
				printf("Error: sourcehits.tmp does not exist!\n");
				exit(EXIT_FAILURE);
			}

			nx = 0;
			while(fgets(strRefLine, LONG_LINE_LENGTH, fpTemp) != NULL)
			{
				StrTrimLeft(strRefLine);
				StrTrimRight(strRefLine);
				if(strRefLine[0] == '\0')
					continue;

				fprintf(fpMapOut, "#%d\t0\n", nx);
				sscanf(strRefLine, "%d", &ny);
				fprintf(fpMapOut, "%s\t0\t0\t", strSourceSpecies);
				RefGeneWrite(vSourceRefGene[ny], fpMapOut);
				/* fprintf(fpMapOut, "%s\t0\t0\t---\n", strDestSpecies); */
				fprintf(fpMapOut, "\n");

				nx++;
			}

			fclose(fpTemp);
			if(nx != nSourceHitNum)
			{
				printf("Error: source hit number not match!\n");
				exit(EXIT_FAILURE);
			}

		}
		else
		{
			pHitLimit = NULL;
			pHitLimit = IMLOAD("maphits.tmp");
			if(pHitLimit == NULL)
			{
				printf("Error: cannot reload maphits.tmp!\n");
				exit(EXIT_FAILURE);
			}
			if(pHitLimit->nHeight != nMapHitNum)
			{
				printf("Error: map hit number not match!\n");
				exit(EXIT_FAILURE);
			}

			fprintf(fpMapOut, ">%s\n", strRefName);
			fpTemp = NULL;
			fpTemp = fopen("sourcehits.tmp", "r");
			if(fpTemp == NULL)
			{
				printf("Error: sourcehits.tmp does not exist!\n");
				exit(EXIT_FAILURE);
			}

			nx = 0;
			while(fgets(strRefLine, LONG_LINE_LENGTH, fpTemp) != NULL)
			{
				StrTrimLeft(strRefLine);
				StrTrimRight(strRefLine);
				if(strRefLine[0] == '\0')
					continue;

				fprintf(fpMapOut, "#%d\t%d\n", nx, nMapHitNum);
				sscanf(strRefLine, "%d", &ny);
				fprintf(fpMapOut, "%s\t0\t0\t", strSourceSpecies);
				RefGeneWrite(vSourceRefGene[ny], fpMapOut);
				
				chp = strchr(strRefLine, '\t');
				chp++;

				/* load string */
				for(nz=0; nz<ORTHOLOG_COLINEAR_NUM; nz++)
				{
					chp2 = strchr(chp, '\t');
					if(chp2 != NULL)
					{
						*chp2 = '\0';
						strcpy(vRefColinear[nz], chp);
						chp = chp2+1;
					}
					else
					{
						strcpy(vRefColinear[nz], chp);
						if( (nz+1) != ORTHOLOG_COLINEAR_NUM)
						{
							printf("Error: source hit number not match!\n");
							exit(EXIT_FAILURE);
						}
					}
				}

				/* match */
				for(nz=0; nz<nMapHitNum; nz++)
				{
					nMatchNum = 0;
					nOrientationNum = 0;
					nPos1 = IMGETAT(pHitLimit, nz, 0);
					nPos2 = IMGETAT(pHitLimit, nz, 1);
					nPos3 = IMGETAT(pHitLimit, nz, 2);

					for(nk=0; nk<ORTHOLOG_COLINEAR_NUM; nk++)
					{
						nMatched = 0;
						vRefSourceIdx[nk] = 0;
						for(ny=nPos1; ny<nPos2; ny++)
						{
							if(strcmp(vRefColinear[nk], vMapRefGene[ny]->strName) == 0)
							{
								nMatched = 1;
								nMatchNum++;
								vRefSourceIdx[nk] = 1;
								break;
							}
						}
						if(nMatched == 1)
							continue;

						if(strcmp(vRefColinear[nk], vMapRefGene[ny]->strName) == 0)
						{
							nMatchNum++;
							continue;
						}

						ny++;

						for(; ny<=nPos3; ny++)
						{
							if(strcmp(vRefColinear[nk], vMapRefGene[ny]->strName) == 0)
							{
								nMatched = 1;
								nMatchNum++;
								vRefSourceIdx[nk] = -1;
								break;
							}
						}
					}

					nflank = (ORTHOLOG_COLINEAR_NUM-1)/2;
					for(nk=0; nk<nflank; nk++)
					{
						for(nj=0; nj<nflank; nj++)
						{
							if(nk != nj)
							{
								if(vRefSourceIdx[nk]*vRefSourceIdx[nj] > 0)
								{
									nOrientationNum++;
								}
							}
						}
						nj++;
						for(; nj<ORTHOLOG_COLINEAR_NUM; nj++)
						{
							if(vRefSourceIdx[nk]*vRefSourceIdx[nj] < 0)
							{
								nOrientationNum++;
							}
						}
					}

					nk++;
					for(; nk<ORTHOLOG_COLINEAR_NUM; nk++)
					{
						for(nj=0; nj<nflank; nj++)
						{
							if(vRefSourceIdx[nk]*vRefSourceIdx[nj] < 0)
							{
								nOrientationNum++;
							}
						}
						nj++;
						for(; nj<ORTHOLOG_COLINEAR_NUM; nj++)
						{
							if(nj != nk)
							{
								if(vRefSourceIdx[nk]*vRefSourceIdx[nj] > 0)
								{
									nOrientationNum++;
								}
							}
						}
					}

					fprintf(fpMapOut, "%s\t%d\t%d\t", strDestSpecies, nMatchNum, nOrientationNum);
					RefGeneWrite(vMapRefGene[nPos2], fpMapOut);
				}

				fprintf(fpMapOut, "\n");
				nx++;
			}

			fclose(fpTemp);
			if(nx != nSourceHitNum)
			{
				printf("Error: source hit number not match!\n");
				exit(EXIT_FAILURE);
			}

			DestroyIntMatrix(pHitLimit);
		}
	}

	fclose(fpTemp1);
	fclose(fpTemp2);

	/* delete temp files */
	if(strcmp(OS_SYSTEM, "UNIX") == 0)
	{
		sprintf(strLine, "rm sourcehits.tmp");
		system(strLine);
		sprintf(strLine, "rm maphits.tmp");
		system(strLine);
		sprintf(strLine, "rm sourcehits2.tmp");
		system(strLine);
		sprintf(strLine, "rm maphits2.tmp");
		system(strLine);
	}
	else
	{
		sprintf(strLine, "del sourcehits.tmp");
		system(strLine);
		sprintf(strLine, "del maphits.tmp");
		system(strLine);
		sprintf(strLine, "del sourcehits2.tmp");
		system(strLine);
		sprintf(strLine, "del maphits2.tmp");
		system(strLine);
	}	

	/* close files */
	fclose(fpIn);
	fclose(fpMapOut);
	if(nDestMap == 1)
	{
		fclose(fpDestOut);
	}


	/* release memory */
	for(ni=0; ni<nSourceRefNum; ni++)
	{
		RefGeneDestroy(vSourceRefGene[ni]);
		vSourceRefGene[ni] = NULL;
	}
	free(vSourceRefGene);

	for(ni=0; ni<nMapRefNum; ni++)
	{
		RefGeneDestroy(vMapRefGene[ni]);
		vMapRefGene[ni] = NULL;
	}
	free(vMapRefGene);

	if(nDestMap == 1)
	{
		for(ni=0; ni<nDestRefNum; ni++)
		{
			RefGeneDestroy(vDestRefGene[ni]);
			vDestRefGene[ni] = NULL;
		}
		free(vDestRefGene);
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGene_MergeMultiOrtholog():                                          */
/*  Merge ortholog genes from multiple pairwise ortholog mapping.          */
/*  Return the number of ortholog groups.                                  */
/* ----------------------------------------------------------------------- */ 
int RefGene_MergeMultiOrtholog(int nSpeciesNum, struct tagString **vSpeciesName,
							   struct tagString **vOrthologPath, char strOutPath[])
{
	/* define */
	struct DOUBLEMATRIX *pMatchCount;
	struct DOUBLEMATRIX *pOrienCount;
	struct tagRefGene **vRefGene;
	struct tagRefGene *pRefGene;
	double dMatchSum,dOrienSum;
	int nMatchCount, nOrienCount;

	struct DOUBLEMATRIX *pOptMatchCount;
	struct DOUBLEMATRIX *pOptOrienCount;
	struct tagRefGene **vOptRefGene;
	double dOptMatchSum,dOptOrienSum;
	int nInitMap;
	int nOrthoGroupNum;

	/* file */
	FILE **vfpIn;
	FILE *fpOut;
	int ni,nj;
	char strRefLine0[LONG_LINE_LENGTH];
	char strRefLine1[LONG_LINE_LENGTH];
	char strRefName0[LINE_LENGTH];
	char strRefName1[LINE_LENGTH];

	char strLine[LINE_LENGTH];
	char strRefLine[LONG_LINE_LENGTH];
	char strSpecies[LINE_LENGTH];

	int nSourceId,nDestNum;
	char *chp,*chp2;
	int nAddNum,nSwitchon;

	/* init */
	pOptMatchCount = NULL;
	pOptOrienCount = NULL;
	vOptRefGene = NULL;
	dOptMatchSum = 0.0;
	dOptOrienSum = 0.0;
	nInitMap = 0;
	nOrthoGroupNum = 0;


	/* open files */
	if(nSpeciesNum <= 1)
	{
		printf("Warning: RefGene_MergeMultiOrtholog, no files to merge!\n");
		return PROC_FAILURE;
	}

	vfpIn = NULL;
	vfpIn = (FILE **)calloc(nSpeciesNum, sizeof(FILE *));
	if(vfpIn == NULL)
	{
		printf("Error: RefGene_MergeMultiOrtholog, cannot match ortholog!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nSpeciesNum; ni++)
	{
		sprintf(strLine, "%s.pair", vOrthologPath[ni]->m_pString);
		vfpIn[ni] = fopen(strLine, "r");
		if(vfpIn[ni] == NULL)
		{
			printf("Error: RefGene_MergeMultiOrtholog, cannot open file!\n");
			exit(EXIT_FAILURE);
		}
	}

	fpOut = NULL;
	fpOut = fopen(strOutPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: RefGene_MergeMultiOrtholog, cannot open file!\n");
		exit(EXIT_FAILURE);
	}

	/* get the best ortholog */
	while(fgets(strRefLine0, LONG_LINE_LENGTH, vfpIn[0]) != NULL)
	{
		StrTrimLeft(strRefLine0);
		StrTrimRight(strRefLine0);
		if(strRefLine0[0] == '\0')
			continue;

		if(strRefLine0[0] == '>')
		{
			/* write last gene */
			if(nInitMap == 1)
			{
				fprintf(fpOut, ">%s\n", strRefName0);
				for(ni=0; ni<nSpeciesNum; ni++)
				{
					if((pOptMatchCount != NULL) && (pOptOrienCount != NULL))
					{
						fprintf(fpOut, "%s\t%d\t%d\t", vSpeciesName[ni]->m_pString, (int)(pOptMatchCount->pMatElement[ni]),
							(int)(pOptOrienCount->pMatElement[ni]));
					}
					else
					{
						fprintf(fpOut, "%s\t0\t0\t", vSpeciesName[ni]->m_pString);
					}
					if(vOptRefGene[ni] != NULL)
					{
						RefGeneWrite(vOptRefGene[ni], fpOut);
					}
					else
					{
						fprintf(fpOut, "---\n");
					}
				}
				fprintf(fpOut, "\n");
				nOrthoGroupNum++;

				DestroyDoubleMatrix(pOptMatchCount);
				DestroyDoubleMatrix(pOptOrienCount);
				for(ni=0; ni<nSpeciesNum; ni++)
				{
					RefGeneDestroy(vOptRefGene[ni]);
					vOptRefGene[ni] = NULL;
				}
				free(vOptRefGene);
			}

			/* init map */
			nInitMap = 1;
			nAddNum = 0;
			nSwitchon = 0;
			sprintf(strRefName0, "%s", strRefLine0+1);

			/* synchronize all files */
			for(ni=1; ni<nSpeciesNum; ni++)
			{
				while(fgets(strRefLine1, LONG_LINE_LENGTH, vfpIn[ni]) != NULL)
				{
					StrTrimLeft(strRefLine1);
					StrTrimRight(strRefLine1);
					if(strRefLine1[0] == '\0')
						continue;

					if(strRefLine1[0] == '>')
					{
						sprintf(strRefName1, "%s", strRefLine1+1);
						break;
					}
				}

				if(strcmp(strRefName0, strRefName1) != 0)
				{
					printf("Error: RefGene_MergeMultiOrtholog, reference refgene id not match!\n");
					exit(EXIT_FAILURE);
				}
			}

			/* process one by one */
			pOptMatchCount = NULL;
			pOptOrienCount = NULL;
			vOptRefGene = NULL;
			vOptRefGene = (struct tagRefGene **)calloc(nSpeciesNum, sizeof(struct tagRefGene*));
			if(vOptRefGene == NULL)
			{
				printf("Error: RefGene_MergeMultiOrtholog, cannot create refgene vector!\n");
				exit(EXIT_FAILURE);
			}
			dOptMatchSum = 0.0; 
			dOptOrienSum = 0.0;				
		}
		else if(strRefLine0[0] == '#')
		{
			/* init */
			pMatchCount = NULL;
			pMatchCount = CreateDoubleMatrix(1, nSpeciesNum);
			if(pMatchCount == NULL)
			{
				printf("Error: RefGene_MergeMultiOrtholog, cannot create matchcount matrix!\n");
				exit(EXIT_FAILURE);
			}
			pOrienCount = NULL;
			pOrienCount = CreateDoubleMatrix(1, nSpeciesNum);
			if(pOrienCount == NULL)
			{
				printf("Error: RefGene_MergeMultiOrtholog, cannot create oriencount matrix!\n");
				exit(EXIT_FAILURE);
			}
			vRefGene = NULL;
			vRefGene = (struct tagRefGene **)calloc(nSpeciesNum, sizeof(struct tagRefGene*));
			if(vRefGene == NULL)
			{
				printf("Error: RefGene_MergeMultiOrtholog, cannot create refgene vector!\n");
				exit(EXIT_FAILURE);
			}
			dMatchSum = 0.0;
			dOrienSum = 0.0;

			/* process the first file */
			sprintf(strLine, "%s", strRefLine0+1);
			sscanf(strLine, "%d %d", &nSourceId, &nDestNum);

			fgets(strRefLine0, LONG_LINE_LENGTH, vfpIn[0]);
			StrTrimLeft(strRefLine0);
			StrTrimRight(strRefLine0);

			chp = strRefLine0;
			chp2 = strchr(chp, '\t');
			*chp2 = '\0';
			strcpy(strSpecies, chp);
			chp = chp2+1;
			chp2 = strchr(chp, '\t');
			*chp2 = '\0';
			nMatchCount = atoi(chp);
			chp = chp2+1;
			chp2 = strchr(chp, '\t');
			*chp2 = '\0';
			nOrienCount = atoi(chp);
			chp = chp2+1;
			strcpy(strRefLine, chp);

			if(strcmp(strSpecies, vSpeciesName[0]->m_pString) != 0)
			{
				printf("Error: RefGene_MergeMultiOrtholog, species not match!\n");
				exit(EXIT_FAILURE);
			}
			pMatchCount->pMatElement[0] = nMatchCount;
			pOrienCount->pMatElement[0] = nOrienCount;
			pRefGene = NULL;
			pRefGene = RefGeneCreate();
			RefGeneInit_FromGenomeLabFormat(pRefGene, strRefLine, strSpecies);
			if(strcmp(strRefName0, pRefGene->strName) != 0)
			{
				printf("Error: RefGene_MergeMultiOrtholog, refgene not match!\n");
				exit(EXIT_FAILURE);
			}
			vRefGene[0] = pRefGene;

			for(nj=0; nj<nDestNum; nj++)
			{
				fgets(strRefLine0, LONG_LINE_LENGTH, vfpIn[0]);
				StrTrimLeft(strRefLine0);
				StrTrimRight(strRefLine0);
				
				chp = strRefLine0;
				chp2 = strchr(chp, '\t');
				*chp2 = '\0';
				strcpy(strSpecies, chp);
				chp = chp2+1;
				chp2 = strchr(chp, '\t');
				*chp2 = '\0';
				nMatchCount = atoi(chp);
				chp = chp2+1;
				chp2 = strchr(chp, '\t');
				*chp2 = '\0';
				nOrienCount = atoi(chp);
				chp = chp2+1;
				strcpy(strRefLine, chp);

				if(strcmp(strSpecies, vSpeciesName[0]->m_pString) != 0)
				{
					printf("Error: RefGene_MergeMultiOrtholog, species not match!\n");
					exit(EXIT_FAILURE);
				}
				if((nMatchCount+ORTHOLOG_COMBCOUNT_COEF*nOrienCount) >= (pMatchCount->pMatElement[0]+ORTHOLOG_COMBCOUNT_COEF*(pOrienCount->pMatElement[0])))
				{
					pMatchCount->pMatElement[0] = nMatchCount;
					pOrienCount->pMatElement[0] = nOrienCount;
					/* pRefGene = NULL;
					pRefGene = RefGeneCreate();
					RefGeneInit_FromGenomeLabFormat(pRefGene, strRefLine, strSpecies); */
				}
			}

			/* process remaining files */
			for(ni=1; ni<nSpeciesNum; ni++)
			{
				while(fgets(strRefLine1, LONG_LINE_LENGTH, vfpIn[ni]) != NULL)
				{
					StrTrimLeft(strRefLine1);
					StrTrimRight(strRefLine1);
					if(strRefLine1[0] == '\0')
						continue;

					if(strRefLine1[0] == '#')
						break;
				}
				
				sprintf(strLine, "%s", strRefLine1+1);
				sscanf(strLine, "%d %d", &nSourceId, &nDestNum);

				fgets(strRefLine1, LONG_LINE_LENGTH, vfpIn[ni]);
				StrTrimLeft(strRefLine1);
				StrTrimRight(strRefLine1);
				chp = strRefLine1;
				chp2 = strchr(chp, '\t');
				*chp2 = '\0';
				strcpy(strSpecies, chp);
				chp = chp2+1;
				chp2 = strchr(chp, '\t');
				*chp2 = '\0';
				nMatchCount = atoi(chp);
				chp = chp2+1;
				chp2 = strchr(chp, '\t');
				*chp2 = '\0';
				nOrienCount = atoi(chp);
				chp = chp2+1;
				strcpy(strRefLine, chp);
				if(strcmp(strSpecies, vSpeciesName[0]->m_pString) != 0)
				{
					printf("Error: RefGene_MergeMultiOrtholog, species not match!\n");
					exit(EXIT_FAILURE);
				}

				pRefGene = NULL;
				pRefGene = RefGeneCreate();
				RefGeneInit_FromGenomeLabFormat(pRefGene, strRefLine, strSpecies);
				if(strcmp(strRefName0, pRefGene->strName) != 0)
				{
					printf("Error: RefGene_MergeMultiOrtholog, refgene not match!\n");
					exit(EXIT_FAILURE);
				}
				RefGeneDestroy(pRefGene);
				
				for(nj=0; nj<nDestNum; nj++)
				{
					fgets(strRefLine1, LONG_LINE_LENGTH, vfpIn[ni]);
					StrTrimLeft(strRefLine1);
					StrTrimRight(strRefLine1);
					chp = strRefLine1;
					chp2 = strchr(chp, '\t');
					*chp2 = '\0';
					strcpy(strSpecies, chp);
					chp = chp2+1;
					chp2 = strchr(chp, '\t');
					*chp2 = '\0';
					nMatchCount = atoi(chp);
					chp = chp2+1;
					chp2 = strchr(chp, '\t');
					*chp2 = '\0';
					nOrienCount = atoi(chp);
					chp = chp2+1;
					strcpy(strRefLine, chp);
					if(strcmp(strSpecies, vSpeciesName[ni]->m_pString) != 0)
					{
						printf("Error: RefGene_MergeMultiOrtholog, species not match!\n");
						exit(EXIT_FAILURE);
					}
					if((nMatchCount+ORTHOLOG_COMBCOUNT_COEF*nOrienCount) >= (pMatchCount->pMatElement[ni]+ORTHOLOG_COMBCOUNT_COEF*pOrienCount->pMatElement[ni]))
					{
						pMatchCount->pMatElement[ni] = nMatchCount;
						pOrienCount->pMatElement[ni] = nOrienCount;
						pRefGene = NULL;
						pRefGene = RefGeneCreate();
						RefGeneInit_FromGenomeLabFormat(pRefGene, strRefLine, strSpecies);
						if(vRefGene[ni] == NULL)
						{
							vRefGene[ni] = pRefGene;
						}
						else
						{
							RefGeneDestroy(vRefGene[ni]);
							vRefGene[ni] = pRefGene;
						}
					}
				}

			}

			
			/* compute total score */
			dMatchSum = 0.0;
			dOrienSum = 0.0;
			for(ni=1; ni<nSpeciesNum; ni++)
			{
				if( (pMatchCount->pMatElement[ni] >= ORTHOLOG_MATCHCOUNT_LIM) &&
					(pOrienCount->pMatElement[ni] >= ORTHOLOG_ORIENCOUNT_LIM) )
				{
					dMatchSum += pMatchCount->pMatElement[ni];
					dOrienSum += pOrienCount->pMatElement[ni];
				}
			}
			
			/* update optimal mapping */
			if((dMatchSum+ORTHOLOG_COMBCOUNT_COEF*dOrienSum) > (dOptMatchSum+ORTHOLOG_COMBCOUNT_COEF*dOptOrienSum))
			{
				dOptMatchSum = dMatchSum;
				dOptOrienSum = dOrienSum;
				if(pOptMatchCount != NULL)
					DestroyDoubleMatrix(pOptMatchCount);
				pOptMatchCount = pMatchCount;
				if(pOptOrienCount != NULL)
					DestroyDoubleMatrix(pOptOrienCount);
				pOptOrienCount = pOrienCount;

				for(ni=0; ni<nSpeciesNum; ni++)
				{
					if(vOptRefGene[ni] != NULL)
						RefGeneDestroy(vOptRefGene[ni]);
					vOptRefGene[ni] = vRefGene[ni];
				}
				free(vRefGene);
				nSwitchon = 1;
			}
			else if(nSwitchon == 0)
			{
				if(nAddNum == 0)
				{
					pOptMatchCount = pMatchCount;
					pOptOrienCount = pOrienCount;

					vOptRefGene[0] = vRefGene[0];
					for(ni=1; ni<nSpeciesNum; ni++)
					{
						pOptMatchCount->pMatElement[ni] = 0.0;
						pOptOrienCount->pMatElement[ni] = 0.0;
						RefGeneDestroy(vRefGene[ni]);
					}
					free(vRefGene);
				}
				else
				{
					if((pMatchCount->pMatElement[0]+ORTHOLOG_COMBCOUNT_COEF*(pOrienCount->pMatElement[0])) > (pOptMatchCount->pMatElement[0]+ORTHOLOG_COMBCOUNT_COEF*(pOptOrienCount->pMatElement[0])))
					{
						pOptMatchCount->pMatElement[0] = pMatchCount->pMatElement[0];
						pOptOrienCount->pMatElement[0] = pOrienCount->pMatElement[0];
						pRefGene = vOptRefGene[0];
						vOptRefGene[0] = vRefGene[0];
						vRefGene[0] = pRefGene;
					}
					DestroyDoubleMatrix(pMatchCount);
					DestroyDoubleMatrix(pOrienCount);
					for(ni=0; ni<nSpeciesNum; ni++)
					{
						RefGeneDestroy(vRefGene[ni]);
						vRefGene[ni] = NULL;
					}
					free(vRefGene);
				}
			}
			else
			{
				DestroyDoubleMatrix(pMatchCount);
				DestroyDoubleMatrix(pOrienCount);
				for(ni=0; ni<nSpeciesNum; ni++)
				{
					RefGeneDestroy(vRefGene[ni]);
					vRefGene[ni] = NULL;
				}
				free(vRefGene);
			}
			nAddNum++;
		}
	}

	if(nInitMap == 1)
	{
		fprintf(fpOut, ">%s\n", strRefName0);
		for(ni=0; ni<nSpeciesNum; ni++)
		{
			if((pOptMatchCount != NULL) && (pOptOrienCount != NULL))
			{
				fprintf(fpOut, "%s\t%d\t%d\t", vSpeciesName[ni]->m_pString, (int)(pOptMatchCount->pMatElement[ni]),
					(int)(pOptOrienCount->pMatElement[ni]));
			}
			else
			{
				fprintf(fpOut, "%s\t0\t0\t", vSpeciesName[ni]->m_pString);
			}
			if(vOptRefGene[ni] != NULL)
			{
				RefGeneWrite(vOptRefGene[ni], fpOut);
			}
			else
			{
				fprintf(fpOut, "---\n");
			}
		}
		fprintf(fpOut, "\n");
		nOrthoGroupNum++;

		DestroyDoubleMatrix(pOptMatchCount);
		DestroyDoubleMatrix(pOptOrienCount);
		for(ni=0; ni<nSpeciesNum; ni++)
		{
			RefGeneDestroy(vOptRefGene[ni]);
			vOptRefGene[ni] = NULL;
		}
		free(vOptRefGene);
	}

	/* close files */
	for(ni=0; ni<nSpeciesNum; ni++)
	{
		fclose(vfpIn[ni]);
		vfpIn[ni] = NULL;
	}
	free(vfpIn);
	fclose(fpOut);

	/* return */
	return nOrthoGroupNum;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGene_GetMultiOrtholog1way_Main()                                    */
/*  Get ortholog genes from 1way refgene mapping.                          */
/* ----------------------------------------------------------------------- */ 
int RefGene_GetMultiOrtholog1way_Main(char strTargetPath[], char strParamPath[], 
								  char strOutPath[])
{
	/* define */
	int nMapSpeciesNum = 0;
	struct tagString **vMapSpeciesName;
	struct tagString **vRefSpeciesName;
	struct tagString **vRefSpeciesMapbase;
	struct tagString **vMapSpeciesMapbase;
	struct tagString **vMapSpeciesDatabase;

	FILE *fpIn;
	char strLine[LONG_LINE_LENGTH];
	char strInFileName[MED_LINE_LENGTH];
	char strOutFileName[MED_LINE_LENGTH];

	/* others */
	char *chp;
	int ni;
	
	/* init */
	vMapSpeciesName = NULL;
	vRefSpeciesName = NULL;
	vRefSpeciesMapbase = NULL;
	vMapSpeciesMapbase = NULL;
	vMapSpeciesDatabase = NULL;

	/* load parameters */
	fpIn = NULL;
	fpIn = fopen(strParamPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: cannot load parameters!\n");
		exit(EXIT_FAILURE);
	}

	/* load line by line */
	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strstr(strLine, "[Map Species Number]") == strLine)
		{
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);
			nMapSpeciesNum = atoi(chp);
			if(nMapSpeciesNum <= 0)
			{
				printf("Error: RefGene_GetMultiOrtholog1way_Main, map species number <= 0!\n");
				exit(EXIT_FAILURE);
			}

			vRefSpeciesName = (struct tagString **)calloc(nMapSpeciesNum, sizeof(struct tagString*));
			if(vRefSpeciesName == NULL)
			{
				printf("Error: RefGene_GetMultiOrtholog1way_Main, cannot create database\n");
				exit(EXIT_FAILURE);
			}

			vMapSpeciesName = (struct tagString **)calloc(nMapSpeciesNum, sizeof(struct tagString*));
			if(vMapSpeciesName == NULL)
			{
				printf("Error: RefGene_GetMultiOrtholog1way_Main, cannot create database\n");
				exit(EXIT_FAILURE);
			}

			vRefSpeciesMapbase = (struct tagString **)calloc(nMapSpeciesNum, sizeof(struct tagString*));
			if(vRefSpeciesMapbase == NULL)
			{
				printf("Error: RefGene_GetMultiOrtholog1way_Main, cannot create database\n");
				exit(EXIT_FAILURE);
			}

			vMapSpeciesDatabase = (struct tagString **)calloc(nMapSpeciesNum, sizeof(struct tagString*));
			if(vMapSpeciesDatabase == NULL)
			{
				printf("Error: RefGene_GetMultiOrtholog1way_Main, cannot create database\n");
				exit(EXIT_FAILURE);
			}

			vMapSpeciesMapbase = (struct tagString **)calloc(nMapSpeciesNum, sizeof(struct tagString*));
			if(vMapSpeciesMapbase == NULL)
			{
				printf("Error: RefGene_GetMultiOrtholog1way_Main, cannot create database\n");
				exit(EXIT_FAILURE);
			}
		}
		else if(strstr(strLine, "[Map Species Name]") == strLine)
		{
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);
			StringAddTail(vMapSpeciesName+ni, chp);

			/* load ref species name */
			fgets(strLine, LONG_LINE_LENGTH, fpIn);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "[Ref Species Name]") != strLine)
			{
				printf("Error: RefGene_GetMultiOrtholog1way_Main, parameter wrong\n");
				exit(EXIT_FAILURE);
			}
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);
			StringAddTail(vRefSpeciesName+ni, chp);

			/* load ref species mapbase */
			fgets(strLine, LONG_LINE_LENGTH, fpIn);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "[Ref Species RefGeneMap]") != strLine)
			{
				printf("Error: RefGene_GetMultiOrtholog1way_Main, parameter wrong\n");
				exit(EXIT_FAILURE);
			}
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);
			StringAddTail(vRefSpeciesMapbase+ni, chp);

			/* load map species mapbase */
			fgets(strLine, LONG_LINE_LENGTH, fpIn);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "[Map Species RefGeneMap]") != strLine)
			{
				printf("Error: RefGene_GetMultiOrtholog1way_Main, parameter wrong\n");
				exit(EXIT_FAILURE);
			}
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);
			StringAddTail(vMapSpeciesMapbase+ni, chp);

			/* load map species database */
			fgets(strLine, LONG_LINE_LENGTH, fpIn);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "[Map Species RefGene]") != strLine)
			{
				printf("Error: RefGene_GetMultiOrtholog1way_Main, parameter wrong\n");
				exit(EXIT_FAILURE);
			}
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);
			StringAddTail(vMapSpeciesDatabase+ni, chp);

			/* add the species number */
			ni++;
		}
		else
		{
			printf("Error: RefGene_GetMultiOrtholog1way_Main, unknown parameters\n");
			exit(EXIT_FAILURE);
		}
	}

	if(ni != nMapSpeciesNum)
	{
		printf("Error: RefGene_GetMultiOrtholog1way_Main, species number not match\n");
		exit(EXIT_FAILURE);
	}

	/* close files */
	fclose(fpIn);

	/* get orthologs from species to species */
	if(nMapSpeciesNum == 1)
	{
		RefGene_GetMultiOrtholog1way(strTargetPath, strOutPath,
			vRefSpeciesName[0]->m_pString, vMapSpeciesName[0]->m_pString,
			vRefSpeciesMapbase[0]->m_pString, vMapSpeciesMapbase[0]->m_pString,
			vMapSpeciesDatabase[0]->m_pString);
	}
	else
	{
		sprintf(strOutFileName, "%s0.tmp", strOutPath);
		RefGene_GetMultiOrtholog1way(strTargetPath, strOutFileName,
			vRefSpeciesName[0]->m_pString, vMapSpeciesName[0]->m_pString,
			vRefSpeciesMapbase[0]->m_pString, vMapSpeciesMapbase[0]->m_pString,
			vMapSpeciesDatabase[0]->m_pString);
		
		for(ni=1; ni<(nMapSpeciesNum-1); ni++)
		{
			sprintf(strInFileName, "%s%d.tmp", strOutPath, (ni-1));
			sprintf(strOutFileName, "%s%d.tmp", strOutPath, ni);
			RefGene_GetMultiOrtholog1way(strInFileName, strOutFileName,
				vRefSpeciesName[ni]->m_pString, vMapSpeciesName[ni]->m_pString,
				vRefSpeciesMapbase[ni]->m_pString, vMapSpeciesMapbase[ni]->m_pString,
				vMapSpeciesDatabase[ni]->m_pString);
		}

		sprintf(strInFileName, "%s%d.tmp", strOutPath, (ni-1));
		RefGene_GetMultiOrtholog1way(strInFileName, strOutPath,
			vRefSpeciesName[ni]->m_pString, vMapSpeciesName[ni]->m_pString,
			vRefSpeciesMapbase[ni]->m_pString, vMapSpeciesMapbase[ni]->m_pString,
			vMapSpeciesDatabase[ni]->m_pString);

		/* delete temp files */
		if(strcmp(OS_SYSTEM, "UNIX") == 0)
		{
			sprintf(strLine, "rm %s*.tmp", strOutPath);
			system(strLine);
		}
		else
		{
			sprintf(strLine, "del %s*.tmp", strOutPath);
			system(strLine);
		}
	}

	/* release memory */
	for(ni=0; ni<nMapSpeciesNum; ni++)
	{
		DeleteString(vRefSpeciesName[ni]);
		DeleteString(vMapSpeciesName[ni]);
		DeleteString(vRefSpeciesMapbase[ni]);
		DeleteString(vMapSpeciesMapbase[ni]);
		DeleteString(vMapSpeciesDatabase[ni]);
	}
	free(vRefSpeciesName);
	free(vMapSpeciesName);
	free(vRefSpeciesMapbase);
	free(vMapSpeciesMapbase);
	free(vMapSpeciesDatabase);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGene_GetMultiOrtholog1way()                                         */
/*  Get ortholog genes from transcript mapping. This mapping does not      */
/*  require colinearity as getmultiortholog and is therefore much less     */
/*  stringent.                                                             */
/* ----------------------------------------------------------------------- */ 
int RefGene_GetMultiOrtholog1way(char strInPath[], char strOutPath[],
			char strRefSpeciesName[], char strMapSpeciesName[],
			char strRefSpeciesMapbase[], char strMapSpeciesMapbase[],
			char strMapSpeciesDatabase[])
{
	/* define */
	FILE *fpSourceRefGene;
	struct tagRefGene *pSourceRefGeneList;
	struct tagRefGene **vSourceRefGene;
	int nSourceRefNum;

	FILE *fpMapRefGene;
	struct tagRefGene *pMapRefGeneList;
	struct tagRefGene **vMapRefGene;
	int nMapRefNum;

	int nDestMap;
	FILE *fpDestRefGene;
	struct tagRefGene *pDestRefGeneList;
	struct tagRefGene **vDestRefGene;
	int nDestRefNum;

	struct tagRefGene *pRefGene, *pCurrentRefGene;
	char strRefLine[LONG_LINE_LENGTH];
	int ni;

	/* for target search */
	FILE *fpIn;
	FILE *fpOut;
	int nFirstRefGene = 1;
	char strSpecies[LINE_LENGTH];
	int nSeedMatchId,nOrthoMatchId,nDestMatchId;
	struct tagRefGene *pSeedRefGene,*pRefMapRefGene,*pMapMapRefGene,*pMapDataRefGene;
	double dOptOrienCount,dOptMatchCount;
	
	char *chp,*chp2;


	/* ----------------------------- */
	/* init databases                */
	/* ----------------------------- */
	if(strcmp(strMapSpeciesDatabase, "NULL") == 0)
	{
		nDestMap = 0;
	}
	else
	{
		nDestMap = 1;
	}

	pSourceRefGeneList = NULL;
	pMapRefGeneList = NULL;
	pDestRefGeneList = NULL;
	nSourceRefNum = 0;
	nMapRefNum = 0;
	nDestRefNum = 0;
	vSourceRefGene = NULL;
	vMapRefGene = NULL;
	vDestRefGene = NULL;

	/* load source refgene */
	fpSourceRefGene = NULL;
	fpSourceRefGene = fopen(strRefSpeciesMapbase, "r");
	if(fpSourceRefGene == NULL)
	{
		printf("Error: RefGene_GetMultiOrtholog1way, cannot open source refgene file!\n");
		exit(EXIT_FAILURE);
	}

	pCurrentRefGene = NULL;
	while(fgets(strRefLine, LONG_LINE_LENGTH, fpSourceRefGene) != NULL)
	{
		StrTrimLeft(strRefLine);
		StrTrimRight(strRefLine);
		if(strRefLine[0] == '\0')
			continue;

		pRefGene = NULL;
		pRefGene = RefGeneCreate();
		RefGeneInit_FromGenomeLabFormat(pRefGene, strRefLine, strRefSpeciesName);
		if(pSourceRefGeneList == NULL)
		{
			pSourceRefGeneList = pRefGene;
			pCurrentRefGene = pRefGene;
		}
		else
		{
			pCurrentRefGene->pNext = pRefGene;
			pCurrentRefGene = pRefGene;
		}
		nSourceRefNum++;
	}

	fclose(fpSourceRefGene);

	vSourceRefGene = (struct tagRefGene **)calloc(nSourceRefNum, sizeof(struct tagRefGene*));
	if(vSourceRefGene == NULL)
	{
		printf("Error: RefGene_GetMultiOrtholog1way, organize refgene into array!\n");
		exit(EXIT_FAILURE);
	}

	ni=0;
	while(pSourceRefGeneList != NULL)
	{
		pRefGene = pSourceRefGeneList;
		pSourceRefGeneList = pRefGene->pNext;
		pRefGene->pNext = NULL;
		vSourceRefGene[ni] = pRefGene;
		ni++;
	}

	if(ni != nSourceRefNum)
	{
		printf("Error: RefGene_GetMultiOrtholog1way, refgene number not match!\n");
		exit(EXIT_FAILURE);
	}


	/* load map refgene */
	fpMapRefGene = NULL;
	fpMapRefGene = fopen(strMapSpeciesMapbase, "r");
	if(fpMapRefGene == NULL)
	{
		printf("Error: RefGene_GetMultiOrtholog1way, cannot open map refgene file!\n");
		exit(EXIT_FAILURE);
	}

	pCurrentRefGene = NULL;
	while(fgets(strRefLine, LONG_LINE_LENGTH, fpMapRefGene) != NULL)
	{
		StrTrimLeft(strRefLine);
		StrTrimRight(strRefLine);
		if(strRefLine[0] == '\0')
			continue;

		pRefGene = NULL;
		pRefGene = RefGeneCreate();
		RefGeneInit_FromGenomeLabFormat(pRefGene, strRefLine, strMapSpeciesName);
		if(pMapRefGeneList == NULL)
		{
			pMapRefGeneList = pRefGene;
			pCurrentRefGene = pRefGene;
		}
		else
		{
			pCurrentRefGene->pNext = pRefGene;
			pCurrentRefGene = pRefGene;
		}
		nMapRefNum++;
	}

	fclose(fpMapRefGene);

	vMapRefGene = (struct tagRefGene **)calloc(nMapRefNum, sizeof(struct tagRefGene*));
	if(vMapRefGene == NULL)
	{
		printf("Error: RefGene_GetMultiOrtholog1way, cannot organize refgene into array!\n");
		exit(EXIT_FAILURE);
	}

	ni=0;
	while(pMapRefGeneList != NULL)
	{
		pRefGene = pMapRefGeneList;
		pMapRefGeneList = pRefGene->pNext;
		pRefGene->pNext = NULL;
		vMapRefGene[ni] = pRefGene;
		ni++;
	}

	if(ni != nMapRefNum)
	{
		printf("Error: RefGene_GetMultiOrtholog1way, refgene number not match!\n");
		exit(EXIT_FAILURE);
	}


	/* load dest refgene */
	if(nDestMap == 1)
	{
		fpDestRefGene = NULL;
		fpDestRefGene = fopen(strMapSpeciesDatabase, "r");
		if(fpDestRefGene == NULL)
		{
			printf("Error: RefGene_GetMultiOrtholog1way, cannot open dest refgene file!\n");
			exit(EXIT_FAILURE);
		}

		pCurrentRefGene = NULL;
		while(fgets(strRefLine, LONG_LINE_LENGTH, fpDestRefGene) != NULL)
		{
			StrTrimLeft(strRefLine);
			StrTrimRight(strRefLine);
			if(strRefLine[0] == '\0')
				continue;

			pRefGene = NULL;
			pRefGene = RefGeneCreate();
			RefGeneInit_FromGenomeLabFormat(pRefGene, strRefLine, strMapSpeciesName);
			if(pDestRefGeneList == NULL)
			{
				pDestRefGeneList = pRefGene;
				pCurrentRefGene = pRefGene;
			}
			else
			{
				pCurrentRefGene->pNext = pRefGene;
				pCurrentRefGene = pRefGene;
			}
			nDestRefNum++;
		}

		fclose(fpDestRefGene);

		vDestRefGene = (struct tagRefGene **)calloc(nDestRefNum, sizeof(struct tagRefGene*));
		if(vDestRefGene == NULL)
		{
			printf("Error: RefGene_GetMultiOrtholog1way, organize refgene into array!\n");
			exit(EXIT_FAILURE);
		}

		ni=0;
		while(pDestRefGeneList != NULL)
		{
			pRefGene = pDestRefGeneList;
			pDestRefGeneList = pRefGene->pNext;
			pRefGene->pNext = NULL;
			vDestRefGene[ni] = pRefGene;
			ni++;
		}

		if(ni != nDestRefNum)
		{
			printf("Error: RefGene_GetMultiOrtholog1way, refgene number not match!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* ----------------------------- */
	/* locate ortholog               */
	/* ----------------------------- */
	pSeedRefGene = NULL;
	pRefMapRefGene = NULL;
	pMapMapRefGene = NULL;
	pMapDataRefGene = NULL;

	fpIn = NULL;
	fpIn = fopen(strInPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: RefGene_GetMultiOrtholog1way, cannot open input map file!\n");
		exit(EXIT_FAILURE);
	}
	fpOut = NULL;
	fpOut = fopen(strOutPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: RefGene_GetMultiOrtholog1way, cannot open output map file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strRefLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strRefLine);
		StrTrimRight(strRefLine);
		if(strRefLine[0] == '\0')
			continue;

		if(strRefLine[0] == '>')
		{
			if(pSeedRefGene == NULL)
			{
				if(nFirstRefGene == 1)
				{
					nFirstRefGene = 0;
				}
				else
				{
					fprintf(fpOut, "%s\t0\t0\t---\n", strMapSpeciesName);
					fprintf(fpOut, "\n");
				}
			}
			else
			{
				/* --------------------- */
				/* step I: --> refmap    */
				/* --------------------- */
				nSeedMatchId = RefGene_GetBestOverlap(pSeedRefGene, vSourceRefGene, nSourceRefNum);
				if(nSeedMatchId < 0)
				{
					pRefMapRefGene = NULL;
				}
				else
				{
					pRefMapRefGene = vSourceRefGene[nSeedMatchId];
				}

				/* --------------------- */
				/* step II: --> mapmap   */
				/* --------------------- */
				if(nSeedMatchId >= 0)
				{
					nOrthoMatchId = RefGene_GetSeededRelaxedOrtholog(nSeedMatchId, vSourceRefGene, nSourceRefNum, 
						vMapRefGene, nMapRefNum, &dOptOrienCount, &dOptMatchCount);
					if(nOrthoMatchId < 0)
					{
						pMapMapRefGene = NULL;
					}
					else
					{
						pMapMapRefGene = vMapRefGene[nOrthoMatchId];
					}
				}

				/* --------------------- */
				/* step III: --> mapdata */
				/* --------------------- */
				if((nDestMap == 1) && (nOrthoMatchId >= 0))
				{
					nDestMatchId = RefGene_GetBestOverlap(pMapMapRefGene, vDestRefGene, nDestRefNum);
					if(nDestMatchId < 0)
					{
						pMapDataRefGene = NULL;
					}
					else
					{
						pMapDataRefGene = vDestRefGene[nDestMatchId];
					}
				}

				/* --------------------- */
				/* step IV: --> export   */
				/* --------------------- */
				fprintf(fpOut, "%s\t%d\t%d\t", strMapSpeciesName, (int)dOptMatchCount,
					(int)dOptOrienCount);
				if( (nDestMap == 1) && (pMapDataRefGene != NULL) )
				{
					RefGeneWrite(pMapDataRefGene, fpOut);
				}
				else
				{
					if(pMapMapRefGene == NULL)
					{
						fprintf(fpOut, "---\n");
					}
					else
					{
						RefGeneWrite(pMapMapRefGene, fpOut);
					}
				}

				fprintf(fpOut, "\n");
			}

			/* --------------------- */
			/* step V: --> reinit    */
			/* --------------------- */
			nSeedMatchId = -1;
			nOrthoMatchId = -1;
			nDestMatchId = -1;
			RefGeneDestroy(pSeedRefGene);
			pSeedRefGene = NULL;
			pRefMapRefGene = NULL;
			pMapMapRefGene = NULL;
			pMapDataRefGene = NULL;
			dOptOrienCount = 0.0;
			dOptMatchCount = 0.0;
		}
		else
		{
			/* ------------------ */
			/* get ref species id */
			/* ------------------ */
			sscanf(strRefLine, "%s", strSpecies);
			if(strcmp(strSpecies, strRefSpeciesName) == 0)
			{
				chp = strchr(strRefLine, '\t');
				chp++;
				
				chp2 = strchr(chp, '\t');
				chp = chp2+1;
				chp2 = strchr(chp, '\t');
				chp = chp2+1;

				if(strcmp(chp, "---") != 0)
				{
					if(pSeedRefGene != NULL)
					{
						printf("Error: RefGene_GetMultiOrtholog1way, there are more than one seed refgene!\n");
						exit(EXIT_FAILURE);
					}

					pSeedRefGene = RefGeneCreate();
					RefGeneInit_FromGenomeLabFormat(pSeedRefGene, chp, strRefSpeciesName);
				}
			}
		}

		fprintf(fpOut, "%s\n", strRefLine);
	}

	/* -------------------- */
	/* process the last one */
	/* -------------------- */
	if(pSeedRefGene == NULL)
	{
		if(nFirstRefGene == 1)
		{
			nFirstRefGene = 0;
		}
		else
		{
			fprintf(fpOut, "%s\t0\t0\t---\n", strMapSpeciesName);
			fprintf(fpOut, "\n");
		}
	}
	else
	{
		/* --------------------- */
		/* step I: --> refmap    */
		/* --------------------- */
		nSeedMatchId = RefGene_GetBestOverlap(pSeedRefGene, vSourceRefGene, nSourceRefNum);
		if(nSeedMatchId < 0)
		{
			pRefMapRefGene = NULL;
		}
		else
		{
			pRefMapRefGene = vSourceRefGene[nSeedMatchId];
		}

		/* --------------------- */
		/* step II: --> mapmap   */
		/* --------------------- */
		if(nSeedMatchId >= 0)
		{
			nOrthoMatchId = RefGene_GetSeededRelaxedOrtholog(nSeedMatchId, vSourceRefGene, nSourceRefNum, 
				vMapRefGene, nMapRefNum, &dOptOrienCount, &dOptMatchCount);
			if(nOrthoMatchId < 0)
			{
				pMapMapRefGene = NULL;
			}
			else
			{
				pMapMapRefGene = vMapRefGene[nOrthoMatchId];
			}
		}

		/* --------------------- */
		/* step III: --> mapdata */
		/* --------------------- */
		if((nDestMap == 1) && (nOrthoMatchId >= 0))
		{
			nDestMatchId = RefGene_GetBestOverlap(pMapMapRefGene, vDestRefGene, nDestRefNum);
			if(nDestMatchId < 0)
			{
				pMapDataRefGene = NULL;
			}
			else
			{
				pMapDataRefGene = vDestRefGene[nDestMatchId];
			}
		}

		/* --------------------- */
		/* step IV: --> export   */
		/* --------------------- */
		fprintf(fpOut, "%s\t%d\t%d\t", strMapSpeciesName, (int)dOptMatchCount,
			(int)dOptOrienCount);
		if( (nDestMap == 1) && (pMapDataRefGene != NULL) )
		{
			RefGeneWrite(pMapDataRefGene, fpOut);
		}
		else
		{
			if(pMapMapRefGene == NULL)
			{
				fprintf(fpOut, "---\n");
			}
			else
			{
				RefGeneWrite(pMapMapRefGene, fpOut);
			}
		}

		fprintf(fpOut, "\n");
		RefGeneDestroy(pSeedRefGene);
		pSeedRefGene = NULL;
	}

	fclose(fpIn);
	fclose(fpOut);

	/* ----------------------------- */
	/* release memory                */
	/* ----------------------------- */
	for(ni=0; ni<nSourceRefNum; ni++)
	{
		RefGeneDestroy(vSourceRefGene[ni]);
		vSourceRefGene[ni] = NULL;
	}
	free(vSourceRefGene);

	for(ni=0; ni<nMapRefNum; ni++)
	{
		RefGeneDestroy(vMapRefGene[ni]);
		vMapRefGene[ni] = NULL;
	}
	free(vMapRefGene);

	if(nDestMap == 1)
	{
		for(ni=0; ni<nDestRefNum; ni++)
		{
			RefGeneDestroy(vDestRefGene[ni]);
			vDestRefGene[ni] = NULL;
		}
		free(vDestRefGene);
	}


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGene_GetBestOverlap()                                               */
/*  Get best overlap from a refgene database.                              */
/* ----------------------------------------------------------------------- */ 
int RefGene_GetBestOverlap(struct tagRefGene *pSeedRefGene, 
				struct tagRefGene **vRefGeneDatabase, int nRefGeneNum)
{
	/* define */
	int nBestId;
	struct tagRefGene *pBestRefGene;
	double dBestScore,dScore;
	int ni;
	int nCDSL1,nCDSL2,nCDSL,nCDSS,nCDSE;
	double dCDSR1,dCDSR2;

	/* search for best match */
	nBestId = -1;
	pBestRefGene = NULL;
	dBestScore = 0.0;

	for(ni=0; ni<nRefGeneNum; ni++)
	{
		dScore = RefGene_Match(vRefGeneDatabase[ni], pSeedRefGene);
		if(dScore > dBestScore)
		{
			if(dScore > REFGENE_MATCH_TH)
			{
				nBestId = ni;
				dBestScore = dScore;
				pBestRefGene = vRefGeneDatabase[ni];
			}
			else
			{
				nCDSS = pSeedRefGene->nCdsStart;
				if(vRefGeneDatabase[ni]->nCdsStart > nCDSS)
					nCDSS = vRefGeneDatabase[ni]->nCdsStart;
				nCDSE = pSeedRefGene->nCdsEnd;
				if(vRefGeneDatabase[ni]->nCdsEnd < nCDSE)
					nCDSE = vRefGeneDatabase[ni]->nCdsEnd;

				nCDSL1 = pSeedRefGene->nCdsEnd-pSeedRefGene->nCdsStart+1;
				nCDSL2 = vRefGeneDatabase[ni]->nCdsEnd-vRefGeneDatabase[ni]->nCdsStart+1;
				nCDSL = nCDSE-nCDSS+1;

				dCDSR1 = (double)nCDSL/(double)nCDSL1;
				dCDSR2 = (double)nCDSL/(double)nCDSL2;
				if(dCDSR2<dCDSR1)
					dCDSR1 = dCDSR2;

				if(dCDSR1 > REFGENE_MATCH_TH)
				{
					nBestId = ni;
					dBestScore = dScore;
					pBestRefGene = vRefGeneDatabase[ni];
				}
			}
		}
	}

	/* return */
	return nBestId;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGene_GetSeededRelaxedOrtholog()                                     */
/*  Get ortholog genes for a seed refgene. No colinearity is required.     */
/* ----------------------------------------------------------------------- */ 
int RefGene_GetSeededRelaxedOrtholog(int nSeedMatchId, 
					struct tagRefGene **vSourceRefGene, int nSourceRefNum, 
					struct tagRefGene **vMapRefGene, int nMapRefNum,
					double *pOptOrienCount, double *pOptMatchCount)
{
	/* define */

	/* for database */
	FILE *fpTemp;
	int nMapHitNum;
	struct INTMATRIX *pHitLimit;

	int nPos1,nPos2,nPos3;
	int nMatchNum;
	int nOrientationNum;
	int nMatched;

	/* best id */
	int nBestId;
	double dBestScore,dScore;
	int nBestMatchNum;
	int nBestOrienNum;

	char strLine[MED_LINE_LENGTH];

	char strRefName[LINE_LENGTH];
	int vRefSourceIdx[ORTHOLOG_COLINEAR_NUM];
	char vRefColinear[ORTHOLOG_COLINEAR_NUM][LINE_LENGTH];
	int ni,nj,nk,nx,ny,nz;
	int nflank,nlastid;

	/* init */
	
	/* look up source database */
	nflank = (ORTHOLOG_COLINEAR_NUM-1)/2;
	for(nk=0; nk<ORTHOLOG_COLINEAR_NUM; nk++)
	{
		vRefSourceIdx[nk] = -1;
		strcpy(vRefColinear[nk], "---");
	}

	strcpy(strRefName, vSourceRefGene[nSeedMatchId]->strName);
	vRefSourceIdx[nflank] = nSeedMatchId;
	strcpy(vRefColinear[nflank], vSourceRefGene[nSeedMatchId]->strName);

	/* get upstream flank */
	nk = 0;
	nlastid = nSeedMatchId;
	for(nj=nSeedMatchId-1; nj>=0; nj--)
	{
		if(vSourceRefGene[nSeedMatchId]->nChrom == vSourceRefGene[nj]->nChrom)
		{
			if((vSourceRefGene[nSeedMatchId]->nTxStart-vSourceRefGene[nj]->nTxEnd) <= ORTHOLOG_COLINEAR_WIN)
			{
				if(vSourceRefGene[nj]->nTxEnd < vSourceRefGene[nlastid]->nTxStart)
				{
					nk++;
					vRefSourceIdx[nflank-nk] = nj;
					nlastid = nj;
					strcpy(vRefColinear[nflank-nk], vSourceRefGene[nj]->strName);
				}

				if(nk == nflank)
					break;
			}
			else
			{
				break;
			}

		}
		else
		{
			break;
		}
	}

	/* get downstream flank */
	nk = 0;
	nlastid = nSeedMatchId;
	for(nj=nSeedMatchId+1; nj<nSourceRefNum; nj++)
	{
		if(vSourceRefGene[nSeedMatchId]->nChrom == vSourceRefGene[nj]->nChrom)
		{
			if((vSourceRefGene[nj]->nTxStart-vSourceRefGene[nSeedMatchId]->nTxEnd) <= ORTHOLOG_COLINEAR_WIN)
			{
				if( vSourceRefGene[nlastid]->nTxEnd < vSourceRefGene[nj]->nTxStart)
				{
					nk++;
					vRefSourceIdx[nflank+nk] = nj;
					nlastid = nj;
					strcpy(vRefColinear[nflank+nk], vSourceRefGene[nj]->strName);
				}

				if(nk == nflank)
					break;
			}
			else
			{
				break;
			}
		}
		else
		{
			break;
		}
	}
				

	/* remove redundancy */
	for(nk=0; nk<ORTHOLOG_COLINEAR_NUM; nk++)
	{
		if(vRefSourceIdx[nk] != -1)
		{
			for(nx=0; nx<nk; nx++)
			{
				if(vRefSourceIdx[nx] != -1)
				{
					if(strcmp(vRefColinear[nk], vRefColinear[nx]) == 0)
					{
						strcpy(vRefColinear[nk], "---");
						vRefSourceIdx[nk] = -1;
						break;
					}
				}
			}
		}
	}
	
	/* look up map database */
	nMapHitNum = 0;
	fpTemp = NULL;
	fpTemp = fopen("maphits.tmp", "w");
	if(fpTemp == NULL)
	{
		printf("Error: cannot create temp file!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nMapRefNum; ni++)
	{
		/* if found the location in the map organism */
		if(strcmp(strRefName, vMapRefGene[ni]->strName) == 0)
		{
			nMapHitNum++;

			/* get upstream flank */
			nk = 0;
			nlastid = ni;
			for(nj=ni-1; nj>=0; nj--)
			{
				if(vMapRefGene[ni]->nChrom == vMapRefGene[nj]->nChrom)
				{
					if( (vMapRefGene[ni]->nTxStart-vMapRefGene[nj]->nTxEnd) <= ORTHOLOG_COLINEAR_WIN )
					{
						if(vMapRefGene[nj]->nTxEnd < vMapRefGene[nlastid]->nTxStart)
						{
							nk++;
							if(nk > ORTHOLOG_MAP_COLINEAR_NUM)
							{
								nlastid = nj+1;
								break;
							}
							nlastid = nj;
						}
					}
					else
					{
						nlastid = nj+1;
						break;
					}
				}
				else
				{
					nlastid = nj+1;
					break;
				}
			}

			fprintf(fpTemp, "%d\t", nlastid);
			fprintf(fpTemp, "%d\t", ni);

			/* get downstream flank */
			nk = 0;
			nlastid = ni;
			for(nj=ni+1; nj<nMapRefNum; nj++)
			{
				if(vMapRefGene[ni]->nChrom == vMapRefGene[nj]->nChrom)
				{
					if( (vMapRefGene[nj]->nTxStart-vMapRefGene[ni]->nTxEnd) <= ORTHOLOG_COLINEAR_WIN )
					{
						if(vMapRefGene[nlastid]->nTxEnd < vMapRefGene[nj]->nTxStart)
						{
							nk++;
							if(nk > ORTHOLOG_MAP_COLINEAR_NUM)
							{
								nlastid = nj-1;
								break;
							}
							nlastid = nj;
						}
					}
					else
					{
						nlastid = nj-1;
						break;
					}
				}
				else
				{
					nlastid = nj-1;
					break;
				}
			}

			fprintf(fpTemp, "%d\n", nlastid);
		}
	}

	fclose(fpTemp);
	

	/* verify collinearity */
	if(nMapHitNum == 0)
	{
		return -1;
	}
	else
	{
		nBestId = -1;
		nBestMatchNum = 0;
		nBestOrienNum = 0;
		dBestScore = 0.0;

		pHitLimit = NULL;
		pHitLimit = IMLOAD("maphits.tmp");
		if(pHitLimit == NULL)
		{
			printf("Error: cannot reload maphits.tmp!\n");
			exit(EXIT_FAILURE);
		}
		if(pHitLimit->nHeight != nMapHitNum)
		{
			printf("Error: map hit number not match!\n");
			exit(EXIT_FAILURE);
		}

		/* match */
		for(nz=0; nz<nMapHitNum; nz++)
		{
			nMatchNum = 0;
			nOrientationNum = 0;
			nPos1 = IMGETAT(pHitLimit, nz, 0);
			nPos2 = IMGETAT(pHitLimit, nz, 1);
			nPos3 = IMGETAT(pHitLimit, nz, 2);

			for(nk=0; nk<ORTHOLOG_COLINEAR_NUM; nk++)
			{
				nMatched = 0;
				vRefSourceIdx[nk] = 0;
				for(ny=nPos1; ny<nPos2; ny++)
				{
					if(strcmp(vRefColinear[nk], vMapRefGene[ny]->strName) == 0)
					{
						nMatched = 1;
						nMatchNum++;
						vRefSourceIdx[nk] = 1;
						break;
					}
				}
				if(nMatched == 1)
					continue;

				if(strcmp(vRefColinear[nk], vMapRefGene[ny]->strName) == 0)
				{
					nMatchNum++;
					continue;
				}

				ny++;
				for(; ny<=nPos3; ny++)
				{
					if(strcmp(vRefColinear[nk], vMapRefGene[ny]->strName) == 0)
					{
						nMatched = 1;
						nMatchNum++;
						vRefSourceIdx[nk] = -1;
						break;
					}
				}
			}

			nflank = (ORTHOLOG_COLINEAR_NUM-1)/2;
			for(nk=0; nk<nflank; nk++)
			{
				for(nj=0; nj<nflank; nj++)
				{
					if(nk != nj)
					{
						if(vRefSourceIdx[nk]*vRefSourceIdx[nj] > 0)
						{
							nOrientationNum++;
						}
					}
				}
				nj++;
				for(; nj<ORTHOLOG_COLINEAR_NUM; nj++)
				{
					if(vRefSourceIdx[nk]*vRefSourceIdx[nj] < 0)
					{
						nOrientationNum++;
					}
				}
			}

			nk++;
			for(; nk<ORTHOLOG_COLINEAR_NUM; nk++)
			{
				for(nj=0; nj<nflank; nj++)
				{
					if(vRefSourceIdx[nk]*vRefSourceIdx[nj] < 0)
					{
						nOrientationNum++;
					}
				}
				nj++;
				for(; nj<ORTHOLOG_COLINEAR_NUM; nj++)
				{
					if(nj != nk)
					{
						if(vRefSourceIdx[nk]*vRefSourceIdx[nj] > 0)
						{
							nOrientationNum++;
						}
					}
				}
			}

			dScore = (double)nMatchNum+ORTHOLOG_COMBCOUNT_COEF*(double)nOrientationNum;
			if(dScore > dBestScore)
			{
				nBestId = nPos2;
				dBestScore = dScore;
				nBestMatchNum = nMatchNum;
				nBestOrienNum = nOrientationNum;
			}
			else if(fabs(dScore-dBestScore) < ZERO_BOUND)
			{
				if( rand_u() < 0.5 )
				{
					nBestId = nPos2;
					dBestScore = dScore;
					nBestMatchNum = nMatchNum;
					nBestOrienNum = nOrientationNum;
				}
			}
		}

		DestroyIntMatrix(pHitLimit);
	}

	/* delete temp files */
	if(strcmp(OS_SYSTEM, "UNIX") == 0)
	{
		sprintf(strLine, "rm maphits.tmp");
		system(strLine);
	}
	else
	{
		sprintf(strLine, "del maphits.tmp");
		system(strLine);
	}

	*pOptOrienCount = nBestOrienNum;
	*pOptMatchCount = nBestMatchNum;

	/* return */
	return nBestId;
}

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
					char strParamPath[], char strOutPath[])
{
	/* define */
	int nSpeciesNum;
	int nOrthoGroupNum;
	int nRefIdConvert;
	struct tagString **vSpeciesName;
	struct tagString **vSpeciesMapbase;
	struct tagString **vSpeciesDatabase;
	struct tagString **vPairwisePath;
	struct tagRefGene **vOriginRefGene = NULL;
	int nOriginRefNum = 0;

	FILE *fpIn;
	char strLine[MED_LINE_LENGTH];
	char *chp;
	int ni;
	char strTempOutPath[LINE_LENGTH];
	char strNROutPath[LINE_LENGTH];
	int nNeedMapBack;

	/* init */
	vSpeciesName = NULL;
	vSpeciesMapbase = NULL;
	vSpeciesDatabase = NULL;
	vPairwisePath = NULL;

	/* load parameters */
	fpIn = NULL;
	fpIn = fopen(strParamPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: cannot load parameters!\n");
		exit(EXIT_FAILURE);
	}

	/* load line by line */
	ni = 0;
	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strstr(strLine, "[Species Number]") == strLine)
		{
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);
			nSpeciesNum = atoi(chp);
			if(nSpeciesNum <= 1)
			{
				printf("Error: RefFlex_GetMultiOrtholog_Main, species number < 2!\n");
				exit(EXIT_FAILURE);
			}

			vSpeciesName = (struct tagString **)calloc(nSpeciesNum, sizeof(struct tagString*));
			if(vSpeciesName == NULL)
			{
				printf("Error: RefFlex_GetMultiOrtholog_Main, cannot create database\n");
				exit(EXIT_FAILURE);
			}

			vSpeciesMapbase = (struct tagString **)calloc(nSpeciesNum, sizeof(struct tagString*));
			if(vSpeciesMapbase == NULL)
			{
				printf("Error: RefFlex_GetMultiOrtholog_Main, cannot create database\n");
				exit(EXIT_FAILURE);
			}

			vSpeciesDatabase = (struct tagString **)calloc(nSpeciesNum, sizeof(struct tagString*));
			if(vSpeciesDatabase == NULL)
			{
				printf("Error: RefFlex_GetMultiOrtholog_Main, cannot create database\n");
				exit(EXIT_FAILURE);
			}

			vPairwisePath = (struct tagString **)calloc(nSpeciesNum, sizeof(struct tagString*));
			if(vPairwisePath == NULL)
			{
				printf("Error: RefFlex_GetMultiOrtholog_Main, cannot create database\n");
				exit(EXIT_FAILURE);
			}
		}
		else if(strstr(strLine, "[Reference Species]") == strLine)
		{
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);
			StringAddTail(vSpeciesName+ni, chp);

			fgets(strLine, MED_LINE_LENGTH, fpIn);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "[RefId From RefGeneMap(0) or RefGene(1)]") != strLine)
			{
				printf("Error: RefFlex_GetMultiOrtholog_Main, parameter wrong\n");
				exit(EXIT_FAILURE);
			}
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);
			nRefIdConvert = atoi(chp);

			fgets(strLine, MED_LINE_LENGTH, fpIn);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "[Reference RefGeneMap]") != strLine)
			{
				printf("Error: RefFlex_GetMultiOrtholog_Main, parameter wrong\n");
				exit(EXIT_FAILURE);
			}
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);
			StringAddTail(vSpeciesMapbase+ni, chp);

			fgets(strLine, MED_LINE_LENGTH, fpIn);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "[Reference RefGene]") != strLine)
			{
				printf("Error: RefFlex_GetMultiOrtholog_Main, parameter wrong\n");
				exit(EXIT_FAILURE);
			}
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);
			StringAddTail(vSpeciesDatabase+ni, chp);

			ni++;
		}
		else if(strstr(strLine, "[Map Species Name]") == strLine)
		{
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);
			StringAddTail(vSpeciesName+ni, chp);

			fgets(strLine, MED_LINE_LENGTH, fpIn);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "[Map Species RefGeneMap]") != strLine)
			{
				printf("Error: RefFlex_GetMultiOrtholog_Main, parameter wrong\n");
				exit(EXIT_FAILURE);
			}
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);
			StringAddTail(vSpeciesMapbase+ni, chp);

			fgets(strLine, MED_LINE_LENGTH, fpIn);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "[Map Species RefGene]") != strLine)
			{
				printf("Error: RefFlex_GetMultiOrtholog_Main, parameter wrong\n");
				exit(EXIT_FAILURE);
			}
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);
			StringAddTail(vSpeciesDatabase+ni, chp);

			ni++;
		}
		else
		{
			printf("Error: RefFlex_GetMultiOrtholog_Main, unknown parameters\n");
			exit(EXIT_FAILURE);
		}
	}

	if(ni != nSpeciesNum)
	{
		printf("Error: RefFlex_GetMultiOrtholog_Main, species number not match\n");
		exit(EXIT_FAILURE);
	}

	/* close files */
	fclose(fpIn);

	/* convert refernce map */
	vOriginRefGene = RefFlex_LoadOriginRefGene(strTargetPath, nColumn, 
		vSpeciesName[0]->m_pString, &nOriginRefNum);
	if(nOriginRefNum <= 0)
	{
		printf("Warning: RefFlex_GetMultiOrtholog_Main, no gene needs to be mapped!\n"); 
		return PROC_SUCCESS;
	}
	
	/* pairwise mapping */
	for(ni=0; ni<nSpeciesNum; ni++)
	{
		sprintf(strTempOutPath, "%s_%s", strOutPath, vSpeciesName[ni]->m_pString);
		StringAddTail(vPairwisePath+ni, strTempOutPath);
		RefFlex_GetOrtholog(vOriginRefGene, nOriginRefNum, 
			vSpeciesName[0]->m_pString, vSpeciesMapbase[0]->m_pString,
			vSpeciesName[ni]->m_pString, vSpeciesMapbase[ni]->m_pString,
			"NULL", strTempOutPath);
	}
	
	/* merge multiple species */
	sprintf(strTempOutPath, "%s.nstmap", strOutPath);
	nOrthoGroupNum = RefFlex_MergeMultiOrtholog(nSpeciesNum, vSpeciesName,
			strTargetPath, vOriginRefGene, nOriginRefNum, vPairwisePath, strTempOutPath);
	
	/* map back to the refgene database of the original species */
	nNeedMapBack = 0;
	for(ni=0; ni<nSpeciesNum; ni++)
	{
		if(strcmp(vSpeciesDatabase[ni]->m_pString, "NULL") != 0)
		{
			nNeedMapBack = 1;
			break;
		}
	}
	if(nNeedMapBack == 1)
	{
		sprintf(strTempOutPath, "%s.nsomap", strOutPath);
		sprintf(strNROutPath, "%s.nstmap", strOutPath);
		RefFlex_MapMultiOrthologToOriginalSpecies(nSpeciesNum, vSpeciesName, vSpeciesMapbase,
			vSpeciesDatabase, vOriginRefGene, nOriginRefNum, strNROutPath, strTempOutPath);
	}
	
	/* release memory */
	RefGene_ClearDatabase(&vOriginRefGene, nOriginRefNum);
	
	for(ni=0; ni<nSpeciesNum; ni++)
	{
		DeleteString(vSpeciesName[ni]);
		DeleteString(vSpeciesMapbase[ni]);
		DeleteString(vSpeciesDatabase[ni]);
		DeleteString(vPairwisePath[ni]);
	}
	free(vSpeciesName);
	free(vSpeciesMapbase);
	free(vSpeciesDatabase);
	free(vPairwisePath);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefFlex_LoadOriginRefGene()                                            */
/*  load refgene from the nColumn-th column of a file.                     */
/* ----------------------------------------------------------------------- */ 
struct tagRefGene **RefFlex_LoadOriginRefGene(char strInPath[], int nColumn, 
						char strSpecies[], int *pRefNum)
{
	/* define */
	FILE *fpRefGene;
	struct tagRefGene **vSourceRefGene = NULL;
	int nSourceRefNum = 0;
	char strRefLine[LONG_LINE_LENGTH];
	int ni,nj;
	char *chp1,*chp2;

	/* get refgene number */
	fpRefGene = NULL;
	fpRefGene = fopen(strInPath, "r");
	if(fpRefGene == NULL)
	{
		printf("Error: RefFlex_LoadOriginRefGene, cannot open source refgene file!\n");
		exit(EXIT_FAILURE);
	}

	nSourceRefNum = 0;
	while(fgets(strRefLine, LONG_LINE_LENGTH, fpRefGene) != NULL)
	{
		StrTrimLeft(strRefLine);
		StrTrimRight(strRefLine);
		if(strRefLine[0] == '\0')
			continue;

		nSourceRefNum++;
	}

	fclose(fpRefGene);
	if(nSourceRefNum <= 0)
	{
		printf("Error: RefFlex_LoadOriginRefGene, no gene available!\n");
		*pRefNum = nSourceRefNum;
		return NULL;
	}

	/* load */
	vSourceRefGene = NULL;
	vSourceRefGene = (struct tagRefGene **)calloc(nSourceRefNum, sizeof(struct tagRefGene*));
	if(vSourceRefGene == NULL)
	{
		printf("Error: RefFlex_LoadOriginRefGene, cannot allocate memory for loading refgene!\n");
		exit(EXIT_FAILURE);
	}

	fpRefGene = NULL;
	fpRefGene = fopen(strInPath, "r");
	if(fpRefGene == NULL)
	{
		printf("Error: RefFlex_LoadOriginRefGene, cannot open source refgene file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strRefLine, LONG_LINE_LENGTH, fpRefGene) != NULL)
	{
		StrTrimLeft(strRefLine);
		StrTrimRight(strRefLine);
		if(strRefLine[0] == '\0')
			continue;

		chp1 = strRefLine;
		for(nj=0; nj<nColumn; nj++)
		{
			chp2 = strchr(chp1, '\t');
			if(chp2 == NULL)
				break;

			chp1 = chp2+1;
		}

		if(nj != nColumn)
		{
			vSourceRefGene[ni] = NULL;
		}
		else
		{
			if( (strcmp(chp1, "---") == 0) || (strcmp(chp1, "NA") == 0) )
			{
				vSourceRefGene[ni] = NULL;
			}
			else
			{
				vSourceRefGene[ni] = RefGeneCreate();
				RefGeneInit_FromGenomeLabFormat(vSourceRefGene[ni], chp1, strSpecies);
			}
		}
			
		ni++;
	}
	fclose(fpRefGene);

	
	if(ni != nSourceRefNum)
	{
		printf("Error: RefFlex_LoadOriginRefGene, refgene number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* return */
	*pRefNum = nSourceRefNum;
	return vSourceRefGene;
}

/* ----------------------------------------------------------------------- */ 
/*  RefFlex_GetOrtholog():                                                 */
/*  get ortholog genes according to colinear refgene structure & alignmetn */
/* ----------------------------------------------------------------------- */ 
int RefFlex_GetOrtholog(struct tagRefGene **vOriginRefGene, int nOriginRefNum, 
			char strSourceSpecies[], char strSourceRefGenePath[],
			char strDestSpecies[], char strMapRefGenePath[],
			char strDestRefGenePath[], 
			char strOutPath[])
{
	/* define */

	/* for database */
	struct tagRefGene **vSourceRefGene;
	int nSourceRefNum;

	struct tagRefGene **vMapRefGene;
	int nMapRefNum;

	int nDestMap;
	struct tagRefGene **vDestRefGene;
	int nDestRefNum;

	int nGeneId,nSourceId;

	int ni,nj,nk,nflank;
	int nx,ny,nz,nlastid;

	/* for target search */
	FILE *fpMapOut;
	FILE *fpDestOut;
	
	int vRefSourceIdx[ORTHOLOG_COLINEAR_NUM];
	char vRefColinear[ORTHOLOG_COLINEAR_NUM][LINE_LENGTH];

	char strLine[LINE_LENGTH];
	char strRefName[LINE_LENGTH];

	FILE *fpTemp;
	FILE *fpTemp1,*fpTemp2;
	int nSourceHitNum,nMapHitNum;
	struct INTMATRIX *pHitLimit;

	int nPos1,nPos2,nPos3;
	int nMatchNum;
	int nOrientationNum;
	int nMatched;

	/* init */
	if(strcmp(strDestRefGenePath, "NULL") == 0)
	{
		nDestMap = 0;
	}
	else
	{
		nDestMap = 1;
	}

	nSourceRefNum = 0;
	nMapRefNum = 0;
	nDestRefNum = 0;
	vSourceRefGene = NULL;
	vMapRefGene = NULL;
	vDestRefGene = NULL;

	/* load source refgene */
	vSourceRefGene = RefGene_LoadDatabase(strSourceRefGenePath, 0, strSourceSpecies, &nSourceRefNum);

	/* load map refgene */
	vMapRefGene = RefGene_LoadDatabase(strMapRefGenePath, 0, strDestSpecies, &nMapRefNum);

	/* load dest refgene */
	if(nDestMap == 1)
	{
		vDestRefGene = RefGene_LoadDatabase(strDestRefGenePath, 0, strDestSpecies, &nDestRefNum);
	}


	/* locate ortholog */
	sprintf(strLine, "%s.pair", strOutPath);
	fpMapOut = NULL;
	fpMapOut = fopen(strLine, "w");
	if(fpMapOut == NULL)
	{
		printf("Error: RefFlex_GetOrtholog, cannot open map output file!\n");
		exit(EXIT_FAILURE);
	}

	if(nDestMap == 1)
	{
		sprintf(strLine, "%s.map", strOutPath);
		fpDestOut = NULL;
		fpDestOut = fopen(strLine, "w");
		if(fpDestOut == NULL)
		{
			printf("Error: RefFlex_GetOrtholog, cannot open dest output file!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* process target refgene one by one */
	fpTemp1 = NULL;
	fpTemp1 = fopen("sourcehits2.tmp", "w");
	if(fpTemp1 == NULL)
	{
		printf("Error: cannot create temp file!\n");
		exit(EXIT_FAILURE);
	}
	fpTemp2 = NULL;
	fpTemp2 = fopen("maphits2.tmp", "w");
	if(fpTemp2 == NULL)
	{
		printf("Error: cannot create temp file!\n");
		exit(EXIT_FAILURE);
	}

	for(nGeneId=0; nGeneId<nOriginRefNum; nGeneId++)
	{
		if(vOriginRefGene[nGeneId] == NULL)
		{
			fprintf(fpMapOut, ">%d\t0\n", nGeneId);
			fprintf(fpMapOut, "%s\t0\t0\t", strSourceSpecies);
			fprintf(fpMapOut, "---\n");
			fprintf(fpMapOut, "\n");

			continue;
		}

		/* look up source database */
		strcpy(strRefName, vOriginRefGene[nGeneId]->strName);
		ni = RefGene_GetBestOverlap(vOriginRefGene[nGeneId], vSourceRefGene, nSourceRefNum);
		nflank = (ORTHOLOG_COLINEAR_NUM-1)/2;

		if(ni < 0)
		{
			fprintf(fpMapOut, ">%d\t0\n", nGeneId);
			fprintf(fpMapOut, "%s\t0\t0\t", strSourceSpecies);
			fprintf(fpMapOut, "---\n");
			fprintf(fpMapOut, "\n");

			continue;
		}

		nSourceHitNum = 1;
		strcpy(strRefName, vSourceRefGene[ni]->strName);
		
		/* init flank */
		for(nk=0; nk<ORTHOLOG_COLINEAR_NUM; nk++)
		{
			vRefSourceIdx[nk] = -1;
			strcpy(vRefColinear[nk], "---");
		}

		nSourceId = ni;
		vRefSourceIdx[nflank] = ni;
		strcpy(vRefColinear[nflank], strRefName);

		/* get upstream flank */
		nk = 0;
		nlastid = ni;
		for(nj=ni-1; nj>=0; nj--)
		{
			if(vSourceRefGene[ni]->nChrom == vSourceRefGene[nj]->nChrom)
			{
				if((vSourceRefGene[ni]->nTxStart-vSourceRefGene[nj]->nTxEnd) <= ORTHOLOG_COLINEAR_WIN)
				{
					if(vSourceRefGene[nj]->nTxEnd < vSourceRefGene[nlastid]->nTxStart)
					{
						nk++;
						vRefSourceIdx[nflank-nk] = nj;
						nlastid = nj;
						strcpy(vRefColinear[nflank-nk], vSourceRefGene[nj]->strName);
					}

					if(nk == nflank)
						break;
				}
				else
				{
					break;
				}

			}
			else
			{
				break;
			}
		}

		/* get downstream flank */
		nk = 0;
		nlastid = ni;
		for(nj=ni+1; nj<nSourceRefNum; nj++)
		{
			if(vSourceRefGene[ni]->nChrom == vSourceRefGene[nj]->nChrom)
			{
				if((vSourceRefGene[nj]->nTxStart-vSourceRefGene[ni]->nTxEnd) <= ORTHOLOG_COLINEAR_WIN)
				{
					if( vSourceRefGene[nlastid]->nTxEnd < vSourceRefGene[nj]->nTxStart)
					{
						nk++;
						vRefSourceIdx[nflank+nk] = nj;
						nlastid = nj;
						strcpy(vRefColinear[nflank+nk], vSourceRefGene[nj]->strName);
					}

					if(nk == nflank)
						break;
				}
				else
				{
					break;
				}
			}
			else
			{
				break;
			}
		}
				
		fprintf(fpTemp1, "%d\t", ni);
				
		/* remove redundancy */
		for(nk=0; nk<ORTHOLOG_COLINEAR_NUM; nk++)
		{
			if(vRefSourceIdx[nk] != -1)
			{
				for(nx=0; nx<nk; nx++)
				{
					if(vRefSourceIdx[nx] != -1)
					{
						if(strcmp(vRefColinear[nk], vRefColinear[nx]) == 0)
						{
							strcpy(vRefColinear[nk], "---");
							vRefSourceIdx[nk] = -1;
							break;
						}
					}
				}
			}
			fprintf(fpTemp1, "%s\t", vRefColinear[nk]);
		}
		fprintf(fpTemp1, "\n");
		
		/* look up map database */
		nMapHitNum = 0;
		fpTemp = NULL;
		fpTemp = fopen("maphits.tmp", "w");
		if(fpTemp == NULL)
		{
			printf("Error: cannot create temp file!\n");
			exit(EXIT_FAILURE);
		}

		for(ni=0; ni<nMapRefNum; ni++)
		{
			/* if found the location in the map organism */
			if(strcmp(strRefName, vMapRefGene[ni]->strName) == 0)
			{
				nMapHitNum++;

				/* get upstream flank */
				nk = 0;
				nlastid = ni;
				for(nj=ni-1; nj>=0; nj--)
				{
					if(vMapRefGene[ni]->nChrom == vMapRefGene[nj]->nChrom)
					{
						if( (vMapRefGene[ni]->nTxStart-vMapRefGene[nj]->nTxEnd) <= ORTHOLOG_COLINEAR_WIN )
						{
							if(vMapRefGene[nj]->nTxEnd < vMapRefGene[nlastid]->nTxStart)
							{
								nk++;
								if(nk > ORTHOLOG_MAP_COLINEAR_NUM)
								{
									nlastid = nj+1;
									break;
								}
								nlastid = nj;
							}
						}
						else
						{
							nlastid = nj+1;
							break;
						}
					}
					else
					{
						nlastid = nj+1;
						break;
					}
				}

				fprintf(fpTemp, "%d\t", nlastid);
				fprintf(fpTemp, "%d\t", ni);

				fprintf(fpTemp2, "%d\t", nlastid);
				fprintf(fpTemp2, "%d\t", ni);

				/* get downstream flank */
				nk = 0;
				nlastid = ni;
				for(nj=ni+1; nj<nMapRefNum; nj++)
				{
					if(vMapRefGene[ni]->nChrom == vMapRefGene[nj]->nChrom)
					{
						if( (vMapRefGene[nj]->nTxStart-vMapRefGene[ni]->nTxEnd) <= ORTHOLOG_COLINEAR_WIN )
						{
							if(vMapRefGene[nlastid]->nTxEnd < vMapRefGene[nj]->nTxStart)
							{
								nk++;
								if(nk > ORTHOLOG_MAP_COLINEAR_NUM)
								{
									nlastid = nj-1;
									break;
								}
								nlastid = nj;
							}
						}
						else
						{
							nlastid = nj-1;
							break;
						}
					}
					else
					{
						nlastid = nj-1;
						break;
					}
				}

				fprintf(fpTemp, "%d\n", nlastid);
				fprintf(fpTemp2, "%d\n", nlastid);
			}
		}

		fclose(fpTemp);
		

		/* verify collinearity */
		if(nMapHitNum == 0)
		{
			fprintf(fpMapOut, ">%d\t0\n", nGeneId);
			fprintf(fpMapOut, "%s\t0\t0\t", strSourceSpecies);
			RefGeneWrite(vSourceRefGene[nGeneId], fpMapOut);
			fprintf(fpMapOut, "\n");
		}
		else
		{
			pHitLimit = NULL;
			pHitLimit = IMLOAD("maphits.tmp");
			if(pHitLimit == NULL)
			{
				printf("Error: cannot reload maphits.tmp!\n");
				exit(EXIT_FAILURE);
			}
			if(pHitLimit->nHeight != nMapHitNum)
			{
				printf("Error: map hit number not match!\n");
				exit(EXIT_FAILURE);
			}

			fprintf(fpMapOut, ">%d\t%d\n", nGeneId, nMapHitNum);
			fprintf(fpMapOut, "%s\t0\t0\t", strSourceSpecies);
			RefGeneWrite(vSourceRefGene[nSourceId], fpMapOut);
				
			/* match */
			for(nz=0; nz<nMapHitNum; nz++)
			{
				nMatchNum = 0;
				nOrientationNum = 0;
				nPos1 = IMGETAT(pHitLimit, nz, 0);
				nPos2 = IMGETAT(pHitLimit, nz, 1);
				nPos3 = IMGETAT(pHitLimit, nz, 2);

				for(nk=0; nk<ORTHOLOG_COLINEAR_NUM; nk++)
				{
					nMatched = 0;
					vRefSourceIdx[nk] = 0;
					for(ny=nPos1; ny<nPos2; ny++)
					{
						if(strcmp(vRefColinear[nk], vMapRefGene[ny]->strName) == 0)
						{
							nMatched = 1;
							nMatchNum++;
							vRefSourceIdx[nk] = 1;
							break;
						}
					}
					if(nMatched == 1)
						continue;

					if(strcmp(vRefColinear[nk], vMapRefGene[ny]->strName) == 0)
					{
						nMatchNum++;
						continue;
					}

					ny++;

					for(; ny<=nPos3; ny++)
					{
						if(strcmp(vRefColinear[nk], vMapRefGene[ny]->strName) == 0)
						{
							nMatched = 1;
							nMatchNum++;
							vRefSourceIdx[nk] = -1;
							break;
						}
					}
				}

				nflank = (ORTHOLOG_COLINEAR_NUM-1)/2;
				for(nk=0; nk<nflank; nk++)
				{
					for(nj=0; nj<nflank; nj++)
					{
						if(nk != nj)
						{
							if(vRefSourceIdx[nk]*vRefSourceIdx[nj] > 0)
							{
								nOrientationNum++;
							}
						}
					}
					nj++;
					for(; nj<ORTHOLOG_COLINEAR_NUM; nj++)
					{
						if(vRefSourceIdx[nk]*vRefSourceIdx[nj] < 0)
						{
							nOrientationNum++;
						}
					}
				}

				nk++;
				for(; nk<ORTHOLOG_COLINEAR_NUM; nk++)
				{
					for(nj=0; nj<nflank; nj++)
					{
						if(vRefSourceIdx[nk]*vRefSourceIdx[nj] < 0)
						{
							nOrientationNum++;
						}
					}
					nj++;
					for(; nj<ORTHOLOG_COLINEAR_NUM; nj++)
					{
						if(nj != nk)
						{
							if(vRefSourceIdx[nk]*vRefSourceIdx[nj] > 0)
							{
								nOrientationNum++;
							}
						}
					}
				}

				fprintf(fpMapOut, "%s\t%d\t%d\t", strDestSpecies, nMatchNum, nOrientationNum);
				RefGeneWrite(vMapRefGene[nPos2], fpMapOut);
			}

			fprintf(fpMapOut, "\n");
			DestroyIntMatrix(pHitLimit);
		}
	}

	fclose(fpTemp1);
	fclose(fpTemp2);

	/* delete temp files */
	RemoveFiles("maphits.tmp");
	RemoveFiles("sourcehits2.tmp");
	RemoveFiles("maphits2.tmp");
	
	/* close files */
	fclose(fpMapOut);
	if(nDestMap == 1)
	{
		fclose(fpDestOut);
	}


	/* release memory */
	RefGene_ClearDatabase(&vSourceRefGene, nSourceRefNum);
	RefGene_ClearDatabase(&vMapRefGene, nMapRefNum);
	if(nDestMap == 1)
	{
		RefGene_ClearDatabase(&vDestRefGene, nDestRefNum);
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefFlex_MergeMultiOrtholog():                                          */
/*  Merge ortholog genes from multiple pairwise ortholog mapping.          */
/*  Return the number of ortholog groups.                                  */
/* ----------------------------------------------------------------------- */ 
int RefFlex_MergeMultiOrtholog(int nSpeciesNum, struct tagString **vSpeciesName,
			char strOriginPath[], struct tagRefGene **vOriginRefGene, int nOriginRefNum,
			struct tagString **vOrthologPath, char strOutPath[])
{
	/* define */
	struct DOUBLEMATRIX *pMatchCount;
	struct DOUBLEMATRIX *pOrienCount;
	struct tagRefGene **vRefGene;
	struct tagRefGene *pRefGene;
	int nMatchCount, nOrienCount;


	/* file */
	FILE **vfpIn;
	FILE *fpOri;
	FILE *fpOut;
	int nGeneId;
	int ni,nj;
	char strRefLine0[LONG_LINE_LENGTH];
	char strRefLine1[LONG_LINE_LENGTH];

	char strLine[LINE_LENGTH];
	char strRefLine[LONG_LINE_LENGTH];
	char strSpecies[LINE_LENGTH];

	int nDestNum;
	char *chp,*chp2;

	int nED1,nED2,nLD1,nLD2;

	/* init */
	
	/* open files */
	if(nSpeciesNum <= 1)
	{
		printf("Warning: RefFlex_MergeMultiOrtholog, no files to merge!\n");
		return PROC_FAILURE;
	}

	vfpIn = NULL;
	vfpIn = (FILE **)calloc(nSpeciesNum, sizeof(FILE *));
	if(vfpIn == NULL)
	{
		printf("Error: RefFlex_MergeMultiOrtholog, cannot match ortholog!\n");
		exit(EXIT_FAILURE);
	}

	fpOri = fopen(strOriginPath, "r");
	if(fpOri == NULL)
	{
		printf("Error: RefFlex_MergeMultiOrtholog, cannot open file!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nSpeciesNum; ni++)
	{
		sprintf(strLine, "%s.pair", vOrthologPath[ni]->m_pString);
		vfpIn[ni] = fopen(strLine, "r");
		if(vfpIn[ni] == NULL)
		{
			printf("Error: RefFlex_MergeMultiOrtholog, cannot open file!\n");
			exit(EXIT_FAILURE);
		}
	}

	fpOut = NULL;
	fpOut = fopen(strOutPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: RefFlex_MergeMultiOrtholog, cannot open file!\n");
		exit(EXIT_FAILURE);
	}

	/* get the best ortholog */
	nGeneId = 0;
	while(fgets(strRefLine0, LONG_LINE_LENGTH, fpOri) != NULL)
	{
		StrTrimLeft(strRefLine0);
		StrTrimRight(strRefLine0);
		if(strRefLine0[0] == '\0')
			continue;

		/* output the first line */
		fprintf(fpOut, ">%d\n", nGeneId);
		fprintf(fpOut, "%s\n", strRefLine0);

		/* initialize ortholog settings */
		pMatchCount = NULL;
		pMatchCount = CreateDoubleMatrix(1, nSpeciesNum);
		if(pMatchCount == NULL)
		{
			printf("Error: RefFlex_MergeMultiOrtholog, cannot create matchcount matrix!\n");
			exit(EXIT_FAILURE);
		}
		pOrienCount = NULL;
		pOrienCount = CreateDoubleMatrix(1, nSpeciesNum);
		if(pOrienCount == NULL)
		{
			printf("Error: RefFlex_MergeMultiOrtholog, cannot create oriencount matrix!\n");
			exit(EXIT_FAILURE);
		}
		vRefGene = NULL;
		vRefGene = (struct tagRefGene **)calloc(nSpeciesNum, sizeof(struct tagRefGene*));
		if(vRefGene == NULL)
		{
			printf("Error: RefFlex_MergeMultiOrtholog, cannot create refgene vector!\n");
			exit(EXIT_FAILURE);
		}
			
		/* process all files */
		for(ni=0; ni<nSpeciesNum; ni++)
		{
			/* synchronize */
			while(fgets(strRefLine1, LONG_LINE_LENGTH, vfpIn[ni]) != NULL)
			{
				StrTrimLeft(strRefLine1);
				StrTrimRight(strRefLine1);
				if(strRefLine1[0] == '\0')
					continue;

				if(strRefLine1[0] == '>')
				{
					break;
				}
			}

			/* load head info */
			chp = strRefLine1+1;
			sscanf(chp, "%d %d", &nj, &nDestNum);
			if(nj != nGeneId)
			{
				printf("Error: RefFlex_MergeMultiOrtholog, file not synchronized!\n");
				exit(EXIT_FAILURE);
			}

			/* load the first sequence */
			fgets(strRefLine1, LONG_LINE_LENGTH, vfpIn[ni]);
						
			/* load ortholog sequences */
			for(nj=0; nj<nDestNum; nj++)
			{
				fgets(strRefLine1, LONG_LINE_LENGTH, vfpIn[ni]);
				StrTrimLeft(strRefLine1);
				StrTrimRight(strRefLine1);
				chp = strRefLine1;
				chp2 = strchr(chp, '\t');
				*chp2 = '\0';
				strcpy(strSpecies, chp);
				chp = chp2+1;
				chp2 = strchr(chp, '\t');
				*chp2 = '\0';
				nMatchCount = atoi(chp);
				chp = chp2+1;
				chp2 = strchr(chp, '\t');
				*chp2 = '\0';
				nOrienCount = atoi(chp);
				chp = chp2+1;
				if(strcmp(strSpecies, vSpeciesName[ni]->m_pString) != 0)
				{
					printf("Error: RefFlex_MergeMultiOrtholog, species not match!\n");
					exit(EXIT_FAILURE);
				}
				sscanf(chp, "%s", strRefLine);
				if(strcmp(strRefLine, "---") == 0)
				{
					pRefGene = NULL;
				}
				else
				{
					strcpy(strRefLine, chp);
					pRefGene = NULL;
					pRefGene = RefGeneCreate();
					RefGeneInit_FromGenomeLabFormat(pRefGene, strRefLine, strSpecies);
				}
				
				if(pRefGene != NULL)
				{
					if((nMatchCount+ORTHOLOG_COMBCOUNT_COEF*nOrienCount) > (pMatchCount->pMatElement[ni]+ORTHOLOG_COMBCOUNT_COEF*pOrienCount->pMatElement[ni]))
					{
						pMatchCount->pMatElement[ni] = nMatchCount;
						pOrienCount->pMatElement[ni] = nOrienCount;
						
						if(vRefGene[ni] != NULL)
						{
							RefGeneDestroy(vRefGene[ni]);
						}
						vRefGene[ni] = pRefGene;
					}
					else if( fabs( (nMatchCount+ORTHOLOG_COMBCOUNT_COEF*nOrienCount) - (pMatchCount->pMatElement[ni]+ORTHOLOG_COMBCOUNT_COEF*pOrienCount->pMatElement[ni]) ) < ZERO_BOUND )
					{
						if(vRefGene[ni] == NULL)
						{
							pMatchCount->pMatElement[ni] = nMatchCount;
							pOrienCount->pMatElement[ni] = nOrienCount;
							vRefGene[ni] = pRefGene;
						}
						else
						{
							nED1 = (int)fabs(vOriginRefGene[nGeneId]->nExonCount-pRefGene->nExonCount);
							nED2 = (int)fabs(vOriginRefGene[nGeneId]->nExonCount-vRefGene[ni]->nExonCount);
							nLD1 = (int)fabs(pRefGene->nTxEnd-pRefGene->nTxStart-vOriginRefGene[nGeneId]->nTxEnd+vOriginRefGene[nGeneId]->nTxStart);
							nLD2 = (int)fabs(vRefGene[ni]->nTxEnd-vRefGene[ni]->nTxStart-vOriginRefGene[nGeneId]->nTxEnd+vOriginRefGene[nGeneId]->nTxStart);
							if( nED1 < nED2 )
							{
								pMatchCount->pMatElement[ni] = nMatchCount;
								pOrienCount->pMatElement[ni] = nOrienCount;
								RefGeneDestroy(vRefGene[ni]);
								vRefGene[ni] = pRefGene;
							}
							else if( nLD1 < nLD2 )
							{
								pMatchCount->pMatElement[ni] = nMatchCount;
								pOrienCount->pMatElement[ni] = nOrienCount;
								RefGeneDestroy(vRefGene[ni]);
								vRefGene[ni] = pRefGene;
							}
							else
							{
								RefGeneDestroy(pRefGene);
							}
						}
					}
					else
					{
						RefGeneDestroy(pRefGene);
					}
				}
			}
		}
		
		for(ni=0; ni<nSpeciesNum; ni++)
		{
			if((pMatchCount != NULL) && (pOrienCount != NULL))
			{
				fprintf(fpOut, "%s\t%d\t%d\t", vSpeciesName[ni]->m_pString, (int)(pMatchCount->pMatElement[ni]),
					(int)(pOrienCount->pMatElement[ni]));
			}
			else
			{
				fprintf(fpOut, "%s\t0\t0\t", vSpeciesName[ni]->m_pString);
			}
			
			if(vRefGene[ni] != NULL)
			{
				RefGeneWrite(vRefGene[ni], fpOut);
			}
			else
			{
				fprintf(fpOut, "---\n");
			}
		}
		fprintf(fpOut, "\n");
		nGeneId++;

		DestroyDoubleMatrix(pMatchCount);
		DestroyDoubleMatrix(pOrienCount);
		for(ni=0; ni<nSpeciesNum; ni++)
		{
			RefGeneDestroy(vRefGene[ni]);
			vRefGene[ni] = NULL;
		}
		free(vRefGene);
	}

	/* close files */
	for(ni=0; ni<nSpeciesNum; ni++)
	{
		fclose(vfpIn[ni]);
		vfpIn[ni] = NULL;
	}
	free(vfpIn);
	fclose(fpOut);
	fclose(fpOri);

	if(nOriginRefNum != nGeneId)
	{
		printf("Error: RefFlex_MergeMultiOrtholog, gene numbers do not match!\n");
		exit(EXIT_FAILURE);
	}

	/* return */
	return nGeneId;
}

/* ----------------------------------------------------------------------- */ 
/*  RefFlex_MapMultiOrthologToOriginalSpecies()                            */
/*  Map refgene orthologs back to the original species.                    */
/* ----------------------------------------------------------------------- */ 
int RefFlex_MapMultiOrthologToOriginalSpecies(int nSpeciesNum, struct tagString **vSpeciesName, 
			struct tagString **vSpeciesMapbase, struct tagString **vSpeciesDatabase,
			struct tagRefGene **vOriginRefGene, int nOriginRefNum,
			char strInPath[], char strOutPath[])
{
	/* define */
	struct tagRefGene ***vvDatabase;
	struct INTMATRIX *pDataNum;

	int ni,nj,nk;
	char *chp,*chp2;
	char strRefLine[LONG_LINE_LENGTH];
	char strRefLine0[LONG_LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	int nMatchCount,nOrienCount;
	struct tagRefGene *pRefGene;
	
	FILE *fpIn;
	FILE *fpOut;

	/* init */
	if(nSpeciesNum <= 0)
	{
		printf("Warning: no species available!\n");
		return PROC_FAILURE;
	}
	if(nOriginRefNum <= 0)
	{
		printf("Warning: no genes available!\n");
		return PROC_FAILURE;
	}


	/* load database */
	vvDatabase = NULL;
	vvDatabase = (struct tagRefGene ***)calloc(nSpeciesNum, sizeof(struct tagRefGene**));
	if(vvDatabase == NULL)
	{
		printf("Error: RefFlex_MapMultiOrthologToOriginalSpecies, cannot create species database!\n");
		exit(EXIT_FAILURE);
	}

	pDataNum = NULL;
	pDataNum = CreateIntMatrix(1, nSpeciesNum);
	if(pDataNum == NULL)
	{
		printf("Error: RefFlex_MapMultiOrthologToOriginalSpecies, cannot create species database!\n");
		exit(EXIT_FAILURE);
	}

	for(nj=0; nj<nSpeciesNum; nj++)
	{
		/* if no need to map, then keep the original one */
		if(strcmp(vSpeciesDatabase[nj]->m_pString, "NULL") == 0)
			continue;

		/* if need map, create database first */
		strcpy(strSpecies, vSpeciesName[nj]->m_pString);
		vvDatabase[nj] = RefGene_LoadDatabase(vSpeciesDatabase[nj]->m_pString, 0, 
			strSpecies, pDataNum->pMatElement+nj);

		if(vvDatabase[nj] == NULL)
		{
			printf("Warning: RefFlex_MapMultiOrthologToOriginalSpecies, cannot open dest refgene database file!\n");
		}
	}


	/* open files */
	fpIn = NULL;
	fpIn = fopen(strInPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: RefFlex_MapMultiOrthologToOriginalSpecies, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: RefFlex_MapMultiOrthologToOriginalSpecies, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}


	/* process one by one */
	ni = 0;
	while(fgets(strRefLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strRefLine);
		StrTrimRight(strRefLine);
		if(strRefLine[0] == '\0')
			continue;

		if(strRefLine[0] == '>')
		{
			/* head info */
			fprintf(fpOut, "%s\n", strRefLine);
			
			/* first line */
			fgets(strRefLine, LONG_LINE_LENGTH, fpIn);
			StrTrimLeft(strRefLine);
			StrTrimRight(strRefLine);
			fprintf(fpOut, "%s\n", strRefLine);

			/* load species */
			for(nj=0; nj<nSpeciesNum; nj++)
			{
				fgets(strRefLine, LONG_LINE_LENGTH, fpIn);
				StrTrimLeft(strRefLine);
				StrTrimRight(strRefLine);
				
				chp = strRefLine;
				chp2 = strchr(chp, '\t');
				*chp2 = '\0';
				strcpy(strSpecies, chp);
				chp = chp2+1;
				chp2 = strchr(chp, '\t');
				*chp2 = '\0';
				nMatchCount = atoi(chp);
				chp = chp2+1;
				chp2 = strchr(chp, '\t');
				*chp2 = '\0';
				nOrienCount = atoi(chp);
				fprintf(fpOut, "%s\t%d\t%d\t", strSpecies, nMatchCount, nOrienCount);

				chp = chp2+1;
				strcpy(strRefLine0, chp);
				
				if(strcmp(strSpecies, vSpeciesName[nj]->m_pString) != 0)
				{
					printf("Error: RefFlex_MapMultiOrthologToOriginalSpecies, species not match!\n");
					exit(EXIT_FAILURE);
				}
				if(strstr(strRefLine0, "---") == strRefLine0)
				{
					fprintf(fpOut, "---\n");
				}
				else
				{
					pRefGene = NULL;
					pRefGene = RefGeneCreate();
					RefGeneInit_FromGenomeLabFormat(pRefGene, strRefLine0, strSpecies);
					
					if(vvDatabase[nj] == NULL)
					{
						RefGeneWrite(pRefGene, fpOut);
					}
					else
					{
						nk = RefGene_GetBestOverlap(pRefGene, vvDatabase[nj], pDataNum->pMatElement[nj]);

						if(nk >= 0)
						{
							RefGeneWrite(vvDatabase[nj][nk], fpOut);
						}
						else
						{
							RefGeneWrite(pRefGene, fpOut);
						}
					}

					RefGeneDestroy(pRefGene);
				}
			}

			fprintf(fpOut, "\n");
			ni++;
		}
	}

	if(ni != nOriginRefNum)
	{
		printf("Error: RefFlex_MapMultiOrthologToOriginalSpecies, loading error!\n");
		exit(EXIT_FAILURE);
	}

	/* close files */
	fclose(fpIn);
	fclose(fpOut);

	/* release memory */
	for(nj=0; nj<nSpeciesNum; nj++)
	{
		if(vvDatabase[nj] != NULL)
		{
			RefGene_ClearDatabase(vvDatabase+nj, pDataNum->pMatElement[nj]);
		}
	}
	free(vvDatabase);
	DestroyIntMatrix(pDataNum);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefFlex_GetMultiOrtholog1way_Main()                                    */
/*  Get ortholog genes from 1way refgene mapping.                          */
/* ----------------------------------------------------------------------- */ 
int RefFlex_GetMultiOrtholog1way_Main(char strTargetPath[], char strParamPath[], 
								  char strOutPath[])
{
	/* define */
	int nMapSpeciesNum = 0;
	struct tagString **vMapSpeciesName;
	struct tagString **vRefSpeciesName;
	struct tagString **vRefSpeciesMapbase;
	struct tagString **vMapSpeciesMapbase;
	struct tagString **vMapSpeciesDatabase;

	FILE *fpIn;
	char strLine[LONG_LINE_LENGTH];
	char strInFileName[MED_LINE_LENGTH];
	char strOutFileName[MED_LINE_LENGTH];

	/* others */
	char *chp;
	int ni;
	
	/* init */
	vMapSpeciesName = NULL;
	vRefSpeciesName = NULL;
	vRefSpeciesMapbase = NULL;
	vMapSpeciesMapbase = NULL;
	vMapSpeciesDatabase = NULL;

	/* load parameters */
	fpIn = NULL;
	fpIn = fopen(strParamPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: cannot load parameters!\n");
		exit(EXIT_FAILURE);
	}

	/* load line by line */
	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strstr(strLine, "[Map Species Number]") == strLine)
		{
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);
			nMapSpeciesNum = atoi(chp);
			if(nMapSpeciesNum <= 0)
			{
				printf("Error: RefFlex_GetMultiOrtholog1way_Main, map species number <= 0!\n");
				exit(EXIT_FAILURE);
			}

			vRefSpeciesName = (struct tagString **)calloc(nMapSpeciesNum, sizeof(struct tagString*));
			if(vRefSpeciesName == NULL)
			{
				printf("Error: RefFlex_GetMultiOrtholog1way_Main, cannot create database\n");
				exit(EXIT_FAILURE);
			}

			vMapSpeciesName = (struct tagString **)calloc(nMapSpeciesNum, sizeof(struct tagString*));
			if(vMapSpeciesName == NULL)
			{
				printf("Error: RefFlex_GetMultiOrtholog1way_Main, cannot create database\n");
				exit(EXIT_FAILURE);
			}

			vRefSpeciesMapbase = (struct tagString **)calloc(nMapSpeciesNum, sizeof(struct tagString*));
			if(vRefSpeciesMapbase == NULL)
			{
				printf("Error: RefFlex_GetMultiOrtholog1way_Main, cannot create database\n");
				exit(EXIT_FAILURE);
			}

			vMapSpeciesDatabase = (struct tagString **)calloc(nMapSpeciesNum, sizeof(struct tagString*));
			if(vMapSpeciesDatabase == NULL)
			{
				printf("Error: RefFlex_GetMultiOrtholog1way_Main, cannot create database\n");
				exit(EXIT_FAILURE);
			}

			vMapSpeciesMapbase = (struct tagString **)calloc(nMapSpeciesNum, sizeof(struct tagString*));
			if(vMapSpeciesMapbase == NULL)
			{
				printf("Error: RefFlex_GetMultiOrtholog1way_Main, cannot create database\n");
				exit(EXIT_FAILURE);
			}
		}
		else if(strstr(strLine, "[Map Species Name]") == strLine)
		{
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);
			StringAddTail(vMapSpeciesName+ni, chp);

			/* load ref species name */
			fgets(strLine, LONG_LINE_LENGTH, fpIn);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "[Ref Species Name]") != strLine)
			{
				printf("Error: RefFlex_GetMultiOrtholog1way_Main, parameter wrong\n");
				exit(EXIT_FAILURE);
			}
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);
			StringAddTail(vRefSpeciesName+ni, chp);

			/* load ref species mapbase */
			fgets(strLine, LONG_LINE_LENGTH, fpIn);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "[Ref Species RefGeneMap]") != strLine)
			{
				printf("Error: RefFlex_GetMultiOrtholog1way_Main, parameter wrong\n");
				exit(EXIT_FAILURE);
			}
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);
			StringAddTail(vRefSpeciesMapbase+ni, chp);

			/* load map species mapbase */
			fgets(strLine, LONG_LINE_LENGTH, fpIn);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "[Map Species RefGeneMap]") != strLine)
			{
				printf("Error: RefFlex_GetMultiOrtholog1way_Main, parameter wrong\n");
				exit(EXIT_FAILURE);
			}
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);
			StringAddTail(vMapSpeciesMapbase+ni, chp);

			/* load map species database */
			fgets(strLine, LONG_LINE_LENGTH, fpIn);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "[Map Species RefGene]") != strLine)
			{
				printf("Error: RefFlex_GetMultiOrtholog1way_Main, parameter wrong\n");
				exit(EXIT_FAILURE);
			}
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);
			StringAddTail(vMapSpeciesDatabase+ni, chp);

			/* add the species number */
			ni++;
		}
		else
		{
			printf("Error: RefFlex_GetMultiOrtholog1way_Main, unknown parameters\n");
			exit(EXIT_FAILURE);
		}
	}

	if(ni != nMapSpeciesNum)
	{
		printf("Error: RefFlex_GetMultiOrtholog1way_Main, species number not match\n");
		exit(EXIT_FAILURE);
	}

	/* close files */
	fclose(fpIn);

	/* get orthologs from species to species */
	if(nMapSpeciesNum == 1)
	{
		RefFlex_GetMultiOrtholog1way(strTargetPath, strOutPath,
			vRefSpeciesName[0]->m_pString, vMapSpeciesName[0]->m_pString,
			vRefSpeciesMapbase[0]->m_pString, vMapSpeciesMapbase[0]->m_pString,
			vMapSpeciesDatabase[0]->m_pString);
	}
	else
	{
		sprintf(strOutFileName, "%s0.tmp", strOutPath);
		RefFlex_GetMultiOrtholog1way(strTargetPath, strOutFileName,
			vRefSpeciesName[0]->m_pString, vMapSpeciesName[0]->m_pString,
			vRefSpeciesMapbase[0]->m_pString, vMapSpeciesMapbase[0]->m_pString,
			vMapSpeciesDatabase[0]->m_pString);
		
		for(ni=1; ni<(nMapSpeciesNum-1); ni++)
		{
			sprintf(strInFileName, "%s%d.tmp", strOutPath, (ni-1));
			sprintf(strOutFileName, "%s%d.tmp", strOutPath, ni);
			RefFlex_GetMultiOrtholog1way(strInFileName, strOutFileName,
				vRefSpeciesName[ni]->m_pString, vMapSpeciesName[ni]->m_pString,
				vRefSpeciesMapbase[ni]->m_pString, vMapSpeciesMapbase[ni]->m_pString,
				vMapSpeciesDatabase[ni]->m_pString);
		}

		sprintf(strInFileName, "%s%d.tmp", strOutPath, (ni-1));
		RefFlex_GetMultiOrtholog1way(strInFileName, strOutPath,
			vRefSpeciesName[ni]->m_pString, vMapSpeciesName[ni]->m_pString,
			vRefSpeciesMapbase[ni]->m_pString, vMapSpeciesMapbase[ni]->m_pString,
			vMapSpeciesDatabase[ni]->m_pString);

		/* delete temp files */
		if(strcmp(OS_SYSTEM, "UNIX") == 0)
		{
			sprintf(strLine, "rm %s*.tmp", strOutPath);
			system(strLine);
		}
		else
		{
			sprintf(strLine, "del %s*.tmp", strOutPath);
			system(strLine);
		}
	}

	/* release memory */
	for(ni=0; ni<nMapSpeciesNum; ni++)
	{
		DeleteString(vRefSpeciesName[ni]);
		DeleteString(vMapSpeciesName[ni]);
		DeleteString(vRefSpeciesMapbase[ni]);
		DeleteString(vMapSpeciesMapbase[ni]);
		DeleteString(vMapSpeciesDatabase[ni]);
	}
	free(vRefSpeciesName);
	free(vMapSpeciesName);
	free(vRefSpeciesMapbase);
	free(vMapSpeciesMapbase);
	free(vMapSpeciesDatabase);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefFlex_GetMultiOrtholog1way()                                         */
/*  Get ortholog genes from transcript mapping. This mapping does not      */
/*  require colinearity as getmultiortholog and is therefore much less     */
/*  stringent.                                                             */
/* ----------------------------------------------------------------------- */ 
int RefFlex_GetMultiOrtholog1way(char strInPath[], char strOutPath[],
			char strRefSpeciesName[], char strMapSpeciesName[],
			char strRefSpeciesMapbase[], char strMapSpeciesMapbase[],
			char strMapSpeciesDatabase[])
{
	/* define */
	FILE *fpSourceRefGene;
	struct tagRefGene *pSourceRefGeneList;
	struct tagRefGene **vSourceRefGene;
	int nSourceRefNum;

	FILE *fpMapRefGene;
	struct tagRefGene *pMapRefGeneList;
	struct tagRefGene **vMapRefGene;
	int nMapRefNum;

	int nDestMap;
	FILE *fpDestRefGene;
	struct tagRefGene *pDestRefGeneList;
	struct tagRefGene **vDestRefGene;
	int nDestRefNum;

	struct tagRefGene *pRefGene, *pCurrentRefGene;
	char strRefLine[LONG_LINE_LENGTH];
	int ni;

	/* for target search */
	FILE *fpIn;
	FILE *fpOut;
	int nFirstRefGene = 1;
	char strSpecies[LINE_LENGTH];
	int nSeedMatchId,nOrthoMatchId,nDestMatchId;
	struct tagRefGene *pSeedRefGene,*pRefMapRefGene,*pMapMapRefGene,*pMapDataRefGene;
	double dOptOrienCount,dOptMatchCount;
	
	char *chp,*chp2;


	/* ----------------------------- */
	/* init databases                */
	/* ----------------------------- */
	if(strcmp(strMapSpeciesDatabase, "NULL") == 0)
	{
		nDestMap = 0;
	}
	else
	{
		nDestMap = 1;
	}

	pSourceRefGeneList = NULL;
	pMapRefGeneList = NULL;
	pDestRefGeneList = NULL;
	nSourceRefNum = 0;
	nMapRefNum = 0;
	nDestRefNum = 0;
	vSourceRefGene = NULL;
	vMapRefGene = NULL;
	vDestRefGene = NULL;

	/* load source refgene */
	fpSourceRefGene = NULL;
	fpSourceRefGene = fopen(strRefSpeciesMapbase, "r");
	if(fpSourceRefGene == NULL)
	{
		printf("Error: RefFlex_GetMultiOrtholog1way, cannot open source refgene file!\n");
		exit(EXIT_FAILURE);
	}

	pCurrentRefGene = NULL;
	while(fgets(strRefLine, LONG_LINE_LENGTH, fpSourceRefGene) != NULL)
	{
		StrTrimLeft(strRefLine);
		StrTrimRight(strRefLine);
		if(strRefLine[0] == '\0')
			continue;

		pRefGene = NULL;
		pRefGene = RefGeneCreate();
		RefGeneInit_FromGenomeLabFormat(pRefGene, strRefLine, strRefSpeciesName);
		if(pSourceRefGeneList == NULL)
		{
			pSourceRefGeneList = pRefGene;
			pCurrentRefGene = pRefGene;
		}
		else
		{
			pCurrentRefGene->pNext = pRefGene;
			pCurrentRefGene = pRefGene;
		}
		nSourceRefNum++;
	}

	fclose(fpSourceRefGene);

	vSourceRefGene = (struct tagRefGene **)calloc(nSourceRefNum, sizeof(struct tagRefGene*));
	if(vSourceRefGene == NULL)
	{
		printf("Error: RefFlex_GetMultiOrtholog1way, organize refgene into array!\n");
		exit(EXIT_FAILURE);
	}

	ni=0;
	while(pSourceRefGeneList != NULL)
	{
		pRefGene = pSourceRefGeneList;
		pSourceRefGeneList = pRefGene->pNext;
		pRefGene->pNext = NULL;
		vSourceRefGene[ni] = pRefGene;
		ni++;
	}

	if(ni != nSourceRefNum)
	{
		printf("Error: RefFlex_GetMultiOrtholog1way, refgene number not match!\n");
		exit(EXIT_FAILURE);
	}


	/* load map refgene */
	fpMapRefGene = NULL;
	fpMapRefGene = fopen(strMapSpeciesMapbase, "r");
	if(fpMapRefGene == NULL)
	{
		printf("Error: RefFlex_GetMultiOrtholog1way, cannot open map refgene file!\n");
		exit(EXIT_FAILURE);
	}

	pCurrentRefGene = NULL;
	while(fgets(strRefLine, LONG_LINE_LENGTH, fpMapRefGene) != NULL)
	{
		StrTrimLeft(strRefLine);
		StrTrimRight(strRefLine);
		if(strRefLine[0] == '\0')
			continue;

		pRefGene = NULL;
		pRefGene = RefGeneCreate();
		RefGeneInit_FromGenomeLabFormat(pRefGene, strRefLine, strMapSpeciesName);
		if(pMapRefGeneList == NULL)
		{
			pMapRefGeneList = pRefGene;
			pCurrentRefGene = pRefGene;
		}
		else
		{
			pCurrentRefGene->pNext = pRefGene;
			pCurrentRefGene = pRefGene;
		}
		nMapRefNum++;
	}

	fclose(fpMapRefGene);

	vMapRefGene = (struct tagRefGene **)calloc(nMapRefNum, sizeof(struct tagRefGene*));
	if(vMapRefGene == NULL)
	{
		printf("Error: RefFlex_GetMultiOrtholog1way, cannot organize refgene into array!\n");
		exit(EXIT_FAILURE);
	}

	ni=0;
	while(pMapRefGeneList != NULL)
	{
		pRefGene = pMapRefGeneList;
		pMapRefGeneList = pRefGene->pNext;
		pRefGene->pNext = NULL;
		vMapRefGene[ni] = pRefGene;
		ni++;
	}

	if(ni != nMapRefNum)
	{
		printf("Error: RefFlex_GetMultiOrtholog1way, refgene number not match!\n");
		exit(EXIT_FAILURE);
	}


	/* load dest refgene */
	if(nDestMap == 1)
	{
		fpDestRefGene = NULL;
		fpDestRefGene = fopen(strMapSpeciesDatabase, "r");
		if(fpDestRefGene == NULL)
		{
			printf("Error: RefFlex_GetMultiOrtholog1way, cannot open dest refgene file!\n");
			exit(EXIT_FAILURE);
		}

		pCurrentRefGene = NULL;
		while(fgets(strRefLine, LONG_LINE_LENGTH, fpDestRefGene) != NULL)
		{
			StrTrimLeft(strRefLine);
			StrTrimRight(strRefLine);
			if(strRefLine[0] == '\0')
				continue;

			pRefGene = NULL;
			pRefGene = RefGeneCreate();
			RefGeneInit_FromGenomeLabFormat(pRefGene, strRefLine, strMapSpeciesName);
			if(pDestRefGeneList == NULL)
			{
				pDestRefGeneList = pRefGene;
				pCurrentRefGene = pRefGene;
			}
			else
			{
				pCurrentRefGene->pNext = pRefGene;
				pCurrentRefGene = pRefGene;
			}
			nDestRefNum++;
		}

		fclose(fpDestRefGene);

		vDestRefGene = (struct tagRefGene **)calloc(nDestRefNum, sizeof(struct tagRefGene*));
		if(vDestRefGene == NULL)
		{
			printf("Error: RefFlex_GetMultiOrtholog1way, organize refgene into array!\n");
			exit(EXIT_FAILURE);
		}

		ni=0;
		while(pDestRefGeneList != NULL)
		{
			pRefGene = pDestRefGeneList;
			pDestRefGeneList = pRefGene->pNext;
			pRefGene->pNext = NULL;
			vDestRefGene[ni] = pRefGene;
			ni++;
		}

		if(ni != nDestRefNum)
		{
			printf("Error: RefFlex_GetMultiOrtholog1way, refgene number not match!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* ----------------------------- */
	/* locate ortholog               */
	/* ----------------------------- */
	pSeedRefGene = NULL;
	pRefMapRefGene = NULL;
	pMapMapRefGene = NULL;
	pMapDataRefGene = NULL;

	fpIn = NULL;
	fpIn = fopen(strInPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: RefFlex_GetMultiOrtholog1way, cannot open input map file!\n");
		exit(EXIT_FAILURE);
	}
	fpOut = NULL;
	fpOut = fopen(strOutPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: RefFlex_GetMultiOrtholog1way, cannot open output map file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strRefLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strRefLine);
		StrTrimRight(strRefLine);
		if(strRefLine[0] == '\0')
			continue;

		if(strRefLine[0] == '>')
		{
			if(pSeedRefGene == NULL)
			{
				if(nFirstRefGene == 1)
				{
					nFirstRefGene = 0;
				}
				else
				{
					fprintf(fpOut, "%s\t0\t0\t---\n", strMapSpeciesName);
					fprintf(fpOut, "\n");
				}
			}
			else
			{
				/* --------------------- */
				/* step I: --> refmap    */
				/* --------------------- */
				nSeedMatchId = RefGene_GetBestOverlap(pSeedRefGene, vSourceRefGene, nSourceRefNum);
				if(nSeedMatchId < 0)
				{
					pRefMapRefGene = NULL;
				}
				else
				{
					pRefMapRefGene = vSourceRefGene[nSeedMatchId];
				}

				/* --------------------- */
				/* step II: --> mapmap   */
				/* --------------------- */
				if(nSeedMatchId >= 0)
				{
					nOrthoMatchId = RefGene_GetSeededRelaxedOrtholog(nSeedMatchId, vSourceRefGene, nSourceRefNum, 
						vMapRefGene, nMapRefNum, &dOptOrienCount, &dOptMatchCount);
					if(nOrthoMatchId < 0)
					{
						pMapMapRefGene = NULL;
					}
					else
					{
						pMapMapRefGene = vMapRefGene[nOrthoMatchId];
					}
				}

				/* --------------------- */
				/* step III: --> mapdata */
				/* --------------------- */
				if((nDestMap == 1) && (nOrthoMatchId >= 0))
				{
					nDestMatchId = RefGene_GetBestOverlap(pMapMapRefGene, vDestRefGene, nDestRefNum);
					if(nDestMatchId < 0)
					{
						pMapDataRefGene = NULL;
					}
					else
					{
						pMapDataRefGene = vDestRefGene[nDestMatchId];
					}
				}

				/* --------------------- */
				/* step IV: --> export   */
				/* --------------------- */
				fprintf(fpOut, "%s\t%d\t%d\t", strMapSpeciesName, (int)dOptMatchCount,
					(int)dOptOrienCount);
				if( (nDestMap == 1) && (pMapDataRefGene != NULL) )
				{
					RefGeneWrite(pMapDataRefGene, fpOut);
				}
				else
				{
					if(pMapMapRefGene == NULL)
					{
						fprintf(fpOut, "---\n");
					}
					else
					{
						RefGeneWrite(pMapMapRefGene, fpOut);
					}
				}

				fprintf(fpOut, "\n");
			}

			/* --------------------- */
			/* step V: --> reinit    */
			/* --------------------- */
			nSeedMatchId = -1;
			nOrthoMatchId = -1;
			nDestMatchId = -1;
			RefGeneDestroy(pSeedRefGene);
			pSeedRefGene = NULL;
			pRefMapRefGene = NULL;
			pMapMapRefGene = NULL;
			pMapDataRefGene = NULL;
			dOptOrienCount = 0.0;
			dOptMatchCount = 0.0;

			fprintf(fpOut, "%s\n", strRefLine);
			while(fgets(strRefLine, LONG_LINE_LENGTH, fpIn) != NULL)
			{
				StrTrimLeft(strRefLine);
				StrTrimRight(strRefLine);
				if(strRefLine[0] != '\0')
					break;
			}
		}
		else
		{
			/* ------------------ */
			/* get ref species id */
			/* ------------------ */
			sscanf(strRefLine, "%s", strSpecies);
			if(strcmp(strSpecies, strRefSpeciesName) == 0)
			{
				chp = strchr(strRefLine, '\t');
				chp++;
				
				chp2 = strchr(chp, '\t');
				chp = chp2+1;
				chp2 = strchr(chp, '\t');
				chp = chp2+1;

				if(strcmp(chp, "---") != 0)
				{
					if(pSeedRefGene != NULL)
					{
						printf("Error: RefFlex_GetMultiOrtholog1way, there are more than one seed refgene!\n");
						exit(EXIT_FAILURE);
					}

					pSeedRefGene = RefGeneCreate();
					RefGeneInit_FromGenomeLabFormat(pSeedRefGene, chp, strRefSpeciesName);
				}
			}
		}

		fprintf(fpOut, "%s\n", strRefLine);
	}

	/* -------------------- */
	/* process the last one */
	/* -------------------- */
	if(pSeedRefGene == NULL)
	{
		if(nFirstRefGene == 1)
		{
			nFirstRefGene = 0;
		}
		else
		{
			fprintf(fpOut, "%s\t0\t0\t---\n", strMapSpeciesName);
			fprintf(fpOut, "\n");
		}
	}
	else
	{
		/* --------------------- */
		/* step I: --> refmap    */
		/* --------------------- */
		nSeedMatchId = RefGene_GetBestOverlap(pSeedRefGene, vSourceRefGene, nSourceRefNum);
		if(nSeedMatchId < 0)
		{
			pRefMapRefGene = NULL;
		}
		else
		{
			pRefMapRefGene = vSourceRefGene[nSeedMatchId];
		}

		/* --------------------- */
		/* step II: --> mapmap   */
		/* --------------------- */
		if(nSeedMatchId >= 0)
		{
			nOrthoMatchId = RefGene_GetSeededRelaxedOrtholog(nSeedMatchId, vSourceRefGene, nSourceRefNum, 
				vMapRefGene, nMapRefNum, &dOptOrienCount, &dOptMatchCount);
			if(nOrthoMatchId < 0)
			{
				pMapMapRefGene = NULL;
			}
			else
			{
				pMapMapRefGene = vMapRefGene[nOrthoMatchId];
			}
		}

		/* --------------------- */
		/* step III: --> mapdata */
		/* --------------------- */
		if((nDestMap == 1) && (nOrthoMatchId >= 0))
		{
			nDestMatchId = RefGene_GetBestOverlap(pMapMapRefGene, vDestRefGene, nDestRefNum);
			if(nDestMatchId < 0)
			{
				pMapDataRefGene = NULL;
			}
			else
			{
				pMapDataRefGene = vDestRefGene[nDestMatchId];
			}
		}

		/* --------------------- */
		/* step IV: --> export   */
		/* --------------------- */
		fprintf(fpOut, "%s\t%d\t%d\t", strMapSpeciesName, (int)dOptMatchCount,
			(int)dOptOrienCount);
		if( (nDestMap == 1) && (pMapDataRefGene != NULL) )
		{
			RefGeneWrite(pMapDataRefGene, fpOut);
		}
		else
		{
			if(pMapMapRefGene == NULL)
			{
				fprintf(fpOut, "---\n");
			}
			else
			{
				RefGeneWrite(pMapMapRefGene, fpOut);
			}
		}

		fprintf(fpOut, "\n");
		RefGeneDestroy(pSeedRefGene);
		pSeedRefGene = NULL;
	}

	fclose(fpIn);
	fclose(fpOut);

	/* ----------------------------- */
	/* release memory                */
	/* ----------------------------- */
	for(ni=0; ni<nSourceRefNum; ni++)
	{
		RefGeneDestroy(vSourceRefGene[ni]);
		vSourceRefGene[ni] = NULL;
	}
	free(vSourceRefGene);

	for(ni=0; ni<nMapRefNum; ni++)
	{
		RefGeneDestroy(vMapRefGene[ni]);
		vMapRefGene[ni] = NULL;
	}
	free(vMapRefGene);

	if(nDestMap == 1)
	{
		for(ni=0; ni<nDestRefNum; ni++)
		{
			RefGeneDestroy(vDestRefGene[ni]);
			vDestRefGene[ni] = NULL;
		}
		free(vDestRefGene);
	}


	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGene_GetMultiOrtholog_Main()                                        */
/*  Get ortholog genes from multiple pairwise ortholog mapping.            */
/* ----------------------------------------------------------------------- */ 
int RefGene_GetMultiOrtholog_Main(char strTargetPath[], char strParamPath[], 
								  char strOutPath[])
{
	/* define */
	int nSpeciesNum;
	int nOrthoGroupNum;
	int nRefIdConvert;
	struct tagString **vSpeciesName;
	struct tagString **vSpeciesMapbase;
	struct tagString **vSpeciesDatabase;
	struct tagString **vPairwisePath;

	FILE *fpIn;
	char strLine[MED_LINE_LENGTH];
	char *chp;
	int ni;
	char strTempOutPath[LINE_LENGTH];
	char strNROutPath[LINE_LENGTH];
	char strRefIdPath[LINE_LENGTH];
	int nNeedMapBack;

	/* init */
	vSpeciesName = NULL;
	vSpeciesMapbase = NULL;
	vSpeciesDatabase = NULL;
	vPairwisePath = NULL;

	/* load parameters */
	fpIn = NULL;
	fpIn = fopen(strParamPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: cannot load parameters!\n");
		exit(EXIT_FAILURE);
	}

	/* load line by line */
	ni = 0;
	while(fgets(strLine, MED_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strstr(strLine, "[Species Number]") == strLine)
		{
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);
			nSpeciesNum = atoi(chp);
			if(nSpeciesNum <= 1)
			{
				printf("Error: RefGene_GetMultiOrtholog_Main, species number < 2!\n");
				exit(EXIT_FAILURE);
			}

			vSpeciesName = (struct tagString **)calloc(nSpeciesNum, sizeof(struct tagString*));
			if(vSpeciesName == NULL)
			{
				printf("Error: RefGene_GetMultiOrtholog_Main, cannot create database\n");
				exit(EXIT_FAILURE);
			}

			vSpeciesMapbase = (struct tagString **)calloc(nSpeciesNum, sizeof(struct tagString*));
			if(vSpeciesMapbase == NULL)
			{
				printf("Error: RefGene_GetMultiOrtholog_Main, cannot create database\n");
				exit(EXIT_FAILURE);
			}

			vSpeciesDatabase = (struct tagString **)calloc(nSpeciesNum, sizeof(struct tagString*));
			if(vSpeciesDatabase == NULL)
			{
				printf("Error: RefGene_GetMultiOrtholog_Main, cannot create database\n");
				exit(EXIT_FAILURE);
			}

			vPairwisePath = (struct tagString **)calloc(nSpeciesNum, sizeof(struct tagString*));
			if(vPairwisePath == NULL)
			{
				printf("Error: RefGene_GetMultiOrtholog_Main, cannot create database\n");
				exit(EXIT_FAILURE);
			}
		}
		else if(strstr(strLine, "[Reference Species]") == strLine)
		{
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);
			StringAddTail(vSpeciesName+ni, chp);

			fgets(strLine, MED_LINE_LENGTH, fpIn);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "[RefId From RefGeneMap(0) or RefGene(1)]") != strLine)
			{
				printf("Error: RefGene_GetMultiOrtholog_Main, parameter wrong\n");
				exit(EXIT_FAILURE);
			}
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);
			nRefIdConvert = atoi(chp);

			fgets(strLine, MED_LINE_LENGTH, fpIn);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "[Reference RefGeneMap]") != strLine)
			{
				printf("Error: RefGene_GetMultiOrtholog_Main, parameter wrong\n");
				exit(EXIT_FAILURE);
			}
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);
			StringAddTail(vSpeciesMapbase+ni, chp);

			fgets(strLine, MED_LINE_LENGTH, fpIn);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "[Reference RefGene]") != strLine)
			{
				printf("Error: RefGene_GetMultiOrtholog_Main, parameter wrong\n");
				exit(EXIT_FAILURE);
			}
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);
			StringAddTail(vSpeciesDatabase+ni, chp);

			ni++;
		}
		else if(strstr(strLine, "[Map Species Name]") == strLine)
		{
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);
			StringAddTail(vSpeciesName+ni, chp);

			fgets(strLine, MED_LINE_LENGTH, fpIn);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "[Map Species RefGeneMap]") != strLine)
			{
				printf("Error: RefGene_GetMultiOrtholog_Main, parameter wrong\n");
				exit(EXIT_FAILURE);
			}
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);
			StringAddTail(vSpeciesMapbase+ni, chp);

			fgets(strLine, MED_LINE_LENGTH, fpIn);
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strstr(strLine, "[Map Species RefGene]") != strLine)
			{
				printf("Error: RefGene_GetMultiOrtholog_Main, parameter wrong\n");
				exit(EXIT_FAILURE);
			}
			chp = strchr(strLine, '=');
			chp++;
			StrTrimLeft(chp);
			StringAddTail(vSpeciesDatabase+ni, chp);

			ni++;
		}
		else
		{
			printf("Error: RefGene_GetMultiOrtholog_Main, unknown parameters\n");
			exit(EXIT_FAILURE);
		}
	}

	if(ni != nSpeciesNum)
	{
		printf("Error: RefGene_GetMultiOrtholog_Main, species number not match\n");
		exit(EXIT_FAILURE);
	}

	/* close files */
	fclose(fpIn);

	/* convert refernce map */
	if( nRefIdConvert == 1 )
	{
		sprintf(strRefIdPath, "%s_refmap.txt", strOutPath);
		RefGene_MapBetweenDatabase(strTargetPath, strRefIdPath,
			vSpeciesName[0]->m_pString, vSpeciesDatabase[0]->m_pString, vSpeciesMapbase[0]->m_pString);
	}
	else
	{
		strcpy(strRefIdPath, strTargetPath);
	}

	/* remove redundancy */
	sprintf(strNROutPath, "%s_nr.txt", strOutPath);
	RefGene_RemoveRedundancy(strRefIdPath, vSpeciesName[0]->m_pString, vSpeciesMapbase[0]->m_pString, strNROutPath);


	/* pairwise mapping */
	for(ni=0; ni<nSpeciesNum; ni++)
	{
		sprintf(strTempOutPath, "%s_%s", strOutPath, vSpeciesName[ni]->m_pString);
		StringAddTail(vPairwisePath+ni, strTempOutPath);
		RefGene_GetOrtholog(strNROutPath, 
			vSpeciesName[0]->m_pString, vSpeciesMapbase[0]->m_pString,
			vSpeciesName[ni]->m_pString, vSpeciesMapbase[ni]->m_pString,
			"NULL", strTempOutPath);
	}
	
	/* merge multiple species */
	sprintf(strTempOutPath, "%s.mstmap", strOutPath);
	nOrthoGroupNum = RefGene_MergeMultiOrtholog(nSpeciesNum, vSpeciesName,
							   vPairwisePath, strTempOutPath);

	/* map back to the refgene database of the original species */
	nNeedMapBack = 0;
	for(ni=0; ni<nSpeciesNum; ni++)
	{
		if(strcmp(vSpeciesDatabase[ni]->m_pString, "NULL") != 0)
		{
			nNeedMapBack = 1;
			break;
		}
	}
	if(nNeedMapBack == 1)
	{
		sprintf(strTempOutPath, "%s.msomap", strOutPath);
		sprintf(strNROutPath, "%s.mstmap", strOutPath);
		RefGene_MapMultiOrthologToOriginalSpecies(nSpeciesNum, vSpeciesName, vSpeciesMapbase,
			vSpeciesDatabase, nOrthoGroupNum, strNROutPath, strTempOutPath);
	}

	/* release memory */
	for(ni=0; ni<nSpeciesNum; ni++)
	{
		DeleteString(vSpeciesName[ni]);
		DeleteString(vSpeciesMapbase[ni]);
		DeleteString(vSpeciesDatabase[ni]);
		DeleteString(vPairwisePath[ni]);
	}
	free(vSpeciesName);
	free(vSpeciesMapbase);
	free(vSpeciesDatabase);
	free(vPairwisePath);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGene_MapBetweenDatabase()                                           */
/*  map a list of refid from one database to another database.             */
/* ----------------------------------------------------------------------- */ 
int RefGene_MapBetweenDatabase(char strInPath[], char strOutPath[],
			char strSpecies[], char strFromDatabase[], char strToDatabase[])
{
	/* define */

	/* for database */
	FILE *fpSourceRefGene;
	FILE *fpDestRefGene;
	struct tagRefGene *pSourceList;
	struct tagRefGene *pTargetList;

	struct tagRefGene *pDestDatabase;
	struct tagRefGene **vDestDatabase;
	int nDestRefNum;
	int ni;

	struct tagRefGene *pRefGene,*pCurrentRefGene,*pOptRefGene,*pPrev;
	double dOptScore,dScore;
	char strRefLine[LONG_LINE_LENGTH];
	int nInserted;

	int nCDSL1,nCDSL2;
	int nCDSS,nCDSE,nCDSL;
	double dCDSR1,dCDSR2;

	/* for target search */
	FILE *fpIn;
	FILE *fpOut;
	
	char strLine[LINE_LENGTH];

	/* init */
	pSourceList = NULL;
	pTargetList = NULL;
	
	/* load source refgene */
	fpSourceRefGene = NULL;
	fpSourceRefGene = fopen(strFromDatabase, "r");
	if(fpSourceRefGene == NULL)
	{
		printf("Error: RefGene_MapBetweenDatabase, cannot open source refgene file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strRefLine, LONG_LINE_LENGTH, fpSourceRefGene) != NULL)
	{
		StrTrimLeft(strRefLine);
		StrTrimRight(strRefLine);
		if(strRefLine[0] == '\0')
			continue;

		pRefGene = NULL;
		pRefGene = RefGeneCreate();
		RefGeneInit_FromGenomeLabFormat(pRefGene, strRefLine, strSpecies);
		nInserted = 0;

		fpIn = NULL;
		fpIn = fopen(strInPath, "r");
		if(fpIn == NULL)
		{
			printf("Error: RefGene_MapBetweenDatabase, cannot open target refgene file!\n");
			exit(EXIT_FAILURE);
		}

		while(fgets(strLine, LINE_LENGTH, fpIn) != NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;
			
			if(strcmp(strLine, pRefGene->strName) == 0)
			{
				RefGeneInsert(&pSourceList, pRefGene);
				nInserted = 1;
				break;
			}
		}

		fclose(fpIn);

		if(nInserted == 0)
		{
			RefGeneDestroy(pRefGene);
		}
	}

	fclose(fpSourceRefGene);


	/* load dest database */
	fpDestRefGene = NULL;
	fpDestRefGene = fopen(strToDatabase, "r");
	if(fpDestRefGene == NULL)
	{
		printf("Error: RefGene_MapBetweenDatabase, cannot open dest refgene file!\n");
		exit(EXIT_FAILURE);
	}

	pDestDatabase = NULL;
	nDestRefNum = 0;
	pCurrentRefGene = NULL;
	while(fgets(strRefLine, LONG_LINE_LENGTH, fpDestRefGene) != NULL)
	{
		StrTrimLeft(strRefLine);
		StrTrimRight(strRefLine);
		if(strRefLine[0] == '\0')
			continue;

		pRefGene = NULL;
		pRefGene = RefGeneCreate();
		RefGeneInit_FromGenomeLabFormat(pRefGene, strRefLine, strSpecies);
		if(pDestDatabase == NULL)
		{
			pDestDatabase = pRefGene;
			pCurrentRefGene = pRefGene;
		}
		else
		{
			pCurrentRefGene->pNext = pRefGene;
			pCurrentRefGene = pRefGene;
		}
		nDestRefNum++;
	}

	fclose(fpDestRefGene);

	vDestDatabase = (struct tagRefGene **)calloc(nDestRefNum, sizeof(struct tagRefGene*));
	if(vDestDatabase == NULL)
	{
		printf("Error: RefGene_MapBetweenDatabase, organize refgene into array!\n");
		exit(EXIT_FAILURE);
	}

	ni=0;
	while(pDestDatabase != NULL)
	{
		pRefGene = pDestDatabase;
		pDestDatabase = pRefGene->pNext;
		pRefGene->pNext = NULL;
		vDestDatabase[ni] = pRefGene;
		ni++;
	}

	if(ni != nDestRefNum)
	{
		printf("Error: RefGene_MapBetweenDatabase, refgene number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* TODO: map refgene */
	pPrev = NULL;
	pRefGene = pSourceList;
	while(pRefGene != NULL)
	{
		dOptScore = 0.0;
		pOptRefGene = NULL;
		for(ni=0; ni<nDestRefNum; ni++)
		{
			if(vDestDatabase[ni]->nChrom < pRefGene->nChrom)
			{
				continue;
			}
			else if(vDestDatabase[ni]->nChrom == pRefGene->nChrom)
			{
				if(vDestDatabase[ni]->nTxEnd < pRefGene->nTxStart)
				{
					continue;
				}
				else if(pRefGene->nTxEnd < vDestDatabase[ni]->nTxStart)
				{
					break;
				}
				else
				{
					dScore = RefGene_Match(vDestDatabase[ni], pRefGene);
					
					/* NOTE: below is for one-multi mapping */
					/*
					if(dScore > REFGENE_MATCH_TH)
					{						 
						if(pTargetList == NULL)
						{
							pTargetList = vDestDatabase[ni];
							pPrev = vDestDatabase[ni];
						}
						else
						{
							pPrev->pNext = vDestDatabase[ni];
							pPrev = vDestDatabase[ni];
						}
					}
					else
					{
						nCDSS = pRefGene->nCdsStart;
						if(vDestDatabase[ni]->nCdsStart > nCDSS)
							nCDSS = vDestDatabase[ni]->nCdsStart;
						nCDSE = pRefGene->nCdsEnd;
						if(vDestDatabase[ni]->nCdsEnd < nCDSE)
							nCDSE = vDestDatabase[ni]->nCdsEnd;

						nCDSL1 = pRefGene->nCdsEnd-pRefGene->nCdsStart+1;
						nCDSL2 = vDestDatabase[ni]->nCdsEnd-vDestDatabase[ni]->nCdsStart+1;
						nCDSL = nCDSE-nCDSS+1;

						dCDSR1 = (double)nCDSL/(double)nCDSL1;
						dCDSR2 = (double)nCDSL/(double)nCDSL2;
						if(dCDSR2<dCDSR1)
							dCDSR1 = dCDSR2;

						if(dCDSR1 > REFGENE_MATCH_TH)
						{
							if(pTargetList == NULL)
							{
								pTargetList = vDestDatabase[ni];
								pPrev = vDestDatabase[ni];
							}
							else
							{
								pPrev->pNext = vDestDatabase[ni];
								pPrev = vDestDatabase[ni];
							}
						}
					}
					*/
					/* END NOTE */

					/* NOTE: below is for one-one mapping */
					if(dScore > dOptScore)
					{
						if(dScore > REFGENE_MATCH_TH)
						{
							dOptScore = dScore;
							pOptRefGene = vDestDatabase[ni];
						}
						else
						{
							nCDSS = pRefGene->nCdsStart;
							if(vDestDatabase[ni]->nCdsStart > nCDSS)
								nCDSS = vDestDatabase[ni]->nCdsStart;
							nCDSE = pRefGene->nCdsEnd;
							if(vDestDatabase[ni]->nCdsEnd < nCDSE)
								nCDSE = vDestDatabase[ni]->nCdsEnd;

							nCDSL1 = pRefGene->nCdsEnd-pRefGene->nCdsStart+1;
							nCDSL2 = vDestDatabase[ni]->nCdsEnd-vDestDatabase[ni]->nCdsStart+1;
							nCDSL = nCDSE-nCDSS+1;

							dCDSR1 = (double)nCDSL/(double)nCDSL1;
							dCDSR2 = (double)nCDSL/(double)nCDSL2;
							if(dCDSR2<dCDSR1)
								dCDSR1 = dCDSR2;

							if(dCDSR1 > REFGENE_MATCH_TH)
							{
								dOptScore = dScore;
								pOptRefGene = vDestDatabase[ni];
							}
						}
					}
					/* END NOTE */
				}
			}
			else
			{
				break;
			}
		}

		/* NOTE: below is for one-one mapping */
		if(pOptRefGene != NULL)
		{
			if(pTargetList == NULL)
			{
				pTargetList = pOptRefGene;
				pPrev = pOptRefGene;
			}
			else
			{
				pPrev->pNext = pOptRefGene;
				pPrev = pOptRefGene;
			}
		}
		/* END NOTE */

		pRefGene = pRefGene->pNext;
	}

	/* export */
	fpOut = NULL;
	fpOut = fopen(strOutPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: RefGene_MapBetweenDatabase, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	pRefGene = pTargetList;
	while(pRefGene != NULL)
	{
		fprintf(fpOut, "%s\n", pRefGene->strName);
		pRefGene = pRefGene->pNext;
	}
	
	fclose(fpOut);

	/* release memory */
	while(pTargetList != NULL)
	{
		pRefGene = pTargetList;
		pTargetList = pRefGene->pNext;
		pRefGene->pNext = NULL;
	}

	for(ni=0; ni<nDestRefNum; ni++)
	{
		RefGeneDestroy(vDestDatabase[ni]);
		vDestDatabase[ni] = NULL;
	}
	free(vDestDatabase);

	
	while(pSourceList != NULL)
	{
		pRefGene = pSourceList;
		pSourceList = pRefGene->pNext;
		pRefGene->pNext = NULL;
		RefGeneDestroy(pRefGene);
	}

	/* return */
	return PROC_SUCCESS;
}


/* ----------------------------------------------------------------------- */ 
/*  RefGene_RemoveRedundancy()                                             */
/*  Remove redundancy in the target list.                                  */
/* ----------------------------------------------------------------------- */ 
int RefGene_RemoveRedundancy(char strTargetPath[], char strSpecies[], 
							 char strRefGenePath[], char strOutPath[])
{
	/* define */
	
	/* for database */
	FILE *fpSourceRefGene;
	struct tagRefGene *pSourceList;
	struct tagRefGene *pTargetList;
	
	struct tagRefGene *pRefGene,*pPrev,*pCurr;
	char strRefLine[LONG_LINE_LENGTH];
	int nRed;
	int nInserted;

	/* for target search */
	FILE *fpIn;
	FILE *fpOut;
	
	char strLine[LINE_LENGTH];

	/* init */
	pSourceList = NULL;
	pTargetList = NULL;
	
	/* load source refgene */
	fpSourceRefGene = NULL;
	fpSourceRefGene = fopen(strRefGenePath, "r");
	if(fpSourceRefGene == NULL)
	{
		printf("Error: RemoveRedundancy, cannot open source refgene file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strRefLine, LONG_LINE_LENGTH, fpSourceRefGene) != NULL)
	{
		StrTrimLeft(strRefLine);
		StrTrimRight(strRefLine);
		if(strRefLine[0] == '\0')
			continue;

		pRefGene = NULL;
		pRefGene = RefGeneCreate();
		RefGeneInit_FromGenomeLabFormat(pRefGene, strRefLine, strSpecies);
		nInserted = 0;

		fpIn = NULL;
		fpIn = fopen(strTargetPath, "r");
		if(fpIn == NULL)
		{
			printf("Error: RemoveRedundancy, cannot open target refgene file!\n");
			exit(EXIT_FAILURE);
		}

		while(fgets(strLine, LINE_LENGTH, fpIn) != NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;
			
			if(strcmp(strLine, pRefGene->strName) == 0)
			{
				RefGeneInsert(&pSourceList, pRefGene);
				nInserted = 1;
				break;
			}

		}

		fclose(fpIn);

		if(nInserted == 0)
		{
			RefGeneDestroy(pRefGene);
		}
	}

	fclose(fpSourceRefGene);


	/* remove overlap */
	pPrev = NULL;
	pCurr = NULL;
	while(pSourceList != NULL)
	{
		pRefGene = pSourceList;
		pSourceList = pRefGene->pNext;
		pRefGene->pNext = NULL;

		if(pTargetList == NULL)
		{
			pTargetList = pRefGene;
			pPrev = NULL;
			pCurr = pRefGene;
		}
		else
		{
			/* if overlap */
			if(RefGene_Overlap(pCurr, pRefGene) == 1)
			{
				if( (pCurr->nTxEnd-pCurr->nTxStart) >= (pRefGene->nTxEnd-pRefGene->nTxStart) )
				{
					RefGeneDestroy(pRefGene);
				}
				else
				{
					if(pPrev == NULL)
					{
						pTargetList = pRefGene;
						RefGeneDestroy(pCurr);
						pCurr = pRefGene;
					}
					else
					{
						pPrev->pNext = pRefGene;
						RefGeneDestroy(pCurr);
						pCurr = pRefGene;
					}
				}

			}
			/* if not overlap */
			else
			{
				pCurr->pNext = pRefGene;
				pPrev = pCurr;
				pCurr = pRefGene;
			}
		}
	}
	
	/* remove common name & write */
	fpOut = NULL;
	fpOut = fopen(strOutPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: RefGene_RemoveRedundancy, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	pRefGene = pTargetList;
	while(pRefGene != NULL)
	{
		nRed = 0;
		pPrev = pTargetList;
		while(pPrev != pRefGene)
		{
			if(strcmp(pPrev->strName, pRefGene->strName) == 0)
			{
				nRed = 1;
				break;
			}
			/* get next */
			pPrev = pPrev->pNext;
		}

		if(nRed == 0)
		{
			fprintf(fpOut, "%s\n", pRefGene->strName);
		}

		/* get next */
		pRefGene = pRefGene->pNext;
	}

	fclose(fpOut);

	/* release memory */
	while(pTargetList != NULL)
	{
		pRefGene = pTargetList;
		pTargetList = pRefGene->pNext;
		pRefGene->pNext = NULL;
		RefGeneDestroy(pRefGene);
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefFlat_RemoveRedundancy()                                             */
/*  Remove redundancy in the target list.                                  */
/* ----------------------------------------------------------------------- */ 
int RefFlat_RemoveRedundancy(char strTargetPath[], char strSpecies[], 
							 char strRefGenePath[], char strOutPath[])
{
	/* define */
	
	/* for database */
	FILE *fpSourceRefGene;
	struct tagRefGene *pSourceList;
	struct tagRefGene *pTargetList;
	
	struct tagRefGene *pRefGene,*pPrev,*pCurr;
	char strRefLine[LONG_LINE_LENGTH];
	int nRed;
	int nInserted;

	/* for target search */
	FILE *fpIn;
	FILE *fpOut;
	
	char strLine[LINE_LENGTH];

	/* init */
	pSourceList = NULL;
	pTargetList = NULL;
	
	/* load source refgene */
	fpSourceRefGene = NULL;
	fpSourceRefGene = fopen(strRefGenePath, "r");
	if(fpSourceRefGene == NULL)
	{
		printf("Error: RefFlat_RemoveRedundancy, cannot open source refgene file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strRefLine, LONG_LINE_LENGTH, fpSourceRefGene) != NULL)
	{
		StrTrimLeft(strRefLine);
		StrTrimRight(strRefLine);
		if(strRefLine[0] == '\0')
			continue;

		pRefGene = NULL;
		pRefGene = RefGeneCreate();
		RefFlatInit_FromGenomeLabFormat(pRefGene, strRefLine, strSpecies);
		nInserted = 0;

		fpIn = NULL;
		fpIn = fopen(strTargetPath, "r");
		if(fpIn == NULL)
		{
			printf("Error: RefGene_RemoveRedundancy, cannot open target refgene file!\n");
			exit(EXIT_FAILURE);
		}

		while(fgets(strLine, LINE_LENGTH, fpIn) != NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;
			
			if(strcmp(strLine, pRefGene->strName) == 0)
			{
				RefGeneInsert(&pSourceList, pRefGene);
				nInserted = 1;
				break;
			}

		}

		fclose(fpIn);

		if(nInserted == 0)
		{
			RefGeneDestroy(pRefGene);
		}
	}

	fclose(fpSourceRefGene);


	/* remove overlap */
	pPrev = NULL;
	pCurr = NULL;
	while(pSourceList != NULL)
	{
		pRefGene = pSourceList;
		pSourceList = pRefGene->pNext;
		pRefGene->pNext = NULL;

		if(pTargetList == NULL)
		{
			pTargetList = pRefGene;
			pPrev = NULL;
			pCurr = pRefGene;
		}
		else
		{
			/* if overlap */
			if(RefGene_Overlap(pCurr, pRefGene) == 1)
			{
				if( (pCurr->nTxEnd-pCurr->nTxStart) >= (pRefGene->nTxEnd-pRefGene->nTxStart) )
				{
					RefGeneDestroy(pRefGene);
				}
				else
				{
					if(pPrev == NULL)
					{
						pTargetList = pRefGene;
						RefGeneDestroy(pCurr);
						pCurr = pRefGene;
					}
					else
					{
						pPrev->pNext = pRefGene;
						RefGeneDestroy(pCurr);
						pCurr = pRefGene;
					}
				}

			}
			/* if not overlap */
			else
			{
				pCurr->pNext = pRefGene;
				pPrev = pCurr;
				pCurr = pRefGene;
			}
		}
	}
	
	/* remove common name & write */
	fpOut = NULL;
	fpOut = fopen(strOutPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: RefGene_RemoveRedundancy, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	pRefGene = pTargetList;
	while(pRefGene != NULL)
	{
		nRed = 0;
		pPrev = pTargetList;
		while(pPrev != pRefGene)
		{
			if(strcmp(pPrev->strName, pRefGene->strName) == 0)
			{
				nRed = 1;
				break;
			}
			/* get next */
			pPrev = pPrev->pNext;
		}

		if(nRed == 0)
		{
			fprintf(fpOut, "%s\n", pRefGene->strName);
		}

		/* get next */
		pRefGene = pRefGene->pNext;
	}

	fclose(fpOut);

	/* release memory */
	while(pTargetList != NULL)
	{
		pRefGene = pTargetList;
		pTargetList = pRefGene->pNext;
		pRefGene->pNext = NULL;
		RefGeneDestroy(pRefGene);
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGene_Overlap()                                                      */
/*  Judge if two refgene entry are redundant. If yes, return 1; else 0.    */
/* ----------------------------------------------------------------------- */ 
int RefGene_Overlap(struct tagRefGene *pRefGene1, struct tagRefGene *pRefGene2)
{
	/* define */
	int nOverlap;
	int nRegionOverlap;
	int ni,nj;
	int nE1,nE2,nD1,nD2;
	int nL1,nL2,nML;
	int nUpdate1,nUpdate2;
	int nMinEnd,nMaxStart;
	double dR1,dR2;

	/* init */
	nOverlap = 0;

	/* judge */
	if((pRefGene1 == NULL) || (pRefGene2 == NULL))
	{
		return nOverlap;
	}

	/* same chromosome ? */
	if(pRefGene1->nChrom != pRefGene2->nChrom)
	{
		return nOverlap;
	}
	/* same strand ? */
	if(pRefGene1->chStrand != pRefGene2->chStrand)
	{
		return nOverlap;
	}

	/* region overlap ? */
	nRegionOverlap = 0;
	if( (pRefGene1->nTxStart >= pRefGene2->nTxStart) && (pRefGene1->nTxStart <= pRefGene2->nTxEnd) )
	{
		nRegionOverlap = 1;
	}
	else if((pRefGene2->nTxStart >= pRefGene1->nTxStart) && (pRefGene2->nTxStart <= pRefGene1->nTxEnd))
	{
		nRegionOverlap = 1;
	}

	if(nRegionOverlap == 0)
	{
		return nOverlap;
	}

	/* structure overlap ? */
	ni = 0;
	nj = 0;
	nL1 = 0;
	nL2 = 0;
	nML = 0;
	nUpdate1 = 1;
	nUpdate2 = 1;
	
	while((ni<pRefGene1->nExonCount) && (nj<pRefGene2->nExonCount))
	{
		if(nUpdate1 == 1)
		{
			nD1 = IMGETAT(pRefGene1->pmatExonStartsEnds, ni, 0);
			nD2 = IMGETAT(pRefGene1->pmatExonStartsEnds, ni, 1);
			nL1 += nD2-nD1+1;
		}
		if(nUpdate2 == 1)
		{
			nE1 = IMGETAT(pRefGene2->pmatExonStartsEnds, nj, 0);
			nE2 = IMGETAT(pRefGene2->pmatExonStartsEnds, nj, 1);
			nL2 += nE2-nE1+1;
		}

		if(nD1 > nE1)
		{
			nMaxStart = nD1;
		}
		else
		{
			nMaxStart = nE1;
		}

		if(nD2 < nE2)
		{
			nMinEnd = nD2;
		}
		else
		{
			nMinEnd = nE2;
		}

		if(nMinEnd >= nMaxStart)
		{
			nML += nMinEnd-nMaxStart+1;
		}

		if(nD2 > nE2)
		{
			nj++;
			nUpdate1 = 0;
			nUpdate2 = 1;
		}
		else if(nD2 < nE2)
		{
			ni++;
			nUpdate1 = 1;
			nUpdate2 = 0;
		}
		else
		{
			ni++;
			nj++;
			nUpdate1 = 1;
			nUpdate2 = 1;
		}
	}

	if( (ni == pRefGene1->nExonCount) && (nj < pRefGene2->nExonCount) )
	{
		if(nUpdate2 == 0)
			nj++;

		while(nj<pRefGene2->nExonCount)
		{
			nE1 = IMGETAT(pRefGene2->pmatExonStartsEnds, nj, 0);
			nE2 = IMGETAT(pRefGene2->pmatExonStartsEnds, nj, 1);
			nL2 += nE2-nE1+1;
			nj++;
		}
	}
	else if( (ni < pRefGene1->nExonCount) && (nj == pRefGene2->nExonCount) )
	{
		if(nUpdate1 == 0)
			ni++;

		while(ni<pRefGene1->nExonCount)
		{
			nD1 = IMGETAT(pRefGene1->pmatExonStartsEnds, ni, 0);
			nD2 = IMGETAT(pRefGene1->pmatExonStartsEnds, ni, 1);
			nL1 += nD2-nD1+1;
			ni++;
		}
	}

	dR1 = (double)nML/(double)nL1;
	dR2 = (double)nML/(double)nL2;

	if( ( dR1 >= REFGENE_OVERLAP_TH ) && ( dR2 >= REFGENE_OVERLAP_TH) )
	{
		nOverlap = 1;
	}

	/* return */
	return nOverlap;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGene_Match()                                                        */
/*  Evaluate how two refgene entries are similar to each other.            */
/* ----------------------------------------------------------------------- */ 
double RefGene_Match(struct tagRefGene *pRefGene1, struct tagRefGene *pRefGene2)
{
	/* define */
	double dScore;
	int nRegionOverlap;
	int ni,nj;
	int nE1,nE2,nD1,nD2;
	int nL1,nL2,nML;
	int nUpdate1,nUpdate2;
	int nMinEnd,nMaxStart;
	double dR1,dR2;

	/* init */
	dScore = 0.0;

	/* judge */
	if((pRefGene1 == NULL) || (pRefGene2 == NULL))
	{
		return dScore;
	}

	/* same chromosome ? */
	if(pRefGene1->nChrom != pRefGene2->nChrom)
	{
		return dScore;
	}
	/* same strand ? */
	if(pRefGene1->chStrand != pRefGene2->chStrand)
	{
		return dScore;
	}

	/* region overlap ? */
	nRegionOverlap = 0;
	if( (pRefGene1->nTxStart >= pRefGene2->nTxStart) && (pRefGene1->nTxStart <= pRefGene2->nTxEnd) )
	{
		nRegionOverlap = 1;
	}
	else if((pRefGene2->nTxStart >= pRefGene1->nTxStart) && (pRefGene2->nTxStart <= pRefGene1->nTxEnd))
	{
		nRegionOverlap = 1;
	}

	if(nRegionOverlap == 0)
	{
		return dScore;
	}

	/* structure overlap ? */
	ni = 0;
	nj = 0;
	nL1 = 0;
	nL2 = 0;
	nML = 0;
	nUpdate1 = 1;
	nUpdate2 = 1;
	
	while((ni<pRefGene1->nExonCount) && (nj<pRefGene2->nExonCount))
	{
		if(nUpdate1 == 1)
		{
			nD1 = IMGETAT(pRefGene1->pmatExonStartsEnds, ni, 0);
			nD2 = IMGETAT(pRefGene1->pmatExonStartsEnds, ni, 1);
			nL1 += nD2-nD1+1;
		}
		if(nUpdate2 == 1)
		{
			nE1 = IMGETAT(pRefGene2->pmatExonStartsEnds, nj, 0);
			nE2 = IMGETAT(pRefGene2->pmatExonStartsEnds, nj, 1);
			nL2 += nE2-nE1+1;
		}

		if(nD1 > nE1)
		{
			nMaxStart = nD1;
		}
		else
		{
			nMaxStart = nE1;
		}

		if(nD2 < nE2)
		{
			nMinEnd = nD2;
		}
		else
		{
			nMinEnd = nE2;
		}

		if(nMinEnd >= nMaxStart)
		{
			nML += nMinEnd-nMaxStart+1;
		}

		if(nD2 > nE2)
		{
			nj++;
			nUpdate1 = 0;
			nUpdate2 = 1;
		}
		else if(nD2 < nE2)
		{
			ni++;
			nUpdate1 = 1;
			nUpdate2 = 0;
		}
		else
		{
			ni++;
			nj++;
			nUpdate1 = 1;
			nUpdate2 = 1;
		}
	}

	if( (ni == pRefGene1->nExonCount) && (nj < pRefGene2->nExonCount) )
	{
		if(nUpdate2 == 0)
			nj++;

		while(nj<pRefGene2->nExonCount)
		{
			nE1 = IMGETAT(pRefGene2->pmatExonStartsEnds, nj, 0);
			nE2 = IMGETAT(pRefGene2->pmatExonStartsEnds, nj, 1);
			nL2 += nE2-nE1+1;
			nj++;
		}
	}
	else if( (ni < pRefGene1->nExonCount) && (nj == pRefGene2->nExonCount) )
	{
		if(nUpdate1 == 0)
			ni++;

		while(ni<pRefGene1->nExonCount)
		{
			nD1 = IMGETAT(pRefGene1->pmatExonStartsEnds, ni, 0);
			nD2 = IMGETAT(pRefGene1->pmatExonStartsEnds, ni, 1);
			nL1 += nD2-nD1+1;
			ni++;
		}
	}

	dR1 = (double)nML/(double)nL1;
	dR2 = (double)nML/(double)nL2;

	if(dR1 < dR2)
		dScore = dR1;
	else
		dScore = dR2;

	/* return */
	return dScore;
}


/* ----------------------------------------------------------------------- */ 
/*  RefGene_MapMultiOrthologToOriginalSpecies()                            */
/*  Map refgene orthologs back to the original species.                    */
/* ----------------------------------------------------------------------- */ 
int RefGene_MapMultiOrthologToOriginalSpecies(int nSpeciesNum, struct tagString **vSpeciesName, 
							struct tagString **vSpeciesMapbase, struct tagString **vSpeciesDatabase, 
							int nOrthoGroupNum, char strInPath[], char strOutPath[])
{
	/* define */
	struct DOUBLEMATRIX *pCountMat;
	struct DOUBLEMATRIX *pOrienMat;
	struct tagRefGene **vRefGene;
	struct tagString **vRefName;
	int nTotRefNum;
	int ni,nj,nk,nx;
	char *chp,*chp2;
	char strRefLine[LONG_LINE_LENGTH];
	char strRefLine0[LONG_LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	int nMatchCount,nOrienCount;
	struct tagRefGene *pRefGene,*pCurrentRefGene,*pBestRefGene;
	double dBestScore,dScore;
	int nCDSL1,nCDSL2;
	int nCDSS,nCDSE,nCDSL;
	double dCDSR1,dCDSR2;

	FILE *fpIn;
	FILE *fpOut;

	/* for creating maps */
	FILE *fpDestRefGene;
	struct tagRefGene *pDestRefGeneList;
	int nDestRefNum;

	/* init */
	if(nSpeciesNum <= 0)
	{
		printf("Warning: no species available!\n");
		return PROC_FAILURE;
	}
	if(nOrthoGroupNum <= 0)
	{
		printf("Warning: no genes available!\n");
		return PROC_FAILURE;
	}

	pCountMat = NULL;
	pCountMat = CreateDoubleMatrix(nOrthoGroupNum, nSpeciesNum);
	if(pCountMat == NULL)
	{
		printf("Error: RefGene_MapMultiOrthologToOriginalSpecies, cannot create matrix!\n");
		exit(EXIT_FAILURE);
	}
	pOrienMat = NULL;
	pOrienMat = CreateDoubleMatrix(nOrthoGroupNum, nSpeciesNum);
	if(pOrienMat == NULL)
	{
		printf("Error: RefGene_MapMultiOrthologToOriginalSpecies, cannot create matrix!\n");
		exit(EXIT_FAILURE);
	}
	nTotRefNum = nSpeciesNum*nOrthoGroupNum;
	vRefGene = NULL;
	vRefGene = (struct tagRefGene **)calloc(nTotRefNum, sizeof(struct tagRefGene *));
	if(vRefGene == NULL)
	{
		printf("Error: RefGene_MapMultiOrthologToOriginalSpecies, cannot create matrix!\n");
		exit(EXIT_FAILURE);
	}
	vRefName = NULL;
	vRefName = (struct tagString **)calloc(nOrthoGroupNum, sizeof(struct tagString *));
	if(vRefName == NULL)
	{
		printf("Error: RefGene_MapMultiOrthologToOriginalSpecies, cannot create matrix!\n");
		exit(EXIT_FAILURE);
	}

	/* load orthologs */
	fpIn = NULL;
	fpIn = fopen(strInPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: RefGene_MapMultiOrthologToOriginalSpecies, cannot load refgene!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	nk = 0;
	while(fgets(strRefLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strRefLine);
		StrTrimRight(strRefLine);
		if(strRefLine[0] == '\0')
			continue;

		if(strRefLine[0] == '>')
		{
			chp = strRefLine+1;
			StringAddTail(vRefName+ni, chp);

			for(nj=0; nj<nSpeciesNum; nj++)
			{
				fgets(strRefLine, LONG_LINE_LENGTH, fpIn);
				StrTrimLeft(strRefLine);
				StrTrimRight(strRefLine);
				chp = strRefLine;
				chp2 = strchr(chp, '\t');
				*chp2 = '\0';
				strcpy(strSpecies, chp);
				chp = chp2+1;
				chp2 = strchr(chp, '\t');
				*chp2 = '\0';
				nMatchCount = atoi(chp);
				chp = chp2+1;
				chp2 = strchr(chp, '\t');
				*chp2 = '\0';
				nOrienCount = atoi(chp);

				DMSETAT(pCountMat, ni, nj, (double)nMatchCount);
				DMSETAT(pOrienMat, ni, nj, (double)nOrienCount);

				chp = chp2+1;
				strcpy(strRefLine0, chp);
				
				if(strcmp(strSpecies, vSpeciesName[nj]->m_pString) != 0)
				{
					printf("Error: RefGene_MapMultiOrthologToOriginalSpecies, species not match!\n");
					exit(EXIT_FAILURE);
				}

				if(strstr(strRefLine0, "---") == strRefLine0)
				{
					nk++;
				}
				else
				{
					pRefGene = NULL;
					pRefGene = RefGeneCreate();
					RefGeneInit_FromGenomeLabFormat(pRefGene, strRefLine0, strSpecies);
					vRefGene[nk] = pRefGene;
					nk++;
				}
			}
			ni++;
		}
	}

	fclose(fpIn);

	if((ni != nOrthoGroupNum) || (nk != nTotRefNum))
	{
		printf("Error: RefGene_MapMultiOrthologToOriginalSpecies, loading error!\n");
		exit(EXIT_FAILURE);
	}

	/* map back */
	for(nj=0; nj<nSpeciesNum; nj++)
	{
		/* if no need to map, then keep the original one */
		if(strcmp(vSpeciesDatabase[nj]->m_pString, "NULL") == 0)
			continue;

		/* if need map, create database first */
		strcpy(strSpecies, vSpeciesName[nj]->m_pString);
		fpDestRefGene = NULL;
		fpDestRefGene = fopen(vSpeciesDatabase[nj]->m_pString, "r");
		if(fpDestRefGene == NULL)
		{
			printf("Error: RefGene_MapMultiOrthologToOriginalSpecies, cannot open dest refgene file!\n");
			exit(EXIT_FAILURE);
		}

		nDestRefNum = 0;
		pDestRefGeneList = NULL;
		while(fgets(strRefLine, LONG_LINE_LENGTH, fpDestRefGene) != NULL)
		{
			StrTrimLeft(strRefLine);
			StrTrimRight(strRefLine);
			if(strRefLine[0] == '\0')
				continue;

			pRefGene = NULL;
			pRefGene = RefGeneCreate();
			RefGeneInit_FromGenomeLabFormat(pRefGene, strRefLine, strSpecies);
			if(pDestRefGeneList == NULL)
			{
				pDestRefGeneList = pRefGene;
				pCurrentRefGene = pRefGene;
			}
			else
			{
				pCurrentRefGene->pNext = pRefGene;
				pCurrentRefGene = pRefGene;
			}
			nDestRefNum++;
		}

		fclose(fpDestRefGene);

		/* remap */
		for(ni=0; ni<nOrthoGroupNum; ni++)
		{
			nk = ni*nSpeciesNum+nj;
			if(vRefGene[nk] != NULL)
			{
				pBestRefGene = NULL;
				dBestScore = 0.0;

				pRefGene = pDestRefGeneList;
				nx = 0;
				while(pRefGene != NULL)
				{				
					dScore = RefGene_Match(vRefGene[nk], pRefGene);
					if(dScore > dBestScore)
					{
						if(dScore > REFGENE_MATCH_TH)
						{
							dBestScore = dScore;
							pBestRefGene = pRefGene;
						}
						else
						{
							nCDSS = pRefGene->nCdsStart;
							if(vRefGene[nk]->nCdsStart > nCDSS)
								nCDSS = vRefGene[nk]->nCdsStart;
							nCDSE = pRefGene->nCdsEnd;
							if(vRefGene[nk]->nCdsEnd < nCDSE)
								nCDSE = vRefGene[nk]->nCdsEnd;

							nCDSL1 = pRefGene->nCdsEnd-pRefGene->nCdsStart+1;
							nCDSL2 = vRefGene[nk]->nCdsEnd-vRefGene[nk]->nCdsStart+1;
							nCDSL = nCDSE-nCDSS+1;

							dCDSR1 = (double)nCDSL/(double)nCDSL1;
							dCDSR2 = (double)nCDSL/(double)nCDSL2;
							if(dCDSR2<dCDSR1)
								dCDSR1 = dCDSR2;

							if(dCDSR1 > REFGENE_MATCH_TH)
							{
								dBestScore = dScore;
								pBestRefGene = pRefGene;
							}
						}
					}

					pRefGene = pRefGene->pNext;

					/* nx++;
					if(nx == nDestRefNum)
					{
						nx = nx;
					} */
				}

				if(pBestRefGene != NULL)
				{
					RefGeneDestroy(vRefGene[nk]);
					vRefGene[nk] = RefGeneClone(pBestRefGene);
				}
			}
		}

		/* clear data base */
		while(pDestRefGeneList != NULL)
		{
			pRefGene = pDestRefGeneList;
			pDestRefGeneList = pRefGene->pNext;
			pRefGene->pNext = NULL;
			RefGeneDestroy(pRefGene);
		}

	}

	/* write */
	fpOut = NULL;
	fpOut = fopen(strOutPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: RefGene_MapMultiOrthologToOriginalSpecies, cannot write refgene!\n");
		exit(EXIT_FAILURE);
	}

	nk = 0;
	for(ni=0; ni<nOrthoGroupNum; ni++)
	{
		fprintf(fpOut, ">%s\n", vRefName[ni]->m_pString);
		for(nj=0; nj<nSpeciesNum; nj++)
		{
			nMatchCount = (int)DMGETAT(pCountMat, ni, nj);
			nOrienCount = (int)DMGETAT(pOrienMat, ni, nj);
			fprintf(fpOut, "%s\t%d\t%d\t", vSpeciesName[nj]->m_pString,
				nMatchCount, nOrienCount);
			if(vRefGene[nk] != NULL)
			{
				RefGeneWrite(vRefGene[nk], fpOut);
			}
			else
			{
				fprintf(fpOut, "---\n");
			}
			nk++;
		}
		fprintf(fpOut, "\n");
	}

	fclose(fpOut);

	/* release memory */
	DestroyDoubleMatrix(pCountMat);
	DestroyDoubleMatrix(pOrienMat);
	for(ni=0; ni<nTotRefNum; ni++)
	{
		if(vRefGene[ni] != NULL)
		{
			RefGeneDestroy(vRefGene[ni]);
			vRefGene[ni] = NULL;
		}
	}
	free(vRefGene);
	for(ni=0; ni<nOrthoGroupNum; ni++)
	{
		DeleteString(vRefName[ni]);
		vRefName[ni] = NULL;
	}
	free(vRefName);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGene_MultiOrthologGetTSSAround()                                    */
/*  Get TSS UP and Down sequence for multiortholog.                        */
/* ----------------------------------------------------------------------- */ 
int RefGene_MultiOrthologGetTSSAround(char strInPath[], int nTSSUP, int nTSSDOWN,
									  char strParamPath[], char strOutPath[],
									  char strSeqFile[])
{
	/* define */
	
	/* environment */
	int nSpeciesNum;
	struct tagString **vSpeciesName;
	struct tagString **vSpeciesSeqbase;
	struct tagString **vSpeciesConsbase;
	struct tagString **vSpeciesAnnotbase;
	struct INTMATRIX **vChrLen;
	struct tagString **vSpeciesCord;
	struct tagString **vSpeciesOut;
	int ni,nj;
	int nPos1,nPos2;
	int nChrid,nChrLen;

	/* ortholog */
	int nMatchCount,nOrienCount;
	struct tagRefGene *pRefGene;
	char strSpecies[LINE_LENGTH];

	char strRefName[LINE_LENGTH];
	char strStrandType[LINE_LENGTH];
	int nStrandType;
	char strConservationType[LINE_LENGTH];

	char *chp,*chp2;
	char strRefLine[LONG_LINE_LENGTH];
	char strRefLine0[LONG_LINE_LENGTH];
	
	FILE *fpIn;
	FILE **vfpOut;

	/* init environment */
	fpIn = NULL;
	fpIn = fopen(strParamPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: RefGene_MultiOrthologGetTSSAround, cannot load environment settings!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strRefLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strRefLine);
		StrTrimRight(strRefLine);
		if(strRefLine[0] == '\0')
			continue;

		if(strstr(strRefLine, "[Strand Type]") == strRefLine)
		{
			chp = strchr(strRefLine, '=');
			chp++;
			StrTrimLeft(chp);
			sprintf(strStrandType, "%s", chp);
			StrMakeUpper(strStrandType);
			if(strcmp(strStrandType, "GENEWISE") == 0)
			{
				nStrandType = 1;
			}
			else
			{
				nStrandType = 0;
			}
		}
		else if(strstr(strRefLine, "[Conservation Type]") == strRefLine)
		{
			chp = strchr(strRefLine, '=');
			chp++;
			StrTrimLeft(chp);
			sprintf(strConservationType, "%s", chp);
		}
		else if(strstr(strRefLine, "[Species Number]") == strRefLine)
		{
			chp = strchr(strRefLine, '=');
			chp++;
			StrTrimLeft(chp);
			nSpeciesNum = atoi(chp);

			if(nSpeciesNum <= 0)
			{
				printf("Warning: RefGene_MultiOrthologGetTSSAround, no species!\n");
				return PROC_SUCCESS;
			}

			vSpeciesName = NULL;
			vSpeciesName = (struct tagString **)calloc(nSpeciesNum, sizeof(struct tagString *));
			if(vSpeciesName == NULL)
			{
				printf("Error: RefGene_MultiOrthologGetTSSAround, cannot load environment settings!\n");
				exit(EXIT_FAILURE);
			}

			vSpeciesSeqbase = NULL;
			vSpeciesSeqbase = (struct tagString **)calloc(nSpeciesNum, sizeof(struct tagString *));
			if(vSpeciesSeqbase == NULL)
			{
				printf("Error: RefGene_MultiOrthologGetTSSAround, cannot load environment settings!\n");
				exit(EXIT_FAILURE);
			}

			vSpeciesConsbase = NULL;
			vSpeciesConsbase = (struct tagString **)calloc(nSpeciesNum, sizeof(struct tagString *));
			if(vSpeciesConsbase == NULL)
			{
				printf("Error: RefGene_MultiOrthologGetTSSAround, cannot load environment settings!\n");
				exit(EXIT_FAILURE);
			}

			vSpeciesAnnotbase = NULL;
			vSpeciesAnnotbase = (struct tagString **)calloc(nSpeciesNum, sizeof(struct tagString *));
			if(vSpeciesAnnotbase == NULL)
			{
				printf("Error: RefGene_MultiOrthologGetTSSAround, cannot load environment settings!\n");
				exit(EXIT_FAILURE);
			}

			vChrLen = NULL;
			vChrLen = (struct INTMATRIX **)calloc(nSpeciesNum, sizeof(struct INTMATRIX *));
			if(vChrLen == NULL)
			{
				printf("Error: RefGene_MultiOrthologGetTSSAround, cannot load environment settings!\n");
				exit(EXIT_FAILURE);
			}
		}
		else if(strstr(strRefLine, "[Species Name]") == strRefLine)
		{
			chp = strchr(strRefLine, '=');
			chp++;
			StrTrimLeft(chp);
			StringAddTail(vSpeciesName+ni, chp);
		}
		else if(strstr(strRefLine, "[Species Genome]") == strRefLine)
		{
			chp = strchr(strRefLine, '=');
			chp++;
			StrTrimLeft(chp);
			StringAddTail(vSpeciesSeqbase+ni, chp);
			sprintf(strRefLine0, "%schrlen.txt", chp);
			vChrLen[ni] = IMLOAD(strRefLine0);
			if(vChrLen[ni] == NULL)
			{
				printf("Error: RefGene_MultiOrthologGetTSSAround, cannot load environment settings!\n");
				exit(EXIT_FAILURE);
			}
		}
		else if(strstr(strRefLine, "[Species Conservation]") == strRefLine)
		{
			chp = strchr(strRefLine, '=');
			chp++;
			StrTrimLeft(chp);
			StringAddTail(vSpeciesConsbase+ni, chp);
		}
		else if(strstr(strRefLine, "[Species Annotation]") == strRefLine)
		{
			chp = strchr(strRefLine, '=');
			chp++;
			StrTrimLeft(chp);
			StringAddTail(vSpeciesAnnotbase+ni, chp);
			ni++;
		}
		else
		{
			printf("Error: RefGene_MultiOrthologGetTSSAround, unknown environment settings!\n");
			exit(EXIT_FAILURE);
		}

	}

	fclose(fpIn);

	if(ni != nSpeciesNum)
	{
		printf("Error: RefGene_MultiOrthologGetTSSAround, species number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* prepare coordinates */
	vSpeciesCord = NULL;
	vSpeciesCord = (struct tagString **)calloc(nSpeciesNum, sizeof(struct tagString *));
	if(vSpeciesCord == NULL)
	{
		printf("Error: RefGene_MultiOrthologGetTSSAround, cannot prepare environment settings!\n");
		exit(EXIT_FAILURE);
	}

	vSpeciesOut = NULL;
	vSpeciesOut = (struct tagString **)calloc(nSpeciesNum, sizeof(struct tagString *));
	if(vSpeciesOut == NULL)
	{
		printf("Error: RefGene_MultiOrthologGetTSSAround, cannot prepare environment settings!\n");
		exit(EXIT_FAILURE);
	}

	vfpOut = NULL;
	vfpOut = (FILE **)calloc(nSpeciesNum, sizeof(FILE *));
	if(vfpOut == NULL)
	{
		printf("Error: RefGene_MultiOrthologGetTSSAround, cannot prepare environment settings!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nSpeciesNum; ni++)
	{
		sprintf(strRefLine, "%s%s_%s.cod", strOutPath, strSeqFile, vSpeciesName[ni]->m_pString);
		StringAddTail(vSpeciesCord+ni, strRefLine);
		vfpOut[ni] = fopen(strRefLine, "w");
		if(vfpOut[ni] == NULL)
		{
			printf("Error: RefGene_MultiOrthologGetTSSAround, cannot prepare environment settings!\n");
			exit(EXIT_FAILURE);
		}
		sprintf(strRefLine, "%s_%s", strSeqFile, vSpeciesName[ni]->m_pString);
		StringAddTail(vSpeciesOut+ni, strRefLine);
	}

	/* load orthologs */
	fpIn = NULL;
	fpIn = fopen(strInPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: RefGene_MultiOrthologGetTSSAround, cannot load refgene!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strRefLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strRefLine);
		StrTrimRight(strRefLine);
		if(strRefLine[0] == '\0')
			continue;

		if(strRefLine[0] == '>')
		{
			chp = strRefLine+1;
			sprintf(strRefName, "%s", chp);

			for(nj=0; nj<nSpeciesNum; nj++)
			{
				fgets(strRefLine, LONG_LINE_LENGTH, fpIn);
				StrTrimLeft(strRefLine);
				StrTrimRight(strRefLine);
				chp = strRefLine;
				chp2 = strchr(chp, '\t');
				*chp2 = '\0';
				strcpy(strSpecies, chp);
				chp = chp2+1;
				chp2 = strchr(chp, '\t');
				*chp2 = '\0';
				nMatchCount = atoi(chp);
				chp = chp2+1;
				chp2 = strchr(chp, '\t');
				*chp2 = '\0';
				nOrienCount = atoi(chp);

				chp = chp2+1;
				strcpy(strRefLine0, chp);
				
				if(strcmp(strSpecies, vSpeciesName[nj]->m_pString) != 0)
				{
					printf("Error: RefGene_MultiOrthologGetTSSAround, species not match!\n");
					exit(EXIT_FAILURE);
				}

				if(strstr(strRefLine0, "---") == strRefLine0)
				{
				}
				else
				{
					pRefGene = NULL;
					pRefGene = RefGeneCreate();
					RefGeneInit_FromGenomeLabFormat(pRefGene, strRefLine0, strSpecies);
					
					if(pRefGene->chStrand == '-')
					{
						nPos1 = pRefGene->nTxEnd-nTSSDOWN;
						nPos2 = pRefGene->nTxEnd-nTSSUP;
					}
					else
					{
						nPos1 = pRefGene->nTxStart+nTSSUP;
						nPos2 = pRefGene->nTxStart+nTSSDOWN;
					}
					if(nPos1 < 0)
						nPos1 = 0;
					if(nPos2 < 0)
						nPos2 = 0;
				
					nChrid = pRefGene->nChrom-1;
					nChrLen = vChrLen[nj]->pMatElement[nChrid];
					if(nPos1 >= nChrLen)
						nPos1 = nChrLen-1;
					if(nPos2 >= nChrLen)
						nPos2 = nChrLen-1;

					if(nPos1 > nPos2)
					{
						printf("Error: RefGene_GetTargetTSSAround, direction not correct (start>end?)!\n");
						exit(EXIT_FAILURE);
					}

					fprintf(vfpOut[nj], "%s_%s_%s\t%d\t%d\t%d\t%c\n", strRefName, vSpeciesName[nj]->m_pString, 
						pRefGene->strName, pRefGene->nChrom, nPos1, nPos2, pRefGene->chStrand);

					RefGeneDestroy(pRefGene);
				}
			}
			ni++;
		}
	}

	fclose(fpIn);

	/* close coordinate file */
	for(ni=0; ni<nSpeciesNum; ni++)
	{
		fclose(vfpOut[ni]);
	}
	free(vfpOut);

	/* get sequence */
	for(ni=0; ni<nSpeciesNum; ni++)
	{
		/* sequence only */
		if(strcmp(vSpeciesConsbase[ni]->m_pString, "NULL") == 0)
		{
			sprintf(strRefLine, "%s%s.fa", strOutPath, vSpeciesOut[ni]->m_pString);
			Genome_Code_4bit_GetSeq_Main(vSpeciesSeqbase[ni]->m_pString, vSpeciesName[ni]->m_pString, 
				vSpeciesCord[ni]->m_pString, strRefLine, nStrandType);
		}
		/* sequence + conservation */
		else
		{
			Genome_Code_4bit_GetSeqCS_Main(vSpeciesSeqbase[ni]->m_pString, vSpeciesConsbase[ni]->m_pString,
				vSpeciesName[ni]->m_pString, vSpeciesCord[ni]->m_pString, strOutPath, 
				vSpeciesOut[ni]->m_pString, nStrandType, strConservationType);
		}
	}

	/* release memory */
	for(ni=0; ni<nSpeciesNum; ni++)
	{
		DeleteString(vSpeciesName[ni]);
		DeleteString(vSpeciesSeqbase[ni]);
		DeleteString(vSpeciesConsbase[ni]);
		DeleteString(vSpeciesAnnotbase[ni]);
		DestroyIntMatrix(vChrLen[ni]);
		DeleteString(vSpeciesCord[ni]);
		DeleteString(vSpeciesOut[ni]);
	}
	free(vSpeciesName);
	free(vSpeciesSeqbase);
	free(vSpeciesConsbase);
	free(vSpeciesAnnotbase);
	free(vChrLen);
	free(vSpeciesCord);
	free(vSpeciesOut);
	

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGene_MultiOrthologGetTSSAroundExonMasked()                          */
/*  Get TSS UP and Down sequence for multiortholog, exons will be masked   */
/*  in little cases.                                                       */
/* ----------------------------------------------------------------------- */ 
int RefGene_MultiOrthologGetTSSAroundExonMasked(char strInPath[], int nTSSUP, int nTSSDOWN,
									  char strParamPath[], char strOutPath[],
									  char strSeqFile[])
{
	/* define */
	
	/* environment */
	int nSpeciesNum;
	struct tagString **vSpeciesName;
	struct tagString **vSpeciesSeqbase;
	struct tagString **vSpeciesConsbase;
	struct tagString **vSpeciesAnnotbase;
	struct INTMATRIX **vChrLen;
	struct tagString **vSpeciesCord;
	struct tagString **vSpeciesOut;
	int ni,nj;
	int nPos1,nPos2;
	int nChrid,nChrLen;

	/* ortholog */
	int nMatchCount,nOrienCount;
	struct tagRefGene *pRefGene;
	char strSpecies[LINE_LENGTH];

	char strRefName[LINE_LENGTH];
	char strStrandType[LINE_LENGTH];
	int nStrandType;
	char strConservationType[LINE_LENGTH];

	char *chp,*chp2;
	char strRefLine[LONG_LINE_LENGTH];
	char strRefLine0[LONG_LINE_LENGTH];
	
	FILE *fpIn;
	FILE **vfpOut;
	FILE **vfpSeqOut;

	/* get seq & cs */
	char strSeqInFile[LINE_LENGTH];
	char strConsFile[LINE_LENGTH];
	char strChrName[LINE_LENGTH];
	char strAlias[LINE_LENGTH];
	
	int nStart,nEnd;
	struct tagSequence *pSeq;
	struct BYTEMATRIX *pCS;
	int nx,ny,nz,nl;
	int nRealStart;


	/* init environment */
	fpIn = NULL;
	fpIn = fopen(strParamPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: RefGene_MultiOrthologGetTSSAround, cannot load environment settings!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strRefLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strRefLine);
		StrTrimRight(strRefLine);
		if(strRefLine[0] == '\0')
			continue;

		if(strstr(strRefLine, "[Strand Type]") == strRefLine)
		{
			chp = strchr(strRefLine, '=');
			chp++;
			StrTrimLeft(chp);
			sprintf(strStrandType, "%s", chp);
			StrMakeUpper(strStrandType);
			if(strcmp(strStrandType, "GENEWISE") == 0)
			{
				nStrandType = 1;
			}
			else
			{
				nStrandType = 0;
			}
		}
		else if(strstr(strRefLine, "[Conservation Type]") == strRefLine)
		{
			chp = strchr(strRefLine, '=');
			chp++;
			StrTrimLeft(chp);
			sprintf(strConservationType, "%s", chp);
		}
		else if(strstr(strRefLine, "[Species Number]") == strRefLine)
		{
			chp = strchr(strRefLine, '=');
			chp++;
			StrTrimLeft(chp);
			nSpeciesNum = atoi(chp);

			if(nSpeciesNum <= 0)
			{
				printf("Warning: RefGene_MultiOrthologGetTSSAround, no species!\n");
				return PROC_SUCCESS;
			}

			vSpeciesName = NULL;
			vSpeciesName = (struct tagString **)calloc(nSpeciesNum, sizeof(struct tagString *));
			if(vSpeciesName == NULL)
			{
				printf("Error: RefGene_MultiOrthologGetTSSAround, cannot load environment settings!\n");
				exit(EXIT_FAILURE);
			}

			vSpeciesSeqbase = NULL;
			vSpeciesSeqbase = (struct tagString **)calloc(nSpeciesNum, sizeof(struct tagString *));
			if(vSpeciesSeqbase == NULL)
			{
				printf("Error: RefGene_MultiOrthologGetTSSAround, cannot load environment settings!\n");
				exit(EXIT_FAILURE);
			}

			vSpeciesConsbase = NULL;
			vSpeciesConsbase = (struct tagString **)calloc(nSpeciesNum, sizeof(struct tagString *));
			if(vSpeciesConsbase == NULL)
			{
				printf("Error: RefGene_MultiOrthologGetTSSAround, cannot load environment settings!\n");
				exit(EXIT_FAILURE);
			}

			vSpeciesAnnotbase = NULL;
			vSpeciesAnnotbase = (struct tagString **)calloc(nSpeciesNum, sizeof(struct tagString *));
			if(vSpeciesAnnotbase == NULL)
			{
				printf("Error: RefGene_MultiOrthologGetTSSAround, cannot load environment settings!\n");
				exit(EXIT_FAILURE);
			}

			vChrLen = NULL;
			vChrLen = (struct INTMATRIX **)calloc(nSpeciesNum, sizeof(struct INTMATRIX *));
			if(vChrLen == NULL)
			{
				printf("Error: RefGene_MultiOrthologGetTSSAround, cannot load environment settings!\n");
				exit(EXIT_FAILURE);
			}
		}
		else if(strstr(strRefLine, "[Species Name]") == strRefLine)
		{
			chp = strchr(strRefLine, '=');
			chp++;
			StrTrimLeft(chp);
			StringAddTail(vSpeciesName+ni, chp);
		}
		else if(strstr(strRefLine, "[Species Genome]") == strRefLine)
		{
			chp = strchr(strRefLine, '=');
			chp++;
			StrTrimLeft(chp);
			StringAddTail(vSpeciesSeqbase+ni, chp);
			sprintf(strRefLine0, "%schrlen.txt", chp);
			vChrLen[ni] = IMLOAD(strRefLine0);
			if(vChrLen[ni] == NULL)
			{
				printf("Error: RefGene_MultiOrthologGetTSSAround, cannot load environment settings!\n");
				exit(EXIT_FAILURE);
			}
		}
		else if(strstr(strRefLine, "[Species Conservation]") == strRefLine)
		{
			chp = strchr(strRefLine, '=');
			chp++;
			StrTrimLeft(chp);
			StringAddTail(vSpeciesConsbase+ni, chp);
		}
		else if(strstr(strRefLine, "[Species Annotation]") == strRefLine)
		{
			chp = strchr(strRefLine, '=');
			chp++;
			StrTrimLeft(chp);
			StringAddTail(vSpeciesAnnotbase+ni, chp);
			ni++;
		}
		else
		{
			printf("Error: RefGene_MultiOrthologGetTSSAround, unknown environment settings!\n");
			exit(EXIT_FAILURE);
		}

	}

	fclose(fpIn);

	if(ni != nSpeciesNum)
	{
		printf("Error: RefGene_MultiOrthologGetTSSAround, species number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* prepare coordinates */
	vSpeciesCord = NULL;
	vSpeciesCord = (struct tagString **)calloc(nSpeciesNum, sizeof(struct tagString *));
	if(vSpeciesCord == NULL)
	{
		printf("Error: RefGene_MultiOrthologGetTSSAround, cannot prepare environment settings!\n");
		exit(EXIT_FAILURE);
	}

	vSpeciesOut = NULL;
	vSpeciesOut = (struct tagString **)calloc(nSpeciesNum, sizeof(struct tagString *));
	if(vSpeciesOut == NULL)
	{
		printf("Error: RefGene_MultiOrthologGetTSSAround, cannot prepare environment settings!\n");
		exit(EXIT_FAILURE);
	}

	vfpOut = NULL;
	vfpOut = (FILE **)calloc(nSpeciesNum, sizeof(FILE *));
	if(vfpOut == NULL)
	{
		printf("Error: RefGene_MultiOrthologGetTSSAround, cannot prepare environment settings!\n");
		exit(EXIT_FAILURE);
	}

	vfpSeqOut = NULL;
	vfpSeqOut = (FILE **)calloc(nSpeciesNum, sizeof(FILE *));
	if(vfpSeqOut == NULL)
	{
		printf("Error: RefGene_MultiOrthologGetTSSAround, cannot prepare environment settings!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nSpeciesNum; ni++)
	{
		sprintf(strRefLine, "%s%s_%s.cod", strOutPath, strSeqFile, vSpeciesName[ni]->m_pString);
		StringAddTail(vSpeciesCord+ni, strRefLine);
		vfpOut[ni] = fopen(strRefLine, "w");
		if(vfpOut[ni] == NULL)
		{
			printf("Error: RefGene_MultiOrthologGetTSSAround, cannot prepare environment settings!\n");
			exit(EXIT_FAILURE);
		}

		sprintf(strRefLine, "%s_%s", strSeqFile, vSpeciesName[ni]->m_pString);
		StringAddTail(vSpeciesOut+ni, strRefLine);
		sprintf(strRefLine, "%s%s.fa", strOutPath, vSpeciesOut[ni]->m_pString);
		vfpSeqOut[ni] = fopen(strRefLine, "w");
		if(vfpSeqOut[ni] == NULL)
		{
			printf("Error: RefGene_MultiOrthologGetTSSAround, cannot prepare environment settings!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* load orthologs */
	fpIn = NULL;
	fpIn = fopen(strInPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: RefGene_MultiOrthologGetTSSAround, cannot load refgene!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strRefLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strRefLine);
		StrTrimRight(strRefLine);
		if(strRefLine[0] == '\0')
			continue;

		if(strRefLine[0] == '>')
		{
			chp = strRefLine+1;
			sprintf(strRefName, "%s", chp);

			for(nj=0; nj<nSpeciesNum; nj++)
			{
				fgets(strRefLine, LONG_LINE_LENGTH, fpIn);
				StrTrimLeft(strRefLine);
				StrTrimRight(strRefLine);
				chp = strRefLine;
				chp2 = strchr(chp, '\t');
				*chp2 = '\0';
				strcpy(strSpecies, chp);
				chp = chp2+1;
				chp2 = strchr(chp, '\t');
				*chp2 = '\0';
				nMatchCount = atoi(chp);
				chp = chp2+1;
				chp2 = strchr(chp, '\t');
				*chp2 = '\0';
				nOrienCount = atoi(chp);

				chp = chp2+1;
				strcpy(strRefLine0, chp);
				
				if(strcmp(strSpecies, vSpeciesName[nj]->m_pString) != 0)
				{
					printf("Error: RefGene_MultiOrthologGetTSSAround, species not match!\n");
					exit(EXIT_FAILURE);
				}

				if(strstr(strRefLine0, "---") == strRefLine0)
				{
					fprintf(vfpOut[nj], "%s_%s_NA\t-1\t-1\t-1\t?\n", strRefName, vSpeciesName[nj]->m_pString);
					sprintf(strAlias, "%s_%s_NA", strRefName, vSpeciesName[nj]->m_pString);
					fprintf(vfpSeqOut[nj], ">%s\n", strAlias);
					fprintf(vfpSeqOut[nj], "N\n");

					if(strcmp(vSpeciesConsbase[nj]->m_pString, "NULL") != 0)
					{
						pCS = NULL;
						pCS = CreateByteMatrix(1,1);
						if(pCS == NULL)
							continue;

						pCS->pMatElement[0] = 0;
						
						if(strcmp(strConservationType, "cs") == 0)
						{
							sprintf(strConsFile, "%s%s_%d_%s.cs", strOutPath, vSpeciesOut[nj]->m_pString, ni, strAlias);
							ConsScoreWriteToBinaryFile_ByStrand(pCS, strConsFile, '+');
						}
						else if(strcmp(strConservationType, "txt") == 0)
						{
							sprintf(strConsFile, "%s%s_%d_%s.txt", strOutPath, vSpeciesOut[nj]->m_pString, ni, strAlias);
							ConsScoreWriteToTextFile_ByStrand(pCS, strConsFile, '+');
						}
						else if(strcmp(strConservationType, "bed") == 0)
						{
						}

						DestroyByteMatrix(pCS);
					}
				}
				else
				{
					pRefGene = NULL;
					pRefGene = RefGeneCreate();
					RefGeneInit_FromGenomeLabFormat(pRefGene, strRefLine0, strSpecies);

					
					if(pRefGene->chStrand == '-')
					{
						nPos1 = pRefGene->nTxEnd-nTSSDOWN;
						nPos2 = pRefGene->nTxEnd-nTSSUP;
						nRealStart = pRefGene->nTxEnd;
					}
					else
					{
						nPos1 = pRefGene->nTxStart+nTSSUP;
						nPos2 = pRefGene->nTxStart+nTSSDOWN;
						nRealStart = pRefGene->nTxStart;
					}

					nChrid = pRefGene->nChrom-1;
					nChrLen = vChrLen[nj]->pMatElement[nChrid];

					if(nPos1 < 0)
						nPos1 = 0;
					if(nPos2 < 0)
						nPos2 = 0;
					/* if(nRealStart < 0)
						nRealStart = 0;
					*/

					if(nPos1 >= nChrLen)
						nPos1 = nChrLen-1;
					if(nPos2 >= nChrLen)
						nPos2 = nChrLen-1;
					/* if(nRealStart >= nChrLen)
						nRealStart = nChrLen-1;
					*/

					nRealStart -= nPos1;

					
					if(nPos1 > nPos2)
					{
						printf("Error: RefGene_GetTargetTSSAround, direction not correct (start>end?)!\n");
						exit(EXIT_FAILURE);
					}

					sprintf(strAlias, "%s_%s_%s", strRefName, vSpeciesName[nj]->m_pString, pRefGene->strName);
					fprintf(vfpOut[nj], "%s\t%d\t%d\t%d\t%c\n",  strAlias, pRefGene->nChrom, nPos1, nPos2, pRefGene->chStrand);

					nStart = nPos1;
					nEnd = nPos2;

					Genome_Index_To_ChromosomeName(strChrName, vSpeciesName[nj]->m_pString, pRefGene->nChrom);
					sprintf(strSeqInFile, "%s%s.sq", vSpeciesSeqbase[nj]->m_pString, strChrName);
					pSeq = NULL;
					pSeq = Genome_Code_4bit_GetSeq(strSeqInFile, nStart, nEnd);
					if(pSeq == NULL)
						continue;

					/* TODO: adjust exons */
					if(pRefGene->chStrand == '-')
					{
						for(nx=pRefGene->nExonCount-1; nx>=0; nx--)
						{
							ny = IMGETAT(pRefGene->pmatExonStartsEnds, nx, 0)-pRefGene->nTxEnd;
							nz = IMGETAT(pRefGene->pmatExonStartsEnds, nx, 1)-pRefGene->nTxEnd;
							
							ny = ny+nRealStart;
							nz = nz+nRealStart;

							if(ny < 0)
								ny = 0;
							if(nz < 0)
								nz = -1;

							if(ny >= pSeq->m_nLength)
								ny = pSeq->m_nLength;
							if(nz >= pSeq->m_nLength)
								nz = pSeq->m_nLength-1;

							for(nl=ny; nl<=nz; nl++)
							{
								switch(pSeq->m_pSequence->m_pString[nl])
								{
									case 'A': pSeq->m_pSequence->m_pString[nl] = 'a';
										break;
									case 'C': pSeq->m_pSequence->m_pString[nl] = 'c';
										break;
									case 'G': pSeq->m_pSequence->m_pString[nl] = 'g';
										break;
									case 'T': pSeq->m_pSequence->m_pString[nl] = 't';
										break;
									default: break;
								}
							}
						}
					}
					else
					{
						for(nx=0; nx<pRefGene->nExonCount; nx++)
						{
							ny = IMGETAT(pRefGene->pmatExonStartsEnds, nx, 0)-pRefGene->nTxStart;
							nz = IMGETAT(pRefGene->pmatExonStartsEnds, nx, 1)-pRefGene->nTxStart;
							
							ny = ny+nRealStart;
							nz = nz+nRealStart;
							
							if(ny < 0)
								ny = 0;
							if(nz < 0)
								nz = -1;

							if(ny >= pSeq->m_nLength)
								ny = pSeq->m_nLength;
							if(nz >= pSeq->m_nLength)
								nz = pSeq->m_nLength-1;

							for(nl=ny; nl<=nz; nl++)
							{
								switch(pSeq->m_pSequence->m_pString[nl])
								{
									case 'A': pSeq->m_pSequence->m_pString[nl] = 'a';
										break;
									case 'C': pSeq->m_pSequence->m_pString[nl] = 'c';
										break;
									case 'G': pSeq->m_pSequence->m_pString[nl] = 'g';
										break;
									case 'T': pSeq->m_pSequence->m_pString[nl] = 't';
										break;
									default: break;
								}
							}
						}
					}

					pSeq->m_nIndex = ni;
					strcpy(pSeq->m_strAlias, strAlias);
					if(nStrandType == 1)
						SequenceWriteToFasta_ByStrand(pSeq, vfpSeqOut[nj], pRefGene->chStrand, 1);
					else
						SequenceWriteToFasta_ByStrand(pSeq, vfpSeqOut[nj], '+', 1);
					SequenceDelete(pSeq);

					if(strcmp(vSpeciesConsbase[nj]->m_pString, "NULL") != 0)
					{
						sprintf(strConsFile, "%s%s.cs", vSpeciesConsbase[nj]->m_pString, strChrName);
						pCS = NULL;
						pCS = Genome_GetConsScore(strConsFile, nStart, nEnd);
						if(pCS == NULL)
							continue;

				
						if(nStrandType == 1)
						{
							if(strcmp(strConservationType, "cs") == 0)
							{
								sprintf(strConsFile, "%s%s_%d_%s.cs", strOutPath, vSpeciesOut[nj]->m_pString, ni, strAlias);
								ConsScoreWriteToBinaryFile_ByStrand(pCS, strConsFile, pRefGene->chStrand);
							}
							else if(strcmp(strConservationType, "txt") == 0)
							{
								sprintf(strConsFile, "%s%s_%d_%s.txt", strOutPath, vSpeciesOut[nj]->m_pString, ni, strAlias);
								ConsScoreWriteToTextFile_ByStrand(pCS, strConsFile, pRefGene->chStrand);
							}
							else if(strcmp(strConservationType, "bed") == 0)
							{
								sprintf(strConsFile, "%s%s_%d_%s.bed", strOutPath, vSpeciesOut[nj]->m_pString, ni, strAlias);
								ConsScoreWriteToBedFile_ByStrand(strChrName, nStart, nEnd, pCS, strConsFile, pRefGene->chStrand);
							}
						}
						else
						{
							if(strcmp(strConservationType, "cs") == 0)
							{
								sprintf(strConsFile, "%s%s_%d_%s.cs", strOutPath, vSpeciesOut[nj]->m_pString, ni, strAlias);
								ConsScoreWriteToBinaryFile_ByStrand(pCS, strConsFile, '+');
							}
							else if(strcmp(strConservationType, "txt") == 0)
							{
								sprintf(strConsFile, "%s%s_%d_%s.txt", strOutPath, vSpeciesOut[nj]->m_pString, ni, strAlias);
								ConsScoreWriteToTextFile_ByStrand(pCS, strConsFile, '+');
							}
							else if(strcmp(strConservationType, "bed") == 0)
							{
								sprintf(strConsFile, "%s%s_%d_%s.bed", strOutPath, vSpeciesOut[nj]->m_pString, ni, strAlias);
								ConsScoreWriteToBedFile_ByStrand(strChrName, nStart, nEnd, pCS, strConsFile, '+');
							}
						}

						DestroyByteMatrix(pCS);
					}

					
					/* destroy refgene */
					RefGeneDestroy(pRefGene);
				}
			}
			ni++;
		}
	}

	fclose(fpIn);

	/* close coordinate file */
	for(ni=0; ni<nSpeciesNum; ni++)
	{
		fclose(vfpOut[ni]);
		fclose(vfpSeqOut[ni]);
	}
	free(vfpOut);
	free(vfpSeqOut);

	/* release memory */
	for(ni=0; ni<nSpeciesNum; ni++)
	{
		DeleteString(vSpeciesName[ni]);
		DeleteString(vSpeciesSeqbase[ni]);
		DeleteString(vSpeciesConsbase[ni]);
		DeleteString(vSpeciesAnnotbase[ni]);
		DestroyIntMatrix(vChrLen[ni]);
		DeleteString(vSpeciesCord[ni]);
		DeleteString(vSpeciesOut[ni]);
	}
	free(vSpeciesName);
	free(vSpeciesSeqbase);
	free(vSpeciesConsbase);
	free(vSpeciesAnnotbase);
	free(vChrLen);
	free(vSpeciesCord);
	free(vSpeciesOut);
	

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGene_PickSpeciesSpecific()                                          */
/*  Pick up specified refgenes from a database.                            */
/*  This can be used to create species-species map using xenoRefGene map   */
/* ----------------------------------------------------------------------- */ 
int RefGene_PickSpeciesSpecific(char strTargetPath[], char strDatabasePath[], 
								char strOutPath[])
{
	/* define */
	char strTempPath[MED_LINE_LENGTH];
	char strTempPath2[MED_LINE_LENGTH];
	char strTempPath3[MED_LINE_LENGTH];
	FILE *fpIn;
	FILE *fpOut;
	FILE *fpTemp;
	char strLine[LONG_LINE_LENGTH];
	char strLine2[LONG_LINE_LENGTH];
	char strCommand[MED_LINE_LENGTH];
		
	int nWrite,nType;
	char strAlias[LINE_LENGTH];
	char strRefAlias[LINE_LENGTH];
	char *chSep;

	/* Step I: Merge database */
	sprintf(strTempPath, "%s.tmp", strOutPath);
	fpOut = NULL;
	fpOut = fopen(strTempPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: RefGene_PickSpeciesSpecific, cannot open temporary merge file!\n");
		exit(EXIT_FAILURE);
	}
	
	fpIn = NULL;
	fpIn = fopen(strTargetPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: RefGene_PickSpeciesSpecific, cannot open target list!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		fprintf(fpOut, "1\t%s\n", strLine);
	}

	fclose(fpIn);

	fpIn = NULL;
	fpIn = fopen(strDatabasePath, "r");
	if(fpIn == NULL)
	{
		printf("Error: RefGene_PickSpeciesSpecific, cannot open database list!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		fprintf(fpOut, "2\t%s\n", strLine);
	}

	fclose(fpIn);
	fclose(fpOut);

	/* Step II: Merge database */
	sprintf(strTempPath2, "%s.sort", strOutPath);
	sprintf(strCommand, "sort -o %s -d +1 %s\n", strTempPath2, strTempPath);
	system(strCommand);

	/* step III: pick up refgene */
	fpIn = NULL;
	fpIn = fopen(strTempPath2, "r");
	if(fpIn == NULL)
	{
		printf("Error: RefGene_PickSpeciesSpecific, cannot open sorted and combined refgene list!\n");
		exit(EXIT_FAILURE);
	}
	
	sprintf(strTempPath3, "%s.out", strOutPath);
	fpOut = NULL;
	fpOut = fopen(strTempPath3, "w");
	if(fpOut == NULL)
	{
		printf("Error: RefGene_PickSpeciesSpecific, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fpTemp = NULL;
	fpTemp = fopen(strTempPath, "w");
	if(fpTemp == NULL)
	{
		printf("Error: RefGene_PickSpeciesSpecific, cannot open temp output file!\n");
		exit(EXIT_FAILURE);
	}

	/* write */
	nWrite = 0;
	strcpy(strRefAlias, "");
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%d %s", &nType, strAlias);
		chSep = strchr(strLine, '\t');
		chSep++;
		
		if(strcmp(strAlias, strRefAlias) == 0)
		{
			/* do nothing */
		}
		else
		{
			/* write last */
			fclose(fpTemp);

			if(nWrite == 1)
			{
				fpTemp = NULL;
				fpTemp = fopen(strTempPath, "r");
				if(fpTemp == NULL)
				{
					printf("Error: RefGene_PickSpeciesSpecific, cannot open temp output file!\n");
					exit(EXIT_FAILURE);
				}

				while(fgets(strLine2, LONG_LINE_LENGTH, fpTemp) != NULL)
				{
					fprintf(fpOut, "%s", strLine2);
				}

				fclose(fpTemp);
			}

			fpTemp = NULL;
			fpTemp = fopen(strTempPath, "w");
			if(fpTemp == NULL)
			{
				printf("Error: RefGene_PickSpeciesSpecific, cannot open temp output file!\n");
				exit(EXIT_FAILURE);
			}
			nWrite = 0;
			strcpy(strRefAlias, strAlias);	
		}

		if(nType == 1)
		{
			nWrite = 1;
		}
		else
		{
			fprintf(fpTemp, "%s\n", chSep);
		}
	}

	if(nWrite == 1)
	{
		fpTemp = NULL;
		fpTemp = fopen(strTempPath, "r");
		if(fpTemp == NULL)
		{
			printf("Error: RefGene_PickSpeciesSpecific, cannot open temp output file!\n");
			exit(EXIT_FAILURE);
		}

		while(fgets(strLine2, LONG_LINE_LENGTH, fpTemp) != NULL)
		{
			fprintf(fpOut, "%s", strLine2);
		}

		fclose(fpTemp);
	}

	fclose(fpIn);
	fclose(fpOut);

	sprintf(strCommand, "sort -o %s -n +1 +3 +4 +5 +6 %s", strOutPath, strTempPath3);
	system(strCommand);

	sprintf(strCommand, "rm %s\n", strTempPath);
	system(strCommand);
	sprintf(strCommand, "rm %s\n", strTempPath2);
	system(strCommand);
	sprintf(strCommand, "rm %s\n", strTempPath3);
	system(strCommand);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGene_GetAffy_Main()                                                 */
/*  Link refgene id to affy probeset id.                                   */
/* ----------------------------------------------------------------------- */ 
int RefGene_GetAffy_Main(char strDatabasePath[], char strInputPath[], 
			int nColumn, char strOutputPath[])
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	char strLine[LONG_LINE_LENGTH];
	char strRefId[LINE_LENGTH];
	char *chp1,*chp2;
	int nj;
	int nPairNum = 0;
	int nFind;
	struct tagStringPair **vDataMap = NULL;

	/* load database */
	vDataMap = Affy_LoadDatabase_Affy2Refid(strDatabasePath, &nPairNum);

	/* init */
	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: RefGene_GetAffy_Main, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: RefGene_GetAffy_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	/* process one by one */
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		/* find refgene id */
		chp1 = strLine;
		for(nj=0; nj<nColumn; nj++)
		{
			chp2 = strchr(chp1, '\t');
			if(chp2 == NULL)
				break;

			chp1 = chp2+1;
		}

		if(nj != nColumn)
		{
			strcpy(strRefId, "---");	
		}
		else
		{
			sscanf(chp1, "%s", strRefId);

			if( (strcmp(strRefId, "---") == 0) || (strcmp(strRefId, "NA") == 0) )
			{
				strcpy(strRefId, "---");
			}
			else
			{
			}
		}

		/* search for affy */
		if(strcmp(strRefId, "---") == 0)
		{
			fprintf(fpOut, "---\t%s\n", strLine);
		}
		else
		{
			nFind = 0;
			for(nj=0; nj<nPairNum; nj++)
			{
				if(vDataMap[nj] != NULL)
				{
					if(strcmp(strRefId, vDataMap[nj]->m_pStr2->m_pString) == 0)
					{
						nFind = 1;
						fprintf(fpOut, "%s\t%s\n", vDataMap[nj]->m_pStr1->m_pString, strLine);
					}
				}
			}

			if(nFind == 0)
			{
				fprintf(fpOut, "---\t%s\n", strLine);
			}
		}
	}

	/* close files */
	fclose(fpIn);
	fclose(fpOut);

	/* clear memory */
	StringPair_ClearDatabase(&vDataMap, nPairNum);
	
	/* return */
	return PROC_SUCCESS;
}

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
			char strNetworkPath[], char strNetworkAnnotationPath[], int nNetDepth)
{
	/* define */
	int nRefGeneNum = 0;
	struct tagRefGene **vRefGeneDatabase = NULL;
	struct INTMATRIX *pMap = NULL;
	struct tagString **vValidationMap = NULL;
	struct BYTEMATRIX *pHasValue = NULL;
	struct DOUBLEMATRIX *pSigValue = NULL;

	int nValueNum = 0;
	int nInconsistencyNum = 0;
	struct DOUBLEMATRIX *pValue = NULL;
	struct tagString **vTag = NULL;
	double dMean,dSD;
	int nTransformType = 0;

	/* Network */
	struct NETWORKPHYSICAL *pNet = NULL;

	FILE *fpIn;
	FILE *fpOut;
	char strLine[LONG_LINE_LENGTH];
	int ni,nj,nId;
	char strTemp[MED_LINE_LENGTH];
	double dTemp;
	char strLocTransform[LINE_LENGTH];

	/* load refgene database */
	vRefGeneDatabase = RefGene_LoadDatabase(strDatabasePath, 2, strSpecies, &nRefGeneNum);
	if( (vRefGeneDatabase == NULL) || (nRefGeneNum <= 0) )
	{
		printf("Error: RefGene_Database_AssignValues_Main, cannot load refgene database!\n");
		exit(EXIT_FAILURE);
	}

	/* load map index */
	pMap = IMLOAD(strMapPath);
	if(pMap == NULL)
	{
		printf("Error: RefGene_Database_AssignValues_Main, cannot load map file!\n");
		exit(EXIT_FAILURE);
	}
	if(pMap->nHeight != nRefGeneNum)
	{
		printf("Error: RefGene_Database_AssignValues_Main, map/database dimensions do not match!\n");
		exit(EXIT_FAILURE);
	}

	fpIn = NULL;
	fpIn = fopen(strMapValidationPath, "r");
	if(fpIn != NULL)
	{
		vValidationMap = (struct tagString **)calloc(nRefGeneNum, sizeof(struct tagString *));
		if(vValidationMap == NULL)
		{
			printf("Error: RefGene_Database_AssignValues_Main, cannot create memory for loading validation index!\n");
			exit(EXIT_FAILURE);
		}

		ni = 0;
		while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;

			if(ni >= nRefGeneNum)
			{
				printf("Error: RefGene_Database_AssignValues_Main, validationmap/database dimensions do not match!\n");
				exit(EXIT_FAILURE);
			}

			StringAddTail(vValidationMap+ni, strLine);

			ni++;
		}

		if(ni != nRefGeneNum)
		{
			printf("Error: RefGene_Database_AssignValues_Main, validationmap/database dimensions do not match!\n");
			exit(EXIT_FAILURE);
		}

		fclose(fpIn);
	}
	else
	{
		printf("Proceed with no index map validation\n"); 
	}

	/* load significance values */
	fpIn = NULL;
	fpIn = fopen(strValuePath, "r");
	if(fpIn == NULL)
	{
		printf("Error: RefGene_Database_AssignValues_Main, cannot open value file!\n");
		exit(EXIT_FAILURE);
	}

	nValueNum = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		nValueNum++;
	}
	fclose(fpIn);

	if(nValueNum <= 0)
	{
		printf("Error: RefGene_Database_AssignValues_Main, empty value file!\n");
		exit(EXIT_FAILURE);
	}

	vTag = (struct tagString **)calloc(nValueNum, sizeof(struct tagString *));
	if(vTag == NULL)
	{
		printf("Error: RefGene_Database_AssignValues_Main, cannot create memory for loading values!\n");
		exit(EXIT_FAILURE);
	}

	pValue = CreateDoubleMatrix(1, nValueNum);
	if(pValue == NULL)
	{
		printf("Error: RefGene_Database_AssignValues_Main, cannot create memory for loading values!\n");
		exit(EXIT_FAILURE);
	}

	fpIn = NULL;
	fpIn = fopen(strValuePath, "r");
	if(fpIn == NULL)
	{
		printf("Error: RefGene_Database_AssignValues_Main, cannot open value file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(ni >= nValueNum)
		{
			printf("Error: RefGene_Database_AssignValues_Main, value dimensions do not match!\n");
			exit(EXIT_FAILURE);
		}

		sscanf(strLine, "%s %lf", strTemp, &dTemp);
		StringAddTail(vTag+ni, strTemp);
		pValue->pMatElement[ni] = dTemp;

		ni++;
	}
	fclose(fpIn);

	if(ni != nValueNum)
	{
		printf("Error: RefGene_Database_AssignValues_Main, value dimensions do not match!\n");
		exit(EXIT_FAILURE);
	}

	/* preprocessing */
	nTransformType = 0;
	strcpy(strLocTransform, strTransform);
	StrMakeUpper(strLocTransform);
	if(strcmp(strLocTransform, "LOGIT") == 0)
		nTransformType = 1;
		
	/* preprocessing values */
	for(ni=0; ni<nValueNum; ni++)
	{
		if(pValue->pMatElement[ni] < dTruncateLowerBound)
			pValue->pMatElement[ni] = dTruncateLowerBound;
		if(pValue->pMatElement[ni] > dTruncateUpperBound)
			pValue->pMatElement[ni] = dTruncateUpperBound;
		if(nTransformType == 1)
		{
			pValue->pMatElement[ni] = log(pValue->pMatElement[ni]/(1.0-pValue->pMatElement[ni]));
		}
	}


	/* normalize values */
	if(nNormalize == 1)
	{
		dMean = 0.0;
		dSD = 0.0;
		for(ni=0; ni<nValueNum; ni++)
		{
			dMean += pValue->pMatElement[ni];
		}


		dMean /= (double)nValueNum;
		for(ni=0; ni<nValueNum; ni++)
		{
			dSD += (pValue->pMatElement[ni]-dMean)*(pValue->pMatElement[ni]-dMean);
		}
		if(nValueNum > 1)
		{
			dSD /= (double)(nValueNum-1);
			dSD = sqrt(dSD);
		}
		else
		{
			dSD = 1e-6;
		}

		for(ni=0; ni<nValueNum; ni++)
		{
			pValue->pMatElement[ni] = (pValue->pMatElement[ni]-dMean)/dSD;
		}
	}

	/* post processing */
	nTransformType = 0;
	strcpy(strLocTransform, strPostTransform);
	StrMakeUpper(strLocTransform);
	if(strcmp(strLocTransform, "-1") == 0)
		nTransformType = 2;

	for(ni=0; ni<nValueNum; ni++)
	{
		if(pValue->pMatElement[ni] < dPostLowerBound)
			pValue->pMatElement[ni] = dPostLowerBound;
		if(pValue->pMatElement[ni] > dPostUpperBound)
			pValue->pMatElement[ni] = dPostUpperBound;
		if(nTransformType == 2)
		{
			pValue->pMatElement[ni] = -pValue->pMatElement[ni];
		}
	}

	/* assign values */
	pHasValue = CreateByteMatrix(nRefGeneNum,1);
	pSigValue = CreateDoubleMatrix(nRefGeneNum,1);
	if( (pHasValue == NULL) || (pSigValue == NULL) )
	{
		printf("Error: RefGene_Database_AssignValues_Main, cannot create memory for assigning values!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nRefGeneNum; ni++)
	{
		nId = pMap->pMatElement[ni];
		if( (nId<0) || (nId>=nValueNum) )
		{
			pHasValue->pMatElement[ni] = 0;
		}
		else
		{
			pHasValue->pMatElement[ni] = 1;
			pSigValue->pMatElement[ni] = pValue->pMatElement[nId];
			if(vValidationMap != NULL)
			{
				if(strcmp(vValidationMap[ni]->m_pString, vTag[nId]->m_pString) != 0)
				{
					printf("Warning: %s (%s) mapping inconsistency: Lineid = %d, Mapid = %d, Maptag = %s, Realtag = %s\n",
						vRefGeneDatabase[ni]->strGene, vRefGeneDatabase[ni]->strName,
						ni, nId, vValidationMap[ni]->m_pString, vTag[nId]->m_pString);
					pHasValue->pMatElement[ni] = 0;
					nInconsistencyNum++;
				}
			}
		}
	}

	/* adjust values based on network structure */
	pNet = Network_LoadFromSIF(strNetworkPath);
	if(pNet != NULL)
	{
		/* annotate the network */
		Network_LoadAnnotation(pNet, strNetworkAnnotationPath);

		/* network prepare values */
		Network_InitNodeDV(pNet, 1, (nNetDepth+1));
		Network_InitNodeBV(pNet, 1, (nNetDepth+1));
		Network_InitNodeSV(pNet, (nNetDepth+1));

		/* network assign values */
		Network_LoadValue_FromRefGene(pNet, 0, 0, vRefGeneDatabase, nRefGeneNum, pHasValue, pSigValue, nTakeAbsoluteValue, 1);

		/* network find biggest value with N steps */
		Network_DepthSearch_GetMaxScore(pNet, 0, 0, nNetDepth);
	}

	/* export values */
	fpOut = NULL;
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: RefGene_Database_AssignValues_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nRefGeneNum; ni++)
	{
		fprintf(fpOut, "%d\t%f", pHasValue->pMatElement[ni], pSigValue->pMatElement[ni]);
		if(pNet != NULL)
		{
			nId = vRefGeneDatabase[ni]->nGeneID;
			if( (nId<0) || (nId>pNet->nMaxNodeNum) )
			{
				for(nj=0; nj<nNetDepth; nj++)
				{
					fprintf(fpOut, "\t---\t---");
				}
			}
			else if(pNet->vNodes[nId] == NULL)
			{
				for(nj=0; nj<nNetDepth; nj++)
				{
					fprintf(fpOut, "\t---\t---");
				}
			}
			else
			{
				for(nj=0; nj<nNetDepth; nj++)
				{
					if(BMGETAT(pNet->vNodes[nId]->pBV, 0, (nj+1)) == 0)
					{
                        fprintf(fpOut, "\t---\t---");
					}
					else
					{
						fprintf(fpOut, "\t%f\t%s", 
							DMGETAT(pNet->vNodes[nId]->pDV, 0, (nj+1))/(double)(nj+2),
							pNet->vNodes[nId]->vSV[nj+1]->m_pString);
					}
				}
			}
		}

		fprintf(fpOut, "\n");
	}

	fclose(fpOut);

	/* clear memory */
	RefGene_ClearDatabase(&vRefGeneDatabase, nRefGeneNum);
	DestroyIntMatrix(pMap);
	if(vValidationMap != NULL)
	{
		for(ni=0; ni<nRefGeneNum; ni++)
		{
			DeleteString(vValidationMap[ni]);
			vValidationMap[ni] = NULL;
		}
		free(vValidationMap);
	}
	DestroyDoubleMatrix(pValue);
	for(ni=0; ni<nValueNum; ni++)
	{
		DeleteString(vTag[ni]);
		vTag[ni] = NULL;
	}
	free(vTag);
	DestroyByteMatrix(pHasValue);
	DestroyDoubleMatrix(pSigValue);
	/* destroy net */
	if(pNet != NULL)
	{
		NETWORKPHYSICALDESTROY(&pNet);
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  RefGene_Database_AssignShortest_Main()                                 */
/*  Search shorted path to target genes for individual genes in a refgene  */
/*  database.                                                              */
/* ----------------------------------------------------------------------- */ 
int RefGene_Database_AssignShortest_Main(char strDatabasePath[], char strSpecies[], 
			char strOutputPath[], char strValuePath[], 
			char strMapPath[], char strMapValidationPath[],
			char strNetworkPath[], char strNetworkAnnotationPath[], 
			char strTargetPath[], int nMaxIterNum)
{
	/* define */
	int nRefGeneNum = 0;
	struct tagRefGene **vRefGeneDatabase = NULL;
	struct INTMATRIX *pMap = NULL;
	struct tagString **vValidationMap = NULL;
	struct BYTEMATRIX *pHasValue = NULL;
	struct DOUBLEMATRIX *pSigValue = NULL;

	int nValueNum = 0;
	int nInconsistencyNum = 0;
	struct DOUBLEMATRIX *pValue = NULL;
	struct tagString **vTag = NULL;
	double dMean,dSD;

	/* Network */
	struct NETWORKPHYSICAL *pNet = NULL;
	struct INTMATRIX *pTargetNodes = NULL;
	int nMaxIter;
	struct NETWORKSHORTESTPATH *pPath = NULL;
	double dMaxNodeValue = -log(2.0*normcdf(0.0, 1.0, 1.0)-1.0)/log(10.0);

	FILE *fpIn;
	FILE *fpOut;
	char strLine[LONG_LINE_LENGTH];
	int ni,nj,nk,nId,nNetId,nMinIndex;
	char strTemp[MED_LINE_LENGTH];
	double dTemp;
	double dMinDist;
	int nTemp;

	/* load refgene database */
	vRefGeneDatabase = RefGene_LoadDatabase(strDatabasePath, 2, strSpecies, &nRefGeneNum);
	if( (vRefGeneDatabase == NULL) || (nRefGeneNum <= 0) )
	{
		printf("Error: RefGene_Database_AssignValues_Main, cannot load refgene database!\n");
		exit(EXIT_FAILURE);
	}

	/* load map index */
	pMap = IMLOAD(strMapPath);
	if(pMap == NULL)
	{
		printf("Error: RefGene_Database_AssignValues_Main, cannot load map file!\n");
		exit(EXIT_FAILURE);
	}
	if(pMap->nHeight != nRefGeneNum)
	{
		printf("Error: RefGene_Database_AssignValues_Main, map/database dimensions do not match!\n");
		exit(EXIT_FAILURE);
	}

	fpIn = NULL;
	fpIn = fopen(strMapValidationPath, "r");
	if(fpIn != NULL)
	{
		vValidationMap = (struct tagString **)calloc(nRefGeneNum, sizeof(struct tagString *));
		if(vValidationMap == NULL)
		{
			printf("Error: RefGene_Database_AssignValues_Main, cannot create memory for loading validation index!\n");
			exit(EXIT_FAILURE);
		}

		ni = 0;
		while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;

			if(ni >= nRefGeneNum)
			{
				printf("Error: RefGene_Database_AssignValues_Main, validationmap/database dimensions do not match!\n");
				exit(EXIT_FAILURE);
			}

			StringAddTail(vValidationMap+ni, strLine);

			ni++;
		}

		if(ni != nRefGeneNum)
		{
			printf("Error: RefGene_Database_AssignValues_Main, validationmap/database dimensions do not match!\n");
			exit(EXIT_FAILURE);
		}

		fclose(fpIn);
	}
	else
	{
		printf("Proceed with no index map validation\n"); 
	}

	/* load significance values */
	fpIn = NULL;
	fpIn = fopen(strValuePath, "r");
	if(fpIn == NULL)
	{
		printf("Error: RefGene_Database_AssignValues_Main, cannot open value file!\n");
		exit(EXIT_FAILURE);
	}

	nValueNum = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		nValueNum++;
	}
	fclose(fpIn);

	if(nValueNum <= 0)
	{
		printf("Error: RefGene_Database_AssignValues_Main, empty value file!\n");
		exit(EXIT_FAILURE);
	}

	vTag = (struct tagString **)calloc(nValueNum, sizeof(struct tagString *));
	if(vTag == NULL)
	{
		printf("Error: RefGene_Database_AssignValues_Main, cannot create memory for loading values!\n");
		exit(EXIT_FAILURE);
	}

	pValue = CreateDoubleMatrix(1, nValueNum);
	if(pValue == NULL)
	{
		printf("Error: RefGene_Database_AssignValues_Main, cannot create memory for loading values!\n");
		exit(EXIT_FAILURE);
	}

	fpIn = NULL;
	fpIn = fopen(strValuePath, "r");
	if(fpIn == NULL)
	{
		printf("Error: RefGene_Database_AssignValues_Main, cannot open value file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(ni >= nValueNum)
		{
			printf("Error: RefGene_Database_AssignValues_Main, value dimensions do not match!\n");
			exit(EXIT_FAILURE);
		}

		sscanf(strLine, "%s %lf", strTemp, &dTemp);
		StringAddTail(vTag+ni, strTemp);
		pValue->pMatElement[ni] = dTemp;

		ni++;
	}
	fclose(fpIn);

	if(ni != nValueNum)
	{
		printf("Error: RefGene_Database_AssignValues_Main, value dimensions do not match!\n");
		exit(EXIT_FAILURE);
	}

	
	/* normalize values */
	dMean = 0.0;
	dSD = 0.0;
	for(ni=0; ni<nValueNum; ni++)
	{
		dMean += pValue->pMatElement[ni];
	}
	dMean /= (double)nValueNum;
	for(ni=0; ni<nValueNum; ni++)
	{
		dSD += (pValue->pMatElement[ni]-dMean)*(pValue->pMatElement[ni]-dMean);
	}
	if(nValueNum > 1)
	{
		dSD /= (double)(nValueNum-1);
		dSD = sqrt(dSD);
	}
	else
	{
		dSD = 1e-6;
	}
	for(ni=0; ni<nValueNum; ni++)
	{
		pValue->pMatElement[ni] = (pValue->pMatElement[ni]-dMean)/dSD;
	}

	/* assign values */
	pHasValue = CreateByteMatrix(nRefGeneNum,1);
	pSigValue = CreateDoubleMatrix(nRefGeneNum,1);
	if( (pHasValue == NULL) || (pSigValue == NULL) )
	{
		printf("Error: RefGene_Database_AssignValues_Main, cannot create memory for assigning values!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nRefGeneNum; ni++)
	{
		nId = pMap->pMatElement[ni];
		if( (nId<0) || (nId>=nValueNum) )
		{
			pHasValue->pMatElement[ni] = 0;
		}
		else
		{
			pHasValue->pMatElement[ni] = 1;
			pSigValue->pMatElement[ni] = pValue->pMatElement[nId];
			if(vValidationMap != NULL)
			{
				if(strcmp(vValidationMap[ni]->m_pString, vTag[nId]->m_pString) != 0)
				{
					printf("Warning: %s (%s) mapping inconsistency: Lineid = %d, Mapid = %d, Maptag = %s, Realtag = %s\n",
						vRefGeneDatabase[ni]->strGene, vRefGeneDatabase[ni]->strName,
						ni, nId, vValidationMap[ni]->m_pString, vTag[nId]->m_pString);
					pHasValue->pMatElement[ni] = 0;
					nInconsistencyNum++;
				}
			}
		}
	}

	/* adjust values based on network structure */
	pNet = Network_LoadFromSIF(strNetworkPath);
	if(pNet != NULL)
	{
		nMaxIter = pNet->nRealNodeNum;
		if(nMaxIter > nMaxIterNum)
			nMaxIter = nMaxIterNum;

		/* annotate the network */
		Network_LoadAnnotation(pNet, strNetworkAnnotationPath);

		/* network prepare values */
		Network_InitNodeDV(pNet, 1, 1);
		Network_InitNodeBV(pNet, 1, 1);
		
		/* network assign values */
		Network_LoadValue_FromRefGene_ForShortestPath(pNet, 0, 0, vRefGeneDatabase, nRefGeneNum, pHasValue, pSigValue, 1, 0, dMaxNodeValue);

		/* network find shortest path */
		pTargetNodes = IMLOAD(strTargetPath);
		if(pTargetNodes == NULL)
		{
			printf("Error: cannot load target nodes for shortest path search!\n");
			exit(EXIT_FAILURE);
		}
		if( (pTargetNodes->nHeight > 1) && (pTargetNodes->nWidth == 1) )
		{
			nTemp = pTargetNodes->nHeight;
			pTargetNodes->nHeight = pTargetNodes->nWidth;
			pTargetNodes->nWidth = nTemp;
		}
		pPath = Network_FindShortestPath_V(pNet, pTargetNodes, 0, 0, dMaxNodeValue, nMaxIter);
		if(pPath == NULL)
		{
			printf("Error: cannot successfully perform shortest path search!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* export values */
	fpOut = NULL;
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: RefGene_Database_AssignValues_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nRefGeneNum; ni++)
	{
		fprintf(fpOut, "%d\t%f", pHasValue->pMatElement[ni], pSigValue->pMatElement[ni]);
		if(pNet != NULL)
		{
			nId = vRefGeneDatabase[ni]->nGeneID;
			if( (nId>=0) && (nId<pNet->nMaxNodeNum) )
			{
				if(pNet->vNodes[nId] != NULL)
				{
					nNetId = pNet->vNodes[nId]->nNetIndex;
					dMinDist = dMaxNodeValue*pNet->nRealNodeNum+1.0;
					nMinIndex = -1;
					for(nk=0; nk<pPath->vShortDist[nNetId]->nWidth; nk++)
					{
						if(pPath->vShortDist[nNetId]->pMatElement[nk] < dMinDist)
						{
							dMinDist = pPath->vShortDist[nNetId]->pMatElement[nk];
							nMinIndex = nk;
						}
					}
					if( (nMinIndex < 0) || (dMinDist >=dMaxNodeValue*pNet->nRealNodeNum) )
					{
						fprintf(fpOut, "\t-1\t---");
					}
					else
					{
						fprintf(fpOut, "\t%f\t", dMinDist);
						if(pNet->vNodes[nId]->nTag >= 0)
						{
							fprintf(fpOut, "(%s,%f)", pNet->vNodes[nId]->strName->m_pString,
								DMGETAT(pNet->vNodes[nId]->pDV, 0, 0));
						}
						else
						{
							fprintf(fpOut, "%s,%f", pNet->vNodes[nId]->strName->m_pString,
								DMGETAT(pNet->vNodes[nId]->pDV, 0, 0));
						}

						nj = pPath->vNeighborNode[nNetId]->pMatElement[nMinIndex];
						while(nj != -1)
						{
							if(pNet->vNodes[nj]->nTag >= 0)
							{
								fprintf(fpOut, "-->(%s,%f)", pNet->vNodes[nj]->strName->m_pString,
									DMGETAT(pNet->vNodes[nj]->pDV, 0, 0));
							}
							else
							{
								fprintf(fpOut, "-->%s,%f", pNet->vNodes[nj]->strName->m_pString,
									DMGETAT(pNet->vNodes[nj]->pDV, 0, 0));
							}

							nNetId = pNet->vNodes[nj]->nNetIndex;
							nj = pPath->vNeighborNode[nNetId]->pMatElement[nMinIndex];
						}
					}
				}
				else
				{
					fprintf(fpOut, "\t-1\t---");
				}
			}
			else
			{
				fprintf(fpOut, "\t-1\t---");
			}
		}

		fprintf(fpOut, "\n");
	}

	fclose(fpOut);

	/* clear memory */
	RefGene_ClearDatabase(&vRefGeneDatabase, nRefGeneNum);
	DestroyIntMatrix(pMap);
	if(vValidationMap != NULL)
	{
		for(ni=0; ni<nRefGeneNum; ni++)
		{
			DeleteString(vValidationMap[ni]);
			vValidationMap[ni] = NULL;
		}
		free(vValidationMap);
	}
	DestroyDoubleMatrix(pValue);
	for(ni=0; ni<nValueNum; ni++)
	{
		DeleteString(vTag[ni]);
		vTag[ni] = NULL;
	}
	free(vTag);
	DestroyByteMatrix(pHasValue);
	DestroyDoubleMatrix(pSigValue);
	/* destroy net */
	if(pNet != NULL)
	{
		NETWORKPHYSICALDESTROY(&pNet);
		NETWORKSHORTESTPATHDESTROY(&pPath);
		DestroyIntMatrix(pTargetNodes);
	}
	
	/* return */
	return PROC_SUCCESS;
}



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
			char strNetworkPath[], char strNetworkAnnotationPath[], int nOperationType)
{
	/* define */
	int nRefGeneNum = 0;
	struct tagRefGene **vRefGeneDatabase = NULL;
	struct INTMATRIX *pMap = NULL;
	struct tagString **vValidationMap = NULL;
	struct BYTEMATRIX *pHasValue = NULL;
	struct DOUBLEMATRIX *pSigValue = NULL;

	int nValueNum = 0;
	int nInconsistencyNum = 0;
	struct DOUBLEMATRIX *pValue = NULL;
	struct tagString **vTag = NULL;
	double dMean,dSD;
	int nTransformType = 0;

	/* Network */
	struct NETWORKPHYSICAL *pNet = NULL;
	struct DOUBLEMATRIX *pEdgeScore = NULL;
	int nNodeID1,nNodeID2;
	int nHasValue;

	FILE *fpIn;
	FILE *fpOut;
	char strLine[LONG_LINE_LENGTH];
	int ni,nId;
	char strTemp[MED_LINE_LENGTH];
	double dTemp;
	char strLocTransform[LINE_LENGTH];

	/* load refgene database */
	vRefGeneDatabase = RefGene_LoadDatabase(strDatabasePath, 2, strSpecies, &nRefGeneNum);
	if( (vRefGeneDatabase == NULL) || (nRefGeneNum <= 0) )
	{
		printf("Error: RefGene_Database_AssignValues_Main, cannot load refgene database!\n");
		exit(EXIT_FAILURE);
	}

	/* load map index */
	pMap = IMLOAD(strMapPath);
	if(pMap == NULL)
	{
		printf("Error: RefGene_Database_AssignValues_Main, cannot load map file!\n");
		exit(EXIT_FAILURE);
	}
	if(pMap->nHeight != nRefGeneNum)
	{
		printf("Error: RefGene_Database_AssignValues_Main, map/database dimensions do not match!\n");
		exit(EXIT_FAILURE);
	}

	fpIn = NULL;
	fpIn = fopen(strMapValidationPath, "r");
	if(fpIn != NULL)
	{
		vValidationMap = (struct tagString **)calloc(nRefGeneNum, sizeof(struct tagString *));
		if(vValidationMap == NULL)
		{
			printf("Error: RefGene_Database_AssignValues_Main, cannot create memory for loading validation index!\n");
			exit(EXIT_FAILURE);
		}

		ni = 0;
		while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;

			if(ni >= nRefGeneNum)
			{
				printf("Error: RefGene_Database_AssignValues_Main, validationmap/database dimensions do not match!\n");
				exit(EXIT_FAILURE);
			}

			StringAddTail(vValidationMap+ni, strLine);

			ni++;
		}

		if(ni != nRefGeneNum)
		{
			printf("Error: RefGene_Database_AssignValues_Main, validationmap/database dimensions do not match!\n");
			exit(EXIT_FAILURE);
		}

		fclose(fpIn);
	}
	else
	{
		printf("Proceed with no index map validation\n"); 
	}

	/* load significance values */
	fpIn = NULL;
	fpIn = fopen(strValuePath, "r");
	if(fpIn == NULL)
	{
		printf("Error: RefGene_Database_AssignValues_Main, cannot open value file!\n");
		exit(EXIT_FAILURE);
	}

	nValueNum = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		nValueNum++;
	}
	fclose(fpIn);

	if(nValueNum <= 0)
	{
		printf("Error: RefGene_Database_AssignValues_Main, empty value file!\n");
		exit(EXIT_FAILURE);
	}

	vTag = (struct tagString **)calloc(nValueNum, sizeof(struct tagString *));
	if(vTag == NULL)
	{
		printf("Error: RefGene_Database_AssignValues_Main, cannot create memory for loading values!\n");
		exit(EXIT_FAILURE);
	}

	pValue = CreateDoubleMatrix(1, nValueNum);
	if(pValue == NULL)
	{
		printf("Error: RefGene_Database_AssignValues_Main, cannot create memory for loading values!\n");
		exit(EXIT_FAILURE);
	}

	fpIn = NULL;
	fpIn = fopen(strValuePath, "r");
	if(fpIn == NULL)
	{
		printf("Error: RefGene_Database_AssignValues_Main, cannot open value file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(ni >= nValueNum)
		{
			printf("Error: RefGene_Database_AssignValues_Main, value dimensions do not match!\n");
			exit(EXIT_FAILURE);
		}

		sscanf(strLine, "%s %lf", strTemp, &dTemp);
		StringAddTail(vTag+ni, strTemp);
		pValue->pMatElement[ni] = dTemp;

		ni++;
	}
	fclose(fpIn);

	if(ni != nValueNum)
	{
		printf("Error: RefGene_Database_AssignValues_Main, value dimensions do not match!\n");
		exit(EXIT_FAILURE);
	}

	/* preprocessing */
	nTransformType = 0;
	strcpy(strLocTransform, strTransform);
	StrMakeUpper(strLocTransform);
	if(strcmp(strLocTransform, "LOGIT") == 0)
		nTransformType = 1;
		
	/* preprocessing values */
	for(ni=0; ni<nValueNum; ni++)
	{
		if(pValue->pMatElement[ni] < dTruncateLowerBound)
			pValue->pMatElement[ni] = dTruncateLowerBound;
		if(pValue->pMatElement[ni] > dTruncateUpperBound)
			pValue->pMatElement[ni] = dTruncateUpperBound;
		if(nTransformType == 1)
		{
			pValue->pMatElement[ni] = log(pValue->pMatElement[ni]/(1.0-pValue->pMatElement[ni]));
		}
	}


	/* normalize values */
	if(nNormalize == 1)
	{
		dMean = 0.0;
		dSD = 0.0;
		for(ni=0; ni<nValueNum; ni++)
		{
			dMean += pValue->pMatElement[ni];
		}


		dMean /= (double)nValueNum;
		for(ni=0; ni<nValueNum; ni++)
		{
			dSD += (pValue->pMatElement[ni]-dMean)*(pValue->pMatElement[ni]-dMean);
		}
		if(nValueNum > 1)
		{
			dSD /= (double)(nValueNum-1);
			dSD = sqrt(dSD);
		}
		else
		{
			dSD = 1e-6;
		}

		for(ni=0; ni<nValueNum; ni++)
		{
			pValue->pMatElement[ni] = (pValue->pMatElement[ni]-dMean)/dSD;
		}
	}

	/* post processing */
	nTransformType = 0;
	strcpy(strLocTransform, strPostTransform);
	StrMakeUpper(strLocTransform);
	if(strcmp(strLocTransform, "-1") == 0)
		nTransformType = 2;

	for(ni=0; ni<nValueNum; ni++)
	{
		if(pValue->pMatElement[ni] < dPostLowerBound)
			pValue->pMatElement[ni] = dPostLowerBound;
		if(pValue->pMatElement[ni] > dPostUpperBound)
			pValue->pMatElement[ni] = dPostUpperBound;
		if(nTransformType == 2)
		{
			pValue->pMatElement[ni] = -pValue->pMatElement[ni];
		}
	}

	/* assign values */
	pHasValue = CreateByteMatrix(nRefGeneNum,1);
	pSigValue = CreateDoubleMatrix(nRefGeneNum,1);
	if( (pHasValue == NULL) || (pSigValue == NULL) )
	{
		printf("Error: RefGene_Database_AssignValues_Main, cannot create memory for assigning values!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nRefGeneNum; ni++)
	{
		nId = pMap->pMatElement[ni];
		if( (nId<0) || (nId>=nValueNum) )
		{
			pHasValue->pMatElement[ni] = 0;
		}
		else
		{
			pHasValue->pMatElement[ni] = 1;
			pSigValue->pMatElement[ni] = pValue->pMatElement[nId];
			if(vValidationMap != NULL)
			{
				if(strcmp(vValidationMap[ni]->m_pString, vTag[nId]->m_pString) != 0)
				{
					printf("Warning: %s (%s) mapping inconsistency: Lineid = %d, Mapid = %d, Maptag = %s, Realtag = %s\n",
						vRefGeneDatabase[ni]->strGene, vRefGeneDatabase[ni]->strName,
						ni, nId, vValidationMap[ni]->m_pString, vTag[nId]->m_pString);
					pHasValue->pMatElement[ni] = 0;
					nInconsistencyNum++;
				}
			}
		}
	}

	/* adjust values based on network structure */
	pNet = Network_LoadFromSIF(strNetworkPath);
	if(pNet == NULL)
	{
		printf("Error: RefGene_Database_AssignValues_Main, cannot load network!\n");
		exit(EXIT_FAILURE);
	}

	/* annotate the network */
	Network_LoadAnnotation(pNet, strNetworkAnnotationPath);

	/* network prepare values */
	Network_InitNodeDV(pNet, 1, 1);
	Network_InitNodeBV(pNet, 1, 1);
	
	/* network assign values */
	Network_LoadValue_FromRefGene(pNet, 0, 0, vRefGeneDatabase, nRefGeneNum, pHasValue, pSigValue, nTakeAbsoluteValue, 2);

	/* export values */
	fpOut = NULL;
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: RefGene_Database_AssignValues_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<pNet->nEdgeNum; ni++)
	{
		nNodeID1 = pNet->vEdges[ni]->nNodeID1;
		nNodeID2 = pNet->vEdges[ni]->nNodeID2;
		if(pNet->vEdges[ni]->nInteractionType == 1)
		{
			strcpy(strTemp, "pp");
		}
		else if(pNet->vEdges[ni]->nInteractionType == 2)
		{
			strcpy(strTemp, "pd");
		}
		else if(pNet->vEdges[ni]->nInteractionType == 3)
		{
			strcpy(strTemp, "dp");
		}
		else if(pNet->vEdges[ni]->nInteractionType == 4)
		{
			strcpy(strTemp, "pr");
		}
		else if(pNet->vEdges[ni]->nInteractionType == 5)
		{
			strcpy(strTemp, "rp");
		}
		else
		{
			strcpy(strTemp, "--");
		}

		
		/* multiply */
		if(nNodeID1 != nNodeID2)
		{
			if(nOperationType == 0)
			{
				if( (BMGETAT(pNet->vNodes[nNodeID1]->pBV, 0, 0) == 0) ||
					(BMGETAT(pNet->vNodes[nNodeID2]->pBV, 0, 0) == 0) )
				{
					nHasValue = 0;
					dTemp = 0.0;
				}
				else
				{
					nHasValue = 1;
					dTemp = DMGETAT(pNet->vNodes[nNodeID1]->pDV, 0, 0)*DMGETAT(pNet->vNodes[nNodeID2]->pDV, 0, 0);
				}
			}
			/* add */
			else
			{
				if( (BMGETAT(pNet->vNodes[nNodeID1]->pBV, 0, 0) == 0) ||
					(BMGETAT(pNet->vNodes[nNodeID2]->pBV, 0, 0) == 0) )
				{
					nHasValue = 0;
					dTemp = 0.0;
				}
				else
				{
					nHasValue = 1;
					dTemp = DMGETAT(pNet->vNodes[nNodeID1]->pDV, 0, 0)+DMGETAT(pNet->vNodes[nNodeID2]->pDV, 0, 0);
				}
			}
		}
		else
		{
			if( BMGETAT(pNet->vNodes[nNodeID1]->pBV, 0, 0) == 0)
			{
				nHasValue = 0;
				dTemp = 0.0;
			}
			else
			{
				nHasValue = 1;
				dTemp = fabs(DMGETAT(pNet->vNodes[nNodeID1]->pDV, 0, 0));
			}
		}

		if(nNodeID1 != nNodeID2)
		{
			fprintf(fpOut, "%s\t%s\t%s\t%d\t%f\t%d\t%f\t%d\t%f\n", 
				pNet->vNodes[nNodeID1]->strName->m_pString, strTemp,
				pNet->vNodes[nNodeID2]->strName->m_pString, nHasValue, dTemp,
				BMGETAT(pNet->vNodes[nNodeID1]->pBV, 0, 0),
				DMGETAT(pNet->vNodes[nNodeID1]->pDV, 0, 0),
				BMGETAT(pNet->vNodes[nNodeID2]->pBV, 0, 0),
				DMGETAT(pNet->vNodes[nNodeID2]->pDV, 0, 0));
		}
	}

	fclose(fpOut);

	/* clear memory */
	RefGene_ClearDatabase(&vRefGeneDatabase, nRefGeneNum);
	DestroyIntMatrix(pMap);
	if(vValidationMap != NULL)
	{
		for(ni=0; ni<nRefGeneNum; ni++)
		{
			DeleteString(vValidationMap[ni]);
			vValidationMap[ni] = NULL;
		}
		free(vValidationMap);
	}
	DestroyDoubleMatrix(pValue);
	for(ni=0; ni<nValueNum; ni++)
	{
		DeleteString(vTag[ni]);
		vTag[ni] = NULL;
	}
	free(vTag);
	DestroyByteMatrix(pHasValue);
	DestroyDoubleMatrix(pSigValue);
	/* destroy net */
	if(pNet != NULL)
	{
		NETWORKPHYSICALDESTROY(&pNet);
	}

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Genome_CreateHash_Main                                                 */
/*  Create Hash Table for genome.                                          */
/* ----------------------------------------------------------------------- */ 
int Genome_CreateHash_Main(char strGenomePath[], char strOutPath[], 
						   int nKeyLen)
{
	/* define */
	int nBaseTypeNum = 4;
	char strChr[MED_LINE_LENGTH];
	char strGenomePathC[MED_LINE_LENGTH];
	char strOutPathC[MED_LINE_LENGTH];
	char strFileName[MED_LINE_LENGTH];
	char strLine[MED_LINE_LENGTH];
	FILE *fpIn;
	struct INTMATRIX *pChrLen;
	struct tagString **vChrName = NULL; 
	int nChrNum = 0;
	int ni;

	/* init */
	if(nKeyLen<=0)
	{
		printf("Error: Genome_CreateHash_Main, key length should be >0!\n");
		exit(EXIT_FAILURE);
	}
	strcpy(strGenomePathC, strGenomePath);
	strcpy(strOutPathC, strOutPath);
	AdjustDirectoryPath(strGenomePathC);
	AdjustDirectoryPath(strOutPathC);

	sprintf(strFileName, "%schrlen.txt", strGenomePathC);
	pChrLen = NULL;
	pChrLen = IMLOAD(strFileName);
	if(pChrLen == NULL)
	{
		printf("Error: Genome_CreateHash_Main, cannot load chromosome length!\n");
		exit(EXIT_FAILURE);
	}

	nChrNum = pChrLen->nHeight;
	if(nChrNum <= 0)
	{
		printf("Error: Genome_CreateHash_Main, cannot load chromosome length!\n");
		exit(EXIT_FAILURE);
	}
	vChrName = NULL;
	vChrName = (struct tagString **)calloc(nChrNum, sizeof(struct tagString *));
	if(vChrName == NULL)
	{
		printf("Error: Genome_CreateHash_Main, cannot load chromosome name!\n");
		exit(EXIT_FAILURE);
	}
	
	sprintf(strFileName, "%schrlist.txt", strGenomePathC);
	fpIn = NULL;
	fpIn = fopen(strFileName, "r");
	if(fpIn == NULL)
	{
		printf("Error: Genome_CreateHash_Main, cannot load chromosome name!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, MED_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(ni >= nChrNum)
		{
			printf("Error: Genome_CreateHash_Main, chromosome number not match!\n");
			exit(EXIT_FAILURE);
		}

		StringAddTail(vChrName+ni, strLine);
		ni++;
	}

	fclose(fpIn);

    if(ni != nChrNum)
	{
		printf("Error: Genome_CreateHash_Main, chromosome number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* process one by one */
	for(ni=0; ni<nChrNum; ni++)
	{
		printf("Processing %s ...\n", vChrName[ni]->m_pString);
		Genome_CreateHash_Chr(strGenomePathC, vChrName[ni]->m_pString, 
			strOutPathC, pChrLen->pMatElement[ni], nBaseTypeNum, nKeyLen);
	}

	/* release memory */
	DestroyIntMatrix(pChrLen);
	for(ni=0; ni<nChrNum; ni++)
	{
		DeleteString(vChrName[ni]);
		vChrName[ni] = NULL;
	}
	free(vChrName);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Genome_CreateHash_Chr_v0                                                  */
/*  Create Hash Table for a chromosome.                                    */
/* ----------------------------------------------------------------------- */ 
int Genome_CreateHash_Chr_v0(char strGenomePath[],	char strChr[], char strOutPath[],
				int nChrLen, int nBaseTypeNum, int nKey1Len, int nKey2Len)
{
	/* define */
	int **vGenomeHash = NULL;
	int nDim1,nDim2;
	int nDS1,nDS2;
	int ni,nj,nk,nz;
	int nM1,nM2;
	int nSeqCount = 0;
	struct FLEXSEQMOTIF *pSeqMtf;
	int nActualStart,nActualEnd,nTo;
	int nMaxHitNum = 0;

	int nWID1[2],nWID2[2];
	int nWBad1[2],nWBad2[2];
	int nSwitch = 0;
	unsigned char *pB1,*pB2;

	/* hash index */
	int nTemp;
	int nTotalWord = 0;
	double dTotalSearch = 0.0;
	int nTotalIndex = 0;
	FILE *fpOut;
	FILE *fpIn;
	char strFileName[MED_LINE_LENGTH];
	int nHashID;
	int nLine1,nLine2;

	/* step1: initialize */
	nDim1 = (int)pow((double)nBaseTypeNum, (double)nKey1Len);
	nDim2 = (int)pow((double)nBaseTypeNum, (double)nKey2Len);
	if(nKey1Len == 0)
		nDS1 = 0;
	else
		nDS1 = (int)pow((double)nBaseTypeNum, (double)(nKey1Len-1));
	if(nKey2Len == 0)
		nDS2 = 0;
	else
		nDS2 = (int)pow((double)nBaseTypeNum, (double)(nKey2Len-1));
	vGenomeHash = (int **)calloc(nDim1, sizeof(int *));
	if(vGenomeHash == NULL)
	{
		printf("Error: Genome_CreateHash_Chr, cannot create hash table!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nDim1; ni++)
	{
		vGenomeHash[ni] = (int *)calloc(nDim2, sizeof(int));
		if(vGenomeHash[ni] == NULL)
		{
			printf("Error: Genome_CreateHash_Chr, cannot create hash table!\n");
			exit(EXIT_FAILURE);
		}
	}

	/* step2: load sequence and compute hash table size */
	nM1 = 0;
	while(nM1 < nChrLen)
	{
		nM2 = nM1+GENOME_CONTIG_LEN-1;
		if(nM2 >= nChrLen)
			nM2 = nChrLen-1;

		nActualStart = nM1;
		nActualEnd = nM2+2*(nKey1Len+nKey2Len-1);
		if(nActualEnd >= nChrLen)
		{
			nActualEnd = nChrLen-1;
		}
		if( (nActualEnd-2*(nKey1Len+nKey2Len-1))< nM1)
		{
			break;
		}
		nTo = nActualEnd-2*(nKey1Len+nKey2Len-1);

		/* #################################### */ 
		/* load sequence and score              */
		/* #################################### */
		pSeqMtf = NULL;
		pSeqMtf = FLEXSEQMTFLOADSEQ_FOR_HASHGENOME(nSeqCount,
				strGenomePath, strChr, nActualStart, nActualEnd, 
				0, "NULL");
		if(pSeqMtf == NULL)
		{
			break;
		}

			
		/* #################################### */
		/* hash sequence                        */
		/* #################################### */
		for(nj=0; nj<2; nj++)
		{
			nWID1[nj] = 0;
			nWID2[nj] = 0;
			nWBad1[nj] = 0;
			nWBad2[nj] = 0;
		}
		nSwitch = 0;
		pB1 = pSeqMtf->vSeq[0]->pMatElement;
		pB2 = pB1+2*nKey1Len;
		
		for(ni=0; ni<nKey1Len; ni++)
		{
			nz = 2*ni;
			if(pB1[nz] < nBaseTypeNum)
			{
				nWID1[nSwitch] = nBaseTypeNum*nWID1[nSwitch]+pB1[nz];
			}
			else
			{
				nWID1[nSwitch] = nBaseTypeNum*nWID1[nSwitch];
				nWBad1[nSwitch] += 1;
			}
		}

		for(ni=0; ni<nKey2Len; ni++)
		{
			nz = 2*ni;
			if(pB2[nz] < nBaseTypeNum)
			{
				nWID2[nSwitch] = nBaseTypeNum*nWID2[nSwitch]+pB2[nz];
			}
			else
			{
				nWID2[nSwitch] = nBaseTypeNum*nWID2[nSwitch];
				nWBad2[nSwitch] += 1;
			}
		}

		nk = nActualStart;

		if( (nWBad1[nSwitch] == 0) && (nWBad2[nSwitch] == 0) )
		{
			/* Add to Hash */
			nz = nWID1[nSwitch];
			*(vGenomeHash[nz]+nWID2[nSwitch]) += 1;
			if(*(vGenomeHash[nz]+nWID2[nSwitch]) > nMaxHitNum)
				nMaxHitNum = *(vGenomeHash[nz]+nWID2[nSwitch]);
		}

		nSwitch = 1-nSwitch;
		
		if(nk < nTo)
		{
			pB1++;
			pB2++;
		
			for(ni=0; ni<nKey1Len; ni++)
			{
				nz = 2*ni;
				if(pB1[nz] < nBaseTypeNum)
				{
					nWID1[nSwitch] = nBaseTypeNum*nWID1[nSwitch]+pB1[nz];
				}
				else
				{
					nWID1[nSwitch] = nBaseTypeNum*nWID1[nSwitch];
					nWBad1[nSwitch] += 1;
				}
			}

			for(ni=0; ni<nKey2Len; ni++)
			{
				nz = 2*ni;
				if(pB2[nz] < nBaseTypeNum)
				{
					nWID2[nSwitch] = nBaseTypeNum*nWID2[nSwitch]+pB2[nz];
				}
				else
				{
					nWID2[nSwitch] = nBaseTypeNum*nWID2[nSwitch];
					nWBad2[nSwitch] += 1;
				}
			}

			if( (nWBad1[nSwitch] == 0) && (nWBad2[nSwitch] == 0) )
			{
				/* Add to Hash */
				nz = nWID1[nSwitch];
				*(vGenomeHash[nz]+nWID2[nSwitch]) += 1;
				if(*(vGenomeHash[nz]+nWID2[nSwitch]) > nMaxHitNum)
					nMaxHitNum = *(vGenomeHash[nz]+nWID2[nSwitch]);
			}

			nk++;
			nSwitch = 1-nSwitch;

		}

		while(nk < nTo)
		{
			pB1++;
			pB2++;
			
			/* update key value */
			if( (nKey1Len > 0) && (nKey2Len > 0) )
			{
				if( *(pB1-2) < nBaseTypeNum)
				{
					nWID1[nSwitch] -= nDS1*(*(pB1-2));
				}
				else
				{
					nWBad1[nSwitch] -= 1;
				}
			

				if(*(pB2-2) < nBaseTypeNum)
				{
					nWID2[nSwitch] -= nDS2*(*(pB2-2));
					nWID1[nSwitch] = nBaseTypeNum*nWID1[nSwitch]+(*(pB2-2));
				}
				else
				{
                    nWBad2[nSwitch] -= 1;
					nWID1[nSwitch] = nBaseTypeNum*nWID1[nSwitch];
					nWBad1[nSwitch] += 1;
				}

				nz = 2*nKey2Len-2;
				if(pB2[nz] < nBaseTypeNum)
				{
					nWID2[nSwitch] = nBaseTypeNum*nWID2[nSwitch]+pB2[nz];
				}
				else
				{
					nWID2[nSwitch] = nBaseTypeNum*nWID2[nSwitch];
					nWBad2[nSwitch] += 1;
				}
			}
			else if(nKey1Len > 0)
			{
				printf("Error: Genome_CreateHash_Chr, key2len must be > 0!\n");
				exit(EXIT_FAILURE);
			}
			else if(nKey2Len > 0)
			{
				if( *(pB2-2) < nBaseTypeNum)
				{
					nWID2[nSwitch] -= nDS2*(*(pB2-2));
				}
				else
				{
					nWBad2[nSwitch] -= 1;
				}

				nz = 2*nKey2Len-2;
				if(pB2[nz] < nBaseTypeNum)
				{
					nWID2[nSwitch] = nBaseTypeNum*nWID2[nSwitch]+pB2[nz];
				}
				else
				{
					nWID2[nSwitch] = nBaseTypeNum*nWID2[nSwitch];
					nWBad2[nSwitch] += 1;
				}
			}
			else
			{
				printf("Error: Genome_CreateHash_Chr, cannot create hash key!\n");
				exit(EXIT_FAILURE);
			}

			if( (nWBad1[nSwitch] == 0) && (nWBad2[nSwitch] == 0) )
			{
				/* Add to Hash */
				nz = nWID1[nSwitch];
				*(vGenomeHash[nz]+nWID2[nSwitch]) +=1;
				if(*(vGenomeHash[nz]+nWID2[nSwitch]) > nMaxHitNum)
					nMaxHitNum = *(vGenomeHash[nz]+nWID2[nSwitch]);
			}

			nk++;
			nSwitch = 1-nSwitch;
		}

		/* #################################### */
		/* release memory                       */
		/* #################################### */
		FLEXSEQMTFDESTROY(pSeqMtf);
		
		nM1 = nM2+1;
		nSeqCount++;
	}

	/* step3: create hash index */
	sprintf(strFileName, "%s%s.hashidx", strOutPath, strChr);
	fpOut = NULL;
	fpOut = fopen(strFileName, "wb");
	if(fpOut == NULL)
	{
		printf("Error: Genome_CreateHash_Chr, cannot save hash index!\n");
		exit(EXIT_FAILURE);
	}
	nTotalWord = 0;
	dTotalSearch = 0.0;
	nTotalIndex = 0;
	fwrite(&nTotalWord, sizeof(int), 1, fpOut);
	for(ni=0; ni<nDim1; ni++)
	{
		for(nj=0; nj<nDim2; nj++)
		{
			nTemp = *(vGenomeHash[ni]+nj);
			if(nTemp > 0)
			{
				nTotalWord += nTemp;
				dTotalSearch += nTemp*nTemp;
				nTotalIndex += 1;
			}
			fwrite(&nTotalWord, sizeof(int), 1, fpOut);
		}
	}
	fclose(fpOut);

	dTotalSearch /= (double)nTotalWord;
	
	/* first release of memory */
	for(ni=0; ni<nDim1; ni++)
	{
		free(vGenomeHash[ni]);
	}
	free(vGenomeHash);


	/* step4: save hash map */
	vGenomeHash = NULL;
	vGenomeHash = (int **)calloc(nDim1, sizeof(int *));
	if(vGenomeHash == NULL)
	{
		printf("Error: Genome_CreateHash_Chr, cannot create hash table!\n");
		exit(EXIT_FAILURE);
	}
	
	for(ni=0; ni<nDim1; ni++)
	{
		vGenomeHash[ni] = (int *)calloc(nDim2, sizeof(int));
		if(vGenomeHash[ni] == NULL)
		{
			printf("Error: Genome_CreateHash_Chr, cannot create hash table!\n");
			exit(EXIT_FAILURE);
		}
	}
    
	nMaxHitNum = 0;
	fpIn = NULL;
	fpIn = fopen(strFileName, "rb");
	if(fpIn == NULL)
	{
		printf("Error: Genome_CreateHash_Chr, cannot save hash index!\n");
		exit(EXIT_FAILURE);
	}

	sprintf(strFileName, "%s%s.hashmap", strOutPath, strChr);
	fpOut = NULL;
	fpOut = fopen(strFileName, "wb");
	if(fpOut == NULL)
	{
		printf("Error: Genome_CreateHash_Chr, cannot save hash index!\n");
		exit(EXIT_FAILURE);
	}
	
	nz = 0;
	for(ni=0; ni<nTotalWord; ni++)
	{
		fwrite(&nz, sizeof(int), 1, fpOut);
	}

	fclose(fpOut);

	fpOut = NULL;
	fpOut = fopen(strFileName, "r+b");
	if(fpOut == NULL)
	{
		printf("Error: Genome_CreateHash_Chr, cannot save hash index!\n");
		exit(EXIT_FAILURE);
	}

	nM1 = 0;
	while(nM1 < nChrLen)
	{
		printf("%d...\n", nM1);
		nM2 = nM1+GENOME_CONTIG_LEN-1;
		if(nM2 >= nChrLen)
			nM2 = nChrLen-1;

		nActualStart = nM1;
		nActualEnd = nM2+2*(nKey1Len+nKey2Len-1);
		if(nActualEnd >= nChrLen)
		{
			nActualEnd = nChrLen-1;
		}
		if( (nActualEnd-2*(nKey1Len+nKey2Len-1))< nM1)
		{
			break;
		}
		nTo = nActualEnd-2*(nKey1Len+nKey2Len-1);

		/* #################################### */ 
		/* load sequence and score              */
		/* #################################### */
		pSeqMtf = NULL;
		pSeqMtf = FLEXSEQMTFLOADSEQ_FOR_HASHGENOME(nSeqCount,
				strGenomePath, strChr, nActualStart, nActualEnd, 
				0, "NULL");
		if(pSeqMtf == NULL)
		{
			break;
		}

			
		/* #################################### */
		/* hash sequence                        */
		/* #################################### */
		for(nj=0; nj<2; nj++)
		{
			nWID1[nj] = 0;
			nWID2[nj] = 0;
			nWBad1[nj] = 0;
			nWBad2[nj] = 0;
		}
		nSwitch = 0;
		pB1 = pSeqMtf->vSeq[0]->pMatElement;
		pB2 = pB1+2*nKey1Len;
		
		for(ni=0; ni<nKey1Len; ni++)
		{
			nz = 2*ni;
			if(pB1[nz] < nBaseTypeNum)
			{
				nWID1[nSwitch] = nBaseTypeNum*nWID1[nSwitch]+pB1[nz];
			}
			else
			{
				nWID1[nSwitch] = nBaseTypeNum*nWID1[nSwitch];
				nWBad1[nSwitch] += 1;
			}
		}

		for(ni=0; ni<nKey2Len; ni++)
		{
			nz = 2*ni;
			if(pB2[nz] < nBaseTypeNum)
			{
				nWID2[nSwitch] = nBaseTypeNum*nWID2[nSwitch]+pB2[nz];
			}
			else
			{
				nWID2[nSwitch] = nBaseTypeNum*nWID2[nSwitch];
				nWBad2[nSwitch] += 1;
			}
		}

		nk = nActualStart;

		if( (nWBad1[nSwitch] == 0) && (nWBad2[nSwitch] == 0) )
		{
			/* Add to Hash */
			nz = nWID1[nSwitch];
			nHashID = nz*nDim2+nWID2[nSwitch];
			
            if( fseek( fpIn, nHashID*sizeof(int), SEEK_SET ) != 0 )
			{
				printf("Error: Genome_CreateHash_Chr, cannot locate the required hash index!\n");
				exit(EXIT_FAILURE);
			}

			fread(&nLine1, sizeof(int), 1, fpIn);
			fread(&nLine2, sizeof(int), 1, fpIn);

			if(*(vGenomeHash[nz]+nWID2[nSwitch]) >= (nLine2-nLine1))
			{
				printf("Error: Genome_CreateHash_Chr, hash index not match!\n");
				exit(EXIT_FAILURE);
			}

			nLine1 += *(vGenomeHash[nz]+nWID2[nSwitch]);
			if( fseek( fpOut, nLine1*sizeof(int), SEEK_SET ) != 0 )
			{
				printf("Error: Genome_CreateHash_Chr, cannot locate the required hash map!\n");
				exit(EXIT_FAILURE);
			}

			fwrite(&nk, sizeof(int), 1, fpOut);
		
			*(vGenomeHash[nz]+nWID2[nSwitch]) += 1;
			if(*(vGenomeHash[nz]+nWID2[nSwitch]) > nMaxHitNum)
				nMaxHitNum = *(vGenomeHash[nz]+nWID2[nSwitch]);
		}

		nSwitch = 1-nSwitch;
		
		if(nk < nTo)
		{
			pB1++;
			pB2++;
		
			for(ni=0; ni<nKey1Len; ni++)
			{
				nz = 2*ni;
				if(pB1[nz] < nBaseTypeNum)
				{
					nWID1[nSwitch] = nBaseTypeNum*nWID1[nSwitch]+pB1[nz];
				}
				else
				{
					nWID1[nSwitch] = nBaseTypeNum*nWID1[nSwitch];
					nWBad1[nSwitch] += 1;
				}
			}

			for(ni=0; ni<nKey2Len; ni++)
			{
				nz = 2*ni;
				if(pB2[nz] < nBaseTypeNum)
				{
					nWID2[nSwitch] = nBaseTypeNum*nWID2[nSwitch]+pB2[nz];
				}
				else
				{
					nWID2[nSwitch] = nBaseTypeNum*nWID2[nSwitch];
					nWBad2[nSwitch] += 1;
				}
			}

			nk++;
			if( (nWBad1[nSwitch] == 0) && (nWBad2[nSwitch] == 0) )
			{
				/* Add to Hash */
				nz = nWID1[nSwitch];
				nHashID = nz*nDim2+nWID2[nSwitch];
				
				if( fseek( fpIn, nHashID*sizeof(int), SEEK_SET ) != 0 )
				{
					printf("Error: Genome_CreateHash_Chr, cannot locate the required hash index!\n");
					exit(EXIT_FAILURE);
				}

				fread(&nLine1, sizeof(int), 1, fpIn);
				fread(&nLine2, sizeof(int), 1, fpIn);

				if(*(vGenomeHash[nz]+nWID2[nSwitch]) >= (nLine2-nLine1))
				{
					printf("Error: Genome_CreateHash_Chr, hash index not match!\n");
					exit(EXIT_FAILURE);
				}

				nLine1 += *(vGenomeHash[nz]+nWID2[nSwitch]);
				if( fseek( fpOut, nLine1*sizeof(int), SEEK_SET ) != 0 )
				{
					printf("Error: Genome_CreateHash_Chr, cannot locate the required hash map!\n");
					exit(EXIT_FAILURE);
				}

				fwrite(&nk, sizeof(int), 1, fpOut);
			
				*(vGenomeHash[nz]+nWID2[nSwitch]) += 1;
				if(*(vGenomeHash[nz]+nWID2[nSwitch]) > nMaxHitNum)
					nMaxHitNum = *(vGenomeHash[nz]+nWID2[nSwitch]);
			}

			nSwitch = 1-nSwitch;
		}

		while(nk < nTo)
		{
			pB1++;
			pB2++;
			
			/* update key value */
			if( (nKey1Len > 0) && (nKey2Len > 0) )
			{
				if( *(pB1-2) < nBaseTypeNum)
				{
					nWID1[nSwitch] -= nDS1*(*(pB1-2));
				}
				else
				{
					nWBad1[nSwitch] -= 1;
				}
			

				if(*(pB2-2) < nBaseTypeNum)
				{
					nWID2[nSwitch] -= nDS2*(*(pB2-2));
					nWID1[nSwitch] = nBaseTypeNum*nWID1[nSwitch]+(*(pB2-2));
				}
				else
				{
                    nWBad2[nSwitch] -= 1;
					nWID1[nSwitch] = nBaseTypeNum*nWID1[nSwitch];
					nWBad1[nSwitch] += 1;
				}

				nz = 2*nKey2Len-2;
				if(pB2[nz] < nBaseTypeNum)
				{
					nWID2[nSwitch] = nBaseTypeNum*nWID2[nSwitch]+pB2[nz];
				}
				else
				{
					nWID2[nSwitch] = nBaseTypeNum*nWID2[nSwitch];
					nWBad2[nSwitch] += 1;
				}
			}
			else if(nKey1Len > 0)
			{
				printf("Error: Genome_CreateHash_Chr, key2len must be > 0!\n");
				exit(EXIT_FAILURE);
			}
			else if(nKey2Len > 0)
			{
				if( *(pB2-2) < nBaseTypeNum)
				{
					nWID2[nSwitch] -= nDS2*(*(pB2-2));
				}
				else
				{
					nWBad2[nSwitch] -= 1;
				}

				nz = 2*nKey2Len-2;
				if(pB2[nz] < nBaseTypeNum)
				{
					nWID2[nSwitch] = nBaseTypeNum*nWID2[nSwitch]+pB2[nz];
				}
				else
				{
					nWID2[nSwitch] = nBaseTypeNum*nWID2[nSwitch];
					nWBad2[nSwitch] += 1;
				}
			}
			else
			{
				printf("Error: Genome_CreateHash_Chr, cannot create hash key!\n");
				exit(EXIT_FAILURE);
			}

			nk++;
			if( (nWBad1[nSwitch] == 0) && (nWBad2[nSwitch] == 0) )
			{
				/* Add to Hash */
				nz = nWID1[nSwitch];
				nHashID = nz*nDim2+nWID2[nSwitch];
				
				if( fseek( fpIn, nHashID*sizeof(int), SEEK_SET ) != 0 )
				{
					printf("Error: Genome_CreateHash_Chr, cannot locate the required hash index!\n");
					exit(EXIT_FAILURE);
				}

				fread(&nLine1, sizeof(int), 1, fpIn);
				fread(&nLine2, sizeof(int), 1, fpIn);

				if(*(vGenomeHash[nz]+nWID2[nSwitch]) >= (nLine2-nLine1))
				{
					printf("Error: Genome_CreateHash_Chr, hash index not match!\n");
					exit(EXIT_FAILURE);
				}

				nLine1 += *(vGenomeHash[nz]+nWID2[nSwitch]);
				if( fseek( fpOut, nLine1*sizeof(int), SEEK_SET ) != 0 )
				{
					printf("Error: Genome_CreateHash_Chr, cannot locate the required hash map!\n");
					exit(EXIT_FAILURE);
				}

				fwrite(&nk, sizeof(int), 1, fpOut);
			
				*(vGenomeHash[nz]+nWID2[nSwitch]) += 1;
				if(*(vGenomeHash[nz]+nWID2[nSwitch]) > nMaxHitNum)
					nMaxHitNum = *(vGenomeHash[nz]+nWID2[nSwitch]);
			}
			nSwitch = 1-nSwitch;
		}

		/* #################################### */
		/* release memory                       */
		/* #################################### */
		FLEXSEQMTFDESTROY(pSeqMtf);
		
		nM1 = nM2+1;
		nSeqCount++;
	}

	fclose(fpIn);
	fclose(fpOut);

	/* final release of memory */
	for(ni=0; ni<nDim1; ni++)
	{
		free(vGenomeHash[ni]);
	}
	free(vGenomeHash);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*  Genome_CreateHash_Chr                                                  */
/*  Create Hash Table for a chromosome.                                    */
/* ----------------------------------------------------------------------- */ 
int Genome_CreateHash_Chr(char strGenomePath[],	char strChr[], char strOutPath[],
				int nChrLen, int nBaseTypeNum, int nKeyLen)
{
	/* define */
	int *vGenomeHash = NULL;
	int nDim;
	int nWID[2];
	int nWBad[2];
	int nSwitch = 0;
	unsigned char *pB;

	int nDS;
	int ni,nj,nk,nz;
	int nM1,nM2;
	int nSeqCount = 0;
	struct FLEXSEQMOTIF *pSeqMtf;
	int nActualStart,nActualEnd,nTo;
	int nMaxHitNum = 0;

	/* hash index */
	int *vGenomeHashIndex = NULL;
	int *vGenomeHashTable = NULL;

	int nTemp;
	int nTotalWord = 0;
	double dTotalSearch = 0.0;
	int nTotalIndex = 0;
	FILE *fpOut;
	FILE *fpIn;
	char strFileName[MED_LINE_LENGTH];
	int nHashID;
	int nLine1,nLine2;

	/* step1: initialize */
	nDS = (int)pow((double)nBaseTypeNum, (double)(nKeyLen-1));
	nDim = (int)pow((double)nBaseTypeNum, (double)nKeyLen);
	vGenomeHash = (int *)calloc(nDim, sizeof(int));
	if(vGenomeHash == NULL)
	{
		printf("Error: Genome_CreateHash_Chr, cannot create hash table!\n");
		exit(EXIT_FAILURE);
	}

	/* step2: load sequence and compute hash table size */
	nM1 = 0;
	while(nM1 < nChrLen)
	{
		nM2 = nM1+GENOME_CONTIG_LEN-1;
		if(nM2 >= nChrLen)
			nM2 = nChrLen-1;

		nActualStart = nM1;
		nActualEnd = nM2+2*(nKeyLen-1);
		if(nActualEnd >= nChrLen)
		{
			nActualEnd = nChrLen-1;
		}
		if( (nActualEnd-2*(nKeyLen-1))< nM1)
		{
			break;
		}

		nTo = nActualEnd-2*(nKeyLen-1);

		/* #################################### */ 
		/* load sequence and score              */
		/* #################################### */
		pSeqMtf = NULL;
		pSeqMtf = FLEXSEQMTFLOADSEQ_FOR_HASHGENOME(nSeqCount,
				strGenomePath, strChr, nActualStart, nActualEnd, 
				0, "NULL");
		if(pSeqMtf == NULL)
		{
			break;
		}

			
		/* #################################### */
		/* hash sequence                        */
		/* #################################### */
		for(nj=0; nj<2; nj++)
		{
			nWID[nj] = 0;
			nWBad[nj] = 0;
		}
		nSwitch = 0;
		pB = pSeqMtf->vSeq[0]->pMatElement;
		
		for(ni=0; ni<nKeyLen; ni++)
		{
			nz = 2*ni;
			if(pB[nz] < nBaseTypeNum)
			{
				nWID[nSwitch] = nBaseTypeNum*nWID[nSwitch]+pB[nz];
			}
			else
			{
				nWID[nSwitch] = nBaseTypeNum*nWID[nSwitch];
				nWBad[nSwitch] += 1;
			}
		}

		nk = nActualStart;

		if(nWBad[nSwitch] == 0)
		{
			/* Add to Hash */
			nz = nWID[nSwitch];
			vGenomeHash[nz] += 1;
			if(vGenomeHash[nz] > nMaxHitNum)
				nMaxHitNum = vGenomeHash[nz];
		}

		nSwitch = 1-nSwitch;
		
		if(nk < nTo)
		{
			pB++;
			
			for(ni=0; ni<nKeyLen; ni++)
			{
				nz = 2*ni;
				if(pB[nz] < nBaseTypeNum)
				{
					nWID[nSwitch] = nBaseTypeNum*nWID[nSwitch]+pB[nz];
				}
				else
				{
					nWID[nSwitch] = nBaseTypeNum*nWID[nSwitch];
					nWBad[nSwitch] += 1;
				}
			}

			if(nWBad[nSwitch] == 0)
			{
				/* Add to Hash */
				nz = nWID[nSwitch];
				vGenomeHash[nz] += 1;
				if(vGenomeHash[nz] > nMaxHitNum)
					nMaxHitNum = vGenomeHash[nz];
			}

			nk++;
			nSwitch = 1-nSwitch;
		}

		while(nk < nTo)
		{
			pB++;
			
			/* update key value */
			if( *(pB-2) < nBaseTypeNum)
			{
				nWID[nSwitch] -= nDS*(*(pB-2));
			}
			else
			{
				nWBad[nSwitch] -= 1;
			}
			
			nz = 2*nKeyLen-2;
			if(pB[nz] < nBaseTypeNum)
			{
				nWID[nSwitch] = nBaseTypeNum*nWID[nSwitch]+pB[nz];
			}
			else
			{
				nWID[nSwitch] = nBaseTypeNum*nWID[nSwitch];
				nWBad[nSwitch] += 1;
			}
			

			if(nWBad[nSwitch] == 0)
			{
				/* Add to Hash */
				nz = nWID[nSwitch];
				vGenomeHash[nz] +=1;
				if(vGenomeHash[nz] > nMaxHitNum)
					nMaxHitNum = vGenomeHash[nz];
			}

			nk++;
			nSwitch = 1-nSwitch;
		}

		/* #################################### */
		/* release memory                       */
		/* #################################### */
		FLEXSEQMTFDESTROY(pSeqMtf);
		
		nM1 = nM2+1;
		nSeqCount++;
	}

	/* step3: create hash index */
	vGenomeHashIndex = (int *)calloc((nDim+1), sizeof(int));
	if(vGenomeHashIndex == NULL)
	{
		printf("Error: Genome_CreateHash_Chr, cannot create hash index!\n");
		exit(EXIT_FAILURE);
	}

	nTotalWord = 0;
	dTotalSearch = 0.0;
	nTotalIndex = 0;
	vGenomeHashIndex[0] = nTotalWord;
	for(ni=0; ni<nDim; ni++)
	{
		nTotalWord += vGenomeHash[ni];
		vGenomeHashIndex[ni+1] = nTotalWord;
		if(vGenomeHash[ni] > 0)
		{
			dTotalSearch += vGenomeHash[ni]*vGenomeHash[ni];
			nTotalIndex += 1;
		}
	}
	dTotalSearch /= (double)nTotalWord;
	
	sprintf(strFileName, "%s%s.hashidx", strOutPath, strChr);
	fpOut = NULL;
	fpOut = fopen(strFileName, "wb");
	if(fpOut == NULL)
	{
		printf("Error: Genome_CreateHash_Chr, cannot save hash index!\n");
		exit(EXIT_FAILURE);
	}
	fwrite(vGenomeHashIndex, sizeof(int), (nDim+1), fpOut);
	fclose(fpOut);

	/* first release of memory */
	free(vGenomeHash);

	/* step4: save hash map */
	vGenomeHash = NULL;
	vGenomeHash = (int *)calloc(nDim, sizeof(int));
	if(vGenomeHash == NULL)
	{
		printf("Error: Genome_CreateHash_Chr, cannot create hash table!\n");
		exit(EXIT_FAILURE);
	}

	vGenomeHashTable = NULL;
	vGenomeHashTable = (int *)calloc(nTotalWord, sizeof(int));
	if(vGenomeHashTable == NULL)
	{
		printf("Error: Genome_CreateHash_Chr, cannot create hash map!\n");
		exit(EXIT_FAILURE);
	}
	
	nMaxHitNum = 0;
	
	nM1 = 0;
	while(nM1 < nChrLen)
	{
		/* printf("%d...\n", nM1); */
		nM2 = nM1+GENOME_CONTIG_LEN-1;
		if(nM2 >= nChrLen)
			nM2 = nChrLen-1;

		nActualStart = nM1;
		nActualEnd = nM2+2*(nKeyLen-1);
		if(nActualEnd >= nChrLen)
		{
			nActualEnd = nChrLen-1;
		}
		if( (nActualEnd-2*(nKeyLen-1))< nM1)
		{
			break;
		}
		nTo = nActualEnd-2*(nKeyLen-1);

		/* #################################### */ 
		/* load sequence and score              */
		/* #################################### */
		pSeqMtf = NULL;
		pSeqMtf = FLEXSEQMTFLOADSEQ_FOR_HASHGENOME(nSeqCount,
				strGenomePath, strChr, nActualStart, nActualEnd, 
				0, "NULL");
		if(pSeqMtf == NULL)
		{
			break;
		}

			
		/* #################################### */
		/* hash sequence                        */
		/* #################################### */
		for(nj=0; nj<2; nj++)
		{
			nWID[nj] = 0;
			nWBad[nj] = 0;
		}
		nSwitch = 0;
		pB = pSeqMtf->vSeq[0]->pMatElement;
		
		for(ni=0; ni<nKeyLen; ni++)
		{
			nz = 2*ni;
			if(pB[nz] < nBaseTypeNum)
			{
				nWID[nSwitch] = nBaseTypeNum*nWID[nSwitch]+pB[nz];
			}
			else
			{
				nWID[nSwitch] = nBaseTypeNum*nWID[nSwitch];
				nWBad[nSwitch] += 1;
			}
		}

		nk = nActualStart;

		if( nWBad[nSwitch] == 0 )
		{
			/* Add to Hash */
			nHashID = nWID[nSwitch];
			nLine1 = vGenomeHashIndex[nHashID];
			nLine2 = vGenomeHashIndex[nHashID+1];

			if( vGenomeHash[nHashID] >= (nLine2-nLine1))
			{
				printf("Error: Genome_CreateHash_Chr, hash index not match!\n");
				exit(EXIT_FAILURE);
			}

			nLine1 += vGenomeHash[nHashID];
			vGenomeHashTable[nLine1] = nk;
			vGenomeHash[nHashID] += 1;

			if(vGenomeHash[nHashID] > nMaxHitNum)
				nMaxHitNum = vGenomeHash[nHashID];
		}

		nSwitch = 1-nSwitch;
		
		if(nk < nTo)
		{
			pB++;
			
			for(ni=0; ni<nKeyLen; ni++)
			{
				nz = 2*ni;
				if(pB[nz] < nBaseTypeNum)
				{
					nWID[nSwitch] = nBaseTypeNum*nWID[nSwitch]+pB[nz];
				}
				else
				{
					nWID[nSwitch] = nBaseTypeNum*nWID[nSwitch];
					nWBad[nSwitch] += 1;
				}
			}

			nk++;
			if( nWBad[nSwitch] == 0 )
			{
				/* Add to Hash */
				nHashID = nWID[nSwitch];
				nLine1 = vGenomeHashIndex[nHashID];
				nLine2 = vGenomeHashIndex[nHashID+1];

				if( vGenomeHash[nHashID] >= (nLine2-nLine1))
				{
					printf("Error: Genome_CreateHash_Chr, hash index not match!\n");
					exit(EXIT_FAILURE);
				}

				nLine1 += vGenomeHash[nHashID];
				vGenomeHashTable[nLine1] = nk;
				vGenomeHash[nHashID] += 1;

				if(vGenomeHash[nHashID] > nMaxHitNum)
					nMaxHitNum = vGenomeHash[nHashID];
			}

			nSwitch = 1-nSwitch;
		}

		while(nk < nTo)
		{
			pB++;
			
			/* update key value */
			if( *(pB-2) < nBaseTypeNum)
			{
				nWID[nSwitch] -= nDS*(*(pB-2));
			}
			else
			{
				nWBad[nSwitch] -= 1;
			}
		
			nz = 2*nKeyLen-2;
			if(pB[nz] < nBaseTypeNum)
			{
				nWID[nSwitch] = nBaseTypeNum*nWID[nSwitch]+pB[nz];
			}
			else
			{
				nWID[nSwitch] = nBaseTypeNum*nWID[nSwitch];
				nWBad[nSwitch] += 1;
			}
			
			nk++;

			if( nWBad[nSwitch] == 0 )
			{
				/* Add to Hash */
				nHashID = nWID[nSwitch];
				nLine1 = vGenomeHashIndex[nHashID];
				nLine2 = vGenomeHashIndex[nHashID+1];

				if( vGenomeHash[nHashID] >= (nLine2-nLine1))
				{
					printf("Error: Genome_CreateHash_Chr, hash index not match!\n");
					exit(EXIT_FAILURE);
				}

				nLine1 += vGenomeHash[nHashID];
				vGenomeHashTable[nLine1] = nk;
				vGenomeHash[nHashID] += 1;

				if(vGenomeHash[nHashID] > nMaxHitNum)
					nMaxHitNum = vGenomeHash[nHashID];
			}
			nSwitch = 1-nSwitch;
		}

		/* #################################### */
		/* release memory                       */
		/* #################################### */
		FLEXSEQMTFDESTROY(pSeqMtf);
		
		nM1 = nM2+1;
		nSeqCount++;
	}

	sprintf(strFileName, "%s%s.hashmap", strOutPath, strChr);
	fpOut = NULL;
	fpOut = fopen(strFileName, "wb");
	if(fpOut == NULL)
	{
		printf("Error: Genome_CreateHash_Chr, cannot save hash map!\n");
		exit(EXIT_FAILURE);
	}
	fwrite(vGenomeHashTable, sizeof(int), nTotalWord, fpOut);
	fclose(fpOut);

	/* final release of memory */
	free(vGenomeHash);
	free(vGenomeHashIndex);
	free(vGenomeHashTable);

	/* return */
	return PROC_SUCCESS;
}



/* ----------------------------------------------------------------------- */ 
/*  Genome_RegionExtend_Main: extend genomic regions                       */
/* ----------------------------------------------------------------------- */ 
int Genome_RegionExtend_Main(char strInputPath[], char strGenomePath[], 
			char strOutputPath[], char strSpecies[], int nL, int nR, 
			int nCN, int nA, int nUseStrand)
{
	/* define */
	FILE *fpIn;
	FILE *fpOut;
	struct INTMATRIX *pChrLen;
	int nChr,nStart,nEnd;
	char chStrand;
	char strChr[LINE_LENGTH];
	char strAlias[LINE_LENGTH];

	char strLine[LONG_LINE_LENGTH];
	char strFileName[MED_LINE_LENGTH];
	int ni;

	/* init */
	StrMakeUpper(strSpecies);
	AdjustDirectoryPath(strGenomePath);
	
	/* get chromosome length */
	sprintf(strFileName, "%schrlen.txt", strGenomePath);
	pChrLen = NULL;
	pChrLen = IMLOAD(strFileName);
	if(pChrLen == NULL)
	{
		printf("Error: Genome_RegionExtend_Main, cannot load chromosome length!\n");
		exit(EXIT_FAILURE);
	}

	/* load file */
	fpIn = NULL;
	fpIn = fopen(strInputPath, "r");
	if(fpIn == NULL)
	{
		printf("Error: Genome_RegionExtend_Main, cannot open input file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutputPath, "w");
	if(fpOut == NULL)
	{
		printf("Error: Genome_RegionExtend_Main, cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	fprintf(fpOut, "#seq_id\tchromosome\tstart\tend\tstrand\n");

	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		sscanf(strLine, "%s %s %d %d %c", strAlias, strChr, &nStart, &nEnd, &chStrand);

		nChr = Genome_ChromosomeName_To_Index(strChr, strSpecies);
		if(nChr <= 0)
			continue;

		if(nUseStrand == 1)
		{
			if(chStrand == '-')
			{
				nStart -= nR;
				nEnd += nL;
			}
			else
			{
				nStart -= nL;
				nEnd += nR;
			}
		}
		else
		{
			nStart -= nL;
			nEnd += nR;
		}

		if(nStart < 0)
			nStart = 0;
		
		if(nEnd >= pChrLen->pMatElement[nChr-1])
			nEnd = pChrLen->pMatElement[nChr-1]-1;

		if(nA == 0)
		{
			if(nCN == 1)
			{
				fprintf(fpOut, "%s\t%d\t%d\t%d\t%c\n", strAlias, nChr, nStart, nEnd, chStrand);
			}
			else
			{
				fprintf(fpOut, "%s\t%s\t%d\t%d\t%c\n", strAlias, strChr, nStart, nEnd, chStrand);
			}
		}
		else
		{
			if(nCN == 1)
			{
				fprintf(fpOut, "%d\t%d\t%d\t%d\t%c\n", ni, nChr, nStart, nEnd, chStrand);
			}
			else
			{
				fprintf(fpOut, "%d\t%s\t%d\t%d\t%c\n", ni, strChr, nStart, nEnd, chStrand);
			}
		}

		ni++;
	}

	fclose(fpIn);
	fclose(fpOut);

	/* release memory */
	DestroyIntMatrix(pChrLen);

	/* return */
	return PROC_SUCCESS;
}