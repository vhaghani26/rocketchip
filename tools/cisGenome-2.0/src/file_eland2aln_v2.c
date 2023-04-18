#include "stdlib.h"
#include "stdio.h"
#include "string.h"
#include "math.h"
#include "limits.h"
#include "time.h"

#define LONG_LINE_LENGTH 65535
#define MED_LINE_LENGTH 1024
#define LINE_LENGTH 255
#define SAMP_NUM 31

int File_Eland2Aln_Main(char strInputFile[], char strOutputFile[], int nS);
char* StrTrimLeft(char strLongLine[]);
char* StrTrimRight(char strLongLine[]);

int main(int argv, char **argc)
{
	/* define */
	char strInputPath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	int nS = 12;
	int nResult;
	int iOK;
	int oOK;
	int ni;

	/* ------------------------------- */
	/*       file_eland2aln            */
	/* -i input                        */
	/* -o output                       */
	/* -s shift                        */
	/* ------------------------------- */
	if(argv < 1)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("          file_eland2aln           \n");
		printf(" -i input      \n");
		printf(" -o output     \n");
		printf(" -s number of bp to shift the start position (default: 12bp) \n");
		printf(" example: \n");
		printf("   file_eland2aln -i seq.eland.out -o seq.aln \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	iOK = 0;
	oOK = 0;

	ni = 1;
	while(ni < argv)
	{
		if(strcmp(argc[ni], "-i") == 0)
		{
			ni++;
			strcpy(strInputPath, argc[ni]);
			iOK = 1;
		}
		else if(strcmp(argc[ni], "-o") == 0)
		{
			ni++;
			strcpy(strOutputPath, argc[ni]);
			oOK = 1;
		}
		else if(strcmp(argc[ni], "-s") == 0)
		{
			ni++;
			nS = atoi(argc[ni]);
		}
		else 
		{
			printf("Error: unknown parameters!\n");
			exit(EXIT_FAILURE);
		}

		ni++;
	}

	if((iOK == 0) || (oOK == 0) )
	{
		printf("Error: Input Parameter not correct!\n");
		exit(EXIT_FAILURE);
	}
	else
	{
		nResult = File_Eland2Aln_Main(strInputPath, strOutputPath, nS);
	}

	return nResult;
}

int File_Eland2Aln_Main(char strInputFile[], char strOutputFile[], int nS)
{
	/* define */
	char strLine[LONG_LINE_LENGTH];
	FILE *fpIn;
	FILE *fpOut;
	
	char strID[MED_LINE_LENGTH];
	char strSeq[MED_LINE_LENGTH];
	char strQual[MED_LINE_LENGTH];
	char strP1[MED_LINE_LENGTH];
	char strP2[MED_LINE_LENGTH];
	char strP3[MED_LINE_LENGTH];
	char strP4[MED_LINE_LENGTH];
	char strP5[MED_LINE_LENGTH];
	char strP6[MED_LINE_LENGTH];
	char strP7[MED_LINE_LENGTH];
	char chStrand;
	char strChr[MED_LINE_LENGTH];
	/* char strC[MED_LINE_LENGTH]; */
	int nPos;
	char *chp2;
	int nSave;
	/* int len_strP2; */

	fpIn = NULL;
	fpIn = fopen(strInputFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: cannot open input file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutputFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		nSave = 1;
		if (13 == sscanf(strLine, "%s %s %s %s %s %s %s %s %s %s %s %d %c", strID, strP1, strP2, strP3, strP4, strP5, strP6, strP7, strSeq, strQual, strChr, &nPos, &chStrand)) {
			chp2 = strrchr(strChr, '.');
			if(chp2 != NULL)
				*chp2 = '\0';
			/* strcpy(strC, strChr+1);
			sprintf(strChr, "chr%s", strC); */
		
			fprintf(fpOut, "%s\t%d\t%c\n", strChr, nPos+nS, chStrand);
		}
	}

	fclose(fpIn);
	fclose(fpOut);

	return 0;
}


/* ----------------------------------------------------------------------- */ 
/*                         char* StrTrimRight()                            */
/* This function is used for deleting '\n' ' ' '\t' '\r' on the right of a */
/* string.                                                                 */
/* ----------------------------------------------------------------------- */ 
char* StrTrimRight(char strLongLine[])
{
	char *endp;
	if(strLongLine != NULL)
	{
		endp=strchr(strLongLine, '\0');
		while(endp != strLongLine)
		{
			endp--;
			if((*endp == ' ') || (*endp == '\t') || (*endp == '\n') || (*endp == '\r'))
			{
				*endp = '\0';
				endp=strchr(strLongLine, '\0');
			}
			else
			{
				break;
			}
		}
	}

	return strLongLine;
}

/* ----------------------------------------------------------------------- */ 
/*                         char* StrTrimLeft()                             */
/* This function is used for deleting '\n' ' ' '\t' '\r' on the left of a  */
/* string.                                                                 */
/* ----------------------------------------------------------------------- */ 
char* StrTrimLeft(char strLongLine[])
{
	char *startp;
	long nLength,ni;

	if(strLongLine != NULL)
	{
		nLength = strlen(strLongLine);
		if(nLength>0)
		{
			startp=strLongLine;
			
			/* search for start point */
			ni = 0;
			while(*startp != '\0')
			{
				if((*startp == ' ') || (*startp == '\t') || (*startp == '\n') || (*startp == '\r'))
				{
					ni++;
					startp++;
				}
				else
				{
					break;
				}
			}

			/* memmove */
			if(ni == nLength)
			{
				strLongLine[0] = '\0';
			}
			else if(ni != 0)
			{
				memmove(strLongLine, strLongLine+ni, nLength-ni+1);
			}
		}
	}

	return strLongLine;
}
