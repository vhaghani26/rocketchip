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
	char strP1[MED_LINE_LENGTH];
	char strP2[MED_LINE_LENGTH];
	char strP3[MED_LINE_LENGTH];
	char strP4[MED_LINE_LENGTH];
	char chStrand;
	char strChr[MED_LINE_LENGTH];
	char strC[MED_LINE_LENGTH];
	int nPos;
	char *chp2;
	int nSave;
	int len_strP2;

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
		if (6 <= sscanf(strLine, "%s %s %s %s %s %s %s %d %c", strID, strSeq, strP1, strP2, strP3, strP4, strChr, &nPos, &chStrand)) {
			chp2 = strrchr(strChr, '.');
			if(chp2 != NULL)
				*chp2 = '\0';
			strcpy(strC, strChr);
		
			if(strP1[0] == 'U')
				fprintf(fpOut, "%s\t%d\t%c\n", strC, nPos+nS, chStrand);
		} else if (4 == sscanf(strLine, "%s %s %s %s", strID, strSeq, strP1, strP2)) {
			if (strlen(strP1)==5 && strP1[1] == ':' && strP1[3] == ':') {
				if (strP1[0] == '1' && strP1[2] == '0' 
				|| strP1[0] == '0' && strP1[2] == '1' && strP1[4] == '0'
				|| strP1[0] == '0' && strP1[2] == '0' && strP1[4] == '1') {					
					len_strP2 = strlen(strP2);
					if (len_strP2 <= 2) continue;
					chStrand = strP2[len_strP2-2];
					strP2[len_strP2-2]= '\0';
					chp2 = strrchr(strP2, ':');
					if (chp2 == NULL) continue;
					sscanf(chp2+1, "%d", &nPos);
					*chp2 = '\0';
					chp2 = strrchr(strP2, '.');
					if(chp2 != NULL) 
						*chp2 = '\0';
					strcpy(strC, strP2);
					fprintf(fpOut, "%s\t%d\t%c\n", strC, nPos+nS, chStrand);
				}
			}
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
