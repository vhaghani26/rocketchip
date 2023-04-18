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

int ESRemoveRedundancy(int argv, char **argc);


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
	ESRemoveRedundancy(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int ESRemoveRedundancy(int argv, char **argc)
{
	/* define */
	char strInFile[LINE_LENGTH];
	char strOutFile[LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	FILE *fpIn;
	FILE *fpOut;
	char strChr[LINE_LENGTH];
	char strLastChr[LINE_LENGTH];
	int nStart,nEnd;
	int nLastStart = -1;
	int nLastEnd = -1;
	int nId;
	char *chp1,*chp2;
	int ni;
	
	/* init */
	if(argv != 3)
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	strcpy(strInFile, argc[1]);
	strcpy(strOutFile, argc[2]);

	strcpy(strLastChr, "NA");

	fpIn = NULL;
	fpIn = fopen(strInFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: cannot open file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: cannot open file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%d %s %d %d", &nId, strChr, &nStart, &nEnd);

		if(strcmp(strChr, strLastChr) == 0)
		{
			if(nStart > nLastEnd)
			{
				fprintf(fpOut, "%s\n", strLine);
				nLastStart = nStart;
				nLastEnd = nEnd;
			}
			else
			{
				nStart = nLastEnd+1;
				if(nEnd > nLastEnd)
				{
					fprintf(fpOut, "%d\t%s\t%d\t%d\t", nId, strChr, nStart, nEnd);
					
					chp1 = strchr(strLine, '\t');
					chp2 = chp1+1;
					chp1 = strchr(chp2, '\t');
					chp2 = chp1+1;
					chp1 = strchr(chp2, '\t');
					chp2 = chp1+1;
					chp1 = strchr(chp2, '\t');
					chp2 = chp1+1;

					fprintf(fpOut, "%s\n", chp2);

					nLastStart = nStart;
					nLastEnd = nEnd;
				}
				else
				{
					/* ignore */
				}
			}
		}
		else
		{
			fprintf(fpOut, "%s\n", strLine);
			strcpy(strLastChr, strChr);
			nLastStart = nStart;
			nLastEnd = nEnd;
		}

	}

	fclose(fpIn);
	fclose(fpOut);


	/* return */
	return PROC_SUCCESS;
}