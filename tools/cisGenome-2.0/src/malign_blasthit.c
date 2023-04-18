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

int menu_malign_blasthit(int argv, char **argc);

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
	menu_malign_blasthit(argv, argc);

	/* exit */
	exit(EXIT_SUCCESS);
}

int menu_malign_blasthit(int argv, char **argc)
{
	/* define */
	char strParamPath[LINE_LENGTH];	

	int nResult;
	
	/* ------------------------------- */
	/*    malign_blasthit              */
	/* ------------------------------- */
	if((argv < 1) || (argv > 2) )
	{
		printf("Error: parameter wrong!\n");
		exit(EXIT_FAILURE);
	}
	else if(argv == 1)
	{
		printf("/* ----------------------------- */\n");
		printf("    malign_blasthit                \n");
		printf(" example: \n");
		printf("    malign_blasthit malign_blasthit_arg.txt \n");
		printf("/* ----------------------------- */\n");
		exit(EXIT_SUCCESS);
	}

	strcpy(strParamPath, argc[1]);
	nResult = MAlign_BlastHit_Main(strParamPath);
	
	/* return */
	return nResult;
}