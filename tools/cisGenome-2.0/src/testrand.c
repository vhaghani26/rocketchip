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
#include "WorkLib.h"

int menu_flexmodule(int argv, char **argc);

int main(int argv, char **argc)
{
	int nLen;
	int nseed;
	int ni;
	int nRand;

	/* init rand */
	srand( (unsigned)time( NULL ) );
	nRand = rand();	
	nseed = nRand%1000;
	
	printf("Rand = %d\n", nRand);
	printf("Seed = %d\n", nseed);
	
	
	rand_u_init(nseed);


	/* ---- */
	/* menu */
	/* ---- */
	printf("Seed = %d\n", nseed);
	for(ni=0; ni<10; ni++)
	{
		printf("%f\n", rand_u());
	}

	/* exit */
	exit(EXIT_SUCCESS);
}