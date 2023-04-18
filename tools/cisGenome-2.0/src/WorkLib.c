/* ----------------------------------------------------------------------- */
/*  WorkLib.h : interface of the research projects                         */
/*  Author : Ji HongKai ; Time: 2004.10                                    */
/* ----------------------------------------------------------------------- */

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

/* affy tiling match */
int menu_affychr2122_probematch()
{
	FILE *fpIn1;
	FILE *fpIn2;
	FILE *fpOut;
	char strFileName[LINE_LENGTH];
	char strDataPath[LINE_LENGTH];
	char strLine1[LINE_LENGTH];
	char strLine2[LINE_LENGTH];
	int nPos1,nPos2;

	strcpy(strDataPath, "C:\\Projects\\research_harvard\\affy_project\\analysis\\tiling_paper\\");
	
	sprintf(strFileName, "%scMyc_A_samp2gtrans01_pvalue.txt", strDataPath);
	fpIn1 = NULL;
	fpIn1 = fopen(strFileName, "r");

	sprintf(strFileName, "%scMyc_A_ref3naive.ori", strDataPath);
	fpIn2 = NULL;
	fpIn2 = fopen(strFileName, "r");

	sprintf(strFileName, "%scMyc_A_ref3naive_matchaffy.ori", strDataPath);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");


	while(fgets(strLine1, LINE_LENGTH, fpIn1) != NULL)
	{
		StrTrimLeft(strLine1);
		StrTrimRight(strLine1);
		if(strLine1[0] == '\0')
			continue;

		sscanf(strLine1, "%d ", &nPos1);
		while(fgets(strLine2, LINE_LENGTH, fpIn2) != NULL)
		{
			StrTrimLeft(strLine2);
			StrTrimRight(strLine2);
			if(strLine2[0] == '\0')
				continue;

			sscanf(strLine2, "%d ", &nPos2);
			if(nPos2 == nPos1)
				break;
			else
			{
				printf("%d-%d\n", nPos1, nPos2);
			}
		}

		fprintf(fpOut, "%s\n", strLine2);
	}



	fclose(fpIn1);
	fclose(fpIn2);
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}

/* sonic hedgehog insitu */
int menu_shhinsitu()
{
	char strInSituPath[LINE_LENGTH];
	char strAnnotPath[LINE_LENGTH];
	char strOutPath[LINE_LENGTH];
	char strProbe[LINE_LENGTH];
	char strProbe1[LINE_LENGTH];
	int nOn;
	int ni;
	char strLine[LONG_LINE_LENGTH];
	FILE *fpIn1;
	FILE *fpIn2;
	FILE *fpOut;

	/* match */
	strcpy(strInSituPath, "C:\\Projects\\research_harvard\\genomelab_project\\Projects\\shh\\insitu.txt");
	/* strcpy(strAnnotPath, "C:\\Projects\\research_harvard\\genomelab_project\\Projects\\shh\\shh_U74_geneinfo.txt");
	strcpy(strOutPath, "C:\\Projects\\research_harvard\\genomelab_project\\Projects\\shh\\shh_U74_insitu.txt");*/

	strcpy(strAnnotPath, "C:\\Projects\\research_harvard\\genomelab_project\\Projects\\shh\\shhdata1.txt");
	strcpy(strOutPath, "C:\\Projects\\research_harvard\\genomelab_project\\Projects\\shh\\shh_all_insitu.txt");
	
	fpIn1 = NULL;
	fpIn1 = fopen(strInSituPath, "rt");
	if(fpIn1 == NULL)
	{
		printf("Error!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutPath, "wt");
	if(fpOut == NULL)
	{
		printf("Error!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn1)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%s %d", strProbe, &nOn);

		fpIn2 = NULL;
		fpIn2 = fopen(strAnnotPath, "rt");
		if(fpIn2 == NULL)
		{
			printf("Error!\n");
			exit(EXIT_FAILURE);
		}

		fgets(strLine, LONG_LINE_LENGTH, fpIn2);
		ni=1;
		while(fgets(strLine, LONG_LINE_LENGTH, fpIn2)!=NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;

			sscanf(strLine, "%s ", strProbe1);
			if(strcmp(strProbe, strProbe1) == 0)
			{
				fprintf(fpOut, "%d\t%d\n", ni, nOn);
			}
			ni++;
		}

		fclose(fpIn2);
	}
	fclose(fpIn1);
	fclose(fpOut);

	return PROC_SUCCESS;
}

/* mapping RPL RNA to genome */
int menu_MapRPL(char strAlnPath[], char strRefGenePath[], char strOutPath[], char strSpecies[])
{
	/* define */
	int nRefGeneNum;
	struct tagRefGene **vRefGene;
	FILE *fpIn;
	FILE *fpOut;
	char strLine[LONG_LINE_LENGTH];

	int nLastChrom;
	int nChromNum;
	int vChromDelimits[MAX_CHROMOSOME_NUM][2];
	int ni,nj,nk,nD1;

	char strTemp1[LINE_LENGTH];
	char strTemp2[LINE_LENGTH];
	char strTemp3[LINE_LENGTH];
	int nAlnNum,nQS,nQE,nQLen,nTS,nTE,nTLen;
	char strTChr[LINE_LENGTH];
	char chTStrand;
	char strIdentity[LINE_LENGTH];
	double dIdentity;
	int nTChr;

	int npp1,npp2;
	int nMapResult;
	int nExonOverLap, nIntronOverLap, nReverseStrand;
	int nz1,nz2,nz3;
	int nlen;

	/* ------------------------ */
	/* load refgene annotations */
	/* ------------------------ */

	/* get refgene number */
	fpIn = NULL;
	fpIn = fopen(strRefGenePath, "rt");
	if(fpIn == NULL)
	{
		printf("Error: menu_MapRPL, cannot open the refgene file!\n");
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
		printf("Error: menu_MapRPL, cannot create the memory for refgenes!\n");
		exit(EXIT_FAILURE);
	}

	/* load all refgene annotations */
	fpIn = NULL;
	fpIn = fopen(strRefGenePath, "rt");
	if(fpIn == NULL)
	{
		printf("Error: menu_MapRPL, cannot open the refgene file!\n");
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
		printf("Error: menu_MapRPL, refgene number not match!\n");
		exit(EXIT_FAILURE);
	}

	/* ------------------------ */
	/* map RPL annotations */
	/* ------------------------ */
	fpIn = NULL;
	fpIn = fopen(strAlnPath, "rt");
	if(fpIn == NULL)
	{
		printf("Error: menu_MapRPL, cannot load alignment!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutPath, "wt");
	if(fpOut == NULL)
	{
		printf("Error: menu_MapRPL, cannot output!\n");
		exit(EXIT_FAILURE);
	}


	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		/*if(strcmp(strLine, "browser details 33866947-33867249(303)   163    21   256   303  85.6%     3   +  196389893 196390119    227") == 0)
		{
			nAlnNum = 0;
		}*/

		sscanf(strLine, "%s %s %s %d %d %d %d %s %s %c %d %d %d",
			strTemp1, strTemp2, strTemp3, &nAlnNum, &nQS, &nQE, 
			&nQLen, strIdentity, strTChr, &chTStrand, &nTS, &nTE, &nTLen);

		sprintf(strTemp1, "chr%s", strTChr);
		nTChr = Genome_ChromosomeName_To_Index(strTemp1, strSpecies);
		nlen = strlen(strIdentity);
		strIdentity[nlen-1] = '\0';
		dIdentity = atof(strIdentity)/100.0;

		if((nAlnNum < 20) || ( dIdentity < 0.7 ) )
			continue;

		npp1 = vChromDelimits[nTChr-1][0];
		npp2 = vChromDelimits[nTChr-1][1];

		for(nj=npp1; nj<=npp2; nj++)
		{
			
			/*if(strcmp(vRefGene[nj]->strName, "NM_012287") == 0)
			{
				nMapResult = 0;
			}*/

			nMapResult = 0;

			if(vRefGene[nj]->nChrom != nTChr)
				continue;

			if( ((vRefGene[nj]->nTxStart <= nTS) && (nTS <= vRefGene[nj]->nTxEnd)) || 
				((vRefGene[nj]->nTxStart <= nTE) && (nTE <= vRefGene[nj]->nTxEnd)) )
			{
				if(vRefGene[nj]->chStrand == chTStrand)
				{
					nReverseStrand = 0;
				}
				else
				{
					nReverseStrand = 1;
				}

				nExonOverLap = 0;
				nIntronOverLap = 0;
				for(nk=0; nk<vRefGene[nj]->nExonCount; nk++)
				{
					nz1 = IMGETAT(vRefGene[nj]->pmatExonStartsEnds, nk, 0);
					nz2 = IMGETAT(vRefGene[nj]->pmatExonStartsEnds, nk, 1);
					if( ((nz1<=nTS) && (nTS<=nz2)) || ((nz1<=nTE) && (nTE<=nz2)) )
					{
						nExonOverLap = 1;
					}

					if(nk != 0)
					{
						nz3 = IMGETAT(vRefGene[nj]->pmatExonStartsEnds, (nk-1), 1);
						if( ((nz3<=nTS) && (nTS<=nz1)) || ((nz3<=nTE) && (nTE<=nz1)) )
						{
							nIntronOverLap = 1;
						}
					}
				}

				fprintf(fpOut, "%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
					vRefGene[nj]->strName, vRefGene[nj]->nChrom, 
					vRefGene[nj]->nTxStart, vRefGene[nj]->nTxEnd,
					nReverseStrand, nExonOverLap, nIntronOverLap,
					nQS, nQE, nTS, nTE);
			}
		}
	}

	fclose(fpIn);
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
	
	
	/* return */
	return PROC_SUCCESS;
}

/* map cMyc binding site to chromosome 21 and 22 */
int Map_cMyc_Chr2122()
{
	/* define */
	char strGenomePath[LINE_LENGTH];
	char strSeqFile[LINE_LENGTH];
	char strOutFile[LINE_LENGTH];
	FILE *fpSeq;
	FILE *fpOut;
	int nChrLen;
	struct INTMATRIX *pChrLen;
	int nEffecLen = 0;
	int nSiteNum = 0;
	struct SEQMOTIF *pSeqMtf;
	int nMotifLen = 6;
	int nChr = 21;
	unsigned char *pBase;
	int nMisNum; 
	int nBadLen;
	char strSite[LINE_LENGTH];
	int nInitLen,nP1,nP2,nRemainLen;
	int ni,nj;

	/* init */
	sprintf(strGenomePath, "C:\\Projects\\research_harvard\\genomelab_project\\genomes\\human\\b34_hg16\\");

	pChrLen = NULL;
	sprintf(strSeqFile, "%schrlen.txt", strGenomePath);
	pChrLen = IMLOAD(strSeqFile);
	if(pChrLen == NULL)
	{
		printf("cannot load chromosome length\n");
		exit(EXIT_FAILURE);
	}

	sprintf(strSeqFile, "%schr%d.sq", strGenomePath, nChr);
	sprintf(strOutFile, "cMyc_chr%d.txt", nChr);
	nChrLen = pChrLen->pMatElement[nChr-1];


	pSeqMtf = NULL;
	nInitLen = GENOME_MED_CONTIG_LEN+nMotifLen-1;

	/* create 1 sequence, 1 background scores, 1 filter information */
	pSeqMtf = SEQMTFCREATE(0, 1, 0, 0);
	if(pSeqMtf == NULL)
	{
		printf("Error: cannot create sequence/motif complex!\n");
		exit(EXIT_FAILURE);
	}
	/* sequence */
	SEQMTFCREATESEQ(pSeqMtf, 0, nInitLen);
	/* forward mc background */
	/* SEQMTFCREATESCORE(pSeqMtf, 0, nInitLen); */
	/* backward mc background */
	/* SEQMTFCREATESCORE(pSeqMtf, 1, nInitLen); */


	/* open sequence file and get the sequence */
	fpSeq = NULL;
	fpSeq = fopen(strSeqFile, "rb");
	if(fpSeq == NULL)
	{
		printf("Warning: cannot open sequence file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Warning: cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	nEffecLen = 0;

	/* process contig by contig */
	nP1 = 0;
	nP2 = nP1+GENOME_MED_CONTIG_LEN-1;
	nRemainLen = nChrLen;
	if(nRemainLen < (GENOME_MED_CONTIG_LEN-nMotifLen+1))
	{
		nP2 = nP1+nRemainLen-nMotifLen;
	}

	while(nRemainLen > 0)
	{
		/* load sequence */
		MapcMyc_LOADSEQFROMGENOME(pSeqMtf, 0, fpSeq, nP1, (nP2+nMotifLen-1), '+');
		pBase = pSeqMtf->ppSeq[0]->pMatElement;
		for(ni=nP1; ni<=nP2; ni++)
		{
			nBadLen = 0;
			for(nj=0; nj<nMotifLen; nj++)
			{
				if(pBase[nj] >= 4)
				{
					nBadLen++;
					break;
				}
			}

			if(nBadLen == 0)
			{
				nEffecLen++;
			}
			else
			{
				pBase++;
				continue;
			}


			/* check pos strand */
			nMisNum = 0;
			
			if(pBase[0] != 1)
			{
				nMisNum = nMotifLen;
			}
			if(pBase[1] != 0)
			{
				nMisNum = nMotifLen;
			}
			if(pBase[2] != 1) 
			{
				if(pBase[2] == 3)
				{
					nMisNum++;
				}
				else
				{
					nMisNum = nMotifLen;
				}
			}
			if(pBase[3] != 2)
			{
				nMisNum = nMotifLen;
			}
			if(pBase[4] != 3)
			{
				if(pBase[4] == 1)
				{
					nMisNum++;
				}
				else
				{
					nMisNum = nMotifLen;
				}
			}
			if(pBase[5] != 2)
			{
				nMisNum = nMotifLen;
			}

			if(nMisNum <= 1)
			{
				/* write site */
				for(nj=0; nj<nMotifLen; nj++)
				{
					switch(pBase[nj])
					{
						case 0: strSite[nj] = 'A';
							break;
						case 1: strSite[nj] = 'C';
							break;
						case 2: strSite[nj] = 'G';
							break;
						case 3: strSite[nj] = 'T';
							break;
						default: strSite[nj] = 'N';
					}
				}
				strSite[nj] = '\0';

				fprintf(fpOut, "%d\t%d\t+\t%s\n", nChr, ni, strSite);

				nSiteNum++;
				pBase++;
				continue;
			}

			/* check neg strand */
			nMisNum = 0;
			
			if(pBase[5] != 2)
			{
				nMisNum = nMotifLen;
			}
			if(pBase[4] != 3)
			{
				nMisNum = nMotifLen;
			}
			if(pBase[3] != 2)
			{
				if(pBase[3] == 0)
				{
					nMisNum++;
				}
				else
				{
					nMisNum = nMotifLen;
				}
			}
			if(pBase[2] != 1)
			{
				nMisNum = nMotifLen;
			}
			if(pBase[1] != 0)
			{
				if(pBase[1] == 2)
				{
					nMisNum++;
				}
				else
				{
					nMisNum = nMotifLen;
				}
			}
			if(pBase[0] != 1)
			{
				nMisNum = nMotifLen;
			}
			if(nMisNum <= 1)
			{
				/* write site */
				for(nj=0; nj<nMotifLen; nj++)
				{
					switch(pBase[nMotifLen-1-nj])
					{
						case 0: strSite[nj] = 'T';
							break;
						case 1: strSite[nj] = 'G';
							break;
						case 2: strSite[nj] = 'C';
							break;
						case 3: strSite[nj] = 'A';
							break;
						default: strSite[nj] = 'N';
					}
				}
				strSite[nj] = '\0';

				fprintf(fpOut, "%d\t%d\t-\t%s\n", nChr, ni, strSite);

				nSiteNum++;
			}

			pBase++;
		}

		
		/* count processed bases */
		nRemainLen -= (nP2-nP1+1);
		if(nRemainLen < nMotifLen)
		{
			break;
		}

		/* get next */
		nP1 = nP2+1;
		if(nRemainLen < (GENOME_MED_CONTIG_LEN-nMotifLen+1))
		{
			nP2 = nP1+nRemainLen-nMotifLen;
		}
		else
		{
			nP2 = nP1+GENOME_MED_CONTIG_LEN-1;
		}
	}

	fprintf(fpOut, "\n===\n");
	fprintf(fpOut, "%d\t%d\n", nSiteNum, nEffecLen);

	/* close files */
	fclose(fpSeq);
	fclose(fpOut);

	/* release memory */
	DestroyIntMatrix(pChrLen);
	SEQMTFDESTROY(pSeqMtf);

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*      MapcMyc_LOADSEQFROMGENOME: load sequences from genome.             */
/* ----------------------------------------------------------------------- */ 
int MapcMyc_LOADSEQFROMGENOME(struct SEQMOTIF *pSeqMtf, int nSeqId, 
							 FILE *fpSeq, int nStart, int nEnd, char chStrand)
{
	/* define */
	unsigned char *pBase;
	int nLen;
	int nP1,nP2,nR1,nR2;
	int numread;
	int ni,nk;
	unsigned char bBase,bChar;

	/* init */
	if( (pSeqMtf == NULL) || (fpSeq == NULL) )
	{
		printf("Warning: MapcMyc_LOADSEQFROMGENOME, null seqmotif or seqfile pointer!\n");
		return PROC_FAILURE;
	}
	if( (nSeqId<0) || (nSeqId >= pSeqMtf->nSeqNum) )
	{
		printf("Warning: MapcMyc_LOADSEQFROMGENOME, sequence id out of range!\n");
		return PROC_FAILURE;
	}
	nLen = (nEnd-nStart+1);
	if( (nLen<0) || (nLen >  pSeqMtf->ppSeq[nSeqId]->nWidth ) )
	{
		printf("Warning: MapcMyc_LOADSEQFROMGENOME, the length of the sequence to be retrieved is larger than the length allowed by seqmotif!\n");
		return PROC_FAILURE;
	}
	if((chStrand != '+') && (chStrand != '-'))
	{
		printf("Warning: MapcMyc_LOADSEQFROMGENOME, strand information not specified!\n");
	}


	/* get seq */
	pBase = pSeqMtf->ppSeq[nSeqId]->pMatElement;
	nP1 = (int)(nStart/2);
	nR1 = nStart%2;
	nP2 = (int)(nEnd/2);
	nR2 = nEnd%2;
	nk = 0;

	if( fseek( fpSeq, nP1, SEEK_SET ) != 0 )
	{
		printf("Error: MapcMyc_LOADSEQFROMGENOME, cannot locate the required sequence!\n");
		exit(EXIT_FAILURE);
	}


	/* first base */
	numread = fread( &bChar, sizeof( unsigned char ), 1, fpSeq );
	if(nR1 == 0)
	{
		bBase = bChar >> 4;
		switch(bBase)
		{
			case 0: *pBase = 0;
				break;
			case 1: *pBase = 1;
				break;
			case 2: *pBase = 2;
				break;
			case 3: *pBase = 3;
				break;
			case 4: *pBase = 0;
				break;
			case 5: *pBase = 1;
				break;
			case 6: *pBase = 2;
				break;
			case 7: *pBase = 3;
				break;
			default: *pBase = 4;
		}
		pBase++;
		nk++;
	}

	bBase = bChar & 0x0F;
	switch(bBase)
	{
		case 0: *pBase = 0;
			break;
		case 1: *pBase = 1;
			break;
		case 2: *pBase = 2;
			break;
		case 3: *pBase = 3;
			break;
		case 4: *pBase = 0;
			break;
		case 5: *pBase = 1;
			break;
		case 6: *pBase = 2;
			break;
		case 7: *pBase = 3;
			break;
		default: *pBase = 4;
	}
	pBase++;
	nk++;

	/* middle bases */
	for(ni=(nP1+1); ni<nP2; ni++)
	{
		numread = fread( &bChar, sizeof( unsigned char ), 1, fpSeq );
		bBase = bChar >> 4;
		switch(bBase)
		{
			case 0: *pBase = 0;
				break;
			case 1: *pBase = 1;
				break;
			case 2: *pBase = 2;
				break;
			case 3: *pBase = 3;
				break;
			case 4: *pBase = 0;
				break;
			case 5: *pBase = 1;
				break;
			case 6: *pBase = 2;
				break;
			case 7: *pBase = 3;
				break;
			default: *pBase = 4;
		}
		pBase++;
		nk++;

		bBase = bChar & 0x0F;
		switch(bBase)
		{
			case 0: *pBase = 0;
				break;
			case 1: *pBase = 1;
				break;
			case 2: *pBase = 2;
				break;
			case 3: *pBase = 3;
				break;
			case 4: *pBase = 0;
				break;
			case 5: *pBase = 1;
				break;
			case 6: *pBase = 2;
				break;
			case 7: *pBase = 3;
				break;
			default: *pBase = 4;
		}
		pBase++;
		nk++;
	}
	
	/* last base */
	numread = fread( &bChar, sizeof( unsigned char ), 1, fpSeq );
	bBase = bChar >> 4;
	switch(bBase)
	{
		case 0: *pBase = 0;
			break;
		case 1: *pBase = 1;
			break;
		case 2: *pBase = 2;
			break;
		case 3: *pBase = 3;
			break;
		case 4: *pBase = 0;
			break;
		case 5: *pBase = 1;
			break;
		case 6: *pBase = 2;
			break;
		case 7: *pBase = 3;
			break;
		default: *pBase = 4;
	}
	pBase++;
	nk++;

	if(nR2 == 1)
	{
		bBase = bChar & 0x0F;
		switch(bBase)
		{
			case 0: *pBase = 0;
				break;
			case 1: *pBase = 1;
				break;
			case 2: *pBase = 2;
				break;
			case 3: *pBase = 3;
				break;
			case 4: *pBase = 0;
				break;
			case 5: *pBase = 1;
				break;
			case 6: *pBase = 2;
				break;
			case 7: *pBase = 3;
				break;
			default: *pBase = 4;
		}
		pBase++;
		nk++;
	}

	if(nk != nLen)
	{
		printf("Error: MapcMyc_LOADSEQFROMGENOME, sequence length not match!\n");
		exit(EXIT_FAILURE);
	}

	/* get reverse complement if necessary */
	if(chStrand == '-')
	{
		pBase = pSeqMtf->ppSeq[nSeqId]->pMatElement;
		ni = 0;
		nk = nLen-1;
		while(nk>=ni)
		{
			bBase = pBase[ni];
			if(pBase[nk]<4)
			{
				pBase[ni] = 3-pBase[nk];
			}
			else
			{
				pBase[ni] = 4;
			}

			if(bBase<4)
			{
				pBase[nk] = 3-bBase;
			}
			else
			{
				pBase[nk] = 4;
			}

			ni++;
			nk--;
		}
	}

	/* return */
	return PROC_SUCCESS;
}



/* count cMyc in target region */
int Count_cMyc_In_TargetRegion_Main()
{
	/* define */
	char strGenomePath[LINE_LENGTH];
	char strRegionFile[LINE_LENGTH];
	char strSiteFile[LINE_LENGTH];
	char strJobName[LINE_LENGTH];
	char strWorkPath[LINE_LENGTH];
	char strOutPath[LINE_LENGTH];
	int ni;
	int nChr = 21;
	int nMotifLen = 6;
	

	/* sprintf(strGenomePath, "C:\\Projects\\research_harvard\\genomelab_project\\genomes\\human\\b34_hg16\\");
	sprintf(strSiteFile, "C:\\Projects\\research_harvard\\affy_project\\analysis\\tiling_paper\\cMycSites\\cMyc_chr%d.txt", nChr);
	strcpy(strWorkPath, "C:\\Projects\\research_harvard\\affy_project\\analysis\\tiling_paper\\"); */

	sprintf(strGenomePath, "/data/genomes/human/b33_hg15/");
	sprintf(strSiteFile, "cMyc_chr%d.txt", nChr);
	strcpy(strWorkPath, "./");
	
	/* naive */
	for(ni=1; ni<=3; ni++)
	{
		sprintf(strJobName, "cMyc_A_anova_sample2r2naive_%d_f_mean", ni);
		sprintf(strRegionFile, "%s%s.rstat", strWorkPath, strJobName);
		sprintf(strOutPath, "%s%s.site", strWorkPath, strJobName);
		Count_cMyc_In_TargetRegion(strGenomePath, nChr, nMotifLen,
			strSiteFile, strRegionFile, strOutPath);
	}
	for(ni=1; ni<=3; ni++)
	{
		sprintf(strJobName, "cMyc_A_anova_sample3r2naive_%d_f_mean", ni);
		sprintf(strRegionFile, "%s%s.rstat", strWorkPath, strJobName);
		sprintf(strOutPath, "%s%s.site", strWorkPath, strJobName);
		Count_cMyc_In_TargetRegion(strGenomePath, nChr, nMotifLen,
			strSiteFile, strRegionFile, strOutPath);
	}
	for(ni=1; ni<=3; ni++)
	{
		sprintf(strJobName, "cMyc_A_anova_sample3r4naive_%d_f_mean", ni);
		sprintf(strRegionFile, "%s%s.rstat", strWorkPath, strJobName);
		sprintf(strOutPath, "%s%s.site", strWorkPath, strJobName);
		Count_cMyc_In_TargetRegion(strGenomePath, nChr, nMotifLen,
			strSiteFile, strRegionFile, strOutPath);
	}

	/* varsh */
	for(ni=1; ni<=3; ni++)
	{
		sprintf(strJobName, "cMyc_A_anova_sample2r2varsh_%d_f_mean", ni);
		sprintf(strRegionFile, "%s%s.rstat", strWorkPath, strJobName);
		sprintf(strOutPath, "%s%s.site", strWorkPath, strJobName);
		Count_cMyc_In_TargetRegion(strGenomePath, nChr, nMotifLen,
			strSiteFile, strRegionFile, strOutPath);
	}
	for(ni=1; ni<=3; ni++)
	{
		sprintf(strJobName, "cMyc_A_anova_sample3r2varsh_%d_f_mean", ni);
		sprintf(strRegionFile, "%s%s.rstat", strWorkPath, strJobName);
		sprintf(strOutPath, "%s%s.site", strWorkPath, strJobName);
		Count_cMyc_In_TargetRegion(strGenomePath, nChr, nMotifLen,
			strSiteFile, strRegionFile, strOutPath);
	}
	for(ni=1; ni<=3; ni++)
	{
		sprintf(strJobName, "cMyc_A_anova_sample3r4varsh_%d_f_mean", ni);
		sprintf(strRegionFile, "%s%s.rstat", strWorkPath, strJobName);
		sprintf(strOutPath, "%s%s.site", strWorkPath, strJobName);
		Count_cMyc_In_TargetRegion(strGenomePath, nChr, nMotifLen,
			strSiteFile, strRegionFile, strOutPath);
	}

	/* hmm */
	for(ni=1; ni<=3; ni++)
	{
		sprintf(strJobName, "cMyc_A_anova_sample2r2varsh_%d_f_refbind_pl", ni);
		sprintf(strRegionFile, "%s%s.rstat", strWorkPath, strJobName);
		sprintf(strOutPath, "%s%s.site", strWorkPath, strJobName);
		Count_cMyc_In_TargetRegion(strGenomePath, nChr, nMotifLen,
			strSiteFile, strRegionFile, strOutPath);
	}
	for(ni=1; ni<=3; ni++)
	{
		sprintf(strJobName, "cMyc_A_anova_sample3r2varsh_%d_f_refbind_pl", ni);
		sprintf(strRegionFile, "%s%s.rstat", strWorkPath, strJobName);
		sprintf(strOutPath, "%s%s.site", strWorkPath, strJobName);
		Count_cMyc_In_TargetRegion(strGenomePath, nChr, nMotifLen,
			strSiteFile, strRegionFile, strOutPath);
	}
	for(ni=1; ni<=3; ni++)
	{
		sprintf(strJobName, "cMyc_A_anova_sample3r4varsh_%d_f_refbind_pl", ni);
		sprintf(strRegionFile, "%s%s.rstat", strWorkPath, strJobName);
		sprintf(strOutPath, "%s%s.site", strWorkPath, strJobName);
		Count_cMyc_In_TargetRegion(strGenomePath, nChr, nMotifLen,
			strSiteFile, strRegionFile, strOutPath);
	}


	/* gtrans */
	for(ni=1; ni<=3; ni++)
	{
		sprintf(strJobName, "cMyc_A_mix_sample2_gtrans_%d_pvalue_f", ni);
		sprintf(strRegionFile, "%s%s.rstat", strWorkPath, strJobName);
		sprintf(strOutPath, "%s%s.site", strWorkPath, strJobName);
		Count_cMyc_In_TargetRegion(strGenomePath, nChr, nMotifLen,
			strSiteFile, strRegionFile, strOutPath);
	}
	for(ni=1; ni<=3; ni++)
	{
		sprintf(strJobName, "cMyc_A_mix_sample3r2_gtrans_%d_maxpvalue_f", ni);
		sprintf(strRegionFile, "%s%s.rstat", strWorkPath, strJobName);
		sprintf(strOutPath, "%s%s.site", strWorkPath, strJobName);
		Count_cMyc_In_TargetRegion(strGenomePath, nChr, nMotifLen,
			strSiteFile, strRegionFile, strOutPath);
	}
	for(ni=1; ni<=3; ni++)
	{
		sprintf(strJobName, "cMyc_A_mix_sample3r4_gtrans_%d_maxpvalue_f", ni);
		sprintf(strRegionFile, "%s%s.rstat", strWorkPath, strJobName);
		sprintf(strOutPath, "%s%s.site", strWorkPath, strJobName);
		Count_cMyc_In_TargetRegion(strGenomePath, nChr, nMotifLen,
			strSiteFile, strRegionFile, strOutPath);
	}
	
	/* return */
	return PROC_SUCCESS;
}


int Count_cMyc_In_TargetRegion(char strGenomePath[], int nChr, int nMotifLen, 
							   char strSiteFile[], 
							   char strRegionFile[], char strOutFile[])
{
	/* define */
	FILE *fpIn;
	FILE *fpSeq;
	FILE *fpOut;
	FILE *fpSite;

	char strSeqFile[LINE_LENGTH];
	char strLine[MED_LINE_LENGTH];
	
	int nRegionNum = 0;
	struct DOUBLEMATRIX *pRegionStat;

	int nChrLen;
	struct INTMATRIX *pChrLen;
	int nEffecLen = 0;
	int nSiteNum = 0;
	struct SEQMOTIF *pSeqMtf;
	
	unsigned char *pBase;
	int nBadLen;
	int nRemainLen;
	int ni,nj,nk;
	double dStart, dEnd, dMin, dMean;
	int nStart, nEnd;
	int nSiteChr,nSitePos;
	double dTemp;

	/* init */
	pChrLen = NULL;
	sprintf(strSeqFile, "%schrlen.txt", strGenomePath);
	pChrLen = IMLOAD(strSeqFile);
	if(pChrLen == NULL)
	{
		printf("cannot load chromosome length\n");
		exit(EXIT_FAILURE);
	}

	sprintf(strSeqFile, "%schr%d.sq", strGenomePath, nChr);
	nChrLen = pChrLen->pMatElement[nChr-1];

	/* get effec len region by region */
	fpIn = NULL;
	fpIn = fopen(strRegionFile, "r");
	if(fpIn == NULL)
	{
		printf("cannot load chromosome length\n");
		exit(EXIT_FAILURE);
	}
	while(fgets(strLine, MED_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		nRegionNum++;
	}
	fclose(fpIn);

	pRegionStat = NULL;
	pRegionStat = CreateDoubleMatrix(nRegionNum, 6);
	if(pRegionStat == NULL)
	{
		printf("cannot get region stat\n");
		exit(EXIT_FAILURE);
	}

	/* open sequence file and get the sequence */
	fpSeq = NULL;
	fpSeq = fopen(strSeqFile, "rb");
	if(fpSeq == NULL)
	{
		printf("Warning: cannot open sequence file!\n");
		exit(EXIT_FAILURE);
	}

	fpIn = NULL;
	fpIn = fopen(strRegionFile, "r");
	if(fpIn == NULL)
	{
		printf("cannot load chromosome length\n");
		exit(EXIT_FAILURE);
	}
	ni = 0;
	while(fgets(strLine, MED_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		
		sscanf(strLine, "%lf %lf %lf %lf", &dStart, &dEnd, &dMin, &dMean);
		
		nStart = (int)dStart;
		nEnd = (int)dEnd;
		nRemainLen = 1000-(nEnd-nStart+1);
		if(nRemainLen > 0)
		{
			nStart -= nRemainLen/2;
			nEnd += nRemainLen/2;
		}

		dStart = (double)nStart;
		dEnd = (double)nEnd;
		DMSETAT(pRegionStat, ni, 0, dStart);
		DMSETAT(pRegionStat, ni, 1, dEnd);
		DMSETAT(pRegionStat, ni, 2, dMin);
		DMSETAT(pRegionStat, ni, 3, dMean);


		pSeqMtf = NULL;
		pSeqMtf = SEQMTFCREATE(0, 1, 0, 0);
		if(pSeqMtf == NULL)
		{
			printf("Error: cannot create sequence/motif complex!\n");
			exit(EXIT_FAILURE);
		}
		/* sequence */
		SEQMTFCREATESEQ(pSeqMtf, 0, nEnd-nStart+nMotifLen);

		MapcMyc_LOADSEQFROMGENOME(pSeqMtf, 0, fpSeq, nStart, nEnd+nMotifLen-1, '+');
		pBase = pSeqMtf->ppSeq[0]->pMatElement;
		nEffecLen = 0;
		for(nk=nStart; nk<=nEnd; nk++)
		{
			nBadLen = 0;
			for(nj=0; nj<nMotifLen; nj++)
			{
				if(pBase[nj] >= 4)
				{
					nBadLen++;
					break;
				}
			}

			if(nBadLen == 0)
			{
				nEffecLen++;
			}

			pBase++;
		}


		DMSETAT(pRegionStat, ni, 4, (double)(nEffecLen));

		SEQMTFDESTROY(pSeqMtf);

		ni++;
	}
	fclose(fpIn);
	fclose(fpSeq);


	fpSite = NULL;
	fpSite = fopen(strSiteFile, "r");
	if(fpSite == NULL)
	{
		printf("cannot open file\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, MED_LINE_LENGTH, fpSite)!= NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strstr(strLine, "===") == strLine)
			break;

		sscanf(strLine, "%d %d", &nSiteChr, &nSitePos);

		for(ni=0; ni<nRegionNum; ni++)
		{
			nStart = (int)(DMGETAT(pRegionStat, ni, 0));
			nEnd = (int)(DMGETAT(pRegionStat, ni, 1));
			if( (nStart <= nSitePos) && (nSitePos <= nEnd) )
			{
				dTemp = DMGETAT(pRegionStat, ni, 5)+1.0;
				DMSETAT(pRegionStat, ni, 5, dTemp);
			}
		}
	}


	fclose(fpSite);

	/* write */
	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("cannot load chromosome length\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nRegionNum; ni++)
	{
		fprintf(fpOut, "%d %d % 9.7e % 9.7e %d %d\n",
			(int)(DMGETAT(pRegionStat, ni, 0)),
			(int)(DMGETAT(pRegionStat, ni, 1)),
			DMGETAT(pRegionStat, ni, 2),
			DMGETAT(pRegionStat, ni, 3),
			(int)(DMGETAT(pRegionStat, ni, 4)),
			(int)(DMGETAT(pRegionStat, ni, 5)));
	}

	fclose(fpOut);

	/* release memory */
	DestroyIntMatrix(pChrLen);
	DestroyDoubleMatrix(pRegionStat);

	/* return */
	return PROC_SUCCESS;
}


/* map cMyc binding site to chromosome 21 and 22 */
int Map_cMyc_Chr2122_Conserved()
{
	/* define */
	char strGenomePath[LINE_LENGTH];
	char strCSPath[LINE_LENGTH];
	char strSeqFile[LINE_LENGTH];
	char strCSFile[LINE_LENGTH];
	char strOutFile[LINE_LENGTH];
	FILE *fpSeq;
	FILE *fpCS;
	FILE *fpOut;
	int nChrLen;
	struct INTMATRIX *pChrLen;
	int nEffecLen = 0;
	int nSiteNum = 0;
	struct SEQMOTIF *pSeqMtf;
	int nMotifLen = 6;
	int nChr = 21;
	unsigned char *pBase,*pC;
	int nMisNum;
	int nBadLen;
	char strSite[LINE_LENGTH];
	int nInitLen,nP1,nP2,nRemainLen;
	int ni,nj;
	double dScore;
	double dCutoff = CMYC_CONS_CUTOFF;
	/* double dCutoff = 18.0; */

	/* init */
	/* sprintf(strGenomePath, "C:\\Projects\\research_harvard\\genomelab_project\\genomes\\human\\b34_hg16\\");
	sprintf(strCSPath, "C:\\Projects\\research_harvard\\genomelab_project\\genomes\\crossspecies\\hg16Mm3Rn3\\hg_cs\\");
	*/

	sprintf(strGenomePath, "/data/genomes/human/b33_hg15/");
	sprintf(strCSPath, "/data/genomes/human/b33_hg15/conservation/cs/");
	
	pChrLen = NULL;
	sprintf(strSeqFile, "%schrlen.txt", strGenomePath);
	pChrLen = IMLOAD(strSeqFile);
	if(pChrLen == NULL)
	{
		printf("cannot load chromosome length\n");
		exit(EXIT_FAILURE);
	}

	sprintf(strSeqFile, "%schr%d.sq", strGenomePath, nChr);
	sprintf(strCSFile, "%schr%d.cs", strCSPath, nChr);
	sprintf(strOutFile, "cMyc_chr%d.txt", nChr);
	nChrLen = pChrLen->pMatElement[nChr-1];


	pSeqMtf = NULL;
	nInitLen = GENOME_MED_CONTIG_LEN+nMotifLen-1;

	/* create 1 sequence, 1 background scores, 1 filter information */
	pSeqMtf = SEQMTFCREATE(0, 2, 0, 0);
	if(pSeqMtf == NULL)
	{
		printf("Error: cannot create sequence/motif complex!\n");
		exit(EXIT_FAILURE);
	}
	/* sequence */
	SEQMTFCREATESEQ(pSeqMtf, 0, nInitLen);
	/* forward mc background */
	SEQMTFCREATESEQ(pSeqMtf, 1, nInitLen);
	/* backward mc background */
	/* SEQMTFCREATESCORE(pSeqMtf, 1, nInitLen); */


	/* open sequence file and get the sequence */
	fpSeq = NULL;
	fpSeq = fopen(strSeqFile, "rb");
	if(fpSeq == NULL)
	{
		printf("Warning: cannot open sequence file!\n");
		exit(EXIT_FAILURE);
	}
	fpCS = NULL;
	fpCS = fopen(strCSFile, "rb");
	if(fpCS == NULL)
	{
		printf("Warning: cannot open conservation file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Warning: cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	nEffecLen = 0;

	/* process contig by contig */
	nP1 = 0;
	nP2 = nP1+GENOME_MED_CONTIG_LEN-1;
	nRemainLen = nChrLen;
	if(nRemainLen < (GENOME_MED_CONTIG_LEN-nMotifLen+1))
	{
		nP2 = nP1+nRemainLen-nMotifLen;
	}

	while(nRemainLen > 0)
	{
		/* load sequence */
		MapcMyc_LOADSEQFROMGENOME_Conserved(pSeqMtf, 0, fpSeq, fpCS, nP1, (nP2+nMotifLen-1), '+');
		
		pBase = pSeqMtf->ppSeq[0]->pMatElement;
		pC = pSeqMtf->ppSeq[1]->pMatElement;

		for(ni=nP1; ni<=nP2; ni++)
		{
			nBadLen = 0;
			dScore = 0.0;
			for(nj=0; nj<nMotifLen; nj++)
			{
				if(pBase[nj] >= 4)
				{
					nBadLen++;
					break;
				}
				dScore += pC[nj];
			}
			dScore /= (double)nMotifLen;

			if(nBadLen == 0)
			/* if((nBadLen == 0) && (dScore >= dCutoff)) */
			{
				nEffecLen++;
			}
			else
			{
				pBase++;
				pC++;
				continue;
			}


			/* check pos strand */
			nMisNum = 0;
			
			if(pBase[0] != 1)
			{
				nMisNum = nMotifLen;
			}
			if(pBase[1] != 0)
			{
				nMisNum = nMotifLen;
			}
			if(pBase[2] != 1) 
			{
				if(pBase[2] == 3)
				{
					nMisNum++;
				}
				else
				{
					nMisNum = nMotifLen;
				}
			}
			if(pBase[3] != 2)
			{
				nMisNum = nMotifLen;
			}
			if(pBase[4] != 3)
			{
				if(pBase[4] == 1)
				{
					nMisNum++;
				}
				else
				{
					nMisNum = nMotifLen;
				}
			}
			if(pBase[5] != 2)
			{
				nMisNum = nMotifLen;
			}

			if( (nMisNum <= 2) && (dScore >= dCutoff) )
			{
				/* write site */
				for(nj=0; nj<nMotifLen; nj++)
				{
					switch(pBase[nj])
					{
						case 0: strSite[nj] = 'A';
							break;
						case 1: strSite[nj] = 'C';
							break;
						case 2: strSite[nj] = 'G';
							break;
						case 3: strSite[nj] = 'T';
							break;
						default: strSite[nj] = 'N';
					}
				}
				strSite[nj] = '\0';

				fprintf(fpOut, "%d\t%d\t+\t%s\n", nChr, ni, strSite);

				nSiteNum++;
				pBase++;
				pC++;
				continue;
			}

			/* check neg strand */
			nMisNum = 0;
			
			if(pBase[5] != 2)
			{
				nMisNum = nMotifLen;
			}
			if(pBase[4] != 3)
			{
				nMisNum = nMotifLen;
			}
			if(pBase[3] != 2)
			{
				if(pBase[3] == 0)
				{
					nMisNum++;
				}
				else
				{
					nMisNum = nMotifLen;
				}
			}
			if(pBase[2] != 1)
			{
				nMisNum = nMotifLen;
			}
			if(pBase[1] != 0)
			{
				if(pBase[1] == 2)
				{
					nMisNum++;
				}
				else
				{
					nMisNum = nMotifLen;
				}
			}
			if(pBase[0] != 1)
			{
				nMisNum = nMotifLen;
			}
			if( (nMisNum <= 2) && (dScore >= dCutoff) )
			{
				/* write site */
				for(nj=0; nj<nMotifLen; nj++)
				{
					switch(pBase[nMotifLen-1-nj])
					{
						case 0: strSite[nj] = 'T';
							break;
						case 1: strSite[nj] = 'G';
							break;
						case 2: strSite[nj] = 'C';
							break;
						case 3: strSite[nj] = 'A';
							break;
						default: strSite[nj] = 'N';
					}
				}
				strSite[nj] = '\0';

				fprintf(fpOut, "%d\t%d\t-\t%s\n", nChr, ni, strSite);

				nSiteNum++;
			}

			pBase++;
			pC++;
		}

		
		/* count processed bases */
		nRemainLen -= (nP2-nP1+1);
		if(nRemainLen < nMotifLen)
		{
			break;
		}

		/* get next */
		nP1 = nP2+1;
		if(nRemainLen < (GENOME_MED_CONTIG_LEN-nMotifLen+1))
		{
			nP2 = nP1+nRemainLen-nMotifLen;
		}
		else
		{
			nP2 = nP1+GENOME_MED_CONTIG_LEN-1;
		}
	}

	fprintf(fpOut, "\n===\n");
	fprintf(fpOut, "%d\t%d\n", nSiteNum, nEffecLen);

	/* close files */
	fclose(fpSeq);
	fclose(fpOut);
	fclose(fpCS);

	/* release memory */
	DestroyIntMatrix(pChrLen);
	SEQMTFDESTROY(pSeqMtf);

	/* return */
	return PROC_SUCCESS;
}

int MapcMyc_LOADSEQFROMGENOME_Conserved(struct SEQMOTIF *pSeqMtf, int nSeqId, 
							 FILE *fpSeq, FILE *fpCS, int nStart, int nEnd, char chStrand)
{
	/* define */
	unsigned char *pBase,*pC;
	int nLen;
	int nP1,nP2,nR1,nR2;
	int numread;
	int ni,nk;
	unsigned char bBase,bChar;

	/* init */
	if( (pSeqMtf == NULL) || (fpSeq == NULL) )
	{
		printf("Warning: MapcMyc_LOADSEQFROMGENOME, null seqmotif or seqfile pointer!\n");
		return PROC_FAILURE;
	}
	if( (nSeqId<0) || (nSeqId >= pSeqMtf->nSeqNum) )
	{
		printf("Warning: MapcMyc_LOADSEQFROMGENOME, sequence id out of range!\n");
		return PROC_FAILURE;
	}
	nLen = (nEnd-nStart+1);
	if( (nLen<0) || (nLen >  pSeqMtf->ppSeq[nSeqId]->nWidth ) )
	{
		printf("Warning: MapcMyc_LOADSEQFROMGENOME, the length of the sequence to be retrieved is larger than the length allowed by seqmotif!\n");
		return PROC_FAILURE;
	}
	if((chStrand != '+') && (chStrand != '-'))
	{
		printf("Warning: MapcMyc_LOADSEQFROMGENOME, strand information not specified!\n");
	}


	/* get cs */
	pC = pSeqMtf->ppSeq[1]->pMatElement;
	if( fseek( fpCS, nStart, SEEK_SET ) != 0 )
	{
		printf("Error: MapcMyc_LOADSEQFROMGENOME, cannot locate the required sequence!\n");
		exit(EXIT_FAILURE);
	}


	/* first base */
	numread = fread(pC, sizeof(unsigned char), (nEnd-nStart+1), fpCS );

	/* get seq */
	pBase = pSeqMtf->ppSeq[nSeqId]->pMatElement;
	nP1 = (int)(nStart/2);
	nR1 = nStart%2;
	nP2 = (int)(nEnd/2);
	nR2 = nEnd%2;
	nk = 0;

	if( fseek( fpSeq, nP1, SEEK_SET ) != 0 )
	{
		printf("Error: MapcMyc_LOADSEQFROMGENOME, cannot locate the required sequence!\n");
		exit(EXIT_FAILURE);
	}


	/* first base */
	numread = fread( &bChar, sizeof( unsigned char ), 1, fpSeq );
	if(nR1 == 0)
	{
		bBase = bChar >> 4;
		switch(bBase)
		{
			case 0: *pBase = 0;
				break;
			case 1: *pBase = 1;
				break;
			case 2: *pBase = 2;
				break;
			case 3: *pBase = 3;
				break;
			case 4: *pBase = 0;
				break;
			case 5: *pBase = 1;
				break;
			case 6: *pBase = 2;
				break;
			case 7: *pBase = 3;
				break;
			default: *pBase = 4;
		}
		pBase++;
		nk++;
	}

	bBase = bChar & 0x0F;
	switch(bBase)
	{
		case 0: *pBase = 0;
			break;
		case 1: *pBase = 1;
			break;
		case 2: *pBase = 2;
			break;
		case 3: *pBase = 3;
			break;
		case 4: *pBase = 0;
			break;
		case 5: *pBase = 1;
			break;
		case 6: *pBase = 2;
			break;
		case 7: *pBase = 3;
			break;
		default: *pBase = 4;
	}
	pBase++;
	nk++;

	/* middle bases */
	for(ni=(nP1+1); ni<nP2; ni++)
	{
		numread = fread( &bChar, sizeof( unsigned char ), 1, fpSeq );
		bBase = bChar >> 4;
		switch(bBase)
		{
			case 0: *pBase = 0;
				break;
			case 1: *pBase = 1;
				break;
			case 2: *pBase = 2;
				break;
			case 3: *pBase = 3;
				break;
			case 4: *pBase = 0;
				break;
			case 5: *pBase = 1;
				break;
			case 6: *pBase = 2;
				break;
			case 7: *pBase = 3;
				break;
			default: *pBase = 4;
		}
		pBase++;
		nk++;

		bBase = bChar & 0x0F;
		switch(bBase)
		{
			case 0: *pBase = 0;
				break;
			case 1: *pBase = 1;
				break;
			case 2: *pBase = 2;
				break;
			case 3: *pBase = 3;
				break;
			case 4: *pBase = 0;
				break;
			case 5: *pBase = 1;
				break;
			case 6: *pBase = 2;
				break;
			case 7: *pBase = 3;
				break;
			default: *pBase = 4;
		}
		pBase++;
		nk++;
	}
	
	/* last base */
	numread = fread( &bChar, sizeof( unsigned char ), 1, fpSeq );
	bBase = bChar >> 4;
	switch(bBase)
	{
		case 0: *pBase = 0;
			break;
		case 1: *pBase = 1;
			break;
		case 2: *pBase = 2;
			break;
		case 3: *pBase = 3;
			break;
		case 4: *pBase = 0;
			break;
		case 5: *pBase = 1;
			break;
		case 6: *pBase = 2;
			break;
		case 7: *pBase = 3;
			break;
		default: *pBase = 4;
	}
	pBase++;
	nk++;

	if(nR2 == 1)
	{
		bBase = bChar & 0x0F;
		switch(bBase)
		{
			case 0: *pBase = 0;
				break;
			case 1: *pBase = 1;
				break;
			case 2: *pBase = 2;
				break;
			case 3: *pBase = 3;
				break;
			case 4: *pBase = 0;
				break;
			case 5: *pBase = 1;
				break;
			case 6: *pBase = 2;
				break;
			case 7: *pBase = 3;
				break;
			default: *pBase = 4;
		}
		pBase++;
		nk++;
	}

	if(nk != nLen)
	{
		printf("Error: MapcMyc_LOADSEQFROMGENOME, sequence length not match!\n");
		exit(EXIT_FAILURE);
	}

	/* get reverse complement if necessary */
	if(chStrand == '-')
	{
		
	}

	/* return */
	return PROC_SUCCESS;
}

/* count cMyc in target region */
int Count_cMyc_In_ConservedTargetRegion_Main()
{
	/* define */
	char strGenomePath[LINE_LENGTH];
	char strCSPath[LINE_LENGTH];
	char strRegionFile[LINE_LENGTH];
	char strSiteFile[LINE_LENGTH];
	char strJobName[LINE_LENGTH];
	char strWorkPath[LINE_LENGTH];
	char strOutPath[LINE_LENGTH];
	int ni;
	int nChr = 21;
	int nMotifLen = 6;
	double dCutoff = CMYC_CONS_CUTOFF;
	/* double dCutoff = 18.0; */
	

	/* sprintf(strGenomePath, "C:\\Projects\\research_harvard\\genomelab_project\\genomes\\human\\b34_hg16\\");
	sprintf(strCSPath, "C:\\Projects\\research_harvard\\genomelab_project\\genomes\\crossspecies\\hg16Mm3Rn3\\hg_cs\\");
	sprintf(strSiteFile, "C:\\Projects\\research_harvard\\affy_project\\analysis\\tiling_paper\\cMycSites\\cMyc_chr%d.txt", nChr);
	strcpy(strWorkPath, "C:\\Projects\\research_harvard\\affy_project\\analysis\\tiling_paper\\");
	*/

	sprintf(strGenomePath, "/data/genomes/human/b33_hg15/");
	sprintf(strCSPath, "/data/genomes/human/b33_hg15/conservation/cs/");
	sprintf(strSiteFile, "cMyc_chr%d.txt", nChr);
	strcpy(strWorkPath, "./");
	
	/* naive */
	for(ni=1; ni<=3; ni++)
	{
		sprintf(strJobName, "cMyc_A_anova_sample2r2naive_%d_f_mean", ni);
		sprintf(strRegionFile, "%s%s.rstat", strWorkPath, strJobName);
		sprintf(strOutPath, "%s%s.site", strWorkPath, strJobName);
		Count_cMyc_In_ConservedTargetRegion(strGenomePath, strCSPath, nChr, nMotifLen,
			strSiteFile, strRegionFile, strOutPath, dCutoff);
	}
	for(ni=1; ni<=3; ni++)
	{
		sprintf(strJobName, "cMyc_A_anova_sample3r2naive_%d_f_mean", ni);
		sprintf(strRegionFile, "%s%s.rstat", strWorkPath, strJobName);
		sprintf(strOutPath, "%s%s.site", strWorkPath, strJobName);
		Count_cMyc_In_ConservedTargetRegion(strGenomePath, strCSPath, nChr, nMotifLen,
			strSiteFile, strRegionFile, strOutPath, dCutoff);
	}
	for(ni=1; ni<=3; ni++)
	{
		sprintf(strJobName, "cMyc_A_anova_sample3r4naive_%d_f_mean", ni);
		sprintf(strRegionFile, "%s%s.rstat", strWorkPath, strJobName);
		sprintf(strOutPath, "%s%s.site", strWorkPath, strJobName);
		Count_cMyc_In_ConservedTargetRegion(strGenomePath, strCSPath, nChr, nMotifLen,
			strSiteFile, strRegionFile, strOutPath, dCutoff);
	}

	sprintf(strJobName, "cMyc_A_anova_refnaive_f_mean");
	sprintf(strRegionFile, "%s%s.rstat", strWorkPath, strJobName);
	sprintf(strOutPath, "%s%s.site", strWorkPath, strJobName);
	Count_cMyc_In_ConservedTargetRegion(strGenomePath, strCSPath, nChr, nMotifLen,
		strSiteFile, strRegionFile, strOutPath, dCutoff);

	/* varsh */
	for(ni=1; ni<=3; ni++)
	{
		sprintf(strJobName, "cMyc_A_anova_sample2r2varsh_%d_f_mean", ni);
		sprintf(strRegionFile, "%s%s.rstat", strWorkPath, strJobName);
		sprintf(strOutPath, "%s%s.site", strWorkPath, strJobName);
		Count_cMyc_In_ConservedTargetRegion(strGenomePath, strCSPath, nChr, nMotifLen,
			strSiteFile, strRegionFile, strOutPath, dCutoff);
	}
	for(ni=1; ni<=3; ni++)
	{
		sprintf(strJobName, "cMyc_A_anova_sample3r2varsh_%d_f_mean", ni);
		sprintf(strRegionFile, "%s%s.rstat", strWorkPath, strJobName);
		sprintf(strOutPath, "%s%s.site", strWorkPath, strJobName);
		Count_cMyc_In_ConservedTargetRegion(strGenomePath, strCSPath, nChr, nMotifLen,
			strSiteFile, strRegionFile, strOutPath, dCutoff);
	}
	for(ni=1; ni<=3; ni++)
	{
		sprintf(strJobName, "cMyc_A_anova_sample3r4varsh_%d_f_mean", ni);
		sprintf(strRegionFile, "%s%s.rstat", strWorkPath, strJobName);
		sprintf(strOutPath, "%s%s.site", strWorkPath, strJobName);
		Count_cMyc_In_ConservedTargetRegion(strGenomePath, strCSPath, nChr, nMotifLen,
			strSiteFile, strRegionFile, strOutPath, dCutoff);
	}
	sprintf(strJobName, "cMyc_A_anova_refvarsh_f_mean");
	sprintf(strRegionFile, "%s%s.rstat", strWorkPath, strJobName);
	sprintf(strOutPath, "%s%s.site", strWorkPath, strJobName);
	Count_cMyc_In_ConservedTargetRegion(strGenomePath, strCSPath, nChr, nMotifLen,
		strSiteFile, strRegionFile, strOutPath, dCutoff);

	/* hmm */
	for(ni=1; ni<=3; ni++)
	{
		sprintf(strJobName, "cMyc_A_anova_sample2r2varsh_%d_f_refbind_pl", ni);
		sprintf(strRegionFile, "%s%s.rstat", strWorkPath, strJobName);
		sprintf(strOutPath, "%s%s.site", strWorkPath, strJobName);
		Count_cMyc_In_ConservedTargetRegion(strGenomePath, strCSPath, nChr, nMotifLen,
			strSiteFile, strRegionFile, strOutPath, dCutoff);
	}
	for(ni=1; ni<=3; ni++)
	{
		sprintf(strJobName, "cMyc_A_anova_sample3r2varsh_%d_f_refbind_pl", ni);
		sprintf(strRegionFile, "%s%s.rstat", strWorkPath, strJobName);
		sprintf(strOutPath, "%s%s.site", strWorkPath, strJobName);
		Count_cMyc_In_ConservedTargetRegion(strGenomePath, strCSPath, nChr, nMotifLen,
			strSiteFile, strRegionFile, strOutPath, dCutoff);
	}
	for(ni=1; ni<=3; ni++)
	{
		sprintf(strJobName, "cMyc_A_anova_sample3r4varsh_%d_f_refbind_pl", ni);
		sprintf(strRegionFile, "%s%s.rstat", strWorkPath, strJobName);
		sprintf(strOutPath, "%s%s.site", strWorkPath, strJobName);
		Count_cMyc_In_ConservedTargetRegion(strGenomePath, strCSPath, nChr, nMotifLen,
			strSiteFile, strRegionFile, strOutPath, dCutoff);
	}
	sprintf(strJobName, "cMyc_A_anova_refvarsh_f_refbind_pl");
	sprintf(strRegionFile, "%s%s.rstat", strWorkPath, strJobName);
	sprintf(strOutPath, "%s%s.site", strWorkPath, strJobName);
	Count_cMyc_In_ConservedTargetRegion(strGenomePath, strCSPath, nChr, nMotifLen,
		strSiteFile, strRegionFile, strOutPath, dCutoff);


	/* gtrans */
	for(ni=1; ni<=3; ni++)
	{
		sprintf(strJobName, "cMyc_A_mix_sample2_gtrans_%d_pvalue_f", ni);
		sprintf(strRegionFile, "%s%s.rstat", strWorkPath, strJobName);
		sprintf(strOutPath, "%s%s.site", strWorkPath, strJobName);
		Count_cMyc_In_ConservedTargetRegion(strGenomePath, strCSPath, nChr, nMotifLen,
			strSiteFile, strRegionFile, strOutPath, dCutoff);
	}
	for(ni=1; ni<=3; ni++)
	{
		sprintf(strJobName, "cMyc_A_mix_sample3r2_gtrans_%d_maxpvalue_f", ni);
		sprintf(strRegionFile, "%s%s.rstat", strWorkPath, strJobName);
		sprintf(strOutPath, "%s%s.site", strWorkPath, strJobName);
		Count_cMyc_In_ConservedTargetRegion(strGenomePath, strCSPath, nChr, nMotifLen,
			strSiteFile, strRegionFile, strOutPath, dCutoff);
	}
	for(ni=1; ni<=3; ni++)
	{
		sprintf(strJobName, "cMyc_A_mix_sample3r4_gtrans_%d_maxpvalue_f", ni);
		sprintf(strRegionFile, "%s%s.rstat", strWorkPath, strJobName);
		sprintf(strOutPath, "%s%s.site", strWorkPath, strJobName);
		Count_cMyc_In_ConservedTargetRegion(strGenomePath, strCSPath, nChr, nMotifLen,
			strSiteFile, strRegionFile, strOutPath, dCutoff);
	}
	sprintf(strJobName, "cMyc_A_mix_refmax_gtrans_pvalue_f");
	sprintf(strRegionFile, "%s%s.rstat", strWorkPath, strJobName);
	sprintf(strOutPath, "%s%s.site", strWorkPath, strJobName);
	Count_cMyc_In_ConservedTargetRegion(strGenomePath, strCSPath, nChr, nMotifLen,
		strSiteFile, strRegionFile, strOutPath, dCutoff);

	/* return */
	return PROC_SUCCESS;
}


int Count_cMyc_In_ConservedTargetRegion(char strGenomePath[], char strCSPath[],
							   int nChr, int nMotifLen, 
							   char strSiteFile[], char strRegionFile[], 
							   char strOutFile[], double dCutoff)
{
	/* define */
	FILE *fpIn;
	FILE *fpSeq;
	FILE *fpCS;
	FILE *fpOut;
	FILE *fpSite;

	char strSeqFile[LINE_LENGTH];
	char strCSFile[LINE_LENGTH];
	char strLine[MED_LINE_LENGTH];
	
	int nRegionNum = 0;
	struct DOUBLEMATRIX *pRegionStat;

	int nChrLen;
	struct INTMATRIX *pChrLen;
	int nEffecLen = 0;
	int nSiteNum = 0;
	struct SEQMOTIF *pSeqMtf;
	
	unsigned char *pBase,*pC;
	int nMisNum; 
	int nBadLen;
	char strSite[LINE_LENGTH];
	int nInitLen,nP1,nP2,nRemainLen;
	int ni,nj,nk;
	double dStart, dEnd, dMin, dMean;
	int nStart, nEnd;
	int nSiteChr,nSitePos;
	double dTemp;
	double dScore;

	/* init */
	pChrLen = NULL;
	sprintf(strSeqFile, "%schrlen.txt", strGenomePath);
	pChrLen = IMLOAD(strSeqFile);
	if(pChrLen == NULL)
	{
		printf("cannot load chromosome length\n");
		exit(EXIT_FAILURE);
	}

	sprintf(strSeqFile, "%schr%d.sq", strGenomePath, nChr);
	sprintf(strCSFile, "%schr%d.cs", strCSPath, nChr);
	nChrLen = pChrLen->pMatElement[nChr-1];

	/* get effec len region by region */
	fpIn = NULL;
	fpIn = fopen(strRegionFile, "r");
	if(fpIn == NULL)
	{
		printf("cannot load chromosome length\n");
		exit(EXIT_FAILURE);
	}
	while(fgets(strLine, MED_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		nRegionNum++;
	}
	fclose(fpIn);

	pRegionStat = NULL;
	pRegionStat = CreateDoubleMatrix(nRegionNum, 6);
	if(pRegionStat == NULL)
	{
		printf("cannot get region stat\n");
		exit(EXIT_FAILURE);
	}

	/* open sequence file and get the sequence */
	fpSeq = NULL;
	fpSeq = fopen(strSeqFile, "rb");
	if(fpSeq == NULL)
	{
		printf("Warning: cannot open sequence file!\n");
		exit(EXIT_FAILURE);
	}

	fpCS = NULL;
	fpCS = fopen(strCSFile, "rb");
	if(fpCS == NULL)
	{
		printf("Warning: cannot open conservation file!\n");
		exit(EXIT_FAILURE);
	}

	fpIn = NULL;
	fpIn = fopen(strRegionFile, "r");
	if(fpIn == NULL)
	{
		printf("cannot load chromosome length\n");
		exit(EXIT_FAILURE);
	}
	ni = 0;
	while(fgets(strLine, MED_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		
		sscanf(strLine, "%lf %lf %lf %lf", &dStart, &dEnd, &dMin, &dMean);
		
		nStart = (int)dStart;
		nEnd = (int)dEnd;
		nRemainLen = CMYC_WIN_SIZE-(nEnd-nStart+1);
		if(nRemainLen > 0)
		{
			nStart -= nRemainLen/2;
			nEnd += nRemainLen/2;
		}

		dStart = (double)nStart;
		dEnd = (double)nEnd;
		DMSETAT(pRegionStat, ni, 0, dStart);
		DMSETAT(pRegionStat, ni, 1, dEnd);
		DMSETAT(pRegionStat, ni, 2, dMin);
		DMSETAT(pRegionStat, ni, 3, dMean);


		pSeqMtf = NULL;
		pSeqMtf = SEQMTFCREATE(0, 2, 0, 0);
		if(pSeqMtf == NULL)
		{
			printf("Error: cannot create sequence/motif complex!\n");
			exit(EXIT_FAILURE);
		}
		/* sequence */
		SEQMTFCREATESEQ(pSeqMtf, 0, nEnd-nStart+nMotifLen);
		SEQMTFCREATESEQ(pSeqMtf, 1, nEnd-nStart+nMotifLen);

		MapcMyc_LOADSEQFROMGENOME_Conserved(pSeqMtf, 0, fpSeq, fpCS, nStart,
							nEnd+nMotifLen-1, '+');
			
		pBase = pSeqMtf->ppSeq[0]->pMatElement;
		pC = pSeqMtf->ppSeq[1]->pMatElement;

		nEffecLen = 0;
		for(nk=nStart; nk<=nEnd; nk++)
		{
			dScore = 0.0;
			nBadLen = 0;
			for(nj=0; nj<nMotifLen; nj++)
			{
				if(pBase[nj] >= 4)
				{
					nBadLen++;
					break;
				}
				dScore += pC[nj];
			}
			dScore /= (double)nMotifLen;

			/* if((nBadLen == 0) && (dScore >= dCutoff)) */
			if(nBadLen == 0)
			{
				nEffecLen++;
			}

			pBase++;
			pC++;
		}


		DMSETAT(pRegionStat, ni, 4, (double)(nEffecLen));

		SEQMTFDESTROY(pSeqMtf);

		ni++;
	}
	fclose(fpIn);
	fclose(fpSeq);
	fclose(fpCS);


	fpSite = NULL;
	fpSite = fopen(strSiteFile, "r");
	if(fpSite == NULL)
	{
		printf("cannot open file\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, MED_LINE_LENGTH, fpSite)!= NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strstr(strLine, "===") == strLine)
			break;

		sscanf(strLine, "%d %d", &nSiteChr, &nSitePos);

		for(ni=0; ni<nRegionNum; ni++)
		{
			nStart = (int)(DMGETAT(pRegionStat, ni, 0));
			nEnd = (int)(DMGETAT(pRegionStat, ni, 1));
			if( (nStart <= nSitePos) && (nSitePos <= nEnd) )
			{
				dTemp = DMGETAT(pRegionStat, ni, 5)+1.0;
				DMSETAT(pRegionStat, ni, 5, dTemp);
			}
		}
	}


	fclose(fpSite);

	/* write */
	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("cannot load chromosome length\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nRegionNum; ni++)
	{
		fprintf(fpOut, "%d %d % 9.7e % 9.7e %d %d\n",
			(int)(DMGETAT(pRegionStat, ni, 0)),
			(int)(DMGETAT(pRegionStat, ni, 1)),
			DMGETAT(pRegionStat, ni, 2),
			DMGETAT(pRegionStat, ni, 3),
			(int)(DMGETAT(pRegionStat, ni, 4)),
			(int)(DMGETAT(pRegionStat, ni, 5)));
	}

	fclose(fpOut);

	/* release memory */
	DestroyIntMatrix(pChrLen);
	DestroyDoubleMatrix(pRegionStat);

	/* return */
	return PROC_SUCCESS;
}

/* for the analysis of rich young's data */
int GetEBChIPCod()
{
	/* define */
	struct INTMATRIX *pRegions;
	char strRegionPath[LINE_LENGTH];
	char strRefGenePath[LINE_LENGTH];
	char strOutputPath[LINE_LENGTH];
	char strFileName[LINE_LENGTH];
	
	/* for handling refgene */
	char strSourceSpecies[LINE_LENGTH];
	FILE *fpRefGene;
	struct tagRefGene *pRefGene,*pCurrentRefGene;
	struct tagRefGene *pSourceRefGeneList;
	struct tagRefGene **vSourceRefGene;
	int nSourceRefNum;
	char strRefLine[LONG_LINE_LENGTH];

	/* for regions */
	int nRegionNum;
	struct tagRefGene **vTarRefGene;
	int nChr,nPos;
	int nDist;
	struct INTMATRIX *pMinDist;

	/* negative controls */
	struct DOUBLEMATRIX *pRand;
	struct DOUBLEMATRIX *pSortRand;
	struct LONGMATRIX *pSortIndex;
	int nP1, nP2;


	/* out */
	FILE *fpOut;
	int ni,nj;

	/* init */
	strcpy(strSourceSpecies, "human");
	strcpy(strRegionPath, "C:\\Projects\\research_harvard\\stemcell_project\\ChIP-chip\\RichYoungChIP\\YoungLab-ChIP_region_all\\Nanog_region.txt");
	/* strcpy(strRefGenePath, "..\\genomes\\human\\b35_hg17\\refGene_sorted.txt"); */
	strcpy(strRefGenePath, "C:\\Data\\genomes\\human\\b35_hg17\\annotation\\ensGene_sorted.txt");
	strcpy(strOutputPath, "C:\\Projects\\research_harvard\\stemcell_project\\ChIP-chip\\RichYoungChIP\\YoungLab-ChIP_region_all\\Nanog_ens_cod");

	/* load source refgene */
	fpRefGene = NULL;
	fpRefGene = fopen(strRefGenePath, "r");
	if(fpRefGene == NULL)
	{
		printf("Error: RefGene_GetOrtholog, cannot open source refgene file!\n");
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
	fclose(fpRefGene);

	vSourceRefGene = NULL;
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

	fclose(fpRefGene);


	/* load regions */
	pRegions = NULL;
	pRegions = IMLOAD(strRegionPath);
	if(pRegions == NULL)
	{
		printf("Error: RefGene_GetOrtholog, cannot load regions!\n");
		exit(EXIT_FAILURE);
	}

	nRegionNum = pRegions->nHeight;
	vTarRefGene = NULL;
	vTarRefGene = (struct tagRefGene **)calloc(nRegionNum, sizeof(struct tagRefGene*));
	if(vTarRefGene == NULL)
	{
		printf("Error: RefGene_GetOrtholog, organize refgene into array!\n");
		exit(EXIT_FAILURE);
	}

	pMinDist = NULL;
	pMinDist = CreateIntMatrix(nRegionNum, 1);
	if(pMinDist == NULL)
	{
		printf("Error: RefGene_GetOrtholog, organize refgene into array!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nRegionNum; ni++)
	{
		nChr = IMGETAT(pRegions, ni, 0);
		nPos = (IMGETAT(pRegions, ni, 1)+IMGETAT(pRegions, ni, 2))/2;
		pMinDist->pMatElement[ni] = 1000000;

		for(nj=0; nj<nSourceRefNum; nj++)
		{
			if(vSourceRefGene[nj]->nChrom != nChr)
				continue;

			if(vSourceRefGene[nj]->chStrand == '+')
			{
				nDist = vSourceRefGene[nj]->nTxStart-nPos;
			}
			else
			{
				nDist = nPos-vSourceRefGene[nj]->nTxEnd;
			}

			if( (nDist > 8500) || (nDist < -2500) )
				continue;

			if( fabs(nDist) < fabs(pMinDist->pMatElement[ni]) )
			{
				pMinDist->pMatElement[ni] = nDist;
				vTarRefGene[ni] = vSourceRefGene[nj];
			}
		}
	}

	/* export */
	sprintf(strFileName, "%s.txt", strOutputPath);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: cannot open output files!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nRegionNum; ni++)
	{
		if(vTarRefGene[ni] != NULL)
		{
			fprintf(fpOut, "%s\t%s\t%d\t%d\t%c\t%d\n",
				vTarRefGene[ni]->strName, vTarRefGene[ni]->strChrom, 
				IMGETAT(pRegions, ni, 1)-200, IMGETAT(pRegions, ni, 2)+200,
				vTarRefGene[ni]->chStrand, pMinDist->pMatElement[ni]);
		}
	}

	fclose(fpOut);

	/* get random negative control */
	sprintf(strFileName, "%s.neg", strOutputPath);
	fpOut = NULL;
	fpOut = fopen(strFileName, "w");
	if(fpOut == NULL)
	{
		printf("Error: cannot open output files!\n");
		exit(EXIT_FAILURE);
	}
	pRand = DMRANDU(1, nSourceRefNum);
	DMSORTMERGEA_0(pRand, &pSortRand, &pSortIndex);
	
	for(ni=0; ni<nRegionNum; ni++)
	{
		nj = pSortIndex->pMatElement[ni];

		if(vTarRefGene[ni] != NULL)
		{
			if(vSourceRefGene[nj]->chStrand == '+')
			{
				/* nP1 = vSourceRefGene[nj]->nTxStart-pMinDist->pMatElement[ni];
				nP2 = nP1+500;
				nP1 = nP1-500; */
				nP1 = vSourceRefGene[nj]->nTxStart-8000;
				nP2 = vSourceRefGene[nj]->nTxStart+2000;

			}
			else
			{
				/* nP1 = vSourceRefGene[nj]->nTxEnd+pMinDist->pMatElement[ni];
				nP2 = nP1+500;
				nP1 = nP1-500; */
				nP1 = vSourceRefGene[nj]->nTxEnd-2000;
				nP2 = vSourceRefGene[nj]->nTxEnd+8000;

			}

			if(nP1 < 0)
				nP1 = 0;
			if(nP2 < 0)
				nP2 = 0;

			fprintf(fpOut, "%s\t%s\t%d\t%d\t%c\t%d\n",
					vSourceRefGene[nj]->strName, vSourceRefGene[nj]->strChrom, 
					nP1, nP2,
					vSourceRefGene[nj]->chStrand, pMinDist->pMatElement[ni]);
		}
	}
	fclose(fpOut);


	/* release memory */
	DestroyDoubleMatrix(pRand);
	DestroyDoubleMatrix(pSortRand);
	DestroyLongMatrix(pSortIndex);

	DestroyIntMatrix(pRegions);
	DestroyIntMatrix(pMinDist);
	
	for(ni=0; ni<nRegionNum; ni++)
	{
		vTarRefGene[ni] = NULL;
	}
	free(vTarRefGene);

	for(ni=0; ni<nSourceRefNum; ni++)
	{
		RefGeneDestroy(vSourceRefGene[ni]);
		vSourceRefGene[ni] = NULL;
	}
	free(vSourceRefGene);

	/* return */
	return PROC_SUCCESS;
}

int ConnectSitewithOrtho(int argv, char **argc)
{
	/* environment */
	int nSpeciesNum;
	int nGeneNum;

	char vSpeciesName[3][LINE_LENGTH];
	struct tagRefGene **vRefGene[3];

	struct INTMATRIX **vStat;
	int ni,nj;
	char *chp1,*chp2;

	/* file */
	char strLine[LONG_LINE_LENGTH];
	char strSiteFile[LINE_LENGTH];
	char strMapFile[LINE_LENGTH];
	char strOutFile[LINE_LENGTH];
	char strRefId[LINE_LENGTH];
	FILE *fpIn;
	FILE *fpOut1;
	FILE *fpOut2;
	char strFileName[LINE_LENGTH];
	int nTemp;

	/* sites */
	char strAlias[LINE_LENGTH];
	char strChr[LINE_LENGTH];
	int nStart,nEnd;
	char chStrand;
	int nRelPos;
	int nOK;

	/* init */
	/* strcpy(strSiteFile, argc[1]);
	strcpy(strMapFile, argc[2]);
	strcpy(strOutFile, argc[3]); */

	strcpy(strSiteFile, "C:\\Projects\\research_harvard\\embroynicstemcell_project\\RichYoungChIP\\YoungLab-ChIP_region_all\\Oct4-Sox2_cod.neg");
	strcpy(strMapFile, "C:\\Projects\\research_harvard\\embroynicstemcell_project\\RichYoungChIP\\YoungLab-ChIP_region_all\\Oct4-Sox2_neg_ortho.msomap");
	strcpy(strOutFile, "C:\\Projects\\research_harvard\\embroynicstemcell_project\\RichYoungChIP\\YoungLab-ChIP_region_all\\Oct4-Sox2_neg");

	nSpeciesNum = 3;
	nGeneNum = 0;

	/* get gene number */
	fpIn = NULL;
	fpIn = fopen(strMapFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: cannot open input file\n");
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strLine[0] == '>')
			nGeneNum++;
	}

	fclose(fpIn);

	/* load ortholog */
	for(ni=0; ni<nSpeciesNum; ni++)
	{
		vRefGene[ni] = NULL;
		vRefGene[ni] = (struct tagRefGene **)calloc(nGeneNum, sizeof(struct tagRefGene *));
	}
	vStat = NULL;
	vStat = (struct INTMATRIX **)calloc(nGeneNum, sizeof(struct INTMATRIX *));

	fpIn = NULL;
	fpIn = fopen(strMapFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: cannot open input file\n");
	}

	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(strLine[0] == '>')
		{
			vStat[ni] = CreateIntMatrix(nSpeciesNum, 2);

			for(nj=0; nj<nSpeciesNum; nj++)
			{
				fgets(strLine, LONG_LINE_LENGTH, fpIn);
				StrTrimLeft(strLine);
				StrTrimRight(strLine);

				chp2 = strchr(strLine, '\t');
				*chp2 = '\0';
				strcpy(vSpeciesName[nj], strLine);

				chp1 = chp2+1;
				chp2 = strchr(chp1, '\t');
				*chp2 = '\0';
				nTemp = atoi(chp1);
				IMSETAT(vStat[ni], nj, 0, nTemp);
				chp1 = chp2+1;
				chp2 = strchr(chp1, '\t');
				*chp2 = '\0';
				nTemp = atoi(chp1);
				IMSETAT(vStat[ni], nj, 1, nTemp);
				chp1 = chp2+1;
				sscanf(chp1, "%s", strRefId);
				
				if(strcmp(strRefId, "---") == 0)
				{
					continue;
				}
				
				vRefGene[nj][ni] = RefGeneCreate();
				RefGeneInit_FromGenomeLabFormat(vRefGene[nj][ni], chp1, vSpeciesName[nj]);
			}
			ni++;
		}
	}

	fclose(fpIn);
	if(ni != nGeneNum)
	{
		printf("Error: cannot open input file\n");
	}

	/* load regions */
	fpIn = NULL;
	fpIn = fopen(strSiteFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: cannot open input file\n");
	}

	sprintf(strFileName, "%s.site", strOutFile);
	fpOut1 = NULL;
	fpOut1 = fopen(strFileName, "w");
	if(fpOut1 == NULL)
	{
		printf("Error: cannot open input file\n");
	}

	sprintf(strFileName, "%s.ortho", strOutFile);
	fpOut2 = NULL;
	fpOut2 = fopen(strFileName, "w");
	if(fpOut2 == NULL)
	{
		printf("Error: cannot open input file\n");
	}

	
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%s %s %d %d %c %d", 
			strAlias, strChr, &nStart, &nEnd, &chStrand, &nRelPos);
		nOK = 0;

		for(ni=0; ni<nGeneNum; ni++)
		{
			if(strcmp(strAlias, vRefGene[0][ni]->strName) == 0)
			{
				nOK = 1;
				break;
			}
		}

		if(nOK == 1)
		{
			fprintf(fpOut1, "%s_%s_%s\t%s\t%d\t%d\t%c\t%d\n", 
				strAlias, vSpeciesName[0], vRefGene[0][ni]->strName,
				strChr, nStart, nEnd, chStrand, nRelPos);
			fprintf(fpOut2, ">%s\n", strAlias);
			for(nj=0; nj<nSpeciesNum; nj++)
			{
				if(vRefGene[nj][ni] == NULL)
				{
					fprintf(fpOut2, "%s\t0\t0\t---\n", vSpeciesName[nj]);
				}
				else
				{
					fprintf(fpOut2, "%s\t%d\t%d\t", vSpeciesName[nj],
						IMGETAT(vStat[ni], nj, 0), IMGETAT(vStat[ni], nj, 1));
					RefGeneWrite(vRefGene[nj][ni], fpOut2);
				}
			}
			fprintf(fpOut2, "\n");
		}
	}


	fclose(fpIn);
	fclose(fpOut1);
	fclose(fpOut2);


	/* release memory */
	for(ni=0; ni<nSpeciesNum; ni++)
	{
		for(nj=0; nj<nGeneNum; nj++)
		{
			RefGeneDestroy(vRefGene[ni][nj]);
			vRefGene[ni][nj] = NULL;
		}
		free(vRefGene[ni]);
	}

	for(ni=0; ni<nGeneNum; ni++)
	{
		DestroyIntMatrix(vStat[ni]);
		vStat[ni] = NULL;
	}
	free(vStat);


	/* return */
	return PROC_SUCCESS;
}


int HG17GetTSS_Main(int argv, char **argc)
{
	FILE *fpIn;
	FILE *fpOut;
	struct tagRefGene **vRefGeneDatabase;
	int nRefNum = 0;
	char strSpecies[LINE_LENGTH];
	char strGenomePath[LINE_LENGTH];
	char strDatabasePath[LINE_LENGTH];
	char strFileName[LINE_LENGTH];
	char strOutFile[LINE_LENGTH];
	struct INTMATRIX *pChrLen;
	int ni, nStart, nEnd, nChr;

	/* init */
	strcpy(strOutFile, argc[1]);
	strcpy(strSpecies, "human");
	strcpy(strGenomePath, "/data/genomes/human/b35_hg17/");
	strcpy(strDatabasePath, "/data/genomes/human/b35_hg17/annotation/refFlat_sorted.txt");

	sprintf(strFileName, "%schrlen.txt", strGenomePath);
	pChrLen = NULL;
	pChrLen = IMLOAD(strFileName);
	if(pChrLen == NULL)
	{
		printf("Error!\n");
		exit(EXIT_FAILURE);
	}

	/* open file */
	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error!\n");
		exit(EXIT_FAILURE);
	}

	vRefGeneDatabase = NULL;
	vRefGeneDatabase = RefGene_LoadDatabase(strDatabasePath, 1, strSpecies, &nRefNum);

	/* write */
	for(ni=0; ni<nRefNum; ni++)
	{
		nChr = vRefGeneDatabase[ni]->nChrom;
		nStart = vRefGeneDatabase[ni]->nTxStart-50000;
		nEnd = vRefGeneDatabase[ni]->nTxEnd+50000;

		if(nStart < 0)
			nStart = 0;
		if(nEnd >= pChrLen->pMatElement[nChr-1])
			nEnd = pChrLen->pMatElement[nChr-1]-1;

		fprintf(fpOut, "%d\t%d\t%d\n", nChr, nStart, nEnd);
	}

	/* close files */
	RefGene_ClearDatabase(&vRefGeneDatabase, nRefNum);
	fclose(fpOut);
	DestroyIntMatrix(pChrLen);

	/* return */
	return PROC_SUCCESS;
}

int MM6GetGeneCover_Main(int argv, char **argc)
{
	return 1;
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

/* ----------------------------------------------------------------------- */ 
/*  EEL_GetReverseComplement_Main.                                         */
/* ----------------------------------------------------------------------- */ 
int EEL_GetReverseComplement_Main()
{
	/* define */
	struct tagSequence *pSeq = NULL;
	char strFileName[LINE_LENGTH];
	char strWorkPath[LINE_LENGTH];
	char strOutFile[LINE_LENGTH];
	int nClustNum = 249;
	int ni;
	char strSpecies[LINE_LENGTH];
	FILE *fpOut;

	/* init */
	strcpy(strSpecies, "zebrafish");
	strcpy(strWorkPath, "C:\\Projects\\research_harvard\\hedgehog_project\\ChIP-chip\\Agilent-200508-Gli1EB\\analysis_paper_eel\\seq\\");

	/* process one by one */
	for(ni=0; ni<nClustNum; ni++)
	{	
		sprintf(strFileName, "%s%d_%s.fa", strWorkPath, ni, strSpecies);
		pSeq = NULL;
		LoadFullSequenceList(strFileName, &pSeq);

		sprintf(strOutFile, "%s%d_%s_rev.fa", strWorkPath, ni, strSpecies);
		fpOut = NULL;
		fpOut = fopen(strOutFile, "w");
		if(fpOut == NULL)
		{
			printf("Error: cannot open output file!\n");
			exit(EXIT_FAILURE);
		}

		if(pSeq != NULL)
		{
			SequenceWriteToFasta_ByStrand(pSeq, fpOut, '-', 1);
		}

		fclose(fpOut);
		SequenceListClear(&pSeq);
	}

	/* return */
	return PROC_SUCCESS;
}

int GliArrayGenerateRand()
{
	char strInFile[LINE_LENGTH];
	char strOutFile[LINE_LENGTH];
	FILE *fpIn;
	FILE *fpOut;
	double dRand;

	char strLine[LONG_LINE_LENGTH];

	strcpy(strInFile, "C:\\Projects\\research_harvard\\hedgehog_project\\ChIP-chip\\Agilent-200604-Promoter\\tier8_r250rand.txt");
	strcpy(strOutFile, "C:\\Projects\\research_harvard\\hedgehog_project\\ChIP-chip\\Agilent-200604-Promoter\\tier8_r250rand_r.txt");

	fpIn = NULL;
	fpIn = fopen(strInFile, "r");
	if(fpIn == NULL)
	{
		printf("Error\n");
		exit(EXIT_FAILURE);
	}
	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!= NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		dRand = rand_u();
		fprintf(fpOut, "%s\t%f\n", strLine, dRand);
	}

	fclose(fpIn);
	fclose(fpOut);


	return PROC_SUCCESS;
}

/* Gene set analysis */
int menu_geneset(int argv, char **argc)
{
	/* define */
	char strGeneSetFile[MED_LINE_LENGTH];
	char strGeneLabFile[MED_LINE_LENGTH];
	char strOutFile[MED_LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	struct tagString ** vGeneName;
	struct tagString ** vRefName;
	struct BYTEMATRIX *vLabel;
	char strGeneName[LINE_LENGTH];
	char strRefGeneName[LINE_LENGTH];
	int nLabel;
	int nGeneNum = 0;
	int nTotNum,nHitNum;
	char *chp,*chp2;

	FILE *fpIn;
	FILE *fpOut;
	int ni;

	/* init */
	if(argv < 4)
	{
		printf("parmater error!\n");
		exit(EXIT_FAILURE);
	}

	strcpy(strGeneLabFile, argc[1]);
	strcpy(strGeneSetFile, argc[2]);
	strcpy(strOutFile, argc[3]);

	/* load label */
	fpIn = NULL;
	fpIn = fopen(strGeneLabFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: cannot load gene labels!\n");
		exit(EXIT_FAILURE);
	}

	fgets(strLine, LONG_LINE_LENGTH, fpIn);
	sscanf(strLine, "%d", &nGeneNum);
	vGeneName = NULL;
	vRefName = NULL;
	vLabel = NULL;
	vGeneName = (struct tagString **)calloc(nGeneNum, sizeof(struct tagString *));
	vRefName = (struct tagString **)calloc(nGeneNum, sizeof(struct tagString *));
	vLabel = CreateByteMatrix(1, nGeneNum);
	if( (vGeneName == NULL) || (vRefName == NULL) || (vLabel == NULL) ) 
	{
		printf("Error: cannot create memory!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(ni>=nGeneNum)
		{
			printf("Error: gene number not correct!\n");
			exit(EXIT_FAILURE);
		}

		sscanf(strLine, "%s %d %s", strGeneName, &nLabel, strRefGeneName);
		StringAddTail(vGeneName+ni, strGeneName);
		StringAddTail(vRefName+ni, strRefGeneName);
		vLabel->pMatElement[ni] = (unsigned char)nLabel;

		ni++;
	}

	if(ni!=nGeneNum)
	{
		printf("Error: gene number not correct!\n");
		exit(EXIT_FAILURE);
	}

	fclose(fpIn);

	/* create overlap mapping */
	fpIn = NULL;
	fpIn = fopen(strGeneSetFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: cannot load gene sets!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: cannot open output files!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		nTotNum = 0;
		nHitNum = 0;

		chp = strchr(strLine, '\t');
		chp2 = chp+1;

		chp = strchr(chp2, '\t');
		if(chp != NULL)
		{
			*chp = '\0';
			chp2 = chp+1;
			chp = strchr(chp2, '\t');
		}
		else
		{
			chp2 = NULL;
		}

		fprintf(fpOut, "%s\t", strLine);

		while(chp != NULL)
		{
			nTotNum += 1;
			*chp = '\0';
			strcpy(strGeneName, chp2);
			chp2 = chp+1;

			StrTrimRight(strGeneName);
			StrTrimLeft(strGeneName);
			StrMakeUpper(strGeneName);

			for(ni=0; ni<nGeneNum; ni++)
			{
				if(strcmp(strGeneName, vGeneName[ni]->m_pString) == 0)
				{
					nHitNum += (int)(vLabel->pMatElement[ni]);
					break;
				}
			}

			chp = strchr(chp2, '\t');
		}

		if(chp2 != NULL)
		{
			nTotNum += 1;
			strcpy(strGeneName, chp2);
			StrTrimRight(strGeneName);
			StrTrimLeft(strGeneName);
			StrMakeUpper(strGeneName);

			for(ni=0; ni<nGeneNum; ni++)
			{
				if(strcmp(strGeneName, vGeneName[ni]->m_pString) == 0)
				{
					nHitNum += (int)(vLabel->pMatElement[ni]);
					break;
				}
			}
		}

		fprintf(fpOut, "%d\t%d\n", nTotNum, nHitNum);
	}


	fclose(fpIn);
	fclose(fpOut);

	/* release memory */
	DestroyByteMatrix(vLabel);
	for(ni=0; ni<nGeneNum; ni++)
	{
		DeleteString(vGeneName[ni]);
		vGeneName[ni] = NULL;
		DeleteString(vRefName[ni]);
		vRefName[ni] = NULL;
	}
	free(vGeneName);
	free(vRefName);

	return PROC_SUCCESS;
}

int MapHomoloGene()
{
	/* define */
	int nDatabaseNum = 174762;
	struct INTMATRIX *pHomoloID = NULL;
	struct INTMATRIX *pSpecies = NULL;
	struct INTMATRIX *pEntrezID = NULL;
	struct tagString **vGeneName = NULL;
	int nSource,nDest;
	int ni,nj,nk,nOK;
	FILE *fpIn;
	FILE *fpOut;
	char strDatabasePath[LINE_LENGTH];
	char strInputFile[LINE_LENGTH];
	char strOutputFile[LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	int nHID,nSID,nGID;
	char strGeneName[LINE_LENGTH];

	/* init species */
	nSource = 7227;
	nDest = 10090;

	/* init data */
	strcpy(strDatabasePath, "C:\\Data\\homologene\\build48.1\\homologene.data");
	strcpy(strInputFile, "C:\\Projects\\research_harvard\\hedgehog_project\\RNAi\\fly_gene_2.txt");
	strcpy(strOutputFile, "C:\\Projects\\research_harvard\\hedgehog_project\\RNAi\\fly-mouse_GeneID.txt");

	/* load database */
	pHomoloID = CreateIntMatrix(1,nDatabaseNum);
	pSpecies = CreateIntMatrix(1,nDatabaseNum);
	pEntrezID = CreateIntMatrix(1,nDatabaseNum);
	vGeneName = (struct tagString **)calloc(nDatabaseNum, sizeof(struct tagString*));

	if( (pHomoloID == NULL) || (pSpecies == NULL) ||
		(pEntrezID == NULL) || (vGeneName == NULL))
	{
		printf("error: cannot create memory!\n");
		exit(EXIT_FAILURE);
	}

	fpIn = NULL;
	fpIn = fopen(strDatabasePath, "r");
	if(fpIn == NULL)
	{
		printf("error: cannot open file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%d %d %d %s", &nHID, &nSID, &nGID, strGeneName); 
		pHomoloID->pMatElement[ni] = nHID;
		pSpecies->pMatElement[ni] = nSID;
		pEntrezID->pMatElement[ni] = nGID;
		StrMakeUpper(strGeneName);
		StringAddTail(vGeneName+ni, strGeneName);
		ni++;
	}

	fclose(fpIn);

	if(ni != nDatabaseNum)
	{
		printf("error: database size not match!\n");
		exit(EXIT_FAILURE);
	}

	/* map homolo gene */
	fpIn = NULL;
	fpIn = fopen(strInputFile, "r");
	if(fpIn == NULL)
	{
		printf("error: cannot open file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutputFile, "w");
	if(fpOut == NULL)
	{
		printf("error: cannot open file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%s", strGeneName); 
		StrMakeUpper(strGeneName);
		for(ni=0; ni<nDatabaseNum; ni++)
		{
			if( (pSpecies->pMatElement[ni] == nSource) && (strcmp(strGeneName, vGeneName[ni]->m_pString) == 0) )
			{
				break;
			}
		}

		/* if found */
		if(ni<nDatabaseNum)
		{
			nOK = 0;
			for(nj=ni-1; nj>=0; nj--)
			{
				if(pHomoloID->pMatElement[nj] != pHomoloID->pMatElement[ni])
				{
					break;
				}

				if(pSpecies->pMatElement[nj] == nDest)
				{
					nOK = 1;
					fprintf(fpOut, "%d\t%d\t%d\t%s\t%s\n", pEntrezID->pMatElement[nj],
						pEntrezID->pMatElement[ni], pHomoloID->pMatElement[ni],
						vGeneName[nj]->m_pString, vGeneName[ni]->m_pString);
					break;
				}
			}

			if(nOK == 0)
			{
				for(nj=ni+1; nj<nDatabaseNum; nj++)
				{
					if(pHomoloID->pMatElement[nj] != pHomoloID->pMatElement[ni])
					{
						break;
					}

					if(pSpecies->pMatElement[nj] == nDest)
					{
						nOK = 1;
						fprintf(fpOut, "%d\t%d\t%d\t%s\t%s\n", pEntrezID->pMatElement[nj],
							pEntrezID->pMatElement[ni], pHomoloID->pMatElement[ni],
							vGeneName[nj]->m_pString, vGeneName[ni]->m_pString);
						break;
					}
				}
			}
		}

	}


	fclose(fpIn);
	fclose(fpOut);


	/* release memory */
	DestroyIntMatrix(pHomoloID);
	DestroyIntMatrix(pSpecies);
	DestroyIntMatrix(pEntrezID);
	for(ni=0; ni<nDatabaseNum; ni++)
	{
		DeleteString(vGeneName[ni]);
		vGeneName[ni] = NULL;
	}
	free(vGeneName);

	return PROC_SUCCESS;
}

int MapHomoloGene_ByHID()
{
	/* define */
	int nDatabaseNum = 174762;
	struct INTMATRIX *pHomoloID = NULL;
	struct INTMATRIX *pSpecies = NULL;
	struct INTMATRIX *pEntrezID = NULL;
	struct tagString **vGeneName = NULL;
	int nSource,nDest;
	int ni,nj,nk,nOK;
	FILE *fpIn;
	FILE *fpOut;
	char strDatabasePath[LINE_LENGTH];
	char strInputFile[LINE_LENGTH];
	char strOutputFile[LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	int nHID,nSID,nGID;
	char strGeneName[LINE_LENGTH];

	/* init species */
	nSource = 7227;
	nDest = 10090;

	/* init data */
	strcpy(strDatabasePath, "C:\\Data\\homologene\\build48.1\\homologene.data");
	strcpy(strInputFile, "C:\\Projects\\research_harvard\\hedgehog_project\\RNAi\\fly_homologene.txt");
	strcpy(strOutputFile, "C:\\Projects\\research_harvard\\hedgehog_project\\RNAi\\fly-mouse_GeneID.txt");

	/* load database */
	pHomoloID = CreateIntMatrix(1,nDatabaseNum);
	pSpecies = CreateIntMatrix(1,nDatabaseNum);
	pEntrezID = CreateIntMatrix(1,nDatabaseNum);
	vGeneName = (struct tagString **)calloc(nDatabaseNum, sizeof(struct tagString*));

	if( (pHomoloID == NULL) || (pSpecies == NULL) ||
		(pEntrezID == NULL) || (vGeneName == NULL))
	{
		printf("error: cannot create memory!\n");
		exit(EXIT_FAILURE);
	}

	fpIn = NULL;
	fpIn = fopen(strDatabasePath, "r");
	if(fpIn == NULL)
	{
		printf("error: cannot open file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%d %d %d %s", &nHID, &nSID, &nGID, strGeneName); 
		pHomoloID->pMatElement[ni] = nHID;
		pSpecies->pMatElement[ni] = nSID;
		pEntrezID->pMatElement[ni] = nGID;
		StrMakeUpper(strGeneName);
		StringAddTail(vGeneName+ni, strGeneName);
		ni++;
	}

	fclose(fpIn);

	if(ni != nDatabaseNum)
	{
		printf("error: database size not match!\n");
		exit(EXIT_FAILURE);
	}

	/* map homolo gene */
	fpIn = NULL;
	fpIn = fopen(strInputFile, "r");
	if(fpIn == NULL)
	{
		printf("error: cannot open file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutputFile, "w");
	if(fpOut == NULL)
	{
		printf("error: cannot open file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%d", &nk); 
		for(ni=0; ni<nDatabaseNum; ni++)
		{
			if( (pHomoloID->pMatElement[ni] == nk) )
			{
				break;
			}
		}

		/* if found */
		if(ni<nDatabaseNum)
		{
			nOK = 0;
			for(nj=ni-1; nj>=0; nj--)
			{
				if(pHomoloID->pMatElement[nj] != pHomoloID->pMatElement[ni])
				{
					break;
				}

				if(pSpecies->pMatElement[nj] == nDest)
				{
					nOK = 1;
					fprintf(fpOut, "%d\t%d\t%d\t%s\n", pEntrezID->pMatElement[nj],
						pSpecies->pMatElement[nj], pHomoloID->pMatElement[ni],
						vGeneName[nj]->m_pString);
					break;
				}
			}

			if(nOK == 0)
			{
				for(nj=ni+1; nj<nDatabaseNum; nj++)
				{
					if(pHomoloID->pMatElement[nj] != pHomoloID->pMatElement[ni])
					{
						break;
					}

					if(pSpecies->pMatElement[nj] == nDest)
					{
						nOK = 1;
						fprintf(fpOut, "%d\t%d\t%d\t%s\n", pEntrezID->pMatElement[nj],
							pSpecies->pMatElement[nj], pHomoloID->pMatElement[ni],
							vGeneName[nj]->m_pString);
						break;
					}
				}
			}

			nOK = 0;
			for(nj=ni-1; nj>=0; nj--)
			{
				if(pHomoloID->pMatElement[nj] != pHomoloID->pMatElement[ni])
				{
					break;
				}

				if(pSpecies->pMatElement[nj] == nSource)
				{
					nOK = 1;
					fprintf(fpOut, "%d\t%d\t%d\t%s\n", pEntrezID->pMatElement[nj],
						pSpecies->pMatElement[nj], pHomoloID->pMatElement[ni],
						vGeneName[nj]->m_pString);
					break;
				}
			}

			if(nOK == 0)
			{
				for(nj=ni+1; nj<nDatabaseNum; nj++)
				{
					if(pHomoloID->pMatElement[nj] != pHomoloID->pMatElement[ni])
					{
						break;
					}

					if(pSpecies->pMatElement[nj] == nSource)
					{
						nOK = 1;
						fprintf(fpOut, "%d\t%d\t%d\t%s\n", pEntrezID->pMatElement[nj],
							pSpecies->pMatElement[nj], pHomoloID->pMatElement[ni],
							vGeneName[nj]->m_pString);
						break;
					}
				}
			}
		}

	}


	fclose(fpIn);
	fclose(fpOut);


	/* release memory */
	DestroyIntMatrix(pHomoloID);
	DestroyIntMatrix(pSpecies);
	DestroyIntMatrix(pEntrezID);
	for(ni=0; ni<nDatabaseNum; ni++)
	{
		DeleteString(vGeneName[ni]);
		vGeneName[ni] = NULL;
	}
	free(vGeneName);

	return PROC_SUCCESS;
}

int MapGeneID2RefGene()
{
	/* define */
	int nDatabaseNum = 50751;
	struct INTMATRIX *pSpecies = NULL;
	struct INTMATRIX *pEntrezID = NULL;
	struct tagString **vGeneName = NULL;
	int nSource,nDest;
	int ni,nj,nk,nOK;
	FILE *fpIn;
	FILE *fpOut;
	char strDatabasePath[LINE_LENGTH];
	char strInputFile[LINE_LENGTH];
	char strOutputFile[LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	int nHID,nSID,nGID;
	char strGeneName[LINE_LENGTH];
	char *chp;

	/* init species */
	nSource = 10090;

	/* init data */
	strcpy(strDatabasePath, "C:\\Projects\\research_harvard\\hedgehog_project\\RNAi\\mouse_geneid2refseq.txt");
	strcpy(strInputFile, "C:\\Projects\\research_harvard\\hedgehog_project\\RNAi\\fly-mouse_GeneID.txt");
	strcpy(strOutputFile, "C:\\Projects\\research_harvard\\hedgehog_project\\RNAi\\RNAi_mouse_refid.txt");

	/* load database */
	pSpecies = CreateIntMatrix(1,nDatabaseNum);
	pEntrezID = CreateIntMatrix(1,nDatabaseNum);
	vGeneName = (struct tagString **)calloc(nDatabaseNum, sizeof(struct tagString*));

	if( (pSpecies == NULL) || (pEntrezID == NULL) || (vGeneName == NULL))
	{
		printf("error: cannot create memory!\n");
		exit(EXIT_FAILURE);
	}

	fpIn = NULL;
	fpIn = fopen(strDatabasePath, "r");
	if(fpIn == NULL)
	{
		printf("error: cannot open file!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%d %d %s", &nSID, &nGID, strGeneName); 
		pSpecies->pMatElement[ni] = nSID;
		pEntrezID->pMatElement[ni] = nGID;
		StrMakeUpper(strGeneName);
		chp = strrchr(strGeneName, '.');
		if(chp != NULL)
			*chp = '\0';
		StringAddTail(vGeneName+ni, strGeneName);
		ni++;
	}

	fclose(fpIn);

	if(ni != nDatabaseNum)
	{
		printf("error: database size not match!\n");
		exit(EXIT_FAILURE);
	}

	/* map homolo gene */
	fpIn = NULL;
	fpIn = fopen(strInputFile, "r");
	if(fpIn == NULL)
	{
		printf("error: cannot open file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutputFile, "w");
	if(fpOut == NULL)
	{
		printf("error: cannot open file!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%d %d %d %s", &nGID, &nSID, &nHID, strGeneName); 
		if(nSID != nSource)
			continue;

		for(ni=0; ni<nDatabaseNum; ni++)
		{
			if( ( pSpecies->pMatElement[ni] == nSID ) && ( pEntrezID->pMatElement[ni] == nGID) )
			{
				fprintf(fpOut, "%s\t%d\t%d\t%d\t%s\n", vGeneName[ni]->m_pString,
					nGID, nSID, nHID, strGeneName);
			}
		}
	}


	fclose(fpIn);
	fclose(fpOut);


	/* release memory */
	DestroyIntMatrix(pSpecies);
	DestroyIntMatrix(pEntrezID);
	for(ni=0; ni<nDatabaseNum; ni++)
	{
		DeleteString(vGeneName[ni]);
		vGeneName[ni] = NULL;
	}
	free(vGeneName);

	return PROC_SUCCESS;
}

int Feinberg_CollectGenoType()
{
	/* define */
	FILE *fpCEPH;
	FILE *fpData;
	FILE *fpOut;
	struct DOUBLEMATRIX **vData;
	struct BYTEMATRIX **vInfo;
	struct tagString **vSNPID;
	struct tagString **vSampleID;

	struct INTMATRIX *pCol;
	int nDataNum,nSampleNum,nInfoNum,nCEUNum,nCEUInfoNum;
	char strLine[LONG_LINE_LENGTH];
	char strDataFile[MED_LINE_LENGTH];
	char strOutFile[MED_LINE_LENGTH];
	char strCEPHFolder[MED_LINE_LENGTH];
	char strCEPHFile[MED_LINE_LENGTH];
	char strAlias[LINE_LENGTH];
	char strAllels[LINE_LENGTH];
	char strChr[LINE_LENGTH];
	int nPos;
	char strCEUSNP[1024][3];
	char strCEUName[1024][10];
	int ni,nj,nk,nl,nu;
	int nCmpResult;
	int nx,ny,nz;
	char *chp1,*chp2;
	int nFirstLine = 1;
	int nSNPNum;
	int nType;
	int nChr;

	/* init */
	nCEUInfoNum = 11;
	nInfoNum = 3;
	nSampleNum = 0;
	/*strcpy(strDataFile, "D:\\Projects\\alleotyping_project\\data\\Madison_20071109\\Madison_20071109.dat");
	strcpy(strCEPHFolder, "D:\\Projects\\alleotyping_project\\data\\CEPH\\");
	strcpy(strOutFile, "D:\\Projects\\alleotyping_project\\data\\Madison_20071109\\Madison_20071109_genotype.dat");
	*/
	
	strcpy(strDataFile, "/home/bst2/faculty/hji/projects/alleotyping_project/Human_ASE_NG067_011408/Human_ASE_NG067_011408_log2norm_intron.dat");
	strcpy(strCEPHFolder, "/home/bst2/faculty/hji/projects/alleotyping_project/CEPHdata/");
	strcpy(strOutFile, "/home/bst2/faculty/hji/projects/alleotyping_project/Human_ASE_NG067_011408/Human_ASE_NG067_011408_intron_genotype_n.dat");

	

	/* load data */
	fpData = NULL;
	fpData = fopen(strDataFile, "r");
	if(fpData == NULL)
	{
		printf("Error: cannot load data!\n");
		exit(EXIT_FAILURE);
	}
	
	nFirstLine = 1;
	nDataNum = 0;
	nSampleNum = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpData) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(nFirstLine == 1)
		{
			chp1 = strLine;
			chp2 = strchr(chp1, '\t');
			while(chp2 != NULL)
			{
				nSampleNum++;
				chp1 = chp2+1;
				chp2 = strchr(chp1, '\t');
			}
			nSampleNum++;

			nFirstLine = 0;
		}

		if(strLine[0] == '#')
			continue;

		nDataNum++;
	}

	fclose(fpData);
	nSampleNum = nSampleNum-2-nInfoNum; 
	printf("sample = %d\n", nSampleNum);
	printf("data = %d\n", nDataNum);

	/* create memory */
	vData = NULL;
	vData = (struct DOUBLEMATRIX **)calloc(nSampleNum, sizeof(struct DOUBLEMATRIX *));
	if(vData == NULL)
	{
		printf("Error: out of memory!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nSampleNum; ni++)
	{
		vData[ni] = CreateDoubleMatrix(1, nDataNum);
		if(vData[ni] == NULL)
		{
			printf("Error: out of memory!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	vInfo = NULL;
	vInfo = (struct BYTEMATRIX **)calloc(nInfoNum, sizeof(struct BYTEMATRIX *));
	if(vInfo == NULL)
	{
		printf("Error: out of memory!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nInfoNum; ni++)
	{
		vInfo[ni] = CreateByteMatrix(1, nDataNum);
		if(vInfo[ni] == NULL)
		{
			printf("Error: out of memory!\n");
			exit(EXIT_FAILURE);
		}
	}

	vSNPID = NULL;
	vSNPID = (struct tagString **)calloc(nDataNum, sizeof(struct tagString *));
	if(vSNPID == NULL)
	{
		printf("Error: out of memory!\n");
		exit(EXIT_FAILURE);
	}

	vSampleID = NULL;
	vSampleID = (struct tagString **)calloc(nSampleNum, sizeof(struct tagString *));
	if(vSampleID == NULL)
	{
		printf("Error: out of memory!\n");
		exit(EXIT_FAILURE);
	}

	/* load data */
	fpData = NULL;
	fpData = fopen(strDataFile, "r");
	if(fpData == NULL)
	{
		printf("Error: cannot load data!\n");
		exit(EXIT_FAILURE);
	}
	
	ni = 0;
	nFirstLine = 1;
	while(fgets(strLine, LONG_LINE_LENGTH, fpData) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(nFirstLine == 1)
		{
			chp1 = strLine;
			for(nj=0; nj<(2+nInfoNum); nj++)
			{
				chp2 = strchr(chp1, '\t');
				chp1 = chp2+1;
			}

			for(nj=0; nj<nSampleNum; nj++)
			{
				chp2 = strchr(chp1, '\t');
				if(chp2 != NULL)
				{
					*chp2 = '\0';
				}
				StringAddTail(vSampleID+nj, chp1);
				chp1 = chp2+1;
			}

			nFirstLine = 0;
		}

		if(strLine[0] == '#')
			continue;

		chp1 = strLine;
		chp2 = strchr(chp1, '\t');
		chp1 = chp2+1;

		chp2 = strchr(chp1, '\t');
		*chp2 = '\0';
		StringAddTail(vSNPID+ni, chp1);
		chp1 = chp2+1;

		for(nj=0; nj<nInfoNum; nj++)
		{
			chp2 = strchr(chp1, '\t');
			*chp2 = '\0';
			vInfo[nj]->pMatElement[ni] = (unsigned char)(atoi(chp1));
			chp1 = chp2+1;
		}

		for(nj=0; nj<nSampleNum; nj++)
		{
			chp2 = strchr(chp1, '\t');
			if(chp2 != NULL)
				*chp2 = '\0';
			
			vData[nj]->pMatElement[ni] = atof(chp1);
			chp1 = chp2+1;
		}

		ni++;
	}

	fclose(fpData);
	if(ni != nDataNum)
	{
		printf("Error: data not match!\n");
		exit(EXIT_FAILURE);
	}

	/* map  file */
	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: cannot write data!\n");
		exit(EXIT_FAILURE);
	}

	nSNPNum = 0;
	for(ni=1; ni<=24; ni++)
	{
		if(ni == 23)
		{
			sprintf(strCEPHFile, "%sgenotypes_chrX_CEU.txt", strCEPHFolder);
		}
		else if(ni == 24)
		{
			sprintf(strCEPHFile, "%sgenotypes_chrY_CEU.txt", strCEPHFolder);
		}
		else
		{
			sprintf(strCEPHFile, "%sgenotypes_chr%d_CEU.txt", strCEPHFolder, ni);
		}

		fpCEPH = NULL;
		fpCEPH = fopen(strCEPHFile, "r");
		if(fpCEPH == NULL)
		{
			printf("Error: cannot load data!\n");
			exit(EXIT_FAILURE);
		}

		pCol = NULL;
		pCol = CreateIntMatrix(1, nSampleNum);
		if(pCol == NULL)
		{
			printf("Error: cannot load data!\n");
			exit(EXIT_FAILURE);
		}
		for(nx=0; nx<nSampleNum; nx++)
			pCol->pMatElement[nx] = -1;


		fgets(strLine, LONG_LINE_LENGTH, fpCEPH);
		nCEUNum = 0;
		chp1 = strLine;
		chp2 = strpbrk(chp1, " \t");
		while(chp2 != NULL)
		{
			*chp2 = '\0';

			if(nCEUNum>=nCEUInfoNum)
			{
				strcpy(strCEUName[nCEUNum-nCEUInfoNum], chp1);
			}

			nCEUNum++;
			chp1 = chp2+1;
			chp2 = strpbrk(chp1, " \t");
		}
		if(nCEUNum>=nCEUInfoNum)
		{
			strcpy(strCEUName[nCEUNum-nCEUInfoNum], chp1);
		}
		nCEUNum++;
		nCEUNum -= nCEUInfoNum;
		
		for(nx=0; nx<nSampleNum; nx++)
		{
			for(ny=0; ny<nCEUNum; ny++)
			{
				if(strcmp(vSampleID[nx]->m_pString, strCEUName[ny]) == 0)
				{
					pCol->pMatElement[nx] = ny;
					break;
				}
			}
		}

		for(nx=0; nx<nSampleNum; nx++)
		{
			if(pCol->pMatElement[nx] == -1)
				printf("no matching found for column %s !\n", vSampleID[nx]->m_pString);
		}

		while(fgets(strLine, LONG_LINE_LENGTH, fpCEPH) != NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;
			if(strLine[0] == '#')
				continue;


			sscanf(strLine, "%s %s %s %d", strAlias, strAllels, strChr, &nPos);

			nx = 0;
			ny = nDataNum-1;
			while((ny-nx)>1)
			{
				nz = (nx+ny)/2;
				nCmpResult = strcmp(strAlias, vSNPID[nz]->m_pString);

				if(nCmpResult == 0)
					break;
				else if(nCmpResult < 0)
					ny = nz;
				else
					nx = nz;
			}

			if( strcmp(strAlias, vSNPID[nz]->m_pString) != 0)
				continue;

			chp1 = strLine;
			for(nx=0; nx<nCEUInfoNum; nx++)
			{
				chp2 = strpbrk(chp1, " \t");
				chp1 = chp2+1;
			}
			for(nx=0; nx<nCEUNum; nx++)
			{
				chp2 = strpbrk(chp1, " \t");
				if(chp2 != NULL)
					*chp2 = '\0';

				strcpy(strCEUSNP[nx], chp1);
				chp1 = chp2+1;
			}

			for(nx=nz; nx>=0; nx--)
			{
				if(strcmp(strAlias, vSNPID[nx]->m_pString) != 0)
				{
					break;
				}
			}
			nx++;

			for(ny=nz; ny<nDataNum; ny++)
			{
				if(strcmp(strAlias, vSNPID[ny]->m_pString) != 0)
				{
					break;
				}
			}
			ny--;

			strChr[0] = 'c';
			if(strcmp(strChr, "chrX") == 0)
				nChr = 23;
			else if(strcmp(strChr, "chrY") == 0)
				nChr = 24;
			else
				nChr = atoi(strChr+3);

			for(nk=nx; nk<=ny; nk++)
			{
				/* fprintf(fpOut, "%d\t%s\t%s\t%d", nSNPNum, strAlias, strChr, nPos); */
				/* fprintf(fpOut, "%d", nSNPNum); */
				fprintf(fpOut, "%d\t%s\t%d\t%d", nSNPNum, strAlias+2, nChr, nPos);
				
				/* possible alleles */
				nType = 0;
				if(strAllels[0] == 'A')
				{
					if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 0) )
					{
						nType += 1;
					}
					else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 3) )
					{
						nType += 1;
					}
				}
				else if(strAllels[0] == 'C')
				{
					if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 1) )
					{
						nType += 1;
					}
					else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 2) )
					{
						nType += 1;
					}
				}
				else if(strAllels[0] == 'G')
				{
					if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 2) )
					{
						nType += 1;
					}
					else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 1) )
					{
						nType += 1;
					}
				}
				else if(strAllels[0] == 'T')
				{
					if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 3) )
					{
						nType += 1;
					}
					else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 0) )
					{
						nType += 1;
					}
				}

				if(strAllels[2] == 'A')
				{
					if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 0) )
					{
						nType += 2;
					}
					else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 3) )
					{
						nType += 2;
					}
				}
				else if(strAllels[2] == 'C')
				{
					if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 1) )
					{
						nType += 2;
					}
					else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 2) )
					{
						nType += 2;
					}
				}
				else if(strAllels[2] == 'G')
				{
					if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 2) )
					{
						nType += 2;
					}
					else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 1) )
					{
						nType += 2;
					}
				}
				else if(strAllels[2] == 'T')
				{
					if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 3) )
					{
						nType += 2;
					}
					else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 0) )
					{
						nType += 2;
					}
				}

				for(nl=0; nl<nInfoNum; nl++)
				{
					fprintf(fpOut, "\t%d", vInfo[nl]->pMatElement[nk]);
				}
				fprintf(fpOut, "\t%d", nType);

				for(nl=0; nl<nSampleNum; nl++)
				{
					nu = pCol->pMatElement[nl];
					nType = 0;

					if(nu>=0)
					{
						if(strCEUSNP[nu][0] == 'A')
						{
							if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 0) )
							{
								nType += 1;
							}
							else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 3) )
							{
								nType += 1;
							}
						}
						else if(strCEUSNP[nu][0] == 'C')
						{
							if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 1) )
							{
								nType += 1;
							}
							else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 2) )
							{
								nType += 1;
							}
						}
						else if(strCEUSNP[nu][0] == 'G')
						{
							if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 2) )
							{
								nType += 1;
							}
							else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 1) )
							{
								nType += 1;
							}
						}
						else if(strCEUSNP[nu][0] == 'T')
						{
							if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 3) )
							{
								nType += 1;
							}
							else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 0) )
							{
								nType += 1;
							}
						}

						if(strCEUSNP[nu][1] == 'A')
						{
							if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 0) )
							{
								nType += 2;
							}
							else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 3) )
							{
								nType += 2;
							}
						}
						else if(strCEUSNP[nu][1] == 'C')
						{
							if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 1) )
							{
								nType += 2;
							}
							else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 2) )
							{
								nType += 2;
							}
						}
						else if(strCEUSNP[nu][1] == 'G')
						{
							if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 2) )
							{
								nType += 2;
							}
							else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 1) )
							{
								nType += 2;
							}
						}
						else if(strCEUSNP[nu][1] == 'T')
						{
							if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 3) )
							{
								nType += 2;
							}
							else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 0) )
							{
								nType += 2;
							}
						}
					}

					fprintf(fpOut, "\t%f\t%d", vData[nl]->pMatElement[nk], nType);
				}
			
				fprintf(fpOut, "\n");
			}
			
			nSNPNum++;
		}

		fclose(fpCEPH);

		DestroyIntMatrix(pCol);
	}

	fclose(fpOut);

	/* free data matrix */
	for(ni=0; ni<nSampleNum; ni++)
	{
		DestroyDoubleMatrix(vData[ni]);
		vData[ni] = NULL;
	}
	free(vData);

	for(ni=0; ni<nInfoNum; ni++)
	{
		DestroyByteMatrix(vInfo[ni]);
		vInfo[ni] = NULL;
	}
	free(vInfo);

	for(ni=0; ni<nSampleNum; ni++)
	{
		DeleteString(vSampleID[ni]);
		vSampleID[ni] = NULL;
	}
	free(vSampleID);

	for(ni=0; ni<nDataNum; ni++)
	{
		DeleteString(vSNPID[ni]);
		vSNPID[ni] = NULL;
	}
	free(vSNPID);

	/* return */
	return PROC_SUCCESS;
}

int Feinberg_CollectGenoType_mouse()
{
	/* define */
	FILE *fpCEPH;
	FILE *fpData;
	FILE *fpOut;
	struct DOUBLEMATRIX **vData;
	struct BYTEMATRIX **vInfo;
	struct tagString **vSNPID;
	struct tagString **vSampleID;

	int nDataNum,nSampleNum,nInfoNum,nCEUInfoNum;
	char strLine[LONG_LINE_LENGTH];
	char strDataFile[MED_LINE_LENGTH];
	char strOutFile[MED_LINE_LENGTH];
	char strCEPHFolder[MED_LINE_LENGTH];
	char strCEPHFile[MED_LINE_LENGTH];
	char strAlias[LINE_LENGTH];
	char strAllels[LINE_LENGTH];
	char strCEUSNP[1024][3];
	int ni,nj,nk,nl,nu,ne;
	int nCmpResult;
	int nx,ny,nz;
	char *chp1,*chp2;
	int nFirstLine = 1;
	int nSNPNum;
	int nType;

	/* init */
	nCEUInfoNum = 11;
	nInfoNum = 3;
	nSampleNum = 0;
	/*strcpy(strDataFile, "D:\\Projects\\alleotyping_project\\data\\Madison_20071109\\Madison_20071109.dat");
	strcpy(strCEPHFolder, "D:\\Projects\\alleotyping_project\\data\\CEPH\\");
	strcpy(strOutFile, "D:\\Projects\\alleotyping_project\\data\\Madison_20071109\\Madison_20071109_genotype.dat");
	*/
	
	/* strcpy(strDataFile, "/home/bst2/faculty/hji/projects/alleotyping_project/Mousedata/NG067_20071023_exon.dat");
	strcpy(strCEPHFolder, "/home/bst2/faculty/hji/projects/alleotyping_project/MouseGenotype/");
	strcpy(strOutFile, "/home/bst2/faculty/hji/projects/alleotyping_project/Mousedata/NG067_20071023_exon_genotype_X_n.dat");
	*/

	strcpy(strDataFile, "/home/bst2/faculty/hji/projects/alleotyping_project/Mousedata/NG067_20080109_exon.dat");
	strcpy(strCEPHFolder, "/home/bst2/faculty/hji/projects/alleotyping_project/MouseGenotype/");
	strcpy(strOutFile, "/home/bst2/faculty/hji/projects/alleotyping_project/Mousedata/NG067_20080109_exon_genotype_autosome_n.dat");

	

	/* load data */
	fpData = NULL;
	fpData = fopen(strDataFile, "r");
	if(fpData == NULL)
	{
		printf("Error: cannot load data!\n");
		exit(EXIT_FAILURE);
	}
	
	nFirstLine = 1;
	nDataNum = 0;
	nSampleNum = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpData) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(nFirstLine == 1)
		{
			chp1 = strLine;
			chp2 = strchr(chp1, '\t');
			while(chp2 != NULL)
			{
				nSampleNum++;
				chp1 = chp2+1;
				chp2 = strchr(chp1, '\t');
			}
			nSampleNum++;

			nFirstLine = 0;
		}

		if(strLine[0] == '#')
			continue;

		nDataNum++;
	}

	fclose(fpData);
	nSampleNum = nSampleNum-2-nInfoNum; 

	/* create memory */
	vData = NULL;
	vData = (struct DOUBLEMATRIX **)calloc(nSampleNum, sizeof(struct DOUBLEMATRIX *));
	if(vData == NULL)
	{
		printf("Error: out of memory!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nSampleNum; ni++)
	{
		vData[ni] = CreateDoubleMatrix(1, nDataNum);
		if(vData[ni] == NULL)
		{
			printf("Error: out of memory!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	vInfo = NULL;
	vInfo = (struct BYTEMATRIX **)calloc(nInfoNum, sizeof(struct BYTEMATRIX *));
	if(vInfo == NULL)
	{
		printf("Error: out of memory!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nInfoNum; ni++)
	{
		vInfo[ni] = CreateByteMatrix(1, nDataNum);
		if(vInfo[ni] == NULL)
		{
			printf("Error: out of memory!\n");
			exit(EXIT_FAILURE);
		}
	}

	vSNPID = NULL;
	vSNPID = (struct tagString **)calloc(nDataNum, sizeof(struct tagString *));
	if(vSNPID == NULL)
	{
		printf("Error: out of memory!\n");
		exit(EXIT_FAILURE);
	}

	vSampleID = NULL;
	vSampleID = (struct tagString **)calloc(nSampleNum, sizeof(struct tagString *));
	if(vSampleID == NULL)
	{
		printf("Error: out of memory!\n");
		exit(EXIT_FAILURE);
	}

	/* load data */
	fpData = NULL;
	fpData = fopen(strDataFile, "r");
	if(fpData == NULL)
	{
		printf("Error: cannot load data!\n");
		exit(EXIT_FAILURE);
	}
	
	ni = 0;
	nFirstLine = 1;
	while(fgets(strLine, LONG_LINE_LENGTH, fpData) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(nFirstLine == 1)
		{
			chp1 = strLine;
			for(nj=0; nj<(2+nInfoNum); nj++)
			{
				chp2 = strchr(chp1, '\t');
				chp1 = chp2+1;
			}

			for(nj=0; nj<nSampleNum; nj++)
			{
				chp2 = strchr(chp1, '\t');
				if(chp2 != NULL)
				{
					*chp2 = '\0';
				}
				StringAddTail(vSampleID+nj, chp1);
				chp1 = chp2+1;
			}

			nFirstLine = 0;
		}

		if(strLine[0] == '#')
			continue;

		chp1 = strLine;
		chp2 = strchr(chp1, '\t');
		chp1 = chp2+1;

		chp2 = strchr(chp1, '\t');
		*chp2 = '\0';
		StringAddTail(vSNPID+ni, chp1);
		chp1 = chp2+1;

		for(nj=0; nj<nInfoNum; nj++)
		{
			chp2 = strchr(chp1, '\t');
			*chp2 = '\0';
			vInfo[nj]->pMatElement[ni] = (unsigned char)(atoi(chp1));
			chp1 = chp2+1;
		}

		for(nj=0; nj<nSampleNum; nj++)
		{
			chp2 = strchr(chp1, '\t');
			if(chp2 != NULL)
				*chp2 = '\0';
			
			vData[nj]->pMatElement[ni] = atof(chp1);
			chp1 = chp2+1;
		}

		ni++;
	}

	fclose(fpData);
	if(ni != nDataNum)
	{
		printf("Error: data not match!\n");
		exit(EXIT_FAILURE);
	}

	/* map  file */
	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: cannot write data!\n");
		exit(EXIT_FAILURE);
	}

	nSNPNum = 0;
	for(ni=1; ni<=19; ni++)
	{
		for(ne=0; ne<2; ne++)
		{
			if(ne == 0)
			{
				if(ni == 20)
				{
					sprintf(strCEPHFile, "%sChrX_CAST_C57B6_coding_genotype.txt", strCEPHFolder);
				}
				else
				{
					sprintf(strCEPHFile, "%sChr%d_CAST_C57B6_coding_genotype.txt", strCEPHFolder, ni);
				}
			}
			else
			{
				if(ni == 20)
				{
					sprintf(strCEPHFile, "%sChrX_CAST_C57B6_intron_genotype.txt", strCEPHFolder);
				}
				else
				{
					sprintf(strCEPHFile, "%sChr%d_CAST_C57B6_intron_genotype.txt", strCEPHFolder, ni);
				}
			}

			fpCEPH = NULL;
			fpCEPH = fopen(strCEPHFile, "r");
			if(fpCEPH == NULL)
			{
				printf("Error: cannot load data!\n");
				exit(EXIT_FAILURE);
			}

			
			fgets(strLine, LONG_LINE_LENGTH, fpCEPH);
			
			while(fgets(strLine, LONG_LINE_LENGTH, fpCEPH) != NULL)
			{
				StrTrimLeft(strLine);
				StrTrimRight(strLine);
				if(strLine[0] == '\0')
					continue;
				if(strLine[0] == '#')
					continue;


				sscanf(strLine, "%s %s", strAlias, strAllels);
				strcpy(strCEUSNP[0], strAllels);

				nx = 0;
				ny = nDataNum-1;
				while((ny-nx)>1)
				{
					nz = (nx+ny)/2;
					nCmpResult = strcmp(strAlias, vSNPID[nz]->m_pString);

					if(nCmpResult == 0)
						break;
					else if(nCmpResult < 0)
						ny = nz;
					else
						nx = nz;
				}

				if( strcmp(strAlias, vSNPID[nz]->m_pString) != 0)
					continue;

				for(nx=nz; nx>=0; nx--)
				{
					if(strcmp(strAlias, vSNPID[nx]->m_pString) != 0)
					{
						break;
					}
				}
				nx++;

				for(ny=nz; ny<nDataNum; ny++)
				{
					if(strcmp(strAlias, vSNPID[ny]->m_pString) != 0)
					{
						break;
					}
				}
				ny--;

				for(nk=nx; nk<=ny; nk++)
				{
					fprintf(fpOut, "%d\t%s", nSNPNum, strAlias+2);
					/* fprintf(fpOut, "%d", nSNPNum); */
					/* fprintf(fpOut, "%s\t%s\t%d", strAlias, strChr, nPos); */
					
					for(nl=0; nl<nInfoNum; nl++)
					{
						fprintf(fpOut, "\t%d", vInfo[nl]->pMatElement[nk]);
					}

					for(nl=0; nl<nSampleNum; nl++)
					{
						nu = 0;
						nType = 0;
						if(strCEUSNP[nu][0] == 'A')
						{
							if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 0) )
							{
								nType += 1;
							}
							else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 3) )
							{
								nType += 1;
							}
						}
						else if(strCEUSNP[nu][0] == 'C')
						{
							if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 1) )
							{
								nType += 1;
							}
							else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 2) )
							{
								nType += 1;
							}
						}
						else if(strCEUSNP[nu][0] == 'G')
						{
							if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 2) )
							{
								nType += 1;
							}
							else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 1) )
							{
								nType += 1;
							}
						}
						else if(strCEUSNP[nu][0] == 'T')
						{
							if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 3) )
							{
								nType += 1;
							}
							else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 0) )
							{
								nType += 1;
							}
						}

						if(strCEUSNP[nu][1] == 'A')
						{
							if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 0) )
							{
								nType += 2;
							}
							else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 3) )
							{
								nType += 2;
							}
						}
						else if(strCEUSNP[nu][1] == 'C')
						{
							if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 1) )
							{
								nType += 2;
							}
							else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 2) )
							{
								nType += 2;
							}
						}
						else if(strCEUSNP[nu][1] == 'G')
						{
							if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 2) )
							{
								nType += 2;
							}
							else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 1) )
							{
								nType += 2;
							}
						}
						else if(strCEUSNP[nu][1] == 'T')
						{
							if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 3) )
							{
								nType += 2;
							}
							else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 0) )
							{
								nType += 2;
							}
						}

						fprintf(fpOut, "\t%f\t%d", vData[nl]->pMatElement[nk], nType);
					}
				
					fprintf(fpOut, "\n");
				}
				
				nSNPNum++;
			}

			fclose(fpCEPH);

		}
	}

	fclose(fpOut);

	/* free data matrix */
	for(ni=0; ni<nSampleNum; ni++)
	{
		DestroyDoubleMatrix(vData[ni]);
		vData[ni] = NULL;
	}
	free(vData);

	for(ni=0; ni<nInfoNum; ni++)
	{
		DestroyByteMatrix(vInfo[ni]);
		vInfo[ni] = NULL;
	}
	free(vInfo);

	for(ni=0; ni<nSampleNum; ni++)
	{
		DeleteString(vSampleID[ni]);
		vSampleID[ni] = NULL;
	}
	free(vSampleID);

	for(ni=0; ni<nDataNum; ni++)
	{
		DeleteString(vSNPID[ni]);
		vSNPID[ni] = NULL;
	}
	free(vSNPID);

	/* return */
	return PROC_SUCCESS;
}

int Feinberg_CollectGenoType_ASM()
{
	/* define */
	FILE *fpCEPH;
	FILE *fpData;
	FILE *fpOut;
	struct DOUBLEMATRIX **vData;
	struct BYTEMATRIX **vInfo;
	struct tagString **vSNPID;
	struct tagString **vSampleID;

	struct INTMATRIX *pCol;
	int nDataNum,nSampleNum,nInfoNum,nCEUNum,nCEUInfoNum;
	char strLine[LONG_LINE_LENGTH];
	char strDataFile[MED_LINE_LENGTH];
	char strOutFile[MED_LINE_LENGTH];
	char strCEPHFolder[MED_LINE_LENGTH];
	char strCEPHFile[MED_LINE_LENGTH];
	char strAlias[LINE_LENGTH];
	char strAllels[LINE_LENGTH];
	char strChr[LINE_LENGTH];
	int nPos;
	char strCEUSNP[1024][3];
	char strCEUName[1024][10];
	int ni,nj,nk,nl,nu;
	int nCmpResult;
	int nx,ny,nz;
	char *chp1,*chp2;
	int nFirstLine = 1;
	int nSNPNum;
	int nType;
	int nChr;

	/* init */
	nCEUInfoNum = 11;
	nInfoNum = 3;
	nSampleNum = 0;
	/*strcpy(strDataFile, "D:\\Projects\\alleotyping_project\\data\\Madison_20071109\\Madison_20071109.dat");
	strcpy(strCEPHFolder, "D:\\Projects\\alleotyping_project\\data\\CEPH\\");
	strcpy(strOutFile, "D:\\Projects\\alleotyping_project\\data\\Madison_20071109\\Madison_20071109_genotype.dat");
	*/
	
	strcpy(strDataFile, "/home/bst2/faculty/hji/projects/alleotyping_project/Human_ASM_NG084_20080510/Human_ASM_NG084_20080510_log2norm_data.txt");
	strcpy(strCEPHFolder, "/home/bst2/faculty/hji/projects/alleotyping_project/CEPHdata/");
	strcpy(strOutFile, "/home/bst2/faculty/hji/projects/alleotyping_project/Human_ASM_NG084_20080510/Human_ASM_NG084_20080510_log2norm_data_genotype_n.dat");
	

	/* load data */
	fpData = NULL;
	fpData = fopen(strDataFile, "r");
	if(fpData == NULL)
	{
		printf("Error: cannot load data!\n");
		exit(EXIT_FAILURE);
	}
	
	nFirstLine = 1;
	nDataNum = 0;
	nSampleNum = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpData) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(nFirstLine == 1)
		{
			chp1 = strLine;
			chp2 = strchr(chp1, '\t');
			while(chp2 != NULL)
			{
				nSampleNum++;
				chp1 = chp2+1;
				chp2 = strchr(chp1, '\t');
			}
			nSampleNum++;

			nFirstLine = 0;
		}

		if(strLine[0] == '#')
			continue;

		nDataNum++;
	}

	fclose(fpData);
	nSampleNum = nSampleNum-2-nInfoNum; 
	printf("sample = %d\n", nSampleNum);
	printf("data = %d\n", nDataNum);

	/* create memory */
	vData = NULL;
	vData = (struct DOUBLEMATRIX **)calloc(nSampleNum, sizeof(struct DOUBLEMATRIX *));
	if(vData == NULL)
	{
		printf("Error: out of memory!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nSampleNum; ni++)
	{
		vData[ni] = CreateDoubleMatrix(1, nDataNum);
		if(vData[ni] == NULL)
		{
			printf("Error: out of memory!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	vInfo = NULL;
	vInfo = (struct BYTEMATRIX **)calloc(nInfoNum, sizeof(struct BYTEMATRIX *));
	if(vInfo == NULL)
	{
		printf("Error: out of memory!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nInfoNum; ni++)
	{
		vInfo[ni] = CreateByteMatrix(1, nDataNum);
		if(vInfo[ni] == NULL)
		{
			printf("Error: out of memory!\n");
			exit(EXIT_FAILURE);
		}
	}

	vSNPID = NULL;
	vSNPID = (struct tagString **)calloc(nDataNum, sizeof(struct tagString *));
	if(vSNPID == NULL)
	{
		printf("Error: out of memory!\n");
		exit(EXIT_FAILURE);
	}

	vSampleID = NULL;
	vSampleID = (struct tagString **)calloc(nSampleNum, sizeof(struct tagString *));
	if(vSampleID == NULL)
	{
		printf("Error: out of memory!\n");
		exit(EXIT_FAILURE);
	}

	/* load data */
	fpData = NULL;
	fpData = fopen(strDataFile, "r");
	if(fpData == NULL)
	{
		printf("Error: cannot load data!\n");
		exit(EXIT_FAILURE);
	}
	
	ni = 0;
	nFirstLine = 1;
	while(fgets(strLine, LONG_LINE_LENGTH, fpData) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(nFirstLine == 1)
		{
			chp1 = strLine;
			for(nj=0; nj<(2+nInfoNum); nj++)
			{
				chp2 = strchr(chp1, '\t');
				chp1 = chp2+1;
			}

			for(nj=0; nj<nSampleNum; nj++)
			{
				chp2 = strchr(chp1, '\t');
				if(chp2 != NULL)
				{
					*chp2 = '\0';
				}
				StringAddTail(vSampleID+nj, chp1);
				chp1 = chp2+1;
			}

			nFirstLine = 0;
		}

		if(strLine[0] == '#')
			continue;

		chp1 = strLine;
		chp2 = strchr(chp1, '\t');
		chp1 = chp2+1;

		chp2 = strchr(chp1, '\t');
		*chp2 = '\0';
		StringAddTail(vSNPID+ni, chp1);
		chp1 = chp2+1;

		for(nj=0; nj<nInfoNum; nj++)
		{
			chp2 = strchr(chp1, '\t');
			*chp2 = '\0';
			vInfo[nj]->pMatElement[ni] = (unsigned char)(atoi(chp1));
			chp1 = chp2+1;
		}

		for(nj=0; nj<nSampleNum; nj++)
		{
			chp2 = strchr(chp1, '\t');
			if(chp2 != NULL)
				*chp2 = '\0';
			
			vData[nj]->pMatElement[ni] = atof(chp1);
			chp1 = chp2+1;
		}

		ni++;
	}

	fclose(fpData);
	if(ni != nDataNum)
	{
		printf("Error: data not match!\n");
		exit(EXIT_FAILURE);
	}

	/* map  file */
	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: cannot write data!\n");
		exit(EXIT_FAILURE);
	}

	nSNPNum = 0;
	for(ni=1; ni<=24; ni++)
	{
		if(ni == 23)
		{
			sprintf(strCEPHFile, "%sgenotypes_chrX_CEU.txt", strCEPHFolder);
		}
		else if(ni == 24)
		{
			sprintf(strCEPHFile, "%sgenotypes_chrY_CEU.txt", strCEPHFolder);
		}
		else
		{
			sprintf(strCEPHFile, "%sgenotypes_chr%d_CEU.txt", strCEPHFolder, ni);
		}

		fpCEPH = NULL;
		fpCEPH = fopen(strCEPHFile, "r");
		if(fpCEPH == NULL)
		{
			printf("Error: cannot load data!\n");
			exit(EXIT_FAILURE);
		}

		pCol = NULL;
		pCol = CreateIntMatrix(1, nSampleNum);
		if(pCol == NULL)
		{
			printf("Error: cannot load data!\n");
			exit(EXIT_FAILURE);
		}
		for(nx=0; nx<nSampleNum; nx++)
			pCol->pMatElement[nx] = -1;


		fgets(strLine, LONG_LINE_LENGTH, fpCEPH);
		nCEUNum = 0;
		chp1 = strLine;
		chp2 = strpbrk(chp1, " \t");
		while(chp2 != NULL)
		{
			*chp2 = '\0';

			if(nCEUNum>=nCEUInfoNum)
			{
				strcpy(strCEUName[nCEUNum-nCEUInfoNum], chp1);
			}

			nCEUNum++;
			chp1 = chp2+1;
			chp2 = strpbrk(chp1, " \t");
		}
		if(nCEUNum>=nCEUInfoNum)
		{
			strcpy(strCEUName[nCEUNum-nCEUInfoNum], chp1);
		}
		nCEUNum++;
		nCEUNum -= nCEUInfoNum;
		
		for(nx=0; nx<nSampleNum; nx++)
		{
			for(ny=0; ny<nCEUNum; ny++)
			{
				if(strcmp(vSampleID[nx]->m_pString, strCEUName[ny]) == 0)
				{
					pCol->pMatElement[nx] = ny;
					break;
				}
			}
		}

		for(nx=0; nx<nSampleNum; nx++)
		{
			if(pCol->pMatElement[nx] == -1)
				printf("no matching found for column %s !\n", vSampleID[nx]->m_pString);
		}

		while(fgets(strLine, LONG_LINE_LENGTH, fpCEPH) != NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;
			if(strLine[0] == '#')
				continue;


			sscanf(strLine, "%s %s %s %d", strAlias, strAllels, strChr, &nPos);

			nx = 0;
			ny = nDataNum-1;
			while((ny-nx)>1)
			{
				nz = (nx+ny)/2;
				nCmpResult = strcmp(strAlias, vSNPID[nz]->m_pString);

				if(nCmpResult == 0)
					break;
				else if(nCmpResult < 0)
					ny = nz;
				else
					nx = nz;
			}

			if( strcmp(strAlias, vSNPID[nz]->m_pString) != 0)
				continue;

			chp1 = strLine;
			for(nx=0; nx<nCEUInfoNum; nx++)
			{
				chp2 = strpbrk(chp1, " \t");
				chp1 = chp2+1;
			}
			for(nx=0; nx<nCEUNum; nx++)
			{
				chp2 = strpbrk(chp1, " \t");
				if(chp2 != NULL)
					*chp2 = '\0';

				strcpy(strCEUSNP[nx], chp1);
				chp1 = chp2+1;
			}

			for(nx=nz; nx>=0; nx--)
			{
				if(strcmp(strAlias, vSNPID[nx]->m_pString) != 0)
				{
					break;
				}
			}
			nx++;

			for(ny=nz; ny<nDataNum; ny++)
			{
				if(strcmp(strAlias, vSNPID[ny]->m_pString) != 0)
				{
					break;
				}
			}
			ny--;

			strChr[0] = 'c';
			if(strcmp(strChr, "chrX") == 0)
				nChr = 23;
			else if(strcmp(strChr, "chrY") == 0)
				nChr = 24;
			else
				nChr = atoi(strChr+3);

			for(nk=nx; nk<=ny; nk++)
			{
				/* fprintf(fpOut, "%d\t%s\t%s\t%d", nSNPNum, strAlias, strChr, nPos); */
				/* fprintf(fpOut, "%d", nSNPNum); */
				fprintf(fpOut, "%d\t%s", nSNPNum, strAlias+2);
				
				/* possible alleles */
				nType = 0;
				if(strAllels[0] == 'A')
				{
					if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 0) )
					{
						nType += 1;
					}
					else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 3) )
					{
						nType += 1;
					}
				}
				else if(strAllels[0] == 'C')
				{
					if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 1) )
					{
						nType += 1;
					}
					else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 2) )
					{
						nType += 1;
					}
				}
				else if(strAllels[0] == 'G')
				{
					if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 2) )
					{
						nType += 1;
					}
					else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 1) )
					{
						nType += 1;
					}
				}
				else if(strAllels[0] == 'T')
				{
					if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 3) )
					{
						nType += 1;
					}
					else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 0) )
					{
						nType += 1;
					}
				}

				if(strAllels[2] == 'A')
				{
					if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 0) )
					{
						nType += 2;
					}
					else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 3) )
					{
						nType += 2;
					}
				}
				else if(strAllels[2] == 'C')
				{
					if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 1) )
					{
						nType += 2;
					}
					else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 2) )
					{
						nType += 2;
					}
				}
				else if(strAllels[2] == 'G')
				{
					if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 2) )
					{
						nType += 2;
					}
					else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 1) )
					{
						nType += 2;
					}
				}
				else if(strAllels[2] == 'T')
				{
					if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 3) )
					{
						nType += 2;
					}
					else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 0) )
					{
						nType += 2;
					}
				}

				for(nl=0; nl<nInfoNum; nl++)
				{
					fprintf(fpOut, "\t%d", vInfo[nl]->pMatElement[nk]);
				}
				fprintf(fpOut, "\t%d", nType);

				for(nl=0; nl<nSampleNum; nl++)
				{
					nu = pCol->pMatElement[nl];
					nType = 0;

					if(nu>=0)
					{
						if(strCEUSNP[nu][0] == 'A')
						{
							if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 0) )
							{
								nType += 1;
							}
							else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 3) )
							{
								nType += 1;
							}
						}
						else if(strCEUSNP[nu][0] == 'C')
						{
							if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 1) )
							{
								nType += 1;
							}
							else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 2) )
							{
								nType += 1;
							}
						}
						else if(strCEUSNP[nu][0] == 'G')
						{
							if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 2) )
							{
								nType += 1;
							}
							else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 1) )
							{
								nType += 1;
							}
						}
						else if(strCEUSNP[nu][0] == 'T')
						{
							if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 3) )
							{
								nType += 1;
							}
							else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 0) )
							{
								nType += 1;
							}
						}

						if(strCEUSNP[nu][1] == 'A')
						{
							if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 0) )
							{
								nType += 2;
							}
							else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 3) )
							{
								nType += 2;
							}
						}
						else if(strCEUSNP[nu][1] == 'C')
						{
							if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 1) )
							{
								nType += 2;
							}
							else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 2) )
							{
								nType += 2;
							}
						}
						else if(strCEUSNP[nu][1] == 'G')
						{
							if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 2) )
							{
								nType += 2;
							}
							else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 1) )
							{
								nType += 2;
							}
						}
						else if(strCEUSNP[nu][1] == 'T')
						{
							if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 3) )
							{
								nType += 2;
							}
							else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 0) )
							{
								nType += 2;
							}
						}

						fprintf(fpOut, "\t%f\t%d", vData[nl]->pMatElement[nk], nType);
					}
					else
					{
						fprintf(fpOut, "\t%f", vData[nl]->pMatElement[nk]);
					}				
				}
			
				fprintf(fpOut, "\n");
			}
			
			nSNPNum++;
		}

		fclose(fpCEPH);

		DestroyIntMatrix(pCol);
	}

	fclose(fpOut);

	/* free data matrix */
	for(ni=0; ni<nSampleNum; ni++)
	{
		DestroyDoubleMatrix(vData[ni]);
		vData[ni] = NULL;
	}
	free(vData);

	for(ni=0; ni<nInfoNum; ni++)
	{
		DestroyByteMatrix(vInfo[ni]);
		vInfo[ni] = NULL;
	}
	free(vInfo);

	for(ni=0; ni<nSampleNum; ni++)
	{
		DeleteString(vSampleID[ni]);
		vSampleID[ni] = NULL;
	}
	free(vSampleID);

	for(ni=0; ni<nDataNum; ni++)
	{
		DeleteString(vSNPID[ni]);
		vSNPID[ni] = NULL;
	}
	free(vSNPID);

	/* return */
	return PROC_SUCCESS;
}

int Feinberg_CollectGenoType_Tumor()
{
	/* define */
	FILE *fpCEPH;
	FILE *fpData;
	FILE *fpOut;
	struct DOUBLEMATRIX **vData;
	struct BYTEMATRIX **vInfo;
	struct tagString **vSNPID;
	struct tagString **vSampleID;

	/* struct INTMATRIX *pCol; */
	int nDataNum,nSampleNum,nInfoNum,nCEUNum,nCEUInfoNum;
	char strLine[LONG_LINE_LENGTH];
	char strDataFile[MED_LINE_LENGTH];
	char strOutFile[MED_LINE_LENGTH];
	char strCEPHFolder[MED_LINE_LENGTH];
	char strCEPHFile[MED_LINE_LENGTH];
	char strAlias[LINE_LENGTH];
	char strAllels[LINE_LENGTH];
	char strChr[LINE_LENGTH];
	int nPos;
	/* char strCEUSNP[1024][3]; */
	char strCEUName[1024][10];
	int ni,nj,nk,nl,nu;
	int nCmpResult;
	int nx,ny,nz;
	char *chp1,*chp2;
	int nFirstLine = 1;
	int nSNPNum;
	int nType;
	int nChr;

	/* init */
	nCEUInfoNum = 11;
	nInfoNum = 3;
	nSampleNum = 0;
		
	/* strcpy(strDataFile, "/nexsan/bst2/cisreg/feinberg/NG087_20080613/analysis/Human_ASE_NG087_061908_log2norm_exon.txt");
	strcpy(strCEPHFolder, "/nexsan/bst2/cisreg/feinberg/HapMap/genotypes/");
	strcpy(strOutFile, "/nexsan/bst2/cisreg/feinberg/NG087_20080613/analysis/Human_ASE_NG087_061908_exon_genotype_n.dat");
	*/

	/* strcpy(strDataFile, "/nexsan/bst2/cisreg/feinberg/NG089_20080627/analysis/Human_ASE_NG089_norm_exon.txt");
	strcpy(strCEPHFolder, "/nexsan/bst2/cisreg/feinberg/HapMap/genotypes/");
	strcpy(strOutFile, "/nexsan/bst2/cisreg/feinberg/NG089_20080627/analysis/Human_ASE_NG089_norm_exon_genotype_n.txt");
	*/

	strcpy(strDataFile, "/nexsan/bst2/cisreg/feinberg/NG107_20091022/analysis/Human_ASE_NG107_norm_exon.txt");
	strcpy(strCEPHFolder, "/nexsan/bst2/cisreg/feinberg/HapMap/genotypes/");
	strcpy(strOutFile, "/nexsan/bst2/cisreg/feinberg/NG107_20091022/analysis/Human_ASE_NG107_norm_exon_genotype_n.txt");
	

	/* load data */
	fpData = NULL;
	fpData = fopen(strDataFile, "r");
	if(fpData == NULL)
	{
		printf("Error: cannot load data!\n");
		exit(EXIT_FAILURE);
	}
	
	nFirstLine = 1;
	nDataNum = 0;
	nSampleNum = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpData) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(nFirstLine == 1)
		{
			chp1 = strLine;
			chp2 = strchr(chp1, '\t');
			while(chp2 != NULL)
			{
				nSampleNum++;
				chp1 = chp2+1;
				chp2 = strchr(chp1, '\t');
			}
			nSampleNum++;

			nFirstLine = 0;
		}

		if(strLine[0] == '#')
			continue;

		nDataNum++;
	}

	fclose(fpData);
	nSampleNum = nSampleNum-2-nInfoNum; 
	printf("sample = %d\n", nSampleNum);
	printf("data = %d\n", nDataNum);

	/* create memory */
	vData = NULL;
	vData = (struct DOUBLEMATRIX **)calloc(nSampleNum, sizeof(struct DOUBLEMATRIX *));
	if(vData == NULL)
	{
		printf("Error: out of memory!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nSampleNum; ni++)
	{
		vData[ni] = CreateDoubleMatrix(1, nDataNum);
		if(vData[ni] == NULL)
		{
			printf("Error: out of memory!\n");
			exit(EXIT_FAILURE);
		}
	}
	
	vInfo = NULL;
	vInfo = (struct BYTEMATRIX **)calloc(nInfoNum, sizeof(struct BYTEMATRIX *));
	if(vInfo == NULL)
	{
		printf("Error: out of memory!\n");
		exit(EXIT_FAILURE);
	}
	for(ni=0; ni<nInfoNum; ni++)
	{
		vInfo[ni] = CreateByteMatrix(1, nDataNum);
		if(vInfo[ni] == NULL)
		{
			printf("Error: out of memory!\n");
			exit(EXIT_FAILURE);
		}
	}

	vSNPID = NULL;
	vSNPID = (struct tagString **)calloc(nDataNum, sizeof(struct tagString *));
	if(vSNPID == NULL)
	{
		printf("Error: out of memory!\n");
		exit(EXIT_FAILURE);
	}

	vSampleID = NULL;
	vSampleID = (struct tagString **)calloc(nSampleNum, sizeof(struct tagString *));
	if(vSampleID == NULL)
	{
		printf("Error: out of memory!\n");
		exit(EXIT_FAILURE);
	}

	/* load data */
	fpData = NULL;
	fpData = fopen(strDataFile, "r");
	if(fpData == NULL)
	{
		printf("Error: cannot load data!\n");
		exit(EXIT_FAILURE);
	}
	
	ni = 0;
	nFirstLine = 1;
	while(fgets(strLine, LONG_LINE_LENGTH, fpData) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		if(nFirstLine == 1)
		{
			chp1 = strLine;
			for(nj=0; nj<(2+nInfoNum); nj++)
			{
				chp2 = strchr(chp1, '\t');
				chp1 = chp2+1;
			}

			for(nj=0; nj<nSampleNum; nj++)
			{
				chp2 = strchr(chp1, '\t');
				if(chp2 != NULL)
				{
					*chp2 = '\0';
				}
				StringAddTail(vSampleID+nj, chp1);
				chp1 = chp2+1;
			}

			nFirstLine = 0;
		}

		if(strLine[0] == '#')
			continue;

		chp1 = strLine;
		chp2 = strchr(chp1, '\t');
		chp1 = chp2+1;

		chp2 = strchr(chp1, '\t');
		*chp2 = '\0';
		StringAddTail(vSNPID+ni, chp1);
		chp1 = chp2+1;

		for(nj=0; nj<nInfoNum; nj++)
		{
			chp2 = strchr(chp1, '\t');
			*chp2 = '\0';
			vInfo[nj]->pMatElement[ni] = (unsigned char)(atoi(chp1));
			chp1 = chp2+1;
		}

		for(nj=0; nj<nSampleNum; nj++)
		{
			chp2 = strchr(chp1, '\t');
			if(chp2 != NULL)
				*chp2 = '\0';
			
			vData[nj]->pMatElement[ni] = atof(chp1);
			chp1 = chp2+1;
		}

		ni++;
	}

	fclose(fpData);
	if(ni != nDataNum)
	{
		printf("Error: data not match!\n");
		exit(EXIT_FAILURE);
	}

	/* map  file */
	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: cannot write data!\n");
		exit(EXIT_FAILURE);
	}

	nSNPNum = 0;
	for(ni=1; ni<=24; ni++)
	{
		if(ni == 23)
		{
			sprintf(strCEPHFile, "%sgenotypes_chrX_CEU_r23a_nr.b36_fwd.txt", strCEPHFolder);
		}
		else if(ni == 24)
		{
			sprintf(strCEPHFile, "%sgenotypes_chrY_CEU_r23a_nr.b36_fwd.txt", strCEPHFolder);
		}
		else
		{
			sprintf(strCEPHFile, "%sgenotypes_chr%d_CEU_r23a_nr.b36_fwd.txt", strCEPHFolder, ni);
		}

		fpCEPH = NULL;
		fpCEPH = fopen(strCEPHFile, "r");
		if(fpCEPH == NULL)
		{
			printf("Error: cannot load data!\n");
			exit(EXIT_FAILURE);
		}

		
		fgets(strLine, LONG_LINE_LENGTH, fpCEPH);
		nCEUNum = 0;
		chp1 = strLine;
		chp2 = strpbrk(chp1, " \t");
		while(chp2 != NULL)
		{
			*chp2 = '\0';

			if(nCEUNum>=nCEUInfoNum)
			{
				strcpy(strCEUName[nCEUNum-nCEUInfoNum], chp1);
			}

			nCEUNum++;
			chp1 = chp2+1;
			chp2 = strpbrk(chp1, " \t");
		}
		if(nCEUNum>=nCEUInfoNum)
		{
			strcpy(strCEUName[nCEUNum-nCEUInfoNum], chp1);
		}
		nCEUNum++;
		nCEUNum -= nCEUInfoNum;
		
		while(fgets(strLine, LONG_LINE_LENGTH, fpCEPH) != NULL)
		{
			StrTrimLeft(strLine);
			StrTrimRight(strLine);
			if(strLine[0] == '\0')
				continue;
			if(strLine[0] == '#')
				continue;


			sscanf(strLine, "%s %s %s %d", strAlias, strAllels, strChr, &nPos);

			nx = 0;
			ny = nDataNum-1;
			while((ny-nx)>1)
			{
				nz = (nx+ny)/2;
				nCmpResult = strcmp(strAlias, vSNPID[nz]->m_pString);

				if(nCmpResult == 0)
					break;
				else if(nCmpResult < 0)
					ny = nz;
				else
					nx = nz;
			}

			if( strcmp(strAlias, vSNPID[nz]->m_pString) != 0)
				continue;

			for(nx=nz; nx>=0; nx--)
			{
				if(strcmp(strAlias, vSNPID[nx]->m_pString) != 0)
				{
					break;
				}
			}
			nx++;

			for(ny=nz; ny<nDataNum; ny++)
			{
				if(strcmp(strAlias, vSNPID[ny]->m_pString) != 0)
				{
					break;
				}
			}
			ny--;

			strChr[0] = 'c';
			if(strcmp(strChr, "chrX") == 0)
				nChr = 23;
			else if(strcmp(strChr, "chrY") == 0)
				nChr = 24;
			else
				nChr = atoi(strChr+3);

			for(nk=nx; nk<=ny; nk++)
			{
				/* fprintf(fpOut, "%d\t%s\t%s\t%d", nSNPNum, strAlias, strChr, nPos); */
				/* fprintf(fpOut, "%d", nSNPNum); */
				fprintf(fpOut, "%d\t%s\t%d\t%d", nSNPNum, strAlias+2, nChr, nPos);
				
				/* possible alleles */
				nType = 0;
				if(strAllels[0] == 'A')
				{
					if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 0) )
					{
						nType += 1;
					}
					else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 3) )
					{
						nType += 1;
					}
				}
				else if(strAllels[0] == 'C')
				{
					if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 1) )
					{
						nType += 1;
					}
					else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 2) )
					{
						nType += 1;
					}
				}
				else if(strAllels[0] == 'G')
				{
					if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 2) )
					{
						nType += 1;
					}
					else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 1) )
					{
						nType += 1;
					}
				}
				else if(strAllels[0] == 'T')
				{
					if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 3) )
					{
						nType += 1;
					}
					else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 0) )
					{
						nType += 1;
					}
				}

				if(strAllels[2] == 'A')
				{
					if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 0) )
					{
						nType += 2;
					}
					else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 3) )
					{
						nType += 2;
					}
				}
				else if(strAllels[2] == 'C')
				{
					if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 1) )
					{
						nType += 2;
					}
					else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 2) )
					{
						nType += 2;
					}
				}
				else if(strAllels[2] == 'G')
				{
					if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 2) )
					{
						nType += 2;
					}
					else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 1) )
					{
						nType += 2;
					}
				}
				else if(strAllels[2] == 'T')
				{
					if( (vInfo[0]->pMatElement[nk] == 0) && (vInfo[2]->pMatElement[nk] == 3) )
					{
						nType += 2;
					}
					else if( (vInfo[0]->pMatElement[nk] == 1) && (vInfo[2]->pMatElement[nk] == 0) )
					{
						nType += 2;
					}
				}

				for(nl=0; nl<nInfoNum; nl++)
				{
					fprintf(fpOut, "\t%d", vInfo[nl]->pMatElement[nk]);
				}
				fprintf(fpOut, "\t%d", nType);

				for(nl=0; nl<nSampleNum; nl++)
				{
					fprintf(fpOut, "\t%f", vData[nl]->pMatElement[nk]);
				}
			
				fprintf(fpOut, "\n");
			}
			
			nSNPNum++;
		}

		fclose(fpCEPH);

	}

	fclose(fpOut);

	/* free data matrix */
	for(ni=0; ni<nSampleNum; ni++)
	{
		DestroyDoubleMatrix(vData[ni]);
		vData[ni] = NULL;
	}
	free(vData);

	for(ni=0; ni<nInfoNum; ni++)
	{
		DestroyByteMatrix(vInfo[ni]);
		vInfo[ni] = NULL;
	}
	free(vInfo);

	for(ni=0; ni<nSampleNum; ni++)
	{
		DeleteString(vSampleID[ni]);
		vSampleID[ni] = NULL;
	}
	free(vSampleID);

	for(ni=0; ni<nDataNum; ni++)
	{
		DeleteString(vSNPID[ni]);
		vSNPID[ni] = NULL;
	}
	free(vSNPID);

	/* return */
	return PROC_SUCCESS;
}

int HTS_selectread(char strInputFile[], char strOutputFile[], double dR)
{
	FILE *fpIn = NULL;
	FILE *fpOut = NULL;
	char strLine[LONG_LINE_LENGTH];
	double dRand;

	fpIn = fopen(strInputFile, "r");
	if(fpIn == NULL)
	{
		printf("Error: cannot open file!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = fopen(strOutputFile, "w");
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
		if(strLine[0] == '#')
			continue;

		dRand = rand_u();

		if(dRand <= dR)
			fprintf(fpOut, "%s\n", strLine);
	}


	fclose(fpIn);
	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;

}

int IQR_ProbeMatch()
{
	FILE *fpIn;
	FILE *fpIn2;
	FILE *fpOut;
	FILE *fpOut2;
	int ni,nj,nk,nz;
	struct tagString **vExon;
	struct INTMATRIX *pExonID;
	int nExonProbeNum = 0;
	char strLine[LONG_LINE_LENGTH];
	char strLine2[LONG_LINE_LENGTH];
	int nID;
	char strExonFile[LINE_LENGTH];
	char strU133File[LINE_LENGTH];
	char strMapFile[LINE_LENGTH];
	char strExonOutFile[LINE_LENGTH];
	char strU133OutFile[LINE_LENGTH];
	char strProbe[LINE_LENGTH];
	char strProbe2[LINE_LENGTH];
	char strProbeF[LINE_LENGTH];
	char strProbeF2[LINE_LENGTH];
	char *chp;
	int nEID;
	int n1end,n2end;
	int nCmp,nCmp1;

	/* load exon array data */
	strcpy(strExonFile, "G:\\Projects\\IQR_project\\summary.selected");
	strcpy(strU133File, "G:\\Projects\\IQR_project\\apt-probeset-summarize-results\\rma.summary.s.txt");
	strcpy(strMapFile, "G:\\Projects\\IQR_project\\HGU133Plus2_0-HuEx1_0st-transcript-cluster-mapping.txt");
	strcpy(strExonOutFile, "G:\\Projects\\IQR_project\\exonmatch.txt");
	strcpy(strU133OutFile, "G:\\Projects\\IQR_project\\u133match.txt");

	fpIn = NULL;
	fpIn = fopen(strExonFile, "r");
	if(fpIn == NULL)
	{
		printf("Error!\n");
		exit(EXIT_FAILURE);
	}

	fgets(strLine, LONG_LINE_LENGTH, fpIn);
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		nExonProbeNum++;
	}

	fclose(fpIn);

	pExonID = NULL;
	pExonID = CreateIntMatrix(1, nExonProbeNum);
	vExon = NULL;
	vExon = (struct tagString **)calloc(nExonProbeNum, sizeof(struct tagString *));
	if(pExonID == NULL || vExon == NULL)
	{
		printf("Error!\n");
		exit(EXIT_FAILURE);
	}

	fpIn = NULL;
	fpIn = fopen(strExonFile, "r");
	if(fpIn == NULL)
	{
		printf("Error!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	fgets(strLine, LONG_LINE_LENGTH, fpIn);
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;

		sscanf(strLine, "%d", &nID);
		pExonID->pMatElement[ni] = nID;
		StringAddTail(vExon+ni, strLine);

		ni++;
	}

	fclose(fpIn);


	/* match */
	fpIn = NULL;
	fpIn = fopen(strU133File, "r");
	if(fpIn == NULL)
	{
		printf("Error!\n");
		exit(EXIT_FAILURE);
	}

	fpIn2 = NULL;
	fpIn2 = fopen(strMapFile, "r");
	if(fpIn2 == NULL)
	{
		printf("Error!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strU133OutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error!\n");
		exit(EXIT_FAILURE);
	}

	fpOut2 = NULL;
	fpOut2 = fopen(strExonOutFile, "w");
	if(fpOut2 == NULL)
	{
		printf("Error!\n");
		exit(EXIT_FAILURE);
	}

	fgets(strLine, LONG_LINE_LENGTH, fpIn);
	fgets(strLine2, LONG_LINE_LENGTH, fpIn2);
	n1end = 0;
	n2end = 0;
	fgets(strLine, LONG_LINE_LENGTH, fpIn);
	fgets(strLine2, LONG_LINE_LENGTH, fpIn2);
	sscanf(strLine, "%s", strProbe);
	sscanf(strLine2, "%s %d", strProbe2, &nEID);
	while((n1end == 0) && (n2end == 0))
	{
		strcpy(strProbeF, strProbe);
		chp = strchr(strProbeF, '_');
		if(chp != NULL)
			*chp = '\0';

		strcpy(strProbeF2, strProbe2);
		chp = strchr(strProbeF2, '_');
		if(chp != NULL)
			*chp = '\0';

		nCmp1 = strcmp(strProbeF, strProbeF2);
		if(nCmp1 == 0)
		{
			nCmp = strcmp(strProbe, strProbe2);
		}
		else
		{
			nCmp = nCmp1;
		}


		/* nCmp = strcmp(strProbe, strProbe2); */
		if(nCmp == 0)
		{
			nj=0;
			nk=nExonProbeNum-1;
			while( (nk-nj)>1 )
			{
				nz = (int)((nk+nj)/2);
				if(pExonID->pMatElement[nz] > nEID)
				{
					nk = nz;
				}
				else
				{
					nj = nz;
				}
			}
			if(pExonID->pMatElement[nj] == nEID)
			{
				nz = nj;
				fprintf(fpOut, "%s", strLine);
				fprintf(fpOut2, "%s\n", vExon[nz]->m_pString);
			}
			else if(pExonID->pMatElement[nk] == nEID)
			{
				nz = nk;
				fprintf(fpOut, "%s", strLine);
				fprintf(fpOut2, "%s\n", vExon[nz]->m_pString);
			}

			if( fgets(strLine, LONG_LINE_LENGTH, fpIn) == NULL)
			{
				n1end = 1;
				break;
			}
			if( fgets(strLine2, LONG_LINE_LENGTH, fpIn2) == NULL)
			{
				n2end = 1;
				break;
			}

			sscanf(strLine, "%s", strProbe);
			sscanf(strLine2, "%s %d", strProbe2, &nEID);
		}
		else if(nCmp < 0)
		{
			if( fgets(strLine, LONG_LINE_LENGTH, fpIn) == NULL)
			{
				n1end = 1;
				break;
			}

			sscanf(strLine, "%s", strProbe);
		}
		else
		{
			if( fgets(strLine2, LONG_LINE_LENGTH, fpIn2) == NULL)
			{
				n2end = 1;
				break;
			}
			sscanf(strLine2, "%s %d", strProbe2, &nEID);
		}
	}


	fclose(fpIn);
	fclose(fpIn2);
	fclose(fpOut);
	fclose(fpOut2);

	/* release memory */
	DestroyIntMatrix(pExonID);
	for(ni=0; ni<nExonProbeNum; ni++)
	{
		DeleteString(vExon[ni]);
		vExon[ni] = NULL;
	}
	free(vExon);

	/* return */
	return PROC_SUCCESS;
}

int cMyc_AffyMatch()
{
	struct tagRefGene **vGene = NULL;
	int nGeneNum = 0;
	char strDatabasePath[LINE_LENGTH];
	char strSpecies[LINE_LENGTH];
	struct DOUBLEMATRIX *pMaxT;
	struct INTMATRIX *pLink;
	struct tagString **vData;
	int nDataNum;
	char strInFile[LINE_LENGTH];
	char strOutFile[LINE_LENGTH];
	char strLine[LONG_LINE_LENGTH];
	FILE *fpIn;
	FILE *fpOut;
	int ni,nj,nk,nx,ny,nz,nu;
	int nInfoNum1 = 5;
	int nInfoNum2 = 18;
	int nInfoNum3 = 5;
	char *chp,*chp2,*chp3;
	int nMaxLoc = 1024;
	int nLocusID[1024];
	char strRefID[1024][50];
	int nLocusNum,nRefNum;
	double dT;
	int nHit;


	/* init database */
	strcpy(strInFile, "G:\\Projects\\cMyc_Project\\analysis\\Expression\\MycRMALimma.txt");
	strcpy(strOutFile, "G:\\Projects\\cMyc_Project\\analysis\\Expression\\hg18_MycRMALimma.txt");
	strcpy(strDatabasePath, "G:\\Projects\\cMyc_Project\\analysis\\Expression\\refLocus_hg18.txt");
	strcpy(strSpecies, "human");

	vGene = RefGene_LoadDatabase(strDatabasePath, 2, strSpecies, &nGeneNum);
	pMaxT = NULL;
	pMaxT = CreateDoubleMatrix(1, nGeneNum);
	if(pMaxT == NULL)
	{
		printf("Error!\n");
		exit(EXIT_FAILURE);
	}

	pLink = NULL;
	pLink = CreateIntMatrix(1, nGeneNum);
	if(pLink == NULL)
	{
		printf("Error!\n");
		exit(EXIT_FAILURE);
	}


	/* load data */
	fpIn = NULL;
	fpIn = fopen(strInFile, "r");
	if(fpIn == NULL)
	{
		printf("Error!\n");
		exit(EXIT_FAILURE);
	}

	nDataNum = 0;
	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		nDataNum++;
	}

	fclose(fpIn);

	vData = NULL;
	vData = (struct tagString **)calloc(nDataNum, sizeof(struct tagString *));
	if(vData == NULL)
	{
		printf("Error!\n");
		exit(EXIT_FAILURE);
	}

	ni = 0;
	fpIn = NULL;
	fpIn = fopen(strInFile, "r");
	if(fpIn == NULL)
	{
		printf("Error!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
			continue;

		StringAddTail(vData+ni, strLine);

		chp = strchr(strLine, '\t');
		chp2 = chp+1;
		chp = strchr(chp2, '\t');
		*chp = '\0';

		nLocusNum = 0;
		chp3 = chp2;
		chp2 = chp+1;

		StrTrimLeft(chp3);
		StrTrimRight(chp3);
		if(strcmp(chp3, "---") == 0)
		{
			nLocusID[nLocusNum] = -1;
			nLocusNum = 1;
		}
		else
		{
			chp = strstr(chp3, "///");
			while(chp != NULL)
			{
				*chp = '\0';
				StrTrimLeft(chp3);
				StrTrimRight(chp3);
				nLocusID[nLocusNum] = atoi(chp3);
				chp3 = chp+3;
				chp = strstr(chp3, "///");
				nLocusNum++;
				if(nLocusNum >= nMaxLoc)
				{
					printf("Error!\n");
					exit(EXIT_FAILURE);
				}
			}
			StrTrimLeft(chp3);
			StrTrimRight(chp3);
			nLocusID[nLocusNum] = atoi(chp3);
			nLocusNum++;
			if(nLocusNum >= nMaxLoc)
			{
				printf("Error!\n");
				exit(EXIT_FAILURE);
			}
		}

		chp = strchr(chp2, '\t');
		chp2 = chp+1;
		chp = strchr(chp2, '\t');
		*chp = '\0';

		nRefNum = 0;
		chp3 = chp2;
		chp2 = chp+1;
		StrTrimLeft(chp3);
		StrTrimRight(chp3);
		if(strcmp(chp3, "---") == 0)
		{
			nRefNum = 0;
		}
		else
		{
			chp = strstr(chp3, "///");
			while(chp != NULL)
			{
				*chp = '\0';
				StrTrimLeft(chp3);
				StrTrimRight(chp3);
				strcpy(strRefID[nRefNum], chp3);
				chp3 = chp+3;
				chp = strstr(chp3, "///");
				nRefNum++;
				if(nRefNum >= nMaxLoc)
				{
					printf("Error!\n");
					exit(EXIT_FAILURE);
				}
			}
			StrTrimLeft(chp3);
			StrTrimRight(chp3);
			strcpy(strRefID[nRefNum], chp3);
			nRefNum++;
			if(nRefNum >= nMaxLoc)
			{
				printf("Error!\n");
				exit(EXIT_FAILURE);
			}
		}

		chp = strchr(chp2, '\t');
		chp2 = chp+1;
		chp = strchr(chp2, '\t');
		chp2 = chp+1;
		chp = strchr(chp2, '\t');
		chp2 = chp+1;
		chp = strchr(chp2, '\t');
		*chp = '\0';
		dT = fabs(atof(chp2));

		for(nj=0; nj<nLocusNum; nj++)
		{
			nx = 0;
			ny = nGeneNum-1;
			while(ny-nx > 1)
			{
				nz = (nx+ny)/2;
				if(nLocusID[nj] < vGene[nz]->nGeneID)
				{
					ny = nz;
				}
				else if(nLocusID[nj] > vGene[nz]->nGeneID)
				{
					nx = nz;
				}
				else
				{
					break;
				}
			}

			if(nLocusID[nj] != vGene[nz]->nGeneID)
			{
				if(nLocusID[nj] == vGene[nx]->nGeneID)
					nz = nx;
				else if(nLocusID[nj] == vGene[ny]->nGeneID)
					nz = ny;
				else
					continue;
			}

			nx = nz;
			ny = nz;
			while(nx>0)
			{
				if(nLocusID[nj] == vGene[nx-1]->nGeneID)
					nx--;
				else
					break;
			}

			while(ny<(nGeneNum-1))
			{
				if(nLocusID[nj] == vGene[ny+1]->nGeneID)
					ny++;
				else
					break;
			}

			for(nz=nx; nz<=ny; nz++)
			{
				if(nRefNum == 0)
					nHit = 0;
				else
				{
					nHit = 0;
					for(nu=0; nu<nRefNum; nu++)
					{
						if(strcmp(strRefID[nu], vGene[nz]->strName) == 0)
						{
							nHit = 1;
							break;
						}
					}
				}

				if(nHit == 1)
				{
					if(dT > pMaxT->pMatElement[nz])
					{
						pMaxT->pMatElement[nz] = dT;
						pLink->pMatElement[nz] = ni+1;
					}
				}
			}

		}

		ni++;
	}

	fclose(fpIn);
	if(ni != nDataNum)
	{
		printf("Error!\n");
		exit(EXIT_FAILURE);
	}

	/* export */
	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<nGeneNum; ni++)
	{
		nj = pLink->pMatElement[ni]-1;
		if(nj>=0)
		{
			fprintf(fpOut, "%s\n", vData[nj]->m_pString);
		}
		else
		{
			fprintf(fpOut, "---");
			for(nk=1; nk<nInfoNum1; nk++)
			{
				fprintf(fpOut, "\t---");
			}
			for(nk=0; nk<nInfoNum2; nk++)
			{
				fprintf(fpOut, "\t-1000");
			}
			for(nk=0; nk<nInfoNum3; nk++)
			{
				fprintf(fpOut, "\t---");
			}
			fprintf(fpOut, "\n");
		}
	}

	fclose(fpOut);

	/* clear memory */
	RefGene_ClearDatabase(&vGene, nGeneNum);
	for(ni=0; ni<nDataNum; ni++)
	{
		DeleteString(vData[ni]);
		vData[ni] = NULL;
	}
	free(vData);
	DestroyDoubleMatrix(pMaxT);
	DestroyIntMatrix(pLink);

	/* return */
	return PROC_SUCCESS;
}

int cMyc_AgilentMatch()
{
	struct tagString **strGeneLine;
	int nGeneNum = 0;
	char strRefID[255];
	char strCID[255];
	double dT,dMaxT;
	char strInFile[255];
	char strOutFile[255];
	char strSpecies[255];
	char strLine[LONG_LINE_LENGTH];
	char strDataLine[LONG_LINE_LENGTH];
	char strTemp[LONG_LINE_LENGTH];
	char strLine0[LONG_LINE_LENGTH];
	char *chp1,*chp2;
	int ni;
	FILE *fpIn,*fpOut;


	/* init database */
	strcpy(strInFile, "G:\\Projects\\cMyc_Project\\analysis\\Expression\\AgilentLimma_hg18presort.txt");
	strcpy(strOutFile, "G:\\Projects\\cMyc_Project\\analysis\\Expression\\hg18_AgilentLimma.txt");
	strcpy(strSpecies, "human");
	strcpy(strRefID, "---");
	strcpy(strDataLine, "---");

	strGeneLine = NULL;
	strGeneLine = (struct tagString **)calloc(1024, sizeof(struct tagString *));
	if(strGeneLine == NULL)
	{
		printf("Error!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<14; ni++)
	{
		strcpy(strTemp, strDataLine);
		sprintf(strDataLine, "%s\t-1000", strTemp);
	}

	/* load data */
	fpIn = NULL;
	fpIn = fopen(strInFile, "r");
	if(fpIn == NULL)
	{
		printf("Error!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		
		if(strstr(strLine, "###") == strLine)
		{
			chp1 = strchr(strLine, '\t');
			chp2 = chp1+1;
			chp1 = strchr(chp2, '\t');
			chp2 = chp1+1;
			strcpy(strLine0, chp2);
			chp1 = strchr(chp2, '\t');
			*chp1 = '\0';
			strcpy(strCID, chp2);
			for(ni=0; ni<5; ni++)
			{
				chp2 = chp1+1;
				chp1 = strchr(chp2, '\t');
			}
			*chp1 = '\0';
			dT = fabs(atof(chp2));

			if(strcmp(strCID, strRefID) != 0)
			{
				if(strcmp(strRefID, "---") != 0)
				{
					for(ni=0; ni<nGeneNum; ni++)
					{
						fprintf(fpOut, "%s\t%s\n", strGeneLine[ni]->m_pString, strDataLine);
						DeleteString(strGeneLine[ni]);
						strGeneLine[ni] = NULL;
					}
				}

				dMaxT = dT;
				nGeneNum = 0;
				strcpy(strRefID, strCID);
				strcpy(strDataLine, strLine0);
			}
			else
			{
				if(dT > dMaxT)
				{
					dMaxT = dT;
					strcpy(strDataLine, strLine0);
				}
			}
		}
		else
		{
			strcpy(strLine0, strLine);
			chp1 = strchr(strLine, '\t');
			chp2 = chp1+1;
			chp1 = strchr(chp2, '\t');
			chp2 = chp1+1;
			chp1 = strchr(chp2, '\t');
			*chp1 = '\0';
			strcpy(strCID, chp2);
			
			if(strcmp(strCID, strRefID) != 0)
			{
				if(strcmp(strRefID, "---") != 0)
				{
					for(ni=0; ni<nGeneNum; ni++)
					{
						fprintf(fpOut, "%s\t%s\n", strGeneLine[ni]->m_pString, strDataLine);
						DeleteString(strGeneLine[ni]);
						strGeneLine[ni] = NULL;
					}
					nGeneNum = 0;
				}

				strcpy(strDataLine, "---");
				for(ni=0; ni<14; ni++)
				{
					strcpy(strTemp, strDataLine);
					sprintf(strDataLine, "%s\t-1000", strTemp);
				}
				dMaxT = 0.0;
				StringAddTail(strGeneLine+nGeneNum, strLine0);
				nGeneNum++;
				if(nGeneNum >= 1024)
				{
					printf("Error!\n");
					exit(EXIT_FAILURE);
				}
				strcpy(strRefID, strCID);
			}
			else
			{
				StringAddTail(strGeneLine+nGeneNum, strLine0);
				nGeneNum++;
				if(nGeneNum >= 1024)
				{
					printf("Error!\n");
					exit(EXIT_FAILURE);
				}
			}

		}
		
	}

	if(strcmp(strRefID, "---") != 0)
	{
		for(ni=0; ni<nGeneNum; ni++)
		{
			fprintf(fpOut, "%s\t%s\n", strGeneLine[ni]->m_pString, strDataLine);
			DeleteString(strGeneLine[ni]);
			strGeneLine[ni] = NULL;
		}
	}

	fclose(fpIn);
	fclose(fpOut);
	free(strGeneLine);

	/* return */
	return PROC_SUCCESS;
}

int cMyc_ExonMatch()
{
	struct tagString **strGeneLine;
	int nGeneNum = 0;
	char strRefID[255];
	char strCID[255];
	double dT,dMaxT,dU;
	char strInFile[255];
	char strOutFile[255];
	char strSpecies[255];
	char strLine[LONG_LINE_LENGTH];
	char strDataLine[LONG_LINE_LENGTH];
	char strTemp[LONG_LINE_LENGTH];
	char strLine0[LONG_LINE_LENGTH];
	char *chp1,*chp2;
	int ni;
	FILE *fpIn,*fpOut;


	/* init database */
	dMaxT = -1000.0;
	strcpy(strInFile, "G:\\Projects\\cMyc_Project\\analysis\\ExonArray\\MycExonLimma_AllwAnnot_presort.txt");
	strcpy(strOutFile, "G:\\Projects\\cMyc_Project\\analysis\\ExonArray\\hg18_ExonLimma_All.txt");
	strcpy(strSpecies, "human");
	strcpy(strRefID, "---");
	strcpy(strDataLine, "---\t---");

	strGeneLine = NULL;
	strGeneLine = (struct tagString **)calloc(1024, sizeof(struct tagString *));
	if(strGeneLine == NULL)
	{
		printf("Error!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<7; ni++)
	{
		strcpy(strTemp, strDataLine);
		sprintf(strDataLine, "%s\t-1000", strTemp);
	}

	/* load data */
	fpIn = NULL;
	fpIn = fopen(strInFile, "r");
	if(fpIn == NULL)
	{
		printf("Error!\n");
		exit(EXIT_FAILURE);
	}

	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error!\n");
		exit(EXIT_FAILURE);
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn)!=NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		
		if(strstr(strLine, "###") == strLine)
		{
			chp1 = strchr(strLine, '\t');
			chp2 = chp1+1;
			strcpy(strLine0, chp2);
			chp1 = strchr(chp2, '\t');
			chp2 = chp1+1;
			chp1 = strchr(chp2, '\t');
			*chp1 = '\0';
			strcpy(strCID, chp2);
			for(ni=0; ni<4; ni++)
			{
				chp2 = chp1+1;
				chp1 = strchr(chp2, '\t');
			}
			*chp1 = '\0';
			dT = atof(chp2);

			chp2 = chp1+1;
			chp1 = strchr(chp2, '\t');
			*chp1 = '\0';
			dU = atof(chp2);

			if(dU > -500)
				dT = fabs(dT);
			else
				dT = -500.0;
							

			if(strcmp(strCID, strRefID) != 0)
			{
				if(strcmp(strRefID, "---") != 0)
				{
					for(ni=0; ni<nGeneNum; ni++)
					{
						fprintf(fpOut, "%s\t%s\n", strGeneLine[ni]->m_pString, strDataLine);
						DeleteString(strGeneLine[ni]);
						strGeneLine[ni] = NULL;
					}
				}

				dMaxT = dT;
				nGeneNum = 0;
				strcpy(strRefID, strCID);
				strcpy(strDataLine, strLine0);
			}
			else
			{
				if(dT > dMaxT)
				{
					dMaxT = dT;
					strcpy(strDataLine, strLine0);
				}
			}
		}
		else
		{
			strcpy(strLine0, strLine);
			chp1 = strchr(strLine, '\t');
			chp2 = chp1+1;
			chp1 = strchr(chp2, '\t');
			chp2 = chp1+1;
			chp1 = strchr(chp2, '\t');
			*chp1 = '\0';
			strcpy(strCID, chp2);
			
			if(strcmp(strCID, strRefID) != 0)
			{
				if(strcmp(strRefID, "---") != 0)
				{
					for(ni=0; ni<nGeneNum; ni++)
					{
						fprintf(fpOut, "%s\t%s\n", strGeneLine[ni]->m_pString, strDataLine);
						DeleteString(strGeneLine[ni]);
						strGeneLine[ni] = NULL;
					}
					nGeneNum = 0;
				}

				strcpy(strDataLine, "---\t---");
				for(ni=0; ni<7; ni++)
				{
					strcpy(strTemp, strDataLine);
					sprintf(strDataLine, "%s\t-1000", strTemp);
				}
				dMaxT = -1000.0;
				StringAddTail(strGeneLine+nGeneNum, strLine0);
				nGeneNum++;
				if(nGeneNum >= 1024)
				{
					printf("Error!\n");
					exit(EXIT_FAILURE);
				}
				strcpy(strRefID, strCID);
			}
			else
			{
				StringAddTail(strGeneLine+nGeneNum, strLine0);
				nGeneNum++;
				if(nGeneNum >= 1024)
				{
					printf("Error!\n");
					exit(EXIT_FAILURE);
				}
			}

		}
		
	}

	if(strcmp(strRefID, "---") != 0)
	{
		for(ni=0; ni<nGeneNum; ni++)
		{
			fprintf(fpOut, "%s\t%s\n", strGeneLine[ni]->m_pString, strDataLine);
			DeleteString(strGeneLine[ni]);
			strGeneLine[ni] = NULL;
		}
	}

	fclose(fpIn);
	fclose(fpOut);
	free(strGeneLine);

	/* return */
	return PROC_SUCCESS;
}


int book_norm()
{
	char strFileName[4][LINE_LENGTH];
	char strOutFile[LINE_LENGTH];
	char strTag[4][LINE_LENGTH];
	struct tagBARData *pBARData;
	FILE *fpOut;
	int ni;
	int nj,nk,nl;

	/* init */
	
	/* 
	strcpy(strOutFile, "G:\\Projects\\temp\\Raw.txt");
	strcpy(strFileName[0], "G:\\Projects\\temp\\060504_ChIP_GliAct_EB.CEL.bar");
	strcpy(strFileName[1], "G:\\Projects\\temp\\060516_ChIP_Gli3T_5_10.CEL.bar");
	strcpy(strFileName[2], "G:\\Projects\\temp\\060504_Input_EB.CEL.bar");
	strcpy(strFileName[3], "G:\\Projects\\temp\\060516_Input_Gli3T_5_11.CEL.bar");
	*/

	/* 
	strcpy(strOutFile, "G:\\Projects\\temp\\QN.txt");
	strcpy(strFileName[0], "D:\\Projects\\hedgehog_project\\ChIP-chip\\Gli_Promoter\\analysis_mm8\\step0_signal\\060504_ChIP_GliAct_EB.CEL.bar");
	strcpy(strFileName[1], "D:\\Projects\\hedgehog_project\\ChIP-chip\\Gli_Promoter\\analysis_mm8\\step0_signal\\060516_ChIP_Gli3T_5_10.CEL.bar");
	strcpy(strFileName[2], "D:\\Projects\\hedgehog_project\\ChIP-chip\\Gli_Promoter\\analysis_mm8\\step0_signal\\060504_Input_EB.CEL.bar");
	strcpy(strFileName[3], "D:\\Projects\\hedgehog_project\\ChIP-chip\\Gli_Promoter\\analysis_mm8\\step0_signal\\060516_Input_Gli3T_5_11.CEL.bar");
	*/

	strcpy(strOutFile, "G:\\Projects\\temp\\MAT.txt");
	strcpy(strFileName[0], "G:\\Projects\\TileProbe_Project\\mymat\\060504_ChIP_GliAct_EB.bar");
	strcpy(strFileName[1], "G:\\Projects\\TileProbe_Project\\mymat\\060516_ChIP_Gli3T_5_10.bar");
	strcpy(strFileName[2], "G:\\Projects\\TileProbe_Project\\mymat\\060504_Input_EB.bar");
	strcpy(strFileName[3], "G:\\Projects\\TileProbe_Project\\mymat\\060516_Input_Gli3T_5_11.bar");

	strcpy(strTag[0], "IP1");
	strcpy(strTag[1], "IP2");
	strcpy(strTag[2], "CT1");
	strcpy(strTag[3], "CT2");

	/* process */
	fpOut = NULL;
	fpOut = fopen(strOutFile, "w");
	if(fpOut == NULL)
	{
		printf("Error: cannot open output file!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<4; ni++)
	{
		pBARData = NULL;
		pBARData = Affy_LoadBAR_Fast(strFileName[ni]);
		if(pBARData == NULL)
		{
			printf("Error: cannot load bar data!\n");
			exit(EXIT_FAILURE);
		}

		nl = 0;
		for(nj=0; nj<pBARData->nSeqNum; nj++)
		{
			for(nk=0; nk<pBARData->vSeqData[nj]->nDataNum; nk++)
			{
				fprintf(fpOut, "%f\t%s\n", pBARData->vSeqData[nj]->vData[1]->pMatElement[nk], strTag[ni]);
				nl++;
				if(nl == 10000)
					break;
			}

			if(nl == 10000)
				break;
		}

		Affy_BARData_Destroy(&pBARData);
	}

	fclose(fpOut);

	/* return */
	return PROC_SUCCESS;
}

int shhmbtest(int argv, char **argc)
{
	/* define */
	struct INTMATRIX *pIP;
	struct INTMATRIX *pCT;
	char strIPFile[LINE_LENGTH];
	char strCTFile[LINE_LENGTH];
	char strSampleFile[LINE_LENGTH];
	char strLine[MED_LINE_LENGTH];
	char strArgFile[LINE_LENGTH];
	int ni,nj,nk,nx,ny;
	char strSampleName[6][LINE_LENGTH];
	char strTitle[LINE_LENGTH];
	FILE *fpIn;
	FILE *fpOut;
	int nResult;

	/* init */
	if(argv != 5)
	{
		printf("Parameter error!\n");
		exit(EXIT_FAILURE);
	}
	strcpy(strSampleFile, argc[1]);
	strcpy(strIPFile, argc[2]);
	strcpy(strCTFile, argc[3]);
	strcpy(strTitle, argc[4]);

	/* load data */
	pIP = NULL;
	pIP = IMLOAD(strIPFile);
	if(pIP == NULL)
	{
		printf("Error!\n");
		exit(EXIT_FAILURE);
	}
	pCT = NULL;
	pCT = IMLOAD(strCTFile);
	if(pCT == NULL)
	{
		printf("Error!\n");
		exit(EXIT_FAILURE);
	}

	fpIn = NULL;
	fpIn = fopen(strSampleFile, "r");
	if(fpIn == NULL)
	{
		printf("Error!\n");
		exit(EXIT_FAILURE);
	}

	for(ni=0; ni<6; ni++)
	{
		fgets(strLine, MED_LINE_LENGTH, fpIn);
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		strcpy(strSampleName[ni], strLine);
	}

	fclose(fpIn);

	/* process */
	nk = 0;
	for(ni=0; ni<pIP->nHeight; ni++)
	{
		for(nj=0; nj<pCT->nHeight; nj++)
		{
			/* write argument file */
			sprintf(strArgFile, "%s_arg.txt", strTitle);
			fpOut = NULL;
			fpOut = fopen(strArgFile, "w");
			if(fpOut == NULL)
			{
				printf("Error!\n");
				exit(EXIT_FAILURE);
			}

			fprintf(fpOut, "[Comparison Type] (1: One Sample; 2: Two Sample; 3: Multiple Sample) = 2\n");
			fprintf(fpOut, "[Working Directory] = .\n");
			fprintf(fpOut, "[Project Title] = %s_%d\n", strTitle, nk);
			fprintf(fpOut, "[No. of Libraries] = 1\n");
			fprintf(fpOut, "[No. of Samples] = %d\n", pIP->nWidth+pCT->nWidth);
			fprintf(fpOut, "[No. of Groups] = 2\n");
			fprintf(fpOut, "[Data]\n");
			for(nx=0; nx<pIP->nWidth; nx++)
			{
				fprintf(fpOut, "1->IP%d\n", nx+1);
				ny = IMGETAT(pIP, ni, nx)-1;
				fprintf(fpOut, "%s\n", strSampleName[ny]);
			}
			for(nx=0; nx<pCT->nWidth; nx++)
			{
				fprintf(fpOut, "2->CT%d\n", nx+1);
				ny = IMGETAT(pCT, nj, nx)-1;
				fprintf(fpOut, "%s\n", strSampleName[ny]);
			}
			fprintf(fpOut, "[Patterns of Interest]\n");
			fprintf(fpOut, "1>2\n");
			fprintf(fpOut, "[Masking Bad Data Points] (0:No, 1:Yes) = 0\n");
			fprintf(fpOut, "[Truncation Lower Bound] = -1000000.000000\n");
			fprintf(fpOut, "[Truncation Upper Bound] = 1000000.000000\n");
			fprintf(fpOut, "[Transformation] (0: Identity; 1: log2; 2: logit; 3: exp(x)/1+exp(x)) = 0\n");
			fprintf(fpOut, "[Common Variance Groups] = 1\n");
			fprintf(fpOut, "1 2\n");
			fprintf(fpOut, "[Method to Combine Neighboring Probes] (0:HMM, 1:MA) = 1\n");
			fprintf(fpOut, "[Method to Compute FDR] (0: Estimate from Left Tail; 1: Permutation Test; 2: UMS; 3: No FDR) = 0\n");
			fprintf(fpOut, "[W] = 5\n");
			fprintf(fpOut, "[Window Boundary] = 300\n");
			fprintf(fpOut, "[Standardize MA Statistics] (0:No; 1:Yes) = 1\n");
			fprintf(fpOut, "[Region Boundary Cutoff, MA>] = 2.000000\n");
			fprintf(fpOut, "[Expected Hybridization Length] = 28\n");
			fprintf(fpOut, "[Posterior Probability Cutoff, P>] = 0.500000\n");
			fprintf(fpOut, "[G0 Selection Criteria, p%] = 0.010000\n");
			fprintf(fpOut, "[G1 Selection Criteria, q%] = 0.050000\n");
			fprintf(fpOut, "[Selection Offset] = 6\n");
			fprintf(fpOut, "[Grid Size] = 1000\n");
			fprintf(fpOut, "[Number of Permutations] = 10\n");
			fprintf(fpOut, "[Exchangeable Groups] = 1\n");
			fprintf(fpOut, "1 2\n");
			fprintf(fpOut, "[Max Gap within a Region] = 300\n");
			fprintf(fpOut, "[Max Run of Insignificant Probes within a Region] = 5\n");
			fprintf(fpOut, "[Min Region Length] = 100\n");
			fprintf(fpOut, "[Min No. of Significant Probes within a Region] = 5\n");

			/* close file */
			fclose(fpOut);

			/* run tilemap */
			sprintf(strLine, "/home/bst/faculty/hji/projects/cisgenome_project/bin/tilemapv2 %s", strArgFile);
			nResult = system(strLine);
			nk++;
		}
	}

	/* release memory */
	DestroyIntMatrix(pIP);
	DestroyIntMatrix(pCT);

	/* return */
	return PROC_SUCCESS;
}