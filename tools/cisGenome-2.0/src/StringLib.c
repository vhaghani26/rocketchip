/* ----------------------------------------------------------------------- */
/*  stringprocessing.c : process strings.                                  */
/*  Author : Ji HongKai ; Time: 2001.12                                    */
/* ----------------------------------------------------------------------- */

/* ----------------------------------------------------------------------- */ 
/*                                Include                                  */
/* ----------------------------------------------------------------------- */ 
#include "stdio.h" 
#include "stdlib.h" 
#include "string.h"
#include "limits.h"

#include "StringLib.h"

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

/* ----------------------------------------------------------------------- */ 
/*                         char* StrTrimSpace()                            */
/* This function is used for deleting ' '\t' within a string.              */
/* ----------------------------------------------------------------------- */
char* StrTrimSpace(char strLongLine[])
{
	char *startp;
	char separators[]="\t ";
	long nLen;

	if(strLongLine != NULL)
	{
		startp = strLongLine;
		while((startp = strpbrk(startp, separators)) != NULL)
		{
			nLen = strlen(startp);
			memmove(startp, startp+1, nLen);
		}
	}

	return strLongLine;
}


/* ----------------------------------------------------------------------- */ 
/*                         void StrMakeUpper()                             */
/* This function is used for converting lower case to upper case.          */
/* ----------------------------------------------------------------------- */ 
void StrMakeUpper(char strLongLine[])
{
	char *startp;

	startp = strLongLine;
	while(*startp != '\0')
	{
		if((*startp >= 'a')&&(*startp <= 'z'))
			*startp = *startp-32;
		startp++;
	}
}

/* ----------------------------------------------------------------------- */ 
/*                         void StrMakeLower()                             */
/* This function is used for converting upper case to lower case.          */
/* ----------------------------------------------------------------------- */ 
void StrMakeLower(char strLongLine[])
{
	char *startp;

	startp = strLongLine;
	while(*startp != '\0')
	{
		if((*startp >= 'A')&&(*startp <= 'Z'))
			*startp = *startp+32;
		startp++;
	}
}


/* ----------------------------------------------------------------------- */ 
/*                    struct tagString *CreateString()                     */
/* This function is used for creating a string.                            */
/* ----------------------------------------------------------------------- */ 
struct tagString *CreateString(long nLength)
{
	struct tagString *pNewString;

	/* Init */
	if(nLength < 0)
		return NULL;
	pNewString = NULL;

	/* new */
	pNewString = (struct tagString *)calloc(1, sizeof(struct tagString));
	if(pNewString == NULL)
		return NULL;
	pNewString->m_nLength = nLength;
	pNewString->m_pString = (char *)calloc((nLength+1), sizeof(char));
	if(pNewString->m_pString == NULL)
	{
		free(pNewString);
		return NULL;
	}
	*(pNewString->m_pString) = '\0';

	/* return */
	return pNewString;
}


/* ----------------------------------------------------------------------- */ 
/*                          void DeleteString()                            */
/* This function is used for deleting a string.                            */
/* ----------------------------------------------------------------------- */ 
void DeleteString(struct tagString *pDelString)
{
	char *pDelPoint;

	if( pDelString != NULL )
	{
		if( pDelString->m_pString != NULL )
		{
			pDelPoint = pDelString->m_pString;
			pDelString->m_pString = NULL;
			pDelPoint = realloc(pDelPoint, sizeof(char));
			free(pDelPoint);
		}
		free(pDelString);
	}	
}

/* ----------------------------------------------------------------------------- */ 
/*                    struct tagString *StringAddTail()                          */
/* This function is used for adding a segment to the tail of a string.           */
/* If the process is successful, it will return the pointer to the string;       */
/* else it will return NULL.                                                     */            
/* ----------------------------------------------------------------------------- */ 
struct tagString *StringAddTail(struct tagString **pString, char strLine[])
{
	/* nLen is used for recording the length of a string */
	int nLen;
	/* pString is used for manipulating the string */
	char *pvec,*pvectemp;

	/* check parameter */
	nLen = (int)strlen(strLine);
	if(nLen == 0)
		return *pString;

	if(*pString == NULL)
	{
		*pString = CreateString(nLen);
		if(*pString == NULL)
		{
			printf("Error in StringAddTail: cannot create a new string...\n");
			return NULL;
		}
		strcpy((*pString)->m_pString, strLine);
	}
	else
	{
		/* add */
		pvec = (*pString)->m_pString;
		pvectemp = (char *)malloc(sizeof(char)*((*pString)->m_nLength+nLen+1));
		if(pvectemp == NULL)
		{
			printf("Error in StringAddTail: cannot create a new string...\n");
			return NULL;
		}
		strcpy(pvectemp, pvec);
		strcat(pvectemp, strLine);
		(*pString)->m_pString = pvectemp;
		free(pvec);
		(*pString)->m_nLength += nLen;
	}
	
	/* return */
	return *pString;
}

/* ----------------------------------------------------------------------------- */ 
/*                   struct tagStringPair *StringPairCreate()                    */
/* This function is used for creating a string pair.                             */
/* If the process is successful, it will return the pointer to the string pair.  */
/* else it will return NULL.                                                     */
/* ----------------------------------------------------------------------------- */ 
struct tagStringPair *StringPairCreate()
{
	/* define */
	struct tagStringPair *pNewPair;

	/* create */
	pNewPair = NULL;
	pNewPair = (struct tagStringPair *)calloc(1, sizeof(struct tagStringPair));
	if(pNewPair == NULL)
	{
		printf("Error: cannot create string pair!\n");
		return NULL;
	}

	pNewPair->m_pNext = NULL;
	pNewPair->m_pStr1 = NULL;
	pNewPair->m_pStr2 = NULL;

	/* return */
	return pNewPair;
}

/* ----------------------------------------------------------------------------- */ 
/*                             void StringPairDestroy()                          */
/* This function is used for destroy a string pair.                              */
/* ----------------------------------------------------------------------------- */ 
void StringPairDestroy(struct tagStringPair **pStringPair)
{
	if(pStringPair != NULL)
	{
		if(*pStringPair != NULL)
		{
			DeleteString( (*pStringPair)->m_pStr1 );
			DeleteString( (*pStringPair)->m_pStr2 );
			free(*pStringPair);
			*pStringPair = NULL;
		}
	}
}

/* ----------------------------------------------------------------------------- */ 
/*                           void StringPairClearList()                          */
/* This function is used for destroy a list of string pairs.                     */
/* ----------------------------------------------------------------------------- */ 
void StringPairClearList(struct tagStringPair **pStringPairList)
{
	/* define */
	struct tagStringPair *pPair;
	
	if(pStringPairList != NULL)
	{
		while(*pStringPairList != NULL)
		{
			pPair = *pStringPairList;
			*pStringPairList = pPair->m_pNext;
			pPair->m_pNext = NULL;
			StringPairDestroy(&pPair);
		}
	}
}

/* ----------------------------------------------------------------------- */ 
/*  StringPair_ClearDatabase:                                              */
/*  clear StringPair database.                                             */
/* ----------------------------------------------------------------------- */ 
int StringPair_ClearDatabase(struct tagStringPair ***vSourceStringPair, int nPairNum)
{
	/* define */
	int ni;
	struct tagStringPair **vDatabase;

	/* clear */
	vDatabase = *vSourceStringPair;
	for(ni=0; ni<nPairNum; ni++)
	{
		StringPairDestroy(vDatabase+ni);
	}
	free(*vSourceStringPair);
	*vSourceStringPair = NULL;

	/* return */
	return PROC_SUCCESS;
}

/* ----------------------------------------------------------------------- */ 
/*                   struct tagWString *CreateWString()                    */
/* This function is used for creating a wchar_t string.                    */
/* ----------------------------------------------------------------------- */ 
struct tagWString *CreateWString(long nLength)
{
	struct tagWString *pNewString;

	/* Init */
	if(nLength < 0)
		return NULL;
	pNewString = NULL;

	/* new */
	pNewString = (struct tagWString *)calloc(1, sizeof(struct tagWString));
	if(pNewString == NULL)
		return NULL;
	pNewString->m_nLength = nLength;
	pNewString->m_pWString = (wchar_t *)calloc(nLength+1, sizeof(wchar_t));
	if(pNewString->m_pWString == NULL)
	{
		free(pNewString);
		return NULL;
	}
	
	/* return */
	return pNewString;
}


/* ----------------------------------------------------------------------- */ 
/*                         void DeleteWString()                            */
/* This function is used for deleting a wchar_t string.                    */
/* ----------------------------------------------------------------------- */ 
void DeleteWString(struct tagWString **pDelWString)
{
	wchar_t *pDelPoint;

	if(pDelWString == NULL)
		return;

	if( *pDelWString != NULL )
	{
		if( (*pDelWString)->m_pWString != NULL )
		{
			pDelPoint = (*pDelWString)->m_pWString;
			(*pDelWString)->m_pWString = NULL;
			free(pDelPoint);
		}
		free(*pDelWString);
		*pDelWString = NULL;
	}	
}