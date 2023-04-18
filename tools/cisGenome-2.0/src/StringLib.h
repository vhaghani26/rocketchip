/* ----------------------------------------------------------------------- */
/*                                 Macro                                   */
/* ----------------------------------------------------------------------- */
#define LINE_LENGTH 255
#define MED_LINE_LENGTH 1024
#define MEDLONG_LINE_LENGTH 20000
#define LONG_LINE_LENGTH 65535
#define PROC_SUCCESS 1
#define PROC_FAILURE 0
#define WORD_SEPARATORS "\t "

/* ----------------------------------------------------------------------- */ 
/*                                 Struct                                  */
/* ----------------------------------------------------------------------- */ 
/* tagString is used for string processing. */
struct tagString
{
	/* pString is used for recording string. */
	char *m_pString;
	/* m_nLength is used for recording the length of string. */
	long m_nLength;
};

/* tagString is used for linking a pair of strings. */
struct tagStringPair
{
	/* m_pStr is used for recording string. */
	struct tagString *m_pStr1;
	struct tagString *m_pStr2;
	/* m_pNext is used for constructing linear list */
	struct tagStringPair *m_pNext;
};

/* tagWString is used for wchar_t string processing. */
struct tagWString
{
	/* pWString is used for recording string. */
	wchar_t *m_pWString;
	/* m_nLength is used for recording the length of string. */
	long m_nLength;
};

/* ----------------------------------------------------------------------- */
/*                              Declaration                                */
/* ----------------------------------------------------------------------- */
char* StrTrimRight(char strLongLine[]);
char* StrTrimLeft(char strLongLine[]);
char* StrTrimSpace(char strLongLine[]);
void StrMakeUpper(char strLongLine[]);
void StrMakeLower(char strLongLine[]);
struct tagString *CreateString(long nLength);
void DeleteString(struct tagString *pDelString);
struct tagString *StringAddTail(struct tagString **pString, char strLine[]);

struct tagWString *CreateWString(long nLength);
void DeleteWString(struct tagWString **pDelWString);


/* ----------------------------------------------------------------------------- */ 
/*                   struct tagStringPair *StringPairCreate()                    */
/* This function is used for creating a string pair.                             */
/* If the process is successful, it will return the pointer to the string pair.  */
/* else it will return NULL.                                                     */
/* ----------------------------------------------------------------------------- */ 
struct tagStringPair *StringPairCreate();

/* ----------------------------------------------------------------------------- */ 
/*                             void StringPairDestroy()                          */
/* This function is used for destroy a string pair.                              */
/* ----------------------------------------------------------------------------- */ 
void StringPairDestroy(struct tagStringPair **pStringPair);

/* ----------------------------------------------------------------------------- */ 
/*                           void StringPairClearList()                          */
/* This function is used for destroy a list of string pairs.                     */
/* ----------------------------------------------------------------------------- */ 
void StringPairClearList(struct tagStringPair **pStringPairList);

/* ----------------------------------------------------------------------- */ 
/*  StringPair_ClearDatabase:                                              */
/*  clear StringPair database.                                             */
/* ----------------------------------------------------------------------- */ 
int StringPair_ClearDatabase(struct tagStringPair ***vSourceStringPair, int nPairNum);

