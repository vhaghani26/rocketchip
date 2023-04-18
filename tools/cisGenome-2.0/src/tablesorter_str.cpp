#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

#define LONG_LINE_LENGTH 65535

using namespace std;

bool tablefile_sort(string infile, string outfile);
char* StrTrimLeft(char strLongLine[]);
char* StrTrimRight(char strLongLine[]);

class table_item
{
public:
	char* chr;
	int coord;
	char* intensities;

};

bool itemsmaller(const table_item &item1, const table_item &item2)
{
   if(string(item1.chr) < string(item2.chr))
	   return true;

	if(string(item1.chr) == string(item2.chr) && item1.coord < item2.coord)
		return true;

	return false;
}

int main(int argc, char *argv[])
{
	string infile;
	string outfile;

	if (argc == 2)
	{
		infile = string(argv[1]);
		outfile = string(argv[1])+".sort";
	}
	else if (argc == 3) 
	{
		infile = string(argv[1]);
		outfile = string(argv[2]);
	}
	else
	{
		cout << "tablesorter <src txt file>" << endl;
		cout << "tablesorter <src txt file> <dest txt file>" << endl;
		return 1;

	}

	tablefile_sort(infile, outfile);

	return 0;
}

bool tablefile_sort(string infile, string outfile)
{
	vector<table_item> table_data;
	char strLine[LONG_LINE_LENGTH];
	int i;
	int line_num = 0;
	int field_num = -1;
	int data_start = 0;
	char *chp1,*chp2;

	FILE *fpIn = NULL;
	FILE *fpOut = NULL;

	fpIn = fopen(infile.c_str(), "r");
	if(fpIn == NULL)
	{
		cout << "Error: cannot open input file file!" << endl;
		return false;
	}

	fpOut = fopen(outfile.c_str(), "w");
	if(fpOut == NULL)
	{
		cout << "Error: cannot open output file!" << endl;
		return false;
	}

	while(fgets(strLine, LONG_LINE_LENGTH, fpIn) != NULL)
	{
		StrTrimLeft(strLine);
		StrTrimRight(strLine);
		if(strLine[0] == '\0')
			continue;
		if(strLine[0] == '#')
		{
			if(data_start == 0)
				fprintf(fpOut, "%s\n", strLine);
			continue;
		}

		data_start = 1;

		table_item newitem;
		newitem.coord = 0;
		string intensities = "";

		i = 0;
		chp1 = strLine;
		chp2 = strpbrk(chp1, " \t");
		while(chp2 != NULL)
		{
			*chp2 = '\0';
			if(i == 0)
			{
				newitem.chr = new char[strlen(chp1)];
				memcpy(newitem.chr, chp1, strlen(chp1)+1);
			}
			else if(i == 1)
			{
				newitem.coord = atoi(chp1);
			}
			else
			{
				intensities += string(chp1) + "\t";
			}

			i++;

			chp1 = chp2+1;
			while( *chp1 == ' ' || *chp1 == '\t')
				chp1++;
			chp2 = strpbrk(chp1, " \t");
		}

		if(i == 0)
		{
			newitem.chr = new char[strlen(chp1)];
			memcpy(newitem.chr, chp1, strlen(chp1)+1);
		}
		else if(i == 1)
		{
			newitem.coord = atoi(chp1);
		}
		else
		{
			intensities += string(chp1);
		}
		i++;

		newitem.intensities = new char[intensities.length()];
		memcpy(newitem.intensities, intensities.c_str(), intensities.length()+1);


		if(field_num < 0)
			field_num = i;
		else if(field_num != i)
		{
			cout << "Error: bad source file, column number not consistent!" << endl;
			table_data.clear();
			return false;
		}

		table_data.push_back(newitem);
		line_num++;
	}

	field_num -= 2;
	printf("%d data units loaded.\n", table_data.size());
	printf("Sorting ...\n");
	sort(table_data.begin(), table_data.end(), itemsmaller);
	for(i=0; i<line_num; i++)
	{
		fprintf(fpOut, "%s", table_data[i].chr);

		if(field_num >= 0)
			fprintf(fpOut, "\t%d", table_data[i].coord);

		fprintf(fpOut, "\t%s", table_data[i].intensities);

		fprintf(fpOut, "\n");
	}

	/* close file */
	fclose(fpIn);
	fclose(fpOut);
	
	for (i=0; i < (int)table_data.size(); i++) {
		delete[] table_data[i].chr;
		delete[] table_data[i].intensities;
	}
	table_data.clear();

	/* return */
	return true;
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

