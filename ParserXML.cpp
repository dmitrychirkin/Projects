// ParserXML.cpp : Defines the entry point for the console application.
//

#include <Windows.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <io.h>

using namespace std;

struct Item
{
	string tagName;
	vector<string> tagValue;
};

class ParserFile
{
	public:
		ParserFile();
		~ParserFile();
		int parseFile(string fileNameIn, string filenameOut);
	protected:
		vector <string> inStr;
		vector <string> outStr;
		map <int, Item> tagmultiMap;
		int parseString(int numBer, vector<string> inStr);
};

class ListFolder
{
protected:
	map<string, vector<string>>	fileNames;

	vector<string> fillFolderTree(string folderName);

public:
	ListFolder();
	~ListFolder();

	map<string,vector<string>> getFolderContent(string folderName);
};

ListFolder::ListFolder()
{
	fileNames.clear();
}

ListFolder::~ListFolder()
{
	fileNames.clear();
}

vector<string> ListFolder::fillFolderTree(string folderName)
{
	_finddata32_t findData;
	long int hFile = 0;
	string mask = ".xml";
	string fileName = folderName+"/*.*";
	vector <string> filesInDir;
	int files_count = 1;

	try
	{
		memset(&findData,0,sizeof(findData));
		hFile = _findfirst32(fileName.c_str(),&findData);
		if ( !hFile )
			throw - 1;
		do
		{
			fileName = findData.name;
			if(fileName == "." || fileName == "..")
			{
			}
			else
			{
				if(findData.attrib == _A_SUBDIR)
				{
					vector <string> filesList = fillFolderTree(folderName+"/"+fileName);
					fileNames.insert( pair<string, vector<string>>(folderName+"/"+fileName,filesList) );
				}
				else if(fileName.rfind(mask) != string::npos)
					filesInDir.push_back(fileName);
			}
		} while(_findnext32(hFile,&findData) == 0);

	}
	catch(...)
	{
	}
	if(hFile)
		_findclose(hFile);
	return filesInDir;
}

map<string,vector<string>> ListFolder::getFolderContent(string folderName)
{
	fileNames.clear();
	vector <string> filesList = fillFolderTree(folderName);
	fileNames.insert( pair<string, vector<string>>(folderName,filesList) );
	return fileNames;
}

int main(int argc, char* argv[])
{
	string filenameIn = "";
	string filenameIn8 = "";
	string filenameOut = "";
	ParserFile m_parse;

	ListFolder m_listFolder;
	map<string,vector<string>> fileName = m_listFolder.getFolderContent("H:/Projects");

	map<string,vector<string>>::iterator mapiterator = fileName.begin();
	int count_of_elements = fileName.size();

	for(int i = 0; i < count_of_elements; i++)
	{
		printf("%s\n",mapiterator->first.c_str());
		for(int j = 0; j < mapiterator->second.size(); j++)
		{
			printf("\t%s\n", mapiterator->second[j].c_str());
		}
		mapiterator++;
	}

	filenameIn8 = "H:/Projects/test_utf8.xml";//"D:\\Authorise_1\\ParserXML\\test.xml";
	filenameOut = "H:/Projects/test_utf8.csv";//"D:\\Authorise_1\\ParserXML\\test.csv";
	m_parse.parseFile(filenameIn8, filenameOut);
	return 0;
}

ParserFile::ParserFile()
{
}

ParserFile::~ParserFile()
{
}

int ParserFile::parseFile(string fileNameIn, string filenameOut)
{
	ifstream m_fileIn;
	ofstream m_fileOut;
	FILE* m_ff = 0;
	string temp;

	try
	{
		m_fileIn.open(fileNameIn.c_str(),std::ios_base::in);
		if( !m_fileIn.is_open() )
			throw fileNameIn.c_str();
		inStr.resize(0);
		while ( getline(m_fileIn,temp) )
		{
			inStr.push_back(temp);
		}
		
		parseString(0, inStr);

		map<int, Item>::iterator mapiterator = tagmultiMap.begin();
		int count_of_elements = tagmultiMap.size();

		m_ff = fopen(filenameOut.c_str(),"w");
		for(int i = 0; i < count_of_elements+1; i++)
			for(int j = 0; j < 16; j++)
				fprintf(m_ff,"-");
		fprintf(m_ff,"\n|\t");
		for( int i = 0; i < count_of_elements; i++)
		{
			fprintf(m_ff,"%-16s|",tagmultiMap[i].tagName.c_str());
//			cout << "\t" << tagmultiMap[i].tagName ;
		}
		fprintf(m_ff,"\n|");
		for(int i = 0; i < count_of_elements+1; i++)
			for(int j = 0; j < 16; j++)
				fprintf(m_ff,"-");
		fprintf(m_ff,"\n|\t");
		for( int i= 0; i < count_of_elements; i++)
		{
			fprintf(m_ff,"%-16s|", tagmultiMap[i].tagValue[0].c_str());
			size_t size_vec = tagmultiMap[i].tagValue.size();
			for( unsigned int j = 1; j < size_vec; j++)
			{
				fprintf(m_ff,"\n|\t");
				for( int k = 0; k < i; k++ )
					fprintf(m_ff,"%-16s|","");
				fprintf(m_ff,"%-16s|",tagmultiMap[i].tagValue[j].c_str());
			}
		}
		fprintf(m_ff,"\n");
		for(int i = 0; i < count_of_elements+1; i++)
			for(int j = 0; j < 16; j++)
				fprintf(m_ff,"-");
	
	}
	catch(char* mess)
	{
		printf("Error file open : %s\n",mess);
	}
	catch(...)
	{
	}
	if(m_fileIn.is_open())
		m_fileIn.close();
	if(m_fileOut.is_open())
		m_fileOut.close();
	if(m_ff)
		fclose(m_ff), m_ff = 0;
	return 0;
}

int ParserFile::parseString(int number, vector<string> inStrings)
{

	int count = inStrings.size();
	string inString;

	int found_count = 0;

	map <string,int> testmap;

	testmap.clear();
	for (int i = 0; i < count; i++)
	{
		inString = inStrings[i];
		int lengthStr = inString.length();
		size_t tagName_begin = inString.find('<', 0);
		size_t tagName_end = inString.find('>', tagName_begin);
		string tagName = inString.substr(tagName_begin + 1, tagName_end - 1);
		if (lengthStr == (tagName_end + 1))
			continue;

		size_t tagName_begin_begin = inString.find('<', tagName_end);
		size_t tagName_end_end = inString.find('>', tagName_begin_begin);
		string tagValue = inString.substr(tagName_end + 1, tagName_begin_begin - (tagName_end + 1));

		int prev_size = testmap.size();
		testmap.insert(pair<string,int>(tagName,i));
//		testmap.insert(tagName, tagName);
	
		if(testmap.size() == prev_size)
		{
			for(int i = 0; i < prev_size; i++)
			{
				if(tagName == tagmultiMap[i].tagName)
				{
					tagmultiMap[i].tagValue.push_back(tagValue);
					break;
				}
			}
		}
		else
		{
			Item iItem;
			iItem.tagName = tagName;
			iItem.tagValue.push_back(tagValue);
			tagmultiMap.insert(pair<int, Item>(found_count, iItem));
			found_count++;
		}
	}

	return 0;
}