#include <dirent.h>
#include "structure.h"
#include "../version.h"

time_t StartProcessTime;
vector<string> ReadLibraryVec;

float minSeqRatio;
bool bDebugMode, bFastMode;
vector<string> IndexPrefixVec;
map<int64_t, int> RefSeqLocMap;
char *IndexFileName, *OutputFilename;
string NodesDumpFilePath = "taxonomy/nodes.dmp";
string MergedDumpFilePath = "taxonomy/merged.dmp";
int iThreadNum, iRefSeqNum, MaxMismatchNum, minFrequency, minDepth;

void ShowProgramUsage(const char* program)
{
	fprintf(stderr, "\n");
	fprintf(stderr, "%s v%s\n", program, VERSION);
	fprintf(stderr, "Usage: %s -i Index_Prefix -f <ReadFile_1 ReadFile_2 ...> -o OutputPrefix\n\n", program);
	fprintf(stderr, "Options: IndexPrefix can be either an index prefix or a directory of multiple indexes\n");
	fprintf(stderr, "         -t     INT     number of threads [%d]\n", iThreadNum);
	fprintf(stderr, "         -n     INT     maximal mismatches in an alignment [%d]\n", MaxMismatchNum);
	fprintf(stderr, "\n");
}

bool CheckReadFiles()
{
	struct stat s;
	bool bRet = true;

	for (vector<string>::iterator iter = ReadLibraryVec.begin(); iter != ReadLibraryVec.end(); iter++)
	{
		if (stat(iter->c_str(), &s) == -1)
		{
			bRet = false;
			fprintf(stderr, "Cannot access file: %s\n", (char*)iter->c_str());
		}
	}
	return bRet;
}

bool CheckIndexFileName()
{
	DIR *dir;

	dir = opendir(IndexFileName);
	if (dir == NULL)
	{
		if (CheckBWAIndexFiles(IndexFileName)) IndexPrefixVec.push_back(IndexFileName);
	}
	else
	{
		int p;
		struct dirent *ds;
		string path, filename;

		path = IndexFileName; if (*path.rbegin() != '/') path.push_back('/');
		while ((ds = readdir(dir)) != NULL)
		{
			filename = ds->d_name; p = filename.find_last_of('.');
			if (filename.substr(p + 1) == "bwt")
			{
				filename = path + filename.substr(0, p);
				if (CheckBWAIndexFiles((char*)filename.c_str())) IndexPrefixVec.push_back(filename);
			}
		}
		closedir(dir);
	}
	sort(IndexPrefixVec.begin(), IndexPrefixVec.end());
	//for (vector<string>::iterator iter = IndexPrefixVec.begin(); iter != IndexPrefixVec.end(); iter++) printf("%s\n", iter->c_str());

	return IndexPrefixVec.size() > 0 ? true : false;
}

bool CheckRequirementPaths()
{
	struct stat s;

	if (stat(NodesDumpFilePath.c_str(), &s) == -1)
	{
		fprintf(stderr, "Cannot access file: %s\n", (char*)NodesDumpFilePath.c_str());
		return false;
	}
	if (stat(MergedDumpFilePath.c_str(), &s) == -1)
	{
		fprintf(stderr, "Cannot access file: %s\n", (char*)MergedDumpFilePath.c_str());
		return false;
	}
	return true;
}

int main(int argc, char* argv[])
{
	int i;
	string parameter, str;

	iThreadNum = 16;
	minSeqRatio = 0.2;
	bFastMode = false;
	MaxMismatchNum = 3;
	bDebugMode = false;
	FastQFormat = true;
	RefSequence = NULL;
	minFrequency = 50;
	OutputFilename = (char*)"output.res";

	if (argc == 1 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0)
	{
		ShowProgramUsage(argv[0]);
		exit(0);
	}
	else
	{
		for (i = 1; i < argc; i++)
		{
			parameter = argv[i];

			if (parameter == "-i" && ++i < argc) IndexFileName = argv[i];
			else if (parameter == "-o" && ++i < argc) OutputFilename = argv[i];
			else if (parameter == "-f")
			{
				while (++i < argc && argv[i][0] != '-') ReadLibraryVec.push_back(argv[i]);
				i--;
			}
			else if (parameter == "-mis" && ++i < argc) MaxMismatchNum = atoi(argv[i]);
			else if (parameter == "-freq" && ++i < argc) minFrequency = atoi(argv[i]);
			else if (parameter == "-depth" && ++i < argc) minSeqRatio = (atoi(argv[i]) / 101.0);
			else if (parameter == "-d" || parameter == "-debug") bDebugMode = true;
			//else if (parameter == "-fast") bFastMode = true;
			else if (parameter == "-t" && ++i < argc)
			{
				if ((iThreadNum = atoi(argv[i])) < 0)
				{
					fprintf(stderr, "Warning! Thread number must be positive!\n");
					iThreadNum = 16;
				}
			}
			else fprintf(stderr, "Warning! Unknow parameter: %s\n", argv[i]);
		}
		StartProcessTime = time(NULL);
		if (ReadLibraryVec.size() == 0)
		{
			fprintf(stderr, "Error! Please specify a valid read input!\n");
			ShowProgramUsage(argv[0]);
			exit(1);
		}
		else if (CheckReadFiles() == false) exit(1);

		if (IndexFileName == NULL || CheckIndexFileName() == false)
		{
			fprintf(stderr, "Error! Please specify a valid reference index or direcotry with indices!\n");
			ShowProgramUsage(argv[0]);
			exit(1);
		}
		if (CheckRequirementPaths() == false) exit(1);

		//ShowTaxSize();
		GetTaxInfomation(); MetaMapping(); MetaTyping();
		fprintf(stderr, "\nDone! It took %lld seconds.\n", (long long)(time(NULL) - StartProcessTime));
	}
	return 0;
}
