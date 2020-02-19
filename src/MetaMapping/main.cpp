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
int iThreadNum, iRefSeqNum, MaxMismatchNum, minFrequency, minDepth;

string StrainProDir, TaxonomyDir;

void ShowProgramUsage(const char* program)
{
	fprintf(stderr, "\n");
	fprintf(stderr, "%s v%s\n", program, VERSION);
	fprintf(stderr, "Usage: %s -i Index_Prefix -f <ReadFile_1 ReadFile_2 ...> -o OutputPrefix\n\n", program);
	fprintf(stderr, "Options: IndexPrefix can be either an index prefix or a directory of multiple indexes\n");
	fprintf(stderr, "         -t           INT     number of threads [%d]\n", iThreadNum);
	fprintf(stderr, "         -n           INT     maximal mismatches in an alignment [%d]\n", MaxMismatchNum);
	fprintf(stderr, "         -depth       INT     minimal read depth of sequenced data [20]\n");
	fprintf(stderr, "         -freq        INT     minimal read count for identified genomes [%d]\n", minFrequency);
	fprintf(stderr, "         -dump        STRING  dump file path\n");
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

int is_regular_file(const char *path)
{
	struct stat path_stat;
	stat(path, &path_stat);
	return S_ISREG(path_stat.st_mode);
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

	if (stat(StrainProDir.c_str(), &s) == -1)
	{
		fprintf(stderr, "Cannot access %s\n", (char*)StrainProDir.c_str());
		return false;
	}
	if (stat(TaxonomyDir.c_str(), &s) == -1)
	{
		fprintf(stderr, "Cannot access %s\n", (char*)TaxonomyDir.c_str());
		return false;
	}
	return true;
}

void ParseEnvSetting(char * envp[])
{
	string str;

	for (int i = 0; envp[i] != NULL; i++)
	{
		str = envp[i];
		if (str.find("StrainPro_DIR") == 0)
		{
			StrainProDir = str.substr(14);
			if (StrainProDir[StrainProDir.length() - 1] == '/') StrainProDir.resize(StrainProDir.length() - 1);
			TaxonomyDir = StrainProDir.substr(0, StrainProDir.find_last_of('/')) + "/taxonomy";
		}
	}
}

int main(int argc, char* argv[], char* envp[])
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
	OutputFilename = (char*)"output.tsv";

	ParseEnvSetting(envp);

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
			else if (parameter == "-t" && ++i < argc)
			{
				if ((iThreadNum = atoi(argv[i])) < 0)
				{
					fprintf(stderr, "Warning! Thread number must be positive!\n");
					iThreadNum = 16;
				}
			}
			else if (parameter == "-dump" && i + 1 < argc) TaxonomyDir = argv[++i];
			else fprintf(stderr, "Warning! Unknow parameter: %s\n", argv[i]);
		}
		if (StrainProDir == "")
		{
			fprintf(stderr, "Cannot find StrainPro_bin directory, please add the environment variable StrainPro_DIR to the StrainPro's bin directory!\n");
			exit(0);
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
