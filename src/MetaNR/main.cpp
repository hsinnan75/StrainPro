#include "structure.h"
#include <sys/stat.h>
#include "../version.h"

time_t StartProcessTime;

bwt_t *Refbwt;
int iThreadNum;
bwaidx_t *RefIdx;
string StrainProDir, TaxonomyDir;
int64_t RefSeqSize, DoubleRefSeqSize;
char *RefSequence, *IndexPrefix, *OutputFASTA;

void ShowProgramUsage(const char* program)
{
	fprintf(stderr, "\n");
	fprintf(stderr, "%s v%s\n", program, VERSION);
	fprintf(stderr, "Usage: %s -i IndexPrefix\n\n", program);
	fprintf(stderr, "Options: -t     INT     number of threads [%d]\n", iThreadNum);
	fprintf(stderr, "         -o     STR     output fasta [%s]\n", OutputFASTA);
	fprintf(stderr, "         -dump  STRING  dump file path\n");
	fprintf(stderr, "\n");
}

int is_regular_file(const char *path)
{
	struct stat path_stat;
	stat(path, &path_stat);
	return S_ISREG(path_stat.st_mode);
}

bool CheckBWAIndexFiles()
{
	char fn[1024];
	bool bRet = true;
	struct stat buffer;

	sprintf(fn, "%s.ann", IndexPrefix); if (stat(fn, &buffer) != 0) bRet = false, fprintf(stderr, "Cannot access file: %s\n", fn);
	sprintf(fn, "%s.amb", IndexPrefix); if (stat(fn, &buffer) != 0) bRet = false, fprintf(stderr, "Cannot access file: %s\n", fn);
	sprintf(fn, "%s.pac", IndexPrefix); if (stat(fn, &buffer) != 0) bRet = false, fprintf(stderr, "Cannot access file: %s\n", fn);
	sprintf(fn, "%s.bwt", IndexPrefix); if (stat(fn, &buffer) != 0) bRet = false, fprintf(stderr, "Cannot access file: %s\n", fn);
	sprintf(fn, "%s.sa",  IndexPrefix); if (stat(fn, &buffer) != 0) bRet = false, fprintf(stderr, "Cannot access file: %s\n", fn);

	return bRet;
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

int main(int argc, char* argv[], char* envp[])
{
	int i;
	string parameter, str;

	iThreadNum = 16;
	OutputFASTA = (char*)"output.fasta";

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

			if (parameter == "-i" && i + 1 < argc) IndexPrefix = argv[++i];
			else if (parameter == "-o" && i + 1 < argc) OutputFASTA = argv[++i];
			else if (parameter == "-t" && i + 1 < argc)
			{
				if ((iThreadNum = atoi(argv[++i])) < 0) iThreadNum = 16;
			}
			else if (parameter == "-dump" && i + 1 < argc) TaxonomyDir = argv[++i];
			else fprintf(stderr, "Warning! Unknow parameter: %s\n", argv[i]);
		}
	}
	if (StrainProDir == "")
	{
		fprintf(stderr, "Cannot find StrainPro_bin directory, please add the environment variable StrainPro_DIR to the StrainPro's bin directory!\n");
		exit(0);
	}
	if (CheckRequirementPaths() == false) exit(1);

	StartProcessTime = time(NULL);

	if (IndexPrefix == NULL || CheckBWAIndexFiles() == false)
	{
		fprintf(stderr, "\n\nError! Index files are corrupt!\n");
		exit(1);
	}
	GetTaxInfomation();
	RefIdx = bwa_idx_load(IndexPrefix);
	if (RefIdx == 0)
	{
		fprintf(stderr, "\n\nError! Index files are corrupt!\n");
		exit(1);
	}
	else
	{
		Refbwt = RefIdx->bwt;
		RestoreReferenceInfo();
	}
	SeqTaxing();
	fprintf(stderr, "Done! It took %lld seconds.\n\n", (long long)(time(NULL) - StartProcessTime));
	
	bwa_idx_destroy(RefIdx);

	return 0;
}
