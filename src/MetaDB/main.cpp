#include "../version.h"
#include "structure.h"

time_t StartProcessTime;

int iThreadNum;
vector<SeqInfo_t> SeqVec;
int64_t TotalSeqSize = 0, NRS_Size = 0;
string StrainProDir, TaxonomyDir, ReferenceFilename, OutputFolder;

void ShowProgramUsage(const char* program)
{
	fprintf(stderr, "\n");
	fprintf(stderr, "%s v%s\n", program, VERSION);
	fprintf(stderr, "Usage: %s -r DB_SeqFile[fa] -o OutputFolder\n\n", program);
	fprintf(stderr, "Options: -t     INT     number of threads [%d]\n", iThreadNum);
	fprintf(stderr, "         -dump  STRING  dump file path\n");
	fprintf(stderr, "\n");
}

int is_regular_file(const char *path)
{
	struct stat path_stat;
	stat(path, &path_stat);
	return S_ISREG(path_stat.st_mode);
}

bool CheckRequirementPaths()
{
	struct stat s;

	if (OutputFolder == "")
	{
		fprintf(stderr, "Please specify an output folder\n");
		return false;
	}
	if (stat(OutputFolder.c_str(), &s) != -1)
	{
		fprintf(stderr, "OutputFolder %s exists\n", (char*)OutputFolder.c_str());
		return false;
	}
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
	if (stat(ReferenceFilename.c_str(), &s) == -1)
	{
		fprintf(stderr, "Cannot access file: %s\n", (char*)ReferenceFilename.c_str());
		return false;
	}
	return true;
}

int64_t GetDBseq(const char* filename)
{
	int p1, p2;
	fstream file;
	int64_t total_n;
	SeqInfo_t SeqInfo;
	string str, seq, header, tax;

	fprintf(stderr, "Get all sequences...\n");
	file.open(filename, ios_base::in); total_n = 0;
	while (!file.eof())
	{
		getline(file, str); if (str == "") continue;
		if (str[0] == '>')
		{
			if (SeqInfo.seq != "")
			{
				total_n += SeqInfo.seq.length();
				SeqVec.push_back(SeqInfo);
			}
			//reset seq_item
			if ((p1 = str.find("taxid|")) == -1)
			{
				fprintf(stderr, "Warning! [%s] does not contain the taxid label\n", str.c_str());
				continue;
			}
			p1 += 6; p2 = str.find_first_of('|', p1); tax = str.substr(p1, p2 - p1); 
			SeqInfo.taxid = atoi(tax.c_str());
			SeqInfo.header = str.substr(1); SeqInfo.seq.clear();
		}
		else SeqInfo.seq += str;
	}
	if (SeqInfo.seq != "")
	{
		total_n += SeqInfo.seq.length();
		SeqVec.push_back(SeqInfo);
	}
	return total_n;
}

void InitializeClusterPath()
{
	DIR *dir;

	dir = opendir(OutputFolder.c_str());
	if (dir != NULL)
	{
		int p;
		struct dirent *ds;
		string path, filename;

		while ((ds = readdir(dir)) != NULL)
		{
			filename = ds->d_name; p = filename.find_last_of('.');
			if (filename.substr(p + 1) == "fna")
			{
				filename = OutputFolder + "/" + filename;
				ClusterSeqPathVec.push_back(make_pair(filename, 0));
			}
		}
		closedir(dir);
	}
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

int main(int argc, char* argv[], char * envp[])
{
	int i;
	string parameter, str;

	iThreadNum = 16;

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

			if (parameter == "-o" && i + 1 < argc) OutputFolder = argv[++i];
			else if (parameter == "-r" && i + 1 < argc) ReferenceFilename = argv[++i];
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

	StartProcessTime = time(NULL);

	if (CheckRequirementPaths() == false) exit(1);
	else system(("mkdir -p " + OutputFolder).c_str());

	//clustering all sequences by taxonomy
	GetTaxInfomation();
	TotalSeqSize = GetDBseq(ReferenceFilename.c_str());
	Seq_Clustering("clr1"); SeqVec.clear(); Make_DB_Index();

	//making non-redundant sequences
	Make_nrDB(); Remove_BWT_Files();
	NRS_Size = Load_All_NRS();
	Seq_Clustering("clr2"); SeqVec.clear(); Make_DB_Index();

	//removing temporary files
	RemoveSeqFiles();

	fprintf(stderr, "Original db_size = %dMb, representative sequence db_size = %dMb\n", (int)(TotalSeqSize/1000000), (int)(NRS_Size/1000000));
	fprintf(stderr, "Done! It took %lld seconds.\n\n", (long long)(time(NULL) - StartProcessTime));

	return 0;
}
