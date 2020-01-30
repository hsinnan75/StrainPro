#include "../version.h"
#include "structure.h"

time_t StartProcessTime;

int iThreadNum;
string OutputFolder;
vector<SeqInfo_t> SeqVec;
string ReferenceFilename;
string MetaNR_Path = "bin/StrainPro-rep";
int64_t TotalSeqSize = 0, NRS_Size = 0;
string NodesDumpFilePath = "taxonomy/nodes.dmp";
string MergedDumpFilePath = "taxonomy/merged.dmp";

void ShowProgramUsage(const char* program)
{
	fprintf(stderr, "\n");
	fprintf(stderr, "%s v%s\n", program, VERSION);
	fprintf(stderr, "Usage: %s -r DB_SeqFile[fa] -o OutputFolder\n\n", program);
	fprintf(stderr, "Options: -t     INT     number of threads [%d]\n", iThreadNum);
	fprintf(stderr, "         -dump  STRING  dump file path\n");
	fprintf(stderr, "\n");
}

void LoadDumpFilePath(const char* filename)
{
	fstream file, f;
	stringstream ss;
	string str, s1, s2;

	file.open(filename, ios_base::in);
	if (!file.is_open())
	{
		fprintf(stderr, "cannot open file %s\n", filename);
		exit(1);
	}
	while (!file.eof())
	{
		getline(file, str); if (str == "") continue;
		ss.clear(); ss.str(str); ss >> s1 >> s2;
		if (s1 == "NodesDumpFilePath")
		{
			f.close(); f.open(s2.c_str());
			if (!f.is_open())
			{
				fprintf(stderr, "Error! File (%s) is not accessible\n", s2.c_str());
				exit(1);
			}
			else NodesDumpFilePath = s2;
		}
		else if (s1 == "MergedDumpFilePath")
		{
			f.close(); f.open(s2.c_str());
			if (!f.is_open())
			{
				fprintf(stderr, "Error! File (%s) is not accessible\n", s2.c_str());
				exit(1);
			}
			else MergedDumpFilePath = s2;
		}
	}
	file.close();

	fprintf(stderr, "Use the dump files: %s and %s\n", NodesDumpFilePath.c_str(), MergedDumpFilePath.c_str());
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
	if (stat(ReferenceFilename.c_str(), &s) == -1)
	{
		fprintf(stderr, "Cannot access file: %s\n", (char*)ReferenceFilename.c_str());
		return false;
	}
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
	if (stat(MetaNR_Path.c_str(), &s) == -1)
	{
		fprintf(stderr, "Cannot access file: %s\n", (char*)MetaNR_Path.c_str());
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

int main(int argc, char* argv[])
{
	int i;
	string parameter, str;

	iThreadNum = 16;

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
			else if (parameter == "-r" && i+1<argc) ReferenceFilename = argv[++i];
			else if (parameter == "-t" && i + 1 < argc)
			{
				if ((iThreadNum = atoi(argv[++i])) < 0) iThreadNum = 16;
			}
			else if (parameter == "-dump" && ++i < argc) LoadDumpFilePath(argv[i]);
			else fprintf(stderr, "Warning! Unknow parameter: %s\n", argv[i]);
		}
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
