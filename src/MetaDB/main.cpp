#include "structure.h"

time_t StartProcessTime;
const char* VersionStr = "0.9.0";

int iThreadNum;
string OutputFolder;
vector<SeqInfo_t> SeqVec;
string ReferenceFilename;
string MetaNR_Path = "./StrainPro-rep";
int64_t TotalSeqSize = 0, NRS_Size = 0;
string NodesDumpFilePath = "taxonomy/nodes.dmp";
string MergedDumpFilePath = "taxonomy/merged.dmp";

void ShowProgramUsage(const char* program)
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage: %s -r DB_SeqFile[fa] -o OutputFolder\n\n", program);
	fprintf(stderr, "Options: -t     INT     number of threads [%d]\n", iThreadNum);
	fprintf(stderr, "         -v             development version\n");
	fprintf(stderr, "\n");
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
	int p;
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
			p = str.find_first_of('|', 7); tax = str.substr(7, p - 7); 
			//taxid = atoi(tax.c_str()); while (TaxMap[taxid].rank < 10 && TaxMap[TaxMap[taxid].parent_taxid].rank < 10) taxid = TaxMap[taxid].parent_taxid;
			//if (taxid != atoi(tax.c_str())) printf("%s --> %d\n", tax.c_str(), taxid);
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

	if (argc == 1 || strcmp(argv[1], "-h") == 0)
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
			else if (parameter == "-v")
			{
				fprintf(stderr, "%s v%s\n", argv[0], VersionStr);
				exit(0);
			}
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

	fprintf(stderr, "Original db_size = %lld, representative sequence db_size = %lld\n", (long long)TotalSeqSize, (long long)NRS_Size);
	fprintf(stderr, "Done! It took %lld seconds. (MemUsage: %d MB)\n\n", (long long)(time(NULL) - StartProcessTime), CheckMemoryUsage());

	return 0;
}
