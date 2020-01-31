#include "structure.h"

vector<string> RepDB_vec;
//static pthread_mutex_t Lock;

void Make_nrDB()
{
	char cmd[1024];
	string IdxPrefix, output_fn;

	for (vector<pair<string, int64_t> >::iterator iter = ClusterSeqPathVec.begin(); iter != ClusterSeqPathVec.end(); iter++)
	{
		IdxPrefix = iter->first.substr(0, iter->first.find_last_of('.'));
		output_fn = IdxPrefix + ".rep"; RepDB_vec.push_back(output_fn);

		sprintf(cmd, "%s/StrainPro-rep -i %s -t %d -o %s", StrainProDir.c_str(), IdxPrefix.c_str(), iThreadNum, output_fn.c_str());
		fprintf(stderr, "cmd=%s\n", cmd); system(cmd);
	}
}

int64_t Load_All_NRS()
{
	fstream file;
	int64_t total_n;
	int taxid, p1, p2;
	SeqInfo_t SeqInfo;
	string str, header, tax;

	fprintf(stderr, "Get all representative sequence segments...\n"); taxid = 0; total_n = 0;
	for (vector<string>::iterator iter = RepDB_vec.begin(); iter != RepDB_vec.end(); iter++)
	{
		file.clear(); file.open(iter->c_str(), ios_base::in);
		while (!file.eof())
		{
			getline(file, str); if (str == "") continue;
			getline(file, SeqInfo.seq); total_n += SeqInfo.seq.length();
			p1 = str.find("taxid|") + 6; p2 = str.find_first_of('|', p1);
			tax = str.substr(p1, p2 - p1); SeqInfo.taxid = atoi(tax.c_str());
			SeqInfo.header = str.substr(1); SeqVec.push_back(SeqInfo);
		}
		file.close();
	}
	return total_n;
}
