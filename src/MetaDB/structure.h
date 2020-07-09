#include <iostream>
#include <string>
#include <stdio.h>
#include <cstring>
#include <fstream>
#include <sstream>
#include <map>
#include <cstdio>
#include <cmath>
#include <vector>
#include <iterator>
#include <algorithm>
#include <ctime>
#include <pthread.h>
#include <sys/time.h>
#include <inttypes.h>
#include <dirent.h>
#include <sys/stat.h>

using namespace std;

typedef struct
{
	int parent_taxid;
	uint8_t rank; //0:sub-species,1:species,2:genus,3:family,4:order,5:class,6:phylum,7:kingdom
} TaxItem_t;

typedef struct
{
	int taxid;
	string header;
	string seq;
} SeqInfo_t;

// Global variables
extern int iThreadNum;
extern const char* TaxLevel;
extern time_t StartProcessTime;
extern vector<SeqInfo_t> SeqVec;
extern map<int, TaxItem_t> TaxMap;
extern map<string, int> TaxRankMap;
extern vector<string> ClusterSeqPath;
extern string StrainProDir, TaxonomyDir, OutputFolder;
extern vector<pair<string, int64_t> > ClusterSeqPathVec;

// nr_db.cpp
extern void Make_nrDB();
extern void Make_nrDB_Index();
extern int64_t Load_All_NRS();

// seq_clustering.cpp
extern void Make_DB_Index();
extern void Seq_Clustering(const char* clr_suffix);

// taxonomy.cpp
extern void GetTaxInfomation();

// tools.cpp
extern void NR_Analysis();
extern void RemoveSeqFiles();
extern void Remove_BWT_Files();
extern int CheckMemoryUsage();
