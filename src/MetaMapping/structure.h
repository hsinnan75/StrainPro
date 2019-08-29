#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include <iterator>
#include <algorithm>
#include <ctime>
#include <zlib.h>
#include <inttypes.h>
#include <pthread.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <unistd.h>

using namespace std;

#define ReadLength 101
#define MaxPosDiffGap 3
#define MinSeedLength 16
#define ReadChunkSize 10000

typedef uint64_t bwtint_t;
typedef unsigned char ubyte_t;

typedef struct {
	bwtint_t primary; // S^{-1}(0), or the primary index of BWT
	bwtint_t L2[5]; // C(), cumulative count
	bwtint_t seq_len; // sequence length
	bwtint_t bwt_size; // size of bwt, about seq_len/4
	uint32_t *bwt; // BWT
	uint32_t cnt_table[256];
	int sa_intv;
	bwtint_t n_sa;
	bwtint_t *sa;
} bwt_t;

typedef struct {
	int64_t offset;
	int32_t len;
	int32_t n_ambs;
	uint32_t gi;
	char *name, *anno;
} bntann1_t;

typedef struct {
	int64_t offset;
	int32_t len;
	char amb;
} bntamb1_t;

typedef struct {
	int64_t l_pac;
	int32_t n_seqs;
	uint32_t seed;
	bntann1_t *anns; // n_seqs elements
	int32_t n_holes;
	bntamb1_t *ambs; // n_holes elements
	FILE *fp_pac;
} bntseq_t;

typedef struct {
	bwt_t    *bwt; // FM-index
	bntseq_t *bns; // information on the reference sequences
	uint8_t  *pac; // the actual 2-bit encoded reference sequences with 'N' converted to a random base
} bwaidx_t;

typedef struct
{
	bwtint_t x[3];
} bwtintv_t;

typedef struct
{
	int len;
	int freq;
	bwtint_t* LocArr;
} bwtSearchResult_t;

typedef struct
{
	int parent_taxid;
	uint8_t rank; //0:sub-species,1:species,2:genus,3:family,4:order,5:class,6:phylum,7:kingdom
} TaxItem_t;

typedef struct
{
	bool bSimple;
	int rPos; // read position
	int64_t gPos; // genome position
	int rLen; // read block size
	int gLen; // genome block size
	int64_t PosDiff; // gPos-rPos
	string aln1; // read fragment alignment
	string aln2; // genomic fragment alignment
} FragPair_t;

typedef struct
{
	int score;
	vector<FragPair_t> FragPairVec;
} AlnCan_t;

typedef struct
{
	int score;
	int taxid;
	int ref_idx;
} AlnSummary_t;

typedef struct
{
	int rlen;
	char* seq;
	//char* header;
	AlnSummary_t AlnSummary;
	vector<AlnCan_t> AlnCanVec;
} ReadItem_t;

typedef struct
{
	int taxid;
	int ref_idx;
	unsigned short score;
	unsigned short idx_id;
} MapInfo_t;

// Global variables
extern bwt_t *Refbwt;
extern bwaidx_t *RefIdx; 
extern const char* VersionStr;
extern time_t StartProcessTime;
extern map<int, TaxItem_t> TaxMap;
extern map<string, int> TaxRankMap;
extern unsigned char nst_nt4_table[256];

extern float minSeqRatio;
extern vector<int> RefTaxVec;
extern vector<MapInfo_t> MapInfoVec;
extern vector<string> IndexPrefixVec;
extern vector<string> ReadLibraryVec;
extern map<int64_t, int> RefSeqLocMap;
extern string NodesDumpFilePath, MergedDumpFilePath;
extern char *RefSequence, *IndexFileName, *OutputFilename;
extern int64_t RefSeqSize, DoubleRefSeqSize, iTotalReadNum;
extern bool bDebugMode, bFastMode, FastQFormat, gzCompressed;
extern int iThreadNum, iRefSeqNum, MaxMismatchNum, IndexID, IndexNum, minFrequency, minDepth;

// bwt_index.cpp
extern void bns_destroy(bntseq_t *bns);
extern void RestoreReferenceSequences();
extern void bwa_idx_destroy(bwaidx_t *idx);
extern bwaidx_t *bwa_idx_load(const char* IndexPrefix);
extern bntseq_t *bns_restore(const char* IndexPrefix);
extern void InitializeReferenceData(const char* IndexPrefix);

// bwt_search.cpp
extern bwtSearchResult_t BWT_Search(char* seq, int start, int stop);

// Mapping.cpp
extern void MetaMapping();

// Typing.cpp
extern void MetaTyping();

// ReadAlignment.cpp
extern bool ProduceReadAlignment(ReadItem_t& read);

// GetData.cpp
extern bool CheckBWAIndexFiles(char* IndexPrefix);
extern bool CheckReadFormat(const char* filename);
extern int GetNextChunk(FILE *file, ReadItem_t* ReadArr);
extern int ConvertReadLibrary(char* ReadLibraryFileName);
extern int GetMappingChunk(FILE *file, ReadItem_t* ReadArr);
extern int gzGetNextChunk(gzFile file, ReadItem_t* ReadArr);

// tools.cpp
extern void ShowTaxSize();
extern int CheckMemoryUsage();
extern int64_t FindMappedReadNum();
extern string ShowTaxRank(int rank);
//extern string GenerateRandomString(int len);
extern void ShowFragPair(FragPair_t& FragPair);
extern void ShowFragPairInfo(vector<FragPair_t>& FragPairVec);
extern void GetComplementarySeq(int len, char* seq, char* rseq);

// nw_alignment.cpp
extern void nw_alignment(int m, string& s1, int n, string& s2);

// taxonomy.cpp
extern void GetTaxInfomation();
extern int Pairwise_LCA(int taxid1, int taxid2);
