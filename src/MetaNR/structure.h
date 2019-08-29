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

#define ReadLength 101

using namespace std;

typedef unsigned char ubyte_t;
typedef unsigned long long bwtint_t;

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
	bwt_t    *bwt; // FM-index
	bntseq_t *bns; // information on the reference sequences
	uint8_t  *pac; // the actual 2-bit encoded reference sequences with 'N' converted to a random base
} bwaidx_t;

typedef struct
{
	int parent_taxid;
	uint8_t rank; //0:sub-species,1:species,2:genus,3:family,4:order,5:class,6:phylum,7:kingdom
} TaxItem_t;

// Global variables
extern bwt_t *Refbwt;
extern int *RefTaxArr;
extern bwaidx_t *RefIdx;
extern time_t StartProcessTime;
extern map<int, TaxItem_t> TaxMap;
extern map<int64_t, int> RefPosTaxMap;
extern unsigned char nst_nt4_table[256];
extern int64_t RefSeqSize, DoubleRefSeqSize;
extern int iThreadNum;
extern char *RefSequence, *IndexPrefix, *OutputFASTA;
extern string NodesDumpFilePath, MergedDumpFilePath, ClusterFolder;

// bwt_index.cpp
extern void RestoreReferenceInfo();
extern void bwa_idx_destroy(bwaidx_t *idx);
extern bwaidx_t *bwa_idx_load(const char *hint);

// bwt_search.cpp
extern void BWT_Search(int64_t myPos, char *seq);

// SeqTaxing.cpp
extern void SeqTaxing();

// tools.cpp
extern int CheckMemoryUsage();
extern string ShowTaxRank(int rank);
extern void ShowReadSequence(int64_t gPos);

// taxonomy.cpp
extern void GetTaxInfomation();
extern int LCA(int taxid1, int taxid2);
extern int Find_LCA(vector<int>& TaxVec);
