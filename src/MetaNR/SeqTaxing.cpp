#include "structure.h"

#define JobBlockSize 10000

int* RefTaxArr;
int64_t cur_Ref_gPos;
map<int64_t, int> RefPosTaxMap;
map<int, int64_t> TaxRankSizeMap;
static pthread_mutex_t MutexLock;

void InitializeRefPosTaxMap()
{
	int64_t gPos;
	string str, tax;
	int i, p, taxid;

	for (gPos = 0, i = 0; i < RefIdx->bns->n_seqs; i++)
	{
		str = RefIdx->bns->anns[i].name; p = str.find_first_of('|', 6); tax = str.substr(6, p - 6); taxid = atoi(tax.c_str());
		gPos += RefIdx->bns->anns[i].len; RefPosTaxMap.insert(make_pair(gPos - 1, taxid));
		//printf("[%lld - %lld]: %s (%d)\n", gPos - RefIdx->bns->anns[i].len, gPos - 1, RefIdx->bns->anns[i].name, taxid);
	}
}

void InitializeRefTaxArr()
{
	int i;
	int64_t iTotalBase, gPos;

	iTotalBase = 0;
	for (i = 0; i < RefIdx->bns->n_seqs; i++)
	{
		gPos = iTotalBase; iTotalBase += RefIdx->bns->anns[i].len; 
		memset(RefTaxArr + gPos, -1, sizeof(int)*RefIdx->bns->anns[i].len);
		//for (gPos = iTotalBase - RefIdx->bns->anns[i].len; gPos < read_end_pos; gPos++) RefTaxArr[gPos] = -1;
		// Mask read sequences that span two distinct chromosomes
		for (gPos = iTotalBase - ReadLength + 1; gPos < iTotalBase; gPos++) RefTaxArr[gPos] = 0;
	}
}

void *PseudoReadMapping(void *arg)
{
	int64_t gPos, start_pos, stop_pos;

	while (true)
	{
		pthread_mutex_lock(&MutexLock);
		start_pos = cur_Ref_gPos; stop_pos = (cur_Ref_gPos += JobBlockSize);
		if (stop_pos > RefSeqSize) stop_pos = cur_Ref_gPos = RefSeqSize;
		fprintf(stderr, "\33[2K\r\tScan %lld nucleobases (%d%%)...", (long long)start_pos, (int)(100 * (1.0*start_pos / RefSeqSize)));
		pthread_mutex_unlock(&MutexLock);

		if (start_pos == RefSeqSize) break;
		for (gPos = start_pos; gPos < stop_pos; gPos++)
		{
			if (RefTaxArr[gPos] == -1) BWT_Search(gPos, RefSequence + gPos);
		}
	}
	return (void*)(1);
}

void OutputAnnotatedSequences()
{
	int64_t gPos;
	FILE *out_file;
	int n, frag_len, taxid;

	out_file = fopen(OutputFASTA, "w");
	for (taxid = 0, n = 0, gPos = 0; gPos < RefSeqSize; gPos++)
	{
		if (RefTaxArr[gPos] == taxid) n++;
		else
		{
			if (taxid > 0)
			{
				frag_len = n + ReadLength - 1;
				TaxRankSizeMap[TaxMap[taxid].rank] += frag_len;
				//printf("%lld - %lld: tax=%d (rank=%d), len=%d \n", gPos - n, gPos - n + frag_len - 1, taxid, TaxMap[taxid].rank, frag_len);
				fprintf(out_file, ">taxid|%d|size=%d\n%.*s\n", taxid, frag_len, frag_len, (RefSequence + gPos - n));
			}
			taxid = RefTaxArr[gPos]; n = 1;
		}
	}
	//for (map<int, int64_t>::iterator iter = TaxRankSizeMap.begin(); iter != TaxRankSizeMap.end(); iter++) fprintf(stderr, "TaxRank=%s(%d), size=%lld (%.5f)\n", ShowTaxRank(iter->first).c_str(), iter->first, iter->second, 1.*iter->second / RefSeqSize);
	fclose(out_file);
}

void SeqTaxing()
{
	int i;

	fprintf(stderr, "Step4. Annotate all sequence fragments in the reference sequences.\n");
	
	RefTaxArr = new int[RefSeqSize];
	//for (cur_Ref_gPos = 0; cur_Ref_gPos < RefSeqSize; cur_Ref_gPos++) RefTaxArr[cur_Ref_gPos] = -1;
	InitializeRefTaxArr(); InitializeRefPosTaxMap();

	////iThreadNum = 1;
	cur_Ref_gPos = 0; pthread_t *ThreadArr = new pthread_t[iThreadNum];
	for (i = 0; i < iThreadNum; i++) pthread_create(&ThreadArr[i], NULL, PseudoReadMapping, NULL);
	for (i = 0; i < iThreadNum; i++) pthread_join(ThreadArr[i], NULL);
	delete[] ThreadArr;
	fprintf(stderr, "\nStep5. Output annotated sequence fragments.\n"); OutputAnnotatedSequences();

	delete[] RefTaxArr;
}
