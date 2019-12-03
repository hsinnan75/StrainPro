#include <unistd.h>
#include <sys/stat.h>
#include "structure.h"

int parseLine(char* line) {
	// This assumes that a digit will be found and the line ends in " Kb".
	int i = strlen(line);
	const char* p = line;
	while (*p <'0' || *p > '9') p++;
	line[i - 3] = '\0';
	i = atoi(p);
	return i;
}

void ShowFragmentPair(char* ReadSeq, FragPair_t& fp)
{
	string frag1, frag2;
	frag1.resize(fp.rLen); strncpy((char*)frag1.c_str(), ReadSeq + fp.rPos, fp.rLen);
	frag2.resize(fp.gLen); strncpy((char*)frag2.c_str(), RefSequence + fp.gPos, fp.gLen);
	printf("FragmentPair:\n%s #read[%d-%d]=%d\n%s #chr[%lld-%lld]=%d\n\n", frag1.c_str(), fp.rPos, fp.rPos + fp.rLen - 1, fp.rLen, frag2.c_str(), (long long)fp.gPos, (long long)(fp.gPos + fp.gLen - 1), fp.gLen);
}

void ShowFragPairInfo(vector<FragPair_t>& FragPairVec)
{
	for (vector<FragPair_t>::const_iterator iter = FragPairVec.begin(); iter != FragPairVec.end(); iter++)
	{
		if (iter->rLen == 0 && iter->gLen == 0) continue;
		else
		{
			//printf("\t\tFragPair#%d: R[%d-%d]=%d G[%lld-%lld]=%d %s\n", (int)(iter - FragPairVec.begin() + 1), iter->rPos, iter->rPos + iter->rLen - 1, iter->rLen, (long long)iter->gPos, (long long)(iter->gPos + iter->gLen - 1), iter->gLen, (iter->bSimple ? "Simple" : "Normal"));
			if (iter->gPos < RefSeqSize) printf("\t\tFragPair#%d: R[%d-%d]=%d G[%lld-%lld]=%d %s\n", (int)(iter - FragPairVec.begin() + 1), iter->rPos, iter->rPos + iter->rLen - 1, iter->rLen, (long long)iter->gPos, (long long)(iter->gPos + iter->gLen - 1), iter->gLen, (iter->bSimple ? "Simple" : "Normal"));
			else printf("\t\tFragPair#%d: R[%d-%d]=%d G[%lld-%lld]=%d %s (Rev)\n", (int)(iter - FragPairVec.begin() + 1), iter->rPos, iter->rPos + iter->rLen - 1, iter->rLen, (long long)(DoubleRefSeqSize - (iter->gPos + iter->gLen)), (long long)(DoubleRefSeqSize - 1 - iter->gPos), iter->gLen, (iter->bSimple ? "Simple" : "Normal"));
			//if (iter->bSimple)
			//{
			//	//char* str = new char[iter->gLen + 1]; str[iter->gLen] = '\0';
			//	//strncpy(str, RefSequence + iter->gPos, iter->gLen);
			//	////if (iter->gPos >= RefSeqSize) SelfComplementarySeq(iter->gLen, str);
			//	////printf("\t\t%s\n", str);
			//	//delete[] str;
			//}
			//else if (iter->rLen > 0 && iter->gLen > 0)
			//{
			//	printf("\t\t%s\n\t\t%s\n", iter->aln1.c_str(), iter->aln2.c_str());
			//}
		}
	}
	printf("\n"); fflush(stdout);
}

void ShowMappedRegion(vector<FragPair_t>& FragPairVec)
{
	map<int64_t, int>::iterator iter = RefSeqLocMap.lower_bound(FragPairVec.begin()->PosDiff);

	printf("gPos = %lld, ref_key=%lld, ref_idx=%d\n", (long long)FragPairVec.begin()->PosDiff, (long long)iter->first, iter->second);
	if (iter != RefSeqLocMap.end()) printf("Mapped Strain = %s\n", RefIdx->bns->anns[iter->second].name);
	else printf("Error!\n");
}

void ShowFragPairCluster(vector<AlnCan_t>& AlnCanVec)
{
	int i, num = (int)AlnCanVec.size();

	//printf("AlnCan#=%d\n", num);
	for (i = 0; i < num; i++)
	{
		if (AlnCanVec[i].score == 0) continue;
		printf("AlnCan#%d: score = %d\n", i + 1, AlnCanVec[i].score);
		ShowMappedRegion(AlnCanVec[i].FragPairVec);
		//ShowFragPairInfo(AlnCanVec[i].FragPairVec);
		printf("\n");
	}
}

string ShowTaxRank(int rank)
{
	switch (rank)
	{
	case 5: return "subspecies";
	case 10: return "species";
	case 15: return "subgenus";
	case 20: return "genus";
	case 25: return "subfamily";
	case 30: return "family";
	case 35: return "suborder";
	case 40: return "order";
	case 45: return "subclass";
	case 50: return "calss";
	case 55: return "subphylum";
	case 60: return "phylum";
	case 65: return "subkingdom";
	case 70: return "kingdom";
	default: return "root";
	}
}

int64_t FindMappedReadNum()
{
	int64_t rid, iMappedReadNum = 0;

	for (rid = 0; rid < iTotalReadNum; rid++)
	{
		if (MapInfoVec[rid].score > 0) iMappedReadNum++;
	}
	//fprintf(stderr, "\t%lld reads (%.4f%%) are mapped to reference sequences\n", iMappedReadNum, 100 * (1.0*iMappedReadNum / iTotalReadNum));
	return iMappedReadNum;
}

void ShowTaxSize()
{
	int i, taxid;
	map<int, int> TSmap;

	IndexNum = (int)IndexPrefixVec.size();
	for (IndexID = 0; IndexID < IndexNum; IndexID++)
	{
		fprintf(stderr, "Load BWT index: %s (%d / %d)", IndexPrefixVec[IndexID].c_str(), IndexID + 1, IndexNum); fflush(stderr);
		InitializeReferenceData(IndexPrefixVec[IndexID].c_str());
		for (i = 0; i < RefIdx->bns->n_seqs; i++)
		{
			taxid = RefTaxVec[i];
			TSmap[taxid] += RefIdx->bns->anns[i].len - 100;
		}
	}
	for (map<int, int>::iterator iter = TSmap.begin(); iter != TSmap.end(); iter++) printf("%d %d\n", iter->first, iter->second);
}
