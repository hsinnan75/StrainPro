#include "structure.h"

//static pthread_mutex_t Lock;

float CalFragAlnSeqIdy(string& aln1, string& aln2)
{
	int i, len = (int)aln1.length(), n = 0, mis = 0;

	for (i = 0; i < len; i++)
	{
		if (aln1[i] != '-' && aln2[i] != '-')
		{
			n++;
			if (aln1[i] != aln2[i]) mis++;
		}
	}
	if (n != len && mis > 0) return 0;
	else if(n > 0) return 1.0*(n - mis) / n;
	else return 1;
}

bool CompByReadPos(const FragPair_t& p1, const FragPair_t& p2)
{
	if (p1.rPos == p2.rPos) return p1.gPos < p2.gPos;
	else return p1.rPos < p2.rPos;
}

void RemoveNullFragPairs(vector<FragPair_t>& FragPairVec)
{
	for (vector<FragPair_t>::iterator iter = FragPairVec.begin(); iter != FragPairVec.end();)
	{
		if (iter->rLen == 0) iter = FragPairVec.erase(iter);
		else iter++;
	}
}

bool RemoveOverlaps(vector<FragPair_t>& FragPairVec)
{
	bool bOverlap = false;
	int i, j, overlap_size, num = (int)FragPairVec.size();

	for (i = 0, j = 1; j < num; i++, j++)
	{
		if (FragPairVec[i].rPos == FragPairVec[j].rPos)
		{
			bOverlap = true;
			FragPairVec[i].rLen = FragPairVec[i].gLen = 0;
		}
		else if (FragPairVec[i].gPos >= FragPairVec[j].gPos || (FragPairVec[i].gPos + FragPairVec[i].gLen) > FragPairVec[j].gPos)
		{
			// shrink FragPairVec[i]
			bOverlap = true;
			overlap_size = FragPairVec[i].gPos + FragPairVec[i].gLen - FragPairVec[j].gPos;
			if ((FragPairVec[i].rLen -= overlap_size) < 0) FragPairVec[i].rLen = 0;
			if ((FragPairVec[i].gLen -= overlap_size) < 0) FragPairVec[i].gLen = 0;
		}
	}
	//if (bOverlap && bDebugMode)
	//{
	//	printf("After RemoveOverlaps\n");
	//	ShowSimplePairInfo(FragPairVec);
	//}
	return bOverlap;
}

void IdentifyNormalPairs(int rlen, vector<FragPair_t>& FragPairVec)
{
	FragPair_t FragPair;
	int i, j, rGaps, gGaps, num = (int)FragPairVec.size();

	FragPair.bSimple = false;

	for (i = 0, j = 1; j < num; i++, j++)
	{
		if ((rGaps = FragPairVec[j].rPos - (FragPairVec[i].rPos + FragPairVec[i].rLen)) < 0) rGaps = 0;
		if ((gGaps = FragPairVec[j].gPos - (FragPairVec[i].gPos + FragPairVec[i].gLen)) < 0) gGaps = 0;

		if (rGaps > 0 || gGaps > 0)
		{
			FragPair.rPos = FragPairVec[i].rPos + FragPairVec[i].rLen;
			FragPair.gPos = FragPairVec[i].gPos + FragPairVec[i].gLen;
			FragPair.PosDiff = FragPair.gPos - FragPair.rPos;
			FragPair.rLen = rGaps; FragPair.gLen = gGaps;
			FragPairVec.push_back(FragPair);
			//if (bDebugMode) printf("insert a normal pair: r[%d-%d] g[%lld-%lld] and r[%d-%d] g[%lld-%lld]: r[%d-%d] g[%lld-%lld]\n", FragPairVec[i].rPos, FragPairVec[i].rPos + FragPairVec[i].rLen - 1, FragPairVec[i].gPos, FragPairVec[i].gPos + FragPairVec[i].gLen - 1, FragPairVec[j].rPos, FragPairVec[j].rPos + FragPairVec[j].rLen - 1, FragPairVec[j].gPos, FragPairVec[j].gPos + FragPairVec[j].gLen - 1, FragPair.rPos, FragPair.rPos + FragPair.rLen - 1, FragPair.gPos, FragPair.gPos + FragPair.gLen - 1);
		}
	}
	if ((int)FragPairVec.size() > num) inplace_merge(FragPairVec.begin(), FragPairVec.begin() + num, FragPairVec.end(), CompByReadPos);

	// Check missing blocks at both ends
	if (FragPairVec[0].rPos > 0)
	{
		FragPair.rPos = 0;
		FragPair.gPos = FragPair.PosDiff = FragPairVec[0].PosDiff;
		FragPair.rLen = FragPair.gLen = FragPairVec[0].rPos;
		FragPairVec.insert(FragPairVec.begin(), FragPair);
	}
	num = (int)FragPairVec.size();
	if ((FragPairVec[num - 1].rPos + FragPairVec[num - 1].rLen) < rlen)
	{
		FragPair.rPos = FragPairVec[num - 1].rPos + FragPairVec[num - 1].rLen;
		FragPair.gPos = FragPairVec[num - 1].gPos + FragPairVec[num - 1].gLen;
		FragPair.PosDiff = FragPairVec[num - 1].PosDiff;
		FragPair.rLen = FragPair.gLen = rlen - FragPair.rPos;
		FragPairVec.push_back(FragPair);
	}
}

void RemoveEmptyFragPairs(vector<FragPair_t>& FragPairVec)
{
	for (vector<FragPair_t>::iterator iter = FragPairVec.begin(); iter != FragPairVec.end();)
	{
		if (iter->rLen == 0 && iter->gLen == 0) iter = FragPairVec.erase(iter);
		else iter++;
	}
}

bool CheckAlnCanCoverage(int rlen, vector<FragPair_t>& FragPairVec)
{
	bool bChecked = true;
	int TotalCovLength = 0;

	for (vector<FragPair_t>::iterator iter = FragPairVec.begin(); iter != FragPairVec.end(); iter++)
	{
		if (iter->rLen > 0) TotalCovLength += iter->rLen;
	}
	if (TotalCovLength != rlen) bChecked = false;

	return bChecked;
}

int CalFragPairMismatches(int len, string& str1, string& str2)
{
	int i, mismatch;

	for (mismatch = i = 0; i < len; i++)
	{
		if (str1[i] != str2[i]) mismatch++;
	}
	return mismatch;
}

int CalFragPairMatches(int len, string& str1, string& str2)
{
	int i, mismatch;

	for (mismatch = i = 0; i < len; i++)
	{
		if (str1[i] != str2[i]) mismatch++;
	}
	return (len - mismatch);
}

void ProcessNormalPair(char* seq, FragPair_t& fp)
{
	if (fp.rLen > 0 && fp.gLen > 0)
	{
		int mismatch = 0;

		fp.aln1.resize(fp.rLen); strncpy((char*)fp.aln1.c_str(), seq + fp.rPos, fp.rLen);
		fp.aln2.resize(fp.gLen); strncpy((char*)fp.aln2.c_str(), RefSequence + fp.gPos, fp.gLen);

		if (fp.rLen != fp.gLen || ((mismatch = CalFragPairMismatches(fp.rLen, fp.aln1, fp.aln2) > 3) && mismatch > (int)(fp.rLen*0.2)))
		{
			nw_alignment(fp.rLen, fp.aln1, fp.gLen, fp.aln2);
		}
	}
}

int EvaluateAlignmentScore(vector<FragPair_t>& FragPairVec)
{
	int len, score = 0;

	for (vector<FragPair_t>::iterator iter = FragPairVec.begin(); iter != FragPairVec.end(); iter++)
	{
		if (iter->bSimple) score += iter->rLen;
		else if ((len = (int)iter->aln1.length()) > 0) score += CalFragPairMatches(len, iter->aln1, iter->aln2);
	}
	return score;
}

int CheckAlignmentBounary(vector<FragPair_t>& FragPairVec)
{
	map<int64_t, int>::iterator iter1, iter2;

	iter1 = RefSeqLocMap.lower_bound(FragPairVec.begin()->gPos);
	iter2 = RefSeqLocMap.lower_bound((FragPairVec.rbegin()->gPos + FragPairVec.rbegin()->gLen - 1));

	if (iter1->second != iter2->second) return -1;
	else return iter1->second;
}

bool ProduceReadAlignment(ReadItem_t& read)
{
	vector<AlnCan_t>::iterator iter;
	int i, thr, ref_idx, FragPairNum;

	read.AlnSummary.score = read.AlnSummary.taxid = 0;
	thr = read.rlen - MaxMismatchNum;
	for (iter = read.AlnCanVec.begin(); iter != read.AlnCanVec.end(); iter++)
	{
		if (iter->score == 0) continue;

		sort(iter->FragPairVec.begin(), iter->FragPairVec.end(), CompByReadPos);
		if (RemoveOverlaps(iter->FragPairVec)) RemoveNullFragPairs(iter->FragPairVec);
		IdentifyNormalPairs(read.rlen, iter->FragPairVec);
		if ((ref_idx = CheckAlignmentBounary(iter->FragPairVec)) == -1) continue;
		for (FragPairNum = (int)iter->FragPairVec.size(), i = 0; i < FragPairNum; i++) if (!iter->FragPairVec[i].bSimple) ProcessNormalPair(read.seq, iter->FragPairVec[i]);
		//if (bDebugMode) ShowFragPairInfo(iter->FragPairVec);
		if ((iter->score = EvaluateAlignmentScore(iter->FragPairVec)) >= thr)
		{
			if (iter->score > read.AlnSummary.score)
			{
				read.AlnSummary.score = iter->score;
				read.AlnSummary.ref_idx = ref_idx;
				read.AlnSummary.taxid = RefTaxVec[ref_idx];
			}
			else if (iter->score == read.AlnSummary.score)
			{
				read.AlnSummary.ref_idx = -1;
				if (read.AlnSummary.taxid != RefTaxVec[ref_idx]) read.AlnSummary.taxid = Pairwise_LCA(read.AlnSummary.taxid, RefTaxVec[ref_idx]);
			}
		}
		else iter->score = 0;
	}
	return read.AlnSummary.score > 0;
}
