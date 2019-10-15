#include "structure.h"

bwt_t *Refbwt;
bwaidx_t *RefIdx;
char *RefSequence;
time_t MappingTime;
int IndexID, IndexNum;
FILE *ReadFileHandler;
gzFile gzReadFileHandler;
vector<string> ReadSeqVec;
static pthread_mutex_t Lock;
vector<MapInfo_t> MapInfoVec;
bool FastQFormat, gzCompressed;
int64_t ReadSeqVecSize, iTotalReadNum, RefSeqSize, DoubleRefSeqSize;

bool CompByPosDiff(const FragPair_t& p1, const FragPair_t& p2)
{
	if (p1.PosDiff == p2.PosDiff) return p1.rPos < p2.rPos;
	else return p1.PosDiff < p2.PosDiff;
}

bool CompByAlnCanScore(const AlnCan_t& p1, const AlnCan_t& p2)
{
	if (p1.score == p2.score) return p1.FragPairVec[0].PosDiff < p2.FragPairVec[0].PosDiff;
	else return p1.score > p2.score;
}

void CopyReadSeq(int ReadNum, ReadItem_t* ReadArr)
{
	for (int i = 0; i < ReadNum; i++) ReadSeqVec.push_back(ReadArr[i].seq);
}

int GetNextChunkFromMemory(int64_t idx, ReadItem_t* ReadArr)
{
	int iCount = 0;
	while(idx < ReadSeqVecSize)
	{
		ReadArr[iCount].AlnSummary.score = 0;
		ReadArr[iCount].rlen = ReadSeqVec[idx].length();
		ReadArr[iCount].seq = (char*)ReadSeqVec[idx].c_str();
		idx++; if (++iCount == ReadChunkSize) break;
	}
	return iCount;
}

vector<FragPair_t> IdentifySimplePairs(int rlen, char* seq)
{
	int i, pos;
	FragPair_t FragPair;
	vector<FragPair_t> SPvec;
	bwtSearchResult_t bwtSearchResult;

	FragPair.bSimple = true; pos = 0;
	while (pos < rlen)
	{
		if (nst_nt4_table[(uint8_t)seq[pos]] > 3) pos++;
		else
		{
			bwtSearchResult = BWT_Search(seq, pos, rlen);
			if (bwtSearchResult.freq > 0)
			{
				FragPair.rPos = pos; FragPair.rLen = FragPair.gLen = bwtSearchResult.len;
				for (i = 0; i != bwtSearchResult.freq; i++)
				{
					FragPair.PosDiff = (FragPair.gPos = bwtSearchResult.LocArr[i]) - FragPair.rPos;
					if (FragPair.PosDiff > 0) SPvec.push_back(FragPair);
				}
				delete[] bwtSearchResult.LocArr;
			}
			pos += (bwtSearchResult.len + 1);
		}
	}
	sort(SPvec.begin(), SPvec.end(), CompByPosDiff);

	FragPair.rLen = FragPair.gLen = 0; FragPair.gPos = FragPair.PosDiff = DoubleRefSeqSize;
	SPvec.push_back(FragPair); // add a terminal fragment pair

	return SPvec;
}

vector<AlnCan_t> SimplePairClustering(int rlen, vector<FragPair_t>& SimplePairVec)
{
	AlnCan_t AlnCan;
	vector<AlnCan_t> AlnCanVec;
	int i, j, score_thr, HeadIdx, num;

	if ((num = (int)SimplePairVec.size()) == 1) return AlnCanVec;
	score_thr = (rlen >> 1); HeadIdx = 0; AlnCan.score = SimplePairVec[0].rLen;
	for (i = 0, j = 1; j < num; i++, j++)
	{
		if (SimplePairVec[j].PosDiff - SimplePairVec[i].PosDiff > MaxPosDiffGap)
		{
			if (AlnCan.score > score_thr)
			{
				copy(SimplePairVec.begin() + HeadIdx, SimplePairVec.begin() + j, back_inserter(AlnCan.FragPairVec));
				AlnCanVec.push_back(AlnCan); AlnCan.FragPairVec.clear();
			}
			HeadIdx = j; AlnCan.score = SimplePairVec[j].rLen;
		}
		else AlnCan.score += SimplePairVec[j].rLen;
	}
	return AlnCanVec;
}

void RemoveRedundantAlnCan(vector<AlnCan_t>& AlnCanVec)
{
	if ((int)AlnCanVec.size() > 1)
	{
		int maxScore = 0;
		vector<AlnCan_t>::iterator iter;

		for (iter = AlnCanVec.begin(); iter != AlnCanVec.end(); iter++) if (iter->score > maxScore) maxScore = iter->score;
		maxScore -= 20;
		for (iter = AlnCanVec.begin(); iter != AlnCanVec.end(); iter++) if (iter->score < maxScore) iter->score = 0;
	}
}

//char* FindOriginalReadHeader(char* header)
//{
//	int i, len;
//
//	len = strlen(header);
//	for (i = len - 1; i > 0; i--)
//	{
//		if (header[i] == ' ') return header + i + 1;
//	}
//	return header;
//}

void UpdateMappingInfo(int64_t rid_base, int ReadNum, ReadItem_t* ReadArr)
{
	int64_t rid;

	for (int i = 0; i != ReadNum; i++)
	{
		if (ReadArr[i].AlnSummary.score == 0) continue;

		rid = rid_base + i;
		if (ReadArr[i].AlnSummary.score > MapInfoVec[rid].score)
		{
			MapInfoVec[rid].score = ReadArr[i].AlnSummary.score;
			MapInfoVec[rid].taxid = ReadArr[i].AlnSummary.taxid;
			MapInfoVec[rid].idx_id = IndexID; MapInfoVec[rid].ref_idx = ReadArr[i].AlnSummary.ref_idx;
		}
		else if (ReadArr[i].AlnSummary.score == MapInfoVec[rid].score)
		{
			MapInfoVec[rid].taxid = Pairwise_LCA(MapInfoVec[rid].taxid, ReadArr[i].AlnSummary.taxid); // multiple hits
			MapInfoVec[rid].idx_id = -1; MapInfoVec[rid].ref_idx = -1;
		}
	}
}

void *ReadMapping(void *arg)
{
	int i, ReadNum;
	int64_t rid_base;
	ReadItem_t* ReadArr = NULL;
	vector<FragPair_t> SimplePairVec;
	map<int64_t, int>::iterator RefIter;

	ReadArr = new ReadItem_t[ReadChunkSize];
	while (true)
	{
		pthread_mutex_lock(&Lock);
		if (bFastMode && IndexID > 0) ReadNum = GetNextChunkFromMemory(iTotalReadNum, ReadArr); // memory
		else if (gzCompressed) ReadNum = gzGetNextChunk(gzReadFileHandler, ReadArr); // gz file
		else ReadNum = GetNextChunk(ReadFileHandler, ReadArr); //fq or fa file

		rid_base = iTotalReadNum; iTotalReadNum += ReadNum;
		if (IndexID == 0 && ReadNum > 0) MapInfoVec.resize(iTotalReadNum);
		if (IndexID == 0 && bFastMode) CopyReadSeq(ReadNum, ReadArr);
		fprintf(stderr, "\33[2K\r\t%lld reads have been aligned in %lld seconds", (long long)rid_base, (long long)(time(NULL) - MappingTime));
		pthread_mutex_unlock(&Lock);

		if (ReadNum == 0) break;
		for (i = 0; i != ReadNum; i++)
		{
			SimplePairVec = IdentifySimplePairs(ReadArr[i].rlen, ReadArr[i].seq);
			ReadArr[i].AlnCanVec = SimplePairClustering(ReadArr[i].rlen, SimplePairVec);
			RemoveRedundantAlnCan(ReadArr[i].AlnCanVec);
			//if (ReadArr[i].AlnCanVec.size() > 0) printf("%s\nReadCluster\n%s\n", string().assign(50, '*').c_str(), ReadArr[i].seq), ShowFragPairCluster(ReadArr[i].AlnCanVec);
			ProduceReadAlignment(ReadArr[i]);
			//printf("%d: score=%d, best_hit_idx=%d\n", i, ReadArr[i].AlnSummary.score, ReadArr[i].AlnSummary.best_hit_idx);
		}
		if (IndexID == 0)
		{
			pthread_mutex_lock(&Lock);
			UpdateMappingInfo(rid_base, ReadNum, ReadArr); // update MapInfoArr
			pthread_mutex_unlock(&Lock);
		}
		else
		{
			UpdateMappingInfo(rid_base, ReadNum, ReadArr); // update MapInfoArr
		}
		//for (i = 0; i != ReadNum; i++) delete[] ReadArr[i].header;
		if (!bFastMode || IndexID == 0) for (i = 0; i != ReadNum; i++) delete[] ReadArr[i].seq;
	}
	delete[] ReadArr;

	return (void*)(1);
}

void MetaMapping()
{
	int i, LibraryID;
	pthread_t *ThreadArr = new pthread_t[iThreadNum];

	IndexNum = (int)IndexPrefixVec.size(); MapInfoVec.clear();
	for (IndexID = 0; IndexID < IndexNum; IndexID++)
	{
		fprintf(stderr, "Load BWT index: %s (%d/%d)", IndexPrefixVec[IndexID].c_str(), IndexID + 1, IndexNum); 
		ReadFileHandler = NULL; gzReadFileHandler = NULL; iTotalReadNum = 0; MappingTime = time(NULL);
		InitializeReferenceData(IndexPrefixVec[IndexID].c_str()); RestoreReferenceSequences();
		fprintf(stderr, "\n"); fflush(stderr);
		for (LibraryID = 0; LibraryID < (int)ReadLibraryVec.size(); LibraryID++)
		{
			if (IndexID == 0 || bFastMode == false)
			{
				if (ReadLibraryVec[LibraryID].substr(ReadLibraryVec[LibraryID].find_last_of('.') + 1) == "gz") gzCompressed = true;
				else gzCompressed = false;
				FastQFormat = CheckReadFormat(ReadLibraryVec[LibraryID].c_str());
				if (gzCompressed) gzReadFileHandler = gzopen(ReadLibraryVec[LibraryID].c_str(), "rb");
				else ReadFileHandler = fopen(ReadLibraryVec[LibraryID].c_str(), "r");
			}
			for (i = 0; i < iThreadNum; i++) pthread_create(&ThreadArr[i], NULL, ReadMapping, NULL);
			for (i = 0; i < iThreadNum; i++) pthread_join(ThreadArr[i], NULL);

			if (IndexID == 0 || bFastMode == false)
			{
				if (gzCompressed) gzclose(gzReadFileHandler); else fclose(ReadFileHandler);
			}
		}
		if (bFastMode && IndexID == 0) ReadSeqVecSize = iTotalReadNum;
		fprintf(stderr, "\r\t%lld reads have been processed (%.2f%% mapped) in %lld seconds\n", (long long)iTotalReadNum, 100*(1.0*FindMappedReadNum()/ iTotalReadNum), (long long)(time(NULL) - MappingTime)); fflush(stderr);
		bwa_idx_destroy(RefIdx); delete[] RefSequence;
	}
	delete[] ThreadArr;
}
