#include "structure.h"

typedef struct
{
	int freq;
	int depth_size;
} TaxMappedInfo_t;

map<int, int64_t> TaxLevelFreqMap;
map<int, TaxMappedInfo_t> TaxMappedInfoMap;

void EstTaxDepth()
{
	int64_t rid;
	map<int, int> RefIdxFreqMap; //first=ref_idx, second=frequency
	map<int, int>::iterator iter;
	int i, ref_idx, taxid, size, thr;

	fprintf(stderr, "Get reference seqeunce's taxnomic information\n");
	for (IndexID = 0; IndexID < IndexNum; IndexID++)
	{
		RefIdxFreqMap.clear();
		fprintf(stderr, "\33[2K\r\tLoad BWT index: %s (%d/%d)", IndexPrefixVec[IndexID].c_str(), IndexID + 1, IndexNum);
		InitializeReferenceData(IndexPrefixVec[IndexID].c_str());

		for (rid = 0; rid < iTotalReadNum; rid++)
		{
			if (MapInfoVec[rid].score > 0 && MapInfoVec[rid].idx_id == IndexID && MapInfoVec[rid].ref_idx != -1) RefIdxFreqMap[MapInfoVec[rid].ref_idx]++;
		}
		for (iter = RefIdxFreqMap.begin(); iter != RefIdxFreqMap.end(); iter++)
		{
			taxid = RefTaxVec[iter->first]; size = (RefIdx->bns->anns[iter->first].len - 100); 
			if (iter->second > (int)(size*minSeqRatio))
			{
				TaxMappedInfoMap[taxid].freq += iter->second;
				TaxMappedInfoMap[taxid].depth_size += size;
			}
		}
		bwa_idx_destroy(RefIdx);
	}
	fprintf(stderr, "\n");
}

int FindNextLevelTaxID(int taxid, int TaxLevel)
{
	int target_taxLevel, parentid;
	map<int, TaxItem_t>::iterator iter;
	
	if ((iter = TaxMap.find(taxid)) == TaxMap.end()) return 1;
	else parentid = TaxMap[taxid].parent_taxid;

	if (TaxLevel == 70) return 1;
	else if (TaxLevel % 10 == 5) target_taxLevel += 5;
	else target_taxLevel += 10;
	
	while(TaxMap[parentid].rank % 10 == 5 || TaxMap[parentid].rank < target_taxLevel) parentid = TaxMap[parentid].parent_taxid;

	return parentid;
}

void TaxTreeTraversal()
{
	map<int, TaxMappedInfo_t>::iterator iter;
	int taxid, parentid, freq, size, taxLevel;

	for (taxLevel = 5; taxLevel <= 70; taxLevel+=5)
	{
		for (iter = TaxMappedInfoMap.begin(); iter != TaxMappedInfoMap.end(); iter++)
		{
			if (iter->second.freq >= minFrequency && TaxMap[(taxid = iter->first)].rank == taxLevel)
			{
				freq = iter->second.freq; size = iter->second.depth_size;
				if (taxLevel == 5) TaxLevelFreqMap[taxLevel] += (int)(101 * (1.0*freq / size));
				else TaxLevelFreqMap[taxLevel] += freq;
				parentid = FindNextLevelTaxID(taxid, taxLevel); TaxMappedInfoMap[parentid].freq += freq; TaxMappedInfoMap[parentid].depth_size += size;
			}
		}
	}
}

void CreateMappingReport()
{
	float conf, abundance;
	FILE *OutputFileHandler;
	char OutputFileName[256];
	map<int, TaxMappedInfo_t>::iterator iter;
	int taxid, parentid, freq, size, depth, taxLevel;

	OutputFileHandler = fopen(OutputFilename, "w");
	fprintf(OutputFileHandler, "%-20s\t%-20s\t%-20s\t%-20s\t%-20s\n", "#TaxID", "#Read_count", "#Depth", "#Relative_abundance", "#Confidence_score");
	for (taxLevel = 5; taxLevel <= 70;)
	{
		fprintf(OutputFileHandler, "@TaxRank:%s\n", ShowTaxRank(taxLevel).c_str());
		for (iter = TaxMappedInfoMap.begin(); iter != TaxMappedInfoMap.end(); iter++)
		{
			if (iter->second.freq >= minFrequency && TaxMap[(taxid = iter->first)].rank == taxLevel)
			{
				freq = iter->second.freq; size = iter->second.depth_size; if (size < (1000 * (taxLevel / 10))) size = 1000 * (taxLevel / 10);
				depth = (int)(101 * (1.0*freq / size));
				conf = taxLevel == 5 ? 1.0*freq / size : (taxLevel / 10.0) * freq / size; if (conf > 1.0) conf = 1.0;
				abundance = 100*(taxLevel == 5 ? (float)depth / TaxLevelFreqMap[taxLevel] : (float)freq / TaxLevelFreqMap[taxLevel]);
				//if (conf > 0.5 || (freq >= minFrequency && depth >= minDepth)) fprintf(OutputFileHandler, "%-20d\t%-20d\t%-20d\t%-20f\n", taxid, freq, depth, conf); //taxid, frequency, est_depth, confidence_score
				//if((taxLevel == 5 && freq >= minFrequency && depth >= minDepth) || (taxLevel >= 10 && freq >= 1000 && conf > 0.01))
				fprintf(OutputFileHandler, "%-20d\t%-20d\t%-20d\t%-20.2f\t%-20f\n", taxid, freq, depth, abundance, conf);
			}
		}
		if (taxLevel == 5) taxLevel = 10;
		else taxLevel += 10;
	}
	fclose(OutputFileHandler);
	fprintf(stderr, "Write typing result to [%s]\n", OutputFilename);
}

void MetaTyping()
{
	EstTaxDepth();
	TaxTreeTraversal();
	CreateMappingReport();
}
