#include "structure.h"

#define MinClusterSize 1000000000

static pthread_mutex_t Lock;
map<int, int64_t> ClusterSizeMap;
static int ClusterSize, ClusterID;
vector<pair<string, int64_t> > ClusterSeqPathVec;

bool CompByClusterSize(const pair<string, int64_t>& p1, const pair<string, int64_t>& p2)
{
	return p1.second > p2.second;
}

int GetTaxIdByTaxLevel(int taxid)
{
	map<int, TaxItem_t>::iterator iter;
	
	iter = TaxMap.find(taxid);
	if (iter == TaxMap.end()) return 0;
	else
	{
		while (iter->second.rank % 10 != 0)
		{
			iter = TaxMap.find(iter->second.parent_taxid);
		}
		return iter->first;
	}
}

void Seq_Clustering(const char* clr_suffix)
{
	FILE* cluster_fp;
	char cluster_fn[24];
	string str, tax, tmp, cmd;
	vector<pair<int, int> > vec;
	int i, sid, cluster_id, taxid;
	map<int, vector<int> > ClusterMemberMap;

	ClusterSizeMap.clear();
	for (i = 0; i < (int)SeqVec.size(); i++)
	{
		if ((taxid = SeqVec[i].taxid) != 0)
		{
			cluster_id = GetTaxIdByTaxLevel(taxid);
			ClusterSizeMap[cluster_id] += (int)SeqVec[i].seq.length();
		}
	}
	while (true)
	{
		vec.clear();
		for (map<int, int64_t>::iterator iter = ClusterSizeMap.begin(); iter != ClusterSizeMap.end(); iter++)
		{
			if (iter->first == 0 || iter->second == 0) continue;
			if (iter->second < MinClusterSize && TaxMap[iter->first].rank <= 70)
			{
				vec.push_back(make_pair(iter->first, iter->second));
				iter->second = 0;
			}
		}
		if (vec.size() == 0) break;
		for (i = 0; i < (int)vec.size(); i++) 
		{
			taxid = GetTaxIdByTaxLevel(TaxMap[vec[i].first].parent_taxid);
			ClusterSizeMap[taxid] += vec[i].second;
		}
	}
	//for (map<int, int64_t>::iterator iter = ClusterSizeMap.begin(); iter != ClusterSizeMap.end(); iter++)
	//{
	//	if (iter->second == 0) continue;
	//	printf("cluster %d: %d M\n", iter->first, iter->second / 1000000);
	//}
	ClusterSeqPathVec.clear(); ClusterMemberMap.clear();
	fprintf(stderr, "Cluster %d sequences...\n", (int)SeqVec.size());
	for (i = 0; i < (int)SeqVec.size(); i++)
	{
		if ((taxid = SeqVec[i].taxid) != 0)
		{
			cluster_id = GetTaxIdByTaxLevel(taxid);
			while (ClusterSizeMap[cluster_id] == 0) cluster_id = GetTaxIdByTaxLevel(TaxMap[cluster_id].parent_taxid);

			ClusterMemberMap[cluster_id].push_back(i);
		}
	}
	for (map<int, vector<int> >::iterator iter = ClusterMemberMap.begin(); iter != ClusterMemberMap.end(); iter++)
	{
		cluster_id = iter->first; 
		sprintf(cluster_fn, "%s/%d.%s", OutputFolder.c_str(), cluster_id, clr_suffix);
		ClusterSeqPathVec.push_back(make_pair(cluster_fn, ClusterSizeMap[cluster_id])); cluster_fp = fopen(cluster_fn, "w");
		for (i = 0; i < (int)iter->second.size(); i++)
		{
			sid = iter->second[i];
			fprintf(cluster_fp, ">%s\n%s\n", SeqVec[sid].header.c_str(), SeqVec[sid].seq.c_str());
		}
		fclose(cluster_fp);
	}
	sort(ClusterSeqPathVec.begin(), ClusterSeqPathVec.end(), CompByClusterSize);
}

static void *Make_BWT_index(void *arg)
{
	int job_id;
	char cmd[256];
	string fn, IdxPrefix;

	while (true)
	{
		pthread_mutex_lock(&Lock);
		job_id = ClusterID++;
		pthread_mutex_unlock(&Lock);

		if (job_id >= ClusterSize) break;
		fn = ClusterSeqPathVec[job_id].first; IdxPrefix = fn.substr(0, fn.find_last_of('.'));
		pthread_mutex_lock(&Lock);
		fprintf(stderr, "Build index for %s (%d Mbp) (%d/%d)\n", fn.c_str(), (int)(ClusterSeqPathVec[job_id].second / 1000000), job_id + 1, ClusterSize);
		pthread_mutex_unlock(&Lock);
		sprintf(cmd, "./bwt_index %s %s > /dev/null", fn.c_str(), IdxPrefix.c_str()); system(cmd);
	}
	return (void*)(1);
}

void Make_DB_Index()
{
	int i;

	ClusterID = 0; ClusterSize = (int)ClusterSeqPathVec.size();
	pthread_t *ThreadArr = new pthread_t[iThreadNum];
	for (i = 0; i < iThreadNum; i++) pthread_create(&ThreadArr[i], NULL, Make_BWT_index, NULL);
	for (i = 0; i < iThreadNum; i++) pthread_join(ThreadArr[i], NULL);
	fprintf(stderr, "\n");
}
