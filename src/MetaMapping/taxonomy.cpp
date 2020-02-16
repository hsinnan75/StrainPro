#include "structure.h"

map<int, TaxItem_t> TaxMap;
map<string, int> TaxRankMap;

void InitializeTaxRankMap()
{
	TaxRankMap.insert(make_pair("species", 10));
	TaxRankMap.insert(make_pair("genus", 20));
	TaxRankMap.insert(make_pair("family", 30));
	TaxRankMap.insert(make_pair("order", 40));
	TaxRankMap.insert(make_pair("class", 50));
	TaxRankMap.insert(make_pair("phylum", 60));
	TaxRankMap.insert(make_pair("superkingdom", 70));
}

int FindTaxRank(string rank)
{
	map<string, int>::iterator iter;

	iter = TaxRankMap.find(rank);
	if (iter == TaxRankMap.end()) return 0;
	else return iter->second;
}

int FindParentTaxidRank(int parentid)
{
	map<int, TaxItem_t>::iterator iter;

	while (true)
	{
		iter = TaxMap.find(parentid);
		if (iter->second.rank != 0 && iter->second.rank % 10 == 0) return iter->second.rank;
		else parentid = iter->second.parent_taxid;
	}
}

void GetTaxInfomation()
{
	fstream file;
	TaxItem_t TaxItem;
	string fn, str, tax, rank, tmp;
	map<int, TaxItem_t>::iterator iter;
	int p1, p2, rankid, taxid, parentid;

	InitializeTaxRankMap();

	fprintf(stderr, "Load taxonomy information.\n");
	fn = TaxonomyDir + "/nodes.dmp"; file.open(fn.c_str(), ios_base::in); // tax_id, parent tax_id, rank
	if (!file.is_open())
	{
		fprintf(stderr, "Error! File: %s is not accesible\n", fn.c_str());
		exit(1);
	}
	while (!file.eof())
	{
		getline(file, str); if (str == "") break;

		p1 = 0; p2 = str.find_first_of('|', p1 + 1); tmp = str.substr(p1, p2 - p1 - 1); taxid = atoi(tmp.c_str());
		p1 = p2 + 2; p2 = str.find_first_of('|', p1); tmp = str.substr(p1, p2 - 1 - p1); parentid = atoi(tmp.c_str());
		p1 = p2 + 2; p2 = str.find_first_of('|', p1); rank = str.substr(p1, p2 - 1 - p1);
		//printf("%s\ntaxid=[%d], pid=[%d], rank=[%s]\n", str.substr(0, p2).c_str(), taxid, parentid, rank.c_str());
		TaxItem.parent_taxid = parentid;
		if (taxid == parentid) TaxItem.rank = 80; else TaxItem.rank = FindTaxRank(rank);
		TaxMap.insert(make_pair(taxid, TaxItem));
	}
	file.close();
	for (iter = TaxMap.begin(); iter != TaxMap.end(); iter++)
	{
		if (iter->second.rank == 0)
		{
			if (TaxMap.find(iter->second.parent_taxid) == TaxMap.end())
			{
				fprintf(stderr, "\nError! Cannot find taxid:%d in the taxonomy dump files\n", iter->second.parent_taxid);
				exit(1);
			}
			rankid = FindParentTaxidRank(iter->second.parent_taxid);
			if (rankid > 0) iter->second.rank = rankid - 5;
			else iter->second.rank = 0;
			//printf("%d: %d\n", iter->first, iter->second.rank);
		}
	}
	fn = TaxonomyDir + "/merged.dmp"; file.clear(); file.open(fn.c_str(), ios_base::in); // tax_id, parent tax_id, rank
	if (!file.is_open())
	{
		fprintf(stderr, "Error! File: %s is not accesible\n", fn.c_str());
		exit(1);
	}
	while (!file.eof())
	{
		getline(file, str); if (str == "") break;

		p1 = 0; p2 = str.find_first_of('|', p1 + 1); tmp = str.substr(p1, p2 - p1 - 1); taxid = atoi(tmp.c_str());
		p1 = p2 + 2; p2 = str.find_first_of('|', p1); tmp = str.substr(p1, p2 - 1 - p1); parentid = atoi(tmp.c_str());
		if ((iter = TaxMap.find(parentid)) != TaxMap.end()) TaxMap.insert(make_pair(taxid, iter->second));
	}
	file.close();
}

int Pairwise_LCA(int taxid1, int taxid2)
{
	map<int, TaxItem_t>::iterator iter1, iter2;

	if ((iter1 = TaxMap.find(taxid1)) == TaxMap.end())
	{
		fprintf(stderr, "Cannot find taxid:%d in the taxonomy dump files\n", taxid1);
		exit(1);
	}
	if ((iter2 = TaxMap.find(taxid2)) == TaxMap.end())
	{
		fprintf(stderr, "Cannot find taxid:%d in the taxonomy dump files\n", taxid2);
		exit(1);
	}
	while (true)
	{
		while (iter1->second.rank != iter2->second.rank)
		{
			while (iter1->second.rank < iter2->second.rank)
			{
				taxid1 = iter1->second.parent_taxid;
				if ((iter1 = TaxMap.find(taxid1)) == TaxMap.end())
				{
					fprintf(stderr, "Cannot find taxid:%d in the taxonomy dump files\n", taxid1);
					exit(1);
				}
			}
			while (iter2->second.rank < iter1->second.rank)
			{
				taxid2 = iter2->second.parent_taxid;
				if ((iter2 = TaxMap.find(taxid2)) == TaxMap.end())
				{
					fprintf(stderr, "Cannot find taxid:%d in the taxonomy dump files\n", taxid2);
					exit(1);
				}
			}
		}
		if (iter1 != iter2)
		{
			taxid1 = iter1->second.parent_taxid;
			if ((iter1 = TaxMap.find(taxid1)) == TaxMap.end())
			{
				fprintf(stderr, "Cannot find taxid:%d in the taxonomy dump files\n", taxid1);
				exit(1);
			}
			taxid2 = iter2->second.parent_taxid;
			if ((iter2 = TaxMap.find(taxid2)) == TaxMap.end())
			{
				fprintf(stderr, "Cannot find taxid:%d in the taxonomy dump files\n", taxid2);
				exit(1);
			}
		}
		else break;
	}
	//printf("LCA=%d\n\n", taxid1);
	return taxid1;
}


int Find_LCA(vector<int>& TaxVec)
{
	int LCA = 0;
	vector<int>::iterator iter;

	for (iter = TaxVec.begin(); iter != TaxVec.end(); iter++)
	{
		//printf("taxid=%d, parentid=%d, rank=%d\n", *iter, TaxMap[*iter].parent_taxid, TaxMap[*iter].rank);
		if (*iter == 0) continue;
		else if (LCA == 0) LCA = *iter;
		else LCA = Pairwise_LCA(LCA, *iter);
	}
	return LCA;
}
