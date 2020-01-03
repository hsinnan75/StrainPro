#include "structure.h"

vector<int> RefTaxVec;

static bwtint_t fread_fix(FILE *fp, bwtint_t size, void *a)
{
	const unsigned int bufsize = 0x1000000; // 16M block
	bwtint_t offset = 0;
	while (size) {
		bwtint_t x = bufsize < size ? bufsize : size;
		x = fread(((char*)a + offset), 1, x, fp);
		size -= x; offset += x;
	}
	return offset;
}

void bwt_restore_sa(const char *fn, bwt_t *bwt)
{
	char skipped[256];
	FILE *fp;
	bwtint_t primary;

	fp = fopen(fn, "rb");
	fread(&primary, sizeof(bwtint_t), 1, fp);
	//xassert(primary == bwt->primary, "SA-BWT inconsistency: primary is not the same.");
	fread(skipped, sizeof(bwtint_t), 4, fp); // skip
	fread(&bwt->sa_intv, sizeof(bwtint_t), 1, fp);
	fread(&primary, sizeof(bwtint_t), 1, fp);
	//xassert(primary == bwt->seq_len, "SA-BWT inconsistency: seq_len is not the same.");

	bwt->n_sa = (bwt->seq_len + bwt->sa_intv) / bwt->sa_intv;
	bwt->sa = (bwtint_t*)calloc(bwt->n_sa, sizeof(bwtint_t));
	bwt->sa[0] = -1;

	fread_fix(fp, sizeof(bwtint_t) * (bwt->n_sa - 1), bwt->sa + 1);
	fclose(fp);
}

bntseq_t *bns_restore_core(const char *ann_filename, const char* amb_filename, const char* pac_filename)
{
	char str[1024];
	FILE *fp;
	const char *fname;
	bntseq_t *bns;
	long long xx;
	int i;
	int scanres;
	bns = (bntseq_t*)calloc(1, sizeof(bntseq_t));
	{ // read .ann
		fp = fopen(fname = ann_filename, "r");
		scanres = fscanf(fp, "%lld%d%u", &xx, &bns->n_seqs, &bns->seed);
		bns->l_pac = xx;
		bns->anns = (bntann1_t*)calloc(bns->n_seqs, sizeof(bntann1_t));
		for (i = 0; i < bns->n_seqs; ++i) {
			bntann1_t *p = bns->anns + i;
			char *q = str;
			int c;
			// read gi and sequence name
			scanres = fscanf(fp, "%u%s", &p->gi, str);
			p->name = strdup(str);
			// read fasta comments 
			while (str - q < (int)(sizeof(str) - 1) && (c = fgetc(fp)) != '\n' && c != EOF) *q++ = c;
			while (c != '\n' && c != EOF) c = fgetc(fp);

			*q = 0;
			if (q - str > 1) p->anno = strdup(str + 1); // skip leading space
			else p->anno = strdup("");
			// read the rest
			scanres = fscanf(fp, "%lld%d%d", &xx, &p->len, &p->n_ambs);
			p->offset = xx;
		}
		fclose(fp);
	}
	{ // read .amb
		int64_t l_pac;
		int32_t n_seqs;
		fp = fopen(fname = amb_filename, "r");
		scanres = fscanf(fp, "%lld%d%d", &xx, &n_seqs, &bns->n_holes);
		l_pac = xx;
		//xassert(l_pac == bns->l_pac && n_seqs == bns->n_seqs, "inconsistent .ann and .amb files.");
		bns->ambs = bns->n_holes? (bntamb1_t*)calloc(bns->n_holes, sizeof(bntamb1_t)) : 0;
		for (i = 0; i < bns->n_holes; ++i) {
			bntamb1_t *p = bns->ambs + i;
			scanres = fscanf(fp, "%lld%d%s", &xx, &p->len, str);
			p->offset = xx;
			p->amb = str[0];
		}
		fclose(fp);
	}
	{ // open .pac
		bns->fp_pac = fopen(pac_filename, "rb");
	}
	return bns;
}

void bwt_gen_cnt_table(bwt_t *bwt)
{
	int i, j;
	for (i = 0; i != 256; ++i) {
		uint32_t x = 0;
		for (j = 0; j != 4; ++j)
			x |= (((i&3) == j) + ((i>>2&3) == j) + ((i>>4&3) == j) + (i>>6 == j)) << (j<<3);
		bwt->cnt_table[i] = x;
	}
}

bwt_t *bwt_restore_bwt(const char *fn)
{
	bwt_t *bwt;
	FILE *fp;

	bwt = (bwt_t*)calloc(1, sizeof(bwt_t));
	fp = fopen(fn, "rb");
	fseek(fp, 0, SEEK_END);
	bwt->bwt_size = (ftell(fp) - sizeof(bwtint_t) * 5) >> 2;
	bwt->bwt = (uint32_t*)calloc(bwt->bwt_size, 4);
	fseek(fp, 0, SEEK_SET);
	fread(&bwt->primary, sizeof(bwtint_t), 1, fp);
	fread(bwt->L2 + 1, sizeof(bwtint_t), 4, fp);
	fread_fix(fp, bwt->bwt_size << 2, bwt->bwt);
	bwt->seq_len = bwt->L2[4];
	fclose(fp);
	bwt_gen_cnt_table(bwt);

	return bwt;
}

bwt_t *bwa_idx_load_bwt(const char* IndexPrefix)
{
	char *tmp;
	bwt_t *bwt;

	tmp = (char*)calloc(strlen(IndexPrefix) + 5, 1);
	strcat(strcpy(tmp, IndexPrefix), ".bwt"); // FM-index
	bwt = bwt_restore_bwt(tmp);
	strcat(strcpy(tmp, IndexPrefix), ".sa");  // partial suffix array (SA)
	bwt_restore_sa(tmp, bwt);
	free(tmp);

	return bwt;
}

bntseq_t *bns_restore(const char* IndexPrefix)
{  
	char ann_filename[256], amb_filename[256], pac_filename[256];
	strcat(strcpy(ann_filename, IndexPrefix), ".ann");
	strcat(strcpy(amb_filename, IndexPrefix), ".amb");
	strcat(strcpy(pac_filename, IndexPrefix), ".pac");
	return bns_restore_core(ann_filename, amb_filename, pac_filename);
}

bwaidx_t *bwa_idx_load(const char* IndexPrefix)
{
	bwaidx_t *idx;

	//fprintf(stderr, "\tLoad the reference index files...");
	idx = (bwaidx_t*)calloc(1, sizeof(bwaidx_t));
	idx->bwt = bwa_idx_load_bwt(IndexPrefix);
	idx->bns = bns_restore(IndexPrefix);
	idx->pac = (uint8_t*)calloc(idx->bns->l_pac/4+1, 1);
	//fprintf(stderr, "Done!\n");

	return idx;
}

void bwt_destroy(bwt_t *bwt)
{
	if (bwt == 0) return;
	free(bwt->sa); free(bwt->bwt);
	free(bwt);
}

void bns_destroy(bntseq_t *bns)
{
	if (bns == 0) return;
	else {
		int i;
		if (bns->fp_pac) fclose(bns->fp_pac);
		free(bns->ambs);
		for (i = 0; i < bns->n_seqs; ++i) {
			free(bns->anns[i].name);
			free(bns->anns[i].anno);
		}
		free(bns->anns);
		free(bns);
	}
}

void bwa_idx_destroy(bwaidx_t *idx)
{
	if (idx == 0) return;
	if (idx->bwt) bwt_destroy(idx->bwt);
	if (idx->bns) bns_destroy(idx->bns);
	if (idx->pac) free(idx->pac);
	free(idx);
}

int GetTaxFromHeader(string header)
{
	string tax;
	int p1, p2, taxid;

	p1 = header.find("taxid|") + 6; p2 = header.find_first_of('|', p1);
	tax = header.substr(p1, p2 - p1); taxid = atoi(tax.c_str());

	return taxid;
}

int GetSizeFromHeader(string header)
{
	string str;

	str = header.substr(str.find_last_of('=') + 1);
	return atoi(str.c_str());
}

void GetRefSeqTaxInfo()
{
	RefTaxVec.resize(RefIdx->bns->n_seqs);
	for (int i = 0; i < RefIdx->bns->n_seqs; i++) RefTaxVec[i] = GetTaxFromHeader(RefIdx->bns->anns[i].name);
}

void RestoreReferenceInfo()
{
	int i;
	int64_t gPos = 0;

	//fprintf(stderr, "\tLoad the %d reference sequences (%lldM bp)", iRefSeqNum, RefSeqSize/1000000);
	fseek(RefIdx->bns->fp_pac, 0, SEEK_SET);
	fread(RefIdx->pac, 1, RefSeqSize / 4 + 1, RefIdx->bns->fp_pac);

	RefSeqLocMap.clear(); 
	for (i = 0; i < iRefSeqNum; i++)
	{
		RefSeqLocMap.insert(make_pair((DoubleRefSeqSize - 1 - gPos), i));
		RefSeqLocMap.insert(make_pair(((gPos += RefIdx->bns->anns[i].len) - 1), i));
	}
	//fprintf(stderr, "\n");
}

void InitializeReferenceData(const char* IndexPrefix)
{
	if ((RefIdx = bwa_idx_load(IndexPrefix)) == 0)
	{
		fprintf(stderr, "Error! Index files are corrupt!\n");
		exit(1);
	}
	else
	{
		//fprintf(stderr, "\tLoading index files..."); fflush(stderr);
		Refbwt = RefIdx->bwt; iRefSeqNum = RefIdx->bns->n_seqs;
		RefSeqSize = RefIdx->bns->l_pac; DoubleRefSeqSize = (RefSeqSize << 1);
		RestoreReferenceInfo();
		GetRefSeqTaxInfo();
	}
}

void *IdvLoadReferenceSequences(void *arg)
{
	int base, *my_id;
	int64_t fPos, rPos;

	my_id = (int*)arg;
	for (fPos = *my_id, rPos = DoubleRefSeqSize - fPos - 1; fPos < RefSeqSize; fPos += iThreadNum, rPos -= iThreadNum)
	{
		base = RefIdx->pac[fPos >> 2] >> ((~fPos & 3) << 1) & 3;
		switch (base)
		{
		case 0: RefSequence[fPos] = 'A'; RefSequence[rPos] = 'T'; break;
		case 1: RefSequence[fPos] = 'C'; RefSequence[rPos] = 'G'; break;
		case 2: RefSequence[fPos] = 'G'; RefSequence[rPos] = 'C'; break;
		case 3: RefSequence[fPos] = 'T'; RefSequence[rPos] = 'A'; break;
		default:RefSequence[fPos] = RefSequence[rPos] = 'N';
		}
	}
	return (void*)(1);
}

void RestoreReferenceSequences()
{
	int i, *JobIDArr = new int[iThreadNum];

	RefSequence = new char[DoubleRefSeqSize + 1]; RefSequence[DoubleRefSeqSize] = '\0';
	pthread_t *ThreadArr = new pthread_t[iThreadNum];
	for (i = 0; i < iThreadNum; i++)
	{
		JobIDArr[i] = i;
		pthread_create(&ThreadArr[i], NULL, IdvLoadReferenceSequences, JobIDArr + i);
	}
	for (i = 0; i < iThreadNum; i++) pthread_join(ThreadArr[i], NULL);

	delete[] ThreadArr; delete[] JobIDArr;
}
