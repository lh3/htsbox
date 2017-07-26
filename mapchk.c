#include <limits.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "sam.h"
#include "faidx.h"

void *bed_read(const char *fn);
int bed_overlap(const void *_h, const char *chr, int beg, int end);
void bed_destroy(void *_h);

typedef struct {
	int64_t q[94][5][8]; // qual, read base, ref base (5=ins, 6=del, 7=clip)
} errstat_t;

static uint8_t seq_nt16to4_table[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };

static uint8_t seq_nt6_table[256] = {
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 1, 5, 2,  5, 5, 5, 3,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  4, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,

    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,
    5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5,  5, 5, 5, 5
};

static void print_stat(errstat_t *e, int pos, int qthres)
{
	int i, j, k;
	int64_t sum[2] = {0, 0}, same[2] = {0, 0}, c[2][4][7];
	memset(c, 0, 2 * 4 * 7 * 8);
	for (i = 0; i < 94; ++i) {
		int x = i < qthres? 0 : 1;
		for (j = 0; j < 4; ++j) {
			for (k = 0; k < 7; ++k)
				if (k != 4)
					sum[x] += e->q[i][j][k], c[x][j][k] += e->q[i][j][k];
			same[x] += e->q[i][j][j];
		}
	}
	if (pos <= 0) printf("ALL");
	else printf("%d", pos);
	printf("\tQ%d", (int)(-4.343 * log((sum[0]+sum[1]-same[0]-same[1]+1e-6) / (sum[0]+sum[1]+1e-6))));
	for (i = 0; i < 2; ++i) {
		printf("\t%lld", (long long)sum[i]);
		for (j = 0; j < 4; ++j) {
			int64_t s = 0;
			putchar('\t');
			for (k = 0; k < 7; ++k)
				if (j != k && k != 4) s += c[i][j][k];
			for (k = 0; k < 4; ++k) {
				if (k) putchar(':');
				if (j == k) printf("Q%d", (int)(-4.343 * log((s+1e-6) / (c[i][j][j]+s+1e-6))));
				else printf("%.2d", (int)(100. * (c[i][j][k]+.25e-6) / (s+1e-6) + .499));
			}
			printf(":%.2d", (int)(100. * (c[i][j][5]+c[i][j][6]+.25e-6) / (s+1e-6) + .499));
		}
	}
	putchar('\n');
}

typedef struct {
	BGZF *fp;
	hts_itr_t *itr;
	bam_hdr_t *h;
	void *bed;
} aux_t;

static int read1_plp(void *data, bam1_t *b)
{
	int ret;
	aux_t *a = (aux_t*)data;
	ret = a->itr? bam_itr_next(a->fp, a->itr, b) : bam_read1(a->fp, b);
	if ((b->core.flag & (BAM_FSECONDARY|BAM_FSUPP)) || b->core.tid < 0)
		b->core.flag |= BAM_FUNMAP;
	else if (a->bed) {
		int tlen;
		const char *chr = a->h->target_name[b->core.tid];
		tlen = bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b));
		if (!bed_overlap(a->bed, chr, b->core.pos, b->core.pos + tlen))
			b->core.flag |= BAM_FUNMAP;
	}
	return ret;
}

int main_mapchk(int argc, char *argv[])
{
	int c, last_tid = -1, tid, pos, ref_len = 0, max_len = 0, max_alloc = 0, qthres = 20, min_sys_dp = 0, n_plp, excl_flag = 0;
	double fthres = 0.35;
	const bam_pileup1_t *plp;
	bam_mplp_t mplp;
	char *ref = 0, *reg = 0;
	faidx_t *fai;
	bam_hdr_t *h;
	aux_t aux = {0,0,0,0}, *auxp = &aux;
	errstat_t all, *e = 0;

	while ((c = getopt(argc, argv, "r:q:f:b:d:12")) >= 0) {
		if (c == 'r') reg = optarg;
		else if (c == 'q') qthres = atoi(optarg);
		else if (c == 'f') fthres = atof(optarg);
		else if (c == 'b') aux.bed = bed_read(optarg);
		else if (c == 'd') min_sys_dp = atoi(optarg);
		else if (c == '1') excl_flag |= 0x40;
		else if (c == '2') excl_flag |= 0x80;
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   htsbox mapchk [options] <aln.bam> <ref.fa>\n\n");
		fprintf(stderr, "Options: -r STR       region [null]\n");
		fprintf(stderr, "         -q INT       threshold for HIGH quality [%d]\n", qthres);
		fprintf(stderr, "         -f FLOAT     skip sites with excessive non-ref bases [%.2f]\n", fthres);
		fprintf(stderr, "         -b FILE      BED file to include []\n");
		fprintf(stderr, "         -d INT       min non-ref count [0]\n");
		fprintf(stderr, "         -1           exclude read1\n");
		fprintf(stderr, "         -2           exclude read2\n");
		fprintf(stderr, "\n");
		return 1;
	}
	aux.fp = bgzf_open(argv[optind], "r");
	fai = fai_load(argv[optind+1]);
	aux.h = h = bam_hdr_read(aux.fp);
	if (reg) {
		hts_idx_t *idx;
		idx = bam_index_load(argv[optind]);
		aux.itr = bam_itr_querys(idx, h, reg);
		hts_idx_destroy(idx);
	}
	mplp = bam_mplp_init(1, read1_plp, (void**)&auxp);
	while (bam_mplp_auto(mplp, &tid, &pos, &n_plp, &plp) > 0) {
		int i, r[2], n_var, n_high, cnt_high[8];
		if (n_plp == 0) continue;
		if (aux.bed && !bed_overlap(aux.bed, h->target_name[tid], pos, pos + 1)) continue;
		// get the reference sequence
		if (tid != last_tid) {
			free(ref);
			ref = faidx_fetch_seq(fai, h->target_name[tid], 0, INT_MAX, &ref_len);
			for (i = 0; i < ref_len; ++i)
				ref[i] = seq_nt6_table[(int)ref[i]] - 1;
			last_tid = tid;
		}
		r[0] = ref[pos], r[1] = 3 - r[0];
		if (r[0] >= 4) continue; // ignore the rest if the reference base is "N"
		// compute n_var
		memset(cnt_high, 0, 8 * sizeof(int));
		for (i = n_var = n_high = 0; i < n_plp; ++i) {
			bam_pileup1_t *p = (bam_pileup1_t*)&plp[i];
			const uint8_t *seq = bam_get_seq(p->b), *qual = seq + ((p->b->core.l_qseq + 1) >> 1);
			int b = seq_nt16to4_table[bam_seqi(seq, p->qpos)], q = qual[p->qpos];
			if (p->is_del) continue;
			p->aux = b<<8 | q; // cache the base and the quality
			if (q >= qthres) {
				++n_high;
				++cnt_high[b];
				if (p->indel > 0) ++cnt_high[5];
				else if (p->indel < 0) ++cnt_high[6];
				if (p->indel != 0 || b != r[0]) ++n_var;
			}
		}
		if (n_high >= 3 && (double)n_var / n_high > fthres) continue;
		// expand $e when necessary
		for (i = 0; i < n_plp; ++i)
			max_len = max_len > plp[i].b->core.l_qseq? max_len : plp[i].b->core.l_qseq;
		if (max_len > max_alloc) {
			int old_max = max_alloc;
			max_alloc = max_len;
			kroundup32(max_alloc);
			e = (errstat_t*)realloc(e, max_alloc * sizeof(errstat_t));
			memset(&e[old_max], 0, (max_alloc - old_max) * sizeof(errstat_t));
		}
		// fill $e
		for (i = 0; i < n_plp; ++i) {
			const bam_pileup1_t *p = &plp[i];
			int x = p->qpos, b = p->aux>>8, b0 = b, q = p->aux&0xff, is_rev = !!bam_is_rev(p->b);
			if (p->is_del) continue;
			if (p->b->core.flag & excl_flag) continue;
			if (is_rev) {
				x = p->b->core.l_qseq - 1 - p->qpos;
				b = b < 4? 3 - b : 4;
			}
			if (cnt_high[b0] >= min_sys_dp) ++e[x].q[q][b][r[is_rev]];
			if (p->indel > 0 && cnt_high[5] >= min_sys_dp) ++e[x].q[q][b][5];
			if (p->indel < 0 && cnt_high[6] >= min_sys_dp) ++e[x].q[q][b][6];
		}
	}
	bam_mplp_destroy(mplp);
	if (aux.itr) hts_itr_destroy(aux.itr);
	bam_hdr_destroy(h);
	free(ref);
	fai_destroy(fai);
	bgzf_close(aux.fp);
	if (aux.bed) bed_destroy(aux.bed);

	memset(&all, 0, sizeof(errstat_t));
	{
		int i, j, k, l; 
		for (l = 0; l < max_len; ++l)
			for (i = 0; i < 94; ++i)
				for (j = 0; j < 5; ++j)
					for (k = 0; k < 8; ++k)
						all.q[i][j][k] += e[l].q[i][j][k];
		print_stat(&all, 0, qthres);
		for (l = 0; l < max_len; ++l)
			print_stat(&e[l], l + 1, qthres);
	}
	free(e);
	return 0;
}
