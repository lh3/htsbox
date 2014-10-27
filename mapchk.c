#include <limits.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "sam.h"
#include "faidx.h"

typedef struct {
	int64_t q[94][5][8]; // qual, read base, ref base (5=ins, 6=del, 7=clip)
} errstat_t;

static uint8_t seq_nt16comp_table[] = { 0, 8, 4, 12, 2, 10, 9, 14, 1, 6, 5, 13, 3, 11, 7, 15 };
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
	int64_t sum[2], c[2][4][7];
	sum[0] = sum[1] = 0;
	memset(c, 0, 2 * 4 * 7 * 8);
	for (i = 0; i < 94; ++i) {
		int x = i < qthres? 0 : 1;
		for (j = 0; j < 4; ++j)
			for (k = 0; k < 7; ++k)
				if (k != 4)
					sum[x] += e->q[i][j][k], c[x][j][k] += e->q[i][j][k];
	}
	if (pos <= 0) printf("ALL");
	else printf("%d", pos);
	for (i = 0; i < 2; ++i) {
		printf("\t%lld", (long long)sum[i]);
		for (j = 0; j < 4; ++j) {
			int64_t s = 0;
			putchar('\t');
			for (k = 0; k < 7; ++k)
				if (j != k && k != 4) s += c[i][j][k];
			for (k = 0; k < 7; ++k) {
				if (k != 4) {
					if (k) putchar(':');
					if (j == k) printf("Q%d", (int)(-4.343 * log((s+1e-6) / (c[i][j][j]+s+1e-6))));
					else printf("%.2d", (int)(100. * (c[i][j][k]+.25e-6) / (s+1e-6) + .499));
				}
			}
		}
	}
	putchar('\n');
}

static inline int read1(BGZF *fp, bam1_t *b, hts_itr_t *itr)
{
	if (itr == 0) return bam_read1(fp, b);
	else return bam_itr_next(fp, itr, b);
}

int main_mapchk2(int argc, char *argv[])
{
	int c, tid = -1, ref_len = 0, max_len = 0, max_alloc = 0, qthres = 20;
	char *ref = 0, *reg = 0;
	faidx_t *fai;
	bam_hdr_t *h;
	hts_itr_t *itr = 0;
	bam1_t *b;
	BGZF *fp;
	errstat_t all, *e = 0;

	while ((c = getopt(argc, argv, "r:q:")) >= 0) {
		if (c == 'r') reg = optarg;
		else if (c == 'q') qthres = atoi(optarg);
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   htsbox mapchk [options] <aln.bam> <ref.fa>\n\n");
		fprintf(stderr, "Options: -r STR       region [null]\n");
		fprintf(stderr, "         -q INT       threshold for HIGH quality [%d]\n", qthres);
		fprintf(stderr, "\n");
		return 1;
	}
	fp = bgzf_open(argv[optind], "r");
	fai = fai_load(argv[optind+1]);
	h = bam_hdr_read(fp);
	if (reg) {
		hts_idx_t *idx;
		idx = bam_index_load(argv[optind]);
		itr = bam_itr_querys(idx, h, reg);
		hts_idx_destroy(idx);
	}
	b = bam_init1();
	while (read1(fp, b, itr) >= 0) {
		const uint32_t *cigar;
		int k, qlen, is_rev, x, y;
		const uint8_t *seq, *qual;
		if ((b->core.flag & (BAM_FUNMAP|BAM_FSECONDARY|BAM_FSUPP)) || b->core.tid < 0)
			continue;
		// get the reference sequence
		if (b->core.tid != tid) {
			free(ref);
			ref = faidx_fetch_seq(fai, h->target_name[b->core.tid], 0, INT_MAX, &ref_len);
			for (k = 0; k < ref_len; ++k)
				ref[k] = seq_nt6_table[(int)ref[k]] - 1;
			tid = b->core.tid;
		}
		// get the entire query length (not the length in the alignment)
		cigar = bam_get_cigar(b);
		for (k = qlen = 0; k < b->core.n_cigar; ++k) {
			int op = bam_cigar_op(cigar[k]), len = bam_cigar_oplen(cigar[k]);
			if (op == BAM_CMATCH || op == BAM_CINS || op == BAM_CSOFT_CLIP)
				qlen += len;
			else if (op == BAM_CHARD_CLIP) break;
		}
		if (k < b->core.n_cigar) continue; // the primary record should not have hard clipping
		// expand $e when necessary
		max_len = max_len > qlen? max_len : qlen;
		if (max_len > max_alloc) {
			int old_max = max_alloc;
			max_alloc = max_len;
			kroundup32(max_alloc);
			e = (errstat_t*)realloc(e, max_alloc * sizeof(errstat_t));
			memset(&e[old_max], 0, (max_alloc - old_max) * sizeof(errstat_t));
		}
		// count errors
		is_rev = bam_is_rev(b);
		seq = bam_get_seq(b);
		qual = bam_get_qual(b);
		for (k = x = 0, y = b->core.pos; k < b->core.n_cigar; ++k) {
			int i, op = bam_cigar_op(cigar[k]), len = bam_cigar_oplen(cigar[k]);
			if (op == BAM_CMATCH) { // FIXME: check BAM_CEQUAL and BAM_CDIFF
				for (i = 0; i < len; ++i) { // NB: "i<len && i+x<qlen && i+y<ref_len" would be safer
					int c16 = bam_seqi(seq, i + x), b = ref[i+y];
					if (is_rev) ++e[qlen-1-i-x].q[qual[i+x]][seq_nt16to4_table[seq_nt16comp_table[c16]]][b<4? 3-b : 4];
					else ++e[i+x].q[qual[i+x]][seq_nt16to4_table[c16]][b];
				}
				x += len, y += len;
			} else if (op == BAM_CSOFT_CLIP || op == BAM_CINS) {
				int r = op == BAM_CSOFT_CLIP? 7 : 5;
				for (i = 0; i < len; ++i) {
					int c16 = bam_seqi(seq, i + x);
					if (is_rev) ++e[qlen-1-i-x].q[qual[i+x]][seq_nt16to4_table[seq_nt16comp_table[c16]]][r];
					else ++e[i+x].q[qual[i+x]][seq_nt16to4_table[c16]][r];
				}
				x += len;
			} else if (op == BAM_CDEL) {
				if (x > 0) { // ignore deletion at the beginning of the alignment
					int c16 = bam_seqi(seq, x - 1);
					if (is_rev) ++e[qlen-x].q[qual[x-1]][seq_nt16to4_table[seq_nt16comp_table[c16]]][6];
					else ++e[x-1].q[qual[x-1]][seq_nt16to4_table[c16]][6];
				}
				y += len;
			}
		}
	}
	bam_destroy1(b);
	if (itr) hts_itr_destroy(itr);
	bam_hdr_destroy(h);
	free(ref);
	fai_destroy(fai);
	bgzf_close(fp);

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

typedef struct {
	BGZF *fp;
	hts_itr_t *itr;
} aux_t;

static int read1_plp(void *data, bam1_t *b)
{
	int ret;
	aux_t *a = (aux_t*)data;
	ret = a->itr? bam_itr_next(a->fp, a->itr, b) : bam_read1(a->fp, b);
	if ((b->core.flag & (BAM_FSECONDARY|BAM_FSUPP)) || b->core.tid < 0)
		b->core.flag |= BAM_FUNMAP;
	return ret;
}

int main_mapchk(int argc, char *argv[])
{
	int c, last_tid = -1, tid, pos, ref_len = 0, max_len = 0, max_alloc = 0, qthres = 20, n_plp;
	double fthres = 0.4;
	const bam_pileup1_t *plp;
	bam_mplp_t mplp;
	char *ref = 0, *reg = 0;
	faidx_t *fai;
	bam_hdr_t *h;
	aux_t aux = {0,0}, *auxp = &aux;
	errstat_t all, *e = 0;

	while ((c = getopt(argc, argv, "r:q:f:")) >= 0) {
		if (c == 'r') reg = optarg;
		else if (c == 'q') qthres = atoi(optarg);
		else if (c == 'f') fthres = atof(optarg);
	}
	if (optind + 2 > argc) {
		fprintf(stderr, "\n");
		fprintf(stderr, "Usage:   htsbox mapchk [options] <aln.bam> <ref.fa>\n\n");
		fprintf(stderr, "Options: -r STR       region [null]\n");
		fprintf(stderr, "         -q INT       threshold for HIGH quality [%d]\n", qthres);
		fprintf(stderr, "         -f FLOAT     skip sites with excessive non-ref bases [%.2f]\n", fthres);
		fprintf(stderr, "\n");
		return 1;
	}
	aux.fp = bgzf_open(argv[optind], "r");
	fai = fai_load(argv[optind+1]);
	h = bam_hdr_read(aux.fp);
	if (reg) {
		hts_idx_t *idx;
		idx = bam_index_load(argv[optind]);
		aux.itr = bam_itr_querys(idx, h, reg);
		hts_idx_destroy(idx);
	}
	mplp = bam_mplp_init(1, read1_plp, (void**)&auxp);
	while (bam_mplp_auto(mplp, &tid, &pos, &n_plp, &plp) > 0) {
		int i, r[2], n_var, n_high;
		if (n_plp == 0) continue;
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
		for (i = n_var = n_high = 0; i < n_plp; ++i) {
			bam_pileup1_t *p = (bam_pileup1_t*)&plp[i];
			const uint8_t *seq = bam_get_seq(p->b), *qual = seq + ((p->b->core.l_qseq + 1) >> 1);
			int b = seq_nt16to4_table[bam_seqi(seq, p->qpos)], q = qual[p->qpos];
			p->aux = b<<8 | q; // cache the base and the quality
			if (q >= qthres) ++n_high;
			if (q >= qthres && (p->indel != 0 || b != r[0])) ++n_var;
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
			int x = p->qpos, b = p->aux>>8, q = p->aux&0xff, is_rev = !!bam_is_rev(p->b);
			if (is_rev) {
				x = p->b->core.l_qseq - 1 - p->qpos;
				b = b < 4? 3 - b : 4;
			}
			++e[x].q[q][b][r[is_rev]];
			if (p->indel > 0) ++e[x].q[q][b][5];
			if (p->indel < 0) ++e[x].q[q][b][6];
		}
	}
	bam_mplp_destroy(mplp);
	if (aux.itr) hts_itr_destroy(aux.itr);
	bam_hdr_destroy(h);
	free(ref);
	fai_destroy(fai);
	bgzf_close(aux.fp);

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
