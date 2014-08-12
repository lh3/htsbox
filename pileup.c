// This piece of code is modified from samtools/bam2depth.c
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <stdio.h>
#include "sam.h"
#include "faidx.h"
#include "ksort.h"

const char *hts_parse_reg(const char *s, int *beg, int *end);

typedef struct {     // auxiliary data structure
	BGZF *fp;        // the file handler
	hts_itr_t *itr;  // NULL if a region not specified
	int min_mapQ, min_len; // mapQ filter; length filter
} aux_t;

// This function reads a BAM alignment from one BAM file.
static int read_bam(void *data, bam1_t *b) // read level filters better go here to avoid pileup
{
	aux_t *aux = (aux_t*)data; // data in fact is a pointer to an auxiliary structure
	int ret = aux->itr? bam_itr_next(aux->fp, aux->itr, b) : bam_read1(aux->fp, b);
	if (!(b->core.flag&BAM_FUNMAP)) {
		if ((int)b->core.qual < aux->min_mapQ) b->core.flag |= BAM_FUNMAP;
		else if (aux->min_len > 0) {
			int k, l;
			const uint32_t *cigar = bam_get_cigar(b);
			for (k = l = 0; k < b->core.n_cigar; ++k) { // compute the query length in the alignment
				int op = bam_cigar_op(cigar[k]);
				if ((bam_cigar_type(op)&1) && op != BAM_CSOFT_CLIP)
					l += bam_cigar_oplen(cigar[k]);
			}
			if (l < aux->min_len) b->core.flag |= BAM_FUNMAP;
		}
	}
	return ret;
}

typedef struct {
	uint32_t is_skip:1, is_rev:1, b:4, q:8, k:18; // b=base, q=quality, k=allele id
	int indel; // <0: deleteion; >0: insertion
	uint64_t hash;
	uint64_t pos; // i<<32|j: j-th read of the i-th sample
} allele_t;

#define allele_lt(a, b) ((a).hash < (b).hash || ((a).hash == (b).hash && (a).indel < (b).indel))
KSORT_INIT(allele, allele_t, allele_lt)

static inline allele_t pileup2allele(const bam_pileup1_t *p, int min_baseQ, uint64_t pos, int ref)
{ // collect allele information given a pileup1 record
	allele_t a;
	int i;
	const uint8_t *seq = bam_get_seq(p->b);
	a.q = bam_get_qual(p->b)[p->qpos];
	a.is_rev = bam_is_rev(p->b);
	a.is_skip = (p->is_del || p->is_refskip || a.q < min_baseQ);
	a.indel = p->indel;
	a.b = a.hash = bam_seqi(seq, p->qpos);
	a.pos = pos;
	if (p->indel > 0) // compute the hash for the insertion
		for (i = 0; i < p->indel; ++i)
			a.hash = (a.hash<<4) + a.hash + bam_seqi(seq, p->qpos + i + 1);
	a.hash = a.hash << 1 >> 1;
	if (p->indel != 0 || a.b != ref || ref == 15) // the highest bit tells whether it is a reference allele or not
		a.hash |= 1ULL<<63;
	return a;
}

static inline void print_allele(const bam_pileup1_t *p, int l_ref, const char *ref, int pos, int max_del, int is_vcf)
{ // print the allele. The format depends on is_vcf.
	const uint8_t *seq = bam_get_seq(p->b);
	int i, rest = max_del;
	putchar(seq_nt16_str[bam_seqi(seq, p->qpos)]);
	if (p->indel > 0) {
		if (!is_vcf) printf("+%d", p->indel);
		for (i = 1; i <= p->indel; ++i)
			putchar(seq_nt16_str[bam_seqi(seq, p->qpos + i)]);
	} else if (p->indel < 0) {
		if (!is_vcf) {
			printf("%d", p->indel);
			for (i = 1; i <= -p->indel; ++i)
				putchar(pos + i < l_ref? toupper(ref[pos+i]) : 'N');
		} else rest -= -p->indel, pos += -p->indel;
	}
	if (is_vcf)
		for (i = 1; i <= rest; ++i)
			putchar(pos + i < l_ref? toupper(ref[pos+i]) : 'N');
}

typedef struct {
	int n_a, n_alleles, max_del; // n_a: #reads used to compute quality sum; max_del: max deletion length
	int tot_dp, max_dp, n_cnt, max_cnt;
	allele_t *a; // allele of each read, of size tot_dp
	int *cnt_b, *cnt_q; // cnt_b: count of supporting reads on both strands; cnt_q: sum of quality
	int *sum_q; // sum of qual across entire _a_. It points to the last "row" of cnt_q.
} paux_t;

static void count_alleles(paux_t *pa, int n)
{
	allele_t *a = pa->a;
	int i, j;
	a[0].k = 0; // the first allele is given allele id 0
	for (i = pa->n_alleles = 1, pa->max_del = 0; i < pa->n_a; ++i) {
		if (a[i].indel != a[i-1].indel || a[i].hash != a[i-1].hash) // change of allele
			++pa->n_alleles;
		a[i].k = pa->n_alleles - 1;
		pa->max_del = pa->max_del > -a[i].indel? pa->max_del : -a[i].indel; // max deletion
	}
	// collect per-BAM counts
	pa->n_cnt = pa->n_alleles * (n + 1);
	if (pa->n_cnt > pa->max_cnt) { // expand the arrays if necessary
		pa->max_cnt = pa->n_cnt;
		kroundup32(pa->max_cnt);
		pa->cnt_b = (int*)realloc(pa->cnt_b, pa->max_cnt * 2 * sizeof(int));
		pa->cnt_q = (int*)realloc(pa->cnt_q, pa->max_cnt * sizeof(int));
	}
	memset(pa->cnt_b, 0, pa->n_cnt * 2 * sizeof(int));
	memset(pa->cnt_q, 0, pa->n_cnt * sizeof(int));
	pa->sum_q = pa->cnt_q + pa->n_alleles * n; // points to the last row of cnt_q
	for (i = 0; i < pa->n_a; ++i) { // compute counts and sums of qualities
		j = (a[i].pos>>32)*pa->n_alleles + a[i].k;
		++pa->cnt_b[j<<1|a[i].is_rev];
		pa->cnt_q[j] += a[i].q;
		pa->sum_q[a[i].k] += a[i].q;
	}
}

int main_pileup(int argc, char *argv[])
{
	int i, j, n, tid, beg, end, pos, last_tid, *n_plp, baseQ = 0, mapQ = 0, min_len = 0, l_ref = 0, depth_only = 0, min_sum_q = 0, is_vcf = 0, var_only = 0;
	const bam_pileup1_t **plp;
	char *ref = 0, *reg = 0, *chr_end; // specified region
	faidx_t *fai = 0;
	bam_hdr_t *h = 0; // BAM header of the 1st input
	aux_t **data;
	paux_t aux;
	bam_mplp_t mplp;

	// parse the command line
	while ((n = getopt(argc, argv, "r:q:Q:l:f:dvcs:")) >= 0) {
		if (n == 'f') fai = fai_load(optarg);
		else if (n == 'l') min_len = atoi(optarg); // minimum query length
		else if (n == 'r') reg = strdup(optarg);   // parsing a region requires a BAM header
		else if (n == 'Q') baseQ = atoi(optarg);   // base quality threshold
		else if (n == 'q') mapQ = atoi(optarg);    // mapping quality threshold
		else if (n == 's') min_sum_q = atoi(optarg);
		else if (n == 'd') depth_only = 1;
		else if (n == 'v') var_only = 1;
		else if (n == 'c') is_vcf = var_only = 1;
	}
	if (optind == argc) {
        fprintf(stderr, "\n");
        fprintf(stderr, "Usage:   pileup [options] in1.bam [in2.bam [...]]\n\n");
        fprintf(stderr, "Options: -f FILE    reference genome [null]\n");
		fprintf(stderr, "         -l INT     minimum query length [%d]\n", min_len);
		fprintf(stderr, "         -q INT     minimum mapping quality [%d]\n", mapQ);
		fprintf(stderr, "         -Q INT     minimum base quality [%d]\n", baseQ);
		fprintf(stderr, "         -s INT     drop alleles with sum of quality below INT [%d]\n", min_sum_q);
		fprintf(stderr, "         -r STR     region [null]\n");
		fprintf(stderr, "         -v         show variants only\n");
		fprintf(stderr, "         -c         output in the VCF format (force -v)\n");
		fprintf(stderr, "\n");
		return 1;
	}

	// initialize the auxiliary data structures
	n = argc - optind; // the number of BAMs on the command line
	data = (aux_t**)calloc(n, sizeof(aux_t*)); // data[i] for the i-th input
	beg = 0; end = 1<<30; tid = -1;  // set the default region
	if (reg) {
		chr_end = (char*)hts_parse_reg(reg, &beg, &end);
		ref = fai? fai_fetch(fai, reg, &l_ref) : 0;
	} else chr_end = 0;

	// load the index or put the file position at the right place
	last_tid = -1;
	for (i = 0; i < n; ++i) {
		bam_hdr_t *htmp;
		data[i] = (aux_t*)calloc(1, sizeof(aux_t));
		data[i]->fp = bgzf_open(argv[optind+i], "r"); // open BAM
		data[i]->min_mapQ = mapQ;                     // set the mapQ filter
		data[i]->min_len  = min_len;                  // set the qlen filter
		htmp = bam_hdr_read(data[i]->fp);             // read the BAM header
		if (i == 0 && chr_end) {
			char c = *chr_end;
			*chr_end = 0;
			last_tid = tid = bam_name2id(htmp, reg);
			*chr_end = c;
		}
		if (i) bam_hdr_destroy(htmp); // if not the 1st BAM, trash the header
		else h = htmp; // keep the header of the 1st BAM
		if (tid >= 0) { // if a region is specified and parsed successfully
			hts_idx_t *idx = bam_index_load(argv[optind+i]); // load the index
			data[i]->itr = bam_itr_queryi(idx, tid, beg, end); // set the iterator
			hts_idx_destroy(idx); // the index is not needed any more; phase out of the memory
		}
	}

	// the core multi-pileup loop
	mplp = bam_mplp_init(n, read_bam, (void**)data); // initialization
	n_plp = (int*)calloc(n, sizeof(int)); // n_plp[i] is the number of covering reads from the i-th BAM
	plp = (const bam_pileup1_t**)calloc(n, sizeof(const bam_pileup1_t*)); // plp[i] points to the array of covering reads (internal in mplp)
	memset(&aux, 0, sizeof(paux_t));
	if (is_vcf) {
		puts("##fileformat=VCFv4.1");
		puts("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
		puts("##FORMAT=<ID=SQ,Number=A,Type=Integer,Description=\"Sum of quality for each allele\">");
		puts("##FORMAT=<ID=FC,Number=A,Type=Integer,Description=\"Number of supporting reads on the forward strand\">");
		puts("##FORMAT=<ID=RC,Number=A,Type=Integer,Description=\"Number of supporting reads on the reverse strand\">");
	}
	while (bam_mplp_auto(mplp, &tid, &pos, n_plp, plp) > 0) { // come to the next covered position
		if (pos < beg || pos >= end) continue; // out of range; skip
		for (i = aux.tot_dp = 0; i < n; ++i) aux.tot_dp += n_plp[i];
		if (last_tid != tid && fai) { // switch of chromosomes
			free(ref);
			ref = fai_fetch(fai, h->target_name[tid], &l_ref);
			last_tid = tid; // note that last_tid is only used when fai is not NULL
		}
		if (aux.tot_dp == 0) continue; // well, this should not happen
		if (depth_only) { // only print read depth; no allele information
			fputs(h->target_name[tid], stdout); printf("\t%d", pos+1); // a customized printf() would be faster
			if (ref == 0 || pos >= l_ref + beg) printf("\tN");
			else printf("\t%c", ref[pos - beg]);
			for (i = 0; i < n; ++i) { // base level filters have to go here
				int m = 0;
				for (j = 0; j < n_plp[i]; ++j) {
					const bam_pileup1_t *p = plp[i] + j; // DON'T modfity plp[][] unless you really know
					if (p->is_del || p->is_refskip) ++m; // having dels or refskips at tid:pos
					else if (bam_get_qual(p->b)[p->qpos] < baseQ) ++m; // low base quality
				}
				printf("\t%d", n_plp[i] - m); // this the depth to output
			}
		} else { // detailed summary of each allele
			int k, r = 15, a1, a2, shift = 0, qual;
			allele_t *a;
			if (aux.tot_dp + 1 > aux.max_dp) { // expand array
				aux.max_dp = aux.tot_dp + 1;
				kroundup32(aux.max_dp);
				aux.a = (allele_t*)realloc(aux.a, aux.max_dp * sizeof(allele_t));
			}
			a = aux.a;
			// collect alleles
			r = (ref && pos - beg < l_ref)? seq_nt16_table[(int)ref[pos - beg]] : 15; // the reference allele
			for (i = aux.n_a = 0; i < n; ++i) {
				for (j = 0; j < n_plp[i]; ++j) {
					allele_t aa;
					aa = pileup2allele(&plp[i][j], baseQ, (uint64_t)i<<32 | j, r);
					if (!aa.is_skip) a[aux.n_a++] = aa;
				}
			}
			if (aux.n_a == 0) continue; // no reads are good enough; zero effective coverage
			// count alleles
			ks_introsort(allele, aux.n_a, aux.a);
			count_alleles(&aux, n);
			// squeeze out weak alleles
			if (min_sum_q > 0) {
				for (i = k = 0; i < aux.n_a; ++i)
					if (aux.sum_q[a[i].k] >= min_sum_q)
						a[k++] = a[i];
				if (k < aux.n_a) {
					aux.n_a = k;
					count_alleles(&aux, n);
				}
			}
			if (var_only && aux.n_alleles == 1 && a[0].hash>>63 == 0) continue; // var_only mode, but no ALT allele; skip
			// identify alleles
			if (aux.n_alleles >= 2) {
				int max1 = -1, max2 = -1;
				for (i = 0; i < aux.n_alleles; ++i)
					if (aux.sum_q[i] > max1) max2 = max1, a2 = a1, max1 = aux.sum_q[i], a1 = i;
					else if (aux.sum_q[i] > max2) max2 = aux.sum_q[i], a2 = i;
				qual = (a1 == 0 && a[0].hash>>63 == 0)? max2 : max1;
			} else a1 = a2 = 0, qual = aux.sum_q[0];
			// print
			fputs(h->target_name[tid], stdout); printf("\t%d", pos+1);
			if (is_vcf) {
				fputs("\t.\t", stdout);
				for (i = 0; i <= aux.max_del; ++i) // print the reference allele up to the longest deletion
					putchar(ref && pos + i < l_ref + beg? ref[pos + i - beg] : 'N');
				putchar('\t');
			} else printf("\t%c\t", ref && pos < l_ref + beg? ref[pos - beg] : 'N'); // print a single reference base
			// print alleles
			if (!is_vcf || a[0].hash>>63) { // print if there is no reference allele
				print_allele(&plp[a[0].pos>>32][(uint32_t)a[0].pos], l_ref, ref, pos - beg, aux.max_del, is_vcf);
				if (aux.n_alleles > 1) putchar(',');
			}
			for (i = k = 1; i < aux.n_a; ++i)
				if (a[i].indel != a[i-1].indel || a[i].hash != a[i-1].hash) {
					print_allele(&plp[a[i].pos>>32][(uint32_t)a[i].pos], l_ref, ref, pos - beg, aux.max_del, is_vcf);
					if (++k != aux.n_alleles) putchar(',');
				}
			if (is_vcf && aux.n_alleles == 1 && a[0].hash>>63 == 0) putchar('.'); // print placeholder if there is only the reference allele
			if (is_vcf) printf("\t%d\t.\t.\tGT:SQ:FC:RC", qual);
			// print counts
			shift = (is_vcf && a[0].hash>>63); // in VCF, if there is no ref allele, we need to shift the allele number
			for (i = k = 0; i < n; ++i, k += aux.n_alleles) {
				printf("\t%d/%d:", a1 + shift, a2 + shift);
				if (shift) fputs("0,", stdout);
				for (j = 0; j < aux.n_alleles; ++j) {
					if (j) putchar(',');
					printf("%d", aux.cnt_q[k+j]);
				}
				putchar(':');
				if (shift) fputs("0,", stdout);
				for (j = 0; j < aux.n_alleles; ++j) {
					if (j) putchar(',');
					printf("%d", aux.cnt_b[(k+j)<<1]);
				}
				putchar(':');
				if (shift) fputs("0,", stdout);
				for (j = 0; j < aux.n_alleles; ++j) {
					if (j) putchar(',');
					printf("%d", aux.cnt_b[(k+j)<<1|1]);
				}
			}
		}
		putchar('\n');
	}
	free(n_plp); free(plp);
	bam_mplp_destroy(mplp);

	bam_hdr_destroy(h);
	for (i = 0; i < n; ++i) {
		bgzf_close(data[i]->fp);
		if (data[i]->itr) bam_itr_destroy(data[i]->itr);
		free(data[i]);
	}
	if (ref) free(ref);
	if (fai) fai_destroy(fai);
	free(aux.cnt_b); free(aux.cnt_q); free(aux.a);
	free(data); free(reg);
	return 0;
}
