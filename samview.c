#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "sam.h"
#include "kstring.h"

void *bed_read(const char *fn);
int bed_overlap(const void *_h, const char *chr, int beg, int end);

static void print_pas(const bam_hdr_t *h, const bam1_t *b, kstring_t *buf)
{
	buf->l = 0;
	kputsn(bam_get_qname(b), b->core.l_qname - 1, buf);
	if (b->core.flag>>6&3) {
		kputc('/', buf);
		kputc("0123"[b->core.flag>>6&3], buf);
	}
	kputc('\t', buf);
	if ((b->core.flag & BAM_FUNMAP) || b->core.tid < 0 || b->core.n_cigar == 0) {
		kputw(b->core.l_qseq, buf);
		kputs("\t0\t0\t*\t*\t0\t0\t0\t1.0000\t0", buf);
	} else {
		const uint32_t *cigar = bam_get_cigar(b);
		const uint8_t *tag;
		int clip[2], k, m = 0, nm = 0, od = 0, oi = 0, os = 0, ns = 0, nd = 0, ni = 0, ql = 0, tl = 0, is_rev = !!bam_is_rev(b);
		double diff;
		clip[0] = clip[1] = 0;
		for (k = 0; k < b->core.n_cigar; ++k) {
			int op = bam_cigar_op(cigar[k]);
			int l = bam_cigar_oplen(cigar[k]);
			if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) m += l;
			else if (op == BAM_CINS) ++oi, ni += l;
			else if (op == BAM_CDEL) ++od, oi += l;
			else if (op == BAM_CREF_SKIP) ++os, ns += l;
			else if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) clip[!!k] = l;
		}
		tag = bam_aux_get(b, "NM");
		nm = tag? bam_aux2i(tag) : 0;
		if (nm < ni + oi) nm = ni + oi;

		ql = clip[0] + m + ni + clip[1];
		tl = h->target_len[b->core.tid];

		kputw(ql, buf); kputc('\t', buf);
		kputw(clip[is_rev], buf); kputc('\t', buf);
		kputw(ql - clip[!is_rev], buf); kputc('\t', buf);
		kputc("+-"[is_rev], buf); kputc('\t', buf);
		kputs(h->target_name[b->core.tid], buf); kputc('\t', buf);
		kputw(tl, buf); kputc('\t', buf);
		kputw(b->core.pos, buf); kputc('\t', buf);
		kputw(b->core.pos + m + nd + ns, buf); kputc('\t', buf);
		diff = (double)nm / (m + ni + nd);
		ksprintf(buf, "%.4f\t%d\tmm:%d;ins:%d,%d;del:%d,%d", diff, b->core.qual, nm - ni - nd, oi, ni, od, nd);
		tag = bam_aux_get(b, "AS");
		if (tag) ksprintf(buf, ";score:%d", bam_aux2i(tag));
	}
	puts(buf->s);
}

static void print_line(int flag, htsFile *out, kstring_t *buf, const bam_hdr_t *h, bam1_t *b, const void *bed)
{
	int i;
	if (!(flag&4)) {
		if (bed) {
			int n[16], k;
			const uint32_t *cigar;
			const bam1_core_t *c = &b->core;
			if (b->core.tid < 0) return;
			memset(n, 0, 16 * sizeof(int));
			cigar = bam_get_cigar(b);
			for (k = 0; k < c->n_cigar; ++k)
				n[bam_cigar_op(cigar[k])] += bam_cigar_oplen(cigar[k]);
			if (!bed_overlap(bed, h->target_name[c->tid], c->pos, c->pos + n[0] + n[2] + n[3] + n[7] + n[8]))
				return;
		}
		if (flag&8) {
			uint8_t *qual = bam_get_qual(b);
			for (i = 0; i < b->core.l_qseq; ++i)
				qual[i] = qual[i] >= 20? 30 : 10;
		}
		if (flag&16) {
			uint8_t *qual = bam_get_qual(b);
			for (i = 0; i < b->core.l_qseq; ++i) qual[i] = 25;
		}
		sam_write1(out, h, b);
	} else print_pas(h, b, buf);
}

int main_samview(int argc, char *argv[])
{
	samFile *in;
	char *fn_ref = 0;
	int flag = 0, c, clevel = -1, ignore_sam_err = 0;
	char moder[8];
	void *bed = 0;
	kstring_t buf = {0, 0, 0};
	bam_hdr_t *h;
	bam1_t *b;

	while ((c = getopt(argc, argv, "IbpSl:t:QUL:")) >= 0) {
		switch (c) {
		case 'S': flag |= 1; break;
		case 'b': flag |= 2; break;
		case 'p': flag |= 4; break;
		case 'Q': flag |= 8; break; // 1-bit quality (<Q20 to Q10; >=Q20 to Q30)
		case 'U': flag |= 16; break; // no quality (to Q25)
		case 'l': clevel = atoi(optarg); flag |= 2; break;
		case 't': fn_ref = optarg; break;
		case 'I': ignore_sam_err = 1; break;
		case 'L': bed = bed_read(optarg); break;
		}
	}
	if (argc == optind) {
		fprintf(stderr, "Usage: samview [-bSIp] [-L reg.bed] [-l level] <in.bam>|<in.sam> [region]\n");
		return 1;
	}
	strcpy(moder, "r");
	if ((flag&1) == 0) strcat(moder, "b");

	in = sam_open(argv[optind], moder, fn_ref);
	h = sam_hdr_read(in);
	h->ignore_sam_err = ignore_sam_err;
	b = bam_init1();

	{ // SAM/BAM output
		htsFile *out = 0;
		char modew[8];
		strcpy(modew, "w");
		if (clevel >= 0 && clevel <= 9) sprintf(modew + 1, "%d", clevel);
		if (flag&2) strcat(modew, "b");
		if (!(flag&4)) {
			out = hts_open("-", modew, 0);
			sam_hdr_write(out, h);
		}
		if (optind + 1 < argc && !(flag&1)) { // BAM input and has a region
			int i;
			hts_idx_t *idx;
			if ((idx = bam_index_load(argv[optind])) == 0) {
				fprintf(stderr, "[E::%s] fail to load the BAM index\n", __func__);
				return 1;
			}
			for (i = optind + 1; i < argc; ++i) {
				hts_itr_t *iter;
				if ((iter = bam_itr_querys(idx, h, argv[i])) == 0) {
					fprintf(stderr, "[E::%s] fail to parse region '%s'\n", __func__, argv[i]);
					continue;
				}
				while (bam_itr_next((BGZF*)in->fp, iter, b) >= 0)
					print_line(flag, out, &buf, h, b, bed);
				hts_itr_destroy(iter);
			}
			hts_idx_destroy(idx);
		} else {
			while (sam_read1(in, h, b) >= 0)
				print_line(flag, out, &buf, h, b, bed);
		}
		if (out) sam_close(out);
	}

	free(buf.s);
	bam_destroy1(b);
	bam_hdr_destroy(h);
	sam_close(in);
	return 0;
}
