#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "sam.h"

typedef struct {
	uint64_t sum_all, sum_high, sum_all2;
	int start, cnt;
} counter_t;

static counter_t null_count = { 0, 0, 0, 0, 0 };

static inline void print_counts(counter_t *cnt, const char *ctg)
{
	if (cnt->cnt)
		printf("%s\t%d\t%d\t%.2f\t%.2f\n", ctg, cnt->start, cnt->cnt,
			   (double)cnt->sum_all / cnt->cnt, (double)cnt->sum_high / cnt->cnt);
	*cnt = null_count;
}

int main_depth(int argc, char *argv[])
{
	int c, qthres = 20;
	int n, tid, pos, upper = -1, lower = -1, last_tid = -1, last_pos = -2, i;
	double pc = 0.2;
	counter_t cnt;
	BGZF *fp;
	bam_hdr_t *h;
	bam_plp_t plp;
	const bam_pileup1_t *p;
	while ((c = getopt(argc, argv, "q:p:")) >= 0) {
		if (c == 'q') qthres = atoi(optarg);
		else if (c == 'p') pc = atof(optarg);
	}
	if (optind == argc) {
		fprintf(stderr, "Usage: htsbox depth [-q %d] [-p %.2f] <in.bam>\n", qthres, pc);
		return 1;
	}
	fp = bgzf_open(argv[optind], "r");
	h = bam_hdr_read(fp);
	plp = bam_plp_init((bam_plp_auto_f)bam_read1, fp);
	bam_plp_set_maxcnt(plp, 1000000);
	cnt = null_count;
	while ((p = bam_plp_auto(plp, &tid, &pos, &n)) != 0) {
		if (tid != last_tid || pos != last_pos + 1 || n < lower || n > upper) {
			if (last_tid >= 0)
				print_counts(&cnt, h->target_name[last_tid]);
			cnt.start = pos;
			lower = (int)(n * (1. - pc) + .499);
			upper = (int)(n * (1. + pc) + .499);
		}
		for (i = 0; i < n; ++i)
			if (p[i].b->core.qual >= qthres)
				++cnt.sum_high;
		cnt.sum_all += n, cnt.sum_all2 += n * n, ++cnt.cnt;
		last_tid = tid, last_pos = pos;
	}
	if (last_tid >= 0)
		print_counts(&cnt, h->target_name[last_tid]);
	bam_plp_destroy(plp);
	bam_hdr_destroy(h);
	bgzf_close(fp);
	return 0;
}
