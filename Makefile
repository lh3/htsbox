CC=			gcc
CFLAGS=		-g -Wall -Wc++-compat -Wno-unused-function -O2
CPPFLAGS=	-Ihtslib
OBJS=		main.o samview.o vcfview.o bamidx.o bcfidx.o bamshuf.o bam2fq.o tabix.o \
			abreak.o bam2bed.o razf.o razip.o faidx.o bedidx.o pileup.o mapchk.o depth.o genreg.o \
			kthread.o qualbin.o
PROG=		htsbox

.SUFFIXES:.c .o
.PHONY:all lib

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@

all:$(PROG)

lib:
		cd htslib; $(MAKE) CC="$(CC)" CFLAGS="$(CFLAGS)" libhts.a || exit 1; cd ..

qualbin.o:qualbin.c
		$(CC) -c $(CFLAGS) $(CPPFLAGS) -DBGZF_MT $< -o $@

htsbox:lib $(OBJS)
		$(CC) $(CFLAGS) -o $@ $(OBJS) -Lhtslib -lhts -lpthread -lz -lm

clean:
		rm -fr gmon.out *.o a.out *.dSYM *~ $(PROG); cd htslib; $(MAKE) clean; cd ..
