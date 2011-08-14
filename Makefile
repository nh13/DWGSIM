PACKAGE_VERSION="0.0.3"
CC=			gcc
CFLAGS=		-g -Wall -O2 #-m64 #-arch ppc
DFLAGS=		-D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -DPACKAGE_VERSION=\\\"${PACKAGE_VERSION}\\\"
PACKAGE_VERSION=0.0.3
DWGSIM_AOBJS = dwgsim.o
DWGSIM_EVAL_AOBJS = dwgsim_eval.o \
					samtools/knetfile.o \
					samtools/bgzf.o samtools/kstring.o samtools/bam_aux.o samtools/bam.o samtools/bam_import.o samtools/sam.o samtools/bam_index.o \
					samtools/bam_pileup.o samtools/bam_lpileup.o samtools/bam_md.o samtools/razf.o samtools/faidx.o samtools/bedidx.o \
					samtools/bam_sort.o samtools/sam_header.o samtools/bam_reheader.o samtools/kprobaln.o samtools/bam_cat.o

PROG=		dwgsim dwgsim_eval
INCLUDES=	-I.
SUBDIRS=	samtools . 
LIBPATH=

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all-recur lib-recur clean-recur cleanlocal-recur install-recur:
		@target=`echo $@ | sed s/-recur//`; \
		wdir=`pwd`; \
		list='$(SUBDIRS)'; for subdir in $$list; do \
			cd $$subdir; \
			$(MAKE) CC="$(CC)" DFLAGS="$(DFLAGS)" CFLAGS="$(CFLAGS)" \
				INCLUDES="$(INCLUDES)" LIBPATH="$(LIBPATH)" $$target || exit 1; \
			cd $$wdir; \
		done;

all:$(PROG)

.PHONY:all lib clean cleanlocal
.PHONY:all-recur lib-recur clean-recur cleanlocal-recur install-recur

dwgsim:lib-recur $(DWGSIM_AOBJS)
	$(CC) $(CFLAGS) -o $@ $(DWGSIM_AOBJS) -lm -lz

dwgsim_eval:lib-recur $(DWGSIM_EVAL_AOBJS)
	$(CC) $(CFLAGS) -o $@ $(DWGSIM_EVAL_AOBJS) -Lsamtools -lm -lz

cleanlocal:
		rm -fr gmon.out *.o a.out *.exe *.dSYM razip bgzip $(PROG) *~ *.a *.so.* *.so *.dylib

clean:cleanlocal-recur
