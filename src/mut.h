#ifndef MUT_H
#define MUT_H

#include "contigs.h"
#include "mut_txt.h"
#include "mut_bed.h"
#include "dwgsim_opt.h"
#include "dwgsim.h"

/* FASTA parser, copied from seq.c */
typedef struct {
    int l, m; /* length and maximum buffer size */
    unsigned char *s; /* sequence */
} seq_t;

#define INIT_SEQ(seq) (seq).s = 0; (seq).l = (seq).m = 0

void 
seq_set_block_size(int size);

int 
seq_read_fasta(FILE *fp, seq_t *seq, char *locus, char *comment);

enum muttype_t {
    NOCHANGE = 0, 
    INSERT = 0x10, 
    SUBSTITUTE = 0x20, 
    DELETE = 0x30
};

typedef uint64_t mut_t;
extern mut_t mutmsk;
extern mut_t muttype_shift; 
extern mut_t ins_length_shift;
extern mut_t ins_length_mask;
extern mut_t ins_length_max;
extern mut_t ins_mask;

typedef struct {
    int l, m; /* length and maximum buffer size */
    mut_t *s; /* sequence */
} mutseq_t;

void 
mut_init();

void 
mut_diref(dwgsim_opt_t *opt, const seq_t *seq, mutseq_t *hap1, mutseq_t *hap2, int32_t contig_i, muts_txt_t *muts_txt, muts_bed_t *muts_bed);

// Columns:
// 1 - chromosome name
// 2 - position (one-based)
// 3 - reference base (dash if there was an insertion)
// 4 - variant allele (IUPAC code or insertion base(s))
// 5 - '-' for homozygous, '+' for heterozygous
void 
mut_print(const char *name, const seq_t *seq, mutseq_t *hap1, mutseq_t *hap2, FILE *fpout);

#endif
