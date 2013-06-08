#ifndef MUT_H
#define MUT_H

#include "contigs.h"
#include "mut_txt.h"
#include "mut_bed.h"
#include "mut_input.h"
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
extern mut_t mut_and_type_mask;
extern mut_t muttype_shift; 
extern mut_t ins_length_shift;
extern mut_t ins_length_mask;
extern mut_t ins_length_max;
extern mut_t ins_long_length_max;
extern mut_t ins_mask;

typedef struct {
    int l, m; /* length and maximum buffer size */
    mut_t *s; /* sequence */
    uint8_t **ins; /* long insertions */
    int ins_l, ins_m; /* length and maximum buffer size for long insertions */
} mutseq_t;

void
mutseq_init_bounds();

mutseq_t *
mutseq_init();

void
mutseq_destroy(mutseq_t *seq);

inline int32_t
mut_get_ins_bytes(int32_t n);

inline uint8_t*
mut_get_ins_long_n(uint8_t *ins, uint32_t *n);

mut_t
mut_get_ins_length(mutseq_t *seq, int32_t i);

int32_t
mut_get_ins(mutseq_t *seq, int32_t i, mut_t *n, mut_t *ins);

void 
mut_diref(dwgsim_opt_t *opt, const seq_t *seq, mutseq_t *hap1, mutseq_t *hap2, int32_t contig_i, muts_input_t *muts_input);

void mut_print_header_pre(FILE *fpout_vcf);
void mut_print_header_post(FILE *fpout_vcf);
void mut_print_header_contig(FILE *fpout_vcf, const char *name, int32_t length);

// Columns:
// 1 - chromosome name
// 2 - position (one-based)
// 3 - reference base (dash if there was an insertion)
// 4 - variant allele (IUPAC code or insertion base(s))
// 5 - '-' for homozygous, '+' for heterozygous
void 
mut_print(const char *name, const seq_t *seq, mutseq_t *hap1, mutseq_t *hap2, FILE *fpout_txt, FILE *fpout_vcf);

// 0 - 0
// 1 - 1
// 5 - 2
// 9 - 3
// ...
#define mut_packed_len(_n) ((_n + 3) >> 2)

#endif
