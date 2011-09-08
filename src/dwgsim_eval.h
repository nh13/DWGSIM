#ifndef DWGSIM_EVAL_H_
#define DWGSIM_EVAL_H_

#define DWGSIM_EVAL_MAXQ 255
#define DWGSIM_EVAL_MINAS -5000

typedef struct {
    int32_t a; // alignment score or not
    int32_t b; // bwa or not
    int32_t c; // color space or not
    int32_t d; // divide by factor
    int32_t e; // print only alignments with # of errors
    int32_t g; // gap "wiggle"
    int32_t i; // indel only
    int32_t m; // multi-mapped
    int32_t n; // # of pe alignments
    int32_t p; // print incorrect alignments or not
    int32_t q; // consider only alignments with this mapping quality or greater
    int32_t s; // print only alignments with # of SNPs 
    int32_t z; // input reads are single end
    int32_t S; // input reads are in text SAM format
} dwgsim_eval_args_t;

// Actual value
enum {
    DWGSIM_EVAL_MAPPABLE   =0,
    DWGSIM_EVAL_UNMAPPABLE =1
};

// Prediction
enum {
    DWGSIM_EVAL_MAPPED_CORRECTLY   =0,
    DWGSIM_EVAL_MAPPED_INCORRECTLY =1,
    DWGSIM_EVAL_UNMAPPED           =2
};

typedef struct {
    int32_t *mc; // mappable && mapped correctly
    int32_t *mi; // mappable && mapped incorrectly
    int32_t *mu; // mappable && unmapped
    int32_t *um; // unmappable && mapped
    int32_t *uu; // unmappable && unmapped
    int32_t min_score, max_score;
} dwgsim_eval_counts_t;

void 
run(dwgsim_eval_args_t *args,
         int32_t num_files,
         char *files[]);
void
process_bam(dwgsim_eval_counts_t *counts,
                    dwgsim_eval_args_t *args,
                    bam_header_t *header,
                    bam1_t *b,
                    samfile_t *fp_out);
dwgsim_eval_counts_t *
dwgsim_eval_counts_init();
void
dwgsim_eval_counts_destroy(dwgsim_eval_counts_t *counts);
void 
dwgsim_eval_counts_add(dwgsim_eval_counts_t *counts, int32_t score, int32_t actual_value, int32_t predicted_value);
void 
dwgsim_eval_counts_print(dwgsim_eval_counts_t *counts, int32_t a, int32_t d, int32_t n);

#endif
