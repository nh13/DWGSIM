#ifndef DWGSIM_OPT_H
#define DWGSIM_OPT_H

#include "zlib.h"

#define ERROR_RATE_NUM_RANDOM_READS 1000000

#define OUTPUT_TYPE_ALL 0
#define OUTPUT_TYPE_BWA 1
#define OUTPUT_TYPE_BFAST 2


typedef struct {
    double start, by, end;
} error_t;

typedef struct {
    error_t e[2];
    int32_t is_inner;
    int32_t dist;
    double std_dev;
    int64_t N;
    double C;
    int32_t length[2];
    double mut_rate;
    double mut_freq;
    double indel_frac;
    double indel_extend;
    int32_t indel_min;
    double rand_read;
    int32_t max_n;
    int32_t data_type;
    int32_t strandedness;
    int32_t read_one_strand;
    int8_t *flow_order;
    int32_t flow_order_len;
    int32_t use_base_error;
    int32_t is_hap;
    int32_t seed;
    int32_t muts_only;
    char *fixed_quality;
    double quality_std;
    char *fn_muts_input;
    int32_t fn_muts_input_type;
    char *fn_regions_bed;
    FILE *fp_mut;
    FILE *fp_vcf;
    gzFile *fp_bfast;
    gzFile *fp_bwa1;
    gzFile *fp_bwa2;
    FILE *fp_fa;
    FILE *fp_fai;
    char *read_prefix;
    int32_t output_type;
} dwgsim_opt_t;

dwgsim_opt_t* 
dwgsim_opt_init();

void 
dwgsim_opt_destroy(dwgsim_opt_t *opt);

#define __check_option(_val, _min, _max, _opt) \
  if(_val < _min || _max < _val) { \
      fprintf(stderr, "Error: command line option %s was out of range\n", _opt); \
      return 0; \
  } 

int 
dwgsim_opt_usage(dwgsim_opt_t *opt);

int32_t
dwgsim_opt_parse(dwgsim_opt_t *opt, int argc, char *argv[]); 

#endif
