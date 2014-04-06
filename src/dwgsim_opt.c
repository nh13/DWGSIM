/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
   */

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <stdio.h>
#include <unistd.h>
#include <stdint.h>
#include <ctype.h>
#include <string.h>
#include "mut.h"
#include "dwgsim.h"
#include "dwgsim_opt.h"

dwgsim_opt_t* dwgsim_opt_init()
{
  dwgsim_opt_t *opt;
  opt = calloc(1, sizeof(dwgsim_opt_t));
  opt->e[0].start = opt->e[0].end = opt->e[1].start = opt->e[1].end = 0.02;
  opt->e[0].by = opt->e[1].by = 0;
  opt->is_inner = 0;
  opt->dist = 500;
  opt->std_dev = 50;
  opt->N = -1;
  opt->C = 100;
  opt->length[0] = opt->length[1] = 70;
  opt->mut_rate = 0.001;
  opt->mut_freq = 0.5;
  opt->indel_frac = 0.1;
  opt->indel_extend = 0.3;
  opt->indel_min = 1;
  opt->rand_read = 0.05;
  opt->data_type = ILLUMINA;
  opt->strandedness = 0;
  opt->max_n = 0;
  opt->flow_order = NULL;
  opt->flow_order_len = 0;
  opt->use_base_error = 0;
  opt->seed = -1;
  opt->muts_only = 0;
  opt->fixed_quality = NULL;
  opt->quality_std = 2.0;
  opt->fn_muts_input = NULL;
  opt->fn_muts_input_type = -1;
  opt->fn_regions_bed = NULL;
  opt->fp_mut = opt->fp_bfast = opt->fp_bwa1 = opt->fp_bwa2 = NULL;
  opt->fp_fa = opt->fp_fai = NULL;
  opt->read_prefix = NULL;

  return opt;
}

void dwgsim_opt_destroy(dwgsim_opt_t *opt)
{
  free(opt->fixed_quality);
  free(opt->fn_muts_input);
  free(opt->fn_regions_bed);
  free(opt->flow_order);
  free(opt->read_prefix);
  free(opt);
}

int dwgsim_opt_usage(dwgsim_opt_t *opt)
{
  mutseq_init_bounds();
  fprintf(stderr, "\n");
  fprintf(stderr, "Program: dwgsim (short read simulator)\n");
  fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
  fprintf(stderr, "Contact: Nils Homer <dnaa-help@lists.sourceforge.net>\n\n");
  fprintf(stderr, "Usage:   dwgsim [options] <in.ref.fa> <out.prefix>\n\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "         -e FLOAT      per base/color/flow error rate of the first read [from %.3f to %.3f by %.3f]\n", opt->e[0].start, opt->e[0].end, opt->e[0].by);
  fprintf(stderr, "         -E FLOAT      per base/color/flow error rate of the second read [from %.3f to %.3f by %.3f]\n", opt->e[1].start, opt->e[1].end, opt->e[1].by);
  fprintf(stderr, "         -i            use the inner distance instead of the outer distance for pairs [%s]\n", __IS_TRUE(opt->is_inner));
  fprintf(stderr, "         -d INT        %s distance between the two ends for pairs [%d]\n", 
          (0 == opt->is_inner) ? "outer" : "inner", opt->dist);
  fprintf(stderr, "         -s INT        standard deviation of the distance for pairs [%.3f]\n", opt->std_dev);
  fprintf(stderr, "         -N INT        number of read pairs (-1 to disable) [%lld]\n", (signed long long int)opt->N);
  fprintf(stderr, "         -C FLOAT      mean coverage across available positions (-1 to disable) [%.2lf]\n", opt->C);
  fprintf(stderr, "         -1 INT        length of the first read [%d]\n", opt->length[0]);
  fprintf(stderr, "         -2 INT        length of the second read [%d]\n", opt->length[1]);
  fprintf(stderr, "         -r FLOAT      rate of mutations [%.4f]\n", opt->mut_rate);
  fprintf(stderr, "         -F FLOAT      frequency of given mutation to simulate low fequency somatic mutations [%.4f]\n", opt->mut_freq);
  fprintf(stderr, "                           NB: freqeuncy F refers to the first strand of mutation, therefore mutations \n"); 
  fprintf(stderr, "                           on the second strand occour with a frequency of 1-F \n");
  fprintf(stderr, "         -R FLOAT      fraction of mutations that are indels [%.2f]\n", opt->indel_frac);
  fprintf(stderr, "         -X FLOAT      probability an indel is extended [%.2f]\n", opt->indel_extend);
  fprintf(stderr, "         -I INT        the minimum length indel [%d]\n", opt->indel_min);
  fprintf(stderr, "         -y FLOAT      probability of a random DNA read [%.2f]\n", opt->rand_read);
  fprintf(stderr, "         -n INT        maximum number of Ns allowed in a given read [%d]\n", opt->max_n);
  fprintf(stderr, "         -c INT        generate reads for [%d]:\n", opt->data_type);
  fprintf(stderr, "                           0: Illumina\n");
  fprintf(stderr, "                           1: SOLiD\n");
  fprintf(stderr, "                           2: Ion Torrent\n");
  fprintf(stderr, "         -S INT        generate reads [%d]:\n", opt->strandedness);
  fprintf(stderr, "                           0: default (opposite strand for Illumina, same strand for SOLiD/Ion Torrent)\n");
  fprintf(stderr, "                           1: same strand (mate pair)\n");
  fprintf(stderr, "                           2: opposite strand (paired end)\n");
  fprintf(stderr, "         -f STRING     the flow order for Ion Torrent data [%s]\n", (char*)opt->flow_order);
  fprintf(stderr, "         -B            use a per-base error rate for Ion Torrent data [%s]\n", __IS_TRUE(opt->use_base_error));
  fprintf(stderr, "         -H            haploid mode [%s]\n", __IS_TRUE(opt->is_hap));
  fprintf(stderr, "         -z INT        random seed (-1 uses the current time) [%d]\n", opt->seed);
  fprintf(stderr, "         -M            generate a mutations file only [%s]\n", __IS_TRUE(opt->muts_only));
  fprintf(stderr, "         -m FILE       the mutations txt file to re-create [%s]\n", (MUT_INPUT_TXT != opt->fn_muts_input_type) ? "not using" : opt->fn_muts_input);
  fprintf(stderr, "         -b FILE       the bed-like file set of candidate mutations [%s]\n", (MUT_INPUT_BED == opt->fn_muts_input_type) ? "not using" : opt->fn_muts_input);
  fprintf(stderr, "         -v FILE       the vcf file set of candidate mutations (use pl tag for strand) [%s]\n", (MUT_INPUT_VCF == opt->fn_muts_input_type) ? "not using" : opt->fn_muts_input);
  fprintf(stderr, "         -x FILE       the bed of regions to cover [%s]\n", (NULL == opt->fn_regions_bed) ? "not using" : opt->fn_regions_bed);
  fprintf(stderr, "         -P STRING     a read prefix to prepend to each read name [%s]\n", (NULL == opt->read_prefix) ? "not using" : opt->read_prefix);
  fprintf(stderr, "         -q STRING     a fixed base quality to apply (single character) [%s]\n", (NULL == opt->fixed_quality) ? "not using" : opt->fixed_quality);
  fprintf(stderr, "         -Q FLOAT      standard deviation of the base quality scores [%.2lf]\n", (NULL == opt->fixed_quality) ? opt->quality_std : 0.0);
  fprintf(stderr, "         -s INT        standard deviation of the distance for pairs [%.3f]\n", opt->std_dev);
  fprintf(stderr, "         -h            print this message\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Note: For SOLiD mate pair reads and BFAST, the first read is F3 and the second is R3. For SOLiD mate pair reads\n");
  fprintf(stderr, "and BWA, the reads in the first file are R3 the reads annotated as the first read etc.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Note: The longest supported insertion is %u.\n", UINT32_MAX);
  fprintf(stderr, "\n");
  return 1;
}

static void get_error_rate(const char *str, error_t *e)
{
  int32_t i;

  e->start = atof(str);
  for(i=0;i<strlen(str);i++) {
      if(',' == str[i] || '-' == str[i]) {
          break;
      }
  }
  if(i<strlen(str)-1) {
      i++;
      e->end = atof(str+i);
  }
  else {
      e->end = e->start;
  }
}

int32_t
dwgsim_opt_is_int(char *optarg)
{
  int32_t i;
  for(i=0;i<strlen(optarg);i++) {
      if(!isdigit(optarg[i])) return 0;
  }
  return 1;
}

int32_t
dwgsim_atoi(char *optarg, char flag)
{
  if(0 == dwgsim_opt_is_int(optarg)) {
      fprintf(stderr, "Error: command line option -%c is not a number [%s]\n", flag, optarg);
      exit(1);
  }
  return atoi(optarg);
}

int32_t
dwgsim_opt_parse(dwgsim_opt_t *opt, int argc, char *argv[]) 
{
  int32_t i;
  int c;
  int muts_input_type = 0;
  
  while ((c = getopt(argc, argv, "id:s:N:C:1:2:e:E:r:F:R:X:I:c:S:n:y:BHf:z:Mm:b:v:x:P:q:Q:h")) >= 0) {
      switch (c) {
        case 'i': opt->is_inner = 1; break;
        case 'd': opt->dist = dwgsim_atoi(optarg, 'd'); break;
        case 's': opt->std_dev = atof(optarg); break;
        case 'N': opt->N = dwgsim_atoi(optarg, 'N'); opt->C = -1; break;
        case 'C': opt->C = atof(optarg); opt->N = -1; break;
        case '1': opt->length[0] = dwgsim_atoi(optarg, '1'); break;
        case '2': opt->length[1] = dwgsim_atoi(optarg, '2'); break;
        case 'e': get_error_rate(optarg, &opt->e[0]); break;
        case 'E': get_error_rate(optarg, &opt->e[1]); break;
        case 'r': opt->mut_rate = atof(optarg); break;
        case 'F': opt->mut_freq = atof(optarg); break;
        case 'R': opt->indel_frac = atof(optarg); break;
        case 'X': opt->indel_extend = atof(optarg); break;
        case 'I': opt->indel_min = dwgsim_atoi(optarg, 'I'); break;
        case 'c': opt->data_type = dwgsim_atoi(optarg, 'c'); break;
        case 'S': opt->strandedness = dwgsim_atoi(optarg, 'S'); break;
        case 'n': opt->max_n = dwgsim_atoi(optarg, 'n'); break;
        case 'y': opt->rand_read = atof(optarg); break;
        case 'f': 
                  if(NULL != opt->flow_order) free(opt->flow_order);
                  opt->flow_order = (int8_t*)strdup(optarg);
                  break;
        case 'B': opt->use_base_error = 1; break;
        case 'H': opt->is_hap = 1; break;
        case 'h': return 0;
        case 'z': opt->seed = dwgsim_atoi(optarg, 'n'); break;
        case 'M': opt->muts_only = 1; break;
        case 'm': free(opt->fn_muts_input); opt->fn_muts_input = strdup(optarg); opt->fn_muts_input_type = MUT_INPUT_TXT; muts_input_type |= 0x1; break;
        case 'b': free(opt->fn_muts_input); opt->fn_muts_input = strdup(optarg); opt->fn_muts_input_type = MUT_INPUT_BED; muts_input_type |= 0x2; break;
        case 'v': free(opt->fn_muts_input); opt->fn_muts_input = strdup(optarg); opt->fn_muts_input_type = MUT_INPUT_VCF; muts_input_type |= 0x4; break;
        case 'x': free(opt->fn_regions_bed); opt->fn_regions_bed = strdup(optarg); break;
        case 'P': free(opt->read_prefix); opt->read_prefix = strdup(optarg); break;
        case 'q': opt->fixed_quality = strdup(optarg); break;
        case 'Q': opt->quality_std = atof(optarg); break;
        default: fprintf(stderr, "Unrecognized option: -%c\n", c); return 0;
      }
  }
  if (argc - optind < 2) return 0;

  __check_option(opt->is_inner, 0, 1, "-i");
  __check_option(opt->dist, 0, INT32_MAX, "-d");
  __check_option(opt->std_dev, 0, INT32_MAX, "-s");
  if(opt->N < 0 && opt->C < 0) {
      fprintf(stderr, "Must use one of -N or -C");
      return 0;
  }
  else if(0 < opt->N && 0 < opt->C) {
      fprintf(stderr, "Cannot use both -N or -C");
      return 0;
  }
  else if(0 < opt->N) {
      __check_option(opt->N, 1, INT32_MAX, "-N");
      __check_option(opt->C, INT32_MIN, -1, "-C");
  }
  else {
      __check_option(opt->N, INT32_MIN, -1, "-N");
      __check_option(opt->C, 0, INT32_MAX, "-C");
  }
  __check_option(opt->length[0], 1, INT32_MAX, "-1");
  __check_option(opt->length[1], 0, INT32_MAX, "-2");
  // error rate
  for(i=0;i<2;i++) {
      if(opt->e[i].start < 0.0 || 1.0 < opt->e[i].start) {
          fprintf(stderr, "End %s: the start error is out of range (-e)\n", (0 == i) ? "one" : "two");
          return 0;
      }
      if(opt->e[i].end < 0.0 || 1.0 < opt->e[i].end) {
          fprintf(stderr, "End %s: the end error is out of range (-e)\n", (0 == i) ? "one" : "two");
          return 0;
      }
      if(IONTORRENT == opt->data_type) {
          if(opt->e[i].end != opt->e[i].start) {
              fprintf(stderr, "End %s: a uniform error rate must be given for Ion Torrent data\n", (0 == i) ? "one" : "two");
              return 0;
          }
      }
  }
  __check_option(opt->mut_rate, 0, 1.0, "-r");
  __check_option(opt->indel_frac, 0, 1.0, "-R");
  __check_option(opt->indel_extend, 0, 1.0, "-X");
  __check_option(opt->indel_min, 1, INT32_MAX, "-I");
  __check_option(opt->data_type, 0, 2, "-c");
  __check_option(opt->strandedness, 0, 2, "-S");
  __check_option(opt->max_n, 0, INT32_MAX, "-n");
  __check_option(opt->rand_read, 0, 1.0, "-y");
  if(IONTORRENT == opt->data_type && NULL == opt->flow_order) {
      fprintf(stderr, "Error: command line option -f is required\n");
      return 0;
  }
  __check_option(opt->use_base_error, 0, 1, "-B");
  __check_option(opt->is_hap, 0, 1, "-H");

  if(NULL != opt->fixed_quality && 1 != strlen(opt->fixed_quality)) {
      fprintf(stderr, "Error: command line option -q requires one character\n");
      return 0;
  }
  __check_option(opt->quality_std, 0, INT32_MAX, "-Q");

  if(NULL != opt->read_prefix) {
      fprintf(stderr, "Warning: remember to use the -P option with dwgsim_eval\n");
  }

  if(0 < opt->length[1]){ //paired end / mate pair
    int s = opt->length[0] + opt->length[1] + opt->is_inner;
    double r = (s-opt->dist)*1.0/opt->std_dev;
    if(r>6){
      fprintf(stderr, "Error: command line option -d is too small for the read length (%d)\n",opt->dist);
      return 0;
    }
    if(r>4){
      fprintf(stderr, "Warning: command line option -d is small for the read length (%d). Generation speed could be affected.\n",opt->dist);
    }
  }

  switch(muts_input_type) {
    case 0x0:
    case 0x1:
    case 0x2:
    case 0x4:
      break;
    default:
      fprintf(stderr, "Error: -m/-b/-v cannot be used together\n");
      return 0;
      break;
  }
  
  // random seed
  srand48((-1 == opt->seed) ? time(0) : opt->seed);

  if(IONTORRENT == opt->data_type) {
      if(NULL != opt->flow_order) {
          // uniform error rates only (so far)
          if(opt->e[0].start != opt->e[0].end || opt->e[1].start != opt->e[1].end) {
              fprintf(stderr, "Error: non-uniform error rate not support for Ion Torrent data\n");
              return 0;
          }
          // update flow order
          opt->flow_order_len = strlen((char*)opt->flow_order);
          for(i=0;i<opt->flow_order_len;i++) {
              opt->flow_order[i] = nst_nt4_table[opt->flow_order[i]];
          }
      }
      else {
          fprintf(stderr, "Error: -f must be given for Ion Torrent data\n");
          return 0;
      }
  }
  // use base error rate
  if(IONTORRENT == opt->data_type && NULL != opt->flow_order && 1 == opt->use_base_error) {
      uint8_t *tmp_seq=NULL;
      uint8_t *tmp_seq_flow_mask=NULL;
      int32_t tmp_seq_mem, s, cur_n_err, n_err, counts;
      int32_t j, k;
      double sf = 0.0;
      for(i=0;i<2;i++) {
          if(opt->length[i] <= 0) continue;
          fprintf(stderr, "[dwgsim_core] Updating error rate for end %d\n", i+1);
          if(0 < i && opt->length[i] == opt->length[1-i]) {
              opt->e[i].start = opt->e[1-i].start;
              opt->e[i].end = opt->e[1-i].end;
              opt->e[i].by = opt->e[1-i].by;
              fprintf(stderr, "[dwgsim_core] Using scaling factor from previous end\n[dwgsim_core] Updated with scaling factor %.5lf\n", sf);
              continue;
          }
          tmp_seq_mem = opt->length[i]+2;
          tmp_seq = (uint8_t*)calloc(tmp_seq_mem, 1);
          tmp_seq_flow_mask = (uint8_t*)calloc(tmp_seq_mem, 1);
          n_err = counts = 0;
          for(j=0;j<ERROR_RATE_NUM_RANDOM_READS;j++) {
              if(0 == (j % 10000)) {
                  fprintf(stderr, "\r[dwgsim_core] %d", j);
              }
              for(k=0;k<opt->length[i];k++) {
                  tmp_seq[k] = (int)(drand48() * 4.0) & 3;
              }
              cur_n_err = 0;
              s = opt->length[i];
              s = generate_errors_flows(opt, &tmp_seq, &tmp_seq_flow_mask, &tmp_seq_mem, s, 0, opt->e[i].start, &cur_n_err);
              n_err += cur_n_err;
              counts += s;
          }
          //fprintf(stderr, "before %lf,%lf,%lf\n", opt->e[i].start, opt->e[i].by, opt->e[i].end); 
          sf = opt->e[i].start / (n_err / (1.0 * counts));
          opt->e[i].start = opt->e[i].end *= sf;
          opt->e[i].by = (opt->e[i].end - opt->e[i].start) / opt->length[i];
          //fprintf(stderr, "after %lf,%lf,%lf\n", opt->e[i].start, opt->e[i].by, opt->e[i].end); 
          free(tmp_seq);
          free(tmp_seq_flow_mask);
          fprintf(stderr, "\r[dwgsim_core] %d\n[dwgsim_core] Updated with scaling factor %.5lf!\n", j, sf);
      }
  }
  else {
      opt->e[0].by = (opt->e[0].end - opt->e[0].start) / opt->length[0];
      opt->e[1].by = (opt->e[1].end - opt->e[1].start) / opt->length[1];
  }
  
  __check_option(opt->muts_only, 0, 1, "-M");

  return 1;
}
