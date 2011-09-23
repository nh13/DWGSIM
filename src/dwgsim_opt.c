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
#include "dwgsim.h"
#include "dwgsim_opt.h"

dwgsim_opt_t* dwgsim_opt_init()
{
  dwgsim_opt_t *opt;
  opt = calloc(1, sizeof(dwgsim_opt_t));
  opt->e[0].start = opt->e[0].end = opt->e[1].start = opt->e[1].end = 0.02;
  opt->e[0].by = opt->e[1].by = 0;
  opt->dist = 500;
  opt->std_dev = 50;
  opt->N = -1;
  opt->C = 100;
  opt->length[0] = opt->length[1] = 70;
  opt->mut_rate = 0.001;
  opt->indel_frac = 0.1;
  opt->indel_extend = 0.3;
  opt->rand_read = 0.05;
  opt->data_type = ILLUMINA;
  opt->max_n = 0;
  opt->flow_order = NULL;
  opt->flow_order_len = 0;
  opt->use_base_error = 0;
  opt->seed = -1;
  opt->fn_muts_txt = NULL;
  opt->fn_muts_bed = NULL;
  opt->fn_regions_bed = NULL;
  opt->fp_mut = opt->fp_bfast = opt->fp_bwa1 = opt->fp_bwa2 = NULL;
  opt->fp_fa = opt->fp_fai = NULL;

  return opt;
}

void dwgsim_opt_destroy(dwgsim_opt_t *opt)
{
  free(opt->fn_muts_txt);
  free(opt->fn_muts_bed);
  free(opt->fn_regions_bed);
  free(opt->flow_order);
  free(opt);
}

int dwgsim_opt_usage(dwgsim_opt_t *opt)
{
  fprintf(stderr, "\n");
  fprintf(stderr, "Program: dwgsim (short read simulator)\n");
  fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
  fprintf(stderr, "Contact: Nils Homer <dnaa-help@lists.sourceforge.net>\n\n");
  fprintf(stderr, "Usage:   dwgsim [options] <in.ref.fa> <out.prefix>\n\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "         -e FLOAT      per base/color/flow error rate of the first read [from %.3f to %.3f by %.3f]\n", opt->e[0].start, opt->e[0].end, opt->e[0].by);
  fprintf(stderr, "         -E FLOAT      per base/color/flow error rate of the second read [from %.3f to %.3f by %.3f]\n", opt->e[1].start, opt->e[1].end, opt->e[1].by);
  fprintf(stderr, "         -d INT        inner distance between the two ends [%d]\n", opt->dist);
  fprintf(stderr, "         -s INT        standard deviation [%.3f]\n", opt->std_dev);
  fprintf(stderr, "         -N INT        number of read pairs (-1 to disable) [%lld]\n", opt->N);
  fprintf(stderr, "         -C FLOAT      mean coverage across available positions (-1 to disable) [%.2lf]\n", opt->C);
  fprintf(stderr, "         -1 INT        length of the first read [%d]\n", opt->length[0]);
  fprintf(stderr, "         -2 INT        length of the second read [%d]\n", opt->length[1]);
  fprintf(stderr, "         -r FLOAT      rate of mutations [%.4f]\n", opt->mut_rate);
  fprintf(stderr, "         -R FLOAT      fraction of mutations that are indels [%.2f]\n", opt->indel_frac);
  fprintf(stderr, "         -X FLOAT      probability an indel is extended [%.2f]\n", opt->indel_extend);
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
  fprintf(stderr, "         -m FILE       the mutations txt file to re-create [%s]\n", (NULL == opt->fn_muts_txt) ? "not using" : opt->fn_muts_txt);
  fprintf(stderr, "         -b FILE       the bed-like set of candidate mutations [%s]\n", (NULL == opt->fn_muts_bed) ? "not using" : opt->fn_muts_bed);
  fprintf(stderr, "         -x FILE       the bed of regions to cover [%s]\n", (NULL == opt->fn_regions_bed) ? "not using" : opt->fn_regions_bed);
  fprintf(stderr, "         -h            print this message\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Note: For SOLiD mate pair reads and BFAST, the first read is F3 and the second is R3. For SOLiD mate pair reads\n");
  fprintf(stderr, "and BWA, the reads in the first file are R3 the reads annotated as the first read etc.\n");
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
dwgsim_opt_parse(dwgsim_opt_t *opt, int argc, char *argv[]) 
{
  int c;
  while ((c = getopt(argc, argv, "d:s:N:C:1:2:e:E:r:R:X:c:S:n:y:BHf:z:m:b:x:h")) >= 0) {
      switch (c) {
        case 'd': opt->dist = atoi(optarg); break;
        case 's': opt->std_dev = atof(optarg); break;
        case 'N': opt->N = atoi(optarg); opt->C = -1; break;
        case 'C': opt->C = atof(optarg); opt->N = -1; break;
        case '1': opt->length[0] = atoi(optarg); break;
        case '2': opt->length[1] = atoi(optarg); break;
        case 'e': get_error_rate(optarg, &opt->e[0]); break;
        case 'E': get_error_rate(optarg, &opt->e[1]); break;
        case 'r': opt->mut_rate = atof(optarg); break;
        case 'R': opt->indel_frac = atof(optarg); break;
        case 'X': opt->indel_extend = atof(optarg); break;
        case 'c': opt->data_type = atoi(optarg); break;
        case 'S': opt->strandedness = atoi(optarg); break;
        case 'n': opt->max_n = atoi(optarg); break;
        case 'y': opt->rand_read = atof(optarg); break;
        case 'f': 
                  if(NULL != opt->flow_order) free(opt->flow_order);
                  opt->flow_order = (int8_t*)strdup(optarg);
                  break;
        case 'B': opt->use_base_error = 1; break;
        case 'H': opt->is_hap = 1; break;
        case 'h': return 0;
        case 'z': opt->seed = atoi(optarg); break;
        case 'm': free(opt->fn_muts_txt); opt->fn_muts_txt = strdup(optarg); break;
        case 'b': free(opt->fn_muts_bed); opt->fn_muts_bed = strdup(optarg); break;
        case 'x': free(opt->fn_regions_bed); opt->fn_regions_bed = strdup(optarg); break;
        default: fprintf(stderr, "Unrecognized option: -%c\n", c); return 0;
      }
  }
  if (argc - optind < 2) return 0;

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
      __check_option(opt->C, 1, INT32_MAX, "-C");
  }
  __check_option(opt->length[0], 1, INT32_MAX, "-1");
  __check_option(opt->length[1], 0, INT32_MAX, "-2");
  // error rate
  if(IONTORRENT == opt->data_type) {
      if(opt->e[0].end != opt->e[0].start) {
          fprintf(stderr, "End one: a uniform error rate must be given for Ion Torrent data");
          return 0;
      }
      if(opt->e[1].end != opt->e[1].start) {
          fprintf(stderr, "End two: a uniform error rate must be given for Ion Torrent data");
          return 0;
      }
  }
  __check_option(opt->mut_rate, 0, 1.0, "-r");

  return 1;
}
