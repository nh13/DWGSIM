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
#include "contigs.h"
#include "mut.h"
#include "dwgsim.h"
#include "mut_txt.h"
#include "mut_bed.h"
#include "mut_vcf.h"
#include "mut_input.h"

muts_input_t *muts_input_init(FILE *fp, contigs_t *c, int32_t type)
{
  muts_input_t *m = calloc(1, sizeof(muts_input_t));
  m->type = type;
  switch(type) {
    case MUT_INPUT_BED:
      m->data.bed = muts_bed_init(fp, c);
      break;
    case MUT_INPUT_TXT:
      m->data.txt = muts_txt_init(fp, c);
      break;
    case MUT_INPUT_VCF:
      m->data.vcf = muts_vcf_init(fp, c);
      break;
    default:
      fprintf(stderr, "Error: mutation input type unrecognized!\n");
      exit(1);
  }
  return m;
}

void muts_input_destroy(muts_input_t *m)
{
  switch(m->type) {
    case MUT_INPUT_BED:
      muts_bed_destroy(m->data.bed);
      break;
    case MUT_INPUT_TXT:
      muts_txt_destroy(m->data.txt);
      break;
    case MUT_INPUT_VCF:
      muts_vcf_destroy(m->data.vcf);
      break;
    default:
      fprintf(stderr, "Error: mutation input type unrecognized!\n");
      exit(1);
  }
  free(m);
}
