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

#ifndef MUT_INPUT_H
#define MUT_INPUT_H

#include "mut_bed.h"
#include "mut_txt.h"
#include "mut_vcf.h"

enum {
    MUT_INPUT_BED = 0,
    MUT_INPUT_TXT = 1,
    MUT_INPUT_VCF = 2
};

typedef struct {
    int32_t type;
    union {
        muts_bed_t *bed;
        muts_txt_t *txt;
        muts_vcf_t *vcf;
    } data;
} muts_input_t;

muts_input_t *muts_input_init(FILE *fp, contigs_t *c, int32_t type);

void muts_input_destroy(muts_input_t *m);

#endif
