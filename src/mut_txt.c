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

muts_txt_t *muts_txt_init(FILE *fp, contigs_t *c)
{
  muts_txt_t *m = NULL;
  int32_t i;
  char name[1024];
  uint32_t pos, prev_pos = 0, is_hap;
  char ref, mut[1024];

  m = calloc(1, sizeof(muts_txt_t));
  m->n = 0;
  m->mem = 1024;
  m->muts = malloc(m->mem * sizeof(mut_txt_t));

  i = 0;
  while(0 < fscanf(fp, "%s\t%u\t%c\t%s\t%d", name, &pos, &ref, mut, &is_hap)) {
      // find the contig
      while(i < c->n && 0 != strcmp(name, c->contigs[i].name)) {
          i++;
          prev_pos = 0;
      }
      if(c->n == i) {
          fprintf(stderr, "Error: mutation contig not found or out of order [%s]\n", name);
          exit(1);
      }
      else if(pos <= 0 || c->contigs[i].len < pos) {
          fprintf(stderr, "Error: start out of range [%s,%u]\n", name, pos);
          exit(1);
      }
      else if(pos < prev_pos) {
          fprintf(stderr, "Error: out of order [%s,%u]\n", name, pos);
          exit(1);
      }

      if(m->n == m->mem) {
          m->mem <<= 1;
          m->muts = realloc(m->muts, m->mem * sizeof(mut_txt_t)); 
      }
      m->muts[m->n].contig = i;
      m->muts[m->n].pos = pos;
      m->muts[m->n].bases = strdup(mut);

      if('-' == ref && '-' != mut[0]) {
          m->muts[m->n].type = INSERT;
      }
      else if('-' != ref && '-' == mut[0]) {
          m->muts[m->n].type = DELETE;
      }
      else if('-' != ref && '-' != mut[0]) {
          m->muts[m->n].type = SUBSTITUTE;
          if(is_hap < 3) { // heterozygous
              if(nst_nt4_table[(int)m->muts[m->n].bases[0]] < 4) {
                  fprintf(stderr, "Error: heterozygous bases must be in IUPAC form\n");
                  exit(1);
              }
              m->muts[m->n].bases[0] = iupac_and_base_to_mut(m->muts[m->n].bases[0], ref);
              if('X' == m->muts[m->n].bases[0]) {
                  fprintf(stderr, "Error: out of range\n");
                  exit(1);
              }
              m->muts[m->n].bases[1] = '\0';
          }
      }
      else {
          fprintf(stderr, "Error: out of range\n");
          exit(1);
      }
      
      m->muts[m->n].is_hap = is_hap;

      m->n++;
  }
  if(m->n < m->mem) {
      m->mem = m->n;
      m->muts = realloc(m->muts, m->mem * sizeof(mut_txt_t)); 
  }

  return m;
}

void muts_txt_destroy(muts_txt_t *m)
{
  int32_t i;
  for(i=0;i<m->n;i++) {
      free(m->muts[i].bases);
  }
  free(m->muts);
  free(m);
}
