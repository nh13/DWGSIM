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
#include <stdio.h>
#include <assert.h>
#include <unistd.h>
#include <stdint.h>
#include <string.h>
#include "contigs.h"
#include "mut.h"
#include "dwgsim.h"
#include "mut_bed.h"

muts_bed_t *muts_bed_init(FILE *fp, contigs_t *c)
{
  muts_bed_t *m = NULL;
  int32_t i;
  char name[1024];
  uint32_t start, end;
  char type[1024];
  char bases[1024];
  uint32_t prev_contig=0, max_end=0;

  m = calloc(1, sizeof(muts_bed_t));
  m->n = 0;
  m->mem = 1024;
  m->muts = malloc(m->mem * sizeof(mut_bed_t));

  i = 0;
  // start is zero-based
  // one is one-based
  while(0 < fscanf(fp, "%s\t%u\t%u\t%s\t%s", name, &start, &end, bases, type)) {
      // find the contig
      while(i < c->n && 0 != strcmp(name, c->contigs[i].name)) {
          i++;
      }
      if(c->n == i) {
          fprintf(stderr, "Error: contig not found [%s]\n", name);
          exit(1);
      }
      else if(c->contigs[i].len <= start) {
          fprintf(stderr, "Error: start out of range [%s,%u]\n", name, start);
          exit(1);
      }
      else if(c->contigs[i].len < end) {
          fprintf(stderr, "Error: end out of range [%s,%u]\n", name, end);
          exit(1);
      }
      else if(end <= start) {
          fprintf(stderr, "Error: end <= start [%s,%u,%u]\n", name, start, end);
          exit(1);
      }
      else if(0 != strcmp("*", bases) && (end - start) != strlen(bases)) {
          fprintf(stderr, "Error: bases did not match start and end [%s,%u,%u,%s]\n", name, start, end, bases);
          exit(1);
      }
      else if(prev_contig == i && start+1 <= max_end) {
          fprintf(stderr, "Warning: overlapping entries, ignoring entry [%s\t%u\t%u\t%s\t%s]\n", name, start, end, bases, type);
          continue;
      }

      if(prev_contig != i || max_end < end) {
          prev_contig = i;
          max_end = end;
      }

      if(m->n == m->mem) {
          m->mem <<= 1;
          m->muts = realloc(m->muts, m->mem * sizeof(mut_bed_t)); 
      }
      m->muts[m->n].contig = i;
      m->muts[m->n].start = start;
      m->muts[m->n].end = end; // one-based

      m->muts[m->n].type = get_muttype(type);
      switch(m->muts[m->n].type) {
        case SUBSTITUTE:
          break;
        case INSERT:
          if(ins_length_max < end - start) {
              fprintf(stderr, "Error: insertion of length %d exceeded the maximum supported length of %d\n",
                      end - start,
                      (int32_t)ins_length_max);
              exit(1);
          }
          break;
        case DELETE:
          break;
        default:
          // error
          fprintf(stderr, "Error: mutation type unrecognized [%s]\n", type);
          exit(1);
      }
      m->muts[m->n].bases = strdup(bases);
      m->n++;
  }
  if(m->n < m->mem) {
      m->mem = m->n;
      m->muts = realloc(m->muts, m->mem * sizeof(mut_bed_t)); 
  }

  return m;
}

void muts_bed_destroy(muts_bed_t *m)
{
  int32_t i;
  for(i=0;i<m->n;i++) {
      free(m->muts[i].bases);
  }
  free(m->muts);
  free(m);
}
