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
#include "regions_bed.h"

regions_bed_txt *regions_bed_init(FILE *fp, contigs_t *c)
{
  regions_bed_txt *r = NULL;
  char name[1024];
  uint32_t start, end, len;
  int32_t i, prev_contig, prev_start, prev_end, b;

  r = calloc(1, sizeof(regions_bed_txt));
  r->n = 0;
  r->mem = 4;
  r->contig = malloc(r->mem * sizeof(uint32_t));
  r->start = malloc(r->mem * sizeof(uint32_t));
  r->end = malloc(r->mem * sizeof(uint32_t));
  
  i = 0;
  prev_contig = prev_start = prev_end = -1;
  while(0 < fscanf(fp, "%s\t%u\t%u", name, &start, &end)) {
      len = end - start + 1;
      // find the contig
      while(i < c->n && 0 != strcmp(name, c->contigs[i].name)) {
          i++;
      }
      if(c->n == i) {
          fprintf(stderr, "Error: contig not found [%s].  Are you sure your BED is coordinate sorted?\n", name);
          exit(1);
      }
      else if(c->contigs[i].len < start) {
          fprintf(stderr, "Error: start out of range [%s,%u]\n", name, start);
          exit(1);
      }
      else if(c->contigs[i].len < end) {
          fprintf(stderr, "Error: end out of range [%s,%u]\n", name, end);
          exit(1);
      }
      else if(end < start) {
          fprintf(stderr, "Error: end < start [%s,%u,%u]\n", name, start, end);
          exit(1);
      }
      else if(prev_contig == i && start < prev_start) {
          fprintf(stderr, "Error: the input was not sorted [%s,%u,%u,%u]\n", name, start, end, len);
          exit(1);
      }
      
      if(end - start + 1 != len) {
          fprintf(stderr, "Warning: len != end - start + 1 [%s,%u,%u,%u]\n", name, start, end, len);
      }
      
      if(prev_contig == i && start <= prev_end && prev_start <= start) {
          if(prev_end < end) {
              r->end[r->n-1] = end;
              prev_end = end;
          }
      }
      else {
          prev_contig = i;
          prev_start = start;
          prev_end = end;
          while(r->mem <= r->n) {
              r->mem <<= 1;
              r->contig = realloc(r->contig, r->mem * sizeof(uint32_t));
              r->start = realloc(r->start, r->mem * sizeof(uint32_t));
              r->end = realloc(r->end, r->mem * sizeof(uint32_t));
          }
          r->contig[r->n] = i;
          r->start[r->n] = start;
          r->end[r->n] = end;
          r->n++;
      }
      // move to the end of the line
      while(EOF != (b = fgetc(fp))) {
          if('\n' == b || '\r' == b) break;
      }
  }
  return r;
}

void regions_bed_destroy(regions_bed_txt *r)
{
  free(r->contig);
  free(r->start);
  free(r->end);
  free(r);
}

int32_t regions_bed_query(regions_bed_txt *r, uint32_t contig, uint32_t start, uint32_t end) 
{
  int32_t low, high, mid;
  if(NULL == r) return 1;

  low = 0;
  high = r->n-1;
  
  while(low <= high) {
      mid = (low + high) / 2;
      if(contig < r->contig[mid] ||
         (contig == r->contig[mid] && start < r->start[mid])) {
          high = mid - 1;
      }
      else if(r->contig[mid] < contig ||
              (r->contig[mid] == contig && r->end[mid] < end)) {
          low = mid + 1;
      }
      else if(r->contig[mid] == contig && r->start[mid] <= start && end <= r->end[mid]) {
          return 1;
      }
      else {
          break;
      }
  }
  return 0;
}
