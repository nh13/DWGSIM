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
#include "mut_vcf.h"

#define BUFFER_L 1048576

muts_vcf_t *muts_vcf_init(FILE *fp, contigs_t *c)
{
  static int32_t warned = 0;
  muts_vcf_t *m = NULL;
  int32_t i, j;
  char name[1024]; // 1. #CHROM
  uint32_t pos; // 2. POS
  char id[1024]; // 3. ID
  char ref[1024]; // 4. REF
  char alt[1025]; // 5. ALT
  size_t n_read = -1;
  uint32_t prev_pos = 0, is_hap;
  int32_t ref_l, alt_l;
  char buffer[BUFFER_L];
  int32_t s, n, len;

  m = calloc(1, sizeof(muts_vcf_t));
  m->n = 0;
  m->mem = 1024;
  m->muts = malloc(m->mem * sizeof(mut_vcf_t));
  
  // read in some buffer
  n_read = fread(buffer, sizeof(char), BUFFER_L, fp);

  i = s = n = 0;
  while(0 != s || 0 != n_read) {
      len = n_read + s;
      s = 0;
      while(s < len) { // while more characters in the buffer
          // find the next newline
          for(n=s;n<len;n++) {
              if('\n' == buffer[n] || '\r' == buffer[n]) break;
          }
          if(n == len && 0 < n_read) { // there are more characters to read
              break;
          }
          else if(s == n) { // newline found, skip
              s++; 
              continue;
          }
          else if(buffer[s] == '#') { // skip headers
              s = n;
              continue;
          }

          // process
          if(EOF == sscanf(buffer+s, "%s\t%u\t%s\t%s\t%s", name, &pos, id, ref, alt)) {
              fprintf(stderr, "Error: VCF parsing error\n"); 
              exit(1);
          }
          // find ploidy in the "pl" (lowercase) tag
          is_hap = 4;
          while(s+4 < n) { // [\t;]pl=[1-3]
              if(('\t' == buffer[s] || ';' == buffer[s]) && 'p' == buffer[s+1] && 'l' == buffer[s+2] && '=' == buffer[s+3]) {
                  s += 4; // skip over '[\t;]pl='
                  switch(buffer[s]) {
                    case '1':
                      is_hap = 1;
                      break;
                    case '2':
                      is_hap = 2;
                      break;
                    case '3':
                      is_hap = 3;
                      break;
                    default:
                      fprintf(stderr, "Error: Could not determine the strand of the mutation from the 'pl' tag.\n");
                      exit(1);
                      break;
                  }
                  break;
              }
              s++;
          }
          if(4 == is_hap && 0 == warned) {
              fprintf(stderr, "Warning: strand of the mutation not found; please use the 'pl' tag.\n"); 
              warned = 1;
              is_hap = 3;
          }

          s = n+1; // move past the next newline

          // find the contig
          while(i < c->n && 0 != strcmp(name, c->contigs[i].name)) {
              i++;
              prev_pos = 0;
          }
          if(c->n == i) {
              fprintf(stderr, "Error: contig not found [%s]\n", name);
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

          ref_l = strlen(ref);
          alt_l = strlen(alt);

          if(1 == ref_l && ref[0] == '.') {
              ref[0] = '\0';
              ref_l = 0;
          }
          if(1 == alt_l && alt[0] == '.') {
              alt[0] = '\0';
              alt_l = 0;
          }
          if(0 == alt_l && 0 == ref_l) {
              fprintf(stderr, "Error: empty alleles\n");
              exit(1);
          }

          // TODO: support multiple alleles
          for(j=0;j<alt_l;j++) {
              if(',' == alt[j]) {
                  fprintf(stderr, "Error: multiple alleles are not supported\n");
                  exit(1);
              }
          }
          
          // to upper, and check for non-ACGT bases
          for(j=0;j<ref_l;j++) {
              ref[j] = "ACGTN"[nst_nt4_table[(int)ref[j]]];
              if('N' == ref[j]) {
                  fprintf(stderr, "Error: non-ACGT base found\n");
                  exit(1);
              }
          }
          for(j=0;j<alt_l;j++) {
              alt[j] = "ACGTN"[nst_nt4_table[(int)alt[j]]];
              if('N' == alt[j]) {
                  fprintf(stderr, "Error: non-ACGT base found\n");
                  exit(1);
              }
          }

          if(ref_l == alt_l) { // SNP
              while(m->mem <= m->n + ref_l - 1) {
                  m->mem <<= 1;
                  m->muts = realloc(m->muts, m->mem * sizeof(mut_vcf_t)); 
              }
              for(j=0;j<ref_l;j++) { // multiple snps
                  m->muts[m->n].contig = i;
                  m->muts[m->n].pos = pos + j;
                  m->muts[m->n].type = SUBSTITUTE;
                  m->muts[m->n].is_hap = is_hap;
                  m->muts[m->n].bases = malloc(sizeof(char) * 2);
                  m->muts[m->n].bases[0] = alt[j];
                  m->muts[m->n].bases[1] = '\0';
                  m->n++;
              }
          }
          else if(ref_l < alt_l) { // INS
              if(m->n == m->mem) {
                  m->mem <<= 1;
                  m->muts = realloc(m->muts, m->mem * sizeof(mut_vcf_t)); 
              }
              // skip over match bases
              for(j=0;j<ref_l;j++,pos++) {
                  if(ref[j] != alt[j]) break;
              }
              //pos--; // make 0-based for INS
              m->muts[m->n].contig = i;
              m->muts[m->n].pos = pos;
              m->muts[m->n].type = INSERT;
              m->muts[m->n].is_hap = is_hap;
              m->muts[m->n].bases = strdup(alt+j);
              m->n++;
          }
          else { // DEL
              // skip over match bases
              for(j=0;j<alt_l;j++,pos++) {
                  if(ref[j] != alt[j]) break;
              }
              if(j == ref_l) {
                  fprintf(stderr, "Error: no deleted bases\n");
                  exit(1);
              }
              while(m->mem < m->n + (ref_l - j)) {
                  m->mem <<= 1;
                  m->muts = realloc(m->muts, m->mem * sizeof(mut_vcf_t)); 
              }
              for(;j<ref_l;j++,pos++) {
                  m->muts[m->n].contig = i;
                  m->muts[m->n].pos = pos;
                  m->muts[m->n].type = DELETE;
                  m->muts[m->n].is_hap = is_hap;
                  m->muts[m->n].bases = NULL;
                  m->n++;
              }
          }

          prev_pos = pos;

          // move to 'n'
          s = n;
      }
      n++; // move past last character
      // shift unused characters down 
      for(s = 0; n < len; n++, s++) {
          buffer[s] = buffer[n];
      }
      // read in some buffer
      n_read = fread(buffer + s, sizeof(char), BUFFER_L - s, fp);
  }
  if(m->n < m->mem) {
      m->mem = m->n;
      m->muts = realloc(m->muts, m->mem * sizeof(mut_vcf_t)); 
  }

  return m;
}

void muts_vcf_destroy(muts_vcf_t *m)
{
  int32_t i;
  for(i=0;i<m->n;i++) {
      free(m->muts[i].bases);
  }
  free(m->muts);
  free(m);
}
