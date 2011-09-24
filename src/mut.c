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
#include "mut_txt.h"
#include "mut_bed.h"
#include "regions_bed.h"
#include "dwgsim_opt.h"
#include "mut.h"

static int SEQ_BLOCK_SIZE = 512;

void seq_set_block_size(int size)
{
    SEQ_BLOCK_SIZE = size;
}

int seq_read_fasta(FILE *fp, seq_t *seq, char *locus, char *comment)
{
  int c, l, max;
  char *p;

  c = 0;
  while (!feof(fp) && fgetc(fp) != '>');
  if (feof(fp)) return -1;
  p = locus;
  while (!feof(fp) && (c = fgetc(fp)) != ' ' && c != '\t' && c != '\n')
    if (c != '\r') *p++ = c;
  *p = '\0';
  if (comment) {
      p = comment;
      if (c != '\n') {
          while (!feof(fp) && ((c = fgetc(fp)) == ' ' || c == '\t'));
          if (c != '\n') {
              *p++ = c;
              while (!feof(fp) && (c = fgetc(fp)) != '\n')
                if (c != '\r') *p++ = c;
          }
      }
      *p = '\0';
  } else if (c != '\n') while (!feof(fp) && fgetc(fp) != '\n');
  l = 0; max = seq->m;
  while (!feof(fp) && (c = fgetc(fp)) != '>') {
      if (isalpha(c) || c == '-' || c == '.') {
          if (l + 1 >= max) {
              max += SEQ_BLOCK_SIZE;
              seq->s = (unsigned char*)realloc(seq->s, sizeof(char) * max);
          }
          seq->s[l++] = (unsigned char)c;
      }
  }
  if (c == '>') ungetc(c,fp);
  seq->s[l] = 0;
  seq->m = max; seq->l = l;
  return l;
} 

mut_t mutmsk;
mut_t muttype_shift;
mut_t ins_length_shift;
mut_t ins_length_mask;
mut_t ins_length_max;
mut_t ins_mask;

void mut_init()
{
  int32_t i;

  mutmsk = (mut_t)0x30;
  muttype_shift = 6; // bits 5-6 store the mutation type
  ins_length_shift = 59; // bits 60-64 store the insertion length
  ins_length_mask = 0x1F; // bits 60-64 store the insertion length
  ins_length_max = ((ins_length_shift - muttype_shift) >> 1);
  if(ins_length_mask < ins_length_max) { // we exceed storing the length of the indel
      ins_length_max = ins_length_mask;
  }
  // create the insertion mask based on the maximum length insertion
  ins_mask = 0;
  for(i=0;i<ins_length_max;i++) {
      ins_mask = (ins_mask << 2) | 0x3; 
  }
}

void mut_diref(dwgsim_opt_t *opt, const seq_t *seq, mutseq_t *hap1, mutseq_t *hap2, int32_t contig_i, muts_txt_t *muts_txt, muts_bed_t *muts_bed)
{
  int32_t i, j, deleting = 0, deletion_length = 0;
  mutseq_t *ret[2];

  ret[0] = hap1; ret[1] = hap2;
  ret[0]->l = seq->l; ret[1]->l = seq->l;
  ret[0]->m = seq->m; ret[1]->m = seq->m;
  ret[0]->s = (mut_t *)calloc(seq->m, sizeof(mut_t));
  ret[1]->s = (mut_t *)calloc(seq->m, sizeof(mut_t));

  if(NULL == muts_bed && NULL == muts_txt) {
      for (i = 0; i != seq->l; ++i) {
          mut_t c;
          c = ret[0]->s[i] = ret[1]->s[i] = (mut_t)nst_nt4_table[(int)seq->s[i]];
          if (deleting) {
              if (deletion_length < opt->indel_min || drand48() < opt->indel_extend) {
                  if (deleting & 1) ret[0]->s[i] |= DELETE|c;
                  if (deleting & 2) ret[1]->s[i] |= DELETE|c;
                  deletion_length++;
                  continue;
              } else deleting = deletion_length = 0;
          }
          if (c < 4 && drand48() < opt->mut_rate) { // mutation
              if (drand48() >= opt->indel_frac) { // substitution
                  double r = drand48();
                  c = (c + (mut_t)(r * 3.0 + 1)) & 3;
                  if (opt->is_hap || drand48() < 0.333333) { // hom
                      ret[0]->s[i] = ret[1]->s[i] = SUBSTITUTE|c;
                  } else { // het
                      ret[drand48()<0.5?0:1]->s[i] = SUBSTITUTE|c;
                  }
              } else { // indel
                  if (drand48() < 0.5) { // deletion
                      if (opt->is_hap || drand48() < 0.3333333) { // hom-del
                          ret[0]->s[i] = ret[1]->s[i] = DELETE|c;
                          deleting = 3;
                      } else { // het-del
                          deleting = drand48()<0.5?1:2;
                          ret[deleting-1]->s[i] = DELETE|c;
                      }
                      deletion_length = 1;
                  } else { // insertion
                      mut_t num_ins = 0, ins = 0;
                      do {
                          num_ins++;
                          ins = (ins << 2) | (mut_t)(drand48() * 4.0);
                      } while (num_ins < ins_length_max && (num_ins < opt->indel_min || drand48() < opt->indel_extend));
                      assert(0 < num_ins);

                      if (opt->is_hap || drand48() < 0.333333) { // hom-ins
                          ret[0]->s[i] = ret[1]->s[i] = (num_ins << ins_length_shift) | (ins << muttype_shift) | INSERT | c;
                      } else { // het-ins
                          ret[drand48()<0.5?0:1]->s[i] = (num_ins << ins_length_shift) | (ins << muttype_shift) | INSERT | c;
                      }
                  }
              }
          }
      }
  }
  else if(NULL != muts_txt) {
      // seed
      for (i = 0; i != seq->l; ++i) {
          ret[0]->s[i] = ret[1]->s[i] = (mut_t)nst_nt4_table[(int)seq->s[i]];
      }
      for (i =0; i < muts_txt->n; ++i) {
          if (muts_txt->muts[i].contig == contig_i) {
              int8_t type = muts_txt->muts[i].type;
              uint32_t pos = muts_txt->muts[i].pos;
              int8_t is_hap = muts_txt->muts[i].is_hap;
              mut_t c = (mut_t)nst_nt4_table[(int)seq->s[pos-1]];

              if (DELETE == type) {
                  if (is_hap & 1) ret[0]->s[pos-1] |= DELETE|c;
                  if (is_hap & 2) ret[1]->s[pos-1] |= DELETE|c;
              }
              else if (SUBSTITUTE == type) {
                  if (is_hap & 1) ret[0]->s[pos-1] = SUBSTITUTE|nst_nt4_table[(int)muts_txt->muts[i].bases[0]];
                  if (is_hap & 2) ret[1]->s[pos-1] = SUBSTITUTE|nst_nt4_table[(int)muts_txt->muts[i].bases[0]];
              }
              else if (INSERT == type) {
                  mut_t num_ins = 0, ins = 0;
                  for (j = strlen(muts_txt->muts[i].bases)-1; 0 <= j; --j) {
                      if(ins_length_max <= num_ins) break;
                      num_ins++;
                      ins = (ins << 2) | nst_nt4_table[(int)muts_txt->muts[i].bases[j]];
                  } 
                  assert(0 < num_ins);

                  if (is_hap & 1) ret[0]->s[pos-1] = (num_ins << ins_length_shift) | (ins << muttype_shift) | INSERT | c;
                  if (is_hap & 2) ret[1]->s[pos-1] = (num_ins << ins_length_shift) | (ins << muttype_shift) | INSERT | c;
              }
          }
      }
  }
  else { // BED
      // seed
      for (i = 0; i != seq->l; ++i) {
          ret[0]->s[i] = ret[1]->s[i] = (mut_t)nst_nt4_table[(int)seq->s[i]];
      }
      // mutates exactly based on a BED file
      for (i = 0; i < muts_bed->n; ++i) {
          if (muts_bed->muts[i].contig == contig_i) {
              int32_t has_bases = 1, is_hom = 0;
              int32_t which_hap = 0;
              mut_t c;

              // does this mutation have random bases?
              if (0 == strcmp("*", muts_bed->muts[i].bases)) has_bases = 0; // random bases

              // het or hom?
              if (opt->is_hap || drand48() < 0.333333) {
                  is_hom = 1; // hom
              }
              else {
                  which_hap = drand48()<0.5?0:1;
              }

              // mutate
              if (SUBSTITUTE == muts_bed->muts[i].type) {
                  for (j = muts_bed->muts[i].start; j < muts_bed->muts[i].end; ++j) { // for each base
                      c = (mut_t)nst_nt4_table[(int)seq->s[j]];
                      if (0 == has_bases) { // random DNA base
                          double r = drand48();
                          c = (c + (mut_t)(r * 3.0 + 1)) & 3;
                      }
                      else {
                          c = (mut_t)nst_nt4_table[(int)muts_bed->muts[i].bases[j - muts_bed->muts[i].start]]; 
                      }
                      if (1 == is_hom) {
                          ret[0]->s[j] = ret[1]->s[j] = SUBSTITUTE|c;
                      } else { // het
                          ret[which_hap]->s[j] = SUBSTITUTE|c;
                      }
                  }
              }
              else if (DELETE == muts_bed->muts[i].type) {
                  for (j = muts_bed->muts[i].start; j < muts_bed->muts[i].end; ++j) { // for each base
                      c = (mut_t)nst_nt4_table[(int)seq->s[j]];
                      if (1 == is_hom) {
                          ret[0]->s[j] = ret[1]->s[j] = DELETE|c;
                      } else { // het-del
                          ret[which_hap]->s[j] = DELETE|c;
                      }
                  }
              }
              else if (INSERT == muts_bed->muts[i].type) {
                  mut_t num_ins = 0, ins = 0;

                  c = (mut_t)nst_nt4_table[(int)seq->s[muts_bed->muts[i].start]];
                  for (j = muts_bed->muts[i].end-1;
                       muts_bed->muts[i].start <= j && num_ins < ins_length_max; 
                       --j) { // for each base
                      num_ins++;
                      if(0 == has_bases) {
                          ins = (ins << 2) | (mut_t)(drand48() * 4.0);
                      }
                      else {
                          ins = (ins << 2) | (mut_t)(nst_nt4_table[(int)muts_bed->muts[i].bases[j - muts_bed->muts[i].start]]);
                      }
                  } while (num_ins < ins_length_max && (1 == has_bases || drand48() < opt->indel_extend));
                  assert(0 < num_ins);

                  j = muts_bed->muts[i].start;
                  if (1 == is_hom) {
                      ret[0]->s[j] = ret[1]->s[j] = (num_ins << ins_length_shift) | (ins << muttype_shift) | INSERT | c;
                  } else { // het-ins
                      ret[which_hap]->s[j] = (num_ins << ins_length_shift) | (ins << muttype_shift) | INSERT | c;
                  }
              }
          }
          else if (contig_i < muts_bed->muts[i].contig) {
              break;
          }
      }
  }
  
  // DEBUG
  for (i = 0; i != seq->l; ++i) {
      mut_t c[3];
      c[0] = nst_nt4_table[(mut_t)seq->s[i]];
      c[1] = hap1->s[i]; c[2] = hap2->s[i];
      if (c[0] >= 4) continue;
      if ((c[1] & mutmsk) != NOCHANGE || (c[2] & mutmsk) != NOCHANGE) {
          if (c[1] == c[2]) { // hom
              if ((c[1]&mutmsk) == SUBSTITUTE) { // substitution
                  continue;
              } else if ((c[1]&mutmsk) == DELETE) { // del
                  continue;
              } else if ((c[1] & mutmsk) == INSERT) { // ins
                  mut_t n = (c[1] >> ins_length_shift) & ins_length_mask;
                  assert(n > 0);
              }  else assert(0);
          } else { // het
              if ((c[1]&mutmsk) == SUBSTITUTE || (c[2]&mutmsk) == SUBSTITUTE) { // substitution
                  continue;
              } else if ((c[1]&mutmsk) == DELETE) {
                  continue;
              } else if ((c[2]&mutmsk) == DELETE) {
                  continue;
              } else if ((c[1]&mutmsk) == INSERT) { // ins 1
                  mut_t n = (c[1] >> ins_length_shift) & ins_length_mask;
                  assert(n > 0);
              } else if ((c[2]&mutmsk) == INSERT) { // ins 2
                  mut_t n = (c[2] >> ins_length_shift) & ins_length_mask;
                  assert(n > 0);
              } else assert(0);
          }
      }
  }
  // left-justify all the insertions and deletions
  int del_length;
  int prev_del[2] = {0, 0};
  for (i = 0; i != seq->l; ++i) {
      mut_t c[3];
      c[0] = nst_nt4_table[(mut_t)seq->s[i]];
      c[1] = hap1->s[i]; c[2] = hap2->s[i];
      if (c[0] >= 4) continue;
      if ((c[1] & mutmsk) != NOCHANGE || (c[2] & mutmsk) != NOCHANGE) {
          if (c[1] == c[2]) { // hom
              // TODO: code re-use
              if ((c[1]&mutmsk) == SUBSTITUTE) { // substitution
                  prev_del[0] = prev_del[1] = 0;
                  continue;
              } else if ((c[1]&mutmsk) == DELETE) { // del
                  if(prev_del[0] == 1 || prev_del[1] == 1) continue;
                  prev_del[0] = prev_del[1] = 1;
                  for(j=i+1,del_length=1;j<seq->l && (hap1->s[j]&mutmsk) == DELETE;j++) { // get del length
                      del_length++;
                  }
                  if(seq->l <= i+del_length) continue;
                  // left-justify
                  for(j=i-1;0<=j;j--) {
                      if(INSERT != (hap1->s[j]&mutmsk) && INSERT != (hap2->s[j]&mutmsk) // no insertion
                         && DELETE != (hap1->s[j]&mutmsk) && DELETE != (hap2->s[j]&mutmsk) // no deletion 
                         && (hap1->s[j]&3) == (hap1->s[j+del_length]&3)  // hap1 bases match
                         && (hap2->s[j]&3) == (hap2->s[j+del_length]&3)) { // hap2 bases match
                          // shift and make it NOCHANGE
                          mut_t t;
                          t = hap1->s[j]; hap1->s[j] = hap1->s[j+del_length]; hap1->s[j+del_length] = (t | mutmsk) ^ mutmsk;
                          t = hap2->s[j]; hap2->s[j] = hap2->s[j+del_length]; hap2->s[j+del_length] = (t | mutmsk) ^ mutmsk;
                      }
                      else {
                          break;
                      }
                  }
              } else if ((c[1] & mutmsk) == INSERT) { // ins
                  prev_del[0] = prev_del[1] = 0;
                  mut_t n = (c[1] >> ins_length_shift) & ins_length_mask, ins = (c[1] >> muttype_shift) & ins_mask;
                  assert(n > 0);
                  j=i;
                  while(0 < j
                        && INSERT != (hap1->s[j-1]&mutmsk) && INSERT != (hap2->s[j-1]&mutmsk) // no insertion
                        && DELETE != (hap1->s[j-1]&mutmsk) && DELETE != (hap2->s[j-1]&mutmsk) // no deletion 
                        && ((ins >> ((n-1) << 1)) & 3) == (hap1->s[j-1]&3)) { // end of insertion matches previous base
                      // update ins
                      ins = (ins | (3 << ((n-1) << 1))) ^ ((n-1) << 1); // zero out last base
                      ins <<= 2; // make room for the first base
                      ins |= (hap1->s[j-1]&3); // insert the first base
                      hap1->s[j] = hap2->s[j] = (hap1->s[j]&3); // make it NOCHANGE
                      j--;
                  }
                  hap1->s[j] = hap2->s[j] = (n << ins_length_shift) | (ins << muttype_shift) | INSERT | (hap1->s[j]&3); // re-insert
              }  else assert(0);
          } else { // het
              if ((c[1]&mutmsk) == SUBSTITUTE || (c[2]&mutmsk) == SUBSTITUTE) { // substitution
                  prev_del[0] = prev_del[1] = 0;
                  continue;
              } else if ((c[1]&mutmsk) == DELETE) {
                  if(prev_del[0] == 1) continue;
                  prev_del[0] = 1;
                  for(j=i+1,del_length=1;j<seq->l && (hap1->s[j]&mutmsk) == DELETE;j++) { // get del length
                      del_length++;
                  }
                  if(seq->l <= i+del_length) continue;
                  // left-justify
                  for(j=i-1;0<=j;j--) {
                      if(INSERT != (hap1->s[j]&mutmsk) // no insertion
                         && DELETE != (hap1->s[j]&mutmsk) // no deletion 
                         && (hap1->s[j]&3) == (hap1->s[j+del_length]&3))  { // hap1 bases match
                          // shift and make it NOCHANGE
                          mut_t t;
                          t = hap1->s[j]; hap1->s[j] = hap1->s[j+del_length]; hap1->s[j+del_length] = (t | mutmsk) ^ mutmsk;
                      }
                      else {
                          break;
                      }
                  }
              } else if ((c[2]&mutmsk) == DELETE) {
                  if(prev_del[1] == 1) continue;
                  prev_del[1] = 1;
                  for(j=i+1,del_length=1;j<seq->l && (hap2->s[j]&mutmsk) == DELETE;j++) { // get del length
                      del_length++;
                  }
                  if(seq->l <= i+del_length) continue;
                  // left-justify
                  for(j=i-1;0<=j;j--) {
                      if(INSERT != (hap2->s[j]&mutmsk) // no insertion
                         && DELETE != (hap2->s[j]&mutmsk) // no deletion 
                         && (hap2->s[j]&3) == (hap2->s[j+del_length]&3))  { // hap2 bases match
                          // shift and make it NOCHANGE
                          mut_t t;
                          t = hap2->s[j]; hap2->s[j] = hap2->s[j+del_length]; hap2->s[j+del_length] = (t | mutmsk) ^ mutmsk;
                      }
                      else {
                          break;
                      }
                  }
              } else if ((c[1]&mutmsk) == INSERT) { // ins 1
                  prev_del[0] = prev_del[1] = 0;
                  mut_t n = (c[1] >> ins_length_shift) & ins_length_mask, ins = (c[1] >> muttype_shift) & ins_mask;
                  assert(n > 0);
                  j=i;
                  while(0 < j
                        && INSERT != (hap1->s[j-1]&mutmsk) // no insertion
                        && DELETE != (hap1->s[j-1]&mutmsk) // no deletion 
                        && ((ins >> ((n-1) << 1)) & 3) == (hap1->s[j-1]&3)) { // end of insertion matches previous base
                      // update ins
                      ins = (ins | (3 << ((n-1) << 1))) ^ ((n-1) << 1); // zero out last base
                      ins <<= 2;
                      ins |= (hap1->s[j-1]&3); // insert the first base
                      hap1->s[j] = (hap1->s[j]&3); // make it NOCHANGE
                      j--;
                  }
                  hap1->s[j] = (n << ins_length_shift) | (ins << muttype_shift) | INSERT | (hap1->s[j]&3); // re-insert
              } else if ((c[2]&mutmsk) == INSERT) { // ins 2
                  prev_del[0] = prev_del[1] = 0;
                  mut_t n = (c[2] >> ins_length_shift) & ins_length_mask, ins = (c[2] >> muttype_shift) & ins_mask;
                  assert(n > 0);
                  j=i;
                  while(0 < j
                        && INSERT != (hap2->s[j-1]&mutmsk) // no insertion
                        && DELETE != (hap2->s[j-1]&mutmsk) // no deletion 
                        && ((ins >> ((n-1) << 1)) & 3) == (hap2->s[j-1]&3)) { // end of insertion matches previous base
                      // update ins
                      ins = (ins | (3 << ((n-1) << 1))) ^ ((n-1) << 1); // zero out last base
                      ins <<= 2;
                      ins |= (hap2->s[j-1]&3); // insert the first base
                      hap2->s[j] = (hap2->s[j]&3); // make it NOCHANGE
                      j--;
                  }
                  hap2->s[j] = (n << ins_length_shift) | (ins << muttype_shift) | INSERT | (hap2->s[j]&3); // re-insert
              } else assert(0);
          }
      }
      else {
          prev_del[0] = prev_del[1] = 0;
      }
  }
}

void mut_print(const char *name, const seq_t *seq, mutseq_t *hap1, mutseq_t *hap2, FILE *fpout)
{
  int32_t i, hap;
  for (i = 0; i != seq->l; ++i) {
      mut_t c[3];
      c[0] = nst_nt4_table[(int)seq->s[i]];
      c[1] = hap1->s[i]; c[2] = hap2->s[i];
      if (c[0] >= 4) continue;
      if ((c[1] & mutmsk) != NOCHANGE || (c[2] & mutmsk) != NOCHANGE) {
          fprintf(fpout, "%s\t%d\t", name, i+1);
          if (c[1] == c[2]) { // hom
              if ((c[1]&mutmsk) == SUBSTITUTE) { // substitution
                  fprintf(fpout, "%c\t%c\t3\n", "ACGTN"[c[0]], "ACGTN"[c[1]&0xf]);
              } else if ((c[1]&mutmsk) == DELETE) { // del
                  fprintf(fpout, "%c\t-\t3\n", "ACGTN"[c[0]]);
              } else if ((c[1] & mutmsk) == INSERT) { // ins
                  fprintf(fpout, "-\t");
                  mut_t n = (c[1] >> ins_length_shift) & ins_length_mask, ins = (c[1] >> muttype_shift) & ins_mask;
                  assert(n > 0);
                  while(n > 0) {
                      fputc("ACGTN"[ins & 0x3], fpout);
                      ins >>= 2;
                      n--;
                  }
                  fprintf(fpout, "\t3\n");
              }  else assert(0);
          } else { // het
              if ((c[1]&mutmsk) == SUBSTITUTE || (c[2]&mutmsk) == SUBSTITUTE) { // substitution
                  hap = ((c[1]&mutmsk) == SUBSTITUTE) ? 1 : 2;
                  fprintf(fpout, "%c\t%c\t%d\n", "ACGTN"[c[0]], "XACMGRSVTWYHKDBN"[1<<(c[1]&0x3)|1<<(c[2]&0x3)], hap);
              } else if ((c[1]&mutmsk) == DELETE) {
                  fprintf(fpout, "%c\t-\t1\n", "ACGTN"[c[0]]);
              } else if ((c[2]&mutmsk) == DELETE) {
                  fprintf(fpout, "%c\t-\t2\n", "ACGTN"[c[0]]);
              } else if ((c[1]&mutmsk) == INSERT) { // ins 1
                  fprintf(fpout, "-\t");
                  mut_t n = (c[1] >> ins_length_shift) & ins_length_mask, ins = (c[1] >> muttype_shift) & ins_mask;
                  assert(n > 0);
                  while (n > 0) {
                      fputc("ACGTN"[ins & 0x3], fpout);
                      ins >>= 2;
                      n--;
                  }
                  fprintf(fpout, "\t1\n");
              } else if ((c[2]&mutmsk) == INSERT) { // ins 2
                  fprintf(fpout, "-\t");
                  mut_t n = (c[2] >> ins_length_shift) & ins_length_mask, ins = (c[2] >> muttype_shift) & ins_mask;
                  assert(n > 0);
                  while (n > 0) {
                      fputc("ACGTN"[ins & 0x3], fpout);
                      ins >>= 2;
                      n--;
                  }
                  fprintf(fpout, "\t2\n");
              } else assert(0);
          }
      }
  }
}
