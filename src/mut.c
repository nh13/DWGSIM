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
mut_t mut_and_type_mask;
mut_t muttype_shift;
mut_t ins_length_shift;
mut_t ins_length_mask;
mut_t ins_length_max;
mut_t ins_long_length_max;
mut_t ins_mask;

void
mutseq_init_bounds()
{
  int32_t i;
  mutmsk = (mut_t)0x30;
  mut_and_type_mask = (mut_t)0x3F;
  muttype_shift = 6; // bits 5-6 store the mutation type
  ins_length_shift = 59; // bits 60-64 store the insertion length
  ins_length_mask = 0x1F; // bits 60-64 store the insertion length
  ins_length_max = ((ins_length_shift - muttype_shift) >> 1);
  ins_long_length_max = UINT32_MAX;
  if(ins_length_mask < ins_length_max) { // we exceed storing the length of the indel
      ins_length_max = ins_length_mask;
  }
  // create the insertion mask based on the maximum length insertion
  ins_mask = 0;
  for(i=0;i<ins_length_max;i++) {
      ins_mask = (ins_mask << 2) | 0x3; 
  }
}

mutseq_t *
mutseq_init()
{
  mutseq_t *seq = NULL;
  seq = calloc(1, sizeof(mutseq_t));
  mutseq_init_bounds();
  return seq;
}

void
mutseq_destroy(mutseq_t *seq)
{
  int32_t i;
  for(i=0;i<seq->ins_l;i++) {
      free(seq->ins[i]);
  }
  free(seq->ins);
  free(seq->s);
  free(seq);
}

mut_t
mut_get_ins_length(mutseq_t *seq, int32_t i)
{
  mut_t m, n, index;
  assert(0 <= i && i < seq->l);
  m = seq->s[i];
  assert(INSERT == (m & mutmsk));
  n = (m >> ins_length_shift) & ins_length_mask;
  if(0 == n) {
      uint32_t l = 0;
      index = (m >> muttype_shift) & ins_mask;
      assert(0 <= index && index < seq->ins_l);
      mut_get_ins_long_n(seq->ins[index], &l);
      return (mut_t)l;
  }
  else {
      return n;
  }
}

// returns 0 if it is a long insertion, 1 otherwise
// if 0 is returned, the index into the long insertion list
// is returned in "ins" and "n" is 0.
int32_t
mut_get_ins(mutseq_t *seq, int32_t i, mut_t *n, mut_t *ins)
{
  mut_t m;
  assert(0 <= i && i < seq->l);
  m = seq->s[i];
  assert(INSERT == (m & mutmsk));
  (*n) = (m >> ins_length_shift) & ins_length_mask;
  (*ins) = (m >> muttype_shift) & ins_mask;
  if(0 == (*n)) {
      assert(0 <= (*ins) && (*ins) < seq->ins_l);
      return 0;
  }
  else {
      assert(0 < (*n));
      assert((*n) <= ins_length_max);
      return 1;
  }
}

inline int32_t
mut_get_ins_bytes(int32_t n)
{
  int32_t num_bytes = 1; // one-byte for the type
  // variable number of bytes to store the length
  if(n <= UINT8_MAX) num_bytes += sizeof(uint8_t); 
  else if(n <= UINT16_MAX) num_bytes += sizeof(uint16_t); 
  else if(n <= UINT32_MAX) num_bytes += sizeof(uint32_t); 
  else {
      fprintf(stderr, "UINT32_MAX < n");
      exit(1);
  }
  // store the insertion
  num_bytes += mut_packed_len(n);
  return num_bytes;
}

uint8_t*
mut_get_ins_long_n(uint8_t *ins, uint32_t *n)
{
  switch(ins[0]) {
    case sizeof(uint8_t):
      (*n) = (uint32_t)(ins[1]);
      return ins + 2;
      break;
    case sizeof(uint16_t):
      (*n) = (uint32_t)(((uint16_t*)(ins+1))[0]);
      return ins + 3;
      break;
    case sizeof(uint32_t):
    default:
      (*n) = (uint32_t)(((uint32_t*)(ins+1))[0]);
      return ins + 5;
      break;
  }
  return NULL;
}

inline uint8_t*
mut_set_ins_long_n(uint8_t *ins, uint32_t n)
{
  // set type
  if(n <= UINT8_MAX) ins[0] = (uint8_t)sizeof(uint8_t);
  else if(n <= UINT16_MAX) ins[0] = (uint8_t)sizeof(uint16_t);
  else ins[0] = (uint8_t)sizeof(uint32_t);
      
  // set length
  switch(ins[0]) {
    case sizeof(uint8_t):
      ins[1] = (uint8_t)n;
      return ins + 2;
      break;
    case sizeof(uint16_t):
      (*((uint16_t*)(ins+1))) = (uint16_t)n;
      return ins + 3;
      break;
    case sizeof(uint32_t):
    default:
      (*((uint32_t*)(ins+1))) = n;
      return ins + 5;
      break;
  }
  return NULL;
}


static inline void
mut_print_ins(FILE *fp, mutseq_t *seq, int32_t i)
{
  mut_t n, ins;
  if(1 == mut_get_ins(seq, i, &n, &ins)) {
      while (n > 0) {
          fputc("ACGTN"[ins & 0x3], fp);
          ins >>= 2;
          n--;
      }
  }
  else { // long insertion
      int32_t byte_index, bit_index;
      uint32_t num_ins = 0;
      uint8_t *insertion = NULL;
      assert(NULL != seq->ins[ins]);
      insertion = mut_get_ins_long_n(seq->ins[ins], &num_ins);
      // reverse order
      byte_index = mut_packed_len(num_ins)-1;
      bit_index = (num_ins+3) & 3; // % 4
      while(0 < num_ins) {
          fputc("ACGTN"[(insertion[byte_index] >> (bit_index << 1)) & 0x3], fp);
          bit_index--;
          if(bit_index < 0) {
              bit_index = 3;
              byte_index--;
          }
          num_ins--;
      }
  }
}

// bases is NULL if we are to randomly simulate the bases
void mut_add_ins(dwgsim_opt_t *opt, mutseq_t *hap1, mutseq_t *hap2, int32_t i, int32_t c, int8_t hap, char *bases, mut_t num_ins)
{
  mut_t ins = 0;
  int64_t j;

  if (NULL == bases) {
      if(num_ins == 0) {
          // get the new insertion length
          do {
              num_ins++;
          } while (num_ins < ins_long_length_max && (num_ins < opt->indel_min || drand48() < opt->indel_extend));
      }
  } else {
      num_ins = strlen(bases); // ignores num_ins
  }
  if (ins_long_length_max < num_ins) num_ins = ins_long_length_max;

  if (hap < 0) {
      // set ploidy
      if (opt->is_hap || drand48() < 0.333333) { // hom-ins
          hap = 3;
      } else if (drand48() < 0.5) {
          hap = 1;
      } else {
          hap = 2;
      }
  }

  if (num_ins <= ins_length_max) { // short
      // generate the insertion
      if (NULL == bases) {
          for (j=0;j<num_ins;j++) {
              ins = (ins << 2) | (mut_t)(drand48() * 4.0);
          }
      } else {
          for (j = num_ins; 0 <= j; --j) {
              ins = (ins << 2) | nst_nt4_table[(int)bases[j]];
          } 
      }
      // store
      if (hap & 0x1) hap1->s[i] = (num_ins << ins_length_shift) | (ins << muttype_shift) | INSERT | c;
      if (hap & 0x2) hap2->s[i] = (num_ins << ins_length_shift) | (ins << muttype_shift) | INSERT | c;
  } else { // long
      int32_t byte_index, bit_index;
      uint8_t *hap1_ins = NULL, *hap2_ins = NULL;
      if (hap & 0x1) {
          while (hap1->ins_m <= hap1->ins_l) { // realloc
              hap1->ins_m = (hap1->ins_m < 16) ? 16 : (hap1->ins_m << 1); 
              hap1->ins = realloc(hap1->ins, sizeof(uint8_t*) * hap1->ins_m);
          }
          hap1->ins[hap1->ins_l] = calloc(mut_get_ins_bytes(num_ins), sizeof(uint8_t));
          hap1_ins = mut_set_ins_long_n(hap1->ins[hap1->ins_l], num_ins); 
      }
      if (hap & 0x2) {
          while (hap2->ins_m <= hap2->ins_l) { // realloc
              hap2->ins_m = (hap2->ins_m < 16) ? 16 : (hap2->ins_m << 1); 
              hap2->ins = realloc(hap2->ins, sizeof(uint8_t*) * hap2->ins_m);
          }
          hap2->ins[hap2->ins_l] = calloc(mut_get_ins_bytes(num_ins), sizeof(uint8_t));
          hap2_ins = mut_set_ins_long_n(hap2->ins[hap2->ins_l], num_ins); 
      }
      byte_index = 0;
      bit_index = 0;
      while(0 < num_ins) {
          uint8_t b;
          if (NULL == bases) {
              b = ((uint8_t)(drand48() * 4.0)) << (bit_index << 1);
          } else {
              b = nst_nt4_table[(int)bases[num_ins-1]] << (bit_index << 1);
          }
          if (hap & 0x1) hap1_ins[byte_index] |= b;
          if (hap & 0x2) hap2_ins[byte_index] |= b;
          bit_index++;
          if(4 == bit_index) {
              bit_index=0;
              byte_index++;
          }
          num_ins--;
      }
      if (hap & 0x1) {
          assert(hap1->ins_l <= ins_mask);
          hap1->s[i] = (hap1->ins_l << muttype_shift) | INSERT | c;
          hap1->ins_l++;
      }
      if (hap & 0x2) {
          assert(hap2->ins_l <= ins_mask);
          hap2->s[i] = (hap2->ins_l << muttype_shift) | INSERT | c;
          hap2->ins_l++;
      }
  }
}

void mut_debug(const seq_t *seq, mutseq_t *hap1, mutseq_t *hap2)
{
  int32_t i;
  // DEBUG
  for (i = 0; i < seq->l; ++i) {
      mut_t c[3];
      c[0] = nst_nt4_table[(mut_t)seq->s[i]];
      c[1] = hap1->s[i]; c[2] = hap2->s[i];
      if (c[0] >= 4) continue;
      if ((c[1] & mutmsk) != NOCHANGE || (c[2] & mutmsk) != NOCHANGE) {
          if ((c[1] & mut_and_type_mask) == (c[2] & mut_and_type_mask)) { // hom
              if ((c[1]&mutmsk) == SUBSTITUTE) { // substitution
                  assert((c[1]&0x3) == (c[2]&0x3));
                  assert((c[0]&0x3) != (c[1]&0x3));
                  continue;
              } else if ((c[1]&mutmsk) == DELETE) { // del
                  continue;
              } else if ((c[1] & mutmsk) == INSERT) { // ins
                  mut_t n1 = mut_get_ins_length(hap1, i);
                  mut_t n2 = mut_get_ins_length(hap1, i);
                  assert(n1 > 0);
                  assert(n1 <= ins_long_length_max);
                  assert(n1 == n2);
              } else assert(0);
          } else { // het
              if ((c[1]&mutmsk) == SUBSTITUTE || (c[2]&mutmsk) == SUBSTITUTE) { // substitution
                  assert((c[1]&0x3) != (c[2]&0x3));
                  assert((c[0]&0x3) == (c[1]&0x3) || (c[0]&0x3) == (c[2]&0x3));
                  assert((c[0]&0x3) != (c[1]&0x3) || (c[0]&0x3) != (c[2]&0x3));
                  continue;
              } else if ((c[1]&mutmsk) == DELETE) {
                  continue;
              } else if ((c[2]&mutmsk) == DELETE) {
                  continue;
              } else if ((c[1]&mutmsk) == INSERT) { // ins 1
                  mut_t n = mut_get_ins_length(hap1, i);
                  assert(n > 0);
                  assert(n <= ins_long_length_max);
              } else if ((c[2]&mutmsk) == INSERT) { // ins 2
                  mut_t n = mut_get_ins_length(hap2, i);
                  assert(n > 0);
                  assert(n <= ins_long_length_max);
              } else assert(0);
          }
      }
  }
}

static void
mut_left_justify_ins(mutseq_t *hap1, int32_t i)
{
  // NB: should we also be checking both haps for homozygous cases
  mut_t n, ins, j;
  if(1 == mut_get_ins(hap1, i, &n, &ins)) { // short
      assert(n > 0);
      j=i;
      while(0 < j
            && INSERT != (hap1->s[j-1]&mutmsk) // no insertion
            && SUBSTITUTE != (hap1->s[j-1]&mutmsk) // no substitution
            && DELETE != (hap1->s[j-1]&mutmsk) // no deletion 
            && ((ins >> ((n-1) << 1)) & 3) == (hap1->s[j-1]&3)) { // end of insertion matches previous base
          // update ins
          ins = (ins & ~((mut_t)3 << ((n-1) << 1))); // zero out the last base
          ins <<= 2; // shift over
          ins |= (hap1->s[j-1]&3); // insert the first base
          hap1->s[j] = (hap1->s[j]&3); // make it NOCHANGE
          j--;
      }
      hap1->s[j] = (n << ins_length_shift) | (ins << muttype_shift) | INSERT | (hap1->s[j]&3); // re-insert
  } else { // long
      int32_t bit_index, byte_index;
      uint32_t num_ins;
      uint8_t *insertion = NULL;
      insertion = mut_get_ins_long_n(hap1->ins[ins], &num_ins);
      assert(num_ins > 0);
      j=i;
      while(0 < j
            && INSERT != (hap1->s[j-1]&mutmsk) // no insertion
            && SUBSTITUTE != (hap1->s[j-1]&mutmsk) // no substitution
            && DELETE != (hap1->s[j-1]&mutmsk) // no deletion 
            // NB: last base of the insertion has a bit index of zero and byte index of zero
            && ((insertion[0] >> (0 << 1)) & 3) == (hap1->s[j-1]&3)) { // end of insertion matches previous base
          // start with the last base
          for (byte_index = 0; byte_index < mut_packed_len(num_ins); byte_index++) {
              insertion[byte_index] >>= 2; // shift down, get rid of first base (last two bits)
              if (byte_index + 1 < mut_packed_len(num_ins)) {
                  uint8_t b = (insertion[byte_index+1] >> (0 << 1)) & 3; // last base of the next byte (first two bits)
                  insertion[byte_index] |= b << 6; // move the new base to the last two bits
              }
          }
          // add a base to the front
          bit_index = (num_ins+3) & 3; // % 4
          byte_index = mut_packed_len(num_ins)-1;
          insertion[byte_index] |= (hap1->s[j-1]&3) << (bit_index << 1);
          hap1->s[j] = (hap1->s[j]&3); // make it NOCHANGE
          j--;
      }
      hap1->s[j] = (ins << muttype_shift) | INSERT | (hap1->s[j]&3); // re-insert
  }
}

// left-justify all the insertions and deletions
static void
mut_left_justify(const seq_t *seq, mutseq_t *hap1, mutseq_t *hap2)
{
  mut_t i, j;
  int del_length;
  int prev_del[2] = {0, 0};
  assert(hap1->l == seq->l);
  assert(hap2->l == seq->l);
  for (i = 0; i < seq->l; ++i) {
      mut_t c[3];
      c[0] = nst_nt4_table[(mut_t)seq->s[i]];
      c[1] = hap1->s[i]; c[2] = hap2->s[i];
      if (c[0] >= 4) continue;
      if ((c[1] & mutmsk) != NOCHANGE || (c[2] & mutmsk) != NOCHANGE) {
          if ((c[1] & mut_and_type_mask) == (c[2] & mut_and_type_mask)) { // hom
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
                  for(j=i-1;0 < i && 0<=j;j--) {
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
                      if(0 == j) break;
                  }
              } else if ((c[1] & mutmsk) == INSERT) { // ins
                  prev_del[0] = prev_del[1] = 0;
                  mut_left_justify_ins(hap1, i);
                  mut_left_justify_ins(hap2, i);
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
                  for(j=i-1;0 < i && 0<=j;j--) {
                      if(INSERT != (hap1->s[j]&mutmsk) // no insertion
                         && SUBSTITUTE != (hap1->s[j]&mutmsk) // no substitution
                         && DELETE != (hap1->s[j]&mutmsk) // no deletion 
                         && (hap1->s[j]&3) == (hap1->s[j+del_length]&3))  { // hap1 bases match
                          // shift and make it NOCHANGE
                          mut_t t;
                          t = hap1->s[j]; hap1->s[j] = hap1->s[j+del_length]; hap1->s[j+del_length] = (t | mutmsk) ^ mutmsk;
                      }
                      else {
                          break;
                      }
                      if(0 == j) break;
                  }
              } else if ((c[2]&mutmsk) == DELETE) {
                  if(prev_del[1] == 1) continue;
                  prev_del[1] = 1;
                  for(j=i+1,del_length=1;j<seq->l && (hap2->s[j]&mutmsk) == DELETE;j++) { // get del length
                      del_length++;
                  }
                  if(seq->l <= i+del_length) continue;
                  // left-justify
                  for(j=i-1;0 < i && 0<=j;j--) {
                      if(INSERT != (hap2->s[j]&mutmsk) // no insertion
                         && SUBSTITUTE != (hap2->s[j]&mutmsk) // no substitution
                         && DELETE != (hap2->s[j]&mutmsk) // no deletion 
                         && (hap2->s[j]&3) == (hap2->s[j+del_length]&3))  { // hap2 bases match
                          // shift and make it NOCHANGE
                          mut_t t;
                          t = hap2->s[j]; hap2->s[j] = hap2->s[j+del_length]; hap2->s[j+del_length] = (t | mutmsk) ^ mutmsk;
                      }
                      else {
                          break;
                      }
                      if(0 == j) break;
                  }
              } else if ((c[1]&mutmsk) == INSERT) { // ins 1
                  prev_del[0] = prev_del[1] = 0;
                  mut_left_justify_ins(hap1, i);
              } else if ((c[2]&mutmsk) == INSERT) { // ins 2
                  prev_del[0] = prev_del[1] = 0;
                  mut_left_justify_ins(hap2, i);
              } else assert(0);
          }
      }
      else {
          prev_del[0] = prev_del[1] = 0;
      }
  }
}

void mut_diref(dwgsim_opt_t *opt, const seq_t *seq, mutseq_t *hap1, mutseq_t *hap2, 
               int32_t contig_i, muts_input_t *muts_input)
{
  int32_t i, j, deleting = 0, deletion_length = 0;
  mutseq_t *ret[2];

  ret[0] = hap1; ret[1] = hap2;
  ret[0]->l = seq->l; ret[1]->l = seq->l;
  ret[0]->m = seq->m; ret[1]->m = seq->m;
  ret[0]->s = (mut_t *)calloc(seq->m, sizeof(mut_t));
  ret[1]->s = (mut_t *)calloc(seq->m, sizeof(mut_t));
  ret[0]->ins = NULL; ret[1]->ins = NULL;
  ret[0]->ins_l = 0; ret[1]->ins_l = 0;
  ret[0]->ins_m = 0; ret[1]->ins_m = 0;

  if(NULL == muts_input) {
      for (i = 0; i < seq->l; ++i) {
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
                      mut_add_ins(opt, ret[0], ret[1], i, c, -1, NULL, 0);
                  }
              }
          }
      }
  }
  else {
      if(MUT_INPUT_BED == muts_input->type) { // BED
          muts_bed_t *muts_bed = muts_input->data.bed; 
          // seed
          for (i = 0; i < seq->l; ++i) {
              ret[0]->s[i] = ret[1]->s[i] = (mut_t)nst_nt4_table[(int)seq->s[i]];
          }
          // mutates exactly based on a BED file
          for (i = 0; i < muts_bed->n; ++i) {
              if (muts_bed->muts[i].contig == contig_i) {
                  int32_t has_bases = 1, is_hom = 0, hap = 0;
                  int32_t which_hap = 0;
                  mut_t c;

                  // does this mutation have random bases?
                  if (0 == strcmp("*", muts_bed->muts[i].bases)) has_bases = 0; // random bases

                  // het or hom?
                  if (opt->is_hap || drand48() < 0.333333) {
                      is_hom = 1; // hom
                      hap = 3;
                  }
                  else {
                      which_hap = drand48()<0.5?0:1;
                      hap = 1 << which_hap;
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
                      c = (mut_t)nst_nt4_table[(int)seq->s[muts_bed->muts[i].start]];
                      if (0 == has_bases) {
                          mut_add_ins(opt, ret[0], ret[1], muts_bed->muts[i].start, c, hap, NULL, muts_bed->muts[i].end - muts_bed->muts[i].start);
                      } else {
                          mut_add_ins(opt, ret[0], ret[1], muts_bed->muts[i].start, c, hap, muts_bed->muts[i].bases, 0);
                      }
                  }
              }
              else if (contig_i < muts_bed->muts[i].contig) {
                  break;
              }
          }
      }
      else if(MUT_INPUT_TXT == muts_input->type || MUT_INPUT_VCF == muts_input->type) { // TXT or VCF
          muts_txt_t *muts_txt = NULL;
          if(MUT_INPUT_TXT == muts_input->type) {
              muts_txt = muts_input->data.txt; 
          }
          else {
              muts_txt = muts_input->data.vcf; 
          }
          // seed
          for (i = 0; i < seq->l; ++i) {
              ret[0]->s[i] = ret[1]->s[i] = (mut_t)nst_nt4_table[(int)seq->s[i]];
          }
          for (i = 0; i < muts_txt->n; ++i) {
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
                      mut_add_ins(opt, ret[0], ret[1], pos-1, c, is_hap, muts_txt->muts[i].bases, 0);
                  }
              }
          }
      }
      else {
          fprintf(stderr, "Error: unknown mutation input type!\n");
          exit(1);
      }
  }

  // DEBUG
  mut_debug(seq, hap1, hap2);
  
  // DEBUG
  mut_left_justify(seq, hap1, hap2);
  mut_debug(seq, hap1, hap2);
}

void mut_print_header_pre(FILE *fpout_vcf)
{
  fprintf(fpout_vcf, "##fileformat=VCFv4.1\n");
}

void mut_print_header_post(FILE *fpout_vcf)
{
  fprintf(fpout_vcf, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n");
  fprintf(fpout_vcf, "##INFO=<ID=pl,Number=1,Type=Integer,Description=\"Phasing: 1 - HET contig 1, #2 - HET contig #2, 3 - HOM both contigs\">\n"); 
  fprintf(fpout_vcf, "##INFO=<ID=mt,Number=1,Type=String,Description=\"Variant Type: SUBSTITUTE/INSERT/DELETE\">\n");
  fprintf(fpout_vcf, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
}

void mut_print_header_contig(FILE *fpout_vcf, const char *name, int32_t length)
{
  // TODO: does mutating the reference change the sequence length?  If so, this
  // is not correct.
  fprintf(fpout_vcf, "##contig=<ID=%s,length=%d>\n", name, length);
}

// TODO: add AD, DP, GQ, GT, PL
void mut_print(const char *name, const seq_t *seq, mutseq_t *hap1, mutseq_t *hap2, FILE *fpout_txt, FILE *fpout_vcf)
{
  int32_t i, j, hap;
  
  // body
  mut_t mut_prev[2] = {NOCHANGE,NOCHANGE}; // for deletions
  for (i = 0; i < seq->l; ++i) {
      mut_t c[3];
      c[0] = nst_nt4_table[(int)seq->s[i]];
      c[1] = hap1->s[i]; c[2] = hap2->s[i];
      if (c[0] >= 4) {
          // do nothing
      } else if ((c[1] & mutmsk) != NOCHANGE || (c[2] & mutmsk) != NOCHANGE) {
          fprintf(fpout_txt, "%s\t%d\t", name, i+1);
          if ((c[1] & mut_and_type_mask) == (c[2] & mut_and_type_mask)) { // hom
              if ((c[1]&mutmsk) == SUBSTITUTE) { // substitution
                  fprintf(fpout_txt, "%c\t%c\t3\n", "ACGTN"[c[0]], "ACGTN"[c[1]&0xf]);
                  fprintf(fpout_vcf, "%s\t%d\t.\t%c\t%c\t100\tPASS\tAF=1.0;pl=3;mt=SUBSTITUTE\n", name, i+1, "ACGTN"[c[0]], "ACGTN"[c[1]&0xf]);
              } else if ((c[1]&mutmsk) == DELETE) { // del
                  fprintf(fpout_txt, "%c\t-\t3\n", "ACGTN"[c[0]]);
                  if(0 == mut_prev[0] || 0 == mut_prev[1]) {
                      fprintf(fpout_vcf, "%s\t%d\t.\t", name, i);
                      if (0 < i) fputc("ACGTN"[nst_nt4_table[(int)seq->s[i-1]]], fpout_vcf);
                      for (j = i; j < seq->l && (c[1] & mut_and_type_mask) == (c[2] & mut_and_type_mask) && (c[1]&mutmsk) == DELETE; ++j) {
                          fputc("ACGTN"[c[0]], fpout_vcf);
                          // NB: this modifies 'c'
                          if (j+1 < seq->l) { 
                              c[0] = nst_nt4_table[(int)seq->s[j+1]];
                              c[1] = hap1->s[j+1]; c[2] = hap2->s[j+1];
                          }
                      }
                      if (0 < i) fprintf(fpout_vcf, "\t%c", "ACGTN"[nst_nt4_table[(int)seq->s[i-1]]]);
                      else fprintf(fpout_vcf, "\t.");
                      fprintf(fpout_vcf, "\t100\tPASS\tAF=1.0;pl=3;mt=DELETE\n"); 
                  }
                  // NB: convert back 'c'
                  c[0] = nst_nt4_table[(int)seq->s[i]];
                  c[1] = hap1->s[i]; c[2] = hap2->s[i];
              } else if ((c[1] & mutmsk) == INSERT) { // ins
                  fprintf(fpout_txt, "-\t");
                  mut_print_ins(fpout_txt, hap1, i);
                  fprintf(fpout_txt, "\t3\n");
                  fprintf(fpout_vcf, "%s\t%d\t.\t%c\t%c", name, i+1, "ACGTN"[c[0]], "ACGTN"[c[0]]);
                  mut_print_ins(fpout_vcf, hap1, i);
                  fprintf(fpout_vcf, "\t100\tPASS\tAF=1.0;pl=3;mt=INSERT\n");
              }  else assert(0);
          } else { // het
              if ((c[1]&mutmsk) == SUBSTITUTE || (c[2]&mutmsk) == SUBSTITUTE) { // substitution
                  hap = ((c[1]&mutmsk) == SUBSTITUTE) ? 1 : 2;
                  fprintf(fpout_txt, "%c\t%c\t%d\n", "ACGTN"[c[0]], "XACMGRSVTWYHKDBN"[1<<(c[1]&0x3)|1<<(c[2]&0x3)], hap);
                  if(1 == hap) fprintf(fpout_vcf, "%s\t%d\t.\t%c\t%c\t100\tPASS\tAF=0.5;pl=1;mt=SUBSTITUTE\n", name, i+1, "ACGTN"[c[0]], "ACGTN"[c[1]&0xf]);
                  else fprintf(fpout_vcf, "%s\t%d\t.\t%c\t%c\t100\tPASS\tAF=0.5;pl=2;mt=SUBSTITUTE\n", name, i+1, "ACGTN"[c[0]], "ACGTN"[c[2]&0xf]);
              } else if ((c[1]&mutmsk) == DELETE) {
                  fprintf(fpout_txt, "%c\t-\t1\n", "ACGTN"[c[0]]);
                  if(0 == mut_prev[0]) {
                      fprintf(fpout_vcf, "%s\t%d\t.\t", name, i);
                      if (0 < i) fputc("ACGTN"[nst_nt4_table[(int)seq->s[i-1]]], fpout_vcf);
                      for (j = i; j < seq->l && (c[1] & mut_and_type_mask) != (c[2] & mut_and_type_mask) && (c[1]&mutmsk) == DELETE; ++j) {
                          fputc("ACGTN"[c[0]], fpout_vcf);
                          // NB: this modifies 'c'
                          if (j+1 < seq->l) { 
                              c[0] = nst_nt4_table[(int)seq->s[j+1]];
                              c[1] = hap1->s[j+1]; c[2] = hap2->s[j+1];
                          }
                      }
                      if (0 < i) fprintf(fpout_vcf, "\t%c", "ACGTN"[nst_nt4_table[(int)seq->s[i-1]]]);
                      else fprintf(fpout_vcf, "\t.");
                      fprintf(fpout_vcf, "\t100\tPASS\tAF=0.5;pl=1;mt=DELETE\n"); 
                  }
                  // NB: convert back 'c'
                  c[0] = nst_nt4_table[(int)seq->s[i]];
                  c[1] = hap1->s[i]; c[2] = hap2->s[i];
              } else if ((c[2]&mutmsk) == DELETE) {
                  fprintf(fpout_txt, "%c\t-\t2\n", "ACGTN"[c[0]]);
                  if(0 == mut_prev[1]) {
                      fprintf(fpout_vcf, "%s\t%d\t.\t", name, i);
                      if (0 < i) fputc("ACGTN"[nst_nt4_table[(int)seq->s[i-1]]], fpout_vcf);
                      for (j = i; j < seq->l && (c[1] & mut_and_type_mask) != (c[2] & mut_and_type_mask) && (c[2]&mutmsk) == DELETE; ++j) {
                          fputc("ACGTN"[c[0]], fpout_vcf);
                          // NB: this modifies 'c'
                          if (j+1 < seq->l) { 
                              c[0] = nst_nt4_table[(int)seq->s[j+1]];
                              c[1] = hap1->s[j+1]; c[2] = hap2->s[j+1];
                          }
                      }
                      if (0 < i) fprintf(fpout_vcf, "\t%c", "ACGTN"[nst_nt4_table[(int)seq->s[i-1]]]);
                      else fprintf(fpout_vcf, "\t.");
                      fprintf(fpout_vcf, "\t100\tPASS\tAF=0.5;pl=2;mt=DELETE\n"); 
                  }
                  // NB: convert back 'c'
                  c[0] = nst_nt4_table[(int)seq->s[i]];
                  c[1] = hap1->s[i]; c[2] = hap2->s[i];
              } else if ((c[1]&mutmsk) == INSERT) { // ins 1
                  fprintf(fpout_txt, "-\t");
                  mut_print_ins(fpout_txt, hap1, i);
                  fprintf(fpout_txt, "\t1\n");
                  fprintf(fpout_vcf, "%s\t%d\t.\t%c\t%c", name, i+1, "ACGTN"[c[0]], "ACGTN"[c[0]]);
                  mut_print_ins(fpout_vcf, hap1, i);
                  fprintf(fpout_vcf, "\t100\tPASS\tAF=0.5;pl=1;mt=INSERT\n");
              } else if ((c[2]&mutmsk) == INSERT) { // ins 2
                  fprintf(fpout_txt, "-\t");
                  mut_print_ins(fpout_txt, hap2, i);
                  fprintf(fpout_txt, "\t2\n");
                  fprintf(fpout_vcf, "%s\t%d\t.\t%c\t%c", name, i+1, "ACGTN"[c[0]], "ACGTN"[c[0]]);
                  mut_print_ins(fpout_vcf, hap2, i);
                  fprintf(fpout_vcf, "\t100\tPASS\tAF=0.5;pl=2;mt=INSERT\n");
              } else assert(0);
          }
      }
      mut_prev[0] = (c[1] & mutmsk) != NOCHANGE;
      mut_prev[1] = (c[2] & mutmsk) != NOCHANGE;
  }
}
