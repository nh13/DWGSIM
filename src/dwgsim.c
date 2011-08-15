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

/* Original Contact: Heng Li <lh3@sanger.ac.uk> */

/* This program is separated from maq's read simulator with Colin
 * Hercus' modification to allow longer indels. Colin is the chief
 * developer of novoalign. */

/* This program was modified by Nils Homer for inclusion in the 
 * DNAA package */

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <stdio.h>
#include <unistd.h>
#include <stdint.h>
#include <ctype.h>
#include <string.h>
//#include <config.h>

#define __gf_add(_x, _y) ((_x >= 4 || _y >= 4) ? 4 : (_x ^ _y))

const uint8_t nst_nt4_table[256] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

#define __gen_read(x, start, iter) do {									\
    for (i = (start), k = 0, ext_coor[x] = -10; i >= 0 && i < seq.l && k < s[x]; iter) {	\
        mut_t c = target[i], mut_type = c & mutmsk;			\
        if (ext_coor[x] < 0) {								\
            if (mut_type != NOCHANGE && mut_type != SUBSTITUTE) continue; \
            ext_coor[x] = i;								\
            if(1 == strand[x]) ext_coor[x] -= s[x]-1; \
        }													\
        if (mut_type == DELETE) { \
            ++n_indel[x];				\
            if(1 == strand[x]) ext_coor[x]--; \
            if(0 == k) n_indel_first[x]++; \
        } \
        else if (mut_type == NOCHANGE || mut_type == SUBSTITUTE) { \
            tmp_seq[x][k++] = c & 0xf;						\
            if (mut_type == SUBSTITUTE) { \
                ++n_sub[x];			\
                if(0 == k) n_sub_first[x]++; \
            } 												\
        } else {											\
            mut_t n, ins;										\
            assert(mut_type == INSERT); \
            ++n_indel[x];									\
            n_indel_first[x]++;							\
            if(0 == strand[x]) { \
                n = (c>>ins_length_shift) & ins_length_mask; \
                ins = c>>muttype_shift; \
                while(n > 0 && k < s[x]) { \
                    tmp_seq[x][k++] = ins & 0x3;                \
                    --n, ins >>= 2; \
                } \
                if(k < s[x]) tmp_seq[x][k++] = c & 0xf;						\
            } else { \
                tmp_seq[x][k++] = c & 0xf;						\
                n = (c>>ins_length_shift) & ins_length_mask; \
                ins = c>>muttype_shift; \
                while(n > 0 && k < s[x]) { \
                    ext_coor[x]++; \
                    tmp_seq[x][k++] = (ins >> ((n-1) << 1) & 0x3);                \
                    --n; \
                } \
            } \
        }													\
    }														\
    if (k != s[x]) ext_coor[x] = -10;						\
    if (1 == strand[x]) { \
        for (k = 0; k < s[x]; ++k) tmp_seq[x][k] = tmp_seq[x][k] < 4? 3 - tmp_seq[x][k] : 4; \
    } 														\
} while (0)

/* Simple normal random number generator, copied from genran.c */

double ran_normal()
{ 
  static int iset = 0; 
  static double gset; 
  double fac, rsq, v1, v2; 
  if (iset == 0) {
      do { 
          v1 = 2.0 * drand48() - 1.0;
          v2 = 2.0 * drand48() - 1.0; 
          rsq = v1 * v1 + v2 * v2;
      } while (rsq >= 1.0 || rsq == 0.0);
      fac = sqrt(-2.0 * log(rsq) / rsq); 
      gset = v1 * fac; 
      iset = 1;
      return v2 * fac;
  } else {
      iset = 0;
      return gset;
  }
}

/* FASTA parser, copied from seq.c */

typedef struct {
    int l, m; /* length and maximum buffer size */
    unsigned char *s; /* sequence */
} seq_t;

#define INIT_SEQ(seq) (seq).s = 0; (seq).l = (seq).m = 0

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

/* Error-checking open, copied from utils.c */

#define xopen(fn, mode) err_xopen_core(__func__, fn, mode)

FILE *err_xopen_core(const char *func, const char *fn, const char *mode)
{
  FILE *fp = 0;
  if (strcmp(fn, "-") == 0)
    return (strstr(mode, "r"))? stdin : stdout;
  if ((fp = fopen(fn, mode)) == 0) {
      fprintf(stderr, "[%s] fail to open file '%s'. Abort!\n", func, fn);
      abort();
  }
  return fp;
}

/* dwgsim */

#define ERR_RATE 0.02

typedef struct {
    double start, by, end;
} error_t;

enum muttype_t {NOCHANGE = 0, INSERT = 0x10, SUBSTITUTE = 0x20, DELETE = 0x30};
typedef uint64_t mut_t;
static mut_t mutmsk = (mut_t)0x30;

typedef struct {
    int l, m; /* length and maximum buffer size */
    mut_t *s; /* sequence */
} mutseq_t;

//static mut_t base_shift = 4; // lower 4-bits store the base
static mut_t muttype_shift = 6; // bits 5-6 store the mutation type
// bits 7-28 store the insertion
static mut_t ins_length_shift = 60; // bits 61-64 store the insertion length
static mut_t ins_length_mask = 0xF; // bits 61-64 store the insertion length
static double MUT_RATE = 0.001;
static double INDEL_FRAC = 0.1;
static double INDEL_EXTEND = 0.3;
static double RAND_READ = 0.1;
enum {
    ILLUMINA=0,
    SOLID=1,
    IONTORRENT=2
};
static int DATA_TYPE = ILLUMINA;
static int MAX_N = 0;
static uint8_t FLOW_ORDER[1024]="TACG";
static int FLOW_ORDER_LEN = 4;

#define __gen_errors_mismatches(_cur_seq, _j, _start, _iter, _len) do { \
    for (i = (_start); 0 <= i && i < _len; _iter) { \
        mut_t c = _cur_seq[_j][i]; \
        if (c >= 4) c = 4; \
        else if(drand48() < e[_j]->start + e[_j]->by*i) { \
            c = (c + (mut_t)(drand48() * 3.0 + 1)) & 3; \
            ++n_err[_j]; \
            if(0 == i) ++n_err_first[_j]; \
        } \
        _cur_seq[_j][i] = c; \
    } \
} while(0)

static int32_t
generate_errors_flows(uint8_t **seq, int32_t *mem, int32_t len, uint8_t strand, double e, int32_t *_n_err)
{
  int32_t i, j, k, hp_l, flow_i, n_err;
  uint8_t prev_c, c;

  // remove Ns
  for(i=0;i<len;i++) {
      if(4 <= (*seq)[i]) {
          (*seq)[i] = 0;
      }
  }

  if(1 == strand) {
      for(i=0;i<len>>1;i++) {
          c = (*seq)[i];
          (*seq)[i] = (*seq)[len-i-1];
          (*seq)[len-i-1] = c;
      }
  }

  // skip forward to the first flow
  for(i=0;i<FLOW_ORDER_LEN;i++) {
      c = (4 <= (*seq)[0]) ? 0 : (*seq)[0];
      if(c == FLOW_ORDER[i]) {
          break;
      }
  }
  assert(FLOW_ORDER_LEN != i);
  flow_i = i;
  // first pass: add errors in the current sequence (non-empty flows)
  prev_c = 4;
  for(i=0;i<len;i++) {
      c = (4 <= (*seq)[i]) ? 0 : (*seq)[i];
      while(c != FLOW_ORDER[flow_i]) { // skip paste empty flows
          flow_i = (flow_i + 1) % FLOW_ORDER_LEN;
      }
      if(prev_c != c) { // new hp
          n_err = 0;
          while(drand48() < e) { // how many bases should we insert/delete
              n_err++;
          }
          if(0 < n_err) {
              if(drand48() < 0.5) { // insert
                  // more memory
                  while((*mem) <= len + n_err) {
                      (*mem) <<= 1; // double
                      (*seq) = realloc((*seq), sizeof(uint8_t) * (*mem));
                  }
                  // shift up
                  for(j=len-1;i<=j;j--) {
                      (*seq)[j+n_err] = (*seq)[j];
                  }
                for(j=i;j<i+n_err;j++) {
                      (*seq)[j] = c;
                  }
                  len += n_err;
              }
              else {
                  int32_t next_c = 4;
                  // get the hp length
                  for(j=i,hp_l=0;j<len;j++,hp_l++) {
                      next_c = (4 <= (*seq)[j]) ? 0 : (*seq)[j];
                      if(c != next_c) {
                          break;
                      }
                  }
                  // bound errors by the hp length
                  n_err = (hp_l < n_err) ? hp_l : n_err;
                  // shift down
                  for(j=i;j<len-n_err;j++) {
                      (*seq)[j] = (*seq)[j+n_err];
                  }
                  len -= n_err;
                  // Note: we need to make sure that if we delete all the bases,
                  // that the neighboring bases are not the same.  If they are,
                  // we need to simulate a "dot-fill", whereby we add a
                  // non-incorporating flow.
                  if(n_err == hp_l && (0 == i || prev_c == next_c)) {
                      j = 0;
                      // get the number of subsequent flows until the next base
                      while(next_c != FLOW_ORDER[(flow_i + j) % FLOW_ORDER_LEN]) { // skip paste empty flows
                          j++;
                      }
                      assert(0 < j);
                      // pick one to fill in
                      k = (int)(drand48() * j);
                      // shift up
                      for(j=len-1;i<=j;j--) {
                          (*seq)[j+1] = (*seq)[j];
                      }
                      // add the filled in base
                      (*seq)[i] = FLOW_ORDER[(flow_i + k) % FLOW_ORDER_LEN];
                      len++;
                  }
              }
              (*_n_err) += n_err;
          }
          prev_c = c;
      }
  }

  // second pass: add empty flows
  prev_c = 4;
  for(i=0;i<len;i++) {
      c = (4 <= (*seq)[i]) ? 0 : (*seq)[i];
      while(c != FLOW_ORDER[flow_i]) {
          n_err = 0;
          while(drand48() < e) {
              n_err++;
          }
          if(0 < n_err) {  // insert
              // more memory
              while((*mem) <= len + n_err) {
                  (*mem) <<= 1; // double
                  (*seq) = realloc((*seq), sizeof(uint8_t) * (*mem));
              }
              // shift up
              for(j=len-1;i<=j;j--) {
                  (*seq)[j+n_err] = (*seq)[j];
              }
              // copy
              for(j=i;j<i+n_err;j++) {
                  (*seq)[j] = FLOW_ORDER[flow_i];
              }
              len += n_err;
              (*_n_err) += n_err;
          }
          flow_i = (flow_i + 1) % FLOW_ORDER_LEN;
      }
  }

  if(1 == strand) {
      for(i=0;i<len>>1;i++) {
          c = (*seq)[i];
          (*seq)[i] = (*seq)[len-i-1];
          (*seq)[len-i-1] = c;
      }
  }

  return len;
}

void maq_mut_diref(const seq_t *seq, int is_hap, mutseq_t *hap1, mutseq_t *hap2)
{
  int i, deleting = 0;
  mutseq_t *ret[2];

  ret[0] = hap1; ret[1] = hap2;
  ret[0]->l = seq->l; ret[1]->l = seq->l;
  ret[0]->m = seq->m; ret[1]->m = seq->m;
  ret[0]->s = (mut_t *)calloc(seq->m, sizeof(mut_t));
  ret[1]->s = (mut_t *)calloc(seq->m, sizeof(mut_t));
  for (i = 0; i != seq->l; ++i) {
      mut_t c;
      c = ret[0]->s[i] = ret[1]->s[i] = (mut_t)nst_nt4_table[(int)seq->s[i]];
      if (deleting) {
          if (drand48() < INDEL_EXTEND) {
              if (deleting & 1) ret[0]->s[i] |= DELETE|c;
              if (deleting & 2) ret[1]->s[i] |= DELETE|c;
              continue;
          } else deleting = 0;
      }
      if (c < 4 && drand48() < MUT_RATE) { // mutation
          if (drand48() >= INDEL_FRAC) { // substitution
              double r = drand48();
              c = (c + (mut_t)(r * 3.0 + 1)) & 3;
              if (is_hap || drand48() < 0.333333) { // hom
                  ret[0]->s[i] = ret[1]->s[i] = SUBSTITUTE|c;
              } else { // het
                  ret[drand48()<0.5?0:1]->s[i] = SUBSTITUTE|c;
              }
          } else { // indel
              if (drand48() < 0.5) { // deletion
                  if (is_hap || drand48() < 0.333333) { // hom-del
                      ret[0]->s[i] = ret[1]->s[i] = DELETE|c;
                      deleting = 3;
                  } else { // het-del
                      deleting = drand48()<0.5?1:2;
                      ret[deleting-1]->s[i] = DELETE|c;
                  }
              } else { // insertion
                  mut_t num_ins = 0, ins = 0;
                  do {
                      num_ins++;
                      ins = (ins << 2) | (mut_t)(drand48() * 4.0);
                  } while (num_ins < ((ins_length_shift - muttype_shift) >> 1) && drand48() < INDEL_EXTEND);
                  assert(0 < num_ins);

                  if (is_hap || drand48() < 0.333333) { // hom-ins
                      ret[0]->s[i] = ret[1]->s[i] = (num_ins << ins_length_shift) | (ins << muttype_shift) | INSERT | c;
                  } else { // het-ins
                      ret[drand48()<0.5?0:1]->s[i] = (num_ins << ins_length_shift) | (ins << muttype_shift) | INSERT | c;
                  }
              }
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
  int j, del_length;
  int prev_del[2] = {0, 0};
  for (i = 0; i != seq->l; ++i) {
      mut_t c[3];
      c[0] = nst_nt4_table[(mut_t)seq->s[i]];
      c[1] = hap1->s[i]; c[2] = hap2->s[i];
      if (c[0] >= 4) continue;
      if ((c[1] & mutmsk) != NOCHANGE || (c[2] & mutmsk) != NOCHANGE) {
          if (c[1] == c[2]) { // hom
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
                  mut_t n = (c[1] >> ins_length_shift) & ins_length_mask, ins = c[1] >> muttype_shift;
                  assert(n > 0);
                  j=i;
                  while(0 < j
                        && INSERT != (hap1->s[j-1]&mutmsk) && INSERT != (hap2->s[j-1]&mutmsk) // no insertion
                        && DELETE != (hap1->s[j-1]&mutmsk) && DELETE != (hap2->s[j-1]&mutmsk) // no deletion 
                        && ((ins >> ((n-1) << 1)) & 3) == (hap1->s[j-1]&3)) { // end of insertion matches previous base
                      // update ins
                      ins = (ins | (3 << ((n-1) << 1))) ^ ((n-1) << 1); // zero out last base
                      ins <<= 2;
                      ins |= (hap1->s[j-1]&mutmsk);
                      hap1->s[j] = hap2->s[j] = (hap1->s[j] | mutmsk) ^ mutmsk; // make it NOCHANGE
                      j--;
                  }
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
                  mut_t n = (c[1] >> ins_length_shift) & ins_length_mask, ins = c[1] >> muttype_shift;
                  assert(n > 0);
                  j=i;
                  while(0 < j
                        && INSERT != (hap1->s[j-1]&mutmsk) // no insertion
                        && DELETE != (hap1->s[j-1]&mutmsk) // no deletion 
                        && ((ins >> ((n-1) << 1)) & 3) == (hap1->s[j-1]&3)) { // end of insertion matches previous base
                      // update ins
                      ins = (ins | (3 << ((n-1) << 1))) ^ ((n-1) << 1); // zero out last base
                      ins <<= 2;
                      ins |= (hap1->s[j-1]&mutmsk);
                      hap1->s[j] = (hap1->s[j] | mutmsk) ^ mutmsk; // make it NOCHANGE
                      j--;
                  }
              } else if ((c[2]&mutmsk) == INSERT) { // ins 2
                  prev_del[0] = prev_del[1] = 0;
                  mut_t n = (c[2] >> ins_length_shift) & ins_length_mask, ins = c[2] >> muttype_shift;
                  assert(n > 0);
                  j=i;
                  while(0 < j
                        && INSERT != (hap2->s[j-1]&mutmsk) // no insertion
                        && DELETE != (hap2->s[j-1]&mutmsk) // no deletion 
                        && ((ins >> ((n-1) << 1)) & 3) == (hap2->s[j-1]&3)) { // end of insertion matches previous base
                      // update ins
                      ins = (ins | (3 << ((n-1) << 1))) ^ ((n-1) << 1); // zero out last base
                      ins <<= 2;
                      ins |= (hap2->s[j-1]&mutmsk);
                      hap2->s[j] = (hap2->s[j] | mutmsk) ^ mutmsk; // make it NOCHANGE
                      j--;
                  }
              } else assert(0);
          }
      }
      else {
          prev_del[0] = prev_del[1] = 0;
      }
  }
}

// Columns:
// 1 - chromosome name
// 2 - position (one-based)
// 3 - reference base (dash if there was an insertion)
// 4 - variant allele (IUPAC code or insertion base(s))
// 5 - '-' for homozygous, '+' for heterozygous
void maq_print_mutref(const char *name, const seq_t *seq, mutseq_t *hap1, mutseq_t *hap2, FILE *fpout)
{
  int i;
  for (i = 0; i != seq->l; ++i) {
      mut_t c[3];
      c[0] = nst_nt4_table[(int)seq->s[i]];
      c[1] = hap1->s[i]; c[2] = hap2->s[i];
      if (c[0] >= 4) continue;
      if ((c[1] & mutmsk) != NOCHANGE || (c[2] & mutmsk) != NOCHANGE) {
          fprintf(fpout, "%s\t%d\t", name, i+1);
          if (c[1] == c[2]) { // hom
              if ((c[1]&mutmsk) == SUBSTITUTE) { // substitution
                  fprintf(fpout, "%c\t%c\t-\n", "ACGTN"[c[0]], "ACGTN"[c[1]&0xf]);
              } else if ((c[1]&mutmsk) == DELETE) { // del
                  fprintf(fpout, "%c\t-\t-\n", "ACGTN"[c[0]]);
              } else if ((c[1] & mutmsk) == INSERT) { // ins
                  fprintf(fpout, "-\t");
                  mut_t n = (c[1] >> ins_length_shift) & ins_length_mask, ins = c[1] >> muttype_shift;
                  assert(n > 0);
                  while(n > 0) {
                      fputc("ACGTN"[ins & 0x3], fpout);
                      ins >>= 2;
                      n--;
                  }
                  fprintf(fpout, "\t-\n");
              }  else assert(0);
          } else { // het
              if ((c[1]&mutmsk) == SUBSTITUTE || (c[2]&mutmsk) == SUBSTITUTE) { // substitution
                  fprintf(fpout, "%c\t%c\t+\n", "ACGTN"[c[0]], "XACMGRSVTWYHKDBN"[1<<(c[1]&0x3)|1<<(c[2]&0x3)]);
              } else if ((c[1]&mutmsk) == DELETE) {
                  fprintf(fpout, "%c\t-\t+\n", "ACGTN"[c[0]]);
              } else if ((c[2]&mutmsk) == DELETE) {
                  fprintf(fpout, "%c\t-\t+\n", "ACGTN"[c[0]]);
              } else if ((c[1]&mutmsk) == INSERT) { // ins 1
                  fprintf(fpout, "-\t");
                  mut_t n = (c[1] >> ins_length_shift) & ins_length_mask, ins = c[1] >> muttype_shift;
                  assert(n > 0);
                  while (n > 0) {
                      fputc("ACGTN"[ins & 0x3], fpout);
                      ins >>= 2;
                      n--;
                  }
                  fprintf(fpout, "\t+\n");
              } else if ((c[2]&mutmsk) == INSERT) { // ins 2
                  fprintf(fpout, "-\t");
                  mut_t n = (c[2] >> ins_length_shift) & ins_length_mask, ins = c[2] >> muttype_shift;
                  assert(n > 0);
                  while (n > 0) {
                      fputc("ACGTN"[ins & 0x3], fpout);
                      ins >>= 2;
                      n--;
                  }
                  fprintf(fpout, "\t+\n");
              } else assert(0);
          }
      }
  }
}

void dwgsim_core(FILE *fpout0, FILE *fpout1, FILE *fpout2, FILE *fpout3, FILE *fp_fa, FILE *fp_fai, 
                 error_t *e1, error_t *e2, int is_hap, uint64_t N, int dist, int std_dev, 
                 int size_l, int size_r, int max_n, int strandedness)
{
  seq_t seq;
  mutseq_t rseq[2];
  uint64_t tot_len, ii=0, ctr=0;
  int i, l, n_ref;
  char name[256], *qstr;
  int size[2], prev_skip=0, qstr_l=0;
  int num_n[2];
  uint8_t *tmp_seq[2]={NULL,NULL};
  int32_t tmp_seq_mem[2]={0,0};
  uint64_t n_sim = 0;
  mut_t *target;
  error_t *e[2]={NULL,NULL};

  e[0] = e1; e[1] = e2;

  INIT_SEQ(seq);
  srand48(time(0));
  seq_set_block_size(0x1000000);
  l = size_l > size_r? size_l : size_r;
  qstr_l = l;
  qstr = (char*)calloc(qstr_l+1, 1);
  tmp_seq[0] = (uint8_t*)calloc(l+2, 1);
  tmp_seq[1] = (uint8_t*)calloc(l+2, 1);
  tmp_seq_mem[0] = tmp_seq_mem[1] = l+2;
  size[0] = size_l; size[1] = size_r;

  tot_len = n_ref = 0;
  if(NULL != fp_fai) {
      int dummy_int[3];
      while(0 < fscanf(fp_fai, "%s\t%d\t%d\t%d\t%d", name, &l, &dummy_int[0], &dummy_int[1], &dummy_int[2])) {
          fprintf(stderr, "[dwgsim_core] %s length: %d\n", name, l);
          tot_len += l;
          ++n_ref;
      }
  }
  else {
      while ((l = seq_read_fasta(fp_fa, &seq, name, 0)) >= 0) {
          fprintf(stderr, "[dwgsim_core] %s length: %d\n", name, l);
          tot_len += l;
          ++n_ref;
      }
  }
  fprintf(stderr, "[dwgsim_core] %d sequences, total length: %llu\n", n_ref, (long long)tot_len);
  rewind(fp_fa);

  fprintf(stderr, "[dwgsim_core] Currently on: \n0");
  while ((l = seq_read_fasta(fp_fa, &seq, name, 0)) >= 0) {
      uint64_t n_pairs;
      n_ref--;
      if(0 == n_ref) {
          n_pairs = N - n_sim;
      }
      else {
          n_pairs = (uint64_t)((long double)l / tot_len * N + 0.5);
      }
      if (0 < size_r && l < dist + 3 * std_dev) {
          if(0 == prev_skip) fprintf(stderr, "\n");
          prev_skip = 1;
          fprintf(stderr, "[dwgsim_core] skip sequence '%s' as it is shorter than %d!\n", name, dist + 3 * std_dev);
          continue;
      }
      prev_skip = 0;

      // generate mutations and print them out
      maq_mut_diref(&seq, is_hap, rseq, rseq+1);
      maq_print_mutref(name, &seq, rseq, rseq+1, fpout0);

      for (ii = 0; ii != n_pairs; ++ii, ++ctr) { // the core loop
          if(0 == (ctr % 10000)) {
              fprintf(stderr, "\r[dwgsim_core] %llu",
                      (unsigned long long int)ctr);
          }
          double ran;
          int d, pos, s[2], strand[2];
          int n_sub[2], n_indel[2], n_err[2], ext_coor[2]={0,0}, j, k;
          int n_sub_first[2], n_indel_first[2], n_err_first[2]; // need this for SOLID data
          int c1, c2, c;
              
          s[0] = size[0]; s[1] = size[1];

          if(RAND_READ < drand48()) { 

              do { // avoid boundary failure
                  ran = ran_normal();
                  ran = ran * std_dev + dist;
                  d = (int)(ran + 0.5);
                  pos = (int)((l - d + 1) * drand48());
              } while (pos < 0 || pos >= seq.l || pos + d - 1 >= seq.l);

              if(2 == strandedness || (0 == strandedness && ILLUMINA == DATA_TYPE)) {
                  // opposite strand by default for Illumina
                  strand[0] = 0; strand[1] = 1; 
              }
              else if(1 == strandedness || (0 == strandedness && (SOLID == DATA_TYPE || IONTORRENT == DATA_TYPE))) {
                  // same strands by default for SOLiD
                  strand[0] = 0; strand[1] = 0; 
              }
              else {
                  // should not reach here
                  assert(1 == 0);
              }

              if (drand48() < 0.5) { // which strand ?
                  // Flip strands 
                  strand[0] = (1 + strand[0]) % 2;
                  strand[1] = (1 + strand[1]) % 2;
              }

              // generate the read sequences
              target = rseq[drand48()<0.5?0:1].s; // haplotype from which the reads are generated
              n_sub[0] = n_sub[1] = n_indel[0] = n_indel[1] = n_err[0] = n_err[1] = 0;
              n_sub_first[0] = n_sub_first[1] = n_indel_first[0] = n_indel_first[1] = n_err_first[0] = n_err_first[1] = 0;
              num_n[0]=num_n[1]=0;
                  
              // generate the reads in base space
              if(0 == strand[0]) {
                  if(0 < s[0]) { // F[FR]
                      __gen_read(0, pos, ++i); // + strand
                  }
                  if(0 < s[1]) {
                      if(1 == strand[1]) { // FR
                          __gen_read(1, pos + s[0] + s[1] + d - 1, --i); // - strand
                      }
                      else { // FF
                          __gen_read(1, pos + d + s[0], ++i); // + strand
                      }
                  }
              }
              else { // (1 == strand[0]) 
                  if(0 < s[0]) { // R[FR]
                      __gen_read(0, pos + s[0] + s[1] + d - 1, --i); // - strand
                  }
                  if(0 < s[1]) {
                      if(0 == strand[1]) { // RF
                          __gen_read(1, pos, ++i); // + strand
                      }
                      else {
                          __gen_read(1, pos + s[1] - 1, --i); // - strand
                      }
                  }
              }

              // Count # of Ns
              for (j = 0; j < 2; ++j) {
                  num_n[j]=0;
                  if(0 < s[j]) {
                      for (i = 0; i < s[j]; ++i) {
                          if(tmp_seq[j][i] == 4) num_n[j]++;
                      }
                  }
              }

              if (ext_coor[0] < 0 || ext_coor[1] < 0 || max_n < num_n[0] || max_n < num_n[1]) { // fail to generate the read(s)
                  --ii;
                  --ctr;
                  continue;
              }

              if(SOLID == DATA_TYPE) {
                  // Convert to color sequence, use the first base as the adaptor
                  for (j = 0; j < 2; ++j) {
                      if(0 < s[j]) {
                          c1 = 0; // adaptor 
                          for (i = 0; i < s[j]; ++i) {
                              c2 = tmp_seq[j][i]; // current base
                              c = __gf_add(c1, c2);
                              tmp_seq[j][i] = c;
                              c1 = c2; // save previous base
                          }
                      }
                  }
              }

              // generate sequencing errors
              if(IONTORRENT == DATA_TYPE) {
                  s[0] = generate_errors_flows(&tmp_seq[0], &tmp_seq_mem[0], s[0], strand[0], e[0]->start, &n_err[0]);
                  s[1] = generate_errors_flows(&tmp_seq[1], &tmp_seq_mem[1], s[1], strand[1], e[1]->start, &n_err[1]);
              }
              else { // Illumina/SOLiD
                  if(0 < s[0]) {
                      if(0 == strand[0]) { 
                          __gen_errors_mismatches(tmp_seq, 0, 0, ++i, s[0]); 
                      }
                      else { 
                          __gen_errors_mismatches(tmp_seq, 0, s[0]-1, --i, s[0]); 
                      }
                  }
                  if(0 < s[1]) {
                      if(0 == strand[1]) { 
                          __gen_errors_mismatches(tmp_seq, 1, 0, ++i, s[1]); 
                      }
                      else { 
                          __gen_errors_mismatches(tmp_seq, 1, s[1]-1, --i, s[1]); 
                      }
                  }
              }

              // print
              for (j = 0; j < 2; ++j) {
                  if(s[j] <= 0) {
                      continue;
                  }
                  if(IONTORRENT == DATA_TYPE && qstr_l < s[j]) {
                      qstr_l = s[j];
                      qstr = realloc(qstr, (1+qstr_l) * sizeof(char));
                  }
                  for (i = 0; i < s[j]; ++i) {
                      qstr[i] = (int)(-10.0 * log(e[j]->start + e[j]->by*i) / log(10.0) + 0.499) + 33;
                  }
                  qstr[i] = 0;
                  // BWA
                  FILE *fpo = (0 == j) ? fpout2 : fpout3;
                  if(ILLUMINA == DATA_TYPE || IONTORRENT == DATA_TYPE) {
                      fprintf(fpo, "@%s_%u_%u_%1u_%1u_%1u_%1u_%d:%d:%d_%d:%d:%d_%llx/%d\n", 
                              name, ext_coor[0]+1, ext_coor[1]+1, strand[0], strand[1], 0, 0,
                              n_err[0] - n_err_first[0], 
                              n_sub[0] - n_sub_first[0], 
                              n_indel[0] - n_indel_first[0], 
                              n_err[1] - n_err_first[1], 
                              n_sub[1] - n_sub_first[1], 
                              n_indel[1] - n_indel_first[1],
                              (long long)ii,
                              j+1);
                      for (i = 0; i < s[j]; ++i)
                        fputc("ACGTN"[(int)tmp_seq[j][i]], fpo);
                      fprintf(fpo, "\n+\n%s\n", qstr);
                  }
                  else {
                      // Note: BWA ignores the adapter and the first color, so this is a misrepresentation 
                      // in samtools.  We must first skip the first color.  Basically, a 50 color read is a 
                      // 49 color read for BWA.
                      //
                      // Note: BWA outputs F3 to read1, annotated as read "2", and outputs R3 to read2,
                      // annotated as read "1".
                      fprintf(fpo, "@%s_%u_%u_%1u_%1u_%1u_%1u_%d:%d:%d_%d:%d:%d_%llx/%d\n", 
                              name, ext_coor[0]+1, ext_coor[1]+1, strand[0], strand[1], 0, 0,
                              n_err[0], n_sub[0], n_indel[0], n_err[1], n_sub[1], n_indel[1],
                              (long long)ii, 2 - j);
                      //fputc('A', fpo);
                      for (i = 1; i < s[j]; ++i)
                        fputc("ACGTN"[(int)tmp_seq[j][i]], fpo);
                      fprintf(fpo, "\n+\n");
                      for (i = 1; i < s[j]; ++i) 
                        fputc(qstr[i], fpo);
                      fprintf(fpo, "\n");
                  }

                  // BFAST output
                  fprintf(fpout1, "@%s_%u_%u_%1u_%1u_%1u_%1u_%d:%d:%d_%d:%d:%d_%llx\n", 
                          name, ext_coor[0]+1, ext_coor[1]+1, strand[0], strand[1], 0, 0,
                          n_err[0], n_sub[0], n_indel[0], n_err[1], n_sub[1], n_indel[1],
                          (long long)ii);
                  if(ILLUMINA == DATA_TYPE || IONTORRENT == DATA_TYPE) {
                      for (i = 0; i < s[j]; ++i)
                        fputc("ACGTN"[(int)tmp_seq[j][i]], fpout1);
                      fprintf(fpout1, "\n+\n%s\n", qstr);
                  }
                  else {
                      fputc('A', fpout1);
                      for (i = 0; i < s[j]; ++i)
                        fputc("01234"[(int)tmp_seq[j][i]], fpout1);
                      fprintf(fpout1, "\n+\n");
                      for (i = 0; i < s[j]; ++i) 
                        fputc(qstr[i], fpout1);
                      fprintf(fpout1, "\n");
                  }
              }
              n_sim++;
          }
          else { // random DNA read
              for(j=0;j<2;j++) {
                  if(s[j] <= 0) {
                     continue;
                  } 
                  if(IONTORRENT == DATA_TYPE && qstr_l < s[j]) {
                      qstr_l = s[j];
                      qstr = realloc(qstr, (1+qstr_l) * sizeof(char));
                  }
                  // get random sequence
                  for(i=0;i<s[j];i++) {
                      tmp_seq[j][i] = (int)(drand48() * 4.0) & 3;
                      qstr[i] = (int)(-10.0 * log(e[j]->start + e[j]->by*i) / log(10.0) + 0.499) + 33;
                  }
                  qstr[i] = 0;
                  if(SOLID == DATA_TYPE) { // convert to color space
                      if(0 < s[j]) {
                          c1 = 0; // adaptor 
                          for (i = 0; i < s[j]; ++i) {
                              c2 = tmp_seq[j][i]; // current base
                              c = __gf_add(c1, c2);
                              tmp_seq[j][i] = c;
                              c1 = c2; // save previous base
                          }
                      }
                  }
                  // BWA
                  FILE *fpo = (0 == j) ? fpout2 : fpout3;
                  if(ILLUMINA == DATA_TYPE || IONTORRENT == DATA_TYPE) {
                      fprintf(fpo, "@%s_%u_%u_%1u_%1u_%1u_%1u_%d:%d:%d_%d:%d:%d_%llx/%d\n", 
                              "rand", 0, 0, 0, 0, 1, 1,
                              0, 0, 0, 0, 0, 0,
                              (long long)ii,
                              j+1);
                      for (i = 0; i < s[j]; ++i)
                        fputc("ACGTN"[(int)tmp_seq[j][i]], fpo);
                      fprintf(fpo, "\n+\n%s\n", qstr);
                  }
                  else {
                      // Note: BWA ignores the adapter and the first color, so this is a misrepresentation 
                      // in samtools.  We must first skip the first color.  Basically, a 50 color read is a 
                      // 49 color read for BWA.
                      //
                      // Note: BWA outputs F3 to read1, annotated as read "2", and outputs R3 to read2,
                      // annotated as read "1".
                      fprintf(fpo, "@%s_%u_%u_%1u_%1u_%1u_%1u_%d:%d:%d_%d:%d:%d_%llx/%d\n", 
                              "rand", 0, 0, 0, 0, 1, 1,
                              0, 0, 0, 0, 0, 0,
                              (long long)ii, 2 - j);
                      //fputc('A', fpo);
                      for (i = 1; i < s[j]; ++i)
                        fputc("ACGTN"[(int)tmp_seq[j][i]], fpo);
                      fprintf(fpo, "\n+\n");
                      for (i = 1; i < s[j]; ++i) 
                        fputc(qstr[i], fpo);
                      fprintf(fpo, "\n");
                  }

                  // BFAST output
                  fprintf(fpout1, "@%s_%u_%u_%1u_%1u_%1u_%1u_%d:%d:%d_%d:%d:%d_%llx\n", 
                          "rand", 0, 0, 0, 0, 1, 1,
                          0, 0, 0, 0, 0, 0,
                          (long long)ii);
                  if(ILLUMINA == DATA_TYPE || IONTORRENT == DATA_TYPE) {
                      for (i = 0; i < s[j]; ++i)
                        fputc("ACGTN"[(int)tmp_seq[j][i]], fpout1);
                      fprintf(fpout1, "\n+\n%s\n", qstr);
                  }
                  else {
                      fputc('A', fpout1);
                      for (i = 0; i < s[j]; ++i)
                        fputc("01234"[(int)tmp_seq[j][i]], fpout1);
                      fprintf(fpout1, "\n+\n");
                      for (i = 0; i < s[j]; ++i) 
                        fputc(qstr[i], fpout1);
                      fprintf(fpout1, "\n");
                  }
              }
              n_sim++;
          }
      }
      free(rseq[0].s); free(rseq[1].s);
      fprintf(stderr, "\r[dwgsim_core] %llu",
              (unsigned long long int)ctr);
  }
  fprintf(stderr, "\r[dwgsim_core] Complete!\n");
  free(seq.s); free(qstr);
  free(tmp_seq[0]); free(tmp_seq[1]);
}

void get_error_rate(const char *str, error_t *e)
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

static void check_option_int(int32_t val, int32_t min, int32_t max, char *opt)
{
  if(val < min || max < val) {
      fprintf(stderr, "Error: command line option %s was out of range\n", opt);
      exit(1);
  }
}

static int simu_usage()
{
  fprintf(stderr, "\n");
  fprintf(stderr, "Program: dwgsim (short read simulator)\n");
  fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
  fprintf(stderr, "Contact: Nils Homer <dnaa-help@lists.sourceforge.net>\n\n");
  fprintf(stderr, "Usage:   dwgsim [options] <in.ref.fa> <out.prefix>\n\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "         -e FLOAT      base/color error rate of the first read [%.3f]\n", ERR_RATE);
  fprintf(stderr, "         -E FLOAT      base/color error rate of the second read [%.3f]\n", ERR_RATE);
  fprintf(stderr, "         -d INT        inner distance between the two ends [500]\n");
  fprintf(stderr, "         -s INT        standard deviation [50]\n");
  fprintf(stderr, "         -N INT        number of read pairs [1000000]\n");
  fprintf(stderr, "         -1 INT        length of the first read [70]\n");
  fprintf(stderr, "         -2 INT        length of the second read [70]\n");
  fprintf(stderr, "         -r FLOAT      rate of mutations [%.4f]\n", MUT_RATE);
  fprintf(stderr, "         -R FLOAT      fraction of indels [%.2f]\n", INDEL_FRAC);
  fprintf(stderr, "         -X FLOAT      probability an indel is extended [%.2f]\n", INDEL_EXTEND);
  fprintf(stderr, "         -y FLOAT      probability of a random DNA read [%.2f]\n", RAND_READ);
  fprintf(stderr, "         -n INT        maximum number of Ns allowed in a given read [%d]\n", MAX_N);
  fprintf(stderr, "         -c INT        generate reads for [%d]:\n", DATA_TYPE);
  fprintf(stderr, "                           0: Illumina\n");
  fprintf(stderr, "                           1: SOLiD\n");
  fprintf(stderr, "                           2: Ion Torrent\n");
  fprintf(stderr, "         -S INT        generate reads [%d]:\n", DATA_TYPE);
  fprintf(stderr, "                           0: default (opposite strand for Illumina, same strand for SOLiD/Ion Torrent)\n");
  fprintf(stderr, "                           1: same strand (mate pair)\n");
  fprintf(stderr, "                           2: opposite strand (paired end)\n");
  fprintf(stderr, "         -f            the flow order for Ion Torrent data [%s]\n", FLOW_ORDER);
  fprintf(stderr, "         -H            haploid mode\n");
  fprintf(stderr, "         -h            print this message\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Note: For SOLiD mate pair reads and BFAST, the first read is F3 and the second is R3. For SOLiD mate pair reads\n");
  fprintf(stderr, "and BWA, the reads in the first file are R3 the reads annotated as the first read etc.\n");
  return 1;
}

int main(int argc, char *argv[])
{
  int64_t N;
  int i, dist, std_dev, c, size_l, size_r, is_hap = 0, strandedness = 0;
  FILE *fpout0, *fpout1, *fpout2, *fpout3, *fp_fa, *fp_fai;
  error_t e1, e2;
  char fn_fai[1024]="\0";
  char fn_tmp[1024]="\0";
  int max_n;

  e1.start = e1.end = e2.start = e2.end = ERR_RATE;

  N = 1000000; dist = 500; std_dev = 50;
  size_l = size_r = 70;
  max_n = MAX_N;
  strandedness = 0;
  while ((c = getopt(argc, argv, "d:s:N:1:2:e:E:r:R:X:c:S:n:y:Hf:h")) >= 0) {
      switch (c) {
        case 'd': dist = atoi(optarg); break;
        case 's': std_dev = atoi(optarg); break;
        case 'N': N = atoi(optarg); break;
        case '1': size_l = atoi(optarg); break;
        case '2': size_r = atoi(optarg); break;
        case 'e': get_error_rate(optarg, &e1); break;
        case 'E': get_error_rate(optarg, &e2); break;
        case 'r': MUT_RATE = atof(optarg); break;
        case 'R': INDEL_FRAC = atof(optarg); break;
        case 'X': INDEL_EXTEND = atof(optarg); break;
        case 'c': DATA_TYPE = atoi(optarg); break;
        case 'S': strandedness = atoi(optarg); break;
        case 'n': max_n = atoi(optarg); break;
        case 'y': RAND_READ = atof(optarg); break;
        case 'f': strcpy((char*)FLOW_ORDER, optarg); break;
        case 'H': is_hap = 1; break;
        case 'h': return simu_usage();
        default: fprintf(stderr, "Unrecognized option: -%c\n", c); return 1;
      }
  }
  if (argc - optind < 1) return simu_usage();

  check_option_int(DATA_TYPE, 0, 2, "-c");
  check_option_int(strandedness, 0, 2, "-s");

  FLOW_ORDER_LEN = strlen((char*)FLOW_ORDER);
  for(i=0;i<FLOW_ORDER_LEN;i++) {
      FLOW_ORDER[i] = nst_nt4_table[FLOW_ORDER[i]];
  }

  e1.by = (e1.end - e1.start) / size_l;
  e2.by = (e2.end - e2.start) / size_r;
  if(IONTORRENT == DATA_TYPE) {
      if(e1.end != e1.start) {
          fprintf(stderr, "End one: a uniform error rate must be given for Ion Torrent data");
          return 1;
      }
      if(e2.end != e2.start) {
          fprintf(stderr, "End two: a uniform error rate must be given for Ion Torrent data");
          return 1;
      }
  }

  // Open files
  fp_fa =	xopen(argv[optind+0], "r");
  strcpy(fn_fai, argv[optind+0]); strcat(fn_fai, ".fai");
  fp_fai = fopen(fn_fai, "r"); // depends on returning NULL;
  strcpy(fn_tmp, argv[optind+1]); strcat(fn_tmp, ".mutations.txt");
  fpout0 = xopen(fn_tmp, "w");
  strcpy(fn_tmp, argv[optind+1]); strcat(fn_tmp, ".bfast.fastq");
  fpout1 = xopen(fn_tmp, "w");
  strcpy(fn_tmp, argv[optind+1]); strcat(fn_tmp, ".bwa.read1.fastq");
  fpout2 = xopen(fn_tmp, "w");
  strcpy(fn_tmp, argv[optind+1]); strcat(fn_tmp, ".bwa.read2.fastq");
  fpout3 = xopen(fn_tmp, "w");

  // Run simulation
  dwgsim_core(fpout0, fpout1, fpout2, fpout3, fp_fa, fp_fai, &e1, &e2, is_hap, N, dist, std_dev, size_l, size_r, max_n, strandedness);

  // Close files
  fclose(fp_fa); fclose(fpout1); fclose(fpout2); fclose(fpout3); 
  if(NULL != fp_fai) fclose(fp_fai);

  return 0;
}
