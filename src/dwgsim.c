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

// TODO: refactor a lot of this code into separate source files

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
#define __IS_TRUE(_val) ((_val == 1) ? "True" : "False")

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

enum data_type_t {
    ILLUMINA=0,
    SOLID=1,
    IONTORRENT=2
};

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

typedef struct {
    double start, by, end;
} error_t;

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

typedef struct {
    error_t e1;
    error_t e2;
    int32_t dist;
    double std_dev;
    int64_t N;
    int32_t length1;
    int32_t length2;
    double mut_rate;
    double indel_frac;
    double indel_extend;
    double rand_read;
    int32_t max_n;
    int32_t data_type;
    int32_t strandedness;
    int8_t *flow_order;
    int32_t flow_order_len;
    int32_t is_hap;
    int32_t seed;
    char *fn_muts_bed;
    FILE *fp_mut;
    FILE *fp_bfast;
    FILE *fp_bwa1;
    FILE *fp_bwa2;
    FILE *fp_fa;
    FILE *fp_fai;
} dwgsim_opt_t;

dwgsim_opt_t* dwgsim_opt_init()
{
  dwgsim_opt_t *opt;
  opt = calloc(1, sizeof(dwgsim_opt_t));
  opt->e1.start = opt->e1.end = opt->e2.start = opt->e2.end = 0.02;
  opt->e1.by = opt->e2.by = 0;
  opt->dist = 500;
  opt->N = 1000000;
  opt->length1 = opt->length2 = 70;
  opt->mut_rate = 0.001;
  opt->indel_frac = 0.1;
  opt->indel_extend = 0.3;
  opt->rand_read = 0.05;
  opt->data_type = ILLUMINA;
  opt->max_n = 0;
  opt->flow_order = NULL;
  opt->flow_order_len = 0;
  opt->seed = -1;
  opt->fn_muts_bed = NULL;
  opt->fp_mut = opt->fp_bfast = opt->fp_bwa1 = opt->fp_bwa2 = NULL;
  opt->fp_fa = opt->fp_fai = NULL;

  return opt;
}

void dwgsim_opt_destroy(dwgsim_opt_t *opt)
{
  free(opt->fn_muts_bed);
  free(opt->flow_order);
  free(opt);
}

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

int32_t ran_num(double prob, int32_t n)
{
  int32_t i, r;
  for(i = r = 0; i < n; i++) {
      if(drand48() < prob) r++;
  }
  return r;
}

typedef struct {
    char *name;
    int32_t len;
} contig_t;

typedef struct {
    contig_t *contigs;
    int32_t n;
} contigs_t;

contigs_t* contigs_init()
{
  return calloc(1, sizeof(contigs_t));
}

void contigs_add(contigs_t *c, char *name, uint32_t len)
{
  c->n++;
  c->contigs = realloc(c->contigs, c->n * sizeof(contig_t));
  c->contigs[c->n-1].name = strdup(name);
  c->contigs[c->n-1].len = len;
}

void contigs_destroy(contigs_t *c)
{
  int32_t i;
  for(i=0;i<c->n;i++) {
      free(c->contigs[i].name);
  }
  free(c->contigs);
  free(c);
} 

int32_t get_muttype(char *str)
{
  int32_t i;
  for(i=0;i<strlen(str);i++) {
      str[i] = tolower(str[i]);
  }
  if(0 == strcmp("snp", str) || 0 == strcmp("substitute", str) || 0 == strcmp("sub", str) || 0 == strcmp("s", str)) {
      return SUBSTITUTE;
  }
  else if(0 == strcmp("insertion", str) || 0 == strcmp("insert", str) || 0 == strcmp("ins", str) || 0 == strcmp("i", str)) {
      return INSERT;
  }
  else if(0 == strcmp("deletion", str) || 0 == strcmp("delet", str) || 0 == strcmp("del", str) || 0 == strcmp("d", str)) {
      return DELETE;
  }
  return -1;
}

typedef struct {
    uint32_t contig; // zero-based
    uint32_t start; // zero-based
    uint32_t end; // zero-based end position
    uint8_t type; // mutation type
    char *bases; // mutation bases
} mut_bed_t;

typedef struct {
    mut_bed_t *muts;
    int32_t n;
    int32_t mem;
} muts_bed_t;

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
      else if(prev_contig == i && start <= max_end) {
          fprintf(stderr, "Error: overlapping entries\n");
          exit(1);
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
          if(((ins_length_shift - muttype_shift) >> 1) < end - start) {
              fprintf(stderr, "Error: insertion of length %d exceeded the maximum supported length of %d\n",
                      end - start,
                      (int32_t)((ins_length_shift - muttype_shift) >> 1));
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
generate_errors_flows(dwgsim_opt_t *opt, uint8_t **seq, int32_t *mem, int32_t len, uint8_t strand, double e, int32_t *_n_err)
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
  for(i=0;i<opt->flow_order_len;i++) {
      c = (4 <= (*seq)[0]) ? 0 : (*seq)[0];
      if(c == opt->flow_order[i]) {
          break;
      }
  }
  assert(opt->flow_order_len != i);
  flow_i = i;
  // first pass: add errors in the current sequence (non-empty flows)
  prev_c = 4;
  for(i=0;i<len;i++) {
      c = (4 <= (*seq)[i]) ? 0 : (*seq)[i];
      while(c != opt->flow_order[flow_i]) { // skip paste empty flows
          flow_i = (flow_i + 1) % opt->flow_order_len;
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
                      while(next_c != opt->flow_order[(flow_i + j) % opt->flow_order_len]) { // skip paste empty flows
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
                      (*seq)[i] = opt->flow_order[(flow_i + k) % opt->flow_order_len];
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
      while(c != opt->flow_order[flow_i]) {
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
                  (*seq)[j] = opt->flow_order[flow_i];
              }
              len += n_err;
              (*_n_err) += n_err;
          }
          flow_i = (flow_i + 1) % opt->flow_order_len;
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

void maq_mut_diref(dwgsim_opt_t *opt, const seq_t *seq, mutseq_t *hap1, mutseq_t *hap2, int32_t contig_i, muts_bed_t *muts_bed)
{
  int i, deleting = 0;
  mutseq_t *ret[2];

  ret[0] = hap1; ret[1] = hap2;
  ret[0]->l = seq->l; ret[1]->l = seq->l;
  ret[0]->m = seq->m; ret[1]->m = seq->m;
  ret[0]->s = (mut_t *)calloc(seq->m, sizeof(mut_t));
  ret[1]->s = (mut_t *)calloc(seq->m, sizeof(mut_t));

  if(NULL == muts_bed) {
      for (i = 0; i != seq->l; ++i) {
          mut_t c;
          c = ret[0]->s[i] = ret[1]->s[i] = (mut_t)nst_nt4_table[(int)seq->s[i]];
          if (deleting) {
              if (drand48() < opt->indel_extend) {
                  if (deleting & 1) ret[0]->s[i] |= DELETE|c;
                  if (deleting & 2) ret[1]->s[i] |= DELETE|c;
                  continue;
              } else deleting = 0;
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
                  } else { // insertion
                      mut_t num_ins = 0, ins = 0;
                      do {
                          num_ins++;
                          ins = (ins << 2) | (mut_t)(drand48() * 4.0);
                      } while (num_ins < ((ins_length_shift - muttype_shift) >> 1) && drand48() < opt->indel_extend);
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
  else {
      // seed
      for (i = 0; i != seq->l; ++i) {
          ret[0]->s[i] = ret[1]->s[i] = (mut_t)nst_nt4_table[(int)seq->s[i]];
      }
      // mutates exactly based on a BED file
      for (i = 0; i < muts_bed->n; ++i) {
          if (muts_bed->muts[i].contig == contig_i) {
              int32_t j;
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
                       muts_bed->muts[i].start <= j && num_ins < ((ins_length_shift - muttype_shift) >> 1); 
                       --j) { // for each base
                      num_ins++;
                      if(0 == has_bases) {
                          ins = (ins << 2) | (mut_t)(drand48() * 4.0);
                      }
                      else {
                          ins = (ins << 2) | (mut_t)(nst_nt4_table[(int)muts_bed->muts[i].bases[j - muts_bed->muts[i].start]]);
                      }
                  } while (num_ins < ((ins_length_shift - muttype_shift) >> 1) && drand48() < opt->indel_extend);
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
  int j, del_length;
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
                  mut_t n = (c[1] >> ins_length_shift) & ins_length_mask, ins = c[1] >> muttype_shift;
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
                      ins |= (hap1->s[j-1]&3); // insert the first base
                      hap1->s[j] = (hap1->s[j]&3); // make it NOCHANGE
                      j--;
                  }
                  hap1->s[j] = (n << ins_length_shift) | (ins << muttype_shift) | INSERT | (hap1->s[j]&3); // re-insert
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

// Columns:
// 1 - chromosome name
// 2 - position (one-based)
// 3 - reference base (dash if there was an insertion)
// 4 - variant allele (IUPAC code or insertion base(s))
// 5 - '-' for homozygous, '+' for heterozygous
void maq_print_mutref(const char *name, const seq_t *seq, mutseq_t *hap1, mutseq_t *hap2, FILE *fpout)
{
  int32_t i;
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

void dwgsim_core(dwgsim_opt_t * opt)
{
  seq_t seq;
  mutseq_t rseq[2];
  uint64_t tot_len, ii=0, ctr=0;
  int i, l, n_ref, contig_i;
  char name[256], *qstr;
  int size[2], prev_skip=0, qstr_l=0;
  int num_n[2];
  uint8_t *tmp_seq[2]={NULL,NULL};
  int32_t tmp_seq_mem[2]={0,0};
  uint64_t n_sim = 0;
  mut_t *target;
  error_t *e[2]={NULL,NULL};
  FILE *fp_muts_bed = NULL;
  muts_bed_t *muts_bed = NULL;
  contigs_t *contigs = NULL;

  e[0] = &opt->e1; e[1] = &opt->e2;

  INIT_SEQ(seq);
  seq_set_block_size(0x1000000);
  l = opt->length1 > opt->length2? opt->length1 : opt->length2;
  qstr_l = l;
  qstr = (char*)calloc(qstr_l+1, 1);
  tmp_seq[0] = (uint8_t*)calloc(l+2, 1);
  tmp_seq[1] = (uint8_t*)calloc(l+2, 1);
  tmp_seq_mem[0] = tmp_seq_mem[1] = l+2;
  size[0] = opt->length1; size[1] = opt->length2;
  
  if(NULL != opt->fn_muts_bed) {
      fp_muts_bed = xopen(opt->fn_muts_bed, "r");
      contigs = contigs_init();
  }
  
  tot_len = n_ref = 0;
  if(NULL != opt->fp_fai) {
      int dummy_int[3];
      while(0 < fscanf(opt->fp_fai, "%s\t%d\t%d\t%d\t%d", name, &l, &dummy_int[0], &dummy_int[1], &dummy_int[2])) {
          fprintf(stderr, "[dwgsim_core] %s length: %d\n", name, l);
          tot_len += l;
          ++n_ref;
          if(NULL != contigs) {
              contigs_add(contigs, name, l);
          }
      }
  }
  else {
      while ((l = seq_read_fasta(opt->fp_fa, &seq, name, 0)) >= 0) {
          fprintf(stderr, "[dwgsim_core] %s length: %d\n", name, l);
          tot_len += l;
          ++n_ref;
          if(NULL != contigs) {
              contigs_add(contigs, name, l);
          }
      }
  }
  fprintf(stderr, "[dwgsim_core] %d sequences, total length: %llu\n", n_ref, (long long)tot_len);
  rewind(opt->fp_fa);

  if(NULL != opt->fn_muts_bed) {
      muts_bed = muts_bed_init(fp_muts_bed, contigs); // read in the BED
      contigs_destroy(contigs);
      contigs = NULL;
  }
  
  fprintf(stderr, "[dwgsim_core] Currently on: \n0");
  contig_i = 0;
  while ((l = seq_read_fasta(opt->fp_fa, &seq, name, 0)) >= 0) {
      uint64_t n_pairs;
      n_ref--;
      if(0 == n_ref) {
          n_pairs = opt->N - n_sim;
      }
      else {
          n_pairs = (uint64_t)((long double)l / tot_len * opt->N + 0.5);
      }
      if (0 < opt->length2 && l < opt->dist + 3 * opt->std_dev) {
          if(0 == prev_skip) fprintf(stderr, "\n");
          prev_skip = 1;
          fprintf(stderr, "[dwgsim_core] skip sequence '%s' as it is shorter than %f!\n", name, opt->dist + 3 * opt->std_dev);
          continue;
      }
      prev_skip = 0;

      // generate mutations and print them out
      maq_mut_diref(opt, &seq, rseq, rseq+1, contig_i, muts_bed);
      maq_print_mutref(name, &seq, rseq, rseq+1, opt->fp_mut);

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

          if(opt->rand_read < drand48()) { 

              do { // avoid boundary failure
                  ran = ran_normal();
                  ran = ran * opt->std_dev + opt->dist;
                  d = (int)(ran + 0.5);
                  pos = (int)((l - d + 1) * drand48());
              } while (pos < 0 || pos >= seq.l || pos + d - 1 >= seq.l);

              if(2 == opt->strandedness || (0 == opt->strandedness && ILLUMINA == opt->data_type)) {
                  // opposite strand by default for Illumina
                  strand[0] = 0; strand[1] = 1; 
              }
              else if(1 == opt->strandedness || (0 == opt->strandedness && (SOLID == opt->data_type || IONTORRENT == opt->data_type))) {
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

              if (ext_coor[0] < 0 || ext_coor[1] < 0 || opt->max_n < num_n[0] || opt->max_n < num_n[1]) { // fail to generate the read(s)
                  --ii;
                  --ctr;
                  continue;
              }

              if(SOLID == opt->data_type) {
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
              if(IONTORRENT == opt->data_type) {
                  s[0] = generate_errors_flows(opt, &tmp_seq[0], &tmp_seq_mem[0], s[0], strand[0], e[0]->start, &n_err[0]);
                  s[1] = generate_errors_flows(opt, &tmp_seq[1], &tmp_seq_mem[1], s[1], strand[1], e[1]->start, &n_err[1]);
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
                  if(IONTORRENT == opt->data_type && qstr_l < s[j]) {
                      qstr_l = s[j];
                      qstr = realloc(qstr, (1+qstr_l) * sizeof(char));
                  }
                  for (i = 0; i < s[j]; ++i) {
                      qstr[i] = (int)(-10.0 * log(e[j]->start + e[j]->by*i) / log(10.0) + 0.499) + 33;
                  }
                  qstr[i] = 0;
                  // BWA
                  FILE *fpo = (0 == j) ? opt->fp_bwa1: opt->fp_bwa2;
                  if(ILLUMINA == opt->data_type || IONTORRENT == opt->data_type) {
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
                  fprintf(opt->fp_bfast, "@%s_%u_%u_%1u_%1u_%1u_%1u_%d:%d:%d_%d:%d:%d_%llx\n", 
                          name, ext_coor[0]+1, ext_coor[1]+1, strand[0], strand[1], 0, 0,
                          n_err[0], n_sub[0], n_indel[0], n_err[1], n_sub[1], n_indel[1],
                          (long long)ii);
                  if(ILLUMINA == opt->data_type || IONTORRENT == opt->data_type) {
                      for (i = 0; i < s[j]; ++i)
                        fputc("ACGTN"[(int)tmp_seq[j][i]], opt->fp_bfast);
                      fprintf(opt->fp_bfast, "\n+\n%s\n", qstr);
                  }
                  else {
                      fputc('A', opt->fp_bfast);
                      for (i = 0; i < s[j]; ++i)
                        fputc("01234"[(int)tmp_seq[j][i]], opt->fp_bfast);
                      fprintf(opt->fp_bfast, "\n+\n");
                      for (i = 0; i < s[j]; ++i) 
                        fputc(qstr[i], opt->fp_bfast);
                      fprintf(opt->fp_bfast, "\n");
                  }
              }
              n_sim++;
          }
          else { // random DNA read
              for(j=0;j<2;j++) {
                  if(s[j] <= 0) {
                     continue;
                  } 
                  if(IONTORRENT == opt->data_type && qstr_l < s[j]) {
                      qstr_l = s[j];
                      qstr = realloc(qstr, (1+qstr_l) * sizeof(char));
                  }
                  // get random sequence
                  for(i=0;i<s[j];i++) {
                      tmp_seq[j][i] = (int)(drand48() * 4.0) & 3;
                      qstr[i] = (int)(-10.0 * log(e[j]->start + e[j]->by*i) / log(10.0) + 0.499) + 33;
                  }
                  qstr[i] = 0;
                  if(SOLID == opt->data_type) { // convert to color space
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
                  FILE *fpo = (0 == j) ? opt->fp_bwa1: opt->fp_bwa2;
                  if(ILLUMINA == opt->data_type || IONTORRENT == opt->data_type) {
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
                  fprintf(opt->fp_bfast, "@%s_%u_%u_%1u_%1u_%1u_%1u_%d:%d:%d_%d:%d:%d_%llx\n", 
                          "rand", 0, 0, 0, 0, 1, 1,
                          0, 0, 0, 0, 0, 0,
                          (long long)ii);
                  if(ILLUMINA == opt->data_type || IONTORRENT == opt->data_type) {
                      for (i = 0; i < s[j]; ++i)
                        fputc("ACGTN"[(int)tmp_seq[j][i]], opt->fp_bfast);
                      fprintf(opt->fp_bfast, "\n+\n%s\n", qstr);
                  }
                  else {
                      fputc('A', opt->fp_bfast);
                      for (i = 0; i < s[j]; ++i)
                        fputc("01234"[(int)tmp_seq[j][i]], opt->fp_bfast);
                      fprintf(opt->fp_bfast, "\n+\n");
                      for (i = 0; i < s[j]; ++i) 
                        fputc(qstr[i], opt->fp_bfast);
                      fprintf(opt->fp_bfast, "\n");
                  }
              }
              n_sim++;
          }
      }
      free(rseq[0].s); free(rseq[1].s);
      fprintf(stderr, "\r[dwgsim_core] %llu",
              (unsigned long long int)ctr);
      contig_i++;
  }
  fprintf(stderr, "\r[dwgsim_core] Complete!\n");
  free(seq.s); free(qstr);
  free(tmp_seq[0]); free(tmp_seq[1]);
  if(NULL != opt->fn_muts_bed) {
      muts_bed_destroy(muts_bed);
      fclose(fp_muts_bed);
  }
}

#define __check_option(_val, _min, _max, _opt) \
  if(_val < _min || _max < _val) { \
      fprintf(stderr, "Error: command line option %s was out of range\n", _opt); \
      exit(1); \
  } 

static int simu_usage(dwgsim_opt_t *opt)
{
  fprintf(stderr, "\n");
  fprintf(stderr, "Program: dwgsim (short read simulator)\n");
  fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
  fprintf(stderr, "Contact: Nils Homer <dnaa-help@lists.sourceforge.net>\n\n");
  fprintf(stderr, "Usage:   dwgsim [options] <in.ref.fa> <out.prefix>\n\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "         -e FLOAT      base/color error rate of the first read [from %.3f to %.3f by %.3f]\n", opt->e1.start, opt->e1.end, opt->e1.by);
  fprintf(stderr, "         -E FLOAT      base/color error rate of the second read [from %.3f to %.3f by %.3f]\n", opt->e2.start, opt->e2.end, opt->e2.by);
  fprintf(stderr, "         -d INT        inner distance between the two ends [%d]\n", opt->dist);
  fprintf(stderr, "         -s INT        standard deviation [%.3f]\n", opt->std_dev);
  fprintf(stderr, "         -N INT        number of read pairs [%lld]\n", opt->N);
  fprintf(stderr, "         -1 INT        length of the first read [%d]\n", opt->length1);
  fprintf(stderr, "         -2 INT        length of the second read [%d]\n", opt->length2);
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
  fprintf(stderr, "         -H            haploid mode [%s]\n", __IS_TRUE(opt->is_hap));
  fprintf(stderr, "         -z INT        random seed (-1 uses the current time) [%d]\n", opt->seed);
  fprintf(stderr, "         -b FILE       the bed-like set of candidate mutations [%s]\n", (NULL == opt->fn_muts_bed) ? "not using" : opt->fn_muts_bed);
  fprintf(stderr, "         -h            print this message\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Note: For SOLiD mate pair reads and BFAST, the first read is F3 and the second is R3. For SOLiD mate pair reads\n");
  fprintf(stderr, "and BWA, the reads in the first file are R3 the reads annotated as the first read etc.\n");
  return 1;
}

int main(int argc, char *argv[])
{
  dwgsim_opt_t *opt = NULL;

  opt = dwgsim_opt_init();

  int i, c;
  char fn_fai[1024]="\0";
  char fn_tmp[1024]="\0";

  while ((c = getopt(argc, argv, "d:s:N:1:2:e:E:r:R:X:c:S:n:y:Hf:z:b:h")) >= 0) {
      switch (c) {
        case 'd': opt->dist = atoi(optarg); break;
        case 's': opt->std_dev = atof(optarg); break;
        case 'N': opt->N = atoi(optarg); break;
        case '1': opt->length1 = atoi(optarg); break;
        case '2': opt->length2 = atoi(optarg); break;
        case 'e': get_error_rate(optarg, &opt->e1); break;
        case 'E': get_error_rate(optarg, &opt->e2); break;
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
        case 'H': opt->is_hap = 1; break;
        case 'h': return simu_usage(opt);
        case 'z': opt->seed = atoi(optarg); break;
        case 'b': free(opt->fn_muts_bed); opt->fn_muts_bed = strdup(optarg); break;
        default: fprintf(stderr, "Unrecognized option: -%c\n", c); return 1;
      }
  }
  if (argc - optind < 1) return simu_usage(opt);

  __check_option(opt->dist, 0, INT32_MAX, "-d");
  __check_option(opt->std_dev, 0, INT32_MAX, "-s");
  __check_option(opt->N, 1, INT32_MAX, "-N");
  __check_option(opt->length1, 1, INT32_MAX, "-1");
  __check_option(opt->length2, 0, INT32_MAX, "-2");
  // error rate
  opt->e1.by = (opt->e1.end - opt->e1.start) / opt->length1;
  opt->e2.by = (opt->e2.end - opt->e2.start) / opt->length2;
  if(IONTORRENT == opt->data_type) {
      if(opt->e1.end != opt->e1.start) {
          fprintf(stderr, "End one: a uniform error rate must be given for Ion Torrent data");
          return 1;
      }
      if(opt->e2.end != opt->e2.start) {
          fprintf(stderr, "End two: a uniform error rate must be given for Ion Torrent data");
          return 1;
      }
  }
  __check_option(opt->mut_rate, 0, 1.0, "-r");
  __check_option(opt->indel_frac, 0, 1.0, "-R");
  __check_option(opt->indel_extend, 0, 1.0, "-X");
  __check_option(opt->data_type, 0, 2, "-c");
  __check_option(opt->strandedness, 0, 2, "-S");
  __check_option(opt->max_n, 0, INT32_MAX, "-n");
  __check_option(opt->rand_read, 0, 1.0, "-y");
  if(IONTORRENT == opt->data_type && NULL == opt->flow_order) {
      fprintf(stderr, "Error: command line option -f is required\n");
      return 1;
  }
  __check_option(opt->is_hap, 0, 1, "-H");
  
  // random seed
  srand48((-1 == opt->seed) ? time(0) : opt->seed);

  // update flow order
  if(IONTORRENT == opt->data_type && NULL != opt->flow_order) {
      opt->flow_order_len = strlen((char*)opt->flow_order);
      for(i=0;i<opt->flow_order_len;i++) {
          opt->flow_order[i] = nst_nt4_table[opt->flow_order[i]];
      }
  }

  // Open files
  opt->fp_fa =	xopen(argv[optind+0], "r");
  strcpy(fn_fai, argv[optind+0]); strcat(fn_fai, ".fai");
  opt->fp_fai = fopen(fn_fai, "r"); // NB: depends on returning NULL;
  strcpy(fn_tmp, argv[optind+1]); strcat(fn_tmp, ".mutations.txt");
  opt->fp_mut = xopen(fn_tmp, "w");
  strcpy(fn_tmp, argv[optind+1]); strcat(fn_tmp, ".bfast.fastq");
  opt->fp_bfast = xopen(fn_tmp, "w");
  strcpy(fn_tmp, argv[optind+1]); strcat(fn_tmp, ".bwa.read1.fastq");
  opt->fp_bwa1 = xopen(fn_tmp, "w");
  strcpy(fn_tmp, argv[optind+1]); strcat(fn_tmp, ".bwa.read2.fastq");
  opt->fp_bwa2 = xopen(fn_tmp, "w");

  // Run simulation
  dwgsim_core(opt);

  // Close files
  fclose(opt->fp_fa); fclose(opt->fp_bfast); fclose(opt->fp_bwa1); fclose(opt->fp_bwa2); 
  if(NULL != opt->fp_fai) fclose(opt->fp_fai);

  dwgsim_opt_destroy(opt);

  return 0;
}
