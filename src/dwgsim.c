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
#include "contigs.h"
#include "mut.h"
#include "mut_txt.h"
#include "mut_bed.h"
#include "regions_bed.h"
#include "dwgsim_opt.h"
#include "dwgsim.h"
//#include <config.h>

#define QUAL_MAX 40

uint8_t nst_nt4_table[256] = {
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
        mut_t c = currseq->s[i], mut_type = c & mutmsk;			\
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
            if(1 == mut_get_ins(currseq, i, &n, &ins)) { \
                if(0 == strand[x]) { \
                    while(n > 0 && k < s[x]) { \
                        tmp_seq[x][k++] = ins & 0x3;                \
                        --n, ins >>= 2; \
                    } \
                    if(k < s[x]) tmp_seq[x][k++] = c & 0xf;						\
                } else { \
                    tmp_seq[x][k++] = c & 0xf;						\
                    while(n > 0 && k < s[x]) { \
                        ext_coor[x]++; \
                        tmp_seq[x][k++] = (ins >> ((n-1) << 1) & 0x3);                \
                        --n; \
                    } \
                } \
            } else { \
                int32_t byte_index, bit_index; \
                uint32_t num_ins; \
                uint8_t *insertion = NULL; \
                insertion = mut_get_ins_long_n(currseq->ins[ins], &num_ins); \
                if(0 == strand[x]) { \
                    byte_index = mut_packed_len(num_ins) - 1; bit_index = 3 - (num_ins & 3); \
                    while(num_ins > 0 && k < s[x]) { \
                        assert(0 <= byte_index); \
                        tmp_seq[x][k++] = (insertion[byte_index] >> (bit_index << 1)) & 0x3;                \
                        --num_ins; \
                        bit_index--; \
                        if (bit_index < 0) { \
                            bit_index = 3; \
                            byte_index--; \
                        } \
                    } \
                    if(k < s[x]) tmp_seq[x][k++] = c & 0xf;						\
                } else { \
                    tmp_seq[x][k++] = c & 0xf;						\
                    byte_index = 0; bit_index = 0; \
                    while(num_ins > 0 && k < s[x]) { \
                        ext_coor[x]++; \
                        tmp_seq[x][k++] = (insertion[byte_index] >> (bit_index << 1)) & 0x3;                \
                        --num_ins; \
                        bit_index++; \
                        if (4 == bit_index) { \
                            bit_index = 0; \
                            byte_index++; \
                        } \
                    } \
                } \
            }													\
        } \
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

char iupac_and_base_to_mut(char iupac, char base)
{
  static char *codes = "XACMGRSVTWYHKDBN";
  int32_t i;
  int32_t b = nst_nt4_table[(int)base];
  for(i=0;i<4;i++) { 
      if(codes[1<<(b&0x3) | 1<<(i&0x3)] == iupac) {
          return "ACGTN"[i];
      }
  }
  return 'X';
}

char bases_to_iupac(char b1, char b2)
{
  int32_t a1, a2, n;
  a1 = nst_nt4_table[(int)b1];
  a2 = nst_nt4_table[(int)b2];
  if(4 == a1 || 4 == a2) return 'X';
  if(a1 == b1) return 'X';
  if(a1 < a2) {
      n = a1 + (a2 << 2);
  }
  else {
      n = a2 + (a1 << 2);
  }
  switch(n) {
    case 4: return 'M'; // 0 + 4*1 = M
    case 8: return 'R'; // 0 + 4*2 = R
    case 12: return 'W'; // 0 + 4*3 = W
    case 9: return 'S'; // 1 + 4*2 = S
    case 13: return 'Y'; // 1 + 4*3 = Y
    case 14: return 'K'; // 2 + 4*3 = K
    default: break;
  }
  return 'X';
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

int32_t
generate_errors_flows(dwgsim_opt_t *opt, uint8_t **seq, uint8_t **mask, int32_t *mem, int32_t len, uint8_t strand, double e, int32_t *_n_err)
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
      (*mask)[i] = 0;
  }
  assert(opt->flow_order_len != i);
  flow_i = i;
  prev_c = 4;
  for(i=0;i<len;i++) {
      c = (4 <= (*seq)[i]) ? 0 : (*seq)[i];
      while(c != opt->flow_order[flow_i]) { // skip paste empty flows
          (*mask)[flow_i] = 0;
          flow_i = (flow_i + 1) % opt->flow_order_len;
      }
      if(prev_c != c) { // new hp
          (*mask)[flow_i] = 0;
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
                      (*mask) = realloc((*mask), sizeof(uint8_t) * (*mem));
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
                  (*mask)[flow_i] = 1;
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
          if(0 == (*mask)[flow_i] && 0 < n_err) {  // insert
              // more memory
              while((*mem) <= len + n_err) {
                  (*mem) <<= 1; // double
                  (*seq) = realloc((*seq), sizeof(uint8_t) * (*mem));
                  (*mask) = realloc((*mask), sizeof(uint8_t) * (*mem));
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

void dwgsim_core(dwgsim_opt_t * opt)
{
  seq_t seq;
  mutseq_t *mutseq[2]={NULL,NULL};
  uint64_t tot_len, ii=0, ctr=0;
  int i, l, m, n_ref, contig_i;
  char name[1024], *qstr;
  int32_t name_len_max=0;
  int size[2], prev_skip=0, qstr_l=0;
  int num_n[2];
  uint8_t *tmp_seq[2]={NULL,NULL};
  uint8_t *tmp_seq_flow_mask[2]={NULL,NULL};
  int32_t tmp_seq_mem[2]={0,0};
  int64_t n_sim = 0;
  error_t *e[2]={NULL,NULL};
  FILE *fp_muts_input = NULL;
  FILE *fp_regions_bed = NULL;
  muts_input_t *muts_input = NULL;
  regions_bed_txt *regions_bed = NULL;
  contigs_t *contigs = NULL;

  e[0] = &opt->e[0]; e[1] = &opt->e[1];

  INIT_SEQ(seq);
  seq_set_block_size(0x1000000);
  l = opt->length[0] > opt->length[1]? opt->length[0] : opt->length[1];
  qstr_l = l;
  qstr = (char*)calloc(qstr_l+1, 1);
  tmp_seq[0] = (uint8_t*)calloc(l+2, 1);
  tmp_seq[1] = (uint8_t*)calloc(l+2, 1);
  if(IONTORRENT == opt->data_type) {
      tmp_seq_flow_mask[0] = (uint8_t*)calloc(l+2, 1);
      tmp_seq_flow_mask[1] = (uint8_t*)calloc(l+2, 1);
  }
  tmp_seq_mem[0] = tmp_seq_mem[1] = l+2;
  size[0] = opt->length[0]; size[1] = opt->length[1];
  
  if(0 <= opt->fn_muts_input_type) {
      contigs = contigs_init();
  }
  
  if(NULL != opt->fn_regions_bed) {
      fp_regions_bed = xopen(opt->fn_regions_bed, "r");
      if(NULL == contigs) contigs = contigs_init();
  }

  tot_len = n_ref = 0;
  mut_print_header_pre(opt->fp_vcf);
  if(NULL != opt->fp_fai) {
      int dummy_int[3];
      while(0 < fscanf(opt->fp_fai, "%s\t%d\t%d\t%d\t%d", name, &l, &dummy_int[0], &dummy_int[1], &dummy_int[2])) {
          fprintf(stderr, "[dwgsim_core] %s length: %d\n", name, l);
          tot_len += l;
          ++n_ref;
          mut_print_header_contig(opt->fp_vcf, name, l);
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
          mut_print_header_contig(opt->fp_vcf, name, l);
          if(NULL != contigs) {
              contigs_add(contigs, name, l);
          }
      }
  }
  fprintf(stderr, "[dwgsim_core] %d sequences, total length: %llu\n", n_ref, (long long)tot_len);
  rewind(opt->fp_fa);
  mut_print_header_post(opt->fp_vcf);

  if(0 <= opt->fn_muts_input_type) {
      fp_muts_input = xopen(opt->fn_muts_input, "r");
      muts_input = muts_input_init(fp_muts_input, contigs, opt->fn_muts_input_type); // read in the mutation file
  }
  
  if(NULL != opt->fn_regions_bed) {
      regions_bed = regions_bed_init(fp_regions_bed, contigs);
      // recalculate the total length
      tot_len = 0;
      for(i=0;i<regions_bed->n;i++) {
          tot_len += regions_bed->end[i] - regions_bed->start[i] + 1;
      }
  }
  if(NULL != contigs) {
      contigs_destroy(contigs);
      contigs = NULL;
  }

  if(0 == opt->muts_only) {
      fprintf(stderr, "[dwgsim_core] Currently on: \n0");
  }
  else {
      fprintf(stderr, "[dwgsim_core] Currently on:");
  }
  contig_i = 0;
  while ((l = seq_read_fasta(opt->fp_fa, &seq, name, 0)) >= 0) {
      int64_t n_pairs = 0;
      n_ref--;
      
      if(1 == opt->muts_only) {
          fprintf(stderr, "\r[dwgsim_core] Currently on: %s", name);
          if(name_len_max < strlen(name)) {
              name_len_max = strlen(name);
          } 
          else {
              for(i=0;i<name_len_max-strlen(name);i++) {
                  fputc(' ', stderr);
              }
          }
      }

      if(0 == opt->muts_only) {
          if(0 == n_ref && opt->C < 0) {
              n_pairs = opt->N - n_sim;
          }
          else {
              if(NULL != regions_bed) {
                  // recalculate l
                  m = 0;
                  for(i=0;i<regions_bed->n;i++) {
                      if(contig_i == regions_bed->contig[i]) {
                          m += regions_bed->end[i] - regions_bed->start[i] + 1;
                      }
                  }
                  if(0 == m) {
                      fprintf(stderr, "[dwgsim_core] #0 skip sequence '%s' as it is not in the targeted region\n", name);
                      contig_i++;
                      continue; // skip this region
                  }
                  l = m;

                  int num_n = 0;
                  for(i=0;i<regions_bed->n;i++) {           
                      if(contig_i == regions_bed->contig[i]) {                    
                          int m;                                                                            
                          for(m=regions_bed->start[i];m<=regions_bed->end[i];m++) {                                               
                              switch (seq.s[m-1]) {                                                                                                             
                                case 'a':
                                case 'A':
                                case 'c':
                                case 'C':
                                case 'g':
                                case 'G':
                                case 't':
                                case 'T':
                                  break;                                                                                                                                                                         default:     
                                    num_n++;                                                                                                                                                                         break;     
                              }                                                                                                                                                                            }                  
                      }                                                                                                                                                                            } 
                  if(0.95 < num_n / (double)l) { // TODO: arbitrary cutoff
                      fprintf(stderr, "[dwgsim_core] #1 skip sequence '%s' as %d out of %d bases are non-ACGT\n", name, num_n, l);
                      contig_i++;
                      continue;
                  }
              }
              if(0 < opt->N) {
                  // based on -N
                  n_pairs = (uint64_t)((long double)l / tot_len * opt->N + 0.5);
                  if(opt->N - n_sim < n_pairs) n_pairs = opt->N - n_sim; // make sure we don't simulate too many reads
              }
              else {
                  // based on coverage, with added random reads
                  n_pairs = (uint64_t)(l * opt->C / ((long double)(size[0] + size[1])) / (1.0 - opt->rand_read) + 0.5);
              }
          }

          // for paired end/mate pair, make sure we have enough bases in this
          // sequence
          if (0 < opt->length[1] && l < opt->dist + 3 * opt->std_dev) {
              if(0 == prev_skip) fprintf(stderr, "\n");
              prev_skip = 1;
              fprintf(stderr, "[dwgsim_core] #2 skip sequence '%s' as it is shorter than %f!\n", name, opt->dist + 3 * opt->std_dev);
              contig_i++;
              continue;
          }
          else if (l < opt->length[0] || (0 < opt->length[1] && l < opt->length[1])) {
              if(0 == prev_skip) fprintf(stderr, "\n");
              prev_skip = 1;
              fprintf(stderr, "[dwgsim_core] #3 skip sequence '%s' as it is shorter than %d!\n", name, (l < opt->length[0]) ? opt->length[0] : opt->length[1]);
              contig_i++;
              continue;
          }
          else if (n_pairs < 0) { // NB: this should not happen
              // not enough pairs
              continue;
          }
          prev_skip = 0;
      }

      // generate mutations and print them out
      mutseq[0] = mutseq_init(); mutseq[1] = mutseq_init();
      mut_diref(opt, &seq, mutseq[0], mutseq[1], contig_i, muts_input);
      mut_print(name, &seq, mutseq[0], mutseq[1], opt->fp_mut, opt->fp_vcf);

      if(0 == opt->muts_only) {
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

                  if(NULL == regions_bed) {
                      do { // avoid boundary failure
                          if(0 < s[1]) { // paired end/mate pair
                              ran = ran_normal();
                              ran = ran * opt->std_dev + opt->dist;
                              d = (int)(ran + 0.5);
                          }
                          else {
                              d = 0;
                          }
                          pos = (int)((l - d + 1) * drand48());
                      } while (pos < 0 
                               || pos >= seq.l 
                               || pos + d - 1 >= seq.l 
                               || (0 < s[1] && 0 == opt->is_inner && ((0 < s[0] && d <= s[1]) || (d <= s[0] && 0 < s[1]))));
                  } 
                  else {
                      do { // avoid boundary failure
                          if(0 < s[1]) {
                              ran = ran_normal();
                              ran = ran * opt->std_dev + opt->dist;
                              d = (int)(ran + 0.5);
                          }
                          else {
                              d = 0;
                          }
                          pos = (int)((l - d + 1) * drand48());
                          // convert in the bed file
                          for(i=0;i<regions_bed->n;i++) { // TODO: regions are in sorted order... so optimize
                              if(contig_i == regions_bed->contig[i]) {
                                  j = regions_bed->end[i] - regions_bed->start[i] + 1;
                                  if(pos < j) {
                                      pos = regions_bed->start[i] + pos - 1; // zero-based
                                      break;
                                  }
                                  else {
                                      pos -= j;
                                  }
                              }
                          }
                      } while (pos < 0 
                               || pos >= seq.l 
                               || pos + d - 1 >= seq.l 
                               || (0 < s[1] && 0 == opt->is_inner && ((0 < s[0] && d <= s[1]) || (d <= s[0] && 0 < s[1])))
                               || 0 == regions_bed_query(regions_bed, contig_i, pos, pos + s[0] + s[1] + d - 1));
                  }

                  // generate the read sequences
                  mutseq_t *currseq = mutseq[drand48()<opt->mut_freq?0:1]; // haplotype from which the reads are generated
                  n_sub[0] = n_sub[1] = n_indel[0] = n_indel[1] = n_err[0] = n_err[1] = 0;
                  n_sub_first[0] = n_sub_first[1] = n_indel_first[0] = n_indel_first[1] = n_err_first[0] = n_err_first[1] = 0;
                  num_n[0]=num_n[1]=0;

                  // strand
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

                  // generate the reads in base space
                  if(0 < s[1]) { // paired end or mate pair
                      if(strand[0] == strand[1]) { // same strand
                          if(0 == strand[0]) { // + strand
                              /*
                               * 5' E2 -----> .... E1 -----> 3'
                               * 3'           ....           5'
                               */
                              if(0 == opt->is_inner) {
                                  __gen_read(0, pos + d - s[0], ++i); 
                              }
                              else {
                                  __gen_read(0, pos + s[1] + d, ++i); 
                              }
                              __gen_read(1, pos, ++i);
                          }
                          else { // - strand
                              /*
                               * 3'           ....            5'
                               * 5' <----- E1 .... <----- E2  3'
                               */
                              __gen_read(0, pos + s[0], --i);
                              if(0 == opt->is_inner) {
                                  __gen_read(1, pos + d, --i);
                              }
                              else {
                                  __gen_read(1, pos + s[0] + d + s[1], --i);
                              }
                          }
                      }
                      else { // opposite strand
                          if(0 == strand[0]) { // + strand
                              /*
                               * 5' E1 -----> ....           3'
                               * 3'           .... <----- E2 5'
                               */
                              __gen_read(0, pos, ++i);
                              if(0 == opt->is_inner) {
                                  __gen_read(1, pos + d, --i);
                              }
                              else {
                                  __gen_read(1, pos + s[0] + d + s[1], --i);
                              }
                          }
                          else { // - strand
                              /*
                               * 5' E2 -----> ....           3'
                               * 3'           .... <----- E1 5'
                               */
                              if(0 == opt->is_inner) {
                                  __gen_read(0, pos + d, --i);
                              }
                              else {
                                  __gen_read(0, pos + s[1] + d + s[0], --i); 
                              }
                              __gen_read(1, pos, i++);
                          }
                      }
                  }
                  else { // fragment
                      if(0 == strand[0]) {
                          __gen_read(0, pos, ++i); // + strand
                      }
                      else {
                          __gen_read(0, pos + s[0] - 1, --i); // - strand
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
                      s[0] = generate_errors_flows(opt, &tmp_seq[0], &tmp_seq_flow_mask[0], &tmp_seq_mem[0], s[0], strand[0], e[0]->start, &n_err[0]);
                      s[1] = generate_errors_flows(opt, &tmp_seq[1], &tmp_seq_flow_mask[1], &tmp_seq_mem[1], s[1], strand[1], e[1]->start, &n_err[1]);
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
                      if(NULL != opt->fixed_quality) {
                          for (i = 0; i < s[j]; ++i) {
                              qstr[i] = opt->fixed_quality[0];
                          }
                      }
                      else {
                          for (i = 0; i < s[j]; ++i) {
                              if (e[j]->start+e[j]->by*i>0) {
                                  qstr[i] = (int)(-10.0 * log(e[j]->start + e[j]->by*i) / log(10.0) + 0.499) + '!';
                              } else {
                                  qstr[i] = QUAL_MAX + '!';
                              }
                              if(0 < opt->quality_std) {
                                  qstr[i] += (int)((ran_normal() * opt->quality_std) + 0.5);
                              }
                              if(qstr[i] < '!') qstr[i] = '!';
                              if(QUAL_MAX + '!' < qstr[i]) qstr[i] = QUAL_MAX + '!';
                          }
                      }
                      qstr[i] = 0;
                      // BWA
                      FILE *fpo = (0 == j) ? opt->fp_bwa1: opt->fp_bwa2;
                      if(ILLUMINA == opt->data_type || IONTORRENT == opt->data_type) {
                          fprintf(fpo, "@%s%s%s_%u_%u_%1u_%1u_%1u_%1u_%d:%d:%d_%d:%d:%d_%llx/%d\n", 
                                  (NULL == opt->read_prefix) ? "" : opt->read_prefix,
                                  (NULL == opt->read_prefix) ? "" : "_",
                                  name, ext_coor[0]+1, ext_coor[1]+1, strand[0], strand[1], 0, 0,
                                  n_err[0], n_sub[0], n_indel[0],
                                  n_err[1], n_sub[1],n_indel[1],
                                  (long long)ii, j+1);
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
                          fprintf(fpo, "@%s%s%s_%u_%u_%1u_%1u_%1u_%1u_%d:%d:%d_%d:%d:%d_%llx/%d\n", 
                                  (NULL == opt->read_prefix) ? "" : opt->read_prefix,
                                  (NULL == opt->read_prefix) ? "" : "_",
                                  name, ext_coor[0]+1, ext_coor[1]+1, strand[0], strand[1], 0, 0,
                                  n_err[0] - n_err_first[0], n_sub[0] - n_sub_first[0], n_indel[0] - n_indel_first[0], 
                                  n_err[1] - n_err_first[1], n_sub[1] - n_sub_first[1], n_indel[1] - n_indel_first[1],
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
                      fprintf(opt->fp_bfast, "@%s%s%s_%u_%u_%1u_%1u_%1u_%1u_%d:%d:%d_%d:%d:%d_%llx\n", 
                              (NULL == opt->read_prefix) ? "" : opt->read_prefix,
                              (NULL == opt->read_prefix) ? "" : "_",
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
                      for (i = 0; i < s[j]; ++i) {
                          tmp_seq[j][i] = (int)(drand48() * 4.0) & 3;
                      }
                      if(NULL != opt->fixed_quality) {
                          for (i = 0; i < s[j]; ++i) {
                              qstr[i] = opt->fixed_quality[0];
                          }
                      }
                      else {
                          for (i = 0; i < s[j]; ++i) {
                              if (e[j]->start+e[j]->by*i>0) {
                                  qstr[i] = (int)(-10.0 * log(e[j]->start + e[j]->by*i) / log(10.0) + 0.499) + '!';
                              } else {
                                  qstr[i] = QUAL_MAX + '!';
                              }
                              if(0 < opt->quality_std) {
                                  qstr[i] += (int)((ran_normal() * opt->quality_std) + 0.5);
                              }
                              if(qstr[i] < '!') qstr[i] = '!';
                              if(QUAL_MAX + '!' < qstr[i]) qstr[i] = '!';
                          }
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
                          fprintf(fpo, "@%s%s%s_%u_%u_%1u_%1u_%1u_%1u_%d:%d:%d_%d:%d:%d_%llx/%d\n", 
                                  (NULL == opt->read_prefix) ? "" : opt->read_prefix,
                                  (NULL == opt->read_prefix) ? "" : "_",
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
                          fprintf(fpo, "@%s%s%s_%u_%u_%1u_%1u_%1u_%1u_%d:%d:%d_%d:%d:%d_%llx/%d\n", 
                                  (NULL == opt->read_prefix) ? "" : opt->read_prefix,
                                  (NULL == opt->read_prefix) ? "" : "_",
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
                      fprintf(opt->fp_bfast, "@%s%s%s_%u_%u_%1u_%1u_%1u_%1u_%d:%d:%d_%d:%d:%d_%llx\n", 
                              (NULL == opt->read_prefix) ? "" : opt->read_prefix,
                              (NULL == opt->read_prefix) ? "" : "_",
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
          fprintf(stderr, "\r[dwgsim_core] %llu",
                  (unsigned long long int)ctr);
      }
      mutseq_destroy(mutseq[0]);
      mutseq_destroy(mutseq[1]);
      contig_i++;
  }
  fprintf(stderr, "\n[dwgsim_core] Complete!\n");

  free(seq.s); free(qstr);
  free(tmp_seq[0]); free(tmp_seq[1]);
  if(IONTORRENT == opt->data_type) {
      free(tmp_seq_flow_mask[0]); free(tmp_seq_flow_mask[1]);
  }
  if(0 <= opt->fn_muts_input_type) {
      muts_input_destroy(muts_input);
  }
  if(NULL != opt->fn_regions_bed) {
      regions_bed_destroy(regions_bed);
      fclose(fp_regions_bed);
  }
}

int main(int argc, char *argv[])
{
  dwgsim_opt_t *opt = NULL;

  // update the mutant sequence bounds
  mutseq_init_bounds();

  opt = dwgsim_opt_init();

  char fn_fai[1024]="\0";
  char fn_tmp[1024]="\0";

  if(0 == dwgsim_opt_parse(opt, argc, argv)) {
      return dwgsim_opt_usage(opt);
  }

  // Open files
  opt->fp_fa =	xopen(argv[optind+0], "r");
  strcpy(fn_fai, argv[optind+0]); strcat(fn_fai, ".fai");
  opt->fp_fai = fopen(fn_fai, "r"); // NB: depends on returning NULL;
  strcpy(fn_tmp, argv[optind+1]); strcat(fn_tmp, ".mutations.txt");
  opt->fp_mut = xopen(fn_tmp, "w");
  strcpy(fn_tmp, argv[optind+1]); strcat(fn_tmp, ".mutations.vcf");
  opt->fp_vcf = xopen(fn_tmp, "w");
  if(0 == opt->muts_only) {
      strcpy(fn_tmp, argv[optind+1]); strcat(fn_tmp, ".bfast.fastq");
      opt->fp_bfast = xopen(fn_tmp, "w");
      strcpy(fn_tmp, argv[optind+1]); strcat(fn_tmp, ".bwa.read1.fastq");
      opt->fp_bwa1 = xopen(fn_tmp, "w");
      strcpy(fn_tmp, argv[optind+1]); strcat(fn_tmp, ".bwa.read2.fastq");
      opt->fp_bwa2 = xopen(fn_tmp, "w");
  }

  // Run simulation
  dwgsim_core(opt);

  // Close files
  if(0 == opt->muts_only) {
      fclose(opt->fp_fa); fclose(opt->fp_bfast); fclose(opt->fp_bwa1); fclose(opt->fp_bwa2); 
  }
  if(NULL != opt->fp_fai) fclose(opt->fp_fai);
  fclose(opt->fp_mut);
  fclose(opt->fp_vcf);

  dwgsim_opt_destroy(opt);

  return 0;
}
