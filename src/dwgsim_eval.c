#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <unistd.h>
#include <float.h>
#include <sys/resource.h>
#include <math.h>
#include "samtools/bam.h"
#include "samtools/sam.h"
#include "dwgsim_eval.h"

#define __IS_TRUE(_val) ((_val == 1) ? "True" : "False")

int 
print_usage(dwgsim_eval_args_t *args)
{
  fprintf(stderr, "\n");
  fprintf(stderr, "Program: dwgsim_eval (short read simulation evaluator)\n");
  fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
  fprintf(stderr, "Contact: Nils Homer <dnaa-help@lists.sourceforge.net>\n\n");
  fprintf(stderr, "Usage: dwgsim_eval [options] <in.sam/in.bam>\n\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "\t-a\tINT\tsplit by [%d]:\n", args->a);
  fprintf(stderr, "\t\t\t\t\t0: by mapping quality\n");
  fprintf(stderr, "\t\t\t\t\t1: by alignment score\n");
  fprintf(stderr, "\t\t\t\t\t2: by suboptimal alignment score\n");
  fprintf(stderr, "\t\t\t\t\t3: by alignment score - suboptimal alignment score\n");
  fprintf(stderr, "\t-b\t\talignments are from BWA (for SOLiD data only) [%s]\n", __IS_TRUE(args->b));
  fprintf(stderr, "\t-c\t\tcolor space alignments [%s]\n", __IS_TRUE(args->c));
  fprintf(stderr, "\t-d\tINT\tdivide quality/alignment score by this factor [%d]\n", args->d);
  fprintf(stderr, "\t-g\t\tgap \"wiggle\" [%d]\n", args->g);
  fprintf(stderr, "\t-m\t\tconsecutive alignments with the same name (and end for multi-ends) should be treated as multi-mapped reads [%s]\n", __IS_TRUE(args->m));
  fprintf(stderr, "\t-n\tINT\tnumber of raw input paired-end reads (otherwise, inferred from all SAM records present) [%d]\n", args->n);
  fprintf(stderr, "\t-q\tINT\tconsider only alignments with this mapping quality or greater [%d]\n", args->q);
  fprintf(stderr, "\t-z\t\tinput contains only single end reads [%s]\n", __IS_TRUE(args->z));
  fprintf(stderr, "\t-S\t\tinput is SAM [%s]\n", __IS_TRUE(args->S));
  fprintf(stderr, "\t-p\t\tprint incorrect alignments [%s]\n", __IS_TRUE(args->p));
  fprintf(stderr, "\t-s\tINT\tconsider only alignments with the number of specified SNPs [%d]\n", args->s);
  fprintf(stderr, "\t-e\tINT\tconsider only alignments with the number of specified errors [%d]\n", args->e);
  fprintf(stderr, "\t-i\t\tconsider only alignments with indels [%s]\n", __IS_TRUE(args->i));
  fprintf(stderr, "\t-P\tSTRING\ta read prefix that was prepended to each read name [%s]\n", (NULL == args->P) ? "not using" : args->P);
  fprintf(stderr, "\t-h\t\tprint this help message\n");
  return 1;
}

/* Action */
enum {Exit, Warn, LastActionType};
/* Type */
enum {  
    Dummy,
    OutOfRange, /* e.g. command line args */
    InputArguments, 
    IllegalFileName,   
    IllegalPath,
    OpenFileError,
    EndOfFile,
    ReallocMemory,
    MallocMemory,
    ThreadError,
    ReadFileError,
    WriteFileError,
    DeleteFileError,
    LastErrorType,
};       
#define BREAK_LINE "************************************************************\n"
void
dwgsim_eval_print_error(char* FunctionName, char *VariableName, char* Message, int Action, int type)
{
  static char ErrorString[][20]=
    { "\0", "OutOfRange", "InputArguments", "IllegalFileName", "IllegalPath", "OpenFileError", "EndOfFile", "ReallocMemory", "MallocMemory", "ThreadError", "ReadFileError", "WriteFileError", "DeleteFileError"};
  static char ActionType[][20]={"Fatal Error", "Warning"};
  fprintf(stderr, "%s\rIn function \"%s\": %s[%s]. ",
          BREAK_LINE, FunctionName, ActionType[Action], ErrorString[type]);

  /* Only print variable name if is available */
  if(VariableName) {
      fprintf(stderr, "Variable/Value: %s.\n", VariableName);
  }
  /* Only print message name if is available */
  if(Message) {
      fprintf(stderr, "Message: %s.\n", Message);
  }
  if(type == ReadFileError ||
     type == OpenFileError ||
     type == WriteFileError) {
      perror("The file stream error was:");
  }

  switch(Action) {
    case Exit:
      fprintf(stderr, " ***** Exiting due to errors *****\n");
      fprintf(stderr, "%s", BREAK_LINE);
      exit(EXIT_FAILURE);
      break; /* Not necessary actually! */
    case Warn:
      fprintf(stderr, " ***** Warning *****\n");
      fprintf(stderr, "%s", BREAK_LINE);
      break;
    default:
      fprintf(stderr, "Trouble!!!\n");
      fprintf(stderr, "%s", BREAK_LINE);
  }
}


int 
main(int argc, char *argv[])
{
  char c;
  dwgsim_eval_args_t args;

  args.a = args.b = args.c = args.i = args.m = args.n = args.p = args.q = args.z = 0; 
  args.d = 1;
  args.e = -1;
  args.g = 5;
  args.s = -1;
  args.S = 0;
  args.P = NULL;

  while(0 <= (c = getopt(argc, argv, "a:d:e:g:m:n:q:s:bchimpzSP:"))) {
      switch(c) {
        case 'a': args.a = atoi(optarg); break;
        case 'b': args.b = 1; break;
        case 'c': args.c = 1; break;
        case 'd': args.d = atoi(optarg); break;
        case 'g': args.g = atoi(optarg); break;
        case 'm': args.m = 1; break;
        case 'h': return print_usage(&args); break;
        case 'n': args.n = atoi(optarg); break;
        case 'q': args.q = atoi(optarg); break;
        case 'z': args.z = 1; break;
        case 'S': args.S = 1; break;
        case 'p': args.p = 1; break;
        case 's': args.s = atoi(optarg); break;
        case 'e': args.e = atoi(optarg); break;
        case 'i': args.i = 1; break;
        case 'P': free(args.P); args.P = strdup(optarg); break;
        default: fprintf(stderr, "Unrecognized option: -%c\n", c); return 1;
      }
  }

  if(argc == optind) {
      return print_usage(&args);
  }

  run(&args, argc - optind, argv + optind);

  free(args.P);

  return 0;
}

void 
run(dwgsim_eval_args_t *args,
    int32_t num_files,
    char *files[]) 
{
  char *FnName="run";
  int32_t i, n = 0;
  samfile_t *fp_in = NULL;
  samfile_t *fp_out = NULL;
  bam1_t *b=NULL;
  dwgsim_eval_counts_t *counts;
  char *prev_qname=NULL;
  int32_t prev_end=-1;

  // initialize counts
  counts = dwgsim_eval_counts_init();

  fprintf(stderr, "Analyzing...\nCurrently on:\n0");
  for(i=0;i<num_files;i++) {
      // Open the file
      fp_in = samopen(files[i], (1 == args->S) ? "r" : "rb", 0); 
      if(NULL == fp_in) {
          dwgsim_eval_print_error(FnName, files[i], "Could not open file for reading", Exit, OpenFileError);
      }

      if(0 == i && 1 == args->p) {
          fp_out = samopen("-", "wh", fp_in->header);
          if(NULL == fp_out) {
              dwgsim_eval_print_error(FnName, "stdout", "Could not open file stream for writing", Exit, OpenFileError);
          }
      }

      b = bam_init1();
      while(0 < samread(fp_in, b)) {
          if(1 == args->m &&
             NULL != prev_qname && 
             prev_end == (BAM_FREAD1 & b->core.flag) &&
             0 == strcmp(prev_qname, bam1_qname(b))) {
              // do nothing
          }
          else {
              if(1 == args->m) {
                  free(prev_qname);
                  prev_qname = strdup(bam1_qname(b));
                  prev_end = (BAM_FREAD1 & b->core.flag);
              }

              process_bam(counts, args, fp_in->header, b, fp_out);

              if((BAM_FPAIRED & b->core.flag)) { // paired end
                  if(1 == args->z) { // expect single end
                      dwgsim_eval_print_error(FnName, NULL, "Found a read that was paired end", Exit, OutOfRange);
                  }
                  if((BAM_FREAD1 & b->core.flag)) { // count # of pairs
                      n++;
                  }
              }
              else { // single end
                  if(0 == args->z) { // expect paired end
                      dwgsim_eval_print_error(FnName, NULL, "Found a read that was not paired", Exit, OutOfRange);
                  }
                  n++;
              }


              if(0 == (n % 10000)) {
                  fprintf(stderr, "\r%lld", (long long int)n);
              }
          }

          bam_destroy1(b);
          b = bam_init1();
      }
      bam_destroy1(b);

      // Close the file
      samclose(fp_in);
  }
  free(prev_qname);

  if(1 == args->p) samclose(fp_out);

  fprintf(stderr, "\r%lld\n", (long long int)n);

  if(0 < args->n) {
      if(n != args->n) {
          fprintf(stderr, "(-n)=%d\tn=%d\n", args->n, n);
          dwgsim_eval_print_error(FnName, NULL, "Number of reads found differs from the number specified (-n)", Warn, OutOfRange);
      }
  }
  if(0 == args->z) {
      dwgsim_eval_counts_print(counts, args->a, args->d, (0 < args->n) ? 2*args->n : 2*n);
  }
  else {
      dwgsim_eval_counts_print(counts, args->a, args->d, (0 < args->n) ? args->n : n);
  }

  dwgsim_eval_counts_destroy(counts);

  fprintf(stderr, "Analysis complete.\n");
}

uint32_t bam_calclip(const bam1_t *b)
{
  const bam1_core_t *c = &b->core;
  const uint32_t *cigar = bam1_cigar(b);
  uint32_t k, end = 0;
  for (k = 0; k < c->n_cigar; ++k) {
      int op = cigar[k] & BAM_CIGAR_MASK;
      if(op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) {
        end += cigar[k] >> BAM_CIGAR_SHIFT;
      }
      else {
        break;
      }
  }
  return end;
}


void
process_bam(dwgsim_eval_counts_t *counts,
            dwgsim_eval_args_t *args,
            bam_header_t *header,
            bam1_t *b,
            samfile_t *fp_out)
{
  char *FnName="process_bam";
  int32_t left, metric=INT_MIN;
  char *chr=NULL;
  char *name=NULL, *ptr=NULL;
  char chr_name[1028]="\0";
  char read_num[1028]="\0";
  int32_t pos_1, pos_2, str_1, str_2, rand_1, rand2; 
  int32_t n_err_1, n_sub_1, n_indel_1, n_err_2, n_sub_2, n_indel_2;
  int32_t pos, str, rand;
  int32_t i, j, tmp;
  int32_t predicted_value, actual_value;
  int32_t clip;

  // mapping quality threshold
  if(b->core.qual < args->q) return;

  // parse read name
  name = strdup(bam1_qname(b));
  ptr = name; // save to be freed
  char *to_rm="_::_::_______"; // to remove
  for(i=b->core.l_qname-1,j=0;0<=i && j<13;i--) { // replace with spaces 
      if(name[i] == to_rm[j]) {
          name[i] = ' '; j++; 
      }
  }
  // check for the prefix
  if(NULL != args->P) {
      j = strlen(name);
      tmp = strlen(args->P);
      if(j < tmp || 0 != strncmp(args->P, name, tmp)) {
          dwgsim_eval_print_error(FnName, name, "[dwgsim_eval] could not match read name with given read name prefix (-P)", Exit, OutOfRange);
          free(ptr);
          return;
      }
      name += tmp + 1;
  }
  if(14 != sscanf(name, "%s %d %d %1d %1d %1d %1d %d %d %d %d %d %d %s",
                  chr_name, &pos_1, &pos_2, &str_1, &str_2, &rand_1, &rand2,
                  &n_err_1, &n_sub_1, &n_indel_1,
                  &n_err_2, &n_sub_2, &n_indel_2,
                  read_num)) {
      dwgsim_eval_print_error(FnName, name, "[dwgsim_eval] read was not generated by dwgsim?", Exit, OutOfRange);
      free(ptr);
      return;
  }
  // check for a prefix, and make sure it was removed correctly
  if(1 == args->z || (b->core.flag & BAM_FREAD1)) {
      rand = rand_1;
  }
  else {
      rand = rand2;
  }
  if(0 == rand) {
      for(j=0;j<header->n_targets;j++) {
          i = strlen(name);
          tmp = strlen(header->target_name[j]); 
          i = (i < tmp) ? i : tmp;
          if(0 == strncmp(name, header->target_name[j], i)) {
              break;
          }
      }
      if(j == header->n_targets) {
          dwgsim_eval_print_error(FnName, name, "[dwgsim_eval] the mapped contig does not exist in the SAM header; perhaps you have a read name prefix?", Exit, OutOfRange);
      }
  }
  free(ptr);
  ptr = name = NULL;

  // get metric value
  if(0 == args->a) {
      metric = (b->core.qual / args->d); 
      if(DWGSIM_EVAL_MAXQ < metric) metric = DWGSIM_EVAL_MAXQ;
  }
  else if((BAM_FUNMAP & b->core.flag) || 0 == b->core.qual) { // unmapped or zero quality
      metric = DWGSIM_EVAL_MINAS;
  }
  else {
      uint8_t *aux_AS=NULL, *aux_XS=NULL;
      metric = DWGSIM_EVAL_MINAS+1;
      if(1 == args->a || 3 == args->a) {
          aux_AS = bam_aux_get(b, "AS");
          if(NULL == aux_AS) {
              metric = DWGSIM_EVAL_MINAS;
          }
      }
      if(2 == args->a || 3 == args->a) {
          aux_XS = bam_aux_get(b, "XS");
          if(NULL == aux_XS) {
              metric = DWGSIM_EVAL_MINAS;
          }
      }
      if(metric != DWGSIM_EVAL_MINAS) {
          switch(args->a) {
            case 1:
              metric = bam_aux2i(aux_AS);
              break;
            case 2:
              metric = bam_aux2i(aux_XS);
              break;
            case 3:
              metric = bam_aux2i(aux_AS) - bam_aux2i(aux_XS);
              break;
            default:
              metric = DWGSIM_EVAL_MINAS;
              break;
          }
      }
  }
  metric /= args->d;
  if(metric < DWGSIM_EVAL_MINAS) metric = DWGSIM_EVAL_MINAS;

  if(1 == args->i) { // indels only
      if(1 == args->z || (b->core.flag & BAM_FREAD1)) {
          if(0 == n_indel_1) return;
      }
      else {
          if(0 == n_indel_2) return;
      }
  }
  else if(0 <= args->e && n_err_1 !=  args->e) { // # of errors
      return;
  }
  else if(0 <= args->s && n_sub_1 !=  args->s) { // # of snps
      return;
  }

  if(1 == args->c && 1 == args->b) { // SOLiD and BWA
      // Swap 1 and 2
      tmp=n_err_1; n_err_1=n_err_2; n_err_2=tmp;
      tmp=n_sub_1; n_sub_1=n_sub_2; n_sub_2=tmp;
      tmp=n_indel_1; n_indel_1=n_indel_2; n_indel_2=tmp;
  }

  // copy data
  if(1 == args->z || (b->core.flag & BAM_FREAD1)) {
      pos = pos_1; str = str_1; rand = rand_1;
  }
  else {
      pos = pos_2; str = str_2; rand = rand2;
  }

  // get the actual value 
  if(1 == rand) {
      actual_value = DWGSIM_EVAL_UNMAPPABLE;
  }
  else {
      actual_value = DWGSIM_EVAL_MAPPABLE;
  }

  // get the predicted value
  if((BAM_FUNMAP & b->core.flag)) { // unmapped
      predicted_value = DWGSIM_EVAL_UNMAPPED;
  }
  else { // mapped (correctly?)
      clip = bam_calclip(b); 
      chr = header->target_name[b->core.tid];
      left = b->core.pos - clip;

      if(1 == rand // should not map 
         || str != bam1_strand(b) // different strand
         || 0 != strcmp(chr, chr_name)  // different chromosome
         || args->g < fabs(pos - left)) { // out of bounds (positionally) 
          predicted_value = DWGSIM_EVAL_MAPPED_INCORRECTLY;
      }
      else {
          predicted_value = DWGSIM_EVAL_MAPPED_CORRECTLY;
      }
  }

  dwgsim_eval_counts_add(counts, metric, actual_value, predicted_value);

  // print incorrect alignments
  if(1 == args->p && DWGSIM_EVAL_MAPPED_INCORRECTLY == predicted_value) {
      if(samwrite(fp_out, b) <= 0) {
          dwgsim_eval_print_error(FnName, "stdout", "Could not write to stream", Exit, WriteFileError);
      }
  }
}

dwgsim_eval_counts_t *
dwgsim_eval_counts_init()
{
  dwgsim_eval_counts_t *counts;

  counts = malloc(sizeof(dwgsim_eval_counts_t));

  counts->min_score = counts->max_score = 0;

  counts->mc = malloc(sizeof(int32_t)); assert(NULL != counts->mc);
  counts->mi = malloc(sizeof(int32_t)); assert(NULL != counts->mi);
  counts->mu = malloc(sizeof(int32_t)); assert(NULL != counts->mu);
  counts->um = malloc(sizeof(int32_t)); assert(NULL != counts->um);
  counts->uu = malloc(sizeof(int32_t)); assert(NULL != counts->uu);

  counts->mc[0] = counts->mi[0] = counts->mu[0] = 0;
  counts->um[0] = counts->uu[0] = 0;

  return counts;
}

void
dwgsim_eval_counts_destroy(dwgsim_eval_counts_t *counts)
{
  free(counts->mc);
  free(counts->mi);
  free(counts->mu);
  free(counts->um);
  free(counts->uu);
  free(counts);
}

void 
dwgsim_eval_counts_add(dwgsim_eval_counts_t *counts, int32_t score, int32_t actual_value, int32_t predicted_value)
{
  char *FnName="dwgsim_eval_counts_add";
  int32_t i, m, n;
  if(counts->max_score < score) {
      m = score - counts->min_score + 1;
      n = counts->max_score - counts->min_score + 1;

      counts->mc = realloc(counts->mc, sizeof(int32_t)*m); assert(NULL != counts->mc);
      counts->mi = realloc(counts->mi, sizeof(int32_t)*m); assert(NULL != counts->mi);
      counts->mu = realloc(counts->mu, sizeof(int32_t)*m); assert(NULL != counts->mu);
      counts->um = realloc(counts->um, sizeof(int32_t)*m); assert(NULL != counts->um);
      counts->uu = realloc(counts->uu, sizeof(int32_t)*m); assert(NULL != counts->uu);

      // initialize to zero
      for(i=n;i<m;i++) {
          counts->mc[i] = counts->mi[i] = counts->mu[i] = 0;
          counts->um[i] = counts->uu[i] = 0;
      }
      counts->max_score = score;
  }
  else if(score < counts->min_score) {
      m = counts->max_score - score + 1;
      n = counts->max_score - counts->min_score + 1;

      counts->mc = realloc(counts->mc, sizeof(int32_t)*m); assert(NULL != counts->mc);
      counts->mi = realloc(counts->mi, sizeof(int32_t)*m); assert(NULL != counts->mi);
      counts->mu = realloc(counts->mu, sizeof(int32_t)*m); assert(NULL != counts->mu);
      counts->um = realloc(counts->um, sizeof(int32_t)*m); assert(NULL != counts->um);
      counts->uu = realloc(counts->uu, sizeof(int32_t)*m); assert(NULL != counts->uu);

      // shift up
      for(i=m-1;m-n<=i;i--) {
          counts->mc[i] = counts->mc[i-(m-n)]; 
          counts->mi[i] = counts->mi[i-(m-n)]; 
          counts->mu[i] = counts->mu[i-(m-n)]; 
          counts->um[i] = counts->um[i-(m-n)]; 
          counts->uu[i] = counts->uu[i-(m-n)]; 
      }
      // initialize to zero
      for(i=0;i<m-n;i++) {
          counts->mc[i] = counts->mi[i] = counts->mu[i] = 0;
          counts->um[i] = counts->uu[i] = 0;
      }
      counts->min_score = score;
  }

  // check actual value
  switch(actual_value) {
    case DWGSIM_EVAL_MAPPABLE:
    case DWGSIM_EVAL_UNMAPPABLE:
      break;
    default:
      dwgsim_eval_print_error(FnName, "actual_value", "Could not understand actual value", Exit, OutOfRange);
  }

  // check predicted value
  switch(predicted_value) {
    case DWGSIM_EVAL_MAPPED_CORRECTLY:
    case DWGSIM_EVAL_MAPPED_INCORRECTLY:
    case DWGSIM_EVAL_UNMAPPED:
      break;
    default:
      dwgsim_eval_print_error(FnName, "predicted_value", "Could not understand predicted value", Exit, OutOfRange);
  }

  switch(actual_value) {
    case DWGSIM_EVAL_MAPPABLE:
      switch(predicted_value) {
        case DWGSIM_EVAL_MAPPED_CORRECTLY:
          counts->mc[score-counts->min_score]++; break;
        case DWGSIM_EVAL_MAPPED_INCORRECTLY:
          counts->mi[score-counts->min_score]++; break;
        case DWGSIM_EVAL_UNMAPPED:
          counts->mu[score-counts->min_score]++; break;
        default:
          break; // should not reach here
      }
      break;
    case DWGSIM_EVAL_UNMAPPABLE:
      switch(predicted_value) {
        case DWGSIM_EVAL_MAPPED_CORRECTLY:
          dwgsim_eval_print_error(FnName, "predicted_value", "predicted value cannot be mapped correctly when the read is unmappable", Exit, OutOfRange); break;
        case DWGSIM_EVAL_MAPPED_INCORRECTLY:
          counts->um[score-counts->min_score]++; break;
        case DWGSIM_EVAL_UNMAPPED:
          counts->uu[score-counts->min_score]++; break;
        default:
          break; // should not reach here
      }
      break;
    default:
      break; // should not reach here
  }
}

void 
dwgsim_eval_counts_print(dwgsim_eval_counts_t *counts, int32_t a, int32_t d, int32_t n)
{
  int32_t i;
  int32_t max = 0;
  int32_t mc_sum, mi_sum, mu_sum, um_sum, uu_sum;
  int32_t m_total, mm_total, u_total;
  char format[1024]="\0";

  mc_sum = mi_sum = mu_sum = um_sum = uu_sum = 0;
  m_total = mm_total = u_total = 0;

  // create the format string
  for(i=counts->max_score - counts->min_score;0<=i;i--) {
      m_total += counts->mc[i] + counts->mi[i] + counts->mu[i];
      u_total += counts->um[i] + counts->uu[i];
      max += counts->mc[i] + counts->mi[i] + counts->mu[i] + counts->um[i] + counts->uu[i];
  }
  max = 1 + log10(max);
  strcat(format, "%.2d ");
  for(i=0;i<12;i++) {
      sprintf(format + (int)strlen(format), "%%%dd ", max);
  }
  strcat(format + (int)strlen(format), "%.3e %.3e %.3e %.3e %.3e %.3e\n");

  // header
  fprintf(stdout, "# thr | the minimum %s threshold\n", (0 == a) ? "mapping quality" : "alignment score");

  fprintf(stdout, "# mc | the number of correctly mapped reads that should be mapped at the threshold\n");
  fprintf(stdout, "# mi | the number of incorrectly mapped reads that should be mapped at the threshold\n");
  fprintf(stdout, "# mu | the number of unmapped reads that should be mapped at the threshold\n");
          
  fprintf(stdout, "# um | the number of mapped reads that should be unmapped at the threshold\n");
  fprintf(stdout, "# uu | the number of unmapped reads that should be unmapped at the threshold\n");

  fprintf(stdout, "# mc + mi + mu + um + uu | the total number of reads at the threshold\n");

  fprintf(stdout, "# mc' | the number of correctly mapped reads that should be mapped at or greater than that threshold\n");
  fprintf(stdout, "# mi' | the number of incorrectly mapped reads that should be mapped at or greater than that threshold\n");
  fprintf(stdout, "# mu' | the number of unmapped reads that should be mapped at or greater than that threshold\n");
          
  fprintf(stdout, "# um' | the number of mapped reads that should be unmapped at or greater than that threshold\n");
  fprintf(stdout, "# uu' | the number of unmapped reads that should be unmapped at or greater than that threshold\n");

  fprintf(stdout, "# mc' + mi' + mu' + um' + uu' | the total number of reads at or greater than the threshold\n");
          
  fprintf(stdout, "# (mc / (mc' + mi' + mu')) | sensitivity: the fraction of mappable reads that are mapped correctly at the threshold\n");
  fprintf(stdout, "# (mc / (mc' + mi')) | positive predictive value: the fraction of mapped mappable reads that are mapped correctly at the threshold\n");
  fprintf(stdout, "# (um / (um' + uu')) | false discovery rate: the fraction of random reads that are mapped at the threshold\n");
  fprintf(stdout, "# (mc' / (mc' + mi' + mu')) | sensitivity: the fraction of mappable reads that are mapped correctly at or greater than the threshold\n");
  fprintf(stdout, "# (mc' / (mc' + mi')) | positive predictive value: the fraction of mapped mappable reads that are mapped correctly at or greater than the threshold\n");
  fprintf(stdout, "# (um' / (um' + uu')) | false discovery rate: the fraction of random reads that are mapped at or greater than the threshold\n");


  // print
  for(i=counts->max_score - counts->min_score;0<=i;i--) {
      double num, den;

      mc_sum += counts->mc[i];
      mi_sum += counts->mi[i];
      mu_sum += counts->mu[i];
      um_sum += counts->um[i];
      uu_sum += counts->uu[i];
      mm_total += counts->mc[i] + counts->mi[i];

      /* Notes:
       *  notice that the denominator for sensitivity (and fdr) for the "ge" 
       *  (greater than or equal) threshold is the "total", while the denominator 
       *  for ppv is "@ >= Q".  The reasoning behind this is ppv is a measure of
       *  the quality of mappings that will be returned when using a Q threshold
       *  while the sensitivity want to measure the fraction of mappings that
       *  will be returned compared to the maximum.  Basically, if we accept
       *  only mappings at a given threshold, and call the rest unmapped, what
       *  happens?
       *  - sensitivity tells us the # of correct mappings out of the total possible
       *  mappings.
       *  - ppv tells us the # of correct mappings out of the total mappings.
       *  - fdr tells us the # of random mappings out of the total unmappable.
       */
       
      // "at" sensitivity: mapped correctly @ Q / mappable @ Q
      num = counts->mc[i];
      den = counts->mc[i] + counts->mi[i] + counts->mu[i];
      double sens_at_thr = (0 == den) ? 0. : (num / (double)den);  
      // "ge" sensitivity: mapped correctly @ >= Q / total mappable
      double sens_ge_thr = (0 == m_total) ? 0. : (mc_sum / (double)m_total);
      
      // "at" positive predictive value: mapped correctly @ Q / mappable and mapped @ Q
      num = counts->mc[i];
      den = counts->mc[i] + counts->mi[i];
      double ppv_at_thr = (0 == den) ? 0. : (num / (double)den);
      // "ge" positive predictive value: mapped correctly @ >= Q / mappable and mapped @ >= Q
      double ppv_ge_thr = (0 == mm_total) ? 0. : (mc_sum / (double)mm_total);

      // "at" false discovery rate: unmappable and mapped @ Q / unmappable @ Q
      num = counts->um[i];
      den = counts->um[i] + counts->uu[i];
      double fdr_at_thr = (0 == den) ? 0. : (num / (double)den);
      // "ge" false discovery rate: unmappable and mapped @ >= Q / unmappable @ >= Q
      double fdr_ge_thr = (0 == u_total) ? 0. : (um_sum / (double)u_total);
      
      fprintf(stdout, format,
              (i + counts->min_score)*d,
              counts->mc[i], counts->mi[i], counts->mu[i], counts->um[i], counts->uu[i],
              counts->mc[i] + counts->mi[i] + counts->mu[i] + counts->um[i] + counts->uu[i],
              mc_sum, mi_sum, mu_sum, um_sum, uu_sum,
              mc_sum + mi_sum + mu_sum + um_sum + uu_sum,
              sens_at_thr, ppv_at_thr, fdr_at_thr,
              sens_ge_thr, ppv_ge_thr, fdr_ge_thr);
  }
}
