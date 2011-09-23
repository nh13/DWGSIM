#ifndef CONTIGS_H_
#define CONTIGS_H_

typedef struct {
    char *name;
    int32_t len;
} contig_t;

typedef struct {
    contig_t *contigs;
    int32_t n;
} contigs_t;

contigs_t* 
contigs_init();

void 
contigs_add(contigs_t *c, char *name, uint32_t len);

void 
contigs_destroy(contigs_t *c);

#endif
