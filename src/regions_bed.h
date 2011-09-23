#ifndef REGIONS_BED_H
#define REGIONS_BED_H

typedef struct {
    uint32_t *contig;
    uint32_t *start;
    uint32_t *end;
    uint32_t n;
    uint32_t mem;
} regions_bed_txt;

regions_bed_txt *
regions_bed_init(FILE *fp, contigs_t *c);

void 
regions_bed_destroy(regions_bed_txt *r);

int32_t 
regions_bed_query(regions_bed_txt *r, uint32_t contig, uint32_t start, uint32_t end); 

#endif
