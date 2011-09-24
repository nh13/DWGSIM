#ifndef DWGSIM_H
#define DWGSIM_H

#include "dwgsim_opt.h"

#define __gf_add(_x, _y) ((_x >= 4 || _y >= 4) ? 4 : (_x ^ _y))
#define __IS_TRUE(_val) ((_val == 1) ? "True" : "False")

extern uint8_t nst_nt4_table[256];

enum data_type_t {
    ILLUMINA=0,
    SOLID=1,
    IONTORRENT=2
};

char iupac_and_base_to_mut(char iupac, char base);

int32_t get_muttype(char *str);

int32_t
generate_errors_flows(dwgsim_opt_t *opt, uint8_t **seq, uint8_t **mask, int32_t *mem, int32_t len, uint8_t strand, double e, int32_t *_n_err);

#endif
