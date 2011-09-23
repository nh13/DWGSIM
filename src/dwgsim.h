#ifndef DWGSIM_H
#define DWGSIM_H

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

#endif
