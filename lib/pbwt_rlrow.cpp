//
// Created by dlcgold on 28/10/21.
//

#include "../include/pbwt_rlrow.h"

pbwt_rlrow::pbwt_rlrow(unsigned int p, unsigned int permP, unsigned int nextPerm) : p(p),
                                                         perm_p(permP),
                                                         next_perm(nextPerm) {}

pbwt_rlrow::~pbwt_rlrow() = default;

int pbwt_rlrow::lf_mapping(unsigned int i) const {
    return perm_p + i - p;
}
