//
// Created by dlcgold on 28/10/21.
//

#include "../include/pbwt_rlrow.h"

pbwt_rlrow::pbwt_rlrow(int p, int permP, int nextPerm) : p(p),
                                                         perm_p(permP),
                                                         next_perm(nextPerm) {}

pbwt_rlrow::~pbwt_rlrow() {

}

int pbwt_rlrow::lf_mapping(unsigned int i) {
    return perm_p + i - p;
}
