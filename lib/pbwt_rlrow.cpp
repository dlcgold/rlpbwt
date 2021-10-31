//
// Created by dlcgold on 28/10/21.
//

#include "../include/pbwt_rlrow.h"

pbwt_rlrow::pbwt_rlrow(unsigned int p, unsigned int permP,
                       unsigned int nextPerm, unsigned int threshold) : p(p),
                                                                        perm_p(permP),
                                                                        next_perm(
                                                                                nextPerm),
                                                                        threshold(
                                                                                threshold) {}

pbwt_rlrow::~pbwt_rlrow() = default;

unsigned int pbwt_rlrow::lf_mapping(unsigned int i) const {
    return perm_p + i - p;
}

std::ostream &operator<<(std::ostream &os, const pbwt_rlrow &rlrow) {
    os << rlrow.p << "\t" << rlrow.perm_p << "\t"
       << rlrow.next_perm << "\t" << rlrow.threshold;
    return os;
}

bool pbwt_rlrow::operator==(const pbwt_rlrow &rhs) const {
    return p == rhs.p &&
           perm_p == rhs.perm_p &&
           next_perm == rhs.next_perm &&
           threshold == rhs.threshold;
}

bool pbwt_rlrow::operator!=(const pbwt_rlrow &rhs) const {
    return !(rhs == *this);
}
