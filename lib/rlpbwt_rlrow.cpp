//
// Created by dlcgold on 28/10/21.
//

#include "../include/rlpbwt_rlrow.h"

rlpbwt_rlrow::rlpbwt_rlrow(unsigned int p, unsigned int permP,
                       unsigned int nextPerm,
                       unsigned int threshold) : p(p),
                                                 perm_p(permP),
                                                 next_perm(nextPerm),
                                                 threshold(threshold) {}

rlpbwt_rlrow::~rlpbwt_rlrow() = default;

unsigned int rlpbwt_rlrow::lf_mapping(unsigned int i) const {
    return perm_p + i - p;
}

std::ostream &operator<<(std::ostream &os, const rlpbwt_rlrow &rlrow) {
//    os << rlrow.p << "\t" << rlrow.perm_p << "\t"
//       << rlrow.next_perm << "\t" << rlrow.threshold;
    os << rlrow.p << "\t" << rlrow.perm_p << "\t"
       << rlrow.next_perm;
    return os;
}

bool rlpbwt_rlrow::operator==(const rlpbwt_rlrow &rhs) const {
    return p == rhs.p &&
           perm_p == rhs.perm_p &&
           next_perm == rhs.next_perm &&
           threshold == rhs.threshold;
}

bool rlpbwt_rlrow::operator!=(const rlpbwt_rlrow &rhs) const {
    return !(rhs == *this);
}
