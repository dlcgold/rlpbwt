//
// Created by dlcgold on 17/01/22.
//

#include "../include/rlpbwtm.h"

rlpbwtm::rlpbwtm(unsigned int begin, unsigned int end, unsigned int nhaplo)
        : begin(begin), end(end), nhaplo(nhaplo) {}

rlpbwtm::~rlpbwtm() {

}

std::ostream &operator<<(std::ostream &os, const rlpbwtm &rlpbwtm) {
    os << "match in [" << rlpbwtm.begin << ", " << rlpbwtm.end << "] with "
       << rlpbwtm.nhaplo << " haplotypes";
    return os;
}
