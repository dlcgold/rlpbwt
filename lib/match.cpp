//
// Created by dlcgold on 17/01/22.
//

#include "../include/match.h"

match::match(unsigned int begin, unsigned int end, unsigned int nhaplo)
        : begin(begin), end(end), nhaplo(nhaplo) {}

match::~match() = default;

std::ostream &operator<<(std::ostream &os, const match &rlpbwtm) {
    os << "match in [" << rlpbwtm.begin << ", " << rlpbwtm.end << "] with "
       << rlpbwtm.nhaplo << " haplotypes";
    return os;
}

bool match::operator==(const match &rhs) const {
    return begin == rhs.begin &&
           end == rhs.end &&
           nhaplo == rhs.nhaplo;
}

bool match::operator!=(const match &rhs) const {
    return !(rhs == *this);
}
