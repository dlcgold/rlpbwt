//
// Created by dlcgold on 17/01/22.
//

#include "../include/rlrow.h"

rlrow::rlrow(unsigned int p, unsigned int uv) : p(p), uv(uv) {}

std::ostream &operator<<(std::ostream &os, const rlrow &rlrow) {
    os << rlrow.p << "\t" << rlrow.uv;
    return os;
}

rlrow::~rlrow() = default;
