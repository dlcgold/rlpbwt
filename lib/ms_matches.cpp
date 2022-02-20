//
// Created by dlcgold on 20/02/22.
//

#include "../include/ms_matches.h"

ms_matches::ms_matches() = default;


std::ostream &operator<<(std::ostream &os, const ms_matches &ms) {
    os << "\nmatches (<ending>, <length>, [haplotypes]):\n";
    bool haplo = ms.haplos.empty();
    for (unsigned int i = 0; i < ms.basic_matches.size(); i++) {
        os << std::get<2>(ms.basic_matches[i]) << ", " << std::get<1>(ms.basic_matches[i]);
        if (!haplo) {
            os << ", [";
            for (auto h: ms.haplos[i]) {
                os << h << " ";
            }
            os << "]";
        }
        os << "\n";
    }
    return os;
}