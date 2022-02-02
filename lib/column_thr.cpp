//
// Created by dlcgold on 02/02/22.
//

#include "../include/column_thr.h"

#include <utility>

column_thr::column_thr(bool zeroFirst, unsigned int count0,
                       const sdsl::bit_vector &runs, const sdsl::bit_vector &u,
                       const sdsl::bit_vector &v, const sdsl::bit_vector &thr,
                       std::vector<std::pair<unsigned int, unsigned int>> rows)
        : zero_first(zeroFirst), count_0(count0), rows(std::move(rows)) {
    this->runs = sdsl::sd_vector<>(runs);
    this->u = sdsl::sd_vector<>(u);
    this->v = sdsl::sd_vector<>(v);
    this->thr = sdsl::sd_vector<>(thr);
}

std::ostream &operator<<(std::ostream &os, const column_thr &thr) {
    auto yesno = "yes";
    if (!thr.zero_first) {
        yesno = "no";
    }
    os << "zero_first: " << yesno << " c: " << thr.count_0
       << "\nruns: " << thr.runs << "\nu: " << thr.u << "\nv: " << thr.v
       << "\nthr: " << thr.thr << "\n";
    for (auto e: thr.rows) {
        std::cout << "(" << e.first << ", " << e.second << ") ";
    }
    std::cout << "\n";
    return os;
}

column_thr::~column_thr() = default;

column_thr::column_thr() = default;

