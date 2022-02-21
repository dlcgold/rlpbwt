//
// Created by dlcgold on 21/02/22.
//

#ifndef RLPBWT_COLUMN_NAIVE_H
#define RLPBWT_COLUMN_NAIVE_H

#include <sdsl/int_vector.hpp>
#include <ostream>
#include <utility>

class column_naive {
public:
    bool zero_first{};
    unsigned int count_0{};
    sdsl::int_vector<> p;
    sdsl::int_vector<> uv;
    sdsl::int_vector<> lcp;

    column_naive(bool zeroFirst, unsigned int count0,
                 sdsl::int_vector<> p, sdsl::int_vector<> uv,
                 sdsl::int_vector<> lcp);

    column_naive();

    virtual ~column_naive();

    friend std::ostream &
    operator<<(std::ostream &os, const column_naive &naive);
    unsigned long long size_in_bytes(bool verbose = false) const;
    double size_in_mega_bytes(bool verbose = false) const;
    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr,
                     const std::string &name = "");

    void load(std::istream &in);
};


#endif //RLPBWT_COLUMN_NAIVE_H
