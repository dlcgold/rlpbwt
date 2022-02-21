//
// Created by dlcgold on 21/02/22.
//

#ifndef RLPBWT_MATCHES_NAIVE_H
#define RLPBWT_MATCHES_NAIVE_H


#include <vector>
#include <tuple>
#include <ostream>

class matches_naive {
public:
    std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> basic_matches;
    unsigned long long size_in_bytes() const;
    double size_in_mega_bytes() const;
    friend std::ostream &operator<<(std::ostream &os, const matches_naive &matches);
    matches_naive();
};


#endif //RLPBWT_MATCHES_NAIVE_H
