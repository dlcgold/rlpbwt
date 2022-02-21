//
// Created by dlcgold on 19/02/22.
//

#ifndef RLPBWT_MS_H
#define RLPBWT_MS_H


#include <utility>
#include <iostream>
#include <vector>
#include <tuple>
#include <ostream>
#include <sdsl/int_vector.hpp>

class ms {
public:
    std::vector<unsigned int> row;
    std::vector<unsigned int> len;
    ms();

    ms(std::vector<unsigned int> pos, std::vector<unsigned int> len);

    explicit ms(unsigned int size);
    unsigned long long size_in_bytes() const;
    double size_in_mega_bytes() const;

    friend std::ostream &operator<<(std::ostream &os, const ms &ms);
};


#endif //RLPBWT_MS_H
