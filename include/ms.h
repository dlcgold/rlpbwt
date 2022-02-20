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

class ms {
public:
    std::vector<unsigned int> row;
    std::vector<unsigned int> len;

    ms();

    ms(std::vector<unsigned int> pos, std::vector<unsigned int> len);

    explicit ms(unsigned int size);

    friend std::ostream &operator<<(std::ostream &os, const ms &ms);
};


#endif //RLPBWT_MS_H
