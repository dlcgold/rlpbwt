//
// Created by dlcgold on 20/02/22.
//

#ifndef RLPBWT_MS_MATCHES_H
#define RLPBWT_MS_MATCHES_H

#include <utility>
#include <iostream>
#include <vector>
#include <tuple>
#include <ostream>

class ms_matches {
public:
    std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> basic_matches;
    std::vector<std::vector<unsigned int>> haplos;

    friend std::ostream &operator<<(std::ostream &os, const ms_matches &ms);
    ms_matches();
};


#endif //RLPBWT_MS_MATCHES_H
