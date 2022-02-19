//
// Created by dlcgold on 19/02/22.
//

#ifndef RLPBWT_MS_H
#define RLPBWT_MS_H


#include <vector>

class ms {
public:
    std::vector<unsigned int> pos;
    std::vector<unsigned int> len;
    std::vector<std::pair<unsigned int, unsigned int>> matches;
    std::vector<std::vector<unsigned int>> haplos;

    ms();

    ms(std::vector<unsigned int> pos,
       std::vector<unsigned int> len,
       std::vector<std::pair<unsigned int, unsigned int>> matches);
};


#endif //RLPBWT_MS_H
