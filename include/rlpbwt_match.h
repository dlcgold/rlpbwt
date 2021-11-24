//
// Created by dlcgold on 10/11/21.
//

#ifndef RLPBWT_RLPBWT_MATCH_H
#define RLPBWT_RLPBWT_MATCH_H


#include <vector>

class rlpbwt_match {
public:
    rlpbwt_match(unsigned int begin, unsigned int length,
          std::vector<unsigned int> rows);

    virtual ~rlpbwt_match();

    unsigned int begin;
    unsigned int length;
    std::vector<unsigned int> rows;
};


#endif //RLPBWT_RLPBWT_MATCH_H
