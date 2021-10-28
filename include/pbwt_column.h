//
// Created by dlcgold on 28/10/21.
//

#ifndef RLPBWT_PBWT_COLUMN_H
#define RLPBWT_PBWT_COLUMN_H

#include <vector>
#include "pbwt_rlrow.h"

class pbwt_column {
public:
    bool zero_first;
    std::vector<pbwt_rlrow> rows;

    pbwt_column(bool zeroFirst, const std::vector<pbwt_rlrow> &rows);
    virtual ~pbwt_column();

};


#endif //RLPBWT_PBWT_COLUMN_H
