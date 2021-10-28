//
// Created by dlcgold on 28/10/21.
//

#ifndef RLPBWT_RLPBWT_H
#define RLPBWT_RLPBWT_H

#include <vector>
#include "pbwt_column.h"

class rlpbwt {
public:
    std::vector<pbwt_column> cols;

    rlpbwt(const std::vector<pbwt_column> &cols);
    virtual ~rlpbwt();
};


#endif //RLPBWT_RLPBWT_H
