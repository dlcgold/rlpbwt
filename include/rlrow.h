//
// Created by dlcgold on 17/01/22.
//

#ifndef RLPBWT_RLROW_H
#define RLPBWT_RLROW_H


#include <ostream>

class rlrow {
public:
    rlrow(unsigned int p, unsigned int uv);

    unsigned int p;
    unsigned int uv;
public:
    friend std::ostream &operator<<(std::ostream &os, const rlrow &rlrow);

public:
    virtual ~rlrow();
};


#endif //RLPBWT_RLROW_H
