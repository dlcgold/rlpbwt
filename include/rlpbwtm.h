//
// Created by dlcgold on 17/01/22.
//

#ifndef RLPBWT_RLPBWTM_H
#define RLPBWT_RLPBWTM_H


#include <ostream>

class rlpbwtm {
private:
    unsigned int begin;
    unsigned int end;
    unsigned int nhaplo;
public:
    bool operator==(const rlpbwtm &rhs) const;

    bool operator!=(const rlpbwtm &rhs) const;

    friend std::ostream &operator<<(std::ostream &os, const rlpbwtm &rlpbwtm);

    rlpbwtm(unsigned int begin, unsigned int end, unsigned int nhaplo);

    virtual ~rlpbwtm();

};


#endif //RLPBWT_RLPBWTM_H
