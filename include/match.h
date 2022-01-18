//
// Created by dlcgold on 17/01/22.
//

#ifndef RLPBWT_MATCH_H
#define RLPBWT_MATCH_H


#include <ostream>

class match {
private:
    unsigned int begin;
    unsigned int end;
    unsigned int nhaplo;
public:
    bool operator==(const match &rhs) const;

    bool operator!=(const match &rhs) const;

    friend std::ostream &operator<<(std::ostream &os, const match &rlpbwtm);

    match(unsigned int begin, unsigned int end, unsigned int nhaplo);

    virtual ~match();

};


#endif //RLPBWT_MATCH_H
