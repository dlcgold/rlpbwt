//
// Created by dlcgold on 28/10/21.
//

#ifndef RLPBWT_PBWT_RLROW_H
#define RLPBWT_PBWT_RLROW_H


class pbwt_rlrow {
public:
    unsigned int p;
    unsigned int perm_p;
    unsigned int next_perm;

    pbwt_rlrow(unsigned int p, unsigned int permP, unsigned int nextPerm);
    virtual ~pbwt_rlrow();
    int lf_mapping(unsigned int) const;
};


#endif //RLPBWT_PBWT_RLROW_H
