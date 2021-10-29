//
// Created by dlcgold on 28/10/21.
//

#ifndef RLPBWT_PBWT_RLROW_H
#define RLPBWT_PBWT_RLROW_H


class pbwt_rlrow {
public:
    int p;
    int perm_p;
    int next_perm;

    pbwt_rlrow(int p, int permP, int nextPerm);
    virtual ~pbwt_rlrow();
    int lf_mapping(unsigned int);
};


#endif //RLPBWT_PBWT_RLROW_H
