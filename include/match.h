//
// Created by dlcgold on 17/01/22.
//

#ifndef RLPBWT_MATCH_H
#define RLPBWT_MATCH_H


#include <ostream>

/**
 * @brief class to rapresent a maximal match
 */
class match {
public:
    /**
     * @brief index on columns for the begin of a match
     */
    unsigned int begin;

    /**
     * @brief index on columns for the end of a match
     */
    unsigned int end;

    /**
     * @brief number of haplotypes in the match
     */
    unsigned int nhaplo;

    bool operator==(const match &rhs) const;

    bool operator!=(const match &rhs) const;

    friend std::ostream &operator<<(std::ostream &os, const match &rlpbwtm);

    match(unsigned int begin, unsigned int end, unsigned int nhaplo);

    virtual ~match();

};


#endif //RLPBWT_MATCH_H
