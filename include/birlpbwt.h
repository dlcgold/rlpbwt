//
// Created by dlcgold on 20/01/22.
//

#ifndef RLPBWT_BIRLPBWT_H
#define RLPBWT_BIRLPBWT_H

#include <iostream>
#include <vector>
#include <list>
#include "rlpbwt.h"
#include "match.h"

class birlpbwt {
public:
    rlpbwt frlpbwt = rlpbwt();
    rlpbwt brlpbwt = rlpbwt();

    explicit birlpbwt(const char *filename, bool verbose = false);

    /**
     * @brief function to compute matches between the panel and a new query
     *
     * @param query an haplotype string of the same length of the panel
     * @param verbose bool for extra print
     * @return a vector of matches (begin, end, number of matches)
     */
    std::vector<match>
    external_match(const std::string &query, bool verbose = false);

};


#endif //RLPBWT_BIRLPBWT_H
