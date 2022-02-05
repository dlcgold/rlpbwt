//
// Created by dlcgold on 20/01/22.
//

#ifndef RLPBWT_BIRLPBWT_H
#define RLPBWT_BIRLPBWT_H

#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include "rlpbwt.h"
#include "match.h"

/**
 * @brief class to rappresent a double RLPBWT in order to identify maximal
 * matches without saving lcp arrays
 */
class birlpbwt {
public:

    /**
     * @brief regular RLPBWT (forward) builded from left to right
     */
    rlpbwt frlpbwt = rlpbwt();

    /**
     * @brief reverdse RLPBWT (backward) builded from right to left
     * (this will be queried with rhe reverse of the query)
     */
    rlpbwt brlpbwt = rlpbwt();

    /**
     * Constructor for BIRLPWBT
     *
     * @param filename file for the input
     * @param vcf bool to identify that the input file is a vcf file
     * @param verbose bool for extra print
     */
    explicit birlpbwt(const char *filename, bool vcf, bool verbose = false);

    /**
     * @brief function to compute matches between the panel and a new query
     *
     * @param query an haplotype string of the same length of the panel
     * @param verbose bool for extra print
     * @return a vector of matches (begin, end, number of matches)
     */
    std::vector<match>
    external_match(const std::string &query, unsigned int min_len = 1,
                   bool verbose = false);

    /**
     * @brief function to compute matches between the panel and a new query
     * from a vcf file
     * @param query an haplotype string of the same length of the panel
     * @param min_len minimum length of a match
     * @param verbose bool for extra print
     * @return a vector of matches (begin, end, number of matches)
     */
    void external_match_vcf(const char *filename, unsigned int min_len = 1,
                            bool verbose = false);


    /**
     * utility to print both forward and backward RLPBWT
     */
    void print();
};


#endif //RLPBWT_BIRLPBWT_H
