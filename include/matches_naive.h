//
// Created by dlcgold on 21/02/22.
//

#ifndef RLPBWT_MATCHES_NAIVE_H
#define RLPBWT_MATCHES_NAIVE_H


#include <vector>
#include <tuple>
#include <ostream>

/**
 * @brief class to represent naive matches
 */
class matches_naive {
public:
    /**
     * @brief vector of tuple that represent a naive match. In order:
     * - number of haplotypes that are matching
     * - length of the match
     * - ending column of the match
     */
    std::vector<std::tuple<unsigned int, unsigned int, unsigned int>> basic_matches;

    /**
     * @brief default constructor
     */
    matches_naive();

    /**
    * @brief function to obtain size in bytes of the mathces
    * @param verbose bool for extra prints
    * @return size in bytes
    */
    unsigned long long size_in_bytes() const;

    /**
     * @brief function to obtain size in megabytes of the matches
     * @param verbose bool for extra prints
     * @return size in megabytes
     */
    double size_in_mega_bytes() const;

    friend std::ostream &
    operator<<(std::ostream &os, const matches_naive &matches);

};


#endif //RLPBWT_MATCHES_NAIVE_H
