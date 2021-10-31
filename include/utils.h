//
// Created by dlcgold on 29/10/21.
//

#ifndef RLPBWT_UTILS_H
#define RLPBWT_UTILS_H

#include <string>
#include <vector>
#include <climits>
#include "pbwt_column.h"

/**
 * @brief function to get pref and div arrays at column i+1 from the values in
 * column i
 * @param column current column
 * @param pref prefix array for previous column
 * @param div divergence array for previous column
 */
void
update(std::string &column, std::vector<unsigned int> &pref,
       std::vector<unsigned int> &div);

/**
 * @brief function to obtain the struct for the run-length encoded PBWT column,
 * except for next_perm values
 * @param column current column
 * @param pref current prefix array
 * @param div current divergence array
 * @return the struct for the run-length encoded PBWT column
 */
pbwt_column
build_column(std::string &column, std::vector<unsigned int> &pref,
             std::vector<unsigned int> &div);

/**
 * @brief fucntion to obtain the next_perm values for the ith run-length encoded
 * PBWT column from the i+1th one
 * @param prev ith run-length encoded PBWT column, the prevoius one
 * @param curr i+1th run-length encoded PBWT column, the current one
 */
void build_next_perm(pbwt_column &prev, pbwt_column &curr);

/**
 * @brief function to get the char (in bialleic case with 0 and 1) at certain
 * run index
 * @param zero_first bool to check first value of the column
 * @param run run index in the run-length encoded PBWT column
 * @return the char at the queried run
 */
char get_next_char(bool zero_first, unsigned int run);

#endif //RLPBWT_UTILS_H
