//
// Created by dlcgold on 29/10/21.
//

#ifndef RLPBWT_UTILS_H
#define RLPBWT_UTILS_H

#include <string>
#include <vector>
#include <climits>
#include "rlpbwt_column.h"

/**
 * @brief function to get the char (in bialleic case with 0 and 1) at certain
 * run index
 * @param zero_first bool to check first value of the column
 * @param index_run run index in the run-length encoded PBWT column
 * @return the char at the queried run
 */
char get_next_char(bool zero_first, unsigned int index_run);

template<typename T> double vectorsizeof(const typename std::vector<T> &vec) {
    return (sizeof(T) * vec.size()) * 0.000001;
}
#endif //RLPBWT_UTILS_H
