//
// Created by dlcgold on 02/02/22.
//

#ifndef RLPBWT_RLPBWT_THR_H
#define RLPBWT_RLPBWT_THR_H


#include <vector>
#include <string>
#include <iostream>
#include <climits>
#include <cassert>
#include <utility>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <sdsl/bit_vectors.hpp>
#include "column_thr.h"
#include "utils.h"
#include "exceptions.h"

class rlpbwt_thr {
public:
    std::vector<column_thr> cols;
    std::vector<std::string> panel;

    /**
     * @brief constructor for run-length encoded PBWT matrix
     * @param filename path of the file with the panel, every row of the file is
     * one column of the panel
     * @param verbose bool for extra print
     */
    explicit rlpbwt_thr(const char *filename, bool vcf, bool verbose = false);

    /**
     * @brief function to obtain the struct for the run-length encoded PBWT
     * column, except for next_perm values
     * @param column current column
     * @param pref current prefix array
     * @param div current divergence array
     * @return the struct for the run-length encoded PBWT column
     */
    static column_thr
    build_column(std::string &column, std::vector<unsigned int> &pref,
                 sdsl::int_vector<> &div);

    /**
     * @brief utility to compute prefix and divergence array
     * @param column the current column
     * @param pref the previous prefix array
     * @param div the previous divergence array
     */
    static void
    update(std::string &column, std::vector<unsigned int> &pref,
           sdsl::int_vector<> &div);

    /**
    * @brief function to compute the lf mapping, w(i, s) function in Durbin
    *
    * @param col_index index of the column
    * @param row_index index of the run
    * @param symbol symbol s
    * @param offset offset to correctly calculate the mapping
    * @param verbose bool for extra print
    * @return the index computed with the lf-mapping
    */
    unsigned int
    lf(unsigned int col_index, unsigned int index, char symbol,
       bool verbose = false) const;

    /**
    * @brief trick to extract u and v value from a run in rlpbwt column
    * @param col_index index of the column
    * @param row_index index of the run
    * @return a std::pair with u as first and v as second
    */
    std::pair<unsigned int, unsigned int>
    uvtrick(unsigned int col_index, unsigned int index) const;
    void match_thr(const std::string &query, bool verbose = false);
};


#endif //RLPBWT_RLPBWT_THR_H
