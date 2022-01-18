//
// Created by dlcgold on 04/01/22.
//

#ifndef RLPBWT_RLPBWT_H
#define RLPBWT_RLPBWT_H


#include <vector>
#include <string>
#include <iostream>
#include <climits>
#include <cassert>
#include <utility>
#include <fstream>
#include <algorithm>
#include "col.h"
#include "utils.h"
#include "exceptions.h"
#include "match.h"

/**
 * @brief class to rappresent the run-length encoded PBWT matrix
 */
class rlpbwt {
public:
    /**
     * @brief vector with all the structs for every column in run-length encoded
     * PBWT matrix
     */
    std::vector<col> cols;

    /**
     * @brief heigth of the original panel
     */
    unsigned int heigth;

    /**
     * @brief width of the original panel
     */
    unsigned int width;

    /**
     * @brief constructor for run-length encoded PBWT matrix
     * @param filename path of the file with the panel, every row of the file is
     * one column of the panel
     * @param verbose bool for extra print
     */
    explicit rlpbwt(const char *filename, bool verbose = false);

    /**
     * @brief default destructor
     */
    virtual ~rlpbwt();

    std::vector<match> ematch(const std::string &query, bool verbose = false);

    unsigned int
    prev_run(unsigned int col_index, unsigned int index, bool verbose) const;

private:


    /**
     * @brief function to obtain the struct for the run-length encoded PBWT
     * column, except for next_perm values
     * @param column current column
     * @param pref current prefix array
     * @param div current divergence array
     * @return the struct for the run-length encoded PBWT column
     */
    static col
    build_column(std::string &column, std::vector<unsigned int> &pref,
                 sdsl::int_vector<> &div);


    unsigned int
    occ(unsigned int col_index, unsigned int row_index, char symbol,
        unsigned int offset, bool verbose = false) const;


    unsigned int index_to_run(unsigned int index, unsigned int col_index) const;

    static void
    update(std::string &column, std::vector<unsigned int> &pref,
           sdsl::int_vector<> &div);

    std::pair<unsigned int, unsigned int>
    uvtrick(unsigned int col_index, unsigned int row_index) const;

};


#endif //RLPBWT_RLPBWT_H

