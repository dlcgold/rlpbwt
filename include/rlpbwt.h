//
// Created by dlcgold on 28/10/21.
//

#ifndef RLPBWT_RLPBWT_H
#define RLPBWT_RLPBWT_H

#include <vector>
#include <string>
#include <iostream>
#include <climits>
#include <cassert>
#include "rlpbwt_column.h"
#include "rlpbwt_match.h"

/**
 * @brief class to rappresent the run-length encoded PBWT matrix
 */
class rlpbwt {
public:
    /**
     * @brief vector with all the structs for every column in run-length encoded
     * PBWT matrix
     */
    std::vector<rlpbwt_column> cols;

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
     * @brief function to extract a row of the original panel from the
     * run-length encoded PBWT matrix
     * @param i index of the query row
     * @param verbose bool for extra print
     * @return the queried row in a std::string
     */
    std::string search_row(unsigned int i, bool verbose = false);

    std::vector<rlpbwt_match> external_match(const std::string &query) const;


    /**
     * @brief default destructor
     */
    virtual ~rlpbwt();

private:

    /**
     * @brief fucntion to obtain the next_perm values for the ith run-length
     * encoded PBWT column from the i+1th one
     * @param prev ith run-length encoded PBWT column, the prevoius one
     * @param curr i+1th run-length encoded PBWT column, the current one
     */
    static void build_next_perm(rlpbwt_column &prev, rlpbwt_column &curr);

    /**
     * @brief function to obtain the struct for the run-length encoded PBWT
     * column, except for next_perm values
     * @param column current column
     * @param pref current prefix array
     * @param div current divergence array
     * @return the struct for the run-length encoded PBWT column
     */
    static rlpbwt_column
    build_column(std::string &column, std::vector<unsigned int> &pref,
                 std::vector<unsigned int> &div);


    /**
     * @brief function to get pref and div arrays at column i+1 from the values
     * in column i
     * @param column current column
     * @param pref prefix array for previous column
     * @param div divergence array for previous column
     */
    __attribute__((unused)) static void
    update_old(std::string &column, std::vector<unsigned int> &pref,
               std::vector<unsigned int> &div, unsigned int k);

    unsigned int
    next_run(unsigned int col_index, unsigned int start, unsigned int end,
             bool verbose) const;

    unsigned int
    occ(unsigned int col_index, unsigned int row_index, char symbol) const;

    unsigned int index_to_run(unsigned int index) const;
    unsigned int get_run_length(unsigned int col_index, unsigned int run_index) const;
    static void
    update(std::string &column, std::vector<unsigned int> &pref,
           std::vector<unsigned int> &div);

    std::vector<unsigned int>
    update_external(unsigned int index, unsigned int e, unsigned int f,
                    unsigned int g, const std::string &query);

};


#endif //RLPBWT_RLPBWT_H
