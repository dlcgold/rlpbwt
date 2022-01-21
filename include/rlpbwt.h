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
#include "column.h"
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
    std::vector<column> cols;


    /**
     * @brief heigth of the original panel
     */
    unsigned int heigth{};

    /**
     * @brief width of the original panel
     */
    unsigned int width{};

    /**
     * @brief constructor for run-length encoded PBWT matrix
     * @param filename path of the file with the panel, every row of the file is
     * one column of the panel
     * @param verbose bool for extra print
     */
    explicit rlpbwt(const char *filename, bool verbose = false);

    /**
     * @brief default constructor
     */
    rlpbwt();

    /**
     * @brief default destructor
     */
    virtual ~rlpbwt();

    /**
     * @brief function to compute matches between the panel and a new query
     *
     * @param query an haplotype string of the same length of the panel
     * @param verbose bool for extra print
     * @return a vector of matches (begin, end, number of matches)
     */
    std::vector<match>
    external_match(const std::string &query, bool verbose = false);

    /**
     * @brief function to get the run in previous column which come from the
     * current run
     *
     * @param col_index index of the column
     * @param index virtual index of the row of the original panel
     * @param verbose bool for extra print
     * @return
     */
    unsigned int
    prev_run(unsigned int col_index, unsigned int index, bool verbose) const;

    /**
     * @brief function to obtain the struct for the run-length encoded PBWT
     * column, except for next_perm values
     * @param column current column
     * @param pref current prefix array
     * @param div current divergence array
     * @return the struct for the run-length encoded PBWT column
     */
    static column
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
     * @brief function to compute end matches between the panel and a new query
     *
     * @param query an haplotype string of the same length of the panel
     * @param verbose bool for extra print
     * @return a vector of matches (begin, end, number of matches)
     */
    std::vector<match>
    end_external_match(const std::string &query, bool forward = true,
                       bool verbose = false);

private:

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
    lf(unsigned int col_index, unsigned int row_index, char symbol,
       unsigned int offset, bool verbose = false) const;


    /**
     * @brief function to map an index to the correct run in a column
     * @param index index to map
     * @param col_index column index
     * @return run index
     */
    unsigned int index_to_run(unsigned int index, unsigned int col_index) const;


    /**
     * @brief trick to extract u and v value from a run in rlpbwt column
     * @param col_index index of the column
     * @param row_index index of the run
     * @return a std::pair with u as first and v as second
     */
    std::pair<unsigned int, unsigned int>
    uvtrick(unsigned int col_index, unsigned int row_index) const;

};


#endif //RLPBWT_RLPBWT_H

