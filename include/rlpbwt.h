//
// Created by dlcgold on 28/10/21.
//

#ifndef RLPBWT_RLPBWT_H
#define RLPBWT_RLPBWT_H

#include <vector>
#include <string>
#include <iostream>
#include <climits>
#include "rlpbwt_column.h"

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
     */
    explicit rlpbwt(const char *filename);

    /**
     * @brief function to extract a row of the original panel from the
     * run-length encoded PBWT matrix
     * @param i index of the query row
     * @return the queried row in a std::string
     */
    std::string search_row(unsigned int i);

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
    static void
    update(std::string &column, std::vector<unsigned int> &pref,
           std::vector<unsigned int> &div);
};


#endif //RLPBWT_RLPBWT_H
