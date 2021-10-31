//
// Created by dlcgold on 28/10/21.
//

#ifndef RLPBWT_PBWT_COLUMN_H
#define RLPBWT_PBWT_COLUMN_H

#include <vector>
#include <ostream>
#include "pbwt_rlrow.h"

/**
 * @brief class to rappresent every column in run-length encoded
 * PBWT matrix
 */
class pbwt_column {
public:
    /**
     * @brief bool to check first value of the column in PBWT matrix
     * (assuming biallelic)
     */
    bool zero_first;

    /**
     * @brief vector with the quadruple for every run in the column in PBWT
     * matrix (assuming biallelic)
     */
    std::vector<pbwt_rlrow> rows;

    /**
     * @brief default constructor
     */
    pbwt_column();

    /**
     * @brief constructor of a run-length encoded PBWT column
     *
     * @param zeroFirst bool to check first value of the column
     * @param rows vector with every quadruple for every run
     */
    pbwt_column(bool zeroFirst, std::vector<pbwt_rlrow> rows);

    /**
     * @brief default destructor
     */
    virtual ~pbwt_column();

    /**
     * @brief ostream overload to print the struct for a column
     */
    friend std::ostream &
    operator<<(std::ostream &os, const pbwt_column &column);

};


#endif //RLPBWT_PBWT_COLUMN_H
