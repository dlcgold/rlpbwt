//
// Created by dlcgold on 28/10/21.
//

#ifndef RLPBWT_RLPBWT_H
#define RLPBWT_RLPBWT_H

#include <vector>
#include "pbwt_column.h"

/**
 * @brief class to rappresent the run-length encoded PBWT matrix
 */
class rlpbwt {
public:
    /**
     * @brief vector with all the structs for every column in run-length encoded
     * PBWT matrix
     */
    std::vector<pbwt_column> cols;

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
    explicit rlpbwt(char *filename);

    /**
     * @brief fucntion to extract a row of the original panel from the
     * run-length encoded PBWT matrix
     * @param i index of the query row
     * @return the queried row in a std::string
     */
    std::string search_row(unsigned int i);

    /**
     * @brief default destructor
     */
    virtual ~rlpbwt();
};


#endif //RLPBWT_RLPBWT_H
