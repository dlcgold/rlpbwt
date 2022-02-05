//
// Created by dlcgold on 28/10/21.
//

#ifndef RLPBWT_COLUMN_H
#define RLPBWT_COLUMN_H

#include <vector>
#include <ostream>
//#include <sdsl/vectors.hpp>
#include <sdsl/structure_tree.hpp>
#include <sdsl/util.hpp>
#include "utils.h"
#include "rlrow.h"

/**
 * @brief class to rappresent every column in run-length encoded
 * PBWT matrix
 */
class column {
public:
    /**
     * @brief bool to check first value of the column in PBWT matrix
     * (assuming biallelic)
     */
    bool zero_first;

    /**
     * @brief total number of zeros in the column
     */
    unsigned int count_0;

    /**
     * @brief vector with the quadruple for every run in the column in PBWT
     * matrix (assuming biallelic)
     */
    std::vector<rlrow> rows;

    // TODO think about lcp saved as a sdsl::int_vector<> with compression
    /**
     * @brief lcp array
     */
    std::vector<unsigned int> lcp;

    /**
     * @brief default constructor
     */
    column();

    /**
     * @brief constructor of a run-length encoded PBWT column
     *
     * @param zeroFirst bool to check first value of the column
     * @param rows vector with every quadruple for every run
     */
    column(bool zeroFirst, std::vector<rlrow> rows,
           unsigned int count0);

    /**
     * @brief default destructor
     */
    virtual ~column();

    /**
     * @brief ostream overload to print the struct for a column
     */
    friend std::ostream &
    operator<<(std::ostream &os, const column &column);
    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr,
                     const std::string &name = "");

    void load(std::istream &in);
};


#endif //RLPBWT_COLUMN_H
