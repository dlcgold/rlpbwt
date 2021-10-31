//
// Created by dlcgold on 28/10/21.
//

#ifndef RLPBWT_RLPBWT_RLROW_H
#define RLPBWT_RLPBWT_RLROW_H

#include <ostream>

/**
 * @brief class to rappresent every row of a column in run-length encoded
 * PBWT matrix
 */
class rlpbwt_rlrow {
public:
    /**
     * @brief position of i-th run in the j-th column of the PBWT column
     */
    unsigned int p;

    /**
     * @brief permutation pi_j(p)
     */
    unsigned int perm_p;

    /**
     * @brief index of the run containing pi_j(p) in j+1th column of the
     * PBWT column
     */
    unsigned int next_perm;

    /**
     * @brief threshold of the run
     */
    unsigned int threshold;

    /**
     * @brief constructor of the quadruple for every row
     *
     * @param p position of i-th run in the j-th column of the PBWT column
     * @param perm_p permutation pi_j(p)
     * @param next_perm index of the run containing pi_j(p) in j+1th column of
     * the PBWT column
     * @param threshold threshold of the run
     */
    rlpbwt_rlrow(unsigned int p, unsigned int perm_p, unsigned int next_perm,
               unsigned int threshold);

    /**
     * @brief default destructor
     */
    virtual ~rlpbwt_rlrow();

    /**
     * @brief function to find the row for the run containing the bit
     * obtained with previous permutations
     *
     * @param i result of the previous permutation: i, p_1(i), p_2(p_1(i))...
     * @return the row for the run containing the bit
     * obtained with previous permutations
     */
    unsigned int lf_mapping(unsigned int i) const;

    /**
     * @brief ostream overload to print the struct for a row with the 4 values
     */
    friend std::ostream &operator<<(std::ostream &os, const rlpbwt_rlrow &rlrow);

    /**
     * @brief Equality operator for a row of a run-length encoded PBWT column
     */
    bool operator==(const rlpbwt_rlrow &rhs) const;

    /**
     * @brief Inequality operator for a row of a run-length encoded PBWT column
     */
    bool operator!=(const rlpbwt_rlrow &rhs) const;
};


#endif //RLPBWT_RLPBWT_RLROW_H
