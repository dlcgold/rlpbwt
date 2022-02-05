//
// Created by dlcgold on 17/01/22.
//

#ifndef RLPBWT_RLROW_H
#define RLPBWT_RLROW_H

#include <ostream>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/vector.hpp>

/**
 * @brief class to rappresent a single row in a column of RLPBWT
 */
class rlrow {
public:
    /**
     * @brief index of the head of the run
     */
    unsigned int p;

    /**
    * @brief value of u or v (in biallelic case we can extract both based on
    * the status of the column (first symbol and even/odd index for the row)
    */
    unsigned int uv;

    /**
     * Constrcutor for a single row
     *
     * @param p index of head run
     * @param uv value for u or v
     */
    rlrow(unsigned int p, unsigned int uv);

    /**
     * @brief ostream overload to print the struct for a row
     */
    friend std::ostream &operator<<(std::ostream &os, const rlrow &rlrow);

    /**
    * @brief default destructor
    */
    virtual ~rlrow();
};


#endif //RLPBWT_RLROW_H
