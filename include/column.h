//
// Created by dlcgold on 28/10/21.
//

#ifndef RLPBWT_COLUMN_H
#define RLPBWT_COLUMN_H

#include <vector>
#include <ostream>
#include <sdsl/vectors.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/vector.hpp>
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

    /**
     * @brief lcp array saved as a sdsl::int_vector
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

private:
    friend class boost::serialization::access;
};

namespace boost {
    namespace serialization {
        template<class Archive>
        void serialize(Archive &a, column &e,
                       const unsigned version){
            a & e.count_0 & e.zero_first & e.lcp & e.rows;
        }
    }
}

#endif //RLPBWT_COLUMN_H
