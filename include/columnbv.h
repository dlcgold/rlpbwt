//
// Created by dlcgold on 28/01/22.
//

#ifndef RLPBWT_COLUMNBV_H
#define RLPBWT_COLUMNBV_H

#include <vector>
#include <sdsl/int_vector.hpp>
#include <sdsl/bit_vectors.hpp>

class columnbv {
public:
    /**
    * @brief bool to check first value of the column in PBWT matrix
    * (assuming biallelic)
    */
    bool zero_first{};

    /**
     * @brief total number of zeros in the column
     */
    unsigned int count_0{};

    sdsl::sd_vector<> runs;
    sdsl::sd_vector<>::rank_1_type rank_runs;
    sdsl::sd_vector<>::select_1_type select_runs;

    sdsl::sd_vector<> u;
    sdsl::sd_vector<>::rank_1_type rank_u;
    sdsl::sd_vector<>::select_1_type select_u;

    sdsl::sd_vector<> v;
    sdsl::sd_vector<>::rank_1_type rank_v;
    sdsl::sd_vector<>::select_1_type select_v;

    sdsl::int_vector<> lcp;

/*
    columnbv(bool zeroFirst, unsigned int count0,
             sdsl::sd_vector<> runs,
             sdsl::sd_vector<>::rank_1_type rankRuns,
             sdsl::sd_vector<>::select_1_type selectRuns,
             sdsl::sd_vector<> u,
             sdsl::sd_vector<>::rank_1_type rankU,
             sdsl::sd_vector<>::select_1_type selectU,
             sdsl::sd_vector<> v,
             sdsl::sd_vector<>::rank_1_type rankV,
             sdsl::sd_vector<>::select_1_type selectV,
             std::vector<unsigned int> lcp);
    */
    columnbv(bool zeroFirst, unsigned int count0,
             const sdsl::bit_vector& runs,
             const sdsl::bit_vector& u,
             const sdsl::bit_vector& v,
             sdsl::int_vector<> lcp);

    /**
     * @brief default constructor
     */
    columnbv();

    virtual ~columnbv();
};


#endif //RLPBWT_COLUMNBV_H
