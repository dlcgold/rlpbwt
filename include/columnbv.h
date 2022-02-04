//
// Created by dlcgold on 28/01/22.
//

#ifndef RLPBWT_COLUMNBV_H
#define RLPBWT_COLUMNBV_H

#include <vector>
#include <sdsl/int_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include <ostream>

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

    friend std::ostream &operator<<(std::ostream &os, const columnbv &columnbv);

    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr,
                     const std::string& name = "");
    void load(std::istream &in);
};


#endif //RLPBWT_COLUMNBV_H
