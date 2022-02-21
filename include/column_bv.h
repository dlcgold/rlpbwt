//
// Created by dlcgold on 28/01/22.
//

#ifndef RLPBWT_COLUMN_BV_H
#define RLPBWT_COLUMN_BV_H

#include <vector>
#include <sdsl/int_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include <ostream>

class column_bv {
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

    column_bv(bool zeroFirst, unsigned int count0,
             const sdsl::bit_vector& runs,
             const sdsl::bit_vector& u,
             const sdsl::bit_vector& v,
             sdsl::int_vector<> lcp);

    /**
     * @brief default constructor
     */
    column_bv();

    virtual ~column_bv();

    friend std::ostream &operator<<(std::ostream &os, const column_bv &columnbv);
    unsigned long long size_in_bytes(bool verbose = false) const;
    double size_in_mega_bytes(bool verbose = false) const;
    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr,
                     const std::string& name = "");
    void load(std::istream &in);
};


#endif //RLPBWT_COLUMN_BV_H
