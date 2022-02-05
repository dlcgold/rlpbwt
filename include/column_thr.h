//
// Created by dlcgold on 02/02/22.
//

#ifndef RLPBWT_COLUMN_THR_H
#define RLPBWT_COLUMN_THR_H

#include <vector>
#include <ostream>
#include <utility>
#include <sdsl/int_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include "utils.h"

class column_thr {
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

    sdsl::sd_vector<> thr;
    sdsl::sd_vector<>::rank_1_type rank_thr;
    sdsl::sd_vector<>::select_1_type select_thr;

    //std::vector<std::pair<unsigned int, unsigned int>> rows;
    sdsl::int_vector<> sample_beg;
    sdsl::int_vector<> sample_end;

    column_thr(bool zeroFirst, unsigned int count0,
               const sdsl::bit_vector &runs,
               const sdsl::bit_vector &u,
               const sdsl::bit_vector &v,
               const sdsl::bit_vector &thr,
               sdsl::int_vector<> sample_beg,
               sdsl::int_vector<> sample_end);

    /**
     * @brief default constructor
     */
    column_thr();

    virtual ~column_thr();

    friend std::ostream &operator<<(std::ostream &os, const column_thr &thr);
    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr,
                     const std::string& name = "");
    void load(std::istream &in);
};


#endif //RLPBWT_COLUMN_THR_H
