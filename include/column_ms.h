//
// Created by dlcgold on 02/02/22.
//

#ifndef RLPBWT_COLUMN_MS_H
#define RLPBWT_COLUMN_MS_H

#include <vector>
#include <ostream>
#include <utility>
#include <sdsl/int_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include "utils.h"

/**
 * @brief class to represent a column that supports matching statistics
 */
class column_ms {
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

    /**
     * @brief sparse bit vector representing run
     */
    sdsl::sd_vector<> runs;

    /**
     * @brief rank support for runs sparse bitvector
     */
    sdsl::sd_vector<>::rank_1_type rank_runs;

    /**
     * @brief select support for runs sparse bitvector
     */
    sdsl::sd_vector<>::select_1_type select_runs;

    /**
     * @brief sparse bit vector representing u (as in Durbin)
     */
    sdsl::sd_vector<> u;

    /**
    * @brief rank support for u sparse bitvector
    */
    sdsl::sd_vector<>::rank_1_type rank_u;

    /**
    * @brief select support for u sparse bitvector
    */
    sdsl::sd_vector<>::select_1_type select_u;

    /**
     * @brief sparse bit vector representing v (as in Durbin)
     */
    sdsl::sd_vector<> v;

    /**
     * @brief rank support for v sparse bitvector
    */
    sdsl::sd_vector<>::rank_1_type rank_v;

    /**
     * @brief select support for v sparse bitvector
     */
    sdsl::sd_vector<>::select_1_type select_v;

    /**
     * @brief sparse bit vector representing the thresholds (min lcp in a run)
     */
    sdsl::sd_vector<> thr;

    /**
    * @brief rank support for thresholds sparse bitvector
    */
    sdsl::sd_vector<>::rank_1_type rank_thr;

    /**
     * @brief select support for thresholds sparse bitvector
     */
    sdsl::sd_vector<>::select_1_type select_thr;

    /**
     * @brief sdsl compressed int vector for runs-head prefix array values
     */
    sdsl::int_vector<> sample_beg;

    /**
     * @brief sdsl compressed int vector for runs-tail prefix array values
     */
    sdsl::int_vector<> sample_end;

    /**
     * @brief constructor for a column supporting matching statistics
     * use
     *
     * @param zeroFirst
     * @param count0
     * @param runs
     * @param u
     * @param v
     * @param thr
     * @param sample_beg
     * @param sample_end
     */
    column_ms(bool zeroFirst, unsigned int count0,
               const sdsl::bit_vector &runs,
               const sdsl::bit_vector &u,
               const sdsl::bit_vector &v,
               const sdsl::bit_vector &thr,
               sdsl::int_vector<> sample_beg,
               sdsl::int_vector<> sample_end);

    /**
     * @brief default constructor
     */
    column_ms();

    /**
     * @brief default destructor
     */
    virtual ~column_ms();

    friend std::ostream &operator<<(std::ostream &os, const column_ms &thr);

    /**
     * @brief function to obtain size in bytes of the column
     * @param verbose bool for extra prints
     * @return size in bytes
     */
    unsigned long long size_in_bytes(bool verbose = false) const;

    /**
     * @brief function to obtain size in megabytes of the column
     * @param verbose bool for extra prints
     * @return size in megabytes
     */
    double size_in_mega_bytes(bool verbose = false) const;

    /**
     * @brief function to serialize the column object
     * @param out std::ostream object to stream the serialization
     * @return size of the serialization
     */
    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr,
                     const std::string &name = "");

    /**
     * @brief function to load the column object
     * @param in std::istream object from which load the column
     */
    void load(std::istream &in);
};


#endif //RLPBWT_COLUMN_MS_H
