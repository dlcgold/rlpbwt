//
// Created by dlcgold on 28/01/22.
//

#ifndef RLPBWT_RLPBWT_BV_H
#define RLPBWT_RLPBWT_BV_H

#include <vector>
#include <string>
#include <iostream>
#include <climits>
#include <cassert>
#include <utility>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <list>
#include <sdsl/bit_vectors.hpp>
#include "utils.h"
#include "exceptions.h"
#include "column_bv.h"
#include "matches_naive.h"

/**
 * @brief data structure for RLPBWT using sparse bitvectors
 */
class rlpbwt_bv {
public:
    /**
     * @brief vector of bitvectors columns
     */
    std::vector<column_bv> cols;

    /**
     * @brief height of the panel
     */
    unsigned int height{};

    /**
     * @brief width of the panel
     */
    unsigned int width{};

    /**
     * @brief costructor of the bitvectors RLPBWT
     * @param filename file with the original panel
     * @param verbose bool for extra prints
     */
    explicit rlpbwt_bv(const char *filename, bool verbose = false);

    /**
     * @brief default constructor
     */
    rlpbwt_bv();

    /**
     * @brief function to obtain size in bytes of the bitvectors RLPBWT
     * @param verbose bool for extra prints
     * @return size in bytes
    */
    unsigned long long size_in_bytes(bool verbose = false);

    /**
     * @brief function to obtain size in megabytes of the bitvectors RLPBWT
     * @param verbose bool for extra prints
     * @return size in megabytes
     */
    double size_in_mega_bytes(bool verbose = false);

    /**
     * @brief function to produce naive matches from a query haplotype
     * @param query std::string with the query
     * @param verbose bool for extra prints
     * @return a matches_naive object with naive matches
     */
    matches_naive
    external_match(const std::string &query, unsigned int min_len = 1,
                   bool verbose = false);

    //void external_match_vcf(const char *filename, unsigned int min_len = 1,
    //                        bool verbose = false);


    /**
     * @brief function to compute queries from a tsv file and output them
     * on a file
     * @param filename queries file
     * @param out output file
     * @param verbose bool for extra prints
     */
    void
    match_tsv(const char *filename, const char *out,
              bool verbose = false);

    /**
     * @brief function to compute queries from a transposed tsv file and output
     * them on a file
     * @param filename queries file
     * @param out output file
     * @param verbose bool for extra prints
     */
    void
    match_tsv_tr(const char *filename, const char *out, bool verbose = false);

    /**
     * @brief function to serialize the bitvector RLPBWT
     * @param out std::ostream object to stream the serialization
     * @return size of the serialization
     */
    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr,
                     const std::string &name = "");

    /**
     * @brief function to load the bitvector RLPBWT object
     * @param in std::istream object from which load the bitvector RLPBWT
     * structure object
     */
    void load(std::istream &in);

private:

    /**
     * @brief function to build a sparse bitvectors column for a column of the
     * panel
     * @param column column of the panel
     * @param pref current prefix array
     * @param div current divergence array (as lcp array)
     * @return the new naive column
     */
    static column_bv
    build_column(std::string &column, std::vector<unsigned int> &pref,
                 sdsl::int_vector<> &div);

    /**
     * @brief function to update prefix/divergence (lcp) arrays as in Durbin
     * @param column column of the panel
     * @param pref current prefix array
     * @param div current divergence array (as lcp array)
     */
    static void
    update(std::string &column, std::vector<unsigned int> &pref,
           sdsl::int_vector<> &div);

    /**
     * @brief function to compute the lf mapping, w(i, s) function in Durbin
     * @param col_index index of the column
     * @param index index of the row
     * @param symbol symbol s
     * @param verbose bool for extra print
     * @return the index computed with the lf-mapping
     */
    unsigned int
    lf(unsigned int col_index, unsigned int index, char symbol,
       bool verbose = false) const;


    /**
     * @brief trick to extract u and v value from a run in rlpbwt column
     * @param col_index index of the column
     * @param index virtual index of the row of the original panel
     * @return a std::pair with u as first and v as second
     */
    std::pair<unsigned int, unsigned int>
    uvtrick(unsigned int col_index, unsigned int index) const;

    /**
     * @brief function to get the run in previous column from which come the
     * current run (a sort of reverse lf mapping)
     * @param col_index index of the column
     * @param index virtual index of the row of the original panel
     * @param verbose bool for extra print
     * @return
     */
    unsigned int
    reverse_lf(unsigned int col_index, unsigned int index,
               bool verbose) const;
};


#endif //RLPBWT_RLPBWT_BV_H
