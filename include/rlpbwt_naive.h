//
// Created by dlcgold on 21/02/22.
//

#ifndef RLPBWT_RLPBWT_NAIVE_H
#define RLPBWT_RLPBWT_NAIVE_H


#include <vector>
#include <climits>
#include <ostream>
#include "column_naive.h"
#include "matches_naive.h"
#include "utils.h"


/**
 * @brief data structure for naive RLPBWT
 */
class rlpbwt_naive {
public:
    /**
     * @brief vector of naive columns
     */
    std::vector<column_naive> cols;

    /**
     * @brief height of the panel
     */
    unsigned int height{};

    /**
     * @brief width of the panel
     */
    unsigned int width{};

    /**
     * @brief costructor of the naive RLPBWT
     * @param filename file with the original panel
     * @param verbose bool for extra prints
     */
    explicit rlpbwt_naive(const char *filename, bool verbose = false);

    /**
     * @brief default constructor
     */
    rlpbwt_naive();

    /**
     * @brief default destructor
     */
    virtual ~rlpbwt_naive();

    /**
     * @brief function to obtain size in bytes of the naive RLPBWT
     * @param verbose bool for extra prints
     * @return size in bytes
    */
    unsigned long long size_in_bytes(bool verbose = false);

    /**
    * @brief function to obtain size in megabytes of the naive RLPBWT
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
    external_match(const std::string &query, bool verbose = false);

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

    void
    match_tsv_conc(const char *filename, const char *out, bool verbose = false);

    friend std::ostream &
    operator<<(std::ostream &os, const rlpbwt_naive &naive);

    /**
     * function to get the total number of runs in the RLPBWT
     * @return total number of run
     */
    unsigned int get_run_number();

    /**
     * @brief function to serialize the naive RLPBWT
     * @param out std::ostream object to stream the serialization
     * @return size of the serialization
     */
    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr,
                     const std::string &name = "");

    /**
     * @brief function to load the naive RLPBWT object
     * @param in std::istream object from which load the naive RLPBWT
     * structure object
     */
    void load(std::istream &in);

private:

    /**
     * @brief function to build a naive column for a column of the panel
     * @param column column of the panel
     * @param pref current prefix array
     * @param div current divergence array (as lcp array)
     * @return the new naive column
     */
    static column_naive
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
     * @param row_index index of the run
     * @param symbol symbol s
     * @param offset offset to correctly calculate the mapping
     * @param verbose bool for extra print
     * @return the index computed with the lf-mapping
     */
    unsigned int
    lf(unsigned int col_index, unsigned int row_index, char symbol,
       unsigned int offset, bool verbose = false) const;


    /**
     * @brief function to get the run in previous column from which come the
     * current run (a sort of reverse lf mapping)
     * @param col_index index of the column
     * @param index virtual index of the row of the original panel
     * @param verbose bool for extra print
     * @return
     */
    unsigned int
    reverse_lf(unsigned int col_index, unsigned int index, bool verbose) const;

    /**
     * @brief function to map an index to the correct run in a column
     * @param index index to map
     * @param col_index column index
     * @return run index
     */
    unsigned int index_to_run(unsigned int index, unsigned int col_index) const;

    /**
     * @brief trick to extract u and v value from a run in rlpbwt column
     * @param col_index index of the column
     * @param row_index index of the run
     * @return a std::pair with u as first and v as second
     */
    std::pair<unsigned int, unsigned int>
    uvtrick(unsigned int col_index, unsigned int row_index) const;
};


#endif //RLPBWT_RLPBWT_NAIVE_H
