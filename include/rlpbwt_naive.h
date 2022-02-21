//
// Created by dlcgold on 21/02/22.
//

#ifndef RLPBWT_RLPBWT_NAIVE_H
#define RLPBWT_RLPBWT_NAIVE_H


#include <vector>
#include <climits>
#include <ostream>
#include "column_naive.h"
#include "utils.h"
#include "matches_naive.h"

class rlpbwt_naive {
public:
    std::vector<column_naive> cols;
    unsigned int height{};
    unsigned int width{};

    explicit rlpbwt_naive(const char *filename, bool verbose = false);

    rlpbwt_naive();

    virtual ~rlpbwt_naive();

    static column_naive
    build_column(std::string &column, std::vector<unsigned int> &pref,
                 sdsl::int_vector<> &div);

    static void
    update(std::string &column, std::vector<unsigned int> &pref,
           sdsl::int_vector<> &div);
    unsigned long long size_in_bytes(bool verbose = false);
    double size_in_mega_bytes(bool verbose = false);
    matches_naive
    external_match(const std::string &query, bool verbose = false);

    void
    match_tsv(const char *filename, const char *out,
              bool verbose = false);

    void
    match_tsv_tr(const char *filename, const char *out, bool verbose = false);

    friend std::ostream &
    operator<<(std::ostream &os, const rlpbwt_naive &naive);


private:

    unsigned int
    lf(unsigned int col_index, unsigned int row_index, char symbol,
       unsigned int offset, bool verbose = false) const;

    unsigned int
    reverse_lf(unsigned int col_index, unsigned int index, bool verbose) const;

    unsigned int index_to_run(unsigned int index, unsigned int col_index) const;

    std::pair<unsigned int, unsigned int>
    uvtrick(unsigned int col_index, unsigned int row_index) const;
};


#endif //RLPBWT_RLPBWT_NAIVE_H
