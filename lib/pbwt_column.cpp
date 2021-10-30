//
// Created by dlcgold on 28/10/21.
//

#include "../include/pbwt_column.h"

pbwt_column::pbwt_column(bool zeroFirst, const std::vector<pbwt_rlrow> &rows)
        : zero_first(zeroFirst), rows(rows) {}
pbwt_column::pbwt_column()
        : zero_first(true), rows() {}
pbwt_column::~pbwt_column() = default;
