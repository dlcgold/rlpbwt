//
// Created by dlcgold on 28/10/21.
//

#include "../include/pbwt_column.h"

#include <utility>

pbwt_column::pbwt_column(bool zeroFirst, std::vector<pbwt_rlrow> rows)
        : zero_first(zeroFirst), rows(std::move(rows)) {}
pbwt_column::pbwt_column()
        : zero_first(true), rows() {}
pbwt_column::~pbwt_column() = default;
