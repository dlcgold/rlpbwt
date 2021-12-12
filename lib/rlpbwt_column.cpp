//
// Created by dlcgold on 28/10/21.
//

#include "../include/rlpbwt_column.h"

#include <utility>

rlpbwt_column::rlpbwt_column(bool zeroFirst, std::vector<rlpbwt_rlrow> rows,
                             unsigned int count0)
        : zero_first(zeroFirst), rows(std::move(rows)), count_0(count0) {}
rlpbwt_column::rlpbwt_column()
        : zero_first(), rows() {}

std::ostream &operator<<(std::ostream &os, const rlpbwt_column &column) {
    std::string value;
    if(column.zero_first){
        value = "0";
    }else{
        value = "1";
    }
    os << "first value: " << value << "\n";
    for(const auto& e: column.rows){
        os << e << "\n";
    }
    return os;
}

rlpbwt_column::~rlpbwt_column() = default;
