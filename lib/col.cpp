//
// Created by dlcgold on 28/10/21.
//

#include "../include/col.h"

#include <utility>

col::col(bool zeroFirst, std::vector<rlrow> rows,
         unsigned int count0)
        : zero_first(zeroFirst), count_0(count0), rows(std::move(rows)) {}
col::col()
        : zero_first(), count_0(), rows() {}

std::ostream &operator<<(std::ostream &os, const col &column) {
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

col::~col() = default;
