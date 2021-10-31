//
// Created by dlcgold on 28/10/21.
//

#include "../include/pbwt_column.h"

#include <utility>

pbwt_column::pbwt_column(bool zeroFirst, std::vector<pbwt_rlrow> rows)
        : zero_first(zeroFirst), rows(std::move(rows)) {}
pbwt_column::pbwt_column()
        : zero_first(), rows() {}

std::ostream &operator<<(std::ostream &os, const pbwt_column &column) {
    //os << "zero_first: " << column.zero_first << " rows: " << column.rows;
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

pbwt_column::~pbwt_column() = default;
