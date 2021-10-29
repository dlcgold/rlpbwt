//
// Created by dlcgold on 29/10/21.
//

#ifndef RLPBWT_UTILS_H
#define RLPBWT_UTILS_H

#include <string>
#include <vector>
#include <climits>
#include "pbwt_column.h"

void update(std::string, std::vector<int>&, std::vector<int>&);
pbwt_column build_column(std::string, std::vector<int>, std::vector<int>);
void build_next_perm(std::vector<pbwt_column>&, int);
#endif //RLPBWT_UTILS_H
