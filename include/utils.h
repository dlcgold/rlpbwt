//
// Created by dlcgold on 29/10/21.
//

#ifndef RLPBWT_UTILS_H
#define RLPBWT_UTILS_H

#include <string>
#include <vector>
#include <climits>
#include "pbwt_column.h"

void update(std::string, std::vector<unsigned int>&, std::vector<unsigned int>&);
pbwt_column build_column(std::string, std::vector<unsigned int>);
void build_next_perm(std::vector<pbwt_column>&, unsigned int);
char get_next_char(bool, unsigned int);

#endif //RLPBWT_UTILS_H
