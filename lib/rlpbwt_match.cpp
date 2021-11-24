//
// Created by dlcgold on 10/11/21.
//

#include "../include/rlpbwt_match.h"

#include <utility>

rlpbwt_match::rlpbwt_match(unsigned int begin, unsigned int length,
             std::vector<unsigned int> rows) : begin(begin),
                                                      length(length),
                                                      rows(std::move(rows)) {}

rlpbwt_match::~rlpbwt_match() {

}
