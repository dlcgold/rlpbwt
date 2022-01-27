//
// Created by dlcgold on 28/01/22.
//

#include "../include/columnbv.h"

#include <utility>


columnbv::columnbv() = default;

columnbv::~columnbv() = default;

columnbv::columnbv(bool zeroFirst, unsigned int count0, const sdsl::bit_vector& runs,
                   const sdsl::bit_vector& u, const sdsl::bit_vector& v,
                   sdsl::int_vector<> lcp) : zero_first(zeroFirst),
                                                    count_0(count0),
                                                    lcp(std::move(lcp)) {
    this->runs = sdsl::sd_vector<>(runs);
    //std::cout << this->runs;
    //this->rank_runs = sdsl::sd_vector<>::rank_1_type(&this->runs);
    //this->select_runs = sdsl::sd_vector<>::select_1_type(&this->runs);
    this->u = sdsl::sd_vector<>(u);
    //this->rank_u = sdsl::sd_vector<>::rank_1_type(&this->u);
    //this->select_u = sdsl::sd_vector<>::select_1_type(&this->u);
    this->v = sdsl::sd_vector<>(v);
    //this->rank_v = sdsl::sd_vector<>::rank_1_type(&this->v);
    //this->select_v = sdsl::sd_vector<>::select_1_type(&this->v);

}
