//
// Created by dlcgold on 28/01/22.
//

#include "../include/column_bv.h"

#include <utility>


column_bv::column_bv() = default;

column_bv::~column_bv() = default;

column_bv::column_bv(bool zeroFirst, unsigned int count0,
                   const sdsl::bit_vector &runs,
                   const sdsl::bit_vector &u, const sdsl::bit_vector &v,
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

size_t column_bv::serialize(std::ostream &out, sdsl::structure_tree_node *v,
                           const std::string &name) {
    sdsl::structure_tree_node *child =
            sdsl::structure_tree::add_child(v, name,
                                            sdsl::util::class_name(
                                                    *this));
    size_t written_bytes = 0;

    out.write((char *) &this->zero_first, sizeof(this->zero_first));
    written_bytes += sizeof(this->zero_first);

    out.write((char *) &this->count_0, sizeof(this->count_0));
    written_bytes += sizeof(this->count_0);

    written_bytes += this->runs.serialize(out, child, "runs");
    written_bytes += this->u.serialize(out, child, "u");
    written_bytes += this->v.serialize(out, child, "v");
    written_bytes += this->lcp.serialize(out, child, "lcp");

    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

void column_bv::load(std::istream &in) {
    in.read((char *) &this->zero_first, sizeof(this->zero_first));
    in.read((char *) &this->count_0, sizeof(this->count_0));
    this->runs.load(in);
    this->u.load(in);
    this->v.load(in);
    this->lcp.load(in);
    this->rank_runs = sdsl::sd_vector<>::rank_1_type(&this->runs);
    this->select_runs = sdsl::sd_vector<>::select_1_type(&this->runs);
    this->rank_u = sdsl::sd_vector<>::rank_1_type(&this->u);
    this->select_u = sdsl::sd_vector<>::select_1_type(&this->u);
    this->rank_v = sdsl::sd_vector<>::rank_1_type(&this->v);
    this->select_v = sdsl::sd_vector<>::select_1_type(&this->v);
}

std::ostream &operator<<(std::ostream &os, const column_bv &columnbv) {
    auto yesno = "yes";
    if (!columnbv.zero_first) {
        yesno = "no";
    }
    os << "zero_first: " << yesno << " c: " << columnbv.count_0
       << "\nruns: " << columnbv.runs << "\nu: " << columnbv.u << "\nv: "
       << columnbv.v
       << "\nlcp: " << columnbv.lcp << "\n";
    return os;
}

unsigned long long column_bv::size_in_bytes(bool verbose) const {
    unsigned long long size = 0;
    size += sizeof(bool);
    size += sizeof(unsigned int);
    auto size_run = sdsl::size_in_bytes(this->runs) +
                    sdsl::size_in_bytes(this->rank_runs) +
                    sdsl::size_in_bytes(this->select_runs);
    auto size_u = sdsl::size_in_bytes(this->u) +
                  sdsl::size_in_bytes(this->rank_u) +
                  sdsl::size_in_bytes(this->select_u);
    auto size_v = sdsl::size_in_bytes(this->v) +
                  sdsl::size_in_bytes(this->rank_v) +
                  sdsl::size_in_bytes(this->select_v);
    auto size_lcp = sdsl::size_in_bytes(this->lcp);
    size += size_run;
    size += size_u;
    size += size_v;
    size += size_lcp;
    if (verbose) {
        std::cout << "run: " << size_run << " bytes\n";
        std::cout << "u: " << size_u << " bytes\n";
        std::cout << "v: " << size_v << " bytes\n";
        std::cout << "lcp: " << size_lcp << " bytes\n";
    }
    return size;
}

double column_bv::size_in_mega_bytes(bool verbose) const {
    double size = 0;
    double to_mega = ((double) 1 / (double) 1024) / (double) 1024;
    size += (double) (sizeof(bool) * to_mega);
    size += (double) (sizeof(unsigned int) * to_mega);
    auto size_run = sdsl::size_in_mega_bytes(this->runs) +
                    sdsl::size_in_mega_bytes(this->rank_runs) +
                    sdsl::size_in_mega_bytes(this->select_runs);
    auto size_u = sdsl::size_in_mega_bytes(this->u) +
                  sdsl::size_in_mega_bytes(this->rank_u) +
                  sdsl::size_in_mega_bytes(this->select_u);
    auto size_v = sdsl::size_in_mega_bytes(this->v) +
                  sdsl::size_in_mega_bytes(this->rank_v) +
                  sdsl::size_in_mega_bytes(this->select_v);
    auto size_lcp = sdsl::size_in_mega_bytes(this->lcp);
    size += size_run;
    size += size_u;
    size += size_v;
    size += size_lcp;
    if (verbose) {
        std::cout << "run: " << size_run << " megabytes\n";
        std::cout << "u: " << size_u << " megabytes\n";
        std::cout << "v: " << size_v << " megabytes\n";
        std::cout << "lcp: " << size_lcp << " megabytes\n";
    }
    return size;
}
