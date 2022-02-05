//
// Created by dlcgold on 02/02/22.
//

#include "../include/column_thr.h"

column_thr::column_thr(bool zeroFirst, unsigned int count0,
                       const sdsl::bit_vector &runs, const sdsl::bit_vector &u,
                       const sdsl::bit_vector &v, const sdsl::bit_vector &thr,
                       sdsl::int_vector<> sample_beg,
                       sdsl::int_vector<> sample_end) : zero_first(zeroFirst),
                                                        count_0(count0),
                                                        sample_beg(std::move(
                                                                sample_beg)),
                                                        sample_end(std::move(
                                                                sample_end)) {
    this->runs = sdsl::sd_vector<>(runs);
    this->u = sdsl::sd_vector<>(u);
    this->v = sdsl::sd_vector<>(v);
    this->thr = sdsl::sd_vector<>(thr);
}

std::ostream &operator<<(std::ostream &os, const column_thr &thr) {
    auto yesno = "yes";
    if (!thr.zero_first) {
        yesno = "no";
    }
    os << "zero_first: " << yesno << " c: " << thr.count_0
       << "\nruns: " << thr.runs << "\nu: " << thr.u << "\nv: " << thr.v
       << "\nthr: " << thr.thr << "\n";
    for (unsigned int i = 0; i < thr.sample_beg.size(); i++) {
        std::cout << "(" << thr.sample_beg[i] << ", " << thr.sample_end[i] << ") ";
    }
    std::cout << "\n";
    return os;
}

size_t column_thr::serialize(std::ostream &out, sdsl::structure_tree_node *v,
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
    written_bytes += this->thr.serialize(out, child, "thr");

    written_bytes += this->sample_beg.serialize(out, child, "s_beg");
    written_bytes += this->sample_end.serialize(out, child, "s_end");

    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

void column_thr::load(std::istream &in) {
    in.read((char *) &this->zero_first, sizeof(this->zero_first));
    in.read((char *) &this->count_0, sizeof(this->count_0));
    this->runs.load(in);
    this->u.load(in);
    this->v.load(in);
    this->thr.load(in);
    this->sample_beg.load(in);
    this->sample_end.load(in);
    this->rank_runs = sdsl::sd_vector<>::rank_1_type(&this->runs);
    this->select_runs = sdsl::sd_vector<>::select_1_type(&this->runs);
    this->rank_u = sdsl::sd_vector<>::rank_1_type(&this->u);
    this->select_u = sdsl::sd_vector<>::select_1_type(&this->u);
    this->rank_v = sdsl::sd_vector<>::rank_1_type(&this->v);
    this->select_v = sdsl::sd_vector<>::select_1_type(&this->v);
    this->rank_thr = sdsl::sd_vector<>::rank_1_type(&this->thr);
    this->select_thr = sdsl::sd_vector<>::select_1_type(&this->thr);
}


column_thr::~column_thr() = default;

column_thr::column_thr() = default;

