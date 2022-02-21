//
// Created by dlcgold on 21/02/22.
//

#include "../include/column_naive.h"


column_naive::column_naive(bool zeroFirst, unsigned int count0,
                           sdsl::int_vector<> p,
                           sdsl::int_vector<> uv,
                           sdsl::int_vector<> lcp) : zero_first(
        zeroFirst), count_0(count0), p(std::move(p)), uv(std::move(uv)),
                                                     lcp(std::move(lcp)) {}


column_naive::column_naive() = default;

column_naive::~column_naive() = default;

std::ostream &operator<<(std::ostream &os, const column_naive &naive) {
    std::string value;
    if (naive.zero_first) {
        value = "0";
    } else {
        value = "1";
    }
    os << "first value: " << value << "\n";
    os << "c value: " << naive.count_0 << "\n";
    os << "lcp: " << naive.lcp << "\n";
    for (unsigned int i = 0; i < naive.p.size(); i++) {
        os << naive.p[i] << "\t" << naive.uv[i] << "\n";
    }
    return os;
}

size_t column_naive::serialize(std::ostream &out, sdsl::structure_tree_node *v,
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

    this->p.serialize(out, child, "p");
    this->uv.serialize(out, child, "uv");
    this->lcp.serialize(out, child, "lcp");
    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

void column_naive::load(std::istream &in) {
    in.read((char *) &this->zero_first, sizeof(this->zero_first));
    in.read((char *) &this->count_0, sizeof(this->count_0));
    this->p.load(in);
    this->uv.load(in);
    this->lcp.load(in);
}

unsigned long long column_naive::size_in_bytes(bool verbose) const {
    unsigned long long size = 0;
    size += sizeof(bool);
    size += sizeof(unsigned int);
    auto p_size = sdsl::size_in_bytes(this->p);
    auto uv_size = sdsl::size_in_bytes(this->uv);
    auto lcp_size = sdsl::size_in_bytes(this->lcp);
    size += p_size;
    size += uv_size;
    size += lcp_size;
    if (verbose) {
        std::cout << "p: " << p_size << " bytes\n";
        std::cout << "uv: " << uv_size << " bytes\n";
        std::cout << "lcp: " << lcp_size << " bytes\n";
    }
    return size;
}

double column_naive::size_in_mega_bytes(bool verbose) const {
    double size = 0;
    double to_mega = ((double) 1 / (double) 1024) / (double) 1024;
    size += (double) (sizeof(bool) * to_mega);
    size += (double) (sizeof(unsigned int) * to_mega);
    auto p_size = sdsl::size_in_mega_bytes(this->p);
    auto uv_size = sdsl::size_in_mega_bytes(this->uv);
    auto lcp_size = sdsl::size_in_mega_bytes(this->lcp);
    size += p_size;
    size += uv_size;
    size += lcp_size;
    if (verbose) {
        std::cout << "p: " << p_size << " megabytes\n";
        std::cout << "uv: " << uv_size << " megabytes\n";
        std::cout << "lcp: " << lcp_size << " megabytes\n";
    }
    return size;
}
