//
// Created by dlcgold on 28/10/21.
//

#include "../include/column.h"

#include <utility>

column::column(bool zeroFirst, std::vector<rlrow> rows,
               unsigned int count0)
        : zero_first(zeroFirst), count_0(count0), rows(std::move(rows)) {}

column::column()
        : zero_first(), count_0(), rows() {}

std::ostream &operator<<(std::ostream &os, const column &column) {
    std::string value;
    if (column.zero_first) {
        value = "0";
    } else {
        value = "1";
    }
    os << "first value: " << value << "\n";
    for (const auto &e: column.rows) {
        os << e << "\n";
    }
    return os;
}

size_t column::serialize(std::ostream &out, sdsl::structure_tree_node *v,
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

    std::vector<unsigned int> first(this->rows.size());
    std::vector<unsigned int> second(this->rows.size());
    for (unsigned int i = 0; i < this->rows.size(); i++) {
        first[i] = rows[i].p;
        second[i] = rows[i].uv;
    }
    written_bytes += my_serialize(first, out, child, "rows_p");
    written_bytes += my_serialize(second, out, child, "rows_uv");
    written_bytes += my_serialize(this->lcp, out, child, "lcp");
    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

void column::load(std::istream &in) {
    in.read((char *) &this->zero_first, sizeof(this->zero_first));
    in.read((char *) &this->count_0, sizeof(this->count_0));
    std::vector<unsigned int> first;
    std::vector<unsigned int> second;
    my_load(first, in);
    my_load(second, in);
    std::vector<rlrow> pairs;
    for (unsigned int i = 0; i < first.size(); i++) {
        pairs.emplace_back(rlrow(first[i], second[i]));
    }
    this->rows = pairs;
    pairs.clear();
    my_load(this->lcp, in);
}

column::~column() = default;
