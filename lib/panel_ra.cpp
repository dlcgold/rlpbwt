//
// Created by dlcgold on 03/02/22.
//

#include "../include/panel_ra.h"


panel_ra::panel_ra(unsigned int h, unsigned int w) : h(h), w(w) {
    this->panel = std::vector<sdsl::bit_vector>(w, sdsl::bit_vector(h, 0));
}

char panel_ra::getElem(unsigned int i, unsigned int j) const {
    if (this->panel[j][i]) {
        return '1';
    } else {
        return '0';
    }
}

std::ostream &operator<<(std::ostream &os, const panel_ra &ra) {
    os << "h: " << ra.h << " w: " << ra.w << "\n";
    for (unsigned int i = 0; i < ra.h; ++i) {
        for (unsigned int j = 0; j < ra.w; ++j) {
            os << ra.getElem(i, j) << " ";
        }
        os << "\n";
    }
    return os;
}

panel_ra::panel_ra() = default;


bool
panel_ra::lceToRCheck(unsigned int col, unsigned int f_r, unsigned int s_r,
                      unsigned int length) const {
    unsigned int start = 0;
    if ((int) col - ((int) length - 1) >= 0) {
        start = col - (length - 1);
    }
    bool check = true;
    for (unsigned int i = start; i <= col; i++) {
        if (this->getElem(f_r, i) != this->getElem(s_r, i)) {
            check = false;
            break;
        }
    }
    return check;
}

unsigned int
panel_ra::lceToR(unsigned int col, unsigned int f_r, unsigned int s_r) const {
    int tmp_col = (int) col;
    bool extend = true;
    while (extend) {
        if (this->getElem(f_r, tmp_col) != this->getElem(s_r, tmp_col)) {
            extend = false;
        }
        if (tmp_col == 0) {
            extend = false;
            tmp_col--;
        }
        if (extend) {
            tmp_col--;
        }
    }

    return col - tmp_col;
}

size_t panel_ra::serialize(std::ostream &out, sdsl::structure_tree_node *v,
                           const std::string &name) {
    sdsl::structure_tree_node *child =
            sdsl::structure_tree::add_child(v, name,
                                            sdsl::util::class_name(
                                                    *this));
    size_t written_bytes = 0;
    out.write((char *) &this->h, sizeof(this->h));
    written_bytes += sizeof(this->h);

    out.write((char *) &this->w, sizeof(this->w));
    written_bytes += sizeof(this->w);
    for (const auto &bv: panel) {
        bv.serialize(out, child, "panel");
    }
    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

void panel_ra::load(std::istream &in) {
    in.read((char *) &this->h, sizeof(this->h));
    in.read((char *) &this->w, sizeof(this->w));
    this->panel = std::vector<sdsl::bit_vector>(w, sdsl::bit_vector(h, 0));
    for (unsigned int i = 0; i < this->w; i++) {
        this->panel[i].load(in);
    }
}


unsigned long long panel_ra::size_in_bytes(bool verbose) {
    unsigned long long size = 0;

    for (const auto &bv: this->panel) {
        size += sdsl::size_in_bytes(bv);
    }
    if(verbose){
        std::cout << "panel: " << size << " bytes\n";
    }
    size += (sizeof(unsigned int) * 2);

    return size;
}

double panel_ra::size_in_mega_bytes(bool verbose) {
    double size = 0;
    double to_mega = ((double) 1 / (double) 1024) / (double) 1024;

    for (const auto &bv: this->panel) {
        size += sdsl::size_in_mega_bytes(bv);
    }
    if(verbose){
        std::cout << "panel: " << size << " megabytes\n";
    }
    size += (sizeof(unsigned int) * 2 * to_mega);
    return size;
}

panel_ra::~panel_ra() = default;
