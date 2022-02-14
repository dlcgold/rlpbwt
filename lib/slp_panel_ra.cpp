//
// Created by dlcgold on 09/02/22.
//

#include "../include/slp_panel_ra.h"

slp_panel_ra::slp_panel_ra(const char *filename, unsigned int h, unsigned int w)
        : h(h),
          w(w) {
    std::ifstream in(filename);
    if (in.is_open()) {
        this->panel.load(in);
        in.close();
    } else {
        throw FileNotFoundException{};
    }
}


slp_panel_ra::slp_panel_ra() = default;

char slp_panel_ra::getElem(unsigned int i, unsigned int j) const {
    unsigned rev_col = (this->w - 1) - j;
    return this->panel.charAt(rev_col + (i * this->w));
}

std::ostream &operator<<(std::ostream &os, const slp_panel_ra &ra) {
    os << "h: " << ra.h << " w: " << ra.w << "\n";
    for (unsigned int i = 0; i < ra.h; i++) {
        for (unsigned int j = 0; j < ra.w; j++) {
            os << ra.getElem(i, j) << " ";
        }
        os << "\n";
    }
    return os;
}

/*
size_t slp_panel_ra::serialize(std::ostream &out, sdsl::structure_tree_node *v,
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
    written_bytes += this->panel.serialize(out, child, "panel");
    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

void slp_panel_ra::load(std::istream &in) {
    in.read((char *) &this->h, sizeof(this->h));
    in.read((char *) &this->w, sizeof(this->w));
    this->panel.load(in);
}
 */

slp_panel_ra::~slp_panel_ra() = default;
