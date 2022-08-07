//
// Created by dlcgold on 09/02/22.
//

#include "../include/slp_panel_ra.h"

slp_panel_ra::slp_panel_ra(const char *filename, unsigned int h, unsigned int w)
        : slp_file(filename),
          h(h),
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
    auto w_l = (unsigned long long int) this->w;
    auto i_l = (unsigned long long int) i;
    auto j_l = (unsigned long long int) j;
    auto rev_col = (w_l - (unsigned long long int) 1) - j_l;
    unsigned long long int pos = rev_col + (i_l * w_l);
    return this->panel.charAt(pos);
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
    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

void slp_panel_ra::load(std::istream &in, const char *slp_filename) {
    in.read((char *) &this->h, sizeof(this->h));
    in.read((char *) &this->w, sizeof(this->w));
    std::ifstream slp_in;
    slp_in.open(slp_filename);
    this->panel.load(slp_in);
    slp_in.close();
    this->slp_file = slp_filename;
}

unsigned long long slp_panel_ra::size_in_bytes(bool verbose) {
    unsigned long long size = 0;
    std::filesystem::path slp{this->slp_file};
    size += std::filesystem::file_size(slp);
    if (verbose) {
        std::cout << "slp: " << size << " bytes\n";
    }
    size += (sizeof(unsigned int) * 2);

    return size;
}

double slp_panel_ra::size_in_mega_bytes(bool verbose) {
    double size = 0;
    std::filesystem::path slp{this->slp_file};
    double to_mega = ((double) 1 / (double) 1024) / (double) 1024;
    size += ((double) std::filesystem::file_size(slp) * to_mega);
    if (verbose) {
        std::cout << "slp: " << size << " megabytes\n";
    }
    size += (sizeof(unsigned int) * 2);
    return size * to_mega;
}


// slp_panel_ra::~slp_panel_ra() {
//   delete slp_file;
// };
