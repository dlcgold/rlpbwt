//
// Created by dlcgold on 03/02/22.
//

#include "../include/panel_ra.h"


panel_ra::panel_ra(unsigned int h, unsigned int w) : h(h), w(w) {
    this->panel = sdsl::bit_vector(h * w, 0);
}

char panel_ra::getElem(unsigned int i, unsigned int j) const {
    if (this->panel[j + (i * w)]) {
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
