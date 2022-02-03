//
// Created by dlcgold on 03/02/22.
//

#ifndef RLPBWT_PANEL_RA_H
#define RLPBWT_PANEL_RA_H

#include <sdsl/bit_vectors.hpp>
#include <ostream>


class panel_ra {
public:

    unsigned int h{};
    unsigned int w{};

    sdsl::bit_vector panel;

    panel_ra(unsigned int h, unsigned int w);
    panel_ra();

    char getElem(unsigned int i, unsigned int j) const;

    friend std::ostream &operator<<(std::ostream &os, const panel_ra &ra);
};


#endif //RLPBWT_PANEL_RA_H
