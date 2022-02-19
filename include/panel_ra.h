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

    //sdsl::bit_vector panel;
    std::vector<sdsl::bit_vector> panel;

    panel_ra(unsigned int h, unsigned int w);

    panel_ra();

    virtual ~panel_ra();

    char getElem(unsigned int i, unsigned int j) const;

    friend std::ostream &operator<<(std::ostream &os, const panel_ra &ra);

    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr,
                     const std::string &name = "");

    unsigned int lceToR(unsigned int col, unsigned int f_r, unsigned int s_r) const;
    bool lceToRCheck(unsigned int col, unsigned int f_r, unsigned int s_r,
                     unsigned int length) const;
    void load(std::istream &in);
};


#endif //RLPBWT_PANEL_RA_H
