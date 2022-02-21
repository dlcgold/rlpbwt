//
// Created by dlcgold on 09/02/22.
//

#ifndef RLPBWT_SLP_PANEL_RA_H
#define RLPBWT_SLP_PANEL_RA_H

#include <fstream>
#include <filesystem>
#include <sdsl/int_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include "exceptions.h"
#include "SelfShapedSlp.hpp"
#include "DirectAccessibleGammaCode.hpp"
#include "SelectType.hpp"

class slp_panel_ra {
private:
    const char *slp_file{};
public:

    unsigned int h{};
    unsigned int w{};

    using SelSd = SelectSdvec<>;
    using DagcSd = DirectAccessibleGammaCode<SelSd>;
    using shaped_slp_t = SelfShapedSlp<uint32_t, DagcSd, DagcSd, SelSd>;
    shaped_slp_t panel;

    slp_panel_ra(const char *filename, unsigned int h, unsigned int w);

    slp_panel_ra();

    virtual ~slp_panel_ra();

    friend std::ostream &operator<<(std::ostream &os, const slp_panel_ra &ra);

    char getElem(unsigned int i, unsigned int j) const;

    unsigned long long size_in_bytes(bool verbose = false);
    double size_in_mega_bytes(bool verbose = false);

    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr,
                     const std::string &name = "");

    void load(std::istream &in, const char *slp_filename);
};


#endif //RLPBWT_SLP_PANEL_RA_H
