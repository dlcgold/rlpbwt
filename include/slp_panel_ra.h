//
// Created by dlcgold on 09/02/22.
//

#ifndef RLPBWT_SLP_PANEL_RA_H
#define RLPBWT_SLP_PANEL_RA_H

#include <fstream>
#include <filesystem>
#include <string>
#include <sdsl/int_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include "exceptions.h"
#include "SelfShapedSlp.hpp"
#include "DirectAccessibleGammaCode.hpp"
#include "SelectType.hpp"

/**
 * @brief class to represent panel with logarithmic random access using an SLP
 */
class slp_panel_ra {
private:
    /**
     * @brief slp filename
     */
  //const char *slp_file{};
  std::string slp_file;
public:
    /**
     * @brief height of the panel
     */
    unsigned int h{};

    /**
    * @brief width of the panel
    */
    unsigned int w{};

    using SelSd = SelectSdvec<>;
    using DagcSd = DirectAccessibleGammaCode<SelSd>;
    using shaped_slp_t = SelfShapedSlp<uint32_t, DagcSd, DagcSd, SelSd>;

    /**
     * @brief panel as shaped SLP
     */
    shaped_slp_t panel;

    /**
     * @brief constructor from an height and a width
     * @param filename slp filename
     * @param h
     * @param w
     */
    slp_panel_ra(const char *filename, unsigned int h, unsigned int w);

    /**
    * @brief default constructor
    */
    slp_panel_ra();

    /**
     * @brief default destructor
     */
     virtual ~slp_panel_ra() = default;
  //~slp_panel_ra();

    friend std::ostream &operator<<(std::ostream &os, const slp_panel_ra &ra);

    /**
    * @brief function to return an element of the panel
    * @param i row index
    * @param j col index
    * @return value (0 or 1) as char
    */
    char getElem(unsigned int i, unsigned int j) const;

    /**
     * @brief function to obtain size in bytes of the panel structure
     * @param verbose bool for extra prints
     * @return size in bytes
    */
    unsigned long long size_in_bytes(bool verbose = false);

    /**
     * @brief function to obtain size in megabytes of the panel structure
     * @param verbose bool for extra prints
     * @return size in megabytes
     */
    double size_in_mega_bytes(bool verbose = false);

    /**
     * @brief function to serialize the panel structure
     * @param out std::ostream object to stream the serialization
     * @return size of the serialization
     */
    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr,
                     const std::string &name = "");

    /**
     * @brief function to load the panel object
     * @param in std::istream object from which load the panel
     * structure object
     * @param slp_filename slp filename
     */
    void load(std::istream &in, const char *slp_filename);
};


#endif //RLPBWT_SLP_PANEL_RA_H
