//
// Created by dlcgold on 03/02/22.
//

#ifndef RLPBWT_PANEL_RA_H
#define RLPBWT_PANEL_RA_H

#include <sdsl/bit_vectors.hpp>
#include <ostream>

/**
 * @brief class to represent panel with constant random access using vector of
 * bit vectors
 */
class panel_ra {
public:
    /**
     * @brief height of the panel
     */
    unsigned int h{};

    /**
    * @brief width of the panel
    */
    unsigned int w{};

    /**
     * @brief panel of bitvectors
     */
    std::vector<sdsl::bit_vector> panel;

    /**
     * @brief constructor from an height and a width
     * @param h
     * @param w
     */
    panel_ra(unsigned int h, unsigned int w);

    /**
     * @brief default constructor
     */
    panel_ra();

    /**
     * @brief default destructor
     */
    virtual ~panel_ra();

    /**
     * @brief function to return an element of the panel
     * @param i row index
     * @param j col index
     * @return value (0 or 1) as char
     */
    char getElem(unsigned int i, unsigned int j) const;

    friend std::ostream &operator<<(std::ostream &os, const panel_ra &ra);


    /**
     * @brief function to compute longest common extension between two rows
     * starting from a column
     * @param col starting column
     * @param f_r first row
     * @param s_r second row
     * @return size of the lce
     */
    unsigned int
    lceToR(unsigned int col, unsigned int f_r, unsigned int s_r) const;

    /**
     * @brief function to check if longest common extension between two rows
     * starting from a column is equal or greater than a bound
     * @param col starting column
     * @param f_r first row
     * @param s_r second row
     * @param length length bound of the lce
     * @return true if lce is equal or greater than the bound
     */
    bool lceToRCheck(unsigned int col, unsigned int f_r, unsigned int s_r,
                     unsigned int length) const;


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
     */
    void load(std::istream &in);
};


#endif //RLPBWT_PANEL_RA_H
