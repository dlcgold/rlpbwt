//
// Created by dlcgold on 17/02/22.
//

#ifndef RLPBWT_PHI_SUPPORT_H
#define RLPBWT_PHI_SUPPORT_H


#include <vector>
#include <sdsl/int_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include <optional>
#include "panel_ra.h"
#include "slp_panel_ra.h"
#include "column_ms.h"

/**
 * @brief class to rapresent the additional data structure for phi and phi_inv
 * support
 * @tparam ra_t type of the random access used in RLPBWT
 */
template<typename ra_t>
class phi_support {
private:
    /**
     * @brief default value (used to return null from the two functions)
     */
    unsigned int def{};
public:
    /**
     * @brief panel of sparse bitvectors for phi function
     */
    std::vector<sdsl::sd_vector<>> phi_vec;

    /**
     * @brief panel of sparse bitvectors for phi_inv function
     */
    std::vector<sdsl::sd_vector<>> phi_inv_vec;

    /**
     * @brief panel of rank support for phi function
     */
    std::vector<sdsl::sd_vector<>::rank_1_type> phi_rank;

    /**
     * @brief panel of rank support for phi_inv function
     */
    std::vector<sdsl::sd_vector<>::rank_1_type> phi_inv_rank;

    /**
     * @brief panel of select support for phi function
     */
    std::vector<sdsl::sd_vector<>::select_1_type> phi_select;

    /**
     * @brief panel of select support for phi_inv function
     */
    std::vector<sdsl::sd_vector<>::select_1_type> phi_inv_select;

    /**
     * @brief compressed int vector for prefix samples used by phi function
     */
    std::vector<sdsl::int_vector<>> phi_supp;

    /**
     * @brief compressed int vector for prefix samples used by phi_inv function
     */
    std::vector<sdsl::int_vector<>> phi_inv_supp;

    /**
     * @brief default constructor
     */
    phi_support() = default;

    /**
     * @brief default destructor
     */
    virtual ~phi_support() = default;

    /**
     * @brief constructor of the phi/phi_inv support data structure
     * @param cols vector of the columns of the RLPBWT
     * @param panelbv random access data structure for the panel of RLPBWT
     * @param last_pref last prefix array of the PBWT
     * @param verbose bool for extra prints
     */
    explicit phi_support(std::vector<column_ms> &cols, ra_t *panelbv,
                         sdsl::int_vector<> &last_pref, bool verbose = false) {
        // TODO select are useless at the moment
        // default value is the panel height
        this->def = panelbv->h;
        // initialize temporary panel of no-sparse bitvectors
        auto phi_tmp = std::vector<sdsl::bit_vector>(panelbv->h,
                                                     sdsl::bit_vector(
                                                             panelbv->w + 1,
                                                             0));
        auto phi_inv_tmp = std::vector<sdsl::bit_vector>(panelbv->h,
                                                         sdsl::bit_vector(
                                                                 panelbv->w + 1,
                                                                 0));
        // initialize panels and vectors of the data structure
        this->phi_vec = std::vector<sdsl::sd_vector<>>(panelbv->h);
        this->phi_inv_vec = std::vector<sdsl::sd_vector<>>(panelbv->h);
        this->phi_rank = std::vector<sdsl::sd_vector<>::rank_1_type>(
                panelbv->h);
        this->phi_inv_rank = std::vector<sdsl::sd_vector<>::rank_1_type>(
                panelbv->h);
        this->phi_select = std::vector<sdsl::sd_vector<>::select_1_type>(
                panelbv->h);
        this->phi_inv_select = std::vector<sdsl::sd_vector<>::select_1_type>(
                panelbv->h);
        this->phi_supp = std::vector<sdsl::int_vector<>>(panelbv->h);
        this->phi_inv_supp = std::vector<sdsl::int_vector<>>(panelbv->h);

        // remporary vector for supports
        std::vector<std::vector<unsigned int>> phi_supp_tmp(panelbv->h);
        std::vector<std::vector<unsigned int>> phi_inv_supp_tmp(panelbv->h);

        // iterate over every column
        for (unsigned int i = 0; i < cols.size(); i++) {
            // iterate over every run index
            for (unsigned int j = 0; j < cols[i].sample_beg.size(); j++) {
                // use sample beg to compute phi panel
                phi_tmp[cols[i].sample_beg[j]][i] = true;
                // use sample_end to compute
                // support phi panel (if we are in the first run we use default
                // value)
                if (j == 0) {
                    phi_supp_tmp[cols[i].sample_beg[j]].push_back(panelbv->h);
                } else {
                    phi_supp_tmp[cols[i].sample_beg[j]].push_back(
                            cols[i].sample_end[j - 1]);
                }

                // use sample end to compute phi_inv panel
                phi_inv_tmp[cols[i].sample_end[j]][i] = true;

                // use sample_beg to compute
                // support phi panel (if we are in the last run we use default
                // value)
                if (j == cols[i].sample_beg.size() - 1) {
                    phi_inv_supp_tmp[cols[i].sample_end[j]].push_back(
                            panelbv->h);
                } else {
                    phi_inv_supp_tmp[cols[i].sample_end[j]].push_back(
                            cols[i].sample_beg[j + 1]);
                }
            }
        }

        // use the last prefix array to compute the remain values for the
        // phi support data structure (with the same "rules" of the previous
        // case)
        for (unsigned int j = 0; j < phi_supp_tmp.size(); j++) {
            if (!phi_tmp[j][phi_tmp[j].size() - 1]) {
                phi_tmp[j][phi_tmp[j].size() - 1] = true;
            }
            if (j == 0) {
                if (phi_supp_tmp[last_pref[j]].empty() ||
                    phi_supp_tmp[last_pref[j]].back() != panelbv->h) {
                    phi_supp_tmp[last_pref[j]].push_back(panelbv->h);
                }
            } else {
                if (phi_supp_tmp[last_pref[j]].empty() ||
                    phi_supp_tmp[last_pref[j]].back() != last_pref[j - 1]) {
                    phi_supp_tmp[last_pref[j]].push_back(last_pref[j - 1]);
                }
            }
            if (!phi_inv_tmp[j][phi_inv_tmp[j].size() - 1]) {
                phi_inv_tmp[j][phi_inv_tmp[j].size() - 1] = true;
            }
            if (j == phi_supp_tmp.size() - 1) {
                if (phi_inv_supp_tmp[last_pref[j]].empty() ||
                    phi_inv_supp_tmp[last_pref[j]].back() != panelbv->h) {
                    phi_inv_supp_tmp[last_pref[j]].push_back(panelbv->h);
                }
            } else {
                if (phi_inv_supp_tmp[last_pref[j]].empty() ||
                    phi_inv_supp_tmp[last_pref[j]].back() != last_pref[j + 1]) {
                    phi_inv_supp_tmp[last_pref[j]].push_back(last_pref[j + 1]);
                }
            }
        }
        // compress and push the support sdsl int_vectors
        for (unsigned int i = 0; i < phi_supp_tmp.size(); i++) {
            sdsl::int_vector<> tmp(phi_supp_tmp[i].size());
            for (unsigned int j = 0; j < phi_supp_tmp[i].size(); j++) {
                tmp[j] = phi_supp_tmp[i][j];
            }
            sdsl::util::bit_compress(tmp);
            this->phi_supp[i] = tmp;

            sdsl::int_vector<> tmp_inv(phi_inv_supp_tmp[i].size());
            for (unsigned int j = 0; j < phi_inv_supp_tmp[i].size(); j++) {
                tmp_inv[j] = phi_inv_supp_tmp[i][j];
            }
            sdsl::util::bit_compress(tmp_inv);
            this->phi_inv_supp[i] = tmp_inv;
        }
        // create sparse bit vector and relative rank/select for phi panel
        unsigned int count = 0;
        for (auto &i: phi_tmp) {
            this->phi_vec[count] = sdsl::sd_vector<>(i);
            this->phi_rank[count] = sdsl::sd_vector<>::rank_1_type(
                    &this->phi_vec[count]);
            this->phi_select[count] = sdsl::sd_vector<>::select_1_type(
                    &this->phi_vec[count]);
            count++;
        }
        // create sparse bit vector and relative rank/select for phi_inv panel
        count = 0;
        for (auto &i: phi_inv_tmp) {
            this->phi_inv_vec[count] = sdsl::sd_vector<>(i);
            this->phi_inv_rank[count] = sdsl::sd_vector<>::rank_1_type(
                    &this->phi_inv_vec[count]);
            this->phi_inv_select[count] = sdsl::sd_vector<>::select_1_type(
                    &this->phi_inv_vec[count]);
            count++;
        }


        if (verbose) {
            unsigned int index = 0;
            for (auto &col: cols) {
                std::cout << index << ":\t" << col.sample_beg
                          << "\n";
                std::cout << index << ":\t" << col.sample_end
                          << "\n";
                std::cout << "\n";
                index++;
            }
            std::cout << "\n\n";
            index = 0;
            for (auto &i: this->phi_vec) {
                std::cout << index << ":\t" << i << "\n";
                index++;
            }
            std::cout << "----------\n";
            index = 0;
            for (auto &i: this->phi_inv_vec) {
                std::cout << index << ":\t" << i << "\n";
                index++;
            }
            std::cout << "\n\n";
            index = 0;
            for (auto &i: this->phi_supp) {
                std::cout << index << ":\t" << i << "\n";
                index++;
            }
            std::cout << "----------\n";
            index = 0;
            for (auto &i: this->phi_inv_supp) {
                std::cout << index << ":\t" << i << "\n";
                index++;
            }
        }
    }

    /**
     * @brief phi function that return an optional
     * @param pref prefix array value
     * @param col current column index
     * @return previous prefix array value at current column (if exists)
     */
    std::optional<unsigned int> phi(unsigned int pref, unsigned int col) {
        auto res = static_cast<unsigned int>(this->phi_supp[pref][this->phi_rank[pref](
                col)]);
        if (res == this->def) {
            return std::nullopt;
        } else {
            return res;
        }
    }

    /**
     * @brief phi_inv function that return an optional
     * @param pref prefix array value
     * @param col current column index
     * @return next prefix array value at current column (if exists)
    */
    std::optional<unsigned int> phi_inv(unsigned int pref, unsigned int col) {
        auto res = static_cast<unsigned int>(this->phi_inv_supp[pref][this->phi_inv_rank[pref](
                col)]);
        if (res == this->def) {
            return std::nullopt;
        } else {
            return res;
        }
    }

    /**
     * @brief function to obtain size in bytes of the phi/phi_inv support data
     * structure
     * @param verbose bool for extra prints
     * @return size in bytes
    */
    unsigned long long size_in_bytes(bool verbose = false) {
        unsigned long long size = 0;
        for (unsigned int i = 0; i < this->phi_vec.size(); ++i) {
            size += sdsl::size_in_bytes(phi_vec[i]);
            size += sdsl::size_in_bytes(phi_inv_vec[i]);
            size += sdsl::size_in_bytes(phi_rank[i]);
            size += sdsl::size_in_bytes(phi_select[i]);
            size += sdsl::size_in_bytes(phi_inv_rank[i]);
            size += sdsl::size_in_bytes(phi_inv_select[i]);
            size += sdsl::size_in_bytes(phi_supp[i]);
            size += sdsl::size_in_bytes(phi_inv_supp[i]);
        }
        if (verbose) {
            std::cout << "phi: " << size << " bytes\n";
        }
        return size;
    }

    /**
     * @brief function to obtain size in megabytes of the phi/phi_inv support data
     * structure
     * @param verbose bool for extra prints
     * @return size in megabytes
     */
    double size_in_mega_bytes(bool verbose = false) {
        double size = 0;
        for (unsigned int i = 0; i < this->phi_vec.size(); ++i) {
            size += sdsl::size_in_mega_bytes(phi_vec[i]);
            size += sdsl::size_in_mega_bytes(phi_inv_vec[i]);
            size += sdsl::size_in_mega_bytes(phi_rank[i]);
            size += sdsl::size_in_mega_bytes(phi_select[i]);
            size += sdsl::size_in_mega_bytes(phi_inv_rank[i]);
            size += sdsl::size_in_mega_bytes(phi_inv_select[i]);
            size += sdsl::size_in_mega_bytes(phi_supp[i]);
            size += sdsl::size_in_mega_bytes(phi_inv_supp[i]);
        }
        if (verbose) {
            std::cout << "phi: " << size << " megabytes\n";
        }
        return size;
    }

    /**
     * @brief function to serialize the phi/phi_inv data structure object
     * @param out std::ostream object to stream the serialization
     * @return size of the serialization
     */
    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr,
                     const std::string &name = "") {
        sdsl::structure_tree_node *child =
                sdsl::structure_tree::add_child(v, name,
                                                sdsl::util::class_name(
                                                        *this));
        size_t written_bytes = 0;
        out.write((char *) &this->def, sizeof(this->def));
        written_bytes += sizeof(this->def);
        for (unsigned int i = 0; i < this->phi_vec.size(); i++) {
            std::string label = "phi_vec_" + std::to_string(i);
            written_bytes += this->phi_vec[i].serialize(out, child, label);
        }
        for (unsigned int i = 0; i < this->phi_inv_vec.size(); i++) {
            std::string label = "phi_inv_vec_" + std::to_string(i);
            written_bytes += this->phi_inv_vec[i].serialize(out, child, label);
        }
        for (unsigned int i = 0; i < this->phi_supp.size(); i++) {
            std::string label = "phi_supp_" + std::to_string(i);
            written_bytes += this->phi_supp[i].serialize(out, child, label);
        }
        for (unsigned int i = 0; i < this->phi_inv_supp.size(); i++) {
            std::string label = "phi_inv_supp_" + std::to_string(i);
            written_bytes += this->phi_inv_supp[i].serialize(out, child,
                                                             label);
        }

        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    /**
     * @brief function to load the phi/phi_inv data structure object
     * @param in std::istream object from which load the phi/phi_inv data
     * structure object
     */
    void load(std::istream &in) {
        in.read((char *) &this->def, sizeof(this->def));

        for (unsigned int i = 0; i < this->def; i++) {
            auto s = new sdsl::sd_vector<>();
            s->load(in);
            this->phi_vec.emplace_back(*s);
        }
        for (unsigned int i = 0; i < this->def; i++) {
            auto s = new sdsl::sd_vector<>();
            s->load(in);
            this->phi_inv_vec.emplace_back(*s);
        }
        for (unsigned int i = 0; i < this->def; i++) {
            auto s = new sdsl::int_vector<>();
            s->load(in);
            this->phi_supp.emplace_back(*s);
        }
        for (unsigned int i = 0; i < this->def; i++) {
            auto s = new sdsl::int_vector<>();
            s->load(in);
            this->phi_inv_supp.emplace_back(*s);
        }

        for (unsigned int i = 0; i < this->phi_vec.size(); i++) {
            this->phi_rank.emplace_back(sdsl::sd_vector<>::rank_1_type(
                    &this->phi_vec[i]));
            this->phi_select.emplace_back(sdsl::sd_vector<>::select_1_type(
                    &this->phi_vec[i]));
        }
        for (unsigned int i = 0; i < this->phi_inv_vec.size(); i++) {
            this->phi_inv_rank.emplace_back(sdsl::sd_vector<>::rank_1_type(
                    &this->phi_inv_vec[i]));
            this->phi_inv_select.emplace_back(sdsl::sd_vector<>::select_1_type(
                    &this->phi_inv_vec[i]));
        }
    }
};


#endif //RLPBWT_PHI_SUPPORT_H
