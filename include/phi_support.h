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
        this->def = panelbv->h;
        auto phi_tmp = std::vector<sdsl::bit_vector>(panelbv->h,
                                                     sdsl::bit_vector(
                                                             panelbv->w,
                                                             0));
        auto phi_inv_tmp = std::vector<sdsl::bit_vector>(panelbv->h,
                                                         sdsl::bit_vector(
                                                                 panelbv->w,
                                                                 0));
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
        this->phi_supp = std::vector(panelbv->h,
                                        sdsl::int_vector(panelbv->w));
        this->phi_inv_supp = std::vector(panelbv->h,
                                         sdsl::int_vector(
                                                    panelbv->w));

        std::vector<std::pair<unsigned int, unsigned int>> counts(
                panelbv->h, std::make_pair(0, 0));

        for (unsigned int i = 0; i < cols.size(); i++) {
            for (unsigned int j = 0;
                 j < cols[i].sample_beg.size(); j++) {
                phi_tmp[cols[i].sample_beg[j]][i] = true;
                if (j == 0) {
                    this->phi_supp[cols[i].sample_beg[j]]
                    [counts[cols[i].sample_beg[j]].first] =
                            panelbv->h;
                } else {
                    this->phi_supp[cols[i].sample_beg[j]]
                    [counts[cols[i].sample_beg[j]].first] =
                            cols[i].sample_end[j - 1];
                }
                counts[cols[i].sample_beg[j]].first++;
                phi_inv_tmp[cols[i].sample_end[j]][i] = true;
                if (j == cols[i].sample_beg.size() - 1) {
                    this->phi_inv_supp[cols[i].sample_end[j]]
                    [counts[cols[i].sample_end[j]].second] =
                            panelbv->h;
                } else {
                    this->phi_inv_supp[cols[i].sample_end[j]]
                    [counts[cols[i].sample_end[j]].second] =
                            cols[i].sample_beg[j + 1];
                }
                counts[cols[i].sample_end[j]].second++;
            }
        }
        for (unsigned int j = 0; j < counts.size(); j++) {
            if (!phi_tmp[j][phi_tmp[j].size() - 1]) {
                phi_tmp[j][phi_tmp[j].size() - 1] = true;
                if (j == 0) {
                    this->phi_supp[last_pref[j]]
                    [counts[last_pref[j]].first] = panelbv->h;
                } else {
                    this->phi_supp[last_pref[j]]
                    [counts[last_pref[j]].first] =
                            last_pref[j - 1];
                }
                counts[last_pref[j]].first++;
            }

            if (!phi_inv_tmp[j][phi_tmp[j].size() - 1]) {
                phi_inv_tmp[j][phi_inv_tmp[j].size() - 1] = true;
                if (j == counts.size() - 1) {
                    this->phi_inv_supp[last_pref[j]]
                    [counts[last_pref[j]].second] = panelbv->h;
                } else {
                    this->phi_inv_supp[last_pref[j]]
                    [counts[last_pref[j]].second] =
                            last_pref[j + 1];
                }
                counts[last_pref[j]].second++;
            }
        }
        for (unsigned int i = 0; i < counts.size(); i++) {
            this->phi_supp[i].resize(counts[i].first);
            sdsl::util::bit_compress(this->phi_supp[i]);
            this->phi_inv_supp[i].resize(counts[i].second);
            sdsl::util::bit_compress(this->phi_inv_supp[i]);
        }
        unsigned int count = 0;
        for (auto &i: phi_tmp) {
            this->phi_vec[count] = sdsl::sd_vector<>(i);
            this->phi_rank[count] = sdsl::sd_vector<>::rank_1_type(
                    &this->phi_vec[count]);
            this->phi_select[count] = sdsl::sd_vector<>::select_1_type(
                    &this->phi_vec[count]);
            count++;
        }
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
        if(verbose){
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
        if(verbose){
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
