//
// Created by dlcgold on 10/02/22.
//

#ifndef RLPBWT_RLPBWT_MS_H
#define RLPBWT_RLPBWT_MS_H

#include <vector>
#include <type_traits>
#include "column_ms.h"
#include "exceptions.h"
#include "panel_ra.h"
#include "slp_panel_ra.h"
#include "phi_support.h"
#include "ms.h"
#include "ms_matches.h"

/**
 * @brief data structure for matching-statistics supported RLPBWT
 * @tparam ra_t type of panel (panel_ra or slp_panel_ra)
 */
template<typename ra_t>
class rlpbwt_ms {
public:
    /**
     * @brief bool to check if the RLPBWT support the use of thresholds
     */
    bool is_thr_enabled{};

    /**
     * @brief bool to check if the RLPBWT is extended with phi/phi_inv structure
     */
    bool is_extended{};

    /**
    * @brief vector of matcing statistics supported columns
    */
    std::vector<column_ms> cols;

    /**
     * @brief panel (panel_ra or slp_panel_ra)
     */

    ra_t *panel;
    /**
     * @brief phi/phi_inv support data structure
     * @tparam type of panel (panel_ra or slp_panel_ra)
     */
    phi_support<ra_t> *phi;

private:

    /**
     * @brief compressed int vector for last prefix array
     */
    sdsl::int_vector<> last_pref;

    /**
     * @brief function to build a naive column for a column of the panel
     * @param column column of the panel
     * @param pref current prefix array
     * @param div current divergence array (as lcp array)
     * @return the new naive column
     */
    column_ms
    build_column(std::string &column, std::vector<unsigned int> &pref,
                 sdsl::int_vector<> &div) {
        unsigned int height_tmp = pref.size();
        // variable for "c" value
        unsigned int count0 = 0;
        // variables for compute "u" and "v" values
        unsigned int u = 0;
        unsigned int v = 0;
        // temporary variables for compute "u" and "v" values
        unsigned int count0tmp = 0;
        unsigned int count1 = 0;
        // variable to compute curr lcs in order to eventually compute thresholds
        // unsigned int lcs = 0;
        // bool to check first symbol fo the column
        bool start = true;

        std::vector<std::pair<unsigned int, unsigned int>> rows;
        unsigned int tmp_beg = 0;
        unsigned int thr_tmp = 0;
        unsigned int lcs = 0;
        // update start and "c" value
        for (unsigned int i = 0; i < height_tmp; i++) {
            if (i == 0 && column[pref[i]] == '1') {
                start = false;
            }
            if (column[i] == '0') {
                count0++;
            }
        }

        // initialize the three bitvectors
        sdsl::bit_vector runvec(height_tmp + 1, 0);
        sdsl::bit_vector thr(height_tmp, 0);
        sdsl::bit_vector zerovec(count0, 0);
        sdsl::bit_vector onevec(height_tmp - count0, 0);

        // bools to check if we are at the beginning of a run and if we have to swap
        // the counting of zeros and ones
        bool begrun = true;
        bool pushz = false;
        bool pusho = false;

        // first swap according to first value of the column
        if (start) {
            pusho = true;
        } else {
            pushz = true;
        }

        // iteration over the entire column
        for (unsigned int i = 0; i < height_tmp; i++) {
            // if we are at the beginning of a run we save previous temporary values
            // for "u" and "v"
            if (begrun) {
                tmp_beg = pref[i];
                u = count0tmp;
                v = count1;
                begrun = false;
            }

            // increment temporary variables
            if (column[pref[i]] == '1') {
                count1++;
            } else {
                count0tmp++;
            }


            // stuff for thresholds
            if (this->is_thr_enabled) {
                if ((i == 0) || (column[pref[i]] != column[pref[i - 1]])) {
                    thr_tmp = i;
                    lcs = div[i];
                }
                if (div[i] < lcs) {
                    thr_tmp = i;
                    lcs = div[i];
                }
            }

            if ((i == height_tmp - 1) ||
                (column[pref[i]] != column[pref[i + 1]])) {
                // 1 in bitvectors for runs et every end of a run
                runvec[i] = true;
                // 1 in bitvectors for thresholds index
                if (this->is_thr_enabled) {
                    thr[thr_tmp] = true;
                }
                rows.emplace_back(tmp_beg, pref[i]);
                // update bitvectors for "u" and "v" and swap the case to study in
                // next run
                if (pusho) {
                    if (v != 0) {
                        onevec[v - 1] = true;
                    }
                    std::swap(pusho, pushz);
                } else {
                    if (u != 0) {
                        zerovec[u - 1] = true;
                    }
                    std::swap(pusho, pushz);
                }
                begrun = true;
            }
        }

        // set last bit to one in bitvectors "u" and "v"
        if (!zerovec.empty()) {
            zerovec[zerovec.size() - 1] = true;
        }
        if (!onevec.empty()) {
            onevec[onevec.size() - 1] = true;
        }

        // compress div array
        sdsl::util::bit_compress(div);
        sdsl::int_vector<> sample_beg(rows.size());
        sdsl::int_vector<> sample_end(rows.size());
        for (unsigned int i = 0; i < rows.size(); i++) {
            sample_beg[i] = rows[i].first;
            sample_end[i] = rows[i].second;
        }
        sdsl::util::bit_compress(sample_beg);
        sdsl::util::bit_compress(sample_end);
        if (!this->is_thr_enabled) {
            thr = sdsl::bit_vector(0);
        }
        return {start, count0, runvec, zerovec, onevec, thr, sample_beg,
                sample_end};
    }

    /**
     * @brief function to update prefix/divergence (lcp) arrays as in Durbin
     * @param column column of the panel
     * @param pref current prefix array
     * @param div current divergence array (as lcp array)
     */
    static void
    update(std::string &column, std::vector<unsigned int> &pref,
           sdsl::int_vector<> &div) {
        unsigned int height = pref.size();
        std::vector<unsigned int> new_pref(height);
        sdsl::int_vector<> new_div(height);
        unsigned int count0 = 0;
        unsigned int lcs = -1;

        for (unsigned int i = 0; i < height; i++) {
            lcs = std::min(lcs, static_cast<unsigned int> (div[i]));
            if (column[pref[i]] == '0') {
                new_pref[count0] = pref[i];
                new_div[count0] = lcs + 1;
                count0++;
                lcs = INT_MAX;
            }
        }

        int count1 = 0;
        lcs = -1;
        for (unsigned int i = 0; i < height; i++) {
            lcs = std::min(lcs, static_cast<unsigned int> (div[i]));
            if (column[pref[i]] == '1') {
                new_pref[count0 + count1] = pref[i];
                new_div[count0 + count1] = lcs + 1;
                count1++;
                lcs = INT_MAX;
            }
        }
        new_div[0] = 0;
        if (count0 != height) {
            new_div[count0] = 0;
        }
        div = new_div;
        pref = new_pref;
    }

    /**
     * @brief function to compute the lf mapping, w(i, s) function in Durbin
     * @param col_index index of the column
     * @param index virtual index of the row of the original panel
     * @param symbol symbol s
     * @param verbose bool for extra print
     * @return the index computed with the lf-mapping
     */
    unsigned int
    lf(unsigned int col_index, unsigned int index, char symbol,
       bool verbose = false) const {
        // obtain "u" and "v"
        auto uv = uvtrick(col_index, index);
        if (verbose) {
            std::cout << "at col: " << col_index << " and index " << index
                      << " with symbol " << symbol << "\n";
            std::cout << uv.first << " " << uv.second << " with c: "
                      << this->cols[col_index].count_0 << "\n";
        }
        // computer lf-mapping as Durbin's w(i, s)
        if (symbol == '0') {
            return uv.first;
        } else {
            return this->cols[col_index].count_0 + uv.second;
        }
    }

    /**
     * @brief trick to extract u and v value from a run in rlpbwt column
     * @param col_index index of the column
     * @param index virtual index of the row of the original panel
     * @return a std::pair with u as first and v as second
     */
    std::pair<unsigned int, unsigned int>
    uvtrick(unsigned int col_index, unsigned int index) const {
        // for index 0 u = v = 0
        if (index == 0) {
            return {0, 0};
        }

        // number of previous zeros and ones
        unsigned int prev_u = 0;
        unsigned int prev_v = 0;
        // runs of previous zeros and ones in respective bitvectors
        unsigned int run_prev_u = 0;
        unsigned int run_prev_v = 0;
        // offset to finally consider inside the run
        unsigned int offset = 0;

        // compute the run of the index using rank over runs bitvector
        unsigned run = this->cols[col_index].rank_runs(index);
        // two main cases:
        // - run > 1
        // - run = 0 or run = 1 (we separate this two cases due to some errors using
        //   select)
        if (run > 1) {
            // if the run is an even index run it means that before we have the same
            // amount of zeros runs and ones runs
            if (run % 2 == 0) {
                // we compute the number of the zeros and ones runs before the
                // current run
                run_prev_u = run / 2;
                run_prev_v = run / 2;

                // we compute the number of zeros and ones before the current run
                // using select over the two bitvectors of zeros and ones
                prev_u = this->cols[col_index].select_u(run_prev_u) + 1;
                prev_v = this->cols[col_index].select_v(run_prev_v) + 1;

                // we compute the offset
                offset = index - (this->cols[col_index].select_runs(run) + 1);

                // depending on the first symbol of the column the offset increase
                // the number of zeros or the number of ones (they increase
                // alternating)
                if (this->cols[col_index].zero_first) {
                    return {prev_u + offset, prev_v};
                } else {
                    return {prev_u, prev_v + offset};
                }
            } else {
                // if the run is an odd index run it means that before we have x
                // runs of zeros and x+1 runs of ones, or vice-versa,
                // depending on the first symbol of the column
                if (this->cols[col_index].zero_first) {
                    run_prev_u = (run / 2) + 1;
                    run_prev_v = run / 2;
                } else {
                    run_prev_u = run / 2;
                    run_prev_v = (run / 2) + 1;
                }
                // we compute the number of zeros and ones before the current run
                // using select over the two bitvectors of zeros and ones
                prev_u = this->cols[col_index].select_u(run_prev_u) + 1;
                prev_v = this->cols[col_index].select_v(run_prev_v) + 1;

                // we compute the offset
                offset = index - (this->cols[col_index].select_runs(run) + 1);

                // depending on the first symbol of the column the offset increase
                // the number of zeros or the number of ones (they increase
                // alternating)
                if (this->cols[col_index].zero_first) {
                    return {prev_u, prev_v + offset};
                } else {
                    return {prev_u + offset, prev_v};
                }
            }
        } else {
            if (run == 0) {
                // if we are in the first run (of index 0) simply v is always 0 if
                // we are in a zero beginning column and vice-versa for u if we are
                // in a one beginning column. The non-zero value is the same of
                // index because it's the first run
                if (this->cols[col_index].zero_first) {
                    return {index, 0};
                } else {
                    return {0, index};
                }
            } else {
                // if we are in the second run (of index 1) we fix one of the two
                // values (based practically on the length of the first run),
                // according to the first symbol of the column, and we
                // increase the other one according to the index
                if (this->cols[col_index].zero_first) {
                    return {this->cols[col_index].select_runs(run) + 1,
                            index -
                            (this->cols[col_index].select_runs(run) + 1)};
                } else {
                    return {index -
                            (this->cols[col_index].select_runs(run) + 1),
                            this->cols[col_index].select_runs(run) + 1};
                }
            }
        }
    }

    /**
     * @brief function to extend matching statistics matches with matching
     * rows indices
     * @param ms_matches matching statistics matches that will be extended
     */
    void extend_haplos(ms_matches &ms_matches) {
        for (unsigned int i = 0; i < ms_matches.basic_matches.size(); i++) {
            std::vector<unsigned int> haplos;
            unsigned int start_row = std::get<0>(
                    ms_matches.basic_matches[i]);
            haplos.emplace_back(start_row);
            unsigned int curr_len = std::get<1>(
                    ms_matches.basic_matches[i]);
            unsigned int curr_col = std::get<2>(
                    ms_matches.basic_matches[i]);
            bool check_down = true;
            unsigned int down_row = 0;
            bool check_up = true;
            unsigned int up_row = 0;
            //down
            if constexpr (std::is_same_v<ra_t, slp_panel_ra>) {
                while (check_down) {
                    auto phi_res = this->phi->phi_inv(start_row, curr_col);
                    if (!phi_res.has_value()) {
                        break;
                    }
                    down_row = phi_res.value();
                    if (lceBound(curr_col, start_row, down_row,
                                 curr_len)) {
                        haplos.emplace_back(down_row);
                        start_row = down_row;
                    } else {
                        check_down = false;
                    }
                }
                //up
                while (check_up) {
                    auto phi_res = this->phi->phi(start_row, curr_col);
                    if (!phi_res.has_value()) {
                        break;
                    }
                    up_row = phi_res.value();

                    if (lceBound(curr_col, start_row, up_row, curr_len)) {
                        haplos.emplace_back(up_row);
                        start_row = up_row;
                    } else {
                        check_up = false;
                    }
                }

            } else {
                // TODO use panelbv
                while (check_down) {
                    auto phi_res = this->phi->phi_inv(start_row, curr_col);
                    if (!phi_res.has_value()) {
                        break;
                    }
                    down_row = phi_res.value();
                    if (this->panel->lceToRCheck(curr_col, start_row, down_row,
                                                 curr_len)) {
                        haplos.emplace_back(down_row);
                        start_row = down_row;
                    } else {
                        check_down = false;
                    }
                }
                //up
                while (check_up) {
                    auto phi_res = this->phi->phi(start_row, curr_col);
                    if (!phi_res.has_value()) {
                        break;
                    }
                    up_row = phi_res.value();
                    if (this->panel->lceToRCheck(curr_col, start_row, up_row,
                                                 curr_len)) {
                        haplos.emplace_back(up_row);
                        start_row = up_row;
                    } else {
                        check_up = false;
                    }
                }
            }

            ms_matches.haplos.emplace_back(haplos);
        }
    }

public:
    /**
    * @brief height of the panel
    */
    unsigned int width{};

    /**
     * @brief width of the panel
     */
    unsigned int height{};

    /**
     * @brief default constructor
     */
    rlpbwt_ms() = default;

    /**
    * @brief default destructor
    */
    virtual ~rlpbwt_ms() = default;

    /**
     * @brief costructor of a RLPBWT that support matching statistics
     * @param filename file with the panel
     * @param thr bool to enable thresholds computation
     * @param verbose bool fro extra prints
     * @param slp_filename file with the slp of the panel
     */
    explicit rlpbwt_ms(const char *filename, bool thr = false,
                       bool verbose = false, const char *slp_filename = "") {
        std::ifstream input_matrix(filename);
        if constexpr (!std::is_same_v<ra_t, panel_ra>) {
            if (std::string(slp_filename).empty()) {
                throw SlpNotFoundException{};
            }
        }
        if (input_matrix.is_open()) {
            this->is_thr_enabled = thr;
            std::string header1;
            std::string header2;
            std::string line;
            std::string garbage;
            std::string new_column;

            getline(input_matrix, header1);
            getline(input_matrix, header2);
            getline(input_matrix, line);
            std::istringstream is(line);
            is >> garbage >> garbage >> garbage >> garbage >> new_column;
            unsigned int tmp_height = new_column.size();
            std::cout << "h: " << tmp_height << "\n";
            unsigned int tmp_width = std::count(
                    std::istreambuf_iterator<char>(input_matrix),
                    std::istreambuf_iterator<char>(), '\n') + 1;
            std::cout << "w: " << tmp_width << "\n";
            input_matrix.clear();
            input_matrix.seekg(0, std::ios::beg);
            this->cols = std::vector<column_ms>(tmp_width);
            std::vector<unsigned int> pref(tmp_height);
            std::vector<unsigned int> lpref(tmp_height);
            sdsl::int_vector<> div(tmp_height);
            this->last_pref.resize(tmp_height);
            for (unsigned int i = 0; i < tmp_height; i++) {
                pref[i] = i;
                div[i] = 0;
            }
            if constexpr (std::is_same_v<ra_t, panel_ra>) {
                this->panel = new ra_t(tmp_height, tmp_width);
            } else if constexpr (std::is_same_v<ra_t, slp_panel_ra>) {
                this->panel = new ra_t(slp_filename, tmp_height, tmp_width);
            } else {
                throw WrongRaTypeException{};
            }
            unsigned int count = 0;
            std::string last_col;
            getline(input_matrix, line);
            getline(input_matrix, line);

            while (getline(input_matrix, line) && !line.empty()) {
                std::cout << count << "\r";
                std::istringstream is_col(line);
                is_col >> garbage;
                if (garbage == "TOTAL_SAMPLES:") {
                    break;
                }
                is_col >> garbage >> garbage >> garbage >> new_column;
                if (verbose) {
                    std::cout << "\nnew_column " << count << "\n";
                    std::cout << new_column << "\n" << this->cols[count].runs
                              << "\n"
                              << this->cols[count].u << "\n"
                              << this->cols[count].v
                              << "\n-------------------------------\n";
                }
                auto col = rlpbwt_ms::build_column(new_column, pref, div);
                if constexpr (std::is_same_v<ra_t, panel_ra>) {
                    for (unsigned int k = 0; k < new_column.size(); k++) {
                        if (new_column[k] != '0') {
                            this->panel->panel[count][k] = true;
                        }
                    }
                }
                this->cols[count] = col;
                this->cols[count].rank_runs = sdsl::sd_vector<>::rank_1_type(
                        &this->cols[count].runs);
                this->cols[count].select_runs = sdsl::sd_vector<>::select_1_type(
                        &this->cols[count].runs);
                this->cols[count].rank_u = sdsl::sd_vector<>::rank_1_type(
                        &this->cols[count].u);
                this->cols[count].select_u = sdsl::sd_vector<>::select_1_type(
                        &this->cols[count].u);
                this->cols[count].rank_v = sdsl::sd_vector<>::rank_1_type(
                        &this->cols[count].v);
                this->cols[count].select_v = sdsl::sd_vector<>::select_1_type(
                        &this->cols[count].v);
                this->cols[count].rank_thr = sdsl::sd_vector<>::rank_1_type(
                        &this->cols[count].thr);
                this->cols[count].select_thr = sdsl::sd_vector<>::select_1_type(
                        &this->cols[count].thr);
                if (count == tmp_width - 1) {
                    lpref = pref;
                }
                rlpbwt_ms::update(new_column, pref, div);
                count++;
            }
            for (unsigned int i = 0; i < lpref.size(); i++) {
                this->last_pref[i] = lpref[i];
            }
            sdsl::util::bit_compress(this->last_pref);
            this->is_extended = false;
            this->width = tmp_width;
            this->height = tmp_height;
            input_matrix.close();
        } else {
            throw FileNotFoundException{};
        }
    }

    /**
     * @brief function to extend the RLPBWT with the phi/phi_inv structure
     * @param verbose bool for extra prints
     */
    void extend(bool verbose = false) {
        if (!this->is_extended) {
            this->phi = new phi_support<ra_t>(this->cols, this->panel,
                                              this->last_pref, verbose);
            this->is_extended = true;
        }
    }

    /**
     * @brief function to delete the phi/phi_inv structure
     */
    void unextend() {
        if (this->is_extended) {
            this->phi = nullptr;
            this->is_extended = false;
        }
    }

    /**
     * @brief function to check if longest common extension between two rows is
     * equal or greater than a bound ending at a given column
     * @param col ending column column
     * @param curr current row
     * @param other other row
     * @param bound boiund to check
     * @param verbose bool for extra prints
     * @return true if ongest common extension between two rows is
     * equal or greater than the bound
     * @attention this function is enabled iff the panel is an SLP
     */
    template<typename U = ra_t>
    std::enable_if_t<sizeof(U) && (!std::is_same<ra_t, panel_ra>::value), bool>
    lceBound(unsigned int col, unsigned int curr, unsigned int other,
             unsigned int bound, bool verbose = false) {

        if (col == 0) {
            return false;
        }
        unsigned int rev_col = (this->panel->w - 1) - col + 1;

        unsigned int pos_curr = rev_col + ((this->panel->w) * curr);
        unsigned int pos_other = rev_col + ((this->panel->w) * other);
        if (verbose) {
            std::cout << "at " << rev_col << ": " << pos_curr << ", "
                      << pos_other << "\n";
        }
//        unsigned int lcp_pair_naive = lceToR_NaiveBounded(this->panel->panel, pos_other,
//                                              pos_curr, bound);
        // TODO check the bound in no naive version, it doesn't work
        unsigned int lcp_pair = lceToRBounded(this->panel->panel, pos_other,
                                              pos_curr, bound);
        // TODO if the bound wotks it should be lcp_pair == bound
        if (lcp_pair >= bound) {
            return true;
        }
        return false;
    }

    /**
     * @brief cunction to compute longest common extension between two rows
     * ending at a given column
     * @param col ending column
     * @param curr current row
     * @param other other row
     * @param verbose
     * @return the lce value and the other row index
     * @attention this function is enabled iff the panel is an SLP
     */
    template<typename U = ra_t>
    std::enable_if_t<sizeof(U) && (!std::is_same<ra_t, panel_ra>::value),
            std::pair<unsigned int, unsigned int>>
    lce(unsigned int col, unsigned int curr, unsigned int other,
        bool verbose = false) {

        if (col == 0) {
            return std::make_pair(other, 0);
        }
        unsigned int rev_col = (this->panel->w - 1) - col + 1;

        unsigned int pos_curr = rev_col + ((this->panel->w) * curr);
        unsigned int pos_other = rev_col + ((this->panel->w) * other);
        if (verbose) {
            std::cout << "at " << rev_col << ": " << pos_curr << ", "
                      << pos_other << "\n";
        }
        unsigned int lcp_other = lceToR(this->panel->panel, pos_other,
                                        pos_curr);
        if (lcp_other >= col) {
            lcp_other = col;
        }
        return std::make_pair(other, lcp_other);
    }

    /**
     * @brief function to compute longest common extension between a row and
     * other two rows ending at a given column
     * @param col ending column
     * @param curr current row
     * @param prev first row to compare
     * @param next second row to compare
     * @param verbose bool
     * @return the best lce value and the best row index between curr row and
     * the other two
     * @attention this function is enabled iff the panel is an SLP
     */
    template<typename U = ra_t>
    std::enable_if_t<sizeof(U) && (!std::is_same<ra_t, panel_ra>::value),
            std::pair<unsigned int, unsigned int>>
    lce_pair(unsigned int col, unsigned int curr, unsigned int prev,
             unsigned int next, bool verbose = false) {

        if (col == 0) {
            return std::make_pair(next, 0);
        }
        unsigned int rev_col = (this->panel->w - 1) - col + 1;

        unsigned int pos_curr = rev_col + ((this->panel->w) * curr);
        unsigned int pos_prev = rev_col + ((this->panel->w) * prev);
        unsigned int pos_next = rev_col + ((this->panel->w) * next);
        if (verbose) {
            std::cout << "at " << rev_col << ": " << pos_curr << ", "
                      << pos_prev
                      << ", " << pos_next << "\n";
        }
        unsigned int lcp_prev = lceToR(this->panel->panel, pos_prev,
                                       pos_curr);
        if (lcp_prev >= col) {
            lcp_prev = col;
        }
        unsigned int lcp_next = lceToR(this->panel->panel, pos_curr,
                                       pos_next);
        if (lcp_next >= col) {
            lcp_next = col;
        }
        if (lcp_prev > lcp_next) {
            return std::make_pair(prev, lcp_prev);
        } else {
            return std::make_pair(next, lcp_next);
        }
    }

    /**
     * @brief function to compute matching statistics matches with a given query
     * usning lce queries
     * @param query haplotype query as std::string
     * @param extend_matches bool to check if extend matching statistics matches
     * with rows
     * @param verbose bool for extra prints
     * @return matching statistics matches
     * @attention this function is enabled iff the panel is an SLP
     */
    template<typename U = ra_t>
    std::enable_if_t<sizeof(U) && (!std::is_same<ra_t, panel_ra>::value),
            ms_matches>
    match_lce(const std::string &query, bool extend_matches = false,
              bool verbose = false) {
        if (query.size() != this->panel->w) {
            throw NotEqualLengthException{};
        }
        if (extend_matches && !this->is_extended) {
            this->extend();
        }

        ms ms(query.size());
        auto curr_pos = static_cast<unsigned int>(this->cols[0].sample_end[
                this->cols[0].sample_end.size() - 1]);
        auto curr_index = curr_pos;
        unsigned int curr_run = this->cols[0].rank_runs(curr_index);
        char symbol = get_next_char(this->cols[0].zero_first, curr_run);

        for (unsigned int i = 0; i < query.size(); i++) {
            std::cout << "processed " << i << "\r";
            if (verbose) {
                std::cout << "at " << i << ": " << curr_run << " "
                          << this->cols[i].rank_thr(curr_index) << "\n";
                std::cout << curr_index << " " << curr_run << " " << curr_pos
                          << " "
                          << symbol << "\n";
            }
            if (query[i] == symbol) {
                if (verbose) {
                    std::cout << "match:\n";
                }
                // save in matching statistics row vector
                ms.row[i] = curr_pos;
                if (i != 0) {
                    ms.len[i] = ms.len[i - 1] + 1;
                } else {
                    ms.len[i] = 1;
                }
                if (i != query.size() - 1) {
                    // update index, run, symbol
                    // TODO maybe make a function for this update
                    curr_index = lf(i, curr_index, query[i]);
                    curr_run = this->cols[i + 1].rank_runs(curr_index);
                    symbol = get_next_char(this->cols[i + 1].zero_first,
                                           curr_run);
                    if (verbose) {
                        std::cout << "new: " << curr_index << " " << curr_run
                                  << " "
                                  << curr_pos << " "
                                  << symbol << "\n";
                    }
                }
            } else {
                if (this->cols[i].sample_beg.size() == 1) {
                    if (verbose) {
                        std::cout << "complete mismatch\n";
                    }
                    ms.row[i] = this->panel->h;
                    ms.len[i] = 0;
                    if (i != query.size() - 1) {
                        curr_pos = static_cast<unsigned int>(this->cols[0].sample_end[
                                this->cols[0].sample_end.size() - 1]);
                        curr_index = curr_pos;
                        curr_run = this->cols[i + 1].rank_runs(curr_index);
                        symbol = get_next_char(this->cols[i + 1].zero_first,
                                               curr_run);
                        if (verbose) {
                            std::cout << "update: " << curr_index << " "
                                      << curr_pos
                                      << " "
                                      << symbol << "\n";
                        }
                    }
                } else {
                    if (curr_run == this->cols[i].sample_beg.size() - 1) {
                        curr_index =
                                (this->cols[i].select_runs(curr_run) + 1) - 1;
                        auto prev_pos = static_cast<unsigned int>(this->cols[i].sample_end[
                                curr_run - 1]);
                        if (verbose) {
                            std::cout << "end_run with " << curr_pos << ", "
                                      << prev_pos << "\n";
                        }
                        auto lce_value = lce(i, curr_pos, prev_pos, false);
                        ms.row[i] = prev_pos;
                        if (i == 0) {
                            ms.len[i] = 1;
                        } else {
                            ms.len[i] =
                                    std::min(ms.len[i - 1],
                                             lce_value.second) +
                                    1;
                        }
                        curr_pos = prev_pos;
                        if (verbose) {
                            std::cout << "update: " << curr_index << " "
                                      << curr_pos
                                      << " "
                                      << symbol << "\n";
                        }

                        if (i != query.size() - 1) {
                            curr_index = lf(i, curr_index, query[i]);
                            curr_run = this->cols[i + 1].rank_runs(curr_index);
                            symbol = get_next_char(this->cols[i + 1].zero_first,
                                                   curr_run);
                            if (verbose) {
                                std::cout << "new: " << curr_index << " "
                                          << curr_run
                                          << " "
                                          << curr_pos << " " << symbol << "\n";
                            }
                        }
                    } else if (curr_run == 0) {
                        curr_index = (this->cols[i].select_runs(curr_run + 1) +
                                      1);
                        auto next_pos = static_cast<unsigned int>(this->cols[i].sample_beg[
                                curr_run + 1]);
                        if (verbose) {
                            std::cout << "first_run with " << curr_pos << ", "
                                      << next_pos << "\n";
                        }
                        auto lce_value = lce(i, curr_pos, next_pos, false);
                        ms.row[i] = next_pos;
                        if (i == 0) {
                            ms.len[i] = 1;
                        } else {
                            ms.len[i] =
                                    std::min(ms.len[i - 1],
                                             lce_value.second) +
                                    1;
                        }
                        curr_pos = next_pos;
                        if (verbose) {
                            std::cout << "update: " << curr_index << " "
                                      << curr_pos
                                      << " "
                                      << symbol << "\n";
                        }
                        if (i != query.size() - 1) {
                            curr_index = lf(i, curr_index, query[i]);
                            curr_run = this->cols[i + 1].rank_runs(curr_index);
                            symbol = get_next_char(this->cols[i + 1].zero_first,
                                                   curr_run);
                            if (verbose) {
                                std::cout << "new: " << curr_index << " "
                                          << curr_run
                                          << " "
                                          << curr_pos << " " << symbol << "\n";
                            }
                        }
                    } else {
                        auto prev_pos = static_cast<unsigned int>(this->cols[i].sample_end[
                                curr_run - 1]);
                        auto next_pos = static_cast<unsigned int>(this->cols[i].sample_beg[
                                curr_run + 1]);
                        if (verbose) {
                            std::cout << "curr: " << curr_pos << ", prev: "
                                      << prev_pos << ", next: " << next_pos
                                      << " -> ";
                        }
                        auto lce = this->lce_pair(i, curr_pos, prev_pos,
                                                  next_pos,
                                                  false);
                        if (verbose) {
                            std::cout << lce.first << ", " << lce.second
                                      << "\n";
                        }
                        if (lce.first == next_pos) {
                            curr_pos = next_pos;
                            ms.row[i] = curr_pos;
                            if (i == 0) {
                                ms.len[i] = 1;
                            } else {
                                ms.len[i] =
                                        std::min(ms.len[i - 1], lce.second) +
                                        1;
                            }
                            curr_index = (
                                    this->cols[i].select_runs(curr_run + 1) +
                                    1);
                            if (i != query.size() - 1) {
                                curr_index = lf(i, curr_index, query[i]);
                                curr_run = this->cols[i + 1].rank_runs(
                                        curr_index);
                                symbol = get_next_char(
                                        this->cols[i + 1].zero_first,
                                        curr_run);
                                if (verbose) {
                                    std::cout << "new: " << curr_index << " "
                                              << curr_run
                                              << " "
                                              << curr_pos << " " << symbol
                                              << "\n";
                                }
                            }
                        } else {
                            curr_pos = prev_pos;
                            ms.row[i] = curr_pos;
                            if (i == 0) {
                                ms.len[i] = 1;
                            } else {
                                ms.len[i] =
                                        std::min(ms.len[i - 1], lce.second) +
                                        1;
                            }
                            curr_index =
                                    (this->cols[i].select_runs(curr_run) + 1) -
                                    1;
                            if (i != query.size() - 1) {
                                curr_index = lf(i, curr_index, query[i]);
                                curr_run = this->cols[i + 1].rank_runs(
                                        curr_index);
                                symbol = get_next_char(
                                        this->cols[i + 1].zero_first,
                                        curr_run);
                                if (verbose) {
                                    std::cout << "new: " << curr_index << " "
                                              << curr_run
                                              << " "
                                              << curr_pos << " " << symbol
                                              << "\n";
                                }
                            }
                        }
                    }
                }
            }
        }
        ms_matches ms_matches;
        for (unsigned int i = 0; i < ms.len.size(); i++) {
            // TODO check if for last match
            if ((ms.len[i] != 0 && ms.len[i] > ms.len[i + 1]) ||
                (i == ms.len.size() - 1 && ms.len[i] != 0)) {
                ms_matches.basic_matches.emplace_back(ms.row[i], ms.len[i], i);
            } else if (ms.len[i] != 0 && i < ms.len.size() - 1 &&
                       ms.len[i] == ms.len[i + 1]) {
                unsigned int pos = 0;
                for (unsigned int j = i + 1; j < ms.len.size(); j++) {
                    if (ms.len[j] > ms.len[j + 1]) {
                        pos = j + 1;
                        break;
                    } else {
                        pos = i + 1;
                        break;
                    }
                }
                if (i + 1 != pos) {
                    for (unsigned int j = i; j < pos; j++) {
                        ms_matches.basic_matches.emplace_back(ms.row[j],
                                                              ms.len[j], j);
                    }
                    i = pos;
                }
                if (pos == ms.len.size() - 1) {
                    break;
                }
            }
        }

        if (extend_matches) {
            if (verbose) {
                std::cout << "extending\n";
            }
            extend_haplos(ms_matches);
        }
        if (verbose) {
            std::cout << ms << "\n";
            std::cout << ms_matches << "\n";
        }

        return ms_matches;
    }

    /**
     * @brief function to compute matching statistics matches with a given query
     * using thresholds
     * @param query haplotype query as std::string
     * @param extend_matches bool to check if extend matching statistics matches
     * with rows
     * @param verbose bool for extra prints
     * @return matching statistics matches
     * @attention use this function is enabled iff thresholds are calculated
     */
    ms_matches
    match_thr(const std::string &query, bool extend_matches = false,
              bool verbose = false) {
        if (!this->is_thr_enabled) {
            throw NoThrException{};
        }
        if (query.size() != this->panel->w) {
            throw NotEqualLengthException{};
        }
        if (extend_matches && !this->is_extended) {
            this->extend();
        }

        ms ms(query.size());
        auto curr_pos = static_cast<unsigned int>(this->cols[0].sample_end[
                this->cols[0].sample_end.size() - 1]);
        auto curr_index = curr_pos;
        unsigned int curr_run = this->cols[0].rank_runs(curr_index);
        char symbol = get_next_char(this->cols[0].zero_first, curr_run);

        for (unsigned int i = 0; i < query.size(); i++) {
            std::cout << "processed " << i << "\r";
            if (verbose) {
                std::cout << "at " << i << ": " << curr_run << " "
                          << this->cols[i].rank_thr(curr_index) << "\n";
                std::cout << curr_index << " " << curr_run << " " << curr_pos
                          << " "
                          << symbol << "\n";
            }
            if (query[i] == symbol) {
                if (verbose) {
                    std::cout << "match:\n";
                }
                // save in matching statistics row vector
                ms.row[i] = curr_pos;
                if (i != query.size() - 1) {
                    // update index, run, symbol
                    // TODO maybe make a function for this update
                    curr_index = lf(i, curr_index, query[i]);
                    curr_run = this->cols[i + 1].rank_runs(curr_index);
                    symbol = get_next_char(this->cols[i + 1].zero_first,
                                           curr_run);
                    if (verbose) {
                        std::cout << "new: " << curr_index << " " << curr_run
                                  << " "
                                  << curr_pos << " "
                                  << symbol << "\n";
                    }
                }
            } else {
                auto thr = this->cols[i].rank_thr(curr_index);
                /* TODO if index is in the position of a thresholds
                 * understand what to do
                 */
                bool in_thr = false;
                if (this->cols[i].sample_beg[curr_run] ==
                    this->cols[i].sample_end[curr_run]) {
                    in_thr = true;
                }

                if (this->cols[i].sample_beg.size() == 1) {
                    if (verbose) {
                        std::cout << "complete mismatch\n";
                    }
                    ms.row[i] = this->panel->h;
                    if (i != query.size() - 1) {
                        curr_pos = static_cast<unsigned int>(this->cols[0].sample_end[
                                this->cols[0].sample_end.size() - 1]);
                        curr_index = curr_pos;
                        curr_run = this->cols[i + 1].rank_runs(curr_index);
                        symbol = get_next_char(this->cols[i + 1].zero_first,
                                               curr_run);
                        if (verbose) {
                            std::cout << "update: " << curr_index << " "
                                      << curr_pos
                                      << " "
                                      << symbol << "\n";
                        }
                    }
                } else if ((curr_run != 0 && !in_thr && curr_run == thr) ||
                           curr_run == this->cols[i].sample_beg.size() - 1) {
                    if (verbose) {
                        std::cout << "mismatch_up: ";
                    }
                    // threshold below index so we go up
                    curr_index = (this->cols[i].select_runs(curr_run) + 1) - 1;
                    //curr_prefs = this->cols[i].rows[curr_run - 1];
                    curr_pos = static_cast<unsigned int>(this->cols[i].sample_end[
                            curr_run - 1]);
                    if (verbose) {
                        std::cout << "update: " << curr_index << " " << curr_pos
                                  << " "
                                  << symbol << "\n";
                    }
                    ms.row[i] = curr_pos;
                    if (i != query.size() - 1) {
                        curr_index = lf(i, curr_index, query[i]);
                        curr_run = this->cols[i + 1].rank_runs(curr_index);
                        symbol = get_next_char(this->cols[i + 1].zero_first,
                                               curr_run);
                        if (verbose) {
                            std::cout << "new: " << curr_index << " "
                                      << curr_run
                                      << " "
                                      << curr_pos << " " << symbol << "\n";
                        }
                    }
                } else {
                    if (verbose) {
                        std::cout << "mismatch_down: ";
                    }
                    // threshold above index so we go down
                    curr_index = (this->cols[i].select_runs(curr_run + 1) + 1);
                    curr_pos = static_cast<unsigned int>(this->cols[i].sample_beg[
                            curr_run + 1]);
                    ms.row[i] = curr_pos;
                    if (verbose) {
                        std::cout << "update: " << curr_index << " " << curr_pos
                                  << " "
                                  << symbol << "\n";
                    }
                    if (i != query.size() - 1) {
                        curr_index = lf(i, curr_index, query[i]);
                        curr_run = this->cols[i + 1].rank_runs(curr_index);
                        symbol = get_next_char(this->cols[i + 1].zero_first,
                                               curr_run);
                        if (verbose) {
                            std::cout << "new: " << curr_index << " "
                                      << curr_run
                                      << " "
                                      << curr_pos << " " << symbol << "\n";
                        }
                    }
                }
            }
        }

        int tmp_index = 0;

        for (int i = (int) ms.row.size() - 1; i >= 0; i--) {
            if (ms.row[i] == this->panel->h) {
                ms.len[i] = 0;
                continue;
            }
            if (i != (int) ms.row.size() - 1 && ms.row[i] == ms.row[i + 1]) {
                ms.len[i] = ms.len[i + 1] - 1;
                continue;
            }
            tmp_index = i;
            while (tmp_index >= 0 &&
                   query[tmp_index] == panel->getElem(ms.row[i], tmp_index)) {
                tmp_index--;
            }
            ms.len[i] = i - tmp_index;
        }
        ms_matches ms_matches;
        for (unsigned int i = 0; i < ms.len.size(); i++) {
            if ((ms.len[i] != 0 && ms.len[i] > ms.len[i + 1]) ||
                (i == ms.len.size() - 1 && ms.len[i] != 0)) {
                ms_matches.basic_matches.emplace_back(ms.row[i], ms.len[i], i);
            } else if (ms.len[i] != 0 && i < ms.len.size() - 1 &&
                       ms.len[i] == ms.len[i + 1]) {
                unsigned int pos = 0;
                for (unsigned int j = i + 1; j < ms.len.size(); j++) {
                    if (ms.len[j] > ms.len[j + 1]) {
                        pos = j + 1;
                        break;
                    } else {
                        pos = i + 1;
                        break;
                    }
                }
                if (i + 1 != pos) {
                    for (unsigned int j = i; j < pos; j++) {
                        ms_matches.basic_matches.emplace_back(ms.row[j],
                                                              ms.len[j], j);
                    }
                    i = pos;
                }
                if (pos == ms.len.size() - 1) {
                    break;
                }
            }
        }

        if (extend_matches) {
            if (verbose) {
                std::cout << "\nextending\n";
            }
            extend_haplos(ms_matches);
        }
        if (verbose) {
            std::cout << ms << "\n";
            std::cout << ms_matches << "\n";
        }
        return ms_matches;
    }

    /**
     * @brief function to compute queries with lce from a tsv file and
     * output them on a file
     * @param filename queries file
     * @param out output file
     * @param extend_matches bool to extende mathc with rows values
     * @param verbose bool for extra prints
     */
    template<typename U = ra_t>
    std::enable_if_t<sizeof(U) && (!std::is_same<ra_t, panel_ra>::value),
            void>
    match_tsv_lce(const char *filename, const char *out,
                  bool extend_matches = false, bool verbose = false) {
        std::ifstream input_matrix(filename);
        std::ofstream out_match(out);
        if (input_matrix.is_open()) {
            std::string header1;
            std::string header2;
            std::string line;
            std::string garbage;
            std::string new_column;
            getline(input_matrix, line);
            getline(input_matrix, line);
            std::vector<std::string> queries_panel;
            while (getline(input_matrix, line) && !line.empty()) {
                std::istringstream is_col(line);
                is_col >> garbage;
                if (garbage == "TOTAL_SAMPLES:") {
                    break;
                }
                is_col >> garbage >> garbage >> garbage >> new_column;
                queries_panel.push_back(new_column);
            }
            input_matrix.close();
            std::string query;
            if (out_match.is_open()) {
                for (unsigned int i = 0; i < queries_panel.size(); i++) {
                    query = queries_panel[i];
                    ms_matches matches;
                    matches = this->match_lce(query, extend_matches, verbose);
                    if (verbose) {
                        std::cout << i << ": ";
                    }
                    out_match << i << ": ";

                    if (verbose) {
                        std::cout << matches;
                    }
                    out_match << matches;

                    if (verbose) {
                        std::cout << "\n";
                    }
                    out_match << "\n";
                    query.clear();
                }
                out_match.close();
            } else {
                throw FileNotFoundException{};
            }

        } else {
            throw FileNotFoundException{};
        }
    }

    /**
     * @brief function to compute queries with lce from a transposed tsv
     * file and output them on a file
     * @param filename queries file
     * @param out output file
     * @param extend_matches bool to extende mathc with rows values
     * @param verbose bool for extra prints
     */
    template<typename U = ra_t>
    std::enable_if_t<sizeof(U) && (!std::is_same<ra_t, panel_ra>::value),
            void>
    match_tsv_tr_lce(const char *filename, const char *out,
                     bool extend_matches = false,
                     bool verbose = false) {
        std::ifstream input_matrix(filename);
        std::ofstream out_match(out);
        if (input_matrix.is_open()) {
            std::string header1;
            std::string header2;
            std::string line;
            std::string garbage;
            std::string new_column;
            getline(input_matrix, line);
            getline(input_matrix, line);
            std::vector<std::string> queries_panel;
            while (getline(input_matrix, line) && !line.empty()) {
                std::istringstream is_col(line);
                is_col >> garbage;
                if (garbage == "TOTAL_SAMPLES:") {
                    break;
                }
                is_col >> garbage >> garbage >> garbage >> new_column;
                queries_panel.push_back(new_column);
            }
            input_matrix.close();
            std::string query = "";
            if (out_match.is_open()) {
                for (unsigned int i = 0; i < queries_panel[0].size(); i++) {
                    for (auto &j: queries_panel) {
                        query.push_back(j[i]);
                    }
                    if (verbose) {
                        std::cout << query << "\n";
                    }
                    ms_matches matches;

                    matches = this->match_lce(query, extend_matches, verbose);
                    if (verbose) {
                        std::cout << i << ": ";
                    }
                    out_match << i << ": ";

                    if (verbose) {
                        std::cout << matches;
                    }
                    out_match << matches;

                    if (verbose) {
                        std::cout << "\n";
                    }
                    out_match << "\n";
                    query.clear();
                }
                out_match.close();
            } else {
                throw FileNotFoundException{};
            }

        } else {
            throw FileNotFoundException{};
        }
    }

    /**
     * @brief function to compute queries with thresholds  from a tsv file and
     * output them on a file
     * @param filename queries file
     * @param out output file
     * @param extend_matches bool to extende mathc with rows values
     * @param verbose bool for extra prints
     */
    void
    match_tsv_thr(const char *filename, const char *out,
                  bool extend_matches = false, bool verbose = false) {
        std::ifstream input_matrix(filename);
        std::ofstream out_match(out);
        if (input_matrix.is_open()) {
            std::string header1;
            std::string header2;
            std::string line;
            std::string garbage;
            std::string new_column;
            getline(input_matrix, line);
            getline(input_matrix, line);
            std::vector<std::string> queries_panel;
            while (getline(input_matrix, line) && !line.empty()) {
                std::istringstream is_col(line);
                is_col >> garbage;
                if (garbage == "TOTAL_SAMPLES:") {
                    break;
                }
                is_col >> garbage >> garbage >> garbage >> new_column;
                queries_panel.push_back(new_column);
            }
            input_matrix.close();
            std::string query;
            if (out_match.is_open()) {
                for (unsigned int i = 0; i < queries_panel.size(); i++) {
                    query = queries_panel[i];
                    ms_matches matches;
                    matches = this->match_thr(query, extend_matches, verbose);
                    if (verbose) {
                        std::cout << i << ": ";
                    }
                    out_match << i << ": ";

                    if (verbose) {
                        std::cout << matches;
                    }
                    out_match << matches;

                    if (verbose) {
                        std::cout << "\n";
                    }
                    out_match << "\n";
                    query.clear();
                }
                out_match.close();
            } else {
                throw FileNotFoundException{};
            }

        } else {
            throw FileNotFoundException{};
        }
    }

    /**
     * @brief function to compute queries with thresholds from a transposed tsv
     * file and output them on a file
     * @param filename queries file
     * @param out output file
     * @param extend_matches bool to extende mathc with rows values
     * @param verbose bool for extra prints
     */
    void
    match_tsv_tr_thr(const char *filename, const char *out,
                     bool extend_matches = false,
                     bool verbose = false) {
        std::ifstream input_matrix(filename);
        std::ofstream out_match(out);
        if (input_matrix.is_open()) {
            std::string header1;
            std::string header2;
            std::string line;
            std::string garbage;
            std::string new_column;
            getline(input_matrix, line);
            getline(input_matrix, line);
            std::vector<std::string> queries_panel;
            while (getline(input_matrix, line) && !line.empty()) {
                std::istringstream is_col(line);
                is_col >> garbage;
                if (garbage == "TOTAL_SAMPLES:") {
                    break;
                }
                is_col >> garbage >> garbage >> garbage >> new_column;
                queries_panel.push_back(new_column);
            }
            input_matrix.close();
            std::string query = "";
            if (out_match.is_open()) {
                for (unsigned int i = 0; i < queries_panel[0].size(); i++) {
                    for (auto &j: queries_panel) {
                        query.push_back(j[i]);
                    }
                    if (verbose) {
                        std::cout << query << "\n";
                    }
                    ms_matches matches;

                    matches = this->match_thr(query, extend_matches, verbose);
                    if (verbose) {
                        std::cout << i << ": ";
                    }
                    out_match << i << ": ";

                    if (verbose) {
                        std::cout << matches;
                    }
                    out_match << matches;

                    if (verbose) {
                        std::cout << "\n";
                    }
                    out_match << "\n";
                    query.clear();
                }
                out_match.close();
            } else {
                throw FileNotFoundException{};
            }

        } else {
            throw FileNotFoundException{};
        }
    }

    /**
     * @brief function to obtain size in bytes of the matching statistics
     * supported RLPBWT
     * @param verbose bool for extra prints
     * @return size in bytes
    */
    unsigned long long size_in_bytes(bool verbose = false) {
        unsigned long long size = 0;
        unsigned long long size_run = 0;
        unsigned long long size_thr = 0;
        unsigned long long size_u = 0;
        unsigned long long size_v = 0;
        unsigned long long size_samples = 0;
        auto lp_size = sdsl::size_in_bytes(this->last_pref);
        size += lp_size;
        for (unsigned int i = 0; i < this->cols.size(); ++i) {
            size += this->cols[i].size_in_bytes();
            size_run += sdsl::size_in_bytes(this->cols[i].runs) +
                        sdsl::size_in_bytes(this->cols[i].rank_runs) +
                        sdsl::size_in_bytes(this->cols[i].select_runs);
            size_thr += sdsl::size_in_bytes(this->cols[i].thr) +
                        sdsl::size_in_bytes(this->cols[i].rank_thr) +
                        sdsl::size_in_bytes(this->cols[i].select_thr);
            size_u += sdsl::size_in_bytes(this->cols[i].u) +
                      sdsl::size_in_bytes(this->cols[i].rank_u) +
                      sdsl::size_in_bytes(this->cols[i].select_u);
            size_v += sdsl::size_in_bytes(this->cols[i].v) +
                      sdsl::size_in_bytes(this->cols[i].rank_v) +
                      sdsl::size_in_bytes(this->cols[i].select_v);
            size_samples += sdsl::size_in_bytes(this->cols[i].sample_beg) +
                            sdsl::size_in_bytes(this->cols[i].sample_end) +
                            lp_size;
        }
        if (verbose) {
            std::cout << "run: " << size_run << " bytes\n";
            std::cout << "thr: " << size_thr << " bytes\n";
            std::cout << "u: " << size_u << " bytes\n";
            std::cout << "v: " << size_v << " bytes\n";
            std::cout << "samples: " << size_samples << " bytes\n";
            std::cout
                    << "rlpbwt (with also c values and other support variables): "
                    << size << " bytes\n";
        }
        size += (sizeof(bool) * 2);
        size += (sizeof(unsigned int) * 2);
        size += this->panel->size_in_bytes(verbose);
        if (this->is_extended) {
            size += this->phi->size_in_bytes(verbose);
        }

        return size;
    }

    /**
     * @brief function to obtain size in megabytes of the matching statistics
     * supported RLPBWT
     * @param verbose bool for extra prints
     * @return size in megabytes
    */
    double size_in_mega_bytes(bool verbose = true) {
        double size = 0;
        double to_mega = ((double) 1 / (double) 1024) / (double) 1024;
        double size_run = 0;
        double size_thr = 0;
        double size_u = 0;
        double size_v = 0;
        double size_samples = 0;
        auto lp_size = sdsl::size_in_mega_bytes(this->last_pref);
        size += lp_size;
        for (unsigned int i = 0; i < this->cols.size(); ++i) {
            size += this->cols[i].size_in_mega_bytes();
            size_run += sdsl::size_in_mega_bytes(this->cols[i].runs) +
                        sdsl::size_in_mega_bytes(this->cols[i].rank_runs) +
                        sdsl::size_in_mega_bytes(this->cols[i].select_runs);
            size_thr += sdsl::size_in_mega_bytes(this->cols[i].thr) +
                        sdsl::size_in_mega_bytes(this->cols[i].rank_thr) +
                        sdsl::size_in_mega_bytes(this->cols[i].select_thr);
            size_u += sdsl::size_in_mega_bytes(this->cols[i].u) +
                      sdsl::size_in_mega_bytes(this->cols[i].rank_u) +
                      sdsl::size_in_mega_bytes(this->cols[i].select_u);
            size_v += sdsl::size_in_mega_bytes(this->cols[i].v) +
                      sdsl::size_in_mega_bytes(this->cols[i].rank_v) +
                      sdsl::size_in_mega_bytes(this->cols[i].select_v);
            size_samples += sdsl::size_in_mega_bytes(this->cols[i].sample_beg) +
                            sdsl::size_in_mega_bytes(this->cols[i].sample_end) +
                            lp_size;
        }
        if (verbose) {
            std::cout << "run: " << size_run << " megabytes\n";
            std::cout << "thr: " << size_thr << " megabytes\n";
            std::cout << "u: " << size_u << " megabytes\n";
            std::cout << "v: " << size_v << " megabytes\n";
            std::cout << "samples: " << size_samples << " megabytes\n";
            std::cout
                    << "rlpbwt (with also c values and other support variables): "
                    << size << " megabytes\n";
        }
        size += (sizeof(bool) * 2 * to_mega);
        size += (sizeof(unsigned int) * 2 * to_mega);
        size += this->panel->size_in_mega_bytes(verbose);
        if (this->is_extended) {
            size += this->phi->size_in_mega_bytes(verbose);
        }

        return size;
    }

    /**
     * @brief function to serialize the matching statistics supported RLPBWT
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
        out.write((char *) &this->is_thr_enabled, sizeof(this->is_thr_enabled));
        written_bytes += sizeof(this->is_thr_enabled);
        out.write((char *) &this->is_extended, sizeof(this->is_extended));
        written_bytes += sizeof(this->is_extended);

        written_bytes += this->panel->serialize(out, child, "panel");
        for (unsigned int i = 0; i < this->cols.size(); i++) {
            std::string label = "col_" + std::to_string(i);
            written_bytes += this->cols[i].serialize(out, child, label);
        }
        written_bytes += this->last_pref.serialize(out, child, "last_pref");
        if (this->is_extended) {
            written_bytes += this->phi->serialize(out, child, "phi");
        }
        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    /**
     * @brief function to load the matching statistics supported RLPBWT object
     * @param in std::istream object from which load the matching statistics
     * supported RLPBWT structure object
     */
    void load(std::istream &in, const char *slp_filename = "") {
        in.read((char *) &this->is_thr_enabled, sizeof(this->is_thr_enabled));
        in.read((char *) &this->is_extended, sizeof(this->is_extended));
        if constexpr (!std::is_same_v<ra_t, panel_ra>) {
            if (std::string(slp_filename).empty()) {
                throw SlpNotFoundException{};
            }
        }
        if constexpr (!std::is_same_v<ra_t, panel_ra>) {
            auto _panel = new ra_t();
            _panel->load(in, slp_filename);
            this->panel = _panel;
        } else {
            auto _panel = new ra_t();
            _panel->load(in);
            this->panel = _panel;
        }
        for (unsigned int i = 0; i < this->panel->w; i++) {
            auto c = new column_ms();
            c->load(in);
            this->cols.emplace_back(*c);
        }
        this->last_pref.load(in);
        if (this->is_extended) {
            auto _phi = new phi_support<ra_t>;
            _phi->load(in);
            this->phi = _phi;
        }
    }
};


#endif //RLPBWT_RLPBWT_MS_H
