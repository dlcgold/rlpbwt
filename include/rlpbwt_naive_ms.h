//
// Created by dlcgold on 16/08/22.
//

#ifndef RLPBWT_RLPBWT_NAIVE_MS_H
#define RLPBWT_RLPBWT_NAIVE_MS_H

#include <vector>
#include <type_traits>
#include "column_naive_ms.h"
#include "exceptions.h"
#include "panel_ra.h"
#include "slp_panel_ra.h"
#include "phi_support.h"
#include "ms.h"
#include "ms_matches.h"


/**
 * @brief data structure for matching-statistics supported RLPBWT naive
 * @tparam ra_t type of panel (panel_ra or slp_panel_ra)
 */
template<typename ra_t>
class rlpbwt_naive_ms {
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
    * @brief vector of matcing statistics supported naive columns
    */
    std::vector<column_naive_ms> cols;

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
    column_naive_ms
    build_column(std::string &column, std::vector<unsigned int> &pref,
                 sdsl::int_vector<> &div) {
        unsigned int height = pref.size();
        // variable for "c" value
        unsigned int count0 = 0;
        // variables for compute "u" and "v" values
        unsigned int u = 0;
        unsigned int v = 0;
        // temporary variables for compute "u" and "v" values
        unsigned int count0tmp = 0;
        unsigned int count1 = 0;
        // bool to check first symbol fo the column
        bool start = true;

        // update start and "c" value
        for (unsigned int i = 0; i < height; i++) {
            if (i == 0 && column[pref[i]] == '1') {
                start = false;
            }
            if (column[i] == '0') {
                count0++;
            }
        }

        // temporary variable to compute thresholds
        unsigned int lcs = 0;

        // initialize a vector of pair in order to build final sdsl int_vector for
        // p and u/v
        std::vector<std::pair<unsigned int, unsigned int>> rows;
        std::vector<unsigned int> thr;
        // support vector to store prefix array samples
        std::vector<std::pair<unsigned int, unsigned int>> samples;

        // temporary variable for p
        unsigned int p_tmp = 0;

        // temporary variable for prefix array value at the begin of a run
        unsigned int tmp_beg = 0;

        // temporary variable for thresholds
        unsigned int tmp_thr = 0;

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
        for (unsigned int i = 0; i < height; i++) {
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

            // updating thresholds (iff thresholds are required)
            if (this->is_thr_enabled) {
                if ((i == 0) || (column[pref[i]] != column[pref[i - 1]])) {
                    tmp_thr = i;
                    lcs = div[i];
                }
                if (div[i] < lcs) {
                    tmp_thr = i;
                    lcs = div[i];
                }
            }

            // record starting position of a run
            if ((i == 0) || (column[pref[i]] != column[pref[i - 1]])) {
                p_tmp = i;
            }

            // do stuff at run change
            if ((i == height - 1) || (column[pref[i]] != column[pref[i + 1]])) {
                // update vector for p and u/v and swap the case to study in
                // next run
                if (pusho) {
                    rows.emplace_back(p_tmp, v);
                    std::swap(pusho, pushz);
                } else {
                    rows.emplace_back(p_tmp, u);
                    std::swap(pusho, pushz);
                }
                if (this->is_thr_enabled) {
                    if (i + 1 != height && div[i + 1] < div[tmp_thr]) {
                        thr.push_back(i + 1);
                    } else {
                        thr.push_back(tmp_thr);
                    }
                }
                samples.emplace_back(tmp_beg, pref[i]);
                begrun = true;
            }
        }
        // create and compress sdsl int_vector
        if (this->is_thr_enabled) {
            sdsl::int_vector<> t_vec(rows.size());
            sdsl::int_vector<> p_vec(rows.size());
            sdsl::int_vector<> uv_vec(rows.size());
            sdsl::int_vector<> sb_vec(rows.size());
            sdsl::int_vector<> se_vec(rows.size());
            for (unsigned int i = 0; i < rows.size(); i++) {
                p_vec[i] = rows[i].first;
                uv_vec[i] = rows[i].second;
                t_vec[i] = thr[i];
                sb_vec[i] = samples[i].first;
                se_vec[i] = samples[i].second;
            }
            sdsl::util::bit_compress(p_vec);
            sdsl::util::bit_compress(uv_vec);
            sdsl::util::bit_compress(t_vec);
            sdsl::util::bit_compress(sb_vec);
            sdsl::util::bit_compress(se_vec);
            // return the column
            return {start, count0, p_vec, uv_vec, t_vec, sb_vec, se_vec};
        } else {
            sdsl::int_vector<> p_vec(rows.size());
            sdsl::int_vector<> uv_vec(rows.size());
            sdsl::int_vector<> sb_vec(rows.size());
            sdsl::int_vector<> se_vec(rows.size());
            for (unsigned int i = 0; i < rows.size(); i++) {
                p_vec[i] = rows[i].first;
                uv_vec[i] = rows[i].second;
                sb_vec[i] = samples[i].first;
                se_vec[i] = samples[i].second;
            }
            sdsl::util::bit_compress(p_vec);
            sdsl::util::bit_compress(uv_vec);
            sdsl::util::bit_compress(sb_vec);
            sdsl::util::bit_compress(se_vec);
            sdsl::int_vector<> t_vec;
            sdsl::util::bit_compress(t_vec);
            // return the column
            return {start, count0, p_vec, uv_vec, t_vec, sb_vec, se_vec};
        }
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
        // get run
        unsigned int run_index = index_to_run(index, col_index);
        //get offset
        unsigned int offset = index - this->cols[col_index].p[run_index];
        // undoing the offsets when they are wrong/useless
        if ((symbol == '0' && get_next_char(this->cols[col_index].zero_first,
                                            run_index) == '1') ||
            (symbol == '1' && get_next_char(this->cols[col_index].zero_first,
                                            run_index) == '0')) {
            offset = 0;
        }
        // obtain "u" and "v"
        auto uv = uvtrick(col_index, run_index);
        if (verbose) {
            std::cout << uv.first << ", " << uv.second << "\n";
        }
        // fix for the last index that's "outside" the column
        if (this->cols[col_index].p[run_index] + offset == this->height) {
            if (get_next_char(this->cols[col_index].zero_first, run_index) ==
                '0') {
                uv.second--;
            } else {
                uv.first--;
            }
        }
        // computer lf-mapping as Durbin's w(i, s), eventually using offset
        if (symbol == '0') {
            return uv.first + offset;
        } else {
            return this->cols[col_index].count_0 + uv.second + offset;
        }
    }

    /**
     * @brief function to map an index to the correct run in a column
     * @param index index to map
     * @param col_index column index
     * @return run index
     */
    unsigned int
    index_to_run(unsigned int index, unsigned int col_index) const {
        // if requested index is equal or greater than the p value of the last run
        // return the index of the last run
        if (index >=
            this->cols[col_index].p[this->cols[col_index].p.size() - 1]) {
            return this->cols[col_index].p.size() - 1;
        }

        // binary search to compute run index
        unsigned int bi = 0;
        unsigned int e = this->cols[col_index].p.size();
        unsigned int pos = (e - bi) / 2;
        while (pos != e && this->cols[col_index].p[pos] != index) {
            if (index < (unsigned int) this->cols[col_index].p[pos]) {
                e = pos;
            } else {
                if (pos + 1 == e ||
                    (unsigned int) this->cols[col_index].p[pos + 1] > index) {
                    break;
                }
                bi = pos + 1;
            }
            pos = bi + (e - bi) / 2;
        }
        return pos;
    }

    /**
     * @brief trick to extract u and v value from a run in rlpbwt column
     * @param col_index index of the column
     * @param index virtual index of the row of the original panel
     * @return a std::pair with u as first and v as second
     */
    std::pair<unsigned int, unsigned int>
    uvtrick(unsigned int col_index, unsigned int run_index) const {
        unsigned int u;
        unsigned int v;
        // if run index is 0 u = v = 0
        // in other case, based on first symbol of the column
        // we have u/v in the same row and v/u in the previous one
        if (run_index == 0) {
            u = 0;
            v = 0;
        } else if (run_index % 2 == 0) {
            u = this->cols[col_index].uv[run_index - 1];
            v = this->cols[col_index].uv[run_index];
            if (!this->cols[col_index].zero_first) {
                std::swap(u, v);
            }
        } else {
            u = this->cols[col_index].uv[run_index];
            v = this->cols[col_index].uv[run_index - 1];
            if (!this->cols[col_index].zero_first) {
                std::swap(u, v);
            }
        }
        return {u, v};
    }

    /**
     * @brief function to extend matching statistics matches with matching
     * rows indices
     * @param ms_matches matching statistics matches that will be extended
     */
    void extend_haplos(ms_matches &ms_matches) {
        // iterate over every basic match
        for (unsigned int i = 0; i < ms_matches.basic_matches.size(); i++) {
            // initialize the vector that will contain row indices
            std::vector<unsigned int> haplos;

            // extract information from the current basic match
            unsigned int start_row = std::get<0>(
                    ms_matches.basic_matches[i]);
            haplos.emplace_back(start_row);
            unsigned int curr_len = std::get<1>(
                    ms_matches.basic_matches[i]);
            unsigned int curr_col = std::get<2>(
                    ms_matches.basic_matches[i]);

            // initialize boolean and temporary variables for go up/down in
            // search of matching rows
            bool check_down = true;
            unsigned int down_row = 0;
            bool check_up = true;
            unsigned int up_row = 0;
            // go down/up and add row that has a lce of at least current len in
            // matching statistics (if panel is not saved as SLP the computation
            // of lce is simulated with random access to the panel)
            if constexpr (std::is_same_v<ra_t, slp_panel_ra>) {
                // down
                while (check_down) {
                    auto phi_res = this->phi->phi_inv(start_row, curr_col + 1);
                    if (!phi_res.has_value()) {
                        check_down = false;
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
                start_row = std::get<0>(ms_matches.basic_matches[i]);
                // up
                while (check_up) {
                    auto phi_res = this->phi->phi(start_row, curr_col + 1);
                    if (!phi_res.has_value()) {
                        check_up = false;
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
                // down
                while (check_down) {
                    auto phi_res = this->phi->phi_inv(start_row, curr_col + 1);
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
                start_row = std::get<0>(ms_matches.basic_matches[i]);
                // up
                while (check_up) {
                    auto phi_res = this->phi->phi(start_row, curr_col + 1);
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
            // record the haplotypes
            ms_matches.haplos.emplace_back(haplos);
        }
    }

/**
     * @brief function to check if longest common extension between two rows is
     * equal or greater than a bound ending at a given column
     * @param col ending column column
     * @param curr current row
     * @param other other row
     * @param bound bound to check
     * @param verbose bool for extra prints
     * @return true if longest common extension between two rows is
     * equal or greater than the bound
     * @attention this function is enabled iff the panel is an SLP
     */
    template<typename U = ra_t>
    std::enable_if_t<sizeof(U) && (!std::is_same<ra_t, panel_ra>::value), bool>
    lceBound(unsigned int col, unsigned int curr, unsigned int other,
             unsigned int bound, bool verbose = false) {
        // by design if we are at first column we don't compute any lce
        if (col == 0) {
            return false;
        }
        auto w_l = (unsigned long long int) this->panel->w;
        // obtain the column in the reverse order (as the SLP is saved)
        unsigned long long int rev_col =
                (w_l - (unsigned long long int) 1) - col;
        // obtain indices of the symbols required in the SLP (built from the
        // reverse matrix saved as a single line string)
        unsigned long long int pos_curr = rev_col + (w_l * curr);
        unsigned long long int pos_other = rev_col + (w_l * other);
        if (verbose) {
            std::cout << "at " << rev_col << ": " << pos_curr << ", "
                      << pos_other << "\n";
        }
        // compute lce and return true if equal (or greater) than the bound
//        auto lcp_pair_naive = lceToR_NaiveBounded(this->panel->panel,
//                                                          pos_other,
//                                                          pos_curr, bound);
        // TODO check the bound in no naive version, it doesn't work
        auto lcp_pair = lceToRBounded(this->panel->panel,
                                      pos_other,
                                      pos_curr, bound);
        // TODO if the bound works it should be lcp_pair == bound
        if (lcp_pair >= bound) {
            return true;
        }
        return false;
    }

    /**
     * @brief function to compute longest common extension between two rows
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
        // by design if we are at first column we don't compute any lce
        if (col == 0) {
            return std::make_pair(other, 0);
        }
        auto w_l = (unsigned long long int) this->panel->w;
        // obtain the column in the reverse order (as the SLP is saved)
        // the previous one in order to compute lce for MS
        unsigned long long int rev_col =
                (w_l - (unsigned long long int) 1) - col +
                (unsigned long long int) 1;
        // obtain indices of the symbols required in the SLP (built from the
        // reverse matrix saved as a single line string)
        unsigned long long int pos_curr = rev_col + (w_l * curr);
        unsigned long long int pos_other = rev_col + (w_l * other);
        if (verbose) {
            std::cout << "at " << rev_col << ": " << pos_curr << ", "
                      << pos_other << "\n";
        }
        // compute lce (eventually bounded with the column index value) and
        // return the length (and the prefix value of the other position)
        auto lcp_other = lceToR(this->panel->panel, pos_other, pos_curr);
        if (verbose) {
            std::cout << "at " << rev_col << " lce is : " << lcp_other << "\n";
        }
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
        // by design if we are at first column we don't compute any lce
        if (col == 0) {
            return std::make_pair(next, 0);
        }

        auto w_l = (unsigned long long int) this->panel->w;
        // obtain the column in the reverse order (as the SLP is saved)
        // the previous one in order to compute lce for MS
        unsigned long long int rev_col =
                (w_l - (unsigned long long int) 1) - col +
                (unsigned long long int) 1;

        // obtain indices of the symbols required in the SLP (built from the
        // reverse matrix saved as a single line string)
        unsigned long long int pos_curr = rev_col + (w_l * curr);
        unsigned long long int pos_prev = rev_col + (w_l * prev);
        unsigned long long int pos_next = rev_col + (w_l * next);
        if (verbose) {
            std::cout << "at " << rev_col << ": " << pos_curr << ", "
                      << pos_prev
                      << ", " << pos_next << "\n";
        }
        // compute lce for both requested indices (eventually bounded with the
        // column index value)
        auto lcp_prev = lceToR(this->panel->panel, pos_prev,
                               pos_curr);
        if (lcp_prev >= col) {
            lcp_prev = col;
        }
        auto lcp_next = lceToR(this->panel->panel, pos_curr,
                               pos_next);
        if (lcp_next >= col) {
            lcp_next = col;
        }
        // select the greater lce and return the length plus the relative prefix
        // array index
        if (lcp_prev > lcp_next) {
            return std::make_pair(prev, lcp_prev);
        } else {
            return std::make_pair(next, lcp_next);
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
    rlpbwt_naive_ms() = default;

    /**
     * @brief default destructor
     */
    ~rlpbwt_naive_ms() {
        delete panel;
        if (this->is_extended) {
            delete phi;
        }
    }

    /**
     * @brief constructor of a RLPBWT that support matching statistics
     * @param filename file with the panel
     * @param thr bool to enable thresholds computation
     * @param verbose bool fro extra prints
     * @param slp_filename file with the slp of the panel
     */
    explicit rlpbwt_naive_ms(const char *filename, bool thr = false,
                             bool verbose = false,
                             const char *slp_filename = "") {
        std::ifstream input_matrix(filename);
        if constexpr (!std::is_same_v<ra_t, panel_ra>) {
            if (std::string(slp_filename).empty()) {
                throw SlpNotFoundException{};
            }
        }
        if constexpr (std::is_same_v<ra_t, panel_ra>) {
            thr = true;
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
            /*unsigned int tmp_width = std::count(
                    std::istreambuf_iterator<char>(input_matrix),
                    std::istreambuf_iterator<char>(), '\n');
                    */
            auto tmp_width = 1;
            while (getline(input_matrix, line) && !line.empty()) {
                std::istringstream is_col(line);
                is_col >> garbage;
                if (garbage == "TOTAL_SAMPLES:") {
                    break;
                }
                tmp_width++;
            }
            std::cout << "w: " << tmp_width << "\n";
            this->width = tmp_width;
            this->height = tmp_height;
            input_matrix.clear();
            input_matrix.seekg(0, std::ios::beg);
            this->cols = std::vector<column_naive_ms>(tmp_width + 1);
            std::vector<unsigned int> pref(tmp_height);
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
            std::string last_column;
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
                    std::cout << new_column << "\n" << this->cols[count]
                              << "\n-------------------------------\n";
                }
                auto col = rlpbwt_naive_ms::build_column(new_column, pref, div);
                if constexpr (std::is_same_v<ra_t, panel_ra>) {
                    for (unsigned int k = 0; k < new_column.size(); k++) {
                        if (new_column[k] != '0') {
                            this->panel->panel[count][k] = true;
                        }
                    }
                }
                this->cols[count] = col;
                rlpbwt_naive_ms::update(new_column, pref, div);
                last_col = new_column;
                count++;
            }
            for (unsigned int i = 0; i < pref.size(); i++) {
                this->last_pref[i] = pref[i];
            }

            auto col = rlpbwt_naive_ms::build_column(last_col, pref, div);
            this->cols[count] = col;
            sdsl::util::bit_compress(this->last_pref);
            this->is_extended = false;
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
     * @brief function to compute matching statistics matches with a given query
     * using lce queries
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
        // compute the match iff |query| is equal to the width of the panel
        if (query.size() != this->panel->w) {
            std::cout << query.size() << " != " << this->panel->w << "\n";
            throw NotEqualLengthException{};
        }
        // if required extend with the phi support struct (iff not already
        // extended)
        if (extend_matches && !this->is_extended) {
            this->extend();
        }

        // initialize matching statistics
        //ms ms(query.size());
        unsigned int p = 0;
        unsigned int l = 0;
        unsigned int prep = 0;
        unsigned int prel = 0;

        // algorithm begin from the last row of the first column
        // so we obtain the prefix array value (from the samples), the run index
        // and the relative symbol
        auto curr_pos = static_cast<unsigned int>(this->cols[0].sample_end[
                this->cols[0].sample_end.size() - 1]);
        auto curr_index = curr_pos;
        unsigned int curr_run = index_to_run(curr_index, 0);
        char symbol = get_next_char(this->cols[0].zero_first, curr_run);


        // initialize struct for matches
        ms_matches ms_matches;

        // iterate over every query's symbol/column index
        for (unsigned int i = 0; i < query.size(); i++) {
            //std::cout << "processed " << i << "\r";
            if (verbose) {
                std::cout << "at " << i << ": " << curr_run << "\n";
                std::cout << curr_index << " " << curr_run << " " << curr_pos
                          << " "
                          << symbol << "\n";
            }
            // a lot of cases:
            // - if the pointer in the RLPBWT match the symbol in the query
            //   we proceed simply using lf-mapping (if we are not at the end)
            // - if the pointer in the RLPBWT mismatch the symbol in the query
            //   and we have only that symbol in the column we restart from the
            //   next column at last index whit the relative prefix array value
            // - otherwise we proceeed using lce queries to select the best
            //   symbol, between the previous and next good symbol (if the
            //   exists), to jump
            if (query[i] == symbol) {
                if (verbose) {
                    std::cout << "match:\n";
                }
                // report in matching statistics row/len vector
                //ms.row[i] = curr_pos;
                p = curr_pos;
                /*if (i != 0 && ms.row[i - 1] == ms.row[i]) {
                    ms.len[i] = ms.len[i - 1] + 1;

                } else {
                    ms.len[i] = 1;
                }*/
                if (i != 0 && prep == p) {
                    l = prel + 1;
                } else {
                    l = 1;
                }
                // update index, run, symbol if we are not at the end
                if (i != query.size() - 1) {
                    curr_index = lf(i, curr_index, query[i]);
                    curr_run = index_to_run(curr_index, i + 1);
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
                    // report in matching statistics row vector using panel
                    // height as sentinel
                    //ms.row[i] = this->panel->h;
                    p = this->panel->h;
                    // report in matching statistics len vector
                    l = 0;
                    // ms.len[i] = 0;
                    // update index, run, symbol (as explained before) if we are
                    // not at the end
                    if (i != query.size() - 1) {
                        curr_pos = static_cast<unsigned int>(this->cols[i +
                                                                        1].sample_end[
                                this->cols[i + 1].sample_end.size() - 1]);
                        curr_index = this->panel->h - 1;

                        curr_run = index_to_run(curr_index, i + 1);
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
                    // if we are in the last run we check lce only with the
                    // previous correct symbol
                    if (curr_run == this->cols[i].sample_beg.size() - 1) {
                        // select symbol
                        curr_index = this->cols[i].p[curr_run] - 1;
                        auto prev_pos = static_cast<unsigned int>(this->cols[i].sample_end[
                                curr_run - 1]);
                        if (verbose) {
                            std::cout << "end_run with " << curr_pos << ", "
                                      << prev_pos << "\n";
                        }
                        // compute lce
                        auto lce_value = lce(i, curr_pos, prev_pos, false);
                        // report in matching statistics row/len vector
                        p = prev_pos;
                        if (i == 0) {
                            l = 1;
                        } else {
                            if (prel == 0) {
                                l = 1;
                            } else {
                                l = std::min(prel, lce_value.second) + 1;
                            }
                        }
                        // update current position
                        curr_pos = prev_pos;
                        if (verbose) {
                            std::cout << "update: " << curr_index << " "
                                      << curr_pos
                                      << " "
                                      << symbol << "\n";
                        }
                        // update index, run, symbol if we are not at the end
                        if (i != query.size() - 1) {
                            curr_index = lf(i, curr_index, query[i]);
                            curr_run = index_to_run(curr_index, i + 1);
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
                        // if we are in the first run we check lce only with the
                        // next correct symbol

                        // select index
                        curr_index = this->cols[i].p[curr_run + 1];
                        auto next_pos = static_cast<unsigned int>(this->cols[i].sample_beg[
                                curr_run + 1]);
                        if (verbose) {
                            std::cout << "first_run with " << curr_pos << ", "
                                      << next_pos << "\n";
                        }
                        // compute lce
                        auto lce_value = lce(i, curr_pos, next_pos, false);
                        // report in matching statistics row/len vector
                        p = next_pos;
                        if (i == 0) {
                            l = 1;
                        } else {
                            if (prel == 0) {
                                l = 1;
                            } else {
                                l = std::min(prel, lce_value.second) + 1;
                            }
                        }
                        // update current position
                        curr_pos = next_pos;
                        if (verbose) {
                            std::cout << "update: " << curr_index << " "
                                      << curr_pos
                                      << " "
                                      << symbol << "\n";
                        }
                        // update index, run, symbol if we are not at the end
                        if (i != query.size() - 1) {
                            curr_index = lf(i, curr_index, query[i]);
                            curr_run = index_to_run(curr_index, i + 1);
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
                        // in every other case we have to choice the best
                        // solution

                        // select the two indices of the previous/next correct
                        // symbol
                        auto prev_pos = static_cast<unsigned int>(this->cols[i].sample_end[
                                curr_run - 1]);
                        auto next_pos = static_cast<unsigned int>(this->cols[i].sample_beg[
                                curr_run + 1]);
                        if (verbose) {
                            std::cout << "curr: " << curr_pos << ", prev: "
                                      << prev_pos << ", next: " << next_pos
                                      << " -> ";
                        }
                        // compute lce
                        auto lce = this->lce_pair(i, curr_pos, prev_pos,
                                                  next_pos,
                                                  false);
                        if (verbose) {
                            std::cout << lce.first << ", " << lce.second
                                      << "\n";
                        }

                        // select best choice, report in matching statistics
                        // row/len vector and update index, run, symbol if we
                        // are not at the end
                        if (lce.first == next_pos) {
                            curr_pos = next_pos;
                            p = curr_pos;
                            if (i == 0) {
                                l = 1;
                            } else {
                                if (prel == 0) {
                                    l = 1;
                                } else {
                                    l = std::min(prel, lce.second) + 1;
                                }
                            }
                            curr_index = this->cols[i].p[curr_run + 1];
                            if (i != query.size() - 1) {
                                curr_index = lf(i, curr_index, query[i]);
                                curr_run = index_to_run(curr_index, i + 1);
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
                            p = curr_pos;
                            if (i == 0) {
                                l = 1;
                            } else {
                                if (prel == 0) {
                                    l = 1;
                                } else {
                                    l = std::min(prel, lce.second) + 1;
                                }
                            }
                            curr_index = this->cols[i].p[curr_run] - 1;
                            if (i != query.size() - 1) {
                                curr_index = lf(i, curr_index, query[i]);
                                curr_run = index_to_run(curr_index, i + 1);
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
            if (i != 0 && (prel > 0 && prel >= l)) {
                ms_matches.basic_matches.emplace_back(prep, prel, i - 1);
            }
            if (i == query.size() - 1 && l != 0) {
                ms_matches.basic_matches.emplace_back(p, l, i);
            }
            prep = p;
            prel = l;
        }
        // compute every row that are matching if required
        if (extend_matches) {
            if (verbose) {
                std::cout << "extending\n";
            }
            extend_haplos(ms_matches);
        }
        if (verbose) {
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
        // TODO fix some bugs using panel_ra
        // if RLPBWT is built without thresholds throw relative exception
        if (!this->is_thr_enabled) {
            throw NoThrException{};
        }
        // compute the match iff |query| is equal to the width of the panel
        if (query.size() != this->panel->w) {
            throw NotEqualLengthException{};
        }
        // if required extend with the phi support struct (iff not already
        // extended)
        if (extend_matches && !this->is_extended) {
            this->extend();
        }
        // initialize matching statistics
        ms ms(query.size());

        // algorithm begin from the last row of the first column
        // so we obtain the prefix array value (from the samples), the run index
        // and the relative symbol
        auto curr_pos = static_cast<unsigned int>(this->cols[0].sample_end[
                this->cols[0].sample_end.size() - 1]);
        auto curr_index = curr_pos;
        unsigned int curr_run = index_to_run(curr_index, 0);
        char symbol = get_next_char(this->cols[0].zero_first, curr_run);
        // iterate over every query's symbol/column index
        // iterate over every query's symbol/column index
        for (unsigned int i = 0; i < query.size(); i++) {
            //std::cout << "processed " << i << "\r";
            if (verbose) {
                std::cout << "at " << i << ": " << curr_run << " "
                          << this->cols[i].t[curr_run] << "\n";
                std::cout << curr_index << " " << curr_run << " " << curr_pos
                          << " "
                          << symbol << "\n";
            }
            // a lot of cases:
            // - if the pointer in the RLPBWT match the symbol in the query
            //   we proceed simply using lf-mapping (if we are not at the end)
            // - if the pointer in the RLPBWT mismatch the symbol in the query
            //   and we have only that symbol in the column we restart from  the
            //   next column at last index whit the relative prefix array value
            // - otherwise we proceeed using thresholds to select the best
            //   symbol, between the previous and next good symbol (if the
            //   exists), to jump
            if (query[i] == symbol) {
                if (verbose) {
                    std::cout << "match:\n";
                }
                // report in matching statistics row vector
                ms.row[i] = curr_pos;
                // update index, run, symbol if we are not at the end
                if (i != query.size() - 1) {
                    curr_index = lf(i, curr_index, query[i]);
                    curr_run = index_to_run(curr_index, i + 1);
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
                // get threshold
                auto thr = this->cols[i].t[curr_run];

                if (this->cols[i].sample_beg.size() == 1) {
                    if (verbose) {
                        std::cout << "complete mismatch\n";
                    }
                    // report in matching statistics row vector using panel
                    // height as sentinel
                    ms.row[i] = this->panel->h;
                    // update index, run, symbol (as explained before) if we are
                    // not at the end
                    if (i != query.size() - 1) {
                        curr_pos = static_cast<unsigned int>(this->cols[i +
                                                                        1].sample_end[
                                this->cols[i + 1].sample_end.size() - 1]);
                        curr_index = this->panel->h - 1;
                        curr_run = index_to_run(curr_index, i + 1);
                        symbol = get_next_char(this->cols[i + 1].zero_first,
                                               curr_run);
                        if (verbose) {
                            std::cout << "update: " << curr_index << " "
                                      << curr_pos
                                      << " "
                                      << symbol << "\n";
                        }
                    }
                } else if (curr_run != 0 &&
                           ((curr_index < thr) ||
                            curr_run == this->cols[i].sample_beg.size() - 1)) {
                    // if we are above the threshold we go up (if we are not in
                    // the first run). We also go up if we are in the last run
                    if (verbose) {
                        std::cout << "mismatch_up: ";
                    }
                    curr_index = this->cols[i].p[curr_run] - 1;
                    curr_pos = static_cast<unsigned int>(this->cols[i].sample_end[
                            curr_run - 1]);
                    if (verbose) {
                        std::cout << "update: " << curr_index << " " << curr_pos
                                  << " "
                                  << symbol << "\n";
                    }
                    // report in matching statistics row vector
                    ms.row[i] = curr_pos;
                    // update index, run, symbol if we are not at the end
                    if (i != query.size() - 1) {
                        curr_index = lf(i, curr_index, query[i]);
                        curr_run = index_to_run(curr_index, i + 1);
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
                    // we are below threshold  so we go down
                    if (verbose) {
                        std::cout << "mismatch_down: ";
                    }
                    curr_index = this->cols[i].p[curr_run + 1];
                    curr_pos = static_cast<unsigned int>(this->cols[i].sample_beg[
                            curr_run + 1]);

                    // report in matching statistics row vector
                    ms.row[i] = curr_pos;
                    if (verbose) {
                        std::cout << "update: " << curr_index << " " << curr_pos
                                  << " "
                                  << symbol << "\n";
                    }
                    // update index, run, symbol if we are not at the end
                    if (i != query.size() - 1) {
                        curr_index = lf(i, curr_index, query[i]);
                        curr_run = index_to_run(curr_index, i + 1);
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
        
        // compute the len vector using random access on the panel, proceeding
        // from left to right
	for (unsigned int i = 0; i < ms.len.size(); i++) {
            if (ms.row[i] == this->panel->h) {
                // if we have the sentinel in row vector than the length is 0
                ms.len[i] = 0;
            } else if (i != 0 && ms.row[i] == ms.row[i - 1] &&
                       ms.len[i - 1] != 0) {
                ms.len[i] = ms.len[i - 1] + 1;
            } else {
                int tmp_index = (int) i;
                unsigned int len = 0;
                while (tmp_index >= 0 &&
                       query[tmp_index] == this->panel->getElem(ms.row[i],
                                                                tmp_index)) {
                    tmp_index--;
                    len++;
                }
                ms.len[i] = len;
            }
        }

        // initialize struct for matches
        ms_matches ms_matches;
        // save every match from matching statistics (when we have a "peak" in
        // ms len vector)
        for (unsigned int i = 0; i < ms.len.size(); i++) {
            if ((i != ms.len.size() - 1 && ms.len[i] > 0 &&
                 ms.len[i] >= ms.len[i + 1]) ||
                (i == ms.len.size() - 1 && ms.len[i] != 0)) {
                ms_matches.basic_matches.emplace_back(ms.row[i], ms.len[i], i);
            }
        }
        // compute every row that are matching if required
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
            std::string query;
            if (out_match.is_open()) {
                for (unsigned int i = 0; i < queries_panel[0].size(); i++) {
                    if (verbose) {
                        std::cout << i << ": \n";
                    }
                    for (auto &j: queries_panel) {
                        query.push_back(j[i]);
                    }
                    //std::cout << i << "\n";
                    /*if (i == 167) {
                        verbose = true;
                    }*/
                    ms_matches matches;
                    matches = this->match_lce(query, extend_matches, verbose);
                    /*
                    if (i == 167) {
                        exit(-1);
                    }
                     */
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

    template<typename U = ra_t>
    std::enable_if_t<sizeof(U) && (!std::is_same<ra_t, panel_ra>::value),
            void>
    match_tsv_conc_lce(const char *filename, const char *out,
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
            std::string query;
            std::vector<std::string> queries;
            if (out_match.is_open()) {
                for (unsigned int i = 0; i < queries_panel[0].size(); i++) {
                    if (verbose) {
                        std::cout << i << ": \n";
                    }
                    for (auto &j: queries_panel) {
                        query.push_back(j[i]);
                    }
                    queries.push_back(query);
                    query.clear();
                }

                auto n_queries = queries.size();
                std::vector<ms_matches> matches_vec(n_queries);
#pragma omp parallel for default(none) shared(queries, matches_vec, n_queries, extend_matches, verbose)
                for (unsigned int i = 0; i < n_queries; i++) {
                    //std::cout << i << "\n";
                    matches_vec[i] = this->match_lce(queries[i], extend_matches,
                                                     verbose);
                }

                if (extend_matches) {
                    for (unsigned int i = 0; i < queries.size(); i++) {
                        if (!matches_vec[i].haplos.empty()) {
                            for (unsigned int j = 0;
                                 j < matches_vec[i].basic_matches.size(); j++) {
                                auto len = std::get<1>(
                                        matches_vec[i].basic_matches[j]);
                                auto end = std::get<2>(
                                        matches_vec[i].basic_matches[j]);
                                for (unsigned int k = 0;
                                     k < matches_vec[i].haplos[j].size(); k++) {
                                    out_match << "MATCH\t" << i << "\t" <<
                                              matches_vec[i].haplos[j][k]
                                              << "\t"
                                              << end - (len - 1) << "\t" << end
                                              << "\t" << len << "\n";
                                }
                            }
                        }
                    }
                } else {
                    for (unsigned int i = 0; i < queries.size(); i++) {
                        for (unsigned int j = 0;
                             j < matches_vec[i].basic_matches.size(); j++) {
                            auto len = std::get<1>(
                                    matches_vec[i].basic_matches[j]);
                            auto pos = std::get<0>(
                                    matches_vec[i].basic_matches[j]);
                            auto end = std::get<2>(
                                    matches_vec[i].basic_matches[j]);
                            out_match << "MATCH\t" << i << "\t" << pos << "\t"
                                      << end - (len - 1) << "\t" << end
                                      << "\t" << len << "\n";


                        }
                    }
                    /*for (unsigned int i = 0; i < queries.size(); i++) {
                        if (verbose) {
                            std::cout << i << ": ";
                        }
                        out_match << i << ": ";

                        if (verbose) {
                            std::cout << matches_vec[i];
                        }
                        out_match << matches_vec[i];

                        if (verbose) {
                            std::cout << "\n";
                        }
                        out_match << "\n";
                    }*/
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
            std::string query;
            if (out_match.is_open()) {
                for (unsigned int i = 0; i < queries_panel[0].size(); i++) {
                    for (auto &j: queries_panel) {
                        query.push_back(j[i]);
                    }
                    /*if (i == 167) {
                        verbose = true;
                    }*/

                    ms_matches matches;

                    matches = this->match_thr(query, extend_matches, verbose);
                    /*
                    if (i == 167) {
                        exit(-1);
                    }
                     */
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

    void
    match_tsv_conc_thr(const char *filename, const char *out,
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
            std::string query;
            std::vector<std::string> queries;
            if (out_match.is_open()) {
                for (unsigned int i = 0; i < queries_panel[0].size(); i++) {
                    if (verbose) {
                        std::cout << i << ": \n";
                    }
                    for (auto &j: queries_panel) {
                        query.push_back(j[i]);
                    }
                    queries.push_back(query);
                    query.clear();
                }

                auto n_queries = queries.size();
                std::vector<ms_matches> matches_vec(n_queries);
#pragma omp parallel for default(none) shared(queries, matches_vec, n_queries, extend_matches, verbose)
                for (unsigned int i = 0; i < n_queries; i++) {
                    //std::cout << i << "\n";
                    matches_vec[i] = this->match_thr(queries[i], extend_matches,
                                                     verbose);
                }
                if (extend_matches) {
                    for (unsigned int i = 0; i < queries.size(); i++) {
                        if (!matches_vec[i].haplos.empty()) {
                            for (unsigned int j = 0;
                                 j < matches_vec[i].basic_matches.size(); j++) {
                                auto len = std::get<1>(
                                        matches_vec[i].basic_matches[j]);
                                auto end = std::get<2>(
                                        matches_vec[i].basic_matches[j]);
                                for (unsigned int k = 0;
                                     k < matches_vec[i].haplos[j].size(); k++) {
                                    out_match << "MATCH\t" << i << "\t" <<
                                              matches_vec[i].haplos[j][k]
                                              << "\t"
                                              << end - (len - 1) << "\t" << end
                                              << "\t" << len << "\n";
                                }
                            }
                        }
                    }
                } else {
                    for (unsigned int i = 0; i < queries.size(); i++) {
                        for (unsigned int j = 0;
                             j < matches_vec[i].basic_matches.size(); j++) {
                            auto len = std::get<1>(
                                    matches_vec[i].basic_matches[j]);
                            auto pos = std::get<0>(
                                    matches_vec[i].basic_matches[j]);
                            auto end = std::get<2>(
                                    matches_vec[i].basic_matches[j]);
                            out_match << "MATCH\t" << i << "\t" << pos << "\t"
                                      << end - (len - 1) << "\t" << end
                                      << "\t" << len << "\n";


                        }
                    }
//                    for (unsigned int i = 0; i < queries.size(); i++) {
//                        if (verbose) {
//                            std::cout << i << ": ";
//                        }
//                        out_match << i << ": ";
//
//                        if (verbose) {
//                            std::cout << matches_vec[i];
//                        }
//                        out_match << matches_vec[i];
//
//                        if (verbose) {
//                            std::cout << "\n";
//                        }
//                        out_match << "\n";
//                    }
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
     * function to get the total number of runs in the RLPBWT
     * @return total number of run
     */
    unsigned int get_run_number() {
        unsigned int count_run = 0;
        for (unsigned int i = 0; i < this->cols.size(); ++i) {
            count_run += cols[i].sample_beg.size();
        }
        return count_run;
    }

    /**
     * function to print in runs.txt the number of run in every
     * column
    */
    void get_run_col() {
        unsigned int count_run = 0;
        std::ofstream myfile;
        myfile.open("runs.txt");
        unsigned int min = 1;
        unsigned int max = 1;
        for (unsigned int i = 0; i < this->cols.size(); ++i) {
            auto x = cols[i].sample_beg.size();
            myfile << x << " ";
            if (x < min) {
                min = x;
            }
            if (x > max) {
                max = x;
            }
        }
        myfile << "\n" << min << " " << max;
    }

    /**
    * function to get the total number of phi/phi_inv element() in the RLPBWT
    * @return total number of run
    */
    std::pair<unsigned int, unsigned int> get_phi_number() {
        if (this->phi) {
            unsigned int count_phi = 0;
            for (unsigned int i = 0; i < this->phi->phi_supp.size(); ++i) {
                count_phi += this->phi->phi_supp[i].size();
            }
            unsigned int count_phi_inv = 0;
            for (unsigned int i = 0; i < this->phi->phi_inv_supp.size(); ++i) {
                count_phi_inv += this->phi->phi_inv_supp[i].size();
            }
            return std::make_pair(count_phi, count_phi_inv);
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
        unsigned long long size_uv = 0;
        auto lp_size = sdsl::size_in_bytes(this->last_pref);
        unsigned long long size_samples = lp_size;
        size += lp_size;
        for (unsigned int i = 0; i < this->cols.size(); ++i) {
            size += this->cols[i].size_in_bytes();
            size_run += sdsl::size_in_bytes(this->cols[i].p);
            if (this->is_thr_enabled) {
                size_thr += sdsl::size_in_bytes(this->cols[i].t);
            }
            size_uv += sdsl::size_in_bytes(this->cols[i].uv);
            size_samples += sdsl::size_in_bytes(this->cols[i].sample_beg) +
                            sdsl::size_in_bytes(this->cols[i].sample_end);
        }
        if (verbose) {
            std::cout << "run: " << size_run << " bytes\n";
            if (this->is_thr_enabled) {
                std::cout << "thr: " << size_thr << " bytes\n";
            }
            std::cout << "uv: " << size_uv << " bytes\n";
            std::cout << "samples: " << size_samples << " bytes\n";
            std::cout
                    << "rlpbwt (with also c values and other support variables): "
                    << size << " bytes\n";
        }
        size += (sizeof(bool) * 2);
        size += (sizeof(unsigned int) * 2);
        size += this->panel->size_in_bytes(verbose);
        if (this->is_extended) {
            auto size_phi = this->phi->size_in_bytes(verbose);
            size += size_phi;
            //std::cout << "phi support: " << size_phi << " bytes\n";
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
        double size_uv = 0;
        auto lp_size = sdsl::size_in_mega_bytes(this->last_pref);
        double size_samples = lp_size;
        size += lp_size;
        for (unsigned int i = 0; i < this->cols.size(); ++i) {
            size += this->cols[i].size_in_mega_bytes();
            size_run += sdsl::size_in_mega_bytes(this->cols[i].p);
            if (this->is_thr_enabled) {
                size_thr += sdsl::size_in_mega_bytes(this->cols[i].t);
            }
            size_uv += sdsl::size_in_mega_bytes(this->cols[i].uv);
            size_samples += sdsl::size_in_mega_bytes(this->cols[i].sample_beg) +
                            sdsl::size_in_mega_bytes(this->cols[i].sample_end);
        }
        if (verbose) {
            std::cout << "run: " << size_run << " megabytes\n";
            if (this->is_thr_enabled) {
                std::cout << "thr: " << size_thr << " megabytes\n";
            }
            std::cout << "uv: " << size_uv << " megabytes\n";
            std::cout << "samples: " << size_samples << " megabytes\n";
            std::cout
                    << "rlpbwt (with also c values and other support variables): "
                    << size << " megabytes\n";
        }
        size += (sizeof(bool) * 2 * to_mega);
        size += (sizeof(unsigned int) * 2 * to_mega);
        size += this->panel->size_in_mega_bytes(verbose);
        if (this->is_extended) {
            auto size_phi = this->phi->size_in_mega_bytes(verbose);
            size += size_phi;
            //std::cout << "phi support: " << size_phi << " megabytes\n";
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
        out.write((char *) &this->width, sizeof(this->width));
        written_bytes += sizeof(this->width);
        out.write((char *) &this->height, sizeof(this->height));
        written_bytes += sizeof(this->height);
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
        in.read((char *) &this->width, sizeof(this->width));
        in.read((char *) &this->height, sizeof(this->height));
        in.read((char *) &this->is_thr_enabled, sizeof(this->is_thr_enabled));
        in.read((char *) &this->is_extended, sizeof(this->is_extended));
        this->cols = std::vector<column_naive_ms>(this->width + 1);
        if constexpr (!std::is_same_v<ra_t, panel_ra>) {
            if (std::string(slp_filename).empty()) {
                throw SlpNotFoundException{};
            }
        }
        if constexpr (!std::is_same_v<ra_t, panel_ra>) {
            this->panel = new ra_t();
            this->panel->load(in, slp_filename);
        } else {
            this->panel = new ra_t();
            this->panel->load(in);
        }

        for (unsigned int i = 0; i <= this->panel->w; i++) {
            this->cols[i].load(in);
        }
        this->last_pref.load(in);
        if (this->is_extended) {
            this->phi = new phi_support<ra_t>();
            this->phi->load(in);
        }
        if (this->cols[0].t.size() > 0) {
            this->is_thr_enabled = true;
        }
    }
};

#endif //RLPBWT_RLPBWT_NAIVE_MS_H
