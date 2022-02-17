//
// Created by dlcgold on 10/02/22.
//

#ifndef RLPBWT_RLPBWT_RA_H
#define RLPBWT_RLPBWT_RA_H

#include <vector>
#include <omp.h>
#include "column_thr.h"
#include "exceptions.h"
#include "panel_ra.h"
#include "slp_panel_ra.h"

template<typename ra_t>
class rlpbwt_ra {
public:
    std::vector<column_thr> cols;
    // panel saved as requested
    ra_t *panelbv;

private:
    static column_thr
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
        // variable to compute curr lcs in order to eventually compute thresholds
        // unsigned int lcs = 0;
        // bool to check first symbol fo the column
        bool start = true;

        std::vector<std::pair<unsigned int, unsigned int>> rows;
        unsigned int tmp_beg = 0;
        unsigned int thr_tmp = 0;
        unsigned int lcs = 0;
        // update start and "c" value
        for (unsigned int i = 0; i < height; i++) {
            if (i == 0 && column[pref[i]] == '1') {
                start = false;
            }
            if (column[i] == '0') {
                count0++;
            }
        }

        // initialize the three bitvectors
        sdsl::bit_vector runvec(height + 1, 0);
        sdsl::bit_vector thr(height, 0);
        sdsl::bit_vector zerovec(count0, 0);
        sdsl::bit_vector onevec(height - count0, 0);

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


            // stuff for thresholds
            if ((i == 0) || (column[pref[i]] != column[pref[i - 1]])) {
                thr_tmp = i;
                lcs = div[i];
            }
            if (div[i] < lcs) {
                thr_tmp = i;
                lcs = div[i];
            }

            if ((i == height - 1) || (column[pref[i]] != column[pref[i + 1]])) {
                // 1 in bitvectors for runs et every end of a run
                runvec[i] = true;
                // 1 in bitvectors for thresholds index
                thr[thr_tmp] = true;

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
        return {start, count0, runvec, zerovec, onevec, thr, sample_beg,
                sample_end};
    }

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


public:
    rlpbwt_ra() = default;

    rlpbwt_ra(const char *filename, bool verbose,
              const char *slp_filename = "") {

        std::ifstream input_matrix(filename);
        if constexpr (!std::is_same_v<ra_t, panel_ra>) {
            if (std::string(slp_filename).empty()) {
                throw SlpNotFoundException{};
            }
        }
        if (input_matrix.is_open()) {
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
            this->cols = std::vector<column_thr>(tmp_width);
            std::vector<unsigned int> pref(tmp_height);
            sdsl::int_vector<> div(tmp_height);
            for (unsigned int i = 0; i < tmp_height; i++) {
                pref[i] = i;
                div[i] = 0;
            }
            if constexpr (std::is_same_v<ra_t, panel_ra>) {
                this->panelbv = new ra_t(tmp_height, tmp_width);
            } else if constexpr (std::is_same_v<ra_t, slp_panel_ra>) {
                this->panelbv = new ra_t(slp_filename, tmp_height, tmp_width);
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
                auto col = rlpbwt_ra::build_column(new_column, pref, div);
                if constexpr (std::is_same_v<ra_t, panel_ra>) {
                    for (unsigned int k = 0; k < new_column.size(); k++) {
                        if (new_column[k] != '0') {
                            this->panelbv->panel[count][k] = true;
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
                rlpbwt_ra::update(new_column, pref, div);
                last_col = new_column;
                count++;
            }
            input_matrix.close();
        } else {
            throw FileNotFoundException{};
        }
    }

    std::pair<unsigned int, unsigned int>
    lce(unsigned int col, unsigned int curr, unsigned int prev,
        bool verbose = false) {
        if (col == 0) {
            return std::make_pair(prev, 0);
        }
        unsigned int rev_col = (this->panelbv->w - 1) - col + 1;

        unsigned int pos_curr = rev_col + ((this->panelbv->w) * curr);
        unsigned int pos_prev = rev_col + ((this->panelbv->w) * prev);
        if (verbose) {
            std::cout << "at " << rev_col << ": " << pos_curr << ", "
                      << pos_prev << "\n";
        }
        unsigned int lcp_prev = lceToR(this->panelbv->panel, pos_prev,
                                       pos_curr);
        if (lcp_prev >= col) {
            lcp_prev = col;
        }
        return std::make_pair(prev, lcp_prev);
    }

    std::pair<unsigned int, unsigned int>
    lce_pair(unsigned int col, unsigned int curr, unsigned int prev,
             unsigned int next, bool verbose = false) {
        if (col == 0) {
            return std::make_pair(next, 0);
        }
        unsigned int rev_col = (this->panelbv->w - 1) - col + 1;

        unsigned int pos_curr = rev_col + ((this->panelbv->w) * curr);
        unsigned int pos_prev = rev_col + ((this->panelbv->w) * prev);
        unsigned int pos_next = rev_col + ((this->panelbv->w) * next);
        if (verbose) {
            std::cout << "at " << rev_col << ": " << pos_curr << ", "
                      << pos_prev
                      << ", " << pos_next << "\n";
        }
        unsigned int lcp_prev = lceToR(this->panelbv->panel, pos_prev,
                                       pos_curr);
        if (lcp_prev >= col) {
            lcp_prev = col;
        }
        unsigned int lcp_next = lceToR(this->panelbv->panel, pos_curr,
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

    std::vector<std::pair<unsigned int, unsigned int>>
    match_thr(const std::string &query, bool verbose = false) {
        if (query.size() != this->panelbv->w) {
            throw NotEqualLengthException{};
        }
        std::vector<unsigned int> ms_row(query.size(), 0);
        std::vector<unsigned int> ms_len(query.size(), 0);
        //auto curr_prefs = this->cols[0].rows[this->cols[0].rows.size() - 1];
        //auto curr_pos = curr_prefs.second;
        auto curr_pos = static_cast<unsigned int>(this->cols[0].sample_end[
                this->cols[0].sample_end.size() - 1]);
        auto curr_index = curr_pos;
        unsigned int curr_run = this->cols[0].rank_runs(curr_index);
        char symbol = get_next_char(this->cols[0].zero_first, curr_run);

        for (unsigned int i = 0; i < query.size(); i++) {
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
                // save in matching statistics pos vector
                ms_row[i] = curr_pos;
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
                    ms_row[i] = this->panelbv->h;
                    //curr_index = curr_pos;

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
                    ms_row[i] = curr_pos;
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
                    ms_row[i] = curr_pos;
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
        if (verbose) {
            std::cout << "ind:\t";
            for (unsigned int i = 0; i < ms_row.size(); i++) {
                std::cout << i << "\t";
            }
            std::cout << "\npos:\t";
            for (auto e: ms_row) {
                std::cout << e << "\t";
            }
            std::cout << "\nlen:\t";
        }
        int tmp_index = 0;

        for (int i = (int) ms_row.size() - 1; i >= 0; i--) {
            if (ms_row[i] == this->panelbv->h) {
                ms_len[i] = 0;
                continue;
            }
            tmp_index = i;
            while (tmp_index >= 0 &&
                   query[tmp_index] == panelbv->getElem(ms_row[i], tmp_index)) {
                tmp_index--;
            }
            ms_len[i] = i - tmp_index;
        }

        std::vector<std::pair<unsigned int, unsigned int>> ms_match;
        for (unsigned int i = 0; i < ms_len.size(); i++) {
            if (verbose) {
                std::cout << ms_len[i] << "\t";
            }
            // TODO check if for last match
            /*if (i > 0 && ms_len[i] >= ms_len[i - 1] &&
                (i == ms_len.size() - 1 || ms_len[i] >= ms_len[i + 1])) {
                ms_match.emplace_back(i, ms_len[i]);
            }*/

            if ((ms_len[i] != 0 && ms_len[i] > ms_len[i + 1]) ||
                (i == ms_len.size() - 1 && ms_len[i] != 0)) {
                ms_match.emplace_back(i, ms_len[i]);
            } else if (ms_len[i] != 0 && i < ms_len.size() - 1 &&
                       ms_len[i] == ms_len[i + 1]) {
                unsigned int pos = 0;
                for (unsigned int j = i + 1; j < ms_len.size(); j++) {
                    if (ms_len[j] > ms_len[j + 1]) {
                        pos = j + 1;
                        break;
                    } else {
                        pos = i + 1;
                        break;
                    }
                }
                if (i + 1 != pos) {
                    if (verbose) {
                        for (unsigned int j = i + 1; j <= pos; j++) {
                            std::cout << ms_len[j] << "\t";
                        }
                    }
                    for (unsigned int j = i; j < pos; j++) {
                        ms_match.emplace_back(j, ms_len[j]);
                    }
                    i = pos;
                }
                if (pos == ms_len.size() - 1) {
                    break;
                }
            }
        }
        if (verbose) {
            std::cout << "\nmatches:\t";
            for (auto e: ms_match) {
                std::cout << "(col: " << e.first << ", len: " << e.second
                          << ")\t";
            }
            std::cout << "\n";
        }
        return ms_match;
    }

    std::vector<std::pair<unsigned int, unsigned int>>
    match_lce(const std::string &query, bool verbose = false) {
        if (query.size() != this->panelbv->w) {
            throw NotEqualLengthException{};
        }
        std::vector<unsigned int> ms_row(query.size(), 0);
        std::vector<unsigned int> ms_len(query.size(), 0);
        auto curr_pos = static_cast<unsigned int>(this->cols[0].sample_end[
                this->cols[0].sample_end.size() - 1]);
        auto curr_index = curr_pos;
        unsigned int curr_run = this->cols[0].rank_runs(curr_index);
        char symbol = get_next_char(this->cols[0].zero_first, curr_run);

        for (unsigned int i = 0; i < query.size(); i++) {
            std::cout << "\t" << i << "\n";
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
                // save in matching statistics pos vector
                ms_row[i] = curr_pos;
                if (i != 0) {
                    ms_len[i] = ms_len[i - 1] + 1;
                } else {
                    ms_len[i] = 1;
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
                    ms_row[i] = this->panelbv->h;
                    ms_len[i] = 0;
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
                        ms_row[i] = prev_pos;
                        if (i == 0) {
                            ms_len[i] = 1;
                        } else {
                            ms_len[i] =
                                    std::min(ms_len[i - 1], lce_value.second) +
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
                        ms_row[i] = next_pos;
                        if (i == 0) {
                            ms_len[i] = 1;
                        } else {
                            ms_len[i] =
                                    std::min(ms_len[i - 1], lce_value.second) +
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
                            ms_row[i] = curr_pos;
                            if (i == 0) {
                                ms_len[i] = 1;
                            } else {
                                ms_len[i] =
                                        std::min(ms_len[i - 1], lce.second) + 1;
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
                            ms_row[i] = curr_pos;
                            if (i == 0) {
                                ms_len[i] = 1;
                            } else {
                                ms_len[i] =
                                        std::min(ms_len[i - 1], lce.second) + 1;
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
        if (verbose) {
            std::cout << "ind:\t";
            for (unsigned int i = 0; i < ms_row.size(); i++) {
                std::cout << i << "\t";
            }
            std::cout << "\npos:\t";
            for (auto e: ms_row) {
                std::cout << e << "\t";
            }
            std::cout << "\nlen:\t";
            for (auto e: ms_len) {
                std::cout << e << "\t";
            }
        }
        //int tmp_index = 0;

        std::vector<std::pair<unsigned int, unsigned int>> ms_match;
        for (unsigned int i = 0; i < ms_len.size(); i++) {
            // TODO check if for last match

            if ((ms_len[i] != 0 && ms_len[i] > ms_len[i + 1]) ||
                (i == ms_len.size() - 1 && ms_len[i] != 0)) {
                ms_match.emplace_back(i, ms_len[i]);
            } else if (ms_len[i] != 0 && i < ms_len.size() - 1 &&
                       ms_len[i] == ms_len[i + 1]) {
                unsigned int pos = 0;
                for (unsigned int j = i + 1; j < ms_len.size(); j++) {
                    if (ms_len[j] > ms_len[j + 1]) {
                        pos = j + 1;
                        break;
                    } else {
                        pos = i + 1;
                        break;
                    }
                }
                if (i + 1 != pos) {
                    for (unsigned int j = i; j < pos; j++) {
                        ms_match.emplace_back(j, ms_len[j]);
                    }
                    i = pos;
                }
                if (pos == ms_len.size() - 1) {
                    break;
                }
            }
        }
        if (verbose) {
            std::cout << "\nmatches:\t";
            for (auto e: ms_match) {
                std::cout << "(col: " << e.first << ", len: " << e.second
                          << ")\t";
            }
            std::cout << "\n";
        }
        return ms_match;
    }

    void
    match_tsv(const char *filename, const char *out, bool verbose = false) {
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
                    auto matches = this->match_thr(query, verbose);
                    if (verbose) {
                        std::cout << i << ": ";
                    }
                    out_match << i << ": ";
                    for (auto m: matches) {
                        if (verbose) {
                            std::cout << "(col: " << m.first << ", len: "
                                      << m.second << ") ";
                        }
                        out_match << "(col: " << m.first << ", len: "
                                  << m.second << ") ";
                    }
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
    match_tsv_tr(const char *filename, const char *out, bool verbose = false) {
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
            if (out_match.is_open()) {
                for (unsigned int i = 0; i < queries_panel.size(); i++) {
                    auto matches = this->match_thr(queries_panel[i], verbose);
                    if (verbose) {
                        std::cout << i << ": ";
                    }
                    out_match << i << ": ";
                    for (auto m: matches) {
                        if (verbose) {
                            std::cout << "(col: " << m.first << ", len: "
                                      << m.second << ") ";
                        }
                        out_match << "(col: " << m.first << ", len: "
                                  << m.second << ") ";
                    }
                    if (verbose) {
                        std::cout << "\n";
                    }
                    out_match << "\n";
                }
                out_match.close();
            } else {
                throw FileNotFoundException{};
            }

        } else {
            throw FileNotFoundException{};
        }
    }

    size_t serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr,
                     const std::string &name = "") {
        sdsl::structure_tree_node *child =
                sdsl::structure_tree::add_child(v, name,
                                                sdsl::util::class_name(
                                                        *this));
        size_t written_bytes = 0;
        written_bytes += this->panelbv.serialize(out, child, "panel");

        for (unsigned int i = 0; i < this->cols.size(); i++) {
            std::string label = "col_" + std::to_string(i);
            written_bytes += this->cols[i].serialize(out, child, label);
        }

        sdsl::structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    void load(std::istream &in) {
        this->panelbv.load(in);
        for (unsigned int i = 0; i < this->panelbv.w; i++) {
            auto c = new column_thr();
            c->load(in);
            this->cols.emplace_back(*c);
        }
    }
};


#endif //RLPBWT_RLPBWT_RA_H