//
// Created by dlcgold on 02/02/22.
//

#include <list>
#include "../include/rlpbwt_thr.h"

rlpbwt_thr::rlpbwt_thr(const char *filename, bool vcf, bool verbose) {
    std::ifstream input_matrix(filename);
    if (input_matrix.is_open()) {
        std::string new_column;
        getline(input_matrix, new_column);
        new_column.erase(
                std::remove(new_column.begin(), new_column.end(), ' '),
                new_column.end());
        const unsigned int tmp_height = new_column.size();
        unsigned int tmp_width = std::count(
                std::istreambuf_iterator<char>(input_matrix),
                std::istreambuf_iterator<char>(), '\n') + 1;
        std::vector<column_thr> tmp_cols(tmp_width);
        this->cols = std::vector<column_thr>(tmp_width + 1);
        input_matrix.clear();
        input_matrix.seekg(0, std::ios::beg);
        std::vector<unsigned int> pref(tmp_height);
        sdsl::int_vector<> div(tmp_height);
        for (unsigned int i = 0; i < tmp_height; i++) {
            pref[i] = i;
            div[i] = 0;
        }
        unsigned int count = 0;
        std::string last_col;
        while (getline(input_matrix, new_column)) {
            if (verbose) {
                std::cout << "\nnew_column " << count << "\n";
                std::cout << new_column << "\n" << this->cols[count].runs
                          << "\n"
                          << this->cols[count].u << "\n"
                          << this->cols[count].v
                          << "\n-------------------------------\n";
            }
            new_column.erase(
                    std::remove(new_column.begin(), new_column.end(), ' '),
                    new_column.end());
            auto col = rlpbwt_thr::build_column(new_column, pref, div);
            this->panel.push_back(new_column);
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
            rlpbwt_thr::update(new_column, pref, div);
            count++;
            last_col = new_column;
        }

        auto col = rlpbwt_thr::build_column(last_col, pref, div);
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

        input_matrix.close();
    } else {
        throw FileNotFoundException{};
    }
}

column_thr
rlpbwt_thr::build_column(std::string &column, std::vector<unsigned int> &pref,
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
    return {start, count0, runvec, zerovec, onevec, thr, rows};
}

void rlpbwt_thr::update(std::string &column, std::vector<unsigned int> &pref,
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
    new_div[count0] = 0;
    div = new_div;
    pref = new_pref;
}

unsigned int
rlpbwt_thr::lf(unsigned int col_index, unsigned int index, char symbol,
               bool verbose) const {
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

// TODO this function is written very bad
std::pair<unsigned int, unsigned int>
rlpbwt_thr::uvtrick(unsigned int col_index, unsigned int index) const {
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
                        index - (this->cols[col_index].select_runs(run) + 1)};
            } else {
                return {index - (this->cols[col_index].select_runs(run) + 1),
                        this->cols[col_index].select_runs(run) + 1};
            }
        }
    }
}

void rlpbwt_thr::match_thr(const std::string &query, bool verbose) {
    std::vector<unsigned int> ms_row(query.size(), 0);
    std::vector<unsigned int> ms_len(query.size(), 0);
    auto curr_prefs = this->cols[0].rows[this->cols[0].rows.size() - 1];
    auto curr_pos = curr_prefs.second;
    auto curr_index = curr_pos;
    unsigned int curr_run = this->cols[0].rank_runs(curr_index);
    char symbol = get_next_char(this->cols[0].zero_first, curr_run);
    for (unsigned int i = 0; i < query.size(); i++) {

        std::cout << i << ": " << curr_run << " "
                  << this->cols[i].rank_thr(curr_index) << "\n";
        std::cout << curr_index << " " << curr_run << " " << curr_pos << " "
                  << symbol << "\n";

        if (query[i] == symbol) {
            std::cout << "match: ";
            // save in matching statistics pos vector
            ms_row[i] = curr_pos;
            // update index, run, symbol
            curr_index = lf(i, curr_index, query[i]);
            curr_run = this->cols[i + 1].rank_runs(curr_index);
            symbol = get_next_char(this->cols[i + 1].zero_first, curr_run);
            std::cout << "new: " << curr_index << " " << curr_run << " "
                      << curr_pos << " "
                      << symbol << "\n";
        } else {
            //auto run = this->cols[i].rank_runs(curr_index);
            auto thr = this->cols[i].rank_thr(curr_index);
            bool single = false;
            if (this->cols[i].rows[curr_run].first ==
                this->cols[i].rows[curr_run].second) {
                single = true;
            }
            if ((curr_run != 0 && !single && curr_run == thr) ||
                curr_run == this->cols[i].rows.size()) {
                std::cout << "mismatch_up: ";
                // threshold below index so we go up
                curr_index = (this->cols[i].select_runs(curr_run) + 1) - 1;
                curr_prefs = this->cols[i].rows[curr_run - 1];
                curr_pos = curr_prefs.second;
                std::cout << "update: " << curr_index << " " << curr_pos << " "
                          << symbol << "\n";
                ms_row[i] = curr_pos;
                curr_index = lf(i, curr_index, query[i]);
                curr_run = this->cols[i + 1].rank_runs(curr_index);
                symbol = get_next_char(this->cols[i + 1].zero_first, curr_run);
                std::cout << "new: " << curr_index << " " << curr_run << " "
                          << curr_pos << " " << symbol << "\n";

            } else {
                std::cout << "mismatch_up: ";
                // threshold above index so we go down
                curr_index = (this->cols[i].select_runs(curr_run + 1) + 1);
                curr_prefs = this->cols[i].rows[curr_run + 1];
                curr_pos = curr_prefs.first;
                ms_row[i] = curr_pos;
                std::cout << "update: " << curr_index << " " << curr_pos << " "
                          << symbol << "\n";
                curr_index = lf(i, curr_index, query[i]);
                curr_run = this->cols[i + 1].rank_runs(curr_index);
                symbol = get_next_char(this->cols[i + 1].zero_first, curr_run);
                std::cout << "new: " << curr_index << " " << curr_run << " "
                          << curr_pos << " " << symbol << "\n";

            }
        }
    }
    std::cout << "ind:\t";
    for (unsigned int i = 0; i < ms_row.size(); i++) {
        std::cout << i << "\t";
    }
    std::cout << "\npos:\t";
    for (auto e: ms_row) {
        std::cout << e << "\t";
    }
    std::cout << "\nlen:\t";
    int tmp_index = 0;
    for (int i = (int) ms_row.size() - 1; i >= 0; i--) {
        /*if (i != (int) ms_row.size() - 1 && ms_row[i] == ms_row[i + 1]) {
            if (ms_len[i + 1] != 0) {
                ms_len[i] = ms_len[i + 1] - 1;
            } else {
                ms_len[i] = 0;
            }
        } else {*/
        tmp_index = i;
        while (tmp_index > 0 &&
               query[tmp_index] == panel[tmp_index][ms_row[i]]) {
            tmp_index--;
        }
        if (tmp_index == 0) {
            tmp_index--;
        }
        ms_len[i] = i - tmp_index;
    }
    //}
    for (auto e: ms_len) {
        std::cout << e << "\t";
    }
    std::cout << "\n";

}
