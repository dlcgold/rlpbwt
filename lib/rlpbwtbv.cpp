//
// Created by dlcgold on 28/01/22.
//

#include "../include/rlpbwtbv.h"

rlpbwtbv::rlpbwtbv(const char *filename, bool vcf, bool verbose) {
    // two cases:
    // - build from vcf file
    // - build from raw matrix file (every row is a column of the matrix)
    if (vcf) {
        // read vcf file from https://github.com/ZhiGroup/Syllable-PBWT
        std::string line;
        std::ifstream input(filename);
        if (!input.is_open()) {
            throw FileNotFoundException{};
        }
        while (std::getline(input, line)) {
            if (line[0] != '#' || line[1] != '#') {
                break;
            }
        }
        std::stringstream ss(line);
        int tmp_height = -9;
        int tmp_width = 0;
        while (getline(ss, line, '\t')) {
            tmp_height++;
        }
        tmp_height <<= 1;
        while (std::getline(input, line)) {
            tmp_width++;
        }
        input.clear();
        input.seekg(0);
        while (std::getline(input, line)) {
            if (line[0] != '#' || line[1] != '#') {
                break;
            }
        }
        ss = std::stringstream(line);
        for (int i = 0; i < 9; i++) {
            getline(ss, line, '\t');
        }
        /*std::vector<std::string> IDs(tmp_height);
        for (int i = 0; i < tmp_height; i += 2) {
            getline(ss, IDs[i], '\t');
            IDs[i + 1] = IDs[i] + "-1";
            IDs[i] += "-0";
            std::cout << IDs[i] << ", " << IDs[i+1] << "\n";
        }
        */

        // initialize prefix and divergence arrays (latter as sdsl int_vector)
        std::vector<unsigned int> pref(tmp_height);
        sdsl::int_vector<> div(tmp_height);
        for (int i = 0; i < tmp_height; i++) {
            pref[i] = i;
            div[i] = 0;
        }
        std::string new_column;
        // initialize vector for the column
        this->cols = std::vector<columnbv>(tmp_width + 1);
        int k = 0;
        for (k = 0; k < tmp_width; k++) {
            if (verbose) {
                std::cout << "\nnew_column " << k << "\n";
                for (auto e: pref) {
                    std::cout << e << " ";
                }
                std::cout << "\n";
                for (auto e: div) {
                    std::cout << e << " ";
                }
                std::cout << "\n";
            }

            // read column as in https://github.com/ZhiGroup/Syllable-PBWT
            new_column.clear();
            std::getline(input, line);
            ss = std::stringstream(line);
            for (int i = 0; i < 9; i++) {
                getline(ss, line, '\t');
            }
            //int index = 0;
            while (getline(ss, line, '\t')) {
                new_column.push_back(line[0]);
                new_column.push_back(line[2]);
            }

            // create column with the bitvectors (but not rank/select) and add
            // to columns vector
            auto col = rlpbwtbv::build_column(new_column, pref, div);
            this->cols[k] = col;

            // create rank/select (here because of some problems using sdsl)
            this->cols[k].rank_runs = sdsl::sd_vector<>::rank_1_type(
                    &this->cols[k].runs);
            this->cols[k].select_runs = sdsl::sd_vector<>::select_1_type(
                    &this->cols[k].runs);
            this->cols[k].rank_u = sdsl::sd_vector<>::rank_1_type(
                    &this->cols[k].u);
            this->cols[k].select_u = sdsl::sd_vector<>::select_1_type(
                    &this->cols[k].u);
            this->cols[k].rank_v = sdsl::sd_vector<>::rank_1_type(
                    &this->cols[k].v);
            this->cols[k].select_v = sdsl::sd_vector<>::select_1_type(
                    &this->cols[k].v);
            rlpbwtbv::update(new_column, pref, div);
        }
        // build another column with last prefix and divergence arrays
        auto col = rlpbwtbv::build_column(new_column, pref, div);
        this->cols[k] = col;
        this->cols[k].rank_runs = sdsl::sd_vector<>::rank_1_type(
                &this->cols[k].runs);
        this->cols[k].select_runs = sdsl::sd_vector<>::select_1_type(
                &this->cols[k].runs);
        this->cols[k].rank_u = sdsl::sd_vector<>::rank_1_type(
                &this->cols[k].u);
        this->cols[k].select_u = sdsl::sd_vector<>::select_1_type(
                &this->cols[k].u);
        this->cols[k].rank_v = sdsl::sd_vector<>::rank_1_type(
                &this->cols[k].v);
        this->cols[k].select_v = sdsl::sd_vector<>::select_1_type(
                &this->cols[k].v);
        this->width = tmp_width;
        this->height = tmp_height;
    } else {
        // same as for vcf files but using raw matrix input files
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
            std::vector<columnbv> tmp_cols(tmp_width);
            this->cols = std::vector<columnbv>(tmp_width + 1);
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
                auto col = rlpbwtbv::build_column(new_column, pref, div);
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

                rlpbwtbv::update(new_column, pref, div);
                count++;
                last_col = new_column;
            }

            auto col = rlpbwtbv::build_column(last_col, pref, div);
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

            this->width = tmp_width;
            this->height = tmp_height;
            input_matrix.close();
        } else {
            throw FileNotFoundException{};
        }
    }
}


columnbv
rlpbwtbv::build_column(std::string &column, std::vector<unsigned int> &pref,
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

        /*
        // stuff for thresholds
        if ((i == 0) || (column[pref[i]] != column[pref[i - 1]])) {
            lcs = div[i];
        }
        if (div[i] < lcs) {
            lcs = div[i];

        }*/

        if ((i == height - 1) || (column[pref[i]] != column[pref[i + 1]])) {
            // 1 in bitvectors for runs et every end of a run
            runvec[i] = true;

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
    return {start, count0, runvec, zerovec, onevec, div};
}


// algorithm 2 from Durbin's paper
void rlpbwtbv::update(std::string &column, std::vector<unsigned int> &pref,
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

std::vector<match>
rlpbwtbv::external_match(const std::string &query, unsigned int min_len,
                         bool verbose) {
    // query allowed iff |query| is uqual to RLPBWT width
    if (query.size() != this->width) {
        throw NotEqualLengthException{};
    }
    // initialize matches vector
    std::vector<match> matches;

    // initialize all variables to support the computation
    // "curr" variables refers to "f" and "end" to "g" in Durbin's paper
    // TODO maybe some of them are useless
    unsigned int curr_run = 0;
    unsigned int end_run = 0;
    unsigned int curr_index = 0;
    unsigned int end_index = 0;
    unsigned int curr_tmp = 0;
    unsigned int end_tmp = 0;

    // current amount of rows that are matching
    unsigned int curr_len = 0;

    // current begin of a match, as "e" in Durbin's paper
    unsigned int curr_beg = 0;

    // iterate over every column
    for (unsigned i = 0; i < query.size(); i++) {
        if (verbose) {
            std::cout << "before at " << i << " from " << curr_index << " to "
                      << end_index
                      << "\n";
        }

        // compute run of the current indices using rank over the runs bv
        curr_run = this->cols[i].rank_runs(curr_index);
        end_run = this->cols[i].rank_runs(end_index);

        // update current number of rows that are matching
        curr_len = (end_index - curr_index);

        // lf-mapping of the two indices
        curr_tmp = lf(i, curr_index, query[i]);
        end_tmp = lf(i, end_index, query[i]);

        if (verbose) {
            std::cout << "middle at " << i << " from " << curr_tmp << " to "
                      << end_tmp
                      << "\n";
        }

        // two cases:
        // - we could continue to next column because match is extendable to
        //   the right
        // - we can't continue so at previous column is ending a match. In this
        //   case we have to update all the indices
        if (curr_tmp < end_tmp) {
            // first case
            if (verbose) {
                std::cout << "case 1\n";
            }
            curr_index = curr_tmp;
            end_index = end_tmp;
        } else {
            // second case
            if (verbose) {
                std::cout << "f: " << curr_tmp << ", g: " << end_tmp << "\n";
            }

            // report matches if longer than a minimum length
            // TODO this is a row solution regarding the minimum length aspect
            if (i > 0) {
                if ((i - 1) - curr_beg >= min_len) {
                    if (verbose) {
                        std::cout << "match at (" << curr_beg << ", " << i - 1
                                  << ") with " << curr_len << " haplotypes \n";
                    }
                    matches.emplace_back(curr_beg, i - 1, curr_len);
                }
            }

            // update current begin index as in Durbin using LCP array and curr
            // index, eventually using a little trick (the first case) if curr
            // is at the end of the column
            if (curr_tmp == this->cols[i + 1].lcp.size()) {
                curr_beg = i + 1;
            } else {
                curr_beg = i - this->cols[i + 1].lcp[curr_tmp];
            }
            if (verbose) {
                std::cout << "before curr beg: " << curr_beg << "\n";
            }

            // now two cases in order to update the indices:
            // - we read a 0 from query and curr index is not at the beginning
            //   of the column or curr index is at the end of the column.
            //   In this case end index won't be updated while curr index will
            //   be decreased
            // - every other situation. In this case curr index won't be updated
            //   while end index will be inreased
            if ((query[curr_beg] == '0' && curr_tmp > 0) ||
                curr_tmp == this->height) {
                if (verbose) {
                    std::cout << "begin case 2 curr run " << curr_run
                              << ", curr index " << curr_index << "\n";
                    std::cout << "currtmp: " << curr_tmp << ", endtmp: "
                              << end_tmp
                              << "\n";
                }
                // first update to curr index
                curr_tmp = end_tmp - 1;

                // we look back to update index of begin of matches
                // (iff it's not already at the minimum possible, 0)
                // we don't have both the panel and the prefix array in
                // memory, so we have to proceed reversing the permutations
                if (curr_beg >= 1) {
                    // as in Durbin we have to look at next column
                    unsigned int i_tmp = i + 1;
                    // we use a temporary index
                    unsigned int curr_rev = curr_tmp;
                    // reach (curr_beg - 1) column using the "reverse
                    // lf-mapping" and we calculate a temporary new index
                    while (i_tmp != curr_beg - 1) {
                        curr_rev = reverse_lf(i_tmp, curr_rev, false);
                        i_tmp--;
                    }

                    // we calculate the run of the temporary index
                    // and the symbol in that run
                    unsigned int tmp_run = this->cols[i_tmp].rank_runs(
                            curr_rev);
                    char curr_elem = get_next_char(this->cols[i_tmp].zero_first,
                                                   tmp_run);

                    // we virtually follow the row of the temporary index
                    // in order to update the beginning index
                    while (curr_beg > 0 && query[curr_beg - 1] == curr_elem) {
                        curr_beg -= 1;
                        if (curr_beg > 0) {
                            curr_rev = reverse_lf(curr_beg, curr_rev,
                                                  false);
                            tmp_run = this->cols[curr_beg - 1].rank_runs(
                                    curr_rev);
                            curr_elem = get_next_char(
                                    this->cols[curr_beg - 1].zero_first,
                                    tmp_run);
                        }
                    }
                }
                // after the computation of the new beginning of a match we use
                // lcp array in order to calculate how many rows are matching at
                // this point of the query, updating the curr index
                while (curr_tmp > 0 &&
                       (i + 1) - this->cols[i + 1].lcp[curr_tmp] <= curr_beg) {
                    curr_tmp--;
                }
                // update in first case is complete
                curr_index = curr_tmp;
                end_index = end_tmp;
                if (verbose) {
                    std::cout << "end case 2 curr curr index " << curr_index
                              << ", end index " << end_index << "\n";
                }
            } else {
                if (verbose) {
                    std::cout << "begin case 3 end run " << end_run
                              << ", end index " << end_index << "\n";
                    std::cout << "currtmp: " << curr_tmp << ", endtmp: "
                              << end_tmp
                              << "\n";
                }
                // first update to end index
                end_tmp = curr_tmp + 1;

                // we look back to update index of begin of matches
                // (iff it's not already at the minimum possible, 0)
                // we don't have both the panel and the prefix array in
                // memory, so we have to proceed reversing the permutations
                if (curr_beg >= 1) {
                    // as in Durbin we have to look at next column
                    unsigned int i_tmp = i + 1;
                    // we use a temporary index
                    unsigned int curr_rev = curr_tmp;
                    // reach (curr_beg - 1) column using the "reverse
                    // lf-mapping" and we calculate a temporary new index
                    while (i_tmp != curr_beg - 1) {
                        curr_rev = reverse_lf(i_tmp, curr_rev, false);
                        i_tmp--;
                    }

                    // we calculate the run of the temporary index
                    // and the symbol in that run
                    unsigned int tmp_run = this->cols[i_tmp].rank_runs(
                            curr_rev);
                    char curr_elem = get_next_char(
                            this->cols[i_tmp].zero_first,
                            tmp_run);

                    // we virtually follow the row of the temporary index
                    // in order to update the beginning index
                    while (curr_beg > 0 && query[curr_beg - 1] == curr_elem) {
                        curr_beg -= 1;
                        if (curr_beg > 0) {
                            curr_rev = reverse_lf(curr_beg, curr_rev,
                                                  false);
                            tmp_run = this->cols[curr_beg - 1].rank_runs(
                                    curr_rev);
                            curr_elem = get_next_char(
                                    this->cols[curr_beg - 1].zero_first,
                                    tmp_run);
                        }
                    }
                }

                if (verbose) {
                    std::cout << "end curr beg: " << curr_beg << "\n";
                }
                // after the computation of the new beginnoing of a match we use
                // lcp array in order to calculate how many rows are matching at
                // this point of the query, updating the end index
                while (end_tmp < this->height &&
                       (i + 1) - this->cols[i + 1].lcp[end_tmp] <= curr_beg) {
                    end_tmp++;
                }

                // update in second case is complete
                curr_index = curr_tmp;
                end_index = end_tmp;
                if (verbose) {
                    std::cout << "end case 3 curr index " << curr_index
                              << ", end index " << end_index << "\n";
                }
            }
        }
        if (verbose) {
            std::cout << "after at " << i << " from " << curr_index << " to "
                      << end_index
                      << "\n";
        }
    }
    // check matches at the end
    if (curr_index < end_index) {
        curr_len = end_index - curr_index;
        if (verbose) {
            std::cout << "match at (" << curr_beg << ", " << query.size() - 1
                      << ") with " << curr_len << " haplotypes \n";
        }
        // TODO this is a row solution regarding the minimum length aspect
        if ((query.size() - 1) - curr_beg >= min_len) {
            matches.emplace_back(curr_beg, query.size() - 1, curr_len);
        }
    }
    return matches;
}

unsigned int
rlpbwtbv::lf(unsigned int col_index, unsigned int index, char symbol,
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
rlpbwtbv::uvtrick(unsigned int col_index, unsigned int index) const {
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

unsigned int
rlpbwtbv::reverse_lf(unsigned int col_index, unsigned int index,
                     bool verbose) const {
    // by design if we try to work on first column the function return 0
    if (col_index == 0) {
        return 0;
    }
    // we extract the "c" value from the previous column
    unsigned int prev_col = col_index - 1;
    unsigned int c = this->cols[prev_col].count_0;
    // index of the run in the previous column
    unsigned int run = 0;

    // two cases:
    // - if index is less than previous "c" it means that it comes from a zero
    //   element, so we will use "u" bitvector
    // - otherwise it means that it comes from a one element, so we will use "v"
    //   bitvector
    if (index < c) {
        // computing the run using "u" bitvector according to the first
        // symbol of the column
        run = this->cols[prev_col].rank_u(index) * 2;
        if (!this->cols[prev_col].zero_first) {
            run++;
        }
        if (verbose) {
            std::cout << "run1: " << run << "\n";
        }

        // compute the index of the run
        unsigned int index_run = 0;
        if (run != 0) {
            index_run = this->cols[prev_col].select_runs(run) + 1;
        }

        // compute the previous zeros before the run to calculate the offset
        unsigned int zeros = uvtrick(prev_col, index_run).first;
        if (verbose) {
            std::cout << "case 1: " << index_run
                      << " + " << index << " z: " << zeros << "\n";
        }

        // return the index of the run plus the offset
        return index_run + (index - zeros);
    } else {
        // computing the run using "v" bitvector according to the first
        // symbol of the column
        run = this->cols[prev_col].rank_v(index - c) * 2;
        if (this->cols[prev_col].zero_first) {
            run++;
        }
        if (verbose) {
            std::cout << "run2: " << run << "\n";
        }

        // compute the index of the run
        unsigned int index_run = 0;
        if (run != 0) {
            index_run = this->cols[prev_col].select_runs(run) + 1;
        }

        // compute the previous ones before the run to calculate the offset
        unsigned int ones = uvtrick(prev_col, index_run).second;
        if (verbose) {
            std::cout << "case 2: " << index_run
                      << " + " << index - c << "\n";
        }

        // return the index of the run plus the offset
        // (calculated using also "c")
        return index_run + (index - (c + ones));
    }
}

void rlpbwtbv::external_match_vcf(const char *filename, unsigned int min_len,
                                  bool verbose) {
    std::ifstream in(filename);
    if (in.fail()) {
        throw FileNotFoundException{};
    }

    std::string line;
    while (getline(in, line)) {
        if (line.size() < 2u) {
            throw FileNotGoodException{};
        }
        if (line[0] != '#' || line[1] != '#') {
            break;
        }
    }
    int Q = -9;
    unsigned int tmp_width = 0;
    std::stringstream ss(line);
    while (getline(ss, line, '\t')) {
        Q++;
    }
    if (Q < 1) {
        throw FileNotGoodException{};
    }
    while (std::getline(in, line)) {
        tmp_width++;
    }
    if (tmp_width != this->width) {
        throw NotEqualLengthException{};
    }
    in.clear();
    in.seekg(0);

    while (getline(in, line)) {
        if (line[0] != '#' || line[1] != '#') break;
    }
    ss = std::stringstream(line);

    std::vector<std::string> qIDs;

    for (int i = 0; i < 9; i++) {
        getline(ss, line, '\t');
    }
    std::string tmp;
    for (int i = 0; i < Q; i++) {
        getline(ss, tmp, '\t');
        qIDs.push_back(tmp);
        qIDs.push_back(tmp);
    }
    Q <<= 1;
    std::vector<std::vector<char>> queries_tmp(tmp_width,
                                               std::vector<char>(Q));
    // read query panel
    for (unsigned int k = 0; k < this->width; k++) {
        getline(in, line);
        ss = std::stringstream(line);
        for (int i = 0; i < 9; i++) {
            getline(ss, line, '\t');
        }
        int i = 0;
        while (getline(ss, line, '\t')) {
            queries_tmp[k][i++] = line[0];
            queries_tmp[k][i++] = line[2];
        }
    }
    in.close();
    std::vector<std::string> queries;
    std::string tmpq;
    std::cout << queries_tmp.size() << " " << queries_tmp[0].size() << "\n";
    for (unsigned int i = 0; i < queries_tmp[0].size(); i++) {
        for (auto &j: queries_tmp) {
            tmpq.push_back(j[i]);
        }
        queries.push_back(tmpq);
        tmpq.clear();
    }
    unsigned count = 0;
    for (const auto &s: queries) {
        if (verbose) {
            std::cout << count << "(" << s.size() << "): " << s << "\n";
        }

        auto matches = external_match(s, min_len, verbose);
        if (!matches.empty()) {
            std::cout << "matches with " << count << " " << qIDs[count]
                      << "\n";
            for (const auto &m: matches) {
                std::cout << m << "\n";
            }
        }
        count++;
    }
}
