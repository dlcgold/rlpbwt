//
// Created by dlcgold on 21/02/22.
//

#include "../include/rlpbwt_naive.h"
#include "../include/exceptions.h"


rlpbwt_naive::rlpbwt_naive() = default;

rlpbwt_naive::~rlpbwt_naive() = default;

rlpbwt_naive::rlpbwt_naive(const char *filename, bool verbose) {
    std::ifstream input_matrix(filename);
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
        // unsigned int tmp_width = std::count(
        //         std::istreambuf_iterator<char>(input_matrix),
        //         std::istreambuf_iterator<char>(), '\n');
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
        input_matrix.clear();
        input_matrix.seekg(0, std::ios::beg);
        this->cols = std::vector<column_naive>(tmp_width + 1);
        std::vector<unsigned int> pref(tmp_height);
        sdsl::int_vector<> div(tmp_height);
        for (unsigned int i = 0; i < tmp_height; i++) {
            pref[i] = i;
            div[i] = 0;
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
            auto col = rlpbwt_naive::build_column(new_column, pref, div);

            this->cols[count] = col;

            rlpbwt_naive::update(new_column, pref, div);
            last_col = new_column;
            count++;
        }
        auto col = rlpbwt_naive::build_column(last_col, pref, div);
        this->cols[count] = col;
        this->height = tmp_height;
        this->width = tmp_width;
        input_matrix.close();
    } else {
        throw FileNotFoundException{};
    }
}

column_naive
rlpbwt_naive::build_column(std::string &column, std::vector<unsigned int> &pref,
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

    // initialize a vector of pair in order to build final sdsl int_vector for
    // p and u/v
    std::vector<std::pair<unsigned int, unsigned int>> rows;
    // temporary variable for p
    unsigned int p_tmp = 0;
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
            begrun = true;
        }
    }
    // create and compress the two sdsl int_vector for p and u/v
    sdsl::int_vector<> p_vec(rows.size());
    sdsl::int_vector<> uv_vec(rows.size());
    for (unsigned int i = 0; i < rows.size(); i++) {
        p_vec[i] = rows[i].first;
        uv_vec[i] = rows[i].second;
    }
    sdsl::util::bit_compress(p_vec);
    sdsl::util::bit_compress(uv_vec);

    // compress also the divergence/lcp vector
    sdsl::util::bit_compress(div);
    return {start, count0, p_vec, uv_vec, div};
}

// algorithm 2 from Durbin's paper to update prefix and divergence arrays
void rlpbwt_naive::update(std::string &column, std::vector<unsigned int> &pref,
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

std::ostream &operator<<(std::ostream &os, const rlpbwt_naive &naive) {
    os << " height: " << naive.height << " width: " << naive.width << ":\n";
    for (unsigned int i = 0; i < naive.cols.size(); i++) {
        std::cout << "column " << i << ":\n" << naive.cols[i]
                  << "-----------------\n";
    }
    return os;
}

unsigned int
rlpbwt_naive::lf(unsigned int col_index, unsigned int row_index, char symbol,
                 unsigned int offset, bool verbose) const {
    // obtain "u" and "v"
    auto uv = uvtrick(col_index, row_index);
    if (verbose) {
        std::cout << uv.first << ", " << uv.second << "\n";
    }
    // fix for the last index that's "outside" the column
    if (this->cols[col_index].p[row_index] + offset == this->height) {
        if (get_next_char(this->cols[col_index].zero_first, row_index) == '0') {
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

unsigned int
rlpbwt_naive::reverse_lf(unsigned int col_index, unsigned int index,
                         bool verbose) const {
    // by design if we try to work on first column the function return 0
    if (col_index == 0) {
        return 0;
    }
    // we extract the "c" value from the previous column
    col_index = col_index - 1;
    unsigned int c = this->cols[col_index].count_0;
    // initialize u/v, offset and the new run index
    unsigned int u = 0;
    unsigned int v = 0;
    unsigned int offset = 0;
    unsigned int pos = 0;

    // bool to interrupt the search
    bool found = false;
    if (verbose) {
        std::cout << "c: " << c << "\n";
    }
    // two cases:
    // - if index is less than previous "c" it means that it comes from a zero
    //   element, so we will search the correct u value
    // - otherwise it means that it comes from a one element, so we will search
    //   the correct v value
    if (index < c) {
        u = index;
        if (verbose) {
            std::cout << "u: " << u << "\n";
        }
        // iteration over the u values to find the correct one
        unsigned int prevu = 0;
        unsigned int nextu = 0;
        for (unsigned int i = 0;
             i < this->cols[col_index].p.size() - 1; i++) {
            prevu = uvtrick(col_index, i).first;
            nextu = uvtrick(col_index, i + 1).first;
            if (prevu <= u && u < nextu) {
                pos = i;
                found = true;
                break;
            }
        }
        // if not found we are at the last run
        if (!found) {
            pos = this->cols[col_index].p.size() - 1;
        }
        if (verbose) {
            std::cout << "row: " << pos << "\n";
        }

        // using the correct u value in previous column to obtain the previous
        // index using offset
        unsigned int curru = uvtrick(col_index, pos).first;
        offset = u - curru;
        if (verbose) {
            std::cout << "offset: " << offset << "\n";
        }
        return this->cols[col_index].p[pos] + offset;
    } else {
        unsigned int prevv = 0;
        unsigned int nextv = 0;
        v = index - c;
        if (verbose) {
            std::cout << "v: " << v << "\n";
        }
        // iteration over the v values to find the correct one
        for (unsigned int i = 0;
             i < this->cols[col_index].p.size() - 1; ++i) {
            prevv = uvtrick(col_index, i).second;
            nextv = uvtrick(col_index, i + 1).second;
            if (prevv <= v && v < nextv) {
                pos = i;
                found = true;
                break;
            }
        }
        // if not found we are at the last run
        if (!found) {
            pos = this->cols[col_index].p.size() - 1;
        }
        if (verbose) {
            std::cout << "row: " << pos << "\n";
        }
        // using the correct v value in previous column to obtain the previous
        // index using offset
        unsigned int currv = uvtrick(col_index, pos).second;
        offset = v - currv;
        if (verbose) {
            std::cout << "offset: " << offset << "\n";
        }
        return this->cols[col_index].p[pos] + offset;
    }
    return 0;
}

// TODO optimize this function
unsigned int
rlpbwt_naive::index_to_run(unsigned int index, unsigned int col_index) const {
    unsigned int pos = 0;
    bool found_first = false;
    // if requested index is equal or greater than the p value of the last run
    // return the index of the last run
    if (index >= this->cols[col_index].p[this->cols[col_index].p.size() - 1]) {
        return this->cols[col_index].p.size() - 1;
    }

    // iterate over the runs in order to find the correct run that contain the
    // index using a run and the next one
    for (unsigned int i = 0; i < this->cols[col_index].p.size() - 1; i++) {
        if (this->cols[col_index].p[i] <= index &&
            index < this->cols[col_index].p[i + 1]) {
            pos = i;
            found_first = true;
            break;
        }
    }

    // if not found in previous iteration the index is in the last run
    if (!found_first) {
        pos = this->cols[col_index].p.size() - 1;
    }
    return pos;
}

std::pair<unsigned int, unsigned int>
rlpbwt_naive::uvtrick(unsigned int col_index, unsigned int run_index) const {
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

matches_naive
rlpbwt_naive::external_match(const std::string &query, bool verbose) {
    // query allowed iff |query| is equal to RLPBWT width
    if (query.size() != this->width) {
        throw NotEqualLengthException{};
    }
    // initialize basic_matches vector
    matches_naive matches;

    // initialize all variables to support the computation
    // "curr" variables refers to "f" and "end" to "g" in Durbin's paper
    unsigned int curr_run = 0;
    unsigned int end_run = 0;
    unsigned int curr_index = 0;
    unsigned int end_index = 0;
    unsigned int curr_tmp = 0;
    unsigned int end_tmp = 0;

    // initialize offsets to obrain real indices
    unsigned int curr_offset = 0;
    unsigned int end_offset = 0;

    // current amount of rows that are matching
    unsigned int curr_len = 0;

    // current begin of a match, as "e" in Durbin's paper
    unsigned int curr_beg = 0;

    // iterate over every column
    for (unsigned i = 0; i < query.size(); i++) {
        //std::cout << "processed " << i << "\r";
        if (verbose) {
            std::cout << "before at " << i << " from " << curr_index << " to "
                      << end_index
                      << "\n";
        }

        // compute run of the current indices using rank over the runs bv
        curr_run = index_to_run(curr_index, i);
        end_run = index_to_run(end_index, i);

        // update current number of rows that are matching
        curr_len = (end_index - curr_index);

        // compute the two offsets
        curr_offset = curr_index - this->cols[i].p[curr_run];
        end_offset = end_index - this->cols[i].p[end_run];

        // undoing the offsets when they are wrong/useless
        if (query[i] == '0') {
            if (get_next_char(this->cols[i].zero_first,
                              curr_run) == '1') {
                curr_offset = 0;
            }
            if (get_next_char(this->cols[i].zero_first,
                              end_run) == '1') {
                end_offset = 0;
            }
        }
        if (query[i] == '1') {
            if (get_next_char(this->cols[i].zero_first,
                              curr_run) == '0') {
                curr_offset = 0;
            }
            if (get_next_char(this->cols[i].zero_first,
                              end_run) == '0') {
                end_offset = 0;
            }
        }

        // lf-mapping of the two indices
        curr_tmp = lf(i, curr_run, query[i], curr_offset);
        end_tmp = lf(i, end_run, query[i], end_offset);

        // fix at column end
        if (curr_tmp > this->height) {
            curr_tmp -= curr_offset;
        }
        if (end_tmp > this->height) {
            end_tmp -= end_offset;
        }
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
            if (verbose) {
                std::cout << "f: " << curr_tmp << ", g: " << end_tmp << "\n";
            }
            // save basic_matches if longer than 1
            if (i > 0) {
                if (i - curr_beg > 1) {
                    if (verbose) {
                        std::cout << "match at (" << curr_beg << ", " << i - 1
                                  << ") with " << curr_len << " haplotypes \n";
                    }
                    matches.basic_matches.emplace_back(curr_len, i - curr_beg,
                                                       i - 1);
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
            //   while end index will be increased
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
                // we look back to update index of begin of basic_matches
                // (iff it's not already at the minimum possible, 0)
                // we don't have both the panel and the prefix array in
                // memory, so we have to proceed reversing the permutations
                if (curr_beg > 1) {
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
                    unsigned int tmp_run = index_to_run(curr_rev, i_tmp);
                    char curr_elem = get_next_char(
                            this->cols[i_tmp].zero_first,
                            tmp_run);

                    // we virtually follow the row of the temporary index
                    // in order to update the beginning index
                    while (curr_beg > 0 && query[curr_beg - 1] == curr_elem) {
                        curr_beg -= 1;
                        if (curr_beg > 0) {
                            curr_rev = reverse_lf(curr_beg, curr_rev, false);
                            tmp_run = index_to_run(curr_rev, curr_beg - 1);
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

                // we look back to update index of begin of basic_matches
                // (iff it's not already at the minimum possible, 0)
                // we don't have both the panel and the prefix array in
                // memory, so we have to proceed reversing the permutations
                if (curr_beg > 1) {
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
                    unsigned int tmp_run = index_to_run(curr_rev, i_tmp);
                    char curr_elem = get_next_char(
                            this->cols[i_tmp].zero_first,
                            tmp_run);

                    // we virtually follow the row of the temporary index
                    // in order to update the beginning index
                    while (curr_beg > 0 && query[curr_beg - 1] == curr_elem) {
                        curr_beg -= 1;
                        if (curr_beg > 0) {
                            curr_rev = reverse_lf(curr_beg, curr_rev, false);
                            tmp_run = index_to_run(curr_rev, curr_beg - 1);
                            curr_elem = get_next_char(
                                    this->cols[curr_beg - 1].zero_first,
                                    tmp_run);
                        }
                    }
                }
                if (verbose) {
                    std::cout << "end curr beg: " << curr_beg << "\n";
                }
                // after the computation of the new beginning of a match we use
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
    // check/save matches at the end (iff length > 1)
    if (curr_index < end_index) {
        curr_len = end_index - curr_index;
        if (verbose) {
            std::cout << "match at (" << curr_beg << ", " << query.size() - 1
                      << ") with " << curr_len << " haplotypes \n";
        }
        if (query.size() - curr_beg > 1) {
            matches.basic_matches.emplace_back(curr_len,
                                               query.size() - curr_beg,
                                               query.size() - 1);
        }
    }
    return matches;
}

void
rlpbwt_naive::match_tsv(const char *filename, const char *out, bool verbose) {
    std::ifstream input_matrix(filename);
    std::ofstream out_match(out);
    if (input_matrix.is_open()) {
        std::string header1;
        std::string header2;
        std::string line;
        std::string garbage;
        std::string new_column;
        getline(input_matrix, header1);
        getline(input_matrix, header2);
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
                matches_naive matches;
                if (verbose) {
                    std::cout << query << "\n";
                }
                matches = this->external_match(query, verbose);
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

void rlpbwt_naive::match_tsv_tr(const char *filename, const char *out,
                                bool verbose) {
    std::ifstream input_matrix(filename);
    std::ofstream out_match(out);
    if (input_matrix.is_open()) {
        std::string header1;
        std::string header2;
        std::string line;
        std::string garbage;
        std::string new_column;
        getline(input_matrix, header1);
        getline(input_matrix, header2);
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
                matches_naive matches;
                matches = this->external_match(query, verbose);
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
rlpbwt_naive::match_tsv_conc(const char *filename, const char *out, bool verbose) {
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
            std::vector<matches_naive> matches_vec(n_queries);
#pragma omp parallel for default(none) shared(queries, matches_vec, n_queries, verbose)
            for (unsigned int i = 0; i < n_queries; i++) {
                //std::cout << i << "\n";
                matches_vec[i] = this->external_match(queries[i], verbose);
            }
            for (unsigned int i = 0; i < queries.size(); i++) {
                for (unsigned int j = 0;
                     j < matches_vec[i].basic_matches.size(); j++) {
                    auto len = std::get<1>(
                            matches_vec[i].basic_matches[j]);
                    auto end = std::get<2>(
                            matches_vec[i].basic_matches[j]);
                    for (unsigned int k = 0;
                         k < std::get<0>(
                                 matches_vec[i].basic_matches[j]); k++) {
                        out_match << "MATCH\t" << i << "\t?\t"
                                  << end - (len - 1) << "\t" << end
                                  << "\t" << len << "\n";
                    }

                }
            }
            /*
            for (unsigned int i = 0; i < queries.size(); i++) {
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
            }
             */
            out_match.close();
        } else {
            throw FileNotFoundException{};
        }

    } else {
        throw FileNotFoundException{};
    }
}

unsigned long long rlpbwt_naive::size_in_bytes(bool verbose) {
    unsigned long long size = 0;
    unsigned long long size_p = 0;
    unsigned long long size_uv = 0;
    unsigned long long size_lcp = 0;
    for (const auto &c: this->cols) {
        size += c.size_in_bytes();
        size_p += sdsl::size_in_bytes(c.p);
        size_uv += sdsl::size_in_bytes(c.uv);
        size_lcp += sdsl::size_in_bytes(c.lcp);
    }
    if (verbose) {
        std::cout << "p: " << size_p << " bytes\n";
        std::cout << "uv: " << size_uv << " bytes\n";
        std::cout << "lcp: " << size_lcp << " bytes\n";
        std::cout << "rlpbwt (with also c values and other support variables): "
                  << size << " bytes\n";
    }
    size += sizeof(bool);
    size += sizeof(unsigned int);

    return size;
}

double rlpbwt_naive::size_in_mega_bytes(bool verbose) {
    double size = 0;
    double size_p = 0;
    double size_uv = 0;
    double size_lcp = 0;
    double to_mega = ((double) 1 / (double) 1024) / (double) 1024;
    for (const auto &c: this->cols) {
        size += c.size_in_mega_bytes();
        size_p += sdsl::size_in_mega_bytes(c.p);
        size_uv += sdsl::size_in_mega_bytes(c.uv);
        size_lcp += sdsl::size_in_mega_bytes(c.lcp);
    }
    if (verbose) {
        std::cout << "p: " << size_p << " megabytes\n";
        std::cout << "uv: " << size_uv << " megabytes\n";
        std::cout << "lcp: " << size_lcp << " megabytes\n";
        std::cout << "rlpbwt (with also c values and other support variables): "
                  << size << " megabytes\n";
    }
    size += (double) (sizeof(bool) * to_mega);
    size += (double) (sizeof(unsigned int) * to_mega);

    return size;
}

size_t rlpbwt_naive::serialize(std::ostream &out, sdsl::structure_tree_node *v,
                               const std::string &name) {
    sdsl::structure_tree_node *child =
            sdsl::structure_tree::add_child(v, name,
                                            sdsl::util::class_name(
                                                    *this));
    size_t written_bytes = 0;
    out.write((char *) &this->height, sizeof(this->height));
    written_bytes += sizeof(this->height);

    out.write((char *) &this->width, sizeof(this->width));
    written_bytes += sizeof(this->width);

    for (unsigned int i = 0; i < this->cols.size(); i++) {
        std::string label = "col_" + std::to_string(i);
        written_bytes += this->cols[i].serialize(out, child, label);
    }

    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

void rlpbwt_naive::load(std::istream &in) {
    in.read((char *) &this->height, sizeof(this->height));
    in.read((char *) &this->width, sizeof(this->width));
    //auto c = new column_naive();
    this->cols = std::vector<column_naive>(this->width + 1);
    for (unsigned int i = 0; i <= this->width; i++) {
      /*
        c->load(in);
        this->cols.emplace_back(*c);
	*/
	this->cols[i].load(in);
    }
}

unsigned int rlpbwt_naive::get_run_number() {
    unsigned int count_run = 0;
    for (auto &col: this->cols) {
        count_run += col.p.size();
    }
    return count_run;
}



