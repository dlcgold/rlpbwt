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
        unsigned int tmp_width = std::count(
                std::istreambuf_iterator<char>(input_matrix),
                std::istreambuf_iterator<char>(), '\n') + 1;
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
    unsigned int count0 = 0;
    unsigned int u = 0;
    unsigned int v = 0;
    unsigned int count0tmp = 0;
    unsigned int count1 = 0;
    bool start = true;

    for (unsigned int i = 0; i < height; i++) {
        if (i == 0 && column[pref[i]] == '1') {
            start = false;
        }
        if (column[i] == '0') {
            count0++;
        }
    }

    std::vector<std::pair<unsigned int, unsigned int>> rows;
    unsigned int p_tmp = 0;
    bool begrun = true;
    bool pushz = false;
    bool pusho = false;
    if (start) {
        pusho = true;
    } else {
        pushz = true;
    }
    for (unsigned int i = 0; i < height; i++) {
        if (begrun) {
            u = count0tmp;
            v = count1;
            begrun = false;
        }
        if (column[pref[i]] == '1') {
            count1++;
        } else {
            count0tmp++;
        }
        if ((i == 0) || (column[pref[i]] != column[pref[i - 1]])) {
            p_tmp = i;
        }
        if ((i == height - 1) || (column[pref[i]] != column[pref[i + 1]])) {
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
    sdsl::int_vector<> p_vec(rows.size());
    sdsl::int_vector<> uv_vec(rows.size());
    for (unsigned int i = 0; i < rows.size(); i++) {
        p_vec[i] = rows[i].first;
        uv_vec[i] = rows[i].second;
    }
    sdsl::util::bit_compress(p_vec);
    sdsl::util::bit_compress(uv_vec);
    sdsl::util::bit_compress(div);
    return {start, count0, p_vec, uv_vec, div};
}

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
    if (symbol == '0') {
        return uv.first + offset;
    } else {
        return this->cols[col_index].count_0 + uv.second + offset;
    }
}

unsigned int
rlpbwt_naive::reverse_lf(unsigned int col_index, unsigned int index,
                         bool verbose) const {
    col_index = col_index - 1;
    unsigned int c = this->cols[col_index].count_0;
    unsigned int u = 0;
    unsigned int v = 0;
    unsigned int pos = 0;
    unsigned int offset = 0;
    bool found = false;
    if (verbose) {
        std::cout << "c: " << c << "\n";
    }
    if (index < c) {
        u = index;
        if (verbose) {
            std::cout << "u: " << u << "\n";
        }
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
        if (!found) {
            pos = this->cols[col_index].p.size() - 1;
        }
        if (verbose) {
            std::cout << "row: " << pos << "\n";
        }
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
        if (!found) {
            pos = this->cols[col_index].p.size() - 1;
        }
        if (verbose) {
            std::cout << "row: " << pos << "\n";
        }
        unsigned int currv = uvtrick(col_index, pos).second;
        offset = v - currv;
        if (verbose) {
            std::cout << "offset: " << offset << "\n";
        }
        return this->cols[col_index].p[pos] + offset;
    }
    return 0;
}

unsigned int
rlpbwt_naive::index_to_run(unsigned int index, unsigned int col_index) const {
    unsigned int pos = 0;
    bool found_first = false;
    if (index >= this->cols[col_index].p[this->cols[col_index].p.size() - 1]) {
        return this->cols[col_index].p.size() - 1;
    }
    for (unsigned int i = 0; i < this->cols[col_index].p.size() - 1; i++) {
        if (this->cols[col_index].p[i] <= index &&
            index < this->cols[col_index].p[i + 1]) {
            pos = i;
            found_first = true;
            break;
        }
    }
    if (!found_first) {
        pos = this->cols[col_index].p.size() - 1;
    }
    return pos;
}

std::pair<unsigned int, unsigned int>
rlpbwt_naive::uvtrick(unsigned int col_index, unsigned int row_index) const {
    unsigned int u;
    unsigned int v;
    if (row_index == 0) {
        u = 0;
        v = 0;
    } else if (row_index % 2 == 0) {
        u = this->cols[col_index].uv[row_index - 1];
        v = this->cols[col_index].uv[row_index];
        if (!this->cols[col_index].zero_first) {
            std::swap(u, v);
        }
    } else {
        u = this->cols[col_index].uv[row_index];
        v = this->cols[col_index].uv[row_index - 1];
        if (!this->cols[col_index].zero_first) {
            std::swap(u, v);
        }
    }
    return {u, v};
}

matches_naive
rlpbwt_naive::external_match(const std::string &query, bool verbose) {
    if (query.size() != this->width) {
        throw NotEqualLengthException{};
    }
    matches_naive matches;
    unsigned int curr_run = 0;
    unsigned int end_run = 0;
    unsigned int curr_index = 0;
    unsigned int end_index = 0;
    unsigned int curr_tmp = 0;
    unsigned int end_tmp = 0;
    unsigned int curr_offset = 0;
    unsigned int end_offset = 0;
    unsigned int curr_len = 0;
    unsigned int curr_beg = 0;
    for (unsigned i = 0; i < query.size(); i++) {
        if (verbose) {
            std::cout << "before at " << i << " from " << curr_index << " to "
                      << end_index
                      << "\n";
        }
        curr_run = index_to_run(curr_index, i);
        end_run = index_to_run(end_index, i);
        curr_len = (end_index - curr_index);
        curr_offset = curr_index - this->cols[i].p[curr_run];
        end_offset = end_index - this->cols[i].p[end_run];

        // check to undoing the offsets
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
        curr_tmp = lf(i, curr_run, query[i], curr_offset);
        end_tmp = lf(i, end_run, query[i], end_offset);

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
        if (curr_tmp < end_tmp) {
            if (verbose) {
                std::cout << "case 1\n";
            }
            curr_index = curr_tmp;
            end_index = end_tmp;
        } else {
            if (verbose) {
                std::cout << "f: " << curr_tmp << ", g: " << end_tmp << "\n";
            }
            if (i > 0) {
                if (verbose) {
                    std::cout << "match at (" << curr_beg << ", " << i - 1
                              << ") with " << curr_len << " haplotypes \n";
                }
                matches.basic_matches.emplace_back(curr_len, i - curr_beg,
                                                   i - 1);
            }

            // update e
            if (curr_tmp == this->cols[i + 1].lcp.size()) {
                curr_beg = i + 1;
            } else {
                curr_beg = i - this->cols[i + 1].lcp[curr_tmp];
            }
            if (verbose) {
                std::cout << "before curr beg: " << curr_beg << "\n";
            }

            if ((query[curr_beg] == '0' && curr_tmp > 0) ||
                curr_tmp == this->height) {
                if (verbose) {
                    std::cout << "begin case 2 curr run " << curr_run
                              << ", curr index " << curr_index << "\n";
                    std::cout << "currtmp: " << curr_tmp << ", endtmp: "
                              << end_tmp
                              << "\n";
                }
                curr_tmp = end_tmp - 1;
                if (curr_beg > 1) {
                    unsigned int i_tmp = i + 1;
                    unsigned int curr_rev = curr_tmp;
                    // reach curr_neg -1 column
                    while (i_tmp != curr_beg - 1) {
                        curr_rev = reverse_lf(i_tmp, curr_rev, false);
                        i_tmp--;
                    }
                    // first step to reverse procede with lf
                    unsigned int tmp_run = index_to_run(curr_rev, i_tmp);
                    char curr_elem = get_next_char(
                            this->cols[i_tmp].zero_first,
                            tmp_run);
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
                while (curr_tmp > 0 &&
                       (i + 1) - this->cols[i + 1].lcp[curr_tmp] <= curr_beg) {
                    curr_tmp--;
                }
                curr_index = curr_tmp;
                end_index = end_tmp;
                /*curr_run = index_to_run(curr_index, i);
                end_run = index_to_run(end_index, i);*/
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
                end_tmp = curr_tmp + 1;
                if (curr_beg > 1) {
                    unsigned int i_tmp = i + 1;
                    unsigned int curr_rev = curr_tmp;
                    // reach curr_neg -1 column
                    while (i_tmp != curr_beg - 1) {
                        curr_rev = reverse_lf(i_tmp, curr_rev, false);
                        i_tmp--;
                    }

                    // first step to reverse procede with lf
                    unsigned int tmp_run = index_to_run(curr_rev, i_tmp);
                    char curr_elem = get_next_char(
                            this->cols[i_tmp].zero_first,
                            tmp_run);
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
                while (end_tmp < this->height &&
                       (i + 1) - this->cols[i + 1].lcp[end_tmp] <= curr_beg) {
                    end_tmp++;
                }
                curr_index = curr_tmp;
                end_index = end_tmp;
                /*curr_run = index_to_run(curr_index, i);
                end_run = index_to_run(end_index, i);*/
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
    if (curr_index < end_index) {
        curr_len = end_index - curr_index;
        //curr_beg = (query.size() - 1) - this->cols[query.size()].lcp[curr_index];
        if (verbose) {
            std::cout << "match at (" << curr_beg << ", " << query.size() - 1
                      << ") with " << curr_len << " haplotypes \n";
        }
        matches.basic_matches.emplace_back(curr_len, query.size() - curr_beg,
                                           query.size() - 1);

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
                matches_naive matches;

                matches = this->external_match(queries_panel[i], verbose);
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
            }
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



