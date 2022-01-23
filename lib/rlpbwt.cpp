//
// Created by dlcgold on 04/01/22.
//

#include "../include/rlpbwt.h"

rlpbwt::rlpbwt() : cols(), heigth(0), width(0) {}

rlpbwt::rlpbwt(const char *filename, bool vcf, bool verbose) {
    if (vcf) {
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
        std::vector<unsigned int> pref(tmp_height);
        sdsl::int_vector<> div(tmp_height);
        for (int i = 0; i < tmp_height; i++) {
            pref[i] = i;
            div[i] = 0;
        }
        std::string new_column;
        std::vector<column> tmp_cols(tmp_width);
        for (int k = 0; k < tmp_width; k++) {
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
            new_column.clear();
            std::getline(input, line);
            ss = std::stringstream(line);
            for (int i = 0; i < 9; i++) {
                getline(ss, line, '\t');
            }
            int index = 0;
            while (getline(ss, line, '\t')) {
                new_column.push_back(line[0]);
                new_column.push_back(line[2]);
            }
            //std::cout << new_column << "\n";
            auto col = rlpbwt::build_column(new_column, pref, div);
            sdsl::util::bit_compress(div);
            col.lcp = div;
            tmp_cols[k] = col;
            rlpbwt::update(new_column, pref, div);
        }
        auto col = rlpbwt::build_column(new_column, pref, div);
        sdsl::util::bit_compress(div);
        col.lcp = div;
        tmp_cols.push_back(col);
        this->cols = tmp_cols;
        this->width = tmp_width;
        this->heigth = tmp_height;
    } else {
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
            std::vector<column> tmp_cols(tmp_width);
            input_matrix.clear();
            input_matrix.seekg(0, std::ios::beg);
            std::vector<unsigned int> pref(tmp_height);
            sdsl::int_vector<> div(tmp_height);
            for (unsigned int i = 0; i < tmp_height; i++) {
                pref[i] = i;
                div[i] = 0;
            }
            unsigned int count = 0;
            while (getline(input_matrix, new_column)) {
                if (verbose) {
                    std::cout << "\nnew_column " << count << "\n";
                    for (auto e: pref) {
                        std::cout << e << " ";
                    }
                    std::cout << "\n";
                    for (auto e: div) {
                        std::cout << e << " ";
                    }
                    std::cout << "\n";
                }
                new_column.erase(
                        std::remove(new_column.begin(), new_column.end(), ' '),
                        new_column.end());
                auto col = rlpbwt::build_column(new_column, pref, div);
                //sdsl::util::bit_compress(div);
                col.lcp = div;
                tmp_cols[count] = col;
                rlpbwt::update(new_column, pref, div);
                count++;
            }
            auto col = rlpbwt::build_column(new_column, pref, div);
            //sdsl::util::bit_compress(div);
            col.lcp = div;
            tmp_cols.push_back(col);
            this->cols = tmp_cols;
            this->width = tmp_width;
            this->heigth = tmp_height;
            input_matrix.close();
        } else {
            throw FileNotFoundException{};
        }
    }
}

rlpbwt::~rlpbwt() = default;

column
rlpbwt::build_column(std::string &column, std::vector<unsigned int> &pref,
                     sdsl::int_vector<> &div) {
    unsigned int height = pref.size();
    unsigned int count0 = 0;
    unsigned int u = 0;
    unsigned int v = 0;
    unsigned int count0tmp = 0;
    unsigned int count1 = 0;
    unsigned int lcs = 0;
    bool start = true;

    for (unsigned int i = 0; i < height; i++) {
        if (i == 0 && column[pref[i]] == '1') {
            start = false;
        }
        if (column[i] == '0') {
            count0++;
        }
    }

    std::vector<rlrow> rows;
    unsigned int p_tmp = 0;
    bool begrun = true;
    bool pushz = false;
    bool pusho = false;
    if (start) {
        pusho = true;
    } else {
        pushz = true;
    }
    std::vector<unsigned int> uv;
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
            lcs = div[i];
        }

        if (div[i] < lcs) {
            lcs = div[i];
        }
        if ((i == height - 1) ||
            (column[pref[i]] != column[pref[i + 1]])) {
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
    return {start, rows, count0};
}


void rlpbwt::update(std::string &column, std::vector<unsigned int> &pref,
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
rlpbwt::lf(unsigned int col_index, unsigned int row_index, char symbol,
           unsigned int offset, bool verbose) const {
    auto uv = uvtrick(col_index, row_index);
    if (verbose) {
        std::cout << uv.first << ", " << uv.second << "\n";
    }
    if (symbol == '0') {
        return uv.first + offset;
    } else {
        return this->cols[col_index].count_0 + uv.second + offset;
    }
}

unsigned int
rlpbwt::index_to_run(unsigned int index, unsigned int col_index) const {
    unsigned int pos = 0;
    bool found_first = false;
    if (index >= this->cols[col_index].rows.back().p) {
        return this->cols[col_index].rows.size() - 1;
    }
    for (unsigned int i = 0; i < this->cols[col_index].rows.size() - 1; i++) {
        if (this->cols[col_index].rows[i].p <= index &&
            index < this->cols[col_index].rows[i + 1].p) {
            pos = i;
            found_first = true;
            break;
        }
    }
    if (!found_first) {
        pos = this->cols[col_index].rows.size() - 1;
    }
    return pos;
}

std::vector<match>
rlpbwt::external_match(const std::string &query, unsigned int min_len,
                       bool verbose) {
    if (query.size() != this->width) {
        throw NotEqualLengthException{};
    }
    std::vector<match> matches;
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
        curr_offset =
                curr_index - this->cols[i].rows[curr_run].p;
        end_offset =
                end_index - this->cols[i].rows[end_run].p;

        if (query[i] == '0') {
            if (curr_run == 0 &&
                get_next_char(this->cols[i].zero_first,
                              curr_run) == '1') {
                curr_offset = 0;
            } else if (end_run == 0 &&
                       get_next_char(this->cols[i].zero_first,
                                     end_run) == '1') {
                end_offset = 0;
            }
        }
        if (query[i] == '1') {
            if (curr_run == 1 &&
                get_next_char(this->cols[i].zero_first,
                              curr_run) == '0') {
                curr_offset = 0;
            } else if (end_run == 1 &&
                       get_next_char(this->cols[i].zero_first,
                                     end_run) == '0') {
                end_offset = 0;
            }
        }
        curr_tmp = lf(i, curr_run, query[i], curr_offset);
        end_tmp = lf(i, end_run, query[i], end_offset);
        if (curr_tmp > this->heigth) {
            curr_tmp -= curr_offset;
        }
        if (end_tmp > this->heigth) {
            end_tmp -= end_offset;
        }
        /*curr_run = index_to_run(curr_tmp, i);
        end_run = index_to_run(end_tmp, i);
        curr_index = this->cols[i].rows[curr_run].p;
        end_index = this->cols[i].rows[end_run].p;*/
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
                if ((i - 1) - curr_beg >= min_len) {
                    matches.emplace_back(curr_beg, i - 1, curr_len);
                }
            }

            // update e
            // TODO check correctness
            if (curr_tmp == this->cols[i + 1].lcp.size()) {
                curr_beg = i + 1;
            } else if ((int) i - (int) this->cols[i + 1].lcp[curr_tmp] < 0) {
                curr_beg = 0;
            }else {
                curr_beg = i - this->cols[i + 1].lcp[curr_tmp];
            }
            if (verbose) {
                std::cout << "before curr beg: " << curr_beg << "\n";
            }

            if ((query[curr_beg] == '0' && curr_tmp > 0) ||
                curr_tmp == this->heigth) {
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
                        curr_rev = prev_run(i_tmp, curr_rev, false);
                        i_tmp--;
                    }
                    // first step to reverse procede with lf
                    unsigned int tmp_run = index_to_run(curr_rev, i_tmp);
                    char curr_elem = get_next_char(
                            this->cols[i_tmp].zero_first,
                            tmp_run);
                    while (curr_beg > 0 && query[curr_beg - 1] == curr_elem) {
                        curr_beg -= 1;
                        if(curr_beg >0) {
                            curr_rev = prev_run(curr_beg, curr_rev, false);
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
                        curr_rev = prev_run(i_tmp, curr_rev, false);
                        i_tmp--;
                    }

                    // first step to reverse procede with lf
                    unsigned int tmp_run = index_to_run(curr_rev, i_tmp);
                    char curr_elem = get_next_char(
                            this->cols[i_tmp].zero_first,
                            tmp_run);
                    while (curr_beg > 0 && query[curr_beg - 1] == curr_elem) {
                        curr_beg -= 1;
                        if(curr_beg > 0) {
                            curr_rev = prev_run(curr_beg, curr_rev, false);
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
                while (end_tmp < this->heigth &&
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
        curr_beg = (query.size() - 1) - this->cols[query.size()].lcp[curr_tmp];
        if (verbose) {
            std::cout << "match at (" << curr_beg << ", " << query.size() - 1
                      << ") with " << curr_len << " haplotypes \n";
        }
        if ((query.size() - 1) - curr_beg >= min_len) {
            matches.emplace_back(curr_beg, query.size() - 1, curr_len);
    }
    }
    return matches;
}

unsigned int rlpbwt::prev_run(unsigned int col_index, unsigned int index,
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
             i < this->cols[col_index].rows.size() - 1; i++) {
            prevu = uvtrick(col_index, i).first;
            nextu = uvtrick(col_index, i + 1).first;
            if (prevu <= u && u < nextu) {
                pos = i;
                found = true;
                break;
            }
        }
        if (!found) {
            pos = this->cols[col_index].rows.size() - 1;
        }
        if (verbose) {
            std::cout << "pos: " << pos << "\n";
        }
        unsigned int curru = uvtrick(col_index, pos).first;

        offset = u - curru;
        if (verbose) {
            std::cout << "offset: " << offset << "\n";
        }
        return this->cols[col_index].rows[pos].p + offset;
    } else {
        unsigned int prevv = 0;
        unsigned int nextv = 0;
        v = index - c;
        if (verbose) {
            std::cout << "v: " << v << "\n";
        }
        for (unsigned int i = 0;
             i < this->cols[col_index].rows.size() - 1; ++i) {
            prevv = uvtrick(col_index, i).second;
            nextv = uvtrick(col_index, i + 1).second;
            if (prevv <= v && v < nextv) {
                pos = i;
                found = true;
                break;
            }
        }
        if (!found) {
            pos = this->cols[col_index].rows.size() - 1;
        }
        if (verbose) {
            std::cout << "pos: " << pos << "\n";
        }
        unsigned int currv = uvtrick(col_index, pos).second;
        offset = v - currv;
        if (verbose) {
            std::cout << "offset: " << offset << "\n";
        }
        return this->cols[col_index].rows[pos].p + offset;
    }
    return 0;
}

std::pair<unsigned int, unsigned int>
rlpbwt::uvtrick(unsigned int col_index, unsigned int row_index) const {
    unsigned int u;
    unsigned int v;
    if (row_index == 0) {
        u = 0;
        v = 0;
    } else if (row_index % 2 == 0) {
        u = this->cols[col_index].rows[row_index - 1].uv;
        v = this->cols[col_index].rows[row_index].uv;
        if (!this->cols[col_index].zero_first) {
            std::swap(u, v);
        }
    } else {
        u = this->cols[col_index].rows[row_index].uv;
        v = this->cols[col_index].rows[row_index - 1].uv;
        if (!this->cols[col_index].zero_first) {
            std::swap(u, v);
        }
    }
    return {u, v};
}

std::vector<match>
rlpbwt::end_external_match(const std::string &query, bool forward,
                           bool verbose) {
    if (query.size() != this->width) {
        throw NotEqualLengthException{};
    }
    std::vector<match> matches;
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
    unsigned int i;
    for (i = 0; i < query.size(); i++) {
        if (verbose) {
            std::cout << "before at " << i << " from " << curr_index << " to "
                      << end_index
                      << "\n";
        }
        curr_run = index_to_run(curr_index, i);
        end_run = index_to_run(end_index, i);
        curr_len = (end_index - curr_index);
        curr_offset = curr_index - this->cols[i].rows[curr_run].p;
        end_offset = end_index - this->cols[i].rows[end_run].p;

        if (query[i] == '0') {
            if (curr_run == 0 &&
                get_next_char(this->cols[i].zero_first,
                              curr_run) == '1') {
                curr_offset = 0;
            } else if (end_run == 0 &&
                       get_next_char(this->cols[i].zero_first,
                                     end_run) == '1') {
                end_offset = 0;
            }
        }
        if (query[i] == '1') {
            if (curr_run == 1 &&
                get_next_char(this->cols[i].zero_first,
                              curr_run) == '0') {
                curr_offset = 0;
            } else if (end_run == 1 &&
                       get_next_char(this->cols[i].zero_first,
                                     end_run) == '0') {
                end_offset = 0;
            }
        }

        curr_tmp = lf(i, curr_run, query[i], curr_offset);
        end_tmp = lf(i, end_run, query[i], end_offset);

        if (curr_tmp > this->heigth) {
            curr_tmp -= curr_offset;
        }
        if (end_tmp > this->heigth) {
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
                if (forward) {
                    matches.emplace_back(curr_beg, i - 1, curr_len);
                } else {
                    unsigned int rev_beg = query.size() - (i - 1) - 1;
                    unsigned int rev_end = query.size() - curr_beg - 1;
                    matches.emplace_back(rev_beg, rev_end, curr_len);
                }
            }
            curr_beg = i;
            if (verbose) {
                std::cout << "before curr beg: " << curr_beg << "\n";
            }

            if (query[curr_beg] == '0') {
                curr_index = 0;
                end_index = this->cols[i].count_0;
            } else {
                curr_index = this->cols[i].count_0;
                end_index = this->heigth;
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
        //curr_beg = (query.size() - 1) - this->cols[query.size()].lcp[curr_tmp];
        //curr_beg = i;
        if (verbose) {
            std::cout << "match at (" << curr_beg << ", " << query.size() - 1
                      << ") with " << curr_len << " haplotypes \n";
        }
        if (forward) {
            if(curr_len>0) {
                matches.emplace_back(curr_beg, query.size() - 1, curr_len);
            }
        } else {
            if(curr_len>0) {
                unsigned int rev_beg = query.size() - (query.size() - 1) - 1;
                unsigned int rev_end = query.size() - curr_beg;
                matches.emplace_back(rev_beg, rev_end, curr_len);
            }
        }
    }
    if (!forward) {
        std::reverse(matches.begin(), matches.end());
    }
    return matches;
}

void rlpbwt::print() {
    int count = 0;
    for (const auto &c: this->cols) {
        std::string z;
        if (c.zero_first) {
            z = "yes";
        } else {
            z = "no";
        }
        std::cout << "column: " << count << "\nstart with 0? "
                  << z << ", c: " << c.count_0 << "\n";

        for (const auto &r: c.rows) {
            std::cout << r << "\n";
        }
        for (auto d: c.lcp) {
            std::cout << d << " ";
        }
        count++;
        if (count != (int) this->cols.size())
            std::cout << "\n-------------- \n";
    }
}

void rlpbwt::external_match_vcf(const char *filename, unsigned int min_len,
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
    std::cout  << queries_tmp.size() << " " << queries_tmp[0].size() << "\n";
    for (unsigned int i = 0; i < queries_tmp[0].size(); i++) {
        for (unsigned int j = 0; j < queries_tmp.size(); j++) {
            tmpq.push_back(queries_tmp[j][i]);
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
            std::cout << "matches with " << qIDs[count] << "\n";
            for (const auto &m: matches) {
                std::cout << m << "\n";
            }
        }
        count++;
    }
}








