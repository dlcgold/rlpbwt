//
// Created by dlcgold on 04/01/22.
//

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include "../include/utils.h"
#include "../include/rlpbwtc.h"
#include "../include/exceptions.h"

rlpbwtc::rlpbwtc(const char *filename, bool verbose) {
    std::ifstream input_matrix(filename);
    if (input_matrix.is_open()) {
        std::string column;
        getline(input_matrix, column);
        column.erase(std::remove(column.begin(), column.end(), ' '),
                     column.end());
        const unsigned int tmp_height = column.size();
        unsigned int tmp_width = std::count(
                std::istreambuf_iterator<char>(input_matrix),
                std::istreambuf_iterator<char>(), '\n') + 1;
        std::vector<rlpbwt_column> tmp_cols(tmp_width);
        input_matrix.clear();
        input_matrix.seekg(0, std::ios::beg);
        std::vector<unsigned int> pref(tmp_height);
        std::vector<unsigned int> div(tmp_height);
        for (unsigned int i = 0; i < tmp_height; i++) {
            pref[i] = i;
            div[i] = 0;
        }
        unsigned int count = 0;
        while (getline(input_matrix, column)) {
            if (verbose) {
                std::cout << "\ncolumn " << count << "\n";
                for (auto e: pref) {
                    std::cout << e << " ";
                }
                std::cout << "\n";
                for (auto e: div) {
                    std::cout << e << " ";
                }
                std::cout << "\n";
            }
            column.erase(std::remove(column.begin(), column.end(), ' '),
                         column.end());
            auto col = rlpbwtc::build_column(column, pref, div);
            tmp_cols[count] = col;

            //rlpbwtc::update_old(column, pref, div, count);
            std::cout << "build at column " << count << "\n";
            for (auto e: div) {
                std::cout << e << " ";
            }
            std::cout << std::endl;
            rlpbwtc::update(column, pref, div);
            count++;
        }
        this->cols = tmp_cols;
        this->width = tmp_width;
        this->heigth = tmp_height;
        input_matrix.close();
    } else {
        throw FileNotFoundException{};
    }
}

rlpbwtc::~rlpbwtc() = default;

std::string rlpbwtc::search_row(unsigned int row_index, bool verbose) {
    unsigned int start;
    unsigned int pos = 0;
    std::string row;
    bool found_first = false;
    for (unsigned int i = 0; i < this->cols[0].rows.size() - 1; i++) {
        if (this->cols[0].rows[i].p <= row_index &&
            row_index < this->cols[0].rows[i + 1].p) {
            pos = i;
            row.push_back(get_next_char(this->cols[0].zero_first, pos));
            found_first = true;
            break;
        }
    }
    if (!found_first) {
        pos = this->cols[0].rows.size() - 1;
        row.push_back(get_next_char(this->cols[0].zero_first, pos));
    }

    start = this->cols[0].rows[pos].next_perm;
    unsigned int end = this->cols[0].rows[pos].lf_mapping(row_index);

    for (unsigned int i = 1; i < this->cols.size(); i++) {
        if (start == this->cols[i].rows.size() - 1) {
            row.push_back(get_next_char(this->cols[i].zero_first, start));
            end = this->cols[i].rows[start].lf_mapping(end);
            start = this->cols[i].rows[start].next_perm;
            if (verbose) {
                std::cout << "column " << i << ": " << start << ", " << end
                          << "\n";
            }
        } else {
            bool found = false;
            for (unsigned int j = start; j < cols[i].rows.size() - 1; j++) {
                if (cols[i].rows[j].p <= end && end < cols[i].rows[j + 1].p) {
                    row.push_back(get_next_char(this->cols[i].zero_first, j));
                    found = true;
                    end = this->cols[i].rows[j].lf_mapping(end);
                    start = this->cols[i].rows[j].next_perm;
                    if (verbose) {
                        std::cout << "column " << i << ": " << start << ", "
                                  << end
                                  << "\n";
                    }
                    break;
                }
            }
            if (!found) {
                unsigned int endrow = this->cols[i].rows.size() - 1;
                row.push_back(get_next_char(this->cols[i].zero_first, endrow));
                end = this->cols[i].rows[endrow].lf_mapping(end);
                start = this->cols[i].rows[endrow].next_perm;
                if (verbose) {
                    std::cout << "column " << i << ": " << start << ", " << end
                              << "\n";
                }
            }
        }
    }
    return row;
}

rlpbwt_column
rlpbwtc::build_column(std::string &column, std::vector<unsigned int> &pref,
                      std::vector<unsigned int> &div) {
    unsigned int height = pref.size();
    unsigned int count0 = 0;
    unsigned int count0tmp = 0;
    unsigned int count1tmp = 0;
    unsigned int count0tmptmp = 0;
    unsigned int count1 = 0;
    unsigned int threshold = 0;
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

    std::vector<rlpbwt_rlrow> rows;
    unsigned int p_tmp = 0;
    unsigned int perm_tmp = 0;
    bool begrun = true;
    for (unsigned int i = 0; i < height; i++) {
        if (begrun) {
            count0tmp = count0tmptmp;
            count1tmp = count1;
            begrun = false;
        }
        if (column[pref[i]] == '1') {
            count1++;
        } else {
            count0tmptmp++;
        }

        if ((i == 0) || (column[pref[i]] != column[pref[i - 1]])) {
            p_tmp = i;
            if (column[pref[i]] == '0') {
                perm_tmp = i - count1;
            } else {
                perm_tmp = count0 + count1 - 1;
            }
            threshold = i;
            lcs = div[i];
        }

        if (div[i] < lcs) {
            threshold = i;
            lcs = div[i];
        }
        if ((i == height - 1) ||
            (column[pref[i]] != column[pref[i + 1]])) {
            if (p_tmp == 0) {
                rows.emplace_back(p_tmp, 0, 0, threshold);
            } else {
                rows.emplace_back(p_tmp, count0tmp, count1tmp, threshold);
            }
            begrun = true;
        }
    }
    return {start, rows, count0};
}


void rlpbwtc::update(std::string &column, std::vector<unsigned int> &pref,
                     std::vector<unsigned int> &div) {
    unsigned int height = pref.size();
    std::vector<unsigned int> new_pref(height);
    std::vector<unsigned int> new_div(height);
    unsigned int count0 = 0;
    unsigned int lcs = -1;

    for (unsigned int i = 0; i < height; i++) {
        lcs = std::min(lcs, div[i]);
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
        lcs = std::min(lcs, div[i]);
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
rlpbwtc::next_run(unsigned int col_index, unsigned int start, unsigned int end,
                  bool verbose) const {
    if (start == cols[col_index].rows.size() - 1) {
        end = this->cols[col_index].rows[start].lf_mapping(end);
        start = this->cols[col_index].rows[start].next_perm;
        if (verbose) {
            std::cout << "column " << col_index << ": " << start << ", " << end
                      << "\n";
        }
    } else {
        bool found = false;
        for (unsigned int j = start; j < cols[col_index].rows.size() - 1; j++) {
            if (cols[col_index].rows[j].p <= end &&
                end < cols[col_index].rows[j + 1].p) {
                found = true;
                end = this->cols[col_index].rows[j].lf_mapping(end);
                start = this->cols[col_index].rows[j].next_perm;
                if (verbose) {
                    std::cout << "column " << col_index << ": " << start << ", "
                              << end
                              << "\n";
                }
                break;
            }
        }
        if (!found) {
            unsigned int endrow = this->cols[col_index].rows.size() - 1;
            end = this->cols[col_index].rows[endrow].lf_mapping(end);
            start = this->cols[col_index].rows[endrow].next_perm;
            if (verbose) {
                std::cout << "column " << col_index << ": " << start << ", "
                          << end
                          << "\n";
            }
        }
    }
    return end;
}

unsigned int
rlpbwtc::occ(unsigned int col_index, unsigned int row_index, char symbol,
             unsigned int offset, bool verbose) const {
    //unsigned int u = 0;
    //unsigned int v = 0;
    //unsigned int c = this->cols[col_index].count_0;
//    for (unsigned int i = 0; i < row_index; i++) {
//        if (verbose && false) {
//            std::cout << this->cols[col_index].rows[i + 1].p << " - "
//                      << this->cols[col_index].rows[i].p << "\n";
//        }
//        if (i % 2 == 0) {
//            if (i != this->cols[col_index].rows.size() - 1) {
//                u += this->cols[col_index].rows[i + 1].p -
//                     this->cols[col_index].rows[i].p;
//            } else {
//                u += this->heigth - this->cols[col_index].rows[i].p;
//            }
//        } else {
//            if (i != this->cols[col_index].rows.size() - 1) {
//                v += this->cols[col_index].rows[i + 1].p -
//                     this->cols[col_index].rows[i].p;
//            } else {
//                v += this->heigth - this->cols[col_index].rows[i].p;
//            }
//        }
//    }
//    if (!this->cols[col_index].zero_first) {
//        std::swap(u, v);
//    }
//    if (verbose) {
//        std::cout << "symbol " << symbol << " at col " << col_index
//                  << " in run "
//                  << row_index << " c: " << c << ", u: " << u << ", v: " << v
//                  << "\n";
//    }
//    if (get_next_char(this->cols[col_index].zero_first, row_index) == '1'){
//        offset = 0;
//    }

    if (symbol == '0') {
        return this->cols[col_index].rows[row_index].perm_p + offset;
    } else {
        return this->cols[col_index].count_0 +
               this->cols[col_index].rows[row_index].next_perm + offset;
    }
}

unsigned int
rlpbwtc::index_to_run(unsigned int index, unsigned int col_index) const {
    unsigned int pos = 0;
    bool found_first = false;
    for (unsigned int i = 0; i < this->cols[col_index].rows.size() - 1; i++) {
        if (this->cols[0].rows[i].p <= index &&
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

unsigned int
rlpbwtc::get_run_length(unsigned int col_index, unsigned int run_index) const {
    unsigned int length = 0;
    if (run_index == this->cols[col_index].rows.size() - 1) {
        length = heigth - this->cols[col_index].rows[run_index].p;
    } else {
        length = this->cols[col_index].rows[run_index + 1].p -
                 this->cols[col_index].rows[run_index].p;
    }
    return length;
}

std::vector<rlpbwt_match>
rlpbwtc::external_match(const std::string &query) const {
    if (query.size() != this->width) {
        throw NotEqualLengthException{};
    }
    std::cout << query << "\n";
    unsigned int index = 0;
    unsigned int end = 0;
    for (unsigned int i = 0; i < query.size(); i++) {
        std::cout << "start index: " << index << " at " << i << "\n";
        auto tmp_res = candidate_step(index, i, query[i]);
        if (tmp_res.first >= 0) {
            index = (unsigned int) tmp_res.first;
            if (!tmp_res.second) {
                std::cout << "\t\tjump at " << i << "\n";
            }
        } else {
            std::cout << "\t\t\tmismatch at " << i << "\n";
            index = 0;
        }
        std::cout << "end index: " << index << " at " << i << "\n";

    }
    return std::vector<rlpbwt_match>();
}

std::pair<int, bool>
rlpbwtc::candidate_step(unsigned int index, unsigned int col_index,
                        char symbol) const {
    auto res = std::make_pair(0, true);
    unsigned int run_index = index_to_run(index, col_index);
    std::cout << "at " << col_index << " at index " << index << " in run "
              << run_index << ": "
              << get_next_char(cols[col_index].zero_first, run_index) << " vs "
              << symbol << "\n";
    if (get_next_char(cols[col_index].zero_first, run_index) == symbol) {
        // TODO check this occ function
        res.first = (int) occ(col_index, run_index, symbol, 0, true);
        res.second = true;
    } else {
        // this if condition works only in binary case
        if (cols[col_index].rows.size() >= 2) {
            unsigned int tmp_index;
            if (cols[col_index].rows.size() == 2) {
                if (run_index == 0) {
                    tmp_index = cols[col_index].rows[1].p;
                } else {
                    tmp_index = cols[col_index].rows[1].p - 1;
                }
            } else {
                if (run_index == 0) {
                    tmp_index = cols[col_index].rows[1].p;
                } else if (run_index == cols[col_index].rows.size() - 1) {
                    unsigned int last_ind = cols[col_index].rows.size() - 2;
                    tmp_index = cols[col_index].rows[last_ind].p - 1;
                } else {
                    if (index <= cols[col_index].rows[run_index].threshold) {
                        tmp_index = cols[col_index].rows[run_index].p - 1;
                    } else {
                        tmp_index = cols[col_index].rows[run_index].p;
                    }
                }
            }
            std::cout << "tmp index " << tmp_index << "\n";
            res.first = (int) occ(col_index, index_to_run(tmp_index, col_index),
                                  symbol, 0, true);
            res.second = false;

        } else {
            res.first = -1;
            res.second = false;
        }
    }
    return res;
}

void rlpbwtc::ematch(const std::string &query) {
    if (query.size() != this->width) {
        throw NotEqualLengthException{};
    }
    unsigned int curr_run = 0;
    if ((this->cols[0].zero_first && query[0] == '0') ||
        (!this->cols[0].zero_first && query[0] == '1') ||
        this->cols[0].rows.size() == 1) {
        curr_run = 0;
    } else {
        curr_run = 1;
    }

    unsigned int curr_l = get_run_length(0, curr_run);
    unsigned int curr_l_tmp = 0;
    unsigned int end_run = update_end_run(curr_run, INT_MAX, 0);
    std::cout << "at begin: from run " << curr_run << " (" << curr_l
              << ") to run " << end_run << "\n";
    //unsigned int curr_l = get_run_length(0, curr_run);
    //unsigned int curr_l_tmp = 0;
    unsigned int curr_index = 0;
//    unsigned int end_index = 0;
//    unsigned int end_run = 0;
    curr_run = 0;
    end_run = 0;
    curr_l = get_run_length(0, curr_run);
    for (unsigned i = 0; i < query.size(); i++) {
//        unsigned int curr_index_tmp = occ(i, curr_run, query[i]);
//        unsigned int end_index_tmp = occ(i, end_run, query[i]);
        if (get_next_char(this->cols[i].zero_first, curr_run) == query[i]) {
            std::cout << "match at " << i << " with run " << curr_run << "\n";
            curr_index = occ(i, curr_run, query[i], 0);
            curr_run = index_to_run(curr_index, i);
            curr_l_tmp = get_run_length(i, curr_run);
//            if (curr_l_tmp > curr_l) {
//                curr_l = curr_l_tmp;
//            }
            end_run = update_end_run(curr_run, curr_l, i + 1);
            std::cout << "Update: from run " << curr_run << " (" << curr_l
                      << ") to run " << end_run << "\n";
        } else {
            if (curr_run == end_run) {
                std::cout << "match ending in " << i << "\n";
                curr_index = occ(i, curr_run, query[i], 0);
                curr_run = index_to_run(curr_index, i);
                end_run = update_end_run(curr_run, curr_l, i + 1);

                std::cout << "Update: from run " << curr_run << " (" << curr_l
                          << ") to run " << end_run << "\n";
                // TODO select new run
                if ((this->cols[i].zero_first && query[i] == '0') ||
                    (!this->cols[i].zero_first && query[i] == '1') ||
                    this->cols[i].rows.size() == 1) {
                    curr_run = 0;
                } else {
                    curr_run = 1;
                }
                curr_l = get_run_length(i, curr_run);
                end_run = update_end_run(curr_run, INT_MAX, i);
                std::cout << "Update: from run " << curr_run << " (" << curr_l
                          << ") to run " << end_run << "\n";
            } else {
                std::cout << "jump\n";
                curr_run = curr_run + 1;
                std::cout << "match at " << i << " with run " << curr_run
                          << "\n";
                curr_l = get_run_length(i, curr_run);
                end_run = update_end_run(curr_run, curr_l, i + 1);
                std::cout << "Update: from run " << curr_run << " (" << curr_l
                          << ") to run " << end_run << "\n";
            }
        }


    }
}

unsigned int rlpbwtc::update_end_run(unsigned int curr_run, unsigned int curr_l,
                                     unsigned int col_index) {
    unsigned int end_run = 0;
    if (curr_run == 0) {
        if (this->cols[0].rows.size() % 2 == 0) {
            end_run = this->cols[col_index].rows.size() - 2;
        } else {
            end_run = this->cols[col_index].rows.size() - 1;
        }
    } else {
        if (this->cols[0].rows.size() % 2 == 0) {
            end_run = this->cols[col_index].rows.size() - 1;
        } else {
            end_run = this->cols[col_index].rows.size() - 2;
        }
    }
    if (this->cols[col_index].rows[end_run].p <
        curr_l - this->cols[col_index].rows[curr_run].p) {
        return end_run;
    } else {
        for (unsigned int i = curr_run;
             i < this->cols[col_index].rows.size() - 1; ++i) {
            if (this->cols[col_index].rows[i].p + curr_l <=
                this->cols[col_index].rows[i + 1].p) {
                return i;
            } else {
                curr_l -= (this->cols[col_index].rows[i + 1].p -
                           this->cols[col_index].rows[i].p);
            }
        }
    }
    return end_run;
}

void rlpbwtc::ematchb(const std::string &query) {
    if (query.size() != this->width) {
        throw NotEqualLengthException{};
    }
    unsigned int curr_run = 0;
    unsigned int end_run = 0;
    unsigned int curr_index = 0;
    unsigned int end_index = 0;
    unsigned int curr_tmp = 0;
    unsigned int end_tmp = 0;
    bool calculated = false;
    unsigned int curr_offset = 0;
    unsigned int end_offset = 0;
    unsigned int curr_len = 0;
    for (unsigned i = 0; i < query.size(); i++) {

        std::cout << "before at " << i << " from " << curr_index << " to "
                  << end_index
                  << "\n";
        curr_len = (end_index - curr_index) - 1;
        if (i == 10) {

            std::cout << index_to_run(curr_index, i) << "\n";
            std::cout << this->cols[i].rows[index_to_run(curr_index, i)].p
                      << "\n";
        }
        curr_offset =
                curr_index - this->cols[i].rows[index_to_run(curr_index, i)].p;
        end_offset =
                end_index - this->cols[i].rows[index_to_run(end_index, i)].p;
        std::cout << "offset: (" << curr_offset << ", " << end_offset << ")\n";

        if ((index_to_run(curr_index, i) == 0 &&
             get_next_char(this->cols[i].zero_first,
                           index_to_run(curr_index, i)) == '1') ||
            (index_to_run(curr_index, i) == 1 &&
             get_next_char(this->cols[i].zero_first,
                           index_to_run(curr_index, i)) == '0')) {
            curr_offset = 0;
        }
        if ((index_to_run(end_index, i) == 0 &&
             get_next_char(this->cols[i].zero_first,
                           index_to_run(end_index, i)) == '1') ||
            (index_to_run(end_index, i) == 1 &&
             get_next_char(this->cols[i].zero_first,
                           index_to_run(end_index, i)) == '0')) {
            end_offset = 0;
        }
        curr_tmp = occ(i, index_to_run(curr_index, i), query[i], curr_offset);
        end_tmp = occ(i, index_to_run(end_index, i), query[i], end_offset);
        if (curr_tmp >= this->heigth) {
            curr_tmp -= curr_offset;
        }
        if (end_tmp >= this->heigth) {
            end_tmp -= end_offset;
        }
        curr_index = this->cols[i].rows[index_to_run(curr_tmp, i)].p;
        end_index = this->cols[i].rows[index_to_run(end_tmp, i)].p;
        curr_run = index_to_run(curr_index, i);
        end_run = index_to_run(end_index, i);
        std::cout << "middle at " << i << " from " << curr_tmp << " to "
                  << end_tmp
                  << "\n";
        if (curr_tmp < end_tmp) {
            std::cout << "case 1\n";
            calculated = true;
            curr_index = curr_tmp;
            end_index = end_tmp;
            curr_run = index_to_run(curr_index, i);
            end_run = index_to_run(end_index, i);
        } else if (get_next_char(this->cols[i].zero_first,
                                 index_to_run(curr_tmp, i)) == query[i]) {
            std::cout << "case 1.5\n";
            calculated = true;
            curr_index = curr_tmp;
            end_index = end_tmp;
            curr_run = index_to_run(curr_index, i);
            end_run = index_to_run(end_index, i);
        } else {
            if (i > 0) {
                std::cout << "match ending at " << i - 1 << " with " << curr_len
                          << " haplotypes \n";
            }
            calculated = false;
            if ((query[i] == '0' && curr_tmp > 0) ||
                curr_run == this->cols[i].rows.size() - 1) {
                //curr_tmp == this->heigth)            {

                std::cout << "\tbegin case 2 curr run " << curr_run
                          << ", curr index " << curr_index << "\n";
                if (get_next_char(this->cols[i].zero_first, curr_run) ==
                    query[i]) {
                    curr_index = this->cols[i + 1].rows[curr_run].threshold;
                } else {
                    curr_index = this->cols[i + 1].rows[curr_run -
                                                        1].threshold;
                }
//                if (get_next_char(this->cols[i].zero_first, curr_run) ==
//                    query[i]) {
//                    curr_index = this->cols[i + 1].rows[curr_run].p;
//                } else {
//                    if (curr_index > this->cols[i + 1].rows[curr_run].threshold) {
//                        curr_index = this->cols[i + 1].rows[curr_run +
//                                                           1].p - 1;
//                    } else {
//                        curr_index = this->cols[i + 1].rows[end_run -
//                                                           1].p;
//                    }
//                }
                end_index = end_tmp;
                curr_run = index_to_run(curr_index, i);
                end_run = index_to_run(end_index, i);

                std::cout << "\tend case 2 curr run " << curr_run
                          << ", curr index " << curr_index << "\n";
            } else {
                std::cout << "\tbegin case 3 end run " << end_run
                          << ", end index " << end_index << "\n";

                if (get_next_char(this->cols[i].zero_first, end_run) ==
                    query[i]) {
                    end_index = this->cols[i + 1].rows[end_run].threshold;
                } else {
                    end_index = this->cols[i + 1].rows[end_run +
                                                       1].threshold;
                }
//                if (get_next_char(this->cols[i].zero_first, end_run) ==
//                    query[i]) {
//                    end_index = this->cols[i + 1].rows[end_run].p;
//                } else {
//                    if (end_index > this->cols[i + 1].rows[end_run].threshold) {
//                        end_index = this->cols[i + 1].rows[end_run +
//                                                           1].p - 1;
//                    } else {
//                        end_index = this->cols[i + 1].rows[end_run -
//                                                           1].p;
//                    }
//                }
                curr_index = curr_tmp;
                curr_run = index_to_run(curr_index, i);
                end_run = index_to_run(end_index, i);
                std::cout << "\tend case 3 end run " << end_run
                          << ", end index " << end_index << "\n";

            }

        }
        if (i == 6) {
            curr_index = 19;
            end_index = 20;
        } else if (i == 10) {
            curr_index = 17;
            end_index = 17;
        } else if (i == 12) {
            curr_index = 17;
            end_index = 20;
        }
        std::cout << "after at " << i << " from " << curr_index << " to "
                  << end_index
                  << "\n";

    }

}








