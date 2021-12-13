//
// Created by dlcgold on 28/10/21.
//

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include "../include/rlpbwt.h"
#include "../include/utils.h"
#include "../include/exceptions.h"

rlpbwt::rlpbwt(const char *filename, bool verbose) {
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
            auto col = rlpbwt::build_column(column, pref, div);
            tmp_cols[count] = col;
            if (count != 0) {
                rlpbwt::build_next_perm(tmp_cols[count - 1], tmp_cols[count]);
            }
            //rlpbwt::update_old(column, pref, div, count);
            rlpbwt::update(column, pref, div);
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

rlpbwt::~rlpbwt() = default;

std::string rlpbwt::search_row(unsigned int row_index, bool verbose) {
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


void rlpbwt::build_next_perm(rlpbwt_column &prev, rlpbwt_column &curr) {
    for (auto &row: prev.rows) {
        bool found = false;
        for (unsigned int j = 0; j < curr.rows.size() - 1; j++) {
            if (curr.rows[j].p <= row.perm_p &&
                row.perm_p < curr.rows[j + 1].p) {
                row.next_perm = j;
                found = true;
                break;
            }
        }
        if (!found) {
            row.next_perm = curr.rows.size() - 1;
        }
    }
}

rlpbwt_column
rlpbwt::build_column(std::string &column, std::vector<unsigned int> &pref,
                     std::vector<unsigned int> &div) {
    unsigned int height = pref.size();
    unsigned int count0 = 0;
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
    for (unsigned int i = 0; i < height; i++) {
        if (column[pref[i]] == '1') {
            count1++;
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

        if ((i == height - 1) || (column[pref[i]] != column[pref[i + 1]])) {
            rows.emplace_back(p_tmp, perm_tmp, 0, threshold);
        }
    }
    return {start, rows, count0};
}

__attribute__((unused)) void
rlpbwt::update_old(std::string &column, std::vector<unsigned int> &pref,
                   std::vector<unsigned int> &div, unsigned int k) {
    unsigned int height = pref.size();
    std::vector<unsigned int> new_pref(height);
    std::vector<unsigned int> new_div(height);
    unsigned int u = 0;
    unsigned int v = 0;
    unsigned int p = k + 1;
    unsigned int q = k + 1;
    std::vector<unsigned int> a(height);
    std::vector<unsigned int> b(height);
    std::vector<unsigned int> d(height);
    std::vector<unsigned int> e(height);
    for (unsigned int i = 0; i < height; i++) {
        if (div[i] > p) {
            p = div[i];
        }
        if (div[i] > q) {
            q = div[i];
        }
        if (column[pref[i]] == '0') {
            a[u] = pref[i];
            d[u] = p;
            u++;
            p = 0;
        } else {
            b[v] = pref[i];
            e[v] = q;
            v++;
            q = 0;
        }
    }
    unsigned int offset = 0;
    for (unsigned int l = 0; l < u; l++) {
        new_pref[offset] = a[l];
        new_div[offset] = k - d[l];
        offset++;
    }
    for (unsigned int l = 0; l < v; l++) {
        new_pref[offset] = b[l];
        new_div[offset] = k - e[l];
        offset++;
    }
    new_div[0] = 0;
    new_div[u] = 0;
    div = new_div;
    pref = new_pref;

}

void rlpbwt::update(std::string &column, std::vector<unsigned int> &pref,
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
rlpbwt::next_run(unsigned int col_index, unsigned int start, unsigned int end,
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
rlpbwt::occ(unsigned int col_index, unsigned int row_index, char symbol) const {
    unsigned int u = 0;
    unsigned int v = 0;
    unsigned int c = this->cols[col_index].count_0;
    for (unsigned int i = 0; i < row_index; i++) {
        std::cout << this->cols[col_index].rows[i + 1].p << " - "
                  << this->cols[col_index].rows[i].p << "\n";
        if (i % 2 == 0) {
            if (row_index != this->cols[col_index].rows.size() &&
                i != row_index - 1) {
                u += this->cols[col_index].rows[i + 1].p -
                     this->cols[col_index].rows[i].p;
            } else {
                u += this->heigth - this->cols[col_index].rows[i].p;
            }
        } else {
            if (row_index != this->cols[col_index].rows.size() &&
                i != row_index - 1) {
                v += this->cols[col_index].rows[i + 1].p -
                     this->cols[col_index].rows[i].p;
            } else {
                v += this->heigth - this->cols[col_index].rows[i].p;
            }
        }
    }
    if (!this->cols[col_index].zero_first) {
        std::swap(u, v);
    }

    std::cout << symbol << " c: " << c << ", u: " << u << ", v: " << v << "\n";
    if (symbol == '0') {
        return u;
    } else {
        return c + v;
    }
}


std::vector<rlpbwt_match>
rlpbwt::external_match(const std::string &query) const {
    if (query.size() != this->width) {
        throw NotEqualLengthException{};
    }
    std::cout << query << "\n";
    unsigned int curr_r = 0;
    unsigned int curr_l = 0;
    unsigned int curr_bit = 0;
    unsigned int curr_start = 0;

    // search longest run of first char in first column
    unsigned int max_r = 0;
    unsigned int ind_max = 0;
    //std::cout << this->cols[0].zero_first << " " << query[0] << "\n";
    if ((this->cols[0].zero_first && query[0] == '0') ||
        (!this->cols[0].zero_first && query[0] == '1')) {
        for (unsigned int i = 0; i < this->cols[0].rows.size(); i++) {
            unsigned int run_l = 0;
            if (i != this->cols[0].rows.size() - 1) {
                run_l = this->cols[0].rows[i + 1].p - this->cols[0].rows[i].p;
            } else {
                run_l = this->heigth - this->cols[0].rows[i].p;
            }
            //std::cout << run_l << ": " << this->cols[0].rows[i] << "\n";
            if (max_r < run_l) {
                max_r = run_l;
                ind_max = i;
            }
            if (i <= this->cols[0].rows.size() - 2) {
                i++;
            } else {
                break;
            }
        }
    } else {
        for (unsigned int i = 1; i < this->cols[0].rows.size(); i++) {
            unsigned int run_l = 0;
            if (i != this->cols[0].rows.size() - 1) {
                run_l = this->cols[0].rows[i + 1].p - this->cols[0].rows[i].p;
            } else {
                run_l = this->heigth - this->cols[0].rows[i].p;
            }
            //std::cout << run_l << ": " << this->cols[0].rows[i] << "\n";
            if (max_r < run_l) {
                max_r = run_l;
                ind_max = i;
            }
            if (i <= this->cols[0].rows.size() - 2) {
                i++;
            } else {
                break;
            }
        }
    }
    std::cout << ind_max << "\n";
    //unsigned int start = occ(0, ind_max, query[1]);
    //unsigned int end = occ(0, ind_max + 1, query[1]);

    //unsigned int start = 0;
    //unsigned int end = 0;
    //std::cout << "(" << start << ", " << end << ")\n";
    curr_r = ind_max;
    unsigned int curr_row = this->cols[0].rows[curr_r].p;
    unsigned int start = this->cols[0].rows[curr_r].next_perm;
    unsigned int end = this->cols[0].rows[curr_r].lf_mapping(curr_row);
    bool verbose = true;
    char curr_c;
    unsigned int curr_e = 1;
    curr_l = get_run_length(0, curr_r);
    std::cout << "run: " << curr_r << ", row " << curr_row << "\n";
    std::cout << "(" << start << ", " << end << ")\n";
    std::cout << "current_length: " << curr_l << "\n";
    for (unsigned int o = 1; o < this->cols.size(); o++) {
        bool check = true;
        unsigned int i = curr_e;
        while (check) {
            if (start == this->cols[i].rows.size() - 1) {
                curr_c = get_next_char(this->cols[i].zero_first, start);
                end = this->cols[i].rows[start].lf_mapping(end);
                start = this->cols[i].rows[start].next_perm;
                if (verbose) {
                    std::cout << "column " << i << ": " << start << ", " << end
                              << "\n";
                }
                curr_r = start;
            } else {
                bool found = false;
                for (unsigned int j = start; j < cols[i].rows.size() - 1; j++) {
                    if (cols[i].rows[j].p <= end &&
                        end < cols[i].rows[j + 1].p) {
                        curr_c = get_next_char(this->cols[i].zero_first, j);
                        found = true;
                        end = this->cols[i].rows[j].lf_mapping(end);
                        start = this->cols[i].rows[j].next_perm;
                        if (verbose) {
                            std::cout << "column " << i << ": " << start << ", "
                                      << end
                                      << "\n";
                        }
                        curr_r = j;
                        break;
                    }
                }
                if (!found) {
                    unsigned int endrow = this->cols[i].rows.size() - 1;
                    curr_c = get_next_char(this->cols[i].zero_first, endrow);
                    end = this->cols[i].rows[endrow].lf_mapping(end);
                    start = this->cols[i].rows[endrow].next_perm;
                    if (verbose) {
                        std::cout << "column " << i << ": " << start << ", "
                                  << end
                                  << "\n";
                    }
                    curr_r = endrow;
                }
            }
            curr_row = this->cols[i].rows[curr_r].p;
            unsigned int tmp_length = get_run_length(i, curr_r);
            std::cout << "i: " << i << ", o: " << o << "\n";
            std::cout << "current char " << curr_c << " vs " << query[o]
                      << "\n";
            std::cout << "current run " << curr_r << " current row " << curr_row
                      << "\n";
            std::cout << "prev length: " << curr_l << ", current length: "
                      << tmp_length << "\n";
            if (o == this->width - 1) {
                std::cout << "match ending in " << o << "\n";
                break;
            }
            if (curr_c != query[o]) {
                if (curr_l > tmp_length &&
                    curr_r != this->cols[i].rows.size()) {
                    curr_r = curr_r + 1;
                    tmp_length = get_run_length(i, curr_r);
                    curr_row = this->cols[i].rows[curr_r].p;
                } else {
                    check = false;
                    if (o != this->width) {
                        std::cout << "match ending in " << o << "\n";
                    } else {
                        std::cout << "match ending in " << o - 1 << "\n";
                    }
                }
            }
            curr_l = tmp_length;
            i++;
            o++;
        }
        if (o == this->width) {
            break;
        }
        std::cout << "old row: " << curr_row << "\n";
        std::cout << "old run: " << curr_r << "\n";
        if (o <= this->cols[o].rows[curr_r].threshold) {
            if (curr_row - 1 >= 0) {
                curr_row = curr_row - 1;
            } else {
                curr_row = 0;
            }
        } else {
            curr_row = 0;
        }
        std::cout << "new row: " << curr_row << "\n";
        curr_r = index_to_run(curr_row);
        start = this->cols[o].rows[curr_r].next_perm;
        end = this->cols[o].rows[curr_r].lf_mapping(curr_row);
        curr_e = o;
    }
    return std::vector<rlpbwt_match>();
}

std::vector<unsigned int>
rlpbwt::update_external(unsigned int index, unsigned int e, unsigned int f,
                        unsigned int g, const std::string &query) {

}

unsigned int rlpbwt::index_to_run(unsigned int index) const {
    unsigned int pos = 0;
    bool found_first = false;
    for (unsigned int i = 0; i < this->cols[0].rows.size() - 1; i++) {
        if (this->cols[0].rows[i].p <= index &&
            index < this->cols[0].rows[i + 1].p) {
            pos = i;
            found_first = true;
            break;
        }
    }
    if (!found_first) {
        pos = this->cols[0].rows.size() - 1;
    }
    return pos;
}

unsigned int
rlpbwt::get_run_length(unsigned int col_index, unsigned int run_index) const {
    unsigned int length = 0;
    if (run_index == this->cols[col_index].rows.size() - 1) {
        length = heigth - this->cols[col_index].rows[run_index].p;
    } else {
        length = this->cols[col_index].rows[run_index + 1].p -
                 this->cols[col_index].rows[run_index].p;
    }
    return length;
}




