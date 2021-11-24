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

rlpbwt::rlpbwt(const char *filename) {
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
            column.erase(std::remove(column.begin(), column.end(), ' '),
                         column.end());
            auto col = rlpbwt::build_column(column, pref, div);
            tmp_cols[count] = col;
            if (count != 0) {
                rlpbwt::build_next_perm(tmp_cols[count - 1], tmp_cols[count]);
            }
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

std::string rlpbwt::search_row(unsigned int row_index) {
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
        if (start == cols[i].rows.size() - 1) {
            row.push_back(get_next_char(this->cols[i].zero_first, start));
            end = this->cols[i].rows[start].lf_mapping(end);
            start = this->cols[i].rows[start].next_perm;
        } else {
            bool found = false;
            for (unsigned int j = start; j < cols[i].rows.size() - 1; j++) {
                if (cols[i].rows[j].p <= end && end < cols[i].rows[j + 1].p) {
                    row.push_back(get_next_char(this->cols[i].zero_first, j));
                    found = true;
                    end = this->cols[i].rows[j].lf_mapping(end);
                    start = this->cols[i].rows[j].next_perm;
                    break;
                }
            }
            if (!found) {
                unsigned int endrow = this->cols[i].rows.size() - 1;
                row.push_back(get_next_char(this->cols[i].zero_first, endrow));
                end = this->cols[i].rows[endrow].lf_mapping(end);
                start = this->cols[i].rows[endrow].next_perm;
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
    return {start, rows};
}

void rlpbwt::update(std::string &column, std::vector<unsigned int> &pref,
                    std::vector<unsigned int> &div) {
    unsigned int height = pref.size();
    std::vector<unsigned int> new_pref(height);
    std::vector<unsigned int> new_div(height);
    unsigned int count0 = 0;
    unsigned int lcs = INT_MAX;

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
    lcs = INT_MAX;

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

std::vector<rlpbwt_match> rlpbwt::external_match(const std::string &query) {
    if(query.size() != this->width){
        throw NotEqualLengthException{};
    }
    std::vector<unsigned int> f_arr(this->width + 1);
    std::vector<unsigned int> g_arr(this->width + 1);
    std::vector<unsigned int> e_arr(this->width + 1);

    return std::vector<rlpbwt_match>();
}

std::vector<unsigned int>
rlpbwt::update_external(unsigned int index, unsigned int e, unsigned int f,
                        unsigned int g, const std::string &query) {

}
