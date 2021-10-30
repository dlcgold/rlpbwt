//
// Created by dlcgold on 29/10/21.
//


#include <iostream>
#include "../include/utils.h"


void update(std::string column, std::vector<unsigned int> &pref,
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

pbwt_column
build_column(std::string column, std::vector<unsigned int> pref) {
    unsigned int count0 = 0;
    unsigned int count1 = 0;
    //int threshold;
    // int lcs;
    unsigned int height = pref.size();
    bool start = true;
    for (unsigned int i = 0; i < height; i++) {
        if (i == 0 && column[pref[i]] == '1') {
            start = false;
        }
        if (column[i] == '0') {
            count0++;
        }
    }

    std::vector<pbwt_rlrow> rows;
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
        }
        if ((i == height - 1) || (column[pref[i]] != column[pref[i + 1]])) {
            rows.emplace_back(p_tmp, perm_tmp, 0);
        }
    }
    return pbwt_column(start, rows);
}

void build_next_perm(std::vector<pbwt_column> &cols, unsigned int index) {
    for (unsigned int i = 0; i < cols[index - 1].rows.size(); i++) {
        bool found = false;
        for (unsigned int j = 0; j < cols[index].rows.size() - 1; j++) {
            if (cols[index].rows[j].p <= cols[index - 1].rows[i].perm_p &&
                cols[index - 1].rows[i].perm_p < cols[index].rows[j + 1].p) {
                cols[index - 1].rows[i].next_perm = j;
                found = true;
                break;
            }
        }
        if (!found) {
            cols[index - 1].rows[i].next_perm =
                    cols[index].rows.size() - 1;
        }
    }

}

char get_next_char(bool zero_first, unsigned int pos) {
    if (zero_first) {
        if (pos % 2 == 0) {
            return '0';
        } else {
            return '1';
        }
    } else {
        if (pos % 2 == 0) {
            return '1';
        } else {
            return '0';
        }
    }
}
