//
// Created by dlcgold on 29/10/21.
//


#include <iostream>
#include "../include/utils.h"


void update(std::string column, std::vector<int> &pref, std::vector<int> &div) {
    unsigned int height = pref.size();
    std::vector<int> new_pref(height);
    std::vector<int> new_div(height);
    int count0 = 0;
    int lcs = INT_MAX;

    for (int i = 0; i < height; i++) {
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

    for (int i = 0; i < height; i++) {
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
build_column(std::string column, std::vector<int> pref, std::vector<int> div) {
    int count0 = 0;
    int count1 = 0;
    int threshold;
    int lcs;
    unsigned int height = pref.size();
    bool start = true;
    for (int i = 0; i < height; i++) {
        if (i == 0 && column[i] == '1') {
            start = false;
        }
        if (column[i] == '0') {
            count0++;
        }
    }

    std::vector<pbwt_rlrow> rows;
    int p_tmp = 0;
    int perm_tmp = 0;
    int next_perm_tmp = 0;
    int run0 = 0;
    int run1 = 0;
    for (int i = 0; i < height; i++) {
        if (column[pref[i]] == '1') {
            count1++;
        }

        if ((i == 0) || (column[pref[i]] != column[pref[i - 1]])) {
            //fprintf(filePtr, "\t%i\t%i", i, column[pref[i]]);
            p_tmp = i;
            if (column[pref[i]] == '0') {
                run0++;
            } else {
                run1++;
            }
            if (column[pref[i]] == '0') {
                //fprintf(filePtr, "\t%i", i - count1);
                perm_tmp = i - count1;
            } else {
                //fprintf(filePtr, "\t%i", count0 + count1 - 1);
                perm_tmp = count0 + count1 - 1;
            }
            //next_perm_tmp = pref[i];
            //fprintf(filePtr, "\t%i", pref[i]);
            threshold = i;
            lcs = div[i];
        }

        if (div[i] < lcs) {
            threshold = i;
            lcs = div[i];
        }

        if ((i == height - 1) || (column[pref[i]] != column[pref[i + 1]])) {
            //fprintf(filePtr, "\t%i\t%i\n", pref[i], threshold);
            rows.emplace_back(p_tmp, perm_tmp, next_perm_tmp);
        }
    }
    /*for(int i = 0; i < rows.size(); i++){
        if(!start && ){
            rows[i].next_perm=run0+
        }
    }*/
    return pbwt_column(start, rows);
}

void build_next_perm(std::vector<pbwt_column> &cols, int index) {
    for (int i = 0; i < cols[index - 1].rows.size(); i++) {
        for (int j = 0; j < cols[index].rows.size() - 1; j++) {
            if (cols[index].rows[j].p < cols[index - 1].rows[i].perm_p &&
                cols[index - 1].rows[i].perm_p <= cols[index].rows[j + 1].p) {
                cols[index - 1].rows[i].next_perm = j + 1;
            }
        }
    }
}
