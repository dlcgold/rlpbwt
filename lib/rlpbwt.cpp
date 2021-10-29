//
// Created by dlcgold on 28/10/21.
//
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include "../include/rlpbwt.h"
#include "../include/pbwt_column.h"
#include "../include/pbwt_rlrow.h"
#include "../include/utils.h"

rlpbwt::rlpbwt(char *filename) {
    std::ifstream input_matrix(filename);
    if (input_matrix.is_open()) {
        std::string column;
        getline(input_matrix, column);
        column.erase(std::remove(column.begin(), column.end(), ' '),
                     column.end());
        const unsigned int height = column.size();
        unsigned int width = std::count(
                std::istreambuf_iterator<char>(input_matrix),
                std::istreambuf_iterator<char>(), '\n') + 1;
        std::vector<pbwt_column> cols(width);
        input_matrix.clear();
        input_matrix.seekg(0, std::ios::beg);
        std::vector<int> pref(height);
        std::vector<int> div(height);
        for (int i = 0; i < height; i++) {
            pref[i] = i;
            div[i] = 0;
        }
        unsigned int count = 0;
        while (getline(input_matrix, column)) {
            column.erase(std::remove(column.begin(), column.end(), ' '),
                         column.end());
            auto col = build_column(column, pref, div);
            cols[count] = col;
            if (count != 0){
                build_next_perm(cols, count);
            }
            update(column, pref, div);
            count++;
        }
        this->cols = cols;
    } else {
        throw "File not found";
    }
}

rlpbwt::~rlpbwt() {}

std::string rlpbwt::search_row(unsigned int row_index) {
    int tmp = 0;
    int start = 0;
    int pos = 0;
    std::string row = "";
    for (int i = 0; i < this->cols[0].rows.size() - 1; i++) {
        if (this->cols[0].rows[i].p < row_index &&
            row_index <= this->cols[0].rows[i + 1].p) {
            start = this->cols[0].rows[i + 1].perm_p;
            pos = i + 1;
            if (this->cols[0].zero_first) {
                if (i + 1 % 2 == 0) {
                    row.push_back('0');
                } else {
                    row.push_back('1');
                }
            } else {
                if (i + 1 % 2 == 0) {
                    row.push_back('1');
                } else {
                    row.push_back('0');
                }
            }
            break;
        }
    }
    int end = this->cols[0].rows[pos].lf_mapping(row_index);
    std::cout << "start: " << start << " " << "end: " << end << "\n";

    for (int i = 1; i < this->cols.size(); i++) {
        std::cout << "start: " << start << " " << "end: " << end << "\n";
        for (int j = start; j < this->cols[i].rows.size() - 1; j++) {
            if ((this->cols[i].rows[j].p < end &&
                 end <= this->cols[i].rows[j + 1].p) ||
                j == this->cols[i].rows.size() - 1) {
                if (this->cols[i].zero_first) {
                    if (j + 1 % 2 == 0) {
                        row.push_back('0');
                    } else {
                        row.push_back('1');
                    }
                } else {
                    if (j + 1 % 2 == 0) {
                        row.push_back('1');
                    } else {
                        row.push_back('0');
                    }
                }
                end = this->cols[i].rows[start].lf_mapping(end);
                start = this->cols[i].rows[j + 1].perm_p;

                break;
            }
        }
    }
    return row;
}
