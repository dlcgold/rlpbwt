//
// Created by dlcgold on 20/01/22.
//

#include "../include/birlpbwt.h"

birlpbwt::birlpbwt(const char *filename, bool verbose) {
    std::ifstream input_matrix(filename);
    std::list<std::string> revrows;
    if (input_matrix.is_open()) {
        std::string new_column;
        getline(input_matrix, new_column);
        new_column.erase(std::remove(new_column.begin(), new_column.end(), ' '),
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
                std::cout << new_column << "\n";

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
            revrows.push_front(new_column);
            auto col = rlpbwt::build_column(new_column, pref, div);
            sdsl::util::bit_compress(div);
            col.lcp = div;
            tmp_cols[count] = col;
            rlpbwt::update(new_column, pref, div);
            count++;
        }
        auto col = rlpbwt::build_column(new_column, pref, div);
        sdsl::util::bit_compress(div);
        col.lcp = div;
        tmp_cols.push_back(col);
        this->frlpbwt.cols = tmp_cols;
        this->frlpbwt.width = tmp_width;
        this->frlpbwt.heigth = tmp_height;
        input_matrix.close();

        std::vector<column> tmp_colsb(tmp_width);
        input_matrix.clear();
        input_matrix.seekg(0, std::ios::beg);
        for (unsigned int i = 0; i < tmp_height; i++) {
            pref[i] = i;
            div[i] = 0;
        }
        count = 0;
        for(auto& new_columnb : revrows){
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
            auto colf = rlpbwt::build_column(new_columnb, pref, div);
            sdsl::util::bit_compress(div);
            colf.lcp = div;
            tmp_colsb[count] = col;
            rlpbwt::update(new_columnb, pref, div);
            count++;
        }
        this->brlpbwt.cols = tmp_colsb;
        this->brlpbwt.width = tmp_width;
        this->brlpbwt.heigth = tmp_height;
    } else {
        throw FileNotFoundException{};
    }
}

std::vector<match>
birlpbwt::external_match(const std::string &query, bool verbose) {
    std::vector<match> fmatch;
    std::vector<match> bmatch;
    std::vector<match> matches;

    return matches;
}
