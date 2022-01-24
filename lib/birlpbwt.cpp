//
// Created by dlcgold on 20/01/22.
//

#include "../include/birlpbwt.h"

birlpbwt::birlpbwt(const char *filename, bool vcf, bool verbose) {
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
        std::vector<unsigned int> pref(tmp_height);
        std::vector<unsigned int> div(tmp_height);
        for (int i = 0; i < tmp_height; i++) {
            pref[i] = i;
            div[i] = 0;
        }
        std::string new_column;
        std::vector<column> tmp_cols(tmp_width);
        std::list<std::string> revrows;

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
            while (getline(ss, line, '\t')) {
                new_column.push_back(line[0]);
                new_column.push_back(line[2]);
            }
            revrows.push_front(new_column);
            auto col = rlpbwt::build_column(new_column, pref, div);
            tmp_cols[k] = col;
            rlpbwt::update(new_column, pref, div);
        }
        auto col = rlpbwt::build_column(new_column, pref, div);
        tmp_cols.push_back(col);
        this->frlpbwt.cols = tmp_cols;
        this->frlpbwt.width = tmp_width;
        this->frlpbwt.heigth = tmp_height;
        std::vector<column> tmp_colsb(tmp_width);
        for (int i = 0; i < tmp_height; i++) {
            pref[i] = i;
            div[i] = 0;
        }
        unsigned int count = 0;
        for (auto &new_columnb: revrows) {
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
            auto colb = rlpbwt::build_column(new_columnb, pref, div);
            tmp_colsb[count] = colb;
            rlpbwt::update(new_columnb, pref, div);
            count++;
        }
        auto colb = rlpbwt::build_column(revrows.back(), pref, div);
        tmp_colsb.push_back(colb);
        this->brlpbwt.cols = tmp_colsb;
        this->brlpbwt.width = tmp_width;
        this->brlpbwt.heigth = tmp_height;
    } else {
        std::ifstream input_matrix(filename);
        std::list<std::string> revrows;
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
            std::vector<unsigned int> div(tmp_height);
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
                tmp_cols[count] = col;
                rlpbwt::update(new_column, pref, div);
                count++;
            }
            auto col = rlpbwt::build_column(new_column, pref, div);
            tmp_cols.push_back(col);
            this->frlpbwt.cols = tmp_cols;
            this->frlpbwt.width = tmp_width;
            this->frlpbwt.heigth = tmp_height;
            input_matrix.close();

            std::vector<column> tmp_colsb(tmp_width);
            for (unsigned int i = 0; i < tmp_height; i++) {
                pref[i] = i;
                div[i] = 0;
            }
            count = 0;
            for (auto &new_columnb: revrows) {
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
                auto colb = rlpbwt::build_column(new_columnb, pref, div);
                tmp_colsb[count] = colb;
                rlpbwt::update(new_columnb, pref, div);
                count++;
            }
            auto colb = rlpbwt::build_column(revrows.back(), pref, div);
            tmp_colsb.push_back(colb);
            this->brlpbwt.cols = tmp_colsb;
            this->brlpbwt.width = tmp_width;
            this->brlpbwt.heigth = tmp_height;
        } else {
            throw FileNotFoundException{};
        }
    }
}


std::vector<match>
birlpbwt::external_match(const std::string &query, unsigned int min_len,
                         bool verbose) {
    std::vector<match> matches;

    std::vector<match> fm = this->frlpbwt.end_external_match(query, true,
                                                             verbose);
    if (verbose) {
        std::cout << "forward matches:\n";
        for (const auto &m: fm) {
            std::cout << m << "\n";
        }
    }
    std::string query_rev(query.rbegin(), query.rend());
    std::vector<match> bm = this->brlpbwt.end_external_match(query_rev, false,
                                                             verbose);

    if (verbose) {
        std::cout << "backward matches:\n";
        for (const auto &m: bm) {
            std::cout << m << "\n";
        }
    }

    // TODO this not work
    unsigned int begin = 0;
    unsigned int end = 0;
    unsigned int nhaplo = 0;
    for (unsigned int i = 0; i < fm.size(); ++i) {
        begin = std::min(fm[i].begin, bm[i].begin);
        end = std::max(fm[i].end, bm[i].end);
        nhaplo = std::min(fm[i].nhaplo, bm[i].nhaplo);
        if (end - begin != 0 && end - begin >= min_len && nhaplo > 0) {
            matches.emplace_back(begin, end, nhaplo);
        }
    }
    return matches;
}


void birlpbwt::external_match_vcf(const char *filename, unsigned int min_len,
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
    if (tmp_width != this->frlpbwt.width) {
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
    for (unsigned int k = 0; k < this->frlpbwt.width; k++) {
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
    std::cout << queries_tmp.size() << " " << queries_tmp[0].size() << "\n";
    std::ofstream out;
    for (unsigned int i = 0; i < queries_tmp[0].size(); i++) {
        for (auto & j : queries_tmp) {
            tmpq.push_back(j[i]);
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
            std::cout << "matches with " << count << " " << qIDs[count] << "\n";
            for (const auto &m: matches) {
                std::cout << m << "\n";
            }
        }
        count++;
    }

}

void birlpbwt::print() {
    std::cout << "\tforward\n";
    this->frlpbwt.print();
    std::cout << "\t--------------------------------------\n";
    std::cout << "\ttackward\n";
    this->brlpbwt.print();
}

