//
// Created by dlcgold on 02/02/22.
//

#include <list>
#include "../include/rlpbwt_thr.h"


rlpbwt_thr::rlpbwt_thr(const char *filename, unsigned int w, unsigned int h,
                       bool verbose) {
    std::ifstream input_matrix(filename);
    if (input_matrix.is_open()) {
        std::string header1;
        std::string header2;
        std::string line;
        std::string garbage;
        std::string new_column;
        /*
        getline(input_matrix, header1);
        getline(input_matrix, header2);
        getline(input_matrix, line);
        std::istringstream is(line);
        is >> garbage >> garbage >> garbage >> garbage >> new_column;
        unsigned int tmp_height = new_column.size();
        unsigned int tmp_width = std::count(
                std::istreambuf_iterator<char>(input_matrix),
                std::istreambuf_iterator<char>(), '\n') + 1;
        if(tmp_width >= w){
            tmp_width = w;
        }else{
            std::cout << "warn: set width to max\n";
        }
        if(tmp_height >= h){
            tmp_height = h;
        }else{
            std::cout << "warn: set height to max\n";
        }
         input_matrix.clear();
        input_matrix.seekg(0, std::ios::beg);
         */
        unsigned int tmp_width = w;
        unsigned int tmp_height = h;
        this->cols = std::vector<column_thr>(tmp_width + 1);
        std::vector<unsigned int> pref(tmp_height);
        sdsl::int_vector<> div(tmp_height);
        for (unsigned int i = 0; i < tmp_height; i++) {
            pref[i] = i;
            div[i] = 0;
        }

        sdsl::bit_vector bv(tmp_height * tmp_width, 0);
        this->panelbv = panel_ra(tmp_height, tmp_width);
        unsigned int count = 0;
        std::string last_col;
        getline(input_matrix, line);
        getline(input_matrix, line);
        while (getline(input_matrix, line) && !line.empty()) {
            std::istringstream is_col(line);
            is_col >> garbage >> garbage >> garbage >> garbage >> new_column;
            if (verbose) {
                std::cout << "\nnew_column " << count << "\n";
                std::cout << new_column << "\n" << this->cols[count].runs
                          << "\n"
                          << this->cols[count].u << "\n"
                          << this->cols[count].v
                          << "\n-------------------------------\n";
            }
            auto col = rlpbwt_thr::build_column(new_column, pref, div);
            for (unsigned int k = 0; k < tmp_height; k++) {
                if (new_column[k] != '0') {
                    //this->panelbv.panel[count + (k * tmp_width)] = true;
                    this->panelbv.panel[count][k] = true;
                }
            }
            this->cols[count] = col;
            this->cols[count].rank_runs = sdsl::sd_vector<>::rank_1_type(
                    &this->cols[count].runs);
            this->cols[count].select_runs = sdsl::sd_vector<>::select_1_type(
                    &this->cols[count].runs);
            this->cols[count].rank_u = sdsl::sd_vector<>::rank_1_type(
                    &this->cols[count].u);
            this->cols[count].select_u = sdsl::sd_vector<>::select_1_type(
                    &this->cols[count].u);
            this->cols[count].rank_v = sdsl::sd_vector<>::rank_1_type(
                    &this->cols[count].v);
            this->cols[count].select_v = sdsl::sd_vector<>::select_1_type(
                    &this->cols[count].v);
            this->cols[count].rank_thr = sdsl::sd_vector<>::rank_1_type(
                    &this->cols[count].thr);
            this->cols[count].select_thr = sdsl::sd_vector<>::select_1_type(
                    &this->cols[count].thr);
            rlpbwt_thr::update(new_column, pref, div);
            last_col = new_column;
            if (count == tmp_width) {
                break;
            }
            count++;

        }

        auto col = rlpbwt_thr::build_column(last_col, pref, div);

        this->cols[count] = col;
        this->cols[count].rank_runs = sdsl::sd_vector<>::rank_1_type(
                &this->cols[count].runs);
        this->cols[count].select_runs = sdsl::sd_vector<>::select_1_type(
                &this->cols[count].runs);
        this->cols[count].rank_u = sdsl::sd_vector<>::rank_1_type(
                &this->cols[count].u);
        this->cols[count].select_u = sdsl::sd_vector<>::select_1_type(
                &this->cols[count].u);
        this->cols[count].rank_v = sdsl::sd_vector<>::rank_1_type(
                &this->cols[count].v);
        this->cols[count].select_v = sdsl::sd_vector<>::select_1_type(
                &this->cols[count].v);
        this->cols[count].rank_thr = sdsl::sd_vector<>::rank_1_type(
                &this->cols[count].thr);
        this->cols[count].select_thr = sdsl::sd_vector<>::select_1_type(
                &this->cols[count].thr);
        input_matrix.close();
    } else {
        throw FileNotFoundException{};
    }
}

rlpbwt_thr::rlpbwt_thr(const char *filename, bool vcf, bool verbose) {
    if (vcf) {
        // read vcf file from https://github.com/ZhiGroup/Syllable-PBWT
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
        // initialize prefix and divergence arrays (latter as sdsl int_vector)
        std::vector<unsigned int> pref(tmp_height);
        sdsl::int_vector<> div(tmp_height);
        for (int i = 0; i < tmp_height; i++) {
            pref[i] = i;
            div[i] = 0;
        }
        std::string new_column;
        this->panelbv = panel_ra(tmp_height, tmp_width);
        // initialize vector for the column
        this->cols = std::vector<column_thr>(tmp_width + 1);
        int k = 0;
        for (k = 0; k < tmp_width; k++) {
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

            // read column as in https://github.com/ZhiGroup/Syllable-PBWT
            new_column.clear();
            std::getline(input, line);
            ss = std::stringstream(line);
            for (int i = 0; i < 9; i++) {
                getline(ss, line, '\t');
            }
            //int index = 0;
            while (getline(ss, line, '\t')) {
                new_column.push_back(line[0]);
                new_column.push_back(line[2]);
            }

            // create column with the bitvectors (but not rank/select) and add
            // to columns vector
            auto col = rlpbwt_thr::build_column(new_column, pref, div);
            for (unsigned int j = 0; j < new_column.size(); j++) {
                if (new_column[j] != '0') {
                    //this->panelbv.panel[k + (j * tmp_width)] = true;
                    this->panelbv.panel[k][j] = true;

                }
            }
            this->cols[k] = col;
            // create rank/select (here because of some problems using sdsl)
            this->cols[k].rank_runs = sdsl::sd_vector<>::rank_1_type(
                    &this->cols[k].runs);
            this->cols[k].select_runs = sdsl::sd_vector<>::select_1_type(
                    &this->cols[k].runs);
            this->cols[k].rank_u = sdsl::sd_vector<>::rank_1_type(
                    &this->cols[k].u);
            this->cols[k].select_u = sdsl::sd_vector<>::select_1_type(
                    &this->cols[k].u);
            this->cols[k].rank_v = sdsl::sd_vector<>::rank_1_type(
                    &this->cols[k].v);
            this->cols[k].select_v = sdsl::sd_vector<>::select_1_type(
                    &this->cols[k].v);
            rlpbwt_thr::update(new_column, pref, div);
        }
        // build another column with last prefix and divergence arrays
        auto col = rlpbwt_thr::build_column(new_column, pref, div);
        this->cols[k] = col;
        this->cols[k].rank_runs = sdsl::sd_vector<>::rank_1_type(
                &this->cols[k].runs);
        this->cols[k].select_runs = sdsl::sd_vector<>::select_1_type(
                &this->cols[k].runs);
        this->cols[k].rank_u = sdsl::sd_vector<>::rank_1_type(
                &this->cols[k].u);
        this->cols[k].select_u = sdsl::sd_vector<>::select_1_type(
                &this->cols[k].u);
        this->cols[k].rank_v = sdsl::sd_vector<>::rank_1_type(
                &this->cols[k].v);
        this->cols[k].select_v = sdsl::sd_vector<>::select_1_type(
                &this->cols[k].v);
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
            std::cout << tmp_height << ", " << tmp_width << "\n";
            std::vector<column_thr> tmp_cols(tmp_width);
            this->cols = std::vector<column_thr>(tmp_width + 1);
            input_matrix.clear();
            input_matrix.seekg(0, std::ios::beg);
            std::vector<unsigned int> pref(tmp_height);
            sdsl::int_vector<> div(tmp_height);
            for (unsigned int i = 0; i < tmp_height; i++) {
                pref[i] = i;
                div[i] = 0;
            }
            sdsl::bit_vector bv(tmp_height * tmp_width, 0);
            this->panelbv = panel_ra(tmp_height, tmp_width);
            unsigned int count = 0;
            std::string last_col;
            while (getline(input_matrix, new_column)) {
                if (verbose) {
                    std::cout << "\nnew_column " << count << "\n";
                    std::cout << new_column << "\n" << this->cols[count].runs
                              << "\n"
                              << this->cols[count].u << "\n"
                              << this->cols[count].v
                              << "\n-------------------------------\n";
                }
                new_column.erase(
                        std::remove(new_column.begin(), new_column.end(), ' '),
                        new_column.end());
                auto col = rlpbwt_thr::build_column(new_column, pref, div);
                for (unsigned int k = 0; k < new_column.size(); k++) {
                    if (new_column[k] != '0') {
                        //this->panelbv.panel[count + (k * tmp_width)] = true;
                        this->panelbv.panel[count][k] = true;
                    }
                }
                this->cols[count] = col;
                this->cols[count].rank_runs = sdsl::sd_vector<>::rank_1_type(
                        &this->cols[count].runs);
                this->cols[count].select_runs = sdsl::sd_vector<>::select_1_type(
                        &this->cols[count].runs);
                this->cols[count].rank_u = sdsl::sd_vector<>::rank_1_type(
                        &this->cols[count].u);
                this->cols[count].select_u = sdsl::sd_vector<>::select_1_type(
                        &this->cols[count].u);
                this->cols[count].rank_v = sdsl::sd_vector<>::rank_1_type(
                        &this->cols[count].v);
                this->cols[count].select_v = sdsl::sd_vector<>::select_1_type(
                        &this->cols[count].v);
                this->cols[count].rank_thr = sdsl::sd_vector<>::rank_1_type(
                        &this->cols[count].thr);
                this->cols[count].select_thr = sdsl::sd_vector<>::select_1_type(
                        &this->cols[count].thr);
                rlpbwt_thr::update(new_column, pref, div);
                count++;
                last_col = new_column;
            }

            auto col = rlpbwt_thr::build_column(last_col, pref, div);
            this->cols[count] = col;

            this->cols[count].rank_runs = sdsl::sd_vector<>::rank_1_type(
                    &this->cols[count].runs);
            this->cols[count].select_runs = sdsl::sd_vector<>::select_1_type(
                    &this->cols[count].runs);
            this->cols[count].rank_u = sdsl::sd_vector<>::rank_1_type(
                    &this->cols[count].u);
            this->cols[count].select_u = sdsl::sd_vector<>::select_1_type(
                    &this->cols[count].u);
            this->cols[count].rank_v = sdsl::sd_vector<>::rank_1_type(
                    &this->cols[count].v);
            this->cols[count].select_v = sdsl::sd_vector<>::select_1_type(
                    &this->cols[count].v);
            this->cols[count].rank_thr = sdsl::sd_vector<>::rank_1_type(
                    &this->cols[count].thr);
            this->cols[count].select_thr = sdsl::sd_vector<>::select_1_type(
                    &this->cols[count].thr);
            input_matrix.close();
        } else {
            throw FileNotFoundException{};
        }
    }
}

rlpbwt_thr::rlpbwt_thr(const char *filename, bool verbose) {
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
        this->cols = std::vector<column_thr>(tmp_width + 1);
        std::vector<unsigned int> pref(tmp_height);
        sdsl::int_vector<> div(tmp_height);
        for (unsigned int i = 0; i < tmp_height; i++) {
            pref[i] = i;
            div[i] = 0;
        }
        //sdsl::bit_vector bv(tmp_height * tmp_width, 0);
        this->panelbv = panel_ra(tmp_height, tmp_width);
        unsigned int count = 0;
        std::string last_col;
        getline(input_matrix, line);
        getline(input_matrix, line);
        while (getline(input_matrix, line) && !line.empty()) {
            std::cout << count << "\r";
            std::istringstream is_col(line);
            is_col >> garbage >> garbage >> garbage >> garbage >> new_column;
            if (verbose) {
                std::cout << "\nnew_column " << count << "\n";
                std::cout << new_column << "\n" << this->cols[count].runs
                          << "\n"
                          << this->cols[count].u << "\n"
                          << this->cols[count].v
                          << "\n-------------------------------\n";
            }
            auto col = rlpbwt_thr::build_column(new_column, pref, div);
            for (unsigned int k = 0; k < new_column.size(); k++) {
                if (new_column[k] != '0') {
                    //this->panelbv.panel[count + (k * tmp_width)] = true;
                    this->panelbv.panel[count][k] = true;

                }
            }

            this->cols[count] = col;
            this->cols[count].rank_runs = sdsl::sd_vector<>::rank_1_type(
                    &this->cols[count].runs);
            this->cols[count].select_runs = sdsl::sd_vector<>::select_1_type(
                    &this->cols[count].runs);
            this->cols[count].rank_u = sdsl::sd_vector<>::rank_1_type(
                    &this->cols[count].u);
            this->cols[count].select_u = sdsl::sd_vector<>::select_1_type(
                    &this->cols[count].u);
            this->cols[count].rank_v = sdsl::sd_vector<>::rank_1_type(
                    &this->cols[count].v);
            this->cols[count].select_v = sdsl::sd_vector<>::select_1_type(
                    &this->cols[count].v);
            this->cols[count].rank_thr = sdsl::sd_vector<>::rank_1_type(
                    &this->cols[count].thr);
            this->cols[count].select_thr = sdsl::sd_vector<>::select_1_type(
                    &this->cols[count].thr);
            rlpbwt_thr::update(new_column, pref, div);
            last_col = new_column;
            count++;

        }
        auto col = rlpbwt_thr::build_column(last_col, pref, div);
        this->cols[count] = col;
        this->cols[count].rank_runs = sdsl::sd_vector<>::rank_1_type(
                &this->cols[count].runs);
        this->cols[count].select_runs = sdsl::sd_vector<>::select_1_type(
                &this->cols[count].runs);
        this->cols[count].rank_u = sdsl::sd_vector<>::rank_1_type(
                &this->cols[count].u);
        this->cols[count].select_u = sdsl::sd_vector<>::select_1_type(
                &this->cols[count].u);
        this->cols[count].rank_v = sdsl::sd_vector<>::rank_1_type(
                &this->cols[count].v);
        this->cols[count].select_v = sdsl::sd_vector<>::select_1_type(
                &this->cols[count].v);
        this->cols[count].rank_thr = sdsl::sd_vector<>::rank_1_type(
                &this->cols[count].thr);
        this->cols[count].select_thr = sdsl::sd_vector<>::select_1_type(
                &this->cols[count].thr);
        input_matrix.close();
    } else {
        throw FileNotFoundException{};
    }
}

column_thr
rlpbwt_thr::build_column(std::string &column, std::vector<unsigned int> &pref,
                         sdsl::int_vector<> &div) {
    unsigned int height = pref.size();
    // variable for "c" value
    unsigned int count0 = 0;
    // variables for compute "u" and "v" values
    unsigned int u = 0;
    unsigned int v = 0;
    // temporary variables for compute "u" and "v" values
    unsigned int count0tmp = 0;
    unsigned int count1 = 0;
    // variable to compute curr lcs in order to eventually compute thresholds
    // unsigned int lcs = 0;
    // bool to check first symbol fo the column
    bool start = true;

    std::vector<std::pair<unsigned int, unsigned int>> rows;
    unsigned int tmp_beg = 0;
    unsigned int thr_tmp = 0;
    unsigned int lcs = 0;
    // update start and "c" value
    for (unsigned int i = 0; i < height; i++) {
        if (i == 0 && column[pref[i]] == '1') {
            start = false;
        }
        if (column[i] == '0') {
            count0++;
        }
    }

    // initialize the three bitvectors
    sdsl::bit_vector runvec(height + 1, 0);
    sdsl::bit_vector thr(height, 0);
    sdsl::bit_vector zerovec(count0, 0);
    sdsl::bit_vector onevec(height - count0, 0);

    // bools to check if we are at the beginning of a run and if we have to swap
    // the counting of zeros and ones
    bool begrun = true;
    bool pushz = false;
    bool pusho = false;

    // first swap according to first value of the column
    if (start) {
        pusho = true;
    } else {
        pushz = true;
    }

    // iteration over the entire column
    for (unsigned int i = 0; i < height; i++) {
        // if we are at the beginning of a run we save previous temporary values
        // for "u" and "v"
        if (begrun) {
            tmp_beg = pref[i];
            u = count0tmp;
            v = count1;
            begrun = false;
        }

        // increment temporary variables
        if (column[pref[i]] == '1') {
            count1++;
        } else {
            count0tmp++;
        }


        // stuff for thresholds
        if ((i == 0) || (column[pref[i]] != column[pref[i - 1]])) {
            thr_tmp = i;
            lcs = div[i];
        }
        if (div[i] < lcs) {
            thr_tmp = i;
            lcs = div[i];
        }

        if ((i == height - 1) || (column[pref[i]] != column[pref[i + 1]])) {
            // 1 in bitvectors for runs et every end of a run
            runvec[i] = true;
            // 1 in bitvectors for thresholds index
            thr[thr_tmp] = true;

            rows.emplace_back(tmp_beg, pref[i]);
            // update bitvectors for "u" and "v" and swap the case to study in
            // next run
            if (pusho) {
                if (v != 0) {
                    onevec[v - 1] = true;
                }
                std::swap(pusho, pushz);
            } else {
                if (u != 0) {
                    zerovec[u - 1] = true;
                }
                std::swap(pusho, pushz);
            }
            begrun = true;
        }
    }

    // set last bit to one in bitvectors "u" and "v"
    if (!zerovec.empty()) {
        zerovec[zerovec.size() - 1] = true;
    }
    if (!onevec.empty()) {
        onevec[onevec.size() - 1] = true;
    }

    // compress div array
    sdsl::util::bit_compress(div);
    return {start, count0, runvec, zerovec, onevec, thr, rows};
}

void rlpbwt_thr::update(std::string &column, std::vector<unsigned int> &pref,
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

unsigned int
rlpbwt_thr::lf(unsigned int col_index, unsigned int index, char symbol,
               bool verbose) const {
    // obtain "u" and "v"
    auto uv = uvtrick(col_index, index);
    if (verbose) {
        std::cout << "at col: " << col_index << " and index " << index
                  << " with symbol " << symbol << "\n";
        std::cout << uv.first << " " << uv.second << " with c: "
                  << this->cols[col_index].count_0 << "\n";
    }
    // computer lf-mapping as Durbin's w(i, s)
    if (symbol == '0') {
        return uv.first;
    } else {
        return this->cols[col_index].count_0 + uv.second;
    }
}

// TODO this function is written very bad
std::pair<unsigned int, unsigned int>
rlpbwt_thr::uvtrick(unsigned int col_index, unsigned int index) const {
    // for index 0 u = v = 0
    if (index == 0) {
        return {0, 0};
    }

    // number of previous zeros and ones
    unsigned int prev_u = 0;
    unsigned int prev_v = 0;
    // runs of previous zeros and ones in respective bitvectors
    unsigned int run_prev_u = 0;
    unsigned int run_prev_v = 0;
    // offset to finally consider inside the run
    unsigned int offset = 0;

    // compute the run of the index using rank over runs bitvector
    unsigned run = this->cols[col_index].rank_runs(index);
    // two main cases:
    // - run > 1
    // - run = 0 or run = 1 (we separate this two cases due to some errors using
    //   select)
    if (run > 1) {
        // if the run is an even index run it means that before we have the same
        // amount of zeros runs and ones runs
        if (run % 2 == 0) {
            // we compute the number of the zeros and ones runs before the
            // current run
            run_prev_u = run / 2;
            run_prev_v = run / 2;

            // we compute the number of zeros and ones before the current run
            // using select over the two bitvectors of zeros and ones
            prev_u = this->cols[col_index].select_u(run_prev_u) + 1;
            prev_v = this->cols[col_index].select_v(run_prev_v) + 1;

            // we compute the offset
            offset = index - (this->cols[col_index].select_runs(run) + 1);

            // depending on the first symbol of the column the offset increase
            // the number of zeros or the number of ones (they increase
            // alternating)
            if (this->cols[col_index].zero_first) {
                return {prev_u + offset, prev_v};
            } else {
                return {prev_u, prev_v + offset};
            }
        } else {
            // if the run is an odd index run it means that before we have x
            // runs of zeros and x+1 runs of ones, or vice-versa,
            // depending on the first symbol of the column
            if (this->cols[col_index].zero_first) {
                run_prev_u = (run / 2) + 1;
                run_prev_v = run / 2;
            } else {
                run_prev_u = run / 2;
                run_prev_v = (run / 2) + 1;
            }
            // we compute the number of zeros and ones before the current run
            // using select over the two bitvectors of zeros and ones
            prev_u = this->cols[col_index].select_u(run_prev_u) + 1;
            prev_v = this->cols[col_index].select_v(run_prev_v) + 1;

            // we compute the offset
            offset = index - (this->cols[col_index].select_runs(run) + 1);

            // depending on the first symbol of the column the offset increase
            // the number of zeros or the number of ones (they increase
            // alternating)
            if (this->cols[col_index].zero_first) {
                return {prev_u, prev_v + offset};
            } else {
                return {prev_u + offset, prev_v};
            }
        }
    } else {
        if (run == 0) {
            // if we are in the first run (of index 0) simply v is always 0 if
            // we are in a zero beginning column and vice-versa for u if we are
            // in a one beginning column. The non-zero value is the same of
            // index because it's the first run
            if (this->cols[col_index].zero_first) {
                return {index, 0};
            } else {
                return {0, index};
            }
        } else {
            // if we are in the second run (of index 1) we fix one of the two
            // values (based practically on the length of the first run),
            // according to the first symbol of the column, and we
            // increase the other one according to the index
            if (this->cols[col_index].zero_first) {
                return {this->cols[col_index].select_runs(run) + 1,
                        index - (this->cols[col_index].select_runs(run) + 1)};
            } else {
                return {index - (this->cols[col_index].select_runs(run) + 1),
                        this->cols[col_index].select_runs(run) + 1};
            }
        }
    }
}

void rlpbwt_thr::match_thr(const std::string &query, bool verbose) {
    if (query.size() != this->panelbv.w) {
        throw NotEqualLengthException{};
    }
    std::vector<unsigned int> ms_row(query.size(), 0);
    std::vector<unsigned int> ms_len(query.size(), 0);
    auto curr_prefs = this->cols[0].rows[this->cols[0].rows.size() - 1];
    auto curr_pos = curr_prefs.second;
    auto curr_index = curr_pos;
    unsigned int curr_run = this->cols[0].rank_runs(curr_index);
    char symbol = get_next_char(this->cols[0].zero_first, curr_run);
    for (unsigned int i = 0; i < query.size(); i++) {
        if (verbose) {
            std::cout << "at " << i << ": " << curr_run << " "
                      << this->cols[i].rank_thr(curr_index) << "\n";
            std::cout << curr_index << " " << curr_run << " " << curr_pos << " "
                      << symbol << "\n";
        }
        if (query[i] == symbol) {
            if (verbose) {
                std::cout << "match: ";
            }
            // save in matching statistics pos vector
            ms_row[i] = curr_pos;
            // update index, run, symbol
            curr_index = lf(i, curr_index, query[i]);
            curr_run = this->cols[i + 1].rank_runs(curr_index);
            symbol = get_next_char(this->cols[i + 1].zero_first, curr_run);
            if (verbose) {
                std::cout << "new: " << curr_index << " " << curr_run << " "
                          << curr_pos << " "
                          << symbol << "\n";
            }
        } else {
            auto thr = this->cols[i].rank_thr(curr_index);
            /* TODO if index is in the position of a thresholds
             * understand what to do
             */
            bool in_thr = false;
            if (this->cols[i].rows[curr_run].first ==
                this->cols[i].rows[curr_run].second) {
                in_thr = true;
            }
            if (this->cols[i].rows.size() == 1) {
                if (verbose) {
                    std::cout << "complete mismatch\n";
                }
                ms_row[i] = this->panelbv.h;
                curr_prefs = this->cols[i].rows[
                        this->cols[i].rows.size() - 1];
                curr_pos = curr_prefs.second;
                //curr_index = curr_pos;

                curr_prefs = this->cols[0].rows[this->cols[0].rows.size() - 1];
                curr_pos = curr_prefs.second;
                curr_index = curr_pos;
                curr_run = this->cols[i + 1].rank_runs(curr_index);
                symbol = get_next_char(this->cols[i + 1].zero_first, curr_run);
                if (verbose) {
                    std::cout << "update: " << curr_index << " " << curr_pos
                              << " "
                              << symbol << "\n";
                }
            } else if ((curr_run != 0 && !in_thr && curr_run == thr) ||
                       curr_run == this->cols[i].rows.size() - 1) {
                if (verbose) {
                    std::cout << "mismatch_up: ";
                }
                // threshold below index so we go up
                curr_index = (this->cols[i].select_runs(curr_run) + 1) - 1;
                curr_prefs = this->cols[i].rows[curr_run - 1];
                curr_pos = curr_prefs.second;
                if (verbose) {
                    std::cout << "update: " << curr_index << " " << curr_pos
                              << " "
                              << symbol << "\n";
                }
                ms_row[i] = curr_pos;
                curr_index = lf(i, curr_index, query[i]);
                curr_run = this->cols[i + 1].rank_runs(curr_index);
                symbol = get_next_char(this->cols[i + 1].zero_first,
                                       curr_run);
                if (verbose) {
                    std::cout << "new: " << curr_index << " " << curr_run
                              << " "
                              << curr_pos << " " << symbol << "\n";
                }
            } else {
                if (verbose) {
                    std::cout << "mismatch_down: ";
                }
                // threshold above index so we go down
                curr_index = (this->cols[i].select_runs(curr_run + 1) + 1);
                curr_prefs = this->cols[i].rows[curr_run + 1];
                curr_pos = curr_prefs.first;
                ms_row[i] = curr_pos;
                if (verbose) {
                    std::cout << "update: " << curr_index << " " << curr_pos
                              << " "
                              << symbol << "\n";
                }
                curr_index = lf(i, curr_index, query[i]);
                curr_run = this->cols[i + 1].rank_runs(curr_index);
                symbol = get_next_char(this->cols[i + 1].zero_first,
                                       curr_run);
                if (verbose) {
                    std::cout << "new: " << curr_index << " " << curr_run
                              << " "
                              << curr_pos << " " << symbol << "\n";
                }

            }
        }
    }
    std::cout << "ind:\t";
    for (unsigned int i = 0; i < ms_row.size(); i++) {
        std::cout << i << "\t";
    }
    std::cout << "\npos:\t";
    for (auto e: ms_row) {
        std::cout << e << "\t";
    }
    std::cout << "\nlen:\t";
    int tmp_index = 0;
    for (int i = (int) ms_row.size() - 1; i >= 0; i--) {
        if (ms_row[i] == this->panelbv.h) {
            ms_len[i] = 0;
            continue;
        }
        tmp_index = i;
        while (tmp_index > 0 &&
               query[tmp_index] == panelbv.getElem(ms_row[i], tmp_index)) {
            tmp_index--;
        }
        if (tmp_index == 0) {
            tmp_index--;
        }
        ms_len[i] = i - tmp_index;
    }
    std::vector<std::pair<unsigned int, unsigned int>> ms_match;
    for (unsigned int i = 0; i < ms_len.size(); i++) {
        std::cout << ms_len[i] << "\t";
        // TODO check if for last match
        if (i > 0 && ms_len[i] >= ms_len[i - 1] &&
            (i == ms_len.size() - 1 || ms_len[i] >= ms_len[i + 1])) {
            ms_match.emplace_back(i, ms_len[i]);
        }
    }
    std::cout << "\nmatches:\t";
    for (auto e: ms_match) {
        std::cout << "(" << e.first << ", " << e.second << ")\t";
    }
    std::cout << "\n";

}

size_t rlpbwt_thr::serialize(std::ostream &out, sdsl::structure_tree_node *v,
                             const std::string &name) {
    sdsl::structure_tree_node *child =
            sdsl::structure_tree::add_child(v, name,
                                            sdsl::util::class_name(
                                                    *this));
    size_t written_bytes = 0;
    written_bytes += this->panelbv.serialize(out, child, "panel");

    for (unsigned int i = 0; i < this->cols.size(); i++) {
        std::string label = "col_" + std::to_string(i);
        written_bytes += this->cols[i].serialize(out, child, label);
    }

    sdsl::structure_tree::add_size(child, written_bytes);
    return written_bytes;
}

void rlpbwt_thr::load(std::istream &in) {
    this->panelbv.load(in);
    auto c = new column_thr();
    for (unsigned int i = 0; i <= this->panelbv.w; i++) {
        c->load(in);
        this->cols.emplace_back(*c);
    }
}


rlpbwt_thr::rlpbwt_thr() = default;
