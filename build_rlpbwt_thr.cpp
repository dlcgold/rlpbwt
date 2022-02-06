#include <iostream>
#include <sdsl/int_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include "include/exceptions.h"
#include "include/rlpbwt_thr.h"
#include <getopt.h>

void printHelp() {
    std::cout << "Usage: RLPBTW [options]\n"
              << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -i, --input <path>" << std::endl;
    std::cout << "  -o, --output <path_folder>" << std::endl;
    std::cout << "  -l, --length <value> height of the panel" << std::endl;
    std::cout << "  -w, --width <value> width of the panel" << std::endl;
    std::cout << "  -q, --query <path>" << std::endl;
    std::cout << "  -v, --verbose: extra prints" << std::endl;
    std::cout << "  -t, --test: run tests" << std::endl;
    std::cout << "  -h, --help: show this help message and exit" << std::endl;

}

template<typename rlpbwt_t>
void build(std::string in_filename, const std::string &out_filename) {
    rlpbwt_t rlpbwt(in_filename.c_str(), false);
    clock_t START = clock();
    auto matches = rlpbwt.match_thr("010010100011101", false);
    //auto matches = rlpbwt.match_thr("111111111111111", true);
    //rlpbwt.match_tsv_tr("../input/query_tr.txt", "../output/query_tr_out.txt");
    //rlpbwt.match_tsv("../input/query.txt", "../output/query_out.txt");
    std::cout << clock() - START << " time\n";
    /*
    for (auto m: matches) {
        std::cout << "(col: " << m.first << ", len:" << m.second << ") ";
    }*/
    std::cout << "\n";
    std::ofstream file_o(out_filename);
    size_t size = rlpbwt.serialize(file_o);
    std::cout << "total size: " << size << "\n";
    file_o.close();
}

void print_size(const std::string &out_filename) {
    auto rlpbwt = new rlpbwt_thr();
    std::ifstream in;
    in.open(out_filename);
    rlpbwt->load(in);
    in.close();
    double thr = 0;
    double run = 0;
    double u = 0;
    double v = 0;
    double panel = 0;
    double sample = 0;
    for (auto &i: rlpbwt->panelbv.panel) {
        panel += sdsl::size_in_mega_bytes(i);
        std::cout << i << "\n";
    }
    unsigned long long thr_count = 0;
    unsigned long long run_count = 0;
    std::cout << "\n";
    for (auto &c: rlpbwt->cols) {
        thr += sdsl::size_in_mega_bytes(c.thr);
        run += sdsl::size_in_mega_bytes(c.runs);
        u += sdsl::size_in_mega_bytes(c.u);
        v += sdsl::size_in_mega_bytes(c.v);
        sample += sdsl::size_in_mega_bytes(c.sample_beg);
        sample += sdsl::size_in_mega_bytes(c.sample_end);
        std::cout << c;
        thr_count += c.rank_thr(c.runs.size()-1);
        run_count += c.sample_beg.size();
        std::cout << "runs: " << run_count <<"\nthrs: " << thr_count << "\n-------\n";
    }

    std::cout << "panel size: " << panel << " mb\n";
    std::cout << "run size: " << run << " mb\n";
    std::cout << "thr size: " << thr << " mb\n";
    std::cout << "u size: " << u << " mb\n";
    std::cout << "v size: " << v << " mb\n";
    std::cout << "sample size: " << sample << " mb\n";
    std::cout << "total size (excluded single value in structure): "
              << panel + run + thr + u + v + sample
              << " mb\n";
    std::cout << "runs: " << run_count <<"\nthrs: " << thr_count << "\n";

}

int main(int argc, char **argv) {
    /*
    bool verbose = false;
    std::string matrix_input;
    std::string out_dir;
    std::string query_input;
    unsigned int height;
    unsigned int width;
    int c;
    while (true) {
        static struct option long_options[] =
                {
                        {"input",   required_argument, nullptr, 'i'},
                        {"output",   required_argument, nullptr, 'o'},
                        {"length",  required_argument, nullptr, 'l'},
                        {"width",   required_argument, nullptr, 'w'},
                        {"query",   required_argument, nullptr, 'q'},
                        {"verbose", no_argument,       nullptr, 'v'},
                        {"test",    no_argument,       nullptr, 't'},
                        {"help",    no_argument,       nullptr, 'h'},
                        {nullptr, 0,                   nullptr, 0}
                };

        int option_index = 0;
        c = getopt_long(argc, argv, "i:l:w:q:vth", long_options,
                        &option_index);

        if (c == -1) {
            break;
        }

        switch (c) {
            case 'i':
                matrix_input = optarg;
                break;
            case 'o':
                out_dir = optarg;
                break;
            case 'l':
                height = std::stoi(optarg);
                break;
            case 'w':
                width = std::stoi(optarg);
                break;
            case 'q':
                query_input = optarg;
                break;
            case 'v':
                verbose = true;
                break;
            case 't':
                std::cout << "current not available\n";
                //::testing::InitGoogleTest(&argc, argv);
                //::testing::GTEST_FLAG(filter) = "BuildRlpbwtTest*";
                //::testing::GTEST_FLAG(filter) = "BuildRlpbwtSerOrig*";
                //::testing::GTEST_FLAG(filter) = "BuildBiRlpbwtTest*";
                //::testing::GTEST_FLAG(filter) = "BuildRlpbwtVCF*";
                //::testing::GTEST_FLAG(filter) = "BuildBiRlpbwtVCF*";
                //::testing::GTEST_FLAG(filter) = "BuildRlpbwtSerBV*";
                //::testing::GTEST_FLAG(filter) = "BuildRlpbwtBVtest*";
                //::testing::GTEST_FLAG(filter) = "BuildRlpbwtBVVCF*";
                //::testing::GTEST_FLAG(filter) = "BuildRlpbwtThr*";
                //::testing::GTEST_FLAG(filter) = "BuildRlpbwtSerThr*";
                //::testing::GTEST_FLAG(filter) = "BuildRlpbwtNewThr*";
                //return RUN_ALL_TESTS();
                break;
            case 'h':
                printHelp();
                exit(EXIT_SUCCESS);
            default:
                printHelp();
                exit(EXIT_FAILURE);
        }
    }*/

    std::string in_filename("../input/sample_new.txt");
    std::string out_filename("../output/samplenew.txt.pbwt");

    build<rlpbwt_thr>(in_filename, out_filename);
    print_size(out_filename);
    return 0;
}