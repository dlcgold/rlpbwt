#include <iostream>
#include <getopt.h>
#include "include/rlpbwt_ms.h"
#include "include/rlpbwt_bv.h"
#include "include/rlpbwt_naive.h"

void printHelp() {
    std::cout << "Usage: RLPBWT [options]\n"
              << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -i, --input_file <path>\t macs file for panel" << std::endl;
    std::cout << "  -m, --memorize <path>\t  path to save serialization " << std::endl;
    std::cout << "  -l, --load <path>\t path to load serialization" << std::endl;
    std::cout << "  -s, --input_slp <path>\t path to SLP file" << std::endl;
    std::cout << "  -o, --output <path>\t path to query output" << std::endl;
    std::cout << "  -q, --query <path>\t path to macs query file" << std::endl;
    std::cout << "  -N, --Naive\t naive RLPBWT (only one mode allowed)"
              << std::endl;
    std::cout << "  -B, --Bv\t bitvectors RLPBWT (only one mode allowed)"
              << std::endl;
    std::cout << "  -S, --Slp\t slp RLPBWT "
              << "(only one mode allowed | slp file required)"
              << std::endl;
    std::cout << "  -P, --Panel\t panel RLPBWT (only one mode allowed)"
              << std::endl;
    std::cout << "  -t, --thresholds\t enable thresholds (slp/panel mode only)"
              << std::endl;
    std::cout << "  -e, --extend\t extend matches (slp/panel mode only)"
              << std::endl;
    std::cout << "  -v, --verbose\t extra prints" << std::endl;
    std::cout << "  -V, --fverbose\t extra prints for functions (cautions)" << std::endl;
    std::cout << "  -h, --help\t show this help message and exit" << std::endl;
}

/*
template<typename rlpbwt_t>
void build(std::string in_filename, const std::string &out_filename) {
    rlpbwt_t rlpbwt(in_filename.c_str(), true);
    clock_t START = clock();
    auto matches = rlpbwt.match_thr("010010100011101", false);
    //auto matches = rlpbwt.match_thr("111111111111111", true);
    //rlpbwt.match_tsv_tr("../input/query_tr.txt", "../output/query_tr_out.txt");
    //rlpbwt.match_tsv("../input/query.txt", "../output/query_out.txt");
    //std::cout << clock() - START << " time\n";

    for (auto m: matches) {
        std::cout << "(col: " << m.first << ", len:" << m.second << ") ";
    }
    std::cout << "\n";
    std::ofstream file_o(out_filename);
    size_t size = rlpbwt.serialize(file_o);
    std::cout << "total size: " << size << "\n";
    file_o.close();
}

template<typename rlpbwt_t>
void print_size(const std::string &out_filename) {
    auto rlpbwt = new rlpbwt_t();
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
    double sparse = 0;
    for (auto &i: rlpbwt->panel->panel) {
        panel += sdsl::size_in_mega_bytes(i);
        auto s = sdsl::sd_vector<>(i);
        sparse += sdsl::size_in_mega_bytes(s);
    }
    unsigned long long thr_count = 0;
    unsigned long long run_count = 0;
    for (auto &c: rlpbwt->cols) {
        thr += sdsl::size_in_mega_bytes(c.thr);
        run += sdsl::size_in_mega_bytes(c.runs);
        u += sdsl::size_in_mega_bytes(c.u);
        v += sdsl::size_in_mega_bytes(c.v);
        sample += sdsl::size_in_mega_bytes(c.sample_beg);
        sample += sdsl::size_in_mega_bytes(c.sample_end);
        thr_count += c.rank_thr(c.runs.size() - 1);
        run_count += c.sample_beg.size();
    }

    std::cout << "panel size: " << panel << " mb\n";
    std::cout << "panel size sparse: " << sparse << " mb\n";
    std::cout << "run size: " << run << " mb\n";
    std::cout << "thr size: " << thr << " mb\n";
    std::cout << "u size: " << u << " mb\n";
    std::cout << "v size: " << v << " mb\n";
    std::cout << "sample size: " << sample << " mb\n";
    std::cout << "total size (excluded single value in structure): "
              << panel + run + thr + u + v + sample
              << " mb\n";
    std::cout << "runs: " << run_count << "\nthrs: " << thr_count << "\n";

}
*/

int main(int argc, char **argv) {
    if (argc == 1) {
        printHelp();
        exit(EXIT_SUCCESS);
    }
    bool verbose = false;
    bool print_verbose = false;
    bool naive = false;
    bool bv = false;
    bool slp = false;
    bool panel = false;
    bool thr = false;
    bool extend = false;
    bool query = false;
    std::string matrix_input = "";
    std::string memorize_file = "";
    std::string load_file = "";
    std::string slp_input = "";
    std::string output = "";
    std::string query_input = "";
    int c;

    while (true) {
        static struct option long_options[] =
                {
                        {"input",      required_argument, nullptr, 'i'},
                        {"memorize",   required_argument, nullptr, 'm'},
                        {"load",       required_argument, nullptr, 'l'},
                        {"slp",        required_argument, nullptr, 's'},
                        {"output",     required_argument, nullptr, 'o'},
                        {"query",      required_argument, nullptr, 'q'},
                        {"Naive",      no_argument,       nullptr, 'N'},
                        {"Bv",         no_argument,       nullptr, 'B'},
                        {"Slp",        no_argument,       nullptr, 'S'},
                        {"Panel",      no_argument,       nullptr, 'P'},
                        {"thresholds", no_argument,       nullptr, 't'},
                        {"extend",     no_argument,       nullptr, 'e'},
                        {"fverbose",   no_argument,       nullptr, 'V'},
                        {"verbose",    no_argument,       nullptr, 'v'},
                        {"help",       no_argument,       nullptr, 'h'},
                        {nullptr, 0,                      nullptr, 0}
                };

        int option_index = 0;
        c = getopt_long(argc, argv, "i:m:l:s:o:q:NBSPtevVh", long_options,
                        &option_index);

        if (c == -1) {
            break;
        }

        switch (c) {
            case 'i':
                matrix_input = optarg;
                break;
            case 'm':
                memorize_file = optarg;
                break;
            case 'l':
                load_file = optarg;
                break;
            case 's':
                slp_input = optarg;
                break;
            case 'o':
                output = optarg;
                break;
            case 'q':
                query_input = optarg;
                break;
            case 'N':
                naive = true;
                break;
            case 'B':
                bv = true;
                break;
            case 'S':
                slp = true;
                break;
            case 'P':
                panel = true;
                break;
            case 't':
                thr = true;
                break;
            case 'e':
                extend = true;
                break;
            case 'v':
                print_verbose = true;
                break;
            case 'V':
                verbose = true;
                break;
            case 'h':
                printHelp();
                exit(EXIT_SUCCESS);
            default:
                printHelp();
                exit(EXIT_FAILURE);
        }
    }

    if (!query_input.empty()) {
        query = true;
    }
    if (!(naive || bv || slp || panel)) {
        std::cerr << "Error: one mode required\n";
        printHelp();
        exit(EXIT_FAILURE);
    }
    if (memorize_file.empty() && (query_input.empty() || output.empty())) {
        std::cerr << "Error: nothing to do\n";
        printHelp();
        exit(EXIT_FAILURE);
    }
    if (query && output.empty()) {
        std::cerr << "Error: output file required if query requested\n";
        printHelp();
        exit(EXIT_FAILURE);
    }

    if ((naive && (bv || slp || panel)) || (bv && (naive || slp || panel)) ||
        (slp && (naive || bv || panel)) || (panel && (naive || bv || slp))) {
        std::cerr << "Error: only one mode allowed\n";
        exit(EXIT_FAILURE);
    }


    if (naive) {
        if (matrix_input.empty() && load_file.empty()) {
            std::cerr << "Error: input or load file required\n";
            exit(EXIT_FAILURE);
        }
        if (query) {
            if (output.empty()) {
                std::cerr << "Error: output file required\n";
                exit(EXIT_FAILURE);
            }
            if (query_input.empty()) {
                std::cerr << "Error: query file required\n";
                exit(EXIT_FAILURE);
            }
        }

        if (extend) {
            std::cout << "matches will not be extended\n";
        }
        if (thr) {
            std::cout << "thresholds will not be used\n";
        }
        clock_t START = clock();
        if (load_file.empty()) {
            rlpbwt_naive rlpbwt(matrix_input.c_str(), verbose);
            if (!memorize_file.empty()) {
                std::ofstream outstream;
                outstream.open(memorize_file.c_str());
                rlpbwt.serialize(outstream);
                outstream.close();
            }
            auto time_build = (float) (clock() - START) / CLOCKS_PER_SEC;
            if (print_verbose) {
                std::cout << "rlpbwt: " << rlpbwt.size_in_bytes(verbose)
                          << " bytes\n";
                std::cout << "rlpbwt: " << rlpbwt.size_in_mega_bytes(verbose)
                          << " megabytes\n----\n";
                std::cout << "estimated dense size: "
                          << dense_size_megabyte(rlpbwt.height, rlpbwt.width)
                          << " megabytes\n----\n";
            }
            std::cout << "built/loaded in " << time_build << " s\n";

            if (query) {
                START = clock();
                rlpbwt.match_tsv_tr(query_input.c_str(), output.c_str());
                auto time_query = (float) (clock() - START) / CLOCKS_PER_SEC;
                std::cout << "queried in " << time_query << " s\n";
            }
        } else {
            rlpbwt_naive rlpbwt;
            std::ifstream load;
            load.open(load_file.c_str());
            rlpbwt.load(load);
            load.close();

            auto time_build = (float) (clock() - START) / CLOCKS_PER_SEC;
            if (print_verbose) {
                std::cout << "rlpbwt: " << rlpbwt.size_in_bytes(verbose)
                          << " bytes\n";
                std::cout << "rlpbwt: " << rlpbwt.size_in_mega_bytes(verbose)
                          << " megabytes\n----\n";
                std::cout << "estimated dense size: "
                          << dense_size_megabyte(rlpbwt.height, rlpbwt.width)
                          << " megabytes\n----\n";
            }
            std::cout << "built/loaded in " << time_build << " s\n";

            if (query) {
                START = clock();
                rlpbwt.match_tsv_tr(query_input.c_str(), output.c_str());
                auto time_query = (float) (clock() - START) / CLOCKS_PER_SEC;
                std::cout << "queried in " << time_query << " s\n";
            }
        }
    }
    if (bv) {
        if (matrix_input.empty() && load_file.empty()) {
            std::cerr << "Error: input or load file required\n";
            exit(EXIT_FAILURE);
        }
        if (query) {
            if (output.empty()) {
                std::cerr << "Error: output file required\n";
                exit(EXIT_FAILURE);
            }
            if (query_input.empty()) {
                std::cerr << "Error: query file required\n";
                exit(EXIT_FAILURE);
            }
        }
        if (extend) {
            std::cout << "matches will not be extended\n";
        }
        if (thr) {
            std::cout << "thresholds will not be used\n";
        }
        clock_t START = clock();

        //rlpbwt_bv rlpbwt(matrix_input.c_str(), verbose);
        if (load_file.empty()) {
            rlpbwt_bv rlpbwt(matrix_input.c_str(), verbose);
            if (!memorize_file.empty()) {
                std::ofstream outstream;
                outstream.open(memorize_file.c_str());
                rlpbwt.serialize(outstream);
                outstream.close();
            }
            auto time_build = (float) (clock() - START) / CLOCKS_PER_SEC;
            if (print_verbose) {
                std::cout << "rlpbwt: " << rlpbwt.size_in_bytes(verbose)
                          << " bytes\n";
                std::cout << "rlpbwt: " << rlpbwt.size_in_mega_bytes(verbose)
                          << " megabytes\n----\n";
                std::cout << "estimated dense size: "
                          << dense_size_megabyte(rlpbwt.height, rlpbwt.width)
                          << " megabytes\n----\n";
            }
            std::cout << "built/loaded in " << time_build << " s\n";
            if (query) {
                START = clock();
                rlpbwt.match_tsv_tr(query_input.c_str(), output.c_str());
                auto time_query = (float) (clock() - START) / CLOCKS_PER_SEC;
                std::cout << "queried in " << time_query << " s\n";
            }
        } else {
            rlpbwt_bv rlpbwt;
            std::ifstream load;
            load.open(load_file.c_str());
            rlpbwt.load(load);
            load.close();
            auto time_build = (float) (clock() - START) / CLOCKS_PER_SEC;
            if (print_verbose) {
                std::cout << "rlpbwt: " << rlpbwt.size_in_bytes(verbose)
                          << " bytes\n";
                std::cout << "rlpbwt: " << rlpbwt.size_in_mega_bytes(verbose)
                          << " megabytes\n----\n";
                std::cout << "estimated dense size: "
                          << dense_size_megabyte(rlpbwt.height, rlpbwt.width)
                          << " megabytes\n----\n";
            }
            std::cout << "built/loaded in " << time_build << " s\n";
            if (query) {
                START = clock();
                rlpbwt.match_tsv_tr(query_input.c_str(), output.c_str());
                auto time_query = (float) (clock() - START) / CLOCKS_PER_SEC;
                std::cout << "queried in " << time_query << " s\n";
            }
        }
    }

    if (slp) {
        if (matrix_input.empty() && load_file.empty()) {
            std::cerr << "Error: input or load file required\n";
            exit(EXIT_FAILURE);
        }
        if (query) {
            if (output.empty()) {
                std::cerr << "Error: output file required\n";
                exit(EXIT_FAILURE);
            }
            if (query_input.empty()) {
                std::cerr << "Error: query file required\n";
                exit(EXIT_FAILURE);
            }
        }
        if (slp_input.empty()) {
            std::cerr << "Error: slp required\n";
            exit(EXIT_FAILURE);
        }

        if (thr) {
            std::cout << "thresholds enabled\n";
        } else {
            std::cout << "thresholds are not enabled\n";
        }
        clock_t START = clock();
        if (load_file.empty()) {
            rlpbwt_ms<slp_panel_ra> rlpbwt(matrix_input.c_str(), thr, verbose,
                                           slp_input.c_str());
            if (extend) {
                rlpbwt.extend();
            }
            if (!memorize_file.empty()) {
                std::ofstream outstream;
                outstream.open(memorize_file.c_str());
                rlpbwt.serialize(outstream);
                outstream.close();
            }

            auto time_build = (float) (clock() - START) / CLOCKS_PER_SEC;
            if (print_verbose) {
                std::cout << "rlpbwt: " << rlpbwt.size_in_bytes(verbose)
                          << " bytes\n";
                std::cout << "rlpbwt: " << rlpbwt.size_in_mega_bytes(verbose)
                          << " megabytes\n----\n";
                std::cout << "estimated dense size: "
                          << dense_size_megabyte(rlpbwt.height, rlpbwt.width)
                          << " megabytes\n----\n";
            }
            std::cout << "built/loaded in " << time_build << " s\n";
            if (query) {
                START = clock();
                if (thr) {
                    rlpbwt.match_tsv_tr_thr(query_input.c_str(), output.c_str(),
                                            extend);
                } else {
                    rlpbwt.match_tsv_tr_lce(query_input.c_str(), output.c_str(),
                                            extend);
                }
                auto time_query = (float) (clock() - START) / CLOCKS_PER_SEC;
                std::cout << "queried in " << time_query << " s\n";
            }
        } else {
            rlpbwt_ms<slp_panel_ra> rlpbwt;
            std::ifstream load;
            load.open(load_file.c_str());
            rlpbwt.load(load, slp_input.c_str());
            load.close();
            if (extend) {
                rlpbwt.extend();
            }
            auto time_build = (float) (clock() - START) / CLOCKS_PER_SEC;
            if (print_verbose) {
                std::cout << "rlpbwt: " << rlpbwt.size_in_bytes(verbose)
                          << " bytes\n";
                std::cout << "rlpbwt: " << rlpbwt.size_in_mega_bytes(verbose)
                          << " megabytes\n----\n";
                std::cout << "estimated dense size: "
                          << dense_size_megabyte(rlpbwt.height, rlpbwt.width)
                          << " megabytes\n----\n";
            }
            std::cout << "built/loaded in " << time_build << " s\n";
            if (query) {
                std::cout << "here\n";
                std::cout << query_input << "\n" << output << "\n";
                START = clock();
                if (thr) {
                    rlpbwt.match_tsv_tr_thr(query_input.c_str(), output.c_str(),
                                            extend);
                } else {
                    std::cout << "l ce\n";
                    rlpbwt.match_tsv_tr_lce(query_input.c_str(), output.c_str(),
                                            extend);
                }
                auto time_query = (float) (clock() - START) / CLOCKS_PER_SEC;
                std::cout << "queried in " << time_query << " s\n";
            }
        }
    }

    if (panel) {
        if (matrix_input.empty() && load_file.empty()) {
            std::cerr << "Error: input or load file required\n";
            exit(EXIT_FAILURE);
        }
        if (query) {
            if (output.empty()) {
                std::cerr << "Error: output file required\n";
                exit(EXIT_FAILURE);
            }
            if (query_input.empty()) {
                std::cerr << "Error: query file required\n";
                exit(EXIT_FAILURE);
            }
        }
        if (!thr) {
            std::cerr << "thresholds required\n";
            exit(EXIT_FAILURE);
        }
        if (!slp_input.empty()) {
            std::cout << "slp will not be used\n";
        }
        clock_t START = clock();
        if (load_file.empty()) {
            rlpbwt_ms<panel_ra> rlpbwt(matrix_input.c_str(), thr, verbose);
            if (extend) {
                rlpbwt.extend();
            }
            if (!memorize_file.empty()) {
                std::ofstream outstream;
                outstream.open(memorize_file.c_str());
                rlpbwt.serialize(outstream);
                outstream.close();
            }
            auto time_build = (float) (clock() - START) / CLOCKS_PER_SEC;
            if (print_verbose) {
                std::cout << "rlpbwt: " << rlpbwt.size_in_bytes(verbose)
                          << " bytes\n";
                std::cout << "rlpbwt: " << rlpbwt.size_in_mega_bytes(verbose)
                          << " megabytes\n----\n";
                std::cout << "estimated dense size: "
                          << dense_size_megabyte(rlpbwt.height, rlpbwt.width)
                          << " megabytes\n----\n";
            }
            std::cout << "built/loaded in " << time_build << " s\n";
            if (query) {
                START = clock();
                rlpbwt.match_tsv_tr_thr(query_input.c_str(), output.c_str(),
                                        extend);
                auto time_query = (float) (clock() - START) / CLOCKS_PER_SEC;
                std::cout << "queried in " << time_query << " s\n";
            }
        } else {
            rlpbwt_ms<panel_ra> rlpbwt;
            std::ifstream load;
            load.open(load_file.c_str());
            rlpbwt.load(load);
            load.close();
            auto time_build = (float) (clock() - START) / CLOCKS_PER_SEC;
            if (extend) {
                rlpbwt.extend();
            }
            if (print_verbose) {
                std::cout << "rlpbwt: " << rlpbwt.size_in_bytes(verbose)
                          << " bytes\n";
                std::cout << "rlpbwt: " << rlpbwt.size_in_mega_bytes(verbose)
                          << " megabytes\n----\n";
                std::cout << "estimated dense size: "
                          << dense_size_megabyte(rlpbwt.height, rlpbwt.width)
                          << " megabytes\n----\n";
            }
            std::cout << "built/loaded in " << time_build << " s\n";
            if (query) {
                START = clock();
                rlpbwt.match_tsv_tr_thr(query_input.c_str(), output.c_str(),
                                        extend);
                auto time_query = (float) (clock() - START) / CLOCKS_PER_SEC;
                std::cout << "queried in " << time_query << " s\n";
            }
        }
    }

    /*std::string in_filename("../input/sample_new.txt");
    std::string out_filename("../output/samplenew.txt.pbwt");

    build<rlpbwt_ms<panel_ra>>(in_filename, out_filename);
    print_size<rlpbwt_ms<panel_ra>>(out_filename);
    */
    return 0;
}
