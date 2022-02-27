#include <iostream>
#include <gtest/gtest.h>
#include <gperftools/heap-profiler.h>
#include "benchmark_test.h"
#include "../include/utils.h"
#include "../include/exceptions.h"
#include "../include/rlpbwt_bv.h"
#include "../include/panel_ra.h"
#include "../include/rlpbwt_ms.h"
#include "../include/rlpbwt_naive.h"

/*
TEST (BuildRlpbwtSerBV, TestSer) {
    auto rlpbwt_t = new rlpbwt_bv("../input/sample2.txt", false, false);

    auto filename = "../output/sample2bv_panel.ser";
    std::ofstream file_o;
    std::ifstream file_i;
    unsigned int size = 0;
    file_o.open(filename);
    size = rlpbwt_t->cols[0].serialize(file_o);
    std::cout << size << "\n";
    file_o.close();
    file_i.open(filename);
    auto col = new column_bv();
    col->load(file_i);
    std::cout << col->runs << "\n";
    std::cout << col->rank_runs(18) << " " << col->select_runs(2) << "\n";
    file_i.close();

    filename = "../output/sample2bv_pbwt.ser";
    file_o.open(filename);
    size = rlpbwt_t->serialize(file_o);
    std::cout << size << "\n";
    file_o.close();
    file_i.open(filename);
    auto rlpbwt = new rlpbwt_bv();
    rlpbwt->load(file_i);
    std::cout << rlpbwt->cols[13] << "\n";
    std::cout << rlpbwt->cols.size() << " vs " << rlpbwt_t->cols.size() << "\n";
    file_i.close();
}


TEST (BuildRlpbwtBVtest, TestBuildQuery) {
    rlpbwt_bv rlpbwt_bv("../input/sample2.txt", false);

    auto matches = rlpbwt_bv.external_match("010010100011101", 1);
    for (const auto &m: matches) {
        std::cout << m << "\n";
    }
    auto match0 = match(0, 5, 4);
    auto match1 = match(2, 8, 1);
    auto match2 = match(3, 9, 1);
    auto match3 = match(7, 11, 1);
    auto match4 = match(11, 14, 3);
    EXPECT_EQ(matches[0], match0);
    EXPECT_EQ(matches[1], match1);
    EXPECT_EQ(matches[2], match2);
    EXPECT_EQ(matches[3], match3);
    EXPECT_EQ(matches[4], match4);
}

TEST (BuildRlpbwtBVVCF, TestBuildQuery) {
    rlpbwt_bv rlpbwt_bv("../input/sample_panel.vcf", true);
    EXPECT_EQ(rlpbwt_bv.height, 900);
    EXPECT_EQ(rlpbwt_bv.width, 500);
    clock_t START = clock();
    // TODO add test check for basic_matches
    rlpbwt_bv.external_match_vcf("../input/sample_query.vcf", 255, false);
    std::cout << clock() - START << " time\n";
}*/

TEST(RlpbwtRaTest, TestBuildQuery) {
    rlpbwt_ms<panel_ra> rlpbwtRa("../input/sample_new.txt", true);
    auto matches = rlpbwtRa.match_thr("010010100011101", true);
    std::cout << matches << "\n";
    rlpbwt_ms<slp_panel_ra> rlpbwtSlp("../input/sample_new.txt", true,
                                      false, "../input/sample.slp");

    matches = rlpbwtRa.match_thr("010010100011101", false);
    std::cout << matches << "\n";
    rlpbwtSlp.extend();

    matches = rlpbwtSlp.match_thr("010010100011101", false);
    std::cout << matches << "\n";
    std::cout << rlpbwtSlp.is_extended << "\n";
    auto filename = "../output/phi_pbwt.ser";
    std::ofstream file_o;
    file_o.open(filename);
    auto size = rlpbwtSlp.serialize(file_o);
    std::cout << size << "\n";
    file_o.close();
    std::ifstream file_i;
    file_i.open(filename);
    auto rlpbwt = new rlpbwt_ms<slp_panel_ra>();
    rlpbwt->load(file_i, "../input/sample.slp");
    file_i.close();
    for (unsigned int i = 0; i < rlpbwt->panel->h; i++) {
        std::cout << i << ": ";
        for (unsigned int j = 0; j < rlpbwt->panel->w; j++) {
            std::cout << rlpbwt->phi->phi(i, j).value_or(
                    rlpbwt->panel->h) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "----------\n";
    for (unsigned int i = 0; i < rlpbwt->panel->h; i++) {
        std::cout << i << ": ";
        for (unsigned int j = 0; j < rlpbwt->panel->w; j++) {
            std::cout << rlpbwt->phi->phi_inv(i, j).value_or(
                    rlpbwt->panel->h) << " ";
        }
        std::cout << "\n";
    }

}

TEST(Lce, Test) {
    rlpbwt_ms<slp_panel_ra> rlpbwtSlp("../input/sample_new.txt", true, false,
                                      "../input/sample.slp");
    rlpbwtSlp.extend();
    for (unsigned int i = 0; i < rlpbwtSlp.panel->h; i++) {
        std::cout << i << ": ";
        for (unsigned int j = 0; j < rlpbwtSlp.panel->w; j++) {
            std::cout << rlpbwtSlp.phi->phi(i, j).value_or(
                    rlpbwtSlp.panel->h) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "----------\n";
    for (unsigned int i = 0; i < rlpbwtSlp.panel->h; i++) {
        std::cout << i << ": ";
        for (unsigned int j = 0; j < rlpbwtSlp.panel->w; j++) {
            std::cout << rlpbwtSlp.phi->phi_inv(i, j).value_or(
                    rlpbwtSlp.panel->h) << " ";
        }
        std::cout << "\n";
    }
    auto matches = rlpbwtSlp.match_lce("010010100011101", true, false);
    std::cout << matches;
    matches = rlpbwtSlp.match_lce("000000000000000", true, false);
    std::cout << matches;
}

TEST(RlpbwtNaive, BuildQuery) {
    rlpbwt_naive rlpbwt("../input/sample_new.txt");
    std::cout << rlpbwt;
    auto matches = rlpbwt.external_match("010010100011101");
    std::cout << matches;

    rlpbwt_bv rlpbwtbv("../input/sample_new.txt");
    std::cout << rlpbwt;
    matches = rlpbwt.external_match("010010100011101");
    std::cout << matches;
}

TEST(MixRlpbwt, TestBuildSize) {
    bool size_verbose = false;
    bool match_verbose = false;
    auto input = "../input/sample_new.txt";
    auto slp = "../input/sample.slp";

    rlpbwt_naive rlpbwt(input);
    std::cout << "rlpbwt:\n" << rlpbwt.size_in_mega_bytes(size_verbose)
              << " megabytes\n----\n";

    rlpbwt_bv rlpbwtbv(input);
    std::cout << "rlpbwt_bv:\n" << rlpbwtbv.size_in_mega_bytes(size_verbose)
              << " megabytes\n----\n";

    rlpbwt_ms<slp_panel_ra> rlpbwtSlp(input, true,
                                      false, slp);
    std::cout << "rlpbwt_slp:\n" << rlpbwtSlp.size_in_mega_bytes(size_verbose)
              << " megabytes\n----\n";
    rlpbwtSlp.extend();
    std::cout << "rlpbwt_slp_ext:\n"
              << rlpbwtSlp.size_in_mega_bytes(size_verbose)
              << " megabytes\n----\n";

    rlpbwt_ms<panel_ra> rlpbwtPan(input, true, false);
    std::cout << "rlpbwt_pa:\n" << rlpbwtPan.size_in_mega_bytes(size_verbose)
              << " megabytes\n----\n";
    rlpbwtPan.extend();
    std::cout << "rlpbwt_pa_ext:\n"
              << rlpbwtPan.size_in_mega_bytes(size_verbose)
              << " megabytes\n----\n";
    clock_t START = clock();
    auto matches_naive = rlpbwt.external_match("010010100011101");
    if (match_verbose) {
        std::cout << matches_naive;
    }
    std::cout << "rlpbwt_naive time: "
              << (float) (clock() - START) / CLOCKS_PER_SEC << " s\n";
    START = clock();
    matches_naive = rlpbwtbv.external_match("010010100011101");
    if (match_verbose) {
        std::cout << matches_naive;
    }
    std::cout << "rlpbwt_bv time: "
              << (float) (clock() - START) / CLOCKS_PER_SEC << " s\n";
    START = clock();
    auto matches = rlpbwtSlp.match_lce("010010100011101", true);
    if (match_verbose) {
        std::cout << matches;
    }
    std::cout << "slp_lce_ext time: "
              << (float) (clock() - START) / CLOCKS_PER_SEC << " s\n";
    START = clock();
    matches = rlpbwtSlp.match_lce("010010100011101", false);
    if (match_verbose) {
        std::cout << matches;
    }
    std::cout << "slp_lce_noext time: "
              << (float) (clock() - START) / CLOCKS_PER_SEC << " s\n";
    START = clock();
    matches = rlpbwtSlp.match_thr("010010100011101", true);
    if (match_verbose) {
        std::cout << matches;
    }
    std::cout << "slp_thr_ext time: "
              << (float) (clock() - START) / CLOCKS_PER_SEC << " s\n";
    START = clock();
    matches = rlpbwtSlp.match_thr("010010100011101", false);
    if (match_verbose) {
        std::cout << matches;
    }
    std::cout << "slp_thr_noext time: "
              << (float) (clock() - START) / CLOCKS_PER_SEC << " s\n";
    START = clock();
    matches = rlpbwtPan.match_thr("010010100011101", true);
    if (match_verbose) {
        std::cout << matches;
    }
    std::cout << "panel_thr_ext time: "
              << (float) (clock() - START) / CLOCKS_PER_SEC << " s\n";
    START = clock();
    matches = rlpbwtPan.match_thr("010010100011101", false);
    if (match_verbose) {
        std::cout << matches;
    }
    std::cout << "panel_thr_noext time: "
              << (float) (clock() - START) / CLOCKS_PER_SEC << " s\n";
}

TEST(Benchmark, Query) {
    bool size_verbose = false;
    bool match_verbose = false;
    auto input = "../input/sample_new.txt";
    auto slp = "../input/sample.slp";
    auto queries = "../input/query.txt";
    auto queries_tr = "../input/query_tr.txt";

    rlpbwt_naive rlpbwt(input);
    std::cout << "rlpbwt:\n" << rlpbwt.size_in_mega_bytes(size_verbose)
              << " megabytes\n----\n";

    rlpbwt_bv rlpbwtbv(input);
    std::cout << "rlpbwt_bv:\n" << rlpbwtbv.size_in_mega_bytes(size_verbose)
              << " megabytes\n----\n";

    rlpbwt_ms<slp_panel_ra> rlpbwtSlp(input, true,
                                      false, slp);
    std::cout << "rlpbwt_slp:\n" << rlpbwtSlp.size_in_mega_bytes(size_verbose)
              << " megabytes\n----\n";

    rlpbwt_ms<panel_ra> rlpbwtPan(input, true, false);
    std::cout << "rlpbwt_pa:\n" << rlpbwtPan.size_in_mega_bytes(size_verbose)
              << " megabytes\n----\n";

    std::cout << "estimated dense size: "
              << dense_size_megabyte(rlpbwt.height, rlpbwt.width)
              << " megabytes\n----\n";

    rlpbwt.match_tsv(queries, "../output/rlpbwt_queries.txt", false);
    rlpbwt.match_tsv_tr(queries_tr, "../output/rlpbwt_queries_tr.txt");


    rlpbwtbv.match_tsv(queries, "../output/rlpbwtbv_queries.txt", false);
    rlpbwtbv.match_tsv_tr(queries_tr, "../output/rlpbwtbv_queries_tr.txt");


    rlpbwtSlp.match_tsv_thr(queries,
                            "../output/rlpbwtSlp_thr_noext_queries.txt",
                            false);
    rlpbwtSlp.match_tsv_lce(queries,
                            "../output/rlpbwtSlp_lce_noext_queries.txt",
                            false);
    rlpbwtSlp.match_tsv_thr(queries, "../output/rlpbwtSlp_thr_ext_queries.txt",
                            true);
    rlpbwtSlp.match_tsv_lce(queries, "../output/rlpbwtSlp_lce_ext_queries.txt",
                            true);
    rlpbwtSlp.match_tsv_tr_thr(queries_tr,
                               "../output/rlpbwtSlp_thr_noext_queries_tr.txt",
                               false);
    rlpbwtSlp.match_tsv_tr_lce(queries_tr,
                               "../output/rlpbwtSlp_lce_noext_queries_tr.txt",
                               false);
    rlpbwtSlp.match_tsv_tr_thr(queries_tr,
                               "../output/rlpbwtSlp_thr_ext_queries_tr.txt",
                               true);
    rlpbwtSlp.match_tsv_tr_lce(queries_tr,
                               "../output/rlpbwtSlp_lce_ext_queries_tr.txt",
                               true);


    rlpbwtPan.match_tsv_thr(queries,
                            "../output/rlpbwtPan_thr_noext_queries.txt",
                            false);
    rlpbwtPan.match_tsv_thr(queries, "../output/rlpbwtPan_thr_ext_queries.txt",
                            true);
    rlpbwtPan.match_tsv_tr_thr(queries_tr,
                               "../output/rlpbwtPan_thr_noext_queries_tr.txt",
                               false);
    rlpbwtPan.match_tsv_tr_thr(queries_tr,
                               "../output/rlpbwtPan_thr_ext_queries_tr.txt",
                               true);

    std::cout << "\n\n";
    std::cout << "runs rlpbwt: " << rlpbwt.get_run_number() << "\n";
    std::cout << "runs rlpbwtbv: " << rlpbwtbv.get_run_number() << "\n";
    std::cout << "runs rlpbwtSlp: " << rlpbwtSlp.get_run_number() << "\n";
    std::cout << "runs rlpbwtPan: " << rlpbwtPan.get_run_number() << "\n";
    //BENCHMARK(BM_Query)->ArgName("test") ;
    //BENCHMARK_MAIN();
}

TEST(BigMatrix, Build) {
//    auto input = "../input_big/pbwt_matrix_10";
//    rlpbwt_ms<slp_panel_ra> rlpbwt(input, true, false,
//                                   "../output_big/pbwt_matrix_10.slp");
//    std::cout << "rlpbwt:\n" << rlpbwt.size_in_mega_bytes(true)
//              << " megabytes\n----\n";
//    std::cout << "extending...\n";
//    rlpbwt.extend();
//    std::cout << "extended...\n";
//    std::cout << "rlpbwt:\n" << rlpbwt.size_in_mega_bytes(true)
//              << " megabytes\n----\n";
//    auto filename = "../output_big/pbwt_matrix_10_s.ser";
//    std::ofstream file_o;
//    unsigned int size = 0;
//    file_o.open(filename);
//    rlpbwt.serialize(file_o);
//    file_o.close();
//    rlpbwt.match_tsv_tr_lce("../output_big/pbwt_matrix_10_query3",
//                             "../output_big/pbwt_matrix_10_matches3",
//                             false);
    std::ifstream file_i;
    file_i.open("../output_big/pbwt_matrix_10_s.ser");
    auto rlpbwt = new rlpbwt_ms<slp_panel_ra>();
    rlpbwt->load(file_i, "../output_big/pbwt_matrix_10.slp");
    file_i.close();

    rlpbwt->match_tsv_tr_lce("../output_big/pbwt_matrix_10_query3",
                             "../output_big/pbwt_matrix_10_le_matches3",
                             false, false);


    /*
    rlpbwt_naive rlpbwtn(input);
    auto filenamen = "../output_big/pbwt_matrix_10_n.ser";
    std::ofstream file_n;
    file_n.open(filenamen);
    rlpbwtn.serialize(file_n);
    file_n.close();
    rlpbwtn.match_tsv_tr("../output_big/pbwt_matrix_10_query3",
                        "../output_big/pbwt_matrix_10_n_matches3",
                        false);
    */
    /*
    rlpbwt_ms<panel_ra> rlpbwtp(input, true);
    auto filenamep = "../output_big/pbwt_matrix_10_p.ser";
    std::ofstream file_p;
    file_p.open(filenamep);
    rlpbwtp.serialize(file_p);
    file_p.close();
    rlpbwtp.match_tsv_tr_thr("../output_big/pbwt_matrix_10_query3",
                         "../output_big/pbwt_matrix_10_p_matches3",
                         false);
    rlpbwtp.match_tsv_tr_thr("../output_big/pbwt_matrix_10_query3",
                             "../output_big/pbwt_matrix_10_pe_matches3",
                             true);
    */
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    //::testing::GTEST_FLAG(filter) = "BuildRlpbwtSerBV*";
    //::testing::GTEST_FLAG(filter) = "BuildRlpbwtBVtest*";
    //::testing::GTEST_FLAG(filter) = "BuildRlpbwtBVVCF*";
    //::testing::GTEST_FLAG(filter) = "BuildRlpbwtThr*";
    //::testing::GTEST_FLAG(filter) = "BuildRlpbwtSerThr*";
    //::testing::GTEST_FLAG(filter) = "BuildRlpbwtNewThr*";
    //::testing::GTEST_FLAG(filter) = "RlpbwtRaTest*";
    //::testing::GTEST_FLAG(filter) = "Lce*";
    //::testing::GTEST_FLAG(filter) = "RlpbwtNaive*";
    //::testing::GTEST_FLAG(filter) = "MixRlpbwt*";
    //::testing::GTEST_FLAG(filter) = "Benchmark*";
    ::testing::GTEST_FLAG(filter) = "BigMatrix*";
    return RUN_ALL_TESTS();
}