#include <iostream>
#include <gtest/gtest.h>
#include <gperftools/heap-profiler.h>
#include "../include/exceptions.h"
#include "../include/rlpbwt_bv.h"
#include "../include/panel_ra.h"
#include "../include/rlpbwt_ra.h"
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
    rlpbwt_ra<panel_ra> rlpbwtRa("../input/sample_new.txt", true);
    auto matches = rlpbwtRa.match_thr("010010100011101", true);
    std::cout << matches << "\n";
    rlpbwt_ra<slp_panel_ra> rlpbwtSlp("../input/sample_new.txt", true,
                                      false, "../input/sample.slp");

    matches = rlpbwtRa.match_thr("010010100011101", false);
    std::cout << matches << "\n";
    rlpbwtSlp.extend();

    matches = rlpbwtSlp.match_thr("010010100011101", false);
    std::cout << matches << "\n";
    std::cout << rlpbwtSlp.extended << "\n";
    auto filename = "../output/phi_pbwt.ser";
    std::ofstream file_o;
    file_o.open(filename);
    auto size = rlpbwtSlp.serialize(file_o);
    std::cout << size << "\n";
    file_o.close();
    std::ifstream file_i;
    file_i.open(filename);
    auto rlpbwt = new rlpbwt_ra<slp_panel_ra>();
    rlpbwt->load(file_i, "../input/sample.slp");

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
    rlpbwt_ra<slp_panel_ra> rlpbwtSlp("../input/sample_new.txt", true, false,
                                      "../input/sample.slp");
    //rlpbwtSlp.extend();
    auto matches = rlpbwtSlp.match_thr("010010100011101", true, false);
    std::cout << matches;
    std::cout << rlpbwtSlp.panel->size_in_bytes() << " "
              << rlpbwtSlp.panel->size_in_mega_bytes();
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
    rlpbwt_naive rlpbwt("../input/sample_new.txt");
    rlpbwt_bv rlpbwtbv("../input/sample_new.txt");
    rlpbwt_ra<slp_panel_ra> rlpbwtSlp("../input/sample_new.txt", true,
                                      false, "../input/sample.slp");
    rlpbwt_ra<panel_ra> rlpbwtPan("../input/sample_new.txt", true, false);
    std::cout << "rlpbwt:\n"<< rlpbwt.size_in_mega_bytes(true) << " megabytes\n----\n";
    std::cout << "rlpbwt_bv:\n"<< rlpbwtbv.size_in_mega_bytes(true) << " megabytes\n----\n";
    std::cout << "rlpbwt_slp:\n"<< rlpbwtSlp.size_in_mega_bytes(true) << " megabytes\n----\n";
    std::cout << "rlpbwt_pa:\n"<< rlpbwtPan.size_in_mega_bytes(true) << " megabytes\n----\n";

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
    ::testing::GTEST_FLAG(filter) = "MixRlpbwt*";
    return RUN_ALL_TESTS();
}