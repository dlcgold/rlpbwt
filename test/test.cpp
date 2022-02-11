#include <iostream>
#include <gtest/gtest.h>
#include <gperftools/heap-profiler.h>
#include <sdsl/int_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include <SelfShapedSlp.hpp>
#include <DirectAccessibleGammaCode.hpp>
#include <SelectType.hpp>
#include "../include/exceptions.h"
#include "../include/rlpbwt.h"
#include "../include/rlpbwtbv.h"
#include "../include/birlpbwt.h"
#include "../include/rlpbwt_thr.h"
#include "../include/panel_ra.h"
//#include "../include/slp_panel_ra.h"


TEST (BuildRlpbwtTest, TestBuildAndQuery) {
    //HeapProfilerStart("heap.prof");
    //std::cout << IsHeapProfilerRunning() << "\n";
    rlpbwt rlpbwt("../input/sample.txt", false);
    //HeapProfilerDump("end construction");
    //HeapProfilerStop();
    std::cout << rlpbwt.cols[rlpbwt.cols.size() - 1] << "\n";
    EXPECT_EQ(rlpbwt.height, 20);
    EXPECT_EQ(rlpbwt.width, 15);
    std::cout << rlpbwt.height << " " << rlpbwt.width << "\n";
    bool verbose = false;
    if (verbose) {
        rlpbwt.print();
    }
    auto matches = rlpbwt.external_match("010010100011101", true);
    for (const auto &m: matches) {
        std::cout << m << "\n";
    }
    auto match0 = match(0, 5, 4);
    auto match1 = match(3, 9, 1);
    auto match2 = match(7, 11, 1);
    auto match3 = match(11, 14, 3);
    EXPECT_EQ(matches[0], match0);
    EXPECT_EQ(matches[1], match1);
    EXPECT_EQ(matches[2], match2);
    EXPECT_EQ(matches[3], match3);

    auto rlsize = sizeof(rlpbwt.width) * 10E-6;
    rlsize += sizeof(rlpbwt.height) * 10E-6;
    auto rlsizeb = rlsize;
    for (const auto &c: rlpbwt.cols) {
        rlsize += sizeof(bool) * 10E-6;
        rlsize += sizeof(unsigned int) * 10E-6;
        rlsize += sizeof(unsigned int) * (double) c.rows.size() * 10E-6;
        rlsize += sizeof(unsigned int) * (double) c.rows.size() * 10E-6;
        rlsizeb += sizeof(bool) * 10E-6;
        rlsizeb += sizeof(unsigned int) * 10E-6;
        rlsizeb += sizeof(unsigned int) * (double) c.rows.size() * 10E-6;
        rlsizeb += sizeof(unsigned int) * (double) c.rows.size() * 10E-6;
        rlsize += sdsl::size_in_mega_bytes(c.lcp);
    }
    double nrlsize = sizeof(unsigned int) *
                     (double) (rlpbwt.height * rlpbwt.width * 5) * 10E-6;
    nrlsize += sizeof(unsigned int) * (double) rlpbwt.width * 10E-6;
    std::cout << rlsize << " vs " << rlsizeb * 2 << " vs " << nrlsize << "\n";
}

TEST (BuildBiRlpbwtTest, TestBuildAndQuery) {
    birlpbwt birlpbwt("../input/sample2.txt", false);
    bool verbose = false;
    if (verbose) {
        birlpbwt.print();
    }
    auto matches = birlpbwt.external_match("010010100011101");
    for (const auto &m: matches) {
        std::cout << m << "\n";
    }
    /*auto match0 = match(0, 5, 4);
    auto match1 = match(3, 9, 1);
    auto match2 = match(7, 11, 1);
    auto match3 = match(11, 14, 3);
    EXPECT_EQ(matches[0], match0);
    EXPECT_EQ(matches[1], match1);
    EXPECT_EQ(matches[2], match2);
    EXPECT_EQ(matches[3], match3);
    auto rlsizeb = sizeof(birlpbwt.frlpbwt.width) * 10E-6;
    rlsizeb += sizeof(birlpbwt.frlpbwt.height) * 10E-6;
    rlsizeb *= 2;
    for (const auto &c: birlpbwt.frlpbwt.cols) {
        rlsizeb += sizeof(bool) * 10E-6;
        rlsizeb += sizeof(unsigned int) * 10E-6;
        rlsizeb += sizeof(unsigned int) * (double) c.rows.size() * 10E-6;
        rlsizeb += sizeof(unsigned int) * (double) c.rows.size() * 10E-6;
    }
    for (const auto &c: birlpbwt.brlpbwt.cols) {
        rlsizeb += sizeof(bool) * 10E-6;
        rlsizeb += sizeof(unsigned int) * 10E-6;
        rlsizeb += sizeof(unsigned int) * (double) c.rows.size() * 10E-6;
        rlsizeb += sizeof(unsigned int) * (double) c.rows.size() * 10E-6;
    }
    double nrlsize = sizeof(unsigned int) * (double) (birlpbwt.frlpbwt.height *
                                                      birlpbwt.frlpbwt.width *
                                                      5) * 10E-6;
    std::cout << rlsizeb << " vs " << nrlsize << "\n";*/
}

TEST (BuildRlpbwtVCF, TestBuildAndQuery) {
    //HeapProfilerStart("heaprlp2.prof");
    //std::cout << IsHeapProfilerRunning() << "\n";
    rlpbwt rlpbwt("../input/sample_panel.vcf", true);
    //HeapProfilerDump("end construction");
    //HeapProfilerStop();

    std::cout << rlpbwt.width << " " << rlpbwt.height << "\n";
    //rlpbwt.print();
    EXPECT_EQ(rlpbwt.height, 900);
    EXPECT_EQ(rlpbwt.width, 500);
    clock_t START = clock();
    // TODO add test check for matches
    rlpbwt.external_match_vcf("../input/sample_query.vcf", 255);
    std::cout << clock() - START << " time\n";
}

TEST (BuildRlpbwtSerOrig, TestSer) {
    auto rlpbwt_t = new rlpbwt("../input/sample2.txt", false, false);

    auto filename = "../output/sample2orig_panel.ser";
    std::ofstream file_o;
    std::ifstream file_i;
    unsigned int size = 0;
    file_o.open(filename);
    size = rlpbwt_t->cols[14].serialize(file_o);
    std::cout << size << "\n";
    file_o.close();
    file_i.open(filename);
    auto col = new column();
    col->load(file_i);
    for (const auto &r: col->rows) {
        std::cout << r.p << ", " << r.uv << "\n";
    }
    for (const auto &e: col->lcp) {
        std::cout << e << " ";
    }
    std::cout << "\n";
    file_i.close();


    filename = "../output/sample2orig_pbwt.ser";
    file_o.open(filename);
    size = rlpbwt_t->serialize(file_o);
    std::cout << size << "\n";
    file_o.close();
    // TODO FIX THIS
    //file_i.open(filename);
    //auto rlpbwt = new rlpbwtbv();
    //rlpbwt->load(file_i);
    /*std::cout << rlpbwt->cols[13] << "\n";
    std::cout << rlpbwt->cols.size() << " vs " << rlpbwt_t->cols.size() << "\n";
    file_i.close();*/
}

TEST (BuildBiRlpbwtVCF, TestBuild) {
    //HeapProfilerStart("birlpbwtheap.prof");
    //std::cout << IsHeapProfilerRunning() << "\n";
    birlpbwt birlpbwt("../input/sample_panel.vcf", true);
    //HeapProfilerDump("end construction");
    //HeapProfilerStop();

    std::cout << birlpbwt.frlpbwt.width << " " << birlpbwt.frlpbwt.height
              << "\n";
    EXPECT_EQ(birlpbwt.frlpbwt.height, 900);
    EXPECT_EQ(birlpbwt.frlpbwt.width, 500);
    //birlpbwt.print()
    unsigned int size = 0;
    unsigned int sizeb = 0;
    unsigned int obj = 0;
    unsigned int objb = 0;
    auto rlsizeb = sizeof(birlpbwt.frlpbwt.width) * 10E-6;
    rlsizeb += sizeof(birlpbwt.frlpbwt.height) * 10E-6;
    rlsizeb *= 2;
    for (const auto &c: birlpbwt.frlpbwt.cols) {
        rlsizeb += sizeof(bool) * 10E-6;
        rlsizeb += sizeof(unsigned int) * 10E-6;
        rlsizeb += sizeof(unsigned int) * (double) c.rows.size() * 10E-6;
        rlsizeb += sizeof(unsigned int) * (double) c.rows.size() * 10E-6;
        size += c.rows.size();
        obj += (c.rows.size() * 2 + 2);
    }
    for (const auto &c: birlpbwt.brlpbwt.cols) {
        rlsizeb += sizeof(bool) * 10E-6;
        rlsizeb += sizeof(unsigned int) * 10E-6;
        rlsizeb += sizeof(unsigned int) * (double) c.rows.size() * 10E-6;
        rlsizeb += sizeof(unsigned int) * (double) c.rows.size() * 10E-6;
        sizeb += c.rows.size();
        objb += (c.rows.size() * 2 + 2);
    }
    double nrlsize = sizeof(unsigned int) * (double) (birlpbwt.frlpbwt.height *
                                                      birlpbwt.frlpbwt.width *
                                                      5) * 10E-6;
    std::cout << rlsizeb << " vs " << nrlsize << "\n";
    std::cout << "saved: " << obj + objb << " vs "
              << (birlpbwt.frlpbwt.height * birlpbwt.frlpbwt.width * 5) + 1
              << " (~"
              << ((birlpbwt.frlpbwt.height * birlpbwt.frlpbwt.width * 5) + 1) /
                 (obj + objb) << " times less) and avg height "
              << (size + sizeb) / (birlpbwt.frlpbwt.width * 2) << "\n";
}

TEST (BuildBiRlpbwtVCF, TestQuery) {
    birlpbwt birlpbwt("../input/sample_panel.vcf", true);
    EXPECT_EQ(birlpbwt.frlpbwt.height, 900);
    EXPECT_EQ(birlpbwt.frlpbwt.width, 500);
    clock_t START = clock();
    // TODO add test check for matches
    birlpbwt.external_match_vcf("../input/sample_query.vcf", 255, false);
    std::cout << clock() - START << " time\n";
}

TEST (BuildRlpbwtSerBV, TestSer) {
    auto rlpbwt_t = new rlpbwtbv("../input/sample2.txt", false, false);

    auto filename = "../output/sample2bv_panel.ser";
    std::ofstream file_o;
    std::ifstream file_i;
    unsigned int size = 0;
    file_o.open(filename);
    size = rlpbwt_t->cols[0].serialize(file_o);
    std::cout << size << "\n";
    file_o.close();
    file_i.open(filename);
    auto col = new columnbv();
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
    auto rlpbwt = new rlpbwtbv();
    rlpbwt->load(file_i);
    std::cout << rlpbwt->cols[13] << "\n";
    std::cout << rlpbwt->cols.size() << " vs " << rlpbwt_t->cols.size() << "\n";
    file_i.close();
}


TEST (BuildRlpbwtBVtest, TestBuildQuery) {
    rlpbwtbv rlpbwtbv("../input/sample2.txt", false);

    auto matches = rlpbwtbv.external_match("010010100011101", 1);
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
    rlpbwtbv rlpbwtbv("../input/sample_panel.vcf", true);
    EXPECT_EQ(rlpbwtbv.height, 900);
    EXPECT_EQ(rlpbwtbv.width, 500);
    clock_t START = clock();
    // TODO add test check for matches
    rlpbwtbv.external_match_vcf("../input/sample_query.vcf", 255, false);
    std::cout << clock() - START << " time\n";
}

TEST (BuildRlpbwtNewThr, TestBuildQuery) {
    rlpbwt_thr rlpbwt_thr("../input/sample_new.txt", 15, 19, false);
    std::cout << rlpbwt_thr.cols.size() << "\n";
    unsigned int count = 0;
    /*for (const auto &c: rlpbwt_thr.cols) {
        std::cout << "column " << count << ":\n";
        std::cout << c;
        count++;
        std::cout << "--------------------\n";
    }
     */
    std::cout << rlpbwt_thr.panelbv;

    rlpbwt_thr.match_thr("010010100011101", false);

}

TEST (BuildRlpbwtThr, TestBuildQuery) {
    rlpbwt_thr rlpbwt_thr("../input/sample2.txt", false, false);
    std::cout << rlpbwt_thr.panelbv.w << " ," << rlpbwt_thr.panelbv.h << "\n";
    std::cout << rlpbwt_thr.cols.size();
    /*
     unsigned int count = 0;
     for (const auto &c: rlpbwt_thr.cols) {
         std::cout << "column " << count << ":\n";
         std::cout << c;
         count++;
         std::cout << "--------------------\n";
     }
     */
    rlpbwt_thr.match_thr("010010100011101", false);

}

TEST (BuildRlpbwtSerThr, TestSer) {
    auto rlpbwt_t = new rlpbwt_thr("../input/sample2.txt", false, false);

    auto filename = "../output/sample2_panel.ser";
    std::ofstream file_o;
    file_o.open(filename);
    auto size = rlpbwt_t->panelbv.serialize(file_o);
    std::cout << size << "\n";
    file_o.close();
    auto panel = new panel_ra();
    std::ifstream file_i;
    file_i.open(filename);
    panel->load(file_i);
    std::cout << *panel << "\n";
    file_i.close();

    filename = "../output/sample2_col.ser";
    file_o.open(filename);
    size = rlpbwt_t->cols[0].serialize(file_o);
    std::cout << size << "\n";
    file_o.close();
    file_i.open(filename);
    auto col = new column_thr();
    col->load(file_i);
    std::cout << col->runs << "\n";
    std::cout << col->rank_runs(18) << " " << col->select_runs(2) << "\n";
    file_i.close();

    filename = "../output/sample2_pbwt.ser";
    file_o.open(filename);
    size = rlpbwt_t->serialize(file_o);
    std::cout << size << "\n";
    file_o.close();
    file_i.open(filename);
    auto rlpbwt = new rlpbwt_thr();
    rlpbwt->load(file_i);
    std::cout << rlpbwt->panelbv << "\n";
    std::cout << rlpbwt->cols[13] << "\n";
    std::cout << rlpbwt->cols.size() << " vs " << rlpbwt_t->cols.size() << "\n";
    file_i.close();
}

TEST(Slp, TestLoad) {
    //auto slp = new slp_panel_ra("../input/sample.slp", 20, 15);
   // std::cout << slp << "\n";
    using SelSd = SelectSdvec<>;
    using DagcSd = DirectAccessibleGammaCode<SelSd>;
    using shaped_slp_t = SelfShapedSlp<uint32_t, DagcSd, DagcSd, SelSd>;
    shaped_slp_t panel;
    std::ifstream in("../input/sample.slp");
    panel.load(in);
    std::cout << panel.charAt(0) << "\n";
    in.close();
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
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
    ::testing::GTEST_FLAG(filter) = "Slp*";
    return RUN_ALL_TESTS();
}