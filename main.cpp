#include <iostream>
#include <gtest/gtest.h>
#include <gperftools/heap-profiler.h>
#include <boost/archive/text_oarchive.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include "include/exceptions.h"
#include "include/rlpbwt.h"
#include "include/rlpbwtbv.h"
#include "include/birlpbwt.h"


TEST (BuildRlpbwtTest, TestBuildAndQuery) {
    //HeapProfilerStart("heap.prof");
    //std::cout << IsHeapProfilerRunning() << "\n";
    rlpbwt rlpbwt("../input/sample.txt", false);
    //HeapProfilerDump("end construction");
    //HeapProfilerStop();
    std::cout << rlpbwt.cols[rlpbwt.cols.size() - 1] << "\n";
    EXPECT_EQ(rlpbwt.heigth, 20);
    EXPECT_EQ(rlpbwt.width, 15);
    std::cout << rlpbwt.heigth << " " << rlpbwt.width << "\n";
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
    rlsize += sizeof(rlpbwt.heigth) * 10E-6;
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
                     (double) (rlpbwt.heigth * rlpbwt.width * 5) * 10E-6;
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
    rlsizeb += sizeof(birlpbwt.frlpbwt.heigth) * 10E-6;
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
    double nrlsize = sizeof(unsigned int) * (double) (birlpbwt.frlpbwt.heigth *
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
    std::ofstream outfile("rlpbwt.ser");
    boost::archive::text_oarchive archive(outfile);
    archive << rlpbwt;

    std::cout << rlpbwt.width << " " << rlpbwt.heigth << "\n";
    //rlpbwt.print();
    EXPECT_EQ(rlpbwt.heigth, 900);
    EXPECT_EQ(rlpbwt.width, 500);
    clock_t START = clock();
    // TODO add test check for matches
    rlpbwt.external_match_vcf("../input/sample_query.vcf", 255);
    std::cout << clock() - START << " time\n";
}


TEST (BuildBiRlpbwtVCF, TestBuild) {
    //HeapProfilerStart("birlpbwtheap.prof");
    //std::cout << IsHeapProfilerRunning() << "\n";
    birlpbwt birlpbwt("../input/sample_panel.vcf", true);
    //HeapProfilerDump("end construction");
    //HeapProfilerStop();

    std::cout << birlpbwt.frlpbwt.width << " " << birlpbwt.frlpbwt.heigth
              << "\n";
    EXPECT_EQ(birlpbwt.frlpbwt.heigth, 900);
    EXPECT_EQ(birlpbwt.frlpbwt.width, 500);
    //birlpbwt.print()
    unsigned int size = 0;
    unsigned int sizeb = 0;
    unsigned int obj = 0;
    unsigned int objb = 0;
    auto rlsizeb = sizeof(birlpbwt.frlpbwt.width) * 10E-6;
    rlsizeb += sizeof(birlpbwt.frlpbwt.heigth) * 10E-6;
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
    double nrlsize = sizeof(unsigned int) * (double) (birlpbwt.frlpbwt.heigth *
                                                      birlpbwt.frlpbwt.width *
                                                      5) * 10E-6;
    std::cout << rlsizeb << " vs " << nrlsize << "\n";
    std::cout << "saved: " << obj + objb << " vs "
              << (birlpbwt.frlpbwt.heigth * birlpbwt.frlpbwt.width * 5) + 1
              << " (~"
              << ((birlpbwt.frlpbwt.heigth * birlpbwt.frlpbwt.width * 5) + 1) /
                 (obj + objb) << " times less) and avg height "
              << (size + sizeb) / (birlpbwt.frlpbwt.width * 2) << "\n";
}

TEST (BuildBiRlpbwtVCF, TestQuery) {
    birlpbwt birlpbwt("../input/sample_panel.vcf", true);
    std::ofstream outfile("birlpbwt.ser");
    boost::archive::text_oarchive archive(outfile);
    archive << birlpbwt;
    EXPECT_EQ(birlpbwt.frlpbwt.heigth, 900);
    EXPECT_EQ(birlpbwt.frlpbwt.width, 500);
    clock_t START = clock();
    // TODO add test check for matches
    birlpbwt.external_match_vcf("../input/sample_query.vcf", 255, false);
    std::cout << clock() - START << " time\n";
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
    EXPECT_EQ(rlpbwtbv.heigth, 900);
    EXPECT_EQ(rlpbwtbv.width, 500);
    clock_t START = clock();
    // TODO add test check for matches
    rlpbwtbv.external_match_vcf("../input/sample_query.vcf", 255, false);
    std::cout << clock() - START << " time\n";
}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    //::testing::GTEST_FLAG(filter) = "BuildRlpbwtTest*";
    //::testing::GTEST_FLAG(filter) = "BuildBiRlpbwtTest*";
    //::testing::GTEST_FLAG(filter) = "BuildRlpbwtVCF*";
    //::testing::GTEST_FLAG(filter) = "BuildBiRlpbwtVCF*";
    ::testing::GTEST_FLAG(filter) = "BuildRlpbwtBVtest*";
    //::testing::GTEST_FLAG(filter) = "BuildRlpbwtBVVCF*";
    return RUN_ALL_TESTS();
}