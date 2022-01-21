#include <iostream>
#include <gtest/gtest.h>
#include <gperftools/heap-profiler.h>
#include "include/exceptions.h"
#include "include/rlpbwt.h"
#include "include/birlpbwt.h"


TEST (BuildRlpbwtTest, TestBuildAndQuery) {
    //HeapProfilerStart("heap.prof");
    //std::cout << IsHeapProfilerRunning() << "\n";
    rlpbwt rlpbwt("../input/sample.txt");
    //HeapProfilerDump("end construction");
    //HeapProfilerStop();
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
        rlsize += sizeof(c.zero_first) * 10E-6;
        rlsize += sizeof(c.count_0) * 10E-6;
        rlsize += sizeof(unsigned int) * (double) c.rows.size() * 10E-6;
        rlsize += sizeof(unsigned int) * (double) c.rows.size() * 10E-6;
        rlsizeb += sizeof(c.zero_first) * 10E-6;
        rlsizeb += sizeof(c.count_0) * 10E-6;
        rlsizeb += sizeof(unsigned int) * (double) c.rows.size() * 10E-6;
        rlsizeb += sizeof(unsigned int) * (double) c.rows.size() * 10E-6;
        rlsize += sdsl::size_in_mega_bytes(c.lcp);
    }
    double nrlsize = sizeof(unsigned int) *
                     (double) (rlpbwt.heigth * rlpbwt.width * 5) * 10E-6;
    nrlsize += sizeof(unsigned int) * (double) rlpbwt.width * 10E-6;
    std::cout << rlsize << " vs " << rlsizeb << " vs " << nrlsize << "\n";
}

TEST (BuildBiRlpbwtTest, TestBuildAndQuery) {
    birlpbwt birlpbwt("../input/sample.txt");
    bool verbose = false;
    if (verbose) {
        birlpbwt.print();
    }
    auto matches = birlpbwt.external_match("010010100011101", false);
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
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    //::testing::GTEST_FLAG(filter) = "BuildBiRlpbwtTest*";
    //::testing::GTEST_FLAG(filter) = "BuildRlpbwtTest*";

    return RUN_ALL_TESTS();
}