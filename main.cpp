#include <iostream>
#include <gtest/gtest.h>
#include <gperftools/heap-profiler.h>
#include "include/exceptions.h"
#include "include/rlpbwt.h"


TEST (BuildRlpbwtTest, TestBuildAndQuery) {
    HeapProfilerStart("heap.prof");
    std::cout << IsHeapProfilerRunning() << "\n";
    rlpbwt rlpbwt("../input/sample.txt");
    HeapProfilerDump("end construction");
    HeapProfilerStop();
    EXPECT_EQ(rlpbwt.heigth, 20);
    EXPECT_EQ(rlpbwt.width, 15);
    std::cout << rlpbwt.heigth << " " << rlpbwt.width << "\n";
    bool verbose = true;
    if (verbose) {
        int count = 0;
        for (const auto &c: rlpbwt.cols) {
            std::string z;
            if (c.zero_first) {
                z = "yes";
            } else {
                z = "no";
            }
            std::cout << "column: " << count << "\nstart with 0? "
                      << z << ", c: " << c.count_0 << "\n";

            for (const auto &r: c.rows) {
                std::cout << r << "\n";
            }
            for (auto d: c.div) {
                std::cout << d << " ";
            }
            count++;
            std::cout << "\n--------------\n";
        }
    }
    auto matches = rlpbwt.external_match("010010100011101");
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
    for (const auto &c: rlpbwt.cols) {
        rlsize += sizeof(c.zero_first) * 10E-6;
        rlsize += sizeof(c.count_0) * 10E-6;
        rlsize += sizeof(unsigned int) * (double) c.rows.size() * 10E-6;
        rlsize += sizeof(unsigned int) * (double) c.rows.size() * 10E-6;
        rlsize += sdsl::size_in_mega_bytes(c.div);
    }
    double nrlsize = sizeof(unsigned int) *
                     (double) (rlpbwt.heigth * rlpbwt.width * 5) * 10E-6;
    nrlsize += sizeof(unsigned int) * (double) rlpbwt.width * 10E-6;
    std::cout << rlsize << " vs " << nrlsize << "\n";
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}