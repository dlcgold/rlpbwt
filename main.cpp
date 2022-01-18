
#include "include/rlpbwt_rlrow.h"
#include "include/exceptions.h"
#include "gtest/gtest.h"
#include "include/rlpbwtc.h"
#include <gperftools/heap-profiler.h>

//TEST (BuildRlpbwtTest, TestSize) {
//    rlpbwt rlpbwt("../input/matrix3.txt");
//    EXPECT_EQ (20, rlpbwt.heigth);
//    EXPECT_EQ (12, rlpbwt.width);
//    EXPECT_EQ (12, rlpbwt.cols.size());
//}
//
//TEST (BuildRlpbwtTest, TestColumnFour) {
//    rlpbwt rlpbwt("../input/matrix3.txt");
//    auto row1 = rlpbwt_rlrow(0, 0, 0, 0);
//    auto row2 = rlpbwt_rlrow(11, 17, 5, 11);
//    auto row3 = rlpbwt_rlrow(14, 11, 5, 19);
//    EXPECT_TRUE(rlpbwt.cols[4].zero_first);
//    EXPECT_EQ (row1, rlpbwt.cols[4].rows[0]);
//    EXPECT_EQ (row2, rlpbwt.cols[4].rows[1]);
//    EXPECT_EQ (row3, rlpbwt.cols[4].rows[2]);
//}
//
//TEST (QueryRlpbwtTest, TestRowNine) {
//    rlpbwt rlpbwt("../input/matrix3.txt");
//    EXPECT_EQ ("100001000011", rlpbwt.search_row(9));
//}
//
//TEST (BuildRlpbwtTestSimple, TestBuild) {
//    rlpbwt rlpbwt("../input/matrix_simple.txt");
//    //std::cout << rlpbwt.heigth << " " << rlpbwt.width;
//    bool verbose = false;
//    if (verbose) {
//        int count = 0;
//        for (const auto& c: rlpbwt.cols) {
//            std::cout << "column: " << count << ", start with 0? "
//                      << c.zero_first << "\n";
//            for (const auto& r: c.rows) {
//                std::cout << r << "\n";
//            }
//            count++;
//            std::cout << "--------------\n";
//        }
//    }
//    EXPECT_EQ (7, rlpbwt.heigth);
//    EXPECT_EQ (20, rlpbwt.width);
//    EXPECT_EQ("11010011111001001001", rlpbwt.search_row(0));
//    EXPECT_EQ("01000011111001110010", rlpbwt.search_row(1));
//    EXPECT_EQ("00010000011110001010", rlpbwt.search_row(2));
//    EXPECT_EQ("10011010100011100010", rlpbwt.search_row(3));
//    EXPECT_EQ("01101111100100111100", rlpbwt.search_row(4));
//    EXPECT_EQ("11001010101010001111", rlpbwt.search_row(5));
//    EXPECT_EQ("00010111111100100011", rlpbwt.search_row(6));
//    rlpbwt.external_match("00101111111001000111");
//}

TEST (BuildRlpbwtTestSample, TestBuild) {

    //std::cout << rlpbwt.search_row(4);
    //rlpbwt.external_match("010010100011101");
    HeapProfilerStart("heap.prof");
    std::cout << IsHeapProfilerRunning() << "\n";
    rlpbwtc rlpbwtc("../input/sample.txt");
    HeapProfilerDump("end construction");
    HeapProfilerStop();
    std::cout << rlpbwtc.heigth << " " << rlpbwtc.width << "\n";
    bool verbose = true;
    if (verbose) {
        int count = 0;
        for (const auto &c: rlpbwtc.cols) {
            std::cout << "column: " << count << ", start with 0? "
                      << c.zero_first << ", zeros? " << c.count_0 << "\n";

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
    auto matches = rlpbwtc.ematch("010010100011101");
    for (const auto &m: matches) {
        std::cout << m << "\n";
    }
    auto rlsize = sizeof(rlpbwtc.width) * 10E-6;
    rlsize += sizeof(rlpbwtc.heigth) * 10E-6;
    for (auto c: rlpbwtc.cols) {
        rlsize += sizeof(c.zero_first) * 10E-6;
        rlsize += sizeof(c.count_0) * 10E-6;
        rlsize += sizeof(unsigned int) * (double) c.rows.size() * 10E-6;
        rlsize += sizeof(unsigned int) * (double) c.rows.size() * 10E-6;
        rlsize += sdsl::size_in_mega_bytes(c.div);
    }
    double nrlsize = sizeof(unsigned int) *
                     (double) (rlpbwtc.heigth * rlpbwtc.width * 5) * 10E-6;
    std::cout << rlsize << " vs " << nrlsize << "\n";
}

//TEST (BuildRlpbwtTestSimple, TestSizeZero) {
//    rlpbwt rlpbwt("../input/simple_matrix_zero.txt");
//    EXPECT_EQ (7, rlpbwt.heigth);
//    EXPECT_EQ (22, rlpbwt.width);
//    EXPECT_EQ("1101001111100100100010", rlpbwt.search_row(0));
//    EXPECT_EQ("0100001111100111000100", rlpbwt.search_row(1));
//    EXPECT_EQ("0001000001111000100100", rlpbwt.search_row(2));
//    EXPECT_EQ("1001101010001110000100", rlpbwt.search_row(3));
//    EXPECT_EQ("0110111110010011110000", rlpbwt.search_row(4));
//    EXPECT_EQ("1100101010101000110110", rlpbwt.search_row(5));
//    EXPECT_EQ("0001011111110010000110", rlpbwt.search_row(6));
//}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}