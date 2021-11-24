#include "include/rlpbwt.h"
#include "include/rlpbwt_rlrow.h"
#include "include/exceptions.h"
#include "gtest/gtest.h"


TEST (BuildRlpbwtTest, TestSize) {
    rlpbwt rlpbwt("../input/matrix3.txt");
    EXPECT_EQ (20, rlpbwt.heigth);
    EXPECT_EQ (12, rlpbwt.width);
    EXPECT_EQ (12, rlpbwt.cols.size());
}

TEST (BuildRlpbwtTest, TestColumnFour) {
    rlpbwt rlpbwt("../input/matrix3.txt");
    auto row1 = rlpbwt_rlrow(0, 0, 0, 0);
    auto row2 = rlpbwt_rlrow(11, 17, 5, 11);
    auto row3 = rlpbwt_rlrow(14, 11, 5, 19);
    EXPECT_TRUE(rlpbwt.cols[4].zero_first);
    EXPECT_EQ (row1, rlpbwt.cols[4].rows[0]);
    EXPECT_EQ (row2, rlpbwt.cols[4].rows[1]);
    EXPECT_EQ (row3, rlpbwt.cols[4].rows[2]);
}

TEST (QueryRlpbwtTest, TestRowNine) {
    rlpbwt rlpbwt("../input/matrix3.txt");
    EXPECT_EQ ("100001000011", rlpbwt.search_row(9));
}

TEST (BuildRlpbwtTestSimple, TestBuild) {
    rlpbwt rlpbwt("../input/matrix_simple.txt");
    //std::cout << rlpbwt.heigth << " " << rlpbwt.width;
    //int count = 0;
    for(auto c: rlpbwt.cols){
        //std::cout << "column: " << count << ", start with 0? " << c.zero_first << "\n";
        for(auto r: c.rows){
            std::cout << r << "\n";
        }
        //count ++;
        std::cout << "--------------\n";
    }
    EXPECT_EQ (7, rlpbwt.heigth);
    EXPECT_EQ (20, rlpbwt.width);
}

int main(int argc, char **argv) {
    /*rlpbwt rlpbwt("../input/matrix3.txt");
    int count = 0;
    for (const auto &e: rlpbwt.cols) {
        std::cout << "Table " << count << ":\n";
        std::cout << e;
        std::cout << "------------------\n";
        count++;
    }
    std::cout << rlpbwt.search_row(9);
    for(int i = 0; i < 20; i++){
        std::cout << rlpbwt.search_row(i);
        std::cout << "\n";
    }*/
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}