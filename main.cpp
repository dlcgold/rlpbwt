#include <iostream>
#include "include/rlpbwt.h"
#include "include/exceptions.h"

int main(int argc, char **argv) {
    if (argc != 2) {
        throw FileNotFoundException{};
    }
    rlpbwt rlpbwt(argv[1]);
    int count = 0;
    for (const auto& e: rlpbwt.cols) {
        std::cout << "Table " << count << " start with 0? " << e.zero_first
                  << "\n";
        for (const auto& r: e.rows) {
            std::cout << r.p << " " << r.perm_p << " " << r.next_perm << "\n";
        }
        std::cout << "------------------\n";
        count++;
    }
    std::cout << rlpbwt.search_row(9);
    /*for(int i = 0; i < 20; i++){
        std::cout << rlpbwt.search_row(i);
        std::cout << "\n";
    }*/
}