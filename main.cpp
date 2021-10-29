#include <iostream>
#include "include/rlpbwt.h"

int main(int argc, char** argv) {
    rlpbwt rlpbwt(argv[1]);
    int count = 0;
    for(auto e: rlpbwt.cols){
        std::cout << "Table " << count << " start with 0? "<< e.zero_first << "\n";
        for(auto r: e.rows){
            std::cout<< r.p << " " <<r.perm_p<<" "<<r.next_perm <<"\n";
        }
        std::cout << "------------------\n";
        count++;
    }
    //std::cout << rlpbwt.search_row(9);

}