#include <iostream>
#include <sdsl/int_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include "include/exceptions.h"
#include "include/rlpbwt_thr.h"

template<typename rlpbwt_t>
void build(std::string in_filename, const std::string &out_filename) {

    rlpbwt_t rlpbwt(in_filename.c_str(), false);
    std::ofstream file_o(out_filename);
    size_t size = rlpbwt.serialize(file_o);
    std::cout << size << "\n";
    file_o.close();
}

int main(int argc, char **argv) {

    std::string in_filename("../input/sample_new.txt");
    std::string out_filename("../output/sample2.txt.pbwt");
    

    build<rlpbwt_thr>(in_filename, out_filename);

    return 0;
}