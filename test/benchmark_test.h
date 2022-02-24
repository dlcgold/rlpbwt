//
// Created by dlcgold on 22/02/22.
//

#ifndef RLPBWT_BENCHMARK_TEST_H
#define RLPBWT_BENCHMARK_TEST_H

#include <benchmark/benchmark.h>
#include "../include/exceptions.h"
#include "../include/rlpbwt_bv.h"
#include "../include/panel_ra.h"
#include "../include/rlpbwt_ms.h"
#include "../include/rlpbwt_naive.h"

static void BM_Query(benchmark::State& state) {
    std::cout << "bench\n";
    auto input = "../input/sample_new.txt";
    auto slp = "../input/sample.slp";
    rlpbwt_naive rlpbwt(input);
    // This code gets timed
    for (auto _: state) {
        rlpbwt.external_match("010010100011101");
    }
}


#endif //RLPBWT_BENCHMARK_TEST_H
