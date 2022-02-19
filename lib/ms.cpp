//
// Created by dlcgold on 19/02/22.
//

#include "../include/ms.h"

#include <utility>

ms::ms() = default;


ms::ms(std::vector<unsigned int> pos,
       std::vector<unsigned int> len,
       std::vector<std::pair<unsigned int, unsigned int>> matches) : pos(std::move(
        pos)), len(std::move(len)), matches(std::move(matches)) {}
