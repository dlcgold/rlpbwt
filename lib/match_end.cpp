//
// Created by dlcgold on 27/01/22.
//

#include "../include/match_end.h"

#include <utility>

match_end::~match_end() = default;

match_end::match_end(unsigned int begin, unsigned int anEnd,
                     std::pair<unsigned int, unsigned int> interval)
        : begin(begin), end(anEnd), interval(std::move(interval)) {}

bool match_end::operator==(const match_end &rhs) const {
    return begin == rhs.begin &&
           end == rhs.end &&
           interval == rhs.interval;
}

bool match_end::operator!=(const match_end &rhs) const {
    return !(rhs == *this);
}

std::ostream &operator<<(std::ostream &os, const match_end &anEnd) {
    os << "begin: " << anEnd.begin << " end: " << anEnd.end <<
       " interval: [" << anEnd.interval.first << ", "
       << anEnd.interval.second << ")";
    return os;
}

unsigned int match_end::get_length() {
    return this->interval.second - this->interval.first;
}
