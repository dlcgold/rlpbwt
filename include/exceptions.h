//
// Created by dlcgold on 30/10/21.
//

#ifndef RLPBWT_EXCEPTIONS_H
#define RLPBWT_EXCEPTIONS_H

#include <exception>

/**
 * @brief exception file not found
 */
class FileNotFoundException: public std::exception {
  const char* what() const noexcept override {
    return "file not found";
  }
};

#endif //RLPBWT_EXCEPTIONS_H
