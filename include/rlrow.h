//
// Created by dlcgold on 17/01/22.
//

#ifndef RLPBWT_RLROW_H
#define RLPBWT_RLROW_H


#include <ostream>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/vector.hpp>

class rlrow {
public:
    rlrow(unsigned int p, unsigned int uv);

    unsigned int p;
    unsigned int uv;
    friend std::ostream &operator<<(std::ostream &os, const rlrow &rlrow);

    virtual ~rlrow();

private:
    friend class boost::serialization::access;

};

namespace boost {
    namespace serialization {
        template<class Archive>
        void serialize(Archive &a, rlrow &e,
                       const unsigned version){
            a & e.p & e.uv;
        }
    }
}

#endif //RLPBWT_RLROW_H
