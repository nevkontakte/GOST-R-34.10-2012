#ifndef NAF_H
#define NAF_H

#include <boost/multiprecision/cpp_int.hpp>

namespace gost_ecc {

namespace mp = ::boost::multiprecision;

template<unsigned window, typename integer_type>
class naf {
protected:
    integer_type n;
    const unsigned mask;
    const unsigned shift;

public:
    naf(const integer_type n)
        :n(n), mask( (1 << window) - 1 ), shift( 1 << window )
    {}

    int next() {
        if (n > 0) {
            int ki;
            if (mp::bit_test(this->n, 0)) {
                ki = static_cast<int>(this->n & this->mask);
                if (ki > (1 << (window - 1))) {
                    ki -= shift;
                }
                this->n = static_cast<integer_type>(this->n - ki);
            } else {
                ki = 0;
            }
            this->n >>= 1;

            return ki;
        } else {
            return 0;
        }
    }
};

}

#endif // NAF_H
