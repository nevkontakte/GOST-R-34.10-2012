#ifndef NAF_H
#define NAF_H

#include <boost/multiprecision/cpp_int.hpp>

namespace gost_ecc {

namespace mp = ::boost::multiprecision;

template<unsigned window, typename integer_type>
unsigned naf(integer_type n, short table[]) {
    const unsigned short mask = (1 << window) - 1;
    const unsigned short shift = 1 << window;

    unsigned i;
    for (i = 0; n > 0; i++) {
        if (mp::bit_test(n, 0)) {
            table[i] = static_cast<short>(n & mask);
            if (table[i] > (1 << (window - 1))) {
                table[i] -= shift;
            }
            n = static_cast<integer_type>(n - table[i]);
        } else {
            table[i] = 0;
        }
        n >>= 1;
    }

    return i;
}

}

#endif // NAF_H
