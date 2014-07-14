#ifndef SIGNATURE_H
#define SIGNATURE_H

#include <sign_engine.h>

#include <elliptic_curve.h>
#include <cstdint>

namespace gost_ecc {

using byte = unsigned char;

class signature
{
    using ec = elliptic_curve<mp::uint512_t, mp::uint1024_t, cpp_int_fixed<528>>;
    using pf = prime_field<mp::uint512_t, mp::uint1024_t>;

    static const unsigned comb_window = 10;
    static const unsigned dynamic_naf_window = 6;
    static const unsigned static_naf_window = 10;

    ec curve;
    pf subgroup;
    ec::point basePoint;
    ec::jacobian_point basePointTable[1 << comb_window];
    ec::jacobian_point basePointNafTable[1 << (static_naf_window - 2)];

public:
    signature(u_int64_t (&modulus)[8], u_int64_t (&a)[8], u_int64_t (&b)[8],
              u_int64_t (&subgroupModulus)[8],
              u_int64_t (&base_x)[8], u_int64_t (&base_y)[8]);

    Gost12S512Status sign(const byte* private_key, const byte* rand, const byte* hash, byte* signature);
    Gost12S512Status verify(const byte* public_key_x, const byte* public_key_y, const byte* hash, const byte* signature);
};

}

#endif // SIGNATURE_H
