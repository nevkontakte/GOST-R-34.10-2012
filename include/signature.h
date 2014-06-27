#ifndef SIGNATURE_H
#define SIGNATURE_H

#include <sign_engine.h>

#include <cryptopp/ecp.h>
#include <cstdint>

namespace gost_ecc {

class signature
{
    ::CryptoPP::ECP curve;
    ::CryptoPP::ModularArithmetic subgroup;
    ::CryptoPP::ECPPoint basePoint;
public:
    signature(u_int64_t (&modulus)[8], u_int64_t (&a)[8], u_int64_t (&b)[8],
              u_int64_t (&subgroupModulus)[8],
              u_int64_t (&base_x)[8], u_int64_t (&base_y)[8]);

    Gost12S512Status sign(const byte* private_key, const byte* rand, const byte* hash, byte* signature);
    Gost12S512Status verify(const byte* public_key_x, const byte* public_key_y, const byte* hash, const byte* signature);
};

}

#endif // SIGNATURE_H
