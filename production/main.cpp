#include <curve.h>
#include <signature.h>

#include <iostream>
#include <algorithm>
#include <sign_engine.h>
#include <cryptopp/ecp.h>

using namespace CryptoPP;

::gost_ecc::signature* s;

Gost12S512Status Gost12S512Init() {
    s = new ::gost_ecc::signature(::gost_ecc::p, ::gost_ecc::a, ::gost_ecc::b, ::gost_ecc::q, ::gost_ecc::x0, ::gost_ecc::y0);
    return kStatusOk;
}

Gost12S512Status Gost12S512Sign( const char* privateKey,
                                 const char* rand,
                                 const char* hash,
                                 char* signature ) {
    return s->sign(reinterpret_cast<const byte*>(privateKey),
                   reinterpret_cast<const byte*>(rand),
                   reinterpret_cast<const byte*>(hash),
                   reinterpret_cast<byte*>(signature));
}

Gost12S512Status Gost12S512Verify( const char* publicKeyX,
                                   const char* publicKeyY,
                                   const char* hash,
                                   const char* signature ) {
    return s->verify(reinterpret_cast<const byte*>(publicKeyX),
                   reinterpret_cast<const byte*>(publicKeyY),
                   reinterpret_cast<const byte*>(hash),
                   reinterpret_cast<const byte*>(signature));
}
