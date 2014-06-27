#include <iostream>
#include <sign_engine.h>

Gost12S512Status Gost12S512Sign( const char* privateKey,
                                 const char* rand,
                                 const char* hash,
                                 char* signature ) {
    return kStatusInternalError;
}

Gost12S512Status Gost12S512Verify( const char* publicKeyX,
                                   const char* publicKeyY,
                                   const char* hash,
                                   const char* signature ) {
    return kStatusInternalError;
}
