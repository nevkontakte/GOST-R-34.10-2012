#ifndef LOAD_DL_H
#define LOAD_DL_H

#include <signature_contest/sign_engine.h>

#define INIT_FUNCTION_NAME "Gost12S512Init"
#define SIGN_FUNCTION_NAME "Gost12S512Sign"
#define VERIFY_FUNCTION_NAME "Gost12S512Verify"

typedef Gost12S512Status (*Gost12S512Init_t)();
typedef Gost12S512Status (*Gost12S512Sign_t)( const char* privateKey,
                                              const char* rand,
                                              const char* hash,
                                              char* signature );
typedef Gost12S512Status (*Gost12S512Verify_t)( const char* publicKeyX,
                                                const char* publicKeyY,
                                                const char* hash,
                                                const char* signature );

#ifdef __cplusplus
extern "C" {
#endif //__cplusplus

void* Gost12S512DlLoad( const char* libpath,
                        Gost12S512Init_t* init,
                        Gost12S512Sign_t* sign,
                        Gost12S512Verify_t* verify );

void Gost12S512DlCleanup( void* handle );

#ifdef __cplusplus
}
#endif //__cplusplus

#endif /* LOAD_DL_H */

