#include <signature_contest/core/src/windows/targetver.h>

#include <windows.h>

#include <signature_contest/core/src/load_dl.h>

void* Gost12S512DlLoad( const char* libpath,
                        Gost12S512Init_t* init,
                        Gost12S512Sign_t* sign,
                        Gost12S512Verify_t* verify )
{
     HMODULE handle;

     if( libpath == NULL || init == NULL || sign == NULL || verify == NULL )
     {
          return NULL;
     }

     handle = LoadLibraryA( libpath );
     if( NULL == handle )
     {
          return NULL;
     }

     *(void**) init = GetProcAddress( handle, INIT_FUNCTION_NAME );

     *(void**) sign = GetProcAddress( handle, SIGN_FUNCTION_NAME );
     if ( NULL == *sign )
     {
          FreeLibrary( handle );
          return NULL;
     }

     *(void**) verify = GetProcAddress( handle, VERIFY_FUNCTION_NAME );
     if( NULL == *verify )
     {
          FreeLibrary( handle );
          return NULL;
     }

     return handle;
}

void Gost12S512DlCleanup( void* handle )
{
     if( handle )
     {
          FreeLibrary( (HMODULE)handle );
     }
}

