#include <stdlib.h>
#include <dlfcn.h>

#include <signature_contest/core/src/load_dl.h>

void* Gost12S512DlLoad( const char* libpath,
                        Gost12S512Init_t* init,
                        Gost12S512Sign_t* sign,
                        Gost12S512Verify_t* verify )
{
     void* handle;

     if( libpath == NULL || init == NULL || sign == NULL || verify == NULL )
     {
          return NULL;
     }

     handle = dlopen( libpath, RTLD_NOW );
     if( handle == NULL )
     {
          return NULL;
     }

     dlerror();     /* Clear any existing error */

     *(void**) init = dlsym( handle, INIT_FUNCTION_NAME );
     dlerror();     /* Init is optional */

     *(void**) sign = dlsym( handle, SIGN_FUNCTION_NAME );
     if( dlerror() != NULL )
     {
          dlclose( handle );
          return NULL;
     }

     *(void**) verify = dlsym( handle, VERIFY_FUNCTION_NAME );
     if( dlerror() != NULL )
     {
          dlclose( handle );
          return NULL;
     }

     return handle;
}

void Gost12S512DlCleanup( void* handle )
{
     if( handle )
     {
          dlclose( handle );
     }
}

