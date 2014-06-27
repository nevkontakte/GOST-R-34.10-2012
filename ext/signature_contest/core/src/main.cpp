#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <string.h>

#include <signature_contest/core/src/load_dl.h>

using namespace std;
using namespace std::chrono;

#define GOST12S512_PRIVATE_KEY_SIZE 64
#define GOST12S512_PUBLIC_KEY_SIZE 64
#define GOST12S512_HASH_SIZE 64
#define GOST12S512_RAND_SIZE 64
#define GOST12S512_SIGN_SIZE 128

static const char* strstatus( Gost12S512Status status )
{
     switch( status )
     {
          case kStatusOk:
               return "kStatusOk";
          case kStatusWrongSignature:
               return "kStatusWrongSignature";
          case kStatusBadInput:
               return "kStatusBadInput";
          case kStatusInternalError:
               return "kStatusInternalError";
     }
     return 0;
}

static void print_sign_error( const char* context, Gost12S512Status status )
{
     cerr << "[" << context << "]: ";

     const char* str = strstatus( status );
     if( str == 0 )
     {
          cerr << "unknown error status (" << status << ")" << endl;
     }
     else
     {
          cerr << str << endl;
     }
}

typedef struct
{
     char privateKey[ GOST12S512_PRIVATE_KEY_SIZE ];
     char publicKeyX[ GOST12S512_PUBLIC_KEY_SIZE ];
     char publicKeyY[ GOST12S512_PUBLIC_KEY_SIZE ];
     char rand[ GOST12S512_RAND_SIZE ];
     char hash[ GOST12S512_HASH_SIZE ];
     char sign[ GOST12S512_SIGN_SIZE ];
} TestEntry;

#define DUMP_BYTES_AT_ROW 16

static void dump_array( const char* arr, size_t size )
{
     size_t i, j;
     size_t ni = (size + DUMP_BYTES_AT_ROW - 1) / DUMP_BYTES_AT_ROW;

     for( i = 0; i < ni; i++ )
     {
          size_t nj = (i < ni - 1) ? DUMP_BYTES_AT_ROW : (size - i * DUMP_BYTES_AT_ROW);

          cerr << "\t";
          for( j = 0; j < nj; j++ )
          {
               cerr << setfill('0') << setw(2) << hex
                    << +static_cast< unsigned char >( arr[ i * DUMP_BYTES_AT_ROW + j ] ) << " ";
          }
          cerr << "\n";
     }
}

static void dump_input_data( bool sign, const TestEntry& entry )
{
     cerr << "Dump of input data:" << endl << endl;
     if( sign )
     {
          cerr << "Private key:" << endl;
          dump_array( entry.privateKey, sizeof(entry.privateKey) );
          cerr << endl;

          cerr << "Random:" << endl;
          dump_array( entry.rand, sizeof(entry.rand) );
          cerr << endl;
     }
     else
     {
          cerr << "Public key (X):" << endl;
          dump_array( entry.publicKeyX, sizeof(entry.publicKeyX) );
          cerr << endl;

          cerr << "Public key (Y):" << endl;
          dump_array( entry.publicKeyY, sizeof(entry.publicKeyY) );
          cerr << endl;
     }

     cerr << "Hash:" << endl;
     dump_array( entry.hash, sizeof(entry.hash) );
     cerr << endl;

     if( !sign )
     {
          cerr << "Generated signature:" << endl;
     }
     else
     {
          cerr << "Expected signature:" << endl;
     }
     dump_array( entry.sign, sizeof(entry.sign) );
     cerr << endl;
}

static bool get_next_entry( ifstream& in_file, TestEntry& entry )
{
     in_file.read( entry.privateKey, sizeof(entry.privateKey) );
     if( !in_file )
     {
          return false;
     }

     in_file.read( entry.publicKeyX, sizeof(entry.publicKeyX) );
     if( !in_file )
     {
          return false;
     }

     in_file.read( entry.publicKeyY, sizeof(entry.publicKeyY) );
     if( !in_file )
     {
          return false;
     }

     in_file.read( entry.rand, sizeof(entry.rand) );
     if( !in_file )
     {
          return false;
     }

     in_file.read( entry.hash, sizeof(entry.hash) );
     if( !in_file )
     {
          return false;
     }

     in_file.read( entry.sign, sizeof(entry.sign) );
     if( !in_file )
     {
          return false;
     }

     return true;
}

static void usage( const char* exec_name )
{
     cerr << "Usage: " << exec_name << " <library> <input file>" << endl;
}

#define MEASURED_CALL( func, ret, dur, ... )                          \
     do                                                               \
     {                                                                \
          high_resolution_clock::time_point before, after;            \
          before = high_resolution_clock::now();                      \
          ret = func( __VA_ARGS__ );                                  \
          after = high_resolution_clock::now();                       \
          dur = duration_cast< microseconds >( after - before );      \
     } while( 0 )

int main( int argc, char* argv[] )
{
     const char* libpath = argv[1];
     const char* infile  = argv[2];

     if( argc != 3 )
     {
          usage( argv[0] );
          return 1;
     }

     ifstream in_file( infile );
     if( !in_file.is_open() )
     {
          cerr << "Can't open file: " << infile << endl;
          return 2;
     }

     Gost12S512Init_t Gost12S512Init;
     Gost12S512Sign_t Gost12S512Sign;
     Gost12S512Verify_t Gost12S512Verify;
     void* dl_handle;

     dl_handle = Gost12S512DlLoad( libpath, &Gost12S512Init, &Gost12S512Sign, &Gost12S512Verify );
     if( dl_handle == 0 )
     {
          cerr << "Incorrect library: " << libpath << endl;
          in_file.close();
          return 3;
     }

     Gost12S512Status status;

     if( Gost12S512Init )
     {
          status = Gost12S512Init();
          if( status != kStatusOk )
          {
               print_sign_error( INIT_FUNCTION_NAME, status );
               Gost12S512DlCleanup( dl_handle );
               in_file.close();
               return 4;
          }
     }

     char signature[ GOST12S512_SIGN_SIZE ];

     microseconds sign_duration( 0 );
     microseconds verify_true_duration( 0 );
     microseconds verify_false_duration( 0 );
     unsigned int test_count = 0;
     TestEntry entry;

     while( get_next_entry( in_file, entry ) )
     {
          microseconds duration;

          // Sign test
          MEASURED_CALL( Gost12S512Sign, status, duration,
                         entry.privateKey, entry.rand, entry.hash, signature );
          if( status != kStatusOk ||
              memcmp( signature, entry.sign, sizeof(entry.sign) ) != 0 )
          {
               if( status == kStatusOk )
               {
                    cerr << "Incorrect signature generated" << endl;
//                    memcpy( entry.sign, signature, sizeof(entry.sign) );
               }
               print_sign_error( SIGN_FUNCTION_NAME, status );
               cerr << endl;
               dump_input_data( true, entry );
               cerr << "Actual sugnature: " << std::endl;
               dump_array(signature, sizeof(entry.sign));

               Gost12S512DlCleanup( dl_handle );
               in_file.close();
               return 5;
          }
          sign_duration += duration;

          // Verify test (with correct data)
          MEASURED_CALL( Gost12S512Verify, status, duration,
                         entry.publicKeyX, entry.publicKeyY, entry.hash, entry.sign );
          if( status == kStatusWrongSignature )
          {
               cerr << "Correct signature was treated as incorrect" << endl;
          }
          if( status != kStatusOk )
          {
               print_sign_error( VERIFY_FUNCTION_NAME, status );
               cerr << endl;
               dump_input_data( false, entry );
               Gost12S512DlCleanup( dl_handle );
               in_file.close();
               return 6;
          }
          verify_true_duration += duration;

          // Verify test (with incorrect data)
          entry.hash[ test_count % sizeof(entry.hash) ] ^= 0xFF;

          MEASURED_CALL( Gost12S512Verify, status, duration,
                         entry.publicKeyX, entry.publicKeyY, entry.hash, entry.sign );
          if( status == kStatusOk )
          {
               cerr << "Incorrect signature was treated as correct" << endl;
          }
          if( status != kStatusWrongSignature )
          {
               print_sign_error( VERIFY_FUNCTION_NAME, status );
               cerr << endl;
               dump_input_data( false, entry );
               Gost12S512DlCleanup( dl_handle );
               in_file.close();
               return 7;
          }
          verify_false_duration += duration;

          test_count++;
     }

     Gost12S512DlCleanup( dl_handle );

     if( !in_file.eof() )
     {
          cerr << "Input file has incorrect format" << endl;
          in_file.close();
          return 8;
     }

     in_file.close();

     cout << "Average sign time: " << sign_duration.count() / test_count
          << " microseconds" << endl;
     cout << "Average correct signature verification time "
          << verify_true_duration.count() / test_count
          << " microseconds" << endl;
     cout << "Average incorrect signature verification time "
          << verify_false_duration.count() / test_count
          << " microseconds" << endl;

     return 0;
}

