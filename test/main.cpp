#include <curve.h>
#include <signature.h>
#include <prime_field.h>

#include <iostream>

#define ASSERT_TRUE(expr) \
    if(!(expr)) {throw std::logic_error("Assertion failed in " + std::string(__FILE__) + " at line " + std::to_string(__LINE__));}

using namespace std;
using namespace gost_ecc;

const uint64_t d[8] =
{0x1D19CE9891EC3B28,
 0x1B60961F49397EEE,
 0x10ED359DD39A72C1,
 0x7A929ADE789BB9BE,
 0x0000000000000000,
 0x0000000000000000,
 0x0000000000000000,
 0x0000000000000000};

const uint64_t alpha[8] =
{0x67ECE6672B043EE5,
 0xCE52032AB1022E8E,
 0x88C09C52E0EEC61F,
 0x2DFBC1B372D89A11,
 0x0000000000000000,
 0x0000000000000000,
 0x0000000000000000,
 0x0000000000000000};

const uint64_t rnd[8] =
{0x4FED924594DCEAB3,
 0x6DE33814E95B7FE6,
 0x2823C8CF6FCC7B95,
 0x77105C9B20BCD312,
 0x0000000000000000,
 0x0000000000000000,
 0x0000000000000000,
 0x0000000000000000};

inline const byte* to_bytes(const uint64_t (&data)[8]) {
    return reinterpret_cast<const byte*>(&data);
}

using namespace mp;

int main() {

    {
        typedef prime_field<mp::uint256_t, mp::uint512_t> pf;

        mp::uint256_t left = pf::import_bytes(to_bytes(alpha));
        mp::uint256_t right = pf::import_bytes(to_bytes(rnd));

        mp::uint256_t modulus = pf::import_bytes(to_bytes(p));

        std::cout << std::hex << left << std::endl << right << std::endl << modulus << std::endl;

        pf field(modulus);

        const mp::uint256_t sum("0x250c1e4e93956d23b0e4652250bb41b53c353b3f9a5dae74b7da78acbfe12567");
        ASSERT_TRUE(sum == field.add(left, right));

        const mp::uint256_t product("0x42a551bb2d40345366cbaae3a4d42b269ec7f48355533a13c0ad3953a31a55e8");
        ASSERT_TRUE(product == field.mul(left, right));

        const mp::uint256_t diff("0x36eb6518521bc6ff609cd38371224a8a606ecb15c7a6aea817ff542196275863");
        ASSERT_TRUE(diff == field.sub(left, right));

        const mp::uint256_t inv = field.mul_inverse(left);
        ASSERT_TRUE(1 == field.mul(left, inv));
    }

    signature s(p, a, b, q, x, y);
    byte result[64*2];
    std::cout << ">> Signing << " << std::endl;
    s.sign(to_bytes(d), to_bytes(rnd), to_bytes(alpha), result);
}
