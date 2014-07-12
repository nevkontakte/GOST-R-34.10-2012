#include <curve.h>
#include <signature.h>
#include <prime_field.h>
#include <elliptic_curve.h>

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

const uint64_t x_q[8] =
{0x6689dbd8e56fd80b,
 0x8585ba1d4e9b788f,
 0xd8595bec458b50c5,
 0x7f2b49e270db6d90};

const uint64_t y_q[8] =
{0xdffb101a87ff77da,
 0xaf64d1c593d26627,
 0x85c8413a977b3cbb,
 0x26f1b489d6701dd1};


template <unsigned n>
inline const byte* to_bytes(const uint64_t (&data)[n]) {
    return reinterpret_cast<const byte*>(&data);
}

template <unsigned n>
inline byte* to_bytes(uint64_t (&data)[n]) {
    return reinterpret_cast<byte*>(&data);
}

using namespace mp;

int main() {

    {
        typedef prime_field<mp::uint256_t, mp::uint512_t> pf;

        mp::uint256_t left = pf::import_bytes(to_bytes(alpha));
        mp::uint256_t right = pf::import_bytes(to_bytes(rnd));

        mp::uint256_t modulus = pf::import_bytes(to_bytes(p));

        ASSERT_TRUE(modulus == mp::uint256_t("0x8000000000000000000000000000000000000000000000000000000000000431"));

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

    {
        typedef prime_field<mp::uint512_t, mp::uint1024_t> pf;

        pf::integer_type modulus ("0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFDC7");
        pf field(modulus); // 2^512 - 569
        const pf::double_integer_type n ("0x67D4C9B41F674164FA6D707CC77486E71A2ED5149776F56514827F17A7990E30FAD580DFDFF4BE1E81975E3057D97AB43AC359645A234581805C4FE6FE51ED1BA241600EE7716F6C73D7EE518FAE0FEB0F89E9FBA5F9B3FD6781F6B5A006C5E7342DF0D67A31279EC44C26AE4D0D537F7AA592E39DEB5DF61E296EF605FA8BA9");
        const pf::integer_type expected ("0x6A35B168B3F1C8DD1116F3A8E1ADE79441A184C04D6121A7FD8E7249233750C4B8B9626F412BB96CCDBC80218F6F0E1616D742EFF650DACC6B5707614A148E1B");
        ASSERT_TRUE(expected == field.reduce(n))
    }

    {
        typedef elliptic_curve<mp::uint256_t, mp::uint512_t> ec;

        ec curve(17, 17-3, 2);

        ec::point left(16, 13);
        ec::point right(3, 1);

        ec::point result = curve.add(left, right);
        ASSERT_TRUE(ec::point(7, 11) == result);
        ASSERT_TRUE(ec::point(7, 11) == curve.add(ec::jacobian_point(left), right).to_affine(curve));

        ASSERT_TRUE(ec::point(2, 4) == curve.twice(left)); // 2 4
        ASSERT_TRUE(ec::point(2, 4) == curve.twice(ec::jacobian_point(left)).to_affine(curve));
        ASSERT_TRUE(ec::point(2, 4) == curve.repeated_twice(ec::jacobian_point(left), 1).to_affine(curve));
        ASSERT_TRUE(curve.twice(ec::point(2, 4)) == curve.repeated_twice(ec::jacobian_point(left), 2).to_affine(curve));

        ASSERT_TRUE(ec::point::inf == curve.add(ec::point::inf, ec::point::inf));
        ASSERT_TRUE(right == curve.add(ec::point::inf, right));
        ASSERT_TRUE(left == curve.add(left, ec::point::inf));
        ASSERT_TRUE(ec::point::inf == curve.twice(ec::point::inf));

        ASSERT_TRUE(ec::point(14, 16) == curve.mul_scalar(ec::point(0, 6), 6));

        ec::point table[256];
        curve.comb_precompute(left, table);
        ASSERT_TRUE(table[0] == ec::point::inf);
        ASSERT_TRUE(table[1] == left);
        ASSERT_TRUE(table[128] == curve.repeated_twice(left, ec::field_type::bits - ec::field_type::bits/8).to_affine(curve));

        const ec::integer_type multiplier("0x2DFBC1B372D89A1188C09C52E0EEC61FCE52032AB1022E8E67ECE6672B043EE5");
        const ec::point expected = curve.mul_scalar(left, multiplier);
        ASSERT_TRUE(expected == curve.mul_scalar(table, multiplier));
    }

    {
        typedef prime_field<mp::uint512_t, mp::uint1024_t> pf;

        uint64_t m[8] = {
            0x0000000000000431,
            0x0000000000000000,
            0x0000000000000000,
            0x8000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
            0x0000000000000000,
        };
        pf field(pf::import_bytes(m));
        auto result = field.mul_inverse(pf::integer_type("8037948113078075006670898845874119551271478779811090161381958730426863132560"));
        ASSERT_TRUE(pf::integer_type("10297018950366695783893339349991847804652198502529290759859202820577721280788") == result);
    }

    std::cout << "General test passed, testing signature..." << std::endl;

    signature s(p, a, b, q, x, y);
    uint64_t result[8 * 2];
    s.sign(to_bytes(d), to_bytes(rnd), to_bytes(alpha), to_bytes(result));

    uint64_t expected[8 * 2] = {
        0x3ad043fd39dc0493,
        0x74053554a42767b8,
        0x80cd9ed56feda419,
        0x41aa28d2f1ab1482,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x928014f6c5bf9c40,
        0xbcd6d3f746b631df,
        0x653c235a98a60249,
        0x01456c64ba4642a1,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
        0x0000000000000000,
    };

    ASSERT_TRUE(std::equal(std::begin(expected), std::end(expected), std::begin(result)));

    ASSERT_TRUE(s.verify(to_bytes(x_q), to_bytes(y_q), to_bytes(alpha), to_bytes(result)) == kStatusOk);

    std::cout << "All tests passed!" << std::endl;
}
