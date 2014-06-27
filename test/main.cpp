#include <curve.h>
#include <signature.h>

#include <iostream>

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

int main() {
    signature s(p, a, b, q, x, y);
    byte result[64*2];
    std::cout << ">> Signing << " << std::endl;
    s.sign(to_bytes(d), to_bytes(rnd), to_bytes(alpha), result);
}
