#include <signature.h>

#include <iostream>

using namespace CryptoPP;

namespace gost_ecc {

static const std::size_t private_key_size = 64;
static const std::size_t rand_size = 64;
static const std::size_t hash_size = 64;
static const std::size_t signature_size = 64 * 2;

// Convert Little-Endian representation into CryptoPP Integer object
template<std::size_t size>
inline Integer import_integer(byte* data) {
    std::reverse(data, data + size);
    return Integer(data, size);
}

template<std::size_t size>
inline Integer import_integer(const byte* data) {
    byte bytes[size];
    std::copy(data, data + size, &bytes[0]);
    return import_integer<size>(bytes);
}

inline Integer import_integer(u_int64_t (&data)[8]) {
    byte* bytes = reinterpret_cast<byte*>(&data);
    return import_integer<sizeof(data)>(bytes);
}

template<std::size_t size>
inline void export_integer(const Integer& n, byte* output) {
    n.Encode(output, size);
    std::reverse(output, output + size);
}

signature::signature(u_int64_t (&modulus)[8], u_int64_t (&a)[8], u_int64_t (&b)[8],
                     u_int64_t (&subgroupModulus)[8],
                     u_int64_t (&base_x)[8], u_int64_t (&base_y)[8])
    :curve(import_integer(modulus), import_integer(a), import_integer(b)),
      subgroup(import_integer(subgroupModulus)),
      basePoint(import_integer(base_x), import_integer(base_y))
{
#ifdef DEBUG
    std::cout << this->curve.GetField().GetModulus() << std::endl << std::endl
                 << this->curve.GetA() << std::endl << this->curve.GetB() << std::endl << std::endl
                    << this->subgroup.GetModulus() << std::endl << std::endl
                       << this->basePoint.x << std::endl << this->basePoint.y << std::endl << std::endl;
#endif
}

Gost12S512Status signature::sign(const byte* private_key, const byte* rand, const byte* hash, byte* signature) {
    Integer alpha = import_integer<hash_size>(hash);
    Integer d = import_integer<private_key_size>(private_key);

#ifdef DEBUG
    std::cout << "alpha: " << alpha << std::endl;
    std::cout << "d: " << d << std::endl;
#endif

    Integer e = this->subgroup.ConvertIn(alpha);
    if (e == 0) {
        e = 1;
    }
#ifdef DEBUG
    std::cout << "e: " << e << std::endl;
#endif

    Integer k = import_integer<rand_size>(rand);
    if(k >= this->subgroup.GetModulus()) {
        return kStatusBadInput;
    }

#ifdef DEBUG
    std::cout << "k: " << k << std::endl;
#endif

    ECPPoint C = this->curve.ScalarMultiply(this->basePoint, k);

#ifdef DEBUG
    std::cout << "x_c: " << C.x << std::endl << "y_c: " << C.y << std::endl;
#endif

    Integer r = this->subgroup.ConvertIn(C.x);

#ifdef DEBUG
    std::cout << "r: " << r << std::endl;
#endif

    if (r == 0) {
        return kStatusBadInput;
    }

    Integer rd = this->subgroup.Multiply(r, d);
    Integer ke = this->subgroup.Multiply(k, e);
    Integer s = this->subgroup.Add(rd, ke);

#ifdef DEBUG
    std::cout << "rd: " << rd << std::endl;
    std::cout << "ke: " << ke << std::endl;

    std::cout << "s:  " << s << std::endl;
#endif

    if (s == 0) {
        return kStatusBadInput;
    }

    export_integer<signature_size / 2>(r, signature);
    export_integer<signature_size / 2>(s, signature + signature_size / 2);

#ifdef DEBUG
    std::cout << "Done!" << std::endl;
#endif

    return kStatusOk;
}

}
