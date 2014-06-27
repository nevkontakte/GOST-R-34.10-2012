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
}

Gost12S512Status signature::sign(const byte* private_key, const byte* rand, const byte* hash, byte* signature) {
    Integer alpha = import_integer<hash_size>(hash);
    Integer d = import_integer<private_key_size>(private_key);

    Integer e = this->subgroup.ConvertIn(alpha);
    if (e == 0) {
        e = 1;
    }

    Integer k = import_integer<rand_size>(rand);
    if(k >= this->subgroup.GetModulus()) {
        return kStatusBadInput;
    }

    ECPPoint C = this->curve.ScalarMultiply(this->basePoint, k);

    Integer r = this->subgroup.ConvertIn(C.x);

    if (r == 0) {
        return kStatusBadInput;
    }

    Integer s = this->subgroup.Add(
                this->subgroup.Multiply(r, d),
                this->subgroup.Multiply(k, e)
                );

    if (s == 0) {
        return kStatusBadInput;
    }

    export_integer<signature_size / 2>(r, signature);
    export_integer<signature_size / 2>(r, signature + signature_size / 2);

    std::cout << "Done!" << std::endl;

    return kStatusOk;
}

}
