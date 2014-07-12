#include <signature.h>

#include <iostream>

namespace gost_ecc {

const unsigned signature_size = 64 * 2;

signature::signature(u_int64_t (&modulus)[8], u_int64_t (&a)[8], u_int64_t (&b)[8],
                     u_int64_t (&subgroupModulus)[8],
                     u_int64_t (&base_x)[8], u_int64_t (&base_y)[8])
    :curve(pf::import_bytes(modulus), pf::import_bytes(a), pf::import_bytes(b)),
      subgroup(pf::import_bytes(subgroupModulus)),
      basePoint(pf::import_bytes(base_x), pf::import_bytes(base_y))
{
#ifdef DEBUG
    std::cout << "p: " << this->curve.field.modulus << std::endl
                 << "a: " << this->curve.a << std::endl << "b: " << this->curve.b << std::endl
                    << "m: " << this->subgroup.modulus << std::endl
                       << "x_p: " << this->basePoint.x << std::endl << "y_p: " << this->basePoint.y << std::endl << std::endl;
#endif
    this->curve.comb_precompute(this->basePoint, this->basePointTable);
}

Gost12S512Status signature::sign(const byte* private_key, const byte* rand, const byte* hash, byte* signature) {
    pf::integer_type alpha = pf::import_bytes(hash);
    pf::integer_type d = pf::import_bytes(private_key);

#ifdef DEBUG
    std::cout << "alpha: " << alpha << std::endl;
    std::cout << "d: " << d << std::endl;
#endif

    pf::integer_type e = this->subgroup.acquire(alpha);
    if (e == 0) {
        e = 1;
    }
#ifdef DEBUG
    std::cout << "e: " << e << std::endl;
#endif

    pf::integer_type k = pf::import_bytes(rand);
    if(k >= this->subgroup.modulus) {
        return kStatusBadInput;
    }

#ifdef DEBUG
    std::cout << "k: " << k << std::endl;
#endif

    ec::point C = this->curve.mul_scalar(this->basePointTable, k);

#ifdef DEBUG
    std::cout << "x_c: " << C.x << std::endl << "y_c: " << C.y << std::endl;
#endif

    pf::integer_type r = this->subgroup.acquire(C.x);

#ifdef DEBUG
    std::cout << "r: " << r << std::endl;
#endif

    if (r == 0) {
        return kStatusBadInput;
    }

    pf::integer_type rd = this->subgroup.mul(r, d);
    pf::integer_type ke = this->subgroup.mul(k, e);
    pf::integer_type s = this->subgroup.add(rd, ke);

#ifdef DEBUG
    std::cout << "rd: " << rd << std::endl;
    std::cout << "ke: " << ke << std::endl;

    std::cout << "s:  " << s << std::endl;
#endif

    if (s == 0) {
        return kStatusBadInput;
    }

    pf::export_bytes(r, signature);
    pf::export_bytes(s, signature + signature_size / 2);

#ifdef DEBUG
    std::cout << std::hex << r << std::endl << s << std::endl;

    std::cout << "Done!" << std::endl;
#endif

    return kStatusOk;
}

Gost12S512Status signature::verify(const byte* public_key_x, const byte* public_key_y, const byte* hash, const byte* signature) {
    pf::integer_type r = pf::import_bytes(signature);
    pf::integer_type s = pf::import_bytes(signature + signature_size / 2);

    pf::integer_type alpha = pf::import_bytes(hash);

    pf::integer_type e = this->subgroup.acquire(alpha);
    if (e == 0) {
        e = 1;
    }

    pf::integer_type v = this->subgroup.mul_inverse(e);

    pf::integer_type z_1 = this->subgroup.mul(s, v);
    pf::integer_type z_2 = this->subgroup.mul(r, v);
    z_2 = this->subgroup.inverse(z_2);

    ec::point Q(pf::import_bytes(public_key_x), pf::import_bytes(public_key_y));

    ec::point C = this->curve.add(this->curve.mul_scalar(this->basePointTable, z_1), this->curve.mul_scalar(Q, z_2));

    pf::integer_type R = this->subgroup.acquire(C.x);

    if (R == r) {
        return kStatusOk;
    } else {
        return kStatusWrongSignature;
    }
}

}
