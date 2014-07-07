#ifndef PRIME_FIELD_H
#define PRIME_FIELD_H

#include <cyclic_array.h>

#include <boost/multiprecision/cpp_int.hpp>
#include <array>

namespace gost_ecc {

namespace mp = ::boost::multiprecision;

template <typename _integer_type, typename _double_integer_type>
class prime_field {
public:
    typedef _integer_type integer_type;
    typedef _double_integer_type double_integer_type;

    static const std::size_t bits = integer_type::backend_type::internal_limb_count * integer_type::backend_type::limb_bits;

    const integer_type modulus;

    prime_field(integer_type modulus)
        :modulus(modulus)
    {}

    integer_type acquire(const integer_type& n) const {
        return n % this->modulus;
    }

    integer_type add(const integer_type& left, const integer_type& right) const {
        double_integer_type sum;
        mp::add(sum, left, right);

        return static_cast<integer_type>(sum % this->modulus);
    }

    integer_type sub(const integer_type& left, const integer_type& right) const {
        double_integer_type sum = left;

        if (left < right) {
            sum += this->modulus;
        }
        sum -= right;

        return static_cast<integer_type>(sum % this->modulus);
    }

    integer_type inverse(const integer_type& n) const {
        return (n == 0) ? n : (this->modulus - n);
    }

    integer_type mul(const integer_type& left, const integer_type& right) const {
        double_integer_type sum;
        mp::multiply(sum, left, right);

        return static_cast<integer_type>(sum % this->modulus);
    }

    integer_type mul_inverse(const integer_type& n) const {
        std::cout << n << " mod " << this->modulus << std::endl;

        cyclic_array<integer_type, 200> s{1, 0};
        cyclic_array<integer_type, 200> r{n, this->modulus};
        integer_type quotient, remainder;

        std::size_t i;

        for (i = 1; r[i] != 0; i++) {
            std::cout << r[i-1] << " / " <<  r[i] << std::endl;
            mp::divide_qr(r[i-1], r[i], quotient, remainder);
            r[i+1] = this->sub(r[i-1], this->mul(quotient, r[i]));
            s[i+1] = this->sub(s[i-1], this->mul(quotient, s[i]));;
        }

        if (r[i-1] != 1) {
            throw std::invalid_argument("Provided number isn't inversible by specified modulus.");
        }

        return s[i-1];
    }

    template<typename T>
    static integer_type import_bytes(const T* data) {
        const mp::limb_type* src = reinterpret_cast<const mp::limb_type*>(data);

        integer_type val;
        val.backend().resize(bits,bits);
        std::copy(src, src + val.backend().size(), val.backend().limbs());
        return val;
    }

    template<typename T>
    static T* export_bytes(integer_type val, T* data) {
        mp::limb_type* dst = reinterpret_cast<mp::limb_type*>(data);

        std::copy(val.backend().limbs(), val.backend().limbs() + val.backend().size(), dst);
        return data;
    }
};

template <unsigned bits>
using cpp_int_fixed = mp::number<mp::cpp_int_backend<bits, bits, mp::unsigned_magnitude, mp::unchecked, void> >;

}
#endif // PRIME_FIELD_H
