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
//        std::cout << n << " mod " << this->modulus << std::endl;

        integer_type s0 = 1, s1 = 0;
        integer_type r0 = n, r1 = this->modulus;
        integer_type quotient, remainder;

        std::size_t i;

        for (i = 1; r1 != 0; i++) {
//            std::cout << r0 << " / " <<  r1 << std::endl;
            mp::divide_qr(r0, r1, quotient, remainder);
            r0 = this->sub(r0, this->mul(quotient, r1));
            s0 = this->sub(s0, this->mul(quotient, s1));;

            std::swap(r0, r1);
            std::swap(s0, s1);
        }

        if (r0 != 1) {
            throw std::invalid_argument("Provided number isn't inversible by specified modulus.");
        }

        return s0;
    }

    template<typename T>
    static integer_type import_bytes(const T* data) {
        const mp::limb_type* src = reinterpret_cast<const mp::limb_type*>(data);

        integer_type val = 0;
        for (unsigned i = 0; i < bits / (sizeof(mp::limb_type) * 8); i++) {
            val += integer_type(src[i]) << (i * sizeof(mp::limb_type) * 8);
        }

        return val;
    }

    template<typename T>
    static T* export_bytes(integer_type val, T* data) {
        mp::limb_type* dst = reinterpret_cast<mp::limb_type*>(data);

        mp::limb_type* limbs = val.backend().limbs();
        unsigned limb_number = val.backend().size();
        unsigned limb_limit = bits / (sizeof(mp::limb_type) * 8);

        std::copy(limbs, limbs + limb_number, dst);
        std::fill_n(dst + limb_number, limb_limit - limb_number, 0);
        return data;
    }
};

template <unsigned bits>
using cpp_int_fixed = mp::number<mp::cpp_int_backend<bits, bits, mp::unsigned_magnitude, mp::unchecked, void> >;

}
#endif // PRIME_FIELD_H
