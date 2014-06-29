#ifndef PRIME_FIELD_H
#define PRIME_FIELD_H

#include <boost/multiprecision/cpp_int.hpp>

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

    integer_type add(const integer_type& left, const integer_type& right) {
        double_integer_type sum;
        mp::add(sum, left, right);

        return static_cast<integer_type>(sum % this->modulus);
    }

    integer_type sub(const integer_type& left, const integer_type& right) {
        double_integer_type sum = left;

        if (left < right) {
            sum += this->modulus;
        }
        sum -= right;

        return static_cast<integer_type>(sum % this->modulus);
    }

    integer_type mul(const integer_type& left, const integer_type& right) {
        double_integer_type sum;
        mp::multiply(sum, left, right);

        return static_cast<integer_type>(sum % this->modulus);
    }

    integer_type mul_inverse(const integer_type& n) {
        return mp::powm(n, this->modulus - 2, this->modulus);
    }

    static integer_type import_bytes(const byte* data) {
        const mp::limb_type* src = reinterpret_cast<const mp::limb_type*>(data);

        integer_type val;
        val.backend().resize(bits,bits);
        std::copy(src, src + val.backend().size(), val.backend().limbs());
        return val;
    }
};

}
#endif // PRIME_FIELD_H
