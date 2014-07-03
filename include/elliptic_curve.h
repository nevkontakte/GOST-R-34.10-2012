#ifndef ELLIPTIC_CURVE_H
#define ELLIPTIC_CURVE_H

#include <prime_field.h>

#include <ostream>

namespace gost_ecc {

template <typename _integer_type, typename _double_integer_type>
class elliptic_curve {
public:
    using field_type = prime_field<_integer_type, _double_integer_type>;
    using integer_type = typename field_type::integer_type;

    struct point {
        integer_type x;
        integer_type y;

        bool is_infinity;

        point()
            :point(0, 0)
        {}

        point(integer_type x, integer_type y)
            :x(x), y(y), is_infinity(false)
        {}

        bool operator ==(const point& that) const {
            return (this->is_infinity == false && that.is_infinity == false && this->x == that.x && this->y == that.y) ||
                    (this->is_infinity == true && that.is_infinity == true);
        }

        static const point inf;

    private:
        point(bool infinity)
            :is_infinity(infinity)
        {}
    };

    const field_type field;
    const integer_type a;
    const integer_type b;

    elliptic_curve(integer_type modulus, integer_type a, integer_type b)
        :field(modulus), a(a), b(b)
    {}

    point add(const point& left, const point& right) const {
        if (left == point::inf) {
            return right;
        } else if (right == point::inf) {
            return left;
        } else if (left == right) {
            return this->twice(left);
        }

        const field_type& f = this->field;
        integer_type delta_y = f.sub(right.y, left.y);
        integer_type delta_x = f.sub(right.x, left.x);
        integer_type lambda = f.mul(delta_y, f.mul_inverse(delta_x));

        point result;
        result.x = f.sub(f.sub(f.mul(lambda, lambda), left.x), right.x);
        result.y = f.sub(f.mul(lambda, f.sub(left.x, result.x)), left.y);

        return result;
    }

    point twice(const point& p) const {
        if (p == point::inf) {
            return p;
        }

        const field_type& f = this->field;
        integer_type lambda = f.add(f.mul(3, f.mul(p.x, p.x)), this->a);
        lambda = f.mul(lambda, f.mul_inverse(f.mul(2, p.y)));

        point result;
        result.x = f.sub(f.mul(lambda, lambda), f.mul(2, p.x));
        result.y = f.sub(f.mul(lambda, f.sub(p.x, result.x)), p.y);

        return result;
    }

    point mulScalar(const point& p, const integer_type& multiplier) const {
        point result = point::inf;

        unsigned total_bits = multiplier.backend().size() * sizeof(mp::limb_type) * 8;

        for (unsigned i = total_bits; i > 0; i--) {
            result = this->twice(result);
            if (mp::bit_test(multiplier, i - 1)) {
                result = this->add(result, p);
            }
        }

        return result;
    }

    friend std::ostream& operator<<(std::ostream& out, point& p) {
        if (p == point::inf) {
            out << "(inf, inf)";
        } else {
            out << '(' << p.x << ", " << p.y << ')';
        }
        return out;
    }
};

template <typename _integer_type, typename _double_integer_type>
const typename elliptic_curve<_integer_type, _double_integer_type>::point elliptic_curve<_integer_type, _double_integer_type>::point::inf(true);
}

#endif // ELLIPTIC_CURVE_H
