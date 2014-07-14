#ifndef ELLIPTIC_CURVE_H
#define ELLIPTIC_CURVE_H

#include <prime_field.h>
#include <naf.h>

#include <ostream>
#include <cstdlib>

namespace gost_ecc {

template <typename _integer_type, typename _double_integer_type, typename _pm_integer_type = _double_integer_type>
class elliptic_curve {
public:
    using field_type = prime_field<_integer_type, _double_integer_type, _pm_integer_type>;
    using integer_type = typename field_type::integer_type;
    using double_integer_type = typename field_type::double_integer_type;

    struct point {
        integer_type x;
        integer_type y;

        point()
            :point(0, 0)
        {}

        point(integer_type x, integer_type y)
            :x(x), y(y)
        {}

        bool operator ==(const point& that) const {
            return (this->x == that.x && this->y == that.y);
        }

        static const point inf;
    };

    struct jacobian_point {
        integer_type x;
        integer_type y;
        integer_type z;

        jacobian_point()
            :x(0),y(0),z(0)
        {}

        jacobian_point(const integer_type& x, const integer_type& y, const integer_type& z = 1)
            :x(x),y(y),z(z)
        {}

        explicit jacobian_point(const point& affine)
        {
            if (affine == point::inf) {
                *this = jacobian_point::inf;
            } else {
                x = affine.x;
                y = affine.y;
                z = 1;
            }
        }

        point to_affine(const elliptic_curve& curve) const {
            if (*this == jacobian_point::inf) {
                return point::inf;
            }
            const field_type& f = curve.field;
            integer_type inv_z = f.mul_inverse(this->z); // z^-1
            integer_type inv_zz = f.mul(inv_z, inv_z); // z^-2
            integer_type inv_zzz = f.mul(inv_zz, inv_z); // z^-3

            return point(f.mul(this->x, inv_zz), f.mul(this->y, inv_zzz));
        }

        bool operator ==(const jacobian_point& that) const {
            return (this->x == that.x && this->y == that.y && this->z == that.z);
        }

        const static jacobian_point inf;
    };

    const field_type field;
    const integer_type a;
    const integer_type b;

protected:

    const bool allow_jacobian;

    const integer_type inv_2;

public:

    elliptic_curve(integer_type modulus, integer_type a, integer_type b)
        :field(modulus), a(a), b(b), allow_jacobian(this->field.inverse(a) == 3),
          inv_2(field.mul_inverse(2))
    {
    }

    point negate(const point& p) const {
        return point(p.x, this->field.inverse(p.y));
    }

    jacobian_point negate(const jacobian_point& p) const {
        return jacobian_point(p.x, this->field.inverse(p.y), p.z);
    }

    point add(const point& left, const point& right) const {
        if (left == point::inf) {
            return right;
        } else if (right == point::inf) {
            return left;
        } else if (left == right) {
            return this->twice(left);
        } else if (left == this->negate(right)) {
            return point::inf;
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

    /**
     * @brief Point addition for mixed Jacobian-affine coordinates.
     *
     * 8M + 3S
     * See: Hankerson, D., Vanstone, S., & Menezes, A. (2004). Guide to elliptic curve cryptography.
     * Page 91, alg. 3.22.
     * @param left
     * @param right
     * @return
     */
    jacobian_point add(const jacobian_point& left, const point& right) const {
        if (!this->allow_jacobian) {
            throw std::invalid_argument("Parameter a for curve must be -3");
        }

        if (left == jacobian_point::inf) {
            return jacobian_point(right);
        } else if (right == point::inf) {
            return left;
        }

        const field_type& f = this->field;

        jacobian_point result;
        integer_type t1, t2, t3, t4;

        t1 = f.mul(left.z, left.z);
        t2 = f.mul(t1, left.z);
        t1 = f.mul(t1, right.x);
        t2 = f.mul(t2, right.y);
        t1 = f.sub(t1, left.x);
        t2 = f.sub(t2, left.y);

        if (t1 == 0) {
            if (t2 == 0) {
                return this->twice(jacobian_point(right));
            } else {
                return jacobian_point::inf;
            }
        }

        result.z    = f.mul(left.z, t1);
        t3          = f.mul(t1, t1);
        t4          = f.mul(t3, t1);
        t3          = f.mul(t3, left.x);
        t1          = f.mul2(t3);
        result.x    = f.mul(t2, t2);
        result.x    = f.sub(result.x, t1);
        result.x    = f.sub(result.x, t4);
        t3          = f.sub(t3, result.x);
        t3          = f.mul(t3, t2);
        t4          = f.mul(t4, left.y);
        result.y    = f.sub(t3, t4);

        return result;
    }

    /**
     * @brief See: http://en.wikibooks.org/wiki/Cryptography/Prime_Curve/Jacobian_Coordinates
     *
     * 12M + 4S
     * @param left
     * @param right
     * @return
     */
    jacobian_point add(const jacobian_point &left, const jacobian_point &right) const {
        if (left == jacobian_point::inf) {
            return right;
        } else if (right == jacobian_point::inf) {
            return left;
        }

        const field_type& f = this->field;

        jacobian_point result;

        integer_type u1, u2, s1, s2, h, r, left_z_squared, right_z_squared, h_squared;

        left_z_squared = f.mul(left.z, left.z);
        right_z_squared = f.mul(right.z, right.z);

        u1 = f.mul(left.x, right_z_squared);
        u2 = f.mul(right.x, left_z_squared);

        left_z_squared = f.mul(left_z_squared, left.z);
        right_z_squared = f.mul(right_z_squared, right.z);

        s1 = f.mul(left.y, right_z_squared);
        s2 = f.mul(right.y, left_z_squared);

        if (u1 == u2) {
            if (s1 != s2) {
                return jacobian_point::inf;
            } else {
                return this->twice(left);
            }
        }

        h = f.sub(u2, u1);
        r = f.sub(s2, s1);

        result.z = f.mul(h, left.z, right.z);
        h_squared = f.mul(h, h);

        result.y = f.mul(u1, h_squared);
        result.x = f.mul2(result.y);

        h_squared = f.mul(h, h_squared);

        result.x = f.sub(f.mul(r, r), h_squared, result.x);
        result.y = f.sub(result.y, result.x);
        result.y = f.mul(r, result.y);
        result.y = f.sub(result.y, f.mul(s1, h_squared));

        return result;
    }

    template<typename P1, typename P2>
    P1 sub(const P1& left, const P2& right) const {
        return this->add(left, this->negate(right));
    }

    point twice(const point& p) const {
        if (p == point::inf) {
            return p;
        }

        const field_type& f = this->field;
        integer_type lambda = f.add(f.mul3(f.mul(p.x, p.x)), this->a);
        lambda = f.mul(lambda, f.mul_inverse(f.mul2(p.y)));

        point result;
        result.x = f.sub(f.mul(lambda, lambda), f.mul2(p.x));
        result.y = f.sub(f.mul(lambda, f.sub(p.x, result.x)), p.y);

        return result;
    }

    /**
     * @brief Point doubling for Jacobian coordinates.
     *
     * See: Hankerson, D., Vanstone, S., & Menezes, A. (2004). Guide to elliptic curve cryptography.
     * Page 91, alg. 3.21.
     * @param p
     * @return
     */
    jacobian_point twice(const jacobian_point& p) const {
        if (!this->allow_jacobian) {
            throw std::invalid_argument("Parameter a for curve must be -3");
        }

        if (p == jacobian_point::inf) {
            return p;
        }

        const field_type& f = this->field;
        jacobian_point result;
        integer_type t1, t2, t3;
        t1          = f.mul(p.z, p.z);
        t2          = f.sub(p.x, t1);
        t1          = f.add(p.x, t1);
        t2          = f.mul(t1, t2);
        t2          = f.mul3(t2); // A = 3 (X1 - Z1^2) (X1 + Z1^2)
        result.y    = f.mul2(p.y); // B = 2Y1
        result.z    = f.mul(result.y, p.z);
        result.y    = f.mul(result.y, result.y); // C = B^2
        t3          = f.mul(result.y, p.x); // D = CX1
        result.y    = f.mul(result.y, result.y);
        result.y    = f.mul(result.y, this->inv_2);
        result.x    = f.mul(t2, t2);
        t1          = f.mul2(t3);
        result.x    = f.sub(result.x, t1);
        t1          = f.sub(t3, result.x);
        t1          = f.mul(t1, t2);
        result.y    = f.sub(t1, result.y);

        return result;
    }


    /**
     * @brief Repeated doubling algorithm.
     *
     * See: Hankerson, D., Vanstone, S., & Menezes, A. (2004). Guide to elliptic curve cryptography.
     * Page 93, alg. 3.23.
     * @param p
     * @param count
     * @return
     */
    jacobian_point repeated_twice(const jacobian_point& p, unsigned count) const {
        if (p == jacobian_point::inf) {
            return jacobian_point::inf;
        }

        jacobian_point result = p;
        const field_type& f = this->field;

        integer_type a, b, w, y_squared, t1, t2;

        result.y    = f.mul2(result.y); // Y <- 2Y
        w           = f.mul(result.z, result.z);
        w           = f.mul(w, w); // W <- Z^4

        while (count > 0) {
            a           = f.mul(result.x, result.x); // a = X^2
            a           = f.sub(a, w); // a = X^2 - W
            a           = f.mul3(a); // a = 3 (X^2 - W)

            y_squared   = f.mul(result.y, result.y);
            b           = f.mul(result.x, y_squared); // B = X Y^2
            result.x    = f.sub(f.mul(a, a), f.mul2(b)); // X = A^2 - 2B
            result.z    = f.mul(result.z, result.y); // Z = ZY

            count--;

            y_squared   = f.mul(y_squared, y_squared); // y_squared = Y^4

            if (count > 0) {
                w = f.mul(w, y_squared); // W = W Y^4
            }

            result.y = f.sub(b, result.x); // B - X
            result.y = f.mul(a, result.y); // A(B - X)
            result.y = f.mul2(result.y); // 2A (B - X)
            result.y = f.sub(result.y, y_squared); // Y = 2A (B - X) - Y^4
        }

        result.y = f.mul(result.y, this->inv_2);

        return result;
    }

    point mul_scalar(const point& p, const integer_type& multiplier) const {
        jacobian_point result = jacobian_point::inf;

        unsigned total_bits = multiplier.backend().size() * sizeof(mp::limb_type) * 8;

        for (unsigned i = total_bits; i > 0; i--) {
            result = this->twice(result);

            if (mp::bit_test(multiplier, i - 1)) {
                result = this->add(result, p);
            }

        }

        return result.to_affine(*this);
    }

    template<unsigned win_left = 8>
    void comb_precompute(const point& base, point (&table)[1 << win_left]) const {
        const unsigned d = field_type::bits / win_left + ((field_type::bits % win_left) ? 1 : 0);

        point pow2[win_left];
        pow2[0] = base;

        for (unsigned i = 1; i < win_left; i++) {
            pow2[i] = this->repeated_twice(jacobian_point(pow2[i-1]), d).to_affine(*this);
        }

        for (unsigned i = 0; i < (1 << win_left); i++) {
            jacobian_point p = jacobian_point::inf;

            for (unsigned offset = 0; offset < win_left; offset++) {
                if ((i & (1 << offset)) != 0) {
                    p = this->add(p, pow2[offset]);
                }
            }

            table[i] = p.to_affine(*this);
        }
    }

    template<unsigned win_left = 8>
    point mul_scalar(const point (&comb_table)[1 << win_left], integer_type multiplier) const {
        const unsigned d = field_type::bits / win_left + ((field_type::bits % win_left) ? 1 : 0);

        integer_type mask = static_cast<integer_type>((double_integer_type(1) << d) - 1);

        integer_type chunks[win_left];
        for (unsigned i = 0; i < win_left; i++) {
            chunks[i] = multiplier & mask;
            multiplier >>= d;
        }

        jacobian_point result = jacobian_point::inf;

        for (unsigned i = d; i > 0; i--) {
            result = this->twice(result);

            unsigned key = 0;
            for (unsigned j = 0; j < win_left; j++) {
                key += mp::bit_test(chunks[j], i-1) << j;
            }

            result = this->add(result, comb_table[key]);
        }

        return result.to_affine(*this);
    }

    template<unsigned win_left = 4>
    void naf_precompute(const point& base, jacobian_point (&table)[1 << (win_left - 2)]) const {
        jacobian_point base_doubled = this->twice(jacobian_point(base));
        table[0] = jacobian_point(base);

        for (unsigned i = 3; i < (1 << (win_left - 1)); i += 2) {
            unsigned index = i / 2;
            table[index] = this->add(base_doubled, table[index - 1]);
        }

    }

    template<unsigned win_left = 4>
    jacobian_point mul_scalar(jacobian_point (&p)[1 << (win_left - 2)], const integer_type& multiplier) const {
        jacobian_point result = jacobian_point::inf;

        short naf_table[field_type::bits + 1];
        unsigned naf_length = naf<win_left, integer_type>(multiplier, naf_table);

        for (unsigned i = naf_length; i > 0; i--) {
            result = this->twice(result);

            short ki = naf_table[i-1];
            if (ki != 0) {
                if (ki > 0) {
                    result = this->add(result, p[ki/2]);
                } else {
                    result = this->sub(result, p[-ki/2]);
                }
            }
        }

        return result;
    }

    template<unsigned win_left = 4, unsigned win_right = 4>
    jacobian_point add_mul(
            jacobian_point (&left)[1 << (win_left - 2)], const integer_type& mul_left,
            jacobian_point (&right)[1 << (win_right - 2)], const integer_type& mul_right
    ) const {
        jacobian_point result = jacobian_point::inf;

        short naf_table_left[field_type::bits + 1];
        short naf_table_right[field_type::bits + 1];

        std::fill(std::begin(naf_table_left), std::end(naf_table_left), 0);
        std::fill(std::begin(naf_table_right), std::end(naf_table_right), 0);

        unsigned naf_length = naf<win_left, integer_type>(mul_left, naf_table_left);
        naf_length = std::max(naf<win_right, integer_type>(mul_right, naf_table_right));

        for (unsigned i = naf_length; i > 0; i--) {
            result = this->twice(result);

            short ki = naf_table_left[i-1];
            if (ki != 0) {
                if (ki > 0) {
                    result = this->add(result, left[ki/2]);
                } else {
                    result = this->sub(result, left[-ki/2]);
                }
            }

            ki = naf_table_right[i-1];
            if (ki != 0) {
                if (ki > 0) {
                    result = this->add(result, right[ki/2]);
                } else {
                    result = this->sub(result, right[-ki/2]);
                }
            }
        }

        return result;
    }


    friend std::ostream& operator<<(std::ostream& out, const point& p) {
        if (p == point::inf) {
            out << "(inf, inf)";
        } else {
            out << '(' << p.x << ", " << p.y << ')';
        }
        return out;
    }

    friend std::ostream& operator<<(std::ostream& out, const jacobian_point& p) {
        if (p == jacobian_point::inf) {
            out << "(inf, inf, inf)";
        } else {
            out << '(' << p.x << ", " << p.y << ", " << p.z << ')';
        }
        return out;
    }
};

template <typename _integer_type, typename _double_integer_type, typename _pm_integer_type>
const typename elliptic_curve<_integer_type, _double_integer_type, _pm_integer_type>::point elliptic_curve<_integer_type, _double_integer_type, _pm_integer_type>::point::inf(1, 0);

template <typename _integer_type, typename _double_integer_type, typename _pm_integer_type>
const typename elliptic_curve<_integer_type, _double_integer_type, _pm_integer_type>::jacobian_point elliptic_curve<_integer_type, _double_integer_type, _pm_integer_type>::jacobian_point::inf(1, 1, 0);
}

#endif // ELLIPTIC_CURVE_H
