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

        jacobian_point(const point& affine)
            :x(affine.x),y(affine.y),z(1)
        {}

        point to_affine(const elliptic_curve& curve) const {
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
        t2          = f.mul(3, t2); // A = 3 (X1 - Z1^2) (X1 + Z1^2)
        result.y    = f.mul(2, p.y); // B = 2Y1
        result.z    = f.mul(result.y, p.z);
        result.y    = f.mul(result.y, result.y); // C = B^2
        t3          = f.mul(result.y, p.x); // D = CX1
        result.y    = f.mul(result.y, result.y);
        result.y    = f.mul(result.y, this->inv_2);
        result.x    = f.mul(t2, t2);
        t1          = f.mul(2, t3);
        result.x    = f.sub(result.x, t1);
        t1          = f.sub(t3, result.x);
        t1          = f.mul(t1, t2);
        result.y    = f.sub(t1, result.y);

        return result;
    }

    /**
     * @brief Point addition for mixed Jacobian-affine coordinates.
     *
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
            return right;
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
        t1          = f.mul(2, t3);
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
     * @brief Repeated doubling algorithm.
     *
     * See: Hankerson, D., Vanstone, S., & Menezes, A. (2004). Guide to elliptic curve cryptography.
     * Page 93, alg. 3.23.
     * @param p
     * @param count
     * @return
     */
    jacobian_point repeated_twice(const jacobian_point& p, unsigned count) {
        if (p == jacobian_point::inf) {
            return jacobian_point::inf;
        }

        jacobian_point result = p;
        const field_type& f = this->field;

        integer_type a, b, w, y_squared, t1, t2;

        result.y    = f.mul(2, result.y); // Y <- 2Y
        w           = f.mul(result.z, result.z);
        w           = f.mul(w, w); // W <- Z^4

        while (count > 0) {
            a           = f.mul(result.x, result.x); // a = X^2
            a           = f.sub(a, w); // a = X^2 - W
            a           = f.mul(3, a); // a = 3 (X^2 - W)
            std::cout << "a: " << a << std::endl;

            y_squared   = f.mul(result.y, result.y);
            b           = f.mul(result.x, y_squared); // B = X Y^2
            std::cout << "b: " << a << std::endl;
            result.x    = f.sub(f.mul(a, a), f.mul(2, b)); // X = A^2 - 2B
            result.z    = f.mul(result.z, result.y); // Z = ZY
            std::cout << "x: " << result.x << std::endl;
            std::cout << "z: " << result.z << std::endl;

            count--;

            y_squared   = f.mul(y_squared, y_squared); // y_squared = Y^4

            if (count > 0) {
                w = f.mul(w, y_squared); // W = W Y^4
            }

            result.y = f.sub(b, result.x); // B - X
            result.y = f.mul(a, result.y); // A(B - X)
            result.y = f.mul(2, result.y); // 2A (B - X)
            result.y = f.sub(result.y, y_squared); // Y = 2A (B - X) - Y^4

            std::cout << "y: " << result.y << std::endl;
        }

        result.y = f.mul(result.y, this->inv_2);

        return result;
    }

    point mulScalar(const point& p, const integer_type& multiplier) const {
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

template <typename _integer_type, typename _double_integer_type>
const typename elliptic_curve<_integer_type, _double_integer_type>::point elliptic_curve<_integer_type, _double_integer_type>::point::inf(1, 0);

template <typename _integer_type, typename _double_integer_type>
const typename elliptic_curve<_integer_type, _double_integer_type>::jacobian_point elliptic_curve<_integer_type, _double_integer_type>::jacobian_point::inf(1, 1, 0);
}

#endif // ELLIPTIC_CURVE_H
