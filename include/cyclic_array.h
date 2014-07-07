#ifndef CYCLIC_ARRAY_H
#define CYCLIC_ARRAY_H

#include <array>

namespace gost_ecc {

/**
 * @brief Auxiliary class, which encapsulates array with indexing wrapping around N.
 *
 * Extends std::array class from C++11 standard library.
 *
 * Mostly useful in recurrent elgorithms which require access to N previous values of some sequence,
 * such as Euclidean algorithm.
 */
template <typename T, std::size_t N>
class cyclic_array {
    typedef std::array<T, N> array;
	typedef typename array::size_type size_type;
	typedef typename array::reference reference;
	typedef typename array::const_reference const_reference;

	array _data;
public:
    /**
     * Curly braces constructor (C++11).
     *
     * Usage: cyclic_array<int, 3>{1, 2, 3}
     */
    template<typename... _T>
    cyclic_array(const _T& ... values)
        :_data({{values...}})
    {}

    reference
    operator[](size_type __n)
    {
        __n = (__n >= 0) ? (__n % N) : ( (__n % N) + N );
        return _data[__n];
    }

    constexpr const_reference
    operator[](size_type __n) const
    {
        __n = (__n > 0) ? (__n % N) : ( (__n % N) + N );
		return _data[__n];
    }

};
}

#endif // CYCLIC_ARRAY_H
