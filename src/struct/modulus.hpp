#pragma once
#include <cstdlib>
#include <type_traits>
#include <concepts>
#include <armadillo>

// like std::div with rem wrapped to [0, base)
// Should always satisfy x = base * res.quot + rem
// Assumes positive base; disable asserion using -DNDEBUG
template<typename V>
requires std::signed_integral<V>
std::lldiv_t moddiv(V x, V base) {
	std::lldiv_t res = std::lldiv(x, base);
	assert(base > 0);
	if (res.rem<0) {
		// moddiv(-1, 5)
		// -> quot: 0, rem: -1, 
		res.rem += base;
		--res.quot;
	}
	return res;
}


template <typename T>
concept has_static_size = requires(T) {
        { T::size() } -> std::convertible_to<size_t>;
    };

template<typename T>
requires has_static_size<T>
inline std::pair<T, T> moddiv(T x, T base){
	std::pair<T, T> retval;
	for (size_t i=0; i<x.size(); i++){
		std::lldiv_t tmp = moddiv(x[i], base[i]);
		retval.first[i]= tmp.quot; 
		retval.second[i] = tmp.rem;
	}
	return retval;
}


template<typename V>
requires std::signed_integral<V>
V mod(V x, V base) {
	assert(base > 0);
	V res = x % base;
	return res >= 0 ? res : res + base;

    //return (a % b + b) % b;
}



// Explicit functions preserving intness
inline arma::sword det2(arma::imat33 a){
	return a(0,0)*a(1,1) - a(0,1)*a(1,0);
}

inline arma::sword det3(arma::imat33 a){
    return (a(0,0) * (a(1,1) * a(2,2) - a(2,1) * a(1,2))
           -a(1,0) * (a(0,1) * a(2,2) - a(2,1) * a(0,2))
           +a(2,0) * (a(0,1) * a(1,2) - a(1,1) * a(0,2)));
}


