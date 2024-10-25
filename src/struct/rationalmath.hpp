#pragma once
#include "modulus.hpp"
#include "vec3.hpp"
#include "modulus.hpp"

#include <array>
#include <cassert>
#include <concepts>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <nlohmann/json.hpp>
namespace rational {


struct Rational {
	int64_t num;
	int64_t denom;

	Rational() {};

	Rational(int64_t _num) : num(_num), denom(1){};

	Rational(int64_t _num, int64_t _denom) : num(_num), denom(_denom) {
	if (denom == 0) {
		throw std::domain_error("Attempting to construct rational with zero denominator");
	}
	};

	Rational& operator+=(const Rational& other);

	Rational& operator-=(const Rational& other);

	void simplify();

	void make_denom_positive();

	Rational& operator/=(int64_t x);

	Rational& operator/=(const Rational& x); 

	Rational& operator*=(int64_t x) ;

	Rational& operator*=(const Rational& x); 

	Rational operator/(int64_t x) const; 

	bool operator==(int64_t x) const {
		return num == x * denom;
	}

	bool operator==(const Rational& r) const {
		return static_cast<int64_t>(num) * r.denom == static_cast<int64_t>(r.num)*denom;
	}
};

inline Rational operator*(int64_t x, const Rational& f){
	return Rational(x*f.num, f.denom);
}


inline Rational operator*(const Rational& x, const Rational& f){
	return Rational(x.num*f.num, x.denom*f.denom);
}

inline Rational operator+(const Rational& x, const Rational& y){
	Rational retval(x);
	retval += y;
	return retval;
}

inline Rational operator-(const Rational& x, const Rational& y){
	Rational retval(x);
	retval -= y;
	return retval;
}

inline Rational operator-(const Rational& x){
	return Rational(-x.num, x.denom);
}


typedef vector3::vec3<Rational> rvec3;
typedef vector3::mat33<Rational> rmat33;


//TODO template this BS
void rswap(rmat33& A, int row_i, int row_j, rvec3& b);


// Performs A[row_i] <- a*A[row_j] - A[row_i]
void rsub(rmat33& A, Rational a_j, int row_j, int row_i, rvec3& b);

void rmult(rmat33& A, int row, Rational a, rvec3& b);

inline Rational inv(const Rational& a){	
	if (a.denom == 0) {
		throw std::domain_error("Attempted inverse of 0");
	}
	return Rational(a.denom, a.num);
}

// Solves Ax = b using Gaussian elimination
void rlinsolve(rvec3& x, const rmat33& A, const rvec3& b);

// Computes the closed-form, simplified matrix inverse
rmat33 inv(const rmat33& A);

inline std::ostream& operator<<(std::ostream& os, const Rational& r){
	os << r.num << "/"<<r.denom;
	return os;
}

inline std::istream& operator>>(std::istream& is, Rational& r){
	is >> r.num;
	if (is.get() != '/'){
		r.denom=1;
	} else {
		is >> r.denom;
	}
	return is;
}

// Returns an integer Z such that final state X
// x = Z + X
// where X has +ve denominator, and x.num is in [0, x.denom)
int64_t make_proper(Rational& x);

template<typename S>
requires std::signed_integral<S>
inline rmat33 operator*(const rmat33& a, const vector3::mat33<S> b){
	// suboptimal but idc
	return a * rmat33::from_other(b);
}


template<typename S>
requires std::signed_integral<S>
inline rvec3 operator*(const rmat33& a, const vector3::vec3<S> b){
	rvec3 res;
	for (int i=0; i<3; i++){
		res[i] = a(i, 0)*b(0) + a(i,1)*b(1) + a(i,2)*b(2);
	}
	return res;
}


/// Approximant finding (thanks chatgpt)
///

// Function to find the lower and upper rational approximations
Rational find_nearest_rational(double v, int64_t max_denominator=1000, int64_t max_iters=100);

/////////JSON conversions
using json=nlohmann::json;
inline void to_json(json& j, const Rational& r) {
	j = static_cast<double>( r.num ) / r.denom ;
}

inline void from_json(const json& j, Rational& r) {
	double tmp;
	j.get_to(tmp);
	r = find_nearest_rational(tmp, 1000);
}

};



// HASH ability
//
template<>
struct std::hash<rational::Rational>
{
  std::size_t operator()(const rational::Rational& k) const
  {
    using std::size_t;
    using std::hash;

    // Compute individual hash values for first,
    // second and third and combine them using XOR
    // and bit shifting:

    return hash<int64_t>()(k.num * 0x10000 / k.denom) ;
  }
};
