#pragma once
#include "modulus.hpp"
#include "vec3.hpp"
#include "modulus.hpp"

#include <array>
#include <cassert>
#include <concepts>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <nlohmann/json.hpp>
namespace rational {


struct Rational {
	int64_t num;
	int64_t denom;

	Rational(int64_t _num) : num(_num), denom(1){};

	Rational(int64_t _num, int64_t _denom) : num(_num), denom(_denom) {
	if (denom == 0) {
		throw std::domain_error("Attempting to construct rational with zero denominator");
	}
	};

	Rational& operator+=(const Rational& other){
		int64_t tmp = static_cast<int64_t>(num)*other.denom + static_cast<int64_t>(other.num)*denom;
		denom = denom * other.denom;
		num = tmp;
		return *this;
	}

	Rational& operator-=(const Rational& other){
		int64_t tmp = static_cast<int64_t>(num)*other.denom - static_cast<int64_t>(other.num)*denom;
		denom = denom * other.denom;
		num = tmp;
		return *this;
	}

	void simplify(){
		int64_t this_gcd = std::gcd(num, denom);
		num /= this_gcd;
		denom /= this_gcd;
	}

	void make_denom_positive(){
		if (denom < 0){
			denom = -denom;
			num = -num;
		}
	}

	Rational& operator/=(int64_t x) {
		if (x == 0) {
			throw std::domain_error("Attempted division by zero");
		}
		denom *= x;
		return *this;
	}


	Rational& operator/=(const Rational& x) {
		if (x.num == 0) {
			throw std::domain_error("Attempted division by zero");
		}
		denom *= x.num;
		num *= x.denom;
		return *this;
	}

	Rational& operator*=(int64_t x) {
		num *= x;
		return *this;
	}

	Rational& operator*=(const Rational& x) {
		num *= x.num;
		denom *= x.denom;
		return *this;
	}

	Rational operator/(int64_t x) const {
		if (x == 0) {
			throw std::domain_error("Attempted division by zero");
		}
		return Rational(num, denom*x);
	}

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
inline void rswap(rmat33& A, int row_i, int row_j, rvec3& b){
	for (int j=0; j<3; j++){
		std::swap(A(row_i, j), A(row_j, j));
	}
	std::swap(b(row_i), b(row_j));
}

// Performs A[row_i] <- a*A[row_j] - A[row_i]
inline void rsub(rmat33& A, Rational a_j, int row_j, int row_i, rvec3& b){
	for (int j=0; j<3; j++){
		A(row_i, j) -= a_j*A(row_j, j);
	}
	b(row_i) -= a_j*b(row_j);
}

inline void rmult(rmat33& A, int row, Rational a, rvec3& b){
	for (int j=0; j<3; j++){
		A(row, j) *= a;
	}
	b(row) *= a;
}

inline Rational inv(const Rational& a){	
	if (a.denom == 0) {
		throw std::domain_error("Attempted inverse of 0");
	}
	return Rational(a.denom, a.num);
}

// Solves Ax = b using Gaussian elimination
inline void rlinsolve(rvec3& x, const rmat33& A, const rvec3& b){
	rmat33 B(A);
	x = b;
	
	for (int col = 0; col<3; col++){
		// pemute rows until A(col, col) != 0
		int row = col + 1;
		while ((B(col,col) == 0)) {
			if (row > 3){
				throw std::invalid_argument("Matrix is singular");
			}

			rswap(B, col, row, x);
			row += 1;
		}

		// rmult to make B(col, col) == 1
		rmult(B, col, inv(B(col,col)), x);
		// simplify the row
		for (int ci=col; ci<3; ci++){
			B(col, ci).simplify();
		}

		// delete the lower part
		for (row = col + 1; row < 3; row++){
			rsub(B, B(row, col), col, row, x);
		}

		// Simplify
		for (int i=0; i<9; i++){
			B[i].simplify();
		}
	}


	// matrix shold now be upper triangular with 1s on diagonal
	// Remove top part and simplify answer
	for (int col=2; col>=0; col--){
		for (int row=0; row<col; row++){
			rsub(B, B(row,col), col, row, x);
		}
		x[col].simplify();
	}

	assert(A * x == b);

}


// Computes the closed-form, simplified matrix inverse
inline rmat33 inv(const rmat33& A){
	auto denom = vector3::det(A);
	std::array<Rational, 3> b0 = {-A(1, 2) * A(2, 1) + A(1, 1) * A(2, 2),
		A(0, 2) * A(2, 1) - A(0, 1) * A(2, 2),
		-A(0, 2) * A(1, 1) + A(0, 1) * A(1, 2)};
	std::array<Rational, 3> b1 = {A(1, 2) * A(2, 0) - A(1, 0) * A(2, 2),
		-A(0, 2) * A(2, 0) + A(0, 0) * A(2, 2),
		A(0, 2) * A(1, 0) - A(0, 0) * A(1, 2)};
	std::array<Rational, 3> b2 = {-A(1, 1) * A(2, 0) + A(1, 0) * A(2, 1),
		A(0, 1) * A(2, 0) - A(0, 0) * A(2, 1),
		-A(0, 1) * A(1, 0) + A(0, 0) * A(1, 1)};

	auto retval =  inv(denom) * rmat33::from_rows(b0,b1,b2);
	for (int i=0; i<9; i++){
		retval[i].simplify();

#ifndef NDEBUG
	if (retval[i].denom == 0) {
		throw std::domain_error("Inverse of singular matrix");
	}
#endif

	}

	return retval;
}


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
inline int64_t make_proper(Rational& x){
	if (x.denom < 0) {
		x.denom *= -1;
		x.num *= -1;
	}
#ifdef DEBUG
	if (x.denom == 0) {
		throw std::domain_error("Divide by zero in make_proper");
	}
#endif
	auto [quot, rem] = moddiv(x.num, x.denom);
	assert(rem >= 0);
	assert(quot*x.denom + rem == x.num);
	x.num = rem;
	return quot;
}



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


static_assert(std::convertible_to<int, Rational>, "int cannot be implicitly converted");




/////////JSON conversions
using json=nlohmann::json;

inline void to_json(json& j, const Rational& r) {
	Rational r2 = r;
	r2.simplify();
	r2.make_denom_positive();
	j = json{ {r.num, r.denom} };
}

inline void from_json(const json& j, Rational& r) {
	j.at(0).get_to(r.num);
	j.at(1).get_to(r.denom);
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

    return hash<int64_t>()(k.num / k.denom) ;
  }
};
