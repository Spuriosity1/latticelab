#include "rationalmath.hpp"
#include <limits>
#include <stdexcept>
namespace rational {


	Rational& Rational::operator+=(const Rational& other){
		int64_t tmp = static_cast<int64_t>(num)*other.denom + static_cast<int64_t>(other.num)*denom;
		denom = denom * other.denom;
		num = tmp;
		return *this;
	}

	Rational& Rational::operator-=(const Rational& other){
		int64_t tmp = static_cast<int64_t>(num)*other.denom - static_cast<int64_t>(other.num)*denom;
		denom = denom * other.denom;
		num = tmp;
		return *this;
	}

	void Rational::simplify(){
		int64_t this_gcd = std::gcd(num, denom);
		num /= this_gcd;
		denom /= this_gcd;
	}

	void Rational::make_denom_positive(){
		if (denom < 0){
			denom = -denom;
			num = -num;
		}
	}

	Rational& Rational::operator/=(int64_t x) {
		if (x == 0) {
			throw std::domain_error("Attempted division by zero");
		}
		denom *= x;
		return *this;
	}


	Rational& Rational::operator/=(const Rational& x) {
		if (x.num == 0) {
			throw std::domain_error("Attempted division by zero");
		}
		denom *= x.num;
		num *= x.denom;
		return *this;
	}

	Rational& Rational::operator*=(int64_t x) {
		num *= x;
		return *this;
	}

	Rational& Rational::operator*=(const Rational& x) {
		num *= x.num;
		denom *= x.denom;
		return *this;
	}

	Rational Rational::operator/(int64_t x) const {
		if (x == 0) {
			throw std::domain_error("Attempted division by zero");
		}
		return Rational(num, denom*x);
	}



typedef vector3::vec3<Rational> rvec3;
typedef vector3::mat33<Rational> rmat33;

// LINEAR ALGEBRA
void rswap(rmat33& A, int row_i, int row_j, rvec3& b){
	for (int j=0; j<3; j++){
		std::swap(A(row_i, j), A(row_j, j));
	}
	std::swap(b(row_i), b(row_j));
}

// Performs A[row_i] <- a*A[row_j] - A[row_i]
void rsub(rmat33& A, Rational a_j, int row_j, int row_i, rvec3& b){
	for (int j=0; j<3; j++){
		A(row_i, j) -= a_j*A(row_j, j);
	}
	b(row_i) -= a_j*b(row_j);
}

void rmult(rmat33& A, int row, Rational a, rvec3& b){
	for (int j=0; j<3; j++){
		A(row, j) *= a;
	}
	b(row) *= a;
}


// Solves Ax = b using Gaussian elimination
void rlinsolve(rvec3& x, const rmat33& A, const rvec3& b){
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
rmat33 inv(const rmat33& A){
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



// Returns an integer Z such that final state X
// x = Z + X
// where X has +ve denominator, and x.num is in [0, x.denom)
int64_t make_proper(Rational& x){
	if (x.denom < 0) {
		x.denom *= -1;
		x.num *= -1;
	}
#ifdef DEBUG
	if (x.denom == 0) {
		throw std::domain_error("Divide by zero in make_proper");
	}
#endif
	auto res = moddiv(x.num, x.denom);
	assert(res.rem >= 0);
	assert(res.quot*x.denom + res.rem == x.num);
	x.num = res.rem;
	return res.quot;
}




/// Approximant finding (thanks chatgpt)
///

// Function to find the lower and upper rational approximations
Rational find_nearest_rational(double v, int64_t max_denominator, int64_t max_iters) {
	int sign = (v > 0) ? 1 : -1;
	v *= sign;

	int64_t Aprev[2] = {1,0};
	int64_t Bprev[2] = {0,1};

	double x = v;
	int64_t flx = floor(x);
	x -= flx;

	int n_iter=0;
	while(n_iter < max_iters && std::abs(x) > std::numeric_limits<double>::epsilon()){
		x = 1/x;

		int64_t A = Aprev[0] + flx * Aprev[1];
		int64_t B = Bprev[0] + flx * Bprev[1];
		
		if (A > max_denominator){
			break;
		}
		Aprev[0] = Aprev[1];
		Aprev[1] = A;
		Bprev[0] = Bprev[1];
		Bprev[1] = B;

		flx = floor(x);
		x -= flx;

		n_iter++;
	}
	if (n_iter == max_iters) {
		throw std::domain_error("Rational approximant algorithm failed to converge");
	}

	return Rational(Bprev[1]*sign, Aprev[1]);
}



};

