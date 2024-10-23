#include <cstdlib>
#include <gtest/gtest.h>

#include <stdexcept>
#include <struct/rationalmath.hpp>

#ifndef DEBUG
#warning "Compiling without DEBUG defined; bounds checking will not throw"
#endif


using namespace rational;

TEST(RationalTest, AddSubOpsCorrect){
	Rational r(-5,3);
	r += 1;
	EXPECT_EQ(r, Rational(-2, 3));
	r -= 1;
	EXPECT_EQ(r, Rational(-5, 3));
	r += Rational(3,2);
	EXPECT_EQ(r, Rational(-1, 6));
	EXPECT_EQ(r, Rational(2, -12));
	r -= Rational(3,2);
	EXPECT_EQ(r, Rational(-5, 3));
}

TEST(RationalTest, MulOpsCorrect){
	Rational r(-5,3);
	r *= 3;
	EXPECT_EQ(r, Rational(-5));
	r *= Rational(10,-11);
	EXPECT_EQ(r, Rational(50,11));
}

TEST(RationalTest, DivOpsCorrect){
	Rational r(12, -17);
	r /= 3;
	EXPECT_EQ(r, Rational(4,-17));
	r /= Rational(-1, 17);
	EXPECT_EQ(r, Rational(4));
}

TEST(RationalTest, DivExceptions){
	Rational r(0,6);
	EXPECT_THROW(r/=0, std::domain_error);
}


