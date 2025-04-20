#include "nlohmann/json.hpp"
#include <cstdlib>
#include <gtest/gtest.h>

#include <limits>
#include <stdexcept>
#include <rationalmath.hpp>

#ifndef DEBUG
#warning "Compiling without DEBUG defined; bounds checking will not throw"
#endif


using json = nlohmann::json;

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

TEST(RationalTest,  NearestRational){
	double x = 0.333;
	EXPECT_EQ(find_nearest_rational(x,100), Rational(1,3));
	EXPECT_EQ(find_nearest_rational(x,999), Rational(1,3));
	EXPECT_EQ(find_nearest_rational(x,1000), Rational(333,1000));


	x = -0.05882;
	EXPECT_EQ(find_nearest_rational(x,100), Rational(-1,17));
	EXPECT_EQ(find_nearest_rational(x,1000), Rational(-1,17));
	EXPECT_EQ(find_nearest_rational(x,10000), Rational(-1,17));
	EXPECT_EQ(find_nearest_rational(x,100000), Rational(-5882,100000));

	x = sqrt(2);
	EXPECT_EQ(find_nearest_rational(x,100), Rational(99,70));
	EXPECT_EQ(find_nearest_rational(x,10000), Rational(8119,5741));
	EXPECT_THROW(find_nearest_rational(x,std::numeric_limits<int64_t>::max()), std::domain_error);
}

TEST(RationalTest, JSONify){
	Rational r(-100, 501);
	json j = r;
	std::cerr << j;
	auto r2 = j.template get<rational::Rational>();
	EXPECT_EQ(r,r2);
	EXPECT_EQ(j, json(-0.1996007984031936));
}
