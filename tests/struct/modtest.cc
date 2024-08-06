#include <cstdlib>
#include <gtest/gtest.h>
#include <struct/modulus.hpp>

std::lldiv_t moddiv_result(long long quot, long long rem){
	std::lldiv_t res;
	res.quot = quot;
	res.rem = rem;
	return res;
}

bool operator==(const std::lldiv_t& x1, const std::lldiv_t x2){
	return x1.quot==x2.quot && x1.rem==x2.rem;
}

TEST(ModTest, ModDivCorrectnessLL) { 
	EXPECT_EQ(moddiv(10ll, 3ll),  moddiv_result(3, 1));
	EXPECT_EQ(moddiv(-10ll, 3ll),  moddiv_result(-4, 2));
	EXPECT_EQ(moddiv(-10ll, 10ll),  moddiv_result(-1, 0));
	EXPECT_EQ(moddiv(-20ll, 10ll),  moddiv_result(-2, 0));
}
