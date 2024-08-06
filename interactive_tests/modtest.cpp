#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <struct/modulus.hpp>

int main (int argc, char *argv[]) {
	assert(argc == 3);
	long long num = atoll(argv[1]);
	long long denom = atoll(argv[2]);
	auto res = std::lldiv(num, denom);
	auto res2 = moddiv(num, denom);
	printf("%lld std::lldiv %lld  -> quot: %lld  rem: %lld\n", num, denom, res.quot, res.rem);
	printf("%lld moddiv %lld  -> quot: %lld  rem: %lld\n", num, denom, res2.quot, res2.rem);
	return 0;
}
