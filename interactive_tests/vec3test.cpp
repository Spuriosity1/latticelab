#include <iostream>
#include <sstream>
#include "rationalmath.hpp"


using namespace rational;

int main (int argc, char *argv[]) {
	if (argc < 3){
		std::cout << "USAGE: vec3test r1 r2 r3, where r# are in format e.g. 11/7\n";
		return 1;
	}
	
	Rational r[3] = { 0,0,0};
	for (int i=0; i<3; i++){
		std::istringstream iss(argv[i+1]);
		iss >> r[i];
	}
	std::cout << "read in " << r[0] <<" "<< r[1] <<" "<< r[2] <<std::endl;
	rvec3 rv(r[0], r[1], r[2]);

	std::cout << "read in " << rv[0] <<" "<< rv[1] <<" "<< rv[2] <<std::endl;
	
	return 0;
}
