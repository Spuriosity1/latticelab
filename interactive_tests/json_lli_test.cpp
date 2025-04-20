#include "nlohmann/json.hpp"
#include "vec3.hpp"
#include <iostream>
#include <sstream>


using namespace vector3;
using namespace nlohmann;

// Tests long long int deserailisation (in hex)
int main (int argc, char *argv[]) {
	if (argc < 4){
		std::cout << "USAGE: "<<argv[0]<<" v1 v2 v3\n";
	}

	vec3<long long int> v;


	for (int i=0; i<3; i++){
		std::stringstream ss;
		ss << std::hex;
		ss << argv[i+1];
		ss >> v[i];
	}

	json j(v);
	std::cout << "raw :" << v << "\n";
	std::cout << "from json: " << j << "\n";
	
	return 0;
}
