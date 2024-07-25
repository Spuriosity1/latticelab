#include <cstdio>
#include <cstdlib>
#include <armadillo>
#include <vector>
#include "struct/geometry.hpp"
using namespace std;



int main (int argc, char *argv[]) {
	if (argc < 5 ){
		printf("Usage: %s L1 L2 L3 max_num_spinons", argv[0]);
	}
	 
	pyro_unitcell cell(atoi(argv[1]),atoi(argv[2]),atoi(argv[3]));

	int max_num_spinons = atoi(argv[4]);

	// Setup complete
	



	return 0;
}
