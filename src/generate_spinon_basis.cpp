#include <cstdio>
#include <cstdlib>
#include <complex>
//#include <vector>
#include "struct/cell_geometry.hpp"
#include "struct/preset_cellspecs.hpp"
using namespace std;

#ifdef __SIZEOF_INT128__
    typedef __uint128_t state_t;
#else
    // do some fallback stuff here
#error 128-bit precision not supported, some hand-optimisation is needed
#endif

using namespace CellGeometry;

struct Tetra : public Cell<0> {
    double Q;
};

inline double abs2(const std::complex<double> c){
    return c.real()*c.real() + c.imag()*c.imag();
}

struct PyroSite : public Cell<1> {
    std::complex<double> state[2]; 
    double Sz() const { 
        return 0.5*(abs2(state[0]) - abs2(state[1]));
    }
};

// For dealing with geometry and computing the connectivity graph
// Does not store the full state
typedef PeriodicLinkLattice<Tetra, PyroSite> PeriodicLattice;



// setting up the state map: interpret 0-> down, 1-> up
// Maps to a product state, just for visualisation purposes
void setstate(PeriodicLattice& lat, state_t state){
    for (size_t i=0; i<lat.links.size(); i++){
        if (state& ((state_t)1 << i == 1)){
            lat.links[i].state[0] = 1.0;
            lat.links[i].state[1] = 0.0;
        } else {
            lat.links[i].state[0] = 0.0;
            lat.links[i].state[1] = 1.0;
        }
    }

    for (auto& t : lat.points){
        t.Q = 0;
        for (auto& [pyro, mul] : t.coboundary) {
            const PyroSite* ps = static_cast<const PyroSite*>(pyro);
            t.Q += ps->Sz();
        }
    }
}


int main (int argc, char *argv[]) {
	if (argc < 2 ){
		printf("Usage: %s max_num_spinons", argv[0]);
	}
    // Import and construct the cell specifier (this can be done statically
    // in principle, but there is very little point optimising it)
    const PrimitiveSpecifers::Pyrochlore pyrospec;

    // Construct the periodic lattice extension using these extensions
    PeriodicLattice lat(pyrospec, imat33_t::from_cols({2,0,0},{0,2,0},{0,0,2}));

	int max_num_spinons = atoi(argv[1]);

	// Idea: Work in local sigma-Z basis
	// Represent states as 128-bit integers
	


	return 0;
}
