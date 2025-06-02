#include <cell_geometry.hpp>
#include <preset_cellspecs.hpp>
#include <chrono>

using namespace CellGeometry;
using namespace std;


using time_pt=std::chrono::time_point<std::chrono::steady_clock>;


static void clobber() {
  asm volatile("" : : : "memory");
}


void print_dt(time_pt start, time_pt end, size_t n_samples){
	auto res = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
	cout << "\tElapsed time: " << res*1.0/n_samples<< "ns per unit cell\n";
}


template<typename T>
void construct_bench( int L){
	imat33_t supercell_spec = imat33_t::from_cols(
			{L,-L,-L},{-L,L,-L},{-L,-L,L});

	const size_t n_samples=10;
	cout << "Constructing supercell specified by" << supercell_spec <<std::endl;
	const auto spec = PrimitiveSpecifiers::CubicSpec();
	auto start = chrono::steady_clock::now();
	for (size_t i=0; i<n_samples; i++){
		T point_lat(spec, supercell_spec);
		clobber();
	}
	auto end = chrono::steady_clock::now();
	auto res = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
	cout << "Elapsed time: " << res*1.0/n_samples/L/L/L<< "ns per unit cell\n";
}


template <typename Container>
void benchmark(const string& label, const Container& container, size_t n_samples) {
    cout << label << endl;
    auto start = chrono::steady_clock::now();
    for (size_t i = 0; i < n_samples; ++i) {
        for ([[maybe_unused]] auto [i, p] : container) {
            clobber();
        }
    }
    auto end = chrono::steady_clock::now();
    print_dt(start, end, n_samples * container.size());
}


template <typename Container>
void benchmark_boundary(const string& label, const Container& container, size_t n_samples) {
    cout << label << endl;
    auto start = chrono::steady_clock::now();
    for (size_t i = 0; i < n_samples; ++i) {
        for (auto [i, p] : container) {
			for ([[maybe_unused]] auto [q, m] : p->boundary){
				clobber();
			}
        }
    }
    auto end = chrono::steady_clock::now();
    print_dt(start, end, n_samples * container.size());
}


template <typename Container>
void benchmark_coboundary(const string& label, const Container& container, size_t n_samples) {
    cout << label << endl;
    auto start = chrono::steady_clock::now();
    for (size_t i = 0; i < n_samples; ++i) {
        for (auto [i, p] : container) {
			for ([[maybe_unused]]auto [q, m] : p->coboundary){
				clobber();
			}
        }
    }
    auto end = chrono::steady_clock::now();
    print_dt(start, end, n_samples * container.size());
}



void iter_bench( int L){

	const size_t n_samples=100;
	imat33_t supercell_spec = imat33_t::from_cols(
			{L,-L,-L},{-L,L,-L},{-L,-L,L});

	const auto spec = PrimitiveSpecifiers::CubicSpec();
	auto lat =  PeriodicVolLattice<Cell<0>,Cell<1>,Cell<2>,Cell<3>>(spec, supercell_spec);
	
	benchmark("Point enumeration", lat.points, n_samples);
	benchmark("Link enumeration", lat.links, n_samples);
	benchmark("Plaq enumeration", lat.plaqs, n_samples);
	benchmark("Vol enumeration", lat.vols, n_samples);
	cout<<"\n";

	benchmark_boundary("Link boundary enumeration", lat.links, n_samples);
	benchmark_coboundary("Link coboundary enumeration", lat.links, n_samples);


}

int main (int argc, char *argv[]) {
	assert(argc >= 2);
	int L = atoi(argv[1]);
	cout << "AbstractLat\n";
	construct_bench<PeriodicAbstractLattice>(L);
	cout << "PointLat\n";
	construct_bench<PeriodicPointLattice<Cell<0>>>(L);
	cout << "LinkLat\n";
	construct_bench<PeriodicLinkLattice<Cell<0>,Cell<1>>>(L);
	cout << "PlaqLat\n";
	construct_bench<PeriodicPlaqLattice<Cell<0>,Cell<1>,Cell<2>>>(L);
	cout << "VolLat\n";
	construct_bench<PeriodicVolLattice<Cell<0>,Cell<1>,Cell<2>,Cell<3>>>(L);

	iter_bench(L);
	return 0;
}
