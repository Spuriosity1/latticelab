#include <SortedVectorMap.hpp>
#include <iostream>
#include <map>
#include <unordered_map>
#include <chrono>
#include <vector>

using namespace std;

struct Spin{
	int Sz;
};

template <typename T> 
size_t time_enumeration(const T& container){
	auto start = chrono::steady_clock::now();
	for (const auto& [J, s] : container){
		s->Sz*=-1;
	}
	auto end = chrono::steady_clock::now();
	return chrono::duration_cast<chrono::nanoseconds>(end - start).count();
}


size_t time_vector_enumeration(const std::vector<Spin*>& container){
	auto start = chrono::steady_clock::now();
	for (size_t i=0; i<container.size(); i++){
		container[i]->Sz*=-1;
	}
	auto end = chrono::steady_clock::now();
	return chrono::duration_cast<chrono::nanoseconds>(end - start).count();
}


template <typename T> 
size_t time_creation(T& container, size_t n_sites){
	auto start = chrono::steady_clock::now();
	for (size_t i=0; i<n_sites; i++){
		container[i] = new Spin;
	}
	auto end = chrono::steady_clock::now();
	return chrono::duration_cast<chrono::nanoseconds>(end - start).count();
}


size_t time_vector_creation(std::vector<Spin*>& container, size_t n_sites){
	auto start = chrono::steady_clock::now();
	container.resize(n_sites);
	for (size_t i=0; i<n_sites; i++){
		container[i] = new Spin;
	}
	auto end = chrono::steady_clock::now();
	return chrono::duration_cast<chrono::nanoseconds>(end - start).count();
}


int main (int argc, char *argv[]) {
	// set it up
	size_t n_sites = atoi(argv[1]);

	std::map<size_t, Spin*> map_v;
	std::unordered_map<size_t, Spin*> umap_v;
	SortedVectorMap<size_t, Spin*> svm_v;
	std::vector<Spin*> v;

	cout<<"Creation:\n";
	cout<<"\tMap:       "<<time_creation(map_v, n_sites) <<" ns\n";
	cout<<"\tUMap:      "<<time_creation(umap_v, n_sites) <<" ns\n";
	cout<<"\tSortedVec: "<<time_creation(svm_v, n_sites) <<" ns\n";
	cout<<"\tVec:       "<<time_vector_creation(v, n_sites) <<" ns\n";


	cout<<"Enumeration:\n";
	cout<<"\tMap:       "<<time_enumeration(map_v) <<" ns\n";
	cout<<"\tUMap:      "<<time_enumeration(umap_v) <<" ns\n";
	cout<<"\tSortedVec: "<<time_enumeration(svm_v) <<" ns\n";
	cout<<"\tVec:       "<<time_vector_enumeration(v) <<" ns\n";

	
	return 0;
}


//   Output: (n=10 000)
//////////////////////////////////////
// Creation:
//         Map:       1166500 ns
//         UMap:      665041 ns
//         SortedVec: 786167 ns
//         Vec:       232166 ns
// Enumeration:
//         Map:       62291 ns
//         UMap:      25291 ns
//         SortedVec: 12459 ns
//         Vec:       11250 ns
//                                                                                                    
		
