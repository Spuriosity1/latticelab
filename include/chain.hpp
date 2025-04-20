#pragma once
#include <ios>
#include <map>
#include <vector>
#include <concepts>
/*
#include <armadillo>
typedef arma::ivec3 ipos_t;
typedef arma::imat33 imat33_t;
*/
#include "vec3.hpp" 
#include "rationalmath.hpp"

typedef vector3::vec3<rational::Rational> ipos_t;
typedef vector3::vec3<long long int> ivec3_t;
typedef vector3::mat33<long long int> imat33_t;
typedef unsigned int idx_t;

template <int order>
struct Cell;

// Implementation of an n-chain: a std::map realises a sparse vector of int
template <int order>
using Chain = std::map<Cell<order>*, int>;

template <int order>
inline void cleanup_chain(Chain<order>& c){
	// delete any canceled cells
	for (auto it = c.cbegin(); it != c.cend(); ) {
		if ( it->second == 0) {
			c.erase(it++);
		} else {
			++it;
		}
	}
}

template<int order>
inline Chain<order> operator+(const Chain<order> c1, const Chain<order> c2){
	auto retval = c1;
	for (const auto& cell : c2){
		if(c1.find(cell) == c1.end()){
			retval[cell] = c2[cell];
		} else {
			retval[cell] += c2[cell];
		}
	}
	cleanup_chain(retval);
		
	return retval;
}

template<int order>
inline bool eq_superset(const Chain<order>& c1, const Chain<order>& c2){
	for (const auto& [key, value] : c1){
		if (value == 0) continue;
		auto it = c2.find(key);
		if (it == c2.end() || it->second == 0) return false;
		if (it->second != value) return false;
	}
	return true;
}

template<int order>
inline bool operator==(const Chain<order>& c1, const Chain<order>& c2){
	// actually very annoying to do
	return eq_superset(c1, c2) && eq_superset(c2, c1); 
}


template<int order>
inline Chain<order> operator-(const Chain<order>& c1, const Chain<order>& c2){
	auto retval = c1;
	for (const auto& cell : c2){
		if(c1.find(cell) == c1.end()){
			retval[cell] = c2[cell];
		} else {
			retval[cell] -= c2[cell];
		}
	}
	cleanup_chain(retval);	
	return retval;
}

template<int order>
inline Chain<order> operator+=(Chain<order>& c, const Cell<order>& cell){
	if (c.find(&cell) == c.end()){
		c[&cell] = 1;
	} else {
		c[&cell] += 1;
	}
	cleanup_chain(c);
	return c;
}

template<int order>
inline Chain<order> operator-=(Chain<order>& c, const Cell<order>& cell){
	if (c.find(&cell) == c.end()){
		c[&cell] = -1;
	} else {
		c[&cell] -= 1;
	}

	cleanup_chain(c);
	return c;
}

template<int order>
Chain<order> operator*(int x, const Chain<order> c){
	if (x == 0) {
		return Chain<order>{};
	}
	auto retval = c;
	for (auto& [cell, m] : retval) {
		retval[cell] *= x;
	}

	return retval;
}

template<int order>
std::ostream &operator<<(std::ostream &stream, const Chain<order> &c) {
	for (const auto& [cell, m] : c){
		if (m==0) continue;
		stream<< std::showpos << m <<" "<<cell->position;
	}
	return stream;
}

// Data Storage class (inherit from these for physical simulations)
struct GeometricObject{
	ipos_t position;
};

// The cells over which r-chains are defined in 3D
// Boundary gives an [order-1]-chain
// Coboundary gives an [order+1]-chain corresponding to all order+1
// chains that this cell participates in
template <int order>
struct Cell : public GeometricObject {
	Chain<order - 1> boundary;
	Chain<order + 1> coboundary;
};

template <>
struct Cell<0> : public GeometricObject {
	Chain<1> coboundary;
};

template <>
struct Cell<3> : public GeometricObject {
	Chain<2> boundary;
};

template <int order>
requires (order > 0)
Chain<order-1> d(const Chain<order>& chain) {
	// computes a sum over the cells
	Chain<order-1> retval;
	for (const auto& [cell, mult] : chain){
		for (const auto& [cell_b, mult_b] : cell->boundary) {
			retval[cell_b] += mult_b*mult;
		}
	}
	return retval;
}

template <int order>
requires (order < 3)
Chain<order+1> co_d(const Chain<order>& chain) {
	Chain<order+1> retval;
	for (const auto& [cell, mult] : chain){
		for (const auto& [cell_b, mult_b] : cell->coboundary) {
			retval[cell_b] += mult_b*mult;
		}
	}
	return retval;
}


// An ordered pair of an integer multiplier and a relative position
// (pointing to the n-1 cell)
struct VectorSignPair {
	int multiplier;
	ipos_t relative_position;
};

/** Cell-specifier class: A lightweight contianer for only a position 
 * and a list of boundary members.
 * The user should avoid overriding this class.
*/
template<int _order>
struct CellSpecifier : public GeometricObject {
	// Coordinates of the objects on the boundary relative to `position`.
	static constexpr int order = _order;
	std::vector<VectorSignPair> boundary;
};


template <typename T, int order>
concept CellLike = std::derived_from<T, Cell<order>>;

