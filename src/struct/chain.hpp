#pragma once
#include <map>
#include <vector>

typedef unsigned int idx_t;
/*
#include <armadillo>

typedef arma::ivec3 ipos_t;
typedef arma::imat33 imat33_t;
*/
#include "vec3.hpp"
typedef vector3::vec3<long long int> ipos_t;
typedef vector3::mat33<long long int> imat33_t;

template <int order>
struct Cell;

template <int order>
using Chain =  std::map<Cell<order>*, int>;

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
requires (order > 1)
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



