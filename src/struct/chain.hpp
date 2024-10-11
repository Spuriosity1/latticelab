#pragma once
#include <concepts>
#include <map>
#include <type_traits>
#include <vector>

/*
#include <armadillo>
typedef arma::ivec3 ipos_t;
typedef arma::imat33 imat33_t;
*/
#include "vec3.hpp"
typedef vector3::vec3<long long int> ipos_t;
typedef vector3::mat33<long long int> imat33_t;
typedef unsigned int idx_t;


// Data Storage class (inherit from these for physical simulations)
struct GeometricObject{
	ipos_t position;
};



template<typename T>
using Chain = std::map<T*, int>;

// The cells over which r-chains are defined in 3D
// Boundary gives an [order-1]-chain
// Coboundary gives an [order+1]-chain corresponding to all order+1
// chains that this cell participates in
template <int order, typename boundary_T, typename coboundary_T>
struct Cell : public GeometricObject {
	static_assert(boundary_T::order == order-1, "boundary must be an n-1 cell");
	static_assert(coboundary_T::order == order+1, "coboundary must be an n+1 cell");

	Chain<boundary_T*>  boundary;
	Chain<coboundary_T*> coboundary;
};

template <int _order>
struct NullCell {
	const static int order=_order;
};

template <typename cell_t, typename boundary_cell_t>
requires (cell_t::order > 1) && (boundary_cell_t::order == cell_t::order-1)
Chain<boundary_cell_t> d(const Chain<cell_t>& chain) {
	// computes a sum over the cells
	Chain<boundary_cell_t> retval;
	for (const auto& [cell, mult] : chain){
		for (const auto& [cell_b, mult_b] : cell->boundary) {
			retval[cell_b] += mult_b*mult;
		}
	}
	return retval;
}

template <typename cell_t, typename cobounday_cell_t>
requires (cell_t::order > 1) && (cobounday_cell_t::order == cell_t::order+1)
Chain<cobounday_cell_t> co_d(const Chain<cell_t>& chain) {
	Chain<cobounday_cell_t> retval;
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

template <typename T, int _order>
concept CellLike = requires {
	{ T::order } -> std::same_as<const int>;
	requires T::order == _order;
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



