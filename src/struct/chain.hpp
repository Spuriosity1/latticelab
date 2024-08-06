#pragma once


#include <armadillo>

typedef unsigned int idx_t;
typedef arma::ivec3 ipos_t;

// Cursed chicken-egg shenanigans:
// Need to forward declare this
template <int order>
struct CellMultPair;

template <int order>
using Chain =  std::vector<CellMultPair<order>>;

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
struct CellMultPair{
	int multiplier; // generically +- 1
	const Cell<order>* cell;
};

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



