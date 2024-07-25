#pragma once 


#include <armadillo>
#include <vector>
#include <type_traits>


typedef unsigned int idx_t;
typedef arma::ivec3 ipos_t;

// forward declaration
template <int order>
struct Chain;


// The cells over which r-chains are defined in 3D
// Boundary gives an [order-1]-chain
// Coboundary gives an [order+1]-chain corresponding to all order+1
// chains that this cell participates in
template <int order>
struct Cell {
	ipos_t position;
	Chain<order - 1> boundary;
	Chain<order + 1> coboundary;
};

template <int order>
struct Chain {
	std::vector<const Cell<order>*> members;
	std::vector<int> coeff;
};

template <>
struct Cell<0> {
	ipos_t position;
	Chain<1> coboundary;
};

template <>
struct Cell<3> {
	ipos_t position;
	Chain<2> boundary;
};

/**
 * This arcane magic actually has a simple purpose - 
 * it insures that whatever pointlike or plaquettelike objects are present,
 * they always have position, boundary and coboundary operators.
 * The programmer is invited to override these as needed to store additional
 * information. For example, giving an Ising defree of freedom to the points
 * may be done with
 * struct Spin : Point {
 *     int ising = 1;
 * };
 *
 * PrimitiveCell<Spin> spincell; 
 *
 */
template<
	std::derived_from<Cell<3>> Vol   = Cell<3>,
	std::derived_from<Cell<2>> Plaq  = Cell<2>,
	std::derived_from<Cell<1>> Link  = Cell<1>,
	std::derived_from<Cell<0>> Point = Cell<0> 
	>
struct PrimitiveCell {	
	arma::imat33 primitive_vectors;
	std::vector<Point> points;
	std::vector<Link> links;
	std::vector<Plaq> plaqs;
	std::vector<Vol> volumes;
};


template<
	std::derived_from<Cell<0>> Point = Cell<0>, 
	std::derived_from<Cell<1>> Link  = Cell<1>,
	std::derived_from<Cell<2>> Plaq  = Cell<2>,
	std::derived_from<Cell<3>> Vol   = Cell<3>
	>
struct PeriodicLattice : PrimitiveCell<Point, Link, Plaq, Vol> {
	PeriodicLattice(const PrimitiveCell<Point, Link, Plaq, Vol>& primitive, const arma::imat33& supercell);

	private:
	// The primitive cell used after SNF decomposition.	
	PrimitiveCell<Point, Link, Plaq, Vol>& primitive_idx;

};


