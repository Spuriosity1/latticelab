#pragma once


#include "chain.hpp"
#include "matrix.hpp"
#include "smithNormalForm.hpp"
#include "rationalmath.hpp"

namespace CellGeometry {


/////////////////////////////////////////////
/////////////////////////////////////////////
///  CONVENIENCE TYPEDEFS
/////////////////////////////////////////////
/////////////////////////////////////////////


typedef SmithNormalFormCalculator::Matrix<long long int> snfmat;

typedef  CellSpecifier<0> PointSpec ;
typedef  CellSpecifier<1> LinkSpec  ;
typedef  CellSpecifier<2> PlaqSpec  ;
typedef  CellSpecifier<3> VolSpec   ;

typedef int sl_t;

/////////////////////////////////////////////
/////////////////////////////////////////////
///  FORWARD DECLARATIONS
/////////////////////////////////////////////
/////////////////////////////////////////////
///

// type converters
template<typename T>
inline vector3::mat33<T> from_snfmat(SmithNormalFormCalculator::Matrix<T> m){
	vector3::mat33<T> out;
	for (int i=0; i<3; i++)
		for (int j=0; j<3; j++)
			out(i,j) = m[i][j];
	return out;
}

template<typename T>
inline SmithNormalFormCalculator::Matrix<T> to_snfmat(vector3::mat33<T>m){
	SmithNormalFormCalculator::Matrix<T> out(3,3);
	for (int i=0; i<3; i++)
		for (int j=0; j<3; j++)
			out[i][j] = m(i,j);
	return out;
}

struct SNF_decomp {
	SNF_decomp(
			SmithNormalFormCalculator::SmithNormalFormDecomposition<long long int>decomp ) : 
		L(from_snfmat(decomp.L)),
		Linv(from_snfmat(SmithNormalFormCalculator::inverse(decomp.L))),
		D(from_snfmat(decomp.D).diagonal()),
		R(from_snfmat(decomp.R)),
		Rinv(from_snfmat(SmithNormalFormCalculator::inverse(decomp.R)))
	{}

	const imat33_t L;
	const imat33_t Linv;
	const ivec3_t D;
	const imat33_t R;
	const imat33_t Rinv;
};

struct UnitCellSpecifier {	
	UnitCellSpecifier(const rational::rmat33& lattice_vectors_);
	/**
	 * Creates a new unitcell from the points of the old one, with
	 * new primitive_vectors such that
	 * `this->primitive_vectors = other.primitive_vectors * cellspec`
	 *
	 *
	 * @param other -> another UnitCellSpecifier to build from
	 * @param cellspec -> the change of basis matrix
	 */
	UnitCellSpecifier(const UnitCellSpecifier& other, const imat33_t& cellspec);

	// The lattice vectors (columns are interpreted as vectors)
	// [ a1  a2  a3]
	const rational::rmat33 lattice_vectors;
	const rational::rmat33 lattice_vectors_inverse; 
	// Smith decomposition of lattice_vectors
	//const SNF_decomp UPV;

	// Creation methods
	// Use these over manipulating the arrays directly
	void add_point(const PointSpec& p);
	void add_link(const LinkSpec& p);
	void add_plaq(const PlaqSpec& p);
	void add_vol(const VolSpec& v);

	// Wraps the vector X in-place to within this unit cell, i.e. such that
	// lattice_vectors_inverse * X in [0,1)^3
	void wrap(ipos_t& X) const;
	void wrap(GeometricObject& X) const;
	ipos_t wrap_copy(const ipos_t& X) const {
		ipos_t Y(X);
		wrap(Y);
		return Y;
	}

	// Const access methods 
	// with bounds check (consider if needed)
	const PointSpec& point_no(sl_t sl) const { return points.at(sl); }
	const LinkSpec&  link_no (sl_t sl) const { return links.at(sl); }
	const PlaqSpec&  plaq_no(sl_t sl) const { return plaqs.at(sl); }
	const VolSpec&   vol_no(sl_t sl) const { return vols.at(sl); }

	// Sublattice counters
	sl_t num_point_sl() const { return points.size(); }
	sl_t num_link_sl() const { return links.size(); }
	sl_t num_plaq_sl() const { return plaqs.size(); }
	sl_t num_vol_sl() const { return vols.size(); }

	// Sublattice index access from a physical position (real units)
	// Linear search suboptimal here -- optimisation pointless, 
	// pointers to everything already stored.
	sl_t sl_of_point(const ipos_t& R) const;
	sl_t sl_of_link(const ipos_t& R) const;
	sl_t sl_of_plaq(const ipos_t& R) const;
	sl_t sl_of_vol(const ipos_t& R) const;

	bool is_point(const ipos_t& R);
	bool is_link(const ipos_t& R);
	bool is_plaq(const ipos_t& R);
	bool is_vol(const ipos_t& R);



protected:
	std::vector<PointSpec> points;
	std::vector<LinkSpec> links;
	std::vector<PlaqSpec> plaqs;
	std::vector<VolSpec> vols;

	bool is_valid_position(const ipos_t& R){
		return (R[0].denom != 0) && (R[1].denom != 0) && (R[2].denom != 0);
	}

	void assert_valid_position(const ipos_t& R){
#ifdef DEBUG
		if ( !is_valid_position(R) ) {
			throw std::invalid_argument("spec has invalid location");
		}
#endif
	}


	/*
	// An index of hashmap-like things. Thought of as a function, it is
	// index : [0,D[0]) x [0,D[1]) x [0,D[2]) -> sl_t
	// point_index[ {D[0], 0, 0} ] or similar will throw.
	// It is an index of U*R's.
	point_idx_t point_index;	
	point_idx_t link_index;	
	point_idx_t plaq_index;	
	point_idx_t vol_index;
	*/
/*
    // attempts to insert the position `x` into the pointMap hmap, throwing
    // exceptions if `x` is outside the index scheme.
	// @param obj_name only used to make the error message easier to read. 
    inline void insert_position(point_idx_t& hmap, ipos_t x, sl_t sl,
		const char* obj_name="object");

    // attempts to read the value at position `x` of pointMap hmap, 
	// no bounds check.
	// @param obj_name only used to make the error message easier to read. 
    sl_t sl_from_position(point_idx_t& hmap, ipos_t x);

	// Checks if the point 'x' in the index scheme is initialised in hmap.
	// If x is an invalid coordinate, an exception is thrown.
	bool is_initialised(point_idx_t& hmap,ipos_t x, const char* obj_name="object");
*/
};

}; // end of namespace
