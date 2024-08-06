#pragma once

#include <smithNormalForm.hpp>
#include <concepts>
#include <stdexcept>
#include <unordered_map>

#include "chain.hpp"
#include "matrix.hpp"
#include "modulus.hpp"

namespace CellGeometry {


/////////////////////////////////////////////
/////////////////////////////////////////////
///  CONVENIENCE TYPEDEFS
/////////////////////////////////////////////
/////////////////////////////////////////////


typedef SmithNormalFormCalculator::Matrix<arma::sword> snfmat;

typedef  CellSpecifier<3> VolSpec   ;
typedef  CellSpecifier<2> PlaqSpec  ;
typedef  CellSpecifier<1> LinkSpec  ;
typedef  CellSpecifier<0> PointSpec ;

typedef int sl_t;

/////////////////////////////////////////////
/////////////////////////////////////////////
///  FORWARD DECLARATIONS
/////////////////////////////////////////////
/////////////////////////////////////////////
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
	SmithNormalFormCalculator::Matrix<T> out;
	for (int i=0; i<3; i++)
		for (int j=0; j<3; j++)
			out[i][j] = m(i,j);
	return out;
}

// TODO move to the subproject
struct ArmaSmithDecomposition {
	ArmaSmithDecomposition(
			SmithNormalFormCalculator::SmithNormalFormDecomposition<arma::sword>decomp ) : 
		L(from_snfmat(decomp.L)),
		Linv(from_snfmat(SmithNormalFormCalculator::inverse(decomp.L))),
		D(from_snfmat(decomp.D).diagonal()),
		R(from_snfmat(decomp.R)),
		Rinv(from_snfmat(SmithNormalFormCalculator::inverse(decomp.R)))
	{}

	const imat33_t L;
	const imat33_t Linv;
	const ipos_t D;
	const imat33_t R;
	const imat33_t Rinv;
};

/** A specialised "hashmap" implementation
 * exploiting the knowledge that operator[] takes in an integer vector of the 
 * form nx, ny, nz, entrywise 0 <= nx < Lx
 *
 */
struct pointMap {
	pointMap(const ipos_t& L, const size_t max_memory=2000):
	L(L),
	max_memory(max_memory)
	{
		assert(L[0] > 0);
		assert(L[1] > 0);
		assert(L[2] > 0);
		this->values.resize(L[0]*L[1]*L[2]);
		for (auto& v : values){
			v = -1;
		}
		if ( static_cast<size_t>(L[0]*L[1]*L[2]) > max_memory ) {
			throw std::out_of_range("Unit cell is to large");
		}
	}

	/**
	 * Attempts to insert sublattice 'value' at position 'key'.
	 * Returns 'true' if the value was overwritten
	 */
	bool insert(ipos_t& key, sl_t value){
		assert_inbounds(key);
		sl_t& v = values[(key[2]*L[1]+key[1])*L[0]+key[0]];
		if (v == -1){
			v = value;
			return true;
		} else {
			// Something already there
			return false;
		}
	}

	sl_t operator[](const ipos_t& key){
#ifdef DEBUG
		assert_inbounds(key);
		if (values[(key[2]*L[1]+key[1])*L[0]+key[0]] == -1){
			throw std::out_of_range("Sublattice is not initialised");
		}
#endif
		return values[(key[2]*L[1]+key[1])*L[0]+key[0]];
	}

	const ipos_t L;
	const size_t max_memory;

	private:
/*	TODO reindex everything if there is a common factor of 2^n
 *	to save memory
 *	auto gcd(const ipos_t& x){
		return std::gcd(std::gcd(x[0],x[1]),x[2]);
	}*/

	void assert_inbounds(const ipos_t& x) const {
		if (x[0] < 0 || x[0] >= L[0] || x[1] < 0 
				|| x[1] >= L[1] || x[2] < 0 || x[2]>= L[2]){
			std::stringstream s;
			s << "Point x = " << x << " is outside of the SL index scheme";
			s << "(Expected 0 <= x < " << L << " )\n";
			throw std::out_of_range(s.str());
		}
	}

	std::vector<sl_t> values;
};

struct UnitCellSpecifier {	
	UnitCellSpecifier(const imat33_t& lattice_vectors_);
	UnitCellSpecifier(const snfmat& lattice_vectors_);/**
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
	const imat33_t lattice_vectors;
	// Smith decomposition of lattice_vectors
	const ArmaSmithDecomposition UPV;

	// Creation methods
	// Use these over manipulating the arrays directly
	void add_point(const PointSpec& p);
	void add_link(const LinkSpec& p);
	void add_plaq(const PlaqSpec& p);
	void add_vol(const VolSpec& v);

	// Wraps the vector X in-place to within this unit cell.
	void wrap(ipos_t& X) const;
	void wrap(GeometricObject& X) const;
	ipos_t wrap_copy(const ipos_t& X) const { ipos_t Y(X); wrap(Y); return Y; }

	// Const access methods 
	// (with bounds checking since this shouldn't be performance critical)
	const PointSpec& point_no(sl_t sl) const { return points.at(sl); }
	const LinkSpec&  link_no (sl_t sl) const { return links.at(sl); }
	const PlaqSpec&  plaq_no(sl_t sl) const { return plaqs.at(sl); }
	const VolSpec&   vol_no(sl_t sl) const { return vols.at(sl); }

	// Sublattice counters
	size_t num_point_sl() const { return points.size(); }
	size_t num_link_sl() const { return links.size(); }
	size_t num_plaq_sl() const { return plaqs.size(); }
	size_t num_vol_sl() const { return vols.size(); }

	// Sublattice index access from a position (no bounds check)
	sl_t get_sl_of_point(const ipos_t& R){return point_index[R];}
	sl_t get_sl_of_link(const ipos_t& R){return link_index[R];}
	sl_t get_sl_of_plaq(const ipos_t& R){return plaq_index[R];}
	sl_t get_sl_of_vol(const ipos_t& R){return vol_index[R];}

	// Checking correctness
	bool is_point(const ipos_t& R){return point_index[R] != -1;}
	bool is_link(const ipos_t& R){return link_index[R] != -1;}
	bool is_plaq(const ipos_t& R){return plaq_index[R] != -1;}
	bool is_vol(const ipos_t& R){return vol_index[R] != -1;}

	protected:
	std::vector<PointSpec> points;
	std::vector<LinkSpec> links;
	std::vector<PlaqSpec> plaqs;
	std::vector<VolSpec> vols;

	//typedef std::unordered_map<ipos_t,sl_t,iposHasher,decltype(iposEqual)> point_idx_t;
	typedef pointMap point_idx_t;

	point_idx_t point_index;	
	point_idx_t link_index;	
	point_idx_t plaq_index;	
	point_idx_t vol_index;	


	void try_insert_position(
			UnitCellSpecifier::point_idx_t hmap,
			ipos_t x, sl_t sl,
			const char* obj_name);

	

};

}; // end of namespace
