#pragma once 

#include <armadillo>
#include <array>
#include <vector>
#include <cstdlib>
#include <smithNormalForm.hpp>
#include <concepts>
#include <type_traits>

#include "modulus.hpp"
#include "UnitCellSpecifier.hpp"


 
namespace CellGeometry {



/////////////////////////////////////////////
/////////////////////////////////////////////
///  FORWARD DECLARATIONS
/////////////////////////////////////////////
/////////////////////////////////////////////



/*
 * Represents a periodic torus of tiled primitive cells.
 * @tparam Point the class representing pointlike objects.
 * @tparam Link the class representing bondlike objects.
 * @tparam Plaq the class representing plaquette-like objects.
 * @tparam Vol the class representing volume-like objects.
 */
template<
	std::derived_from<Cell<0>> Point = Cell<0>, 
	std::derived_from<Cell<1>> Link  = Cell<1>,
	std::derived_from<Cell<2>> Plaq  = Cell<2>,
	std::derived_from<Cell<3>> Vol   = Cell<3>
	>
struct PeriodicLattice {
	/*
	 * @param primitive   a primitive UnitCellSpecifier to tile over space.
	 * @param supercell   Specifies the supercell transformation matrix, such
	 *  that the supercell translation vectors (i.e. the ones that becoem equivalent to 0)
	 *  are given by 
	 *  `supercell_vectors = primitive_vectors * supercell
	 */
	PeriodicLattice(
			const UnitCellSpecifier& primitive,
			const arma::imat33& supercell
			);

	// Const properties
	const arma::imat33 cell_vectors;

	// Object access
	Point* get_point_at(const ipos_t &R);
	Point* get_point_at(const arma::ivec3& I, int sl);

	Link* get_link_at(const ipos_t &R);
	Link* get_link_at(const arma::ivec3& I, int sl);

	Plaq* get_plaq_at(const ipos_t &R);
	Plaq* get_plaq_at(const arma::ivec3& I, int sl);

	Vol* get_vol_at(const ipos_t &R);
	Vol* get_vol_at(const arma::ivec3& I, int sl);

	// Size of the supercell
	arma::ivec3 size() const {return LDW.D;} 
	int size(int idx) const {
		assert(idx >= 0 && idx < 3);
		return LDW.D[idx];
	}

	int num_primitive() const {
		return LDW.D[0]*LDW.D[1]*LDW.D[2];
	}

	///////////////////////////////////////////////////////
	protected:
	std::vector<Point> points;
	std::vector<Link>  links;
	std::vector<Plaq>  plaqs;
	std::vector<Vol>   vols;

	// Populates all cells of the lattice object at their correct positions
	void populate_cells();
	// Links (n)-cells with their boundary (n-1)-cells
	void link_boundaries();
	// Links (n)-cells with their coboundary (n+1)-cells
	void link_coboundaries();

	// The primitive cell used after SNF decomposition.	
	const UnitCellSpecifier primitive_spec;
	// The Smith decompositions of the basis vectors, for indexing purposes
	const ArmaSmithDecomposition LDW;
//	const ArmaSmithDecomposition UPV;
};


/////////////////////////////////////////////
/////////////////////////////////////////////
///  IMPLEMENTATION
/////////////////////////////////////////////
/////////////////////////////////////////////

template<
	std::derived_from<Cell<0>> Point , 
	std::derived_from<Cell<1>> Link  ,
	std::derived_from<Cell<2>> Plaq  ,
	std::derived_from<Cell<3>> Vol   
	>
PeriodicLattice<Point, Link, Plaq, Vol>::PeriodicLattice(
	const UnitCellSpecifier& specified_primitive,
	const arma::imat33& supercell
) : 
	// Store the supercell's lattice vectors here by calling the base class constructor
	cell_vectors(specified_primitive.lattice_vectors * supercell),
	// Smith decopose the supercell spec to find a primitive cell that aligns nicely with the supercell
	LDW(ComputeSmithNormalForm( snfmat(supercell))),
	// Store the new primitve cell
	primitive_spec( specified_primitive,  LDW.Linv )
{


}


template<
	std::derived_from<Cell<0>> Point , 
	std::derived_from<Cell<1>> Link  ,
	std::derived_from<Cell<2>> Plaq  ,
	std::derived_from<Cell<3>> Vol   
	>
void PeriodicLattice<Point, Link, Plaq, Vol>::populate_cells(){ 
	// First Pass: Populate all cells
	arma::ivec3 IDX = {0,0,0};
	this->points.resize(num_primitive()*this->primitive_spec.points.size());
	this->links.resize(num_primitive()*this->primitive_spec.links.size());
	this->plaqs.resize(num_primitive()*this->primitive_spec.plaqs.size());
	this->vols.resize(num_primitive()*this->primitive_spec.vols.size());
	for (IDX[0]=0; IDX[0]<this->size(0); IDX[0]++){
	for (IDX[1]=1; IDX[1]<this->size(1); IDX[1]++){
	for (IDX[2]=2; IDX[2]<this->size(2); IDX[2]++){
		// Iterate over all of the cell specifiers
		// I am cure that some tempalte magic could shorten this,
		// but a) idk how and b) that's probably horrendous to read
		for (int sl=0; sl<this->primitive_spec.num_point_sl(); sl++){
			const PointSpec& spec = primitive_spec.point_no(sl);
			Point* tmp = this->get_point_at(IDX, sl);
			tmp->position = spec.position;
			tmp->position += primitive_spec.lattice_vectors * IDX;
		}

		for (int sl=0; sl<this->primitive_spec.num_link_sl(); sl++){
			const LinkSpec& spec = primitive_spec.link_no(sl);
			Link* tmp = this->get_link_at(IDX, sl);
			tmp->position = spec.position;
			tmp->position += primitive_spec.lattice_vectors * IDX;
		}

		for (int sl=0; sl<this->primitive_spec.num_link_sl(); sl++){
			const PlaqSpec& spec = primitive_spec.plaq_no(sl);
			Plaq* tmp = this->get_plaq_at(IDX, sl);
			tmp->position = spec.position;
			tmp->position += primitive_spec.lattice_vectors * IDX;
		}

		for (int sl=0; sl<this->primitive_spec.num_vol_sl(); sl++){
			const VolSpec& spec = primitive_spec.vol_no(sl);
			Vol* tmp = this->get_vol_at(IDX, sl);
			tmp->position = spec.position;
			tmp->position += primitive_spec.lattice_vectors * IDX;
		}
	}}}	
}

// NOTE:
// Index convention ->
//
// J =  sl*L2*L1*L0  + (i2*L1 + i1)*L0 +  i0
//
// Pros (+) and cons (-) ->
// + Easy vectorisation for adding sublattices
// + Convenient for calculating SL-dependent FT
// - physically nearby spins have very distant memory addresses (unavoidable)
// 

template<
	std::derived_from<Cell<0>> Point , 
	std::derived_from<Cell<1>> Link  ,
	std::derived_from<Cell<2>> Plaq  ,
	std::derived_from<Cell<3>> Vol   
	>
void PeriodicLattice<Point, Link, Plaq, Vol>::link_boundaries(){ 
	// Second Pass: Stitch together the boundaries		
	arma::ivec3 IDX = {0,0,0};
	for (IDX[0]=0; IDX[0]<this->size(0); IDX[0]++){
	for (IDX[1]=1; IDX[1]<this->size(1); IDX[1]++){
	for (IDX[2]=2; IDX[2]<this->size(2); IDX[2]++){
		// Iterate over link sublattices
		for (int sl=0; sl<primitive_spec.links.size(); sl++){
			const auto& linkspec = primitive_spec.links[sl]; // link reference
			for (auto& bp : linkspec.boundary) {
				CellMultPair<0> tmp;
				tmp.multiplier = bp.multiplier;
				tmp.cell = this->get_point_at(linkspec.position + bp.relative_position);
				this->link_at(IDX, sl)->boundary.push_back(tmp);
			}
		}

		// Iterate over link sublattices
		for (int sl=0; sl<primitive_spec.plaqs.size(); sl++){
			const auto& plaqspec = primitive_spec.links[sl]; // link reference
			for (auto& bp : plaqspec.boundary) {
				CellMultPair<1> tmp;
				tmp.multiplier = bp.multiplier;
				tmp.cell = this->get_point_at(plaqspec.position + bp.relative_position);
				this->plaq_at(IDX, sl)->boundary.push_back(tmp);
			}
		}
	}}}
}


// If you think I've lost it, you're correct. 
// I have this kind of an interface in mind:
/*
struct spin: public Cell<1> {
	double heis[3];

};

double plaq_ring(const Cell<2>& plaq){
	double real = 1;
	double imag = 0;
	double tmp;
	for (auto& [mul, link] : plaq.boundary){
		// real + i imag *= heis[0] + i mul * heis[1]
		tmp = real * heis[0] - imag*mul*heis[1];
		imag = real *mul*heis[1] + imag*heis[0];
		real = tmp;	
	}
	return retval.real();
}
*/


/*
# R = a[Z n + x] + r
# U R = P V^{-1} [D m + L x ] + U r
UR = (U * R)
Y = V *  fldiv.( UR, diag[P] ) # Floor divide to snap to correct unit cell
Lx = modulo.(Y, D) # Entrywise float valued modulo wraps all entries to [0, D[i])
# sublattice info is contained here
r = Uinv * modulo,(UR, diag[P]) 
*/
template<
	std::derived_from<Cell<0>> Point , 
	std::derived_from<Cell<1>> Link  ,
	std::derived_from<Cell<2>> Plaq  ,
	std::derived_from<Cell<3>> Vol   
	>
Point* PeriodicLattice<Point, Link, Plaq, Vol>::get_point_at(const ipos_t &R){
	arma::ivec3 UR = primitive_spec.UPV.L * R;
	auto [quot, rem] = moddiv(UR, primitive_spec.UPV.D);
	quot = primitive_spec.UPV.R  * quot;
	rem =  primitive_spec.UPV.Linv * rem;
	// quot = V * ( UR // diag[P] ) 
	// rem  = Uinv * mod( UR, diag[P] )
	//
	// wrap quot around to yield the index
	for (int i=0; i<3; i++){
		quot[i] = mod(quot[i], LDW.D[i]);
	}
	
	
}


}; // End of namespace CellGeometry
