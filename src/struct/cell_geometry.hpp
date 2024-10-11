#pragma once 

#include <array>
#include <cstddef>
#include <vector>
#include <cstdlib>
#include <smithNormalForm.hpp>
#include <cassert>


#include "chain.hpp"
#include "modulus.hpp"
#include "vec3.hpp"
#include "UnitCellSpecifier.hpp"


 
namespace CellGeometry {



/////////////////////////////////////////////
/////////////////////////////////////////////
///  FORWARD DECLARATIONS
/////////////////////////////////////////////
/////////////////////////////////////////////


typedef ipos_t idx3_t;


// Does the main part of the 3d indexing work
struct PeriodicAbstractLattice {
	PeriodicAbstractLattice(
			const UnitCellSpecifier& specified_primitive,
			const imat33_t& supercell
			) : 
	// Store the supercell's lattice vectors here by calling the base class constructor
	cell_vectors(specified_primitive.lattice_vectors * supercell),
	// Smith decopose the supercell spec to find a primitive cell that aligns nicely with the supercell
	LDW(ComputeSmithNormalForm( to_snfmat(supercell))),
	// Store the new primitve cell
	primitive_spec( specified_primitive,  LDW.Linv ),
	num_primitive(LDW.D[0]*LDW.D[1]*LDW.D[2])
	{}

	// 3-vectors, arranged columnwise, corresponding to supercell lengths 
	const imat33_t cell_vectors;

	// Size of the supercell in units of primitive cells
	inline ipos_t size() const {return LDW.D;} 
	inline int size(int idx) const {
		assert(idx >= 0 && idx < 3);
		return LDW.D[idx];
	}

	idx3_t get_supercell_IDX(ipos_t&R);

	///////////////////////////////////////////////////////
protected:
	// The Smith decompositions of the basis vectors, for indexing purposes
	const SNF_decomp LDW;
	// The primitive cell used after SNF decomposition
	const UnitCellSpecifier primitive_spec;
	//	const ArmaSmithDecomposition UPV;
	inline size_t idx_from_idx3(const idx3_t&I){
		return (I[2]*LDW.D[1] + I[1])*LDW.D[0] + I[0];
	}
public:
	const int num_primitive;
};


/**
 * Wraps r to primitive cell, and returns the primitive-cell index
 *
 * Given a 3D unit cell, seek an index I in [0,D_0) x [0, D_2) x [0, D_3)
 * such that R = b * I + A N + r
 * for some int vector N we do not care about
 * where the columns of A are the cell_vectors, 
 * cols of b are the primitive_lattice vectors, and 
 * r is within the primitive cell specified by b.
 *
 * Letting a be the specified cell_vectors, constructor set up
 * A = a Z = a L^-1 D W ^-1
 * b = a L^-1
 *
 * R = b * (I + D W^-1 N) + r
 */
inline idx3_t PeriodicAbstractLattice::get_supercell_IDX(ipos_t& R)
{	
	ipos_t UR = primitive_spec.UPV.L * R;
	auto [quot, rem] = moddiv(UR, primitive_spec.UPV.D);
	quot = primitive_spec.UPV.R  * quot;
	R =  primitive_spec.UPV.Linv * rem;
	// quot = V * ( UR // diag[P] ) 
	// R  = Uinv * mod( UR, diag[P] )
	// the sl position
	
	// wrap quot around supercell to  yield the index
	return mod(quot, LDW.D);
}


/* for reference:
 *  specified_unitcell:
 *      primitive_vectors <- a
 *      U, P, V s.t. U a V = P (ignore me)
 *
 *  primitive_vectors:
 *		primitive_vectors <- b
 *		U, P, V s.t. U b V = P
 *
 *		
 */   


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////// POINTS

template<
	CellLike<0> Point 
	>
	requires (Point::order == 0)
struct PeriodicPointLattice : public PeriodicAbstractLattice {
	PeriodicPointLattice(
			const UnitCellSpecifier& specified_primitive,
			const imat33_t& supercell
			) : PeriodicAbstractLattice(specified_primitive, supercell)
	{
		initialise_points();
	}
	// Object access
	// For all below:
	// R -> a realspace position in the units originally provided
	// I -> a vector in [0,L) (enable bounds check with #define DEBUG)
	// sl -> a sublattice index
	//
	inline Point& get_point_at(const ipos_t &R){
		ipos_t r(R);
		const idx3_t& I = get_supercell_IDX(r); // r now contains the sublattice index
		sl_t sl=primitive_spec.sl_of_point(r);
		return get_point_at(I, sl);
	}
	inline Point& get_point_at(const idx3_t& I, sl_t sl){
		return points[idx_from_idx3(I) + sl*num_primitive];
	}
	// const accessors	
	inline const Point& get_point_at(const ipos_t& R) const { return get_point_at(R); }
	inline const Point& get_point_at(const idx3_t& I, sl_t sl) const { return get_point_at(I, sl);}
	// Contains the 'point' geometric objects
	std::vector<Point> points;
private:
	void initialise_points(){	
		idx3_t IDX = {0,0,0};
		// ensure all have the right number of spaces
		this->points.resize(this->num_primitive*this->primitive_spec.num_point_sl());
	
		for (IDX[0]=0; IDX[0]<this->size(0); IDX[0]++){
		for (IDX[1]=0; IDX[1]<this->size(1); IDX[1]++){
		for (IDX[2]=0; IDX[2]<this->size(2); IDX[2]++){
			for (sl_t sl=0; sl<(sl_t) this->primitive_spec.num_point_sl(); sl++){
				const PointSpec& spec = primitive_spec.point_no(sl);
				Point& tmp = this->get_point_at(IDX, sl);
				tmp.position = spec.position;
				tmp.position += primitive_spec.lattice_vectors * IDX;
			}
		}}}	
	}
};

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////// LINKS
///
template<
	CellLike<0> Point,
	CellLike<1> Link
	>
struct PeriodicLinkLattice : public PeriodicPointLattice<Point>
{
	PeriodicLinkLattice(
			const UnitCellSpecifier& primitive,
			const imat33_t& supercell
			) :
		PeriodicPointLattice<Point>(primitive, supercell)
	{
		initialise_links();
		connect_boundaries();
	}

	// Object access
	// For all below:
	// R -> a realspace position in the units originally provided
	// I -> a vector in [0,L) (enable bounds check with #define DEBUG)
	// sl -> a sublattice index
	//
	inline Link& get_link_at(const ipos_t &R){
		ipos_t r(R);
		const idx3_t& I = this->get_supercell_IDX(r); // r now contains sl index
		sl_t sl=this->primitive_spec.sl_of_link(r);
		return get_link_at(I, sl);
	}
	inline Link& get_link_at(const idx3_t& I, sl_t sl){
		return links[this->idx_from_idx3(I) + sl*this->num_primitive];
	}
	inline const Link& get_link_at(const ipos_t& R) const { return get_link_at(R); }
	inline const Link& get_link_at(const idx3_t& I, sl_t sl) const { return get_link_at(I, sl);}
	std::vector<Link> links;

	private:
	void initialise_links(){
		idx3_t IDX = {0,0,0};
		// ensure all have the right number of spaces
		this->links.resize(this->num_primitive*this->primitive_spec.num_link_sl());
	
		for (IDX[0]=0; IDX[0]<this->size(0); IDX[0]++){
		for (IDX[1]=0; IDX[1]<this->size(1); IDX[1]++){
		for (IDX[2]=0; IDX[2]<this->size(2); IDX[2]++){
			for (int sl=0; sl<(sl_t)this->primitive_spec.num_link_sl(); sl++){
				const LinkSpec& spec = this->primitive_spec.link_no(sl);
				Link& tmp = this->get_link_at(IDX, sl);
				tmp.position = spec.position;
				tmp.position += this->primitive_spec.lattice_vectors * IDX;
			}
		}}}
	}
	
	void connect_boundaries(){ 
		// Stitch together the boundaries and coboundaries
		for (Link& l : links){
			// Iterate over link sublattices
			sl_t sl = this->primitive_spec.sl_of_link(l.position);
			const auto& linkspec = this->primitive_spec.link_no(sl); // link reference
			for (auto& bp : linkspec.boundary) {
				auto point = this->get_point_at(linkspec.position + bp.relative_position);
				l.boundary[&point] = bp.multiplier;
				point.coboundary[&l] = bp.multiplier;
			}
		}
	}

};



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////// PLAQS
///
template<
	CellLike<0> Point,
	CellLike<1> Link,
	CellLike<2> Plaq
	>
struct PeriodicPlaqLattice : public PeriodicLinkLattice<Point,Link>
{
	PeriodicPlaqLattice(
			const UnitCellSpecifier& primitive,
			const imat33_t& supercell
			) :
		PeriodicLinkLattice<Point,Link>(primitive, supercell)
	{
		initialise_plaqs();
		connect_boundaries();
	}

	// Object access
	// For all below:
	// R -> a realspace position in the units originally provided
	// I -> a vector in [0,L) (enable bounds check with #define DEBUG)
	// sl -> a sublattice index
	//
	inline Plaq& get_plaq_at(const ipos_t &R){
		ipos_t r(R);
		const idx3_t& I = this->get_supercell_IDX(r); // r now contains sl index
		sl_t sl=this->primitive_spec.sl_of_plaq(r);
		return get_plaq_at(I, sl);
	}
	inline Plaq& get_plaq_at(const idx3_t& I, sl_t sl){
		return plaqs[this->idx_from_idx3(I) + sl*this->num_primitive];
	}

	// const accessors	
	inline const Plaq& get_plaq_at(const ipos_t& R) const { return get_plaq_at(R); }
	inline const Plaq& get_plaq_at(const idx3_t& I, sl_t sl) const { return get_plaq_at(I, sl);}

	std::vector<Plaq> plaqs;

	private:
	void initialise_plaqs(){
		idx3_t IDX = {0,0,0};
		// ensure all have the right number of spaces
		this->plaqs.resize(this->num_primitive*this->primitive_spec.num_plaq_sl());
	
		for (IDX[0]=0; IDX[0]<this->size(0); IDX[0]++){
		for (IDX[1]=0; IDX[1]<this->size(1); IDX[1]++){
		for (IDX[2]=0; IDX[2]<this->size(2); IDX[2]++){
			for (int sl=0; sl<(sl_t)this->primitive_spec.num_plaq_sl(); sl++){
				const PlaqSpec& spec = this->primitive_spec.plaq_no(sl);
				Plaq& tmp = this->get_plaq_at(IDX, sl);
				tmp.position = spec.position;
				tmp.position += this->primitive_spec.lattice_vectors * IDX;
			}
		}}}
	}
	
	void connect_boundaries(){ 
		// Stitch together the boundaries and coboundaries
		for (Plaq& l : plaqs){
			// Iterate over plaq sublattices
			int sl = this->primitive_spec.sl_of_plaq(l.position);
			const auto& plaqspec = this->primitive_spec.plaq_no(sl); // plaq reference
			for (auto& bp : plaqspec.boundary) {
				auto link = this->get_link_at(plaqspec.position + bp.relative_position);
				l.boundary[&link] = bp.multiplier;
				link.coboundary[&l] = bp.multiplier;
			}
		}
	}
};


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////// VOLUMES
///
template<
	CellLike<0> Point,
	CellLike<1> Link,
	CellLike<2> Plaq,
	CellLike<3> Vol
	>
struct PeriodicVolLattice : public PeriodicPlaqLattice<Point,Link,Plaq>
{
	PeriodicVolLattice(
			const UnitCellSpecifier& primitive,
			const imat33_t& supercell
			) :
		PeriodicPlaqLattice<Point,Link,Plaq>(primitive, supercell)
	{
		initialise_vols();
		connect_boundaries();
	}

	// Object access
	// For all below:
	// R -> a realspace position in the units originally provided
	// I -> a vector in [0,L) (enable bounds check with #define DEBUG)
	// sl -> a sublattice index
	//
	inline Vol& get_vol_at(const ipos_t &R){
		ipos_t r(R);
		const idx3_t& I = this->get_supercell(r); // r now contains sl index
		sl_t sl=this->primitive_spec.sl_of_vol(r);
		return get_vol_at(I, sl);
	}
	inline Vol& get_vol_at(const idx3_t& I, sl_t sl){
		return &vols[this->idx_from_idx3(I) + sl*this->num_primitive()];
	}

	// const accessors	
	inline const Vol& get_vol_at(const ipos_t& R) const { return get_vol_at(R); }
	inline const Vol& get_vol_at(const idx3_t& I, sl_t sl) const { return get_vol_at(I, sl);}

	std::vector<Vol> vols;

	private:
	void initialise_vols(){
		idx3_t IDX = {0,0,0};
		// ensure all have the right number of spaces
		this->vols.resize(this->num_primitive()*this->primitive_spec.num_vol_sl());
	
		for (IDX[0]=0; IDX[0]<this->size(0); IDX[0]++){
		for (IDX[1]=0; IDX[1]<this->size(1); IDX[1]++){
		for (IDX[2]=0; IDX[2]<this->size(2); IDX[2]++){
			for (int sl=0; sl<this->primitive_spec.num_vol_sl(); sl++){
				const VolSpec& spec = this->primitive_spec.vol_no(sl);
				Vol* tmp = this->get_vol_at(IDX, sl);
				tmp->position = spec.position;
				tmp->position += this->primitive_spec.lattice_vectors * IDX;
			}
		}}}
	}
	
	void connect_boundaries(){ 
		// Stitch together the boundaries and coboundaries
		for (Vol& l : vols){
			// Iterate over vol sublattices
			int sl = sl_of_vol(l.position);
			const auto& volspec = this->primitive_spec.vol_no(sl); // vol reference
			for (auto& bp : volspec.boundary) {
				auto plaq = this->get_plaq_at(volspec.position + bp.relative_position);
				l.boundary[plaq] = bp.multiplier;
				plaq.coboundary[&l] = bp.multiplier;
			}
		}
	}
};





}; // End of namespace CellGeometry
