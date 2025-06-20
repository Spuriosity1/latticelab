#pragma once 

#include <cstddef>
#include <memory>
#include <unordered_map>
#include <vector>
#include <cstdlib>
#include <smithNormalForm.hpp>
#include <cassert>

#include "chain.hpp"
#include "modulus.hpp"
#include "vec3.hpp"
#include "UnitCellSpecifier.hpp"
#include "SortedVectorMap.hpp"


 
namespace CellGeometry {



/////////////////////////////////////////////
/////////////////////////////////////////////
///  FORWARD DECLARATIONS
/////////////////////////////////////////////
/////////////////////////////////////////////


typedef ivec3_t idx3_t;

template<class Key, class Tp>
using SparseMap = std::unordered_map<Key, Tp>;
// using SparseMap = SortedVectorMap<Key, Tp>;
//using SparseMap = FilteredVector<Key, Tp>;

// Does the main part of the 3d indexing work
// Represents a periodic region of space with nothing filling it
struct PeriodicAbstractLattice {
	PeriodicAbstractLattice(
			const UnitCellSpecifier& specified_primitive,
			const imat33_t& supercell
			) : 
	// Smith decopose the supercell spec to find a primitive cell that aligns 
	// nicely with the supercell
	LDW(ComputeSmithNormalForm( to_snfmat(supercell))),
	// Store the reparameterised supercell
	cell_vectors(specified_primitive.latvecs 
			* imat33_t::from_other(supercell)),
	// Cell vectors only used for indexing
	index_cell_vectors(specified_primitive.latvecs 
			* imat33_t::from_other(supercell * LDW.R)),
	num_primitive(LDW.D[0]*LDW.D[1]*LDW.D[2]),
	// Store the new primitve cell
	primitive_spec( specified_primitive,  LDW.Linv )
	{
	}

	// Size of the supercell in units of primitive cells
	inline ivec3_t size() const {return LDW.D;} 
	inline int size(int idx) const {
		assert(idx >= 0 && idx < 3);
		return LDW.D[idx];
	}

	// Converts a position R of a supercell point into a three-tuple lying in
	// [0, D[0]) x [0, D[1]) x [0, D[2]) \subset Z^3
	// Modifies its argument, leaving remainder there
	idx3_t get_supercell_IDX(ipos_t&R);



protected:
	// The Smith decompositions of the supercell spec Z, for indexing purposes
	const SNF_decomp LDW;
	inline size_t idx_from_idx3(const idx3_t&I){
		return (I[2]*LDW.D[1] + I[1])*LDW.D[0] + I[0];
	}
public:
	///////////////////////////////////////////////////////
	// 3-vectors, arranged columnwise, corresponding to supercell lengths 
	// i.e. j'th vector is cell_vectors[:, j]
	// a0 b0 c0
	// a1 b1 c1
	// a2 b2 c2
	const imat33_t cell_vectors; // = specified primitive * supercell
	const imat33_t index_cell_vectors;

	// The number of primitive cells
	const int num_primitive;

	// The primitive cell used after SNF decomposition
	const UnitCellSpecifier primitive_spec;
	//const rational::rmat33 primitive_cell_vectors;
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
//inline idx3_t PeriodicAbstractLattice::get_supercell_IDX(ipos_t& R)
//{	
//	ipos_t UR = primitive_spec.UPV.L * R;
//	auto [quot, rem] = moddiv(UR, primitive_spec.UPV.D);
//	quot = primitive_spec.UPV.R  * quot;
//	R =  primitive_spec.UPV.Linv * rem;
//	// quot = V * ( UR // diag[P] ) 
//	// R  = Uinv * mod( UR, diag[P] )
//	// the sl position
//	
//	// wrap quot around supercell to  yield the index
//	return mod(quot, LDW.D);
//}



/*
 * Wraps r to primitive cell, and returns the primitive-cell index
 *
 * Given a 3D unit cell, seek an index I in [0,D_0) x [0, D_1) x [0, D_2)
 * such that R = b * (I + D N) + r
 * for some N in Z3
 * where b is primitive_spec.lattice_vectors
 * mutating R to now contain the remainder r
*/
inline idx3_t PeriodicAbstractLattice::get_supercell_IDX(ipos_t& R) {
	// b^-1 R  = I + D N + b^-1 r
	ipos_t x = this->primitive_spec.latvecs_unnormed_inverse * R;
	idx3_t I;
	for (int n=0; n<3; n++){
		auto res = moddiv(x[n], primitive_spec.abs_det_latvecs);
		x[n] = res.rem;
		I[n] = res.quot;
		I[n] = mod(I[n], LDW.D[n]);
	}
	R = this->primitive_spec.latvecs * x;
	for (int n=0; n<3; n++){
		R[n] /= primitive_spec.abs_det_latvecs;
	}
	return I;
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


template<typename T>
void display_boundary(T boundary, unsigned verbosity){
	if (verbosity < 2) return;
	for (const auto& [p, _] : boundary) {
		std::cout << "d   boundary cell at " << p << "; coboundary size " << p->coboundary.size() << std::endl;
		if (verbosity >= 3){
			for (const auto& [link_ptr, coeff] : p->coboundary) {
				std::cout << "Dd      at " << link_ptr << " coeff " << coeff << std::endl;
			}
		}
	}
}

template <typename T>
void display_coboundary(T coboundary, unsigned verbosity){
	if (verbosity < 2) return;
	for (const auto& [p, _] : coboundary) {
		std::cout << "d   coboundary Plaq at " << p << "; boundary size " << p->boundary.size() << std::endl;
		if (verbosity >= 3){
			for (const auto& [link_ptr, coeff] : p->boundary) {
				std::cout << "Dd      at " << link_ptr << " coeff " << coeff << std::endl;
			}
		}
	}
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
///////// POINTS

template<
	CellLike<0> Point 
	>
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
		return *points.at(get_point_idx_at(R));
	}

	// const accessors	
	inline const Point& get_point_at(const ipos_t& R) const { return get_point_at(R); }
//	inline const Point& get_point_at(const idx3_t& I, sl_t sl) const { return get_point_at(I, sl);}

	// Tests if point exists in the map (slow; in principle can do in log time if we know spin* are sorted)
	bool has_point(const Point* point_it) const {
		for (const auto& [_, p] : points){
			if (p == point_it) return true;
		}
		return false;
	}

	// Deletes a point and all references to it
	void erase_point(Point* point_it){
		points.erase(this->get_point_idx_at(point_it->position));
		delete point_it;
	}

	void print_state(unsigned verbosity=3){
		if (verbosity == 0) {
			std::cout<< points.size() << " points" << std::endl;
			return;
		}
		for (auto [R, x] : points){
			std::cout<< "Point " << x <<" at "<<x->position<<std::endl;
			display_coboundary(x->coboundary, verbosity);
		}
	}

/*
    auto get_points() {
        return points 
             | std::views::values 
             | std::views::transform([](Point* p) -> Point& { return *p; });
    }

    auto get_points() const {
        return points 
             | std::views::values 
             | std::views::transform([](Point* p) -> Point const& { return *p; });
    } 
	*/

	// Contains the 'point' geometric objects
	SparseMap<sl_t, Point*> points;

private:

	inline sl_t get_point_idx_at(const ipos_t& R){	
		ipos_t r(R);
		const idx3_t& I = get_supercell_IDX(r); // r now contains the sublattice index
		sl_t sl = this->primitive_spec.sl_of_point(r);
		return this->idx_from_idx3(I)  + sl * this->num_primitive;
	}

	void initialise_points(){	
		idx3_t IDX = {0,0,0};
		// Allocate memory for all of the points we want
	
		for (IDX[0]=0; IDX[0]<this->size(0); IDX[0]++){
		for (IDX[1]=0; IDX[1]<this->size(1); IDX[1]++){
		for (IDX[2]=0; IDX[2]<this->size(2); IDX[2]++){
			for (sl_t sl=0; sl< this->primitive_spec.num_point_sl(); sl++){
				const PointSpec& spec = primitive_spec.point_no(sl);
				Point* tmp = new Point();
				tmp->position = spec.position + primitive_spec.latvecs * IDX;
				points[get_point_idx_at(tmp->position)] = tmp;
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
		connect_link_boundaries();
	}

	// Object access
	// For all below:
	// R -> a realspace position in the units originally provided
	// I -> a vector in [0,L) (enable bounds check with #define DEBUG)
	// sl -> a sublattice index
	//
	inline Link& get_link_at(const ipos_t &R){
		return *links.at(get_link_idx_at(R));
	}

	inline const Link& get_link_at(const ipos_t& R) const { return get_link_at(R); }

	// Deletes a link (and erases corresponding coboundary terms in point)
	void erase_link(Link* link_ptr){
		for (auto [p, _] : link_ptr->boundary){
			p->coboundary.erase(link_ptr);
			// silently fails if link_it not in the coboundary
		}
		// remove from the index
		links.erase(get_link_idx_at(link_ptr->position));
		delete link_ptr;
	}


	// Tests if link exists in the map (slow; in principle can do in log time if we know spin* are sorted)
	bool has_link(const Link* link_it) const {
		for (const auto& [_, p] : links){
			if (p == link_it) return true;
		}
		return false;
	}

	// Deletes a point (and connected links)
	void erase_point(Point* point_ptr) {
		// remove connected links
		std::vector<Link*> to_purge;
		for (auto [link_ptr, _] : point_ptr->coboundary){
			to_purge.push_back(static_cast<Link*>(link_ptr));
		}
		// note that erase_link modifies the point coboundary,
		// invalidating the iterator
		for (auto& link_ptr : to_purge) {
			erase_link( link_ptr );
		}
		PeriodicPointLattice<Point>::erase_point(point_ptr);
	}

	SparseMap<sl_t, Link*> links;

	void print_state(unsigned verbosity=3){
		PeriodicPointLattice<Point>::print_state(verbosity);
		if (verbosity == 0) {
			std::cout<< links.size() << " links" << std::endl;
			return;
		}
		for (auto [R, l] : links){
			std::cout<< "Link " << l <<" at "<<l->position<<std::endl;

			display_boundary(l->boundary, verbosity);
			display_coboundary(l->coboundary, verbosity);
			
		}
	}

private:

	inline sl_t get_link_idx_at(const ipos_t& R){	
		ipos_t r(R);
		const idx3_t& I = this->get_supercell_IDX(r); // r now contains the sublattice index
		sl_t sl = this->primitive_spec.sl_of_link(r);
		return this->idx_from_idx3(I)  + sl * this->num_primitive;
	}


	void initialise_links(){
		idx3_t IDX = {0,0,0};
		// Aloocate memory for the links
		// Place all of the links	
		for (IDX[0]=0; IDX[0]<this->size(0); IDX[0]++){
		for (IDX[1]=0; IDX[1]<this->size(1); IDX[1]++){
		for (IDX[2]=0; IDX[2]<this->size(2); IDX[2]++){
			for (sl_t sl=0; sl<this->primitive_spec.num_link_sl(); sl++){
				const LinkSpec& spec = this->primitive_spec.link_no(sl);
				auto tmp = new Link();
				tmp->position = spec.position + this->primitive_spec.latvecs * IDX;
				links[get_link_idx_at(tmp->position)] = tmp;
			}
		}}}
	}
	
	void connect_link_boundaries(){ 
		// Stitch together the boundaries and coboundaries
		for (auto [R, l] : links){
			// Iterate over link sublattices
			sl_t sl = this->primitive_spec.sl_of_link(l->position);
			const auto& linkspec = this->primitive_spec.link_no(sl); // link reference
			for (const auto& bp : linkspec.boundary) {
				auto& point = this->get_point_at(l->position + bp.relative_position);
				l->boundary[std::addressof(point)] = bp.multiplier;
				point.coboundary[l] = bp.multiplier;
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
		connect_plaq_boundaries();
	}

	// Object access
	// For all below:
	// R -> a realspace position in the units originally provided
	// I -> a vector in [0,L) (enable bounds check with #define DEBUG)
	// sl -> a sublattice index
	//
	inline Plaq& get_plaq_at(const ipos_t &R){
		return *plaqs.at(get_plaq_idx_at(R));
	}

	// const accessors	
	inline const Plaq& get_plaq_at(const ipos_t& R) const { return get_plaq_at(R); }


	// Tests if plaq exists in the map (slow; in principle can do in log time if we know spin* are sorted)
	bool has_plaq(const Plaq* plaq_it) const {
		for (const auto& [_, p] : plaqs){
			if (p == plaq_it) return true;
		}
		return false;
	}

	//Deletes a plaquette
	void erase_plaq(Plaq* plaq_ptr){
		// remove coreferences
		for (auto [l, _] : plaq_ptr->boundary){
			l->coboundary.erase(plaq_ptr);
			// silently fails if plaq_ptr not in the coboundary
		}
		// remove from index
		plaqs.erase(get_plaq_idx_at(plaq_ptr->position));
		delete plaq_ptr;
	}

	// Deletes a link (and associated points, plaqs...)
	void erase_link(Link* link_ptr){
		// remove connected plaqs 
		// note that erase_plaq invalidates the coboundary iterator, 
		// so we must iterate then store
		std::vector<Plaq*> to_purge;
		for (auto [plaq_ptr, _]: link_ptr->coboundary){
			to_purge.push_back(static_cast<Plaq*>(plaq_ptr));
		}
		for (auto plaq_ptr : to_purge){
			erase_plaq(plaq_ptr);
		}
		PeriodicLinkLattice<Point, Link>::erase_link(link_ptr);
	}

	// cascades up - deletes all connected plaqs too
	void erase_point(const Point* point_ptr){
		// remove connected links
		std::vector<Link*> to_purge; 
		for (auto [link_ptr, _] : point_ptr->coboundary){
			to_purge.push_back(static_cast<Link*>(link_ptr));
		}
		for (auto link_ptr : to_purge){
			erase_link(link_ptr);
		}
		PeriodicPointLattice<Point>::erase_point(point_ptr);
	}


	SparseMap<sl_t, Plaq*> plaqs;


	void print_state(unsigned verbosity =3){
		PeriodicLinkLattice<Point, Link>::print_state(verbosity);	
		if (verbosity == 0) {
			std::cout<< plaqs.size() << " plaqs" << std::endl;
			return;
		}

		for (auto [R, pl] : plaqs){
			std::cout<< "Plaq " << pl <<" at " <<pl->position<<std::endl;
			display_boundary(pl->boundary, verbosity);
			display_coboundary(pl->coboundary, verbosity);
		}

	}

private:
	inline sl_t get_plaq_idx_at(const ipos_t& R){	
		ipos_t r(R);
		const idx3_t& I = this->get_supercell_IDX(r); // r now contains the sublattice index
		sl_t sl = this->primitive_spec.sl_of_plaq(r);
		return this->idx_from_idx3(I)  + sl * this->num_primitive;
	}

	void initialise_plaqs(){
		idx3_t IDX = {0,0,0};
		// ensure all have the right number of spaces
	
		for (IDX[0]=0; IDX[0]<this->size(0); IDX[0]++){
		for (IDX[1]=0; IDX[1]<this->size(1); IDX[1]++){
		for (IDX[2]=0; IDX[2]<this->size(2); IDX[2]++){
			for (sl_t sl=0; sl<this->primitive_spec.num_plaq_sl(); sl++){
				const PlaqSpec& spec = this->primitive_spec.plaq_no(sl);
				auto tmp = new Plaq();
				tmp->position = spec.position + this->primitive_spec.latvecs * IDX;
				plaqs[get_plaq_idx_at(tmp->position)] = tmp;
			}
		}}}
	}
	
	void connect_plaq_boundaries(){ 
		// Stitch together the boundaries and coboundaries
		for (auto [R, pl] : plaqs){
			// Iterate over plaq sublattices
			int sl = this->primitive_spec.sl_of_plaq(pl->position);
			const auto& plaqspec = this->primitive_spec.plaq_no(sl); // plaq reference
			for (const auto& bp : plaqspec.boundary) {
				auto& link = this->get_link_at(pl->position + bp.relative_position);
				pl->boundary[std::addressof(link)] = bp.multiplier;
				link.coboundary[pl] = bp.multiplier;
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
		connect_vol_boundaries();
	}

	// Object access
	// For all below:
	// R -> a realspace position in the units originally provided
	// I -> a vector in [0,L) (enable bounds check with #define DEBUG)
	// sl -> a sublattice index
	//
	inline Vol& get_vol_at(const ipos_t &R){
		return *vols.at(get_vol_idx_at(R));
	}

	// const accessors	
	inline const Vol& get_vol_at(const ipos_t& R) const { return get_vol_at(R); }


	// Tests if vol exists in the map (slow; in principle can do in log time if we know spin* are sorted)
	bool has_vol(const Vol* vol_it) const {
		for (const auto& [_, p] : vols){
			if (p == vol_it) return true;
		}
		return false;
	}


	void erase_vol(Vol* vol_ptr){
		// remove coreferences
		for (auto [p, _] : vol_ptr->boundary){
			p->coboundary.erase(vol_ptr);
		}
		// remove from index
		vols.erase(get_vol_idx_at(vol_ptr->position));
		delete vol_ptr;
	}


	//Deletes a plaquette
	void erase_plaq(Plaq* plaq_ptr){
		// delete volumes
		std::vector<Vol*> to_purge;
		for (auto [v, _] : plaq_ptr->coboundary){
			to_purge.push_back(static_cast<Vol*>(v));
		}
		for (auto v : to_purge) { erase_vol(v); }

		PeriodicPlaqLattice<Point, Link, Plaq>::erase_plaq(plaq_ptr);
	}

	// Deletes a link (and associated points, plaqs...)
	void erase_link(Link* link_ptr){
		// remove connected plaqs (merging too weird to think about...)
		std::vector<Plaq*> to_purge;
		for (auto [plaq_ptr, _]: link_ptr->coboundary){
			to_purge.push_back(static_cast<Plaq*>(plaq_ptr));
		}
		for (auto plaq_ptr : to_purge){
			erase_plaq(plaq_ptr);
		}
		PeriodicLinkLattice<Point, Link>::erase_link(link_ptr);
	}

	// cascades up - deletes all connected plaqs too
	void erase_point(Point* point_ptr){
		// remove connected links

		std::vector<Link*> to_purge;
		for (auto [c1_ptr, _] : point_ptr->coboundary){
			to_purge.push_back(static_cast<Link*>(c1_ptr));
		}
		for (auto link_ptr: to_purge){
			erase_link(link_ptr);
		}
		PeriodicPointLattice<Point>::erase_point(point_ptr);
	}

	SparseMap<sl_t, Vol*> vols;


	void print_state(unsigned verbosity=3){
		PeriodicPlaqLattice<Point, Link, Plaq>::print_state(verbosity);
		if (verbosity == 0) {
			std::cout<< vols.size() << " vols" << std::endl;
			return;
		}
		for (auto [R, pl] : vols){
			std::cout<< "Vol " << pl <<" at " <<pl->position<<std::endl;
			display_boundary(pl->boundary, verbosity);
		}


	}


private:
	inline sl_t get_vol_idx_at(const ipos_t& R){
		ipos_t r(R);
		const idx3_t& I = this->get_supercell_IDX(r); // r now contains the sublattice index
		sl_t sl = this->primitive_spec.sl_of_vol(r);
		return this->idx_from_idx3(I)  + sl * this->num_primitive;
	}

	void initialise_vols(){
		idx3_t IDX = {0,0,0};
	
		for (IDX[0]=0; IDX[0]<this->size(0); IDX[0]++){
		for (IDX[1]=0; IDX[1]<this->size(1); IDX[1]++){
		for (IDX[2]=0; IDX[2]<this->size(2); IDX[2]++){
			for (sl_t sl=0; sl<this->primitive_spec.num_vol_sl(); sl++){
				const VolSpec& spec = this->primitive_spec.vol_no(sl);
				auto tmp = new Vol();
				tmp->position = spec.position + this->primitive_spec.latvecs * IDX;
				vols[get_vol_idx_at(tmp->position)] = tmp;
			} 
		}}}
	}

	void connect_vol_boundaries(){ 
		// Stitch together the boundaries and coboundaries
		for (auto [R, v]: vols){
			// Iterate over vol sublattices
			int sl = this->primitive_spec.sl_of_vol(v->position);
			const auto& volspec = this->primitive_spec.vol_no(sl); // vol reference
			for (const auto& bp : volspec.boundary) {
				auto& plaq = this->get_plaq_at(v->position + bp.relative_position);
				v->boundary[std::addressof(plaq)] = bp.multiplier;
				plaq.coboundary[v] = bp.multiplier;
			}
		}
	}
};



template <typename T>
requires (CellOfAnyOrder<T> && T::order > 0)
std::vector<T*> get_neighbours(T* x0){
    std::vector<T*> retval;
    for (auto& [pl, m] : x0->boundary){
        for (auto& [x1, _] : pl->coboundary){
            if (x1 != x0) retval.push_back(static_cast<T*>(x1));
        }
    }
    return retval;
}


template <typename T, int order>
requires (CellOfAnyOrder<T> && T::order < 3)
std::vector<T*> get_coneighbours(T* x0){
    std::vector<T*> retval;
    for (auto& [pl, m] : x0->coboundary){
        for (auto& [x1, _] : pl->boundary){
            if (x1 != x0) retval.push_back(static_cast<T*>(x1));
        }
    }
    return retval;
}




}; // End of namespace CellGeometry
