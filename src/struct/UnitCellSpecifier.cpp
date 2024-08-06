#include "UnitCellSpecifier.hpp"
#include "chain.hpp"
#include <sstream>
#include <stdexcept>
#include <string>
#include <format>

namespace CellGeometry {

std::string pprint(ipos_t x){
	return std::format("[{} {} {}]",x[0],x[1],x[2]);
}

// utility
void UnitCellSpecifier::try_insert_position(
		point_idx_t hmap,
		ipos_t x, sl_t sl,
		const char* obj_name){
	bool status = hmap.insert(x,sl);
	if (status == false){
		std::stringstream s;
		s<<"Failed to add "<<obj_name<<" at "<<pprint(x)<<": ";
		s<<"There is already a(n) "<<obj_name<<" there\n";
		throw std::runtime_error(s.str());	
	}
}

template<int _order>
void throw_bad_boundary(CellSpecifier<_order> p, ipos_t missing_point) {
	std::stringstream s;
	s << "Problem with boundary of "<<_order<<"-cell specifier ";
	s << "at "<<p.position<<": ";
	s << "There is no "<<_order-1<<"-cell at " << missing_point;
	throw std::out_of_range(s.str());
}

// constructors
UnitCellSpecifier::UnitCellSpecifier(const imat33_t& lattice_vectors_) :
	lattice_vectors(lattice_vectors_),
	UPV(ComputeSmithNormalForm(to_snfmat(lattice_vectors))),
	point_index(UPV.D),link_index(UPV.D),plaq_index(UPV.D),vol_index(UPV.D)
{}

UnitCellSpecifier::UnitCellSpecifier(const snfmat& lattice_vectors_) :
	lattice_vectors(from_snfmat(lattice_vectors_)),
	UPV(ComputeSmithNormalForm(lattice_vectors_)),
	point_index(UPV.D),link_index(UPV.D),plaq_index(UPV.D),vol_index(UPV.D)
{}

UnitCellSpecifier::UnitCellSpecifier(
		const UnitCellSpecifier& other,
		const imat33_t& cellspec) : 
	lattice_vectors(other.lattice_vectors * cellspec),
	UPV(ComputeSmithNormalForm(to_snfmat(lattice_vectors))),
	point_index(UPV.D),link_index(UPV.D),plaq_index(UPV.D),vol_index(UPV.D)
{
	// populate the cells
	for (const auto& p : other.points){	add_point(p); }
	for (const auto& p : other.links){	add_link(p); }
	for (const auto& p : other.plaqs){	add_plaq(p); }
	for (const auto& p : other.vols){	add_vol(p); }
}

void UnitCellSpecifier::wrap(ipos_t& X) const {
	auto [quot, rem] = moddiv(UPV.L * X, UPV.D);
	X = UPV.Linv * rem;
}


void UnitCellSpecifier::wrap(GeometricObject& X) const {
	wrap(X.position);
}


void UnitCellSpecifier::add_point(const PointSpec& p){
	try_insert_position(point_index, p.position, points.size(), "point");	
	points.push_back(p);
	wrap(points.back());
}

void UnitCellSpecifier::add_link(const LinkSpec& p){
	try_insert_position(link_index, p.position, links.size(), "link"); // store index now
	links.push_back(p);
	wrap(links.back());
	// check that boundary is resolvable	
	for (const auto& dp : p.boundary){
		ipos_t x = wrap_copy(p.position + dp.relative_position);
		if (!is_point(x)){ throw_bad_boundary(p,x); }
	}
}

void UnitCellSpecifier::add_plaq(const PlaqSpec& p){
	try_insert_position(plaq_index, p.position, plaqs.size(), "plaq"); // store index now
	plaqs.push_back(p);
	wrap(plaqs.back());
	// check that boundary is resolvable	
	for (const auto& dp : p.boundary){
		ipos_t x = wrap_copy(p.position + dp.relative_position);
		if (!is_link(x)){ throw_bad_boundary(p,x); }
	}
}



void UnitCellSpecifier::add_vol(const VolSpec& p){
	try_insert_position(vol_index, p.position, vols.size(), "vol"); // store index now
	vols.push_back(p);
	wrap(vols.back());
	// check that boundary is resolvable	
	for (const auto& dp : p.boundary){
		ipos_t x = wrap_copy(p.position + dp.relative_position);
		if (!is_vol(x)){ throw_bad_boundary(p,x); }
	}
}


};
