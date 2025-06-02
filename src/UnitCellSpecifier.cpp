#include "UnitCellSpecifier.hpp"
#include "chain.hpp"
#include "modulus.hpp"
#include <sstream>
#include <stdexcept>


#ifdef DEBUG
#define ASSERT_VALID_POS(R) assert_valid_position(R);
#else
#define ASSERT_VALID_POS(R);
#endif


namespace CellGeometry {

// #define PPRINT(x) "["<<x[0]<<" "<<x[1]<<" "<<x[2]<<"]"



// Sublattice index access from a physical position (real units)
// Linear search suboptimal here, but this fn should not be performance
// critical.
sl_t UnitCellSpecifier::sl_of_point(const ipos_t& R_) const {
	auto R = wrap_copy(R_);
	for (size_t i=0; i<points.size(); i++) { if (points[i].position == R) return i; }
	std::stringstream s("No point found at "); s << R;
	throw std::out_of_range(s.str()); 
}
sl_t UnitCellSpecifier::sl_of_link(const ipos_t& R_) const {
	auto R = wrap_copy(R_);
	for (size_t i=0; i<links.size(); i++) { if (links[i].position == R) return i; }
	std::stringstream s("No link found at "); s << R;
	throw std::out_of_range(s.str()); 
}
sl_t UnitCellSpecifier::sl_of_plaq(const ipos_t& R_) const {
	auto R = wrap_copy(R_);
	for (size_t i=0; i<plaqs.size(); i++) { if (plaqs[i].position == R) return i; }
	std::stringstream s("No plaq found at "); s << R;
	throw std::out_of_range(s.str()); 
}
sl_t UnitCellSpecifier::sl_of_vol(const ipos_t& R_) const {
	auto R = wrap_copy(R_);
	for (size_t i=0; i<vols.size(); i++) { if (vols[i].position == R) return i; }
	std::stringstream s("No vol found at "); s << R;
	throw std::out_of_range(s.str()); 
}

// testers -- linear search, since these lists should be very small
// Consider upgrading to binary search or better
bool UnitCellSpecifier::is_point(const ipos_t& R_){
	auto R = wrap_copy(R_);
	ASSERT_VALID_POS(R);
	for (size_t i=0; i<points.size(); i++) { if (points[i].position == R) return true; }
	return false;
}
bool UnitCellSpecifier::is_link(const ipos_t& R_){
	auto R = wrap_copy(R_);
	ASSERT_VALID_POS(R);
	for (size_t i=0; i<links.size(); i++) { if (links[i].position == R) return true; }
	return false;
}
bool UnitCellSpecifier::is_plaq(const ipos_t& R_){
	auto R = wrap_copy(R_);
	ASSERT_VALID_POS(R);
	for (size_t i=0; i<plaqs.size(); i++) { if (plaqs[i].position == R) return true; }
	return false;
}
bool UnitCellSpecifier::is_vol(const ipos_t& R_){
	auto R = wrap_copy(R_);
	ASSERT_VALID_POS(R);
	for (size_t i=0; i<vols.size(); i++) { if (vols[i].position == R) return true; }
	return false;
}


template<int _order>
void throw_bad_boundary(CellSpecifier<_order> p, ipos_t missing_point) {
	std::stringstream s;
	s << "Problem with boundary of "<<_order<<"-cell specifier ";
	s << "at "<<p.position<<": ";
	s << "There is no "<<_order-1<<"-cell at " << missing_point;
	throw std::out_of_range(s.str());
}


// Computes the closed-form, simplified matrix inverse
constexpr imat33_t unnormed_inverse(const imat33_t& A){
	std::array<int64_t, 3> b0 = {-A(1, 2) * A(2, 1) + A(1, 1) * A(2, 2),
		A(0, 2) * A(2, 1) - A(0, 1) * A(2, 2),
		-A(0, 2) * A(1, 1) + A(0, 1) * A(1, 2)};
	std::array<int64_t, 3> b1 = {A(1, 2) * A(2, 0) - A(1, 0) * A(2, 2),
		-A(0, 2) * A(2, 0) + A(0, 0) * A(2, 2),
		A(0, 2) * A(1, 0) - A(0, 0) * A(1, 2)};
	std::array<int64_t, 3> b2 = {-A(1, 1) * A(2, 0) + A(1, 0) * A(2, 1),
		A(0, 1) * A(2, 0) - A(0, 0) * A(2, 1),
		-A(0, 1) * A(1, 0) + A(0, 0) * A(1, 1)};
	return imat33_t::from_rows(b0,b1,b2);
}


constexpr imat33_t make_positive(const imat33_t& mat){
	if (det(mat) < 0) { return -1*mat; }
	return mat;
}

// constructors
//
UnitCellSpecifier::UnitCellSpecifier(const imat33_t& lattice_vectors_) :
	latvecs(make_positive(lattice_vectors_)),
	latvecs_unnormed_inverse( unnormed_inverse(latvecs) ),
	abs_det_latvecs(det(latvecs))
	//UPV(ComputeSmithNormalForm(to_snfmat(lattice_vectors)))
//	point_index(UPV.D),link_index(UPV.D),plaq_index(UPV.D),vol_index(UPV.D)
{}

UnitCellSpecifier::UnitCellSpecifier(
		const UnitCellSpecifier& other,
		const imat33_t& cellspec) : 
	latvecs(make_positive(other.latvecs * imat33_t::from_other(cellspec))),
	latvecs_unnormed_inverse( unnormed_inverse(latvecs) ),
	abs_det_latvecs(det(latvecs))
	// UPV(ComputeSmithNormalForm(to_snfmat(lattice_vectors)))
	// point_index(UPV.D),link_index(UPV.D),plaq_index(UPV.D),vol_index(UPV.D)
{
	if (abs(det(cellspec)) != 1) {
		throw std::invalid_argument(
				"cellspec of re-parameterised unit cell must have determinant 1");
	}
	// populate the cells
	for (const auto& p : other.points){	add_point(p); }
	for (const auto& p : other.links){	add_link(p); }
	for (const auto& p : other.plaqs){	add_plaq(p); }
	for (const auto& p : other.vols){	add_vol(p); }
}

void UnitCellSpecifier::wrap(ipos_t& X) const {
	ipos_t x = latvecs_unnormed_inverse * X; // this / det(A) is the true x
	for (int i=0; i<3; i++){
		x[i] = mod(x[i], abs_det_latvecs);
	}
	X = latvecs*x;
	for (int i=0; i<3; i++){
		assert(X[i] % abs_det_latvecs == 0);
		X[i] /= abs_det_latvecs;
	}
}


void UnitCellSpecifier::wrap(GeometricObject& X) const {
	wrap(X.position);
}


void UnitCellSpecifier::add_point(const PointSpec& p){
	ASSERT_VALID_POS(p.position);
	points.push_back(p);
	wrap(points.back());
}

void UnitCellSpecifier::add_link(const LinkSpec& p){
	ASSERT_VALID_POS(p.position);
	links.push_back(p);
	wrap(links.back());
	// check that boundary is resolvable	
	for (const auto& dp : p.boundary){
		ipos_t x = wrap_copy(p.position + dp.relative_position);
		if (!is_point(x)){ throw_bad_boundary(p,x); }
	}
}

void UnitCellSpecifier::add_plaq(const PlaqSpec& p){
	ASSERT_VALID_POS(p.position);
	plaqs.push_back(p);
	wrap(plaqs.back());
	// check that boundary is resolvable	
	for (const auto& dp : p.boundary){
		ipos_t x = wrap_copy(p.position + dp.relative_position);
		if (!is_link(x)){ throw_bad_boundary(p,x); }
	}
}

void UnitCellSpecifier::add_vol(const VolSpec& p){
	ASSERT_VALID_POS(p.position);
	vols.push_back(p);
	wrap(vols.back());
	// check that boundary is resolvable	
	for (const auto& dp : p.boundary){
		ipos_t x = wrap_copy(p.position + dp.relative_position);
		if (!is_plaq(x)){ throw_bad_boundary(p,x); }
	}
}


};
