#pragma once

#include "UnitCellSpecifier.hpp"

namespace CellGeometry {
namespace PrimitiveSpecifers {

struct Pyrochlore : public UnitCellSpecifier {
	Pyrochlore() : UnitCellSpecifier(
		imat33_t::from_cols({0, 4, 4}, {4, 0, 4}, {4, 4, 0}) )
	{
		setup_points();
		setup_links();
		setup_plaqs();
		setup_vols();
	}
private:
	void setup_points();
	void setup_links();
	void setup_plaqs();
	void setup_vols();
};


/*
struct Cubic : public UnitCellSpecifier {
	Cubic() : UnitCellSpecifier(
		imat33_t::from_cols({2,0,0},{0,2,0},{0,0,2}) )
	{
		setup_points();
		setup_links();
		setup_plaqs();
		setup_vols();
	}
private:
	void setup_points();
	void setup_links();
	void setup_plaqs();
	void setup_vols();
};
*/

};
};


