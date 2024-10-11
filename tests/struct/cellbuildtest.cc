#include "struct/chain.hpp"
#include <cstdlib>
#include <gtest/gtest.h>
#include <stdexcept>
#include <struct/cell_geometry.hpp>

#ifndef DEBUG
#warning "Compiling without DEBUG defined; bounds checking will not throw"
#endif


using namespace CellGeometry;

struct CustomPoint : public Cell<0> {
	double data1;
	double data2;
};



class CustomPointsTest : public testing::Test {
	protected:
		CustomPointsTest()
			: cell(imat33_t::from_cols({2, 0, 0}, {0, 2, 0}, {0, 0, 2})) {
				PointSpec pointspec;
				pointspec.position = {0, 0, 0};
				cell.add_point(pointspec);
			}
		UnitCellSpecifier cell;
};

TEST_F(CustomPointsTest, ConstructPointLattice){
	PeriodicPointLattice<CustomPoint> pointlat(cell, 
			imat33_t::from_cols({4,-2,0},{-4,-3,4},{0,2,-1})
			);	
	auto p = pointlat.get_point_at({0,0,0});
}



