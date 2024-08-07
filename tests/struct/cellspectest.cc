#include "struct/chain.hpp"
#include <cstdlib>
#include <gtest/gtest.h>
#include <stdexcept>
#include <struct/UnitCellSpecifier.hpp>

#ifndef DEBUG
#warning "Compiling without DEBUG defined; bounds checking will not throw"
#endif

using namespace CellGeometry;

TEST(PointMap, slMap) {
	pointMap pmap(ipos_t({3,4,5}));
	ipos_t x({0,0,0});
	ipos_t y({2,3,4});
	pmap.insert(x, 4);
	pmap.insert(y, 5);
	EXPECT_EQ(pmap[x], 4);
	EXPECT_EQ(pmap[y], 5);
};


TEST(PointMap, slMapDoubling) {
	pointMap pmap(ipos_t({3,4,5}));
	ipos_t x({0,0,0});
	pmap.insert(x, 4);
	EXPECT_EQ(pmap.insert(x, 5), false);
};

TEST(PointMap, slMapBadAccess) {
	pointMap pmap(ipos_t({3,4,5}));
	EXPECT_THROW(pmap[ipos_t({-1,0,0})],std::out_of_range);
	EXPECT_THROW(pmap[ipos_t({0,-5,0})],std::out_of_range);
	EXPECT_THROW(pmap[ipos_t({0,0,-5})],std::out_of_range);
	EXPECT_THROW(pmap[ipos_t({3,0,0})],std::out_of_range);
	EXPECT_THROW(pmap[ipos_t({0,4,0})],std::out_of_range);
	EXPECT_THROW(pmap[ipos_t({0,0,5})],std::out_of_range);
};

TEST(UnitCellSpecifier, CellInsert) {
	UnitCellSpecifier cell(imat33_t::from_cols(
				{2,0,0},
				{0,2,0},
				{0,0,2}
				));
	PointSpec pointspec;
	pointspec.position = {0,0,0};	
	cell.add_point(pointspec);

	EXPECT_EQ(cell.num_point_sl(), 1);
	EXPECT_TRUE(cell.is_point({0,0,0}));
	EXPECT_FALSE(cell.is_point({1,0,0}));
	EXPECT_FALSE(cell.is_point({0,1,0}));
	EXPECT_FALSE(cell.is_point({0,0,1}));	

	LinkSpec linkspec;
	linkspec.position = {1,0,0};
	linkspec.boundary = { {1, {-1,0,0}}, {-1, {1,0,0}} };
	cell.add_link(linkspec);

	linkspec.position = {0,1,0};
	linkspec.boundary = { {1, {0,-1,0}}, {-1, {0,1,0}} };
	cell.add_link(linkspec);

	linkspec.position = {0,0,1};
	linkspec.boundary = { {1, {0,0,-1}}, {-1, {0,0,1}} };
	cell.add_link(linkspec);

	EXPECT_EQ(cell.num_link_sl(), 3);
	EXPECT_TRUE(cell.is_link({1,0,0}));
	EXPECT_TRUE(cell.is_link({0,1,0}));
	EXPECT_TRUE(cell.is_link({0,0,1}));

	EXPECT_FALSE(cell.is_link({0,1,1}));
	EXPECT_FALSE(cell.is_link({1,1,0}));
	EXPECT_FALSE(cell.is_link({1,0,1}));

	PlaqSpec plaqspec;
	plaqspec.position = {0,1,1};
	plaqspec.boundary = {
		{1, {0,0,-1}}, {1, {0,1,0}}, {-1, {0,0,-1}}, {-1, {0,-1,0}}
	};
	cell.add_plaq(plaqspec);

	plaqspec.position = {1,0,1};
	plaqspec.boundary = {
		{1, {-1,0,0}}, {1, {0,0,1}}, {-1, {-1,0,0}}, {-1, {0,0,-1}}
	};
	cell.add_plaq(plaqspec);

	plaqspec.position = {1,1,0};
	plaqspec.boundary = {
		{1, {0,-1,0}}, {1, {1,0,0}}, {-1, {0,-1,0}}, {-1, {-1,0,0}}
	};
	cell.add_plaq(plaqspec);

	EXPECT_EQ(cell.num_plaq_sl(), 3);

	EXPECT_TRUE(cell.is_plaq({0,1,1}));
	EXPECT_TRUE(cell.is_plaq({1,1,0}));
	EXPECT_TRUE(cell.is_plaq({1,0,1}));

	VolSpec volspec;
	volspec.position = {1,1,1};
	volspec.boundary = {
		{1, {1,0,0}}, {-1, {-1,0,0}},
		{1, {0,1,0}}, {-1, {0,-1,0}},
		{1, {0,0,1}}, {-1, {0,0,-1}}
	};	
	cell.add_vol(volspec);
	EXPECT_EQ(cell.num_vol_sl(), 1);
	EXPECT_TRUE(cell.is_vol({1,1,1}));
}
