#include "struct/chain.hpp"
#include <cstdlib>
#include <gtest/gtest.h>
#include <struct/UnitCellSpecifier.hpp>
#include <struct/cell_geometry.hpp>
#include <unordered_set>
#include <vector>

#ifndef DEBUG
#warning "Compiling without DEBUG defined; bounds checking will not throw"
#endif

using namespace CellGeometry;
/*
typedef pointMap<ipos_t, sl_t> pmap_t;


TEST(PointMap, slMap) {
	pmap_t pmap(ipos_t({3,4,5}));
	ipos_t x({0,0,0});
	ipos_t y({2,3,4});
	pmap.insert(x, 4);
	pmap.insert(y, 5);
	EXPECT_EQ(pmap[x], 4);
	EXPECT_EQ(pmap[y], 5);
};


TEST(PointMap, slMapDoubling) {
	pmap_t pmap(ipos_t({3,4,5}));
	ipos_t x({0,0,0});
	pmap.insert(x, 4);
	EXPECT_EQ(pmap.insert(x, 5), false);
};

TEST(PointMap, slMapBadAccess) {
	pmap_t pmap(ipos_t({3,4,5}));
	EXPECT_THROW(pmap[ipos_t({-1,0,0})],std::out_of_range);
	EXPECT_THROW(pmap[ipos_t({0,-5,0})],std::out_of_range);
	EXPECT_THROW(pmap[ipos_t({0,0,-5})],std::out_of_range);
	EXPECT_THROW(pmap[ipos_t({3,0,0})],std::out_of_range);
	EXPECT_THROW(pmap[ipos_t({0,4,0})],std::out_of_range);
	EXPECT_THROW(pmap[ipos_t({0,0,5})],std::out_of_range);
};
*/


typedef PeriodicPointLattice<Cell<0>> PeriodicPointLattice_std;
typedef PeriodicLinkLattice<Cell<0>,Cell<1>> PeriodicLinkLattice_std;
typedef PeriodicPlaqLattice<Cell<0>,Cell<1>,Cell<2>> PeriodicPlaqLattice_std;
typedef PeriodicVolLattice<Cell<0>,Cell<1>,Cell<2>,Cell<3>> PeriodicVolLattice_std;

class CubicPointsTest : public testing::Test {
	protected:
		CubicPointsTest()
			: cell(imat33_t::from_cols({2, 0, 0}, {0, 2, 0}, {0, 0, 2})) {
				PointSpec pointspec;
				pointspec.position = {0, 0, 0};
				cell.add_point(pointspec);
			}
		UnitCellSpecifier cell;
};

class CubicLinksTest : public CubicPointsTest {
	protected:
		CubicLinksTest() {
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
		}

};

class CubicPlaqsTest : public CubicLinksTest {
	protected:
		CubicPlaqsTest(){
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
		}
};

class CubicVolTest : public CubicPlaqsTest {
	protected:
		CubicVolTest(){
			VolSpec volspec;
			volspec.position = {1,1,1};
			volspec.boundary = {
				{1, {1,0,0}}, {-1, {-1,0,0}},
				{1, {0,1,0}}, {-1, {0,-1,0}},
				{1, {0,0,1}}, {-1, {0,0,-1}}
			};	
			cell.add_vol(volspec);
		}
};

TEST_F(CubicPointsTest, AddPointWorks) {
	EXPECT_EQ(cell.num_point_sl(), 1);
	EXPECT_TRUE(cell.is_point({0, 0, 0}));
	EXPECT_TRUE(cell.is_point({2, 0, 0}));
	EXPECT_TRUE(cell.is_point({-2, 0, 0}));
	EXPECT_TRUE(cell.is_point({-2, -2, 16}));
	EXPECT_FALSE(cell.is_point({1, 0, 0}));
	EXPECT_FALSE(cell.is_point({0, 1, 0}));
	EXPECT_FALSE(cell.is_point({0, 0, 1}));
};

TEST_F(CubicLinksTest, AddLinkWorks) {
	EXPECT_EQ(cell.num_link_sl(), 3);
	EXPECT_TRUE(cell.is_link({1,0,0}));
	EXPECT_TRUE(cell.is_link({0,1,0}));
	EXPECT_TRUE(cell.is_link({0,0,1}));

	EXPECT_FALSE(cell.is_link({0,1,1}));
	EXPECT_FALSE(cell.is_link({1,1,0}));
	EXPECT_FALSE(cell.is_link({1,0,1}));
};

TEST_F(CubicPlaqsTest, AddPlaqWorks) {
	EXPECT_EQ(cell.num_plaq_sl(), 3);

	EXPECT_TRUE(cell.is_plaq({0,1,1}));
	EXPECT_TRUE(cell.is_plaq({1,1,0}));
	EXPECT_TRUE(cell.is_plaq({1,0,1}));


	EXPECT_FALSE(cell.is_link({0,1,1}));
	EXPECT_FALSE(cell.is_link({1,1,0}));
	EXPECT_FALSE(cell.is_link({1,0,1}));
};

TEST_F(CubicVolTest, AddVolWorks) {
	EXPECT_EQ(cell.num_vol_sl(), 1);
	EXPECT_TRUE(cell.is_vol({1,1,1}));
};

// Pyrochlore fixture
//

class PyroPointsTest : public testing::Test {
	protected:
		PyroPointsTest()
			: cell(imat33_t::from_cols({0,4,4},	{4,0,4},{4,4,0})) {			
				PointSpec pointspec;
				for (const auto& x : point_positions){
					pointspec.position = x;
					cell.add_point(pointspec);
				}
			}
		UnitCellSpecifier cell;

		const std::vector<ipos_t> point_positions = {
			{0,0,0},{2,2,2}
		};
};

class PyroLinksTest : public PyroPointsTest {
	protected:
		PyroLinksTest() {
			LinkSpec linkspec;
			for (const auto& x : link_positions){
				linkspec.position = x;
				linkspec.boundary = {{1, x}, {-1, -x}};
				cell.add_link(linkspec);
			}
		}


		const std::vector<ipos_t> link_positions = {
			{1,1,1}, {1,-1,-1}, {-1,1,-1}, {-1,-1,1}
		};	

};

class PyroPlaqsTest : public PyroLinksTest {
	protected:
		PyroPlaqsTest(){
			PlaqSpec plaqspec;
			for (int mu=0; mu<4; mu++){
				plaqspec.position = plaq_positions[mu];
				plaqspec.boundary = plaq_boundaries[mu];
				cell.add_plaq(plaqspec);
			}
		}

		const ipos_t plaq_positions[4]  = {
			{-3,-3,-3}, {-3,-1,-1}, {-1,-3,-1}, {-1,-1,-3}
			//{-1,-1,-1}, {-1, 1, 1}, { 1,-1, 1}, { 1, 1,-1}
		};

		const std::vector<VectorSignPair> plaq_boundaries[4] = {
			{
				{ 1,{0, -2, 2}},    {-1, {2, -2, 0}},
				{ 1,{2, 0, -2}},    {-1, {0, 2, -2}}, 
				{ 1,{-2, 2, 0}},    {-1, {-2, 0, 2}}
			},
			{
				{ 1,{ 0, 2,-2}},    {-1,{ 2, 2, 0}},
				{ 1,{ 2, 0, 2}},    {-1,{ 0,-2, 2}},
				{ 1,{-2,-2, 0}},    {-1,{-2, 0,-2}}
			},
			{
				{ 1,{ 0,-2,-2}},	{-1,{-2,-2, 0}},
				{ 1,{-2, 0, 2}},	{-1,{ 0, 2, 2}},
				{ 1,{ 2, 2, 0}},	{-1,{ 2, 0,-2}}
			},
			{
				{ 1,{ 0, 2, 2}},	{-1,{-2, 2, 0}},
				{ 1,{-2, 0,-2}},	{-1,{ 0,-2,-2}},
				{ 1,{ 2,-2, 0}},	{-1,{ 2, 0, 2}}
			}
		};

};

class PyroVolTest : public PyroPlaqsTest {
	protected:
		PyroVolTest(){
			VolSpec volspec;
			volspec.position = vol_positions[0];
			volspec.boundary = {
				{1, { 1, 1, 1}}, {1, { 1,-1,-1}},
				{1, {-1, 1,-1}}, {1, {-1,-1, 1}},
			};	
			cell.add_vol(volspec);
			volspec.position = vol_positions[1];
			volspec.boundary = {
				{-1, {-1,-1,-1}}, {-1, {-1, 1, 1}},
				{-1, { 1,-1, 1}}, {-1, { 1, 1,-1}},
			};
			cell.add_vol(volspec);
		}


		const std::vector<ipos_t> vol_positions = {
			{4,4,4},{6,6,6}
		};
};

TEST_F(PyroPointsTest, AddPointWorks){
	for (const auto& x : point_positions){
		EXPECT_TRUE(cell.is_point(x));
	}
}


TEST_F(PyroLinksTest, AddLinkWorks) {
	for (const auto& x : link_positions){
		EXPECT_TRUE(cell.is_link(x));
	}
}

TEST_F(PyroPlaqsTest, AddPlaqWorks){
	for (const auto& x : plaq_positions){
		EXPECT_TRUE(cell.is_plaq(x));
	}
}

TEST_F(PyroVolTest, AddVolWorks){
	EXPECT_EQ(cell.num_vol_sl(), 2);
	EXPECT_TRUE(cell.is_vol(vol_positions[0]));
	EXPECT_TRUE(cell.is_vol(vol_positions[1]));
}


// specifier is probably OK, move on
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/////// POINT LATTICE   ////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
TEST_F(PyroPointsTest, ConstructAbstractPrimitive){
	PeriodicAbstractLattice lat(cell, 
			imat33_t::from_cols({1,0,0},{0,1,0},{0,0,1})
			);
	EXPECT_EQ(lat.num_primitive, 1);
	std::vector<ipos_t> R = {
		{0,0,0},{1,1,1},{2,2,2},{-1,0,0},{-1,-1,2},
		{4,5,6},{-4,1,3}
	};
	for (auto& r : R){
		EXPECT_EQ(lat.get_supercell_IDX(r), idx3_t({0,0,0}));
	}
}


TEST_F(PyroPointsTest, ConstructAbstractPyro){
	PeriodicAbstractLattice lat(cell, 
			imat33_t::from_cols({-1,1,1},{1,-1,1},{1,1,-1})
			);
	EXPECT_EQ(lat.num_primitive, 4);
}

TEST_F(PyroPointsTest, ConstructPointPrim){
	PeriodicPointLattice_std lat(cell, 
			imat33_t::from_cols({1,0,0},{0,1,0},{0,0,1})
			);
	std::unordered_set<ipos_t> positions;
	for (size_t j =0; j < lat.points.size(); j++ ) {
		auto& R = lat.points[j].position;
		positions.insert(R);
		EXPECT_EQ(&lat.get_point_at(R), &lat.points[j]);
		EXPECT_EQ(&lat.get_point_at(R - lat.cell_vectors * idx3_t({-1,2,-3}) ),
				&lat.points[j]);
	}
	ASSERT_EQ(positions.size(), lat.points.size());
}


TEST_F(PyroPointsTest, ConstructPointCube){
	PeriodicPointLattice_std lat(cell, 
			imat33_t::from_cols({-1,1,1},{1,-1,1},{1,1,-1})
			);
	std::unordered_set<ipos_t> positions;
	for (size_t j =0; j < lat.points.size(); j++ ) {
		auto& R = lat.points[j].position;
		positions.insert(R);
		EXPECT_EQ(&lat.get_point_at(R), &lat.points[j]);
		EXPECT_EQ(&lat.get_point_at(R - lat.cell_vectors * idx3_t({-1,2,-3}) ),
				&lat.points[j]);
	}
	ASSERT_EQ(positions.size(), lat.points.size());
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/////// Link LATTICE   ////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///
///


TEST_F(PyroLinksTest, ConstructLinkPrim){
	PeriodicLinkLattice_std lat(cell, 
			imat33_t::from_cols({1,0,0},{0,1,0},{0,0,1})
			);
	std::unordered_set<ipos_t> positions;
	for (size_t j =0; j < lat.links.size(); j++ ) {
		auto& R = lat.links[j].position;
		positions.insert(R);
		EXPECT_EQ(&lat.get_link_at(R), &lat.links[j]);
		EXPECT_EQ(&lat.get_link_at(R - lat.cell_vectors * idx3_t({-1,2,-3}) ),
				&lat.links[j]);
	}
	ASSERT_EQ(positions.size(), lat.links.size());
}


TEST_F(PyroLinksTest, ConstructLinkCube){
	PeriodicLinkLattice_std lat(cell, 
			imat33_t::from_cols({-1,1,1},{1,-1,1},{1,1,-1})
			);
	std::unordered_set<ipos_t> positions;
	for (size_t j =0; j < lat.links.size(); j++ ) {
		auto& R = lat.links[j].position;
		positions.insert(R);
		EXPECT_EQ(&lat.get_link_at(R), &lat.links[j]);
		EXPECT_EQ(&lat.get_link_at(R - lat.cell_vectors * idx3_t({-1,2,-3}) ),
				&lat.links[j]);
	}
	ASSERT_EQ(positions.size(), lat.links.size());
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/////// Plaq LATTICE   ////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/// These check that the lattice indexor finctions behave as expected


TEST_F(PyroPlaqsTest, ConstructPlaqPrim){
	PeriodicPlaqLattice_std lat(cell, 
			imat33_t::from_cols({1,0,0},{0,1,0},{0,0,1})
			);
	std::unordered_set<ipos_t> positions;
	for (size_t j =0; j < lat.plaqs.size(); j++ ) {
		auto& R = lat.plaqs[j].position;
		positions.insert(R);
		EXPECT_EQ(&lat.get_plaq_at(R), &lat.plaqs[j]);
		EXPECT_EQ(&lat.get_plaq_at(R - lat.cell_vectors * idx3_t({-1,2,-3}) ),
				&lat.plaqs[j]);
	}
	ASSERT_EQ(positions.size(), lat.plaqs.size());
}


TEST_F(PyroPlaqsTest, ConstructPlaqCube){
	PeriodicPlaqLattice_std lat(cell, 
			imat33_t::from_cols({-1,1,1},{1,-1,1},{1,1,-1})
			);
	std::unordered_set<ipos_t> positions;
	for (size_t j =0; j < lat.plaqs.size(); j++ ) {
		auto& R = lat.plaqs[j].position;
		positions.insert(R);
		EXPECT_EQ(&lat.get_plaq_at(R), &lat.plaqs[j]);
		EXPECT_EQ(&lat.get_plaq_at(R - lat.cell_vectors * idx3_t({-1,2,-3}) ),
				&lat.plaqs[j]);
	}
	ASSERT_EQ(positions.size(), lat.plaqs.size());
}
