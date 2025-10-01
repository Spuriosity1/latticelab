#include <gtest/gtest.h>
#include "UnitCellSpecifier.hpp"
#include "cell_geometry.hpp"

using namespace CellGeometry;


TEST(AbstractTest, DistanceCubicSmoke) {
    UnitCellSpecifier cell(imat33_t::from_cols({1, 0, 0}, {0, 1, 0}, {0, 0, 1}));
    PeriodicAbstractLattice lat(cell, imat33_t::from_cols({8, 0, 0}, {0, 8, 0}, {0, 0, 8}));

    ipos_t x0{0,0,0};
    ASSERT_EQ(lat.distance({1,2,3}, x0), ipos_t(1,2,3));
    ASSERT_EQ(lat.distance(x0, {1,2,3}), ipos_t(-1,-2,-3));

    ASSERT_EQ(lat.distance({9,-6,3}, x0), ipos_t(1,2,3));
}


TEST(AbstractTest, DistanceCubic2D) {
    UnitCellSpecifier cell(imat33_t::from_cols({1, 0, 0}, {0, 1, 0}, {0, 0, 1}));
    PeriodicAbstractLattice lat(cell, imat33_t::from_cols({4, -4, 0}, {4, 4, 0}, {0, 0, 8}));

    ASSERT_EQ(lat.distance({2,0,0}, {5,3,0}), ipos_t(1,1,0));
    ASSERT_EQ(lat.distance({5,3,0}, {2,0,0}), ipos_t(-1,-1,0));

    ASSERT_EQ(lat.distance({2,0,6}, {5,3,0}), ipos_t(1,1,-2));
    ASSERT_EQ(lat.distance({5,3,6}, {2,0,0}), ipos_t(-1,-1,-2));
}


