#include <gtest/gtest.h>
#include "UnitCellSpecifier.hpp"
#include "cell_geometry.hpp"
#include <random>

using namespace CellGeometry;


TEST(AbstractTest, DetWorks) {
    UnitCellSpecifier cell(imat33_t::from_cols({4,-4,4},{8,-5,5}, {5,-3,5}));
    PeriodicAbstractLattice lat(cell, imat33_t::from_cols({3,-3,5},{3,-3,2},{8,4,-4}));

    auto primitive = lat.primitive_spec;

    ASSERT_EQ(primitive.latvecs_unnormed_inverse * primitive.latvecs, primitive.abs_det_latvecs
            *imat33_t::eye());
}


TEST(AbstractTest, IndexVecsExpected) {
    
    std::mt19937 rng(12345); // fixed seed for reproducibility
    std::uniform_int_distribution<int> dist(-5, 5);

    imat33_t a;
    imat33_t z;

    for (int test=0; test<1000; test++){

        for (int j=0; j<9; j++){
            a[j] = dist(rng);
            z[j] = dist(rng);
        }
        if (det(a) == 0 || det(z) == 0) continue;
        if (det(a) < 0) a *= -1;

        UnitCellSpecifier cell(a);
        PeriodicAbstractLattice lat(cell, z);

        auto primitive = lat.primitive_spec;
        imat33_t D;

        for (int i=0; i<3; i++) 
            D(i,i) = lat.size(i);

        EXPECT_EQ(lat.index_cell_vectors, lat.primitive_spec.latvecs * D);
    }
}

TEST(AbstractTest, DistanceCubicSmoke) {
    UnitCellSpecifier cell(imat33_t::from_cols({1, 0, 0}, {0, 1, 0}, {0, 0, 1}));
    PeriodicAbstractLattice lat(cell, imat33_t::from_cols({8, 0, 0}, {0, 8, 0}, {0, 0, 8}));

    ipos_t x0{0,0,0};
    EXPECT_EQ(lat.distance({1,2,3}, x0), ipos_t(1,2,3));
    EXPECT_EQ(lat.distance(x0, {1,2,3}), ipos_t(-1,-2,-3));

    EXPECT_EQ(lat.distance({1,2,3}, {5,6,7}), ipos_t(4,4,4));
    EXPECT_EQ(lat.distance({-1,-2,-3}, {3,2,1}), ipos_t(4,4,4));

    EXPECT_EQ(lat.distance({9,-6,3}, x0), ipos_t(1,2,3));
}


TEST(AbstractTest, DistanceCubic211) {
    UnitCellSpecifier cell(imat33_t::from_cols({2, 0, 0}, {0, 1, 0}, {0, 0, 1}));
    PeriodicAbstractLattice lat(cell, imat33_t::from_cols({4, 0, 0}, {0, 8, 0}, {0, 0, 8}));

    ipos_t x0{0,0,0};
    EXPECT_EQ(lat.distance({1,2,3}, x0), ipos_t(1,2,3));
    EXPECT_EQ(lat.distance(x0, {1,2,3}), ipos_t(-1,-2,-3));

    EXPECT_EQ(lat.distance({1,2,3}, {5,6,7}), ipos_t(4,4,4));
    EXPECT_EQ(lat.distance({-1,-2,-3}, {3,2,1}), ipos_t(4,4,4));

    EXPECT_EQ(lat.distance({9,-6,3}, x0), ipos_t(1,2,3));
}


TEST(AbstractTest, DistanceCubic2D) {
    UnitCellSpecifier cell(imat33_t::from_cols({2, 0, 0}, {0, 2, 0}, {0, 0, 2}));
    PeriodicAbstractLattice lat(cell, imat33_t::from_cols({2, -2, 0}, {2, 2, 0}, {0, 0, 4}));

    EXPECT_EQ(lat.distance({2,0,0}, {5,3,0}), ipos_t(1,1,0));
    EXPECT_EQ(lat.distance({5,3,0}, {2,0,0}), ipos_t(-1,-1,0));

    EXPECT_EQ(lat.distance({2,0,6}, {5,3,0}), ipos_t(1,1,-2));
    EXPECT_EQ(lat.distance({5,3,6}, {2,0,0}), ipos_t(-1,-1,-2));
}


TEST(AbstractTest, EdgeCases1) {
    UnitCellSpecifier cell(imat33_t::from_cols({1, 0, 0}, {0, 1, 0}, {0, 0, 1}));
    PeriodicAbstractLattice lat(cell, imat33_t::from_cols({1,-1, 0}, {2,2, 0}, {0, 0, 8}));

    EXPECT_EQ(lat.distance({2,1,0}, {1,-1,0}), ipos_t(0,-1,0));
    EXPECT_EQ(lat.distance({2,1,0}, {3,0,0}), ipos_t(0,0,0));
    EXPECT_EQ(lat.distance({0,-2,0}, {-1,1,0}), ipos_t(2,0,0));
}



TEST(AbstractTest, EdgeCases2) {
    UnitCellSpecifier cell(imat33_t::from_cols({0,4,4},{4,0,4},{4,4,0}));
    PeriodicAbstractLattice lat(cell, imat33_t::from_cols({-5,5,5},{5,-5,5},{5,5,-5}));

    EXPECT_EQ(lat.distance({0,3,8}, {0,3,-8}), ipos_t(0,0,16));
}

//size_t d2(const ipos_t& x){
//    size_t res=0;
//    for (int j=0; j<3; j++){
//        res += x[j]*x[j];
//    }
//    return res;
//}

size_t bruteforce_d2(const PeriodicAbstractLattice& lat, const ipos_t& x, const ipos_t& y){
    size_t min_dist = std::numeric_limits<size_t>::max();
    for (int ix=-1; ix<=1; ix++){
    for (int iy=-1; iy<=1; iy++){
    for (int iz=-1; iz<=1; iz++){
        auto dist = d2_raw(x - y - lat.cell_vectors * ipos_t({ix,iy,iz}));
        if (dist < min_dist){
            min_dist = dist;
        }
    }
    }
    }
    return min_dist;
}

TEST(AbstractTest, Correctness1) {
    UnitCellSpecifier cell(imat33_t::from_cols({1,-1,0},{2,2,0}, {0,0,1}));
    PeriodicAbstractLattice lat(cell, imat33_t::from_cols({3,-3,0},{3,3,0},{0,0,2}));

    std::mt19937 rng(12345); // fixed seed for reproducibility
    std::uniform_int_distribution<int> dist(-3, 3);

    lat.print_diagnostics();

    for (int trial = 0; trial < 200; ++trial) {
        ipos_t x(dist(rng), dist(rng), dist(rng));
        ipos_t y(dist(rng), dist(rng), dist(rng));
        
        auto delta = lat.distance(x, y);

        auto d_normal = d2_raw(delta);
        auto d_brute  = bruteforce_d2(lat, x, y);

        EXPECT_EQ(d_normal, d_brute) 
            << "Mismatch for x=" << x << " y=" << y << " ... delta = " << delta;
    }
}



TEST(AbstractTest, Correctness2) {
    UnitCellSpecifier cell(imat33_t::from_cols({0,4,4},{4,0,4}, {4,4,0}));
    PeriodicAbstractLattice lat(cell, imat33_t::from_cols({-4,4,4},{4,-4,4},{4,4,-4}));

    std::mt19937 rng(12345); // fixed seed for reproducibility
    std::uniform_int_distribution<int> dist(-3, 3);

    lat.print_diagnostics();

    for (int trial = 0; trial < 2; ++trial) {
        ipos_t x(dist(rng), dist(rng), dist(rng));
        ipos_t y(dist(rng), dist(rng), dist(rng));

        auto delta = lat.distance(x, y);

        auto d_normal = d2_raw(delta);
        auto d_brute  = bruteforce_d2(lat, x, y);


        EXPECT_EQ(d_normal, d_brute) 
            << "Mismatch for x=" << x << " y=" << y  << " ... delta = " << delta;
    }
}


