#include <gtest/gtest.h>
#include <struct/geometry.hpp>



TEST(IndexTest, SpinIndexConsistency) {
	pyro_unitcell cell(2,2,2);
	for (idx_t i=0; i<cell.num_spins(); i++){
	EXPECT_EQ(i, cell.spin_idx_from_coord(cell.spin_sites[i].R));
	}
}

TEST(IndexTest, SpinIndexShiftConsistency) {
	pyro_unitcell cell(2,2,2);
	for (idx_t i=0; i<cell.num_spins(); i++){
		EXPECT_EQ(i,cell.spin_idx_from_coord(cell.spin_sites[i].R + cell.supercell_vectors() * arma::ivec3({-1, 0, 0}) ) );
		EXPECT_EQ(i,cell.spin_idx_from_coord(cell.spin_sites[i].R + cell.supercell_vectors() * arma::ivec3({-1,-1, 0}) ) );
		EXPECT_EQ(i,cell.spin_idx_from_coord(cell.spin_sites[i].R + cell.supercell_vectors() * arma::ivec3({-1,-1,-1}) ) );
	}
}

