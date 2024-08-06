#pragma once 

#include <armadillo>
#include <vector>
#include <type_traits>


typedef unsigned int idx_t;
typedef arma::ivec3 ipos_t;
struct tetra;

struct spin {
	ipos_t R;
	tetra* neighbours[2];
};

struct tetra {
	ipos_t R;
	spin* spins[4];
};

struct pyro_unitcell {
	const arma::ivec3 L;
	std::vector<spin> spin_sites;

	// constructors
	pyro_unitcell(const arma::ivec3& _L) : L(_L) {
		construct_spins();
	}
	pyro_unitcell(int Lx, int Ly, int Lz) :
		L( arma::ivec3({Lx,Ly,Lz}) ) {
			construct_spins();
		}

	spin& spin_at(const arma::ivec3& I, idx_t spin_sl);

	idx_t num_spins() const;
	idx_t num_tetras() const;

	// Spin index finders
	idx_t spin_idx_from_cell(const arma::ivec3& I, idx_t spin_sl) const;
	idx_t spin_idx_from_coord(const ipos_t& r);

	idx_t tetra_idx_from_cell(const arma::ivec3& I, idx_t tetra_sl) const;
	idx_t tetra_idx_from_coord(const ipos_t& r);
	
	/**
	 * Computes and returns the cell index (a three-component integer vector 
	 * of nonnegative indices) of a particular point in space
	 *
	 */
	arma::ivec3 cell_IDX_from_coord(const ipos_t& R);


	const arma::imat33 primitive_vectors = {{0,8,8},{8,0,8},{8,8,0}};
	arma::imat33 supercell_vectors() {
		arma::imat33 retval(primitive_vectors);
		retval.col(0) *= L[0];
		retval.col(1) *= L[1];
		retval.col(2) *= L[2];
		return retval;
	}



	private:

	void construct_spins();

	// The SNF decomposition of the matrix [0 8 8; 8 0 8; 8 8 0]
	
};
