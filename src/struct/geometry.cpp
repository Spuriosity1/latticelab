#include "geometry.hpp"



// Utility functions
const ipos_t sl_vectors[4] = {
		{0,0,0},
		{0,4,4},
		{4,0,4},
		{4,4,0}
};

bool modeq(ipos_t x, int b, int v){
	return (x[0] % b == v) && (x[1] % b == v) && (x[2] % b == v);
}

// Wraps correctly for negative ints
arma::sword modulo(arma::sword x, arma::sword base){
	return (x%base + base) % base;
}

// entrywise modulo
arma::ivec3 modulo(arma::ivec3 x, arma::ivec3 base){
	return arma::ivec3({
			modulo(x[0], base[0]), 
			modulo(x[1], base[1]), 
			modulo(x[2], base[2])});
}

idx_t pyro_unitcell::num_spins() const {
	return L[0]*L[1]*L[2]*4;
}

spin& pyro_unitcell::spin_at(const arma::ivec3& I, idx_t spin_sl) {
	return spin_sites[spin_idx_from_cell(I, spin_sl)];
}

void pyro_unitcell::construct_spins() {

	spin_sites.resize(num_spins());
	for (int i0=0; i0<L[0]; i0++){	
		for (int i1=0; i1<L[0]; i1++){
			for (int i2=0; i2<L[0]; i2++){
				// The unitcell index
				arma::ivec3 I = {i0, i1, i2};

				// Make the spins
				for (idx_t mu=0; mu<4; mu++){
					spin_at(arma::ivec3({i0,i1,i2}), mu).R = 
						primitive_vectors * I + sl_vectors[mu];
				}
			}
		}
	}
}


idx_t pyro_unitcell::spin_idx_from_cell(const arma::ivec3& I, idx_t spin_sl) const {
	assert(spin_sl < 4);
	return ((I[0]*L[1]+I[1])*L[2] + I[2])*4 + spin_sl;
}


arma::ivec3 pyro_unitcell::cell_IDX_from_coord(const ipos_t& R){
	// Exploit the SNF of the primitve vectors [0 8 8; 8 0 8; 8 8 0]
	// Note that in armadillo % is entrywise product
	const arma::imat33 u = {{ 0, 1, 0,},{ 1, 0, 0}, {1, 1, -1}};
	const arma::imat33 v = {{ 1, 0, -1}, {0, 1, -1}, { 0, 0, 1 }};
	const arma::ivec3  d = { 8, 8, 16 };

	return   modulo( v* ( (u * R) / d), L);
	// Explanation: 
	// R = primitive_vectors * I + pyro_sl
	//  where pyro_sl all lie strictly within the unit cell
	//   = u^-1 * diag(d) * v^-1 * I
	//  u * R + u*pyro_sl = diag(d) * v^-1 * I
	//  Now do entrywise integer (i.e. floor) division (this works even if we are off a perfect site)
	//  (u * R) / d + (u*pyro_sl)/d = v^-1 I
	//  and finally change basis and wrap the unit cell
	//  v * ((u*R) / d) = I
	
}


idx_t pyro_unitcell::spin_idx_from_coord(const ipos_t& R){
	int spin_sl = 0;
	ipos_t tmp;
	for (; spin_sl < 4; spin_sl++){
		tmp = R - sl_vectors[spin_sl];
		if (modeq(tmp, 8, 0)) {
			break;
		}
	}
	// tmp now stores the primitive cell origin
	//
	auto I = cell_IDX_from_coord(R);
	return (((I[0])*L[1] + I[1])*L[2] + I[2])*4 + spin_sl;
}





