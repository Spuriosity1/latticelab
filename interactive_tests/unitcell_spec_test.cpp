#include <iostream>
#include <struct/UnitCellSpecifier.hpp>
#include <cstdlib>
#include <sstream>
#include "struct/cell_geometry.hpp"
#include "struct/chain.hpp"
#include "struct/preset_cellspecs.hpp"
#include "admin/basic_parser.hh"
#include "struct/lattice_IO.hpp"

using namespace CellGeometry;

void parse_supercell_spec(imat33_t& supercell_spec,
        std::string L1_s, std::string L2_s, std::string L3_s){

    std::stringstream L1_ss(L1_s);
    std::stringstream L2_ss(L2_s);
    std::stringstream L3_ss(L3_s);
    for (int row=0; row<3; row++){
        L1_ss >> supercell_spec(row,0);
        L2_ss >> supercell_spec(row,1);
        L3_ss >> supercell_spec(row,2);
    }
}


int main (int argc, const char *argv[]) {
    std::string L1_s,L2_s,L3_s;
    std::filesystem::path outpath;

    basic_parser::Parser args(1,1);
    args.declare("Z1",&L1_s);
    args.declare("Z2",&L2_s);
    args.declare("Z3",&L3_s);

    args.from_file(argv[1]);
    outpath = argv[2];
    args.from_argv(argc, argv, 3);

    args.assert_initialised();

    // parse L1 L2 L3
    imat33_t supercell_spec;
    parse_supercell_spec(supercell_spec, L1_s, L2_s, L3_s);

    std::cout<<"Constructing supercell of dimensions \n"<<supercell_spec<<std::endl;

    // const PrimitiveSpecifers::Diamond spec;
    const PrimitiveSpecifers::Cubic spec;

    // Construct the periodic lattice extension using these extensions
    
    PeriodicLinkLattice<Cell<0>,Cell<1>> point_lat(spec, supercell_spec);

    save(point_lat, outpath/"link_lat.json");


    
}

