#include <iostream>
#include <UnitCellSpecifier.hpp>
#include <cstdlib>
#include <sstream>
#include "cell_geometry.hpp"
#include "chain.hpp"
#include "preset_cellspecs.hpp"
#include "basic_parser.hh"
#include "lattice_IO.hpp"

using namespace CellGeometry;

void parse_supercell_spec(imat33_t& supercell_spec,
        std::string Z1_s, std::string Z2_s, std::string Z3_s){

    std::stringstream Z1_ss(Z1_s);
    std::stringstream Z2_ss(Z2_s);
    std::stringstream Z3_ss(Z3_s);
    for (int row=0; row<3; row++){
        Z1_ss >> supercell_spec(row,0);
        Z2_ss >> supercell_spec(row,1);
        Z3_ss >> supercell_spec(row,2);
    }
}


int main (int argc, const char *argv[]) {
    std::string Z1_s,Z2_s,Z3_s;
    std::filesystem::path outpath;

    basic_parser::Parser args(1,1);
    args.declare("Z1",&Z1_s);
    args.declare("Z2",&Z2_s);
    args.declare("Z3",&Z3_s);

    args.from_file(argv[1]);
    outpath = argv[2];
    args.from_argv(argc, argv, 3);

    args.assert_initialised();

    // parse L1 L2 L3
    imat33_t supercell_spec;
    parse_supercell_spec(supercell_spec, Z1_s, Z2_s, Z3_s);

    std::cout<<"Constructing supercell of dimensions \n"<<supercell_spec<<std::endl;

    // const PrimitiveSpecifers::Diamond spec;
    const auto spec = PrimitiveSpecifiers::CubicSpec();

    // Construct the periodic lattice extension using these extensions
    
    PeriodicLinkLattice<Cell<0>,Cell<1>> point_lat(spec, supercell_spec);

    save(point_lat, outpath/"link_lat.json");


    
}

