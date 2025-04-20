#include "cell_geometry.hpp"
#include "basic_parser.hh"
#include "lattice_IO.hpp"
#include "preset_cellspecs.hpp"
#include <UnitCellSpecifier.hpp>
#include <algorithm>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

using namespace CellGeometry;


void parse_supercell_spec(imat33_t& supercell_spec, std::string Z1_s, std::string Z2_s, std::string Z3_s){

    std::stringstream Z1_ss(Z1_s);
    std::stringstream Z2_ss(Z2_s);
    std::stringstream Z3_ss(Z3_s);
    for (int row=0; row<3; row++){
        Z1_ss >> supercell_spec(row,0);
        Z2_ss >> supercell_spec(row,1);
        Z3_ss >> supercell_spec(row,2);
    }
}



inline std::string hash_parameters(
        const std::string& Z1_s, 
        const std::string& Z2_s, 
        const std::string& Z3_s,
        double cell_disorder[]
        ){
    // hashing the arguments 
    std::ostringstream name;
    name << "Z1="+Z1_s+";Z2="+Z2_s+";Z3="+Z3_s+";";
    for (int r = 0; r<4; r++){
        if (cell_disorder[r]>0){
            name << "d0="<< std::fixed << std::setprecision(4) <<
                    cell_disorder[r] << ";";
        }
    }
    auto s = name.str();
    std::replace(s.begin(), s.end(), ' ', ',');  // replace space by commas
    return s;
}


int main (int argc, const char *argv[]) {
	std::string Z1_s, Z2_s, Z3_s;
	double cell_disorder[4];

    std::filesystem::path outpath;

    basic_parser::Parser args(1,1);
    args.declare("Z1", &Z1_s);
    args.declare("Z2", &Z2_s);
    args.declare("Z3", &Z3_s);

	
    args.declare_optional("cell0_disorder",cell_disorder, 0.);
    args.declare_optional("cell1_disorder",cell_disorder+1, 0.);
    args.declare_optional("cell2_disorder",cell_disorder+2, 0.);
    args.declare_optional("cell3_disorder",cell_disorder+3, 0.);


    if (argc == 1){
        std::cerr<<"Usage: "<<argv[0]<<" <output_path> <options...>\n";
        return argc;
    }
    outpath = argv[1];
    args.from_argv(argc, argv, 2);

    args.assert_initialised();

    // parse L1 L2 L3
    imat33_t supercell_spec;
    parse_supercell_spec(supercell_spec, Z1_s, Z2_s, Z3_s);

    std::cout<<"Constructing supercell of dimensions \n"<<supercell_spec<<std::endl;

    const PrimitiveSpecifers::Diamond spec;
    auto name = hash_parameters(Z1_s, Z2_s, Z3_s, cell_disorder);



    // Construct the periodic lattice extension using these extensions
    PeriodicPointLattice<Cell<0>> lat0(spec, supercell_spec);
    PeriodicLinkLattice<Cell<0>, Cell<1>> lat1(spec, supercell_spec);
    PeriodicPlaqLattice<Cell<0>, Cell<1>, Cell<2>> lat2(spec, supercell_spec);
    PeriodicVolLattice<Cell<0>, Cell<1>, Cell<2>, Cell<3>> lat3(spec, supercell_spec);

    auto lat = lat3;

    // debug
    lat.print_state();

    
    if (cell_disorder[0]>1e-16){
        std::cout<<"Deleting "<<cell_disorder[0]*100<<"% of 0-cells"<<std::endl;
        for (auto p: lat.get_points()){
            if (rand() > cell_disorder[0]){
                lat.erase_point(&p);
            }

        }
    }

    if (cell_disorder[1]>1e-16){
        std::cout<<"Deleting "<<cell_disorder[1]*100<<"% of 1-cells"<<std::endl;
        for (auto l : lat.get_links()){
            if (rand() > cell_disorder[1]){
                lat1.erase_link(&l);
            }
        }
    }   


    save(lat, outpath/(name+".lat.json"));

    return 0;
}


