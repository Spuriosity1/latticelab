#include "cell_geometry.hpp"
#include "basic_parser.hh"
#include "chain.hpp"
#include "lattice_IO.hpp"
#include "preset_cellspecs.hpp"
#include <UnitCellSpecifier.hpp>
#include <algorithm>
#include <random>
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
            name << "dis"<<r<<"="<< std::fixed << std::setprecision(4) <<
                    cell_disorder[r] << ";";
        }
    }
    auto s = name.str();
    std::replace(s.begin(), s.end(), ' ', ',');  // replace space by commas
    return s;
}


template <int Dim, CellLike<0> Point, CellLike<1> Link, CellLike<2> Plaq, CellLike<3> Vol>
void delete_sites(PeriodicVolLattice<Point, Link, Plaq, Vol>& lat, double disorder_rate) {
}

typedef Cell<0> Point;
typedef Cell<1> Link;
typedef Cell<2> Plaq;
typedef Cell<3> Vol;



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

/*
    PeriodicPointLattice<Cell<0>> lat0(spec, supercell_spec);
    PeriodicLinkLattice<Cell<0>, Cell<1>> lat1(spec, supercell_spec);
    PeriodicPlaqLattice<Cell<0>, Cell<1>, Cell<2>> lat2(spec, supercell_spec);
    */
    PeriodicVolLattice<Cell<0>, Cell<1>, Cell<2>, Cell<3>> lat(spec, supercell_spec);





   
    lat.print_state(3);

    
    std::random_device rd;
    std::mt19937 gen(rd());
    
    for (int dim=0; dim<4; dim++){
        if (cell_disorder[dim] > 1e-6){
            std::cout << "Deleting " << cell_disorder[dim] * 100 << "% of " << dim << "-cells" << std::endl;
        }
    }

    std::bernoulli_distribution d0(cell_disorder[0]);
    std::bernoulli_distribution d1(cell_disorder[1]);
    std::bernoulli_distribution d2(cell_disorder[2]);
    std::bernoulli_distribution d3(cell_disorder[3]);

    {
        std::vector<Point*> to_delete;
        for (const auto& [_, p] : lat.points) {
            if (d0(gen)) to_delete.push_back(p);
        }
        for (const auto& p : to_delete){
            lat.erase_point(p);
        }
    }
    
    {
        std::vector<Link*> to_delete;
        for (const auto& [_, x] : lat.links) {
            if (d1(gen))  to_delete.push_back(x);
        }
        for (const auto& x: to_delete){
            lat.erase_link(x);
        }
    }
    
    {
        std::vector<Plaq*> to_delete;
        for (const auto& [_, x] : lat.plaqs) {
            if (d2(gen))  to_delete.push_back(x);
        }
        for (const auto& x: to_delete){
            lat.erase_plaq(x);
        }
    }
    
    {
        std::vector<Vol*> to_delete;
        for (const auto& [_, x] : lat.vols) {
            if (d3(gen))  to_delete.push_back(x);
        }
        for (const auto& x: to_delete){
            lat.erase_vol(x);
        }
    }

    


    // debug
    lat.print_state(0);

    save(lat, outpath/(name+".lat.json"));

    return 0;
}


