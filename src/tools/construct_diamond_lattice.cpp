#include "cell_geometry.hpp"
#include "argparse/argparse.hpp"
#include "chain.hpp"
#include "lattice_IO.hpp"
#include "preset_cellspecs.hpp"
#include <UnitCellSpecifier.hpp>
#include <algorithm>
#include <random>
#include <sstream>
#include <stack>
#include <stdexcept>
#include <string>
#include <vector>
/**
 * dmndlat, a utility for generating diamond lattices precisely
 */

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
        const std::string& Z3_s
        ){
    // hashing the arguments 
    std::ostringstream name;
    name << "Z1="+Z1_s+";Z2="+Z2_s+";Z3="+Z3_s+";";
    auto s = name.str();
    std::replace(s.begin(), s.end(), ' ', ',');  // replace space by commas
    return s;
}

inline std::string hash_deletions(
    const std::array<std::vector<int>, 4>& del_list
        ){
    std::stringstream name_ss;
    for (int order=0; order<4; order++){
        auto& dlist = del_list[order];
        if (dlist.size() > 0){
            // hash the setup for naming
            name_ss << "d"<<order<<"=";
            auto sep = "";
            for (const auto& j : dlist){
                name_ss << sep << j;
                sep = ",";
            }
            name_ss << ";";
        }
    }
    return name_ss.str();
}


struct Point : public Cell<0> {
};

struct Link : public Cell<1> {
    bool visited = false;
};

struct Plaq : public Cell<2> {
};

struct Vol : public Cell<3> {
    bool visited = false;
};

typedef PeriodicVolLattice<Point, Link, Plaq, Vol> Lattice;

using namespace std;

struct bfs_node {
    const Cell<0>* point;
    Chain<1> path;
};


inline std::vector<Chain<1>> find_paths_neighbours(Lattice& lat, Point* origin, Point* finish, unsigned len){
    // finds all chains of length 'len' connecting p0 to p1
    std::vector<Chain<1>> res;
    for (auto& l : lat.get_links()){
        l.visited = false;
    }
    std::stack<bfs_node> to_visit;
    to_visit.push(bfs_node(origin, Chain<1>()));

    Chain<0> expected_boundary;
    expected_boundary[origin] = -1;
    expected_boundary[finish] = 1;

    while (!to_visit.empty()){
        auto curr = to_visit.top();
        to_visit.pop();
        if (curr.path.size() == len){
            if (curr.point == finish) {
                assert( d(curr.path) == expected_boundary );
                res.push_back(curr.path);
            }
        }
        for (const auto& [l, m] : curr.point->coboundary){
            const Link* ln = static_cast<const Link*>(l);
            if (ln->visited) continue;
            Chain<1> this_ln; this_ln[l]=m;
            for (auto& [p2, n] : d(this_ln)){
                if (p2 != curr.point){ 
                    to_visit.push(bfs_node(p2,curr.path + this_ln));
                }
            }
        }
    }
    return res;
}


inline std::vector<std::set<Vol*>> find_connected_components(Lattice& lat){
    // greedy algorithm -- traverse via all boundaries recursively
    std::stack<Vol*> vols;

    for (auto& [_,v] : lat.vols){
        v->visited = false;
    }

    std::vector<std::set<Vol*>> retval;

    for (auto& [_,v] : lat.vols){
        if (!v->visited){
            retval.push_back({});
            // start a DFS
            vols.push(v);
            while(!vols.empty()){
                auto curr =vols.top();
                retval.back().insert(curr);
                vols.pop();
                curr->visited = true;
                for(auto v2 :  get_neighbours<Vol>(curr)){
                    if (! v2->visited ){
                        vols.push(v2);
                    }
                }
            }
        }
    }

    return retval;
}


int main (int argc, const char *argv[]) {

    argparse::ArgumentParser prog("dmndlat");

    prog.add_argument("outpath")
        .help("Path to store output");
    prog.add_argument("Z1")
        .help("First lattice vector in primitive units (as space separated string)");
    prog.add_argument("Z2")
        .help("First lattice vector in primitive units (as space separated string)");
    prog.add_argument("Z3")
        .help("First lattice vector in primitive units (as space separated string)");

    
    prog.add_argument("--verbosity")
        .scan<'i', int>()
        .default_value(0);

    std::array<std::vector<int>, 4> del_list;
    prog.add_argument("--d0")
        .help("Specifies index of 0-forms to delete")
        .default_value<std::vector<int>>({})
        .append()
        .store_into(del_list[0]);
    prog.add_argument("--d1")
        .help("Specifies index of 1-forms to delete")
        .default_value<std::vector<int>>({})
        .append()
        .store_into(del_list[1]);
    prog.add_argument("--d2")
        .help("Specifies index of 2-forms to delete")
        .default_value<std::vector<int>>({})
        .append()
        .store_into(del_list[2]);
    prog.add_argument("--d3")
        .help("Specifies index of 3-forms to delete")
        .default_value<std::vector<int>>({})
        .append()
        .store_into(del_list[3]);
        

    try {
        prog.parse_args(argc, argv);
    } catch (const std::exception& err){
        cerr << err.what() << endl;
        cerr << prog;
        std::exit(1);
    }



    std::filesystem::path outpath(prog.get<string>("outpath"));
    string Z1_s = prog.get<string>("Z1");
    string Z2_s = prog.get<string>("Z2");
    string Z3_s = prog.get<string>("Z3");

    // parse L1 L2 L3
    imat33_t supercell_spec;
    parse_supercell_spec(supercell_spec, Z1_s, Z2_s, Z3_s);

    std::cout<<"Constructing supercell of dimensions \n"<<supercell_spec<<std::endl;

    const PrimitiveSpecifers::Diamond spec;
    auto name = hash_parameters(Z1_s, Z2_s, Z3_s);

    Lattice lat(spec, supercell_spec);


    std::stringstream name_ss;
    for (int order=0; order<4; order++){
        auto& dlist = del_list[order];

        // Sort all of these so we don't muck up any indexing
        std::sort(dlist.begin(), dlist.end());
        // silently remove any duplicates
        dlist.erase( unique( dlist.begin(), dlist.end() ), dlist.end() );
        
    }
    name += hash_deletions(del_list);

    // Perform the deletions starting from vols (which have no dise effects)
    // in order to make sure that any deleted things still exist
    try {
        for (auto i=del_list[3].rbegin(); i<del_list[3].rend(); i++){
            lat.erase_vol(lat.vols.at(*i)); 
        }
        for (auto i=del_list[2].rbegin(); i<del_list[2].rend(); i++){
            lat.erase_plaq(lat.plaqs.at(*i));
        }
        for (auto i=del_list[1].rbegin(); i<del_list[1].rend(); i++){
            lat.erase_link(lat.links.at(*i));
        }
        for (auto i=del_list[0].rbegin(); i<del_list[0].rend(); i++){
            lat.erase_point(lat.points.at(*i));
        }
    } catch (const std::out_of_range& e){
        cout << "Tried to delete object at nonexistent index " <<std::endl;
        return 1;
    }

    auto verbosity = prog.get<int>("--verbosity");
    lat.print_state(verbosity);

    save(lat, outpath/(name+".lat.json"));

    return 0;
}


