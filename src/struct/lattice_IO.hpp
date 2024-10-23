#pragma once

#include "cell_geometry.hpp"
#include "chain.hpp"
#include "nlohmann/json_fwd.hpp"
#include <concepts>
#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>


namespace CellGeometry {
	template<typename lat_t>
		requires std::derived_from<lat_t, PeriodicAbstractLattice>
	inline bool save(const lat_t& lat, const std::filesystem::path& out_path){
		using namespace nlohmann;
		json j = {};
		write_data(lat, j);
		std::ofstream of(out_path);
		if (of.is_open()){
			of << j;
			of.close();
			return true;
		}
		return false;
	}

	inline void write_data(const PeriodicAbstractLattice& lat, nlohmann::json& j){
		j["cell_vectors"] = lat.cell_vectors;
		j["index_cell_vectors"] = lat.index_cell_vectors;
		j["primitive_cell_vectors"] = lat.primitive_spec.lattice_vectors;
	}

	template<int order, CellLike<order> T>
	requires (order < 3)
	inline void store_coboundary(const T& obj, nlohmann::json& pt){
		pt["coboundary"] = {};
		for (const auto& [cellptr, mult] : obj.coboundary){
			pt["coboundary"].push_back({
					{"pos",cellptr->position},
					{"mult", mult}
					});
		}
	}

	template<int order, CellLike<order> T>
	requires (order > 0)
	inline void store_boundary(const T& obj, nlohmann::json& pt){
		pt["boundary"] = {};
		for (const auto& [cellptr, mult] : obj.boundary){
			
			pt["boundary"].push_back({
					{"pos",cellptr->position},
					{"mult", mult}
					});
		}
	}

	template<CellLike<0> Point>
	inline void write_data(const PeriodicPointLattice<Point>& lat, nlohmann::json& j){
		write_data(static_cast<const PeriodicAbstractLattice&>(lat), j);
		j["points"] = {};
		for (const auto& point : lat.points){
			nlohmann::json pt = {{"pos", point.position}};	
			store_coboundary<0>(point, pt);
			j["points"].push_back(pt);
		}
	}

	template<
		CellLike<0> Point,
		CellLike<1> Link
	>
	inline void write_data(const PeriodicLinkLattice<Point, Link>& lat, nlohmann::json& j){
		write_data(static_cast<const PeriodicPointLattice<Point>&>(lat), j);
		j["links"] = {};
		for (const auto& link : lat.links){
			nlohmann::json ln = {};
			ln["pos"] = link.position;
			store_coboundary<1>(link, ln);
			store_boundary<1>(link, ln);
			j["links"].push_back(ln);
		}
	}

	template<
		CellLike<0> Point,
		CellLike<1> Link,
		CellLike<2> Plaq
	>
	inline void write_data(const PeriodicPlaqLattice<Point, Link, Plaq>& lat,
			nlohmann::json& j){
		write_data(static_cast<const PeriodicLinkLattice<Point, Link>&>(lat), j);
		j["plaqs"] = {};
		for (const auto& plaq : lat.plaqs){
			nlohmann::json ln = {{"pos", plaq.position}};
			store_coboundary<2>(plaq, ln);
			store_boundary<2>(plaq, ln);
			j["plaqs"].push_back(plaq);
		}
	}

	template<
		CellLike<0> Point,
		CellLike<1> Link,
		CellLike<2> Plaq,
		CellLike<3> Vol
	>
	inline void write_data(const PeriodicVolLattice<Point, Link, Plaq, Vol>& lat,
			nlohmann::json& j){
		write_data(static_cast<const PeriodicPlaqLattice<Point, Link, Plaq>&>(lat), j);
		j["volumes"] = {};
		for (const auto& vol : lat.vols){
			nlohmann::json v = {{"pos", vol.position}};
			store_boundary<3>(vol, v);
			j["volumes"].push_back(v);
		}
	}


};
