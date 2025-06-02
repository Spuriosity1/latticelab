#pragma once

#include "cell_geometry.hpp"
#include <concepts>
#include <filesystem>
#include <highfive/highfive.hpp>






namespace CellGeometry {
	template<typename lat_t>
		requires std::derived_from<lat_t, PeriodicAbstractLattice>
	inline herr_t save(const lat_t& lat, const std::filesystem::path& out_path){
		using namespace HighFive;
		File file(out_path, File::Truncate);
		write_data(lat, file);
	}

	inline void write_data(const PeriodicAbstractLattice& lat, HighFive::File& file){
		auto g = file.createGroup("/geometry");
		g.createAttribute<std::vector<int>>("version",{1,0});

		auto dset = file.createDataSet<int64_t>(
				"/geometry/cell_vectors", HighFive::DataSpace(3,3));
		int64_t data[3][3];
		for (int i=0; i<3; i++){ for (int j=0; j<3; j++) {
			data[i][j] = lat.cell_vectors(i,j);
		}}
		dset.write(data);
	}

	template<CellLike<0> Point>
	inline void write_data(const PeriodicPointLattice<Point>& lat, HighFive::File& file){
		write_data(static_cast<const PeriodicAbstractLattice&>(lat), file);

		size_t num_points = lat.points.size();
		auto loc_dset = file.createDataSet<int64_t>(
				"/geometry/point_loc", HighFive::DataSpace(num_points,3)
				);

		auto data = new int64_t[num_points][3];
		for (size_t row=0; row<num_points; row++) {
			for (size_t j=0; j<3; j++) {
				data[row][j] = lat.points[row].position[j];
			}
		}
		loc_dset.write(data);
		delete[] data;
	}

template<
	CellLike<0> Point,
	CellLike<1> Link
	>
	inline void write_data(const PeriodicLinkLattice<Point, Link>& lat, HighFive::File& file){
		write_data(static_cast<const PeriodicPointLattice<Point>&>(lat), file);

		size_t num_links = lat.links.size();

		auto loc_dset = file.createDataSet<int64_t>(
				"/geometry/link_loc", HighFive::DataSpace(num_links,3)
				);
		
		size_t max_boundary_size = 0;
		for (const auto& l : lat.links) {
			max_boundary_size = max(l.boundary.size(), max_boundary_size);
		}

		auto boundary_dset = file.createDataSet<int64_t>(
				"/geometry/link_boundary", HighFive::DataSpace(num_links,2,3)
				);


		auto data = new int64_t[num_links][3];
		auto boundary = new int64_t[num_links][3];
		for (size_t row=0; row<num_links; row++) {
			auto l = lat.links[row];
			for (int j=0; j<3; j++) {
				data[row][j] = l.position[j];
				for (auto& [point_ptr, mult] : lat){
				}
			}
		}
		loc_dset.write(data);
		delete[] data;
		delete[] boundary;
	}


};
