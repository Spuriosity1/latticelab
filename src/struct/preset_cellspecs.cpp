#include "preset_cellspecs.hpp"


namespace CellGeometry {
namespace PrimitiveSpecifers {


	void Pyrochlore::setup_points(){
		const std::vector<ipos_t> point_positions = {{0, 0, 0}, {2, 2, 2}};
		// set up the points
		CellGeometry::PointSpec pointspec;
		for (const auto &x : point_positions) {
			pointspec.position = x;
			this->add_point(pointspec);
		}
	}

	void Pyrochlore::setup_links(){
		// set up the links
		const std::vector<ipos_t> link_positions = {
			{1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}};
		CellGeometry::LinkSpec linkspec;
		for (const auto &x : link_positions) {
			linkspec.position = x;
			linkspec.boundary = {{1, x}, {-1, -x}};
			this->add_link(linkspec);
		}
	}

	void Pyrochlore::setup_plaqs(){
		// set up the plaqs	
		const ipos_t plaq_positions[4]  = {
			{-3,-3,-3}, {-3,-1,-1}, {-1,-3,-1}, {-1,-1,-3}
			//{-1,-1,-1}, {-1, 1, 1}, { 1,-1, 1}, { 1, 1,-1}
		};
		const std::vector<VectorSignPair> plaq_boundaries[4] = {
			{
				{ 1,{0, -2, 2}},    {-1, {2, -2, 0}},
				{ 1,{2, 0, -2}},    {-1, {0, 2, -2}}, 
				{ 1,{-2, 2, 0}},    {-1, {-2, 0, 2}}
			},
			{
				{ 1,{ 0, 2,-2}},    {-1,{ 2, 2, 0}},
				{ 1,{ 2, 0, 2}},    {-1,{ 0,-2, 2}},
				{ 1,{-2,-2, 0}},    {-1,{-2, 0,-2}}
			},
			{
				{ 1,{ 0,-2,-2}},	{-1,{-2,-2, 0}},
				{ 1,{-2, 0, 2}},	{-1,{ 0, 2, 2}},
				{ 1,{ 2, 2, 0}},	{-1,{ 2, 0,-2}}
			},
			{
				{ 1,{ 0, 2, 2}},	{-1,{-2, 2, 0}},
				{ 1,{-2, 0,-2}},	{-1,{ 0,-2,-2}},
				{ 1,{ 2,-2, 0}},	{-1,{ 2, 0, 2}}
			}
		};

		CellGeometry::PlaqSpec plaqspec;
		for (int mu=0; mu<4; mu++){
			plaqspec.position = plaq_positions[mu];
			plaqspec.boundary = plaq_boundaries[mu];
			this->add_plaq(plaqspec);
		}
	}

	void Pyrochlore::setup_vols(){
		const std::vector<ipos_t> vol_positions = {
			{4,4,4},{6,6,6}
		};

		CellGeometry::VolSpec volspec;
		volspec.position = vol_positions[0];
		volspec.boundary = {
			{1, { 1, 1, 1}}, {1, { 1,-1,-1}},
			{1, {-1, 1,-1}}, {1, {-1,-1, 1}},
		};	
		this->add_vol(volspec);
		volspec.position = vol_positions[1];
		volspec.boundary = {
			{-1, {-1,-1,-1}}, {-1, {-1, 1, 1}},
			{-1, { 1,-1, 1}}, {-1, { 1, 1,-1}},
		};
		this->add_vol(volspec);
	}





};
};


