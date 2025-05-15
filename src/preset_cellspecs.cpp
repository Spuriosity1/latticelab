#include "preset_cellspecs.hpp"
#include "rationalmath.hpp"


namespace CellGeometry {
namespace PrimitiveSpecifiers {


void Diamond::setup_points(){
	const std::vector<ipos_t> point_positions = {{0, 0, 0}, {2, 2, 2}};
	// set up the points
	CellGeometry::PointSpec pointspec;
	for (const auto &x : point_positions) {
		pointspec.position = x;
		this->add_point(pointspec);
	}
}

void Diamond::setup_links(){
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

void Diamond::setup_plaqs(){
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

void Diamond::setup_vols(){
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


void Cubic::setup_points(){
	PointSpec pointspec;
	pointspec.position = {0, 0, 0};
	this->add_point(pointspec);
}

void Cubic::setup_links(){
	LinkSpec linkspec;
	linkspec.position = {1,0,0};
	linkspec.boundary = { {1, {-1,0,0}}, {-1, {1,0,0}} };
	this->add_link(linkspec);

	linkspec.position = {0,1,0};
	linkspec.boundary = { {1, {0,-1,0}}, {-1, {0,1,0}} };
	this->add_link(linkspec);

	linkspec.position = {0,0,1};
	linkspec.boundary = { {1, {0,0,-1}}, {-1, {0,0,1}} };
	this->add_link(linkspec);
}

void Cubic::setup_plaqs(){
	PlaqSpec plaqspec;
	plaqspec.position = {0,1,1};
	plaqspec.boundary = {
		{1, {0,0,-1}}, {1, {0,1,0}}, {-1, {0,0,-1}}, {-1, {0,-1,0}}
	};
	this->add_plaq(plaqspec);

	plaqspec.position = {1,0,1};
	plaqspec.boundary = {
		{1, {-1,0,0}}, {1, {0,0,1}}, {-1, {-1,0,0}}, {-1, {0,0,-1}}
	};
	this->add_plaq(plaqspec);

	plaqspec.position = {1,1,0};
	plaqspec.boundary = {
		{1, {0,-1,0}}, {1, {1,0,0}}, {-1, {0,-1,0}}, {-1, {-1,0,0}}
	};
	this->add_plaq(plaqspec);
}

void Cubic::setup_vols(){
	VolSpec volspec;
	volspec.position = {1,1,1};
	volspec.boundary = {
		{1, {1,0,0}}, {-1, {-1,0,0}},
		{1, {0,1,0}}, {-1, {0,-1,0}},
		{1, {0,0,1}}, {-1, {0,0,-1}}
	};	
	this->add_vol(volspec);
}


const UnitCellSpecifier DiamondSpec() { return Diamond(); }
const UnitCellSpecifier CubicSpec()   { return Cubic(); }

};
};


