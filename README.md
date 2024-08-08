# LatticeLab (working title)

Any condensed matter physicist will be familiar with the pain of setting up 
and using an indexing scheme for points on a lattice. This becomes even more
complicated when one wishes to also consider more complicated geometrical
objects - links, elementary plaquettes, and elementary cells.

This project aims to provide a **single, universal** solution to this problem,
defining a flexible, templated class `Lattice` which handles all the
geometric pain for you.


# Usage Examples

## 3D Cubic lattice
```c++

	UnitCellSpecifier cell(imat33_t::from_cols({
				{2,0,0},
				{0,2,0},
				{0,0,2}
				}));
	PointSpec pointspec;
	pointspec.position = {0,0,0};	
	cell.add_point(pointspec);

	LinkSpec linkspec;
	linkspec.position = {1,0,0};
	linkspec.boundary = { {1, {-1,0,0}}, {-1, {1,0,0}} };
	cell.add_link(linkspec);

	linkspec.position = {0,1,0};
	linkspec.boundary = { {1, {0,-1,0}}, {-1, {0,1,0}} };
	cell.add_link(linkspec);

	linkspec.position = {0,0,1};
	linkspec.boundary = { {1, {0,0,-1}}, {-1, {0,0,1}} };
	cell.add_link(linkspec);

	PlaqSpec plaqspec;
	plaqspec.position = {0,1,1};
	plaqspec.boundary = {
		{1, {1,0,0}}, {1, {2,1,0}}, {-1, {1,2,0}}, {-1, {0,1,0}}
	};
	cell.add_plaq(plaqspec);

	plaqspec.position = {1,0,1};
	plaqspec.boundary = {
		{1, {0,1,0}}, {1, {0,2,1}}, {-1, {0,1,2}}, {-1, {0,0,1}}
	};
	cell.add_plaq(plaqspec);

	plaqspec.position = {1,1,0};
	plaqspec.boundary = {
		{1, {0,0,1}}, {1, {1,0,2}}, {-1, {2,0,1}}, {-1, {1,0,0}}
	};
	cell.add_plaq(plaqspec);

	VolSpec volspec;
	volspec.position = {1,1,1};
	volspec.boundary = {
		{1, {0,1,1}}, {-1, {2,1,1}},
		{1, {1,0,1}}, {-1, {1,2,1}},
		{1, {1,1,0}}, {-1, {1,1,2}},
	};	
	cell.add_vol(volspec);
```
