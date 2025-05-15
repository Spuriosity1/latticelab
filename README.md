# Latlib

[![Build and Test](https://github.com/Spuriosity1/latticelab/actions/workflows/main.yml/badge.svg)](https://github.com/Spuriosity1/latticelab/actions/workflows/main.yml)

Any condensed matter physicist will be familiar with the pain of setting up 
and using a good indexing scheme for points on a periodic lattice. It's tough to get boundary conditions right, and even tougher to consider *disorder* (i.e. deletions) in a way that doesn't feel, well, hacky.

This becomes even more complicated when one wishes to also consider higher-dimensional geometrical objects than simple points, such as links, elementary plaquettes, and elementary cells.

This project aims to provide a **single, universal** solution to this problem,
defining a flexible, templated class `Lattice` which allows the user to quickly and easily define general chain complexes in 3D (or lower).


# INSTALLING

```bash
git clone https://github.com/Spuriosity1/latticelab.git && cd latticelab
meson setup build 
ninja -C build install
```

If you can't touch `/usr/local`, instead do
```bash
meson setup build -Dprefix="/your/install/prefix..."
```

# Philosophy

The `UnitCellSpecifier` is essentially a plan for how to build a regular lattice.


# Usage Examples

## Basic Geometry

### 3D Cubic lattice
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

## Extending the Cell classes

Index convention: `J =  sl*L2*L1*L0  + (i2*L1 + i1)*L0 + i0`

Pros (+) and cons (-):
+ + Easy vectorisation for adding sublattices
+ + Convenient for calculating SL-dependent FT
+ - physically nearby spins have very distant memory addresses (unavoidable)

```c++


/*
struct spin: public Cell<1> {
	double heis[3];

};

double plaq_ring(const Cell<2>& plaq){
	double real = 1;
	double imag = 0;
	double tmp;
	for (auto& [mul, link] : plaq.boundary){
		// real + i imag *= heis[0] + i mul * heis[1]
		tmp = real * heis[0] - imag*mul*heis[1];
		imag = real *mul*heis[1] + imag*heis[0];
		real = tmp;	
	}
	return retval.real();
}
```


