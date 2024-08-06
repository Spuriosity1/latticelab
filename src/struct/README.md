# Simplical Complexes

This library models unit cells in the language of __r-chains__, an idea from algebraic geometry. See Weber 2013, https://doi.org/10.1007/JHEP05(2013)157 for a detailed explanation.

The idea here is to formalise the notion of an '$n$-dimensional region of space'. 

The singular simplices $\Sigma_{r}$ are the maximally-simple elements of the cell complex - single points, single links, and individual plaquettes.  Specifically, we care about the situation in $n(+1)$ dimensional lattice gauge theory, in which we naturally think about four kinds of objects:

$$
\{ \rm Volumes \} \xrightarrow{\partial} \{\rm Plaquettes\} \xrightarrow{\partial} \{\rm Links \} \xrightarrow{\partial} \{\rm Vertices \}
$$
where $\partial$ takes the boundary of given object. In the language of algebraic geometry, we identify cells, plaquettes, links and vertices / sites with 3,2,1, and 0-singular simplices respectively.
## $r$-chains

In fact, we should encode more information than just the topology - we also care about _orientation_. Formally, the set $\{\rm Links \}$ contains two copies of every physical location of a cell, each having opposite orientations. Denote $-c$ for an object $c$ with opposite orientation.

A finite exact chain complex is a chain of sets
$$
\begin{align}
&C_n \xrightarrow{\partial} C_{n-1} \xrightarrow{\partial} ... \xrightarrow{\partial} C_0
\end{align}
$$
with maps $\partial_n : C_n \to C_{n-1}$ such that 
1. (chain) Each map, $\partial$, acts like taking the boundary.
2. (exact) $\partial^2 = 0$

The objects residing in the set $C_r$ are the $r$-chains, formal $\mathbb{Z}$-linear combinations of $r$-dimensional singular simplices which we think about as oriented regions of space.
$$C_{r}= \left\{\sum\limits m_{i}\sigma_{i} | m_{i}\in \mathbb{Z} \right\}$$
Note that these contain the orientation information - the $\Sigma_{r}$ are not _a priori_ oriented. However, this coordinate expression for an $r$ chain implicitly assumes that an orientation convention has been chosen.

Now, let $c_{r,}c_{r+1}$ be simplices, (i.e. 'single' plaquettes/links/cells/vertices), and define the _incidence numbers_
$$
\iota_{c_{r}c_{r+1}} = \begin{cases}1 & {c_{r} \subset c_{r+1}}\\ 
-1 & {-c_{r} \subset c_{r+1}}\\ 
0 & \text{else} \end{cases}
$$
which is to say that $\iota_{c_{r},c_{r+1}}$ gives 1 if $c_r$ is contained in $c_{r+1}$ with the correct orientation, -1 if it is contained with the wrong orientation, and zero otherwise. 

One can then see that $\iota$ is a coordinate representation of $\partial$  in terms of the elementary cells $c_r$ :
$$\partial F_{c_{r-1}} = \sum\limits_{c_{r}} \iota_{c_{r-1},c_r}F_{c_{r}}$$

the intermediate sum running over all elementary $r+1$-cells.

Exactness implies the identity
$$\sum\limits_{c_{r}}\iota_{c_{r-1}, c_{r}}\iota_{c_{r}, c_{r+1}} = 0 \Leftrightarrow \partial^2=0$$

Dual (in the Hodge sense) to the boundary operation is the _coboundary operator_ $\delta$,  $\delta_n: C_{n-1} \to C_{n}$, which can be defined in a very sophisticated way using chain cohomology. Instead, we will be physicists about it and choose "coordinates" with which to understand it. 
$$\delta F_{c_{r+1}} = \sum\limits_{c_{r}} \iota_{c_{r},c_{r+1}}F_{c_{r}}$$
The coboundary $\delta$ is therefore understood as the 'transpose' of the boundary. 


