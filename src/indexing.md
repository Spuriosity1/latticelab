# Index convention

The problems to solve:
1. An efficient location-to-index method.
2. A generalisable and scalable indexing convention.

The approach:
- The primitive unit cell has some basis vectors $[a_1~a_2~a_3] =: a$
- The enlarged unit cell has some basis vectors $[A_1~A_2~A_3] =: A$, with the constraint that each $A_i = [a_1~a_2~a_3] * Z_i$ for $Z_i \in \mathbb{Z}^3 \setminus 0$.
- 

Given $r\in \mathbb{Z}^3$, wish to solve the generalised diophantine equation
$$
R = A n + a x + r 
$$
where $x, n \in \mathbb{Z}^3$ and $r$ is within a primitive cell.
We don't care bout the value of $n$, since $r + A_i \equiv r$.
Clearly, this is not unique - any translation $x\mapsto x + z_i$ 
yields another solution with a different $n$:

$$
R = a (Z n + x) + r
$$

and so we must consider $x$ solutions "modulo Z". This is awkward to implement directly.


Instead, at initialisation we compute the Smith normal form -
$$
L Z W = D
$$
Where $L$ and $W$ are invertible integer valued matrices, and $D$ a diagonal integer matrix.

This tells us that there are exactly $D_{11}D_{22}D_{33}$ primitive unit cells in the supercell.

$$
R = a L^{-1} (D W^{-1} n + L x) + r
$$

As an equivalence class on the synthetic, infinite lattice this is
$$ 
\{ a L^{-1} (D W^{-1} n + L x) + b | n \int \mathbb{Z}^3 \} 
= \{ a L^{-1} (D m + L x) + b | m \int \mathbb{Z}^3 \} 
$$
(equality follows from invertibility of $W$)

We can therefore straightforwardly find a unique solution for Lx by considering entries modulo D.

It is convenient to define 
 - $Lx \equiv y \in \mathbb{Z}^3$,
 - $a L^{-1} \equiv b$

 to get the simpler $R = b(D m + y) + r$
 where $y_i \in \{0, ..., D_{ii} - 1\}$

# WHAT WE DO HERE

The user specifies a primitive cell, $a$=`specified_primitive.cell_vectors`, 
with its cell vectors understood as columns, as well as a supercell $Z$ such that
$$A = a Z\,.$$

We want a way to number the copies of cells, i.e. Bravais points, consistently and uniquely.
This means solving $$X = a n$$ modulo $A n \equiv 0$. This is clunky -- in general, it's not just $n=[0,1,2,Lx] x [0,1,2,Ly]$, but something more confusing. For example, think of the equivalence classes when $a = [1 0; 0 1], Z = [2 -2; 2 2]$.

Instead, for indexing we change the definition of the primitive cell to something aligned with $A$. We define the alternative vectors

$a' = a L^-1 $

where $L, W$ are invertible int matrices such that $LZW = D$, a diagonal matrix. 

so we solve $X = a' n$ for $n$ modulo $a' D n = 0$


# APPROACHES

## Float math
Suppose a, A are float-value. Writing the remainder $r$ in the $a$ basis,
$r = a L^{-1} L q, q\in \mathbb{R}$ (where $q$ is from from a rescaled cell of unit volume).

One can then solve
$ (a L^{-1}) Y = R $ and then consider the three entries of $Y = D m + L x + L q$ modulo $D$.

```julia
Y = modulo.(Y, D) # Entrywise float valued modulo wraps all entries to [0, D[i])
Lx = floor.(Y) # Lx is integer valued, and we choose the convention such that Lb is entrywise in [0,1]
Lq = Y - Lx
```

The trouble: Floating-point roundoff error will cause problems for any sites located on the zone boundary.

## Integer math
Instead, we demand that all geometric objects have _integer-valued_ positions. (This is not as restrictive as one might think- if the exact positions are important, one can add by hand a float-valued "realpos" field)

Compute yet another SNF, this time for $a L^{-1}$:

$U (a L^{-1}) V = P$ for diagonal $P$.

```julia
# R = a[Z n + x] + r
# U R = P V^{-1} [D m + L x ] + U r
UR = (U * R)
Y = V *  fldiv.( UR, diag[P] ) # Floor divide to snap to correct unit cell
Lx = modulo.(Y, D) # Entrywise float valued modulo wraps all entries to [0, D[i])
# sublattice info is contained here
r = Uinv * modulo,(UR, diag[P]) 
```


