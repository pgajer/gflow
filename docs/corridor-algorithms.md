# Graph corridor constructions

Let a weighted graph have nonnegative edge lengths and shortest-path distance
`d`. For endpoints `u` and `v` in the same connected component, let
`gamma` be a shortest path and `L = d(u, v)` its length. The constructions below
define different neighborhoods of that path; their parameters have different
geometric meanings.

## Segment detour corridor

For each consecutive pair `(a, b)` on `gamma`, retain a vertex `x` when

$$
d(x,a) + d(x,b) \le d(a,b)(1+\mathrm{rel.tol}) + \mathrm{abs.tol}.
$$

Take the union over path segments and include the path vertices. Here
`rel.tol` is dimensionless and `abs.tol` has the units of edge length. Small
tolerances restrict the corridor to vertices reached with little extra travel
between the segment endpoints; larger tolerances admit alternative routes and
side branches. This is a detour criterion, rather than a fixed-radius tube.

## Path-relative tube intersected with an endpoint lens

Choose a nonnegative, dimensionless radius fraction `rho` and set `r = rho * L`. Let
`d_gamma(x)` be the shortest distance from `x` to any vertex on `gamma`.
The tube-and-lens construction is

$$
C(u,v) = \{x: d_\gamma(x)\le r,\ d(u,x)\le L,\ d(v,x)\le L\}.
$$

It can be computed with two endpoint Dijkstra searches truncated at `L` and
one multi-source Dijkstra search from all path vertices truncated at `r`.
The endpoint searches must extend to `L`, not `r`: a small radius would
otherwise exclude path vertices far from either endpoint. A vertex absent
from either endpoint search is excluded. Non-strict endpoint inequalities
retain `u` and `v`; strict inequalities exclude them unless they are added
explicitly.

The endpoint lens limits the tube in terms of distance to both endpoints. It
can still include side pockets and cul-de-sacs close to the path, so membership
does not imply that a vertex lies on a nearly shortest route between `u` and
`v`.

## Additional geodesic-excess restriction

A more restrictive corridor also requires

$$
d(u,x) + d(x,v) - L \le \tau,
$$

where `tau` is a nonnegative allowed excess length, for example a chosen multiple of
`r`. This condition bounds the extra travel through `x` between the two
endpoints. It distinguishes a broad tube-and-lens neighborhood from a
near-geodesic corridor while using the same three distance searches.

These definitions describe the geometric choices. They do not establish
biological arm membership or empirical accuracy for a particular dataset.
