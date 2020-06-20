# Overview

This document describes the data formats and directory structure of the example
files for

   Yu, Schumacher, Crane
   _"Repulsive Curves"_
   **(submitted to SIGGRAPH 2020)**

# Directory Structure

Example files are contained in `self-avoiding/Media`.  Each example is contained
in a subdirectory with an identical structure:

```
self-avoiding/Media/MyExample:
| scene.txt
| curve.obj
| [surface.obj]
| flow/step%05d.obj
| frames/frame%05d.obj
```

The file `scene.txt`, detailed below, contains the "scene" describing the
example.  The file `curve.obj` is the input curve network, as an OBJ file
(also described below).  The optional file `surface.obj` is a standard
OBJ file that describes any surface mesh used to define a potential in
the scene.  The subdirectory `flow/` contains the curve at each time step
of the optimization; each step is an OBJ file with the same format as the
input curve.  Finally, the file `frames/` contains rendered images
corresponding to each of the files in the `flow/` subdirectory.

**IMPORTANT:** The input scene, curve, and optional surface file must be
constructed in such a way that all constraints are already satisfied (or
_very_ nearly satisfied) by the input data.  Otherwise, there is no
guarantee that the solver will find a feasible initial state for the curve,
and may produce garbage output (or terminate without computing the flow).

# Curve files

Curve files are stored as Wavefront OBJ files, but unlike standard surface
mesh files do _not_ contain any triangles.  Instead, they contain only a
list of vertices

```
v x1 y1 z1
v x2 y2 z2
v x3 y3 z3
...
```

as well as a list of edges, given by lines

```
l i1 j1
l i2 j2
l i3 j3
...
```

where `in` and `jn` are 1-based indices into the vertex list.

# Scene files

Each line of a scene files describes either

1. scene geometry,
2. a constraint, or
3. an objective term.

Note that all file paths in `scene.txt` should be relative to the location of `scene.txt`.

## Scene geometry

**Curve.** Each scene file must specify one and only one curve network, via a line

`curve path/to/curve.obj`

where `curve` is a fixed keyword, and a path is given to the curve OBJ file.
These paths are _relative to the path of the scene file_.

## Constraints

Several constraints are supported by solver; each constraint is specified by
an initial keyword followed by some number of parameters (possibly zero).
**Note:** if _no_ constraints are specified, a barycenter constraint will be
added by default (in order to make the flow well-defined).

### Barycenter

```
fix_barycenter
```

Constrains the center of mass of the curve to remain fixed.  No arguments.

### Total length

```
fix_length [scale]
```

Constrains the total length of the curve to remain fixed.  The optional argument
`scale` is a positive number, giving a target scale factor for the total length
relative to the input length; the solver will try to gradually scale up the curve
length to satisfy this constraint.

### Edge length

```
fix_edgelengths [scale]
```

Constrains each individual edge length of the curve to remain fixed.  The
optional argument `scale` is a positive number, giving a target scale factor
for these edge lengths relative to the input edge lengths; the solver will try
to gradually scale up edge lengths to satisfy these constraints.

### Pinned vertices

```
fix_vertex i
```

Fixes the ith vertex of the curve to its input location.  Any number of these
constraints can be specified (on separate lines).  Note that these indices are
_0-indexed_, as opposed to vertices in the OBJ file, which are 1-indexed.

### Tangent constraints

```
fix_tangent i
```

Fixes the tangent of vertex `i`, which is again 0-indexed. The vertex and its neighbors
are free to move, but the discrete tangent at vertex `i` (computed as the normalized average
of neighboring edge tangents) must remain unchanged. Behavior is undefined if vertex `i` has
degree 3 or more. Multiple tangent constraints can be used.

### Convenience pin options

```
fix_special_vertices
```

Equivalent to writing `fix_vertex i` for every vertex `i` with degree not equal to 2,
meaning it is either an endpoint with only 1 edge, or a juncture of at least 3 edges.

```
fix_endpoint_vertices
```

Equivalent to writing `fix_vertex i` for every vertex `i` with degree equal to 1,
meaning it is an endpoint of a curve.

```
fix_special_tangents
```

Causes `fix_special_vertices` or `fix_endpoint_vertices` to add tangent constraints
as well, equivalent to writing `fix_tangent i` for all affected vertices `i`.


### Constrained vertices

It is also possible to constrain vertices to slide around on a surface.  The
surface must first be specified via a line

```
constraint_surface surfaceType [params]
```

Here `constraint_surface` is a fixed keyword, and `surfaceType` is the name of
an implicit surface supported by the solver.  A given surface type may have
several parameters---currently supported names are given in the list below.
Note that multiple constraint surfaces can be specified; the overall constraint
surface will then be the disjoint union of these surfaces.  (The behavior for
overlapping constraint surfaces is undefined.)

To specify which vertices are constrained to the surface, one can either use
the constraints

```
constrain_vertex i
```

for each vertex i that should be constrained, or (for convenience)

```
constrain_all
```

which has no parameters, and just constrains the entire curve to the implicit
surface.

_Available constraint surfaces:_

- `sphere x y z r` --- a sphere with center (x,y,z) and radius r; defaults to
(0, 0, 0) and 1 if the parameters are left unspecified
- `torus x y z r1 r2` --- a torus with center (x,y,z), major radius r1,
and minor radius r2; defaults to (0, 0, 0), 1, and 0.25 if left unspecified
- `doubletorus` --- a topological double torus with fixed origin and size
- `yplane` --- the infinite plane y = 0

## Objectives

Like constraints, each objective is specified via a single line with a keyword
followed by a list of parameters.

### Curve repulsion

```
repel_curve [alpha] [beta] [weight]
```

Adds a term that prevents curve self-collision.  The optional parameters `alpha` and
`beta` control the strength of repulsion.  (See paper for details.)  Reasonable
defaults will be used if alpha and beta are not specified.

The optional parameter `weight` controls the relative strength of this term.
The default weight for all objective terms is 1.

### Surface repulsion

```
repel_surface path/to/surface.obj [weight]
```

Adds a term that prevents the curve from colliding with a fixed surface, which is
the one found at the given path (again relative to the scene file).

The optional parameter `weight` controls the relative strength of this term.

### Plane repulsion

```
repel_plane c_x c_y c_z n_x n_y n_z [weight]
```

Adds a term that prevents the curve from passing through an infinite plane passing
through the point (c_x, c_y, c_z), with normal (n_x, n_y, n_z).

The optional parameter `weight` controls the relative strength of this term.

### Total length

```
optimize_length [weight]
```

Adds a term that penalizes the total curve length.

The optional parameter `weight` controls the relative strength of this term.
(Note that the weight can be negative.)

### Length difference

```
optimize_length_diffs [weight]
```

Adds a term that penalizes differences in lengths of adjacent edges. This is
meant to encourage uniformity of lengths without enforcing a hard constraint.

The optional parameter `weight` controls the relative strength of this term.
(Note that the weight can be negative.)

### Vector field

```
optimize_field fieldType [weight]
```

Adds a term that encourages the curve to be tangent to a fixed vector field.
The parameter `fieldType` is a string that specifies which field should be
used.  These strings come from a predefined list of available vector fields:

- `circular` --- A vector field that circulates around the origin in the XZ plane
- `interesting` --- A vector field on the sphere that involves some more singularities

## Other options

For convenience, we allow a few other settings to be specified in the scene file.
These settings are all optional.

### Iteration limit

```
iteration_limit n
```

Automatically halts the flow after `n` steps. If `n` is 0, or if this option is
omitted, then no limit is used.


### Subdivision limit

```
subdivide_limit n
```

Allows the system to subdivide the curve automatically once the average edge length
has exceeded twice the original average edge length. This would occur, for instance,
as a result of using the optional scale arguments to fix_edgelengths or fix_length.

The argument `n` specifies how many times the curve may be subdivided.
By default, the limit is 0, meaning no subdivision will occur.


## Example scene file

Here's a small example of a scene file that sets up a repulsive curve subject to
fixed barycenter and edge length constraints, which is also affected by a surface penalty:

```
curve curve.obj
repel_curve 3 6
fix_barycenter
fix_edgelengths
repel_surface surface.obj
```

