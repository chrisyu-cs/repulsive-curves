# Repulsive Curves
Christopher Yu, Henrik Schumacher, Keenan Crane
ACM Transactions on Graphics 2020 (accepted)

## Quick setup instructions

First, clone the project and all its dependencies:
```
git clone --recursive https://github.com/icethrush/repulsive-curves.git
```

If the recursive flag was not used to clone, then one can also get the dependencies by running:
```
git submodule update --init --recursive
```

From there, the project can be built using CMake.
```
cd repulsive-curves
mkdir build
cd build
cmake ..
make -j4
```
We highly recommend using Clang to build the project. Building with GCC/G++ is possible, but will require a different set of warnings to be suppressed.

The code can then be run:
```
./bin/rcurves_app path/to/scene.txt
```

For best performance, you should make sure that OpenMP is supported on your system.

Note that the file `scene.txt` has a particular format that describes where to find the curve data, as well as what constraints will be used. See `scenes/FORMATS.md` for details.

## Using the project

The important options for manipulating curves are all under the "Curve options" panel in the system. These options are:

+ Run TPE: While checked, the system will run the gradient flow of the tangent-point energy.
+ Output frames: If checked, screenshots of every frame of the gradient flow will be saved as PNG images in the `./frames` directory; note that this is relative to the working directory from which the executable is run.
+ Normalize view: If checked, the objects will be visually rescaled to fit within the camera frame every timestep. This rescaling is purely visual and does not affect the flow.
+ Output OBJs: If checked, OBJs of the curve on every frame will be output to the `./objs` directory.
+ Use Sobolev: If checked, our fractional Sobolev preconditioner will be used. If unchecked, the L2 flow is used instead.
+ Use backprojection: If checked, the system will perform a projection step to enforce hard constraints, correctin for drift. If unchecked, no such step is performed.
+ Use Barnes-Hut: If checked, hierarchical Barnes-Hut approximation is used for energy and gradient evaluations. If unchecked, the energy and gradient are evaluated exactly (and slowly).
+ Use multigrid: If checked, multigrid is used to perform linear solves. If unchecked, dense linear solves are performed.

