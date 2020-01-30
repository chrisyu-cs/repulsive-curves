# Repulsive Curves
Christopher Yu, Henrik Schumacher, Keenan Crane
SIGGRAPH 2020 (submitted)

## Quick setup instructions

First, clone the project and all its dependencies:
```
git clone --recursive https://github.com/icethrush/self-avoiding-flow.git
```

If the recursive flag was not used to clone, then one can also get the dependencies by running:
```
git submodule update --init --recursive
```

From there, the project can be built using CMake.
```
cd self-avoiding-flow
mkdir build
cd build
cmake ..
make -j4
```
We highly recommend using Clang to build the project. Building with GCC/G++ is possible, but will require a different set of warnings to be suppressed.

The code can then be run:
```
./bin/lws path/to/scene.txt
```

Note that the file `scene.txt` has a particular format that describes where to find the curve data, as well as what constraints will be used. See FORMATS.md for details.
