# gc-polyscope-project-template
A template project to get started with geometry-central and Polyscope.

TODO instructions for IDEs and Windows

TODO this repo is set up as a template, but it seems the submodules don't get copied to the new project?

### Get the code
Clone the project 
```
git clone --recursive https://github.com/nmwsharp/gc-polyscope-project-template.git
```

### Build the code

Configure (with cmake) and compile
```
cd gc-polyscope-project-template
mkdir build
cd build
cmake ..
make -j6
```

### Run the code
```
./bin/gc_project /path/to/a/mesh
```
