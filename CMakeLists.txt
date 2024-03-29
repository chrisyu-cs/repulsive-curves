cmake_minimum_required(VERSION 3.10.0)

project(self_avoiding)

### Configure output locations
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_BUILD_TYPE "Release")

# Print the build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release" FORCE)
endif()
message(STATUS "cmake build type: ${CMAKE_BUILD_TYPE}")

### Configure the compiler
# This is a basic, decent setup that should do something sane on most compilers

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

  # using Clang (linux or apple) or GCC
  message("Using clang/gcc compiler flags")
  SET(BASE_CXX_FLAGS "-std=c++11 -Wall -Wextra -g3 -fPIC -fopenmp")
  SET(DISABLED_WARNINGS " -Wno-unused-parameter -Wno-return-std-move -Wno-deprecated-copy -Wno-unused-variable -Wno-unused-function -Wno-deprecated-declarations -Wno-missing-braces -Wno-unused-private-field")
  SET(TRACE_INCLUDES " -H -Wno-error=unused-command-line-argument")

  if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    message("Setting clang-specific options")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -ferror-limit=3 -fcolor-diagnostics")
    SET(CMAKE_CXX_FLAGS_DEBUG          "-fsanitize=address -fno-limit-debug-info")
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message("Setting gcc-specific options")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -fmax-errors=5")
    SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} -Wno-maybe-uninitialized -Wno-format-zero-length -Wno-unused-but-set-parameter -Wno-unused-but-set-variable")
  endif()
  
  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  SET(CMAKE_CXX_FLAGS_DEBUG          "-g3")
  SET(CMAKE_CXX_FLAGS_RELEASE        "-O3 -march=native -DNDEBUG")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
  message("Using Visual Studio compiler flags")
  set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} /W4")
  set(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} /MP") # parallel build
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4267\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4244\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4305\"")  # ignore truncation on initialization
  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")

  add_definitions(/D "_CRT_SECURE_NO_WARNINGS")
  add_definitions (-DNOMINMAX)
else()
  # unrecognized
  message( FATAL_ERROR "Unrecognized compiler [${CMAKE_CXX_COMPILER_ID}]" )
endif()


# == Deps
add_subdirectory(deps/geometry-central)
add_subdirectory(deps/polyscope)

# == Build our project stuff

set(SRCS
  src/circle_search.cpp
  src/curve_io.cpp
  src/extra_potentials.cpp
  src/implicit_surface.cpp
  src/lws_options.cpp
  src/poly_curve_network.cpp
  src/scene_file.cpp
  src/sobo_slobo.cpp
  src/tpe_energy_sc.cpp
  src/tpe_flow_sc.cpp
  src/utils.cpp
  src/vert_jacobian.cpp
  src/applications/pathplanning.cpp
  src/flow/constraint_functions.cpp
  src/flow/gradient_constraint_enum.cpp
  src/marchingcubes/CIsoSurface.cpp
  src/marchingcubes/Vectors.cpp
  src/obstacles/obstacle.cpp
  src/obstacles/mesh_obstacle.cpp
  src/obstacles/plane_obstacle.cpp
  src/obstacles/sphere_obstacle.cpp
  src/product/block_cluster_tree.cpp
  src/product/dense_matrix.cpp
  src/product/test_matrices.cpp
  src/spatial/spatial_tree.cpp
  src/spatial/tpe_bvh.cpp
  src/spatial/vertex_body.cpp
	# add any other source files here
)

# add_executable(rcurves_app src/lws_app.cpp "${SRCS}")
# target_include_directories(rcurves_app PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
# target_include_directories(rcurves_app PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/libgmultigrid/include/")
# target_link_libraries(rcurves_app geometry-central polyscope)

add_library(rcurves STATIC "${SRCS}")
target_include_directories(rcurves PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/")
target_include_directories(rcurves PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/deps/libgmultigrid/include/")
target_link_libraries(rcurves geometry-central polyscope)
target_compile_options(rcurves PUBLIC -fvisibility=hidden)

add_executable(rcurves_app src/lws_app.cpp)
target_link_libraries(rcurves_app rcurves)

add_library(rcurves_shared SHARED src/export/mvproduct.cpp)
target_link_libraries(rcurves_shared rcurves)
target_compile_options(rcurves_shared PUBLIC -fvisibility=default)

