#pragma once

#include <vector>
#include <array>

namespace LWS {
class PolyCurveNetwork;
class BVHNode3D;
class BlockClusterTree;

// Arguments:
// positions - std::vector of size-3 arrays of doubles, each one containing the position of that vertex,
// edges - std::vector of size-2 arrays of size_t, each one with a pair of indices indicating an edge.
// Returns:
// a pointer to a new PolyCurveNetwork object. 
__attribute__ ((visibility ("default")))
PolyCurveNetwork* createCurveNetwork(std::vector<std::array<double, 3>> &positions, std::vector<std::array<size_t, 2>> &edges);

// Arguments:
// curve - pointer to the PolyCurveNetwork that you want the BVH to be built from
// Returns:
// a pointer to a new BVH, which can be used to evaluate energy and gradient with Barnes-Hut.
__attribute__ ((visibility ("default")))
BVHNode3D* createBVHForEnergy(PolyCurveNetwork *curve);

// Arguments:
// curve - pointer to the PolyCurveNetwork that you want the BVH to be built from
// sep - parameter controlling how far away clusters must be to be "admissible." The value we use in the project is 1.0.
// alpha, beta - parameters for the exponents of the tangent-point energy.
// Returns:
// a pointer to the new BlockClusterTree.
__attribute__ ((visibility ("default")))
BlockClusterTree* createBlockClusterTree(PolyCurveNetwork *curve, double sep, double alpha, double beta);

// Arguments:
// tree - a pointer to the BlockClusterTree that you are using to multiply the metric
// vec - a vector of length m (# edges) containing the input data
// output - a vector of length m, where the result of the product will be written
// Returns:
// nothing, but the product will be written into the vector "output".
__attribute__ ((visibility ("default")))
void multiplyMetricWithVector(BlockClusterTree* tree, std::vector<double> &vec, std::vector<double> &output);

// Arguments:
// curve - pointer to the PolyCurveNetwork that you want to compute the energy of
// root - pointer to the BVH you want to use to do the approximation
// alpha, beta - parameters for the exponents of the tangent-point energy.
// Returns:
// the value of the tangent-point energy on the curve. 
__attribute__ ((visibility ("default")))
double evaluateEnergy(PolyCurveNetwork* curve, BVHNode3D *root, double alpha, double beta);

// Arguments:
// curve - pointer to the PolyCurveNetwork that you want to compute the energy of
// root - pointer to the BVH you want to use to do the approximation
// out - a vector of size m x 3, where m is the number of edges. The output will be written here.
// alpha, beta - parameters for the exponents of the tangent-point energy.
// Returns:
// nothing, but the gradient will be written into the vector "out".
__attribute__ ((visibility ("default")))
void evaluateGradient(PolyCurveNetwork* curve, BVHNode3D *root, std::vector<std::array<double, 3>> &out, double alpha, double beta);

}
