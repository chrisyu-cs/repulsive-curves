#include "export/mvproduct.h"

#include <Eigen/Core>
#include "geometrycentral/utilities/vector3.h"
#include "poly_curve_network.h"
#include "spatial/tpe_bvh.h"
#include "product/block_cluster_tree.h"

using namespace geometrycentral;

namespace LWS
{
PolyCurveNetwork *createCurveNetwork(std::vector<std::array<double, 3>> &positions, std::vector<std::array<size_t, 2>> &edges)
{
    std::vector<Vector3> vectors(positions.size());
    for (size_t i = 0; i < positions.size(); i++)
    {
        vectors[i] = Vector3{positions[i][0], positions[i][1], positions[i][2]};
    }
    return new PolyCurveNetwork(vectors, edges);
}

BVHNode3D *createBVHForEnergy(PolyCurveNetwork *curves)
{
    return CreateBVHFromCurve(curves);
}

BlockClusterTree *createBlockClusterTree(PolyCurveNetwork *curves, double sep, double alpha, double beta)
{
    BVHNode3D *edgeBVH = CreateEdgeBVHFromCurve(curves);
    BlockClusterTree *tree = new BlockClusterTree(curves, edgeBVH, sep, alpha, beta, 0);
    // This block cluster tree is only going to multiply the dense upper-left block, with no duplication of entries.
    tree->SetBlockTreeMode(BlockTreeMode::MatrixOnly);
    return tree;
}

void multiplyMetricWithVector(BlockClusterTree *tree, std::vector<double> &vec, std::vector<double> &output)
{
    // Copy input to an Eigen matrix
    Eigen::VectorXd in(vec.size());
    for (size_t i = 0; i < vec.size(); i++)
    {
        in(i) = vec[i];
    }
    // Set up an Eigen matrix as output
    Eigen::VectorXd out(vec.size());
    out.setZero();

    // Call the block cluster tree routine
    tree->Multiply(in, out);

    // Copy result to the output
    for (size_t i = 0; i < vec.size(); i++)
    {
        output[i] = out(i);
    }
}

double evaluateEnergy(PolyCurveNetwork *curve, BVHNode3D *root, double alpha, double beta)
{
    return SpatialTree::TPEnergyBH(curve, root, alpha, beta);
}

void evaluateGradient(PolyCurveNetwork *curve, BVHNode3D *root, std::vector<std::array<double, 3>> &out, double alpha, double beta)
{
    // Set up an Eigen matrix for the computation to use
    Eigen::MatrixXd grad(out.size(), 3);
    grad.setZero();

    // Use the BVH routine
    SpatialTree::TPEGradientBarnesHut(curve, root, grad, alpha, beta);

    // Copy result to the output
    for (size_t i = 0; i < out.size(); i++)
    {
        out[i][0] = grad(i, 0);
        out[i][1] = grad(i, 1);
        out[i][2] = grad(i, 2);
    }
}

} // namespace LWS
