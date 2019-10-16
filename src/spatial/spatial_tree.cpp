#include "spatial/spatial_tree.h"

namespace LWS {
    SpatialTree::~SpatialTree() {}


    void SpatialTree::TPEGradientBarnesHut(PolyCurveNetwork* curveNetwork, SpatialTree *root, Eigen::MatrixXd &output, double alpha, double beta) {
        // The single energy term (i, j) affects six vertices:
        // (i_prev, i, i_next, j_prev, j, j_next).
        // We can restructure the computation as follows:
        // for each single 1-ring (i, i_prev, i_next), accumulate the
        // contributions from the gradients of both terms (i, j) and (j, i).
        int nVerts = curveNetwork->NumVertices();
        output.setZero();

        for (int i = 0; i < nVerts; i++) {
            CurveVertex* i_pt = curveNetwork->GetVertex(i);
            root->accumulateTPEGradient(output, i_pt, curveNetwork, alpha, beta);
        }
    }
}
