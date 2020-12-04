#include "spatial/spatial_tree.h"
#include <omp.h>

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
        Eigen::MatrixXd partialOutput = output;

        #pragma omp parallel firstprivate(partialOutput) shared(root, output)
        {
            #pragma omp for
            for (int i = 0; i < nVerts; i++)
            {
                CurveVertex* i_pt = curveNetwork->GetVertex(i);
                root->accumulateTPEGradient(partialOutput, i_pt, curveNetwork, alpha, beta);
            }

            #pragma omp critical
            {
                output += partialOutput;
            }
        }
    }

    double SpatialTree::TPEnergyBH(PolyCurveNetwork* curveNetwork, SpatialTree *root, double alpha, double beta) {
        int nVerts = curveNetwork->NumVertices();
        double fullSum = 0;
        
        #pragma omp parallel for reduction(+ : fullSum) shared(root)
        // Loop over all vertices and add up energy contributions
        for (int i = 0; i < nVerts; i++) {
            CurveVertex* i_pt = curveNetwork->GetVertex(i);
            double vertSum = 0;
            root->accumulateVertexEnergy(vertSum, i_pt, curveNetwork, alpha, beta);
            fullSum += vertSum;
        }
        return fullSum;
    }
}
