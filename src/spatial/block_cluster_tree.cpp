#include "spatial/block_cluster_tree.h"
#include "utils.h"
#include "sobo_slobo.h"

namespace LWS {

    BlockClusterTree::BlockClusterTree(PolyCurveGroup* cg, BVHNode3D* tree, double sepCoeff, double a, double b) {
        curves = cg;
        alpha = a;
        beta = b;
        separationCoeff = sepCoeff;

        tree_root = tree;
        ClusterPair pair{tree, tree};
        unresolvedPairs.push_back(pair);

        while (unresolvedPairs.size() > 0) {
            splitInadmissibleNodes();
        }
    }

    void BlockClusterTree::splitInadmissibleNodes() {
        std::vector<ClusterPair> nextPairs;

        for (ClusterPair pair : unresolvedPairs) {
            if (isPairAdmissible(pair, separationCoeff)) {
                // If the pair is admissible, mark it as such and leave it
                admissiblePairs.push_back(pair);
            }
            else if (isPairSmallEnough(pair)) {
                imadmissiblePairs.push_back(pair);
            }
            else {
                // Otherwise, subdivide it into child pairs
                for (size_t i = 0; i < pair.cluster1->children.size(); i++) {
                    for (size_t j = 0; j < pair.cluster2->children.size(); j++) {
                        ClusterPair pair_ij{pair.cluster1->children[i], pair.cluster2->children[j]};
                        nextPairs.push_back(pair_ij);
                    }
                }
            }
        }
        // Replace the inadmissible pairs by the next set
        unresolvedPairs.clear();
        unresolvedPairs = nextPairs;
    }

    bool BlockClusterTree::isPairSmallEnough(ClusterPair pair) {
        int s1 = pair.cluster1->NumElements();
        int s2 = pair.cluster2->NumElements();
        return (s1 <= 1) || (s2 <= 1) || (s1 + s2 <= 8);
    }

    bool BlockClusterTree::isPairAdmissible(ClusterPair pair, double theta) {
        if (pair.cluster1 == pair.cluster2) return false;

        Vector3 max1 = pair.cluster1->maxBound().position;
        Vector3 max2 = pair.cluster2->maxBound().position;
        Vector3 min1 = pair.cluster1->minBound().position;
        Vector3 min2 = pair.cluster2->minBound().position;

        // Compute diagonals for the cells
        Vector3 diag1 = max1 - min1;
        Vector3 diag2 = max2 - min2;

        // Compute centers of the bounding boxes
        Vector3 center1 = (max1 + min1) / 2;
        Vector3 center2 = (max2 + min2) / 2;

        // Compute lines between centers of boxes
        Vector3 centerLine = center2 - center1;
        double centerDistance = norm(centerLine);
        centerLine /= centerDistance;

        // Compute the distance between the edges of boxes
        double diam1 = dot(diag1, centerLine);
        double diam2 = dot(diag2, centerLine);
        double boxDistance = centerDistance - (diam1 / 2) - (diam2 / 2);

        return fmax(diam1, diam2) <= theta * boxDistance;
    }

    void BlockClusterTree::PrintData() {
        std::cout << admissiblePairs.size() << " admissible pairs" << std::endl;
        std::cout << imadmissiblePairs.size() << " inadmissible pairs" << std::endl;
    }

    void BlockClusterTree::MultiplyVector(std::vector<double> &v, std::vector<double> &b) {
        std::vector<Vector3> v_hat(v.size());
        SobolevCurves::ApplyDf(curves, v, v_hat);

        std::vector<Vector3> b_hat(v.size());

        ClusterPair top{tree_root, tree_root};

        for (ClusterPair pair : imadmissiblePairs) {
            AfFullProduct_hat(pair, v_hat, b_hat);
        }
        for (ClusterPair pair : admissiblePairs) {
            AfFullProduct_hat(pair, v_hat, b_hat);
        }

        SobolevCurves::ApplyDfTranspose(curves, b_hat, b);
    }

    void BlockClusterTree::AfFullProduct_hat(ClusterPair pair, std::vector<Vector3> &v_hat, std::vector<Vector3> &result)
    {
        std::vector<VertexBody6D> children1;
        pair.cluster1->accumulateChildren(children1);
        std::vector<VertexBody6D> children2;
        pair.cluster2->accumulateChildren(children2);
        
        double pow_s = (beta - alpha);

        std::vector<double> a_times_one(children1.size());
        std::vector<Vector3> a_times_v(children1.size());
        
        for (size_t i = 0; i < children1.size(); i++) {
            VertexBody6D e1 = children1[i];

            for (size_t j = 0; j < children2.size(); j++) {
                VertexBody6D e2 = children2[j];
                
                double mass = 1;
                double af_ij = (e1.vertIndex1 == e2.vertIndex1) ? 0 :
                    mass / pow(norm(e1.pt.position - e2.pt.position), pow_s);

                // We dot this row of Af(i, j) with the all-ones vector, which means we
                // just add up all entries of that row.
                a_times_one[i] += af_ij;

                // We also dot it with v_hat(J).
                a_times_v[i] += af_ij * v_hat[e2.vertIndex1];
            }

            // We've computed everything from row i now, so add to the results vector
            result[e1.vertIndex1] += a_times_one[i] * v_hat[e1.vertIndex1] - a_times_v[i];
        }
    }
    
    void BlockClusterTree::AfApproxProduct_hat(ClusterPair pair, std::vector<Vector3> &v_hat, std::vector<Vector3> &result)
    {
        std::vector<VertexBody6D> children1;
        pair.cluster1->accumulateChildren(children1);
        std::vector<VertexBody6D> children2;
        pair.cluster2->accumulateChildren(children2);
        
        double pow_s = beta - alpha;

        std::vector<double> wf_i(children1.size());
        std::vector<double> wf_j(children2.size());

        for (size_t i = 0; i < children1.size(); i++) {
            wf_i[i] = children1[i].mass;
        }

        for (size_t j = 0; j < children2.size(); j++) {
            wf_j[j] = children2[j].mass;
        }

        double a_IJ = 1.0 / pow(norm(pair.cluster1->centerOfMass - pair.cluster2->centerOfMass), pow_s);
        // Evaluate a(I,J) * w_f(J)^T * 1(J)
        double a_wf_1 = a_IJ * std_vector_sum_entries(wf_j);

        // Evaluate a(I,J) * w_f(J)^T * v_hat(J)
        Vector3 a_wf_J{0, 0, 0};
        // Dot w_f(J) with v_hat(J)
        for (size_t j = 0; j < children2.size(); j++) {
            a_wf_J += wf_j[j] * v_hat[children2[j].vertIndex1];
        }
        a_wf_J *= a_IJ;

        // Add in the results
        for (size_t i = 0; i < children1.size(); i++) {
            // TODO: why is it off by a factor of 2?
            result[children1[i].vertIndex1] += (wf_i[i] * a_wf_1 * v_hat[children1[i].vertIndex1] - wf_i[i] * a_wf_J);
        }
    }
}