#include "product/block_cluster_tree.h"
#include "utils.h"

#include <omp.h>
#include <thread>
#include <fstream>
#include <sstream>
#include <queue>
#include <map>

namespace LWS {

    long BlockClusterTree::illSepTime = 0;
    long BlockClusterTree::wellSepTime = 0;
    long BlockClusterTree::traversalTime = 0;

    BlockClusterTree::BlockClusterTree(PolyCurveNetwork* cg, BVHNode3D* tree, double sepCoeff, double a, double b, double e) {
        curves = cg;
        alpha = a;
        beta = b;
        separationCoeff = sepCoeff;
        epsilon = e;

        // std::cout << "Using " << nThreads << " threads." << std::endl;

        tree_root = tree;
        ClusterPair pair(tree, tree);
        unresolvedPairs.push_back(pair);

        nVerts = curves->NumVertices();
        constraintsSet = false;

        admissibleByCluster.resize(tree->numNodes);
        int depth = 0;
        while (unresolvedPairs.size() > 0) {
            splitInadmissibleNodes(depth);
            depth++;
        }

#ifdef DUMP_BCT_VISUALIZATION
        writeVisualization();
#endif

        // Premultiply A_f * 1, for reuse in later multiplications with G_f
        Af_1.setOnes(tree->numElements);
        Af_1 = MultiplyAf(Af_1);
        Af_1_low.setOnes(tree->numElements);
        Af_1_low = MultiplyAfLow(Af_1_low);

        mode = BlockTreeMode::MatrixOnly;
    }

    BlockClusterTree::~BlockClusterTree() {
        //delete threadpool;
    }

    void BlockClusterTree::splitInadmissibleNodes(int depth) {
        std::vector<ClusterPair> nextPairs;

        for (ClusterPair pair : unresolvedPairs) {
           pair.depth = depth;
            if (pair.cluster1->NumElements() == 0 || pair.cluster2->NumElements() == 0) {
                // Drop pairs where one of the sides has 0 vertices
                continue;
            }
            else if (pair.cluster1->NumElements() == 1 && pair.cluster2->NumElements() == 1) {
                // If this is two singleton vertices, put in the inadmissible list
                // so they get multiplied accurately
                inadmissiblePairs.push_back(pair);
            }
            else if (isPairAdmissible(pair, separationCoeff)) {
                // If the pair is admissible, mark it as such and leave it
                admissiblePairs.push_back(pair);
                admissibleByCluster[pair.cluster1->thisNodeID].push_back(pair);
            }
            else if (isPairSmallEnough(pair)) {
                inadmissiblePairs.push_back(pair);
            }
            else {
                // Otherwise, subdivide it into child pairs
                for (size_t i = 0; i < pair.cluster1->children.size(); i++) {
                    for (size_t j = 0; j < pair.cluster2->children.size(); j++) {
                        ClusterPair pair_ij(pair.cluster1->children[i], pair.cluster2->children[j]);
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
        double distance = norm(pair.cluster1->centerOfMass - pair.cluster2->centerOfMass);
        double radius1 = norm(pair.cluster1->maxBound().position - pair.cluster1->minBound().position) / 2;
        double radius2 = norm(pair.cluster2->maxBound().position - pair.cluster2->minBound().position) / 2;
        if (distance < radius1 || distance < radius2) return false;

        double ratio1 = pair.cluster1->nodeRatio(distance - radius2);
        double ratio2 = pair.cluster2->nodeRatio(distance - radius1);

        bool isAdm = fmax(ratio1, ratio2) < theta;
        return isAdm;
    }

    void BlockClusterTree::PrintData() {
        std::cout << admissiblePairs.size() << " admissible pairs" << std::endl;
        std::cout << inadmissiblePairs.size() << " inadmissible pairs" << std::endl;
    }

    void BlockClusterTree::PrintAdmissibleClusters(std::ofstream &stream) {
        for (ClusterPair p : admissiblePairs) {
            stream << p.cluster1->thisNodeID << ", " << p.cluster2->thisNodeID << std::endl;
        }
    }

    void BlockClusterTree::PrintInadmissibleClusters(std::ofstream &stream) {
        for (ClusterPair p : inadmissiblePairs) {
            stream << p.cluster1->thisNodeID << ", " << p.cluster2->thisNodeID << std::endl;
        }
    }

    void BlockClusterTree::SetBlockTreeMode(BlockTreeMode m) {
        mode = m;
    }

    void BlockClusterTree::sum_AIJ_VJ() const {
        // First accumulate the sums of a_IJ * V_J from admissible cluster pairs
        for (size_t i = 0; i < admissibleByCluster.size(); i++) {
            if (admissibleByCluster[i].size() > 0) {
                for (auto &pair : admissibleByCluster[i]) {
                    double a_IJ = SobolevCurves::MetricDistanceTerm(alpha, beta,
                        pair.cluster1->centerOfMass, pair.cluster2->centerOfMass,
                        pair.cluster1->averageTangent, pair.cluster2->averageTangent);
                    pair.cluster1->aIJ_VJ += a_IJ * pair.cluster2->V_I;
                }
            }
        }
    }

    void BlockClusterTree::sum_AIJ_VJ_Parallel() const {
        // First accumulate the sums of a_IJ * V_J from admissible cluster pairs
        #pragma omp parallel for shared(admissibleByCluster)
        for (size_t i = 0; i < admissibleByCluster.size(); i++) {
            if (admissibleByCluster[i].size() > 0) {
                for (auto &pair : admissibleByCluster[i]) {
                    double a_IJ = SobolevCurves::MetricDistanceTerm(alpha, beta,
                        pair.cluster1->centerOfMass, pair.cluster2->centerOfMass,
                        pair.cluster1->averageTangent, pair.cluster2->averageTangent);
                    pair.cluster1->aIJ_VJ += a_IJ * pair.cluster2->V_I;
                }
            }
        }
    }

    void BlockClusterTree::sum_AIJ_VJ_Low() const {
        // First accumulate the sums of a_IJ * V_J from admissible cluster pairs
        for (size_t i = 0; i < admissibleByCluster.size(); i++) {
            if (admissibleByCluster[i].size() > 0) {
                for (auto &pair : admissibleByCluster[i]) {
                    double a_IJ = SobolevCurves::MetricDistanceTermLow(alpha, beta,
                        pair.cluster1->centerOfMass, pair.cluster2->centerOfMass,
                        pair.cluster1->averageTangent, pair.cluster2->averageTangent);
                    pair.cluster1->aIJ_VJ += a_IJ * pair.cluster2->V_I;
                }
            }
        }
    }

    void BlockClusterTree::sum_AIJ_VJ_Low_Parallel() const {
        // First accumulate the sums of a_IJ * V_J from admissible cluster pairs
        #pragma omp parallel for shared(admissibleByCluster)
        for (size_t i = 0; i < admissibleByCluster.size(); i++) {
            if (admissibleByCluster[i].size() > 0) {
                for (auto &pair : admissibleByCluster[i]) {
                    double a_IJ = SobolevCurves::MetricDistanceTermLow(alpha, beta,
                        pair.cluster1->centerOfMass, pair.cluster2->centerOfMass,
                        pair.cluster1->averageTangent, pair.cluster2->averageTangent);
                    pair.cluster1->aIJ_VJ += a_IJ * pair.cluster2->V_I;
                }
            }
        }
    }

    void BlockClusterTree::MultiplyAdmissibleFast(const Eigen::MatrixXd &v_hat, Eigen::MatrixXd &b_hat) const {
        for (int i = 0; i < 3; i++) {
            b_hat.col(i) = MultiplyAf(v_hat.col(i));
        }

        b_hat = 2 * (Af_1.asDiagonal() * v_hat - b_hat);
    }

    void BlockClusterTree::MultiplyAdmissibleLowFast(const Eigen::VectorXd &v_mid, Eigen::VectorXd &b_mid) const {
        b_mid = MultiplyAfLow(v_mid);
        b_mid = 2 * (Af_1_low.asDiagonal() * v_mid - b_mid);
    }

    void BlockClusterTree::MultiplyInadmissibleLowParallel(const Eigen::VectorXd &v_mid, Eigen::VectorXd &b_mid) const {
        Eigen::VectorXd partialOutput = b_mid;
        partialOutput.setZero();
        
        #pragma omp parallel firstprivate(partialOutput) shared(v_mid, b_mid, inadmissiblePairs)
        {
            #pragma omp for
            for (size_t i = 0; i < inadmissiblePairs.size(); i++) {
                AfFullProductLow(inadmissiblePairs[i], v_mid, partialOutput);
            }
            #pragma omp critical
            {
                b_mid += partialOutput;
            }
        }
    }

    void BlockClusterTree::MultiplyInadmissibleParallel(const Eigen::MatrixXd &v_hat, Eigen::MatrixXd &b_hat) const {
        Eigen::MatrixXd partialOutput = b_hat;
        partialOutput.setZero();

        #pragma omp parallel firstprivate(partialOutput) shared(v_hat, b_hat, inadmissiblePairs)
        {
            #pragma omp for
            for (size_t i = 0; i < inadmissiblePairs.size(); i++) {
                AfFullProduct(inadmissiblePairs[i], v_hat, partialOutput);
            }
            #pragma omp critical
            {
                b_hat += partialOutput;
            }
        }
    }

    void BlockClusterTree::AfFullProduct(ClusterPair pair, const Eigen::MatrixXd &v_hat, Eigen::MatrixXd &result) const
    {
        std::vector<double> a_times_one(pair.cluster1->clusterIndices.size());
        std::vector<Vector3> a_times_v(pair.cluster1->clusterIndices.size());
        
        for (size_t i = 0; i < pair.cluster1->clusterIndices.size(); i++) {
            int e1index = pair.cluster1->clusterIndices[i];
            CurveEdge* p1 = curves->GetEdge(e1index);
            double l1 = tree_root->bvhRoot->fullMasses(e1index);

            for (size_t j = 0; j < pair.cluster2->clusterIndices.size(); j++) {
                int e2index = pair.cluster2->clusterIndices[j];
                CurveEdge* p2 = curves->GetEdge(e2index);
                bool isNeighbors = (p1 == p2 || p1->IsNeighbors(p2));
                
                Vector3 mid1 = p1->Midpoint();
                Vector3 mid2 = p2->Midpoint();
                Vector3 tan1 = p1->Tangent();
                Vector3 tan2 = p2->Tangent();

                double l2 = tree_root->bvhRoot->fullMasses(e2index);

                // Save on a few operations by only multiplying l2 now,
                // and multiplying l1 only once, after inner loop
                double distTerm = SobolevCurves::MetricDistanceTerm(alpha, beta, mid1, mid2, tan1, tan2);
                double af_ij = (isNeighbors) ? 0 : l2  * distTerm;

                // We dot this row of Af(i, j) with the all-ones vector, which means we
                // just add up all entries of that row.
                a_times_one[i] += af_ij;

                // We also dot it with v_hat(J).
                a_times_v[i] += af_ij * SelectRow(v_hat, e2index);
            }

            a_times_one[i] *= l1;
            a_times_v[i] *= l1;

            // We've computed everything from row i now, so add to the results vector
            Vector3 toAdd = 2 * (a_times_one[i] * SelectRow(v_hat, e1index) - a_times_v[i]);
            AddToRow(result, e1index, toAdd);
        }
    }

    void BlockClusterTree::AfFullProductLow(ClusterPair pair, const Eigen::VectorXd &v_mid, Eigen::VectorXd &result) const
    {
        std::vector<double> a_times_one(pair.cluster1->clusterIndices.size());
        std::vector<double> a_times_v(pair.cluster1->clusterIndices.size());

        for (size_t i = 0; i < pair.cluster1->clusterIndices.size(); i++) {
            int e1index = pair.cluster1->clusterIndices[i];
            CurveEdge* p1 = curves->GetEdge(e1index);
            double l1 = tree_root->bvhRoot->fullMasses(e1index);

            for (size_t j = 0; j < pair.cluster2->clusterIndices.size(); j++) {
                int e2index = pair.cluster2->clusterIndices[j];
                CurveEdge* p2 = curves->GetEdge(e2index);
                bool isNeighbors = (p1 == p2 || p1->IsNeighbors(p2));
                
                Vector3 mid1 = p1->Midpoint();
                Vector3 mid2 = p2->Midpoint();
                Vector3 tan1 = p1->Tangent();
                Vector3 tan2 = p2->Tangent();

                double l2 = tree_root->bvhRoot->fullMasses(e2index);

                // Save on a few operations by only multiplying l2 now,
                // and multiplying l1 only once, after inner loop
                double distTerm = SobolevCurves::MetricDistanceTermLow(alpha, beta, mid1, mid2, tan1, tan2);

                double af_ij = (isNeighbors) ? 0 : l2  * distTerm;

                // We dot this row of Af(i, j) with the all-ones vector, which means we
                // just add up all entries of that row.
                a_times_one[i] += af_ij;

                // We also dot it with v_hat(J).
                a_times_v[i] += af_ij * v_mid(e2index);
            }

            a_times_one[i] *= l1;
            a_times_v[i] *= l1;

            // We've computed everything from row i now, so add to the results vector
            double toAdd = 2 * (a_times_one[i] * v_mid(e1index) - a_times_v[i]);
            result(e1index) += toAdd;
        }
    }
    
    void BlockClusterTree::AfApproxProduct(ClusterPair pair, const Eigen::MatrixXd &v_hat, Eigen::MatrixXd &result) const
    {
        long traversalStart = Utils::currentTimeMilliseconds();
        Eigen::VectorXd wf_i;
        pair.cluster1->fillClusterMassVector(wf_i);
        Eigen::VectorXd wf_j;
        pair.cluster2->fillClusterMassVector(wf_j);
        long traversalEnd = Utils::currentTimeMilliseconds();

        traversalTime += (traversalEnd - traversalStart);
        
        double a_IJ = SobolevCurves::MetricDistanceTerm(alpha, beta,
            pair.cluster1->centerOfMass, pair.cluster2->centerOfMass,
            pair.cluster1->averageTangent, pair.cluster2->averageTangent);
        // Evaluate a(I,J) * w_f(J)^T * 1(J)
        double a_wf_1 = a_IJ * wf_j.sum();

        // Evaluate a(I,J) * w_f(J)^T * v_hat(J)
        Vector3 a_wf_J{0, 0, 0};
        // Dot w_f(J) with v_hat(J)
        for (int j = 0; j < wf_j.rows(); j++) {
            a_wf_J += wf_j(j) * SelectRow(v_hat, pair.cluster2->clusterIndices[j]);
        }
        a_wf_J *= a_IJ;

        // Add in the results
        for (int i = 0; i < wf_i.rows(); i++) {
            Vector3 toAdd = wf_i[i] * 2 * (a_wf_1 * SelectRow(v_hat, pair.cluster1->clusterIndices[i]) - a_wf_J);
            CurveEdge* e_i = curves->GetEdge(pair.cluster1->clusterIndices[i]);
            AddToRow(result, pair.cluster1->clusterIndices[i], toAdd);
        }
    }
    
    void BlockClusterTree::AfApproxProductLow(ClusterPair pair, const Eigen::VectorXd &v_mid, Eigen::VectorXd &result) const
    {
        long traversalStart = Utils::currentTimeMilliseconds();
        Eigen::VectorXd wf_i;
        pair.cluster1->fillClusterMassVector(wf_i);
        Eigen::VectorXd wf_j;
        pair.cluster2->fillClusterMassVector(wf_j);
        long traversalEnd = Utils::currentTimeMilliseconds();

        traversalTime += (traversalEnd - traversalStart);
        
        double a_IJ = SobolevCurves::MetricDistanceTermLow(alpha, beta,
            pair.cluster1->centerOfMass, pair.cluster2->centerOfMass,
            pair.cluster1->averageTangent, pair.cluster2->averageTangent);

        // Evaluate a(I,J) * w_f(J)^T * 1(J)
        double a_wf_1 = a_IJ * wf_j.sum();

        // Evaluate a(I,J) * w_f(J)^T * v_hat(J)
        double a_wf_J = 0;
        // Dot w_f(J) with v_hat(J)
        for (int j = 0; j < wf_j.rows(); j++) {
            a_wf_J += wf_j(j) * v_mid(pair.cluster2->clusterIndices[j]);
        }
        a_wf_J *= a_IJ;

        // Add in the results
        for (int i = 0; i < wf_i.rows(); i++) {
            double toAdd = wf_i[i] * 2 * (a_wf_1 * v_mid(pair.cluster1->clusterIndices[i]) - a_wf_J);
            CurveEdge* e_i = curves->GetEdge(pair.cluster1->clusterIndices[i]);
            result(pair.cluster1->clusterIndices[i]) += toAdd;
        }
    }
}
