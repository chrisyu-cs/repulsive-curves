#include "product/block_cluster_tree.h"
#include "utils.h"

#include <thread>
#include <fstream>
#include <sstream>
#include <queue>
#include <map>

namespace LWS {

    long BlockClusterTree::illSepTime = 0;
    long BlockClusterTree::wellSepTime = 0;
    long BlockClusterTree::traversalTime = 0;

    // Set up thread pool here, since we need it for premultiplying Af
    int nMachineThreads = std::thread::hardware_concurrency();
    int BlockClusterTree::nThreads = std::max( 2, nMachineThreads );
    progschj::ThreadPool* BlockClusterTree::threadpool = new progschj::ThreadPool(nThreads);

    BlockClusterTree::BlockClusterTree(PolyCurveNetwork* cg, BVHNode3D* tree, double sepCoeff, double a, double b, double e) {
        curves = cg;
        alpha = a;
        beta = b;
        separationCoeff = sepCoeff;
        epsilon = e;

        std::cout << "Using " << nThreads << " threads." << std::endl;

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
        Af_1.setOnes(nVerts);
        Af_1 = MultiplyAf(Af_1);
        Af_1_low.setOnes(nVerts);
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
        for (size_t i = 0; i < admissibleByCluster.size(); i++) {
            if (admissibleByCluster[i].size() > 0) {
                auto multiply = [&, i]() {
                    for (auto &pair : admissibleByCluster[i]) {
                        double a_IJ = SobolevCurves::MetricDistanceTerm(alpha, beta,
                            pair.cluster1->centerOfMass, pair.cluster2->centerOfMass,
                            pair.cluster1->averageTangent, pair.cluster2->averageTangent);
                        pair.cluster1->aIJ_VJ += a_IJ * pair.cluster2->V_I;
                    }
                };
                threadpool->enqueue(multiply);
                // multiply();
            }
        }
        
        threadpool->wait_until_empty();
        threadpool->wait_until_nothing_in_flight();
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
        for (size_t i = 0; i < admissibleByCluster.size(); i++) {
            if (admissibleByCluster[i].size() > 0) {
                auto multiply = [&, i]() {
                    for (auto &pair : admissibleByCluster[i]) {
                        double a_IJ = SobolevCurves::MetricDistanceTermLow(alpha, beta,
                            pair.cluster1->centerOfMass, pair.cluster2->centerOfMass,
                            pair.cluster1->averageTangent, pair.cluster2->averageTangent);
                        pair.cluster1->aIJ_VJ += a_IJ * pair.cluster2->V_I;
                    }
                };
                threadpool->enqueue(multiply);
                // multiply();
            }
        }

        threadpool->wait_until_empty();
        threadpool->wait_until_nothing_in_flight();
    }

    void BlockClusterTree::MultiplyAdmissible(Eigen::MatrixXd &v_hat, Eigen::MatrixXd &b_hat) const {
        for (ClusterPair pair : admissiblePairs) {
            AfApproxProduct(pair, v_hat, b_hat);
        }
    }

    void BlockClusterTree::MultiplyAdmissibleExact(Eigen::MatrixXd &v_hat, Eigen::MatrixXd &b_hat) const {
        for (ClusterPair pair : admissiblePairs) {
            AfFullProduct(pair, v_hat, b_hat);
        }
    }

    void BlockClusterTree::MultiplyAdmissibleFast(const Eigen::MatrixXd &v_hat, Eigen::MatrixXd &b_hat) const {
        for (int i = 0; i < 3; i++) {
            b_hat.col(i) = MultiplyAf(v_hat.col(i));
        }

        b_hat = 2 * (Af_1.asDiagonal() * v_hat - b_hat);
    }

    void BlockClusterTree::MultiplyAdmissibleLowExact(const Eigen::VectorXd &v_mid, Eigen::VectorXd &b_mid) const {
        for (auto &pair : admissiblePairs) {
            AfFullProductLow(pair, v_mid, b_mid);
        }
    }

    void BlockClusterTree::MultiplyAdmissibleLow(const Eigen::VectorXd &v_mid, Eigen::VectorXd &b_mid) const {
        for (auto &pair : admissiblePairs) {
            AfApproxProductLow(pair, v_mid, b_mid);
        }
    }

    void BlockClusterTree::MultiplyAdmissibleLowFast(const Eigen::VectorXd &v_mid, Eigen::VectorXd &b_mid) const {
        b_mid = MultiplyAfLow(v_mid);
        b_mid = 2 * (Af_1_low.asDiagonal() * v_mid - b_mid);

        // for (auto &pair : admissiblePairs) {
        //     AfApproxProductLow(pair, v_mid, b_mid);
        // }

        // for (auto &pair : admissiblePairs) {
        //     AfFullProductLow(pair, v_mid, b_mid);
        // }
    }

    void BlockClusterTree::MultiplyInadmissible(const Eigen::MatrixXd &v_hat, Eigen::MatrixXd &b_hat, int startIndex, int endIndex) const {
        for (int i = startIndex; i < endIndex; i++) {
            AfFullProduct(inadmissiblePairs[i], v_hat, b_hat);
        }
    }

    void BlockClusterTree::MultiplyInadmissibleLow(const Eigen::VectorXd &v_mid, Eigen::VectorXd &b_mid, int startIndex, int endIndex) const {
        for (int i = startIndex; i < endIndex; i++) {
            AfFullProductLow(inadmissiblePairs[i], v_mid, b_mid);
        }
    }

    void BlockClusterTree::MultiplyInadmissibleLowParallel(const Eigen::VectorXd &v_mid, Eigen::VectorXd &b_mid) const {
        int nPairs = inadmissiblePairs.size();
        // How many would be left if we rounded down to the nearest multiple of nThreads
        int leftovers = nPairs % nThreads;
        std::vector<int> counts(nThreads);
        // Split into nThreads roughly equal groups
        for (int i = 0; i < nThreads; i++) {
            counts[i] = nPairs / nThreads;
            // Add the extras from rounding down, so that the totals match up
            if (i < leftovers) counts[i]++;
        }

        std::vector<Eigen::VectorXd> out_bs(nThreads);

        for (int i = 0; i < nThreads; i++) {
            out_bs[i].setZero(b_mid.rows());
        }

        // std::vector<std::thread> threads;

        int startIndex = 0;
        for (int i = 0; i < nThreads; i++) {
            int endIndex = startIndex + counts[i];
            
            auto multiply = [&, i, startIndex, endIndex]() {
                MultiplyInadmissibleLow(v_mid, out_bs[i], startIndex, endIndex);
            };
            // threads.push_back(std::thread(multiply));
            threadpool->enqueue(multiply);

            startIndex = endIndex;
        }

        // for (int i = 0; i < nThreads; i++) {
        //     threads[i].join();
        // }
        threadpool->wait_until_empty();
        threadpool->wait_until_nothing_in_flight();

        // Add up all the result vectors
        b_mid = out_bs[0];
        for (int i = 1; i < nThreads; i++) {
            b_mid += out_bs[i];
        }
    }

    void BlockClusterTree::MultiplyInadmissibleParallel(const Eigen::MatrixXd &v_hat, Eigen::MatrixXd &b_hat) const {
        int nPairs = inadmissiblePairs.size();
        // How many would be left if we rounded down to the nearest multiple of nThreads
        int leftovers = nPairs % nThreads;
        std::vector<int> counts(nThreads);
        // Split into nThreads roughly equal groups
        for (int i = 0; i < nThreads; i++) {
            counts[i] = nPairs / nThreads;
            // Add the extras from rounding down, so that the totals match up
            if (i < leftovers) counts[i]++;
        }

        std::vector<Eigen::MatrixXd> out_bs(nThreads);

        for (int i = 0; i < nThreads; i++) {
            out_bs[i].setZero(b_hat.rows(), b_hat.cols());
        }

        // std::vector<std::thread> threads;

        int startIndex = 0;
        for (int i = 0; i < nThreads; i++) {
            int endIndex = startIndex + counts[i];
            
            auto multiply = [&, i, startIndex, endIndex]() {
                MultiplyInadmissible(v_hat, out_bs[i], startIndex, endIndex);
            };
            // threads.push_back(std::thread(multiply));
            threadpool->enqueue(multiply);

            startIndex = endIndex;
        }

        // for (int i = 0; i < nThreads; i++) {
        //     threads[i].join();
        // }
        threadpool->wait_until_empty();
        threadpool->wait_until_nothing_in_flight();

        // Add up all the result vectors
        b_hat = out_bs[0];
        for (int i = 1; i < nThreads; i++) {
            b_hat += out_bs[i];
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

    void BlockClusterTree::AfFullProductLowVerts(ClusterPair pair, const Eigen::VectorXd &v_in, Eigen::VectorXd &result) const
    {
        for (size_t i = 0; i < pair.cluster1->clusterIndices.size(); i++) {
            int e1index = pair.cluster1->clusterIndices[i];
            CurveEdge* p1 = curves->GetEdge(e1index);
            double len1 = p1->Length();
            double l1 = tree_root->bvhRoot->fullMasses(e1index);

            for (size_t j = 0; j < pair.cluster2->clusterIndices.size(); j++) {
                int e2index = pair.cluster2->clusterIndices[j];
                CurveEdge* p2 = curves->GetEdge(e2index);
                double len2 = p2->Length();
                bool isNeighbors = (p1 == p2 || p1->IsNeighbors(p2));
                if (isNeighbors) continue;

                CurveVertex* endpoints[4] = {p1->prevVert, p1->nextVert, p2->prevVert, p2->nextVert};

                double kf_st = SobolevCurves::MetricDistanceTermLow(alpha, beta,
                    p1->Midpoint(), p2->Midpoint(), p1->Tangent(), p2->Tangent());

                for (CurveVertex* u : endpoints) {
                    for (CurveVertex* v : endpoints) {
                        double u_s = SobolevCurves::HatMidpoint(p1, u);
                        double u_t = SobolevCurves::HatMidpoint(p2, u);
                        double v_s = SobolevCurves::HatMidpoint(p1, v);
                        double v_t = SobolevCurves::HatMidpoint(p2, v);

                        double numer = (u_s - u_t) * (v_s - v_t);
                        int index_u = u->GlobalIndex();
                        int index_v = v->GlobalIndex();

                        result(index_v) += v_in(index_u) * numer * kf_st * len1 * len2;
                    }
                }
            }
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

    void BlockClusterTree::CompareBlocks() {
        Eigen::MatrixXd fullBlock;
        Eigen::MatrixXd approxBlock;

        double totalError = 0;
        double totalNorm = 0;

        for (ClusterPair pair : inadmissiblePairs) {
            // Inadmissible blocks are computed exactly, so no error contribution
            fullBlock = AfFullBlock(pair);

            double normFull = fullBlock.norm();
            totalNorm += normFull * normFull;
        }
        for (ClusterPair pair : admissiblePairs) {
            fullBlock = AfFullBlock(pair);
            approxBlock = AfApproxBlock(pair);

            double normFull = fullBlock.norm();
            double normDiff = (fullBlock - approxBlock).norm();
            double relative = 100 * normDiff / normFull;

            if (relative > 50) {
                std::cout << "(" << pair.cluster1->NumElements() << ", " << pair.cluster2->NumElements() << ")" << std::endl;
                std::cout << "Full:\n" << fullBlock << std::endl;
                std::cout << "Approx:\n" << approxBlock << std::endl;
                std::cout << "Error: " << normDiff << " (" << relative << " percent)" << std::endl;
            }

            totalNorm += normFull * normFull;
            totalError += normDiff * normDiff;
        }

        totalError = sqrt(totalError);
        totalNorm = sqrt(totalNorm);
        double totalRelative = 100 * totalError / totalNorm;

        std::cout << "Total error = " << totalError << " (" << totalRelative << " percent; total norm = " << totalNorm << ")" << std::endl;
    }

    Eigen::MatrixXd BlockClusterTree::AfFullBlock(ClusterPair pair) {
        std::vector<VertexBody6D> children1;
        pair.cluster1->accumulateChildren(children1);
        std::vector<VertexBody6D> children2;
        pair.cluster2->accumulateChildren(children2);
        
        Eigen::MatrixXd block;
        block.setZero(children1.size(), children2.size());

        for (size_t i = 0; i < children1.size(); i++) {
            for (size_t j = 0; j < children2.size(); j++) {
                CurveEdge* p1 = curves->GetEdge(children1[i].elementIndex);
                CurveEdge* p2 = curves->GetEdge(children2[j].elementIndex);

                bool isNeighbors = p1 == p2 || p1->IsNeighbors(p2);

                double w_i = children1[i].mass;
                double w_j = children2[j].mass;

                Vector3 c_i = children1[i].pt.position;
                Vector3 c_j = children2[j].pt.position;
                Vector3 t_i = children1[i].pt.tangent;
                Vector3 t_j = children2[j].pt.tangent;

                block(i, j) = (isNeighbors) ? 0 : -w_i * w_j * SobolevCurves::MetricDistanceTerm(alpha, beta, c_i, c_j, t_i, t_j);
            }
        }

        return block;
    }

    Eigen::MatrixXd BlockClusterTree::AfApproxBlock(ClusterPair pair) {
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
        double a_IJ = SobolevCurves::MetricDistanceTerm(alpha, beta, 
            pair.cluster1->centerOfMass, pair.cluster2->centerOfMass,
            pair.cluster1->averageTangent, pair.cluster2->averageTangent);
        
        Eigen::MatrixXd block;
        block.setZero(children1.size(), children2.size());

        for (size_t i = 0; i < children1.size(); i++) {
            for (size_t j = 0; j < children2.size(); j++) {
                block(i, j) = -wf_i[i] * a_IJ * wf_j[j];
            }
        }

        return block;
    }

#ifdef DUMP_BCT_VISUALIZATION
    void BlockClusterTree::writeVisualization()
    {
       int nE = curves->NumEdges();

       // create filenames that depend on the number of vertices
       // (so that we get different files for each level of the
       // multigrid hierarchy)
       std::string nLeavesString = std::to_string( nE );
       std::stringstream fnK, fnPK, fnPK0, fnBCT;
       fnK   << "K_"   << nLeavesString << ".csv";
       fnPK  << "PK_"  << nLeavesString << ".csv";
       fnPK0 << "PK0_" << nLeavesString << ".csv";
       fnBCT << "BCT_" << nLeavesString << ".svg";

       // write a CSV file that stores a dense kernel matrix, with
       // row/column order determined by the blocking used in the BCT
       std::ofstream outK( fnK.str().c_str() );
       std::ofstream outPK( fnPK.str().c_str() );
       std::ofstream outPK0( fnPK0.str().c_str() );
       std::ofstream outBCT( fnBCT.str().c_str() );

       // compute start and end indices for each internal node
       // of the BVH that span the indices of their children
       tree_root->indexNodes( 0 );

       // allocate storage for the full, reordered, and hierarchical matrices
       std::vector< std::vector<double> > K( nE ), PK( nE ), PK0( nE );
       for( int i = 0; i < nE; i++ )
       {
            K[i].resize( nE );
           PK[i].resize( nE );
          PK0[i].resize( nE );
       }
       
       // construct a map from curve vertex
       // indices to BVH leaf indices
       std::map<int,int> element2leafIndex;
       std::queue<BVHNode3D*> Q;
       Q.push( tree_root );
       while( !Q.empty() )
       {
          BVHNode3D* node = Q.front(); Q.pop();
          for( BVHNode3D* child : node->children )
          {
             Q.push( child );

             if( child->isLeaf )
             {
                int elementIndex = child->body.elementIndex;
                int leafIndex = child->indexStart;
                element2leafIndex[elementIndex] = leafIndex;
             }
          }
       }

       // iterate over all pairs and build the dense matrices
       for( int i = 0; i < nE; i++ )
       {
          CurveEdge* ei = curves->GetEdge( i );
          int Pi = element2leafIndex[i];
          Vector3 xi = ei->Midpoint();
          Vector3 Ti = ei->Tangent();

          for( int j = 0; j < nE; j++ )
          {
             CurveEdge* ej = curves->GetEdge( j );
             int Pj = element2leafIndex[j];
             Vector3 xj = ej->Midpoint();
             Vector3 Tj = ej->Tangent();

             //double kij = SobolevCurves::MetricDistanceTerm( alpha, beta, xi, xj, Ti, Tj );
             double kij = TPESC::tpe_pair_pts( xi, xj, Ti, 1., 1., alpha, beta );
             if( isnan(kij) || isinf(kij)) kij = 0.;

             K[i][j] = kij;
             PK[Pi][Pj] = kij;
             PK0[Pi][Pj] = kij;
          }
       }

       // evaluate the energy approximation in each admissible block,
       // and write it to all the entries in that block
       outBCT << "<svg>" << std::endl;
       outBCT << "<style type=\"text/css\">" << std::endl;
       outBCT << ".st0{fill:none;stroke:#000000;stroke-width:32.0;}" << std::endl;
       outBCT << ".st1{fill:none;stroke:#000000;stroke-width:16.0;}" << std::endl;
       outBCT << ".st2{fill:none;stroke:#000000;stroke-width:8.0;}" << std::endl;
       outBCT << ".st3{fill:none;stroke:#000000;stroke-width:4.0;}" << std::endl;
       outBCT << ".st4{fill:none;stroke:#000000;stroke-width:2.0;}" << std::endl;
       outBCT << ".st5{fill:none;stroke:#000000;stroke-width:1.0;}" << std::endl;
       outBCT << ".st6{fill:none;stroke:#000000;stroke-width:0.5;}" << std::endl;
       outBCT << ".st7{fill:none;stroke:#000000;stroke-width:0.25;}" << std::endl;
       outBCT << ".st8{fill:none;stroke:#000000;stroke-width:0.125;}" << std::endl;
       outBCT << ".st9{fill:none;stroke:#000000;stroke-width:0.0625;}" << std::endl;
       outBCT << "</style>" << std::endl;

       for( const ClusterPair& AB : admissiblePairs )
       {
          BVHNode3D* A = AB.cluster1;
          BVHNode3D* B = AB.cluster2;

          int A0 = A->indexStart;
          int A1 = A->indexEnd;
          int B0 = B->indexStart;
          int B1 = B->indexEnd;

          Vector3 xA = A->centerOfMass;
          Vector3 xB = B->centerOfMass;
          Vector3 TA = A->averageTangent;
          Vector3 TB = B->averageTangent;

          //double kAB = SobolevCurves::MetricDistanceTerm( alpha, beta, xA, xB, TA, TB );
          double kAB = TPESC::tpe_pair_pts( xA, xB, TA, 1., 1., alpha, beta );
          if( isnan(kAB) || isinf(kAB)) kAB = 0.;

          for( int Pi = A0; Pi < A1; Pi++ )
          for( int Pj = B0; Pj < B1; Pj++ )
          {
             PK0[Pi][Pj] = kAB;
          }

          // use depth in tree to determine line width
          int d = std::min( AB.depth, 9 );

          outBCT << "<rect ";
          outBCT << "class=\"st" << d << "\" ";
          outBCT << "x=\"" << A0 << "\" ";
          outBCT << "y=\"" << B0 << "\" ";
          outBCT << "width=\"" << A1 - A0 << "\" ";
          outBCT << "height=\"" << B1 - B0 << "\"";
          outBCT << "/>" << std::endl;
       }
       outBCT << "</svg>" << std::endl;

       // write matrices to disk
       for( int i = 0; i < nE; i++ )
       {
          for( int j = 0; j < nE; j++ )
          {
             outK   <<   K[i][j];
             outPK  <<  PK[i][j];
             outPK0 << PK0[i][j];

             if( j != nE-1 )
             {
                outK   << ",";
                outPK  << ",";
                outPK0 << ",";
             }
          }
          outK   << std::endl;
          outPK  << std::endl;
          outPK0 << std::endl;
       }
    }
#endif
}
