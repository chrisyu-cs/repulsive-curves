#include "lws_cluster.h"

#include <Eigen/Eigenvalues>
#include "vertexderivatives.h" 

using namespace std;
using namespace geometrycentral::surface;

namespace LWS {

    bool CompareEntries::operator()(const ClusterQueueEntry &x, const ClusterQueueEntry &y) {
        return x.priority > y.priority;
    }

    // ====================== Aggregate clustering functions ======================

    LWSClustering::LWSClustering() {
        queue = new std::priority_queue<ClusterQueueEntry, std::vector<ClusterQueueEntry>, CompareEntries>();
    }

    LWSClustering::LWSClustering(polyscope::SurfaceMesh* s, HalfedgeMesh* m, VertexPositionGeometry* g) {
        mesh = m;
        geom = g;
        surface = s;
        precomputeScalars(mesh, geom);
        vertexOwnership = VertexData<int>(*mesh);
        queue = new std::priority_queue<ClusterQueueEntry, std::vector<ClusterQueueEntry>, CompareEntries>();
    }

    void LWSClustering::precomputeScalars(HalfedgeMesh* mesh, VertexPositionGeometry* geom) {
        areas = VertexData<double>(*mesh);
        meanCurvatures = VertexData<double>(*mesh);
        gaussCurvatures = VertexData<double>(*mesh);

        for (Vertex v : mesh->vertices()) {
            areas[v] = VertexArea(geom, v);
            meanCurvatures[v] = VertexMeanCurvature(geom, v);
            gaussCurvatures[v] = VertexGaussCurvature(geom, v);
        }
    }

    void LWSClustering::initializeClusters(int numClusters) {
        // Delete all existing clusters
        for (LWSCluster* cluster : clusters) {
            delete cluster;
        }
        clusters.clear();
        claimedVertices.clear();

        // Repopulate with new clusters
        for (int i = 0; i < numClusters; i++) {
            LWSCluster* c = new LWSCluster(mesh, geom);
            c->index = i;
            clusters.push_back(c);
        }

        // Seed each cluster with an initial vertex
        int index = randomInt(0, mesh->nVertices() - 1);
        Vertex v = mesh->vertex(index);
        for (LWSCluster* cluster : clusters) {
            // Rejection sampling until we find a free vertex
            while (claimedVertices.find(v) != claimedVertices.end()) {
                index = randomInt(0, mesh->nVertices() - 1);
                v = mesh->vertex(index);
            }
            // Add it to the vertex
            cluster->vertices.push_back(v);
            claimedVertices.insert(v);
            vertexOwnership[v] = cluster->index;
        }
    }

    void LWSClustering::printClusters() {
        int sumClusters = 0;

        for (LWSCluster* cluster : clusters) {
            cout << "Cluster " << cluster->index << " (" << cluster->a
                << ", " << cluster->b << ", " << cluster->c <<  ") has size " << cluster->vertices.size() << endl;
            sumClusters += cluster->vertices.size();
        }

        cout << "Size of claimed vertices = " << claimedVertices.size() << endl;
        cout << "Sum = " << sumClusters << ", mesh has " << mesh->nVertices() << endl;
    }

    void LWSClustering::colorByClusters() {
        VertexData<Vector3> colors(*mesh);

        for (Vertex v : mesh->vertices()) {
            int i = vertexOwnership[v];
            LWSCluster* cluster = clusters[i];
            Vector3 col{cluster->a, cluster->b, cluster->c};
            col = (col + Vector3{1, 1, 1}) / 2;
            colors[v] = col;
        }

        cout << "Adding color quantity" << endl;
        surface->addVertexColorQuantity("colors", colors);
    }

    void LWSClustering::priorityFillClusters() {
        claimedVertices.clear();
        if (queue) delete queue;
        queue = new std::priority_queue<ClusterQueueEntry, std::vector<ClusterQueueEntry>, CompareEntries>();
        vertexOwnership.fill(-1);
    
        // Reset each cluster to its best (i.e. lowest energy) vertex
        for (LWSCluster* cluster : clusters) {
            Vertex v = cluster->minEnergyVertex();
            cluster->vertices.clear();
            cluster->vertices.push_back(v);
            claimedVertices.insert(v);
            vertexOwnership[v] = cluster->index;
        }

        for (LWSCluster* cluster : clusters) {
            // Seed initial priority queue with neighbors
            for (Vertex v : cluster->vertices) {
                for (Vertex neighbor : v.adjacentVertices()) {
                    // If the neighbor hasn't been claimed already, enqueue it
                    if (claimedVertices.find(neighbor) == claimedVertices.end()) {
                        double energy = cluster->vertexEnergyPre(neighbor, areas, meanCurvatures, gaussCurvatures);
                        ClusterQueueEntry entry{neighbor, cluster, energy};
                        queue->push(entry);
                    }
                }
            }
        }

        // Now pull nodes from the queue and flood-fill according to priorities
        while (claimedVertices.size() < mesh->nVertices()) {
            if (queue->empty()) {
                cout << "Ran out of neighbor pairs before including all vertices" << endl;
                break;
            }

            ClusterQueueEntry entry = queue->top();
            queue->pop();

            // If the popped entry hasn't already been absorbed by another cluster,
            // absorb it now.
            if (claimedVertices.find(entry.vertex) == claimedVertices.end()) {
                entry.cluster->vertices.push_back(entry.vertex);
                claimedVertices.insert(entry.vertex);
                vertexOwnership[entry.vertex] = entry.cluster->index;
                // If any neighbor of the just-added vertex is unclaimed, enqueue it for this cluster.
                for (Vertex neighbor : entry.vertex.adjacentVertices()) {
                    if (claimedVertices.find(neighbor) == claimedVertices.end()) {
                        double energy = entry.cluster->vertexEnergyPre(neighbor, areas, meanCurvatures, gaussCurvatures);
                        ClusterQueueEntry newEntry{neighbor, entry.cluster, energy};
                        queue->push(newEntry);
                    }
                }
            }
        }

        // Every vertex should now be clustered
        // Recompute best coefficients for each cluster
        for (LWSCluster* cluster : clusters) {
            cluster->computeCoefficients(areas, meanCurvatures, gaussCurvatures);
        }

        cout << endl;
        printClusters();
        //colorByClusters();
    }


    // ====================== Individual cluster functions ======================
    
    LWSCluster::LWSCluster() {}

    LWSCluster::LWSCluster(HalfedgeMesh* m, VertexPositionGeometry* g) {
        mesh = m;
        geom = g;
    }

    double LWSCluster::vertexEnergy(Vertex v) {
        double A = VertexArea(geom, v);
        double H = VertexMeanCurvature(geom, v);
        double K = VertexGaussCurvature(geom, v);

        return a * A + b * H + c * K;
    }

    double LWSCluster::vertexEnergyPre(Vertex v, VertexData<double> &As, VertexData<double> &Hs, VertexData<double> &Ks) {
        double A = As[v];
        double H = Hs[v];
        double K = Ks[v];

        return a * A + b * H + c * K;
    }

    Vertex LWSCluster::minEnergyVertex() {
        if (vertices.size() == 0) return Vertex();
        
        Vertex minV = vertices[0];
        double minE = vertexEnergy(minV);
        for (Vertex v : vertices) {
            double e = vertexEnergy(v);
            if (e < minE) {
                minV = v;
                minE = e;
            }
        }

        return minV;
    }

    Eigen::Matrix3d LWSCluster::clusterQuadric(VertexData<double> &As, VertexData<double> &Hs, VertexData<double> &Ks) {
        Eigen::Matrix3d M;
        M << 0, 0, 0,
            0, 0, 0,
            0, 0, 0;

        for (Vertex v : vertices) {
            double A = As[v];
            double H = Hs[v];
            double K = Ks[v];

            double H2 = H * H;
            double K2 = K * K;
            double HK = H * K;

            M(0, 0) += A; M(0, 1) += H; M(0, 2) += K;
            M(1, 0) += H; M(1, 1) += H2 / A; M(1, 2) += HK / A;
            M(2, 0) += K; M(2, 1) += HK / A; M(2, 2) += K2 / A;
        }

        //cout << M << endl;

        return M;
    }

    void LWSCluster::computeCoefficients(VertexData<double> &As, VertexData<double> &Hs, VertexData<double> &Ks) {
        Eigen::Matrix3d quadric = clusterQuadric(As, Hs, Ks);
        Eigen::EigenSolver<Eigen::Matrix3d> solver(quadric, true);

        auto evals = solver.eigenvalues();
        auto evecs = solver.eigenvectors();

        int lowest = 0;
        double lowestVal = evals(0).real();
        
        // Find the smallest eigenvalue
        for (int i = 1; i < 3; i++) {
            if (evals(i).real() < lowestVal) {
                lowest = i;
                lowestVal = evals(i).real();
            }
        }

        //cout << "Eigenvalues are " << evals << endl;

        // Get the coresponding eigenvector
        auto col = evecs.col(lowest);
        a = col(0).real();
        b = col(1).real();
        c = col(2).real();

        //cout << "Lowest eigenvector is " << col << ", corresponding to eigenvalue " << evals(lowest) << endl;
    }
}
