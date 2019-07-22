#include "lws_flow.h"

#include <iostream>
#include <vector>
#include <fstream>
#include "lws_options.h"
#include "mesh_helpers.h"
#include "vertexderivatives.h"

#include <chrono>

using namespace std;
using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace LWS {

    GradientSolver::GradientSolver() {
        cholmod_start(context);

        laplacian = 0;
        bilaplacian = 0;
        saddleMatrix = 0;
        symbolicFactor = 0;
        vertexGradients = 0;
        transformedGradients = 0;

        originalPositions = 0;

        workspaceE = 0;
        workspaceY = 0;
        doublingLimit = 1;

        coeffA = 1;
        coeffB = -2;
        coeffC = 1;

        workspaceE = 0;
        workspaceY = 0;

        context.setSimplicial();

        cout << "Started cholmod context for solver" << endl;
    }

    GradientSolver::~GradientSolver() {
        cholmod_finish(context);
        cout << "Finished cholmod context for solver" << endl;
    }

    CholmodContext* GradientSolver::GetContext() {
        return &context;
    }

    double vectorNorm(cholmod_dense* vec) {
        double sum = 0;
        for (size_t i = 0; i < vec->nrow; i++) {
            double entry = ((double*)vec->x)[i];
            sum += entry * entry;
        }
        return sqrt(sum);
    }

    Vector3 vertexValue(cholmod_dense* vec, int vertIndex) {
        int base = 3 * vertIndex;
        double x = ((double*)vec->x)[base + 0];
        double y = ((double*)vec->x)[base + 1];
        double z = ((double*)vec->x)[base + 2];
        return Vector3{x, y, z};
    }

    void setVertexValue(cholmod_dense* vec, int vertIndex, Vector3 value) {
        int base = 3 * vertIndex;
        ((double*)vec->x)[base + 0] = value.x;
        ((double*)vec->x)[base + 1] = value.y;
        ((double*)vec->x)[base + 2] = value.z;
    }

    void GradientSolver::setSurface(HalfedgeMesh* m, VertexPositionGeometry* g) {
        mesh = m;
        geom = g;
        if (originalPositions) delete originalPositions;
        originalPositions = new VertexData<Vector3>(*mesh);
    }

    void GradientSolver::setEnergy(LWSEnergyFunction* lws) {
        lwsFunc = lws;
    }

    void addTriplet(cholmod_triplet* M, int rowInd, int colInd, double value) {
        // Add to each of the row/col/value arrays
        ((int*)M->i)[M->nnz] = rowInd;
        ((int*)M->j)[M->nnz] = colInd;
        ((double*)M->x)[M->nnz] = value;
        // Increment the number of nonzeroes
        M->nnz++;
    }

    void GradientSolver::printMatrix(cholmod_dense* matrix, std::string name) {
        ofstream file;
        file.open(name);
        for (size_t i = 0; i < matrix->nrow; i++) {
            std::string line = "";
            for (size_t j = 0; j < matrix->ncol; j++) {
                size_t index = i * matrix->ncol + j;
                if (j < matrix->ncol - 1) {
                    line += to_string(((double*)matrix->x)[index]) + ", ";
                }
                else {
                    line += to_string(((double*)matrix->x)[index]);
                }
            }
            file << line << endl;
        }
        file.close();
        cout << "Wrote matrix to " << name << endl;
    }

    void GradientSolver::printMatrix(cholmod_sparse* matrix, std::string name) {
        cholmod_dense* A = cholmod_sparse_to_dense(matrix, context);
        printMatrix(A, name);
    }

    void GradientSolver::computeCotanLaplacian() {
        if (laplacian) {
            cholmod_free_sparse(&laplacian, context);
        }
        if (bilaplacian) {
            cholmod_free_sparse(&bilaplacian, context);
        }

        // Allocate a sparse triplet matrix
        int numVertices = mesh->nVertices();
        int nzmax = 12 * numVertices;
        int matrixType = 0;
        // TODO: better max number of non-zeroes?
        cholmod_triplet* L_triplets = cholmod_allocate_triplet(numVertices, numVertices, nzmax,
            matrixType, CHOLMOD_REAL, context);
        
        cholmod_triplet* Minv_L_triplets = cholmod_allocate_triplet(numVertices, numVertices, nzmax,
            matrixType, CHOLMOD_REAL, context);

        VertexData<size_t> indices = mesh->getVertexIndices();

        // Build the Laplacian row by row
        for (Vertex v : mesh->vertices()) {
            double invDualWeight = 1.0 / VertexArea(geom, v);
            int rowInd = indices[v];
            double rowSum = 0;

            // Laplacian boundary conditions
            if (v.isBoundary()) {
                addTriplet(L_triplets, rowInd, rowInd, 1);
                continue;
            }

            // Add up the cotan weights for all neighbors
            for (Halfedge he : v.outgoingHalfedges()) {
                // Skip boundary edges, since one of the faces is imaginary
                if (he.edge().isBoundary()) {
                    continue;
                }
                Vertex neighbor = he.twin().vertex();
                int colInd = indices[neighbor];

                double cotWeight = 0.5 * (geom->halfedgeCotanWeight(he) + geom->halfedgeCotanWeight(he.twin()));
                if (matrixType == 0 || colInd < rowInd) {
                    addTriplet(L_triplets, rowInd, colInd, -cotWeight);
                    addTriplet(Minv_L_triplets, rowInd, colInd, -cotWeight * invDualWeight);
                }
                // Accumulate the row sum
                rowSum += cotWeight;
            }
            addTriplet(L_triplets, rowInd, rowInd, rowSum);
        }

        cholmod_sparse* L = cholmod_triplet_to_sparse(L_triplets, nzmax, context);
        laplacian = L;
        cholmod_sparse* Minv_L = cholmod_triplet_to_sparse(Minv_L_triplets, nzmax, context);
        bilaplacian = cholmod_ssmult(L, Minv_L, 0, 1, 1, context);

        cholmod_free_triplet(&L_triplets, context);
        cholmod_free_triplet(&Minv_L_triplets, context);

        cholmod_free_sparse(&Minv_L, context);
    }

    EigenSparse GradientSolver::CotanLaplacianEigen() {
        VertexData<size_t> indices = mesh->getVertexIndices();

        // TODO: assemble cotan laplacian in an Eigen sparse matrix
        std::vector<EigenTriplet> triplets;
        int numVertices = mesh->nVertices();
        triplets.reserve(numVertices * 6);

        for (Vertex v : mesh->vertices()) {
            int rowInd = indices[v];

            if (v.isBoundary()) {
                triplets.push_back(EigenTriplet(rowInd, rowInd, 1.0));
                continue;
            }

            double rowSum = 0;
            for (Halfedge he : v.outgoingHalfedges()) {
                Vertex neighbor = he.twin().vertex();
                int colInd = indices[neighbor];
                double cotWeight = 0.5 * (geom->halfedgeCotanWeight(he) + geom->halfedgeCotanWeight(he.twin()));

                triplets.push_back(EigenTriplet(rowInd, colInd, -cotWeight));
                rowSum += cotWeight;
            }

            triplets.push_back(EigenTriplet(rowInd, rowInd, rowSum));
        }

        EigenSparse mat(numVertices, numVertices);
        mat.setFromTriplets(triplets.begin(), triplets.end());
        return mat;
    }

    void GradientSolver::computeSaddleMatrix() {
        if (saddleMatrix) {
            cholmod_free_sparse(&saddleMatrix, context);
        }

        VertexData<size_t> indices = mesh->getVertexIndices();

        // Get the triplets back
        cholmod_triplet* J_triplets = cholmod_sparse_to_triplet(bilaplacian, context);
        // Adding the barycenter constraint expands the matrix by 3 in each dimension
        int rowBase = 3 * J_triplets->nrow;
        int nRows = 3 * J_triplets->nrow + 3;
        // We need to expand the laplacian, replacing every entry with 3x3 diagonal blocks
        // We also need room for the barycenter rows and columns at the end.
        int newMax = 3 * J_triplets->nzmax + mesh->nVertices() * 6;

        const int SYMMETRY = 1;

        // Make a triplet matrix storing only the upper triangle
        cholmod_triplet* A_triplets = cholmod_allocate_triplet(nRows, nRows, newMax,
            SYMMETRY, CHOLMOD_REAL, context);
        
        for (size_t i = 0; i < J_triplets->nnz; i++) {
            int r = ((int*)J_triplets->i)[i];
            int c = ((int*)J_triplets->j)[i];
            double val = ((double*)J_triplets->x)[i];
            // Add only the part above the diagonal, meaning row index <= column index
            if (r <= c) {
                // Every entry goes to a 3x3 diagonal block with the value duplicated on the diagonal
                addTriplet(A_triplets, 3 * r + 0, 3 * c + 0, val);
                addTriplet(A_triplets, 3 * r + 1, 3 * c + 1, val);
                addTriplet(A_triplets, 3 * r + 2, 3 * c + 2, val);
            }
        }

        // Sum the total area of the surface
        double sumArea = 0;
        for (Face f : mesh->faces()) {
            sumArea += geom->faceArea(f);
        }

        // Add rows and columns of barycenter constraint
        for (Vertex v : mesh->vertices()) {
            // The gradient of the barycenter w.r.t. any vertex is just the all-ones vector,
            // multiplied by the proportion of the total area affected by that vertex,
            // i.e. the dual area.
            double areaWeight = VertexArea(geom, v) / sumArea;
            int vBase = 3 * indices[v];

            // Fill the bottom-left part based on symmetry
            if (SYMMETRY <= 0) {
                addTriplet(A_triplets, rowBase + 0, vBase + 0, areaWeight);
                addTriplet(A_triplets, rowBase + 1, vBase + 1, areaWeight);
                addTriplet(A_triplets, rowBase + 2, vBase + 2, areaWeight);
            }
            // Fill the top-right part based on symmetry
            if (SYMMETRY >= 0) {
                addTriplet(A_triplets, vBase + 0, rowBase + 0, areaWeight);
                addTriplet(A_triplets, vBase + 1, rowBase + 1, areaWeight);
                addTriplet(A_triplets, vBase + 2, rowBase + 2, areaWeight);
            }
        }

        saddleMatrix = cholmod_triplet_to_sparse(A_triplets, newMax, context);

        cholmod_free_triplet(&J_triplets, context);
        cholmod_free_triplet(&A_triplets, context);
    }

    void GradientSolver::symbolicFactorization(cholmod_sparse* matrix) {
        if (symbolicFactor) {
            cholmod_free_factor(&symbolicFactor, context);
        }
        symbolicFactor = cholmod_analyze(matrix, context);
        cout << "Computed symbolic factorization" << endl;
    }

    void GradientSolver::numericFactorization(cholmod_sparse* matrix) {
        cholmod_factorize(matrix, symbolicFactor, context);
    }

    Vector3 GradientSolver::vertexGradient(GradientFunction gradFunc, Vertex v) {
        // Start with gradient of self wrt v
        Vector3 sumGradient = gradFunc(geom, v, v);
        // Add up gradient of all neighbors wrt v
        for (Vertex neighbor : v.adjacentVertices()) {
            sumGradient += gradFunc(geom, neighbor, v);
        }
        return sumGradient;
    }

    void GradientSolver::computeVertexGradients(GradientFunction gradFunc) {
        VertexData<size_t> indices = mesh->getVertexIndices();

        int nRows = 3 * mesh->nVertices() + 3;
        
        // Only reallocate the vertex gradients if they are unallocated
        if (!vertexGradients)
            vertexGradients = cholmod_zeros(nRows, 1, CHOLMOD_REAL, context);

        for (Vertex v : mesh->vertices()) {
            int vBase = indices[v] * 3;
            Vector3 sumGradient = vertexGradient(gradFunc, v);

            ((double*)vertexGradients->x)[vBase + 0] = sumGradient.x;
            ((double*)vertexGradients->x)[vBase + 1] = sumGradient.y;
            ((double*)vertexGradients->x)[vBase + 2] = sumGradient.z;
        }
    }

    cholmod_dense* GradientSolver::transformGradient(cholmod_dense* rhs) {
        cholmod_solve2(CHOLMOD_A, symbolicFactor, rhs, 0, &transformedGradients, 0, &workspaceY, &workspaceE, context);
        return transformedGradients;
    }

    cholmod_dense* GradientSolver::copyGradient(cholmod_dense* rhs) {
        cholmod_copy_dense2(rhs, transformedGradients, context);
        return transformedGradients;
    }

    void GradientSolver::projectOntoNormal(cholmod_dense* gradient) {
        VertexData<size_t> indices = mesh->getVertexIndices();

        for (Vertex v : mesh->vertices()) {
            Vector3 current = vertexValue(gradient, indices[v]);
            Vector3 normal = geom->vertexNormals[v];
            current = dot(current, normal) * normal;
            setVertexValue(gradient, indices[v], current);
        }
    }

    void GradientSolver::backprojectConstraints(cholmod_factor* factor) {
        VertexData<size_t> indices = mesh->getVertexIndices();

        Vector3 barycenter = Barycenter(mesh, geom);

        int nRows = 3 * mesh->nVertices() + 3;
        // Reuse the already-allocated vertex gradients matrix
        cholmod_dense* data = vertexGradients;

        for (int i = 0; i < nRows - 3; i++) {
            ((double*)data->x)[i] = 0;
        }
        
        ((double*)data->x)[nRows - 3] = -barycenter.x;
        ((double*)data->x)[nRows - 2] = -barycenter.y;
        ((double*)data->x)[nRows - 1] = -barycenter.z;

        cholmod_solve2(CHOLMOD_A, factor, data, 0, &transformedGradients, 0, &workspaceY, &workspaceE, context);

        for (Vertex v : mesh->vertices()) {
            Vector3 correction = vertexValue(transformedGradients, indices[v]);
            geom->vertexPositions[v] = geom->vertexPositions[v] + correction;
        }
    }

    void GradientSolver::FixedStep(cholmod_dense* gradient, double delta) {
        VertexData<size_t> indices = mesh->getVertexIndices();
        // Take a fixed step of size delta along the gradient
        for (Vertex v : mesh->vertices()) {
            int vBase = indices[v] * 3;
            Vector3 offset = vertexValue(gradient, indices[v]);
            geom->vertexPositions[v] = geom->vertexPositions[v] - delta * offset;
        }
    }

    double GradientSolver::EvaluateGradientStep(EnergyFunction energy, cholmod_dense* gradient, double delta) {
        // Write the new vertex positions to the mesh
        VertexData<size_t> indices = mesh->getVertexIndices();
        // Step every vertex by the gradient times delta
        for (Vertex v : mesh->vertices()) {
            Vector3 vertGrad = vertexValue(gradient, indices[v]);
            geom->vertexPositions[v] = (*originalPositions)[v] - delta * vertGrad;
        }
        // Compute the new energy
        return MeshEnergy(geom, mesh, energy);
    }

    bool GradientSolver::LineSearchStep(EnergyFunction energy, cholmod_dense* gradient, double initGuess) {
        VertexData<size_t> indices = mesh->getVertexIndices();

        double delta = initGuess;

        // Save initial positions
        for (Vertex v : mesh->vertices()) {
            (*originalPositions)[v] = geom->vertexPositions[v];
            Vector3 grad = vertexValue(gradient, indices[v]);
        }

        double initialEnergy = MeshEnergy(geom, mesh, energy);
        double gradNorm = vectorNorm(gradient);
        int numBacktracks = 0, numDoubles = 0;
        double sigma = 0.01f;
        double newEnergy = initialEnergy;

        while (delta > 1e-10) {
            newEnergy = EvaluateGradientStep(energy, gradient, delta);
            double decrease = initialEnergy - newEnergy;
            double targetDecrease = sigma * initialEnergy * delta * gradNorm;

            // If the energy hasn't decreased enough to meet the Armijo condition,
            // halve the step size.
            if (decrease < targetDecrease) {
                delta /= 2;
                numBacktracks++;
            }
            else if (decrease >= targetDecrease && numBacktracks == 0 && numDoubles < doublingLimit) {
                delta *= 2;
                numDoubles++;
            }
            // Otherwise, accept the current step.
            else {
                cout << "Took step of size " << delta << " after " << numBacktracks
                    << " backtracks, " << numDoubles << " doubles" << endl;
                break;
            }
        }

        if (delta <= 1e-10) {
            cout << "Failed to find a non-trivial step after " << numBacktracks << " backtracks" << endl;
            for (Vertex v : mesh->vertices()) {
                // Restore initial positions if step size goes to 0
                geom->vertexPositions[v] = (*originalPositions)[v];
            }
            return false;
        }
        else {
            cout << "Energy: " << initialEnergy << " -> " << newEnergy << endl;
            return true;
        }
    }
    
    void GradientSolver::stepProjectedGradient(EnergyFunction energyFunc, GradientFunction gradFunc, double h) {
        chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();

        // Evaluate the new gradient
        lwsFunc->PrecomputeEnergies(mesh, geom);
        computeVertexGradients(gradFunc);

        if (LWSOptions::useSobalev) {
            // Reassemble the saddle matrix after each step
            computeCotanLaplacian();
            computeSaddleMatrix();

            // We can reuse the symbolic factorization from before.
            numericFactorization(saddleMatrix);
            transformGradient(vertexGradients);

            // Take a step
            double gradNorm = vectorNorm(transformedGradients);
            if (gradNorm > 0) {
                LineSearchStep(energyFunc, transformedGradients, 1.0 / sqrt(gradNorm));
                // Correct for any drift with backprojection
                backprojectConstraints(symbolicFactor);
            }
            else {
                cout << "Gradient is 0" << endl;
            }
        }

        else {
            // Just use the L2 gradient
            double gradNorm = vectorNorm(vertexGradients);
            LineSearchStep(energyFunc, vertexGradients, 1.0 / gradNorm);
            // No backprojection in this case, since we never computed the saddle matrix
        }

        //test_quantities(surface);

        chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();

        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        double fps = 1000.0 / duration;

        cout << "Took time step in " << duration << " ms (" << fps << " FPS)" << endl;
        cout << endl;
    }

    void print_deriv_comparison(VertexPositionGeometry* geom, Vertex v1, Vertex v2) {
        Vector3 areaAnalytic = GradientArea(geom, v1, v2);
        Vector3 meanAnalytic = GradientMeanCurvature(geom, v1, v2);
        Vector3 gaussAnalytic = GradientGaussCurvature(geom, v1, v2);

        Vector3 areaNumerical = GradientAreaNumerical(geom, v1, v2);
        Vector3 meanNumerical = GradientMeanNumerical(geom, v1, v2);
        Vector3 gaussNumerical = GradientGaussNumerical(geom, v1, v2);

        double areaScale = norm(areaAnalytic) / norm(areaNumerical);
        double meanScale = norm(meanAnalytic) / norm(meanNumerical);
        double gaussScale = norm(gaussAnalytic) / norm(gaussNumerical);

        double areaDot = dot(areaAnalytic, areaNumerical) / (norm(areaAnalytic) * norm(areaNumerical));
        double meanDot = dot(meanAnalytic, meanNumerical) / (norm(meanAnalytic) * norm(meanNumerical));
        double gaussDot = dot(gaussAnalytic, gaussNumerical) / (norm(gaussAnalytic) * norm(gaussNumerical));

        cout << "Area:  " << areaAnalytic << " (dot product " << areaDot << ")" << endl;
        cout << "       " << areaNumerical << " (magnitude ratio " << areaScale << ")" << endl;

        cout << "Mean:  " << meanAnalytic << " (dot product " << meanDot << ")" << endl;
        cout << "       " << meanNumerical << " (magnitude ratio " << meanScale << ")" << endl;
        
        cout << "Gauss: " << gaussAnalytic << " (dot product " << gaussDot << ")" << endl;
        cout << "       " << gaussNumerical << " (magnitude ratio " << gaussScale << ")\n" << endl;
    }

    void test_derivatives(HalfedgeMesh* mesh, VertexPositionGeometry* geom) {
        double totalK = 0;
        
        for (Vertex v : mesh->vertices()) {
            for (Vertex v2 : v.adjacentVertices()) {
                print_deriv_comparison(geom, v, v2);
            }
        }
    }

    void test_quantities(HalfedgeMesh* mesh, VertexPositionGeometry* geom) {
        double totalH = 0;
        double totalK = 0;
        
        for (Vertex v : mesh->vertices()) {
            double h = VertexMeanCurvature(geom, v);
            totalH += h;
            totalK += VertexGaussCurvature(geom, v);
        }

        cout << "Total mean curvature = " << totalH << endl;
        cout << "Total Gauss curvature = " << totalK << endl;
    }

    void test_dihedral(HalfedgeMesh* mesh, VertexPositionGeometry* geom) {
        for (Edge e : mesh->edges()) {
            Vertex p1 = e.halfedge().vertex();
            Vertex p2 = e.halfedge().twin().vertex();
            Vertex p3 = e.halfedge().next().next().vertex();
            Vertex p4 = e.halfedge().twin().next().next().vertex();

            Vector3 deriv1a = GradientDihedralAngle(geom, e, p1);
            Vector3 deriv2a = GradientDihedralAngle(geom, e, p2);
            Vector3 deriv3a = GradientDihedralAngle(geom, e, p3);
            Vector3 deriv4a = GradientDihedralAngle(geom, e, p4);

            Vector3 derivsA[4] = { deriv1a, deriv2a, deriv3a, deriv4a };

            Vector3 deriv1n = GradientDihedralNumerical(geom, e, p1, 0.001);
            Vector3 deriv2n = GradientDihedralNumerical(geom, e, p2, 0.001);
            Vector3 deriv3n = GradientDihedralNumerical(geom, e, p3, 0.001);
            Vector3 deriv4n = GradientDihedralNumerical(geom, e, p4, 0.001);

            Vector3 derivsN[4] = { deriv1n, deriv2n, deriv3n, deriv4n };

            for (int i = 0; i < 4; i++) {
                double scaleDiff = norm(derivsA[i]) / norm(derivsN[i]);
                double angleDiff = dot(derivsA[i], derivsN[i]) / (norm(derivsA[i]) * norm(derivsN[i]));

                if (abs(angleDiff - 1) > 0.01 && abs(scaleDiff - 1) > 0.01) {
                    cout << "Analytic: " << derivsA[i] << endl;
                    cout << "Numeric:  " << derivsN[i] << endl;
                    cout << "    Dot product = " << angleDiff << ", scale diff = " << scaleDiff << "\n" << endl;
                }
            }
        }
    }

    void step_vertex_flow(HalfedgeMesh* mesh, VertexPositionGeometry* geom, EnergyFunction energyFunc, GradientFunction gradFunc, double h) {
        VertexData<size_t> indices = mesh->getVertexIndices();

        VertexData<Vector3> gradients(*mesh, Vector3{0, 0, 0});

        for (Vertex v : mesh->vertices()) {
            // For each vertex v, we need to add up the gradient terms:
            // gradient of own energy wrt v, and gradient of all neighbors' energies wrt v.
            // This gives the gradient of the global energy wrt v.

            // Start with gradient of self wrt v
            Vector3 sumGradient = energyFunc(geom, v) * gradFunc(geom, v, v);
            // Add up gradient of all neighbors wrt v
            for (Vertex neighbor : v.adjacentVertices()) {
                sumGradient += energyFunc(geom, neighbor) * gradFunc(geom, neighbor, v);
            }
            gradients[v] = sumGradient;
        }

        // Now we displace the mesh by h * vertex gradients.
        for (Vertex v : mesh->vertices()) {
            geom->vertexPositions[v] = geom->vertexPositions[v] - h * gradients[v];
        }
    }
}