#include "boundary_loop.h"
#include "vertexderivatives.h"
#include "boundary_derivatives.h"

#include "lws_flow.h"

using namespace std;

namespace LWS {

    using namespace geometrycentral::surface;

    LWSBoundaryLoop::LWSBoundaryLoop(surface::HalfedgeMesh* m, surface::VertexPositionGeometry* g) {
        cholmod_start(context);
        mesh = m;
        geom = g;

        VertexData<size_t> indices = m->getVertexIndices();

        if (mesh->nBoundaryLoops() > 0) {
            boundary = mesh->boundaryLoop(0);
            Halfedge start = boundary.halfedge();
            Halfedge he = boundary.halfedge();

            do {
                vertices.push_back(he.vertex());
                size_t global_i = indices[he.vertex()];
                loopToGlobal.push_back(global_i);
                he = he.next();
            }
            while (he != start);

            valid = true;

            for (Vertex v : boundary.adjacentVertices()) {
                double sum = 2 * M_PI;
                for (Corner c : v.adjacentCorners()) {
                    if (c.face().isBoundaryLoop()) {
                        sum -= geom->cornerAngle(c);
                    }
                }
                interiorAngles.push_back(sum / (2 * M_PI));
            }

            ReparameterizeArc();
            cout << "Extracted boundary loop of " << NumVertices() << " vertices" << endl;
        }
        else {
            cout << "Surface does not have boundary" << endl;
            valid = false;
        }

        gradients = 0;
        workspaceE = 0;
        workspaceY = 0;
        symbolicFactor = 0;
        boundary_bilaplacian = 0;
        boundary_saddle = 0;
        fullRHS = 0;

        CreateNodesAndElements();
    }

    bool LWSBoundaryLoop::isValid() {
        return valid;
    }

    Vertex LWSBoundaryLoop::GetVertex(int index) {
        return vertices[index];
    }

    BoundaryElement* LWSBoundaryLoop::GetElement(int index) {
        return bem_elements[index];
    }

    size_t LWSBoundaryLoop::GlobalIndex(int index) {
        return loopToGlobal[index];
    }

    int LWSBoundaryLoop::PrevVertex(int index) {
        if (index == 0) {
            return NumVertices() - 1;
        }
        else return index - 1;
    }
    
    int LWSBoundaryLoop::NextVertex(int index) {
        if (index == NumVertices() - 1) {
            return 0;
        }
        else return index + 1;
    }

    Vector3 LWSBoundaryLoop::Position(int index) {
        return geom->vertexPositions[GetVertex(index)];
    }

    void LWSBoundaryLoop::PrintVertices2D() {
        for (int i = 0; i < NumVertices(); i++) {
            Vector3 p = Position(i);
            std::cout << p.x << ", " << p.z << std::endl;
        }
    }

    Vector3 LWSBoundaryLoop::BisectorNormal(int index) {
        Vertex center = GetVertex(index);
        Vertex prev = GetVertex(PrevVertex(index));
        Vertex next = GetVertex(NextVertex(index));

        Vector3 v1 = geom->vertexPositions[prev];
        Vector3 v2 = geom->vertexPositions[center];
        Vector3 v3 = geom->vertexPositions[next];

        Vector3 e1 = (v2 - v1);
        Vector3 e2 = (v3 - v2);

        Vector3 up{0, 1, 0};

        Vector3 n1 = cross(up, e1);
        n1.normalize();
        Vector3 n2 = cross(up, e2);
        n2.normalize();

        // Return the angle bisector normal
        Vector3 normal = n1 + n2;
        normal.normalize();
        return normal;
    }

    Vector3 LWSBoundaryLoop::VertexTangent(int index) {
        Vector3 prev = Position(PrevVertex(index));
        Vector3 next = Position(NextVertex(index));
        Vector3 center = Position(index);

        Vector3 prevT = center - prev;
        Vector3 nextT = next - center;
        prevT.normalize();
        nextT.normalize();

        Vector3 tangent = prevT + nextT;
        tangent.normalize();
        return tangent;
    }

    Vector3 LWSBoundaryLoop::LengthWeightedNormal(int index) {
        Vertex center = GetVertex(index);
        Vertex prev = GetVertex(PrevVertex(index));
        Vertex next = GetVertex(NextVertex(index));

        Vector3 v1 = geom->vertexPositions[prev];
        Vector3 v2 = geom->vertexPositions[center];
        Vector3 v3 = geom->vertexPositions[next];

        Vector3 e1 = (v2 - v1);
        Vector3 e2 = (v3 - v2);
        double w1 = norm(e1);
        double w2 = norm(e2);

        Vector3 up{0, 1, 0};

        Vector3 n1 = cross(up, e1);
        n1.normalize();
        Vector3 n2 = cross(up, e2);
        n2.normalize();

        // Weight normals by length
        Vector3 normal = (w1 * n1 + w2 * n2) / (w1 + w2);
        normal.normalize();
        return normal;
    }

    Vector3 LWSBoundaryLoop::EdgeNormal(int index) {
        Vertex center = GetVertex(index);
        Vertex next = GetVertex(NextVertex(index));

        Vector3 edgeVec = geom->vertexPositions[next] - geom->vertexPositions[center];
        Vector3 up{0, 1, 0};
        Vector3 normal = cross(up, edgeVec);
        normal.normalize();
        return normal;
    }

    double LWSBoundaryLoop::EdgeLength(int index) {
        Vertex center = GetVertex(index);
        Vertex next = GetVertex(NextVertex(index));
        Vector3 edgeVec = geom->vertexPositions[next] - geom->vertexPositions[center];
        return norm(edgeVec);
    }

    double LWSBoundaryLoop::LoopCurvature(int index) {
        Vertex center = GetVertex(index);
        Vertex prev = GetVertex(PrevVertex(index));
        Vertex next = GetVertex(NextVertex(index));

        Vector3 v1 = geom->vertexPositions[prev];
        Vector3 v2 = geom->vertexPositions[center];
        Vector3 v3 = geom->vertexPositions[next];
        
        Vector3 e1 = (v2 - v1);
        Vector3 e2 = (v3 - v2);
        e1.normalize();
        e2.normalize();
        Vector3 up{0, 1, 0};
        double sinT = dot(up, cross(e1, e2));
        return asin(sinT);
    }
    
    void LWSBoundaryLoop::MoveAlongNormal(double h) {
        for (int i = 0; i < NumVertices(); i++) {
            geom->vertexPositions[vertices[i]] += h * BisectorNormal(i);
        }
    }

    void LWSBoundaryLoop::ScaleByHalf() {
        for (int i = 0; i < NumVertices(); i++) {
            Vector3 p = geom->vertexPositions[vertices[i]];
            geom->vertexPositions[vertices[i]] = p / 2;
        }
    }
    
    int LWSBoundaryLoop::NumVertices() {
        return vertices.size();
    }

    void LWSBoundaryLoop::StepBoundaryGradient(double areaCoeff, double lenCoeff, double h) {
        // Assemble required matrices
        BoundaryLaplacian();
        BoundarySaddle();
        // Compute and step gradients
        ComputeGradients(areaCoeff, lenCoeff);
        StepGradients(h);
    }

    double VertexMass(Vector3 centerPos, Vector3 prevPos, Vector3 nextPos) {
        double length1 = norm(centerPos - prevPos);
        double length2 = norm(centerPos - nextPos);
        double invMass = 1 / ((length1 + length2) / 2);
        return invMass;
    }

    void LWSBoundaryLoop::BoundaryLaplacian() {

        if (boundary_bilaplacian) {
            cholmod_free_sparse(&boundary_bilaplacian, context);
        }

        int numVerts = NumVertices();
        int nzmax = 3 * numVerts + 3;

        cholmod_triplet* Minv_triplets = cholmod_allocate_triplet(numVerts, numVerts, nzmax,
            0, CHOLMOD_REAL, context);
        cholmod_triplet* L_triplets = cholmod_allocate_triplet(numVerts, numVerts, nzmax,
            0, CHOLMOD_REAL, context);

        for (int row = 0; row < numVerts; row++) {
            int prevRow = (row - 1 + numVerts) % numVerts;
            int nextRow = (row + 1 + numVerts) % numVerts;

            // Add the three entries in the row
            addTriplet(L_triplets, row, prevRow, -1);
            addTriplet(L_triplets, row, nextRow, -1);
            addTriplet(L_triplets, row, row, 2);

            // Compute mass matrix weights
            Vector3 centerPos = geom->vertexPositions[GetVertex(row)];
            Vector3 prevPos = geom->vertexPositions[GetVertex(prevRow)];
            Vector3 nextPos = geom->vertexPositions[GetVertex(nextRow)];

            double invMass = VertexMass(centerPos, prevPos, nextPos);

            // Add the same three entries to the weighted matrix
            addTriplet(Minv_triplets, row, prevRow, -invMass);
            addTriplet(Minv_triplets, row, nextRow, -invMass);
            addTriplet(Minv_triplets, row, row, 2 * invMass);
        }

        cholmod_sparse* L = cholmod_triplet_to_sparse(L_triplets, nzmax, context);
        cholmod_sparse* Minv_L = cholmod_triplet_to_sparse(Minv_triplets, nzmax, context);

        boundary_bilaplacian = cholmod_ssmult(L, Minv_L, 0, 1, 1, context);

        cholmod_free_sparse(&L, context);
        cholmod_free_sparse(&Minv_L, context);
        cholmod_free_triplet(&Minv_triplets, context);
        cholmod_free_triplet(&L_triplets, context);
    }

    void LWSBoundaryLoop::BoundarySaddle() {
        if (boundary_saddle) {
            cholmod_free_sparse(&boundary_saddle, context);
        }
        // Get the triplets back
        cholmod_triplet* J_triplets = cholmod_sparse_to_triplet(boundary_bilaplacian, context);
        int numVerts = NumVertices();
        int nRows = NumVertices() * 3;
        int newMax = 3 * J_triplets->nzmax + nRows * 2;

        const int SYMMETRY = 1;

        // Make a triplet matrix storing only the upper triangle
        cholmod_triplet* A_triplets = cholmod_allocate_triplet(nRows, nRows, newMax,
            SYMMETRY, CHOLMOD_REAL, context);
        
        for (size_t i = 0; i < J_triplets->nnz; i++) {
            int r = ((int*)J_triplets->i)[i];
            int c = ((int*)J_triplets->j)[i];
            double val = ((double*)J_triplets->x)[i];
            // Add a little of the mass matrix to the diagonal
            if (r == c) {
                int prevRow = (r - 1 + numVerts) % numVerts;
                int nextRow = (r + 1 + numVerts) % numVerts;

                Vector3 centerPos = geom->vertexPositions[GetVertex(r)];
                Vector3 prevPos = geom->vertexPositions[GetVertex(prevRow)];
                Vector3 nextPos = geom->vertexPositions[GetVertex(nextRow)];

                double mass = VertexMass(centerPos, prevPos, nextPos);
                val += 0.001 * mass;
            }

            // Add only the part above the diagonal, meaning row index <= column index
            if (r <= c) {
                // Every entry goes to a 3x3 diagonal block with the value duplicated on the diagonal
                addTriplet(A_triplets, 3 * r + 0, 3 * c + 0, val);
                addTriplet(A_triplets, 3 * r + 1, 3 * c + 1, val);
                addTriplet(A_triplets, 3 * r + 2, 3 * c + 2, val);
            }
        }

        boundary_saddle = cholmod_triplet_to_sparse(A_triplets, newMax, context);

        cholmod_free_triplet(&J_triplets, context);
        cholmod_free_triplet(&A_triplets, context);
    }

    void LWSBoundaryLoop::symbolicFactorization(cholmod_sparse* matrix) {
        if (symbolicFactor) {
            cholmod_free_factor(&symbolicFactor, context);
        }
        symbolicFactor = cholmod_analyze(matrix, context);
    }

    void LWSBoundaryLoop::numericFactorization(cholmod_sparse* matrix) {
        cholmod_factorize(matrix, symbolicFactor, context);
    }

    void LWSBoundaryLoop::ComputeGradients(double areaCoeff, double bLengthCoeff) {
        int nRows = 3 * NumVertices();

        cholmod_dense* temp_gradients = cholmod_zeros(nRows, 1, CHOLMOD_REAL, context);
        
        for (int i = 0; i < NumVertices(); i++) {
            Vertex vert = GetVertex(i);
            Vector3 areaGrad = GradientArea(geom, vert, vert);
            Vector3 lengthGrad = GradientBLength(geom, vert, vert);

            Vector3 sumGrad = areaCoeff * areaGrad + bLengthCoeff * lengthGrad;

            int baseIndex = 3 * i;
            ((double*)temp_gradients->x)[baseIndex + 0] = sumGrad.x;
            ((double*)temp_gradients->x)[baseIndex + 1] = sumGrad.y;
            ((double*)temp_gradients->x)[baseIndex + 2] = sumGrad.z;
        }
        /*
        ((double*)gradients->x)[nRows - 3] = 0;
        ((double*)gradients->x)[nRows - 2] = 0;
        ((double*)gradients->x)[nRows - 1] = 0;
        */

        symbolicFactorization(boundary_saddle);
        numericFactorization(boundary_saddle);

        gradients = cholmod_solve(CHOLMOD_A, symbolicFactor, temp_gradients, context);
        cholmod_free_dense(&temp_gradients, context);
    }

    void LWSBoundaryLoop::StepGradients(double h) {
        for (int i = 0; i < NumVertices(); i++) {
            int baseIndex = 3 * i;
            Vector3 gradient{
                ((double*)gradients->x)[baseIndex + 0],
                ((double*)gradients->x)[baseIndex + 1],
                ((double*)gradients->x)[baseIndex + 2]
            };

            Vertex vert = GetVertex(i);
            geom->vertexPositions[vert] = geom->vertexPositions[vert] + h * gradient;
        }

        ReparameterizeArc();
    }

    void LWSBoundaryLoop::ReparameterizeArc() {
        double currentLength = 0;
        arcLength.clear();
        int numVerts = NumVertices();
        // Record the arc length at each vertex of the loop
        for (int i = 0; i < numVerts; i++) {
            arcLength.push_back(currentLength);
            
            int next = (i + 1) % numVerts;
            Vector3 curPos = geom->vertexPositions[GetVertex(i)];
            Vector3 nextPos = geom->vertexPositions[GetVertex(next)];
            double distance = norm(curPos - nextPos);
            currentLength += distance;
        }
        // Record the total (cyclical) length of the loop
        totalLength = currentLength;
        cout << "Length of curve = " << totalLength << endl;
    }

    void LWSBoundaryLoop::CreateNodesAndElements() {
        // Create the nodes first
        for (size_t i = 0; i < vertices.size(); i++) {
            BoundaryNode* node_i = new BoundaryNode(i, geom->vertexPositions[vertices[i]], BisectorNormal(i));
            bem_nodes.push_back(node_i);
        }
        // Create elements connecting adjacent nodes
        for (size_t i = 0; i < vertices.size(); i++) {
            size_t i_plus = (i + 1) % vertices.size();
            BoundaryNode* node_i = bem_nodes[i];
            BoundaryNode* node_plus = bem_nodes[i_plus];
            BoundaryElement* element = new BoundaryElement(i, node_i, node_plus);
            bem_elements.push_back(element);

            node_i->element2 = element;
            node_plus->element1 = element;
        }

        boundary_u.setZero(vertices.size());
        boundary_v.setZero(vertices.size());
    }

    // Set up the right-hand vector for the system Au = Bq.
    // In this case, we provide the Dirichlet boundary condition u = 1
    // (the all-ones vector), so we compute A*1.
    void LWSBoundaryLoop::AllOnesRHS() {
        int rows = NumVertices();
        boundary_u.setZero(rows);
        for (int i = 0; i < rows; i++) {
            boundary_u(i) = 1;
        }
        rhs_vec = bemMatrixA * boundary_u;
    }
    
    // Set up and solve the BEM system Au = Bq for some Dirichlet boundary
    // conditions u.
    void LWSBoundaryLoop::LaplaceTestRHS() {
        int rows = NumVertices();
        boundary_u.setZero(rows);
        for (int i = 0; i < rows; i++) {
            Vector3 pos = geom->vertexPositions[vertices[i]];
            double value = geometrycentral::clamp(pos.x, -2.0, 2.0);
            if (value > 1) value = 1;
            else if (value < -0.5) value = -1;
            else value = 0;
            boundary_u(i) = value;
        }

        CreateBemMatrix();
        rhs_vec = bemMatrixA * boundary_u;
        SolveV();
    }

    // Get the vector of Dirichlet boundary conditions u as a dense vector.
    cholmod_dense* LWSBoundaryLoop::GetFullRHS(CholmodContext *c) {
        if (fullRHS) {
            cholmod_free_dense(&fullRHS, *c);
        }
        int rows = mesh->nVertices();
        VertexData<size_t> indices = mesh->getVertexIndices();
        // Insert the boundary entries in the appropriate places of the global vector
        fullRHS = cholmod_zeros(rows, 1, CHOLMOD_REAL, *c);
        for (int i = 0; i < NumVertices(); i++) {
            int vInd = indices[vertices[i]];
            ((double*)fullRHS->x)[vInd] = boundary_u(i);
        }
        return fullRHS;
    }

    void LWSBoundaryLoop::SolveV() {
        boundary_v = bemMatrixB.partialPivLu().solve(rhs_vec);
    }

    void LWSBoundaryLoop::MoveByV(double h) {
        for (int i = 0; i < NumVertices(); i++) {
            geom->vertexPositions[vertices[i]] += h * boundary_v(i) * BisectorNormal(i);
        }
        cout << "Moved by normal derivative" << endl;
    }

    void LWSBoundaryLoop::CreateBemMatrix() {
        bemMatrixA.setZero(NumVertices(), NumVertices());
        bemMatrixB.setZero(NumVertices(), NumVertices());

        for (BoundaryElement* element : bem_elements) {
            int colInd1 = element->node1->nodeID;
            int colInd2 = element->node2->nodeID;

            // Iterate over all collocation points (i.e. all nodes)
            for (int row = 0; row < NumVertices(); row++) {
                // Add the two columns for the two nodes to the global matrix
                double a_ia1 = element->CoeffIntA_ia(bem_nodes[row], bem_nodes[colInd1]);
                double a_ia2 = element->CoeffIntA_ia(bem_nodes[row], bem_nodes[colInd2]);
                bemMatrixA(row, colInd1) += a_ia1;
                bemMatrixA(row, colInd2) += a_ia2;

                double b_ia1 = element->CoeffIntB_ia(bem_nodes[row], bem_nodes[colInd1]);
                double b_ia2 = element->CoeffIntB_ia(bem_nodes[row], bem_nodes[colInd2]);
                bemMatrixB(row, colInd1) += b_ia1;
                bemMatrixB(row, colInd2) += b_ia2;
            }

            for (int row = 0; row < NumVertices(); row++) {
                // Add in the c_i term for each u_i
                bemMatrixA(row, row) += 0.5;
            }
        }
    }

    double LWSBoundaryLoop::DualLength(int x) {
        Vertex v = GetVertex(x);
        double length = 0;
        for (Edge e : v.adjacentEdges()) {
            if (e.isBoundary()) {
                length += geom->edgeLength(e);
            }
        }
        return length / 2;
    }

    double LWSBoundaryLoop::TotalLength() {
        return totalLength;
    }

    Vector3 LWSBoundaryLoop::Barycenter() {
        Vector3 center{0, 0, 0};
        double sumLength = 0;

        for (int i = 0; i < NumVertices(); i++) {
            Vertex v_i = GetVertex(i);
            double length = DualLength(i);
            center += length * geom->vertexPositions[v_i];
            sumLength += length;
        }

        center /= sumLength;
        return center;
    }

    double LWSBoundaryLoop::ValueAtInterior(Vector3 point) {
        // TODO: integrate u and v over boundary
        // Sum, for each element j:
        //   Sum over each node alpha:
        //     v(alpha) * j.coeff_b(pos, alpha) + u(alpha) * j.coeff_a(pos, alpha)
        double total = 0;
        
        for (BoundaryElement* &j : bem_elements) {
            int n1 = j->node1->nodeID;
            int n2 = j->node2->nodeID;

            total += (boundary_v(n1) * j->CoeffIntB_ia(point, j->node1) - boundary_u(n1) * j->CoeffIntA_ia(point, j->node1));
            total += (boundary_v(n2) * j->CoeffIntB_ia(point, j->node2) - boundary_u(n2) * j->CoeffIntA_ia(point, j->node2));
        }

        return total;
    }
}
