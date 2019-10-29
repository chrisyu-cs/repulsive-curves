#include "spatial/tpe_bvh.h"
#include <algorithm>

namespace LWS {

    int BVHNode3D::globalID = 0;

    inline double GetCoordFromBody(VertexBody6D body, int axis) {
        switch (axis) {
            case 0: return body.pt.position.x;
            case 1: return body.pt.position.y;
            case 2: return body.pt.position.z;
            case 3: return body.pt.tangent.x;
            case 4: return body.pt.tangent.y;
            case 5: return body.pt.tangent.z;
            default:
            std::cout << "Invalid axis passed to GetCoordFromBody: " << axis << std::endl;
            return body.pt.position.x;
        }
    }

    inline PosTan ReplaceCoord(PosTan body, int axis, double value) {
        switch (axis) {
            case 0: body.position.x = value; break;
            case 1: body.position.y = value; break;
            case 2: body.position.z = value; break;
            case 3: body.tangent.x = value; break;
            case 4: body.tangent.y = value; break;
            case 5: body.tangent.z = value; break;
            default:
            std::cout << "Invalid axis passed to ReplaceCoord: " << axis << std::endl; break;
        }
        return body;
    }

    inline int NextAxis(int axis) {
        switch (axis) {
            case 0: return 1;
            case 1: return 2;
            case 2: return 3;
            case 3: return 4;
            case 4: return 5;
            case 5: return 0;
            default: return 0;
        }
    }

    inline VertexBody6D vertToBody(PolyCurveNetwork* curves, int i) {
        CurveVertex* p = curves->GetVertex(i);
        Vector3 pos = p->Position();
        Vector3 tangent = p->Tangent();

        PosTan ptan{pos, tangent};
        double mass = p->DualLength();
        int globalIndex = p->GlobalIndex();

        return VertexBody6D{ptan, mass, globalIndex, BodyType::Vertex};
    }

    inline VertexBody6D edgeToBody(PolyCurveNetwork* curves, int i) {
        CurveEdge* edge = curves->GetEdge(i);

        Vector3 pos = edge->Midpoint();
        double mass = edge->Length();
        Vector3 tangent = edge->Tangent();

        PosTan ptan{pos, tangent};
        int gIndex1 = edge->GlobalIndex();

        return VertexBody6D{ptan, mass, gIndex1, BodyType::Edge};
    }

    BVHNode3D* CreateBVHFromCurve(PolyCurveNetwork *curves) {
        int nVerts = curves->NumVertices();
        std::vector<VertexBody6D> verts(nVerts);

        // Loop over all the vertices
        for (int i = 0; i < nVerts; i++) {
            VertexBody6D curBody = vertToBody(curves, i);
            // Put vertex body into full list
            verts[curBody.elementIndex] = curBody;
        }

        BVHNode3D* tree = new BVHNode3D(verts, 0, 0);
        tree->recomputeCentersOfMass(curves);
        BVHNode3D::globalID = 0;
        tree->recursivelyAssignIDs();

        std::cout << "Created vertex BVH with " << BVHNode3D::globalID << " nodes" << std::endl;
        tree->numNodes = BVHNode3D::globalID;

        return tree;
    }

    BVHNode3D* CreateEdgeBVHFromCurve(PolyCurveNetwork *curves) {
        int nEdges = curves->NumEdges();
        std::vector<VertexBody6D> verts(nEdges);

        // Loop over all the vertices
        for (int i = 0; i < nEdges; i++) {
            VertexBody6D curBody = edgeToBody(curves, i);
            // Put vertex body into full list
            verts[curBody.elementIndex] = curBody;
        }

        BVHNode3D* tree = new BVHNode3D(verts, 0, 0);
        tree->recomputeCentersOfMass(curves);
        BVHNode3D::globalID = 0;
        tree->recursivelyAssignIDs();

        std::cout << "Created edge BVH with " << BVHNode3D::globalID << " nodes" << std::endl;
        tree->numNodes = BVHNode3D::globalID;

        return tree;
    }

    BVHNode3D::BVHNode3D(std::vector<VertexBody6D> &points, int axis, BVHNode3D* root) {
        // Split the points into sets somehow
        thresholdTheta = 0.25;
        splitAxis = axis;
        zeroMVFields();

        if (!root) bvhRoot = this;
        else bvhRoot = root;

        if (points.size() == 0) {
            isLeaf = false;
            isEmpty = true;
            totalMass = 0;
            numElements = 0;
            // This is just to find bugs more easily
            body.elementIndex = -999;
        }
        // If there's only one point, then we don't need to do anything
        // except set the fields from that one point
        else if (points.size() == 1) {
            body = points[0];
            isLeaf = true;
            isEmpty = false;
            totalMass = body.mass;
            centerOfMass = body.pt.position;
            averageTangent = body.pt.tangent;
            minCoords = body.pt;
            maxCoords = body.pt;
            numElements = 1;

            clusterIndices.push_back(body.elementIndex);
        }
        else {
            // Reserve space for splitting the points into lesser and greater
            int nPoints = points.size();
            std::vector<VertexBody6D> lesserPoints;
            lesserPoints.reserve(nPoints / 2 + 1);
            std::vector<VertexBody6D> greaterPoints;
            greaterPoints.reserve(nPoints / 2 + 1);

            // Compute the plane over which to split the points
            splitPoint = AxisSplittingPlane(points, axis);

            // Split the points over the median
            for (int i = 0; i < nPoints; i++) {
                double coord = GetCoordFromBody(points[i], axis);

                if (coord <= splitPoint) {
                    lesserPoints.push_back(points[i]);
                }
                else {
                    greaterPoints.push_back(points[i]);
                }
            }

            // Compute the middle extents for the two halves

            //std::cout << "Split " << nPoints << " into " << lesserPoints.size() << " and " << greaterPoints.size() << std::endl;

            // Recursively construct children
            int nextAxis = NextAxis(axis);

            BVHNode3D* nextRoot = (root) ? root : this;
            BVHNode3D* lesserNode = new BVHNode3D(lesserPoints, nextAxis, nextRoot);
            BVHNode3D* greaterNode = new BVHNode3D(greaterPoints, nextAxis, nextRoot);

            children.push_back(lesserNode);
            children.push_back(greaterNode);

            clusterIndices.insert(clusterIndices.end(), lesserNode->clusterIndices.begin(), lesserNode->clusterIndices.end());
            clusterIndices.insert(clusterIndices.end(), greaterNode->clusterIndices.begin(), greaterNode->clusterIndices.end());

            if (!root) {
                fullMasses.setZero(points.size());
                for (size_t i = 0; i < points.size(); i++) {
                    fullMasses(i) = points[i].mass;
                }
            }

            isLeaf = false;
            isEmpty = false;
        }
    }

    BVHNode3D::~BVHNode3D() {
        for (size_t i = 0; i < children.size(); i++) {
            if (children[i]) delete children[i];
        }
    }

    void BVHNode3D::findCurveSegments(std::vector<VertexBody6D> &points, PolyCurveNetwork* curves) {
        // Trace out all of the connected curve segments contained
        // within the given point set.
        std::cerr << "findCurveSegments not implemented" << std::endl;
        throw 1;
    }

    double BVHNode3D::AxisSplittingPlane(std::vector<VertexBody6D> &points, int axis) {
        size_t nPoints = points.size();
        std::vector<double> coords(nPoints);

        for (size_t i = 0; i < nPoints; i++) {
            coords[i] = GetCoordFromBody(points[i], axis);
        }

        std::sort(coords.begin(), coords.end());

        size_t splitIndex = -1;
        double minWidths = INFINITY;

        for (size_t i = 0; i < nPoints; i++) {
            double width1 = coords[i] - coords[0];
            double width2 = (i == nPoints - 1) ? 0 : coords[nPoints - 1] - coords[i + 1];

            double sumSquares = width1 * width1 + width2 * width2;
            if (sumSquares < minWidths) {
                minWidths = sumSquares;
                splitIndex = i;
            }
        }

        double splitPoint = (coords[splitIndex] + coords[splitIndex + 1]) / 2;
        return splitPoint;
    }

    void BVHNode3D::setLeafData(PolyCurveNetwork* curves) {
        if (body.type == BodyType::Vertex) {
            CurveVertex* p = curves->GetVertex(body.elementIndex);
            body.mass = p->DualLength();
            body.pt.position = p->Position();
            body.pt.tangent = p->Tangent();
        }
        else if (body.type == BodyType::Edge) {
            CurveEdge* p1 = curves->GetEdge(body.elementIndex);

            // Mass of an edge is its length
            body.mass = p1->Length();
            // Use midpoint as center of mass
            body.pt.position = p1->Midpoint();
            // Tangent direction is normalized edge vector
            body.pt.tangent = p1->Tangent();
        }

        totalMass = body.mass;
        centerOfMass = body.pt.position;
        averageTangent = body.pt.tangent;

        minCoords = PosTan{body.pt.position, body.pt.tangent};
        maxCoords = minCoords;
    }

    void BVHNode3D::refreshWeightsVector(PolyCurveNetwork* curves, BodyType bType) {
        if (bType == BodyType::Vertex) {
            int nVerts = curves->NumVertices();
            for (int i = 0; i < nVerts; i++) {
                CurveVertex* v = curves->GetVertex(i);
                fullMasses(v->id) = v->DualLength();
            }
        }
        else if (bType == BodyType::Edge) {
            int nEdges = curves->NumEdges();
            for (int i = 0; i < nEdges; i++) {
                CurveEdge* e = curves->GetEdge(i);
                fullMasses(e->id) = e->Length();
            }
        }
    }

    void BVHNode3D::recomputeCentersOfMass(PolyCurveNetwork* curves) {
        if (isEmpty) {
            totalMass = 0;
            numElements = 0;
        }
        // For a leaf, just set centers and bounds from the one body
        else if (isLeaf) {
            setLeafData(curves);
            numElements = 1;
        }
        else {
            // Recursively compute bounds for all children
            for (size_t i = 0; i < children.size(); i++) {
                children[i]->recomputeCentersOfMass(curves);
            }

            minCoords = children[0]->minCoords;
            maxCoords = children[0]->maxCoords;

            totalMass = 0;
            centerOfMass = Vector3{0, 0, 0};
            averageTangent = Vector3{0, 0, 0};
            
            // Accumulate max/min over all nonempty children
            for (size_t i = 0; i < children.size(); i++) {
                if (!children[i]->isEmpty) {
                    minCoords = postan_min(children[i]->minCoords, minCoords);
                    maxCoords = postan_max(children[i]->maxCoords, maxCoords);

                    totalMass += children[i]->totalMass;
                    centerOfMass += children[i]->centerOfMass * children[i]->totalMass;
                    averageTangent += children[i]->averageTangent * children[i]->totalMass;
                }
            }

            centerOfMass /= totalMass;
            averageTangent /= totalMass;

            averageTangent = averageTangent.normalize();

            numElements = 0;
            for (size_t i = 0; i < children.size(); i++) {
                numElements += children[i]->numElements;
            }
        }
    }

    int BVHNode3D::NumElements() {
        return numElements;
    }

    inline double BVHNode3D::nodeRadius() {
        // Compute diagonal distance from corner to corner
        double diag = norm(maxCoords.position - minCoords.position);
        return diag;
    }

    Vector2 BVHNode3D::viewspaceBounds(Vector3 point) {
        Vector3 center = (maxCoords.position + minCoords.position) / 2;
        Vector3 offset = (maxCoords.position - minCoords.position) / 2;

        Vector3 toCenter = (centerOfMass - point);
        toCenter = toCenter.normalize();

        Vector3 corners[8] = {
            Vector3{center.x - offset.x, center.y - offset.y, center.z - offset.z},
            Vector3{center.x - offset.x, center.y - offset.y, center.z + offset.z},
            Vector3{center.x - offset.x, center.y + offset.y, center.z - offset.z},
            Vector3{center.x - offset.x, center.y + offset.y, center.z + offset.z},
            Vector3{center.x + offset.x, center.y - offset.y, center.z - offset.z},
            Vector3{center.x + offset.x, center.y - offset.y, center.z + offset.z},
            Vector3{center.x + offset.x, center.y + offset.y, center.z - offset.z},
            Vector3{center.x + offset.x, center.y + offset.y, center.z + offset.z}
        };

        // Measure distance along the offset to the point
        double minLinear = dot(toCenter, corners[0] - point);
        double maxLinear = minLinear;

        // Measure distance projected onto the plane orthogonal to the offset from the point
        Vector3 offset0 = corners[0] - point;
        double maxRadial = norm(offset0 - toCenter * dot(toCenter, offset0));

        for (int i = 0; i < 8; i++) {
            double linearDist = dot(toCenter, corners[i] - point);
            minLinear = fmin(minLinear, linearDist);
            maxLinear = fmax(maxLinear, linearDist);

            Vector3 offset = corners[i] - point;
            double radialDist = norm(offset - toCenter * dot(toCenter, offset));
            maxRadial = fmax(maxRadial, radialDist);
        }

        double linearSpread = maxLinear - minLinear;
        return Vector2{maxRadial, linearSpread};
    }

    bool BVHNode3D::shouldUseCell(Vector3 vertPos) {
        double d = norm(centerOfMass - vertPos);
        Vector2 ratios = viewspaceBounds(vertPos) / d;
        // TODO: take into account some tangent-related criteria?
        return fmax(ratios.x, ratios.y) < thresholdTheta;
        // return (nodeRadius() / d) < thresholdTheta;
    }

    void BVHNode3D::accumulateVertexEnergy(double &result, CurveVertex* &i_pt,
    PolyCurveNetwork* curves, double alpha, double beta) {
        if (isEmpty) {
            return;
        }
        else if (isLeaf) {
            // If this is a leaf, then it only has one element in it, so just use it
            if (body.type == BodyType::Vertex) {
                // If the element is a vertex, just get that vertex
                CurveVertex* j = curves->GetVertex(body.elementIndex);
                // Don't sum if it's the same vertex
                if (j == i_pt) return;
                // Add contribution to energy at i from vertex j
                result += TPESC::tpe_pair(i_pt, j, alpha, beta);
            }
            else if (body.type == BodyType::Edge) {
                // Otherwise the element is an edge, so we use its midpoint (stored in the body)
                result += TPESC::tpe_pair_pts(i_pt->Position(), body.pt.position, i_pt->Tangent(), i_pt->DualLength(), body.mass, alpha, beta);
            }
        }
        else {
            if (shouldUseCell(i_pt->Position())) {
                // This cell is far enough away that we can treat it as a single body
                result += bodyEnergyEvaluation(i_pt, alpha, beta);
            }
            else {
                // Otherwise we continue recursively traversing the tree
                for (size_t i = 0; i < children.size(); i++) {
                    if (children[i]) {
                        children[i]->accumulateVertexEnergy(result, i_pt, curves, alpha, beta);
                    }
                }
            }
        }
    }

    double BVHNode3D::bodyEnergyEvaluation(CurveVertex* &i_pt, double alpha, double beta) {
        Vector3 tangent = averageTangent;
        tangent = tangent.normalize();
        return TPESC::tpe_pair_pts(i_pt->Position(), centerOfMass, tangent, i_pt->DualLength(), totalMass, alpha, beta);
    }

    void BVHNode3D::accumulateTPEGradient(Eigen::MatrixXd &gradients, CurveVertex* &i_pt,
    PolyCurveNetwork* curves, double alpha, double beta) {
        if (isEmpty) {
            return;
        }
        else if (isLeaf) {
            // With a vertex, we add gradient terms the same way as usual
            if (body.type == BodyType::Vertex) {
                // If this is a leaf, then it only has one vertex in it, so just use it
                CurveVertex* j_pt = curves->GetVertex(body.elementIndex);
                // Don't sum if it's the same vertex
                if (j_pt == i_pt) return;

                // Add i and neighbors of i
                std::vector<CurveVertex*> i_pts;
                i_pts.push_back(i_pt);
                for (int e = 0; e < i_pt->numEdges(); e++) {
                    i_pts.push_back(i_pt->edge(e)->Opposite(i_pt));
                }

                // Add j and neighbors of j
                std::vector<CurveVertex*> j_pts;
                j_pts.push_back(j_pt);
                for (int e = 0; e < j_pt->numEdges(); e++) {
                    j_pts.push_back(j_pt->edge(e)->Opposite(j_pt));
                }

                for (CurveVertex* i_n : i_pts) {
                    AddToRow(gradients, i_n->GlobalIndex(), TPESC::tpe_grad(i_pt, j_pt, alpha, beta, i_n));
                    bool noOverlap = true;
                    // Avoid double-counting on the reverse terms with these checks
                    for (CurveVertex* j_n : j_pts) {
                        if (i_n == j_n) noOverlap = false;
                    }
                    if (noOverlap) {
                        AddToRow(gradients, i_n->GlobalIndex(), TPESC::tpe_grad(j_pt, i_pt, alpha, beta, i_n));
                    }
                }
            }
            // With an edge, we have to pass in a TangentMassPoint instead
            else if (body.type == BodyType::Edge) {
                Vector3 tangent = body.pt.tangent;
                tangent = tangent.normalize();
                CurveEdge* e = curves->GetEdge(body.elementIndex);
                CurveVertex* j1 = e->prevVert;
                CurveVertex* j2 = e->nextVert;

                // Add i and neighbors of i
                std::vector<CurveVertex*> i_pts;
                i_pts.push_back(i_pt);
                for (int e = 0; e < i_pt->numEdges(); e++) {
                    i_pts.push_back(i_pt->edge(e)->Opposite(i_pt));
                }

                TangentMassPoint jm{tangent, body.mass, body.pt.position, j1, j2};

                if (i_pt != j1 && i_pt != j2) {
                    for (CurveVertex* i_n : i_pts) {
                        AddToRow(gradients, i_n->GlobalIndex(), TPESC::tpe_grad(i_pt, jm, alpha, beta, i_n));
                        AddToRow(gradients, i_n->GlobalIndex(), TPESC::tpe_grad(jm, i_pt, alpha, beta, i_n));
                    }
                }
            }
        }
        else {
            if (shouldUseCell(i_pt->Position())) {
                Vector3 tangent = averageTangent;
                tangent = tangent.normalize();
                // This cell is far enough away that we can treat it as a single body
                TangentMassPoint j{tangent, totalMass, centerOfMass, 0, 0};

                // Add i and neighbors of i
                std::vector<CurveVertex*> i_pts;
                i_pts.push_back(i_pt);
                for (int e = 0; e < i_pt->numEdges(); e++) {
                    i_pts.push_back(i_pt->edge(e)->Opposite(i_pt));
                }

                // Differentiate both terms for previous, middle, and next
                for (CurveVertex* i_n : i_pts) {
                    AddToRow(gradients, i_n->GlobalIndex(), TPESC::tpe_grad(i_pt, j, alpha, beta, i_n));
                    AddToRow(gradients, i_n->GlobalIndex(), TPESC::tpe_grad(j, i_pt, alpha, beta, i_n));
                }
            }
            else {
                // Otherwise we continue recursively traversing the tree
                for (size_t i = 0; i < children.size(); i++) {
                    if (children[i]) {
                        children[i]->accumulateTPEGradient(gradients, i_pt, curves, alpha, beta);
                    }
                }
            }
        }
    }

    Vector3 BVHNode3D::bodyForceEvaluation(CurveVertex* &i_pt, double alpha, double beta) {
        // TODO: placeholder
        return Vector3{0, 0, 0};
    }

    Vector3 BVHNode3D::exactGradient(CurveVertex* i_pt, PolyCurveNetwork* curves, double alpha, double beta) {
        if (isEmpty) {
            return Vector3{0, 0, 0};
        }
        else if (isLeaf) {
            // If this is a leaf, then it only has one vertex in it, so just use it
            CurveVertex* j_pt = curves->GetVertex(body.elementIndex);
            // Don't sum if it's the same vertex
            if (j_pt == i_pt) return Vector3{0, 0, 0};

            Vector3 deriv{0, 0, 0};

            deriv += TPESC::tpe_grad(i_pt, j_pt, alpha, beta, i_pt);
            deriv += TPESC::tpe_grad(j_pt, i_pt, alpha, beta, i_pt);

            return deriv;
        }
        else {
            Vector3 total{0, 0, 0};
            // Otherwise we continue recursively traversing the tree
            for (size_t i = 0; i < children.size(); i++) {
                if (children[i]) {
                    total += children[i]->exactGradient(i_pt, curves, alpha, beta);
                }
            }
            return total;
        }
    }

    PosTan BVHNode3D::minBound() {
        return minCoords;
    }
    
    PosTan BVHNode3D::maxBound() {
        return maxCoords;
    }

    Vector3 BVHNode3D::BoxCenter() {
        return (minCoords.position + maxCoords.position) / 2;
    }

    void BVHNode3D::accumulateChildren(std::vector<VertexBody6D> &result) {
        if (isEmpty) {
            return;
        }
        else if (isLeaf) {
            result.push_back(body);
        }
        else {
            for (BVHNode3D* child : children) {
                child->accumulateChildren(result);
            }
        }
    }
}