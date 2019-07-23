#include "spatial/tpe_bvh.h"
#include <algorithm>

namespace LWS {

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

    inline VertexBody6D vertToBody(PolyCurveGroup* curves, PolyCurve* curve, int i) {
        PointOnCurve p = curve->GetCurvePoint(i);
        Vector3 pos = p.Position();
        Vector3 tangent = p.Tangent();

        PosTan ptan{pos, tangent};
        double mass = p.DualLength();
        int globalIndex = curves->GlobalIndex(p);

        return VertexBody6D{ptan, mass, globalIndex, -1};
    }

    inline VertexBody6D edgeToBody(PolyCurveGroup* curves, PolyCurve* curve, int i) {
        PointOnCurve p1 = curve->GetCurvePoint(i);
        PointOnCurve p2 = p1.Next();

        Vector3 x1 = p1.Position();
        Vector3 x2 = p2.Position();

        Vector3 pos = (x1 + x2) / 2;
        double mass = norm(x2 - x1);
        Vector3 tangent = (x2 - x1) / mass;

        PosTan ptan{pos, tangent};
        int gIndex1 = curves->GlobalIndex(p1);
        int gIndex2 = curves->GlobalIndex(p2);

        return VertexBody6D{ptan, mass, gIndex1, gIndex2};
    }

    BVHNode3D* CreateBVHFromCurve(PolyCurveGroup *curves) {
        std::vector<VertexBody6D> verts(curves->NumVertices());

        // Loop over all the vertices
        for (size_t i = 0; i < curves->curves.size(); i++) {
            PolyCurve* curve = curves->curves[i];
            size_t nVerts = curve->NumVertices();
            for (size_t v = 0; v < nVerts; v++) {
                VertexBody6D curBody = vertToBody(curves, curve, v);
                // Put vertex body into full list
                verts[curBody.vertIndex1] = curBody;
            }
        }

        BVHNode3D* tree = new BVHNode3D(verts, 0);
        tree->recomputeCentersOfMass(curves);

        return tree;
    }

    BVHNode3D* CreateEdgeBVHFromCurve(PolyCurveGroup *curves) {
        std::vector<VertexBody6D> verts(curves->NumVertices());

        // Loop over all the vertices
        for (size_t i = 0; i < curves->curves.size(); i++) {
            PolyCurve* curve = curves->curves[i];
            size_t nVerts = curve->NumVertices();
            for (size_t v = 0; v < nVerts; v++) {
                VertexBody6D curBody = edgeToBody(curves, curve, v);
                // Put vertex body into full list
                verts[curBody.vertIndex1] = curBody;
            }
        }

        BVHNode3D* tree = new BVHNode3D(verts, 0);
        tree->recomputeCentersOfMass(curves);

        return tree;
    }

    BVHNode3D::BVHNode3D(std::vector<VertexBody6D> &points, int axis) {
        // Split the points into sets somehow
        thresholdTheta = 0.25;
        splitAxis = axis;

        if (points.size() == 0) {
        
            isLeaf = true;
            isEmpty = true;
            totalMass = 0;
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
            BVHNode3D* lesserNode = new BVHNode3D(lesserPoints, nextAxis);
            BVHNode3D* greaterNode = new BVHNode3D(greaterPoints, nextAxis);

            children.push_back(lesserNode);
            children.push_back(greaterNode);
            
            isLeaf = false;
            isEmpty = false;
        }
    }

    BVHNode3D::~BVHNode3D() {
        for (size_t i = 0; i < children.size(); i++) {
            if (children[i]) delete children[i];
        }
    }

    void BVHNode3D::findCurveSegments(std::vector<VertexBody6D> &points, PolyCurveGroup* curves) {
        // Trace out all of the connected curve segments contained
        // within the given point set.
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

    void BVHNode3D::setLeafData(PolyCurveGroup* curves) {
        if (body.type() == BodyType::Vertex) {
            PointOnCurve p = curves->GetCurvePoint(body.vertIndex1);

            body.mass = p.DualLength();
            body.pt.position = p.Position();
            body.pt.tangent = p.Tangent();
        }
        else if (body.type() == BodyType::Edge) {
            PointOnCurve p1 = curves->GetCurvePoint(body.vertIndex1);
            PointOnCurve p2 = curves->GetCurvePoint(body.vertIndex2);

            Vector3 x1 = p1.Position();
            Vector3 x2 = p2.Position();

            // Mass of an edge is its length
            body.mass = norm(x2 - x1);
            // Use midpoint as center of mass
            body.pt.position = (x1 + x2) / 2;
            // Tangent direction is normalized edge vector
            body.pt.tangent = (x2 - x1) / body.mass;
        }

        totalMass = body.mass;
        centerOfMass = body.pt.position;
        averageTangent = body.pt.tangent;

        minCoords = PosTan{body.pt.position, body.pt.tangent};
        maxCoords = minCoords;
    }

    void BVHNode3D::recomputeCentersOfMass(PolyCurveGroup* curves) {
        if (isEmpty) {
            totalMass = 0;
        }
        // For a leaf, just set centers and bounds from the one body
        else if (isLeaf) {
            setLeafData(curves);
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
        }
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

    void BVHNode3D::accumulateVertexEnergy(double &result, PointOnCurve &i_pt,
    PolyCurveGroup* curves, double alpha, double beta) {
        if (isEmpty) {
            return;
        }
        else if (isLeaf) {
            // If this is a leaf, then it only has one element in it, so just use it
            if (body.type() == BodyType::Vertex) {
                // If the element is a vertex, just get that vertex
                PointOnCurve j = curves->GetCurvePoint(body.vertIndex1);
                // Don't sum if it's the same vertex
                if (j == i_pt) return;
                // Add contribution to energy at i from vertex j
                result += TPESC::tpe_pair(i_pt, j, alpha, beta);
            }
            else if (body.type() == BodyType::Edge) {
                // Otherwise the element is an edge, so we use its midpoint (stored in the body)
                result += TPESC::tpe_pair_pts(i_pt.Position(), body.pt.position, i_pt.Tangent(), i_pt.DualLength(), body.mass, alpha, beta);
            }
        }
        else {
            if (shouldUseCell(i_pt.Position())) {
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

    double BVHNode3D::bodyEnergyEvaluation(PointOnCurve &i_pt, double alpha, double beta) {
        Vector3 tangent = averageTangent;
        tangent = tangent.normalize();
        return TPESC::tpe_pair_pts(i_pt.Position(), centerOfMass, tangent, i_pt.DualLength(), totalMass, alpha, beta);
    }

    void BVHNode3D::accumulateTPEGradient(std::vector<Vector3> &gradients, PointOnCurve &i_pt,
    PolyCurveGroup* curves, double alpha, double beta) {
        if (isEmpty) {
            return;
        }
        else if (isLeaf) {
            // With a vertex, we add gradient terms the same way as usual
            if (body.type() == BodyType::Vertex) {
                // If this is a leaf, then it only has one vertex in it, so just use it
                PointOnCurve j = curves->GetCurvePoint(body.vertIndex1);
                // Don't sum if it's the same vertex
                if (j == i_pt) return;

                PointOnCurve i_prev = i_pt.Prev();
                PointOnCurve i_next = i_pt.Next();

                PointOnCurve j_prev = j.Prev();
                PointOnCurve j_next = j.Next();

                // Differentiate both terms for previous, middle, and next
                gradients[curves->GlobalIndex(i_prev)] += TPESC::tpe_grad(i_pt, j, alpha, beta, i_prev);
                gradients[curves->GlobalIndex(i_pt)]   += TPESC::tpe_grad(i_pt, j, alpha, beta, i_pt);
                gradients[curves->GlobalIndex(i_next)] += TPESC::tpe_grad(i_pt, j, alpha, beta, i_next);
                // Avoid double-counting on the reverse terms with these checks
                if (i_prev != j_prev && i_prev != j && i_prev != j_next) {
                    gradients[curves->GlobalIndex(i_prev)] += TPESC::tpe_grad(j, i_pt, alpha, beta, i_prev);
                }
                if (i_pt != j_prev && i_pt != j && i_pt != j_next) {
                    gradients[curves->GlobalIndex(i_pt)]   += TPESC::tpe_grad(j, i_pt, alpha, beta, i_pt);
                }
                if (i_next != j_prev && i_next != j && i_next != j_next) {
                    gradients[curves->GlobalIndex(i_next)] += TPESC::tpe_grad(j, i_pt, alpha, beta, i_next);
                }
            }
            // With an edge, we have to pass in a TangentMassPoint instead
            else if (body.type() == BodyType::Edge) {
                Vector3 tangent = body.pt.tangent;
                tangent = tangent.normalize();

                PointOnCurve j1 = curves->GetCurvePoint(body.vertIndex1);
                PointOnCurve j2 = curves->GetCurvePoint(body.vertIndex2);

                PointOnCurve i_prev = i_pt.Prev();
                PointOnCurve i_next = i_pt.Next();

                TangentMassPoint jm{tangent, body.mass, body.pt.position, j1, j2};

                if (i_pt != j1 && i_pt != j2) {
                    gradients[curves->GlobalIndex(i_prev)] += TPESC::tpe_grad(i_pt, jm, alpha, beta, i_prev);
                    gradients[curves->GlobalIndex(i_pt)]   += TPESC::tpe_grad(i_pt, jm, alpha, beta, i_pt);
                    gradients[curves->GlobalIndex(i_next)] += TPESC::tpe_grad(i_pt, jm, alpha, beta, i_next);

                    gradients[curves->GlobalIndex(i_prev)] += TPESC::tpe_grad(jm, i_pt, alpha, beta, i_prev);
                    gradients[curves->GlobalIndex(i_pt)]   += TPESC::tpe_grad(jm, i_pt, alpha, beta, i_pt);
                    gradients[curves->GlobalIndex(i_next)] += TPESC::tpe_grad(jm, i_pt, alpha, beta, i_next);
                }
            }
        }
        else {
            if (shouldUseCell(i_pt.Position())) {
                Vector3 tangent = averageTangent;
                tangent = tangent.normalize();
                // This cell is far enough away that we can treat it as a single body
                TangentMassPoint j{tangent, totalMass, centerOfMass, PointOnCurve{0, 0}, PointOnCurve{0, 0}};

                PointOnCurve i_prev = i_pt.Prev();
                PointOnCurve i_next = i_pt.Next();

                // Differentiate both terms for previous, middle, and next
                gradients[curves->GlobalIndex(i_prev)] += TPESC::tpe_grad(i_pt, j, alpha, beta, i_prev);
                gradients[curves->GlobalIndex(i_prev)] += TPESC::tpe_grad(j, i_pt, alpha, beta, i_prev);
                gradients[curves->GlobalIndex(i_pt)]   += TPESC::tpe_grad(i_pt, j, alpha, beta, i_pt);
                gradients[curves->GlobalIndex(i_pt)]   += TPESC::tpe_grad(j, i_pt, alpha, beta, i_pt);
                gradients[curves->GlobalIndex(i_next)] += TPESC::tpe_grad(i_pt, j, alpha, beta, i_next);
                gradients[curves->GlobalIndex(i_next)] += TPESC::tpe_grad(j, i_pt, alpha, beta, i_next);
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

    Vector3 BVHNode3D::bodyForceEvaluation(PointOnCurve &i_pt, double alpha, double beta) {
        // TODO: placeholder
        return Vector3{0, 0, 0};
    }

    void BVHNode3D::testGradientSingle(std::vector<BHPlotData> &out, PointOnCurve i_pt,
    PolyCurveGroup* curves, double alpha, double beta) {
        if (isEmpty) {
            return;
        }
        else if (isLeaf) {
            // In this case, no error is possible
            return;
        }
        else {
            if (shouldUseCell(i_pt.Position())) {
                Vector3 tangent = averageTangent;
                tangent = tangent.normalize();
                // This cell is far enough away that we can treat it as a single body
                TangentMassPoint j{tangent, totalMass, centerOfMass, PointOnCurve{0, 0}, PointOnCurve{0, 0}};

                PointOnCurve i_prev = i_pt.Prev();
                PointOnCurve i_next = i_pt.Next();

                Vector3 approx{0, 0, 0};
                approx += TPESC::tpe_grad(i_pt, j, alpha, beta, i_pt);
                approx += TPESC::tpe_grad(j, i_pt, alpha, beta, i_pt);

                Vector3 exact = exactGradient(i_pt, curves, alpha, beta);

                double normDiff = norm(approx - exact);
                double pctDiff = 100 * (normDiff / norm(exact));
                double gradNorm = norm(exact);

                Vector2 bounds = viewspaceBounds(i_pt.Position());
                double d = norm(centerOfMass - i_pt.Position());

                double radialTheta = bounds.x / d;
                double linearTheta = bounds.y / d;
                double minWidth = fmin(bounds.x, bounds.y);
                double maxWidth = fmax(bounds.x, bounds.y);

                out.push_back(BHPlotData{fmax(radialTheta, linearTheta), pctDiff, gradNorm, minWidth, maxWidth});
            }
            else {
                // Otherwise we continue recursively traversing the tree
                for (size_t i = 0; i < children.size(); i++) {
                    if (children[i]) {
                        children[i]->testGradientSingle(out, i_pt, curves, alpha, beta);
                    }
                }
            }
        }
    }

    Vector3 BVHNode3D::exactGradient(PointOnCurve i_pt, PolyCurveGroup* curves, double alpha, double beta) {
        if (isEmpty) {
            return Vector3{0, 0, 0};
        }
        else if (isLeaf) {
            // If this is a leaf, then it only has one vertex in it, so just use it
            PointOnCurve j_pt = curves->GetCurvePoint(body.vertIndex1);
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